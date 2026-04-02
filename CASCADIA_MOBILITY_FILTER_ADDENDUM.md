# CASCADIA_MOBILITY_FILTER_ADDENDUM.md
# Mobility-Filtered Spectrum Extraction for Cascadia De Novo Sequencing
# Version: 1.0 | Priority: P1 (modifies Phase 1 of CASCADIA_DENOVO_SPEC.md)

---

## 1. Problem Statement

The current `bruker_augment.py` spec (Phase 1 of CASCADIA_DENOVO_SPEC.md) reads timsTOF
`.d` files via `timsrust_pyo3` and **flattens the ion mobility dimension** — summing all
ions across the full 1/K0 range within each diaPASEF isolation window. This is equivalent
to what mzML conversion does.

The problem: diaPASEF isolation windows are wide (typically 25–50 m/z) and span a
significant 1/K0 range. Summing across the full mobility range produces chimeric MS2
spectra contaminated with fragments from co-isolated precursors at different mobilities.

For **database search** (DIA-NN), this is tolerable — the search engine uses its spectral
library to deconvolve chimeric spectra via scoring. For **de novo sequencing** (Cascadia),
there is no library to rescue chimeric input. The transformer model must predict sequences
directly from the fragment pattern, and chimeric contamination directly degrades prediction
confidence and accuracy.

### The insight

timsTOF's entire value proposition is PASEF — using ion mobility to achieve near-DDA
selectivity in a DIA acquisition. But if we flatten the mobility dimension before
handing spectra to Cascadia, we throw away that selectivity advantage. The spectra
become equivalent to wide-window Orbitrap DIA — noisy and chimeric.

**Spectronaut solves this** by using ion mobilograms (intensity vs. 1/K0 profiles) to
extract mobility-filtered spectra for each precursor. We can do the same thing in the
extraction step, producing cleaner pseudo-DDA-quality spectra for Cascadia.

---

## 2. timsTOF Data Model Refresher

```
Frame (one RT time point, one full TIMS ramp)
  └── Scans (each scan = one 1/K0 slice of the TIMS ramp)
        └── Peaks (m/z, intensity pairs within that scan)
```

- A **frame** is all data acquired during one TIMS ramp cycle (~100 ms). All peaks
  in a frame share the same retention time.
- Each frame contains ~700–1000 **scans**, each at a different 1/K0 value.
  Scan number maps to 1/K0 via a calibration function stored in `analysis.tdf`.
- In diaPASEF, a single frame contains multiple isolation windows at different
  m/z × 1/K0 positions (defined in `DiaFrameMsMsWindows` table).

### Key SQLite tables in `analysis.tdf`

| Table | Contents |
|---|---|
| `Frames` | One row per frame: Id, Time, MsMsType, SummedIntensities, etc. |
| `DiaFrameMsMsWindows` | diaPASEF window definitions: Frame → WindowGroup → IsolationMz, IsolationWidth, ScanNumBegin, ScanNumEnd |
| `GlobalMetadata` | Calibration coefficients for scan → 1/K0 conversion |

The `ScanNumBegin` / `ScanNumEnd` columns in `DiaFrameMsMsWindows` define the
1/K0 boundaries of each isolation window — this is the mobility range that was
actively selected by the quadrupole for that window.

### timsrust access model

timsrust exposes two data types:
- **Spectra**: Pre-summed (mobility-collapsed) — what the current spec uses
- **Frames**: Full scan-level data preserving the 1/K0 dimension

For mobility-filtered extraction, we need **Frame-level access**.

---

## 3. Mobility-Filtered Extraction Strategy

### 3.1 Overview

Instead of summing all scans within an isolation window (current approach), we:

1. Get precursor targets from DIA-NN's `report.tsv` (which includes `IM` = 1/K0)
2. For each precursor, read the raw frame data at the precursor's RT
3. Extract only the scans within a narrow 1/K0 window around the precursor's
   expected mobility
4. Build the MS2 spectrum from those filtered scans only
5. Hand the cleaned spectrum to Cascadia

This is conceptually identical to what Spectronaut does with its mobilogram-based
extraction, but we're using DIA-NN's precursor detection as the front-end instead
of Spectronaut's proprietary algorithm.

### 3.2 Two operating modes

**Mode A — DIA-NN-guided extraction (recommended)**

Uses DIA-NN's `report.tsv` as a precursor target list. DIA-NN has already identified
precursors with their m/z, RT, charge, and 1/K0 values. We extract mobility-filtered
spectra only for these targets.

Advantages:
- Precise 1/K0 targets from DIA-NN's scoring
- Only extracts spectra for real precursors (not noise)
- Natural integration with the parallel DIA-NN + Cascadia pipeline
- Cascadia confirms/refutes DIA-NN's identifications with cleaner input

Disadvantage:
- Requires DIA-NN to finish first (breaks the parallel job model)
- Misses precursors DIA-NN didn't detect (defeats part of the de novo purpose)

**Mode B — Window-based extraction with mobility binning (no DIA-NN dependency)**

For each diaPASEF isolation window, instead of summing all scans, we:
1. Read all scans in the window's `ScanNumBegin:ScanNumEnd` range
2. Build a mobilogram (summed intensity vs. scan number / 1/K0)
3. Detect peaks in the mobilogram (local maxima above noise threshold)
4. For each mobility peak, extract fragments only from scans within ±N scans
   of the peak center
5. Each mobility peak becomes a separate pseudo-DDA spectrum for Cascadia

Advantages:
- Fully independent of DIA-NN (parallel jobs preserved)
- Discovers precursors DIA-NN might miss
- No additional input files needed

Disadvantage:
- Mobilogram peak detection adds complexity
- Some false positive mobility peaks from noise
- Slightly more spectra to process (but Cascadia filters on score anyway)

**Recommendation: Implement Mode B first**, because it preserves the parallel
job architecture and serves the de novo discovery use case. Mode A can be added
later as an optional "targeted refinement" pass.

### 3.3 Mode B: Detailed algorithm

```
For each .d file:
  Read frame metadata from analysis.tdf
  Read DiaFrameMsMsWindows for isolation window definitions
  Build scan → 1/K0 calibration from GlobalMetadata

  For each MS2 frame (MsMsType != 0):
    For each isolation window active in this frame:
      scan_begin, scan_end = window's ScanNumBegin, ScanNumEnd
      
      # Step 1: Build mobilogram
      For each scan in [scan_begin, scan_end]:
        mobilogram[scan] = sum(intensities in this scan within isolation m/z range)
      
      # Step 2: Find mobility peaks
      peaks = find_peaks(mobilogram, prominence=noise_threshold)
      
      # Step 3: Extract filtered spectra
      For each mobility peak at scan_center:
        scan_lo = max(scan_begin, scan_center - half_width)
        scan_hi = min(scan_end, scan_center + half_width)
        
        # Collect fragments only from this narrow mobility slice
        mz_filtered = []
        intensity_filtered = []
        For scan in [scan_lo, scan_hi]:
          mz_filtered.extend(peaks_in_scan.mz)
          intensity_filtered.extend(peaks_in_scan.intensity)
        
        # Centroid / merge duplicate m/z values within tolerance
        spectrum = centroid_merge(mz_filtered, intensity_filtered, tol=0.01)
        
        # Record as a pseudo-DDA spectrum
        yield Spectrum(
          mz = spectrum.mz,
          intensity = spectrum.intensity,
          rt_sec = frame.time * 60,
          precursor_mz = isolation_center_mz,  # approximate
          isolation_width = narrow,  # effective width is much smaller now
          im = scan_to_1k0(scan_center),
          charge = 0,  # unknown; Cascadia handles charge prediction
        )
```

### 3.4 Key parameters

| Parameter | Default | Rationale |
|---|---|---|
| `mobility_half_width` | 3 scans (~0.015 Vs/cm²) | Typical peptide 1/K0 peak width is ~0.02–0.04 Vs/cm². 3 scans captures the core while excluding neighbors. |
| `mobilogram_prominence` | 100 (intensity units) | Filters out noise peaks. Adjustable based on instrument sensitivity. |
| `min_fragments` | 5 | Minimum fragment peaks in the filtered spectrum. Below this, skip — too little info for de novo. |
| `mz_merge_tol` | 0.01 Da | For centroiding fragments from adjacent scans. |
| `max_spectra_per_window` | 5 | Safety cap — a single wide window shouldn't produce dozens of spectra. |

### 3.5 Expected impact on spectrum quality

**Before (mobility-summed):**
- Isolation window: 400–440 m/z, scans 200–600 (wide 1/K0 range)
- All fragments from all precursors in this m/z × 1/K0 box summed together
- Typical chimeric contamination: 3–10 co-isolated precursors contributing fragments
- Cascadia sees a noisy mess → low confidence scores → fewer usable calls

**After (mobility-filtered):**
- Same isolation window, but fragments extracted from scans 350–356 only
- Contains fragments primarily from one precursor (or maybe 1–2 near-mobility neighbors)
- Chimeric contamination reduced by ~5–10×
- Cascadia sees a clean spectrum → higher confidence → more usable de novo calls

Conservative estimate: **2–3× more high-confidence de novo peptides** from the same
data from mobility filtering alone. With the IM-aware transformer (fine-tuned on
timsTOF data), potentially **3–5×** improvement as the model learns to use 1/K0
as a discriminating feature for sequence prediction.

---

## 4a. Integration with IM-Aware Transformer (Claude Code Fork)

### Background

A separate Claude Code session has already modified Cascadia's transformer architecture
to accept ion mobility as a 5th input feature. The key changes:

1. **ASF format**: Peaks go from 4 values `(m/z, intensity, RT_offset, MS_level)` to
   5 values `(m/z, intensity, RT_offset, MS_level, 1/K0)`. The parser auto-detects
   4 vs 5 column files.

2. **AugmentedPeakEncoder**: Gains a 4th encoding channel:
   `im_encoder = FloatEncoder(d_model, min_wavelength=0.001, max_wavelength=2.0)`
   using the same sinusoidal encoding as RT. Combiner linear layer expands from
   `Linear(3*d_model, d_model)` to `Linear(4*d_model, d_model)`.

3. **Checkpoint compatibility**: `use_ion_mobility` parameter (default `False`) on both
   `AugmentedPeakEncoder` and `AugmentedSpec2Pep`. `from_pretrained()` loads old 3-channel
   checkpoints into the new 4-channel model by zero-initializing the IM combiner columns —
   producing identical output when IM=0. `save_hyperparameters(ignore=['tokenizer'])`
   persists the flag in checkpoints.

4. **Config**: `use_ion_mobility: bool = False`, `min_im_wavelength: float = 0.001`,
   `max_im_wavelength: float = 2.0`. HDF5 storage includes `im_array` field.

### Why this matters for mobility-filtered extraction

Without the IM-aware model, mobility filtering helps by providing **cleaner spectra**
(less chimeric contamination). The model still only sees m/z/intensity/RT.

With the IM-aware model, mobility filtering provides **both**:
1. Cleaner spectra (same as above)
2. An accurate 1/K0 value that the transformer can use as a **discriminating feature**

The 1/K0 value constrains the sequence search space — ions at 1/K0 = 0.8 have different
likely sequences than ions at 1/K0 = 1.3 (because CCS correlates with peptide length,
charge, and amino acid composition). The model can learn these correlations during
fine-tuning.

### Three-phase deployment

| Phase | Model | Extraction | Expected gain |
|---|---|---|---|
| Phase A | Original checkpoint, `use_ion_mobility=False` | Mobility-filtered, 5-col ASF (IM ignored) | 2–3× from cleaner spectra |
| Phase B | Original checkpoint, `use_ion_mobility=True` | Mobility-filtered, 5-col ASF (IM zero-init) | ~2–3× (IM has minimal effect before fine-tuning) |
| Phase C | Fine-tuned checkpoint, `use_ion_mobility=True` | Mobility-filtered, 5-col ASF (IM fully used) | 3–5× (cleaner spectra + learned IM features) |

Phase A and B can be tested immediately. Phase C requires a fine-tuning run on
timsTOF data with mobility-filtered extraction — this is the dataset that would
teach the model what 1/K0 values mean for sequence prediction.

---

## 4. Implementation

### 4.1 Modified `bruker_augment.py`

This replaces the `_read_bruker_spectra()` function in the original spec. The new
version uses Frame-level access instead of pre-summed Spectra.

```python
"""
bruker_augment_mobility.py — Mobility-filtered Bruker .d loader for Cascadia.

Extracts mobility-resolved pseudo-DDA spectra from diaPASEF data by:
1. Reading raw frame data (preserving scan/1K0 dimension)
2. Building per-window mobilograms
3. Detecting mobility peaks
4. Extracting fragment spectra from narrow mobility slices

This produces dramatically cleaner input for de novo sequencing compared to
the mobility-summed approach.
"""

import numpy as np
import sqlite3
from pathlib import Path
from scipy.signal import find_peaks

try:
    import timsrust_pyo3 as timsrust
except ImportError:
    raise ImportError("pip install timsrust_pyo3")


class MobilityFilteredExtractor:
    """Extract mobility-filtered pseudo-DDA spectra from a timsTOF .d folder."""
    
    def __init__(
        self,
        d_folder: str,
        mobility_half_width: int = 3,       # scans
        mobilogram_prominence: float = 100,  # intensity units
        min_fragments: int = 5,
        mz_merge_tol: float = 0.01,          # Da
        max_spectra_per_window: int = 5,
    ):
        self.d_folder = Path(d_folder)
        self.mobility_half_width = mobility_half_width
        self.mobilogram_prominence = mobilogram_prominence
        self.min_fragments = min_fragments
        self.mz_merge_tol = mz_merge_tol
        self.max_spectra_per_window = max_spectra_per_window
        
        # Load calibration and window definitions from analysis.tdf
        self.tdf_path = self.d_folder / "analysis.tdf"
        if not self.tdf_path.exists():
            raise FileNotFoundError(f"No analysis.tdf in {d_folder}")
        
        self._load_metadata()
    
    def _load_metadata(self):
        """Load scan↔1/K0 calibration and diaPASEF window definitions."""
        con = sqlite3.connect(str(self.tdf_path))
        
        # Scan-to-1/K0 calibration
        # GlobalMetadata stores calibration coefficients
        meta = dict(con.execute(
            "SELECT Key, Value FROM GlobalMetadata"
        ).fetchall())
        
        # The exact calibration varies by TDF version, but the key fields are:
        # OneOverK0Begin, OneOverK0End, and the number of scans
        # For now, store raw metadata — the calibration function is:
        #   1/K0 = OneOverK0End - (scan / NumScans) * (OneOverK0End - OneOverK0Begin)
        # (Higher scan numbers = lower 1/K0 = smaller CCS)
        self.one_over_k0_begin = float(meta.get('OneOverK0Begin', 0.6))
        self.one_over_k0_end = float(meta.get('OneOverK0End', 1.6))
        
        # Frame info
        self.frames_df = con.execute("""
            SELECT Id, Time, MsMsType, NumScans
            FROM Frames
            ORDER BY Id
        """).fetchall()
        
        # diaPASEF window definitions
        # This table maps frames to their isolation windows with mobility bounds
        try:
            self.dia_windows = con.execute("""
                SELECT Frame, WindowGroup, ScanNumBegin, ScanNumEnd,
                       IsolationMz, IsolationWidth
                FROM DiaFrameMsMsWindows
                ORDER BY Frame, WindowGroup
            """).fetchall()
        except sqlite3.OperationalError:
            # Fallback for older TDF versions or ddaPASEF
            # Try PasefFrameMsMsInfo for DDA data
            self.dia_windows = con.execute("""
                SELECT Frame, 0 as WindowGroup, ScanNumber as ScanNumBegin,
                       ScanNumber as ScanNumEnd, IsolationMz, IsolationWidth
                FROM PasefFrameMsMsInfo
                ORDER BY Frame
            """).fetchall()
        
        con.close()
    
    def scan_to_1k0(self, scan_num: int, num_scans: int) -> float:
        """Convert scan number to 1/K0 value."""
        # Linear interpolation between calibration endpoints
        # Scan 0 = highest 1/K0 (largest ions), scan N = lowest 1/K0 (smallest)
        if num_scans <= 1:
            return (self.one_over_k0_begin + self.one_over_k0_end) / 2
        frac = scan_num / (num_scans - 1)
        return self.one_over_k0_end - frac * (self.one_over_k0_end - self.one_over_k0_begin)
    
    def _build_mobilogram(
        self, frame_data, scan_begin: int, scan_end: int,
        iso_mz: float, iso_width: float
    ) -> np.ndarray:
        """
        Build an intensity-vs-scan mobilogram for a given isolation window.
        
        frame_data: timsrust Frame object with per-scan peak lists
        scan_begin/end: mobility boundaries of the isolation window
        iso_mz/iso_width: m/z boundaries (for filtering fragment m/z range)
        
        Returns: 1D array of summed intensity per scan in [scan_begin, scan_end]
        """
        n_scans = scan_end - scan_begin + 1
        mobilogram = np.zeros(n_scans, dtype=np.float64)
        
        for scan_offset in range(n_scans):
            scan_idx = scan_begin + scan_offset
            # Get peaks for this specific scan from the frame
            # timsrust Frame API: frame.scan_offsets gives boundaries,
            # frame.mz_values and frame.intensities are flat arrays
            # indexed by scan via offsets
            scan_start = frame_data.scan_offsets[scan_idx]
            scan_end_idx = frame_data.scan_offsets[scan_idx + 1]
            
            if scan_start >= scan_end_idx:
                continue
            
            intensities = frame_data.intensities[scan_start:scan_end_idx]
            mobilogram[scan_offset] = np.sum(intensities)
        
        return mobilogram
    
    def _extract_scan_spectrum(
        self, frame_data, scan_lo: int, scan_hi: int
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Extract and merge fragment peaks from a narrow scan range.
        
        Returns: (mz_array, intensity_array) after centroid merging.
        """
        all_mz = []
        all_int = []
        
        for scan_idx in range(scan_lo, scan_hi + 1):
            scan_start = frame_data.scan_offsets[scan_idx]
            scan_end_idx = frame_data.scan_offsets[scan_idx + 1]
            
            if scan_start >= scan_end_idx:
                continue
            
            all_mz.append(frame_data.mz_values[scan_start:scan_end_idx])
            all_int.append(frame_data.intensities[scan_start:scan_end_idx])
        
        if not all_mz:
            return np.array([]), np.array([])
        
        mz = np.concatenate(all_mz)
        intensity = np.concatenate(all_int)
        
        if len(mz) == 0:
            return mz, intensity
        
        # Centroid merge: sort by m/z, merge peaks within tolerance
        order = np.argsort(mz)
        mz = mz[order]
        intensity = intensity[order]
        
        merged_mz = []
        merged_int = []
        i = 0
        while i < len(mz):
            j = i + 1
            while j < len(mz) and (mz[j] - mz[i]) < self.mz_merge_tol:
                j += 1
            # Intensity-weighted average m/z
            cluster_mz = mz[i:j]
            cluster_int = intensity[i:j]
            total_int = np.sum(cluster_int)
            avg_mz = np.average(cluster_mz, weights=cluster_int)
            merged_mz.append(avg_mz)
            merged_int.append(total_int)
            i = j
        
        return np.array(merged_mz, dtype=np.float32), np.array(merged_int, dtype=np.float32)
    
    def extract_filtered_spectra(self) -> list[dict]:
        """
        Main extraction method.
        
        Returns a list of pseudo-DDA spectrum dicts suitable for Cascadia:
        {
            'mz': np.ndarray,
            'intensity': np.ndarray,
            'rt_sec': float,
            'precursor_mz': float,
            'isolation_width': float,
            'im': float,          # 1/K0 at mobility peak center
            'charge': int,        # 0 = unknown
            'scan_idx': int,      # for SSL output
            'mobility_filtered': True,
        }
        """
        reader = timsrust.FileReader(str(self.d_folder))
        
        # Read ALL frames (not spectra — we need scan-level resolution)
        frames = reader.read_all_frames()
        
        # Build frame index: frame_id → frame object
        frame_index = {}
        for frame in frames:
            frame_index[frame.frame_id] = frame
        
        # Build frame metadata index
        frame_meta = {}
        for (fid, time, msms_type, num_scans) in self.frames_df:
            frame_meta[fid] = {
                'time': time,           # seconds
                'msms_type': msms_type,
                'num_scans': num_scans,
            }
        
        spectra = []
        
        for (frame_id, window_group, scan_begin, scan_end,
             iso_mz, iso_width) in self.dia_windows:
            
            if frame_id not in frame_index:
                continue
            
            frame = frame_index[frame_id]
            meta = frame_meta.get(frame_id, {})
            rt_sec = meta.get('time', 0)  # analysis.tdf Time is in seconds
            num_scans = meta.get('num_scans', 1000)
            
            # Step 1: Build mobilogram for this window
            mobilogram = self._build_mobilogram(
                frame, scan_begin, scan_end, iso_mz, iso_width
            )
            
            if np.max(mobilogram) == 0:
                continue
            
            # Step 2: Find mobility peaks
            peak_indices, properties = find_peaks(
                mobilogram,
                prominence=self.mobilogram_prominence,
                distance=self.mobility_half_width,  # minimum separation
            )
            
            if len(peak_indices) == 0:
                # No clear mobility peaks — fall back to summing the whole window
                # (this handles edge cases like very low-abundance windows)
                mz, intensity = self._extract_scan_spectrum(
                    frame, scan_begin, scan_end
                )
                if len(mz) >= self.min_fragments:
                    spectra.append({
                        'mz': mz,
                        'intensity': intensity,
                        'rt_sec': rt_sec,
                        'precursor_mz': iso_mz,
                        'isolation_width': iso_width,
                        'im': self.scan_to_1k0(
                            (scan_begin + scan_end) // 2, num_scans
                        ),
                        'charge': 0,
                        'scan_idx': len(spectra),
                        'mobility_filtered': False,
                    })
                continue
            
            # Cap the number of spectra per window
            if len(peak_indices) > self.max_spectra_per_window:
                # Keep the most prominent peaks
                top_idx = np.argsort(properties['prominences'])[::-1]
                peak_indices = peak_indices[top_idx[:self.max_spectra_per_window]]
            
            # Step 3: Extract filtered spectrum for each mobility peak
            for peak_idx in peak_indices:
                scan_center = scan_begin + peak_idx
                scan_lo = max(scan_begin, scan_center - self.mobility_half_width)
                scan_hi = min(scan_end, scan_center + self.mobility_half_width)
                
                mz, intensity = self._extract_scan_spectrum(
                    frame, scan_lo, scan_hi
                )
                
                if len(mz) < self.min_fragments:
                    continue
                
                im_value = self.scan_to_1k0(scan_center, num_scans)
                
                spectra.append({
                    'mz': mz,
                    'intensity': intensity,
                    'rt_sec': rt_sec,
                    'precursor_mz': iso_mz,
                    'isolation_width': iso_width,
                    'im': im_value,
                    'charge': 0,
                    'scan_idx': len(spectra),
                    'mobility_filtered': True,
                })
        
        return spectra


def augment_spectra_bruker_mobility(
    d_folder: str,
    output_asf_path: str,
    augmentation_width: int = 5,
    mobility_half_width: int = 3,
    mobilogram_prominence: float = 100,
    min_fragments: int = 5,
) -> int:
    """
    Main entry point. Replaces augment_spectra_bruker() from the original spec.
    
    Reads a Bruker .d folder with mobility-filtered extraction and writes an
    augmented spectrum file (.asf) for Cascadia's transformer.
    
    Returns: number of spectra written.
    """
    extractor = MobilityFilteredExtractor(
        d_folder,
        mobility_half_width=mobility_half_width,
        mobilogram_prominence=mobilogram_prominence,
        min_fragments=min_fragments,
    )
    
    spectra = extractor.extract_filtered_spectra()
    
    # Log stats
    n_filtered = sum(1 for s in spectra if s.get('mobility_filtered', False))
    n_fallback = len(spectra) - n_filtered
    print(f"[MobilityFilter] {len(spectra)} total spectra: "
          f"{n_filtered} mobility-filtered, {n_fallback} fallback (summed)")
    
    # Convert to Cascadia's 5-column ASF format (IM-aware model)
    #
    # The modified Cascadia transformer (AugmentedPeakEncoder) expects ASF peaks
    # with 5 values per peak: (m/z, intensity, RT_offset, MS_level, 1/K0)
    #
    # When use_ion_mobility=True, the model uses the 5th column via a dedicated
    # im_encoder (FloatEncoder with sinusoidal encoding, wavelength 0.001–2.0).
    # The combiner linear layer is 4*d_model → d_model (vs 3*d_model for 4-col).
    #
    # When IM is unavailable (e.g., Orbitrap data), the 5th value is 0.0 and the
    # model produces identical output to 4-column mode via zero-initialized
    # combiner weights.
    #
    # For mobility-filtered spectra, we pass the actual 1/K0 at the mobility peak
    # center. For fallback (summed) spectra, we pass the midpoint 1/K0 of the
    # isolation window (less precise, but better than 0.0).
    
    records = []
    for spec in spectra:
        if len(spec['mz']) < min_fragments:
            continue
        n_peaks = len(spec['mz'])
        records.append({
            'mz': spec['mz'],
            'intensity': spec['intensity'],
            'rt': spec['rt_sec'],
            'precursor_mz': spec['precursor_mz'],
            'charge': spec['charge'],
            'im': spec['im'],  # 1/K0 in Vs/cm², typically 0.6–1.6
        })
    
    # Write as 5-column ASF (auto-detected by Cascadia's parser)
    # The parser checks column count: 4 = legacy, 5 = IM-aware
    import pickle
    with open(output_asf_path, "wb") as f:
        pickle.dump(records, f, protocol=4)
    
    return len(records)
```

### 4.2 Routing change in `cascadia.py`

Update the routing in CASCADIA_DENOVO_SPEC.md Phase 1 to use the mobility-filtered
extractor by default for Bruker data:

```python
if input_path.suffix == ".d" or (input_path.is_dir() and
        (input_path / "analysis.tdf").exists()):
    # Bruker timsTOF .d directory — use mobility-filtered native loader
    from cascadia.bruker_augment_mobility import augment_spectra_bruker_mobility
    n_spectra = augment_spectra_bruker_mobility(
        str(input_path),
        str(asf_output_path),
        augmentation_width=config.augmentation_width,
        mobility_half_width=config.mobility_half_width,   # new config param
        mobilogram_prominence=config.mobilogram_prominence,  # new config param
    )
```

### 4.3 New Cascadia config parameters

Add to Cascadia's config:

```yaml
# Mobility-filtered extraction (Bruker timsTOF only)
mobility_half_width: 3        # Scans on each side of mobility peak (default: 3)
mobilogram_prominence: 100    # Min prominence for mobility peak detection (default: 100)
min_fragments: 5              # Min fragment peaks per spectrum (default: 5)
mobility_filter: true         # Set false to revert to summed extraction

# IM-aware transformer (from Claude Code fork)
use_ion_mobility: true        # Enable 4-channel AugmentedPeakEncoder (im_encoder)
min_im_wavelength: 0.001      # FloatEncoder sinusoidal encoding range
max_im_wavelength: 2.0        # Covers 1/K0 range of 0.6–1.6 Vs/cm²
```

### 4.4 How the pieces connect (end-to-end data flow)

This is the critical integration point between three separate development efforts:

```
bruker_augment_mobility.py          Cascadia model (Claude Code fork)
──────────────────────────          ─────────────────────────────────
                                    
1. Read .d frames                   
2. Build mobilograms                
3. Detect mobility peaks            
4. Extract filtered spectra         
5. Write 5-column ASF:              
   (m/z, intensity, RT, MS_level,   → Parser auto-detects 5 columns
    1/K0)                           → AugmentedPeakEncoder uses im_encoder
                                    → FloatEncoder(d_model, 0.001, 2.0)
   Each peak's 1/K0 = actual        → Sinusoidal encoding of 1/K0 value
   mobility at extraction center    → Combiner: Linear(4*d_model, d_model)
                                    → Transformer predicts sequence
                                    → SSL output (unchanged)
```

The 1/K0 value written by the extractor is the **actual ion mobility at the
mobilogram peak center** — not an average, not a window midpoint. This is
exactly the signal the IM-aware transformer was designed to consume.

For fallback (summed) spectra where no clear mobility peak was detected, the
1/K0 is the midpoint of the isolation window's scan range. This is less precise
but still informative — the model will weight it lower via the learned combiner.

### 4.5 Checkpoint compatibility

The modified Cascadia model supports both old and new checkpoints:

- **Old 3-channel checkpoint + `use_ion_mobility=False`**: Original behavior.
  4-column ASF input, 3-channel encoder. Ignores any 5th column if present.
  
- **Old 3-channel checkpoint + `use_ion_mobility=True`**: `from_pretrained()`
  zero-initializes the IM combiner columns. The model produces identical output
  when IM=0.0, and begins using IM signal immediately (even before fine-tuning)
  via the untrained im_encoder.
  
- **New 4-channel checkpoint + `use_ion_mobility=True`**: Full IM-aware model.
  Best performance with mobility-filtered extraction.

**For initial testing**, you can use the zero-initialized path (old checkpoint +
`use_ion_mobility=True`) to verify the pipeline works end-to-end before committing
to a fine-tuning run. The IM signal will have minimal impact until fine-tuned, but
the mobility-filtered extraction alone (cleaner spectra) should still improve results.

---

## 5. Validation Plan

### 5.1 Quick sanity check (~30 min)

Before running full Cascadia, validate the extraction itself:

```python
"""
test_mobility_extraction.py — Validate mobility-filtered extraction on a real .d file.
Run on Hive in the cascadia conda env.
"""
from bruker_augment_mobility import MobilityFilteredExtractor
import numpy as np

d_file = "/path/to/test_sample.d"  # one of your timsTOF DIA files
extractor = MobilityFilteredExtractor(d_file)
spectra = extractor.extract_filtered_spectra()

# Basic stats
n_total = len(spectra)
n_filtered = sum(1 for s in spectra if s['mobility_filtered'])
n_fallback = n_total - n_filtered

frag_counts = [len(s['mz']) for s in spectra]
im_values = [s['im'] for s in spectra if s['mobility_filtered']]

print(f"Total spectra: {n_total}")
print(f"Mobility-filtered: {n_filtered} ({100*n_filtered/n_total:.1f}%)")
print(f"Fallback (summed): {n_fallback}")
print(f"Fragment count: median={np.median(frag_counts):.0f}, "
      f"mean={np.mean(frag_counts):.1f}, range={min(frag_counts)}-{max(frag_counts)}")
print(f"1/K0 range: {min(im_values):.3f} - {max(im_values):.3f}")
print(f"Memory: ~{sum(s['mz'].nbytes + s['intensity'].nbytes for s in spectra) / 1e6:.1f} MB")

# Sanity checks
assert n_total > 10000, f"Too few spectra ({n_total}) — extraction may be broken"
assert n_filtered > n_total * 0.5, f"Too few filtered spectra — mobilogram peak detection may be too strict"
assert np.median(frag_counts) >= 5, "Median fragment count too low"
assert min(im_values) > 0.5 and max(im_values) < 2.0, "1/K0 values out of expected range"
print("All sanity checks passed!")
```

### 5.2 Head-to-head comparison (~2 hr)

Run Cascadia on the same `.d` file three ways:

| Run | Extraction | Model | Config |
|---|---|---|---|
| A | Mobility-summed (original) | Original checkpoint | `use_ion_mobility: false` |
| B | Mobility-filtered | Original checkpoint | `use_ion_mobility: false` |
| C | Mobility-filtered | Original checkpoint | `use_ion_mobility: true` (zero-init) |

Run B vs A isolates the extraction improvement. Run C vs B shows whether the
zero-initialized IM encoder helps even before fine-tuning.

Compare:

```python
"""
compare_extraction_modes.py — Compare summed vs filtered vs IM-aware Cascadia output.
"""
import pandas as pd

runs = {
    "A_summed":       "output_summed/sample.ssl",
    "B_filtered":     "output_filtered/sample.ssl",
    "C_filtered_im":  "output_filtered_im/sample.ssl",
}

ssls = {name: pd.read_csv(path, sep="\t") for name, path in runs.items()}

# 1. Total yield at different score thresholds
for threshold in [0.5, 0.7, 0.8, 0.9]:
    counts = {name: (ssl['score'] >= threshold).sum() for name, ssl in ssls.items()}
    baseline = max(counts["A_summed"], 1)
    print(f"Score >= {threshold}:")
    for name, n in counts.items():
        print(f"  {name}: {n} ({n/baseline:.2f}x vs summed)")

# 2. Score distribution
print(f"\nScore distribution:")
for name, ssl in ssls.items():
    print(f"  {name}: median={ssl['score'].median():.3f}, "
          f"mean={ssl['score'].mean():.3f}")

# 3. Sequence overlap at score >= 0.8
seqs = {name: set(ssl.loc[ssl['score'] >= 0.8, 'sequence']) for name, ssl in ssls.items()}

print(f"\nAt score >= 0.8:")
for name, s in seqs.items():
    print(f"  {name}: {len(s)} unique sequences")

overlap_ab = seqs["A_summed"] & seqs["B_filtered"]
filtered_gain = seqs["B_filtered"] - seqs["A_summed"]
im_gain = seqs["C_filtered_im"] - seqs["B_filtered"]
print(f"\n  A∩B overlap: {len(overlap_ab)}")
print(f"  B-only (filtering gain): +{len(filtered_gain)}")
print(f"  C-only (IM model gain): +{len(im_gain)}")
```

**Expected results:**
- 1.5–3× more high-confidence (≥0.8) de novo calls from mobility-filtered extraction
- Higher median score across all predictions
- Most summed-mode hits are a subset of filtered-mode hits (filtered finds everything
  summed finds, plus more)
- The filtered-only sequences should cross-reference well against DIA-NN's report.tsv
  (confirming they're real peptides, not noise-derived false positives)

### 5.3 Cross-validation against DIA-NN (~30 min)

The ultimate test: do the mobility-filtered de novo sequences match DIA-NN's database
search results better?

```python
"""
crossref_diann.py — Compare de novo sequences to DIA-NN identifications.
Uses the classify_denovo_peptides() logic from CASCADIA_DENOVO_SPEC Phase 7.
"""
# Load DIA-NN report
diann_peptides = set(diann_report['Stripped.Sequence'].unique())

# Normalize I/L
def normalize_il(seq):
    return seq.replace('I', 'L').replace('[+', '[').upper()

diann_normalized = {normalize_il(p) for p in diann_peptides}

for label, ssl_path in [("summed", ssl_summed_path), ("filtered", ssl_filtered_path)]:
    ssl = pd.read_csv(ssl_path, sep="\t")
    ssl_good = ssl[ssl['score'] >= 0.8]
    
    denovo_seqs = {normalize_il(s.split('[')[0]) for s in ssl_good['sequence']}
    # Simplified — strip mods for matching
    
    confirmed = denovo_seqs & diann_normalized
    novel = denovo_seqs - diann_normalized
    
    print(f"{label}: {len(confirmed)} confirmed, {len(novel)} novel "
          f"({100*len(confirmed)/max(len(denovo_seqs),1):.1f}% confirmation rate)")
```

**Expected:** Mobility-filtered extraction should have a **higher confirmation rate**
(more of its de novo sequences match DIA-NN), because cleaner spectra → more accurate
sequence predictions → more matches to the database search ground truth.

---

## 6. Risks and Mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| timsrust Frame API doesn't expose per-scan peak lists | Medium | Blocking | Fall back to direct SQLite reading of the binary TDF blob (opentimsr approach). OR use `rustims` (imspy) which explicitly supports frame-level scan access. |
| Mobilogram peak detection is too sensitive/insensitive | Medium | Quality | Make prominence threshold a tunable parameter. Start conservative (high threshold), lower if too few spectra extracted. |
| More spectra → longer Cascadia runtime | Low | Minor | Typically 1.5–2× more spectra, but they're individually smaller. Net GPU time increase ~30%. Acceptable. |
| Cascadia IM-aware model not yet fine-tuned on timsTOF data | Medium | Moderate | Zero-initialized combiner means the model works immediately (IM has no effect until fine-tuned). Mobility-filtered extraction still helps via cleaner spectra alone. Fine-tuning on a timsTOF dataset is the next step for full benefit. |
| 1/K0 calibration inaccuracy (timsrust vs Bruker SDK) | Low | Negligible | Calibration only needs to be accurate enough to separate mobility peaks. ±0.01 Vs/cm² is plenty. We're not doing CCS prediction. |
| Scan numbering/offset mismatch between analysis.tdf and timsrust | Medium | Blocking | Validate scan offsets against known precursors from DIA-NN report (which includes IM values). |

### Biggest risk: timsrust Frame API coverage

The timsrust GitHub README confirms that Frames are a first-class data type. But the
exact Python API for accessing per-scan data within a frame (scan_offsets, per-scan
mz/intensity arrays) needs verification on Hive. The recon step should include:

```python
import timsrust_pyo3 as timsrust

reader = timsrust.FileReader("/path/to/test.d")
frames = reader.read_all_frames()

frame = frames[100]  # pick an MS2 frame
print(f"Frame attributes: {dir(frame)}")
print(f"Frame ID: {frame.frame_id if hasattr(frame, 'frame_id') else 'N/A'}")
print(f"Has scan_offsets: {hasattr(frame, 'scan_offsets')}")
print(f"Has tof_indices: {hasattr(frame, 'tof_indices')}")

# If scan_offsets exists, we're good.
# If not, check for alternative scan-level access patterns.
```

**Fallback if timsrust doesn't expose scan-level access:**

Use `opentimsr` (Python) or direct SQLite + binary blob reading. The TDF binary
format stores compressed scan data in the `analysis.tdf_bin` file, indexed by the
`TimsId` column in the `Frames` table. Libraries like OpenTIMS/TimsPy handle the
decompression. This is more complex but definitely works — it's how AlphaTims and
DIA-NN itself read the data.

Alternative fallback: Use `rustims` (`imspy-core`) instead of `timsrust_pyo3`.
The rustims framework explicitly supports frame-level scan access with its
`TimsFrame` class and has Python bindings via PyO3.

---

## 7. Integration with Existing Specs

### CASCADIA_DENOVO_SPEC.md

This addendum **replaces Phase 1 Section 5** ("Bruker Native Loader"). Everything
downstream of spectrum extraction is unchanged:
- SSL output format: unchanged
- SLURM job scripts: unchanged (just swap the import)
- R-side integration (Phase 3+): unchanged
- DIA-NN cross-referencing: unchanged (and will benefit from higher-quality de novo calls)

### EVIDENCE_DEPTH_SPEC.md

Mobility-filtered de novo confirmation is a stronger evidence signal than summed-mode
de novo confirmation. Consider adding a sub-tier:
- **Strong+**: DIA-NN identification + mobility-filtered de novo confirmation
- **Strong**: DIA-NN identification + summed-mode de novo confirmation

### XIC Viewer

The mobilogram data extracted during this process could feed a "Mobilogram" pane
in the XIC Viewer — showing the 1/K0 profile for a selected precursor alongside
its fragment XICs. This would give DE-LIMP a Spectronaut-like mobilogram visualization
essentially for free, since the extraction code already computes it.

---

## 8. Claude Code Implementation Plan

### Step 0: Recon — verify timsrust Frame API (~15 min)

```
Run on Hive in cascadia conda env:
python -c "
import timsrust_pyo3 as tr
r = tr.FileReader('/path/to/test.d')
frames = r.read_all_frames()
f = frames[100]
print(dir(f))
print(type(f))
# Check for scan_offsets, mz_values, intensities, tof_indices
"
```

If `scan_offsets` or equivalent exists → proceed with this spec.
If not → evaluate `rustims` / `imspy-core` as alternative, or fall back to
OpenTIMS/TimsPy.

### Step 1: Implement MobilityFilteredExtractor (~1 hr)

Tell Claude Code:
> "Create `cascadia/bruker_augment_mobility.py` from the spec in
> CASCADIA_MOBILITY_FILTER_ADDENDUM.md section 4.1. Read the existing
> `bruker_augment.py` first to understand the output format. The key change
> is reading Frames instead of Spectra and doing per-scan mobility filtering.
> Adapt the Frame API calls based on what we found in Step 0 recon."

### Step 2: Run validation test (~30 min)

Tell Claude Code:
> "Run test_mobility_extraction.py from section 5.1 on [test .d file].
> Report the stats. If sanity checks fail, adjust mobilogram_prominence
> or mobility_half_width parameters."

### Step 3: Head-to-head Cascadia comparison (~2 hr)

Tell Claude Code:
> "Submit three Cascadia jobs on the same .d file:
> (A) original bruker_augment.py + use_ion_mobility=false
> (B) bruker_augment_mobility.py + use_ion_mobility=false
> (C) bruker_augment_mobility.py + use_ion_mobility=true
> Use the sbatch script from CASCADIA_DENOVO_SPEC Phase 6.
> After all complete, run compare_extraction_modes.py from section 5.2."

### Step 4: Wire up the routing (~15 min)

Tell Claude Code:
> "Update cascadia.py routing to use bruker_augment_mobility by default
> for .d input. Add mobility_half_width and mobilogram_prominence to the
> Cascadia config. Add a --mobility-filter / --no-mobility-filter CLI flag.
> Ensure the ASF writer outputs 5-column format with the im field."

### Step 5: Fine-tuning dataset preparation (future, after validation)

Once Steps 2–4 confirm the pipeline works, prepare a fine-tuning dataset:

Tell Claude Code:
> "Extract mobility-filtered spectra from [N timsTOF DIA files] using
> bruker_augment_mobility.py. Cross-reference against DIA-NN report.tsv
> to get ground-truth peptide labels for each spectrum. Write as HDF5
> with im_array field. This will be the fine-tuning set for the IM-aware
> transformer."

This is Phase C — teaching the model what 1/K0 means for sequence prediction.

---

*Spec version 1.0 — Brett Phinney / UC Davis Proteomics Core*
*Addendum to CASCADIA_DENOVO_SPEC.md*
