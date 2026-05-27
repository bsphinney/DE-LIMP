# Cascadia De Novo Sequencing — Bruker timsTOF Integration

## Overview

This document describes how DE-LIMP integrates Cascadia (a transformer-based de novo peptide sequencer) with native Bruker timsTOF `.d` files, eliminating the need for mzML conversion. The implementation replaces Cascadia's default mzML-based augment.py with a timsrust_pyo3-based pipeline that reads `.d` files directly.

**Authors**: Brett Phinney, UC Davis Proteomics Core  
**Date**: April 2026  
**Status**: In development (feature/cascadia-denovo branch)

## Architecture

```
Standard Cascadia Pipeline:
  .d file → msconvert → .mzML → augment.py → .asf → model → .ssl

Native Bruker Pipeline (this work):
  .d file → bruker_patch.py (timsrust_pyo3) → .asf → model → .ssl
```

The native pipeline skips mzML conversion (45–90 min per file) and reads directly from the Bruker TDF format (~2–5 min per file).

## Key Components

### 1. `bruker_patch.py` — Drop-in replacement for augment.py

**Location**: Install into `cascadia/bruker_patch.py` in the Cascadia package directory.

**Functions**:
- `get_centers_bruker(d_folder)` — Replaces `get_centers()`. Reads MS2 DIA windows from `.d` via `timsrust.read_all_spectra()`. Returns the same `(f_to_mzrt_to_pep, max_mz, window_size, cycle_time)` tuple.
- `extract_spectra_bruker(d_folder, ...)` — Replaces `extract_spectra()`. Reads MS2 fragments AND MS1 survey frames, applies the same normalization as augment.py, returns `prec_to_spec` dict compatible with `write_asf()`.

**Critical implementation details**:

#### MS1 Frame Reading
timsrust's `read_all_spectra()` returns only MS2 (DIA window) spectra — no MS1. MS1 survey scans must be read separately via `FrameReader.read_ms1_frames()`, which returns `Frame` objects with raw `tof_indices` (not m/z values).

#### TOF-to-m/z Conversion (the Sage formula)
Frame objects have `tof_indices` but NOT `mz_values`. The conversion uses the same formula as Sage/timsrust internally (from `src/converters.rs`, `Tof2MzConverter`):

```python
import numpy as np
import sqlite3

def tof_to_mz(tof_indices, d_folder):
    """Convert raw TOF indices to m/z using the timsrust/Sage formula."""
    # Read calibration from GlobalMetadata table
    tdf_path = os.path.join(d_folder, "analysis.tdf")
    conn = sqlite3.connect(f"file:{tdf_path}?mode=ro", uri=True)
    meta = dict(conn.execute("SELECT Key, Value FROM GlobalMetadata").fetchall())
    conn.close()
    
    mz_min = float(meta['MzAcqRangeLower'])   # e.g., 99.994
    mz_max = float(meta['MzAcqRangeUpper'])   # e.g., 1700.0
    dns = int(meta['DigitizerNumSamples'])      # e.g., 631025
    
    # Sage formula (timsrust src/converters.rs):
    #   tof_intercept = sqrt(mz_min)
    #   tof_slope = (sqrt(mz_max) - sqrt(mz_min)) / DigitizerNumSamples
    #   mz = (tof_intercept + tof_slope * tof_index) ^ 2
    tof = np.array(tof_indices, dtype=np.float64)
    tof_intercept = np.sqrt(mz_min)
    tof_slope = (np.sqrt(mz_max) - tof_intercept) / dns
    return (tof_intercept + tof_slope * tof) ** 2
```

**Why not use `Metadata.mz_converter`?** The timsrust `Metadata()` class has a `mz_converter` attribute that does this internally, but it fails with `OSError: unable to open database file` on network-mounted storage (quobyte, NFS, SMB). Reading the 3 values from `GlobalMetadata` via SQLite with `?mode=ro` works reliably on all filesystems.

**WARNING**: A naive quadratic approximation (`mz = mz_min + (mz_max - mz_min) * (tof/dns)^2`) is wrong by up to 155 Da in the mid-range. The Sage formula is a linear interpolation in sqrt(m/z) space, not a quadratic in m/z space. Use the exact formula above.

#### Intensity Normalization (must match augment.py exactly)

| Scan Type | augment.py normalization | Description |
|-----------|------------------------|-------------|
| MS1 | `sqrt(sqrt(intensity))`, then `/ max` | Fourth root, then unit normalize |
| MS2 | `sqrt(intensity)`, then `/ max` | Square root, then unit normalize |

Both also: sort by intensity (top 150), re-sort by m/z.

### 2. `run_bruker_e2e.py` — End-to-end inference script

**Location**: `/quobyte/proteomics-grp/de-limp/cascadia/run_bruker_e2e.py`

Usage:
```bash
python3 run_bruker_e2e.py <d_folder> <model_checkpoint> <output_dir> [score_threshold] [batch_size]
```

Steps:
1. `get_centers_bruker()` — scan MS2 windows, estimate cycle time
2. `extract_spectra_bruker()` — build augmented spectrum file with MS1+MS2 context
3. Load ASF into Cascadia's `AnnotatedSpectrumDataset`
4. Run GPU inference with `AugmentedSpec2Pep` model
5. Write SSL output via Cascadia's `write_results()`

**Important parameters**:
- `time_width = (scan_width + 1) * cycle_time` — must match augment.py
- `score_threshold = 0.5` for testing, `0.8` for production (Cascadia default)
- `batch_size = 64` on A100 80GB

### 3. `prepare_training_data.py` — Training data from DIA-NN results

**Location**: `/quobyte/proteomics-grp/de-limp/cascadia/prepare_training_data.py`

Creates labeled ASF files for Cascadia fine-tuning using DIA-NN search results as ground truth labels.

Usage:
```bash
python3 prepare_training_data.py \
    --report /path/to/report.parquet \
    --raw-dir /path/to/raw_files/ \
    --output /path/to/training.asf \
    --q-threshold 0.01
```

**How it works**:
1. Load DIA-NN `report.parquet` — filter to Q.Value ≤ 0.01
2. For each `.d` file in the report:
   a. Read MS2 spectra via `timsrust.read_all_spectra()` — normalize, index by RT/mz
   b. Read MS1 frames via `FrameReader.read_ms1_frames()` — convert tof to m/z
   c. Match DIA-NN peptide-spectrum matches to raw spectra by RT and precursor m/z
   d. Write labeled ASF (peptide sequence in SEQ field instead of dummy 'K' * 30)
3. If no MS1 frame matches a precursor, generate a synthetic isotope pattern as fallback

**DIA-NN 2.0 note**: The `Run` column (not `File.Name`) contains bare filenames without `.d` extension. The script appends `.d` and joins with `--raw-dir`.

## HPC Environment Setup

### Conda Environment
```bash
# Environment: /quobyte/proteomics-grp/envs/cascadia5/
conda activate /quobyte/proteomics-grp/envs/cascadia5

# Key packages:
#   cascadia 0.0.7
#   pytorch 2.0.1 (CUDA)
#   timsrust_pyo3 (for native .d reading)
#   pyarrow (for DIA-NN report.parquet)
```

### Model Checkpoint
```
/quobyte/proteomics-grp/de-limp/cascadia/models/cascadia.ckpt  (558 MB)
```
This is the pre-trained Cascadia model from the original paper. Fine-tuned models will be saved alongside it.

### SLURM Job Templates

**Inference (single file)**:
```bash
#SBATCH --partition=gpu-a100
#SBATCH --account=genome-center-grp
#SBATCH --qos=genome-center-grp-gpu-a100-qos
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
```

**Training (fine-tuning)**:
```bash
#SBATCH --partition=gpu-a100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
```

## Gotchas & Troubleshooting

| Problem | Solution |
|---------|----------|
| `lxml.etree.XMLSyntaxError: Resource limit exceeded` on large mzML | Add `huge_tree=True` to `mzml.read()` calls in `cascadia/augment.py` |
| `Metadata()` fails with `OSError: unable to open database file` | Network storage (quobyte/NFS) issue. Read GlobalMetadata via SQLite `?mode=ro` instead |
| `read_all_spectra()` returns 0 MS1 spectra | By design — timsrust only returns MS2 DIA spectra. Use `FrameReader.read_ms1_frames()` for MS1 |
| `Frame` objects have `tof_indices` not `mz_values` | Apply Sage tof-to-mz formula (see above). DO NOT use naive quadratic approximation |
| `write_asf()` produces 0-byte file | Missing MS1 context — `prec_to_spec` entries without `ms1_scans` are skipped |
| `score >= 0.5` TypeError in SSL writing | Score field read as string. Use Cascadia's `write_results()` which handles conversion internally |
| Training produces 0 labeled spectra | RT unit mismatch (minutes vs seconds) OR missing MS1 frames. Check both |
| `pyarrow` not found | `pip install pyarrow` in the cascadia5 env |
| Cascadia `train` CLI rejects `--batch-size` | Use short flags: `-b 64 -e 10 -lr 1e-5 -c 4` (underscores in long form: `--batch_size`) |
| Older Berger DIA data has isolation window centers as precursor m/z | Standard DIA (not diaPASEF) — precursor m/z values like 412.5, 612.5 are window centers, not real precursors. Training data from these files may not be suitable for diaPASEF fine-tuning |
| `read_all_spectra()` returns spectra without DIA window info | Use `SpectrumReader.new_with_span_step(d_folder, 1, 1)` to get spectra with `precursor_mz` and `isolation_width` fields. ~7.7x more spectra (174k vs 21k) due to per-window expansion |
| Training data matches 0 PSMs to spectra | DIA window matching must use actual `isolation_width` from span-step reader, not a fixed tolerance. Also check RT units (minutes vs seconds) |
| IM-enhanced model incompatible with old checkpoint | Zero-initialize the 4th embedding weight column after loading. Ensures exact match with 3-channel model when IM=0 |
| ASF file has wrong number of columns | IM-enhanced ASF uses 5 columns (mz, intensity, charge_flag, peak_type, im_value). Original uses 4. Check `n_features` parameter in model config |

## Validation

### mzML vs Native .d Comparison

Both paths process the same `.d` file (`03232025__100SPD_DIA-LV100_S3-D4_1_20588.d`, 12 GB bovine liver):

| Metric | mzML Path | Native .d Path (v3) |
|--------|-----------|----------------------|
| Augmented spectra | 90,264 | 90,264 |
| Inference batches | 252 (batch=64) | ~1,411 (batch=64) |
| GPU time | ~4.5 min | ~22 min |
| Preprocessing time | 45-90 min (msconvert) + ~10 min (augment) | ~2 min (timsrust) |
| Total time | 60-105 min | ~24 min |
| Peptides found (threshold 0.5) | 44 | 738 |

**v3 results**: The native `.d` path now produces **738 peptides** (16.8x more than the 44 from mzML). The dramatic improvement is due to: (1) correct MS1 frame reading with Sage tof-to-m/z formula, (2) proper DIA window matching using `SpectrumReader.new_with_span_step()` which expands MS2 spectra with DIA window metadata (174k spectra from 21k raw spectra, ~7.7x expansion), (3) correct b/y ion filtering. The mzML path likely underperformed because `msconvert` loses some timsTOF-specific metadata during conversion.

### Training Data Sources

#### DIA Training Data (self-supervised from DIA-NN results)

| Dataset | Species | Files | Report |
|---------|---------|-------|--------|
| Berger DIA | Porcine | 87 .d (older DIA, not diaPASEF) | `/service/on_campus/Berger/dia/out/report.parquet` |
| Kim | Human | 65 .d | `/service/on_campus/Kim-Jinhwan/out/report.parquet` |
| Zhao | California mouse | 43 .d (validation set) | `/service/on_campus/Zhao_Gemma/.../report.parquet` |
| Liver (bovine) | Bovine | 60 .d | `.../Vinning-paul/Liver/Bovine_Liver__20260330_1149/report.parquet` |
| Muscle (bovine) | Bovine | 60 .d | `.../Vinning-paul/Muscle/Liver_Muscle_research_20260330_1735/report.parquet` |

#### ddaPASEF Pre-Training Data (public repositories)

| Dataset | PXD ID | Species | Description | Status |
|---------|--------|---------|-------------|--------|
| Meier 2018 | PXD010012 | Human (HeLa) | Original ddaPASEF paper, 200ng HeLa on timsTOF Pro | Downloading |
| Prianichnikov 2020 | PXD014777 | Human (HeLa) | MaxQuant.Live ddaPASEF, 200ng HeLa, various gradient lengths | Downloading |

ddaPASEF data provides clean, isolated precursor spectra with unambiguous peptide assignments from database search, making it ideal for pre-training or fine-tuning the Cascadia model specifically for timsTOF fragmentation patterns. Unlike DIA data where multiple precursors co-fragment in wide isolation windows, ddaPASEF isolates individual precursors.

### SpectrumReader.new_with_span_step() — DIA Window Expansion

A critical discovery: `timsrust.read_all_spectra()` returns raw MS2 spectra without DIA window metadata (precursor m/z and isolation width). The `SpectrumReader.new_with_span_step()` constructor provides this by reading the DIA method definition and expanding each raw spectrum into multiple entries, one per DIA window that covers it.

```python
from timsrust_pyo3 import SpectrumReader

# Standard reader: ~21k spectra (raw MS2 frames)
reader_basic = SpectrumReader(d_folder)

# Span-step reader: ~174k spectra (expanded with DIA window info)
reader_dia = SpectrumReader.new_with_span_step(d_folder, 1, 1)
```

The 7.7x expansion (21k to 174k) occurs because each MS2 frame covers multiple DIA windows, and the span-step reader creates one spectrum per window. Each expanded spectrum has correct `precursor_mz` and `isolation_width` fields, which are essential for matching DIA-NN PSMs to raw spectra during training data preparation.

### Training Pipeline Details

#### DIA Window Matching Fix
The original `prepare_training_data.py` matched PSMs to spectra using a fixed 25 Da tolerance centered on the precursor m/z. This failed because DIA windows are defined by their center m/z and width (e.g., center=500, width=25 means 487.5-512.5). The fix uses the actual `isolation_width` from `new_with_span_step()` spectra:

```python
window_half = spectrum.isolation_width / 2
if abs(spectrum.precursor_mz - psm_mz) <= window_half:
    # This spectrum's DIA window covers this precursor
```

#### b/y Ion Quality Filter
Training data quality is critical for de novo model performance. A 10% relative intensity threshold filters out noise peaks that would otherwise confuse the model:

```python
# Keep only b/y ions above 10% of the base peak
max_intensity = max(spectrum.intensities)
threshold = 0.10 * max_intensity
filtered = [(mz, i) for mz, i in zip(mzs, intensities) if i >= threshold]
```

The 10% threshold is appropriate for timsTOF data which has higher baseline noise than Orbitrap. Lower thresholds (1-5%) retain too many noise peaks; higher thresholds (20%+) lose real fragment ions.

## Ion Mobility Enhancement

### Motivation
timsTOF instruments measure ion mobility (1/K0) for every ion, providing an additional dimension of separation. The base Cascadia model ignores this information because it was trained on Orbitrap data. Adding ion mobility as a 4th input channel should improve peptide identification on timsTOF data by leveraging the additional selectivity.

### Model Architecture Changes

The Cascadia `AugmentedSpec2Pep` model uses a multi-channel embedding approach. Each spectrum position has 3 input channels: m/z, intensity, and a channel flag (0=MS1, 1=MS2). The IM enhancement adds a 4th channel for ion mobility.

**Key files modified**:
- `cascadia/models.py` — `AugmentedSpec2Pep.__init__()` and `forward()`
- `cascadia/augment.py` — ASF writing (5-column format)
- `bruker_patch.py` — IM extraction from timsrust spectra

**Architecture detail**: The embedding layer changes from `nn.Linear(3, d_model)` to `nn.Linear(4, d_model)`. The 4th dimension receives the mean 1/K0 value for each spectrum peak (0.0 for Orbitrap/mzML data where IM is unavailable).

### Zero-Initialization Trick
To maintain backward compatibility with the pre-trained checkpoint (trained without IM), the 4th column of the embedding weight matrix is initialized to zero:

```python
# After loading checkpoint:
with torch.no_grad():
    model.embedding.weight[:, 3] = 0.0  # or model.embedding.bias if present
```

This ensures that when IM=0.0 (Orbitrap data), the model produces **exactly** the same output as the original 3-channel model. The IM pathway only contributes when IM > 0 (timsTOF data). Verified by test: predictions match to floating-point precision with zero-init.

### 5-Column ASF Format
The augmented spectrum file format extends from 4 columns to 5:

```
# Original 4-column: mz, intensity, charge_flag, peak_type
# IM-enhanced 5-column: mz, intensity, charge_flag, peak_type, im_value

# Example lines:
500.2847  0.8432  2  1  0.7821    # MS2 peak with IM
500.2847  0.8432  2  0  0.0000    # MS1 peak (no IM available)
```

The `im_value` is the mean 1/K0 across all scans contributing to that spectrum peak. For MS1 frames, IM is set to 0.0 (MS1 survey scans aggregate across all mobility).

### IM Test Results
All 5 unit tests pass:
1. **4-channel embedding creation** -- model initializes with `nn.Linear(4, d_model)`
2. **Forward pass shape** -- output tensor has correct sequence length and vocabulary size
3. **Zero-init backward compatibility** -- with IM=0 input, output matches 3-channel model exactly
4. **IM sensitivity** -- non-zero IM values produce different predictions than IM=0
5. **Mixed batch handling** -- batch with some IM=0 and some IM>0 rows processes correctly

## DE-LIMP R Integration

The R-side integration (feature/cascadia-denovo branch) includes:

- `R/helpers_denovo.R` — `parse_cascadia_ssl()`, `classify_denovo_peptides()`, `run_diamond_blast()`, `generate_cascadia_sbatch()`
- `R/server_denovo.R` — Full server module (SSL upload, confirmed/novel peptide tables, DIAMOND BLAST, job submission, AI context)
- `R/ui.R` — "De Novo" dropdown in navbar with SSL upload in sidebar

### SSL File Format
```
file    scan    charge    sequence    score-type    score    retention-time    start-time    end-time
sample.d    1234    2    PEPTIDER    UNKNOWN    0.85    45.2    43.1    47.3
```

### Peptide Classification
- **Confirmed**: Peptide sequence matches a DIA-NN identification (database search agreement)
- **Novel**: Peptide not found in DIA-NN results (potential novel sequence)
- Classification uses `Stripped.Sequence` matching between SSL and DIA-NN report

## References

- Cascadia: Noble Lab, University of Washington — https://github.com/Noble-Lab/cascadia
- timsrust: Lazear — https://github.com/lazear/timsrust
- Sage: Lazear — https://github.com/lazear/sage
- DIA-NN: Demichev et al. (2020) Nature Methods 17:41-44
- DIAMOND: Buchfink et al. (2015) Nature Methods 12:59-60
