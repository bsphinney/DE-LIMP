# Cascadia Native Bruker .d Support — Proposed Contribution

## Summary

Native Bruker timsTOF `.d` file support for Cascadia, eliminating the need for mzML conversion. Uses `timsrust_pyo3` (open-source Rust-based TDF reader) as a drop-in replacement for the mzML parser front-end.

**Impact**: 10-20x faster data loading (4 min vs 45-90 min per file), zero disk overhead from mzML intermediates, no ProteoWizard/Bruker SDK dependency.

## Motivation

Cascadia currently requires mzML input, which means Bruker timsTOF users must convert `.d` files before running de novo sequencing. This conversion step:
- Takes 45-90 minutes per file via msconvert/TIMSCONVERT
- Creates 3-8x disk overhead (temporary mzML files)
- Requires ProteoWizard or Bruker SDK installation
- Is the primary barrier to adopting Cascadia in timsTOF workflows

With timsTOF being one of the dominant DIA platforms, native support would significantly expand Cascadia's user base.

## Implementation

### Approach

Replace only the mzML parsing front-end (`get_centers()` and `extract_spectra()` in `augment.py`) with timsrust equivalents. The transformer model, ASF format, and all downstream code are unchanged.

### New Files

**`cascadia/bruker_patch.py`** (~120 lines)
- `get_centers_bruker(d_folder)` — Reads spectra via `timsrust_pyo3.read_all_spectra()`, builds the same `f_to_mzrt_to_pep` dict that `get_centers()` returns from mzML
- `extract_spectra_bruker(d_folder, ...)` — Extracts MS1/MS2 pairs into the same `prec_to_spec` format for `write_asf()`

### Routing Change

6-line modification to `cascadia/cascadia.py`:

```python
# In sequence() function, before the augment call:
from pathlib import Path
input_path = Path(spectrum_file)

if input_path.suffix == ".d" or (input_path.is_dir() and (input_path / "analysis.tdf").exists()):
    # Bruker timsTOF — use native loader
    from cascadia.bruker_patch import get_centers_bruker, extract_spectra_bruker
    f_to_mzrt_to_pep, max_mz, window_size, cycle_time = get_centers_bruker(spectrum_file)
    # ... rest uses existing write_asf + inference pipeline
else:
    # Existing mzML path
    asf_file, isolation_window_size, cycle_time = augment_spectra(spectrum_file, temp_path, ...)
```

### Optional Dependency

```toml
# pyproject.toml
[project.optional-dependencies]
bruker = ["timsrust_pyo3>=0.4.0"]
```

Install with: `pip install cascadia[bruker]`

## Benchmarks

Tested on UC Davis HIVE HPC (A100 80GB GPU node):

| Step | mzML Path | Native .d Path |
|------|-----------|----------------|
| **File reading** | N/A (pre-converted) | 51 seconds |
| **Augmentation** | ~45-90 min (conversion) | 4 min 5 sec |
| **Inference** | ~26 min (same) | ~26 min (same) |
| **Total** | ~70-120 min | **~31 min** |
| **Disk overhead** | 3-8x file size | Zero |
| **Memory** | Similar | 34 GB peak |

Test file: timsTOF HT DIA-PASEF, 100 SPD, ~22,500 MS2 spectra, bovine liver tissue.

## timsTOF-Specific Considerations

1. **Ion mobility dimension**: timsrust_pyo3 handles IM by returning spectra already flattened within each isolation window — same as mzML conversion does
2. **diaPASEF window scheme**: Variable windows handled transparently via `precursor.mz` per spectrum
3. **RT units**: timsrust returns minutes; conversion to seconds applied (`rt * 60`)
4. **No proprietary dependencies**: timsrust_pyo3 uses open Rust implementations, not Bruker SDK

## Why timsrust_pyo3

- Open-source (MIT license), no proprietary Bruker SDK needed
- Python bindings via PyO3: `pip install timsrust_pyo3`
- Already used by the Sage search engine for timsTOF support
- Actively maintained by the Lazear lab
- Reads TDF/miniTDF formats natively
- Small accuracy difference vs Bruker SDK (acceptable for de novo sequencing where we're doing broad peptide matching)

## Testing

- [x] timsrust_pyo3 reads real timsTOF .d files (22,566 spectra)
- [x] Native augmentation produces valid ASF format
- [x] Cascadia transformer runs inference on native ASF output
- [x] SSL peptide predictions generated from .d input
- [ ] Peptide overlap comparison: native vs mzML path (>90% expected)
- [ ] Multi-file batch processing
- [ ] Memory optimization for large files (currently 34 GB peak)

## Authors

- Brett Phinney — UC Davis Proteomics Core Facility
- Implementation assisted by Claude (Anthropic)

## Related

- [timsrust_pyo3](https://github.com/lazear/timsrust) — Rust-based TDF reader
- [Cascadia](https://github.com/Noble-Lab/cascadia) — DIA de novo sequencing
- [DE-LIMP](https://github.com/bsphinney/DE-LIMP) — Shiny proteomics pipeline integrating Cascadia
