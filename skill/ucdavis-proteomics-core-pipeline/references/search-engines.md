# Search engines reference

Acquisition detection, engine routing, and normalizing each engine's output to the
DE-input contract.

## Acquisition detection (`detect_acquisition.py`)

| format | how DIA/DDA is decided | instrument |
|---|---|---|
| Bruker `.d` | `analysis.tdf` SQLite: `DiaFrameMsMsInfo`/`DiaFrameMsMsWindowGroups` or `Frames.MsMsType==9` → DIA; `PasefFrameMsMsInfo`/`MsMsType==8` → DDA | `GlobalMetadata.InstrumentName` (e.g. "timsTOF Pro") |
| `.mzML[.gz]` | stream MS2 isolation windows: median width ≥3 Da over few centers → DIA; ≤2 Da, many centers → DDA | none (mzML rarely carries model reliably) |
| Thermo `.raw` | needs ThermoRawFileParser → mzML, else `unknown` | ThermoRawFileParser `metadata` model string, if present |
| `.wiff` | convert to mzML first | none |

Confidence is `high`/`medium`/`low`. **Anything not `high`, mixed, or with
disagreeing instruments sets `needs_confirmation` — ask the user.** Instrument may
be null; that's fine (matcher falls back to score 0 and the user confirms).

## Parameters

Search parameters are derived from the data type (instrument + acquisition) by
`estimate_params.py`, not hand-maintained per workflow. → `references/parameters.md`.

## Routing

Default by acquisition: **DIA → DIA-NN, DDA → Sage.** `--engine` overrides. The
bundle's `engine.name` is authoritative when present. FragPipe is opt-in only (the
bundle names it or the user asks) because MSFragger/IonQuant are license-gated.

## Per-engine invocation & output adapter (`run_search.py`)

### DIA-NN (native contract — no adapter)
```
<cmd> --cfg <bundle diann.cfg> --f <file> [--f <file> ...] \
      --fasta <fasta> --out report.parquet --threads N
```
`report.parquet` is already the DE contract — `run_de.R` reads it directly.

### Sage (mzML-first; adapter required)
1. Convert `.d`/`.raw` → mzML with `msconvert` if needed (fails loudly if msconvert
   is absent and inputs aren't mzML).
2. `<cmd> <bundle sage_config.json> -f <fasta> -o <out> --parquet
   --disable-telemetry-i-dont-want-to-improve-sage <mzml...>`
3. **Adapter:** map `lfq.parquet` (protein, filename, intensity) →
   DIA-NN-shaped `report.parquet` with `Run, Protein.Group, PG.MaxLFQ` and zeroed
   Q-value columns (Sage already FDR-filtered). Q-values default to 0 so the
   MaxLFQ DE path keeps every row.
   - v0.14.x has no native protein grouping → may need `sage_protein_groups.py`
     post-hoc (DE-LIMP keeps it on HIVE). v0.15+ has IDPicker grouping.

### FragPipe (opt-in; adapter required)
- Build an `.fp-manifest` from the files + data type (DIA/DDA).
- `<cmd> --headless --workflow <.workflow> --manifest <m> --workdir <out>
  [--config-tools-folder $FRAGPIPE_TOOLS_FOLDER]`
- **Adapter:** `combined_protein.tsv` per-sample `MaxLFQ Intensity` columns →
  DIA-NN-shaped `report.parquet`.
- Needs Java 9+; MSFragger/IonQuant must already be licensed/present.

## The adapters are the test surface

DIA-NN's path is fully wired and native. The **Sage and FragPipe → DE-contract
adapters are the part that needs real-data validation** (flagged in each
workflow's `VALIDATION.md`). Verify the protein×run matrix shape and that
intensities are in the expected (linear, pre-log) scale before trusting DE output.
