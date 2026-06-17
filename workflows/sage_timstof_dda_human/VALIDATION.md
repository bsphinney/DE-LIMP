# Validation — sage_timstof_dda_human

> EXAMPLE provenance. Replace with your real record.

## What this workflow is
Sage database search of human ddaPASEF (timsTOF), label-free quant, then
MaxLFQ + limma DE.

## Pinned software
- Sage 0.15.0 (CONFIRM the version you validated). Telemetry disabled at runtime
  with `--disable-telemetry-i-dont-want-to-improve-sage`.

## Two adapter steps that need real-data testing
This is the honest gap, flagged so it gets tested rather than assumed:

1. **`.d` → mzML.** Sage reads mzML; ddaPASEF `.d` must be converted first
   (msconvert, or a Sage build with the Bruker reader). The skill's search
   router handles conversion when needed, but the conversion settings should be
   validated on your data.
2. **Sage → DE input.** Sage emits `results.sage.parquet` + `lfq.parquet`, not a
   DIA-NN report. The MaxLFQ DE path expects a protein × run matrix, so the
   router maps `lfq.parquet` (protein LFQ per run) into that matrix before
   `run_de.R --method maxlfq`. Confirm column mapping on a known dataset.

## Search parameters
See `sage_config.json`. The FASTA is injected at runtime (`-f`), so the config's
`fasta` field is a placeholder.

## Validation dataset (CONFIRM)
- Sample, expected PSMs / protein groups, acceptance criteria: (fill in).

## Change log
- 2026-05-15: initial validation (EXAMPLE).
