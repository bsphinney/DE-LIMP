# Validation — diann_astral_dia_human

> EXAMPLE provenance. Replace with your real validation record before this
> workflow is trusted in production.

## What this workflow is
Library-free DIA-NN search of human narrow-window DIA from the Orbitrap Astral,
QuantUMS quantification, then limpa DPC-Quant + limma DE.

## Pinned software
- DIA-NN 2.6.0 (Academia), reads Thermo `.raw` natively (2.1+).
- limpa (Bioconductor, R 4.5+) + limma for DE.

## Search parameters
See `diann.cfg`. Runtime adds `--f <files> --fasta <fasta> --out report.parquet
--threads <N>`. Mass accuracy left on auto (`--mass-acc 0 --mass-acc-ms1 0`) so
DIA-NN calibrates per run — appropriate for Astral narrow-window data.

## Validation dataset (CONFIRM)
- Sample: HeLa QC digest, 3 technical replicates.
- Expected: ~9,000–10,000 protein groups at 1% FDR (fill in your observed range).
- CV / ID depth acceptance criteria: (fill in).

## Reproducibility
DE step recorded with the exact recipe:
```
dat <- limpa::readDIANN("report.parquet", format="parquet", q.cutoffs=0.01)
y   <- limpa::dpcQuant(dat, "Protein.Group", dpc=limpa::dpcCN(dat))
fit <- limpa::dpcDE(y, design, plot=FALSE) |> contrasts.fit(...) |> eBayes()
```

## Change log
- 2026-05-01: initial validation (EXAMPLE).
