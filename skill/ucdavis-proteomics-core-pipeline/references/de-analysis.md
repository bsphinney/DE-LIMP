# Differential expression reference

The DE step (`run_de.R` + `build_maxlfq.R`) is a faithful port of what DE-LIMP logs
in `R/server_data.R` / `R/helpers.R`. Two pipelines, picked by the bundle's
`de.method`.

## Invocation
```
Rscript scripts/run_de.R --input report.parquet --metadata conditions.csv \
        --method {dpc|maxlfq} --outdir de_results \
        [--contrasts "B-A,C-A"] [--q-cutoff 0.01] [--logfc 1.0] [--adjp 0.05]
```
`metadata` CSV: `File.Name,Group[,Batch,Covariate1,Covariate2]`. `File.Name` must
match the `Run` / column names in the report. Default contrasts = every group vs
the first factor level.

## `--method dpc` (limpa DPC-Quant + limma) — DE-LIMP default, use with DIA-NN
```r
dat <- limpa::readDIANN("report.parquet", format="parquet", q.cutoffs=0.01)
y   <- limpa::dpcQuant(dat, "Protein.Group", dpc=limpa::dpcCN(dat))
fit <- limpa::dpcDE(y, design, plot=FALSE)        # wraps voomaLmFitWithImputation
fit <- contrasts.fit(fit, makeContrasts(...)) |> eBayes()
topTable(fit, coef=cn, number=Inf, adjust.method="BH")
```
Missing precursors are modelled by the detection-probability curve — **not imputed,
not dropped** — and the imputation uncertainty is propagated into the limma fit.
Needs **R 4.5+ / Bioconductor 3.22+**.

## `--method maxlfq` (MaxLFQ + limma) — use with Sage/FragPipe (or DIA-NN MaxLFQ)
`build_maxlfq.R`: filter `Q/Lib.Q/Lib.PG.Q ≤ q` (+ optional QuantUMS `eQ`/`pgQ`),
pivot `PG.MaxLFQ` to a protein×run matrix, log2, quantile-normalize
(`limma::normalizeBetweenArrays`), then `lmFit → contrasts.fit → eBayes →
topTable(BH)`. NAs are left in place; limma drops them per row.

## DE-input contract (§8.3)
A DIA-NN-shaped report with: `Run, Protein.Group, PG.MaxLFQ, Q.Value, Lib.Q.Value,
Lib.PG.Q.Value` (+ optional `Empirical.Quality, PG.MaxLFQ.Quality, Genes,
Protein.Names`). `run_search.py` produces this for non-DIA-NN engines.

## Design matrix
`~ 0 + groups [+ Batch + Covariate1 + Covariate2]`, colnames = group levels.
**Rank-checked before fitting** (`qr(design)$rank`); fails on confounded covariates
or empty groups. Groups with <2 replicates have no within-group variance — warn the
user at the design step (`collect_conditions.py --validate` flags singletons).

## Provenance (self-describing — DE-LIMP architectural rule #1)
Each method path returns a `descriptor` (pipeline_id, display_label, rollup_method,
de_engine, missing_policy, citation). `methods.txt` is built from it — **never
hardcode a description of what ran**, and hand `methods.txt` to the user verbatim.

## Citations (verified June 2026)
- **limpa / DPC:** Li M, Cobbold SA, Smyth GK (2025) bioRxiv 10.1101/2025.04.28.651125;
  Li M, Smyth GK (2023) Bioinformatics 39(5):btad200. (DE-LIMP's
  `dpc_pipeline_descriptor()` mis-cites this as "Law CW, Smyth GK" — fix upstream.)
- **MaxLFQ path:** DIA-NN MaxLFQ (Demichev et al. 2020, Nat Methods 17:41) +
  limma (Ritchie et al. 2015, NAR 43:e47). (DE-LIMP's "Moschem et al. 2025"
  reference could not be verified — reconcile before publishing methods text.)
