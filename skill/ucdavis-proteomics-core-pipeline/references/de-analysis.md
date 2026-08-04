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

**`--logfc` is a reference line, not a filter.** Significance is `adj.P.Val < --adjp`
(BH) and nothing else. `--logfc` is drawn and labelled on the volcano so a reader can see
effect size, and it is reported descriptively ("N of the significant proteins are also
≥2-fold"), but it never removes a protein from the significant set. This keeps the
reported claim identical to the hypothesis `eBayes`/`topTable` actually tested
(H0: log2FC = 0). A post-hoc `|logFC|` cut would assert "changed by at least Nx" while
the FDR only covers "changed at all" — and because the observed fold change is a point
estimate, proteins whose true effect is below the line pass it routinely. If a genuine
fold-change threshold is ever wanted, use `limma::treat()` + `topTreat()`, whose interval
null tests it properly; note that `treat` is far more stringent than the same cut applied
post hoc, so limma recommends a *small* threshold there (fc 1.1–1.5, not 2).

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

## Method choice — limpa/DPC is the default

`--method dpc` (limpa) is the default and should stay that way. It models the
detection-probability curve and quantifies from precursor intensities, so it uses the
whole measurement instead of a pre-collapsed protein number, and it handles missingness
as information rather than as a hole to filter around.

`--method maxlfq` is for two situations only:

1. **The user asked** — usually because they want QuantUMS quality filtering, or to
   match an earlier MaxLFQ analysis.
2. **The file you point at is protein-level.** `readDIANN()` keys on `Precursor.Id`
   and `Precursor.Normalised`. This is a property of the FILE, not of the engine —
   most engines can feed limpa if you use the right output:

   - **DIA-NN** — `search_out/report.parquet`, native.
   - **FragPipe** — `dia-quant-output/report.tsv` with `--format tsv`. Its DIA route
     bundles DIA-NN, so this is a DIA-NN report and carries `Lib.Q.Value` /
     `Lib.PG.Q.Value` (columns 39-40). Verified end-to-end on the 9-file class data:
     22,425 precursors -> 12,485 proteins, both contrasts written.
   - **Radiant** — `radiant_to_delimp.py` output. Verified twice: 3-run HeLa
     (25,722 precursors x 3) and 18-run Poplar (21,370 precursors x 18).

   The ADAPTED `report.parquet` from `adapt_*` is the exception: it collapses to one
   row per protein x run on purpose, to feed the maxlfq path. `run_de.R` checks the
   columns before limpa is called and names the precursor-level file to use instead.

**q-value columns.** `readDIANN()` prints `Q-value columns <x> not found.` and then
continues with that filter simply not applied — a message, not an error, easy to miss
in a long log. `run_de.R` resolves `q.columns` against the real header, states which
filters it applied, and stops if none are usable, so an unfiltered result can never
look filtered.

Do not read a bundle's `de.method: maxlfq` as a recommendation — it records what that
engine's adapted output can support, not which method is better.

## Coverage filter (maxlfq path)

`--coverage-min` (default 0.5) drops proteins quantified in fewer than that fraction
of samples **before** limma. The MaxLFQ matrix keeps every protein seen in *any* run,
so rows with 1-2 finite values otherwise reach `eBayes`, which then moderates variance
against rows whose variance is barely estimable — destabilising the whole fit, not just
those rows. Dropped proteins are an on/off observation, not a differential-abundance
result, and the count lands in `filters_applied` so the Methods text says so.

Ported from DE-LIMP (`server_data.R`), whose default and rationale come from the
UC Davis Bioinformatics Core. If it leaves <10 testable proteins the run stops rather
than fitting noise — loosen `--coverage-min`, or the QuantUMS cutoffs if those are
what emptied the matrix.

## Provenance (self-describing — DE-LIMP architectural rule #1)
Each method path returns a `descriptor` (pipeline_id, display_label, rollup_method,
de_engine, missing_policy, citation). `methods.txt` is built from it — **never
hardcode a description of what ran**, and hand `methods.txt` to the user verbatim.

## Citations (verified June 2026)
- **limpa / DPC:** Li M, Cobbold SA, Smyth GK (2025) bioRxiv 10.1101/2025.04.28.651125;
  Li M, Smyth GK (2023) Bioinformatics 39(5):btad200. (DE-LIMP's
  `dpc_pipeline_descriptor()` mis-cites this as "Law CW, Smyth GK" — fix upstream.)
- **MaxLFQ path:** DIA-NN MaxLFQ (Demichev et al. 2020, Nat Methods 17:41) +
  limma (Ritchie et al. 2015, NAR 43:e47).
- **QuantUMS quality filtering:** da Cruz Moschem J, Silva Campitelli de Barros BC,
  de Toledo Serrano SM, Chaves AFA (2025) *Decoding the Impact of Isolation Window
  Selection and QuantUMS Filtering in DIA-NN for DIA Quantification of Peptides and
  Proteins.* J Proteome Res 24:3860-3873. doi:10.1021/acs.jproteome.5c00009.
  **VERIFIED 2026-08-04** — an earlier note here said this "could not be verified";
  that was wrong. DE-LIMP cites it as "Moschem et al." (the first author's surname is
  da Cruz Moschem). QuantUMS computes three scores: protein-group MaxLFQ quality,
  empirical quality, and quantity quality, all measuring MS1/MS2 feature agreement.
  The skill filters on the first two (`--pgq-cutoff`, `--eq-cutoff`).
