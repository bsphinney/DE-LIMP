# MaxLFQ + QuantUMS path: what it does, and where it departs from Moschem et al.

**RETRACTION NOTICE.** The first version of this document (and the first two commits on
this branch) claimed the pipeline "does not implement Moschem et al." because the QuantUMS
filter censors cells instead of re-quantifying. **That claim was wrong**, and it is
retracted here. Reading the paper and the authors' published code shows they do the same
thing. The measurements below are unchanged and still correct; only the interpretation was
wrong. What survives are three genuine, smaller departures from the paper, plus two
independent bugs.

---

## What was measured (unchanged, and still true)

On a DIA-NN 2.x report (12,931,480 precursor rows, 276 runs), one row group sampled =
772,285 protein × run cells:

| Check | Result |
|---|---|
| Distinct `PG.MaxLFQ` per (Protein.Group, Run) | **1 in all 772,285 cells** — broadcast confirmed |
| Cells surviving `Empirical.Quality >= 0.75` | 707,369 (91.6%) |
| Of those survivors, fraction whose `PG.MaxLFQ` is unchanged | **100.00%** |

And the two QuantUMS scores sit at different levels:

| Score | Varies within a cell? | Level |
|---|---|---|
| `Empirical.Quality` | yes (470,937 cells) | precursor |
| `PG.MaxLFQ.Quality` | no (0 of 772,285) | protein × run |

So the filter decides **which cells survive**, never their values.

## Why that is NOT a deviation from the paper

Moschem et al., *J. Proteome Res.* 2025, 24, 3860 is about **applying cutoffs to the
QuantUMS scores** — the abstract says filtering "removes proteins with low abundances and
high coefficients of variation". It is a filtering study, not a re-quantification study.

The decisive evidence is the corresponding author's published workflow
(<https://github.com/41ison/limma-for-proteomics>, Alison Felipe Alencar Chaves):

```r
diann_report <- arrow::read_parquet("report.parquet") %>%
    dplyr::filter(Lib.PG.Q.Value <= 0.01 & Lib.Q.Value <= 0.01 & PG.Q.Value <= 0.01) %>%
    dplyr::filter(str_detect(Protein.Ids, "cRAP", negate = TRUE))

prot_mtx <- diann::diann_matrix(diann_report,
    id.header = "Protein.Ids",
    quantity.header = "Genes.MaxLFQ.Unique",
    proteotypic.only = T,
    pg.q = .01)
```

`diann::diann_matrix()` **pivots a precomputed column**. `diann::diann_maxlfq()` is the
function that recomputes MaxLFQ from precursor quantities — the package README documents
both, describing `diann_maxlfq` as being for "MaxLFQ-based protein quantification after
manual precursor-level filtering". **The authors chose the pivot.** So filter-then-pivot is
the reference behaviour, and this pipeline's approach is correct in kind.

### Corollary: post-hoc re-quantification is also the wrong fix

We implemented `iq::fast_MaxLFQ` re-quantification and tested it with the filter **off**,
where it should reproduce DIA-NN's `PG.MaxLFQ`. It does not:

| `iq::fast_MaxLFQ` input | median per-protein r vs `PG.MaxLFQ` | frac r > 0.9 |
|---|---|---|
| `Precursor.Normalised` | 0.598 | 5.5% |
| `Precursor.Quantity` | 0.327 | 0.0% |
| `Precursor.Normalised`, proteotypic | 0.581 | 5.0% |
| `Precursor.Quantity`, proteotypic | 0.325 | 0.0% |

(400 protein groups present in ≥200 of 276 runs; correlation per protein across runs, so a
constant per-protein offset cannot explain the gap.) DIA-NN 2.x computes protein quantities
with **QuantUMS**, not the classic MaxLFQ algorithm — the column keeps the legacy name — so
`iq` is a *different estimator*, not a recomputation. Useful by-product: the ~2× advantage
of `Precursor.Normalised` confirms DIA-NN's protein quantity is built from normalised
precursor quantities.

---

## Where this pipeline genuinely DOES depart from the paper

1. **All peptides vs proteotypic only.** The paper uses `Genes.MaxLFQ.Unique` with
   `proteotypic.only = TRUE`. This pipeline uses `PG.MaxLFQ`, which includes shared
   peptides. Different entity, different quantity.
2. **`PG.Q.Value` was never filtered.** The paper filters `Lib.PG.Q.Value`, `Lib.Q.Value`
   **and `PG.Q.Value`** at 0.01. This pipeline filtered only the first two. *Fixed in this
   PR*, which also adds the experiment-wide `Global.PG.Q.Value` / `Global.Q.Value`.
3. **Contaminants were never removed.** The paper explicitly strips cRAP. This pipeline
   does not. In our 276-run dataset 93 `Cont_` groups reached the DE tables and 47–53 were
   called significant, including pig trypsin `Cont_P00761` (log2FC −7.2) and bovine
   `Cont_P00760` (+3.0) — real signal, but reported as biological prep effects and
   inflating the BH denominator. **Not fixed here**; it needs a decision about where the
   contaminant list lives.

## Two independent bugs (fixed in this PR)

- **Quantile normalisation crashes on a near-empty run.** `limma::normalizeQuantiles()`
  calls `approx()` per column and fails with *"need at least two non-NA values to
  interpolate"* if a run retains <2 quantified proteins. Observed: at eq/pgq ≥ 0.75 one run
  lost every row (so it never became a column and the job survived at 275 runs), but at
  ≥ 0.50 the same run survived with one value and killed the analysis. Degenerate runs are
  now dropped explicitly and named.
- **The `PG.MaxLFQ` broadcast was assumed, not asserted.** A silent failure would bias the
  affected cells upward with nothing to detect it. Now warns.

## Separate issue worth a maintainer's attention: quantile normalisation and run depth

Not a Moschem question, but it bites the same function. `normalizeBetweenArrays(method =
"quantile")` maps every run onto the full common quantile grid. When a run is shallow its
survivors are its most abundant proteins, so they get spread across the whole dynamic
range. In our dataset 16 runs had <2000 quantified proteins (concentrated in one prep
group, so the distortion is group-dependent and enters the contrast). Measured: excluding
the 15 shallow runs that remained after dropping a failed injection moved one contrast's
log2FC by a median of **−0.403 log2** across the 500 most abundant proteins — larger than
that contrast's median |log2FC| of 0.327. Median-centring on proteins present in ≥70% of
runs avoids this without discarding the shallow runs. Worth exposing as an option.
