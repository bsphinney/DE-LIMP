# The MaxLFQ + QuantUMS path does not implement Moschem et al.

**Status:** partially fixed. The FDR gap, the normalisation crash and the missing
broadcast assertion are fixed in this PR. **Re-quantification is NOT implemented** —
it needs a maintainer decision (new dependency), so this document records the defect,
the evidence, and a working reference implementation.

Found 2026-08-04 while running a 276-sample three-way DE comparison on a mouse DIA
dataset (short-course sample-prep study, HUPO 2026 poster).

---

## The defect

`build_maxlfq_pipeline()` in `R/helpers.R` (and its port,
`skill/ucdavis-proteomics-core-pipeline/scripts/build_maxlfq.R`):

1. filters **precursor rows** on `Empirical.Quality >= eq_cutoff`, then
2. rolls up with `summarise(PG.MaxLFQ = max(PG.MaxLFQ, na.rm = TRUE))`, reusing
   **DIA-NN's precomputed protein-group quantity**.

The code's own comment already states the decisive fact — *"DIA-NN broadcasts it across
rows of a PG within a run"*. Because the value is broadcast, **filtering precursor rows
cannot change it.** DIA-NN computed `PG.MaxLFQ` from *all* precursors, including every
one the filter just discarded, so the retained number still embodies the rejected
precursors.

The eQ filter therefore never improves a retained quantity. It only decides whether a
protein × run **cell exists at all**, under an "any precursor passed" rule that bears no
relationship to the quality of the quantity kept.

This matters because `maxlfq_pipeline_descriptor()` labelled the pipeline
*"MaxLFQ + limma (Moschem 2025)"* and cited
[Moschem et al., *J. Proteome Res.* 2025, 24, 3860](https://doi.org/10.1021/acs.jproteome.5c00009).
As implemented it is not that method. This PR corrects the label and citation.

## Evidence

Measured on a DIA-NN 2.x report (12,931,480 precursor rows, 71 columns, 276 runs),
sampling one row group — 2,513,435 rows covering 772,285 protein × run cells:

| Check | Result |
|---|---|
| Distinct `PG.MaxLFQ` values per (Protein.Group, Run) | **1 in all 772,285 cells** — broadcast confirmed |
| Cells surviving `Empirical.Quality >= 0.75` | 707,369 (91.6%) |
| Of those survivors, fraction whose `PG.MaxLFQ` is **unchanged** | **100.00%** |

## Root cause: the two QuantUMS scores are at different levels

| Score | Varies within a (Protein.Group, Run) cell? | Level |
|---|---|---|
| `Empirical.Quality` | **Yes** — 470,937 cells had >1 distinct value | precursor |
| `PG.MaxLFQ.Quality` | **No** — 0 of 772,285 | protein × run |

So `pgQ` is a coherent **cell** filter and works correctly as written. `eQ` is a
**precursor** filter, and filtering precursors is only meaningful if you then
**re-quantify from the survivors**.

## Scope: the limpa / DPC path is unaffected

`filter_quantums_parquet()` hands the filtered parquet to `limpa::readDIANN()`, and
DPC-Quant quantifies proteins **from precursors**. There an eQ filter genuinely changes
the result. The defect is confined to the MaxLFQ path.

## Reference implementation (not yet wired in)

`iq::fast_MaxLFQ` is the reference implementation of the MaxLFQ algorithm and installs
cleanly alongside the existing container R:

```r
# after the eQ / pgQ filters, instead of reusing PG.MaxLFQ:
r  <- rows[is.finite(rows$Precursor.Normalised) & rows$Precursor.Normalised > 0, ]
nd <- list(protein_list = as.character(r$Protein.Group),
           sample_list  = as.character(r$Run),
           id           = as.character(r$Precursor.Id),
           quant        = log2(r$Precursor.Normalised))
E_pre <- iq::fast_MaxLFQ(nd)$estimate      # log2, proteins x runs
```

Requires `Precursor.Id` and `Precursor.Normalised` in the selected columns.

**Suggested validation when wiring this in:** run it with the filter **off** and compare
against DIA-NN's own `PG.MaxLFQ`. They should track closely; if they do not, the
re-quantification is wrong.

**Decision needed from the maintainer:** adding `iq` is a new dependency. The alternative
is to drop `eq_cutoff` from the MaxLFQ path entirely and expose only `pgq_cutoff`, which
is coherent as written — at the cost of not offering the paper's method at all.

## Two further defects in the same function (fixed in this PR)

- **No experiment-wide protein FDR.** Only `Q.Value`, `Lib.Q.Value` and `Lib.PG.Q.Value`
  were filtered. Across 276 runs, a union of run-level-passing IDs sits well above the
  nominal 1% experiment-wide, and the inflated protein list also inflates the family size
  `m` that BH later corrects over. `PG.Q.Value`, `Global.Q.Value` and `Global.PG.Q.Value`
  are now filtered when present.
- **Quantile normalisation crashes on a near-empty run.** `limma::normalizeQuantiles()`
  calls `approx()` per column and fails with *"need at least two non-NA values to
  interpolate"* if a run retains <2 quantified proteins. Observed: at eq/pgq ≥ 0.75 one
  run lost every row (so it never became a column and the job survived at 275 runs), but
  at ≥ 0.50 the same run survived with a single value and killed the analysis. Degenerate
  runs are now dropped explicitly and named in a warning.

## Not fixed here, but worth a look

Contaminants are never removed. In our run 93 `Cont_` protein groups reached the DE
tables and 47–53 were called significant, including pig trypsin `Cont_P00761`
(log2FC −7.2) and bovine `Cont_P00760` (+3.0). Trypsin autolysis genuinely differs
between digestion protocols, so these are real signal — but they are being reported as
biological prep effects and are inflating the BH denominator. In that dataset every
non-mouse entry sat inside a `Cont_` group, so a single filter would address both.
