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

## IMPORTANT: post-hoc re-quantification is NOT equivalent — measured

The obvious fix is "recompute MaxLFQ from the surviving precursors." We implemented that
with `iq::fast_MaxLFQ` and **measured how well it reproduces DIA-NN's `PG.MaxLFQ` with the
filter switched off**, which is the control that has to pass before the filtered version
means anything. It does not pass:

| `iq::fast_MaxLFQ` input | median per-protein r vs DIA-NN `PG.MaxLFQ` | fraction r > 0.9 |
|---|---|---|
| `Precursor.Normalised` | **0.598** | 5.5% |
| `Precursor.Quantity` | 0.327 | 0.0% |
| `Precursor.Normalised`, proteotypic only | 0.581 | 5.0% |
| `Precursor.Quantity`, proteotypic only | 0.325 | 0.0% |

(400 protein groups present in ≥200 of 276 runs; correlation computed per protein across
runs, so a constant per-protein offset cannot explain the gap.)

**Why:** DIA-NN 2.x does not compute `PG.MaxLFQ` with the classic MaxLFQ algorithm. It
computes protein quantities with **QuantUMS**, its uncertainty-minimising estimator; the
column retains the legacy name. So substituting `iq`'s classic MaxLFQ does not "recompute
the same number from fewer precursors" — it swaps in a **different estimator**.

Consequence: a post-hoc `iq` re-quantification arm conflates two changes (the filter *and*
a different quantification algorithm) and cannot be presented as "the Moschem workflow
done correctly". The higher `Precursor.Normalised` score does at least confirm DIA-NN's
protein quantity is built from normalised precursor quantities, not raw ones.

**The only faithful route is to apply the quality filtering inside DIA-NN** so that DIA-NN
recomputes its own QuantUMS protein quantities from the surviving precursors. Any
post-hoc re-quantification in R should be labelled as a different estimator, not as the
published method.

## Reference implementation (for the different-estimator route only)

If a post-hoc route is wanted anyway — clearly labelled as classic MaxLFQ, not QuantUMS —
`iq::fast_MaxLFQ` installs cleanly alongside the existing container R:

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

**Decision needed from the maintainer.** Three options, in order of faithfulness:

1. **Apply the QuantUMS filtering inside DIA-NN** and let it recompute its own protein
   quantities. This is the only route that actually reproduces the published method. It
   means the filter becomes a search-time setting rather than a post-processing toggle.
2. **Drop `eq_cutoff` from the MaxLFQ path** and expose only `pgq_cutoff`, which is
   coherent as written because pgQ is already a per-cell score. Honest, cheap, and loses
   nothing that currently works — but does not offer the paper's method.
3. **Post-hoc `iq` re-quantification**, clearly labelled as *classic MaxLFQ*, accepting
   that it is a different estimator from DIA-NN's QuantUMS (median per-protein r ≈ 0.6
   against `PG.MaxLFQ` even with no filtering) and adds a dependency.

Option 2 is the smallest correct change. Option 1 is the right long-term answer.

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
