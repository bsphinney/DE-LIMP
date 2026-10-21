# Differential expression reference

The DE step (`run_de.R` + `build_maxlfq.R`) is a faithful port of what DE-LIMP logs
in `R/server_data.R` / `R/helpers.R`. Two pipelines, picked by the bundle's
`de.method`.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 30: Invocation
- 96: `--method dpc` (limpa DPC-Quant + limma) — DE-LIMP default, use with DIA-NN
- 115: `--method maxlfq` (MaxLFQ + limma) — use with Sage/FragPipe (or DIA-NN MaxLFQ)
- 160: DE-input contract (§8.3)
- 165: Design matrix
- 216: Paired / repeated designs — `--block <column>`
- 335: Method choice — limpa/DPC is the default
- 371: Normalised or non-normalised quantities: the experiment type and the data check
- 470: Contaminant filter (both paths, on by default)
- 559: Coverage filter (maxlfq path)
- 573: Provenance (self-describing — DE-LIMP architectural rule #1)
- 586: Top-protein violins (`make_figures.R` → `violin_top_<contrast>.png`)
- 617: Citations (verified June 2026)
- 635: Less-used flags
- 639: Step 8: the normalisation check and the final DE
- 732: Step 8: blocks and numeric covariates
- 765: Step 8: keratin samples and contaminants
- 801: Step 8: the method, what `run_de.R` writes, QuantUMS
- 869: 8b. Figures, and what counts as significant (step 8b)
<!-- /contents -->

## Invocation
```
Rscript scripts/run_de.R --input report.parquet [--format parquet|tsv] \
        --metadata conditions.csv --method {dpc|maxlfq} --outdir de_results \
        [--contrasts "B-A,C-A"] [--q-cutoff 0.01] [--eq-cutoff 0] [--pgq-cutoff 0] \
        [--coverage-min 0.5] [--logfc 1.0] [--adjp 0.05] \
        [--fasta-meta search.fasta.meta.json] [--keep-contaminants] \
        [--numeric-covariates Covariate1[,Covariate2]] \
        [--block Mouse [--block-scope within|all] [--block-effect auto|fixed|random]]
```
- `--q-cutoff` is the q-value filter on the report; `--eq-cutoff` / `--pgq-cutoff` the
  empirical-library and protein-group q-value filters (0 = off; `de_provenance.json` records
  each as applied).
- `--coverage-min` (MaxLFQ only): the fraction of samples a protein must be quantified in to
  be tested.
- `--fasta-meta` (`fetch_fasta.py`'s sidecar) tells the contaminant step whether real
  proteins sat in the database only as `Cont_` entries, which the removal would then take out
  too. `--keep-contaminants` quantifies and tests `Cont_` entries with the sample proteins
  instead of removing them first; `de_provenance.json` records which.
- `--numeric-covariates`: covariates that are continuous measurements, fitted as numbers
  (below, "Continuous covariates").
- `--block` / `--block-scope` / `--block-effect`: paired and repeated designs, below.
`metadata` CSV: `File.Name,Group[,Batch,Covariate1,Covariate2][,<block column>]`. `File.Name` must
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

**A protein measured in neither group of a contrast is not tested in it (2.12.0).** With more
than two groups, a protein can be measured only in groups a contrast does not compare (seen in
C and D, never in A or B). For A vs B it has nothing to test: on dpc both sides are the
detection model's inferred value (logFC ≈ 0, p ≈ 1), on maxlfq both are missing. Before
2.12.0 such a protein was labelled `presence call`, tested, and counted in that contrast's
Benjamini–Hochberg correction — one staff analysis had 186 of its 249 p > 0.975 proteins in
one contrast from this, and a spike at p = 1 in every p-value histogram. Now, per contrast:
- its `Evidence` is `not measured in either group` (`presence call` keeps its meaning:
  measured in one compared group, never in another);
- its `logFC`, `t`, `P.Value`, `adj.P.Val` and `B` are NA — the row stays, so nothing
  disappears — and it is never significant;
- BH runs over the proteins tested in that contrast. The rule depends only on the detection
  counts of the two groups together, not on which group is which, so leaving these out is
  independent filtering (Bourgon et al. 2010) and does not bias the error rate. The model fit
  (and limma's variance prior) is one fit for all contrasts and is not changed;
- `de_provenance.json` records the rule (`not_measured`) and, per DE table, `n_proteins`,
  `n_tested` and `n_not_measured`; `methods.txt`, the Methods paragraph, the analysis page's
  counts, the figure captions and `reproducibility_log.R` all follow it. It applies to every
  `--method` and to blocked fits (technical replicates, `--block`).

The DE-LIMP session (`DE-LIMP_session.rds`) follows the same rule (2.12): DE-LIMP draws its
volcano and tables with `topTable()` from the session's fit, so in that fit an untested protein's
logFC (coefficient), `t`, `p.value` and `lods` are NA for the contrast (its moderated F too, when
any contrast does not test it) and `not_measured` lists them. `topTable()` then matches
`DE_*.csv`: the same NA rows, BH over the tested proteins. DE-LIMP leaves NA values off the
volcano and counts them as not significant; its GSEA, which ranks by logFC, drops them -- a fold
change between two groups that never measured the protein is the model's inference, not data.

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
Needs **R 4.5+ / Bioconductor 3.22+**; setup.sh pins **limpa ≥ 1.4.0**.

**The numbers depend on the limpa version.** On > 2,000 precursors `dpcCN()` fits the curve
on a 2,000-row subset: a seeded random sample in limpa 1.2.x, a systematic one (ordered by
missingness, then mean) in 1.4. On the mouse-brain IP cohort v2 that moved |logFC| by at most 0.064 (median
≤ 0.006) and the significant sets by Jaccard 0.958–1.000. `de_provenance.json` records
`limpa_version`: **compare analyses only within one limpa version** (detail and numbers:
`references/environment.md`, "Results depend on the limpa version").

## `--method maxlfq` (MaxLFQ + limma) — use with Sage/FragPipe (or DIA-NN MaxLFQ)
`build_maxlfq.R`: filter `Q/Lib.Q/Lib.PG.Q ≤ q` (+ optional QuantUMS `eQ`/`pgQ`),
pivot `PG.MaxLFQ` to a protein×run matrix, log2, quantile-normalize
(`limma::normalizeBetweenArrays`), then `lmFit → contrasts.fit → eBayes →
topTable(BH)`. NAs are left in place; limma drops them per row.

**A Sage report is not MaxLFQ, and is not described as MaxLFQ.** `adapt_sage` writes one row
per **peptide** × run, and declares that in the report itself (parquet metadata
`delimp.quantity.*`: level `peptide`, what the value is, Sage's FDR, the citation).
`build_maxlfq.R` takes the highest row per (protein, run), so each protein's value is its
single most intense peptide. `maxlfq_descriptor()` turns the declaration into the one
description every output reads:
- pipeline `peptide_max`, "Highest-peptide intensity + limma";
- a `rollup_method` naming that rollup, and an `identification_fdr` naming Sage's own FDR (the
  report's q-columns are 0.0 placeholders);
- a `caveat` that `audit_results.py` turns into a `quantification` WARN (AUDIT.md, then the
  report's Audit & caveats);
- the AI brief's `plain_language` line.

methods.txt, de_provenance.json, the Methods (`make_methods.py`), the AI brief and
reproducibility_log.R all quote these, so none of them says "DIA-NN" or "MaxLFQ" for a Sage
DE. A real peptide → protein rollup for Sage is a backlog item (`docs/SKILL_OPEN_DEFECTS.md`
in the DE-LIMP repo). DIA-NN's own report declares nothing and is described as DIA-NN
PG.MaxLFQ, as before.

**The other adapters declare too, at level `protein`**:
- `adapt_fragpipe_dda`: FragPipe/IonQuant MaxLFQ from `combined_protein.tsv`, with the
  versions and Philosopher filter read from `fragpipe.workflow`;
- `adapt_alphadia`: AlphaDIA's protein-group LFQ, directLFQ or QuantSelect, with the version and
  FDR read from `frozen_config.yaml`;
- `adapt_radiant`: Radiant's Fulcrum `PG.*` quantity and q-value.

For these, the descriptor's `protein` branch passes the engine's quantity through, labelled
e.g. "FragPipe IonQuant MaxLFQ + limma", with no rollup described, since each is already one
value per protein and run. A Fulcrum output with only a `Precursor.*` quantity is collapsed to
the highest precursor by the adapter. That is the peptide-level rule, so it is declared
`peptide` and described (and caveated) like Sage.

Every declaration also says whether the report's q-value columns are `placeholder` (0.0: Sage,
FragPipe-DDA, AlphaDIA filtered upstream) or `real` (Radiant's Fulcrum q). For `placeholder`, the
methods state the engine's own FDR, and the DE's filter line says it kept every row. For `real`,
the numeric cutoff (e.g. "q <= 0.010", "Q.Value ≤ 0.01") is printed beside the declaration. Citations are
Crossref-verified: MSFragger (Kong 2017), IonQuant (Yu 2021), MaxLFQ (Cox 2014), AlphaDIA
(Wallmann 2025), directLFQ (Ammar 2023), Radiant DIA (Just 2026, bioRxiv).

## DE-input contract (§8.3)
A DIA-NN-shaped report with: `Run, Protein.Group, PG.MaxLFQ, Q.Value, Lib.Q.Value,
Lib.PG.Q.Value` (+ optional `Empirical.Quality, PG.MaxLFQ.Quality, Genes,
Protein.Names`). `run_search.py` produces this for non-DIA-NN engines.

## Design matrix
`~ 0 + groups [+ Batch + Covariate1 + Covariate2]`, colnames = group levels.
**Rank-checked before fitting** (`qr(design)$rank`); fails on confounded covariates
or empty groups. Groups with <2 replicates have no within-group variance — warn the
user at the design step (`collect_conditions.py --validate` flags singletons).

### Continuous covariates — `--numeric-covariates Covariate1[,Covariate2]` (2.12)
Every covariate column is a **factor** unless named here: one design column per level. A
continuous measurement (a centred log2 bait intensity, an age, a run depth) as a factor gets a
level per distinct value — a staff analysis's 11-value TurboID intensity became an 11-level
factor and "not full rank". Named in `--numeric-covariates`, it is **one numeric column**: one
slope, shared by all groups. Ask the user which covariates are continuous; never decide it from
the values (a batch numbered 1/2/3 is a factor, and `Batch` always is — only `Covariate1` /
`Covariate2` can be numeric). When a design fails and a covariate holds only numbers with 3 or
more distinct values, the error names the flag.
- **Values** are used as given, not centred, and must be a finite number for every sample: a
  blank, `NA`, text or `Inf` stops the run before quantification, naming the samples (fill them
  in — never a made-up value — or fit the column as a factor). A column with one value stops too.
- **The contrasts stay between groups.** With `~ 0 + groups + Covariate1` the design's columns
  are the group levels plus `Covariate1`; `makeContrasts("B-A", levels = design)` puts weight 0
  on `Covariate1`, so each contrast is the difference between the group coefficients — the
  groups compared at the same covariate value (ANCOVA's parallel slopes). Centring the
  covariate would move the group coefficients, never their difference. The same holds for both
  `--method`s (dpc: `dpcDE`, whose `voomaLmFitWithImputation` takes any design matrix; maxlfq:
  `lmFit`), for a fixed `--block` (its columns are added to this design; a numeric covariate
  constant within each block is absorbed and dropped, as a factor would be) and for a random
  block or technical replicates (`duplicateCorrelation` and the blocked fit use the same
  design). No group × covariate interaction is fitted: a slope that differs between groups is
  not modelled.
- **Confounding is warned about, per contrast** (2.12). A comparison is flagged when, among the
  samples of the groups it compares, the groups explain ≥ 80% of the covariate's variance (R² of
  covariate ~ group within them; variance-inflation factor ≥ 5), when the design alone widens
  its standard error ≥ √5 (2.24x, the same line), or when the compared groups' covariate values
  do not overlap at all (the comparison is then an extrapolation along the slope). Per contrast,
  not overall: with three or more groups a pair the covariate separates hides in a modest
  overall R² (0.595 overall, 0.988 within the two compared groups, in the review's example). The
  run warns, naming each flagged comparison and why; `methods.txt` carries a `CAUTION`, the report
  a warning callout, the brief a Data Quality Note and AUDIT.md a WARN, all naming it: that
  comparison rests largely on the fitted straight line, and if the covariate differs between the
  groups *because of* the grouping, adjusting for it removes part of the effect being tested. Ask
  whether it belongs in the model.
- **Recorded** in `de_provenance.json` `covariates`: each term's `type` (`factor` with
  `n_levels`, or `numeric` with `source: --numeric-covariates`, `n_distinct`, `min`, `max`,
  `group_r2` (overall), `confounded`, `confounded_contrasts`, `se_inflation` — per contrast, the
  standard error with the covariate over without it, from the design alone — and `per_contrast`:
  for each contrast the groups compared, `r2` within them, `overlap`, `se_inflation`,
  `confounded` and `why`), `numeric_requested`, `absorbed_by_block`, `model`, the two thresholds
  and `warnings`. `methods.txt` has a `Covariates` line, the Methods paragraph and the
  analysis brief a sentence (`make_methods.de_covariate_sentence`), the analysis page a callout,
  and `reproducibility_log.R` writes the values as numbers and fits `as.numeric(...)`.

## Paired / repeated designs — `--block <column>`
When several samples come from one source — the IPs cut from one mouse brain, the
biopsies from one patient, before/after samples from one animal — they are correlated,
and fitting them as independent throws the pairing away. Put the unit in its own
`conditions.csv` column (e.g. `Mouse`) and pass `--block Mouse`.

**Fixed or random — `--block-effect auto|fixed|random` (default `auto`, recorded as
`block.effect` / `block.effect_choice`).**
- **Crossed + every contrast within one block → fixed.** Before/after in the same patient,
  treated vs control in the same donor: the block is crossed with the groups (it stays
  full rank as design columns) and every contrast compares samples of one subject. The
  fixed subject effect (`~ 0 + groups + Patient`) is the exact paired analysis. Simulated
  (6 patients × 2 conditions, blocksim/paired.R): fixed type I 0.050 and power 0.95 at
  every per-protein correlation; the random effect's type I drifts 0.078 → 0.010 and its
  power falls to 0.76 as the correlation goes 0.05 → 0.85, because one consensus
  correlation is applied to every protein.
- **Otherwise → random**, the way limma fits multi-level experiments: one consensus
  within-block correlation estimated across all proteins by `duplicateCorrelation()`, used
  by `lmFit(block =, correlation =)`. This is the nested case (mice within age, the mouse-brain IP cohort):
  a fixed mouse term would be aliased with the age groups and between-mouse contrasts need
  the random effect. `--block-effect fixed` on a nested block stops before quantification.
- **Nested in a covariate, not the groups** (each patient's pre/post pair run in one
  `Batch`): the fixed subject effect absorbs that covariate — its columns are sums of the
  subject's — so the fixed fit leaves it out (`block.absorbed_covariates`, a `Dropped from
  the design` line in `methods.txt`, the Methods sentence, the repro script). Every message
  names what the block is aliased with: the groups, or covariate X.

- **dpc**: `limpa::dpcDE(y, design, block = b)`. `dpcDE` passes `...` to
  `voomaLmFitWithImputation()`, which takes `block` natively: it estimates the
  correlation with the vooma precision weights, refits, recomputes the weights and
  re-estimates it (limpa prints "First/Final intra-block correlation"), then fits
  `lmFit(block, correlation, weights)` (limpa 1.4.0 source). **maxlfq**:
  `duplicateCorrelation(E, design, block)` → `lmFit(E, design, block, correlation)`.
- **Nested is fine, fixed-and-random is not.** The block may sit inside a fixed factor
  (mice 1–3 Old, 4–6 Young; groups `Old_JPH3 … Young_IgG`): within-mouse contrasts
  (bait vs IgG) gain power, between-mouse contrasts (Old vs Young) are judged on the
  number of mice — see `--block-scope` below. Do NOT also put the block in the design: a column named
  `Batch`/`Covariate1`/`Covariate2` is a fixed covariate, so `--block Covariate1` stops
  with an error — rename the column (e.g. `Mouse`). `collect_conditions.py --map` keeps
  a Mouse / Animal / Subject / Patient / Donor column under its own name and reports it
  as `block_column` (other extra columns still become `Covariate1/2`). As a fixed
  covariate nested in the groups the subject makes the design rank-deficient; run_de.R
  then stops up front and points at `--block`.
- **`--block-scope within|all` — which contrasts the blocked fit reports (default
  `within`).** One run fits both models on the one quantification (the slow part is not
  repeated) and reports each contrast from its fit:
  - `within`: a contrast **between** blocks that uses **at most one sample per block**
    (Old_JPH3 vs Young_JPH3: 3 mice vs 3 mice, one IP each) comes from the fit with
    samples independent; every other contrast — within-block, partial, and between-block
    contrasts that pool several samples per block — from the blocked fit. If any block
    holds two samples of one group (technical replicates of a mouse), everything comes
    from the blocked fit.
  - **A between-block contrast on the blocked fit is NOT a remedy, only the lesser evil.**
    Pooled over baits (Old vs Young over all IPs) it cannot use the independent fit —
    unblocked it is pseudo-replicated (simulated type I 0.07–0.36) — but the blocked fit
    is itself anti-conservative for proteins with strong mouse effects (type I 0.10 / 0.22
    at per-protein correlation 0.6 / 0.85, nominal 0.05; blocksim/between.R). Every such
    contrast gets a `block.warnings` entry and a `CAUTION` in `methods.txt`. Define
    between-block contrasts one sample per block (per bait) instead. (A block-level
    analysis — one score per mouse, then a two-sample test — cut it to 0.08 at 0.85 in
    simulation; a follow-up.)
  - `all`: every contrast from the blocked fit.

  **Why `within` is the default** (the mouse-brain IP cohort: 6 mice × 5 IPs, consensus correlation
  0.17). A one-sample-per-mouse age contrast compares independent samples, and in a
  balanced design the independent fit's variance is unbiased for every protein. It is not
  exact: the pooled residuals come from the same mice across groups, so they share the mouse
  effect and their degrees of freedom are overstated — mildly (simulated type I 0.047–0.071
  across per-protein correlation 0.05–0.85; a fit on the two groups' samples alone is
  exact, a follow-up). The blocked fit
  applies ONE consensus correlation to all proteins, so it understates the between-mouse
  variance for proteins with strong mouse-to-mouse variation: the blocked/independent SE
  ratio on the age contrasts was 1.04 at per-protein correlation ≤ 0 and 0.86 at > 0.6
  (ideal: 1), matching what the design predicts to within 0.015. The 77 age calls only
  the blocked fit made were those proteins (median per-protein correlation 0.44 vs 0.14
  overall — red-cell, complement, tRNA-synthetase proteins that vary animal to animal).
  On the bait-vs-IgG contrasts the blocked fit is the right model: +14–58% calls, none
  lost, top hits and fold changes unchanged.

  Recorded: `block.effect` (+ `effect_choice`), `block.scope`, `block.contrast_model` (`{"<contrast>": "blocked" |
  "independent"}`, beside `contrast_structure`), `block.contrast_model_rule`;
  `de_tables` (each `DE_*.csv`, its model, its significant count); the `Blocking` lines
  of `methods.txt`; the console line of each contrast (`[blocked fit]`); the Methods
  sentence. `block.applied` stays true while the blocked fit reports ≥ 1 contrast, and
  `block_column` is present exactly then. The DE-LIMP session's `fit` holds each contrast's
  reporting fit column by column (`fit$contrast_model`; moderated F dropped) with
  `fit_independent` beside it; `reproducibility_log.R` fits both and picks per contrast.
- **Stops** (before quantification): column missing or blank for a sample, the column
  is `Group`/`File.Name`/a covariate, the block is encoded in the design (its levels
  coincide with the groups), or a block holds one analysed sample.
- **Warns** (`CAUTION` in `methods.txt`, `block.warnings` in `de_provenance.json`): a
  consensus correlation ≤ 0 (blocking gains nothing — check the assignments), fewer than
  50 proteins with an estimate, the two limpa passes disagreeing by > 0.1, or only 2
  blocks.
- **Recorded**: `de_provenance.json` `block` — column, levels and sizes, consensus
  correlation (and limpa's first-pass one), proteins estimated, per-protein quartiles,
  estimator, fit, and each contrast labelled `within` / `between` / `partial`; the
  `de_engine` label gains `; block = Mouse (…)`; `methods.txt` has a `Blocking` line (an
  unblocked run says `none -- samples modelled as independent`); `make_methods.py`
  writes the sentence from the record; `reproducibility_log.R` refits it.
- With no `--block`, `run_de.R` prints a note when a metadata column recurs across
  groups (the same Mouse in several conditions) — the hint to ask the user.

### Technical replicates — the `Sample` column (automatic)
Runs of ONE biological sample (`core_submission.py locate --reinjections all` keeps every
injection) are technical replicates. `conditions` / `collect_conditions.py` (`--mapping-json`
key `technical_replicates`) write each run's sample into the column **`Sample`**
(`TECH_REPLICATE_COLUMN`; the name is the declaration, so a CSV copied without its
`.decisions.json` still carries it). When a `Sample` value is shared by several runs, `run_de.R`
blocks on it by itself: the sample is a **random** blocking factor (duplicateCorrelation),
**every** contrast is reported from that fit (scope `all` -- the samples are nested in the
groups, so the independent fit would count each injection), and a sample injected once is a
block of one. It stops instead when a `Sample` has runs in two groups, or with `--block <other>`,
`--block-effect fixed` or `--block-scope within`. `de_provenance.json` `block.technical_replicates`,
methods.txt (`Replicates :`) and the Methods sentence name the samples. Averaging the injections
first is not offered: DPC-Quant's per-run standard errors and observation counts would have to
be combined, not just its means. To keep one injection per sample instead: `locate
--reinjections latest`. (2.10 review, MED 7: before this, two injections counted as n = 2.)

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

## Normalised or non-normalised quantities: the experiment type and the data check

Global normalisation — DIA-NN's, and the maxlfq path's quantile step — assumes most proteins
do not change between samples. An IP / pull-down breaks that: its IgG or bead controls carry
little protein, normalisation scales them up several-fold, and the enrichment is hidden or
turned into depletion. So the DE's input is a decision, made once and recorded:

1. **The experiment type** (`experiment_type.py`, step 3): proposed from the CoreOmics
   submission's experiment type and confirmed with the user, or asked; never taken from file
   names. Its default (`experiment_type.DEFAULTS`):

   | type | default quantities | why |
   |---|---|---|
   | whole proteome | normalised | most proteins do not change |
   | IP / AP-MS / pull-down (IgG or bead controls) | **non-normalised** | the controls carry little protein and are not loaded by amount |
   | proximity labelling (TurboID / BioID / APEX) | normalised | the controls carry plenty of signal; samples usually loaded by equal protein |
   | fractionated for depth (offline high-pH, fractions combined per sample) | normalised | every sample is still the whole proteome |
   | separation / profiling fractions COMPARED (SEC, density or sucrose gradient, BN-PAGE complexome, organelle or subcellular) | **non-normalised** | each fraction holds a different part of the proteome. DIA-NN's README (FAQ "What is normalisation and how does it work?"): "If this condition is not satisfied (e.g. when analysing fractions obtained with some separation technique, like SEC, for instance), then normalisation should not be used"; and "Disabling normalisation" lists "any kind of protein fractionation" |
   | secretome, other | normalised | leaning on the data check |

   The submitter's "Normalization" answer (how the samples were loaded) is set beside it; a
   contradiction (normalised quantities for samples normalised by volume; non-normalised for
   eluates loaded by protein amount) is a question for the user.
2. **The quantities** (`run_de.R --quantities`), verified in DIA-NN 2.7.0's README and a 2.7.0
   report: `Precursor.Quantity` is non-normalised, `Precursor.Normalised` normalised, and
   `normalised = Normalisation.Factor × non-normalised`; `PG.MaxLFQ` is normalised unless the
   search ran with `--no-norm`. **dpc** reads either precursor column through limpa's own
   quantity-column argument (limpa ≥ 1.4 `intensity.column`, 1.2.x `qty.column`; limpa adds no
   between-run normalisation of its own) — the cleanest route, with no re-search. **maxlfq**
   `raw` needs a `--no-norm` report (checked from the data: `Precursor.Normalised` equals
   `Precursor.Quantity` on every row) and then applies no quantile step; a normalised report is
   refused, never relabelled. A pull-down on maxlfq used to be normalised twice. **maxlfq
   `normalised`** says what its `PG.MaxLFQ` carries before the quantile step — DIA-NN's
   normalisation, or none for a `--no-norm` search — from the search's provenance
   (`normalization_check.searched_no_norm`: the chain's `no_norm_report.parquet` name, else
   `search_provenance.json`), reading the data only when that records nothing, and NOT RECORDED
   when neither can tell (`de_provenance.json` `dia_nn_normalisation`).
3. **The data check** (`normalization_check.py`), for EVERY type, after the DE has been run
   both ways. Thresholds, each justified in the script and set on the synthetic whole-proteome,
   IP and TurboID cases of `tests/test_normalization_check.py`:
   - **B1 factors vs groups** — the deciding signal: the per-run factor (median log2
     `Precursor.Normalised / Precursor.Quantity`) differs ≥ 2-fold between group medians AND the
     groups separate completely (a permutation p is reported; with 3 vs 3 runs its floor is 0.1,
     which is why separation, not p, decides). **With n = 2 per group it trips more**: two runs
     per group separate by chance a third of the time, and the 2.10 review's simulation of
     whole-proteome 2 vs 2 cohorts tripped B1 in 0.1% / 4% / 13% at a per-run loading SD of
     0.3 / 0.5 / 0.7 log2. A trip only asks, so B1 stays;
   - **B2 shift before normalisation** — a simple contrast's median per-protein log2 shift in
     `Precursor.Quantity` ≥ 1 with ≥ 80% of proteins moving the same way;
   - **B3 identifications** — group medians of precursors per run ≥ 1.5-fold apart;
   - **B4** — `Normalisation.Instability` above 0.3 in any run (see `normalisation_stability`);
   - **volcano**, per contrast, from the DE tables (never an image), both versions: the
     non-significant proteins' median log2FC ≥ 0.5 off zero; more than 5:1 up:down (or
     down:up) among ≥ 20 significant — for an enrichment-vs-control contrast only many proteins
     DOWN counts (background made to look depleted); over half the proteins significant; the
     bait outside the 10 most enriched or not significant; a p-value histogram with an excess
     near 1. In the non-normalised version an enrichment contrast's background sitting above
     zero, and much of it significant, is expected and not counted.

   Separation / profiling fractions are compared without normalisation by design, so their
   non-normalised volcano's shift, lopsidedness and % significant are reported as expected, not
   counted. With a **normalised** default, B1-B4 and the normalised volcano trip it. With a
   **non-normalised** default (IP), it trips on factors that span ≥ 2-fold across runs WITHOUT
   following the groups (loading variation the raw quantities would keep) and on the raw
   volcano. **A trip stops the analysis before the final DE** — and so does a session with no
   experiment type recorded, where nothing is decided for the user: `NORMALIZATION_CHECK.md`
   gives the side-by-side summary and a recommendation (the type's default unless that
   default's own volcano is clearly the worse one, by 2 or more anomalies; the trips are named as
   questions), a Core staff member chooses — never Claude (`decide --by <their HIVE login>
   --reason`, the login checked against the Core's staff list by `staff.py`) — and the final `run_de.R
   --normalization-check` refuses any other quantities, and any report or conditions other than
   the ones checked (sha256 recorded). It is held to the checked model too (2.12): run_de.R
   records the design it was asked for before fitting (`design_request`: method, formula,
   numeric covariates, block column/effect/scope, contrasts, the flags given) and the conditions'
   sha256; `normalization_check.py run` refuses two check DEs whose designs differ or that read
   another report or conditions file, and records `check_design`; the final DE passes its design
   to `gate` (`--design-json`), which refuses another method, formula, numeric covariate or block
   and warns on other contrasts (`normalization_check.design_match` in de_provenance.json). A
   check from before 2.12 has no `check_design`: the final DE runs, says the design was NOT
   CHECKED, and the Methods line carries `[confirm]`. The step-8 commands run_de.R prints carry
   its design flags. No silent switch in either direction. In a session that
   records a type, `run_de.R` refuses a DE with neither `--normalization-check` nor
   `--check-input`; a DE without a decided check is a WARN in AUDIT.md and holds
   `core_submission.py deliver` (the second gate beside the QC one). A DE from before 2.10 (no
   check record) is let through only by `normalization_check.py ack-legacy --session S --by
   <login> --note "<why it stands>"`, tied to that DE record. MaxLFQ's non-normalised DE needs
   a `--no-norm` report, given to the check as `--report-raw` (dpc reads both quantities from
   one report). Proximity labelling also reports how steady the endogenously
   biotinylated carboxylases are under each option (PC, PCCA, MCCC1, ACACA, ACACB — UniProt
   KW-0092, human and mouse): a reference, never a normaliser.
4. **The record**: `de_provenance.json` `normalization_check` (`schema_version` 1, fixed keys:
   the type and its source, the default, every check's numbers, the volcano metrics per
   contrast and version, the carboxylases, the trips, the decision — by role, "Core staff", with
   `statement` the sentence the Methods and AUDIT.md quote; the name and login are only in the
   staff-only `logs/normalization_check.staff.json`), written
   whether or not anything trips so projects can be compared over time; methods.txt's
   `Quantities` / `Normalization` / `Norm. check` lines; the Methods' normalisation sentence;
   AUDIT.md's `normalization_check`. A DE run without a check says NOT RUN.

## Contaminant filter (both paths, on by default)

`run_de.R` drops every precursor that maps to a `Cont_` (or FragPipe `contam_`) entry
before quantification. That is any accession in `Protein.Ids` (`Protein.Group` for an adapted
report), the rule of DIA-NN's own `--cont-quant-exclude`: a peptide shared between a sample
protein and a contaminant entry can carry the contaminant's signal.
- **Where the rule came from is the pipeline's to say.** The record's `rule` text (which
  methods.txt and `make_methods.py` quote) ends with the descriptor's `contaminant_rule_origin`
  (`contaminants.R`): "the rule of DIA-NN's --cont-quant-exclude, applied here" for a DIA-NN
  report (dpc, and maxlfq on DIA-NN's PG.MaxLFQ), "applied here, in the DE step, to the Sage
  report" (or the FragPipe / AlphaDIA / Radiant label the report declares) for an adapted one.
  Before 2.10 every engine's Methods named DIA-NN's flag, which a reader could take to mean
  DIA-NN had run.
- **Every engine's report.** The rule tests bare accessions, as DIA-NN writes them, so every
  adapter writes those through `protein_ids.py`: `group_accessions` for Sage's `proteins` and
  Radiant's groups, and `fragpipe_group` for FragPipe DDA: its `Protein ID` (keeping `contam_`),
  then every protein under `Indistinguishable Proteins`, so a `Cont_` entry listed only there is
  seen, as any member of a DIA-NN group is (before 2.10 only the leading protein was read).
- **Before 2.9 Sage slipped past.** `adapt_sage` wrote Sage's full FASTA IDs
  (`sp|Cont_P02769|ALBU_BOVIN`), so no Sage DE ever had a contaminant removed. A Core member's HeL50
  UnvPe tested 171 `Cont_` groups, while methods.txt said "none".
- **Counts are distinct items in the report's own unit** (`contaminant_unit` in contaminants.R,
  from the declared quantity level): DIA-NN precursors, Sage peptides (`Peptide.Id`, written by
  `adapt_sage`), a protein-level adapter's protein groups. Up to b1b94f5 an adapted report
  counted item × run rows as "precursors": a Sage DE said 2,744 for 686 peptides × 4 runs.
- **Mixed groups.** A removed group that also names a sample accession is counted
  (`n_protein_groups_mixed`). For a peptide-level report, the sample proteins it named are
  counted per accession: those that keep other peptides, and those left with none. A
  protein-level report cannot show that, and its record says "not computed" rather than 0.
- **A keratin sample's kept keratins** are then described per engine, by the pipeline's own
  descriptor (`kept_contaminant_quant`):
  - DIA-NN: what its `--cont-quant-exclude` left out of its own quantities, stated only when
    the run recorded the flag. The source is the command line in the DIA-NN log beside the report,
    else the parameters file the search ran with (`make_methods.diann_cont_quant_exclude`, the
    one reader). DIA-NN logs a `--cfg <file>` unexpanded, so the reader opens that file where
    `--cfg` stands. A relative path resolves where DIA-NN ran. A `--cfg` that can't be read
    never counts as "not set". Otherwise the text says `[not recorded -- confirm]`. FragPipe's
    DIA-NN step may not have been given the flag.
  - Sage: the highest peptide is taken here from the rows kept.
  - FragPipe / AlphaDIA / Radiant: the engine's own quantities, used as reported. Whether
    keratin is under-quantified there is "not known" (JSON null), because the engine's protein
    inference (Philosopher's razor peptides) can move shared peptides.
  - None of the non-DIA-NN routes is told that DIA-NN did anything. The rule, the tag and the counting
live in `scripts/contaminants.R` (the tag mirrors `fetch_fasta.py`'s `CONT_TAG`; a test
asserts it).

Why it exists: DIA-NN's `--cont-quant-exclude` only shapes DIA-NN's OWN quantities, and
`run_de.R` re-quantifies from the report. On the mouse-brain IP cohort (2026-09-24)
all 121 `Cont_` groups — bovine serum proteins from the antibody prep — went into the DE
and came out as hits, while the Methods said they were excluded. With the filter the same
report loses 2,196 precursors: the 121 `Cont_` groups, plus the shared peptides of 47
sample groups (Eno1, Aldoa, Tubb3 …); no sample group disappears.

What it records, so nothing downstream restates the policy (architectural rules 1 and 2):
- `de_provenance.json` → `contaminants`: `policy` (`removed` / `kept` / `none_present` /
  `not_checked`), the rule and column, precursor and group counts, the share summary, and
  the database check. `filters_applied` lists it with the other filters.
- `methods.txt` → a `Contaminants :` block written from that record; `make_methods.py`
  words the DE paragraph from the same record, and the DIA-NN sentence from the search
  parameters (never from the FASTA sidecar's recommendation).
- `contaminants_removed.csv` — every protein group that lost precursors, and whether it
  lost all of them. `audit_results.py` reads it, so a real protein removed as `Cont_`
  is still named.
- `QC_contaminant_share.csv` — per run, the contaminant share of `Precursor.Quantity`
  (measured signal; `Precursor.Normalised` only if the report lacks it — DIA-NN's
  RT-dependent normalisation does not keep a run's signal fractions: 12.7% vs 28.3%
  median on the same report). Written with `--keep-contaminants` too.
- `reproducibility_log.R` removes the same rows.

`--fasta-meta` (default `./search.fasta.meta.json` if present): a database built before
`fetch_fasta.py` removed contaminant entries identical to target proteins (a sidecar with
no `contaminant_target_rule`), one built by the identity rule alone (the rule but no
`min_unique_peptides`: fetch_fasta.py before skill 2.8.0, which kept near-identical entries
such as bovine EEF1A1 and YWHAZ for mouse), or one listing
`contaminants_identical_to_target_kept`, holds real proteins only as `Cont_` groups — the
filter removes those too. The run then prints a `CAUTION` and records `database_risk: true`
with a note naming the set; `audit_results.py --fasta-meta` re-checks the searched FASTA and
names the proteins. With NO sidecar the run checks the FASTA the search itself names (its
`search_provenance.json`, or `--fasta` in `report.log.txt`): `fetch_fasta.py check-db`
re-checks it, or — when it cannot be read here — recognises a known superseded database by
md5 or name (`superseded_databases.json`: the Core's Sep-2025 MRS human FASTA, 161 entries);
only when neither works is the check recorded as not run. The states are `sidecar_state()` in `fetch_fasta.py`, mirrored in
`contaminants.R` (a test keeps them equal). Rebuild the FASTA with this release's
fetch_fasta.py (skill 2.8.0 or later) and re-search; the Core's shared human+contaminant FASTA
was rebuilt with it on 2026-09-25 (`MRS/UP000005640_9606_plus_universal_contam_2026-09.fasta`,
docs/HPC_PATHS.md).
`--keep-contaminants` keeps every `Cont_` group in the DE instead (the true contaminants
are then tested too).

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

`run_de.R` also writes **`reproducibility_log.R`** (via `repro_script.R`): the same
analysis emitted as flat, literal R — the report path, the q-cutoff and the q-columns
it was actually applied to, any QuantUMS pre-filter, the sample→group map, the
covariates, the blocking factor (if any), the design, the contrasts. Runnable with `Rscript`, needing only R and
limpa/limma. It is built from the objects that ran, for the same reason `methods.txt`
is: a hand-written recipe drifts, a generated one can't. Point users at it whenever
they ask what was done or want the code. → `references/reproducibility.md`.

## Top-protein violins (`make_figures.R` → `violin_top_<contrast>.png`)
One figure per contrast: the top `--violin-top` proteins (default 8), ranked by
`adj.P.Val`, ties broken by the raw `P.Value` (the DE table's own order), as small
multiples. Each panel shows the contrast's groups (reference group on the left), one point
per run, the group mean as a bar, and the model's log2FC + adjusted p (`n.s.` when not
significant; non-significant panels are shown as context only). A group with 5 or more runs
also gets a violin; with 4 or fewer only the points are drawn, since a density through 3
points is a shape the data do not have. The arrow is the **model's** log2 fold change, drawn
from the reference mean (red = higher, blue = lower), so it always agrees in sign with the
label; the difference of plain means can disagree once the model weights runs, models
inferred values or removes a block effect. For a within-block contrast (`--block`, recorded
as `within` in `de_provenance.json`) grey lines join each block's runs: the paired changes
the model tested. Other groups are left out on purpose: `heatmap_top.png` already shows the
top proteins across every group.

**Each point is marked measured or not** (DE-LIMP's expression-grid violin):
filled teal = at least one precursor observed in that run; hollow amber = no precursor
observed. Read it from `Detection_Matrix.csv` (Protein.Group + one column per run, named
as in `Expression_Matrix.csv`; values = precursors observed, 0 = not observed). What a 0
means comes from `de_provenance.json` (`detection_matrix.zero_means` if present, else
`pipeline_id`): **Inferred** under DPC-Quant (the value exists but was modelled), **Missing**
under MaxLFQ (no value exists, so nothing is drawn). A label under a group counts its
measured runs, shown only when not every run was measured.

**Say this when you describe the figure:** a group labelled *all inferred* (0 measured
runs) makes that protein's fold change a **detection event**, since it is seen in one group and
not the other, not a measured magnitude. The subtitle names every such protein; confirm
them before building on the size of the change. With no `Detection_Matrix.csv` (older
runs) the points carry no status, the subtitle says detection status was not recorded, and
no status legend is drawn. Never describe those points as measured.

## Citations (verified June 2026)
- **limpa / DPC:** Li M, Cobbold SA, Smyth GK (2025) bioRxiv 10.1101/2025.04.28.651125;
  Li M, Smyth GK (2023) Bioinformatics 39(5):btad200. (DE-LIMP's
  `dpc_pipeline_descriptor()` mis-cites this as "Law CW, Smyth GK" — fix upstream.)
- **MaxLFQ path:** DIA-NN MaxLFQ (Demichev et al. 2020, Nat Methods 17:41) +
  limma (Ritchie et al. 2015, NAR 43:e47).
- **Sage (highest-peptide path):** Lazear MR (2023) J Proteome Res 22(11):3652-3659,
  doi:10.1021/acs.jproteome.3c00486 (Crossref-verified 2026-09-29) + limma.
- **QuantUMS quality filtering:** da Cruz Moschem J, Silva Campitelli de Barros BC,
  de Toledo Serrano SM, Chaves AFA (2025) *Decoding the Impact of Isolation Window
  Selection and QuantUMS Filtering in DIA-NN for DIA Quantification of Peptides and
  Proteins.* J Proteome Res 24:3860-3873. doi:10.1021/acs.jproteome.5c00009.
  **VERIFIED 2026-08-04** — an earlier note here said this "could not be verified";
  that was wrong. DE-LIMP cites it as "Moschem et al." (the first author's surname is
  da Cruz Moschem). QuantUMS computes three scores: protein-group MaxLFQ quality,
  empirical quality, and quantity quality, all measuring MS1/MS2 feature agreement.
  The skill filters on the first two (`--pgq-cutoff`, `--eq-cutoff`).

## Less-used flags
- `make_figures.R --top <n>`: how many proteins the top-protein heatmap shows (default 50).
- `sample_quality.py --z <x>`: the |z| above which a sample is flagged as elevated (default 1.5).

## Step 8: the normalisation check and the final DE

**Then check the normalisation before the final DE** — for every experiment type. Run the DE
both ways into scratch folders (one job on HIVE), then the check. **A DE, set-test or figure job
on HIVE** starts from `python3 ~/proteomics-pipeline/scripts/on_compute.py --header <S>/logs
--name <short name> --time <minutes> --mem <size>G --cpus <n>`: the queue the skill's rule
picks, `--time`/`--mem`/`--cpus` explicit, one log per job (`%j`), `set -euo pipefail`. Append
the commands, `sbatch` it, and when the `ucd-hpc` skill is installed here, lint it first with
its `lint-jobscript.sh` (`U` = that skill's scripts folder, `~/.claude/skills/ucd-hpc/scripts`):
`bash scripts/hive_exec.sh 'bash -s -- <job file>' < "$U"/lint-jobscript.sh`
(it misreads a log path that holds a space as having no `%j`; sbatch reads it right, so that
one message can be ignored).

```
Rscript scripts/run_de.R --input ./search_out/report.parquet --metadata conditions.csv \
    --method <dpc|maxlfq> --fasta-meta ./search.fasta.meta.json [--block Mouse] [--numeric-covariates Covariate1] \
    --quantities normalised --check-input --outdir <S>/output/norm_check/normalised
Rscript scripts/run_de.R ... (the same) --quantities raw --check-input --outdir <S>/output/norm_check/raw
python3 scripts/normalization_check.py run --report ./search_out/report.parquet \
    --conditions conditions.csv --de-normalised <S>/output/norm_check/normalised \
    --de-raw <S>/output/norm_check/raw --session <S>
```

`--quantities normalised` reads DIA-NN's normalised quantities (dpc: `Precursor.Normalised`;
maxlfq: `PG.MaxLFQ` + quantile normalisation); `raw` reads DIA-NN's non-normalised ones (dpc:
`Precursor.Quantity`, through limpa's own quantity-column argument; maxlfq: only a report
searched with `--no-norm` — `diann_parallel.py --no-norm`, or `--no-norm` in a single-shot
cfg — with no quantile step; a normalised report is refused, since its `PG.MaxLFQ` is
normalised). **maxlfq on a normal report: its `--quantities raw` DE is refused, so do one of
these.** (1) Run all three DEs (both check inputs and the final DE) with `--method dpc`, which
reads both quantities from this report, with no new search. (2) To keep MaxLFQ, search again
with `--no-norm` (`diann_parallel.py --no-norm` writes `no_norm_report.parquet`), run the raw
check DE on that report, and pass it to the check as `--report-raw
<no_norm_report.parquet>`. B1–B4 still read `--report`, and a raw decision holds the final DE
to the `--report-raw` report. The check reads the normalisation factors (`Precursor.Normalised /
Precursor.Quantity` = DIA-NN's `Normalisation.Factor`) against the groups, the fold changes
before normalisation, identifications per group, `Normalisation.Instability`, and each
contrast's volcano both ways (centre, up:down, % significant, bait rank, p-value shape); for
proximity labelling also how steady the biotinylated carboxylases are (a reference, never a
normaliser). `--check-input` marks these two DEs as the check's inputs: in a session that
records an experiment type, `run_de.R` refuses any DE that has neither it nor
`--normalization-check`, and prints these exact commands (`normalization_check.py required
--input <report> --metadata <conditions>` prints them too, and the chain's recovery step in
`RECOVERY.md` / `checkpoint.py status` is this check, never the final DE alone). **Exit 0:** the
data agree with the type's default — use it. **Exit 3 (TRIPPED, or no experiment type
recorded): STOP before the final DE.** With no type recorded nothing is decided for the user:
ask what the experiment is (step 3) and which quantities to use. Show the Core staff member
`NORMALIZATION_CHECK.md` (side by side, normalised vs non-normalised, with a recommendation),
ask which to use, and record their answer. **The decision is the person's, never Claude's.**
Never choose, and never fill in `--by` yourself. `--by` is their HIVE login, never a full name.
`decide` checks it against the Core's staff list (`staff.py`: `CORE_STAFF_FILE`, default
`/quobyte/proteomics-grp/.config/core_staff.txt`, trusted only when a Core admin owns it). With
the list, a login on it is accepted and any other is refused. With no list, a login that is an
agent's name (claude, assistant, bot…) is refused, and any other passes with a NOT VERIFIED
warning. Their name and login go only into the
staff-only `logs/normalization_check.staff.json`. The Methods, AUDIT.md, `de_provenance.json` and
the zip say "Core staff chose … because <reason>", so the reason must not name anyone. A
second `decide` replaces the first: its sentence has no history, and the report, the audit and
the podcast state only the decision that stands -- never the earlier reason, nor that it was
"corrected" (the change and why go in `logs/decisions.md`; both decisions stay in the staff-only
record):

```
python3 scripts/normalization_check.py decide --check <S>/output/norm_check/normalization_check.json \
    --quantities <raw|normalised> --by <their HIVE login> --reason '<why, in their words; no quote marks, $ or backticks>'
```

The recommendation keeps the type's default unless that default's own volcano is clearly the
worse one; the checks that tripped are the questions to ask. Never switch silently, either way.
The final DE is held to the decision AND to the data it was made on (it refuses any other
`--quantities`, and any other report or conditions — their sha256 is in the check) AND to the
model: all three DEs carry the same design flags (`--numeric-covariates`, `--block`,
`--block-effect`, `--block-scope`, `--contrasts`). `normalization_check.py run` refuses two check
DEs fitted with different designs or on other files and records their design (`check_design`);
the final DE refuses another method, formula, numeric covariate or block (other contrasts: a
warning; a check from before 2.12 records no design, said as NOT CHECKED — re-run the check to
check it). It records all of it in `de_provenance.json` `normalization_check` and the Methods.
A DE without a decided check is a WARN in AUDIT.md, and `core_submission.py deliver` holds the delivery until step 8 is done
(a delivery meant to go without a DE — a QC report — is `deliver --mode no-de`, step 12b).
A DE from **before 2.10** has no check record at all. Re-run step 8 on it, or, if Core staff
decide it stands as it is, they record that with `python3 scripts/normalization_check.py
ack-legacy --session <S> --by <their HIVE login> --note "<why this DE stands>"`. This is the
same shape as the QC ack and goes through the same staff-list check. It is accepted only for a
DE with no check record and is tied to that DE record. It is written to the staff-only
`logs/normalization_check.staff.json`, and the gate then passes. The final DE:

```
Rscript scripts/run_de.R --input ./search_out/report.parquet \
    --metadata conditions.csv --method <dpc|maxlfq> --outdir ./de_results \
    --fasta-meta ./search.fasta.meta.json [--block Mouse] [--numeric-covariates Covariate1] \
    --quantities <decided> --normalization-check <S>/output/norm_check/normalization_check.json
```

## Step 8: blocks and numeric covariates

**Several samples from one animal / patient → `--block <column>`.** IPs cut from the same
mouse brain, or before/after samples from one subject, are correlated; without `--block`
they are fitted as independent and the pairing is lost. Give the unit its own
`conditions.csv` column (`Mouse` — not Batch/Covariate1/2) and pass `--block Mouse`.
`run_de.R` picks the model (`--block-effect auto`, recorded as `block.effect`):

- **Crossed paired design** (before/after in the same patient; every contrast within one
  subject): a **fixed** subject effect (`~ 0 + groups + Patient`) — the exact paired
  analysis (type I 0.050 and power 0.95 at every per-protein correlation, where the
  random effect drifts 0.078 → 0.010 and loses power to 0.76).
- **Nested / multi-level** (mice within age, as the mouse-brain IP cohort): a **random** effect via limma
  `duplicateCorrelation` (limpa's `dpcDE` takes it natively). `--block-scope within` (the
  default) reports within-mouse contrasts (bait vs IgG) from the blocked fit and
  between-mouse contrasts that use one sample per mouse (Old vs Young for one bait) from
  the fit with samples independent — one run, one record (`block.contrast_model`). The
  single consensus correlation understates between-mouse variance for proteins that vary
  strongly animal to animal (the mouse-brain IP cohort: SE 0.86× where it should be 1×), so a
  between-mouse contrast that must use the blocked fit (pooled over baits, or
  `--block-scope all`) carries a `CAUTION` — define age contrasts per bait instead.
State n per group in mice, not IPs. run_de.R prints a note when a column looks like a
block and `--block` is missing — ask the user.
A `Sample` column whose value repeats (technical replicates from `locate --reinjections all`)
is blocked by run_de.R itself — random effect, every contrast from that fit — so pass no
`--block` for it; see "Technical replicates" above.
**A covariate that is a measurement** (a bait intensity, an age, a run depth) → ask, then pass
`--numeric-covariates Covariate1` to all three DEs (both check inputs and the final one): it is
fitted as one slope shared by all groups, not a factor with a level per value. Never decided
from the values; `Batch` is always a factor. A blank or non-number stops the run; a covariate
that nearly matches the groups of a comparison is a `CAUTION`, naming it, to raise with the user →
"Continuous covariates" above.

## Step 8: keratin samples and contaminants

**Keratin sample (step 3):** `run_de.R` reads it from the sidecar or `search_provenance.json`
and keeps keratin as the analyte — nothing to pass. For a search that records neither (an older
search, one run outside the skill), pass `--keratin-sample`: precursors that map only to
keratin-family contaminant entries are then kept (the entries are found with `fetch_fasta.py
keratin-db`) — another species' keratin (mouse fur, sheep wool in human hair) only when the
precursor also names a sample protein — trypsin, BSA and every other contaminant entry are still
removed, and the run prints a `CAUTION` to rebuild the FASTA with `--keratin-sample` and re-search.
With `--method maxlfq` the kept protein groups stay, but their quantities are DIA-NN's PG.MaxLFQ,
which already left the shared peptides out, so only a rebuild fixes those numbers. `contam_`
entries (a database built with FragPipe/Philosopher's `--contam`, not by the skill) count as
contaminants the same way wherever the tag survives — FragPipe's DDA adapter keeps it. FragPipe's
DIA route drops it (library.tsv), so `run_search.py` refuses a FragPipe search on such a database:
build it with `fetch_fasta.py`. When it cannot
tell whether the samples are keratin, it removes keratin-family entries as for any sample and
says so in its log and in `de_provenance.json` (`contaminants.keratin_sample`), which
`methods.txt`, `make_methods.py` and the report read.
**Contaminants are removed before quantification, on both methods.** Every precursor that
maps to a `Cont_` entry, or to a `contam_` one (a FragPipe/Philosopher `--contam` database, DDA
route) (any accession in `Protein.Ids` — the rule of DIA-NN's `--cont-quant-exclude`, applied by
`run_de.R` to every engine's report), is dropped before limpa/limma sees it. Every adapter writes
bare accessions as DIA-NN does (`protein_ids.py`): Sage's FASTA IDs (`sp|Cont_P02769|ALBU_BOVIN`)
used to slip past the rule, so a Sage DE tested its contaminants. DIA-NN's flag alone
never reached the DE, so contaminants used to be tested and came out as hits (bovine serum
HBB +10.5 log2 in antibody IPs). The counts land in `de_provenance.json` (`contaminants`),
`methods.txt` and `contaminants_removed.csv`; `QC_contaminant_share.csv` gives each run's
contaminant share of the precursor signal — report a high or group-confounded share in the
Data Quality Notes. `--keep-contaminants` opts out (only when the user asks). Pass
`--fasta-meta`: a database built before target-identical contaminants were removed holds
real proteins (ACTB, EEF1A1, keratins) only as `Cont_` entries, and one built by the identity
rule alone (fetch_fasta.py before 2.8.0) still holds near-identical ones (bovine EEF1A1 and
YWHAZ for mouse), so the filter removes those proteins too — the run then prints a `CAUTION`
naming them, AUDIT.md names them from the searched FASTA, and the fix is to rebuild the FASTA
with fetch_fasta.py from 2.8.0 on and re-search. Tell the collaborator which proteins.

## Step 8: the method, what `run_de.R` writes, QuantUMS

**`dpc` (limpa) is THE DEFAULT — use it unless the user asks otherwise or the data
cannot support it.** limpa models the detection-probability curve and quantifies from
precursor intensities directly, so it uses the whole measurement rather than a
pre-collapsed protein number. Don't switch away from it because a bundle happens to
say `maxlfq`; switch only for the two reasons below.

Use `maxlfq` when **either**:

1. **The user asks for it** — e.g. they want QuantUMS quality filtering (below), or
   to match a previous MaxLFQ-based analysis.
2. **The input cannot support limpa.** `readDIANN()` needs PRECURSOR-level columns
   (`Precursor.Id`, `Precursor.Normalised`). What matters is *which file you point at*,
   not which engine ran:

   | engine | file to use for `dpc` | verified |
   |---|---|---|
   | DIA-NN | `search_out/report.parquet` (native) | yes |
   | FragPipe | `dia-quant-output/report.tsv` + `--format tsv` | yes — its DIA route bundles DIA-NN, so this **is** a DIA-NN report, `Lib.Q.Value`/`Lib.PG.Q.Value` included |
   | Radiant | `radiant_to_delimp.py --out <x>.parquet` | yes |

   What does **not** work is the *adapted* `report.parquet` from `adapt_sage` /
   `adapt_fragpipe` / `adapt_radiant` — those deliberately collapse to one row per
   protein × run to feed the maxlfq path. `run_de.R` checks up front and names the
   precursor-level file to use instead.

Writes `DE_<method>_<contrast>.csv` + `Expression_Matrix.csv` +
`methods.txt` + `sessionInfo.txt` + `de_provenance.json` (exact R package versions) +
`QC_contaminant_share.csv` + `contaminants_removed.csv` + `Detection_Matrix.csv` (per protein ×
sample: precursors observed, 0 = inferred by DPC / missing for MaxLFQ, per
`de_provenance.json` `detection_matrix`) + **`reproducibility_log.R`** + (DPC runs)
`DE-LIMP_session.rds`, which loads straight into the DE-LIMP app
(https://delimp.stan-proteomics.org/) for interactive exploration.

`reproducibility_log.R` is the whole analysis as plain, flat R — every value written
out literally (report path, FDR cutoff, sample→group map, design, contrasts), runnable
with `Rscript` using nothing but R and limpa/limma. It is generated from the objects
that actually ran, so it can't drift from the CSVs next to it. **This is what to show
a user who asks "what did you do?" or "can I have the R code?"** — not the
`reproducibility/` bundle, which pins the environment and re-runs the multi-hour
search. Mention it when you hand over results.

**QuantUMS quality filtering — opt-in, `--method maxlfq`, DIA-NN input only.**
Only reach for this when the user asks. DIA-NN's
QuantUMS scores how well MS1 and MS2 agree for each measurement, and filtering on it
trades depth for quantitative reliability:

| flag | column | default | what it does |
|---|---|---|---|
| `--eq-cutoff` | `Empirical.Quality` | 0 (off) | drop precursors whose empirical quality is below this |
| `--pgq-cutoff` | `PG.MaxLFQ.Quality` | 0 (off) | drop protein groups whose MaxLFQ quality is below this |
| `--coverage-min` | — | 0.5 | require a protein in ≥ this fraction of samples before testing |

Both QuantUMS cutoffs are **off by default** — they discard real measurements, so
they are the user's call, not a silent default. Offer them when quantitative
precision matters more than depth (small fold-changes, few replicates), say roughly
how many proteins each costs, and record the choice: the counts land in
`filters_applied` and flow into `methods.txt`. They need DIA-NN input — the columns
don't exist in Sage/FragPipe/Radiant output, and the filter is skipped if absent.

`--coverage-min` is different: it is **on** at 0.5 because a protein seen in one or
two samples is an on/off observation, and letting such rows into `eBayes` destabilises
the variance moderation for *every* protein, not just those rows.

**After DE, log the contrasts** (`log_decision.py`) **and save the conversation**
(`save_transcript.py <S>`; hive_remote: `--hive "$S"`). See SKILL.md and outputs.md, "Keeping the record of the analysis".

## 8b. Figures, and what counts as significant (step 8b)

```
Rscript scripts/make_figures.R --de-dir ./de_results --conditions ./conditions.csv \
    --outdir ./figures
```

Produces publication-quality volcano (per contrast), PCA (each group circled and named),
heatmap of top proteins, per-contrast top-protein violins (`violin_top_<contrast>.png`: each
run's point filled if measured, hollow if inferred; a group never measured makes that fold
change a detection event), the detected-vs-inferred QC plot, and ONE panel of raw p-value
distributions (`qc_pvalue_panel.png`, for the appendix), plus `figures.json` (`figures`:
captions; `failed`: every figure not drawn, with the reason). The per-sample protein-count
plot is drawn only for a matrix with missing values (MaxLFQ); a DPC matrix is complete, so
every bar would be identical — it is skipped and the reason logged. The significance cutoff
comes from the DE run's `de_provenance.json`; `--adjp`/`--logfc` are only a fallback for a
run with no record, and the captions then say so. These get embedded in the report.

**Significance is `adj.P.Val` alone — there is no fold-change cutoff anywhere in this
pipeline.** `--logfc` only positions a labelled reference line on the volcano. Never
re-impose it when you count, rank, or describe significant proteins, and never justify a
result by "it's ≥2-fold". The moderated t-test tested H0: log2FC = 0, so that is the only
claim the FDR covers; filtering the list on the observed fold change afterwards would
report a stronger claim than the error rate supports, because an observed fold change is
a point estimate with no uncertainty attached. (Testing a fold-change threshold properly
needs `limma::treat()`'s interval null — McCarthy & Smyth 2009 — not a post-hoc cut.) You
may report how many significant proteins also exceed the reference line, as a description
of effect size.

**A protein measured in neither group of a contrast is not tested in it.** With three or more
groups, a protein seen only in other groups has `Evidence` = `not measured in either group`
in that contrast's table, no `logFC` or p-value, and is left out of that contrast's BH
correction (this file). The row stays. Count *tested* proteins as the rows
with a `P.Value` (`de_provenance.json` `de_tables` → `n_tested`), never as all rows, and never
describe such a protein as a presence call or a hit.
