# Differential expression reference

The DE step (`run_de.R` + `build_maxlfq.R`) is a faithful port of what DE-LIMP logs
in `R/server_data.R` / `R/helpers.R`. Two pipelines, picked by the bundle's
`de.method`.

## Invocation
```
Rscript scripts/run_de.R --input report.parquet --metadata conditions.csv \
        --method {dpc|maxlfq} --outdir de_results \
        [--contrasts "B-A,C-A"] [--q-cutoff 0.01] [--logfc 1.0] [--adjp 0.05] \
        [--block Mouse [--block-scope within|all]]
```
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
  by `lmFit(block =, correlation =)`. This is the nested case (mice within age, PROT_0756):
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

  **Why `within` is the default** (PROT_0756: 6 mice × 5 IPs, consensus correlation
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

## Contaminant filter (both paths, on by default)

`run_de.R` drops every precursor that maps to a `Cont_` entry before quantification —
any accession in `Protein.Ids` (`Protein.Group` for protein-level input), which is
DIA-NN's own `--cont-quant-exclude` rule: a peptide shared between a sample protein and
a contaminant entry can carry the contaminant's signal. The rule, the tag and the counting
live in `scripts/contaminants.R` (the tag mirrors `fetch_fasta.py`'s `CONT_TAG`; a test
asserts it).

Why it exists: DIA-NN's `--cont-quant-exclude` only shapes DIA-NN's OWN quantities, and
`run_de.R` re-quantifies from the report. On Silva08172026 (mouse brain IPs, 2026-09-24)
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
names the proteins. The states are `sidecar_state()` in `fetch_fasta.py`, mirrored in
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
- **QuantUMS quality filtering:** da Cruz Moschem J, Silva Campitelli de Barros BC,
  de Toledo Serrano SM, Chaves AFA (2025) *Decoding the Impact of Isolation Window
  Selection and QuantUMS Filtering in DIA-NN for DIA Quantification of Peptides and
  Proteins.* J Proteome Res 24:3860-3873. doi:10.1021/acs.jproteome.5c00009.
  **VERIFIED 2026-08-04** — an earlier note here said this "could not be verified";
  that was wrong. DE-LIMP cites it as "Moschem et al." (the first author's surname is
  da Cruz Moschem). QuantUMS computes three scores: protein-group MaxLFQ quality,
  empirical quality, and quantity quality, all measuring MS1/MS2 feature agreement.
  The skill filters on the first two (`--pgq-cutoff`, `--eq-cutoff`).
