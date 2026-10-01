# Protein-set tests (`run_sets.R`)

A comparison with few significant proteins can still hold a real, coordinated shift: twenty
proteins of one complex each up by 20% may pass no single-protein test and still be the
clearest biology in the run. Set tests look for that. They also make it easy to fool
yourself, so the skill runs them one way only.

## What it runs

```
Rscript scripts/run_sets.R --de-dir <run_de outdir> [--sets go|reactome|go,reactome|none]
    [--gmt sets.gmt --gmt-origin "<where it came from, and when>" [--gmt-kind public|own]]
    [--organism mouse] [--ontologies BP,CC,MF] [--min-size 10] [--max-size 500]
    [--contrasts A-B,...] [--camera-cor 0.01] [--fdr 0.05]
    [--ip-map ip_map.csv [--ip-normaliser interactome|bait] [--ip-min-enrichment 2]]
    [--outdir <de-dir>] [--figdir <de-dir>/../figures]
```

`run_de.R` suggests it for every comparison with fewer than 10 significant proteins
(`--sets-trigger`; `de_provenance.json` `set_tests.suggested_for`), and it can always be run.

**Per comparison, two tests side by side:**

| | question | reference |
|---|---|---|
| **camera** | are the set's proteins **more changed than the other proteins**? (competitive) | Wu & Smyth 2012, *NAR* 40:e133 |
| **fry** | are the set's proteins **changed at all**? (self-contained; roast with infinite rotations) | Giner & Smyth 2016, *F1000Res* 5:2605 |

P-values are adjusted across the sets of each comparison (Benjamini–Hochberg), separately
for each test. A set's `Call` joins the two BH lists; nothing is corrected across comparisons.
Then, for every set: the same tests with log2 run depth in the model, the fraction of its
proteins that are DPC presence calls, both tests on measured proteins only, and — for each
fry call — whether it holds on other reference bases. Each comparison gets a **reading**, one
generated sentence the report and the brief quote (below).

## One model — the DE model

`run_de.R` writes `set_test_inputs.rds`: each contrast's moderated t as reported (from the fit
that reported it), the expression matrix, the precision weights of each fit (DPC: limpa's
`voomaLmFitWithImputation` weights), the design, the contrast matrix, the block and its
consensus correlation. `run_sets.R` tests on exactly these. It refits only for the
sensitivity analyses (depth, bait normalisation), with `fit_like_run_de()` — the call
`run_de.R` makes — and first checks that this refit reproduces run_de's t for every tested
contrast (PROT_0756: max |Δt| 3e-12); it stops if it does not.

- **camera** ranks the contrast's moderated t (as z-scores, `zscoreT`, Hill's approximation,
  df = the fit's `df.total`) and takes its inter-protein correlations from the same model's
  standardised residual effects. With no block and no weights this **is** `limma::camera()` —
  the tests check the p-values agree to 1e-12 — and with a preset correlation it is
  `limma::cameraPR()` on the z-scores.
- **fry** gets the expression, weights, design and contrast, and for a random block the block
  and consensus correlation — `fry()`'s own computation (`.lmEffects`, `standardize =
  "posterior.sd"`, `.fryEffects`) on the canonical basis below.

### One canonical residual basis — and what it does not fix

camera's correlation estimate and fry's set statistic average a set's residual effects
coordinate by coordinate, and fry's robust variance takes each protein's largest squared
coordinate. limma takes the residual basis from a Householder QR — of the design, or with
precision weights of each protein's own weighted design — and that basis is **not unique**: on
PROT_0756, R's `qr()` of the same 30 × 10 design matrix came out with different reflection signs
on an AMD EPYC 9554 HIVE node than on EPYC 7532/7763 nodes (floating-point zeros in a 0/1
design landing either side of zero). Residual effects came out rotated, and camera and fry
p-values moved by up to 0.09 between otherwise identical runs. `lm_effects()` fixes the basis:

- **reference**: the polar factor of the residual projector (which depends only on the design's
  column space) applied to a fixed generic matrix (seed 20260928, R's portable RNG, global RNG
  state restored) — the same on every machine and for every way of writing the design;
- **a protein with weights**: the polar factor of the reference projected into its own weighted
  residual space — the orthonormal basis closest to the reference; constant weights give the
  reference.

Residual sums of squares and the contrast effect (limma's `.lmEffects`) are unchanged. What the
basis does and does not change (the statistician's review, rel-stats check1):

- **A change of reference is a common rotation** of every protein's basis — with weights too,
  since the polar factor carries a common rotation through. camera is invariant to it: it does
  not depend on the reference (to 1e-15 with weights, rel-stats check1b).
- **fry depends on the reference** only through its robust variance's largest-coordinate term:
  up to 0.03 in p without weights, 0.07 with (check1b; 1 of 200 calls changed).
- **Apart from the reference**, both differ from limma's own functions where limma's basis
  differs. Without weights camera is `limma::camera()` exactly, and fry is `limma::fry()` to
  7e-14 once that one term is taken from limma's basis (unit tests). With DPC-style weights
  limma's per-protein bases are not aligned, and camera differs from weighted `limma::camera()` by
  up to 0.04 in p, fry from weighted `limma::fry()` by up to 0.1 — a difference from limma, not a
  dependence on the reference.

The fixed reference makes the result reproducible, not right: another reference is equally
valid. So every fry-significant set is re-tested on four other references (seeds 1–4);
`fry_Reference_Stable` says `stable` or on how many it held, and a call that did not hold on
all is flagged "fry call depends on the residual basis" — read it as borderline. (On PROT_0756,
two other seeds changed 267 and 332 of 10,542 fry calls.)

### A random block (`--block`, e.g. several IPs per mouse)

`limma::camera()` has no block argument. Gordon Smyth (support.bioconductor.org/p/78299):
"you can't use duplicateCorrelation() with camera()"; for a blocked design the route is
`cameraPR()` on the moderated t from `lmFit(..., block, correlation)` + `eBayes`, "as long as
you're happy with a preset value for the inter-gene correlation", or `roast`/`fry`, which take
the block directly. `run_sets.R` does both: fry takes the block; camera takes the blocked fit's
statistics, and — so that its correlation need not be preset — estimates it from the residual
effects whitened by the block correlation (`limma:::.lmEffects`, the whitening fry itself
uses). A shared subject effect then stops inflating the correlation: on the synthetic pulldown
whitening halves it (0.85 → 0.43). Contrasts `run_de.R` reports from the fit
with samples independent (`--block-scope within`, between-mouse) use that fit.

### DPC and fry

limpa reduces the residual df of a protein whose whole group was imputed
(`voomaLmFitWithImputation`); that reaches camera through run_de's t, but fry works from the
values and weights, so for those proteins it sees the model's values as data. That is what
the presence-call columns are for: `Presence_Call_Fraction` (from the DE table's `Evidence`),
`Mostly_Presence_Calls` (> 0.5, default), and `measured_*` — both tests with presence-call
proteins removed from the sets and from the background.

## Which sets

A-priori sets only.

- **GO** BP, CC and MF from the organism's Bioconductor OrgDb, propagated to ancestor terms
  (`GO2ALLEGS`, as `limma::goana`), names from GO.db. Offline: `setup.sh` installs
  AnnotationDbi, GO.db, org.Hs.eg.db and org.Mm.eg.db (bioconda, else the Bioconductor source
  packages). The record names the GO release (`GOSOURCEDATE`) and package versions.
- **Reactome** (`--sets reactome`) when `reactome.db` is installed (about 2 GB; not by setup.sh).
- **A GMT file** (`--gmt`), with `--gmt-origin` saying what it is and `--gmt-kind`:
  - `public` — a published collection named with its version (`--gmt-origin "MSigDB Hallmark
    v2024.1"`): recorded as such;
  - `own` (the default, tagged) — a list someone drew up: valid only if drawn up before the
    results. Its sha256 and file time are recorded against the **first** DE run on this input
    (`de_provenance.json` `first_written`, which `run_de.R` carries forward on every re-run of
    the same input — a re-run must not make a late list look early); a file dated after it is
    flagged in the record, the Methods and on stderr. The residual: the date is carried forward
    only within one DE folder, so a run_de into a **new** outdir starts `first_written` again
    (the record's `de_results_first_source` says so). With `public` the date check is skipped
    and the record says that too — a published collection's file date says nothing about when
    its sets were chosen.
- **Never a list drawn up from the results.** In the Dog CSF session (2026-09-21) hand-built
  lists of the proteins that looked changed "confirmed" themselves: testing a list on the data
  that selected it finds what selected it. Those results were retracted.
- **Not gene-permutation GSEA.** `gseGO` / fgsea permute genes, treating a set's proteins as
  independent — the flaw that keeps `geneSetTest` out. A `--gsea` file given to
  `analysis_prompt.py` is presented as exploratory only, never beside these tests as a second one.

Organism: `--organism`, else the taxid in the search FASTA's `.meta.json`; neither → stop. An
organism with no OrgDb installed is mapped to **human** by gene-symbol identity (upper-cased):
a stand-in for 1:1 orthology that misses renamed genes and lineage-expanded families (the Dog
CSF session's DLA and mitochondrial families), recorded as such with its mapping rate.
Symbols map to Entrez IDs, then unambiguous aliases (an alias naming two genes names
neither). A protein group is in a set when any of its genes is. Size limits count the set's
proteins **among the tested ones**: 10–500 by default, and a set must leave at least 10
proteins outside it (a competitive test of the whole universe has nothing to compete with —
it matters for a pulldown's interactome).

## Why camera estimates the correlation (the default)

`camera(inter.gene.cor = 0.01)` is limma's default ("a useful compromise between strict error
rate control and interpretable gene set rankings"); `inter.gene.cor = NA` estimates it per set
("rigorous error rate control for all sample sizes and all gene sets"). Proteins in a set are
often correlated — complexes, co-regulation, and in pulldowns the IP's recovery, which moves
a whole interactome together. `Rscript scripts/sim_set_tests.R --what camera` (300
replicates per row, 1,000 proteins, one set of 20; latent 0.7 = the set shares a sample
factor, within-set correlation ≈ 0.33), rate of p < 0.05:

| runs per group | set | shift | camera, estimated | camera, preset 0.01 | fry |
|---|---|---|---|---|---|
| 3 | independent | 0 | 0.003 | 0.043 | 0.030 |
| 3 | correlated | 0 | 0.060 | **0.527** | 0.047 |
| 3 | independent | 0.4 | 0.203 | 0.520 | 0.373 |
| 3 | correlated | 0.4 | 0.107 | 0.527 | 0.070 |
| 6 | independent | 0 | 0.017 | 0.030 | 0.063 |
| 6 | correlated | 0 | 0.087 | **0.530** | 0.057 |
| 6 | independent | 0.4 | 0.757 | 0.830 | 0.803 |
| 6 | correlated | 0.4 | 0.240 | 0.690 | 0.147 |

With a preset 0.01, half of all correlated null sets come out significant. Estimating the
correlation keeps them near 5% (6–9%), at the price of power and of being conservative for
independent sets at 3 runs per group. `--camera-cor 0.01` restores the preset; the record
tags whichever is used.

### When camera has little power

camera takes a set's estimated correlation as within-set co-regulation and inflates the set's
variance by it. When **whole runs** move together — a pulldown's recovery, a sample's loading —
random protein sets correlate too, and camera's variance inflation is that global correlation,
not the set's own. `sim_set_tests.R --what camera`, part (b): every protein loads on one run
factor (random sets then correlate at 0.22–0.24); rate of p < 0.05:

| runs per group | set | camera, estimated | camera, preset 0.01 | fry |
|---|---|---|---|---|
| 3 | random | 0.000 | 0.037 | 0.047 |
| 3 | own factor, no shift | 0.037 | **0.457** | 0.070 |
| 3 | shifted 0.5 | 0.070 | 0.750 | 0.133 |
| 6 | random | 0.000 | 0.037 | 0.057 |
| 6 | own factor, no shift | 0.057 | **0.553** | 0.060 |
| 6 | shifted 0.5 | 0.180 | 0.940 | 0.223 |

The estimated correlation stays calibrated but finds a real shift 7–18% of the time. On
PROT_0756 random sets correlated at 0.13–0.145, which inflates a 54-protein set's variance about
8-fold: camera missed Kv2.1's own potassium-channel complex (FDR 0.25 / 0.63, where fry gave
9e-6 / 6e-8) and called nothing in any comparison. So `run_sets.R` measures it for every
comparison — `camera_global_correlation`, the mean estimated correlation of 300 random sets of
20, 50 and 100 proteins from a fixed seed, and `camera_vif_50 = 1 + 49 ×` that — and above 0.05
(`camera_low_power`) the reading says camera has little power and its count is not a result;
"significant in both" is then not a useful headline, and the report does not use it as one.

Not fixed here: removing the global run factor before estimating the set correlation. The
naive version — the correlation in excess of the random-set level — is miscalibrated (false
positives 0.10–0.15 on the review's simulation), and a fix goes in only when a simulation
validates it AND it recovers the Kv2.1 channel complex. For now the limitation is stated.

## Run depth

The Dog CSF lesson: a complement set gave p = 2 × 10⁻⁴, and 0.81 once log run depth was a
covariate — the groups differed in depth and the set tracked it. Here depth is
`Detection_Matrix.csv`'s column sums (DPC: precursors observed per run; MaxLFQ: proteins
quantified), log2, centred, added to the design with everything else unchanged
(`fit_like_run_de`). In a pulldown it gets **one slope per run type** (bait IPs, control IPs),
each centred within its runs: a control's depth and a bait IP's depth are not the same thing,
and one shared slope let the bait runs' slope stand in for the controls'. (PROT_0756, IgG Old vs
Young: with one shared slope 124 of 2,407 fry sets "held"; with the IgG runs' own slope, 0 —
the 124 were the largest movers of the same shift, not a signature.)

A set significant without depth and not with it (same test, same direction) is **"not
separable from run depth"**, and its flag gives the FDR and the set's effect with depth
(`depth_fry_FDR`, `depth_Mean_logFC`). That is not proof of an artefact — a depth difference can
be biology — but the set's change cannot be told apart from it. When most sets of a comparison
move together and are not separable from depth, the reading says exactly that and names no set.

Not run for **bait-vs-control** comparisons: there the depth difference is the enrichment (the
bait brings its interactors and their precursors with it), and adjusting for it would remove
what is being tested. The record gives each comparison's depth difference between the groups
(`depth_difference_log2`).

## Blocked designs: pooled between-subject comparisons

A comparison between different subjects (mice) that averages several samples per subject —
Old vs Young over all of a mouse's IPs — can only come from the blocked fit, and `run_de.R`
reports it with a CAUTION: one consensus correlation understates the between-mouse variance of
proteins with strong mouse-to-mouse variation. Sets inherit it (rel-stats check2: a set sharing
a mouse-level factor, no age effect: camera 0.23 and fry 0.16 false positives). Such comparisons
carry `between_block_caution` in the record, a flag on every significant set and a sentence in
the reading. Prefer between-mouse comparisons with one sample per mouse (per bait), which
`--block-scope within` reports from the independent fit.

## Pulldowns (`--ip-map`)

`ip_map.csv`: one row per group — `Group`, `Role` (`bait` | `control`), `Bait`, `Condition`,
optional `Bait_gene` (the bait protein's gene symbol).

**(a) Control vs control** (e.g. `Old_IgG-Young_IgG`, added to run_de's `--contrasts`): the
lysate background every IP carries. A shift here is in everything the controls pull down
nonspecifically, and it sits under every bait-vs-control and between-condition comparison.

**(b) Between conditions, relative to the bait's complex** (`Sets_baitnorm_<contrast>.csv`).
Bait recovery differs between conditions (PROT_0756: the Old IPs brought down 0.7–1.3 log2
units less of each complex), so an unnormalised Old-vs-Young IP comparison is mostly recovery.
For each bait:

1. **The interactome**: proteins enriched over the control in the bait's bait-vs-control
   comparison **pooled over conditions** — (mean of the bait's groups) − (mean of their
   controls), from run_de's model — at the FDR and at least 4-fold (`--ip-min-enrichment 2`,
   log2; tagged when it is the default).
2. Each of the bait's IPs is offset by the **median of the interactome** in that run.
3. Within the interactome, sets are tested between conditions with `fit_like_run_de` refitted
   on the bait's IPs (the between-condition contrast's own model: independent or blocked).
4. The **bait protein's own value** is used as a second reference when `Bait_gene` names it,
   and every significant set says whether it `holds with the other reference` or `depends on
   the reference`.
5. `Weak_Fraction`: the share of a set's proteins in the interactome's **bottom quarter** of
   pooled enrichment over the control; above one half the set is flagged "mostly weak
   interactors". (Until the review this was "within 1 log2 of the minimum", which flagged 73–88%
   of PROT_0756's tested sets and every simulated one, random nulls included — so with the
   defaults a finding could never stand unflagged. It is now relative to the interactome.)

**Why pooled, not "enriched at every condition".** An earlier version of this page claimed the
every-condition rule "cannot favour" either direction. That was wrong. When recovery differs,
"still enriched in Old" keeps a protein that fell in Old only if it did not fall far: PROT_0756's
JPH3 had 225 Old against 974 Young proteins enriched, so its interactome was "the proteins that
did not drop". The pooled comparison is uncorrelated with the between-condition comparison under
the null when the conditions have equal runs (independent filtering, Bourgon et al. 2010), so it
does not select on the answer.

**Why a minimum enrichment.** A weak interactor's IP signal is mostly the control background; if
the background differs between conditions, it reads as a change relative to the complex. The
minimum keeps them out — at a cost, below.

`Rscript scripts/sim_set_tests.R --what ip --reps 1000` simulates a pulldown the way it is
measured: an IP's intensity is the protein's nonspecific background plus, for interactors, its
specific signal (0.3–3 log2 over background), on the linear scale; noise 0.4 log2; Old recovering
the bait 1 log2 unit less; the Old runs' background 0 or +0.6 log2 higher. "kept" is how many of
a planted set's 20 members enter the interactome; "tested" how often that is run_sets.R's
10-member minimum; power and false-positive rates are over **all** replicates (an untested set is
not significant), so they are the whole pipeline's (HIVE job 24182369; fry / camera):

| Old background | set | enriched at every condition | pooled | pooled ≥ 2-fold | pooled ≥ 4-fold (default) |
|---|---|---|---|---|---|
| +0 | 20 fall 0.8 in Old: kept, tested, power | 10.4, 0.63, 0.57 / 0.56 | 19.4, 1.00, 0.96 / 0.96 | 16.3, 1.00, 0.95 / 0.95 | **4.5, 0.012, 0.011 / 0.011** |
| +0 | 20 rise 0.8: kept, tested, power | 18.6, 1.00, 0.99 / 1.00 | 20.0, 1.00, 0.99 / 0.99 | 19.4, 1.00, 1.00 / 1.00 | 9.8, 0.54, 0.51 / 0.53 |
| +0 | null, random | 0.062 / 0.000 | 0.034 / 0.001 | 0.054 / 0.003 | 0.073 / 0.002 |
| +0 | null, the 20 least enriched | **0.58 / 0.47** | 0.29 / 0.15 | 0.31 / 0.16 | 0.087 / 0.005 |
| +0.6 | 20 fall 0.8 in Old: kept, tested, power | **6.3, 0.13, 0.10 / 0.08** | 19.0, 1.00, 0.89 / 0.79 | 15.3, 0.99, 0.87 / 0.80 | **2.9, 0.002, 0.002 / 0.002** |
| +0.6 | 20 rise 0.8: kept, tested, power | 16.2, 1.00, 0.94 / 0.94 | 19.9, 1.00, 0.96 / 0.91 | 18.4, 1.00, 0.97 / 0.97 | 7.7, 0.23, 0.22 / 0.23 |
| +0.6 | null, random | 0.066 / 0.004 | 0.071 / 0.002 | 0.066 / 0.001 | 0.067 / 0.002 |
| +0.6 | null, the 20 least enriched | **0.86 / 0.86** | 0.67 / 0.49 | 0.68 / 0.59 | 0.087 / 0.010 |

The every-condition rule keeps half (or less) of a set that falls; pooling keeps it whole. Only
the 4-fold minimum stops weak interactors reading as a change (0.087 fry, against 0.29–0.86).

**What the minimum costs: losses are harder to detect than gains.** The minimum acts on the
enrichment pooled over conditions, so a protein that falls in one condition — binds the bait
less there, or leaves the complex — has a pooled enrichment lower by half the fall and can drop
below the minimum before any test. In this simulation, whose interactors sit at 0.3–3 log2, it
is severe: at the default a set of 20 that falls 0.8 log2 keeps 4.5 / 2.9 members and reaches
the 10-member minimum in 1.2% / 0.2% of replicates (power 0.011 / 0.002); one that rises 0.8
log2 keeps 9.8 / 7.7 and is tested 54% / 23% of the time. **How much it costs in a real pulldown
depends on how close that interactome sits to the minimum**, so `run_sets.R` measures it per
bait from the interactome itself (`ip_loss_sensitivity`: the pooled enrichment is the mean of
the per-condition bait-vs-control log2 fold changes, and a fall of d in one condition lowers it
by d/2) and each bait's reading says, e.g., "a 0.8 log2 fall in one condition would put about
37% of this interactome below the 4-fold minimum, and about 63% of the sets tested here would
keep enough members to be tested (about 42% for a 1.6 log2 fall)" — the second figure assumes
every protein of a set fell. **"About", and slightly optimistic:** each protein's fitted
enrichment is taken as exact and shifted by d/2, so noise near the 4-fold threshold — proteins
that would have fallen below it by chance, or failed the FDR half of the rule once weakened — is
ignored, and the real cost is somewhat higher. The report's pulldown section, the Methods and the brief say the same, and "no change
detectable" there is not evidence that the complex lost nothing. (At 2-fold the simulated cost
is small — 16.3 / 15.3 members of a falling set kept — but weak interactors then read as
changes, 0.31 / 0.68.) The fix is not a different threshold but a different correction for the
background: see "Known limitations" below.

**The weak flag, relative to the interactome.** Same simulation, pooled ≥ 4-fold: how often a
set is flagged "mostly weak interactors", by the flag run_sets.R uses now (over half of its
proteins in the interactome's bottom quarter of pooled enrichment) and by the one it replaced
(over half within 1 log2 of the minimum); Old background +0 / +0.6:

| set | flagged now | flagged by the old rule | fry significant and not flagged (now) |
|---|---|---|---|
| null, random | 0.003 / 0.002 | 1.000 / 1.000 | 0.073 / 0.067 |
| null, the 20 least enriched | 0.89 / 0.58 | 1.00 / 1.00 | 0.006 / 0.023 |
| 20 rise 0.8 | 0.010 / 0.014 | 0.99 / 1.00 | 0.51 / 0.22 |
| 20 fall 0.8 (when any are kept) | 0.26 / 0.25 | 1.00 / 1.00 | 0.009 / 0.001 |

The old rule flagged every set, so nothing could stand unflagged (a rising set: 0.003 / 0.000
significant and unflagged, rel-stats check6); the new one leaves random sets alone and still
catches the weakest.

**The false-positive rate is not 0.05.** Random null sets of the default interactome run at
0.070 (SE 0.005, 2,000 sets over both backgrounds; this simulation) and 0.0725 (SE 0.005; the
statistician's independent run, check6): mildly above 0.05 — slightly anti-conservative in this
model, not calibrated.

**Why the interactome median, with the bait as a cross-check** (pooled ≥ 4-fold interactome,
Old background +0.6, 1,000 replicates, untested sets not significant; fry / camera):

| scenario | set | interactome median | bait protein |
|---|---|---|---|
| planted only | 20 fall 0.8 | 0.001 / 0.001 | 0.000 / 0.001 |
| planted only | 20 rise 0.8 | 0.19 / 0.19 | 0.08 / 0.08 |
| planted only | null, random | 0.062 / 0.001 | 0.063 / 0.000 |
| broad: half the interactome +0.8 | null, random | **0.66 / 0.56** | 0.060 / 0.040 |
| bait epitope: the bait's read-out −0.8 | null, random | 0.058 / 0.000 | **0.38 / 0.001** |

- **The bait protein** is one noisy measurement (power 0.08 against 0.19 here), and anything
  that changes the bait's read-out but not the complex — an epitope masked with age, antibody
  affinity, a modification — shifts every set.
- **The interactome median** has the power and ignores the bait's own read-out, but it is a
  median: when a sizeable part of the complex really changes, the rest looks changed the other
  way.

So the median is the reference, the bait the cross-check, and a set whose call changes with the
reference is reported as unresolved. The table in the report gives, per bait, the recovery
difference each reference implies. When the interactome is small the empirical-Bayes variance
prior rests on few proteins (PROT_0756 JPH4: 35 at the old rule), and the reading says so below
100.

## Known limitations and backlog

- **Losses from a complex are harder to detect than gains** with the 4-fold minimum (above) —
  nearly undetectable in the simulation, less so in PROT_0756's interactomes. Each bait's reading
  states its own cost; not fixed in 2.9.
- **Random null sets run at about 0.07**, not 0.05, in the pulldown simulation (above).
- **camera has little power in pulldowns** (whole IPs move together); its zero there is not a
  result (see "When camera has little power").
- **A GMT's date check** is against `first_written`, which a run_de into a new folder starts
  again; with `--gmt-kind public` it is skipped (both recorded).
- **Backlog (after 2.9): background subtraction in place of the minimum (`pooled_bgsub`).**
  Subtract each condition's mean control (IgG) background on the linear scale, with no minimum
  enrichment, and floor what falls to or below the background (the statistician's check6 floored
  at 1/16 of the IP value). In that simulation, 1,000 replicates, Old background +0 / +0.6: a
  falling set keeps 19.4 / 19.1 members, fry power 0.98 / 0.96 (camera 0.97 / 0.92); a rising
  set 0.999 / 0.986; random nulls 0.065 / 0.085; the 20 least-enriched interactors 0.094 /
  0.094 (against 0.29 / 0.67 with no correction at all). **But that simulation generates the
  data with the same additive background + signal model the method assumes, so it cannot
  validate the method: validate it on real pulldowns first** (PROT_0756; a pulldown with known
  interactors or a spike-in), and choose the floor there. Tracked in `docs/TODO.md`.

## How to read the results

`sets_provenance.json` `reading` holds the rules, and each comparison's record holds a generated
`reading` sentence; the report, the analysis brief and AGENTS.md quote them from there, so the
wording cannot be strengthened on the way:

- camera is competitive and fry self-contained. A set's Call joins the two BH lists; no
  correction across comparisons. fry alone: the set moved, but no more than proteins in general
  — often part of a broad shift. camera alone: the set stands out, but the evidence it moved at
  all is weaker. (Dog CSF: a CNS set was camera-significant and fry p = 0.53.)
- Where camera has little power, its count is not a result and "significant in both" is not the
  headline.
- A set test says a predefined group of proteins shifted together — not that a pathway is
  activated or inhibited, and not which proteins drive it.
- GO terms nest and overlap: neighbouring significant terms are largely the same proteins.
- **"A broad shift" — most proteins moved together — is said only when the proteins did**: at
  least 60% of them in that direction and their median log2 fold change too. fry sets alone do
  not make one (PROT_0756 Old_JPH3 vs Young_JPH3: 45 fry sets lower in Old, but 53% of proteins
  lower and a median of 0.000 — the complex's recovery, not the proteome). A broad shift that
  cannot be separated from run depth is said to be exactly that. The sets that happen to survive
  depth adjustment are its largest movers, not a named signature.
- "fry call depends on the residual basis", "mostly presence calls", "mostly weak interactors",
  "depends on the reference", the between-block CAUTION: read those calls as unresolved. Calls
  that depend on the reference are given with their direction, their effect range and how far
  apart the two references put the recovery: an effect within that disagreement is within the
  reference uncertainty.
- **"Not detectable in this design"** — with n per side and the residual df — is how a comparison
  with nothing significant reads. It is not evidence that nothing changed.
- **Relative to a bait's complex, losses are harder to detect than gains** at the default 4-fold
  minimum (Pulldowns, above; each bait's reading gives its own numbers). Finding nothing there
  is not evidence that the complex lost nothing.

## Outputs

| file | what |
|---|---|
| `tables/Sets_<contrast>.csv` | per set: `Call`, camera and fry (direction, p, FDR), `fry_Reference_Stable`, the depth versions with `depth_Mean_logFC`, presence-call fraction, measured-only versions, `Flags`; columns described in `sets_provenance.json` `columns` |
| `tables/sets_members.csv` | every tested set's proteins (gene symbols) — the same in every comparison, so written once |
| `tables/Sets_baitnorm_<contrast>.csv` | pulldowns: the same, relative to the interactome median, with the bait-protein reference's calls, `Reference_Robust`, `Weak_Fraction` and `Members` |
| `tables/sets_provenance.json` | sources + versions, mapping rate, settings (each with its source; defaults tagged `DEFAULT -- not user-confirmed`), the model check, per-contrast counts, n per side, residual df, camera's global correlation and power, the between-block CAUTION and the generated `reading`; the pulldown record; the rules (`reading`), `columns`, `methods_paragraph` |
| `tables/sets_methods.txt` | the Methods paragraph + citations |
| `figures/sets_<contrast>.png` | each set's effect (mean log2FC of its proteins) against camera's FDR and fry's, side by side; colour = which test is significant; hollow = flagged |

`make_analysis_html.py` adds a "Protein-set tests" section from the record (a report section of
the same name replaces it); `make_methods.py` adds the paragraph to `methods.md`.

## PROT_0756 v2 (validation, 2026-09-28; after the statistician's review; HIVE job 24182898)

Read-only on the v2 session (its report, conditions and FASTA sidecar); `run_de.R` re-run to a
scratch folder with the 12 v2 contrasts plus `Old_IgG-Young_IgG` (random Mouse block, scope
within; HIVE r46 env: R 4.6.1, limpa 1.4.0, limma 3.68.5). The re-run has v2's 6,231 proteins
and every significance count; its logFC / t differ from v2's tables by at most 0.002 / 0.006
because v2 ran on a zen4 node and the re-run on zen2 — limpa's per-protein BFGS fit stops at a
slightly different point on each CPU family's BLAS kernels (reproduced exactly within a family).
Sets: GO BP/CC/MF (org.Mm.eg.db 3.23.0, GO release 2026-01-23) and Reactome (reactome.db
1.96.0): 5,410 GO and 808 Reactome sets of 10–500 tested proteins; 99.6% of protein groups
with a symbol mapped.

- **camera has little power in every comparison here**: random protein sets correlate at
  0.13–0.15 across runs (a 50-protein set's variance inflated about 8-fold). It called nothing,
  and that is not a result.
- **Bait vs IgG (what co-purifies)**: fry finds 51–3,593 sets per comparison (none for Old
  JPH4); which of them stand out from the rest of each IP cannot be told here.
- **(a) IgG Old vs Young — the lysate background**: fry 2,407 of 6,218 sets, every one higher
  in Old, as are 68% of the proteins (median +0.41 log2) — a broad shift; whether any category
  moved more than the rest cannot be told (camera has little power). It cannot be separated
  from the Old IgG runs' greater depth (0.64 log2): with the control runs' own depth slope, none
  of the 2,407 remains. No named signature.
- **Unnormalised Old vs Young per bait** (not corrected for bait recovery):
  - RyR: 89 of 92 sets higher in Old, as are 70% of the proteins (median +0.41) — a broad
    shift, though the Old RyR IPs brought down 1.31 log2 less of the complex.
  - JPH3: 45 sets, all lower in Old — but only 53% of proteins are lower and the median is
    0.000, so it is **not** a broad shift (review W1: an earlier reading called it one from the
    fry sets alone). 45 of 6,218 sets is about 0.7%; the Old JPH3 IPs brought down 0.94 log2
    less of the complex, which these tests do not correct.
  - JPH4 and Kv2.1: none.
- **(b) Relative to each bait's complex** (interactome chosen on the pooled enrichment, ≥ 4-fold):
  - JPH3 (212 proteins; 478 sets): fry 6 sets lower in Old relative to the interactome median
    (ER membrane terms; −0.12 to −0.09 log2), camera none. None holds with the bait-protein
    reference, whose recovery estimate differs by 0.19 log2 — the effects are within the
    reference uncertainty: unresolved, not a finding. (Under the old weak flag they were also
    "mostly weak interactors", like 73–88% of every bait's tested sets; by the
    interactome-relative flag none of the 6 is, and 0–1.1% of each bait's tested sets are.)
  - JPH4 (69 proteins; 180 sets), Kv2.1 (218; 601), RyR (121; 209): nothing — **not detectable in
    this design** (3 vs 3 mice, 4 residual df per bait); not evidence that the complexes are
    unchanged. JPH4's variance prior rests on 69 proteins.
  - **Losses are harder to see here than gains.** Each interactome's median pooled enrichment
    is 2.4–2.6 log2. A 0.8 log2 fall with age would put about 37% / 36% / 47% / 41% of the JPH3 / JPH4
    / Kv2.1 / RyR interactome below the 4-fold minimum, and about 63% / 63% / 60% / 54% of the tested
    sets would still keep 10 members if all their proteins fell that much (42% / 38% / 35% / 31%
    for a 1.6 log2 fall) — a real cost, milder than the simulation's. Each bait's reading gives
    its numbers; none of these results says a complex lost nothing with age.
  - The Old IPs brought down less of every complex: −0.94 / −1.13 (JPH3), −1.22 / −1.08 (JPH4),
    −0.74 / −0.70 (Kv2.1), −1.31 / −1.26 (RyR) log2 by the interactome / bait references.
  - The old "enriched at every condition" rule rested, for JPH3, on 225 Old against 974 Young
    enriched proteins — the selection bias the pooled rule removes.
