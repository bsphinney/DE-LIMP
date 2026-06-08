# De Novo Species-ID FDR: Validation & Method

**Status:** validated 2026-06-03 on the ocelot forensic case + a HeLa entrapment. This note records
*why* the de novo Target-Decoy pane uses **NovoBoard-style decoy spectra** (not shuffled peptides) and
*why* a spectral-decoy FDR alone is not sufficient for a species call.

## The question
De novo (Casanovo) peptides are BLASTed against NCBI **nr** and assigned a taxon by LCA. How do we put a
*calibrated false-discovery rate* on a species attribution when the organism may be absent from nr
(non-model wildlife forensics)?

## Method: decoy SPECTRA, target–decoy competition (NovoBoard)
Following **Tran et al. 2024, *Mol Cell Proteomics* (doi:10.1016/j.mcpro.2024.100849)**: for every real
MS/MS spectrum, build a **decoy spectrum** by removing a fraction `FRAC` of its peaks and adding back the
same number sampled from the global peak pool (precursor m/z + charge preserved). Run Casanovo on both,
BLAST both against nr, and estimate FDR by competition: `FDR(t) = (# decoy hits ≥ t) / (# real hits ≥ t)`.
This is the de-novo analogue of database target–decoy. Generator: `scripts/denovo_decoy_gen.py`.

This **replaces** the old shuffled-peptide decoy (residues permuted, BLASTed vs nr). Shuffled peptides
test "does a scrambled *sequence* hit by chance"; decoy spectra test the stronger null "does a peptide
*Casanovo would actually emit from realistic spectra that encode no real peptide* hit by chance."

## Finding 1 — the scramble fraction matters (`FRAC ≥ 0.8`)
HeLa decoy peptides hitting nr, by `FRAC`:

| FRAC | decoy peptides | nr hits | mammal hits |
|------|----------------|---------|-------------|
| 0.5  | 45,316         | 194     | **116** (GAPDH, β-actin, tubulin) |
| 0.8  | 45,803         | 3       | 0 |
| 0.9  | 45,898         | 5       | 0 |

At 50% scrambling Casanovo still reconstructs the most abundant real peptides — they **leak** into the
null with top bitscores (60–80, E≈1e-18) and poison the competition. **At ≥80% the null is clean.** The
pipeline uses **FRAC = 0.8**.

## Finding 2 — three different error rates (don't conflate them)
Entrapment on HeLa (Sage 1% FDR vs human proteome = ground truth; scan-level Casanovo-vs-Sage match):

| Error type | Rate | Controlled by |
|------------|------|---------------|
| **Chance / spurious homology** (decoy-spectra FDR) | ~0 once null is clean | decoy-spectra competition |
| **De novo exact-sequence error** | 45% overall; 21% at conf ≥ 0.95; 12% at conf ≥ 0.99 | **Casanovo confidence** |
| **Species-attribution error** | 36% non-mammal by *best-hit*; mostly conserved peptides | **LCA, not best-hit** |

De novo exact accuracy calibrates cleanly with confidence: ≥0.90 → 71%, ≥0.95 → 79%, ≥0.97 → 83%,
≥0.99 → 88% (scan-level, n=22,970 Sage-confident spectra). nr-BLAST recall is **length-limited**: among
*known-correct* peptides, recall is ~0% at ≤12 aa but 71% at 15–20 aa and 97% at 20+ aa — short exact
peptides cannot clear `E≤1` against nr's ~700M sequences. So nr species-ID is a **longer-peptide method**.

Crucially, a wrong *sequence* often still gives the right *species*: of exact-wrong nr-hits, **64% still
land on a mammal** (BLAST mismatch tolerance). That is why the ocelot converged on Felidae despite
imperfect sequences.

## Conclusion for the pipeline
A spectral-decoy FDR is **necessary but not sufficient**. A defensible de novo species call needs **three
controls together**:
1. **decoy-spectra FDR** (FRAC≥0.8) — kills chance homology;
2. **Casanovo confidence gate** — controls de novo sequence correctness (calibrated against the entrapment);
3. **LCA over all hits** (not best-hit) — controls conserved-peptide mis-assignment + the nr
   over-representation artifact.

The Target-Decoy FDR pane reports (1) and (2) and points to the Species (LCA) pane for (3).

## Composite score (research)
A confidence-led composite (de novo confidence × homology evidence, LCA-gated) ranks correct hits far
better than BLAST score alone (AUC vs Sage truth: mean-confidence 0.85, CWI 0.82, **bitscore 0.43** —
bitscore tracks length/detectability, *not* correctness). Bitscore-led composites are therefore avoided.

## Update 2026-06-04 — the composite/mokapot FDR "advantage" was decoy leakage, not real
A learned rescorer (real **mokapot 0.10.0**, plus a hand composite `−log₁₀E + 2·n_hi`) *appeared* to accept
~2× more peptides at 1% decoy-FDR than bitscore on the ocelot — **but that run used the FRAC-0.5 decoy**
(356 decoy nr hits, global FDR 3.6%, bitscore yield@1% = 2,833; mokapot ≈ 6,750). On the **clean FRAC-0.8
decoy** the null is near-empty: ocelot **36** decoy hits / 9,852 real = **0.37%** global FDR; HeLa **3** /
5,068 = 0.06%. With global FDR already < 1%, **every nr hit passes at 1% FDR and bitscore, composite, and
mokapot are identical** (yield = all hits). The "2×" was mokapot separating real peptides from *leaked*
decoys; with a clean decoy there is nothing to separate.

**Consequence:** do **not** calibrate the homology FDR with mokapot/composite — on a clean (FRAC-0.8) decoy it
adds nothing. The composite/confidence belongs on the **correctness** axis (Finding 2; mean-confidence AUC
0.85 vs Sage truth), calibrated on entrapment, not on the homology decoy. The Target-Decoy pane's bitscore
FDR is correct and, with FRAC-0.8 decoys, will (correctly) report ~all nr hits as real homology — chance
homology ≈ 0. mokapot was *fed all hit peptides unfiltered* (≥7 aa, has-a-hit) and still showed no edge on
the clean decoy. Analysis: `/tmp/decoyval/` (run_mokapot.py, hela_truth_mokapot.py, ocelot_decoy_f08_hits).

**Full f08 re-run (2026-06-05)** — the ocelot document (`docs/denovo_decoy_method.html`) was rebuilt on a clean
FRAC-0.8 decoy (the earlier draft used FRAC-0.5, which leaked). Clean-decoy results: **standard nr** real
9,852 / decoy **36** (global FDR 0.27% → all nr hits are real homology); **EXTRA nr-Euk** bitscore-calibrated
**10,404 @1% FDR** (bits≥33) vs standard's 9,852 = **+6%** (the leaky draft's "2×" was an artifact;
enrichment at bits≥36 rose from 96× to 323×); **SwissProt 2×2** decoy-spectra null **0.030%** ≈
decoy-database null 0.023% (**1.3×**, not the f05 "12×"). EXTRA's raw hit rate still saturates at ~99% under
`--min-score 1` regardless of FRAC (not FDR-controlled). EXTRA still outputs valid E-values (`--min-score 1`
only disables the E-value *filter*), and E-value calibration beats bitscore under EXTRA (13,074 vs 10,404
@1%, length-normalized); a confidence+E-value model beats both. Net: with a clean
decoy the homology FDR is ~trivially satisfied; the limits are de novo correctness + LCA.

## Artifacts
HPC: `…/Ocelot_search_2/decoy_validation/` and `…/HeLa_entrapment/` (Sage, FRAC sweep, lineage).
Analysis: `/tmp/decoyval/` (composite_fdr.R, scan_accuracy.R, true_fdr.R). Figures: composite_fdr,
nr_recall_by_length, frac_leakage.
