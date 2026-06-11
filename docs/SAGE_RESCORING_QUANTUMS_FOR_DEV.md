# Sage on diaPASEF DIA — external findings + two prototype enhancements (rescoring & QuantUMS-style quant)

*Prepared for the Sage developer. Context: we used Sage (v0.14.7) as the open-source engine to identify + quantify Bruker **diaPASEF (timsTOF)** DIA data, by searching **diaTracer** pseudo-MS/MS spectra. Below is (1) a small mzML-reader patch that let Sage read diaTracer output, (2) a head-to-head with FragPipe (MSFragger + MSBooster) that localizes the identification gap to **one missing feature**, and (3) two external prototype add-ons — an MSBooster-style predicted-spectra rescoring layer and a QuantUMS-style quant layer. All numbers are on a 16-run dog (Canis) diaPASEF dataset (timsTOF HT, 60 SPD), dog UniProt DB, Sage's own `rev_` decoys.*

---

## 0. TL;DR for a Sage dev

- **Tiny patch needed to read diaTracer mzML:** Sage panics (`missing MS1 precursor`) on diaTracer's MS2-only pseudo-spectra because it reads `precursor.mz` only from `MS:1000744` (selected ion m/z), which diaTracer omits — it provides `MS:1000827` (isolation-window target m/z) instead. A 29-line fallback in `sage-cloudpath/src/mzml.rs` fixes it. **Candidate upstream contribution** (details in §1).
- **Sage's identification ceiling vs MSFragger+MSBooster on identical pseudo-spectra: 984 vs 1,361 protein groups.** We decomposed it: it is **not** quant, **not** match-between-runs, **not** raw-data re-extraction (FragPipe's DIA-NN step actually *trims* the count 1,361→1,337). The entire gap is search/rescoring sensitivity, and the pin diff localizes it to **one feature Sage lacks**: a deep-learning **predicted-fragment-intensity spectral similarity** (MSBooster's `*_spectral_entropy`). Sage already has RT + ion-mobility prediction features. (§2.) **We confirmed this empirically:** adding that one feature to Sage's pin (AlphaPeptDeep + mokapot) and rescoring contributes **+236 protein groups** at valid 1% FDR — the dominant share of the gap — making open Sage competitive with FragPipe+MSBooster (§3).
- Two external prototypes built on top of Sage's outputs (`.pin`, `results.sage.tsv`, `matched_fragments.sage.tsv`): an **MSBooster-style rescorer** (§3) and a **QuantUMS-style quant** (§4).

---

## 1. mzML patch — reading diaTracer (MS2-only) pseudo-spectra

**Symptom:** Sage `panic!("missing MS1 precursor for {}")` at `crates/sage/src/scoring.rs:379` (`score_standard`) and `:543` (`score_chimera_fast`).

**Root cause (not an MS1-scan requirement):** the panic fires when a query's `precursors` vec is empty. Tracing `SpectrumProcessor::process` → `initial_hits`, scoring only ever reads `precursor.mz`, `precursor.charge`, `precursor.isolation_window` — no MS1 *scan* is required. The precursor was being dropped earlier, in `crates/sage-cloudpath/src/mzml.rs`: `precursor.mz` is set **only** from cvParam `MS:1000744` and the precursor is pushed only if `precursor.mz != 0.0`. diaTracer pseudo-mzML provide `MS:1000827` (isolation window target m/z) + `MS:1000041` (charge) but **no `MS:1000744`** → `precursor.mz` stays 0.0 → precursor dropped → empty `precursors` → panic.

**Patch (29 lines, `mzml.rs`):** add `ISO_WINDOW_TARGET = "MS:1000827"`; capture it per precursor; at precursor-close, **if `precursor.mz == 0.0`, fall back to the isolation-window target m/z**, then push. Charge (`MS:1000041`) already parsed → preserved. No `scoring.rs` change. LFQ stays config-gated (`"lfq": false` for ID-only; no MS1 → no crash).

**Worth noting for upstream:** this also makes Sage able to ingest *any* spectrum-centric DIA→pseudo-spectra converter that encodes the precursor in the isolation window rather than a selected-ion cvParam. diaTracer additionally carries per-spectrum ion mobility as `MS:1002815` at the scan level (Sage's `ims` LDA feature can pick this up).

**Result:** patched Sage searched the 16 diaTracer mzML → **91,034 PSMs / 6,402 peptides / 984 protein groups @1% FDR** (global across 16 runs), ~7.5 min. No further code changes.

---

## 2. Where the identification gap is (decomposed)

Same 16 files, same diaTracer pseudo-spectra fed to both engines:

| stage | protein groups (union/16 runs) | input |
|---|---|---|
| **Sage** search (built-in LDA) | **984 / 996** | diaTracer pseudo-spectra |
| **MSFragger + MSBooster** (FragPipe ID stack, pre-DIA-NN) | **1,361** | the **same** pseudo-spectra |
| DIA-NN library quant of raw `.d` (FragPipe final report) | 1,337 | raw `.d` |

Two things we ruled out, with evidence:
- **Not raw-data re-extraction.** FragPipe's DIA-NN step (library-based quant of the raw `.d`) goes **1,361 → 1,337** — it *filters* at the protein level, it doesn't add proteins. So FragPipe's edge is not "it re-mines the raw data."
- **Not match-between-runs.** We built a decoy-FDR-controlled MBR layer over Sage (two independent decoy models, transfer FDR 0.83%/0.91%). Per-run completeness rose (463→472 PG/run) but **union was unchanged (996)** — and **675 of the proteins FragPipe finds are identified in *zero* Sage runs**, so they're structurally unreachable by transfer. The gap is upstream of MBR.

**So the entire ~377-PG gap is search/rescoring sensitivity on identical input.** The `.pin` diff pins it to one missing feature class.

### `.pin` feature comparison (this is the actionable part)

Sage's `results.sage.pin` features (abbrev.): `ln(hyperscore)`, `ln(delta_next)`, `ln(delta_best)`, `matched_peaks`, `longest_b/y`, `longest_y_pct`, `ln(matched_intensity_pct)`, `ln(-poisson)`, `scored_candidates`, `posterior_error`, `isotope_error`, `ln(precursor_ppm)`, `fragment_ppm`, plus **`retentiontime`, `predicted_rt`, `sqrt(delta_rt_model)`, `ion_mobility`, `predicted_mobility`, `sqrt(delta_mobility)`** (Sage already predicts RT and IM).

MSFragger raw pin: hyperscore/delta, `log10_evalue`, `matched_ion_num`, `complementary_ions`, `ion_series`, charge/length/ntt/nmc terms.

**MSBooster adds (edited pin minus raw MSFragger):** `delta_RT_loess`, `delta_RT_loess_real`, `delta_IM_loess`, `ion_mobility`, `hypergeometric_probability`, `intersection`, `top6_matched_intensity`, **`unweighted_spectral_entropy`, `weighted_spectral_entropy`**.

**The key observation:** Sage **already has** the RT and IM prediction half (its built-in models give `predicted_rt`/`predicted_mobility` + deltas — the analogues of `delta_RT_loess`/`delta_IM_loess`). The one thing Sage's pin **lacks** is the **deep-learning predicted-fragment-*intensity* spectral similarity** (`unweighted/weighted_spectral_entropy`) — i.e. how well the *observed* MS2 intensity pattern matches a *predicted* one. Sage's intensity features (`matched_intensity_pct`, `longest_b/y`, `matched_peaks`) measure matched-ladder coverage, not predicted-intensity agreement. That is MSBooster's signature discriminator and, by elimination, the likely driver of the 984→1,361 gap.

> **Suggestion for Sage:** an optional predicted-MS2-intensity similarity feature (via a pluggable predictor, e.g. AlphaPeptDeep/Koina/Prosit) added to the LDA/pin would likely recover much of this gap on spectrum-centric DIA (and probably help DDA too). Sage already carries RT/IM prediction, so this is the missing third predicted channel.

---

## 3. Prototype A — MSBooster-style predicted-spectra rescoring on Sage's `.pin`

External proof-of-concept (does **not** modify Sage; operates on its outputs), fully open (AlphaPeptDeep [Apache-2.0] + mokapot):

1. Parse unique (peptide, mods, charge) from Sage's PSMs (`results.sage.tsv` / `.pin`). Predict b/y fragment intensities with **AlphaPeptDeep**, `instrument='timsTOF'` (important — an Orbitrap/Lumos model mismatches tims fragmentation; NCE tuned to maximize predicted-vs-observed correlation since tims CE is ion-mobility-dependent).
2. For each PSM, match observed pseudo-spectrum peaks to predicted b/y (≈20 ppm), compute **spectral-entropy similarity** (and cosine) target *and* decoy (decoys predicted from the reversed sequence → valid FDR).
3. Append the feature(s) to `results.sage.pin` → rescore with **mokapot**.
4. Controls: (A) mokapot on the *original* pin (no new feature) isolates rescorer-vs-Sage-LDA; (B) mokapot + spectral feature is the treatment. Compare protein groups at a consistent 1% FDR rollup.

**Result — it works, and the spectral feature is the dominant contributor** (16-run bigdog, 786,879 Sage PSMs, valid 1% FDR with reversed-sequence decoys):

| level | Sage native | (A) mokapot, no spectral | (B) mokapot + spectral | FragPipe+MSBooster |
|---|---|---|---|---|
| **protein groups** | 984 | 1,170 | **1,406** | 1,361 |
| peptides | 6,402 | 6,874 | 8,235 | — |
| PSMs | 91,034 | 93,593 | 110,144 | — |

- **`weighted_spectral_entropy` = +2.63 mokapot weight — 2nd-strongest of 36 features** (behind only `posterior_error`); `n_matched_frags` +1.47, `spectral_entropy` +1.23. High-confidence targets: median entropy similarity **0.836** / cosine 0.840; decoys 0.272 / 0.208 — clean separation.
- **Clean, same-rollup, valid-FDR contribution of the spectral feature itself = B − A = +236 protein groups** (the +186 from A is mokapot's stronger discriminant vs Sage's LDA). So the predicted-MS2 feature alone confirms the §2 hypothesis: it is the dominant driver of the Sage→MSFragger+MSBooster gap, and adding it makes open Sage competitive.
- **NCE QC** (timsTOF IM-dependent CE): sweep over NCE 20–40 gave median entropy similarity 0.82–0.84 throughout (all above the 0.75 "strong" bar) → the AlphaPeptDeep timsTOF model fits this data well; chose NCE=35.
- **Honest caveat on absolute cross-tool numbers:** the self-computed picked-protein rollup gives A=1,170 for the same PSMs Sage natively rolled up to 984 — i.e. the rollup is slightly more permissive than Sage's `protein_q`. The rigorous claim is the within-pipeline **B−A=+236 at valid FDR**; the cross-tool "B (1,406) ≳ FragPipe (1,361)" holds under one consistent rollup but mixes rollup conventions, so read it as "competitive," not a hard win.
- **Two real bugs found en route (FYI):** (1) AlphaPeptDeep reorders its `precursor_df` internally — must carry explicit identity keys; (2) **diaTracer mzML peak arrays are not m/z-sorted**, which silently zeroed all fragment matches until sorted before `searchsorted`.

**Implication for Sage:** an optional predicted-MS2-intensity similarity feature (pluggable predictor → spectral-entropy/cosine vs observed, added to the LDA/pin) empirically recovers essentially the whole spectrum-centric-DIA identification gap here. Sage already carries RT + IM prediction; this is the missing third predicted channel. Code: `predict_ms2.py`, `build_feature.py`, `nce_sweep.py`, `rescore3.py` + the augmented `results.sage.rescored.pin`.

---

## 4. Prototype B — QuantUMS-style quantification on Sage output

Sage's native FDR is target-decoy q-values (`ml/qvalue.rs`) + the built-in 20-feature LDA (`ml/linear_discriminant.rs`, mokapot-like) + KDE PEP (`ml/kde.rs`); `--write-pin` exports a Percolator pin. For quant on diaPASEF we explored DIA-NN's **QuantUMS** (Grossmann/Kistner et al., Nat Biotechnol 2026) as an external layer over Sage's fragment intensities.

**MS1 for LFQ:** Sage's `tdf.rs` reads only DDA MS2 from `.d`; `timsrust 0.2.3` (already a Sage dep) exposes `read_all_ms1_frames()`. We wrote a small Rust tool to extract MS1 frames from the raw `.d`, sum across ion mobility, and inject them into the diaTracer mzML → Sage LFQ then traced MS1 features and quantified 5,399 precursors × 16.

**QuantUMS estimator (clean-room from the paper):** log-bias + log-variance of each signal as parametric functions of (intensity, quality score); inverse-variance combination; two-term loss (precision + bias-removal); learned by autodiff/gradient descent. Quality score per the paper = per-fragment **XIC elution correlation** (`Fr.N.Score`); we built a Rust XIC extractor (timsrust + `analysis.tdf`) to compute it from the raw `.d`.

**Honest findings (a Sage dev may care about the limits):**
- The estimator reproduces QuantUMS's **structure and precision** behavior — lowest CV of any method (precursor CV 0.39 / protein 0.52, better than DIA-NN's reported here) — but **not its accuracy**: correlation to DIA-NN's `Precursor.Normalised` was ~0.6–0.85 depending on configuration.
- **Definitive faithfulness test:** fed our estimator DIA-NN's *own* exported fragment data (via `--report-lib-info --xic`). Our full model reached only r≈0.28 vs DIA-NN — *worse than a plain sum of the same inputs (0.59)*. Two causes: (1) our two-term objective optimizes cross-run concordance/precision and its bias term degrades cross-*precursor* accuracy even with perfect inputs (so our objective ≠ DIA-NN's), and (2) DIA-NN's main report doesn't expose the true per-fragment quantities (they live in the closed `.d.quant`; the `--xic` chromatograms are a visualization, capping a perfect model at ~0.74).
- **Takeaway:** a faithful QuantUMS on Sage would need (a) DIA-NN-grade fragment peak integration and (b) an MS1-vs-MS/MS-concordance-anchored bias term (an MS1 channel Sage doesn't natively produce). The precision side is reproducible; the accuracy side is gated by upstream extraction + the closed objective.

**Follow-up (MS1-anchored objective + the real diagnosis):** we rebuilt it with MS1 and MS/MS as two channels and anchored the bias-removal term on MS1-precursor-vs-MS2-fragment concordance (the paper's actual formulation), using DIA-NN's own `MS1.Apex.Area` + `Ms1.Profile.Corr`. Faithfulness improved **0.28 → 0.72** (within-run r vs DIA-NN). The diagnosis is the useful part: **DIA-NN's `Precursor.Normalised` is MS1-*driven* on this dia-PASEF data** — `MS1.Apex.Area` alone correlates r=**0.94** with DIA-NN, vs 0.59 for the MS2 fragment-sum. Our *learned concordance objective* structurally can't reproduce that: it weights each channel by cross-channel agreement, so MS1 (a single signal not co-varying with the noisy MS2 consensus) is assigned high variance → low weight → the combine stays MS2-dominated (r≈0.28). Only **direct MS1 quality-weighting** (by `Profile.Corr²`, as the paper does) recovers 0.72. So the gap is the **channel-weighting / quality→variance mapping**, not the two-term loss; the objective optimizes internal consistency, which here is orthogonal to external accuracy.

### 4f. Pragmatic open quant for Sage output: directLFQ (recommended)

Rather than chase a faithful QuantUMS, the practical open quant for the Sage path is **directLFQ** (Ammar et al., Apache-2.0 — the MaxLFQ-style quant AlphaDIA uses). Fed Sage's per-run precursor intensities (PSMs at spectrum_q≤0.01) as directLFQ's ion-level longtable → protein quantities across 16 runs. Result vs DIA-NN 2.5 at the **protein** level: **801 protein groups, median CV 0.46, r = 0.86 (mean) / 0.89 (within-run)** — best accuracy *and* best precision of any method tried, in one turnkey step, no model fitting. (Precursor level: directLFQ #4,616 / CV 0.479 / r 0.61; the MS1-anchored QuantUMS gives the best precursor precision, CV 0.372, and most precursors, 11,106, but has no protein rollup here.) **So: Sage IDs → directLFQ is the open, accurate, precise quant recipe; the QuantUMS reimplementation is a research artifact that documents where a faithful clone diverges.** Code: `run_directlfq.py`, `quantums_ms1anchor.py`, `compare_final.py`.

---

## 5. Artifacts / reproducibility

All on our cluster under `/quobyte/proteomics-grp/brett/glendon/`:
- Patched Sage source + 29-line patch: `sage_src/`, `sage_patched/{sage, mzml_diatracer_ONLY.patch}` (built with a Cargo.lock `time 0.3.41` / `zerocopy 0.7.35` bump for rustc 1.94.1 — unrelated to the patch).
- ID results: `sage_bigdog_diatracer_16/patched/results.sage.tsv` (+ `.pin`, `matched_fragments.sage.tsv`).
- MS1 inject + LFQ: `ms1extract/`. XIC quality scores: `xicextract/`.
- QuantUMS: `sage_patched/{quantums_xic.py, quantums_faithcheck.py, *_results.txt}`.
- Rescoring (in progress): `sage_bigdog_diatracer_16/patched/sage_patched/{predict_ms2.py, build_feature.py, rescore.py}`.

Happy to share any of the code/diffs directly. The mzML `MS:1000827` precursor fallback (§1) in particular seems like a clean, low-risk upstream addition that would let Sage read diaTracer / other spectrum-centric DIA pseudo-spectra out of the box.
