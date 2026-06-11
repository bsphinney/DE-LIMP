# Open-Source Sage diaPASEF Stack — Technical Handoff

**Purpose:** A fully license-free path to identify + quantify Bruker **diaPASEF (timsTOF)** DIA data using the open-source **Sage** search engine (MIT), via diaTracer pseudo-spectra, plus a **timsrust** MS1-injection LFQ step and a clean-room **QuantUMS** quant prototype. Built and validated on UC Davis HIVE in June 2026. This document is written so another Claude Code (or engineer) can reproduce, extend, or fix it without re-deriving anything.

> TL;DR results on 16 dog diaPASEF files: patched Sage → **91,034 PSMs / 6,402 peptides / 984 proteins @1% FDR**. With timsrust MS1 injection, Sage LFQ quantified **5,399 precursors**; a QuantUMS-style estimator over Sage fragments gave the **lowest CV (0.55 precursor / 0.71 protein)** and correlated r≈0.85 with DIA-NN 2.5 quant. **UPDATE (corrected):** that gap is a **search/identification** gap — NOT match-between-runs (a rigorous decoy-FDR-controlled MBR added 0 to the union; disproven) and NOT raw-`.d` requant (FragPipe's DIA-NN step *filters*, doesn't add). **Open predicted-spectra rescoring closes it:** adding an AlphaPeptDeep predicted-MS2 spectral-similarity feature to Sage's pin + mokapot adds **+236 PG at 1% FDR → ~1,400**, matching FragPipe+MSBooster (1,361). Independently, **AlphaDIA (open, empirical library) = 1,361 PG** on the same 16 files. So two fully-open paths now match the proprietary/free-academic tools.

### Current scoreboard (16 dog diaPASEF files, protein groups @1% FDR)

| Tool | PG | License | note |
|---|---|---|---|
| **Patched Sage + predicted-spectra rescoring** | **~1,400** (+236 from the open feature) | MIT + Apache | §4e |
| **AlphaDIA (empirical library)** | **1,361** | Apache | needs `extraction_backend: python`; library was FragPipe's EasyPQP |
| FragPipe + MSBooster | 1,361 | free-academic | |
| Patched Sage (native, no rescoring) | 984 | MIT | |
| DIA-NN 2.5 (library-free) | 1,306 | proprietary | sensitivity 2×2: autocal `--mass-acc 0` fails here, relaxed-prot-inf is neutral → fixed-15ppm/default already optimal |
| Spectronaut (allDog, 215-run cohort) | 3,048 | proprietary | cohort advantage |
| AlphaDIA transfer-learning (library-free directDIA) | running (very slow: ~8M-precursor predicted lib × python backend × IM modeling × 3 passes) | Apache | paper-acknowledged timsTOF speed limit |

---

## 0. Cluster environment & access

- **SSH:** `ssh -o ConnectTimeout=25 -i ~/.ssh/id_ed25519 brettsp@hive.hpc.ucdavis.edu "bash -l -c '<cmd>'"`
- **SLURM needs a login shell** (`bash -l -c`). **Never run heavy compute on the login node** — submit via `sbatch`/`srun`.
- **Partitions:**
  - `high` — account `genome-center-grp`, non-preemptible, ~64-CPU/user cap. Use for long monolithic jobs.
  - `low` — account `publicgrp`, qos `publicgrp-low-qos`, add `--requeue`. Preemptible, huge capacity. Use for parallel arrays + builds.
  - GPU — `--account=genome-center-grp -p gpu-a100 --qos=genome-center-grp-gpu-a100-qos --gres=gpu:a100:1`
- **Big-mem nodes:** 2 TB nodes exist on both `high` and `low` (request `--mem=1500G`).
- **Shell-quoting gotcha:** parentheses inside `bash -l -c '... echo (x) ...'` break the remote shell. Write scripts to files + `scp`, or avoid `()` in echo strings.

---

## 1. What the stack is (architecture)

```
Bruker .d (diaPASEF)
   │
   ├─(A) diaTracer ── pseudo-MS/MS mzML (MS2-only, precursor in isolation-window) ─┐
   │                                                                                │
   └─(B) timsrust (ms1extract) ── MS1 mzML (IM-flattened) ──► merge_ms1_ms2.py ─────┤
                                                                                    ▼
                                                              PATCHED SAGE (search + FDR + LFQ)
                                                                    │
                                                                    ├─ results.sage.tsv  (IDs, q-values)
                                                                    ├─ matched_fragments.sage.tsv (fragment intensities)
                                                                    ├─ lfq.tsv (MS1 LFQ, needs merged MS1)
                                                                    └─ .pin (Percolator/mokapot input)
                                                                    │
                                                       quantums_proto.py (QuantUMS-style quant)
```

- **(A) diaTracer** turns each diaPASEF `.d` into pseudo-MS/MS (MS2-only). diaTracer is *closed-source/licensed* — used here only to generate spectra. (A separate effort attempts a clean-room reimpl from Li et al., Nat Commun 16:95 (2025).)
- **Patched Sage** searches the pseudo-MS2 (this is the core open engine; the patch is what makes it read diaTracer files).
- **(B) ms1extract + merge** brings MS1 back from the raw `.d` (via `timsrust`, which Sage already depends on) so Sage can do MS1-based LFQ.
- **QuantUMS prototype** is an alternative quant over Sage's fragment intensities.

---

## 2. Component 1 — Patched Sage (identification)

**Source tree:** `/quobyte/proteomics-grp/brett/glendon/sage_src/` (git clone of `lazear/sage` tag **v0.14.7** — note v0.14.7 is the *latest* Sage release).

### Root cause of the original failure
Sage panicked on diaTracer mzML with `panic!("missing MS1 precursor for {}")` at `crates/sage/src/scoring.rs:379` (`score_standard`) and `:543` (`score_chimera_fast`). The panic fires when a query's `precursors` vec is **empty** — NOT because Sage needs an MS1 *scan*. The precursor was being dropped by the mzML parser:

- In `crates/sage-cloudpath/src/mzml.rs`, `precursor.mz` is set **only** from cvParam `MS:1000744` (selected ion m/z), and the precursor is only pushed if `precursor.mz != 0.0`.
- diaTracer pseudo-mzML provide `MS:1000827` (**isolation window target m/z**) + `MS:1000041` (charge) but **omit `MS:1000744`** → `precursor.mz` stays 0.0 → precursor dropped → empty `precursors` → panic.

### The patch (29 lines, `crates/sage-cloudpath/src/mzml.rs`)
- Add constant `ISO_WINDOW_TARGET = "MS:1000827"`.
- Capture it into a per-precursor `iso_window_target` in the cvParam handler.
- At precursor-close: **if `precursor.mz == 0.0`, fall back to `iso_window_target`**, then push. Reset the per-precursor iso-window vars.
- Charge (`MS:1000041`) is already parsed → preserved. **No `scoring.rs` change needed.**
- **LFQ is config-gated** (`main.rs:301`): set `"lfq": false` for ID-only (no MS1 → no crash); no code change.

**Patch artifacts:** `/quobyte/proteomics-grp/brett/glendon/sage_patched/`
- `mzml_diatracer_ONLY.patch` — the 29-line functional change
- `diatracer_ms2only.patch` — full diff incl. the Cargo.lock toolchain bump (below)
- `Cargo.lock.patched`
- `sage` — the built patched binary (md5 `dbeb0e99…`)

### Build notes (toolchain)
Building under the cluster's rustc 1.94.1 hit a **pre-existing** transitive failure (`time v0.3.34`), unrelated to the patch. Fix before `cargo build --release`:
```bash
cargo update -p time --precise 0.3.41
cargo update -p zerocopy --precise 0.7.35
```
Build on a compute node (low partition, ~16 cpu, ~30 min): `cargo build --release` → `target/release/sage`.

### FDR / rescoring (stock Sage — the patch does NOT touch this)
- **q-values:** target-decoy competition, `crates/sage/src/ml/qvalue.rs`. Decoys labeled `-1`; sort by score desc; FDR = decoys/targets; cumulative-min for monotonic q. Computed at **spectrum / peptide / protein** level (`spectrum_q`, `peptide_q`, `protein_q` columns).
- **Built-in semi-supervised LDA** (mokapot-like, native): `crates/sage/src/ml/linear_discriminant.rs` — 20 features (rank, charge, ln1p(hyperscore), delta_next/best, delta_mass_model, isotope_error, average_ppm, poisson, matched_intensity_pct, matched_peaks, longest_b/y, peptide_len, missed_cleavages, **rt, ims, delta_rt_model, delta_ims_model**). Linear discriminant trained on target/decoy.
- **PEP** via KDE: `crates/sage/src/ml/kde.rs::posterior_error`.
- **mokapot/Percolator interop:** `--write-pin` emits a `.pin` (`sage-cli/src/input.rs`).

### Run (ID-only)
Config = v0.14.7 schema (see §5). Example:
```bash
/quobyte/proteomics-grp/brett/glendon/sage_patched/sage \
  --write-pin --output_directory <out> <sage.json> <*_diatracer.mzML ...>
```
**Result (16 dog files):** 91,034 PSMs / 6,402 peptides / **984 proteins @1% FDR** (global across 16 runs), 7m38s.
Output: `/quobyte/proteomics-grp/brett/glendon/sage_bigdog_diatracer_16/patched/results.sage.tsv`

> On **Orbitrap/Astral DDA** (real MS1 present) Sage works *unpatched-relevant* and natively — the human-bone Astral run gave 100,814 PSMs / 623 peptides / ~150 proteins. Dir: `/quobyte/proteomics-grp/brett/glendon/sage_NP_astral_20260609/`. (The patch only matters for MS2-only diaTracer input.)

---

## 3. Component 2 — ms1extract (timsrust MS1 → LFQ)

**Dir:** `/quobyte/proteomics-grp/brett/glendon/ms1extract/` (`ms1extract_main.rs`, built binary `ms1extract`, `merge_ms1_ms2.py`, `merged/`).

- Sage's `crates/sage/src/tdf.rs` only reads **MS2** from `.d` (no MS1 path). But **`timsrust 0.2.3`** (already a Sage dep) exposes `read_all_ms1_frames()`, `get_tof_converter()`, `get_frame_converter()`.
- **`ms1extract`** reads MS1 frames from a `.d`, **sums intensity across the ion-mobility dimension** (IM-flatten), converts TOF→m/z, emits a centroided MS1 mzML. Thresholds: intensity floor 30 + top-2000 peaks/scan → ~26 MB/file (vs 1.9 GB unthresholded).
- **`merge_ms1_ms2.py`** injects the MS1 spectra into the diaTracer MS2 mzML by retention time (handles diaTracer's minute RT units), renumbering spectrum indices. Output: 16 merged mzML (1,146 MS1 + 175,514 MS2 each) in `merged/`.
- Run patched Sage with `"lfq": true` on the **merged** mzML → Sage logs "tracing MS1 features → discovered 5,184 MS1 peaks @5% FDR" → writes `lfq.tsv` = **5,399 quantified precursors × 16 runs**.
- **Caveat:** the injected MS1 + diaTracer MS2 carry no ion mobility in the precursor record, so Sage's LFQ is m/z+RT-based, **not IM-resolved**.

---

## 4. Component 3 — QuantUMS-style quant prototype

**Script:** `quantums_proto.py` (in `sage_patched/`). **Paper:** QuantUMS, Demichev lab, bioRxiv **2023.06.20.545604**.

- Core estimator (clean-room from the paper, NOT from DIA-NN code): each signal gets an estimated **log-variance** decreasing in intensity + a quality score; combine by **inverse-variance-weighted aggregation**: `q̂ = Σ wᵢ·qᵢ / Σ wᵢ`, `wᵢ = 1/σᵢ²` (minimum-variance estimator).
- Operates over Sage's `matched_fragments.sage.tsv` (fragment intensities). Quality proxies: relative fragment intensity + PSM posterior-error confidence (Sage lacks DIA-NN's XIC correlation feature).
- Output: 11,403 precursor × 16 + 2,398 protein × 16 matrices.
- **Scope (heuristic version):** uses the estimator with hand-picked uncertainty proxies — see §4b for the LEARNED upgrade.

### 4b. Learned QuantUMS (the ML fitting, no longer stubbed)

**Script:** `quantums_learned.py` (+ `qlearn_*_quant_learned.tsv`, `qlearn_learned_hyperparams.json`, `compare_quant4.py`, `quantums_learned_results.txt`, all in `sage_patched/`).

We replaced the hand-picked proxies with a **genuinely learned** model fit by gradient descent on the paper's concordance objective:
- `log σ²_i = v·f_i`, `log-bias_i = b·f_i` (per-signal log-variance + log-bias as linear functions of features `f_i`); consensus per (precursor, run) = inverse-variance-weighted mean of bias-corrected log-intensities.
- Loss = Gaussian NLL of each fragment vs the consensus (`Σ resid²/σ² + log σ²`). Fit with Adam, precursor-level 80/20 holdout + L2, plus a shot-noise prior (variance non-increasing in intensity) to partly substitute for the missing XIC quality feature.

**VERIFIED method (from the Nature Biotech 2026 Methods, PDF obtained):**
- Input per precursor: signal intensities of **MS1 precursor + each MS/MS fragment** across every acquisition it was identified in, plus **one quality score per signal**.
- The quality score in DIA-NN = **the correlation-based similarity score comparing a feature's extracted ion chromatogram (XIC) to that of the 'best' fragment ion** for each elution peak. ← this is the per-fragment XIC-correlation feature Sage lacks.
- QuantUMS models **log-bias** and **log-variance** of each signal, *each as a function of (intensity, quality score)* with the algorithm hyperparameters as coefficients.
- Precursor quantity = a formula summarizing the **bias-corrected** intensities while **minimizing the estimated log-variance** of the result (inverse-variance combination).
- Hyperparameters learned from data by **automatic differentiation + gradient descent**. **Loss = pairwise differences between quantities estimated from different features (MS1 vs MS/MS)**, with **two terms — one prioritizing precision, one prioritizing bias removal**; their balance = the high-precision vs high-accuracy modes (configurable). Hyperparameters fit **separately for MS1 and MS/MS** signals.
- **Standalone QuantUMS source (replicates DIA-NN 1.8.2 beta 25):** https://osf.io/yrvuc/?view_only=6a141dd5b8de405f95a6049ebb25913d · figure scripts: https://figshare.com/s/ab9b1cebc31db0e8a1bc

So the right path is NOT to reimplement from scratch — it's to (i) use the **standalone reference source** for the exact objective, and (ii) manufacture the one input Sage lacks: a **per-fragment XIC-correlation quality score** (extend `ms1extract` to extract fragment XICs from the raw `.d` and correlate them). Our `quantums_learned.py` (below) is the from-abstract reconstruction, superseded by this plan.

**Result vs DIA-NN 2.5 (= real QuantUMS):**

| method | prec CV | r vs DIANN | prot CV | r vs DIANN |
|---|---|---|---|---|
| Sage-LFQ | 0.864 | 0.851 | 1.036 | 0.832 |
| QuantUMS heuristic | 0.548 | 0.847 | 0.707 | 0.762 |
| **QuantUMS LEARNED** | **0.396** | 0.668 | **0.516** | 0.740 |

**Honest finding:** the learned model gives the **lowest CV of any method** (best precision — it does exactly what the concordance objective rewards) but **correlates LESS with DIA-NN** (r 0.847→0.668 precursor). This is robust (holds on >20-fragment precursors and across constrained/unconstrained fits), so it's the estimator's central tendency, not a fragment-count artifact: the concordance objective flattens fragment weights (trusts consistent noise-floor fragments) → more precise, but diverges from DIA-NN's top-fragment-driven abundance.

**The capping feature gap (the real frontier):** Sage is spectrum-centric and lacks DIA-NN's per-fragment XIC elution correlation. We first substituted coarser PSM-level proxies — but then built the real feature (see §4c).

### 4c. With the real XIC quality score (executed) — and the true limiting factor

**The OSF standalone is source-less** (private project, view-only token; the 84 MB zip is runtime DLLs + eigen + config, no Linux binary — QuantUMS is closed compiled C++ inside the DIA-NN Windows binary). But its config `diann-output-config-16-threads.txt` gave the **authoritative input contract**: per-fragment quality = `Fr.N.Score` (fragment XIC correlation), MS1 quality = `Ms1.Profile.Corr`, up to 12 fragments, `--high-acc` = bias-removal mode.

We **built the missing feature**: `xicextract` (Rust, timsrust 0.2.3 + direct `analysis.tdf` SQLite for the dia-PASEF window layout) — for each Sage precursor it picks the matching isolation window (m/z + mobility band), builds each matched fragment's XIC across RT, and correlates each fragment's XIC to the precursor consensus = the QuantUMS `Fr.N.Score`. Ran on all 16 `.d` (SLURM array); scores are discriminative (median ~0.90 with a real interfered-fragment low tail). Then `quantums_xic.py` = log-bias + log-variance as functions of (intensity, XIC score), bias-corrected inverse-variance combination, the **exact two-term loss** (precision + bias-removal, `--high-acc` λ balance), autodiff fit.

**6-way comparison (vs DIA-NN 2.5 = real QuantUMS):**

| method | prec CV | prec r vs DIANN | prot CV | prot r vs DIANN |
|---|---|---|---|---|
| Sage-LFQ | 0.864 | 0.851 | 1.036 | 0.832 |
| QuantUMS heuristic | 0.548 | **0.847** | 0.707 | 0.762 |
| QuantUMS learned (no XIC) | 0.396 | 0.668 | 0.516 | 0.740 |
| QuantUMS-XIC (balanced) | **0.394** | 0.606 | 0.520 | 0.680 |
| QuantUMS-XIC (high-acc) | 0.396 | 0.602 | 0.523 | 0.677 |

**The verified diagnosis — the model is faithful; the *quantity extraction* is the cap.** The fit confirmed the feature works exactly as intended: `v[xic_score] = −0.666` (higher XIC correlation → much lower variance; the no-XIC model had ~0 quality coefficients), and precision is excellent (CV ~0.39, better than DIA-NN). But DIA-NN correlation did **not** improve. Root cause, isolated: the **fidelity of our XIC-extracted fragment *areas* (quantity), not the weighting model or the quality feature.** Summing our XIC fragment areas gives r=0.732 vs DIANN; Sage's diaTracer pseudo-spectrum intensities give r=0.847; the two disagree (r=0.773). Our simple XIC engine (fixed ±15 s window, mobility-band sum, 20 ppm, no peak-shape/alignment) yields noisier *abundances* than DIA-NN's tuned integration — capping every XIC variant at r≈0.6–0.73 regardless of weighting. The XIC *correlation* (quality) is high-fidelity; the XIC *area* (quantity) is the weak link.

**Interim conclusion (later REFUTED — see §4d):** we *thought* the algorithm was faithful and only extraction capped accuracy.

### 4d. Definitive faithfulness test (executed) — our QuantUMS is precision-faithful, NOT accuracy-faithful

Fed our estimator **DIA-NN's own fragment data** (re-ran DIA-NN 2.5.1 on the 16 `.d` with `--report-lib-info --xic --xic-theoretical-fr`; reconstructed per-fragment quantity = XIC area, score = consensus XIC correlation, from DIA-NN's own chromatograms; reference = DIA-NN `Precursor.Normalised`). `quantums_faithcheck.py`, results in `sage_patched/` + `diann251_fragexport16/`.

| estimator on DIA-NN's own fragments | r vs DIA-NN Precursor.Normalised |
|---|---|
| **our full QuantUMS model** | **0.28** |
| our model, bias term removed | 0.51 |
| plain SUM of the same fragment areas | 0.59 (within-run median 0.74) |
| top-3 fragments by score | 0.72 |

**Verdict — the test FAILS, and it refutes "the gap was just our XIC extraction."** Given DIA-NN's *own* inputs, our model reaches only r≈0.28 — **worse than trivially summing those same inputs (0.59)**. Two localized problems, #1 dominant:
1. **Our objective ≠ DIA-NN's objective.** Our inverse-variance combination drops a plain sum 0.59→0.51, and our **bias-correction term drops it further to 0.28** (the fit is dominated by `b[xic_score_sq]=+0.63`, a quality-dependent distortion that corrupts the cross-*precursor* abundance axis). Our two-term loss optimizes cross-*run* concordance/precision (hence the excellent CV) but trades away cross-precursor accuracy *even with perfect inputs*. So we reproduced QuantUMS's **structure + precision behavior, not its accuracy behavior**.
2. **Even a perfect model is capped ~0.74 by the exportable inputs.** DIA-NN 2.5's main report does not expose the real per-fragment quantities (`Fr.N.Quantity` is legacy 1.8.x TSV; the true values live in the closed `.d.quant` binary). The `--xic` chromatograms are a downstream visualization, so summing them tops out ~0.74.

**Net:** our reimplementation is a working, decoy/precision-sound QuantUMS-*like* estimator but **not a faithful DIA-NN clone**. Closing it would require the closed QuantUMS source / internal fragment columns AND a bias term anchored to **MS1-vs-MS/MS concordance** (the MS1 channel we lack), not fragment-vs-fragment alone. This is the honest endpoint of the QuantUMS thread.

Artifacts: `quantums_xic.py`, `compare_quant6.py`, `qxic_*` matrices, `quantums_xic_results.txt` (in `sage_patched/`); `xicextract/` (Rust extractor + 16 `xic_*.tsv`).

---

## 5. Sage config (v0.14.7 schema)

Example at `/quobyte/proteomics-grp/brett/glendon/sage_bigdog_diatracer_16/sage_bigdog.json`. Key fields:
```json
{
  "database": {
    "enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 50, "cleave_at": "KR", "restrict": "P"},
    "static_mods": {"C": 57.0215},
    "variable_mods": {"M": [15.9949], "^": [42.0106]},
    "max_variable_mods": 2,
    "decoy_tag": "rev_", "generate_decoys": true,
    "fasta": "/nfs/lssc0/flinders/proteomics/Data/raw_data/tTOF_HT/may26/Ameer_Taha/databases/UP000805418_full_isoforms_2026_06.fasta"
  },
  "precursor_tol": {"ppm": [-20, 20]},
  "fragment_tol": {"ppm": [-20, 20]},
  "isotope_errors": [-1, 3],
  "deisotope": true,
  "chimera": true,
  "min_matched_peaks": 6,
  "quant": {"lfq": true, "lfq_settings": {"peak_scoring": "Hybrid", "integration": "Sum", "spectral_angle": 0.7}}
}
```
**v0.14.7 quant-schema trap:** `quant.lfq` is a **boolean** + separate `quant.lfq_settings{peak_scoring,integration,spectral_angle}`. (v0.15+ removed the boolean and renamed the object `lfq` — emitting that shape against 0.14.7 errors with `invalid type: map, expected a boolean`.) Sage **generates its own `rev_` decoys** — give it a targets-only FASTA. For ID-only on diaTracer mzML, set `"lfq": false`.

---

## 6. Component 4 — quant comparison harness & results

**Script:** `compare_quant2.py` → `quant_comparison_results.txt` (in `sage_patched/`).
Reference = the **full bigDog DIA-NN 2.5 run** (215 files), `report.pr_matrix.tsv` at
`/nfs/lssc0/flinders/proteomics/Data/raw_data/tTOF_HT/may26/Ameer_Taha/Taha_Big_Dog_diann251_fulldb_mc2_acetyl_20260608/`.
Matched cross-engine on stripped sequence (charges collapsed) / bare accession.

**Precursor level (16 runs):**
| Method | # quantified | median CV | r vs DIA-NN |
|---|---|---|---|
| DIA-NN 2.5 (real QuantUMS) | 15,246 | 0.65 | — |
| Sage-LFQ (timsrust MS1) | 4,839 | 0.86 | 0.851 |
| Sage-QuantUMS-proto | 4,872 | **0.55** | 0.847 |

**Protein level:**
| Method | # quantified | median CV | r vs DIA-NN |
|---|---|---|---|
| DIA-NN 2.5 | 1,670 | 0.97 | — |
| Sage-LFQ | 755 | 1.04 | 0.832 |
| Sage-QuantUMS-proto | 864 | **0.71** | 0.762 |

(Sage-LFQ vs QuantUMS-proto: r=0.78 precursor, 0.84 protein.) **Finding:** QuantUMS-proto gives the lowest CV (consistent with its variance-minimization goal) while tracking DIA-NN's quant at r≈0.85. DIA-NN quantifies ~3× more features due to MBR + spectral library + 215-run context.

---

## 7. Data inputs (reference)

- **16 dog diaPASEF `.d`** (timsTOF HT, 60 SPD): paths in `/quobyte/proteomics-grp/brett/glendon/fragpipe_bigdog_diatracer_16/subset16.txt`. Full 215-file list: `…/Taha_Big_Dog_diann251_fulldb_mc2_acetyl_20260608/file_list.txt`.
- **diaTracer pseudo-mzML:** `<base>_diatracer.mzML` written **next to each `.d`** in `/nfs/lssc0/flinders/proteomics/Data/raw_data/tTOF_HT/may26/Ameer_Taha/`.
- **Dog DB:** `…/Ameer_Taha/databases/UP000805418_full_isoforms_2026_06.fasta` (43,670 seqs; Canis lupus familiaris + contaminants; Sage adds `rev_` decoys).
- **diaTracer binary** (closed/licensed, for generating pseudo-spectra): `/quobyte/proteomics-grp/fragpipe24/fragpipe-24.0/tools/diatracer-2.2.1.jar`. **Standalone CLI gotcha:** all numeric params must be explicit or it NPEs:
  ```bash
  java -jar diatracer-2.2.1.jar -d <.d> -t <threads> -w <tmpdir> -dI 0.01 -dR 3 -mC 0.3 -mF 1 -mO 0.1 -rM 500 -r 0
  ```
  Output `<base>_diatracer.mzML` lands next to the input `.d` (not in `-w`).

---

## 8. Licensing

| Component | License | Notes |
|---|---|---|
| Sage + the mzml patch | MIT | fully open; patch is ours |
| timsrust / ms1extract | open (Apache/MIT) | Sage already depends on timsrust |
| QuantUMS prototype | clean-room from the paper | not derived from DIA-NN code |
| diaTracer (pseudo-spectra) | **closed / academic license** | only used to generate input; clean-room reimpl attempted separately |

So everything *we built* is open; the only non-open dependency in this exact pipeline is diaTracer for pseudo-spectrum generation.

---

## 9. Known limitations & next steps

1. **No MBR / no raw-`.d` requant** → this is the whole ~984 (Sage) vs ~1,306 (DIA-NN 2.5) / ~1,337 (FragPipe) protein-group gap on the same 16 files. Single-shot Sage scores each run independently; DIA-NN/FragPipe transfer IDs across runs and requantify against the full raw DIA data.
   - **AlphaDIA — CORRECTED: it DOES work on timsTOF (earlier "ruled out" was wrong).** Per the AlphaDIA paper (Wallmann et al., Nat Biotechnol 2025), AlphaDIA works excellently on timsTOF dia-PASEF **with an empirical library** (Fig 2 ~6,800 PG; Fig 4c 7,649 PG, beating DIA-NN/Spectronaut). Our prior failure was specifically **fully-predicted-library directDIA WITHOUT transfer learning** — a known PeptDeep model-mismatch case; every timsTOF result in the paper uses an empirical library, and the predicted-library benchmarks are all Orbitrap. The documented fix for predicted-library on a mismatched instrument is **DIA transfer learning** (adapts PeptDeep: spectral corr 0.5→0.85, +48% precursors, FDR stays valid). v2.x config notes: `search.extraction_backend: python` (Rust backend doesn't support timsTOF). **RESULT (verified):** AlphaDIA empirical-library on all 16 files = **1,361 PG / 14,560 precursors @1% FDR** (healthy 87/13 target/decoy) — matches FragPipe, beats DIA-NN 2.5 (1,306). Caveat: the empirical library used was FragPipe's EasyPQP `library.tsv` (AlphaDIA rejected DIA-NN's `.parquet`), so the library provenance is FragPipe; AlphaDIA did the DIA search+quant. The fully-self-contained library-free path (transfer learning) is running but very slow (see scoreboard). Artifacts: `alphadia_correct_bigdog/`.

### 4e. Predicted-spectra rescoring of Sage's pin (MSBooster-equivalent, open) — closes the gap

The pin-feature diff (§ comparison) showed Sage lacks one thing vs MSFragger+MSBooster: a deep-learning **predicted-fragment-intensity spectral similarity**. We added it externally (AlphaPeptDeep `instrument='timsTOF'`, NCE=35 chosen by sweep; spectral-entropy similarity to the observed pseudo-spectra; appended to `results.sage.pin`; rescored with mokapot; decoys predicted from reversed sequences → valid FDR):

| level | Sage native | mokapot, no spectral | **mokapot + spectral** | FragPipe+MSBooster |
|---|---|---|---|---|
| protein groups | 984 | 1,170 | **1,406** | 1,361 |
| peptides | 6,402 | 6,874 | 8,235 | — |

`weighted_spectral_entropy` was mokapot's **2nd-strongest** feature (+2.63). **Same-rollup, valid-FDR contribution of the spectral feature alone = +236 PG.** (Cross-tool 1,406-vs-1,361 mixes protein-rollup conventions, so read as "competitive," not a hard win.) Bugs fixed en route: AlphaPeptDeep reorders its precursor_df (carry identity keys); **diaTracer mzML peak arrays are not m/z-sorted** (must sort before `searchsorted` or all matches zero out). Code: `sage_patched/sage_patched/{predict_ms2.py, build_feature.py, nce_sweep.py, rescore3.py}`.
   - **MBR DONE — and it revealed the gap is NOT an MBR gap.** Built `sage_mbr.py`: empirical run-library from Sage's confident IDs, LOWESS RT alignment, transfer at aligned RT+IM+m/z, **confirmed against raw `.d` via `xicextract`**, **transfer FDR validated by TWO independent decoy models** (RT-shuffle 0.83%, fragment-mismatch 0.91%, both ≤1%) → 2,644 FDR-held transfers. Result: mean per-run PG 463→**472**, but **union PG 996→996 (unchanged)**.
   - **The reframing (corrects an earlier mistake in this doc):** DIA-NN's 1,306 and FragPipe's 1,337 are **union** counts across 16 runs; DIA-NN's *mean per-run* PG is only ~850, Sage's union is already 996. So MBR raises per-run completeness toward Sage's *own* union but cannot exceed it. **675 of DIA-NN's 1,306 union PGs are identified in ZERO Sage runs** (631 shared, 365 **Sage-only**) → structurally unreachable by MBR-over-Sage. **The 984→1,300 gap is a SEARCH/IDENTIFICATION gap, not a match-between-runs gap.** The lever to close it is search sensitivity (diaTracer pseudo-spectrum quality + Sage scoring, or a library-based DIA pass — but the open library-DIA option, AlphaDIA, is blocked on timsTOF). Notably Sage is **partially orthogonal** to DIA-NN (365 unique PGs), so open-Sage + commercial together exceed either alone.
   - Artifacts: `sage_mbr.py`, `mbr_*.py`, `MBR_RESULTS.md`, `mbr_accepted_xic.tsv` in `sage_patched/`.
2. **LFQ not IM-resolved** — the merged MS1 carries no ion mobility in the precursor record. A fuller solution would propagate IM from `.d` MS1 into the precursor records.
3. **diaTracer dependency** — for a *fully* open stack, the clean-room diaTracer reimplementation (from Li et al. Nat Commun 16:95, 2025) needs to land; validate it against the real `_diatracer.mzML` outputs (spectrum counts, precursor distributions, peptide-ID overlap).
4. **QuantUMS prototype** is the estimator only — no learned hyperparameters; treat as a proof-of-concept, not a DIA-NN replacement.
5. **Sage v0.14.7** is current latest; if a newer Sage adds native DIA / MBR, re-evaluate (would obsolete much of this).

---

## 10. Quick reproduce checklist

1. `ssh` in (login shell). Ensure `/quobyte/proteomics-grp/brett/glendon/sage_patched/sage` exists (or rebuild from `sage_src/` after the Cargo.lock bump, on a compute node).
2. Generate pseudo-spectra: run diaTracer (jar above) per `.d` — or reuse the existing `*_diatracer.mzML` next to the `.d`.
3. **ID:** write a v0.14.7 `sage.json` (§5, `lfq:false`), run patched `sage --write-pin --output_directory <out> <sage.json> <*_diatracer.mzML>`. Read `results.sage.tsv` (filter `spectrum_q`/`peptide_q`/`protein_q` ≤ 0.01).
4. **LFQ:** `ms1extract` each `.d` → `merge_ms1_ms2.py` → run patched `sage` with `lfq:true` on the merged mzML → `lfq.tsv`.
5. **QuantUMS quant:** `quantums_proto.py` over `matched_fragments.sage.tsv`.
6. **Compare:** `compare_quant2.py` vs a DIA-NN `report.pr_matrix.tsv`.

All scripts/binaries live under `/quobyte/proteomics-grp/brett/glendon/{sage_patched,ms1extract,sage_bigdog_diatracer_16}/`.
