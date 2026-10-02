# Results audit — catching common proteomics mistakes

New users make predictable proteomics mistakes that invalidate the results. The
skill runs `audit_results.py` after DE to catch them deterministically (grounded in
the real data — never guessed) and **surfaces every issue to the user** before they
over-interpret anything. FAILs should stop interpretation; WARNs go in the report's
"Audit & caveats" section.

## What it checks
| check | FAIL | WARN |
|---|---|---|
| `replication` | a group with <2 replicates (no within-group variance — stats invalid) | a group with exactly 2 (low power) |
| `group_balance` | — | group sizes very unequal (≥3× ratio) |
| `confounding` | a covariate (Batch) perfectly confounded with Group — biology can't be separated from batch | — |
| `acquisition_mix` | DIA and DDA files mixed in one analysis | — |
| `instrument_mix` | — | >1 instrument model (batch effect) |
| `id_depth` | — | suspiciously few proteins quantified (< `--min-proteins`, default 500) → check params / organism FASTA |
| `missingness` | — | very high missing fraction (> `--max-missing`, default 0.5) |
| `contamination` | — | keratin/trypsin/casein contaminants present at a meaningful level |
| `de_signal` | — | 0 significant (underpowered) **or** >50% significant (batch/normalization/confounding artefact, not biology) |
| `quantification` | — | the DE pipeline's own caveat about its protein quantities (de_provenance.json `caveat`, set by the pipeline's descriptor), quoted. Today: a Sage DE, whose protein value is its single most intense peptide per run (no protein rollup model) |
| `normalisation_stability` | — | (a DIA-NN report) DIA-NN's per-run `Normalisation.Instability` from the `<report>.stats.tsv` beside the DE input, for the runs in `--conditions`, above `check_report_runs.NORM_INSTABILITY_WARN` (0.3 — DIA-NN documents no cutoff; the skill's sits between a cohort measured at 0.04–0.07 and one at 0.67–1.00 in every run). Names each run and value, and the alternatives DIA-NN documents (`--global-norm`, `--no-norm`, the report's non-normalised `Precursor.Quantity`). PASS with the median and max; INFO when the stats file or its column is missing. Every value is in `AUDIT.json` |
| `sage_lfq` | — | (`--search-out <search out dir>`, Sage only) a run's median precursor mass error is within 2 ppm of, or beyond, Sage's LFQ window (`quant.lfq_settings.ppm_tolerance`, default ±5 ppm), or Sage kept 0 or very few target MS1 peaks at 5% FDR. Identifications are not affected; the quantities are unreliable. The message is `sage_lfq_check.py`'s own record (`sage_lfq_check.json`), quoted as written and with the fix (a wider window, re-run Sage, `--adapt-only`); runs outside the window are also a gate at `--adapt-only` and `run_de.R`, so a DE on them means the user accepted them, which the message then says (who, when, why). PASS when the window fits; INFO when the check could not run |

## How the orchestrator uses it
- Run it with `--conditions`, `--de-dir`, and the `detect_acquisition` JSON
  (`--acquisition-json`) so it can see acquisition/instrument mixing, and `--search-out`
  (the search's output dir) so a Sage search's LFQ check reaches the report.
- **STOP on any `FAIL`** — e.g. a singleton group or a confounded batch means the
  differential results aren't trustworthy; tell the user plainly and how to fix it
  (add replicates, de-confound the design, split by acquisition/instrument).
- **Relay every `WARN`** and fold the findings into the report.

## Why deterministic, not LLM-judged
These are exactly the checks where an over-eager interpreter would hand a new user a
confident but wrong story (e.g. "1,800 proteins changed!" when a batch effect made
half the proteome "significant"). Grounding them in counts from the real files keeps
the skill honest and protects the user from the default-workflow violations that
matter most. It complements DE-LIMP's architectural rules (no fabricated values; the
pipeline self-describes) and its "Error handling & UX audit" review.

## Tuning
`--adjp` should match the DE significance cutoff and `--logfc` the volcano's reference
line (it is counted descriptively, never used to call significance); `--min-proteins` and
`--max-missing` are soft thresholds — adjust for very small samples or enriched
fractions (e.g. a phospho or secretome experiment legitimately quantifies fewer).
