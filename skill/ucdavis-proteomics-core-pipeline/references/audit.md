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

## How the orchestrator uses it
- Run it with `--conditions`, `--de-dir`, and the `detect_acquisition` JSON
  (`--acquisition-json`) so it can see acquisition/instrument mixing.
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
`--adjp` / `--logfc` should match the DE thresholds; `--min-proteins` and
`--max-missing` are soft thresholds — adjust for very small samples or enriched
fractions (e.g. a phospho or secretome experiment legitimately quantifies fewer).
