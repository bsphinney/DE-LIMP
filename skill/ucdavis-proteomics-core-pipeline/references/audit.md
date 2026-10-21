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
| `numeric_covariate` | — | a covariate fitted as a number (`run_de.R --numeric-covariates`) that Group explains most of: de_provenance.json `covariates` marks it `confounded` (R² of covariate ~ Group at or above `confound_r2_threshold`, 0.8). The message is run_de.R's own caution, quoted (composed from the term when the record has none): the comparisons rest on its fitted slope, how much wider each contrast's standard error is (`se_inflation`), and that adjusting for something the groups themselves change removes part of the effect being tested. PASS for a numeric covariate below the threshold; nothing when every covariate is a factor (that is `confounding`'s) |
| `instrument_mix` | — | >1 instrument model (batch effect) |
| `id_depth` | — | suspiciously few proteins quantified (< `--min-proteins`, default 500) → check params / organism FASTA |
| `missingness` | — | very high missing fraction (> `--max-missing`, default 0.5) |
| `contamination` | — | keratin/trypsin/casein contaminants present at a meaningful level |
| `contaminant_evidence` | — | (a DIA-NN report: parquet or tsv) a contaminant protein group every one of whose peptides in the DE's input report — rows kept as run_de.R keeps them (its recorded q-columns and cutoffs) — is also a target protein's: listed with it in `Protein.Ids`, or one of the peptides the build's digest makes from a target entry of the searched FASTA (`fetch_fasta.target_digest_peptides`, the sidecar's `contaminant_digest` — the rule that kept the entry; a string merely inside a target does not count: BSA's LCVLHEK sits in human ALB after Q, not K/R), as written or the same but for I/L (`fetch_fasta.il_equal`; no spectrum tells I from L). The data then cannot attribute the group to the contaminant's species: it is the sample organism's protein — its own, or handling contamination such as skin keratin — and the message says to report it so, never as that species' signal or an added standard. Horse myoglobin `Cont_P68082` in human heart rested on GHHEAELKPLAQSHATK (human MB's GHHEAEIKPLAQSHATK) and ASEDLKK (shared). The build keeps such an entry — it has 9 peptides of its own — so only the report can show it. Names each group, its peptides and the target protein, whether its precursors were removed, and the fix: the same `fetch_fasta.py` command plus `--drop-contaminant <accession>` for each, and a new search. PASS when every contaminant group has a peptide of its own; INFO when the report, pyarrow or the searched FASTA (checked by sha256) is not readable here |
| `contaminant_competition` | — | (a DIA-NN report) the sample genes the contaminant list took away, from `contaminant_competition.py`: the search's last job writes `contaminant_competition.json` and `search_provenance.json` `contaminant_competition`; the audit checks the DE's input report itself when the search left no record of it. A precursor is a contaminant's when any accession in its `Protein.Ids` carries the tag (DIA-NN's `--cont-quant-exclude` rule, which `run_de.R` applies too). A contaminant entry *shows itself* with at least 2 non-overlapping peptides that name no sample protein and are not a sample protein's digest peptide (I = L; the build's own count). **COMPETED OUT**: every precursor of a gene that passes the q-value filter is shared with entries that do not show themselves, so the gene has no quantity and is missing from the results; **SIPHONED**: it kept some, but the removed ones carry most of its signal. Each is a WARN naming the genes, the entries and the fix (the same `fetch_fasta.py` command plus `--drop-contaminant <accession>` for each, and a new search; the list itself stays the Core's universal one), never a FAIL. A gene seen only through peptides shared with an entry that shows itself (bovine albumin in a cell-culture sample) cannot be attributed: AUDIT.md lists it as not counted as lost. Also each run's contaminant share of the precursor signal. PASS: checked, none; INFO, with the reason, when it could not run |
| `de_signal` | — | 0 significant (underpowered) **or** >50% significant (batch/normalization/confounding artefact, not biology) |
| `quantification` | — | the DE pipeline's own caveat about its protein quantities (de_provenance.json `caveat`, set by the pipeline's descriptor), quoted. Today: a Sage DE, whose protein value is its single most intense peptide per run (no protein rollup model) |
| `normalisation_stability` | — | (a DIA-NN report) DIA-NN's per-run `Normalisation.Instability` from the `<report>.stats.tsv` beside the DE input, for the runs in `--conditions`, above `check_report_runs.NORM_INSTABILITY_WARN` (0.3 — DIA-NN documents no cutoff; the skill's sits between a cohort measured at 0.04–0.07 and one at 0.67–1.00 in every run). Names each run and value, and the alternatives DIA-NN documents (`--global-norm`, `--no-norm`, the report's non-normalised `Precursor.Quantity`). PASS with the median and max; INFO when the stats file or its column is missing. Every value is in `AUDIT.json` |
| `sage_lfq` | — | (`--search-out <search out dir>`, Sage only) a run's median precursor mass error is within 2 ppm of, or beyond, Sage's LFQ window (`quant.lfq_settings.ppm_tolerance`, default ±5 ppm), or Sage kept 0 or very few target MS1 peaks at 5% FDR. Identifications are not affected; the quantities are unreliable. The message is `sage_lfq_check.py`'s own record (`sage_lfq_check.json`), quoted as written and with the fix (a wider window, re-run Sage, `--adapt-only`); runs outside the window are also a gate at `--adapt-only` and `run_de.R`, so a DE on them means the user accepted them, which the message then says (who, when, why). PASS when the window fits; INFO when the check could not run |

## What the client reads: a plain-language table

The report's **Audit & caveats** section (`make_analysis_html.py`, from `AUDIT.json`; from
`AUDIT.md`'s lines for a session without it) is a table with one row per check: **Check**
(a plain name), **Result** (✅ Passed, ⚠️ Caution, ⛔ Problem, ℹ️ Note: the worst of that check's
findings) and **What it means for the results**. The words are written beside each check in
`audit_results.py` (`PLAIN`, and a finding's own `plain` where its meaning differs from the
check's usual one: no significant protein vs more than half, the normalisation decision
quoting Core staff's reason by role). A pipeline's own checks carry their label and words
(`run_silac.R`'s `chk(..., plain =)`, `SILAC_CHECK_LABELS`); a check that has none reads as a
generic caution pointing to AUDIT.md, never its technical message. `AUDIT.md` keeps every
check code, number and file name: it is the Core's full record.

A staff report (2.11.2): the section was AUDIT.md as written, check codes, B1-B4 normalisation
strings, Normalisation.Instability values and "implausibly high ... not real biology". The
staff member could not read it, and a client could take "not real biology" to mean the results
are wrong. A report's own `## Audit & caveats` section still replaces the table (write one only
when staff ask, and in plain words).

## How the orchestrator uses it
- Run it with `--conditions`, `--de-dir`, and the `detect_acquisition` JSON
  (`--acquisition-json`) so it can see acquisition/instrument mixing, and `--search-out`
  (the search's output dir) so a Sage search's LFQ check reaches the report.
- **STOP on any `FAIL`** — e.g. a singleton group or a confounded batch means the
  differential results aren't trustworthy; tell the user plainly and how to fix it
  (add replicates, de-confound the design, split by acquisition/instrument).
- **Relay every `WARN`** and fold the findings into the report.
- **A `contaminant_competition` WARN** names proteins the contaminant list took away. Tell the user; a re-search without those entries is Core staff's call (the director chose, 2026-10-20, to keep the universal list and report what it takes).

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
