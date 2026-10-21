# Results audit — catching common proteomics mistakes

New users make predictable proteomics mistakes that invalidate the results. The
skill runs `audit_results.py` after DE to catch them deterministically (grounded in
the real data — never guessed) and **surfaces every issue to the user** before they
over-interpret anything. FAILs should stop interpretation; WARNs go in the report's
"Audit & caveats" section.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 19: What it checks
- 37: What the client reads: a plain-language table
- 56: How the orchestrator uses it
- 69: Why deterministic, not LLM-judged
- 77: Tuning
- 83: 8c. The audit and the sample-quality check (step 8c)
<!-- /contents -->

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

> Superseded by "8c. The audit and the sample-quality check (step 8c)" below.

- Run it with `--conditions`, `--de-dir`, and the `detect_acquisition` JSON
  (`--acquisition-json`) so it can see acquisition/instrument mixing, `--search-out`
  (the search's output dir) so a Sage search's LFQ check reaches the report, and
  `--fasta-meta` (the searched FASTA's `.meta.json`).
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

## 8c. The audit and the sample-quality check (step 8c)

```
python3 scripts/audit_results.py --out AUDIT.md --conditions ./conditions.csv \
    --de-dir ./de_results --acquisition-json /tmp/acq.json --adjp 0.05 --logfc 1 \
    --fasta-meta ./search.fasta.meta.json --search-out ./search_out
```

Checks for the classic new-user pitfalls: too few/no replicates, imbalanced or
**confounded** design, **mixed acquisition or mixed instruments** in one analysis,
suspiciously low ID depth, very high missingness, contaminant dominance, DE
results that are too-empty or implausibly-large (batch/normalization artefacts), for a
Sage search (`--search-out`), an LFQ window that did not fit the runs' MS1 mass error
(`sage_lfq`), and, for a DIA-NN report, DIA-NN's own per-run `Normalisation.Instability`
(`normalisation_stability`, read from the `<report>.stats.tsv` beside the DE input).
**A `normalisation_stability` WARN** names every run above 0.3 (DIA-NN documents no cutoff;
the skill's is between cohorts measured at 0.04–0.07 and 0.67–1.00). The DE used DIA-NN's
normalised quantities, so tell the user, and offer a sensitivity analysis on a
re-quantification with the alternatives DIA-NN documents: `--global-norm` (one global factor
per run instead of RT-dependent normalisation) or `--no-norm` (none; `Precursor.Normalised`
and `PG.MaxLFQ` are then not normalised). `--global-norm` goes in the cfg on either route;
`--no-norm` goes in the cfg of a single-shot search, while the 5-step chain strips it from the
cfg and takes `diann_parallel.py --no-norm` instead (its report is `no_norm_report.parquet`).
Run the DE on both reports into separate output folders, compare, and
report which one the results come from — never swap one for the other silently.
**Surface every `WARN` to the user, and STOP on any `FAIL`** (e.g. a group with no
replicate, or a batch confounded with the biology) until they resolve it — don't
let a new user over-interpret a broken design. The findings also become the
report's "Audit & caveats" section, as a plain-language table (step 9). → detail:
`references/audit.md`.
`target_contaminants` lists quantified proteins whose sequence is also a common contaminant
— **KEPT in quantification**: a possible contamination source, reported rather than
excluded. Name any that are significant in DE, and judge them by type: skin/hair keratins
(KRT1/2/9/10, KRTAPs, FLG) are usually handling contamination; ACTB, EEF1A1, tubulins and
KRT8/18/19/7 in epithelial cells are usually endogenous. `contaminant_overlap` means real
proteins are missing from quant — rebuild the FASTA and re-search. `contaminant_evidence` names a
contaminant protein group whose every peptide is also one of a sample protein's digest peptides
(shared, or the same but for I/L): the data cannot attribute it to the contaminant's species.
Report it as the sample organism's protein — its own, or handling contamination (skin keratin) —
never as that species' signal or an added standard, and tell the user it was removed with the
contaminants. To quantify it as the organism's (the user's call: a new search), rebuild the FASTA
with the same `fetch` command plus the `--drop-contaminant` the finding names, and search again.
`sample_quality.py`'s CONTAMINANT_IDENTICAL panel uses the same list as `target_contaminants`; a
group-confounded result there is a STOP like any other panel. A `numeric_covariate` WARN (a
covariate fitted as a number that the groups nearly determine) quotes run_de.R's caution: tell
the user, and offer the DE without that covariate to compare.

Then run the **biological sample-quality** check (contamination that mimics biology —
hemolysis, tissue cross-contamination, skin) — and **read `references/anomaly-checks.md`
and walk its decision tree** (the judgment branches `audit_results.py` can't automate:
intensity-distribution outliers, replicate/condition ID overlap, mass/RT/CCS shifts,
p-value-distribution shape, largest-FC-in-low-abundance, log2FC>5 artefacts, PCA sample
swaps, GSEA background; plus XL-MS / phospho / non-model branches):

```
python3 scripts/sample_quality.py --matrix ./de_results/Expression_Matrix.csv \
    --conditions ./conditions.csv --report ./search_out/report.parquet --out SAMPLE_QUALITY.md \
    --fasta-meta ./search.fasta.meta.json
```

**If a contamination panel is confounded with a group, STOP** — DE may be contamination,
not biology, and protein-level filtering will not fix it (→ `references/anomaly-checks.md`).
Read the panels against the **stated** sample type (`sample_type.py show --session <session>`),
never one guessed from file names: muscle proteins are the analyte in a muscle sample and a
contaminant in a liver biopsy; hemolysis matters for plasma/serum. When the panels or the most
abundant proteins disagree with the stated type, ask the user — do not change it yourself. With
none stated (Core staff, `--not-asked`), ask what the samples are only when a panel flags
something whose meaning depends on it, and record the answer (`sample_type.py set --stated`).
**If the sample IS keratin** (step 3), both auditors read it from `--fasta-meta` (a database
built with `fetch_fasta.py --keratin-sample`); pass **`--keratin-sample`** to `sample_quality.py`
**and** `audit_results.py` only for a database that does not record it — keratin is the analyte
there, not a contaminant, and must not be flagged. This is a check, not the fix: a keratin
sample must be decided at step 3, before the database is built.
Flag every anomaly in the 5-part format (*what · where · why · likely cause · fix*); the
report **always** gets a **Data Quality Notes** section, even if "nothing anomalous observed."
