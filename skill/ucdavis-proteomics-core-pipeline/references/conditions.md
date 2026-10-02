# Experimental conditions — collection & mapping

The user should never have to hand-fill a grid. They give conditions however is
easiest; the skill maps them onto the real raw filenames and asks only about
genuine ambiguities.

## Two ways the user provides conditions
1. **In words** — "the first three are control, the last three treated", "A* are
   wild-type, B* are knockout". The agent reads the real run list and turns this
   into an intent JSON: `{"groups": {"control": [...], "treated": [...]}}` or
   `{"mapping": {"<sample>": "<group>"}}`. If they said which animal / patient each
   sample came from, add `"subjects": {"<sample>": "<subject>"}, "subject_column": "Mouse"`.
2. **A file** — any CSV/TSV they already have, or an Excel `.xlsx` / `.xlsm` (read with the
   standard library — the pipeline env has no openpyxl; `--sheet <name>` picks a worksheet,
   else the first one that is not hidden, and when it has no sample and group columns the error
   names the sheet read and the others; a merged range reads as Excel shows it, its value in
   every cell (`merged_ranges_filled`); a broken workbook is a clean error, never a traceback;
   an old binary `.xls` is refused: save it as `.xlsx` or CSV). Column names are auto-detected:
   sample column from {File.Name, filename, run, sample, sample name, name, raw,
   id, unique id, sample id} — with several, the one whose values name the most runs
   (`sample_column`, `sample_column_matches`); a **tie is asked** (`sample_column_tie`: no run is
   given a group until `--sample-column <header>` names the column — the sheet's order never
   decides); group from {group, condition, treatment, class,
   type, cohort, phenotype, condition name, group name};
   batch from {batch, plate, run order}; a subject column from {mouse, mice,
   animal, rat, subject, patient, donor, individual, participant, pair, block} (an `id` / `no` /
   `number` suffix allowed: `Mouse ID`, `animal_no`) is kept **under its own name**
   (`Mouse ID` → `Mouse_ID`) — but only when its VALUES look like subjects: ≥ 3 of them,
   not the groups relabelled, and (where they span groups) recurring across them. A
   `Subject` of M/F, a `Patient` of Yes/No or a `Donor` identical to the group stays a
   covariate and is reported as `subject_ambiguous` for the user to confirm
   (`--subject-column <header>`; `--subject-column none` turns detection off). Up to two
   further columns become Covariate1/Covariate2 (`covariate_columns` says which header
   went where). Anything that does not fit is NAMED, never dropped silently:
   `columns_not_written` lists extra columns beyond the two slots, and an ambiguous subject
   column with no free slot says `"written": false` and how to keep it.

## The Core LIMS sample sheet

`PROT_####.samples.<date>.xlsx` (columns `internal_id, internal_notes, sample_name, unique_id,
condition_name, amt_2_inject`) maps as it is: `unique_id` or `sample_name`, whichever names the
runs, is the sample column; `condition_name` the group. `internal_id` (the submission number on
every row), `internal_notes` and `amt_2_inject` are bookkeeping, never covariates: they are listed
in `columns_skipped` with why, as is the identifier column not used.

**`condition_name` names each replicate** (`X_mix_1` … `X_mix_5`): taken literally it is one
singleton group per sample and no DE at all. When every sample has its own label, `--map` reads
the column as **conditions + replicate numbers** — the condition is the label without its
trailing number (after a separator such as `_`, `-`, `.` or a space, or straight after a letter:
`Control1`; a `rep` word before the number goes too) — but **only when that is unambiguous**:
every label ends in a number, they make ≥ 2 conditions of ≥ 2 samples, no two claim the same
replicate of a condition, no two conditions differ only in case or punctuation, and each
condition's numbers run 1..n (or 0..n−1). Then the **proposed** csv has `Group` = the condition,
`Label` = the label as given (`make_figures.R` names samples by it), `replicates` gives each
run's number and `replicate_labels` records the reading — and it is still **asked**
(`replicate_labels_to_confirm`, `needs_confirmation`): numbers can be the conditions themselves,
and `Day1`–`Day3` beside `Ctl1`–`Ctl3` passes every test above as a time course. Otherwise
nothing is collapsed and
`replicate_labels_ambiguous` carries the proposal (`conditions`) and the `reasons` (a gap — a
missing replicate, or not replicate numbers; numbers that are the conditions themselves, such
as `Day1` … `Day10`; …): **ask the user**, then re-run with `--replicate-labels collapse` (apply
the proposal) or `--replicate-labels keep` (the labels are the groups). The same reading applies
to `--mapping-json` intent, which is how `core_submission.py conditions` hands CoreOmics
condition names to `--map`.

## How mapping works (`collect_conditions.py --map`)
Matching is grounded in the actual run names (never guessed):
1. exact match (case-insensitive, punctuation-insensitive), else
2. whole words of the run name, in order (`LRS-96` names `..._DIA-LRS-96_S3`), or the run
   name's words inside the identifier (a full file name names its run). Never a bare substring:
   a word with digits matches only a whole word (`KG1` never names `KG12`, `A1` never `A10`); a
   word of letters names that word plus digits (`Ctrl` names `Ctrl1..3`). **Never the plate
   position**: the slot-well words of a run name (`_S3-A1_`) say where a tube sat on the Core's
   plate, so a sample labelled `A1` does not match the run in well A1.

It then reports, for confirmation:
- `unassigned_runs` — a raw file no condition matched. **Blocking**: every run
  needs a group or it can't be analyzed.
- `conflicting_runs` — a file matched to >1 group. **Blocking**: resolve it.
- `unmatched_identifiers` — the user named something with no matching file (typo or
  a sample that isn't in this batch).
- `multi_match_identifiers` — one label hit several files, and `duplicate_identifiers` — one
  label given for two samples. **Blocking until the user answers**: which run is which sample
  is theirs to say, so the runs are written with no group and listed under
  `awaiting_decision_runs`, and `decisions_required` holds each question with its runs. Ask
  every one, even when another key (a plate well, a run number) seems to settle it — put that
  key in the question. A label the user says names every one of its runs (a replicate prefix
  such as `ctrl`): re-run `--map --confirm-multi '<label>'`; it is recorded in
  `<csv>.decisions.json`, which `provenance.py` copies beside the CSV. One sample: name its run
  in the mapping. Two samples sharing a label: map each by something only its run carries.
- `singleton_groups` — a group with <2 replicates → no within-group variance, so
  no DE for it. Warn.

The agent confirms each with the user, finalizes `conditions.csv`, then runs
`--validate` against the search report before DE.

## The finished metadata
`conditions.csv` columns: `File.Name,Group[,Batch,Covariate1,Covariate2][,<block column>]`.

**Samples from the same animal / patient?** (five IPs from one mouse brain, paired
before/after samples.) Ask. If so, keep that unit in its own column under its own name
(`Mouse`, `Patient`) and run DE with `--block Mouse` (`references/de-analysis.md`,
"Paired / repeated designs"). `run_de.R` then fits it as a **fixed** effect when it is
crossed with the groups and every contrast is within one subject (before/after in the same
patient: the exact paired analysis), and as a **random** effect when it is nested in a
group (mice within age). `--map` keeps a subject column under its own name and reports it
as `block_column`; `block_suggested` is true when every run has a subject, every subject
holds ≥ 2 runs (all singletons = a sample id, not a block) and the values look like
subjects or were confirmed (`subject_confirmed`). Ambiguities add
`subject_conflicting_runs`, `runs_without_subject` and `single_run_subjects`. run_de.R
prints a note when a column recurs across groups and no `--block` was given.
`File.Name` must equal the Run names in the search report. `--validate` checks
column presence, blank groups, singleton groups, and that the report runs and
metadata rows line up exactly.

## Why deterministic matching (not pure LLM)
Filename↔condition matching is where an LLM can hallucinate (assigning a group to a
file that doesn't exist, or mismatching near-identical names). Keeping the match in
`collect_conditions.py`, grounded in the real run list, means the agent does the
language understanding while the assignment is verifiable — and every uncertainty
is surfaced for explicit confirmation rather than guessed.

## Less-used flags
- `--from-report <report.parquet|.tsv>` takes the run names from a search report, and
  `--runs "run1,run2,..."` gives them directly, instead of `--from-dir` + `--glob`.
- `--emit-template ... --covariates Batch,Covariate1`: add covariate columns to the blank sheet.
