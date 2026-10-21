# Experimental conditions — collection & mapping

The user should never have to hand-fill a grid. They give conditions however is
easiest; the skill maps them onto the real raw filenames and asks only about
genuine ambiguities.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 20: Two ways the user provides conditions
- 54: The Core LIMS sample sheet
- 84: A run whose file name holds `--`
- 97: How mapping works (`collect_conditions.py --map`)
- 138: The finished metadata
- 174: Why deterministic matching (not pure LLM)
- 181: Less-used flags
- 186: Step 3: collecting the conditions and settling every ambiguity
- 246: Step 8: check the design against the report; runs left out of the DE
<!-- /contents -->

## Two ways the user provides conditions

> Superseded by "Step 3: collecting the conditions and settling every ambiguity" below.

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

> Superseded by "Step 3: collecting the conditions and settling every ambiguity" below.

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

## A run whose file name holds `--`

> Superseded by "Step 3: collecting the conditions and settling every ambiguity" below.

DIA-NN cannot be given such a path, so `run_search.py --rename-dashed-runs` searches the run
under a new name, each run of dashes made one (`x_DIA--sampA_S3` → `x_DIA-sampA_S3`), and that is
its `Run` in the report. `detect_acquisition.py` lists these runs at step 2 (`dashed_run_names`).
Build the design with the same names: `--rename-dashed-runs` on `--list-runs`, `--map` and
`--emit-template` (with `--from-dir` or `--runs`; a `--from-report` already has them), or on
`core_submission.py conditions`. Matching ignores punctuation, so a user's file naming the raw
name still finds the run. `--validate --against <report>` names any row that still holds the raw
name (`x_DIA--sampA_S3 -> x_DIA-sampA_S3`).

## How mapping works (`collect_conditions.py --map`)
Matching is grounded in the actual run names (never guessed):
1. exact match (case-insensitive, punctuation-insensitive), else
2. the identifier as a whole **token** of the run name (`LRS-96` names `..._DIA-LRS-96_S3`), or
   the run name inside the identifier (a full file name names its run). `-`, `_` and spaces are
   token boundaries, so a delimiter in the identifier counts: `DIA-XY1_` names `DIA-XY1_S3…`
   and never `DIA-XY10_…` (2.9.1 dropped the `_` and matched XY10–XY19 too). A separator may
   come or go only where letters meet digits (`XY-1` names `XY1`, `LRS96` names `LRS-96`). This
   is the same rule `core_submission.py locate` matches sample ids with. Only when no run carries
   the identifier as a token:
3. its words as whole words of the run name, where a word of letters also names that word plus
   digits (`Ctrl` names `Ctrl1..3`, a prefix). A word with digits always matches only a whole
   word (`KG1` never names `KG12`, `A1` never `A10`).

Never a bare substring, and **never the plate position**: the slot-well words of a run name
(`_S3-A1_`) say where a tube sat on the Core's plate, so a sample labelled `A1` does not match
the run in well A1.

> Superseded by "Step 3: collecting the conditions and settling every ambiguity" below.

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

**Is a covariate a measurement?** `run_de.R` fits every covariate as a factor (one level per
value) unless told otherwise. When `Covariate1`/`Covariate2` holds a continuous measurement (a
bait intensity, an age, a run depth), ask the user, and run DE with `--numeric-covariates
Covariate1` (`references/de-analysis.md`, "Continuous covariates"). Never decide it from the
values: a batch numbered 1/2/3 is a factor. A numeric covariate needs a number for every sample.

> Superseded by "Step 3: collecting the conditions and settling every ambiguity" below.

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

> Superseded by "Step 8: check the design against the report; runs left out of the DE" below.

`File.Name` must equal the Run names in the search report. `--validate` checks
column presence, blank groups, singleton groups, and that the report runs and
metadata rows line up exactly. A report run the user decided to leave out of the DE
(a re-run, a failed injection) has no row; record it with `--validate ... --against
<report> --exclude <run> --reason '<why>'` (one `--reason` right after each
`--exclude`). It then passes and goes into `<csv>.decisions.json` under
`excluded_runs`; `run_de.R` carries the reasons into `de_provenance.json`
(`runs_left_out`), and the Methods and the report list them. An excluded run must be
in the report and must not have a row. Any other missing run still fails.

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

## Step 3: collecting the conditions and settling every ambiguity

- **Conditions:** ask the user to either *tell you* the conditions in plain words
  ("the first three are control, the rest treated") **or** *upload a file* (any
  CSV/TSV with a sample column and a group column, however named, or an Excel `.xlsx` —
  the Core's LIMS export `PROT_####.samples.<date>.xlsx` included, read without openpyxl).
  Don't make them hand-fill a template.

First get the real run names the conditions must map to (add `--rename-dashed-runs` here and to
`--map` when step 2 listed `dashed_run_names`, so the CSV carries the names the search will use):

```
python3 scripts/collect_conditions.py --list-runs --from-dir /path/to/raw --glob '*.d'
```

Then map their conditions onto those runs — the script does the fuzzy filename
matching so you don't guess:

```
# (a) user uploaded a file:
python3 scripts/collect_conditions.py --map conditions.csv \
    --from-dir /path/to/raw --glob '*.d' --from-file <their_file>
# (b) user described it in words: turn it into intent, then ground it against the runs:
python3 scripts/collect_conditions.py --map conditions.csv \
    --from-dir /path/to/raw --glob '*.d' \
    --mapping-json '{"groups": {"control": ["..."], "treated": ["..."]}}'
```

Read the returned `ambiguities` and **confirm every one with the user** before
proceeding — `unassigned_runs` (a raw file no condition matched), `conflicting_runs`
(a file matched to two groups), `unmatched_identifiers` (the user named something
with no matching file), and `singleton_groups` (<2 replicates → no within-group variance).
**`decisions_required` is a sample-identity question, and only the user answers it** — a label
that names several runs (`multi_match_identifiers`: two runs numbered 61, re-injections) or
that two samples share (`duplicate_identifiers`). Ask each one as its own question, with its
runs listed, **before** the compute confirmation and never folded into it; never settle one
yourself, not even by another key that looks safe (a plate well, a run number) — propose that
key in the question instead. Those runs are written with **no group** until the answer comes,
so `--validate` refuses the CSV: re-run `--map` with `--confirm-multi '<label>'` when the user
says the label names every one of its runs (a group label such as `ctrl`), or with a mapping
that names the one run. The answer is recorded in `<csv>.decisions.json` (copied with the CSV
by `provenance.py`); also log it (`log_decision.py`). Do **not** start a search while any run is
unassigned, conflicting or awaiting a decision. A subject column (Mouse, Animal, Subject, Patient, Donor…) whose values look
like subjects (≥ 3, not the groups relabelled, recurring across groups) is kept under its
own name; when the output has `"block_suggested": true`, run DE with
`--block <block_column>` (step 8). A subject-named column whose values do not (Subject =
M/F, Patient = Yes/No) stays a covariate and comes back as `subject_ambiguous` — ask, and
re-run with `--subject-column <header>` only if the user confirms it. A condition column that
names each **replicate** (the LIMS `condition_name`: `X_mix_1` … `X_mix_5`) is read as conditions
plus replicate numbers in the proposed CSV (`replicate_labels`; `conditions.csv` keeps the label
in `Label`) when that is unambiguous, and **always asked** (`replicate_labels_to_confirm`: the
numbers may be time points or doses); a doubtful reading is not applied
(`replicate_labels_ambiguous`, with why). Ask, then re-run with `--replicate-labels collapse` (yes)
or `keep` (no). Two identifier columns naming as many runs (`sample_column_tie`) are asked too:
re-run with `--sample-column <header>`.
Finalize the CSV; it is
checked against the search's own run names in step 8, once `report.parquet` exists.
(Fallback: if the user has nothing yet, `--emit-template` writes a blank
File.Name,Group sheet for them to fill.)

## Step 8: check the design against the report; runs left out of the DE

First check the design against the runs the search actually reported — every report run needs
a row and every row a report run, every row a group, every group two samples — so a renamed or
dropped file fails here, not halfway through the DE:

```
python3 scripts/collect_conditions.py --validate conditions.csv --against ./search_out/report.parquet
```

**A searched run left out of the DE on purpose** (a re-run whose first injection is analysed, a
failed injection: one Core submission searched 13 runs and analysed 11) has no row, so `--validate` names it
under "report runs missing from metadata" and exits 1. **Never go past that failure on your own
judgment.** Ask the user; only when they decide the run is left out, record it, one `--exclude`
per run, each with its own `--reason` in their words (it goes into the Methods and the report,
so it names no person):

```
python3 scripts/collect_conditions.py --validate conditions.csv --against ./search_out/report.parquet \
    --exclude <run> --reason '<why>' [--exclude <run> --reason '<why>' ...]
```

(The reason in single quotes; through `hive_exec.sh`, in double quotes without `$` `` ` `` `\` `"` `'`: step 0a.)
The excluded runs then pass and are written to `<csv>.decisions.json` (`excluded_runs`, copied
into the bundle with the CSV). `run_de.R` copies each run and its reason into
`de_provenance.json` (`runs_left_out`) and `methods.txt`, and the Methods and the report list
them. A run left out with no recorded reason is listed as "not recorded", never given a reason
afterwards. An excluded run must be one the report has, and must not keep a row: `run_de.R`
analyses every row. A missing run without `--exclude` still fails. Also log the decision
(`log_decision.py`).
