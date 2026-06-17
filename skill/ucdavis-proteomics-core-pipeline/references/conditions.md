# Experimental conditions — collection & mapping

The user should never have to hand-fill a grid. They give conditions however is
easiest; the skill maps them onto the real raw filenames and asks only about
genuine ambiguities.

## Two ways the user provides conditions
1. **In words** — "the first three are control, the last three treated", "A* are
   wild-type, B* are knockout". The agent reads the real run list and turns this
   into an intent JSON: `{"groups": {"control": [...], "treated": [...]}}` or
   `{"mapping": {"<sample>": "<group>"}}`.
2. **A file** — any CSV/TSV they already have. Column names are auto-detected:
   sample column from {File.Name, filename, run, sample, sample name, name, raw,
   id}; group from {group, condition, treatment, class, type, cohort, phenotype};
   batch from {batch, block, plate, run order}. Up to two further columns become
   Covariate1/Covariate2.

## How mapping works (`collect_conditions.py --map`)
Matching is grounded in the actual run names (never guessed):
1. exact match (case-insensitive, punctuation-insensitive), else
2. substring either way (so `Treated` matches `Treated_rep1..3`, and a full
   filename matches its run name).

It then reports, for confirmation:
- `unassigned_runs` — a raw file no condition matched. **Blocking**: every run
  needs a group or it can't be analyzed.
- `conflicting_runs` — a file matched to >1 group. **Blocking**: resolve it.
- `unmatched_identifiers` — the user named something with no matching file (typo or
  a sample that isn't in this batch).
- `multi_match_identifiers` — one label hit several files. Usually intended (a
  replicate prefix), but confirm it's what they meant.
- `singleton_groups` — a group with <2 replicates → no within-group variance, so
  no DE for it. Warn.

The agent confirms each with the user, finalizes `conditions.csv`, then runs
`--validate` against the search report before DE.

## The finished metadata
`conditions.csv` columns: `File.Name,Group[,Batch,Covariate1,Covariate2]`.
`File.Name` must equal the Run names in the search report. `--validate` checks
column presence, blank groups, singleton groups, and that the report runs and
metadata rows line up exactly.

## Why deterministic matching (not pure LLM)
Filename↔condition matching is where an LLM can hallucinate (assigning a group to a
file that doesn't exist, or mismatching near-identical names). Keeping the match in
`collect_conditions.py`, grounded in the real run list, means the agent does the
language understanding while the assignment is verifiable — and every uncertainty
is surfaced for explicit confirmation rather than guessed.
