# ucdavis-proteomics-core-pipeline — open defects

**Audited:** 2026-08-17 against `c29f269` (`main`, immediately after PR #48 merged).
**Scope:** `skill/ucdavis-proteomics-core-pipeline/`. Every claim below was checked against the
tree, not recalled — file and line references are from that commit, and the audited files are
byte-identical between `main` and the `skill-v2.0.0` branch, so both carry what follows.

This file is for **open** defects. `docs/GOTCHAS.md` is the quick reference for problems that
already have solutions; when something here is fixed, delete its entry and, if the lesson is
reusable, add a row there instead. A shrinking file is the point.

---

## Status at a glance

**No open defects.** Every entry from the 2026-08-17 audit has been fixed; the fixed list is
below and in `docs/GOTCHAS.md`. When a new one is found, add a row here.

---

## ~~1. The precursor *m/z* range is hardcoded~~ — FIXED 2026-08-17

Fixed in `fix/skill-precursor-mz-range`. `detect_acquisition.py` now returns the acquired
bounds it was already reading (`precursor_mz_range`), and `estimate_params.py` takes
`--precursor-mz-range LO HI` and rounds outward. When no range is available the old 380–980 is
still used but tagged **`FALLBACK — acquired range unknown for this input; NOT measured`**
rather than `universal trypsin/LFQ default`, so a reader can tell a measurement from a guess.

Verified against a real `.d` on HIVE: acquired 299.5–1200.5 → emits `--min-pr-mz 299`
/ `--max-pr-mz 1201`. Covered by `skill/ucdavis-proteomics-core-pipeline/tests/`, now run in CI.

Kept here as a stub only until the next edit of this file; the lesson is in `docs/GOTCHAS.md`.

---

## ~~2. The q-value column set is duplicated~~ — FIXED 2026-08-17

Fixed in `fix/skill-q-column-definition`. The set now lives in
`scripts/diann_q_columns.py`, mirrored by `scripts/diann_q_columns.R` (the two R scripts are
standalone Rscripts that treat `jsonlite` as optional, so reading a shared JSON file would put
a hard dependency in front of the identification filter). `tests/test_q_columns.py` parses the R
source and asserts it equals the Python values, so drift is a CI failure rather than two subtly
different FDR filters.

The audit undercounted: `run_search.py` alone held **six** hand-written copies, not one — four
of them in output adapters (AlphaDIA, Sage, FragPipe DIA, FragPipe combined_protein,
Radiant/Fulcrum) that emit the DE contract. A test asserts no call site restates the set inline,
which is what found them.

It also conflated **two distinct concepts**, and merging them would have been a bug:

* `FDR_REQUIRED` + `FDR_OPTIONAL` — the **filter set**: every column ANDed at the q-cutoff
  (`run_de.R`, `build_maxlfq.R`).
* `PROTEIN_Q_PREFERENCE` — a **preference chain**: the first available protein-level q-column
  wins and the rest are ignored (`compare_searches.py`, `run_search.py`).

Applying the preference order as a filter would over-filter; filtering on only the first
available column would under-filter.

---

## ~~3. A single `--q-cutoff` for all six columns~~ — FIXED 2026-08-17

`COLUMN_CUTOFFS` in `scripts/diann_q_columns.py` (mirrored in the `.R`) gives `PG.Q.Value`
DIA-NN's own recommended **0.05**; the other five keep `--q-cutoff`. `limpa::readDIANN()`
recycles `q.cutoffs` against `q.columns` element-wise, verified in
`EListFromLongFormatFile`, so the dpc path passes a vector; `build_maxlfq.R` applies the same map
and labels any column whose cutoff differs as `PG.Q.Value@0.050` in `filters_applied`.

**This is a behaviour change that LOOSENS identification FDR.** Measured on a 2-run HeLa report:
32,559 → 32,741 rows and 2,695 → 2,715 protein groups. `cutoff_for(..., uniform=True)` restores
one cutoff for all six.

Deliberately *not* clever about a tightened `--q-cutoff`: the pipeline default (0.01) is already
stricter than 0.05, so a "never loosen what the caller asked for" rule would stop the
recommendation ever applying. `--q-cutoff` governs the other five.

---

## ~~4. `--min-pr-charge 2` narrows DIA-NN's default without saying why~~ — FIXED 2026-08-17

Tag only; no behaviour change. It now reads: *"z=1 excluded: rarely informative for tryptic
bottom-up. DIA-NN's own default is 1-4; this drops ~19% of the predicted library (measured:
10,899 → 8,805 precursors on a 60-protein FASTA)"*. The narrowing is still applied — it is the
right call — but a reader of the emitted parameters can now see that it was a decision.

---

## Fixed since the last audit — recorded so they are not re-litigated

All verified present in `c29f269`:

- **Experiment-wide FDR columns on the dpc path** — PR #48. `--method dpc` treated
  `Global.Q.Value` / `Global.PG.Q.Value` / `PG.Q.Value` as fallbacks and so never applied them on
  any modern DIA-NN report. Worth 227 protein groups (median one precursor, present in 9.4% of
  runs) on a 373-run report.
- **Executable bits.** Every directly-invoked `scripts/*.sh` is `100755` in the index.
  `diann_release.sh` is `100644` and correctly so — it is `. sourced`, never executed.
- **`--requeue`** is emitted for generated SLURM steps, so a preempted array task on a
  preemptible partition no longer breaks an `afterok` chain hours in.
- **The step-4 shared-library write race.** Each array task now gets a private library copy;
  previously all tasks re-saved the same file and the losers blocked indefinitely.
- **The glob hazard** in generated commands — `shlex.split` / `shlex.quote` now handle
  `--cut K*,R*` and `--var-mod ...,*n`, which would otherwise be expanded by the shell against
  the working directory.
- **The deprecated one-stage DIA-NN invocation.** `run_search.py` detects the library-free flag
  combination, strips it, and emits two dependent SLURM jobs — the two-stage workflow DIA-NN 2.3+
  requires.
