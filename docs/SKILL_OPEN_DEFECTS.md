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

| # | Defect | Severity | Evidence |
|---|---|---|---|
| 1 | The q-value column set is hand-written in five files and they disagree | Medium — this is how defect #48 drifted in | five files, see below |
| 2 | One `--q-cutoff` applied to all six q-columns | Low | `scripts/run_de.R`, `scripts/build_maxlfq.R` |
| 3 | `--min-pr-charge 2` narrows DIA-NN's own default with no stated reason | Low — decision, not a bug | `scripts/estimate_params.py:181` |

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

## 2. The q-value column set is duplicated across five files, and they disagree

**What.** Which DIA-NN q-value columns constitute "FDR filtering" is written out by hand in
five places, with three different answers:

| File | Columns |
|---|---|
| `scripts/run_de.R:166-167` | `Q.Value` / `Lib.Q.Value` / `Lib.PG.Q.Value` + `Global.Q.Value` / `Global.PG.Q.Value` / `PG.Q.Value` |
| `scripts/build_maxlfq.R:29, 48` | the same six |
| `scripts/run_search.py:554` | `col("Global.PG.Q.Value", "Q.Value")` — two, as a fallback chain |
| `scripts/compare_searches.py:72` | `col("Global.PG.Q.Value", "Global.Q.Value", "Lib.PG.Q.Value", "Q.Value")` — four |
| `scripts/radiant_to_delimp.py:180` | maps `Global.PG.Q.Value` onto `Lib.PG.Q.Value` |

The last three are fallback chains for reading a single q-value out of a report rather than full
FDR filters, so they are not all the same kind of thing — but that is precisely the problem: a
reader cannot tell which copy is authoritative, and there is no single definition to consult.

**Why it matters.** This is the mechanism behind **PR #48**. `run_de.R`'s dpc path and
`build_maxlfq.R` disagreed about the same report for months because the definition lived in two
hand-maintained copies and only one was updated. Since then the count has grown to five. The
same failure mode produced a second, unrelated defect in a downstream session, where a control
script copied from its sibling kept the sibling's filter description.

**Proposed fix.** One exported definition — column names, the cutoff each takes, and a one-line
description of what each controls — imported by both R scripts, with the Python readers pointing
at it in a comment even where they cannot import it. The point is that a future change happens
once and is reviewable in one place.

---

## 3. A single `--q-cutoff` is applied to all six q-columns

**What.** Both paths apply `--q-cutoff` (default 0.01) uniformly. DIA-NN's README recommends
`PG.Q.Value` **"at 0.01 to 0.05, typically 0.05 is sufficient"**, and documents it as a
*run-specific* column, unlike the global ones.

**Why it is low severity.** Measured on a 373-run mouse DIA-NN 2.6 report, holding the other
five at 0.01: `PG.Q.Value` at 0.01 vs 0.05 differs by **two protein groups** (6,564 vs 6,566),
20,370 rows out of ~20 million. So the current uniform 0.01 is *conservative*, not wrong, and
tightening beyond DIA-NN's advice costs almost nothing on that dataset.

That is one measurement on one report, which is why this is recorded rather than dismissed. It
is also worth knowing that `PG.Q.Value` is the **largest single contributor** to what the six
columns exclude on that report (22,465 rows, 20,447 of them uniquely, against 15,752 for
`Global.PG.Q.Value` and 7,873 for `Global.Q.Value`) — so a change to its cutoff has more room to
matter than its "low severity" label suggests, on a dataset less forgiving than this one.

**Proposed fix.** A per-column cutoff map, defaulting `PG.Q.Value` to 0.05 and the rest to 0.01.
This is a **behaviour change**, not a bug fix: it would loosen current `maxlfq` behaviour for
every existing user. It should land in its own change, touching both paths together, after
defect #2 so there is one definition to edit. Deliberately not bundled into PR #48.

---

## 4. `--min-pr-charge 2` narrows DIA-NN's own default without saying why

**What.** `scripts/estimate_params.py:181` emits `--min-pr-charge 2`, tagged `UNIV`. DIA-NN's
own default is 1–4, verified by measurement on the 2.6.0 binary with a 60-protein FASTA:

| config | precursors generated |
|---|---|
| no charge flags (DIA-NN default) | 10,899 |
| `--min-pr-charge 1 --max-pr-charge 4` | 10,899 |
| `--min-pr-charge 2 --max-pr-charge 4` | 8,805 |

So the emitted default discards **2,094 precursors, ~19% of the predicted library**. (The DIA-NN
README documents the flags but states no default, which is why this was measured rather than
looked up.)

**This is probably the right choice** for tryptic bottom-up work — singly-charged precursors are
rarely informative — so it is listed as a decision to record, not a bug to fix. The problem is
only that it is tagged `UNIV` alongside genuinely universal settings, so a user reading the
emitted parameters cannot tell that a deliberate narrowing has been applied on their behalf.

**Proposed fix.** Change the tag text to state the rationale and the cost, e.g.
`"z=1 excluded: rarely informative for tryptic bottom-up; DIA-NN's own default is 1-4"`. No
behaviour change.

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
