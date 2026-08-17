#!/usr/bin/env python3
"""
diann_q_columns.py -- ONE definition of DIA-NN's q-value columns for the skill.

Fixes SKILL_OPEN_DEFECTS #2: the column set was hand-written in five files. That
is how the dpc path in run_de.R drifted away from build_maxlfq.R and applied a
different identification FDR to the same report for months (PR #48), and the same
class of defect was later found in the DE-LIMP app itself.

There are TWO distinct concepts here, and conflating them is a bug, not a
simplification:

  FDR_REQUIRED + FDR_OPTIONAL   the FILTER SET. Every one of these is applied,
                               ANDed together, at the q-cutoff. Used by run_de.R
                               (dpc) and build_maxlfq.R.

  PROTEIN_Q_PREFERENCE         a PREFERENCE ORDER for picking ONE protein-level
                               q-column to read for reporting or cross-tool
                               comparison. The first available wins; the rest are
                               ignored. Used by compare_searches.py and
                               run_search.py.

Applying the preference order as a filter would over-filter; filtering on only
the first available column would under-filter. They are not interchangeable.

COLUMN SEMANTICS -- from the DIA-NN 2.6 README, "Main output reference"
(github.com/vdemichev/DiaNN@master; the 1.8 tag serves a much shorter, different
README, so re-check against master):

  Q.Value            L1242  run-specific precursor q-value
  Global.Q.Value     L1244  global (experiment-wide) precursor q-value
  Lib.Q.Value        L1245  library-entry q-value; 'global' when DIA-NN built the
                            library, so NOT reliably run-level
  PG.Q.Value         L1255  RUN-SPECIFIC protein-group q-value -- despite sitting
                            beside the Global.* pair. Do not call it
                            experiment-wide.
  Global.PG.Q.Value  L1259  global protein-group q-value
  Lib.PG.Q.Value     L1260  library-entry protein-group q-value

Run-level q-values alone do not control FDR across an experiment: a union of
run-level-passing IDs over N runs sits above the nominal cutoff and grows with N,
and the inflated protein list also inflates the family size m that BH corrects
over. Hence the Global.* pair is applied whenever present.

Measured caveat on Lib.*: their values are not dependable across reports. On one
2-run HeLa report both were identically 0.0 in 100% of rows (no filtering power
at all); on a 399-run mouse report they filtered 0.56% / 0.19%. So neither "the
Lib columns already give global control" nor "there is no global control" is true
in general -- which is exactly why Global.* is named explicitly rather than
assumed to be covered.

The R mirror is diann_q_columns.R. The two are asserted identical by
tests/test_q_columns.py, so drift is a CI failure rather than a silent
divergence.
"""

# Required by both quant paths; a report lacking these is not a DIA-NN report.
FDR_REQUIRED = ["Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"]

# Applied additionally WHEN PRESENT. Optional only because older reports predate
# them -- not because they are discretionary.
FDR_OPTIONAL = ["PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value"]

# Preference order for reading a SINGLE protein-level q-value. Most global and
# most protein-specific first, degrading toward the run-level precursor q-value.
PROTEIN_Q_PREFERENCE = [
    "Global.PG.Q.Value",
    "Global.Q.Value",
    "Lib.PG.Q.Value",
    "PG.Q.Value",
    "Q.Value",
]


# Per-column cutoffs. DIA-NN's README recommends PG.Q.Value "at 0.01 to 0.05,
# typically 0.05 is sufficient" and documents it as RUN-SPECIFIC, unlike the
# Global.* pair -- so applying one uniform cutoff to all six was tighter than its
# own author advises. Measured on a 373-run mouse report, holding the other five
# at 0.01, PG.Q.Value at 0.01 vs 0.05 differs by two protein groups (6,564 vs
# 6,566). Small there, but PG.Q.Value is the LARGEST single contributor to what
# the six columns exclude on that report (22,465 rows, 20,447 uniquely), so it has
# more room to matter on a less forgiving dataset.
#
# Anything not named here uses the run's --q-cutoff.
COLUMN_CUTOFFS = {"PG.Q.Value": 0.05}


def cutoff_for(column, q_cutoff, uniform=False):
    """The cutoff to apply to `column` given the run's --q-cutoff.

    A column listed in COLUMN_CUTOFFS uses its own value; everything else uses
    q_cutoff. `uniform=True` disables the per-column map entirely, restoring the
    pre-2026-08 behaviour of one cutoff for all six.

    Note this deliberately does NOT try to be clever about a caller who tightened
    --q-cutoff. The pipeline default (0.01) is already stricter than PG.Q.Value's
    recommended 0.05, so any "never loosen what the caller asked for" rule would
    prevent the recommendation from ever applying. --q-cutoff governs the other
    five columns; pass uniform=True to make it govern all six.
    """
    if uniform:
        return q_cutoff
    return COLUMN_CUTOFFS.get(column, q_cutoff)


def fdr_columns(available=None):
    """The columns to FILTER on, restricted to those the report has.

    `available` is any container of column names; None means "assume all present".
    Never returns a column the report lacks -- limpa::readDIANN() errors on an
    unknown name, and arrow silently returns ZERO ROWS when filtering a column
    that was not selected.
    """
    if available is None:
        return list(FDR_REQUIRED) + list(FDR_OPTIONAL)
    have = set(available)
    return [c for c in (FDR_REQUIRED + FDR_OPTIONAL) if c in have]


def protein_q_column(available):
    """The single best protein-level q-column present, or None."""
    have = set(available)
    for c in PROTEIN_Q_PREFERENCE:
        if c in have:
            return c
    return None


if __name__ == "__main__":  # tiny self-report, handy when debugging a report
    import json
    import sys
    cols = sys.argv[1:] or None
    print(json.dumps({
        "fdr_columns": fdr_columns(cols),
        "protein_q_column": protein_q_column(cols) if cols else None,
    }, indent=2))
