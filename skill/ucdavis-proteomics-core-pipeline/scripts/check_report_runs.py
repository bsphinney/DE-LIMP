#!/usr/bin/env python3
"""
check_report_runs.py  --  Fail unless a DIA-NN report holds every input run.

DIA-NN exits 0 when it cannot read a run, and simply leaves that run out of the cross-run
report. Measured on DIA-NN 2.7.0 (HIVE, 2026-09-16): a single-shot search of one readable
Exploris .raw, one truncated .raw and one path that does not exist logged
"ERROR: DIA-NN tried but failed to load the following files: ...", exited 0, and wrote a
report.parquet whose Run column held 1 of the 3 runs. SLURM records COMPLETED, the watcher
reports success, and the DE step quietly analyses fewer samples than were acquired. The
5-step chain catches this at step 5 by counting .quant files; the single-shot search has no
per-file .quant directory to count, so it checks the report itself.

A run counts as present when the report has rows for it. The report is read two ways:

  1. report.parquet's Run column (pyarrow, else polars) -- the file the DE step consumes,
     so this is the authoritative answer.
  2. DIA-NN's own <report>.stats.tsv, with no third-party packages at all, for a job whose
     python has no parquet reader (the stdlib-only CI; HIVE's /usr/bin/python3 has none
     unless the user pip-installed one into ~/.local). DIA-NN lists every run it ATTEMPTED
     there -- in the measurement above the truncated and the nonexistent file both had rows,
     with 0 in every column -- so a row only counts when Precursors.Identified > 0. The
     .pg_matrix.tsv header is no help either: it named all 3 runs.

If neither can be read (no parquet reader and the cfg passed --no-stats) the check FAILS:
an unverified run count reported as success is exactly the failure this exists to stop.

Usage:
  check_report_runs.py --report <out>/report.parquet --files-list <file with one input per line>
Exit 0 when every input is present; 1 otherwise, with the missing runs named.
"""
import argparse
import csv
import os
import sys

# DIA-NN's main-report Run column is "raw file name without file path" (DIA-NN README,
# Main output reference), and without its extension.
_EXTS = (".raw", ".d", ".mzml", ".wiff", ".wiff2", ".dia")


def run_name(path):
    """The name DIA-NN reports a run under: the file name, no folder, no extension."""
    base = os.path.basename(str(path).replace("\\", "/").rstrip("/"))
    stem, ext = os.path.splitext(base)
    return stem if ext.lower() in _EXTS else base


def runs_from_parquet(report):
    """Distinct Run values in report.parquet, or None when no parquet reader is installed.

    A report that exists but cannot be read raises: that is a broken result, not a reason
    to fall back to a weaker check."""
    try:
        import pyarrow.parquet as pq
        col = pq.read_table(report, columns=["Run"]).column("Run")
        return {str(r) for r in col.unique().to_pylist() if r is not None}
    except ImportError:
        pass
    try:
        import polars as pl
        return {str(r) for r in pl.read_parquet(report, columns=["Run"])["Run"].unique()
                .to_list() if r is not None}
    except ImportError:
        return None


def stats_path(report):
    """DIA-NN writes <report stem>.stats.tsv next to <report stem>.parquet."""
    stem = report[:-len(".parquet")] if report.endswith(".parquet") else report
    return stem + ".stats.tsv"


def runs_from_stats(path):
    """Run names DIA-NN's stats file shows with at least one identified precursor, or None
    when there is no stats file to read."""
    if not os.path.exists(path):
        return None
    found = set()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            try:
                n = float(row.get("Precursors.Identified") or 0)
            except ValueError:
                n = 0
            if n > 0 and row.get("File.Name"):
                found.add(run_name(row["File.Name"]))
    return found


def verify(report, files, parquet_runs=runs_from_parquet):
    """(ok, message). `parquet_runs` is injectable so the stdlib route can be tested."""
    n = len(files)
    expected = [run_name(f) for f in files]
    dupes = sorted({e for e in expected if expected.count(e) > 1})
    if dupes:
        return False, (f"FAILED: {n} inputs but only {len(set(expected))} distinct run names -- "
                       f"{', '.join(dupes)} appear more than once. DIA-NN names a run by its file "
                       "name without the folder, so these would be merged into one Run in the "
                       "report. Rename or search them separately.")
    runs = parquet_runs(report)
    source = os.path.basename(report)
    if runs is None:
        sp = stats_path(report)
        runs = runs_from_stats(sp)
        source = f"{os.path.basename(sp)} (runs with identified precursors; no parquet reader)"
        if runs is None:
            return False, (f"FAILED: cannot confirm the report holds all {n} runs: this python "
                           f"({sys.executable}) has neither pyarrow nor polars to read {report}, "
                           f"and DIA-NN wrote no {sp} (was --no-stats in the cfg?). "
                           "Install pyarrow, or check the Run column by hand before using it.")
    # The COUNT decides; names only say which. Every Run in the report comes from one of the
    # inputs, so n distinct Runs from n distinctly named inputs means all are there even if
    # DIA-NN spells some format's name differently from run_name() -- a naming guess must not
    # fail a complete search.
    if len(runs) != n:
        missing = [e for e in expected if e not in runs]
        which = f" Missing: {', '.join(missing)}." if missing else ""
        return False, (f"FAILED: the report holds {len(runs)} of {n} runs ({source}).{which} "
                       "DIA-NN exits 0 when it cannot read a run -- its log names them on the "
                       "line 'ERROR: DIA-NN tried but failed to load the following files'.")
    return True, f"OK: report holds all {n} runs ({source})"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True)
    ap.add_argument("--files-list", required=True, help="one input path per line")
    a = ap.parse_args()
    with open(a.files_list) as fh:
        files = [ln.strip() for ln in fh if ln.strip()]
    try:
        ok, msg = verify(a.report, files)
    except Exception as e:                        # an unreadable report is a failed search
        ok, msg = False, f"FAILED: could not read {a.report}: {type(e).__name__}: {e}"
    print(msg, file=sys.stdout if ok else sys.stderr)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
