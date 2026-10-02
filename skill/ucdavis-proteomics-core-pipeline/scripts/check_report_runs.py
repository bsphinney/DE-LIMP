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
    """DIA-NN writes <report stem>.stats.tsv next to <report stem>.parquet (DIA-NN 1.x:
    next to <report stem>.tsv)."""
    stem = next((report[:-len(e)] for e in (".parquet", ".tsv") if report.endswith(e)), report)
    return stem + ".stats.tsv"


def _num(row, col):
    try:
        return float(row.get(col) or 0)
    except ValueError:
        return 0.0


def stats_rows(path):
    """{run name: row} for every run in DIA-NN's stats file, or None when there is none."""
    if not os.path.exists(path):
        return None
    with open(path, newline="") as fh:
        return {run_name(r["File.Name"]): r for r in csv.DictReader(fh, delimiter="\t")
                if r.get("File.Name")}


def runs_from_stats(path):
    """Run names DIA-NN's stats file shows with at least one identified precursor, or None
    when there is no stats file to read."""
    rows = stats_rows(path)
    if rows is None:
        return None
    return {name for name, r in rows.items() if _num(r, "Precursors.Identified") > 0}


def distinct_inputs(files):
    """(each input once, [repeated paths]) -- THE rule for one file listed more than once, which
    every route that takes a file list follows (run_search.py for every engine, diann_parallel.py,
    radiant_parallel.py, ht_manifest.py, and verify() below): the file is searched ONCE, because
    an engine cannot take it twice, and the repeat is FLAGGED (repeats_text), never refused and
    never dropped silently. "The same file" is judged on the RESOLVED path, so a symlink and its
    target, or `a/./x.d` and `a/x.d`, are one. The first spelling is kept, in first-seen order. A
    repeat is {"path", "times", "also_listed_as": [other spellings]}.

    Why: STAN's HT manifest listed one run 4 times and another twice (100 lines for 96 files) with
    every gate PASS (a staff plate, 2026-10-01); searched as given, those runs would count 2-4
    times in the experiment and the empirical library. Brett (2026-10-01): repeats are okay, but
    they are flagged."""
    order, seen = [], {}
    for f in files:
        key = os.path.realpath(str(f).rstrip("/"))
        if key in seen:
            seen[key].append(f)
        else:
            seen[key] = [f]
            order.append(f)
    repeats = [{"path": v[0], "times": len(v),
                "also_listed_as": sorted({x for x in v[1:] if x != v[0]})}
               for v in seen.values() if len(v) > 1]
    return order, repeats


def repeated_names(files):
    """[{"run_name", "paths"}]: DIFFERENT files (distinct_inputs) that share a run name -- a rerun
    in another folder, one name on two plates. Kept in every list and flagged (repeats_text). An
    engine cannot search them as given: DIA-NN names a run by its file name without the folder,
    so /plate1/s1.raw and /plate2/s1.raw become ONE Run in the report (and, in the 5-step chain,
    two array tasks writing one .quant), and Sage converts both to one mzML name -- so
    run_search.py stops before any engine with this list in its message."""
    by = {}
    for f in distinct_inputs(files)[0]:
        by.setdefault(run_name(f), []).append(f)
    return [{"run_name": n, "paths": v} for n, v in by.items() if len(v) > 1]


def repeats_text(repeated_paths=(), repeated_names_=()):
    """THE wording of the flag, one sentence per kind -- stderr, the HT manifest, and the analysis
    report all print these."""
    out = []
    if repeated_paths:
        out.append(f"{len(repeated_paths)} input file(s) were listed more than once; each is "
                   f"searched once: " + "; ".join(
                       f"{r['path']} (listed {r['times']} times"
                       + (f", also as {', '.join(r['also_listed_as'])}" if r.get("also_listed_as")
                          else "") + ")" for r in repeated_paths))
    if repeated_names_:
        out.append(f"{len(repeated_names_)} run name(s) are shared by different files, all kept: "
                   + "; ".join(f"{r['run_name']}: {', '.join(r['paths'])}" for r in repeated_names_))
    return out


def repeats_note(repeated_paths, prog, repeated_names_=()):
    """repeats_text() as stderr lines."""
    return "".join(f"[{prog}] FLAG: {t}\n" for t in repeats_text(repeated_paths, repeated_names_))


def names_stop(repeated_names_, prog):
    """THE message an engine route stops with when different inputs share a run name: what is
    shared, that the list keeps and flags them, and why the search cannot take them as given."""
    return (f"[{prog}] STOPPED before searching: different input files share a run name -- "
            + "; ".join(f"{r['run_name']}: {', '.join(r['paths'])}" for r in repeated_names_)
            + ". They are kept in the file list and flagged (repeated_names), but a search cannot "
              "take them as given: DIA-NN names a run by its file name without the folder, so they "
              "would become ONE run in the report (two samples merged), and Sage converts both to "
              "one mzML name; every engine route holds to the same rule. Ways out: if they are "
              "re-injections of one sample, keep one (core_submission.py locate --reinjections "
              "latest, or leave the other out of the list); otherwise give one a different name "
              "(a link under another name works), or search them separately. Searching both under "
              "unique run names is a 2.11 item. Nothing was written or submitted.")


# DIA-NN's per-run Normalisation.Instability, a column of <report>.stats.tsv (DIA-NN 2.7.0 writes
# it, between Median.Mass.Acc.MS2.Corrected and Median.RT.Prediction.Acc). Its README documents
# neither what the value means nor a cutoff, so the threshold is this skill's, set between two
# measured cohorts (a staff report, DIA-NN 2.7.0): an E. coli cohort at 0.04-0.07, and a
# timsTOF tissue cohort at 0.67-1.00 in every run, whose Precursor.Normalised /
# Precursor.Quantity ratio then spread ~8-fold between runs and varied with RT inside each run.
# The ONE definition: audit_results.py warns on it and record_run.py records it.
NORM_INSTABILITY_COLUMN = "Normalisation.Instability"
NORM_INSTABILITY_WARN = 0.3
NORM_INSTABILITY_BASIS = (
    f"DIA-NN documents no cutoff for {NORM_INSTABILITY_COLUMN}; {NORM_INSTABILITY_WARN:g} is this "
    "skill's, between a cohort that measured 0.04-0.07 and one that measured 0.67-1.00 in every "
    "run (DIA-NN 2.7.0)")


def normalisation_instability(report, runs=None):
    """DIA-NN's Normalisation.Instability per run, from the stats file beside `report`.

    -> None when there is no stats file (another engine's report, or --no-stats), else
    {"stats_file", "column" (False: this DIA-NN wrote no such column), "per_run" {run: value},
    "median", "max", "flagged" {run: value above NORM_INSTABILITY_WARN}, "threshold", "basis"}.
    A run with no identified precursor is left out: DIA-NN writes 0 for a run it could not load or
    found nothing in, which says nothing about its normalisation. `runs` (run names, as
    run_name() gives them) restricts the answer to the runs an analysis used."""
    path = stats_path(report)
    rows = stats_rows(path)
    if rows is None:
        return None
    rec = {"stats_file": path, "column": False, "per_run": {}, "median": None, "max": None,
           "flagged": {}, "threshold": NORM_INSTABILITY_WARN, "basis": NORM_INSTABILITY_BASIS}
    want = {run_name(r) for r in runs} if runs else None
    for name, r in rows.items():
        if NORM_INSTABILITY_COLUMN not in r:
            return rec
        rec["column"] = True
        if (want is not None and name not in want) or _num(r, "Precursors.Identified") <= 0:
            continue
        rec["per_run"][name] = _num(r, NORM_INSTABILITY_COLUMN)
    vals = sorted(rec["per_run"].values())
    if vals:
        mid = len(vals) // 2
        rec["median"] = vals[mid] if len(vals) % 2 else (vals[mid - 1] + vals[mid]) / 2
        rec["max"] = vals[-1]
    rec["flagged"] = {k: v for k, v in rec["per_run"].items() if v > NORM_INSTABILITY_WARN}
    return rec


def duplicate_run_names(files):
    """The run names repeated_names() finds, sorted -- what a check that only needs the names
    reads (the step-5 report check)."""
    return sorted(r["run_name"] for r in repeated_names(files))


def sharing_run_name(files, names):
    """'s1: /a/s1.raw, /b/s1.raw; ...' -- the distinct inputs behind each shared run name."""
    return "; ".join(f"{r['run_name']}: {', '.join(r['paths'])}" for r in repeated_names(files)
                     if r["run_name"] in names)


def why_missing(name, stats, stats_file):
    """One line saying what the evidence shows about a run that is not in the report.

    A run is absent both when DIA-NN could not load it and when it loaded the run but nothing
    passed the q-value filter (a blank or a failed injection), and the log line naming load
    failures exists only for the first. The stats file separates them only one way round: a
    run with MS1/MS2 signal was read. Measured on DIA-NN 2.7.0 (HIVE srun 23512769): a wash
    injection searched beside a real run had MS1.Signal 3.48e11 and MS2.Signal 5.39e9 with 1
    precursor identified, where the unreadable runs of the earlier measurement had 0. All zeros
    is NOT proof of a load failure, though: two readable runs searched against an empty library
    also had 0 in every column, so that case keeps both explanations."""
    r = (stats or {}).get(name)
    load = ("DIA-NN could not load it (its log names such files after 'ERROR: DIA-NN tried "
            "but failed to load the following files')")
    if r is not None and max(_num(r, "MS1.Signal"), _num(r, "MS2.Signal")) > 0:
        return (f"  {name}: DIA-NN read it ({stats_file}: MS1.Signal {_num(r, 'MS1.Signal'):.3g},"
                f" MS2.Signal {_num(r, 'MS2.Signal'):.3g}) but identified nothing that passed "
                "the report's q-value filter -- a blank, wash or failed injection? If it is "
                "not a sample, leave it out of the inputs.")
    seen = f"all zeros in {stats_file}" if r is not None else "not in the report"
    return f"  {name}: {seen} -- {load}, or it identified nothing in it."


def verify(report, files, parquet_runs=runs_from_parquet):
    """(ok, message). `parquet_runs` is injectable so the stdlib route can be tested."""
    files = distinct_inputs(files)[0]           # a file listed twice was searched once
    n = len(files)
    expected = [run_name(f) for f in files]
    dupes = duplicate_run_names(files)
    if dupes:
        return False, (f"FAILED: {n} inputs but only {len(set(expected))} distinct run names -- "
                       f"different files share a name ({sharing_run_name(files, dupes)}). DIA-NN "
                       "names a run by its file name without the folder, so these would be merged "
                       "into one Run in the report. Rename or search them separately.")
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
        sp = stats_path(report)
        stats = stats_rows(sp)
        lines = [f"FAILED: the report holds {len(runs)} of {n} runs ({source}).{which}"]
        lines += [why_missing(m, stats, os.path.basename(sp)) for m in missing]
        lines.append("DIA-NN exits 0 either way, so this job fails rather than hand the DE "
                     "step fewer samples than were acquired.")
        return False, "\n".join(lines)
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
