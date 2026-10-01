#!/usr/bin/env python3
"""
pass_comparison.py  --  Did the 5-step chain's final pass keep what its first pass found?

The chain's step 3 is a cross-run analysis of the FIRST pass (each run searched once against the
predicted library): step3_assembly.parquet, with its step3_assembly.stats.tsv. Step 5 is the
report of the second pass against the experiment's empirical library: report.parquet. DIA-NN's
README, on its own two-pass MBR: "analyse first with MBR. This will produce both the final MBR
report as well as a report based on the first MBR pass ... Compare the identification numbers and
data completeness in the two reports. If it appears that the performance is higher for the first
pass, just uncheck MBR". Nothing compared them, and the second pass can LOSE identifications.
SET28 (6 Exploris 480 DDA hair runs, DIA-NN 2.7.0, HIVE 2026-09-29), precursors per run, first
pass -> final: 103 -> 10, 566 -> 529, 126 -> 0, 138 -> 7, 269 -> 183, 165 -> 1. The 28-run SET28
chain was lower in 18 of 28 runs, while FragPipe found 400-570 spectra in those same files.

So step 5 runs this after its report is built, for every chain (DIA included: the mechanism is
the same): a per-run table of precursors and protein groups from both passes (DIA-NN's own
<report>.stats.tsv, no parquet reader needed), a WARNING for every run whose precursors fell by
more than half, the table in <out>/pass_comparison.json and in search_provenance.json
(`pass_comparison`), where audit_results.py puts it in AUDIT.md. The final report is not wrong
and the job does not fail over this: the numbers are a decision for the user.

THE STATS COUNTS ARE BEFORE PROTEIN-GROUP FDR, so a flag alone does not say the first pass would
give the DE step more data -- and what it gives depends on the cutoffs. SET28, rows per run
passing the q-columns (HIVE, 2026-09-29):
    at 1% on every column             the four flagged runs 0 in BOTH reports; N5 396 vs 147
    run_de.R's filter (PG.Q.Value 5%)  flagged runs 109-160 in the first pass vs 0-9 in the
                                       final; N5 398 vs 147
The difference is rows whose protein group has a run-specific q between 1% and 5%, which
run_de.R keeps (diann_q_columns.COLUMN_CUTOFFS). So, where pyarrow
is available, each report's rows per run that pass run_de.R's identification filter are counted
too -- the same filter, from diann_q_columns.py (columns, per-column cutoffs, run_de.R's default
--q-cutoff), applied as limpa applies it -- and switching to the first-pass report is recommended
only for runs where it has MORE of those rows. Without pyarrow the table says so, and the advice
is to compare the post-filter rows before switching.

Usage:
  pass_comparison.py --first-report <out>/step3_assembly.parquet --final-report <out>/report.parquet
                     [--out <out>/pass_comparison.json] [--provenance <out>/search_provenance.json]
Exit 0 whenever the table was made (flagged runs or not); 1 when a stats file cannot be read.
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_report_runs import _num, stats_path, stats_rows  # noqa: E402  one stats reader
# run_de.R's identification filter: the ONE definition (DE-LIMP rule 3), never a copy here
from diann_q_columns import DEFAULT_Q_CUTOFF, cutoff_for, fdr_columns  # noqa: E402

# A run is flagged when the final pass kept less than this fraction of its first-pass precursors.
DROP_FLAG = 0.5
BASIS = ("Precursors.Identified and Proteins.Identified per run, from DIA-NN's <report>.stats.tsv "
         "of each pass -- counts BEFORE protein-group FDR")
BEFORE_PG_FDR = ("the stats counts are before protein-group FDR: compare the rows that pass "
                 "run_de.R's q-value filter before switching")


def rows_after_q(report, q_cutoff=DEFAULT_Q_CUTOFF):
    """({run: rows passing run_de.R's q-value filter}, columns used, None), or (None, [], why).

    The filter is run_de.R's: the q-columns diann_q_columns.fdr_columns() finds in the report,
    each at diann_q_columns.cutoff_for() -- applied as limpa's EListFromLongFormatFile does,
    `which(Report[[col]] > cutoff)` dropped, so a row is kept unless a value EXCEEDS its cutoff
    (a missing value is kept)."""
    try:
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ImportError:
        return None, [], "pyarrow is not installed in this python"
    if not os.path.isfile(report):
        return None, [], f"{report} not found"
    cols = fdr_columns(pq.read_schema(report).names)
    if not cols:
        return None, [], f"{os.path.basename(report)} has no q-value column run_de.R filters on"
    t = pq.read_table(report, columns=["Run"] + cols)
    drop = None
    for c in cols:
        over = pc.fill_null(pc.greater(t[c], cutoff_for(c, q_cutoff)), False)
        drop = over if drop is None else pc.or_(drop, over)
    vc = pc.value_counts(t.filter(pc.invert(drop))["Run"])
    return ({str(v): int(n) for v, n in zip(vc.field("values").to_pylist(),
                                              vc.field("counts").to_pylist())}, cols, None)


def compare(first_report, final_report, drop=DROP_FLAG, q_cutoff=DEFAULT_Q_CUTOFF):
    """The per-run comparison record. Raises OSError when either stats file is missing."""
    first_stats, final_stats = stats_path(first_report), stats_path(final_report)
    first, final = stats_rows(first_stats), stats_rows(final_stats)
    for rows, path in ((first, first_stats), (final, final_stats)):
        if rows is None:
            raise OSError(f"no DIA-NN stats file at {path}")
    q_first, cols_first, why_first = rows_after_q(first_report, q_cutoff)
    q_final, cols_final, why_final = rows_after_q(final_report, q_cutoff)
    have_q = q_first is not None and q_final is not None
    runs = []
    for run in sorted(set(first) | set(final)):
        f, g = first.get(run) or {}, final.get(run) or {}
        f_pr, g_pr = _num(f, "Precursors.Identified"), _num(g, "Precursors.Identified")
        change = (g_pr - f_pr) / f_pr if f_pr > 0 else None
        fq, gq = ((q_first.get(run, 0), q_final.get(run, 0)) if have_q else (None, None))
        runs.append({"run": run,
                     "first": {"precursors": int(f_pr), "proteins": int(_num(f, "Proteins.Identified")),
                               "rows_after_q": fq},
                     "final": {"precursors": int(g_pr), "proteins": int(_num(g, "Proteins.Identified")),
                               "rows_after_q": gq},
                     "precursor_change": None if change is None else round(change, 3),
                     "flagged": change is not None and g_pr < (1 - drop) * f_pr,
                     # the only question that decides a switch: does the first pass give the DE
                     # step MORE rows for this run? None when that could not be counted
                     "first_pass_has_more_after_q": (fq > gq) if have_q else None})
    flagged = [r["run"] for r in runs if r["flagged"]]
    better = [r["run"] for r in runs if r["first_pass_has_more_after_q"]]
    post = {"available": have_q, "q_cutoff": q_cutoff,
            "columns": {"first": cols_first, "final": cols_final},
            "cutoffs": {c: cutoff_for(c, q_cutoff) for c in (cols_first or cols_final)},
            "basis": ("rows per run passing run_de.R's q-value filter (diann_q_columns: the columns "
                      "present, cutoff_for() per column, a row dropped when a value exceeds its "
                      "cutoff, as limpa does)"),
            "why_unavailable": None if have_q else "; ".join(
                x for x in (why_first and f"first pass: {why_first}",
                            why_final and f"final pass: {why_final}") if x)}
    return {"first_pass_report": first_report, "first_pass_stats": first_stats,
            "final_report": final_report, "final_stats": final_stats,
            "basis": BASIS, "threshold": f"final precursors below {1 - drop:.0%} of the first pass",
            "n_runs": len(runs), "n_flagged": len(flagged), "flagged": flagged, "runs": runs,
            "post_filter": post, "switch_recommended_for": better if have_q else None,
            "advice": advice(flagged, better, have_q, len(runs), first_report)}


def advice(flagged, better, have_q, n_runs, first_report):
    """What to tell the user. A switch is recommended only for runs where the first pass has
    more rows AFTER run_de.R's q-value filter; the flag itself is a pre-FDR count."""
    name = os.path.basename(first_report)
    what = f"the first-pass report ({name}: each run searched once, no match-between-runs)"
    if not flagged and not better:
        return "no run lost more than half its first-pass precursors"
    head = (f"{len(flagged)} of {n_runs} runs lost more than half their first-pass precursors in "
            "the final pass (" + BEFORE_PG_FDR + "). " if flagged else "")
    if not have_q:
        return (head + f"The post-filter rows could not be counted here; count them (pyarrow) "
                f"before offering {what} -- DIA-NN's README: use the first pass when it "
                "performs better")
    if better:
        return (head + f"After run_de.R's q-value filter the first pass has more rows for "
                f"{len(better)} run(s): {', '.join(better)}. Tell the user and offer {what} for "
                "the DE -- DIA-NN's README: use the first pass when it performs better")
    return (head + "After run_de.R's q-value filter the first pass has no more rows than the final "
            "pass for any run, so switching to it would not recover them. Tell the user; keep the "
            "final report")


def table(rec):
    """The comparison as a markdown table, one row per run. The first four count columns are
    DIA-NN's stats, before protein-group FDR; the last two are rows after run_de.R's filter."""
    lines = ["| Run | First-pass precursors | Final precursors | Change | First-pass protein "
             "groups | Final protein groups | First-pass rows after q | Final rows after q |",
             "|---|---|---|---|---|---|---|---|"]
    for r in rec["runs"]:
        ch = "n/a" if r["precursor_change"] is None else f"{r['precursor_change']:+.0%}"
        fq, gq = r["first"].get("rows_after_q"), r["final"].get("rows_after_q")
        more = " **first pass has more**" if r.get("first_pass_has_more_after_q") else ""
        lines.append(f"| {r['run']} | {r['first']['precursors']} | {r['final']['precursors']} | "
                     f"{ch}{' **FLAGGED**' if r['flagged'] else ''} | {r['first']['proteins']} | "
                     f"{r['final']['proteins']} | {'n/a' if fq is None else fq}{more} | "
                     f"{'n/a' if gq is None else gq} |")
    lines.append(f"Precursor and protein-group counts: {BEFORE_PG_FDR}. Rows after q: "
                 + ("run_de.R's q-value filter on " + ", ".join(
                     f"{c} <= {v:g}" for c, v in rec['post_filter']['cutoffs'].items())
                    if rec["post_filter"]["available"] else
                    f"not counted ({rec['post_filter']['why_unavailable']})") + ".")
    return "\n".join(lines)


def warnings(rec):
    """One WARNING line per flagged run, for the job log."""
    return [f"WARNING: {r['run']}: {r['first']['precursors']} precursors in the first pass, "
            f"{r['final']['precursors']} in the final pass ({r['precursor_change']:+.0%}; "
            "before protein-group FDR"
            + ("" if r["first"].get("rows_after_q") is None else
               f"; rows after run_de.R's q filter {r['first']['rows_after_q']} vs "
               f"{r['final']['rows_after_q']}") + ")"
            for r in rec["runs"] if r["flagged"]]


def record_in_provenance(path, rec):
    """Add the record to search_provenance.json as `pass_comparison`, in place (a .tmp then a
    rename, so a reader never sees half a file). Returns False when there is no such file."""
    if not os.path.isfile(path):
        return False
    with open(path) as fh:
        prov = json.load(fh)
    prov["pass_comparison"] = rec
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(prov, fh, indent=2)
    os.replace(tmp, path)
    return True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--first-report", required=True, help="step 3's first-pass report")
    ap.add_argument("--final-report", required=True, help="step 5's report")
    ap.add_argument("--out", help="write the record here as JSON")
    ap.add_argument("--provenance", help="search_provenance.json to record it in")
    a = ap.parse_args(argv)
    try:
        rec = compare(a.first_report, a.final_report)
    except (OSError, ValueError) as e:
        sys.exit(f"[pass_comparison] could not compare the passes: {e}")
    print("First pass (step 3) vs final pass (step 5), per run:")
    print(table(rec))
    for w in warnings(rec):
        print(w, file=sys.stderr)
    print(("WARNING: " if rec["flagged"] else "") + rec["advice"],
          file=sys.stderr if rec["flagged"] else sys.stdout)
    if a.out:
        with open(a.out, "w") as fh:
            json.dump(rec, fh, indent=2)
    if a.provenance and not record_in_provenance(a.provenance, rec):
        print(f"[pass_comparison] {a.provenance} not found: recorded in {a.out or 'the log'} "
              "only", file=sys.stderr)


if __name__ == "__main__":
    main()
