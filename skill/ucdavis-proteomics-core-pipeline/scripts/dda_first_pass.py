#!/usr/bin/env python3
"""
dda_first_pass.py  --  A DIA-NN --dda search whose SECOND pass crashed falls back to the complete
first-pass report DIA-NN had already written, and says so.

With match-between-runs (--reanalyse) DIA-NN searches every run once, writes that first pass as
<report stem>-first-pass.parquet, and then searches again against the empirical library it built.
A Core submission (2026-10-29, DIA-NN 2.7.0 --dda, timsTOF .d): the second pass crashed, and the search
folder held report-first-pass.parquet (and its .manifest.txt, .pr_matrix.tsv, .stats.tsv) beside
report-lib.parquet -- and no report.parquet, report.pg_matrix.tsv or report.stats.tsv. The job
failed, the skill had nothing to analyse, and the session left the skill's path for a hand-made
report. The first pass is a complete search of every run, only without match-between-runs.

So a --dda search job (run_search.run_diann; job_lines() below) runs DIA-NN without `set -e`, and
when the report is missing or empty AND the first pass is not, `record` makes the first pass the
report:

  1. copies it to <report>.part, then renames that onto the report -- a copy, not a link, which
     would break in a zip or a delivery -- and its .stats.tsv and .pr_matrix.tsv onto the
     report's (COPIED_SUFFIXES): the job's report guard reads the stats when its interpreter has
     no parquet reader, and so do the audit's normalisation check and the Slack summary; a
     delivery takes the precursor matrix under the report's name (core_submission.SEARCH_FILES).
     Both describe exactly the report beside them, and the caveat travels with them in
     search_provenance.json;
  2. writes <out>/dda_first_pass.json: {report_used: "first_pass", exit_code (DIA-NN's status:
     the evidence, kept here and in AUDIT.md, never in the Methods), first_pass_sha256, at,
     caveat, ...};
  3. sets search_provenance.json `report_used: "first_pass"` and `dda_first_pass: {...}`. A job
     whose search_provenance.json is missing fails instead: the Methods, AUDIT.md and the
     watcher read only that file, so a fallback recorded nowhere else must not be one.

The first pass's .manifest.txt keeps its own name. Every first-pass file of an earlier search
goes before the job runs (stale_files()).

The rest of the job is unchanged: the report guard (check_report_runs.py) then checks the copied
report as it checks any other. caution() is the ONE wording of the caveat: the Methods
(make_methods.py), AUDIT.md (audit_results.py), the watcher (watch_run.sh), the Core run
registry (record_run.py), the job-end Slack post (notify_slack.py) and checkpoint.py status all
print it (caution_for() for a folder).

Never a fallback when the final report exists: the second pass wrote it, and its exit status
decides as before. A crash with no first pass fails the job with DIA-NN's own status. A DIA
search's job has none of this.

The record describes THIS run of the job. A --dda job therefore starts, right after it clears
the previous report, with forget_lines(): `report_used` and `dda_first_pass` leave
search_provenance.json (every other key is kept). A job resubmitted after a fallback -- what
the two-job split is for, requeueing against the same library -- whose second pass then
succeeds reads as the ordinary second-pass search it is, in the Methods, AUDIT.md, the watcher
and FRAN's manifest (L8R-F1).

Engine protein quantities (Ruling 53): DIA-NN 2.7.0's --dda first-pass report has PG.MaxLFQ = 0
in every row (two timsTOF DDA-PASEF runs; Precursor.Normalised > 0 everywhere). pg_maxlfq_state()
is the ONE definition of "this report has engine protein quantities" (is any PG.MaxLFQ > 0?):
the record carries it as `pg_maxlfq` ("absent" / "present"); build_maxlfq.R asks it (the
`pg-maxlfq` command) and refuses --method maxlfq on an "absent" report with MAXLFQ_REFUSAL --
it used to drop every run and stop with "fewer than 2 usable runs survived the filters"; and
caution() adds NO_PG_MAXLFQ for an "absent" record, so the Methods and AUDIT.md say the protein
quantities come from the precursors.

Usage:
  dda_first_pass.py record --report <out>/report.parquet --first-pass <out>/report-first-pass.parquet
      --provenance <out>/search_provenance.json --exit-code N
  dda_first_pass.py caution --provenance search_provenance.json   # the caution, or nothing
  dda_first_pass.py forget --provenance search_provenance.json    # a new run of the job
  dda_first_pass.py pg-maxlfq --report report.parquet             # JSON: pg_maxlfq, refusal

Python 3.9, stdlib only (pyarrow, else polars, to read a parquet's PG.MaxLFQ: without either the
state is not checked, and says so).
"""
import argparse
import datetime
import hashlib
import json
import os
import shlex
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from check_report_runs import stats_path  # noqa: E402  one definition of <report>.stats.tsv
from safe_write import write_file  # noqa: E402  no write through a planted link

# DIA-NN 2.7.0 writes <report stem>-first-pass.parquet (that submission's output/search_dda/).
FIRST_PASS_SUFFIX = "-first-pass.parquet"
RECORD = "dda_first_pass.json"
REPORT_USED = "first_pass"
# the job log's line when the fallback ran (grep for it)
MARKER = "DIANN_DDA_FIRST_PASS_FALLBACK"
# No direction is claimed: there is no second-pass report to compare with, and on SET28 (Exploris
# DDA, DIA-NN 2.7.0) the second pass had FEWER precursors in 18 of 28 runs (pass_comparison.py).
# The shell status is in RECORD and AUDIT.md's evidence, not in this Methods sentence (L8R-F4).
CAVEAT = ("DIA-NN's second pass (match-between-runs) failed for this search in DIA-NN's DDA "
          "mode; these results are DIA-NN's complete first-pass report, identified without "
          "match-between-runs: no "
          "identification was transferred between runs, so a precursor is quantified only in the "
          "runs where it was identified on its own, and more values may be missing across runs.")

# Engine protein quantities: pg_maxlfq_state()'s answers, and the words for "absent" (one each).
PG_MAXLFQ_PRESENT, PG_MAXLFQ_ABSENT = "present", "absent"
NO_ENGINE_QUANT = "no engine protein quantities (PG.MaxLFQ is 0 or missing in every row)"
# The search-level fact, and what a DE does about it. The staff views (watcher, Slack, registry,
# checkpoint, AUDIT -- which runs only after a DE, and a DE on this report can only be dpc) say
# both, looking ahead (NO_PG_MAXLFQ). The Methods state the fact in the search paragraph
# (caution(prov, methods=True)) and the DPC-Quant clause only in the DE paragraph, only when the
# DE that ran was dpc (de_sentence): a no-DE Methods once said limpa estimated protein
# quantities that nothing computed (2.14.2 C2 biologist B2 = statistician H2).
NO_PG_MAXLFQ_FACT = f"This first-pass report has {NO_ENGINE_QUANT}."
DPC_FROM_PRECURSORS = ("protein quantities are estimated from its precursors (limpa DPC-Quant), "
                       "never from DIA-NN's MaxLFQ")
NO_PG_MAXLFQ = f"This first-pass report has {NO_ENGINE_QUANT}, so {DPC_FROM_PRECURSORS}."
DE_DPC_SENTENCE = ("The search's first-pass report has no DIA-NN protein (MaxLFQ) quantities, so "
                   "the protein quantities were estimated from its precursors with limpa "
                   "DPC-Quant.")
MAXLFQ_REFUSAL = (f"--method maxlfq: this report has {NO_ENGINE_QUANT}, e.g. a DIA-NN --dda "
                  f"first-pass report. MaxLFQ reads exactly those quantities, so it cannot run "
                  f"on this report. Use --method dpc, which estimates the protein quantities from "
                  f"this same report's precursors (no new search).")


# What DIA-NN writes beside a report as <report stem><suffix>: the fallback copies these from
# the first pass onto the report's names, and a --dda job clears the first pass's before it runs
COPIED_SUFFIXES = (".stats.tsv", ".pr_matrix.tsv")
FIRST_PASS_SIDE_SUFFIXES = COPIED_SUFFIXES + (".manifest.txt",)


def _stem(report):
    return report[:-len(".parquet")] if report.endswith(".parquet") else report


def first_pass_path(report):
    """The first-pass report DIA-NN writes beside `report`."""
    return _stem(report) + FIRST_PASS_SUFFIX


def companion(report, suffix):
    """<report stem><suffix>: the .stats.tsv / .pr_matrix.tsv / .manifest.txt DIA-NN writes beside
    a report; .stats.tsv is check_report_runs.stats_path(), its one definition."""
    return stats_path(report) if suffix == ".stats.tsv" else _stem(report) + suffix


def stale_files(report):
    """What a --dda job clears with the report before it runs, so nothing of an earlier search's
    fallback survives into this one: the first pass and its side files, the report's precursor
    matrix (the fallback's copy writes it), and RECORD."""
    fp = first_pass_path(report)
    return ([fp] + [companion(fp, x) for x in FIRST_PASS_SIDE_SUFFIXES]
            + [companion(report, ".pr_matrix.tsv"), record_path(report)])


def record_path(report):
    """<out>/dda_first_pass.json, beside the report."""
    return os.path.join(os.path.dirname(os.path.abspath(report)), RECORD)


def forget_lines(provenance):
    """Bash for the top of a --dda search job, after its previous report is cleared: this run's
    record starts without an earlier run's fallback (forget())."""
    q = shlex.quote
    return [f'PY={q(sys.executable)}; [ -x "$PY" ] || PY=python3',
            f'"$PY" {q(os.path.abspath(__file__))} forget --provenance {q(provenance)}']


def forget(provenance):
    """Take `report_used` and `dda_first_pass` out of search_provenance.json, keeping every other
    key. True when there was something to take out; False when the file holds neither, or does
    not exist (generation said so when it could not write one). An unreadable file raises."""
    prov = _load(provenance)
    if prov is None or not ({"report_used", "dda_first_pass"} & set(prov)):
        return False
    for k in ("report_used", "dda_first_pass"):
        prov.pop(k, None)
    write_file(provenance, json.dumps(prov, indent=2))
    return True


def job_lines(search_cmd, report, provenance):
    """Bash for a --dda search job: `search_cmd` without `set -e`, then the fallback when the
    report is missing or empty and the first pass is not; otherwise DIA-NN's non-zero status
    ends the job, as `set -e` did. Every path is shell-quoted."""
    q = shlex.quote
    fp = first_pass_path(report)
    return [
        "# DIA-NN --dda: a crashed second pass leaves the complete first-pass report; use it, and",
        "# say so (dda_first_pass.py). A final report, or no first pass: as before.",
        "set +e",
        search_cmd,
        "_dia_rc=$?",
        "set -e",
        f"if [ ! -s {q(report)} ] && [ -s {q(fp)} ]; then",
        f'  PY={q(sys.executable)}; [ -x "$PY" ] || PY=python3',
        f'  "$PY" {q(os.path.abspath(__file__))} record --report {q(report)} --first-pass {q(fp)} '
        f'--provenance {q(provenance)} --exit-code "$_dia_rc"',
        'elif [ "$_dia_rc" -ne 0 ]; then',
        '  exit "$_dia_rc"',
        "fi",
    ]


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _copy(src, dst):
    """`src` onto `dst` through <dst>.part and a rename: a reader sees the old file or the whole
    new one, and a link planted at the .part is removed, never written through."""
    part = dst + ".part"
    if os.path.lexists(part):
        os.unlink(part)
    fd = os.open(part, os.O_CREAT | os.O_EXCL | os.O_WRONLY | getattr(os, "O_NOFOLLOW", 0), 0o644)
    try:
        with os.fdopen(fd, "wb") as out, open(src, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                out.write(chunk)
        os.replace(part, dst)
    except BaseException:
        if os.path.lexists(part):
            os.unlink(part)
        raise


def nonempty(path):
    """`[ -s path ]`, as the job tests it: a file with something in it."""
    return os.path.isfile(path) and os.path.getsize(path) > 0


def _positive(x):
    try:
        return float(x) > 0           # NaN > 0 is False: NaN is missing, as DIA-NN's NA is
    except (TypeError, ValueError):
        return False


def pg_maxlfq_state(report):
    """THE definition of "this report has engine protein quantities": PG_MAXLFQ_PRESENT when any
    PG.MaxLFQ is > 0, else PG_MAXLFQ_ABSENT (0, NaN or null in every row, or no such column).
    A .parquet is read with pyarrow, else polars, one batch at a time (a quantified report
    answers on its first); anything else as a tab-separated report. -> None when there is no
    parquet reader to check with. A report that cannot be read raises."""
    col = "PG.MaxLFQ"
    if not str(report).lower().endswith(".parquet"):
        import csv
        with open(report, newline="") as fh:
            rows = csv.DictReader(fh, delimiter="\t")
            if col not in (rows.fieldnames or ()):
                return PG_MAXLFQ_ABSENT
            return (PG_MAXLFQ_PRESENT if any(_positive(r.get(col)) for r in rows)
                    else PG_MAXLFQ_ABSENT)
    try:
        import pyarrow.compute as pc
        import pyarrow.parquet as pq
    except ImportError:
        pq = None
    if pq is not None:
        pf = pq.ParquetFile(report)
        if col not in pf.schema_arrow.names:
            return PG_MAXLFQ_ABSENT
        for batch in pf.iter_batches(columns=[col], batch_size=1 << 16):
            if pc.any(pc.greater(batch.column(0), 0)).as_py() is True:
                return PG_MAXLFQ_PRESENT
        return PG_MAXLFQ_ABSENT
    try:
        import polars as pl
    except ImportError:
        return None
    if col not in pl.read_parquet_schema(report):
        return PG_MAXLFQ_ABSENT
    hit = (pl.scan_parquet(report).select((pl.col(col) > 0).fill_null(False).any())
           .collect().item())
    return PG_MAXLFQ_PRESENT if hit else PG_MAXLFQ_ABSENT


def pg_maxlfq_answer(report):
    """pg_maxlfq_state() for the record and for build_maxlfq.R: {pg_maxlfq, why, refusal} --
    pg_maxlfq None (and `why` says why) when it could not be checked; `refusal` is
    MAXLFQ_REFUSAL for an "absent" report, else None. Never raises: what failed is said."""
    try:
        state = pg_maxlfq_state(report)
    except Exception as e:  # noqa: BLE001 -- said in `why`, and the record keeps it
        return {"pg_maxlfq": None, "refusal": None,
                "why": f"not checked: {report} could not be read ({type(e).__name__}: {e})"}
    if state is None:
        return {"pg_maxlfq": None, "refusal": None,
                "why": f"not checked: {sys.executable} has neither pyarrow nor polars to read "
                       f"{report}"}
    return {"pg_maxlfq": state, "refusal": MAXLFQ_REFUSAL if state == PG_MAXLFQ_ABSENT else None,
            "why": f"any PG.MaxLFQ > 0 in {os.path.basename(str(report))}: "
                   f"{'yes' if state == PG_MAXLFQ_PRESENT else 'no'}"}


def record(report, first_pass, provenance, exit_code):
    """Make `first_pass` the report, write RECORD, and set search_provenance.json `report_used`
    and `dda_first_pass` (when `provenance` names a file; the inline route has none yet, and
    main() writes the returned record into it). Returns the record. Raises ValueError when the
    final report exists or the first pass is missing or empty -- never a fallback then."""
    if nonempty(report):
        raise ValueError(f"{report} exists: the second pass wrote it, so it is the report")
    if not nonempty(first_pass):
        raise ValueError(f"{first_pass} is missing or empty: there is no first pass to use")
    if provenance is not None and not os.path.isfile(provenance):
        raise ValueError(f"{provenance} does not exist, and the Methods, AUDIT.md and the watcher "
                         f"read only it, so the fallback could not be recorded where it is read")
    # read before anything is copied: an unreadable record stops the fallback, not half of it
    prov = _load(provenance)
    # whether the first pass has engine protein quantities (pg_maxlfq_state, the one definition):
    # a --dda first pass has none, and build_maxlfq.R refuses --method maxlfq on it
    quant = pg_maxlfq_answer(first_pass)
    copied = {os.path.basename(report): os.path.basename(first_pass)}
    _copy(first_pass, report)
    for suffix in COPIED_SUFFIXES:
        src, dst = companion(first_pass, suffix), companion(report, suffix)
        if os.path.isfile(src):
            _copy(src, dst)
            copied[os.path.basename(dst)] = os.path.basename(src)
    rec = {"report_used": REPORT_USED, "exit_code": exit_code,
           "report": os.path.abspath(report), "first_pass": os.path.abspath(first_pass),
           "first_pass_sha256": _sha256(first_pass), "copied": copied,
           "at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
           "pg_maxlfq": quant["pg_maxlfq"], "pg_maxlfq_why": quant["why"],
           "by": "dda_first_pass.py record (the search job)"}
    rec["caveat"] = caution({"report_used": REPORT_USED, "dda_first_pass": rec})
    write_file(record_path(report), json.dumps(rec, indent=2) + "\n")
    rec["provenance_updated"] = record_in_provenance(provenance, rec, prov)
    return rec


def _load(path):
    """search_provenance.json as a dict, or None when there is no such file."""
    if not path or not os.path.isfile(path):
        return None
    with open(path) as fh:
        prov = json.load(fh)
    if not isinstance(prov, dict):
        raise ValueError(f"{path} is not a JSON object")
    return prov


def record_in_provenance(path, rec, prov=None):
    """search_provenance.json: `report_used` and `dda_first_pass`. False when there is no such
    file (the record is still in RECORD)."""
    prov = prov if prov is not None else _load(path)
    if prov is None:
        return False
    prov["report_used"] = rec["report_used"]
    prov["dda_first_pass"] = {k: v for k, v in rec.items() if k != "provenance_updated"}
    write_file(path, json.dumps(prov, indent=2))
    return True


def _no_pg_maxlfq(prov):
    """True when a fallen-back search's record says its first pass has no engine protein
    quantities (pg_maxlfq "absent"); a record that does not say (one from before T23b, or not
    checked) claims nothing."""
    if not isinstance(prov, dict) or prov.get("report_used") != REPORT_USED:
        return False
    rec = prov.get("dda_first_pass") if isinstance(prov.get("dda_first_pass"), dict) else {}
    return rec.get("pg_maxlfq") == PG_MAXLFQ_ABSENT


def caution(prov, methods=False):
    """CAVEAT for a search provenance whose `report_used` is "first_pass", else None (every
    other search, including any from before 2.14.2) -- with NO_PG_MAXLFQ after it when its
    dda_first_pass record says the first pass has no engine protein quantities (pg_maxlfq
    "absent"). `methods`: the Methods' search paragraph, which states only the fact
    (NO_PG_MAXLFQ_FACT); what the DE did about it is de_sentence()'s, in the DE paragraph."""
    if not isinstance(prov, dict) or prov.get("report_used") != REPORT_USED:
        return None
    if not _no_pg_maxlfq(prov):
        return CAVEAT
    return CAVEAT + " " + (NO_PG_MAXLFQ_FACT if methods else NO_PG_MAXLFQ)


def de_sentence(prov, de_prov):
    """The Methods' DE-paragraph sentence for a first-pass report with no engine protein
    quantities, when the DE that ran read its precursors with DPC-Quant (de_provenance.json
    pipeline_id "dpc", the record, never assumed) -- else None: no DE, or another pipeline
    (MaxLFQ refuses such a report), or a search that did not fall back."""
    if not _no_pg_maxlfq(prov) or not isinstance(de_prov, dict):
        return None
    return DE_DPC_SENTENCE if de_prov.get("pipeline_id") == "dpc" else None


def caution_for(out_dir):
    """caution() for a search output folder, from its search_provenance.json -- for the status
    views that have the folder (notify_slack.py's job-end post, checkpoint.py status). None
    when there is no record or nothing fell back; an unreadable record is said, not dropped."""
    path = os.path.join(out_dir, "search_provenance.json")
    try:
        return caution(_load(path))
    except (OSError, ValueError) as e:
        return (f"CAUTION: {path} could not be read ({type(e).__name__}: {e}), so whether this "
                f"--dda search fell back to its first-pass report is not known")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("record", help="make the first pass the report, and record it")
    r.add_argument("--report", required=True)
    r.add_argument("--first-pass", required=True)
    r.add_argument("--provenance", required=True,
                   help="search_provenance.json to record it in (it must exist)")
    r.add_argument("--exit-code", type=int, required=True, help="DIA-NN's exit status")
    c = sub.add_parser("caution", help="print the caution for a search provenance, or nothing")
    c.add_argument("--provenance", required=True)
    f = sub.add_parser("forget", help="a new run of the job: drop an earlier run's fallback "
                                      "from search_provenance.json")
    f.add_argument("--provenance", required=True)
    g = sub.add_parser("pg-maxlfq", help="does a report have engine protein quantities (any "
                                         "PG.MaxLFQ > 0)? JSON {pg_maxlfq, why, refusal}, for "
                                         "build_maxlfq.R")
    g.add_argument("--report", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "pg-maxlfq":
        print(json.dumps(pg_maxlfq_answer(a.report)))
        return 0
    if a.cmd == "forget":
        try:
            if forget(a.provenance):
                print(f"[dda_first_pass] an earlier run of this job fell back to its first pass; "
                      f"that record is taken out of {a.provenance} (this run records its own)")
        except (OSError, ValueError) as e:
            sys.exit(f"[dda_first_pass] FAILED: {a.provenance} could not be read or rewritten "
                     f"({type(e).__name__}: {e}), so this run could not start from a record of "
                     f"its own. Regenerate the job with run_search.py --sbatch.")
        return 0
    if a.cmd == "caution":
        try:
            with open(a.provenance) as fh:
                prov = json.load(fh)
        except FileNotFoundError:
            return 0                     # no record: nothing to say
        except (OSError, ValueError) as e:
            sys.exit(f"[dda_first_pass] {a.provenance} could not be read ({type(e).__name__}: "
                     f"{e})")
        msg = caution(prov)
        if msg:
            print(msg)
        return 0
    try:
        rec = record(a.report, a.first_pass, a.provenance, a.exit_code)
    except (OSError, ValueError) as e:
        sys.exit(f"[dda_first_pass] FAILED: the first-pass report could not be used ({e}). "
                 f"DIA-NN exited {a.exit_code} without {a.report}.")
    print(f"{MARKER}: DIA-NN exited {a.exit_code} without writing {a.report}; its first-pass "
          f"report {a.first_pass} (no match-between-runs) is used instead. Recorded in "
          f"{record_path(a.report)}"
          + (" and search_provenance.json (report_used: first_pass)."
             if rec["provenance_updated"] else
             "; there is no search_provenance.json to record it in."))
    print(f"CAUTION: {rec['caveat']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
