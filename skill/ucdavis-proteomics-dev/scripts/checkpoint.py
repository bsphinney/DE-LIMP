#!/usr/bin/env python3
"""
checkpoint.py -- crash/disconnect recovery for cluster runs.

A search is submitted to SLURM and then takes hours. The jobs survive the user
closing their laptop; the *conversation* does not. This records what was
submitted and what still has to happen, so a fresh session can pick the run back
up instead of starting over or, worse, silently re-running a finished search.

Writes two files in the session directory:
  RECOVERY.md    human-readable: what stage, which jobs, what to do next
  .recovery.json machine-readable: the same, for `status` to re-query

Usage
  # after submitting jobs
  checkpoint.py record --session <dir> --stage search \
      --jobs 19477495,19477496 --desc "DIA-NN 5-step chain" \
      --watch-job 19477500 --watch-log <path> \
      --next "Rscript run_de.R --input <report> --metadata <conditions> --outdir <tables>"

  # on reconnect — queries SLURM live and says what is left
  checkpoint.py status --session <dir>

  # find interrupted runs anywhere under a directory
  checkpoint.py find --base ~/            (or --base /quobyte/proteomics-grp/...)
"""
import argparse, json, os, subprocess, sys, glob, datetime

REC_JSON = ".recovery.json"
REC_MD = "RECOVERY.md"

TERMINAL_OK = {"COMPLETED"}
TERMINAL_BAD = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL",
                "BOOT_FAIL", "DEADLINE", "PREEMPTED"}


def _load(session):
    p = os.path.join(session, REC_JSON)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def _save(session, data):
    with open(os.path.join(session, REC_JSON), "w") as fh:
        json.dump(data, fh, indent=2)
    _write_md(session, data)


def sacct_states(job_ids):
    """Map job id -> state. Tolerates sacct being absent (returns UNKNOWN)."""
    if not job_ids:
        return {}
    try:
        out = subprocess.run(
            ["sacct", "-j", ",".join(str(j) for j in job_ids), "-n", "-P",
             "--format=JobID,State"],
            capture_output=True, text=True, timeout=60).stdout
    except Exception:
        return {str(j): "UNKNOWN" for j in job_ids}
    states = {}
    for line in out.splitlines():
        if "|" not in line:
            continue
        jid, st = line.split("|")[:2]
        if "." in jid:                      # skip .batch / .extern steps
            continue
        base = jid.split("_")[0]
        st = st.split()[0]
        # for an array, roll the tasks up: bad beats running beats completed
        prev = states.get(base)
        if prev is None:
            states[base] = st
        elif prev in TERMINAL_OK and st not in TERMINAL_OK:
            states[base] = st
        elif st in TERMINAL_BAD:
            states[base] = st
    return states


def _rollup(states):
    vals = list(states.values())
    if not vals:
        return "UNKNOWN"
    if any(v in TERMINAL_BAD for v in vals):
        return "FAILED"
    if any(v not in TERMINAL_OK for v in vals):
        return "RUNNING"
    return "COMPLETED"


def _write_md(session, data):
    L = []
    L.append("# Recovery — how to pick this analysis back up\n")
    L.append("You can close your terminal. **Cluster jobs keep running without you.**")
    L.append("When you come back, start Claude Code again in this folder and say:\n")
    L.append("> *resume this analysis*\n")
    L.append("Claude will read this file, check the jobs, and carry on from wherever it got to.")
    L.append("To check by hand instead:\n")
    L.append("```\npython3 checkpoint.py status --session %s\n```\n" % session)
    L.append("| field | value |")
    L.append("|---|---|")
    L.append("| session | `%s` |" % session)
    L.append("| last updated | %s |" % data.get("updated", "?"))
    if data.get("report"):
        L.append("| expected search output | `%s` |" % data["report"])
    L.append("")
    L.append("## Stages\n")
    L.append("| stage | jobs | recorded state | description |")
    L.append("|---|---|---|---|")
    for s in data.get("stages", []):
        L.append("| %s | %s | %s | %s |" % (
            s["stage"], ", ".join(str(j) for j in s.get("jobs", [])) or "—",
            s.get("state", "submitted"), s.get("desc", "")))
    L.append("")
    nxt = [s for s in data.get("stages", []) if s.get("next")]
    if nxt:
        L.append("## What still has to happen\n")
        L.append("Once the jobs above show COMPLETED, run:\n")
        L.append("```bash")
        for s in nxt:
            L.append(s["next"])
        L.append("```")
    L.append("")
    L.append("## If a job failed or hung\n")
    L.append("- A task still RUNNING while its siblings finished long ago is **hung**, not slow.")
    L.append("  Cancel just that one (`scancel <jobid>_<task>`) and resubmit that array index.")
    L.append("- The 5-step chain resumes from the earliest incomplete step; finished `.quant`")
    L.append("  files and the predicted library are reused, so you never redo the whole cohort.")
    L.append("- If one file hangs twice, drop it and continue with the rest — note it in the report.")
    with open(os.path.join(session, REC_MD), "w") as fh:
        fh.write("\n".join(L) + "\n")


def cmd_record(a):
    os.makedirs(a.session, exist_ok=True)
    data = _load(a.session) or {"session": a.session, "stages": []}
    jobs = [j for j in (a.jobs or "").split(",") if j]
    entry = {"stage": a.stage, "jobs": jobs, "desc": a.desc or "",
             "state": "submitted", "next": a.next,
             "watch_job": a.watch_job, "watch_log": a.watch_log}
    data["stages"] = [s for s in data["stages"] if s["stage"] != a.stage] + [entry]
    if a.report:
        data["report"] = a.report
    data["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    _save(a.session, data)
    print(json.dumps({"recorded": a.stage, "jobs": jobs,
                      "recovery_md": os.path.join(a.session, REC_MD)}, indent=2))


def cmd_status(a):
    data = _load(a.session)
    if not data:
        print(json.dumps({"error": "no %s in %s" % (REC_JSON, a.session)}))
        sys.exit(1)
    all_jobs = [j for s in data["stages"] for j in s.get("jobs", [])]
    states = sacct_states(all_jobs)
    out = {"session": a.session, "stages": [], "done": True, "any_failed": False}
    for s in data["stages"]:
        st = _rollup({j: states.get(j, "UNKNOWN") for j in s.get("jobs", [])}) \
             if s.get("jobs") else s.get("state", "submitted")
        s["state"] = st
        if st != "COMPLETED":
            out["done"] = False
        if st == "FAILED":
            out["any_failed"] = True
        out["stages"].append({"stage": s["stage"], "jobs": s.get("jobs", []),
                              "state": st, "next": s.get("next")})
    rep = data.get("report")
    out["report_exists"] = bool(rep and os.path.exists(rep))
    out["report"] = rep
    # the search is only really done when the output file exists — SLURM State
    # alone has lied before (a segfault masked by a trailing echo)
    if rep and not out["report_exists"]:
        out["done"] = False
    data["updated"] = datetime.datetime.now().isoformat(timespec="seconds")
    _save(a.session, data)
    out["next_commands"] = [s["next"] for s in data["stages"]
                            if s.get("next") and s["state"] == "COMPLETED"]
    print(json.dumps(out, indent=2))


def cmd_find(a):
    hits = []
    for p in glob.glob(os.path.join(a.base, "**", REC_JSON), recursive=True):
        try:
            with open(p) as fh:
                d = json.load(fh)
            hits.append({"session": os.path.dirname(p),
                         "updated": d.get("updated"),
                         "stages": [s["stage"] for s in d.get("stages", [])]})
        except Exception:
            continue
    hits.sort(key=lambda h: h.get("updated") or "", reverse=True)
    print(json.dumps({"found": len(hits), "sessions": hits}, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("record", help="record a submitted stage")
    r.add_argument("--session", required=True)
    r.add_argument("--stage", required=True)
    r.add_argument("--jobs", help="comma-separated SLURM job ids")
    r.add_argument("--desc", default="")
    r.add_argument("--next", help="command to run once this stage completes")
    r.add_argument("--report", help="output file that must exist for the stage to count as done")
    r.add_argument("--watch-job"), r.add_argument("--watch-log")
    r.set_defaults(func=cmd_record)

    s = sub.add_parser("status", help="re-query SLURM and report what is left")
    s.add_argument("--session", required=True)
    s.set_defaults(func=cmd_status)

    f = sub.add_parser("find", help="find interrupted runs under a directory")
    f.add_argument("--base", required=True)
    f.set_defaults(func=cmd_find)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
