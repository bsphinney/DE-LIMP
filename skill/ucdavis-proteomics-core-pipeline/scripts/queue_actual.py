#!/usr/bin/env python3
"""
queue_actual.py -- where a search's jobs actually ran, from Slurm, into search_provenance.json.

search_provenance.json `queue` is what was REQUESTED at generation (--partition/--account/--qos,
else the detected rule). A pending job can be moved afterwards (scontrol update, the controller's
queue move to publicgrp/low), and the record then named a queue the search never used (2.15
acceptance). At the search's last job (notify_slack.search_done, the job-end hook) this asks
sacct for the partition, account and QOS of every job of the search and writes them beside the
request as `queue_actual`, never over it:

  "queue_actual": {"jobs": [{"id", "partition", "account", "qos", "state"}],
                   "partition": [..], "account": [..], "qos": [..],
                   "moved": true | false | null,      # null: nothing was requested to compare
                   "source": "sacct", "recorded_at": "..."}
  or {"unavailable": "<why>"} -- said, never left out quietly.

  python3 queue_actual.py --provenance <out>/search_provenance.json [--jobs-file <out>/jobs.txt]
                          [--job ID ...]
Job ids: --jobs-file (submit.sh's jobs.txt), --job, and $SLURM_JOB_ID. Exit 0 always: a record that
cannot be made must never fail a finished search.
"""
import argparse
import datetime
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = "queue_actual"
SACCT_CANDIDATES = ("/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin/sacct",
                    "/usr/bin/sacct", "/usr/local/bin/sacct", "/opt/slurm/bin/sacct")


def _sacct():
    """The sacct to run: $QUEUE_ACTUAL_SACCT (a stand-in), PATH, then HPC@UCD's fixed paths."""
    return (os.environ.get("QUEUE_ACTUAL_SACCT") or shutil.which("sacct")
            or next((c for c in SACCT_CANDIDATES if os.access(c, os.X_OK)), None))


def job_ids(jobs_file=None, extra=()):
    ids = []
    if jobs_file:
        try:
            with open(jobs_file, encoding="utf-8", errors="replace") as fh:
                ids += re.findall(r"\d{5,}", fh.read())
        except OSError:
            pass
    ids += [str(x) for x in extra if x and re.fullmatch(r"\d{5,}", str(x))]
    return list(dict.fromkeys(ids))


def query(ids):
    """([job], None) or (None, why): one row per job (an array collapses to its job)."""
    if not ids:
        return None, "no SLURM job ids known (no jobs.txt, no $SLURM_JOB_ID)"
    exe = _sacct()
    if not exe:
        return None, "sacct is not available here"
    try:
        p = subprocess.run([exe, "-j", ",".join(ids[:300]), "-X", "-n", "-P",
                            "--format=JobID,Partition,Account,QOS,State"],
                           capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=30)
    except (OSError, subprocess.SubprocessError) as e:
        return None, f"sacct failed: {e}"
    jobs, seen = [], set()
    for ln in p.stdout.splitlines():
        f = ln.split("|")
        if len(f) < 5:
            continue
        base = f[0].split("_")[0].split(".")[0]
        if base in seen:
            continue
        seen.add(base)
        jobs.append({"id": base, "partition": f[1] or None, "account": f[2] or None,
                     "qos": f[3] or None, "state": f[4].split(" ")[0]})
    return (jobs, None) if jobs else (None, (p.stderr.strip() or "sacct returned no rows")[:200])


def _requested(queue):
    q = queue if isinstance(queue, dict) else {}
    return {k: q.get(k) for k in ("partition", "account", "qos") if q.get(k)}


def build(jobs, queue):
    """The record for `jobs` against the requested `queue` (search_provenance.json's `queue`)."""
    def distinct(k):
        return sorted({j[k] for j in jobs if j.get(k)})
    rec = {"jobs": jobs, "partition": distinct("partition"), "account": distinct("account"),
           "qos": distinct("qos"), "source": "sacct",
           "recorded_at": datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0)
           .isoformat()}
    want = _requested(queue)
    if want:
        rec["moved"] = any(j.get(k) != v for j in jobs for k, v in want.items() if j.get(k))
        rec["requested"] = want
    else:
        rec["moved"] = None
    return rec


def record(provenance, jobs_file=None, extra=()):
    """Merge `queue_actual` into `provenance` (in place, atomically). -> the record written, or
    None when the file cannot be read or written (said on stderr)."""
    try:
        with open(provenance, encoding="utf-8") as fh:
            prov = json.load(fh)
        if not isinstance(prov, dict):
            raise ValueError("not a JSON object")
    except (OSError, ValueError) as e:
        sys.stderr.write(f"[queue_actual] NOTE: {provenance} cannot be read ({e}); "
                         "the queue the jobs ran in is not recorded\n")
        return None
    jobs, why = query(job_ids(jobs_file, extra))
    rec = build(jobs, prov.get("queue")) if jobs else {"unavailable": why}
    prov[KEY] = rec
    try:
        with open(provenance + ".tmp", "w", encoding="utf-8") as fh:
            json.dump(prov, fh, indent=2)
        os.replace(provenance + ".tmp", provenance)
    except OSError as e:
        sys.stderr.write(f"[queue_actual] NOTE: {provenance} cannot be written ({e})\n")
        return None
    return rec


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--provenance", required=True, help="the search's search_provenance.json")
    ap.add_argument("--jobs-file", help="submit.sh's jobs.txt (every job of the chain)")
    ap.add_argument("--job", action="append", default=[], help="a job id (repeatable)")
    a = ap.parse_args(argv)
    rec = record(a.provenance, a.jobs_file, list(a.job) + [os.environ.get("SLURM_JOB_ID")])
    if rec is not None:
        print(json.dumps(rec))
    return 0


if __name__ == "__main__":
    sys.exit(main())
