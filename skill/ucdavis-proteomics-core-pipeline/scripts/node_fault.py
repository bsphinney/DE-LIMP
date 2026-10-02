#!/usr/bin/env python3
"""
node_fault.py  --  A compute node that cannot reach the storage is a NODE fault, not a search
failure: say so, and retry the step on another node.

On 2026-09-28 the chain's step 3 died in 2 s on hive-dc-7-4-50: `diann-linux: Transport endpoint
is not connected` (that node's /quobyte FUSE mount had dropped), then `FAILED: DIA-NN exited 0 but
did not write the empirical spectral library` -- the message blamed DIA-NN, which never ran. Steps
4 and 5 sat DependencyNeverSatisfied until someone noticed. The same day a job on hive-dc-7-4-54
hung 16 min expanding a raw-data glob that the login node listed instantly; resubmitted with
--exclude it ran in 12 s. So:

  preflight_lines(checks)   bash every generated chain job runs FIRST (diann_parallel.py, through
                            notify_slack.wrap_job_script): can THIS node reach the engine, the
                            FASTA, the raw data and the search folder? Each check is bounded (a
                            hung mount never answers) and a failure prints
                            `NODE_FAULT: node=<n> check=<what> path=<p> detail=<why>` and exits
                            NODE_FAULT_EXIT -- before DIA-NN, so nothing is blamed on it.
  classify(...)             a FAILED job's log tail, sacct state and exit codes -> a node fault,
                            or None. The signatures are I/O errors from the storage, never
                            anything DIA-NN reports about the data. watch_run.sh asks this.
  retry                     resubmit the failed step -- and every step waiting on it, which
                            afterok has left DependencyNeverSatisfied -- with --exclude=<node>,
                            at most MAX_NODE_RETRIES times per step. Recorded in
                            <out>/node_faults.json and search_provenance.json `node_faults`.

What is NOT a node fault here: an ESTALE that step 1b's probe absorbed. probe_window.py rereads
its own log on ESTALE, retries once and then falls back (probe_fallback.py) -- the job succeeds,
and only a FAILED job is classified, so the two never act on the same failure. OOM, a time limit
and a cancel keep their own classes (watch_run.sh); a retry of those on another node would fail
the same way.

Usage:
  node_fault.py classify [--state FAILED] [--exit-codes 75:0,1:0] [--nodes n1,n2] < log_tail
  node_fault.py retry --out <search out dir> [--job <failed job id>] [--node <name>] [--dry-run]
"""
import argparse
import datetime
import json
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import time

# EX_TEMPFAIL (sysexits.h): "a temporary failure ... the user is invited to retry". What a job
# exits with when its node check fails; sacct shows it as ExitCode 75:0 even when the node could
# not write the job's log.
NODE_FAULT_EXIT = 75
# Retries per step. A third failure of the same step on yet another node is not one bad node.
MAX_NODE_RETRIES = 2
# Failed tasks spread over more nodes than this are the storage (or the data), not a node.
MAX_BAD_NODES = 3
# Seconds a node check waits for an answer. A healthy stat takes milliseconds; a dropped FUSE
# mount can block for ever (16 min on hive-dc-7-4-54 before it was cancelled).
NODE_CHECK_S = 120
RECORD = "node_faults.json"

# (name, pattern, what it means). I/O errors from the storage -- never a DIA-NN message.
SIGNATURES = (
    ("preflight", re.compile(r"NODE_FAULT: [^\n]*"),
     "the job's node check failed before the search started"),
    ("enotconn", re.compile(r"[^\n]*Transport endpoint is not connected[^\n]*"),
     "the node's network-filesystem mount had dropped (ENOTCONN)"),
    ("estale", re.compile(r"[^\n]*Stale file handle[^\n]*"),
     "the node's handle on the network filesystem went stale (ESTALE)"),
    ("eio", re.compile(r"[^\n]*Input/output error[^\n]*"),
     "an I/O error reading the storage (EIO)"),
)
# A job SLURM ended for a reason of its own: not retried as a node fault.
NOT_NODE_STATES = ("OUT_OF_MEMORY", "TIMEOUT", "CANCELLED", "PREEMPTED", "DEADLINE")
FAILED_STATES = ("FAILED", "NODE_FAIL", "BOOT_FAIL")
# a line step 1b's probe wrote or relayed (probe_window.py: "[probe_window] ...", "[probe 2] ...";
# diann_parallel.probe_attempts: "[probe] ...")
PROBE_LINE_RE = re.compile(r"\s*\[probe(?:_window| \d+)?\] ")
# What to do about a failed task that is NOT a node fault (the watcher's own classes)
NOT_NODE_FIX = {"OUT_OF_MEMORY": "raise the memory (--mem-per-file / --mem) and resubmit it",
                "TIMEOUT": "raise the time limit (--time-per-file / --time) and resubmit it",
                "CANCELLED": "find out who cancelled it and why, then resubmit it",
                "PREEMPTED": "resubmit it (it was preempted)",
                "DEADLINE": "resubmit it with a later deadline"}


def _ended_badly(r):
    return (r["state"] in FAILED_STATES + NOT_NODE_STATES
            or r["exit"].split(":")[0] == str(NODE_FAULT_EXIT))


def classify(text, state=None, exit_codes=(), nodes=()):
    """A node fault in a failed job, or None. `text` is its log tail, `state` its sacct state,
    `exit_codes` sacct's ExitCode values ("75:0"), `nodes` where it ran. The log's own line is
    the evidence when there is one; NODE_FAULT_EXIT (the node check) and NODE_FAIL (SLURM's own
    verdict) need no log -- a node that lost its mount may not have been able to write one."""
    state = state.split()[0].rstrip("+") if state and state.split() else ""
    nodes = sorted({n for n in nodes if n and n not in ("None", "None assigned")})
    if state in NOT_NODE_STATES:
        return None
    rec = None
    # What step 1b's probe ABSORBED is not this job's failure: it rereads its log on ESTALE,
    # retries, falls back, and tags every line it writes or relays "[probe ...]". A step 1b that
    # did that and then failed for another reason (a run that logged no radius) was read as a
    # node fault and retried on another node (review of 2.10). Its node problems are caught by
    # the node check, which is not tagged.
    text = "\n".join(ln for ln in (text or "").splitlines() if not PROBE_LINE_RE.match(ln))
    for name, rx, meaning in SIGNATURES:
        found = rx.findall(text or "")
        if found:
            rec = {"signature": name, "meaning": meaning, "evidence": found[0].strip()[:300],
                   "count": len(found)}
            break
    if rec is None:
        if any(str(c).split(":")[0] == str(NODE_FAULT_EXIT) for c in exit_codes):
            rec = {"signature": "preflight", "meaning": SIGNATURES[0][2],
                   "evidence": f"exit code {NODE_FAULT_EXIT} (the node check), no log line read"}
        elif state in ("NODE_FAIL", "BOOT_FAIL"):
            rec = {"signature": "slurm_node_fail", "evidence": f"sacct state {state}",
                   "meaning": "SLURM itself marked the node as failed"}
        else:
            return None
    rec["nodes"] = nodes
    rec["say"] = say(rec)
    return rec


def say(rec):
    """The sentence for the user: a node problem, named as one."""
    where = ", ".join(rec.get("nodes") or []) or "the compute node"
    return (f"Node problem, not a DIA-NN or data problem: {where} could not reach the storage "
            f"({rec['meaning']}; {rec['evidence']}). The step is resubmitted on another node "
            f"(node_fault.py retry, at most {MAX_NODE_RETRIES} times per step).")


# --------------------------------------------------------------------------- the node check
def preflight_lines(checks, timeout_s=NODE_CHECK_S):
    """Bash for the start of a job: each (test flag, path, what) in `checks` must answer within
    `timeout_s` on this node, or the job prints NODE_FAULT and exits NODE_FAULT_EXIT.

    Each check runs in a background subshell that this shell polls, instead of `timeout`: a
    process stuck on a dead FUSE mount is in uninterruptible sleep, so `timeout` would wait for it
    as long as the mount does. The job exits without waiting; SLURM reaps what is left. The paths
    were all there on the login node when the job was generated, so one missing HERE is this
    node's view of the storage."""
    q = shlex.quote
    lines = [
        "# Node check (node_fault.py): can THIS node reach what the job needs? A node whose",
        "# storage mount has dropped fails here, named as a node fault (exit "
        f"{NODE_FAULT_EXIT}), not later as a search failure.",
        "_nf_fail() {",
        '  echo "NODE_FAULT: node=$(hostname -s 2>/dev/null || hostname) $*" >&2',
        '  echo "This is a problem with the compute node, not with DIA-NN or the data: resubmit '
        'on another node (node_fault.py retry, or --exclude=<node>)." >&2',
        f"  exit {NODE_FAULT_EXIT}",
        "}",
        "_nf_check() {",
        '  local flag="$1" path="$2" what="$3" t=0 pid err',
        '  err="$(mktemp 2>/dev/null || echo /dev/null)"',
        '  ( ls -ld -- "$path" >/dev/null 2>"$err" && test "$flag" "$path" ) &',
        "  pid=$!",
        '  while kill -0 "$pid" 2>/dev/null; do',
        f'    if [ "$t" -ge {int(timeout_s)} ]; then _nf_fail "check=$what path=$path '
        f'detail=no answer in {int(timeout_s)} s (a hung mount)"; fi',
        "    sleep 1; t=$((t + 1))",
        "  done",
        '  if wait "$pid"; then [ "$err" = /dev/null ] || rm -f "$err"; return 0; fi',
        '  local why; why="$(head -c 300 "$err" 2>/dev/null | tr \'\\n\' \' \')"',
        '  [ "$err" = /dev/null ] || rm -f "$err"',
        '  _nf_fail "check=$what path=$path detail=${why:-fails test $flag here, although it '
        'passed on the login node when this job was generated}"',
        "}",
    ]
    for flag, path, what in checks:
        lines.append(f"_nf_check {q(flag)} {q(path)} {q(what)}")
    lines.append('echo "node check: $(hostname -s 2>/dev/null || hostname) reaches '
                 + ", ".join(w for _, _, w in checks) + '"')
    return lines


def chain_checks(diann, out, fasta, raws, max_dirs=8):
    """What a DIA-NN chain job needs from the storage, as (test flag, path, what) for
    preflight_lines(): every absolute path in the DIA-NN command that exists now (the binary, an
    Apptainer image), the search folder (writable), the FASTA and the raw files' folders. Only
    paths that exist where this runs (the login node) -- a check must not fail on a path that was
    never there."""
    checks, seen = [], set()

    def add(flag, path, what):
        if path and path not in seen and os.path.exists(path):
            seen.add(path)
            checks.append((flag, path, what))
    try:
        words = shlex.split(diann)
    except ValueError:
        words = diann.split()
    for w in words:
        if w.startswith("/") and os.path.isfile(w):
            add("-x" if os.access(w, os.X_OK) else "-r", w, "the DIA-NN engine")
    add("-w", out, "the search folder")
    add("-r", fasta, "the FASTA")
    dirs = []
    for r in raws:
        d = os.path.dirname(os.path.abspath(r.rstrip("/")))
        if d not in dirs:
            dirs.append(d)
    for d in dirs[:max_dirs]:
        add("-r", d, "the raw data")
    return checks


# ------------------------------------------------------------------------------ SLURM helpers
def _run(argv, cwd=None):
    try:
        r = subprocess.run(argv, capture_output=True, text=True, timeout=120, cwd=cwd)
        return r.returncode, r.stdout, r.stderr
    except (OSError, subprocess.SubprocessError) as e:
        return None, "", f"{type(e).__name__}: {e}"


def sacct_rows(job):
    """[{id, task, state, exit, node}] for a job (every array task), from sacct."""
    rc, out, err = _run(["sacct", "-j", str(job), "-X", "-n", "-P",
                         "-o", "JobID,State,ExitCode,NodeList"])
    if rc != 0:
        raise RuntimeError(f"sacct -j {job} failed: {(err or out).strip()[:200]}")
    rows = []
    for ln in out.splitlines():
        f = ln.strip().split("|")
        if len(f) < 4 or "." in f[0]:
            continue
        jid, _, task = f[0].partition("_")
        rows.append({"id": jid, "task": task or None, "state": f[1].split()[0] if f[1] else "",
                     "exit": f[2], "node": f[3]})
    return rows


def parse_submit(path):
    """submit.sh as written by diann_parallel.py (the 5-step chain) or run_search.py (the two-job
    search): its jobs in submission order, [{var, script, deps}], the order of jobs.txt (the
    variables its printf writes) and the folder it submits from (its `cd`)."""
    with open(path) as fh:
        text = fh.read()
    jobs = []
    for m in re.finditer(r"^\s*(\w+)=\$\(sbatch --parsable((?:\s+--\S+)*)\s+(.+?)\)\s*$",
                         text, re.M):
        var, opts, script = m.group(1), m.group(2), m.group(3)
        deps = []
        for d in re.findall(r"--dependency=afterok:(\S+)", opts):
            deps += [x for x in d.split(":") if x]
        jobs.append({"var": var, "script": shlex.split(script)[0],
                     "deps": [x[1:] if x.startswith("$") else x for x in deps]})
    order = None
    for ln in text.splitlines():
        if ln.lstrip().startswith("printf") and "jobs.txt" in ln:
            order = re.findall(r"\$(\w+)", ln.split(">")[0])
    cd = re.search(r'^cd\s+"?([^"\n]+)"?\s*$', text, re.M)
    return jobs, order, (cd.group(1) if cd else None)


def _now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def _write_json(path, data):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
    os.replace(tmp, path)


def _log_tail(out, job, task=None, n=12000):
    """The tail of a failed job's log in the search folder (`<name>_<job>.log`, an array task's
    `<name>_<job>_<task>.log`), or ""."""
    want = f"_{job}_{task}.log" if task else f"_{job}.log"
    try:
        names = [f for f in os.listdir(out) if f.endswith(want)]
    except OSError:
        return ""
    for f in names:
        try:
            with open(os.path.join(out, f), errors="replace") as fh:
                return fh.read()[-n:]
        except OSError:
            continue
    return ""


def _refuse(msg, **extra):
    print(json.dumps(dict({"retried": False, "say": msg}, **extra), indent=2))
    sys.exit(3)


# One retry at a time per search folder. Two at once (two sessions on the shared account, or a
# resumed one beside the original) both read the old jobs.txt and both resubmit the step and its
# dependants: two copies of steps 3-5 writing one folder (review of 2.10). An mkdir lock IN the
# folder, because flock does not lock across HIVE's login nodes on /quobyte. A second retry waits
# for the first, then reads the jobs afresh -- and finds nothing left to retry.
LOCK = ".node_fault_retry.lock"
LOCK_WAIT_S = float(os.environ.get("NODE_FAULT_LOCK_WAIT_S", "180"))   # (tests shorten it)
LOCK_STALE_S = 900            # a retry takes seconds; a lock this old is a dead one's


class RetryLock:
    def __init__(self, out):
        self.path = os.path.join(out, LOCK)

    def __enter__(self):
        deadline = time.time() + LOCK_WAIT_S
        while True:
            try:
                os.mkdir(self.path)
            except FileExistsError:
                try:
                    age = time.time() - os.stat(self.path).st_mtime
                except FileNotFoundError:
                    continue                          # released between the two calls
                if age > LOCK_STALE_S:
                    # a dead retry's lock: moved aside (it is ours once renamed), then removed
                    aside = f"{self.path}.stale-{int(time.time())}-{os.getpid()}"
                    try:
                        os.rename(self.path, aside)
                        shutil.rmtree(aside, ignore_errors=True)
                    except OSError:
                        pass
                    continue
                if time.time() >= deadline:
                    _refuse(f"another node_fault.py retry of this search holds {self.path} "
                            f"(for {age:.0f} s): wait for it, then run watch_run.sh --all")
                time.sleep(1)
                continue
            with open(os.path.join(self.path, "owner"), "w") as fh:
                fh.write(f"{socket.gethostname()} pid {os.getpid()} {_now()}\n")
            return self

    def __exit__(self, *exc):
        shutil.rmtree(self.path, ignore_errors=True)
        return False


def retry(out, job=None, node=None, dry_run=False, force=False):
    """Resubmit the failed step of the search in `out` and the steps waiting on it, excluding
    the node(s) it failed on, holding RetryLock for the whole of it. See the module docstring."""
    out = os.path.abspath(out)
    if dry_run or not os.path.isdir(out):
        return _retry(out, job, node, dry_run, force)
    with RetryLock(out):
        return _retry(out, job, node, dry_run, force)


def _retry(out, job, node, dry_run, force):
    """retry()'s body: everything it reads (jobs.txt, the states) is read under the lock."""
    submit, jobs_txt = os.path.join(out, "submit.sh"), os.path.join(out, "jobs.txt")
    if not (os.path.isfile(submit) and os.path.isfile(jobs_txt)):
        _refuse(f"{out} has no submit.sh + jobs.txt: not a search submitted as a chain. "
                f"Resubmit the one job with `sbatch --exclude=<node> <script>`.")
    jobs, order, cd = parse_submit(submit)
    with open(jobs_txt) as fh:
        ids = [ln.strip() for ln in fh if ln.strip()]
    if not jobs or not order or len(order) != len(ids) or {j["var"] for j in jobs} != set(order):
        _refuse(f"{jobs_txt} does not match {submit} ({len(ids)} ids, "
                f"{len(order or [])} jobs): retry by hand (references/watcher.md, node_fault)")
    id_of = dict(zip(order, ids))
    by_var = {j["var"]: j for j in jobs}

    # which job failed, and where
    states = {}
    for v in order:
        try:
            states[v] = sacct_rows(id_of[v])
        except RuntimeError as e:
            _refuse(f"cannot read the jobs' states: {e}")
    if job:
        base = str(job).split("_")[0].split(".")[0]
        failed_var = next((v for v in order if id_of[v] == base), None)
        if not failed_var:
            _refuse(f"job {job} is not in {jobs_txt} ({', '.join(ids)})")
    else:
        failed_var = next((v for v in order if any(_ended_badly(r) for r in states[v])), None)
        if not failed_var:
            _refuse("no failed job in this search: nothing to retry")
    # EVERY task that ended badly -- an OUT_OF_MEMORY or TIMEOUT task beside a node fault too: the
    # steps after this one need all of its tasks (review of 2.10: an OOM task of the same array
    # was left out of the retry, step 5 ran afterok on the rest and the retry said complete)
    bad = [r for r in states[failed_var] if _ended_badly(r)]
    if not bad:
        _refuse(f"job {id_of[failed_var]} did not fail ({', '.join(sorted({r['state'] for r in states[failed_var]}))})")
    is_array = any(r["task"] for r in states[failed_var])
    busy = [r for r in states[failed_var] if r["state"] in ("RUNNING", "PENDING", "COMPLETING",
                                                              "REQUEUED", "RESIZING")]
    if busy:
        # The steps waiting on it need ALL its tasks: resubmitting the failed ones while the rest
        # still run would let the next step start without them.
        _refuse(f"{len(busy)} task(s) of job {id_of[failed_var]} are still running or queued: "
                "retry once the array has finished, so every failed task is retried together",
                failed_job=id_of[failed_var])

    # Is it a node fault? Judged TASK BY TASK, each on its own sacct row and its own log -- the
    # first bad task's evidence once stood for all of them, and a task that failed with DIA-NN's
    # own error had its node excluded and blamed.
    verdict = {}
    for r in bad:
        verdict[r["task"]] = classify(_log_tail(out, id_of[failed_var], r["task"]), r["state"],
                                      [r["exit"]], [r["node"]])
    faulty = [r for r in bad if verdict[r["task"]]]
    other = [r for r in bad if not verdict[r["task"]]]
    label = (lambda r: f"task {r['task']}" if r["task"] else f"job {r['id']}")
    others = [{"task": r["task"], "state": r["state"], "exit": r["exit"], "node": r["node"],
               "fix": NOT_NODE_FIX.get(r["state"], "read its log; fix the failure itself "
                                                   "(watch_run.sh names its class)")}
              for r in other]
    if other and not force:
        # A node retry cannot fix these, and the steps after this one need them: resubmitting
        # only the node-fault tasks and calling the chain complete would let step 5 build a
        # report without them. Nothing is resubmitted; each failure is said with its fix.
        _refuse(f"job {id_of[failed_var]} has failures a node retry cannot fix: "
                + "; ".join(f"{label(r)} {r['state']} (exit {r['exit']}) on {r['node']}"
                            for r in other)
                + (f". Node faults: {', '.join(label(r) for r in faulty)}" if faulty else "")
                + ". Nothing was resubmitted. Fix each as said under `not_node_faults`, then "
                  "resubmit this step's failed tasks and the steps after it (references/"
                  "watcher.md); --force retries them all as node faults anyway",
                failed_job=id_of[failed_var], not_node_faults=others,
                node_faults=[dict(task=r["task"], node=r["node"], **verdict[r["task"]])
                             for r in faulty])
    retried = bad if force else faulty
    nodes = sorted({r["node"] for r in (faulty or retried)
                    if r["node"] and r["node"] != "None assigned"} | ({node} if node else set()))
    tasks = sorted({r["task"] for r in retried if r["task"]},
                   key=lambda t: int(t) if t.isdigit() else 0)
    ev = next((verdict[r["task"]] for r in faulty), None)
    if len(nodes) > MAX_BAD_NODES:
        _refuse(f"the failed tasks ran on {len(nodes)} different nodes ({', '.join(nodes)}): that "
                "is the storage, or the data, not one bad node. Tell the user; check the HIVE "
                "status page before resubmitting.", nodes=nodes)

    rec_path = os.path.join(out, RECORD)
    try:
        with open(rec_path) as fh:
            record = json.load(fh)
    except FileNotFoundError:
        record = {"retries": []}
    except (OSError, ValueError) as e:
        _refuse(f"{rec_path} cannot be read ({e}), so how often this step was retried is not "
                "known: retry by hand")
    script = by_var[failed_var]["script"]
    done = [r for r in record.get("retries", []) if r.get("step") == script]
    if len(done) >= MAX_NODE_RETRIES:
        _refuse(f"{script} already failed on a node {len(done)} times (nodes "
                f"{', '.join(sorted({n for r in done for n in r.get('nodes', [])}))}) and now on "
                f"{', '.join(nodes) or 'another'}: not one bad node. Tell the user it is a "
                "cluster storage problem; do not resubmit again.", attempts=len(done))
    excluded = sorted({n for r in record.get("retries", []) for n in r.get("nodes", [])}
                      | set(nodes))

    # every job waiting on the failed one, transitively, in submission order
    downstream, grow = [], {failed_var}
    for j in jobs:
        if j["var"] != failed_var and set(j["deps"]) & grow:
            downstream.append(j["var"])
            grow.add(j["var"])
    for v in downstream:
        live = {r["state"] for r in states[v]}
        if live & {"RUNNING", "COMPLETED", "COMPLETING"}:
            _refuse(f"{by_var[v]['script']} (job {id_of[v]}) is {', '.join(sorted(live))} although "
                    f"it waits on the failed {script}: retry by hand")

    plan = []
    excl = f"--exclude={','.join(excluded)}" if excluded else None
    cwd = cd or out

    def argv(v, deps):
        a = ["sbatch", "--parsable"]
        if deps:
            a.append("--dependency=afterok:" + ":".join(deps))
        if excl:
            a.append(excl)
        if v == failed_var and is_array and tasks:
            a.append("--array=" + ",".join(tasks))
        return a + [by_var[v]["script"]]
    plan.append({"script": script, "old": id_of[failed_var], "argv": argv(failed_var, [])})
    if dry_run:
        for v in downstream:
            plan.append({"script": by_var[v]["script"], "old": id_of[v], "cancel": id_of[v],
                         "argv": argv(v, ["<new " + d + ">" if d in grow else id_of.get(d, d)
                                          for d in by_var[v]["deps"]])})
        print(json.dumps({"retried": False, "dry_run": True, "plan": plan, "excluded": excluded,
                          "node_fault": ev}, indent=2))
        return 0

    new_id = {}

    def submit_one(v, deps):
        rc, so, se = _run(argv(v, deps), cwd=cwd)
        nid = (so or "").strip().split(";")[0]
        if rc != 0 or not nid.isdigit():
            raise RuntimeError(f"sbatch {by_var[v]['script']} failed: {(se or so).strip()[:300]}")
        new_id[v] = nid
        return nid

    cancelled, errors = [], []
    try:
        submit_one(failed_var, [])
        for v in downstream:
            # pending for ever on the failed job (DependencyNeverSatisfied); one SLURM already
            # ended (kill_invalid_depend) needs no scancel
            if {r["state"] for r in states[v]} - {"CANCELLED", "FAILED"}:
                rc, so, se = _run(["scancel", id_of[v]])
                (cancelled if rc == 0 else errors).append(
                    id_of[v] if rc == 0 else f"scancel {id_of[v]}: {(se or so).strip()[:200]}")
            submit_one(v, [new_id.get(d) or id_of.get(d, d) for d in by_var[v]["deps"]])
    except RuntimeError as e:
        errors.append(str(e))
    resub = [{"script": by_var[v]["script"], "old": id_of[v], "new": new_id[v]}
             for v in [failed_var] + downstream if v in new_id]
    # jobs.txt: same order, the new ids in place of the old -- what watch_run.sh --all reads
    _ids = [new_id.get(v, id_of[v]) for v in order]
    with open(jobs_txt + ".tmp", "w") as fh:
        fh.write("\n".join(_ids) + "\n")
    os.replace(jobs_txt + ".tmp", jobs_txt)
    entry = {"at": _now(), "step": script, "attempt": len(done) + 1, "max": MAX_NODE_RETRIES,
             "failed_job": id_of[failed_var], "tasks": tasks or None, "nodes": nodes,
             "excluded": excluded, "evidence": ev,
             "evidence_per_task": {str(r["task"]): verdict[r["task"]] for r in retried
                                   if r["task"]} or None,
             "forced": bool(force and other), "resubmitted": resub,
             "cancelled": cancelled, "errors": errors}
    record.setdefault("retries", []).append(entry)
    _write_json(rec_path, record)
    _record_in_provenance(out, record["retries"])
    ck = _update_checkpoint(out, {id_of[v]: new_id[v] for v in new_id})
    # complete: EVERY task of the step that ended badly was resubmitted, and every step after it
    complete = len(resub) == 1 + len(downstream) and not errors and len(retried) == len(bad)
    msg = ((ev or {}).get("say") or "Retried on another node (--force).") + (
        f" Attempt {entry['attempt']} of {MAX_NODE_RETRIES} for {script}: "
        + "; ".join(f"{r['script']} {r['old']} -> {r['new']}" for r in resub)
        + f"; excluding {', '.join(excluded) or 'no node'}.")
    if not complete:
        msg += (" NOT every step could be resubmitted -- the chain is incomplete: "
                + "; ".join(errors) + ". Finish it by hand (references/watcher.md).")
    print(json.dumps({"retried": bool(resub), "complete": complete, "say": msg,
                      "node_fault": ev, "excluded": excluded, "resubmitted": resub,
                      "cancelled": cancelled, "errors": errors, "checkpoint": ck,
                      "record": rec_path, "watch": f"watch_run.sh --all {out}"}, indent=2))
    return 0 if complete else 1


def _record_in_provenance(out, retries):
    """search_provenance.json `node_faults`: the retries, beside what the search ran with."""
    p = os.path.join(out, "search_provenance.json")
    if not os.path.isfile(p):
        return False
    with open(p) as fh:
        prov = json.load(fh)
    prov["node_faults"] = retries
    _write_json(p, prov)
    return True


def _update_checkpoint(out, new_ids):
    """The session's .recovery.json (checkpoint.py) lists the job ids a resumed session checks:
    the old ones would read FAILED for ever. Replaced in place; None when there is none."""
    session = os.path.dirname(os.path.dirname(out))
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        import checkpoint
        data = checkpoint._load(session)
    except (ImportError, OSError, ValueError) as e:
        return f"not updated: {type(e).__name__}: {e}"
    if not data:
        return None
    for s in data.get("stages", []):
        s["jobs"] = [new_ids.get(str(j), j) for j in s.get("jobs", [])]
        for old, new in new_ids.items():
            if str(s.get("watch_job")) == old:
                s["watch_job"] = new
            if s.get("watch_log"):
                s["watch_log"] = s["watch_log"].replace(f"_{old}.log", f"_{new}.log")
    data["updated"] = _now()
    checkpoint._save(session, data)
    return os.path.join(session, checkpoint.REC_JSON)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("classify", help="is a failed job's log tail (stdin) a node fault?")
    c.add_argument("--state", default="")
    c.add_argument("--exit-codes", default="", help="sacct ExitCode values, comma-separated")
    c.add_argument("--nodes", default="", help="sacct NodeList values, comma-separated")
    r = sub.add_parser("retry", help="resubmit the failed step and its dependants elsewhere")
    r.add_argument("--out", required=True, help="the search folder (submit.sh + jobs.txt)")
    r.add_argument("--job", help="the failed job id (default: the first failed one)")
    r.add_argument("--node", help="a node to exclude besides the one(s) sacct names")
    r.add_argument("--dry-run", action="store_true", help="print the plan, submit nothing")
    r.add_argument("--force", action="store_true",
                   help="retry although nothing marks it as a node fault")
    a = ap.parse_args(argv)
    if a.cmd == "classify":
        split = lambda s: [x for x in s.split(",") if x]  # noqa: E731
        rec = classify(sys.stdin.read(), a.state, split(a.exit_codes), split(a.nodes))
        print(json.dumps({"node_fault": bool(rec), **(rec or {})}, indent=2))
        return 0
    return retry(a.out, a.job, a.node, a.dry_run, a.force)


if __name__ == "__main__":
    sys.exit(main())
