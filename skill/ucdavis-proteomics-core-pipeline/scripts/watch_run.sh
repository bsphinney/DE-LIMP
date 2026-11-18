#!/usr/bin/env bash
# =============================================================================
# watch_run.sh  --  One poll of a running search: report its state, detect known
# errors AND stalls, and suggest the fix. The orchestrator MUST watch every search it
# starts: loop this until the run finishes, and on failure OR stall apply the fix and
# resubmit AUTONOMOUSLY (no user prompt) — never leave a multi-hour search unmonitored.
# Detects both hard failures (sacct state) and STALLS (job RUNNING but log frozen,
# which sacct/squeue cannot see) — see --stall-min and references/watcher.md.
#
# Usage:
#   watch_run.sh --log <logfile>                 # local run: the search log
#   watch_run.sh --slurm <jobid> [--log <file>]  # SLURM job
#   add --hive to query HIVE over SSH (needs HIVE_USER/HIVE_KEY; uses hive_exec.sh)
#
# Emits JSON: {state, done, failed, error_class, fix, log_tail}. Loop pattern:
#   while not done: watch_run.sh ...; sleep 300; done  # then act on failed/fix
#   (--hive: 5 min or more -- each poll is several logins and sacct calls on HIVE's shared
#   login node; a local run: 60 s)
#   With --hive and no usable python3 on this computer (Windows' Store stub), the whole poll
#   runs on HIVE instead, through hive_exec.sh -- the skill must be there (--put-skill).
# =============================================================================
set -uo pipefail
ARGS=("$@")
MODE="local"; JOB=""; LOG=""; HIVE=false; STALL_MIN=15; OUTDIR=""; POLL=0; CHAINDIR=""; SESSION=""
while [ $# -gt 0 ]; do case "$1" in
  --slurm)     MODE="slurm"; JOB="$2"; shift 2;;
  --all)       MODE="chain"; CHAINDIR="$2"; shift 2;;   # watch EVERY job in a chain
  --log)       LOG="$2"; shift 2;;
  --hive)      HIVE=true; shift;;
  --stall-min) STALL_MIN="$2"; shift 2;;   # RUNNING + log frozen this long ⇒ stalled
  --out)       OUTDIR="$2"; shift 2;;      # search output dir ⇒ real progress + narration
  --poll)      POLL="$2"; shift 2;;        # poll counter; rotates the note
  --session)   SESSION="$2"; shift 2;;     # the analysis session: where a NEEDS_DECISION marker is
  *) shift;;
esac; done
HERE="$(cd "$(dirname "$0")" && pwd)"
run() { if $HIVE; then bash "$HERE/hive_exec.sh" "$*"; else bash -c "$*"; fi; }

# With --hive every path here is a HIVE path. Git Bash rewrites an argument that looks like a
# POSIX path (/nfs/... -> C:/Program Files/Git/nfs/...) when it starts a Windows program, and the
# watcher then reported no_jobs_file for a search that had one (staff report, skill 2.13.1).
if $HIVE; then export MSYS_NO_PATHCONV=1 MSYS2_ARG_CONV_EXCL='*'; fi

# A JSON string, without python (the case below has none).
json_str() { local s="${1//\\/\\\\}"; s="${s//\"/\\\"}"; printf '"%s"' "$(printf '%s' "$s" | tr -d '\000-\037')"; }
# No usable python3 here -- on Windows it is often the Microsoft Store stub, which prints "Python
# was not found ..." and exits 9009. Every JSON this script builds needs it, and a healthy search
# read as {failed: true, done: true}, which an agent could act on by resubmitting a running
# search (staff report, skill 2.13.1, Git Bash). With --hive the whole poll runs ON HIVE, whose
# python3 is real; never is a missing parser on this side reported as a failed search.
if [ "$(python3 -c 'print("ok")' 2>/dev/null | tr -d '\r')" != ok ]; then     # \r: Windows' own
  if $HIVE; then
    cmd="bash ~/proteomics-pipeline/scripts/watch_run.sh"
    for a in ${ARGS[@]+"${ARGS[@]}"}; do [ "$a" = --hive ] || cmd="$cmd $(printf '%q' "$a")"; done
    ef="$(mktemp 2>/dev/null || echo "${TMPDIR:-/tmp}/watch_run.$$")"
    out="$(bash "$HERE/hive_exec.sh" "$cmd" 2>"$ef")"; rc=$?
    err="$(cat "$ef" 2>/dev/null)"; rm -f "$ef"
    j="$(printf '%s\n' "$out" | sed -n '/^{/,$p')"           # its JSON, one line or indented
    if [ -n "$j" ]; then printf '%s\n' "$j"; exit "$rc"; fi
    echo "{\"mode\":\"$MODE\",\"state\":\"unknown\",\"failed\":false,\"done\":false,\"error_class\":\"watcher_query_failed\",\"detail\":$(json_str "$(printf '%s\n%s' "$err" "$out" | grep -v '^$' | head -n 2 | tr '\n' ' ')"),\"fix\":\"There is no usable python3 on this computer, so this poll was run on HIVE, and HIVE did not answer it -- this is NOT a job failure. Usually HIVE cannot be reached (HIVE_USER/HIVE_KEY, the network), or the skill is not on HIVE yet (bash scripts/hive_exec.sh --put-skill). Re-run the watcher before drawing any conclusion about the run.\"}"
    exit 3
  fi
  echo "{\"mode\":\"$MODE\",\"state\":\"unknown\",\"failed\":false,\"done\":false,\"error_class\":\"watcher_no_python\",\"fix\":\"There is no usable python3 on this computer (on Windows, python3 is often the Microsoft Store stub), so the watcher cannot read the run's state -- this is NOT a job failure. Install Python 3, or for a search on HIVE add --hive, which runs the poll there.\"}"
  exit 3
fi

# What a Slurm state means: slurm_states.py, the skill's one list (ux H1). EREs, each name also
# as sacct cuts it to 9 characters. BAD = ended without finishing (PREEMPTED and not requeued,
# DEADLINE, BOOT_FAIL, REVOKED, ... and CANCELLED); LIVE = not ended (a REQUEUED job runs again).
ST_OK="$(python3 "$HERE/slurm_states.py" regex ok 2>/dev/null)"
ST_BAD="$(python3 "$HERE/slurm_states.py" regex bad 2>/dev/null)"
ST_REQ="$(python3 "$HERE/slurm_states.py" regex requeued 2>/dev/null)"
if [ -z "$ST_OK" ] || [ -z "$ST_BAD" ] || [ -z "$ST_REQ" ]; then
  echo "{\"mode\":\"$MODE\",\"state\":\"unknown\",\"failed\":false,\"done\":false,\"error_class\":\"watcher_no_states\",\"fix\":\"slurm_states.py next to this script could not be run, so the watcher cannot tell a finished job from a failed one -- this is NOT a job failure. Re-copy the skill's scripts folder and re-run the watcher.\"}"
  exit 3
fi
# 1 when every state in a comma/newline list is terminal (OK or BAD); a state the list does not
# know is never read as finished
all_terminal() { ! printf '%s\n' "$1" | tr ',' '\n' | grep -v '^$' | grep -qvE "($ST_OK|$ST_BAD)"; }

# sacct's JobID|State|ExitCode|NodeList|Elapsed rows for a job (every array task); the one
# query (node_fault.SACCT_FORMAT). Elapsed: a black-hole node's job dies in seconds (signal 53).
NF_FORMAT="JobID,State,ExitCode,NodeList,Elapsed"
# node_fault.py classify, from those rows of a failed job and its FAILED tasks' log tails (the
# one definition of a node fault). Prints its JSON; empty when python3 cannot run it.
node_fault_json() {
  local rows="$1" tail="$2" st codes nodes took
  st="$(printf '%s\n' "$rows" | awk -F'|' 'NF>=4 && $2 !~ /COMPLETED|RUNNING|PENDING/{print $2; exit}')"
  codes="$(printf '%s\n' "$rows" | awk -F'|' 'NF>=4 && $2 !~ /COMPLETED|RUNNING|PENDING/{print $3}' | paste -sd, -)"
  nodes="$(printf '%s\n' "$rows" | awk -F'|' 'NF>=4 && $2 !~ /COMPLETED|RUNNING|PENDING/{print $4}' | sort -u | paste -sd, -)"
  took="$(printf '%s\n' "$rows" | awk -F'|' 'NF>=5 && $2 !~ /COMPLETED|RUNNING|PENDING/{print $5}' | paste -sd, -)"
  printf '%s' "$tail" | python3 "$HERE/node_fault.py" classify --state "$st" --exit-codes "$codes" --nodes "$nodes" --elapsed "$took" 2>/dev/null | tr '\n' ' '
}
# A failed job's log tails for node_fault.py, as `run` reads them. A log NOT READ is never a log
# not written: a signal-53 job is a node fault only when it wrote none (L8R-F5). Each sets
# nf_logs (the tails) and nf_unread (why they could not be read, or empty). A folder that
# exists and holds no matching log is a job that wrote none. The folder's logs/ is read too: the
# Skyline jobs write theirs there (skyline_work/logs/sky_build_<id>.log), and a chain watched with
# --out skyline_work read no log tail for them (pm-215x check d).
read_job_logs() {   # <folder> <glob in it>
  local got
  got="$(run "d=$(printf %q "$1"); if [ -d \"\$d\" ] && [ -x \"\$d\" ]; then echo __LOGDIR_READ__; for f in \"\$d\"/$2 \"\$d\"/logs/$2; do [ -e \"\$f\" ] || continue; if [ -r \"\$f\" ]; then tail -n 100 \"\$f\"; else echo __LOG_UNREADABLE__; fi; done; fi" 2>/dev/null)"
  nf_logs=""; nf_unread=""
  if ! printf '%s\n' "$got" | LC_ALL=C grep -aq '^__LOGDIR_READ__$'; then nf_unread="$1 could not be read"
  elif printf '%s\n' "$got" | LC_ALL=C grep -aq '^__LOG_UNREADABLE__$'; then nf_unread="a log in $1 could not be read"
  else nf_logs="$(printf '%s\n' "$got" | LC_ALL=C sed '/^__LOGDIR_READ__$/d')"; fi
}
read_log_file() {   # <file>
  local got
  got="$(run "f=$(printf %q "$1"); if [ -r \"\$f\" ]; then echo __LOG_READ__; tail -n 100 \"\$f\"; fi" 2>/dev/null)"
  nf_logs=""; nf_unread=""
  if printf '%s\n' "$got" | LC_ALL=C grep -aq '^__LOG_READ__$'; then
    nf_logs="$(printf '%s\n' "$got" | LC_ALL=C sed '/^__LOG_READ__$/d')"
  else nf_unread="--log $1 could not be read"; fi
}
# What to run for one: on HIVE itself, or from a laptop through hive_exec.sh. The folder is
# shell-quoted: service folders hold spaces.
node_fault_fix() {
  local o; o="$(printf %q "$1")"
  echo "Retry it on another node: python3 $HERE/node_fault.py retry --out $o --job $2 (from a laptop: bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/node_fault.py retry --out $o --job $2'). It resubmits that step and the steps waiting on it with --exclude=<node>, at most twice per step (node_faults.json), and refuses when any failed task is not a node fault. Tell the user it was a node problem, not DIA-NN and not their data."
}

# A --dda search whose second pass crashed and fell back to its first-pass report
# (dda_first_pass.py): the caveat, worded by dda_first_pass.caution() HERE, from the search's
# provenance as `run` reads it (with --hive the skill's scripts are not at this path on HIVE).
# Empty for every other search, and for a provenance from before 2.14.2 (no report_used).
dda_caution() {
  local f py
  f="$(mktemp 2>/dev/null || echo "${TMPDIR:-/tmp}/watch_run_dfp.$$")"
  run "cat $(printf %q "$1/search_provenance.json") 2>/dev/null" > "$f" 2>/dev/null
  py="$f"; command -v cygpath >/dev/null 2>&1 && py="$(cygpath -m "$f")"
  [ -s "$f" ] && python3 "$HERE/dda_first_pass.py" caution --provenance "$py" 2>/dev/null | tr -d '\r'
  rm -f "$f"
}

# mtime of a file, in the syntax of whichever kernel `run` will actually execute on.
# GNU stat wants -c %Y, BSD/macOS stat wants -f %m. This CANNOT key off the local
# uname alone: with --hive the command is sent to HIVE, which is Linux no matter what
# the laptop is, so a darwin check here would send BSD syntax to a GNU stat. Getting
# this wrong is silent -- stat writes to stderr, the substitution yields "", and stall
# detection below just never fires.
if $HIVE || [ "$(uname -s)" = "Linux" ]; then
  STAT_MTIME="stat -c %Y"; STAT_SIZE_MTIME="stat -c %s:%Y"
else
  STAT_MTIME="stat -f %m"; STAT_SIZE_MTIME="stat -f %z:%m"
fi

# ---- real progress (outdir_progress: both modes) ----------------------------
# Counted from what the 5-step chain actually leaves on disk, not guessed from the log:
# every finished file drops a .quant, so "34 of 66" is a fact. Runs through run() so it
# works over SSH on HIVE exactly as it does locally.
outdir_progress() {
stage="single"; n_total=0; n_done=0
if [ -n "$OUTDIR" ]; then
  q() { run "$*" 2>/dev/null | tr -dc 0-9; }
  n_total="$(q "wc -l < $(printf %q "$OUTDIR/file_list.txt")")"
  q2="$(q "ls -1 $(printf %q "$OUTDIR/quant_step2")/*.quant 2>/dev/null | wc -l")"
  q4="$(q "ls -1 $(printf %q "$OUTDIR/quant_step4")/*.quant 2>/dev/null | wc -l")"
  has() { [ "$(run "test -s $(printf %q "$1") && echo 1 || echo 0" 2>/dev/null | tr -dc 0-9)" = "1" ]; }
  : "${n_total:=0}"; : "${q2:=0}"; : "${q4:=0}"
  if [ "$n_total" -gt 0 ]; then                       # a 5-step chain lives here
    if   has "$OUTDIR/report.parquet";     then stage="step5"; n_done="$n_total"
    elif has "$OUTDIR/empirical.parquet";  then stage="step4"; n_done="$q4"
    elif [ "$q2" -ge "$n_total" ];         then stage="step3"; n_done="$n_total"
    # step 1b (a chain that has one) writes params.resolved.cfg as its last act: until then it
    # is measuring on representative runs, and no file is being searched yet (S6)
    elif has "$OUTDIR/step1.predicted.speclib" && has "$OUTDIR/step1b_window.sbatch" \
         && ! has "$OUTDIR/params.resolved.cfg"; then stage="step1b"; n_done=0
    elif has "$OUTDIR/step1.predicted.speclib"; then stage="step2"; n_done="$q2"
    else stage="step1"; n_done=0
    fi
  fi
fi
}

state="unknown"; done=false; failed=false; q_failed=false; fix_qf=""; exit_code=""

# --out alone (no --slurm, no --all): a search folder that holds SLURM job ids IS a Slurm search.
# A watcher given only --out used to read no state at all ("search running", failed false) for a
# job sacct reported FAILED -- an agent waiting on it never woke (2.15 acceptance, pep24). The ids
# are jobs.txt's (submit.sh: the whole chain) or, for one job submitted by hand, the ids in its
# log names (<step>_<jobid>.log). Only where sacct exists: a local run has none, and stays local.
if [ "$MODE" = "local" ] && [ -n "$OUTDIR" ]; then
  qo="$(printf %q "$OUTDIR")"
  if [ "$(run "command -v sacct >/dev/null 2>&1 && echo yes" 2>/dev/null | tr -d '\r \n')" = yes ]; then
    if [ -n "$(run "grep -E '[0-9]{5,}' $qo/jobs.txt 2>/dev/null | head -n1" 2>/dev/null | tr -d '\r \n')" ]; then
      MODE="chain"; CHAINDIR="$OUTDIR"
    else
      JOB="$(run "ls -1 $qo 2>/dev/null | sed -nE 's/.*_([0-9]{5,})(_[0-9]+)?\\.(log|out|err)\$/\\1/p' | sort -un | tail -n1" 2>/dev/null | tr -dc 0-9)"
      [ -n "$JOB" ] && MODE="slurm"
    fi
  fi
fi

# --all <dir>: watch the WHOLE chain, not just its last link. Watching only the final
# job is why a step-4 array that timed out went unnoticed for hours -- step 5 simply
# sat PENDING on a dependency that could never be satisfied, which reads as "still
# running". jobs.txt is written by submit.sh with every id in the chain.
if [ "$MODE" = "chain" ]; then
  JT="$CHAINDIR/jobs.txt"; QD="$(printf %q "$CHAINDIR")"   # QD: the folder, quoted for run
  probe="$(run "echo __ok__" 2>&1)"
  if ! printf '%s' "$probe" | grep -q __ok__; then
    # failed/done false: an agent acts on `failed`, and resubmitting a running search on it is
    # the one thing this must never cause (staff report, skill 2.13.1)
    echo "{\"mode\":\"chain\",\"failed\":false,\"done\":false,\"err_class\":\"watcher_query_failed\",\"detail\":$(printf '%s' "$probe" | head -2 | python3 -c 'import sys,json;print(json.dumps(sys.stdin.read()))'),\"fix\":\"Cannot reach the cluster at all — this is NOT a job failure and NOT a missing file. Usually HIVE_USER/HIVE_KEY are unset: save them once to ~/.config/ucdavis-proteomics/hive.env. Re-run the watcher before drawing any conclusion about the run.\"}"
    exit 3
  fi
  # the job ids only: a DIA-NN chain's jobs.txt has one id a line, the Skyline check's has
  # "<label><TAB><id>" (build, import, join, import_retry ...), and a label read as a job id
  # gave "build=UNKNOWN" in the chain (pm-215x check d). The first all-digit word of each line.
  ids="$(run "cat $QD/jobs.txt 2>/dev/null" 2>/dev/null | tr -d '\r' | awk '{for (i = 1; i <= NF; i++) if ($i ~ /^[0-9]+$/) {print $i; break}}' | tr '\n' ' ')"
  if [ -z "$ids" ]; then
    echo "{\"mode\":\"chain\",\"failed\":true,\"done\":true,\"err_class\":\"no_jobs_file\",\"fix\":$(json_str "Cluster is reachable but $JT does not exist — this chain was not submitted via submit.sh (older runs predate jobs.txt). Write it by hand with one job id per line, then re-run.")}"
    exit 2
  fi
  worst=""; cancelled=""; anyfail=false; allterm=true; summary=""; unknown_ids=""; requeued_ids=""
  for id in $ids; do
    js="$(run "sacct -j $id -X --noheader -o State 2>/dev/null | tr -d ' ' | sed 's/+$//' | sort -u | tr '\n' ','" 2>/dev/null)"
    # sacct has no row yet (a job submitted seconds ago): squeue's, as single mode does. Neither
    # answering is not "running": it is a question that went unanswered (ux M4, below).
    [ -z "$js" ] && js="$(run "squeue -j $id -h -o %T 2>/dev/null | tr -d ' ' | sort -u | tr '\n' ','" 2>/dev/null)"
    [ -z "$js" ] && { js="UNKNOWN"; unknown_ids="$unknown_ids $id"; }
    summary="$summary $id=${js%,}"
    # one failure list for both modes (slurm_states.py): a PREEMPTED (not requeued), DEADLINE,
    # BOOT_FAIL or REVOKED step read "done, not failed" here before 2.15 (ux H1)
    printf '%s' "$js" | tr ',' '\n' | grep -E "($ST_BAD)" | grep -qv 'CANCELLED' && { anyfail=true; worst="$id"; }
    # a CANCELLED step is a failure too (a chain read "done, ok" with step 5 cancelled), but a
    # step cancelled BEHIND a failed one is a symptom: the failed one is the one to name
    printf '%s' "$js" | grep -q "CANCELLED" && { anyfail=true; [ -z "$worst" ] && cancelled="$id"; }
    printf '%s' "$js" | tr ',' '\n' | grep -qE "($ST_REQ)" && requeued_ids="$requeued_ids $id"
    all_terminal "$js" || allterm=false
  done
  # Step 1b fell back instead of measuring (probe_fallback.py): the chain runs on, and this is
  # the status view the orchestrator reads, so it says so on every poll from then on.
  # From the search's provenance, as probe_fallback.fallback_modes() reads it: a fallback mode,
  # or -- a provenance from before the modes -- a `probe_fallback` record in it (written only
  # when the search fell back); from probe_fallback.json only when there is no provenance. A
  # stale record from an earlier search must not outvote the provenance (dda-review N1). Bash,
  # not the Python reader, because with --hive this runs on HIVE, where the skill's scripts are
  # not at this path; the provenance is written with indent=2, so these are exact.
  fb=""
  if [ "$(run "p=$QD/search_provenance.json; if [ -f \"\$p\" ]; then grep -qE '\"mode\": \"fallback_|\"probe_fallback\": [{]' \"\$p\" && echo yes; else test -s $QD/probe_fallback.json && echo yes; fi" 2>/dev/null)" = yes ]; then
    fb=",\"probe_fallback\":true,\"caution\":\"CAUTION: step 1b FELL BACK -- what it was to measure was NOT measured: the mass accuracy is the documented DEFAULT level for every run (a search from 2.15.0 or earlier: the scan window too, chosen by DIA-NN per run). Tell the user, with the reason in $CHAINDIR/probe_fallback.json, and re-run the search if the cause was transient.\""
  fi
  [ -z "$worst" ] && worst="${cancelled:-}"
  if $anyfail; then
    # A NODE fault (node_fault.py): the node could not reach the storage -- the job's node check
    # exited 75, SLURM said NODE_FAIL, or the job's log has a storage I/O error. Not the search.
    nf_rows="$(run "sacct -j $worst -X --noheader -P -o $NF_FORMAT 2>/dev/null" 2>/dev/null)"
    # the logs of the tasks that FAILED (an array's first task may have succeeded), at most 3
    nf_tasks="$(printf '%s\n' "$nf_rows" | awk -F'|' 'NF>=4 && $2 !~ /COMPLETED|RUNNING|PENDING/{n=split($1,a,"_"); print (n>1 ? a[2] : "-")}' | head -3 | tr '\n' ' ')"
    nf_tail=""; nf_unread_all=""
    for t in $nf_tasks; do
      if [ "$t" = "-" ]; then pat="*_${worst}.log"; else pat="*_${worst}_${t}.log"; fi
      read_job_logs "$CHAINDIR" "$pat"
      [ -n "$nf_unread" ] && nf_unread_all="$nf_unread"
      nf_tail="$nf_tail
$nf_logs"
    done
    # not read is not "no log": classify() then calls no signal-53 job a black hole
    [ -n "$nf_unread_all" ] && nf_tail="$nf_tail
(the job's log could not be read: $nf_unread_all)"
    nf_json="$(node_fault_json "$nf_rows" "$nf_tail")"
    if printf '%s' "$nf_json" | grep -q '"node_fault": true'; then
      nf_say="$(printf '%s' "$nf_json" | python3 -c 'import sys,json; print(json.load(sys.stdin).get("say",""))' 2>/dev/null | tr -d '"\\')"
      # built by json.dumps: the fix quotes the folder for the shell (printf %q), and a backslash
      # pasted into a hand-built JSON string is not JSON
      NF_FIX="$nf_say $(node_fault_fix "$CHAINDIR" "$worst")" NF_DIR="$CHAINDIR" NF_WORST="$worst" \
      NF_CHAIN="${summary# }" NF_JSON="$nf_json" NF_FB="$fb" python3 -c '
import json, os
o = {"mode": "chain", "dir": os.environ["NF_DIR"], "failed": True, "done": True,
     "error_class": "node_fault", "first_failed_job": os.environ["NF_WORST"],
     "chain": os.environ["NF_CHAIN"]}
fb = os.environ.get("NF_FB") or ""
if fb:
    o.update(json.loads("{" + fb.lstrip(",") + "}"))
o["node_fault"] = json.loads(os.environ["NF_JSON"])
o["fix"] = os.environ["NF_FIX"]
print(json.dumps(o))'

      exit 0
    fi
    lu=""; [ -n "$nf_unread_all" ] && lu=",\"log_unread\":$(json_str "the failed job's log could not be read ($nf_unread_all), so whether a node fault ended it is not known: read it before resubmitting")"
    base="{\"mode\":\"chain\",\"dir\":$(json_str "$CHAINDIR"),\"failed\":true,\"done\":true,\"first_failed_job\":\"$worst\",\"chain\":\"${summary# }\"$fb$lu,\"fix\":\"A step FAILED. Inspect it: watch_run.sh --slurm $worst --hive. Downstream steps will sit PENDING with DependencyNeverSatisfied forever until you fix and resubmit them.\"}"
    # the failed step's own diagnosis (its state, exit code, error class and log tail), by the
    # single-job path -- the same words whether the agent watches the chain or the job
    hv=""; $HIVE && hv="--hive"
    sub="$(bash "$0" --slurm "$worst" --out "$CHAINDIR" $hv 2>/dev/null)"
    BASE="$base" SUB="$sub" python3 -c '
import json, os
o = json.loads(os.environ["BASE"])
try:
    sub = json.loads(os.environ.get("SUB") or "{}")
except ValueError:
    sub = {}
for k in ("exit_code", "log_tail", "error_class", "state"):
    if sub.get(k):
        o[k] = sub[k]
if sub.get("error_class") not in (None, "", "unknown_failure"):
    o["fix"] = (sub.get("fix") or "") + " " + o["fix"]
print(json.dumps(o))'
    exit 0
  fi
  if $allterm; then
    # the two-job route's search is watched here: done, with the caveat when it fell back
    dfp="$(dda_caution "$CHAINDIR")"
    if [ -n "$dfp" ]; then
      DFP="$dfp" CDIR="$CHAINDIR" CHAIN="${summary# }" FB="$fb" python3 -c '
import json, os
o = {"mode": "chain", "dir": os.environ["CDIR"], "failed": False, "done": True,
     "chain": os.environ["CHAIN"]}
fb = os.environ.get("FB") or ""
if fb:
    o.update(json.loads("{" + fb.lstrip(",") + "}"))
o["dda_first_pass"] = True
o["caution"] = " ".join(x for x in (o.get("caution"), os.environ["DFP"]) if x)
print(json.dumps(o))'
      exit 0
    fi
    echo "{\"mode\":\"chain\",\"dir\":$(json_str "$CHAINDIR"),\"failed\":false,\"done\":true,\"chain\":\"${summary# }\"$fb}"
    exit 0
  fi
  # Still running: what the chain is doing, as single mode says it (SKILL.md: relay
  # progress.summary, doing and note on every poll) -- a 66-minute two-job search once had
  # nothing to say but RUNNING. Files done from the output folder (outdir_progress), else the
  # running job's own log: DIA-NN's "File #k/N" line.
  # (a chain of two jobs -- library, then search: the two-job route -- has no step folders)
  OUTDIR="${OUTDIR:-$CHAINDIR}"; stage="single"; n_total=0; n_done=0
  set -- $ids; [ "$#" -gt 2 ] && outdir_progress
  rid=""
  for x in $summary; do case "$x" in *=*RUNNING*) rid="${x%%=*}"; break ;; esac; done
  rlog=""; rfile=""
  if [ -n "$rid" ]; then
    # (an array job's tasks log to <job>_<id>_<task>.log)
    rlog="$(run "for f in $QD/*_$rid.log $QD/*_${rid}_*.log; do [ -f \"\$f\" ] && echo \"\${f##*/}\"; done" 2>/dev/null | head -n1)"
    rfile="$(run "for f in $QD/*_$rid.log $QD/*_${rid}_*.log; do [ -f \"\$f\" ] && tail -n 3000 \"\$f\"; done" 2>/dev/null \
             | grep -aoE 'File #[0-9]+/[0-9]+' | tail -n1)"
  fi
  notes="$(python3 "$HERE/pipeline_notes.py" --stage "$stage" --index "$POLL" 2>/dev/null)"
  CDIR="$CHAINDIR" CHAIN="${summary# }" FB="$fb" STAGE="$stage" NDONE="${n_done:-0}" \
  NTOTAL="${n_total:-0}" RID="$rid" RLOG="$rlog" RFILE="$rfile" NOTES="$notes" \
  UNKNOWN_IDS="${unknown_ids# }" N_IDS="$#" REQUEUED_IDS="${requeued_ids# }" python3 - <<'PY'
import json, os, re
STEP1B = ("step 1b/5 measuring the scan window (and any planned mass accuracy) on representative "
          "runs; no file is searched until it is done")
o = {"mode": "chain", "dir": os.environ["CDIR"], "failed": False, "done": False,
     "chain": os.environ["CHAIN"]}
fb = os.environ.get("FB") or ""
if fb:
    o.update(json.loads("{" + fb.lstrip(",") + "}"))
stage = os.environ.get("STAGE") or "single"
n_done, n_total = int(os.environ.get("NDONE") or 0), int(os.environ.get("NTOTAL") or 0)
step = stage[4:] if stage.startswith("step") else None
where = f"step {step}/5" if step else "search"
p = {"stage": stage, "files_done": n_done, "files_total": n_total}
m = re.match(r"File #(\d+)/(\d+)", os.environ.get("RFILE") or "")
if stage == "step1b":
    p["summary"] = STEP1B
elif n_total and step:
    p["percent"] = round(100.0 * n_done / n_total, 1)
    p["summary"] = f"{where}: {n_done}/{n_total} files done ({p['percent']:.0f}%)"
elif m:
    k, n = int(m.group(1)), int(m.group(2))
    p.update(file=k, files=n, source=os.environ.get("RLOG") or "the running job's log")
    p["summary"] = (f"job {os.environ.get('RID')} ({os.environ.get('RLOG') or 'its log'}): "
                    f"DIA-NN is on file {k} of {n}")
elif os.environ.get("RID"):
    p["summary"] = f"job {os.environ['RID']} running; its log names no file yet"
else:
    p["summary"] = "waiting in the SLURM queue"
o["progress"] = p
try:
    o.update({k: v for k, v in json.loads(os.environ.get("NOTES") or "{}").items()
              if k in ("doing", "why", "note", "note_source")})
except ValueError:
    pass
req = (os.environ.get("REQUEUED_IDS") or "").split()
if req:
    p["requeued_jobs"] = req
    p["summary"] = (f"job(s) {', '.join(req)} REQUEUED by Slurm (usually preempted on the low "
                    f"queue); they will run again by themselves -- " + p["summary"])
# ux M4: a job neither sacct nor squeue knows is a question that went unanswered, never "running"
# -- the same watcher_query_failed single mode gives, with the ids, so it cannot read as progress
unk = (os.environ.get("UNKNOWN_IDS") or "").split()
if unk:
    o["unknown_jobs"] = unk
    o["error_class"] = "watcher_query_failed"
    if len(unk) == int(os.environ.get("N_IDS") or 0):
        o["state"] = "query_failed"
        p["summary"] = "the job states could not be read (sacct and squeue returned nothing)"
    o["fix"] = (f"Could not read the state of job(s) {', '.join(unk)}: sacct and squeue both "
                "returned nothing -- this is NOT a job failure, so do not resubmit on it. "
                "Usually SLURM tools are not on PATH (hive_exec.sh must use a LOGIN shell), or "
                "the job id is wrong / aged out of sacct. Verify: hive_exec.sh 'sacct -j <id> "
                "-X'. Until it is read, the chain is not known to be done.")
    o.setdefault("warnings", []).append(o["fix"])
print(json.dumps(o))
PY
  exit 0
fi
if [ "$MODE" = "slurm" ] && [ -n "$JOB" ]; then
  # -X = one row per JOB STEP, not per .batch/.extern subrecord. For a job ARRAY this
  # is many rows; `head -1` would report one arbitrary task's state as the whole job's
  # -- how a 4-COMPLETED / 14-TIMEOUT array read as healthy. Aggregate instead:
  # failed if ANY task failed, done only when ALL tasks are terminal.
  all_states="$(run "sacct -j $JOB -X --noheader -o State 2>/dev/null | tr -d ' ' | sed 's/+$//' | sort | uniq -c" 2>/dev/null)"
  if [ -z "$all_states" ]; then
    live="$(run "squeue -j $JOB -h -o %T 2>/dev/null | tr -d ' ' | sort -u" 2>/dev/null)"
    if [ -n "$live" ]; then
      all_states="$(printf '   1 %s\n' $live)"
    else
      # Neither sacct nor squeue answered. Do NOT call that PENDING -- that is how a
      # dead job passes for a running one. Say the query failed and stop. Not failed either: an
      # agent acts on `failed`, and resubmitting a running search because a question went
      # unanswered is the one thing this must never cause (review of 2.13.2, round 3). HIVE out
      # of reach is told apart from a job SLURM does not know, as chain mode does.
      probe="$(run "echo __ok__" 2>&1)"
      if ! printf '%s' "$probe" | grep -q __ok__; then
        echo "{\"mode\":\"slurm\",\"job\":$(json_str "$JOB"),\"state\":\"unknown\",\"failed\":false,\"done\":false,\"error_class\":\"watcher_query_failed\",\"detail\":$(json_str "$(printf '%s\n' "$probe" | grep -v '^$' | head -n 2 | tr '\n' ' ')"),\"fix\":\"Cannot reach HIVE at all -- this is NOT a job failure. Usually HIVE_USER/HIVE_KEY are unset (save them once to ~/.config/ucdavis-proteomics/hive.env), or the network is down. Re-run the watcher before drawing any conclusion about the run.\"}"
        exit 3
      fi
      state="query_failed"; done=false; failed=false; q_failed=true
      fix_qf="Could not read job state: sacct and squeue both returned nothing -- this is NOT a job failure, so do not resubmit on it. Usually SLURM tools are not on PATH (hive_exec.sh must use a LOGIN shell), or the job id is wrong / aged out of sacct. Verify: hive_exec.sh 'sacct -j <id> -X'."
    fi
  fi
  if [ "${state:-}" != "query_failed" ]; then
    a_total=$(printf '%s\n' "$all_states" | awk '{s+=$1} END{print s+0}')
    # the one state list (slurm_states.py), as chain mode reads it: terminal = OK or BAD
    a_term=$(printf '%s\n' "$all_states" | awk -v re="($ST_OK|$ST_BAD)" '$2 ~ re {s+=$1} END{print s+0}')
    # BAD but not CANCELLED first: an array with a FAILED task and CANCELLED ones names the FAILED
    a_bad=$(printf '%s\n' "$all_states"  | awk -v re="($ST_BAD)" '$2 ~ re && $2 !~ /CANCELLED/ {s+=$1} END{print s+0}')
    a_run=$(printf '%s\n' "$all_states"  | awk '/RUNNING/{s+=$1} END{print s+0}')
    a_req=$(printf '%s\n' "$all_states"  | awk -v re="($ST_REQ)" '$2 ~ re {s+=$1} END{print s+0}')
    # worst state wins, so a partially-failed array is never reported as healthy
    if   [ "$a_bad"  -gt 0 ]; then state="$(printf '%s\n' "$all_states" | awk -v re="($ST_BAD)" '$2 ~ re && $2 !~ /CANCELLED/ {print $2; exit}')"
    elif [ "$a_run"  -gt 0 ]; then state="RUNNING"
    # all terminal: COMPLETED only when every task completed -- a CANCELLED (or PREEMPTED) job
    # counted as terminal read "COMPLETED", success, for a job that never finished
    elif [ "$a_term" -gt 0 ] && [ "$a_term" -eq "$a_total" ]; then
      state="$(printf '%s\n' "$all_states" | awk '!/COMPLETED/{print $2; exit}')"; : "${state:=COMPLETED}"
    # requeued (preempted, then put back in the queue by Slurm): it runs again by itself -- not
    # done, not failed. Without sacct -D this is the job's newest record (slurm_states.py).
    elif [ "$a_req" -gt 0 ]; then state="REQUEUED"
    else state="PENDING"; fi
    array_summary="$(printf '%s\n' "$all_states" | awk '{printf "%s=%s ", $2, $1}')"
  fi
  st="$state"
  reason="$(run "squeue -j $JOB -h -o %R 2>/dev/null" 2>/dev/null)"
  # The job's exit code, as Slurm recorded it (<code>:<signal>): the first task that did not
  # complete, else the first row. A watcher that names FAILED without it leaves the agent to ask.
  if [ "${state:-}" != "query_failed" ]; then
    ex_rows="$(run "sacct -j $JOB -X --noheader -P -o JobID,State,ExitCode 2>/dev/null" 2>/dev/null)"
    exit_code="$(printf '%s\n' "$ex_rows" | awk -F'|' 'NF>=3 && $2 !~ /COMPLETED|RUNNING|PENDING/{print $3; exit}')"
    [ -z "$exit_code" ] && exit_code="$(printf '%s\n' "$ex_rows" | awk -F'|' 'NF>=3{print $3; exit}')"
  fi
  # done/failed from the one list (slurm_states.py): every BAD state, not four of them -- a
  # PREEMPTED job read "waiting / running" forever here before 2.15 (ux H1)
  if [ "$state" = COMPLETED ]; then done=true
  elif printf '%s' "$state" | grep -qE "($ST_BAD)"; then done=true; failed=true; fi
  # An upstream failure leaves a CHAINED job PENDING FOREVER with this reason. A naive
  # "wait until the job leaves the queue" monitor hangs here indefinitely and never fires
  # — so treat it as a terminal FAILURE that needs recovery. (This is exactly what silently
  # ate the first nail semi-tryptic chain overnight.)
  if printf '%s' "$reason" | grep -qiE "DependencyNeverSatisfied|launch failed"; then
    done=true; failed=true; dep_failed=true
  fi
  # A HELD job never starts on its own either, but nothing is wrong with it: not failed (a
  # failed job gets resubmitted, which would duplicate it), yet the orchestrator must act.
  if printf '%s' "$reason" | grep -qE "JobHeld(User|Admin)"; then held=true; fi
fi

tail_txt=""
[ -n "$LOG" ] && tail_txt="$(run "tail -n 100 $(printf %q "$LOG") 2>/dev/null" 2>/dev/null)"
# No --log: the job's own log in --out or its logs/ (<step>_<jobid>.log, an array task's
# <step>_<jobid>_<n>.log; the Skyline jobs write to logs/), the newest -- the error signatures below and the stack-overflow check read it, and so does the
# user (log_tail). Never guessed when there is none.
if [ -z "$tail_txt" ] && [ -n "$OUTDIR" ] && [ "$MODE" = "slurm" ] && [ -n "$JOB" ] && ! $q_failed; then
  qo="$(printf %q "$OUTDIR")"
  tail_txt="$(run "f=\$(ls -t $qo/*_$JOB.log $qo/*_${JOB}_*.log $qo/logs/*_$JOB.log $qo/logs/*_${JOB}_*.log 2>/dev/null | head -n1); [ -n \"\$f\" ] && tail -n 100 \"\$f\"" 2>/dev/null)"
fi
# the Sage job's own record of how it died (run_search.py's job script, sage_failure.json)
sage_rec=""
[ -n "$OUTDIR" ] && $failed && sage_rec="$(run "cat $(printf %q "$OUTDIR/sage_failure.json") 2>/dev/null" 2>/dev/null)"
# a normalisation check that stopped to ask: its job exits 0 and leaves NEEDS_DECISION.json
# (normalization_check.py --job-exit0). In the session --session names, else the one the log's
# folder (<session>/logs) belongs to. Read only for a job that finished without failing.
# The marker records the job that wrote it (job_id): only THAT job's poll says "a decision is
# owed"; any other job of the session gets a warning and keeps its own summary (ux M8). A marker
# that is there but could not be read (or a read that got no answer) is said, never taken as none.
needs_json=""; needs_unread=""
if $done && ! $failed; then
  nd_sess="$SESSION"
  [ -z "$nd_sess" ] && [ -n "$LOG" ] && [ "$(basename "$(dirname "$LOG")")" = logs ] && nd_sess="$(dirname "$(dirname "$LOG")")"
  if [ -n "$nd_sess" ]; then
    nd_f="$nd_sess/output/norm_check/NEEDS_DECISION.json"
    nd_got="$(run "f=$(printf %q "$nd_f"); if [ ! -e \"\$f\" ]; then echo __ND_NONE__; elif [ -r \"\$f\" ]; then echo __ND_READ__; cat \"\$f\"; else echo __ND_UNREADABLE__; fi" 2>/dev/null)"
    if printf '%s\n' "$nd_got" | grep -q '^__ND_READ__$'; then
      needs_json="$(printf '%s\n' "$nd_got" | sed '/^__ND_READ__$/d')"
    elif ! printf '%s\n' "$nd_got" | grep -q '^__ND_NONE__$'; then
      needs_unread="$nd_f could not be read, so whether a normalisation decision is owed in this session is not known"
    fi
  fi
fi

# A NODE fault first (node_fault.py): a failed job whose node could not reach the storage. Only a
# FAILED job -- an ESTALE that step 1b's probe absorbed sits in a job that succeeded -- and never
# over OOM / TIMEOUT / CANCELLED, which classify() leaves to their own classes below.
nf_json=""
if $failed && [ "$MODE" = "slurm" ] && [ -n "$JOB" ] && ! $q_failed; then
  # The log is evidence either way -- a signal-53 job is a node fault only when it wrote none --
  # so: --log; else (or when --log cannot be read) the job's own log in --out; with neither, or
  # when neither can be read, it is said not to be read (read_job_logs), never "no log".
  nf_tail=""; nf_unread_all=""
  if [ -n "$LOG" ]; then read_log_file "$LOG"; nf_tail="$nf_logs"; nf_unread_all="$nf_unread"; fi
  if [ -z "$LOG" ] || [ -n "$nf_unread_all" ]; then
    if [ -n "$OUTDIR" ]; then
      read_job_logs "$OUTDIR" "*_${JOB}.log"; nf_tail="$nf_logs"
      if [ -n "$nf_unread" ]; then nf_unread_all="${nf_unread_all:+$nf_unread_all; }$nf_unread"
      else nf_unread_all=""; fi
    elif [ -z "$LOG" ]; then nf_unread_all="no --log and no --out"; fi
  fi
  [ -n "$nf_unread_all" ] && nf_tail="(the job's log could not be read: $nf_unread_all)"
  nf_json="$(node_fault_json "$(run "sacct -j $JOB -X --noheader -P -o $NF_FORMAT 2>/dev/null" 2>/dev/null)" "$nf_tail")"
fi

# error signatures -> (class, fix). First match wins.
err_class=""; fix=""
hay="$tail_txt
$state"
m() { printf '%s' "$hay" | grep -qiE "$1"; }
if printf '%s' "$nf_json" | grep -q '"node_fault": true'; then
  err_class="node_fault"
  fix="$(printf '%s' "$nf_json" | python3 -c 'import sys,json; print(json.load(sys.stdin).get("say",""))' 2>/dev/null) $(node_fault_fix "${OUTDIR:-<search out dir>}" "$JOB")"
elif $failed && { printf '%s' "$sage_rec" | grep -q 'sage_stack_overflow' || m "fatal runtime error: stack overflow|has overflowed its stack"; };
                                                                                     then err_class="sage_stack_overflow"
  fix="Sage died with a stack overflow (exit 134): an intermittent Rust crash, NOT out of memory and NOT your data. Re-run the same job once, unchanged -- jobs generated by this skill version export RUST_MIN_STACK (64 MiB) before Sage starts, so regenerate it with run_search.py --sbatch if this one was made earlier. Do not change the search settings for it; if it crashes again with RUST_MIN_STACK set, report it."
  h="$(printf '%s' "$sage_rec" | python3 -c 'import sys,json; print(json.load(sys.stdin).get("hint",""))' 2>/dev/null)"
  [ -n "$h" ] && fix="$h $fix"
elif $failed && printf '%s' "$state" | grep -q '^PREEMPTED'; then err_class="preempted"
  fix="Slurm PREEMPTED this job and did not requeue it: a higher-priority job took its node (the preemptible queue, publicgrp/low). It is NOT a fault in the search or the data. Resubmit the step -- on the high queue (not --preemptible, not partition low) so it cannot be preempted again. In a chain, the steps waiting on it never start: resubmit them too, reusing the outputs already finished."
elif $failed && printf '%s' "$state" | grep -qE '^(DEADLINE|BOOT_FAIL|REVOKED|SPECIAL_E)'; then err_class="slurm_ended"
  fix="Slurm ended the job ($state) before it finished: DEADLINE = its --deadline passed; BOOT_FAIL = the node failed to start; REVOKED = a federated copy ran elsewhere; SPECIAL_EXIT = requeued and held. It is NOT a fault in the search or the data. Check it (hive_exec.sh 'sacct -j $JOB -X'), then resubmit the step."
elif m "out.of.memory|oom-kill|OUT_OF_MEMORY|std::bad_alloc|cannot allocate";       then err_class="out_of_memory"; fix="Raise the sbatch --mem (e.g. 64G→128G) and resubmit; for DIA-NN try fewer threads or --min-corr.";
elif m "TIMEOUT|DUE TO TIME LIMIT|CANCELLED.*TIME";                                  then err_class="timeout";       fix="Raise --time in the sbatch (or split the run) and resubmit.";
# CANCELLED (not for time, not for memory: both above): someone ran scancel, or Slurm cancelled it
# because a job it waited on failed. The search's log is not the cause, so no log signature below
# may name one -- a cancelled FragPipe job read "empty_results" from its own counts (pm-215x d).
elif $failed && printf '%s' "$state" | grep -q '^CANCELLED';                        then err_class="cancelled";
  fix="The job was CANCELLED: a person (or a script) ran scancel, or Slurm cancelled it because a job it waited on failed. It is not a fault the search's log shows. Find out which before resubmitting: sacct -j $JOB -X -o JobID,State,Reason (\"CANCELLED by <uid>\" names who). If it was stopped on purpose, nothing is wrong; if a step before it failed, fix that step first."
elif m "dotnet: not found|dotnet: command not found";                               then err_class="diann_no_dotnet";fix="Wrong DIA-NN container (no .NET → .raw silently skipped). Use the HIVE native build (build_<v>/diann-<v>/diann-linux) or a .NET-enabled image.";
elif m "Number of IDs at 0.01 FDR: 0";                                              then err_class="diann_zero_ids";  fix="DIA-NN completed but identified nothing - a SILENT null result, not a crash. PRESERVE report.log.txt before re-running; it is the only evidence. Check in order: (1) the library actually generated precursors, (2) the mzML really are DIA with isolation windows (detect_acquisition.py), (3) FASTA matches the organism, (4) mass accuracy. NOTE: DIA-NN warning about generating the predicted library in a separate step is BENIGN per its author - do not chase it.";
elif m "java.lang.OutOfMemoryError|Java heap space|Answer from Java side is empty";  then err_class="spark_heap";      fix="A JVM heap OOM inside Fulcrum/Spark. RAISING --mem DOES NOT HELP: Spark sizes its driver heap independently of the SLURM allocation. Pass spark_config = {\"spark.driver.memory\" = \"32g\"} in the workflow TOML (Fulcrum forwards it to SparkSession.builder.config).";
# Only for a job that failed, and only the engines' own words: "empty" matched ordinary log text,
# and "0 proteins" matched FragPipe's count line "peptides=1340 proteins=443" (pm-215x check d).
elif $failed && m "(^|[^0-9.,=])0 proteins|No fragment ions|No precursors|no spectra";    then err_class="empty_results"; fix="Check the FASTA matches the organism, the mass-accuracy setting, and that the raw files are the expected acquisition type.";
elif m "CUDA|no kernel image|cuDNN|device-side|GPU.*not";                           then err_class="gpu";           fix="AlphaDIA needs a GPU. Submit to a GPU node (sbatch --gres=gpu:1) or reduce batch size.";
elif m "msconvert.*not found|requires mzML|no mzML|does not match its pinned sha256"; then err_class="sage_no_mzml";  fix="Sage/Radiant read mzML and the conversion in the job failed (run_search.py converts .raw with ThermoRawFileParser and, to quantify, a timsTOF .d with msconvert in the Core's pinned ProteoWizard container -- apptainer + wine, its sha256 checked first -- before the search). Read the converter's own error above the FAILED line: for ThermoRawFileParser usually .NET (bash scripts/ensure_dotnet8.sh); for a .d, a container that does not match its pinned sha256 is not the image the skill was tested with (re-copy it, or pass another with --pwiz-sif), otherwise msconvert's own message; for .raw via msconvert a Linux build without vendor readers. A .d searched for identifications only needs no conversion (estimate_params.py --ids-only). Fix, then regenerate the job with run_search.py --sbatch and resubmit.";
elif m "Disk quota exceeded|No space left";                                         then err_class="disk";          fix="Out of disk/quota. Free space or point --out elsewhere and resubmit.";
elif m "No such file|cannot open|does not exist|not found.*(fasta|\\.d|\\.raw)";     then err_class="missing_input"; fix="An input path is wrong (fasta/raw). Re-check paths (Windows→WSL/HIVE translation) and resubmit.";
elif [ "${dep_failed:-false}" = true ];                                              then err_class="dependency_failed"; fix="An UPSTREAM job in the chain failed, so this one is stuck PENDING with DependencyNeverSatisfied (it will NEVER run and never leave the queue). Find the failed step (sacct -j <arrayjob>), apply that step's fix, and resubmit the downstream steps reusing already-computed outputs (.quant, step1.predicted.speclib) — don't restart the whole chain.";
elif [ "${held:-false}" = true ];                                                   then err_class="held"; fix="The job is HELD (${reason}) and will never start by itself -- do NOT resubmit it (that makes a duplicate). JobHeldUser: release it with scontrol release <jobid> once whatever it was held for is resolved. JobHeldAdmin: ask the HIVE admins why.";
elif $failed;                                                                        then err_class="unknown_failure"; fix="Read the full log; diagnose via references/watcher.md; fix and resubmit.";
fi
# ...and when the log the node-fault check needed could not be read, say so (read_job_logs)
if [ "$err_class" = unknown_failure ] && [ -n "${nf_unread_all:-}" ]; then
  fix="$fix The job's log could not be read here (${nf_unread_all}), so whether a node fault ended it is not known: read it before resubmitting."
fi

# STALL detection: job RUNNING but its log has not advanced in STALL_MIN minutes.
# A hung file (e.g. a pathological DIA-NN run) keeps the job in RUNNING while the log
# freezes — sacct/squeue will NOT flag it. This is the "NA41-class" failure.
stalled=false
if [ -z "$err_class" ] && [ -n "$LOG" ] && printf '%s' "$state" | grep -qiE "RUNNING|^R$"; then
  now="$(run "date +%s" 2>/dev/null)"
  mt="$(run "$STAT_MTIME $(printf %q "$LOG") 2>/dev/null" 2>/dev/null | tr -dc 0-9)"
  if [ -n "$now" ] && [ -n "$mt" ]; then
    age=$(( now - mt ))
    if [ "$age" -ge $(( STALL_MIN * 60 )) ]; then
      stalled=true; err_class="stalled"
      fix="RUNNING but log frozen ${age}s (> ${STALL_MIN} min) — likely a hung file (pathological DIA-NN run). Auto-recover: scancel this task/job, retry it ONCE on a fresh node; if it stalls again, DROP that file and continue (the 5-step chain's step 4 auto-skips a file with no .quant), then resubmit downstream steps reusing the completed .quant. Log the dropped file in Data Quality Notes."
    fi
  fi
fi

# ---- real progress ----------------------------------------------------------
outdir_progress
# ---- Sage LFQ: a finished search can still have unusable quantities -------------------
# sage_lfq_check.py (run in the Sage job and by run_search.py --adapt-only) records whether the
# LFQ window fits the runs' MS1 mass error. A WARNING there is not a failure -- the job succeeded
# and the identifications stand -- so it is reported beside the state, never as error_class.
# Its file has an engine-neutral name and 2.14.x's (sage_lfq_check.RECORD_NAMES, the newer first;
# a test ties these to it).
lfq_rec=""
[ -n "$OUTDIR" ] && lfq_rec="$(run "cat $(printf %q "$OUTDIR/quantification_check.json") 2>/dev/null || cat $(printf %q "$OUTDIR/sage_lfq_check.json") 2>/dev/null" 2>/dev/null)"
# ---- DIA-NN: the sample proteins the contaminant list took (contaminant_competition.py, the
# search's last job). Not a failure either. Read only once the search has finished, and shown
# only when it describes the report as it is now (size and mtime, as the record states them): a
# re-search into the same folder must never show the previous search's genes (review of 2.14).
# Into a file, not a variable: a 400-run record is larger than one environment variable may be.
cc_file=""; cc_py=""; cc_stat=""
if [ -n "$OUTDIR" ] && [ "$done" = true ] && [ "$failed" = false ]; then
  cc_file="$(mktemp 2>/dev/null || echo "${TMPDIR:-/tmp}/watch_run_cc.$$")"
  run "cat $(printf %q "$OUTDIR/contaminant_competition.json") 2>/dev/null" > "$cc_file" 2>/dev/null
  # Git Bash: a native Windows python3 reads C:/..., not /tmp/...
  cc_py="$cc_file"; command -v cygpath >/dev/null 2>&1 && cc_py="$(cygpath -m "$cc_file")"
  cc_rep="$(python3 -c 'import json, sys; print(json.load(open(sys.argv[1])).get("report") or "")' \
            "$cc_py" 2>/dev/null | tr -d '\r')"
  [ -n "$cc_rep" ] && cc_stat="$(run "$STAT_SIZE_MTIME $(printf %q "$cc_rep") 2>/dev/null" \
                                 2>/dev/null | tr -dc '0-9:')"
fi
# ---- DIA-NN --dda: a crashed second pass fell back to the first-pass report (dda_caution). Done,
# not failed: the report is complete, but the user must hear the caveat.
dfp=""
if [ -n "$OUTDIR" ] && [ "$done" = true ] && [ "$failed" = false ]; then
  dfp="$(dda_caution "$OUTDIR")"
fi

# ---- narration: what this stage is doing, plus something to read while waiting
if $q_failed; then err_class="watcher_query_failed"; fix="$fix_qf"; fi

notes="$(python3 "$HERE/pipeline_notes.py" --stage "$stage" --index "$POLL" 2>/dev/null)"

EXITC="$exit_code" NEEDS_JSON="$needs_json" NEEDS_UNREAD="$needs_unread" \
STATE="$state" DONE="$done" FAILED="$failed" STALLED="$stalled" JOB="$JOB" MODE="$MODE" \
ECLASS="$err_class" FIX="$fix" TAIL="$tail_txt" ATASKS="${array_summary:-}" \
STAGE="$stage" NDONE="${n_done:-0}" NTOTAL="${n_total:-0}" NOTES="$notes" \
REASON="${reason:-}" ATERM="${a_term:-0}" LFQ_REC="$lfq_rec" NF_JSON="$nf_json" CC_FILE="$cc_py" \
CC_STAT="$cc_stat" DFP="$dfp" \
python3 - <<'PY'
import os, json, re
n_done, n_total = int(os.environ.get("NDONE") or 0), int(os.environ.get("NTOTAL") or 0)
stage = os.environ.get("STAGE", "single")
out = {
    "mode": os.environ["MODE"], "job": os.environ["JOB"], "state": os.environ["STATE"],
    "done": os.environ["DONE"] == "true", "failed": os.environ["FAILED"] == "true",
    "stalled": os.environ.get("STALLED") == "true",
    "error_class": os.environ["ECLASS"], "fix": os.environ["FIX"],
    **({"array_tasks": os.environ["ATASKS"].strip()} if os.environ.get("ATASKS","").strip() else {}),
}
step_no = stage[4:] if stage.startswith("step") else None
progress = {"stage": stage, "files_done": n_done, "files_total": n_total}
if step_no:
    progress["step"] = f"{step_no}/5"
if n_total:
    progress["percent"] = round(100.0 * n_done / n_total, 1)
# One sentence the orchestrator can hand the user verbatim.
where = f"step {step_no}/5" if step_no else "search"
progress["summary"] = (f"{where}: {n_done}/{n_total} files done"
                       + (f" ({progress['percent']:.0f}%)" if n_total else "")) \
    if n_total else f"{where} running"
if stage == "step1b":
    progress["summary"] = ("step 1b/5 measuring the scan window (and any planned mass accuracy) "
                           "on representative runs; no file is searched until it is done")
# A job that has not started is not searching anything. A held search (JobHeldUser) was
# reported as "search running" / "Searching your files" (HIVE e2e test 2026-09-23). But PENDING
# alone does not mean nothing ran: the chain's step-5 job sits PENDING (Dependency) while
# steps 2-4 do the work, and an array can be half done with the rest queued -- those keep
# their file count (review 2026-09-23).
reason = os.environ.get("REASON", "").strip()
queued = (out["state"] == "PENDING" and not out["done"] and n_done == 0
          and int(os.environ.get("ATERM") or 0) == 0)
if out["state"] == "PENDING" and not out["done"]:
    if "JobHeld" in reason:
        progress["summary"] = f"{where}: HELD ({reason}) -- will not start until released"
    elif queued:
        progress["summary"] = (f"{where}: waiting for earlier steps ({reason})"
                               if reason.startswith("Dependency")
                               else f"{where}: queued, not started yet ({reason or 'PENDING'})")
    elif not n_total:
        progress["summary"] = f"{where}: partly done, remaining work queued ({reason or 'PENDING'})"
# A job that has ended is not "running": its summary says how it ended, with Slurm's exit code
# (2.15 acceptance: a FAILED job read "search running"). A finished job that left a question
# (NEEDS_DECISION.json: a normalisation check) says so -- a question, not a failure.
code = (os.environ.get("EXITC") or "").strip()
if out["failed"] and out["state"] != "query_failed":
    progress["summary"] = (f"{where}: job {out['state']}" + (f" (exit code {code})" if code else "")
                           + " -- it did NOT finish")
    out["doing"] = f"Stopped: the job ended {out['state']}; nothing is being searched now."
elif out["done"] and out["state"] == "COMPLETED" and not n_total:
    progress["summary"] = f"{where}: job finished (COMPLETED)"
elif out["state"] == "REQUEUED":
    progress["summary"] = (f"{where}: REQUEUED by Slurm (usually preempted on the low queue) -- "
                           "it will run again by itself; not finished, not failed")
    out["doing"] = "Requeued -- waiting to run again; nothing is being searched right now."
if code:
    out["exit_code"] = code
out["progress"] = progress
try:
    nd = json.loads(os.environ.get("NEEDS_JSON") or "null")
except ValueError:
    nd = None
nd_job = str(nd.get("job_id") or "") if isinstance(nd, dict) else ""
if isinstance(nd, dict) and nd.get("needs_decision") and nd_job != out["job"]:
    # the marker belongs to another job (or names none): this job's summary is its own (ux M8)
    out.setdefault("warnings", []).append(
        f"A normalisation decision is owed in this session (NEEDS_DECISION.json, written by "
        f"{'job ' + nd_job if nd_job else 'a run with no Slurm job id'}, not job {out['job']}): "
        f"show {nd.get('summary') or 'NORMALIZATION_CHECK.md'} to the user and record their "
        "answer with normalization_check.py decide.")
    nd = None
if os.environ.get("NEEDS_UNREAD"):
    out.setdefault("warnings", []).append(os.environ["NEEDS_UNREAD"])
if isinstance(nd, dict) and nd.get("needs_decision"):
    out["needs_decision"] = True
    shown = nd.get("summary") or "NORMALIZATION_CHECK.md"
    progress["summary"] = (f"{where}: job finished (COMPLETED); a normalisation DECISION is owed "
                           "-- the check stopped to ask, it did not fail")
    out["doing"] = "Waiting for the user's normalisation decision."
    out["fix"] = (f"The normalisation check stopped to ask ({nd.get('why') or 'it needs a decision'}"
                  f"): this is a question, not a failure -- do not resubmit. Show {shown} to the "
                  "user, ask, and record THEIR answer with normalization_check.py decide.")
    out.setdefault("warnings", []).append(out["fix"])
try:                                     # node_fault.py's verdict, when it is a node fault
    nf = json.loads(os.environ.get("NF_JSON") or "null")
    if isinstance(nf, dict) and nf.get("node_fault"):
        out["node_fault"] = nf
except ValueError:
    pass
try:                                     # never over a "doing" said above (a FAILED job
    out.update({k: v for k, v in json.loads(os.environ.get("NOTES") or "{}").items()  # read
                if k in ("doing", "why", "note", "note_source")              # "Searching")
                and not (k == "doing" and "doing" in out)})
except Exception:
    pass
if queued:
    out["doing"] = ("Held -- the job will not start until it is released." if "JobHeld" in reason
                    else "Waiting in the SLURM queue -- the job has not started, so nothing is "
                         "being searched yet.")
# Sage LFQ (sage_lfq_check.py): its record with --out, else what the log tail says. The MS1-peak
# count is Sage's own "discovered N target MS1 peaks at 5% FDR"; 0 means no usable quantities.
tail = os.environ.get("TAIL", "")
try:
    lfq = json.loads(os.environ.get("LFQ_REC") or "null")
except ValueError:
    lfq = None
peaks = re.findall(r"discovered (\d+) target MS1 peaks at 5% FDR", tail)
if isinstance(lfq, dict) and lfq.get("status"):
    out["sage_lfq"] = {k: lfq.get(k) for k in ("status", "target_ms1_peaks_5pct_fdr",
                                               "ppm_tolerance", "runs_outside_lfq_window",
                                               "suggested_ppm_tolerance") if k in lfq}
    if lfq["status"] in ("warn", "unchecked") and lfq.get("message"):
        out.setdefault("warnings", []).append(lfq["message"])
elif peaks:
    out["sage_lfq"] = {"target_ms1_peaks_5pct_fdr": int(peaks[-1]), "source": "log tail"}
if not (isinstance(lfq, dict) and lfq.get("status")):
    for ln in tail.splitlines():
        if "[sage_lfq_check] WARNING:" in ln:
            out.setdefault("warnings", []).append(ln.split("WARNING:", 1)[1].strip())
try:
    with open(os.environ.get("CC_FILE") or "") as fh:
        cc = json.load(fh)
except (OSError, ValueError):
    cc = None
st = cc.get("report_stat") if isinstance(cc, dict) else None
now = os.environ.get("CC_STAT") or ""
fresh = isinstance(st, dict) and f"{st.get('size')}:{st.get('mtime')}" == now
if isinstance(cc, dict) and cc.get("status") and not now:
    out["contaminant_competition"] = {
        "status": "unchecked", "reason": "the report the record names is not there or could not "
                                         "be read, so whether the record describes it is not "
                                         "known (poll again if the connection failed)"}
elif isinstance(cc, dict) and cc.get("status") and not fresh:
    out["contaminant_competition"] = {
        "status": "stale", "reason": "the record does not describe the report as it is now (its "
                                     "size or modification time differ), so it was not read"}
elif isinstance(cc, dict) and cc.get("status"):
    out["contaminant_competition"] = {
        "status": cc["status"], "reason": cc.get("reason"),
        "competed_out": [x.get("gene") for x in cc.get("competed_out") or []],
        "siphoned": [x.get("gene") for x in cc.get("siphoned") or []],
        "drop_contaminant": cc.get("drop_contaminant") or []}
    if (cc.get("competed_out") or cc.get("siphoned") or cc["status"] == "not_run") \
            and cc.get("message"):
        out.setdefault("warnings", []).append(cc["message"])
if os.environ.get("DFP"):
    out["dda_first_pass"] = True
    out["caution"] = os.environ["DFP"]
    out.setdefault("warnings", []).append(os.environ["DFP"])
out["log_tail"] = tail[-1500:]
print(json.dumps(out, indent=2))
PY
rc=$?                                   # the poll's own status, as before the record was read
[ -z "$cc_file" ] || rm -f "$cc_file"
exit "$rc"
