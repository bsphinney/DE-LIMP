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
#   while not done: watch_run.sh ...; sleep 60; done   # then act on failed/fix
# =============================================================================
set -uo pipefail
MODE="local"; JOB=""; LOG=""; HIVE=false; STALL_MIN=15; OUTDIR=""; POLL=0; CHAINDIR=""
while [ $# -gt 0 ]; do case "$1" in
  --slurm)     MODE="slurm"; JOB="$2"; shift 2;;
  --all)       MODE="chain"; CHAINDIR="$2"; shift 2;;   # watch EVERY job in a chain
  --log)       LOG="$2"; shift 2;;
  --hive)      HIVE=true; shift;;
  --stall-min) STALL_MIN="$2"; shift 2;;   # RUNNING + log frozen this long ⇒ stalled
  --out)       OUTDIR="$2"; shift 2;;      # search output dir ⇒ real progress + narration
  --poll)      POLL="$2"; shift 2;;        # poll counter; rotates the note
  *) shift;;
esac; done
HERE="$(cd "$(dirname "$0")" && pwd)"
run() { if $HIVE; then bash "$HERE/hive_exec.sh" "$*"; else bash -c "$*"; fi; }

# sacct's JobID|State|ExitCode|NodeList rows for a job (every array task); the one query.
NF_FORMAT="JobID,State,ExitCode,NodeList"
# node_fault.py classify, from those rows of a failed job and its FAILED tasks' log tails (the
# one definition of a node fault). Prints its JSON; empty when python3 cannot run it.
node_fault_json() {
  local rows="$1" tail="$2" st codes nodes
  st="$(printf '%s\n' "$rows" | awk -F'|' 'NF>=4 && $2 !~ /COMPLETED|RUNNING|PENDING/{print $2; exit}')"
  codes="$(printf '%s\n' "$rows" | awk -F'|' 'NF>=4 && $2 !~ /COMPLETED|RUNNING|PENDING/{print $3}' | paste -sd, -)"
  nodes="$(printf '%s\n' "$rows" | awk -F'|' 'NF>=4 && $2 !~ /COMPLETED|RUNNING|PENDING/{print $4}' | sort -u | paste -sd, -)"
  printf '%s' "$tail" | python3 "$HERE/node_fault.py" classify --state "$st" --exit-codes "$codes" --nodes "$nodes" 2>/dev/null | tr '\n' ' '
}
# What to run for one: on HIVE itself, or from a laptop through hive_exec.sh. The folder is
# shell-quoted: service folders hold spaces.
node_fault_fix() {
  local o; o="$(printf %q "$1")"
  echo "Retry it on another node: python3 $HERE/node_fault.py retry --out $o --job $2 (from a laptop: bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/node_fault.py retry --out $o --job $2'). It resubmits that step and the steps waiting on it with --exclude=<node>, at most twice per step (node_faults.json), and refuses when any failed task is not a node fault. Tell the user it was a node problem, not DIA-NN and not their data."
}

# mtime of a file, in the syntax of whichever kernel `run` will actually execute on.
# GNU stat wants -c %Y, BSD/macOS stat wants -f %m. This CANNOT key off the local
# uname alone: with --hive the command is sent to HIVE, which is Linux no matter what
# the laptop is, so a darwin check here would send BSD syntax to a GNU stat. Getting
# this wrong is silent -- stat writes to stderr, the substitution yields "", and stall
# detection below just never fires.
if $HIVE || [ "$(uname -s)" = "Linux" ]; then STAT_MTIME="stat -c %Y"; else STAT_MTIME="stat -f %m"; fi

state="unknown"; done=false; failed=false; q_failed=false; fix_qf=""

# --all <dir>: watch the WHOLE chain, not just its last link. Watching only the final
# job is why a step-4 array that timed out went unnoticed for hours -- step 5 simply
# sat PENDING on a dependency that could never be satisfied, which reads as "still
# running". jobs.txt is written by submit.sh with every id in the chain.
if [ "$MODE" = "chain" ]; then
  JT="$CHAINDIR/jobs.txt"; QD="$(printf %q "$CHAINDIR")"   # QD: the folder, quoted for run
  probe="$(run "echo __ok__" 2>&1)"
  if ! printf '%s' "$probe" | grep -q __ok__; then
    echo "{\"mode\":\"chain\",\"failed\":true,\"done\":true,\"err_class\":\"watcher_query_failed\",\"detail\":$(printf '%s' "$probe" | head -2 | python3 -c 'import sys,json;print(json.dumps(sys.stdin.read()))'),\"fix\":\"Cannot reach the cluster at all — this is NOT a job failure and NOT a missing file. Usually HIVE_USER/HIVE_KEY are unset: save them once to ~/.config/ucdavis-proteomics/hive.env. Re-run the watcher before drawing any conclusion about the run.\"}"
    exit 3
  fi
  ids="$(run "cat $QD/jobs.txt 2>/dev/null | tr '\n' ' '" 2>/dev/null)"
  if [ -z "$ids" ]; then
    echo "{\"mode\":\"chain\",\"failed\":true,\"done\":true,\"err_class\":\"no_jobs_file\",\"fix\":\"Cluster is reachable but $JT does not exist — this chain was not submitted via submit.sh (older runs predate jobs.txt). Write it by hand with one job id per line, then re-run.\"}"
    exit 2
  fi
  worst=""; anyfail=false; allterm=true; summary=""
  for id in $ids; do
    js="$(run "sacct -j $id -X --noheader -o State 2>/dev/null | tr -d ' ' | sed 's/+$//' | sort -u | tr '\n' ','" 2>/dev/null)"
    [ -z "$js" ] && js="UNKNOWN"
    summary="$summary $id=${js%,}"
    printf '%s' "$js" | grep -qE "FAILED|TIMEOUT|OUT_OF_ME|NODE_FAIL" && { anyfail=true; worst="$id"; }
    printf '%s' "$js" | grep -qE "RUNNING|PENDING|UNKNOWN" && allterm=false
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
    fb=",\"probe_fallback\":true,\"caution\":\"CAUTION: step 1b FELL BACK -- the scan window (and any planned mass accuracy) was NOT measured; DIA-NN chooses it per run. Tell the user, with the reason in $CHAINDIR/probe_fallback.json, and re-run the search if the cause was transient.\""
  fi
  if $anyfail; then
    # A NODE fault (node_fault.py): the node could not reach the storage -- the job's node check
    # exited 75, SLURM said NODE_FAIL, or the job's log has a storage I/O error. Not the search.
    nf_rows="$(run "sacct -j $worst -X --noheader -P -o $NF_FORMAT 2>/dev/null" 2>/dev/null)"
    # the logs of the tasks that FAILED (an array's first task may have succeeded), at most 3
    nf_tasks="$(printf '%s\n' "$nf_rows" | awk -F'|' 'NF>=4 && $2 !~ /COMPLETED|RUNNING|PENDING/{n=split($1,a,"_"); print (n>1 ? a[2] : "-")}' | head -3 | tr '\n' ' ')"
    nf_tail=""
    for t in $nf_tasks; do
      if [ "$t" = "-" ]; then pat="*_${worst}.log"; else pat="*_${worst}_${t}.log"; fi
      nf_tail="$nf_tail
$(run "for f in $QD/$pat; do [ -f \"\$f\" ] && tail -n 100 \"\$f\"; done" 2>/dev/null)"
    done
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
    echo "{\"mode\":\"chain\",\"dir\":\"$CHAINDIR\",\"failed\":true,\"done\":true,\"first_failed_job\":\"$worst\",\"chain\":\"${summary# }\"$fb,\"fix\":\"A step FAILED. Inspect it: watch_run.sh --slurm $worst --hive. Downstream steps will sit PENDING with DependencyNeverSatisfied forever until you fix and resubmit them.\"}"
    exit 0
  fi
  if $allterm; then
    echo "{\"mode\":\"chain\",\"dir\":\"$CHAINDIR\",\"failed\":false,\"done\":true,\"chain\":\"${summary# }\"$fb}"
  else
    echo "{\"mode\":\"chain\",\"dir\":\"$CHAINDIR\",\"failed\":false,\"done\":false,\"chain\":\"${summary# }\"$fb}"
  fi
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
      # dead job passes for a running one. Say the query failed and stop.
      state="query_failed"; done=false; failed=true; q_failed=true
      fix_qf="Could not read job state: sacct and squeue both returned nothing. Usually SLURM tools are not on PATH (hive_exec.sh must use a LOGIN shell), or the job id is wrong / aged out of sacct. Verify: hive_exec.sh 'sacct -j <id> -X'."
    fi
  fi
  if [ "${state:-}" != "query_failed" ]; then
    a_total=$(printf '%s\n' "$all_states" | awk '{s+=$1} END{print s+0}')
    a_term=$(printf '%s\n' "$all_states" | awk '/COMPLETED|FAILED|TIMEOUT|OUT_OF_ME|CANCELLED|NODE_FAIL|PREEMPTED/{s+=$1} END{print s+0}')
    a_bad=$(printf '%s\n' "$all_states"  | awk '/FAILED|TIMEOUT|OUT_OF_ME|NODE_FAIL/{s+=$1} END{print s+0}')
    a_run=$(printf '%s\n' "$all_states"  | awk '/RUNNING/{s+=$1} END{print s+0}')
    # worst state wins, so a partially-failed array is never reported as healthy
    if   [ "$a_bad"  -gt 0 ]; then state="$(printf '%s\n' "$all_states" | awk '/FAILED|TIMEOUT|OUT_OF_ME|NODE_FAIL/{print $2; exit}')"
    elif [ "$a_run"  -gt 0 ]; then state="RUNNING"
    elif [ "$a_term" -gt 0 ] && [ "$a_term" -eq "$a_total" ]; then state="COMPLETED"
    else state="PENDING"; fi
    array_summary="$(printf '%s\n' "$all_states" | awk '{printf "%s=%s ", $2, $1}')"
  fi
  st="$state"
  reason="$(run "squeue -j $JOB -h -o %R 2>/dev/null" 2>/dev/null)"
  case "$state" in
    COMPLETED)                                   done=true;;
    FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)       done=true; failed=true;;
    CANCELLED*)                                   done=true; failed=true;;
  esac
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

# A NODE fault first (node_fault.py): a failed job whose node could not reach the storage. Only a
# FAILED job -- an ESTALE that step 1b's probe absorbed sits in a job that succeeded -- and never
# over OOM / TIMEOUT / CANCELLED, which classify() leaves to their own classes below.
nf_json=""
if $failed && [ "$MODE" = "slurm" ] && [ -n "$JOB" ] && ! $q_failed; then
  nf_json="$(node_fault_json "$(run "sacct -j $JOB -X --noheader -P -o $NF_FORMAT 2>/dev/null" 2>/dev/null)" "$tail_txt")"
fi

# error signatures -> (class, fix). First match wins.
err_class=""; fix=""
hay="$tail_txt
$state"
m() { printf '%s' "$hay" | grep -qiE "$1"; }
if printf '%s' "$nf_json" | grep -q '"node_fault": true'; then
  err_class="node_fault"
  fix="$(printf '%s' "$nf_json" | python3 -c 'import sys,json; print(json.load(sys.stdin).get("say",""))' 2>/dev/null) $(node_fault_fix "${OUTDIR:-<search out dir>}" "$JOB")"
elif m "out.of.memory|oom-kill|OUT_OF_MEMORY|std::bad_alloc|cannot allocate";       then err_class="out_of_memory"; fix="Raise the sbatch --mem (e.g. 64G→128G) and resubmit; for DIA-NN try fewer threads or --min-corr.";
elif m "TIMEOUT|DUE TO TIME LIMIT|CANCELLED.*TIME";                                  then err_class="timeout";       fix="Raise --time in the sbatch (or split the run) and resubmit.";
elif m "dotnet: not found|dotnet: command not found";                               then err_class="diann_no_dotnet";fix="Wrong DIA-NN container (no .NET → .raw silently skipped). Use the HIVE native build (build_<v>/diann-<v>/diann-linux) or a .NET-enabled image.";
elif m "Number of IDs at 0.01 FDR: 0";                                              then err_class="diann_zero_ids";  fix="DIA-NN completed but identified nothing - a SILENT null result, not a crash. PRESERVE report.log.txt before re-running; it is the only evidence. Check in order: (1) the library actually generated precursors, (2) the mzML really are DIA with isolation windows (detect_acquisition.py), (3) FASTA matches the organism, (4) mass accuracy. NOTE: DIA-NN warning about generating the predicted library in a separate step is BENIGN per its author - do not chase it.";
elif m "java.lang.OutOfMemoryError|Java heap space|Answer from Java side is empty";  then err_class="spark_heap";      fix="A JVM heap OOM inside Fulcrum/Spark. RAISING --mem DOES NOT HELP: Spark sizes its driver heap independently of the SLURM allocation. Pass spark_config = {\"spark.driver.memory\" = \"32g\"} in the workflow TOML (Fulcrum forwards it to SparkSession.builder.config).";
elif m "0 proteins|No fragment ions|No precursors|no spectra|empty";                then err_class="empty_results"; fix="Check the FASTA matches the organism, the mass-accuracy setting, and that the raw files are the expected acquisition type.";
elif m "CUDA|no kernel image|cuDNN|device-side|GPU.*not";                           then err_class="gpu";           fix="AlphaDIA needs a GPU. Submit to a GPU node (sbatch --gres=gpu:1) or reduce batch size.";
elif m "msconvert.*not found|requires mzML|no mzML";                                then err_class="sage_no_mzml";  fix="Sage/Radiant read mzML and the conversion in the job failed (run_search.py converts .raw with ThermoRawFileParser, .d with msconvert, before the search). Read the converter's own error above the FAILED line: for ThermoRawFileParser usually .NET (bash scripts/ensure_dotnet8.sh), for .raw via msconvert a Linux build without vendor readers. Fix, then regenerate the job with run_search.py --sbatch and resubmit.";
elif m "Disk quota exceeded|No space left";                                         then err_class="disk";          fix="Out of disk/quota. Free space or point --out elsewhere and resubmit.";
elif m "No such file|cannot open|does not exist|not found.*(fasta|\\.d|\\.raw)";     then err_class="missing_input"; fix="An input path is wrong (fasta/raw). Re-check paths (Windows→WSL/HIVE translation) and resubmit.";
elif [ "${dep_failed:-false}" = true ];                                              then err_class="dependency_failed"; fix="An UPSTREAM job in the chain failed, so this one is stuck PENDING with DependencyNeverSatisfied (it will NEVER run and never leave the queue). Find the failed step (sacct -j <arrayjob>), apply that step's fix, and resubmit the downstream steps reusing already-computed outputs (.quant, step1.predicted.speclib) — don't restart the whole chain.";
elif [ "${held:-false}" = true ];                                                   then err_class="held"; fix="The job is HELD (${reason}) and will never start by itself -- do NOT resubmit it (that makes a duplicate). JobHeldUser: release it with scontrol release <jobid> once whatever it was held for is resolved. JobHeldAdmin: ask the HIVE admins why.";
elif $failed;                                                                        then err_class="unknown_failure"; fix="Read the full log; diagnose via references/watcher.md; fix and resubmit.";
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
# Counted from what the 5-step chain actually leaves on disk, not guessed from the log:
# every finished file drops a .quant, so "34 of 66" is a fact. Runs through run() so it
# works over SSH on HIVE exactly as it does locally.
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
    elif has "$OUTDIR/step1.predicted.speclib"; then stage="step2"; n_done="$q2"
    else stage="step1"; n_done=0
    fi
  fi
fi

# ---- Sage LFQ: a finished search can still have unusable quantities -------------------
# sage_lfq_check.py (run in the Sage job and by run_search.py --adapt-only) records whether the
# LFQ window fits the runs' MS1 mass error. A WARNING there is not a failure -- the job succeeded
# and the identifications stand -- so it is reported beside the state, never as error_class.
lfq_rec=""
[ -n "$OUTDIR" ] && lfq_rec="$(run "cat $(printf %q "$OUTDIR/sage_lfq_check.json") 2>/dev/null" 2>/dev/null)"

# ---- narration: what this stage is doing, plus something to read while waiting
if $q_failed; then err_class="watcher_query_failed"; fix="$fix_qf"; fi

notes="$(python3 "$HERE/pipeline_notes.py" --stage "$stage" --index "$POLL" 2>/dev/null)"

STATE="$state" DONE="$done" FAILED="$failed" STALLED="$stalled" JOB="$JOB" MODE="$MODE" \
ECLASS="$err_class" FIX="$fix" TAIL="$tail_txt" ATASKS="${array_summary:-}" \
STAGE="$stage" NDONE="${n_done:-0}" NTOTAL="${n_total:-0}" NOTES="$notes" \
REASON="${reason:-}" ATERM="${a_term:-0}" LFQ_REC="$lfq_rec" NF_JSON="$nf_json" python3 - <<'PY'
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
out["progress"] = progress
try:                                     # node_fault.py's verdict, when it is a node fault
    nf = json.loads(os.environ.get("NF_JSON") or "null")
    if isinstance(nf, dict) and nf.get("node_fault"):
        out["node_fault"] = nf
except ValueError:
    pass
try:
    out.update({k: v for k, v in json.loads(os.environ.get("NOTES") or "{}").items()
                if k in ("doing", "why", "note", "note_source")})
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
out["log_tail"] = tail[-1500:]
print(json.dumps(out, indent=2))
PY
