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

state="unknown"; done=false; failed=false; q_failed=false; fix_qf=""

# --all <dir>: watch the WHOLE chain, not just its last link. Watching only the final
# job is why a step-4 array that timed out went unnoticed for hours -- step 5 simply
# sat PENDING on a dependency that could never be satisfied, which reads as "still
# running". jobs.txt is written by submit.sh with every id in the chain.
if [ "$MODE" = "chain" ]; then
  JT="$CHAINDIR/jobs.txt"
  probe="$(run "echo __ok__" 2>&1)"
  if ! printf '%s' "$probe" | grep -q __ok__; then
    echo "{\"mode\":\"chain\",\"failed\":true,\"done\":true,\"err_class\":\"watcher_query_failed\",\"detail\":$(printf '%s' "$probe" | head -2 | python3 -c 'import sys,json;print(json.dumps(sys.stdin.read()))'),\"fix\":\"Cannot reach the cluster at all — this is NOT a job failure and NOT a missing file. Usually HIVE_USER/HIVE_KEY are unset: save them once to ~/.config/ucdavis-proteomics/hive.env. Re-run the watcher before drawing any conclusion about the run.\"}"
    exit 3
  fi
  ids="$(run "cat $JT 2>/dev/null | tr '\n' ' '" 2>/dev/null)"
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
  if $anyfail; then
    echo "{\"mode\":\"chain\",\"dir\":\"$CHAINDIR\",\"failed\":true,\"done\":true,\"first_failed_job\":\"$worst\",\"chain\":\"${summary# }\",\"fix\":\"A step FAILED. Inspect it: watch_run.sh --slurm $worst --hive. Downstream steps will sit PENDING with DependencyNeverSatisfied forever until you fix and resubmit them.\"}"
    exit 0
  fi
  if $allterm; then
    echo "{\"mode\":\"chain\",\"dir\":\"$CHAINDIR\",\"failed\":false,\"done\":true,\"chain\":\"${summary# }\"}"
  else
    echo "{\"mode\":\"chain\",\"dir\":\"$CHAINDIR\",\"failed\":false,\"done\":false,\"chain\":\"${summary# }\"}"
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
fi

tail_txt=""
[ -n "$LOG" ] && tail_txt="$(run "tail -n 100 $(printf %q "$LOG") 2>/dev/null" 2>/dev/null)"

# error signatures -> (class, fix). First match wins.
err_class=""; fix=""
hay="$tail_txt
$state"
m() { printf '%s' "$hay" | grep -qiE "$1"; }
if   m "out.of.memory|oom-kill|OUT_OF_MEMORY|std::bad_alloc|cannot allocate";       then err_class="out_of_memory"; fix="Raise the sbatch --mem (e.g. 64G→128G) and resubmit; for DIA-NN try fewer threads or --min-corr.";
elif m "TIMEOUT|DUE TO TIME LIMIT|CANCELLED.*TIME";                                  then err_class="timeout";       fix="Raise --time in the sbatch (or split the run) and resubmit.";
elif m "dotnet: not found|dotnet: command not found";                               then err_class="diann_no_dotnet";fix="Wrong DIA-NN container (no .NET → .raw silently skipped). Use the HIVE native build (build_<v>/diann-<v>/diann-linux) or a .NET-enabled image.";
elif m "Number of IDs at 0.01 FDR: 0";                                              then err_class="diann_zero_ids";  fix="DIA-NN completed but identified nothing - a SILENT null result, not a crash. PRESERVE report.log.txt before re-running; it is the only evidence. Check in order: (1) the library actually generated precursors, (2) the mzML really are DIA with isolation windows (detect_acquisition.py), (3) FASTA matches the organism, (4) mass accuracy. NOTE: DIA-NN warning about generating the predicted library in a separate step is BENIGN per its author - do not chase it.";
elif m "java.lang.OutOfMemoryError|Java heap space|Answer from Java side is empty";  then err_class="spark_heap";      fix="A JVM heap OOM inside Fulcrum/Spark. RAISING --mem DOES NOT HELP: Spark sizes its driver heap independently of the SLURM allocation. Pass spark_config = {\"spark.driver.memory\" = \"32g\"} in the workflow TOML (Fulcrum forwards it to SparkSession.builder.config).";
elif m "0 proteins|No fragment ions|No precursors|no spectra|empty";                then err_class="empty_results"; fix="Check the FASTA matches the organism, the mass-accuracy setting, and that the raw files are the expected acquisition type.";
elif m "CUDA|no kernel image|cuDNN|device-side|GPU.*not";                           then err_class="gpu";           fix="AlphaDIA needs a GPU. Submit to a GPU node (sbatch --gres=gpu:1) or reduce batch size.";
elif m "msconvert.*not found|requires mzML|no mzML";                                then err_class="sage_no_mzml";  fix="Sage needs mzML. Convert .d/.raw with msconvert first (Linux/HIVE), then re-run.";
elif m "Disk quota exceeded|No space left";                                         then err_class="disk";          fix="Out of disk/quota. Free space or point --out elsewhere and resubmit.";
elif m "No such file|cannot open|does not exist|not found.*(fasta|\\.d|\\.raw)";     then err_class="missing_input"; fix="An input path is wrong (fasta/raw). Re-check paths (Windows→WSL/HIVE translation) and resubmit.";
elif [ "${dep_failed:-false}" = true ];                                              then err_class="dependency_failed"; fix="An UPSTREAM job in the chain failed, so this one is stuck PENDING with DependencyNeverSatisfied (it will NEVER run and never leave the queue). Find the failed step (sacct -j <arrayjob>), apply that step's fix, and resubmit the downstream steps reusing already-computed outputs (.quant, step1.predicted.speclib) — don't restart the whole chain.";
elif $failed;                                                                        then err_class="unknown_failure"; fix="Read the full log; diagnose via references/watcher.md; fix and resubmit.";
fi

# STALL detection: job RUNNING but its log has not advanced in STALL_MIN minutes.
# A hung file (e.g. a pathological DIA-NN run) keeps the job in RUNNING while the log
# freezes — sacct/squeue will NOT flag it. This is the "NA41-class" failure.
stalled=false
if [ -z "$err_class" ] && [ -n "$LOG" ] && printf '%s' "$state" | grep -qiE "RUNNING|^R$"; then
  now="$(run "date +%s" 2>/dev/null)"
  mt="$(run "stat -c %Y $(printf %q "$LOG") 2>/dev/null" 2>/dev/null | tr -dc 0-9)"
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

# ---- narration: what this stage is doing, plus something to read while waiting
if $q_failed; then err_class="watcher_query_failed"; fix="$fix_qf"; fi

notes="$(python3 "$HERE/pipeline_notes.py" --stage "$stage" --index "$POLL" 2>/dev/null)"

STATE="$state" DONE="$done" FAILED="$failed" STALLED="$stalled" JOB="$JOB" MODE="$MODE" \
ECLASS="$err_class" FIX="$fix" TAIL="$tail_txt" ATASKS="${array_summary:-}" \
STAGE="$stage" NDONE="${n_done:-0}" NTOTAL="${n_total:-0}" NOTES="$notes" python3 - <<'PY'
import os, json
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
out["progress"] = progress
try:
    out.update({k: v for k, v in json.loads(os.environ.get("NOTES") or "{}").items()
                if k in ("doing", "why", "note", "note_source")})
except Exception:
    pass
out["log_tail"] = os.environ["TAIL"][-1500:]
print(json.dumps(out, indent=2))
PY
