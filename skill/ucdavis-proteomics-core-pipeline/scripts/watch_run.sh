#!/usr/bin/env bash
# =============================================================================
# watch_run.sh  --  One poll of a running search: report its state, detect known
# errors, and suggest the fix. The orchestrator MUST watch every search it starts:
# loop this until the run finishes, and on failure apply the fix (and resubmit) —
# never leave a multi-hour search unmonitored.
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
MODE="local"; JOB=""; LOG=""; HIVE=false
while [ $# -gt 0 ]; do case "$1" in
  --slurm) MODE="slurm"; JOB="$2"; shift 2;;
  --log)   LOG="$2"; shift 2;;
  --hive)  HIVE=true; shift;;
  *) shift;;
esac; done
HERE="$(cd "$(dirname "$0")" && pwd)"
run() { if $HIVE; then bash "$HERE/hive_exec.sh" "$*"; else bash -c "$*"; fi; }

state="unknown"; done=false; failed=false
if [ "$MODE" = "slurm" ] && [ -n "$JOB" ]; then
  st="$(run "sacct -j $JOB --noheader -o State 2>/dev/null | head -1 | tr -d ' '" 2>/dev/null)"
  [ -z "$st" ] && st="$(run "squeue -j $JOB -h -o %T 2>/dev/null" 2>/dev/null)"
  state="${st:-PENDING}"
  case "$state" in
    COMPLETED)                                   done=true;;
    FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL)       done=true; failed=true;;
    CANCELLED*)                                   done=true; failed=true;;
  esac
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
elif m "0 proteins|No fragment ions|No precursors|no spectra|empty";                then err_class="empty_results"; fix="Check the FASTA matches the organism, the mass-accuracy setting, and that the raw files are the expected acquisition type.";
elif m "CUDA|no kernel image|cuDNN|device-side|GPU.*not";                           then err_class="gpu";           fix="AlphaDIA needs a GPU. Submit to a GPU node (sbatch --gres=gpu:1) or reduce batch size.";
elif m "msconvert.*not found|requires mzML|no mzML";                                then err_class="sage_no_mzml";  fix="Sage needs mzML. Convert .d/.raw with msconvert first (Linux/HIVE), then re-run.";
elif m "Disk quota exceeded|No space left";                                         then err_class="disk";          fix="Out of disk/quota. Free space or point --out elsewhere and resubmit.";
elif m "No such file|cannot open|does not exist|not found.*(fasta|\\.d|\\.raw)";     then err_class="missing_input"; fix="An input path is wrong (fasta/raw). Re-check paths (Windows→WSL/HIVE translation) and resubmit.";
elif $failed;                                                                        then err_class="unknown_failure"; fix="Read the full log; diagnose via references/watcher.md; fix and resubmit.";
fi

STATE="$state" DONE="$done" FAILED="$failed" JOB="$JOB" MODE="$MODE" \
ECLASS="$err_class" FIX="$fix" TAIL="$tail_txt" python3 - <<'PY'
import os, json
print(json.dumps({
    "mode": os.environ["MODE"], "job": os.environ["JOB"], "state": os.environ["STATE"],
    "done": os.environ["DONE"] == "true", "failed": os.environ["FAILED"] == "true",
    "error_class": os.environ["ECLASS"], "fix": os.environ["FIX"],
    "log_tail": os.environ["TAIL"][-1500:],
}, indent=2))
PY
