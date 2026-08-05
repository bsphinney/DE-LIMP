#!/usr/bin/env bash
# =============================================================================
# hive_exec.sh  --  Run a command on UC Davis HIVE over SSH using the user's
# private key. Claude Code runs LOCALLY; this is how the HIVE steps execute.
#
#   Set once per session:
#     export HIVE_USER=brettsp
#     export HIVE_KEY=~/.ssh/id_ed25519     # the path the user gave you
#
#   Run a command on HIVE:
#     bash hive_exec.sh 'sbatch ~/run/job.sh'
#     bash hive_exec.sh 'ls -d /quobyte/proteomics-grp/dia-nn/*.sif'
#
#   Copy files to/from HIVE (helpers):
#     bash hive_exec.sh --put  ./local/path   '~/remote/path'
#     bash hive_exec.sh --get  '~/remote/path' ./local/path
#
# Heavy compute must go through SLURM (sbatch), never the login node.
# =============================================================================
set -uo pipefail
# Env vars do not survive between tool calls / shells, so requiring HIVE_USER in the
# environment means every helper that shells out (watch_run.sh especially) silently
# fails with an empty result -- which downstream code reads as "nothing there".
# Persist it once to a config file and every later invocation just works.
CFG="${HIVE_ENV_FILE:-$HOME/.config/ucdavis-proteomics/hive.env}"
if [ -z "${HIVE_USER:-}" ] && [ -f "$CFG" ]; then . "$CFG"; fi
HU="${HIVE_USER:?set HIVE_USER, or save it once: mkdir -p ~/.config/ucdavis-proteomics && printf \'HIVE_USER=<user>\\nHIVE_KEY=~/.ssh/id_ed25519\\n\' > ~/.config/ucdavis-proteomics/hive.env}"
KEY="${HIVE_KEY:?set HIVE_KEY to the private-key path the user gave you}"
KEY="${KEY/#\~/$HOME}"
HOST="${HIVE_HOST:-hive.hpc.ucdavis.edu}"
[ -f "$KEY" ] || { echo "private key not found: $KEY" >&2; exit 2; }
SSH=(ssh -i "$KEY" -o IdentitiesOnly=yes -o ConnectTimeout=20 "$HU@$HOST")

case "${1:-}" in
  --put) shift; rsync -e "ssh -i $KEY -o IdentitiesOnly=yes" -a "$1" "$HU@$HOST:$2" ;;
  --get) shift; rsync -e "ssh -i $KEY -o IdentitiesOnly=yes" -a "$HU@$HOST:$1" "$2" ;;
  "")    echo "usage: hive_exec.sh '<command>' | --put <local> <remote> | --get <remote> <local>" >&2; exit 2 ;;
  # LOGIN shell (bash -l). ssh with a bare command runs a NON-login shell, where
  # HIVE does not put sacct/squeue/sbatch on PATH. Every SLURM query then returned
  # nothing, and watch_run.sh read that emptiness as "PENDING" -- so a job array with
  # 14 TIMEOUT tasks reported as still running. Silence must never look like health.
  *)     "${SSH[@]}" "bash -l -c $(printf '%q' "$*")" ;;
esac
