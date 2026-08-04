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
HU="${HIVE_USER:?set HIVE_USER (ask the user for their HIVE username)}"
KEY="${HIVE_KEY:?set HIVE_KEY to the private-key path the user gave you}"
KEY="${KEY/#\~/$HOME}"
HOST="${HIVE_HOST:-hive.hpc.ucdavis.edu}"
[ -f "$KEY" ] || { echo "private key not found: $KEY" >&2; exit 2; }
SSH=(ssh -i "$KEY" -o IdentitiesOnly=yes -o ConnectTimeout=20 "$HU@$HOST")

case "${1:-}" in
  --put) shift; rsync -e "ssh -i $KEY -o IdentitiesOnly=yes" -a "$1" "$HU@$HOST:$2" ;;
  --get) shift; rsync -e "ssh -i $KEY -o IdentitiesOnly=yes" -a "$HU@$HOST:$1" "$2" ;;
  "")    echo "usage: hive_exec.sh '<command>' | --put <local> <remote> | --get <remote> <local>" >&2; exit 2 ;;
  *)     "${SSH[@]}" "$@" ;;
esac
