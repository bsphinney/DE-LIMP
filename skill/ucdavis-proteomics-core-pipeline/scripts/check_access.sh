#!/usr/bin/env bash
# =============================================================================
# check_access.sh  --  Ground the two onboarding questions by probing what's
# actually reachable, and recommend where to run. Emits JSON on stdout.
#
# The skill asks the user:
#   Q1. Do you have access to UC Davis HIVE (account + SSH private key)?
#   Q2. Are you a member of the UC Davis Proteomics Core?
#
# This script verifies those answers against reality so the orchestrator doesn't
# just take the user's word for it:
#   on_hive               this machine is a HIVE node (sbatch present)
#   proteomics_grp_access /quobyte/proteomics-grp is readable here (Core software)
#   hive_ssh              if not on HIVE and a user is given, can we SSH in? and
#                         does that HIVE account have sbatch + proteomics-grp access?
#
# Model: Claude Code runs LOCALLY; HIVE work is driven over SSH with the user's
# private key. So this tests SSH to HIVE using that key.
#
# Usage: bash check_access.sh [hive_user] [private_key_path]
#        (or set HIVE_USER / HIVE_KEY)
# =============================================================================
set -uo pipefail
have() { command -v "$1" >/dev/null 2>&1; }
HU="${1:-${HIVE_USER:-}}"
KEY="${2:-${HIVE_KEY:-}}"
KEY="${KEY/#\~/$HOME}"   # expand a leading ~
HIVE_HOST="${HIVE_HOST:-hive.hpc.ucdavis.edu}"

ON_HIVE=false; have sbatch && ON_HIVE=true
GRP=false; [ -d /quobyte/proteomics-grp ] && ls /quobyte/proteomics-grp >/dev/null 2>&1 && GRP=true

HIVE_SSH="not_tested"; SSH_SBATCH=false; SSH_GRP=false; KEY_FOUND=true
[ -n "$KEY" ] && [ ! -f "$KEY" ] && KEY_FOUND=false
if ! $ON_HIVE && [ -n "$HU" ] && [ "$KEY_FOUND" = true ]; then
  KEY_OPT=""; [ -n "$KEY" ] && KEY_OPT="-i $KEY"
  out="$(timeout 25 ssh $KEY_OPT -o BatchMode=yes -o ConnectTimeout=12 -o IdentitiesOnly=yes "$HU@$HIVE_HOST" \
        'command -v sbatch >/dev/null 2>&1 && echo HAS_SBATCH; ls -d /quobyte/proteomics-grp 2>/dev/null && echo HAS_GRP' 2>/dev/null)"
  if [ -n "$out" ]; then
    HIVE_SSH="ok"
    echo "$out" | grep -q HAS_SBATCH && SSH_SBATCH=true
    echo "$out" | grep -q proteomics-grp && SSH_GRP=true
  else
    HIVE_SSH="failed"   # key/account/VPN problem, or host unreachable
  fi
fi

# Decide the recommended execution mode + facility-software availability.
# Model: Claude Code is LOCAL; HIVE work is driven over SSH with the key.
HAS_SLURM=$([ "$ON_HIVE" = true ] || [ "$SSH_SBATCH" = true ] && echo true || echo false)
FACILITY_SW=$([ "$GRP" = true ] || [ "$SSH_GRP" = true ] && echo true || echo false)
if   $ON_HIVE;                  then MODE="hive_local"   # already on HIVE -> submit SLURM here
elif [ "$SSH_SBATCH" = true ];  then MODE="hive_remote"  # local Claude Code -> drive HIVE over SSH (the intended HIVE mode)
else                                 MODE="local"; fi     # run on the user's own machine

cat <<JSON
{
  "on_hive": $ON_HIVE,
  "key_path_valid": $KEY_FOUND,
  "proteomics_grp_access": $GRP,
  "hive_ssh": "$HIVE_SSH",
  "hive_ssh_has_sbatch": $SSH_SBATCH,
  "hive_ssh_has_proteomics_grp": $SSH_GRP,
  "can_use_slurm": $HAS_SLURM,
  "facility_software_available": $FACILITY_SW,
  "recommended_mode": "$MODE",
  "notes": [
    "Claude Code runs locally; in hive_remote mode it drives HIVE over SSH with the user's private key (ssh -i <key> <user>@hive).",
    "Proteomics Core members (proteomics_grp_access=true) reuse the software already installed in /quobyte/proteomics-grp (DIA-NN .sif, pre-staged FASTAs).",
    "HIVE users NOT in the Core must rebuild the toolchain in their own HIVE home — see references/access.md 'Rebuild on HIVE'.",
    "No HIVE + no Core is fine: the skill installs its own toolchain locally and uses public engines (DIA-NN Academia, Sage).",
    "hive_ssh='failed' usually means VPN off, wrong key path, or account not set up."
  ]
}
JSON
