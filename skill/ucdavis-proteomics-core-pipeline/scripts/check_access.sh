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
#   on_hive                     this machine is a HIVE node (sbatch present)
#   local_proteomics_grp_access /quobyte/proteomics-grp is readable HERE
#   hive_ssh                    if not on HIVE and a user is given, can we SSH in? and
#                               does that HIVE account have sbatch + proteomics-grp access?
#   core_member                 proteomics-grp is readable here OR over SSH
#   hive_ssh_error              why the SSH probe failed: host_key | permission_denied |
#                               timeout | other, with the line that says so
#   hive_host_key_*             is HIVE in known_hosts yet, and its fingerprint
#   local_python3               is there a python3 here that actually runs?
#
# Model: Claude Code runs LOCALLY; HIVE work is driven over SSH with the user's
# private key. So this tests SSH to HIVE using that key.
#
# Usage: bash check_access.sh [hive_user] [private_key_path]
#        (or set HIVE_USER / HIVE_KEY)
# =============================================================================
set -uo pipefail
have() { command -v "$1" >/dev/null 2>&1; }
js() {  # a JSON string literal
  local s="$1"
  s="${s//\\/\\\\}"; s="${s//\"/\\\"}"; s="${s//$'\t'/\\t}"; s="${s//$'\r'/\\r}"; s="${s//$'\n'/\\n}"
  printf '"%s"' "$s"
}
HU="${1:-${HIVE_USER:-}}"
KEY="${2:-${HIVE_KEY:-}}"
KEY="${KEY/#\~/$HOME}"   # expand a leading ~
HIVE_HOST="${HIVE_HOST:-hive.hpc.ucdavis.edu}"
# HIVE's host keys, from ssh-keyscan 2026-09-23 -- the same list is in references/access.md.
PUBLISHED_FP="SHA256:b5nv86Ciaqg1yrUVai6bZ0Hk4IpzAFLWtIPDBdacbQM SHA256:AmJ+z2miIMXlSAcm7k8YwKlIWk5+VqxyT0R4q3fvpcA"

ON_HIVE=false; have sbatch && ON_HIVE=true
GRP=false; [ -d /quobyte/proteomics-grp ] && ls /quobyte/proteomics-grp >/dev/null 2>&1 && GRP=true

# gabrig 2026-09-23, Git Bash on Windows 11: python3 was the Microsoft Store alias in
# ...\WindowsApps, which opens the Store instead of running anything. hive_remote was fine --
# every script ran on HIVE -- but the orchestrator has to know not to run Python locally.
# A macOS /usr/bin/python3 without the Command Line Tools is the same kind of stub (it pops
# an installer), so neither is executed to find out.
PY="$(command -v python3 2>/dev/null || true)"; PY_OK=false; PY_NOTE=""
case "$PY" in
  "")            PY_NOTE="no python3 on PATH" ;;
  *WindowsApps*) PY_NOTE="Windows Store alias (WindowsApps) -- opens the Store instead of running; not executed" ;;
  /usr/bin/python3)
    if [ "$(uname -s 2>/dev/null)" = Darwin ] && ! xcode-select -p >/dev/null 2>&1; then
      PY_NOTE="macOS stub -- the Command Line Tools are not installed; not executed"
    fi ;;
esac
if [ -n "$PY" ] && [ -z "$PY_NOTE" ]; then
  if "$PY" -c 'import sys; sys.exit(sys.version_info[0] != 3)' >/dev/null 2>&1; then PY_OK=true
  else PY_NOTE="python3 is on PATH but did not run"; fi
fi

# gabrig 2026-09-23: on Windows ssh's default user is the AD login (`AD3+gabrig`), and every
# call that relied on it got "Permission denied". HIVE usernames are plain UC Davis ids, so a
# domain-style name is never right -- say so instead of spending a throttled connection on it.
USER_WARN=""
case "$HU" in
  *[+\\@[:space:]]*)
    guess="${HU##*[+\\]}"; guess="${guess%%@*}"
    USER_WARN="HIVE user '$HU' looks like a Windows domain login, not a HIVE username; use the plain UC Davis login id (probably '$guess')." ;;
esac

HIVE_SSH="not_tested"; SSH_SBATCH=false; SSH_GRP=false; KEY_FOUND=true
ERR_KIND=""; ERR_DETAIL=""
KNOWN=null; FPS=""; FP_MATCH=null
[ -n "$KEY" ] && [ ! -f "$KEY" ] && KEY_FOUND=false
if [ -n "$USER_WARN" ] && ! $ON_HIVE; then
  HIVE_SSH="failed"; ERR_KIND="permission_denied"; ERR_DETAIL="not attempted: $USER_WARN"
elif ! $ON_HIVE && [ -n "$HU" ] && [ "$KEY_FOUND" = true ]; then
  # gabrig 2026-09-23: the first BatchMode connection from a new laptop died "Host key
  # verification failed" -- HIVE was not in known_hosts and BatchMode cannot ask. Report
  # whether it is known and, if not, its fingerprint so the user can compare it with the
  # published one; then connect with accept-new, which records an unknown key but still
  # refuses a CHANGED one.
  fp_lines() {  # "TYPE SHA256:..." from `ssh-keygen -F -l` or `ssh-keygen -lf -` output
    awk '!/^#/ { fp = ""; t = ""
      for (i = 1; i <= NF; i++) if ($i ~ /^SHA256:/) { fp = $i; if (i > 1 && $(i-1) !~ /^[0-9]+$/) t = $(i-1) }
      if ($NF ~ /^\(.*\)$/) { t = $NF; gsub(/[()]/, "", t) }
      if (fp != "") print t " " fp }'
  }
  if have ssh-keygen; then
    # No -f first: ssh-keygen then reads ~/.ssh/known_hosts with ~ resolved the way ssh
    # resolves it (the account's home, not $HOME).
    KNOWN=false
    FPS="$(ssh-keygen -F "$HIVE_HOST" -l 2>/dev/null | fp_lines)"
    [ -n "$FPS" ] || FPS="$(ssh-keygen -F "$HIVE_HOST" -f /etc/ssh/ssh_known_hosts -l 2>/dev/null | fp_lines)"
    [ -n "$FPS" ] && KNOWN=true
    # Scan only when unknown: every connection counts against HIVE's MaxStartups throttle.
    if [ "$KNOWN" = false ] && have ssh-keyscan; then
      FPS="$(ssh-keyscan -T 10 -t ed25519,ecdsa "$HIVE_HOST" 2>/dev/null | ssh-keygen -lf - 2>/dev/null | fp_lines)"
    fi
  fi
  if [ -n "$FPS" ]; then
    FP_MATCH=false
    for fp in $(printf '%s\n' "$FPS" | awk '{print $2}'); do
      case " $PUBLISHED_FP " in *" $fp "*) FP_MATCH=true ;; esac
    done
  fi
  # An unknown key that does NOT match the published fingerprint is exactly the case
  # accept-new must not wave through: probe with strict checking, so it fails as host_key
  # and the user checks with HPC support before anything is trusted.
  HKC=accept-new; [ "$KNOWN" = false ] && [ "$FP_MATCH" = false ] && HKC=yes

  # `timeout` is GNU coreutils and a stock macOS does not have it -- there is no
  # /usr/bin/timeout. Unguarded, the probe below died with "command not found", $out
  # came back empty, and hive_ssh was reported "failed"; that then drops MODE to
  # "local", so a user with a perfectly good HIVE key silently ran the whole analysis
  # on their laptop. Degrading to no wrapper is safe: ssh's own ConnectTimeout already
  # bounds a stalled connect, and this only additionally bounds a hang after the
  # connection is up.
  TMO=()
  if   have timeout;  then TMO=(timeout 25)
  elif have gtimeout; then TMO=(gtimeout 25)     # coreutils installed via Homebrew
  fi
  # Arrays, not a spliced string: a key path containing a space has to reach ssh as
  # ONE argument. ${A[@]+"${A[@]}"} because the system bash on macOS is 3.2, where
  # "${A[@]}" on an empty array is an unbound-variable error under `set -u`.
  KEY_OPT=(); [ -n "$KEY" ] && KEY_OPT=(-i "$KEY")
  # bash -l -c, exactly as hive_exec.sh does. ssh with a BARE command runs a NON-login
  # shell, and HIVE does not put sbatch on PATH there -- so this probe reported
  # HAS_SBATCH missing for an account that has it, can_use_slurm went false, and the
  # recommended mode fell back to "local". A Core user with a working cluster account
  # was told to run a multi-hour search on their laptop. Verified against HIVE: bare
  # command -> no HAS_SBATCH, `bash -l -c` -> HAS_SBATCH.
  # HAS_GRP needs a LISTING, not `ls -d`: /quobyte is mode 777 and proteomics-grp is 2770
  # (checked 2026-09-23), so `ls -d` succeeds for every HIVE account, member or not.
  ERRF="$(mktemp "${TMPDIR:-/tmp}/check_access.XXXXXX")"
  out="$(${TMO[@]+"${TMO[@]}"} ssh ${KEY_OPT[@]+"${KEY_OPT[@]}"} \
        -o BatchMode=yes -o ConnectTimeout=12 -o IdentitiesOnly=yes -o StrictHostKeyChecking=$HKC \
        "$HU@$HIVE_HOST" \
        "bash -l -c 'command -v sbatch >/dev/null 2>&1 && echo HAS_SBATCH; ls /quobyte/proteomics-grp >/dev/null 2>&1 && echo HAS_GRP'" 2>"$ERRF")"
  rc=$?
  if [ -n "$out" ]; then
    HIVE_SSH="ok"
    echo "$out" | grep -q HAS_SBATCH && SSH_SBATCH=true
    echo "$out" | grep -q HAS_GRP && SSH_GRP=true
  else
    HIVE_SSH="failed"   # key/account/VPN problem, or host unreachable
    # The line that names the cause, not ssh's "Permanently added" notice.
    err="$(grep -v '^Warning: Permanently added' "$ERRF" | grep -v '^[[:space:]]*$')"
    pick() { printf '%s\n' "$err" | grep -m1 -E "$1"; }
    if   l="$(pick 'Host key verification failed|REMOTE HOST IDENTIFICATION HAS CHANGED|host key .*(differs|is not known)')"; then ERR_KIND=host_key
    elif l="$(pick 'Permission denied')"; then ERR_KIND=permission_denied
    elif l="$(pick '[Tt]imed out|kex_exchange_identification|Connection closed by')"; then ERR_KIND=timeout
    elif [ $rc -eq 124 ]; then ERR_KIND=timeout; l="no answer within 25 s"
    # Not the FIRST line -- that can be HIVE's pre-auth banner -- and not blindly the last:
    # after "Received disconnect ...: <reason>" OpenSSH always ends with a reasonless
    # "Disconnected from <host> port 22" (review 2026-09-23). The disconnect reason if there
    # is one, else ssh's last line.
    else ERR_KIND=other
         l="$(pick 'Received disconnect')"; [ -n "$l" ] || l="$(printf '%s\n' "$err" | tail -1)"; fi
    ERR_DETAIL="${l:-ssh exited $rc with no message}"
    if [ "$ERR_KIND" = host_key ] && [ "$FP_MATCH" = false ]; then
      ERR_DETAIL="$ERR_DETAIL (HIVE's key does not match the published fingerprint -- confirm with HPC support before trusting it)"
    fi
  fi
  rm -f "$ERRF"
fi

# Decide the recommended execution mode + facility-software availability.
# Model: Claude Code is LOCAL; HIVE work is driven over SSH with the key.
HAS_SLURM=$([ "$ON_HIVE" = true ] || [ "$SSH_SBATCH" = true ] && echo true || echo false)
FACILITY_SW=$([ "$GRP" = true ] || [ "$SSH_GRP" = true ] && echo true || echo false)
if   $ON_HIVE;                  then MODE="hive_local"   # already on HIVE -> submit SLURM here
elif [ "$SSH_SBATCH" = true ];  then MODE="hive_remote"  # local Claude Code -> drive HIVE over SSH (the intended HIVE mode)
else                                 MODE="local"; fi     # run on the user's own machine

ERR_JSON=null
[ -n "$ERR_KIND" ] && ERR_JSON="{\"kind\": $(js "$ERR_KIND"), \"detail\": $(js "$ERR_DETAIL")}"
FP_JSON=""
while IFS= read -r l; do [ -n "$l" ] && FP_JSON="$FP_JSON${FP_JSON:+, }$(js "$l")"; done <<EOF
$FPS
EOF
cat <<JSON
{
  "on_hive": $ON_HIVE,
  "key_path_valid": $KEY_FOUND,
  "local_proteomics_grp_access": $GRP,
  "hive_ssh": "$HIVE_SSH",
  "hive_ssh_error": $ERR_JSON,
  "hive_ssh_has_sbatch": $SSH_SBATCH,
  "hive_ssh_has_proteomics_grp": $SSH_GRP,
  "core_member": $FACILITY_SW,
  "hive_host_key_known": $KNOWN,
  "hive_host_key_fingerprints": [$FP_JSON],
  "hive_host_key_matches_published": $FP_MATCH,
  "hive_user_warning": $([ -n "$USER_WARN" ] && js "$USER_WARN" || echo null),
  "local_python3": {"path": $([ -n "$PY" ] && js "$PY" || echo null), "usable": $PY_OK, "note": $(js "$PY_NOTE")},
  "can_use_slurm": $HAS_SLURM,
  "facility_software_available": $FACILITY_SW,
  "recommended_mode": "$MODE",
  "notes": [
    "Claude Code runs locally; in hive_remote mode it drives HIVE over SSH with the user's private key (ssh -i <key> <user>@hive).",
    "core_member = /quobyte/proteomics-grp is readable, here (local_proteomics_grp_access) or over SSH (hive_ssh_has_proteomics_grp). Core members reuse the software already installed there (DIA-NN builds, pre-staged FASTAs).",
    "HIVE users NOT in the Core must rebuild the toolchain in their own HIVE home — see references/access.md 'Rebuild on HIVE'.",
    "No HIVE + no Core is fine: the skill installs its own toolchain locally and uses public engines (DIA-NN Academia, Sage).",
    "hive_ssh='failed': hive_ssh_error says why. host_key = compare hive_host_key_fingerprints with references/access.md; permission_denied = wrong username or key (on Windows never the DOMAIN+user login name); timeout = VPN off, or HIVE's connection throttle after many quick connections (wait ~20 min).",
    "local_python3.usable=false: run nothing in Python on this machine. In hive_remote every script runs on HIVE through hive_exec.sh, which needs no local Python."
  ]
}
JSON
