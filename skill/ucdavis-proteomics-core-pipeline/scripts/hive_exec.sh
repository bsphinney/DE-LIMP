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
#     bash hive_exec.sh 'ls -d /quobyte/proteomics-grp/dia-nn/build_*/diann-*'
#
#   Copy files to/from HIVE (helpers):
#     bash hive_exec.sh --put  ./local/path   '~/remote/path'
#     bash hive_exec.sh --get  '~/remote/path' ./local/path
#
#   Calls share one SSH connection for 10 minutes (ControlMaster); HIVE_SSH_MUX=0 disables.
#
#   --put refuses a source that is already on HIVE (a mapped drive or SMB mount of a
#   share HIVE mounts -- see hive_path.sh) and prints the HIVE path to use in place;
#   HIVE_PUT_FORCE=1 uploads anyway. With no rsync (Git Bash on Windows has none),
#   --put/--get use scp.
#
# Heavy compute must go through SLURM (sbatch), never the login node.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
# Env vars do not survive between tool calls / shells, so requiring HIVE_USER in the
# environment means every helper that shells out (watch_run.sh especially) silently
# fails with an empty result -- which downstream code reads as "nothing there".
# Persist it once to a config file and every later invocation just works.
CFG="${HIVE_ENV_FILE:-$HOME/.config/ucdavis-proteomics/hive.env}"
if [ -z "${HIVE_USER:-}" ] && [ -f "$CFG" ]; then . "$CFG"; fi
HU="${HIVE_USER:?set HIVE_USER, or save it once: mkdir -p ~/.config/ucdavis-proteomics && printf \'HIVE_USER=<user>\\nHIVE_KEY=~/.ssh/id_ed25519\\n\' > ~/.config/ucdavis-proteomics/hive.env}"
# gabrig 2026-09-23, Windows 11: ssh's default user there is the AD login `AD3+gabrig`, and
# every call that relied on it got "Permission denied". The same name in HIVE_USER fails the
# same way, after spending one of HIVE's throttled connection attempts -- so stop here.
case "$HU" in
  *[+\\@[:space:]]*)
    guess="${HU##*[+\\]}"; guess="${guess%%@*}"
    echo "HIVE_USER='$HU' looks like a Windows domain login, not a HIVE username." >&2
    echo "HIVE usernames are the plain UC Davis login id (probably '$guess'): export HIVE_USER=$guess" >&2
    exit 2 ;;
esac
KEY="${HIVE_KEY:?set HIVE_KEY to the private-key path the user gave you}"
KEY="${KEY/#\~/$HOME}"
HOST="${HIVE_HOST:-hive.hpc.ucdavis.edu}"
[ -f "$KEY" ] || { echo "private key not found: $KEY" >&2; exit 2; }

# Connection reuse. HIVE's sshd throttles rapid NEW connections (MaxStartups): measured
# 2026-09-16, 8 quick separate ssh calls got "kex_exchange_identification: read: Operation
# timed out" for 20+ minutes while the HIVE status page said operational. A multiplexed
# master makes every later call a channel on one connection instead of a new handshake.
# %C (a hash of host/user/port) keeps the socket path short -- a spelled-out ControlPath
# can overflow the ~104-byte unix-socket limit. A stale master prints
# "mux_client_request_session: read from master failed: Broken pipe" and ssh falls back to
# a normal connection, which is harmless. Native Windows OpenSSH has no ControlMaster, so
# it is off there; HIVE_SSH_MUX=0 turns it off anywhere.
MUX=()
MUX_E=""
case "$(uname -s 2>/dev/null)" in MINGW*|MSYS*|CYGWIN*) HIVE_SSH_MUX=0 ;; esac
if [ "${HIVE_SSH_MUX:-1}" != "0" ]; then
  [ -d "$HOME/.ssh" ] || mkdir -m 700 -p "$HOME/.ssh"
  MUX=(-o ControlMaster=auto -o "ControlPath=$HOME/.ssh/cm-%C" -o ControlPersist=10m)
  MUX_E=" -o ControlMaster=auto -o 'ControlPath=$HOME/.ssh/cm-%C' -o ControlPersist=10m"
fi
# StrictHostKeyChecking=accept-new: gabrig 2026-09-23, the first call from a new laptop died
# "Host key verification failed" -- HIVE was not in known_hosts yet and there is no prompt
# to answer in a non-interactive call. accept-new records an UNKNOWN key and still refuses
# a CHANGED one. check_access.sh prints the fingerprint for the user to compare.
HK=(-o StrictHostKeyChecking=accept-new)
# Keep-alives: a laptop that sleeps or a VPN that reconnects inside ControlPersist leaves the
# master on a dead TCP connection, and every later call multiplexes onto it and hangs until
# TCP gives up -- minutes, with ConnectTimeout not applying (review 2026-09-23). 15 s x 3 =
# the dead connection is dropped within ~45 s and the next call opens a fresh one.
HK+=(-o ServerAliveInterval=15 -o ServerAliveCountMax=3)
# ${MUX[@]+...}: an empty array under `set -u` is an "unbound variable" error in bash < 4.4,
# and macOS still ships bash 3.2 as /bin/bash.
SSH=(ssh -i "$KEY" -o IdentitiesOnly=yes -o ConnectTimeout=20 "${HK[@]}" ${MUX[@]+"${MUX[@]}"} "$HU@$HOST")
# $KEY is single-quoted inside -e: rsync parses that string itself and DOES honour
# quotes (verified -- a path with a space arrives as one argv entry), whereas bare
# $KEY gets split and ssh is handed a truncated path. Matches the SSH=() array above.
RSYNC_E="ssh -i '$KEY' -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new -o ServerAliveInterval=15 -o ServerAliveCountMax=3$MUX_E"
SCP=(scp -i "$KEY" -o IdentitiesOnly=yes "${HK[@]}" ${MUX[@]+"${MUX[@]}"} -r)

# rsync and scp both read the "T" of T:\x as a HOST name. On a Windows shell (Git Bash, MSYS2,
# Cygwin) spell a drive-letter path the shell's own way first -- cygpath knows /t/x vs
# /cygdrive/t/x -- for rsync as well as scp (review 2026-09-23: only the scp path did it).
# Anywhere else a path is left exactly as given.
winpath() {
  case "$(uname -s 2>/dev/null)" in
    MINGW*|MSYS*|CYGWIN*)
      case "$1" in
        [A-Za-z]:[\\/]*)
          # A trailing \ (T:\Data\P1\) names the folder in Windows terms, but a trailing /
          # means "copy the CONTENTS" to rsync -- so drop it (review 2026-09-23).
          local p
          if command -v cygpath >/dev/null 2>&1; then p="$(cygpath -u "$1")"
          else p="/$(printf '%s' "${1%%:*}" | tr '[:upper:]' '[:lower:]')/${1:3}"; p="${p//\\//}"; fi
          while [ "${#p}" -gt 3 ] && [ "${p%/}" != "$p" ]; do p="${p%/}"; done
          printf '%s\n' "$p"; return ;;
      esac ;;
  esac
  printf '%s\n' "$1"
}

usage() { echo "usage: hive_exec.sh '<command>' | --put <local> <remote> | --get <remote> <local>" >&2; exit 2; }

# gabrig 2026-09-23: T:\Data\lab\service\... went up by scp (7.1 GB, ~20 min) although T: is
# the Flinders share HIVE mounts at /nfs/lssc0/flinders/proteomics. hive_path.sh makes no
# ssh call for an ordinary local path, so this costs nothing unless the source is on a share.
put_guard() {
  [ "${HIVE_PUT_FORCE:-0}" = 1 ] && return 0
  local j rc
  j="$(bash "$HERE/hive_path.sh" "$1")"; rc=$?
  case $rc in
    0) echo "REFUSING --put: '$1' is already on HIVE -- use the hive_path below in place." >&2
       printf '%s\n' "$j" | grep -E '^  "(hive_path|how)"' >&2
       echo "(HIVE_PUT_FORCE=1 uploads it anyway.)" >&2
       exit 4 ;;
    1) echo "note: '$1' is on a network share but no HIVE copy was verified; uploading:" >&2
       printf '%s\n' "$j" | grep -E '^  "how"' >&2 ;;
    3) ;;
    *) echo "note: hive_path.sh failed (exit $rc); uploading without the on-HIVE check." >&2 ;;
  esac
}

# scp/sftp does not reliably expand a leading ~/ (OpenSSH >= 9 scp speaks SFTP). A relative
# remote path is relative to home under both protocols, so drop it.
rpath() { local r="$1"; case "$r" in "~") r=. ;; "~/"*) r="${r#\~/}"; [ -n "$r" ] || r=. ;; esac; echo "$r"; }

# rsync, which the rest of the skill assumes, is not in Git Bash on Windows (gabrig
# 2026-09-23: "rsync: command not found", then scp commands hand-built by the agent). These
# keep rsync's meaning, which scp -r does not have on its own:
#   - a directory lands INSIDE <dest>, and <dest> is created if missing. scp -r to a
#     missing <dest> would instead make <dest> itself the copy.
#   - a trailing slash on the source ("dir/" = copy the CONTENTS) cannot be said with scp,
#     which copies the directory itself either way -- refused rather than silently changed.
no_contents_copy() {
  echo "hive_exec.sh: '$1' ends in '/', which under rsync copies the folder's CONTENTS." >&2
  echo "Without rsync (scp fallback) that cannot be done; drop the trailing slash to copy" >&2
  echo "the folder itself into the destination, or install rsync." >&2
  exit 2
}
scp_put() {
  local src="$1" dst; dst="$(rpath "$2")"
  case "$src" in */) [ -d "$src" ] && no_contents_copy "$src" ;; esac
  # scp reads the "T" of T:\x as a host name; Git Bash spells that drive /t/x.
  case "$src" in [A-Za-z]:[\\/]*) src="/$(printf '%s' "${src%%:*}" | tr '[:upper:]' '[:lower:]')/${src:3}"; src="${src//\\//}" ;; esac
  if [ -d "$src" ] || [ "${dst%/}" != "$dst" ]; then
    "${SSH[@]}" "mkdir -p -- $(printf '%q' "$dst")" || exit
    dst="${dst%/}/"
  fi
  "${SCP[@]}" "$src" "$HU@$HOST:$dst"
}
scp_get() {
  local src dst="$2" b rc
  case "$1" in */) no_contents_copy "$1" ;; esac
  src="$(rpath "$1")"
  if [ -d "$dst" ] || [ "${dst%/}" != "$dst" ]; then
    mkdir -p "$dst" && "${SCP[@]}" "$HU@$HOST:$src" "${dst%/}/"; return
  elif [ -e "$dst" ]; then
    "${SCP[@]}" "$HU@$HOST:$src" "$dst"; return                    # overwrite a file, as rsync does
  fi
  # A missing <dest>: rsync names a single FILE <dest>, but puts a DIRECTORY at <dest>/<name>.
  # Which it is shows only after the copy, so land it inside a new <dest>/ and unwrap a file.
  mkdir -p "$dst" || exit
  "${SCP[@]}" "$HU@$HOST:$src" "$dst/"; rc=$?
  [ $rc -eq 0 ] || { rmdir "$dst" 2>/dev/null; return $rc; }
  b="$(basename "$src")"
  if [ -f "$dst/$b" ] && [ "$(ls -A "$dst")" = "$b" ]; then
    mv "$dst/$b" "$dst.part.$$" && rmdir "$dst" && mv "$dst.part.$$" "$dst"
  fi
}

case "${1:-}" in
  --put) [ $# -eq 3 ] || usage; shift; put_guard "$1"; src="$(winpath "$1")"
         if command -v rsync >/dev/null 2>&1; then rsync -e "$RSYNC_E" -a "$src" "$HU@$HOST:$2"
         else scp_put "$src" "$2"; fi ;;
  # --get is -rlt, not -a: HIVE's shared dirs are setgid 2775, and -a's permission/group copy
  # makes macOS fail with "fchmodat ... Operation not permitted" (exit 23) AFTER the files
  # arrive -- measured 2026-09-16 -- so a good copy read as a failed one.
  --get) [ $# -eq 3 ] || usage; shift
         dst="$(winpath "$2")"
         if command -v rsync >/dev/null 2>&1; then rsync -e "$RSYNC_E" -rlt "$HU@$HOST:$1" "$dst"
         else scp_get "$1" "$dst"; fi ;;
  "")    usage ;;
  # LOGIN shell (bash -l). ssh with a bare command runs a NON-login shell, where
  # HIVE does not put sacct/squeue/sbatch on PATH. Every SLURM query then returned
  # nothing, and watch_run.sh read that emptiness as "PENDING" -- so a job array with
  # 14 TIMEOUT tasks reported as still running. Silence must never look like health.
  *)     "${SSH[@]}" "bash -l -c $(printf '%q' "$*")" ;;
esac
