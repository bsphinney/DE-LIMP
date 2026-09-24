#!/usr/bin/env bash
# =============================================================================
# report_issue.sh  --  Record a problem with THIS SKILL the moment it happens, so the
# Core can fix it. Not for problems with the user's data (a bad run, a failed
# injection) -- those belong in the analysis report.
#
# Why this exists: on 2026-09-23 a Core member (Windows laptop, hive_remote) hit ten
# skill defects in one session -- a mapped T: drive nobody translated to its HIVE path
# (7.1 GB uploaded that was already on HIVE), no rsync in Git Bash, ThermoRawFileParser
# never found, .NET half-installed, FASTA metadata lost. The only reason any of it was
# fixed is that the agent happened to write a notes file into the user's own folder and
# the user happened to forward it. Most sessions end with nothing written down, and
# whatever went wrong is lost with the conversation.
#
#   bash report_issue.sh --title "ThermoRawFileParser not found" \
#        --what "detect_acquisition.py returned 'not found' for all 15 .raw" \
#        --impact "two wasted srun jobs; search would have used the 380-980 fallback" \
#        [--workaround "set THERMORAWFILEPARSER to the Core copy by hand"] \
#        [--fix "setup.sh should find /quobyte/.../ThermoRawFileParser"] \
#        [--kind bug|docs|missing|agent_mistake] [--severity high|medium|low] \
#        [--step "2 detect acquisition"] [--mode hive_remote] [--session <name>]
#   bash report_issue.sh --where        # where issues go from this machine, and nothing else
#   bash report_issue.sh --list [N]     # the N newest issue files (default 20)
#
# Where it goes (first that works):
#   1. On HIVE with the Core's group folder writable  -> $SKILL_ISSUES_DIR directly.
#   2. Off HIVE with a HIVE login set up (hive.env / HIVE_USER + HIVE_KEY) -> the same
#      folder over SSH, through hive_exec.sh. This is the Windows/Mac hive_remote case.
#   3. Otherwise -> ~/.proteomics-pipeline/issues/ on this machine; the script says so and
#      prints how to send the file to the Core. Never an error: recording a problem must not
#      become a second problem.
# One file per user per day per session (<date>_<user>[_<session>].md), appended to, so a
# session's issues read top to bottom in the order they happened and two users never write
# the same file.
#
# Never put secrets in a report: private keys, passwords, CoreOmics/GitHub/HF tokens. The
# script refuses text that looks like one rather than shipping it to a shared folder.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ISSUES_DIR="${SKILL_ISSUES_DIR:-/quobyte/proteomics-grp/skill_issues}"
LOCAL_DIR="${SKILL_ISSUES_LOCAL_DIR:-$HOME/.proteomics-pipeline/issues}"
HIVE_EXEC="${HIVE_EXEC:-$HERE/hive_exec.sh}"
CFG="${HIVE_ENV_FILE:-$HOME/.config/ucdavis-proteomics/hive.env}"
if [ -z "${HIVE_USER:-}" ] && [ -f "$CFG" ]; then . "$CFG"; fi

die() { echo "report_issue.sh: $*" >&2; exit 2; }
TITLE="" WHAT="" IMPACT="" WORKAROUND="" FIX="" KIND="bug" SEV="medium" STEP="" MODE="" SESSION=""
ACTION="record" LIST_N=20
while [ $# -gt 0 ]; do
  case "$1" in
    --title) TITLE="${2:-}"; shift 2 ;;
    --what) WHAT="${2:-}"; shift 2 ;;
    --impact) IMPACT="${2:-}"; shift 2 ;;
    --workaround) WORKAROUND="${2:-}"; shift 2 ;;
    --fix) FIX="${2:-}"; shift 2 ;;
    --kind) KIND="${2:-}"; shift 2 ;;
    --severity) SEV="${2:-}"; shift 2 ;;
    --step) STEP="${2:-}"; shift 2 ;;
    --mode) MODE="${2:-}"; shift 2 ;;
    --session) SESSION="${2:-}"; shift 2 ;;
    --where) ACTION="where"; shift ;;
    --list) ACTION="list"; case "${2:-}" in ''|-*) shift ;; *) LIST_N="$2"; shift 2 ;; esac ;;
    -h|--help) sed -n '2,40p' "$0"; exit 0 ;;
    *) die "unknown argument: $1 (see --help)" ;;
  esac
done
case "$KIND" in bug|docs|missing|agent_mistake) ;; *) die "--kind must be bug, docs, missing or agent_mistake" ;; esac
case "$SEV" in high|medium|low) ;; *) die "--severity must be high, medium or low" ;; esac

# Which transport reaches the shared folder from here: "direct", "ssh", or "local".
on_hive_with_folder() { [ -d "$ISSUES_DIR" ] && [ -w "$ISSUES_DIR" ]; }
have_hive_login() { [ -n "${HIVE_USER:-}" ] && [ -n "${HIVE_KEY:-}" ] && [ -f "$HIVE_EXEC" ]; }
route() {
  if on_hive_with_folder; then echo direct
  elif [ -d "$(dirname "$ISSUES_DIR")" ]; then echo local   # on HIVE, but not in the Core group
  elif have_hive_login; then echo ssh
  else echo local; fi
}

if [ "$ACTION" = "where" ]; then
  case "$(route)" in
    direct) echo "direct: $ISSUES_DIR" ;;
    ssh)    echo "ssh: $HIVE_USER@hive:$ISSUES_DIR (via hive_exec.sh)" ;;
    local)  echo "local: $LOCAL_DIR (send the file to the Core -- see the note it prints)" ;;
  esac
  exit 0
fi
if [ "$ACTION" = "list" ]; then
  case "$LIST_N" in *[!0-9]*|'') die "--list takes a number" ;; esac
  cmd="ls -1t '$ISSUES_DIR'/*.md 2>/dev/null | grep -v '/README.md\$' | head -n $LIST_N"
  case "$(route)" in
    direct) bash -c "$cmd" ;;
    ssh)    bash "$HIVE_EXEC" "$cmd" ;;
    local)  ls -1t "$LOCAL_DIR"/*.md 2>/dev/null | head -n "$LIST_N" ;;
  esac
  exit 0
fi

[ -n "$TITLE" ] && [ -n "$WHAT" ] || die "--title and --what are required"
[ -n "$IMPACT" ] || IMPACT="(not stated)"
# Secrets must never reach a folder 46 people can read. Refuse rather than redact: a
# redacted report can still be wrong about what it hid, and the agent can rephrase.
ALL="$TITLE $WHAT $IMPACT $WORKAROUND $FIX $STEP"
if printf '%s' "$ALL" | grep -Eq -- '-----BEGIN [A-Z ]*PRIVATE KEY|ghp_[A-Za-z0-9]{20,}|github_pat_|hf_[A-Za-z0-9]{20,}|Authorization: *(Token|Bearer) +[A-Za-z0-9]|[Pp]assword *[:=] *[^ ]'; then
  die "the text looks like it contains a key, token or password -- remove it and re-run"
fi

# ---- context captured automatically, so the agent never has to remember it ----------
PLUGIN_JSON="$HERE/../.claude-plugin/plugin.json"
VER="$(sed -nE 's/^[[:space:]]*"version"[[:space:]]*:[[:space:]]*"([^"]+)".*/\1/p' "$PLUGIN_JSON" 2>/dev/null | head -n1)"
[ -n "$VER" ] || VER="unknown"
LOCAL_USER="$(id -un 2>/dev/null || echo "${USERNAME:-${USER:-unknown}}")"
WHO="${HIVE_USER:-$LOCAL_USER}"
OS="$(uname -sr 2>/dev/null || echo unknown)"
DAY="$(date +%Y-%m-%d)"
NOW="$(date '+%Y-%m-%d %H:%M %Z')"
# Filenames: letters, digits, dot, dash, underscore only -- the name travels through a
# remote shell and lands in a shared folder.
clean() { printf '%s' "$1" | tr -c 'A-Za-z0-9._-' '_' | cut -c1-40; }
NAME="${DAY}_$(clean "$WHO")"
[ -n "$SESSION" ] && NAME="${NAME}_$(clean "$SESSION")"
FILE="$NAME.md"
if [ -z "$MODE" ]; then
  case "$(route)" in direct) MODE="hive_local" ;; ssh) MODE="hive_remote" ;; *) MODE="local" ;; esac
fi

MARK="<!-- end of header: report_issue.sh strips everything above this line when appending -->"
ENTRY="$(mktemp "${TMPDIR:-/tmp}/report_issue.XXXXXX")" || die "cannot create a temp file"
trap 'rm -f "$ENTRY"' EXIT
{
  printf '# Skill issues -- %s -- %s%s\n\n' "$WHO" "$DAY" "${SESSION:+ -- $SESSION}"
  printf -- '- **Skill version:** %s\n- **Mode:** %s\n- **OS:** %s\n- **Local user:** %s\n' \
         "$VER" "$MODE" "$OS" "$LOCAL_USER"
  printf -- '- **HIVE user:** %s\n\n' "${HIVE_USER:-(none)}"
  printf '%s\n' "$MARK"
  printf '\n## %s  [%s, %s]\n\n' "$TITLE" "$KIND" "$SEV"
  printf -- '- **When:** %s%s\n' "$NOW" "${STEP:+ -- step $STEP}"
  printf -- '- **Skill:** v%s, mode %s\n' "$VER" "$MODE"
  printf -- '- **What happened:** %s\n' "$WHAT"
  printf -- '- **Impact:** %s\n' "$IMPACT"
  # `if`, not `[ ] &&`: a block's status is its last command's, and a false test there
  # made every entry without --fix look like a failed write.
  if [ -n "$WORKAROUND" ]; then printf -- '- **Workaround used:** %s\n' "$WORKAROUND"; fi
  if [ -n "$FIX" ]; then printf -- '- **Proposed fix:** %s\n' "$FIX"; fi
} > "$ENTRY" || die "cannot write $ENTRY"

# The same append logic runs on either side of the SSH hop: a new file gets the header,
# an existing one gets only what follows the marker line.
APPEND='f="$1"; if [ -s "$f" ]; then sed "1,/^<!-- end of header/d" >> "$f"; else cat > "$f"; fi'

deliver_local() {
  mkdir -p "$LOCAL_DIR" 2>/dev/null || { echo "report_issue.sh: cannot create $LOCAL_DIR" >&2; return 1; }
  bash -c "$APPEND" _ "$LOCAL_DIR/$FILE" < "$ENTRY" || return 1
  echo "recorded locally: $LOCAL_DIR/$FILE"
  echo "  This machine cannot reach the Core's shared issue folder. Please send that file to the"
  echo "  UC Davis Proteomics Core, or open an issue at https://github.com/bsphinney/DE-LIMP/issues"
}

case "$(route)" in
  direct)
    if bash -c "$APPEND" _ "$ISSUES_DIR/$FILE" < "$ENTRY"; then
      chmod g+r "$ISSUES_DIR/$FILE" 2>/dev/null
      echo "recorded: $ISSUES_DIR/$FILE"
    else deliver_local; fi ;;
  ssh)
    # printf %q quoting: the path crosses a remote shell (hive_exec.sh runs `bash -l -c`).
    remote="set -e; [ -d $(printf '%q' "$ISSUES_DIR") ] && [ -w $(printf '%q' "$ISSUES_DIR") ] || exit 7; bash -c $(printf '%q' "$APPEND") _ $(printf '%q' "$ISSUES_DIR/$FILE"); chmod g+r $(printf '%q' "$ISSUES_DIR/$FILE") 2>/dev/null || true"
    if bash "$HIVE_EXEC" "$remote" < "$ENTRY"; then
      echo "recorded: $HIVE_USER@hive:$ISSUES_DIR/$FILE"
    else
      echo "report_issue.sh: could not write to HIVE ($ISSUES_DIR not writable, or SSH failed); keeping it locally" >&2
      deliver_local
    fi ;;
  local) deliver_local ;;
esac
exit 0
