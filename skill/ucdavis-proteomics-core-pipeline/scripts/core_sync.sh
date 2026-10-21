#!/usr/bin/env bash
# =============================================================================
# core_sync.sh  --  Copy a CoreOmics submission's session between HIVE, which holds the
# session of record, and this computer, where the report is written (SKILL.md step 1c.7).
# The file lists live here, once, so a hand-copied loop cannot drop or mistype one.
#
#   Pull what the local report needs (after the HIVE steps, before step 9):
#     bash core_sync.sh --get --session "$S" --local ~/core/PROT_0807/session
#   Push the finished report files back BEFORE deliver (deliver copies from $S):
#     bash core_sync.sh --put --session "$S" --local ~/core/PROT_0807/session
#   --dry-run prints what would move and moves nothing (no HIVE connection at all).
#
#   --get makes <local>/output/search and <local>/input, then copies each GET item from
#   <session>/<item> into <local>/<the item's folder>/. One listing call on HIVE first says
#   which items exist: an absent one (AUDIT.md before step 8c, no conditions.csv, ...) is
#   "missing", not an error.
#   --put copies each PUT item that exists under <local> into <session>/output/; the rest are
#   "skipped (not here)". Both refuse a --session that is not a directory on HIVE.
#
#   --session  the HIVE session dir ($S): absolute (/quobyte/...) or ~/... . It goes over
#              ssh, so only letters, digits and / . _ - + , = @ % are accepted -- no spaces,
#              quotes, $, backticks, \, ;, |, &, newlines.
#   --local    this computer's copy of the session (a leading ~/ is expanded).
#
# Output: one line per item -- "got <session>/<item>", "missing <session>/<item>",
# "put <local>/<item>", "skipped <local>/<item> (not here)", "failed <path> (...)" -- then
# one JSON line: {"mode": "get"|"put", "session": ..., "local": ..., "done": [...],
# "missing": [...], "failed": [...]}, items relative to the session ("." = the session
# itself, with an "error"); --dry-run adds "dry_run": true and "would": [...].
# Exit: 0 every item moved or was absent; 1 a transfer failed, HIVE could not be reached or
# the session is not there; 2 a usage error.
#
# Transport: hive_exec.sh beside this file (HIVE_EXEC overrides it, as in report_issue.sh):
# rsync, or scp where there is none (Git Bash on Windows). Bash 3.2 and Git Bash safe.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HIVE_EXEC="${HIVE_EXEC:-$HERE/hive_exec.sh}"

# --get, from <session>/<item> into <local>/<folder>/ -- what the local report reads:
# make_analysis_html.py (step 9) the DE tables, figures, AUDIT.md and SAMPLE_QUALITY.md;
# analysis_prompt.py the search report and methods.md.
GET_OUTPUT="output/tables output/figures output/search/report.parquet output/AUDIT.md
            output/SAMPLE_QUALITY.md output/methods.md"
# --get: the session records the report's sections read -- make_analysis_html.py through
# submission_report.py (session.json, the attached submission, input/locate_decisions.json),
# the conditions and the FASTA's metadata, and input/sample_type.json and experiment_type.json
# (sample_type.py / experiment_type.py, read by analysis_prompt.py and session_docs.py).
GET_RECORDS="session.json input/submission.json input/samples.tsv input/raw_files.txt
             input/conditions.csv input/search.fasta.meta.json input/sample_type.json
             input/experiment_type.json input/locate_decisions.json"
# --put, from <local>/<item> into <session>/output/ when it exists here -- what deliver
# (core_submission.py DELIVER_FILES and the podcast) copies from $S. The PDF exists only where
# a headless browser printed it; the podcast only when it was made on this computer.
PUT_OUTPUT="output/AI_Analysis_Report.md output/Analysis_Report.html output/Analysis_Report.md
            output/Analysis_Report.pdf output/Analysis_Report_with_audio.html output/methods.md
            output/methods.docx output/podcast"

USAGE="usage: core_sync.sh --get|--put --session <HIVE session dir> --local <local session dir> [--dry-run]"
die() { echo "core_sync.sh: $*" >&2; echo "$USAGE" >&2; exit 2; }
NL=$'\n'

MODE="" S="" L="" DRY=0 HAVE_S=0 HAVE_L=0
while [ $# -gt 0 ]; do
  case "$1" in
    -h|--help) awk 'NR > 2 && /^# =====/ { exit } NR > 2 { sub(/^# ?/, ""); print }' "$0"
               echo "$USAGE"; exit 0 ;;
    --get|--put) [ -z "$MODE" ] || die "give one of --get or --put, once"
                 MODE="${1#--}"; shift ;;
    --session|--local)
      [ $# -ge 2 ] || die "$1 needs a value"
      if [ "$1" = --session ]; then S="$2"; HAVE_S=1; else L="$2"; HAVE_L=1; fi
      shift 2 ;;
    --session=*) S="${1#*=}"; HAVE_S=1; shift ;;
    --local=*) L="${1#*=}"; HAVE_L=1; shift ;;
    --dry-run) DRY=1; shift ;;
    *) die "unknown argument: $1" ;;
  esac
done
[ -n "$MODE" ] || die "give --get (pull from HIVE) or --put (push to HIVE)"
[ "$HAVE_S" = 1 ] && [ -n "$S" ] || die "--session is required: the session dir on HIVE (\$S)"
[ "$HAVE_L" = 1 ] && [ -n "$L" ] || die "--local is required: this computer's copy of the session"

# --session goes into a command line on HIVE: an allow-list, not a list of bad characters.
# LC_ALL=C (in the subshell only) keeps A-Z an ASCII range in every locale.
session_chars_ok() ( LC_ALL=C; case "${1#\~/}" in *[!A-Za-z0-9/._+,=@%-]*) exit 1 ;; esac; exit 0 )
while [ "${#S}" -gt 1 ] && [ "${S%/}" != "$S" ]; do S="${S%/}"; done
case "$S" in
  /?*|"~/"?*) ;;
  *) die "--session must be the session's absolute path on HIVE (/quobyte/...) or ~/...: '$S'" ;;
esac
session_chars_ok "$S" || die "--session holds a character the remote shell would interpret
(only letters, digits and / . _ - + , = @ % are accepted): '$S'"

case "$L" in "~") L="$HOME" ;; "~/"*) L="$HOME/${L#\~/}" ;; esac
case "$(uname -s 2>/dev/null)" in
  MINGW*|MSYS*|CYGWIN*) case "$L" in
    [A-Za-z]:[\\/]*) command -v cygpath >/dev/null 2>&1 && L="$(cygpath -u "$L")" ;; esac ;;
esac
while [ "${#L}" -gt 1 ] && [ "${L%/}" != "$L" ]; do L="${L%/}"; done

# JSON string, one character at a time: parameter-expansion escaping of \ and " differs
# between bash 3.2 and 5.x.
js() {
  local s="$1" out="" c i
  for ((i = 0; i < ${#s}; i++)); do
    c="${s:i:1}"
    case "$c" in
      '"') out="$out\\\"" ;;
      '\') out="$out\\\\" ;;
      "$NL") out="$out\\n" ;;
      $'\t') out="$out\\t" ;;
      $'\r') out="$out\\r" ;;
      [[:cntrl:]]) out="$out$(printf '\\u%04x' "'$c")" ;;
      *) out="$out$c" ;;
    esac
  done
  printf '"%s"' "$out"
}
DONE_J="" MISS_J="" FAIL_J="" WOULD_J="" ERR=""
summary() {
  printf '{"mode": %s, "session": %s, "local": %s, "done": [%s], "missing": [%s], "failed": [%s]' \
    "$(js "$MODE")" "$(js "$S")" "$(js "$L")" "$DONE_J" "$MISS_J" "$FAIL_J"
  if [ "$DRY" = 1 ]; then printf ', "dry_run": true, "would": [%s]' "$WOULD_J"; fi
  if [ -n "$ERR" ]; then printf ', "error": %s' "$(js "$ERR")"; fi
  printf '}\n'
}

# The remote form of $S: single-quoted (it holds no quote), a leading ~/ left outside the
# quotes so HIVE's shell still expands it.
rq() { case "$1" in "~/"*) printf "~/'%s'" "${1#\~/}" ;; *) printf "'%s'" "$1" ;; esac; }

# ONE call on HIVE: is $S a directory, and which of the items "$@" exist in it. Tells "not
# there" from a transport failure, which hive_exec.sh's exit codes cannot (rsync and scp both
# exit non-zero for a missing source, scp with the same 1 as for any other failure).
# Sets HAVE (marker lines); returns 0 = listed, 1 = no such session, 2 = HIVE not reached.
probe() {
  local cmd out rc
  cmd="cd $(rq "$S") 2>/dev/null || { echo CORE_SYNC_NO_SESSION; exit 0; }"
  [ $# -eq 0 ] || cmd="$cmd; for p in $*; do [ -e \"\$p\" ] && echo \"CORE_SYNC_HAVE \$p\"; done"
  cmd="$cmd; echo CORE_SYNC_END"
  out="$(bash "$HIVE_EXEC" "$cmd" </dev/null)"; rc=$?
  HAVE="$NL${out//$'\r'/}$NL"
  case "$HAVE" in *"${NL}CORE_SYNC_NO_SESSION$NL"*) [ $rc -eq 0 ] && return 1 ;; esac
  case "$HAVE" in *"${NL}CORE_SYNC_END$NL"*) [ $rc -eq 0 ] && return 0 ;; esac
  PROBE_RC=$rc
  return 2
}
probe_failed() {  # $1 = probe's return
  if [ "$1" = 1 ]; then
    ERR="no such session directory on HIVE -- check --session"
  else
    ERR="could not list the session on HIVE (hive_exec.sh exit ${PROBE_RC:-?}); nothing was copied"
  fi
  echo "failed $S ($ERR)"
  FAIL_J="$(js .)"
  summary; exit 1
}
folder() { case "$1" in */*) printf '%s/' "${1%/*}" ;; esac; }

if [ "$MODE" = get ]; then
  ITEMS="$(echo $GET_OUTPUT $GET_RECORDS)"
  if [ "$DRY" = 1 ]; then
    for i in $ITEMS; do
      echo "would get $S/$i -> $L/$(folder "$i")"
      WOULD_J="$WOULD_J${WOULD_J:+, }$(js "$i")"
    done
    summary; exit 0
  fi
  mkdir -p "$L/output/search" "$L/input" || { ERR="could not create $L/output/search and $L/input"
    echo "failed $L ($ERR)"; FAIL_J="$(js .)"; summary; exit 1; }
  probe $ITEMS; p=$?; [ $p -eq 0 ] || probe_failed $p
  for i in $ITEMS; do
    case "$HAVE" in
      *"${NL}CORE_SYNC_HAVE $i$NL"*)
        if bash "$HIVE_EXEC" --get "$S/$i" "$L/$(folder "$i")" </dev/null; then
          echo "got $S/$i"; DONE_J="$DONE_J${DONE_J:+, }$(js "$i")"
        else
          rc=$?; echo "failed $S/$i (hive_exec.sh --get exit $rc)"; FAIL_J="$FAIL_J${FAIL_J:+, }$(js "$i")"
        fi ;;
      *) echo "missing $S/$i"; MISS_J="$MISS_J${MISS_J:+, }$(js "$i")" ;;
    esac
  done
else
  [ -d "$L/output" ] || die "no output/ folder under --local '$L' -- give the session folder --get filled"
  ANY=0
  for i in $PUT_OUTPUT; do
    if [ -e "$L/$i" ]; then
      ANY=1
      if [ "$DRY" = 1 ]; then
        echo "would put $L/$i -> $S/$(folder "$i")"; WOULD_J="$WOULD_J${WOULD_J:+, }$(js "$i")"
      fi
    fi
  done
  if [ "$DRY" != 1 ] && [ "$ANY" = 1 ]; then probe; p=$?; [ $p -eq 0 ] || probe_failed $p; fi
  for i in $PUT_OUTPUT; do
    if [ ! -e "$L/$i" ]; then
      echo "skipped $L/$i (not here)"; MISS_J="$MISS_J${MISS_J:+, }$(js "$i")"
    elif [ "$DRY" != 1 ]; then
      if bash "$HIVE_EXEC" --put "$L/$i" "$S/$(folder "$i")" </dev/null; then
        echo "put $L/$i"; DONE_J="$DONE_J${DONE_J:+, }$(js "$i")"
      else
        rc=$?; echo "failed $L/$i (hive_exec.sh --put exit $rc)"; FAIL_J="$FAIL_J${FAIL_J:+, }$(js "$i")"
      fi
    fi
  done
  [ "$ANY" = 1 ] || echo "core_sync.sh: none of the report files is under $L/output -- nothing was pushed" >&2
fi
summary
[ -z "$FAIL_J" ]
