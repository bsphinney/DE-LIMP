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
#   --get first makes one listing call on HIVE, which says whether the session is there and
#   which items exist: an absent one (AUDIT.md before step 8c, no conditions.csv, ...) is
#   "missing", not an error. Only then does it make <local>/output and <local>/input (a
#   mistyped --session leaves nothing behind here) and copy each GET item from
#   <session>/<item> into <local>/<the item's folder>/. The same call asks
#   `session.py which-search --paths` for the search the DE read -- its search_provenance.json
#   and report, wherever the search went (output/search_sage/ ...), and the workflow manifest
#   (2.14.2, T44: one answer to "which search"). With records and none chosen it pulls the
#   rest, prints STOP and exits 1: choose one (SKILL.md step 12b), then pull again -- a report
#   written without them cannot name the engine. A HIVE copy of the skill from before 2.14.2
#   has no which-search (exit 2): the STOP line says to run step 0 first.
#   --put copies each PUT item that exists under <local> into <session>/output/; the rest are
#   "skipped (not here)". Both refuse a --session that is not a directory on HIVE.
#
#   --session  the HIVE session dir ($S): absolute (/quobyte/...) or ~/... . Any name a lab
#              folder can have is accepted -- spaces, ( ) & ' and the rest (2.14.2 V1; review
#              M5): $S goes to HIVE's shell quoted by printf %q, as hive_exec.sh's rq() quotes
#              it, so nothing in it is interpreted there. Only a newline is refused.
#   --local    this computer's copy of the session (a leading ~/ is expanded).
#
# Output: one line per item -- "got <session>/<item>", "missing <session>/<item>",
# "put <local>/<item>", "skipped <local>/<item> (not here)", "failed <path> (...)" -- then
# one JSON line: {"mode": "get"|"put", "session": ..., "local": ..., "done": [...],
# "missing": [...], "failed": [...]}, items relative to the session ("." = the session
# itself, with an "error"); --dry-run adds "dry_run": true and "would": [...].
# Exit: 0 every item moved or was absent; 1 a transfer failed, HIVE could not be reached, the
# session is not there, no search is chosen (--get: which-search exit 3) or which-search
# failed, or --put found none of its files under --local; 2 a usage error.
#
# Transport: hive_exec.sh beside this file (HIVE_EXEC overrides it, as in report_issue.sh):
# rsync, or scp where there is none (Git Bash on Windows). Bash 3.2 and Git Bash safe.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HIVE_EXEC="${HIVE_EXEC:-$HERE/hive_exec.sh}"

# --get, from <session>/<item> into <local>/<folder>/ -- what the local report reads:
# make_analysis_html.py (step 9) the DE tables, figures, AUDIT.md and SAMPLE_QUALITY.md;
# analysis_prompt.py the search report and methods.md. The search's own files come between
# the two lists, from which-search on HIVE (SESSION_PY), never a fixed output/search path.
GET_OUTPUT="output/tables output/figures"
GET_OUTPUT_AFTER="output/AUDIT.md output/SAMPLE_QUALITY.md output/methods.md"
SESSION_PY="~/proteomics-pipeline/scripts/session.py"    # on HIVE; ~ expands there
# --get: the session records the report's sections read -- make_analysis_html.py through
# submission_report.py (session.json, the attached submission, input/locate_decisions.json),
# the conditions and the FASTA's metadata, and input/sample_type.json and experiment_type.json
# (sample_type.py / experiment_type.py, read by analysis_prompt.py and session_docs.py),
# input/raw_names.tsv (the names as acquired and as delivered: submission_report.stage_map reads
# it when no stage record is above the copy) and logs/core_client.staff.json (a Core client work
# mark counts only with it). This is the one pull list (tests/test_step9_pull.py pins it).
GET_RECORDS="session.json input/submission.json input/samples.tsv input/raw_files.txt
             input/conditions.csv input/search.fasta.meta.json input/sample_type.json
             input/experiment_type.json input/locate_decisions.json input/raw_names.tsv
             logs/core_client.staff.json"
# --put, from <local>/<item> into <session>/output/ when it exists here -- what deliver
# (core_submission.py DELIVER_FILES and the podcast) copies from $S. The PDF exists only where
# a headless browser printed it; the podcast only when it was made on this computer.
# report_provenance.json is the report's record (make_analysis_html.py writes it beside the
# report, make_podcast.py link again): without it finalize, deliver and the podcast on HIVE
# take the report for one made before 2.14.2.
PUT_OUTPUT="output/AI_Analysis_Report.md output/Analysis_Report.html output/Analysis_Report.md
            output/Analysis_Report.pdf output/Analysis_Report_with_audio.html output/methods.md
            output/methods.docx output/report_provenance.json output/podcast"

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

# --session goes into a command line on HIVE quoted by rq() (printf %q), so any character is
# safe there; a newline alone is refused -- it would split this script's one-item-per-line
# lists and HIVE's marker lines.
session_chars_ok() { case "$1" in *"$NL"*) return 1 ;; esac; return 0; }
# a session-relative path which-search names: it goes to hive_exec.sh, which quotes it the same
# way. Refused: empty, absolute, a ".." (outside the session), a newline.
rel_path_ok() { case "$1" in ""|/*|*..*|*"$NL"*) return 1 ;; esac; return 0; }
while [ "${#S}" -gt 1 ] && [ "${S%/}" != "$S" ]; do S="${S%/}"; done
case "$S" in
  /?*|"~/"?*) ;;
  *) die "--session must be the session's absolute path on HIVE (/quobyte/...) or ~/...: '$S'" ;;
esac
session_chars_ok "$S" || die "--session holds a newline: give the session's path on one line"

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
SEARCH_DRY="(the search: session.py which-search --paths)"
summary() {
  printf '{"mode": %s, "session": %s, "local": %s, "done": [%s], "missing": [%s], "failed": [%s]' \
    "$(js "$MODE")" "$(js "$S")" "$(js "$L")" "$DONE_J" "$MISS_J" "$FAIL_J"
  if [ "$DRY" = 1 ]; then printf ', "dry_run": true, "would": [%s]' "$WOULD_J"; fi
  if [ -n "$ERR" ]; then printf ', "error": %s' "$(js "$ERR")"; fi
  printf '}\n'
}

# The remote form of $S, one shell word: printf %q, with a leading ~/ as "$HOME"/ so HIVE's
# shell still expands it -- exactly hive_exec.sh's rq(), which quotes the same path for the
# transfers.
rq() { case "$1" in "~/"*) printf '"$HOME"/%s' "$(printf '%q' "${1#\~/}")" ;; *) printf '%q' "$1" ;; esac; }

# ONE call on HIVE: is $S a directory, and which of the items "$@" exist in it. Tells "not
# there" from a transport failure, which hive_exec.sh's exit codes cannot (rsync and scp both
# exit non-zero for a missing source, scp with the same 1 as for any other failure).
# Sets HAVE (marker lines); returns 0 = listed, 1 = no such session, 2 = HIVE not reached.
# With WITH_SEARCH=1 (--get) the same call asks which-search for the search's files: one
# "CORE_SYNC_SEARCH <path>" line each and "CORE_SYNC_SEARCH_RC <exit>" (its reason, on exit 3,
# goes to HIVE's stderr, which hive_exec.sh passes through).
probe() {
  local cmd out rc
  cmd="cd $(rq "$S") 2>/dev/null || { echo CORE_SYNC_NO_SESSION; exit 0; }"
  [ $# -eq 0 ] || cmd="$cmd; for p in $*; do [ -e \"\$p\" ] && echo \"CORE_SYNC_HAVE \$p\"; done"
  if [ "${WITH_SEARCH:-0}" = 1 ]; then
    cmd="$cmd; P=\$(python3 $SESSION_PY which-search --dir \"\$PWD\" --paths); echo \"CORE_SYNC_SEARCH_RC \$?\""
    cmd="$cmd; printf '%s\\n' \"\$P\" | sed 's/^/CORE_SYNC_SEARCH /'"
  fi
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
  LISTED="$(echo $GET_OUTPUT $GET_OUTPUT_AFTER $GET_RECORDS)"
  if [ "$DRY" = 1 ]; then
    for i in $GET_OUTPUT "$SEARCH_DRY" $GET_OUTPUT_AFTER $GET_RECORDS; do
      if [ "$i" = "$SEARCH_DRY" ]; then
        echo "would get the search the DE read (on HIVE: session.py which-search --paths) -> $L/<its folder>/"
      else
        echo "would get $S/$i -> $L/$(folder "$i")"
      fi
      WOULD_J="$WOULD_J${WOULD_J:+, }$(js "$i")"
    done
    summary; exit 0
  fi
  WITH_SEARCH=1 probe $LISTED; p=$?; [ $p -eq 0 ] || probe_failed $p
  # only now, with the session found: a mistyped --session leaves no empty tree here
  mkdir -p "$L/output" "$L/input" || { ERR="could not create $L/output and $L/input"
    echo "failed $L ($ERR)"; FAIL_J="$(js .)"; summary; exit 1; }
  # the search's files, as which-search named them (each exists on HIVE: it lists only those)
  SEARCH_RC="$(printf '%s\n' "$HAVE" | sed -n 's/^CORE_SYNC_SEARCH_RC \([0-9][0-9]*\)$/\1/p' | head -n 1)"
  SEARCH_ITEMS="" SEARCH_STOP=""
  if [ "$SEARCH_RC" = 0 ]; then
    while IFS= read -r f; do
      [ -n "$f" ] || continue
      if rel_path_ok "$f"; then
        SEARCH_ITEMS="$SEARCH_ITEMS$f$NL"; HAVE="$HAVE${NL}CORE_SYNC_HAVE $f$NL"
      else
        echo "failed $S/$f (which-search named a path this script will not pass to HIVE)"
        FAIL_J="$FAIL_J${FAIL_J:+, }$(js "$f")"
      fi
    done <<EOF_SEARCH
$(printf '%s\n' "$HAVE" | sed -n 's/^CORE_SYNC_SEARCH \(..*\)$/\1/p')
EOF_SEARCH
  elif [ "$SEARCH_RC" = 3 ]; then
    SEARCH_STOP="no search is chosen (session.py which-search exit 3; its reason is above) -- choose one (SKILL.md step 12b), then pull again"
  elif [ "$SEARCH_RC" = 2 ]; then
    # argparse's exit for an unknown command: HIVE's ~/proteomics-pipeline is older than 2.14.2
    SEARCH_STOP="session.py which-search failed on HIVE (exit 2); the search's files were not pulled -- HIVE's skill copy predates 2.14.2: run step 0 (\`skill_version.sh --session-start --mode hive_remote\`), then pull again"
  else
    SEARCH_STOP="session.py which-search failed on HIVE (exit ${SEARCH_RC:-unknown}); the search's files were not pulled"
  fi
  # one item per line: a search folder's name can hold a space (the fixed lists hold none)
  ITEMS="$(printf '%s\n' $GET_OUTPUT)$NL$SEARCH_ITEMS$(printf '%s\n' $GET_OUTPUT_AFTER $GET_RECORDS)"
  while IFS= read -r i; do
    [ -n "$i" ] || continue
    case "$HAVE" in
      *"${NL}CORE_SYNC_HAVE $i$NL"*)
        mkdir -p "$L/$(folder "$i")"
        if bash "$HIVE_EXEC" --get "$S/$i" "$L/$(folder "$i")" </dev/null; then
          echo "got $S/$i"; DONE_J="$DONE_J${DONE_J:+, }$(js "$i")"
        else
          rc=$?; echo "failed $S/$i (hive_exec.sh --get exit $rc)"; FAIL_J="$FAIL_J${FAIL_J:+, }$(js "$i")"
        fi ;;
      *) echo "missing $S/$i"; MISS_J="$MISS_J${MISS_J:+, }$(js "$i")" ;;
    esac
  done <<EOF_ITEMS
$ITEMS
EOF_ITEMS
  if [ -n "$SEARCH_STOP" ]; then
    echo "STOP: $SEARCH_STOP"
    ERR="$SEARCH_STOP"; FAIL_J="$FAIL_J${FAIL_J:+, }$(js "which-search")"
  fi
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
  if [ "$ANY" != 1 ]; then
    # nothing to push is a mistake in the step, not a success: deliver copies from $S
    ERR="none of the report files is under $L/output -- nothing was pushed"
    echo "core_sync.sh: $ERR" >&2
    summary; exit 1
  fi
fi
summary
[ -z "$FAIL_J" ]
