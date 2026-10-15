#!/usr/bin/env bash
# =============================================================================
# coreomics_key_to_hive.sh  --  Put a staff member's CoreOmics API key on HIVE, once, so the
# CoreOmics steps (core_submission.py check / identify / fetch / bioshare / email-draft) can run
# there. It is for hive_remote from a laptop with no usable Python: on Windows `python3` is
# usually the Microsoft Store stub (check_access.sh: local_python3.usable false), so those steps
# cannot run on the laptop, and HIVE had no key. Michelle (msalemi, 2026-10-08) saved her key with
# the skill's own line and her Claude still could not read the submission; Gabriela (gabrig,
# 2026-10-02) first reported it.
#
#   bash coreomics_key_to_hive.sh             THE AGENT runs this. The key the user already saved
#                                             on this computer goes to ~/.coreomics_token on HIVE.
#   bash coreomics_key_to_hive.sh --paste     No key on this computer: THE USER runs this in their
#                                             own Git Bash window (or terminal), never the chat. It
#                                             waits, showing nothing, for the key, then sends it.
#   bash coreomics_key_to_hive.sh --status    Is the key on HIVE, and private? (never the key)
#   add --remove-local                        to the first form: delete this computer's copy once
#                                             HIVE has it. The default keeps it: it is harmless, and
#                                             save_transcript.py blanks that key out of saved
#                                             conversations only when it can read it here.
#
# The key is never printed, never on a command line (the process list, the shell history,
# commands.log), and never copied into a file on this computer: the saved key file itself is
# ssh's standard input, and a pasted key goes from `read -rs` (silent) through printf, a shell
# builtin, into the same pipe. On HIVE it is written into a NEW file made under umask 077 (mode
# 600) beside ~/.coreomics_token and renamed over it. HIVE home folders are drwxrwsr-x -- every
# account can list them and group members can add entries -- so the file's own mode is the
# key's only protection, and a rename replaces anything planted at that name (a link to someone
# else's file) instead of writing the key through it. Before that rename HIVE checks the new
# file with core_submission.py's rule for a usable key; one that cannot be a key (only blanks,
# "Authorization: Token ...") is refused and the key already there kept.
#
# One JSON object on stdout; notes on stderr. Exit 0 done (--status: on HIVE and private);
# 2 nothing to send, a key HIVE refuses, or --status: absent / readable by others; 3 HIVE
# unreachable, the key could not be checked there, or the write failed. The transport is hive_exec.sh beside this file (HIVE_EXEC overrides it, as in
# report_issue.sh), with the login check_access.sh saved.
#
# check_access.sh sources this file for COREOMICS_KEY_PROBE, coreomics_key_state,
# coreomics_local_key, coreomics_route and coreomics_advice -- the one definition of where the
# CoreOmics steps run and what to do about the key. Sourcing it runs nothing.
# =============================================================================
set -uo pipefail
CK_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COREOMICS_KEY_NAME=".coreomics_token"
COREOMICS_KEY_MAX_BYTES=4096           # a key is ~40 characters; anything this big is not one
COREOMICS_CHECK_ON_HIVE="bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/core_submission.py check --json'"

# --- what runs ON HIVE (no single quote inside either: check_access.sh splices the probe into
# --- its own single-quoted `bash -l -c '...'`) -----------------------------------------------
# One line saying what ~/.coreomics_token is. Symlinks are followed, as Python's os.stat does.
COREOMICS_KEY_PROBE='f="$HOME/.coreomics_token"; if [ ! -e "$f" ]; then echo COREOMICS_KEY=absent; elif [ ! -f "$f" ]; then echo COREOMICS_KEY=notfile; elif [ ! -O "$f" ]; then echo COREOMICS_KEY=notmine; elif [ ! -r "$f" ]; then echo COREOMICS_KEY=unreadable; elif [ ! -s "$f" ]; then echo COREOMICS_KEY=empty; else echo "COREOMICS_KEY=mode:$(stat -L -c %a "$f" 2>/dev/null || stat -L -f %Lp "$f" 2>/dev/null)"; fi'
# Is the new file $t a usable key? core_submission.py decides -- read_key_file, then
# key_shape_problem: the one definition of a usable key -- from the skill on HIVE, BEFORE the
# rename, so a blank or malformed laptop file ("Authorization: Token ...") never replaces a
# working key (2.12 review). python3 -I: a module a group member leaves in this (group-writable)
# home folder is never imported. The answer is a fixed reason, never any part of the key.
COREOMICS_KEY_SHAPE='s="$HOME/proteomics-pipeline/scripts"; if [ ! -f "$s/core_submission.py" ]; then echo COREOMICS_KEY_STORE=uncheckable:noskill; exit 3; fi; v="$(python3 -I -c "import sys; sys.path.insert(0, sys.argv[1]); import core_submission as c; k = c.read_key_file(sys.argv[2]); w = c.key_shape_problem(k) if k else \"it is blank: only spaces or line breaks\"; print(\"notakey:\" + w if w else \"keyok\")" "$s" "$t" 2>/dev/null)"; case "$v" in keyok) ;; notakey:*) echo "COREOMICS_KEY_STORE=$v"; exit 2 ;; *) echo COREOMICS_KEY_STORE=uncheckable:python; exit 3 ;; esac; '
# Standard input -> a new mode-600 file -> checked -> renamed over ~/.coreomics_token; then the probe.
COREOMICS_KEY_STORE='umask 077; t="$(mktemp "$HOME/.coreomics_token.XXXXXX")" || { echo COREOMICS_KEY_STORE=failed:mktemp; exit 3; }; trap "rm -f -- \"\$t\"" EXIT; cat > "$t" || { echo COREOMICS_KEY_STORE=failed:write; exit 3; }; if [ ! -s "$t" ]; then echo COREOMICS_KEY_STORE=empty; exit 2; fi; '"$COREOMICS_KEY_SHAPE"'chmod 600 "$t" && mv -f -- "$t" "$HOME/.coreomics_token" || { echo COREOMICS_KEY_STORE=failed:rename; exit 3; }; echo COREOMICS_KEY_STORE=ok; '"$COREOMICS_KEY_PROBE"

ck_js() {  # a JSON string literal
  local s="$1"
  s="${s//\\/\\\\}"; s="${s//\"/\\\"}"; s="${s//$'\t'/\\t}"; s="${s//$'\r'/\\r}"; s="${s//$'\n'/\\n}"
  printf '"%s"' "$s"
}

# The files core_submission.key_files() reads on this computer, in its order: on Windows
# Python's ~ (the profile folder, USERPROFILE) first, then Git Bash's $HOME when that differs.
coreomics_local_key_files() {
  local prof=""
  case "$(uname -s 2>/dev/null)" in
    MINGW*|MSYS*|CYGWIN*)
      if [ -n "${USERPROFILE:-}" ]; then
        if command -v cygpath >/dev/null 2>&1; then prof="$(cygpath -u "$USERPROFILE" 2>/dev/null)"; fi
        if [ -z "$prof" ]; then
          case "$USERPROFILE" in
            [A-Za-z]:[\\/]*) prof="/$(printf '%s' "${USERPROFILE%%:*}" | tr '[:upper:]' '[:lower:]')/${USERPROFILE:3}"
                             prof="${prof//\\//}" ;;
            *) prof="$USERPROFILE" ;;
          esac
        fi
        prof="${prof%/}"
        printf '%s\n' "$prof/$COREOMICS_KEY_NAME"
      fi ;;
  esac
  if [ -n "${HOME:-}" ] && [ "${HOME%/}" != "$prof" ]; then printf '%s\n' "${HOME%/}/$COREOMICS_KEY_NAME"; fi
}

# The first of those that holds something (a regular, non-empty file); exit 1 when none does.
coreomics_local_key() {
  local f
  while IFS= read -r f; do
    if [ -n "$f" ] && [ -f "$f" ] && [ -s "$f" ]; then printf '%s\n' "$f"; return 0; fi
  done <<EOF
$(coreomics_local_key_files)
EOF
  return 1
}

# The fix for a key file on HIVE that is not private: chmod -- unless another account owns it,
# where chmod cannot work and the file is removed and the key sent again (core_submission.py's
# key_mode_fix says the same on HIVE).
COREOMICS_KEY_CHMOD_FIX="bash scripts/hive_exec.sh 'chmod 600 ~/.coreomics_token'"
COREOMICS_KEY_NOTMINE_FIX="it belongs to another account, so chmod cannot fix it: remove it (bash scripts/hive_exec.sh 'rm -f ~/.coreomics_token') and send the key again (bash scripts/coreomics_key_to_hive.sh, or the user runs it with --paste in their own Git Bash window)"

# coreomics_key_state <text holding the probe's line>: sets CK_STATE (true | false | mode-wrong |
# null = not checked), CK_MODE (octal, or ""), CK_DETAIL (plain words) and CK_FIX (for
# mode-wrong: one of the two fixes above; otherwise "").
coreomics_key_state() {
  local line m
  line="$(printf '%s\n' "$1" | LC_ALL=C tr -d '\r' | grep -m1 '^COREOMICS_KEY=' || true)"
  CK_MODE=""; CK_FIX=""
  case "$line" in
    "") CK_STATE=null; CK_DETAIL="HIVE was not asked (no SSH login, or no answer)" ;;
    COREOMICS_KEY=absent) CK_STATE=false; CK_DETAIL="there is no ~/.coreomics_token on HIVE" ;;
    COREOMICS_KEY=empty) CK_STATE=false; CK_DETAIL="~/.coreomics_token on HIVE is empty" ;;
    COREOMICS_KEY=notfile) CK_STATE=false; CK_DETAIL="~/.coreomics_token on HIVE is not a file" ;;
    COREOMICS_KEY=unreadable) CK_STATE=false; CK_DETAIL="~/.coreomics_token on HIVE cannot be read, even by its owner" ;;
    COREOMICS_KEY=notmine) CK_STATE=mode-wrong; CK_FIX="$COREOMICS_KEY_NOTMINE_FIX"
                           CK_DETAIL="~/.coreomics_token on HIVE belongs to another account" ;;
    COREOMICS_KEY=mode:*)
      m="${line#COREOMICS_KEY=mode:}"; CK_MODE="$m"; CK_FIX="$COREOMICS_KEY_CHMOD_FIX"
      case "$m" in
        ""|*[!0-7]*) CK_STATE=mode-wrong; CK_DETAIL="the mode of ~/.coreomics_token on HIVE could not be read" ;;
        # No permission for group or others: the last two octal digits are 0 (600, 400, 700).
        *00) CK_STATE=true; CK_FIX=""; CK_DETAIL="~/.coreomics_token on HIVE, mode $m: only this account can read it" ;;
        *) CK_STATE=mode-wrong; CK_DETAIL="~/.coreomics_token on HIVE is mode $m: other HIVE accounts can read it" ;;
      esac ;;
    *) CK_STATE=null; CK_DETAIL="HIVE's answer was not understood" ;;
  esac
}

# coreomics_route <mode> <local python3 usable: true|false> <key on this computer: true|false>
#                 <key on HIVE: true|false|mode-wrong|null>
# Where core_submission.py check / identify / fetch / bioshare / email-draft run: "hive" (through
# hive_exec.sh, the key in ~/.coreomics_token there) or "this_computer" (the machine Claude Code
# runs on -- which IS HIVE in hive_local). Only hive_remote ever routes to HIVE: always without a
# usable local Python, and with one only when the key is on HIVE and not here.
coreomics_route() {
  if [ "$1" != hive_remote ]; then echo this_computer
  elif [ "$2" != true ]; then echo hive
  elif [ "$3" = true ]; then echo this_computer
  elif [ "$4" = true ] || [ "$4" = mode-wrong ]; then echo hive
  else echo this_computer; fi
}

# coreomics_advice <route> <key on this computer> <key on HIVE> -- one paragraph for the agent.
coreomics_advice() {
  local route="$1" here="$2" hive="$3" fix_mode="" out
  if [ "$hive" = mode-wrong ]; then
    fix_mode="The CoreOmics key on HIVE is not private: ${CK_FIX:-$COREOMICS_KEY_CHMOD_FIX} (core_submission.py refuses it until then). If other accounts could read it for long, the user Regenerates it in CoreOmics and it is put on HIVE again (bash scripts/coreomics_key_to_hive.sh). "
  fi
  if [ "$route" = hive ]; then
    case "$hive" in
      true|mode-wrong) out="Run the CoreOmics steps (core_submission.py check / identify / fetch / bioshare / email-draft) ON HIVE through hive_exec.sh; the key is ~/.coreomics_token there. First: $COREOMICS_CHECK_ON_HIVE" ;;
      *) if [ "$here" = true ]; then
           out="The CoreOmics steps run ON HIVE (no usable Python here), and HIVE has no key yet. Run bash scripts/coreomics_key_to_hive.sh -- it sends the key saved on this computer to HIVE over ssh standard input (never shown, never on a command line), mode 600 -- then: $COREOMICS_CHECK_ON_HIVE"
         else
           out="The CoreOmics steps run ON HIVE (no usable Python here), and there is no CoreOmics key on this computer or on HIVE. When a CoreOmics step comes up, the user makes a key in CoreOmics (Profile > API Key > Create) and runs, in their own Git Bash window, never this chat: bash '$CK_HERE/coreomics_key_to_hive.sh' --paste -- it asks for the key without showing it and saves it only on HIVE. Then: $COREOMICS_CHECK_ON_HIVE"
         fi ;;
    esac
  elif [ "$here" = true ]; then
    out="The CoreOmics steps run on this computer, where the key is: python3 scripts/core_submission.py check --json"
  else
    out="The CoreOmics steps run on this computer, which has no CoreOmics key yet: python3 scripts/core_submission.py check --json says how to save one."
  fi
  printf '%s' "$fix_mode$out"
}

ck_usage() {
  sed -n '3,/^# =====/p' "${BASH_SOURCE[0]}" | sed '$d' | sed 's/^# \{0,1\}//'
  exit "${1:-2}"
}

ck_fail() {  # <exit code> <status> <say> [<fix>]
  printf '{"status": %s, "coreomics_key_on_hive": null, "say": %s, "fix": %s}\n' \
    "$(ck_js "$2")" "$(ck_js "$3")" "$(ck_js "${4:-}")"
  exit "$1"
}

# The line that says why a transport call failed, never one of our own markers.
ck_why() {
  printf '%s\n' "$1" | LC_ALL=C tr -d '\r' | grep -v -e '^COREOMICS_KEY' -e '^Warning: Permanently added' \
    | grep -v '^[[:space:]]*$' | tail -1
}

ck_main() {
  local action=send remove_local=false hx out rc src="" from
  while [ $# -gt 0 ]; do
    case "$1" in
      --paste) action=paste ;;
      --status) action=status ;;
      --remove-local) remove_local=true ;;
      -h|--help) ck_usage 0 ;;
      *) echo "coreomics_key_to_hive.sh: unknown argument '$1'" >&2; ck_usage 2 ;;
    esac
    shift
  done
  if $remove_local && [ "$action" != send ]; then
    ck_fail 2 usage "--remove-local goes with the default form only (the key saved on this computer)."
  fi
  hx="${HIVE_EXEC:-$CK_HERE/hive_exec.sh}"
  [ -f "$hx" ] || ck_fail 3 no_transport "hive_exec.sh is missing beside this script ($hx)." \
    "Reinstall or update the skill, then try again."

  if [ "$action" = status ]; then
    out="$(bash "$hx" "$COREOMICS_KEY_PROBE" 2>&1 </dev/null)"; rc=$?
    coreomics_key_state "$out"
    if [ "$CK_STATE" = null ]; then
      ck_fail 3 unreachable "Could not ask HIVE about the key: $(ck_why "$out") (exit $rc)." \
        "Run check_access.sh to see why the SSH login fails."
    fi
    local fix="" status state_json
    case "$CK_STATE" in
      true) status=ok; state_json=true ;;
      false) status=absent; state_json=false
             fix="bash scripts/coreomics_key_to_hive.sh (a key saved on this computer), or the user runs it with --paste in their own Git Bash window" ;;
      *) status=mode_wrong; state_json='"mode-wrong"'
         fix="${CK_FIX:-$COREOMICS_KEY_CHMOD_FIX}; if other accounts could read it for long, Regenerate the key in CoreOmics and put the new one on HIVE" ;;
    esac
    printf '{"status": %s, "coreomics_key_on_hive": %s, "hive_mode": %s, "say": %s, "fix": %s}\n' \
      "$(ck_js "$status")" "$state_json" "$([ -n "$CK_MODE" ] && ck_js "$CK_MODE" || echo null)" \
      "$(ck_js "$CK_DETAIL")" "$(ck_js "$fix")"
    [ "$CK_STATE" = true ] && exit 0
    exit 2
  fi

  if [ "$action" = send ]; then
    src="$(coreomics_local_key)" || {
      local looked
      looked="$(coreomics_local_key_files | tr '\n' ' ')"
      ck_fail 2 no_key_here "There is no saved CoreOmics key on this computer (looked in: ${looked% })." \
        "The user makes a key in CoreOmics (Profile > API Key > Create) and runs, in their own Git Bash window (never the chat): bash '$CK_HERE/coreomics_key_to_hive.sh' --paste"
    }
    local size
    size="$(wc -c < "$src" | tr -d '[:space:]')"
    if [ "${size:-0}" -gt "$COREOMICS_KEY_MAX_BYTES" ]; then
      ck_fail 2 not_a_key "$src holds $size bytes, far more than a CoreOmics key; it was not sent." \
        "Save the key again on this computer (core_submission.py's line, or --paste)."
    fi
    out="$(bash "$hx" "$COREOMICS_KEY_STORE" 2>&1 < "$src")"; rc=$?
    from="$src"
  else
    # --paste: the user's own terminal. `read -s` shows nothing; printf is a builtin, so the key
    # is on no command line; nothing runs after Ctrl-D or an empty Enter, so a key already on
    # HIVE is never replaced with nothing.
    local TOK=""
    IFS= read -rs -p "Paste your CoreOmics API key (nothing will show), then press Enter: " TOK
    [ -t 0 ] && echo >&2
    TOK="${TOK#"${TOK%%[![:space:]]*}"}"; TOK="${TOK%"${TOK##*[![:space:]]}"}"
    if [ -z "$TOK" ]; then
      unset TOK
      ck_fail 2 nothing_pasted "No key was pasted, so nothing was sent (a key already on HIVE is unchanged)." \
        "The user runs this in their own Git Bash window (not through Claude) and pastes the key at the prompt: bash '$CK_HERE/coreomics_key_to_hive.sh' --paste"
    fi
    out="$(printf '%s\n' "$TOK" | bash "$hx" "$COREOMICS_KEY_STORE" 2>&1)"; rc=$?
    unset TOK
    from="the hidden prompt"
  fi

  coreomics_key_state "$out"
  if ! printf '%s\n' "$out" | LC_ALL=C tr -d '\r' | grep -q '^COREOMICS_KEY_STORE=ok$'; then
    if printf '%s\n' "$out" | grep -q '^COREOMICS_KEY_STORE=empty'; then
      ck_fail 2 empty "Nothing arrived on HIVE, so nothing was saved there (a key already on HIVE is unchanged)." \
        "Save the key again on this computer, or use --paste."
    fi
    # HIVE's key check (COREOMICS_KEY_SHAPE): its reason is fixed text, never part of the key
    local why fix
    why="$(printf '%s\n' "$out" | LC_ALL=C tr -d '\r' | sed -n 's/^COREOMICS_KEY_STORE=notakey://p' | head -1)"
    if [ -n "$why" ]; then
      if [ "$action" = send ]; then
        fix="Save the key again on this computer -- only the key, on one line (python3 scripts/core_submission.py check --json says how) -- or the user runs this with --paste in their own Git Bash window: bash '$CK_HERE/coreomics_key_to_hive.sh' --paste"
      else
        fix="The user runs it again in their own Git Bash window and pastes only the key (CoreOmics: Profile > API Key), nothing before or after it: bash '$CK_HERE/coreomics_key_to_hive.sh' --paste"
      fi
      ck_fail 2 not_a_key "The key from $from is not a usable CoreOmics key: $why. It was not saved on HIVE; a key already there is unchanged." "$fix"
    fi
    if printf '%s\n' "$out" | grep -q '^COREOMICS_KEY_STORE=uncheckable:noskill'; then
      ck_fail 3 not_checked "HIVE could not check the key, because the skill is not on HIVE (no ~/proteomics-pipeline/scripts/core_submission.py), so it was not saved; a key already there is unchanged." \
        "bash scripts/hive_exec.sh --put-skill, then run this again."
    fi
    if printf '%s\n' "$out" | grep -q '^COREOMICS_KEY_STORE=uncheckable:'; then
      ck_fail 3 not_checked "HIVE could not check the key: python3 there could not run the skill's key check (core_submission.py), so it was not saved; a key already there is unchanged." \
        "bash scripts/hive_exec.sh --put-skill (the skill on HIVE may be partial or older than this one), then run this again."
    fi
    ck_fail 3 not_saved "The key was not saved on HIVE: $(ck_why "$out") (exit $rc). A key already there is unchanged." \
      "Run check_access.sh to see whether the SSH login works, then try again."
  fi
  if [ "$CK_STATE" != true ]; then
    ck_fail 3 saved_not_private "The key reached HIVE but is not private there: $CK_DETAIL." \
      "${CK_FIX:-$COREOMICS_KEY_CHMOD_FIX}; then: bash scripts/coreomics_key_to_hive.sh --status"
  fi
  local copy="none"
  if [ "$action" = send ]; then
    copy="kept"
    if $remove_local; then rm -f -- "$src" && copy="removed" || copy="kept (could not remove it)"; fi
  fi
  printf '{"status": "ok", "coreomics_key_on_hive": true, "hive_file": "~/.coreomics_token", "hive_mode": %s, "sent_from": %s, "laptop_copy": %s, "say": %s, "next": %s}\n' \
    "$(ck_js "$CK_MODE")" "$(ck_js "$from")" "$(ck_js "$copy")" \
    "$(ck_js "The CoreOmics key is on HIVE as ~/.coreomics_token (mode $CK_MODE, readable only by you). The CoreOmics steps now run on HIVE.")" \
    "$(ck_js "$COREOMICS_CHECK_ON_HIVE")"
  exit 0
}

if [ "${BASH_SOURCE[0]}" = "$0" ]; then ck_main "$@"; fi
