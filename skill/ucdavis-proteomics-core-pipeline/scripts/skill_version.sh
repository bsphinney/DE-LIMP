#!/usr/bin/env bash
# =============================================================================
# skill_version.sh  --  skill_version.py for bash, and the check that HIVE runs the same skill
# as this computer, and the Core's current release.
#
# THE VERSION. skill_version.py is the ONE reader of .claude-plugin/plugin.json beside scripts/.
# Bash cannot import it, and on Windows `python3` is often the Microsoft Store stub, so this is
# its bash mirror -- the only one: report_issue.sh sources it, and tests/test_skill_version.py
# keeps it equal to the original, fixture for fixture.
#   bash skill_version.sh                 # this copy's version, or the UNKNOWN tag
#   . skill_version.sh                    # skill_version_of / skill_version / skill_label /
#                                         # skill_version_core / skill_version_cmp /
#                                         # skill_installed_version / skill_update_state
#
# DRIFT (SKILL.md step 0, hive_remote). On 2026-09-29 two Core staff were still on 2.6.0, on
# their laptops AND in ~/proteomics-pipeline on HIVE, two releases behind: none of the 2.7/2.8
# fixes had ever run for them. One HIVE copy had no .claude-plugin/ at all, so its sdrf.tsv said
# "v0.0.0". Nothing compared the copies, and nothing said a newer release existed.
#   bash skill_version.sh --check-hive [--mode hive_remote]
#   One SSH call; JSON on stdout; the exit code says what to do:
#     0  HIVE runs this skill. Relay `say` if set: this computer is behind the Core's release.
#        Also 0, silently, with no HIVE login here -- unless --mode hive_remote, where the
#        login is expected: then 5.
#     3  ~/proteomics-pipeline on HIVE is missing, has no .claude-plugin/, is older, or is a
#        different build (its files differ) -> run `next` (hive_exec.sh --put-skill), then
#        check again.
#     4  HIVE's copy is NEWER than this computer's -> stop and relay `say`. Never put the older
#        one over it unasked.
#     5  could not check (no HIVE login in hive_remote, HIVE unreachable, no plugin.json here)
#        -> say so in one line and carry on. The check is never fatal.
#     6  as 3, but SLURM jobs of this user whose batch scripts use that copy are running or
#        queued -- their job-end steps run its scripts hours later -- or squeue answered with an
#        error, so they could not be counted. It is not replaced under them: put it after they
#        finish, or when the user says to.
#
# THE CURRENT RELEASE is $SKILL_RELEASE_DIR/CURRENT_VERSION on HIVE (default
# /quobyte/proteomics-grp/skill_release, beside skill_runs/ and skill_issues/): the version on
# its first line, comments after it. Only the Core group can read it; for anyone else the
# release is not compared and nothing is said. It is written by a maintainer as the last step
# of a release, once the new version is on GitHub main -- never by --put-skill, which also
# carries unreleased test builds:
#   bash skill_version.sh --publish-release [--allow-older]
# /quobyte/proteomics-grp is group-writable and not sticky: any member can rename or replace
# skill_release/. So the file is trusted only while it AND its folder are owned by a Core admin
# (core_admins.txt beside this file, the one list); anything else reads `untrusted`, as if
# nothing were published. Only an admin can publish.
#
# HIVE_EXEC (default: hive_exec.sh beside this file) is the transport, as in report_issue.sh;
# SKILL_SQUEUE (default squeue) is what lists the jobs on HIVE, SKILL_OWNER_CMD (default GNU
# `stat -c %U`, BSD `stat -f %Su` as a fallback) who owns a path there. Tests point all three at
# stand-ins.
# No HIVE login here (hive.env, or HIVE_USER + HIVE_KEY) = no SSH at all.
#
# BEHIND GITHUB MAIN (SKILL.md step 0, every mode). On 2026-10-07 a Windows laptop (Git Bash, no
# GitHub SSH key) ran `claude plugin update ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core`.
# It warned "marketplace not refreshed: SSH host key is not in your known_hosts", then "SSH
# authentication failed", and said "already at the latest version (2.10.0)" while main shipped
# 2.11.2: it had compared against its own stale catalogue. The catalogue (the marketplace clone,
# ~/.claude/plugins/marketplaces/ucdavis-proteomics-core) had an HTTPS remote, and
# `git -C <clone> pull --ff-only` followed by the same update took it to 2.11.2.
#   bash skill_version.sh --check-update     # read-only: this copy vs main's plugin.json
#     0  this copy is what main ships (or newer: a test build). Say nothing.
#     3  behind main -> the update question, then `next` (--update). status reload_needed: the
#        new version is installed already, this session still runs the old one (/reload-plugins).
#     5  main could not be read (offline, no curl) -> one line, carry on. Never "up to date".
#   bash skill_version.sh --update           # does the update, and checks what landed
#     0  updated (the installed plugin is now main's version), or nothing was needed
#     1  the update did NOT happen: `say` says why and gives the by-hand route
#     5  it ran, but could not be checked against main (main unreadable, or the installed
#        version could not be read). Not a confirmation.
#   main's version is read over HTTPS (curl, which Git Bash has; SKILL_MAIN_URL), never SSH, and
#   no Python is needed. --update refreshes a catalogue behind main with
#   `git pull --ff-only` over HTTPS (an SSH remote is pulled from SKILL_REPO_URL instead), runs
#   `claude plugin update` (SKILL_CLAUDE) with CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1, and reads the
#   installed version from `claude plugin list --json`, or installed_plugins.json beside the
#   catalogue. SKILL_MARKETPLACE_DIR overrides where the catalogue is. Tests point all four at
#   stand-ins (a file:// URL, local git repositories, a fake claude); no test uses the network.
# =============================================================================

# The tag for a version that cannot be read -- skill_version.py's UNKNOWN, character for
# character (never a guessed number: CLAUDE.md rule 2).
SKILL_VERSION_UNKNOWN="(unknown — plugin.json not found)"
SV_HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The version in the plugin.json at $1, or SKILL_VERSION_UNKNOWN: the first "version": "<text>"
# in the file -- skill_version.py's answer on every fixture of tests/test_skill_version.py.
skill_version_of() {
  local v
  v="$(sed -nE 's/.*"version"[[:space:]]*:[[:space:]]*"([^"]*)".*/\1/p' "$1" 2>/dev/null \
       | head -n1 | tr -d '[:space:]')"
  if [ -n "$v" ]; then printf '%s\n' "$v"; else printf '%s\n' "$SKILL_VERSION_UNKNOWN"; fi
}

# The version of the skill whose scripts/ folder is $1 (default: this file's).
skill_version() { skill_version_of "${1:-$SV_HERE}/../.claude-plugin/plugin.json"; }

# How a version reads after the skill's name: "v2.8.0", or the UNKNOWN tag as it is.
skill_label() {
  if [ "$1" = "$SKILL_VERSION_UNKNOWN" ]; then printf '%s\n' "$1"; else printf 'v%s\n' "$1"; fi
}

# The release a version belongs to: without a leading "v" and without a pre-release or build
# suffix ("v2.9.0-dev+3" -> "2.9.0"). Anything else as given.
skill_version_core() {
  local v="${1#[vV]}"
  v="${v%%[-+]*}"
  printf '%s\n' "$v"
}

# -1, 0 or 1 as version $1 is older than, the same as, or newer than $2, number by number
# ("2.10.0" is newer than "2.9.0"; "2.9" is "2.9.0"; a number of any length). Prints nothing and
# returns 1 when either is not plain numbers and dots: a suffix (compare skill_version_core's),
# a blank, the UNKNOWN tag.
skill_version_cmp() {
  local a="$1" b="$2" x y v
  for v in "$a" "$b"; do
    case "$v" in ''|.*|*.|*..*|*[!0-9.]*) return 1 ;; esac
  done
  while [ -n "$a$b" ]; do
    x="${a%%.*}"; y="${b%%.*}"
    case "$a" in *.*) a="${a#*.}" ;; *) a="" ;; esac
    case "$b" in *.*) b="${b#*.}" ;; *) b="" ;; esac
    # as text, not as integers: `[ -lt ]` overflows past 19 digits
    x="${x#"${x%%[!0]*}"}"; y="${y#"${y%%[!0]*}"}"
    if [ "${#x}" -ne "${#y}" ]; then
      if [ "${#x}" -lt "${#y}" ]; then echo -1; else echo 1; fi
      return 0
    fi
    if [[ "$x" < "$y" ]]; then echo -1; return 0; fi
    if [[ "$x" > "$y" ]]; then echo 1; return 0; fi
  done
  echo 0
}

# One checksum of a skill copy's files: `cksum` of each of the named files under $1, in the
# order given, then of that list. The same on both sides (POSIX cksum, verified identical on
# macOS and HIVE). A missing file changes it.
skill_files_digest() {
  local root="$1"; shift
  (cd "$root" 2>/dev/null && cksum -- "$@" 2>/dev/null | cksum | tr -s ' \t' ' ')
}

# The Proteomics Core admins, one per line: core_admins.txt ($1, default beside this file)
# without `#` comments, blank lines or surrounding space (a CRLF file reads the same). Nothing
# when the file is missing -- then nobody is an admin, and nothing is trusted.
skill_core_admins() {
  sed -e 's/#.*//' -e 's/[[:space:]]*$//' -e 's/^[[:space:]]*//' \
      "${1:-$SV_HERE/core_admins.txt}" 2>/dev/null | tr -d '\r' | grep -v '^$'
}

# True when $1 is a Core admin.
skill_is_core_admin() {
  [ -n "${1:-}" ] && skill_core_admins "${2:-}" | grep -qxF -- "$1"
}

# How `claude plugin` names this plugin and its marketplace.
SKILL_PLUGIN_ID="ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core"
SKILL_MARKETPLACE="ucdavis-proteomics-core"

# The LOWEST version that the JSON text on stdin gives the plugin $1 (default SKILL_PLUGIN_ID):
# every object that names it, read to the end of that object. That covers
# `claude plugin list --json` ({"id": "<id>", "version": ...}) and installed_plugins.json
# ("<id>": [{..., "version": ...}]), as Claude Code 2.1.295 writes them, pretty-printed or on one
# line. The lowest, so that an install still behind never reads as current. Only plain versions
# count. Prints nothing when no object names the plugin.
skill_installed_version() {
  local key="\"${1:-$SKILL_PLUGIN_ID}\"" rb='}' text chunk v low=""
  text="$(LC_ALL=C tr -d '\r\n')"
  while :; do
    case "$text" in *"$key"*) ;; *) break ;; esac
    text="${text#*"$key"}"
    chunk="${text%%"$rb"*}"
    v="$(printf '%s\n' "$chunk" \
         | LC_ALL=C sed -nE 's/.*"version"[[:space:]]*:[[:space:]]*"([^"]*)".*/\1/p' | head -n1)"
    skill_version_cmp "$(skill_version_core "$v")" 0 >/dev/null || continue
    if [ -z "$low" ] || [ "$(skill_version_cmp "$(skill_version_core "$v")" \
                                               "$(skill_version_core "$low")")" = -1 ]; then
      low="$v"
    fi
  done
  [ -z "$low" ] || printf '%s\n' "$low"
}

# Where this computer stands against GitHub main, for the running copy $1, main's version $2
# and the installed plugin's $3 (may be blank):
#   current        $1 is main's version, or newer (a test build)
#   reload_needed  $1 is behind, but the installed plugin already is main's: this session still
#                  runs the old copy
#   behind         $1 is behind main
#   unknown        $1 or $2 is not a plain version, so nothing can be said
skill_update_state() {
  local here main inst c
  here="$(skill_version_core "$1")"; main="$(skill_version_core "$2")"
  inst="$(skill_version_core "${3:-}")"
  if ! c="$(skill_version_cmp "$here" "$main")"; then echo unknown; return 0; fi
  if [ "$c" != -1 ]; then echo current; return 0; fi
  if c="$(skill_version_cmp "$inst" "$main")" && [ "$c" != -1 ]; then
    echo reload_needed; return 0
  fi
  echo behind
}

# ---- run as a command (not sourced) ------------------------------------------------------
if [ "${BASH_SOURCE[0]}" = "$0" ]; then
  set -uo pipefail
  SKILL_RELEASE_DIR="${SKILL_RELEASE_DIR:-/quobyte/proteomics-grp/skill_release}"
  HIVE_EXEC="${HIVE_EXEC:-$SV_HERE/hive_exec.sh}"
  CFG="${HIVE_ENV_FILE:-$HOME/.config/ucdavis-proteomics/hive.env}"
  if [ -z "${HIVE_USER:-}" ] && [ -f "$CFG" ]; then . "$CFG"; fi
  HIVE_COPY="~/proteomics-pipeline"
  PUT_SKILL="bash scripts/hive_exec.sh --put-skill"
  # Checked on Claude Code 2.1.285: `/plugin marketplace update` refreshes the catalogue only;
  # the installed plugin stays where it was. `claude plugin update` moves it -- but only as far
  # as its catalogue, which it may fail to refresh (over SSH, 2026-10-07), so --update does the
  # refresh over HTTPS itself and checks the version that landed.
  UPDATE_CMD="bash scripts/skill_version.sh --update"
  UPDATE_HOW="$UPDATE_CMD, which refreshes the plugin catalogue over HTTPS when it is behind, runs claude plugin update $SKILL_PLUGIN_ID and checks the version that landed; then /reload-plugins. To stop this recurring: /plugin -> Marketplaces -> $SKILL_MARKETPLACE -> Enable auto-update"
  SKILL_MAIN_URL="${SKILL_MAIN_URL:-https://raw.githubusercontent.com/bsphinney/DE-LIMP/main/skill/ucdavis-proteomics-core-pipeline/.claude-plugin/plugin.json}"
  SKILL_REPO_URL="${SKILL_REPO_URL:-https://github.com/bsphinney/DE-LIMP.git}"
  SKILL_CLAUDE="${SKILL_CLAUDE:-claude}"
  ROOT="$(cd "$SV_HERE/.." && pwd)"

  sv_die() { echo "skill_version.sh: $*" >&2; exit 2; }
  sv_have_login() { [ -n "${HIVE_USER:-}" ] && [ -n "${HIVE_KEY:-}" ] && [ -f "$HIVE_EXEC" ]; }
  # A JSON string, or null for an empty value. Terminal escapes and every other control
  # character are dropped: an ssh error line in colour once made the JSON unparseable.
  sv_js() {
    local s="${1:-}"
    s="$(printf '%s' "$s" | LC_ALL=C sed "s/$(printf '\033')\[[0-9;]*[A-Za-z]//g" \
         | LC_ALL=C tr '\t\n\r' '   ' | LC_ALL=C tr -d '\000-\037\177')"
    [ -n "$s" ] || { printf 'null'; return; }
    s="${s//\\/\\\\}"; s="${s//\"/\\\"}"
    printf '"%s"' "$s"
  }
  sv_is_version() { skill_version_cmp "$1" "$1" >/dev/null; }

  # The files the digest covers, relative to the skill root: plugin.json and every file in this
  # copy's scripts/ (not bytecode). HIVE hashes the same names -- files an older put left there
  # and this copy no longer has are not counted.
  sv_files() {
    local f
    printf '%s\n' ".claude-plugin/plugin.json"
    for f in "$ROOT"/scripts/*; do
      [ -f "$f" ] || continue
      case "${f##*/}" in *.pyc|.*) continue ;; esac
      printf 'scripts/%s\n' "${f##*/}"
    done
  }

  # What HIVE has, in ONE command: the copy's folders, version and digest (read by this
  # file's functions, shipped along -- a copy older than this file has none of its own), the
  # release, and this user's SLURM jobs: all of them, and those whose batch script names the
  # copy's scripts/ in one of the ways a script can spell it. A script that cannot be read
  # counts as one that does (fail closed); a job with no script path (an interactive srun)
  # does not. Not the bare word "proteomics-pipeline": that also matches ~/.proteomics-pipeline
  # (the toolchain) and ~/proteomics-pipeline-280-*. SKILL_VERSION_UNKNOWN is blank there: a
  # copy with no readable version reads "".
  sv_remote_read() {
    local f
    printf 'SKILL_VERSION_UNKNOWN=\nD=%q\nQ=%q\nOWN=%q\nF=(' "$SKILL_RELEASE_DIR" \
      "${SKILL_SQUEUE:-squeue}" "${SKILL_OWNER_CMD:-}"
    while IFS= read -r f; do printf ' %q' "$f"; done <<EOF
$(sv_files)
EOF
    printf ' )\n'
    declare -f skill_version_of skill_files_digest
    cat <<'SH'
P="$HOME/proteomics-pipeline"
if [ -d "$P/scripts" ]; then echo "SKILLCHECK scripts=yes"; else echo "SKILLCHECK scripts=no"; fi
if [ -f "$P/.claude-plugin/plugin.json" ]; then echo "SKILLCHECK plugin=yes"; else echo "SKILLCHECK plugin=no"; fi
echo "SKILLCHECK version=$(skill_version_of "$P/.claude-plugin/plugin.json")"
echo "SKILLCHECK digest=$(skill_files_digest "$P" "${F[@]}")"
if ! command -v "$Q" >/dev/null 2>&1; then echo "SKILLCHECK jobs=absent"
else
  EF="$(mktemp 2>/dev/null || echo /dev/null)"
  if J="$("$Q" -h -u "$(id -un)" -o '%o' 2>"$EF")"; then
    PP="$(cd "$P" 2>/dev/null && pwd -P)"
    PAT=(-e "$P/scripts/" -e '~/proteomics-pipeline/scripts/'
         -e '$HOME/proteomics-pipeline/scripts/' -e '${HOME}/proteomics-pipeline/scripts/')
    [ -n "$PP" ] && PAT+=(-e "$PP/scripts/")
    U=0
    while read -r N S; do
      case "$S" in /*) ;; *) continue ;; esac
      if [ ! -r "$S" ] || grep -qF "${PAT[@]}" -- "$S" 2>/dev/null; then U=$((U + N)); fi
    done <<LIST
$(printf '%s\n' "$J" | grep . | sort | uniq -c)
LIST
    echo "SKILLCHECK jobs=$U"
    echo "SKILLCHECK jobs_total=$(printf '%s\n' "$J" | grep -c .)"
  else
    echo "SKILLCHECK jobs=error"
    echo "SKILLCHECK jobs_error=$(LC_ALL=C grep -av '^[[:space:]]*$' "$EF" | head -n1 | cut -c1-200)"
  fi
  [ "$EF" = /dev/null ] || rm -f "$EF"
fi
G="$(dirname "$D")"
# who owns a path -- the link itself, never what it points to (no -L): a member's symlink to an
# admin's file is still the member's
own() { if [ -n "$OWN" ]; then "$OWN" "$1"
        else stat -c %U -- "$1" 2>/dev/null || stat -f %Su -- "$1" 2>/dev/null; fi; }
echo "SKILLCHECK me=$(id -un)"
if [ -e "$D" ] || [ -L "$D" ]; then echo "SKILLCHECK release_dir_owner=$(own "$D")"; fi
if [ -f "$D/CURRENT_VERSION" ] && [ -r "$D/CURRENT_VERSION" ]; then
  echo "SKILLCHECK release=$(sed -e '/^[[:space:]]*#/d' -e '/^[[:space:]]*$/d' "$D/CURRENT_VERSION" | head -n1 | tr -d '[:space:]')"
  echo "SKILLCHECK release_owner=$(own "$D/CURRENT_VERSION")"
  echo "SKILLCHECK release_state=published"
elif [ -r "$G" ] && [ -x "$G" ]; then echo "SKILLCHECK release_state=not_published"
else echo "SKILLCHECK release_state=no_access"; fi
echo "SKILLCHECK end=1"
SH
  }

  # Run sv_remote_read on HIVE into OUT; false when the answer did not come back whole.
  # LC_ALL=C: under a UTF-8 locale macOS `tr` stops at the first byte that is not UTF-8
  # ("Illegal byte sequence") and an ssh error in another encoding was cut short. GNU grep (Linux,
  # Git Bash) in a UTF-8 locale prints nothing for such a line, only "binary file matches" on
  # stderr, so every grep over HIVE's answer runs as `LC_ALL=C grep -a`.
  sv_read_hive() {
    OUT="$(bash "$HIVE_EXEC" "$(sv_remote_read)" 2>&1 | LC_ALL=C tr -d '\r')"
    printf '%s\n' "$OUT" | LC_ALL=C grep -aq '^SKILLCHECK end=1$'
  }
  sv_field() { printf '%s\n' "$OUT" | sed -n "s/^SKILLCHECK $1=//p" | head -n1; }
  # The release as read from CURRENT_VERSION: a plain version ("v2.9.0" reads 2.9.0), else "".
  # Nothing unless the file and its folder are both owned by a Core admin (sv_release_trusted).
  sv_release() {
    local r; r="$(sv_field release)"; r="${r#[vV]}"
    if sv_release_trusted && sv_is_version "$r"; then printf '%s\n' "$r"; fi
  }
  sv_release_trusted() {
    skill_is_core_admin "$(sv_field release_owner)" \
      && skill_is_core_admin "$(sv_field release_dir_owner)"
  }

  # Text from HIVE's side (an ssh or squeue error) as plain ASCII: a non-UTF-8 byte there
  # would make the JSON unreadable.
  sv_ascii() { LC_ALL=C tr -d '\200-\377' | cut -c1-200; }

  sv_emit() {   # status local hive release release_state behind_release jobs next say
    # (+ JOBS_TOTAL, the user's jobs of any kind, set by sv_check_hive)
    printf '{"status": %s, "local": %s, "hive": %s, "hive_path": %s, "release": %s, ' \
      "$(sv_js "$1")" "$(sv_js "$2")" "$(sv_js "$3")" "$(sv_js "$HIVE_COPY")" "$(sv_js "$4")"
    printf '"release_state": %s, "behind_release": %s, "jobs": %s, "jobs_total": %s, ' \
      "$(sv_js "$5")" "${6:-null}" "${7:-null}" "${JOBS_TOTAL:-null}"
    printf '"next": %s, "say": %s}\n' "$(sv_js "$8")" "$(sv_js "$9")"
  }

  sv_check_hive() {
    local mode="$1" here hive="" rel="" rel_state="not_checked" behind=null status code
    local next="" say="" c why jobs=null jobs_err="" d_here d_hive hc lc f files=()
    here="$(skill_version)"
    if ! sv_have_login; then
      if [ "$mode" = hive_remote ]; then
        sv_emit no_login "$here" "" "" not_checked null null "" \
          "No HIVE login is saved on this computer, so the skill on HIVE was not compared with this one ($here). Run scripts/check_access.sh with the HIVE user and key (it saves them to ~/.config/ucdavis-proteomics/hive.env), then this check again."
        return 5
      fi
      sv_emit skipped "$here" "" "" not_checked null null "" ""   # local mode, or no HIVE
      return 0
    fi
    if [ "$here" = "$SKILL_VERSION_UNKNOWN" ]; then
      sv_emit local_unknown "" "" "" not_checked null null "" \
        "This copy of the skill has no readable .claude-plugin/plugin.json, so its version is unknown and cannot be compared with HIVE's. Reinstall it from the plugin marketplace (references/install.md)."
      return 5
    fi
    if ! sv_read_hive; then
      why="$(printf '%s\n' "$OUT" | LC_ALL=C grep -av '^SKILLCHECK ' | LC_ALL=C grep -av '^[[:space:]]*$' \
             | tail -n1 | sv_ascii)"
      sv_emit unreachable "$here" "" "" not_checked null null "" \
        "Could not reach HIVE to compare skill versions${why:+ ($why)}; carrying on with $here."
      return 5
    fi

    hive="$(sv_field version)"
    rel_state="$(sv_field release_state)"; rel_state="${rel_state:-no_access}"
    rel="$(sv_release)"
    if [ "$rel_state" = published ] && ! sv_release_trusted; then
      rel_state=untrusted                    # not an admin's: as if nothing were published
    elif [ "$rel_state" = published ] && [ -z "$rel" ]; then
      rel_state=unreadable
    fi
    lc="$(skill_version_core "$here")"
    if [ -n "$rel" ] && c="$(skill_version_cmp "$lc" "$rel")"; then
      if [ "$c" = -1 ]; then behind=true; else behind=false; fi
    fi
    case "$(sv_field jobs)" in
      absent) ;;                                          # no squeue there: nothing to count
      error)  jobs_err="$(sv_field jobs_error | sv_ascii)"; jobs_err="${jobs_err:-no message}" ;;
      ''|*[!0-9]*) jobs_err="squeue gave no count" ;;
      *) jobs="$(sv_field jobs)"
         JOBS_TOTAL="$(sv_field jobs_total)"
         case "$JOBS_TOTAL" in ''|*[!0-9]*) JOBS_TOTAL=null ;; esac ;;
    esac
    # an array, not $(sv_files) unquoted: a name with a space ("session 2.py", an iCloud
    # duplicate) split in two, and a file only one side had went unnoticed
    while IFS= read -r f; do files+=("$f"); done <<LIST
$(sv_files)
LIST
    d_here="$(skill_files_digest "$ROOT" ${files[@]+"${files[@]}"})"
    d_hive="$(sv_field digest)"

    hc="$(skill_version_core "$hive")"
    if [ "$(sv_field scripts)" != yes ] || [ "$(sv_field plugin)" != yes ] || [ -z "$hive" ]; then
      status=hive_missing; code=3
      say="HIVE has no complete copy of the skill in $HIVE_COPY (no scripts/, or no readable .claude-plugin/plugin.json, so every record written there would name the version 'unknown')."
    elif c="$(skill_version_cmp "$hc" "$lc")" && [ "$c" = 1 ]; then
      status=hive_ahead; code=4
      case "$behind" in
        true)  say="HIVE has skill $hive in $HIVE_COPY, newer than this computer's $here, and the Core's current release is $rel. Update the plugin on this computer ($UPDATE_HOW), then start again." ;;
        false) say="HIVE has skill $hive in $HIVE_COPY, newer than this computer's $here -- a build not released yet (the Core's current release is $rel), so there is nothing to update to here. Ask the user whether to put $here back on HIVE ($PUT_SKILL) or to leave $hive there." ;;
        *)     say="HIVE has skill $hive in $HIVE_COPY, newer than this computer's $here. If a newer release is out, update the plugin on this computer ($UPDATE_HOW); if $hive is a test build, ask the user whether to put $here back ($PUT_SKILL)." ;;
      esac
    elif [ "$c" = -1 ]; then
      status=hive_behind; code=3
      say="HIVE runs skill $hive from $HIVE_COPY; this computer has $here."
    elif [ "$hive" != "$here" ] || { [ -n "$d_here" ] && [ -n "$d_hive" ] && [ "$d_here" != "$d_hive" ]; }; then
      # the same release, but another build of it (a test build, or files changed by hand):
      # this computer's wins -- its SKILL.md is the one running
      status=hive_differs; code=3
      say="HIVE's copy in $HIVE_COPY ($hive) is not the same build as this computer's ($here): its files differ."
    else
      status=in_step; code=0
    fi

    if [ "$code" = 3 ]; then
      if [ -n "$jobs_err" ]; then
        code=6
        say="$say It is not replaced now: squeue on HIVE answered with an error ($jobs_err), so your SLURM jobs that run scripts from $HIVE_COPY could not be counted. Put this computer's $here there once squeue answers and none are running ($PUT_SKILL), or now if the user says so."
      elif [ "$jobs" != null ] && [ "$jobs" -gt 0 ]; then
        code=6
        say="$say It is not replaced now: $jobs of your ${JOBS_TOTAL/null/?} SLURM job(s) on HIVE use $HIVE_COPY (their batch scripts name its scripts/, or cannot be read), and their end-of-job steps would run the new scripts. Put this computer's $here there after they finish ($PUT_SKILL), or now if the user says so."
      else
        next="$PUT_SKILL"
        say="$say Putting this computer's $here there."
      fi
    fi
    if [ "$behind" = true ] && [ "$status" != hive_ahead ]; then
      say="${say:+$say }This computer has skill $here, but the Core's current release is $rel: the fixes since $here do not run until you update. Update the plugin ($UPDATE_HOW), then start again."
    fi
    sv_emit "$status" "$here" "$hive" "$rel" "$rel_state" "$behind" "$jobs" "$next" "$say"
    return "$code"
  }

  # Write the release file on HIVE, atomically (a reader sees the old file or the new one).
  sv_remote_publish() {
    printf 'D=%q\nV=%q\n' "$SKILL_RELEASE_DIR" "$1"
    cat <<'SH'
G="$(dirname "$D")"
if [ ! -d "$D" ]; then
  if [ -d "$G" ] && [ -w "$G" ] && mkdir "$D" 2>/dev/null; then chmod 2770 "$D" 2>/dev/null
  else echo "SKILLPUB error=cannot create $D (a Proteomics Core account is needed)"; exit 0; fi
fi
[ -w "$D" ] || { echo "SKILLPUB error=cannot write $D (a Proteomics Core account is needed)"; exit 0; }
T="$(mktemp "$D/.CURRENT_VERSION.XXXXXX" 2>/dev/null)" || { echo "SKILLPUB error=cannot write in $D"; exit 0; }
if { printf '%s\n' "$V"
     printf '# The current release of the ucdavis-proteomics-core-pipeline skill (UC Davis Proteomics\n'
     printf '# Core), written by %s on %s with skill_version.sh --publish-release.\n' \
       "$(id -un)" "$(date '+%Y-%m-%d %H:%M %Z')"
     printf '# A laptop older than line 1 is told to update (SKILL.md step 0).\n'
   } > "$T" && chmod 664 "$T" && mv -f "$T" "$D/CURRENT_VERSION"; then
  echo "SKILLPUB ok=1"
else rm -f "$T"; echo "SKILLPUB error=could not write $D/CURRENT_VERSION"; fi
SH
  }

  # The version origin/main ships, as this checkout last fetched it (no network). Nothing, and
  # status 1, when this copy is not a git checkout that has origin/main: the installed plugin.
  sv_main_version() {
    local prefix tmp v=""
    git -C "$ROOT" rev-parse --is-inside-work-tree >/dev/null 2>&1 || return 1
    prefix="$(git -C "$ROOT" rev-parse --show-prefix 2>/dev/null)" || return 1
    tmp="$(mktemp "${TMPDIR:-/tmp}/skill_version_main.XXXXXX")" || return 1
    if git -C "$ROOT" show "origin/main:${prefix}.claude-plugin/plugin.json" >"$tmp" 2>/dev/null; then
      v="$(skill_version_of "$tmp")"
    fi
    rm -f "$tmp"
    [ -n "$v" ] && [ "$v" != "$SKILL_VERSION_UNKNOWN" ] || return 1
    printf '%s\n' "$v"
  }

  sv_publish() {
    local allow_older="$1" here cur c err main main_checked=false me downer
    here="$(skill_version)"; here="${here#[vV]}"
    sv_is_version "$here" || sv_die "this copy's version is '$here' -- only a released version (numbers and dots, from .claude-plugin/plugin.json) can be published"
    if main="$(sv_main_version)"; then
      main_checked=true
      [ "${main#[vV]}" = "$here" ] || sv_die "origin/main (as this checkout last fetched it) ships $main, not $here: publish only a version main ships (git fetch origin, then again)"
    fi
    sv_have_login || sv_die "no HIVE login here (hive.env, or HIVE_USER + HIVE_KEY): CURRENT_VERSION lives on HIVE"
    sv_read_hive || sv_die "could not reach HIVE: $(printf '%s\n' "$OUT" | tail -n1)"
    me="$(sv_field me)"
    skill_is_core_admin "$me" \
      || sv_die "HIVE user '${me:-?}' is not a Proteomics Core admin (scripts/core_admins.txt): only an admin publishes the release"
    downer="$(sv_field release_dir_owner)"
    if [ -n "$downer" ] && ! skill_is_core_admin "$downer"; then
      sv_die "$SKILL_RELEASE_DIR is owned by '$downer', not a Core admin, so what is written there would not be trusted: move it aside and publish again (it is then made anew)"
    fi
    cur="$(sv_release)"
    if [ -n "$cur" ] && c="$(skill_version_cmp "$here" "$cur")" && [ "$c" = -1 ] && [ "$allow_older" != 1 ]; then
      sv_die "the published release is $cur and this copy is $here, older. Going back to it (a withdrawn release) needs --allow-older."
    fi
    OUT="$(bash "$HIVE_EXEC" "$(sv_remote_publish "$here")" 2>&1 | LC_ALL=C tr -d '\r')"
    if ! printf '%s\n' "$OUT" | LC_ALL=C grep -aq '^SKILLPUB ok=1$'; then
      err="$(printf '%s\n' "$OUT" | sed -n 's/^SKILLPUB error=//p' | head -n1)"
      sv_die "not published: ${err:-$(printf '%s\n' "$OUT" | tail -n1)}"
    fi
    printf '{"published": %s, "previous": %s, "file": %s, "main_checked": %s}\n' \
      "$(sv_js "$here")" "$(sv_js "$cur")" "$(sv_js "$SKILL_RELEASE_DIR/CURRENT_VERSION")" \
      "$main_checked"
  }

  # ---- behind GitHub main: --check-update / --update ---------------------------------------

  # The marketplace clone (the plugin catalogue `claude plugin update` reads): beside the
  # plugin cache this copy runs from, else under CLAUDE_CONFIG_DIR or ~/.claude.
  sv_marketplace_dir() {
    if [ -n "${SKILL_MARKETPLACE_DIR:-}" ]; then printf '%s\n' "$SKILL_MARKETPLACE_DIR"; return; fi
    case "$ROOT" in
      */plugins/cache/"$SKILL_MARKETPLACE"/*)
        printf '%s\n' "${ROOT%%/plugins/cache/*}/plugins/marketplaces/$SKILL_MARKETPLACE"; return ;;
    esac
    printf '%s\n' "${CLAUDE_CONFIG_DIR:-$HOME/.claude}/plugins/marketplaces/$SKILL_MARKETPLACE"
  }
  # A path as a user pastes it: "$HOME/..." when under the home folder (spaces in a Windows
  # user name survive), else shell-quoted.
  sv_show_path() {
    case "$1" in
      "$HOME"/*) printf '"$HOME/%s"' "${1#"$HOME"/}" ;;
      *) printf '%q' "$1" ;;
    esac
  }
  # The catalogue's version of the plugin, or "" (no clone, or no readable plugin.json in it).
  sv_catalogue_version() {
    local v
    v="$(skill_version_of "$1/skill/ucdavis-proteomics-core-pipeline/.claude-plugin/plugin.json")"
    [ "$v" = "$SKILL_VERSION_UNKNOWN" ] || printf '%s\n' "$v"
  }
  # PULL_TAIL: what `git pull --ff-only` pulls from. Nothing (the clone's own remote) unless that
  # remote is SSH; then the repository's HTTPS URL and main. The clone Claude Code made on
  # 2026-10-07 had an HTTPS remote, and a plain pull was the route that worked.
  sv_pull_tail() {
    PULL_TAIL=()
    case "$(git -C "$1" remote get-url origin 2>/dev/null)" in
      http://*|https://*|'') ;;
      ssh://*|git@*|*@*:*) PULL_TAIL=("$SKILL_REPO_URL" main) ;;
    esac
  }
  # The update by hand, as one line to paste in Git Bash or a terminal.
  sv_manual() {
    local t=""
    sv_pull_tail "$1"
    [ "${#PULL_TAIL[@]}" -eq 0 ] || t=" ${PULL_TAIL[*]}"
    printf 'git -C %s pull --ff-only%s && claude plugin update %s' \
      "$(sv_show_path "$1")" "$t" "$SKILL_PLUGIN_ID"
  }

  # MAIN: the version GitHub main ships, from its plugin.json over HTTPS. Else "" and MAIN_ERR.
  sv_fetch_main() {
    local tmp err
    MAIN=""; MAIN_ERR=""
    if ! command -v curl >/dev/null 2>&1; then MAIN_ERR="curl is not on PATH"; return 1; fi
    tmp="$(mktemp "${TMPDIR:-/tmp}/skill_main.XXXXXX" 2>/dev/null)" \
      || { MAIN_ERR="no temporary file could be made"; return 1; }
    if err="$(curl -fsSL --max-time 20 --retry 1 -o "$tmp" "$SKILL_MAIN_URL" 2>&1)"; then
      MAIN="$(skill_version_of "$tmp")"
      sv_is_version "$(skill_version_core "$MAIN")" \
        || { MAIN=""; MAIN_ERR="main's plugin.json has no readable version"; }
    else
      MAIN_ERR="$(printf '%s\n' "$err" | LC_ALL=C grep -av '^[[:space:]]*$' | tail -n1 | sv_ascii)"
      MAIN_ERR="${MAIN_ERR:-curl failed}"
    fi
    rm -f "$tmp"
    [ -n "$MAIN" ]
  }

  # The installed plugin's version: `claude plugin list --json`, else installed_plugins.json
  # beside the catalogue (as Claude Code 2.1.295 writes it). "" when neither says.
  sv_installed() {
    local v="" f
    if command -v "$SKILL_CLAUDE" >/dev/null 2>&1; then
      v="$("$SKILL_CLAUDE" plugin list --json </dev/null 2>/dev/null | skill_installed_version)"
    fi
    f="$(dirname "$(dirname "$1")")/installed_plugins.json"
    if [ -z "$v" ] && [ -f "$f" ]; then v="$(skill_installed_version <"$f")"; fi
    printf '%s\n' "$v"
  }

  # Refresh the catalogue: fast-forward only, never a reset (a clone with local changes fails
  # and says so). No prompt can hang it. REFRESH = ok | failed | no_clone, REFRESH_ERR = why.
  sv_refresh() {
    local out
    REFRESH_ERR=""
    if [ ! -e "$1/.git" ]; then
      REFRESH=no_clone; REFRESH_ERR="there is no plugin catalogue (git clone) at $1"; return 1
    fi
    if ! command -v git >/dev/null 2>&1; then
      REFRESH=failed; REFRESH_ERR="git is not on PATH"; return 1
    fi
    sv_pull_tail "$1"
    if out="$(GIT_TERMINAL_PROMPT=0 git -c http.lowSpeedLimit=1000 -c http.lowSpeedTime=60 \
                -C "$1" pull --ff-only ${PULL_TAIL[@]+"${PULL_TAIL[@]}"} </dev/null 2>&1)"; then
      REFRESH=ok; return 0
    fi
    REFRESH=failed
    REFRESH_ERR="$(printf '%s\n' "$out" | LC_ALL=C grep -aE '^(fatal|error):' | head -n1 | sv_ascii)"
    [ -n "$REFRESH_ERR" ] || REFRESH_ERR="$(printf '%s\n' "$out" | LC_ALL=C grep -av '^[[:space:]]*$' \
                                            | tail -n1 | sv_ascii)"
    REFRESH_ERR="${REFRESH_ERR:-git pull failed}"
    return 1
  }

  # `claude plugin update`, with Claude Code's documented switch for a machine without a GitHub
  # SSH key (CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1: an owner/repo marketplace is cloned over HTTPS,
  # no SSH probe). CLAUDE_RC is "" when claude is not on PATH; CLAUDE_SAID its output, on one line.
  sv_claude_update() {
    local out
    CLAUDE_RC=""; CLAUDE_SAID=""
    command -v "$SKILL_CLAUDE" >/dev/null 2>&1 || return 1
    out="$(CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1 "$SKILL_CLAUDE" plugin update "$SKILL_PLUGIN_ID" \
           </dev/null 2>&1)"
    CLAUDE_RC=$?
    CLAUDE_SAID="$(printf '%s\n' "$out" | LC_ALL=C tr -d '\r' | LC_ALL=C tr -d '\200-\377' \
                   | awk 'NF { printf "%s%s", (n++ ? " / " : ""), $0 }' | cut -c1-500)"
  }

  sv_check_update() {
    local here inst="" state say="" next="" behind=null status code mp
    MAIN=""; MAIN_ERR=""
    here="$(skill_version)"
    mp="$(sv_marketplace_dir)"
    if [ "$here" = "$SKILL_VERSION_UNKNOWN" ]; then
      status=local_unknown; code=5; here=""
      say="This copy of the skill has no readable .claude-plugin/plugin.json, so it cannot be compared with GitHub main. Reinstall it from the plugin marketplace (references/install.md)."
    elif ! sv_fetch_main; then
      status=main_unknown; code=5
      say="Could not read the version GitHub main ships ($MAIN_ERR), so this computer's skill ($here) was not compared with it. That is not a confirmation that $here is current; carrying on."
    else
      state="$(skill_update_state "$here" "$MAIN")"
      if [ "$state" = behind ]; then
        inst="$(sv_installed "$mp")"
        state="$(skill_update_state "$here" "$MAIN" "$inst")"
      fi
      case "$state" in
        current) status=current; code=0; behind=false ;;
        reload_needed)
          status=reload_needed; code=3; behind=true
          say="Skill $inst is installed on this computer, but this session still runs $here: type /reload-plugins (or restart Claude Code), then start again." ;;
        behind)
          status=behind; code=3; behind=true; next="$UPDATE_CMD"
          say="This computer has skill $here; GitHub main ships $MAIN, so the fixes since $here do not run here yet. Updating takes a minute: $UPDATE_HOW." ;;
        *) status=main_unknown; code=5
           say="Could not compare this computer's skill ($here) with GitHub main ($MAIN). That is not a confirmation that it is current; carrying on." ;;
      esac
    fi
    printf '{"status": %s, "local": %s, "main": %s, "installed": %s, "behind_main": %s, ' \
      "$(sv_js "$status")" "$(sv_js "$here")" "$(sv_js "$MAIN")" "$(sv_js "$inst")" "$behind"
    printf '"next": %s, "manual": %s, "say": %s}\n' \
      "$(sv_js "$next")" "$(sv_js "$(sv_manual "$mp")")" "$(sv_js "$say")"
    return "$code"
  }

  sv_update() {
    local here mp inst0="" inst1="" cat0 cat1 state=unknown status code say="" why="" c lc
    here="$(skill_version)"; [ "$here" != "$SKILL_VERSION_UNKNOWN" ] || here=""
    mp="$(sv_marketplace_dir)"
    sv_fetch_main
    cat0="$(sv_catalogue_version "$mp")"
    REFRESH=not_needed; REFRESH_ERR=""; CLAUDE_RC=""; CLAUDE_SAID=""
    # the installed version only when this copy is not already main's (claude takes seconds)
    if [ -n "$MAIN" ] && [ -n "$here" ]; then state="$(skill_update_state "$here" "$MAIN")"; fi
    if [ "$state" != current ]; then
      inst0="$(sv_installed "$mp")"
      [ -z "$MAIN" ] || state="$(skill_update_state "${here:-$inst0}" "$MAIN" "$inst0")"
    fi
    case "$state" in
      current)
        status=current; code=0; inst1="$inst0"
        say="This computer's skill (${here:-$inst0}) is what GitHub main ships ($MAIN): nothing to update." ;;
      reload_needed)
        status=reload_needed; code=0; inst1="$inst0"
        say="Skill $inst0 is installed on this computer, but this session still runs ${here:-an older copy}: type /reload-plugins (or restart Claude Code), then start again." ;;
      *)
        # the catalogue first, over HTTPS, unless it already has main's version
        if [ -n "$MAIN" ] && [ -n "$cat0" ] && c="$(skill_version_cmp "$(skill_version_core "$cat0")" \
             "$(skill_version_core "$MAIN")")" && [ "$c" != -1 ]; then
          REFRESH=not_needed
        else
          sv_refresh "$mp"
        fi
        sv_claude_update
        inst1="$(sv_installed "$mp")"
        # Why, in words staff can be told. claude's own "already at the latest version" is never
        # repeated here (it is in claude_said): it compared against a stale catalogue.
        lc="$(printf '%s' "$CLAUDE_SAID" | LC_ALL=C tr '[:upper:]' '[:lower:]')"
        cat1="$(sv_catalogue_version "$mp")"
        if [ -z "$CLAUDE_RC" ]; then
          why="The claude command is not on PATH here, so claude plugin update could not run (in Claude Code: /plugin -> Installed -> ucdavis-proteomics-core-pipeline -> Update now)."
        elif [ "$CLAUDE_RC" != 0 ]; then
          why="claude plugin update failed (exit $CLAUDE_RC)."
        else
          case "$lc" in
            *"not refreshed"*) why="claude plugin update could not refresh its own catalogue, so it compared against the old one (${cat1:-version unknown})." ;;
            *"already at the latest"*) why="claude plugin update compared against its own catalogue, which is at ${cat1:-an unknown version}." ;;
          esac
        fi
        case "$REFRESH" in
          failed|no_clone) why="Refreshing the plugin catalogue over HTTPS failed: ${REFRESH_ERR%.}.${why:+ $why}" ;;
        esac
        if [ -z "$MAIN" ]; then
          status=unverified; code=5
          say="The update ran, but the version GitHub main ships could not be read ($MAIN_ERR), so it is not confirmed: the installed skill is now ${inst1:-unknown}. That is not a confirmation that it is current.${why:+ $why}"
        elif [ -z "$inst1" ]; then
          status=unverified; code=5
          say="The update ran, but the installed version could not be read (claude plugin list --json), so it is not confirmed. Check with claude plugin list that ucdavis-proteomics-core-pipeline shows $MAIN.${why:+ $why}"
        elif c="$(skill_version_cmp "$(skill_version_core "$inst1")" "$(skill_version_core "$MAIN")")" \
             && [ "$c" != -1 ]; then
          status=updated; code=0
          say="Updated the skill from ${inst0:-${here:-an older version}} to $inst1, the version GitHub main ships. Type /reload-plugins (or restart Claude Code), then start again: this session still runs ${here:-the old copy}."
        else
          status=not_updated; code=1
          say="The update did NOT happen: this computer still has skill $inst1, and GitHub main ships $MAIN.${why:+ $why} To do it by hand, in Git Bash or a terminal: $(sv_manual "$mp"), then /reload-plugins."
        fi ;;
    esac
    cat1="$(sv_catalogue_version "$mp")"
    printf '{"status": %s, "local": %s, "main": %s, "installed_before": %s, "installed": %s, ' \
      "$(sv_js "$status")" "$(sv_js "$here")" "$(sv_js "$MAIN")" "$(sv_js "$inst0")" "$(sv_js "$inst1")"
    printf '"catalogue_before": %s, "catalogue": %s, "refresh": %s, "refresh_error": %s, ' \
      "$(sv_js "$cat0")" "$(sv_js "$cat1")" "$(sv_js "$REFRESH")" "$(sv_js "$REFRESH_ERR")"
    printf '"claude_exit": %s, "claude_said": %s, "manual": %s, "say": %s}\n' \
      "${CLAUDE_RC:-null}" "$(sv_js "$CLAUDE_SAID")" "$(sv_js "$(sv_manual "$mp")")" "$(sv_js "$say")"
    return "$code"
  }

  case "${1:-}" in
    "")                [ $# -eq 0 ] || sv_die "usage"; skill_version ;;
    --check-hive)      case "$*" in
                         "--check-hive")                  sv_check_hive "" ;;
                         "--check-hive --mode hive_remote") sv_check_hive hive_remote ;;
                         "--check-hive --mode local")     sv_check_hive local ;;
                         *) sv_die "usage: --check-hive [--mode hive_remote|local]" ;;
                       esac
                       exit $? ;;
    --check-update)    [ $# -eq 1 ] || sv_die "usage: --check-update"
                       sv_check_update; exit $? ;;
    --update)          [ $# -eq 1 ] || sv_die "usage: --update"
                       sv_update; exit $? ;;
    --publish-release) case "${2:-}" in
                         "") sv_publish 0 ;;
                         --allow-older) [ $# -eq 2 ] || sv_die "usage"; sv_publish 1 ;;
                         *) sv_die "unknown argument: $2" ;;
                       esac ;;
    -h|--help)         sed -n '3,77p' "$0" ;;
    *)                 sv_die "unknown argument: $1 (see --help)" ;;
  esac
fi
