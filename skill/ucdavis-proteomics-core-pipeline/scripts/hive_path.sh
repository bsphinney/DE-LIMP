#!/usr/bin/env bash
# =============================================================================
# hive_path.sh  --  Is this local path on a network share that HIVE already mounts?
# If so, print the HIVE path to use IN PLACE instead of uploading it.
#
#   bash hive_path.sh 'T:\Data\lab\service\PROT_0807'           # Windows mapped drive
#   bash hive_path.sh /Volumes/proteomics/Data/lab/service/x     # macOS SMB mount
#
# gabrig 2026-09-23: the raw folder was T:\Data\lab\service\..., and T: is
# \\128.120.208.24\proteomics -- the Flinders share HIVE mounts at
# /nfs/lssc0/flinders/proteomics. Nothing mapped one to the other, so the agent searched
# /quobyte, decided the files were not on HIVE, and uploaded 7.1 GB that already was
# (~20 min). Which share a path lives on is knowable locally; look it up before uploading.
#
# Prints one JSON object:
#   local, unc_or_mount_source, server, share, rest    what the path resolved to
#   candidates  HIVE paths that share should appear at, most likely first
#   hive_path   the candidate HIVE confirmed; null until verified
#   verified    true only when HIVE shows the same thing: the same byte size for a file,
#               every top-level name (the first 50, dotfiles ignored) for a directory.
#               Only a verified path may be used in place.
#   how         how it got there, in one line
# Exit: 0 verified | 1 on a share but nothing verified | 3 not on a network share (no
# ssh was made) | 2 usage.
#
# Plain bash 3.2 + coreutils, no python: gabrig's laptop had no working python3. Tests
# replace the ssh step with HIVE_EXEC=<command> and shim uname, net, mount,
# powershell.exe and cygpath on PATH.
# =============================================================================
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
IN="${1:-}"
[ -n "$IN" ] || { echo "usage: hive_path.sh <local path>" >&2; exit 2; }
NAMES_MAX=50
case "$(uname -s 2>/dev/null)" in MINGW*|MSYS*|CYGWIN*) WIN=true ;; *) WIN=false ;; esac

js() {  # a JSON string literal
  local s="$1"
  s="${s//\\/\\\\}"; s="${s//\"/\\\"}"; s="${s//$'\t'/\\t}"; s="${s//$'\r'/\\r}"; s="${s//$'\n'/\\n}"
  printf '"%s"' "$s"
}
lower() { printf '%s' "$1" | tr '[:upper:]' '[:lower:]'; }
join() { local a="${1%/}" b="${2#/}"; if [ -n "$a" ] && [ -n "$b" ]; then echo "$a/$b"; else echo "$a$b"; fi; }

SRC="" SERVER="" SHARE="" SUB="" REST="" HOW="" LOCAL_FS="$IN"
BASES=() RESTS=()

# \\server\share\a  or  //user@server/share/a  ->  SERVER, SHARE, SUB (the part below the share)
unc() {
  local u="${1//\\//}" host
  u="${u#//}"; host="${u%%/*}"; SERVER="${host##*@}"      # a macOS smbfs source carries user@
  u="${u#"$host"}"; u="${u#/}"
  SHARE="${u%%/*}"; SUB="${u#"$SHARE"}"; SUB="${SUB#/}"
}

# T:\a\b | T:/a/b | T:  ->  DRIVE=T  WREST=a/b
winpath() { DRIVE="${1%%:*}"; WREST="${1#?:}"; WREST="${WREST//\\//}"; WREST="${WREST#/}"; }

# A drive letter's UNC root (\\128.120.208.24\proteomics), or nothing for a local disk.
# `net use` first -- instant, and reachable from WSL as net.exe. PowerShell only when there
# is no net at all: it takes seconds to start, and this runs before every --put.
drive_unc() {
  local d n out
  d="$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]'):"
  for n in net net.exe; do
    command -v "$n" >/dev/null 2>&1 || continue
    out="$("$n" use 2>/dev/null)" || continue          # Samba's `net` has no `use`
    # a row is "<status> T: \\server\share <network>"; the status can be blank
    printf '%s\n' "$out" | tr -d '\r' | awk -v d="$d" '
      { for (i = 1; i < NF; i++) if (toupper($i) == d && substr($(i + 1), 1, 2) == "\\\\") {
          r = substr($0, index($0, "\\\\")); sub(/  +.*/, "", r); sub(/ +$/, "", r); print r; exit } }'
    return 0
  done
  for n in powershell.exe powershell pwsh; do
    command -v "$n" >/dev/null 2>&1 || continue
    "$n" -NoProfile -NonInteractive -Command \
      "(Get-PSDrive -Name ${d%:} -ErrorAction SilentlyContinue).DisplayRoot" 2>/dev/null | tr -d '\r'
    return 0
  done
}

# The path this shell reads for a drive or UNC path (Git Bash: /t/..., WSL: /mnt/t/...).
to_local() {
  if command -v cygpath >/dev/null 2>&1; then cygpath -u "$1"
  elif command -v wslpath >/dev/null 2>&1; then wslpath -u "$1" 2>/dev/null || echo "$1"
  else echo "$1"; fi
}

abspath() {
  local p="${1%/}"; [ -n "$p" ] || p=/
  case "$p" in /*) ;; *) p="$PWD/$p" ;; esac
  if [ -d "$p" ]; then (cd -P "$p" 2>/dev/null && pwd) || echo "$p"
  elif [ -d "$(dirname "$p")" ]; then echo "$(cd -P "$(dirname "$p")" && pwd)/$(basename "$p")"
  else echo "$p"; fi
}

from_drive() {  # $1 drive letter, $2 path below it
  local u; u="$(drive_unc "$1")"
  [ -n "$u" ] || { HOW="$1: is not a mapped network drive"; return 1; }
  SRC="$u"; unc "$u"; REST="$(join "$SUB" "$2")"; HOW="$1: is $u"
}

# The longest network mount holding $1 (smbfs //u@srv/share, cifs //srv/share, nfs
# srv:/export, WSL drvfs T: or \\srv\share). Sets SRC and UNDER (the path below it).
mount_lookup() {
  local line src rest mp best=""
  while IFS= read -r line; do
    src="${line%% on *}"; [ "$src" != "$line" ] || continue
    rest="${line#* on }"
    case "$rest" in *" type "*) mp="${rest%% type *}" ;; *) mp="${rest%% (*}" ;; esac
    mp="$(printf '%b' "$mp")"                         # Linux writes a space as \040
    case "$src" in //*|\\\\*|[A-Za-z]:|[A-Za-z]:[\\/]*|?*:/*) ;; *) continue ;; esac
    case "$1/" in "${mp%/}/"*) [ ${#mp} -gt ${#best} ] && { best="$mp"; SRC="$src"; } ;; esac
  done <<EOF
$(mount 2>/dev/null)
EOF
  [ -n "$best" ] || return 1
  UNDER="${1#"${best%/}"}"; UNDER="${UNDER#/}"
}

# ---- 1. what share is it on? (local only -- no ssh for an ordinary local path) ----
DRIVE="" WREST="" UNDER="" ABS=""
case "$IN" in
  \\\\*)                     IS=unc ;;
  //*)                       if $WIN; then IS=unc; else IS=posix; fi ;;
  [A-Za-z]:|[A-Za-z]:[\\/]*) IS=drive ;;
  *)                         IS=posix ;;
esac
if [ "$IS" = posix ]; then
  ABS="$(abspath "$IN")"; LOCAL_FS="$ABS"
  if $WIN; then        # Git Bash spells T:\x as /t/x (Cygwin: /cygdrive/t/x), and a UNC cwd as //srv/share
    case "$ABS" in
      //*) IS=unc; IN_UNC="$ABS" ;;
      /cygdrive/[A-Za-z]|/cygdrive/[A-Za-z]/*) IS=drive; winpath "${ABS:10:1}:${ABS:11}" ;;
      /[A-Za-z]|/[A-Za-z]/*)                   IS=drive; winpath "${ABS:1:1}:${ABS:2}" ;;
    esac
  fi
else
  LOCAL_FS="$(to_local "$IN")"
fi

ON_SHARE=false
case "$IS" in
  unc)   SRC="${IN_UNC:-$IN}"; unc "$SRC"; REST="$SUB"; HOW="UNC path"; ON_SHARE=true ;;
  drive) [ -n "$DRIVE" ] || winpath "$IN"; from_drive "$DRIVE" "$WREST" && ON_SHARE=true ;;
  posix)
    # Paths HIVE itself uses (a Linux box mounting quobyte/Flinders the same way) map to themselves.
    case "$ABS" in /quobyte/*|/nfs/lssc0/*) BASES+=("$ABS"); RESTS+=(""); HOW="same path as on HIVE"; ON_SHARE=true ;; esac
    if mount_lookup "$ABS"; then
      case "$SRC" in
        //*)   unc "$(printf '%b' "${SRC//\%/\\x}")"         # smbfs URL-encodes: My%20Share
               REST="$(join "$SUB" "$UNDER")"; HOW="mounted from $SRC"; ON_SHARE=true ;;
        \\\\*) unc "$SRC"; REST="$(join "$SUB" "$UNDER")"; HOW="mounted from $SRC"; ON_SHARE=true ;;
        [A-Za-z]:*) winpath "$SRC"; from_drive "$DRIVE" "$(join "$WREST" "$UNDER")" && ON_SHARE=true ;;
        *)     SERVER="${SRC%%:*}"; SHARE="$(basename "${SRC#*:}")"; REST="$UNDER"   # nfs srv:/export
               [ "$SHARE" != / ] || SHARE=""
               HOW="NFS mount of $SRC"; ON_SHARE=true ;;
      esac
    fi ;;
esac

# ---- 2. where HIVE mounts that share ----
add() {
  local i
  for ((i = 0; i < ${#BASES[@]}; i++)); do
    [ "$(join "${BASES[$i]}" "${RESTS[$i]}")" = "$(join "$1" "$2")" ] && return
  done
  BASES+=("$1"); RESTS+=("$2")
}
if [ -n "$SHARE" ]; then
  lc="$(lower "$SHARE")"
  case "$(lower "$SERVER")/$lc" in
    128.120.208.24/proteomics) add /nfs/lssc0/flinders/proteomics "$REST" ;;   # the Flinders proteomics share
    */proteomics-grp)          add /quobyte/proteomics-grp "$REST" ;;          # the Core's quobyte group dir
    *)                         add "/nfs/lssc0/flinders/$lc" "$REST"; add "/quobyte/$lc" "$REST" ;;
  esac
fi

# ---- 3. verify: the same file size, or the same top-level names, on HIVE ----
KIND="" SIZE=0 NAMES=() TOTAL=0
if [ -f "$LOCAL_FS" ]; then
  KIND=f
  SIZE="$(stat -c %s "$LOCAL_FS" 2>/dev/null || stat -f %z "$LOCAL_FS" 2>/dev/null || wc -c < "$LOCAL_FS")"
  SIZE="${SIZE//[!0-9]/}"
elif [ -d "$LOCAL_FS" ]; then
  KIND=d
  listing="$(ls "$LOCAL_FS" 2>/dev/null | LC_ALL=C sort)"          # ls without -A: no dotfiles
  while IFS= read -r n; do
    [ -n "$n" ] || continue
    TOTAL=$((TOTAL + 1)); [ ${#NAMES[@]} -lt $NAMES_MAX ] && NAMES+=("$n")
  done <<EOF
$listing
EOF
fi

HIVE_PATH="" VERIFIED=false
if [ ${#BASES[@]} -eq 0 ]; then
  :
elif [ -z "$KIND" ]; then
  HOW="$HOW; cannot read $LOCAL_FS here, so there is nothing to compare -- not verified"
else
  # Runs on HIVE via `bash -l -c`. One ssh call for every candidate: HIVE throttles rapid new
  # connections (MaxStartups; 8 quick calls locked Brett out for 20 min on 2026-09-16).
  read -r -d '' CALL <<'EOF'
v() {  # v <base> <rest> <f|d> <bytes> [names...]
  cur=$1; rest=$2; kind=$3; size=$4; shift 4
  [ -d "$cur" ] || { printf 'HIVEPATH_NO\t%s\tnot on HIVE\n' "$cur"; return; }
  IFS=/ read -r -a parts <<< "$rest"
  for c in ${parts[@]+"${parts[@]}"}; do
    [ -n "$c" ] || continue
    if [ ! -e "$cur/$c" ]; then
      # Windows and macOS SMB ignore case, so the user may have typed t:\data; HIVE's NFS
      # does not. Take a UNIQUE case-insensitive match, one level at a time.
      m=$(ls -A "$cur" 2>/dev/null | grep -ixF -- "$c")
      [ -n "$m" ] && [ "$(printf '%s\n' "$m" | wc -l)" -eq 1 ] \
        || { printf 'HIVEPATH_NO\t%s\tnot on HIVE\n' "$cur/$c"; return; }
      c=$m
    fi
    cur=$cur/$c
  done
  if [ "$kind" = f ]; then
    s=$(stat -c %s "$cur" 2>/dev/null || stat -f %z "$cur" 2>/dev/null)
    if [ -f "$cur" ] && [ "$s" = "$size" ]; then printf 'HIVEPATH_OK\t%s\tsame size, %s bytes\n' "$cur" "$s"
    else printf 'HIVEPATH_NO\t%s\t%s bytes on HIVE, %s here\n' "$cur" "${s:-no file}" "$size"; fi
  else
    [ -d "$cur" ] || { printf 'HIVEPATH_NO\t%s\tnot a directory on HIVE\n' "$cur"; return; }
    miss=0; for n in "$@"; do [ -e "$cur/$n" ] || miss=$((miss + 1)); done
    if [ $# -gt 0 ] && [ $miss -eq 0 ]; then printf 'HIVEPATH_OK\t%s\tall %s names checked are there\n' "$cur" "$#"
    else printf 'HIVEPATH_NO\t%s\t%s of %s names missing\n' "$cur" "$miss" "$#"; fi
  fi
}
EOF
  args=""
  for n in ${NAMES[@]+"${NAMES[@]}"}; do args="$args $(printf '%q' "$n")"; done
  for ((i = 0; i < ${#BASES[@]}; i++)); do
    CALL="$CALL"$'\n'"v $(printf '%q' "${BASES[$i]}") $(printf '%q' "${RESTS[$i]}") $KIND $SIZE$args"
  done
  if [ -n "${HIVE_EXEC:-}" ]; then out="$($HIVE_EXEC "$CALL" 2>&1)"
  else out="$(bash "$HERE/hive_exec.sh" "$CALL" 2>&1)"; fi
  why="" misses=""
  while IFS=$'\t' read -r tag p w; do
    case "$tag" in
      HIVEPATH_OK) [ -n "$HIVE_PATH" ] || { HIVE_PATH="$p"; why="$w"; } ;;
      HIVEPATH_NO) misses="${misses:+$misses; }$p: $w" ;;
    esac
  done <<EOF
$out
EOF
  if [ -n "$HIVE_PATH" ]; then
    VERIFIED=true; HOW="$HOW; on HIVE at $HIVE_PATH ($why)"
    [ "$KIND" = d ] && [ $TOTAL -gt ${#NAMES[@]} ] && HOW="$HOW -- first ${#NAMES[@]} of $TOTAL"
  elif [ -n "$misses" ]; then HOW="$HOW; not verified -- $misses"
  else HOW="$HOW; could not check on HIVE: $(printf '%s\n' "$out" | grep -v '^[[:space:]]*$' | head -1)"; fi
fi
$ON_SHARE || HOW="${HOW:+$HOW; }not on a network share -- upload is the only way onto HIVE"

cands=""
for ((i = 0; i < ${#BASES[@]}; i++)); do cands="$cands${cands:+, }$(js "$(join "${BASES[$i]}" "${RESTS[$i]}")")"; done
cat <<JSON
{
  "local": $(js "$IN"),
  "unc_or_mount_source": $(js "$SRC"),
  "server": $(js "$SERVER"),
  "share": $(js "$SHARE"),
  "rest": $(js "$REST"),
  "hive_path": $([ -n "$HIVE_PATH" ] && js "$HIVE_PATH" || echo null),
  "candidates": [$cands],
  "verified": $VERIFIED,
  "how": $(js "$HOW")
}
JSON
if $VERIFIED; then exit 0; elif [ ${#BASES[@]} -gt 0 ] || $ON_SHARE; then exit 1; else exit 3; fi
