# =============================================================================
# NO SHEBANG ON PURPOSE: this file is `. sourced` by acquire_tools.sh and
# build_diann_docker.sh, never executed. skill-checks.yml requires every
# shebang-bearing script to be 100755, and marking a sourced library executable
# would misrepresent how it is used.
#
# diann_release.sh -- resolve a DIA-NN "Academia" release asset URL BY VERSION.
#
# Sourced by acquire_tools.sh and build_diann_docker.sh. It is its own file so
# the rule has exactly one definition, because DIA-NN's release layout is a trap:
#
#   EVERY 2.x build is attached as an asset to the SINGLE release tagged "2.0".
#   There is no 2.1 / 2.3 / 2.6 tag. So both obvious approaches fail:
#
#     releases/tags/2.6.1 -> 404. A pinned 2.x version can NEVER download.
#     releases/latest     -> the "2.0" release, whose first Academia-Linux asset
#                            is DIA-NN-2.0 itself. This SILENTLY installs 2.0
#                            and reports it as "latest" -- the dangerous one,
#                            because it succeeds and nothing looks wrong.
#
# So we resolve by ASSET FILENAME across all releases, never by tag.
#
# Verified against the GitHub API on 2026-08-14: assets on tag "2.0" are
# 2.0, 2.0.1, 2.0.2, 2.1.0, 2.2.0, 2.3.0 (-Preview only), 2.3.1, 2.3.2, 2.5.0,
# 2.5.1, 2.6.0, 2.6.1. Releases 1.9.x and older carry their own tags AND a
# different filename scheme (diann-1.9.2.Linux.zip), which this does not cover
# -- the skill pins 2.x.
# =============================================================================

# gh_fetch <api url>: the body of a GitHub API answer on stdout, return 0. On any failure it
# prints nothing, returns 1, and writes the REASON to $GH_ERR_FILE when the caller set one --
# "GitHub could not be reached (curl exit 6)" or "GitHub answered HTTP 403 (a rate limit:
# 60 requests an hour per address without a token; set GITHUB_TOKEN)" -- because the callers
# used to send errors to /dev/null and then say "no build named X exists", which sends the
# reader looking for a wrong version when the network is the problem. GITHUB_TOKEN, when set,
# is sent to api.github.com only.
gh_fetch() {
  local url="$1" tmp code rc=0
  tmp="$(mktemp)" || return 1
  if [ -n "${GITHUB_TOKEN:-}" ] && [ "${url#https://api.github.com/}" != "$url" ]; then
    code="$(curl -sSL -H "Authorization: Bearer $GITHUB_TOKEN" -w '%{http_code}' -o "$tmp" "$url" 2>/dev/null)" || rc=$?
  else
    code="$(curl -sSL -w '%{http_code}' -o "$tmp" "$url" 2>/dev/null)" || rc=$?
  fi
  if [ "$rc" != 0 ]; then
    [ -n "${GH_ERR_FILE:-}" ] && printf 'GitHub could not be reached (curl exit %s) -- a network problem, not a missing version' "$rc" > "$GH_ERR_FILE"
    rm -f "$tmp"; return 1
  fi
  case "$code" in
    2??) cat "$tmp"; rm -f "$tmp"; return 0 ;;
    403|429) [ -n "${GH_ERR_FILE:-}" ] && printf 'GitHub answered HTTP %s -- most likely its rate limit (60 requests an hour per address without a token; set GITHUB_TOKEN, or wait) -- not a missing version' "$code" > "$GH_ERR_FILE" ;;
    404) [ -n "${GH_ERR_FILE:-}" ] && printf 'GitHub answered HTTP 404 -- that release or tag does not exist' > "$GH_ERR_FILE" ;;
    *)   [ -n "${GH_ERR_FILE:-}" ] && printf 'GitHub answered HTTP %s' "$code" > "$GH_ERR_FILE" ;;
  esac
  rm -f "$tmp"; return 1
}

# All Academia asset URLs across every release, one per line.
diann_release_assets() {
  local body
  body="$(gh_fetch "https://api.github.com/repos/vdemichev/DiaNN/releases?per_page=100")" || return 0
  printf '%s\n' "$body" \
    | grep '"browser_download_url"' \
    | sed -E 's/.*"(https[^"]+)".*/\1/'
}

# diann_asset_url <version|latest> [Linux]
# Prints the URL and returns 0, or prints nothing and returns 1.
diann_asset_url() {
  local want="${1:-latest}" kind="${2:-Linux}" all esc
  all="$(diann_release_assets)"
  [ -z "$all" ] && return 1

  if [ "$want" = "latest" ]; then
    # Highest version by FILENAME, not whatever GitHub flags "latest". Preview
    # builds are excluded here; they stay reachable by pinning them explicitly.
    printf '%s\n' "$all" \
      | sed -nE "s#.*/DIA-NN-([0-9][0-9.]*)-Academia-${kind}\.zip#\1 &#p" \
      | sort -V -k1,1 | tail -n1 | cut -d' ' -f2- | grep . && return 0
    return 1
  fi

  esc="$(printf '%s' "$want" | sed 's/\./\\./g')"
  # The literal "-Academia" immediately after the version anchors the match, so
  # a pin of 2.6.1 can never match a future 2.6.10.
  printf '%s\n' "$all" | grep -E "/DIA-NN-${esc}-Academia-${kind}\.zip$" | head -n1 | grep . && return 0
  # Same version, suffixed build (2.3.0 ships ONLY as -Academia-Linux-Preview.zip).
  printf '%s\n' "$all" | grep -E "/DIA-NN-${esc}-Academia-${kind}[^/]*\.zip$" | head -n1 | grep . && return 0
  return 1
}
