#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# DE-LIMP download / usage stats — one-shot summary across all four surfaces.
#
# Sources:
#   1. GitHub repo traffic (clones / views, last 14 days, owner-only)
#   2. GitHub Release asset downloads (all-time, public)
#   3. Hugging Face Space (likes — visible to anyone)
#   4. Docker Hub base image pulls (all-time, public)
#
# Usage:
#   bash scripts/download_stats.sh           # plain output
#   bash scripts/download_stats.sh --json    # machine-readable JSON
#
# Requirements:
#   - `gh` (GitHub CLI), authenticated as repo owner for traffic API
#   - `curl`, `jq`
#
# Notes:
#   - Clones != users (CI runners and mirrors count). Trust `uniques`.
#   - HF "downloads" don't really exist for Spaces — likes is the closest
#     public number. Owner dashboard at huggingface.co/spaces/.../analytics
#     has unique-user counts.
#   - Docker Hub `pull_count` is all-time and public.
# ----------------------------------------------------------------------------
set -euo pipefail

REPO="bsphinney/DE-LIMP"
HF_SPACE="brettsp/de-limp-proteomics"
DOCKER_IMAGE="brettphinney/delimp-base"

OUTPUT_JSON=0
[ "${1:-}" = "--json" ] && OUTPUT_JSON=1

need() {
  command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found on PATH" >&2; exit 1; }
}
need curl
need jq

# ---- 1. GitHub traffic (owner-only) ---------------------------------------
gh_clones_count=null
gh_clones_unique=null
gh_views_count=null
gh_views_unique=null
gh_stars=null
gh_forks=null
gh_open_issues=null

if command -v gh >/dev/null 2>&1; then
  if gh auth status >/dev/null 2>&1; then
    if clones_json=$(gh api "repos/${REPO}/traffic/clones" 2>/dev/null); then
      gh_clones_count=$(echo "$clones_json" | jq -r '.count')
      gh_clones_unique=$(echo "$clones_json" | jq -r '.uniques')
    fi
    if views_json=$(gh api "repos/${REPO}/traffic/views" 2>/dev/null); then
      gh_views_count=$(echo "$views_json" | jq -r '.count')
      gh_views_unique=$(echo "$views_json" | jq -r '.uniques')
    fi
  fi
fi

repo_json=$(curl -s "https://api.github.com/repos/${REPO}" || true)
if [ -n "$repo_json" ]; then
  gh_stars=$(echo "$repo_json"      | jq -r '.stargazers_count // null')
  gh_forks=$(echo "$repo_json"      | jq -r '.forks_count // null')
  gh_open_issues=$(echo "$repo_json"| jq -r '.open_issues_count // null')
fi

# ---- 2. GitHub Release asset downloads ------------------------------------
releases_json='[]'
release_total=0
if command -v gh >/dev/null 2>&1 && gh auth status >/dev/null 2>&1; then
  releases_json=$(gh release list --repo "$REPO" --limit 100 --json tagName,assets 2>/dev/null || echo '[]')
fi
release_total=$(echo "$releases_json" | jq -r '[.[].assets[].downloadCount // 0] | add // 0')

# ---- 3. Hugging Face Space ------------------------------------------------
hf_likes=null
hf_last_modified=null
hf_json=$(curl -s "https://huggingface.co/api/spaces/${HF_SPACE}" || true)
if [ -n "$hf_json" ]; then
  hf_likes=$(echo "$hf_json"        | jq -r '.likes // null')
  hf_last_modified=$(echo "$hf_json"| jq -r '.lastModified // null')
fi

# ---- 4. Docker Hub --------------------------------------------------------
docker_pulls=null
docker_stars=null
docker_json=$(curl -s "https://hub.docker.com/v2/repositories/${DOCKER_IMAGE}/" || true)
if [ -n "$docker_json" ]; then
  docker_pulls=$(echo "$docker_json"| jq -r '.pull_count // null')
  docker_stars=$(echo "$docker_json"| jq -r '.star_count // null')
fi

# ---- Output ---------------------------------------------------------------
if [ "$OUTPUT_JSON" = "1" ]; then
  jq -n \
    --argjson gh_clones_count "$gh_clones_count" \
    --argjson gh_clones_unique "$gh_clones_unique" \
    --argjson gh_views_count "$gh_views_count" \
    --argjson gh_views_unique "$gh_views_unique" \
    --argjson gh_stars "$gh_stars" \
    --argjson gh_forks "$gh_forks" \
    --argjson gh_open_issues "$gh_open_issues" \
    --argjson release_total "$release_total" \
    --argjson hf_likes "$hf_likes" \
    --argjson docker_pulls "$docker_pulls" \
    --argjson docker_stars "$docker_stars" \
    --arg     hf_last_modified "$hf_last_modified" \
    --argjson releases "$releases_json" \
    '{
      generated_at: (now | strftime("%Y-%m-%d %H:%M:%S")),
      github: {
        repo: $ENV.REPO,
        stars: $gh_stars,
        forks: $gh_forks,
        open_issues: $gh_open_issues,
        traffic_14d: {
          clones: $gh_clones_count, clones_unique: $gh_clones_unique,
          views:  $gh_views_count,  views_unique:  $gh_views_unique
        },
        release_downloads_all_time: $release_total,
        per_release: $releases
      },
      huggingface: { space: $ENV.HF_SPACE, likes: $hf_likes, last_modified: $hf_last_modified },
      docker:      { image: $ENV.DOCKER_IMAGE, pull_count: $docker_pulls, stars: $docker_stars }
    }' REPO="$REPO" HF_SPACE="$HF_SPACE" DOCKER_IMAGE="$DOCKER_IMAGE"
  exit 0
fi

date_str=$(date '+%Y-%m-%d %H:%M:%S')
printf '\n=== DE-LIMP usage stats — %s ===\n\n' "$date_str"

printf 'GitHub repo (%s)\n' "$REPO"
printf '  stars            : %s\n'  "$gh_stars"
printf '  forks            : %s\n'  "$gh_forks"
printf '  open issues      : %s\n'  "$gh_open_issues"
printf '  clones (14d)     : %s total, %s unique\n' "$gh_clones_count" "$gh_clones_unique"
printf '  views  (14d)     : %s total, %s unique\n' "$gh_views_count"  "$gh_views_unique"
printf '  release downloads: %s (all-time, sum of all assets across all releases)\n' "$release_total"
echo

if [ "$releases_json" != "[]" ]; then
  printf 'Top release assets (download count, all-time):\n'
  echo "$releases_json" | jq -r '
    [.[] | .tag = .tagName | .assets[] | {tag: .tag, name: .name, dl: .downloadCount}]
    | sort_by(-.dl) | .[0:10]
    | .[] | "  \(.dl|tostring|.[0:7]|gsub(" ";""))  \(.tag)/\(.name)"
  '
  echo
fi

printf 'Hugging Face Space (%s)\n' "$HF_SPACE"
printf '  likes            : %s\n'  "$hf_likes"
printf '  last modified    : %s\n'  "$hf_last_modified"
echo

printf 'Docker Hub (%s)\n' "$DOCKER_IMAGE"
printf '  pulls (all-time) : %s\n'  "$docker_pulls"
printf '  stars            : %s\n'  "$docker_stars"
echo

cat <<EOF
Caveats:
  - GitHub clones include CI runners and mirrors. 'unique' is a better signal.
  - HF Space likes != users; deeper unique-user counts at:
    https://huggingface.co/spaces/${HF_SPACE}/settings/analytics
  - Docker pulls include automated builds and CI; ratio of pulls/stars is a
    rough sanity check.
  - Stars + forks + Discussions activity is the most honest 'active users'
    proxy. Open the About -> Community tab in the running app for the
    rolling daily snapshots produced by .github/workflows/track-stats.yml.
EOF
