#!/usr/bin/env bash
# =============================================================================
# build_diann_docker.sh  --  Build a DIA-NN Docker image from DIA-NN's OWN
# official recipe, so macOS users (no native DIA-NN build exists) can run DIA
# searches without hand-assembling anything. Uses only the free academic
# "Academia-Linux" release + the Dockerfile it ships — no third-party images.
#
# Requires Docker to be installed and running. If it isn't, this prints the one
# step the user must take (install Docker Desktop) and exits non-zero.
#
#   Usage: bash build_diann_docker.sh [version]
#          version defaults to the workflow's pinned DIA-NN version, else latest.
#
# On success: writes the image tag to ~/.proteomics-pipeline/diann_docker_image
# (activate.sh exports it as DIANN_DOCKER_IMAGE, which acquire_tools.sh uses).
# =============================================================================
set -uo pipefail

PP_HOME="${PP_HOME:-$HOME/.proteomics-pipeline}"
VER="${1:-${PIN_VERSION:-latest}}"
TAG="proteomics-pipeline/diann:${VER}"
WORK="$PP_HOME/diann_build"
mkdir -p "$WORK"

if ! command -v docker >/dev/null 2>&1; then
  cat >&2 <<'MSG'
[diann-docker] Docker is not installed.

DIA-NN has no native macOS build, so DIA searches on a Mac need Docker.
ONE step to fix it:

  1. Install Docker Desktop (free):  https://docs.docker.com/desktop/setup/install/mac-install/
  2. Open Docker Desktop once so it finishes setup (whale icon in the menu bar).
  3. Re-run the analysis — the AI will pick up from here automatically.

(DDA data via Sage, and all the differential-expression steps, do NOT need
Docker — only DIA-NN does.)
MSG
  exit 3
fi

if ! docker info >/dev/null 2>&1; then
  echo "[diann-docker] Docker is installed but not running. Open Docker Desktop, wait for the whale icon to go steady, then re-run." >&2
  exit 4
fi

# reuse an already-built image
if docker image inspect "$TAG" >/dev/null 2>&1; then
  echo "$TAG" > "$PP_HOME/diann_docker_image"
  echo "[diann-docker] image already built: $TAG"; exit 0
fi

# resolve the Academia-Linux zip (free academic build) for the requested version
asset_url() { # tag pattern
  local api
  if [ "$1" != "latest" ]; then api="https://api.github.com/repos/vdemichev/DiaNN/releases/tags/$1"
  else api="https://api.github.com/repos/vdemichev/DiaNN/releases/latest"; fi
  curl -fsSL "$api" 2>/dev/null | grep '"browser_download_url"' \
    | sed -E 's/.*"(https[^"]+)".*/\1/' | grep -iE 'Academia-Linux.*zip' | head -n1
}
URL="$(asset_url "$VER")"
if [ -z "$URL" ]; then
  echo "[diann-docker] Could not find the DIA-NN Academia-Linux zip for version '$VER'." >&2
  echo "  Download it manually (free for academics) from https://github.com/vdemichev/DiaNN/releases" >&2
  echo "  unzip it into $WORK, then re-run this script." >&2
  exit 5
fi

echo "[diann-docker] downloading $URL"
curl -fsSL "$URL" -o "$WORK/diann.zip" && (cd "$WORK" && unzip -oq diann.zip)

# DIA-NN's zip ships a Dockerfile (and make_docker.sh). Build from it directly.
DF="$(find "$WORK" -iname Dockerfile -type f | head -n1)"
if [ -z "$DF" ]; then
  echo "[diann-docker] No Dockerfile found in the DIA-NN zip. Check $WORK; the layout may have changed." >&2
  exit 6
fi
CTX="$(dirname "$DF")"
echo "[diann-docker] building $TAG from $DF (this takes a few minutes the first time)..."
if docker build -t "$TAG" -f "$DF" "$CTX"; then
  echo "$TAG" > "$PP_HOME/diann_docker_image"
  echo "[diann-docker] built $TAG  ->  recorded in $PP_HOME/diann_docker_image"
  echo "[diann-docker] re-source activate.sh so DIANN_DOCKER_IMAGE is exported."
else
  echo "[diann-docker] docker build failed. See the output above." >&2
  exit 7
fi
