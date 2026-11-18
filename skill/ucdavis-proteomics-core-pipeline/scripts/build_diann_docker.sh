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

# Resolve the Academia-Linux zip (free academic build) for the requested version.
# Shared with acquire_tools.sh -- see diann_release.sh for why resolution is by
# asset filename and never by release tag.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/diann_release.sh"
. "$SCRIPT_DIR/verify_download.sh"
GH_ERR_FILE="$(mktemp)"; export GH_ERR_FILE
trap 'rm -f "$GH_ERR_FILE"' EXIT
URL="$(diann_asset_url "$VER" Linux)"
if [ -z "$URL" ] && [ "$VER" != "latest" ] && [ -s "$GH_ERR_FILE" ]; then
  # GitHub's API could not answer; a pinned 2.x build has a fixed address (diann_release.sh)
  echo "[diann-docker] $(cat "$GH_ERR_FILE"); using the pinned build's fixed address." >&2
  URL="https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-$VER-Academia-Linux.zip"
fi
if [ -z "$URL" ] && [ -s "$GH_ERR_FILE" ]; then
  echo "[diann-docker] Could not look up DIA-NN '$VER': $(cat "$GH_ERR_FILE")." >&2
  echo "  Nothing was downloaded. Try again, or set GITHUB_TOKEN." >&2
  exit 5
fi
if [ -z "$URL" ]; then
  echo "[diann-docker] No DIA-NN Academia-Linux build named '$VER' exists." >&2
  echo "  Check the version against https://github.com/vdemichev/DiaNN/releases (free for academics)," >&2
  echo "  or download the zip manually, unzip it into $WORK, then re-run this script." >&2
  exit 5
fi
# Tag the image with the version we actually resolved, so "latest" does not
# produce an image whose tag says nothing about what is inside it.
GOT="$(printf '%s' "$URL" | sed -nE 's#.*/DIA-NN-([0-9][0-9.]*)-Academia-.*#\1#p')"
if [ -n "$GOT" ] && [ "$GOT" != "$VER" ]; then
  echo "[diann-docker] '$VER' resolved to DIA-NN $GOT"
  VER="$GOT"; TAG="proteomics-pipeline/diann:${VER}"
  if docker image inspect "$TAG" >/dev/null 2>&1; then
    echo "$TAG" > "$PP_HOME/diann_docker_image"
    echo "[diann-docker] image already built: $TAG"; exit 0
  fi
fi

echo "[diann-docker] downloading $URL"
rm -f "$WORK/diann.zip.part"
if ! curl -fsSL "$URL" -o "$WORK/diann.zip.part"; then
  rm -f "$WORK/diann.zip.part"; echo "[diann-docker] The download failed or was cut short. Nothing was built." >&2; exit 8
fi
# against the sha256 the vendor published for this build (pinned_checksums.txt); a mismatch stops here
verify_download "diann-${VER}-linux" "$WORK/diann.zip.part" "DIA-NN $VER" ; vrc=$?
if [ "$vrc" = 1 ]; then rm -f "$WORK/diann.zip.part"; echo "[diann-docker] The zip does not match the published checksum, so it was deleted and nothing was built." >&2; exit 9; fi
mv -f "$WORK/diann.zip.part" "$WORK/diann.zip" && (cd "$WORK" && unzip -oq diann.zip) || { echo "[diann-docker] The zip could not be unpacked. Nothing was built." >&2; exit 8; }

# DIA-NN's zip ships a Dockerfile (and make_docker.sh). Build from it directly.
DF="$(find "$WORK" -iname Dockerfile -type f | head -n1)"
if [ -z "$DF" ]; then
  echo "[diann-docker] No Dockerfile found in the DIA-NN zip. Check $WORK; the layout may have changed." >&2
  exit 6
fi
CTX="$(dirname "$DF")"
echo "[diann-docker] building $TAG from $DF (this takes a few minutes the first time)..."
# --platform linux/amd64: DIA-NN has no arm64 build; on Apple Silicon the image is built (and run,
# see diann_docker.sh) as amd64 under emulation. Without it the image is built for the host.
if docker build --platform linux/amd64 -t "$TAG" -f "$DF" "$CTX"; then
  echo "$TAG" > "$PP_HOME/diann_docker_image"
  echo "[diann-docker] built $TAG  ->  recorded in $PP_HOME/diann_docker_image"
  echo "[diann-docker] re-source activate.sh so DIANN_DOCKER_IMAGE is exported."
else
  echo "[diann-docker] docker build failed. See the output above." >&2
  exit 7
fi
