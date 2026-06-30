#!/usr/bin/env bash
# =============================================================================
# acquire_tools.sh  --  Resolve (and download if needed) the search engines for
# the detected environment, then write a tool manifest (tools.json) describing
# how to invoke each one. Idempotent: existing installs are reused.
#
#   Usage: bash acquire_tools.sh <platform_class> [tools_root]
#          platform_class = hpc | mac | linux   (from detect_env.sh)
#          tools_root      = where to install    (default ~/.proteomics-pipeline/tools)
#
#   Version pinning (PLAN.md §7c): the orchestrator passes the workflow bundle's
#   engine + version so we install/resolve THAT EXACT version, not GitHub
#   "latest". Pass via env:
#          PIN_ENGINE=diann PIN_VERSION=2.6.0 bash acquire_tools.sh hpc
#   Installs are cached under <root>/<engine>/<version>/ so multiple pinned
#   versions coexist and results stay reproducible.
#
# Engine acquisition matrix (current upstream state, June 2026):
#   Sage     open source, cross-platform binary -> always downloadable.
#   DIA-NN   Win/Linux only; free "Academia" zip on GitHub. No macOS build ->
#            on mac, run via Docker. On HIVE, reuse the existing .sif.
#   FragPipe Java app; MSFragger + IonQuant are LICENSE-GATED and cannot be
#            silently downloaded.
# =============================================================================
set -uo pipefail

CLASS="${1:?platform_class required (hpc|mac|linux)}"
ROOT="${2:-$HOME/.proteomics-pipeline/tools}"
PIN_ENGINE="${PIN_ENGINE:-}"      # diann | sage | fragpipe  (optional)
PIN_VERSION="${PIN_VERSION:-}"    # exact version to honor   (optional)
mkdir -p "$ROOT"
MANIFEST="$ROOT/tools.json"
OS="$(uname -s | tr '[:upper:]' '[:lower:]')"; ARCH="$(uname -m)"
have() { command -v "$1" >/dev/null 2>&1; }

# normalize a version into the tags engines actually use
sage_tag()  { local v="$1"; [ -z "$v" ] && { echo ""; return; }; case "$v" in v*) echo "$v";; *) echo "v$v";; esac; }

# Asset URL from a SPECIFIC release tag (pinned); falls back to latest if tag is blank.
asset_url_tag() { # owner/repo  tag  pattern
  local repo="$1" tag="$2" pat="$3" api
  if [ -n "$tag" ]; then api="https://api.github.com/repos/$repo/releases/tags/$tag"
  else                  api="https://api.github.com/repos/$repo/releases/latest"; fi
  curl -fsSL "$api" 2>/dev/null \
    | grep '"browser_download_url"' | sed -E 's/.*"(https[^"]+)".*/\1/' \
    | grep -iE "$pat" | head -n1
}

DIANN_CMD="null"; SAGE_CMD="null"; FRAGPIPE_CMD="null"; ALPHADIA_CMD="null"; NOTES=()
SAGE_VER="${PIN_VERSION:-latest}"; DIANN_VER="${PIN_VERSION:-latest}"

# Only honor PIN_VERSION for the engine it names.
pin_for() { [ -n "$PIN_ENGINE" ] && [ "$PIN_ENGINE" = "$1" ] && [ -n "$PIN_VERSION" ]; }

# ---------------------------------------------------------------- Sage --------
acquire_sage() {
  # Prefer a Sage already on PATH (setup.sh installs sage-proteomics into the
  # conda env and activate.sh puts it on PATH). No need to re-download.
  if have sage; then SAGE_CMD="$(command -v sage)"; SAGE_VER="env"
    NOTES+=("Sage: using the conda-env build on PATH ($SAGE_CMD)."); return; fi
  local ver tag; ver="latest"; tag=""
  if pin_for sage; then ver="$PIN_VERSION"; tag="$(sage_tag "$PIN_VERSION")"; fi
  SAGE_VER="$ver"
  local dir="$ROOT/sage/$ver" bin
  bin="$(find "$dir" -name sage -type f -perm -u+x 2>/dev/null | head -n1)"
  if [ -n "$bin" ]; then SAGE_CMD="$bin"; return; fi
  mkdir -p "$dir"
  local pat
  case "$OS-$ARCH" in
    linux-x86_64)  pat="x86_64-unknown-linux-gnu.*tar" ;;
    darwin-arm64)  pat="aarch64-apple-darwin.*tar" ;;
    darwin-x86_64) pat="x86_64-apple-darwin.*tar" ;;
    *)             pat="x86_64-unknown-linux-gnu.*tar" ;;
  esac
  local url; url="$(asset_url_tag lazear/sage "$tag" "$pat")"
  if [ -z "$url" ]; then NOTES+=("Sage: could not resolve a release asset for $OS-$ARCH version='$ver'; see https://github.com/lazear/sage/releases"); return; fi
  echo "  Sage: downloading $url"
  curl -fsSL "$url" -o "$dir/sage.tar.gz" && { tar -xzf "$dir/sage.tar.gz" -C "$dir" --strip-components=1 2>/dev/null || tar -xzf "$dir/sage.tar.gz" -C "$dir"; }
  bin="$(find "$dir" -name sage -type f 2>/dev/null | head -n1)"
  [ -n "$bin" ] && chmod +x "$bin" && SAGE_CMD="$bin"
  NOTES+=("Sage $ver: run with --parquet; v0.14.x needs sage_protein_groups.py post-hoc for protein rollup. Telemetry off via --disable-telemetry-i-dont-want-to-improve-sage.")
}

# --------------------------------------------------------------- DIA-NN -------
acquire_diann() {
  local ver; ver="latest"; pin_for diann && ver="$PIN_VERSION"
  DIANN_VER="$ver"
  case "$CLASS" in
    hpc)
      # On HIVE the Proteomics Core keeps DIA-NN as NATIVE builds, e.g. 2.6.0 at
      #   /quobyte/proteomics-grp/dia-nn/build_260/diann-2.6.0/diann-linux
      # (older 2.3.0 is a .sif). Match the PINNED version exactly — never silently
      # substitute a different version (reproducibility).
      local DN=/quobyte/proteomics-grp/dia-nn bin="" sif=""
      if [ "$ver" != "latest" ]; then
        bin="$(ls -1 "$DN"/*/diann-"$ver"/diann-linux 2>/dev/null | head -n1)"
        [ -z "$bin" ] && bin="$(find "$DN" -maxdepth 3 -path "*diann-$ver*/diann-linux" -type f 2>/dev/null | head -n1)"
        [ -z "$bin" ] && sif="$(ls -1 "$DN"/*"$ver"*.sif 2>/dev/null | sort -V | tail -n1)"
      else
        bin="$(find "$DN" -maxdepth 3 -name diann-linux -type f 2>/dev/null | sort -V | tail -n1)"
        [ -z "$bin" ] && sif="$(ls -1 "$DN"/*.sif 2>/dev/null | sort -V | tail -n1)"
      fi
      if [ -n "$bin" ]; then
        DIANN_CMD="$bin"
        NOTES+=("DIA-NN $ver: using the HIVE Proteomics Core native build $bin. Reference invocation: $DN/run_diann_*.sbatch.")
        return
      elif [ -n "$sif" ] && have apptainer && { [ "$ver" = "latest" ] || printf '%s' "$sif" | grep -q "$ver"; }; then
        DIANN_CMD="apptainer exec --bind /quobyte:/quobyte $sif /diann-*/diann-linux"
        NOTES+=("DIA-NN $ver: reusing HIVE container $sif.")
        return
      else
        NOTES+=("DIA-NN $ver NOT found on HIVE under $DN — NOT substituting a different version. Present: $(ls -d "$DN"/*.sif "$DN"/build_*/diann-* 2>/dev/null | xargs -n1 basename 2>/dev/null | tr '\n' ' '). Will fetch the pinned Academia build into your space instead (or build $ver).")
      fi ;;
    mac)
      if [ -n "${DIANN_DOCKER_IMAGE:-}" ] && have docker; then
        DIANN_CMD="docker run --rm -v \$PWD:/data ${DIANN_DOCKER_IMAGE} diann-linux"
        NOTES+=("DIA-NN: using Docker image \$DIANN_DOCKER_IMAGE on macOS (no native mac build exists).")
        return
      fi
      NOTES+=("DIA-NN on macOS: no native build. Set DIANN_DOCKER_IMAGE to a built image (version $ver), or build from the Academia Linux zip's Dockerfile. See references/environment.md.")
      ;;
  esac
  # linux native (or fallback): fetch the pinned Academia Linux zip
  local dir="$ROOT/diann/$ver"; mkdir -p "$dir"
  local existing; existing="$(find "$dir" -name 'diann-linux' -type f 2>/dev/null | head -n1)"
  if [ -n "$existing" ]; then DIANN_CMD="$existing"; return; fi
  # DIA-NN release tags are usually the bare version (e.g. "2.6.0"); blank => latest.
  local tag=""; [ "$ver" != "latest" ] && tag="$ver"
  local url; url="$(asset_url_tag vdemichev/DiaNN "$tag" 'Academia-Linux.*zip')"
  if [ -z "$url" ]; then NOTES+=("DIA-NN: could not resolve Academia Linux asset for version '$ver'; download manually from https://github.com/vdemichev/DiaNN/releases (free for academic use)."); return; fi
  echo "  DIA-NN: downloading $url"
  curl -fsSL "$url" -o "$dir/diann.zip" && (cd "$dir" && unzip -oq diann.zip)
  existing="$(find "$dir" -name 'diann-linux' -type f 2>/dev/null | head -n1)"
  if [ -n "$existing" ]; then chmod +x "$existing"; DIANN_CMD="$existing"
    NOTES+=("DIA-NN $ver: Linux native needs glibc>=Mint 21.2 + .NET 8. If missing, prefer Docker/Apptainer (the zip ships a Dockerfile).")
  fi
}

# ------------------------------------------------------------- FragPipe -------
acquire_fragpipe() {
  local ver tag; ver="latest"; tag=""
  if pin_for fragpipe; then ver="$PIN_VERSION"; tag="$PIN_VERSION"; fi
  local dir="$ROOT/fragpipe/$ver"
  local fp; fp="$(find "$dir" -name 'fragpipe' -type f 2>/dev/null | head -n1)"
  if [ -z "$fp" ]; then
    mkdir -p "$dir"
    local url; url="$(asset_url_tag Nesvilab/FragPipe "$tag" 'FragPipe-[0-9].*zip' | grep -v jre)"
    if [ -z "$url" ]; then NOTES+=("FragPipe: could not resolve release zip for version '$ver'; see https://github.com/Nesvilab/FragPipe/releases"); return; fi
    echo "  FragPipe: downloading $url"
    curl -fsSL "$url" -o "$dir/fragpipe.zip" && (cd "$dir" && unzip -oq fragpipe.zip)
    fp="$(find "$dir" -path '*/bin/fragpipe' -type f 2>/dev/null | head -n1)"
    [ -n "$fp" ] && chmod +x "$fp"
  fi
  if ! have java; then NOTES+=("FragPipe: requires Java 9+ on PATH. Install a JDK first."); fi
  local tools="${FRAGPIPE_TOOLS_FOLDER:-$dir/tools}"
  if [ -n "$fp" ] && ls "$tools"/MSFragger*.jar >/dev/null 2>&1; then
    FRAGPIPE_CMD="$fp"
    NOTES+=("FragPipe $ver: ready (tools folder $tools). run_search.py adds --config-tools-folder via FRAGPIPE_TOOLS_FOLDER.")
  elif [ -n "$fp" ]; then
    FRAGPIPE_CMD="$fp"
    NOTES+=("FragPipe $ver: downloaded, but MSFragger/IonQuant are license-gated and absent. Run the GUI once to accept the academic license, OR set FRAGPIPE_TOOLS_FOLDER. Headless search disabled until then.")
  fi
}

# ------------------------------------------------------------- AlphaDIA -------
# Apache-2.0 (commercial use OK) — the open-source DIA alternative to DIA-NN,
# whose free "Academia" build is academic/non-profit only. pip-installed into the
# active env. Deep-learning based: a GPU is strongly recommended (CPU is slow).
acquire_alphadia() {
  # On HIVE the Proteomics Core keeps an AlphaDIA container.
  if [ "$CLASS" = "hpc" ]; then
    local asif=/quobyte/proteomics-grp/apptainers/alphadia.sif
    if [ -f "$asif" ] && have apptainer; then
      ALPHADIA_CMD="apptainer exec --bind /quobyte:/quobyte $asif alphadia"
      NOTES+=("AlphaDIA: reusing the HIVE container $asif (Apache-2.0, commercial-OK)."); return
    fi
  fi
  if have alphadia; then ALPHADIA_CMD="$(command -v alphadia)"
    NOTES+=("AlphaDIA: using the alphadia on PATH ($ALPHADIA_CMD). Apache-2.0 (commercial use OK)."); return; fi
  local spec="alphadia"; pin_for alphadia && spec="alphadia==$PIN_VERSION"
  if have pip || have pip3; then
    echo "  AlphaDIA: pip install $spec (large download; GPU recommended)..."
    { pip install "$spec" >/dev/null 2>&1 || pip3 install "$spec" >/dev/null 2>&1; }
    if have alphadia; then ALPHADIA_CMD="$(command -v alphadia)"
      NOTES+=("AlphaDIA: installed via pip (Apache-2.0, commercial-OK). GPU strongly recommended; verify with 'alphadia --check'."); return; fi
  fi
  NOTES+=("AlphaDIA: not installed. Install with 'pip install alphadia' in the env (Apache-2.0, commercial-OK; GPU recommended). https://github.com/MannLabs/alphadia")
}

echo "[acquire] platform=$CLASS  os=$OS arch=$ARCH  root=$ROOT  pin=${PIN_ENGINE:-none}/${PIN_VERSION:-none}"
acquire_sage
acquire_diann
acquire_fragpipe
# AlphaDIA: discover the HIVE container automatically; otherwise only acquire when
# pinned/requested or already present (it's a large GPU pip package — don't auto-install).
if [ "${PIN_ENGINE:-}" = "alphadia" ] || have alphadia \
   || { [ "$CLASS" = "hpc" ] && [ -f /quobyte/proteomics-grp/apptainers/alphadia.sif ]; }; then
  acquire_alphadia
fi

# ---- write manifest ----------------------------------------------------------
esc() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }
{
  printf '{\n'
  printf '  "platform_class": "%s",\n' "$CLASS"
  printf '  "tools_root": "%s",\n' "$(esc "$ROOT")"
  printf '  "pinned": {"engine": "%s", "version": "%s"},\n' "${PIN_ENGINE:-}" "${PIN_VERSION:-}"
  printf '  "versions": {"diann": "%s", "sage": "%s"},\n' "$(esc "$DIANN_VER")" "$(esc "$SAGE_VER")"
  printf '  "diann":    %s,\n'  "$( [ "$DIANN_CMD" = null ]    && echo null || printf '"%s"' "$(esc "$DIANN_CMD")" )"
  printf '  "sage":     %s,\n'  "$( [ "$SAGE_CMD" = null ]     && echo null || printf '"%s"' "$(esc "$SAGE_CMD")" )"
  printf '  "fragpipe": %s,\n'  "$( [ "$FRAGPIPE_CMD" = null ] && echo null || printf '"%s"' "$(esc "$FRAGPIPE_CMD")" )"
  printf '  "alphadia": %s,\n'  "$( [ "$ALPHADIA_CMD" = null ] && echo null || printf '"%s"' "$(esc "$ALPHADIA_CMD")" )"
  printf '  "notes": [\n'
  for i in "${!NOTES[@]}"; do
    printf '    "%s"%s\n' "$(esc "${NOTES[$i]}")" "$( [ "$i" -lt $((${#NOTES[@]}-1)) ] && echo , )"
  done
  printf '  ]\n}\n'
} | tee "$MANIFEST"

echo "[acquire] manifest written to $MANIFEST"
