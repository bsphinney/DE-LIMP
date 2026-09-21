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
#          PIN_ENGINE=diann PIN_VERSION=2.6.1 bash acquire_tools.sh hpc   # the pin resolve_defaults.py writes
#   Installs are cached under <root>/<engine>/<version>/ so multiple pinned
#   versions coexist and results stay reproducible.
#
#   tools.json `versions.<engine>` is the build the resolved command IS, never the
#   request: an unpinned ("latest") run records the number it resolved to, and a
#   build whose version cannot be determined records "" rather than echoing the
#   request back (Sage taken from a conda env on PATH records "env", a source rather
#   than a build). run_search.py copies this into search_provenance.json, and FRAN
#   stores it as the engine version -- the word "latest" there names nothing.
#   (Measured 2026-09-16: both HIVE tools.json files in use said "sage": "latest";
#   the cached tarball was sage-v0.14.7.)
#   There is a key for every engine tools.json can name a command for -- diann, sage,
#   radiant, fragpipe, alphadia. A MISSING key is not "no version this time": it makes
#   every search that engine ever runs record `version: null`, because run_search.py has
#   nothing to read. run_search.py accepts only a value SHAPED like a version, so a word
#   written here ("latest", "env", "nightly") is refused there rather than deposited.
#
# Engine acquisition matrix (current upstream state, June 2026):
#   Sage     open source, cross-platform binary -> always downloadable.
#   DIA-NN   Win/Linux only; free "Academia" zip on GitHub. No macOS build ->
#            on mac, run via Docker. On HIVE, reuse the Core's native builds
#            (or its older .sif) under $DIANN_HIVE_DIR.
#   FragPipe Java app; MSFragger + IonQuant are LICENSE-GATED and cannot be
#            silently downloaded.
# =============================================================================
set -uo pipefail

CLASS="${1:?platform_class required (hpc|mac|linux)}"
ROOT="${2:-$HOME/.proteomics-pipeline/tools}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# DIA-NN resolves by asset filename, not by release tag -- see diann_release.sh
# for why (every 2.x build hangs off the single release tagged "2.0").
. "$SCRIPT_DIR/diann_release.sh"
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
RADIANT_CMD="null"; RADIANT_RUNTIME=""; RADIANT_IMAGE=""; RADIANT_VER=""
# Empty until a command is resolved; see "versions" in the header. There is one per engine
# tools.json can name a command for: an engine with no `versions` key records
# `version: null` in search_provenance.json for every search it ever runs, because
# run_search.py has nothing to read.
SAGE_VER=""; DIANN_VER=""; FRAGPIPE_VER=""; ALPHADIA_VER=""

# The Core's shared DIA-NN builds on HIVE. Overridable so the resolution below can be
# tested against a mocked layout (tests/test_acquire_tools_versions.py), and so another
# site can point it at its own shared builds; nobody else needs to set it.
DIANN_HIVE_DIR="${DIANN_HIVE_DIR:-/quobyte/proteomics-grp/dia-nn}"
# The folders searched for the Core's Radiant .sif, colon-separated, in order. Overridable for
# the same two reasons. `-`, not `:-`: an explicitly EMPTY RADIANT_HIVE_DIRS means "this host
# has no shared folder", and used to fall back to the real /quobyte paths instead -- so a test
# (or a site) that set it empty silently searched HIVE's own folders.
RADIANT_HIVE_DIRS="${RADIANT_HIVE_DIRS-/quobyte/proteomics-grp/apptainers:/quobyte/proteomics-grp/radiant}"

# Only honor PIN_VERSION for the engine it names.
pin_for() { [ -n "$PIN_ENGINE" ] && [ "$PIN_ENGINE" = "$1" ] && [ -n "$PIN_VERSION" ]; }

# diann_version_of <path>: the DIA-NN version a build's PATH names, or nothing.
#   .../build_261/diann-2.6.1/diann-linux -> 2.6.1     .../diann_2.3.0.sif -> 2.3.0
# Read from the path, not by running the binary: that would start DIA-NN on what is
# often a login node, and the .sif needs apptainer just to print a banner.
diann_version_of() {
  printf '%s' "$1" | sed -nE 's#.*[Dd][Ii][Aa]-?[Nn][Nn][-_]([0-9]+(\.[0-9]+)+).*#\1#p'
}

# newest_diann <candidate paths on stdin>: the path with the highest VERSION.
# Never `sort -V` the whole paths: that orders by whatever folder sits in front of the
# version, so a 2.5.1 kept under z_archive/ would outrank build_261/diann-2.6.1.
# Candidates whose path names no version are dropped -- "newest" is unknowable for them.
newest_diann() {
  local p v
  while IFS= read -r p; do
    v="$(diann_version_of "$p")"
    [ -n "$v" ] && printf '%s %s\n' "$v" "$p"
  done | sort -V -k1,1 | tail -n1 | cut -d' ' -f2-
}

# fragpipe_version_of <path or asset name>: the release a FragPipe path or release asset
# names, or nothing.
#   .../fragpipe/latest/fragpipe-24.0/bin/fragpipe -> 24.0      FragPipe-24.0.zip -> 24.0
# Read from a NAME, like diann_version_of. The separator has to be `-` or `_`, so the
# request-keyed cache folder (fragpipe/24.0/, fragpipe/latest/) is deliberately not matched:
# it is named after what was ASKED for, and can hold another release.
fragpipe_version_of() {
  printf '%s' "$1" | sed -nE 's#.*[Ff]rag[Pp]ipe[-_]v?([0-9]+(\.[0-9]+)+).*#\1#p'
}

# alphadia_pip_version: the version pip reports for the INSTALLED alphadia distribution, or
# nothing. `pip show` reads the distribution metadata; it does not import alphadia (which
# loads torch -- slow, and not something to do on a login node). Printed VERBATIM: if pip
# says "1.10.0rc1", that is what tools.json records, and run_search.py -- which accepts only
# something shaped like a version -- records no version for the search rather than rounding
# a pre-release down to a release that was not what ran.
alphadia_pip_version() {
  local p out
  for p in pip pip3; do
    have "$p" || continue
    out="$("$p" show alphadia 2>/dev/null | sed -nE 's/^Version: *(.+)$/\1/p' | head -n1)"
    if [ -n "$out" ]; then printf '%s' "$out"; return; fi
  done
}

# radiant_tag_version <image ref>: the release an image TAG names
# (seerbio/radiant-fulcrum:2.3.3 -> 2.3.3), and nothing for `:latest`.
radiant_tag_version() { printf '%s' "$1" | sed -nE 's#.*:v?([0-9]+(\.[0-9]+)+)$#\1#p'; }

# ---------------------------------------------------------------- Sage --------
# sage_release_of <cache dir>: the release a cached Sage came from, read from the tarball
# acquire_sage keeps beside the binary (its top folder is sage-v0.14.7-<target>/). NOT
# `sage --version`: the v0.14.7 release binary prints "sage 0.14.6" (measured on HIVE
# 2026-09-16), while 0.14.7 is what resolve_defaults.py pins and what anyone can download
# again.
#
# Prints NOTHING when there is no such tarball. It used to fall back to the REQUEST
# (PIN_VERSION, i.e. the workflow manifest's own pin) -- so a pinned Sage whose tarball was
# missing recorded the pin as the build that ran. run_search.py could not catch it either:
# the cache path is sage/<pin>/sage, a request-keyed directory it deliberately does not read
# a version from, so it compared the pin with itself and recorded `mismatch: false` -- which
# claims the two were checked and agree. The empty string plus the "cannot tell which release
# this is" note below is the honest answer.
sage_release_of() {
  local dir="$1" got=""
  if [ -f "$dir/sage.tar.gz" ]; then
    got="$(tar -tzf "$dir/sage.tar.gz" 2>/dev/null | head -n1 \
           | sed -nE 's#^(\./)?sage-v?([0-9]+(\.[0-9]+)+)-.*#\2#p')"
  fi
  printf '%s' "$got"
}

acquire_sage() {
  # Prefer a Sage already on PATH (setup.sh installs sage-proteomics into the
  # conda env and activate.sh puts it on PATH). No need to re-download.
  # "env" names a source, not a build; run_search.py treats it as no version.
  if have sage; then SAGE_CMD="$(command -v sage)"; SAGE_VER="env"
    NOTES+=("Sage: using the conda-env build on PATH ($SAGE_CMD)."); return; fi
  local ver tag; ver="latest"; tag=""
  if pin_for sage; then ver="$PIN_VERSION"; tag="$(sage_tag "$PIN_VERSION")"; fi
  local dir="$ROOT/sage/$ver" bin
  bin="$(find "$dir" -name sage -type f -perm -u+x 2>/dev/null | head -n1)"
  if [ -z "$bin" ]; then
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
    [ -n "$bin" ] && chmod +x "$bin"
  fi
  [ -z "$bin" ] && return
  SAGE_CMD="$bin"
  # The cache is keyed on the REQUEST (sage/latest/ is what real installs look like), so
  # the version has to come from what is inside it, not from the folder name.
  SAGE_VER="$(sage_release_of "$dir")"
  if [ -z "$SAGE_VER" ]; then
    NOTES+=("Sage ($bin): cannot tell which release this is (no sage.tar.gz naming one beside it), so tools.json records no version. Delete $dir and re-run to re-download a known release.")
  elif [ "$ver" != "latest" ] && [ "$SAGE_VER" != "${ver#v}" ]; then
    NOTES+=("Sage: pinned '$ver' but the cached build in $dir is release $SAGE_VER.")
  fi
  NOTES+=("Sage ${SAGE_VER:-$ver}: run with --parquet; v0.14.x needs sage_protein_groups.py post-hoc for protein rollup. Telemetry off via --disable-telemetry-i-dont-want-to-improve-sage.")
}

# --------------------------------------------------------------- DIA-NN -------
acquire_diann() {
  local ver; ver="latest"; pin_for diann && ver="$PIN_VERSION"
  case "$CLASS" in
    hpc)
      # On HIVE the Proteomics Core keeps DIA-NN as NATIVE builds at
      #   /quobyte/proteomics-grp/dia-nn/build_<nnn>/diann-<version>/diann-linux
      # (2.5.1, 2.6.0, 2.6.1 on 2026-09-16) plus an older diann_2.3.0.sif. Match the
      # PINNED version exactly -- never silently substitute a different version
      # (reproducibility).
      local DN="$DIANN_HIVE_DIR" bin="" sif="" got=""
      if [ "$ver" != "latest" ]; then
        bin="$(ls -1 "$DN"/*/diann-"$ver"/diann-linux 2>/dev/null | head -n1)"
        [ -z "$bin" ] && bin="$(find "$DN" -maxdepth 3 -path "*diann-$ver*/diann-linux" -type f 2>/dev/null | head -n1)"
        [ -z "$bin" ] && sif="$(ls -1 "$DN"/*"$ver"*.sif 2>/dev/null | sort -V | tail -n1)"
      else
        bin="$(find "$DN" -maxdepth 3 -name diann-linux -type f 2>/dev/null | newest_diann)"
        [ -z "$bin" ] && sif="$(ls -1 "$DN"/*.sif 2>/dev/null | newest_diann)"
      fi
      # The version is read off what was FOUND. The request is only a fallback, and only
      # when it was concrete and the glob that matched it had it in the name.
      if [ -n "$bin" ]; then got="$(diann_version_of "$bin")"
      elif [ -n "$sif" ]; then got="$(diann_version_of "$sif")"; fi
      [ -z "$got" ] && [ "$ver" != "latest" ] && got="$ver"
      if [ -n "$bin" ]; then
        DIANN_CMD="$bin"; DIANN_VER="$got"
        [ "$ver" = "latest" ] && NOTES+=("DIA-NN: 'latest' resolved to $got, the highest version under $DN.")
        [ "$ver" != "latest" ] && [ "$got" != "$ver" ] && NOTES+=("DIA-NN: pinned '$ver' matched build $got ($bin).")
        NOTES+=("DIA-NN $got: using the HIVE Proteomics Core native build $bin. Reference invocation: $DN/run_diann_*.sbatch.")
        return
      elif [ -n "$sif" ] && have apptainer && { [ "$ver" = "latest" ] || printf '%s' "$sif" | grep -q "$ver"; }; then
        DIANN_CMD="apptainer exec --bind /quobyte:/quobyte $sif /diann-*/diann-linux"; DIANN_VER="$got"
        [ "$ver" = "latest" ] && NOTES+=("DIA-NN: 'latest' resolved to $got, the highest-versioned .sif under $DN (no native build found).")
        NOTES+=("DIA-NN $got: reusing HIVE container $sif.")
        return
      else
        NOTES+=("DIA-NN $ver NOT found on HIVE under $DN — NOT substituting a different version. Present: $(ls -d "$DN"/*.sif "$DN"/build_*/diann-* 2>/dev/null | xargs -n1 basename 2>/dev/null | tr '\n' ' '). Will fetch the pinned Academia build into your space instead (or build $ver).")
      fi ;;
    mac)
      if [ -n "${DIANN_DOCKER_IMAGE:-}" ] && have docker; then
        DIANN_CMD="docker run --rm -v \$PWD:/data ${DIANN_DOCKER_IMAGE} diann-linux"
        # build_diann_docker.sh tags the image with the version it RESOLVED
        # (proteomics-pipeline/diann:2.6.1), so the tag is the record. An image tagged
        # otherwise says nothing about what is inside; record nothing rather than a guess.
        DIANN_VER="$(printf '%s' "$DIANN_DOCKER_IMAGE" | sed -nE 's#.*:v?([0-9]+(\.[0-9]+)+)$#\1#p')"
        NOTES+=("DIA-NN: using Docker image \$DIANN_DOCKER_IMAGE on macOS (no native mac build exists).")
        if [ -z "$DIANN_VER" ]; then
          NOTES+=("DIA-NN: Docker image '$DIANN_DOCKER_IMAGE' has no version in its tag, so tools.json records no DIA-NN version. Rebuild with build_diann_docker.sh, which tags the image with the version it contains.")
        elif [ "$ver" != "latest" ] && [ "$DIANN_VER" != "$ver" ]; then
          NOTES+=("DIA-NN: pinned '$ver' but DIANN_DOCKER_IMAGE is $DIANN_VER ('$DIANN_DOCKER_IMAGE'). Build the pinned one: bash scripts/build_diann_docker.sh $ver")
        fi
        return
      fi
      NOTES+=("DIA-NN on macOS: no native build. Set DIANN_DOCKER_IMAGE to a built image (version $ver), or build from the Academia Linux zip's Dockerfile. See references/environment.md.")
      ;;
  esac
  # linux native (or fallback): fetch the pinned Academia Linux zip.
  # Resolve the URL FIRST so the cache is keyed on the real version: with
  # ver=latest, caching under "latest" would both re-download every run and
  # record a version string nobody can reproduce.
  local url; url="$(diann_asset_url "$ver" Linux)"
  if [ -z "$url" ]; then NOTES+=("DIA-NN: no Academia Linux build named '$ver' exists on https://github.com/vdemichev/DiaNN/releases (free for academic use). Check the version, or download manually."); return; fi
  local got; got="$(printf '%s' "$url" | sed -nE 's#.*/DIA-NN-([0-9][0-9.]*)-Academia-.*#\1#p')"
  [ -z "$got" ] && got="$ver"
  [ "$got" != "$ver" ] && NOTES+=("DIA-NN: '$ver' resolved to build $got ($(basename "$url")).")
  local dir="$ROOT/diann/$got"; mkdir -p "$dir"
  local existing; existing="$(find "$dir" -name 'diann-linux' -type f 2>/dev/null | head -n1)"
  if [ -n "$existing" ]; then DIANN_CMD="$existing"; DIANN_VER="$got"; return; fi
  echo "  DIA-NN: downloading $url"
  curl -fsSL "$url" -o "$dir/diann.zip" && (cd "$dir" && unzip -oq diann.zip)
  existing="$(find "$dir" -name 'diann-linux' -type f 2>/dev/null | head -n1)"
  if [ -n "$existing" ]; then chmod +x "$existing"; DIANN_CMD="$existing"; DIANN_VER="$got"
    NOTES+=("DIA-NN $got: Linux native needs glibc>=Mint 21.2 + .NET 8. If missing, prefer Docker/Apptainer (the zip ships a Dockerfile).")
  fi
}

# ------------------------------------------------- Radiant DIA + Fulcrum ------
# Seer ships Radiant only as a container image (multi-arch: linux/amd64 +
# linux/arm64), so there is no native binary to fetch on any platform. We record
# the RUNTIME and IMAGE separately because run_search.py has to inject bind mounts
# for the inputs, and docker (-v) and apptainer (--bind) spell those differently.
acquire_radiant() {
  local ver; ver="latest"; pin_for radiant && ver="$PIN_VERSION"
  local image="seerbio/radiant-fulcrum:$ver"
  # RADIANT_VER stays "" until an image is chosen, and is the release that image IS (see
  # "versions" in the header). An unpinned Docker/Apptainer image is `:latest`, which names
  # no release; finding out which one it points at means pulling ~3 GB, so it records "".
  local unpinned_note="Radiant: the image is seerbio/radiant-fulcrum:latest, which names no release, so tools.json records no Radiant version. To record one, pin it: PIN_ENGINE=radiant PIN_VERSION=<release> (resolve_defaults.py pins 2.3.3)."

  # LICENSE: Apache-2.0 with a MANDATORY GRANT-BACK and the COMMONS CLAUSE, which
  # removes the right to "Sell" the software or a service whose value derives
  # substantially from it. That is a real question for a fee-for-service core --
  # surface it rather than letting a run imply it was cleared.
  NOTES+=("Radiant/Fulcrum LICENSE: Apache-2.0 + Commons Clause + mandatory grant-back (https://github.com/seerbio/radiant-fulcrum-container/blob/main/LICENSE.md). The Commons Clause restricts SELLING a service whose value derives substantially from the software -- confirm with your institution before using it for fee-for-service work. DIA-NN Academia and FragPipe have their own separate terms.")

  case "$CLASS" in
    hpc)
      # Prefer an existing .sif; do NOT auto-build one (a docker->apptainer
      # conversion pulls ~3 GB and needs a writable cache -- not a login-node job).
      # Walked as a string rather than split into an array: on bash 3.2 -- macOS's system
      # bash, and the one these tests run under -- expanding "${arr[@]}" of an EMPTY array
      # under `set -u` is an unbound-variable error that aborts the whole run, which is
      # exactly what an empty RADIANT_HIVE_DIRS now produces.
      local RS sif="" rest="$RADIANT_HIVE_DIRS"
      while [ -n "$rest" ]; do
        RS="${rest%%:*}"
        case "$rest" in *:*) rest="${rest#*:}" ;; *) rest="" ;; esac
        [ -n "$RS" ] && [ -d "$RS" ] || continue
        if [ "$ver" != "latest" ]; then
          sif="$(ls -1 "$RS"/*radiant*"$ver"*.sif 2>/dev/null | sort -V | tail -n1)"
        else
          sif="$(ls -1 "$RS"/*radiant*fulcrum*.sif 2>/dev/null | sort -V | tail -n1)"
        fi
        [ -n "$sif" ] && break
      done
      if [ -n "$sif" ] && have apptainer; then
        RADIANT_CMD="apptainer exec"; RADIANT_RUNTIME="apptainer"; RADIANT_IMAGE="$sif"
        # The .sif's own name, as the build step below writes it (radiant-fulcrum-<ver>.sif);
        # the pin only when the glob that found it had the pin in the name.
        RADIANT_VER="$(basename "$sif" | sed -nE 's#.*[-_]v?([0-9]+(\.[0-9]+)+)\.sif$#\1#p')"
        [ -z "$RADIANT_VER" ] && [ "$ver" != "latest" ] && RADIANT_VER="${ver#v}"
        [ -z "$RADIANT_VER" ] && NOTES+=("Radiant: $sif names no release, so tools.json records no Radiant version.")
        NOTES+=("Radiant ${RADIANT_VER:-$ver}: reusing HIVE container $sif.")
        return
      fi
      local searched="${RADIANT_HIVE_DIRS//:/ or }"
      [ -z "$searched" ] && searched="(none: RADIANT_HIVE_DIRS is empty)"
      NOTES+=("Radiant $ver: no .sif found under $searched. Build one ON A COMPUTE NODE (never the login node): 'srun -c 8 --mem 16G --pty apptainer build radiant-fulcrum-$ver.sif docker://$image', then re-run acquire_tools.sh.")
      ;;
    mac|linux)
      # The release recorded is read from the image REFERENCE -- the thing that is stored and
      # run -- not from $ver. The two are the same string today, but the tag is evidence about
      # what will be pulled, while the pin is only the request, and nothing here can turn one
      # into the other. Nothing has pulled the image yet either, so the tag is the ONLY
      # evidence of what is inside it; the note says so.
      local tag_note="The release recorded is the image TAG; nothing has pulled the image here, so the tag is the only evidence of what is inside it."
      if have docker; then
        RADIANT_CMD="docker run --rm"; RADIANT_RUNTIME="docker"; RADIANT_IMAGE="$image"
        RADIANT_VER="$(radiant_tag_version "$RADIANT_IMAGE")"
        [ -z "$RADIANT_VER" ] && NOTES+=("$unpinned_note")
        NOTES+=("Radiant $ver: using Docker image $image (multi-arch; runs natively on Apple Silicon and x86). First run pulls ~3 GB. $tag_note")
        return
      fi
      if have apptainer; then
        RADIANT_CMD="apptainer exec"; RADIANT_RUNTIME="apptainer"; RADIANT_IMAGE="docker://$image"
        RADIANT_VER="$(radiant_tag_version "$RADIANT_IMAGE")"
        [ -z "$RADIANT_VER" ] && NOTES+=("$unpinned_note")
        NOTES+=("Radiant $ver: no Docker; will let Apptainer pull $image on first use (~3 GB). $tag_note")
        return
      fi
      NOTES+=("Radiant $ver: needs Docker or Apptainer -- Seer ships no native binary. Install Docker (https://docs.docker.com/get-docker/) and re-run.")
      ;;
  esac
}

# ------------------------------------------------------------- FragPipe -------
acquire_fragpipe() {
  local ver tag; ver="latest"; tag=""
  if pin_for fragpipe; then ver="$PIN_VERSION"; tag="$PIN_VERSION"; fi
  local dir="$ROOT/fragpipe/$ver"
  local fp asset=""; fp="$(find "$dir" -name 'fragpipe' -type f 2>/dev/null | head -n1)"
  if [ -z "$fp" ]; then
    mkdir -p "$dir"
    local url; url="$(asset_url_tag Nesvilab/FragPipe "$tag" 'FragPipe-[0-9].*zip' | grep -v jre)"
    if [ -z "$url" ]; then NOTES+=("FragPipe: could not resolve release zip for version '$ver'; see https://github.com/Nesvilab/FragPipe/releases"); return; fi
    echo "  FragPipe: downloading $url"
    # The asset names the release (FragPipe-24.0.zip). Keep it: the file is saved as
    # fragpipe.zip, so after this line the name is the only place the release survives if
    # the zip unpacks into a folder that does not carry it.
    asset="$(basename "$url")"
    curl -fsSL "$url" -o "$dir/fragpipe.zip" && (cd "$dir" && unzip -oq fragpipe.zip)
    fp="$(find "$dir" -path '*/bin/fragpipe' -type f 2>/dev/null | head -n1)"
    [ -n "$fp" ] && chmod +x "$fp"
  fi
  # The release is read from what was FOUND -- the unpacked folder (fragpipe-24.0/bin/fragpipe),
  # else the asset just downloaded. NEVER $ver: that is the request, and "latest" names no
  # release. Without this, tools.json had no `fragpipe` version at all and every FragPipe
  # search recorded `version: null` for ever.
  if [ -n "$fp" ]; then
    FRAGPIPE_VER="$(fragpipe_version_of "$fp")"
    [ -z "$FRAGPIPE_VER" ] && [ -n "$asset" ] && FRAGPIPE_VER="$(fragpipe_version_of "$asset")"
    [ -z "$FRAGPIPE_VER" ] && NOTES+=("FragPipe ($fp): cannot tell which release this is (neither the install path nor the downloaded asset names one), so tools.json records no FragPipe version.")
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
      # Whatever the .sif's NAME carries (alphadia-1.10.0.sif). Not `pip show`: that answers
      # for the login environment's alphadia, not the one inside the image.
      ALPHADIA_VER="$(basename "$asif" | sed -nE 's#.*[-_]v?([0-9]+(\.[0-9]+)+)\.sif$#\1#p')"
      NOTES+=("AlphaDIA: reusing the HIVE container $asif (Apache-2.0, commercial-OK).")
      [ -z "$ALPHADIA_VER" ] && NOTES+=("AlphaDIA: $asif names no release, so tools.json records no AlphaDIA version.")
      return
    fi
  fi
  if have alphadia; then ALPHADIA_CMD="$(command -v alphadia)"
    ALPHADIA_VER="$(alphadia_pip_version)"
    NOTES+=("AlphaDIA: using the alphadia on PATH ($ALPHADIA_CMD). Apache-2.0 (commercial use OK).")
    [ -z "$ALPHADIA_VER" ] && NOTES+=("AlphaDIA ($ALPHADIA_CMD): pip reports no version for an installed alphadia distribution, so tools.json records no AlphaDIA version.")
    return; fi
  local spec="alphadia"; pin_for alphadia && spec="alphadia==$PIN_VERSION"
  if have pip || have pip3; then
    echo "  AlphaDIA: pip install $spec (large download; GPU recommended)..."
    { pip install "$spec" >/dev/null 2>&1 || pip3 install "$spec" >/dev/null 2>&1; }
    if have alphadia; then ALPHADIA_CMD="$(command -v alphadia)"
      ALPHADIA_VER="$(alphadia_pip_version)"
      NOTES+=("AlphaDIA ${ALPHADIA_VER:-(version unknown)}: installed via pip (Apache-2.0, commercial-OK). GPU strongly recommended; verify with 'alphadia --check'."); return; fi
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
# Radiant: a ~3 GB image pull, so only when pinned/requested or already available.
if [ "${PIN_ENGINE:-}" = "radiant" ] || [ "${ACQUIRE_RADIANT:-}" = "1" ]; then
  acquire_radiant
fi

# ---- write manifest ----------------------------------------------------------
esc() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }
{
  printf '{\n'
  printf '  "platform_class": "%s",\n' "$CLASS"
  printf '  "tools_root": "%s",\n' "$(esc "$ROOT")"
  printf '  "pinned": {"engine": "%s", "version": "%s"},\n' "${PIN_ENGINE:-}" "${PIN_VERSION:-}"
  printf '  "versions": {"diann": "%s", "sage": "%s", "radiant": "%s", "fragpipe": "%s", "alphadia": "%s"},\n' "$(esc "$DIANN_VER")" "$(esc "$SAGE_VER")" "$(esc "${RADIANT_VER:-}")" "$(esc "${FRAGPIPE_VER:-}")" "$(esc "${ALPHADIA_VER:-}")"
  printf '  "radiant":  %s,\n'  "$( [ "${RADIANT_CMD:-null}" = null ]  && echo null || printf '"%s"' "$(esc "${RADIANT_CMD}")" )"
  printf '  "radiant_runtime": %s,\n' "$( [ -z "${RADIANT_RUNTIME:-}" ] && echo null || printf '"%s"' "$(esc "${RADIANT_RUNTIME}")" )"
  printf '  "radiant_image":   %s,\n' "$( [ -z "${RADIANT_IMAGE:-}" ]   && echo null || printf '"%s"' "$(esc "${RADIANT_IMAGE}")" )"
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
