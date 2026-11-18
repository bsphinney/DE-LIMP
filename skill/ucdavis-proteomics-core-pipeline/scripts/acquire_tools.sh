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
#          PIN_ENGINE=diann PIN_VERSION=2.7.0 bash acquire_tools.sh hpc   # the pin resolve_defaults.py writes
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
#            silently downloaded. Where this account can run the Core's licensed
#            install (engine_choice.CORE_FRAGPIPE_INSTALL, HIVE), `fragpipe` IS its
#            bin/fragpipe. tools.json `fragpipe_tools` (engine_choice.fragpipe_tools:
#            the tools folder, its jars, Philosopher + bundled DIA-NN of the install
#            that runs, the Java, `licensed`, `why`) is what the search reads.
# =============================================================================
set -uo pipefail

ACQUIRE_ARGS=("$@")              # for a re-run on a compute node (download_here)
CLASS="${1:?platform_class required (hpc|mac|linux)}"
ROOT="${2:-$HOME/.proteomics-pipeline/tools}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# DIA-NN resolves by asset filename, not by release tag -- see diann_release.sh
# for why (every 2.x build hangs off the single release tagged "2.0").
. "$SCRIPT_DIR/diann_release.sh"
# sha256 of what is downloaded, against the vendor's published value (verify_download.sh)
. "$SCRIPT_DIR/verify_download.sh"
PIN_ENGINE="${PIN_ENGINE:-}"      # diann | sage | fragpipe  (optional)
PIN_VERSION="${PIN_VERSION:-}"    # exact version to honor   (optional)
mkdir -p "$ROOT"
MANIFEST="$ROOT/tools.json"
OS="$(uname -s | tr '[:upper:]' '[:lower:]')"; ARCH="$(uname -m)"
have() { command -v "$1" >/dev/null 2>&1; }

# A download (Sage, DIA-NN, FragPipe: a release archive fetched and unpacked) is work for a
# compute node -- HPC@UCD's ucd-hpc rules; finding a build that is already there is not. On a
# cluster login node (on_compute.py --is-login-node) the first download re-runs this whole script
# there through on_compute.py: one blocking srun, and tools.json is written by that run. With no
# python3 + on_compute.py on a cluster, nothing is downloaded and a note says how.
# download_here <what>: returns 0 to download here; 1 = do not.
download_here() {
  if [ -n "${SLURM_JOB_ID:-}" ] || [ -n "${PROTEOMICS_ON_COMPUTE:-}" ]; then return 0; fi
  if have python3 && [ -f "$SCRIPT_DIR/on_compute.py" ]; then
    python3 "$SCRIPT_DIR/on_compute.py" --is-login-node || return 0
    exec python3 "$SCRIPT_DIR/on_compute.py" --time 60 --mem 4G --name pp_acquire \
      --why "acquire_tools.sh downloads $1" -- "${BASH:-bash}" "$SCRIPT_DIR/acquire_tools.sh" "${ACQUIRE_ARGS[@]}"
  fi
  if have sbatch || [ -x "${HPCCF_SLURM_BIN:-/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin}/sbatch" ] \
     || [ -r "${UCD_HPC_CLUSTER_CONF:-/opt/hpccf/etc/cluster_name.conf}" ]; then
    NOTES+=("$1: not downloaded -- this is a cluster login node with no python3 + on_compute.py to move the download to a compute node. Run acquire_tools.sh in a job: srun --account=<yours> --partition=<yours> --time=60 --mem=4G bash $SCRIPT_DIR/acquire_tools.sh ${ACQUIRE_ARGS[*]}")
    return 1
  fi
  return 0
}

# normalize a version into the tags engines actually use
sage_tag()  { local v="$1"; [ -z "$v" ] && { echo ""; return; }; case "$v" in v*) echo "$v";; *) echo "v$v";; esac; }

# Asset URL from a SPECIFIC release tag (pinned); falls back to latest if tag is blank.
# Prints nothing when it cannot answer, and gh_fetch (diann_release.sh) leaves the reason in
# $GH_ERR_FILE -- "GitHub could not be reached" and "no such release asset" are different
# problems and used to read the same.
asset_url_tag() { # owner/repo  tag  pattern
  local repo="$1" tag="$2" pat="$3" api body
  if [ -n "$tag" ]; then api="https://api.github.com/repos/$repo/releases/tags/$tag"
  else                  api="https://api.github.com/repos/$repo/releases/latest"; fi
  body="$(gh_fetch "$api")" || return 0
  printf '%s\n' "$body" \
    | grep '"browser_download_url"' | sed -E 's/.*"(https[^"]+)".*/\1/' \
    | grep -iE "$pat" | head -n1
}

# gh_why: the reason the last GitHub lookup failed ("" when it did not), for a note.
gh_why() { [ -s "${GH_ERR_FILE:-/nonexistent}" ] && cat "$GH_ERR_FILE"; return 0; }
gh_why_paren() { local w; w="$(gh_why)"; [ -n "$w" ] && printf ' (%s)' "$w"; return 0; }

# fetch_archive <url> <dest file> <checksum key> <label>: download to a temporary name, check
# the sha256 against the vendor's published value, and only then rename it into place. A
# download cut short, or one that does not match, never leaves a file at <dest> -- the old
# `curl -o <dest> && unpack` left a half-written archive that the next run reused.
# Returns 0 = in place (verified, or "not verified" said aloud when the vendor publishes no
# checksum); non-zero = nothing at <dest>.
fetch_archive() {
  local url="$1" dest="$2" key="$3" label="$4" part="$2.part" rc
  rm -f "$part"
  if ! curl -fsSL "$url" -o "$part"; then
    rm -f "$part"; NOTES+=("$label: the download of $url failed or was cut short, so nothing was installed. Check the network and re-run."); return 1
  fi
  verify_download "$key" "$part" "$label" || { rc=$?
    if [ "$rc" = 1 ]; then rm -f "$part"; NOTES+=("$label: the download did NOT match the sha256 the vendor published, so it was deleted and not used. Re-run; if it happens again, report it (references/environment.md)."); return 1; fi
    NOTES+=("$label: downloaded but NOT verified (no published checksum is pinned for '$key' in scripts/pinned_checksums.txt).")
  }
  mv -f "$part" "$dest"
}

DIANN_CMD="null"; SAGE_CMD="null"; FRAGPIPE_CMD="null"; ALPHADIA_CMD="null"; NOTES=()
GH_ERR_FILE="$(mktemp 2>/dev/null || true)"; export GH_ERR_FILE
trap '[ -n "$GH_ERR_FILE" ] && rm -f "$GH_ERR_FILE"' EXIT
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
# The Core's AlphaDIA container; overridable like the folders above (tests point it at nothing).
ALPHADIA_HIVE_SIF="${ALPHADIA_HIVE_SIF-/quobyte/proteomics-grp/apptainers/alphadia.sif}"

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

# sage_path_version <sage on PATH>: the release a Sage on PATH IS, or nothing. A conda install
# says it in the package record next to the binary (conda-meta/sage-proteomics-<ver>-<build>.json);
# failing that, `sage --version`. NOT trusted alone: the v0.14.7 release binary prints
# "sage 0.14.6" (measured on HIVE 2026-09-16), so a mismatch here means "cannot confirm the pin",
# not "wrong build" -- the caller then fetches the pinned release instead of guessing.
# sage_conda_version <sage on PATH>: the conda package version recorded beside it
# (conda-meta/sage-proteomics-<ver>-<build>.json), or nothing. The one reader of that record.
# bioconda's sage-proteomics 0.14.7 builds the v0.14.7 tag, whose crates/sage-cli/Cargo.toml says
# 0.14.6 -- so its `sage --version` prints 0.14.6 too (an upstream self-report bug; HIVE F8).
sage_conda_version() {
  local bin real prefix rec v=""
  bin="$1"; real="$bin"
  if have readlink; then real="$(readlink "$bin" 2>/dev/null || true)"; [ -n "$real" ] || real="$bin"
    case "$real" in /*) ;; *) real="$(dirname "$bin")/$real" ;; esac; fi
  prefix="$(cd "$(dirname "$real")/.." 2>/dev/null && pwd)"
  for rec in "$prefix"/conda-meta/sage-proteomics-*.json; do
    [ -f "$rec" ] || continue
    v="$(basename "$rec" | sed -nE 's#^sage-proteomics-([0-9]+(\.[0-9]+)+)-.*#\1#p')"
    [ -n "$v" ] && break
  done
  printf '%s' "$v"
}
sage_path_version() {
  local v; v="$(sage_conda_version "$1")"
  [ -z "$v" ] && v="$("$1" --version 2>/dev/null | sed -nE 's#^[Ss]age[ v]*([0-9]+(\.[0-9]+)+).*#\1#p' | head -n1)"
  printf '%s' "$v"
}

acquire_sage() {
  # Prefer a Sage already on PATH (setup.sh installs sage-proteomics into the
  # conda env and activate.sh puts it on PATH). No need to re-download.
  # "env" names a source, not a build; run_search.py treats it as no version.
  # BUT a pin is a pin: when a version is pinned, a Sage on PATH is used only if it can be shown to
  # be that version (sage_path_version). A different one used to win silently, because this
  # branch returned before the pin was even read, and tools.json recorded "env".
  if have sage; then
    local onpath; onpath="$(command -v sage)"
    if ! pin_for sage; then
      # the package the conda env records is the build it IS (2.15.1, HIVE F8: the setup env's
      # sage-proteomics 0.14.7 was recorded as "env", and its `sage --version` says 0.14.6);
      # with no package record, "env": a source, not a build
      local cv; cv="$(sage_conda_version "$onpath")"
      SAGE_CMD="$onpath"; SAGE_VER="${cv:-env}"
      NOTES+=("Sage${cv:+ $cv (conda package sage-proteomics)}: using the conda-env build on PATH ($SAGE_CMD).")
      return
    fi
    local pv; pv="$(sage_path_version "$onpath")"
    if [ -n "$pv" ] && [ "$pv" = "${PIN_VERSION#v}" ]; then
      SAGE_CMD="$onpath"; SAGE_VER="$pv"
      NOTES+=("Sage $pv: using the build on PATH ($SAGE_CMD); it matches the pin.")
      return
    fi
    NOTES+=("Sage: the sage on PATH ($onpath) ${pv:+is $pv, not the pinned ${PIN_VERSION#v}}${pv:-cannot be shown to be the pinned ${PIN_VERSION#v} (its version is unknown)}, so it is NOT used; the pinned release is fetched into $ROOT/sage/$PIN_VERSION instead.")
  fi
  local ver tag; ver="latest"; tag=""
  if pin_for sage; then ver="$PIN_VERSION"; tag="$(sage_tag "$PIN_VERSION")"; fi
  local dir="$ROOT/sage/$ver" bin
  bin="$(find "$dir" -name sage -type f -perm -u+x 2>/dev/null | head -n1)"
  if [ -z "$bin" ]; then
    mkdir -p "$dir"
    local pat target
    case "$OS-$ARCH" in
      linux-x86_64)  pat="x86_64-unknown-linux-gnu.*tar"; target="x86_64-unknown-linux-gnu" ;;
      darwin-arm64)  pat="aarch64-apple-darwin.*tar";     target="aarch64-apple-darwin" ;;
      darwin-x86_64) pat="x86_64-apple-darwin.*tar";      target="x86_64-apple-darwin" ;;
      *)             pat="x86_64-unknown-linux-gnu.*tar"; target="x86_64-unknown-linux-gnu" ;;
    esac
    : > "$GH_ERR_FILE"
    local url; url="$(asset_url_tag lazear/sage "$tag" "$pat")"
    # A pinned release has a fixed, known address; when GitHub's API cannot answer (offline,
    # rate limited) use it rather than report "no such release".
    if [ -z "$url" ] && [ -n "$tag" ] && [ -s "$GH_ERR_FILE" ]; then
      url="https://github.com/lazear/sage/releases/download/$tag/sage-$tag-$target.tar.gz"
      NOTES+=("Sage: $(gh_why); using the release's fixed address $url instead.")
    fi
    if [ -z "$url" ]; then NOTES+=("Sage: could not resolve a release asset for $OS-$ARCH version='$ver'$(gh_why_paren); see https://github.com/lazear/sage/releases"); return; fi
    download_here Sage || return
    echo "  Sage: downloading $url"
    # Into a staging folder, renamed into place last: a tarball cut short used to leave a
    # half-written `sage` that the next run found and reused. Sage publishes no checksum, so the
    # sha256 step says "not verified" (pinned_checksums.txt has no sage line) -- but the archive
    # is read to its end (tar -tzf) before anything is unpacked.
    local stage="$dir/.stage.$$"
    rm -rf "$stage"; mkdir -p "$stage"
    if fetch_archive "$url" "$stage/sage.tar.gz" "sage-${ver#v}-$target" "Sage ${ver#v}" \
       && tar -tzf "$stage/sage.tar.gz" >/dev/null 2>&1 \
       && { tar -xzf "$stage/sage.tar.gz" -C "$stage" --strip-components=1 2>/dev/null || tar -xzf "$stage/sage.tar.gz" -C "$stage"; }; then
      local got; got="$(find "$stage" -name sage -type f 2>/dev/null | head -n1)"
      if [ -n "$got" ]; then
        chmod +x "$got"
        mv -f "$stage/sage.tar.gz" "$dir/sage.tar.gz"      # the release record first
        mv -f "$got" "$dir/sage"                           # the binary last: it is the cache test
        bin="$dir/sage"
      else
        NOTES+=("Sage: the downloaded archive holds no 'sage' program, so nothing was installed.")
      fi
    elif [ -f "$stage/sage.tar.gz" ]; then
      NOTES+=("Sage: the downloaded archive is damaged (it cannot be read to its end), so nothing was installed.")
    fi
    rm -rf "$stage"
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
      # (2.5.1, 2.6.0, 2.6.1 on 2026-09-16; 2.7.0 added 2026-09-23) plus an older
      # diann_2.3.0.sif. Match the PINNED version exactly -- never silently substitute a
      # different version (reproducibility).
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
        # diann_docker.sh mounts every absolute path of the command at the SAME path inside the
        # container and adds --platform linux/amd64. The old `docker run -v $PWD:/data` mounted one
        # folder at /data, so the absolute --f/--fasta/--out paths the search passes did not exist
        # in the container (audit 2026-11-17, finding 2).
        DIANN_CMD="bash $(printf '%q' "$SCRIPT_DIR/diann_docker.sh") ${DIANN_DOCKER_IMAGE}"
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
  : > "$GH_ERR_FILE"
  local url; url="$(diann_asset_url "$ver" Linux)"
  # GitHub's API could not answer (offline, rate limited): that is not "no such build". A pinned
  # 2.x build hangs off the release tagged 2.0 (diann_release.sh), so it has a fixed address.
  if [ -z "$url" ] && [ "$ver" != "latest" ] && [ -s "$GH_ERR_FILE" ]; then
    url="https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-$ver-Academia-Linux.zip"
    NOTES+=("DIA-NN: $(gh_why); using the fixed address $url for the pinned build (its sha256 is checked).")
  fi
  if [ -z "$url" ]; then
    if [ -s "$GH_ERR_FILE" ]; then NOTES+=("DIA-NN: could not look up '$ver': $(gh_why). Nothing was downloaded; try again, or set GITHUB_TOKEN.")
    else NOTES+=("DIA-NN: no Academia Linux build named '$ver' exists on https://github.com/vdemichev/DiaNN/releases (free for academic use). Check the version, or download manually."); fi
    return
  fi
  local got; got="$(printf '%s' "$url" | sed -nE 's#.*/DIA-NN-([0-9][0-9.]*)-Academia-.*#\1#p')"
  [ -z "$got" ] && got="$ver"
  [ "$got" != "$ver" ] && NOTES+=("DIA-NN: '$ver' resolved to build $got ($(basename "$url")).")
  local dir="$ROOT/diann/$got"; mkdir -p "$dir"
  local existing; existing="$(find "$dir" -name 'diann-linux' -type f 2>/dev/null | head -n1)"
  if [ -n "$existing" ]; then DIANN_CMD="$existing"; DIANN_VER="$got"; return; fi
  download_here DIA-NN || return
  echo "  DIA-NN: downloading $url"
  # verified against the vendor's published sha256, then unpacked into a staging folder and
  # moved into place: an unzip that stops half-way no longer leaves a diann-linux behind.
  local stage="$dir/.stage.$$"; rm -rf "$stage"; mkdir -p "$stage"
  if fetch_archive "$url" "$stage/diann.zip" "diann-$got-linux" "DIA-NN $got" && (cd "$stage" && unzip -oq diann.zip); then
    rm -f "$stage/diann.zip"
    # move whatever it unpacked into place (a folder named for the release, or loose files)
    ( cd "$stage" && for f in * .[!.]*; do [ -e "$f" ] && { rm -rf "${dir:?}/$f"; mv -f "$f" "$dir/"; }; done ) 2>/dev/null
  else
    [ -f "$stage/diann.zip" ] && NOTES+=("DIA-NN $got: the zip could not be unpacked (damaged or incomplete), so nothing was installed.")
  fi
  rm -rf "$stage"
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
      # HPC@UCD's ucd-hpc rules: `apptainer build` once, in a job, every request explicit
      # (on_compute.py spells out the queue), and the ~3 GB pull's cache off the 20 GB home --
      # the job's own $TMPDIR. The .sif goes where this script looks for it.
      local into="" d rest2="$RADIANT_HIVE_DIRS"
      while [ -n "$rest2" ]; do
        d="${rest2%%:*}"; case "$rest2" in *:*) rest2="${rest2#*:}" ;; *) rest2="" ;; esac
        [ -n "$d" ] && [ -d "$d" ] && { into="$d"; break; }
      done
      local again="re-run acquire_tools.sh"
      [ -z "$into" ] && { into="<a folder on your group storage>"; again="re-run acquire_tools.sh with RADIANT_HIVE_DIRS=<that folder>"; }
      NOTES+=("Radiant $ver: no .sif found under $searched. Build one in a job, never on the login node: python3 $SCRIPT_DIR/on_compute.py --time 120 --mem 16G --cpus 8 -- bash -lc 'APPTAINER_CACHEDIR=\"\$TMPDIR/apptainer-cache\" apptainer build $into/radiant-fulcrum-$ver.sif docker://$image', then $again.")
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
  local fp="" asset="" core=""
  # The Core's licensed install, when this account can run it and it is the pinned release
  # (engine_choice.py core-binary: the one definition of where it is). Its own tools/ holds the
  # jars, and FragPipe takes Philosopher from the install that runs -- so on HIVE the home copy
  # (no jars) is never what runs. No python3 = not looked for.
  if have python3; then
    core="$(python3 "$SCRIPT_DIR/engine_choice.py" core-binary ${tag:+--pin "$tag"} 2>/dev/null || true)"
  fi
  if [ -n "$core" ]; then
    fp="$core"
    NOTES+=("FragPipe: using the Core's install $core.")
  else
    fp="$(find "$dir" -name 'fragpipe' -type f 2>/dev/null | head -n1)"
  fi
  if [ -z "$fp" ]; then
    mkdir -p "$dir"
    : > "$GH_ERR_FILE"
    local url; url="$(asset_url_tag Nesvilab/FragPipe "$tag" 'FragPipe-[0-9].*zip' | grep -v jre)"
    if [ -z "$url" ] && [ -n "$tag" ] && [ -s "$GH_ERR_FILE" ]; then
      url="https://github.com/Nesvilab/FragPipe/releases/download/$tag/FragPipe-$tag-linux.zip"
      NOTES+=("FragPipe: $(gh_why); using the release's fixed address $url (its sha256 is checked).")
    fi
    if [ -z "$url" ]; then NOTES+=("FragPipe: could not resolve release zip for version '$ver'$(gh_why_paren); see https://github.com/Nesvilab/FragPipe/releases"); fragpipe_tools_record ""; return; fi
    download_here FragPipe || { fragpipe_tools_record ""; return; }
    echo "  FragPipe: downloading $url"
    # The asset names the release (FragPipe-24.0.zip). Keep it: the file is saved as
    # fragpipe.zip, so after this line the name is the only place the release survives if
    # the zip unpacks into a folder that does not carry it.
    asset="$(basename "$url")"
    # verified against the vendor's published sha256 first (FragPipe-24.0-linux.zip is pinned in
    # pinned_checksums.txt); a download that does not match is deleted, not unpacked.
    local fpkey="fragpipe-${tag}-linux"
    if fetch_archive "$url" "$dir/fragpipe.zip" "$fpkey" "FragPipe ${tag:-latest}"; then
      (cd "$dir" && unzip -oq fragpipe.zip)
    fi
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
  if ! have java; then NOTES+=("FragPipe: requires Java 11+ on PATH (FragPipe 24.0). Install a JDK first."); fi
  fragpipe_tools_record "$fp" "$tag"
  [ -z "$fp" ] && return
  FRAGPIPE_CMD="$fp"
  case "$FRAGPIPE_LICENSED" in
    true)
      NOTES+=("FragPipe $ver: ready (tools folder $FRAGPIPE_FOLDER, tools.json fragpipe_tools.folder; the search passes it to FragPipe as --config-tools-folder).") ;;
    false)
      NOTES+=("FragPipe $ver: downloaded, but MSFragger/IonQuant are license-gated and absent. Run the GUI once to accept the academic license, OR set FRAGPIPE_TOOLS_FOLDER. Headless search disabled until then.") ;;
    *)
      NOTES+=("FragPipe $ver: whether MSFragger/IonQuant are here was not checked (no python3 to run engine_choice.py), so tools.json fragpipe_tools.licensed is null. Re-run acquire_tools.sh where python3 is on PATH.") ;;
  esac
}

# fragpipe_tools_record <bin/fragpipe or ""> [pin]: tools.json `fragpipe_tools`, from
# engine_choice.fragpipe_tools() -- the one place the tools folder, the jars (two levels down,
# any case), Philosopher, the bundled DIA-NN, the Java and `licensed` are worked out. It reads
# $FRAGPIPE_TOOLS_FOLDER itself; the search reads only what is written here (review M4).
FRAGPIPE_TOOLS_JSON=""; FRAGPIPE_LICENSED=""; FRAGPIPE_FOLDER=""
fragpipe_tools_record() {
  local out=""
  if have python3; then
    out="$(python3 "$SCRIPT_DIR/engine_choice.py" fragpipe-tools ${1:+--binary "$1"} ${2:+--pin "$2"} 2>/dev/null || true)"
  fi
  if [ -z "$out" ]; then
    FRAGPIPE_TOOLS_JSON='{"binary": null, "folder": null, "source": null, "jars": null, "philosopher": null, "diann": null, "java": null, "licensed": null, "ready": null, "why": "not checked: python3 is not on PATH where acquire_tools.sh ran", "skipped": []}'
    FRAGPIPE_LICENSED=""; FRAGPIPE_FOLDER=""
    return
  fi
  FRAGPIPE_TOOLS_JSON="$out"
  FRAGPIPE_LICENSED="$(printf '%s\n' "$out" | sed -nE 's/^  "licensed": (true|false|null),?$/\1/p' | head -n1)"
  FRAGPIPE_FOLDER="$(printf '%s\n' "$out" | sed -nE 's/^  "folder": "(.*)",?$/\1/p' | head -n1)"
}

# ------------------------------------------------------------- AlphaDIA -------
# Apache-2.0 (commercial use OK) — the open-source DIA alternative to DIA-NN,
# whose free "Academia" build is academic/non-profit only. pip-installed into the
# active env. Deep-learning based: a GPU is strongly recommended (CPU is slow).
acquire_alphadia() {
  # On HIVE the Proteomics Core keeps an AlphaDIA container.
  if [ "$CLASS" = "hpc" ]; then
    local asif="$ALPHADIA_HIVE_SIF"
    if [ -n "$asif" ] && [ -f "$asif" ] && have apptainer; then
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
  # On a cluster this runs on the login node (hive_exec.sh): a large download there, and a bare
  # pip with the system python installs into ~/.local -- both against HPC@UCD's ucd-hpc rules
  # (never pip into a central python or --user). The note says how: the pipeline env's own pip,
  # in a job.
  if [ "$CLASS" = "hpc" ]; then
    local envpy="${PROTEOMICS_PIPELINE_HOME:-$HOME/.proteomics-pipeline}/micromamba/envs/proteomics-pipeline/bin/python"
    NOTES+=("AlphaDIA: not installed (on a cluster it is never pip-installed where this runs). Install it into the pipeline env, in a job: python3 $SCRIPT_DIR/on_compute.py --time 60 --mem 8G -- $envpy -m pip install '$spec' (Apache-2.0, commercial-OK; GPU recommended). https://github.com/MannLabs/alphadia")
    return
  fi
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
   || { [ "$CLASS" = "hpc" ] && [ -n "$ALPHADIA_HIVE_SIF" ] && [ -f "$ALPHADIA_HIVE_SIF" ]; }; then
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
  printf '  "fragpipe_tools": %s,\n' "$(printf '%s' "${FRAGPIPE_TOOLS_JSON:-null}" | sed '2,$s/^/  /')"
  printf '  "alphadia": %s,\n'  "$( [ "$ALPHADIA_CMD" = null ] && echo null || printf '"%s"' "$(esc "$ALPHADIA_CMD")" )"
  printf '  "notes": [\n'
  for i in "${!NOTES[@]}"; do
    printf '    "%s"%s\n' "$(esc "${NOTES[$i]}")" "$( [ "$i" -lt $((${#NOTES[@]}-1)) ] && echo , )"
  done
  printf '  ]\n}\n'
} | tee "$MANIFEST"

echo "[acquire] manifest written to $MANIFEST"
