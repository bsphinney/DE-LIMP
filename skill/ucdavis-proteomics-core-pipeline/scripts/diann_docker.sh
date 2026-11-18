#!/usr/bin/env bash
# =============================================================================
# diann_docker.sh  --  run DIA-NN inside its Docker image, on a machine with no native build
# (macOS). tools.json records it as the DIA-NN command:
#
#     bash /path/to/diann_docker.sh <image> <DIA-NN arguments...>
#
# What it does that a bare `docker run -v $PWD:/data <image> diann-linux` did not:
#   * Every absolute path in the arguments (--f, --fasta, --cfg, --out, --out-lib, --lib,
#     --temp, and the paths written inside a --cfg file) is made visible INSIDE the container
#     at the SAME path, by mounting the folder that holds it. The search commands pass absolute
#     host paths; with only $PWD mounted at /data, none of them existed in the container and
#     DIA-NN could read neither the raw files nor the FASTA (audit 2026-11-17, finding 2).
#   * `--platform linux/amd64`: DIA-NN ships no arm64 build, so on Apple Silicon the image runs
#     under emulation. Without the flag Docker picks the host's architecture and fails or runs the
#     wrong image.
#   * The working folder (the job's `cd <out>`) is mounted and used as the container's.
#
# Test without Docker:  DIANN_DOCKER_DRYRUN=1 bash diann_docker.sh <image> --f /a/b.d --out /c/r.parquet
# prints the docker command, one word per line. DIANN_DOCKER_BIN names the docker program.
# Works with macOS's bash 3.2.
# =============================================================================
set -uo pipefail

if [ $# -lt 1 ]; then echo "usage: diann_docker.sh <image> [DIA-NN arguments...]" >&2; exit 2; fi
IMAGE="$1"; shift
PLATFORM="${DIANN_DOCKER_PLATFORM:-linux/amd64}"
DOCKER="${DIANN_DOCKER_BIN:-docker}"

MOUNTS=()          # folders to mount, each at its own path
add_mount() {      # add_mount <path>: the nearest EXISTING folder at or above it
  local p="$1" d
  case "$p" in /*) ;; *) return 0 ;; esac
  d="$p"
  while [ -n "$d" ] && [ "$d" != "/" ] && [ ! -d "$d" ]; do d="$(dirname "$d")"; done
  # never the whole disk: a path whose only existing ancestor is / is skipped
  [ -n "$d" ] && [ "$d" != "/" ] && [ -d "$d" ] || return 0
  case "$d" in *,*) echo "diann_docker.sh: a folder name with a comma cannot be mounted: $d" >&2; exit 2 ;; esac
  # a file is mounted through its folder
  local m
  if [ ${#MOUNTS[@]} -gt 0 ]; then
    for m in "${MOUNTS[@]}"; do [ "$m" = "$d" ] && return 0; done
  fi
  MOUNTS+=("$d")
}

mount_for() {      # mount_for <path>: a file is mounted by its folder, a folder by itself
  local p="$1"
  if [ -d "$p" ]; then add_mount "$p"; else add_mount "$(dirname "$p")"; fi
}

mount_for "$PWD"
prev=""
for a in "$@"; do
  case "$a" in
    /*) mount_for "$a" ;;
  esac
  # a --cfg file names more paths inside it (--fasta, --lib, --temp ...)
  if [ "$prev" = "--cfg" ] && [ -f "$a" ]; then
    mount_for "$a"
    # its words as the cfg reader takes them (diann_parallel.cfg_tokens: bash's quoting, so a
    # path with a space is written quoted by write_cfg): xargs honours the quotes, and `read -r`
    # never glob-expands a word (2.15.1 review minor 13: `for w in $(...)` split on spaces and
    # globbed). A word with an unmatched quote ends the list: what was read so far is mounted.
    while IFS= read -r w; do
      case "$w" in /*) mount_for "$w" ;; esac
    done < <(sed 's/^[[:space:]]*#.*//' "$a" 2>/dev/null | xargs -n 1 printf '%s\n' 2>/dev/null | head -n 400)
  fi
  prev="$a"
done

ARGV=("$DOCKER" run --rm --platform "$PLATFORM" -w "$PWD")
if [ ${#MOUNTS[@]} -gt 0 ]; then
  for m in "${MOUNTS[@]}"; do ARGV+=(--mount "type=bind,source=$m,target=$m"); done
fi
ARGV+=("$IMAGE" diann-linux "$@")

if [ -n "${DIANN_DOCKER_DRYRUN:-}" ]; then printf '%s\n' "${ARGV[@]}"; exit 0; fi
exec "${ARGV[@]}"
