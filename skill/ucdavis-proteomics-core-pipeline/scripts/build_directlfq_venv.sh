#!/bin/bash
# build_directlfq_venv.sh -- build directLFQ's own pinned venv, once, in a job.
#
# directLFQ is the protein quantity of a FragPipe DDA search whose IonQuant wrote no MaxLFQ for
# every experiment (2 experiments, or MaxLFQ off; fragpipe_quant.py, Ruling 102). It never goes into the pipeline env: directLFQ 0.3.3 under pandas 3 silently
# skips its between-sample normalisation, so its venv pins pandas (fragpipe_quant.DIRECTLFQ_PINS,
# the one list: `python3 fragpipe_quant.py pins`).
#
#   python3 on_compute.py --time 30 --mem 8G --cpus 2 --why 'build directLFQ venv' -- \
#       bash build_directlfq_venv.sh [--prefix DIR] [--python PY]
#   (or: sbatch --account=... --partition=... --time=30 --mem=8G --cpus-per-task=2 \
#        --wrap "bash <this file's full path> [--prefix DIR]")
#
#   --prefix   where the venv goes; default `fragpipe_quant.py venv --default-prefix`:
#              $DIRECTLFQ_VENV, else the Core's group location on HIVE (fragpipe_quant.GROUP_VENV, a
#              folder named for its pins), else ~/.proteomics-pipeline/venvs/<name> off-cluster. Never the 20 GB home on an
#              HPC@UCD cluster (ucd-hpc): refused there.
#   --python   the base Python (>= 3.9). Default: HIVE's `module load python/3.11.9` (on CVMFS, so
#              every account can run the venv), else the skill env's python3, else python3.
#
# It refuses a login node (a pip solve is compute: ucd-hpc). It writes the marker
# (directlfq_venv.json) LAST, after `fragpipe_quant.py check-env` passes inside the new venv with a
# clean environment: a venv without the marker is never used. An existing complete venv with the
# same pins is left as it is (exit 0); a folder without a marker is refused, never removed.
# pip's cache goes inside the prefix, never ~/.cache. The venv ends group-readable and writable
# by its builder only (never world-writable): one built venv serves every staff account.
#
# Exit: 0 built (or already there); 2 usage / refused before building; 3 built but the check failed.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX=""
PY=""
while [ $# -gt 0 ]; do
  case "$1" in
    --prefix) PREFIX="${2:?--prefix needs a folder}"; shift 2 ;;
    --python) PY="${2:?--python needs a path}"; shift 2 ;;
    -h|--help) sed -n '2,29p' "$0"; exit 0 ;;
    *) echo "[build_directlfq_venv] unknown argument: $1" >&2; exit 2 ;;
  esac
done
say() { echo "[build_directlfq_venv] $*" >&2; }
if [ ! -f "$HERE/fragpipe_quant.py" ]; then
  say "fragpipe_quant.py is not beside this script ($HERE): run it from the skill's scripts/ folder by its full path (sbatch: --wrap \"bash <full path>\")."
  exit 2
fi
if python3 "$HERE/on_compute.py" --is-login-node 2>/dev/null; then
  say "this is a cluster login node, and a pip install is compute (HPC@UCD's ucd-hpc rules). Run:"
  say "  python3 $HERE/on_compute.py --time 30 --mem 8G --cpus 2 --why 'build directLFQ venv' -- bash $HERE/build_directlfq_venv.sh $*"
  exit 2
fi

if [ -z "$PY" ]; then
  if [ -f /etc/profile.d/modules.sh ]; then
    # shellcheck disable=SC1091
    set +u; source /etc/profile.d/modules.sh; set -u
    if module load python/3.11.9 >/dev/null 2>&1; then PY="$(command -v python3)"; fi
  fi
  if [ -z "$PY" ] && [ -x "$HOME/.proteomics-pipeline/micromamba/envs/proteomics-pipeline/bin/python3" ]; then
    PY="$HOME/.proteomics-pipeline/micromamba/envs/proteomics-pipeline/bin/python3"
  fi
  [ -n "$PY" ] || PY="$(command -v python3)"
fi
if ! "$PY" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 9) else 1)'; then
  say "$PY is older than Python 3.9 (directLFQ 0.3.3 needs 3.9 or later); give --python"
  exit 2
fi
[ -n "$PREFIX" ] || PREFIX="$("$PY" -I "$HERE/fragpipe_quant.py" venv --default-prefix)"
case "$PREFIX" in /*) ;; *) PREFIX="$(pwd)/$PREFIX" ;; esac
if [ -f "${UCD_HPC_CLUSTER_CONF:-/opt/hpccf/etc/cluster_name.conf}" ]; then
  case "$PREFIX/" in "$HOME"/*)
    say "REFUSED: $PREFIX is in the 20 GB home on an HPC@UCD cluster; give --prefix on group storage (or set DIRECTLFQ_VENV)."
    exit 2 ;;
  esac
fi
MARKER="$PREFIX/directlfq_venv.json"
if [ -f "$MARKER" ]; then
  if DIRECTLFQ_VENV="$PREFIX" "$PY" -I "$HERE/fragpipe_quant.py" venv >/dev/null 2>&1; then
    say "already built: $PREFIX (marker $MARKER). Nothing done."
    exit 0
  fi
  say "REFUSED: $PREFIX has a marker for other pins than $("$PY" -I "$HERE/fragpipe_quant.py" pins | tr '\n' ' '); move it away or give another --prefix."
  exit 2
fi
if [ -e "$PREFIX" ]; then
  say "REFUSED: $PREFIX exists without $(basename "$MARKER") (a build that did not finish, or another folder). Nothing was removed: move it away or give another --prefix."
  exit 2
fi

say "base python: $PY ($("$PY" -c 'import platform; print(platform.python_version())'))"
say "prefix: $PREFIX"
PARENT="$(dirname "$PREFIX")"
if [ ! -d "$PARENT" ]; then
  mkdir -p "$PARENT"
  chmod g+rx,o-w "$PARENT"
fi
"$PY" -m venv "$PREFIX"
export PIP_CACHE_DIR="$PREFIX/.pip-cache" PYTHONNOUSERSITE=1 PIP_DISABLE_PIP_VERSION_CHECK=1
PINS=()
while IFS= read -r line; do [ -n "$line" ] && PINS+=("$line"); done \
  < <("$PY" -I "$HERE/fragpipe_quant.py" pins)
say "pip install ${PINS[*]}"
"$PREFIX/bin/python" -m pip install --no-cache-dir "${PINS[@]}"
rm -rf "$PREFIX/.pip-cache"
"$PREFIX/bin/python" -m pip freeze > "$PREFIX/pip_freeze.txt"
# the check runs with a clean environment: the venv must work without the module loaded
if ! env -i HOME="$HOME" PATH=/usr/bin:/bin "$PREFIX/bin/python" -I "$HERE/fragpipe_quant.py" check-env > "$PREFIX/check_env.json"; then
  cat "$PREFIX/check_env.json" >&2 || true
  say "FAILED: the new venv does not pass directLFQ's guard checks (above); no marker written, so it is never used."
  exit 3
fi
env -i HOME="$HOME" PATH=/usr/bin:/bin "$PREFIX/bin/python" -I "$HERE/fragpipe_quant.py" write-marker \
  --prefix "$PREFIX" --base-python "$PY" --freeze "$PREFIX/pip_freeze.txt" >/dev/null
chmod -R u+rwX,g+rX,g-w,o-w "$PREFIX"
say "built: $PREFIX"
cat "$PREFIX/check_env.json"
