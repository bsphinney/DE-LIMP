#!/usr/bin/env bash
# ensure_dotnet8.sh -- make a .NET 8 install available that carries BOTH shared runtimes the
# pipeline's Thermo .raw readers need:
#   Microsoft.NETCore.App    8.0.>=17  the DIA-NN 2.6 NATIVE binary reads .raw directly with it
#   Microsoft.AspNetCore.App 8.0.x     a framework-dependent ThermoRawFileParser 2.x needs it
#
# Why NETCore >= 8.0.17: DIA-NN 2.6's .raw reader is a bundled .NET (ThermoFisher.CommonCore)
# component. Without a suitable runtime it aborts with:
#   "ERROR: cannot read .raw files, please download and install .NET Runtime 8:
#    8.0.17 or later"
# and processes 0 files. DIA-NN does a HARD check for >= 8.0.17 -- an older 8.0.x
# (e.g. a cluster's 8.0.4 module) is REJECTED even though the app's runtimeconfig
# only asks for 8.0.0, and a .NET 9 runtime is NOT used (rollForward is LatestMinor,
# it won't cross the major version). So we need an actual 8.0.>=17 runtime.
# When it works, DIA-NN logs: ".NET runtime found, Thermo .raw support enabled".
#
# Why AspNetCore: ThermoRawFileParser 2.0.0.0's framework-dependent build (the Core's shared
# copy on HIVE is one) lists BOTH Microsoft.NETCore.App 8.0.0 and Microsoft.AspNetCore.App
# 8.0.0 in its ThermoRawFileParser.runtimeconfig.json. With this script's old NETCore-only
# install as DOTNET_ROOT it exits 150, "You must install or update .NET to run this
# application." (a Core member, 2026-09-23, 15 Fusion Lumos .raw: every file came back unknown). A
# root with both is ONE root that serves both readers -- HIVE's dotnet-core-sdk/8.0.4 module
# has both runtimes but is 8.0.4, which DIA-NN rejects.
#
# Idempotent: reuses a root that already has both; into its own DEST it installs only the
# runtime that is missing (an existing NETCore-only DEST gets AspNetCore added in place).
# It needs internet only when it has something to install. On a cluster login node that part runs
# on a compute node (on_compute.py, one blocking srun): HPC@UCD's ucd-hpc rules keep downloads and
# installs off login nodes, and HIVE's compute nodes reach Microsoft's downloads (a compute-node
# job, 2026-10-16). It installs to a shared path so compute nodes can read it at run time.
#
# If AspNetCore cannot be installed but a NETCore >= 8.0.17 root exists, that root is still
# printed (exit 0) with a WARNING: DIA-NN's callers (run_search.py, diann_parallel.py) read
# the last line and must not lose a .raw search that worked before AspNetCore was asked for.
# detect_acquisition.py checks the frameworks itself and says what is missing.
#
# Prints DOTNET_ROOT on the LAST stdout line (diagnostics go to stderr):
#   eval "$(bash ensure_dotnet8.sh)"            # not this -- it prints a path, not exports
#   export DOTNET_ROOT="$(bash ensure_dotnet8.sh | tail -1)"; export PATH="$DOTNET_ROOT:$PATH"
#
# Override the install location with PROTEOMICS_DOTNET_DIR. Default: dotnet8 in the skill's
# per-user folder, $PROTEOMICS_PIPELINE_HOME (setup.sh's PP_HOME; activate.sh exports it), else
# ~/.proteomics-pipeline/dotnet8 -- so an install into another prefix never writes to $HOME
# (2.15.1, HIVE F5). A good root at the old default (~/.proteomics-pipeline/dotnet8) is reused.
#
# --print-install-tmp: print the folder Microsoft's installer would be given as TMPDIR, and exit
# (installer_tmp_parent below; the tests use it).
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "$SCRIPT_DIR/verify_download.sh"
MIN_MINOR=17
LEGACY_DEST="$HOME/.proteomics-pipeline/dotnet8"
DEST="${PROTEOMICS_DOTNET_DIR:-${PROTEOMICS_PIPELINE_HOME:-$HOME/.proteomics-pipeline}/dotnet8}"

# Microsoft's dotnet-install.sh unpacks into "${TMPDIR:-/tmp}/dotnet.XXXXXXXXX" and sorts the files
# it copies by the regex '^.*/[0-9]+\.[0-9]+[^/]+/' over their FULL path: a TMPDIR with a
# version-like folder in it (".../2.15.1/...") matches for every file, the copy fails ("cp: cannot
# stat"), and no runtime is installed (2.15.1, HIVE F6). So the installer gets a TMPDIR of its own,
# in the first of these folders whose path has no such component.
VERSION_LIKE_DIR='/[0-9]+\.[0-9]+[^/]+/'
installer_tmp_parent() {
  local d
  for d in "${TMPDIR:-}" "$(dirname "$DEST")" /tmp /var/tmp; do
    [ -n "$d" ] && [ -d "$d" ] && [ -w "$d" ] || continue
    case "$d" in /*) ;; *) d="$(cd "$d" && pwd)" ;; esac
    if ! printf '%s/' "${d%/}" | grep -Eq "$VERSION_LIKE_DIR"; then
      printf '%s\n' "${d%/}"; return 0
    fi
  done
  return 1
}
if [ "${1:-}" = "--print-install-tmp" ]; then installer_tmp_parent; exit $?; fi

runtimes() { [ -x "$1/dotnet" ] && "$1/dotnet" --list-runtimes 2>/dev/null; }
# true if $1/dotnet exposes a Microsoft.NETCore.App 8.0.<patch> with patch >= MIN_MINOR
have_netcore() {
  runtimes "$1" | awk -v m="$MIN_MINOR" '
    $1 == "Microsoft.NETCore.App" && $2 ~ /^8\.0\./ { split($2, a, "."); if (a[3]+0 >= m) ok=1 }
    END { exit ok?0:1 }'
}
# true if $1/dotnet exposes any Microsoft.AspNetCore.App 8.0.x (TRFP asks for 8.0.0 and rolls
# forward on patch)
have_aspnet() {
  runtimes "$1" | awk '$1 == "Microsoft.AspNetCore.App" && $2 ~ /^8\.0\./ { ok=1 }
    END { exit ok?0:1 }'
}
have_ok() { have_netcore "$1" && have_aspnet "$1"; }

sysroot=""
if command -v dotnet >/dev/null 2>&1; then
  sysroot="$(dirname "$(readlink -f "$(command -v dotnet)")")"
fi

# 1) already installed where we put it (or, with no PROTEOMICS_DOTNET_DIR, at the old default)
if have_ok "$DEST"; then echo "reusing .NET 8 at $DEST" >&2; echo "$DEST"; exit 0; fi
if [ -z "${PROTEOMICS_DOTNET_DIR:-}" ] && [ "$LEGACY_DEST" != "$DEST" ] && have_ok "$LEGACY_DEST"; then
  echo "reusing .NET 8 at $LEGACY_DEST" >&2; echo "$LEGACY_DEST"; exit 0
fi

# 2) a system/module dotnet that already has both
if [ -n "$sysroot" ] && have_ok "$sysroot"; then
  echo "system .NET 8 ok at $sysroot" >&2; echo "$sysroot"; exit 0
fi

# 3) install whichever runtime DEST lacks (needs internet). Steps 1-2 above are checks and stay
# where they are; a download on a cluster login node (on_compute.py --is-login-node) is moved to
# a compute node first: this script again, under on_compute.py. srun passes stdout through, so
# the root is still the LAST line. With no python3 + on_compute.py on a cluster, nothing is
# downloaded here, and step 4 still prints a root DIA-NN can use when there is one.
SKIP_INSTALL=false
here="$(cd "$(dirname "$0")" && pwd)"
if [ -z "${SLURM_JOB_ID:-}" ] && [ -z "${PROTEOMICS_ON_COMPUTE:-}" ]; then
  if command -v python3 >/dev/null 2>&1 && [ -f "$here/on_compute.py" ]; then
    if python3 "$here/on_compute.py" --is-login-node; then
      exec python3 "$here/on_compute.py" --time 30 --mem 4G --name pp_dotnet8 \
        --why "ensure_dotnet8.sh downloads and installs .NET 8" -- "${BASH:-bash}" "$here/ensure_dotnet8.sh" "$@"
    fi
  elif command -v sbatch >/dev/null 2>&1 \
       || [ -x "${HPCCF_SLURM_BIN:-/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin}/sbatch" ] \
       || [ -r "${UCD_HPC_CLUSTER_CONF:-/opt/hpccf/etc/cluster_name.conf}" ]; then
    echo "WARNING: .NET 8 needs installing, and this is a cluster login node with no python3 +" \
         "on_compute.py to move the download to a compute node: nothing is downloaded here. Run it" \
         "in a job: srun --account=<yours> --partition=<yours> --time=30 --mem=4G bash" \
         "$here/ensure_dotnet8.sh" >&2
    SKIP_INSTALL=true
  fi
fi
# install_runtime <dotnet|aspnetcore>; returns non-zero instead of exiting, so a failed
# AspNetCore download can still fall back to a NETCore root below
install_runtime() {
  local tmp rc=0 parent work
  echo "installing the .NET 8 '$1' runtime into $DEST ..." >&2
  mkdir -p "$DEST" || return 1
  if ! parent="$(installer_tmp_parent)"; then
    echo "no writable temporary folder without a version-like name (such as 2.15.1) in its path" \
         "(tried \$TMPDIR, $(dirname "$DEST"), /tmp, /var/tmp): Microsoft's installer cannot" \
         "copy from such a folder. Set TMPDIR to one, then run this again." >&2
    return 1
  fi
  work="$(mktemp -d "$parent/ppdotnet.XXXXXX")" || return 1
  tmp="$work/dotnet-install.sh"
  # Bounded: every existing NETCore-only install now reaches this line once, and on a host
  # with no internet (a compute node, another site's login node) curl's default connect
  # timeout would stall each DIA-NN .raw generation for minutes before the fallback below
  # (review 2026-09-23). No internet -> fail in 15 s and degrade.
  # Microsoft's installer is fetched from a FIXED COMMIT of dotnet/install-scripts (the same file
  # https://dot.net/v1/dotnet-install.sh served on 2026-11-17) and checked against its sha256 BEFORE
  # bash runs it -- it used to be piped from a URL whose content changes whenever Microsoft likes.
  # The pin is in pinned_checksums.txt (key dotnet-install-5147e32); to move to a newer script, change
  # the commit below and that line together.
  local commit=5147e32300a8e908f5d737c8cff63a76b4b63531
  if curl -sSL --connect-timeout 15 --max-time 120 \
       "https://raw.githubusercontent.com/dotnet/install-scripts/$commit/src/dotnet-install.sh" -o "$tmp"; then
    if verify_download "dotnet-install-${commit:0:7}" "$tmp" "dotnet-install.sh"; then
      TMPDIR="$work" bash "$tmp" --channel 8.0 --runtime "$1" --install-dir "$DEST" >&2 || rc=$?
    else
      rc=1; echo "the .NET installer script does not match its pinned sha256 (or no tool can check it), so it was NOT run" >&2
    fi
  else
    rc=$?; echo "could not download the .NET installer script (dotnet/install-scripts@${commit:0:7})" >&2
  fi
  rm -rf "$work"
  return "$rc"
}
if ! $SKIP_INSTALL; then
  have_netcore "$DEST" || install_runtime dotnet || true
  have_aspnet "$DEST" || install_runtime aspnetcore || true
fi
if have_ok "$DEST"; then echo "installed .NET 8 at $DEST" >&2; echo "$DEST"; exit 0; fi

# 4) degraded: DIA-NN can still read .raw; a framework-dependent ThermoRawFileParser cannot
for root in "$DEST" "$sysroot"; do
  if [ -n "$root" ] && have_netcore "$root"; then
    echo "WARNING: $root has Microsoft.NETCore.App 8.0.>=$MIN_MINOR (DIA-NN can read .raw)" \
         "but no Microsoft.AspNetCore.App 8.0, and installing it into $DEST failed. A" \
         "framework-dependent ThermoRawFileParser will exit 150 'You must install or update" \
         ".NET'. Re-run this script where there is internet, or use a self-contained parser" \
         "(bioconda thermorawfileparser, which setup.sh installs)." >&2
    echo "$root"; exit 0
  fi
done
echo "ERROR: .NET 8 install did not yield a >= 8.0.$MIN_MINOR runtime in $DEST" >&2
exit 1
