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
# application." (gabrig 2026-09-23, 15 Fusion Lumos .raw: every file came back unknown). A
# root with both is ONE root that serves both readers -- HIVE's dotnet-core-sdk/8.0.4 module
# has both runtimes but is 8.0.4, which DIA-NN rejects.
#
# Idempotent: reuses a root that already has both; into its own DEST it installs only the
# runtime that is missing (an existing NETCore-only DEST gets AspNetCore added in place).
# Run it on a machine with internet (an HPC LOGIN node) -- compute nodes often have none.
# It installs to a shared path so compute nodes can read it at run time.
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
# Override the install location with PROTEOMICS_DOTNET_DIR (default ~/.proteomics-pipeline/dotnet8).
set -euo pipefail
MIN_MINOR=17
DEST="${PROTEOMICS_DOTNET_DIR:-$HOME/.proteomics-pipeline/dotnet8}"

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

# 1) already installed where we put it
if have_ok "$DEST"; then echo "reusing .NET 8 at $DEST" >&2; echo "$DEST"; exit 0; fi

# 2) a system/module dotnet that already has both
if [ -n "$sysroot" ] && have_ok "$sysroot"; then
  echo "system .NET 8 ok at $sysroot" >&2; echo "$sysroot"; exit 0
fi

# 3) install whichever runtime DEST lacks (needs internet -- run on a login node)
# install_runtime <dotnet|aspnetcore>; returns non-zero instead of exiting, so a failed
# AspNetCore download can still fall back to a NETCore root below
install_runtime() {
  local tmp rc=0
  echo "installing the .NET 8 '$1' runtime into $DEST ..." >&2
  mkdir -p "$DEST" || return 1
  tmp="$(mktemp)" || return 1
  # Bounded: every existing NETCore-only install now reaches this line once, and on a host
  # with no internet (a compute node, another site's login node) curl's default connect
  # timeout would stall each DIA-NN .raw generation for minutes before the fallback below
  # (review 2026-09-23). No internet -> fail in 15 s and degrade.
  if curl -sSL --connect-timeout 15 --max-time 120 https://dot.net/v1/dotnet-install.sh -o "$tmp"; then
    bash "$tmp" --channel 8.0 --runtime "$1" --install-dir "$DEST" >&2 || rc=$?
  else
    rc=$?; echo "could not download https://dot.net/v1/dotnet-install.sh" >&2
  fi
  rm -f "$tmp"
  return "$rc"
}
have_netcore "$DEST" || install_runtime dotnet || true
have_aspnet "$DEST" || install_runtime aspnetcore || true
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
