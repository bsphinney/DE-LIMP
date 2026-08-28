#!/usr/bin/env bash
# ensure_dotnet8.sh -- make a .NET 8 runtime (>= 8.0.17) available so the DIA-NN 2.6
# NATIVE binary can read Thermo .raw files directly (no mzML conversion).
#
# Why this exists: DIA-NN 2.6's .raw reader is a bundled .NET (ThermoFisher.CommonCore)
# component. Without a suitable runtime it aborts with:
#   "ERROR: cannot read .raw files, please download and install .NET Runtime 8:
#    8.0.17 or later"
# and processes 0 files. DIA-NN does a HARD check for >= 8.0.17 -- an older 8.0.x
# (e.g. a cluster's 8.0.4 module) is REJECTED even though the app's runtimeconfig
# only asks for 8.0.0, and a .NET 9 runtime is NOT used (rollForward is LatestMinor,
# it won't cross the major version). So we need an actual 8.0.>=17 runtime.
#
# When it works, DIA-NN logs: ".NET runtime found, Thermo .raw support enabled".
#
# Idempotent: reuses an existing good runtime; only downloads if none is found.
# Run it on a machine with internet (an HPC LOGIN node) -- compute nodes often
# have none. It installs to a shared path so compute nodes can read it at run time.
#
# Prints DOTNET_ROOT on the LAST stdout line (diagnostics go to stderr):
#   eval "$(bash ensure_dotnet8.sh)"            # not this -- it prints a path, not exports
#   export DOTNET_ROOT="$(bash ensure_dotnet8.sh | tail -1)"; export PATH="$DOTNET_ROOT:$PATH"
#
# Override the install location with PROTEOMICS_DOTNET_DIR (default ~/.proteomics-pipeline/dotnet8).
set -euo pipefail
MIN_MINOR=17
DEST="${PROTEOMICS_DOTNET_DIR:-$HOME/.proteomics-pipeline/dotnet8}"

# true if $1/dotnet exposes a Microsoft.NETCore.App 8.0.<minor>  with minor >= MIN_MINOR
have_ok() {
  local root="$1"
  [ -x "$root/dotnet" ] || return 1
  "$root/dotnet" --list-runtimes 2>/dev/null | awk -v m="$MIN_MINOR" '
    /Microsoft.NETCore.App 8\.0\./ { split($2, a, "."); if (a[3]+0 >= m) ok=1 }
    END { exit ok?0:1 }'
}

# 1) already installed where we put it
if have_ok "$DEST"; then echo "reusing .NET 8 at $DEST" >&2; echo "$DEST"; exit 0; fi

# 2) a system/module dotnet that already satisfies >= 8.0.17
if command -v dotnet >/dev/null 2>&1; then
  sysroot="$(dirname "$(readlink -f "$(command -v dotnet)")")"
  if have_ok "$sysroot"; then echo "system .NET 8 ok at $sysroot" >&2; echo "$sysroot"; exit 0; fi
fi

# 3) install the latest 8.0 runtime into DEST (needs internet -- run on a login node)
echo "installing .NET 8 runtime (>= 8.0.$MIN_MINOR) into $DEST ..." >&2
mkdir -p "$DEST"
tmp="$(mktemp)"
curl -sSL https://dot.net/v1/dotnet-install.sh -o "$tmp"
bash "$tmp" --channel 8.0 --runtime dotnet --install-dir "$DEST" >&2
rm -f "$tmp"
have_ok "$DEST" || { echo "ERROR: .NET 8 install did not yield a >= 8.0.$MIN_MINOR runtime" >&2; exit 1; }
echo "installed .NET 8 at $DEST" >&2
echo "$DEST"
