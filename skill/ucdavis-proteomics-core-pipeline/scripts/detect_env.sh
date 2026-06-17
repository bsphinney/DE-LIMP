#!/usr/bin/env bash
# =============================================================================
# detect_env.sh  --  Detect the runtime environment and emit a capability map
# as JSON on stdout. The skill reads this to decide how to acquire and run the
# search engines (native binary vs Docker vs Apptainer vs SLURM submission).
#
# platform_class:
#   hpc    a SLURM cluster (sbatch present, or /quobyte visible -> UC Davis HIVE)
#   mac    macOS (DIA-NN has no native build here -> Docker required)
#   linux  a Linux workstation (native binaries OK)
#
# Usage:  bash detect_env.sh  [> env.json]
# =============================================================================
set -euo pipefail

have() { command -v "$1" >/dev/null 2>&1; }
ver()  { "$@" 2>&1 | head -n1 | tr -d '"' | sed 's/[[:cntrl:]]//g' || true; }

OS="$(uname -s | tr '[:upper:]' '[:lower:]')"   # linux | darwin
ARCH="$(uname -m)"

HAS_SLURM=false;     have sbatch     && HAS_SLURM=true
HAS_MODULE=false;    type module >/dev/null 2>&1 && HAS_MODULE=true
HAS_DOCKER=false;    have docker     && HAS_DOCKER=true
HAS_APPTAINER=false; ( have apptainer || have singularity ) && HAS_APPTAINER=true
HAS_CONDA=false;     ( have conda || have mamba || have micromamba ) && HAS_CONDA=true
QUOBYTE=false;       [ -d /quobyte/proteomics-grp ] && QUOBYTE=true

JAVA_VER=null
if have java; then JAVA_VER="\"$(ver java -version)\""; fi
R_VER=null
if have Rscript; then R_VER="\"$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo unknown)\""; fi

# Classify
if $HAS_SLURM || $QUOBYTE; then
  CLASS="hpc"
elif [ "$OS" = "darwin" ]; then
  CLASS="mac"
else
  CLASS="linux"
fi

# Container runtime preference for this class
RUNTIME="native"
if [ "$CLASS" = "hpc" ] && $HAS_APPTAINER; then RUNTIME="apptainer"
elif [ "$CLASS" = "mac" ] && $HAS_DOCKER;   then RUNTIME="docker"
elif [ "$CLASS" = "linux" ];                then RUNTIME="native"
fi

cat <<JSON
{
  "os": "$OS",
  "arch": "$ARCH",
  "platform_class": "$CLASS",
  "container_runtime": "$RUNTIME",
  "has_slurm": $HAS_SLURM,
  "has_module": $HAS_MODULE,
  "has_docker": $HAS_DOCKER,
  "has_apptainer": $HAS_APPTAINER,
  "has_conda": $HAS_CONDA,
  "has_java": $JAVA_VER,
  "has_R": $R_VER,
  "uc_davis_hive": $QUOBYTE
}
JSON
