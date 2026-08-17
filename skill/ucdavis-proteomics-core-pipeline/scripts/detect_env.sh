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

# ------------------------------------------------------ engine capability ----
# Which search engines can actually RUN on this machine, and how. Every claim
# below was checked against the projects' released artifacts on 2026-08-14 --
# none of it is assumed:
#
#   DIA-NN   ZERO macOS assets in ANY release (only *-Academia-Linux.zip and
#            *-Academia.msi). On a Mac it runs only inside a linux/amd64
#            container; on Apple Silicon that means Rosetta emulation. Works,
#            but expect a large slowdown.
#   Sage     Ships native aarch64-apple-darwin AND x86_64-apple-darwin, so it
#            runs natively on Apple Silicon. DDA only.
#   FragPipe 24.0 ships ONLY FragPipe-24.0-installer.exe and -linux.zip. There
#            is NO macOS build -- Intel or ARM. Same for diaTracer.
#   Radiant  Container only, but the image is multi-arch (linux/amd64 +
#            linux/arm64), so it runs NATIVELY on Apple Silicon -- the only DIA
#            engine that does.
APPLE_SILICON=false
[ "$OS" = "darwin" ] && [ "$ARCH" = "arm64" ] && APPLE_SILICON=true
HAS_JAVA=false; have java && HAS_JAVA=true

# NOTE: key these on $OS, not $CLASS. CLASS is decided by SLURM/quobyte BEFORE
# the darwin check, so a Mac with the SLURM client installed comes back as
# "hpc" -- and "there is no macOS build" is a property of the OS, not of
# whether a cluster scheduler happens to be on PATH.
case "$OS" in
  linux) DIANN_OK=true;  DIANN_HOW=native; DIANN_NOTE="native Linux build" ;;
  darwin)
    if $HAS_DOCKER; then
      DIANN_OK=true; DIANN_HOW=docker
      if $APPLE_SILICON; then
        DIANN_NOTE="no macOS build exists; runs as a linux/amd64 container under Rosetta emulation -- correct results, but far slower than native. For a big cohort, run on HIVE instead."
      else
        DIANN_NOTE="no macOS build exists; runs in a linux/amd64 container"
      fi
    else
      DIANN_OK=false; DIANN_HOW=null
      DIANN_NOTE="no macOS build exists and Docker is not installed -- install Docker Desktop (build_diann_docker.sh then builds the image), or run on HIVE"
    fi ;;
esac

SAGE_OK=true; SAGE_HOW=native; SAGE_NOTE="native build for this OS/arch (DDA only)"
$APPLE_SILICON && SAGE_NOTE="native aarch64-apple-darwin build -- runs natively on Apple Silicon (DDA only)"

if [ "$OS" = "darwin" ]; then
  FRAGPIPE_OK=false; FRAGPIPE_HOW=null
  FRAGPIPE_NOTE="FragPipe ships only a Windows installer and a Linux zip -- there is NO macOS build, Intel or Apple Silicon. This route (and diaTracer) is unavailable on any Mac; run it on HIVE."
elif $HAS_JAVA; then
  FRAGPIPE_OK=true; FRAGPIPE_HOW=native
  FRAGPIPE_NOTE="Linux build; MSFragger/IonQuant/diaTracer are licence-gated -- set FRAGPIPE_TOOLS_FOLDER"
else
  FRAGPIPE_OK=false; FRAGPIPE_HOW=null
  FRAGPIPE_NOTE="FragPipe needs Java and none is on PATH"
fi

if $HAS_DOCKER || $HAS_APPTAINER; then
  RADIANT_OK=true
  RADIANT_HOW=docker; $HAS_APPTAINER && [ "$CLASS" = "hpc" ] && RADIANT_HOW=apptainer
  if $APPLE_SILICON; then
    RADIANT_NOTE="multi-arch image runs NATIVELY on Apple Silicon -- the only DIA engine that does. Licence-restricted (Commons Clause) for fee-for-service work."
  else
    RADIANT_NOTE="container only (~3 GB first pull). Licence-restricted (Commons Clause) for fee-for-service work."
  fi
else
  RADIANT_OK=false; RADIANT_HOW=null
  RADIANT_NOTE="Seer ships no native binary -- needs Docker or Apptainer"
fi

j() { [ "$1" = "null" ] && printf null || printf '"%s"' "$1"; }

cat <<JSON
{
  "os": "$OS",
  "arch": "$ARCH",
  "apple_silicon": $APPLE_SILICON,
  "platform_class": "$CLASS",
  "container_runtime": "$RUNTIME",
  "engines": {
    "diann":    {"available": $DIANN_OK,    "how": $(j "$DIANN_HOW"),    "note": "$DIANN_NOTE"},
    "sage":     {"available": $SAGE_OK,     "how": $(j "$SAGE_HOW"),     "note": "$SAGE_NOTE"},
    "fragpipe": {"available": $FRAGPIPE_OK, "how": $(j "$FRAGPIPE_HOW"), "note": "$FRAGPIPE_NOTE"},
    "radiant":  {"available": $RADIANT_OK,  "how": $(j "$RADIANT_HOW"),  "note": "$RADIANT_NOTE"}
  },
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
