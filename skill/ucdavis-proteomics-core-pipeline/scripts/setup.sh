#!/usr/bin/env bash
# =============================================================================
# setup.sh  --  One-shot environment bootstrap for the proteomics-pipeline skill.
#
# GOAL: an average biologist drops this skill into an agentic AI and says
# "analyze my proteomics data". The AI runs THIS first. Everything that can be
# installed WITHOUT admin rights is installed automatically into one
# self-contained conda environment; anything that genuinely needs the user
# (Docker Desktop for DIA-NN on macOS) is reported with exact next steps.
#
# What gets installed (no sudo, into ~/.proteomics-pipeline/):
#   - micromamba          (single static binary; only if no conda/mamba present)
#   - a conda env with:   python + pyarrow + pyyaml
#                         R (>=4.5) + bioconductor-limpa + bioconductor-limma
#                                   + r-arrow + r-dplyr + r-tidyr
#                         sage-proteomics            (the DDA search engine)
#                         proteowizard/msconvert     (LINUX ONLY on bioconda)
#
# What stays special-cased (handled elsewhere / reported):
#   - DIA-NN: license-gated, no conda. Linux -> binary (acquire_tools.sh);
#             HIVE -> existing .sif; macOS -> Docker (see notes below).
#
# Outputs:
#   ~/.proteomics-pipeline/activate.sh   <- source this; puts the env on PATH
#   ~/.proteomics-pipeline/setup.json    <- machine-readable readiness report
#
# Usage:  bash setup.sh            # install/repair everything it can
#         bash setup.sh --check    # report only, install nothing
# =============================================================================
set -uo pipefail

PP_HOME="${PP_HOME:-$HOME/.proteomics-pipeline}"
ENV_NAME="proteomics-pipeline"
MAMBA_ROOT="$PP_HOME/micromamba"
ENV_PREFIX="$MAMBA_ROOT/envs/$ENV_NAME"
ACTIVATE="$PP_HOME/activate.sh"
SETUP_JSON="$PP_HOME/setup.json"
CHECK_ONLY=false; [ "${1:-}" = "--check" ] && CHECK_ONLY=true
mkdir -p "$PP_HOME"

OS="$(uname -s | tr '[:upper:]' '[:lower:]')"   # darwin | linux
ARCH="$(uname -m)"                               # arm64 | x86_64 | aarch64
have() { command -v "$1" >/dev/null 2>&1; }
say()  { printf '%s\n' "$*" >&2; }
NOTES=()

# micromamba platform slug
case "$OS-$ARCH" in
  darwin-arm64)        MM_PLAT="osx-arm64" ;;
  darwin-x86_64)       MM_PLAT="osx-64" ;;
  linux-x86_64)        MM_PLAT="linux-64" ;;
  linux-aarch64|linux-arm64) MM_PLAT="linux-aarch64" ;;
  *)                   MM_PLAT="linux-64" ;;
esac

# ---- 1. find (or install) a conda-family package manager --------------------
CONDA=""
pick_conda() {
  if   have micromamba;                 then CONDA="micromamba"
  elif [ -x "$MAMBA_ROOT/bin/micromamba" ]; then CONDA="$MAMBA_ROOT/bin/micromamba"
  elif have mamba;                      then CONDA="mamba"
  elif have conda;                      then CONDA="conda"
  fi
}
install_micromamba() {
  say "[setup] installing micromamba (no admin needed) for $MM_PLAT ..."
  mkdir -p "$MAMBA_ROOT/bin"
  # official endpoint streams a tarball containing bin/micromamba
  if curl -Ls "https://micro.mamba.pm/api/micromamba/$MM_PLAT/latest" \
       | tar -xj -C "$MAMBA_ROOT" bin/micromamba 2>/dev/null; then
    CONDA="$MAMBA_ROOT/bin/micromamba"
    say "[setup] micromamba installed at $CONDA"
  else
    NOTES+=("Could not download micromamba automatically. Install it manually: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html  then re-run setup.sh.")
  fi
}

pick_conda
if [ -z "$CONDA" ]; then
  if $CHECK_ONLY; then NOTES+=("No conda/micromamba found; run setup.sh (without --check) to install it.")
  else install_micromamba; fi
fi

# create-env helper that works for micromamba OR conda/mamba
create_env() {
  local pkgs=(python=3.11 pyarrow pyyaml
              "r-base>=4.5" bioconductor-limpa bioconductor-limma
              r-arrow r-dplyr r-tidyr sage-proteomics
              pandoc python-docx)   # pandoc + python-docx: Markdown report -> Word .docx
  # proteowizard (msconvert) is bioconda LINUX-only
  if [ "$OS" = "linux" ]; then pkgs+=(proteowizard); fi

  say "[setup] solving + installing the analysis environment (this can take a few minutes)..."
  case "$CONDA" in
    *micromamba)
      "$CONDA" create -y -r "$MAMBA_ROOT" -n "$ENV_NAME" \
        -c conda-forge -c bioconda "${pkgs[@]}" ;;
    *)
      "$CONDA" create -y -p "$ENV_PREFIX" \
        -c conda-forge -c bioconda "${pkgs[@]}" ;;
  esac
}

# ---- 2. ensure the environment exists ---------------------------------------
env_ready() { [ -x "$ENV_PREFIX/bin/python" ] && [ -x "$ENV_PREFIX/bin/Rscript" ]; }

if [ -n "$CONDA" ] && ! env_ready; then
  if $CHECK_ONLY; then NOTES+=("Analysis env not built yet; run setup.sh to create it.")
  else
    create_env || NOTES+=("Environment solve failed. Try: $CONDA create -r $MAMBA_ROOT -n $ENV_NAME -c conda-forge -c bioconda bioconductor-limpa sage-proteomics python pyarrow")
  fi
fi

# limpa is on bioconda, but if the solve dropped it, install via BiocManager.
if env_ready && ! "$ENV_PREFIX/bin/Rscript" -e 'q(status=!requireNamespace("limpa",quietly=TRUE))' 2>/dev/null; then
  if ! $CHECK_ONLY; then
    say "[setup] limpa missing from env; installing via BiocManager..."
    "$ENV_PREFIX/bin/Rscript" -e 'if(!requireNamespace("BiocManager",quietly=TRUE))install.packages("BiocManager",repos="https://cloud.r-project.org");BiocManager::install("limpa",update=FALSE,ask=FALSE)' \
      || NOTES+=("limpa could not be installed. DE --method dpc will be unavailable; --method maxlfq still works (limma only).")
  fi
fi

# ---- 3. resolve tool paths --------------------------------------------------
resolve() { [ -x "$ENV_PREFIX/bin/$1" ] && echo "$ENV_PREFIX/bin/$1" || (have "$1" && command -v "$1" || echo ""); }
PY="$(resolve python)";     [ -z "$PY" ] && PY="$(command -v python3 || true)"
RSCRIPT="$(resolve Rscript)"
SAGE="$(resolve sage)"
MSCONVERT="$(resolve msconvert)"
HAS_DOCKER=false; have docker && HAS_DOCKER=true
HAS_APPTAINER=false; ( have apptainer || have singularity ) && HAS_APPTAINER=true
QUOBYTE=false; [ -d /quobyte/proteomics-grp ] && QUOBYTE=true

# DIA-NN reachability by platform
DIANN_PATH="diann_engine"; DIANN_READY=false; DIANN_NOTE=""
if   $QUOBYTE && $HAS_APPTAINER; then DIANN_READY=true;  DIANN_NOTE="HIVE: reuse existing .sif (acquire_tools.sh resolves it)."
elif [ "$OS" = "linux" ];        then DIANN_READY=true;  DIANN_NOTE="Linux: acquire_tools.sh downloads the free DIA-NN Academia binary."
elif [ "$OS" = "darwin" ] && $HAS_DOCKER; then DIANN_READY=true; DIANN_NOTE="macOS+Docker: build the image with build_diann_docker.sh, then export DIANN_DOCKER_IMAGE."
elif [ "$OS" = "darwin" ];       then DIANN_READY=false; DIANN_NOTE="macOS: DIA-NN has NO native build. Install Docker Desktop (https://docs.docker.com/desktop/setup/install/mac-install/), then re-run setup.sh and build_diann_docker.sh."
fi

# readiness per acquisition type (for the orchestrator to gate on)
DE_READY=false; [ -n "$RSCRIPT" ] && DE_READY=true
DDA_READY=false
if [ -n "$SAGE" ] && [ -n "$RSCRIPT" ]; then DDA_READY=true; fi
DIA_READY=false
if $DIANN_READY && [ -n "$RSCRIPT" ]; then DIA_READY=true; fi

[ -z "$RSCRIPT" ] && NOTES+=("R/Rscript not available — DE cannot run. Re-run setup.sh to install it into the conda env.")
[ -z "$SAGE" ]    && NOTES+=("Sage not found — DDA search unavailable until the conda env is built.")
[ "$OS" = "darwin" ] && [ -z "$MSCONVERT" ] && NOTES+=("msconvert is Linux-only on bioconda. On macOS, Sage can only search files ALREADY in mzML; convert Bruker .d / Thermo .raw elsewhere first, or use DIA-NN (which reads .d/.raw natively) for DIA data.")

# ---- 4. write activate.sh + setup.json --------------------------------------
if ! $CHECK_ONLY || [ ! -f "$ACTIVATE" ]; then
  cat > "$ACTIVATE" <<EOF
# source this to put the proteomics-pipeline environment on PATH
export PROTEOMICS_PIPELINE_HOME="$PP_HOME"
export PATH="$ENV_PREFIX/bin:\$PATH"
[ -f "$PP_HOME/diann_docker_image" ] && export DIANN_DOCKER_IMAGE="\$(cat "$PP_HOME/diann_docker_image")"
EOF
fi

j() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }
{
  printf '{\n'
  printf '  "os": "%s", "arch": "%s",\n' "$OS" "$ARCH"
  printf '  "conda": "%s",\n' "$(j "${CONDA:-}")"
  printf '  "env_prefix": "%s",\n' "$(j "$ENV_PREFIX")"
  printf '  "activate": "%s",\n' "$(j "$ACTIVATE")"
  printf '  "python": "%s",\n'   "$(j "$PY")"
  printf '  "rscript": "%s",\n'  "$(j "$RSCRIPT")"
  printf '  "sage": "%s",\n'     "$(j "$SAGE")"
  printf '  "msconvert": "%s",\n' "$(j "$MSCONVERT")"
  printf '  "has_docker": %s, "has_apptainer": %s, "uc_davis_hive": %s,\n' \
         "$($HAS_DOCKER && echo true || echo false)" \
         "$($HAS_APPTAINER && echo true || echo false)" \
         "$($QUOBYTE && echo true || echo false)"
  printf '  "diann": {"ready": %s, "note": "%s"},\n' "$($DIANN_READY && echo true || echo false)" "$(j "$DIANN_NOTE")"
  printf '  "ready_for": {"de": %s, "dia": %s, "dda": %s},\n' \
         "$($DE_READY && echo true || echo false)" \
         "$($DIA_READY && echo true || echo false)" \
         "$($DDA_READY && echo true || echo false)"
  printf '  "notes": ['
  for i in "${!NOTES[@]}"; do
    printf '%s"%s"' "$( [ "$i" -gt 0 ] && echo ', ' )" "$(j "${NOTES[$i]}")"
  done
  printf ']\n}\n'
} | tee "$SETUP_JSON"

say ""
say "[setup] activate with:  source $ACTIVATE"
say "[setup] readiness report written to $SETUP_JSON"
