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
#   - thermorawfileparser (bioconda, self-contained: reads Thermo .raw for
#                         detect_acquisition.py) + pythonnet (conda-forge: lets
#                         thermo_resolution.py read the Orbitrap resolution) + pandas
#                         -- one separate, non-fatal install that an existing env gets too
#   - .NET 8 (NETCore + AspNetCore) via ensure_dotnet8.sh, into ~/.proteomics-pipeline/dotnet8:
#                         DIA-NN's .raw reader and the resolution reader both run on it
#
# What stays special-cased (handled elsewhere / reported):
#   - DIA-NN: license-gated, no conda. Linux -> binary (acquire_tools.sh);
#             HIVE -> the Core's native builds; macOS -> Docker (see notes below).
#
# Outputs:
#   ~/.proteomics-pipeline/activate.sh   <- source this; puts the env on PATH
#                                           (and ensure_dotnet8.sh's .NET, once it exists)
#   ~/.proteomics-pipeline/setup.json    <- machine-readable readiness report
#                                           (ready_for.thermo_raw + thermo_raw_reader: can
#                                           step 2 read Thermo .raw here, and if not, the fix)
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
# this script's own directory, without dirname (not every PATH this runs under has it)
case "${BASH_SOURCE[0]}" in */*) SCRIPT_DIR="${BASH_SOURCE[0]%/*}" ;; *) SCRIPT_DIR="." ;; esac
SCRIPT_DIR="$(cd "$SCRIPT_DIR" && pwd)"

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
              r-ggplot2 r-ggrepel r-pheatmap r-jsonlite   # publication-quality figures
              pandoc python-docx)   # pandoc + python-docx: Markdown report -> Word .docx
  # proteowizard (msconvert) is bioconda LINUX-only
  if [ "$OS" = "linux" ]; then pkgs+=(proteowizard); fi

  say "[setup] solving + installing the analysis environment (this can take a few minutes)..."
  # -p, not -n: a user .condarc `envs_dirs` outranks -r for a NAMED env. On HIVE (brettsp,
  # 2026-09-23) envs_dirs is the relative "quobyte/proteomics-grp/conda_envs/", so `-r ROOT -n
  # NAME` built the env under $PWD/quobyte/... -- env_ready never found it, and setup.json said
  # rscript "" after a successful create.
  case "$CONDA" in
    *micromamba)
      "$CONDA" create -y -r "$MAMBA_ROOT" -p "$ENV_PREFIX" \
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
    create_env || NOTES+=("Environment solve failed. Try: $CONDA create -r $MAMBA_ROOT -p $ENV_PREFIX -c conda-forge -c bioconda bioconductor-limpa sage-proteomics python pyarrow")
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

# ---- 2b. ThermoRawFileParser + pythonnet + pandas: their own install ---------
# gabrig 2026-09-23 (15 Fusion Lumos .raw, HIVE): nothing had installed a parser, step 2 read
# every file as "not found" while setup.json said ready_for.dia, and the search would have run
# the 380-980 FALLBACK on a method that acquired 357-1105. bioconda's build is the
# self-contained .NET 8 one (run deps icu, libzlib, openssl, wget -- no dotnet), for linux-64,
# osx-64 and osx-arm64. A SEPARATE install, and not fatal: in create_env's one solve a
# thermorawfileparser that does not resolve (no build for this platform, a channel hiccup)
# would take R and limpa down with it.
# pythonnet (conda-forge, noarch; deps clr_loader + cffi, no .NET) goes in the same step:
# thermo_resolution.py uses it to read the Orbitrap MS1/MS2 resolution from the scan trailer
# with the RawFileReader DLLs TRFP ships -- TRFP never outputs it, and without it
# estimate_params.py cannot pin DIA-NN's documented tolerances. Its .NET is ensure_dotnet8.sh's
# (a self-contained TRFP's bundled runtime cannot host it). Installed = its conda-meta record.
# pandas (conda-forge) rides along: gabrig's follow-up (#15) found none in the HIVE env, and
# ad-hoc analysis reaches for it first. Here, not in create_env, because create_env only runs
# for a MISSING env -- this step is how an env that already exists gets it on a re-run.
trfp_in_env() { [ -x "$ENV_PREFIX/bin/ThermoRawFileParser" ] || [ -x "$ENV_PREFIX/bin/thermorawfileparser" ]; }
in_env_meta() {     # in_env_meta PKG: conda's own record of an installed package
  local f; for f in "$ENV_PREFIX"/conda-meta/"$1"-[0-9]*.json; do [ -e "$f" ] && return 0; done
  return 1
}
EXTRA_PKGS=""       # a word list, not an array: bash 3.2 + set -u and empty arrays do not mix
trfp_in_env           || EXTRA_PKGS="$EXTRA_PKGS thermorawfileparser"
in_env_meta pythonnet || EXTRA_PKGS="$EXTRA_PKGS pythonnet"
in_env_meta pandas    || EXTRA_PKGS="$EXTRA_PKGS pandas"
if [ -n "$CONDA" ] && env_ready && [ -n "$EXTRA_PKGS" ] && ! $CHECK_ONLY; then
  say "[setup] installing$EXTRA_PKGS into the env..."
  # shellcheck disable=SC2086  # $EXTRA_PKGS is meant to split into package names
  case "$CONDA" in      # -p for the reason create_env gives
    *micromamba) "$CONDA" install -y -r "$MAMBA_ROOT" -p "$ENV_PREFIX" \
                   -c conda-forge -c bioconda $EXTRA_PKGS >&2 ;;
    *)           "$CONDA" install -y -p "$ENV_PREFIX" \
                   -c conda-forge -c bioconda $EXTRA_PKGS >&2 ;;
  esac || NOTES+=("Could not install$EXTRA_PKGS into the env (the rest of the env is unaffected). Thermo .raw needs thermorawfileparser, and reading the Orbitrap resolution needs pythonnet: see thermo_raw_reader in setup.json.")
fi

# ---- 2c. .NET 8 (ensure_dotnet8.sh) ------------------------------------------
# One root with Microsoft.NETCore.App >= 8.0.17 + AspNetCore serves DIA-NN's .raw reader, a
# framework-dependent ThermoRawFileParser and the resolution reader. Left to a manual step, most
# Thermo users met a not-ready gate at step 2 instead. Idempotent (it reuses a good root at once),
# and not fatal: without it the note says so, and thermo_raw_reader reports what cannot run.
# Not with --check (it downloads), and only on Linux/macOS -- under Git Bash, DIA-NN for Windows
# reads .raw with no .NET step (references/install.md). activate.sh exports the root it made.
# PROTEOMICS_ENSURE_DOTNET8 exists so the tests can stand a stub in for the download.
ENSURE_DOTNET8="${PROTEOMICS_ENSURE_DOTNET8:-$SCRIPT_DIR/ensure_dotnet8.sh}"
DOTNET8_ROOT=""; DOTNET8_NOTE="not run: setup.sh --check installs nothing"
if ! $CHECK_ONLY; then
  case "$OS" in
    linux|darwin)
      say "[setup] making sure .NET 8 is available (ensure_dotnet8.sh)..."
      if DOTNET8_OUT="$("${BASH:-bash}" "$ENSURE_DOTNET8")" && [ -n "$DOTNET8_OUT" ]; then
        DOTNET8_ROOT="${DOTNET8_OUT##*$'\n'}"      # its contract: the root is the LAST line
        DOTNET8_NOTE="ensure_dotnet8.sh: .NET 8 at $DOTNET8_ROOT"
      else
        DOTNET8_NOTE="ensure_dotnet8.sh could not provide .NET 8 (it needs internet the first time; its own messages are on setup.sh's stderr)"
        NOTES+=("$DOTNET8_NOTE. Until it does, DIA-NN cannot read .raw and the Orbitrap resolution cannot be read: re-run \`bash $ENSURE_DOTNET8\` where there is internet (a login node).")
      fi ;;
    *) DOTNET8_NOTE="not run on $OS" ;;
  esac
fi

# ---- 3. resolve tool paths --------------------------------------------------
resolve() { [ -x "$ENV_PREFIX/bin/$1" ] && echo "$ENV_PREFIX/bin/$1" || (have "$1" && command -v "$1" || echo ""); }
PY="$(resolve python)";     [ -z "$PY" ] && PY="$(command -v python3 || true)"
RSCRIPT="$(resolve Rscript)"
SAGE="$(resolve sage)"
MSCONVERT="$(resolve msconvert)"
HAS_DOCKER=false; have docker && HAS_DOCKER=true
HAS_APPTAINER=false; ( have apptainer || have singularity ) && HAS_APPTAINER=true
# The Core's shared folder on HIVE. Overridable so the DIA-NN reachability below can be
# tested against a mocked layout, exactly as acquire_tools.sh's DIANN_HIVE_DIR is; nobody
# else needs to set it.
QUOBYTE_DIR="${QUOBYTE_DIR-/quobyte/proteomics-grp}"
QUOBYTE=false; [ -n "$QUOBYTE_DIR" ] && [ -d "$QUOBYTE_DIR" ] && QUOBYTE=true

# DIA-NN reachability by platform
DIANN_PATH="diann_engine"; DIANN_READY=false; DIANN_NOTE=""
# HIVE's DIA-NN builds are native binaries (build_<nnn>/diann-<version>/diann-linux); only
# the older 2.3.0 is a .sif, so reusing them does not need apptainer -- that is why the
# apptainer test that used to be here is gone. The LINUX test is not optional though:
# /quobyte can be mounted on a Mac, where no diann-linux and no .sif can run at all. On
# `$QUOBYTE` alone such a host reported diann.ready=true and never saw the Docker Desktop
# instructions it actually needs.
if   [ "$OS" = "linux" ] && $QUOBYTE; then DIANN_READY=true;  DIANN_NOTE="HIVE: reuse the Core's native DIA-NN builds under $QUOBYTE_DIR/dia-nn (acquire_tools.sh resolves the pinned version, or downloads it if the Core has none)."
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

# Thermo .raw: can step 2 read it HERE? Asked of detect_acquisition.py itself, so the parser
# search (env var, PATH, the env, the Core's shared copy) and the .NET check are the ones step
# 2 will run, not a second copy of them. Env bin first on PATH, as after activate.sh. It starts
# the parser once (`--version`) and reads no .raw. Exit 0 = ready; stdout = the object.
THERMO_READY=false; THERMO_READER=""
if [ -n "$PY" ]; then
  THERMO_READER="$(PATH="$ENV_PREFIX/bin:$PATH" PROTEOMICS_PIPELINE_HOME="$PP_HOME" \
                   "$PY" "$SCRIPT_DIR/detect_acquisition.py" --check-reader 2>/dev/null)" \
    && THERMO_READY=true
fi
case "$THERMO_READER" in "{"*) ;; *) THERMO_READY=false; THERMO_READER="" ;; esac
$THERMO_READY || NOTES+=("Thermo .raw cannot be read yet (only matters for .raw input): see thermo_raw_reader.note in setup.json for the exact fix.")

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
# .NET 8 from ensure_dotnet8.sh (NETCore >= 8.0.17 + AspNetCore): DIA-NN 2.6 reads .raw with
# it, and so does a framework-dependent ThermoRawFileParser. A DOTNET_ROOT already set is kept.
_pp_dotnet="\${PROTEOMICS_DOTNET_DIR:-\$HOME/.proteomics-pipeline/dotnet8}"
if [ -z "\${DOTNET_ROOT:-}" ] && [ -x "\$_pp_dotnet/dotnet" ]; then
  export DOTNET_ROOT="\$_pp_dotnet"; export PATH="\$DOTNET_ROOT:\$PATH"
fi
unset _pp_dotnet
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
  printf '  "dotnet8": {"root": "%s", "note": "%s"},\n' "$(j "$DOTNET8_ROOT")" "$(j "$DOTNET8_NOTE")"
  printf '  "ready_for": {"de": %s, "dia": %s, "dda": %s, "thermo_raw": %s},\n' \
         "$($DE_READY && echo true || echo false)" \
         "$($DIA_READY && echo true || echo false)" \
         "$($DDA_READY && echo true || echo false)" \
         "$($THERMO_READY && echo true || echo false)"
  if [ -n "$THERMO_READER" ]; then printf '  "thermo_raw_reader": %s,\n' "$THERMO_READER"
  else printf '  "thermo_raw_reader": {"ready": false, "note": "%s"},\n' \
         "$(j "Could not ask detect_acquisition.py whether Thermo .raw can be read (python: ${PY:-none found}). Re-run setup.sh to build the env, then setup.sh --check.")"
  fi
  printf '  "notes": ['
  for i in "${!NOTES[@]}"; do
    printf '%s"%s"' "$( [ "$i" -gt 0 ] && echo ', ' )" "$(j "${NOTES[$i]}")"
  done
  printf ']\n}\n'
} | tee "$SETUP_JSON"

say ""
say "[setup] activate with:  source $ACTIVATE"
say "[setup] readiness report written to $SETUP_JSON"
