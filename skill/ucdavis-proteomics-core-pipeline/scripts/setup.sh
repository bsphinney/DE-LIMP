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
#   - a conda env with:   python + pyarrow + pyyaml + matplotlib-base (the peak check's pictures)
#                         R (>=4.5) + bioconductor-limma + limpa >= 1.4.0 (from
#                                   Bioconductor 3.23 in step 2a: bioconda's only build,
#                                   1.2.5, cannot be solved with the limma pin)
#                                   + r-arrow + r-dplyr + r-tidyr + r-nanoparquet
#                         sage-proteomics            (the DDA search engine)
#                         proteowizard/msconvert     (LINUX ONLY on bioconda)
#   - thermorawfileparser (bioconda, self-contained: reads Thermo .raw for
#                         detect_acquisition.py) + pythonnet (conda-forge: lets
#                         thermo_resolution.py read the Orbitrap resolution) + pandas
#                         + matplotlib-base -- one separate, non-fatal install that an
#                         existing env gets too
#   - AnnotationDbi + GO.db + org.Hs.eg.db + org.Mm.eg.db: offline GO sets for run_sets.R
#                         (AnnotationDbi from bioconda; the annotation packages from Bioconductor's
#                         annotation repository for the env's release; each must load; non-fatal)
#   - .NET 8 (NETCore + AspNetCore) via ensure_dotnet8.sh, into $PP_HOME/dotnet8:
#                         DIA-NN's .raw reader and the resolution reader both run on it
#   - google-genai >= 2.25.0 (conda-forge): the podcast's Gemini voices (make_podcast.py, step
#                         9e) -- its own, non-fatal install
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
# On a cluster login node either one runs itself again on a compute node (on_compute.py: one
# blocking srun) -- never the solve or the R loads on the login node. Exit 2: it could not.
# =============================================================================
set -uo pipefail

PP_HOME="${PP_HOME:-$HOME/.proteomics-pipeline}"
# absolute, so the conda calls below can run from inside it (in_prefix)
mkdir -p "$PP_HOME" && PP_HOME="$(cd "$PP_HOME" && pwd)"
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
# sha256 of what is downloaded, against the vendor's published value (verify_download.sh)
. "$SCRIPT_DIR/verify_download.sh"

# The analysis environment's versions are PINNED, not "whatever is current": two installs a month
# apart used to be able to give different DE numbers. These are the versions the 2.15 results were
# produced and checked with, or -- where the record names only a floor -- what bioconda/conda-forge
# carried for R 4.5 on 2026-11-17 (every pin below has builds for osx-arm64, osx-64 and linux-64):
#   limma 3.68.5   the limma the delivered v2 tables were reproduced with (references/environment.md)
#   limpa 1.4.0    installed in step 2a from Bioconductor 3.23 (bioconda has only 1.2.5)
#   statmod 1.5.1, arrow 24.0.0 (R and Python, one libarrow), sage 0.14.7 (resolve_defaults.py's pin)
# To move to newer ones, change them HERE, re-run the tests on a real dataset, and say so in the
# CHANGELOG. PIN_* may be overridden for a test; nothing else should.
PIN_PYTHON="${PIN_PYTHON:-3.11}"
PIN_LIMMA="${PIN_LIMMA:-3.68.5}"
PIN_STATMOD="${PIN_STATMOD:-1.5.1}"
PIN_ARROW="${PIN_ARROW:-24.0.0}"
PIN_SAGE_CONDA="${PIN_SAGE_CONDA:-0.14.7}"
PIN_LIMPA="${PIN_LIMPA:-1.4.0}"
PIN_MICROMAMBA="${PIN_MICROMAMBA:-2.10.0-0}"

# ---- 0. a compute node, never a cluster login node ---------------------------
# HPC@UCD's ucd-hpc skill (mandatory on Hive, Farm and Franklin since 2026-10-16) names conda
# solves as login-node abuse, and downloads and package builds belong inside a job; even --check
# loads R and a handful of packages and starts the .NET parser. On a login node (on_compute.py
# --is-login-node: no SLURM job, and SLURM or HPCCF's cluster file here) setup.sh runs itself
# again through on_compute.py: one blocking srun with the queue spelled out, on a queue that does
# not preempt, and waits. HIVE's compute nodes reach conda-forge, Bioconductor, CRAN and
# Microsoft's .NET downloads (a compute-node job, 2026-10-16).
# Exit 3: that run did not finish (SLURM refused it, a time limit), so setup.json -- which would
# otherwise be read as current -- was not rewritten. Exit 2: no way to a compute node here.
UCD_SLURM_BIN="${HPCCF_SLURM_BIN:-/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin}"
UCD_CONF="${UCD_HPC_CLUSTER_CONF:-/opt/hpccf/etc/cluster_name.conf}"
slurm_here() { have sbatch || [ -x "$UCD_SLURM_BIN/sbatch" ] || [ -r "$UCD_CONF" ]; }   # a cluster
if [ -z "${SLURM_JOB_ID:-}" ] && [ -z "${PROTEOMICS_ON_COMPUTE:-}" ]; then
  OC_PY="$(command -v python3 2>/dev/null || true)"
  if [ -n "$OC_PY" ] && [ -f "$SCRIPT_DIR/on_compute.py" ]; then
    if "$OC_PY" "$SCRIPT_DIR/on_compute.py" --is-login-node; then
      # --check only reads (setup.json at the end): it may run on low rather than wait behind
      # the user's own searches on high; an install never does
      if $CHECK_ONLY; then OC_ARGS="--time 20 --mem 4G --cpus 1 --preemptible"; OC_WHY="setup.sh --check loads R and its packages"
      else OC_ARGS="--time 120 --mem 8G --cpus 2"; OC_WHY="setup.sh solves and installs the environment"; fi
      # the run writes this id into setup.json (srun passes the environment on): the one that is
      # there afterwards is this run's, or the run did not finish
      PROTEOMICS_SETUP_RUN_ID="$$.${RANDOM:-0}.${RANDOM:-0}"; export PROTEOMICS_SETUP_RUN_ID
      # shellcheck disable=SC2086  # OC_ARGS is meant to split
      "$OC_PY" "$SCRIPT_DIR/on_compute.py" $OC_ARGS --name pp_setup --why "$OC_WHY" \
        -- "${BASH:-bash}" "$SCRIPT_DIR/setup.sh" "$@"
      OC_RC=$?
      OC_WROTE="$(sed -n 's/.*"run_id": "\([^"]*\)".*/\1/p' "$SETUP_JSON" 2>/dev/null)"
      if [ "$OC_WROTE" != "$PROTEOMICS_SETUP_RUN_ID" ]; then
        say "[setup] ERROR: the compute-node run did not finish (exit $OC_RC), and $SETUP_JSON was"
        say "[setup] not rewritten: its ready_for is from an earlier run. Read the srun lines above"
        say "[setup] (an account or partition SLURM refused, a time limit), then run setup.sh again."
        exit 3
      fi
      exit $OC_RC
    fi
  elif slurm_here; then
    say "[setup] REFUSING to run on a cluster login node (SLURM is here, SLURM_JOB_ID is unset),"
    say "[setup] and there is no python3 + on_compute.py here to move it to a compute node. Run it in a"
    say "[setup] job, with your account and partition: srun --account=<yours> --partition=<yours> \\"
    say "[setup]   --time=120 --mem=8G --cpus-per-task=2 bash $SCRIPT_DIR/setup.sh $*"
    exit 2
  fi
fi

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
  say "[setup] installing micromamba $PIN_MICROMAMBA (no admin needed) for $MM_PLAT ..."
  mkdir -p "$MAMBA_ROOT/bin"
  # A PINNED release from the project's GitHub releases, checked against the sha256 GitHub shows
  # for that asset -- not `micro.mamba.pm/.../latest` piped into tar, which ran whatever the
  # endpoint served that day. The asset is the program itself (no archive).
  local url="https://github.com/mamba-org/micromamba-releases/releases/download/$PIN_MICROMAMBA/micromamba-$MM_PLAT"
  local part="$MAMBA_ROOT/bin/micromamba.part"
  rm -f "$part"
  if curl -fsSL "$url" -o "$part"; then
    if verify_download "micromamba-$PIN_MICROMAMBA-$MM_PLAT" "$part" "micromamba $PIN_MICROMAMBA"; then :
    else
      local rc=$?
      if [ "$rc" = 1 ]; then rm -f "$part"; NOTES+=("The micromamba download did not match its published sha256, so it was deleted and not used. Re-run setup.sh; if it happens again, install micromamba yourself: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html")
        return; fi
      NOTES+=("micromamba $PIN_MICROMAMBA was downloaded but NOT verified (no pinned checksum for $MM_PLAT).")
    fi
    chmod +x "$part" && mv -f "$part" "$MAMBA_ROOT/bin/micromamba"
    CONDA="$MAMBA_ROOT/bin/micromamba"
    say "[setup] micromamba installed at $CONDA"
  else
    rm -f "$part"
    NOTES+=("Could not download micromamba automatically. Install it manually: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html  then re-run setup.sh.")
  fi
}

pick_conda
if [ -z "$CONDA" ]; then
  if $CHECK_ONLY; then NOTES+=("No conda/micromamba found; run setup.sh (without --check) to install it.")
  else install_micromamba; fi
fi

# The package cache, an absolute folder under this prefix (2.15.1, HIVE F9): a user ~/.condarc
# with a RELATIVE pkgs_dirs (HIVE, a Core admin's account: "quobyte/proteomics-grp/conda_pkgs/")
# made micromamba create that folder under whatever directory setup.sh was run from. libmamba
# CONCATENATES a list setting from every source (configuration_impl.hpp, Source<std::vector>::
# merge), env first, so CONDA_PKGS_DIRS alone makes the cache absolute but cannot remove the
# rc's entry; in_prefix runs each conda call from inside the prefix, so a relative entry
# resolves there and never in the user's folder. A CONDA_PKGS_DIRS the user set to absolute
# folders is kept.
case ":${CONDA_PKGS_DIRS:-}" in
  ":"|*:[!/]*|*,[!/]*) export CONDA_PKGS_DIRS="$MAMBA_ROOT/pkgs" ;;
esac
in_prefix() { mkdir -p "$MAMBA_ROOT" && ( cd "$MAMBA_ROOT" && "$@" ); }
# create-env helper that works for micromamba OR conda/mamba
create_env() {
  # NO bioconductor-limpa here (2.15.1 review, blocker): bioconda's only build, 1.2.5, requires
  # bioconductor-limma >=3.66.0,<3.67.0, so it and the limma pin cannot be solved together
  # (conda --dry-run: LibMambaUnsatisfiableError on osx-arm64 and linux-64). limpa comes from
  # step 2a only -- Bioconductor's source package, pure R, whose Depends/Imports (limma,
  # data.table, statmod) are all in this list. tests/test_setup_pins_2151.py keeps it so.
  local pkgs=("python=$PIN_PYTHON" "pyarrow=$PIN_ARROW" pyyaml matplotlib-base   # matplotlib: the peak check's pictures
              "r-base>=4.5" "bioconductor-limma=$PIN_LIMMA" "r-statmod=$PIN_STATMOD" r-data.table
              r-nanoparquet "r-arrow=$PIN_ARROW" r-dplyr r-tidyr "sage-proteomics=$PIN_SAGE_CONDA"
              r-ggplot2 r-ggrepel r-pheatmap r-jsonlite   # publication-quality figures
              pandoc python-docx)   # pandoc + python-docx: Markdown report -> Word .docx
  # proteowizard (msconvert) is bioconda LINUX-only
  if [ "$OS" = "linux" ]; then pkgs+=(proteowizard); fi

  say "[setup] solving + installing the analysis environment (this can take a few minutes)..."
  # -p, not -n: a user .condarc `envs_dirs` outranks -r for a NAMED env. On HIVE (a Core admin's account,
  # 2026-09-23) envs_dirs is the relative "quobyte/proteomics-grp/conda_envs/", so `-r ROOT -n
  # NAME` built the env under $PWD/quobyte/... -- env_ready never found it, and setup.json said
  # rscript "" after a successful create.
  # --override-channels: the skill names its own channels. A user .condarc with `defaults`
  # and channel_priority strict (HIVE, a Core admin's account, in a test job) made libmamba consider only
  # defaults' r-statmod, built for R <= 4.3, and bioconductor-limma 3.66 was unsolvable.
  case "$CONDA" in
    *micromamba)
      in_prefix "$CONDA" create -y -r "$MAMBA_ROOT" -p "$ENV_PREFIX" \
        --override-channels -c conda-forge -c bioconda "${pkgs[@]}" ;;
    *)
      in_prefix "$CONDA" create -y -p "$ENV_PREFIX" \
        --override-channels -c conda-forge -c bioconda "${pkgs[@]}" ;;
  esac
}

# ---- 2. ensure the environment exists ---------------------------------------
env_ready() { [ -x "$ENV_PREFIX/bin/python" ] && [ -x "$ENV_PREFIX/bin/Rscript" ]; }

if [ -n "$CONDA" ] && ! env_ready; then
  if $CHECK_ONLY; then NOTES+=("Analysis env not built yet; run setup.sh to create it.")
  else
    create_env || NOTES+=("Environment solve failed. Try: $CONDA create -r $MAMBA_ROOT -p $ENV_PREFIX --override-channels -c conda-forge -c bioconda bioconductor-limma=$PIN_LIMMA sage-proteomics python pyarrow")
  fi
fi

# ---- 2a. limpa >= 1.4.0 -------------------------------------------------------
# run_de.R reads DIA-NN with readDIANN(annotation.columns =): limpa >= 1.4.0, Bioconductor
# 3.23. bioconda's only build is 1.2.5 (Bioconductor 3.22, R 4.5), and a 2.8.0 run_de.R died
# on it (in a test job on HIVE). limpa is pure R -- github.com/bioc/limpa RELEASE_3_23 has no
# src/, Depends limma with no version and no R floor, Imports data.table / statmod -- and
# every limma function it imports (lmFit, duplicateCorrelation, squeezeVar, arrayWeights,
# asMatrixWeights, chooseLowessSpan, loessFit, plotMDS) is in bioconda's limma 3.66. So the
# env's R takes the Bioconductor 3.23 SOURCE package as it is: no compiler, no R 4.6, no
# second Bioconductor stack (conda-forge has r-base 4.6 but no R 4.6 builds of r-statmod,
# r-arrow, ... yet). An env that already has R 4.6 + Bioconductor 3.23 (HIVE's
# proteomics-pipeline-r46) passes the same check. PROTEOMICS_LIMPA_BIOC overrides the release.
LIMPA_MIN="1.4.0"
LIMPA_PIN="$PIN_LIMPA"
LIMPA_BIOC="${PROTEOMICS_LIMPA_BIOC:-3.23}"
LIMPA_REPO="https://bioconductor.org/packages/$LIMPA_BIOC/bioc"
limpa_ok() {
  "$ENV_PREFIX/bin/Rscript" -e "q(status = as.integer(!(requireNamespace('limpa', quietly = TRUE) && utils::packageVersion('limpa') >= '$LIMPA_MIN')))" >/dev/null 2>&1
}
LIMPA_JUST_INSTALLED=false
if env_ready && ! $CHECK_ONLY && ! limpa_ok; then
  LIMPA_JUST_INSTALLED=true
  say "[setup] installing limpa >= $LIMPA_MIN from Bioconductor $LIMPA_BIOC (bioconda has only 1.2.5)..."
  "$ENV_PREFIX/bin/Rscript" -e "install.packages('limpa', repos = c(BioCsoft = '$LIMPA_REPO', CRAN = 'https://cloud.r-project.org'), type = 'source', dependencies = NA)" >&2 \
    || say "[setup] the limpa install failed (above)"
fi
LIMPA_VERSION=""; LIMPA_OK=false
if env_ready; then
  LIMPA_VERSION="$("$ENV_PREFIX/bin/Rscript" -e "cat(tryCatch(as.character(utils::packageVersion('limpa')), error = function(e) ''))" 2>/dev/null)"
  # limpa is PINNED to $LIMPA_PIN (the release the 2.15 results were produced with). The install above
  # takes whatever Bioconductor $LIMPA_BIOC holds now; when that is another version, fetch the pinned
  # source tarball (current release folder, then the archive) and install that one over it. If neither
  # has it, the run goes on with what is there and says so -- limpa_version is recorded either way.
  # Only for a limpa THIS run installed: a working env that already has a newer one is not downgraded.
  if $LIMPA_JUST_INSTALLED && [ -n "$LIMPA_VERSION" ] && [ "$LIMPA_VERSION" != "$LIMPA_PIN" ]; then
    say "[setup] limpa is $LIMPA_VERSION, the pinned release is $LIMPA_PIN: installing the pinned source ..."
    "$ENV_PREFIX/bin/Rscript" -e "for (u in c('$LIMPA_REPO/src/contrib/limpa_$LIMPA_PIN.tar.gz', '$LIMPA_REPO/src/contrib/Archive/limpa/limpa_$LIMPA_PIN.tar.gz')) { ok <- tryCatch({ install.packages(u, repos = NULL, type = 'source'); as.character(utils::packageVersion('limpa')) == '$LIMPA_PIN' }, error = function(e) FALSE, warning = function(w) FALSE); if (isTRUE(ok)) break }" >&2 || true
    LIMPA_VERSION="$("$ENV_PREFIX/bin/Rscript" -e "cat(tryCatch(as.character(utils::packageVersion('limpa')), error = function(e) ''))" 2>/dev/null)"
  fi
  if [ -n "$LIMPA_VERSION" ] && [ "$LIMPA_VERSION" != "$LIMPA_PIN" ]; then
    NOTES+=("The analysis environment has limpa $LIMPA_VERSION, not the pinned $LIMPA_PIN (the release the 2.15 results were produced with). Results can differ slightly; the version is recorded with every analysis (de_provenance.json).")
  fi
  limpa_ok && LIMPA_OK=true
fi
# What THIS limpa reads report.parquet with. Its own DESCRIPTION says: Suggests arrow in the 1.2.x
# line (readDIANN: arrow::read_parquet), nanoparquet from 1.4 (github.com/bioc/limpa RELEASE_3_22
# vs RELEASE_3_23) -- right for a devel limpa between the two as well. The version decides only
# when the DESCRIPTION names neither.
LIMPA_READS=""
if env_ready && { $LIMPA_OK || [ -n "$LIMPA_VERSION" ]; }; then
  LIMPA_READS="$("$ENV_PREFIX/bin/Rscript" -e "d <- tryCatch(utils::packageDescription('limpa')[['Suggests']], error = function(e) NULL); d <- if (is.null(d)) '' else d; cat(if (grepl('nanoparquet', d)) 'nanoparquet' else if (grepl('arrow', d)) 'arrow' else '')" 2>/dev/null)"
  case "$LIMPA_READS" in
    nanoparquet|arrow) ;;
    *) if $LIMPA_OK; then LIMPA_READS=nanoparquet; else LIMPA_READS=arrow; fi ;;
  esac
fi
LIMPA_FIX="bash $SCRIPT_DIR/setup.sh (it installs it; on a cluster it runs itself on a compute node), or by hand -- on a cluster inside a job, never on a login node: $ENV_PREFIX/bin/Rscript -e \"install.packages('limpa', repos = '$LIMPA_REPO', type = 'source')\", then bash setup.sh --check"
LIMPA_STATE="limpa ${LIMPA_VERSION:-is not installed}${LIMPA_VERSION:+ is older than $LIMPA_MIN}"
if env_ready && ! $LIMPA_OK; then
  NOTES+=("$LIMPA_STATE: run_de.R --method dpc needs limpa >= $LIMPA_MIN (it can still read the report with limpa 1.2.x and records that it did, but the pinned stack is >= $LIMPA_MIN). Fix: $LIMPA_FIX.")
fi

# ---- 2b. ThermoRawFileParser + pythonnet + pandas: their own install ---------
# A Core member, 2026-09-23 (15 Fusion Lumos .raw, HIVE): nothing had installed a parser, step 2 read
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
# pandas (conda-forge) rides along: a Core member's follow-up (#15) found none in the HIVE env, and
# ad-hoc analysis reaches for it first. Here, not in create_env, because create_env only runs
# for a MISSING env -- this step is how an env that already exists gets it on a re-run.
trfp_in_env() { [ -x "$ENV_PREFIX/bin/ThermoRawFileParser" ] || [ -x "$ENV_PREFIX/bin/thermorawfileparser" ]; }
in_env_meta() {     # in_env_meta PKG: conda's own record of an installed package
  local f; for f in "$ENV_PREFIX"/conda-meta/"$1"-[0-9]*.json; do [ -e "$f" ] && return 0; done
  return 1
}
nanoparquet_ok() {
  "$ENV_PREFIX/bin/Rscript" -e "q(status = as.integer(!requireNamespace('nanoparquet', quietly = TRUE)))" >/dev/null 2>&1
}
EXTRA_PKGS=""       # a word list, not an array: bash 3.2 + set -u and empty arrays do not mix
trfp_in_env           || EXTRA_PKGS="$EXTRA_PKGS thermorawfileparser"
in_env_meta pythonnet || EXTRA_PKGS="$EXTRA_PKGS pythonnet"
in_env_meta pandas    || EXTRA_PKGS="$EXTRA_PKGS pandas"
# matplotlib-base: skyline_peaks.py draws the peak check's pictures with it (step 8f); an env
# built before 2.14 has none, so it rides along here too
in_env_meta matplotlib-base || EXTRA_PKGS="$EXTRA_PKGS matplotlib-base"
if [ -n "$CONDA" ] && env_ready && [ -n "$EXTRA_PKGS" ] && ! $CHECK_ONLY; then
  say "[setup] installing$EXTRA_PKGS into the env..."
  # shellcheck disable=SC2086  # $EXTRA_PKGS is meant to split into package names
  case "$CONDA" in      # -p for the reason create_env gives
    *micromamba) in_prefix "$CONDA" install -y -r "$MAMBA_ROOT" -p "$ENV_PREFIX" \
                   --override-channels -c conda-forge -c bioconda $EXTRA_PKGS >&2 ;;
    *)           in_prefix "$CONDA" install -y -p "$ENV_PREFIX" \
                   --override-channels -c conda-forge -c bioconda $EXTRA_PKGS >&2 ;;
  esac || NOTES+=("Could not install$EXTRA_PKGS into the env (the rest of the env is unaffected). Thermo .raw needs thermorawfileparser, and reading the Orbitrap resolution needs pythonnet: see thermo_raw_reader in setup.json.")
fi

# ---- 2b1. nanoparquet: its OWN install ----------------------------------------
# limpa reads DIA-NN's report.parquet with it ("nanoparquet package required but is not
# installed", limpa's EListFromLongFormatFile), but it is only in limpa's Suggests -- installing
# limpa never brings it -- and create_env lists it from 2.8.0 only, so an env built before that
# never got it (2.8.0 report). Loadable is the test, not a conda-meta record: an R-installed copy
# counts.
#   * conda-forge builds r-nanoparquet only for r-base < 4.6 (0.5.2; anaconda.org, 2026-10-15).
#     Solved in 2b's call above, without --freeze-installed, an R 4.6.1 env got the plan
#     "r-base 4.6.1 -> 4.5.3" (2.12 review; reproduced with a dry run on HIVE's
#     proteomics-pipeline-r46): R downgraded under limpa and every Bioconductor package. So it is
#     installed alone and --freeze-installed: conda may only add it, or fail (the same dry run
#     frozen: "r-nanoparquet ... is not installable", nothing changed).
#   * R >= 4.6 has no conda build, so it is CRAN's source package, compiled. 2.8.0 ruled CRAN out
#     because install.packages("nanoparquet") failed with "x86_64-conda-linux-gnu-c++: not found":
#     conda's R compiles with the compilers its Makeconf names (R CMD config CC / CXX), which an
#     env has only with conda-forge's c-compiler + cxx-compiler. So CRAN is used only when those
#     run here; otherwise the fix installs them (frozen) first. HIVE's proteomics-pipeline-r46 is
#     built exactly so: conda-forge's c-compiler + cxx-compiler, nanoparquet 0.5.2 from CRAN.
R_VERSION=""
env_ready && R_VERSION="$("$ENV_PREFIX/bin/Rscript" -e 'cat(as.character(getRversion()))' 2>/dev/null)"
NANOPARQUET_CONDA=true      # does conda-forge build r-nanoparquet for this R? (unknown R: try it)
case "$R_VERSION" in ""|[0-3].*|4.[0-5]|4.[0-5].*) ;; *) NANOPARQUET_CONDA=false ;; esac
r_compilers_ok() {          # can R's own C and C++ compilers and make run, the env's bin first?
  local v c
  for v in CC CXX MAKE; do
    c="$(PATH="$ENV_PREFIX/bin:$PATH" "$ENV_PREFIX/bin/R" CMD config "$v" 2>/dev/null)"; c="${c%% *}"
    [ -n "$c" ] && ( PATH="$ENV_PREFIX/bin:$PATH"; command -v "$c" >/dev/null 2>&1 ) || return 1
  done
}
NANOPARQUET_CRAN_R="install.packages('nanoparquet', repos = 'https://cloud.r-project.org', type = 'source')"
NANOPARQUET_CRAN="PATH=\"$ENV_PREFIX/bin:\$PATH\" $ENV_PREFIX/bin/Rscript -e \"$NANOPARQUET_CRAN_R\""
case "$CONDA" in            # -p for the reason create_env gives
  *micromamba) NANOPARQUET_CONDA_CMD="$CONDA install -y -r $MAMBA_ROOT -p $ENV_PREFIX --freeze-installed --override-channels -c conda-forge -c bioconda" ;;
  "")          NANOPARQUET_CONDA_CMD="" ;;
  *)           NANOPARQUET_CONDA_CMD="$CONDA install -y -p $ENV_PREFIX --freeze-installed --override-channels -c conda-forge -c bioconda" ;;
esac
if env_ready && ! nanoparquet_ok && ! $CHECK_ONLY; then
  if $NANOPARQUET_CONDA && [ -n "$CONDA" ]; then
    say "[setup] installing r-nanoparquet into the env (conda-forge; nothing already there changes)..."
    case "$CONDA" in
      *micromamba) in_prefix "$CONDA" install -y -r "$MAMBA_ROOT" -p "$ENV_PREFIX" --freeze-installed \
                     --override-channels -c conda-forge -c bioconda r-nanoparquet >&2 ;;
      *)           in_prefix "$CONDA" install -y -p "$ENV_PREFIX" --freeze-installed \
                     --override-channels -c conda-forge -c bioconda r-nanoparquet >&2 ;;
    esac || say "[setup] conda could not add r-nanoparquet without changing the env (above)"
  fi
  if ! nanoparquet_ok && r_compilers_ok; then
    say "[setup] installing nanoparquet from CRAN's source package (R $R_VERSION, compiled with R's own compilers)..."
    PATH="$ENV_PREFIX/bin:$PATH" "$ENV_PREFIX/bin/Rscript" -e "$NANOPARQUET_CRAN_R" >&2 \
      || say "[setup] the nanoparquet install from CRAN failed (above)"
  fi
fi
NANOPARQUET_OK=false
if env_ready && nanoparquet_ok; then NANOPARQUET_OK=true; fi
# The fix that works for THIS env (setup.json nanoparquet.fix, the note, the ERROR line)
if $NANOPARQUET_CONDA; then
  NANOPARQUET_WHY="conda-forge's binary, added without changing anything else in the env"
  if [ -n "$CONDA" ]; then NANOPARQUET_FIX="$NANOPARQUET_CONDA_CMD r-nanoparquet"
  else NANOPARQUET_FIX="run setup.sh (it installs micromamba), or install conda-forge's r-nanoparquet into $ENV_PREFIX with --freeze-installed"; fi
elif env_ready && r_compilers_ok; then
  NANOPARQUET_WHY="R $R_VERSION: conda-forge has no r-nanoparquet for it, so CRAN's source package, compiled with the env's own compilers"
  NANOPARQUET_FIX="$NANOPARQUET_CRAN"
else
  NANOPARQUET_WHY="R $R_VERSION: conda-forge has no r-nanoparquet for it, and CRAN's source package needs the compilers this R names (R CMD config CC / CXX), which the env does not have -- conda-forge's c-compiler + cxx-compiler provide them"
  if [ -n "$CONDA" ]; then NANOPARQUET_FIX="$NANOPARQUET_CONDA_CMD c-compiler cxx-compiler make && $NANOPARQUET_CRAN"
  else NANOPARQUET_FIX="install conda-forge's c-compiler, cxx-compiler and make into $ENV_PREFIX (with --freeze-installed), then: $NANOPARQUET_CRAN"; fi
fi
# What a missing nanoparquet costs depends on the limpa (DE_READER below): limpa >= 1.4.0 reads
# report.parquet with it; limpa 1.2.x reads with arrow, and only run_de.R's contaminant share
# (Precursor.Quantity, read with nanoparquet) falls back to Precursor.Normalised.
if env_ready && ! $NANOPARQUET_OK; then
  if [ "$LIMPA_READS" = nanoparquet ]; then
    NANOPARQUET_COST="this env's limpa${LIMPA_VERSION:+ ($LIMPA_VERSION)} reads DIA-NN's report.parquet with it, so no DE can run"
  elif [ -z "$LIMPA_READS" ]; then
    NANOPARQUET_COST="limpa >= $LIMPA_MIN reads DIA-NN's report.parquet with it, so no DE can run"
  else
    NANOPARQUET_COST="this env's limpa ($LIMPA_VERSION) reads report.parquet with arrow, so the DE runs, but on the unpinned limpa, and run_de.R's contaminant share is read from Precursor.Normalised instead of Precursor.Quantity (de_provenance.json contaminants.intensity_column); limpa >= $LIMPA_MIN needs it"
  fi
  NOTES+=("nanoparquet cannot be loaded in $ENV_PREFIX: $NANOPARQUET_COST. The route that works here: $NANOPARQUET_WHY. Fix: $NANOPARQUET_FIX.")
fi

# ---- 2b2. protein-set annotation: GO, offline -----------------------------------
# run_sets.R tests Gene Ontology sets. They come from the Bioconductor annotation packages,
# installed here once, so a set test never downloads anything and its record names the GO
# release it used (the package's GOSOURCEDATE). org.Hs.eg.db also serves every organism
# without its own package (human gene symbols; run_sets.R records that route). About 250 MB.
# Not installed: reactome.db (about 2 GB); run_sets.R --sets reactome says how to add it.
# "Installed" is not enough: each package must LOAD. requireNamespace() said FALSE alike for a
# missing package and for conda's bioconductor-org.mm.eg.db / go.db 3.2.x builds, which
# `micromamba install --freeze-installed` resolved into a HIVE env and which did not load
# (staff report 2026-10-08: run_sets.R stopped inside a DE job). So R reports each package as
# ok / missing / broken (installed, does not load, with R's reason), and:
#   - AnnotationDbi (software; its dependencies compile): conda's binary first, as for the env;
#   - GO.db and the org.*.eg.db data packages -- and AnnotationDbi if conda could not -- from
#     Bioconductor's own repositories for the env's release, never conda's annotation builds:
#     install.packages(repos = https://bioconductor.org/packages/<release>/data/annotation) is
#     what fixed that env. The release is the one whose AnnotationDbi the env has (release X's
#     annotation packages need X's AnnotationDbi: GO.db 3.23 wants AnnotationDbi >= 1.73), else
#     its limma's (the env's own Bioconductor), else BiocVersion's -- read from their DESCRIPTION
#     (Repository: Bioconductor X.Y / git_branch RELEASE_X_Y), else $LIMPA_BIOC.
# Not fatal: DE runs without them; setup.json's set_annotation says what loads and the fix, and
# `setup.sh --check` reports it without installing anything.
SET_ANN_PKGS="AnnotationDbi GO.db org.Hs.eg.db org.Mm.eg.db"
# R's answer, one line each: "<pkg> ok <version>" | "<pkg> missing" | "<pkg> broken <reason>", then
# "bioc <X.Y> <package it was read from>" or "bioc -". (No single quote in it: bash quotes it.)
SET_ANN_R='
ok <- character(0)
for (q in strsplit(Sys.getenv("SET_ANN_PKGS"), " ")[[1]]) {
  if (!nzchar(system.file(package = q))) { cat(q, " missing\n", sep = ""); next }
  e <- tryCatch({ loadNamespace(q); "" }, error = function(e) conditionMessage(e))
  if (nzchar(e)) { cat(q, " broken ", gsub("[[:space:]]+", " ", e), "\n", sep = ""); next }
  ok <- c(ok, q)
  cat(q, " ok ", as.character(utils::packageVersion(q)), "\n", sep = "")
}
rel <- function(q) {
  d <- suppressWarnings(utils::packageDescription(q))
  if (!is.list(d)) return(NA)
  r <- sub("^Bioconductor ", "", if (is.null(d$Repository)) "" else d$Repository)
  if (grepl("^[0-9]+[.][0-9]+$", r)) return(r)
  b <- if (is.null(d$git_branch)) "" else d$git_branch
  if (grepl("^RELEASE_[0-9]+_[0-9]+$", b)) return(sub("_", ".", sub("^RELEASE_", "", b)))
  if (q == "BiocVersion") return(sub("^([0-9]+[.][0-9]+).*", "\\1", d$Version))
  NA
}
for (q in c("AnnotationDbi", "limma", "BiocVersion")) {
  if (q == "AnnotationDbi" && !(q %in% ok)) next
  r <- rel(q)
  if (!is.na(r)) { cat("bioc ", r, " ", q, "\n", sep = ""); quit(save = "no") }
}
cat("bioc -\n")
'
# -> SET_ANN_ANSWERED (did R answer), SET_ANN_MISSING (missing or broken, space-separated),
#    SET_ANN_LOADS, SET_ANN_STATUS ("pkg|status" lines), SET_ANN_PROBLEMS (for the note), ENV_BIOC(_FROM)
set_ann_check() {
  SET_ANN_ANSWERED=false; SET_ANN_MISSING=""; SET_ANN_STATUS=""; SET_ANN_PROBLEMS=""
  SET_ANN_LOADS=" "; ENV_BIOC=""; ENV_BIOC_FROM=""
  local out line pkg st rest
  out="$(SET_ANN_PKGS="$SET_ANN_PKGS" "$ENV_PREFIX/bin/Rscript" -e "$SET_ANN_R" 2>/dev/null)"
  while IFS= read -r line; do
    pkg="${line%% *}"; rest=""; case "$line" in *" "*) rest="${line#* }" ;; esac
    st="${rest%% *}"; case "$rest" in *" "*) rest="${rest#* }" ;; *) rest="" ;; esac
    if [ "$pkg" = "bioc" ]; then
      SET_ANN_ANSWERED=true
      [ "$st" != "-" ] && { ENV_BIOC="$st"; ENV_BIOC_FROM="$rest"; }
      continue
    fi
    case " $SET_ANN_PKGS " in *" $pkg "*) ;; *) continue ;; esac
    case "$st" in
      ok)      SET_ANN_LOADS="$SET_ANN_LOADS$pkg "; rest="loads ($rest)" ;;
      missing) rest="not installed" ;;
      broken)  rest="installed but does not load: $rest" ;;
      *)       continue ;;
    esac
    [ "$st" = ok ] || { SET_ANN_MISSING="$SET_ANN_MISSING${SET_ANN_MISSING:+ }$pkg"
                        SET_ANN_PROBLEMS="$SET_ANN_PROBLEMS${SET_ANN_PROBLEMS:+; }$pkg $rest"; }
    SET_ANN_STATUS="$SET_ANN_STATUS$pkg|$rest
"
  done <<SET_ANN_EOF
$out
SET_ANN_EOF
}
set_ann_loads() { case "$SET_ANN_LOADS" in *" $1 "*) return 0 ;; esac; return 1; }
# SET_ANN_STATUS as JSON members, at setup.json time (j is defined there)
set_ann_status_json() {
  local line sep=""
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    printf '%s"%s": "%s"' "$sep" "$(j "${line%%|*}")" "$(j "${line#*|}")"; sep=", "
  done <<SET_ANN_EOF
$SET_ANN_STATUS
SET_ANN_EOF
}
SET_ANN_ANSWERED=false; SET_ANN_MISSING=""; SET_ANN_STATUS=""; SET_ANN_PROBLEMS=""; SET_ANN_LOADS=" "
ENV_BIOC=""; ENV_BIOC_FROM=""
env_ready && set_ann_check
if env_ready && ! $CHECK_ONLY && [ -n "$SET_ANN_MISSING" ]; then
  if ! set_ann_loads AnnotationDbi; then
    say "[setup] installing AnnotationDbi for the protein-set tests (conda)..."
    case "$CONDA" in
      *micromamba) in_prefix "$CONDA" install -y -r "$MAMBA_ROOT" -p "$ENV_PREFIX" --override-channels \
                     -c conda-forge -c bioconda bioconductor-annotationdbi >&2 ;;
      "")          false ;;
      *)           in_prefix "$CONDA" install -y -p "$ENV_PREFIX" --override-channels \
                     -c conda-forge -c bioconda bioconductor-annotationdbi >&2 ;;
    esac || say "[setup] conda could not install AnnotationDbi; trying Bioconductor's source package"
    set_ann_check
  fi
  if [ -n "$SET_ANN_MISSING" ]; then
    SET_ANN_BIOC="${ENV_BIOC:-$LIMPA_BIOC}"
    say "[setup] installing $SET_ANN_MISSING from Bioconductor $SET_ANN_BIOC (annotation repository)..."
    "$ENV_PREFIX/bin/Rscript" -e "install.packages(strsplit('$SET_ANN_MISSING', ' ')[[1]], repos = c(BioCann = 'https://bioconductor.org/packages/$SET_ANN_BIOC/data/annotation', BioCsoft = 'https://bioconductor.org/packages/$SET_ANN_BIOC/bioc', CRAN = 'https://cloud.r-project.org'), dependencies = NA)" >&2 \
      || say "[setup] the annotation install failed (above)"
    set_ann_check
  fi
fi
SET_ANN_BIOC="${ENV_BIOC:-$LIMPA_BIOC}"
if [ -z "$ENV_BIOC" ]; then
  if $SET_ANN_ANSWERED; then ENV_BIOC_FROM="none of AnnotationDbi, limma or BiocVersion names one"
  else ENV_BIOC_FROM="not known: R did not answer"; fi
  ENV_BIOC_FROM="$ENV_BIOC_FROM -- the skill's pinned release (PROTEOMICS_LIMPA_BIOC) instead"
fi
SET_ANN_OK=false; $SET_ANN_ANSWERED && [ -z "$SET_ANN_MISSING" ] && SET_ANN_OK=true
SET_ANN_FIX=""
if env_ready && ! $SET_ANN_ANSWERED; then
  SET_ANN_FIX="re-run bash setup.sh --check; if it says the same, run $ENV_PREFIX/bin/Rscript by hand to see R's error"
  NOTES+=("Protein-set tests (run_sets.R): R in $ENV_PREFIX did not say whether the annotation packages ($SET_ANN_PKGS) load, so they are not known to work. DE does not need them.")
elif env_ready && ! $SET_ANN_OK; then
  SET_ANN_FIX="bash setup.sh (on a cluster it runs itself on a compute node), or by hand -- on a cluster inside a job, never on a login node: $ENV_PREFIX/bin/Rscript -e \"install.packages(c('$(printf '%s' "$SET_ANN_MISSING" | sed "s/ /','/g")'), repos = c('https://bioconductor.org/packages/$SET_ANN_BIOC/data/annotation', 'https://bioconductor.org/packages/$SET_ANN_BIOC/bioc', 'https://cloud.r-project.org'))\""
  NOTES+=("Protein-set tests (run_sets.R) cannot run until these load: $SET_ANN_PROBLEMS. Fix: $SET_ANN_FIX. DE does not need them.")
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
DOTNET8_ROOT=""; DOTNET8_NOTE="not run: setup.sh --check installs nothing"; DOTNET8_STATE=not_run
if ! $CHECK_ONLY; then
  case "$OS" in
    linux|darwin)
      say "[setup] making sure .NET 8 is available (ensure_dotnet8.sh)..."
      # into this prefix (PP_HOME), never $HOME's default when PP_HOME is another (HIVE F5)
      if DOTNET8_OUT="$(PROTEOMICS_PIPELINE_HOME="$PP_HOME" "${BASH:-bash}" "$ENSURE_DOTNET8")" \
          && [ -n "$DOTNET8_OUT" ]; then
        DOTNET8_ROOT="${DOTNET8_OUT##*$'\n'}"      # its contract: the root is the LAST line
        DOTNET8_NOTE="ensure_dotnet8.sh: .NET 8 at $DOTNET8_ROOT"; DOTNET8_STATE=ok
      else
        # its note joins NOTES after the resolution reader is known (section 3, ready_for.dia_raw)
        DOTNET8_NOTE="ensure_dotnet8.sh could not provide .NET 8 (it needs internet the first time; its own messages are on setup.sh's stderr)"; DOTNET8_STATE=failed
      fi ;;
    *) DOTNET8_NOTE="not run on $OS" ;;
  esac
fi

# ---- 2d. google-genai: the podcast's Gemini voices ---------------------------
# A Core run makes the podcast by default (SKILL.md step 9e), and on HIVE its job runs with this
# env's python (make_podcast.py job reads setup.json's `python`). Gemini 3 TTS needs
# client.interactions: google-genai >= GENAI_MIN (make_podcast.GENAI_MIN;
# tests/test_make_podcast_hive.py keeps the two equal). conda-forge's noarch build, in its own
# install like 2b's: only the podcast needs it, so a failure is a note. Installed = its
# conda-meta record (2b's rule); make_podcast.py job asks the python itself before it writes a job.
GENAI_MIN="2.25.0"
genai_version() {   # the version in google-genai's conda-meta record, or nothing
  local f
  for f in "$ENV_PREFIX"/conda-meta/google-genai-[0-9]*.json; do
    [ -e "$f" ] || continue
    f="${f##*/google-genai-}"; printf '%s' "${f%%-*}"; return 0
  done
  return 1
}
version_ge() {      # version_ge A B: is A >= B, on the first three numbers (bash 3.2: no sort -V)
  local a="$1." b="$2." i x y
  for i in 1 2 3; do
    x="${a%%.*}"; y="${b%%.*}"; a="${a#*.}"; b="${b#*.}"
    x="${x//[!0-9]/}"; y="${y//[!0-9]/}"
    [ "${x:-0}" -gt "${y:-0}" ] && return 0
    [ "${x:-0}" -lt "${y:-0}" ] && return 1
  done
  return 0
}
GENAI_VERSION="$(genai_version)"
genai_ok() { [ -n "$GENAI_VERSION" ] && version_ge "$GENAI_VERSION" "$GENAI_MIN"; }
# --freeze-installed, as for nanoparquet (2b): the env is shared by every analysis on the account,
# and the podcast's package must never move a package another analysis runs on (a dry run on
# HIVE, 2026-10: 41 packages added, none changed)
case "$CONDA" in
  *micromamba) GENAI_FIX="$CONDA install -y -r $MAMBA_ROOT -p $ENV_PREFIX --freeze-installed --override-channels -c conda-forge -c bioconda 'google-genai>=$GENAI_MIN'" ;;
  "")          GENAI_FIX="run setup.sh (it installs micromamba), or install conda-forge's google-genai>=$GENAI_MIN into $ENV_PREFIX with --freeze-installed" ;;
  *)           GENAI_FIX="$CONDA install -y -p $ENV_PREFIX --freeze-installed --override-channels -c conda-forge -c bioconda 'google-genai>=$GENAI_MIN'" ;;
esac
if [ -n "$CONDA" ] && env_ready && ! genai_ok && ! $CHECK_ONLY; then
  say "[setup] installing google-genai >= $GENAI_MIN (the podcast's Gemini voices) into the env..."
  case "$CONDA" in      # -p for the reason create_env gives
    *micromamba) in_prefix "$CONDA" install -y -r "$MAMBA_ROOT" -p "$ENV_PREFIX" --freeze-installed \
                   --override-channels -c conda-forge -c bioconda "google-genai>=$GENAI_MIN" >&2 ;;
    *)           in_prefix "$CONDA" install -y -p "$ENV_PREFIX" --freeze-installed \
                   --override-channels -c conda-forge -c bioconda "google-genai>=$GENAI_MIN" >&2 ;;
  esac || say "[setup] the google-genai install failed (above); only the podcast needs it"
  GENAI_VERSION="$(genai_version)"
fi
GENAI_OK=false; genai_ok && GENAI_OK=true
if env_ready && ! $GENAI_OK; then
  NOTES+=("google-genai ${GENAI_VERSION:-is not installed}${GENAI_VERSION:+ is older than $GENAI_MIN}: the podcast's Gemini voices (step 9e; on HIVE, make_podcast.py job) need google-genai >= $GENAI_MIN. Nothing else needs it. Fix: $GENAI_FIX.")
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
# The DE can run when R can read report.parquet the way THIS env's limpa does. run_de.R filters
# the report with arrow on both methods; limpa >= 1.4.0 then reads it with nanoparquet, limpa
# 1.2.x with arrow (readDIANN: arrow::read_parquet, github.com/bioc/limpa RELEASE_3_22). On
# 2026-10-16 an env on limpa 1.2.5 without nanoparquet ran every DE of the week (de_provenance.json:
# limpa 1.2.5) while setup.json said ready_for.de false, "no DE can run". limpa.ok and the exit
# status still hold the env to the pinned stack; this says only whether a DE can run.
# --method maxlfq needs no limpa at all: arrow and limma (de_reader.maxlfq).
r_read_ok() {
  "$ENV_PREFIX/bin/Rscript" -e "q(status = as.integer(!all(vapply(c('arrow', 'limma'), requireNamespace, logical(1), quietly = TRUE))))" >/dev/null 2>&1
}
DE_READER="$LIMPA_READS"; DE_READER_OK=false; R_READ_OK=false
if env_ready; then
  r_read_ok && R_READ_OK=true
  case "$DE_READER" in
    nanoparquet) $R_READ_OK && $NANOPARQUET_OK && DE_READER_OK=true ;;
    arrow)       $R_READ_OK && DE_READER_OK=true ;;
  esac
  $R_READ_OK || NOTES+=("arrow cannot be loaded in $ENV_PREFIX (or limma): run_de.R reads and filters report.parquet with arrow on both methods, so no DE can run. Fix: re-run bash setup.sh (on a cluster it runs itself on a compute node).")
fi
DE_READY=false; [ -n "$RSCRIPT" ] && $DE_READER_OK && DE_READY=true
DDA_READY=false
if [ -n "$SAGE" ] && [ -n "$RSCRIPT" ]; then DDA_READY=true; fi
DIA_READY=false
if $DIANN_READY && [ -n "$RSCRIPT" ]; then DIA_READY=true; fi
# A search can run (an Rscript outside the env counts for it) while the DE cannot: say why, or
# ready_for.de is false beside dia/dda true with no note (2.12 review). With the env built, the
# nanoparquet note above says it.
if ! $DE_READY && ! env_ready && { $DIA_READY || $DDA_READY; }; then
  NOTES+=("ready_for.de is false although a search can run here: the DE runs in the skill's own env ($ENV_PREFIX), which is not built yet, so whether its R can read report.parquet (nanoparquet) is not known (the Rscript found, $RSCRIPT, is not that env's). Run bash setup.sh (without --check) to build it before the DE step.")
fi

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
THERMO_UNASKED="Could not ask detect_acquisition.py whether Thermo .raw can be read (python: ${PY:-none found}). Re-run setup.sh to build the env, then setup.sh --check."
$THERMO_READY || NOTES+=("Thermo .raw cannot be read yet (only matters for .raw input): see thermo_raw_reader.note in setup.json for the exact fix.")

[ -z "$RSCRIPT" ] && NOTES+=("R/Rscript not available — DE cannot run. Re-run setup.sh to install it into the conda env.")
[ -z "$SAGE" ]    && NOTES+=("Sage not found — DDA search unavailable until the conda env is built.")
MSCONVERT_NOTE="msconvert is Linux-only on bioconda (and has no vendor readers). Sage still searches here: Thermo .raw through ThermoRawFileParser, mzML as it is, and a timsTOF .d read directly for identifications only (estimate_params.py --ids-only). To QUANTIFY a .d with Sage, run the search on HIVE, where the job converts it with the Core's pinned ProteoWizard container, or convert it to mzML elsewhere with a vendor-enabled msconvert (e.g. Windows ProteoWizard) first. DIA-NN reads .d/.raw natively."
[ "$OS" = "darwin" ] && [ -z "$MSCONVERT" ] && NOTES+=("$MSCONVERT_NOTE")

# ---- 3b. next_steps: the fix for each thing that is not ready, for THIS platform ------------
# setup.json next_steps {key: one sentence} and a "[setup] next:" line each at the end, so the
# agent relays a fix instead of deriving one from ready_for. Keys de, limpa, dia, dia_raw, dda, thermo_raw,
# resolution_reader; a ready item has none, so all ready is {}. Each sentence says when it matters
# and reuses the fix the notes above give (LIMPA_FIX, NANOPARQUET_FIX, DIANN_NOTE, MSCONVERT_NOTE,
# thermo_raw_reader's own notes) -- never a second wording of it. "A cluster" is setup.sh's own
# test (section 0): in a SLURM job, or SLURM / HPCCF's cluster file here.
ON_CLUSTER=false
if [ -n "${SLURM_JOB_ID:-}" ] || slurm_here; then ON_CLUSTER=true; fi
if $ON_CLUSTER; then
  RERUN_SETUP="bash $SCRIPT_DIR/setup.sh (without --check; from a login node it runs itself on a compute node in one blocking srun, so run it in the background and wait for it, and never install with conda, R or .NET on a login node by hand)"
  BY_HAND="by hand inside a job, never on a login node"
else
  RERUN_SETUP="bash $SCRIPT_DIR/setup.sh (without --check)"
  BY_HAND="by hand"
fi
# thermo_raw_reader's note and resolution_reader's, one line each (whitespace folded)
THERMO_NOTE="$THERMO_UNASKED"; RES_READY=""; RES_NOTE=""
if [ -n "$THERMO_READER" ]; then
  THERMO_FIELDS="$("$PY" -I -c 'import json, sys
d = json.loads(sys.argv[1]); r = d.get("resolution_reader") or {}
ready = {True: "true", False: "false"}.get(r.get("ready"), "")
for v in (d.get("note"), ready, r.get("note")):
    print(" ".join(str(v or "").split()))' "$THERMO_READER" 2>/dev/null)"
  { IFS= read -r THERMO_NOTE; IFS= read -r RES_READY; IFS= read -r RES_NOTE; } <<THERMO_EOF
$THERMO_FIELDS
THERMO_EOF
  [ -n "$THERMO_NOTE" ] || THERMO_NOTE="thermo_raw_reader in setup.json says why"
fi
# ready_for.dia_raw: can DIA-NN read Thermo .raw HERE -- ONE definition, from what this file
# records (2.15.1, HIVE F7: ready_for.dia was true beside a note that DIA-NN could not read .raw).
# true: DIA-NN can run and ensure_dotnet8.sh gave a .NET 8 root (dotnet8.root); false: DIA-NN
# cannot run, or ensure_dotnet8.sh failed; null: not checked (--check runs no ensure_dotnet8.sh;
# on macOS DIA-NN runs in Docker, whose image's own .NET is not checked here). ready_for.dia
# stays "a DIA search can run" (.d and mzML need no .NET), and a .raw search runs
# ensure_dotnet8.sh again itself (run_search.py, diann_parallel.py).
if ! $DIA_READY; then DIA_RAW=false
elif [ "$OS" = "darwin" ]; then DIA_RAW=null
elif [ "$DOTNET8_STATE" = ok ]; then DIA_RAW=true
elif [ "$DOTNET8_STATE" = failed ]; then DIA_RAW=false
else DIA_RAW=null
fi
if [ "$DOTNET8_STATE" = failed ]; then
  # what cannot run without it, as this file records it: the resolution reader may still load a
  # runtime it finds (thermo_raw_reader.resolution_reader.ready), and then it is not listed
  _no_dotnet="DIA-NN cannot read .raw here (ready_for.dia_raw false; a .raw search runs ensure_dotnet8.sh again itself)"
  [ "$RES_READY" = true ] || _no_dotnet="$_no_dotnet, and the Orbitrap resolution cannot be read"
  NOTES+=("$DOTNET8_NOTE. Until it does, $_no_dotnet: re-run bash setup.sh, or \`bash $ENSURE_DOTNET8\` (on a cluster both run themselves on a compute node, which on HIVE reaches Microsoft's downloads).")
  unset _no_dotnet
fi
NEXT_STEPS=""       # "key|sentence" lines (set_ann_status_json's form), in the order above
next_step() { NEXT_STEPS="$NEXT_STEPS$1|$2
"; }
MAXLFQ_TOO=""; $R_READ_OK && MAXLFQ_TOO="; --method maxlfq can still run meanwhile (arrow + limma, no limpa)"
if ! $DE_READY; then
  if ! env_ready; then
    next_step de "No DE can run until the skill's own env ($ENV_PREFIX) is built: run $RERUN_SETUP."
  elif ! $R_READ_OK; then
    next_step de "No DE can run: arrow (or limma) cannot be loaded in $ENV_PREFIX, and run_de.R reads report.parquet with arrow on both methods; run $RERUN_SETUP."
  elif [ "$DE_READER" = nanoparquet ]; then
    next_step de "nanoparquet cannot be loaded in $ENV_PREFIX, and $NANOPARQUET_COST with the default method (dpc); run $RERUN_SETUP, which installs it, or $BY_HAND: $NANOPARQUET_FIX$MAXLFQ_TOO."
  else
    next_step de "No DE with the default method (dpc) can run: $LIMPA_STATE in $ENV_PREFIX; fix: $LIMPA_FIX$MAXLFQ_TOO."
  fi
fi
if ! $LIMPA_OK; then
  if [ -n "$LIMPA_VERSION" ]; then      # an upgrade: other analyses may be running on this limpa
    next_step limpa "$LIMPA_STATE: run_de.R --method dpc needs limpa >= $LIMPA_MIN; fix: $LIMPA_FIX -- but upgrading an env other analyses run on changes their limpa mid-analysis, so Core staff ask the Core's director first."
  elif env_ready; then
    next_step limpa "$LIMPA_STATE: run_de.R --method dpc needs limpa >= $LIMPA_MIN; fix: $LIMPA_FIX."
  else
    next_step limpa "limpa >= $LIMPA_MIN comes with the skill's env ($ENV_PREFIX), which is not built yet: run $RERUN_SETUP."
  fi
fi
if ! $DIA_READY; then
  if $DIANN_READY; then
    next_step dia "Only for DIA data: DIA-NN can run here, but R is not installed yet (${DIANN_NOTE}); run $RERUN_SETUP."
  elif [ "$OS" = "darwin" ]; then
    next_step dia "Only for DIA data, and never by falling back to a DDA engine: run bash $SCRIPT_DIR/build_diann_docker.sh and relay what it prints to the user (install Docker Desktop, open it once) -- $DIANN_NOTE"
  else
    next_step dia "Only for DIA data: setup.sh knows no DIA-NN route on $OS; references/install.md says how to run DIA-NN there (Windows: its native build, from Git Bash or MSYS2)."
  fi
fi
if ! $DDA_READY; then
  DDA_MISSING=""
  [ -z "$SAGE" ] && DDA_MISSING="Sage"
  [ -z "$RSCRIPT" ] && DDA_MISSING="$DDA_MISSING${DDA_MISSING:+ and }R"
  DDA_MAC=""; [ "$OS" = "darwin" ] && [ -z "$MSCONVERT" ] && DDA_MAC=" Even then, $MSCONVERT_NOTE"
  next_step dda "Only for DDA data (DIA data goes to DIA-NN instead): $DDA_MISSING not found here; run $RERUN_SETUP.$DDA_MAC"
fi
if [ "$DIA_RAW" = false ] && $DIA_READY; then
  next_step dia_raw "Only for Thermo .raw input to DIA-NN (not for .d or mzML): no .NET 8 runtime >= 8.0.17 here (dotnet8.note); run $RERUN_SETUP, or bash $ENSURE_DOTNET8."
fi
if ! $THERMO_READY; then
  THERMO_THEN="; then bash $SCRIPT_DIR/setup.sh --check until ready_for.thermo_raw is true."
  [ -n "$THERMO_READER" ] || THERMO_THEN="."      # THERMO_UNASKED says that itself
  next_step thermo_raw "Only if the input has Thermo .raw files (not for .d or mzML), before step 2: ${THERMO_NOTE%.}$THERMO_THEN"
fi
if [ "$RES_READY" = false ]; then
  next_step resolution_reader "Only for Thermo Orbitrap .raw input: step 2 cannot read the Orbitrap resolution here, so it will ask the user for it; to read it instead: ${RES_NOTE%.}."
fi

# ---- 4. write activate.sh + setup.json --------------------------------------
if ! $CHECK_ONLY || [ ! -f "$ACTIVATE" ]; then
  cat > "$ACTIVATE" <<EOF
# source this to put the proteomics-pipeline environment on PATH
export PROTEOMICS_PIPELINE_HOME="$PP_HOME"
export PATH="$ENV_PREFIX/bin:\$PATH"
[ -f "$PP_HOME/diann_docker_image" ] && export DIANN_DOCKER_IMAGE="\$(cat "$PP_HOME/diann_docker_image")"
# .NET 8 from ensure_dotnet8.sh (NETCore >= 8.0.17 + AspNetCore): DIA-NN 2.6 reads .raw with
# it, and so does a framework-dependent ThermoRawFileParser. A DOTNET_ROOT already set is kept.
# The root: \$PROTEOMICS_DOTNET_DIR when set; else the first that exists of the one
# ensure_dotnet8.sh gave this setup (setup.json dotnet8.root), this prefix's dotnet8 and the old
# default ~/.proteomics-pipeline/dotnet8 (ensure_dotnet8.sh's own order).
# (a plain word list, no IFS splitting: activate.sh is also sourced from zsh)
if [ -z "\${DOTNET_ROOT:-}" ]; then
  for _pp_dotnet in "\${PROTEOMICS_DOTNET_DIR:-}" "$DOTNET8_ROOT" "$PP_HOME/dotnet8" "\$HOME/.proteomics-pipeline/dotnet8"; do
    if [ -n "\$_pp_dotnet" ] && [ -x "\$_pp_dotnet/dotnet" ]; then
      export DOTNET_ROOT="\$_pp_dotnet"; export PATH="\$DOTNET_ROOT:\$PATH"; break
    fi
    if [ -n "\${PROTEOMICS_DOTNET_DIR:-}" ]; then break; fi       # set: that root only
  done
fi
unset _pp_dotnet
EOF
fi

j() { printf '%s' "$1" | sed 's/\\/\\\\/g; s/"/\\"/g'; }
next_steps_json() {   # NEXT_STEPS (section 3b) as JSON members
  local line sep=""
  while IFS= read -r line; do
    [ -n "$line" ] || continue
    printf '%s"%s": "%s"' "$sep" "${line%%|*}" "$(j "${line#*|}")"; sep=", "
  done <<NEXT_EOF
$NEXT_STEPS
NEXT_EOF
}
# Which skill's setup.sh wrote this, and when: skill_version.sh --check-hive reads them on HIVE,
# so a setup.json older than the skill there is visible (2026-10: one written for 2.12 said
# ready_for.de true where 2.13's limpa and nanoparquet were missing)
SETUP_SKILL_VERSION="$( (. "$SCRIPT_DIR/skill_version.sh" && skill_version) 2>/dev/null || true)"
{
  printf '{\n'
  printf '  "skill_version": "%s", "written": "%s", "check_only": %s, "run_id": "%s",\n' \
         "$(j "$SETUP_SKILL_VERSION")" "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" "$($CHECK_ONLY && echo true || echo false)" \
         "$(j "${PROTEOMICS_SETUP_RUN_ID:-}")"
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
  printf '  "podcast": {"google_genai": "%s", "required": ">= %s", "ok": %s},\n' \
         "$(j "$GENAI_VERSION")" "$GENAI_MIN" "$($GENAI_OK && echo true || echo false)"
  printf '  "limpa": {"version": "%s", "required": ">= %s", "ok": %s, "source": "Bioconductor %s"},\n' \
         "$(j "$LIMPA_VERSION")" "$LIMPA_MIN" "$($LIMPA_OK && echo true || echo false)" "$LIMPA_BIOC"
  printf '  "de_reader": {"limpa": "%s", "reads_parquet_with": "%s", "ok": %s, "maxlfq": %s},\n' \
         "$(j "$LIMPA_VERSION")" "$DE_READER" "$($DE_READER_OK && echo true || echo false)" \
         "$($R_READ_OK && echo true || echo false)"
  printf '  "nanoparquet": {"ok": %s, "fix": "%s"},\n' \
         "$($NANOPARQUET_OK && echo true || echo false)" "$(j "$($NANOPARQUET_OK || printf '%s' "$NANOPARQUET_FIX")")"
  printf '  "set_annotation": {"ok": %s, "packages": "%s", "missing": "%s", "status": {%s}, "organisms": {"human": %s, "mouse": %s}, "bioconductor": "%s", "bioconductor_from": "%s", "fix": "%s"},\n' \
         "$($SET_ANN_OK && echo true || echo false)" "$SET_ANN_PKGS" "$(j "$SET_ANN_MISSING")" \
         "$(set_ann_status_json)" \
         "$(set_ann_loads AnnotationDbi && set_ann_loads GO.db && set_ann_loads org.Hs.eg.db && echo true || echo false)" \
         "$(set_ann_loads AnnotationDbi && set_ann_loads GO.db && set_ann_loads org.Mm.eg.db && echo true || echo false)" \
         "$(j "$SET_ANN_BIOC")" "$(j "$ENV_BIOC_FROM")" "$(j "$SET_ANN_FIX")"
  printf '  "ready_for": {"de": %s, "dia": %s, "dia_raw": %s, "dda": %s, "thermo_raw": %s},\n' \
         "$($DE_READY && echo true || echo false)" \
         "$($DIA_READY && echo true || echo false)" \
         "$DIA_RAW" \
         "$($DDA_READY && echo true || echo false)" \
         "$($THERMO_READY && echo true || echo false)"
  if [ -n "$THERMO_READER" ]; then printf '  "thermo_raw_reader": %s,\n' "$THERMO_READER"
  else printf '  "thermo_raw_reader": {"ready": false, "note": "%s"},\n' "$(j "$THERMO_UNASKED")"
  fi
  printf '  "next_steps": {%s},\n' "$(next_steps_json)"
  printf '  "notes": ['
  for i in "${!NOTES[@]}"; do
    printf '%s"%s"' "$( [ "$i" -gt 0 ] && echo ', ' )" "$(j "${NOTES[$i]}")"
  done
  printf ']\n}\n'
} | tee "$SETUP_JSON"

say ""
say "[setup] activate with:  source $ACTIVATE"
say "[setup] readiness report written to $SETUP_JSON"

# The pinned R stack is part of a working setup: an env whose limpa is < 1.4.0, or whose R cannot
# load nanoparquet, is a setup failure, said plainly and with the fix -- not a note to be
# scrolled past. Both are said before setup.sh exits.
R_STACK_OK=true
if env_ready && ! $LIMPA_OK; then
  say ""
  say "[setup] ERROR: limpa ${LIMPA_VERSION:-(not installed)} in $ENV_PREFIX -- the skill needs limpa >= $LIMPA_MIN (Bioconductor $LIMPA_BIOC)."
  say "[setup]        Fix: $LIMPA_FIX"
  R_STACK_OK=false
fi
if env_ready && ! $NANOPARQUET_OK; then
  say ""
  say "[setup] ERROR: nanoparquet cannot be loaded in $ENV_PREFIX -- limpa >= $LIMPA_MIN reads report.parquet with it, and run_de.R its contaminant share."
  say "[setup]        Fix: $NANOPARQUET_FIX   ($NANOPARQUET_WHY)"
  R_STACK_OK=false
fi
# setup.json next_steps, one line each
if [ -n "$NEXT_STEPS" ]; then
  say ""
  while IFS= read -r NS_LINE; do
    [ -n "$NS_LINE" ] && say "[setup] next: ${NS_LINE%%|*}: ${NS_LINE#*|}"
  done <<NEXT_EOF
$NEXT_STEPS
NEXT_EOF
fi
if ! $R_STACK_OK && ! $CHECK_ONLY; then exit 1; fi
