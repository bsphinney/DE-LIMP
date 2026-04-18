#!/bin/bash
# =============================================================================
# DE-LIMP WSL Setup & Launch Script (for Ubuntu under Windows WSL2)
# =============================================================================
#
# Usage (inside WSL Ubuntu):
#   bash delimp_wsl_setup.sh install    # One-time install (R + system deps + R packages)
#   bash delimp_wsl_setup.sh update     # git pull + re-check R packages
#   bash delimp_wsl_setup.sh run        # Launch the Shiny app on localhost:3838
#   bash delimp_wsl_setup.sh            # Auto: install if needed, then run
#
# Typically invoked by Launch_DE-LIMP_WSL.bat on the Windows side.
#
# Design notes:
#   - Clones the repo to ~/.delimp/DE-LIMP (native WSL filesystem — fast)
#   - Installs R packages to ~/.delimp/R-lib (native WSL filesystem)
#   - Raw data / FASTA live at /mnt/c/Users/<you>/DE-LIMP/data/ so Windows
#     File Explorer can still drop files in.
#   - Binds 0.0.0.0:3838 so Windows localhost:3838 reaches it via WSL2
#     port forwarding (enabled by default in Windows 10/11).
# =============================================================================

set -e

DELIMP_BASE="${HOME}/.delimp"
REPO_DIR="${DELIMP_BASE}/DE-LIMP"
R_LIB="${DELIMP_BASE}/R-lib"
DIANN_DIR="${DELIMP_BASE}/diann"
DIANN_LICENSE_FLAG="${DELIMP_BASE}/.diann_license_accepted"
DATA_DIR="${DELIMP_DATA_DIR:-${DELIMP_BASE}/data}"
PORT="${DELIMP_PORT:-3838}"
REPO_URL="https://github.com/bsphinney/DE-LIMP.git"
# DIA-NN version. All 2.x releases live under the same GitHub tag ("2.0"),
# but the filename embeds the actual version. Default is pinned to 2.3.2
# (community-validated, matches what's installed on HIVE). Opt in to the
# latest release with DIANN_VERSION=latest, or pin an explicit version
# like DIANN_VERSION=2.5.0.
DIANN_VERSION="${DIANN_VERSION:-2.3.2}"
DIANN_RELEASE_TAG="2.0"   # Static tag used by all 2.x DIA-NN releases

# --- Colors ---
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log()  { echo -e "${BLUE}[delimp]${NC} $*"; }
warn() { echo -e "${YELLOW}[delimp]${NC} $*"; }
err()  { echo -e "${RED}[delimp]${NC} $*" >&2; }
ok()   { echo -e "${GREEN}[delimp]${NC} $*"; }

# -----------------------------------------------------------------------------
# 1. System dependencies (apt)
# -----------------------------------------------------------------------------
# Ubuntu 22.04 ships R 4.1, 24.04 ships R 4.3 — both too old for Bioc 3.22
# (which limpa needs). We add CRAN's Ubuntu repo to get R 4.5+.
install_system_deps() {
    log "Installing system dependencies via apt (may prompt for sudo password)..."

    # Basic tools needed to add the CRAN repo
    sudo apt-get update
    sudo apt-get install -y software-properties-common dirmngr gnupg lsb-release wget

    # Add CRAN repo for latest R (if not already present)
    if ! grep -rq "cloud.r-project.org/bin/linux/ubuntu" /etc/apt/sources.list.d/ 2>/dev/null; then
        log "Adding CRAN repo for latest R..."
        wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
            | sudo tee /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc >/dev/null
        local codename="$(lsb_release -cs)"
        echo "deb https://cloud.r-project.org/bin/linux/ubuntu ${codename}-cran40/" \
            | sudo tee /etc/apt/sources.list.d/cran.list >/dev/null
        sudo apt-get update
    fi

    sudo apt-get install -y \
        r-base r-base-dev \
        build-essential cmake pkg-config \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libcairo2-dev libxt-dev \
        openssh-client git unzip curl \
        python3 python3-pip python3-venv

    local rver=$(R --version 2>/dev/null | head -1)
    ok "System dependencies installed. R: ${rver}"
}

# -----------------------------------------------------------------------------
# 1b. DIA-NN + .NET runtime (for local searches inside WSL)
# -----------------------------------------------------------------------------
# DIA-NN 2.0 ships a Linux binary. Needs:
#   - .NET 8 runtime (to read Thermo .raw via RawFileReader)
#   - libgomp1, libstdc++6 (usually already installed by build-essential)
#   - unzip (to extract the release archive)
#
# License: DIA-NN is free for academic use but proprietary. Users must
# agree to terms at https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md
install_diann() {
    # License check — write a flag file once user agrees, skip on subsequent runs
    if [ ! -f "${DIANN_LICENSE_FLAG}" ]; then
        echo ""
        echo -e "${YELLOW}====================== DIA-NN License ======================${NC}"
        echo "  DIA-NN is developed by Vadim Demichev."
        echo "  Free for academic and non-commercial use."
        echo "  Commercial use requires a separate license from the author."
        echo ""
        echo "  Full terms: https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md"
        echo ""
        echo "  Citation: Demichev V, Messner CB, Vernardis SI, Lilley KS,"
        echo "  Ralser M. DIA-NN. Nature Methods. 2020;17(1):41-44."
        echo -e "${YELLOW}=============================================================${NC}"
        echo ""
        read -p "Do you accept the DIA-NN license terms? [yes/no]: " accept
        if [ "${accept}" != "yes" ]; then
            warn "License not accepted. Skipping DIA-NN install — WSL mode will be HPC-only."
            return 0
        fi
        mkdir -p "${DELIMP_BASE}"
        date > "${DIANN_LICENSE_FLAG}"
    fi

    # .NET 8 runtime — Ubuntu 24.04 has it in the default apt repos;
    # 22.04 needs Microsoft's apt repo.
    if ! command -v dotnet >/dev/null 2>&1; then
        log "Installing .NET 8 runtime..."
        if apt-cache show dotnet-runtime-8.0 >/dev/null 2>&1; then
            sudo apt-get install -y dotnet-runtime-8.0
        else
            log "dotnet-runtime-8.0 not in default repos — adding Microsoft repo..."
            local codename="$(lsb_release -cs)"
            wget -q "https://packages.microsoft.com/config/ubuntu/$(lsb_release -rs)/packages-microsoft-prod.deb" \
                -O /tmp/packages-microsoft-prod.deb
            sudo dpkg -i /tmp/packages-microsoft-prod.deb
            rm /tmp/packages-microsoft-prod.deb
            sudo apt-get update
            sudo apt-get install -y dotnet-runtime-8.0
        fi
    fi

    # Resolve the version to download. "latest" triggers an API lookup for the
    # newest non-Preview Linux zip; anything else is treated as an explicit
    # pin (e.g. "2.3.2", "2.5.0").
    if [ "${DIANN_VERSION}" = "latest" ]; then
        log "Querying GitHub for latest DIA-NN Linux release..."
        local resolved="$(curl -s --max-time 15 \
            "https://api.github.com/repos/vdemichev/DiaNN/releases/tags/${DIANN_RELEASE_TAG}" \
            | grep -oE '"name": "DIA-NN-[0-9.]+-Academia-Linux\.zip"' \
            | grep -v -i 'preview' \
            | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' \
            | sort -V | tail -1)"
        if [ -z "${resolved}" ]; then
            warn "GitHub API unreachable — falling back to pinned version 2.3.2"
            DIANN_VERSION="2.3.2"
        else
            DIANN_VERSION="${resolved}"
            log "Newest DIA-NN stable release: ${DIANN_VERSION}"
        fi
    fi

    # DIA-NN binary — download once, skip if already present
    if [ ! -x "${DIANN_DIR}/diann-linux" ]; then
        log "Downloading DIA-NN ${DIANN_VERSION} (~500 MB)..."
        mkdir -p "${DIANN_DIR}"
        local zipfile="${DIANN_DIR}/diann.zip"
        local url="https://github.com/vdemichev/DiaNN/releases/download/${DIANN_RELEASE_TAG}/DIA-NN-${DIANN_VERSION}-Academia-Linux.zip"
        wget --progress=bar -O "${zipfile}" "${url}"
        if [ ! -s "${zipfile}" ]; then
            err "DIA-NN download failed. Check version ${DIANN_VERSION} exists at https://github.com/vdemichev/DiaNN/releases"
            rm -f "${zipfile}"
            return 1
        fi
        log "Extracting..."
        unzip -q "${zipfile}" -d "${DIANN_DIR}/extract"
        # Binaries sometimes sit in a subdirectory — flatten
        local bin_src="$(find "${DIANN_DIR}/extract" -name diann-linux -type f | head -1)"
        if [ -z "${bin_src}" ]; then
            err "diann-linux not found in extracted archive."
            return 1
        fi
        local bin_dir="$(dirname "${bin_src}")"
        mv "${bin_dir}"/* "${DIANN_DIR}/"
        rm -rf "${DIANN_DIR}/extract" "${zipfile}"
        chmod +x "${DIANN_DIR}/diann-linux"
    fi

    # Create user bin symlink so `diann` is on PATH
    mkdir -p "${HOME}/.local/bin"
    ln -sf "${DIANN_DIR}/diann-linux" "${HOME}/.local/bin/diann"

    # Persistent LD_LIBRARY_PATH and PATH via ~/.delimp/env.sh (sourced by run_app)
    mkdir -p "${DELIMP_BASE}"
    cat > "${DELIMP_BASE}/env.sh" <<EOF
# Auto-generated by delimp_wsl_setup.sh — do not edit.
export LD_LIBRARY_PATH="${DIANN_DIR}:\${LD_LIBRARY_PATH}"
export PATH="\${HOME}/.local/bin:\${PATH}"
EOF

    # Sanity check
    if "${DIANN_DIR}/diann-linux" --help >/dev/null 2>&1; then
        ok "DIA-NN installed: $(${DIANN_DIR}/diann-linux --help 2>&1 | head -1)"
    else
        warn "DIA-NN installed but --help failed. Try running ${DIANN_DIR}/diann-linux --help manually to see the error."
    fi
}

# -----------------------------------------------------------------------------
# 2. Clone or update the repo
# -----------------------------------------------------------------------------
sync_repo() {
    if [ ! -d "${REPO_DIR}/.git" ]; then
        log "Cloning DE-LIMP into ${REPO_DIR}..."
        mkdir -p "${DELIMP_BASE}"
        git clone --depth 1 "${REPO_URL}" "${REPO_DIR}"
        ok "Repo cloned."
    else
        log "Updating DE-LIMP (git pull in ${REPO_DIR})..."
        git -C "${REPO_DIR}" pull --ff-only || warn "git pull failed — continuing with existing code."
    fi
}

# -----------------------------------------------------------------------------
# 3. R packages
# -----------------------------------------------------------------------------
install_r_packages() {
    mkdir -p "${R_LIB}"
    export R_LIBS_USER="${R_LIB}"

    log "Installing R packages into ${R_LIB} (first run: 20-30 min)..."

    R --no-save <<EOF
r_lib <- "${R_LIB}"
if (!dir.exists(r_lib)) dir.create(r_lib, recursive = TRUE)
.libPaths(c(r_lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"),
        Ncpus = max(1, parallel::detectCores() - 1))

cran <- c(
  "bslib", "readr", "tibble", "dplyr", "tidyr", "ggplot2", "httr2",
  "rhandsontable", "DT", "arrow", "shinyjs", "plotly", "stringr", "ggrepel",
  "remotes", "BiocManager", "markdown", "shinyFiles", "jsonlite", "processx",
  "callr", "KSEAapp", "ggseqlogo", "ggdendro", "systemfonts", "gdtools", "Rcpp",
  "ggraph", "graphlayouts", "tidygraph", "scatterpie", "shadowtext", "ggforce",
  "DBI", "RSQLite", "yaml", "uuid", "quarto", "shiny")

missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  cat("[delimp] Installing CRAN:", paste(missing_cran, collapse = ", "), "\n")
  # LIBARROW_MINIMAL=false gives arrow full codec support (zstd needed for parquet)
  Sys.setenv(LIBARROW_MINIMAL = "false")
  install.packages(missing_cran, lib = r_lib)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = r_lib)

bioc <- c(
  "DOSE", "GOSemSim", "yulab.utils",
  "limma", "limpa", "ComplexHeatmap", "AnnotationDbi",
  "org.Hs.eg.db", "org.Mm.eg.db", "ggridges",
  "ggtree", "ggtangle",
  "clusterProfiler", "enrichplot",
  "MOFA2", "basilisk")

missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  cat("[delimp] Installing Bioconductor:", paste(missing_bioc, collapse = ", "), "\n")
  BiocManager::install(missing_bioc, lib = r_lib, ask = FALSE, update = FALSE)
}

# Verify limpa loaded — most likely failure point
if (!requireNamespace("limpa", quietly = TRUE))
  stop("limpa failed to install. Check build logs above.")

cat("[delimp] All R packages OK.\n")
EOF

    ok "R packages installed."
}

# -----------------------------------------------------------------------------
# 4. Run the app
# -----------------------------------------------------------------------------
run_app() {
    if [ ! -f "${REPO_DIR}/app.R" ]; then
        err "app.R not found at ${REPO_DIR}. Run 'install' first."
        exit 1
    fi

    mkdir -p "${DATA_DIR}/raw" "${DATA_DIR}/fasta" "${DATA_DIR}/output" "${DATA_DIR}/ssh"

    export R_LIBS_USER="${R_LIB}"
    export DELIMP_DATA_DIR="${DATA_DIR}"

    # Source DIA-NN PATH/LD_LIBRARY_PATH if installed
    [ -f "${DELIMP_BASE}/env.sh" ] && . "${DELIMP_BASE}/env.sh"

    log "Starting DE-LIMP on http://localhost:${PORT}"
    log "  Repo:  ${REPO_DIR}"
    log "  Data:  ${DATA_DIR}"
    log "  R lib: ${R_LIB}"
    log ""
    log "Press Ctrl+C to stop."

    cd "${REPO_DIR}"
    exec R --no-save -e "shiny::runApp('.', host = '0.0.0.0', port = ${PORT}, launch.browser = FALSE)"
}

# -----------------------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------------------
CMD="${1:-auto}"
case "${CMD}" in
    install)
        install_system_deps
        sync_repo
        install_r_packages
        install_diann
        ok "Install complete. Run: bash $0 run"
        ;;
    update)
        sync_repo
        install_r_packages
        ok "Update complete."
        ;;
    run)
        run_app
        ;;
    diann)
        # Install DIA-NN only (e.g., after declining license on first run)
        install_diann
        ;;
    auto)
        # Install only if missing, then run
        if ! command -v R >/dev/null 2>&1; then
            install_system_deps
        fi
        if [ ! -d "${REPO_DIR}/.git" ]; then
            sync_repo
        fi
        if [ ! -d "${R_LIB}/limpa" ]; then
            install_r_packages
        fi
        if [ ! -x "${DIANN_DIR}/diann-linux" ] && [ ! -f "${DIANN_LICENSE_FLAG}" ]; then
            install_diann
        fi
        run_app
        ;;
    *)
        echo "Usage: bash $0 [install|update|run|diann]"
        echo "  install — install system deps, R packages, and DIA-NN"
        echo "  update  — git pull + refresh R packages"
        echo "  run     — launch the Shiny app on localhost:\${DELIMP_PORT:-3838}"
        echo "  diann   — install DIA-NN only (accepts license on first run)"
        exit 1
        ;;
esac
