#!/bin/bash
# =============================================================================
# DE-LIMP HPC Setup & Launch Script for HIVE (UC Davis)
# =============================================================================
#
# First-time setup:
#   bash hpc_setup.sh install
#
# Launch (interactive):
#   bash hpc_setup.sh run
#
# Update to latest version:
#   bash hpc_setup.sh update
#
# =============================================================================

set -e

# --- Configuration ---
CONTAINER_DIR="${HOME}/containers"
SIF_FILE="${CONTAINER_DIR}/de-limp.sif"
HF_IMAGE="docker://registry.hf.space/brettsp-de-limp-proteomics:latest"
R_USER_LIB="${HOME}/R/delimp-lib"   # Persistent R library for extra packages
PORT=7860
MEM="32GB"
CPUS=8
TIME="8:00:00"
ACCOUNT="genome-center-grp"
PARTITION="high"

# --- Colors ---
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

print_header() {
    echo ""
    echo -e "${BLUE}============================================${NC}"
    echo -e "${BLUE}  DE-LIMP Proteomics — HPC Launcher${NC}"
    echo -e "${BLUE}============================================${NC}"
    echo ""
}

# --- Install / Pull Container ---
cmd_install() {
    print_header
    echo -e "${GREEN}[1/4] Creating directories...${NC}"
    mkdir -p "${CONTAINER_DIR}"
    mkdir -p "${HOME}/data"
    mkdir -p "${HOME}/results"
    mkdir -p "${HOME}/logs"
    mkdir -p "${R_USER_LIB}"

    echo -e "${GREEN}[2/4] Requesting compute node for build...${NC}"
    echo "  (SIF conversion is CPU-intensive — should not run on head node)"
    echo ""

    srun --account=${ACCOUNT} --partition=${PARTITION} --time=1:00:00 --mem=16GB --cpus-per-task=4 --pty bash -c '
        SIF="$1"; IMAGE="$2"

        echo "[3/4] Loading Apptainer..."
        module load apptainer 2>/dev/null || module load singularity 2>/dev/null || {
            echo "ERROR: Neither apptainer nor singularity module found."
            exit 1
        }

        echo "[4/4] Pulling DE-LIMP container from Hugging Face..."
        echo "This may take 10-20 minutes on first pull."
        echo ""
        apptainer pull --force "${SIF}" "${IMAGE}"

        echo ""
        echo "Setup complete!"
        echo ""
        echo "Container: ${SIF}"
        echo "Size: $(du -h "${SIF}" | cut -f1)"
        echo ""
        echo "To launch DE-LIMP, run: bash hpc_setup.sh run"
    ' _ "${SIF_FILE}" "${HF_IMAGE}"
}

# --- Update Container ---
cmd_update() {
    print_header
    echo -e "${GREEN}Requesting compute node for build...${NC}"
    srun --account=${ACCOUNT} --partition=${PARTITION} --time=1:00:00 --mem=16GB --cpus-per-task=4 --pty bash -c '
        SIF="$1"; IMAGE="$2"

        echo "Updating DE-LIMP container..."
        module load apptainer 2>/dev/null || module load singularity 2>/dev/null
        apptainer pull --force "${SIF}" "${IMAGE}"
        echo ""
        echo "Update complete!"
        echo "Size: $(du -h "${SIF}" | cut -f1)"
    ' _ "${SIF_FILE}" "${HF_IMAGE}"
}

# --- Run Interactive ---
cmd_run() {
    print_header

    # Check container exists
    if [ ! -f "${SIF_FILE}" ]; then
        echo -e "${RED}Container not found at ${SIF_FILE}${NC}"
        echo -e "Run ${YELLOW}bash hpc_setup.sh install${NC} first."
        exit 1
    fi

    # Use srun to execute directly on compute node (not login node)
    SRUN_CMD="srun --account=${ACCOUNT} --partition=${PARTITION} --time=${TIME} --mem=${MEM} --cpus-per-task=${CPUS} --pty"

    echo -e "${GREEN}Requesting compute node...${NC}"
    echo "  Account: ${ACCOUNT} | Partition: ${PARTITION} | Time: ${TIME} | Memory: ${MEM} | CPUs: ${CPUS}"
    echo ""

    ${SRUN_CMD} bash -c '
        PORT="$1"; SIF="$2"; HIVE_USER="$3"; REPO_DIR="$4"; R_LIB="$5"
        NODE=$(hostname)

        echo ""
        echo "=========================================="
        echo "  DE-LIMP is starting on: ${NODE}"
        echo "=========================================="
        echo ""
        echo ">>> On your LOCAL computer, run this command:"
        echo ""
        echo "  ssh -L ${PORT}:${NODE}:${PORT} ${HIVE_USER}@hive.hpc.ucdavis.edu"
        echo ""
        echo ">>> Then open your browser to:"
        echo ""
        echo "  http://localhost:${PORT}"
        echo ""
        echo "=========================================="
        echo "Press Ctrl+C to stop DE-LIMP"
        echo ""

        module load apptainer 2>/dev/null || module load singularity 2>/dev/null

        # Bind-mount local code over container code so git pull takes effect
        # without rebuilding the container image.
        # R_LIBS_USER provides a writable library for packages missing from the image
        # (install with: bash hpc_setup.sh packages)
        apptainer exec \
            --env R_LIBS_USER="${R_LIB}" \
            --bind ${HOME}/data:/data \
            --bind ${HOME}/results:/results \
            --bind ${R_LIB}:${R_LIB} \
            --bind ${REPO_DIR}/app.R:/srv/shiny-server/app.R \
            --bind ${REPO_DIR}/R:/srv/shiny-server/R \
            "${SIF}" \
            R -e "shiny::runApp('"'"'/srv/shiny-server/'"'"', host='"'"'0.0.0.0'"'"', port=${PORT})"
    ' _ "${PORT}" "${SIF_FILE}" "${USER}" "$(pwd)" "${R_USER_LIB}"
}

# --- Install Extra R Packages (runs on login node which has internet) ---
cmd_packages() {
    print_header

    if [ ! -f "${SIF_FILE}" ]; then
        echo -e "${RED}Container not found. Run 'bash hpc_setup.sh install' first.${NC}"
        exit 1
    fi

    mkdir -p "${R_USER_LIB}"
    echo -e "${GREEN}Installing missing R packages into ${R_USER_LIB}...${NC}"
    echo "  (This runs on the login node which has internet access)"
    echo ""

    module load apptainer 2>/dev/null || module load singularity 2>/dev/null

    apptainer exec \
        --env R_LIBS_USER="${R_USER_LIB}" \
        --bind "${R_USER_LIB}:${R_USER_LIB}" \
        "${SIF_FILE}" \
        R --no-save -e "
            .libPaths(c('${R_USER_LIB}', .libPaths()))
            options(repos = c(CRAN = 'https://cloud.r-project.org'))
            if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')

            pkgs <- c('clusterProfiler', 'enrichplot', 'org.Hs.eg.db', 'org.Mm.eg.db',
                       'KSEAapp', 'ggseqlogo', 'MOFA2')
            missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

            if (length(missing) > 0) {
                message(paste0('Installing: ', paste(missing, collapse = ', ')))
                BiocManager::install(missing, lib = '${R_USER_LIB}', ask = FALSE, update = FALSE)
            } else {
                message('All packages already installed!')
            }
            message('Done. Library: ${R_USER_LIB}')
            message(paste0('Packages available: ', length(list.dirs('${R_USER_LIB}', recursive = FALSE))))
        "

    echo ""
    echo -e "${GREEN}Package installation complete!${NC}"
    echo "These packages will be available when you run: bash hpc_setup.sh run"
}

# --- Help ---
cmd_help() {
    print_header
    echo "Usage: bash hpc_setup.sh <command>"
    echo ""
    echo "Commands:"
    echo "  install    Pull container from Hugging Face and set up directories"
    echo "  update     Pull the latest container version"
    echo "  packages   Install extra R packages (GSEA, MOFA2) — run once"
    echo "  run        Request a compute node and launch DE-LIMP"
    echo "  help       Show this help message"
    echo ""
    echo "Configuration (edit the top of this script to change):"
    echo "  PORT=${PORT}  MEM=${MEM}  CPUS=${CPUS}  TIME=${TIME}"
    echo ""
    echo "Quick start:"
    echo "  1. bash hpc_setup.sh install"
    echo "  2. bash hpc_setup.sh packages    (one-time, installs GSEA etc.)"
    echo "  3. bash hpc_setup.sh run"
    echo "  4. On your laptop: ssh -L ${PORT}:<node>:${PORT} ${USER}@hive.hpc.ucdavis.edu"
    echo "  5. Open http://localhost:${PORT}"
}

# --- Main ---
case "${1:-help}" in
    install)  cmd_install ;;
    update)   cmd_update ;;
    packages) cmd_packages ;;
    run)      cmd_run ;;
    help)     cmd_help ;;
    *)
        echo -e "${RED}Unknown command: $1${NC}"
        cmd_help
        exit 1
        ;;
esac
