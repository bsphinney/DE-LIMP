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

    echo -e "${GREEN}[2/4] Requesting compute node for build...${NC}"
    echo "  (SIF conversion is CPU-intensive — should not run on head node)"
    echo ""

    salloc --account=${ACCOUNT} --partition=${PARTITION} --time=1:00:00 --mem=16GB --cpus-per-task=4 bash -c '
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
    salloc --account=${ACCOUNT} --partition=${PARTITION} --time=1:00:00 --mem=16GB --cpus-per-task=4 bash -c '
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

    # Build salloc command
    SALLOC_CMD="salloc --account=${ACCOUNT} --partition=${PARTITION} --time=${TIME} --mem=${MEM} --cpus-per-task=${CPUS}"

    echo -e "${GREEN}Requesting compute node...${NC}"
    echo "  Account: ${ACCOUNT} | Partition: ${PARTITION} | Time: ${TIME} | Memory: ${MEM} | CPUs: ${CPUS}"
    echo ""

    ${SALLOC_CMD} bash -c '
        PORT="$1"; SIF="$2"; HIVE_USER="$3"
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
        apptainer exec \
            --bind ${HOME}/data:/data \
            --bind ${HOME}/results:/results \
            "${SIF}" \
            R -e "shiny::runApp('"'"'/srv/shiny-server/'"'"', host='"'"'0.0.0.0'"'"', port=${PORT})"
    ' _ "${PORT}" "${SIF_FILE}" "${USER}"
}

# --- Help ---
cmd_help() {
    print_header
    echo "Usage: bash hpc_setup.sh <command>"
    echo ""
    echo "Commands:"
    echo "  install   Pull container from Hugging Face and set up directories"
    echo "  update    Pull the latest container version"
    echo "  run       Request a compute node and launch DE-LIMP"
    echo "  help      Show this help message"
    echo ""
    echo "Configuration (edit the top of this script to change):"
    echo "  PORT=${PORT}  MEM=${MEM}  CPUS=${CPUS}  TIME=${TIME}"
    echo ""
    echo "Quick start:"
    echo "  1. bash hpc_setup.sh install"
    echo "  2. bash hpc_setup.sh run"
    echo "  3. On your laptop: ssh -L ${PORT}:<node>:${PORT} ${USER}@hive.hpc.ucdavis.edu"
    echo "  4. Open http://localhost:${PORT}"
}

# --- Main ---
case "${1:-help}" in
    install) cmd_install ;;
    update)  cmd_update ;;
    run)     cmd_run ;;
    help)    cmd_help ;;
    *)
        echo -e "${RED}Unknown command: $1${NC}"
        cmd_help
        exit 1
        ;;
esac
