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

    salloc --account=${ACCOUNT} --partition=${PARTITION} --time=1:00:00 --mem=16GB --cpus-per-task=4 bash -c "
        echo -e '${GREEN}[3/4] Loading Apptainer...${NC}'
        module load apptainer 2>/dev/null || module load singularity 2>/dev/null || {
            echo -e '${RED}ERROR: Neither apptainer nor singularity module found.${NC}'
            exit 1
        }

        echo -e '${GREEN}[4/4] Pulling DE-LIMP container from Hugging Face...${NC}'
        echo 'This may take 10-20 minutes on first pull.'
        echo ''
        apptainer pull --force ${SIF_FILE} ${HF_IMAGE}

        echo ''
        echo -e '${GREEN}Setup complete!${NC}'
        echo ''
        echo 'Container: ${SIF_FILE}'
        echo \"Size: \$(du -h ${SIF_FILE} | cut -f1)\"
        echo ''
        echo -e 'To launch DE-LIMP, run: ${YELLOW}bash hpc_setup.sh run${NC}'
    "
}

# --- Update Container ---
cmd_update() {
    print_header
    echo -e "${GREEN}Requesting compute node for build...${NC}"
    salloc --account=${ACCOUNT} --partition=${PARTITION} --time=1:00:00 --mem=16GB --cpus-per-task=4 bash -c "
        echo -e '${GREEN}Updating DE-LIMP container...${NC}'
        module load apptainer 2>/dev/null || module load singularity 2>/dev/null
        apptainer pull --force ${SIF_FILE} ${HF_IMAGE}
        echo ''
        echo -e '${GREEN}Update complete!${NC}'
        echo \"Size: \$(du -h ${SIF_FILE} | cut -f1)\"
    "
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

    ${SALLOC_CMD} bash -c "
        echo ''
        echo -e '${GREEN}==========================================${NC}'
        echo -e '${GREEN}  DE-LIMP is starting on: \$(hostname)${NC}'
        echo -e '${GREEN}==========================================${NC}'
        echo ''
        echo -e '${YELLOW}>>> On your LOCAL computer, run this command:${NC}'
        echo ''
        echo -e '  ssh -L ${PORT}:\$(hostname):${PORT} ${USER}@hive.hpc.ucdavis.edu'
        echo ''
        echo -e '${YELLOW}>>> Then open your browser to:${NC}'
        echo ''
        echo -e '  http://localhost:${PORT}'
        echo ''
        echo -e '${GREEN}==========================================${NC}'
        echo 'Press Ctrl+C to stop DE-LIMP'
        echo ''

        module load apptainer 2>/dev/null || module load singularity 2>/dev/null
        apptainer exec \\
            --bind \${HOME}/data:/data \\
            --bind \${HOME}/results:/results \\
            ${SIF_FILE} \\
            R -e \"shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=${PORT})\"
    "
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
