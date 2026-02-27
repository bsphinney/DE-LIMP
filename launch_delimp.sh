#!/bin/bash
# =============================================================================
# DE-LIMP Employee Launcher (Mac/Linux)
# =============================================================================
# Place this script next to your SSH key and run:
#   bash launch_delimp.sh
#
# It will:
#   1. Find your SSH key
#   2. SSH to HIVE and install the container if needed
#   3. Clone/update the DE-LIMP repo
#   4. Submit the app via sbatch
#   5. Open an SSH tunnel and launch your browser
# =============================================================================

set -euo pipefail

# --- Configuration (edit these if needed) ---
HIVE_HOST="hive.hpc.ucdavis.edu"
PORT=7860
CORE_DIR="/share/genome-center/delimp"
ACCOUNT="genome-center-grp"
PARTITION="high"
CONFIG_FILE=".delimp_config"
SETUP_SCRIPT="hpc_setup.sh"
MAX_WAIT_NODE=600    # seconds to wait for compute node (10 min)
MAX_WAIT_APP=120     # seconds to wait for app to respond (2 min)

# --- Colors ---
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# --- State ---
TUNNEL_PID=""
SLURM_JOB_ID=""
SSH_KEY=""
HIVE_USER=""

# --- Cleanup on exit ---
cleanup() {
    echo ""
    echo -e "${YELLOW}Shutting down...${NC}"

    # Cancel SLURM job
    if [ -n "${SLURM_JOB_ID}" ] && [ -n "${SSH_KEY}" ]; then
        echo "  Cancelling SLURM job ${SLURM_JOB_ID}..."
        ssh -i "${SSH_KEY}" -o StrictHostKeyChecking=accept-new -o ConnectTimeout=5 \
            "${HIVE_USER}@${HIVE_HOST}" \
            "scancel ${SLURM_JOB_ID} 2>/dev/null; rm -f ~/logs/delimp_node_${SLURM_JOB_ID}.txt" \
            2>/dev/null || true
    fi

    # Kill SSH tunnel
    if [ -n "${TUNNEL_PID}" ]; then
        echo "  Closing SSH tunnel (PID ${TUNNEL_PID})..."
        kill "${TUNNEL_PID}" 2>/dev/null || true
    fi

    echo -e "${GREEN}Done. Goodbye!${NC}"
    exit 0
}
trap cleanup INT TERM

print_header() {
    echo ""
    echo -e "${BLUE}============================================${NC}"
    echo -e "${BLUE}  DE-LIMP Proteomics — Employee Launcher${NC}"
    echo -e "${BLUE}============================================${NC}"
    echo ""
}

# --- Step 1: Find SSH key ---
find_ssh_key() {
    echo -e "${GREEN}[1/7] Looking for SSH key...${NC}"

    # Check current directory first
    local SCRIPT_DIR
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    for candidate in "${SCRIPT_DIR}/id_ed25519" "${SCRIPT_DIR}/id_rsa" "${SCRIPT_DIR}"/*.pem; do
        if [ -f "${candidate}" ]; then
            SSH_KEY="${candidate}"
            echo "  Found: ${SSH_KEY}"
            return 0
        fi
    done

    # Fall back to ~/.ssh/
    for candidate in "${HOME}/.ssh/id_ed25519" "${HOME}/.ssh/id_rsa" "${HOME}/.ssh"/*.pem; do
        if [ -f "${candidate}" ]; then
            SSH_KEY="${candidate}"
            echo "  Found: ${SSH_KEY}"
            return 0
        fi
    done

    echo -e "${RED}No SSH key found!${NC}"
    echo "  Place your SSH key (id_ed25519, id_rsa, or *.pem) next to this script"
    echo "  or in ~/.ssh/"
    exit 1
}

# --- Step 2: Get HIVE username ---
get_username() {
    echo -e "${GREEN}[2/7] Getting HIVE username...${NC}"

    local SCRIPT_DIR
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    local CONFIG_PATH="${SCRIPT_DIR}/${CONFIG_FILE}"

    # Try loading from config
    if [ -f "${CONFIG_PATH}" ]; then
        HIVE_USER=$(grep '^HIVE_USER=' "${CONFIG_PATH}" | cut -d= -f2)
    fi

    if [ -n "${HIVE_USER}" ]; then
        echo "  Using saved username: ${HIVE_USER}"
    else
        read -rp "  Enter your HIVE username: " HIVE_USER
        if [ -z "${HIVE_USER}" ]; then
            echo -e "${RED}Username cannot be empty.${NC}"
            exit 1
        fi
        # Save for next time
        echo "HIVE_USER=${HIVE_USER}" > "${CONFIG_PATH}"
        echo "  Saved to ${CONFIG_PATH}"
    fi
}

# --- Helper: run command on HIVE ---
hive_ssh() {
    ssh -i "${SSH_KEY}" -o StrictHostKeyChecking=accept-new -o ConnectTimeout=10 \
        "${HIVE_USER}@${HIVE_HOST}" "$@"
}

# --- Step 3: Check / install container ---
check_container() {
    echo -e "${GREEN}[3/7] Checking container on HIVE...${NC}"

    local HAS_SIF
    HAS_SIF=$(hive_ssh "test -f ~/containers/de-limp.sif && echo yes || echo no")

    if [ "${HAS_SIF}" = "yes" ]; then
        echo "  Container found."
    else
        echo -e "${YELLOW}  Container not found. Installing...${NC}"
        echo "  (This downloads ~5 GB and may take 10-20 minutes)"
        echo ""

        # Copy hpc_setup.sh to HIVE
        local SCRIPT_DIR
        SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        if [ -f "${SCRIPT_DIR}/${SETUP_SCRIPT}" ]; then
            scp -i "${SSH_KEY}" -o StrictHostKeyChecking=accept-new \
                "${SCRIPT_DIR}/${SETUP_SCRIPT}" "${HIVE_USER}@${HIVE_HOST}:~/${SETUP_SCRIPT}"
        else
            echo -e "${RED}Cannot find ${SETUP_SCRIPT} next to this launcher.${NC}"
            exit 1
        fi

        # Run install (needs TTY for srun)
        ssh -t -i "${SSH_KEY}" -o StrictHostKeyChecking=accept-new \
            "${HIVE_USER}@${HIVE_HOST}" "bash ~/${SETUP_SCRIPT} install"

        echo -e "${GREEN}  Container installed!${NC}"
    fi
}

# --- Step 4: Clone/update repo ---
update_repo() {
    echo -e "${GREEN}[4/7] Syncing DE-LIMP repo on HIVE...${NC}"

    hive_ssh "
        if [ -d ~/DE-LIMP/.git ]; then
            cd ~/DE-LIMP && git pull --ff-only 2>&1 | tail -1
        else
            git clone https://github.com/bsphinney/DE-LIMP.git ~/DE-LIMP 2>&1 | tail -1
        fi
    "
}

# --- Step 5: Check R packages ---
check_packages() {
    echo -e "${GREEN}[5/7] Checking R packages...${NC}"

    local PKG_COUNT
    PKG_COUNT=$(hive_ssh "ls -1d ~/R/delimp-lib/*/ 2>/dev/null | wc -l" | tr -d '[:space:]')

    if [ "${PKG_COUNT:-0}" -lt 3 ]; then
        echo -e "${YELLOW}  Missing R packages. Installing (this may take a few minutes)...${NC}"

        # Copy setup script if not already there
        local SCRIPT_DIR
        SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        scp -i "${SSH_KEY}" -o StrictHostKeyChecking=accept-new \
            "${SCRIPT_DIR}/${SETUP_SCRIPT}" "${HIVE_USER}@${HIVE_HOST}:~/${SETUP_SCRIPT}" 2>/dev/null || true

        hive_ssh "bash ~/${SETUP_SCRIPT} packages"
        echo -e "${GREEN}  Packages installed!${NC}"
    else
        echo "  Found ${PKG_COUNT} packages — OK."
    fi
}

# --- Step 6: Submit via sbatch ---
submit_job() {
    echo -e "${GREEN}[6/7] Submitting DE-LIMP to compute node...${NC}"

    # Ensure hpc_setup.sh is on HIVE
    local SCRIPT_DIR
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    scp -i "${SSH_KEY}" -o StrictHostKeyChecking=accept-new \
        "${SCRIPT_DIR}/${SETUP_SCRIPT}" "${HIVE_USER}@${HIVE_HOST}:~/${SETUP_SCRIPT}" 2>/dev/null || true

    local SUBMIT_OUTPUT
    SUBMIT_OUTPUT=$(hive_ssh "bash ~/${SETUP_SCRIPT} sbatch '${CORE_DIR}' ~/DE-LIMP")

    # Parse JOBID:<number>
    SLURM_JOB_ID=$(echo "${SUBMIT_OUTPUT}" | grep '^JOBID:' | cut -d: -f2)

    if [ -z "${SLURM_JOB_ID}" ]; then
        echo -e "${RED}Failed to submit job. Output:${NC}"
        echo "${SUBMIT_OUTPUT}"
        exit 1
    fi

    echo "  Submitted SLURM job: ${SLURM_JOB_ID}"

    # Poll for compute node hostname
    echo "  Waiting for compute node allocation..."
    local ELAPSED=0
    local NODE=""

    while [ ${ELAPSED} -lt ${MAX_WAIT_NODE} ]; do
        NODE=$(hive_ssh "cat ~/logs/delimp_node_${SLURM_JOB_ID}.txt 2>/dev/null" | tr -d '[:space:]')
        if [ -n "${NODE}" ]; then
            break
        fi
        sleep 5
        ELAPSED=$((ELAPSED + 5))
        printf "  Waiting... (%ds / %ds)\r" ${ELAPSED} ${MAX_WAIT_NODE}
    done
    echo ""

    if [ -z "${NODE}" ]; then
        echo -e "${RED}Timed out waiting for compute node (${MAX_WAIT_NODE}s).${NC}"
        echo "  Check job status: ssh ${HIVE_USER}@${HIVE_HOST} squeue -u ${HIVE_USER}"
        cleanup
        exit 1
    fi

    echo -e "${GREEN}  Running on node: ${NODE}${NC}"

    # Open SSH tunnel
    echo "  Opening SSH tunnel (localhost:${PORT} -> ${NODE}:${PORT})..."
    ssh -i "${SSH_KEY}" -o StrictHostKeyChecking=accept-new \
        -N -f -L "${PORT}:${NODE}:${PORT}" "${HIVE_USER}@${HIVE_HOST}"
    TUNNEL_PID=$(lsof -ti :"${PORT}" -sTCP:LISTEN 2>/dev/null | head -1)

    # Poll until app responds
    echo "  Waiting for DE-LIMP to start..."
    ELAPSED=0
    while [ ${ELAPSED} -lt ${MAX_WAIT_APP} ]; do
        if curl -s -o /dev/null -w "%{http_code}" "http://localhost:${PORT}" 2>/dev/null | grep -q "200\|302"; then
            break
        fi
        sleep 3
        ELAPSED=$((ELAPSED + 3))
        printf "  Starting... (%ds / %ds)\r" ${ELAPSED} ${MAX_WAIT_APP}
    done
    echo ""
}

# --- Step 7: Open browser ---
open_browser() {
    echo -e "${GREEN}[7/7] Opening browser...${NC}"

    local URL="http://localhost:${PORT}"

    # macOS
    if command -v open &>/dev/null; then
        open "${URL}"
    # Linux with xdg-open
    elif command -v xdg-open &>/dev/null; then
        xdg-open "${URL}"
    else
        echo "  Open your browser to: ${URL}"
    fi

    echo ""
    echo -e "${BLUE}============================================${NC}"
    echo -e "${BLUE}  DE-LIMP is running!${NC}"
    echo -e "${BLUE}============================================${NC}"
    echo ""
    echo "  URL:       ${URL}"
    echo "  SLURM Job: ${SLURM_JOB_ID}"
    echo "  Tunnel:    localhost:${PORT}"
    echo ""
    echo "  Press Ctrl+C to stop."
    echo ""

    # Wait forever until Ctrl+C
    while true; do
        sleep 60
    done
}

# --- Main ---
print_header
find_ssh_key
get_username
check_container
update_repo
check_packages
submit_job
open_browser
