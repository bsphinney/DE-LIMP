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

    # Auto-update code from GitHub (login node has internet)
    if [ -d "$(pwd)/.git" ]; then
        echo -e "${GREEN}Updating code from GitHub...${NC}"
        git pull --ff-only 2>&1 | tail -1
        echo ""
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

        # --- SLURM Proxy ---
        # The Shiny app runs inside an Apptainer container where sbatch/squeue/sacct
        # are not available. This proxy runs OUTSIDE the container, watches a shared
        # directory for command requests, executes them, and writes results back.
        PROXY_DIR="${HOME}/.delimp_slurm_proxy"
        mkdir -p "${PROXY_DIR}"
        rm -f "${PROXY_DIR}"/cmd_* "${PROXY_DIR}"/result_*

        (
            while true; do
                for cmd_file in "${PROXY_DIR}"/cmd_*; do
                    [ -f "$cmd_file" ] || continue
                    id="${cmd_file##*/cmd_}"
                    result_file="${PROXY_DIR}/result_${id}"
                    cmd=$(cat "$cmd_file")
                    output=$(eval "$cmd" 2>&1)
                    rc=$?
                    printf "%d\n%s\n" "$rc" "$output" > "$result_file"
                    rm -f "$cmd_file"
                done
                sleep 1
            done
        ) &
        PROXY_PID=$!
        echo "SLURM proxy started (PID ${PROXY_PID})"

        # Ensure proxy is cleaned up when container exits
        cleanup_proxy() {
            kill $PROXY_PID 2>/dev/null
            rm -rf "${PROXY_DIR}"
        }
        trap cleanup_proxy EXIT

        # Bind-mount local code over container code so git pull takes effect
        # without rebuilding the container image.
        # R_LIBS_USER provides a writable library for packages missing from the image
        # (install with: bash hpc_setup.sh packages)
        apptainer exec \
            --env R_LIBS_USER="${R_LIB}" \
            --env DELIMP_SLURM_PROXY="${PROXY_DIR}" \
            --bind ${HOME}/data:/data \
            --bind ${HOME}/results:/results \
            --bind ${R_LIB}:${R_LIB} \
            --bind /quobyte/proteomics-grp:/quobyte/proteomics-grp \
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

# --- Clone / Update Repo (for bind-mount code overlay) ---
cmd_repo() {
    local REPO_DIR="${HOME}/DE-LIMP"
    local REPO_URL="https://github.com/bsphinney/DE-LIMP.git"

    if [ -d "${REPO_DIR}/.git" ]; then
        echo -e "${GREEN}Updating DE-LIMP repo...${NC}"
        cd "${REPO_DIR}" && git pull --ff-only
        echo -e "${GREEN}Repo updated: ${REPO_DIR}${NC}"
    else
        echo -e "${GREEN}Cloning DE-LIMP repo...${NC}"
        git clone "${REPO_URL}" "${REPO_DIR}"
        echo -e "${GREEN}Repo cloned: ${REPO_DIR}${NC}"
    fi
}

# --- Submit App via sbatch (non-interactive, for launcher scripts) ---
cmd_sbatch() {
    local CORE_DIR="${2:-}"
    local REPO_DIR="${3:-${HOME}/DE-LIMP}"
    local LOG_DIR="${HOME}/logs"
    local JOB_DIR="${HOME}/jobs"
    local SBATCH_FILE="${JOB_DIR}/delimp_app.sbatch"

    # Ensure directories exist
    mkdir -p "${LOG_DIR}" "${JOB_DIR}" "${R_USER_LIB}"

    # Auto-update code from GitHub
    if [ -d "${REPO_DIR}/.git" ]; then
        (cd "${REPO_DIR}" && git pull --ff-only 2>&1 | tail -1)
    fi

    # Check container exists
    if [ ! -f "${SIF_FILE}" ]; then
        echo "ERROR: Container not found at ${SIF_FILE}" >&2
        echo "Run 'bash hpc_setup.sh install' first." >&2
        exit 1
    fi

    # Build optional Core Facility env
    local ENV_ARGS=""
    if [ -n "${CORE_DIR}" ]; then
        ENV_ARGS="--env DELIMP_CORE_DIR=${CORE_DIR}"
    fi

    # Build optional repo bind-mount
    local REPO_BINDS=""
    if [ -d "${REPO_DIR}/R" ]; then
        REPO_BINDS="--bind ${REPO_DIR}/app.R:/srv/shiny-server/app.R --bind ${REPO_DIR}/R:/srv/shiny-server/R"
    fi

    # Write sbatch script
    cat > "${SBATCH_FILE}" <<SBATCH_EOF
#!/bin/bash
#SBATCH --job-name=delimp
#SBATCH --account=${ACCOUNT}
#SBATCH --partition=${PARTITION}
#SBATCH --time=${TIME}
#SBATCH --mem=${MEM}
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --output=${LOG_DIR}/delimp_%j.out
#SBATCH --error=${LOG_DIR}/delimp_%j.err

# Write hostname sentinel for launcher to discover compute node
echo "\$(hostname)" > ${LOG_DIR}/delimp_node_\${SLURM_JOB_ID}.txt

module load apptainer 2>/dev/null || module load singularity 2>/dev/null

# Start SLURM proxy for container
PROXY_DIR="\${HOME}/.delimp_slurm_proxy"
mkdir -p "\${PROXY_DIR}"
rm -f "\${PROXY_DIR}"/cmd_* "\${PROXY_DIR}"/result_*
(
    while true; do
        for cmd_file in "\${PROXY_DIR}"/cmd_*; do
            [ -f "\$cmd_file" ] || continue
            id="\${cmd_file##*/cmd_}"
            result_file="\${PROXY_DIR}/result_\${id}"
            cmd=\$(cat "\$cmd_file")
            output=\$(eval "\$cmd" 2>&1)
            rc=\$?
            printf "%d\n%s\n" "\$rc" "\$output" > "\$result_file"
            rm -f "\$cmd_file"
        done
        sleep 1
    done
) &
PROXY_PID=\$!
trap "kill \$PROXY_PID 2>/dev/null; rm -rf \${PROXY_DIR}" EXIT

apptainer exec \\
    --env R_LIBS_USER="${R_USER_LIB}" \\
    --env DELIMP_SLURM_PROXY="\${PROXY_DIR}" \\
    ${ENV_ARGS} \\
    --bind ${HOME}/data:/data \\
    --bind ${HOME}/results:/results \\
    --bind ${R_USER_LIB}:${R_USER_LIB} \\
    --bind /quobyte/proteomics-grp:/quobyte/proteomics-grp \\
    ${REPO_BINDS} \\
    "${SIF_FILE}" \\
    R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=${PORT})"
SBATCH_EOF

    # Find sbatch — non-interactive SSH may not have full PATH
    local SBATCH_CMD="sbatch"
    if ! command -v sbatch &>/dev/null; then
        for p in /cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin/sbatch \
                 /usr/bin/sbatch /opt/slurm/bin/sbatch /usr/local/bin/sbatch; do
            if [ -x "$p" ]; then SBATCH_CMD="$p"; break; fi
        done
    fi

    # Submit and parse job ID
    local SUBMIT_OUTPUT
    SUBMIT_OUTPUT=$($SBATCH_CMD "${SBATCH_FILE}" 2>&1)
    local SUBMIT_RC=$?

    if [ ${SUBMIT_RC} -ne 0 ]; then
        echo "ERROR: sbatch submission failed: ${SUBMIT_OUTPUT}" >&2
        exit 1
    fi

    # sbatch outputs "Submitted batch job 12345678"
    local JOB_ID
    JOB_ID=$(echo "${SUBMIT_OUTPUT}" | grep -oE '[0-9]+$')

    if [ -z "${JOB_ID}" ]; then
        echo "ERROR: Could not parse job ID from: ${SUBMIT_OUTPUT}" >&2
        exit 1
    fi

    echo "JOBID:${JOB_ID}"
}

# --- Setup Core Facility Shared Directory ---
cmd_setup_facility() {
    local CF_DIR="${2:-/share/genome-center/delimp}"

    echo -e "${GREEN}Setting up Core Facility directory: ${CF_DIR}${NC}"

    mkdir -p "${CF_DIR}/reports"
    mkdir -p "${CF_DIR}/state"

    # Write template staff.yml if it doesn't exist
    if [ ! -f "${CF_DIR}/staff.yml" ]; then
        cat > "${CF_DIR}/staff.yml" <<'STAFF_EOF'
# DE-LIMP Core Facility Staff Configuration
# Each entry maps a staff member to their SSH and SLURM settings.

staff:
  - name: "Your Name"
    email: "you@ucdavis.edu"
    lab: "Genome Center"
    ssh_host: "hive.hpc.ucdavis.edu"
    ssh_user: "your_username"
    ssh_key: "~/.ssh/id_ed25519"
    slurm_account: "genome-center-grp"
    slurm_partition: "high"

# Instruments tracked for QC
instruments:
  - name: "timsTOF HT"
  - name: "Exploris 480"

# LC methods
lc_methods:
  - name: "30min-Evosep"
  - name: "60min-nanoLC"
  - name: "90min-nanoLC"
STAFF_EOF
        echo -e "${GREEN}Template staff.yml written to ${CF_DIR}/staff.yml${NC}"
        echo "  Edit this file to add your team members."
    else
        echo "  staff.yml already exists — skipping."
    fi

    # Set group-writable permissions
    chmod -R g+w "${CF_DIR}" 2>/dev/null || true

    echo -e "${GREEN}Core Facility setup complete!${NC}"
    echo "  Reports dir: ${CF_DIR}/reports/"
    echo "  State dir:   ${CF_DIR}/state/"
    echo "  Staff config: ${CF_DIR}/staff.yml"
}

# --- Help ---
cmd_help() {
    print_header
    echo "Usage: bash hpc_setup.sh <command>"
    echo ""
    echo "Commands:"
    echo "  install          Pull container from Hugging Face and set up directories"
    echo "  update           Pull the latest container version"
    echo "  packages         Install extra R packages (GSEA, MOFA2) — run once"
    echo "  run              Request a compute node and launch DE-LIMP (interactive)"
    echo "  sbatch           Submit DE-LIMP as a batch job (for launchers)"
    echo "  repo             Clone or update the DE-LIMP GitHub repo"
    echo "  setup-facility   Create Core Facility shared directory structure"
    echo "  help             Show this help message"
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
    install)          cmd_install ;;
    update)           cmd_update ;;
    packages)         cmd_packages ;;
    run)              cmd_run ;;
    sbatch)           cmd_sbatch "$@" ;;
    repo)             cmd_repo ;;
    setup-facility)   cmd_setup_facility "$@" ;;
    help)             cmd_help ;;
    *)
        echo -e "${RED}Unknown command: $1${NC}"
        cmd_help
        exit 1
        ;;
esac
