# DE-LIMP HPC Deployment Guide

This guide covers running DE-LIMP on the UC Davis HIVE HPC cluster. Two approaches are supported:

1. **Docker + SSH (Recommended for Windows)** — Run DE-LIMP locally in Docker, connect to HPC via SSH for DIA-NN searches. See [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md) for Docker setup, then use the SSH file browser and "Load from HPC" features.
2. **Apptainer on HPC (Alternative)** — Run DE-LIMP directly on a compute node via Apptainer container with SLURM proxy. This guide focuses on this approach.

The automated launcher scripts handle everything — container installation, code updates, job submission, SSH tunneling, and cleanup.

### What's Included

- **Differential Expression** — Upload DIA-NN `.parquet` results, run the limpa/limma pipeline, explore volcano plots, heatmaps, and tables
- **DIA-NN Search Integration** — Submit DIA-NN database searches to SLURM directly from the **New Search** tab. Results auto-load when complete.
- **NCBI Proteome Download** — Download FASTA databases from NCBI with automatic gene symbol mapping
- **MOFA2 Multi-Omics Integration** — Combine 2-6 data views for unsupervised factor analysis
- **Phosphoproteomics** — Site-level DE, KSEA kinase activity, motif analysis
- **GSEA** — GO (BP/MF/CC) and KEGG pathways with automatic organism detection
- **Run Comparator** — Cross-tool DE comparison (DE-LIMP vs Spectronaut/FragPipe)
- **AI Chat** — Google Gemini integration for data exploration (requires API key)
- **Chromatography QC** — TIC trace extraction and per-run diagnostic flagging
- **Contaminant Analysis** — Automatic contaminant detection, per-sample breakdown, keratin flagging
- **SSH File Browser** — Visual directory navigation for remote data on HPC

---

## Quick Start (Recommended)

### Windows Users — Two Options

#### Option A: Docker + SSH (Recommended)

Run DE-LIMP locally in Docker and connect to HIVE via SSH for DIA-NN searches. This is the easiest setup for Windows users.

1. **Install Docker Desktop** and follow [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md)
2. **Double-click `Launch_DE-LIMP_Docker.bat`** — it auto-detects your SSH key, starts the Docker container, and opens `http://localhost:3838`
3. The app auto-connects to HIVE via SSH on startup (green "HPC" badge appears)
4. Use the **SSH File Browser** to navigate remote data directories
5. Use **"Load from HPC"** to download and analyze completed search results

> **Shared PC support:** Multiple Windows users on the same PC are handled automatically. The Docker launcher detects each user's SSH key from `C:\Users\{username}\.ssh\`.

#### Option B: Apptainer on HPC

Run DE-LIMP directly on a HIVE compute node via Apptainer.

**Prerequisites:** OpenSSH (built into Windows 10+) and an SSH key for HIVE (see [SSH Key Setup](#ssh-key-setup) below).

1. **Download [`Launch_DE-LIMP.bat`](https://raw.githubusercontent.com/bsphinney/DE-LIMP/main/Launch_DE-LIMP.bat)** into a folder (e.g., your Desktop). That's the only file you need — it auto-downloads everything else.

2. **Place your SSH key** (`id_ed25519` or `id_rsa`) in the same folder, or in `C:\Users\YourName\.ssh\`.

3. **Double-click `Launch_DE-LIMP.bat`**. On first run it will:
   - Download the remaining scripts from GitHub automatically
   - Ask for your HIVE username (saved for next time)
   - Install the container on shared storage (~5 GB, takes 10-20 min the first time)
   - Install any missing R packages to shared storage
   - Submit a SLURM job and wait for a compute node
   - Open an SSH tunnel and launch your browser to `http://localhost:7860`

4. **Press Ctrl+C** in the terminal window to stop. The launcher automatically cancels the SLURM job and closes the SSH tunnel.

> **First launch takes ~20 minutes** (container download + package install). Subsequent launches take ~1-2 minutes.

### Mac / Linux Users — Terminal Launcher

For Apptainer on HPC (Docker + SSH is also available for Mac/Linux — see [MACOS_DOCKER_INSTALL.md](MACOS_DOCKER_INSTALL.md) or [LINUX_DOCKER_INSTALL.md](LINUX_DOCKER_INSTALL.md)):

```bash
# Download the launcher (one time) — it auto-downloads hpc_setup.sh on first run
curl -O https://raw.githubusercontent.com/bsphinney/DE-LIMP/main/launch_delimp.sh
chmod +x launch_delimp.sh

# Launch (does everything automatically)
bash launch_delimp.sh
```

The launcher performs the same 7 steps as the Windows Apptainer version. Press Ctrl+C to stop.

---

## What the Launcher Does (Step by Step)

| Step | What Happens |
|------|-------------|
| 1. Find SSH key | Checks the script directory, then `~/.ssh/` for `id_ed25519`, `id_rsa`, or `*.pem` |
| 2. Get username | Prompts once, saves to `.delimp_config` for future runs |
| 3. Check container | Verifies container exists on shared storage; installs if missing |
| 4. Sync repo | Clones or `git pull`s the DE-LIMP repo on HIVE (for live code updates). Code update detection banner shown in app if behind. |
| 5. Check packages | Installs missing R packages (GSEA, MOFA2, etc.) into shared storage |
| 6. Create user dirs | Creates per-user directories on shared storage (`users/{username}/logs/`, `users/{username}/jobs/`) |
| 7. Submit job | Submits an `sbatch` job with SLURM proxy, polls for the compute node hostname, opens an SSH tunnel |
| 8. Open browser | Opens `http://localhost:7860` and waits until you press Ctrl+C |

On exit (Ctrl+C), the launcher cancels the SLURM job, removes the node sentinel file, and kills the SSH tunnel.

---

## SSH Key Setup

HIVE uses the **HiPPO portal** ([hippo.ucdavis.edu](https://hippo.ucdavis.edu)) for SSH key management.

### 1. Generate an SSH key (if you don't have one)

**Mac/Linux:**
```bash
ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519
```

**Windows (PowerShell):**
```powershell
ssh-keygen -t ed25519 -f $env:USERPROFILE\.ssh\id_ed25519
```

When prompted for a passphrase, **press Enter for no passphrase**. The launcher scripts use non-interactive SSH and cannot prompt for a passphrase.

### 2. Upload your public key to HiPPO

1. Go to [hippo.ucdavis.edu](https://hippo.ucdavis.edu) and log in with your UC Davis CAS credentials
2. Select **HIVE** as your cluster
3. Upload your **public** key (`~/.ssh/id_ed25519.pub`) — you can paste the contents or upload the file
4. Wait a few minutes for the key to propagate

> **Note:** HiPPO accepts one key at a time. The same key works for both Farm and HIVE.

### 3. Test the connection

```bash
ssh -i ~/.ssh/id_ed25519 username@hive.hpc.ucdavis.edu
```

This should log you in without a password prompt. If it asks for a password, your key wasn't uploaded correctly — go back to HiPPO and re-upload.

For more details, see the [UC Davis HPC SSH key guide](https://hpc.ucdavis.edu/faq2/ssh-keypair).

---

## Manual HPC Usage (Without Launcher)

If you prefer to manage things yourself, or are on a non-HIVE cluster:

### Option 1: `hpc_setup.sh` Commands

SSH to HIVE and run these directly:

```bash
# First-time setup (pulls container, creates directories)
bash hpc_setup.sh install

# Install extra R packages (GSEA, MOFA2 — one-time)
bash hpc_setup.sh packages

# Launch interactively (requests compute node via srun)
bash hpc_setup.sh run
# Then on your local machine: ssh -L 7860:<node>:7860 username@hive.hpc.ucdavis.edu
# Open browser to: http://localhost:7860

# Update container to latest version
bash hpc_setup.sh update

# Clone/update the GitHub repo (for code overlay)
bash hpc_setup.sh repo
```

### Option 2: Fully Manual

```bash
# SSH to HIVE
ssh username@hive.hpc.ucdavis.edu

# Pull the container (first time only)
mkdir -p /quobyte/proteomics-grp/de-limp/containers
module load apptainer
apptainer pull /quobyte/proteomics-grp/de-limp/containers/de-limp.sif docker://registry.hf.space/brettsp-de-limp-proteomics:latest

# Request a compute node
salloc --account=genome-center-grp --partition=high --time=8:00:00 --mem=32GB --cpus-per-task=8

# Run DE-LIMP
module load apptainer
apptainer exec /quobyte/proteomics-grp/de-limp/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"

# In a separate terminal on your local machine, set up the tunnel:
ssh -L 7860:<compute-node>:7860 username@hive.hpc.ucdavis.edu

# Open browser to http://localhost:7860
```

---

## Configuration

### Resource Defaults

Edit the top of `hpc_setup.sh` to change defaults:

| Setting | Default | Notes |
|---------|---------|-------|
| `PORT` | 7860 | Change if port conflicts with another user |
| `MEM` | 32GB | Increase to 64GB+ for MOFA2 or large GSEA |
| `CPUS` | 8 | More CPUs = faster pipeline |
| `TIME` | 8:00:00 | Max wall time for the SLURM job |
| `ACCOUNT` | genome-center-grp | Your SLURM account (see note below) |
| `PARTITION` | high | SLURM partition |

> **Non-Genome Center users:** All HIVE users have free access via `publicgrp`, but the `high` partition is limited to **8 CPUs and 128 GB RAM per job**, and `low` partition jobs can be **preempted (killed) without warning**. If you belong to a sponsored group, use your group account (e.g., `mylab-grp`) instead. Check your accounts with: `sacctmgr show assoc user=$USER format="account%20,partition%20"`

### Directory Layout on HIVE

All DE-LIMP files live on shared storage (`/quobyte/proteomics-grp/de-limp/`) to avoid home directory quota limits. Shared resources are used by everyone; per-user directories prevent conflicts when multiple people run simultaneously:

```
/quobyte/proteomics-grp/de-limp/
├── containers/
│   └── de-limp.sif          # Apptainer container (~5 GB, shared)
├── DE-LIMP/                  # Git repo — code overlay (shared)
├── R/
│   └── delimp-lib/           # R packages (shared, installed once)
├── data/                     # Bind-mounted as /data in container
├── results/                  # Bind-mounted as /results in container
└── users/
    ├── brettsp/              # Per-user directories
    │   ├── logs/             #   SLURM job logs + node sentinel files
    │   └── jobs/             #   Generated sbatch scripts
    └── jsmith/
        ├── logs/
        └── jobs/
```

> **Note:** Only `hpc_setup.sh` (~17 KB) and SLURM proxy temp files use your home directory. The container, R packages, repo, and all data stay on shared storage. Multiple users can run DE-LIMP simultaneously without conflicts.

### How Code Updates Work

The launcher uses a **bind-mount overlay** pattern: the latest code from the GitHub repo on shared storage is bind-mounted over the container's `/srv/shiny-server/` files at runtime. This means:

- **Code updates don't require rebuilding the container** — just `git pull`
- The launcher runs `git pull --ff-only` automatically on every launch
- Container rebuilds are only needed when R package dependencies change
- The app detects when local code is behind the repo and shows an update notification banner

### SLURM Proxy

When running inside Apptainer, SLURM commands (`sbatch`, `squeue`, `sacct`, etc.) are not available. The launcher starts a SLURM proxy process outside the container that relays commands via temp files on shared storage. All 9 SLURM command paths are covered: `sbatch`, `squeue`, `scancel`, `sacct`, `sinfo`, `sacctmgr`, `scontrol`, `srun`, `sbatch --test-only`. The cluster monitor dashboard works inside containers via this proxy.

### Shared Storage

The container bind-mounts `/quobyte/proteomics-grp` for access to the shared proteomics group storage. DIA-NN search results, raw data, FASTA databases, and activity logs on this volume are accessible from within the app.

**Pre-staged FASTA files** are available at `/quobyte/proteomics-grp/de-limp/fasta/` — these appear as a dropdown option in the FASTA selector for quick access to commonly used organisms.

---

## Core Facility Mode

For proteomics core facility staff who need report generation and QC tracking:

```bash
# Set up the shared facility directory (one time)
bash hpc_setup.sh setup-facility /share/genome-center/delimp

# Edit staff.yml with your team's SSH/SLURM settings
nano /share/genome-center/delimp/staff.yml

# Launch with Core Facility mode enabled
# The launcher scripts set DELIMP_CORE_DIR automatically
```

Core Facility mode adds:
- Staff selector (auto-fills SSH/SLURM settings per team member)
- QC run tracking (SQLite database)
- Report generation (Quarto HTML reports)
- Job history dashboard

---

## Troubleshooting

### "Container not found"
Run `bash hpc_setup.sh install` on HIVE, or let the launcher do it automatically.

### "Port already in use"
Another user may be running DE-LIMP on the same port. Edit `PORT=7860` in `hpc_setup.sh` to a different number (e.g., 7861, 8080).

### Launcher times out waiting for compute node
The SLURM queue may be full. Check with `squeue -u $USER` on HIVE. Try the `publicgrp/low` partition if `genome-center-grp/high` is busy — but note that **low partition jobs can be killed and requeued** when high-priority jobs need resources, so you may lose your session without warning.

### "Permission denied" on SSH key (Windows)
Windows requires strict file permissions on SSH keys. The PowerShell launcher tries to fix this automatically. If it fails, run:
```powershell
icacls C:\path\to\id_ed25519 /inheritance:r /grant:r "$($env:USERNAME):(R)"
```

### "Load key ... invalid format"
The SSH key file is corrupted or is not actually a key (e.g., it was downloaded as an HTML page). Open it in a text editor — it should start with `-----BEGIN OPENSSH PRIVATE KEY-----`. If not, copy your real key from `C:\Users\YourName\.ssh\id_ed25519` or generate a new one.

### "Disk quota exceeded" during container pull
HIVE home directories have limited quota (~5-10 GB). The launcher stores everything on shared storage, but Apptainer's default cache (`~/.apptainer/cache/`) still uses your home dir. The launcher sets `APPTAINER_CACHEDIR` to redirect this, but if you see this error, clean up your home dir:
```bash
ssh username@hive.hpc.ucdavis.edu
# Check what's using space
du -sh ~/* 2>/dev/null | sort -rh | head -10
# Remove old DE-LIMP files (now on shared storage)
rm -rf ~/containers ~/R/delimp-lib ~/DE-LIMP
# Clear Apptainer cache
rm -rf ~/.apptainer/cache
```

### "srun: command not found"
SLURM commands aren't on the default PATH in non-login SSH sessions. The launcher uses `bash -l` (login shell) to load the environment. If you're running `hpc_setup.sh` manually, use `bash -l hpc_setup.sh install`.

### "remote fsetstat: Failure" (SCP error on Windows)
Newer Windows OpenSSH versions use the SFTP protocol by default, which some servers don't support. The launcher adds the `-O` flag to force legacy SCP protocol. If you're running SCP manually, add `-O`:
```powershell
scp -O -i ~/.ssh/id_ed25519 file.txt username@hive.hpc.ucdavis.edu:/path/
```

### PowerShell script won't run ("not digitally signed")
Don't run `.\launch_delimp.ps1` directly. Use `Launch_DE-LIMP.bat` instead — it runs PowerShell with `-ExecutionPolicy Bypass` automatically.

### File browser only shows "Data"
If running via Docker (`localhost:3838`), the container only sees bind-mounted paths. Switch to **Remote (SSH)** connection mode to browse files on HIVE. If running via Apptainer on HPC (`localhost:7860`), the shared storage is auto-detected and should appear as "Proteomics" in the file picker.

### App shows "vunknown" for version
The container is outdated. Update with:
```bash
bash hpc_setup.sh update
```

### App shows "Local" badge instead of "HPC"
You're running the Docker container on your local machine, not Apptainer on HIVE. Use `Launch_DE-LIMP.bat` to start the HPC version (green badge, port 7860). The Docker version (`update_docker.sh`, red badge, port 3838) runs locally without the cluster.

### Out of memory
Request more memory. MOFA2 and large GSEA analyses need 64GB+:
```bash
# Edit hpc_setup.sh and change MEM="64GB", then relaunch
```

### SLURM proxy errors
The container runs inside Apptainer where SLURM commands (`sbatch`, `squeue`) aren't available. The launcher starts a SLURM proxy process outside the container that relays all 9 SLURM command paths. If DIA-NN search submission fails, check that the proxy is running (`ps aux | grep delimp_slurm_proxy`). The cluster monitor dashboard also relies on this proxy.

---

## Updating DE-LIMP

| What Changed | How to Update |
|--------------|---------------|
| Code only (R files) | Just relaunch — the launcher runs `git pull` automatically |
| R package dependencies | `bash hpc_setup.sh packages` on HIVE |
| Base container (new system libs) | `bash hpc_setup.sh update` on HIVE |

---

**Last updated:** 2026-03-17
**DE-LIMP version:** v3.7.0
**Apptainer/Singularity:** Compatible with versions 3.0+
