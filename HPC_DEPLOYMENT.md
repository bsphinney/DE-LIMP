# DE-LIMP HPC Deployment Guide

This guide covers deploying DE-LIMP on HPC clusters using Apptainer/Singularity containers. While written with UC Davis HPC (FARM/HPC1/HPC2) in mind, these instructions work on any HPC system with Apptainer/Singularity.

### What's Included in v3.0

The DE-LIMP container includes everything you need:

- **Differential Expression** — Upload DIA-NN `.parquet` results, run the limpa/limma pipeline, explore volcano plots, heatmaps, and tables
- **DIA-NN Search Integration** — Submit DIA-NN database searches to your HPC cluster's SLURM scheduler directly from the app's **New Search** tab via SSH. Results auto-load when complete. (Requires DIA-NN installed on your cluster — see [DIA-NN Search on HPC](#-dia-nn-search-on-hpc) below)
- **MOFA2 Multi-Omics Integration** — Combine 2-6 data views (proteomics, phospho, transcriptomics, etc.) for unsupervised factor analysis. Includes example datasets.
- **Phosphoproteomics** — Site-level DE, KSEA kinase activity, motif analysis
- **GSEA** — GO (BP/MF/CC) and KEGG pathways with automatic organism detection
- **XIC Chromatogram Viewer** — Fragment-level validation with MS2 intensity alignment (see [XIC Viewing](#-xic-chromatogram-viewing-on-hpc) below)
- **AI Chat** — Google Gemini integration for data exploration (requires API key)

---

## 🎯 Quick Start (3 Options)

Choose the method that works best for you:

### **Option 1: Pull from Hugging Face** ⭐ **EASIEST - RECOMMENDED**
```bash
# SSH to your HPC cluster
ssh username@your-hpc-cluster.edu

# Load Apptainer (or Singularity)
module load apptainer  # or: module load singularity

# Create directory for containers
mkdir -p ~/containers

# Pull directly from Hugging Face (takes ~5-10 minutes)
apptainer pull ~/containers/de-limp.sif docker://huggingface.co/spaces/brettsp/de-limp-proteomics

# Done! Skip to Part 3 to run it.
```

### **Option 2: Build from GitHub Source**
```bash
# SSH to your HPC cluster
ssh username@your-hpc-cluster.edu

# Clone the repository
git clone https://github.com/bsphinney/DE-LIMP.git
cd DE-LIMP

# Load Apptainer and build (takes 30-45 minutes)
module load apptainer
mkdir -p ~/containers
apptainer build ~/containers/de-limp.sif Dockerfile

# Clean up build directory (optional)
cd ~ && rm -rf ~/DE-LIMP
```

### **Option 3: Build Locally, Transfer to Cluster**
If you have Docker on your local machine and want to build there:

```bash
# On your local machine:
# 1. Clone the repository
git clone https://github.com/bsphinney/DE-LIMP.git
cd DE-LIMP

# 2. Build Docker image (30-45 minutes)
docker build -t de-limp:latest .

# 3. Save as tar archive
docker save de-limp:latest -o de-limp-docker.tar

# 4. Transfer to cluster (adjust username and host)
scp de-limp-docker.tar username@your-hpc-cluster.edu:~/

# 5. SSH to cluster and convert
ssh username@your-hpc-cluster.edu
module load apptainer
mkdir -p ~/containers
apptainer build ~/containers/de-limp.sif docker-archive://de-limp-docker.tar
rm de-limp-docker.tar  # Clean up
```

---

## ⚙️ Part 2: First-Time Setup

After you have the `de-limp.sif` file (from any option above), set up your directories:

```bash
# Create organized directory structure
mkdir -p ~/containers     # For .sif files
mkdir -p ~/jobs           # For SLURM scripts
mkdir -p ~/logs           # For job output logs
mkdir -p ~/data           # For your proteomics data
mkdir -p ~/results        # For analysis results

# Verify the container works
apptainer exec ~/containers/de-limp.sif R --version
# Should show: R version 4.5.0 (2025-04-11)
```

---

## 🖥️ Part 3: Interactive Usage with Port Forwarding

### Method 1: Simple One-Step Launch

**On your local computer (Terminal 1):**
```bash
# SSH with port forwarding (replace with your cluster address)
ssh -L 7860:localhost:7860 username@your-hpc-cluster.edu

# Once connected, request an interactive compute node
# Adjust partition name and resources for your cluster
salloc --time=4:00:00 --mem=32GB --cpus-per-task=8

# Load Apptainer and run DE-LIMP
module load apptainer
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"
```

**On your local computer (Browser):**
```
Open: http://localhost:7860
```

### Method 2: Two-Step Port Forwarding (More Stable for Long Sessions)

**On your local computer (Terminal 1):**
```bash
# SSH to cluster
ssh username@your-hpc-cluster.edu

# Request interactive node
salloc --time=8:00:00 --mem=32GB --cpus-per-task=8

# IMPORTANT: Note the exact compute node name shown
# It will be something like: "compute-0-42" or "node123"
```

**On your local computer (Terminal 2):**
```bash
# Set up port forwarding (replace NODE_NAME and cluster address)
ssh -L 7860:NODE_NAME:7860 username@your-hpc-cluster.edu

# Example:
# ssh -L 7860:compute-0-42:7860 username@farm.hpc.ucdavis.edu
```

**Back in Terminal 1 (on cluster):**
```bash
# Load Apptainer and run
module load apptainer
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"
```

**On your local computer (Browser):**
```
Open: http://localhost:7860
```

### Interactive Session Helper Script

Save as `~/run-delimp-interactive.sh` on cluster:
```bash
#!/bin/bash
# Interactive DE-LIMP launcher with port forwarding instructions

# Get compute node
echo "Requesting compute node..."
salloc --time=8:00:00 --mem=32GB --cpus-per-task=8 << 'SCRIPT'

# Show node info and port forwarding command
CLUSTER_HOST=$(hostname -f | sed 's/.*\.//')  # Extract cluster domain
echo "==================================="
echo "Running on: $(hostname)"
echo ""
echo "To access DE-LIMP, open a NEW terminal and run:"
echo "ssh -L 7860:$(hostname):7860 $USER@$(hostname -f)"
echo ""
echo "Then open browser to: http://localhost:7860"
echo "==================================="
echo ""
echo "Starting DE-LIMP..."

# Load and run
module load apptainer
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"

SCRIPT
```

---

## 🚀 Part 4: Batch Job Submission

### Basic Batch Job

Create `~/jobs/delimp-batch.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=de-limp
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/delimp-%j.log
#SBATCH --error=logs/delimp-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@yourinstitution.edu

# Optional: Specify partition (adjust for your cluster)
# #SBATCH --partition=compute

# Load Apptainer module (name may vary by cluster)
module load apptainer  # or: module load singularity

# Print job info
echo "Job started: $(date)"
echo "Running on: $(hostname)"
echo "Job ID: $SLURM_JOB_ID"
echo "Working directory: $(pwd)"

# Run DE-LIMP with data binding
apptainer exec \
  --bind ${HOME}/data:/data \
  --bind ${HOME}/results:/results \
  ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"

echo "Job finished: $(date)"
```

### Submit and Monitor:
```bash
# Create log directory
mkdir -p ~/logs

# Submit job
sbatch ~/jobs/delimp-batch.slurm

# Check job status
squeue -u $USER

# View output in real-time
tail -f ~/logs/delimp-<jobid>.log

# After job starts, set up port forwarding to access the app
# Find the node name from the log file:
NODE=$(grep "Running on:" ~/logs/delimp-*.log | tail -1 | awk '{print $3}')
# Then on your local computer:
# ssh -L 7860:$NODE:7860 username@your-hpc-cluster.edu
```

---

## 📁 Part 5: Data Binding

### Bind Your Data Directories

```bash
# Example: Bind proteomics data directory
apptainer exec \
  --bind /share/proteomics/data:/data:ro \
  --bind ${HOME}/results:/results:rw \
  ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"
```

**Common bind patterns:**
- `--bind /path/on/host:/data:ro` - Read-only data
- `--bind /path/on/host:/results:rw` - Read-write results
- `--bind ${HOME}/downloads:/downloads` - Downloads folder

### Access in Shiny App

In the running app, your bound directories appear at:
- `/data` - Your bound data directory
- `/results` - Your bound results directory

### 📈 XIC Chromatogram Viewing on HPC

DIA-NN generates `_xic` directories containing per-file `.xic.parquet` files alongside the main report. To use the XIC Viewer on HPC:

1. **Bind the XIC directory** alongside your data:
   ```bash
   apptainer exec \
     --bind /path/to/diann/output:/data:ro \
     ~/containers/de-limp.sif \
     R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"
   ```
   Make sure the bind path includes both the report `.parquet` and the `_xic/` sibling directory.

2. **In the app**, paste the path to the XIC directory in the sidebar under "5. XIC Viewer" (e.g., `/data/report_xic/`), then click "Load XICs".
   - The app auto-detects DIA-NN 1.x vs 2.x format
   - If mobilogram files with non-zero data are found (timsTOF/PASEF), an ion mobility toggle appears

3. **Select a protein** in the DE Dashboard and click "📈 XICs" to view chromatograms.
   - Three display modes available: Facet by sample, Facet by fragment, and **Intensity alignment** (Spectronaut-style fragment ratio consistency check)

**Tip:** XIC files can be large (hundreds of MB per sample). Ensure your `--mem` allocation is sufficient — 32GB+ recommended when viewing XICs for datasets with many samples.

---

## 🔬 DIA-NN Search on HPC

DE-LIMP v3.0 can submit DIA-NN database searches to your HPC cluster's SLURM scheduler directly from the **New Search** tab — no command line needed.

### Prerequisites
- DIA-NN must be installed on your HPC cluster (contact your HPC admin or install in your home directory)
- SSH key-based authentication to the cluster (the app uses `ssh`/`scp` — no passwords stored)

### How It Works
1. Run DE-LIMP on your **local machine** (laptop/desktop) or inside the HPC container
2. Go to the **New Search** tab and select the **HPC (SSH/SLURM)** backend
3. Enter your SSH connection details and click **Test Connection**
4. Browse raw files and FASTA databases on the remote cluster
5. Configure search settings (standard or phosphoproteomics mode) and submit
6. The app generates a SLURM script, uploads it via SCP, and submits via `sbatch`
7. Monitor progress in the job queue — results auto-load when complete

### Notes
- The job queue is non-blocking — submit multiple searches and continue using the app
- Search parameters are automatically captured in the Methodology tab
- Built-in UniProt FASTA downloader and 6 contaminant libraries
- Job queue persists across app restarts; active jobs resume polling

> **Tip:** If running DE-LIMP inside the HPC container, you can still use the SSH backend to submit searches to the same cluster — you just SSH from the container back to the login node.

---

## 🔧 Troubleshooting

### Issue: "Cannot bind mount: directory doesn't exist"
**Solution:** Create directories first
```bash
mkdir -p ~/data ~/results
```

### Issue: "Port already in use"
**Solution:** Use a different port
```bash
# Try port 8080, 8888, or random high port
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=8080)"
```

### Issue: "Out of memory"
**Solution:** Request more memory. MOFA2 training and large GSEA analyses can be memory-intensive.
```bash
salloc --mem=128GB --cpus-per-task=16
```
**Recommended minimums:** 32GB for basic DE analysis, 64GB if using MOFA2 or GSEA with large datasets.

### Issue: Port forwarding not working
**Solution:** Check firewall and try two-step method
```bash
# On Mac, verify forwarding:
netstat -an | grep 7860

# Should show: tcp4  0  0  127.0.0.1.7860  *.*  LISTEN
```

### Issue: Slow performance
**Solution:** Use more CPUs and memory
```bash
#SBATCH --cpus-per-task=32
#SBATCH --mem=256GB
```

---

## 📊 Part 6: Example Workflows

### Workflow 1: Quick Interactive Analysis (Recommended for First-Time Users)
```bash
# 1. Get the container (if you haven't already)
apptainer pull ~/containers/de-limp.sif docker://huggingface.co/spaces/brettsp/de-limp-proteomics

# 2. SSH with port forwarding (replace with your cluster)
ssh -L 7860:localhost:7860 username@your-hpc-cluster.edu

# 3. Start interactive session (adjust resources)
salloc --time=4:00:00 --mem=32GB --cpus-per-task=8

# 4. Run DE-LIMP
module load apptainer
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/', host='0.0.0.0', port=7860)"

# 5. Open browser on your computer: http://localhost:7860
```

### Workflow 2: Long-Running Batch Analysis
```bash
# 1. Create SLURM script (see Part 4)
# 2. Submit job
sbatch ~/jobs/delimp-batch.slurm

# 3. Monitor job status
squeue -u $USER
tail -f ~/logs/delimp-*.log

# 4. When job starts, set up port forwarding in a new terminal
# Find the compute node from the log file:
NODE=$(grep "Running on:" ~/logs/delimp-*.log | tail -1 | awk '{print $3}')

# Then connect with port forwarding:
ssh -L 7860:${NODE}:7860 username@your-hpc-cluster.edu

# 5. Access in browser: http://localhost:7860
```

### Workflow 3: Automated Pipeline (No GUI)
```bash
# For batch processing without interactive Shiny interface
# Create R script with your analysis: ~/scripts/run-analysis.R
# Then execute it:
apptainer exec ~/containers/de-limp.sif \
  Rscript ~/scripts/run-analysis.R
```

---

## 🎓 Cluster-Specific Examples

### UC Davis (FARM, HPC1, HPC2)
```bash
# Login
ssh username@farm.hpc.ucdavis.edu

# Common partitions: high, low, med, long
salloc --partition=high --time=8:00:00 --mem=64GB --cpus-per-task=16

# Module
module load apptainer

# Support: hpc-help@ucdavis.edu
# Docs: https://hpc.ucdavis.edu
```

### General HPC Clusters
```bash
# Check available partitions
sinfo

# Check available modules
module avail

# Check node availability
squeue

# Request resources (adjust for your system)
salloc --time=8:00:00 --mem=64GB --cpus-per-task=16
```

### Common Module Names
- `module load apptainer`
- `module load singularity`
- `module load singularity/3.8`
- Check your cluster documentation for the exact name

### Getting Help
- Check your cluster's documentation website
- Contact your HPC support team
- DE-LIMP issues: https://github.com/bsphinney/DE-LIMP/issues

---

## 💡 Tips & Best Practices

1. **Use screen/tmux** for persistent sessions
   ```bash
   screen -S delimp
   # Run your commands
   # Ctrl+A, D to detach
   # screen -r delimp to reattach
   ```

2. **Monitor resource usage**
   ```bash
   # During interactive session
   top
   htop  # if available
   ```

3. **Save session data** regularly using DE-LIMP's built-in save feature

4. **Keep container updated** - rebuild when app updates

5. **Test with small datasets** first to verify setup

---

**Last updated:** 2026-02-20
**DE-LIMP version:** v3.0.0
**Apptainer/Singularity:** Compatible with versions 3.0+
