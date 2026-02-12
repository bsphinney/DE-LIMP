# DE-LIMP HPC Deployment Guide (UC Davis)

This guide covers deploying DE-LIMP on UC Davis HPC clusters (FARM, HPC1/HPC2) using Apptainer/Singularity containers.

---

## 🏗️ Part 1: Build Docker Image Locally

### Prerequisites
- Docker Desktop installed on Mac
- ~10GB free disk space
- 30-45 minutes for build

### Step 1: Build the Image
```bash
# Navigate to project directory
cd /Users/brettphinney/Documents/claude

# Build Docker image (takes 30-45 minutes)
docker build -t de-limp:latest .

# Verify build
docker images | grep de-limp
# Should show: de-limp   latest   <image-id>   <size>
```

### Step 2: Test Locally (Optional)
```bash
# Run the container locally to verify it works
docker run -p 7860:7860 de-limp:latest

# Open browser to http://localhost:7860
# Ctrl+C to stop
```

### Step 3: Convert to Apptainer Format

**Option A: Direct conversion (if Apptainer installed on Mac)**
```bash
# Install Apptainer on Mac (if not already installed)
brew install apptainer

# Convert Docker image to .sif
apptainer build de-limp.sif docker-daemon://de-limp:latest
```

**Option B: Save for conversion on cluster**
```bash
# Save Docker image as tar archive
docker save de-limp:latest -o de-limp-docker.tar

# The .tar file will be ~3-4GB
ls -lh de-limp-docker.tar
```

---

## 📤 Part 2: Transfer to UC Davis HPC

### Transfer the Image

**If you built .sif file locally:**
```bash
# Transfer to FARM cluster
scp de-limp.sif username@farm.hpc.ucdavis.edu:~/containers/

# Or to HPC1/HPC2
scp de-limp.sif username@hpc1.hpc.ucdavis.edu:~/containers/
```

**If you saved .tar file:**
```bash
# Transfer tar file
scp de-limp-docker.tar username@farm.hpc.ucdavis.edu:~/

# SSH to cluster and convert
ssh username@farm.hpc.ucdavis.edu
mkdir -p ~/containers
module load apptainer
apptainer build ~/containers/de-limp.sif docker-archive://de-limp-docker.tar
rm de-limp-docker.tar  # Clean up
```

---

## 🖥️ Part 3: Interactive Usage with Port Forwarding

### Method 1: Single-Command Interactive Session

**On your Mac (Terminal 1):**
```bash
# SSH with port forwarding
ssh -L 7860:localhost:7860 username@farm.hpc.ucdavis.edu

# Once connected, request interactive node
salloc --partition=high --time=4:00:00 --mem=32GB --cpus-per-task=8

# Note the node name (e.g., "compute-0-42")
# Load Apptainer and run
module load apptainer
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/app.R', host='0.0.0.0', port=7860)"
```

**On your Mac (Browser):**
```
Open: http://localhost:7860
```

### Method 2: Two-Step Port Forwarding (More Stable)

**On your Mac (Terminal 1):**
```bash
# SSH to cluster
ssh username@farm.hpc.ucdavis.edu

# Request interactive node
salloc --partition=high --time=8:00:00 --mem=32GB --cpus-per-task=8

# Note the exact node name shown, e.g., "compute-0-42"
```

**On your Mac (Terminal 2):**
```bash
# Set up port forwarding (replace NODE_NAME with actual node)
ssh -L 7860:compute-0-42:7860 username@farm.hpc.ucdavis.edu
```

**Back in Terminal 1 (on cluster):**
```bash
# Load Apptainer and run
module load apptainer
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/app.R', host='0.0.0.0', port=7860)"
```

**On your Mac (Browser):**
```
Open: http://localhost:7860
```

### Interactive Session Script

Save as `~/run-delimp-interactive.sh` on cluster:
```bash
#!/bin/bash
# Interactive DE-LIMP launcher

# Get compute node
echo "Requesting compute node..."
salloc --partition=high --time=8:00:00 --mem=32GB --cpus-per-task=8 << 'SCRIPT'

# Show node info
echo "==================================="
echo "Running on: $(hostname)"
echo "Setup port forwarding from your Mac:"
echo "ssh -L 7860:$(hostname):7860 $USER@farm.hpc.ucdavis.edu"
echo "==================================="
echo ""
echo "Starting DE-LIMP..."

# Load and run
module load apptainer
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/app.R', host='0.0.0.0', port=7860)"

SCRIPT
```

---

## 🚀 Part 4: Batch Job Submission

### Basic Batch Job

Create `~/jobs/delimp-batch.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=de-limp
#SBATCH --partition=high
#SBATCH --time=12:00:00
#SBATCH --mem=64GB
#SBATCH --cpus-per-task=16
#SBATCH --output=logs/delimp-%j.log
#SBATCH --error=logs/delimp-%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your.email@ucdavis.edu

# Load modules
module load apptainer

# Print job info
echo "Job started: $(date)"
echo "Running on: $(hostname)"
echo "Working directory: $(pwd)"

# Run DE-LIMP
apptainer exec \
  --bind ${HOME}/data:/data \
  --bind ${HOME}/results:/results \
  ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/app.R', host='0.0.0.0', port=7860)"

echo "Job finished: $(date)"
```

### Submit the job:
```bash
# Create log directory
mkdir -p ~/logs

# Submit job
sbatch ~/jobs/delimp-batch.slurm

# Check job status
squeue -u $USER

# View output
tail -f ~/logs/delimp-<jobid>.log
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
  R -e "shiny::runApp('/srv/shiny-server/app.R', host='0.0.0.0', port=7860)"
```

**Common bind patterns:**
- `--bind /path/on/host:/data:ro` - Read-only data
- `--bind /path/on/host:/results:rw` - Read-write results
- `--bind ${HOME}/downloads:/downloads` - Downloads folder

### Access in Shiny App

In the running app, your bound directories appear at:
- `/data` - Your bound data directory
- `/results` - Your bound results directory

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
  R -e "shiny::runApp('/srv/shiny-server/app.R', host='0.0.0.0', port=8080)"
```

### Issue: "Out of memory"
**Solution:** Request more memory
```bash
salloc --mem=128GB --cpus-per-task=16
```

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

### Workflow 1: Quick Interactive Analysis
```bash
# 1. SSH with port forwarding
ssh -L 7860:localhost:7860 username@farm.hpc.ucdavis.edu

# 2. Start interactive session
salloc --partition=high --time=4:00:00 --mem=32GB --cpus-per-task=8

# 3. Run DE-LIMP
module load apptainer
apptainer exec ~/containers/de-limp.sif \
  R -e "shiny::runApp('/srv/shiny-server/app.R', host='0.0.0.0', port=7860)"

# 4. Open browser on Mac: http://localhost:7860
```

### Workflow 2: Long-Running Batch Analysis
```bash
# 1. Create SLURM script (see Part 4)
# 2. Submit job
sbatch ~/jobs/delimp-batch.slurm

# 3. Monitor
squeue -u $USER
tail -f ~/logs/delimp-*.log

# 4. When job starts, set up port forwarding
# Find node in log file
NODE=$(grep "Running on:" ~/logs/delimp-*.log | awk '{print $3}')
ssh -L 7860:${NODE}:7860 username@farm.hpc.ucdavis.edu

# 5. Access in browser: http://localhost:7860
```

### Workflow 3: Automated Pipeline (No GUI)
```bash
# For batch processing without interactive Shiny interface
# Create R script: ~/scripts/run-analysis.R
apptainer exec ~/containers/de-limp.sif \
  Rscript ~/scripts/run-analysis.R
```

---

## 🎓 UC Davis Specific Notes

### FARM Cluster
- **Login node:** farm.hpc.ucdavis.edu
- **Common partitions:** high, low, med, long
- **Max time:** 14 days (long partition)
- **Apptainer module:** `module load apptainer`

### HPC1/HPC2
- **Login node:** hpc1.hpc.ucdavis.edu, hpc2.hpc.ucdavis.edu
- **Check available partitions:** `sinfo`
- **Check node availability:** `squeue`

### Getting Help
- HPC Support: hpc-help@ucdavis.edu
- Proteomics Core: Brett Phinney (your account!)
- Documentation: https://hpc.ucdavis.edu

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

**Last updated:** 2026-02-11
**DE-LIMP version:** v2.0.1
**Apptainer/Singularity:** Compatible with versions 3.0+
