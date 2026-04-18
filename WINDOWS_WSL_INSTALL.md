# Windows WSL Installation Guide

Run DE-LIMP natively on Windows **without Docker**, using Windows Subsystem for Linux (WSL2). Better fit than Docker if:

- Your IT department blocks Docker Desktop
- You've hit the Docker Desktop licensing or performance issues
- You want faster startup (~30s vs Docker's ~3 minutes)
- You already use WSL for other work

**DIA-NN support**: Local searches work. The installer downloads the official DIA-NN Linux binary and .NET 8 runtime into WSL. Thermo `.raw`, Bruker `.d`, and mzML all supported — same as the Docker install. HPC/SSH submit mode also works if you prefer offloading heavy runs.

On first install you're asked to accept the DIA-NN academic license. Decline it if you only want HPC submission or don't need local searches — the rest of the app still works, and you can run `bash ~/delimp_wsl_setup.sh diann` later to change your mind.

---

## Quick Start

1. **Install WSL2** (once per Windows machine, Admin PowerShell):
   ```powershell
   wsl --install
   ```
   Restart Windows when prompted.

2. **Clone DE-LIMP**:
   ```powershell
   cd $env:USERPROFILE
   git clone https://github.com/bsphinney/DE-LIMP.git
   cd DE-LIMP
   ```

3. **Double-click `Launch_DE-LIMP_WSL.bat`.**

   On first run this takes 20–30 minutes to install R + Bioconductor packages inside WSL. Watch the console for progress. When you see `Listening on http://0.0.0.0:3838`, your browser opens automatically.

4. **Subsequent runs**: double-click the launcher again. Takes ~30 seconds.

---

## Where Your Files Go

- **Code & R libraries**: `~/.delimp/` inside WSL (not visible in Windows File Explorer directly — access via `\\wsl.localhost\Ubuntu\home\<you>\.delimp\`)
- **Raw data, FASTA, output, SSH keys**: `~/.delimp/data/` inside WSL by default

**To share `data/` with Windows File Explorer**, set the env var before launching. In the `.bat` you can add (edit and re-run):
```
set DELIMP_DATA_DIR=/mnt/c/Users/%USERNAME%/DE-LIMP/data
```
Then files in `C:\Users\<you>\DE-LIMP\data\raw\` are visible to the app. Note: 9p-mounted files (under `/mnt/c/`) are slower than WSL-native files — fine for small FASTAs, OK for `.raw`, sluggish for very large `.d` directories.

---

## SSH Key Setup (for HPC Submission)

SSH in WSL is much simpler than in Docker — WSL Ubuntu has full-fat SSH with normal Unix permissions.

1. **Generate a key inside WSL** if you don't have one:
   ```bash
   wsl -d Ubuntu
   ssh-keygen -t ed25519 -C "your_email@example.com"
   ```

2. **Copy the public key to your HPC cluster**:
   ```bash
   ssh-copy-id your-hpc-username@hive.hpc.ucdavis.edu
   ```

3. **In the DE-LIMP SSH panel**, point to `~/.ssh/id_ed25519` (or whatever your key is named).

No chmod, no line-ending fixes, no container restart needed.

---

## Troubleshooting

### `wsl --install` says "no access"

Run PowerShell as Administrator. WSL install requires admin rights on the first setup.

### The launcher hangs at "Copying setup script into WSL"

Check that `delimp_wsl_setup.sh` sits next to `Launch_DE-LIMP_WSL.bat` in the same folder. If DE-LIMP is inside OneDrive, try cloning it to `C:\DE-LIMP\` instead — OneDrive paths sometimes confuse `wslpath`.

### R package install fails for `limpa` (Bioconductor)

limpa requires a recent Bioconductor (3.22+). Check the R version inside WSL:
```bash
wsl -d Ubuntu -e R --version
```
If it's older than 4.5, update via:
```bash
wsl -d Ubuntu
sudo apt update
sudo apt install --only-upgrade r-base r-base-dev
```

### Browser opens before app is ready

The launcher waits 90 seconds before opening the browser. If the app isn't up yet, just refresh `http://localhost:3838`.

### App port 3838 already in use

Something else is listening on that port (another DE-LIMP instance, Shiny Server, etc.). Find and stop it, or change the port:
```bash
wsl -d Ubuntu -e bash -c "DELIMP_PORT=3839 bash ~/delimp_wsl_setup.sh run"
```

### Update to the latest DE-LIMP code

```bash
wsl -d Ubuntu -e bash -c "bash ~/delimp_wsl_setup.sh update"
```

---

## Uninstall

Everything lives in `~/.delimp/` inside WSL. To wipe it:

```bash
wsl -d Ubuntu -e bash -c "rm -rf ~/.delimp ~/delimp_wsl_setup.sh"
```

The apt-installed R stays — uninstall via `sudo apt remove r-base r-base-dev` if you want to remove that too.

---

## When to Use This vs Docker vs HPC Apptainer

| Install path | When it fits | Install time | DIA-NN search |
|---|---|---|---|
| **Docker + SSH** (Windows) | You already use Docker, want isolation | ~45 min first run | Local + HPC |
| **WSL** (Windows) | No Docker, lighter-weight, faster | ~30 min first run | Local + HPC |
| **HPC Apptainer** | You work mostly on HIVE directly | ~15 min | HPC only |
| **Hugging Face Spaces** | Zero install, analysis-only (have a report.parquet already) | 0 | None |

WSL is the "middle weight" option — lighter than Docker, heavier than HF Spaces, gives you full UI + HPC integration.
