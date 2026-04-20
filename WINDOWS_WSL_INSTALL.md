# Windows WSL Installation Guide

Run DE-LIMP natively on Windows **without Docker**, using Windows Subsystem for Linux (WSL2). Better fit than Docker if:

- Your IT department blocks Docker Desktop
- You've hit the Docker Desktop licensing or performance issues
- You want faster startup (~30s vs Docker's ~3 minutes)
- You already use WSL for other work

**DIA-NN support**: Local searches work. The installer downloads the official DIA-NN Linux binary and .NET 8 runtime into WSL. Thermo `.raw`, Bruker `.d`, and mzML all supported — same as the Docker install. HPC/SSH submit mode also works if you prefer offloading heavy runs.

On first install you're asked to accept the DIA-NN academic license. Decline it if you only want HPC submission or don't need local searches — the rest of the app still works, and you can run `bash ~/delimp_wsl_setup.sh diann` later to change your mind.

---

## Prerequisites

You need two things installed on Windows **before** running the launcher:

1. **Git for Windows** — needed to `git clone` the repo and pull updates.
   - Download: https://git-scm.com/download/win
   - Install with default options. This gives you the `git` command in PowerShell.
   - Verify: open PowerShell, run `git --version`. Should print something like `git version 2.44.0.windows.1`.

2. **WSL2 with Ubuntu** — the launcher needs this to create the Linux environment.
   - Open PowerShell **as Administrator** and run:
     ```powershell
     wsl --install
     ```
   - Restart Windows when prompted. On first login to Ubuntu you'll be asked to create a username and password — these are **WSL-only credentials**, unrelated to your Windows account or any HPC login.

You do **not** need R, RStudio, or any Bioconductor packages installed on Windows itself — everything lives inside the WSL Ubuntu distro that the launcher sets up.

---

## Quick Start

1. **Clone DE-LIMP** (PowerShell — not admin):
   ```powershell
   cd $env:USERPROFILE
   git clone https://github.com/bsphinney/DE-LIMP.git
   cd DE-LIMP
   ```

2. **Double-click `Launch_DE-LIMP_WSL.bat`** (or run `.\Launch_DE-LIMP_WSL.bat` from PowerShell to see live output).

   > **You may need to click it twice the first time.** If WSL/Ubuntu wasn't fully set up yet, the launcher triggers Ubuntu install and then exits with the message:
   >
   > `Ubuntu install triggered. After it finishes setting up, re-run this launcher.`
   >
   > Complete the Ubuntu first-time setup in the separate window that pops up (create a Linux username/password — these are **WSL-only** credentials, unrelated to your Windows account or HPC login).
   >
   > **⚠️ Write down the Ubuntu password.** The second run of the launcher will prompt for it when installing system packages (`sudo apt-get`). There's no way to recover it if you forget — you'd have to wipe the WSL distro and start over. Ubuntu does **not** echo the password as you type it, which is normal — just type it and press Enter.
   >
   > Close the Ubuntu window when setup is done, then **double-click the launcher again** to start the real install.

   On that second run the installer begins — takes 20–30 minutes to compile R + Bioconductor packages inside WSL. Watch the console for progress. When you see `Listening on http://0.0.0.0:3838`, your browser opens automatically.

3. **Subsequent runs**: double-click the launcher once. Takes ~30 seconds to start.

---

## Where Your Files Go

During the first install the script pauses and asks:

```
======================== Data Directory ========================
  Where should DE-LIMP store raw files, FASTA, and search output?

  Raw mass-spec files can be 5-10 GB each. You probably want this
  on a Windows drive — ideally an internal SSD with plenty of space
  (D:, E:, etc.) so File Explorer can see the output and your WSL
  virtual disk doesn't balloon.

  Enter a Windows path like:   D:\proteomics\delimp-data
  Or leave blank to use:       ~/.delimp/data  (inside WSL)
================================================================

Data directory [leave blank for WSL-internal default]:
```

**Recommended: type a Windows path.** A 30-file experiment is easily 200 GB, and the WSL virtual disk is painful to grow back down if it fills. Keeping data on a real Windows drive also lets you drag-and-drop results back to File Explorer.

**Accepted input forms:**
- Windows path: `F:\DE-LIMP` or `D:\proteomics\delimp-data`
- Linux-style WSL path: `/mnt/f/DE-LIMP`
- Blank (just press Enter): falls back to `~/.delimp/data` inside the WSL virtual disk

Trailing backslashes are stripped automatically, so `F:\DE-LIMP\` and `F:\DE-LIMP` work the same. If you type an invalid path (drive not mounted, parent folder missing) the prompt loops — it won't close the installer.

Your choice is saved to `~/.delimp/data_dir` and reused on every future launch. Three subfolders — `raw/`, `fasta/`, `output/` — are created automatically inside whatever directory you chose.

**To change it later:**
```powershell
wsl -d Ubuntu -e bash ~/delimp_wsl_setup.sh config-data-dir
```

**Performance note on Windows drives**: files under `/mnt/c/`, `/mnt/d/`, etc. are accessed via WSL's 9p bridge — slightly slower than WSL-native paths. In practice this is fine for `.raw` and `.fasta` files (sequential reads), but scanning very large `.d` directories is noticeably slower than doing the same on `/quobyte/` on HPC. If you hit this, the fallback is to keep data inside WSL (blank prompt) but you'll need to manage WSL disk size yourself.

**Where other things live** (inside WSL, not user-configurable):
- Code repo: `~/.delimp/DE-LIMP/`
- R packages: `~/.delimp/R-lib/`
- DIA-NN binary + libs: `~/.delimp/diann/`
- SSH keys: `~/.ssh/` (WSL home — full-fat Unix permissions)

All accessible from Windows File Explorer at `\\wsl.localhost\Ubuntu\home\<you>\` if you need to poke at them.

---

## SSH Key Setup (for HPC Submission)

SSH in WSL is much simpler than in Docker — WSL Ubuntu has full-fat SSH with normal Unix permissions.

You need to run these commands from inside an Ubuntu shell (not PowerShell, not Git Bash — the key's file permissions only get set correctly when it's created inside the WSL filesystem).

**Opening an Ubuntu terminal — pick whichever is easiest:**

- **From the Start menu:** click Start, type `Ubuntu`, press Enter. A black Ubuntu terminal window opens at your home directory.
- **From Windows Terminal** (preinstalled on Win11, free on Win10): click the `▾` next to the tab and pick `Ubuntu` from the dropdown.
- **From any PowerShell window:** run `wsl -d Ubuntu` — that drops the current PowerShell session into the Ubuntu shell.

In any of those you'll see a prompt like `protcore@DESKTOP:~$`. You're in Linux. Now:

1. **Generate a key** if you don't have one:
   ```bash
   ssh-keygen -t ed25519 -C "your_email@example.com"
   ```
   Press Enter at each prompt to accept the defaults (saves to `~/.ssh/id_ed25519`, no passphrase — DE-LIMP can't unlock passphrase-protected keys automatically).

2. **Copy the public key to your HPC cluster:**
   ```bash
   ssh-copy-id your-hpc-username@hive.hpc.ucdavis.edu
   ```
   You'll type your HPC password once; after that SSH uses the key. If your cluster disables password logins entirely, email your cluster admin the public key contents (`cat ~/.ssh/id_ed25519.pub`) and ask them to add it to your `~/.ssh/authorized_keys` on the cluster.

3. **In the DE-LIMP SSH panel** (the app UI), point to `~/.ssh/id_ed25519` and enter your HPC username. Click Test Connection — should succeed without any further setup.

No chmod, no line-ending fixes, no container restart needed — WSL treats this exactly like a real Linux box.

**To get back to PowerShell** when you're done in the Ubuntu shell, type `exit` and press Enter.

### Finding your SSH key on the Windows side

The key lives in WSL's Linux filesystem, but Windows File Explorer can see it through a special network path. Paste this into File Explorer's address bar:

```
\\wsl.localhost\Ubuntu\home\<your-linux-username>\.ssh
```

Replace `<your-linux-username>` with whatever username you created when setting up Ubuntu. If you don't remember it, run `whoami` in the Ubuntu shell.

You'll see two files:

- `id_ed25519` — **private key** (keep secret, never share)
- `id_ed25519.pub` — **public key** (safe to share — this is what goes on your HPC cluster, GitHub, etc.)

**If you just want the public key text to copy-paste** into a web form or email:

```bash
cat ~/.ssh/id_ed25519.pub
```

Prints a single line like `ssh-ed25519 AAAAC3... your_email@example.com`. Select and copy it from the terminal.

**If you want a Windows-side backup copy** (e.g. to share with PuTTY or OpenSSH for Windows):

```bash
cp ~/.ssh/id_ed25519 /mnt/c/Users/<you>/Documents/id_ed25519
cp ~/.ssh/id_ed25519.pub /mnt/c/Users/<you>/Documents/id_ed25519.pub
```

> **Caution:** The copy on the Windows side won't carry Unix 0600 permissions. If you ever try to use it *from* WSL again, SSH will reject it with "bad permissions." The canonical working copy should always stay in `~/.ssh/` inside WSL — the Windows copy is for backup or use by other tools.

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
