# Windows Docker Installation Guide

Run DE-LIMP with full DIA-NN search capability on Windows — no R installation required.

> **DIA-NN License Notice:** DIA-NN is developed by Vadim Demichev and is **free for academic and non-commercial use**. It is **not open source and cannot be redistributed**. The build script in Step 2 downloads DIA-NN directly from the [official GitHub release](https://github.com/vdemichev/DiaNN/releases) and builds a local Docker image on your machine — the binary never leaves your computer. By using this script, you agree to the [DIA-NN license terms](https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md). For commercial use, contact the author directly.

---

## Quick Start: One-Click Launcher

If you already have Docker Desktop installed and the DIA-NN image built, just **double-click `Launch_DE-LIMP_Docker.bat`** in your DE-LIMP folder. It creates the `data\` subfolders, prompts you for your HPC username (optional), starts the container, waits for the app to come up, and opens your browser. Skip to [Step 4](#step-4-run-a-search) once the app loads.

If you're setting up for the first time or prefer manual steps, follow Steps 1–3 below.

---

## What You'll Get

- DE-LIMP running at **http://localhost:3838**
- DIA-NN embedded inside the container, ready for searches
- A shared `data/` folder where you put your files and get results back

## Prerequisites

1. **Docker Desktop for Windows**
   - Download: https://www.docker.com/products/docker-desktop/
   - Install and start Docker Desktop
   - Make sure it's running (whale icon in the system tray)

2. **Git for Windows**
   - Download: https://git-scm.com/download/win
   - Install with default settings

---

## Step 1: Get DE-LIMP

Open **PowerShell** or **Git Bash** and run:

```
git clone https://github.com/bsphinney/DE-LIMP.git
cd DE-LIMP
```

If you already have DE-LIMP:

```
cd C:\path\to\DE-LIMP
git pull
```

## Step 2: Build the DIA-NN Docker Image

DIA-NN is free for academic use but cannot be redistributed, so you build the image locally. This downloads DIA-NN from the official GitHub release.

**Option A — PowerShell** (recommended for Windows):

```powershell
.\build_diann_docker.ps1
```

**Option B — Git Bash**:

```bash
bash build_diann_docker.sh
```

This takes a few minutes. When finished, you'll see:

```
DIA-NN 2.0 Docker image built successfully!
```

**Verify it worked:**

```
docker images diann
```

You should see `diann` with tag `2.0`.

## Step 3: Build and Start DE-LIMP

```
docker compose up --build
```

**First time** this takes 15-30 minutes (downloads R packages, .NET SDK, etc.). All subsequent runs use cached layers and start in seconds.

When you see this line, the app is ready:

```
Listening on http://0.0.0.0:3838
```

Open **http://localhost:3838** in your browser.

---

## How the Data Folder Works

DE-LIMP runs inside a Docker container, which is an isolated environment. It cannot see your Windows files directly. Instead, a single folder — `data/` inside your DE-LIMP directory — is **shared** between Windows and the container.

```
C:\Users\you\DE-LIMP\
  data\
    raw\          ← Put your .d / .raw / .mzML files here
    fasta\        ← Put your .fasta files here
    output\       ← DIA-NN results appear here automatically
```

### Key points:

- **Create the `data/` subfolders first** — before running `docker compose up`, create `data\raw\`, `data\fasta\`, and `data\output\` on your Windows machine (the `Launch_DE-LIMP_Docker.bat` launcher does this for you automatically if you use it). Docker mounts `./data` into the container, so whatever is (or isn't) on your Windows side is what the app sees.
- **Copy your files into `data\raw\` and `data\fasta\`** before or after starting the app — either works
- **This is a live connection** — files you add in Windows File Explorer appear instantly in the app's file browser, and search results written by DIA-NN appear in `data\output\` on your Windows filesystem
- **Files outside `data/` are not visible** to the app. If your raw files are on a different drive (e.g., `D:\proteomics\experiment1\`), you need to copy or move them into `data\raw\` first
- **Results persist** — stopping and restarting the container does not delete anything in `data/`

### Typical workflow:

1. Open File Explorer and navigate to `C:\Users\you\DE-LIMP\data\raw\`
2. Copy/paste your `.raw` or `.d` files into this folder
3. Do the same for your `.fasta` file in `data\fasta\`
4. In the app, click **Browse** — your files appear under `/data/raw/` and `/data/fasta/`
5. After the search completes, results are in `data\output\` — you can open them with any tool

### Starting a new experiment:

You can reuse the same `data/` folder or clean it out between experiments:
- Delete old files from `data\raw\` and `data\output\` in File Explorer
- Copy in new files
- No need to restart the app — just browse the new files in the Search tab

---

## Step 4: Run a Search

1. Copy your raw data files (`.d`, `.raw`, or `.mzML`) into **`data\raw\`**
2. Copy your FASTA database files into **`data\fasta\`**
3. In DE-LIMP, go to the **New Search** tab
4. The backend shows **"Local (Embedded)"** — DIA-NN is ready
5. Click **Browse** to select your raw files from `/data/raw/`
6. Select your FASTA from `/data/fasta/` (or use the built-in UniProt downloader)
7. Configure search settings and click **Submit**
8. Results appear in `data\output\` and auto-load into the DE pipeline

---

## Stopping and Restarting

**Stop the app:**

Press `Ctrl+C` in the terminal, or run:

```
docker compose down
```

**Restart (fast — no rebuild needed):**

```
docker compose up
```

Only use `--build` again if you pull new code changes. Your `data/` folder and all files in it are preserved across restarts.

---

## Troubleshooting

### `image diann:2.0 not found`

You need to build the DIA-NN image first (Step 2). Run `docker images diann` to check if it exists.

### Build fails at R package installation

Usually a transient network issue. Run `docker compose up --build` again — Docker caches completed layers so it picks up where it left off.

### `Listening on http://0.0.0.0:3838` but browser shows nothing

Try **http://localhost:3838** (not 0.0.0.0). If that doesn't work, check that port 3838 isn't used by another app:

```
netstat -an | findstr 3838
```

### Docker Desktop says "WSL 2 not installed"

Follow the Docker Desktop prompt to install WSL 2, or see: https://learn.microsoft.com/en-us/windows/wsl/install

### Files in `data\raw\` don't appear in the file browser

Make sure Docker Desktop has file sharing enabled for your drive:
- Docker Desktop > Settings > Resources > File Sharing
- Add `C:\` (or your drive letter) if not listed

### My files are on a different drive (D:\, E:\, etc.)

The app can only see files inside the `data/` folder. Copy your files into `data\raw\` and `data\fasta\`. You cannot browse arbitrary folders on your computer from inside the container.

### Search runs but is slow

DIA-NN runs under Linux emulation in Docker on Windows. Performance is reasonable but not native speed. For large experiments (50+ files), consider using the HPC/SSH backend to submit to a compute cluster.

### PowerShell script won't run (`execution policy` error)

Run this first:

```powershell
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
```

### Want to update to the latest code?

```
docker compose down
git pull
docker compose up --build
```

Your data files in `data/` are not affected by code updates.

---

## SSH to HPC — Setting Up Your Key

If you want to submit DIA-NN searches to an HPC cluster instead of running them locally on your Windows machine, you need an SSH key configured. The Docker container picks up SSH keys from a specific folder on your Windows host and copies them into itself with the correct permissions on startup.

### Where to put your SSH private key

Put your private key file in **`data\ssh\`** inside your DE-LIMP folder:

```
C:\Users\you\DE-LIMP\
  data\
    ssh\
      brettsp          ← your HPC private key, named after your HPC username
```

**Important conventions:**
- **Name the file after your HPC username** (e.g., `jsmith`, `acme_lab`). This is convention only — the app reads your HPC username from the `DELIMP_SSH_USER` environment variable set by `Launch_DE-LIMP_Docker.bat` (which prompts you for it) or entered in the SSH panel inside the app. The key filename itself doesn't influence authentication.
- **This must be the private key**, not the `.pub` file. Private keys start with `-----BEGIN OPENSSH PRIVATE KEY-----`.
- **Only one key should be in `data\ssh\`** — if there are multiple, only the first one (alphabetical) is used.
- The container auto-copies it to `/tmp/.ssh/` and sets permissions to 600 on startup.
- If you're running `docker compose up` directly instead of the launcher, you must create `data\ssh\` yourself before starting.

### Register the public key on your HPC cluster

The private key alone isn't enough — the public half must be listed in `~/.ssh/authorized_keys` on the HPC cluster. Ask your cluster admin to add it, or add it yourself (one-time):

```powershell
# Print the public key from the private key you already have
docker compose exec delimp ssh-keygen -y -f /tmp/.ssh/<your-username>

# Copy the output (the ssh-ed25519 AAA... line), then from your HPC account:
# $ echo "ssh-ed25519 AAAA... your_comment" >> ~/.ssh/authorized_keys
# $ chmod 600 ~/.ssh/authorized_keys
```

### Restart after adding or changing a key

The entrypoint only reads `data\ssh\` **on container start**. If you add a key after the app is running:

```powershell
docker compose restart delimp
```

### Test the connection

```powershell
docker compose exec delimp ssh -i /tmp/.ssh/<your-username> -o BatchMode=yes <your-username>@<hpc-host> "echo SSH_OK"
```

Should print `SSH_OK`. If it doesn't, see the troubleshooting section below.

---

## SSH Troubleshooting

### `Load key "...": error in libcrypto`

The private key file is malformed. Common causes:

1. **Missing trailing newline on the `-----END OPENSSH PRIVATE KEY-----` line.** The OpenSSH parser requires the file to end with a newline. Fix from PowerShell:

   ```powershell
   $p = (Resolve-Path .\data\ssh\<your-username>).Path
   $b = [IO.File]::ReadAllBytes($p)
   if ($b[-1] -ne 0x0A) {
       [IO.File]::WriteAllBytes($p, $b + [byte]0x0A)
       "Fixed — new size: $((Get-Item $p).Length)"
   } else { "Already OK" }
   docker compose restart delimp
   ```

2. **Windows CRLF line endings** (from pasting into Notepad or using `Set-Content`/`Out-File`). Use `[IO.File]::ReadAllText` / `WriteAllText` to strip `\r`:

   ```powershell
   $p = (Resolve-Path .\data\ssh\<your-username>).Path
   $c = [IO.File]::ReadAllText($p) -replace "`r`n", "`n"
   [IO.File]::WriteAllText($p, $c)
   docker compose restart delimp
   ```

3. **UTF-8 BOM** (3 bytes at the start: `0xEF 0xBB 0xBF`). Same fix as above — use binary-safe `ReadAllBytes` / `WriteAllBytes` and strip the BOM, or regenerate the key from a clean source.

To see the raw bytes and diagnose: `docker compose exec delimp od -c /tmp/.ssh/<your-username> | head -5`.

### `Load key "...": bad permissions`

The key was read from a 9p-mounted path (`/data/...`) where Windows gives everything permissions 0777. SSH refuses keys with broader than 0600 permissions.

**Fix:** put the key in `data\ssh\` so the entrypoint copies it to `/tmp/.ssh/` and sets 600:

```powershell
Move-Item .\data\id_ed25519 .\data\ssh\<your-username>
docker compose restart delimp
```

Don't point SSH at keys living directly under `/data/` — always use `/tmp/.ssh/`.

### `Permission denied (publickey,password)`

The key itself is valid (no libcrypto/permission error), but the cluster doesn't recognize it. Either:

1. The public half isn't in the cluster's `~/.ssh/authorized_keys` for your user — ask your admin to register it, or add it yourself (see above).
2. The file is named wrong. The diag script infers the HPC username from the filename. If your cluster username is `jsmith` but the file is named `id_ed25519`, DE-LIMP will try to log in as `brettsp` (the default fallback).

To check what key fingerprint your container is offering:

```powershell
docker compose exec delimp ssh-keygen -lf /tmp/.ssh/<your-username>
```

Then ask the cluster admin: "is this fingerprint registered on my account?"

### Using the built-in diagnostic

A diagnostic script ships in the container. It checks SSH key format, mounts, DIA-NN install, network, and common failure patterns, then uploads the report to your HPC cluster (so you can share it if you want help):

```powershell
docker compose exec delimp bash /srv/shiny-server/delimp_diag.sh
```

Look for **"11. Auto-diagnosis — likely issues + fixes"** at the bottom — it pattern-matches common problems and prints the exact PowerShell fix for each.

---

## Performance Troubleshooting

### Scans and searches are slow

Windows shares the `data\` folder into the container through a "9p" virtual filesystem. The default chunk size is small (`msize=65536` = 64 KB), which makes large raw files (`.d` directories, Thermo `.raw`) transfer slowly.

**Check the current setting:**

```powershell
docker compose exec delimp bash -c "mount | grep '/data' | grep -oE 'msize=[0-9]+'"
```

**Option A (recommended) — switch Docker Desktop to virtiofs:**

1. Docker Desktop → Settings → General
2. "Choose file sharing implementation for your containers" → select **virtiofs**
3. Apply & Restart

virtiofs doesn't use 9p at all, so there's no msize to tune. Much faster.

**Option B — increase the 9p chunk size (if virtiofs isn't available):**

1. Create or append to `C:\Users\you\.wslconfig`:
   ```
   [wsl2]
   kernelCommandLine = "9p.msize=262144"
   ```
2. Fully quit Docker Desktop (tray icon → Quit) — not just close
3. Run `wsl --shutdown` in PowerShell
4. Verify `wsl --list --verbose` shows all distros as `Stopped`
5. Restart Docker Desktop

**Caveat:** Docker Desktop runs its own WSL2 distribution (`docker-desktop`) on the same WSL2 VM that respects `.wslconfig`. In practice we've seen cases where the `msize` setting doesn't take effect even after a full restart, likely because Docker Desktop auto-restarts its internal distro before the kernel command line is re-read. If `mount | grep /data` still shows `msize=65536` after following the steps above, use Option A (virtiofs) instead — that bypasses 9p entirely and is the supported long-term path.

---

## Citation

If you use DIA-NN in your research, please cite:

> Demichev V, Messner CB, Vernardis SI, Lilley KS, Ralser M. DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nature Methods*. 2020;17(1):41-44.

DIA-NN license: https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md
