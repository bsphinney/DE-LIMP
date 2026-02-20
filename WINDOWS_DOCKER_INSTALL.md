# Windows Docker Installation Guide

Run DE-LIMP with full DIA-NN search capability on Windows — no R installation required.

> **DIA-NN License Notice:** DIA-NN is developed by Vadim Demichev and is **free for academic and non-commercial use**. It is **not open source and cannot be redistributed**. The build script in Step 2 downloads DIA-NN directly from the [official GitHub release](https://github.com/vdemichev/DiaNN/releases) and builds a local Docker image on your machine — the binary never leaves your computer. By using this script, you agree to the [DIA-NN license terms](https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md). For commercial use, contact the author directly.

## What You'll Get

- DE-LIMP running at **http://localhost:3838**
- DIA-NN embedded inside the container, ready for searches
- Your data files shared between Windows and the container via the `data/` folder

## Prerequisites

1. **Docker Desktop for Windows**
   - Download: https://www.docker.com/products/docker-desktop/
   - Install and start Docker Desktop
   - Make sure it's running (whale icon in the system tray)

2. **Git for Windows**
   - Download: https://git-scm.com/download/win
   - Install with default settings

## Step 1: Get DE-LIMP

Open **PowerShell** or **Git Bash** and run:

```
git clone https://github.com/bsphinney/DE-LIMP.git
cd DE-LIMP
git checkout diann-integration
```

If you already have DE-LIMP:

```
cd C:\path\to\DE-LIMP
git checkout diann-integration
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

## Step 4: Run a Search

1. Copy your raw data files (`.d`, `.raw`, or `.mzML`) into the **`data\raw\`** folder
2. Copy your FASTA database files into **`data\fasta\`**
3. In DE-LIMP, go to the **New Search** tab
4. The backend shows **"Local (Embedded)"** — DIA-NN is ready
5. Click **Browse** to select your raw files from `/data/raw/`
6. Select your FASTA from `/data/fasta/` (or use the built-in UniProt downloader)
7. Configure search settings and click **Submit**
8. Results appear in `data\output\` and auto-load into the DE pipeline

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

Only use `--build` again if you pull new code changes.

## Folder Structure

```
DE-LIMP/
  data/
    raw/          <-- Put your .d / .raw / .mzML files here
    fasta/        <-- Put your .fasta files here
    output/       <-- DIA-NN results appear here
```

The `data/` folder is shared between Windows and the container. Files you add on Windows appear instantly inside the app, and search results are written back to your Windows filesystem.

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

### Search runs but is slow

DIA-NN runs under Linux emulation in Docker. Performance is good but not native. For large experiments (50+ files), consider using the HPC/SSH backend to submit to a compute cluster.

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

## Citation

If you use DIA-NN, please cite:

> Demichev V, Messner CB, Vernardis SI, Lilley KS, Ralser M. DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nature Methods*. 2020;17(1):41-44.

DIA-NN license: https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md
