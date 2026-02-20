# Windows Docker Installation Guide

Run DE-LIMP with full DIA-NN search capability on Windows — no R installation required.

> **DIA-NN License Notice:** DIA-NN is developed by Vadim Demichev and is **free for academic and non-commercial use**. It is **not open source and cannot be redistributed**. The build script in Step 2 downloads DIA-NN directly from the [official GitHub release](https://github.com/vdemichev/DiaNN/releases) and builds a local Docker image on your machine — the binary never leaves your computer. By using this script, you agree to the [DIA-NN license terms](https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md). For commercial use, contact the author directly.

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

- **The `data/` folder is created automatically** when you run `docker compose up` — you don't need to create it manually
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

## Citation

If you use DIA-NN in your research, please cite:

> Demichev V, Messner CB, Vernardis SI, Lilley KS, Ralser M. DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. *Nature Methods*. 2020;17(1):41-44.

DIA-NN license: https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md
