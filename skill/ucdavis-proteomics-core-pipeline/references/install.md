# Installation reference — what's automatic, what isn't

Design goal: **a biologist with a fresh laptop runs the skill and it works.**
`setup.sh` does the installing; this doc explains what it does and the one
genuine manual fallback.

## What `setup.sh` installs automatically (no admin / no sudo)

It downloads **micromamba** (a single static binary, ~5 MB) into
`~/.proteomics-pipeline/micromamba/` if no conda/mamba is already present, then
creates one conda environment containing:

| package | why |
|---|---|
| `python=3.11` + `pyarrow` + `pyyaml` | the skill's Python scripts + parquet adapters |
| `r-base` (≥4.5) | the DE step |
| `bioconductor-limpa` (1.2.5, noarch) | the DPC-Quant DE pipeline |
| `bioconductor-limma` | both DE pipelines |
| `r-arrow`, `r-dplyr`, `r-tidyr` | reading parquet + the MaxLFQ matrix builder |
| `sage-proteomics` | the DDA search engine |
| `proteowizard` (msconvert) | `.d`/`.raw` → mzML for Sage — **Linux only on bioconda** |

If the conda solve drops limpa, `setup.sh` installs it via `BiocManager::install("limpa")`.

Everything lands under `~/.proteomics-pipeline/`. `source activate.sh` puts it on
PATH so `Rscript`, `python`, `sage`, `msconvert` resolve to the env. Nothing is
installed system-wide; deleting `~/.proteomics-pipeline/` fully uninstalls it.

`setup.sh --check` reports readiness without installing anything.

## The search engines

- **Sage** (DDA): installed by `setup.sh` from bioconda. `acquire_tools.sh` then
  finds it on PATH — no separate download.
- **DIA-NN** (DIA): license-gated, not on conda. By platform:
  - **Linux:** `acquire_tools.sh` downloads the free DIA-NN *Academia* binary.
  - **UC Davis HIVE:** reuse the existing Apptainer `.sif`.
  - **macOS:** *no native build exists.* See the one manual step below.
- **FragPipe** (opt-in): downloaded on demand; MSFragger/IonQuant are license-gated
  and can't be auto-downloaded — `acquire_tools.sh` says so if you opt in.

## The ONE manual step: Docker on macOS (only for DIA-NN / DIA data)

DIA-NN can't run natively on a Mac. `build_diann_docker.sh` builds an image from
DIA-NN's **own official** Academia zip + Dockerfile (no third-party images), but
that needs Docker:

1. Install Docker Desktop (free): https://docs.docker.com/desktop/setup/install/mac-install/
2. Open it once so setup finishes (steady whale icon in the menu bar).
3. Re-run the analysis — `build_diann_docker.sh` builds the image and records it as
   `DIANN_DOCKER_IMAGE`; the AI continues automatically.

DDA data (Sage) and all DE steps need **no** Docker. A Mac user with DDA data that
is already in mzML never touches Docker at all.

## macOS + Sage + Bruker/Thermo

`msconvert` is Linux-only on bioconda, so on a Mac the skill can't auto-convert
`.d`/`.raw` to mzML for Sage. Options, in order of preference:
1. If the data is DIA, use DIA-NN (reads `.d`/`.raw` natively) — no conversion.
2. Convert to mzML elsewhere (a Windows/Linux box, or ProteoWizard Docker) and point
   the skill at the `.mzML` files.
3. Run the whole skill on HIVE/Linux where msconvert is available.

## Windows

Run inside **WSL2** (the DE-LIMP-recommended path) and follow the Linux flow — the
conda env, Sage, and the DIA-NN Linux binary all work there. Native Windows is not
a target for this skill.
