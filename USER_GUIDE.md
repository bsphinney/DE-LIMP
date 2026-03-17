# 📘 DE-LIMP User Guide

Welcome to **DE-LIMP** (Differential Expression & Limpa Proteomics), your interactive dashboard for analyzing DIA-NN proteomics data. This guide covers the complete workflow, from importing data to discovering biological insights with our integrated AI assistant.

DE-LIMP helps you find which proteins are significantly different between experimental conditions (e.g., treatment vs. control). Upload your DIA-NN output, and the app handles normalization, statistics, and visualization -- including interactive volcano plots, heatmaps, pathway enrichment analysis (GSEA), and AI-powered summaries. No programming required.

---

## 1. Getting Started

### Which Installation Method?

| Platform | Recommended | What You Need | Guide |
|----------|-------------|---------------|-------|
| **Just exploring** | Web browser | Nothing | [Hugging Face](https://huggingface.co/spaces/brettsp/de-limp-proteomics) |
| **Windows** | Docker + SSH to HPC | Docker Desktop + SSH key | [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md) |
| **Mac / Linux** | Native R | R 4.5+ and RStudio | Continue reading below |
| **HPC cluster** | Apptainer (alternative) | Singularity/Apptainer | [HPC_DEPLOYMENT.md](HPC_DEPLOYMENT.md) |

> **Windows users:** R package installation on Windows is often problematic. We strongly recommend the Docker + SSH approach — double-click `Launch_DE-LIMP_Docker.bat` to run DE-LIMP locally in Docker, then connect to your HPC cluster via SSH for DIA-NN searches. Shared PC support is built-in. See **[WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md)** for a step-by-step walkthrough.

### Prerequisites (Native R — Mac / Linux)
* **R & RStudio:** Ensure you have R (version 4.5 or newer) installed.
  * **Important:** The limpa package requires R 4.5+ and Bioconductor 3.22+
  * Download R from: https://cloud.r-project.org/
* **Gemini API Key:** Required for AI Chat features (See below).

### 🔑 1.1 How to Obtain a Free Gemini API Key
To use the "Chat with Data" features, you need a key from Google. It is free for standard use.

1.  Go to **[Google AI Studio](https://aistudio.google.com/)**.
2.  Sign in with your Google Account.
3.  In the top-left corner, click the blue button **"Get API key"**.
4.  Click **"Create API key"**.
    * If asked, select "Create API key in new project".
5.  Copy the long string of text that appears (it starts with `AIza...`).
6.  **Paste this key** into the "Gemini API Key" box in the DE-LIMP sidebar.

### 1.2 Launching the App
1.  Open the DE-LIMP project folder in RStudio (or navigate to it in your R console).
2.  Run the app with: `shiny::runApp('/path/to/de-limp/', port=3838, launch.browser=TRUE)`
3.  The dashboard will launch in your default web browser.

---

## 2. Core Analysis Workflow

The sidebar on the left contains collapsible sections for each stage of the workflow. Click a section header to expand or collapse it. The **Upload Data** section is open by default.

### Step 2.1: Upload Data
You have two options to get started:

**Option A: Load Example Data (Recommended for First-Time Users)**
* Click the **"📊 Load Example Data"** button in the sidebar
* The app will automatically download a demo dataset (Affinisep vs Evosep comparison, 46MB)
* This is the fastest way to explore DE-LIMP's features
* The example data showcases a real proteomics experiment with clear differential expression

**Option B: Upload Your Own Data**
* **Input File:** Click **"Browse..."** and select your DIA-NN report file.
    * *Requirement:* The file must be in **`.parquet`** format (the app only accepts `.parquet` files).
    * *What is report.parquet?* When DIA-NN finishes processing your raw files, it creates a `report.parquet` file in the output directory you specified. This is a compact binary format that loads much faster than the older `.tsv` format. Look for it in your DIA-NN output folder alongside other output files.
    * *Download Example:* Available at [GitHub Releases](https://github.com/bsphinney/DE-LIMP/releases/tag/v1.0)
* **Q-Value Cutoff:** This filters DIA-NN precursors at import time -- only precursors with a confidence score (Q-value) below this threshold are kept. It is *not* the DE significance threshold (that uses adj.P.Val after the pipeline runs). Lower values are more strict, keeping only the most confident identifications. The default of 0.01 (1% FDR) is appropriate for most experiments.

### Step 2.2: Assign Groups & Run Pipeline
This is the most critical step for statistical analysis. The workflow is streamlined into one modal dialog.

> **Replicate guidance:** For reliable statistical results, we recommend at least **3 biological replicates per group**.
> - **n=1 per group (NEW in v3.7)**: The pipeline will complete quantification (normalization, protein-level aggregation) but will skip differential expression analysis -- DE requires replicates. You can still explore the Expression Grid, PCA, and Signal Distribution. An informational message explains what was skipped.
> - **n=2 per group**: The pipeline will run, but you can only detect very large changes (>4-fold) with any confidence. Treat results as exploratory, not publication-ready.
> - **n=3 per group**: The standard minimum for publication-quality results. Limma's empirical Bayes moderation helps compensate for the small sample size.
> - **n=4-6+**: Improved sensitivity for detecting smaller fold changes (e.g., 1.3-fold) and more reliable CV estimates in the CV Analysis tab.
> - **Design caution**: If all treatment samples were run on Day 1 and all controls on Day 2, batch and group are confounded -- the statistics cannot separate them. Add batch as a covariate if possible, or run samples in a balanced design.

1.  Click **"Assign Groups & Run Pipeline"** in the sidebar (or it will auto-open after data upload).
2.  **Auto-Guess Groups (Recommended):**
    * Click the **"🪄 Auto-Guess Groups"** button at the top of the modal
    * The app intelligently detects groups (e.g., "Control", "Treatment", "WT", "KO", "Affinisep", "Evosep") based on your filenames
    * For the example data, this automatically assigns samples to "Affinisep" and "Evosep" groups
3.  **Manual Edit (If Needed):**
    * Click on any cell in the **Group** column to type a custom group name
    * You can also edit **Batch**, **Covariate1**, and **Covariate2** columns
4.  **Template Export/Import (NEW in v2.0.1):**
    * **Export Template**: Click **"📥 Export"** to download current group assignments as CSV
      - Saves all table data: File.Name, Group, Batch, and custom covariates
      - Filename format: `DE-LIMP_group_template_YYYYMMDD_HHMMSS.csv`
      - Use for saving configurations or sharing with collaborators
    * **Import Template**: Click **"📤 Import"** to load previously saved group assignments
      - Opens file picker to select a CSV template
      - Validates columns and matches files by name
      - Perfect for reproducible workflows or applying standard patterns to new data
5.  **Customize Covariate Names (Optional):**
    * Use the text inputs above the table to rename "Covariate1" and "Covariate2"
    * Examples: "Sex", "Diet", "Age", "Time_Point", "Instrument"
    * Check the boxes to include covariates in the statistical model
    * Only covariates with 2+ unique values will be used
    * **When to use covariates**: If your samples were run on different days or instruments, adding "Batch" separates batch effects from your treatment effect. If samples come from male and female animals and sex is not your research question, adding "Sex" removes sex-related variation. Only add covariates you have reason to suspect affect protein levels -- adding too many with few samples can reduce statistical power
6.  **Run the Analysis:**
    * Click the **"▶ Run Pipeline"** button at the top of the modal
    * **What happens?** The app uses the `limpa` package to perform DPC normalization and the `limma` package to fit linear models for differential expression
    * The modal will automatically close and navigate to the QC Plots tab
    * Wait for the status to change to **"✅ Complete!"**

### Step 2.3: Select Comparison
* Use the **"Comparison"** dropdown to select which contrast you want to view (e.g., `Evosep - Affinisep`, `Treatment - Control`).

---

## 3. DIA-NN Database Search

> **Already have a report.parquet file?** If your core facility or collaborator already processed the raw data and gave you a `report.parquet` file, you can skip this entire section. Go directly to **[Section 4: Deep Dive](#4-deep-dive-the-data-overview--grid-view)** -- you only need to upload your file and run the pipeline (Section 2).

The **"New Search"** tab lets you submit DIA-NN database searches directly from DE-LIMP. This is for advanced users who need to process raw mass spectrometry files (.raw, .d, .mzML) into a report. Three backends are supported:

| Backend | When to Use | How It Works |
| :--- | :--- | :--- |
| **Local (Embedded)** | Docker Compose deployment (Windows) | DIA-NN binary runs inside the same container as DE-LIMP |
| **Local (Docker)** | Mac/Linux with Docker installed | DIA-NN runs in a separate Docker container |
| **HPC (SSH/SLURM)** | Access to a compute cluster | Jobs submitted via SLURM; results downloaded via SCP |

> **Note:** The New Search tab only appears when at least one backend is detected. It is not shown on the Hugging Face web version.
>
> **Windows users:** The easiest setup is `docker compose up` which gives you the Local (Embedded) backend with no R installation. See [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md).
>
> **DIA-NN License:** DIA-NN is developed by [Vadim Demichev](https://github.com/vdemichev/DiaNN) and is **free for academic and non-commercial use only**. It cannot be redistributed. DE-LIMP does not bundle DIA-NN — the build scripts download it directly from the [official GitHub release](https://github.com/vdemichev/DiaNN/releases) and create a local Docker image on your machine. By using the DIA-NN search features, you agree to the [DIA-NN license terms](https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md). For commercial use, contact the author directly.

### 🔌 3.1 Backend Selection

At the top of the New Search tab, a **backend selector** lets you choose how DIA-NN runs:

- **Local (Embedded):** DIA-NN is installed inside the DE-LIMP container (Docker Compose deployment). Configure threads and output directory — no other setup needed.
- **Local (Docker):** DIA-NN runs in a separate Docker container. Configure CPU/memory limits via sliders.
- **HPC (SSH/SLURM):** Submit to a SLURM cluster. Choose between local (on-cluster) or remote (SSH) connection mode.

#### SSH Settings (Remote Mode)
When SSH mode is selected, a connection panel appears in the SLURM Resources section:

| Setting | Description |
| :--- | :--- |
| **Hostname** | The HPC login node address (e.g., `hive.hpc.university.edu`) |
| **Username** | Your HPC username |
| **Port** | SSH port (default: 22) |
| **SSH Key Path** | Path to your private key file (e.g., `~/.ssh/id_rsa`) |

- **Key-based authentication only** — no password entry or storage. Your SSH key must not require a passphrase.
- Click **"🔗 Test Connection"** to validate the SSH connection and locate SLURM binaries on the remote system.
- On success, the app caches the full path to `sbatch` (e.g., `/cvmfs/.../slurm/bin/sbatch`) so that all subsequent operations are fast (no login shell overhead).
- If the test reports "sbatch not found," the app automatically probes common HPC paths (`/usr/bin`, `/usr/local/bin`, spack, module directories) to find it.

### 📁 3.2 File Configuration (Panel 1)

The first panel configures your input files: raw data, FASTA database, and optional spectral library.

#### Raw Data Directory
- **Local mode:** Use the file browser to select the directory containing your raw data files.
- **SSH mode:** Click the **"Browse"** button to open the SSH File Browser and navigate to your data directory visually, or type/paste the remote path. Then click **"🔍 Scan Files"**.
- The scan detects mass spectrometry files and displays them with sizes:
  - `.d` directories (Bruker timsTOF)
  - `.raw` files (Thermo)
  - `.mzML` files (open format)
  - `.wiff` files (SCIEX)

#### FASTA Database
Four sources are available:

**📥 Download from UniProt:**
1. Type an organism name (e.g., "Homo sapiens", "Mus musculus") in the search box
2. Select a proteome from the results dropdown
3. Choose the content type:
   - **One protein per gene** (recommended) — canonical isoform only, smallest and cleanest database
   - **Canonical** — all reviewed canonical sequences
   - **Canonical + isoform** — includes splice variants
4. Click **"Download"** — the FASTA is downloaded to the HPC working directory (uploaded via SCP in SSH mode)

**📥 Download from NCBI (NEW in v3.7):**
1. Type an organism name in the search box
2. Select a proteome from the NCBI Datasets results
3. Click **"Download"** — downloads the RefSeq protein FASTA
4. **Gene symbol mapping**: NCBI RefSeq accessions (XP_, NP_, WP_) don't have embedded gene names like UniProt. DE-LIMP automatically runs a batch E-utilities lookup to map accessions to gene symbols. The gene map TSV is cached alongside the FASTA and auto-downloaded via SSH for Docker users.
5. Best for: non-model organisms, organisms with better NCBI than UniProt coverage

**📂 Pre-staged on server:**
- A dropdown of FASTA files already available on the cluster (pre-downloaded to a shared location at `/quobyte/proteomics-grp/de-limp/fasta/`)
- Fastest option for commonly used organisms

**📄 Browse / enter path:**
- **Local mode:** Use the file browser to locate any `.fasta` or `.fa` file
- **SSH mode:** Click the **"Browse"** button to open the SSH File Browser (see [Section 3.7](#-37-ssh-file-browser)), or type/paste the full remote path

#### Contaminant Library
Applies to all FASTA sources. Select from 6 curated contaminant libraries from [HaoGroup-ProtContLib](https://github.com/HaoGroup-ProtContLib):

| Library | Use Case |
| :--- | :--- |
| **Universal** (default) | General-purpose, covers common lab contaminants |
| **Cell Culture** | Optimized for cell line experiments |
| **Mouse Tissue** | Includes mouse-specific environmental contaminants |
| **Rat Tissue** | Includes rat-specific environmental contaminants |
| **Neuron Culture** | Specialized for neuronal cell culture experiments |
| **Stem Cell Culture** | Specialized for stem cell experiments |

The selected contaminant library is passed as a separate `--fasta` flag to DIA-NN, ensuring contaminant proteins are properly identified and can be filtered downstream.

#### Spectral Library (Optional)
- For **library-based** search mode only
- Browse or enter the path to a `.tsv` or `.speclib` spectral library file
- When omitted, DIA-NN runs in library-free mode (generates its own in silico library)

### ⚙️ 3.3 Search Settings (Panel 2)

The second panel configures DIA-NN analysis parameters.

#### Search Mode
| Mode | Description |
| :--- | :--- |
| **Library-free** (default) | DIA-NN generates an in silico spectral library from the FASTA. Best for most experiments. |
| **Library-based** | Uses a provided spectral library for peptide identification. Requires a spectral library file in Panel 1. |
| **Phosphoproteomics** | Auto-configures phospho-specific settings (see below). |

**Phosphoproteomics mode** automatically sets:
- STY phosphorylation variable modification (`UniMod:21` on S, T, Y)
- Maximum 3 variable modifications per peptide
- 2 missed cleavages
- `--phospho-output` flag (generates phosphosite-level output)
- `--report-lib-info` flag (reports library information for site localization)

#### Basic Settings
| Setting | Default | Description |
| :--- | :--- | :--- |
| **Enzyme** | Trypsin/P | Digestion enzyme for in silico digest |
| **Missed cleavages** | 1 | Maximum allowed missed cleavage sites |
| **Mass accuracy** | Auto | MS2 mass accuracy in ppm; "Auto" lets DIA-NN optimize |
| **MS1 mass accuracy** | Auto | MS1 mass accuracy in ppm; "Auto" lets DIA-NN optimize |
| **Max variable mods** | 2 | Maximum variable modifications per peptide |

#### Variable Modifications
| Modification | Default | Description |
| :--- | :--- | :--- |
| **Met oxidation** | ✅ On | Oxidation of methionine (UniMod:35) |
| **N-term acetylation** | Off | Acetylation of protein N-terminus (UniMod:1) |
| **Custom modifications** | — | Add any DIA-NN-compatible modification string |

#### Advanced Settings
| Setting | Default | Description |
| :--- | :--- | :--- |
| **FDR** | 0.01 | False discovery rate threshold (1%) |
| **Scan window** | Auto | Number of scans for chromatographic peak detection |
| **Peptide length** | 7–30 | Min and max peptide length for in silico digest |
| **Precursor m/z** | 300–1800 | Precursor mass-to-charge range |
| **MBR (Match Between Runs)** | ✅ On | Transfer identifications between runs |
| **RT profiling** | ✅ On | Builds a retention time model from your data to improve peptide scoring and quantification accuracy |
| **Normalization** | On | DIA-NN's built-in RT-dependent normalization (DE-LIMP applies DPC-CN on top of this during the pipeline step) |

### 🖥️ 3.4 SLURM Resources (Panel 3)

The third panel configures compute resources for the SLURM job.

| Setting | Default | Description |
| :--- | :--- | :--- |
| **CPUs** | 8 | Number of CPU cores for the DIA-NN search |
| **Memory (GB)** | 64 | RAM allocation (adjust based on FASTA size and file count) |
| **Time limit (hours)** | 24 | Maximum walltime before the job is killed |
| **Partition** | — | SLURM partition/queue to submit to |
| **Account** | — | SLURM account for resource billing |

In SSH mode, the SSH connection panel (hostname, username, port, key path, Test Connection button) also appears in this section.

### 📂 3.5 SSH File Browser (NEW in v3.7)

When using **Remote (SSH)** mode, the "Browse" buttons next to the raw data directory and FASTA path inputs open a visual file browser that lets you navigate your HPC file system without typing paths manually.

**Features:**
- **Clickable breadcrumbs**: Navigate the path hierarchy by clicking any segment in the breadcrumb trail
- **Up / Home buttons**: Go up one directory or jump to your home directory
- **Color-coded entries**: Folders (blue), data files like `.d`, `.raw`, `.parquet` (green), other files (grey)
- **File type filtering**: The browser shows relevant file types based on context (e.g., only `.fasta`/`.fa` when browsing for FASTA files, `.parquet` when loading results)
- **Double-click to enter**: Double-click folders to navigate into them; click "Select" to choose the current directory or file

**Performance notes:**
- The browser uses specific subdirectory roots configured for your HPC (`DELIMP_EXTRA_ROOTS` env var) rather than scanning the entire filesystem
- Directories with thousands of entries are paginated for responsiveness

### 📋 3.6 Job Queue

After clicking **"🚀 Submit Search"**, the job enters the **Job Queue** at the bottom of the New Search tab. You can submit multiple jobs and continue using the rest of DE-LIMP — searches are fully non-blocking.

#### Job Information
Each job in the queue displays:
- **Name:** A descriptive job name
- **SLURM Job ID:** The cluster-assigned job identifier
- **Status badge:** Color-coded status indicator
  - 🟡 **Queued** — Waiting in the SLURM queue
  - 🔵 **Running** — Actively processing on cluster nodes
  - 🟢 **Completed** — Finished successfully
  - 🔴 **Failed** — Exited with an error
  - ⚪ **Cancelled** — Manually cancelled by user
  - ❓ **Unknown** — Status could not be determined (e.g., SLURM purged the record)
- **Elapsed time:** How long the job has been running or total runtime
- **File count:** Number of raw data files in the search

#### Job Actions
| Button | When Available | Action |
| :--- | :--- | :--- |
| **📄 Log** | Always | View the stdout/stderr output from the SLURM job (reads from `{output_dir}/logs/`) |
| **❌ Cancel** | Queued or Running | Send `scancel` to terminate the job on the cluster |
| **📂 Load** | Completed | Download results (SCP for SSH mode) and load into DE-LIMP pipeline |
| **🔄 Refresh** | Unknown status | Re-query SLURM via `sacct` to update the job status |

- **"🔄 Refresh All"** button appears when any jobs have unknown status, refreshing all job statuses at once.
- **Auto-load:** When enabled, completed jobs automatically download results and load them into the pipeline — no manual "Load" click needed. For SSH mode, results are transferred via SCP.
- **Persistence:** The job queue is saved to `~/.delimp_job_queue.rds` and persists across app restarts. Restarting DE-LIMP restores your full job history.
- **Log file organization** (v3.2+): All SLURM `.out` / `.err` files are now stored in a `logs/` subdirectory within the output directory, keeping your results folder clean. The "Log" button checks this location first, with automatic fallback to the old root-level location for backward compatibility.

### 📝 3.7 Search Settings in Methodology

When results are loaded from a DIA-NN search (either via "Load" button or auto-load), DE-LIMP automatically records the search configuration in the **Methodology** tab.

A new **"0. DIA-NN DATABASE SEARCH"** section appears at the top of the methodology, documenting:
- **Raw files:** Count and file type (e.g., "24 Bruker .d files")
- **DIA-NN version:** The version installed on the cluster
- **Search mode:** Library-free, Library-based, or Phosphoproteomics
- **FASTA databases:** Primary database name and source
- **Contaminant library:** Which HaoGroup-ProtContLib contaminant FASTA was used
- **Enzyme:** Human-readable name (e.g., "Trypsin/P" instead of the DIA-NN flag)
- **Modifications:** Human-readable names (e.g., "Methionine oxidation, N-terminal acetylation")
- **FDR:** False discovery rate threshold
- **MBR:** Whether Match Between Runs was enabled
- **Mass accuracy:** Manual values or "auto-determined by DIA-NN"
- **SLURM resources:** CPUs, memory, and time limit used
- **Log files:** Location of SLURM log files (`{output_dir}/logs/`)

This section is **publication-ready** — it uses human-readable names for all parameters and follows standard methods section conventions. The same information is also logged in the **reproducibility code log** for programmatic access.

---

## 4. Deep Dive: The Data Overview & Grid View

### 📊 Data Overview
This is your landing page with 5 sub-tabs:
* **Assign Groups & Run** — Configure experimental groups and run the analysis pipeline
* **Signal Distribution** — Visualizes the dynamic range; automatically colors by DE status with synchronized comparison selector. Checkbox to overlay contaminant proteins in orange (v3.7).
* **Dataset Summary** — QC statistics and DE protein counts per comparison with directional arrows
* **Replicate Consistency** — Average precursor and protein counts per group
* **Expression Grid** — Heatmap-style table with UniProt linking, click-to-plot, and contaminant highlighting (pink/red rows for `Cont_` proteins)
* **Contaminant Analysis** (NEW in v3.7) — Summary cards (contaminant count, % of total, median intensity ratio vs endogenous, keratin count), per-sample stacked bar chart, top contaminants table with keratin flagging, and heatmap of top 20 contaminants by median intensity. Only visible when contaminant proteins (`Cont_` prefix) are detected.
* **Data Explorer** (NEW in v3.7) — Two panels for exploring data without requiring DE analysis:
  * **Abundance Profiles**: Proteins split into intensity quartiles (Q1=highest to Q4=lowest). Heatmap shows top 10 per quartile, colored by per-sample quartile assignment. Proteins that shift 2+ quartiles across samples are flagged as "Variable" — potentially biologically interesting even without replicates.
  * **Sample-Sample Scatter**: Pick any two samples and compare protein intensities. Identity line shows expected correlation. Outliers (>4-fold difference) are labeled with gene names. Shows Pearson correlation, protein count, and number of outliers. Contaminants shown as orange triangles.
* **AI Summary** — Generate AI-powered analysis summaries that analyze all contrasts simultaneously (requires Gemini API key); includes **"Export Report"** for standalone HTML and **"Export for Claude"** for a comprehensive .zip archive (see [Section 8](#8--ai-powered-analysis--export))

### 🔬 The Grid View (New!)
Click the green **"Open Grid View"** button to open the deep-dive table.

#### **Key Features:**
1.  **Bi-Directional Filtering:**
    * If you select proteins in the **Volcano Plot** (DE Dashboard) or if the **AI** selects interesting proteins, the Grid View automatically filters to show *only those proteins*.
    * Click **"Show All / Clear Selection"** in the footer to reset the view.
2.  **Compact Headers:**
    * Columns are labeled with **Run Numbers** (1, 2, 3...) to save space.
    * **Hover** your mouse over a number to see the full **File Name** and **Group**.
    * Headers are **color-coded** by Experimental Group (refer to the Legend at the top).
3.  **Heatmap Coloring:** Cell values (Log2 Intensity) are colored Blue (Low) to Red (High) for identifying patterns at a glance.
4.  **UniProt Integration:** Click any **Protein ID** to open its official UniProt page in a new tab.
5.  **Click-to-Plot:** Click any row in the table to instantly open a **Violin Plot** showing that specific protein's expression across all samples.
6.  **Smart Export:** Click **"Export Full Table"** to download the data as a CSV. The export will use the **Full Filenames** in the header (not the Run Numbers) for publication use.

---

## 5. Visualizing Results

### 📉 DE Dashboard
The DE Dashboard is organized into **four sub-tabs** for a cleaner workflow:

* **Current Comparison Display:** A prominent header banner at the top shows which comparison you're viewing (e.g., "Evosep - Affinisep"). This updates automatically when you change the comparison dropdown.

#### Volcano Sub-tab
* **Volcano Plot:** The volcano plot is the primary way to identify differentially expressed proteins. Each dot represents one protein: the x-axis shows how much the protein changed between conditions (log2 fold change), and the y-axis shows how statistically confident that change is (-log10 P-value). The most interesting candidates -- proteins with large, significant changes -- appear in the upper-left and upper-right corners (colored red). Interactive! Click points to select them. Box-select multiple points to analyze a cluster.
    * **Y-axis:** Shows -log10(raw P-Value) following proteomics best practices — this gives the classic volcano spread shape
    * **Coloring:** Red points indicate FDR-corrected significance (adj.P.Val < 0.05) — logFC vertical lines are visual guides only and do not gate coloring
    * **Threshold line:** The horizontal dashed line is drawn at the raw P.Value corresponding to the FDR boundary (adj.P.Val = 0.05), so the line and coloring always agree
    * **DE protein count:** The info box shows the total number of significant proteins with directional breakdown (e.g., "78 DE proteins (45 up, 33 down)")
    * **Default logFC cutoff:** 0.6 (~1.5-fold change) -- adjustable via the sidebar slider. The slider is in log2 units: 0.6 = ~1.5-fold, 1.0 = 2-fold, 2.0 = 4-fold. The vertical lines are visual guides only and do not gate significance coloring. **Choosing a cutoff**: For dramatic perturbations (knockout vs wild-type), many proteins change >2-fold, so 1.0 is reasonable. For subtle treatments (low-dose drug), even 1.3-fold changes may be biologically meaningful -- try 0.4. When in doubt, start with the default (0.6)
    * **Selection:** Single-click for one protein, box-select for multiple proteins
    * *Sync:* Selecting points here updates the Results Table and the AI context
* **Heatmap:** Displayed directly below the volcano plot. Automatically scales and clusters the top 50 significant proteins (or your specific selection).

#### Results Table Sub-tab
* **DE Results Table:** Shows both raw P-values and FDR-adjusted P-values for transparency, plus an **Avg CV (%)** column for each protein
* **Violin Plots:** Select one or more proteins and click **"📊 Violin Plot"** button
    * Multi-protein support: View multiple proteins in a 2-column grid layout
    * Individual scales: Each protein gets its own Y-axis for better visualization
    * Dynamic height: Adjusts based on number of selected proteins
* **XICs Button:** Click "📈 XICs" to inspect fragment chromatograms (local/HPC only)

#### PCA Sub-tab
* **PCA Plot:** Interactive scatter plot of samples in principal component space
    * **Color by:** Group, Batch, or any covariate column
    * **Axis selection:** Choose which principal components to display (PC1 vs PC2, etc.)
    * Helps identify sample clustering, outliers, and batch effects

#### CV Analysis Sub-tab
The **Coefficient of Variation (CV)** measures how reproducible a protein's measurement is across replicates -- it is the standard deviation divided by the mean, expressed as a percentage. A low CV (e.g., < 20%) means the protein was measured consistently, so you can be more confident that its fold change is real rather than noise. This tab helps you identify the most robust DE findings.

* **CV Analysis scatter plot:** Interactive Plotly scatter plot showing logFC vs Average CV for all significant proteins, color-coded by CV category (Excellent < 10%, Good 10-20%, Moderate 20-30%, High > 30%)
* **Summary stats subtitle:** Per-group median CV and percentage of proteins below 20% CV displayed directly in the plot subtitle
* **CSV export:** Download the full CV analysis data for all significant proteins

### 📐 QC Sample Metrics & Plots
* **Sample Metrics:** A single faceted trend plot showing four key per-run quality metrics stacked vertically:
    * **Precursors:** Number of peptide precursors identified at your Q-value cutoff
    * **Proteins:** Number of protein groups quantified per run
    * **MS1 Signal:** Overall MS1 intensity per run
    * **Data Completeness (%):** Percentage of precursors detected (non-missing) per sample in the raw expression matrix — shown as dots (not bars) to zoom into the relevant range
    * **LOESS Trendline:** A black smoothed trend line on each facet makes injection drift immediately visible — a flat line means stable performance, a downward slope suggests instrument degradation
    * **Group Average Lines:** Dashed horizontal lines show the mean for each experimental group
    * **Sort Order:** Toggle between Run Order (spot acquisition-time drift) and Group (compare conditions side by side)
    * **Fullscreen View:** Click **"🔍 Fullscreen"** to open the plot in a large modal for detailed inspection
* **MDS Plot:** A multidimensional scaling plot to visualize how samples cluster. (Good samples should cluster by Group).

### 📋 Reproducibility & Code Export
DE-LIMP automatically logs every analysis step for complete reproducibility.

**Features:**
* **Automatic Logging:** Every action (upload, pipeline run, contrast change, GSEA) is recorded with timestamps
* **Export Code:** Navigate to **Output > Methods & Code > R Code Log** tab
* **Download Button:** Click **"📥 Download Reproducibility Log"** to save as a timestamped `.R` file
* **Complete Script:** The exported file includes:
  * All analysis steps in executable R code
  * Session info (R version, package versions)
  * Group assignments and model formulas
  * Parameter settings (Q-value cutoffs, covariates)
* **Publication Ready:** Use the exported code to reproduce your analysis or include in Methods sections

**Methodology Summary:**
* View detailed methodology in the **Output > Methods & Code > Methods Summary** tab
* Includes citations for limpa, limma, and DIA-NN
* Explains normalization (DPC-CN), quantification (maxLFQ), and statistics (empirical Bayes)

### 📋 Export Data
The **Export Data** tab (under the Output dropdown) provides one-click access to download your results:

* **Results CSV:** Full DE results table for the selected comparison, including protein IDs, gene symbols, logFC, P-values, adjusted P-values, and per-sample expression values
* **CV Analysis CSV:** Coefficient of variation data for significant proteins across experimental groups

### 📈 XIC Chromatogram Viewer (Local/HPC Only)

The XIC Viewer lets you inspect fragment-level chromatograms for differentially expressed proteins, providing visual validation of quantification quality.

> **Note:** XIC files are generated by DIA-NN alongside the main report and are typically too large for cloud deployment. This feature is available for local and HPC installations only.

> **Hugging Face Users:** The XIC Viewer is not available on the hosted web version. The sidebar section is replaced with a link to download DE-LIMP for local use. To access chromatogram inspection, [download DE-LIMP.R from GitHub](https://github.com/bsphinney/DE-LIMP) and run locally or on your HPC cluster.

#### Setup
1. **Automatic Detection:** When you upload a DIA-NN report, the app automatically checks for a `_xic` directory in the working directory (e.g., `report_xic/` alongside `report.parquet`). If found, the XIC directory path is pre-filled.
2. **Manual Path:** If auto-detection doesn't find your files, paste the path to your XIC directory in the sidebar under **"5. XIC Viewer"** and click **"Load XICs"**.
   - You can also paste the path to the report `.parquet` file — the app will derive the `_xic` directory automatically.
3. The status badge shows the number of XIC files loaded, the detected DIA-NN version (1.x or 2.x), and whether ion mobility data is available.

#### Viewing Chromatograms
1. **Select a protein** in the DE Dashboard (volcano plot click or table row selection).
2. Click the **"📈 XICs"** button in the DE Dashboard results table header (or in the Grid View modal).
3. The XIC modal opens with interactive Plotly chromatograms.

#### Display Controls
- **Display Mode:**
  - *Facet by sample* — Each panel shows one sample with all fragment ions overlaid (color = fragment)
  - *Facet by fragment* — Each panel shows one fragment ion with all samples overlaid (color = group)
  - *Intensity alignment* — Spectronaut-style stacked bar chart showing relative fragment ion proportions per sample. Bars are ordered by experimental group with dashed separators. Automatic inconsistency detection flags samples where fragment ratios deviate significantly (> mean + 2 SD), with green (all consistent) or amber (flagged samples) guidance banners. Tooltips include AUC, proportion, deviation score, and cosine similarity.
- **Show MS1 (split axis):** When checked, the plot splits into two rows:
  - Top row: **MS1 precursor** signal (often much more intense)
  - Bottom row: **MS2 fragment** ions
  - Each row has its own y-axis, preventing the MS1 signal from squishing fragment peaks
- **Precursor Selector:** Choose a specific precursor or view all (top 6 shown for large proteins)
- **Group Filter:** Focus on a specific experimental group
- **Ion Mobility:** When timsTOF/PASEF mobilogram data is detected, a toggle appears with a prominent blue banner indicating ion mobility mode

#### Navigation
- Use **Prev/Next** buttons to step through significant DE proteins
- **Download** button exports the current view as PNG (14×10 inches, 150 DPI)

#### Info Panel
Below the plot, the info panel shows:
- Number of precursors and fragments
- Retention time range
- DE statistics (log2 fold-change and adjusted p-value) for the current comparison

---

### 🧬 Gene Set Enrichment (GSEA)

Gene Set Enrichment Analysis (GSEA) answers the question: "Are my differentially expressed proteins enriched in known biological pathways or functions?" Instead of interpreting proteins one by one, GSEA groups them into predefined sets and tests whether any set is overrepresented. Choose the database that matches your question:
- **Biological Process (BP):** What cellular processes are affected? (e.g., "cell division," "immune response") -- best general-purpose choice
- **Molecular Function (MF):** What molecular activities are changed? (e.g., "kinase activity," "DNA binding")
- **Cellular Component (CC):** Where in the cell are the changes? (e.g., "mitochondria," "nucleus")
- **KEGG:** Which metabolic or signaling pathways are involved? (e.g., "glycolysis," "MAPK signaling")

1.  Select a **database** from the ontology selector dropdown: Biological Process (BP), Molecular Function (MF), Cellular Component (CC), or KEGG pathways.
2.  Click **"Run GSEA"** for the current contrast.
3.  **Automatic organism detection**: The app queries the UniProt API using your protein IDs to determine the correct organism database. For human data this is instant; for non-human data (mouse, rat, etc.) the API lookup runs automatically.
4.  **Per-ontology caching**: Results are cached separately for each database and contrast combination. Switch between BP, MF, CC, and KEGG without re-running the analysis.
5.  **Contrast indicator**: A banner shows which contrast is active. If you change the comparison, a stale-results warning appears prompting you to re-run.
6.  View results as Dot Plots, Enrichment Maps (networks), Ridgeplots, or browse the full Results Table.

### 🔬 Phosphoproteomics

Phosphoproteomics studies how proteins are regulated by phosphorylation -- a chemical modification where a phosphate group is added to specific amino acids (Serine, Threonine, or Tyrosine). Phosphorylation acts as an on/off switch for many cellular processes, including signaling, growth, and metabolism. "Site-level" analysis means DE-LIMP tests each individual phosphorylation site independently, so you can see exactly which positions on which proteins change between your conditions.

The phosphoproteomics module provides site-level analysis of phosphorylation data, available when phospho-enriched data is detected.

#### Auto-Detection
- On file upload, the app scans for phospho modifications (`UniMod:21`) and displays a detection banner if phospho data is present
- The **Phosphoproteomics** tab appears automatically when phospho data is detected

#### Input Paths
- **Site matrix upload** (recommended): Upload a DIA-NN 1.9+ `site_matrix_*.parquet` file directly
- **Parsed from report**: The app extracts phosphosites from `Modified.Sequence` columns in your main report file

#### Phase 1 — Site-Level DE
- **Phospho Volcano Plot**: Interactive volcano plot for phosphosite-level differential expression
- **Site Table**: Full results table with site ID, protein, gene, residue, position, fold-change, and significance
- **Residue Distribution**: Breakdown of Serine/Threonine/Tyrosine phosphorylation frequencies
- **QC Completeness**: Missingness analysis across sites and samples with filtering thresholds

#### Phase 2 — Kinase Activity & Motifs
- **KSEA (Kinase-Substrate Enrichment Analysis)**: Infers upstream kinase activity from phosphosite fold-changes using PhosphoSitePlus and NetworKIN databases
- **Sequence Logo Analysis**: Visualizes enriched amino acid motifs around significant phosphosites (up-regulated vs. down-regulated)

#### Phase 3 — Advanced
- **Protein-level abundance correction**: Subtracts protein-level logFC from phosphosite logFC to isolate changes in phosphorylation stoichiometry (requires matched total proteome and phospho-enriched samples)
- **AI context**: Phosphosite DE results and kinase activities are included in the Gemini chat context when phospho analysis is active

#### References & Methods
The phosphoproteomics module is grounded in the following literature:

**Core Data Processing:**
- **DIA-NN site-level reporting**: DIA-NN 1.9+ natively produces site quantification matrices with localization confidence scores. [github.com/vdemichev/DiaNN](https://github.com/vdemichev/DiaNN)
- Pham TV, Henneman AA, Truong NX, Jimenez CR (2024). "msproteomics sitereport: reporting DIA-MS phosphoproteomics experiments at site level with ease." *Bioinformatics* 40(7):btae432. [doi:10.1093/bioinformatics/btae432](https://doi.org/10.1093/bioinformatics/btae432)
- Bekker-Jensen DB et al. (2020). "Rapid and site-specific deep phosphoproteome profiling by data-independent acquisition without the need for spectral libraries." *Nat Commun* 11:787. [doi:10.1038/s41467-020-14609-1](https://doi.org/10.1038/s41467-020-14609-1)
- Muneer A et al. (2025). "Advancements in Global Phosphoproteomics Profiling: Overcoming Challenges in Sensitivity and Quantification." *PROTEOMICS* 2400087. [doi:10.1002/pmic.202400087](https://doi.org/10.1002/pmic.202400087)

**Kinase Activity Inference:**
- Wiredja DD, Koyutürk M, Chance MR (2017). "The KSEA App: a web-based tool for kinase activity inference from quantitative phosphoproteomics." *Bioinformatics* 33(21):3489–3491. [doi:10.1093/bioinformatics/btx687](https://doi.org/10.1093/bioinformatics/btx687)
- Piersma SR et al. (2024). "Inferring kinase activity from phosphoproteomic data: Tool comparison and recent applications." *Mass Spectrometry Reviews* 43:552–571. [doi:10.1002/mas.21808](https://doi.org/10.1002/mas.21808)
- Kim HJ et al. (2021). "PhosR enables processing and functional analysis of phosphoproteomic data." *Cell Reports* 34(8):108771. [doi:10.1016/j.celrep.2021.108771](https://doi.org/10.1016/j.celrep.2021.108771)

**Motif & Sequence Visualization:**
- Wagih O (2017). "ggseqlogo: a versatile R package for drawing sequence logos." *Bioinformatics* 33(22):3645–3647. [doi:10.1093/bioinformatics/btx469](https://doi.org/10.1093/bioinformatics/btx469)

**DIA Phosphoproteomics Workflows:**
- Skowronek P et al. (2022). "Rapid and In-Depth Coverage of the (Phospho-)Proteome With Deep Libraries and Optimal Window Design for dia-PASEF." *MCP* 21(9):100277.
- Kitata RB et al. (2021). "DIA-based global phosphoproteomics system using hybrid spectral libraries." *Nat Commun* 12:2539. [doi:10.1038/s41467-021-22759-z](https://doi.org/10.1038/s41467-021-22759-z)
- Roßmann K et al. (2024). "Data-Independent Acquisition: A Milestone and Prospect in Clinical Mass Spectrometry–Based Proteomics." *MCP* 23(7):100800. [doi:10.1016/S1535-9476(24)00090-2](https://doi.org/10.1016/S1535-9476(24)00090-2)

**Normalization:**
- Protein-level abundance correction isolates phosphorylation stoichiometry from total protein changes (Piersma 2024; PhosR documentation).
- Tail-based imputation follows the Perseus-style approach: downshifted normal distribution (mean − 1.8 SD, width 0.3 SD) for missing values assumed to be below detection limit (Tyanova et al. 2016).

### 🎓 Education & Resources
Click the **"Education"** tab to access embedded proteomics training materials without leaving the app.
* **UC Davis Proteomics Videos:** Latest YouTube content auto-updates dynamically
* **Hands-On Proteomics Short Course:** Information about UC Davis summer training
* **Core Facility Resources:** Direct links to proteomics.ucdavis.edu
* **Google NotebookLM:** Explore the key citations behind DE-LIMP's methodology (limpa, limma, DIA-NN)
* **Proteomics News Blog:** Stay updated with the latest in the field

### ℹ️ About
Click the **"About"** tab in the navbar to view project information and community activity:
* **Version display**: Shows the current app version, read from the `VERSION` file
* **Community stats cards**: GitHub stars, forks, unique visitors (14-day window), and unique clones — updated daily by a GitHub Actions workflow
* **Trend sparklines**: Interactive 14-day sparkline charts for unique visitors and unique clones, making adoption trends visible at a glance
* **GitHub Discussions feed**: The 10 most recently updated discussions with title, category, author, comment count, and direct link — engage with the community without leaving the app
* **Quick links**: One-click access to the GitHub repository, Hugging Face Space, documentation site, and GitHub Discussions
* **Freshness indicator**: Shows when the stats were last updated

> Community stats are collected by the `track-stats.yml` GitHub Actions workflow, which runs daily and saves data to `stats/community_stats.json`. Stats will appear after the first workflow run.

### 💾 Save & Load Analysis Sessions
DE-LIMP can save your entire analysis state for later use.

**To Save a Session:**
1. Complete your analysis (upload data, run pipeline, explore results)
2. Click the **"Save"** button in the sidebar (below the accordion panels)
3. Choose a filename and location
4. The file saves as `.rds`
5. The session file includes:
   * Raw data
   * Processed data (normalized, quantified)
   * Statistical results (limma fit object)
   * All group assignments and settings
   * App version tag for compatibility tracking

**To Load a Session:**
1. Click the **"Load"** button in the sidebar
2. Select a previously saved `.rds` file
3. The entire analysis state is restored instantly
4. Continue where you left off without re-running the pipeline

**Use Cases:**
* Share analyses with collaborators
* Archive completed projects
* Test different parameters without re-uploading data
* Quick access to previous experiments

---

## 6. Multi-Omics MOFA2

The **Multi-Omics MOFA2** tab provides unsupervised multi-omics integration using MOFA2 (Multi-Omics Factor Analysis). It discovers latent factors -- hidden patterns that explain variation across your datasets.

**When should you use MOFA2?** Use it when you have two or more types of data measured on the same samples and want to understand which patterns are shared versus unique to each data type. For example, if you ran both total proteomics and phosphoproteomics on the same samples, MOFA2 can separate protein abundance changes from true phosphorylation regulation. It is also useful for identifying hidden batch effects or finding biological signals that only become apparent when integrating multiple data types.

### When to Use MOFA2
- **Global + Phospho**: Separate protein abundance effects from true phospho-regulation changes
- **Multiple experiments**: Find what's shared vs unique between different measurements
- **QC discovery**: Identify hidden batch effects or sample outliers across views
- **Multi-omics**: Integrate proteomics with RNA-seq, metabolomics, or other -omics data

### 6.1 Loading Data Views

Each MOFA view is a features × samples matrix. You can load views from:

| Source | How |
| :--- | :--- |
| **Current DE pipeline** | View 1 auto-populates from your loaded proteomics data |
| **Phospho tab** | Click "Use Phospho Data" to add site-level data as a view |
| **File upload** | Upload CSV/TSV/Parquet (first column = feature IDs, remaining = samples) |
| **RDS import** | Upload DE-LIMP session files or limma objects — the smart parser extracts the expression matrix |
| **Example data** | Click "Mouse Brain (2-view)" or "TCGA Breast (3-view)" buttons |

- **2-6 views required** — use the "Add View" button to add more, "Remove" to delete
- **Sample matching** — the app automatically finds common samples across all views and reports overlap statistics

### 6.2 Training Parameters

| Parameter | Description | Default |
| :--- | :--- | :--- |
| **Number of Factors** | Latent factors to discover (auto or manual) | Auto (up to 15) |
| **Convergence Mode** | Fast (~500 iter), Medium (~1000), Thorough (~5000) | Medium |
| **Scale Views** | Equalize contribution of each view (recommended when views have different feature counts) | ON |
| **Min Variance** | Drop factors explaining less than this % of total variance | 1% |
| **Seed** | Random seed for reproducibility | 42 |

Click **"Train MOFA Model"** to start. Training runs in an isolated subprocess and typically takes 1-5 minutes depending on data size.

### 6.3 Results Tabs

After training completes, five results tabs appear:

1. **Variance Explained** — Heatmap showing % variance each factor explains per view. Factors loading heavily on one view indicate view-specific variation; factors loading similarly across views indicate shared biology.
2. **Factor Weights** — Bar chart of top N features driving each factor. Select a view and factor from the dropdowns. High |weight| = strong contributor.
3. **Sample Scores** — Scatter plot of samples in factor space, colored by experimental group. Clustering indicates shared biology captured by those factors.
4. **Top Features** — Sortable table ranking features by absolute weight across all views for a selected factor.
5. **Factor-DE Correlation** — Bar chart showing Pearson correlation between each factor's weights and DE log-fold-changes. Requires the DE pipeline to have been run first.

### 6.4 Example Datasets

Two built-in datasets for testing:

| Dataset | Button | Views | Samples | Groups |
| :--- | :--- | :--- | :--- | :--- |
| **Mouse Brain** | "Mouse Brain (2-view)" | Global Proteomics (10,333) + Phospho (89 sites) | 16 | F_PME, F_PSE, M_PME, M_PSE |
| **TCGA Breast Cancer** | "TCGA Breast (3-view)" | mRNA (~200) + miRNA (184) + Protein (142) | 150 | Basal, Her2, LumA |

---

## 7. 🏢 Core Facility Mode

Core Facility Mode transforms DE-LIMP into a managed proteomics analysis platform for core labs, adding job tracking, instrument QC monitoring, and automated report generation.

> **Note:** Core Facility Mode is optional and not visible unless explicitly activated. Standard users and Hugging Face deployments are unaffected.

### 7.1 Activation

Set the `DELIMP_CORE_DIR` environment variable to a directory containing a `staff.yml` configuration file:

```bash
# Example:
export DELIMP_CORE_DIR=/srv/delimp
# Then launch the app normally
shiny::runApp('.', port=3838)
```

The directory should contain:
- `staff.yml` — Staff member profiles with SSH/SLURM configuration
- `delimp.db` — SQLite database (auto-created on first run)
- `reports/` — Generated HTML reports (auto-created)
- `state/` — Saved analysis state files (auto-created)

### 7.2 Staff Configuration

The `staff.yml` file defines staff members and their HPC credentials:

```yaml
staff:
  - name: "Jane Smith"
    username: "jsmith"
    host: "hpc.university.edu"
    key_path: "~/.ssh/id_rsa"
    account: "proteomics_lab"
    partition: "high"
    lab: "Smith Lab"
```

Selecting a staff member from the dropdown auto-fills SSH host, username, key path, and SLURM account/partition — no manual entry needed.

### 7.3 Search DB Tab

The **Search DB** tab (under the Facility dropdown) provides a searchable history of all DIA-NN searches:

- **6 filters**: Free-text search, lab, status, staff, instrument, and LC method
- **Load Results**: Re-load results from any past search into the analysis pipeline
- **Generate Report**: Create a standalone HTML report from any completed search

### 7.4 Instrument QC Tab

The **Instrument QC** tab monitors instrument performance over time:

- **Trend plots**: Protein count, precursor count, and TIC per QC run
- **Control lines**: Rolling mean ± 2SD for anomaly detection
- **Instrument filter**: Select specific instruments to monitor
- **Date range**: Focus on a specific time period
- **Runs table**: Detailed metrics for each QC run

### 7.5 Report Generation

Click **"Generate Report"** on any completed search to produce a standalone HTML report:

- **Metadata header**: Title, lab, instrument, LC method, project, analyst, date
- **QC bracket**: Comparison with the nearest HeLa QC runs (before and after)
- **Volcano plots**: For each contrast in the analysis
- **DE statistics**: Protein counts by significance threshold
- **Top proteins table**: Most significant differentially expressed proteins
- **Normalization diagnostics**: Pre/post normalization signal distributions

Reports are saved to the `reports/` directory and recorded in the SQLite database for tracking.

---

## 8. 🤖 AI-Powered Analysis & Export

DE-LIMP offers two complementary AI pathways:

- **In-app AI (Google Gemini):** Quick questions and summaries powered by Google Gemini, right inside the app. This includes **AI Summary** (a one-click overview of all comparisons) and **Data Chat** (interactive Q&A about your data). Requires a free Gemini API key.
- **Export for External AI:** Download your complete analysis as a .zip to upload to Claude, ChatGPT, or other AI tools for deeper analysis, manuscript writing, or extended interpretation. No API key needed for the export itself.

### 8.1 Setup — Google Gemini API Key

A free API key from Google is required for all AI features (AI Summary, Data Chat, Auto-Analyze).

1.  Go to **[Google AI Studio](https://aistudio.google.com/)**.
2.  Sign in with your Google Account.
3.  Click **"Get API key"** in the top-left corner.
4.  Click **"Create API key"** (select "Create API key in new project" if prompted).
5.  Copy the key (starts with `AIza...`).
6.  Paste it into the **"Gemini API Key"** box in the DE-LIMP sidebar (under the AI Chat accordion section).
7.  (Optional) Change the Model Name to use a specific Gemini version (default: `gemini-3-flash-preview`).

> **Privacy:**
> - **AI Summary** sends only summary statistics to Gemini: protein names, logFC, adj.P.Val, CV metrics, and dataset dimensions. No raw expression values or sample identifiers are included.
> - **Data Chat** sends per-sample expression values for the top DE proteins and QC metrics that include run identifiers, giving Gemini richer context for interactive Q&A.
> - Neither feature sends file paths or server information.
> - Google retains API data for approximately 48 hours for abuse monitoring. If you are working with clinical or patient-derived data, consult your institutional data governance office before using any AI features.

### 8.2 AI Summary (Data Overview > AI Summary Sub-tab)

The AI Summary analyzes **all contrasts (pairwise comparisons between your experimental groups, e.g., Treatment vs. Control) simultaneously**, not just the currently selected comparison. This provides a global view of your experiment.

**What data is sent to Gemini:**
* Top differentially expressed proteins per comparison (gene names, logFC (log2 fold change -- a value of 1.0 means the protein doubled), adj.P.Val (p-value corrected for multiple testing))
* Cross-comparison biomarkers -- proteins that are significant in two or more contrasts
* CV-based stability metrics -- median coefficient of variation per group, percentage of proteins below 20% CV
* Dataset dimensions (number of proteins, samples, groups, contrasts)

**What is NOT sent:**
* Raw expression values or intensity matrices
* Sample file names or identifiers
* File paths or server information

**The AI generates:**
* Biological interpretation of the top DE proteins in each comparison
* Cross-comparison patterns -- proteins that change consistently across multiple contrasts
* Pathway and functional context for the findings
* Suggestions for follow-up experiments

#### AI Summary HTML Export

Click **"Export Report"** below the AI Summary to download a styled standalone HTML file:
* Gradient header with experiment metadata
* Full AI analysis with markdown formatting preserved
* Print-friendly CSS suitable for sharing with collaborators
* Self-contained -- no external dependencies, opens in any browser

### 8.3 Data Chat (AI Analysis Tab)

The **AI Analysis** tab provides an interactive conversational interface with Google Gemini, where the AI has full awareness of your dataset context.

#### How It Works

When you open the Data Chat, the app automatically sends Gemini:
* QC statistics (protein counts, precursor counts, data completeness per sample)
* Top differentially expressed proteins **for the currently selected comparison** -- change the comparison selector in the DE Dashboard to explore other contrasts with the AI
* Smart data scaling: sends 100-800 proteins depending on dataset size (smaller datasets send more complete results; larger datasets focus on the most significant)
* When phosphoproteomics analysis is active, the top 20 phosphosites and KSEA kinase activity results are automatically included

#### Using Data Chat

* **Ask questions** about your specific data:
    * *"Which group has the highest variance?"*
    * *"Are there any mitochondrial proteins upregulated?"*
    * *"What biological processes are enriched among the top hits?"*
    * *"Generate a figure caption for the volcano plot."*
    * *"Summarize the key findings for a lab meeting presentation."*
* **Auto-Analyze:** Click the **"Auto-Analyze"** button for a one-click comprehensive report. The AI generates a structured analysis covering data quality assessment, top differentially expressed proteins, and biological interpretation in approximately 30-60 seconds -- no manual prompting needed.

#### Interactive AI and Plot Connection

This is one of DE-LIMP's most powerful features -- the AI and the interactive plots are connected in both directions:

**User to AI (select proteins, then ask):**
1. Select proteins on the Volcano Plot (click or box-select) or in the Results Table
2. Ask: *"What are the functions of these selected proteins?"*
3. The AI receives the exact protein list you selected and responds with targeted analysis

**AI to User (AI highlights proteins in plots):**
1. The AI can suggest proteins using a special syntax: `[[SELECT: GAPDH; ENO1; PKM]]`
2. When the AI includes this in a response, DE-LIMP automatically highlights those proteins in the volcano plot and filters the Results Table
3. Example: Ask *"Show me glycolytic enzymes"* -- the AI identifies them and highlights them in your plots

> **Note:** You do not need to type `[[SELECT:]]` yourself. The AI automatically uses this format when it identifies proteins of interest, and DE-LIMP reads it to update your plots.

#### Save Chat History

Click **"Save Chat"** to download the full conversation as a plain text file. The export includes both your messages and all AI responses, with timestamps. Useful for documenting your analytical reasoning or sharing insights with collaborators.

### 8.4 Export for Claude (.zip Archive)

> **Works with any AI:** This export is optimized for Claude but works equally well with ChatGPT, Gemini, Copilot, or any AI assistant that accepts file uploads.

The **"Export for Claude"** button (on the **AI Summary** sub-tab under Data Overview) downloads a comprehensive multi-file package designed for deep analysis with Claude or other external AI systems. While the in-app AI features use Google Gemini, this export creates a portable dataset package optimized for extended conversation-based analysis.

**The .zip archive contains:**

| File | Contents |
| :--- | :--- |
| **`PROMPT.md`** | Full context document explaining the experimental design, statistical methodology, and how to interpret each file -- serves as an instruction manual for the AI |
| **`DE_Results_Full.csv`** | All proteins across all contrasts with logFC, P.Value, adj.P.Val, and expression values |
| **`Expression_Matrix.csv`** | Log2 expression values (rows = proteins, columns = samples) |
| **`QC_Metrics.csv`** | Per-sample quality control statistics (precursor counts, protein counts, MS1 signal, data completeness) |
| **`GSEA_Results.csv`** | Gene set enrichment results across all ontologies (included if GSEA has been run) |
| **`Phospho_DE_Results.csv`** | Site-level phospho differential expression results (included if phospho data is detected) |
| **`Session.rds`** | Full DE-LIMP session state -- can be reloaded into DE-LIMP to restore the exact analysis. **Note:** Contains all raw and processed data, so this file can be very large |
| **`Group_Assignments.csv`** | Sample-to-group mapping table |
| **`Analysis_Parameters.txt`** | Pipeline settings (Q-value cutoff, covariates, normalization method) |
| **`Methods_and_References.txt`** | Statistical methodology text suitable for a paper's Methods section, with citations |
| **`Reproducibility_Code.R`** | Complete R code log with timestamps for every analysis step |

**How to use it:**
1. Click **"Export for Claude"** on the AI Summary sub-tab to download the .zip file
2. Go to [claude.ai](https://claude.ai) (free tier available), [chatgpt.com](https://chatgpt.com), or another AI assistant
3. Start a new conversation and upload the .zip file (or individual files like `PROMPT.md` + the relevant CSVs)
4. Ask questions like *"Summarize the key biological findings"*, *"Help me write a results paragraph for my paper"*, or *"What pathways are most affected?"*

**Use cases:**
* In-depth biological interpretation beyond what the in-app chat provides
* Help writing a methods section or results narrative for a manuscript
* Compare your results against known biology or published datasets
* Generate publication-quality figure descriptions
* Explore specific pathways or protein families in detail

> **Note:** The Export for Claude package is for use with external AI tools (Claude, ChatGPT, etc.). It does not connect to or require any Anthropic API key. The in-app AI features (Sections 8.2 and 8.3) use Google Gemini.

### 8.5 Other Export Features

These additional export options are available from the **Output** dropdown in the navbar:

* **Export Data panel** -- One-click CSV downloads for DE Results and CV Analysis data
* **Reproducibility Code Log** -- Timestamped R script documenting every analysis step (download as `.R` file from the Methods & Code tab)
* **Methods Summary** -- Publication-ready methodology text with citations for limpa, limma, and DIA-NN
* **Session Save/Load** -- Save the full analysis state as `.rds` for later use or sharing (Save/Load buttons are in the sidebar)

---

## 9. Run Comparator

The **Run Comparator** lets you compare two analyses of the same dataset to understand how different tools, settings, or pipelines affect your results. This is essential for benchmarking, validating findings across platforms, and building confidence in your DE protein list.

### 9.1 When to Use It

- You ran the same samples through DE-LIMP twice with different settings (e.g., different normalization, mass accuracy, or FASTA database) and want to know what changed
- Your core facility processed samples with Spectronaut or FragPipe, and you want to compare against your DE-LIMP analysis
- You want to identify which DE proteins are robust (consistent across tools) vs. tool-dependent

> **Protein inference caveat:** Different tools may group shared peptides into different protein groups (e.g., tool A reports P12345 while tool B reports P12345;P67890 as a group). The comparator normalizes protein IDs to bare UniProt accessions, but some mismatches are unavoidable when tools resolve protein ambiguity differently. Proteins unique to one run in the Protein Universe tab may partly reflect these grouping differences rather than true identification failures.

### 9.2 Comparison Modes

| Mode | Run A | Run B | What You Need |
|------|-------|-------|---------------|
| **A: DE-LIMP vs DE-LIMP** | Current session or .rds file | .rds file | Two DE-LIMP session files (or one + current session) |
| **B: DE-LIMP vs Spectronaut** | Current session or .rds file | Spectronaut candidates.tsv export | Spectronaut "Candidates" export with IntensityPG columns |
| **C: DE-LIMP vs FragPipe** | Current session or .rds file | combined_protein.tsv | FragPipe combined_protein.tsv; optionally with FragPipe-Analyst DE stats |

### 9.3 Running a Comparison

1. Navigate to **Analysis > Run Comparator** in the navbar
2. Select the comparison mode (A, B, or C)
3. For **Run A**: Choose "Use current session" (if you have a loaded analysis) or upload a DE-LIMP `.rds` session file
4. For **Run B**: Upload the appropriate file for the chosen mode
5. Select the **contrast** to compare (must exist in both runs -- e.g., "Treatment - Control")
6. Click **Run Comparison**

**Important**: Before comparing, verify that both runs used the same MBR (match-between-runs) setting. MBR can add 10-30% more protein identifications, so comparing MBR-on vs MBR-off will produce large protein universe differences and many "Missing values" hypotheses that reflect the MBR setting rather than genuine analytical disagreement.

**Mode A bonus -- DIA-NN Log Upload**: In the collapsible "DIA-NN Log Files" section, you can optionally upload DIA-NN log files for each run. This enriches the Settings Diff with search-derived parameters like pg-level quantification mode, proteoform detection, library precursor counts, and which pipeline step produced the output. Useful when comparing a first-pass quant vs final assembly output.

### 9.4 The Four Diagnostic Layers

Results appear as sub-tabs, each building on the previous:

**Settings Diff** -- A side-by-side parameter table highlighting differences in red. Covers pipeline settings (normalization, imputation, covariates), DIA-NN search settings (mass accuracy, enzyme, scan window), and DIA-NN log-derived settings (if uploaded). Mismatched rows are highlighted for quick scanning.

**Protein Universe** -- A stacked bar chart showing how many proteins are shared, Run A-only, and Run B-only. Summary cards show total protein counts for each run and the overlap percentage. Large differences here often indicate different FASTA databases, different filtering thresholds, or dramatically different data completeness.

**Quantification** -- Three views of how protein-level intensities compare:
- *Correlation scatter*: Log2 intensity of each protein in Run A vs Run B. Pearson r and systematic bias (median offset) displayed.
- *Per-sample correlation*: Bar chart showing how well each sample's intensities agree between runs. Low-correlation samples may indicate normalization differences or batch effects.
- *Bias density*: Histogram of log2(Run A / Run B) per protein. A symmetric distribution centered at zero means no systematic bias. Shifts indicate normalization or quantification differences.

**DE Concordance** -- The core diagnostic:
- *3x3 concordance matrix*: Each protein classified as Up, Down, or Not Significant in each run. The 9 cells show counts and percentages. Ideally, most proteins fall on the diagonal (agree) rather than off-diagonal (disagree).
- *Volcano overlay*: Both runs' volcano plots superimposed, colored by concordance status.
- *Discordant protein table*: Every protein where the two runs disagree, with per-protein hypothesis explaining why.
- *Hypothesis distribution chart*: Bar chart showing the dominant causes of discordance at a glance.
- *Summary banner*: One-line overview with concordance rate, bias detection badge, and dominant cause badge.

### 9.5 Understanding Hypotheses

The hypothesis engine assigns one of 7 categories to each discordant protein. These are rule-based heuristics (not formal statistical tests) designed to guide your investigation -- they indicate the most likely explanation, not a definitive cause:

| Category | Meaning | Typical Cause |
|----------|---------|---------------|
| **Direction reversal** | One run says Up, the other says Down | Different normalization centering or peptide selection |
| **Normalization offset** | Same direction but one run crosses the significance threshold | One tool normalizes differently (e.g., DIA-NN + DPC-CN vs Spectronaut local regression), shifting all intensities up or down, which pushes borderline proteins above or below the significance threshold |
| **Variance estimation** | Similar fold changes but different significance | The two tools handle measurement noise differently. Limma borrows information across all proteins (empirical Bayes) to stabilize variance estimates; other tools may use per-protein variance only, leading to different p-values for the same fold change |
| **Missing values** | One run has fewer observations | Different MBR (match-between-runs) settings or imputation. MBR can add 10-30% more identifications -- comparing MBR-on vs MBR-off produces many hits in this category |
| **Peptide count** | Different number of supporting peptides | One tool used more peptide measurements for this protein. More peptides generally means a more stable estimate; the tool with fewer peptides may report a noisier fold change |
| **FC magnitude** | Fold change is larger in one run | Different protein quantification methods (DPC-Quant vs MaxLFQ) combine peptide measurements differently, producing different fold-change estimates for the same protein |
| **Borderline** | Both runs have the protein near the significance boundary | Not a true disagreement -- dichotomizing continuous p-values at a fixed threshold (0.05) inevitably creates disagreements for proteins near the boundary. Small perturbations in data processing push them across |

**"Borderline" is the most common and least concerning** -- it means the protein is close to adj.P.Val = 0.05 in both runs and small perturbations push it across the threshold. This is a fundamental limitation of significance thresholds, not a tool-specific problem. Focus your attention on Direction Reversals and Normalization Offsets.

### 9.6 What to Do with the Results

> **Note on concordance rates:** The reported concordance rate includes proteins that are non-significant (NS) in both runs. Since most proteins are NS, the base-rate concordance is naturally high. A 90% concordance rate does not mean 90% of your DE proteins agree -- focus on the 3x3 matrix and discordant protein count for a clearer picture.

- **Concordance rate >80%**: Typical when comparing the same tool with minor parameter changes. Your results are robust.
- **Concordance rate 60-80%**: Common across different tools (e.g., DE-LIMP vs Spectronaut). Focus on proteins that are consistent across both runs -- these are your highest-confidence hits.
- **Concordance rate <60%**: Investigate the dominant cause. Large protein universe differences suggest different FASTA databases or MBR settings. Many "Normalization offset" hits suggest the tools center intensities differently.
- **For publications**: Report "X proteins were significant in both analyses (Y% concordance); Z discordant proteins were predominantly borderline cases." Proteins consistent across tools are your strongest candidates for validation.
- **Which run to trust?** Neither is inherently "correct." If one run used more replicates, better normalization, or a more complete FASTA database, prefer its results. The comparator helps you understand *where* and *why* the tools disagree, so you can make an informed decision.

### 9.7 AI Analysis & Export

- **Gemini Analysis**: Click "Analyze with Gemini" on the AI Analysis sub-tab for a narrative interpretation. The prompt is tool-aware -- it includes context about structural differences between the compared tools.
- **MOFA2 Decomposition**: Click "Run MOFA2" to decompose the joint variance between runs into latent factors. Helps identify whether discordant proteins share hidden biological or technical patterns.
- **Claude ZIP Export**: Download a .zip optimized for deep analysis with Claude or ChatGPT. Includes settings diff, protein universe, DE results, discordant proteins with hypotheses, DIA-NN log parameters (if uploaded), and a structured prompt.

---

## 10. Chromatography QC

The **Chromatography QC** system lets you inspect TIC (Total Ion Current) traces from your timsTOF raw files *before* committing to a DIA-NN search. This catches common problems -- dead injections, sample carryover between runs, retention time shifts, and uneven sample loading -- that would otherwise waste hours of compute time.

> **Currently timsTOF only.** Thermo .raw TIC extraction is planned for a future release.

### 10.1 Extracting TIC Traces

1. **Scan your raw data directory** on the New Search tab (SSH or local)
2. After files are scanned, click the **"Extract TIC"** button that appears
3. The app reads the total ion signal from each raw instrument file to build a chromatographic profile for every run
4. For SSH mode, each file is downloaded temporarily, extracted locally, then cleaned up

### 10.2 QC Tab Views

After extraction, the **QC** tab in the navbar becomes visible with three views:

**Faceted View** -- Each run gets its own panel (4 columns). The blue dashed line is the median trace across all runs. Run traces are colored by diagnostic status: green (pass), yellow (warning), red (fail). Quickly spot runs that deviate from the cohort.

**Overlay View** -- All runs normalized to 0-1 intensity on a single axis. Excellent for spotting outlier shapes -- a run with a very different peak position or width stands out immediately.

**Metrics View** -- Horizontal AUC bar chart plus a detailed diagnostics table with columns: File Size, AUC, Peak RT, Gradient Width, Baseline Ratio, Late Signal, Shape Correlation (r), and Flags. Each flag links to a specific diagnostic:

| Diagnostic | Warning | Fail | What It Means |
|-----------|---------|------|---------------|
| Shape deviation | r < 0.95 | r < 0.90 | Run shape differs from the cohort median |
| RT shift | >2 MAD from median peak RT | >3 MAD | Retention time shifted (column degradation, gradient problem) |
| Loading anomaly | AUC >2x or <0.5x median | >3x or <0.3x | Much more or less sample loaded |
| File size outlier | >2x or <0.5x median | >3x or <0.3x | Acquisition anomaly |
| Late elution | -- | >15% signal in last 20% of gradient | Carryover or column bleed |
| Elevated baseline | -- | >10% of peak intensity | Dirty source or contamination |
| Narrow gradient | -- | <70% of median width | Truncated acquisition |

### 10.3 Interpreting Results

- **All green**: Your data looks good. Proceed to search.
- **Yellow warnings**: Worth investigating. See guidance below for each flag type.
- **Red failures**: Strongly consider excluding these files or re-acquiring them. A dead injection or severe carryover will waste search time and can distort downstream normalization.

**Per-flag guidance:**
- **Shape deviation (yellow)**: Look at the Faceted View for this run. If the peak is slightly shifted but the overall shape is similar, it is usually fine -- DIA-NN corrects for RT drift. If the shape is dramatically different (double peak, flat line, or no clear peak), exclude the file. Note: the first 1-2 runs after column conditioning or blank washes often show lower shape correlation -- consider whether warnings on the first injections are expected.
- **RT shift**: A small shift (yellow) is normal for biological variation. A large shift (red) suggests column degradation, gradient issues, or autosampler problems. Check whether the shift is systematic (all late runs shifted) or isolated.
- **Loading anomaly (yellow)**: Slightly higher loading may reflect real biology (e.g., a tissue sample yielded more protein). Dramatically lower loading (red) usually means a failed injection -- exclude it.
- **Late elution / elevated baseline**: Suggests carryover from a previous highly-loaded sample or dirty LC system. The affected run's quantification may be unreliable.

> **Note on small cohorts**: The MAD-based outlier detection requires 3+ runs for cohort-level checks. With very small batches (3-4 runs), the thresholds may be too permissive. Visually inspect the Overlay View as a sanity check.

TIC traces and metrics are saved with your session, so you can review them later without re-extracting.

---

## 10.5 Load from HPC (NEW in v3.7)

The **"Load from HPC"** button in the sidebar provides a quick way to download and analyze results from completed DIA-NN searches on the cluster without navigating the job queue.

1. Click **"Load from HPC"** in the sidebar (visible when SSH is connected)
2. The SSH File Browser opens, filtered to show only `.parquet` files
3. Navigate to your search output directory and select `report.parquet`
4. The file is downloaded via SCP and automatically loaded through the analysis pipeline

This is especially useful for Docker users who run searches via SSH and want to analyze results locally.

---

## 11. Search & Analysis History

### 11.1 Search History

Navigate to **About > Search History** in the navbar. Every DIA-NN search submitted through DE-LIMP is logged with 26 fields including timestamp, backend, search mode, FASTA files, enzyme, mass accuracy, scan window, MBR, normalization, any additional search parameters, and job status.

**Key features:**
- **Expandable rows**: Click any row to reveal full details (enzyme, mass accuracy, scan window, MBR, normalization, extra flags, output directory, job ID)
- **Import Settings**: Click the yellow "Settings" button to apply that search's parameters to the current search UI -- no need to remember what settings you used last time
- **Import Results**: Click the green "Results" button (only for completed searches) to load the search output (`report.parquet`) directly. Works over SSH for remote files. Auto-runs phospho detection and records to Analysis History.
- **View Log**: Click the log icon to display the `search_info.md` metadata file from the search output directory
- **Cross-reference**: The link icon navigates to the matching entry in Analysis History (if that search output was loaded and analyzed)

In short: **Import Settings** copies the search parameters so you can start a new search with the same configuration. **Import Results** loads the actual output data (report.parquet) from that completed search into DE-LIMP for analysis -- no new search needed.

### 11.2 Analysis History

Navigate to **About > Analysis History** in the navbar. Every pipeline run (from any source -- upload, example data, search import, session load) is logged with file details, protein counts, contrast names, and pipeline settings.

**Key features:**
- **Expandable rows**: Click to reveal full file path, FASTA, output directory, and notes
- **Project assignment**: Click "Assign" on any entry to group it into a project. Use existing project names or create new ones.
- **Project filtering**: Select a project from the dropdown to filter the table. Summary cards appear above showing total analyses, date range, and protein count range.
- **Load button**: Re-load a previous analysis directly from the history table

### 11.3 Cross-Referencing

Search History and Analysis History are linked by the output folder -- the directory where DIA-NN saved its results. When you complete a search and load its results, both tables show a link icon on the matching rows. Click the link icon to jump between them -- useful for tracing which search parameters produced which analysis results.

---

## 12. Accessing DE-LIMP

You have multiple options to access DE-LIMP:

### 🌐 Web Browser (No Installation)
* **Hugging Face Spaces:** https://huggingface.co/spaces/brettsp/de-limp-proteomics
* Run directly in your browser without installing R or any packages
* Perfect for quick analyses or trying out the app
* Note: Limited computational resources compared to local installation

### 💻 Local Installation (Recommended for Regular Use)
* Clone or download the DE-LIMP directory from [GitHub](https://github.com/bsphinney/DE-LIMP)
* Run with: `shiny::runApp('/path/to/de-limp/', port=3838, launch.browser=TRUE)` (directory-based)
* The app uses a modular structure (`app.R` + `R/` directory) — launch the project folder, not a single file
* Full computational power of your machine
* Better for large datasets or multiple analyses

### 🐳 Docker + SSH (Windows — Recommended)
* **No R installation required** — DE-LIMP runs entirely inside Docker
* **One-click launcher**: Double-click `Launch_DE-LIMP_Docker.bat` — handles SSH key detection, container startup, and browser launch
* **Shared PC support**: Multiple Windows users on the same PC are handled automatically
* **SSH to HPC**: Auto-connects to your HPC cluster for DIA-NN searches via the SSH file browser
* **Load from HPC**: Download and analyze completed search results with one click
* Full step-by-step guide: [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md)
* Also works on Mac/Linux, though native R installation is typically easier on those platforms

---

## 13. Troubleshooting

| Issue | Solution |
| :--- | :--- |
| **App crashes on startup** | Ensure R 4.5+ is installed. The `limpa` package requires R 4.5 or newer. Download from: https://cloud.r-project.org/ |
| **"limpa package not found"** | Upgrade to R 4.5+, then run: `BiocManager::install('limpa')`. The app will auto-install missing packages on first run. |
| **"Please select a CRAN mirror"** | This should not happen in the current version. If it does, add `options(repos = c(CRAN = "https://cloud.r-project.org"))` at the top of the script. |
| **GSEA fails** | Ensure you are connected to the internet (it needs to download gene ontologies). |
| **GSEA fails with organism error** | The app now auto-detects organism via UniProt API. Ensure internet connection. If detection fails, check that your protein IDs are valid UniProt accessions (e.g., P12345). |
| **Grid View "Object not found"** | Ensure you have run the pipeline first (click "Assign Groups & Run Pipeline"). The Grid View requires processed data. |
| **AI says "No data"** | Click the "▶ Run Pipeline" button first. The AI needs the statistical results to answer questions. Also verify your Gemini API key is entered correctly. |
| **Example data won't download** | Check your internet connection. The file is 46MB and downloads from GitHub Releases. |
| **Can't select multiple proteins in table** | Use Ctrl+Click (Windows/Linux) or Cmd+Click (Mac) for individual selections. Use Shift+Click for range selections. |
| **Session file won't load** | Ensure the `.rds` file was created by DE-LIMP v2.0+. Older versions may not be compatible. |
| **Covariate columns not showing** | Click "Assign Groups & Run Pipeline" to open the modal. Covariate columns (Batch, Covariate1, Covariate2) are in the table. |
| **SSH connection fails** | Check hostname, username, and SSH key path. The key must be passwordless (key-based auth only). Ensure the remote host is reachable and your key is authorized (`~/.ssh/authorized_keys` on the server). |
| **"sbatch not found on remote PATH"** | Click **"🔗 Test Connection"** — the app probes the login shell and common HPC paths (spack, modules, `/usr/local/bin`) automatically to locate SLURM binaries. |
| **Job shows "Unknown" status** | Click the **"🔄 Refresh"** button on the job. This re-queries SLURM via `sacct`. Jobs older than the SLURM accounting retention period may remain unknown. |
| **DIA-NN search fails with "command not found"** | Usually a line continuation issue in the generated sbatch script. Click **"📄 Log"** on the job to view stdout/stderr for details. |
| **Jobs lost after app restart** | Jobs now persist automatically in `~/.delimp_job_queue.rds`. Restart the app to restore your full job history. If the file is missing or corrupted, jobs from before the persistence feature will not be recoverable. |
| **MallocStackLogging warnings on Mac** | Harmless macOS ARM64 warnings from system libraries. These are suppressed in the latest version and do not affect functionality. |
| **Community stats not showing on About tab** | Stats are populated by the `track-stats.yml` GitHub Actions workflow, which runs daily. They will appear after the first successful workflow run. You can also trigger the workflow manually from the Actions tab on GitHub. |
| **Can't find SLURM log files** | As of v3.2, log files are stored in `{output_dir}/logs/`. For older searches, they may still be in the output directory root. The "Log" button checks both locations automatically. |
| **App shows wrong environment badge** | The colored badge (Docker/HPC/Local/HF) is auto-detected. Docker shows red, Apptainer on HPC shows green, native R shows blue. If Docker shows "Local" instead of "Docker", check that Docker environment variables are set. |
| **SSH auto-connect fails on startup** | SSH auto-connect runs when an SSH key is detected. If it hangs, a stale ControlMaster socket may exist. The app probes with `ssh -O check` and removes dead sockets automatically. Check that your SSH key has no passphrase. |
| **NCBI gene symbols not appearing** | For NCBI FASTA databases, gene symbol mapping requires E-utilities access. Docker users without direct internet access to NCBI get the gene map via SSH from HPC. If symbols show as accessions, check the gene map TSV file alongside the FASTA. |
| **File browser only shows limited directories** | The SSH file browser uses configured root directories (`DELIMP_EXTRA_ROOTS` env var) for performance. Ask your admin to add additional paths if your data is elsewhere. |
| **"Load from HPC" button not visible** | This button only appears when SSH is connected. Click "Test Connection" or wait for SSH auto-connect. |

---

## Glossary

| Term | Definition |
| :--- | :--- |
| **DIA** | Data-Independent Acquisition -- a mass spectrometry method that fragments all ions in wide m/z windows, providing comprehensive peptide coverage |
| **DDA** | Data-Dependent Acquisition -- an older MS method that selects the most abundant ions for fragmentation; DE-LIMP does not support DDA data |
| **DIA-NN** | A software tool (by Vadim Demichev) that processes DIA raw files to identify and quantify peptides and proteins |
| **Parquet** | A columnar file format used by DIA-NN for its output (`report.parquet`). More compact and faster to read than TSV |
| **limpa** | A Bioconductor R package for DIA proteomics data processing. Handles data import from DIA-NN, normalization (DPC-CN), and protein quantification (DPC-Quant, a modified maxLFQ algorithm). limpa reads the raw DIA-NN output; limma then performs the statistical testing |
| **limma** | A Bioconductor R package for linear modeling and differential expression analysis. Uses empirical Bayes moderation (see below) to produce reliable statistics even with small sample sizes. Originally designed for microarrays, now widely used in proteomics |
| **DPC-CN** | Data Point Correspondence - Cyclic Normalization -- a normalization method designed for DIA proteomics that adjusts for systematic intensity differences between runs (e.g., differences in sample loading or instrument performance) so that fold-change comparisons reflect biology, not technical variation. Applied by limpa on top of DIA-NN's built-in RT-dependent normalization |
| **Empirical Bayes moderation** | A statistical technique used by limma that borrows information across all ~3,000+ proteins to produce more stable variance estimates for each individual protein. This is especially helpful with few replicates (n=3-4): rather than relying on noisy per-protein variance from just your replicates, limma combines each protein's variance with a prior estimated from the full dataset |
| **FDR** | False Discovery Rate -- the expected proportion of false positives among all significant results. An FDR threshold of 0.05 means the procedure is calibrated so that, on average, no more than 5% of proteins called significant are expected to be false positives |
| **adj.P.Val** | Adjusted P-value (after FDR correction via Benjamini-Hochberg). This is what DE-LIMP uses to determine significance (default threshold: 0.05) |
| **P.Value** | Raw (unadjusted) P-value from the statistical test. Used on the volcano y-axis for visual spread, but significance is determined by adj.P.Val |
| **logFC** | Log2 fold change between conditions. A logFC of 1.0 means 2-fold higher; -1.0 means 2-fold lower; 0.6 means ~1.5-fold higher |
| **Fold change** | The ratio of expression between two conditions. A 2-fold change means one group has twice the abundance of the other |
| **Volcano plot** | A scatter plot showing fold change (x-axis) vs. statistical significance (y-axis) for every protein. Significant proteins with large changes appear in the upper corners |
| **PCA** | Principal Component Analysis -- a dimensionality reduction method that projects samples onto axes of maximum variance, helping visualize how samples relate to each other based on their overall protein expression profiles. PCA does not perform clustering; visual proximity suggests similarity |
| **GSEA** | Gene Set Enrichment Analysis -- a rank-based method that tests whether predefined sets of genes/proteins (e.g., pathways) tend to cluster toward the top or bottom of a fold-change-ranked list, rather than being randomly distributed. This is different from overrepresentation analysis (ORA), which tests a fixed list of significant hits |
| **GO** | Gene Ontology -- a standardized vocabulary for gene/protein function, organized into three categories: Biological Process (BP), Molecular Function (MF), and Cellular Component (CC) |
| **KEGG** | Kyoto Encyclopedia of Genes and Genomes -- a database of metabolic and signaling pathways |
| **CV** | Coefficient of Variation -- a measure of replicate reproducibility (standard deviation / mean, expressed as %). Computed on linear-scale intensities (back-transformed from log2) within each experimental group, then averaged. Lower CV = more reproducible measurement |
| **MOFA2** | Multi-Omics Factor Analysis -- an unsupervised method that finds latent factors (mathematical patterns) shared across multiple data types (e.g., proteomics + phospho). Factors require biological interpretation -- they may capture genuine biology, batch effects, or technical variation |
| **FASTA** | A text file format containing protein or nucleotide sequences, used as the reference database for peptide identification |
| **Precursor** | A peptide ion at a specific charge state as measured by the mass spectrometer. The same peptide can appear at different charge states (e.g., +2 and +3), creating multiple precursor entries. This is why a dataset may show 40,000 precursors but only 4,000 protein groups |
| **maxLFQ** | A label-free quantification algorithm (Cox et al., 2014) that estimates protein-level abundance from pairwise precursor/peptide ratios, robust to missing values. limpa uses a modified version called DPC-Quant |
| **MBR** | Match Between Runs -- a DIA-NN feature that transfers peptide identifications from runs where a peptide was confidently identified to runs where it was not, increasing data completeness. Can add 10-30% more identifications but introduces more imputed values |
| **TIC** | Total Ion Current -- the summed intensity of all ions detected at each time point during a mass spectrometry run. TIC traces show the chromatographic profile of each injection |
| **MAD** | Median Absolute Deviation -- a robust measure of spread, similar to standard deviation but less sensitive to outliers. "3 MAD from median" means the value is far from the group center |
| **Spectral library** | A collection of previously observed peptide fragmentation patterns used to identify peptides in new DIA data |
| **Covariate** | A known source of variation that is not your variable of interest (e.g., batch, sex, instrument). Including covariates in the statistical model separates their effect from your treatment effect, reducing noise and false positives |

---

## Version History

### What's New in v3.5.0 (March 2026)

**Run Comparator** -- Compare two analyses of the same dataset across tools (DE-LIMP vs DE-LIMP, Spectronaut, or FragPipe). 4 diagnostic layers, 7-rule hypothesis engine, optional MOFA2 decomposition, tool-aware AI integration.

**Search & Analysis History** -- Full audit trail for DIA-NN searches (26 parameters) and pipeline runs. Import Settings or Import Results from past searches. Project-based organization.

**Chromatography QC** -- TIC extraction from timsTOF .d files before search. Three views (Faceted, Overlay, Metrics) with automated per-run diagnostics.

**Smart HPC Job Submission** -- Per-user SLURM CPU limit detection, auto-partition switching, FASTA database library with auto-upload, path validation.

**DIA-NN Log Parser** -- Extended with pg-level, proteoforms, library precursor count, pipeline step detection.

### Previous Releases

**v3.2** -- About tab with community dashboard, search log reorganization, DIA-NN log parser, Claude export enhancements, sacct array fix.

**v3.1.1** -- Volcano plot fixes (P.Value/adj.P.Val handling, DE count annotation), CV Analysis scatter plot redesign, Export Data panel, AI Summary HTML export.

**v3.1** -- UI overhaul (page_navbar, dark navbar, hover dropdowns, accordion sidebar, DE Dashboard sub-tabs). Core Facility Mode (SQLite, staff YAML, QC dashboard, Quarto reports).

**v3.0** -- Multi-Omics MOFA2, DIA-NN Docker backend, phosphoproteomics (KSEA, motifs), GSEA expansion (BP/MF/CC/KEGG), all-contrast AI summary.

See [CHANGELOG.md](CHANGELOG.md) for complete release notes.

### Previous: v3.0--v3.1 Highlights

- **UI Overhaul** (v3.1): Professional dark navbar, hover dropdowns, accordion sidebar, DE Dashboard with four sub-tabs
- **Core Facility Mode** (v3.1): SQLite job tracking, staff auto-configuration, instrument QC dashboard, Quarto report generation
- **Multi-Omics MOFA2** (v3.0): 2-6 view integration with variance heatmap, factor weights, sample scores
- **DIA-NN Docker Search** (v3.0): Three backends (Local/Docker/HPC), Windows Docker Compose deployment
- **Phosphoproteomics** (v2.4--v2.5): Site-level DE, KSEA kinase activity, motif analysis
- **GSEA Expansion** (v2.5): 4 databases (BP/MF/CC/KEGG) with per-ontology caching and automatic organism detection
- **XIC Chromatogram Viewer** (v2.1): Fragment-level inspection with MS2 intensity alignment and ion mobility support

See [CHANGELOG.md](CHANGELOG.md) for complete version history.

---

*Happy analyzing!* 🧬
