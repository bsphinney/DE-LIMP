# 📘 DE-LIMP User Guide

Welcome to **DE-LIMP** (Differential Expression & Limpa Proteomics), your interactive dashboard for analyzing DIA-NN proteomics data. This guide covers the complete workflow, from importing data to discovering biological insights with our integrated AI assistant.

---

## ✨ What's New in v3.1 (February 2026)

### UI Overhaul
- **Professional dark navbar**: Fixed top navigation bar with hover-activated dropdown menus and smooth animations
- **Dropdown section labels**: The Analysis menu is organized with labeled sections (Setup / Results / AI) for easier navigation
- **Collapsible accordion sidebar**: Upload Data, Pipeline Settings, and AI Chat sections collapse independently — less scrolling, more focus
- **Active tab indicator**: Teal underline on the currently active tab for clear visual feedback

### DE Dashboard Restructured
- **Four sub-tabs**: The DE Dashboard now contains **Volcano** (with heatmap below), **Results Table**, **PCA**, and **Robust Changes** as organized sub-tabs
- **PCA moved here**: PCA plot relocated from Data Overview to the DE Dashboard, placing it closer to your differential expression results
- **Cleaner layout**: Each sub-tab has full width for its content — no more side-by-side cramming

### Core Facility Mode *(Optional — for proteomics core labs)*
- **Activated by environment variable**: Set `DELIMP_CORE_DIR` to a directory containing `staff.yml` to enable
- **SQLite job tracking**: Every DIA-NN search recorded with lab, instrument, project, LC method, and analyst metadata
- **Staff auto-configuration**: Selecting a staff member auto-fills SSH host, username, key path, and SLURM account/partition
- **Search DB tab**: Full-width job history table with 6 filters (text search, lab, status, staff, instrument, LC method) — load results or generate reports from any past search
- **Instrument QC dashboard**: Protein count, precursor count, and TIC trend plots per instrument with ±2SD control lines and date range filtering
- **Quarto report generation**: One-click standalone HTML reports with metadata header, QC bracket comparison, volcano plots, DE statistics, and top proteins table
- **Template system**: Save and load search presets for reproducible workflows across staff members

> *Core Facility Mode adds a "Facility" dropdown to the navbar with Search DB and Instrument QC tabs. It is not visible on Hugging Face Spaces or standard (non-facility) installations.*

---

## Previous: v3.0 Highlights

- **Multi-Omics MOFA2**: 2-6 view integration with variance heatmap, factor weights, sample scores, Factor-DE correlation, and built-in example datasets
- **DIA-NN Docker Search**: Three backends (Local/Docker/HPC), Windows Docker Compose deployment
- **Phosphoproteomics**: Site-level DE, KSEA kinase activity, motif analysis, protein-level abundance correction
- **GSEA Expansion**: 4 databases (BP/MF/CC/KEGG) with per-ontology caching and automatic organism detection via UniProt API
- **AI Summary**: All-contrast analysis with cross-comparison biomarker detection
- **DIA-NN Search Integration**: SSH remote submission, non-blocking job queue, auto-load results, UniProt FASTA download, 6 contaminant libraries
- **XIC Chromatogram Viewer** (v2.1): Fragment-level inspection with MS2 intensity alignment and ion mobility support (local/HPC only)
- **Contextual Help System** (v2.2): 15 info modal `?` buttons with actionable guidance

See [CHANGELOG.md](CHANGELOG.md) for complete version history.

---

## 1. Getting Started

### Which Installation Method?

| Platform | Recommended | What You Need | Guide |
|----------|-------------|---------------|-------|
| **Just exploring** | Web browser | Nothing | [Hugging Face](https://huggingface.co/spaces/brettsp/de-limp-proteomics) |
| **Windows** | Docker Compose | Docker Desktop + Git | [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md) |
| **Mac / Linux** | Native R | R 4.5+ and RStudio | Continue reading below |
| **HPC cluster** | Apptainer | Singularity/Apptainer | [HPC_DEPLOYMENT.md](HPC_DEPLOYMENT.md) |

> **Windows users:** R package installation on Windows is often problematic. We strongly recommend the Docker approach — it bundles everything (R, all packages, and DIA-NN) in one container. See **[WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md)** for a step-by-step walkthrough.

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
    * *Requirement:* The file must be in **`.parquet`** format.
    * *Download Example:* Available at [GitHub Releases](https://github.com/bsphinney/DE-LIMP/releases/tag/v1.0)
* **Q-Value Cutoff:** Adjust the slider to set your False Discovery Rate (FDR) threshold (Default: 0.01).

### Step 2.2: Assign Groups & Run Pipeline
This is the most critical step for statistical analysis. The workflow is streamlined into one modal dialog.

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
6.  **Run the Analysis:**
    * Click the **"▶ Run Pipeline"** button at the top of the modal
    * **What happens?** The app uses the `limpa` package to perform DPC normalization and the `limma` package to fit linear models for differential expression
    * The modal will automatically close and navigate to the QC Plots tab
    * Wait for the status to change to **"✅ Complete!"**

### Step 2.3: Select Comparison
* Use the **"Comparison"** dropdown to select which contrast you want to view (e.g., `Evosep - Affinisep`, `Treatment - Control`).

---

## 3. DIA-NN Database Search

The **"New Search"** tab lets you submit DIA-NN database searches directly from DE-LIMP. Three backends are supported:

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
- **SSH mode:** Type or paste the remote path and click **"🔍 Scan Files"**.
- The scan detects mass spectrometry files and displays them with sizes:
  - `.d` directories (Bruker timsTOF)
  - `.raw` files (Thermo)
  - `.mzML` files (open format)
  - `.wiff` files (SCIEX)

#### FASTA Database
Three sources are available:

**📥 Download from UniProt:**
1. Type an organism name (e.g., "Homo sapiens", "Mus musculus") in the search box
2. Select a proteome from the results dropdown
3. Choose the content type:
   - **One protein per gene** (recommended) — canonical isoform only, smallest and cleanest database
   - **Canonical** — all reviewed canonical sequences
   - **Canonical + isoform** — includes splice variants
4. Click **"Download"** — the FASTA is downloaded to the HPC working directory (uploaded via SCP in SSH mode)

**📂 Pre-staged on server:**
- A dropdown of FASTA files already available on the cluster (pre-downloaded to a shared location)
- Fastest option for commonly used organisms

**📄 Browse / enter path:**
- **Local mode:** Use the file browser to locate any `.fasta` or `.fa` file
- **SSH mode:** Type or paste the full remote path to the FASTA file

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
| **RT profiling** | ✅ On | Retention time-based profiling for improved quantification |
| **Normalization** | On | DIA-NN internal normalization |

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

### 📋 3.5 Job Queue

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
| **📄 Log** | Always | View the full stdout/stderr output from the SLURM job |
| **❌ Cancel** | Queued or Running | Send `scancel` to terminate the job on the cluster |
| **📂 Load** | Completed | Download results (SCP for SSH mode) and load into DE-LIMP pipeline |
| **🔄 Refresh** | Unknown status | Re-query SLURM via `sacct` to update the job status |

- **"🔄 Refresh All"** button appears when any jobs have unknown status, refreshing all job statuses at once.
- **Auto-load:** When enabled, completed jobs automatically download results and load them into the pipeline — no manual "Load" click needed. For SSH mode, results are transferred via SCP.
- **Persistence:** The job queue is saved to `~/.delimp_job_queue.rds` and persists across app restarts. Restarting DE-LIMP restores your full job history.

### 📝 3.6 Search Settings in Methodology

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

This section is **publication-ready** — it uses human-readable names for all parameters and follows standard methods section conventions. The same information is also logged in the **reproducibility code log** for programmatic access.

---

## 4. Deep Dive: The Data Overview & Grid View

### 📊 Data Overview
This is your landing page with 5 sub-tabs:
* **Assign Groups & Run** — Configure experimental groups and run the analysis pipeline
* **Signal Distribution** — Visualizes the dynamic range; automatically colors by DE status with synchronized comparison selector
* **Dataset Summary** — QC statistics and DE protein counts per comparison with directional arrows
* **Replicate Consistency** — Average precursor and protein counts per group
* **Expression Grid** — Heatmap-style table with UniProt linking and click-to-plot
* **AI Summary** — Generate AI-powered analysis summaries (requires Gemini API key)

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
* **Volcano Plot:** Interactive! Click points to select them. Box-select multiple points to analyze a cluster.
    * **Y-axis:** Shows -log10(raw P-Value) following proteomics best practices
    * **Coloring:** Red points indicate FDR-corrected significance (adj.P.Val < 0.05)
    * **Selection:** Single-click for one protein, box-select for multiple proteins
    * *Sync:* Selecting points here updates the Results Table and the AI context
* **Heatmap:** Displayed directly below the volcano plot. Automatically scales and clusters the top 50 significant proteins (or your specific selection).

#### Results Table Sub-tab
* **DE Results Table:** Shows both raw P-values and FDR-adjusted P-values for transparency
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

#### Robust Changes Sub-tab
* **High-Consistency Table:** Proteins with the most reproducible differential expression across replicates
* **CV Distribution:** Histogram of coefficient of variation across experimental groups

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
* **Export Code:** Navigate to **"Reproducibility > Code Log"** tab
* **Download Button:** Click **"📥 Download Reproducibility Log"** to save as a timestamped `.R` file
* **Complete Script:** The exported file includes:
  * All analysis steps in executable R code
  * Session info (R version, package versions)
  * Group assignments and model formulas
  * Parameter settings (Q-value cutoffs, covariates)
* **Publication Ready:** Use the exported code to reproduce your analysis or include in Methods sections

**Methodology Summary:**
* View detailed methodology in the **"Reproducibility > Methodology"** tab
* Includes citations for limpa, limma, and DIA-NN
* Explains normalization (DPC-CN), quantification (maxLFQ), and statistics (empirical Bayes)

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
1.  Select a **database** from the ontology selector dropdown: Biological Process (BP), Molecular Function (MF), Cellular Component (CC), or KEGG pathways.
2.  Click **"Run GSEA"** for the current contrast.
3.  **Automatic organism detection**: The app queries the UniProt API using your protein IDs to determine the correct organism database. For human data this is instant; for non-human data (mouse, rat, etc.) the API lookup runs automatically.
4.  **Per-ontology caching**: Results are cached separately for each database and contrast combination. Switch between BP, MF, CC, and KEGG without re-running the analysis.
5.  **Contrast indicator**: A banner shows which contrast is active. If you change the comparison, a stale-results warning appears prompting you to re-run.
6.  View results as Dot Plots, Enrichment Maps (networks), Ridgeplots, or browse the full Results Table.

### 🔬 Phosphoproteomics

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

### 💾 Save & Load Analysis Sessions
DE-LIMP can save your entire analysis state for later use.

**To Save a Session:**
1. Complete your analysis (upload data, run pipeline, explore results)
2. Navigate to the **"Reproducibility"** tab
3. Click **"Save Current Session"**
4. Choose a filename and location (saves as `.rds` file)
5. The session file includes:
   * Raw data
   * Processed data (normalized, quantified)
   * Statistical results (limma fit object)
   * All group assignments and settings

**To Load a Session:**
1. Click **"Load Session (.rds)"** in the sidebar
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

The **Multi-Omics MOFA2** tab provides unsupervised multi-omics integration using MOFA2 (Multi-Omics Factor Analysis). It discovers latent factors — hidden patterns that explain variation across your datasets.

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

## 8. 🤖 AI Chat (Gemini Integration)

DE-LIMP features a context-aware AI assistant.

### Setup
1.  Expand the **AI Chat** section in the sidebar accordion and paste your **Gemini API Key**.
2.  (Optional) Change the Model Name if you want to use a specific version (Default: `gemini-3-flash-preview`).

### "Chat with Your Data"
You aren't just chatting with a bot; you are chatting with **your specific dataset**.
* **Auto-Analyze:** Click this button to generate a comprehensive report summarizing QC quality and the top biological findings.
* **Ask Questions:**
    * *"Which group has the highest variance?"*
    * *"Are there any mitochondrial proteins upregulated?"*
    * *"Generate a figure caption for the volcano plot."*

### Bi-Directional AI Sync
* **User -> AI:** Select points on the Volcano Plot. Then ask: *"What are the functions of these selected proteins?"*. The AI knows exactly which ones you clicked.
* **AI -> User:** If the AI finds interesting proteins (e.g., *"I found several glycolytic enzymes..."*), it will highlight them in your plots automatically.

---

## 9. Accessing DE-LIMP

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

### 🐳 Docker Compose (Windows — Recommended)
* **No R installation required** — DE-LIMP and DIA-NN run entirely inside Docker
* Build DIA-NN image once, then `docker compose up` — that's it
* Includes embedded DIA-NN search capability out of the box
* Full step-by-step guide: [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md)
* Also works on Mac/Linux, though native R installation is typically easier on those platforms

---

## 10. Troubleshooting

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

---

*Happy analyzing!* 🧬
