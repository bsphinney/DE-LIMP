# 📘 DE-LIMP User Guide

Welcome to **DE-LIMP** (Differential Expression & Limpa Proteomics), your interactive dashboard for analyzing DIA-NN proteomics data. This guide covers the complete workflow, from importing data to discovering biological insights with our integrated AI assistant.

---

## ✨ What's New in v2.5 (February 2026)

### GSEA Expansion
- **Four gene-set databases**: Biological Process (BP), Molecular Function (MF), Cellular Component (CC), and KEGG pathways
- **Ontology selector dropdown**: Choose which database to run from the GSEA tab
- **Per-ontology caching**: Results are cached separately for each database and contrast — switch instantly without re-running
- **Automatic organism detection**: The app queries the UniProt API to identify organism from your protein IDs, supporting non-human datasets without manual configuration

### AI Summary — All Comparisons
- New **"Summarize All Comparisons"** mode analyzes every contrast simultaneously
- Identifies **cross-comparison biomarkers** (proteins significant in multiple contrasts) and shared biological themes
- **Info modal** explains exactly what data is sent to the AI and how results are generated

### MDS Plot Coloring
- Color samples by **Group**, **Batch**, or any covariate column
- **Colorblind-friendly palette** applied across all MDS color modes

### Complete Dataset Export
- **Download all results at once** from the Dataset Summary tab
- Export includes: all contrast tables, full expression matrix, and sample metadata

### Phosphoproteomics
- **Phase 1 — Site-Level DE**: Auto-detection of phospho data on upload, site-level differential expression with volcano plot, results table, residue distribution, and QC completeness
- **Phase 2 — Kinase Activity & Motifs**: KSEA kinase activity inference from phosphosite fold-changes, sequence logo motif analysis for enriched residue patterns
- **Phase 3 — Advanced**: Protein-level abundance correction to isolate phosphorylation stoichiometry, AI context integration for phospho results

### DIA-NN HPC Search Integration
- **New "New Search" tab** appears when HPC/SSH mode is available — submit DIA-NN database searches directly from DE-LIMP
- **SSH remote submission**: Connect to your HPC cluster via SSH key-based authentication and submit SLURM jobs without leaving the app
- **Non-blocking job queue**: Submit multiple DIA-NN searches and continue using DE-LIMP while jobs run on the cluster
- **Auto-load results**: Completed jobs automatically download results (via SCP for SSH) and load into the analysis pipeline, navigating to the Assign Groups tab
- **File scanning**: Browse or scan remote directories for raw data files (.d, .raw, .mzML, .wiff) with file sizes
- **FASTA database sources**: Download from UniProt (search by organism), select pre-staged server FASTAs, or browse/enter a path
- **Contaminant libraries**: 6 curated options from HaoGroup-ProtContLib (Universal, Cell Culture, Mouse Tissue, Rat Tissue, Neuron Culture, Stem Cell Culture)
- **Search modes**: Library-free (default), Library-based, and Phosphoproteomics (auto-configures STY mods, --phospho-output)
- **Methodology integration**: Search parameters automatically appear in the Methodology tab as a "0. DIA-NN DATABASE SEARCH" section

### Code Modularization
- App split from a single `DE-LIMP.R` monolith into a modular `R/` directory with `app.R` orchestrator
- 12 focused module files for easier maintenance and development

---

## Previous: v2.1-2.2 Highlights

- **XIC Chromatogram Viewer** (v2.1): Fragment-level chromatogram inspection with DIA-NN 1.x/2.x support, split-axis MS1/MS2 view, ion mobility, and MS2 intensity alignment
- **Comparison Selector Sync** (v2.1): All comparison selectors sync automatically across tabs
- **Enhanced Volcano Plot** (v2.1): Colored threshold lines, significance legend, and box-select filtering
- **P-value Distribution Diagnostic** (v2.1): Automated pattern detection with color-coded guidance banners
- **Contextual Help System** (v2.2): 15 info modal `?` buttons across every major tab with actionable guidance
- **Volcano → Table Filtering** (v2.2): Selecting proteins in the volcano plot filters the results table automatically
- **Responsive Design**: All plots use viewport-relative heights for optimal viewing on any screen size

See [CHANGELOG.md](CHANGELOG.md) for complete version history.

---

## 1. Getting Started

### Prerequisites
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

Follow the sidebar controls on the left to process your data.

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

## 3. DIA-NN Database Search (HPC Integration)

The **"New Search"** tab lets you submit DIA-NN database searches to an HPC cluster directly from DE-LIMP. This feature is available when running locally with access to a SLURM-based HPC — either on the cluster itself or remotely via SSH.

> **Note:** The New Search tab only appears when HPC/SSH mode is available. It is not shown on the Hugging Face web version.

### 🔌 3.1 Connection Mode

At the top of the New Search tab, a **connection mode toggle** lets you choose how to reach the cluster:

- **Local (on HPC):** You are running DE-LIMP directly on the cluster head node. SLURM commands (`sbatch`, `sacct`, `scancel`) are on your PATH. File browsers work natively.
- **Remote (SSH):** You are running DE-LIMP on your local machine (e.g., a Mac or laptop) and connecting to the HPC over SSH.

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
This is your landing page with 6 sub-tabs:
* **Assign Groups & Run** — Configure experimental groups and run the analysis pipeline
* **Signal Distribution** — Visualizes the dynamic range; automatically colors by DE status with synchronized comparison selector
* **Dataset Summary** — QC statistics and DE protein counts per comparison with directional arrows
* **Group QC Summary** — Average precursor and protein counts per group
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
* **Current Comparison Display (NEW in v2.0.1):** A prominent blue header banner at the top shows which comparison you're viewing (e.g., "Evosep - Affinisep"). This updates automatically when you change the comparison dropdown, making it easy to keep track of your current analysis focus.
* **Volcano Plot:** Interactive! Click points to select them. Box-select multiple points to analyze a cluster.
    * **Y-axis:** Shows -log10(raw P-Value) following proteomics best practices
    * **Coloring:** Red points indicate FDR-corrected significance (adj.P.Val < 0.05)
    * **Selection:** Single-click for one protein, box-select for multiple proteins
    * *Sync:* Selecting points here updates the Grid View and the AI context
* **Results Table:** Shows both raw P-values and FDR-adjusted P-values for transparency
* **Violin Plots:** Select one or more proteins and click **"📊 Violin Plot"** button
    * Multi-protein support: View multiple proteins in a 2-column grid layout
    * Individual scales: Each protein gets its own Y-axis for better visualization
    * Dynamic height: Adjusts based on number of selected proteins
* **Heatmap:** Automatically scales and clusters the top 50 significant proteins (or your specific selection).

### 📐 QC Trends & Plots
* **Trends:** Monitor precursors and proteins across run order with automatic group averages
    * **Group Average Lines:** Dashed horizontal lines show the mean for each experimental group
    * Color-coded by group for easy comparison
    * **Fullscreen View:** Click **"🔍 View Fullscreen"** to open the plot in a large modal for detailed inspection
    * Helps spot batch effects, instrument drift, or quality issues
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

## 6. 🤖 AI Chat (Gemini Integration)

DE-LIMP features a context-aware AI assistant.

### Setup
1.  Paste your **Gemini API Key** in the sidebar.
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

## 7. Accessing DE-LIMP

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

### 🐳 Docker Deployment
* For IT administrators or advanced users
* Dockerfile available in the [GitHub repository](https://github.com/bsphinney/DE-LIMP)
* Consistent environment across platforms

---

## 8. Troubleshooting

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
