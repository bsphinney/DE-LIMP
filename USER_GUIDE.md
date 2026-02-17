# 📘 DE-LIMP User Guide

Welcome to **DE-LIMP** (Differential Expression & Limpa Proteomics), your interactive dashboard for analyzing DIA-NN proteomics data. This guide covers the complete workflow, from importing data to discovering biological insights with our integrated AI assistant.

---

## ✨ What's New in v2.1 (February 2026)

### XIC Chromatogram Viewer (NEW)
- **Inspect fragment-level chromatograms** for any differentially expressed protein
- Supports **DIA-NN 1.x and 2.x** XIC formats with automatic detection
- **Split-axis MS1/MS2 view**: MS1 precursor on top, MS2 fragments below with independent y-axes
- **Ion Mobility support**: Mobilogram visualization for timsTOF/PASEF data with prominent mode indicator
- Three display modes: Facet by sample, Facet by fragment, and Intensity alignment
- **MS2 Intensity Alignment**: Spectronaut-style stacked bar chart for fragment ion ratio consistency with automatic inconsistency detection
- **Auto-detection**: XIC directory path auto-populates when data is loaded; XICs auto-load when detected
- Available for local and HPC deployments (XIC files too large for cloud)

### Major UI/UX Enhancements

**Comparison Selector Synchronization**
- All comparison selectors now sync automatically across tabs
- Purple gradient banners show active comparison on Signal Distribution, Expression Grid, and P-value Distribution tabs
- Change any selector → all others update immediately

**Enhanced Visualizations**
- **Volcano Plot**: Now includes colored threshold lines (blue FDR, orange logFC) with significance criteria legend box
- **Dataset Summary**: Shows DE protein counts per comparison with explicit directional arrows (e.g., "↑ 189 proteins higher in Evosep")
- **P-value Distribution**: New diagnostic plot with automated pattern detection and color-coded guidance banners
- **CV Distribution**: New histogram showing protein variability distribution per experimental group
- **Signal Distribution**: Always shows DE coloring automatically (no manual toggle needed)

**AI Summary**
- Moved from modal popup to dedicated sub-tab in Data Overview
- Results display inline for easier workflow

**Responsive Design**
- All plots now use viewport-relative heights for optimal viewing on any screen size

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
1.  Open `DE-LIMP.R` in RStudio.
2.  Click the **"Run App"** button at the top right of the script editor.
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

## 3. Deep Dive: The Data Overview & Grid View

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

## 4. Visualizing Results

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
1.  Click **"Run GSEA"**.
2.  The app automatically detects if your data is Human, Mouse, or Rat based on UniProt suffixes.
3.  View results as Dot Plots, Enrichment Maps (networks), or Ridgeplots to understand biological pathways.

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

## 5. 🤖 AI Chat (Gemini Integration)

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

## 6. Accessing DE-LIMP

You have multiple options to access DE-LIMP:

### 🌐 Web Browser (No Installation)
* **Hugging Face Spaces:** https://huggingface.co/spaces/brettsp/de-limp-proteomics
* Run directly in your browser without installing R or any packages
* Perfect for quick analyses or trying out the app
* Note: Limited computational resources compared to local installation

### 💻 Local Installation (Recommended for Regular Use)
* Download `DE-LIMP.R` from [GitHub Releases](https://github.com/bsphinney/DE-LIMP/releases)
* Run with: `shiny::runApp('DE-LIMP.R', port=3838, launch.browser=TRUE)`
* Full computational power of your machine
* Better for large datasets or multiple analyses

### 🐳 Docker Deployment
* For IT administrators or advanced users
* Dockerfile available in the [GitHub repository](https://github.com/bsphinney/DE-LIMP)
* Consistent environment across platforms

---

## 7. Troubleshooting

| Issue | Solution |
| :--- | :--- |
| **App crashes on startup** | Ensure R 4.5+ is installed. The `limpa` package requires R 4.5 or newer. Download from: https://cloud.r-project.org/ |
| **"limpa package not found"** | Upgrade to R 4.5+, then run: `BiocManager::install('limpa')`. The app will auto-install missing packages on first run. |
| **"Please select a CRAN mirror"** | This should not happen in the current version. If it does, add `options(repos = c(CRAN = "https://cloud.r-project.org"))` at the top of the script. |
| **GSEA fails** | Ensure you are connected to the internet (it needs to download gene ontologies). |
| **Grid View "Object not found"** | Ensure you have run the pipeline first (click "Assign Groups & Run Pipeline"). The Grid View requires processed data. |
| **AI says "No data"** | Click the "▶ Run Pipeline" button first. The AI needs the statistical results to answer questions. Also verify your Gemini API key is entered correctly. |
| **Example data won't download** | Check your internet connection. The file is 46MB and downloads from GitHub Releases. |
| **Can't select multiple proteins in table** | Use Ctrl+Click (Windows/Linux) or Cmd+Click (Mac) for individual selections. Use Shift+Click for range selections. |
| **Session file won't load** | Ensure the `.rds` file was created by DE-LIMP v2.0+. Older versions may not be compatible. |
| **Covariate columns not showing** | Click "Assign Groups & Run Pipeline" to open the modal. Covariate columns (Batch, Covariate1, Covariate2) are in the table. |

---

*Happy analyzing!* 🧬
