# 📘 DE-LIMP User Guide

Welcome to **DE-LIMP** (Differential Expression & Limpa Proteomics), your interactive dashboard for analyzing DIA-NN proteomics data. This guide covers the complete workflow, from importing data to discovering biological insights with our integrated AI assistant.

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
4.  **Customize Covariate Names (Optional):**
    * Use the text inputs above the table to rename "Covariate1" and "Covariate2"
    * Examples: "Sex", "Diet", "Age", "Time_Point", "Instrument"
    * Check the boxes to include covariates in the statistical model
    * Only covariates with 2+ unique values will be used
5.  **Run the Analysis:**
    * Click the **"▶ Run Pipeline"** button at the top of the modal
    * **What happens?** The app uses the `limpa` package to perform DPC normalization and the `limma` package to fit linear models for differential expression
    * The modal will automatically close and navigate to the QC Plots tab
    * Wait for the status to change to **"✅ Complete!"**

### Step 2.3: Select Comparison
* Use the **"Comparison"** dropdown to select which contrast you want to view (e.g., `Evosep - Affinisep`, `Treatment - Control`).

---

## 3. Deep Dive: The Data Overview & Grid View

### 📊 Data Overview
This is your landing page.
* **Signal Plot:** Visualizes the dynamic range of your experiment.
    * Click **"Color by DE Status"** to see which proteins are Up/Down-regulated in your current comparison.
* **Group Summary:** A quick table showing average precursor and protein counts per group.

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
