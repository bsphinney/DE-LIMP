<img src="https://github.com/user-attachments/assets/2aeb5863-2c10-4faa-99e8-25835a9e9330" align="left" width="150" style="margin-right: 20px;" alt="DE-LIMP Logo" />

# DE-LIMP: Differential Expression & Limpa Proteomics

**DE-LIMP**  is a modern, interactive R Shiny application designed for the statistical analysis of proteomics data (specifically DIA-NN outputs). It leverages the **Limpa** and **Limma** packages for robust statistics and integrates **Google Gemini AI** for intelligent, conversational data exploration.

> **Why "DE-LIMP"?**
> Because analyzing differential expression shouldn't make you walk with a limp. 🧬🚶‍♂️

<br clear="left"/>

## 🚀 Key Features

### 1. 📊 Interactive Dashboard
* **Volcano Plots:** Fully interactive (Plotly). Select points to highlight them across all other views.
* **Heatmaps:** Auto-scaled Z-score heatmaps of top differentially expressed proteins.
* **QC Trends:** Monitor precursor and protein counts across run order to spot batch effects.

### 2. 🤖 AI Chat with Data (Powered by Gemini)
* **"Chat with your Data":** The app securely uploads your processed dataset (top significant proteins + QC stats) to Google Gemini's File API.
* **Bi-Directional Sync:**
    * **You Select:** Highlight proteins in the Volcano Plot -> AI knows which ones you are interested in.
    * **AI Selects:** Ask the AI "Show me proteins related to glycolysis" -> The app automatically filters the plots and tables to show those proteins.
* **Auto-Summary:** Generate publication-ready methodology and biological summaries with one click.

### 3. 🔬 Deep-Dive Grid View (New!)
* **Bi-Directional Filtering:** The grid instantly filters to match proteins selected in the Volcano Plot or by the AI.
* **Heatmap-Style Columns:** Expression values are color-coded (Blue-White-Red) for quick visual scanning.
* **Smart Headers:** Column headers use compact Run Numbers to save space; hover over them to see the full filename and experimental group.
* **UniProt Linking:** Click any Protein ID to open its official UniProt page in a new tab.
* **Click-to-Plot:** Click any row to immediately see a **Violin Plot** of that protein's expression across all groups.
* **Smart Export:** Dedicated "Export Full Table" button downloads the data with full original filenames for publication.

## 🛠 Installation

### Prerequisites
You need **R** (version 4.2+) and the following packages. The app will attempt to auto-install missing CRAN/Bioconductor packages on first run.

```r
# Key dependencies
install.packages(c("shiny", "bslib", "plotly", "httr2", "tidyverse"))
BiocManager::install(c("limpa", "limma", "ComplexHeatmap", "clusterProfiler"))
