# DE-LIMP: Differential Expression & Limpa Proteomics

**DE-LIMP** is a modern, interactive R Shiny application designed for the statistical analysis of proteomics data (specifically DIA-NN outputs). It leverages the **Limpa** and **Limma** packages for robust statistics and integrates **Google Gemini AI** for intelligent, conversational data exploration.

> **Why "DE-LIMP"?**
> Because analyzing differential expression shouldn't make you walk with a limp. 🧬🚶‍♂️

<img width="271" height="489" alt="image" src="https://github.com/user-attachments/assets/2aeb5863-2c10-4faa-99e8-25835a9e9330" />


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
* **UniProt Linking:** Click any Protein ID to open its official UniProt page.
* **Click-to-Plot:** Click any row to immediately see a **Violin Plot** of that protein's expression across all groups.
* **Smart Export:** Download the full table (with full original filenames) for your supplementary data.

## 🛠 Installation

### Prerequisites

You need **R** (version 4.2+) and the following packages. The app will attempt to auto-install missing CRAN/Bioconductor packages on first run.

```r
# Key dependencies
install.packages(c("shiny", "bslib", "plotly", "httr2", "tidyverse"))
BiocManager::install(c("limpa", "limma", "ComplexHeatmap", "clusterProfiler"))

```

### Running the App

1. Clone this repository.
2. Open `LIMP-D.R` (or `app.R`) in RStudio.
3. Click **Run App**.

## 🧬 Analysis Workflow

1. **Upload:** Load your `report.parquet` file from DIA-NN.
2. **Assign Groups:** Use the "Auto-Guess" feature to parse groups from filenames, or manually edit them in the Excel-like table.
3. **Run Pipeline:** The app performs DPC normalization (via `limpa`) and linear modeling (via `limma`).
4. **Explore:**
* Check **QC Plots** for outliers.
* Use the **DE Dashboard** to find significant candidates.
* Open **Grid View** to inspect raw values and check specific proteins.


5. **Chat:** Enter your Google Gemini API key to start asking biological questions about your results.

## 🤖 Developer Notes (AI & Contributors)

If you are modifying this code or using an AI assistant to help you, please note the following architectural rules to prevent crashes:

1. **State Management:** All app state is stored in a single `reactiveValues` list named `values`.
2. **Critical Syntax Rule:**
* **CORRECT:** `values$log <- new_data` (Treat it like a list).
* **FORBIDDEN:** `values$log(new_data)` (Do NOT treat it like a function).


3. **Data Flow:** The app uses `httr2` to communicate with the Gemini File API. Large datasets are uploaded as CSVs rather than pasted into the prompt to save token space.

## 📄 License

MIT License. Feel free to fork and limp along with us!
