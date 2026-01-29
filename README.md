# Limpa Proteomics AI Dashboard (v1.0)

A local, interactive R Shiny dashboard for statistical analysis of proteomics data (DIA-NN outputs). This tool leverages the `limpa` package for Data Processing & Cleaning (DPC) and integrates Google's Gemini AI for bidirectional, context-aware biological insights.

## 🚀 Key Features

### 1. Automated Pipeline
* **Input:** Raw DIA-NN report files (`.parquet`).
* **Processing:** Automatic DPC normalization and imputation using the `limpa` Bioconductor package.
* **Statistics:** Differential Expression (DE) analysis via `limma` with empirical Bayes moderation.

### 2. Deep Quality Control (QC)
* **Interactive Trends:** Sort QC metrics (Precursors, Proteins, MS1) by Run Order or Experimental Group.
* **Group Distributions:** Jittered Violin plots to instantly spot batch effects or failing experimental groups.

### 3. AI Data Chat (Gemini Powered)
* **Bidirectional:** Select proteins in the Volcano Plot $\leftrightarrow$ Ask about them in the Chat.
* **Token-Safe Analysis:** Uses a "Smart Context" system that uploads the Top 800 significant proteins + your manual selections directly to Google's File API, bypassing standard context limits.
* **Context Aware:** The AI receives your full QC table and Experimental Design, allowing it to answer technical questions (e.g., "Which group has the lowest MS1 signal?").

### 4. Robust Biomarker Discovery
* **Consistent DE Panel:** Ranks significant proteins by **%CV (Coefficient of Variation)** to identify the most stable, reproducible markers across replicates.
* **Reproducibility:** Automatically generates the R code required to reproduce your specific analysis in a standalone script.

---

## 🛠️ Installation

### Prerequisites
You need **R** and **RStudio** installed.

### Required Packages
Run this command in your R console to install all dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "limpa", "ComplexHeatmap", "shiny", "shinyjs", "plotly", "DT", "tidyr", "tibble", "stringr", "curl", "bslib", "arrow"))