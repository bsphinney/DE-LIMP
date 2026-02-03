# DE-LIMP (Limpa Proteomics Dashboard) v1.0

**DE-LIMP** is a local, interactive R Shiny dashboard for statistical analysis of proteomics data (DIA-NN outputs). This tool leverages the `limpa` package for Data Processing & Cleaning (DPC) and integrates Google's Gemini AI for bidirectional, context-aware biological insights.

## 🚀 Key Features

### 1. Automated Pipeline
* **Input:** Raw DIA-NN report files (`.parquet`).
* **Processing:** Automatic DPC normalization and imputation using the `limpa` Bioconductor package.
* **Statistics:** Differential Expression (DE) analysis via `limma` with empirical Bayes moderation.

### 2. Deep Quality Control (QC)
* **Interactive Trends:** Sort QC metrics (Precursors, Proteins, MS1) by Run Order or Experimental Group.
* **Group Distributions:** Jittered Violin plots to instantly spot batch effects or failing experimental groups.
* **Comprehensive Signal Overview:** Visualize the full dynamic range of protein signals across all groups.
    *   **Dynamic Range Metric:** Quantitative display of the dataset's signal dynamic range in orders of magnitude.
    *   **Interactive DE Coloring:** Toggle coloring of proteins by Differential Expression (DE) status (Up, Down, Not Significant).
* **Grouped QC Summaries:** Instantly view average Precursors, MS1 Signal, and Protein Groups per experimental condition in a dedicated table.

### 4. Robust Biomarker Discovery
* **Consistent DE Panel:** Ranks significant proteins by **%CV (Coefficient of Variation)** to identify the most stable, reproducible markers across replicates.
* **Reproducibility:** Automatically generates the R code required to reproduce your specific analysis in a standalone script.

### 5. Gene Set Enrichment Analysis (GSEA)
* **Powered by `clusterProfiler`:** Performs Gene Ontology (GO) enrichment analysis on ranked gene lists from DE results.
* **Smart Organism Detection:** Automatically identifies the organism (e.g., Human, Mouse, Rat, etc.) based on protein IDs and loads the correct Bioconductor annotation database (`OrgDb`).
* **Interactive Visualizations:** Explore enrichment results with:
    *   **Dot Plot:** Overview of enriched GO terms.
    *   **Enrichment Map:** Network visualization of overlapping enriched terms.
    *   **Ridgeplot:** Shows expression distribution of core enriched genes across terms.

### 6. AI Data Chat (Gemini Powered)
* **Bidirectional:** Select proteins in the Volcano Plot $\leftrightarrow$ Ask about them in the Chat.
* **One-Click Summarization:** New 'Auto-Analyze' button prompts Gemini to summarize the dataset, highlight key biological/technical patterns, and identify interesting proteins.
* **Enhanced Bidirectional Interaction:** AI-selected proteins are now highlighted in both the Volcano Plot and the Signal Distribution plot, creating a seamless feedback loop.

---

## 🛠️ Installation

### Prerequisites
You need **R** and **RStudio** installed.

### Required Packages
Run this command in your R console to install all dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "limpa", "ComplexHeatmap", "shiny", "shinyjs", "plotly", "DT", "tidyr", "tibble", "stringr", "curl", "bslib", "arrow", "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "enrichplot", "ggridges", "ggrepel"))
