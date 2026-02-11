<img src="https://github.com/user-attachments/assets/2aeb5863-2c10-4faa-99e8-25835a9e9330" align="left" width="150" style="margin-right: 20px;" alt="DE-LIMP Logo" />

# DE-LIMP: Differential Expression & Limpa Proteomics

**DE-LIMP** is a modern, interactive R Shiny application for differential expression analysis of DIA-NN proteomics data. Built on **Limpa** and **Limma** for robust statistics, with **Google Gemini AI** integration for intelligent data exploration.

> **Why "DE-LIMP"?**
> Because analyzing differential expression shouldn't make you walk with a limp. 🧬🚶‍♂️

**🌐 Project Website:** [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/)
**📚 Full Documentation:** [USER_GUIDE.md](USER_GUIDE.md) | [CLAUDE.md](CLAUDE.md)

<br clear="left"/>

---

## 🚀 Key Features

### 📊 Interactive Visualizations
- **Volcano Plots** - Fully interactive (Plotly). Click or box-select to highlight across all views
- **Heatmaps** - Auto-scaled Z-score heatmaps of significant proteins
- **QC Trends** - Monitor precursors and proteins across run order with group averages
- **Multi-Protein Violin Plots** - Compare expression distributions across groups
- **MDS & DPC Plots** - Assess sample clustering and normalization quality

### 🤖 AI-Powered Analysis (Gemini Integration)
- **Chat with Your Data** - Upload QC stats + top proteins to Google Gemini
- **Bi-Directional Sync:**
  - Select proteins in plots → AI knows your focus
  - Ask AI questions → Auto-filter plots and tables
- **Auto-Summary** - Generate publication-ready methodology summaries

### 💾 Session Management
- **Save/Load Sessions** - Preserve entire analysis state as .rds files
- **Reproducibility Logging** - Export complete R code for analysis reproduction
- **Example Data** - One-click download of demo dataset (Affinisep vs Evosep comparison)

### 🔬 Advanced Features
- **Multiple Covariates** - Customize covariate names (Batch, Sex, Diet, etc.) and include in models
- **Auto-Guess Groups** - Smart detection of experimental groups from filenames
- **Grid View** - Heatmap-style table with UniProt linking and click-to-plot
- **GSEA Integration** - Gene Ontology enrichment analysis
- **Consistent DE** - Identify highly reproducible significant proteins

### 🎓 Education Tab
- Embedded proteomics resources and training materials
- Latest UC Davis Proteomics YouTube videos
- Links to short courses and methodology citations
- Google NotebookLM for exploring key papers

---

## 🛠 Installation

### System Requirements
- **R 4.5 or newer** (required for limpa package)
- **Bioconductor 3.22+** (automatically configured with R 4.5+)
- Internet connection for package installation

### Quick Start

1. **Download the app:**
   ```bash
   wget https://github.com/bsphinney/DE-LIMP/releases/latest/download/DE-LIMP.R
   # Or download from: https://github.com/bsphinney/DE-LIMP/releases
   ```

2. **Run the app:**
   ```r
   shiny::runApp('DE-LIMP.R', port=3838, launch.browser=TRUE)
   ```

3. **First-time setup:**
   - The app automatically installs missing packages
   - R version check performed on startup
   - Clear upgrade instructions if R 4.5+ is needed

### Dependencies
All dependencies install automatically on first run:
```r
# Core packages
shiny, bslib, plotly, DT, rhandsontable, shinyjs

# Data processing
dplyr, tidyr, stringr, readr, arrow

# Statistics & Bioinformatics
limpa, limma, ComplexHeatmap, clusterProfiler
org.Hs.eg.db, org.Mm.eg.db, AnnotationDbi

# Visualization
ggplot2, ggrepel, ggridges, enrichplot

# AI Integration
httr2, curl
```

---

## 📖 Usage

### 1. Load Data
- **Upload DIA-NN .parquet file** OR
- **Click "Load Example Data"** for instant demo dataset

### 2. Assign Groups & Run Pipeline
- Click "Assign Groups & Run Pipeline"
- Use "Auto-Guess Groups" or manually assign
- Configure covariates (optional): Batch, custom covariate names
- Click "▶ Run Pipeline" to execute limpa analysis

### 3. Explore Results
- **Data Overview** - Signal distributions and QC summaries
- **QC Trends** - Temporal trends with group averages
- **DE Dashboard** - Volcano plots and interactive tables
- **Consistent DE** - High-reproducibility proteins
- **GSEA** - Gene Ontology enrichment
- **Data Chat** - AI-powered exploration

### 4. Export & Reproduce
- Download reproducibility log (.R file)
- Save session (.rds file) for later
- Export tables and plots

---

## 🔬 Methodology

**Normalization:** Data Point Correspondence (DPC-CN)
**Quantification:** maxLFQ algorithm
**Statistics:** limma empirical Bayes moderation
**FDR Correction:** Benjamini-Hochberg

**Key Citations:**
- **limpa:** Bioconductor package for DIA proteomics
- **limma:** Ritchie ME et al. (2015) Nucleic Acids Res
- **DIA-NN:** Demichev V et al. (2020) Nat Methods

Full methodology available in the app's "Reproducibility > Methodology" tab.

---

## 🎓 Resources

- **🌐 Project Website:** [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/)
- **📺 Video Tutorials:** [UC Davis Proteomics YouTube](https://www.youtube.com/channel/UCpulhf8gl-HVxACyJUEFPRw)
- **🏫 Training:** [Hands-On Proteomics Short Course](https://proteomics.ucdavis.edu/events/hands-proteomics-short-course)
- **🏢 Core Facility:** [proteomics.ucdavis.edu](https://proteomics.ucdavis.edu)
- **📰 News:** [Proteomics News Blog](https://proteomicsnews.blogspot.com/)

---

## 📄 License

This project is open source. See repository for license details.

---

## 🤝 Contributing

Issues and pull requests welcome! See [CLAUDE.md](CLAUDE.md) for development documentation.

**Developer:** Brett Phinney, UC Davis Proteomics Core Facility
**Contact:** [GitHub Issues](https://github.com/bsphinney/DE-LIMP/issues)

---

## 📊 Example Data

Demo dataset included in releases: **Affinisep vs Evosep** comparison using 50ng Thermo Hela digest.
Available at: [github.com/bsphinney/DE-LIMP/releases](https://github.com/bsphinney/DE-LIMP/releases)

---

**Built with ❤️ for the proteomics community**
