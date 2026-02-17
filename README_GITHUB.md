<img src="https://github.com/user-attachments/assets/2aeb5863-2c10-4faa-99e8-25835a9e9330" align="left" width="150" style="margin-right: 20px;" alt="DE-LIMP Logo" />

# DE-LIMP: Differential Expression & Limpa Proteomics

**DE-LIMP** is a modern, interactive R Shiny application for differential expression analysis of DIA-NN proteomics data. Built on **Limpa** and **Limma** for robust statistics, with **Google Gemini AI** integration for intelligent data exploration.

> **Why "DE-LIMP"?**
> Because analyzing differential expression shouldn't make you walk with a limp. 🧬🚶‍♂️

**🌐 Project Website:** [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/)
**📚 Full Documentation:** [USER_GUIDE.md](USER_GUIDE.md) | [CLAUDE.md](CLAUDE.md)

<br clear="left"/>

---

## ✨ What's New in v2.1

**XIC Chromatogram Viewer** (February 2026):
- 📈 **Fragment-level XIC inspection** - View chromatograms for any differentially expressed protein
- 🔄 **Dual format support** - Auto-detects DIA-NN 1.x (wide) and 2.x (long) XIC formats
- 📊 **Split-axis MS1/MS2** - MS1 precursor on top, MS2 fragments below with independent y-axes
- ⚡ **Ion mobility support** - Mobilogram visualization for timsTOF/PASEF data with clear mode indicator
- 🔍 **Auto-detection** - XIC directory path auto-populates based on report filename

**Major UI/UX Enhancements**:
- 🎯 **Four-Way Comparison Selector Sync** - Signal Distribution, Expression Grid, P-value Distribution, and DE Dashboard selectors now sync automatically
- 📊 **Enhanced Dataset Summary** - DE protein counts per comparison with explicit directional arrows (↑ higher in X, ↓ higher in Y)
- 🔬 **P-value Distribution Diagnostic** - Automated pattern detection with color-coded guidance (healthy, inflation, low power, model issues)
- 📈 **CV Distribution Histogram** - Visual distribution of protein variability per experimental group
- 🤖 **AI Summary Relocated** - Moved from modal popup to dedicated Data Overview sub-tab with inline display
- 🌋 **Volcano Plot Annotations** - FDR and logFC thresholds clearly labeled with significance criteria legend box
- 🎨 **Signal Distribution Enhancement** - Always shows DE coloring automatically with dedicated comparison selector
- 🎨 **Responsive Plot Heights** - All plots use viewport-relative units for optimal viewing on any screen size

See [CHANGELOG.md](CHANGELOG.md) for full release history.

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
- **Group Assignment Templates** - Export/import group configurations as CSV for reproducible workflows

### 📈 XIC Chromatogram Viewer (Local/HPC)
- **Fragment-level validation** - Inspect chromatograms for any DE protein
- **MS2 Intensity Alignment** - Spectronaut-style stacked bar chart for fragment ion ratio consistency with automatic inconsistency detection
- **Split-axis MS1/MS2** - Independent y-axes prevent MS1 from squishing fragment peaks
- **Ion mobility support** - Mobilogram visualization for timsTOF/PASEF instruments
- **Smart auto-detection** - XIC directory auto-populates and auto-loads when data is uploaded
- **Dual format** - Supports DIA-NN 1.x and 2.x XIC output formats

> *Not available on Hugging Face Spaces (XIC files too large for cloud). Download DE-LIMP.R to use locally.*

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

## 🚀 Deployment & Access

**DE-LIMP is available in three ways:**

### 🌐 Web Access (No Installation Required)
**Try it now:** [huggingface.co/spaces/brettsp/de-limp-proteomics](https://huggingface.co/spaces/brettsp/de-limp-proteomics)
- Run directly in your browser
- No R installation needed
- Perfect for quick analyses and exploring features
- Note: Limited computational resources compared to local installation

### 💻 Local Installation (Recommended for Regular Use)
- Full computational power of your machine
- Better for large datasets and multiple analyses
- See Installation section below for setup instructions
- Download from: [GitHub Releases](https://github.com/bsphinney/DE-LIMP/releases)

### 🖥️ HPC Deployment (High-Performance Computing)
- Deploy on cluster environments (SLURM, PBS, etc.)
- Use Apptainer/Singularity containers
- Three deployment options including direct pull from Hugging Face
- Full guide: [HPC_DEPLOYMENT.md](HPC_DEPLOYMENT.md)

> **Note for Developers:** This repository is deployed to both GitHub (source code + releases) and Hugging Face Spaces (web deployment). The repositories track the same codebase but have different README.md files due to platform requirements.

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
- **DE Dashboard** - Volcano plots and interactive tables (current comparison shown prominently)
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
