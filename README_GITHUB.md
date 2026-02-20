<img src="https://github.com/user-attachments/assets/2aeb5863-2c10-4faa-99e8-25835a9e9330" align="left" width="150" style="margin-right: 20px;" alt="DE-LIMP Logo" />

# DE-LIMP: Differential Expression & Limpa Proteomics

**DE-LIMP** is a modern, interactive R Shiny application for differential expression analysis of DIA-NN proteomics data. Built on **Limpa** and **Limma** for robust statistics, with **Google Gemini AI** integration for intelligent data exploration.

> **Why "DE-LIMP"?**
> Because analyzing differential expression shouldn't make you walk with a limp. 🧬🚶‍♂️

**🌐 Project Website:** [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/)
**📚 Full Documentation:** [USER_GUIDE.md](USER_GUIDE.md) | [CLAUDE.md](CLAUDE.md)

<br clear="left"/>

---

## ✨ What's New in v2.5

**🧬 GSEA Expansion**:
- 4 enrichment databases: GO Biological Process, Molecular Function, Cellular Component, and KEGG Pathways
- Per-ontology caching — switch databases without re-computation
- Automatic organism detection via UniProt REST API (works for any organism)

**🔬 Phosphoproteomics**:
- Auto-detection of phospho-enriched data on upload
- Site-level differential expression via limma
- KSEA kinase activity inference with bar plots and results tables
- Sequence logo motif analysis for regulated phosphosites
- Protein-level abundance correction for stoichiometry changes

**🤖 AI Summary — All Comparisons**:
- Analyzes all contrasts simultaneously with cross-comparison biomarker detection
- Biological insights on high-confidence candidates

**🔬 DIA-NN HPC Search Integration**:
- Submit DIA-NN database searches to an HPC cluster directly from DE-LIMP via SSH
- Non-blocking job queue — submit multiple searches and continue using the app
- Results auto-load when complete, search parameters captured in Methodology tab
- UniProt FASTA download, 6 contaminant libraries, and phosphoproteomics search mode built in

**🧬 Multi-Omics MOFA2**:
- Unsupervised integration of 2-6 data views using MOFA2 (Multi-Omics Factor Analysis)
- Smart data import: RDS (DE-LIMP sessions, limma objects), CSV, TSV, Parquet with auto-log2 detection
- Variance explained heatmap, factor weights browser, sample scores scatter, top features table
- Factor-DE correlation links MOFA factors to differential expression results
- Built-in example datasets: Mouse Brain (2-view proteomics+phospho) and TCGA Breast Cancer (3-view mRNA+miRNA+protein)

**🗺️ MDS Plot Coloring**: Color by Group, Batch, or covariates
**📦 Complete Dataset Export**: Download all contrasts + expression + metadata
**🏗️ Code Modularization**: Split from 5,139-line monolith into app.R + 12 R/ modules

See [CHANGELOG.md](CHANGELOG.md) for full release history.

### Previous: v2.1–2.2 Highlights
- ❓ **Contextual Help System** — 15 info modal buttons across every major tab
- 🌋 **Volcano → Table filtering** and improved DE Dashboard layout
- 📈 **XIC Chromatogram Viewer** with MS2 Intensity Alignment and ion mobility support (Local/HPC)
- 🎯 **Four-Way Comparison Selector Sync** across all tabs
- 🔬 **P-value Distribution Diagnostic** with automated pattern detection

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
- **GSEA Integration** - GO (BP/MF/CC) and KEGG pathways with automatic organism detection
- **Consistent DE** - Identify highly reproducible significant proteins

### 🧪 Phosphoproteomics
- **Auto-Detection** - Phospho-enriched data identified automatically on upload
- **Site-Level DE** - Differential expression at the phosphosite level via limma
- **KSEA Kinase Activity** - Infer upstream kinase activity from phosphosite fold-changes with bar plots and results tables
- **Motif Analysis** - Sequence logo visualization for regulated phosphosites
- **Abundance Correction** - Protein-level correction to isolate stoichiometry changes

### 🧬 Multi-Omics Integration (MOFA2)
- **2-6 data views** — Combine proteomics, phosphoproteomics, transcriptomics, metabolomics, etc.
- **Smart RDS parser** — Import DE-LIMP sessions, limma EList/MArrayLM objects, or plain matrices
- **Sample matching** — Automatic intersection of samples across views with overlap statistics
- **5 results tabs** — Variance Explained heatmap, Factor Weights browser, Sample Scores scatter, Top Features table, Factor-DE Correlation
- **Example datasets** — Mouse Brain (proteomics + phospho) and TCGA Breast Cancer (mRNA + miRNA + protein)
- **Session save/load** — MOFA results persist across sessions; methodology auto-generated

### 🔬 DIA-NN Search Integration
- **Three backends** - Local (embedded), Docker, and HPC (SSH/SLURM) — choose what fits your setup
- **Windows Docker Deployment** - `docker compose up` runs DE-LIMP + DIA-NN with zero R installation. See [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md)
- **HPC Search Submission** - Submit DIA-NN database searches to a remote HPC cluster via SSH without leaving DE-LIMP
- **Non-Blocking Job Queue** - Submit multiple searches and continue using the app; results auto-load when complete
- **UniProt FASTA Download** - Search and download proteome databases directly from UniProt (one-per-gene, reviewed, full, isoforms)
- **Contaminant Libraries** - 6 curated contaminant FASTA libraries bundled (Universal, Cell Culture, Mouse/Rat Tissue, Neuron Culture, Stem Cell Culture)
- **Phosphoproteomics Mode** - Auto-configures DIA-NN for phospho analysis (STY modification, max 3 var mods, `--phospho-output`)
- **Methodology Capture** - Search parameters automatically added to the Methodology tab for publication-ready methods
- **Job Queue Persistence** - Queue survives app restarts; active jobs resume polling automatically

> **DIA-NN License:** DIA-NN is developed by [Vadim Demichev](https://github.com/vdemichev/DiaNN) and is free for academic/non-commercial use. It is **not open source and cannot be redistributed**. DE-LIMP does not bundle DIA-NN — the included build scripts download it directly from the official GitHub release and build a local Docker image on your machine. See the [DIA-NN license](https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md) for full terms.

### 🎓 Education Tab
- Embedded proteomics resources and training materials
- Latest UC Davis Proteomics YouTube videos
- Links to short courses and methodology citations
- Google NotebookLM for exploring key papers

---

## 🚀 Which Installation Should I Use?

| Platform | Recommended Method | DIA-NN Search? | Guide |
|----------|-------------------|----------------|-------|
| **Any (just exploring)** | Web browser | No | [Hugging Face](https://huggingface.co/spaces/brettsp/de-limp-proteomics) |
| **Windows** | Docker Compose | Yes (embedded) | [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md) |
| **Mac** | R/RStudio (native) | Via HPC or Docker | See [Installation](#-installation) below |
| **Linux** | R/RStudio (native) | Via HPC or Docker | See [Installation](#-installation) below |
| **HPC cluster** | Apptainer/Singularity | Via SLURM | [HPC_DEPLOYMENT.md](HPC_DEPLOYMENT.md) |

### 🌐 Web Access (No Installation)
**Try it now:** [huggingface.co/spaces/brettsp/de-limp-proteomics](https://huggingface.co/spaces/brettsp/de-limp-proteomics)
- Run directly in your browser — no R, no Docker, nothing to install
- Perfect for exploring features or quick analyses with small datasets
- Note: No DIA-NN search capability; limited computational resources

### 🐳 Docker Compose (Windows)
- **No R installation required** — DE-LIMP and DIA-NN both run inside Docker
- Build DIA-NN image once, then `docker compose up` — that's it
- Full step-by-step guide: **[WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md)**

> R package installation on Windows is notoriously difficult. Docker sidesteps this entirely.

### 💻 Native R Installation (Mac / Linux)
- Install R 4.5+ and run `shiny::runApp()` — all packages auto-install on first launch
- Full computational power of your machine, best performance
- For DIA-NN searches: connect to an HPC cluster via SSH, or install Docker for local searches
- See [Installation](#-installation) section below

### 🖥️ HPC Deployment
- Deploy on cluster environments (SLURM, PBS, etc.) using Apptainer/Singularity
- **SSH Remote Submission** — Run DE-LIMP on your laptop and submit DIA-NN searches to the cluster
- Full guide: [HPC_DEPLOYMENT.md](HPC_DEPLOYMENT.md)

---

## 🛠 Installation

### System Requirements
- **R 4.5 or newer** (required for limpa package)
- **Bioconductor 3.22+** (automatically configured with R 4.5+)
- Internet connection for package installation

### Quick Start

1. **Download the app:**
   ```bash
   git clone https://github.com/bsphinney/DE-LIMP.git
   cd DE-LIMP
   ```

2. **Run the app:**
   ```r
   shiny::runApp('.', port=3838, launch.browser=TRUE)
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
KSEAapp, ggseqlogo
MOFA2, basilisk, callr

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
- **XIC Viewer** - Click "📈 XICs" on any DE protein to inspect fragment chromatograms (local/HPC only)
- **Consistent DE** - High-reproducibility proteins
- **Phosphoproteomics** - Site-level DE with KSEA kinase activity analysis
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
- **MOFA2:** Argelaguet R et al. (2020) Genome Biology

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
