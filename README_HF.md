---
title: DE-LIMP Proteomics
emoji: 🧬
colorFrom: blue
colorTo: green
sdk: docker
sdk_version: 4.5.0
app_file: app.R
pinned: true
license: mit
tags:
  - proteomics
  - shiny
  - bioinformatics
  - differential-expression
  - mass-spectrometry
  - limma
  - dia-nn
---

# DE-LIMP: Differential Expression & Limpa Proteomics 🧬

An interactive R Shiny application for differential expression analysis of DIA-NN proteomics data. Built on **Limpa** (a Bioconductor package for DIA proteomics normalization and quantification) and **Limma** (a widely-used statistical framework for differential expression), with **Google Gemini AI** integration.

## ✨ What's New in v3.7.0

- 🧪 **Contaminant Analysis**: Summary cards, per-sample bar chart, keratin flagging, contaminant heatmap. Signal Distribution and Expression Grid also highlight contaminants.
- 🔬 **Data Explorer**: Quartile-based abundance profiles and sample-sample scatter plots -- explore data without requiring DE analysis
- 🌐 **NCBI Proteome Download**: Search NCBI Datasets by organism, download RefSeq FASTA with automatic gene symbol mapping (local/HPC only)
- 📂 **SSH File Browser**: Visual directory browser for remote HPC navigation with clickable breadcrumbs and color-coded entries (local/HPC only)
- 🔄 **No-Replicates Mode**: Quantification completes normally with n=1 per group; DE is skipped gracefully
- 🏷️ **Environment Badge**: Colored badge showing deployment mode (Docker/HPC/Local/HF)

**Previous highlights** (v3.5): Run Comparator (cross-tool DE comparison), Search & Analysis History, Chromatography QC, smart HPC partitions

**Earlier** (v3.1): UI overhaul (dark navbar, accordion sidebar, DE Dashboard sub-tabs). (v3.0): Multi-Omics MOFA2 | Phosphoproteomics (site-level DE, KSEA, motifs) | GSEA 4-database | AI Summary | XIC Viewer ([local/HPC only](https://github.com/bsphinney/DE-LIMP))

## 🚀 Features

### 📊 Interactive Analysis
- **Volcano Plots** - Fully interactive (Plotly). Click or box-select to highlight
- **Heatmaps** - Auto-scaled Z-score heatmaps of significant proteins
- **Contaminant Analysis** - Summary cards, bar charts, keratin flagging, and contaminant heatmap
- **Data Explorer** - Quartile abundance profiles and sample-sample scatter plots
- **QC Trends** - Monitor run quality with group averages
- **Multi-Protein Violin Plots** - Compare expression distributions

### 🤖 AI-Powered Analysis (Google Gemini)
> **API Key Required:** You must provide your own free Gemini API key. Get one at [Google AI Studio](https://aistudio.google.com/) and paste it into the sidebar. AI Summary sends only summary statistics (protein names, logFC, adj.P.Val). Data Chat sends per-sample expression data for top DE proteins to enable interactive Q&A.

- **AI Summary** - Analyzes all contrasts at once: top DE proteins, cross-comparison biomarkers, and CV stability metrics. Export as standalone HTML report
- **Interactive Data Chat** - Conversational AI with full dataset context (QC stats, top DE proteins, phospho sites when available). Auto-Analyze button for one-click reports
- **Interactive AI + Plot Connection** - Select proteins in volcano/table to set AI context; AI responses highlight proteins in your plots automatically
- **Export for Claude** - Download your complete analysis as a .zip for deep analysis with Claude, ChatGPT, or other AI assistants (includes DE results, expression matrix, QC metrics, GSEA, methods text, and more)
- **Save Chat History** - Download conversations as plain text

### 💾 Session Management
- **Save/Load Sessions** - Preserve analysis state (.rds files)
- **Reproducibility Logging** - Export complete R code
- **Example Data** - One-click demo dataset (Affinisep vs Evosep)

### 🎓 Education & Resources
- Embedded proteomics training materials
- UC Davis Proteomics video tutorials
- Methodology citations (limpa, limma, DIA-NN)

## 📖 Quick Start

1. **Load Data**: Upload DIA-NN .parquet file or use "Load Example Data"
2. **Assign Groups**: Use auto-guess or manual assignment
3. **Run Pipeline**: Click "▶ Run Pipeline" -- the app normalizes your data, quantifies proteins, and runs statistical tests to identify which proteins differ significantly between your groups
4. **Explore Results**: Interactive plots, tables, GSEA, AI chat

## 🔬 Methodology

- **Normalization**: Data Point Correspondence (DPC-CN)
- **Quantification**: maxLFQ algorithm
- **Statistics**: limma empirical Bayes moderation
- **FDR Correction**: Benjamini-Hochberg

## 📚 Resources

- **GitHub**: [github.com/bsphinney/DE-LIMP](https://github.com/bsphinney/DE-LIMP)
- **Discussions**: [github.com/bsphinney/DE-LIMP/discussions](https://github.com/bsphinney/DE-LIMP/discussions)
- **Website**: [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/)
- **YouTube**: [UC Davis Proteomics](https://www.youtube.com/channel/UCpulhf8gl-HVxACyJUEFPRw)
- **Core Facility**: [proteomics.ucdavis.edu](https://proteomics.ucdavis.edu)

## ⚠️ Hosted Version Limitations

This Hugging Face deployment has some limitations compared to a local installation:
- **No DIA-NN search or XIC viewer** -- these require local file access and are only available on local/HPC installations
- **Upload size limits** -- very large files may fail to upload or cause timeouts
- **Sessions don't persist** -- your analysis is lost when you close the browser tab; use "Save Session" to download an .rds file you can reload later
- **Shared resources** -- large datasets or complex analyses may be slower than on a dedicated machine

For the full experience, [install DE-LIMP locally](https://github.com/bsphinney/DE-LIMP).

## 🛠 System Requirements (Local Installation Only)

The following requirements apply only if you are installing DE-LIMP on your own machine. The hosted version above requires nothing but a web browser.

- R 4.5+
- Bioconductor 3.22+
- All dependencies auto-install on first run

## 👨‍🔬 Developer

**Brett Phinney** - UC Davis Proteomics Core Facility

---

**Built with ❤️ for the proteomics community**
