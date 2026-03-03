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

## ✨ What's New in v3.2.0

- 📊 **About Tab**: Community stats (stars, forks, visitors, clones), trend sparklines, and a live GitHub Discussions feed
- 🔍 **DIA-NN Search Improvements**: Spectral library caching, custom FASTA sequences, cluster resource monitoring
- 📂 **Organized Search Logs**: SLURM logs now stored in a dedicated `logs/` subdirectory

**Previous highlights** (v3.1.1):
- 🌋 Volcano plot fixes: correct significance handling with DE protein counts (up/down annotations)
- 💾 Export Data panel with one-click CSV downloads
- 🤖 AI Summary HTML export for standalone styled reports
- 📈 CV Analysis redesign: interactive scatter plot (logFC vs Avg CV)

**Earlier** (v3.1):
- 🎨 UI overhaul: dark navbar, hover dropdowns, accordion sidebar
- 📊 DE Dashboard sub-tabs: Volcano (+heatmap), Results Table, PCA, CV Analysis

**Earlier** (v3.0): Multi-Omics MOFA2 | Phosphoproteomics (site-level DE, KSEA, motifs) | GSEA 4-database | AI Summary | XIC Viewer ([local/HPC only](https://github.com/bsphinney/DE-LIMP))

## 🚀 Features

### 📊 Interactive Analysis
- **Volcano Plots** - Fully interactive (Plotly). Click or box-select to highlight
- **Heatmaps** - Auto-scaled Z-score heatmaps of significant proteins
- **QC Trends** - Monitor run quality with group averages
- **Multi-Protein Violin Plots** - Compare expression distributions

### 🤖 AI-Powered Exploration
- **Chat with Your Data** - Google Gemini integration
- **Bi-Directional Sync** - Select proteins ↔ AI suggestions
- **Auto-Summary** - Generate draft analysis summaries

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
