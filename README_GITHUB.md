<img src="https://github.com/user-attachments/assets/2aeb5863-2c10-4faa-99e8-25835a9e9330" align="left" width="150" style="margin-right: 20px;" alt="DE-LIMP Logo" />

# DE-LIMP: Differential Expression & Limpa Proteomics

An R Shiny application for differential expression analysis of **label-free DIA proteomics data** (DIA-NN output). Uses the [limpa](https://bioconductor.org/packages/limpa/) pipeline (DPC-CN normalization, modified maxLFQ quantification) with [limma](https://bioconductor.org/packages/limma/) empirical Bayes moderated *t*-statistics and Benjamini-Hochberg FDR correction.

**Input:** DIA-NN `report.parquet` (or `report.tsv`) | **Not for:** DDA data, TMT/iTRAQ, Spectronaut/MaxQuant output

<br clear="left"/>

**Try it now:** [huggingface.co/spaces/brettsp/de-limp-proteomics](https://huggingface.co/spaces/brettsp/de-limp-proteomics) -- no installation required

**Project Website:** [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/) | **Docs:** [USER_GUIDE.md](USER_GUIDE.md) | [CLAUDE.md](CLAUDE.md)

---

## What's New in v3.1

**UI Overhaul** -- Dark navbar with hover-activated dropdowns, collapsible accordion sidebar, DE Dashboard restructured into sub-tabs (Volcano, Results Table, PCA, Robust Changes), PCA moved from Data Overview into DE Dashboard.

**Core Facility Mode** *(optional)* -- Activated by `DELIMP_CORE_DIR` env var. SQLite job tracking with lab/instrument/project metadata. Staff YAML auto-configures SSH/SLURM. Search DB tab with 6-filter job history. Instrument QC dashboard with trend plots and control lines. One-click Quarto HTML report generation. Template system for reproducible presets.

**v3.0 highlights:** Multi-Omics MOFA2 (2-6 views), DIA-NN Docker search backend, phosphoproteomics (site-level DE, KSEA, motif analysis), GSEA expansion (BP/MF/CC/KEGG), all-contrast AI summary.

See [CHANGELOG.md](CHANGELOG.md) for full release history.

---

## Key Features

### Analysis & Visualization
- **Volcano Plots** -- Interactive (Plotly), click or box-select proteins to highlight across all views; all pairwise contrasts available
- **Heatmaps** -- Z-score heatmaps of selected or significant proteins (ComplexHeatmap)
- **QC Sample Metrics** -- Faceted trend plot (Precursors, Proteins, MS1 Signal, Data Completeness) with LOESS smoother for drift detection and group average lines
- **MDS & DPC Plots** -- Sample clustering and normalization diagnostics
- **Covariates** -- Include batch, sex, diet, or custom covariates in the linear model
- **XIC Chromatogram Viewer** -- Fragment-level chromatogram validation, MS2 intensity alignment (Spectronaut-style), ion mobility/mobilogram support for timsTOF, DIA-NN v1/v2 formats (local/HPC only)
- **Robust Changes** -- Identify highly reproducible DE proteins via coefficient of variation analysis across replicates

### Phosphoproteomics
- **Auto-detection** of phospho-enriched data on upload (scans for UniMod:21 in Modified.Sequence)
- **Phosphosite-level DE** via limma (independent from protein-level analysis); supports DIA-NN `site_matrix_*.parquet` or parsed from `report.parquet`
- **KSEA** (Kinase-Substrate Enrichment Analysis) -- infer upstream kinase activity from phosphosite fold-changes using PhosphoSitePlus + NetworKIN databases
- **Motif analysis** -- sequence logos (ggseqlogo) of flanking residues around regulated phosphosites
- **Abundance correction** -- subtract protein-level logFC from site logFC to isolate phosphorylation stoichiometry changes

### Gene Set Enrichment & Multi-Omics
- **GSEA** -- GO (BP/MF/CC) and KEGG pathways via clusterProfiler; per-ontology caching; automatic organism detection (12 species via UniProt REST API or protein ID suffix)
- **MOFA2** (Multi-Omics Factor Analysis) -- unsupervised integration of 2-6 data views (e.g., proteomics + phosphoproteomics + transcriptomics). Import from RDS, CSV, TSV, or Parquet. Variance explained heatmap, factor weights, sample scores, Factor-DE correlation. Built-in example datasets (Mouse Brain, TCGA Breast Cancer)

### AI-Powered Analysis (Google Gemini)
- Query DE results via Google Gemini API -- sends QC stats and top DE proteins as context (requires your own API key)
- Bidirectional protein selection: select in plots to set AI context, AI responses can highlight proteins in the app
- Auto-generated methodology text for methods sections

### DIA-NN Search Integration
- **Three backends** -- Local, Docker, and HPC (SSH/SLURM)
- **Windows Docker** -- `docker compose up` runs DE-LIMP + DIA-NN with zero R installation ([guide](WINDOWS_DOCKER_INSTALL.md))
- **UniProt FASTA download** -- Search and download proteome databases directly; 6 bundled contaminant libraries
- **Non-blocking job queue** -- Submit multiple searches, results auto-load on completion
- **Phospho mode** -- Auto-configures DIA-NN for phospho analysis (STY modification, `--phospho-output`)

> **DIA-NN License:** DIA-NN is developed by [Vadim Demichev](https://github.com/vdemichev/DiaNN) and is free for academic/non-commercial use. It is not open source and cannot be redistributed. DE-LIMP does not bundle DIA-NN. See the [DIA-NN license](https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md).

### Core Facility Mode *(Optional)*
- Staff YAML profiles auto-fill SSH, SLURM, and instrument settings
- SQLite job tracking with searchable history (6 filters), one-click result loading and report generation
- Instrument QC dashboard with protein/precursor/TIC trends and control lines
- Quarto HTML reports with QC bracket, volcanos, DE stats, and top proteins

> *Activated by setting `DELIMP_CORE_DIR`. Not visible on standard installations.*

### Session Management & Education
- Save/load full analysis state as `.rds`; export reproducibility R code log
- One-click example data (Affinisep vs Evosep comparison)
- Group assignment templates (CSV export/import)
- Embedded proteomics resources, UC Davis Proteomics videos, short course links

---

## Which Installation Should I Use?

| Platform | Method | DIA-NN Search? | Guide |
|----------|--------|----------------|-------|
| **Any (just exploring)** | Web browser | No | [Hugging Face](https://huggingface.co/spaces/brettsp/de-limp-proteomics) |
| **Windows** | Docker Compose | Yes (embedded) | [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md) |
| **Mac / Linux** | R/RStudio (native) | Via HPC or Docker | See [Installation](#installation) below |
| **HPC cluster** | Apptainer/Singularity | Via SLURM | [HPC_DEPLOYMENT.md](HPC_DEPLOYMENT.md) |

---

## Installation

**Requirements:** R 4.5+ (for limpa), Bioconductor 3.22+ (auto-configured with R 4.5+)

```bash
git clone https://github.com/bsphinney/DE-LIMP.git
cd DE-LIMP
```

```r
shiny::runApp('.', port=3838, launch.browser=TRUE)
```

All dependencies install automatically on first run:
```r
# Core: shiny, bslib, plotly, DT, rhandsontable, shinyjs
# Data: dplyr, tidyr, stringr, readr, arrow
# Stats: limpa, limma, ComplexHeatmap, clusterProfiler
#        org.Hs.eg.db, org.Mm.eg.db, AnnotationDbi
#        KSEAapp, ggseqlogo, MOFA2, basilisk, callr
# Viz:  ggplot2, ggrepel, ggridges, enrichplot
# AI:   httr2, curl
```

---

## Usage

1. **Load Data** -- Upload a DIA-NN `report.parquet` (or `.tsv`) output file, or click "Load Example Data" for a demo HeLa dataset
2. **Assign Groups & Run** -- Auto-guess groups from filenames or manually assign; optionally add covariates (batch, etc.); click "Run Pipeline" to execute DPC-CN normalization, maxLFQ quantification, and limma DE
3. **Explore Results** -- Data Overview, QC, DE Dashboard (Volcano/Table/PCA/Robust Changes), Phospho, GSEA, MOFA2, AI Analysis, XIC Viewer (local/HPC)
4. **Export** -- Download reproducibility log (.R), save session (.rds), export tables and plots

---

## Methodology

| Step | Method |
|------|--------|
| **Normalization** | Data Point Correspondence - Cyclic Normalization (DPC-CN) via `limpa::dpcCN()` |
| **Quantification** | Modified maxLFQ (precursor-to-protein rollup) via `limpa::dpcQuant()` |
| **DE model** | Linear model fit via `limpa::dpcDE()` + `limma::contrasts.fit()` |
| **Moderation** | Empirical Bayes moderated *t*-statistics via `limma::eBayes()` |
| **FDR** | Benjamini-Hochberg adjusted *p*-values |
| **Phospho DE** | Same limma pipeline at the phosphosite level (independent from protein-level) |

**Key Citations:**
- **limpa** -- Bioconductor package for DIA proteomics ([bioconductor.org/packages/limpa](https://bioconductor.org/packages/limpa/))
- **limma** -- Ritchie ME et al. (2015) *Nucleic Acids Res* 43(7):e47 ([doi:10.1093/nar/gkv007](https://doi.org/10.1093/nar/gkv007))
- **DIA-NN** -- Demichev V et al. (2020) *Nat Methods* 17:41-44 ([doi:10.1038/s41592-019-0638-x](https://doi.org/10.1038/s41592-019-0638-x))
- **MOFA2** -- Argelaguet R et al. (2020) *Genome Biol* 21:111 ([doi:10.1186/s13059-020-02015-1](https://doi.org/10.1186/s13059-020-02015-1))
- **KSEA** -- Wiredja DD et al. (2017) *Bioinformatics* 33:3489-3491; Casado P et al. (2013) *Sci Signaling* 6:rs6
- **clusterProfiler** -- Wu T et al. (2021) *Innovation* 2(3):100141

---

## Resources

- **Project Website:** [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/)
- **Video Tutorials:** [UC Davis Proteomics YouTube](https://www.youtube.com/channel/UCpulhf8gl-HVxACyJUEFPRw)
- **Training:** [Hands-On Proteomics Short Course](https://proteomics.ucdavis.edu/events/hands-proteomics-short-course)
- **Core Facility:** [proteomics.ucdavis.edu](https://proteomics.ucdavis.edu)

---

## License

This project is open source. See repository for license details.

## Contributing

Issues and pull requests welcome! See [CLAUDE.md](CLAUDE.md) for development documentation.

**Developer:** Brett Phinney, UC Davis Proteomics Core Facility | **Contact:** [GitHub Issues](https://github.com/bsphinney/DE-LIMP/issues)

## Example Data

Demo dataset: **Affinisep vs Evosep** SPE column comparison using 50 ng Thermo HeLa protein digest standard (DIA, Orbitrap). Available at [github.com/bsphinney/DE-LIMP/releases](https://github.com/bsphinney/DE-LIMP/releases).
