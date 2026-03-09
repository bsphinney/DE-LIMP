<img src="https://github.com/user-attachments/assets/2aeb5863-2c10-4faa-99e8-25835a9e9330" align="left" width="150" style="margin-right: 20px;" alt="DE-LIMP Logo" />

# DE-LIMP: Differential Expression & Limpa Proteomics

Find which proteins are significantly different between your experimental conditions -- upload a DIA-NN output file and get interactive volcano plots, heatmaps, pathway enrichment, and AI-powered interpretation, all without writing code.

Built on R Shiny with the [limpa](https://bioconductor.org/packages/limpa/) pipeline (DPC-CN normalization, maxLFQ quantification) and [limma](https://bioconductor.org/packages/limma/) empirical Bayes moderated *t*-statistics with Benjamini-Hochberg FDR correction.

**Input:** DIA-NN `report.parquet` | **Not for:** DDA data, TMT/iTRAQ, Spectronaut/MaxQuant output

> **Not sure if your data is DIA?** If your core facility used DIA-NN to process your samples, you have DIA data. Look for a `report.parquet` file in your results folder. If your data was processed with MaxQuant, Spectronaut, or Proteome Discoverer, or if you used isobaric labels (TMT, iTRAQ), DE-LIMP is not the right tool.

<br clear="left"/>

**Try it now:** [huggingface.co/spaces/brettsp/de-limp-proteomics](https://huggingface.co/spaces/brettsp/de-limp-proteomics) -- no installation required

**Project Website:** [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/) | **Docs:** [USER_GUIDE.md](USER_GUIDE.md) | [CLAUDE.md](CLAUDE.md)

---

## What's New in v3.5.0

**Run Comparator** -- Compare two analyses of the same dataset across tools. Supports DE-LIMP vs DE-LIMP, Spectronaut, or FragPipe. Four diagnostic layers (Settings Diff, Protein Universe, Quantification, DE Concordance) with a 7-rule hypothesis engine that explains *why* each discordant protein disagrees. Optional DIA-NN log upload enriches the comparison with search-derived parameters. Optional MOFA2 decomposition identifies hidden variance patterns.

**Search & Analysis History** -- Full audit trail for DIA-NN searches (26 parameters logged) and pipeline runs. Import Settings or Import Results from past searches. Cross-reference links between search and analysis history. Project-based grouping with summary cards.

**Chromatography QC** -- Extract TIC traces from timsTOF .d files before committing to long searches. Faceted, overlay, and metrics views with per-run diagnostics (shape deviation, RT shift, loading anomaly, late elution, and more). Flag dead injections and carryover before wasting HPC hours.

**Smart HPC Job Submission** -- Per-user CPU limit detection (queries SLURM QOS), auto-switches partition when at capacity. FASTA database library with auto-upload to HPC. Path validation prevents local-only FASTA files from reaching HPC jobs.

**Previous highlights:** v3.1 UI overhaul (dark navbar, accordion sidebar, DE Dashboard sub-tabs, Core Facility Mode). v3.0 Multi-Omics MOFA2, Docker search, phosphoproteomics (KSEA, motifs), GSEA (BP/MF/CC/KEGG), all-contrast AI summary.

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
> **Requires a free Gemini API key.** Get one at [Google AI Studio](https://aistudio.google.com/) and paste it into the DE-LIMP sidebar.

- **AI Summary** -- Analyzes all contrasts simultaneously, identifying top DE proteins per comparison, cross-comparison biomarkers, and CV-based stability metrics. AI Summary sends only summary statistics (protein names, logFC, adj.P.Val); Data Chat sends per-sample expression data for top DE proteins to enable interactive Q&A
- **Export for Claude** -- Download your complete analysis as a .zip optimized for deep analysis with Claude, ChatGPT, or other AI assistants (includes DE results, expression matrix, QC metrics, GSEA, methods text, and more)
- **AI Summary HTML Export** -- Styled standalone HTML report with gradient header and markdown formatting, suitable for sharing with collaborators
- **Interactive Data Chat** -- Conversational interface with Google Gemini, auto-injecting QC stats and 100-800 top DE proteins as context. Phospho context (top 20 sites + KSEA kinase results) auto-included when phospho analysis is active
- **Interactive AI + plot connection** -- Select proteins in volcano/table to set AI context; AI can highlight proteins in plots via `[[SELECT: protein1; protein2]]` syntax
- **Auto-Analyze** button for one-click dataset analysis; **Save Chat** to download conversation as plain text
- Auto-generated methodology text for methods sections

### Run Comparator
- **Cross-tool comparison** -- Compare your DE-LIMP analysis against a second DE-LIMP run, Spectronaut export, or FragPipe output to understand how tool choice affects your results
- **4 diagnostic layers** -- Settings Diff (parameter-by-parameter comparison), Protein Universe (overlap analysis), Quantification (log2 intensity correlation, per-sample concordance, systematic bias detection), DE Concordance (3x3 Up/Down/NS matrix, volcano overlay, discordant protein table)
- **7-rule hypothesis engine** -- For each discordant protein, assigns a tool-aware hypothesis explaining *why* the tools disagree (direction reversal, normalization offset, variance estimation, missing values, peptide count, FC magnitude, or borderline significance)
- **Optional DIA-NN log upload** -- Enrich Mode A comparisons with search-derived parameters (pg-level quantification, proteoforms, library precursor counts, pipeline step)
- **Optional MOFA2 decomposition** -- Treats the two runs as views and decomposes joint variance to find hidden patterns among discordant proteins
- **AI integration** -- Tool-aware Gemini prompt and Claude ZIP export for deeper analysis

### Chromatography QC
- **Pre-search quality check** -- Extract TIC traces from timsTOF .d files *before* committing to hours-long DIA-NN searches
- **Three views** -- Faceted panels (per-run with median overlay), Overlay (all runs normalized 0-1 on one axis), Metrics (AUC bar chart + diagnostics table)
- **Automated diagnostics** -- Shape deviation (Pearson r vs median trace), RT shift, loading anomaly (AUC outlier), file size outlier, late elution, elevated baseline, narrow gradient
- **SSH support** -- SCP downloads analysis.tdf from remote .d directories, extracts locally

### DIA-NN Search Integration
- **Three backends** -- Local, Docker, and HPC (SSH/SLURM)
- **Parallel 5-step SLURM pipeline** -- Optimized search with dependency chaining and array jobs for maximum HPC throughput
- **Spectral library caching** -- Reuse predicted libraries across searches to save compute time
- **Custom FASTA sequences** -- Add custom protein sequences inline when submitting searches
- **Smart partition selection** -- Detects per-user SLURM CPU limits, auto-switches to public queue when at capacity
- **FASTA database library** -- Shared catalog with auto-upload to HPC, fragment m/z range tracking, path validation
- **Cluster resource indicator** -- Real-time HPC CPU usage monitoring with traffic-light display (green/yellow/red)
- **Windows Docker** -- `docker compose up` runs DE-LIMP + DIA-NN with zero R installation ([guide](WINDOWS_DOCKER_INSTALL.md))
- **UniProt FASTA download** -- Search and download proteome databases directly; 6 bundled contaminant libraries
- **Non-blocking job queue** -- Submit multiple searches, results auto-load on completion
- **Phospho mode** -- Auto-configures DIA-NN for phospho analysis (STY modification, `--phospho-output`)
- **Organized search logs** -- SLURM `.out`/`.err` and local `.log` files written to `{output_dir}/logs/`

> **DIA-NN License:** DIA-NN is developed by [Vadim Demichev](https://github.com/vdemichev/DiaNN) and is free for academic/non-commercial use. It is not open source and cannot be redistributed. DE-LIMP does not bundle DIA-NN. See the [DIA-NN license](https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md).

### Core Facility Mode *(Optional)*
- Staff YAML profiles auto-fill SSH, SLURM, and instrument settings
- SQLite job tracking with searchable history (6 filters), one-click result loading and report generation
- Instrument QC dashboard with protein/precursor/TIC trends and control lines
- Quarto HTML reports with QC bracket, volcanos, DE stats, and top proteins

> *Activated by setting `DELIMP_CORE_DIR`. Not visible on standard installations.*

### Session Management & History
- **Search History** -- Full audit trail for every DIA-NN search (26 parameters). Import Settings to reuse parameters; Import Results to load completed search output directly. View Log shows search metadata. Cross-reference links to Analysis History.
- **Analysis History & Projects** -- Track every pipeline run with expandable detail rows. Assign analyses to projects for organized grouping with summary cards.
- **About tab** -- Community stats dashboard with GitHub stars, forks, visitors, and clones (14-day trend sparklines), GitHub Discussions feed, version info, and project links
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

1. **Load Data** -- Upload a DIA-NN `report.parquet` output file, or click "Load Example Data" for a demo HeLa dataset
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
- **Discussions:** [github.com/bsphinney/DE-LIMP/discussions](https://github.com/bsphinney/DE-LIMP/discussions) -- Q&A, feature ideas, and announcements
- **Video Tutorials:** [UC Davis Proteomics YouTube](https://www.youtube.com/channel/UCpulhf8gl-HVxACyJUEFPRw)
- **Training:** [Hands-On Proteomics Short Course](https://proteomics.ucdavis.edu/events/hands-proteomics-short-course)
- **Core Facility:** [proteomics.ucdavis.edu](https://proteomics.ucdavis.edu)

---

## License

This project is open source. See repository for license details.

## Contributing

Issues, pull requests, and [Discussions](https://github.com/bsphinney/DE-LIMP/discussions) welcome! See [CLAUDE.md](CLAUDE.md) for development documentation.

**Developer:** Brett Phinney, UC Davis Proteomics Core Facility | **Contact:** [GitHub Issues](https://github.com/bsphinney/DE-LIMP/issues)

## Example Data

Demo dataset: **Affinisep vs Evosep** SPE column comparison using 50 ng Thermo HeLa protein digest standard (DIA, Orbitrap). Available at [github.com/bsphinney/DE-LIMP/releases](https://github.com/bsphinney/DE-LIMP/releases).
