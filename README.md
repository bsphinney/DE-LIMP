<img src="https://github.com/user-attachments/assets/2aeb5863-2c10-4faa-99e8-25835a9e9330" align="left" width="150" style="margin-right: 20px;" alt="DE-LIMP Logo" />

# DE-LIMP: Differential Expression & Limpa Proteomics

Find which proteins are significantly different between your experimental conditions -- upload a DIA-NN output file and get interactive volcano plots, heatmaps, pathway enrichment, and AI-powered interpretation, all without writing code.

Built on R Shiny with the [limpa](https://bioconductor.org/packages/limpa/) pipeline for normalization and protein quantification, and [limma](https://bioconductor.org/packages/limma/) for statistical testing with FDR correction. See [USER_GUIDE.md](USER_GUIDE.md#glossary) for methodology details.

**Input:** DIA-NN `report.parquet` | **Not for:** DDA data, TMT/iTRAQ, Spectronaut/MaxQuant output

> **Not sure if your data is DIA?** If your core facility used DIA-NN to process your samples, you have DIA data. Look for a `report.parquet` file in your results folder. If your data was processed with MaxQuant, Spectronaut, or Proteome Discoverer, or if you used isobaric labels (TMT, iTRAQ), DE-LIMP is not the right tool.

<br clear="left"/>

**Try it now:** [huggingface.co/spaces/brettsp/de-limp-proteomics](https://huggingface.co/spaces/brettsp/de-limp-proteomics) -- no installation required

**Project Website:** [bsphinney.github.io/DE-LIMP](https://bsphinney.github.io/DE-LIMP/) | **Docs:** [USER_GUIDE.md](USER_GUIDE.md) | [CLAUDE.md](CLAUDE.md)

---

## What's New in v3.11.0

**Proteogenomics — Build sample-specific FASTA databases from matched RNA-seq** -- Upload RNA-seq fastq files or SRA/SLIMS accessions, and DE-LIMP runs a complete HPC pipeline (fastp → bowtie2 rRNA filter → STAR → stringtie → merge → gffcompare → gffread → TransDecoder) to predict novel ORFs and alternative splice variants alongside your canonical proteome. Final FASTA auto-assembles with UniProt or NCBI RefSeq, registers in the FASTA library catalog, and appears as "Proteogenomics DBs" in the main search page's FASTA Database dropdown. Reference genomes: Human, Mouse, Pig, Rat, Arabidopsis. Multi-user workflow with "Restore from Hive" and "Discover from Hive" buttons for shared catalog discovery.

**Two analysis pipelines, one app** -- Choose between **DPC-Quant + limma** (limpa's detection-probability model, default) and **MaxLFQ + limma** (paper-faithful Moschem et al. 2025 implementation). The pipeline you pick is recorded in the dataset itself; methods text, AI prompts, exports, and the Reproducibility log all describe whichever pipeline actually ran -- no hardcoded "DPC-Quant" strings anywhere.

**QuantUMS quality filters** -- Optional pre-filtering of DIA-NN precursors by `eQ`, `qQ`, and `pgQ` quality scores (Moschem 2025). Applied at pipeline run-time on the parquet, with a waterfall showing how many precursors and proteins survived each filter. Defaults to off; opt-in checkbox lets you feed filtered data through limpa's DPC-Quant if you want to test the combination.

**On/Off Proteins panel** -- New sub-tab in DE Dashboard surfaces proteins detected in ≥N samples of one condition AND zero samples of the other. limma assigns these `NA` logFC, so they're invisible in the volcano -- the new panel makes them findable.

**Coverage filter** -- Drop proteins with fewer than X non-NA samples before limma fits the model (UC Davis Bioinformatics Core convention). Live waterfall shows protein retention. Available in MaxLFQ + limma mode.

**Run Comparator pipeline-aware** -- Cross-tool DE comparison (DE-LIMP vs DE-LIMP / Spectronaut / FragPipe) now reads the pipeline descriptor on both sides. Hypothesis-engine rules (rollup, normalization, peptide rules) emit the correct contrast based on whether each side used MaxLFQ or DPC-Quant.

**Export Complete Analysis is now a true superset** -- Single ZIP at Output > Export Complete Analysis includes DE results, QC metrics, phospho results (when present), expression matrix, detection matrix, quartile profiles, variable proteins, contaminant summary, search_info.md, methods.txt, parameters.txt, reproducibility_log.R + sessionInfo, session.rds, and a DE-aware **PROMPT.md** for LLM analysis. Every section is wrapped in `safe_section()` and a **MANIFEST.txt** records what was included vs skipped and why. The three redundant "Export for Claude" buttons (Data Explorer, AI Summary, AI Chat) are now consolidated into this single download.

**FASTA picker** -- Scanning a shared FASTA directory (e.g. `/quobyte/proteomics-grp/de-limp/fasta`) used to silently combine every `.fasta` it found. Now: 1 file → use directly; ≥2 files → checkbox modal so you pick exactly which one(s) to use. Local browse and SSH scan both fixed.

**Provenance block** -- Exports include parquet MD5, full sessionInfo, app version, and pipeline label so reanalysis is bit-reproducible.

**Previous highlights** (v3.10): Two analysis pipelines, QuantUMS quality filters, On/Off proteins, Coverage filter, FASTA picker, provenance block. (v3.7): NCBI Proteome Download with gene mapping, Contaminant Analysis with keratin flagging, Data Explorer (quartile + scatter), SSH File Browser, Load from HPC, WSL2 Launcher for Windows, No-Replicates Mode, SSH Auto-Connect, Environment Badge.

**Earlier highlights** (v3.5): Run Comparator, Search & Analysis History, Chromatography QC, smart HPC partitions. (v3.1): UI overhaul, Core Facility Mode. (v3.0): MOFA2, Docker search, phosphoproteomics, GSEA.

See [CHANGELOG.md](CHANGELOG.md) for full release history.

---

## Key Features

### Analysis & Visualization
- **Volcano Plots** -- Interactive (Plotly), click or box-select proteins to highlight across all views; all pairwise contrasts available
- **Heatmaps** -- Z-score heatmaps of selected or significant proteins (ComplexHeatmap)
- **Contaminant Analysis** -- Summary cards, per-sample stacked bar chart, top contaminants table with keratin flagging, and contaminant heatmap; Signal Distribution and Expression Grid also highlight contaminants
- **Data Explorer** -- Quartile-based abundance profiles and sample-sample scatter plots for exploring data without DE analysis
- **QC Sample Metrics** -- Faceted trend plot (Precursors, Proteins, MS1 Signal, Data Completeness) with LOESS smoother for drift detection and group average lines
- **MDS & DPC Plots** -- Sample clustering and normalization diagnostics
- **Covariates** -- Include batch, sex, diet, or custom covariates in the linear model
- **XIC Chromatogram Viewer** -- Fragment-level chromatogram validation, MS2 intensity alignment (Spectronaut-style), ion mobility/mobilogram support for timsTOF, DIA-NN v1/v2 formats (local/HPC only)
- **CV Analysis (Robust Changes)** -- Identify highly reproducible DE proteins via coefficient of variation analysis across replicates

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
- **SSH file browser** -- Visual directory browser for navigating remote HPC filesystems with clickable breadcrumbs, color-coded entries, and file type filtering
- **SSH auto-connect** -- Automatically connects to HPC on startup when an SSH key is detected; environment badge shows deployment mode
- **UniProt FASTA download** -- Search and download proteome databases directly; 6 bundled contaminant libraries
- **NCBI proteome download** -- Download RefSeq protein FASTA from NCBI Datasets with automatic gene symbol mapping for non-model organisms
- **Load from HPC** -- One-click button to browse, download, and analyze completed search results from the cluster
- **Spectral library caching** -- Reuse predicted libraries across searches to save compute time
- **Custom FASTA sequences** -- Add custom protein sequences inline when submitting searches
- **Smart partition selection** -- Detects per-user SLURM CPU limits, auto-switches to public queue when at capacity
- **FASTA database library** -- Shared catalog with auto-upload to HPC, fragment m/z range tracking, path validation
- **Cluster resource indicator** -- Real-time HPC CPU usage monitoring with traffic-light display (green/yellow/red)
- **Windows WSL launcher (recommended)** -- One-click `.bat` runs DE-LIMP + DIA-NN natively in WSL2 Ubuntu, zero R install on Windows ([guide](WINDOWS_WSL_INSTALL.md))
- **Windows Docker launcher (alternative)** -- One-click `.bat` for users who already run Docker Desktop ([guide](WINDOWS_DOCKER_INSTALL.md))
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
- **Unified activity log** -- Single audit trail for all DIA-NN searches and pipeline runs, with remote activity log via SSH for multi-user visibility
- **Search History** -- Full audit trail for every DIA-NN search (26 parameters). Import Settings to reuse parameters; Import Results to load completed search output directly. View Log shows search metadata. Cross-reference links to Analysis History.
- **Analysis History & Projects** -- Track every pipeline run with expandable detail rows. Assign analyses to projects for organized grouping with summary cards.
- **About tab** -- Community stats dashboard with GitHub stars, forks, visitors, and clones (14-day trend sparklines), GitHub Discussions feed, version info, and project links
- **No-replicates mode** -- Quantification without DE for n=1 experiments; PCA, Expression Grid, and Data Explorer still available
- Save/load full analysis state as `.rds`; export reproducibility R code log
- One-click example data (Affinisep vs Evosep comparison)
- Group assignment templates (CSV export/import)
- Embedded proteomics resources, UC Davis Proteomics videos, short course links

---

## Which Installation Should I Use?

| Platform | Method | DIA-NN Search? | Guide |
|----------|--------|----------------|-------|
| **Any (just exploring)** | Web browser | No | [Hugging Face](https://huggingface.co/spaces/brettsp/de-limp-proteomics) |
| **Windows (recommended)** | **WSL2 + native Linux R** | Yes (local + HPC) | **[WINDOWS_WSL_INSTALL.md](WINDOWS_WSL_INSTALL.md)** |
| **Windows (alternative)** | Docker + SSH to HPC | Yes (local + HPC) | [WINDOWS_DOCKER_INSTALL.md](WINDOWS_DOCKER_INSTALL.md) |
| **Mac / Linux** | R/RStudio (native) | Via HPC | See [Installation](#installation) below |
| **HPC cluster** | Apptainer/Singularity | Via SLURM | [HPC_DEPLOYMENT.md](HPC_DEPLOYMENT.md) |

> **Why WSL over Docker on Windows?** WSL2 runs a real Ubuntu environment natively in Windows 10/11. R and Bioconductor packages compile cleanly, SSH keys use normal Unix permissions (no CRLF / 0600 / missing-newline gymnastics), no 9p filesystem perf tax, and no .NET download blocked by corporate firewalls. The Docker path still works and is kept for users who already depend on it, but new installs should start with WSL.

---

## Installation

**Requirements:** R 4.5 or newer (for limpa). Everything else is installed automatically the first time you launch the app.

### Mac -- step-by-step (first-time users)

If this is your first time running an R/Shiny app, do these four steps in order. Skip ahead to **Already have R, RStudio, and Git?** below if you don't need the walkthrough.

**1. Install R**

Download the right installer for your Mac from CRAN: <https://cloud.r-project.org/bin/macosx/>

- Apple Silicon (M1, M2, M3, M4): pick the **`-arm64.pkg`** file.
- Intel Mac: pick the **`-x86_64.pkg`** file (also labelled "older Macs").

Open the `.pkg` file and follow the prompts. Not sure which chip you have? Click the Apple menu → "About This Mac". "Apple M*x*" = arm64; anything starting with "Intel" = x86_64.

**2. Install RStudio (recommended)**

Download from <https://posit.co/download/rstudio-desktop/> and drag RStudio to `/Applications`. RStudio is the friendliest way to run R on a Mac. (VS Code with the R extension also works if you prefer it.)

**3. Install Git**

The simplest way: open the **Terminal** app (Cmd-Space, type "Terminal", Enter) and run:

```bash
git --version
```

If Git is missing, macOS will pop up a dialog asking to install the Xcode Command Line Tools -- click "Install" and wait. That installs Git, no Apple Developer account needed.

**4. Download DE-LIMP and launch it**

Still in Terminal:

```bash
cd ~/Documents
git clone https://github.com/bsphinney/DE-LIMP.git
```

This creates `~/Documents/DE-LIMP/`. Now open RStudio. In the R console (the bottom-left pane), paste:

```r
shiny::runApp('~/Documents/DE-LIMP', port = 3838, launch.browser = TRUE)
```

The first launch takes 5--15 minutes -- DE-LIMP installs ~30 R packages from CRAN and Bioconductor automatically. Subsequent launches start in seconds. When it's done, your browser opens to the DE-LIMP app.

**Updating later:**

```bash
cd ~/Documents/DE-LIMP
git pull
```

Then re-launch from RStudio. Updating the app code never needs you to reinstall R or RStudio.

### Already have R, RStudio, and Git?

```bash
git clone https://github.com/bsphinney/DE-LIMP.git
cd DE-LIMP
```

```r
shiny::runApp('.', port = 3838, launch.browser = TRUE)
```

All dependencies install automatically on first run:

- **Core:** shiny, bslib, plotly, DT, rhandsontable, shinyjs, dplyr, tidyr, stringr, readr, arrow, ggplot2, ggrepel, ggridges
- **Stats:** limpa, limma, ComplexHeatmap, clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, AnnotationDbi, enrichplot, KSEAapp, ggseqlogo, MOFA2, basilisk, callr
- **AI / network:** httr2, curl

### Troubleshooting first launch

- **"Package 'limpa' is missing" or "Missing required packages"** -- DE-LIMP installs missing packages automatically, but right after a brand-new R release (e.g. you just installed R 4.6) the BiocManager helper can lag a few weeks behind in knowing which Bioconductor branch pairs with the new R. DE-LIMP v3.7.3+ falls back to a direct Bioconductor download in that case, so just rerun the launch command once and the install finishes. Make sure Wi-Fi is on; the first run needs internet.
- **"R version: 4.x.y (NEED: 4.5+)"** -- you have an older R. Upgrade from <https://cloud.r-project.org/bin/macosx/> and rerun.
- **Permission denied when installing packages** -- you're trying to install into the system library. Either run RStudio normally (it sets up a personal library automatically) or, in a Terminal R session, type `dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE, showWarnings = FALSE)` once.
- **Port 3838 already in use** -- another R session is already serving the app. In RStudio, click "Session → Restart R", then relaunch. Or change the port: `shiny::runApp(..., port = 3839)`.

---

## Usage

1. **Load Data** -- Upload a DIA-NN `report.parquet` output file, or click "Load Example Data" for a demo HeLa dataset
2. **Assign Groups & Run** -- Auto-guess groups from filenames or manually assign; optionally add covariates (batch, etc.); click "Run Pipeline" to execute DPC-CN normalization, DPC-Quant protein quantification, and limma DE
3. **Explore Results** -- Data Overview, QC, DE Dashboard (Volcano/Table/PCA/CV Analysis), Phospho, GSEA, MOFA2, AI Analysis, XIC Viewer (local/HPC)
4. **Export** -- Download reproducibility log (.R), save session (.rds), export tables and plots

---

## Methodology

| Step | Method |
|------|--------|
| **Normalization** | Data Point Correspondence - Cyclic Normalization (DPC-CN) via `limpa::dpcCN()` |
| **Quantification** | DPC-Quant (Detection Probability Curve Quantification): precursor-to-protein rollup via probabilistic missing-value modelling, via `limpa::dpcQuant()` |
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
