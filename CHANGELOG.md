# Changelog

All notable changes to DE-LIMP will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [3.1.1] - 2026-02-26

### Fixed
- **Volcano plot significance mismatch**: Y-axis used raw P.Value but coloring used adj.P.Val, causing ~746 proteins to appear above the threshold line while colored as "Not Sig". Now the horizontal dashed line is drawn at the raw P.Value corresponding to adj.P.Val = 0.05, so the line and coloring agree visually.
- **Volcano significance coloring**: Removed logFC cutoff from significance determination — all proteins with adj.P.Val < 0.05 are now colored red. logFC vertical lines remain as visual guides.
- **Default logFC cutoff**: Changed from 1.0 (2-fold change) to 0.6 (~1.5-fold change).
- **CV Analysis cards compressed**: Plotly annotation cards above the scatter plot kept overlapping and rendering incorrectly in bslib sub-tabs. Replaced with a ggplot subtitle showing per-group median CV and % proteins below 20% CV.
- **CV Analysis tab layout**: Wrapped in scrollable div with min-height to prevent bslib from compressing the scatter plot. Content stacks vertically with overflow-y scroll.

### Added
- **Volcano DE protein count**: Info box now shows "78 DE proteins (X up, Y down)" so users can immediately see the count.
- **Export Data panel**: New "Export Data" tab in the Output dropdown with prominent buttons for downloading Results CSV and CV Analysis CSV.
- **AI Summary HTML export**: "Export Report" button on AI Summary tab generates a styled standalone HTML report with gradient header, formatted tables, and print-friendly CSS.
- **Docker update scripts**: `update_docker.sh` (bash) and `update_docker.ps1` (PowerShell) for one-command pull + rebuild on Windows/Mac/Linux.

## [3.1.0] - 2026-02-24

### Added
- **UI Overhaul**:
  - `page_navbar()` layout with dark navbar (`#2c3e50`), hover-activated dropdown menus with smooth animations
  - Dropdown section labels (Setup / Results / AI) in the Analysis menu via JS injection
  - Active tab teal underline indicator (`--flatly-success`)
  - `nav_spacer()` + native `nav_item()` for gear icon (replaces JS injection)

- **Accordion Sidebar**:
  - Three collapsible `accordion_panel()` sections: Upload Data (open by default), Pipeline Settings, AI Chat
  - Phospho and XIC sections in separate conditional accordions outside the main group
  - Session buttons and core facility controls outside accordion for always-visible access

- **DE Dashboard Sub-tabs**:
  - Replaced grid+accordion layout with `navset_card_tab(id = "de_dashboard_subtabs")`
  - Four sub-tabs: Volcano (with heatmap below), Results Table, PCA, Robust Changes
  - PCA moved from Data Overview into DE Dashboard

- **Core Facility Mode**:
  - Activated by `DELIMP_CORE_DIR` env var pointing to directory with `staff.yml`
  - SQLite database (`delimp.db`) with 4 tables: searches, qc_runs, reports, templates
  - WAL mode for concurrent SQLite access
  - Staff YAML config auto-fills SSH host, username, key path, SLURM account/partition
  - Search DB tab: full-width job history with 6 filters (text, lab, status, staff, instrument, LC method)
  - Instrument QC dashboard: protein/precursor/TIC trend plots with ±2SD control lines, instrument filter, date range
  - Quarto report generation: standalone HTML with metadata, QC bracket, volcanos, DE stats, top proteins
  - Template system for saving/loading search presets
  - New files: `R/helpers_facility.R`, `R/server_facility.R`, `report_template.qmd`, `seed_test_db.R`

### Changed
- `R/ui.R`: Rewritten outer wrapper from `page_sidebar()` + `navset_card_tab()` to `page_navbar()` with direct nav items
- `R/ui.R`: Sidebar width reduced from 320 to 300px
- `R/ui.R`: Data Overview reduced from 6 to 5 sub-tabs (PCA moved to DE Dashboard)
- `R/ui.R`: Removed `.de-dashboard-grid` CSS class and accordion heatmap wrapper
- `CLAUDE.md`: Condensed from ~500 to ~130 lines; detailed patterns moved to `docs/PATTERNS.md`, TODOs to `docs/TODO.md`
- `README_GITHUB.md`, `README_HF.md`, `USER_GUIDE.md`: Updated for v3.1 features
- App version bumped to v3.1

### Fixed
- **Navbar text invisible on dark background**: Flatly theme renders dark text on `navbar-inverse`. Fixed with CSS `!important` overrides for white text on `.navbar .nav-link` and `.navbar-brand`.
- **Hidden tab fragments visible**: `nav_hide()` on `page_navbar` leaves letter fragments. Fixed with `width: 0 !important; overflow: hidden !important` on hidden nav items.
- **bslib deprecation warning**: `page_navbar(bg=...)` deprecated in bslib 0.9.0+. Changed to `navbar_options = navbar_options(bg = "#2c3e50")`.

## [3.0.0] - 2026-02-20

### Added
- **Multi-Omics MOFA2 Integration**:
  - New **Multi-Omics MOFA2** tab for unsupervised integration of 2-6 data views using MOFA2
  - Dynamic view cards: add/remove views, smart RDS parser (DE-LIMP sessions, limma EList/MArrayLM, matrices, data frames)
  - CSV/TSV/Parquet matrix upload with auto-log2 detection
  - Phospho tab integration as data source
  - Sample matching with overlap statistics and color-coded status
  - MOFA training via `callr::r()` subprocess (isolates basilisk/Python from Shiny's event loop)
  - 5 results tabs: Variance Explained heatmap, Factor Weights browser (plotly), Sample Scores scatter, Top Features table (DT), Factor-DE Correlation
  - **Mouse Brain example dataset** (2-view): proteomics + phosphoproteomics, 16 samples
  - **TCGA Breast Cancer example dataset** (3-view): mRNA + miRNA + protein, 150 samples, 3 subtypes
  - Session save/load, methodology text, reproducibility logging
  - New packages: `MOFA2`, `basilisk`, `callr`

- **DIA-NN Docker Local Backend**:
  - `Dockerfile.search`: Multi-stage build embeds DIA-NN binaries from user's pre-built `diann:2.0` image into DE-LIMP container
  - `docker-compose.yml`: One-command deployment (`docker compose up`) with data volume mounts
  - DIA-NN runs inside the container as a background process — no Docker-in-Docker required
  - Windows Docker Desktop support: intermediate files written to container-internal `/tmp` to avoid FUSE layer issues with large speclib writes
  - `build_diann_docker.sh` / `build_diann_docker.ps1`: User-facing scripts to build DIA-NN image (license compliance)

- **Windows Docker Deployment**:
  - `WINDOWS_DOCKER_INSTALL.md`: Step-by-step guide for Windows users (Docker Desktop + WSL2)
  - Zero R installation required — everything runs in Docker
  - `MACOS_DOCKER_INSTALL.md` and `LINUX_DOCKER_INSTALL.md`: Platform-specific guides

### Changed
- `R/ui.R`: Renamed "Multi-View Integration" tab to "Multi-Omics MOFA2"
- `R/helpers_search.R`: Docker command uses `--entrypoint sh -c` wrapper to write intermediate files to `/tmp/diann_work` then copy final outputs to mounted volume
- `Dockerfile` and `Dockerfile.search`: Added MOFA2, basilisk, callr packages; pre-initialize basilisk Python env at build time
- `docs/index.html`: Comprehensive update for v3.0 features, deployment options, corrected Getting Started
- `README_GITHUB.md`, `README_HF.md`, `USER_GUIDE.md`: Added MOFA2 documentation, fixed Quick Start
- App version bumped to v3.0

### Fixed
- **R session crash during MOFA2 training**: basilisk's Python subprocess conflicts with Shiny's httpuv event loop. Fixed by running entire MOFA pipeline in isolated subprocess via `callr::r()`
- **Variance Explained heatmap empty**: `r2_per_factor` is a matrix (factors × views), not a nested named list. Fixed iteration logic.
- **Factor-DE Correlation tab blank**: Controls required `values$fit` which is NULL when using test data. Now shows info message when no DE results available.
- **Docker DIA-NN "Could not save" on Windows**: Large intermediate `.predicted.speclib` files fail on Windows Docker FUSE layer. Fixed by writing to container-internal `/tmp` then copying final outputs.
- **Variance explained download handler**: Same matrix iteration bug as the heatmap render.

## [2.5.0] - 2026-02-18

### Added
- **GSEA Expansion — Multi-Database Enrichment**:
  - **Four enrichment databases**: GO Biological Process (BP), GO Molecular Function (MF), GO Cellular Component (CC), and KEGG Pathways
  - **Ontology selector**: Dropdown to switch between BP/MF/CC/KEGG on the GSEA tab
  - **Per-ontology caching**: Results cached per database; switching back loads instantly without re-computation
  - **Contrast indicator**: Shows active contrast with stale-results warning when contrast changes
  - **UniProt API organism detection**: Queries `rest.uniprot.org` to determine organism from accession when no suffix is present — works automatically for human, mouse, rat, and 9 other species
  - **Robust ID mapping**: Handles multiple protein ID formats (pipe-separated, isoform suffixes, organism suffixes); fallback from UNIPROT to SYMBOL ID types
  - **KEGG organism mapping**: Supports 11 species with automatic organism code detection
  - **Dynamic plot titles**: All GSEA plots show which database was used
  - Updated info modals with all 4 database descriptions
  - Session save/load for GSEA cache, last contrast, and organism DB

- **AI Summary — All Comparisons Analysis**:
  - AI Summary now analyzes **all contrasts** simultaneously, not just the selected one
  - Cross-comparison biomarker detection: identifies proteins significant in ≥2 comparisons
  - Enhanced prompt with 5 structured sections: Overview, Key Findings Per Comparison, Cross-Comparison Biomarkers, High-Confidence Biomarker Insights, Biological Interpretation
  - Adaptive token budget: top 30/20/10 proteins per contrast (scales with number of contrasts)
  - `?` info modal explaining what data is/isn't sent to Gemini API

- **MDS Plot Coloring**:
  - Color MDS plot by Group, Batch, Covariate1, or Covariate2
  - Colorblind-friendly Okabe-Ito palette
  - Dynamic dropdown updates when metadata changes

- **Complete Dataset Export**:
  - Download button on Dataset Summary tab
  - Exports: protein IDs, gene symbols, DE stats for ALL contrasts (suffixed columns), per-sample expression values, metadata as header comment rows

- **Phosphoproteomics Phase 2 — Kinase Activity & Motif Analysis**:
  - **KSEA kinase-substrate enrichment analysis**: Infers upstream kinase activity from phosphosite fold-changes using `KSEAapp` CRAN package with PhosphoSitePlus + NetworKIN database
  - **KSEA bar plot**: Horizontal bar chart of kinase z-scores (top 15 activated + top 15 inhibited), colored by direction, with substrate count annotations
  - **KSEA results table**: Filterable, sortable DT datatable with kinase gene, z-score, FDR, substrate count. Downloadable as CSV.
  - **Sequence logo motif analysis**: Displays amino acid enrichment around regulated phosphosites using `ggseqlogo`. Separate logos for up-regulated and down-regulated sites.
  - **FASTA upload**: Sidebar file input for protein FASTA. Parses UniProt-format headers to extract accessions. Enables accurate flanking sequence extraction for motif analysis.
  - New "Kinase Activity" and "Motif Analysis" tabs in phospho results navset
  - New packages: `KSEAapp` (CRAN), `ggseqlogo` (CRAN)

- **Phosphoproteomics Phase 3 — Advanced Features**:
  - **Protein-level abundance correction**: Checkbox to subtract protein-level logFC from phosphosite logFC, isolating phosphorylation stoichiometry changes
  - **AI context integration**: Phosphosite DE results and KSEA kinase activities appended to Data Chat context when phospho analysis is active
  - **Session persistence**: KSEA results, FASTA sequences, and all Phase 2/3 state saved/loaded with sessions

### Changed
- `R/server_gsea.R`: Rewritten from 144 to ~400 lines (multi-DB, caching, organism detection)
- `R/server_ai.R`: AI Summary rewritten for all-contrast analysis; send_chat scaled for large datasets
- `R/helpers_phospho.R`: Extended from 210 to ~380 lines (5 new helper functions)
- `R/server_phospho.R`: Extended from 650 to ~950 lines (KSEA, motifs, protein correction)
- `Dockerfile`: Added `KSEAapp` and `ggseqlogo` CRAN package installation
- App version bumped to v2.5

### Fixed
- **Export CSV crash**: "Column name `Protein.Group` must not be duplicated" when limpa's topTable includes Protein.Group column
- **Gemini token limit**: Scale data sent to AI based on sample count; group-level Mean/SD for >100 samples
- **P-value histogram y-axis**: Cap y-axis to show distribution shape when first bin dominates; annotate clipped bin count
- **P-value dropdown clipping**: Removed card_body wrapper and added z-index stacking to prevent plot from overlapping dropdown
- **Comparison dropdown width**: Full-width dropdowns on all comparison banners; buttons moved inline

### Planned — Future
- PhosR integration (RUVphospho normalization, kinase profiling, signalome)
- FASTA-based protein-relative position mapping for Path B sites

## [2.4.0] - 2026-02-17

### Added
- **Phosphoproteomics Tab (Phase 1)**: Site-level differential phosphorylation analysis
  - **Auto-detection**: Scans `Modified.Sequence` for `UniMod:21` on file upload; shows blue banner in Data Overview with "Open Phospho Tab" button
  - **Two input paths**:
    - Path A (recommended): Upload DIA-NN 1.9+ `site_matrix_0.9.parquet` or `site_matrix_0.99.parquet`
    - Path B: Parse phosphosites directly from `report.parquet` with configurable localization confidence threshold (0.5–1.0)
  - **Site extraction algorithm** (Path B): Character-by-character Modified.Sequence parser to locate `(UniMod:21)` positions, expand multiply-phosphorylated peptides, aggregate per SiteID × Run via max intensity ("Top 1" method)
  - **Site-level limma DE**: Filter sites (≥2 non-NA per group), tail-based imputation (Perseus-style: mean − 1.8 SD, width 0.3 SD), optional normalization (none/median/quantile), standard limma pipeline
  - **Phospho Volcano**: ggplot2 volcano with ggrepel labels formatted as "Gene Residue+Position" (e.g., "MAPK1 T185"). Colors: Significant=#E63946, FDR-only=#457B9D, NS=gray70. Downloadable as PDF.
  - **Site Table**: DT datatable with SiteID, Gene, Residue, Position, logFC, adj.P.Val, localization confidence. Filterable, sortable, downloadable as CSV.
  - **Residue Distribution**: Grouped bar chart (S/T/Y) comparing "All quantified" vs "Significant". Subtitle with expected ~85% Ser / ~14% Thr / ~1% Tyr.
  - **QC: Completeness**: Histogram of per-site % samples quantified with red dashed line at 50% threshold.
  - **Sidebar controls**: Conditional on phospho detection — input mode, localization slider, normalization radio, "Run Phosphosite Analysis" button
  - **Normalization warning**: Yellow alert for phospho-enriched data explaining DIA-NN normalization assumptions
  - **Educational expandable**: Explains site-level vs protein-level analysis, localization confidence, imputation approach
  - **Session save/load**: All phospho state persisted and restored
  - **Reproducibility logging**: Pipeline steps logged with parameters
- New files: `R/helpers_phospho.R` (210 lines), `R/server_phospho.R` (650 lines)

## [2.3.0] - 2026-02-17

### Changed
- **Modularization**: Split 5,139-line monolith into `app.R` orchestrator + 12 `R/` module files
- **Upload limit**: Increased from 500 MB to 5 GB
- **Dockerfile**: Updated for directory-based `runApp()` and `COPY R/` directive

## [2.2.0] - 2026-02-17

### Added
- **Contextual Help System**: 15 info modal buttons (`?`) across all major tabs providing in-context guidance
  - QC Plots: Normalization Diagnostic, DPC Fit, MDS Plot, Group Distribution, P-value Distribution
  - Data Overview: Signal Distribution, Expression Grid
  - DE Dashboard: Volcano/table interaction guide with threshold explanation
  - Consistent DE: High-Consistency Table (%CV explained), CV Distribution
  - QC Trends: Metric definitions, sort order explanation, drift detection tips
  - Gene Set Enrichment: GSEA overview + Results Table column definitions (NES, p.adjust, etc.)
  - Reproducibility: Methodology (LIMPA pipeline, limma/eBayes, covariates)
  - Data Chat: Privacy info, API key instructions, selection integration
- **Volcano Plot → Results Table Filtering**: Selecting proteins in the volcano plot now filters the results table to show only selected proteins (bidirectional sync)
- **MDS Plot Legend**: Visible legend at bottom-right with white background box (was being clipped off-screen)
- **Heatmap Expanded by Default**: DE Dashboard heatmap accordion now opens expanded

### Changed
- **Normalization Diagnostic**: "What am I looking at?" moved from in-page expandable to modal dialog (no longer interferes with plot)
- **P-value Distribution**: "How do I interpret this?" moved from in-page expandable to modal dialog; guidance banner moved below plot (no longer overlaps comparison dropdown)

### Fixed
- **Bad Merge Recovery**: Restored 333 lines lost in merge commit `1361b62` including MS2 Intensity Alignment and XIC auto-load features
- **XIC Viewer Facet Error**: Added guard for "Faceting variables must have at least one value" warning when modal first opens
- **DE Table Row Index Mapping**: Fixed row selection observer to correctly index into filtered data when volcano selection is active

## [2.1.1] - 2026-02-16

### Added
- **XIC Chromatogram Viewer**: On-demand fragment-level chromatogram inspection for differentially expressed proteins
  - Sidebar section "5. XIC Viewer" with directory path input and load button
  - XIC buttons on DE Dashboard results table and Grid View modal
  - Three display modes: Facet by sample, Facet by fragment, Intensity alignment
  - Split-axis MS1/MS2 view with independent y-axes (MS1 top, fragments bottom)
- **MS2 Intensity Alignment**: Spectronaut-style stacked bar chart for fragment ion ratio consistency
  - Each bar = one sample, colored segments = relative fragment proportions
  - Automatic inconsistency detection: flags samples with deviation > mean + 2×SD
  - Green/amber guidance banners with sample IDs and possible causes
  - Bars ordered by experimental group with dashed separators
  - Cosine similarity and deviation scores in tooltips
  - Precursor selector, group filter, and MS1 toggle controls
  - Prev/Next protein navigation through significant DE proteins
  - Download handler for PNG export (14×10 inch, 150 DPI)
  - Info panel with protein stats, RT range, and DE statistics
- **DIA-NN 2.x Format Support**: Auto-detects and handles both DIA-NN 1.x (wide) and 2.x (long) XIC formats
- **Ion Mobility / Mobilogram Support**: Detects mobilogram files, checks for non-zero data (timsTOF/PASEF only)
  - Blue gradient toggle with bolt icon; prominent banner when IM mode is active
- **XIC Directory Auto-Population**: Auto-detects `_xic` sibling directory when data is loaded
  - Smart path resolution: accepts `.parquet` file paths or directories without `_xic` suffix
- **Precursor Map from In-Memory Data**: Builds protein→precursor mapping from `values$raw_data` (no file I/O)

- **XIC Auto-Load**: XICs automatically load when `_xic` directory is detected on data upload (no manual button click needed)

### Fixed
- **Assign Groups Layout**: Fixed Run Pipeline button pushed off-screen on MacBook (CSS Grid → Flexbox)
- **Arrow/dplyr Conflicts**: `arrow::select` masking `dplyr::select` — use explicit `dplyr::select()` in XIC code
- **Tidy Evaluation Issues**: `rlang::sym("pr")` and `rename(Precursor.Id = pr)` replaced with base R equivalents
- **Mobilogram Detection Pattern**: Fixed file pattern to match DIA-NN naming convention (`_mobilogram.parquet` not `.mobilogram.parquet`)
- **Windows Path Compatibility**: Fixed regex in XIC smart path resolution to use `basename()` instead of forward-slash pattern
- **HF XIC Visibility**: XIC Viewer sidebar, buttons, and auto-load logic hidden on Hugging Face Spaces (detected via `SPACE_ID` env var); replaced with info note linking to GitHub for local download

## [2.1.0] - 2026-02-13

### Added
- **Four-Way Comparison Selector Synchronization**: Signal Distribution, Expression Grid, P-value Distribution, and DE Dashboard selectors now sync automatically when any one changes
- **P-value Distribution Diagnostic Plot**: New QC Plots sub-tab with automated pattern detection
  - Color-coded guidance banners (healthy, inflation, low power, model issues)
  - Expandable interpretation guide
  - Actionable recommendations for troubleshooting
- **CV Distribution Histogram**: New Consistent DE sub-tab showing protein variability distribution per experimental group
- **Dataset Summary DE Counts**: Shows differential expression protein counts for all comparisons with explicit directional language ("X proteins higher in GroupA")
- **AI Summary Sub-Tab**: Dedicated tab in Data Overview (moved from DE Dashboard button)
  - Inline display instead of modal popup
  - Cleaner workflow for AI-generated summaries
- **Volcano Plot Annotations**:
  - Colored threshold lines (blue FDR, orange logFC)
  - Significance criteria legend box in upper-left corner
  - Uses plotly native annotations for reliable rendering
- **Comparison Selectors**: Purple gradient banners on Signal Distribution, Expression Grid, and P-value Distribution tabs

### Changed
- **Signal Distribution Plot**: Now always shows DE coloring when results are available (removed manual toggle buttons)
- **Signal Distribution Plot**: Uses dedicated `contrast_selector_signal` instead of main selector
- **Responsive Plot Heights**: All plots now use viewport-relative units (`vh`, `calc()`) for optimal viewing on any screen size
- **Assign Groups Interface**: Converted from modal popup to permanent sub-tab in Data Overview
- **DE Dashboard**: Removed AI Summary button (moved to Data Overview tab)

### Fixed
- **Volcano Plot Legend**: Fixed text rendering issues by switching from ggplot annotations to plotly native system
- **Reactive Loops**: Signal Distribution table no longer filters based on selection-derived reactives

### Technical
- All comparison selectors (`contrast_selector`, `contrast_selector_signal`, `contrast_selector_grid`, `contrast_selector_pvalue`) sync bidirectionally
- Plotly annotations use paper coordinates (`xref = "paper", yref = "paper"`) for absolute positioning
- No new R package dependencies (all features use existing packages)

## [2.0.0] - 2025-XX-XX

### Initial Release
- Interactive differential expression analysis for DIA-NN proteomics data
- Limpa pipeline integration (DPC normalization, maxLFQ quantification, limma statistics)
- Google Gemini AI chat integration
- Session save/load functionality
- GSEA integration (GO Biological Process)
- Reproducibility code logging
- Example data with Affinisep vs Evosep comparison
- Education tab with embedded resources
