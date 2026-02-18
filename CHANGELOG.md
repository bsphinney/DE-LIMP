# Changelog

All notable changes to DE-LIMP will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.5.0] - 2026-02-18

### Added
- **Phosphoproteomics Phase 2 — Kinase Activity & Motif Analysis**:
  - **KSEA kinase-substrate enrichment analysis**: Infers upstream kinase activity from phosphosite fold-changes using `KSEAapp` CRAN package with PhosphoSitePlus + NetworKIN database
  - **KSEA bar plot**: Horizontal bar chart of kinase z-scores (top 15 activated + top 15 inhibited), colored by direction, with substrate count annotations
  - **KSEA results table**: Filterable, sortable DT datatable with kinase gene, z-score, FDR, substrate count. Downloadable as CSV.
  - **Sequence logo motif analysis**: Displays amino acid enrichment around regulated phosphosites using `ggseqlogo`. Separate logos for up-regulated and down-regulated sites.
  - **FASTA upload**: Sidebar file input for protein FASTA. Parses UniProt-format headers to extract accessions. Enables accurate flanking sequence extraction for motif analysis.
  - New "Kinase Activity" and "Motif Analysis" tabs in phospho results navset
  - New packages: `KSEAapp` (CRAN), `ggseqlogo` (CRAN)

- **Phosphoproteomics Phase 3 — Advanced Features**:
  - **Protein-level abundance correction**: Checkbox to subtract protein-level logFC from phosphosite logFC, isolating phosphorylation stoichiometry changes. Requires total proteome pipeline to have been run first. Volcano plot title and axis labels update to indicate corrected values.
  - **AI context integration**: Phosphosite DE results (top 20 sites) and KSEA kinase activities (significant kinases) are automatically appended to Data Chat context when phospho analysis is active. Both manual chat and auto-analyze include phospho context.
  - **Session persistence**: KSEA results, FASTA sequences, and all Phase 2/3 state saved/loaded with sessions.
  - New helper functions: `read_fasta_sequences()`, `prepare_ksea_input()`, `extract_flanking_sequences()`, `correct_phospho_for_protein()`, `phospho_ai_context()`

### Changed
- `R/helpers_phospho.R`: Extended from 210 to ~380 lines (5 new helper functions)
- `R/server_phospho.R`: Extended from 650 to ~950 lines (KSEA, motifs, protein correction)
- `Dockerfile`: Added `KSEAapp` and `ggseqlogo` CRAN package installation
- App version bumped to v2.5

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
