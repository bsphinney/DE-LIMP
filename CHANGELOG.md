# Changelog

All notable changes to DE-LIMP will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
