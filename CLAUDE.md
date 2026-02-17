# DE-LIMP Project Context for Claude

## Working Preferences
- **ALWAYS update this file** after implementing changes to the codebase
- Update the "Recent Changes" section with detailed entries including line numbers and functionality
- Update other relevant sections (Main Tab Structure, Key Patterns, TODO list, etc.) as needed
- **When this file exceeds ~400 lines**, proactively recommend using the `claude-md-management:claude-md-improver` skill to audit and refactor
- Proactively suggest updates to this file when new patterns, gotchas, or architectural decisions emerge

## Project Overview
DE-LIMP is a Shiny proteomics data analysis pipeline using the LIMPA R package for differential expression analysis of DIA-NN data.

- **GitHub**: https://github.com/bsphinney/DE-LIMP
- **Hugging Face**: https://huggingface.co/spaces/brettsp/de-limp-proteomics
- **Local URL**: http://localhost:3838

## System Requirements
- **R**: 4.5+ (required for limpa)
- **Bioconductor**: 3.22+ (automatic with R 4.5+)
- **Download R**: https://cloud.r-project.org/

## Key Files

| File | Purpose |
|------|---------|
| `app.R` | Main Shiny app (~2800 lines) - used by HF Spaces |
| `DE-LIMP.R` | Identical copy of app.R - for GitHub releases/local users |
| `Dockerfile` | Docker container for HF Spaces and HPC deployment |
| `UI_LAYOUT_SPEC.md` | Spec for v2.1 responsive UI layout refactoring |
| `NORMALIZATION_DIAGNOSTIC_SPEC.md` | Spec for normalization diagnostic plot feature |
| `HPC_DEPLOYMENT.md` | Guide for HPC cluster deployment with Apptainer/Singularity |
| `USER_GUIDE.md` | End-user documentation |
| `README_GITHUB.md` | **SOURCE** for GitHub README (edit this!) |
| `README_HF.md` | **SOURCE** for HF README (edit this!) |
| `README.md` | Generated - content differs between remotes |
| `Citations/` | Reference PDFs (LIMMA, FDR, DIA-NN papers) |
| `docs/index.html` | GitHub Pages site |
| `.github/workflows/sync-to-hf.yml` | Auto-sync GitHub to HF |

## App Architecture

### Structure (app.R)
```
Lines 1-110:    Auto-installation & setup (R/Bioc version checks, package install)
Lines 111-168:  Server configuration (max upload size, Gemini API)
Lines 169-330:  Helper functions (QC stats, AI chat, organism detection)
Lines 331-850:  UI definition (page_sidebar with 9 main tabs)
Lines 851-4185: Server logic (~80+ reactive elements, incl. XIC viewer)
Line 4186:      shinyApp(ui, server)
```

### Main Tab Structure
- **Data Overview** - 6 sub-tabs: **Assign Groups & Run**, Signal Distribution, Dataset Summary, Group QC Summary, Expression Grid, **AI Summary**
- **QC Trends** - 4 sub-tabs: Precursors, Proteins, MS1 Signal, Stats Table
- **QC Plots** - 5 sub-tabs: Normalization Diagnostic, DPC Fit, MDS Plot, Group Distribution, P-value Distribution
- **DE Dashboard** - Volcano plot, results table, violin plots, interactive comparison selector
- **Consistent DE** - 2 sub-tabs: High-Consistency Table (ranked by %CV), CV Distribution (histogram by group)
- **Reproducibility** - Code Log + Methodology sub-tabs
- **Gene Set Enrichment** - Dot Plot, Enrichment Map, Ridgeplot, Results Table
- **Data Chat** - AI-powered analysis (Google Gemini API)
- **Education** - Learning resources

### Key Reactive Values
- `values$plot_selected_proteins` - Selected proteins from table/volcano (used by all viz)
- `values$fit` - limma fit object from DE analysis
- `values$y_protein` - Protein-level quantification matrix
- `values$repro_log` - Cumulative R code for reproducibility
- `values$xic_dir` - Path to XIC parquet directory
- `values$xic_available` - Whether XIC files were detected
- `values$xic_format` - "v1" (DIA-NN 1.x) or "v2" (DIA-NN 2.x)
- `values$xic_protein` - Currently selected protein for XIC viewing
- `values$xic_data` - Loaded & reshaped XIC data for current protein
- `values$xic_report_map` - Protein → Precursor mapping from report
- `values$uploaded_report_path` - Path to uploaded report.parquet
- `values$original_report_name` - Original filename of uploaded report (for XIC path derivation)
- `values$mobilogram_available` - Whether mobilogram files with non-zero IM data exist

### LIMPA Pipeline Flow
1. `readDIANN()` - Load DIA-NN parquet file
2. `dpcCN()` - Data Point Correspondence normalization
3. `dpcQuant()` - Peptide to protein quantification
4. `dpcDE()` - Differential expression model
5. `contrasts.fit()` + `eBayes()` - Pairwise comparisons

### Selection System
- Volcano plot click/box-select and table row selection both update `values$plot_selected_proteins`
- Bidirectional sync for highlighting

### Comparison Selector Synchronization
- **Four synchronized selectors**: `contrast_selector` (DE Dashboard), `contrast_selector_signal` (Signal Distribution), `contrast_selector_grid` (Expression Grid), `contrast_selector_pvalue` (P-value Distribution)
- **Bidirectional sync**: Changing any selector updates all others automatically
- **Population**: All four selectors populated when pipeline runs and when session loads
- **Independence**: Each plot uses its own dedicated selector (not the main one) to avoid reactive conflicts

### AI Chat Feature
- Uses Google Gemini API, uploads top 800 proteins via File API
- Bidirectional: user selects proteins for AI analysis, AI suggests proteins to highlight
- Selection format: `[[SELECT: P12345; P67890]]`

## Development Workflow

### Running Locally
```r
# In VS Code R terminal (recommended):
shiny::runApp('/Users/brettphinney/Documents/claude/DE-LIMP.r', port=3838, launch.browser=TRUE)

# From command line:
Rscript -e "shiny::runApp('/Users/brettphinney/Documents/claude/DE-LIMP.r', port=3838)" &
```

- **DO NOT** use `source()` - it doesn't work properly in VS Code
- Shiny apps don't hot-reload - must restart after every code change
- Stop: `pkill -f "DE-LIMP.r"` or `Ctrl+C`
- Check if running: `lsof -i :3838`

### Keeping app.R and DE-LIMP.R in Sync
These files must be identical. After editing one: `cp DE-LIMP.R app.R` (or vice versa).

## Deployment

### Three Platforms
1. **GitHub** (`origin` remote) - Source code, releases, download DE-LIMP.R
2. **Hugging Face Spaces** (`hf` remote) - Web app via Docker (app.R)
3. **HPC Clusters** - Apptainer/Singularity containers (see HPC_DEPLOYMENT.md)

### GitHub Actions Auto-Sync (Primary Workflow)
Pushing to `origin/main` automatically syncs to HF within 1-2 minutes.
- Workflow: `.github/workflows/sync-to-hf.yml`
- Syncs all files **except** README.md and .github/
- Uses `README_HF.md` as HF's README (with YAML frontmatter)
- **You only need `git push origin main`** - HF updates automatically

```bash
# Standard workflow for code changes:
git add DE-LIMP.R app.R  # be specific with files
git commit -m "Description"
git push origin main      # HF syncs automatically
```

### README Management (CRITICAL)
The two remotes have **permanently different README.md files**. GitHub Actions handles this automatically, but if you ever need to manually update READMEs:

- Edit `README_GITHUB.md` for GitHub, push to `origin` only
- Edit `README_HF.md` for HF, push to `hf` only
- **NEVER** push README.md changes to both remotes
- **NEVER** use `git add .` when README.md is modified

**Recovery if HF breaks** ("Missing configuration in README" error):
```bash
cp README.md README_GITHUB_BACKUP.md
cp README_HF.md README.md
git add README.md && git commit -m "Fix HF: Restore YAML README" && git push hf main
cp README_GITHUB_BACKUP.md README.md
git add README.md && git commit -m "Restore GitHub README" && git push origin main
rm README_GITHUB_BACKUP.md
```

### Adding New R Packages
When adding features that require new packages:
1. Add `library()` call to app.R/DE-LIMP.R
2. Update Dockerfile to install the package (CRAN in Section 2, Bioconductor in Section 3)
3. Consider system dependencies: graphics packages need `libcairo2-dev`, XML/web need `libxml2-dev`
4. Commit app.R, DE-LIMP.R, and Dockerfile together

### Minimize HF Builds
HF Docker builds take 5-10 min (cached) or 30-45 min (Dockerfile changes). Always test locally first, batch changes, and push only when fully tested.

## Key Patterns & Gotchas

### R Shiny Reactivity
- **DT table row indices** refer to the CURRENTLY DISPLAYED (possibly filtered) data
- **Avoid circular reactivity**: Don't filter data based on a reactive value that the filter updates
- **CRITICAL reactive loop**: If `renderDT` depends on a reactive that uses `values$selection`, AND selecting rows updates `values$selection`, you get an infinite loop (selection triggers re-render which resets selection). **Solution**: Make `renderDT` build data independently, not using reactives that depend on selection state.
- Use `isolate()` to break reactive dependencies when needed
- DT tables with `selection = "multiple"` support Ctrl+Click and Shift+Click

### Package Installation & Loading
- Install packages BEFORE any `library()` calls (attempting to update loaded packages causes unload errors)
- Use `update = FALSE` in `BiocManager::install()` and `upgrade = "never"` in `remotes::install_*()`
- Set CRAN mirror before installs: `options(repos = c(CRAN = "https://cloud.r-project.org"))`

### dplyr/tidyverse
- In `summarise()`, use explicit `{}` braces for multi-line if statements
- Always set `.groups = 'drop'` to avoid grouping warnings

### Reproducibility Logging
- Uses cumulative `add_to_log(action_name, code_lines)` helper function
- Logs with timestamps, appends rather than overwrites
- Include context: which button was clicked, which parameters changed

### UI Design Patterns (v2.1)
- All plot heights use viewport-relative units (`vh`, `calc()`) - no fixed pixel heights
- Multi-plot sections use `navset_card_tab` sub-tabs instead of dropdowns or modals
- DE Dashboard uses responsive grid (`.de-dashboard-grid`) that stacks at <1200px
- Fullscreen modals available as enhancement on all plot panels
- **CRITICAL bslib component issue**: `card()` and `card_body()` components don't render when placed at top level inside `nav_panel()`. Use plain `div()` elements with inline CSS instead for top-level content in nav panels. Cards work fine in the main body content, just not at the very top of a nav_panel.

## Common Issues

| Problem | Solution |
|---------|----------|
| `source()` doesn't start app in VS Code | Use `shiny::runApp()` instead |
| App doesn't reflect code changes | Restart app (no hot-reload) |
| "Please select a CRAN mirror" hangs | Already fixed: CRAN mirror set at line 8 |
| "limpa not found" / R version error | Need R 4.5+ from https://cloud.r-project.org/ |
| "Package cannot be unloaded" during install | Restart R session, packages install before loading |
| Selections disappear after clicking | Reactive loop - table must not depend on selection-derived reactives |
| bslib `card()` doesn't render in nav_panel | Use plain `div()` with inline CSS for top-level nav_panel content |
| HF "Missing configuration in README" | Run README recovery script above |

## Recent Changes (v2.1)

### 2026-02-16: XIC Bug Fixes & Auto-Population
1. **Fixed "object 'pr' not found" error** (load_xic_for_protein helper)
   - `rlang::sym("pr")` filtering didn't work reliably with Arrow/dplyr
   - Replaced with base R subsetting: `df[df$pr %in% target_prs, , drop = FALSE]`
   - Applied to both v2 (`pr` column) and v1 (`Precursor.Id` column)

2. **Fixed precursor map building** (XIC directory loader, lines 3910-3930)
   - Eliminated all report file I/O — builds map directly from `values$raw_data$E` rownames + `values$raw_data$genes$Protein.Group`
   - Previous approach failed due to: temp file deletion, `Run` vs `File.Name` column mismatch, `rlang` rename errors
   - Now just 3 lines of base R — zero file reading needed

3. **Fixed Assign Groups layout** (lines 465-508 CSS)
   - MacBook: Run Pipeline button pushed off-screen, no vertical scroll
   - Changed from rigid CSS Grid (`grid-template-columns: 200px 1fr 250px`) to Flexbox with `flex-wrap: wrap`
   - Controls now wrap responsively on smaller screens

4. **XIC Directory Auto-Population** (lines 1374-1393, 1306-1314, 3875-3900)
   - On data upload: auto-detects `<report_name>_xic/` directory in working directory
   - On example data load: checks for `Affinisep_vs_evosep_noNorm_xic/` in working directory
   - Pre-fills `xic_dir_input` text field with `updateTextInput()`
   - Shows notification when auto-detected
   - Stores `values$original_report_name` for future path derivation

5. **Smart XIC Path Resolution** (XIC directory loader, lines 3875-3900)
   - If user enters a `.parquet` file path → derives `_xic` sibling directory automatically
   - If user enters a path without `_xic` suffix → tries appending `_xic`
   - Updates the text input to show the resolved path

### 2026-02-16: DIA-NN 2.x Format Support & Mobilogram Detection
1. **Dual Format Support** (helper functions lines 1088-1250)
   - `detect_xic_format(xic_dir)`: Auto-detects DIA-NN version from column names
   - **v2** (DIA-NN 2.x): Long format with `pr`, `feature`, `info`, `rt`, `value` columns
   - **v1** (DIA-NN 1.x): Wide format with numbered columns + `Retention.Times`/`Intensities` flag rows
   - `load_xic_for_protein()`: Builds `StrippedSequence+Charge` key for v2, `Precursor.Id` for v1
   - Per-file reading preserves `Source.File` for sample identity (v2 has no File.Name column)
   - `reshape_xic_for_plotting()`: v2 data already long (just rename + join metadata), v1 does full pivot

2. **Mobilogram Detection** (XIC directory loader)
   - Scans for `.mobilogram.parquet` files alongside `.xic.parquet`
   - Samples one file to check if data is non-zero (only timsTOF/IM instruments produce real data)
   - `values$mobilogram_available` flag controls UI toggle visibility
   - Status badge shows "+ IM" when ion mobility data detected
   - "Show Ion Mobility" checkbox conditionally appears in XIC modal controls

3. **Precursor Mapping** (report map loading)
   - Now loads `Precursor.Charge` column from report for v2 key building
   - v2 key: `paste0(Stripped.Sequence, Precursor.Charge)` matches XIC `pr` column
   - Graceful fallback if `Precursor.Charge` column not present in older reports

### 2026-02-15: XIC Viewer for Differentially Expressed Proteins
1. **Sidebar Section "5. XIC Viewer"** (lines 438-447 UI)
   - `textInput` for XIC directory path (HPC/local users paste path)
   - "Load XICs" button with wave-square icon
   - Status badge shows file count, format version, and IM availability
   - Not visible on HF (XIC files too large for cloud deployment)

2. **XIC Trigger Buttons**
   - DE Dashboard results table: "📈 XICs" button (line 864) alongside Reset, Violin, Export
   - Grid View modal footer: "📈 XICs" button (line 1989) alongside Back to Grid

3. **Helper Functions** (lines 1088-1250)
   - `detect_xic_format()`: Auto-detect v1 vs v2 format
   - `load_xic_for_protein(xic_dir, protein_id, report_map, xic_format)`: Arrow predicate pushdown, dual format
   - `reshape_xic_for_plotting(xic_raw, metadata, xic_format)`: Format-aware reshaping

4. **XIC Modal** (lines 3950+)
   - Full-width modal with controls: Display mode, Precursor selector, Group filter, MS1 checkbox, IM toggle
   - Three display modes: Overlay (fragments per sample), Facet by fragment, Facet by sample
   - Plotly interactive chromatograms with tooltips (Sample, Fragment, RT, Intensity)
   - Info panel shows protein stats, RT range, DE stats (logFC, adj.P.Val)
   - Prev/Next protein navigation through significant DE proteins
   - Download handler exports PNG (14x10 inch, 150 DPI)

5. **Data Loading Architecture**
   - On-demand lazy loading: Arrow `open_dataset()` with predicate pushdown
   - Fallback: individual per-file reads if schema mismatch across DIA-NN versions
   - Protein → Precursor mapping loaded from report.parquet on XIC dir load
   - Caps at top 6 precursors for very large proteins (e.g., Titin)

6. **Report Path Storage** (lines 1289, 1234)
   - `values$uploaded_report_path` stored on both user upload and example data load
   - Enables XIC precursor mapping without re-reading the report

### 2026-02-13: Added Comparison Selectors to Data Overview Tab
1. **Signal Distribution sub-tab** (lines 502-523)
   - Purple gradient banner with comparison selector (matches DE Dashboard style)
   - Shows which comparison is being used for coloring
   - Directly selectable - changes sync across all tabs
   - Height adjusted to `calc(100vh - 420px)` to accommodate banner

2. **Expression Grid sub-tab** (lines 524-552)
   - Purple gradient banner with comparison selector
   - Makes it clear which DE results are shown in the grid
   - Selector syncs bidirectionally with main selector and Signal Distribution

3. **Synchronization Logic** (lines 1496-1532)
   - Three comparison selectors stay in sync: main (DE Dashboard), Signal Distribution, Expression Grid
   - Change any selector → all others update automatically
   - Prevents confusion about which comparison is being viewed
   - All three populated when pipeline runs or session loads

### 2026-02-13: Added P-value Distribution Diagnostic to QC Plots
1. **New QC Plots Sub-Tab: P-value Distribution** (lines 704-708 UI, 2180-2333 server)
   - Histogram of raw p-values from differential expression analysis
   - Helps diagnose statistical issues: p-value inflation, lack of power, model problems
   - Red dashed line shows expected uniform distribution under the null hypothesis
   - Fullscreen button for detailed examination

2. **Automated Pattern Detection & Contextual Guidance** (NEW):
   - **Real-time health assessment**: Analyzes p-value distribution automatically
   - **Color-coded guidance banners**: Appear above plot based on detected pattern
     - 🟢 **Green (Healthy)**: Good spike near p=0 with uniform background
     - 🟡 **Yellow (Warning)**: P-value inflation or low statistical power detected
     - 🔴 **Red (Danger)**: U-shaped distribution indicating model issues
     - 🔵 **Blue (Info)**: Completely uniform (no DE signal detected)
   - **Actionable recommendations**: Each banner explains what was detected and what to do
   - **Detection metrics**:
     - Spike detection: Low p-value ratio vs expected (>2x = spike)
     - Inflation detection: Mid-range (0.3-0.7) excess (>1.3x expected)
     - Power assessment: Depletion of small p-values (<0.5x expected)
     - U-shape detection: Both ends elevated (first & last bins >1.5x expected)

3. **What it Detects**:
   - **Healthy pattern**: Mostly flat (uniform) with spike near p=0 (true positives)
   - **P-value inflation**: Too many intermediate values (0.3-0.7) suggests unmodeled variance, batch effects, or small sample size
   - **Low power**: Depletion near zero means test is too conservative or underpowered
   - **Model issues**: U-shaped distribution indicates statistical model problems or normalization issues
   - **No signal**: Completely uniform = no differential expression detected

4. **Built-in Guidance**:
   - **Dynamic banners** with specific suggestions based on detected pattern
   - Expandable "How do I interpret this?" section explains patterns
   - Links to other diagnostic plots (Normalization, MDS, CV Distribution) for troubleshooting
   - Suggests fixes: add covariates, check normalization, verify sample sizes, increase replicates

5. **Technical Details**:
   - Uses raw P.Value (not adj.P.Val) to assess null hypothesis behavior
   - 30 bins (regular view), 40 bins (fullscreen)
   - Comparison-specific (updates when contrast changes)
   - Height: `calc(100vh - 400px)` for responsive sizing
   - Pattern detection runs automatically on each comparison change

### 2026-02-13: Added CV Histogram to Consistent DE Tab
1. **Converted Consistent DE to Sub-Tabs** (lines 769-791 UI, 2296-2445 server)
   - Tab 1: "High-Consistency Table" - existing table ranked by %CV
   - Tab 2: "CV Distribution" - new histogram visualization
   - Uses `navset_card_tab` pattern for easy switching between views

2. **CV Distribution Histogram Features**:
   - Histogram of CV values for all significant proteins (adj.P.Val < 0.05)
   - Faceted by experimental group (separate panel for each group)
   - Dashed vertical line shows average CV for each group
   - Color-coded by group with annotations showing exact average values
   - Height: `calc(100vh - 320px)` for responsive viewport sizing
   - Fullscreen button opens extra-large modal view (700px height, 1000px width)

3. **Implementation Details**:
   - CV calculated in linear space (2^log2 expression) for each group
   - CV = (SD / Mean) × 100
   - Uses `pivot_longer()` to convert wide CV data to long format for ggplot2
   - Histogram bins = 30 (regular view), 40 (fullscreen)
   - `scales = "free_y"` allows each group's histogram to optimize its own Y-axis

4. **Rationale**:
   - Complements the existing table by showing the distribution shape
   - Helps identify outlier proteins with unusually high/low CV
   - Group-level averages show which conditions have more/less variability
   - Lower CV = more stable and reproducible biomarker candidates

### 2026-02-13: Interactive Comparison Selector on DE Dashboard
1. **Moved comparison selector from sidebar to DE Dashboard banner** (lines 654-672)
   - Used plain `div()` with inline CSS (gradient purple background)
   - Dropdown selector embedded directly in banner for contextual access
   - Removes redundancy (was duplicated in sidebar)

### 2026-02-13: Assign Groups Moved to Data Overview Sub-Tab
1. **Converted modal popup to permanent sub-tab** (lines 446-500)
   - Fixes covariate text box layout issues (screenshot from user)
   - Full screen width for metadata table (`calc(100vh - 480px)` height)
   - Proper grid layout for controls: `grid-template-columns: 200px 1fr 250px`
   - Removed modal observer code (~47 lines)
   - Navigation: `nav_select("data_overview_tabs", "Assign Groups & Run")` instead of `click("open_setup")`

### 2026-02-13: Signal Distribution Always Shows DE Coloring
1. **Removed color toggle buttons** (lines 517-521)
   - Deleted "Color by DE" and "Reset" buttons from UI
   - Removed `values$color_plot_by_de` reactive flag
   - Plot automatically colors by DE status when results are available

2. **Updated plot to use dedicated selector** (lines 1757-1784 main, 2158-2180 fullscreen)
   - Changed from `input$contrast_selector` to `input$contrast_selector_signal`
   - Plot now responds to its own comparison selector in the purple banner
   - Always shows DE coloring when fit object exists (no manual toggle needed)
   - Height adjusted to `calc(100vh - 370px)` (removed button row saves 50px)

3. **Removed button observers** (previously lines 1758-1759)
   - Deleted `observeEvent(input$color_de)` and `observeEvent(input$reset_color)`
   - Simplified reactive logic

### 2026-02-13: Four-Way Comparison Selector Synchronization
1. **Added P-value Distribution selector** (lines 729-740 UI)
   - Purple gradient banner with `contrast_selector_pvalue` selector
   - Matches style of Signal Distribution and Expression Grid selectors
   - Plot height adjusted to `calc(100vh - 450px)` to accommodate banner

2. **Updated P-value Distribution to use dedicated selector** (lines 2279-2282, 2435-2495)
   - `assess_pvalue_health()` reactive uses `input$contrast_selector_pvalue`
   - Main histogram plot uses `input$contrast_selector_pvalue`
   - Fullscreen modal uses `input$contrast_selector_pvalue`

3. **Four-way synchronization** (lines 1513-1543)
   - Main selector (DE Dashboard) syncs to: Signal Distribution, Expression Grid, P-value Distribution
   - Each of the 3 sub-selectors syncs back to main selector
   - Change any selector → all 4 update automatically
   - Populated when pipeline runs (lines 1475-1479) and session loads (lines 3028-3032)

### 2026-02-13: Dataset Summary Shows DE Protein Counts
1. **Added Differential Expression Summary section** (lines 1680-1743)
   - Shows DE protein counts for all comparisons
   - Displays total significant proteins (adj.P.Val < 0.05)
   - **Explicit language**: "X proteins higher in [Group A]" and "Y proteins higher in [Group B]"
   - Uses arrows (↑/↓) to indicate direction of change
   - Color-coded: red for Group A (positive logFC), blue for Group B (negative logFC)
   - Bordered cards with purple left border for each comparison
   - Parses comparison names to extract group labels
   - Auto-updates when pipeline runs
   - Example display:
     ```
     🔬 Evosep - Astral (245 significant proteins)
        ↑ 189 proteins higher in Evosep
        ↓ 56 proteins higher in Astral
     ```

### 2026-02-13: AI Summary Moved to Data Overview Tab
1. **New AI Summary sub-tab in Data Overview** (lines 564-600 UI, 1716-1802 server)
   - 6th sub-tab in Data Overview (after Expression Grid)
   - Purple gradient header matching other tabs
   - Large "🤖 Generate AI Summary" button
   - Output displays inline (no modal popup)
   - Uses new button ID: `generate_ai_summary_overview`

2. **Removed from DE Dashboard** (lines 849-857)
   - Deleted "🤖 AI Summary" button from Results Table header
   - Removed old observer `observeEvent(input$generate_ai_summary)` (~31 lines)
   - Cleaner DE Dashboard with just Reset, Violin, and Export buttons

3. **Server logic features**:
   - Gathers top 50 significant proteins with gene symbols
   - Identifies 3 most stable DE proteins (lowest CV across groups)
   - Generates 3-paragraph summary via Gemini API
   - Displays with markdown formatting in bordered output area

### 2026-02-13: Volcano Plot FDR Threshold Annotations & Legend
1. **Added colored threshold lines** (lines 3550-3579 main, 2744-2773 fullscreen)
   - Horizontal FDR line (blue): p-value = 0.05
   - Vertical logFC lines (orange): ±fold-change cutoff (user-adjustable slider)
   - Lines colored instead of plain gray dashed for better visibility

2. **Significance criteria legend box**:
   - Uses **plotly's native `layout()` annotations** (not ggplot - better rendering)
   - White semi-transparent box in upper-left corner
   - Header: "Significant if:"
   - Bullet points showing current thresholds:
     - "FDR-adj. p < 0.05"
     - "|log2FC| > X.XX" (dynamically updates with slider)
   - Positioned using paper coordinates (% from edges) - never overlaps with data
   - HTML formatting supported (`<b>` tags, `<br>` for line breaks)

3. **Technical implementation**:
   - Removed ggplot `annotate()` calls (unreliable in plotly conversion)
   - Added plotly `layout(annotations = list(...), shapes = list(...))` after `ggplotly()`
   - Annotations use `xref = "paper", yref = "paper"` for absolute positioning
   - Legend box: `x0 = 0.01, x1 = 0.28, y0 = 0.88, y1 = 0.99` (paper coords)
   - Text positioned at `x = 0.02, y = 0.98` and `y = 0.93` (paper coords)

4. **Enhanced visual clarity**:
   - Plot title shows current comparison
   - Color scheme: Blue = FDR threshold, Orange = FC thresholds
   - Legend box never obscures data points or threshold lines
   - Both main plot and fullscreen modal have identical annotations
   - Consistent rendering across zoom/pan operations

## Current TODO
- [x] Volcano plot: Add FDR threshold annotation/legend ✅
- [ ] GSEA: Add KEGG/Reactome enrichment; clarify which contrast is used
- [x] DE Dashboard: Make comparison bar clickable to change contrasts ✅
- [ ] Grid View: Open violin plot on protein click with bar plot toggle
- [x] Consistent DE: Add CV histogram by condition ✅
- [x] Data Overview: Show which comparison is used for signal distribution ✅
- [x] Data Overview: Show DE protein counts per comparison ✅
- [x] Move AI Summary to Data Overview tab ✅
- [ ] Publication-quality plot exports (SVG/PNG/TIFF with size controls)
- [ ] Sample correlation heatmap (QC Plots tab)
- [ ] Venn diagram of significant proteins across comparisons
- [x] P-value histogram diagnostic ✅
- [x] P-value Distribution: Add comparison selector ✅
- [ ] Sample CV distribution plots
- [ ] Protein numbers bar plot per sample
- [ ] Absence/presence table for on/off proteins
- [x] XIC Viewer: On-demand chromatogram viewer for DE proteins ✅
