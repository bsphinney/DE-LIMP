# DE-LIMP Project Context for Claude

## Working Preferences
- **Update this file** when new patterns, gotchas, or architectural decisions emerge
- Update TODO list, Main Tab Structure, Key Patterns sections as needed
- For detailed change history, update `CHANGELOG.md` (not this file)
- **When this file exceeds ~300 lines**, proactively recommend using the `claude-md-management:claude-md-improver` skill to audit and refactor

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
| `app.R` | Main Shiny app (~4900 lines) - used by HF Spaces. **Must be identical to DE-LIMP.R** |
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

### Structure (app.R, ~4900 lines)
```
Lines 1-110:     Auto-installation & setup (R/Bioc version checks, package install)
Lines 111-170:   Server configuration (max upload size, HF detection, Gemini API)
Lines 170-330:   Helper functions (QC stats, AI chat, organism detection, XIC helpers)
Lines 331-900:   UI definition (page_sidebar with 9 main tabs, 15 info modal buttons)
Lines 900-4900:  Server logic (~80+ reactive elements, XIC viewer, info modal observers)
Line ~4900:      shinyApp(ui, server)
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
- `values$mobilogram_files_found` - Count of mobilogram files detected (may have zero data)

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

### Hugging Face Environment Detection
- `is_hf_space <- nzchar(Sys.getenv("SPACE_ID", ""))` — set once at startup (line 147)
- XIC sidebar section, XICs buttons (DE Dashboard + Grid View), and auto-load logic all guarded by `if (!is_hf_space)`
- On HF: info note with GitHub download link replaces the XIC sidebar section
- On Local/HPC: full XIC viewer with path input, auto-detect, auto-load

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

### UI Design Patterns
- All plot heights use viewport-relative units (`vh`, `calc()`) - no fixed pixel heights
- Multi-plot sections use `navset_card_tab` sub-tabs instead of dropdowns or modals
- DE Dashboard uses responsive grid (`.de-dashboard-grid`) that stacks at <1200px
- Fullscreen modals available as enhancement on all plot panels
- **CRITICAL bslib component issue**: `card()` and `card_body()` components don't render when placed at top level inside `nav_panel()`. Use plain `div()` elements with inline CSS instead for top-level content in nav panels.
- **Info modal pattern**: `actionButton("[id]_info_btn", icon("question-circle"), class="btn-outline-info btn-sm")` + `observeEvent(input$[id]_info_btn, { showModal(modalDialog(...)) })`. DE Dashboard uses `btn-outline-light` (white on purple banner).
- **Plotly annotations**: Use `layout(annotations = list(...))` with paper coordinates (`xref = "paper"`) for absolute positioning. Don't use ggplot `annotate()` — it renders unreliably in plotly conversion.
- **In-plot legends**: Use `legend("bottomright", bg = "white", box.col = "gray80")` inside the plot area. Never use negative `inset` values — they push legends off-screen.

### XIC Viewer Patterns
- **Arrow masks dplyr**: `arrow::select()` masks `dplyr::select()` — always use `dplyr::select()` explicitly in XIC code
- **Avoid rlang in Arrow context**: `rlang::sym()` and tidy `rename()` fail with Arrow — use base R `df[df$col %in% vals, ]` and `names(df)[...] <- "new"`
- **Precursor map**: Built from in-memory `values$raw_data` (rownames + `$genes$Protein.Group`), not file I/O
- **Dual format**: `detect_xic_format()` auto-detects v1 (wide) vs v2 (long) DIA-NN XIC formats
- **Auto-load**: Uses `shinyjs::delay(500, shinyjs::click("xic_load_dir"))` after auto-detecting `_xic` directory
- **Platform guard**: All XIC UI/logic wrapped in `if (!is_hf_space)` — XIC files too large for cloud

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
| `arrow::select` masks `dplyr::select` | Use `dplyr::select()` explicitly in XIC helper functions |
| `rename(New = old)` fails with "object not found" | Use base R `names(df)[names(df) == "old"] <- "New"` |
| XIC buttons visible on HF but don't work | Already fixed: hidden via `is_hf_space` flag |
| Mobilogram files not detected | Pattern was `.mobilogram.parquet`; DIA-NN uses `_mobilogram.parquet` — already fixed |
| Windows path separator in XIC regex | Already fixed: use `basename()` instead of literal `/` in regex |

## Version History

Current version: **v2.2.0** (2026-02-17). See [CHANGELOG.md](CHANGELOG.md) for detailed release history.

### Key Architecture Decisions (for context)
- **Volcano → Table filtering** (v2.2): `output$de_table` filters by `values$plot_selected_proteins`. Row selection observer uses `isolate()` to index into filtered data. `req(length > 0)` prevents reactive loop when selection resets to NULL.
- **Contextual Help** (v2.2): 15 info modal `?` buttons across all major tabs. See "Info modal pattern" in Key Patterns above.
- **XIC Viewer** (v2.1): Fragment-level chromatogram inspection with 3 display modes (facet by sample, facet by fragment, intensity alignment). Supports DIA-NN v1/v2 formats and ion mobility. Local/HPC only (hidden on HF).
- **MS2 Intensity Alignment** (v2.1): Spectronaut-style stacked bar chart. Inconsistency detection via `xic_alignment_data()` reactive (deviation > mean + 2*SD flags samples).
- **Four-way selector sync** (v2.1): `contrast_selector`, `contrast_selector_signal`, `contrast_selector_grid`, `contrast_selector_pvalue` all sync bidirectionally.
- **P-value Distribution Diagnostic** (v2.1): Automated health assessment with color-coded guidance banners (healthy/warning/danger/info patterns).
- **Assign Groups** (v2.1): Moved from modal to permanent Data Overview sub-tab for better layout.

## Current TODO
- [ ] GSEA: Add KEGG/Reactome enrichment; clarify which contrast is used
- [ ] Grid View: Open violin plot on protein click with bar plot toggle
- [ ] Publication-quality plot exports (SVG/PNG/TIFF with size controls)
- [ ] Sample correlation heatmap (QC Plots tab)
- [ ] Venn diagram of significant proteins across comparisons
- [ ] Sample CV distribution plots
- [ ] Protein numbers bar plot per sample
- [ ] Absence/presence table for on/off proteins
