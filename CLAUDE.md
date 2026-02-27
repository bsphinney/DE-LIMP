# DE-LIMP Project Context for Claude

## Working Preferences
- **Update this file** when new patterns, gotchas, or architectural decisions emerge
- For detailed change history, update `CHANGELOG.md` (not this file)
- Detailed patterns: @docs/PATTERNS.md | TODO list: @docs/TODO.md

## Project Overview
DE-LIMP is a Shiny proteomics data analysis pipeline using the LIMPA R package for differential expression analysis of DIA-NN data.

- **GitHub**: https://github.com/bsphinney/DE-LIMP
- **Hugging Face**: https://huggingface.co/spaces/brettsp/de-limp-proteomics
- **Local URL**: http://localhost:3838
- **R**: 4.5+ required (limpa needs Bioconductor 3.22+)

## Architecture

### Structure (~7,500 lines total)
```
app.R (270 lines):       Package loading, backend detection (Docker/HPC/Core Facility), reactive values, module calls
R/ui.R (~1,500 lines):   build_ui() — page_navbar layout, CSS/JS, accordion sidebar, all tab nav_panels
R/server_*.R (11 files): Server modules, each receives (input, output, session, values, ...)
R/helpers*.R (5 files):  Pure utility functions (no Shiny reactivity)
```

### Key Files

| File | Purpose |
|------|---------|
| `app.R` | Orchestrator — package loading, backend detection, reactive values, module calls |
| `R/ui.R` | `page_navbar` layout, accordion sidebar, all tab definitions |
| `R/server_data.R` | Data upload, example load, group assignment, pipeline execution |
| `R/server_de.R` | Volcano, DE table, heatmap, CV analysis, selection sync |
| `R/server_qc.R` | QC sample metrics (faceted trend plot), diagnostic plots, p-value distribution |
| `R/server_viz.R` | Expression grid, signal distribution, PCA |
| `R/server_gsea.R` | GSEA analysis, multi-DB (BP/MF/CC/KEGG), organism detection |
| `R/server_ai.R` | AI Summary (all contrasts), Data Chat, Gemini integration, HTML report export |
| `R/server_search.R` | Docker/HPC dual backend, SSH, DIA-NN search, job queue |
| `R/server_phospho.R` | Phospho site-level DE, volcano, site table |
| `R/server_mofa.R` | MOFA2 multi-view integration |
| `R/server_facility.R` | Core facility: reports, job history, QC dashboard |
| `R/server_session.R` | Info modals, save/load session, reproducibility |
| `R/helpers_search.R` | `ssh_exec()`, `build_diann_flags()`, `generate_sbatch_script()`, UniProt search |

### Tab Structure (page_navbar)
Navbar: **New Search** (conditional) | **QC** | **Analysis** dropdown | **Output** dropdown (Export Data, Methods & Code) | **Education** | **Facility** dropdown (conditional) | gear icon (far right)

- `page_navbar(id = "main_tabs", navbar_options = navbar_options(bg = "#2c3e50"))` — dark navbar, global sidebar, hover dropdowns
- Dropdown section labels ("Setup"/"Results"/"AI") injected via JS

**Progressive reveal**: `nav_hide()`/`nav_show()` on `"main_tabs"`. Hidden on startup via `session$onFlushed(once=TRUE)`:
- **Always visible**: New Search (if `search_enabled`), Analysis > Data Overview, Education, Facility (if `is_core_facility`)
- **QC**: shown when `values$raw_data` not NULL
- **DE Dashboard, GSEA, MOFA2, AI Analysis, Output**: shown when `values$fit` not NULL
- **Phosphoproteomics**: shown when `values$phospho_detected$detected` is TRUE

**Tab values that MUST NOT change** (used by server nav_select/nav_show/nav_hide):
`"QC"`, `"DE Dashboard"`, `"Gene Set Enrichment"`, `"mofa_tab"`, `"AI Analysis"`, `"Output"`, `"Phosphoproteomics"`, `"Data Overview"`, `"data_overview_tabs"`, `"Assign Groups & Run"`

#### Analysis dropdown
- **Data Overview** — `navset_card_tab(id = "data_overview_tabs")`: Assign Groups & Run, Signal Distribution, Dataset Summary, Replicate Consistency, Expression Grid, AI Summary
- **DE Dashboard** — `navset_card_tab(id = "de_dashboard_subtabs")`: Volcano (+heatmap), Results Table, PCA, CV Analysis. Comparison selector banner above sub-tabs.
- **Phosphoproteomics** — conditional on phospho detection
- **Gene Set Enrichment** — BP/MF/CC/KEGG with per-ontology caching
- **Multi-Omics MOFA2** — `value = "mofa_tab"`
- **AI Analysis** — Gemini chat

### Comparison Selector Sync
Four synchronized selectors: `contrast_selector` (DE Dashboard), `contrast_selector_signal`, `contrast_selector_grid`, `contrast_selector_pvalue`. Bidirectional sync — changing any updates all.

## Development Workflow

### Running Locally
```r
shiny::runApp('/Users/brettphinney/Documents/claude/', port=3838, launch.browser=TRUE)
```
- **DO NOT** use `source()` — it doesn't work properly in VS Code
- No hot-reload — must restart after every code change
- Stop: `pkill -f "shiny::runApp"` | Check: `lsof -i :3838`

## Deployment

### Three Platforms
1. **GitHub** (`origin`) — Source code. `git push origin main` auto-syncs to HF via GitHub Actions.
2. **Hugging Face** (`hf`) — Docker app. Thin `Dockerfile` FROM `brettphinney/delimp-base:v3.1`.
3. **HPC** — Apptainer containers (see `HPC_DEPLOYMENT.md`)

### README Management (CRITICAL)
- Edit `README_GITHUB.md` for GitHub, `README_HF.md` for HF
- **NEVER** push README.md changes to both remotes
- **NEVER** use `git add .` when README.md is modified

### Docker Base Image
- `brettphinney/delimp-base:v3.1` on Docker Hub (public, ~5 GB)
- Adding new R packages requires rebuilding base image on Windows box
- Code-only changes: just `git push origin main`
- **Windows update shortcut**: `bash update_docker.sh` (pulls latest + rebuilds container)

## UI Design Patterns

- **`page_navbar` layout**: Dark navbar with white text (CSS `!important`). `nav_spacer()` + `nav_item()` for gear icon. Hover dropdowns via `.navbar .dropdown:hover > .dropdown-menu { display: block; }`. Active tab gets teal underline.
- **bslib `navbar_options()` required**: bslib 0.9.0+ deprecated `bg` as direct arg to `page_navbar()`. Use `navbar_options = navbar_options(bg = ...)`.
- **Sidebar accordion**: Three collapsible panels: "Upload Data" (open), "Pipeline Settings", "AI Chat". Conditional phospho/XIC sections use separate `accordion()` blocks.
- **DE Dashboard sub-tabs**: `navset_card_tab(id = "de_dashboard_subtabs")` — Volcano+heatmap, Results Table, PCA, CV Analysis.
- **CRITICAL bslib issue**: `card()`/`card_body()` don't render at top level inside `nav_panel()`. Use plain `div()` with inline CSS.
- **CRITICAL bslib sub-tab issue**: `renderUI`/`uiOutput` content disappears inside `navset_card_tab` sub-tabs. `renderPlot` crashes with `invalid quartz() device size` on macOS (0-width hidden container). **Use `plotlyOutput`/`renderPlotly`** — only reliable output type in bslib sub-tabs. See @docs/PATTERNS.md for details.
- **Info modal pattern**: `actionButton("[id]_info_btn", icon("question-circle"), class="btn-outline-info btn-sm")` + `observeEvent(...)`.
- **Plotly annotations**: Use `layout(annotations = ...)` with paper coordinates, not ggplot `annotate()`. For summary stats, prefer ggplot subtitles over plotly annotation cards (more robust in bslib sub-tabs).
- **Scrollable tab content**: Wrap dense sub-tab content in `div(style = "overflow-y: auto; max-height: calc(100vh - 200px);")` with `min-height` on key widgets to prevent bslib compression.
- Plot heights use viewport-relative units (`vh`, `calc()`) — no fixed pixel heights.

## Key Gotchas

| Problem | Solution |
|---------|----------|
| Navbar text invisible on dark bg | Flatly theme needs CSS override: `.navbar .nav-link { color: rgba(255,255,255,0.75) !important; }` |
| Hidden tabs show letter fragments | `.navbar .nav-item[style*='display: none'] { width: 0 !important; overflow: hidden !important; }` |
| `page_navbar(bg=...)` deprecation | Use `navbar_options = navbar_options(bg = ...)` (bslib 0.9.0+) |
| `source()` doesn't start app | Use `shiny::runApp()` instead |
| Selections disappear after clicking | Reactive loop — table must not depend on selection-derived reactives |
| bslib `card()` doesn't render | Use plain `div()` for top-level nav_panel content |
| Volcano P.Value vs adj.P.Val mismatch | Y-axis uses raw P.Value for spread; dashed line at `max(P.Value)` among adj.P.Val < 0.05 proteins |
| `arrow::select` masks `dplyr::select` | Use `dplyr::select()` explicitly |
| Shiny hidden input not registered by JS | Use `div(style="display:none;", radioButtons(...))` for `conditionalPanel` |
| `readDIANN` data.table column error | Must pass `format="parquet"` for .parquet files |
| `return()` inside `withProgress` | Exits `withProgress` not enclosing function. Use flat `tryCatch`. |
| Quarto `output_file` path error | Pass filename only, then `file.rename()` to target dir |
| y_protein `colSums` error | It's limma `EList`. Extract `$E` for expression matrix. |
| SQLite `Parameter N does not have length 1` | Use `NA_character_` instead of `NULL` |
| SSH output encoding crash | `iconv(..., sub="")` in `ssh_exec`/`scp_download`/`scp_upload` |
| R regex `\\s` invalid | Use `[:space:]` in base R regex (POSIX ERE) |

## Version History

Current version: **v3.1.1** (2026-02-26). See [CHANGELOG.md](CHANGELOG.md) for details.

Key decisions: Modularization (v2.3) | XIC Viewer (v2.1) | Phospho Phase 1 (v2.4) | GSEA multi-DB (v2.5) | SSH job submission (v2.5) | Docker backend (v3.0) | MOFA2 (v3.0) | Core Facility (v3.1) | **UI overhaul to page_navbar** (v3.1) | Volcano/CV fixes + Export panel (v3.1.1)
