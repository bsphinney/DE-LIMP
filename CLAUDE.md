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
app.R (~290 lines):      Package loading, backend detection (Docker/HPC/Core Facility), VERSION + stats loading, reactive values, module calls
R/ui.R (~1,500 lines):   build_ui() — page_navbar layout, CSS/JS, accordion sidebar, all tab nav_panels
R/server_*.R (12 files): Server modules, each receives (input, output, session, values, ...)
R/helpers*.R (6 files):  Pure utility functions (no Shiny reactivity)
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
| `R/server_comparator.R` | Run Comparator: cross-tool DE comparison (DE-LIMP vs DE-LIMP/Spectronaut/FragPipe), 4-layer diagnostics, hypothesis engine, MOFA2 decomposition |
| `R/server_facility.R` | Core facility: reports, job history, QC dashboard |
| `R/server_session.R` | Info modals, save/load session, reproducibility, About tab, analysis history, projects |
| `R/helpers_search.R` | `ssh_exec()`, `build_diann_flags()`, `generate_sbatch_script()`, `generate_parallel_scripts()`, `generate_search_info()`, `check_cluster_resources()`, UniProt search, analysis history, projects |
| `R/helpers_instrument.R` | `parse_timstof_metadata()`, `parse_thermo_metadata()`, `parse_raw_file_metadata()`, `extract_tic_timstof()`, `compute_tic_metrics()`, `diagnose_run()`, instrument formatters for Methods/AI |
| `VERSION` | Single-line app version (e.g. `3.1.1`), read at startup into `values$app_version` |
| `stats/community_stats.json` | GitHub traffic data generated daily by `.github/workflows/track-stats.yml` |

### Tab Structure (page_navbar)
Navbar: **New Search** (conditional) | **QC** | **Analysis** dropdown | **Output** dropdown (Export Data, Methods & Code) | **About** dropdown (Community, Search History, Analysis History) | **Education** | **Facility** dropdown (conditional) | gear icon (far right)

- `page_navbar(id = "main_tabs", navbar_options = navbar_options(bg = "#2c3e50"))` — dark navbar, global sidebar, hover dropdowns
- Dropdown section labels ("Setup"/"Results"/"AI") injected via JS

**Progressive reveal**: `nav_hide()`/`nav_show()` on `"main_tabs"`. Hidden on startup via `session$onFlushed(once=TRUE)`:
- **Always visible**: New Search (if `search_enabled`), Analysis > Data Overview, About, Education, Facility (if `is_core_facility`)
- **QC**: shown when `values$raw_data` not NULL OR `values$tic_traces` not NULL
- **DE Dashboard, GSEA, MOFA2, Run Comparator, AI Analysis, Output**: shown when `values$fit` not NULL (Comparator also shown when `values$comparator_results` exists)
- **Phosphoproteomics**: shown when `values$phospho_detected$detected` is TRUE

**Tab values that MUST NOT change** (used by server nav_select/nav_show/nav_hide):
`"QC"`, `"DE Dashboard"`, `"Gene Set Enrichment"`, `"mofa_tab"`, `"comparator_tab"`, `"AI Analysis"`, `"Output"`, `"Phosphoproteomics"`, `"Data Overview"`, `"data_overview_tabs"`, `"Assign Groups & Run"`, `"about_tab"`, `"search_history_tab"`, `"analysis_history_tab"`

#### Analysis dropdown
- **Data Overview** — `navset_card_tab(id = "data_overview_tabs")`: Assign Groups & Run, Signal Distribution, Dataset Summary, Replicate Consistency, Expression Grid, AI Summary
- **DE Dashboard** — `navset_card_tab(id = "de_dashboard_subtabs")`: Volcano (+heatmap), Results Table, PCA, CV Analysis. Comparison selector banner above sub-tabs.
- **Phosphoproteomics** — conditional on phospho detection
- **Gene Set Enrichment** — BP/MF/CC/KEGG with per-ontology caching
- **Multi-Omics MOFA2** — `value = "mofa_tab"`
- **Run Comparator** — `value = "comparator_tab"`, `navset_card_tab(id = "comparator_subtabs")`: Settings Diff, Protein Universe, Quantification, DE Concordance, AI Analysis. Modes: DE-LIMP vs DE-LIMP/Spectronaut/FragPipe.
- **AI Analysis** — Gemini chat

#### About dropdown
- **Community** (`value = "about_tab"`) — Version, GitHub stats cards, trend sparklines, recent discussions, links
- **Search History** (`value = "search_history_tab"`) — DT table with expandable detail rows (enzyme, mass acc, flags), Import Settings + View Log buttons, cross-reference to Analysis History
- **Analysis History** (`value = "analysis_history_tab"`) — DT table with expandable detail rows, project filter, project management, cross-reference to Search History

### Search History
- **Search history CSV**: `search_history_path()` → shared `/Volumes/proteomics-grp/delimp/search_history.csv` or local `~/.delimp_search_history.csv`. Append-only with file locking; updates via `update_search_status()` (rewrite matching row).
- **Record points**: Job submission (`server_search.R` — all 4 backends) and job completion/failure (monitoring observer).
- **Cross-reference**: `output_dir` joins search history ↔ analysis history. Both tables show a link icon to navigate to the matching entry in the other table.
- **Import Settings**: Reads CSV row fields and applies to search UI inputs (mass acc, enzyme, mode, etc.).

### Analysis History & Projects
- **Analysis history CSV**: `analysis_history_path()` → shared `/Volumes/proteomics-grp/delimp/analysis_history.csv` or local `~/.delimp_analysis_history.csv`. Append-only with file locking.
- **Projects JSON**: `projects_path()` → shared or local `delimp_projects.json`. Maps project names to `output_dir` entries.
- **DT expandable rows**: Click a row to show full File, FASTA, Output dir, Notes in a child row. JS callback on `tbody tr td` with `row.child()` toggle. Buttons use `event.stopPropagation()`.
- **n_proteins**: Use `length(unique(raw_data$genes$Protein.Group))` not `nrow(raw_data$E)` — the latter counts precursors (~40k), not protein groups (~3k).
- **Auto-save session**: Pipeline completion in `server_data.R` auto-saves `.rds` to `output_dir` (or `~/.delimp_sessions` fallback). Path stored in `session_file` column of analysis history CSV.
- **Restore button**: Green "Restore" button in history table loads auto-saved `.rds` (full pipeline state). Outline "Load" button loads raw `report.parquet` only. Handler: `observeEvent(input$history_restore_click, ...)` in `server_session.R`.

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

### Version Management
- **Single source of truth**: `VERSION` file in repo root (e.g. `3.1.1`)
- Loaded at startup in `app.R` → stored in `values$app_version` → all modules read from there
- **No hardcoded version strings** — always use `values$app_version`
- Community stats (`stats/community_stats.json`) generated daily by `track-stats.yml` GitHub Action

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
| `uiOutput` vanishes in `navset_card_tab` | Use static HTML + `shinyjs::html("div_id", content)` for dynamic injection. `plotlyOutput` with `req()` is safe. |
| DIA-NN `Genes` column has accessions | Not gene symbols. Validate with length/pattern check. Real genes from `bitr()` UNIPROT → SYMBOL. |
| MOFA2 views need same sample names | Subset to matched pairs, assign common labels (`Sample_1`, `Sample_2`, ...) |
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
| `<<-` inside `withProgress` fails | `withProgress` uses `eval(substitute(expr), env)` — `<<-` can't find parent vars. Use `new.env()` + `<-` instead. |
| SSH rapid connections rejected (255) | HPC `MaxStartups` throttling. Batch operations into fewer SSH calls; use ControlMaster multiplexing. |
| macOS SSH ControlPath too long | Unix domain sockets limited to 104 bytes on macOS. R's `tempdir()` paths are ~105 chars. Use `/tmp/.delimp_<user>_<host>`. |
| `parse_sbatch_output` returns dirty ID | SSH stdout may have trailing `\r`/whitespace. Always `trimws()` parsed job IDs. |
| DIA-NN empirical lib is `.parquet` not `.speclib` | DIA-NN 2.0+ saves empirical libraries in parquet format. Use `empirical.parquet` in `--lib` and `--out-lib`. Predicted libs remain `.predicted.speclib`. |
| DIA-NN `--quant-ori-names` required on ALL steps | Per Vadim (DIA-NN dev): preserves original filenames in `.quant` files. Without it, container bind mount path differences cause naming mismatches between steps. |
| DIA-NN `--fasta-search`/`--predictor` Step 1 only | Including in Steps 2-5 causes full FASTA re-digest. `generate_parallel_scripts()` strips these from `step_flags`. |
| DIA-NN auto mass acc + `--use-quant` | Produces different results. `generate_parallel_scripts()` forces `mass_acc_mode = "manual"`. See @docs/PATTERNS.md for full flag reference. |
| `nrow(raw_data$E)` counts precursors not proteins | Use `length(unique(raw_data$genes$Protein.Group))` for protein group count. `y_protein$E` rows are protein groups (post-pipeline). |
| `sacct` `.extern` step falsely reports COMPLETED | `sacct` includes `.extern`/`.batch` substeps that COMPLETE even when the main job is PENDING/FAILED. `check_slurm_status()` now requests `JobID,State` format and filters out substep lines (those containing `.`). |
| Log import ignores `fr_mz`/`pr_charge` | `parse_diann_log` previously put `--max-fr-mz`, `--min-fr-mz`, `--min-pr-charge`, `--max-pr-charge` into `extra_cli_flags` instead of `params`. Now parsed via `value_map` so they flow properly into `search_params` and `build_diann_flags`. |
| Array progress sacct inflated counts | `sacct -j ARRAY_ID` returns parent job + `.extern`/`.batch` substeps for each task. Filter to only `JOBID_N` format entries: `grepl("_", jid) && !grepl("\\.", jid)`. |
| Docker container name rejected | `analysis_name` with spaces/special chars fails Docker naming rules `[a-zA-Z0-9][a-zA-Z0-9_.-]*`. Sanitize with `gsub("[^a-zA-Z0-9_.-]", "_", name)`. |
| `max_pr_mz` default was 1200 not 1800 | DIA-NN default for `--max-pr-mz` is 1800. UI and all fallbacks were incorrectly set to 1200, causing FASTA library entries and searches to use wrong range when Advanced Options wasn't opened. |
| Parallel search OOM on timsTOF | Default `mem_per_file` was 32 GB, insufficient for timsTOF DIA-PASEF. Now 64 GB. |
| TIC extraction auto-triggered | `observeEvent(list(btn, trigger))` fires when button first renders (NULL→0). Use separate `reactiveVal` trigger pattern instead. |
| Older TDF missing `SummedIntensities` | `extract_tic_timstof()` auto-detects intensity column: `SummedIntensities` → `AccumulatedIntensity` → `MaxIntensity` → any `*ntensit*`. |
| FASTA library `remote_dir` stored local paths | `fasta_library_file_paths()` validates remote paths; auto-uploads via SCP if local-only. Blocks HPC submission with local-only FASTA paths. |
| SLURM limits on QOS not associations | `sacctmgr show assoc` returns empty limits. Use `sacctmgr show qos where name={account}-{partition}-qos` to get `GrpTRES` and `MaxTRESPU`. |
| Per-user CPU limit (not account) is binding | `MaxTRESPU` (e.g., 64 CPUs) constrains individual users. `GrpTRES` (e.g., 616 CPUs) is shared across all lab members. `select_best_partition()` uses per-user limit. |

## Version History

Current version: **v3.5.0** — defined in `VERSION` file. See [CHANGELOG.md](CHANGELOG.md) for details.

Key decisions: Modularization (v2.3) | XIC Viewer (v2.1) | Phospho Phase 1 (v2.4) | GSEA multi-DB (v2.5) | SSH job submission (v2.5) | Docker backend (v3.0) | MOFA2 (v3.0) | Core Facility (v3.1) | **UI overhaul to page_navbar** (v3.1) | Volcano/CV fixes + Export panel (v3.1.1) | **About tab, community stats, docs overhaul** (v3.2.0) | **Search history, log parser, Claude export enhancements, sacct fixes** (v3.2.1) | **Chromatography QC** (v3.3.0) | **Run Comparator** (v3.4.0) | **Run Comparator enhancements, Search/Analysis History, smart partitions, FASTA library fixes** (v3.5.0)
