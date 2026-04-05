# DE-LIMP Project Context for Claude

## Working Preferences
- **Update this file** when new patterns, gotchas, or architectural decisions emerge
- For detailed change history, update `CHANGELOG.md` (not this file)
- **Document as you go**: When the user says "wrap up", "good night", "that's it for now", or asks for a summary — update CLAUDE.md and CHANGELOG.md with all changes from the current work before responding
- **NEVER run heavy computation on HPC login nodes** — always submit via `sbatch` or request an interactive node with `srun`. Login nodes are shared and running CPU/memory-intensive tasks can get the user flagged.
- **Check primary sources before guessing — NEVER guess anything verifiable** — This applies to EVERYTHING: algorithms, formulas, file paths, container locations, module names, binary paths, HPC configurations, API formats, config parameters. If it can be checked, check it FIRST. SSH to the cluster and run `find`/`which`/`ls` for paths. Fetch source code from GitHub for algorithms. Read config files for parameters. Do NOT answer from memory or approximation. Examples of past failures: (1) tof-to-mz formula guessed from first principles was off by 155 Da — correct formula was in timsrust source; (2) said `module load diann` when DIA-NN actually runs as an Apptainer container at `/quobyte/proteomics-grp/apptainers/diann2.3.0.sif`; (3) used depthcharge's default peak filtering thinking it was Cascadia's — Cascadia's actual `train()` function uses different preprocessing.

## Review Agents (spawn before major releases)
After significant changes, spawn these 5 review agents in parallel:
1. **Biological researcher** — workflow intuitiveness, jargon, missing biology features (non-bioinformatician perspective)
2. **Proteomics expert** — DIA-NN integration, QC, core facility readiness, instrument support
3. **Statistician** — statistical validity, multiple testing, no-replicates caveats, methodology
4. **Error handling & UX audit** — silent failures, blank req() screens, missing validation
5. **Documentation audit** — stale references, version mismatches, missing features across all docs
- Detailed patterns: @docs/PATTERNS.md | TODO list: @docs/TODO.md

## Project Overview
DE-LIMP is a Shiny proteomics data analysis pipeline using the LIMPA R package for differential expression analysis of DIA-NN data.

- **GitHub**: https://github.com/bsphinney/DE-LIMP
- **Hugging Face**: https://huggingface.co/spaces/brettsp/de-limp-proteomics
- **Local URL**: http://localhost:3838
- **R**: 4.5+ required (limpa needs Bioconductor 3.22+)

## Architecture

### Structure (~9,000 lines total)
```
app.R (~350 lines):      Package loading, backend detection (Docker/HPC/Core Facility/Apptainer), VERSION + stats loading, reactive values, module calls, SSH auto-connect, container detection
R/ui.R (~1,700 lines):   build_ui() — page_navbar layout, CSS/JS, accordion sidebar, all tab nav_panels, SSH file browser modals, environment badge
R/server_*.R (12 files): Server modules, each receives (input, output, session, values, ...)
R/helpers*.R (6 files):  Pure utility functions (no Shiny reactivity)
```

### Key Files

| File | Purpose |
|------|---------|
| `app.R` | Orchestrator — package loading, backend detection (Docker/HPC/Apptainer/Core Facility), SSH auto-connect, container detection, reactive values, module calls |
| `R/ui.R` | `page_navbar` layout, accordion sidebar, all tab definitions, SSH file browser modals, environment badge |
| `R/server_data.R` | Data upload, example load, group assignment, pipeline execution, contaminant analysis, no-replicates mode |
| `R/server_de.R` | Volcano, DE table, heatmap, CV analysis, selection sync |
| `R/server_qc.R` | QC sample metrics (faceted trend plot), diagnostic plots, p-value distribution, data completeness |
| `R/server_viz.R` | Expression grid (contaminant highlighting), signal distribution (contaminant overlay), PCA |
| `R/server_gsea.R` | GSEA analysis, multi-DB (BP/MF/CC/KEGG), organism detection |
| `R/server_ai.R` | AI Summary (all contrasts), Data Chat, Gemini integration, HTML report export |
| `R/server_search.R` | Docker/HPC dual backend, SSH, DIA-NN search, job queue, SSH file browser, NCBI download, SLURM proxy, Load from HPC |
| `R/server_phospho.R` | Phospho site-level DE, volcano, site table |
| `R/server_mofa.R` | MOFA2 multi-view integration |
| `R/server_comparator.R` | Run Comparator: cross-tool DE comparison (DE-LIMP vs DE-LIMP/Spectronaut/FragPipe), 4-layer diagnostics, 9-rule hypothesis engine (Rule 0: 0-ratio rescue), Spectronaut ZIP parser (TopN/Quant3/RunQC/n_ratios/AnalysisOverview, Spectronaut 20+ key-value RunOverview), contrast mismatch detection, instrument context in AI prompts, DPC-Quant methodology note in Claude export, MOFA2 decomposition |
| `R/helpers_denovo.R` | Cascadia de novo: SSL parsing, peptide classification, DIAMOND BLAST, sbatch generation (feature branch) |
| `R/helpers_dda.R` | DDA pipeline: Sage config, result parsing, Casanovo mztab, classify_dda_denovo, sbatch generation (feature branch) |
| `R/server_dda.R` | DDA search module: Sage SLURM, Casanovo GPU, DIAMOND BLAST, species viz, contaminant filtering, info modals (feature branch) |
| `R/server_denovo_viz.R` | Advanced de novo viz: BLAST alignment, target-decoy FDR, cross-species, protein families, coverage maps (feature branch) |
| `R/server_denovo_controls.R` | De novo controls: confidence slider, manuscript stats, GO annotation, disagreement analysis (feature branch) |
| `R/server_facility.R` | Core facility: reports, job history, QC dashboard |
| `R/server_session.R` | Info modals, save/load session, reproducibility, About tab, unified history, notes, remote history |
| `R/helpers_search.R` | `ssh_exec()`, `build_diann_flags()`, `generate_sbatch_script()`, `generate_parallel_scripts()`, `generate_search_info()`, `check_cluster_resources()`, UniProt/NCBI search, unified activity log, SSH file browser helpers, SLURM proxy |
| `R/helpers_instrument.R` | `parse_timstof_metadata()`, `parse_thermo_metadata()`, `parse_raw_file_metadata()`, `extract_tic_timstof()`, `compute_tic_metrics()`, `diagnose_run()`, instrument formatters for Methods/AI |
| `VERSION` | Single-line app version (e.g. `3.7.0`), read at startup into `values$app_version` |
| `stats/community_stats.json` | GitHub traffic data generated daily by `.github/workflows/track-stats.yml` |
| `Launch_DE-LIMP_Docker.bat` | Windows one-click Docker launcher with SSH key auto-detection and shared PC support |
| `launch_delimp.sh` | Mac/Linux launcher — auto-downloads `hpc_setup.sh`, handles container install + SSH tunnel |
| `hpc_setup.sh` | HPC setup script — container install, R packages, code updates, SLURM proxy, per-user dirs |

### Tab Structure (page_navbar)
Navbar: **New Search** (conditional) | **QC** | **Analysis** dropdown | **Output** dropdown (Export Data, Methods & Code) | **About** dropdown (Community, Search History, Analysis History) | **Education** | **Facility** dropdown (conditional) | gear icon (far right)

- `page_navbar(id = "main_tabs", navbar_options = navbar_options(bg = "#2c3e50"))` — dark navbar, global sidebar, hover dropdowns
- Dropdown section labels ("Setup"/"Results"/"AI") injected via JS

**Progressive reveal**: `nav_hide()`/`nav_show()` on `"main_tabs"`. Hidden on startup via `session$onFlushed(once=TRUE)`:
- **Always visible**: New Search (if `search_enabled`), Analysis > Data Overview, About, Education, Facility (if `is_core_facility`)
- **QC**: shown when `values$raw_data` not NULL OR `values$tic_traces` not NULL. Sub-tabs: TIC traces (faceted/overlay/metrics), Data Completeness (detection vs inferred analysis, Jaccard dendrogram).
- **DE Dashboard, GSEA, MOFA2, Run Comparator, AI Analysis, Output**: shown when `values$fit` not NULL (Comparator also shown when `values$comparator_results` exists). PCA and Expression Grid also work without `values$fit` when quantification is complete (no-replicates mode).
- **Phosphoproteomics**: shown when `values$phospho_detected$detected` is TRUE

**Tab values that MUST NOT change** (used by server nav_select/nav_show/nav_hide):
`"QC"`, `"DE Dashboard"`, `"Gene Set Enrichment"`, `"mofa_tab"`, `"comparator_tab"`, `"AI Analysis"`, `"Output"`, `"Phosphoproteomics"`, `"Data Overview"`, `"data_overview_tabs"`, `"Assign Groups & Run"`, `"about_tab"`, `"history_tab"`

#### Analysis dropdown
- **Data Overview** — `navset_card_tab(id = "data_overview_tabs")`: Assign Groups & Run, Signal Distribution, Dataset Summary, Replicate Consistency, Contaminant Analysis, Expression Grid, Data Explorer, Data Completeness, AI Summary
- **DE Dashboard** — `navset_card_tab(id = "de_dashboard_subtabs")`: Volcano (+heatmap), Results Table, PCA, CV Analysis. Comparison selector banner above sub-tabs.
- **Phosphoproteomics** — conditional on phospho detection
- **Gene Set Enrichment** — BP/MF/CC/KEGG with per-ontology caching
- **Multi-Omics MOFA2** — `value = "mofa_tab"`
- **Run Comparator** — `value = "comparator_tab"`, `navset_card_tab(id = "comparator_subtabs")`: Settings Diff, Protein Universe, Quantification, DE Concordance, AI Analysis. Modes: DE-LIMP vs DE-LIMP/Spectronaut/FragPipe.
- **AI Analysis** — Gemini chat

#### About dropdown
- **Community** (`value = "about_tab"`) — Version, GitHub stats cards, trend sparklines, recent discussions, links
- **History** (`value = "history_tab"`) — Unified activity log (replaces Search History + Analysis History). Single DT table with expandable detail rows, project/status filters, Load/Settings/Notes/Project buttons. Notes modal on search completion.

### Unified Activity Log (v3.6.0)
- **Single CSV**: `activity_log_path()` → shared HPC storage or local `~/.delimp_activity_log.csv`. 33 columns, append-only with file locking. Updates via `update_activity()` (by `output_dir`) or `update_activity_by_id()`. Replaces old search_history.csv + analysis_history.csv + projects.json (migrated to `.bak` on first load).
- **Event types**: `search_submitted`, `search_completed`, `search_failed`, `analysis_completed`, `data_loaded`, `session_restored`. Record points: job submission/completion (`server_search.R`), pipeline completion (`server_data.R`), session upload (`server_session.R`).
- **Projects & Notes**: `project` field in CSV (not separate JSON). Notes modal on search completion; editable anytime via pen icon. Deterministic RDS at `{output_dir}/session.rds`.
- **History UI**: DT expandable rows with Log/Load/Settings/Notes/Project buttons. "Load" tries `session.rds` first, falls back to `report.parquet`.
- **n_proteins**: Use `length(unique(raw_data$genes$Protein.Group))` not `nrow(raw_data$E)` — the latter counts precursors (~40k), not protein groups (~3k).

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

### Four Deployment Modes
1. **GitHub** (`origin`) — Source code. `git push origin main` auto-syncs to HF via GitHub Actions.
2. **Hugging Face** (`hf`) — Docker app. Thin `Dockerfile` FROM `brettphinney/delimp-base:v3.1`.
3. **Docker + SSH** (recommended for Windows) — `Launch_DE-LIMP_Docker.bat` runs DE-LIMP locally in Docker, connects to HPC via SSH for DIA-NN search. Shared PC support with auto SSH key detection.
4. **HPC Apptainer** (alternative) — `launch_delimp.sh` / `Launch_DE-LIMP.bat` launches via Apptainer on HPC with SLURM proxy. See `HPC_DEPLOYMENT.md`.

### Release Checklist
On each version release, do ALL of these:
1. Bump `VERSION` file
2. Update `CHANGELOG.md` with new section
3. Update `README_GITHUB.md` with new features → copy to `README.md`
4. Update `README_HF.md` with new features
5. Update `docs/index.html` version badge and feature cards (GitHub Pages Education site)
6. Create GitHub release: `gh release create vX.Y.Z --title "..." --notes "..."`
7. Run review agents (biologist, proteomics expert, statistician, error audit, docs audit)

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
- **SVG vector export**: Plotly plots have camera icon for SVG download via `config(toImageButtonOptions = list(format = "svg", scale = 2))`. ggplot/ComplexHeatmap plots use `downloadButton` with `ggsave(device = "svg")` or `svg()` device.
- Plot heights use viewport-relative units (`vh`, `calc()`) — no fixed pixel heights.

## Key Gotchas

Full gotchas table: @docs/GOTCHAS.md (70+ items organized by category)
HPC paths & containers: @docs/HPC_PATHS.md

**Most critical gotchas** (ones that waste the most time):
- **Two DIA-NN containers**: `/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif` has .NET (reads .raw). `/apptainers/diann2.3.0.sif` does NOT — .raw files silently skipped.
- **`uiOutput` vanishes in bslib sub-tabs**: Use `plotlyOutput`/`renderPlotly` only. See @docs/PATTERNS.md.
- **NEVER use mounted drives for app state**: All app state in `~/.delimp_*`. Shared data via SSH/SCP.
- **Derived data stays with source data**: session.rds, TIC cache in output dir, not home dir.
- **DIA-NN binary inside container**: `/diann-2.3.0/diann-linux` (NOT just `diann`)

### Queue Switching
Auto-switches parallel jobs between `genome-center-grp/high` and `publicgrp/low` partitions. Steps 2/4 (array) move to low (preemptible); steps 1/3/5 (assembly) stay on high. See @docs/QUEUE_SWITCHING.md for full logic, known issues, and SLURM state mapping.

## Version History

Current version: **v3.7.0** — defined in `VERSION` file. See [CHANGELOG.md](CHANGELOG.md) for details.

Key decisions: Modularization (v2.3) | XIC Viewer (v2.1) | Phospho Phase 1 (v2.4) | GSEA multi-DB (v2.5) | SSH job submission (v2.5) | Docker backend (v3.0) | MOFA2 (v3.0) | Core Facility (v3.1) | **UI overhaul to page_navbar** (v3.1) | Volcano/CV fixes + Export panel (v3.1.1) | **About tab, community stats, docs overhaul** (v3.2.0) | **Search history, log parser, Claude export enhancements, sacct fixes** (v3.2.1) | **Chromatography QC** (v3.3.0) | **Run Comparator** (v3.4.0) | **Run Comparator enhancements, Search/Analysis History, smart partitions, FASTA library fixes** (v3.5.0) | **Spectronaut parsing fixes, TopN scatter fix, sub-tab help modals** (v3.5.1) | **Unified activity log, cluster monitoring, Spectronaut NaN/rescue fixes** (v3.6.0/v3.6.1) | **Docker launcher, NCBI integration, contaminant analysis, SSH file browser, Load from HPC, no-replicates mode** (v3.7.0)

Unreleased (post-v3.7.0): Data Completeness visualization, SVG vector export, NCBI gene symbols in DE/GSEA, auto-queue switching fixes, drift test infrastructure, Cascadia de novo integration (feature branch).
