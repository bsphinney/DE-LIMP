# DE-LIMP Project Context for Claude

> **This file is the always-loaded context.** Keep it lean — it loads into *every* session. Detailed reference lives in `docs/` and is read on demand:
> - **Gotchas / known-issue fixes** → `docs/GOTCHAS.md` (index below)
> - **Subsystem patterns** (Shiny, bslib, SSH, DIA-NN flags, Comparator, Proteogenomics, UI) → `docs/PATTERNS.md`
> - **Proteogenomics architecture** → `docs/PROTEOGENOMICS.md`
> - **HPC paths/containers** → `docs/HPC_PATHS.md` · **Queue switching** → `docs/QUEUE_SWITCHING.md` · **TODO** → `docs/TODO.md`

## Working Preferences
- **Update the right file**: project state/patterns → this file or `docs/`; change history → `CHANGELOG.md`; new gotcha → `docs/GOTCHAS.md`.
- **Bump the patch version after every user-visible fix**: update (a) the `VERSION` file, (b) the `# Version:` line in the `app.R` header comment, and (c) add a `CHANGELOG.md` entry under the new version. (a) drives the runtime console banner; (b) is what the user sees opening `app.R` in the editor. Keep them in sync.
- **NEVER run heavy computation on HPC login nodes** — submit via `sbatch` or request an interactive node with `srun`. Login nodes are shared; CPU/memory-heavy tasks can get the user flagged.
- **Check primary sources before guessing — NEVER guess anything verifiable.** Applies to EVERYTHING: algorithms, formulas, file paths, container locations, module names, binary paths, HPC config, API formats, parameters. If it can be checked, check it FIRST — SSH and run `find`/`which`/`ls`, fetch source from GitHub, read config files. Do NOT answer from memory. Past failures: tof-to-mz formula guessed from first principles was off by 155 Da (correct one was in timsrust source); claimed `module load diann` when DIA-NN is an Apptainer container; used depthcharge's default peak filtering thinking it was Cascadia's.

## Architectural rules (NEVER violate — discovered the hard way in v3.9.x)

These four rules exist because each was violated in early DE-LIMP and produced wrong-but-plausible exports that misled real analyses. Before writing code that touches user-facing exports, methods text, AI prompts, info modals, or reproducibility logs — stop and check this list.

1. **Pipeline objects must self-describe — never hardcode a description of "what we did".** Every quantification path returns a self-describing object; downstream consumers read `$pipeline_id`, `$methods_paragraph`, `$rollup_method` from it. Never write `if (isTRUE(values$pipeline_mode_used == "maxlfq")) ...` in a new file — put that branch in a single helper. Adding a third pipeline must require zero edits to downstream files. (Real users once got reports describing the wrong pipeline because methods/AI/Comparator/repro text all hardcoded `"DPC-Quant... dpcDE()"`.)

2. **`%||%` defaults that flow into user-facing text must be tagged.** For any value that ends up in an export / AI prompt / methods string, either render `(DEFAULT — not user-confirmed)` next to it, or replace the `%||%` with `NA_character_` and have the prompt-builder skip the line. Never silently substitute a fabricated value (e.g. `coalesce_setting(..., "0.05")`) for missing user input — a reviewer reading "FDR threshold: 0.05" wrongly assumes the user set it.

3. **Concepts have one definition.** Every shared concept (detection status, pipeline label, covariate display name, DIA-NN defaults, classifier rules) lives in exactly **one** file as a function/constant, imported elsewhere. If you find yourself copying logic to a second site, refactor instead. ("Detected vs Inferred" once had 4 independent classifiers; the covariate display name is read 22 times.)

4. **Silent catch is banned in export paths.** In every multi-file export bundler use `safe_section(manifest, name, expr)` from `R/helpers.R` — records `[OK]` on success, `[SKIPPED] <name> -- <reason>` to a `MANIFEST.txt` in the ZIP root on failure. Never write `tryCatch(error = function(e) NULL)` around CSV writes / ZIP entries — a downstream user gets a ZIP silently missing whole sections.

When reviewing your own changes: ask whether the export still describes the analysis correctly under *every* pipeline + *every* input shape. If you can't answer "yes" with certainty, branch the code or add a test before shipping.

## Review Agents (spawn before major releases)
After significant changes, spawn these 5 in parallel: (1) **Biological researcher** — workflow intuitiveness, jargon, missing biology features; (2) **Proteomics expert** — DIA-NN integration, QC, core-facility readiness, instrument support; (3) **Statistician** — statistical validity, multiple testing, no-replicates caveats; (4) **Error handling & UX audit** — silent failures, blank `req()` screens, missing validation; (5) **Documentation audit** — stale references, version mismatches, missing features.

## Project Overview
DE-LIMP is a Shiny proteomics data analysis pipeline using the LIMPA R package for differential expression analysis of DIA-NN data.
- **GitHub**: https://github.com/bsphinney/DE-LIMP · **Hugging Face**: https://huggingface.co/spaces/brettsp/de-limp-proteomics
- **Local URL**: http://localhost:3838 · **R**: 4.5+ required (limpa needs Bioconductor 3.22+)

## Architecture

### Structure (~9,000 lines)
```
app.R (~350 lines):      Package loading, backend detection (Docker/HPC/Core Facility/Apptainer), VERSION + stats loading, reactive values, module calls, SSH auto-connect, container detection
R/ui.R (~1,700 lines):   build_ui() — page_navbar layout, CSS/JS, accordion sidebar, all tab nav_panels, SSH file browser modals, environment badge
R/server_*.R:            Server modules, each receives (input, output, session, values, ...)
R/helpers*.R:            Pure utility functions (no Shiny reactivity)
```

### Key Files

| File | Purpose |
|------|---------|
| `app.R` | Orchestrator — package loading, backend detection, SSH auto-connect, container detection, reactive values, module calls |
| `R/ui.R` | `page_navbar` layout, accordion sidebar, all tab definitions, SSH file browser modals, environment badge |
| `R/server_data.R` | Data upload, example load, group assignment, pipeline execution, contaminant analysis, no-replicates mode |
| `R/server_de.R` | Volcano, DE table, heatmap, CV analysis, selection sync |
| `R/server_qc.R` | QC sample metrics (faceted trend plot), diagnostic plots, p-value distribution, data completeness |
| `R/server_viz.R` | Expression grid (contaminant highlighting), signal distribution (contaminant overlay), PCA |
| `R/server_gsea.R` | GSEA, multi-DB (BP/MF/CC/KEGG), organism detection |
| `R/server_ai.R` | AI Summary, Data Chat, Gemini integration, HTML report export |
| `R/server_search.R` | Docker/HPC dual backend, SSH, DIA-NN search, job queue, SSH file browser, NCBI download, SLURM proxy, Load from HPC |
| `R/server_phospho.R` | Phospho site-level DE, volcano, site table |
| `R/server_mofa.R` | MOFA2 multi-view integration |
| `R/server_comparator.R` | Run Comparator: cross-tool DE comparison, 4-layer diagnostics, 9-rule hypothesis engine, Spectronaut ZIP parser, MOFA2 decomposition (details in PATTERNS.md) |
| `R/server_facility.R` | Core facility: reports, job history, QC dashboard |
| `R/server_session.R` | Info modals, save/load session, reproducibility, About tab, unified history, notes |
| `R/server_proteog_builder.R` + `helpers_proteog_assembly.R` / `helpers_rnaseq.R` / `helpers_slims.R` | Proteogenomics RNA-seq → FASTA pipeline (see `docs/PROTEOGENOMICS.md`) |
| `R/helpers_denovo.R` / `helpers_dda.R` / `server_dda.R` / `server_denovo_viz.R` / `server_denovo_controls.R` | Cascadia/Sage/Casanovo de novo + DDA (feature branch) |
| `R/helpers_search.R` | `ssh_exec()`, `build_diann_flags()`, `generate_sbatch_script()`, `generate_parallel_scripts()`, `generate_search_info()`, `check_cluster_resources()`, UniProt/NCBI search, unified activity log, SSH browser helpers, SLURM proxy |
| `R/helpers_instrument.R` | `parse_timstof_metadata()`, `parse_thermo_metadata()`, `extract_tic_timstof()`, `compute_tic_metrics()`, `diagnose_run()`, instrument formatters |
| `VERSION` | Single-line app version, read at startup into `values$app_version` |
| `stats/community_stats.json` | GitHub traffic data generated daily by `.github/workflows/track-stats.yml` |
| `Launch_DE-LIMP_Docker.bat` / `launch_delimp.sh` / `hpc_setup.sh` | Windows Docker launcher / Mac-Linux launcher / HPC setup script |

### Tab Structure (page_navbar)
Navbar: **New Search** (conditional) | **QC** | **Analysis** dropdown | **Output** dropdown | **About** dropdown | **Education** | **Facility** dropdown (conditional) | gear icon.
- `page_navbar(id = "main_tabs", navbar_options = navbar_options(bg = "#2c3e50"))` — dark navbar, global sidebar, hover dropdowns. Dropdown section labels injected via JS.
- **Analysis dropdown**: Data Overview (`navset_card_tab id="data_overview_tabs"`: Assign Groups & Run, Signal Distribution, Dataset Summary, Replicate Consistency, Contaminant Analysis, Expression Grid, Data Explorer, Data Completeness, AI Summary) · DE Dashboard (`id="de_dashboard_subtabs"`: Volcano+heatmap, Results Table, PCA, CV Analysis) · Phosphoproteomics · Gene Set Enrichment · MOFA2 · Run Comparator (`id="comparator_subtabs"`) · AI Analysis.
- **About dropdown**: Community (`about_tab`) · History (`history_tab` — unified activity log).

**Progressive reveal** via `nav_hide()`/`nav_show()` on `"main_tabs"`, hidden on startup via `session$onFlushed(once=TRUE)`:
- **Always visible**: New Search (if `search_enabled`), Analysis > Data Overview, About, Education, Facility (if `is_core_facility`).
- **QC**: when `values$raw_data` OR `values$tic_traces` not NULL.
- **DE Dashboard, GSEA, MOFA2, Comparator, AI Analysis, Output**: when `values$fit` not NULL (Comparator also when `values$comparator_results` exists). PCA + Expression Grid also work without `values$fit` (no-replicates mode).
- **Phosphoproteomics**: when `values$phospho_detected$detected` is TRUE.

**Tab values that MUST NOT change** (used by server `nav_select`/`nav_show`/`nav_hide`):
`"QC"`, `"DE Dashboard"`, `"Gene Set Enrichment"`, `"mofa_tab"`, `"comparator_tab"`, `"AI Analysis"`, `"Output"`, `"Phosphoproteomics"`, `"Data Overview"`, `"data_overview_tabs"`, `"Assign Groups & Run"`, `"about_tab"`, `"history_tab"`, `"search_tab"` (proteogenomics DIA-NN Run Search sub-panel), `"build_database_tab"` (proteogenomics Build Database sub-panel, HPC-gated).

### Unified Activity Log (v3.6.0)
- **Single CSV**: `activity_log_path()` → shared HPC storage or local `~/.delimp_activity_log.csv`. 33 columns, append-only with file locking. Updates via `update_activity()` (by `output_dir`) or `update_activity_by_id()`. Replaces old search/analysis history + projects.json.
- **Event types**: `search_submitted`, `search_completed`, `search_failed`, `analysis_completed`, `data_loaded`, `session_restored`.
- **Projects & Notes**: `project` field in CSV; notes modal on search completion, editable via pen icon. Deterministic RDS at `{output_dir}/session.rds`.
- **History UI**: DT expandable rows with Log/Load/Settings/Notes/Project buttons. "Load" tries `session.rds` first, falls back to `report.parquet`.

### Comparison Selector Sync
Four synchronized selectors (`contrast_selector` on DE Dashboard, plus `_signal`/`_grid`/`_pvalue`). Bidirectional — changing any updates all.

### Proteogenomics workflow (v3.11.0)
RNA-seq → 11-stage SLURM chain → custom FASTA → FASTA-library catalog → main-page picker. Architecture in `docs/PROTEOGENOMICS.md`; implementation patterns in `docs/PATTERNS.md`.

## Development Workflow

### Running Locally
```r
shiny::runApp('/Users/brettphinney/Documents/claude/', port=3838, launch.browser=TRUE)
```
- **DO NOT** use `source()` — doesn't work in VS Code. No hot-reload — restart after every code change.
- Stop: `pkill -f "shiny::runApp"` · Check: `lsof -i :3838`

## Deployment

### Four Deployment Modes
1. **GitHub** (`origin`) — source code; `git push origin main` auto-syncs to HF via GitHub Actions.
2. **Hugging Face** (`hf`) — Docker app; thin `Dockerfile` FROM `brettphinney/delimp-base:v3.1`.
3. **Docker + SSH** (recommended for Windows) — `Launch_DE-LIMP_Docker.bat` runs locally in Docker, connects to HPC via SSH. Shared-PC support with auto SSH key detection.
4. **HPC Apptainer** — `launch_delimp.sh` / `Launch_DE-LIMP.bat` via Apptainer with SLURM proxy. See `HPC_DEPLOYMENT.md`.

### Release Checklist
1. Bump `VERSION` · 2. Update `CHANGELOG.md` · 3. Update `README_GITHUB.md` → copy to `README.md` · 4. Update `README_HF.md` · 5. Update `docs/index.html` version badge + feature cards · 6. `gh release create vX.Y.Z` · 7. Run the 5 review agents.

### README Management (CRITICAL)
- Edit `README_GITHUB.md` for GitHub, `README_HF.md` for HF.
- **NEVER** push README.md changes to both remotes. **NEVER** `git add .` when README.md is modified.

### Docker Base Image
- `brettphinney/delimp-base:v3.1` on Docker Hub (public, ~5 GB). New R packages require rebuilding the base image on the Windows box. Code-only changes: just `git push origin main`. Windows update shortcut: `bash update_docker.sh`.

### Version Management
- **Single source of truth**: `VERSION` file → `app.R` reads at startup → `values$app_version` → all modules. No hardcoded version strings.
- Community stats (`stats/community_stats.json`) generated daily by `track-stats.yml`.

## UI Design Patterns
Critical: bslib `card()` doesn't render at top level of `nav_panel()` (use `div()`); `uiOutput`/`renderPlot` break inside `navset_card_tab` sub-tabs (use `plotlyOutput`/`renderPlotly`, or `shinyjs::html()` for dynamic HTML). Full layout/navbar/SVG-export conventions in `docs/PATTERNS.md` ("UI Design Patterns" + "bslib navset_card_tab Rendering Issues").

## Key Gotchas
Quick-reference fix tables live in **`docs/GOTCHAS.md`**, grouped by subsystem. Read the relevant table before debugging in that area; add a row when you fix a new non-obvious bug. Sections:
- **R Shiny / bslib** — reactive loops, `withProgress` + `<<-`, `observeEvent(list(...))` auto-fire, sub-tab rendering, SQLite/Quarto quirks
- **DIA-NN** — Genes-column accessions, parquet libs, `--quant-ori-names`, mass-acc + `--use-quant`, the two containers (`.raw` support), binary path
- **Data & Columns** — precursors vs protein groups, EList `$E`, arrow/dplyr mask, Linux matrix subsetting, contaminants, no-replicates, NCBI gene mapping
- **SSH / HPC / SLURM** — ControlMaster socket length, `sacct` substeps, QOS vs assoc limits, per-user CPU cap, `nohup ... </dev/null`, mounted-drive ban, derived-data placement
- **Spectronaut** — column-name quirks, Quant3, 0-ratio NaN, RunOverview formats
- **Proteogenomics** — empty-JSON `list()`, no-job_id false-complete, reactiveValues-at-init, Ensembl URL/Plants, per-user catalog
- **Sage / Casanovo / Cascadia** — Sage v0.14.7 schema, `.raw` silent-skip, msconvert bind, Casanovo v4/v5 envs, `pipefail` glob; plus a Training subsection
- HPC paths/containers: `docs/HPC_PATHS.md`. **Queue switching** (genome-center-grp/high ↔ publicgrp/low): `docs/QUEUE_SWITCHING.md`.

## Version History
Current version: **v3.11.0** — see `VERSION` and `CHANGELOG.md` for the full history and per-release details. Unreleased work + code-audit follow-ups are tracked in `docs/TODO.md`. (Cascadia/Sage/Casanovo de novo + DDA is on a feature branch.)
