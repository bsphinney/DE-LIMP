# DE-LIMP Project Context for Claude

## Working Preferences
- **Update this file** when new patterns, gotchas, or architectural decisions emerge
- For detailed change history, update `CHANGELOG.md` (not this file)
- **Bump the patch version after every user-visible fix**: After every fix or small feature, bump (a) the `VERSION` file, (b) the `# Version:` line in the `app.R` header comment, and (c) add a CHANGELOG entry under that new version section. (a) drives the runtime banner in the RStudio console; (b) is what the user sees by just opening `app.R` in the editor without running it. Keep them in sync.
- **NEVER run heavy computation on HPC login nodes** — always submit via `sbatch` or request an interactive node with `srun`. Login nodes are shared and running CPU/memory-intensive tasks can get the user flagged.
- **Check primary sources before guessing — NEVER guess anything verifiable** — This applies to EVERYTHING: algorithms, formulas, file paths, container locations, module names, binary paths, HPC configurations, API formats, config parameters. If it can be checked, check it FIRST. SSH to the cluster and run `find`/`which`/`ls` for paths. Fetch source code from GitHub for algorithms. Read config files for parameters. Do NOT answer from memory or approximation. Examples of past failures: (1) tof-to-mz formula guessed from first principles was off by 155 Da — correct formula was in timsrust source; (2) said `module load diann` when DIA-NN actually runs as an Apptainer container at `/quobyte/proteomics-grp/apptainers/diann2.3.0.sif`; (3) used depthcharge's default peak filtering thinking it was Cascadia's — Cascadia's actual `train()` function uses different preprocessing.

## Architectural rules (NEVER violate — discovered the hard way in v3.9.x)

These four rules exist because all of them were violated in early DE-LIMP and produced wrong-but-plausible exports that misled real analyses. Whenever you're about to write code that touches user-facing exports, methods text, AI prompts, info modals, or reproducibility logs — stop and check this list.

1. **Pipeline objects must self-describe — never hardcode a description of "what we did".**
   The MaxLFQ + limma pipeline (v3.9.0) was added as a parallel branch to the historic DPC-Quant pipeline. Methods text, Claude/AI exports, Comparator AI prompts, and the Reproducibility log all hardcoded `"DPC-CN... DPC-Quant... dpcDE()"` strings regardless of which pipeline ran. Real users got reports describing the wrong pipeline. The fix: every quantification path returns a self-describing object whose downstream consumers read `$pipeline_id`, `$methods_paragraph`, `$rollup_method` from. Never write `if (isTRUE(values$pipeline_mode_used == "maxlfq")) ...` in a new file. If you find yourself wanting to, the right fix is to put that branch in a single helper. Adding a third pipeline (DDA path is in flight) must require zero edits to downstream files.

2. **`%||%` defaults that flow into user-facing text must be tagged.**
   The Comparator's settings diff has a long chain of `coalesce_setting(..., "0.05")` / `%||% "0.6"` calls that fabricate plausible-looking values when the upstream field is `NULL`. A reviewer reading the AI prompt sees `"FDR threshold: 0.05"` and assumes the user set 0.05 — but the value was unrecorded and a hardcoded default got stamped in. **Rule:** for any value that ends up in an export / AI prompt / methods string, either render `(DEFAULT — not user-confirmed)` next to it, or replace the `%||%` with `NA_character_` and have the prompt-builder skip the line. Never silently substitute fabricated values for missing user input.

3. **Concepts have one definition.**
   "Detected vs Inferred" had four independent classifiers across `server_qc.R`, `server_viz.R`, `server_de.R`, `helpers.R` and the same fix had to be applied 5 times for v3.9.x. Covariate display name (`values$cov1_name %||% "Covariate1"`) is read 22 different times — renaming Covariate1 → Year touches 22 lines instead of 1. **Rule:** every shared concept (detection status, pipeline label, covariate display name, DIA-NN defaults, classifier rules) lives in exactly **one** file as a function/constant, imported elsewhere. If you find yourself copying logic to a second site, refactor instead.

4. **Silent catch is banned in export paths.**
   `tryCatch(error = function(e) NULL)` around CSV writes / ZIP entries means a downstream user receives an export ZIP missing entire sections (instrument metadata, QC, phospho summary) and has no idea. **Rule:** in every multi-file export bundler use `safe_section(manifest, name, expr)` from `R/helpers.R`. On success it records `[OK]`; on failure it records `[SKIPPED] <name> -- <reason>` to a `MANIFEST.txt` written into the ZIP root. The user reading the export downstream sees what's missing and why. If you write a new tryCatch with `error = function(e) NULL` in any export path, you owe the next reviewer an explanation — and it had better not just be "graceful degradation."

When reviewing your own changes against these rules: ask yourself whether the export still describes the analysis correctly under every pipeline + every input shape. If you can't answer "yes" with certainty, branch the code or add a test before shipping.

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
`"QC"`, `"DE Dashboard"`, `"Gene Set Enrichment"`, `"mofa_tab"`, `"comparator_tab"`, `"AI Analysis"`, `"Output"`, `"Phosphoproteomics"`, `"Data Overview"`, `"data_overview_tabs"`, `"Assign Groups & Run"`, `"about_tab"`, `"history_tab"`, `"search_tab"` (proteogenomics Phase D — DIA-NN Run Search sub-panel of New Search dropdown), `"build_database_tab"` (proteogenomics Phase D — Build Database sub-panel, HPC-gated)

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

### Proteogenomics workflow (v3.11.0)
- **Module files**: `R/server_proteog_builder.R` (orchestrator + UI), `R/helpers_proteog_assembly.R` (FASTA concat + dedupe), `R/helpers_rnaseq.R` (stage sbatch generators), `R/helpers_slims.R` (SLIMS + ENA download launchers)
- **Public functions**: `submit_proteogenomics_build()`, `submit_assemble_only()`, `poll_proteog_build_status()`, `cancel_proteog_build()`, `assemble_proteogenomics_fasta()`. All accept `ssh_config = NULL`.
- **11-stage SLURM chain**: `fastp → rrna_filter → star → qc_gate → stringtie → merge → gffcompare → gffread → transdecoder → rewrite → assemble`. Each is its own sbatch chained via `--dependency=afterok`. `PROTEOG_STAGE_ORDER` constant is the contract.
- **Per-build output**: `/quobyte/proteomics-grp/de-limp/rnaseq/<project_name>/` with `status.json` (self-describing schema: pipeline_id, project_name, sample_names, reference_key, read_length_tier, stages[], build_metadata), sbatch/ scripts, logs/, per-stage output dirs.
- **Final FASTA output**: `/quobyte/proteomics-grp/de-limp/databases/proteogenomics/<project>_proteogenomics_<YYYY_MM>.fasta`. Predicted ORFs + optional UniProt FASTA concatenated by `assemble` stage on Hive.
- **Reference registry**: `/quobyte/proteomics-grp/de-limp/references/registry.json` — one entry per species with `organism`, `build`, `genome_fasta`, `gtf`, `star_index`, `rrna_index` paths. Add new references via `references/scripts/build_reference_genome.sh` (downloads from Ensembl, builds STAR + bowtie2 rRNA indexes) → writes to `references/registry_pending/` → merge via `references/scripts/merge_registry_pending.sh`.
- **FASTA library catalog auto-registration**: `poll_proteog_build_status()` detects `current_stage == "complete"` and writes a new entry to local `~/.delimp_fasta_library/catalog.rds` with `content_type = "proteogenomics"` + 9 proteog-specific extension fields (`proteog_pipeline_id`, `proteog_project_dir`, `proteog_methods_paragraph`, `proteog_sample_names`, `proteog_reference_key`, `proteog_uniprot_fasta`, `proteog_read_length_tier`, `proteog_status_path`). The `library_entry_id` is stored back in `status.json` so polling doesn't re-register.
- **Shared FASTA + per-user catalog**: FASTAs live on shared Hive storage (anyone can read). Catalog entries live in each user's local `~/.delimp_fasta_library/catalog.rds` (no shared-write conflicts). "Discover from Hive" button on the search page's Proteogenomics DBs modal scans `PROTEOG_RNASEQ_ROOT/*/status.json` and populates the local catalog with any completed build it doesn't already have. Similarly "Restore from Hive" on the Build Database page restores in-progress builds.
- **Active builds persistence**: `~/.delimp_proteog_builds.rds` mirrors `values$proteog_build_jobs`. Save observer is gated against overwriting a non-empty file with an empty list (so a fresh Shiny session doesn't erase prior state before Discover runs). Restore is wrapped in `isolate({})` because reactiveValues access outside a reactive consumer throws.
- **SSH-aware download launchers**: `launch_slims_download()`, `launch_ena_download()`, `poll_download_status()` accept `ssh_config = NULL`. When set, `ssh_exec mkdir`, `scp_upload` the status JSON + script, then `nohup bash ... </dev/null &` over SSH. The `</dev/null` is mandatory — without it the SSH connection waits on the background process's stdin.
- **Deferred build submit for sra/slims modes**: submit handler does NOT call `submit_proteogenomics_build()` immediately; queues in `pending_build_submits` reactiveVal + adds a placeholder row with `download_pending=TRUE` (blue "downloading" badge) to Active Builds. A 15s `reactivePoll` observer watches `download_status.json` and fires the build when state="complete". Skip-if-present: before launching any download, `find`s for `*.fastq.gz` under the target project_dir; if present, treats as local-mode.
- **Per-row Assemble button**: shown in the Active Builds Action column when stage 11 (`assemble`) has `status: "unknown"` and no `job_id`. Opens a modal with the same UniProt source dropdown (None / Download from UniProt / Download from NCBI / Enter path on Hive). When the user picks "Download from UniProt" and the upload completes, the assemble job auto-submits without a second click (gated on `proteog_assemble_target` being non-NULL).
- **UniProt/NCBI download integration**: reuses `search_uniprot_proteomes()`, `download_uniprot_fasta()`, `ncbi_search_assemblies()`, `ncbi_download_proteome()` via proteog-prefixed observers (`proteog_uniprot_state`, `proteog_ncbi_state`). Downloads land in `/quobyte/proteomics-grp/de-limp/databases/uniprot/` on Hive (created via SCP). NCBI also uploads the side-car `_gene_map.tsv` since RefSeq headers don't embed gene symbols.

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
| Spectronaut trailing dots in sample names | Spectronaut appends `.` to labels ending in digits (e.g., "AD12." → "AD12"). `match_samples()` strips with `gsub("\\.$", "", x)`. |
| Spectronaut `PG.UniProtIds` fallback | Some Spectronaut exports lack `PG.ProteinGroups`. Protein column regex includes `UniProtIds` as fallback. Q-value regex includes `Q.Value` variant. |
| Spectronaut Quant3 inflates significance | "Use All MS-Level Quantities" doubles observation count (21v20 → 42v40 in t-test). Detected via `parse_spectronaut_search_settings()`, shown as red "severe" row in settings diff. |
| Spectronaut `Group` not `ProteinGroup` | Candidates.tsv uses `Group` for protein accessions. Regex must include `^Group$`. |
| Spectronaut `Comparison (group1/group2)` format | Comparison column has parenthetical suffix. Regex `^Comparison$` fails — remove `$` anchor. |
| Spectronaut 0-ratio proteins have NaN | Proteins with 0 `# of Ratios` have NaN logFC/Pvalue/Qvalue. `classify_de()` uses `is.finite()`, `assign_hypothesis()` coerces to safe defaults (0 for logFC, 1 for adjP). |
| Spectronaut `AnalyisOverview.txt` typo | Filename may be misspelled by Spectronaut. Regex detector handles both spellings: `analy.?is.?overview`. |
| History tab slow with network CSV | Multiple `activity_log_read()` calls per render cycle. Use `cached_activity_log()` reactive to read once per invalidation. |
| **NEVER use mounted drives for app state** | SMB mounts (`/Volumes/proteomics-grp/`) may be absent, slow, or disappear. All app state files (activity log, cluster usage, lab members) MUST use local paths (`~/.delimp_*`). Cross-user sharing via SSH/SCP sync when connected. |
| **Derived data stays with source data** | Cached/computed data (TIC cache, session.rds, search_info.md) belongs in the raw data or output directory — NOT in a user's home directory. This ensures: (1) any lab member scanning the same directory gets cached results, (2) data is portable (copy dir = copy everything), (3) lifecycle is tied to the dataset (delete data = delete cache). Use SCP for remote directories. Pattern: `.delimp_tic_cache.rds` in raw data dir, `session.rds` in output dir. Only app-level config (activity log, lab members, cluster usage) belongs in `~/.delimp_*`. |
| SSH auto-connect blocks event loop | SSH connection test on startup can take 10-30s. Run via `later::later()` or ensure fast-fail with short timeout. Stale ControlMaster sockets detected with `ssh -O check`. |
| NCBI RefSeq accessions need gene mapping | DIA-NN `Genes` column is empty/accession-only for NCBI FASTA. Use batch E-utilities (`esummary` on protein UIDs) to get gene symbols. Gene map TSV cached alongside FASTA. |
| Contaminant proteins have `Cont_` prefix | After `--fasta` contaminant library, DIA-NN prefixes protein IDs with `Cont_`. Detect via `grepl("^Cont_", Protein.Group)`. Link to NCBI Protein (not UniProt) for these. |
| No-replicates mode skips DE | Groups with <2 replicates get quantification only. `values$fit` remains NULL. Expression Grid, PCA, Signal Distribution still work. Volcano/DE table/GSEA require replicates. |
| SLURM proxy inside Apptainer | SLURM commands not available inside container. Proxy process outside relays commands via temp files. All 9 SLURM command paths covered. |
| SSH file browser Unix socket path | File browser modal creates SSH connections. Same ControlMaster socket path length constraint applies. |
| Docker SSH key permissions | Windows Docker bind mounts lose Unix permissions. `Launch_DE-LIMP_Docker.bat` copies keys to container-internal volume with `chmod 600`. |
| Spectronaut 20+ RunOverview key-value format | Spectronaut 20+ writes RunOverview as 2-column key-value pairs, not wide table. `parse_spectronaut_run_overview()` handles both formats via `ncol()` check. |
| `unlist()` on nested lists causes row mismatch | Nested list elements expand to multiple rows in `data.frame()`. Use `vapply(x, function(v) paste(v, collapse="; "), character(1))` instead. |
| Character matrix subsetting fails on Linux | Empty string rownames cause `mat[protein_ids, ]` subscript errors on Linux (works on macOS). Use numeric indices via `match(protein_ids, rownames(mat))`. |
| Load from HPC needs build-time guard | `conditionalPanel` alone insufficient on HF — button renders briefly before JS hides it. Wrap in `if (!is_hf_space)` in `build_ui()`. |
| Paths with spaces in SLURM scripts | Quote all paths in sbatch: `#SBATCH -o "path"`, `--bind "path":/work`. Launcher uses `q()` helper. |
| Partial retry dependency chain | After retrying failed step 2 tasks, must `scontrol update` step 3's dependency to wait for retry job ID. Otherwise step 3 starts before retries complete. |
| `jsonlite::fromJSON` on `NA_character_` round-trips as `list()` | When `.init_status_json()` serializes `started_at = NA_character_`, jsonlite writes `{}` (empty JSON object). On re-read, `fromJSON(simplifyVector=FALSE)` parses it back as `list()` (length 0). `nzchar(as.character(list()))` returns `logical(0)`, which crashes `if`. **Fix**: use the file-scope `.empty_or_str()` helper in `server_proteog_builder.R` for any value read out of a parsed status.json. |
| Proteogenomics stage with no `job_id` was falsely marked complete | Stage 11 (`assemble`) is non-SLURM (in-process consolidation) — old builds had `status:"unknown"` + `job_id:{}`. The poll loop called `.sacct_state("")` → garbage state → `any_running` stayed FALSE → `current_stage = "complete"`. **Fix**: stages with no SLURM job_id are now treated as still-pending, not complete. |
| `nohup bash ... &` over SSH won't return without `</dev/null` | The remote bash inherits stdin from the SSH session. Even with `nohup`, SSH waits for the inherited stdin to close. The background process runs but `ssh_exec` blocks until timeout. **Fix**: every backgrounded command over SSH needs `< /dev/null` between `bash` and `&`. |
| reactiveValues access outside a reactive consumer crashes the session | `values$proteog_build_jobs` at server-function entry (e.g., restore-from-disk on module init) throws "Can't access reactive value outside of reactive consumer", killing the entire module's observers. **Fix**: wrap all eager `values$*` reads at module init in `isolate({...})`. |
| Persistence observer wrote empty list on startup, clobbering disk | A naive `observe({ saveRDS(values$proteog_build_jobs, ...) })` fires once on startup with the default empty list before any prior state is restored — overwriting the saved file. **Fix**: gate the save: if `length(jobs) == 0 && file.exists(path) && length(readRDS(path)) > 0`, skip. |
| Ensembl URL patterns are species-specific (toplevel vs primary_assembly) | `dna.primary_assembly.fa.gz` only exists for human + mouse. Most other species use `dna.toplevel.fa.gz`. **Fix**: default to `toplevel` for any new reference genome additions to the registry. |
| Ensembl Plants moved to ebi.ac.uk mirror | `ftp.ensemblgenomes.org` returns no response (DNS/connection refused). Use `ftp.ebi.ac.uk/ensemblgenomes/pub/plants/release-<N>/...` instead. The path structure under there mirrors the old layout. |
| NCBI RefSeq accessions lack embedded gene symbols | `XP_*`, `NP_*`, `WP_*` headers have no gene names. DIA-NN's `Genes` column comes out empty/accession-only. **Fix**: `ncbi_download_proteome()` generates a side-car `<basename>_gene_map.tsv`. The proteog NCBI download flow uploads BOTH files to the Hive cache so gene-symbol resolution works at result-load time. |
| Proteog FASTA library catalog is per-user, FASTAs are shared | `~/.delimp_fasta_library/catalog.rds` is local (per CLAUDE.md "never use mounted drives for app state"). The actual FASTAs live on shared Hive storage. **Result**: another lab member's catalog won't auto-populate from your builds. **Fix**: "Discover from Hive" button on the Proteogenomics DBs modal scans `PROTEOG_RNASEQ_ROOT/*/status.json` and registers every completed build that isn't already in their catalog. Idempotent — re-running just updates existing entries. |
| Two DIA-NN containers, only one reads .raw | `/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif` has .NET bundled and reads Thermo `.raw` files. `/quobyte/proteomics-grp/apptainers/diann2.3.0.sif` does NOT have .NET — .raw files silently skipped. Always use the `dia-nn/` version unless you only have `.d`/`.mzML`. |
| DIA-NN binary inside container is `/diann-2.3.0/diann-linux` | Not just `diann` on PATH. Apptainer exec needs the full path: `apptainer exec image.sif /diann-2.3.0/diann-linux ...`. |

Detailed gotchas + HPC paths: @docs/GOTCHAS.md, @docs/HPC_PATHS.md (cascadia-denovo branch extras: ddaPASEF catalogs, Casanovo training, Cascadia model patches).

### Queue Switching
Auto-switches parallel jobs between `genome-center-grp/high` and `publicgrp/low` partitions. Steps 2/4 (array) move to low (preemptible); steps 1/3/5 (assembly) stay on high. See @docs/QUEUE_SWITCHING.md for full logic, known issues, and SLURM state mapping.

## Version History

Current version: **v3.11.0** — defined in `VERSION` file. See [CHANGELOG.md](CHANGELOG.md) for details.

Key decisions: Modularization (v2.3) | XIC Viewer (v2.1) | Phospho Phase 1 (v2.4) | GSEA multi-DB (v2.5) | SSH job submission (v2.5) | Docker backend (v3.0) | MOFA2 (v3.0) | Core Facility (v3.1) | **UI overhaul to page_navbar** (v3.1) | Volcano/CV fixes + Export panel (v3.1.1) | **About tab, community stats, docs overhaul** (v3.2.0) | **Search history, log parser, Claude export enhancements, sacct fixes** (v3.2.1) | **Chromatography QC** (v3.3.0) | **Run Comparator** (v3.4.0) | **Run Comparator enhancements, Search/Analysis History, smart partitions, FASTA library fixes** (v3.5.0) | **Spectronaut parsing fixes, TopN scatter fix, sub-tab help modals** (v3.5.1) | **Unified activity log, cluster monitoring, Spectronaut NaN/rescue fixes** (v3.6.0/v3.6.1) | **Docker launcher, NCBI integration, contaminant analysis, SSH file browser, Load from HPC, no-replicates mode** (v3.7.0) | **Data Completeness visualization, SVG vector export, NCBI gene symbols in DE/GSEA, auto-queue switching fixes** (v3.8.x) | **MaxLFQ + limma alternative pipeline (Moschem 2025), QuantUMS quality filters, On/Off Proteins panel, coverage filter, pipeline descriptor self-describing objects, safe_section export pattern, MANIFEST.txt, R 4.6 / Bioc 3.23 fallback, provenance block (parquet MD5 + sessionInfo)** (v3.9.x) | **Comparator pipeline-aware (descriptor-driven), Methods README pipeline-aware, modal/violin/Methods text pipeline-aware, On/Off panel metadata-alignment fix, Export Complete Analysis consolidated as true superset of Claude exports, FASTA picker modal, On/Off column rename (`detected_g1`/`total_g1`)** (v3.10.x) | **Proteogenomics workflow: RNA-seq → SLURM chain (auto-assemble) → FASTA library catalog → main-page picker, Restore/Discover from Hive for multi-user catalogs, SSH-aware download launchers with deferred-submit on download completion, per-row Assemble button for legacy builds, UniProt+NCBI download integration with auto-submit-after-download** (v3.11.0)

Unreleased (post-v3.10.8): Cascadia de novo integration (feature branch). Audit follow-ups: covariate-name single source of truth (refactor 22 read sites), `<<-` in `add_covariate()` closure refactor, `DIANN_DEFAULTS` constant extraction in `helpers_search.R`, **`detect_organism_db()` silent-fallback refactor** — currently returns `"org.Hs.eg.db"` as a default with no signal to callers; should return `NULL` or `list(db, method)` so callers (PROMPT.md builder, Run Comparator, Explorer prompt, GSEA tab) handle uncertainty explicitly per CLAUDE.md Architectural Rule #2 (caused the v3.10.6 → v3.10.7 "Organism: Human" bug on a Peromyscus dataset). Sweep for similar silent-fallback functions (`coalesce_setting`, default Q-value cutoffs, default mass-acc fallbacks). **Sweep all export bundlers for the `if (!is.null(f)) files_to_zip <- c(...)` and `tryCatch(error = function(e) NULL)` patterns and convert to `safe_section()`** — the v3.10.4 → v3.10.8 hotfix train (5 patches in one day) was almost entirely silent-failure-masquerading-as-success bugs in export paths. Affected files to audit: `R/server_ai.R` (Claude AI export), `R/server_viz.R` (Explorer Claude export), `R/server_phospho.R` (Phospho Claude export), `R/server_mofa.R` (MOFA Claude export), `R/server_comparator.R` (Comparator Claude export). The Complete Analysis bundle in `R/server_session.R` is now clean as of v3.10.8 but the others are still likely-broken in the same ways.
