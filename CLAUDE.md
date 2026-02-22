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
| `app.R` | Thin orchestrator (~270 lines) — package loading, backend detection, reactive values init, calls R/ modules |
| `R/ui.R` | `build_ui(is_hf_space, search_enabled, docker_available, hpc_available, local_sbatch)` — full UI definition (~900 lines) |
| `R/server_data.R` | Data upload, example load, group assignment, pipeline execution (~576 lines) |
| `R/server_de.R` | Volcano, DE table, heatmap, consistent DE, selection sync (~498 lines) |
| `R/server_qc.R` | QC trends, diagnostic plots, p-value distribution (~828 lines) |
| `R/server_viz.R` | Expression grid, signal distribution, AI summary (~381 lines) |
| `R/server_gsea.R` | GSEA analysis, multi-DB (BP/MF/CC/KEGG), organism detection, caching (~400 lines) |
| `R/server_ai.R` | AI Summary (all contrasts), Data Chat, Gemini integration (~280 lines) |
| `R/server_xic.R` | XIC viewer, mobilogram, alignment (~861 lines) |
| `R/server_session.R` | Info modals, save/load session, reproducibility, methodology (~533 lines) |
| `R/helpers.R` | `get_diann_stats_r()`, `cal_z_score()`, `detect_organism_db()` |
| `R/helpers_ai.R` | `list_google_models()`, `upload_csv_to_gemini()`, `ask_gemini_*()` |
| `R/helpers_xic.R` | `detect_xic_format()`, `load_xic_for_protein()`, `reshape_xic_for_plotting()` |
| `R/helpers_phospho.R` | `detect_phospho()`, `parse_phospho_positions()`, `extract_phosphosites()` (~210 lines) |
| `R/server_phospho.R` | Phospho site-level DE, volcano, site table, residue dist, completeness QC (~650 lines) |
| `R/helpers_mofa.R` | `parse_rds_for_mofa()`, `parse_matrix_file()`, `compute_sample_overlap()`, `generate_mofa_code()` (~200 lines) |
| `R/server_mofa.R` | MOFA2 multi-view integration: view management, training, variance heatmap, weights, scores, Factor-DE correlation (~700 lines) |
| `R/helpers_search.R` | `ssh_exec()`, `scp_download/upload()`, `build_diann_flags()`, `build_docker_command()`, `check_docker_*()`, `recover_slurm_jobs()`, `recover_docker_jobs()`, `generate_sbatch_script()`, UniProt proteome search (~1050 lines) |
| `R/server_search.R` | Docker/HPC dual backend, SSH connection, DIA-NN search config, file browsing, job queue (submit/monitor/cancel/recover/load), SCP result download (~1750 lines) |
| `build_diann_docker.sh` | User-facing script to build DIA-NN Docker image locally (license notice, downloads from GitHub, builds x86_64 image) |
| `Dockerfile` | Thin HF/HPC image — `FROM brettphinney/delimp-base:v3.0` + COPY app code |
| `Dockerfile.base` | Pre-built base image with all system deps, R/CRAN/Bioc/MOFA2+basilisk (~5 GB) |
| `HPC_DEPLOYMENT.md` | Guide for HPC cluster deployment with Apptainer/Singularity |
| `USER_GUIDE.md` | End-user documentation |
| `README_GITHUB.md` | **SOURCE** for GitHub README (edit this!) |
| `README_HF.md` | **SOURCE** for HF README (edit this!) |
| `README.md` | Generated - content differs between remotes |
| `Citations/` | Reference PDFs (LIMMA, FDR, DIA-NN papers) |
| `docs/index.html` | GitHub Pages site |
| `R/helpers_facility.R` | `cf_get_bracketing_qc()`, `cf_get_rolling_stats()`, `compute_rolling_cv()`, WAL mode, staff/lab/instrument/LC names (~340 lines) |
| `R/server_facility.R` | Core facility: report generation, job history, QC dashboard, templates, staff auto-config (~660 lines) |
| `report_template.qmd` | Quarto HTML report template — metadata, QC bracket, volcanos, DE stats, top proteins, normalization |
| `seed_test_db.R` | Synthetic test data generator (20 searches, 30 QC runs, 3 reports) |
| `.github/workflows/sync-to-hf.yml` | Auto-sync GitHub to HF |

## App Architecture

### Structure (modular, ~7,100 lines total)
```
app.R (270 lines):      Package loading, auto-install, backend detection (Docker/HPC/Core Facility), reactive values, module calls
R/ui.R:                 build_ui() — CSS/JS, sidebar, all 12+ tab nav_panels (some conditional on core facility mode)
R/server_*.R (11 files): Server modules, each receives (input, output, session, values, ...)
R/helpers*.R (5 files): Pure utility functions (no Shiny reactivity)
```

### Module calling pattern (app.R server function)
```r
server <- function(input, output, session) {
  values <- reactiveValues(...)
  add_to_log <- function(action_name, code_lines) { ... }

  server_data(input, output, session, values, add_to_log, is_hf_space)
  server_de(input, output, session, values, add_to_log)
  server_qc(input, output, session, values)
  server_viz(input, output, session, values, add_to_log, is_hf_space)
  server_gsea(input, output, session, values, add_to_log)
  server_ai(input, output, session, values)
  server_xic(input, output, session, values, is_hf_space)
  server_phospho(input, output, session, values, add_to_log)
  server_search(input, output, session, values, add_to_log,
                search_enabled, docker_available, docker_config, hpc_available, local_sbatch)
  server_mofa(input, output, session, values, add_to_log)
  server_facility(input, output, session, values, add_to_log,
                  is_core_facility, cf_config, search_enabled)
  server_session(input, output, session, values, add_to_log)
}
```

### Main Tab Structure
- **Data Overview** - 6 sub-tabs: **Assign Groups & Run** (+ phospho detection banner), Signal Distribution, Dataset Summary, Group QC Summary, Expression Grid, **AI Summary**
- **QC Trends** - 4 sub-tabs: Precursors, Proteins, MS1 Signal, Stats Table
- **QC Plots** - 5 sub-tabs: Normalization Diagnostic, DPC Fit, MDS Plot, Group Distribution, P-value Distribution
- **DE Dashboard** - Volcano plot, results table, violin plots, interactive comparison selector
- **Consistent DE** - 2 sub-tabs: High-Consistency Table (ranked by %CV), CV Distribution (histogram by group)
- **Reproducibility** - Code Log + Methodology sub-tabs
- **Phosphoproteomics** - Site-level DE: Phospho Volcano, Site Table, Residue Distribution, QC Completeness (conditional on phospho detection)
- **Gene Set Enrichment** - Ontology selector (BP/MF/CC/KEGG), Dot Plot, Enrichment Map, Ridgeplot, Results Table, per-ontology caching
- **Multi-View Integration** - MOFA2: dynamic view cards (2-6), file upload (RDS/CSV/TSV/Parquet), sample matching, training, 5 results tabs (Variance Explained, Factor Weights, Sample Scores, Top Features, Factor-DE Correlation)
- **New Search** - DIA-NN search: dual backend (Docker local + HPC SSH/SLURM), file selection, FASTA database (UniProt download + upload), search settings (standard/phospho), resource controls, job queue (submit/monitor/cancel/recover/load), auto-load results (conditional on search backend availability)
- **Data Chat** - AI-powered analysis (Google Gemini API)
- **Search DB** *(core facility only)* - Full-width job history table, 6 filters (text/lab/status/staff/instrument/LC method), Load Results + Generate Report buttons
- **Instrument QC** *(core facility only)* - QC trend plots (proteins, precursors, TIC), instrument filter, date range, runs table
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
- `values$phospho_detected` - list(detected, n_phospho, n_total, pct_phospho, is_enriched)
- `values$phospho_site_matrix` - numeric matrix: sites x samples (log2)
- `values$phospho_site_info` - data.frame: SiteID, Protein.Group, Genes, Residue, Position
- `values$phospho_fit` - limma MArrayLM object for site-level DE
- `values$phospho_site_matrix_filtered` - after missingness filtering + imputation
- `values$phospho_input_mode` - "site_matrix" or "parsed_report"
- `values$gsea_results_cache` - Named list of cached GSEA results per ontology (BP=res, MF=res, etc.)
- `values$gsea_last_contrast` - Contrast used for cached GSEA results
- `values$gsea_last_org_db` - OrgDb used for cached GSEA results
- `values$ssh_sbatch_path` - Full path to sbatch binary on remote HPC (cached from test connection)
- `values$diann_jobs` - List of DIA-NN job entries (job_id, backend, name, status, output_dir, n_files, log_content, is_ssh, loaded, etc.)
- `values$mofa_view_configs` - List of view config objects (id, name, type, source, matrix, status)
- `values$mofa_views` - Named list of loaded matrices (view_name -> features x samples matrix)
- `values$mofa_view_fits` - Optional limma fit objects from RDS imports (for Factor-DE correlation)
- `values$mofa_object` - Trained MOFA2 model object
- `values$mofa_factors` - Extracted factor matrix (samples x factors)
- `values$mofa_weights` - List of weight matrices per view
- `values$mofa_variance_explained` - Variance explained per view per factor

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

### Core Facility Mode
- **Activation**: `DELIMP_CORE_DIR` env var pointing to a directory containing `staff.yml` (e.g., `/srv/delimp` or `.core_facility_test/`)
- **Detection**: `is_core_facility <- dir.exists(core_facility_config_dir) && file.exists(file.path(core_facility_config_dir, "staff.yml"))`
- **Config object** (`cf_config`): `list(staff, qc, db_path, reports_dir, state_dir, template_qmd)` — loaded once at startup, passed to `build_ui()` and `server_facility()`
- **SQLite database** (`delimp.db`): 4 tables — `searches` (job tracking with lab/instrument/project/LC method), `qc_runs` (HeLa QC metrics per instrument), `reports` (generated report metadata + QC bracket IDs), `templates` (saved search presets/group templates)
- **WAL mode**: All SQLite connections use `PRAGMA journal_mode=WAL` for concurrent read performance
- **Staff selector**: Choosing a staff member auto-fills SSH host, username, key path, SLURM account/partition from YAML config. Auto-tests SSH connection on selection.
- **Report generation flow**: `generate_report_impl()` → save state `.rds` → render Quarto `.qmd` → move HTML to reports dir → attempt HF upload → record in SQLite reports table → show success UI with download button
- **Quarto rendering**: `quarto::quarto_render()` with `output_file` as filename only (not full path — quarto rejects paths). Rendered HTML moved from template dir to reports dir. Requires `quarto` R package + CLI (`brew install quarto`).
- **State object**: Full analysis snapshot (fit, y_protein, metadata, design, phospho, GSEA, MOFA results + metadata_extra with title/lab/instrument/lc_method/project/analyst)
- **y_protein is EList**: When reading state in report template, extract `$E` for the expression matrix (not `as.matrix()` directly)
- **Job history**: `job_history_data()` reactive queries SQLite with 6 filter params. DT table in dedicated Search DB tab (moved out of cramped Panel 3 card).
- **Project field**: `selectizeInput` with `create = TRUE` for autocomplete + new entry. Dropdown filtered by staff + lab selections via `bindEvent(observe(...))`.
- **Conditional UI**: All core facility tabs (`Search DB`, `Instrument QC`) and sidebar sections wrapped in `if (is_core_facility)`. Non-facility users see no difference.
- **Test data**: `seed_test_db.R` generates synthetic SQLite data (20 searches, 30 QC runs across 2 instruments, 3 reports, 2 templates)

## Development Workflow

### Running Locally
```r
# In VS Code R terminal (recommended) — directory-based:
shiny::runApp('/Users/brettphinney/Documents/claude/', port=3838, launch.browser=TRUE)

# From command line:
Rscript -e "shiny::runApp('/Users/brettphinney/Documents/claude/', port=3838)" &
```

- **DO NOT** use `source()` - it doesn't work properly in VS Code
- `app.R` explicitly sources all `R/*.R` files, so it works with both `runApp('.')` and `runApp('app.R')`
- Shiny apps don't hot-reload - must restart after every code change
- **Upload limit**: 5 GB (`options(shiny.maxRequestSize)` in app.R line 145)
- Stop: `pkill -f "shiny::runApp"` or `Ctrl+C`
- Check if running: `lsof -i :3838`


## Deployment

### Three Platforms
1. **GitHub** (`origin` remote) - Source code, releases
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
git add app.R R/           # add orchestrator + modules
git commit -m "Description"
git push origin main       # HF syncs automatically
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

### Docker Base Image (CRITICAL for HF deploys)
HF has a ~30-minute build timeout. All heavy dependencies live in a **pre-built base image** on Docker Hub (`brettphinney/delimp-base:v3.0`). The HF `Dockerfile` is just `FROM` that image + `COPY` app files — builds in under 1 minute.

- **Base image**: `brettphinney/delimp-base:v3.0` on Docker Hub (public)
- **Contents**: rocker/shiny:4.5.0 + all system deps + CRAN + Bioconductor + MOFA2/basilisk
- **Source**: `Dockerfile.base` in repo
- **Size**: ~5 GB
- **Build machine**: Windows box at `~/DE-LIMP` (has Docker Desktop)

**Code-only changes** (most pushes): Just `git push origin main` — HF rebuilds instantly.

**Adding new R packages** (rebuild base image required):
```bash
# 1. Edit Dockerfile.base — add package to appropriate RUN line
# 2. On Windows box:
cd ~/DE-LIMP
git pull origin main
docker build --platform linux/amd64 -t brettphinney/delimp-base:v3.0 -f Dockerfile.base .
docker push brettphinney/delimp-base:v3.0
docker tag brettphinney/delimp-base:v3.0 brettphinney/delimp-base:latest
docker push brettphinney/delimp-base:latest
# 3. Also update Dockerfile.search if it exists (keeps parity)
# 4. Commit Dockerfile.base, app.R, R/ modules together
# 5. git push origin main — HF picks up new base image on next build
```

**If HF build fails with "not found"**: Check Docker Hub repo is set to **Public** at https://hub.docker.com/r/brettphinney/delimp-base. Factory reboot the Space at HF settings if needed.

**Docker Hub credentials**: Username `brettphinney`, log in with `docker login -u brettphinney` on the Windows box.

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

### Phosphoproteomics Patterns
- **Auto-detection**: `detect_phospho()` scans `Modified.Sequence` for `UniMod:21` on file upload; threshold >100 phosphoprecursors
- **Two input paths**: Path A = DIA-NN 1.9+ `site_matrix_*.parquet` upload (recommended); Path B = parse from `report.parquet` via `extract_phosphosites()`
- **Site extraction (Path B)**: Walk `Modified.Sequence` character-by-character to find `(UniMod:21)`, record preceding residue + position. Aggregate per SiteID x Run via max intensity ("Top 1" method).
- **Imputation**: Tail-based / Perseus-style: `rnorm(n, mean = global_mean - 1.8 * global_sd, sd = 0.3 * global_sd)`. `set.seed(42)` for reproducibility.
- **Normalization options**: none (default), median centering, quantile normalization
- **Independent from protein-level**: Phospho pipeline uses its own `phospho_fit` object, separate from `values$fit`. Both can coexist.
- **Conditional UI**: Sidebar controls and tab content use `conditionalPanel(condition = "output.phospho_detected_flag")` with `outputOptions(suspendWhenHidden = FALSE)`
- **Site positions are peptide-relative** in Path B (no FASTA mapping). Labeled as such in the table.
- **Spec documents**: `PHOSPHO_TAB_SPEC.md` (full 3-phase spec with citations), `PHOSPHO_TAB_PHASE1_COMPACT.md` (Phase 1 implementation guide)

### GSEA & Organism Detection
- **Organism detection priority**: 1) Suffix-based (`_HUMAN`, `_MOUSE` in Protein.Group), 2) UniProt REST API lookup (`rest.uniprot.org/uniprotkb/{accession}` → taxonomy ID), 3) Default to human
- **ID mapping fallback**: Try UNIPROT → ENTREZID first, then SYMBOL → ENTREZID
- **ID format handling**: Strips pipe format (`sp|ACC|NAME`), isoform suffixes (`-2`), organism suffixes (`_HUMAN`)
- **Per-ontology caching**: `values$gsea_results_cache` is a named list keyed by ontology (BP/MF/CC/KEGG)
- **Taxonomy mapping**: 12 species mapped from NCBI taxonomy ID to Bioconductor OrgDb name

### SSH Remote Job Submission Patterns
- **`ssh_exec()` is the single chokepoint** — all remote commands flow through one function in `helpers_search.R`
- **Login shell vs non-interactive**: Non-interactive SSH (`ssh host "cmd"`) only sources `~/.bashrc`, NOT `~/.bash_profile`. HPC module paths (spack, conda) require `bash -l -c` (login shell). But login shell on HIVE takes 10-20s (conda init + spack loading), so use it sparingly.
- **Save full SLURM paths at connection test**: `test_ssh_connection()` probes for sbatch once (slow, login shell), saves the full path (e.g., `/cvmfs/.../slurm/bin/sbatch`) in `values$ssh_sbatch_path`. All subsequent operations use the full path directly (fast, no login shell).
- **`processx::run()` for timeouts**: R's `system2()` has no timeout support. `processx::run(timeout=60)` prevents indefinite SSH hangs. Falls back to `system2()` if processx unavailable.
- **Shell quoting layers**: Commands pass through R string → `shQuote()` → SSH arg → remote sshd → bash. Each layer can mangle quotes, globs, and escapes. Key rule: **never `shQuote()` glob patterns** — quote only the directory path, leave `/*.d` unquoted for bash to expand.
- **`find` vs `du` for Bruker .d dirs**: `find -maxdepth 2` recurses INTO `.d` directories (~33k internal files each). Use `du -sm path/*.d` to treat `.d` as a single unit.
- **SCP for script upload**: More robust than heredoc (`cat <<'EOF'`) for writing sbatch scripts to remote — avoids nested quoting issues.
- **FASTA upload**: After UniProt download, check remote existence with `test -f` before SCP upload to avoid overwriting.
- **Connection mode toggle**: `radioButtons("search_connection_mode")` — local users (sbatch on PATH) unaffected; SSH users explicitly opt in.
- **Key-based auth only**: No password storage in memory, no `sshpass` dependency. Uses `StrictHostKeyChecking=accept-new` to avoid interactive prompt on first connect.
- **R regex gotcha**: `\\s` is NOT valid in base R regex (POSIX ERE). Use `[:space:]` in character classes — `[,;[:space:]]+` not `[,;\\s]+`.
- **Shell line continuations**: When building multi-line shell commands from R vectors, use `paste0(flags, " \\")` to append `\` to each flag line — NOT `c(flags, " \\")` which puts `\` on its own line.

### Docker/HPC Dual Backend Patterns
- **Backend detection at startup** (`app.R`): Checks `Sys.which("docker")` + `docker info` for Docker; `Sys.which("sbatch")` + `Sys.which("ssh")` for HPC. Sets `search_enabled`, `docker_available`, `hpc_available` flags passed to UI and server modules.
- **Backend selector UI**: When both backends available, `radioButtons("search_backend", ...)` with Docker/HPC choices. When only one, hidden `radioButtons` in `div(style="display:none;")` — must be a real Shiny input for `conditionalPanel` to work (raw `tags$input(type="hidden")` is NOT registered by Shiny JS).
- **Shared flag builder**: `build_diann_flags()` extracts DIA-NN CLI flag generation from `generate_sbatch_script()`. Called by both HPC (sbatch) and Docker (`build_docker_command()`) paths. Eliminates flag duplication.
- **Docker execution**: `docker run -d --platform linux/amd64 --cpus=N --memory=Xg` returns container ID. Monitoring via `check_docker_container_status()` → `docker inspect` + `docker logs --tail 30`.
- **DIA-NN license compliance**: DO NOT bundle DIA-NN binary. Users build their own image via `build_diann_docker.sh`. UI footer has license attribution + citation link.
- **Job queue persistence**: Saved to `~/.delimp_job_queue.rds`. `observeEvent(values$diann_jobs, {...}, ignoreInit = TRUE)` prevents race condition where save observer fires before load observer restores data.
- **Job entry schema**: `list(job_id, backend="docker"|"hpc", name, status, output_dir, n_files, search_mode, search_settings, auto_load, log_content, submitted_at, completed_at, loaded, is_ssh, container_id)`. Old jobs missing `backend` field default to `"hpc"`.
- **Recover Jobs**: Queries `sacct` for DIA-NN jobs, then enriches with `scontrol show job` (StdOut path), log tailing, file count extraction. Updates existing entries instead of skipping. Docker recovery via `docker ps -a --filter name=delimp_`.
- **SLURM log discovery** (3 strategies): 1) `scontrol show job <id>` → `sed -n 's/.*StdOut=//p'` (most reliable), 2) `sacct SubmitLine` → parse script path, 3) `timeout 10 find <output_base> -maxdepth 4` (targeted, with timeout).
- **File count extraction**: DIA-NN log format is `[HH:MM] N files will be processed`. Use `grep -m1 'files will be processed'` then `regexpr("[0-9]+(?=\\s+files will be processed)", ..., perl=TRUE)` — NOT `^[0-9]+` which matches the `[HH:MM]` timestamp.
- **SSH output encoding**: All `ssh_exec()`/`scp_download()`/`scp_upload()` sanitize output with `iconv(stdout, from="", to="UTF-8", sub="")` to prevent `nchar(msg): invalid multibyte string` crashes from ANSI codes, MOTD banners, or binary data in SSH output.
- **Loading results**: Must pass `format = "parquet"` to `limpa::readDIANN()` — default is `"tsv"` which uses `fread()` (returns data.table), causing `DT[, extra.columns]` syntax error. Wrap in `suppressMessages(suppressWarnings(...))` to prevent crashes from arrow messages with invalid encoding.
- **`return()` inside `withProgress`**: `return()` exits `withProgress` (not the calling function) because `withProgress` uses `eval(expr, envir=env)`. Never rely on `return()` inside `withProgress` to exit an observer handler — restructure as flat `tryCatch` instead.

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
| SSH "sbatch not found on remote PATH" | Non-interactive SSH doesn't load login profile. Fixed: `test_ssh_connection()` uses `bash -l -c` + common path probing |
| SSH submit hangs indefinitely | `bash -l` on some HPCs takes 10-20s (conda/spack init). Fixed: save full sbatch path from test, use directly without login shell |
| SSH file scan returns thousands of spurious files | Bruker `.d` dirs contain ~33k internal files. Fixed: use `du -sm *.d` instead of `find` |
| SSH glob patterns not expanding | `shQuote()` wrapping glob prevents bash expansion. Fixed: quote only directory path, leave `/*.d` unquoted |
| UniProt proteome search returns "No proteomes found" | `proteome_type` removed from UniProt REST API `fields`. Fixed: removed from query params |
| sbatch script "command not found" errors | Line continuation `\` was on its own line. Fixed: `paste0(flags, " \\")` appends backslash directly |
| `nchar(msg): invalid multibyte string` crash | SSH/SCP output contains ANSI codes or non-UTF-8 bytes. Fixed: `iconv(..., sub="")` in `ssh_exec`/`scp_download`/`scp_upload`. Also wrap `readDIANN` in `suppressMessages(suppressWarnings(...))`. |
| Shiny hidden input not registered by JS | `tags$input(type="hidden")` doesn't register with Shiny. Use `div(style="display:none;", radioButtons(...))` for `conditionalPanel` to see the input. |
| Job queue wiped on app restart | `observe()` for save fires on initial empty `reactiveValues`. Fixed: `observeEvent(..., ignoreInit=TRUE)`. |
| `readDIANN` data.table column error | `limpa::readDIANN()` defaults to `format="tsv"`. Must pass `format="parquet"` for .parquet files or `fread()` returns data.table causing `DT[,extra.columns]` error. |
| `return()` inside `withProgress` doesn't exit handler | `withProgress` uses `eval()` — `return()` exits `withProgress`, not the enclosing function. Restructure with flat `tryCatch` + early returns outside `withProgress`. |
| Recovered job file count wrong (e.g. 83 instead of 6) | DIA-NN log lines start with `[HH:MM]`. Regex `^[0-9]+` matches the timestamp. Use lookahead: `[0-9]+(?=\\s+files will be processed)`. |
| Recover button appears to hang | `find` on large shared filesystems (e.g. `/quobyte/`) is very slow. Fixed: use `timeout 10 find <output_base> -maxdepth 4` targeted search. |
| Quarto `output_file` path error | `quarto::quarto_render()` rejects full paths for `output_file`. Pass filename only, then `file.rename()` to target dir. |
| SQLite `Parameter N does not have length 1` | `NULL` has length 0; SQLite params must be length 1. Use `NA_character_` instead of `NULL` for missing values. |
| Generate Report button does nothing | Silent `req(values$fit)` blocks without feedback. Fixed: explicit `showNotification()` warning when no DE results loaded. |
| Quarto template `colSums` error on `y_protein` | `y_protein` is limma `EList`, not a matrix. Extract `$E` for the expression matrix: `if ("E" %in% names(y)) y$E else as.matrix(y)`. |
| DT table invisible in bslib card_body | `DTOutput` defaults to `fill = TRUE`, collapses to zero height in `card_body` with `overflow-y: auto`. Solution: move to standalone tab or set `fill = FALSE`. |

## Version History

Current version: **v3.0.0** (2026-02-20). See [CHANGELOG.md](CHANGELOG.md) for detailed release history.

### Key Architecture Decisions (for context)
- **Modularization** (v2.3): Split 5,139-line monolith into `app.R` orchestrator + 12 `R/` module files. Each server module receives `(input, output, session, values, ...)`. `app.R` explicitly sources `R/` for compatibility with both `runApp('.')` and `runApp('app.R')`.
- **Volcano → Table filtering** (v2.2): `output$de_table` filters by `values$plot_selected_proteins`. Row selection observer uses `isolate()` to index into filtered data. `req(length > 0)` prevents reactive loop when selection resets to NULL.
- **Contextual Help** (v2.2): 15 info modal `?` buttons across all major tabs. See "Info modal pattern" in Key Patterns above.
- **XIC Viewer** (v2.1): Fragment-level chromatogram inspection with 3 display modes (facet by sample, facet by fragment, intensity alignment). Supports DIA-NN v1/v2 formats and ion mobility. Local/HPC only (hidden on HF).
- **MS2 Intensity Alignment** (v2.1): Spectronaut-style stacked bar chart. Inconsistency detection via `xic_alignment_data()` reactive (deviation > mean + 2*SD flags samples).
- **Four-way selector sync** (v2.1): `contrast_selector`, `contrast_selector_signal`, `contrast_selector_grid`, `contrast_selector_pvalue` all sync bidirectionally.
- **P-value Distribution Diagnostic** (v2.1): Automated health assessment with color-coded guidance banners (healthy/warning/danger/info patterns).
- **Assign Groups** (v2.1): Moved from modal to permanent Data Overview sub-tab for better layout.
- **Phosphoproteomics Phase 1** (v2.4): Site-level DE via limma with two input paths (site matrix upload or parsed from report). Auto-detection on file upload. Independent from protein-level pipeline. No new packages required.
- **GSEA Expansion** (v2.5): Multi-database enrichment (BP/MF/CC/KEGG) with per-ontology caching. UniProt REST API organism detection when protein IDs lack species suffixes. Robust ID extraction handles pipe-separated, isoform, and organism suffix formats.
- **AI Summary All-Contrasts** (v2.5): Loops over all `colnames(values$fit$contrasts)`, computes cross-comparison biomarkers (significant in >=2 contrasts), scales token budget by contrast count (30/20/10 top proteins).
- **SSH Remote Job Submission** (v2.5): Run DE-LIMP on local Mac, submit DIA-NN searches to remote HPC via SSH. Non-blocking job queue (submit multiple, continue using app). `ssh_exec()` single chokepoint, full SLURM path caching after connection test, `processx::run()` timeouts, SCP for file transfer. No new R packages (uses system `ssh`/`scp`).
- **Phosphoproteomics Search Mode** (v2.5): DIA-NN search preset for phospho analysis — auto-configures STY modification (UniMod:21), max 3 variable mods, 2 missed cleavages, `--phospho-output` and `--report-lib-info` flags.
- **Docker Local Backend** (v3.0): "One UI, two engines" — same search config UI for both Docker (local) and HPC (SSH/SLURM) backends. Backend detection at startup (`docker info` / `Sys.which("sbatch")`). Shared `build_diann_flags()` eliminates flag duplication. DIA-NN Docker image built separately by users via `build_diann_docker.sh` (license compliance — DIA-NN is proprietary, free for academic use, cannot be redistributed).
- **Job Queue Recovery** (v2.5): "Recover" button queries `sacct`/`scontrol`/`docker ps -a` to repopulate lost job queue. 3-strategy log discovery for HPC (scontrol → sacct SubmitLine → find). Updates existing entries with fresh data.
- **MOFA2 Multi-View Integration** (v3.0): Standalone tab for unsupervised integration of 2-6 data views using MOFA2. Dynamic view cards with add/remove, smart RDS parser (DE-LIMP session, limma objects, matrices, data frames), CSV/TSV/Parquet matrix upload, phospho tab integration, sample matching with overlap stats. Training with configurable factors/convergence/scaling. 5 results tabs: variance explained heatmap, factor weights browser (plotly), sample scores scatter, top features table (DT), Factor-DE correlation. Session save/load, methodology text, reproducibility logging. HF resource limits (max 10 factors, fast convergence).
- **Core Facility Mode** (v3.1, branch `claude/heuristic-wu`): Optional mode activated by `DELIMP_CORE_DIR` env var. SQLite database tracking all DIA-NN searches with metadata (lab, instrument, project, LC method, analyst). Staff YAML config auto-fills SSH/SLURM fields. Instrument QC dashboard with protein/precursor/TIC trend plots and ±2SD control lines. Search DB tab with 6-filter job history table. Quarto report generation producing standalone HTML with metadata, QC bracket comparison, volcano plots, DE stats, top proteins, normalization diagnostics. Report state saved as `.rds` for re-loading. Template system for saving/loading search presets. WAL mode for concurrent SQLite access.

## Current TODO

### Phosphoproteomics — Phase 2 (Kinase Activity & Motifs)
- [ ] **KSEA integration** (`KSEAapp` CRAN package): Infer upstream kinase activity from phosphosite fold-changes using PhosphoSitePlus + NetworKIN database. Horizontal bar plot of kinase z-scores. New package: `KSEAapp`
- [ ] **Sequence logo / Motif analysis** (`ggseqlogo` CRAN package): Extract ±7 flanking residues around significant phosphosites, display as sequence logos (up-regulated vs down-regulated). Requires FASTA upload for protein-relative positions. New package: `ggseqlogo`
- [ ] **Kinase Activity tab** in phospho results navset: Run KSEA button, bar plot, results table
- [ ] **Motif Analysis tab** in phospho results navset: Logos for up/down regulated sites
- [ ] Dockerfile: Add `KSEAapp`, `ggseqlogo` to CRAN install list

### Phosphoproteomics — Phase 3 (Advanced)
- [ ] **Protein-level abundance correction**: Subtract protein logFC from phosphosite logFC to isolate phosphorylation stoichiometry changes (requires total proteome + phospho-enriched from same samples)
- [ ] **PhosR integration** (Bioconductor): RUVphospho normalization using stably phosphorylated sites (SPSs), kinase-substrate scoring, signalome construction. New package: `PhosR`
- [ ] **AI context for phospho**: Append phosphosite DE results and KSEA kinase activities to Gemini chat context when phospho analysis is active
- [ ] **Phospho-specific FASTA upload**: Map peptide-relative positions to protein-relative positions for accurate site IDs and motif extraction

### DIA-NN Local Search (Docker) — IMPLEMENTED
**Spec**: `DIANN_SEARCH_INTEGRATION_SPEC.md` — full design for unified Local Docker + HPC backends ("one UI, two engines")
- [x] **Docker backend detection**: `app.R` detects Docker daemon + DIA-NN image at startup
- [x] **Local compute backend**: Docker option runs DIA-NN via `docker run -d` with resource limits (CPU/RAM sliders)
- [x] **Unified UI**: Backend selector (Docker/HPC), shared search config, backend-specific resource controls
- [x] **`build_diann_docker.sh`**: User builds own DIA-NN Docker image (license compliance)
- [x] **Job queue dual backend**: Docker/HPC badges, backend-specific monitoring/cancel/refresh
- [x] **Job recovery**: Recover button queries sacct (HPC) and docker ps (Docker) to repopulate lost queue
- [x] **Load results**: SCP download for SSH jobs, direct filesystem for Docker/local, auto-load on completion
- [ ] **End-to-end Docker testing**: Test full Docker submit → monitor → auto-load flow with real data

### MOFA2 Multi-View Integration — IMPLEMENTED
**Spec**: `MOFA2_INTEGRATION_SPEC.md`
- [x] Dynamic view cards (2-6 views, add/remove)
- [x] Smart RDS parser (DE-LIMP session, limma EList, MArrayLM, matrix, named list, data.frame)
- [x] CSV/TSV/Parquet matrix upload with auto-log2 detection
- [x] Phospho tab integration as data source
- [x] Sample matching with overlap statistics
- [x] MOFA training (configurable factors, convergence, scaling, seed)
- [x] 5 results visualization tabs
- [x] Factor-DE correlation (links MOFA factors to DE results from any view)
- [x] Session save/load, methodology text, reproducibility logging
- [x] HF Spaces resource limits
- [ ] **MEFISTO integration**: Temporal/spatial MOFA for time-course experiments
- [ ] **Factor annotation**: Link factors to GO terms based on top weights
- [ ] **DIA-NN report processing**: Process raw DIA-NN .parquet as MOFA view via existing pipeline
- [ ] **Dockerfile**: Add MOFA2 + basilisk to Docker image

### Core Facility Mode — Next Steps
- [ ] **QC run ingestion**: Auto-record QC metrics when loading HeLa digest report.parquet (detect by filename pattern or checkbox)
- [ ] **Report template polish**: Add GSEA section, MOFA variance explained, configurable logo/header
- [ ] **Report comparison**: Side-by-side QC bracket + DE summary for two reports on same instrument
- [ ] **HF state upload/download**: Upload `.rds` state to HF Spaces for shareable live links
- [ ] **Template application on search submit**: Auto-apply saved search preset when submitting new DIA-NN job
- [ ] **Audit log**: Track who generated which report, when, with what parameters
- [ ] **Multi-instrument QC alerts**: Flag instruments where protein count drops below rolling mean - 2*SD
- [ ] **End-to-end testing**: Test full flow with real DIA-NN search → QC ingest → report generation

### General
- [ ] Grid View: Open violin plot on protein click with bar plot toggle
- [ ] Publication-quality plot exports (SVG/PNG/TIFF with size controls)
- [ ] Sample correlation heatmap (QC Plots tab)
- [ ] Venn diagram of significant proteins across comparisons
- [ ] Sample CV distribution plots
- [ ] Protein numbers bar plot per sample
- [ ] Absence/presence table for on/off proteins
