# DE-LIMP Detailed Patterns & Gotchas

This file contains detailed patterns that are referenced from CLAUDE.md but too verbose to keep in the main file. Claude loads this on demand when working with the relevant subsystems.

## R Shiny Reactivity

- **DT table row indices** refer to the CURRENTLY DISPLAYED (possibly filtered) data
- **Avoid circular reactivity**: Don't filter data based on a reactive value that the filter updates
- **CRITICAL reactive loop**: If `renderDT` depends on a reactive that uses `values$selection`, AND selecting rows updates `values$selection`, you get an infinite loop. **Solution**: Make `renderDT` build data independently.
- Use `isolate()` to break reactive dependencies when needed
- DT tables with `selection = "multiple"` support Ctrl+Click and Shift+Click
- **CRITICAL `<<-` inside `withProgress`**: `withProgress` uses `eval(substitute(expr), env)` internally. The `<<-` operator does NOT work on lists/scalars inside this eval — it fails with "object not found" because the parent environment chain is broken. **Solution**: Use `new.env(parent = emptyenv())` instead of `list()` for mutable state, and `<-` (not `<<-`) to modify fields. Environments are pass-by-reference, so mutations inside `withProgress` are visible outside. Convert back with `as.list(env)` if needed.

## bslib navset_card_tab Rendering Issues

- **`renderUI`/`uiOutput` inside `navset_card_tab` sub-tabs**: Content renders briefly then disappears. Tried `cancelOutput = TRUE`, `htmlOutput` + `renderText`, `insertUI`/`removeUI` — all fail. bslib clears the DOM for HTML output bindings on tab re-activation.
- **`renderPlot`/`plotOutput` in hidden sub-tabs**: Crashes with `invalid quartz() device size` on macOS. The quartz PNG device gets 0 width/height when the tab isn't visible. Adding `req(input$tab_id == "Tab Name")` does NOT reliably prevent this — the container may still have 0 dimensions even when the tab is "active".
- **`plotlyOutput` works reliably**: Plotly renders via JavaScript in the browser and handles hidden/0-width containers gracefully. Use `renderPlotly`/`plotlyOutput` for any dynamic content inside `navset_card_tab` sub-tabs.
- **First child compression**: The first `plotlyOutput` inside a `navset_card_tab` `nav_panel` can get compressed to a thin strip. Place the primary chart first; secondary/smaller plotly widgets after.
- **Plotly annotation cards are fragile**: `layout(annotations = ...)` with `bgcolor` for card-like summary stats gets compressed, overlaps, and renders inconsistently in bslib sub-tabs. After 15+ attempts, **use ggplot subtitles instead** — a single line of text that's always reliable.
- **Scrollable wrappers prevent compression**: Wrap dense sub-tab content in `div(style = "overflow-y: auto; max-height: calc(100vh - 200px);")` and add `div(style = "min-height: Npx;")` around key widgets.
- **Safe pattern**: `plotlyOutput` > `plotOutput` > `uiOutput` for reliability in bslib sub-tabs.

## Volcano Plot Patterns

- **Y-axis**: Use raw `P.Value` for y-axis (gives classic volcano spread). Do NOT use `adj.P.Val` — it compresses the plot to a flat strip.
- **Threshold line**: Compute as `max(P.Value)` among proteins with `adj.P.Val < 0.05`. This draws the dashed line at the raw P.Value that corresponds to the FDR boundary. BH adjustment is monotonic, so this is exact.
- **Significance coloring**: Base on `adj.P.Val < 0.05` only (not logFC cutoff). logFC vertical lines are visual guides the user can adjust, but don't gate the coloring.
- **Count annotation**: Show "N DE proteins (X up, Y down)" in the info box — users can't easily count overlapping dots.

## Package Installation & Loading

- Install packages BEFORE any `library()` calls (attempting to update loaded packages causes unload errors)
- Use `update = FALSE` in `BiocManager::install()` and `upgrade = "never"` in `remotes::install_*()`
- Set CRAN mirror before installs: `options(repos = c(CRAN = "https://cloud.r-project.org"))`

## Phosphoproteomics Patterns

- **Auto-detection**: `detect_phospho()` scans `Modified.Sequence` for `UniMod:21` on file upload; threshold >100 phosphoprecursors
- **Two input paths**: Path A = DIA-NN 1.9+ `site_matrix_*.parquet` upload (recommended); Path B = parse from `report.parquet` via `extract_phosphosites()`
- **Site extraction (Path B)**: Walk `Modified.Sequence` character-by-character to find `(UniMod:21)`, record preceding residue + position. Aggregate per SiteID x Run via max intensity ("Top 1" method).
- **Imputation**: Tail-based / Perseus-style: `rnorm(n, mean = global_mean - 1.8 * global_sd, sd = 0.3 * global_sd)`. `set.seed(42)` for reproducibility.
- **Normalization options**: none (default), median centering, quantile normalization
- **Independent from protein-level**: Phospho pipeline uses its own `phospho_fit` object, separate from `values$fit`. Both can coexist.
- **Conditional UI**: Sidebar controls and tab content use `conditionalPanel(condition = "output.phospho_detected_flag")` with `outputOptions(suspendWhenHidden = FALSE)`
- **Site positions are peptide-relative** in Path B (no FASTA mapping). Labeled as such in the table.

## GSEA & Organism Detection

- **Organism detection priority**: 1) Suffix-based (`_HUMAN`, `_MOUSE` in Protein.Group), 2) UniProt REST API lookup, 3) Default to human
- **ID mapping fallback**: Try UNIPROT → ENTREZID first, then SYMBOL → ENTREZID
- **ID format handling**: Strips pipe format (`sp|ACC|NAME`), isoform suffixes (`-2`), organism suffixes (`_HUMAN`)
- **Per-ontology caching**: `values$gsea_results_cache` is a named list keyed by ontology (BP/MF/CC/KEGG)
- **Taxonomy mapping**: 12 species mapped from NCBI taxonomy ID to Bioconductor OrgDb name

## SSH Remote Job Submission Patterns

- **`ssh_exec()` is the single chokepoint** — all remote commands flow through one function in `helpers_search.R`
- **Login shell vs non-interactive**: Non-interactive SSH only sources `~/.bashrc`, NOT `~/.bash_profile`. HPC module paths require `bash -l -c` (login shell). But login shell on HIVE takes 10-20s, so use sparingly.
- **Save full SLURM paths at connection test**: `test_ssh_connection()` probes for sbatch once (slow, login shell), saves full path in `values$ssh_sbatch_path`. All subsequent operations use the full path directly (fast, no login shell).
- **`processx::run()` for timeouts**: R's `system2()` has no timeout support. `processx::run(timeout=60)` prevents indefinite SSH hangs. Falls back to `system2()` if processx unavailable.
- **Shell quoting layers**: Commands pass through R string → `shQuote()` → SSH arg → remote sshd → bash. Key rule: **never `shQuote()` glob patterns** — quote only the directory path, leave `/*.d` unquoted.
- **`find` vs `du` for Bruker .d dirs**: `find -maxdepth 2` recurses INTO `.d` directories (~33k internal files). Use `du -sm path/*.d` to treat `.d` as a single unit.
- **SCP for script upload**: More robust than heredoc for writing sbatch scripts to remote.
- **Key-based auth only**: No password storage, no `sshpass`. Uses `StrictHostKeyChecking=accept-new`.
- **R regex gotcha**: `\\s` is NOT valid in base R regex (POSIX ERE). Use `[:space:]` in character classes.
- **Shell line continuations**: Use `paste0(flags, " \\")` to append `\` — NOT `c(flags, " \\")`.
- **SSH ControlMaster multiplexing**: All `ssh_exec()`, `scp_upload()`, `scp_download()` include `-o ControlMaster=auto -o ControlPath=<path> -o ControlPersist=300`. Reuses a single TCP connection across all SSH/SCP calls, eliminating HPC `MaxStartups` throttling. Control socket path must be short (`/tmp/.delimp_<user>_<host>`) — macOS limits Unix domain sockets to 104 bytes.
- **`parse_sbatch_output()` must `trimws()`**: SSH stdout often has trailing `\r` or whitespace. Without trimming, dependency chains like `--dependency=afterok:12345\r` silently fail.
- **`MallocStackLogging` noise**: macOS injects `MallocStackLogging` warnings into stderr. Set `MallocStackLogging = ""` in processx env to suppress.

## Search Output Directory Layout

All DIA-NN search jobs (single, parallel, Docker, local) use this structure:
```
{output_dir}/
├── logs/                          # SLURM .out/.err and local .log files
│   ├── diann_{jobid}.out          # single-job stdout
│   ├── diann_{jobid}.err          # single-job stderr
│   ├── diann_{step}_{jobid}.out   # parallel step stdout (e.g. diann_lib_pred_12345.out)
│   ├── diann_{step}_{jobid}.err   # parallel step stderr
│   └── diann_{name}.log           # local/Docker execution log
├── search_info.md                 # Auto-generated metadata (params, job IDs, file list, logs path)
├── *.sbatch                       # SLURM submission scripts
├── file_list.txt                  # Raw file list (parallel mode)
├── submit_all.sh                  # Launcher script (parallel mode)
├── quant_step2/                   # Per-file quant output (parallel mode)
├── quant_step4/                   # Final quant output (parallel mode)
├── report.parquet                 # Final DIA-NN output
└── *.tsv / *.parquet              # Other DIA-NN outputs
```

- **`logs/` subdir created at submit time**: `mkdir -p` in all paths (SSH, Docker, local)
- **`search_info.md`** includes `**Log files**: \`{output_dir}/logs/\`` for reference
- **Backwards compatibility**: Job monitoring checks `logs/` first, falls back to `output_dir/` root for pre-existing jobs
- **`generate_search_info()`** in `helpers_search.R`: Generates `search_info.md` content with all search metadata

## Parallel DIA-NN Search (5-Step)

- **Architecture**: 5 SLURM steps with dependency chaining: (1) library prediction → (2) first-pass quant array → (3) assemble → (4) second-pass quant array → (5) final assembly
- **3-connection approach**: Batch all remote operations into 3 SSH/SCP calls to avoid HPC `MaxStartups` throttling:
  1. One SSH: `mkdir -p` for output dir + `logs/` + quant subdirs
  2. One SCP: upload all sbatch scripts + file_list.txt + `submit_all.sh` launcher
  3. One SSH: `bash submit_all.sh` — chains all `sbatch` submissions with `--dependency=afterok:$PREV_ID`
- **Launcher script**: Shell script that submits each step, parses job IDs with `grep -oP "[0-9]+$"`, outputs `STEP1:id` through `STEP5:id`. R-side parses with `regexec("^STEP([1-5]):(.+)$", line)`.
- **Resource allocation**: Step 1 (lib prediction) uses `libpred_cpus=16, libpred_mem=64G` — lighter than assembly steps (64 CPUs, 512G) to avoid QOS CPU limits. Steps 2/4 (per-file quant) use `cpus_per_file=16, mem_per_file=32G` as array jobs.
- **`generate_parallel_scripts()`** in `helpers_search.R`: Builds all 5 sbatch scripts + launcher. Parameters include per-step resource controls, max simultaneous array tasks, search params.
- **Cluster resource indicator**: `check_cluster_resources()` queries sacctmgr (group CPU limit), squeue (group CPUs in use), sinfo (partition idle/total). Polled every 60s when SSH connected. Traffic-light UI: green (<50%), yellow (50-80%), red (>80%).

## Docker/HPC Dual Backend Patterns

- **Backend detection at startup** (`app.R`): Checks `Sys.which("docker")` + `docker info` for Docker; `Sys.which("sbatch")` + `Sys.which("ssh")` for HPC.
- **Backend selector UI**: When only one backend, use hidden `radioButtons` in `div(style="display:none;")` — must be a real Shiny input for `conditionalPanel` to work.
- **Shared flag builder**: `build_diann_flags()` called by both HPC (sbatch) and Docker (`build_docker_command()`) paths.
- **DIA-NN license compliance**: DO NOT bundle DIA-NN binary. Users build their own image via `build_diann_docker.sh`.
- **Job queue persistence**: Saved to `~/.delimp_job_queue.rds`. `observeEvent(values$diann_jobs, {...}, ignoreInit = TRUE)` prevents race condition.
- **Job entry schema**: `list(job_id, backend="docker"|"hpc", name, status, output_dir, n_files, ...)`. Old jobs missing `backend` field default to `"hpc"`.
- **Loading results**: Must pass `format = "parquet"` to `limpa::readDIANN()` — default is `"tsv"`.
- **`return()` inside `withProgress`**: Exits `withProgress` not the calling function. Restructure as flat `tryCatch`.
- **SSH output encoding**: All `ssh_exec()`/`scp_download()`/`scp_upload()` sanitize with `iconv(..., sub="")`.

## Instrument Metadata Extraction

- **Auto-extraction at file scan time**: `parse_raw_file_metadata()` dispatches to `parse_timstof_metadata()` or `parse_thermo_metadata()` based on file extension. Only ONE file is parsed (first in list).
- **timsTOF (.d)**: Reads `analysis.tdf` SQLite — `GlobalMetadata` table for instrument model/serial/m/z range/1/K0 range, `Frames` table for RT range/MS1+MS2 spectra counts/dia-PASEF detection/cycle time. Also tries `submethods/*.method` XML for LC method name.
- **Thermo (.raw)**: Uses `rawrr::readFileHeader()` for instrument model, serial, mass range, time range, scan counts, mass resolution, software version. `rawrr` is optional — graceful fallback if not installed.
- **SSH remote timsTOF**: SCP downloads `analysis.tdf` (~few MB) from first `.d` dir to temp, parses locally.
- **SSH remote Thermo (.raw)**: Uses ThermoRawFileParser (`-m=0 -f=0`) on remote system to extract JSON metadata, then SCP downloads the JSON and parses locally via `parse_thermorawfileparser_json()`. ThermoRawFileParser probed via `which` + login shell fallback. Does NOT extract LC method/gradient (GitHub issue #86). Available as BioContainer on HPC systems.
- **timsTOF LC info**: Extracted from `HyStarMetadata.xml` (UTF-8 XML with LC system name, method, runtime, serial). DIA windows from `.m/diaSettings.diasqlite`.
- **Auto-sets m/z range**: Updates `min_pr_mz`/`max_pr_mz` UI inputs from actual instrument range, fixing the 300-1200 default issue.
- **Stored in `values$instrument_metadata`**: Named list, persisted in session save/load.
- **Surfaces in**: mass accuracy hint, raw file summary badge, methodology text, Claude export (CSV + PROMPT.md inline block).

## Chromatography QC — TIC Extraction & Diagnostics

- **Purpose**: Extract TIC traces from raw files at scan time (before search) to flag dead injections, carryover, RT drift, and loading anomalies before committing to hours-long DIA-NN searches.
- **User-triggered**: "Extract TIC" button appears after file scan. Users can skip for verified good files. Not automatic.
- **timsTOF only**: Reads MS1 frames (`MsMsType = 0`) from `analysis.tdf` SQLite. Auto-detects intensity column (`SummedIntensities` → `AccumulatedIntensity` → `MaxIntensity` → any `*ntensit*`) for older TDF compatibility. Also handles missing `MsMsType` column (falls back to `ScanMode` or no filter).
- **SSH mode**: SCP downloads each `analysis.tdf` (~5-30 MB) to temp, extracts locally, deletes temp. Same pattern as instrument metadata extraction.
- **Reactive values**: `values$tic_traces` (named list of data.frames), `values$tic_metrics` (data.frame with per-run metrics, diagnostics, file sizes).
- **QC tab visibility**: QC tab shown when `values$raw_data` OR `values$tic_traces` is not NULL — allows viewing TIC QC before search runs.
- **Per-run metrics** (`compute_tic_metrics()`): Total AUC (trapezoid), peak RT, gradient boundaries (10% threshold), gradient width, baseline ratio, late signal ratio, asymmetry. Smoothed with 2% moving average.
- **Shape similarity** (`compute_shape_similarity()`): Interpolates all traces onto 500-point common RT grid, computes Pearson r vs median trace. Single-file returns r=1.0.
- **Diagnostics** (`diagnose_run()`): MAD-based outlier detection (needs 3+ runs for cohort checks):
  - Shape deviation: r < 0.90 = fail, r < 0.95 = warn
  - RT shift: peak RT > 3 MAD from median
  - Loading anomaly (AUC): >3x/<0.3x median = fail; >2x/<0.5x = warn
  - File size outlier: >3x/<0.3x median = fail; >2x/<0.5x = warn
  - Late elution: >15% signal in last 20% of gradient
  - Elevated baseline: >10% of peak intensity
  - Narrow gradient: <70% of median width
- **Three plot views** (all plotlyOutput for bslib safety):
  - Faceted: Per-run panels (ncol=4) with median trace overlay (blue dashed). Status-colored lines (green/yellow/red).
  - Overlay: All runs normalized 0-1 on one axis. Quick outlier shape detection.
  - Metrics: Horizontal bar chart of AUC per run with DT table below (Size, AUC, Peak RT, Gradient Width, Baseline Ratio, Late Signal, Shape r, Flags).
- **Session persistence**: `tic_traces` and `tic_metrics` saved/loaded with session `.rds`.
- **Auto-trigger bug**: `observeEvent(list(btn, trigger))` fires when button first renders (NULL→0 change). Fixed by using separate `reactiveVal` trigger that both Extract and Re-extract buttons increment.

## XIC Viewer Patterns

- **Arrow masks dplyr**: `arrow::select()` masks `dplyr::select()` — always use `dplyr::select()` explicitly
- **Avoid rlang in Arrow context**: Use base R `df[df$col %in% vals, ]` and `names(df)[...] <- "new"`
- **Precursor map**: Built from in-memory `values$raw_data`, not file I/O
- **Dual format**: `detect_xic_format()` auto-detects v1 (wide) vs v2 (long) DIA-NN XIC formats
- **Platform guard**: All XIC UI/logic wrapped in `if (!is_hf_space)`

## About Dropdown & Community Stats

- **About** is a `nav_menu()` dropdown with two sub-tabs:
  - **Community** (`value = "about_tab"`): Stats, sparklines, discussions, links — narrow layout (`max-width: 900px`)
  - **Analysis History** (`value = "analysis_history_tab"`): Full-width page with project filter, DT table, expandable rows
- **Version source**: `VERSION` file → `app.R` reads at startup → `values$app_version` → all modules use it
- **Community stats**: `stats/community_stats.json` loaded at startup → `values$community_stats`
  - Generated daily by `.github/workflows/track-stats.yml` (GitHub Actions)
  - Contains: stars, forks, 14-day clone/view counts, daily trend arrays
  - Graceful NULL fallback — app works fine without stats file
- **Stats cards**: 4 value boxes (Stars, Forks, Visitors, Clones) rendered via `uiOutput`
- **Trend sparklines**: 2 plotly line charts (views + clones). Uses `plotlyOutput` per bslib safety pattern.
- **Discussions section**: Recent discussions from GitHub Discussions (title, category, author, comment count, link). Fetched via GraphQL in the workflow.
- **Stats freshness**: Displays `updated_at` timestamp from JSON

## Search History

- **Search history CSV**: `search_history_path()` in `helpers_search.R` — shared volume first, local `~/.delimp_search_history.csv` fallback. Same `filelock` pattern as analysis history.
- **Schema (26 fields)**: timestamp, completed_at, user, search_name, backend, search_mode, parallel, n_files, fasta_files, fasta_seq_count, normalization, enzyme, mass_acc_mode, mass_acc, mass_acc_ms1, scan_window, mbr, extra_cli_flags, output_dir, job_id, status, duration_min, speclib_cached, imported_from_log, app_version, notes
- **Record points**: Job submission (`server_search.R` — single shared block after all 4 backends converge at `values$diann_jobs` assignment). Status update in job monitoring observer on completion/failure.
- **`record_search()`**: Append-only, same pattern as `record_analysis()`.
- **`update_search_status()`**: Read-modify-write for status/completed_at/duration_min by `output_dir` match. Not append-only (file rewrite) but file is small.
- **Cross-reference**: `output_dir` joins search history ↔ analysis history. Both DT tables show a link icon button when a matching row exists in the other table. Click navigates to that tab via `nav_select("main_tabs", ...)`.
- **Import Settings button**: Reads CSV row fields, builds settings list, applies to search UI inputs (mass_acc, enzyme, mode, normalization, extra flags). Stores in `values$diann_search_settings` with `imported_from_log = TRUE`.
- **View Log button**: Reads `search_info.md` from `output_dir` (SSH first, local fallback) and displays in modal.
- **Expandable detail rows**: Same JS `row.child()` pattern as analysis history. Shows enzyme, mass acc, scan window, MBR, normalization, extra flags, output_dir, job_id, speclib cached, imported from log, notes.

## Analysis History & Projects

- **Analysis history CSV**: `analysis_history_path()` in `helpers_search.R` — shared volume first, local `~/.delimp_analysis_history.csv` fallback. Append-only with `filelock`.
- **Record points**: Pipeline completion (`server_data.R`), auto-load from search (`server_search.R`), session file load (`server_session.R`).
- **n_proteins pitfall**: `nrow(raw_data$E)` = precursors (~40k). Correct: `length(unique(raw_data$genes$Protein.Group))` for protein groups (~3k). `y_protein$E` rows are protein groups (post-pipeline only).
- **Projects JSON**: `projects_path()` — shared volume first, local `~/.delimp_projects.json` fallback. Schema: `{projects: {name: {created_at, description, entries: [output_dir, ...]}}}`
- **Join key**: `output_dir` links history CSV rows to projects and live job queue statuses.
- **DT expandable rows**: Hidden column 0 contains detail HTML. JS callback on `tbody tr td:not(:last-child)` toggles `row.child()`. Detail div forces `color:#2d3748` to avoid white-on-white with DT row selection highlighting.
- **Action buttons**: Use `event.stopPropagation()` to prevent row expansion when clicking Info/Load/Assign buttons. Pass `output_dir` directly in JS `Shiny.setInputValue()` (not row index, which breaks with project filtering).
- **Project summary cards**: Shown above table when project filter is active. Computed from filtered history rows.
- **Assign modal**: `selectizeInput(create = TRUE)` for existing/new project names. Uses `reactiveVal(assign_od)` to store the target output_dir (not a hidden HTML input, which Shiny can't read).

## Core Facility Mode

- **Activation**: `DELIMP_CORE_DIR` env var pointing to directory containing `staff.yml`
- **Config object** (`cf_config`): `list(staff, qc, db_path, reports_dir, state_dir, template_qmd)`
- **SQLite database** (`delimp.db`): 4 tables — `searches`, `qc_runs`, `reports`, `templates`
- **WAL mode**: All SQLite connections use `PRAGMA journal_mode=WAL`
- **Staff selector**: Auto-fills SSH host, username, key path, SLURM account/partition from YAML
- **Report generation**: `generate_report_impl()` → save state `.rds` → render Quarto `.qmd` → move HTML → record in SQLite
- **Quarto rendering**: `quarto::quarto_render()` with `output_file` as filename only (not full path). Move HTML after render.
- **y_protein is EList**: Extract `$E` for the expression matrix (not `as.matrix()` directly)
- **Project field**: `selectizeInput` with `create = TRUE` for autocomplete + new entry
- **Test data**: `seed_test_db.R` generates synthetic SQLite data

## DIA-NN Parallel Search — Critical Flags & Gotchas

Reference: [DIA-NN Discussion #1414](https://github.com/vdemichev/DiaNN/discussions/1414) (Vadim Demichev)

### 5-Step Pipeline

| Step | Purpose | Key Flags |
|------|---------|-----------|
| 1 | Library prediction (from FASTA) | `--fasta-search --predictor --gen-spec-lib --out-lib step1.speclib` |
| 2 | First-pass per-file quant (array) | `--lib step1.predicted.speclib --gen-spec-lib --quant-ori-names` |
| 3 | Empirical library assembly | `--lib step1.predicted.speclib --use-quant --quant-ori-names --gen-spec-lib --rt-profiling --out-lib empirical.parquet` |
| 4 | Final per-file quant (array) | `--lib empirical.parquet --no-ifs-removal --quant-ori-names` |
| 5 | Cross-run report | `--lib empirical.parquet --use-quant --quant-ori-names --matrices` |

### Critical Rules (per Vadim)

1. **`--quant-ori-names` at ALL steps** — Preserves original filenames in `.quant` files. Without it, quant file naming can mismatch between steps when container bind mount paths differ.
2. **Fixed mass accuracy when using `--use-quant`** — Auto-optimisation + quant reuse produces different results from the original analysis. Always set `--mass-acc`, `--mass-acc-ms1`, and `--window` explicitly. Our `generate_parallel_scripts()` forces `mass_acc_mode = "manual"` for this reason.
3. **`--fasta-search` and `--predictor` in Step 1 ONLY** — Including these in Steps 2-5 causes DIA-NN to re-digest the FASTA from scratch instead of using the predicted/empirical library. This was the original FASTA re-digest bug (commit `d2b9bc6`).
4. **Empirical library is `.parquet`, not `.speclib`** — DIA-NN 2.0+ saves empirical libraries (from `--gen-spec-lib --out-lib`) in Apache Parquet format. Predicted libraries (Step 1) remain `.predicted.speclib`. Both formats are loadable with `--lib`.
5. **Verify quant files before assembly steps** — Array jobs can fail silently (preemption, OOM, CPU mismatch). Steps 3 and 5 include a bash verification block that checks all expected `.quant` files exist before running DIA-NN. With `--quant-ori-names`, quant files are named `BASENAME.quant` (e.g., `sample.raw` → `sample.quant`).
7. **Backup Step 2 quant files before Step 3** — Step 3 uses `--use-quant --temp quant_step2`, which OVERWRITES the Step 2 quant files (same names with `--quant-ori-names`). A `cp -r quant_step2 quant_step2_orig` backup runs before DIA-NN in Step 3. Smart resume from Step 3 restores from this backup. Without it, a Step 3 failure would require re-running Step 2.
6. **`--no-ifs-removal`** — Used in Step 4 (second-pass quant) for high-precision quantification with the empirical library.

### Common Failure Modes

| Failure | Cause | Fix |
|---------|-------|-----|
| "Cannot load spectral library" in Step 4 | Library file referenced as `.speclib` but DIA-NN saved it as `.parquet` | Use `empirical.parquet` in `--lib` and `--out-lib` |
| "Illegal instruction (core dumped)" | DIA-NN binary compiled for newer CPU (e.g., AVX-512) than the assigned HPC node | Add `--constraint` to sbatch or rebuild container for wider CPU compatibility |
| Step 3 re-digests FASTA from scratch | `--fasta-search` or `--predictor` leaked into Steps 2-5 | `step_flags` filtering in `generate_parallel_scripts()` removes these |
| Quant files not found by assembly step | Array tasks failed/were preempted, or naming mismatch from different bind mount paths | Quant verification block + `--quant-ori-names` on all steps |
| DIA-NN warning about mass accuracy + quant reuse | Auto mass accuracy mode with `--use-quant` | Force manual mode: `--mass-acc N --mass-acc-ms1 N --window N` |

### Implementation

- **Script generator**: `generate_parallel_scripts()` in `R/helpers_search.R` (~lines 1467-1815)
- **Step flag filtering**: Lines 1541-1551 strip step-specific flags from `base_flags`
- **Mass accuracy enforcement**: `parallel_sp$mass_acc_mode <- "manual"` (line ~1540)
- **Quant verification**: `quant_verify_block()` helper generates bash check code
- **Resume launcher**: `generate_resume_launcher()` for resubmitting from failed step

## Run Comparator

- **Module**: `R/server_comparator.R` contains all parsers, helpers, and server logic (no separate helpers file)
- **Three modes**: DE-LIMP vs DE-LIMP (Mode A), vs Spectronaut (Mode B), vs FragPipe (Mode C1/C2)
- **Protein ID normalization**: `normalize_protein_id()` strips `sp|...|`, isoform suffixes, semicolon groups to bare UniProt accessions
- **Sample matching**: `match_samples()` normalizes filenames (strip extensions, FragPipe column prefixes) and requires all samples matched before comparison proceeds
- **4 diagnostic layers**: Settings Diff, Protein Universe, Quantification, DE Concordance
- **3x3 concordance matrix**: Improvement over spec's 2x2 — classifies each protein as Up/Down/NS in each run for more nuanced concordance view
- **8-rule hypothesis engine**: `assign_hypothesis()` returns hypothesis text, confidence (High/Medium/Low), and category. Tool-aware: rules include context-specific notes for Spectronaut/FragPipe differences. Rule 3 names Quant3 when detected. Rule 8 checks low Spectronaut n_ratios (underpowered significance).
- **Hypothesis categories**: Direction reversal, Normalization offset, Variance estimation, Missing values, Peptide count, FC magnitude, Borderline, Low n_ratios
- **Protein Universe Venn**: Plotly shape-based Venn diagram with proportional circles. Fixed center positions with minimum separation to prevent label overlap at high Jaccard. Filterable DT table below with All/Shared/A-only/B-only buttons + CSV export.
- **Gene symbol validation**: `is_gene_symbol()` filters out accessions. DIA-NN `Genes` column often has accessions like `A0A075B6K5;P80748`. Real gene symbols from `bitr()` or `sp|ACC|GENE` format parsing.
- **DE Concordance explanatory text**: Blue info banner above matrix explaining Up/Down/NS terminology and concordance percentage.
- **All plotly rendering**: No `renderPlot` in sub-tabs (avoids quartz crash). Venn diagram via plotly shapes. Correlation heatmap replaced with plotly bar chart. Bias density as plotly histogram.
- **AI Analysis tab — static HTML only**: `uiOutput`/`renderUI` vanishes inside `navset_card_tab` sub-tabs. All content is static HTML. Gemini narrative injected via `shinyjs::html("comparator_gemini_container", ...)`. MOFA2 uses static `plotlyOutput` + `DTOutput` with `req()` guards. Button hidden via `shinyjs::hide()` after MOFA2 completes.
- **Optional MOFA2 decomposition**: Treats Run A and Run B as two views, decomposes joint variance. **CRITICAL**: MOFA2 requires all views to have identical column names. Use matched sample pairs with common labels (`Sample_1`, `Sample_2`, ...). Requires >= 4 matched pairs.
- **Spectronaut ZIP upload**: `parse_spectronaut_zip()` handles Report Exporter ZIP containing Pivot report, Candidates.tsv, ConditionSetup.tsv, RunSummaries, ExperimentSetupOverview, AnalysisLog.txt. Returns de_stats, intensities, n_peptides, n_ratios, run_qc, library_info, search_settings, manifest.
- **Spectronaut parsers**: `parse_spectronaut_search_settings()` extracts TopN, normalization, Quant3 from ExperimentSetupOverview. `parse_spectronaut_log()` extracts library size, version from AnalysisLog.txt. Both feed into enriched settings diff.
- **Spectronaut protein ID fallback**: Column regex includes `PG.UniProtIds` fallback. Q-value regex includes `Q.Value` variant. `match_samples()` strips trailing dots (Spectronaut appends `.` to labels ending in digits).
- **TopN quantification analysis**: Settings Diff shows amber banner when >30% proteins limited by TopN cap. Quantification sub-tab shows TopN Effect scatter (n_peptides vs |quant diff|, log scale, LOESS).
- **Quant3 (Use All MS-Level Quantities) detection**: Red "severe" row in settings diff + danger banner when enabled. Effectively doubles observation count in t-test (21v20 → 42v40). Rule 3 in hypothesis engine specifically names Quant3 as cause.
- **Per-Sample QC (Mode B)**: Horizontal bar chart from RunSummaries with outlier detection (<60% median protein count). DT table with Size, Precursors, ProteinGroups, IDRate, CV metrics.
- **Settings diff match styles**: 5 levels — `match` (transparent), `differs` (amber), `unknown` (grey), `structural_difference` (blue #cce5ff), `severe` (red #f8d7da). Blue for inherent tool differences (e.g., TopN vs DPC-Quant), red for problematic settings (e.g., Quant3 enabled).
- **Gemini prompt is data-driven**: `build_gemini_comparator_prompt()` computes pct_limited, mean_pep, n_sig, pct_sig, median_lfc from actual data. Includes Quant3 conditional paragraph. `diff_rows` filter includes `structural_difference` and `severe`.
- **Claude ZIP export**: Settings diff CSV, protein universe CSV, DE results combined CSV, discordant proteins CSV with hypotheses (including n_ratios_B), claude_prompt.md, comparison_context.md. Mode B adds spectronaut_run_qc.csv and spectronaut_library_info.csv.
- **Session persistence**: All comparator state (results, parsed runs, mode, Gemini narrative, MOFA object) saved/loaded with session .rds
- **Progressive reveal**: Tab always visible (not gated behind pipeline completion)
- **No ComplexUpset dependency**: Spec called for it, but Venn diagram is more intuitive for 2-set comparisons
- **Internal reactiveVals**: `comp_run_a`, `comp_run_b`, `comp_results` are local `reactiveVal()` inside the module, synced to `values$comparator_*` for session persistence
- **Spectronaut 0-ratio rescue stats**: Proteins with 0 `# of Ratios` have NaN DE stats (untestable in Spectronaut). `compute_de_concordance()` computes `rescue_stats` (n_zero_ratio, n_rescued_sig, n_low_ratio) and returns them in the result list. Rule 0 in hypothesis engine ("Untestable in Spectronaut") is imputation-aware with three variants (None/enabled/unknown).
- **Contrast mismatch detection**: `compute_de_concordance()` fuzzy-matches Spectronaut conditions against DE-LIMP contrasts via normalized alphanumeric comparison. Amber warning div shown when no contrasts overlap.
- **NaN-safe classify_de()**: Uses `is.finite()` guards for both logFC and adjP. Returns "NS" for any non-finite value.
- **NaN-safe assign_hypothesis()**: Coerces non-finite logFC to 0, non-finite adjP to 1 at function entry, before any conditional logic.
- **Instrument context in AI prompts**: `build_gemini_comparator_prompt()` and `build_claude_comparator_prompt()` accept optional `instrument_meta` parameter. Emits brief "INSTRUMENT: model, LC, SPD, gradient" line when available.

## NCBI Integration Patterns

- **NCBI Proteome Download**: Uses NCBI Datasets API (`datasets.ncbi.nlm.nih.gov/v2alpha/genome/taxon/{organism}/dataset_report`) to find reference proteomes. Downloads protein FASTA via assembly accession.
- **Gene symbol mapping**: NCBI RefSeq accessions (XP_, NP_, WP_) lack embedded gene names unlike UniProt `sp|ACC|GENE` format. Batch E-utilities lookup: `esearch` for protein UIDs, then `esummary` to extract gene name from `DocumentSummarySet`. Results cached as TSV alongside FASTA.
- **Gene map TSV**: Stored as `{fasta_basename}_gene_map.tsv` with columns `accession` and `gene_symbol`. Applied during data loading to populate the `Genes` column in DIA-NN output.
- **Docker gene map download**: Docker users may lack direct E-utilities access. Gene map TSV auto-downloaded from HPC via SSH (`scp_download()`) when not found locally.
- **NCBI protein links**: Proteins with RefSeq accessions link to `https://www.ncbi.nlm.nih.gov/protein/{accession}` instead of UniProt. `Cont_` prefixed proteins detected separately.

## Contaminant Analysis Patterns

- **Detection**: Contaminant proteins identified by `Cont_` prefix in `Protein.Group` column (added by DIA-NN when using `--fasta` contaminant library).
- **Summary cards**: Count, percentage of total proteins, median intensity ratio (contaminant median / endogenous median), keratin count (grep for `KRT|Keratin|keratin` in gene/protein names).
- **Per-sample stacked bar**: Each sample bar shows endogenous (blue) and contaminant (orange) protein counts. Useful for spotting samples with unusual contamination.
- **Top contaminants table**: Sorted by median intensity across samples. Keratin flag column (yellow highlight). Downloadable as CSV.
- **Contaminant heatmap**: Top 20 contaminants by median intensity. Rows = proteins, columns = samples. Color scale: white (low) to orange/red (high).
- **Expression Grid highlighting**: `Cont_` rows get pink background CSS. Helps users visually identify contamination in the full expression matrix.
- **Signal Distribution overlay**: Checkbox adds orange contaminant protein distribution overlaid on the main signal distribution plot. Shows where contaminants fall in the intensity range.

## SSH File Browser Patterns

- **Modal-based**: Opens as a Shiny modal with path display, breadcrumbs, and file list.
- **SSH directory listing**: `ssh_exec()` with `ls -la` to list directory contents. Parsed into data frame with name, type, size.
- **Color coding**: Folders (blue folder icon), data files `.d`/`.raw`/`.parquet`/`.fasta` (green file icon), other (grey).
- **Breadcrumb navigation**: Path split into clickable segments. Each segment navigates to that directory level.
- **Context-specific filtering**: Raw data browser shows `.d`/`.raw`/`.mzML`/`.wiff` directories/files. FASTA browser shows `.fasta`/`.fa`. Results browser shows `.parquet`.
- **Configurable roots**: `DELIMP_EXTRA_ROOTS` env var (comma-separated paths) adds institution-specific root directories to the browser navigation.
- **Performance**: Uses specific subdirectory roots rather than scanning from `/`. Avoids `ls` on directories with >10,000 entries.

## Data Explorer Patterns (v3.7)

- **Location**: Data Overview > Data Explorer subtab in `server_viz.R`
- **No DE required**: Works with just `values$y_protein` — designed for no-replicates workflows
- **Two panels** in scrollable accordion:
  - **Abundance Profiles (Quartile Analysis)**: Splits proteins into Q1-Q4 by average intensity. Plotly heatmap shows top 10 per quartile, colored by per-SAMPLE quartile (not average). Reveals proteins that shift rank across samples. "Variable Proteins" DT table lists proteins with quartile range >= 2.
  - **Sample-Sample Scatter**: Pick two samples, XY scatter with identity line. Points colored by distance from diagonal (grey=similar, red=different). Outliers (>4-fold) labeled with gene names. Contaminants as orange triangles. Stats subtitle: Pearson r, N proteins, N outliers. Built with ggplot + ggplotly.
- **Gene label helper**: `explorer_gene_label()` resolves protein IDs via y_protein$genes > sp|ACC|GENE parsing > truncation fallback
- **Quartile computation**: Per-sample quartile assignment uses each sample's own intensity distribution (not the average). This means a protein can be Q1 in one sample and Q3 in another — the inconsistency IS the signal.
- **Contaminant exclusion**: Both panels have independent "Exclude contaminants" checkboxes (default TRUE)

## No-Replicates Mode (v3.7)

- **Detection**: After group assignment, if any group has <2 samples, DE analysis is skipped.
- **What works**: Data import, normalization (DPC-CN), protein quantification (DPC-Quant), Expression Grid, Signal Distribution, PCA, Dataset Summary, Data Explorer, Contaminant Analysis.
- **What's skipped**: limma model fitting (`values$fit` stays NULL), volcano plot, DE results table, GSEA, CV Analysis. User sees informational message.
- **Expression Grid fallback**: When `values$fit` is NULL, grid renders from `y_protein$E` with a dummy `P.Value = 1` column so the table structure is consistent.
- **PCA without fit**: PCA tab checks for `y_protein` (quantified expression) not `values$fit`. Visible as soon as quantification completes.

## Docker + SSH Deployment Patterns (v3.7)

- **SSH auto-connect**: On startup, `app.R` checks for SSH keys in Docker-mounted volume (`/app/data/ssh/`) and `~/.ssh/`. If found, triggers `test_ssh_connection()` via `later::later()` to avoid blocking the event loop.
- **Stale socket detection**: `ssh -O check` probes ControlMaster socket before auto-connect. Dead sockets (from previous app instances) are removed.
- **Environment badge**: CSS-styled badge in navbar. Detection order: `SPACE_ID` env (HF orange) > `/.singularity.d/` or `APPTAINER_*` env (HPC green) > Docker socket/env (Docker red) > fallback (Local blue).
- **SLURM proxy for Apptainer**: Shell script outside container listens for command files in shared storage, executes SLURM commands, writes results. All 9 SLURM paths: `sbatch`, `squeue`, `scancel`, `sacct`, `sinfo`, `sacctmgr`, `scontrol`, `srun`, `sbatch --test-only`.
- **Container detection**: `is_container()` checks for Apptainer/Docker markers. When true, skips BiocManager package validation (no internet in container = validation hangs).
- **Home directory quota warning**: Startup check via `df -h ~` or `quota` command. Warns if >80% used. Common on HPC where home dirs are 5-10 GB.
- **Auto-adjust search CPUs**: When submitting a DIA-NN search, CPUs automatically capped to per-user SLURM limit (from `check_cluster_resources()`). Prevents job rejection due to QOS limits.
- **Windows Docker SSH key handling**: `Launch_DE-LIMP_Docker.bat` copies SSH keys from `C:\Users\{username}\.ssh\` into Docker-accessible volume (`./data/ssh/`). Sets permissions via `docker exec chmod 600`. Multiple Windows users on same PC each get their own key detected.
