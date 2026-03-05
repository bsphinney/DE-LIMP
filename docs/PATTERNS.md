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
