# DE-LIMP Detailed Patterns & Gotchas

This file contains detailed patterns that are referenced from CLAUDE.md but too verbose to keep in the main file. Claude loads this on demand when working with the relevant subsystems.

## R Shiny Reactivity

- **DT table row indices** refer to the CURRENTLY DISPLAYED (possibly filtered) data
- **Avoid circular reactivity**: Don't filter data based on a reactive value that the filter updates
- **CRITICAL reactive loop**: If `renderDT` depends on a reactive that uses `values$selection`, AND selecting rows updates `values$selection`, you get an infinite loop. **Solution**: Make `renderDT` build data independently.
- Use `isolate()` to break reactive dependencies when needed
- DT tables with `selection = "multiple"` support Ctrl+Click and Shift+Click

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
