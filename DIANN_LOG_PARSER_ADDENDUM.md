# RUN COMPARATOR — DIA-NN LOG PARSER ADDENDUM

> **Adds to:** RUN_COMPARATOR_SPEC_V2.md  
> **Motivation:** Real-world comparison (AD proteomics dataset, March 2026) showed that
> the ZIP export's `settings_diff.csv` reported "unknown" for DIA-NN version, search
> parameters, and library provenance — the exact fields needed to diagnose why two
> DE-LIMP sessions produced 94 vs 7 significant proteins from the same raw files.
> Analysis of the actual DIA-NN log files revealed the root cause: different pipeline
> architectures (all-in-one `--fasta-search` vs multi-step) and `--proteoforms` /
> `--pg-level 0` settings expanding the Run A library to 6M precursors vs 4.6M in
> Run B. The log file is the only place this information lives.

---

## What to Parse

DIA-NN writes a complete record of every search parameter to its log file. The
`report_log.txt` (or `.out` file from a SLURM job) contains the full command line
on line 1, followed by confirmation lines ("Min precursor m/z set to 300", etc.)
and summary statistics per sample. This is the authoritative source for search
provenance.

**Key fields that explain peptide count discrepancies between runs:**

| Field | Flag | Why it matters |
|-------|------|----------------|
| `--pg-level` | 0 = isoform, 1 = gene, 2 = protein | Directly controls protein grouping granularity |
| `--proteoforms` | present/absent | Enables proteoform-aware FDR; massively expands search space |
| `--max-fr-mz` | numeric | Fragment ion ceiling; affects which ions score each precursor |
| `--min-fr-mz` | numeric | Fragment ion floor |
| `--mass-acc` | ppm (MS2) | If different, affects precursor matching sensitivity |
| `--mass-acc-ms1` | ppm (MS1) | As above for MS1 |
| `--window` | scans | RT extraction window; affects co-elution interference |
| `--reanalyse` | present/absent | Two-pass search; improves sensitivity |
| `--fasta` | filename | Library source — different FASTA = different protein universe |
| `--out-lib` | filename | Whether an empirical library was generated and used |
| `--lib` | filename | If a pre-built library was used (vs `--fasta-search`) |
| Precursors in library | from log | Direct measure of search space size |
| Pipeline step | from SLURM job name | Whether log is step 1/5 or the full search |

---

## New UI Element — DIA-NN Log Upload (Mode A only)

Add an **optional** collapsible section to the Run Comparator config card,
shown only in `delimp_delimp` mode. This is optional — the comparison runs
without it, but the settings diff will have fewer "unknown" entries when provided.

```r
# Inside the Mode A config panel in ui_comparator.R
conditionalPanel(
  "input.comparator_mode == 'delimp_delimp'",
  div(
    actionLink("toggle_diann_logs", 
               tagList(icon("chevron-right", id = "diann_log_chevron"),
                       " Attach DIA-NN log files (optional — improves settings diff)")),
    conditionalPanel(
      "input.toggle_diann_logs % 2 == 1",
      div(class = "mt-2",
        div(class = "alert alert-info p-2",
          icon("info-circle"), " Upload ", tags$code("report_log.txt"),
          " or the SLURM ", tags$code(".out"), " file from each DIA-NN search.",
          " Only the first line (command) and summary statistics are read."
        ),
        layout_columns(col_widths = c(6, 6),
          fileInput("comparator_diann_log_a", "Run A — DIA-NN log",
                    accept = c(".txt", ".log", ".out")),
          fileInput("comparator_diann_log_b", "Run B — DIA-NN log",
                    accept = c(".txt", ".log", ".out"))
        ),
        uiOutput("comparator_diann_log_status")   # shows parsed summary once loaded
      )
    )
  )
)
```

**Log status UI** — displayed immediately after upload so user can confirm the
right file was uploaded:

```r
output$comparator_diann_log_status <- renderUI({
  logs <- list()
  if (!is.null(input$comparator_diann_log_a))
    logs$a <- parse_diann_log(input$comparator_diann_log_a$datapath)
  if (!is.null(input$comparator_diann_log_b))
    logs$b <- parse_diann_log(input$comparator_diann_log_b$datapath)
  if (length(logs) == 0) return(NULL)

  items <- lapply(names(logs), function(run) {
    p <- logs[[run]]
    div(class = "small text-muted mt-1",
      tags$b(paste0("Run ", toupper(run), ":")),
      paste0(" DIA-NN ", p$diann_version %||% "unknown",
             " | ", p$n_precursors_library %||% "?", " precursors",
             " | pg-level ", p$pg_level %||% "?",
             if (isTRUE(p$proteoforms)) " | proteoforms" else "",
             if (isTRUE(p$reanalyse))   " | reanalyse"   else "")
    )
  })
  tagList(items)
})
```

---

## New Function — `parse_diann_log()`

Add to `R/server_comparator.R`. Wraps everything in `tryCatch` — returns a
list of NULLs if parsing fails, never errors.

```r
parse_diann_log <- function(log_path) {
  tryCatch({
    lines <- readLines(log_path, warn = FALSE)

    # ── 1. Extract the command line ──────────────────────────────────────────
    # Either the first line (report_log.txt style) or after SLURM prolog
    cmd_line <- lines[grep("^/.*diann|diann-linux|diann\\.exe", lines)[1]]
    tokens   <- if (!is.na(cmd_line)) strsplit(cmd_line, "\\s+")[[1]] else character(0)

    get_flag <- function(flag) {
      idx <- which(tokens == flag)
      if (length(idx) == 0 || idx[1] >= length(tokens)) return(NULL)
      tokens[idx[1] + 1]
    }
    has_flag <- function(flag) any(tokens == flag)

    # ── 2. Extract from confirmation lines (more robust than parsing cmd) ───
    conf <- paste(lines, collapse = "\n")

    extract_conf <- function(pattern) {
      m <- regmatches(conf, regexpr(pattern, conf, perl = TRUE))
      if (length(m) == 0) return(NULL)
      trimws(sub(".*?([0-9.]+)\\s*$", "\\1", m))
    }

    # ── 3. Library size — from "X precursors generated" or "library loaded" ─
    prec_line    <- grep("precursors generated|precursors in", lines, value = TRUE)
    n_precursors <- if (length(prec_line) > 0) {
      as.integer(gsub(",", "", regmatches(
        prec_line[length(prec_line)],
        regexpr("[0-9,]+(?= precursors)", prec_line[length(prec_line)], perl = TRUE)
      )))
    } else NULL

    # ── 4. FASTA files ───────────────────────────────────────────────────────
    fasta_flags <- which(tokens == "--fasta")
    fasta_files <- if (length(fasta_flags) > 0) {
      basename(tokens[fasta_flags + 1])
    } else {
      # fall back to log confirmation lines
      m <- regmatches(lines, gregexpr("(?<=Loading FASTA ).*", lines, perl = TRUE))
      basename(unlist(m[lengths(m) > 0]))
    }

    # ── 5. Library path (--lib vs --fasta-search) ───────────────────────────
    lib_path   <- get_flag("--lib")
    out_lib    <- get_flag("--out-lib")
    fasta_mode <- has_flag("--fasta-search")

    # ── 6. DIA-NN version ────────────────────────────────────────────────────
    ver_line    <- grep("DIA-NN [0-9]", lines, value = TRUE)[1]
    diann_ver   <- if (!is.na(ver_line)) {
      regmatches(ver_line, regexpr("[0-9]+\\.[0-9]+(?:\\.[0-9]+)?", ver_line))
    } else NULL

    # ── 7. Pipeline step detection (SLURM job name) ──────────────────────────
    job_line     <- grep("JobName=|Step [0-9]+/[0-9]+:", lines, value = TRUE)[1]
    pipeline_step <- if (!is.na(job_line)) trimws(job_line) else NULL

    list(
      diann_version      = diann_ver,
      pipeline_step      = pipeline_step,

      # Precursor / library
      n_precursors_library = n_precursors,
      fasta_files          = fasta_files,
      library_path         = lib_path %||% (if (fasta_mode) "(in silico from FASTA)" else NULL),
      out_lib_path         = out_lib,
      fasta_search         = fasta_mode,

      # Search parameters
      pg_level             = as.integer(get_flag("--pg-level") %||%
                               extract_conf("pg.level [0-9]+")),
      proteoforms          = has_flag("--proteoforms"),
      reanalyse            = has_flag("--reanalyse"),
      rt_profiling         = has_flag("--rt-profiling"),
      window               = as.integer(get_flag("--window")),
      mass_acc_ms2         = as.numeric(get_flag("--mass-acc")),
      mass_acc_ms1         = as.numeric(get_flag("--mass-acc-ms1")),
      min_pr_mz            = as.numeric(get_flag("--min-pr-mz")),
      max_pr_mz            = as.numeric(get_flag("--max-pr-mz")),
      min_pr_charge        = as.integer(get_flag("--min-pr-charge")),
      max_pr_charge        = as.integer(get_flag("--max-pr-charge")),
      min_fr_mz            = as.numeric(get_flag("--min-fr-mz")),
      max_fr_mz            = as.numeric(get_flag("--max-fr-mz")),
      min_pep_len          = as.integer(get_flag("--min-pep-len")),
      max_pep_len          = as.integer(get_flag("--max-pep-len")),
      missed_cleavages     = as.integer(get_flag("--missed-cleavages")),
      qvalue               = as.numeric(get_flag("--qvalue")),
      var_mod              = get_flag("--var-mod"),
      ids_to_names         = has_flag("--ids-to-names")
    )
  }, error = function(e) {
    warning("parse_diann_log() failed: ", conditionMessage(e))
    list()   # empty list — caller checks with %||% and NULL guards
  })
}
```

---

## Integration with `build_settings_diff()`

Extend `build_settings_diff()` to incorporate log-derived parameters when
available. The log fields are merged into the existing settings list before
diff construction.

```r
build_settings_diff <- function(comparator_data) {
  s_a <- comparator_data$run_a$settings
  s_b <- comparator_data$run_b$settings

  # Merge log-derived fields into settings (log takes precedence over "unknown")
  if (!is.null(comparator_data$diann_log_a)) {
    log_a <- comparator_data$diann_log_a
    s_a$diann_version    <- log_a$diann_version   %||% s_a$diann_version
    s_a$pg_level         <- log_a$pg_level        %||% s_a$pg_level
    s_a$proteoforms      <- if (isTRUE(log_a$proteoforms)) "yes" else s_a$proteoforms
    s_a$reanalyse        <- if (isTRUE(log_a$reanalyse))   "yes" else s_a$reanalyse
    s_a$window           <- log_a$window          %||% s_a$window
    s_a$mass_acc_ms2_ppm <- log_a$mass_acc_ms2    %||% s_a$mass_acc_ms2_ppm
    s_a$mass_acc_ms1_ppm <- log_a$mass_acc_ms1    %||% s_a$mass_acc_ms1_ppm
    s_a$max_fr_mz        <- log_a$max_fr_mz       %||% s_a$max_fr_mz
    s_a$min_fr_mz        <- log_a$min_fr_mz       %||% s_a$min_fr_mz
    s_a$n_precursors_lib <- log_a$n_precursors_library %||% s_a$n_precursors_lib
    s_a$fasta_files      <- paste(log_a$fasta_files, collapse = "; ") %||% s_a$fasta_files
    s_a$library_source   <- log_a$library_path    %||% s_a$library_source
    s_a$pipeline_step    <- log_a$pipeline_step   %||% "full search"
  }
  if (!is.null(comparator_data$diann_log_b)) {
    log_b <- comparator_data$diann_log_b
    # same pattern as above for s_b
    s_b$diann_version    <- log_b$diann_version   %||% s_b$diann_version
    s_b$pg_level         <- log_b$pg_level        %||% s_b$pg_level
    s_b$proteoforms      <- if (isTRUE(log_b$proteoforms)) "yes" else s_b$proteoforms
    s_b$reanalyse        <- if (isTRUE(log_b$reanalyse))   "yes" else s_b$reanalyse
    s_b$window           <- log_b$window          %||% s_b$window
    s_b$mass_acc_ms2_ppm <- log_b$mass_acc_ms2    %||% s_b$mass_acc_ms2_ppm
    s_b$mass_acc_ms1_ppm <- log_b$mass_acc_ms1    %||% s_b$mass_acc_ms1_ppm
    s_b$max_fr_mz        <- log_b$max_fr_mz       %||% s_b$max_fr_mz
    s_b$min_fr_mz        <- log_b$min_fr_mz       %||% s_b$min_fr_mz
    s_b$n_precursors_lib <- log_b$n_precursors_library %||% s_b$n_precursors_lib
    s_b$fasta_files      <- paste(log_b$fasta_files, collapse = "; ") %||% s_b$fasta_files
    s_b$library_source   <- log_b$library_path    %||% s_b$library_source
    s_b$pipeline_step    <- log_b$pipeline_step   %||% "full search"
  }

  # Build diff rows — existing params first, then log-derived DIA-NN params
  params <- c(
    # ... existing params (Software, DE-LIMP Version, etc.) ...

    # New log-derived block (shown as a separate group in the table)
    "DIA-NN Version"        = list(s_a$diann_version,    s_b$diann_version),
    "Pipeline Step"         = list(s_a$pipeline_step,    s_b$pipeline_step),
    "Library Source"        = list(s_a$library_source,   s_b$library_source),
    "FASTA File(s)"         = list(s_a$fasta_files,      s_b$fasta_files),
    "Precursors in Library" = list(s_a$n_precursors_lib, s_b$n_precursors_lib),
    "pg-level"              = list(s_a$pg_level,         s_b$pg_level),
    "--proteoforms"         = list(s_a$proteoforms,      s_b$proteoforms),
    "--reanalyse"           = list(s_a$reanalyse,        s_b$reanalyse),
    "Scan Window (scans)"   = list(s_a$window,           s_b$window),
    "Mass Acc MS2 (ppm)"    = list(s_a$mass_acc_ms2_ppm, s_b$mass_acc_ms2_ppm),
    "Mass Acc MS1 (ppm)"    = list(s_a$mass_acc_ms1_ppm, s_b$mass_acc_ms1_ppm),
    "Max Fragment m/z"      = list(s_a$max_fr_mz,        s_b$max_fr_mz),
    "Min Fragment m/z"      = list(s_a$min_fr_mz,        s_b$min_fr_mz)
  )

  data.frame(
    Parameter = names(params),
    Run_A     = sapply(params, function(x) as.character(x[[1]] %||% "unknown")),
    Run_B     = sapply(params, function(x) as.character(x[[2]] %||% "unknown")),
    match     = mapply(function(a, b) {
      if (is.null(a) || is.null(b) || a == "unknown" || b == "unknown") "unknown"
      else if (as.character(a) == as.character(b)) "match"
      else "differs"
    }, params),
    stringsAsFactors = FALSE
  )
}
```

**Important:** The "Pipeline Step" row should trigger a **special amber warning banner**
above the settings diff table when the two runs are at different pipeline steps (e.g.,
one is "Step 1/5: Library Prediction" and the other is a full search). This is a
configuration error, not a settings difference:

```r
# Render warning above settings diff if pipeline steps differ
output$comparator_pipeline_warning <- renderUI({
  req(comparator_data())
  log_a <- comparator_data()$diann_log_a
  log_b <- comparator_data()$diann_log_b
  if (is.null(log_a) || is.null(log_b)) return(NULL)

  step_a <- log_a$pipeline_step %||% ""
  step_b <- log_b$pipeline_step %||% ""

  # Detect if one log is a library prediction step only
  is_libpred <- function(s) grepl("libpred|lib_pred|step.*1.*lib|Library Prediction",
                                  s, ignore.case = TRUE)

  if (is_libpred(step_a) || is_libpred(step_b)) {
    div(class = "alert alert-warning",
      icon("triangle-exclamation"), " ",
      tags$b("Log file mismatch: "),
      "One of the uploaded logs appears to be a library prediction step, not a full search. ",
      "The settings comparison may not reflect the actual search parameters. ",
      "Upload the final search step log (e.g., ", tags$code("step3_search.log"),
      " or the full ", tags$code("report_log.txt"), ") for an accurate comparison."
    )
  }
})
```

---

## ZIP Export — Log-Derived Files

When logs are uploaded, add two files to the ZIP:

```
delimp_run_comparison_{timestamp}/
├── ...existing files...
├── diann_params_run_a.txt       # parsed parameter summary — human-readable
├── diann_params_run_b.txt       # same for Run B
└── diann_settings_diff.csv      # just the DIA-NN rows, standalone
```

```r
write_diann_params_txt <- function(log_parsed, run_label, out_path) {
  if (is.null(log_parsed) || length(log_parsed) == 0) {
    writeLines("DIA-NN log not provided.", out_path)
    return(invisible(NULL))
  }
  p <- log_parsed
  lines <- c(
    paste0("DIA-NN Search Parameters — ", run_label),
    paste0(rep("=", 50), collapse = ""),
    paste0("DIA-NN Version:        ", p$diann_version      %||% "unknown"),
    paste0("Pipeline Step:         ", p$pipeline_step      %||% "full search"),
    "",
    "--- Library ---",
    paste0("Library Source:        ", p$library_path       %||% "(in silico from FASTA)"),
    paste0("FASTA File(s):         ", paste(p$fasta_files  %||% "unknown", collapse = "; ")),
    paste0("Precursors in Library: ", p$n_precursors_library %||% "unknown"),
    paste0("Output Library:        ", p$out_lib_path       %||% "none"),
    "",
    "--- Protein Grouping ---",
    paste0("--pg-level:            ", p$pg_level           %||% "unknown"),
    paste0("--proteoforms:         ", if (isTRUE(p$proteoforms)) "yes" else "no"),
    paste0("--ids-to-names:        ", if (isTRUE(p$ids_to_names)) "yes" else "no"),
    "",
    "--- Search Parameters ---",
    paste0("--window:              ", p$window             %||% "unknown"),
    paste0("--mass-acc (MS2 ppm):  ", p$mass_acc_ms2       %||% "unknown"),
    paste0("--mass-acc-ms1 (ppm):  ", p$mass_acc_ms1       %||% "unknown"),
    paste0("--min-pr-mz:           ", p$min_pr_mz          %||% "unknown"),
    paste0("--max-pr-mz:           ", p$max_pr_mz          %||% "unknown"),
    paste0("--min-pr-charge:       ", p$min_pr_charge      %||% "unknown"),
    paste0("--max-pr-charge:       ", p$max_pr_charge      %||% "unknown"),
    paste0("--min-fr-mz:           ", p$min_fr_mz          %||% "unknown"),
    paste0("--max-fr-mz:           ", p$max_fr_mz          %||% "unknown"),
    paste0("--min-pep-len:         ", p$min_pep_len        %||% "unknown"),
    paste0("--max-pep-len:         ", p$max_pep_len        %||% "unknown"),
    paste0("--missed-cleavages:    ", p$missed_cleavages   %||% "unknown"),
    paste0("--qvalue:              ", p$qvalue             %||% "unknown"),
    paste0("--var-mod:             ", p$var_mod            %||% "none"),
    paste0("--reanalyse:           ", if (isTRUE(p$reanalyse))    "yes" else "no"),
    paste0("--rt-profiling:        ", if (isTRUE(p$rt_profiling)) "yes" else "no"),
    paste0("--fasta-search:        ", if (isTRUE(p$fasta_search)) "yes" else "no")
  )
  writeLines(lines, out_path)
}
```

Update `build_claude_prompt()` to include library size when available:

```r
# Add after KEY FINDING block in build_claude_prompt():
diann_section <- if (!is.null(comparator_data$diann_log_a) &&
                     !is.null(comparator_data$diann_log_a$n_precursors_library)) {
  glue::glue("
DIA-NN LIBRARY SIZES:
Run A: {comparator_data$diann_log_a$n_precursors_library} precursors \\
  (pg-level {comparator_data$diann_log_a$pg_level %||% '?'}\\
  {if(isTRUE(comparator_data$diann_log_a$proteoforms)) ', --proteoforms' else ''})
Run B: {comparator_data$diann_log_b$n_precursors_library %||% 'unknown'} precursors \\
  (pg-level {comparator_data$diann_log_b$pg_level %||% '?'}\\
  {if(isTRUE(comparator_data$diann_log_b$proteoforms)) ', --proteoforms' else ''})
")
} else ""
```

---

## Reactive State Additions

```r
# In app.R reactiveValues(), add:
comparator_diann_log_a = NULL,   # parsed output of parse_diann_log() for Run A
comparator_diann_log_b = NULL    # parsed output of parse_diann_log() for Run B
```

Parse on upload:

```r
observeEvent(input$comparator_diann_log_a, {
  req(input$comparator_diann_log_a)
  values$comparator_diann_log_a <- parse_diann_log(
    input$comparator_diann_log_a$datapath
  )
})

observeEvent(input$comparator_diann_log_b, {
  req(input$comparator_diann_log_b)
  values$comparator_diann_log_b <- parse_diann_log(
    input$comparator_diann_log_b$datapath
  )
})
```

Pass into `comparator_results` so all downstream functions can access:

```r
comparator_results <- list(
  # ...existing fields...
  diann_log_a = values$comparator_diann_log_a,
  diann_log_b = values$comparator_diann_log_b
)
```

---

## Hypothesis Engine — New Rule

Add as **Rule 0** (highest priority — checked before all others):

```r
# Rule 0: Library size mismatch → search space difference
if (!is.null(log_a$n_precursors_library) && !is.null(log_b$n_precursors_library)) {
  ratio <- max(log_a$n_precursors_library, log_b$n_precursors_library) /
           min(log_a$n_precursors_library, log_b$n_precursors_library)
  if (ratio > 1.2) {
    return(list(
      hypothesis = glue::glue(
        "Library size differs substantially \\
        ({format(log_a$n_precursors_library, big.mark=',')} vs \\
         {format(log_b$n_precursors_library, big.mark=',')} precursors). \\
        {if(isTRUE(log_a$proteoforms) != isTRUE(log_b$proteoforms))
          'One run used --proteoforms (isoform-aware FDR), the other did not. ' else ''}\\
        {if(!is.null(log_a$pg_level) && !is.null(log_b$pg_level) &&
             log_a$pg_level != log_b$pg_level)
          paste0('pg-level differs (', log_a$pg_level, ' vs ', log_b$pg_level, '). ') else ''}\\
        Different precursor pools affect which peptides are available for \\
        each protein's DPC-Quant rollup."
      ),
      hypothesis_category = "Library / search space",
      confidence          = "High"
    ))
  }
}
```

---

## Session Save / Load

Log data is not saved to session RDS (file paths are ephemeral). On load, 
`comparator_diann_log_a` and `comparator_diann_log_b` will be NULL. The UI
will show the log upload panel again, prompting re-upload if needed.

---

## Testing Checklist

- [ ] `parse_diann_log()` correctly parses `report_log.txt` (all-in-one run)
- [ ] `parse_diann_log()` correctly parses SLURM `.out` file (multi-step run)
- [ ] Pipeline step warning banner appears when one log is a `libpred` step
- [ ] Settings diff table shows DIA-NN block with amber highlighting on `--proteoforms` row when it differs
- [ ] Library size difference triggers Rule 0 hypothesis on a known discordant pair
- [ ] `write_diann_params_txt()` produces readable output for both log types
- [ ] ZIP export includes `diann_params_run_a.txt` when log was uploaded, skips it gracefully when not
- [ ] `parse_diann_log()` returns empty list (not error) on a malformed / unrelated file
- [ ] Log status summary UI renders correctly after upload

