# Local Docker Execution Addendum
## DDA (Sage) + XL-MS (MeroX + xiSearch) — Local Docker Support

> **Extends:** `DDA_SAGE_CASANOVO_SPEC.md`, `XLMS_INTEGRATION_SPEC.md`,
> and `THREE_MODE_SWITCHER_ADDENDUM.md`
>
> **Priority:** P2 — implement after Hive pipeline is confirmed working (Phase 5
> of IMPLEMENTATION_PLAN.md). Estimated effort: 2-3 hrs.
>
> **When implementing:** Read this addendum after the relevant parent spec.
> This addendum defines only what changes for local Docker execution.
> All Hive/SLURM logic in the parent specs is unchanged.

---

## What This Enables

Local Docker users (collaborators, labs without HPC access) can run DDA and
XL-MS searches on their own machines. The workflow is identical to Hive from
the user's perspective — same UI, same settings, same results format, same
Claude export ZIP. The only difference under the hood is serial execution
instead of SLURM array jobs.

**Practical scope for local execution:**

| Experiment type | Local Docker | Hive |
|---|---|---|
| XL-MS, purified complex, 2-4 samples | ✅ ~1-2 hrs | ✅ ~20 min |
| XL-MS, cell lysate, 12+ samples | ⚠️ ~12+ hrs, 64 GB RAM | ✅ ~2 hrs |
| DDA, standard proteomics, 6-12 samples | ✅ ~15-30 min | ✅ ~5 min |
| DDA, whole proteome, 24+ samples | ⚠️ slow, feasible | ✅ fast |

Casanovo (DDA de novo) requires a GPU and is **not supported in local Docker**
unless the user has an NVIDIA GPU and nvidia-container-toolkit installed.
The Casanovo checkbox is hidden in Docker mode.

---

## 1. Execution Mode Detection

Add to `app.R` startup block, after the existing `is_hive` detection:

```r
# Existing
is_hive <- nzchar(Sys.which("sbatch"))

# New — detect local Docker tool availability
is_docker_xlms <- all(file.exists(c(
  getOption("delimp.merox_jar",  "/opt/xlms/MeroX_2025.jar"),
  getOption("delimp.xi_jar",     "/opt/xlms/xiSearch.jar"),
  getOption("delimp.xifdr_jar",  "/opt/xlms/xiFDR.jar"),
  getOption("delimp.timsrust_bin", "/usr/local/bin/timsrust_mgf")
)))
is_docker_dda <- file.exists(
  getOption("delimp.sage_bin", "/usr/local/bin/sage")
)

# Resolved execution mode per pipeline
xlms_exec_mode <- if (is_hive) "hive" else if (is_docker_xlms) "docker" else "unavailable"
dda_exec_mode  <- if (is_hive) "hive" else if (is_docker_dda)  "docker" else "unavailable"
```

Pass these into `reactiveValues()`:

```r
values <- reactiveValues(
  # ... existing values ...
  xlms_exec_mode = xlms_exec_mode,
  dda_exec_mode  = dda_exec_mode
)
```

---

## 2. config.yml — Docker Tool Paths

Add Docker paths alongside existing Hive paths. The `get_tool_path()` helper
(§3 below) selects the right one at runtime.

```yaml
tools:
  # Hive paths (when sbatch available)
  merox_jar_hive:    "/quobyte/proteomics-grp/tools/xlms/MeroX_2025.jar"
  xi_jar_hive:       "/quobyte/proteomics-grp/tools/xlms/xiSearch.jar"
  xifdr_jar_hive:    "/quobyte/proteomics-grp/tools/xlms/xiFDR.jar"
  timsrust_bin_hive: "/quobyte/proteomics-grp/tools/timsrust_mgf"
  sage_bin_hive:     "/quobyte/proteomics-grp/tools/sage"

  # Docker paths (when running locally)
  merox_jar_docker:    "/opt/xlms/MeroX_2025.jar"
  xi_jar_docker:       "/opt/xlms/xiSearch.jar"
  xifdr_jar_docker:    "/opt/xlms/xiFDR.jar"
  timsrust_bin_docker: "/usr/local/bin/timsrust_mgf"
  sage_bin_docker:     "/usr/local/bin/sage"

slurm:
  account:     "genome-center-grp"
  partition:   "high"
  java_module: "openjdk"

docker:
  # Memory cap for local Java processes (GB)
  # Set lower than total system RAM — leave headroom for OS + R
  # Recommended: (total RAM - 4) GB
  # E.g., 16 GB MacBook → 12, 32 GB Mac Studio → 28, 64 GB workstation → 56
  java_mem_gb: 12
```

---

## 3. Tool Path Helper

Add to `R/helpers_hpc.R` (create this file if it doesn't exist — extract
`check_array_job()` here too):

```r
#' Get the correct tool path for current execution environment
#'
#' @param tool_name Base name (e.g. "merox_jar", "sage_bin")
#' @param config Loaded config list
#' @param is_hive Logical — TRUE if running on Hive with sbatch
#' @return Character path to the tool
get_tool_path <- function(tool_name, config, is_hive) {
  suffix <- if (is_hive) "hive" else "docker"
  key    <- paste0(tool_name, "_", suffix)
  path   <- config$tools[[key]]
  if (is.null(path) || !nzchar(path)) {
    stop(sprintf("Tool path not configured: %s (looking for config$tools$%s)",
                 tool_name, key))
  }
  path
}

#' Get Java memory flag for current environment
#'
#' @param config Loaded config list
#' @param is_hive Logical
#' @param override_gb Optional integer to override config
#' @return Character like "-Xmx56g"
get_java_xmx <- function(config, is_hive, override_gb = NULL) {
  mem_gb <- override_gb %||%
    if (is_hive) config$slurm$java_mem_gb %||% 56
    else         config$docker$java_mem_gb %||% 12
  paste0("-Xmx", mem_gb, "g")
}

#' Check if a SLURM array job is complete (Hive only)
check_array_job <- function(jid) {
  if (is.null(jid)) return("not submitted")
  out <- system2("squeue",
    args = c("-j", jid, "-h", "-o", "%T"),
    stdout = TRUE, stderr = FALSE)
  if (length(out) == 0) return("COMPLETED")
  states <- unique(trimws(out))
  if (any(states == "FAILED"))  return("FAILED")
  if (any(states == "RUNNING")) return("RUNNING")
  if (any(states == "PENDING")) return("PENDING")
  "COMPLETED"
}
```

---

## 4. XL-MS Serial Execution (Docker Mode)

Replace SLURM array jobs with blocking `system2()` calls. Add these functions
to `R/server_xlms.R`:

```r
# ── Docker mode: convert one .d file to MGF ──────────────────────────────────
run_timsrust_local <- function(d_path, mgf_out, timsrust_bin) {
  result <- system2(timsrust_bin,
    args   = c(shQuote(d_path), shQuote(mgf_out)),
    stdout = TRUE, stderr = TRUE)
  if (!file.exists(mgf_out)) {
    stop("timsrust conversion failed for: ", d_path, "\n", paste(result, collapse = "\n"))
  }
  mgf_out
}

# ── Docker mode: run MeroX on one MGF file ───────────────────────────────────
run_merox_local <- function(mgf_path, fasta_path, mxf_path,
                             output_dir, merox_jar, mem_gb = 12,
                             progress_callback = NULL) {
  name   <- tools::file_path_sans_ext(basename(mgf_path))
  zhrm   <- file.path(output_dir, paste0(name, ".zhrm"))
  log_f  <- file.path(output_dir, paste0(name, ".log"))

  if (!is.null(progress_callback)) progress_callback(paste("MeroX:", name))

  result <- system2("java",
    args = c(paste0("-Xmx", mem_gb, "g"),
             "-jar", shQuote(merox_jar),
             shQuote(mgf_path),
             shQuote(fasta_path),
             shQuote(mxf_path),
             shQuote(zhrm),
             shQuote(log_f)),
    stdout = TRUE, stderr = TRUE)

  list(zhrm = zhrm, log = log_f, output = result)
}

# ── Docker mode: run xiSearch on one MGF file ────────────────────────────────
run_xi_local <- function(mgf_path, fasta_path, xi_config_path,
                          output_dir, xi_jar, mem_gb = 12,
                          progress_callback = NULL) {
  name   <- tools::file_path_sans_ext(basename(mgf_path))
  out_csv <- file.path(output_dir, paste0(name, "_xi.csv"))

  if (!is.null(progress_callback)) progress_callback(paste("xiSearch:", name))

  system2("java",
    args = c(paste0("-Xmx", mem_gb, "g"),
             "-cp", shQuote(xi_jar),
             "rappsilber.applications.Xi",
             paste0("--config=", xi_config_path),
             paste0("--peaks=", shQuote(mgf_path)),
             paste0("--fasta=", shQuote(fasta_path)),
             paste0("--output=", shQuote(out_csv)),
             "--locale=en"),
    stdout = TRUE, stderr = TRUE)

  out_csv
}

# ── Docker mode: run xiFDR on combined CSV ───────────────────────────────────
run_xifdr_local <- function(combined_csv, fdr_out_dir, xifdr_jar,
                             fdr_threshold, mem_gb = 8) {
  system2("java",
    args = c(paste0("-Xmx", mem_gb, "g"),
             "-jar", shQuote(xifdr_jar),
             paste0("input=",       combined_csv),
             paste0("outputpath=",  fdr_out_dir),
             paste0("psmfdr=",      fdr_threshold),
             paste0("pepfdr=",      fdr_threshold),
             paste0("proteinfdr=",  fdr_threshold),
             "reportfactor=0", "unique=true", "csv=true"),
    stdout = TRUE, stderr = TRUE)
}

# ── Docker mode: full serial pipeline ────────────────────────────────────────
run_xlms_docker <- function(xlms, input, config, session) {
  mem_gb      <- config$docker$java_mem_gb %||% 12
  merox_jar   <- get_tool_path("merox_jar",   config, is_hive = FALSE)
  xi_jar      <- get_tool_path("xi_jar",      config, is_hive = FALSE)
  xifdr_jar   <- get_tool_path("xifdr_jar",   config, is_hive = FALSE)
  timsrust_bin <- get_tool_path("timsrust_bin", config, is_hive = FALSE)

  output_dir    <- xlms$output_dir
  merox_out_dir <- file.path(output_dir, "merox_per_file")
  xi_out_dir    <- file.path(output_dir, "xi_per_file")
  mgf_dir       <- file.path(output_dir, "mgf")
  dir.create(merox_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(xi_out_dir,    recursive = TRUE, showWarnings = FALSE)
  dir.create(mgf_dir,       recursive = TRUE, showWarnings = FALSE)

  n_files <- length(xlms$raw_files)

  withProgress(message = "Running XL-MS search locally...", value = 0, {

    # Stage 0: convert all .d → MGF
    mgf_files <- character(n_files)
    for (i in seq_along(xlms$raw_files)) {
      incProgress(1 / (n_files * 3),
        detail = paste0("Converting (", i, "/", n_files, "): ",
                        basename(xlms$raw_files[i])))
      name         <- tools::file_path_sans_ext(basename(xlms$raw_files[i]))
      mgf_files[i] <- run_timsrust_local(
        xlms$raw_files[i],
        file.path(mgf_dir, paste0(name, ".mgf")),
        timsrust_bin
      )
    }

    # Stages 1a + 1b: MeroX and xiSearch — serial per file
    # (no parallelism locally — Java processes are heavy)
    zhrm_files  <- character(n_files)
    xi_csv_files <- character(n_files)

    for (i in seq_along(mgf_files)) {
      base_progress <- 1/3 + (i - 1) / n_files * (2/3)

      # MeroX
      incProgress(1 / (n_files * 3),
        detail = paste0("MeroX (", i, "/", n_files, "): ", basename(mgf_files[i])))
      merox_result <- run_merox_local(
        mgf_files[i], xlms$fasta_path, xlms$mxf_path,
        merox_out_dir, merox_jar, mem_gb
      )
      zhrm_files[i] <- merox_result$zhrm

      # xiSearch
      incProgress(1 / (n_files * 3),
        detail = paste0("xiSearch (", i, "/", n_files, "): ", basename(mgf_files[i])))
      xi_csv_files[i] <- run_xi_local(
        mgf_files[i], xlms$fasta_path, xlms$xi_config_path,
        xi_out_dir, xi_jar, mem_gb
      )
    }

    # Stage 2a: MeroX merge
    incProgress(0, detail = "Merging MeroX results...")
    merged_merox <- file.path(output_dir, "merox_combined.csv")
    system2("java", args = c(
      "-jar", shQuote(merox_jar),
      "--export-csv", shQuote(merged_merox),
      shQuote(paste(zhrm_files, collapse = " "))
    ))

    # Stage 2b: xiFDR on combined xi CSVs
    incProgress(0, detail = "Running xiFDR...")
    combined_xi <- file.path(output_dir, "xi_combined.csv")
    fdr_out_dir <- file.path(output_dir, "xi_fdr")
    dir.create(fdr_out_dir, recursive = TRUE, showWarnings = FALSE)

    # Concatenate xi CSVs
    header <- readLines(xi_csv_files[1], n = 1)
    writeLines(header, combined_xi)
    for (f in xi_csv_files) {
      lines <- readLines(f)
      write(lines[-1], file = combined_xi, append = TRUE)
    }

    run_xifdr_local(combined_xi, fdr_out_dir, xifdr_jar,
                    as.numeric(input$xlms_fdr), mem_gb = 8)
  })

  # Parse and merge — same as Hive path
  parse_and_merge_results(
    merox_csv = merged_merox,
    xi_fdr_dir = fdr_out_dir
  )
}
```

---

## 5. Main Submit Handler — Mode-Aware Routing

Replace the single `observeEvent(input$xlms_submit)` in `server_xlms.R` with
a mode-aware version:

```r
observeEvent(input$xlms_submit, {
  req(xlms$raw_files, xlms$fasta_path, xlms$output_dir)

  # Generate settings files (same for both modes)
  mxf_path       <- file.path(xlms$output_dir, "merox_settings.mxf")
  xi_config_path <- file.path(xlms$output_dir, "xi_settings.config")
  generate_merox_mxf(reactiveValuesToList(xlms), mxf_path)
  generate_xi_config(reactiveValuesToList(xlms),
                     xlms$fasta_path, character(0), xi_config_path)
  xlms$mxf_path       <- mxf_path
  xlms$xi_config_path <- xi_config_path

  if (values$xlms_exec_mode == "hive") {
    # ── Hive path: SLURM array jobs (non-blocking) ─────────────────────
    # ... existing submit_conversion_array_job() etc. from XLMS spec §7.3 ...
    submit_xlms_hive(xlms, input, config, session)

  } else if (values$xlms_exec_mode == "docker") {
    # ── Docker path: serial blocking execution ──────────────────────────
    run_xlms_docker(xlms, input, config, session)

  } else {
    showNotification(
      "XL-MS search tools not found. Please check your installation.",
      type = "error", duration = 10)
  }
})
```

Apply the same routing pattern to the DDA submit handler in `server_dda.R`:

```r
observeEvent(input$run_dda_search, {
  if (values$dda_exec_mode == "hive") {
    submit_sage_hive(values, input, config, session)
  } else if (values$dda_exec_mode == "docker") {
    run_sage_docker(values, input, config, session)
  } else {
    showNotification("Sage not found. Please check your installation.",
                     type = "error", duration = 10)
  }
})
```

---

## 6. DDA Serial Execution (Docker Mode)

Sage is much simpler — one blocking call:

```r
run_sage_docker <- function(values, input, config, session) {
  sage_bin   <- get_tool_path("sage_bin", config, is_hive = FALSE)
  output_dir <- values$dda$output_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sage_config_path <- generate_sage_config(
    fasta_path        = input$dda_fasta_path,
    raw_paths         = values$dda$raw_files,
    output_dir        = output_dir,
    preset            = input$dda_preset,
    missed_cleavages  = input$dda_missed_cleavages,
    precursor_tol_ppm = input$dda_precursor_tol,
    fragment_tol_da   = input$dda_fragment_tol
  )

  withProgress(message = "Running Sage search locally...", value = 0.1, {
    incProgress(0.1, detail = "Submitting Sage search...")

    result <- system2(sage_bin,
      args   = c(shQuote(sage_config_path),
                 shQuote(values$dda$raw_files)),
      stdout = TRUE, stderr = TRUE)

    incProgress(0.7, detail = "Parsing results...")
    values$dda$sage_output <- result
  })

  # Parse results — same as Hive path
  results_path <- file.path(output_dir, "results.sage.tsv")
  lfq_path     <- file.path(output_dir, "lfq.tsv")

  if (file.exists(results_path)) {
    parsed <- parse_sage_results(results_path, lfq_path)
    values$dda$sage_psms    <- parsed$psms
    values$dda$lfq_wide     <- parsed$lfq_wide
    values$dda$protein_meta <- parsed$protein_meta
    run_dda_pipeline(values, input)
    showNotification("Sage search complete!", type = "message", duration = 5)
  } else {
    showNotification(
      paste("Sage search failed. Output:", paste(result, collapse = "\n")),
      type = "error", duration = 15)
  }
}
```

---

## 7. UI Warnings for Docker Mode

### XL-MS Setup Panel — Docker warning

Add to `xlms_setup_ui()` in `R/ui_xlms.R`, inside the Setup panel:

```r
# Show only in Docker mode
uiOutput("xlms_exec_mode_banner")
```

```r
# In server_xlms.R
output$xlms_exec_mode_banner <- renderUI({
  if (values$xlms_exec_mode == "hive") return(NULL)

  if (values$xlms_exec_mode == "docker") {
    div(class = "alert alert-info mt-2",
      icon("laptop"),
      strong(" Local mode."),
      " Searches run serially on this machine. ",
      "Recommended for small experiments (≤4 samples, purified protein complex). ",
      "Memory available: ", round(memory.limit() / 1024, 0), " GB",
      " — configured limit: ",
      strong(paste0(config$docker$java_mem_gb, " GB")), " per Java process.",
      br(),
      "For large experiments (whole proteome, 12+ samples) use the Hive HPC cluster."
    )
  } else {
    div(class = "alert alert-danger mt-2",
      icon("circle-exclamation"),
      strong(" XL-MS tools not found."),
      " Please ensure MeroX, xiSearch, xiFDR, and timsrust are installed. ",
      "See README for installation instructions."
    )
  }
})
```

### DDA Panel — Casanovo hidden in Docker mode

In `ui_dda` or `server_dda.R`, hide the Casanovo checkbox when not on Hive:

```r
# In DDA setup UI
conditionalPanel(
  condition = "input.acquisition_mode === 'dda'",
  # Casanovo only available on Hive (GPU)
  if (is_hive) {
    checkboxInput("dda_casanovo_enable",
      "Run Casanovo de novo sequencing (GPU required)", FALSE)
  } else {
    div(class = "text-muted small mt-2",
      icon("circle-info"),
      " Casanovo de novo sequencing requires GPU — available on Hive only.")
  }
)
```

---

## 8. Dockerfile Changes

Add to the existing `Dockerfile` (or `Dockerfile.search` if separate):

```dockerfile
# ── Java runtime (MeroX + xiSearch + xiFDR) ───────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    openjdk-17-jre-headless \
    && rm -rf /var/lib/apt/lists/*

# ── XL-MS search engines ────────────────────────────────────────────────────
RUN mkdir -p /opt/xlms
# These JARs must be present in the build context under tools/
# Download them from:
#   MeroX:   https://github.com/Sinz-Lab-CSMS/MeroX
#   xiSearch: https://github.com/Rappsilber-Laboratory/XiSearch
#   xiFDR:   https://github.com/Rappsilber-Laboratory/xiFDR
COPY tools/MeroX_2025.jar /opt/xlms/MeroX_2025.jar
COPY tools/xiSearch.jar   /opt/xlms/xiSearch.jar
COPY tools/xiFDR.jar      /opt/xlms/xiFDR.jar

# ── timsrust MGF converter (pre-compiled static Linux binary from Hive) ─────
# Copy from Hive: scp hive:/quobyte/proteomics-grp/tools/timsrust_mgf tools/
COPY tools/timsrust_mgf /usr/local/bin/timsrust_mgf
RUN chmod +x /usr/local/bin/timsrust_mgf

# ── Sage (DDA) ──────────────────────────────────────────────────────────────
RUN wget -q \
    https://github.com/lazear/sage/releases/latest/download/sage-x86_64-unknown-linux-gnu.tar.gz \
    -O /tmp/sage.tar.gz \
    && tar -xzf /tmp/sage.tar.gz -C /usr/local/bin/ \
    && chmod +x /usr/local/bin/sage \
    && rm /tmp/sage.tar.gz

# ── R packages for XL-MS visualization ─────────────────────────────────────
RUN Rscript -e "install.packages(c('visNetwork', 'ggraph', 'igraph', 'yaml'), \
    repos = 'https://cloud.r-project.org/')"
```

### Build context requirement

The JARs cannot be downloaded automatically (no public direct download URLs
without clicking through license agreements). Add a `tools/` directory to the
repo (gitignored) with a `README_TOOLS.md` explaining where to get them:

```
tools/
├── README_TOOLS.md     ← instructions for downloading JARs
├── MeroX_2025.jar      ← download from github.com/Sinz-Lab-CSMS/MeroX
├── xiSearch.jar        ← download from rappsilberlab.org/software/xisearch
├── xiFDR.jar           ← download from github.com/Rappsilber-Laboratory/xiFDR
└── timsrust_mgf        ← copy from Hive: /quobyte/proteomics-grp/tools/timsrust_mgf
```

Add to `.gitignore`:
```
tools/*.jar
tools/timsrust_mgf
```

Add a `setup_tools.sh` convenience script:

```bash
#!/bin/bash
# setup_tools.sh — download and stage tools for local Docker build
# Run once before building the Docker image

set -e
mkdir -p tools

echo "Downloading Sage..."
wget -q https://github.com/lazear/sage/releases/latest/download/\
sage-x86_64-unknown-linux-gnu.tar.gz -O /tmp/sage.tar.gz
tar -xzf /tmp/sage.tar.gz -C tools/ sage
rm /tmp/sage.tar.gz
echo "  ✓ sage"

echo ""
echo "Manual steps required — cannot auto-download (license agreements):"
echo "  1. MeroX_2025.jar → https://github.com/Sinz-Lab-CSMS/MeroX/releases"
echo "     Save as: tools/MeroX_2025.jar"
echo "  2. xiSearch.jar → https://www.rappsilberlab.org/software/xisearch/"
echo "     Save as: tools/xiSearch.jar"
echo "  3. xiFDR.jar → https://github.com/Rappsilber-Laboratory/xiFDR/releases"
echo "     Save as: tools/xiFDR.jar"
echo "  4. timsrust_mgf → copy from Hive:"
echo "     scp USER@hive.hpc.ucdavis.edu:/quobyte/proteomics-grp/tools/timsrust_mgf tools/"
echo ""
echo "Once all tools are in place, run: docker build -t delimp ."
```

---

## 9. Image Size Impact

| Component | Size |
|---|---|
| openjdk-17-jre-headless | ~180 MB |
| MeroX_2025.jar | ~50 MB |
| xiSearch.jar | ~30 MB |
| xiFDR.jar | ~15 MB |
| timsrust_mgf (static binary) | ~8 MB |
| Sage binary | ~12 MB |
| visNetwork + ggraph + igraph R packages | ~25 MB |
| **Total addition** | **~320 MB** |

Current DE-LIMP image is ~3-4 GB. New image will be ~3.3-4.3 GB — within
HuggingFace Spaces limits. However, XL-MS and DDA tools should only be
included in the full local Docker image, not the HF Spaces image. Use a
build arg to control this:

```dockerfile
ARG INCLUDE_SEARCH_TOOLS=false

RUN if [ "$INCLUDE_SEARCH_TOOLS" = "true" ]; then \
      apt-get install -y openjdk-17-jre-headless; \
    fi

# ... similar for each tool block ...
```

Build commands:
```bash
# HuggingFace Spaces (no search tools — smaller image)
docker build -t delimp:hf .

# Local use (includes DDA + XL-MS tools)
docker build --build-arg INCLUDE_SEARCH_TOOLS=true -t delimp:local .
```

---

## 10. Progress Display — Docker vs Hive

In Docker mode, `withProgress()` gives live feedback since execution is
synchronous in the R session. The job monitor badge display (§7.4 of XLMS spec)
is replaced with a simple progress bar — no polling needed.

In Hive mode, the job monitor badges remain as designed.

```r
output$xlms_job_monitor <- renderUI({
  if (values$xlms_exec_mode == "docker") {
    # Docker: no job monitor needed — withProgress() handles feedback
    # Show a simple "search running" indicator if active
    if (isTRUE(xlms$docker_running)) {
      div(class = "alert alert-info",
        icon("spinner", class = "fa-spin"), " Search running locally...")
    }
  } else {
    # Hive: full five-stage badge display from XLMS spec §7.4
    # ... existing badge renderUI code ...
  }
})
```

---

## 11. Files Modified by This Addendum

| File | Change |
|---|---|
| `app.R` | Add Docker mode detection (`is_docker_xlms`, `is_docker_dda`, exec mode vars) |
| `config.yml` | Add `_docker` tool paths and `docker.java_mem_gb` |
| `R/helpers_hpc.R` | New file: `get_tool_path()`, `get_java_xmx()`, `check_array_job()` |
| `R/server_xlms.R` | Add `run_timsrust_local()`, `run_merox_local()`, `run_xi_local()`, `run_xifdr_local()`, `run_xlms_docker()`; mode-aware submit handler |
| `R/server_dda.R` | Add `run_sage_docker()`; mode-aware submit handler |
| `R/ui_xlms.R` | Add `output$xlms_exec_mode_banner` |
| `Dockerfile` | Add Java, JARs, timsrust, Sage, R packages with `INCLUDE_SEARCH_TOOLS` build arg |
| `tools/README_TOOLS.md` | New file: JAR download instructions |
| `setup_tools.sh` | New file: tool staging script |
| `.gitignore` | Add `tools/*.jar`, `tools/timsrust_mgf`, `tools/sage` |

---

## 12. Implementation Order

Do this after Phase 5 (end-to-end Hive test) passes:

1. Add `helpers_hpc.R` with `get_tool_path()` — extract `check_array_job()` here too
2. Update all tool path references in `server_xlms.R` and `server_dda.R` to use
   `get_tool_path()` instead of hardcoded paths — this is a cleanup that benefits
   both modes
3. Add Docker detection to `app.R`
4. Add Docker execution functions to `server_xlms.R` and `server_dda.R`
5. Update submit handlers to route by exec mode
6. Update `config.yml` with Docker paths
7. Update `Dockerfile` with `INCLUDE_SEARCH_TOOLS` build arg
8. Build local Docker image, test with a small purified complex XL-MS dataset
9. Confirm Hive path still works (regression test)
