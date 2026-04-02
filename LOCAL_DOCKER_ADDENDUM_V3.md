# Local Docker Execution Addendum — v3 (Phase 1)
## DDA (Sage + Casanovo) + XL-MS (MeroX + xiSearch) — Local Docker Support

> **Extends:** `DDA_SAGE_CASANOVO_SPEC.md`, `XLMS_INTEGRATION_SPEC.md`,
> `THREE_MODE_SWITCHER_ADDENDUM.md`
>
> **Scope:** Phase 1 only — CPU-first local execution. No GPU-in-Docker
> complexity. Ship this after Hive pipeline is confirmed working.
> Estimated implementation effort: 2-3 hrs.
>
> **Supersedes:** LOCAL_DOCKER_ADDENDUM.md v1 and v2 — read this instead.
> Previous versions were over-engineered. This is what to actually build.

---

## The One Thing to Understand Before Reading Further

GPU support in Docker is complicated, platform-dependent, and helps
very few of your users. Don't build it yet.

**What actually matters for local Docker:**
- Sage, MeroX, and xiSearch are all CPU-only tools — they don't need a GPU
- Casanovo works on CPU, just slower (~2-4 hrs vs ~10 min with GPU)
- For Casanovo, Apple Silicon Mac users get MPS acceleration — but only
  when running DE-LIMP natively, not inside Docker

For Phase 1, do this: run everything on CPU in Docker, document the native
Mac path for users who want MPS, and move on. GPU-in-Docker is Phase 2,
only if users ask for it.

---

## 1. What Needs a GPU and What Doesn't

| Component | GPU needed? | Local Docker behavior |
|---|---|---|
| **Sage** (DDA search) | ❌ No | Runs fast on CPU — no change needed |
| **MeroX** (XL-MS) | ❌ No | Runs on CPU — no change needed |
| **xiSearch** (XL-MS) | ❌ No | Runs on CPU — no change needed |
| **timsrust** (conversion) | ❌ No | Runs in seconds on CPU |
| **Casanovo** (DDA de novo) | ⚠️ Optional | CPU works, just slower |
| **Cascadia** (DIA de novo) | ⚠️ Optional | CPU works, just slower |

---

## 2. Execution Mode Detection

Add to `app.R` startup, after the existing `is_hive` detection:

```r
# Existing
is_hive <- nzchar(Sys.which("sbatch"))

# New — detect local tool availability
is_docker_xlms <- all(file.exists(c(
  getOption("delimp.merox_jar",    "/opt/xlms/MeroX_2025.jar"),
  getOption("delimp.xi_jar",       "/opt/xlms/xiSearch.jar"),
  getOption("delimp.xifdr_jar",    "/opt/xlms/xiFDR.jar"),
  getOption("delimp.timsrust_bin", "/usr/local/bin/timsrust_mgf")
))) && nzchar(Sys.which("java"))

is_docker_dda <- file.exists(
  getOption("delimp.sage_bin", "/usr/local/bin/sage")
)

# Execution mode per pipeline — used throughout server_xlms.R and server_dda.R
xlms_exec_mode <- if (is_hive) "hive" else if (is_docker_xlms) "docker" else "unavailable"
dda_exec_mode  <- if (is_hive) "hive" else if (is_docker_dda)  "docker" else "unavailable"
```

Add to `reactiveValues()` in `app.R`:

```r
values <- reactiveValues(
  # ... existing values ...
  xlms_exec_mode = xlms_exec_mode,
  dda_exec_mode  = dda_exec_mode
)
```

---

## 3. Accelerator Detection for Casanovo

Simple three-tier detection. Works reliably in all environments.
The only tricky case is Docker on Mac — handle it explicitly.

Add to `R/helpers_hpc.R`:

```r
#' Detect best available compute accelerator for PyTorch (Casanovo/Cascadia)
#'
#' Returns "gpu", "mps", or "cpu".
#' Tests actual availability from within the current process,
#' not just whether hardware exists on the host.
detect_accelerator <- function() {

  # 1. NVIDIA GPU — works if nvidia-container-toolkit is set up on the host
  #    and container was launched with --gpus all, OR running natively on Linux
  nvidia_ok <- tryCatch({
    out <- system2("nvidia-smi",
      args   = c("--query-gpu=name", "--format=csv,noheader"),
      stdout = TRUE, stderr = FALSE)
    length(out) > 0 &&
      nzchar(trimws(out[1])) &&
      !grepl("error|failed|not found", out[1], ignore.case = TRUE)
  }, error = function(e) FALSE)

  if (nvidia_ok) return("gpu")

  # 2. Apple Silicon MPS — native macOS only, never works inside Docker
  #    Docker on Mac runs a Linux VM that cannot access the Metal GPU
  in_docker <- file.exists("/.dockerenv")
  is_macos  <- Sys.info()[["sysname"]] == "Darwin"

  if (is_macos && !in_docker) {
    arch <- tryCatch(
      trimws(system2("uname", "-m", stdout = TRUE)),
      error = function(e) ""
    )
    if (arch == "arm64") return("mps")
  }

  # 3. CPU fallback — always works
  "cpu"
}

# Detect once at session start — cache the result
local_accelerator <- detect_accelerator()
```

---

## 4. Casanovo Configuration by Accelerator

Two things that must be set correctly for Casanovo to work without errors:

**`PYTORCH_ENABLE_MPS_FALLBACK=1`** — Casanovo crashes on Apple Silicon
without this. Set it in every code path that invokes Casanovo.

**`predict_batch_size`** — the default (1024) causes out-of-memory on CPU.
Auto-tune based on accelerator.

Add to `R/helpers_dda.R`:

```r
#' Get Casanovo predict_batch_size for current accelerator
get_casanovo_batch_size <- function(accelerator) {
  switch(accelerator,
    "gpu"  = 1024,   # NVIDIA default — fast
    "mps"  = 512,    # Apple Silicon — slightly reduced for stability
    "cpu"  = 64,     # CPU — must be small to avoid OOM
    64               # safe default
  )
}

#' Generate casanovo config YAML
#' Always sets PYTORCH_ENABLE_MPS_FALLBACK so MPS users don't hit silent crashes
generate_casanovo_config <- function(accelerator, output_path,
                                      n_beams = 1, top_match = 1) {
  config <- list(
    accelerator        = accelerator,
    predict_batch_size = get_casanovo_batch_size(accelerator),
    n_beams            = n_beams,
    top_match          = top_match
  )
  yaml::write_yaml(config, output_path)
  invisible(output_path)
}

#' Build casanovo system2() call with all required env vars
#' @param mgf_path Input MGF file
#' @param output_path Output mzTab file
#' @param config_path Generated casanovo config YAML
#' @param weights_path Pre-cached .ckpt file (avoids GitHub rate limiting)
#' @param accelerator "gpu", "mps", or "cpu"
run_casanovo <- function(mgf_path, output_path, config_path,
                          weights_path, accelerator) {
  # PYTORCH_ENABLE_MPS_FALLBACK=1 must be set before importing torch
  # Otherwise Apple Silicon users hit: NotImplementedError: aten::index.Tensor
  old_env <- Sys.getenv("PYTORCH_ENABLE_MPS_FALLBACK")
  Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK = "1")
  on.exit(Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK = old_env))

  result <- system2(
    "casanovo",
    args = c(
      "sequence",
      shQuote(mgf_path),
      "--model",    shQuote(weights_path),  # always explicit — no GitHub download
      "--config",   shQuote(config_path),
      "--output",   shQuote(output_path)
    ),
    stdout = TRUE,
    stderr = TRUE
  )

  list(
    success = file.exists(output_path),
    output  = result
  )
}
```

---

## 5. Pre-Cache Casanovo Model Weights

Casanovo auto-downloads weights from GitHub at runtime — subject to rate
limits (60 requests/hour/IP). Pre-cache them so they're always available.

**On Hive (run once manually):**

```bash
mkdir -p /share/proteomics/tools/casanovo

# Download weights for the pinned version
CASANOVO_VERSION="4.2.0"
wget \
  "https://github.com/Noble-Lab/casanovo/releases/download/\
v${CASANOVO_VERSION}/casanovo_v${CASANOVO_VERSION//./_}.ckpt" \
  -O "/share/proteomics/tools/casanovo/casanovo_v${CASANOVO_VERSION//./_}.ckpt"

echo "Weights cached at: /share/proteomics/tools/casanovo/"
```

**In `config.yml`:**

```yaml
tools:
  # Casanovo — version pinned (weights are version-specific)
  casanovo_version:      "4.2.0"
  casanovo_weights_hive: "/share/proteomics/tools/casanovo/casanovo_v4_2_0.ckpt"
  casanovo_weights_docker: "/opt/casanovo/weights/casanovo_v4_2_0.ckpt"
```

**Version check helper** — warn if installed Casanovo doesn't match weights:

```r
check_casanovo_version <- function(config) {
  expected <- config$tools$casanovo_version %||% "4.2.0"
  installed <- tryCatch({
    out <- system2("casanovo", "--version", stdout = TRUE, stderr = FALSE)
    regmatches(out, regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", out))
  }, error = function(e) NA_character_)

  if (is.na(installed)) return(invisible(NULL))

  # Check major version match only (minor versions share weights)
  exp_major  <- strsplit(expected,  "\\.")[[1]][1]
  inst_major <- strsplit(installed, "\\.")[[1]][1]

  if (exp_major != inst_major) {
    warning(sprintf(
      paste("Casanovo version mismatch: installed v%s, weights for v%s.",
            "Update weights or pin Casanovo version in config.yml."),
      installed, expected
    ))
  }
}
```

Call at startup: `check_casanovo_version(config)`

---

## 6. MeroX Memory Warning — Based on FASTA Size, Not Sample Count

MeroX memory scales with database size, not sample count. A 5-protein
complex and a full human proteome both have "12 samples" but require
4 GB vs 256 GB respectively.

Add to `R/helpers_hpc.R`:

```r
#' Estimate MeroX RAM requirement from FASTA protein count
#' MeroX search space scales quadratically with database size.
#' Empirical: ~5 MB per protein is a conservative estimate.
estimate_merox_memory_gb <- function(fasta_path) {
  if (is.null(fasta_path) || !nzchar(fasta_path) ||
      !file.exists(fasta_path)) return(NA_integer_)

  n_proteins <- tryCatch({
    as.integer(system2(
      "grep", args = c("-c", "^>", shQuote(fasta_path)),
      stdout = TRUE, stderr = FALSE
    ))
  }, error = function(e) NA_integer_)

  if (is.na(n_proteins) || n_proteins == 0) return(NA_integer_)

  # 5 MB per protein, minimum 4 GB
  as.integer(max(4L, ceiling(n_proteins * 5L / 1000L)))
}
```

Use in the Docker mode warning banner in `server_xlms.R`:

```r
output$xlms_exec_mode_banner <- renderUI({
  req(values$xlms_exec_mode == "docker")

  est_gb  <- estimate_merox_memory_gb(input$xlms_fasta_path)
  limit   <- config$docker$java_mem_gb %||% 12
  n_files <- length(xlms$raw_files)

  # Memory status
  mem_status <- if (!is.na(est_gb)) {
    if (est_gb > limit) {
      div(class = "alert alert-danger small mt-1",
        icon("triangle-exclamation"),
        strong(sprintf(" Estimated RAM needed: ~%d GB", est_gb)),
        sprintf(" — exceeds your configured limit (%d GB).", limit),
        " Increase ", code("docker.java_mem_gb"), " in config.yml,",
        " or use Hive for this experiment."
      )
    } else {
      div(class = "alert alert-success small mt-1",
        icon("check-circle"),
        sprintf(" Estimated RAM: ~%d GB — within your %d GB limit.",
                est_gb, limit)
      )
    }
  } else NULL

  tagList(
    div(class = "alert alert-info mt-2",
      icon("laptop"),
      strong(" Local mode — serial execution."),
      br(),
      "Searches run one file at a time on this machine. ",
      "Best suited for small experiments (purified complex, ≤4 samples). ",
      "For large experiments (whole proteome, 12+ samples) use Hive."
    ),
    mem_status
  )
})
```

---

## 7. De Novo Accelerator Info Panel

Show users what hardware Casanovo will use before they start.
Add to the DDA Setup panel (only visible when Casanovo checkbox is ticked):

```r
output$denovo_accelerator_info <- renderUI({
  req(input$dda_casanovo_enable)

  acc       <- local_accelerator
  in_docker <- file.exists("/.dockerenv")
  is_mac    <- Sys.info()[["sysname"]] == "Darwin"

  # Mac + Docker: MPS unavailable — explain why and offer native path
  if (in_docker && is_mac) {
    return(div(class = "alert alert-info small mt-1",
      icon("apple"),
      strong(" Mac + Docker: GPU not available."),
      " Apple Silicon GPU (MPS) cannot be accessed inside Docker on macOS. ",
      "For faster de novo sequencing, run DE-LIMP natively — ",
      "see ", tags$a("installation guide",
        href = "https://github.com/bsphinney/DE-LIMP#installation",
        target = "_blank"), ".",
      " CPU mode will be used here (~2-4 hrs per dataset)."
    ))
  }

  switch(acc,
    "gpu" = div(class = "alert alert-success small mt-1",
      icon("bolt"),
      strong(" NVIDIA GPU detected."),
      " De novo sequencing: ~5-10 min per dataset."
    ),
    "mps" = div(class = "alert alert-success small mt-1",
      icon("apple"),
      strong(" Apple Silicon GPU (MPS)."),
      " De novo sequencing: ~20-40 min per dataset."
    ),
    "cpu" = div(class = "alert alert-warning small mt-1",
      icon("clock"),
      strong(" No GPU detected — CPU mode."),
      " De novo sequencing: ~2-4 hrs per dataset. ",
      "For faster results: use Hive, or on Apple Silicon ",
      "run DE-LIMP natively outside Docker."
    )
  )
})
```

---

## 8. Serial Execution Functions for Docker Mode

These replace the SLURM array jobs when running locally.
Add to `R/server_xlms.R` and `R/server_dda.R`.

### XL-MS Serial Pipeline

```r
# In server_xlms.R
run_xlms_docker <- function(xlms, input, config, session) {
  mem_gb       <- config$docker$java_mem_gb %||% 12
  merox_jar    <- get_tool_path("merox_jar",    config, is_hive = FALSE)
  xi_jar       <- get_tool_path("xi_jar",       config, is_hive = FALSE)
  xifdr_jar    <- get_tool_path("xifdr_jar",    config, is_hive = FALSE)
  timsrust_bin <- get_tool_path("timsrust_bin", config, is_hive = FALSE)

  output_dir    <- xlms$output_dir
  mgf_dir       <- file.path(output_dir, "mgf")
  merox_out_dir <- file.path(output_dir, "merox_per_file")
  xi_out_dir    <- file.path(output_dir, "xi_per_file")
  fdr_out_dir   <- file.path(output_dir, "xi_fdr")

  dir.create(mgf_dir,       recursive = TRUE, showWarnings = FALSE)
  dir.create(merox_out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(xi_out_dir,    recursive = TRUE, showWarnings = FALSE)
  dir.create(fdr_out_dir,   recursive = TRUE, showWarnings = FALSE)

  n <- length(xlms$raw_files)

  withProgress(message = "Running XL-MS search locally...", value = 0, {

    # ── Stage 0: Convert .d → MGF ────────────────────────────────────────
    mgf_files <- character(n)
    for (i in seq_len(n)) {
      name         <- tools::file_path_sans_ext(basename(xlms$raw_files[i]))
      mgf_files[i] <- file.path(mgf_dir, paste0(name, ".mgf"))
      incProgress(1 / (n * 3),
        detail = sprintf("Converting (%d/%d): %s", i, n, name))
      system2(timsrust_bin,
        args = c(shQuote(xlms$raw_files[i]), shQuote(mgf_files[i])))
    }

    # ── Stage 1: MeroX + xiSearch per file (serial) ──────────────────────
    zhrm_files  <- character(n)
    xi_csv_files <- character(n)

    for (i in seq_len(n)) {
      name <- tools::file_path_sans_ext(basename(mgf_files[i]))

      # MeroX
      incProgress(1 / (n * 3),
        detail = sprintf("MeroX (%d/%d): %s", i, n, name))
      zhrm_files[i] <- file.path(merox_out_dir, paste0(name, ".zhrm"))
      system2("java",
        args = c(paste0("-Xmx", mem_gb - 2, "g"),
                 "-jar", shQuote(merox_jar),
                 shQuote(mgf_files[i]),
                 shQuote(xlms$fasta_path),
                 shQuote(xlms$mxf_path),
                 shQuote(zhrm_files[i]),
                 shQuote(file.path(merox_out_dir, paste0(name, ".log")))))

      # xiSearch
      incProgress(1 / (n * 3),
        detail = sprintf("xiSearch (%d/%d): %s", i, n, name))
      xi_csv_files[i] <- file.path(xi_out_dir, paste0(name, "_xi.csv"))
      system2("java",
        args = c(paste0("-Xmx", mem_gb - 2, "g"),
                 "-cp", shQuote(xi_jar),
                 "rappsilber.applications.Xi",
                 paste0("--config=",   xlms$xi_config_path),
                 paste0("--peaks=",    shQuote(mgf_files[i])),
                 paste0("--fasta=",    shQuote(xlms$fasta_path)),
                 paste0("--output=",   shQuote(xi_csv_files[i])),
                 "--locale=en"))
    }

    # ── Stage 2a: MeroX merge ─────────────────────────────────────────────
    incProgress(0, detail = "Merging MeroX results...")
    merged_merox <- file.path(output_dir, "merox_combined.csv")
    system2("java", args = c(
      "-jar", shQuote(merox_jar),
      "--export-csv", shQuote(merged_merox),
      shQuote(paste(zhrm_files, collapse = " "))
    ))

    # ── Stage 2b: Concatenate xi CSVs → xiFDR ────────────────────────────
    incProgress(0, detail = "Running xiFDR...")
    combined_xi <- file.path(output_dir, "xi_combined.csv")

    # Header from first file, data from all
    header <- readLines(xi_csv_files[1], n = 1)
    writeLines(header, combined_xi)
    for (f in xi_csv_files) {
      data_lines <- readLines(f)[-1]  # drop header
      write(data_lines, file = combined_xi, append = TRUE)
    }

    system2("java",
      args = c(paste0("-Xmx8g"),
               "-jar", shQuote(xifdr_jar),
               paste0("input=",       combined_xi),
               paste0("outputpath=",  fdr_out_dir),
               paste0("psmfdr=",      as.numeric(input$xlms_fdr)),
               paste0("pepfdr=",      as.numeric(input$xlms_fdr)),
               paste0("proteinfdr=",  as.numeric(input$xlms_fdr)),
               "reportfactor=0", "unique=true", "csv=true"))
  })

  # Parse results — same functions as Hive path
  xi_fdr_csv <- list.files(fdr_out_dir, pattern = "\\.csv$",
                            full.names = TRUE)[1]
  parse_and_merge_results(
    merox_csv    = merged_merox,
    xi_fdr_csv   = xi_fdr_csv
  )
}
```

### DDA Serial Pipeline (Sage only — Casanovo optional)

```r
# In server_dda.R
run_sage_docker <- function(values, input, config, session) {
  sage_bin   <- get_tool_path("sage_bin", config, is_hive = FALSE)
  output_dir <- values$dda$output_dir
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  sage_config_path <- generate_sage_config(
    fasta_path        = values$dda$fasta_path,
    raw_paths         = values$dda$raw_files,
    output_dir        = output_dir,
    preset            = input$dda_preset,
    missed_cleavages  = input$dda_missed_cleavages,
    precursor_tol_ppm = input$dda_precursor_tol,
    fragment_tol_da   = input$dda_fragment_tol
  )

  withProgress(message = "Running Sage search locally...", {
    incProgress(0.2, detail = "Searching...")

    system2(sage_bin,
      args   = c(shQuote(sage_config_path),
                 sapply(values$dda$raw_files, shQuote)),
      stdout = file.path(output_dir, "sage.log"),
      stderr = file.path(output_dir, "sage.log"))

    incProgress(0.7, detail = "Parsing results...")
  })

  results_path <- file.path(output_dir, "results.sage.tsv")
  lfq_path     <- file.path(output_dir, "lfq.tsv")

  if (!file.exists(results_path)) {
    showNotification(
      paste("Sage search failed. Check:",
            file.path(output_dir, "sage.log")),
      type = "error", duration = 15)
    return(invisible(NULL))
  }

  parsed <- parse_sage_results(results_path, lfq_path)
  values$dda$sage_psms    <- parsed$psms
  values$dda$lfq_wide     <- parsed$lfq_wide
  values$dda$protein_meta <- parsed$protein_meta
  run_dda_pipeline(values, input)

  showNotification("Sage search complete!", type = "message", duration = 5)
}
```

---

## 9. Mode-Aware Submit Handlers

Replace the single submit observers in each server module with
two-branch versions:

```r
# In server_xlms.R
observeEvent(input$xlms_submit, {
  req(xlms$raw_files, xlms$fasta_path, xlms$output_dir)

  # Generate settings files — same for both modes
  mxf_path        <- file.path(xlms$output_dir, "merox_settings.mxf")
  xi_config_path  <- file.path(xlms$output_dir, "xi_settings.config")
  generate_merox_mxf(reactiveValuesToList(xlms), mxf_path)
  generate_xi_config(reactiveValuesToList(xlms),
                     xlms$fasta_path, character(0), xi_config_path)
  xlms$mxf_path       <- mxf_path
  xlms$xi_config_path <- xi_config_path

  if (values$xlms_exec_mode == "hive") {
    submit_xlms_hive(xlms, input, config, session)    # SLURM arrays
  } else if (values$xlms_exec_mode == "docker") {
    run_xlms_docker(xlms, input, config, session)     # serial blocking
  } else {
    showNotification(
      paste("XL-MS tools not found.",
            "Check config.yml tool paths and verify Java is installed."),
      type = "error", duration = 10)
  }
})

# In server_dda.R
observeEvent(input$run_dda_search, {
  if (values$dda_exec_mode == "hive") {
    submit_sage_hive(values, input, config, session)
  } else if (values$dda_exec_mode == "docker") {
    run_sage_docker(values, input, config, session)
  } else {
    showNotification(
      paste("Sage not found.",
            "Check config.yml sage_bin_docker path."),
      type = "error", duration = 10)
  }
})
```

---

## 10. Dockerfile Changes

Minimal additions — no CUDA, no GPU layers, just the tools:

```dockerfile
# ── Java (MeroX + xiSearch + xiFDR) ───────────────────────────────────────
RUN apt-get update && \
    apt-get install -y --no-install-recommends openjdk-17-jre-headless && \
    rm -rf /var/lib/apt/lists/*

# ── XL-MS JARs ────────────────────────────────────────────────────────────
# Place downloaded JARs in tools/ before building
RUN mkdir -p /opt/xlms
COPY tools/MeroX_2025.jar /opt/xlms/MeroX_2025.jar
COPY tools/xiSearch.jar   /opt/xlms/xiSearch.jar
COPY tools/xiFDR.jar      /opt/xlms/xiFDR.jar

# ── timsrust MGF converter ─────────────────────────────────────────────────
# Copy the x86_64 binary from Hive
# For Mac ARM64 Docker: timsconvert is the fallback (see §11)
COPY tools/timsrust_mgf /usr/local/bin/timsrust_mgf
RUN chmod +x /usr/local/bin/timsrust_mgf

# ── Sage (DDA) ────────────────────────────────────────────────────────────
RUN wget -q \
    https://github.com/lazear/sage/releases/latest/download/\
sage-x86_64-unknown-linux-gnu.tar.gz \
    -O /tmp/sage.tar.gz && \
    tar -xzf /tmp/sage.tar.gz -C /usr/local/bin/ sage && \
    chmod +x /usr/local/bin/sage && \
    rm /tmp/sage.tar.gz

# ── Casanovo model weights ─────────────────────────────────────────────────
# Pre-cached to avoid GitHub rate limiting at runtime
RUN mkdir -p /opt/casanovo/weights
COPY tools/casanovo_weights/casanovo_v4_2_0.ckpt \
     /opt/casanovo/weights/casanovo_v4_2_0.ckpt

# ── R packages for XL-MS visualization ────────────────────────────────────
RUN Rscript -e "install.packages(\
    c('visNetwork', 'ggraph', 'igraph', 'yaml'), \
    repos = 'https://cloud.r-project.org/')"

# ── Environment ────────────────────────────────────────────────────────────
# Required for Casanovo on Apple Silicon MPS
# Safe to set always — no-op on non-Mac platforms
ENV PYTORCH_ENABLE_MPS_FALLBACK=1
```

Build commands:
```bash
# Standard (CPU — works everywhere)
docker build -t delimp:local .
docker run -p 3838:3838 delimp:local

# With NVIDIA GPU (Linux/WSL2 only — requires nvidia-container-toolkit on host)
# Only do this if you specifically need GPU for Casanovo
docker run --gpus all -p 3838:3838 delimp:local
```

---

## 11. Mac ARM64 Docker — timsrust Fallback

The timsrust binary from Hive is x86_64 and won't run inside Docker on
Apple Silicon Macs. Handle this with a runtime fallback to timsconvert,
which is pure Python and installs cleanly on any architecture:

```dockerfile
# Add to Dockerfile
RUN pip install timsconvert --break-system-packages
```

```r
# In R/helpers_hpc.R — used by run_xlms_docker()
get_conversion_command <- function(d_path, mgf_out) {
  timsrust_bin <- "/usr/local/bin/timsrust_mgf"

  # Try timsrust first
  timsrust_works <- tryCatch({
    system2(timsrust_bin, "--help",
            stdout = FALSE, stderr = FALSE) == 0
  }, error = function(e) FALSE)

  if (timsrust_works) {
    list(cmd  = timsrust_bin,
         args = c(shQuote(d_path), shQuote(mgf_out)))
  } else {
    # Fallback: timsconvert → mzML (both MeroX and xiSearch accept mzML)
    mzml_out <- sub("\\.mgf$", ".mzML", mgf_out)
    list(cmd  = "timsconvert",
         args = c("--input", shQuote(d_path),
                  "--output", shQuote(dirname(mzml_out)),
                  "--ms2_only", "True"),
         output_file = mzml_out,   # mzML not MGF
         is_mzml     = TRUE)
  }
}
```

This means Mac Docker users silently get timsconvert → mzML instead of
timsrust → MGF, and MeroX/xiSearch receive mzML files instead. Both
engines accept mzML natively so the output is identical.

---

## 12. config.yml — Complete Tool Paths Section

```yaml
tools:
  # Hive paths (when sbatch available)
  merox_jar_hive:          "/quobyte/proteomics-grp/tools/xlms/MeroX_2025.jar"
  xi_jar_hive:             "/quobyte/proteomics-grp/tools/xlms/xiSearch.jar"
  xifdr_jar_hive:          "/quobyte/proteomics-grp/tools/xlms/xiFDR.jar"
  timsrust_bin_hive:       "/quobyte/proteomics-grp/tools/timsrust_mgf"
  sage_bin_hive:           "/quobyte/proteomics-grp/tools/sage"
  casanovo_weights_hive:   "/share/proteomics/tools/casanovo/casanovo_v4_2_0.ckpt"

  # Docker / local paths
  merox_jar_docker:        "/opt/xlms/MeroX_2025.jar"
  xi_jar_docker:           "/opt/xlms/xiSearch.jar"
  xifdr_jar_docker:        "/opt/xlms/xiFDR.jar"
  timsrust_bin_docker:     "/usr/local/bin/timsrust_mgf"
  sage_bin_docker:         "/usr/local/bin/sage"
  casanovo_weights_docker: "/opt/casanovo/weights/casanovo_v4_2_0.ckpt"

  # Casanovo version — must match weights file
  casanovo_version: "4.2.0"

slurm:
  account:     "genome-center-grp"
  partition:   "high"
  java_module: "openjdk"

docker:
  # Java memory per process (GB)
  # Set to roughly (total system RAM - 4) GB
  # 16 GB laptop → 12, 32 GB workstation → 28, 64 GB workstation → 56
  java_mem_gb: 12
```

---

## 13. setup_tools.sh — Tool Staging Script

Run this once before building the Docker image:

```bash
#!/bin/bash
# setup_tools.sh — stage tools for local Docker build
# Run from the DE-LIMP repo root before: docker build -t delimp:local .
set -e
mkdir -p tools tools/casanovo_weights

# ── Sage (auto-download) ─────────────────────────────────────────────────
echo "Downloading Sage..."
wget -q \
  https://github.com/lazear/sage/releases/latest/download/\
sage-x86_64-unknown-linux-gnu.tar.gz \
  -O /tmp/sage.tar.gz
tar -xzf /tmp/sage.tar.gz -C tools/ sage
rm /tmp/sage.tar.gz
chmod +x tools/sage
echo "  ✓ sage"

# ── timsrust (copy from Hive) ────────────────────────────────────────────
echo ""
echo "Copying timsrust from Hive..."
echo "  Run: scp USER@hive.hpc.ucdavis.edu:\\"
echo "         /quobyte/proteomics-grp/tools/timsrust_mgf tools/timsrust_mgf"
echo "  Then: chmod +x tools/timsrust_mgf"

# ── Casanovo weights (auto-download) ─────────────────────────────────────
echo ""
CKPT="tools/casanovo_weights/casanovo_v4_2_0.ckpt"
if [ ! -f "$CKPT" ]; then
  echo "Downloading Casanovo weights (~500 MB)..."
  wget -q \
    "https://github.com/Noble-Lab/casanovo/releases/download/\
v4.2.0/casanovo_v4_2_0.ckpt" \
    -O "$CKPT"
  echo "  ✓ casanovo_v4_2_0.ckpt"
else
  echo "  ✓ Casanovo weights already present"
fi

# ── JARs — must be downloaded manually ──────────────────────────────────
echo ""
echo "Manual downloads required (license agreements):"
echo ""
echo "  1. MeroX_2025.jar"
echo "     → https://github.com/Sinz-Lab-CSMS/MeroX/releases"
echo "     Save as: tools/MeroX_2025.jar"
echo ""
echo "  2. xiSearch.jar"
echo "     → https://www.rappsilberlab.org/software/xisearch/"
echo "     Save as: tools/xiSearch.jar"
echo ""
echo "  3. xiFDR.jar"
echo "     → https://github.com/Rappsilber-Laboratory/xiFDR/releases"
echo "     Save as: tools/xiFDR.jar"
echo ""
echo "Once all tools are staged, build with:"
echo "  docker build -t delimp:local ."
echo ""
echo "Run with:"
echo "  docker run -p 3838:3838 delimp:local"
echo ""
echo "GPU NOTE: MeroX, xiSearch, and Sage don't need a GPU."
echo "Casanovo de novo sequencing runs on CPU (~2-4 hrs) or GPU (~10 min)."
echo "To use an NVIDIA GPU (Linux/WSL2 only, requires nvidia-container-toolkit):"
echo "  docker run --gpus all -p 3838:3838 delimp:local"
echo "Apple Silicon Mac users: run DE-LIMP natively for GPU acceleration."
```

---

## 14. README Section — User-Facing Documentation

Add this to `README.md` under a "Local Installation" section:

```markdown
### Local Docker Installation

DE-LIMP can run locally using Docker. DDA and XL-MS searches run on CPU —
no GPU required for core functionality.

**Prerequisites:** Docker Desktop, ~5 GB free disk space

**Setup:**
\`\`\`bash
git clone https://github.com/bsphinney/DE-LIMP
cd DE-LIMP
bash setup_tools.sh   # downloads Sage, Casanovo weights; prompts for JARs
docker build -t delimp:local .
docker run -p 3838:3838 delimp:local
# Open http://localhost:3838
\`\`\`

**Performance guide:**

| Task | Local (CPU) | Hive (HPC) |
|---|---|---|
| XL-MS search, purified complex | ~1-2 hrs | ~20 min |
| XL-MS search, cell lysate | ≥12 hrs ⚠️ | ~2 hrs |
| DDA search (Sage) | ~15-30 min | ~5 min |
| De novo sequencing (Casanovo) | ~2-4 hrs | ~10 min |

**Apple Silicon Macs:** For faster de novo sequencing (Casanovo/Cascadia),
run DE-LIMP natively outside Docker to use Apple Silicon GPU acceleration:
\`\`\`bash
Rscript -e "shiny::runApp('app.R', port=3838)"
\`\`\`
MPS acceleration gives ~20-40 min for de novo vs ~2-4 hrs on CPU.

**NVIDIA GPU users (Linux/WSL2):**
If you have nvidia-container-toolkit installed:
\`\`\`bash
docker run --gpus all -p 3838:3838 delimp:local
\`\`\`
This speeds up Casanovo de novo sequencing only — all other tools are unaffected.
```

---

## 15. Files Created or Modified

| File | Change |
|---|---|
| `app.R` | Add Docker mode detection; `local_accelerator <- detect_accelerator()` |
| `config.yml` | Add Docker tool paths; `casanovo_version`; `docker.java_mem_gb` |
| `R/helpers_hpc.R` | New: `detect_accelerator()`, `get_tool_path()`, `get_java_xmx()`, `estimate_merox_memory_gb()`, `get_conversion_command()`, `check_casanovo_version()` |
| `R/helpers_dda.R` | Add `generate_casanovo_config()`, `get_casanovo_batch_size()`, `run_casanovo()` |
| `R/server_xlms.R` | Add `run_xlms_docker()`; mode-aware submit handler |
| `R/server_dda.R` | Add `run_sage_docker()`; mode-aware submit handler |
| `R/ui_xlms.R` | Add `output$xlms_exec_mode_banner` with FASTA-based memory warning |
| `R/ui_dda.R` | Add `output$denovo_accelerator_info` |
| `Dockerfile` | Add Java, JARs, timsrust, Sage, Casanovo weights, timsconvert fallback |
| `setup_tools.sh` | New: tool staging script |
| `tools/README_TOOLS.md` | New: JAR download instructions |
| `.gitignore` | Add `tools/*.jar`, `tools/timsrust_mgf`, `tools/casanovo_weights/` |
| `README.md` | Add Local Docker Installation section |

---

## 16. What This Does NOT Include (Phase 2, Later)

Do not implement these now. Add only if a specific user requests it:

- **NVIDIA GPU in Docker on Windows WSL2** — works in theory, requires IT
  cooperation for Hyper-V/WSL2, low value for academic labs
- **timsrust ARM64 binary for Mac Docker** — timsconvert fallback handles
  Mac Docker users adequately
- **Casanovo weight auto-update system** — manual version bump in config.yml
  is sufficient for now
- **AMD GPU / ROCm support** — uncommon in proteomics labs
- **Multi-GPU** — inference-only workloads don't benefit
