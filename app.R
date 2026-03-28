# ==============================================================================
#  DE-LIMP: Differential Expression & Limpa Proteomics App
#  (Formerly LIMP-D)
#  Status: Production Ready (Hugging Face Compatible v1.2)

# ==============================================================================

# Set CRAN mirror to avoid interactive popup (especially in VS Code)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Suppress Bioconductor internet validation (HPC compute nodes often lack internet)
options(BIOCONDUCTOR_ONLINE_VERSION_DIAGNOSIS = FALSE)

# --- 1. AUTO-INSTALLATION & SETUP ---
# IMPORTANT: Install packages BEFORE loading libraries to avoid conflicts

# Skip package installation inside containers (no internet on HPC compute nodes)
is_container <- nzchar(Sys.getenv("APPTAINER_CONTAINER", "")) ||
                nzchar(Sys.getenv("SINGULARITY_CONTAINER", "")) ||
                file.exists("/.dockerenv")

# Install BiocManager if needed (skip in containers)
if (!is_container && !requireNamespace("BiocManager", quietly = TRUE)) {
  message("BiocManager not found. Installing...")
  install.packages("BiocManager", quiet = TRUE)
}

# Check for limpa and install if needed (skip in containers — already installed)
if (!is_container && !requireNamespace("limpa", quietly = TRUE)) {
  message("Package 'limpa' is missing. Attempting installation...")

  r_version <- getRversion()
  bioc_version <- as.character(BiocManager::version())
  message(paste0("R version: ", r_version, ", Bioconductor version: ", bioc_version))

  # Try installing from Bioconductor (may require devel version)
  limpa_installed <- tryCatch({
    suppressWarnings({
      BiocManager::install("limpa", ask = FALSE, update = FALSE, quiet = TRUE)
    })
    requireNamespace("limpa", quietly = TRUE)
  }, error = function(e) FALSE)

  # If standard Bioconductor failed, try development version
  if (!limpa_installed) {
    message("limpa not found in release version. Trying Bioconductor devel...")
    limpa_installed <- tryCatch({
      suppressWarnings({
        BiocManager::install(version = "devel", ask = FALSE, update = FALSE)
        BiocManager::install("limpa", ask = FALSE, update = FALSE, quiet = TRUE)
      })
      requireNamespace("limpa", quietly = TRUE)
    }, error = function(e) FALSE)
  }

  # Final check
  if (!limpa_installed) {

    # Detect platform for specific instructions
    os_type <- Sys.info()["sysname"]
    download_url <- if (os_type == "Darwin") {
      "https://cloud.r-project.org/bin/macosx/"
    } else if (os_type == "Windows") {
      "https://cloud.r-project.org/bin/windows/base/"
    } else {
      "https://cloud.r-project.org/bin/linux/"
    }

    stop(paste0(
      "\n\n╔════════════════════════════════════════════════════════════════╗\n",
      "║          LIMPA INSTALLATION FAILED - R UPGRADE NEEDED          ║\n",
      "╚════════════════════════════════════════════════════════════════╝\n\n",
      "Current setup:\n",
      "  • R version: ", r_version, " (NEED: 4.5+)\n",
      "  • Bioconductor: ", bioc_version, " (NEED: 3.22+)\n",
      "  • Platform: ", os_type, "\n\n",
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
      "UPGRADE INSTRUCTIONS:\n",
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
      "1. Download R 4.5+ for your platform:\n",
      "   ", download_url, "\n\n",
      if (os_type == "Darwin") {
        "2. Install the .pkg file (R will be upgraded in-place)\n"
      } else if (os_type == "Windows") {
        "2. Run the installer .exe file\n"
      } else {
        "2. Follow platform-specific installation instructions\n"
      },
      "\n3. Restart VSCode/RStudio completely\n",
      "\n4. Verify upgrade by running: R.version.string\n",
      "\n5. Rerun this script - limpa will install automatically\n\n",
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
      "More info: https://bioconductor.org/packages/limpa/\n\n"
    ))
  } else {
    message("✓ limpa installed successfully!")
  }
}

# Required packages (excluding limpa which was handled above)
# Core packages: app won't start without these
core_pkgs <- c("shiny", "bslib", "readr", "tibble", "dplyr", "tidyr",
               "ggplot2", "httr2", "rhandsontable", "DT", "arrow",
               "ComplexHeatmap", "shinyjs", "plotly", "stringr", "limma",
               "AnnotationDbi", "ggridges", "ggrepel", "markdown", "curl",
               "glue", "data.table")

# Optional packages: app runs without them (features disabled gracefully)
optional_pkgs <- c("clusterProfiler", "enrichplot", "org.Hs.eg.db", "org.Mm.eg.db",
                    "KSEAapp", "ggseqlogo", "MOFA2")

# Only install truly missing packages (don't update already-loaded packages)
missing_pkgs <- character(0)
for (pkg in c(core_pkgs, optional_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, pkg)
  }
}

if (length(missing_pkgs) > 0) {
  message(paste0("Installing missing packages: ", paste(missing_pkgs, collapse = ", ")))
  tryCatch(
    BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE, quiet = TRUE),
    error = function(e) {
      message("Note: Could not install packages (no internet?). Checking core dependencies...")
    }
  )
  # Verify core packages are available — these are required
  still_missing_core <- character(0)
  for (pkg in core_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      still_missing_core <- c(still_missing_core, pkg)
    }
  }
  if (length(still_missing_core) > 0) {
    stop(paste0("Missing required packages: ", paste(still_missing_core, collapse = ", "),
                "\nRebuild the container image to include these packages."))
  }
  # Report optional packages that are unavailable
  still_missing_opt <- character(0)
  for (pkg in optional_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      still_missing_opt <- c(still_missing_opt, pkg)
    }
  }
  if (length(still_missing_opt) > 0) {
    message(paste0("Optional packages unavailable (features disabled): ",
                   paste(still_missing_opt, collapse = ", ")))
  }
}

# --- 2. SERVER CONFIGURATION ---
# Set repos (BiocManager::repositories() requires internet — skip in containers)
if (is_container) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
} else {
  options(repos = c(
    tryCatch(BiocManager::repositories(), error = function(e) NULL),
    CRAN = "https://cloud.r-project.org"
  ))
}

library(shiny)
library(bslib)

# Verify bslib version supports responsive UI components
if (packageVersion("bslib") < "0.5.0") {
  stop(paste0(
    "bslib >= 0.5.0 required for responsive UI components.\n",
    "Current version: ", packageVersion("bslib"), "\n",
    "Please upgrade: install.packages('bslib')"
  ))
}

library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(httr2)      # CRITICAL for AI Chat
library(rhandsontable)
library(DT)     
library(arrow)  
library(ComplexHeatmap)
library(shinyjs)
library(plotly)
library(stringr)
library(AnnotationDbi)
library(ggrepel)
library(markdown) # Needed for AI formatting

# Optional packages — load if available, features degrade gracefully
gsea_available <- requireNamespace("clusterProfiler", quietly = TRUE) &&
                  requireNamespace("enrichplot", quietly = TRUE)
if (gsea_available) {
  library(clusterProfiler)
  library(enrichplot)
} else {
  message("Note: clusterProfiler/enrichplot not available — GSEA tab will be disabled")
}

options(shiny.maxRequestSize = 5000 * 1024^2)  # 5 GB upload limit

# Detect Hugging Face Spaces environment (SPACE_ID is set automatically by HF)
is_hf_space <- nzchar(Sys.getenv("SPACE_ID", ""))

# Detect search backends (Docker local + HPC SSH/SLURM)
# Disabled on Hugging Face Spaces — search tab not useful in cloud environment
local_sbatch_path <- Sys.which("sbatch")
if (!nzchar(local_sbatch_path)) {
  # Inside Apptainer container, sbatch may not be on PATH but exists on CVMFS
  cvmfs_sbatch <- c(
    "/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin/sbatch",
    "/usr/bin/sbatch", "/opt/slurm/bin/sbatch", "/usr/local/bin/sbatch"
  )
  for (p in cvmfs_sbatch) {
    if (file.exists(p)) { local_sbatch_path <- p; break }
  }
}
local_sbatch <- nzchar(local_sbatch_path)
hpc_available <- !is_hf_space && (local_sbatch || nzchar(Sys.which("ssh")))

# Docker backend detection
docker_available <- FALSE
docker_config <- list(diann_image = "diann:2.0")
if (!is_hf_space && nzchar(Sys.which("docker"))) {
  docker_available <- tryCatch({
    system2("docker", "info", stdout = TRUE, stderr = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(e) FALSE)
  # Optional config from ~/.delimp_docker.conf
  docker_conf_path <- path.expand("~/.delimp_docker.conf")
  if (docker_available && file.exists(docker_conf_path)) {
    docker_config <- tryCatch(
      jsonlite::fromJSON(docker_conf_path),
      error = function(e) docker_config
    )
  }
}

# Local DIA-NN binary detection (embedded in Docker container or installed on host)
local_diann <- nzchar(Sys.which("diann")) || nzchar(Sys.which("diann-linux"))
delimp_data_dir <- Sys.getenv("DELIMP_DATA_DIR", "")

# Combined flag — at least one backend available
search_enabled <- docker_available || hpc_available || local_diann

# Environment label — helps users tell where the app is running
# Apptainer sets APPTAINER_CONTAINER / SINGULARITY_CONTAINER automatically
is_apptainer <- nzchar(Sys.getenv("APPTAINER_CONTAINER", "")) ||
                nzchar(Sys.getenv("SINGULARITY_CONTAINER", ""))
deploy_env <- if (is_hf_space) {
  "Hugging Face"
} else if (local_diann && nzchar(delimp_data_dir)) {
  "Docker"
} else if (is_apptainer || local_sbatch) {
  "HPC"
} else {
  "Local"
}

# Detect core facility mode — activated by config directory + staff.yml
# The config directory is created by delimp-server setup, never present on
# HF Spaces or regular local installs
core_facility_config_dir <- Sys.getenv("DELIMP_CORE_DIR", "/srv/delimp")
is_core_facility <- dir.exists(core_facility_config_dir) &&
                    file.exists(file.path(core_facility_config_dir, "staff.yml"))

# Load core facility configuration
cf_config <- NULL
if (is_core_facility) {
  # Load required packages for core facility features
  for (pkg in c("DBI", "RSQLite", "yaml", "uuid", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE, repos = "https://cloud.r-project.org")
    }
  }
  library(DBI)
  library(RSQLite)

  cf_config <- list(
    staff    = yaml::read_yaml(file.path(core_facility_config_dir, "staff.yml")),
    qc       = if (file.exists(file.path(core_facility_config_dir, "qc_config.yml")))
                 yaml::read_yaml(file.path(core_facility_config_dir, "qc_config.yml"))
               else NULL,
    db_path  = file.path(core_facility_config_dir, "delimp.db"),
    reports_dir  = file.path(core_facility_config_dir, "reports"),
    state_dir    = file.path(core_facility_config_dir, "state"),
    template_qmd = file.path(core_facility_config_dir, "report_template.qmd")
  )

  # Ensure directories exist
  dir.create(cf_config$reports_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(cf_config$state_dir, showWarnings = FALSE, recursive = TRUE)

  # Initialize SQLite DB if needed
  db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
  cf_init_db(db)  # defined in helpers_facility.R
  DBI::dbDisconnect(db)

  message("Core facility mode: ENABLED (", core_facility_config_dir, ")")
}

# Conditionally load search-related packages
if (search_enabled) {
  for (pkg in c("shinyFiles", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
    }
  }
  library(shinyFiles)
  library(jsonlite)
  # Migrate user-local speclib cache to shared volume if available
  tryCatch(speclib_cache_migrate(), error = function(e) NULL)
}

# Verify Limpa installation
if (!requireNamespace("limpa", quietly = TRUE)) {
  os_type <- Sys.info()["sysname"]
  download_url <- if (os_type == "Darwin") {
    "https://cloud.r-project.org/bin/macosx/"
  } else if (os_type == "Windows") {
    "https://cloud.r-project.org/bin/windows/base/"
  } else {
    "https://cloud.r-project.org/bin/linux/"
  }

  stop(paste0(
    "\n\n╔══════════════════════════════════════════════════════════╗\n",
    "║     CRITICAL: limpa package not found                    ║\n",
    "╚══════════════════════════════════════════════════════════╝\n\n",
    "Your R version: ", getRversion(), " (NEED: 4.5+)\n\n",
    "Upgrade R from: ", download_url, "\n",
    "Then run: BiocManager::install('limpa')\n\n"
  ))
}
library(limpa)

# Load app version from VERSION file
app_version <- tryCatch(
  trimws(readLines("VERSION", n = 1, warn = FALSE)),
  error = function(e) "unknown"
)

# Snapshot source file timestamps at startup for update detection
code_snapshot <- tryCatch({
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  all_files <- c("app.R", r_files)
  sum(file.mtime(all_files[file.exists(all_files)]), na.rm = TRUE)
}, error = function(e) 0)

# Load community stats — static JSON first, then overlay live GitHub API data
community_stats <- tryCatch({
  # Start with static JSON (has traffic trends + discussions from Actions workflow)
  stats_file <- file.path(getwd(), "stats", "community_stats.json")
  stats <- if (file.exists(stats_file)) jsonlite::fromJSON(stats_file) else list()

  # Overlay live stats from GitHub API (public, no auth needed)
  live <- tryCatch({
    repo_data <- jsonlite::fromJSON("https://api.github.com/repos/bsphinney/DE-LIMP")
    list(stars = repo_data$stargazers_count, forks = repo_data$forks_count)
  }, error = function(e) NULL)

  if (!is.null(live)) {
    if (is.null(stats$github)) stats$github <- list()
    stats$github$stars <- live$stars
    stats$github$forks <- live$forks
    stats$updated_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  }

  stats
}, error = function(e) NULL)

# Source R/ modules explicitly — ensures they load whether called via
# runApp('.') (auto-sources R/), runApp('app.R'), or Rscript app.R.
# Re-sourcing already-loaded functions is harmless (just redefines them).
local({
  r_dir <- "R"
  if (!dir.exists(r_dir)) {
    # Try relative to this script's location (e.g., Docker /srv/shiny-server/)
    script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
    r_dir <- file.path(script_dir, "R")
  }
  if (dir.exists(r_dir)) {
    for (f in sort(list.files(r_dir, pattern = "\\.R$", full.names = TRUE))) {
      source(f, local = FALSE)
    }
  }
})

# Clean up stale SSH sockets from previous sessions (prevents zombie mux blocking)
tryCatch(ssh_cleanup_stale_sockets(), error = function(e) NULL)

ui <- build_ui(is_hf_space, search_enabled, docker_available, hpc_available, local_sbatch,
               local_diann, delimp_data_dir,
               is_core_facility, cf_config, deploy_env)

# ==============================================================================
#  SERVER LOGIC — Thin orchestrator calling R/ modules
# ==============================================================================
server <- function(input, output, session) {

  # --- Shared reactive state ---
  values <- reactiveValues(
    raw_data = NULL, metadata = NULL, fit = NULL, y_protein = NULL,
    dpc_fit = NULL, status = "Waiting...", design = NULL, qc_stats = NULL,
    plot_selected_proteins = NULL, chat_history = list(),
    current_file_uri = NULL, gsea_results = NULL,
    gsea_results_cache = list(), gsea_last_contrast = NULL, gsea_last_org_db = NULL,
    repro_log = c(
      "# ==============================================================================",
      "# DE-LIMP Reproducibility Log",
      sprintf("# Session started: %s", Sys.time()),
      "# ==============================================================================",
      "",
      "# --- Load Required Libraries ---",
      "library(limpa); library(limma); library(dplyr); library(stringr); library(ggrepel);"
    ),
    color_plot_by_de = FALSE,
    grid_selected_protein = NULL,
    temp_violin_target = NULL,
    diann_norm_detected = "unknown",
    # XIC Viewer
    xic_dir = NULL, xic_available = FALSE, xic_format = "v2",
    xic_protein = NULL, xic_data = NULL, xic_report_map = NULL,
    uploaded_report_path = NULL, original_report_name = NULL,
    mobilogram_available = FALSE, mobilogram_files_found = 0,
    mobilogram_dir = NULL,
    # Phosphoproteomics
    phospho_detected = NULL,
    phospho_site_matrix = NULL,
    phospho_site_info = NULL,
    phospho_fit = NULL,
    phospho_site_matrix_filtered = NULL,
    phospho_input_mode = NULL,
    # Phospho Phase 2/3
    ksea_results = NULL,
    ksea_last_contrast = NULL,
    phospho_fasta_sequences = NULL,
    phospho_corrected_active = FALSE,
    phospho_annotations = NULL,
    # DIA-NN Search (HPC + Docker backends)
    diann_jobs = list(),
    diann_raw_files = NULL,
    diann_fasta_files = character(),
    diann_speclib = NULL,
    uniprot_results = NULL,
    fasta_info = NULL,
    library_locked = FALSE,
    ssh_connected = FALSE,
    ssh_sbatch_path = NULL,
    cluster_resources = NULL,
    public_resources = NULL,
    auto_partition = NULL,
    diann_search_settings = NULL,
    pending_notes_od = NULL,       # Set after search completes — triggers notes modal
    pending_notes_name = NULL,     # Search name for notes modal title
    instrument_metadata = NULL,    # List from parse_*_metadata() — instrument model, m/z range, etc.
    tic_traces = NULL,             # Named list of data.frames from extract_tic_timstof(), keyed by filename
    tic_metrics = NULL,            # data.frame: run, valid, total_auc, ..., shape_r, status, flags
    excluded_files = NULL,         # data.frame: filename, excluded_at, reason, user_note, source, group
    docker_available = docker_available,
    # Multi-View Integration (MOFA2)
    mofa_view_configs = list(),
    mofa_views = list(),
    mofa_view_fits = list(),
    mofa_sample_metadata = NULL,
    mofa_object = NULL,
    mofa_factors = NULL,
    mofa_weights = list(),
    mofa_variance_explained = NULL,
    mofa_last_run_params = NULL,
    # Run Comparator
    comparator_results          = NULL,
    comparator_run_a            = NULL,
    comparator_run_b            = NULL,
    comparator_mode             = NULL,
    comparator_gemini_narrative = NULL,
    comparator_mofa             = NULL,
    comparator_compare_from_history = NULL,
    comparator_diann_log_a      = NULL,
    comparator_diann_log_b      = NULL,
    per_user_resources          = NULL,
    # App metadata
    app_version = app_version,
    community_stats = community_stats
  )

  # --- Shared helper: append to reproducibility log ---
  add_to_log <- function(action_name, code_lines) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    header <- c("", paste0("# --- ", action_name, " [", timestamp, "] ---"))
    values$repro_log <- c(values$repro_log, header, code_lines)
  }

  # --- Call server modules (defined in R/ directory, auto-sourced by Shiny) ---
  server_data(input, output, session, values, add_to_log, is_hf_space)
  server_de(input, output, session, values, add_to_log)
  server_qc(input, output, session, values)
  server_viz(input, output, session, values, add_to_log, is_hf_space)
  server_gsea(input, output, session, values, add_to_log)
  server_ai(input, output, session, values)
  server_xic(input, output, session, values, is_hf_space)
  server_phospho(input, output, session, values, add_to_log)
  server_search(input, output, session, values, add_to_log,
                search_enabled, docker_available, docker_config, hpc_available, local_sbatch,
                local_diann, delimp_data_dir,
                is_core_facility, cf_config, local_sbatch_path)
  server_mofa(input, output, session, values, add_to_log)
  server_comparator(input, output, session, values, add_to_log)
  server_facility(input, output, session, values, add_to_log,
                  is_core_facility, cf_config, search_enabled)
  server_session(input, output, session, values, add_to_log)

  # --- Home directory quota check (HPC systems often have small quotas) ---
  session$onFlushed(function() {
    tryCatch({
      home <- Sys.getenv("HOME")
      if (!nzchar(home)) return()
      # Try writing a small temp file to detect quota issues
      test_file <- file.path(home, ".delimp_quota_test")
      writeLines("test", test_file)
      unlink(test_file)
    }, error = function(e) {
      if (grepl("quota|permission|read-only|no space", e$message, ignore.case = TRUE)) {
        showNotification(
          tagList(
            icon("exclamation-triangle"),
            tags$strong(" Home directory is full."),
            " Some features (job queue, activity log) may not work.",
            " Free up space in ", tags$code(Sys.getenv("HOME")),
            " or contact your HPC admin."
          ),
          type = "error", duration = NULL, id = "home_quota_warning"
        )
      }
    }, warning = function(w) {
      if (grepl("quota|no space", w$message, ignore.case = TRUE)) {
        showNotification(
          tagList(icon("exclamation-triangle"), " Home directory may be full: ", w$message),
          type = "warning", duration = 30, id = "home_quota_warning"
        )
      }
    })
  }, once = TRUE)

  # --- Code update detection: check if source files changed on disk ---
  observe({
    invalidateLater(30000)  # Check every 30 seconds
    current <- tryCatch({
      r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
      all_files <- c("app.R", r_files)
      sum(file.mtime(all_files[file.exists(all_files)]), na.rm = TRUE)
    }, error = function(e) code_snapshot)

    if (current != code_snapshot) {
      showNotification(
        tagList(
          icon("sync"), " Code update available. ",
          tags$a("Restart the app", href = "javascript:window.location.reload()",
                 style = "color: white; text-decoration: underline; font-weight: bold;"),
          " to apply."
        ),
        type = "warning", duration = NULL, id = "code_update_banner"
      )
    }
  })

  # --- Progressive reveal: hide result-dependent tabs until state exists ---
  session$onFlushed(once = TRUE, function() {
    nav_hide("main_tabs", "QC")
    nav_hide("main_tabs", "DE Dashboard")
    nav_hide("main_tabs", "Gene Set Enrichment")
    nav_hide("main_tabs", "AI Analysis")
    nav_hide("main_tabs", "Output")
    nav_hide("main_tabs", "Phosphoproteomics")
  })

  observe({
    if (!is.null(values$raw_data) || !is.null(values$tic_traces)) {
      nav_show("main_tabs", "QC")
    } else {
      nav_hide("main_tabs", "QC")
    }
  })

  observe({
    if (!is.null(values$fit)) {
      nav_show("main_tabs", "DE Dashboard")
      if (gsea_available) nav_show("main_tabs", "Gene Set Enrichment")
      nav_show("main_tabs", "AI Analysis")
      nav_show("main_tabs", "Output")
    } else if (!is.null(values$y_protein)) {
      # No DE (no replicates) but quantification done — show DE Dashboard for PCA
      nav_show("main_tabs", "DE Dashboard")
      nav_show("main_tabs", "Output")
      nav_hide("main_tabs", "Gene Set Enrichment")
      nav_hide("main_tabs", "AI Analysis")
    } else {
      nav_hide("main_tabs", "DE Dashboard")
      nav_hide("main_tabs", "Gene Set Enrichment")
      nav_hide("main_tabs", "AI Analysis")
      nav_hide("main_tabs", "Output")
    }
  })

  observe({
    phospho <- values$phospho_detected
    if (!is.null(phospho) && isTRUE(phospho$detected)) {
      nav_show("main_tabs", "Phosphoproteomics")
    } else {
      nav_hide("main_tabs", "Phosphoproteomics")
    }
  })
}


shinyApp(ui, server)
