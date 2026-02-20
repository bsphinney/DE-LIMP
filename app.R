# ==============================================================================
#  DE-LIMP: Differential Expression & Limpa Proteomics App
#  (Formerly LIMP-D)
#  Status: Production Ready (Hugging Face Compatible v1.2)
# ==============================================================================

# Set CRAN mirror to avoid interactive popup (especially in VS Code)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- 1. AUTO-INSTALLATION & SETUP ---
# IMPORTANT: Install packages BEFORE loading libraries to avoid conflicts

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("BiocManager not found. Installing...")
  install.packages("BiocManager", quiet = TRUE)
}

# Check for limpa and install if needed
if (!requireNamespace("limpa", quietly = TRUE)) {
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
required_pkgs <- c("shiny", "bslib", "readr", "tibble", "dplyr", "tidyr",
                   "ggplot2", "httr2", "rhandsontable", "DT", "arrow",
                   "ComplexHeatmap", "shinyjs", "plotly", "stringr", "limma",
                   "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
                   "enrichplot", "ggridges", "ggrepel", "markdown", "curl",
                   "KSEAapp", "ggseqlogo")

# Only install truly missing packages (don't update already-loaded packages)
missing_pkgs <- character(0)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, pkg)
  }
}

if (length(missing_pkgs) > 0) {
  message(paste0("Installing missing packages: ", paste(missing_pkgs, collapse = ", ")))
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE, quiet = TRUE)
}

# --- 2. SERVER CONFIGURATION ---
options(repos = c(BiocManager::repositories(), CRAN = "https://cloud.r-project.org"))

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
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(ggrepel)
library(markdown) # Needed for AI formatting

options(shiny.maxRequestSize = 5000 * 1024^2)  # 5 GB upload limit

# Detect Hugging Face Spaces environment (SPACE_ID is set automatically by HF)
is_hf_space <- nzchar(Sys.getenv("SPACE_ID", ""))

# Detect search backends (Docker local + HPC SSH/SLURM)
# Disabled on Hugging Face Spaces — search tab not useful in cloud environment
local_sbatch <- nzchar(Sys.which("sbatch"))
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

# Combined flag — at least one backend available
search_enabled <- docker_available || hpc_available

# Conditionally load search-related packages
if (search_enabled) {
  for (pkg in c("shinyFiles", "jsonlite")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, quiet = TRUE)
    }
  }
  library(shinyFiles)
  library(jsonlite)
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

ui <- build_ui(is_hf_space, search_enabled, docker_available, hpc_available, local_sbatch)

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
    ssh_connected = FALSE,
    diann_search_settings = NULL,
    docker_available = docker_available
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
                search_enabled, docker_available, docker_config, hpc_available, local_sbatch)
  server_session(input, output, session, values, add_to_log)
}


shinyApp(ui, server)
