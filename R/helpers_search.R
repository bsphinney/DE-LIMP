# helpers_search.R — Pure helper functions for DIA-NN Search Integration
# No Shiny reactivity. All functions are testable standalone.
# Supports both HPC (SSH/SLURM) and Local Docker backends.

# =============================================================================
# UniProt API Functions
# =============================================================================

#' Search UniProt for reference proteomes by organism name
#' @param query Character string — organism common or scientific name
#' @return data.frame with proteome ID, organism, protein count, type
search_uniprot_proteomes <- function(query) {
  url <- paste0(
    "https://rest.uniprot.org/proteomes/search?",
    "query=", utils::URLencode(paste0("(", query, ") AND (proteome_type:1)")),
    "&format=json",
    "&fields=upid,organism,organism_id,protein_count",
    "&size=25"
  )

  tryCatch({
    resp <- httr2::request(url) |>
      httr2::req_headers(Accept = "application/json") |>
      httr2::req_timeout(30) |>
      httr2::req_perform()

    data <- httr2::resp_body_json(resp)

    if (length(data$results) == 0) {
      return(data.frame(
        upid = character(), organism = character(),
        common_name = character(), taxonomy_id = integer(),
        protein_count = integer(), proteome_type = character(),
        stringsAsFactors = FALSE
      ))
    }

    data.frame(
      upid = vapply(data$results, function(r) r$id %||% "", character(1)),
      organism = vapply(data$results, function(r) {
        r$taxonomy$scientificName %||% ""
      }, character(1)),
      common_name = vapply(data$results, function(r) {
        r$taxonomy$commonName %||% ""
      }, character(1)),
      taxonomy_id = vapply(data$results, function(r) {
        as.integer(r$taxonomy$taxonId %||% 0L)
      }, integer(1)),
      protein_count = vapply(data$results, function(r) {
        as.integer(r$proteinCount %||% 0L)
      }, integer(1)),
      proteome_type = vapply(data$results, function(r) {
        pt <- r$proteomeType %||% ""
        if (grepl("Reference", pt, ignore.case = TRUE)) "Reference" else "Other"
      }, character(1)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(sprintf("[DE-LIMP Search] UniProt proteome search failed: %s", e$message))
    data.frame(
      upid = character(), organism = character(),
      common_name = character(), taxonomy_id = integer(),
      protein_count = integer(), proteome_type = character(),
      stringsAsFactors = FALSE
    )
  })
}

#' Download FASTA from UniProt for a given proteome
#' @param proteome_id Character — UniProt proteome ID (e.g., "UP000005640")
#' @param content_type Character — "one_per_gene", "reviewed", "full", "full_isoforms"
#' @param output_path Character — full path where FASTA will be saved
#' @return List with success status, path, sequence count, file size
download_uniprot_fasta <- function(proteome_id, content_type, output_path) {
  # Build query based on content type
  base_query <- sprintf("(proteome:%s)", proteome_id)
  query <- switch(content_type,
    "one_per_gene" = base_query,
    "reviewed"     = paste0(base_query, " AND (reviewed:true)"),
    "full"         = base_query,
    "full_isoforms" = base_query,
    base_query
  )

  include_isoform <- content_type == "full_isoforms"
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "query=", utils::URLencode(query),
    "&format=fasta",
    "&compressed=false",
    if (include_isoform) "&includeIsoform=true" else ""
  )

  tryCatch({
    # Download FASTA via httr2
    tmp_file <- tempfile(fileext = ".fasta")

    resp <- httr2::request(url) |>
      httr2::req_headers(Accept = "text/plain") |>
      httr2::req_timeout(300) |>
      httr2::req_perform(path = tmp_file)

    if (!file.exists(tmp_file) || file.size(tmp_file) < 100) {
      stop("Download failed or returned empty file")
    }

    # Count sequences
    n_seqs <- sum(grepl("^>", readLines(tmp_file, warn = FALSE)))

    # Move to final location
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    file.copy(tmp_file, output_path, overwrite = TRUE)
    unlink(tmp_file)

    list(
      success = TRUE,
      path = output_path,
      n_sequences = n_seqs,
      file_size = file.size(output_path),
      url = url
    )
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}

#' Get path to a bundled contaminant FASTA file
#' @param library_name Character — one of: "universal", "cell_culture", etc.
#' @param app_dir Character — app root directory (where contaminants/ lives)
#' @return List with success, path, n_sequences, file_size
get_contaminant_fasta <- function(library_name, app_dir = ".") {
  lib_map <- c(
    universal         = "Universal_Contaminants.fasta",
    cell_culture      = "Cell_Culture_Contaminants.fasta",
    mouse_tissue      = "Mouse_Tissue_Contaminants.fasta",
    rat_tissue        = "Rat_Tissue_Contaminants.fasta",
    neuron_culture    = "Neuron_Culture_Contaminants.fasta",
    stem_cell_culture = "Stem_Cell_Culture_Contaminants.fasta"
  )
  fname <- lib_map[[library_name]]
  if (is.null(fname)) return(list(success = FALSE, error = "Unknown library"))

  local_path <- file.path(app_dir, "contaminants", fname)
  if (!file.exists(local_path)) {
    return(list(success = FALSE, error = paste("File not found:", local_path)))
  }
  n_seqs <- sum(grepl("^>", readLines(local_path, warn = FALSE)))
  list(success = TRUE, path = local_path, n_sequences = n_seqs,
       file_size = file.size(local_path))
}

#' Generate a descriptive FASTA filename
generate_fasta_filename <- function(proteome_id, organism_name, content_type) {
  safe_org <- tolower(gsub("[^A-Za-z0-9]", "_", organism_name))
  safe_org <- gsub("_+", "_", safe_org)
  safe_org <- substr(safe_org, 1, 30)

  type_suffix <- switch(content_type,
    "one_per_gene"  = "opg",
    "reviewed"      = "sprot",
    "full"          = "full",
    "full_isoforms" = "full_iso",
    "custom"
  )

  release <- format(Sys.Date(), "%Y_%m")
  sprintf("%s_%s_%s_%s.fasta", proteome_id, safe_org, type_suffix, release)
}

# =============================================================================
# File Discovery Functions
# =============================================================================

#' Scan a directory for MS raw data files
#' @param dir_path Character — path to scan
#' @return data.frame with filename, size_mb, type columns
scan_raw_files <- function(dir_path) {
  if (!dir.exists(dir_path)) {
    return(data.frame(filename = character(), size_mb = numeric(),
                      type = character(), stringsAsFactors = FALSE))
  }

  # .d directories are Bruker raw data (special handling)
  d_dirs <- list.dirs(dir_path, recursive = FALSE, full.names = TRUE)
  d_dirs <- d_dirs[grepl("\\.d$", d_dirs, ignore.case = TRUE)]

  # .raw and .mzML are regular files
  raw_files <- list.files(dir_path, pattern = "\\.(raw|mzML)$",
                          ignore.case = TRUE, full.names = TRUE)

  all_files <- c(d_dirs, raw_files)

  if (length(all_files) == 0) {
    return(data.frame(filename = character(), size_mb = numeric(),
                      type = character(), stringsAsFactors = FALSE))
  }

  # Get sizes (for .d dirs, sum contents)
  sizes <- vapply(all_files, function(f) {
    if (dir.exists(f)) {
      files_in <- list.files(f, recursive = TRUE, full.names = TRUE)
      sum(file.size(files_in), na.rm = TRUE) / 1e6
    } else {
      file.size(f) / 1e6
    }
  }, numeric(1))

  types <- vapply(all_files, function(f) {
    if (dir.exists(f) && grepl("\\.d$", f, ignore.case = TRUE)) return(".d")
    ext <- tools::file_ext(f)
    paste0(".", tolower(ext))
  }, character(1))

  data.frame(
    filename = basename(all_files),
    full_path = all_files,
    size_mb = round(sizes, 1),
    type = types,
    stringsAsFactors = FALSE
  )
}

#' Scan a directory for pre-staged FASTA databases
#' @param fasta_dir Character — path to scan
#' @return Named character vector suitable for selectInput choices
scan_prestaged_databases <- function(fasta_dir) {
  if (!dir.exists(fasta_dir)) return(character())

  fasta_files <- list.files(fasta_dir, pattern = "\\.(fasta|fa)$",
                            ignore.case = TRUE, full.names = TRUE)
  if (length(fasta_files) == 0) return(character())

  # Build display names from filenames
  display_names <- vapply(fasta_files, function(f) {
    bn <- basename(f)
    size_mb <- round(file.size(f) / 1e6, 1)
    sprintf("%s (%s MB)", bn, size_mb)
  }, character(1))

  stats::setNames(fasta_files, display_names)
}

# =============================================================================
# DIA-NN Flag Building (shared by HPC and Docker backends)
# =============================================================================

#' Build DIA-NN CLI flags from search parameters
#' Returns a character vector of flags (without --f, --fasta, --out, --threads).
#' Used by both generate_sbatch_script() and build_docker_command().
#' @param search_params List of search parameters (qvalue, enzyme, mods, etc.)
#' @param search_mode Character: "libfree", "library", or "phospho"
#' @param normalization Character: "on" or "off"
#' @param speclib_mount Character or NULL: container-internal path to spectral library
#' @return Character vector of DIA-NN CLI flags
build_diann_flags <- function(search_params = list(), search_mode = "libfree",
                              normalization = "on", speclib_mount = NULL) {
  # Defaults for search params
  sp <- list(
    qvalue = 0.01, max_var_mods = 1, scan_window = 6,
    mass_acc_mode = "auto", mass_acc = 14, mass_acc_ms1 = 14,
    unimod4 = TRUE, met_excision = TRUE,
    min_pep_len = 7, max_pep_len = 30,
    min_pr_mz = 300, max_pr_mz = 1200,
    min_pr_charge = 1, max_pr_charge = 4,
    min_fr_mz = 200, max_fr_mz = 1200,
    enzyme = "K*,R*", missed_cleavages = 1,
    mbr = TRUE, rt_profiling = TRUE, xic = TRUE,
    mod_met_ox = TRUE, mod_nterm_acetyl = FALSE,
    extra_var_mods = "", extra_cli_flags = ""
  )
  for (nm in names(search_params)) sp[[nm]] <- search_params[[nm]]

  flags <- c()
  is_phospho <- identical(search_mode, "phospho")

  # Variable modification flags
  if (isTRUE(sp$mod_met_ox)) flags <- c(flags, "--var-mod UniMod:35,15.994915,M")
  if (isTRUE(sp$mod_nterm_acetyl)) flags <- c(flags, "--var-mod UniMod:1,42.010565,*n")
  if (nzchar(sp$extra_var_mods)) {
    extra_lines <- trimws(strsplit(sp$extra_var_mods, "\n")[[1]])
    for (mod in extra_lines) {
      if (nzchar(mod)) flags <- c(flags, sprintf("--var-mod %s", mod))
    }
  }

  # Core shared flags
  flags <- c(flags,
    "--out-lib /work/out/report-lib.parquet",
    "--matrices",
    "--gen-spec-lib",
    sprintf("--qvalue %s", sp$qvalue),
    "--verbose 1",
    sprintf("--var-mods %d", sp$max_var_mods)
  )

  if (isTRUE(sp$xic)) flags <- c(flags, "--xic")
  if (isTRUE(sp$unimod4)) flags <- c(flags, "--unimod4")

  # Library mode
  if (search_mode == "library" && !is.null(speclib_mount)) {
    flags <- c(flags,
      sprintf("--lib %s", speclib_mount),
      sprintf("--window %d", sp$scan_window),
      "--use-quant"
    )
  }

  # Library-free mode (and phospho)
  if (search_mode != "library") {
    flags <- c(flags,
      "--fasta-search",
      "--predictor",
      sprintf("--cut %s", sp$enzyme),
      sprintf("--missed-cleavages %d", sp$missed_cleavages),
      sprintf("--min-pep-len %d", sp$min_pep_len),
      sprintf("--max-pep-len %d", sp$max_pep_len),
      sprintf("--min-pr-mz %d", sp$min_pr_mz),
      sprintf("--max-pr-mz %d", sp$max_pr_mz),
      sprintf("--min-pr-charge %d", sp$min_pr_charge),
      sprintf("--max-pr-charge %d", sp$max_pr_charge),
      sprintf("--min-fr-mz %d", sp$min_fr_mz),
      sprintf("--max-fr-mz %d", sp$max_fr_mz)
    )
    if (isTRUE(sp$met_excision)) flags <- c(flags, "--met-excision")
  }

  # Mass accuracy
  if (sp$mass_acc_mode == "manual") {
    flags <- c(flags,
      sprintf("--window %d", sp$scan_window),
      sprintf("--mass-acc %s", sp$mass_acc),
      sprintf("--mass-acc-ms1 %s", sp$mass_acc_ms1)
    )
  }

  # Toggles
  if (isTRUE(sp$mbr)) flags <- c(flags, "--reanalyse")
  if (isTRUE(sp$rt_profiling)) flags <- c(flags, "--rt-profiling")
  if (normalization == "off") flags <- c(flags, "--no-norm")

  # Phospho-specific
  if (is_phospho) {
    flags <- c(flags, "--phospho-output", "--report-lib-info")
  }

  # Extra CLI flags
  if (nzchar(sp$extra_cli_flags)) {
    flags <- c(flags, trimws(sp$extra_cli_flags))
  }

  flags
}

# =============================================================================
# sbatch Script Generation (HPC backend)
# =============================================================================

#' Generate a complete sbatch script for DIA-NN search
#' @return Character string: complete sbatch script content
generate_sbatch_script <- function(
  analysis_name, raw_files, fasta_files, speclib_path = NULL,
  output_dir, diann_sif, normalization = "on", search_mode = "libfree",
  cpus = 64, mem_gb = 512, time_hours = 12,
  partition = "high", account = "genome-center-grp",
  search_params = list()
) {
  # Determine output filename
  report_name <- if (normalization == "off") "no_norm_report.parquet" else "report.parquet"

  # Determine unique directories for data and fasta
  data_dirs <- unique(dirname(raw_files))
  fasta_dirs <- unique(dirname(fasta_files))

  # Build bind mount string — handle multiple FASTA directories
  fasta_bind_parts <- if (length(fasta_dirs) == 1) {
    sprintf("%s:/work/fasta", fasta_dirs[1])
  } else {
    sprintf("%s:/work/fasta%d", fasta_dirs, seq_along(fasta_dirs))
  }
  bind_parts <- c(
    sprintf("%s:/work/data", data_dirs[1]),
    fasta_bind_parts,
    sprintf("%s:/work/out", output_dir)
  )
  if (!is.null(speclib_path) && nzchar(speclib_path)) {
    bind_parts <- c(bind_parts, sprintf("%s:/work/lib", dirname(speclib_path)))
  }
  bind_mount <- paste(bind_parts, collapse = ",")

  # Build --f flags for raw files
  run_flags <- paste(sprintf("    --f /work/data/%s", basename(raw_files)),
                     collapse = " \\\n")

  # Build --fasta flags — map each file to its mount point
  fasta_mount_map <- if (length(fasta_dirs) == 1) {
    rep("/work/fasta", length(fasta_files))
  } else {
    sprintf("/work/fasta%d", match(dirname(fasta_files), fasta_dirs))
  }
  fasta_flags <- paste(sprintf("    --fasta %s/%s", fasta_mount_map, basename(fasta_files)),
                       collapse = " \\\n")

  # Get shared DIA-NN flags via build_diann_flags()
  speclib_mount <- if (!is.null(speclib_path) && nzchar(speclib_path)) {
    sprintf("/work/lib/%s", basename(speclib_path))
  } else NULL
  shared_flags <- build_diann_flags(search_params, search_mode, normalization, speclib_mount)

  # Build DIA-NN command for apptainer
  diann_cmd_parts <- c(
    sprintf("apptainer exec --bind %s %s /diann-2.3.0/diann-linux \\", bind_mount, diann_sif),
    paste0(run_flags, " \\"),
    paste0(fasta_flags, " \\"),
    sprintf("    --out /work/out/%s \\", report_name),
    sprintf("    --threads %d \\", cpus),
    paste0("    ", shared_flags)
  )

  # Remove NULLs and trailing backslash on last line
  diann_cmd_parts <- Filter(Negate(is.null), diann_cmd_parts)
  # Add line continuations to all but last flag line
  for (i in seq_along(diann_cmd_parts)) {
    if (i < length(diann_cmd_parts) && !grepl(" \\\\$", diann_cmd_parts[i])) {
      diann_cmd_parts[i] <- paste0(diann_cmd_parts[i], " \\")
    }
  }
  # Ensure last line has no trailing backslash
  last <- length(diann_cmd_parts)
  diann_cmd_parts[last] <- sub(" \\\\$", "", diann_cmd_parts[last])
  diann_cmd <- paste(diann_cmd_parts, collapse = "\n")

  # Assemble full sbatch script
  script <- paste0(
    '#!/bin/bash -l\n',
    sprintf('#SBATCH --job-name=diann_%s\n', analysis_name),
    sprintf('#SBATCH --cpus-per-task=%d\n', cpus),
    sprintf('#SBATCH --mem=%dG\n', mem_gb),
    sprintf('#SBATCH -o %s/diann_%%j.out\n', output_dir),
    sprintf('#SBATCH -e %s/diann_%%j.err\n', output_dir),
    sprintf('#SBATCH --account=%s\n', account),
    sprintf('#SBATCH --time=%d:00:00\n', time_hours),
    sprintf('#SBATCH --partition=%s\n', partition),
    '\n',
    'module load apptainer\n',
    '\n',
    sprintf('echo "DIA-NN search: %s"\n', analysis_name),
    'echo "Started: $(date)"\n',
    sprintf('echo "Output: %s"\n', output_dir),
    '\n',
    diann_cmd, '\n',
    '\n',
    'EXIT_CODE=$?\n',
    'echo ""\n',
    'echo "DIA-NN finished with exit code: $EXIT_CODE"\n',
    'echo "Completed: $(date)"\n',
    'exit $EXIT_CODE\n'
  )

  return(script)
}

# =============================================================================
# Docker Helper Functions (Local backend)
# =============================================================================

#' Check if Docker is installed and daemon is running
#' @return list(available, daemon_running, error)
check_docker_available <- function() {
  if (!nzchar(Sys.which("docker"))) {
    return(list(available = FALSE, daemon_running = FALSE,
                error = "Docker CLI not found on PATH"))
  }
  daemon_ok <- tryCatch({
    out <- system2("docker", "info", stdout = TRUE, stderr = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(e) FALSE)

  list(available = TRUE, daemon_running = daemon_ok,
       error = if (!daemon_ok) "Docker daemon not running" else NULL)
}

#' Check if a DIA-NN Docker image exists locally
#' @param image_name Character — Docker image name (e.g., "diann:2.3.0")
#' @return list(exists, image_name, error)
check_diann_image <- function(image_name = "diann:2.3.0") {
  exists <- tryCatch({
    out <- system2("docker", c("image", "inspect", image_name),
                   stdout = TRUE, stderr = TRUE)
    TRUE
  }, error = function(e) FALSE, warning = function(e) FALSE)

  list(exists = exists, image_name = image_name,
       error = if (!exists) paste("Image not found:", image_name) else NULL)
}

#' Detect host machine CPU and memory resources
#' @return list(cpus, memory_gb)
get_host_resources <- function() {
  cpus <- tryCatch(parallel::detectCores(), error = function(e) 4L)
  if (is.na(cpus)) cpus <- 4L

  mem_gb <- tryCatch({
    os <- Sys.info()[["sysname"]]
    if (os == "Darwin") {
      # macOS: sysctl hw.memsize returns bytes
      raw <- system2("sysctl", c("-n", "hw.memsize"), stdout = TRUE, stderr = TRUE)
      as.integer(as.numeric(raw) / 1024^3)
    } else if (os == "Linux") {
      raw <- readLines("/proc/meminfo", n = 1)
      kb <- as.numeric(gsub("[^0-9]", "", raw))
      as.integer(kb / 1024^2)
    } else {
      # Windows or unknown
      64L
    }
  }, error = function(e) 64L)

  list(cpus = cpus, memory_gb = mem_gb)
}

#' Build Docker command arguments for running DIA-NN locally
#' @param raw_files Character vector — local paths to raw data files/dirs
#' @param fasta_files Character vector — local paths to FASTA files
#' @param output_dir Character — local output directory
#' @param image_name Character — Docker image name
#' @param diann_flags Character vector — flags from build_diann_flags()
#' @param cpus Integer — CPU limit
#' @param mem_gb Integer — memory limit (GB)
#' @param container_name Character — name for the container
#' @param speclib_path Character or NULL — local path to spectral library
#' @param report_name Character — output report filename
#' @return Character vector suitable for system2("docker", args)
build_docker_command <- function(raw_files, fasta_files, output_dir, image_name,
                                 diann_flags, cpus, mem_gb, container_name,
                                 speclib_path = NULL, report_name = "report.parquet") {
  # Identify unique data and fasta directories
  data_dirs <- unique(dirname(raw_files))
  fasta_dirs <- unique(dirname(fasta_files))

  # Build volume mounts
  volumes <- c()

  # Data directory mount (read-only)
  if (length(data_dirs) == 1) {
    volumes <- c(volumes, "-v", sprintf("%s:/work/data:ro", data_dirs[1]))
  } else {
    for (i in seq_along(data_dirs)) {
      volumes <- c(volumes, "-v", sprintf("%s:/work/data%d:ro", data_dirs[i], i))
    }
  }

  # FASTA directory mount(s) (read-only)
  if (length(fasta_dirs) == 1) {
    volumes <- c(volumes, "-v", sprintf("%s:/work/fasta:ro", fasta_dirs[1]))
  } else {
    for (i in seq_along(fasta_dirs)) {
      volumes <- c(volumes, "-v", sprintf("%s:/work/fasta%d:ro", fasta_dirs[i], i))
    }
  }

  # Output directory mount (read-write)
  volumes <- c(volumes, "-v", sprintf("%s:/work/out", output_dir))

  # Spectral library mount
  if (!is.null(speclib_path) && nzchar(speclib_path)) {
    volumes <- c(volumes, "-v", sprintf("%s:/work/lib:ro", dirname(speclib_path)))
  }

  # Build --f flags for raw files (mapped to container paths)
  f_flags <- c()
  if (length(data_dirs) == 1) {
    f_flags <- sprintf("--f /work/data/%s", basename(raw_files))
  } else {
    data_map <- match(dirname(raw_files), data_dirs)
    f_flags <- sprintf("--f /work/data%d/%s", data_map, basename(raw_files))
  }

  # Build --fasta flags (mapped to container paths)
  fasta_flags <- c()
  if (length(fasta_dirs) == 1) {
    fasta_flags <- sprintf("--fasta /work/fasta/%s", basename(fasta_files))
  } else {
    fasta_map <- match(dirname(fasta_files), fasta_dirs)
    fasta_flags <- sprintf("--fasta /work/fasta%d/%s", fasta_map, basename(fasta_files))
  }

  # Build the DIA-NN command that runs inside the container.
  # CRITICAL: DIA-NN writes large intermediate files (.predicted.speclib,
  # .quant files) to the same directory as --out. On Windows Docker Desktop,
  # the FUSE layer for bind-mounted volumes can't handle multi-GB writes,
  # causing "Could not save" errors. Fix: run DIA-NN with output to
  # container-internal /tmp, then copy only the final reports to /work/out.
  # Redirect --out-lib from /work/out/ to /tmp/diann_work/
  diann_flags_local <- gsub("--out-lib /work/out/", "--out-lib /tmp/diann_work/", diann_flags)

  diann_shell_cmd <- paste0(
    "mkdir -p /tmp/diann_work && ",
    paste(c(
      "diann-linux",
      f_flags,
      fasta_flags,
      sprintf("--out /tmp/diann_work/%s", report_name),
      sprintf("--threads %d", cpus),
      diann_flags_local
    ), collapse = " "),
    # Copy final outputs to the mounted volume (semicolons: run all even if some globs miss)
    " && { cp /tmp/diann_work/*.parquet /work/out/ 2>/dev/null;",
    " cp /tmp/diann_work/*.tsv /work/out/ 2>/dev/null;",
    " cp -r /tmp/diann_work/*_xic /work/out/ 2>/dev/null;",
    " true; }"
  )

  # Assemble full docker run command args
  args <- c(
    "run", "--rm", "-d",
    "--platform", "linux/amd64",
    "--name", container_name,
    sprintf("--cpus=%d", cpus),
    sprintf("--memory=%dg", mem_gb),
    volumes,
    "--entrypoint", "sh",
    image_name,
    "-c", diann_shell_cmd
  )

  args
}

#' Check Docker container status
#' @param container_id Character — Docker container ID or name
#' @return list(status, exit_code, log_tail)
check_docker_container_status <- function(container_id) {
  # Get container state
  state <- tryCatch({
    out <- system2("docker", c("inspect", "--format", "{{.State.Status}}",
                               container_id), stdout = TRUE, stderr = TRUE)
    trimws(out[1])
  }, error = function(e) "unknown",
     warning = function(e) "unknown")

  exit_code <- NA_integer_
  if (state == "exited") {
    exit_code <- tryCatch({
      out <- system2("docker", c("inspect", "--format", "{{.State.ExitCode}}",
                                 container_id), stdout = TRUE, stderr = TRUE)
      as.integer(trimws(out[1]))
    }, error = function(e) NA_integer_,
       warning = function(e) NA_integer_)
  }

  # Map Docker state to DE-LIMP job status
  status <- switch(state,
    "running"    = "running",
    "created"    = "queued",
    "exited"     = if (!is.na(exit_code) && exit_code == 0) "completed" else "failed",
    "dead"       = "failed",
    "removing"   = "running",
    "unknown"
  )

  # Tail logs
  log_tail <- tryCatch({
    out <- system2("docker", c("logs", "--tail", "30", container_id),
                   stdout = TRUE, stderr = TRUE)
    paste(out, collapse = "\n")
  }, error = function(e) "",
     warning = function(e) "")

  list(status = status, exit_code = exit_code, log_tail = log_tail)
}

# =============================================================================
# Local DIA-NN Execution (embedded binary — no Docker/SLURM)
# =============================================================================

#' Launch DIA-NN as a background process via processx
#' @param raw_files Character vector of raw file paths
#' @param fasta_files Character vector of FASTA file paths
#' @param output_dir Output directory path
#' @param diann_flags Character vector of DIA-NN CLI flags (from build_diann_flags)
#' @param threads Number of threads
#' @param log_file Path to write stdout+stderr log
#' @param speclib_path Optional spectral library path
#' @param report_name Output report filename (default: report.parquet)
#' @return list(process, pid, log_file)
run_local_diann <- function(raw_files, fasta_files, output_dir,
                             diann_flags, threads, log_file,
                             speclib_path = NULL, report_name = "report.parquet") {
  diann_bin <- Sys.which("diann")
  if (!nzchar(diann_bin)) diann_bin <- Sys.which("diann-linux")
  if (!nzchar(diann_bin)) stop("DIA-NN binary not found on PATH")

  # Build argument vector — each flag is a separate element
  args <- c()
  for (f in raw_files) args <- c(args, "--f", f)
  for (f in fasta_files) args <- c(args, "--fasta", f)
  args <- c(args, "--out", file.path(output_dir, report_name))
  args <- c(args, "--threads", as.character(threads))

  if (!is.null(speclib_path) && nzchar(speclib_path)) {
    args <- c(args, "--lib", speclib_path)
  }

  # Add DIA-NN flags (each may be "--flag value" — split on first space)
  for (flag in diann_flags) {
    parts <- strsplit(flag, " ", fixed = TRUE)[[1]]
    args <- c(args, parts)
  }

  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Launch as background process
  proc <- processx::process$new(
    command = diann_bin,
    args = args,
    stdout = log_file,
    stderr = log_file,
    cleanup_tree = TRUE
  )

  list(process = proc, pid = proc$get_pid(), log_file = log_file)
}

#' Check status of a locally running DIA-NN process
#' @param proc processx::process object
#' @param log_file Path to the log file
#' @return list(status, exit_code, log_tail)
check_local_diann_status <- function(proc, log_file) {
  alive <- tryCatch(proc$is_alive(), error = function(e) FALSE)
  exit_code <- if (!alive) {
    tryCatch(proc$get_exit_status(), error = function(e) NA_integer_)
  } else NA_integer_

  status <- if (alive) "running"
            else if (!is.na(exit_code) && exit_code == 0) "completed"
            else if (!is.na(exit_code)) "failed"
            else "unknown"

  log_tail <- tryCatch({
    lines <- readLines(log_file, warn = FALSE)
    iconv(paste(tail(lines, 30), collapse = "\n"), from = "", to = "UTF-8", sub = "")
  }, error = function(e) "")

  list(status = status, exit_code = exit_code, log_tail = log_tail)
}

# =============================================================================
# Job Recovery Functions
# =============================================================================

#' Recover DIA-NN jobs from SLURM accounting (sacct)
#' @param ssh_config SSH config list, or NULL for local
#' @param sbatch_path Full path to sbatch binary (used to find sacct)
#' @param days_back Integer — how many days back to search (default 7)
#' @return data.frame with job_id, name, state, elapsed, or empty data.frame
recover_slurm_jobs <- function(ssh_config = NULL, sbatch_path = NULL, days_back = 7) {
  empty <- data.frame(job_id = character(), name = character(),
                      state = character(), elapsed = character(),
                      stringsAsFactors = FALSE)

  # Build sacct path from sbatch path if available
  sacct_bin <- if (!is.null(sbatch_path) && nzchar(sbatch_path)) {
    file.path(dirname(sbatch_path), "sacct")
  } else {
    "sacct"
  }

  # Query sacct for recent jobs with "diann" in the name
  # Include WorkDir so we can find log files and output
  cmd <- paste0(
    sacct_bin,
    " --starttime=$(date -d '", days_back, " days ago' +%Y-%m-%d 2>/dev/null || ",
    "date -v-", days_back, "d +%Y-%m-%d)",
    " --format=JobID%20,JobName%50,State%20,Elapsed%15,WorkDir%120",
    " --parsable2 --noheader",
    " 2>/dev/null | grep -i diann | grep -v '\\.' "
  )

  result <- if (!is.null(ssh_config)) {
    ssh_exec(ssh_config, cmd, login_shell = is.null(sbatch_path))
  } else if (slurm_proxy_available()) {
    slurm_proxy_exec(cmd, timeout = 30)
  } else {
    tryCatch({
      stdout <- system2("bash", args = c("-c", cmd), stdout = TRUE, stderr = TRUE)
      list(status = 0, stdout = stdout)
    }, error = function(e) list(status = 1, stdout = character()))
  }

  if (result$status != 0 || length(result$stdout) == 0) return(empty)

  lines <- result$stdout[nzchar(result$stdout)]
  if (length(lines) == 0) return(empty)

  parsed <- strsplit(lines, "\\|")
  parsed <- parsed[vapply(parsed, length, integer(1)) >= 4]
  if (length(parsed) == 0) return(empty)

  df <- data.frame(
    job_id = vapply(parsed, `[`, character(1), 1),
    name = trimws(vapply(parsed, `[`, character(1), 2)),
    state = trimws(vapply(parsed, `[`, character(1), 3)),
    elapsed = trimws(vapply(parsed, `[`, character(1), 4)),
    stringsAsFactors = FALSE
  )

  # WorkDir is field 5 (may be missing for some jobs)
  df$work_dir <- vapply(parsed, function(p) {
    if (length(p) >= 5) trimws(p[5]) else ""
  }, character(1))

  df
}

#' Recover DIA-NN jobs from Docker containers
#' @return data.frame with container_id, name, state, created
recover_docker_jobs <- function() {
  empty <- data.frame(container_id = character(), name = character(),
                      state = character(), created = character(),
                      stringsAsFactors = FALSE)

  result <- tryCatch({
    out <- system2("docker", c("ps", "-a",
      "--filter", "name=delimp_",
      "--format", "{{.ID}}\t{{.Names}}\t{{.Status}}\t{{.CreatedAt}}"),
      stdout = TRUE, stderr = TRUE)
    list(status = 0, stdout = out)
  }, error = function(e) list(status = 1, stdout = character()))

  if (result$status != 0 || length(result$stdout) == 0) return(empty)

  lines <- result$stdout[nzchar(result$stdout)]
  if (length(lines) == 0) return(empty)

  parsed <- strsplit(lines, "\t")
  parsed <- parsed[vapply(parsed, length, integer(1)) >= 4]
  if (length(parsed) == 0) return(empty)

  data.frame(
    container_id = vapply(parsed, `[`, character(1), 1),
    name = vapply(parsed, `[`, character(1), 2),
    state = vapply(parsed, `[`, character(1), 3),
    created = vapply(parsed, `[`, character(1), 4),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# SSH Helper Functions
# =============================================================================

#' Execute a command on a remote host via SSH
#' @param ssh_config list(host, user, port, key_path) or NULL for local
#' @param command Character — command to execute remotely
#' @return list(status, stdout) — status is exit code, stdout is character vector
ssh_exec <- function(ssh_config, command, login_shell = FALSE, timeout = 60) {
  # Optionally wrap in login shell so .bash_profile / module paths are loaded
  # Prepend module loads if specified
  if (login_shell) {
    modules <- ssh_config$modules %||% ""
    mod_cmd <- if (nzchar(modules)) {
      mod_names <- trimws(strsplit(modules, "[,;[:space:]]+")[[1]])
      mod_names <- mod_names[nzchar(mod_names)]
      if (length(mod_names) > 0) {
        paste0(paste("module load", mod_names, "2>/dev/null;"), collapse = " ")
      } else ""
    } else ""
    full_cmd <- if (nzchar(mod_cmd)) paste(mod_cmd, command) else command
    remote_cmd <- paste0("bash -l -c ", shQuote(full_cmd))
  } else {
    remote_cmd <- command
  }
  args <- c(
    "-i", ssh_config$key_path,
    "-p", as.character(ssh_config$port %||% 22),
    "-o", "StrictHostKeyChecking=accept-new",
    "-o", "ConnectTimeout=10",
    "-o", "ServerAliveInterval=5",
    "-o", "ServerAliveCountMax=6",
    "-o", "BatchMode=yes",
    paste0(ssh_config$user, "@", ssh_config$host),
    remote_cmd
  )
  # Use processx for timeout support if available, else system2
  # Suppress macOS ARM64 MallocStackLogging warnings via environment
  stdout <- tryCatch({
    if (requireNamespace("processx", quietly = TRUE)) {
      res <- processx::run("ssh", args = args, timeout = timeout,
                           error_on_status = FALSE,
                           env = c("current", MallocStackLogging = "0"))
      out <- strsplit(res$stdout, "\n")[[1]]
      if (res$status != 0) attr(out, "status") <- res$status
      out
    } else {
      system2("ssh", args = args, stdout = TRUE, stderr = TRUE)
    }
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("timeout", msg, ignore.case = TRUE)) {
      msg <- paste("Command timed out after", timeout, "seconds")
    }
    structure(msg, status = 124L)
  })
  status <- attr(stdout, "status") %||% 0L
  # Sanitize output to valid UTF-8 (SSH may return ANSI codes, MOTD banners, etc.)
  stdout <- iconv(stdout, from = "", to = "UTF-8", sub = "")
  list(status = status, stdout = stdout)
}

#' Download a file from remote host via SCP
#' @param ssh_config list(host, user, port, key_path)
#' @param remote_path Character — full path on remote
#' @param local_path Character — full path on local machine
#' @return list(status, stdout)
scp_download <- function(ssh_config, remote_path, local_path) {
  args <- c(
    "-i", ssh_config$key_path,
    "-P", as.character(ssh_config$port %||% 22),
    "-o", "StrictHostKeyChecking=accept-new",
    "-o", "BatchMode=yes",
    paste0(ssh_config$user, "@", ssh_config$host, ":", remote_path),
    local_path
  )
  stdout <- tryCatch(
    system2("scp", args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      structure(conditionMessage(e), status = 1L)
    }
  )
  status <- attr(stdout, "status") %||% 0L
  # Sanitize output to valid UTF-8 (SSH may return ANSI codes, MOTD banners, etc.)
  stdout <- iconv(stdout, from = "", to = "UTF-8", sub = "")
  list(status = status, stdout = stdout)
}

#' Upload a local file to remote host via SCP
#' @param ssh_config list(host, user, port, key_path)
#' @param local_path Character — full path on local machine
#' @param remote_path Character — full path on remote
#' @return list(status, stdout)
scp_upload <- function(ssh_config, local_path, remote_path) {
  args <- c(
    "-i", ssh_config$key_path,
    "-P", as.character(ssh_config$port %||% 22),
    "-o", "StrictHostKeyChecking=accept-new",
    "-o", "BatchMode=yes",
    local_path,
    paste0(ssh_config$user, "@", ssh_config$host, ":", remote_path)
  )
  stdout <- tryCatch(
    system2("scp", args = args, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      structure(conditionMessage(e), status = 1L)
    }
  )
  status <- attr(stdout, "status") %||% 0L
  # Sanitize output to valid UTF-8 (SSH may return ANSI codes, MOTD banners, etc.)
  stdout <- iconv(stdout, from = "", to = "UTF-8", sub = "")
  list(status = status, stdout = stdout)
}

#' Test SSH connection and verify sbatch is available
#' @param ssh_config list(host, user, port, key_path)
#' @return list(success, message, sbatch_path)
test_ssh_connection <- function(ssh_config) {
  if (is.null(ssh_config$host) || !nzchar(ssh_config$host)) {
    return(list(success = FALSE, message = "No hostname specified", sbatch_path = NULL))
  }
  if (!file.exists(ssh_config$key_path %||% "")) {
    return(list(success = FALSE,
                message = paste("SSH key not found:", ssh_config$key_path),
                sbatch_path = NULL))
  }

  # Step 1: Test basic SSH connectivity (no login shell wrapper)
  result <- ssh_exec(ssh_config, "echo SSH_OK", login_shell = FALSE)
  if (!any(grepl("SSH_OK", result$stdout))) {
    msg <- paste(result$stdout, collapse = " ")
    if (!nzchar(msg)) msg <- paste("Exit code", result$status)
    return(list(success = FALSE,
                message = paste("SSH connection failed:", msg),
                sbatch_path = NULL))
  }

  # Step 2: Probe for sbatch — try multiple approaches
  sbatch_path <- NULL

  # Try 1: login shell with modules
  result2 <- ssh_exec(ssh_config, "which sbatch 2>/dev/null",
                       login_shell = TRUE)
  sbatch_line <- grep("^/", result2$stdout, value = TRUE)
  if (length(sbatch_line) > 0) sbatch_path <- sbatch_line[1]

  # Try 2: common HPC paths
  if (is.null(sbatch_path)) {
    result3 <- ssh_exec(ssh_config,
      "for p in /usr/bin/sbatch /usr/local/bin/sbatch /opt/slurm/bin/sbatch /cm/shared/apps/slurm/current/bin/sbatch; do [ -x \"$p\" ] && echo \"$p\" && break; done",
      login_shell = FALSE)
    sbatch_line <- grep("^/", result3$stdout, value = TRUE)
    if (length(sbatch_line) > 0) sbatch_path <- sbatch_line[1]
  }

  # Try 3: use locate or find
  if (is.null(sbatch_path)) {
    result4 <- ssh_exec(ssh_config,
      "command -v sbatch 2>/dev/null || type -P sbatch 2>/dev/null",
      login_shell = FALSE)
    sbatch_line <- grep("^/", result4$stdout, value = TRUE)
    if (length(sbatch_line) > 0) sbatch_path <- sbatch_line[1]
  }

  if (is.null(sbatch_path)) {
    return(list(success = TRUE,
                message = paste0("Connected to ", ssh_config$host,
                                 " but sbatch not found. Check 'Modules to Load' or contact HPC admin."),
                sbatch_path = NULL))
  }

  list(success = TRUE,
       message = sprintf("Connected to %s (sbatch: %s)", ssh_config$host, sbatch_path),
       sbatch_path = sbatch_path)
}

#' Scan raw files on a remote host via SSH
#' @param ssh_config list(host, user, port, key_path)
#' @param dir_path Character — remote directory path
#' @return data.frame(filename, size_mb, type) or empty data.frame
ssh_scan_raw_files <- function(ssh_config, dir_path) {
  empty_df <- data.frame(filename = character(), size_mb = numeric(),
                         type = character(), stringsAsFactors = FALSE)

  # du -sm with globs — no recursion into .d directories
  # Quote the directory path (may contain spaces) but leave glob unquoted for expansion
  qdir <- shQuote(dir_path)
  cmd <- paste0(
    "du -sm ", qdir, "/*.d ", qdir, "/*.raw ", qdir, "/*.mzML ", qdir, "/*.wiff",
    " 2>/dev/null; true"
  )
  result <- ssh_exec(ssh_config, cmd)

  lines <- result$stdout[nzchar(result$stdout)]
  if (length(lines) == 0) return(empty_df)

  parsed <- strsplit(lines, "\t")
  # Filter out malformed lines
  parsed <- parsed[vapply(parsed, length, integer(1)) >= 2]
  if (length(parsed) == 0) return(empty_df)

  data.frame(
    filename = vapply(parsed, function(x) basename(trimws(x[2])), character(1)),
    size_mb  = as.numeric(vapply(parsed, function(x) trimws(x[1]), character(1))),
    type     = vapply(parsed, function(x) {
      f <- x[2]
      if (grepl("\\.d/?$", f)) "Bruker .d"
      else if (grepl("\\.raw$", f, ignore.case = TRUE)) "Thermo .raw"
      else if (grepl("\\.wiff$", f, ignore.case = TRUE)) "SCIEX .wiff"
      else "mzML"
    }, character(1)),
    stringsAsFactors = FALSE
  )
}

#' Scan FASTA files on a remote host via SSH
#' @param ssh_config list(host, user, port, key_path)
#' @param fasta_dir Character — remote directory path
#' @return Named character vector (display name -> full path)
ssh_scan_fasta_files <- function(ssh_config, fasta_dir) {
  qdir <- shQuote(fasta_dir)
  cmd <- paste0("ls -1d ", qdir, "/*.fasta ", qdir, "/*.fa 2>/dev/null; true")
  result <- ssh_exec(ssh_config, cmd)

  paths <- result$stdout[nzchar(result$stdout)]
  # Filter out lines that are literal unexpanded globs (no matches)
  paths <- paths[!grepl("\\*", paths)]
  if (length(paths) == 0) return(character())

  names(paths) <- basename(paths)
  paths
}

# =============================================================================
# SLURM Helper Functions
# =============================================================================

#' Execute a command via the SLURM proxy (for running sbatch/squeue/sacct from
#' inside an Apptainer container). The proxy runs host-side and watches a shared
#' directory for command requests.
#' @param cmd Character — shell command to execute (e.g. "sbatch /path/to/script.sh")
#' @param timeout Numeric — seconds to wait for result (default 30)
#' @return list(status, stdout) where status is exit code and stdout is character vector
slurm_proxy_exec <- function(cmd, timeout = 30) {
  proxy_dir <- Sys.getenv("DELIMP_SLURM_PROXY", "")
  if (!nzchar(proxy_dir) || !dir.exists(proxy_dir)) {
    return(list(status = 1, stdout = "SLURM proxy not available"))
  }

  id <- paste0(as.integer(Sys.time()), "_", sample(1000:9999, 1))
  cmd_file <- file.path(proxy_dir, paste0("cmd_", id))
  result_file <- file.path(proxy_dir, paste0("result_", id))

  writeLines(cmd, cmd_file)

  # Poll for result
  start <- Sys.time()
  while (as.numeric(difftime(Sys.time(), start, units = "secs")) < timeout) {
    if (file.exists(result_file)) {
      lines <- readLines(result_file, warn = FALSE)
      unlink(result_file)
      rc <- as.integer(lines[1])
      stdout <- if (length(lines) > 1) lines[-1] else character(0)
      return(list(status = rc, stdout = stdout))
    }
    Sys.sleep(0.5)
  }

  # Timeout — clean up
  unlink(cmd_file)
  list(status = 1, stdout = "SLURM proxy timed out")
}

#' Check if the SLURM proxy is available (inside Apptainer container)
#' @return logical
slurm_proxy_available <- function() {
  proxy_dir <- Sys.getenv("DELIMP_SLURM_PROXY", "")
  nzchar(proxy_dir) && dir.exists(proxy_dir)
}

#' Check SLURM job status (local or remote via SSH)
#' @param job_id Character — SLURM job ID
#' @param ssh_config list(host, user, port, key_path) or NULL for local
#' @param sbatch_path Character — full path to sbatch (to derive squeue/sacct paths)
#' @return Character: "queued", "running", "completed", "failed", "cancelled", "unknown"
check_slurm_status <- function(job_id, ssh_config = NULL, sbatch_path = NULL) {
  # Derive squeue/sacct/scancel paths from sbatch path
  slurm_cmd <- function(cmd) {
    if (!is.null(sbatch_path)) {
      file.path(dirname(sbatch_path), cmd)
    } else {
      cmd
    }
  }

  # First try squeue (for active jobs)
  if (!is.null(ssh_config)) {
    squeue_result <- ssh_exec(ssh_config,
      sprintf("%s --job %s --format=%%T --noheader 2>/dev/null",
              slurm_cmd("squeue"), job_id),
      login_shell = is.null(sbatch_path))
    status_output <- if (squeue_result$status == 0) squeue_result$stdout else character(0)
  } else if (slurm_proxy_available()) {
    # Inside Apptainer container — use SLURM proxy
    squeue_cmd <- sprintf("%s --job %s --format=%%T --noheader 2>/dev/null",
                          slurm_cmd("squeue"), job_id)
    proxy_result <- slurm_proxy_exec(squeue_cmd, timeout = 15)
    status_output <- if (proxy_result$status == 0) proxy_result$stdout else character(0)
  } else {
    status_output <- tryCatch({
      system2(slurm_cmd("squeue"),
        args = c("--job", job_id, "--format=%T", "--noheader"),
        stdout = TRUE, stderr = TRUE)
    }, error = function(e) character(0))
  }

  if (length(status_output) > 0 && nzchar(trimws(status_output[1]))) {
    state <- toupper(trimws(status_output[1]))
    return(switch(state,
      "PENDING"   = "queued",
      "RUNNING"   = "running",
      "COMPLETING" = "running",
      tolower(state)
    ))
  }

  # Job not in queue — check sacct for final state
  if (!is.null(ssh_config)) {
    sacct_result <- ssh_exec(ssh_config,
      sprintf("%s -j %s --format=State --noheader --parsable2 2>/dev/null",
              slurm_cmd("sacct"), job_id),
      login_shell = is.null(sbatch_path))
    sacct_output <- if (sacct_result$status == 0) sacct_result$stdout else "UNKNOWN"
  } else if (slurm_proxy_available()) {
    sacct_cmd <- sprintf("%s -j %s --format=State --noheader --parsable2 2>/dev/null",
                         slurm_cmd("sacct"), job_id)
    proxy_result <- slurm_proxy_exec(sacct_cmd, timeout = 15)
    sacct_output <- if (proxy_result$status == 0) proxy_result$stdout else "UNKNOWN"
  } else {
    sacct_output <- tryCatch({
      system2(slurm_cmd("sacct"),
        args = c("-j", job_id, "--format=State", "--noheader", "--parsable2"),
        stdout = TRUE, stderr = TRUE)
    }, error = function(e) "UNKNOWN")
  }

  # Filter empty lines and look for meaningful state
  sacct_output <- trimws(sacct_output)
  sacct_output <- sacct_output[nzchar(sacct_output)]
  if (length(sacct_output) == 0) return("unknown")

  # Check all returned states (sacct may return multiple lines for job + steps)
  states <- toupper(sacct_output)
  if (any(grepl("COMPLETED", states))) return("completed")
  if (any(grepl("FAILED|TIMEOUT|OUT_OF_ME", states))) return("failed")
  if (any(grepl("CANCELLED", states))) return("cancelled")
  if (any(grepl("RUNNING", states))) return("running")
  if (any(grepl("PENDING", states))) return("queued")
  return("unknown")
}

#' Parse job ID from sbatch stdout
#' @param sbatch_stdout Character vector — stdout from system2("sbatch", ...)
#' @return Character: job ID, or NULL if parsing fails
parse_sbatch_output <- function(sbatch_stdout) {
  match_line <- grep("Submitted batch job", sbatch_stdout, value = TRUE)
  if (length(match_line) == 0) return(NULL)
  gsub(".*job\\s+", "", match_line[1])
}

#' Estimate search time for display
#' @return Character: human-readable estimate
estimate_search_time <- function(n_files, search_mode = "libfree", cpus = 64) {
  if (n_files == 0) return("")

  # Rough heuristics (minutes per file at 64 cores)
  min_per_file <- if (search_mode == "libfree") 45 else 20
  max_per_file <- if (search_mode == "libfree") 60 else 30

  # Scale by CPU count (assume ~linear scaling from 64)
  scale <- 64 / max(cpus, 4)
  min_per_file <- min_per_file * scale
  max_per_file <- max_per_file * scale

  total_min <- n_files * min_per_file
  total_max <- n_files * max_per_file

  format_time <- function(minutes) {
    if (minutes < 60) return(sprintf("%.0f min", minutes))
    hours <- minutes / 60
    if (hours < 24) return(sprintf("%.0f hours", ceiling(hours)))
    sprintf("%.1f days", hours / 24)
  }

  sprintf("~%s to %s for %d files", format_time(total_min), format_time(total_max), n_files)
}
