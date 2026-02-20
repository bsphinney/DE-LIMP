# helpers_search.R — Pure helper functions for DIA-NN HPC Search Integration
# No Shiny reactivity. All functions are testable standalone.

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
# sbatch Script Generation
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
  # Override defaults with provided params
  for (nm in names(search_params)) sp[[nm]] <- search_params[[nm]]

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

  # Build variable modification flags
  var_mod_flags <- c()
  if (isTRUE(sp$mod_met_ox)) {
    var_mod_flags <- c(var_mod_flags, "    --var-mod UniMod:35,15.994915,M")
  }
  if (isTRUE(sp$mod_nterm_acetyl)) {
    var_mod_flags <- c(var_mod_flags, "    --var-mod UniMod:1,42.010565,*n")
  }
  if (nzchar(sp$extra_var_mods)) {
    extra_lines <- trimws(strsplit(sp$extra_var_mods, "\n")[[1]])
    for (mod in extra_lines) {
      if (nzchar(mod)) var_mod_flags <- c(var_mod_flags, sprintf("    --var-mod %s", mod))
    }
  }
  var_mod_str <- if (length(var_mod_flags) > 0) {
    paste(var_mod_flags, collapse = " \\\n")
  } else ""

  # Phospho-specific flags
  is_phospho <- identical(search_mode, "phospho")

  # Build DIA-NN command
  if (search_mode == "library" && !is.null(speclib_path) && nzchar(speclib_path)) {
    # Library mode
    diann_cmd_parts <- c(
      sprintf("apptainer exec --bind %s %s /diann-2.3.0/diann-linux \\", bind_mount, diann_sif),
      paste0(run_flags, " \\"),
      paste0(fasta_flags, " \\"),
      sprintf("    --lib /work/lib/%s \\", basename(speclib_path)),
      sprintf("    --threads %d \\", cpus),
      "    --verbose 1 \\",
      sprintf("    --out /work/out/%s \\", report_name),
      sprintf("    --qvalue %s \\", sp$qvalue),
      "    --matrices \\",
      "    --out-lib /work/out/report-lib.parquet \\",
      "    --gen-spec-lib \\",
      if (isTRUE(sp$xic)) "    --xic \\" else NULL,
      if (isTRUE(sp$unimod4)) "    --unimod4 \\" else NULL,
      sprintf("    --var-mods %d \\", sp$max_var_mods),
      sprintf("    --window %d \\", sp$scan_window),
      if (sp$mass_acc_mode == "manual") c(
        sprintf("    --mass-acc %s \\", sp$mass_acc),
        sprintf("    --mass-acc-ms1 %s \\", sp$mass_acc_ms1)
      ) else NULL,
      if (isTRUE(sp$mbr)) "    --reanalyse \\" else NULL,
      if (isTRUE(sp$rt_profiling)) "    --rt-profiling \\" else NULL,
      if (normalization == "off") "    --no-norm \\" else NULL,
      "    --use-quant",
      if (nzchar(var_mod_str)) var_mod_str else NULL
    )
  } else {
    # Library-free mode
    diann_cmd_parts <- c(
      sprintf("apptainer exec --bind %s %s /diann-2.3.0/diann-linux \\", bind_mount, diann_sif),
      paste0(run_flags, " \\"),
      paste0(fasta_flags, " \\"),
      sprintf("    --out /work/out/%s \\", report_name),
      "    --out-lib /work/out/report-lib.parquet \\",
      "    --matrices \\",
      "    --fasta-search \\",
      sprintf("    --qvalue %s \\", sp$qvalue),
      "    --gen-spec-lib \\",
      "    --verbose 1 \\",
      "    --predictor \\",
      if (isTRUE(sp$xic)) "    --xic \\" else NULL,
      sprintf("    --threads %d \\", cpus),
      if (isTRUE(sp$met_excision)) "    --met-excision \\" else NULL,
      sprintf("    --min-pep-len %d \\", sp$min_pep_len),
      sprintf("    --max-pep-len %d \\", sp$max_pep_len),
      sprintf("    --min-pr-mz %d \\", sp$min_pr_mz),
      sprintf("    --max-pr-mz %d \\", sp$max_pr_mz),
      sprintf("    --min-pr-charge %d \\", sp$min_pr_charge),
      sprintf("    --max-pr-charge %d \\", sp$max_pr_charge),
      sprintf("    --min-fr-mz %d \\", sp$min_fr_mz),
      sprintf("    --max-fr-mz %d \\", sp$max_fr_mz),
      sprintf("    --cut %s \\", sp$enzyme),
      sprintf("    --missed-cleavages %d \\", sp$missed_cleavages),
      if (isTRUE(sp$unimod4)) "    --unimod4 \\" else NULL,
      sprintf("    --var-mods %d \\", sp$max_var_mods),
      if (sp$mass_acc_mode == "manual") c(
        sprintf("    --window %d \\", sp$scan_window),
        sprintf("    --mass-acc %s \\", sp$mass_acc),
        sprintf("    --mass-acc-ms1 %s \\", sp$mass_acc_ms1)
      ) else NULL,
      if (isTRUE(sp$mbr)) "    --reanalyse \\" else NULL,
      if (isTRUE(sp$rt_profiling)) "    --rt-profiling \\" else NULL,
      if (normalization == "off") "    --no-norm \\" else NULL,
      if (is_phospho) "    --phospho-output \\" else NULL,
      if (is_phospho) "    --report-lib-info \\" else NULL,
      if (nzchar(var_mod_str)) var_mod_str else NULL
    )
  }

  # Remove NULLs and trailing backslash on last line
  diann_cmd_parts <- Filter(Negate(is.null), diann_cmd_parts)
  last <- length(diann_cmd_parts)
  diann_cmd_parts[last] <- sub(" \\\\$", "", diann_cmd_parts[last])
  diann_cmd <- paste(diann_cmd_parts, collapse = "\n")

  # Extra CLI flags
  extra_cli <- ""
  if (nzchar(sp$extra_cli_flags)) {
    extra_cli <- paste0(" \\\n    ", trimws(sp$extra_cli_flags))
  }

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
    diann_cmd, extra_cli, '\n',
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
  } else {
    status_output <- tryCatch({
      system2("squeue",
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
  } else {
    sacct_output <- tryCatch({
      system2("sacct",
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
