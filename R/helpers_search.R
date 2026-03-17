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

  # One-per-gene: use FTP (REST API &onePerGene=true is silently ignored)
  if (content_type == "one_per_gene") {
    return(download_uniprot_fasta_ftp(proteome_id, output_path))
  }

  # All other content types: use REST API
  base_query <- sprintf("(proteome:%s)", proteome_id)
  query <- switch(content_type,
    "reviewed"          = paste0(base_query, " AND (reviewed:true)"),
    "reviewed_isoforms" = paste0(base_query, " AND (reviewed:true)"),
    "full"              = base_query,
    "full_isoforms"     = base_query,
    base_query
  )

  include_isoform <- content_type %in% c("full_isoforms", "reviewed_isoforms")
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "query=", utils::URLencode(query),
    "&format=fasta",
    "&compressed=false",
    if (include_isoform) "&includeIsoform=true" else ""
  )

  tryCatch({
    tmp_file <- tempfile(fileext = ".fasta")

    resp <- httr2::request(url) |>
      httr2::req_headers(Accept = "text/plain") |>
      httr2::req_timeout(300) |>
      httr2::req_perform(path = tmp_file)

    if (!file.exists(tmp_file) || file.size(tmp_file) < 100) {
      stop("Download failed or returned empty file")
    }

    n_seqs <- sum(grepl("^>", readLines(tmp_file, warn = FALSE)))

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

#' Download one-per-gene FASTA from UniProt FTP (reference proteomes)
#'
#' The REST API &onePerGene=true parameter is silently ignored, so we use the
#' FTP reference proteome files which are the true canonical one-per-gene sets.
#' URL pattern: https://ftp.uniprot.org/pub/databases/uniprot/current_release/
#'   knowledgebase/reference_proteomes/{Kingdom}/{UPID}/{UPID}_{TAXID}.fasta.gz
#'
#' Falls back to REST API (full proteome) if FTP download fails (e.g., organism
#' is not a reference proteome or FTP file doesn't exist).
#'
#' @param proteome_id Character — UniProt proteome ID (e.g., "UP000005640")
#' @param output_path Character — full path where FASTA will be saved
#' @return List with success status, path, sequence count, file size
download_uniprot_fasta_ftp <- function(proteome_id, output_path) {
  ftp_base <- "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes"

  # Map superkingdom to FTP directory name
  kingdom_map <- c(
    eukaryota = "Eukaryota",
    bacteria  = "Bacteria",
    archaea   = "Archaea",
    viruses   = "Viruses"
  )

  tryCatch({
    # Query proteomes API for superkingdom + taxonomy ID
    meta_url <- sprintf("https://rest.uniprot.org/proteomes/%s?format=json", proteome_id)
    meta_resp <- httr2::request(meta_url) |>
      httr2::req_headers(Accept = "application/json") |>
      httr2::req_timeout(30) |>
      httr2::req_perform()

    meta <- httr2::resp_body_json(meta_resp)
    superkingdom <- tolower(meta$superkingdom %||% "")
    taxon_id <- meta$taxonomy$taxonId %||% 0L

    if (!nzchar(superkingdom) || taxon_id == 0L) {
      stop("Could not determine superkingdom or taxonomy ID from proteomes API")
    }

    kingdom_dir <- kingdom_map[[superkingdom]]
    if (is.null(kingdom_dir)) {
      stop(sprintf("Unknown superkingdom '%s' — cannot map to FTP directory", superkingdom))
    }

    # Build FTP URL: {Kingdom}/{UPID}/{UPID}_{TAXID}.fasta.gz
    ftp_url <- sprintf("%s/%s/%s/%s_%s.fasta.gz",
                       ftp_base, kingdom_dir, proteome_id, proteome_id, taxon_id)

    # Download compressed FASTA
    tmp_gz <- tempfile(fileext = ".fasta.gz")
    tmp_fasta <- tempfile(fileext = ".fasta")

    resp <- httr2::request(ftp_url) |>
      httr2::req_timeout(300) |>
      httr2::req_perform(path = tmp_gz)

    if (!file.exists(tmp_gz) || file.size(tmp_gz) < 100) {
      stop("FTP download returned empty file")
    }

    # Decompress .gz → .fasta using text-mode gzfile connection
    fasta_lines <- readLines(gzfile(tmp_gz), warn = FALSE)

    if (length(fasta_lines) < 2) {
      stop("Decompression produced empty file")
    }

    writeLines(fasta_lines, tmp_fasta)
    n_seqs <- sum(grepl("^>", fasta_lines))

    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    file.copy(tmp_fasta, output_path, overwrite = TRUE)
    unlink(c(tmp_gz, tmp_fasta))

    list(
      success = TRUE,
      path = output_path,
      n_sequences = n_seqs,
      file_size = file.size(output_path),
      url = ftp_url,
      source = "ftp"
    )
  }, error = function(e) {
    message(sprintf("[DE-LIMP] FTP one-per-gene download failed: %s. Falling back to REST API.", e$message))
    # Fallback: REST API full proteome with warning
    fallback_result <- download_uniprot_fasta_rest_fallback(proteome_id, output_path)
    if (fallback_result$success) {
      fallback_result$warning <- paste0(
        "One-per-gene FASTA not available via FTP for this organism. ",
        "Downloaded full proteome instead (may contain multiple isoforms per gene).")
    }
    fallback_result
  })
}

#' REST API fallback for one-per-gene when FTP is unavailable
#' @keywords internal
download_uniprot_fasta_rest_fallback <- function(proteome_id, output_path) {
  url <- paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "query=", utils::URLencode(sprintf("(proteome:%s)", proteome_id)),
    "&format=fasta",
    "&compressed=false"
  )
  tryCatch({
    tmp_file <- tempfile(fileext = ".fasta")
    resp <- httr2::request(url) |>
      httr2::req_headers(Accept = "text/plain") |>
      httr2::req_timeout(300) |>
      httr2::req_perform(path = tmp_file)

    if (!file.exists(tmp_file) || file.size(tmp_file) < 100) {
      stop("Fallback download failed or returned empty file")
    }

    n_seqs <- sum(grepl("^>", readLines(tmp_file, warn = FALSE)))
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    file.copy(tmp_file, output_path, overwrite = TRUE)
    unlink(tmp_file)

    list(
      success = TRUE,
      path = output_path,
      n_sequences = n_seqs,
      file_size = file.size(output_path),
      url = url,
      source = "rest_fallback"
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
    "one_per_gene"      = "opg",
    "reviewed"          = "sprot",
    "reviewed_isoforms" = "sprot_iso",
    "full"              = "full",
    "full_isoforms"     = "full_iso",
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

# =============================================================================
# NCBI Datasets API — Search and Download Proteomes
# =============================================================================

#' Search NCBI for genome assemblies with protein annotations
#' @param query Character — organism name (e.g., "Peromyscus californicus")
#' @return data.frame with accession, organism, assembly_level, protein_count, annotation_name
ncbi_search_assemblies <- function(query) {
  url <- sprintf(
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon/%s/dataset_report",
    utils::URLencode(query, reserved = TRUE)
  )

  resp <- tryCatch(
    httr2::request(url) |>
      httr2::req_headers(Accept = "application/json") |>
      httr2::req_timeout(30) |>
      httr2::req_perform(),
    error = function(e) {
      message("[NCBI] Search failed: ", e$message)
      return(NULL)
    }
  )
  if (is.null(resp)) return(data.frame())

  body <- httr2::resp_body_json(resp)
  reports <- body$reports
  if (length(reports) == 0) return(data.frame())

  rows <- lapply(reports, function(r) {
    acc <- r$accession %||% ""
    paired <- r$paired_accession %||% ""
    org <- r$organism$organism_name %||% ""
    level <- r$assembly_info$assembly_level %||% ""
    refseq_cat <- r$assembly_info$refseq_category %||% ""
    ann <- r$annotation_info
    ann_name <- ann$name %||% ""
    prot_count <- tryCatch(
      ann$stats$gene_counts$protein_coding %||% 0L,
      error = function(e) 0L
    )
    # Prefer RefSeq (GCF_) accession which has the annotation
    best_acc <- if (nzchar(paired) && grepl("^GCF_", paired)) paired else acc
    data.frame(
      accession = best_acc,
      genbank = acc,
      organism = org,
      assembly_level = level,
      refseq_category = refseq_cat,
      protein_count = as.integer(prot_count),
      annotation = ann_name,
      stringsAsFactors = FALSE
    )
  })

  df <- do.call(rbind, rows)
  # Keep only annotated assemblies (protein_count > 0), deduplicate by accession
  df <- df[df$protein_count > 0, , drop = FALSE]
  df <- df[!duplicated(df$accession), , drop = FALSE]
  # Sort: reference genome first, then by protein count descending
  df <- df[order(df$refseq_category == "reference genome", -df$protein_count,
                 decreasing = c(TRUE, FALSE)), , drop = FALSE]
  rownames(df) <- NULL
  df
}

#' Download protein FASTA from NCBI for a genome assembly
#' @param accession Character — RefSeq accession (e.g., "GCF_007827085.1")
#' @param output_dir Character — directory to save the FASTA file
#' @return Path to the downloaded .fasta file, or NULL on failure
ncbi_download_proteome <- function(accession, output_dir) {
  url <- sprintf(
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/%s/download?include_annotation_type=PROT_FASTA",
    accession
  )

  zip_path <- file.path(output_dir, paste0(accession, "_protein.zip"))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # Download ZIP — use download.file for reliability across environments
  message("[NCBI] Downloading proteome from: ", url)
  message("[NCBI] Saving to: ", zip_path)
  dl_result <- tryCatch(
    download.file(url, zip_path, mode = "wb", quiet = TRUE, timeout = 300),
    error = function(e) {
      message("[NCBI] Download failed: ", e$message)
      return(1L)
    }
  )
  if (dl_result != 0 || !file.exists(zip_path) || file.size(zip_path) < 1000) {
    message("[NCBI] Download failed or file too small: ",
            if (file.exists(zip_path)) paste0(file.size(zip_path), " bytes") else "file missing")
    unlink(zip_path)
    return(NULL)
  }
  message("[NCBI] Downloaded: ", round(file.size(zip_path) / 1e6, 1), " MB")

  # Extract protein.faa from ZIP
  fasta_entry <- tryCatch({
    entries <- utils::unzip(zip_path, list = TRUE)
    faa <- entries$Name[grepl("protein\\.faa$", entries$Name)]
    if (length(faa) == 0) return(NULL)
    faa[1]
  }, error = function(e) NULL)

  if (is.null(fasta_entry)) {
    unlink(zip_path)
    return(NULL)
  }

  utils::unzip(zip_path, files = fasta_entry, exdir = output_dir, junkpaths = TRUE)
  unlink(zip_path)

  fasta_path <- file.path(output_dir, basename(fasta_entry))

  # Rename to a descriptive filename
  final_name <- paste0(accession, "_protein.fasta")
  final_path <- file.path(output_dir, final_name)
  file.rename(fasta_path, final_path)

  final_path
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
    min_pr_mz = 300, max_pr_mz = 1800,
    min_pr_charge = 1, max_pr_charge = 4,
    min_fr_mz = 200, max_fr_mz = 1800,
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
    sprintf("--var-mods %d", as.integer(sp$max_var_mods))
  )

  if (isTRUE(sp$xic)) flags <- c(flags, "--xic")
  if (isTRUE(sp$unimod4)) flags <- c(flags, "--unimod4")

  # Library mode
  if (search_mode == "library" && !is.null(speclib_mount)) {
    flags <- c(flags,
      sprintf("--lib %s", speclib_mount),
      sprintf("--window %d", as.integer(sp$scan_window)),
      "--use-quant"
    )
  }

  # Library-free mode (and phospho)
  if (search_mode != "library") {
    flags <- c(flags,
      "--fasta-search",
      "--predictor",
      sprintf("--cut %s", sp$enzyme),
      sprintf("--missed-cleavages %d", as.integer(sp$missed_cleavages)),
      sprintf("--min-pep-len %d", as.integer(sp$min_pep_len)),
      sprintf("--max-pep-len %d", as.integer(sp$max_pep_len)),
      sprintf("--min-pr-mz %d", as.integer(sp$min_pr_mz)),
      sprintf("--max-pr-mz %d", as.integer(sp$max_pr_mz)),
      sprintf("--min-pr-charge %d", as.integer(sp$min_pr_charge)),
      sprintf("--max-pr-charge %d", as.integer(sp$max_pr_charge)),
      sprintf("--min-fr-mz %d", as.integer(sp$min_fr_mz)),
      sprintf("--max-fr-mz %d", as.integer(sp$max_fr_mz))
    )
    if (isTRUE(sp$met_excision)) flags <- c(flags, "--met-excision")
  }

  # Mass accuracy
  if (sp$mass_acc_mode == "manual") {
    flags <- c(flags,
      sprintf("--window %d", as.integer(sp$scan_window)),
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
# DIA-NN Log File Parsing (inverse of build_diann_flags)
# =============================================================================

#' Parse a DIA-NN log file and extract search parameters
#' @param log_path Path to DIA-NN log file (.log, .txt, .out)
#' @return List with success, message, params, search_mode, normalization,
#'   version, fasta_files, n_raw_files, command_line
parse_diann_log <- function(log_path) {
  fail <- function(msg) list(success = FALSE, message = msg,
    params = list(), search_mode = NULL, normalization = NULL,
    version = NULL, fasta_files = character(), n_raw_files = 0L,
    command_line = NULL)

  if (!file.exists(log_path)) return(fail("File not found."))

  lines <- tryCatch(readLines(log_path, warn = FALSE), error = function(e) NULL)
  if (is.null(lines) || length(lines) == 0) return(fail("Could not read file or file is empty."))

  # Strip \r from Windows-style logs
  lines <- gsub("\r", "", lines)

  # Extract DIA-NN version from early log lines
  version <- NULL
  ver_match <- regmatches(lines, regexpr("DIA-NN\\s+([0-9]+\\.[0-9]+\\.?[0-9]*)", lines))
  if (length(ver_match) > 0 && any(nzchar(ver_match))) {
    ver_line <- ver_match[nzchar(ver_match)][1]
    version <- sub("DIA-NN\\s+", "", ver_line)
  }

  # Find command line: first line matching diann binary followed by flags
  cmd_idx <- grep("^(diann[^[:space:]]*|.*/diann[^[:space:]]*)\\s+--", lines)
  if (length(cmd_idx) == 0) return(fail("No DIA-NN command line found in log file."))
  cmd_line <- lines[cmd_idx[1]]

  # Split into flag chunks: first remove the binary name, then split on --
  binary_end <- regexpr("\\s+--", cmd_line)
  if (binary_end < 1) return(fail("Could not parse command line."))
  flags_str <- substring(cmd_line, binary_end + 1)

  # Split on whitespace-preceded -- to get individual flag chunks
  # Each chunk is like "flag-name value" or just "flag-name"
  chunks <- strsplit(flags_str, "\\s+--")[[1]]
  # First chunk still has leading --, strip it
  chunks[1] <- sub("^--", "", chunks[1])

  # Parse each chunk into flag name + value
  parsed <- list()
  for (chunk in chunks) {
    chunk <- trimws(chunk)
    if (!nzchar(chunk)) next
    # Split flag name from value (value may contain spaces, e.g. var-mod)
    parts <- regmatches(chunk, regexec("^([^[:space:]]+)(\\s+(.*))?$", chunk))[[1]]
    flag_name <- parts[2]
    flag_value <- if (length(parts) >= 4 && nzchar(parts[4])) trimws(parts[4]) else NULL
    parsed <- c(parsed, list(list(flag = flag_name, value = flag_value)))
  }

  # Operational flags to skip (these don't map to user-facing settings)
  # Note: "lib" and "out-lib" are NOT here — handled separately
  skip_flags <- c("f", "out", "temp", "threads", "verbose",
                  "matrices", "gen-spec-lib", "use-quant", "no-ifs-removal",
                  "report-lib-info", "quant-ori-names")


  # Valid enzyme values for the dropdown

  valid_enzymes <- c("K*,R*", "K,R", "K", "F,W,Y,L", "-")

  # Initialize output params with NULLs (only set what the log contains)
  params <- list()
  fasta_files <- character()
  n_raw_files <- 0L
  var_mods <- character()
  extra_cli_parts <- character()
  has_mass_acc <- FALSE
  has_fasta_search <- FALSE
  has_phospho <- FALSE
  has_lib <- FALSE
  lib_path <- NULL
  out_lib_path <- NULL
  has_no_norm <- FALSE

  # Value flag mapping: DIA-NN flag -> params key + type
  value_map <- list(
    "qvalue"           = list(key = "qvalue", type = "numeric"),
    "var-mods"         = list(key = "max_var_mods", type = "integer"),
    "window"           = list(key = "scan_window", type = "integer"),
    "mass-acc"         = list(key = "mass_acc", type = "numeric"),
    "mass-acc-ms1"     = list(key = "mass_acc_ms1", type = "numeric"),
    "cut"              = list(key = "enzyme", type = "character"),
    "missed-cleavages" = list(key = "missed_cleavages", type = "integer"),
    "min-pep-len"      = list(key = "min_pep_len", type = "integer"),
    "max-pep-len"      = list(key = "max_pep_len", type = "integer"),
    "min-pr-mz"        = list(key = "min_pr_mz", type = "integer"),
    "max-pr-mz"        = list(key = "max_pr_mz", type = "integer"),
    "min-pr-charge"    = list(key = "min_pr_charge", type = "integer"),
    "max-pr-charge"    = list(key = "max_pr_charge", type = "integer"),
    "min-fr-mz"        = list(key = "min_fr_mz", type = "integer"),
    "max-fr-mz"        = list(key = "max_fr_mz", type = "integer"),
    "pg-level"         = list(key = "pg_level", type = "integer")
  )

  # Boolean flag mapping: DIA-NN flag -> params key
  bool_map <- list(
    "reanalyse"    = "mbr",
    "rt-profiling" = "rt_profiling",
    "xic"          = "xic",
    "unimod4"      = "unimod4",
    "met-excision" = "met_excision",
    "proteoforms"  = "proteoforms"
  )

  # Track which boolean flags we see
  seen_bools <- character()

  for (item in parsed) {
    fl <- item$flag
    val <- item$value

    # Count raw files
    if (fl == "f") {
      n_raw_files <- n_raw_files + 1L
      next
    }

    # Collect FASTA paths
    if (fl == "fasta") {
      if (!is.null(val)) fasta_files <- c(fasta_files, val)
      next
    }

    # Detect library mode and store path
    if (fl == "lib") {
      has_lib <- TRUE
      lib_path <- val
      next
    }

    # Store output library path
    if (fl == "out-lib") {
      out_lib_path <- val
      next
    }

    # Skip operational flags
    if (fl %in% skip_flags) next

    # Special: --fasta-search
    if (fl == "fasta-search") {
      has_fasta_search <- TRUE
      next
    }

    # Special: --predictor (skip, implied by fasta-search)
    if (fl == "predictor") next

    # Special: --phospho-output
    if (fl == "phospho-output") {
      has_phospho <- TRUE
      next
    }

    # Special: --no-norm
    if (fl == "no-norm") {
      has_no_norm <- TRUE
      next
    }

    # Special: --var-mod
    if (fl == "var-mod") {
      if (!is.null(val)) {
        if (grepl("^UniMod:35", val)) {
          params$mod_met_ox <- TRUE
        } else if (grepl("^UniMod:1", val)) {
          params$mod_nterm_acetyl <- TRUE
        } else {
          var_mods <- c(var_mods, val)
        }
      }
      next
    }

    # Value flags
    if (fl %in% names(value_map)) {
      mapping <- value_map[[fl]]
      if (!is.null(val)) {
        converted <- switch(mapping$type,
          "numeric" = suppressWarnings(as.numeric(val)),
          "integer" = suppressWarnings(as.integer(val)),
          "character" = val
        )
        # Enzyme validation
        if (fl == "cut") {
          if (!(val %in% valid_enzymes)) {
            extra_cli_parts <- c(extra_cli_parts, paste0("--", fl, " ", val))
            next
          }
        }
        if (!is.na(converted)) params[[mapping$key]] <- converted
        if (fl %in% c("mass-acc", "mass-acc-ms1")) has_mass_acc <- TRUE
      }
      next
    }

    # Boolean flags
    if (fl %in% names(bool_map)) {
      params[[bool_map[[fl]]]] <- TRUE
      seen_bools <- c(seen_bools, fl)
      next
    }

    # Unrecognized flag → extra_cli_flags
    flag_str <- paste0("--", fl)
    if (!is.null(val)) flag_str <- paste(flag_str, val)
    extra_cli_parts <- c(extra_cli_parts, flag_str)
  }

  # Boolean flags not seen → FALSE
  for (fl in names(bool_map)) {
    key <- bool_map[[fl]]
    if (!(fl %in% seen_bools)) params[[key]] <- FALSE
  }

  # Set mod defaults if not seen
  if (is.null(params$mod_met_ox)) params$mod_met_ox <- FALSE
  if (is.null(params$mod_nterm_acetyl)) params$mod_nterm_acetyl <- FALSE

  # Extra var mods
  if (length(var_mods) > 0) {
    params$extra_var_mods <- paste(var_mods, collapse = "\n")
  }

  # Mass accuracy mode
  if (has_mass_acc) {
    params$mass_acc_mode <- "manual"
  } else {
    params$mass_acc_mode <- "auto"
  }

  # Extra CLI flags
  if (length(extra_cli_parts) > 0) {
    params$extra_cli_flags <- paste(extra_cli_parts, collapse = " ")
  }

  # Search mode
  # If --fasta + --cut are present alongside --lib, this is likely a second pass
  # of a library-free search (DIA-NN uses --lib for the predicted speclib).
  # Prefer libfree so the user can reproduce the full workflow.
  has_fasta_with_enzyme <- length(fasta_files) > 0 && !is.null(params$enzyme)
  search_mode <- if (has_phospho) {
    "phospho"
  } else if (has_fasta_search || has_fasta_with_enzyme) {
    "libfree"
  } else if (has_lib) {
    "library"
  } else {
    "libfree"
  }

  # Normalization
  normalization <- if (has_no_norm) "off" else "on"

  # Extract library precursor count from log body
  # Lines like "6028174 precursors generated" or "4612085 precursors in the library"
  n_precursors_library <- tryCatch({
    prec_lines <- grep("precursors (generated|in)", lines, value = TRUE)
    if (length(prec_lines) > 0) {
      # Take the last match (final library size after all steps)
      m <- regmatches(prec_lines[length(prec_lines)],
                      regexpr("[0-9,]+(?=\\s+precursors)", prec_lines[length(prec_lines)], perl = TRUE))
      if (length(m) > 0) as.integer(gsub(",", "", m)) else NULL
    } else NULL
  }, error = function(e) NULL)

  # Pipeline step detection from SLURM job name or log content
  pipeline_step <- tryCatch({
    job_lines <- grep("JobName=|Step [0-9]+/[0-9]+:", lines, value = TRUE)
    if (length(job_lines) > 0) trimws(job_lines[1]) else NULL
  }, error = function(e) NULL)

  list(
    success = TRUE,
    message = NULL,
    params = params,
    search_mode = search_mode,
    normalization = normalization,
    version = version,
    fasta_files = fasta_files,
    n_raw_files = n_raw_files,
    command_line = cmd_line,
    # Comparator-relevant fields
    n_precursors_library = n_precursors_library,
    lib_path = lib_path,
    out_lib_path = out_lib_path,
    pipeline_step = pipeline_step
  )
}

# =============================================================================
# sbatch Script Generation (HPC backend)
# =============================================================================

#' Generate a complete sbatch script for DIA-NN search
#' @return Character string: complete sbatch script content
generate_sbatch_script <- function(
  analysis_name, raw_files, fasta_files, speclib_path = NULL,
  output_dir, diann_sif, normalization = "on", search_mode = "libfree",
  cpus = 64, mem_gb = 128, time_hours = 12,
  partition = "high", account = "genome-center-grp",
  search_params = list(), requeue = FALSE
) {
  # Determine output filename
  report_name <- if (normalization == "off") "no_norm_report.parquet" else "report.parquet"

  # Determine unique directories for data and fasta
  data_dirs <- unique(dirname(raw_files))
  has_fasta <- length(fasta_files) > 0 && any(nzchar(fasta_files))

  # Build bind mount string — handle multiple data and FASTA directories
  if (length(data_dirs) == 1) {
    data_bind_parts <- sprintf("%s:/work/data", data_dirs[1])
    data_mount_map <- rep("/work/data", length(raw_files))
  } else {
    data_bind_parts <- sprintf("%s:/work/data%d", data_dirs, seq_along(data_dirs))
    data_mount_map <- sprintf("/work/data%d", match(dirname(raw_files), data_dirs))
  }
  bind_parts <- data_bind_parts
  if (has_fasta) {
    fasta_dirs <- unique(dirname(fasta_files))
    fasta_bind_parts <- if (length(fasta_dirs) == 1) {
      sprintf("%s:/work/fasta", fasta_dirs[1])
    } else {
      sprintf("%s:/work/fasta%d", fasta_dirs, seq_along(fasta_dirs))
    }
    bind_parts <- c(bind_parts, fasta_bind_parts)
  }
  bind_parts <- c(bind_parts, sprintf("%s:/work/out", output_dir))
  if (!is.null(speclib_path) && nzchar(speclib_path)) {
    bind_parts <- c(bind_parts, sprintf("%s:/work/lib", dirname(speclib_path)))
  }
  bind_mount <- paste(bind_parts, collapse = ",")

  # Build --f flags for raw files — map each file to its mount point
  run_flags <- paste(sprintf("    --f %s/%s", data_mount_map, basename(raw_files)),
                     collapse = " \\\n")

  # Build --fasta flags — map each file to its mount point (skip if library-only)
  fasta_flags <- NULL
  if (has_fasta) {
    fasta_mount_map <- if (length(fasta_dirs) == 1) {
      rep("/work/fasta", length(fasta_files))
    } else {
      sprintf("/work/fasta%d", match(dirname(fasta_files), fasta_dirs))
    }
    fasta_flags <- paste(sprintf("    --fasta %s/%s", fasta_mount_map, basename(fasta_files)),
                         collapse = " \\\n")
  }

  # Get shared DIA-NN flags via build_diann_flags()
  speclib_mount <- if (!is.null(speclib_path) && nzchar(speclib_path)) {
    sprintf("/work/lib/%s", basename(speclib_path))
  } else NULL
  shared_flags <- build_diann_flags(search_params, search_mode, normalization, speclib_mount)

  # Build DIA-NN command for apptainer
  diann_cmd_parts <- c(
    sprintf("apptainer exec --bind %s %s /diann-2.3.0/diann-linux \\", bind_mount, diann_sif),
    paste0(run_flags, " \\"),
    if (!is.null(fasta_flags)) paste0(fasta_flags, " \\"),
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
    sprintf('#SBATCH -o %s/logs/diann_%%j.out\n', output_dir),
    sprintf('#SBATCH -e %s/logs/diann_%%j.err\n', output_dir),
    sprintf('#SBATCH --account=%s\n', account),
    sprintf('#SBATCH --time=%d:00:00\n', time_hours),
    sprintf('#SBATCH --partition=%s\n', partition),
    if (isTRUE(requeue)) '#SBATCH --requeue\n' else '',
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
  # Include WorkDir and StdOut so we can find log files and output
  cmd <- paste0(
    sacct_bin,
    " --starttime=$(date -d '", days_back, " days ago' +%Y-%m-%d 2>/dev/null || ",
    "date -v-", days_back, "d +%Y-%m-%d)",
    " --format=JobID%20,JobName%50,State%20,Elapsed%15,WorkDir%120,StdOut%300",
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

  # StdOut is field 6 — contains log path template (may have %j, %A placeholders)
  df$std_out <- vapply(parsed, function(p) {
    if (length(p) >= 6) trimws(p[6]) else ""
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
# SSH ControlMaster multiplexing: reuse one TCP connection for all SSH/SCP calls.
# First call creates the socket; subsequent calls reuse it (no new TCP handshake).
# Eliminates HPC MaxStartups throttling from rapid SSH connections.
ssh_control_path <- function(ssh_config) {
  # Unix domain sockets have max 104-byte path on macOS — keep it short.
  file.path("/tmp", sprintf(".delimp_%s_%s",
    ssh_config$user, gsub("[^a-zA-Z0-9]", "", ssh_config$host)))
}

ssh_mux_args <- function(ssh_config) {
  c("-o", sprintf("ControlPath=%s", ssh_control_path(ssh_config)),
    "-o", "ControlMaster=auto",
    "-o", "ControlPersist=300")
}

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
    ssh_mux_args(ssh_config),
    paste0(ssh_config$user, "@", ssh_config$host),
    remote_cmd
  )
  # Use processx for timeout support if available, else system2
  stdout <- tryCatch({
    if (requireNamespace("processx", quietly = TRUE)) {
      res <- processx::run("ssh", args = args, timeout = timeout,
                           error_on_status = FALSE,
                           env = c("current", MallocStackLogging = ""))
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
    "-o", "ConnectTimeout=10",
    "-o", "BatchMode=yes",
    ssh_mux_args(ssh_config),
    paste0(ssh_config$user, "@", ssh_config$host, ":", remote_path),
    local_path
  )
  stdout <- tryCatch({
    if (requireNamespace("processx", quietly = TRUE)) {
      res <- processx::run("scp", args = args, timeout = 60,
                           error_on_status = FALSE,
                           env = c("current", MallocStackLogging = ""))
      out <- paste0(res$stdout, res$stderr)
      if (res$status != 0) attr(out, "status") <- res$status
      out
    } else {
      system2("scp", args = args, stdout = TRUE, stderr = TRUE)
    }
  }, error = function(e) {
    structure(conditionMessage(e), status = 1L)
  })
  status <- attr(stdout, "status") %||% 0L
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
    "-o", "ConnectTimeout=10",
    "-o", "BatchMode=yes",
    ssh_mux_args(ssh_config),
    local_path,
    paste0(ssh_config$user, "@", ssh_config$host, ":", remote_path)
  )
  stdout <- tryCatch({
    if (requireNamespace("processx", quietly = TRUE)) {
      res <- processx::run("scp", args = args, timeout = 60,
                           error_on_status = FALSE,
                           env = c("current", MallocStackLogging = ""))
      out <- paste0(res$stdout, res$stderr)
      if (res$status != 0) attr(out, "status") <- res$status
      out
    } else {
      system2("scp", args = args, stdout = TRUE, stderr = TRUE)
    }
  }, error = function(e) {
    structure(conditionMessage(e), status = 1L)
  })
  status <- attr(stdout, "status") %||% 0L
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
  # Include JobID so we can filter out .extern/.batch substeps
  sacct_fmt <- "JobID,State"
  if (!is.null(ssh_config)) {
    sacct_result <- ssh_exec(ssh_config,
      sprintf("%s -j %s --format=%s --noheader --parsable2 2>/dev/null",
              slurm_cmd("sacct"), job_id, sacct_fmt),
      login_shell = is.null(sbatch_path))
    sacct_output <- if (sacct_result$status == 0) sacct_result$stdout else "UNKNOWN"
  } else if (slurm_proxy_available()) {
    sacct_cmd <- sprintf("%s -j %s --format=%s --noheader --parsable2 2>/dev/null",
                         slurm_cmd("sacct"), job_id, sacct_fmt)
    proxy_result <- slurm_proxy_exec(sacct_cmd, timeout = 15)
    sacct_output <- if (proxy_result$status == 0) proxy_result$stdout else "UNKNOWN"
  } else {
    sacct_output <- tryCatch({
      system2(slurm_cmd("sacct"),
        args = c("-j", job_id, paste0("--format=", sacct_fmt), "--noheader", "--parsable2"),
        stdout = TRUE, stderr = TRUE)
    }, error = function(e) "UNKNOWN")
  }

  # Filter empty lines
  sacct_output <- trimws(sacct_output)
  sacct_output <- sacct_output[nzchar(sacct_output)]
  if (length(sacct_output) == 0) return("unknown")

  # Parse JobID|State lines — filter out .extern/.batch substeps
  # sacct returns lines like "12345|COMPLETED", "12345.extern|COMPLETED", "12345.batch|COMPLETED"
  # The .extern step always COMPLETED immediately (even for PENDING/FAILED jobs),
  # so we must only look at the main job line (no dot in the JobID).
  states <- character(0)
  for (line in sacct_output) {
    parts <- strsplit(line, "\\|")[[1]]
    if (length(parts) >= 2) {
      jid <- trimws(parts[1])
      st <- toupper(trimws(parts[2]))
      # Keep only main job line (no .extern, .batch, .0, etc.)
      if (!grepl("\\.", jid)) states <- c(states, st)
    } else {
      # Fallback for lines without JobID (shouldn't happen with new format)
      states <- c(states, toupper(trimws(line)))
    }
  }

  if (length(states) == 0) return("unknown")

  # Check states — order matters: FAILED before COMPLETED
  if (any(grepl("FAILED|TIMEOUT|OUT_OF_ME", states))) return("failed")
  if (any(grepl("CANCELLED", states))) return("cancelled")
  if (any(grepl("RUNNING", states))) return("running")
  if (any(grepl("PENDING", states))) return("queued")
  if (any(grepl("COMPLETED", states))) return("completed")
  return("unknown")
}

#' Get failed array task indices from a SLURM array job
#' @param array_job_id Character — SLURM array job ID (parent, e.g. "9591429")
#' @param ssh_config SSH config list or NULL
#' @param sbatch_path Full path to sbatch
#' @return list(failed_tasks = integer vector of 0-based task IDs,
#'              reasons = character vector, max_rss_gb = numeric) or NULL
get_failed_array_tasks <- function(array_job_id, ssh_config = NULL, sbatch_path = NULL) {
  slurm_cmd <- function(cmd) {
    if (!is.null(sbatch_path)) file.path(dirname(sbatch_path), cmd) else cmd
  }

  sacct_cmd <- sprintf(
    "%s -j %s --format=JobID,State,MaxRSS --noheader --parsable2 2>/dev/null",
    slurm_cmd("sacct"), array_job_id)

  result <- if (!is.null(ssh_config)) {
    ssh_exec(ssh_config, sacct_cmd, login_shell = is.null(sbatch_path), timeout = 15)
  } else {
    tryCatch({
      out <- system2(slurm_cmd("sacct"),
        args = c("-j", array_job_id, "--format=JobID,State,MaxRSS",
                 "--noheader", "--parsable2"),
        stdout = TRUE, stderr = TRUE)
      list(status = 0, stdout = out)
    }, error = function(e) list(status = 1, stdout = character()))
  }

  if (result$status != 0 || length(result$stdout) == 0) return(NULL)

  failed_tasks <- integer(0)
  reasons <- character(0)
  max_rss_bytes <- 0
  total_tasks <- 0L

  for (line in result$stdout) {
    parts <- strsplit(trimws(line), "\\|")[[1]]
    if (length(parts) < 2) next
    jid <- trimws(parts[1])
    st <- toupper(trimws(parts[2]))
    rss <- if (length(parts) >= 3) trimws(parts[3]) else ""

    # Only array task entries: contain "_" but no "." (excludes .batch/.extern)
    if (!grepl("_", jid) || grepl("\\.", jid)) next

    total_tasks <- total_tasks + 1L
    # Extract task index (e.g., "9591429_16" -> 16)
    task_idx <- as.integer(sub(".*_", "", jid))

    if (grepl("FAILED|TIMEOUT|OUT_OF_ME", st)) {
      failed_tasks <- c(failed_tasks, task_idx)
      reasons <- c(reasons, st)
    }

    # Parse MaxRSS (e.g., "67102388K" or "51254.50M")
    if (nzchar(rss)) {
      rss_val <- tryCatch({
        num <- as.numeric(gsub("[^0-9.]", "", rss))
        if (grepl("K", rss, ignore.case = TRUE)) num / (1024 * 1024)  # KB -> GB
        else if (grepl("M", rss, ignore.case = TRUE)) num / 1024      # MB -> GB
        else if (grepl("G", rss, ignore.case = TRUE)) num             # GB
        else num / (1024^3)                                            # bytes -> GB
      }, error = function(e) 0)
      max_rss_bytes <- max(max_rss_bytes, rss_val)
    }
  }

  if (length(failed_tasks) == 0) return(NULL)

  list(
    failed_tasks = sort(failed_tasks),
    reasons = reasons,
    max_rss_gb = round(max_rss_bytes, 1),
    n_failed = length(failed_tasks),
    n_total = total_tasks
  )
}

#' Get SLURM estimated start time for a queued job
#' @param job_id SLURM job ID
#' @param ssh_config SSH config list (NULL for local)
#' @param sbatch_path Full path to sbatch (used to derive squeue path)
#' @return List with est_start (character or NULL) and priority (integer or NULL)
get_slurm_start_time <- function(job_id, ssh_config = NULL, sbatch_path = NULL) {
  squeue_cmd <- if (!is.null(sbatch_path)) {
    file.path(dirname(sbatch_path), "squeue")
  } else "squeue"

  # Query est start, priority, and reason in one call
  cmd <- sprintf('%s --job %s --format="%%S|%%Q|%%r" --noheader 2>/dev/null', squeue_cmd, job_id)

  output <- if (!is.null(ssh_config)) {
    res <- ssh_exec(ssh_config, cmd, login_shell = is.null(sbatch_path), timeout = 10)
    if (res$status == 0) res$stdout else character(0)
  } else {
    tryCatch(system2(squeue_cmd,
      args = c("--job", job_id, '--format="%S|%Q|%r"', "--noheader"),
      stdout = TRUE, stderr = TRUE), error = function(e) character(0))
  }

  output <- trimws(output)
  output <- output[nzchar(output)]
  if (length(output) == 0) return(list(est_start = NULL, priority = NULL, reason = NULL))

  # Parse first line: "est_start|priority|reason"
  parts <- strsplit(output[1], "\\|")[[1]]
  est_start <- if (length(parts) >= 1 && !parts[1] %in% c("N/A", "")) parts[1] else NULL
  priority <- if (length(parts) >= 2) suppressWarnings(as.integer(parts[2])) else NULL
  reason <- if (length(parts) >= 3 && nzchar(trimws(parts[3]))) trimws(parts[3]) else NULL

  list(est_start = est_start, priority = priority, reason = reason)
}

#' Parse job ID from sbatch stdout
#' @param sbatch_stdout Character vector — stdout from system2("sbatch", ...)
#' @return Character: job ID, or NULL if parsing fails
parse_sbatch_output <- function(sbatch_stdout) {
  match_line <- grep("Submitted batch job", sbatch_stdout, value = TRUE)
  if (length(match_line) == 0) return(NULL)
  trimws(gsub(".*job[[:space:]]+", "", match_line[1]))
}

#' Estimate search time for display
#' @return Character: human-readable estimate
estimate_search_time <- function(n_files, search_mode = "libfree", cpus = 64,
                                  parallel = FALSE, jobs = NULL) {
  if (n_files == 0) return("")

  format_time <- function(minutes) {
    if (minutes < 60) return(sprintf("%.0f min", minutes))
    hours <- minutes / 60
    if (hours < 24) return(sprintf("%.1f hours", hours))
    sprintf("%.1f days", hours / 24)
  }

  # Try to compute from historical job data
  hist_rates <- NULL
  if (!is.null(jobs) && length(jobs) > 0) {
    rates <- c()
    for (j in jobs) {
      if (!identical(j$status, "completed")) next
      if (is.null(j$submitted_at) || is.null(j$completed_at)) next
      n <- j$n_files %||% 0
      if (n < 2) next
      elapsed <- as.numeric(difftime(j$completed_at, j$submitted_at, units = "mins"))
      if (elapsed < 1) next
      is_par <- isTRUE(j$parallel) || grepl("parallel", j$name, ignore.case = TRUE)
      if (is_par == parallel) rates <- c(rates, elapsed / n)
    }
    if (length(rates) >= 1) hist_rates <- rates
  }

  if (!is.null(hist_rates)) {
    med_rate <- median(hist_rates)
    lo <- min(hist_rates) * 0.9
    hi <- max(hist_rates) * 1.1
    # Ensure lo < hi
    if (lo >= hi) { lo <- med_rate * 0.8; hi <- med_rate * 1.2 }
    total_lo <- n_files * lo
    total_hi <- n_files * hi
    source_note <- sprintf(" (based on %d past %s)",
      length(hist_rates), if (length(hist_rates) == 1) "search" else "searches")
  } else {
    # Fallback heuristics (minutes per file)
    if (parallel) {
      min_per_file <- if (search_mode == "libfree") 3 else 2
      max_per_file <- if (search_mode == "libfree") 8 else 5
    } else {
      min_per_file <- if (search_mode == "libfree") 15 else 8
      max_per_file <- if (search_mode == "libfree") 30 else 15
      # Scale by CPU count (assume ~linear scaling from 64)
      scale <- 64 / max(cpus, 4)
      min_per_file <- min_per_file * scale
      max_per_file <- max_per_file * scale
    }
    total_lo <- n_files * min_per_file
    total_hi <- n_files * max_per_file
    source_note <- ""
  }

  sprintf("~%s to %s for %d files%s",
    format_time(total_lo), format_time(total_hi), n_files, source_note)
}

# =============================================================================
# Search Info Archive
# =============================================================================

#' Generate search_info.md content for archiving search metadata
#'
#' Creates a markdown file with all search parameters, job IDs, and file paths
#' so that search history is preserved even if SLURM purges its records.
#'
#' @param analysis_name Character: name of the analysis
#' @param output_dir Character: output directory path
#' @param raw_files Character vector: raw file paths
#' @param fasta_files Character vector: FASTA file paths
#' @param search_params List: search parameters
#' @param search_mode Character: "libfree", "lib", or "phospho"
#' @param normalization Character: "on" or "off"
#' @param sif_path Character: path to DIA-NN Apptainer SIF
#' @param job_ids Named list or character: job ID(s). For parallel: list(step1="id", ...).
#'   For single: character scalar.
#' @param parallel Logical: whether this is a parallel (5-step) search
#' @param resources Named list: CPU/memory/time per step
#' @param partition Character: SLURM partition
#' @param account Character: SLURM account
#' @return Character: markdown content
generate_search_info <- function(analysis_name, output_dir, raw_files, fasta_files,
                                  search_params, search_mode = "libfree",
                                  normalization = "on", sif_path = "",
                                  job_ids = NULL, parallel = FALSE,
                                  resources = list(), partition = "", account = "",
                                  cached_speclib = NULL, custom_fasta_sequences = NULL,
                                  instrument_metadata = NULL) {

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")

  # Job IDs section
  if (parallel && is.list(job_ids)) {
    step_labels <- c(step1 = "Library Prediction", step2 = "First-pass Quant (array)",
                     step3 = "Empirical Library Assembly", step4 = "Final-pass Quant (array)",
                     step5 = "Cross-run Report")
    job_lines <- vapply(names(job_ids), function(s) {
      sprintf("- **%s** (%s): `%s`", step_labels[s] %||% s, s, job_ids[[s]])
    }, character(1))
    job_section <- paste(c("### Job IDs (5-Step Parallel)", job_lines), collapse = "\n")
  } else {
    job_section <- sprintf("### Job ID\n- `%s`", as.character(job_ids)[1])
  }

  # Resources section
  res_lines <- character()
  if (length(resources) > 0) {
    for (nm in names(resources)) {
      r <- resources[[nm]]
      res_lines <- c(res_lines, sprintf("- **%s**: %d CPUs, %dG memory, %dh walltime",
        nm, r$cpus %||% 0, r$mem %||% 0, r$time %||% 0))
    }
  }
  res_section <- if (length(res_lines) > 0) {
    paste(c("### Resources", res_lines), collapse = "\n")
  } else ""

  # Search params section
  param_lines <- vapply(names(search_params), function(nm) {
    val <- search_params[[nm]]
    if (is.null(val)) return("")
    sprintf("- **%s**: `%s`", nm, paste(as.character(val), collapse = ", "))
  }, character(1))
  param_lines <- param_lines[nzchar(param_lines)]

  # Raw files
  file_lines <- sprintf("- `%s`", raw_files)
  if (length(file_lines) > 50) {
    file_lines <- c(file_lines[1:50], sprintf("- ... and %d more files", length(file_lines) - 50))
  }

  # FASTA files
  fasta_lines <- if (length(fasta_files) > 0 && any(nzchar(fasta_files))) {
    sprintf("- `%s`", fasta_files[nzchar(fasta_files)])
  } else "- (none)"

  # Cached spectral library section
  cache_section <- if (!is.null(cached_speclib)) {
    paste(c(
      "### Cached Spectral Library",
      sprintf("- **Source**: `%s` (from analysis: %s)",
              cached_speclib$speclib_path, cached_speclib$analysis_name),
      sprintf("- **Cache key**: `%s`", cached_speclib$key)
    ), collapse = "\n")
  } else ""

  # Custom FASTA sequences section
  custom_fasta_section <- if (!is.null(custom_fasta_sequences) &&
                               nzchar(trimws(custom_fasta_sequences))) {
    paste(c(
      "### Custom Protein Sequences",
      "```fasta",
      trimws(custom_fasta_sequences),
      "```"
    ), collapse = "\n")
  } else ""

  # Instrument metadata section
  inst_section <- ""
  if (!is.null(instrument_metadata)) {
    inst_lines <- c("### Instrument & Acquisition")
    if (!is.null(instrument_metadata$instrument_model))
      inst_lines <- c(inst_lines, sprintf("- **Instrument**: %s", instrument_metadata$instrument_model))
    if (!is.null(instrument_metadata$instrument_serial))
      inst_lines <- c(inst_lines, sprintf("- **Serial**: %s", instrument_metadata$instrument_serial))
    if (!is.null(instrument_metadata$acquisition_mode) && instrument_metadata$acquisition_mode != "unknown")
      inst_lines <- c(inst_lines, sprintf("- **Acquisition mode**: %s", instrument_metadata$acquisition_mode))
    if (!is.null(instrument_metadata$lc_system))
      inst_lines <- c(inst_lines, sprintf("- **LC system**: %s", instrument_metadata$lc_system))
    if (!is.null(instrument_metadata$lc_method))
      inst_lines <- c(inst_lines, sprintf("- **LC method**: %s", instrument_metadata$lc_method))
    if (!is.null(instrument_metadata$dia_windows))
      inst_lines <- c(inst_lines, sprintf("- **DIA windows**: %d", instrument_metadata$dia_windows))
    if (!is.na(instrument_metadata$mz_range_low %||% NA) && !is.na(instrument_metadata$mz_range_high %||% NA))
      inst_lines <- c(inst_lines, sprintf("- **m/z range**: %.0f-%.0f",
                                           instrument_metadata$mz_range_low, instrument_metadata$mz_range_high))
    if (length(inst_lines) > 1) inst_section <- paste(inst_lines, collapse = "\n")
  }

  paste(c(
    sprintf("# DIA-NN Search: %s", analysis_name),
    "",
    sprintf("**Submitted**: %s", timestamp),
    sprintf("**Output directory**: `%s`", output_dir),
    sprintf("**Log files**: `%s/logs/`", output_dir),
    sprintf("**Search mode**: %s", search_mode),
    sprintf("**Normalization**: %s", normalization),
    sprintf("**DIA-NN container**: `%s`", sif_path),
    sprintf("**Partition**: %s | **Account**: %s", partition, account),
    "",
    if (nzchar(inst_section)) c(inst_section, ""),
    job_section,
    "",
    res_section,
    "",
    "### Search Parameters",
    param_lines,
    "",
    sprintf("### Raw Files (%d)", length(raw_files)),
    file_lines,
    "",
    "### FASTA Files",
    fasta_lines,
    if (nzchar(cache_section)) c("", cache_section),
    if (nzchar(custom_fasta_section)) c("", custom_fasta_section)
  ), collapse = "\n")
}

# =============================================================================
# Parallel DIA-NN Search (5-Step Workflow)
# =============================================================================

#' Write a file_list.txt for SLURM array jobs
#' @param raw_files Character vector of raw file paths
#' @param output_dir Character: directory to write file_list.txt
#' @return Character: path to the written file
write_file_list <- function(raw_files, output_dir) {
  file_list_path <- file.path(output_dir, "file_list.txt")
  writeLines(raw_files, file_list_path)
  file_list_path
}

#' Generate 5 sbatch scripts for parallel DIA-NN search
#'
#' Implements the canonical parallel workflow:
#' Step 1: Library prediction (single job, no raw files)
#' Step 2: First-pass per-file (SLURM array, predicted library)
#' Step 3: Empirical library assembly (single job, --use-quant)
#' Step 4: Final per-file pass (SLURM array, empirical library)
#' Step 5: Cross-run report (single job, --use-quant, matrices)
#'
#' @param analysis_name Character: sanitized analysis name
#' @param raw_files Character vector: full paths to raw files
#' @param fasta_files Character vector: full paths to FASTA files
#' @param speclib_path Character or NULL: path to user-provided spectral library
#' @param output_dir Character: remote output directory
#' @param diann_sif Character: path to Apptainer SIF
#' @param normalization Character: "on" or "off"
#' @param search_mode Character: "libfree" or "phospho"
#' @param cpus_per_file Integer: CPUs per array task
#' @param mem_per_file Integer: GB per array task
#' @param time_per_file Integer: hours per array task
#' @param assembly_cpus Integer: CPUs for assembly/report steps
#' @param assembly_mem Integer: GB for assembly/report steps
#' @param assembly_time Integer: hours for assembly/report steps
#' @param partition Character: SLURM partition
#' @param account Character: SLURM account
#' @param search_params List: search parameters for build_diann_flags()
#' @param max_simultaneous Integer: max concurrent array tasks
#' @return Named list of 5 script strings
generate_parallel_scripts <- function(
  analysis_name, raw_files, fasta_files, speclib_path = NULL,
  output_dir, diann_sif, normalization = "on", search_mode = "libfree",
  cpus_per_file = 16, mem_per_file = 64, time_per_file = 2,
  libpred_cpus = 16, libpred_mem = 64, libpred_time = 4,
  assembly_cpus = 64, assembly_mem = 128, assembly_time = 12,
  partition = "high", account = "genome-center-grp",
  search_params = list(), max_simultaneous = 20,
  array_partition = NULL, array_account = NULL,
  assembly_partition = NULL, assembly_account = NULL
) {
  n_files <- length(raw_files)
  has_fasta <- length(fasta_files) > 0 && any(nzchar(fasta_files))
  has_speclib <- !is.null(speclib_path) && nzchar(speclib_path)
  report_name <- if (normalization == "off") "no_norm_report.parquet" else "report.parquet"

  # --- Shared bind mount computation ---
  data_dirs <- unique(dirname(raw_files))
  # Handle multiple data directories (e.g., files from different instruments/runs)
  if (length(data_dirs) == 1) {
    data_bind_parts <- sprintf("%s:/work/data", data_dirs[1])
    data_mount_map <- rep("/work/data", length(raw_files))
  } else {
    data_bind_parts <- sprintf("%s:/work/data%d", data_dirs, seq_along(data_dirs))
    data_mount_map <- sprintf("/work/data%d", match(dirname(raw_files), data_dirs))
  }

  fasta_bind_parts <- character()
  fasta_mount_map <- character()
  if (has_fasta) {
    fasta_dirs <- unique(dirname(fasta_files))
    if (length(fasta_dirs) == 1) {
      fasta_bind_parts <- sprintf("%s:/work/fasta", fasta_dirs[1])
      fasta_mount_map <- rep("/work/fasta", length(fasta_files))
    } else {
      fasta_bind_parts <- sprintf("%s:/work/fasta%d", fasta_dirs, seq_along(fasta_dirs))
      fasta_mount_map <- sprintf("/work/fasta%d", match(dirname(fasta_files), fasta_dirs))
    }
  }

  speclib_bind <- if (has_speclib) sprintf("%s:/work/lib", dirname(speclib_path)) else NULL
  out_bind <- sprintf("%s:/work/out", output_dir)

  # Full bind mount for assembly/report steps (need all data dirs)
  all_binds <- c(data_bind_parts, fasta_bind_parts, out_bind, speclib_bind)
  all_binds <- Filter(Negate(is.null), all_binds)
  full_bind_mount <- paste(all_binds, collapse = ",")

  # Per-file bind mount (single file dir + output)
  # Array jobs bind $FILE_DIR dynamically
  perfile_bind_parts <- c("${FILE_DIR}:/work/data", fasta_bind_parts, out_bind, speclib_bind)
  perfile_bind_parts <- Filter(Negate(is.null), perfile_bind_parts)
  perfile_bind_mount <- paste(perfile_bind_parts, collapse = ",")

  # FASTA flags
  fasta_flags_str <- ""
  if (has_fasta) {
    fasta_flags_str <- paste(sprintf("    --fasta %s/%s", fasta_mount_map, basename(fasta_files)),
                              collapse = " \\\n")
  }

  # All --f flags (for assembly/report steps) — map each file to its mount point
  all_f_flags <- paste(sprintf("    --f %s/%s", data_mount_map, basename(raw_files)),
                        collapse = " \\\n")

  # --- Base search params for build_diann_flags ---
  # Override settings for parallel mode
  parallel_sp <- search_params
  parallel_sp$mbr <- FALSE        # No MBR in parallel (5-step replaces it)
  parallel_sp$rt_profiling <- FALSE  # Controlled per-step
  # Force manual mass accuracy — DIA-NN warns about auto-optimisation when
  # reusing .quant files (Steps 3/5 use --use-quant). Auto-calibration in the
  # assembly step may differ from per-file calibration in Step 2, producing
  # inconsistent results. Per recommendation from Vadim Demichev (DIA-NN dev).
  parallel_sp$mass_acc_mode <- "manual"

  speclib_mount <- if (has_speclib) sprintf("/work/lib/%s", basename(speclib_path)) else NULL
  base_flags <- build_diann_flags(parallel_sp, search_mode, "on", speclib_mount)

  # Remove flags that are step-specific (we add them per-step)
  # --fasta-search and --predictor belong in Step 1 only; including them
  # in Steps 2-5 causes DIA-NN to re-digest the FASTA instead of using
  # the predicted/empirical library
  remove_patterns <- c("^--out-lib ", "^--matrices$", "^--gen-spec-lib$",
                        "^--reanalyse$", "^--rt-profiling$", "^--no-norm$",
                        "^--xic$", "^--lib ", "^--fasta-search$", "^--predictor$")
  step_flags <- base_flags
  for (pat in remove_patterns) {
    step_flags <- step_flags[!grepl(pat, step_flags)]
  }

  # Resolve per-step partition/account overrides (NULL = use default)
  arr_part <- array_partition %||% partition
  arr_acct <- array_account %||% account
  asm_part <- assembly_partition %||% partition
  asm_acct <- assembly_account %||% account
  arr_requeue <- tolower(arr_part) == "low"  # preemptible partition
  asm_requeue <- tolower(asm_part) == "low"

  # --- SBATCH header helper ---
  sbatch_header <- function(job_suffix, cpus, mem_gb, time_hours,
                            array_spec = NULL, step_partition = partition,
                            step_account = account, requeue = FALSE) {
    lines <- c(
      '#!/bin/bash -l',
      sprintf('#SBATCH --job-name=diann_%s_%s', analysis_name, job_suffix),
      sprintf('#SBATCH --cpus-per-task=%d', cpus),
      sprintf('#SBATCH --mem=%dG', mem_gb),
      sprintf('#SBATCH -o %s/logs/diann_%s_%%j.out', output_dir, job_suffix),
      sprintf('#SBATCH -e %s/logs/diann_%s_%%j.err', output_dir, job_suffix),
      sprintf('#SBATCH --account=%s', step_account),
      sprintf('#SBATCH --time=%d:00:00', time_hours),
      sprintf('#SBATCH --partition=%s', step_partition)
    )
    if (isTRUE(requeue)) {
      lines <- c(lines, '#SBATCH --requeue')
    }
    if (!is.null(array_spec)) {
      lines <- c(lines, sprintf('#SBATCH --array=%s', array_spec))
    }
    paste(lines, collapse = "\n")
  }

  # --- Apptainer exec prefix ---
  apptainer_cmd <- function(bind_mount) {
    sprintf("apptainer exec --bind %s %s /diann-2.3.0/diann-linux", bind_mount, diann_sif)
  }

  # --- Quant file verification block (bash) ---
  # Generates bash code to verify all expected .quant files exist before
  # running the assembly/report step. Prevents silent failures when array
  # jobs crash without producing quant files. (Per Vadim Demichev recommendation)
  # With --quant-ori-names on all steps, quant files are named BASENAME.quant
  # (e.g., sample.raw → sample.quant, sample.d → sample.quant)
  quant_verify_block <- function(quant_subdir, prev_step, max_missing_pct = 5) {
    paste0(
      sprintf('# Verify quant files from Step %d — log missing files without mutating file_list.txt\n', prev_step),
      sprintf('echo "Verifying Step %d quant files..."\n', prev_step),
      'MISSING=0\n',
      'TOTAL=0\n',
      sprintf('EXCLUDED_FILE="%s/.excluded_step%d.txt"\n', output_dir, prev_step),
      '> "$EXCLUDED_FILE"\n',
      sprintf('while IFS= read -r RAW_FILE; do\n'),
      '  TOTAL=$((TOTAL + 1))\n',
      '  BASENAME=$(basename "$RAW_FILE")\n',
      '  # --quant-ori-names: quant files use original basename with .quant extension\n',
      '  QUANT_NAME="${BASENAME%.*}.quant"\n',
      sprintf('  if [ ! -f "%s/%s/${QUANT_NAME}" ]; then\n', output_dir, quant_subdir),
      '    echo "MISSING: ${QUANT_NAME} (from ${RAW_FILE})"\n',
      '    echo "$RAW_FILE" >> "$EXCLUDED_FILE"\n',
      '    MISSING=$((MISSING + 1))\n',
      '  fi\n',
      sprintf('done < "%s/file_list.txt"\n', output_dir),
      'if [ $MISSING -gt 0 ]; then\n',
      sprintf('  MAX_MISSING=$(( TOTAL * %d / 100 ))\n', max_missing_pct),
      '  if [ $MAX_MISSING -lt 3 ]; then MAX_MISSING=3; fi\n',
      '  if [ $MISSING -le $MAX_MISSING ]; then\n',
      sprintf('    echo "WARNING: ${MISSING} of ${TOTAL} quant files missing from Step %d (within %d%% tolerance)."\n',
              prev_step, max_missing_pct),
      '    echo "Excluded files logged to ${EXCLUDED_FILE}. Continuing with available quant files."\n',
      '  else\n',
      sprintf('    echo "ERROR: ${MISSING} of ${TOTAL} quant files missing from Step %d — exceeds %d%% tolerance. Aborting."\n',
              prev_step, max_missing_pct),
      '    exit 1\n',
      '  fi\n',
      'else\n',
      '  rm -f "$EXCLUDED_FILE"\n',
      '  echo "All ${TOTAL} quant files verified."\n',
      'fi\n\n'
    )
  }

  # =========================================================================
  # Step 1 — Library Prediction (single job, no raw files)
  # =========================================================================
  step1_script <- if (has_speclib) {
    # Skip step 1 if user provides a spectral library
    NULL
  } else {
    step1_flags <- c(
      step_flags,
      "--fasta-search",
      "--predictor",
      "--gen-spec-lib",
      sprintf("--out-lib /work/out/step1.speclib"),
      sprintf("--threads %d", libpred_cpus)
    )
    # Remove --fasta-search and --predictor from step_flags if already present (avoid dups)
    step1_flags <- unique(step1_flags)

    step1_cmd_parts <- c(
      sprintf("%s \\", apptainer_cmd(full_bind_mount)),
      if (nzchar(fasta_flags_str)) paste0(fasta_flags_str, " \\"),
      sprintf("    --out /work/out/step1_lib.parquet \\"),
      paste0("    ", step1_flags)
    )
    step1_cmd_parts <- Filter(Negate(is.null), step1_cmd_parts)
    # Add line continuations
    for (i in seq_along(step1_cmd_parts)) {
      if (i < length(step1_cmd_parts) && !grepl(" \\\\$", step1_cmd_parts[i])) {
        step1_cmd_parts[i] <- paste0(step1_cmd_parts[i], " \\")
      }
    }
    last <- length(step1_cmd_parts)
    step1_cmd_parts[last] <- sub(" \\\\$", "", step1_cmd_parts[last])

    paste0(
      sbatch_header("s1_libpred", libpred_cpus, libpred_mem, libpred_time,
                    step_partition = arr_part, step_account = arr_acct,
                    requeue = arr_requeue), "\n\n",
      "module load apptainer\n\n",
      sprintf('echo "Step 1/5: Library Prediction for %s"\n', analysis_name),
      'echo "Started: $(date)"\n\n',
      paste(step1_cmd_parts, collapse = "\n"), "\n\n",
      'EXIT_CODE=$?\n',
      'echo "Step 1 finished with exit code: $EXIT_CODE"\n',
      'echo "Completed: $(date)"\n',
      'exit $EXIT_CODE\n'
    )
  }

  # =========================================================================
  # Step 2 — First-pass per-file (SLURM array)
  # =========================================================================
  predicted_lib <- if (has_speclib) speclib_mount else "/work/out/step1.predicted.speclib"
  array_spec_2 <- sprintf("0-%d%%%d", n_files - 1, max_simultaneous)

  step2_script <- paste0(
    sbatch_header("s2_firstpass", cpus_per_file, mem_per_file, time_per_file,
                  array_spec = array_spec_2, step_partition = arr_part,
                  step_account = arr_acct, requeue = arr_requeue), "\n\n",
    "module load apptainer\n\n",
    sprintf('echo "Step 2/5: First-pass file ${SLURM_ARRAY_TASK_ID} of %d"\n', n_files),
    'echo "Started: $(date)"\n\n',
    '# Read file path from file list\n',
    sprintf('FILE_LIST="%s/file_list.txt"\n', output_dir),
    'RAW_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$FILE_LIST")\n',
    'FILE_DIR=$(dirname "$RAW_FILE")\n',
    'FILE_BASE=$(basename "$RAW_FILE")\n\n',
    'if [ -z "$RAW_FILE" ]; then\n',
    '  echo "ERROR: No file found for array task $SLURM_ARRAY_TASK_ID"\n',
    '  exit 1\n',
    'fi\n\n',
    sprintf('echo "Processing: $RAW_FILE"\n\n'),
    sprintf('%s \\\n', apptainer_cmd(perfile_bind_mount)),
    sprintf('    --f /work/data/$FILE_BASE \\\n'),
    sprintf('    --lib %s \\\n', predicted_lib),
    sprintf('    --temp /work/out/quant_step2 \\\n'),
    '    --rt-profiling \\\n',
    '    --gen-spec-lib \\\n',
    '    --quant-ori-names \\\n',
    sprintf('    --threads %d \\\n', cpus_per_file),
    if (nzchar(fasta_flags_str)) paste0(fasta_flags_str, " \\\n"),
    paste0("    ", paste(step_flags, collapse = " \\\n    ")), "\n\n",
    'EXIT_CODE=$?\n',
    'echo "Step 2 task ${SLURM_ARRAY_TASK_ID} finished with exit code: $EXIT_CODE"\n',
    'echo "Completed: $(date)"\n',
    'exit $EXIT_CODE\n'
  )

  # =========================================================================
  # Step 3 — Empirical Library Assembly (single job)
  # =========================================================================
  step3_script <- paste0(
    sbatch_header("s3_assembly", assembly_cpus, assembly_mem, assembly_time,
                  step_partition = asm_part, step_account = asm_acct,
                  requeue = asm_requeue), "\n\n",
    "module load apptainer\n\n",
    sprintf('echo "Step 3/5: Empirical Library Assembly for %s"\n', analysis_name),
    'echo "Started: $(date)"\n\n',
    quant_verify_block("quant_step2", 2),
    # Backup Step 2 quant files before assembly — Step 3 overwrites them
    # (same --temp dir + --quant-ori-names = same filenames). Backup enables
    # smart resume from Step 3 without re-running Step 2.
    sprintf('echo "Backing up Step 2 quant files..."\n'),
    sprintf('cp -r %s/quant_step2 %s/quant_step2_orig\n', output_dir, output_dir),
    sprintf('echo "Backup saved to %s/quant_step2_orig/"\n\n', output_dir),
    sprintf('%s \\\n', apptainer_cmd(full_bind_mount)),
    paste0(all_f_flags, " \\\n"),
    if (nzchar(fasta_flags_str)) paste0(fasta_flags_str, " \\\n"),
    sprintf('    --lib %s \\\n', predicted_lib),
    '    --use-quant \\\n',
    '    --quant-ori-names \\\n',
    '    --rt-profiling \\\n',
    '    --gen-spec-lib \\\n',
    sprintf('    --out-lib /work/out/empirical.parquet \\\n'),
    sprintf('    --temp /work/out/quant_step2 \\\n'),
    sprintf('    --out /work/out/step3_assembly.parquet \\\n'),
    sprintf('    --threads %d \\\n', assembly_cpus),
    paste0("    ", paste(step_flags, collapse = " \\\n    ")), "\n\n",
    'EXIT_CODE=$?\n',
    'echo "Step 3 finished with exit code: $EXIT_CODE"\n',
    'echo "Completed: $(date)"\n',
    'exit $EXIT_CODE\n'
  )

  # =========================================================================
  # Step 4 — Final per-file pass (SLURM array)
  # =========================================================================
  array_spec_4 <- sprintf("0-%d%%%d", n_files - 1, max_simultaneous)

  step4_script <- paste0(
    sbatch_header("s4_finalpass", cpus_per_file, mem_per_file,
                  time_per_file,
                  array_spec = array_spec_4, step_partition = arr_part,
                  step_account = arr_acct, requeue = arr_requeue), "\n\n",
    "module load apptainer\n\n",
    sprintf('echo "Step 4/5: Final-pass file ${SLURM_ARRAY_TASK_ID} of %d"\n', n_files),
    'echo "Started: $(date)"\n\n',
    '# Read file path from file list\n',
    sprintf('FILE_LIST="%s/file_list.txt"\n', output_dir),
    'RAW_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$FILE_LIST")\n',
    'FILE_DIR=$(dirname "$RAW_FILE")\n',
    'FILE_BASE=$(basename "$RAW_FILE")\n\n',
    'if [ -z "$RAW_FILE" ]; then\n',
    '  echo "ERROR: No file found for array task $SLURM_ARRAY_TASK_ID"\n',
    '  exit 1\n',
    'fi\n\n',
    '# Skip files that failed Step 2 — they will likely fail again\n',
    'QUANT_NAME="${FILE_BASE%.*}.quant"\n',
    sprintf('if [ ! -f "%s/quant_step2/${QUANT_NAME}" ]; then\n', output_dir),
    '  echo "SKIPPED: No Step 2 quantification for ${FILE_BASE} (${QUANT_NAME} not found)"\n',
    '  exit 0\n',
    'fi\n\n',
    sprintf('echo "Processing: $RAW_FILE"\n\n'),
    sprintf('%s \\\n', apptainer_cmd(perfile_bind_mount)),
    sprintf('    --f /work/data/$FILE_BASE \\\n'),
    sprintf('    --lib /work/out/empirical.parquet \\\n'),
    sprintf('    --temp /work/out/quant_step4 \\\n'),
    '    --no-ifs-removal \\\n',
    '    --quant-ori-names \\\n',
    sprintf('    --threads %d \\\n', cpus_per_file),
    if (nzchar(fasta_flags_str)) paste0(fasta_flags_str, " \\\n"),
    paste0("    ", paste(step_flags, collapse = " \\\n    ")), "\n\n",
    'EXIT_CODE=$?\n',
    'echo "Step 4 task ${SLURM_ARRAY_TASK_ID} finished with exit code: $EXIT_CODE"\n',
    'echo "Completed: $(date)"\n',
    'exit $EXIT_CODE\n'
  )

  # =========================================================================
  # Step 5 — Cross-run Report (single job)
  # =========================================================================
  norm_flag <- if (normalization == "off") "--no-norm" else ""

  step5_script <- paste0(
    sbatch_header("s5_report", assembly_cpus, assembly_mem, assembly_time,
                  step_partition = asm_part, step_account = asm_acct,
                  requeue = asm_requeue), "\n\n",
    "module load apptainer\n\n",
    sprintf('echo "Step 5/5: Cross-run Report for %s"\n', analysis_name),
    'echo "Started: $(date)"\n\n',
    quant_verify_block("quant_step4", 4),
    sprintf('%s \\\n', apptainer_cmd(full_bind_mount)),
    paste0(all_f_flags, " \\\n"),
    if (nzchar(fasta_flags_str)) paste0(fasta_flags_str, " \\\n"),
    '    --lib /work/out/empirical.parquet \\\n',
    '    --use-quant \\\n',
    '    --quant-ori-names \\\n',
    sprintf('    --temp /work/out/quant_step4 \\\n'),
    '    --matrices \\\n',
    sprintf('    --out /work/out/%s \\\n', report_name),
    sprintf('    --threads %d \\\n', assembly_cpus),
    if (nzchar(norm_flag)) paste0("    ", norm_flag, " \\\n"),
    paste0("    ", paste(step_flags, collapse = " \\\n    ")), "\n\n",
    'EXIT_CODE=$?\n',
    'echo "Step 5 finished with exit code: $EXIT_CODE"\n',
    'echo "Completed: $(date)"\n',
    'exit $EXIT_CODE\n'
  )

  list(
    step1_library   = step1_script,
    step2_firstpass = step2_script,
    step3_assembly  = step3_script,
    step4_finalpass = step4_script,
    step5_report    = step5_script
  )
}

# =============================================================================
# Resume Launcher for Failed Parallel Search
# =============================================================================

#' Generate a resume launcher script for a failed parallel DIA-NN search
#'
#' When a parallel search fails at Step N, this generates a bash script that
#' skips Steps 1..(N-1) and submits only Steps N..5 with dependency chaining.
#' Skipped steps output "STEP<n>:skipped" for the R-side parser.
#'
#' @param resume_from Integer 1-5 — which step to resume from
#' @param sbatch_bin Character — full path to sbatch binary
#' @param step_script_paths Character vector of length 5 — remote paths to sbatch scripts
#' @return Character — bash script content
generate_resume_launcher <- function(resume_from, sbatch_bin, step_script_paths) {
  resume_from <- as.integer(resume_from)
  stopifnot(length(resume_from) == 1, resume_from >= 1L, resume_from <= 5L,
            length(step_script_paths) == 5)
  lines <- c("#!/bin/bash", "set -e", "")

  # Mark skipped steps
  for (s in seq_len(resume_from - 1)) {
    lines <- c(lines, sprintf('echo "STEP%d:skipped"', s))
  }

  # Submit remaining steps with dependency chaining
  prev_var <- NULL
  for (s in resume_from:5) {
    var <- sprintf("JOB%d", s)
    if (is.null(prev_var)) {
      # First submitted step — no dependency
      lines <- c(lines,
        sprintf('%s=$(%s %s 2>&1)', var, sbatch_bin, step_script_paths[s]),
        sprintf('%s_ID=$(echo "$%s" | grep -oP "[0-9]+$")', var, var),
        sprintf('echo "STEP%d:$%s_ID"', s, var), "")
    } else {
      # Chain to previous — use afterany for Step 2→3 and Step 4→5 so a few
      # OOM/timeout tasks don't collapse the pipeline (verify blocks handle it)
      dep_type <- if ((s == 3 && prev_var == "JOB2") ||
                     (s == 5 && prev_var == "JOB4")) "afterany" else "afterok"
      lines <- c(lines,
        sprintf('%s=$(%s --dependency=%s:$%s_ID --kill-on-invalid-dep=yes %s 2>&1)',
                var, sbatch_bin, dep_type, prev_var, step_script_paths[s]),
        sprintf('%s_ID=$(echo "$%s" | grep -oP "[0-9]+$")', var, var),
        sprintf('echo "STEP%d:$%s_ID"', s, var), "")
    }
    prev_var <- var
  }

  paste(lines, collapse = "\n")
}

# =============================================================================
# Cluster Resource Monitoring
# =============================================================================

#' Check SLURM cluster resource utilization
#'
#' Queries SLURM for group CPU allocation, usage, and partition stats.
#' Supports both SSH and local SLURM access.
#'
#' @param ssh_config list(host, user, port, key_path, modules) or NULL for local
#' @param account Character — SLURM account name (e.g., "genome-center-grp")
#' @param partition Character — SLURM partition name (e.g., "high")
#' @param sbatch_path Character or NULL — full path to sbatch (to derive other SLURM tool paths)
#' @return list(success, group_limit, group_used, group_available,
#'              partition_idle, partition_total, error)
check_cluster_resources <- function(ssh_config, account, partition, sbatch_path = NULL) {
  # Derive SLURM command paths from sbatch path
  slurm_cmd <- function(cmd) {
    if (!is.null(sbatch_path) && nzchar(sbatch_path)) {
      file.path(dirname(sbatch_path), cmd)
    } else {
      cmd
    }
  }

  run_cmd <- function(command) {
    if (!is.null(ssh_config)) {
      ssh_exec(ssh_config, command, login_shell = is.null(sbatch_path), timeout = 15)
    } else if (slurm_proxy_available()) {
      slurm_proxy_exec(command, timeout = 15)
    } else {
      parts <- strsplit(command, " ")[[1]]
      stdout <- tryCatch(
        system2(parts[1], args = parts[-1], stdout = TRUE, stderr = TRUE),
        error = function(e) structure(e$message, status = 1L)
      )
      list(status = attr(stdout, "status") %||% 0L, stdout = stdout)
    }
  }

  # Get username for per-user queries
  username <- if (!is.null(ssh_config)) ssh_config$user else Sys.info()[["user"]]

  result <- list(
    success = FALSE,
    # Account-level
    group_limit = NA_integer_, group_used = NA_integer_, group_available = NA_integer_,
    # Per-user (the real constraint)
    user_limit = NA_integer_, user_used = NA_integer_, user_available = NA_integer_,
    # Partition
    partition_idle = NA_integer_, partition_total = NA_integer_,
    error = NULL
  )

  # --- 1. Account + per-user CPU limits via sacctmgr ---
  # Limits are set on the QOS (e.g. genome-center-grp-high-qos), not on associations
  # QOS name follows pattern: {account}-{partition}-qos
  tryCatch({
    qos_name <- sprintf("%s-%s-qos", account, partition)
    sacctmgr_cmd <- sprintf(
      "%s show qos where name=%s format=GrpTRES%%80,MaxTRESPU%%80 --noheader --parsable2",
      slurm_cmd("sacctmgr"), qos_name
    )
    res <- run_cmd(sacctmgr_cmd)
    if (res$status == 0 && length(res$stdout) > 0) {
      all_lines <- trimws(res$stdout)
      all_lines <- all_lines[nzchar(all_lines)]

      # Format: "cpu=616,gres/gpu=0,mem=9856G|cpu=64,gres/gpu=0,mem=1T"
      #          GrpTRES                       | MaxTRESPU
      for (line in all_lines) {
        fields <- strsplit(line, "\\|")[[1]]
        # GrpTRES is first field — account-level limit
        if (length(fields) >= 1 && is.na(result$group_limit)) {
          grp_match <- regmatches(fields[1], regexpr("cpu=[0-9]+", fields[1]))
          if (length(grp_match) > 0 && nzchar(grp_match)) {
            result$group_limit <- as.integer(sub("cpu=", "", grp_match))
          }
        }
        # MaxTRESPU is second field — per-user limit
        if (length(fields) >= 2 && is.na(result$user_limit)) {
          user_match <- regmatches(fields[2], regexpr("cpu=[0-9]+", fields[2]))
          if (length(user_match) > 0 && nzchar(user_match)) {
            result$user_limit <- as.integer(sub("cpu=", "", user_match))
          }
        }
      }
    }
  }, error = function(e) NULL)

  # --- 2a. Account-wide CPUs in use via squeue ---
  tryCatch({
    squeue_cmd <- sprintf(
      "%s -A %s -t RUNNING -o \"%%C\" --noheader",
      slurm_cmd("squeue"), account
    )
    res <- run_cmd(squeue_cmd)
    if (res$status == 0 && length(res$stdout) > 0) {
      cpu_counts <- suppressWarnings(as.integer(trimws(res$stdout)))
      cpu_counts <- cpu_counts[!is.na(cpu_counts)]
      result$group_used <- if (length(cpu_counts) > 0) sum(cpu_counts) else 0L
    } else {
      result$group_used <- 0L
    }
  }, error = function(e) {
    result$group_used <<- 0L
  })

  # --- 2b. THIS USER's CPUs in use on this account ---
  tryCatch({
    squeue_cmd <- sprintf(
      "%s -u %s -A %s -t RUNNING -o \"%%C\" --noheader",
      slurm_cmd("squeue"), username, account
    )
    res <- run_cmd(squeue_cmd)
    if (res$status == 0 && length(res$stdout) > 0) {
      cpu_counts <- suppressWarnings(as.integer(trimws(res$stdout)))
      cpu_counts <- cpu_counts[!is.na(cpu_counts)]
      result$user_used <- if (length(cpu_counts) > 0) sum(cpu_counts) else 0L
    } else {
      result$user_used <- 0L
    }
  }, error = function(e) {
    result$user_used <<- 0L
  })

  # --- 3. Partition stats via sinfo ---
  tryCatch({
    sinfo_cmd <- sprintf(
      "%s -p %s -o \"%%C\" --noheader",
      slurm_cmd("sinfo"), partition
    )
    res <- run_cmd(sinfo_cmd)
    if (res$status == 0 && length(res$stdout) > 0) {
      # sinfo %C format: Allocated/Idle/Other/Total
      line <- trimws(res$stdout[1])
      parts <- strsplit(line, "/")[[1]]
      if (length(parts) == 4) {
        result$partition_idle <- as.integer(parts[2])
        result$partition_total <- as.integer(parts[4])
      }
    }
  }, error = function(e) NULL)

  # Compute available CPUs
  if (!is.na(result$group_limit) && !is.na(result$group_used)) {
    result$group_available <- result$group_limit - result$group_used
  }
  if (!is.na(result$user_limit) && !is.na(result$user_used)) {
    result$user_available <- result$user_limit - result$user_used
  }

  # --- 4. Queue wait time: average wait for PENDING jobs (excluding dependency) ---
  result$pending_count <- 0L
  result$avg_wait_min <- NA_real_
  result$max_wait_min <- NA_real_
  tryCatch({
    # %V = submit time, %r = reason; filter to current user, exclude Dependency
    squeue_cmd <- sprintf(
      "%s -u %s -A %s -p %s -t PENDING -o \"%%V|%%r\" --noheader",
      slurm_cmd("squeue"), username, account, partition
    )
    res_q <- run_cmd(squeue_cmd)
    if (res_q$status == 0 && length(res_q$stdout) > 0) {
      lines <- trimws(res_q$stdout)
      lines <- lines[nzchar(lines)]
      if (length(lines) > 0) {
        # Parse submit_time|reason, keep only non-dependency pending
        now <- Sys.time()
        wait_mins <- numeric()
        for (ln in lines) {
          parts <- strsplit(ln, "\\|")[[1]]
          if (length(parts) >= 2) {
            reason <- trimws(parts[2])
            # Skip dependency-pending jobs
            if (grepl("Depend", reason, ignore.case = TRUE)) next
            submit_str <- trimws(parts[1])
            submit_time <- tryCatch(
              as.POSIXct(submit_str, format = "%Y-%m-%dT%H:%M:%S"),
              error = function(e) NA)
            if (!is.na(submit_time)) {
              wait_mins <- c(wait_mins, as.numeric(difftime(now, submit_time, units = "mins")))
            }
          }
        }
        result$pending_count <- length(wait_mins)
        if (length(wait_mins) > 0) {
          result$avg_wait_min <- mean(wait_mins)
          result$max_wait_min <- max(wait_mins)
        }
      }
    }
  }, error = function(e) NULL)

  # Success if we got at least squeue data
  if (!is.na(result$group_used) || !is.na(result$user_used)) {
    result$success <- TRUE
  } else {
    result$error <- "Could not query SLURM cluster resources"
  }

  result
}

#' Move a pending SLURM job to a different account/partition via scontrol
#'
#' @param job_id SLURM job ID (or array job ID)
#' @param new_account Target account
#' @param new_partition Target partition
#' @param ssh_config SSH config list or NULL
#' @param sbatch_path Path to sbatch (used to derive scontrol path)
#' @return list(success, message)
slurm_move_job <- function(job_id, new_account, new_partition,
                            ssh_config = NULL, sbatch_path = NULL) {
  scontrol_cmd <- if (!is.null(sbatch_path) && nzchar(sbatch_path)) {
    file.path(dirname(sbatch_path), "scontrol")
  } else "scontrol"

  cmd <- sprintf('%s update jobid=%s Account=%s Partition=%s',
    scontrol_cmd, job_id, new_account, new_partition)

  res <- if (!is.null(ssh_config)) {
    ssh_exec(ssh_config, cmd, login_shell = is.null(sbatch_path), timeout = 15)
  } else {
    parts <- strsplit(cmd, " ")[[1]]
    stdout <- tryCatch(
      system2(parts[1], args = parts[-1], stdout = TRUE, stderr = TRUE),
      error = function(e) structure(e$message, status = 1L))
    list(status = attr(stdout, "status") %||% 0L, stdout = stdout)
  }

  if (res$status == 0) {
    list(success = TRUE, message = sprintf("Job %s moved to %s/%s", job_id, new_account, new_partition))
  } else {
    list(success = FALSE, message = sprintf("Failed to move job %s: %s", job_id,
      paste(res$stdout, collapse = " ")))
  }
}

#' Select the best SLURM account/partition based on current cluster utilization
#'
#' @param lab_resources Result from check_cluster_resources() for lab group
#' @param public_resources Result from check_cluster_resources() for publicgrp
#' @param peak_cpus Numeric — estimated peak CPU need for the job
#' @return list(account, partition, reason)
select_best_partition <- function(lab_resources, public_resources, peak_cpus = 64) {
  # Per-user limits are the real constraint (e.g. 64 CPUs per user on genome-center-grp/high)
  # Account limits (e.g. 616 CPUs) are shared across all lab members
  user_limit <- NA_integer_
  user_used <- 0L
  user_available <- NA_integer_
  group_limit <- NA_integer_
  group_used <- 0L

  if (!is.null(lab_resources) && isTRUE(lab_resources$success)) {
    user_limit <- lab_resources$user_limit %||% NA_integer_
    user_used <- lab_resources$user_used %||% 0L
    user_available <- lab_resources$user_available %||% NA_integer_
    group_limit <- lab_resources$group_limit %||% NA_integer_
    group_used <- lab_resources$group_used %||% 0L
  }

  pub_idle <- 0L
  if (!is.null(public_resources) && isTRUE(public_resources$success)) {
    pub_idle <- public_resources$partition_idle %||% 0L
  }

  # Use per-user limit if available, otherwise fall back to group limit
  effective_limit <- if (!is.na(user_limit)) user_limit else group_limit
  effective_used <- if (!is.na(user_limit)) user_used else group_used
  effective_available <- if (!is.na(user_available)) user_available
    else if (!is.na(effective_limit)) effective_limit - effective_used
    else NA_integer_

  min_useful_cpus <- min(peak_cpus, 16L)  # at least 1 array task worth

  has_limit_info <- !is.na(effective_limit) && effective_limit > 0
  has_capacity <- has_limit_info && !is.na(effective_available) && effective_available >= min_useful_cpus
  at_limit <- has_limit_info && !is.na(effective_available) && effective_available < min_useful_cpus
  pub_has_idle <- pub_idle >= min_useful_cpus

  # Format the limit label for messages
  limit_label <- if (!is.na(user_limit)) {
    sprintf("Your usage: %d/%d CPUs", user_used, user_limit)
  } else if (!is.na(group_limit)) {
    sprintf("Group: %d/%d CPUs", group_used, group_limit)
  } else "no limit info"

  if (has_capacity) {
    list(account = "genome-center-grp", partition = "high",
         reason = sprintf("%s (%d available)", limit_label, effective_available))
  } else if (at_limit && pub_has_idle) {
    list(account = "publicgrp", partition = "low",
         reason = sprintf("%s \u2014 at capacity. Public has %d idle CPUs, faster start",
                          limit_label, pub_idle))
  } else if (at_limit) {
    list(account = "genome-center-grp", partition = "high",
         reason = sprintf("%s \u2014 at capacity. Public also busy. Using priority queue",
                          limit_label))
  } else {
    list(account = "genome-center-grp", partition = "high",
         reason = "Lab group (no limit info available)")
  }
}

# =============================================================================
# Predicted Spectral Library Caching
# =============================================================================

#' Compute a cache key for a predicted spectral library
#'
#' The key is an MD5 hash of the FASTA file basenames (sorted, path-independent)
#' plus the subset of search parameters that affect library prediction.
#'
#' @param fasta_files Character vector of FASTA file paths
#' @param search_params Named list of search parameters
#' @param search_mode Character: "libfree", "phospho", etc.
#' @param custom_fasta_text Character or NULL: raw custom FASTA sequences text
#'   (included in hash since the filename is always custom_proteins.fasta)
#' @return Character: MD5 hex digest
speclib_cache_key <- function(fasta_files, search_params, search_mode,
                               custom_fasta_text = NULL,
                               fasta_seq_count = NULL) {
  # Parameters that affect library prediction
  libpred_params <- c("enzyme", "missed_cleavages", "min_pep_len", "max_pep_len",
                       "min_pr_mz", "max_pr_mz", "min_pr_charge", "max_pr_charge",
                       "min_fr_mz", "max_fr_mz", "met_excision", "mod_met_ox",
                       "mod_nterm_acetyl", "extra_var_mods", "max_var_mods", "unimod4")
  sp_subset <- search_params[intersect(names(search_params), libpred_params)]
  # Sort for deterministic ordering
  sp_subset <- sp_subset[sort(names(sp_subset))]

  canonical <- list(
    fasta = sort(basename(fasta_files)),
    search_mode = search_mode,
    params = sp_subset,
    custom_fasta = if (!is.null(custom_fasta_text) && nzchar(trimws(custom_fasta_text))) {
      trimws(custom_fasta_text)
    },
    # Include sequence count to distinguish FASTAs with same name but different content
    fasta_seq_count = as.integer(fasta_seq_count)
  )
  digest::digest(canonical, algo = "md5")
}

#' Path to the speclib cache file
#'
#' Checks for shared proteomics volume first, falls back to user-local.
#' @return Character: file path to speclib_cache.rds
speclib_cache_path <- function() {
  # Allow env var override
  env_path <- Sys.getenv("DELIMP_SPECLIB_CACHE", "")
  if (nzchar(env_path)) return(env_path)

  # Shared proteomics volume (macOS local mount)
  shared_dir <- "/Volumes/proteomics-grp/dia-nn"
  if (dir.exists(shared_dir)) {
    return(file.path(shared_dir, "speclib_cache.rds"))
  }

  # HPC mount (Quobyte)
  hpc_dir <- "/quobyte/proteomics-grp/dia-nn"
  if (dir.exists(hpc_dir)) {
    return(file.path(hpc_dir, "speclib_cache.rds"))
  }

  # Fallback to user-local
  file.path(Sys.getenv("HOME"), ".delimp_speclib_cache.rds")
}

#' Check if the speclib cache is on a shared volume
#' @return Logical
speclib_cache_is_shared <- function() {
  path <- speclib_cache_path()
  grepl("^/Volumes/proteomics-grp/|^/quobyte/proteomics-grp/", path)
}

#' Load the speclib cache registry
#' @return List of cache entries
speclib_cache_load <- function() {
  path <- speclib_cache_path()
  if (file.exists(path)) {
    tryCatch(readRDS(path), error = function(e) {
      message("[DE-LIMP] Failed to read speclib cache: ", e$message)
      list()
    })
  } else {
    list()
  }
}

#' Save the speclib cache registry with file locking
#'
#' Uses a lockfile to prevent concurrent write corruption when multiple
#' users access the shared volume simultaneously.
#' @param cache List of cache entries
#' @return Logical: TRUE on success
speclib_cache_save <- function(cache) {
  path <- speclib_cache_path()
  lock_path <- paste0(path, ".lock")

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  # File-based locking: retry up to 10 times with 0.5s delay
  for (attempt in seq_len(10)) {
    if (!file.exists(lock_path)) {
      result <- tryCatch({
        writeLines(as.character(Sys.getpid()), lock_path)
        # Write atomically: temp file then rename
        tmp_path <- paste0(path, ".tmp.", Sys.getpid())
        saveRDS(cache, tmp_path)
        file.rename(tmp_path, path)
        unlink(lock_path)
        TRUE
      }, error = function(e) {
        unlink(lock_path)
        message("[DE-LIMP] speclib cache save failed: ", e$message)
        FALSE
      })
      return(result)
    }
    Sys.sleep(0.5)
  }

  message("[DE-LIMP] speclib cache save failed: could not acquire lock after 5 seconds")
  FALSE
}

#' Migrate user-local speclib cache entries to shared volume
#'
#' If a shared volume is available and a local cache exists, merge local entries
#' into the shared cache and remove the local file. Called once at startup.
#' @return Invisible NULL
speclib_cache_migrate <- function() {
  local_path <- file.path(Sys.getenv("HOME"), ".delimp_speclib_cache.rds")
  if (!file.exists(local_path)) return(invisible(NULL))
  if (!speclib_cache_is_shared()) return(invisible(NULL))

  tryCatch({
    local_cache <- readRDS(local_path)
    if (length(local_cache) == 0) {
      unlink(local_path)
      return(invisible(NULL))
    }

    shared_cache <- speclib_cache_load()
    shared_keys <- vapply(shared_cache, function(e) e$key %||% "", character(1))

    # Merge: add local entries not already in shared cache
    n_added <- 0L
    for (entry in local_cache) {
      if (!is.null(entry$key) && !(entry$key %in% shared_keys)) {
        shared_cache <- c(shared_cache, list(entry))
        n_added <- n_added + 1L
      }
    }

    if (n_added > 0L) {
      speclib_cache_save(shared_cache)
      message(sprintf("[DE-LIMP] Migrated %d speclib cache entries to shared volume", n_added))
    }

    # Remove local file after successful migration
    unlink(local_path)
  }, error = function(e) {
    message("[DE-LIMP] speclib cache migration failed: ", e$message)
  })
  invisible(NULL)
}

#' Register a predicted spectral library in the cache
#'
#' @param fasta_files Character vector of FASTA paths
#' @param search_params Named list of search parameters
#' @param search_mode Character: search mode
#' @param speclib_path Character: full path to the .predicted.speclib file
#' @param analysis_name Character: name of the analysis that produced this library
#' @param output_dir Character: output directory of the analysis
speclib_cache_register <- function(fasta_files, search_params, search_mode,
                                    speclib_path, analysis_name, output_dir,
                                    custom_fasta_text = NULL,
                                    fasta_seq_count = NULL) {
  key <- speclib_cache_key(fasta_files, search_params, search_mode,
                           custom_fasta_text, fasta_seq_count)

  entry <- list(
    key = key,
    speclib_path = speclib_path,
    fasta_files = basename(fasta_files),
    search_mode = search_mode,
    params_subset = search_params[intersect(names(search_params),
      c("enzyme", "missed_cleavages", "min_pep_len", "max_pep_len",
        "min_pr_mz", "max_pr_mz", "min_pr_charge", "max_pr_charge",
        "min_fr_mz", "max_fr_mz", "met_excision", "mod_met_ox",
        "mod_nterm_acetyl", "extra_var_mods", "max_var_mods", "unimod4"))],
    analysis_name = analysis_name,
    created_at = Sys.time(),
    output_dir = output_dir,
    fasta_seq_count = as.integer(fasta_seq_count),
    registered_by = Sys.getenv("USER", "unknown")
  )

  # Re-read cache inside save to minimize race window on shared volumes
  cache <- speclib_cache_load()
  # Deduplicate on key (newer wins)
  cache <- Filter(function(e) e$key != key, cache)
  cache <- c(cache, list(entry))
  speclib_cache_save(cache)
}

#' Look up a cached predicted spectral library
#'
#' @param fasta_files Character vector of FASTA paths
#' @param search_params Named list of search parameters
#' @param search_mode Character: search mode
#' @return Matching cache entry (list) or NULL
speclib_cache_lookup <- function(fasta_files, search_params, search_mode,
                                  custom_fasta_text = NULL,
                                  fasta_seq_count = NULL) {
  key <- speclib_cache_key(fasta_files, search_params, search_mode,
                           custom_fasta_text, fasta_seq_count)
  cache <- speclib_cache_load()
  for (entry in cache) {
    if (identical(entry$key, key)) return(entry)
  }
  NULL
}

# =============================================================================
# Shared FASTA Database Library
# =============================================================================

#' Get the base path for the shared FASTA library
#'
#' Checks for the shared proteomics volume first, falls back to a local
#' user directory if the shared volume is not mounted.
#'
#' @return Character: directory path to the FASTA library root
fasta_library_path <- function() {
  # Allow env var override for custom deployments

  env_path <- Sys.getenv("DELIMP_FASTA_LIBRARY", "")
  if (nzchar(env_path) && dir.exists(env_path)) return(env_path)

  # Shared proteomics volume (macOS local mount)
  shared_path <- "/Volumes/proteomics-grp/dia-nn/fasta_library"
  if (dir.exists(shared_path)) return(shared_path)

  # HPC mount (Quobyte)
  hpc_path <- "/quobyte/proteomics-grp/dia-nn/fasta_library"
  if (dir.exists(hpc_path)) return(hpc_path)

  # Fallback to user-local directory

  local_path <- file.path(Sys.getenv("HOME"), ".delimp_fasta_library")
  if (!dir.exists(local_path)) {
    dir.create(local_path, recursive = TRUE, showWarnings = FALSE)
  }
  local_path
}

#' Check if the shared FASTA library volume is available
#' @return Logical: TRUE if on shared volume, FALSE if using local fallback
fasta_library_is_shared <- function() {
  env_path <- Sys.getenv("DELIMP_FASTA_LIBRARY", "")
  if (nzchar(env_path) && dir.exists(env_path)) return(TRUE)

  dir.exists("/Volumes/proteomics-grp/dia-nn/fasta_library") ||
    dir.exists("/quobyte/proteomics-grp/dia-nn/fasta_library")
}

#' Get the HPC-equivalent remote path for a library entry
#'
#' Translates the local macOS mount path to the HPC Quobyte path.
#' @param local_dir Character: local path (e.g., /Volumes/proteomics-grp/...)
#' @return Character: HPC path (e.g., /quobyte/proteomics-grp/...)
fasta_library_remote_path <- function(local_dir) {
  # macOS mount -> HPC Quobyte
  gsub("^/Volumes/proteomics-grp/", "/quobyte/proteomics-grp/", local_dir)
}

#' Load the FASTA library catalog
#'
#' Reads catalog.rds from the library path. Returns an empty list if the
#' file doesn't exist or is corrupted.
#'
#' @return List of catalog entries (each is a named list per the schema)
fasta_library_load <- function() {
  catalog_path <- file.path(fasta_library_path(), "catalog.rds")
  if (!file.exists(catalog_path)) return(list())

  tryCatch({
    catalog <- readRDS(catalog_path)
    if (!is.list(catalog)) return(list())
    catalog
  }, error = function(e) {
    message(sprintf("[DE-LIMP] FASTA library catalog load failed: %s", e$message))
    list()
  })
}

#' Save the FASTA library catalog with file locking
#'
#' Uses a lockfile to prevent concurrent write corruption when multiple
#' users access the shared volume simultaneously.
#'
#' @param catalog List of catalog entries
#' @return Logical: TRUE on success
fasta_library_save <- function(catalog) {
  lib_path <- fasta_library_path()
  catalog_path <- file.path(lib_path, "catalog.rds")
  lock_path <- paste0(catalog_path, ".lock")

  # Ensure library directory exists
  if (!dir.exists(lib_path)) {
    dir.create(lib_path, recursive = TRUE, showWarnings = FALSE)
  }

  # Simple file-based locking: create lockfile, write, remove lock
  # Retry up to 10 times with 0.5s delay if lock is held
  for (attempt in seq_len(10)) {
    if (!file.exists(lock_path)) {
      result <- tryCatch({
        # Create lock
        writeLines(as.character(Sys.getpid()), lock_path)

        # Write catalog atomically: write to temp, then rename
        tmp_path <- paste0(catalog_path, ".tmp")
        saveRDS(catalog, tmp_path)
        file.rename(tmp_path, catalog_path)

        # Remove lock
        unlink(lock_path)
        TRUE
      }, error = function(e) {
        unlink(lock_path)
        message(sprintf("[DE-LIMP] FASTA library save failed: %s", e$message))
        FALSE
      })
      return(result)
    }
    Sys.sleep(0.5)
  }

  message("[DE-LIMP] FASTA library save failed: could not acquire lock after 5 seconds")
  FALSE
}

#' Add an entry to the FASTA library catalog
#'
#' Deduplicates on entry id. If an entry with the same id already exists,
#' it is replaced with the new one.
#'
#' @param entry Named list: catalog entry conforming to the schema
#' @return Logical: TRUE on success
fasta_library_add <- function(entry) {
  if (is.null(entry$id) || !nzchar(entry$id)) {
    message("[DE-LIMP] FASTA library add: entry missing id")
    return(FALSE)
  }

  catalog <- fasta_library_load()

  # Remove any existing entry with the same id (dedup)
  catalog <- Filter(function(e) !identical(e$id, entry$id), catalog)
  catalog <- c(catalog, list(entry))

  fasta_library_save(catalog)
}

#' Update fields on an existing FASTA library catalog entry
#'
#' @param id Character: entry id to update
#' @param updates Named list of fields to merge into the entry
#' @return Logical: TRUE on success
fasta_library_update_entry <- function(id, updates) {
  catalog <- fasta_library_load()
  idx <- which(vapply(catalog, function(e) identical(e$id, id), logical(1)))
  if (length(idx) == 0) return(FALSE)
  for (nm in names(updates)) catalog[[idx[1]]][[nm]] <- updates[[nm]]
  fasta_library_save(catalog)
}

#' Parse DIA-NN log output to extract actual flags used
#'
#' Reads a DIA-NN .out log file and extracts the command-line flags that
#' were actually used. Returns a named list of search parameters.
#'
#' @param log_lines Character vector of log file lines
#' @return Named list of verified search parameters
parse_diann_log_flags <- function(log_lines) {
  if (length(log_lines) == 0) return(list())

  all_text <- paste(log_lines, collapse = " ")

  # DIA-NN log has TWO sources of truth:
  # 1. The apptainer/diann command line in the sbatch echo (has --flag syntax)
  # 2. DIA-NN's own confirmation output (has "Max fragment m/z set to 1800" syntax)
  # We parse BOTH — DIA-NN's confirmation is authoritative

  # Helper: extract from DIA-NN confirmation lines like "Max fragment m/z set to 1800"
  extract_set_to <- function(pattern) {
    m <- regmatches(all_text, regexpr(sprintf("%s set to\\s+([0-9.]+)", pattern), all_text))
    if (length(m) > 0) as.numeric(sub(".*set to\\s+", "", m[1])) else NULL
  }

  # Helper: extract from CLI flags like --max-fr-mz 1800
  extract_cli_int <- function(flag) {
    m <- regmatches(all_text, regexpr(sprintf("%s\\s+(\\d+)", flag), all_text))
    if (length(m) > 0) as.integer(sub(sprintf("^%s\\s+", flag), "", m[1])) else NULL
  }

  has_text <- function(pattern) grepl(pattern, all_text, fixed = TRUE)

  list(
    enzyme = {
      # "In silico digest will involve cuts at K*,R*"
      m <- regmatches(all_text, regexpr("cuts at\\s+(\\S+)", all_text))
      if (length(m) > 0) sub("^cuts at\\s+", "", m[1])
      else {
        m2 <- regmatches(all_text, regexpr("--cut\\s+(\\S+)", all_text))
        if (length(m2) > 0) sub("^--cut\\s+", "", m2[1]) else NULL
      }
    },
    missed_cleavages = as.integer(extract_set_to("Maximum number of missed cleavages") %||%
                                    extract_cli_int("--missed-cleavages")),
    min_pep_len = as.integer(extract_set_to("Min peptide length") %||%
                               extract_cli_int("--min-pep-len")),
    max_pep_len = as.integer(extract_set_to("Max peptide length") %||%
                               extract_cli_int("--max-pep-len")),
    min_pr_mz = as.integer(extract_set_to("Min precursor m/z") %||%
                             extract_cli_int("--min-pr-mz")),
    max_pr_mz = as.integer(extract_set_to("Max precursor m/z") %||%
                             extract_cli_int("--max-pr-mz")),
    min_fr_mz = as.integer(extract_set_to("Min fragment m/z") %||%
                             extract_cli_int("--min-fr-mz")),
    max_fr_mz = as.integer(extract_set_to("Max fragment m/z") %||%
                             extract_cli_int("--max-fr-mz")),
    min_pr_charge = as.integer(extract_set_to("Min precursor charge") %||%
                                 extract_cli_int("--min-pr-charge")),
    max_pr_charge = as.integer(extract_set_to("Max precursor charge") %||%
                                 extract_cli_int("--max-pr-charge")),
    mass_acc = {
      m <- regmatches(all_text, regexpr("--mass-acc\\s+([0-9.]+)", all_text))
      if (length(m) > 0 && !grepl("--mass-acc-ms1", m[1]))
        as.numeric(sub("^--mass-acc\\s+", "", m[1])) else NULL
    },
    mass_acc_ms1 = {
      m <- regmatches(all_text, regexpr("--mass-acc-ms1\\s+([0-9.]+)", all_text))
      if (length(m) > 0) as.numeric(sub("^--mass-acc-ms1\\s+", "", m[1])) else NULL
    },
    qvalue = {
      m <- regmatches(all_text, regexpr("filtered at\\s+([0-9.]+)\\s+FDR", all_text))
      if (length(m) > 0) {
        as.numeric(sub("\\s+FDR.*", "", sub(".*filtered at\\s+", "", m[1])))
      } else {
        m2 <- regmatches(all_text, regexpr("--qvalue\\s+([0-9.]+)", all_text))
        if (length(m2) > 0) as.numeric(sub("^--qvalue\\s+", "", m2[1])) else NULL
      }
    },
    max_var_mods = as.integer(extract_set_to("Maximum number of variable modifications") %||%
                                extract_cli_int("--var-mods")),
    scan_window = as.integer(extract_set_to("Scan window radius") %||%
                               extract_cli_int("--window")),
    mod_met_ox = has_text("UniMod:35") || has_text("Modification UniMod:35"),
    mod_nterm_acetyl = has_text("UniMod:1,") || has_text("UniMod:1 "),
    unimod4 = has_text("carbamidomethylation enabled") || has_text("--unimod4"),
    met_excision = has_text("methionine excision enabled") || has_text("--met-excision"),
    fasta_search = has_text("FASTA digest") || has_text("--fasta-search"),
    gen_spec_lib = has_text("spectral library will be generated") || has_text("--gen-spec-lib"),
    # Library generation stats
    n_precursors = {
      m <- regmatches(all_text, regexpr("([0-9]+)\\s+precursors generated", all_text))
      if (length(m) > 0) as.integer(sub("\\s+precursors.*", "", m[1])) else NULL
    },
    n_proteins_lib = {
      m <- regmatches(all_text, regexpr("Library contains\\s+([0-9]+)\\s+proteins", all_text))
      if (length(m) > 0) as.integer(sub(".*contains\\s+", "", sub("\\s+proteins.*", "", m[1]))) else NULL
    },
    n_genes_lib = {
      m <- regmatches(all_text, regexpr("and\\s+([0-9]+)\\s+genes", all_text))
      if (length(m) > 0) as.integer(sub(".*and\\s+", "", sub("\\s+genes.*", "", m[1]))) else NULL
    }
  )
}

#' Remove an entry from the FASTA library catalog
#'
#' Optionally deletes the associated files on disk.
#'
#' @param id Character: entry id to remove
#' @param delete_files Logical: also delete FASTA files from disk (default FALSE)
#' @return Logical: TRUE on success
fasta_library_remove <- function(id, delete_files = FALSE) {
  catalog <- fasta_library_load()

  # Find the entry before removing
  entry <- NULL
  for (e in catalog) {
    if (identical(e$id, id)) {
      entry <- e
      break
    }
  }

  if (is.null(entry)) return(FALSE)

  # Optionally delete files
  if (delete_files && !is.null(entry$fasta_dir)) {
    dir_path <- file.path(fasta_library_path(), entry$fasta_dir)
    if (dir.exists(dir_path)) {
      unlink(dir_path, recursive = TRUE)
    }
  }

  # Remove from catalog
  catalog <- Filter(function(e) !identical(e$id, id), catalog)
  fasta_library_save(catalog)
}

#' Check the age/freshness status of a FASTA library entry
#'
#' Based on `created_at` in the catalog entry. Uses 6-month default
#' (configurable via DELIMP_FASTA_MAX_AGE_DAYS env var).
#'
#' @param entry Named list: catalog entry
#' @return Character: "fresh" (< 5 months), "expiring" (5-6 months), "expired" (> 6 months)
fasta_library_check_age <- function(entry) {
  max_age_days <- as.integer(Sys.getenv("DELIMP_FASTA_MAX_AGE_DAYS", "180"))
  warning_days <- max_age_days - 30L  # 1 month before expiry

  created <- tryCatch(
    as.POSIXct(entry$created_at),
    error = function(e) NA
  )
  if (is.na(created)) return("expired")

  age_days <- as.numeric(difftime(Sys.time(), created, units = "days"))

  if (age_days > max_age_days) return("expired")
  if (age_days > warning_days) return("expiring")
  "fresh"
}

#' Compute human-readable age string for a library entry
#'
#' @param entry Named list: catalog entry
#' @return Character: e.g., "2 months", "15 days", "8 months"
fasta_library_age_label <- function(entry) {
  created <- tryCatch(
    as.POSIXct(entry$created_at),
    error = function(e) NA
  )
  if (is.na(created)) return("Unknown")

  age_days <- as.numeric(difftime(Sys.time(), created, units = "days"))

  if (age_days < 1) return("Today")
  if (age_days < 30) return(sprintf("%.0f days", age_days))
  months <- floor(age_days / 30)
  if (months == 1) return("1 month")
  sprintf("%d months", months)
}

#' Look up a library entry by id
#'
#' @param id Character: entry id
#' @return Named list (catalog entry) or NULL
fasta_library_lookup <- function(id) {
  catalog <- fasta_library_load()
  for (entry in catalog) {
    if (identical(entry$id, id)) return(entry)
  }
  NULL
}

#' Verify that FASTA files for a library entry actually exist on disk
#'
#' @param entry Named list: catalog entry
#' @return Logical: TRUE if all FASTA files exist
fasta_library_verify_files <- function(entry) {
  if (is.null(entry$fasta_dir) || is.null(entry$fasta_files)) return(FALSE)

  dir_path <- file.path(fasta_library_path(), entry$fasta_dir)
  if (!dir.exists(dir_path)) return(FALSE)

  all(vapply(entry$fasta_files, function(f) {
    file.exists(file.path(dir_path, f))
  }, logical(1)))
}

#' Get absolute file paths for a library entry's FASTA files
#'
#' @param entry Named list: catalog entry
#' @param use_remote Logical: return HPC remote paths instead of local
#' @return Character vector of absolute FASTA file paths
fasta_library_file_paths <- function(entry, use_remote = FALSE) {
  if (use_remote) {
    # Try remote_dir first, then translate local paths
    rd <- entry$remote_dir
    if (!is.null(rd)) {
      # Translate macOS shared volume mount → HPC Quobyte path
      if (grepl("^/Volumes/proteomics-grp/", rd)) {
        rd <- fasta_library_remote_path(rd)
      }
      # If path is still local-only (e.g. /Users/...), can't use it
      if (!grepl("^/quobyte/|^/share/|^/home/", rd)) {
        rd <- NULL
      }
    }
    if (!is.null(rd)) {
      return(file.path(rd, entry$fasta_files))
    }
    # Fallback: return local paths (caller must handle upload)
  }
  base_dir <- file.path(fasta_library_path(), entry$fasta_dir)
  file.path(base_dir, entry$fasta_files)
}

#' Check if FASTA paths are reachable on remote HPC
#' @return TRUE if all paths are valid HPC paths, FALSE if any are local-only
fasta_paths_are_remote <- function(fasta_files) {
  all(grepl("^/quobyte/|^/share/|^/home/", fasta_files))
}

#' Build a catalog entry from a UniProt download result
#'
#' Convenience function that creates a properly structured catalog entry
#' from the results of download_uniprot_fasta() plus UniProt metadata.
#'
#' @param download_result List from download_uniprot_fasta()
#' @param uniprot_row data.frame row from search_uniprot_proteomes()
#' @param content_type Character: "one_per_gene", "reviewed", "full", "full_isoforms"
#' @param contam_info List or NULL: contaminant FASTA info (from get_contaminant_fasta)
#' @param contam_name Character: contaminant library name (e.g., "universal")
#' @param custom_sequences Character or NULL: custom FASTA text
#' @param search_params List: search settings used
#' @param speclib_path Character or NULL: path to predicted spectral library
#' @param created_by Character: username
#' @param notes Character: optional notes
#' @return Named list: catalog entry
fasta_library_build_entry <- function(download_result, uniprot_row,
                                       content_type,
                                       contam_info = NULL,
                                       contam_name = "none",
                                       custom_sequences = NULL,
                                       search_params = list(),
                                       speclib_path = NULL,
                                       created_by = Sys.info()[["user"]],
                                       notes = "") {
  # Generate a unique id
  id <- paste0(
    tolower(gsub("[^A-Za-z0-9]", "", uniprot_row$organism)),
    "_", content_type, "_",
    format(Sys.time(), "%Y%m%d_%H%M%S")
  )

  # Build name
  type_label <- switch(content_type,
    "one_per_gene"  = "OPG",
    "reviewed"      = "Swiss-Prot",
    "full"          = "Full",
    "full_isoforms" = "Full+Iso",
    "Custom"
  )

  common <- uniprot_row$common_name %||% ""
  name_parts <- c(
    if (nzchar(common)) common else uniprot_row$organism,
    type_label
  )
  if (!is.null(contam_info) && contam_name != "none") {
    name_parts <- c(name_parts, paste0("+ ", tools::toTitleCase(gsub("_", " ", contam_name))))
  }
  if (!is.null(custom_sequences) && nzchar(trimws(custom_sequences %||% ""))) {
    n_custom <- sum(grepl("^>", strsplit(custom_sequences, "\n")[[1]]))
    if (n_custom > 0) name_parts <- c(name_parts, sprintf("+ %d custom", n_custom))
  }
  name <- paste(name_parts, collapse = " ")

  # Collect FASTA file basenames
  fasta_files <- basename(download_result$path)
  if (!is.null(contam_info) && isTRUE(contam_info$success)) {
    fasta_files <- c(fasta_files, basename(contam_info$path))
  }

  # Build directory name
  safe_org <- tolower(gsub("[^A-Za-z0-9]", "_", uniprot_row$organism))
  safe_org <- gsub("_+", "_", safe_org)
  safe_org <- substr(safe_org, 1, 20)
  fasta_dir <- sprintf("%s_%s_%s", safe_org, content_type, format(Sys.Date(), "%Y_%m"))

  # Custom sequence count
  custom_count <- 0L
  if (!is.null(custom_sequences) && nzchar(trimws(custom_sequences %||% ""))) {
    custom_count <- as.integer(sum(grepl("^>", strsplit(custom_sequences, "\n")[[1]])))
  }

  # Build HPC remote dir
  lib_path <- fasta_library_path()
  remote_dir <- fasta_library_remote_path(file.path(lib_path, fasta_dir))

  list(
    id = id,
    name = name,
    organism = uniprot_row$organism %||% "",
    organism_common = uniprot_row$common_name %||% "",
    proteome_id = uniprot_row$upid %||% "",
    content_type = content_type,
    protein_count = as.integer(download_result$n_sequences %||% 0L),
    file_size_bytes = as.integer(download_result$file_size %||% 0L),
    contaminant_library = if (contam_name != "none") tools::toTitleCase(gsub("_", " ", contam_name)) else NULL,
    contaminant_count = as.integer(contam_info$n_sequences %||% 0L),
    custom_sequences = custom_sequences,
    custom_sequence_count = custom_count,
    fasta_files = fasta_files,
    fasta_dir = fasta_dir,
    search_settings = list(
      enzyme = search_params$enzyme %||% "K*,R*",
      missed_cleavages = as.integer(search_params$missed_cleavages %||% 1L),
      var_mods = paste(Filter(nzchar, c(
        if (isTRUE(search_params$mod_met_ox)) "UniMod:35 (Met oxidation)" else NULL,
        if (isTRUE(search_params$mod_nterm_acetyl)) "UniMod:1 (N-term acetylation)" else NULL,
        if (nzchar(search_params$extra_var_mods %||% "")) search_params$extra_var_mods else NULL
      )), collapse = "; "),
      fixed_mods = if (isTRUE(search_params$unimod4)) "UniMod:4 (Carbamidomethylation)" else "",
      min_pep_len = as.integer(search_params$min_pep_len %||% 7L),
      max_pep_len = as.integer(search_params$max_pep_len %||% 30L),
      min_pr_mz = as.numeric(search_params$min_pr_mz %||% 300),
      max_pr_mz = as.numeric(search_params$max_pr_mz %||% 1800),
      min_fr_mz = as.numeric(search_params$min_fr_mz %||% 200),
      max_fr_mz = as.numeric(search_params$max_fr_mz %||% 1800)
    ),
    speclib_path = speclib_path,
    speclib_search_mode = if (!is.null(speclib_path)) "libfree" else NULL,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    created_by = created_by,
    notes = notes,
    remote_dir = remote_dir
  )
}

#' Build a display data.frame from the FASTA library catalog
#'
#' Creates a data.frame suitable for rendering in a DT table.
#'
#' @param catalog List of catalog entries
#' @return data.frame with display columns
fasta_library_display_df <- function(catalog) {
  if (length(catalog) == 0) {
    return(data.frame(
      Name = character(), Organism = character(), Proteins = integer(),
      Age = character(), Status = character(), `Created By` = character(),
      id = character(), stringsAsFactors = FALSE, check.names = FALSE
    ))
  }

  data.frame(
    Name = vapply(catalog, function(e) e$name %||% "", character(1)),
    Organism = vapply(catalog, function(e) e$organism %||% "", character(1)),
    Proteins = vapply(catalog, function(e) {
      n <- as.integer(e$protein_count %||% 0L)
      contam <- as.integer(e$contaminant_count %||% 0L)
      n + contam
    }, integer(1)),
    Age = vapply(catalog, fasta_library_age_label, character(1)),
    Status = vapply(catalog, fasta_library_check_age, character(1)),
    `Created By` = vapply(catalog, function(e) e$created_by %||% "", character(1)),
    id = vapply(catalog, function(e) e$id %||% "", character(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

# =============================================================================
# Unified Activity Log — replaces search_history + analysis_history + projects
# =============================================================================

activity_log_headers <- c(
  "id", "event_type", "timestamp", "user", "search_name", "project", "notes",
  "backend", "search_mode", "parallel", "n_files", "fasta_files", "fasta_seq_count",
  "normalization", "enzyme", "mass_acc_mode", "mass_acc", "mass_acc_ms1",
  "scan_window", "mbr", "extra_cli_flags",
  "output_dir", "job_id", "status", "duration_min",
  "n_proteins", "n_samples", "n_contrasts", "n_de_proteins",
  "session_file", "speclib_cached", "app_version", "source_type"
)

#' Get path for the unified activity log CSV
#' Always local — never depends on mounted network drives (may be absent/slow).
#' Cross-user sharing is handled via SSH/SCP sync when connected to HPC.
activity_log_path <- function() {
  file.path(Sys.getenv("HOME"), ".delimp_activity_log.csv")
}

#' Read the activity log CSV
activity_log_read <- function(path = activity_log_path()) {
  if (!file.exists(path)) return(data.frame())
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) data.frame())
}

#' Generate a simple unique ID (no uuid dependency)
generate_activity_id <- function() {
  paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample(999999, 1)))
}

#' Record an activity event (append-only)
record_activity <- function(entry, path = activity_log_path()) {
  if (is.null(entry$id)) entry$id <- generate_activity_id()
  row <- as.data.frame(
    lapply(activity_log_headers, function(h) entry[[h]] %||% NA),
    stringsAsFactors = FALSE
  )
  names(row) <- activity_log_headers

  lock_path <- paste0(path, ".lock")
  lock <- filelock::lock(lock_path, timeout = 5000)
  on.exit(filelock::unlock(lock), add = TRUE)

  needs_header <- !file.exists(path)
  tryCatch({
    suppressWarnings(
      write.table(row, file = path, append = TRUE, sep = ",", row.names = FALSE,
        col.names = needs_header, quote = TRUE)
    )
  }, error = function(e) message("[DE-LIMP] Failed to write activity log: ", e$message))
}

#' Update an activity log row by output_dir (read-modify-write)
update_activity <- function(output_dir, updates, event_type_filter = NULL,
                            path = activity_log_path()) {
  if (!file.exists(path)) return(invisible(NULL))

  lock_path <- paste0(path, ".lock")
  lock <- filelock::lock(lock_path, timeout = 5000)
  on.exit(filelock::unlock(lock), add = TRUE)

  tryCatch({
    log <- read.csv(path, stringsAsFactors = FALSE)
    if (nrow(log) == 0) return(invisible(NULL))

    idx <- which(log$output_dir == output_dir)
    if (!is.null(event_type_filter))
      idx <- idx[log$event_type[idx] == event_type_filter]
    if (length(idx) == 0) return(invisible(NULL))
    idx <- idx[length(idx)]  # most recent match

    for (nm in names(updates)) {
      if (nm %in% names(log)) log[[nm]][idx] <- updates[[nm]]
    }
    write.csv(log, file = path, row.names = FALSE, quote = TRUE)
  }, error = function(e) message("[DE-LIMP] Failed to update activity log: ", e$message))
}

#' Update an activity log row by id
update_activity_by_id <- function(id, updates, path = activity_log_path()) {
  if (!file.exists(path)) return(invisible(NULL))

  lock_path <- paste0(path, ".lock")
  lock <- filelock::lock(lock_path, timeout = 5000)
  on.exit(filelock::unlock(lock), add = TRUE)

  tryCatch({
    log <- read.csv(path, stringsAsFactors = FALSE)
    if (nrow(log) == 0) return(invisible(NULL))

    idx <- which(log$id == id)
    if (length(idx) == 0) return(invisible(NULL))
    idx <- idx[1]

    for (nm in names(updates)) {
      if (nm %in% names(log)) log[[nm]][idx] <- updates[[nm]]
    }
    write.csv(log, file = path, row.names = FALSE, quote = TRUE)
  }, error = function(e) message("[DE-LIMP] Failed to update activity by id: ", e$message))
}

#' Get unique project names from activity log
get_projects <- function(path = activity_log_path()) {
  log <- activity_log_read(path)
  if (nrow(log) == 0 || !"project" %in% names(log)) return(character(0))
  sort(unique(na.omit(log$project[nzchar(log$project)])))
}

#' Set project for all rows matching an output_dir
set_project <- function(output_dir, project_name, path = activity_log_path()) {
  if (!file.exists(path)) return(invisible(NULL))

  lock_path <- paste0(path, ".lock")
  lock <- filelock::lock(lock_path, timeout = 5000)
  on.exit(filelock::unlock(lock), add = TRUE)

  tryCatch({
    log <- read.csv(path, stringsAsFactors = FALSE)
    if (nrow(log) == 0) return(invisible(NULL))
    idx <- which(!is.na(log$output_dir) & log$output_dir == output_dir)
    if (length(idx) > 0) {
      log$project[idx] <- project_name
      write.csv(log, file = path, row.names = FALSE, quote = TRUE)
    }
  }, error = function(e) message("[DE-LIMP] Failed to set project: ", e$message))
}

#' Migrate old search_history + analysis_history CSVs + projects.json to unified activity log
migrate_to_activity_log <- function() {
  new_path <- activity_log_path()
  if (file.exists(new_path)) return(invisible(FALSE))  # already migrated

  # Old paths (local only — no mounted drive dependency)
  sh_local <- file.path(Sys.getenv("HOME"), ".delimp_search_history.csv")
  ah_local <- file.path(Sys.getenv("HOME"), ".delimp_analysis_history.csv")
  pj_local <- file.path(Sys.getenv("HOME"), ".delimp_projects.json")

  sh_path <- if (file.exists(sh_local)) sh_local else NULL
  ah_path <- if (file.exists(ah_local)) ah_local else NULL
  pj_path <- if (file.exists(pj_local)) pj_local else NULL

  if (is.null(sh_path) && is.null(ah_path)) return(invisible(FALSE))  # nothing to migrate

  # Build project lookup: output_dir -> project_name
  proj_map <- list()
  if (!is.null(pj_path)) {
    tryCatch({
      pj <- jsonlite::fromJSON(pj_path, simplifyVector = FALSE)
      for (pname in names(pj$projects)) {
        for (od in unlist(pj$projects[[pname]]$entries)) {
          proj_map[[od]] <- pname
        }
      }
    }, error = function(e) NULL)
  }

  rows <- list()

  # Migrate search history rows
  if (!is.null(sh_path)) {
    tryCatch({
      sh <- read.csv(sh_path, stringsAsFactors = FALSE)
      for (i in seq_len(nrow(sh))) {
        evt <- if (!is.na(sh$status[i]) && sh$status[i] == "completed") "search_completed"
               else if (!is.na(sh$status[i]) && sh$status[i] == "failed") "search_failed"
               else "search_submitted"
        od <- sh$output_dir[i] %||% NA
        rows[[length(rows) + 1]] <- list(
          id = generate_activity_id(),
          event_type = evt,
          timestamp = sh$timestamp[i],
          user = sh$user[i] %||% NA,
          search_name = sh$search_name[i] %||% NA,
          project = if (!is.na(od) && !is.null(proj_map[[od]])) proj_map[[od]] else NA,
          notes = sh$notes[i] %||% NA,
          backend = sh$backend[i] %||% NA,
          search_mode = sh$search_mode[i] %||% NA,
          parallel = sh$parallel[i] %||% NA,
          n_files = sh$n_files[i] %||% NA,
          fasta_files = sh$fasta_files[i] %||% NA,
          fasta_seq_count = sh$fasta_seq_count[i] %||% NA,
          normalization = sh$normalization[i] %||% NA,
          enzyme = sh$enzyme[i] %||% NA,
          mass_acc_mode = sh$mass_acc_mode[i] %||% NA,
          mass_acc = sh$mass_acc[i] %||% NA,
          mass_acc_ms1 = sh$mass_acc_ms1[i] %||% NA,
          scan_window = sh$scan_window[i] %||% NA,
          mbr = sh$mbr[i] %||% NA,
          extra_cli_flags = sh$extra_cli_flags[i] %||% NA,
          output_dir = od,
          job_id = sh$job_id[i] %||% NA,
          status = sh$status[i] %||% NA,
          duration_min = sh$duration_min[i] %||% NA,
          speclib_cached = sh$speclib_cached[i] %||% NA,
          app_version = sh$app_version[i] %||% NA,
          source_type = "search"
        )
        Sys.sleep(0.001)  # ensure unique IDs
      }
    }, error = function(e) message("[DE-LIMP] Search history migration error: ", e$message))
  }

  # Migrate analysis history rows
  if (!is.null(ah_path)) {
    tryCatch({
      ah <- read.csv(ah_path, stringsAsFactors = FALSE)
      for (i in seq_len(nrow(ah))) {
        od <- ah$output_dir[i] %||% NA
        rows[[length(rows) + 1]] <- list(
          id = generate_activity_id(),
          event_type = "analysis_completed",
          timestamp = ah$timestamp[i],
          user = ah$user[i] %||% NA,
          search_name = ah$source_file[i] %||% NA,
          project = if (!is.na(od) && !is.null(proj_map[[od]])) proj_map[[od]] else NA,
          notes = ah$notes[i] %||% NA,
          output_dir = od,
          n_proteins = ah$n_proteins[i] %||% NA,
          n_samples = ah$n_samples[i] %||% NA,
          n_contrasts = ah$n_contrasts[i] %||% NA,
          n_de_proteins = ah$n_de_proteins[i] %||% NA,
          session_file = ah$session_file[i] %||% NA,
          fasta_files = ah$fasta_file[i] %||% NA,
          fasta_seq_count = ah$fasta_seq_count[i] %||% NA,
          app_version = ah$app_version[i] %||% NA,
          source_type = ah$source_type[i] %||% NA
        )
        Sys.sleep(0.001)
      }
    }, error = function(e) message("[DE-LIMP] Analysis history migration error: ", e$message))
  }

  if (length(rows) == 0) return(invisible(FALSE))

  # Write unified log — ensure every row has all headers
  row_dfs <- lapply(rows, function(r) {
    vals <- lapply(activity_log_headers, function(h) {
      v <- r[[h]]
      if (is.null(v)) NA else v
    })
    df <- as.data.frame(vals, stringsAsFactors = FALSE)
    names(df) <- activity_log_headers
    df
  })
  df <- do.call(rbind, row_dfs)
  df <- df[order(df$timestamp, decreasing = FALSE), ]

  tryCatch({
    write.csv(df, file = new_path, row.names = FALSE, quote = TRUE)
    message(sprintf("[DE-LIMP] Migrated %d entries to activity log: %s", nrow(df), new_path))

    # Rename old files to .bak
    for (f in c(sh_path, ah_path, pj_path)) {
      if (!is.null(f) && file.exists(f)) {
        tryCatch(file.rename(f, paste0(f, ".bak")), error = function(e) NULL)
      }
    }
  }, error = function(e) message("[DE-LIMP] Migration write failed: ", e$message))

  invisible(TRUE)
}

#' Backfill activity log from existing job queue RDS
backfill_activity_log <- function(jobs, path = activity_log_path(),
                                   app_version = "unknown") {
  if (length(jobs) == 0) return(invisible(NULL))

  existing <- activity_log_read(path)
  existing_ods <- if (nrow(existing) > 0) existing$output_dir else character(0)

  n_added <- 0
  for (j in jobs) {
    if (j$output_dir %in% existing_ods) next

    ss <- j$search_settings
    sp <- if (!is.null(ss)) ss$search_params else list()

    dur <- if (!is.null(j$submitted_at) && !is.null(j$completed_at)) {
      round(as.numeric(difftime(j$completed_at, j$submitted_at, units = "mins")), 1)
    } else NA

    record_activity(list(
      event_type = "search_submitted",
      timestamp = if (!is.null(j$submitted_at)) format(j$submitted_at, "%Y-%m-%d %H:%M:%S") else NA,
      user = Sys.info()[["user"]],
      search_name = j$name %||% NA,
      backend = j$backend %||% "hpc",
      search_mode = ss$search_mode %||% j$search_mode %||% NA,
      parallel = isTRUE(j$parallel),
      n_files = j$n_files %||% ss$n_raw_files %||% NA,
      fasta_files = if (!is.null(ss$fasta_files)) paste(basename(ss$fasta_files), collapse = ", ") else NA,
      fasta_seq_count = ss$fasta_seq_count %||% NA,
      normalization = ss$normalization %||% NA,
      enzyme = sp$enzyme %||% NA,
      mass_acc_mode = sp$mass_acc_mode %||% NA,
      mass_acc = sp$mass_acc %||% NA,
      mass_acc_ms1 = sp$mass_acc_ms1 %||% NA,
      scan_window = sp$scan_window %||% NA,
      mbr = isTRUE(sp$mbr),
      extra_cli_flags = sp$extra_cli_flags %||% NA,
      output_dir = j$output_dir,
      job_id = j$job_id %||% NA,
      status = j$status %||% "unknown",
      duration_min = dur,
      speclib_cached = isTRUE(j$speclib_cached),
      app_version = app_version,
      notes = "Backfilled from job queue",
      source_type = "search"
    ), path = path)
    n_added <- n_added + 1
  }

  if (n_added > 0) message(sprintf("[DE-LIMP] Backfilled %d entries into activity log", n_added))
  invisible(n_added)
}

# Legacy aliases — keep temporarily for any straggling references
analysis_history_path <- activity_log_path
search_history_path <- activity_log_path

#' Count DE proteins across all contrasts
#'
#' @param fit limma fit object with contrasts
#' @param alpha Numeric: FDR threshold (default 0.05)
#' @return Integer: total DE proteins across all contrasts
count_de_proteins <- function(fit, alpha = 0.05) {
  tryCatch({
    n <- 0
    for (coef in colnames(fit$contrasts)) {
      tt <- limma::topTable(fit, coef = coef, number = Inf, sort.by = "none")
      n <- n + sum(tt$adj.P.Val < alpha, na.rm = TRUE)
    }
    n
  }, error = function(e) NA)
}

# Legacy function stubs for core facility module compatibility
record_search <- function(entry, path = activity_log_path()) {
  entry$event_type <- "search_submitted"
  entry$source_type <- "search"
  record_activity(entry, path)
}

update_search_status <- function(output_dir, status, completed_at = NA,
                                  duration_min = NA,
                                  path = activity_log_path()) {
  updates <- list(status = status)
  if (!is.na(completed_at)) updates$event_type <- "search_completed"
  if (!is.na(duration_min)) updates$duration_min <- duration_min
  update_activity(output_dir, updates, event_type_filter = "search_submitted", path = path)
}

record_analysis <- function(entry, path = activity_log_path()) {
  entry$event_type <- "analysis_completed"
  record_activity(entry, path)
}

search_history_read <- function(path = activity_log_path()) activity_log_read(path)
analysis_history_read <- function(path = activity_log_path()) activity_log_read(path)
backfill_search_history <- function(jobs, path = activity_log_path(), app_version = "unknown") {
  backfill_activity_log(jobs, path, app_version)
}

# =============================================================================
# Cluster Usage History — persistent resource monitoring for grant reporting
# =============================================================================

cluster_usage_headers <- c(
  "timestamp", "user", "account", "partition",
  "group_limit", "group_used", "group_available",
  "user_limit", "user_used", "user_available",
  "partition_idle", "partition_total",
  "pending_count", "avg_wait_min", "max_wait_min",
  "auto_selected"
)

cluster_usage_history_path <- function() {
  file.path(Sys.getenv("HOME"), ".delimp_cluster_usage_history.csv")
}

#' Append a cluster resource snapshot to the usage history CSV
#'
#' @param lab_res Result from check_cluster_resources() for genome-center-grp
#' @param pub_res Result from check_cluster_resources() for publicgrp
#' @param auto_partition Result from select_best_partition()
record_cluster_snapshot <- function(lab_res, pub_res, auto_partition,
                                     path = cluster_usage_history_path()) {
  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  usr <- Sys.info()[["user"]]
  auto_sel <- if (!is.null(auto_partition))
    paste0(auto_partition$account, "/", auto_partition$partition) else NA

  rows <- list()

  # Record genome-center-grp row

  if (!is.null(lab_res) && isTRUE(lab_res$success)) {
    rows[[1]] <- data.frame(
      timestamp = ts, user = usr,
      account = "genome-center-grp", partition = "high",
      group_limit = lab_res$group_limit %||% NA,
      group_used = lab_res$group_used %||% NA,
      group_available = lab_res$group_available %||% NA,
      user_limit = lab_res$user_limit %||% NA,
      user_used = lab_res$user_used %||% NA,
      user_available = lab_res$user_available %||% NA,
      partition_idle = lab_res$partition_idle %||% NA,
      partition_total = lab_res$partition_total %||% NA,
      pending_count = lab_res$pending_count %||% NA,
      avg_wait_min = round(lab_res$avg_wait_min %||% NA, 1),
      max_wait_min = round(lab_res$max_wait_min %||% NA, 1),
      auto_selected = auto_sel,
      stringsAsFactors = FALSE
    )
  }

  # Record publicgrp row
  if (!is.null(pub_res) && isTRUE(pub_res$success)) {
    rows[[length(rows) + 1]] <- data.frame(
      timestamp = ts, user = usr,
      account = "publicgrp", partition = "low",
      group_limit = pub_res$group_limit %||% NA,
      group_used = pub_res$group_used %||% NA,
      group_available = pub_res$group_available %||% NA,
      user_limit = pub_res$user_limit %||% NA,
      user_used = pub_res$user_used %||% NA,
      user_available = pub_res$user_available %||% NA,
      partition_idle = pub_res$partition_idle %||% NA,
      partition_total = pub_res$partition_total %||% NA,
      pending_count = pub_res$pending_count %||% NA,
      avg_wait_min = round(pub_res$avg_wait_min %||% NA, 1),
      max_wait_min = round(pub_res$max_wait_min %||% NA, 1),
      auto_selected = auto_sel,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) return(invisible(NULL))
  combined <- do.call(rbind, rows)

  lock_path <- paste0(path, ".lock")
  lock <- filelock::lock(lock_path, timeout = 5000)
  on.exit(filelock::unlock(lock), add = TRUE)

  needs_header <- !file.exists(path)
  tryCatch({
    suppressWarnings(
      write.table(combined, file = path, append = TRUE, sep = ",",
        row.names = FALSE, col.names = needs_header, quote = TRUE)
    )
  }, error = function(e) message("[DE-LIMP] Failed to write cluster usage: ", e$message))
}

#' Record per-job wait time for grant reporting
#' Appends to ~/.delimp_job_wait_log.csv
#' @param job Job list entry with job_id, submitted_at, wait_min, n_files, etc.
record_job_wait <- function(job) {
  path <- file.path(Sys.getenv("HOME"), ".delimp_job_wait_log.csv")
  # For parallel jobs, cpus/mem_gb/partition are in search_settings$slurm
  slurm_ss <- job$search_settings$slurm %||% list()
  par_ss <- job$search_settings$parallel %||% list()
  cpus <- job$cpus %||% slurm_ss$cpus %||% NA_integer_
  mem_gb <- job$mem_gb %||% slurm_ss$mem_gb %||% NA_integer_
  partition <- job$partition %||% slurm_ss$partition %||% NA_character_
  # For parallel: also record per-file resources
  cpus_per_file <- par_ss$cpus_per_file %||% NA_integer_
  row <- data.frame(
    timestamp    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    job_id       = job$job_id %||% NA_character_,
    name         = job$name %||% "unnamed",
    backend      = job$backend %||% "hpc",
    partition    = partition,
    n_files      = job$n_files %||% 0L,
    cpus         = cpus,
    mem_gb       = mem_gb,
    parallel     = isTRUE(job$parallel),
    cpus_per_file = cpus_per_file,
    priority     = job$priority %||% NA_integer_,
    submitted_at = if (!is.null(job$submitted_at)) format(job$submitted_at, "%Y-%m-%dT%H:%M:%S") else NA,
    started_at   = if (!is.null(job$started_at)) format(job$started_at, "%Y-%m-%dT%H:%M:%S") else NA,
    wait_min     = job$wait_min %||% NA_real_,
    est_start    = job$est_start %||% NA_character_,
    stringsAsFactors = FALSE
  )
  needs_header <- !file.exists(path)
  # If existing file has fewer columns (schema upgrade), rewrite with new header
  if (!needs_header) {
    existing <- tryCatch(read.csv(path, nrows = 1, stringsAsFactors = FALSE),
      error = function(e) NULL)
    if (!is.null(existing) && ncol(existing) < ncol(row)) {
      # Re-read all rows, add missing columns, rewrite
      all_rows <- tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
      if (!is.null(all_rows)) {
        for (col in setdiff(names(row), names(all_rows))) all_rows[[col]] <- NA
        all_rows <- rbind(all_rows[, names(row)], row)
        write.csv(all_rows, file = path, row.names = FALSE, quote = TRUE)
        return(invisible(NULL))
      }
    }
  }
  tryCatch({
    write.table(row, file = path, append = TRUE, sep = ",",
      row.names = FALSE, col.names = needs_header, quote = TRUE)
  }, error = function(e) message("[DE-LIMP] Failed to write job wait log: ", e$message))
}

#' Read job wait log for grant reporting
#' @return Data frame with all recorded job wait times
read_job_wait_log <- function() {
  path <- file.path(Sys.getenv("HOME"), ".delimp_job_wait_log.csv")
  if (!file.exists(path)) return(data.frame())
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) data.frame())
}

#' Read cluster usage history, optionally filtered by time
#'
#' @param since POSIXct timestamp — only return rows newer than this
#' @param account Filter to specific account (NULL for all)
cluster_usage_history_read <- function(path = cluster_usage_history_path(),
                                        since = NULL, account = NULL) {
  if (!file.exists(path)) return(data.frame())
  tryCatch({
    # Detect and fix header/column mismatch (schema evolved from 13 to 16 cols)
    lines <- readLines(path, n = 2)
    if (length(lines) >= 2) {
      header_cols <- length(strsplit(lines[1], ",")[[1]])
      data_cols <- length(strsplit(lines[2], ",")[[1]])
      if (header_cols < data_cols) {
        message(sprintf("[DE-LIMP] Fixing cluster usage CSV header: %d -> %d columns", header_cols, data_cols))
        all_lines <- readLines(path)
        all_lines[1] <- paste0('"', paste(cluster_usage_headers, collapse = '","'), '"')
        writeLines(all_lines, path)
      }
    }
    df <- read.csv(path, stringsAsFactors = FALSE)
    if (nrow(df) == 0) return(df)
    df$timestamp <- as.POSIXct(df$timestamp, format = "%Y-%m-%dT%H:%M:%S")
    if (!is.null(since)) df <- df[!is.na(df$timestamp) & df$timestamp >= since, ]
    if (!is.null(account)) df <- df[df$account == account, ]
    df
  }, error = function(e) data.frame())
}

#' Summarize cluster usage for grant reporting (hourly aggregation)
#'
#' @param df Data frame from cluster_usage_history_read()
#' @return Data frame with hourly summary statistics
# =============================================================================
# Per-user resource tracking — CPU + memory for each user in the group
# =============================================================================

per_user_usage_headers <- c(
  "timestamp", "account", "partition", "username",
  "cpus_running", "mem_gb_running", "n_jobs_running",
  "cpus_pending", "n_jobs_pending"
)

per_user_usage_path <- function() {
  file.path(Sys.getenv("HOME"), ".delimp_per_user_usage.csv")
}

#' Get path for lab members config JSON
lab_members_path <- function() {
  file.path(Sys.getenv("HOME"), ".delimp_lab_members.json")
}

#' Read lab member HPC usernames
#' @param ssh_user Current SSH username (always included in result)
get_lab_members <- function(ssh_user = NULL) {
  path <- lab_members_path()
  members <- character(0)
  if (file.exists(path)) {
    tryCatch({
      cfg <- jsonlite::fromJSON(path)
      members <- as.character(cfg$members)
    }, error = function(e) NULL)
  }
  # Always include current SSH user
  if (!is.null(ssh_user) && nzchar(ssh_user)) {
    members <- unique(c(members, ssh_user))
  }
  if (length(members) == 0) members <- ssh_user %||% Sys.info()[["user"]]
  members
}

#' Query per-user CPU and memory usage for lab members
#'
#' @param ssh_config SSH config list or NULL for local
#' @param account SLURM account name
#' @param partition SLURM partition name
#' @param sbatch_path Path to sbatch (used to derive squeue path)
#' @param members Character vector of usernames to track
#' @return Data frame with one row per user
check_per_user_resources <- function(ssh_config, account, partition, sbatch_path = NULL,
                                      members = get_lab_members()) {
  slurm_cmd <- function(cmd) {
    if (!is.null(sbatch_path) && nzchar(sbatch_path)) {
      file.path(dirname(sbatch_path), cmd)
    } else cmd
  }

  run_cmd <- function(command) {
    if (!is.null(ssh_config)) {
      ssh_exec(ssh_config, command, login_shell = is.null(sbatch_path), timeout = 15)
    } else if (slurm_proxy_available()) {
      slurm_proxy_exec(command, timeout = 15)
    } else {
      parts <- strsplit(command, " ")[[1]]
      stdout <- tryCatch(
        system2(parts[1], args = parts[-1], stdout = TRUE, stderr = TRUE),
        error = function(e) structure(e$message, status = 1L))
      list(status = attr(stdout, "status") %||% 0L, stdout = stdout)
    }
  }

  ts <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  rows <- list()

  # Build user filter: only query specified lab members
  user_filter <- paste(members, collapse = ",")

  # Query running jobs: user, CPUs, memory per job
  tryCatch({
    # %u=user, %C=CPUs, %m=min_memory
    squeue_cmd <- sprintf(
      '%s -u %s -A %s -p %s -t RUNNING -o "%%u|%%C|%%m" --noheader',
      slurm_cmd("squeue"), user_filter, account, partition
    )
    res <- run_cmd(squeue_cmd)
    if (res$status == 0 && length(res$stdout) > 0) {
      lines <- trimws(res$stdout)
      lines <- lines[nzchar(lines)]
      if (length(lines) > 0) {
        # Parse each line: user|cpus|memory
        parsed <- lapply(lines, function(line) {
          parts <- strsplit(line, "\\|")[[1]]
          if (length(parts) < 3) return(NULL)
          mem_str <- trimws(parts[3])
          # Parse memory: "64G", "65536M", "64000", etc.
          mem_gb <- tryCatch({
            if (grepl("G$", mem_str, ignore.case = TRUE)) {
              as.numeric(sub("[Gg]$", "", mem_str))
            } else if (grepl("M$", mem_str, ignore.case = TRUE)) {
              as.numeric(sub("[Mm]$", "", mem_str)) / 1024
            } else if (grepl("T$", mem_str, ignore.case = TRUE)) {
              as.numeric(sub("[Tt]$", "", mem_str)) * 1024
            } else {
              as.numeric(mem_str) / 1024  # assume MB
            }
          }, error = function(e) 0)
          list(user = trimws(parts[1]), cpus = as.integer(parts[2]), mem_gb = mem_gb)
        })
        parsed <- Filter(Negate(is.null), parsed)

        # Aggregate by user
        user_data <- list()
        for (p in parsed) {
          u <- p$user
          if (is.null(user_data[[u]])) user_data[[u]] <- list(cpus = 0L, mem_gb = 0, n = 0L)
          user_data[[u]]$cpus <- user_data[[u]]$cpus + p$cpus
          user_data[[u]]$mem_gb <- user_data[[u]]$mem_gb + p$mem_gb
          user_data[[u]]$n <- user_data[[u]]$n + 1L
        }

        for (u in names(user_data)) {
          rows[[length(rows) + 1]] <- data.frame(
            timestamp = ts, account = account, partition = partition, username = u,
            cpus_running = user_data[[u]]$cpus,
            mem_gb_running = round(user_data[[u]]$mem_gb, 1),
            n_jobs_running = user_data[[u]]$n,
            cpus_pending = 0L, n_jobs_pending = 0L,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }, error = function(e) NULL)

  # Query pending jobs per user
  # Include job ID (%i) to detect array jobs and compute real CPU demand
  tryCatch({
    squeue_cmd <- sprintf(
      '%s -u %s -A %s -p %s -t PENDING -o "%%i|%%u|%%C" --noheader',
      slurm_cmd("squeue"), user_filter, account, partition
    )
    res <- run_cmd(squeue_cmd)
    if (res$status == 0 && length(res$stdout) > 0) {
      lines <- trimws(res$stdout)
      lines <- lines[nzchar(lines)]
      if (length(lines) > 0) {
        pending <- list()
        for (line in lines) {
          parts <- strsplit(line, "\\|")[[1]]
          if (length(parts) < 3) next
          job_id <- trimws(parts[1])
          u <- trimws(parts[2])
          cpus_per_task <- as.integer(parts[3])
          # For array jobs like "12345_[0-231%200]", compute max concurrent CPUs
          # n_tasks = range size, max_simultaneous = throttle (%N)
          array_match <- regexec("_\\[(\\d+)-(\\d+)(%\\d+)?\\]$", job_id)
          if (array_match[[1]][1] != -1) {
            parts_m <- regmatches(job_id, array_match)[[1]]
            lo <- as.integer(parts_m[2])
            hi <- as.integer(parts_m[3])
            n_tasks <- hi - lo + 1L
            max_simul <- if (nzchar(parts_m[4])) as.integer(sub("^%", "", parts_m[4])) else n_tasks
            cpus <- cpus_per_task * n_tasks
            n_jobs <- as.integer(n_tasks)
          } else {
            cpus <- cpus_per_task
            n_jobs <- 1L
          }
          if (is.null(pending[[u]])) pending[[u]] <- list(cpus = 0L, n = 0L)
          pending[[u]]$cpus <- pending[[u]]$cpus + cpus
          pending[[u]]$n <- pending[[u]]$n + n_jobs
        }
        # Merge into existing rows or add new
        existing_users <- vapply(rows, function(r) r$username, character(1))
        for (u in names(pending)) {
          idx <- which(existing_users == u)
          if (length(idx) > 0) {
            rows[[idx[1]]]$cpus_pending <- pending[[u]]$cpus
            rows[[idx[1]]]$n_jobs_pending <- pending[[u]]$n
          } else {
            rows[[length(rows) + 1]] <- data.frame(
              timestamp = ts, account = account, partition = partition, username = u,
              cpus_running = 0L, mem_gb_running = 0, n_jobs_running = 0L,
              cpus_pending = pending[[u]]$cpus, n_jobs_pending = pending[[u]]$n,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }, error = function(e) NULL)

  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

#' Record per-user resource snapshot to CSV
record_per_user_snapshot <- function(user_df, path = per_user_usage_path()) {
  if (is.null(user_df) || nrow(user_df) == 0) return(invisible(NULL))

  lock_path <- paste0(path, ".lock")
  lock <- filelock::lock(lock_path, timeout = 5000)
  on.exit(filelock::unlock(lock), add = TRUE)

  needs_header <- !file.exists(path)
  tryCatch({
    suppressWarnings(
      write.table(user_df, file = path, append = TRUE, sep = ",",
        row.names = FALSE, col.names = needs_header, quote = TRUE))
  }, error = function(e) message("[DE-LIMP] Failed to write per-user usage: ", e$message))
}

#' Read per-user usage history
per_user_usage_read <- function(path = per_user_usage_path(), since = NULL, account = NULL) {
  if (!file.exists(path)) return(data.frame())
  tryCatch({
    df <- read.csv(path, stringsAsFactors = FALSE)
    if (nrow(df) == 0) return(df)
    df$timestamp <- as.POSIXct(df$timestamp, format = "%Y-%m-%dT%H:%M:%S")
    if (!is.null(since)) df <- df[!is.na(df$timestamp) & df$timestamp >= since, ]
    if (!is.null(account)) df <- df[df$account == account, ]
    df
  }, error = function(e) data.frame())
}

cluster_usage_grant_summary <- function(df) {
  if (nrow(df) == 0) return(data.frame())

  # Filter to genome-center-grp only
  df <- df[df$account == "genome-center-grp", ]
  if (nrow(df) == 0) return(data.frame())

  df$date <- as.Date(df$timestamp)
  df$hour <- as.integer(format(df$timestamp, "%H"))

  # Aggregate by date + hour
  result <- do.call(rbind, lapply(split(df, paste(df$date, df$hour)), function(chunk) {
    data.frame(
      date = chunk$date[1],
      hour = chunk$hour[1],
      account = "genome-center-grp",
      n_snapshots = nrow(chunk),
      avg_group_used = round(mean(chunk$group_used, na.rm = TRUE), 1),
      max_group_used = max(chunk$group_used, na.rm = TRUE),
      avg_user_used = round(mean(chunk$user_used, na.rm = TRUE), 1),
      max_user_used = max(chunk$user_used, na.rm = TRUE),
      group_limit = chunk$group_limit[1],
      user_limit = chunk$user_limit[1],
      pct_group_utilization = round(mean(chunk$group_used, na.rm = TRUE) /
        max(chunk$group_limit[1], 1, na.rm = TRUE) * 100, 1),
      pct_user_at_capacity = round(
        sum(!is.na(chunk$user_available) & chunk$user_available < 64) /
        max(nrow(chunk), 1) * 100, 1),
      avg_pending_count = if ("pending_count" %in% names(chunk))
        round(mean(chunk$pending_count, na.rm = TRUE), 1) else NA_real_,
      avg_wait_min = if ("avg_wait_min" %in% names(chunk))
        round(mean(chunk$avg_wait_min, na.rm = TRUE), 1) else NA_real_,
      max_wait_min = if ("max_wait_min" %in% names(chunk))
        round(max(chunk$max_wait_min, na.rm = TRUE), 1) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result[order(result$date, result$hour), ]
}
