# helpers_slims.R — Ingestion helpers for the Proteogenomics builder.
# No Shiny reactivity. Pure functions.
#
# Two ingestion modes:
#   Mode A — SLIMS URL (UC Davis DNA Tech Core delivery)
#   Mode B — SRA/ENA accession (re-analysis / external collaborators)
#
# Both modes share the same downstream pipeline. This file handles only the
# pre-pipeline parts: URL scanning, accession metadata verification, download
# launchers (login-node `nohup`, not sbatch — compute nodes may lack outbound HTTP).
#
# Also hosts the reference-registry loaders since they're conceptually about
# "where do data sources live."

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# =============================================================================
# Registry I/O — read-only loaders. Writers live in helpers_proteog_assembly.R.
# =============================================================================

.reference_registry_path <- function() {
  Sys.getenv(
    "DELIMP_REFERENCE_REGISTRY",
    unset = "/quobyte/proteomics-grp/de-limp/references/registry.json"
  )
}

#' Load the reference-genomes registry as a named list
#'
#' Empty/missing registry returns `list()`. Never throws.
#'
#' When `ssh_config` is non-NULL and the user is connected to Hive,
#' fetches the registry via `ssh_exec("cat <path>")` instead of a local
#' filesystem read. This is how DE-LIMP-on-Mac with SSH-to-Hive picks
#' up the registry on /quobyte/.
#'
#' @param ssh_config NULL (local read) OR ssh_config list with $host/$user/$key
#' @return named list keyed by reference_key (e.g. "mm39_GRCm39")
load_reference_registry <- function(ssh_config = NULL) {
  path <- .reference_registry_path()

  raw <- if (!is.null(ssh_config) && exists("ssh_exec")) {
    res <- tryCatch(
      ssh_exec(ssh_config, sprintf("cat %s", shQuote(path)),
               login_shell = FALSE, timeout = 15),
      error = function(e) NULL
    )
    if (is.null(res) || !identical(res$status, 0L) ||
        length(res$stdout) == 0) {
      return(list())
    }
    paste(res$stdout, collapse = "\n")
  } else {
    if (!file.exists(path)) return(list())
    tryCatch(paste(readLines(path, warn = FALSE), collapse = "\n"),
             error = function(e) NULL)
  }

  if (is.null(raw) || !nzchar(raw)) return(list())
  parsed <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE),
                     error = function(e) NULL)
  if (is.null(parsed)) {
    warning("load_reference_registry(): could not parse registry at ", path)
    return(list())
  }
  parsed
}

# =============================================================================
# Mode A — SLIMS URL scanning
# =============================================================================

#' Validate that a URL is a SLIMS data delivery URL
#'
#' SLIMS URLs look like:
#'   http://slimsdata.genomecenter.ucdavis.edu/Data/<random_id>/Unaligned/
#'
#' @param slims_url character
#' @return logical
is_slims_url <- function(slims_url) {
  isTRUE(nzchar(slims_url) &&
         grepl("^https?://slimsdata\\.genomecenter\\.ucdavis\\.edu/Data/[a-z0-9]+/",
               slims_url, ignore.case = TRUE))
}

#' Scan a SLIMS URL for sample files
#'
#' Returns a structured list with $success (logical), $error (if !success),
#' and on success: $url, $n_samples, $sample_names, $files, $is_paired, $has_md5.
#'
#' The HTTP fetch runs synchronously and respects a short timeout — this is
#' called from the UI to validate the URL before the user clicks Build.
#'
#' @param slims_url character
#' @param timeout_sec integer — HTTP timeout (default 20s)
#' @return list (see above)
scan_slims_url <- function(slims_url, timeout_sec = 20L) {
  if (!is_slims_url(slims_url)) {
    return(list(
      success = FALSE,
      error   = "URL does not match SLIMS format. Expected http(s)://slimsdata.genomecenter.ucdavis.edu/Data/<id>/..."
    ))
  }

  index_html <- tryCatch({
    if (requireNamespace("httr2", quietly = TRUE)) {
      resp <- httr2::request(slims_url) |>
        httr2::req_timeout(timeout_sec) |>
        httr2::req_perform()
      httr2::resp_body_string(resp)
    } else if (requireNamespace("httr", quietly = TRUE)) {
      httr::content(
        httr::GET(slims_url, httr::timeout(timeout_sec)),
        as = "text", encoding = "UTF-8"
      )
    } else {
      stop("Neither httr2 nor httr is available; cannot fetch URL")
    }
  }, error = function(e) NULL)

  if (is.null(index_html) || !nzchar(index_html)) {
    return(list(
      success = FALSE,
      error   = paste0(
        "Could not reach SLIMS from this host. Possible causes: ",
        "(1) URL expired (SLIMS retains data ~1 month); ",
        "(2) network/firewall block from this node; ",
        "(3) URL typo. Verify the link in your SLIMS email and retry."
      )
    ))
  }

  file_pattern <- "[A-Za-z0-9_.-]+\\.fastq\\.gz"
  files <- unique(regmatches(index_html, gregexpr(file_pattern, index_html))[[1]])

  if (length(files) == 0) {
    return(list(
      success = FALSE,
      error   = paste0(
        "No .fastq.gz files found at this URL. ",
        "If this URL was issued more than ~1 month ago, the data may have been purged. ",
        "Contact DNA Technologies Core to confirm."
      )
    ))
  }

  r1_files  <- grep("_R1[._]", files, value = TRUE)
  r2_files  <- grep("_R2[._]", files, value = TRUE)
  is_paired <- length(r1_files) > 0 && length(r2_files) > 0

  sample_names <- if (is_paired) {
    gsub("_R1[._].*", "", r1_files)
  } else {
    gsub("\\.fastq\\.gz$", "", files)
  }
  sample_names <- unique(sample_names)

  list(
    success      = TRUE,
    url          = slims_url,
    n_samples    = length(sample_names),
    sample_names = sample_names,
    files        = files,
    is_paired    = is_paired,
    has_md5      = grepl("checksums\\.md5", index_html, ignore.case = TRUE)
  )
}

# =============================================================================
# Mode B — SRA/ENA accession metadata verification
# =============================================================================

# Library strategies that are NOT suitable for proteogenomics novel-ORF discovery.
# Validation surfaced this gate when an "RNA-Seq" claim turned out to be Tag-Seq.
.UNSUITABLE_LIBRARY_STRATEGIES <- c("Tag-Seq", "miRNA-Seq", "OTHER", "small-RNA")

#' Fetch and parse ENA metadata for an SRA/ENA accession
#'
#' Critical: must run BEFORE download, not after. Validation surfaced 30 min
#' of wasted compute when an SRR claimed to be K562 human turned out to be
#' Mus musculus. A 2-second metadata call would have caught it.
#'
#' Uses ENA's PORTAL API (returns plain TSV) instead of the BROWSER API (XML).
#' The TSV endpoint is more robust across R/xml2 versions — no XML parsing,
#' no encoding quirks — and lets us request exactly the fields we need.
#'
#' @param accession character — e.g., "SRR1303776"
#' @param timeout_sec integer
#' @return list with $success, $error (on failure), and on success:
#'   accession, scientific_name, library_strategy, library_source,
#'   library_selection, instrument, layout, suitable (logical),
#'   suitable_warning (NULL or character — non-fatal warning)
verify_sra_accession <- function(accession, timeout_sec = 15L) {
  if (!nzchar(accession) || !grepl("^[A-Z]{3}[0-9]+$", accession)) {
    return(list(
      success = FALSE,
      error   = paste0("Invalid accession format: ", accession,
                       " (expected e.g. SRR1303776 / ERR12345)")
    ))
  }

  fields <- "run_accession,scientific_name,library_strategy,library_source,library_selection,instrument_model,library_layout"
  url <- sprintf(
    "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=%s&result=read_run&fields=%s&format=tsv",
    accession, fields
  )

  raw <- tryCatch({
    if (requireNamespace("httr2", quietly = TRUE)) {
      resp <- httr2::request(url) |>
        httr2::req_timeout(timeout_sec) |>
        httr2::req_perform()
      httr2::resp_body_string(resp)
    } else if (requireNamespace("httr", quietly = TRUE)) {
      r <- httr::GET(url, httr::timeout(timeout_sec))
      if (httr::status_code(r) != 200) stop("HTTP ", httr::status_code(r))
      httr::content(r, as = "text", encoding = "UTF-8")
    } else {
      stop("Neither httr2 nor httr available")
    }
  }, error = function(e) conditionMessage(e))

  if (is.null(raw) || !is.character(raw) || !nzchar(raw)) {
    return(list(
      success = FALSE,
      error   = paste0("Could not reach ENA metadata API for ", accession,
                       ". Check network access from this host.")
    ))
  }

  lines <- strsplit(raw, "\n", fixed = TRUE)[[1]]
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2) {
    return(list(
      success = FALSE,
      error   = paste0("ENA returned no record for ", accession,
                       " — accession may not exist, or may not have public reads yet.")
    ))
  }

  header <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  fields_row <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
  # Pad short rows so column indexing doesn't blow up on a truncated response
  if (length(fields_row) < length(header)) {
    fields_row <- c(fields_row, rep(NA_character_, length(header) - length(fields_row)))
  }
  rec <- setNames(as.list(fields_row), header)

  .get <- function(k) {
    v <- rec[[k]]
    if (is.null(v) || is.na(v) || !nzchar(v)) NA_character_ else v
  }

  scientific_name   <- .get("scientific_name")
  library_strategy  <- .get("library_strategy")
  library_source    <- .get("library_source")
  library_selection <- .get("library_selection")
  instrument        <- .get("instrument_model")
  layout_raw        <- .get("library_layout")
  layout <- if (is.na(layout_raw)) "unknown"
            else if (toupper(layout_raw) == "PAIRED") "paired"
            else if (toupper(layout_raw) == "SINGLE") "single"
            else tolower(layout_raw)

  is_unsuitable <- !is.na(library_strategy) &&
    library_strategy %in% .UNSUITABLE_LIBRARY_STRATEGIES

  list(
    success           = TRUE,
    accession         = accession,
    scientific_name   = scientific_name,
    library_strategy  = library_strategy,
    library_source    = library_source,
    library_selection = library_selection,
    instrument        = instrument,
    layout            = layout,
    suitable          = !is_unsuitable,
    suitable_warning  = if (is_unsuitable) {
      sprintf(
        "Library strategy '%s' is not suitable for proteogenomics novel-ORF discovery (does not cover full transcripts). Use a different ingestion workflow.",
        library_strategy
      )
    } else NULL
  )
}

# =============================================================================
# Download launchers — login-node nohup, NOT sbatch
# =============================================================================
# Compute nodes may lack outbound HTTP. Spec §4 explicitly requires the
# download step to run on the login node so the SLIMS HTTP server and the
# ENA FTP/HTTP servers are reachable.

#' Sanitize a string for use as a filesystem path component
#'
#' Allows [A-Za-z0-9_.-] only; replaces everything else with `_`.
#' Used for project_name validation throughout the proteogenomics builder.
#'
#' @param name character
#' @return character — sanitized name; throws if input is empty
sanitize_project_name <- function(name) {
  if (!nzchar(name)) stop("sanitize_project_name(): name must be non-empty")
  gsub("[^A-Za-z0-9_.-]", "_", name)
}

#' Launch a SLIMS download as a background process on the local host
#'
#' Mirrors the SLIMS tree under `<rnaseq_root>/<project_name>/`, then verifies
#' md5 checksums. Writes a status file the UI can poll. Runs via `nohup` so
#' the Shiny session can return immediately.
#'
#' When `ssh_config` is non-NULL, the project dir is created on Hive via
#' ssh_exec, the status JSON is written there via scp_upload, and the
#' wget+md5 script runs on the Hive login node via ssh_exec with nohup +
#' </dev/null detachment (which lets the SSH connection return immediately
#' while the background bash keeps running).
#'
#' @param slims_url    character — URL from scan_slims_url()
#' @param project_name character — sanitized via sanitize_project_name()
#' @param rnaseq_root  character — base dir (default: spec'd path)
#' @param ssh_config   NULL (local) or ssh_config list (run on Hive via SSH)
#' @return list with $project_dir, $status_file, $download_log
launch_slims_download <- function(slims_url,
                                  project_name,
                                  rnaseq_root = "/quobyte/proteomics-grp/de-limp/rnaseq",
                                  ssh_config = NULL) {
  if (!is_slims_url(slims_url)) {
    stop("launch_slims_download(): URL does not match SLIMS format")
  }
  project_name <- sanitize_project_name(project_name)
  project_dir  <- file.path(rnaseq_root, project_name)

  status_file  <- file.path(project_dir, "download_status.json")
  download_log <- file.path(project_dir, "download.log")
  md5_log      <- file.path(project_dir, "md5_verify.log")

  status_json <- jsonlite::toJSON(
    list(state = "running",
         started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
         mode = "slims", url = slims_url),
    auto_unbox = TRUE, pretty = TRUE
  )

  # Inner command — runs the wget mirror then md5 verification, then writes
  # the final status JSON. Quote everything carefully.
  inner_cmd <- sprintf(
    paste0(
      "cd %s && ",
      "wget -q -r -nH --cut-dirs=3 -nc -R 'index.html*' %s > %s 2>&1; rc=$?; ",
      "if [ $rc -eq 0 ]; then ",
      "  if [ -f checksums.md5 ]; then md5sum -c checksums.md5 > %s 2>&1; mc=$?; else mc=0; fi; ",
      "  if [ $mc -eq 0 ]; then state=complete; else state=md5_failed; fi; ",
      "else state=download_failed; fi; ",
      "printf '{\"state\":\"%%s\",\"finished_at\":\"%%s\"}\\n' \"$state\" \"$(date -Iseconds)\" > %s"
    ),
    shQuote(project_dir),
    shQuote(slims_url),
    shQuote(download_log),
    shQuote(md5_log),
    shQuote(status_file)
  )

  if (is.null(ssh_config)) {
    # ── Local: dir.create + jsonlite::write_json + system(nohup ...) ──────
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(as.character(status_json), status_file)
    outer <- sprintf("nohup bash -c %s > %s 2>&1 &",
                     shQuote(inner_cmd), shQuote(download_log))
    system(outer)
  } else {
    # ── SSH: mkdir on Hive, scp status.json, ssh_exec nohup bash ──────────
    res_mk <- ssh_exec(ssh_config,
                       sprintf("mkdir -p %s", shQuote(project_dir)),
                       login_shell = FALSE, timeout = 15)
    if (!identical(res_mk$status, 0L)) {
      stop("launch_slims_download (ssh): mkdir failed: ",
           paste(res_mk$stderr %||% character(), collapse = "; "))
    }
    tmp <- tempfile(pattern = "slims_status_")
    on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
    writeLines(as.character(status_json), tmp)
    scp_res <- scp_upload(ssh_config, tmp, status_file, timeout = 60)
    if (!identical(scp_res$status, 0L)) {
      stop("launch_slims_download (ssh): scp status.json failed: ",
           paste(scp_res$stderr %||% character(), collapse = "; "))
    }
    # `< /dev/null` is mandatory — without it the SSH stays attached
    # to the bash stdin and won't return until the download finishes.
    bg <- sprintf("nohup bash -c %s > %s 2>&1 < /dev/null & echo BG_PID:$!",
                  shQuote(inner_cmd), shQuote(download_log))
    res_bg <- ssh_exec(ssh_config, bg, login_shell = FALSE, timeout = 15)
    if (!identical(res_bg$status, 0L)) {
      stop("launch_slims_download (ssh): nohup launch failed: ",
           paste(res_bg$stderr %||% character(), collapse = "; "))
    }
  }

  list(
    project_dir  = project_dir,
    status_file  = status_file,
    download_log = download_log
  )
}

#' Launch an ENA accession download as a background process
#'
#' Builds ENA FTP URLs from accession IDs (ENA mirrors essentially all of SRA).
#' Supports an optional read-pair subsample for quick testing.
#'
#' Set +o pipefail inside the wget|zcat|head|gzip pipeline so the SIGPIPE that
#' `head` produces when it closes upstream doesn't fail the script (validation
#' bug — fixed in spec §4).
#'
#' When `ssh_config` is non-NULL, the project dir + shell script are
#' created on Hive via ssh_exec / scp_upload, and the download runs on the
#' Hive login node via ssh_exec with nohup detachment.
#'
#' @param accessions      character vector — e.g., c("SRR1303776", "SRR1303777")
#' @param project_name    character
#' @param subsample_reads integer or NULL — if set, stream-subsample N read pairs
#' @param rnaseq_root     character
#' @param ssh_config      NULL (local) or ssh_config list (remote via SSH)
#' @return list (same shape as launch_slims_download)
launch_ena_download <- function(accessions,
                                project_name,
                                subsample_reads = NULL,
                                rnaseq_root = "/quobyte/proteomics-grp/de-limp/rnaseq",
                                ssh_config = NULL) {
  if (length(accessions) == 0) stop("launch_ena_download(): no accessions provided")
  if (length(accessions) > 24) {
    stop("launch_ena_download(): too many accessions (limit 24 per build for sanity)")
  }
  bad <- accessions[!grepl("^[A-Z]{3}[0-9]+$", accessions)]
  if (length(bad) > 0) {
    stop("launch_ena_download(): invalid accession format: ",
         paste(bad, collapse = ", "))
  }

  project_name <- sanitize_project_name(project_name)
  project_dir  <- file.path(rnaseq_root, project_name)

  status_file  <- file.path(project_dir, "download_status.json")
  download_log <- file.path(project_dir, "download.log")

  status_json <- jsonlite::toJSON(
    list(state = "running",
         started_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
         mode = "ena", accessions = as.list(accessions),
         subsample_reads = subsample_reads %||% NA),
    auto_unbox = TRUE, pretty = TRUE
  )

  # Build per-accession URL pairs. ENA FTP URL convention:
  #   https://ftp.sra.ebi.ac.uk/vol1/fastq/<SRR130>/<00X>/<SRR1303776>/SRR1303776_{1,2}.fastq.gz
  #   <00X> = "00" + last digit of accession when accession is 7 digits; for >7 digits use modulo (see ENA docs).
  # Use the simplest form that matches the common case; users can override
  # by giving a SLIMS URL instead for non-standard cases.
  build_urls <- function(acc) {
    n <- nchar(sub("^[A-Z]{3}", "", acc))
    prefix3 <- substr(acc, 1, 6)  # e.g. "SRR130"
    bucket <- if (n <= 6) {
      ""  # no bucket subdir
    } else if (n == 7) {
      sprintf("00%s/", substr(acc, nchar(acc), nchar(acc)))
    } else if (n == 8) {
      sprintf("0%s/", substr(acc, nchar(acc) - 1, nchar(acc)))
    } else {
      sprintf("%s/", substr(acc, nchar(acc) - 2, nchar(acc)))
    }
    base <- sprintf("https://ftp.sra.ebi.ac.uk/vol1/fastq/%s/%s%s",
                    prefix3, bucket, acc)
    c(
      r1 = sprintf("%s/%s_1.fastq.gz", base, acc),
      r2 = sprintf("%s/%s_2.fastq.gz", base, acc)
    )
  }

  # Compose a single shell script that downloads (or stream-subsamples) every
  # accession sequentially. set +o pipefail around the head pipe to tolerate SIGPIPE.
  per_acc_lines <- character()
  for (acc in accessions) {
    urls <- build_urls(acc)
    if (is.null(subsample_reads)) {
      per_acc_lines <- c(per_acc_lines, sprintf(
        "wget -q -O %s/%s_R1.fastq.gz %s && wget -q -O %s/%s_R2.fastq.gz %s || { echo FAILED_%s; exit 1; }",
        shQuote(project_dir), acc, shQuote(urls["r1"]),
        shQuote(project_dir), acc, shQuote(urls["r2"]),
        acc
      ))
    } else {
      # Stream-subsample. Each read is 4 lines, so 4*N lines per FASTQ.
      n_lines <- as.integer(subsample_reads) * 4L
      per_acc_lines <- c(per_acc_lines, sprintf(
        "( set +o pipefail; curl -fsSL %s | zcat | head -n %d | gzip > %s/%s_R1.fastq.gz ) && ( set +o pipefail; curl -fsSL %s | zcat | head -n %d | gzip > %s/%s_R2.fastq.gz ) || { echo FAILED_%s; exit 1; }",
        shQuote(urls["r1"]), n_lines, shQuote(project_dir), acc,
        shQuote(urls["r2"]), n_lines, shQuote(project_dir), acc,
        acc
      ))
    }
  }

  finalize_line <- sprintf(
    "printf '{\"state\":\"%%s\",\"finished_at\":\"%%s\"}\\n' \"$1\" \"$(date -Iseconds)\" > %s",
    shQuote(status_file)
  )

  script <- paste(c(
    "#!/bin/bash",
    "set -uo pipefail",
    paste(per_acc_lines, collapse = " && \\\n"),
    "rc=$?",
    sprintf("state=$([ $rc -eq 0 ] && echo complete || echo download_failed)"),
    sprintf("printf '{\"state\":\"%%s\",\"finished_at\":\"%%s\"}\\n' \"$state\" \"$(date -Iseconds)\" > %s",
            shQuote(status_file))
  ), collapse = "\n")

  script_path <- file.path(project_dir, "ena_download.sh")

  if (is.null(ssh_config)) {
    # ── Local ────────────────────────────────────────────────────────────
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
    writeLines(as.character(status_json), status_file)
    writeLines(script, script_path)
    Sys.chmod(script_path, "755")
    outer <- sprintf("nohup bash %s > %s 2>&1 &",
                     shQuote(script_path), shQuote(download_log))
    system(outer)
  } else {
    # ── SSH: mkdir, scp status+script to Hive, nohup over ssh ────────────
    res_mk <- ssh_exec(ssh_config,
                       sprintf("mkdir -p %s", shQuote(project_dir)),
                       login_shell = FALSE, timeout = 15)
    if (!identical(res_mk$status, 0L)) {
      stop("launch_ena_download (ssh): mkdir failed: ",
           paste(res_mk$stderr %||% character(), collapse = "; "))
    }
    tmp_status <- tempfile(pattern = "ena_status_")
    tmp_script <- tempfile(pattern = "ena_script_", fileext = ".sh")
    on.exit({
      if (file.exists(tmp_status)) file.remove(tmp_status)
      if (file.exists(tmp_script)) file.remove(tmp_script)
    }, add = TRUE)
    writeLines(as.character(status_json), tmp_status)
    writeLines(script, tmp_script)
    up1 <- scp_upload(ssh_config, tmp_status, status_file, timeout = 60)
    if (!identical(up1$status, 0L)) {
      stop("launch_ena_download (ssh): scp status.json failed: ",
           paste(up1$stderr %||% character(), collapse = "; "))
    }
    up2 <- scp_upload(ssh_config, tmp_script, script_path, timeout = 60)
    if (!identical(up2$status, 0L)) {
      stop("launch_ena_download (ssh): scp script failed: ",
           paste(up2$stderr %||% character(), collapse = "; "))
    }
    res_chmod <- ssh_exec(ssh_config,
                          sprintf("chmod 755 %s", shQuote(script_path)),
                          login_shell = FALSE, timeout = 15)
    if (!identical(res_chmod$status, 0L)) {
      stop("launch_ena_download (ssh): chmod failed: ",
           paste(res_chmod$stderr %||% character(), collapse = "; "))
    }
    # `< /dev/null` is mandatory — without it the SSH connection waits
    # on the background bash's stdin and won't return.
    bg <- sprintf("nohup bash %s > %s 2>&1 < /dev/null & echo BG_PID:$!",
                  shQuote(script_path), shQuote(download_log))
    res_bg <- ssh_exec(ssh_config, bg, login_shell = FALSE, timeout = 15)
    if (!identical(res_bg$status, 0L)) {
      stop("launch_ena_download (ssh): nohup launch failed: ",
           paste(res_bg$stderr %||% character(), collapse = "; "))
    }
  }

  list(
    project_dir  = project_dir,
    status_file  = status_file,
    download_log = download_log,
    script       = script_path
  )
}

#' Poll the status JSON for a download
#'
#' When `ssh_config` is non-NULL, fetches the JSON from Hive via ssh_exec
#' (cat) — short file, no need to scp.
#'
#' @param project_dir character — from launch_*_download() return
#' @param ssh_config  NULL (local) or ssh_config list (remote via SSH)
#' @return list with $state ("running"|"complete"|"download_failed"|"md5_failed"|"missing"), $finished_at
poll_download_status <- function(project_dir, ssh_config = NULL) {
  status_file <- file.path(project_dir, "download_status.json")
  raw <- NULL
  if (is.null(ssh_config)) {
    if (!file.exists(status_file)) {
      return(list(state = "missing", finished_at = NA_character_))
    }
    raw <- tryCatch(jsonlite::read_json(status_file), error = function(e) NULL)
  } else {
    res <- ssh_exec(ssh_config,
                    sprintf("cat %s 2>/dev/null || echo MISSING_STATUS",
                            shQuote(status_file)),
                    login_shell = FALSE, timeout = 15)
    txt <- paste(res$stdout %||% character(), collapse = "\n")
    if (!identical(res$status, 0L) ||
        identical(trimws(txt), "MISSING_STATUS") ||
        !nzchar(trimws(txt))) {
      return(list(state = "missing", finished_at = NA_character_))
    }
    raw <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE),
                    error = function(e) NULL)
  }
  if (is.null(raw)) {
    return(list(state = "missing", finished_at = NA_character_))
  }
  list(
    state       = raw$state %||% "unknown",
    started_at  = raw$started_at %||% NA_character_,
    finished_at = raw$finished_at %||% NA_character_
  )
}
