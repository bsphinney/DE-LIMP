#!/usr/bin/env Rscript
# ==============================================================================
#  qc_monitor.R
#  DE-LIMP QC Monitor — CLI script for automated ingestion of QC run results
#  into the Core Facility SQLite database.
#
#  This script is designed to run via cron on the HPC cluster. It does NOT
#  depend on Shiny — it uses only base R, DBI, RSQLite, arrow, and jsonlite.
#
#  Usage:
#    Rscript qc_monitor.R \
#      --core-dir /share/genome-center/delimp \
#      --watch-dir /share/genome-center/qc_outputs \
#      --instrument "timsTOF HT"
#
#  Options:
#    --core-dir   (required) Path to Core Facility directory containing delimp.db
#    --watch-dir  (required) Directory to scan for new QC report.parquet files
#    --instrument (required) Instrument name to record in qc_runs
#    --pattern    (optional) Filename pattern to match, default "report.parquet"
#    --dry-run    (optional) Show what would be ingested without writing to DB
#
#  Cron example (run every hour, log output):
#    0 * * * * Rscript /path/to/qc_monitor.R \
#      --core-dir /share/genome-center/delimp \
#      --watch-dir /share/genome-center/qc_outputs \
#      --instrument "timsTOF HT" \
#      >> /var/log/delimp_qc.log 2>&1
#
#  Apptainer note:
#    If running inside the DE-LIMP Apptainer container on HIVE, you may need
#    to set R_LIBS_USER so that R can find installed packages. For example:
#      export R_LIBS_USER=/share/genome-center/delimp/R_libs
#    Or wrap the cron entry:
#      0 * * * * R_LIBS_USER=/share/genome-center/delimp/R_libs \
#        apptainer exec /path/to/delimp.sif \
#        Rscript /path/to/qc_monitor.R --core-dir ... --watch-dir ... --instrument ...
#
# ==============================================================================

# Null-coalescing operator (standalone — no rlang dependency)
`%||%` <- function(x, y) if (is.null(x)) y else x

# ==============================================================================
#  1. Parse command-line arguments
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

# Simple argument parser — no optparse dependency
parse_args <- function(args) {
  parsed <- list(
    core_dir   = NULL,
    watch_dir  = NULL,
    instrument = NULL,
    pattern    = "report.parquet",
    dry_run    = FALSE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--core-dir" && i < length(args)) {
      parsed$core_dir <- args[i + 1]
      i <- i + 2
    } else if (arg == "--watch-dir" && i < length(args)) {
      parsed$watch_dir <- args[i + 1]
      i <- i + 2
    } else if (arg == "--instrument" && i < length(args)) {
      parsed$instrument <- args[i + 1]
      i <- i + 2
    } else if (arg == "--pattern" && i < length(args)) {
      parsed$pattern <- args[i + 1]
      i <- i + 2
    } else if (arg == "--dry-run") {
      parsed$dry_run <- TRUE
      i <- i + 1
    } else if (arg == "--help" || arg == "-h") {
      cat("Usage: Rscript qc_monitor.R --core-dir DIR --watch-dir DIR --instrument NAME [--pattern PAT] [--dry-run]\n")
      cat("\n")
      cat("Options:\n")
      cat("  --core-dir DIR      Path to Core Facility directory (contains delimp.db)\n")
      cat("  --watch-dir DIR     Directory to scan for QC report files\n")
      cat("  --instrument NAME   Instrument name (e.g., 'timsTOF HT')\n")
      cat("  --pattern PAT       Filename to match (default: report.parquet)\n")
      cat("  --dry-run           Show what would be ingested without writing to DB\n")
      cat("  --help, -h          Show this help message\n")
      quit(status = 0)
    } else {
      cat(sprintf("WARNING: Unknown argument '%s', ignoring.\n", arg))
      i <- i + 1
    }
  }

  parsed
}

opts <- parse_args(args)

# Validate required arguments
if (is.null(opts$core_dir)) {
  stop("ERROR: --core-dir is required. Use --help for usage.")
}
if (is.null(opts$watch_dir)) {
  stop("ERROR: --watch-dir is required. Use --help for usage.")
}
if (is.null(opts$instrument)) {
  stop("ERROR: --instrument is required. Use --help for usage.")
}

# Validate paths exist
if (!dir.exists(opts$core_dir)) {
  stop(sprintf("ERROR: Core directory does not exist: %s", opts$core_dir))
}
if (!dir.exists(opts$watch_dir)) {
  stop(sprintf("ERROR: Watch directory does not exist: %s", opts$watch_dir))
}

db_path <- file.path(opts$core_dir, "delimp.db")
if (!file.exists(db_path)) {
  stop(sprintf("ERROR: Database not found: %s", db_path))
}


# ==============================================================================
#  2. Load required packages
# ==============================================================================

required_pkgs <- c("DBI", "RSQLite", "arrow", "jsonlite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1),
                                       quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(sprintf("ERROR: Missing required packages: %s\nInstall with: install.packages(c(%s))",
               paste(missing_pkgs, collapse = ", "),
               paste(sprintf("'%s'", missing_pkgs), collapse = ", ")))
}

suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
  library(arrow)
  library(jsonlite)
})


# ==============================================================================
#  3. Scan watch directory for parquet files
# ==============================================================================

timestamp <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

cat(sprintf("%s QC Monitor started\n", timestamp()))
cat(sprintf("%s   Core dir:   %s\n", timestamp(), opts$core_dir))
cat(sprintf("%s   Watch dir:  %s\n", timestamp(), opts$watch_dir))
cat(sprintf("%s   Instrument: %s\n", timestamp(), opts$instrument))
cat(sprintf("%s   Pattern:    %s\n", timestamp(), opts$pattern))
if (opts$dry_run) {
  cat(sprintf("%s   *** DRY RUN MODE — no database writes ***\n", timestamp()))
}

# Recursively find matching files up to 2 levels deep.
# list.files recursive = TRUE has no depth limit, so we search each level
# explicitly to enforce the 2-level cap.
find_files <- function(base_dir, pattern, max_depth = 2) {
  found <- character(0)
  dirs_to_scan <- base_dir

  for (depth in 0:max_depth) {
    for (d in dirs_to_scan) {
      # Check for matching files in this directory
      matches <- list.files(d, pattern = paste0("^", pattern, "$"),
                            full.names = TRUE, recursive = FALSE)
      # Only keep actual files (not directories)
      matches <- matches[file.info(matches)$isdir == FALSE]
      found <- c(found, matches)
    }

    # Go one level deeper — list subdirectories
    if (depth < max_depth) {
      next_dirs <- character(0)
      for (d in dirs_to_scan) {
        subdirs <- list.dirs(d, full.names = TRUE, recursive = FALSE)
        next_dirs <- c(next_dirs, subdirs)
      }
      dirs_to_scan <- next_dirs
      if (length(dirs_to_scan) == 0) break
    }
  }

  # Normalize to absolute paths for consistent DB storage
  normalizePath(found, mustWork = FALSE)
}

all_files <- find_files(opts$watch_dir, opts$pattern)

cat(sprintf("%s Found %d '%s' file(s) in watch directory\n",
            timestamp(), length(all_files), opts$pattern))

if (length(all_files) == 0) {
  cat(sprintf("%s Nothing to do. Exiting.\n", timestamp()))
  quit(status = 0)
}


# ==============================================================================
#  4. Filter out files already in the database
# ==============================================================================

db <- dbConnect(SQLite(), db_path)
# Enable WAL mode (matches pattern in helpers_facility.R)
dbExecute(db, "PRAGMA journal_mode=WAL")

# Get all file_path values already recorded
existing_paths <- tryCatch({
  result <- dbGetQuery(db, "SELECT file_path FROM qc_runs WHERE file_path IS NOT NULL")
  normalizePath(result$file_path, mustWork = FALSE)
}, error = function(e) {
  cat(sprintf("%s WARNING: Could not query existing paths: %s\n", timestamp(), e$message))
  character(0)
})

new_files <- setdiff(all_files, existing_paths)

cat(sprintf("%s %d file(s) already in database, %d new file(s) to process\n",
            timestamp(), length(all_files) - length(new_files), length(new_files)))

if (length(new_files) == 0) {
  dbDisconnect(db)
  cat(sprintf("%s Nothing new to ingest. Exiting.\n", timestamp()))
  quit(status = 0)
}


# ==============================================================================
#  5. Extract metrics from each new parquet file and insert into DB
# ==============================================================================

#' Extract QC metrics from a DIA-NN report.parquet file.
#'
#' Handles varying column names across DIA-NN versions:
#'   - Protein.Group / Protein.Names / Protein.Ids
#'   - Stripped.Sequence / Modified.Sequence
#'   - Precursor.Id / Precursor.Quantity
#'   - Ms1.Area / MS1.Area
#'
#' @param file_path Absolute path to the parquet file
#' @return Named list of metrics, or NULL on failure
extract_metrics <- function(file_path) {
  tryCatch({
    # Read the parquet file
    df <- arrow::read_parquet(file_path)
    col_names <- names(df)

    # --- n_proteins ---
    protein_col <- intersect(
      c("Protein.Group", "Protein.Names", "Protein.Ids"),
      col_names
    )
    n_proteins <- if (length(protein_col) > 0) {
      length(unique(df[[protein_col[1]]]))
    } else {
      NA_integer_
    }

    # --- n_peptides ---
    peptide_col <- intersect(
      c("Stripped.Sequence", "Modified.Sequence"),
      col_names
    )
    n_peptides <- if (length(peptide_col) > 0) {
      length(unique(df[[peptide_col[1]]]))
    } else {
      NA_integer_
    }

    # --- n_precursors ---
    precursor_col <- intersect(c("Precursor.Id"), col_names)
    n_precursors <- if (length(precursor_col) > 0) {
      length(unique(df[[precursor_col[1]]]))
    } else {
      # Fallback: total rows (each row is typically one precursor measurement)
      nrow(df)
    }

    # --- median_ms1_tic ---
    ms1_col <- intersect(c("Ms1.Area", "MS1.Area"), col_names)
    median_ms1_tic <- if (length(ms1_col) > 0) {
      vals <- as.numeric(df[[ms1_col[1]]])
      vals <- vals[!is.na(vals) & vals > 0]
      if (length(vals) > 0) median(vals) else NA_real_
    } else {
      NA_real_
    }

    # --- Identify runs (samples) in the file ---
    run_col <- intersect(c("Run", "File.Name"), col_names)
    runs <- if (length(run_col) > 0) unique(df[[run_col[1]]]) else character(0)

    # --- Compute protein-level intensities and CV ---
    # Use Precursor.Quantity or similar for quantification
    quant_col <- intersect(
      c("Precursor.Quantity", "Precursor.Normalised", "Precursor.Translated",
        "Intensity", "Ms1.Area", "MS1.Area"),
      col_names
    )
    intensity_col <- if (length(quant_col) > 0) quant_col[1] else NULL

    # Initialize outputs
    median_cv <- NA_real_
    protein_intensities_json <- NA_character_

    if (!is.null(intensity_col) && length(protein_col) > 0) {
      pcol <- protein_col[1]
      icol <- intensity_col

      # Aggregate: sum precursor intensities per protein per run
      if (length(run_col) > 0 && length(runs) > 1) {
        # Multiple runs: compute CV across runs
        rcol <- run_col[1]

        # Use base R aggregate for reliability (no dplyr dependency)
        agg <- aggregate(
          as.numeric(df[[icol]]),
          by = list(protein = df[[pcol]], run = df[[rcol]]),
          FUN = sum, na.rm = TRUE
        )
        names(agg) <- c("protein", "run", "intensity")

        # Remove zero/NA intensities
        agg <- agg[!is.na(agg$intensity) & agg$intensity > 0, ]

        # Mean intensity per protein (for ranking and JSON output)
        mean_by_protein <- aggregate(
          agg$intensity,
          by = list(protein = agg$protein),
          FUN = mean, na.rm = TRUE
        )
        names(mean_by_protein) <- c("protein", "mean_intensity")

        # Top 500 by mean intensity
        mean_by_protein <- mean_by_protein[order(-mean_by_protein$mean_intensity), ]
        top_n <- min(500, nrow(mean_by_protein))
        top_proteins <- mean_by_protein$protein[seq_len(top_n)]

        # CV across runs for top proteins
        top_agg <- agg[agg$protein %in% top_proteins, ]
        cv_by_protein <- aggregate(
          top_agg$intensity,
          by = list(protein = top_agg$protein),
          FUN = function(x) {
            if (length(x) < 2) return(NA_real_)
            sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
          }
        )
        names(cv_by_protein) <- c("protein", "cv")
        cv_vals <- cv_by_protein$cv[!is.na(cv_by_protein$cv)]
        median_cv <- if (length(cv_vals) > 0) median(cv_vals) else NA_real_

        # Build JSON for top 500 protein intensities
        top_intensities <- mean_by_protein[seq_len(top_n), ]
        int_list <- as.list(setNames(
          round(top_intensities$mean_intensity, 2),
          top_intensities$protein
        ))
        protein_intensities_json <- jsonlite::toJSON(int_list, auto_unbox = TRUE)

      } else {
        # Single run (or no run column): compute mean per protein, CV = NA
        agg <- aggregate(
          as.numeric(df[[icol]]),
          by = list(protein = df[[pcol]]),
          FUN = sum, na.rm = TRUE
        )
        names(agg) <- c("protein", "intensity")
        agg <- agg[!is.na(agg$intensity) & agg$intensity > 0, ]
        agg <- agg[order(-agg$intensity), ]

        top_n <- min(500, nrow(agg))
        top_intensities <- agg[seq_len(top_n), ]

        int_list <- as.list(setNames(
          round(top_intensities$intensity, 2),
          top_intensities$protein
        ))
        protein_intensities_json <- jsonlite::toJSON(int_list, auto_unbox = TRUE)

        # median_cv stays NA for single-run files
      }
    }

    # --- run_date: file modification time ---
    run_date <- as.character(file.info(file_path)$mtime)

    # --- run_name: parent directory basename ---
    run_name <- basename(dirname(file_path))

    list(
      run_name            = run_name,
      run_date            = run_date,
      file_path           = file_path,
      n_proteins          = as.integer(n_proteins),
      n_peptides          = as.integer(n_peptides),
      n_precursors        = as.integer(n_precursors),
      median_ms1_tic      = median_ms1_tic,
      median_cv           = median_cv,
      protein_intensities = protein_intensities_json,
      n_runs              = length(runs)
    )
  }, error = function(e) {
    cat(sprintf("%s ERROR reading %s: %s\n", timestamp(), file_path, e$message))
    NULL
  })
}


# Process each new file
ingested <- 0
skipped  <- 0

for (fp in new_files) {
  cat(sprintf("%s Processing: %s\n", timestamp(), fp))

  metrics <- extract_metrics(fp)
  if (is.null(metrics)) {
    cat(sprintf("%s   SKIPPED (could not extract metrics)\n", timestamp()))
    skipped <- skipped + 1
    next
  }

  # Print summary
  cat(sprintf("%s   run_name:     %s\n", timestamp(), metrics$run_name))
  cat(sprintf("%s   run_date:     %s\n", timestamp(), metrics$run_date))
  cat(sprintf("%s   n_proteins:   %s\n", timestamp(),
              ifelse(is.na(metrics$n_proteins), "NA", as.character(metrics$n_proteins))))
  cat(sprintf("%s   n_peptides:   %s\n", timestamp(),
              ifelse(is.na(metrics$n_peptides), "NA", as.character(metrics$n_peptides))))
  cat(sprintf("%s   n_precursors: %s\n", timestamp(),
              ifelse(is.na(metrics$n_precursors), "NA", as.character(metrics$n_precursors))))
  cat(sprintf("%s   median_ms1:   %s\n", timestamp(),
              ifelse(is.na(metrics$median_ms1_tic), "NA", signif(metrics$median_ms1_tic, 4))))
  cat(sprintf("%s   median_cv:    %s\n", timestamp(),
              ifelse(is.na(metrics$median_cv), "NA (single run or insufficient data)",
                     sprintf("%.1f%%", metrics$median_cv))))
  cat(sprintf("%s   n_runs:       %d\n", timestamp(), metrics$n_runs))

  if (opts$dry_run) {
    cat(sprintf("%s   [DRY RUN] Would insert into qc_runs\n", timestamp()))
    ingested <- ingested + 1
    next
  }

  # Insert into database
  insert_result <- tryCatch({
    dbExecute(db, "
      INSERT INTO qc_runs (run_name, instrument, run_date, file_path,
                            n_proteins, n_peptides, n_precursors,
                            median_ms1_tic, median_cv, protein_intensities)
      VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)",
      params = list(
        metrics$run_name,
        opts$instrument,
        metrics$run_date,
        metrics$file_path,
        metrics$n_proteins %||% NA_integer_,
        metrics$n_peptides %||% NA_integer_,
        metrics$n_precursors %||% NA_integer_,
        metrics$median_ms1_tic %||% NA_real_,
        metrics$median_cv %||% NA_real_,
        metrics$protein_intensities %||% NA_character_
      )
    )
    TRUE
  }, error = function(e) {
    cat(sprintf("%s   ERROR inserting: %s\n", timestamp(), e$message))
    FALSE
  })

  if (isTRUE(insert_result)) {
    cat(sprintf("%s   Inserted into qc_runs OK\n", timestamp()))
    ingested <- ingested + 1
  } else {
    skipped <- skipped + 1
  }
}


# ==============================================================================
#  6. Summary
# ==============================================================================

dbDisconnect(db)

cat(sprintf("%s ---- Summary ----\n", timestamp()))
cat(sprintf("%s   Files found:    %d\n", timestamp(), length(all_files)))
cat(sprintf("%s   Already in DB:  %d\n", timestamp(), length(all_files) - length(new_files)))
cat(sprintf("%s   Newly ingested: %d\n", timestamp(), ingested))
cat(sprintf("%s   Skipped/errors: %d\n", timestamp(), skipped))
if (opts$dry_run) {
  cat(sprintf("%s   (Dry run — no actual writes performed)\n", timestamp()))
}
cat(sprintf("%s QC Monitor finished\n", timestamp()))
