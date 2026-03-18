# ==============================================================================
#  server_comparator.R
#  Run Comparator — Compare two analyses of the same dataset
#  Modes: DE-LIMP vs DE-LIMP, vs Spectronaut, vs FragPipe
#  Called from app.R as: server_comparator(input, output, session, values, add_to_log)
# ==============================================================================

# --- Pure helper functions (no Shiny reactivity) ---

#' Normalize protein IDs to bare UniProt accessions
normalize_protein_id <- function(ids) {
  # sp|P12345|GENE_HUMAN  ->  P12345
  ids <- gsub(".*\\|([A-Z][0-9][A-Z0-9]{3}[0-9])(-[0-9]+)?\\|.*", "\\1", ids)
  # P12345-2  ->  P12345  (strip isoform suffix)
  ids <- gsub("-[0-9]+$", "", ids)
  # P12345;Q67890  ->  P12345  (take first in group)
  ids <- gsub(";.*", "", ids)
  trimws(ids)
}

#' Detect Spectronaut version from export columns
detect_spectronaut_version <- function(df) {
  if ("EG.SpectroVersion" %in% names(df)) return(df[["EG.SpectroVersion"]][1])
  "Unknown"
}

#' Detect FragPipe version from export
detect_fp_version <- function(df) {
  if ("FragPipe.Version" %in% names(df)) return(df[["FragPipe.Version"]][1])
  "Unknown"
}

#' Detect intensity columns in FragPipe exports
detect_fp_intensity_cols <- function(df) {
  cols <- grep("MaxLFQ Intensity$| Intensity$", names(df), value = TRUE)
  cols[!grepl("^(Protein|Gene|Description|log|adj|p|Total|Unique|Combined)", cols, ignore.case = TRUE)]
}

#' Safe log2 transformation (0 -> NA)
log2_safe <- function(mat) {
  mat[mat == 0] <- NA
  log2(mat)
}

#' Coalesce setting value — returns val if non-null/non-"Not available", else fallback
coalesce_setting <- function(val, fallback = "Not available") {
  if (is.null(val) || identical(val, "Not available") || identical(val, "")) fallback else val
}

#' Format integer with commas or return "N/A"
format_or_na <- function(val) {
  if (is.null(val) || is.na(val)) return("N/A")
  format(as.integer(val), big.mark = ",")
}

#' Compute TopN effect summary — whether TopN drives quantification divergence
compute_topn_effect_summary <- function(joined, topn) {
  below <- mean(joined$quant_diff_abs[joined$n_peptides <= topn], na.rm = TRUE)
  above <- mean(joined$quant_diff_abs[joined$n_peptides >  topn], na.rm = TRUE)
  ratio <- above / below

  if (is.nan(ratio) || is.na(ratio)) return(NULL)

  if (ratio > 1.3) {
    interpretation <- sprintf(
      "Proteins with >%d peptides have %.1fx larger quantification divergence than proteins with <=%d peptides (%.2f vs %.2f log2 units). Spectronaut's TopN=%d cap is a primary driver of quantification differences.",
      topn, ratio, topn, above, below, topn)
    severity <- "warning"
  } else {
    interpretation <- sprintf(
      "Quantification divergence does not increase markedly above the TopN=%d threshold (%.2f vs %.2f log2 units). TopN may not be a major driver of discordance here.",
      topn, above, below)
    severity <- "info"
  }
  list(text = interpretation, severity = severity, ratio = ratio)
}

# --- Parsers ---

#' Parse a DE-LIMP session (from .rds file or current session reactive values)
parse_delimp_session <- function(rds_path = NULL, values = NULL) {
  if (!is.null(rds_path)) {
    session_obj <- readRDS(rds_path)
  } else if (!is.null(values)) {
    session_obj <- reactiveValuesToList(values)
  } else {
    stop("Either rds_path or values must be provided")
  }

  if (is.null(session_obj$fit) || is.null(session_obj$y_protein)) {
    stop("Session does not contain a completed DE-LIMP analysis (fit + y_protein required)")
  }

  contrasts <- colnames(session_obj$fit$contrasts)
  expr_mat <- as.data.frame(session_obj$y_protein$E)

  # Build protein ID column
  if ("Protein.Group" %in% colnames(session_obj$y_protein$genes)) {
    protein_ids <- normalize_protein_id(session_obj$y_protein$genes$Protein.Group)
  } else {
    protein_ids <- normalize_protein_id(rownames(expr_mat))
  }

  # Peptide counts
  n_pep <- NULL
  if ("NPeptides" %in% colnames(session_obj$y_protein$genes)) {
    n_pep <- data.frame(
      protein_id = protein_ids,
      n_peptides = as.integer(session_obj$y_protein$genes$NPeptides),
      stringsAsFactors = FALSE
    )
  }

  # Per-contrast DE stats
  de_stats <- lapply(contrasts, function(ct) {
    tbl <- limma::topTable(session_obj$fit, coef = ct, number = Inf)
    tbl <- as.data.frame(tbl)
    if (!"Protein.Group" %in% colnames(tbl)) {
      tbl$Protein.Group <- rownames(tbl)
    }
    tbl$protein_id <- normalize_protein_id(tbl$Protein.Group)
    tbl$contrast <- ct
    tbl
  })
  names(de_stats) <- contrasts

  # Extract search params if available
  # If diann_search_settings is NULL/empty, try to find settings from the most recent completed job
  ss <- session_obj$diann_search_settings
  if (is.null(ss) || length(ss) == 0) {
    jobs <- session_obj$diann_jobs
    if (!is.null(jobs) && length(jobs) > 0) {
      # Find most recent completed job with search_settings
      for (j in rev(seq_along(jobs))) {
        if (!is.null(jobs[[j]]$search_settings) && length(jobs[[j]]$search_settings) > 0 &&
            identical(jobs[[j]]$status, "completed")) {
          ss <- jobs[[j]]$search_settings
          message("[DE-LIMP] parse_delimp_session: recovered search_settings from job '",
                  jobs[[j]]$name %||% "?", "'")
          break
        }
      }
    }
  }
  sp <- ss$search_params

  # Extract DIA-NN version: prefer search_settings, then parse from job log_content, then SIF name
  diann_version <- ss$diann_version
  if (is.null(diann_version) || !nzchar(diann_version %||% "") || diann_version == "unknown") {
    # Try parsing from job log_content (authoritative — DIA-NN prints its own version)
    jobs <- session_obj$diann_jobs
    if (!is.null(jobs)) {
      for (j in rev(seq_along(jobs))) {
        log_text <- jobs[[j]]$log_content %||% ""
        if (nzchar(log_text)) {
          log_lines <- strsplit(log_text, "\n")[[1]]
          ver_hit <- regmatches(log_lines,
            regexpr("DIA-NN\\s+([0-9]+\\.[0-9]+\\.?[0-9]*)", log_lines))
          ver_hit <- ver_hit[nzchar(ver_hit)]
          if (length(ver_hit) > 0) {
            diann_version <- sub("DIA-NN\\s+", "", ver_hit[1])
            message("[DE-LIMP] parse_delimp_session: parsed DIA-NN version ", diann_version, " from job log")
            break
          }
        }
      }
    }
  }
  # Last resort: extract from SIF filename (less reliable)
  if (is.null(diann_version) || !nzchar(diann_version %||% "") || diann_version == "unknown") {
    sif <- ss$diann_sif %||% ""
    ver_match <- regmatches(sif, regexec("(\\d+\\.\\d+\\.\\d+)", sif))[[1]]
    if (length(ver_match) > 1) {
      diann_version <- ver_match[2]
      message("[DE-LIMP] parse_delimp_session: extracted DIA-NN version ", diann_version, " from SIF name (fallback)")
    }
  }

  # Data-level stats from raw_data
  n_precursors <- nrow(session_obj$raw_data$E) %||% NA
  n_samples    <- ncol(session_obj$raw_data$E) %||% NA
  n_proteins   <- length(protein_ids)

  # Precursor m/z range from raw_data$genes
  mz_range <- tryCatch({
    if (!is.null(session_obj$raw_data$genes) && "Precursor.Mz" %in% colnames(session_obj$raw_data$genes)) {
      mz <- session_obj$raw_data$genes$Precursor.Mz
      paste0(round(min(mz, na.rm = TRUE), 1), " - ", round(max(mz, na.rm = TRUE), 1))
    } else NA
  }, error = function(e) NA)

  # Charge range
  charge_range <- tryCatch({
    if (!is.null(session_obj$raw_data$genes) && "Precursor.Charge" %in% colnames(session_obj$raw_data$genes)) {
      ch <- session_obj$raw_data$genes$Precursor.Charge
      paste0(min(ch, na.rm = TRUE), " - ", max(ch, na.rm = TRUE))
    } else NA
  }, error = function(e) NA)

  # Dynamic range from raw (precursor-level) intensities — log2 scale
  dynamic_range <- tryCatch({
    raw_mat <- session_obj$raw_data$E
    if (!is.null(raw_mat) && length(raw_mat) > 0) {
      # Per-sample dynamic range
      per_sample <- apply(raw_mat, 2, function(col) {
        v <- col[!is.na(col)]
        if (length(v) < 10) return(list(min = NA, max = NA, median = NA, iqr = NA, orders = NA, n = 0, pct_missing = 100))
        list(
          min = min(v), max = max(v), median = median(v), iqr = IQR(v),
          orders = round((max(v) - min(v)) / log2(10), 2),
          n = length(v),
          pct_missing = round(100 * mean(is.na(col)), 1)
        )
      })
      # Global summary
      all_vals <- raw_mat[!is.na(raw_mat)]
      list(
        global_min = min(all_vals), global_max = max(all_vals),
        global_median = median(all_vals), global_iqr = IQR(all_vals),
        global_orders = round((max(all_vals) - min(all_vals)) / log2(10), 2),
        n_precursors = nrow(raw_mat),
        pct_missing = round(100 * mean(is.na(raw_mat)), 1),
        per_sample = per_sample
      )
    } else NULL
  }, error = function(e) NULL)

  # FASTA info
  fasta_entries <- ss$fasta_seq_count %||% session_obj$fasta_info$n_sequences %||% NA
  fasta_file    <- if (length(ss$fasta_files %||% character()) > 0) paste(basename(ss$fasta_files), collapse = ", ") else NA

  list(
    source      = "delimp",
    de_stats    = de_stats,
    intensities = expr_mat,
    protein_ids = protein_ids,
    n_peptides  = n_pep,
    missing_pct = setNames(
      apply(session_obj$y_protein$E, 1, function(x) mean(is.na(x))),
      protein_ids
    ),
    metadata    = session_obj$metadata,
    dynamic_range = dynamic_range,
    settings    = list(
      software        = "DE-LIMP",
      delimp_version  = session_obj$app_version %||% "unknown",
      limpa_version   = tryCatch(as.character(packageVersion("limpa")), error = function(e) "unknown"),
      limma_version   = tryCatch(as.character(packageVersion("limma")), error = function(e) "unknown"),
      diann_version   = diann_version %||% session_obj$diann_version %||% "unknown",
      de_significance = session_obj$fdr_threshold %||% "0.05",
      identification_fdr = "0.01 (DIA-NN default)",
      lfc_threshold   = session_obj$lfc_threshold %||% "0.6",
      min_peptides    = session_obj$min_peptides %||% "unknown",
      normalization   = "DIA-NN RT-dependent (Precursor.Normalised)",
      rollup_method   = "DPC-Quant (empirical Bayes precursor aggregation)",
      de_engine       = "limma moderated t-test",
      covariates      = paste(session_obj$covariates %||% "none", collapse = ", "),
      contrast_string = session_obj$contrast_string %||% "unknown",
      # Search & data stats
      n_precursors    = as.character(n_precursors %||% "unknown"),
      n_proteins_total = as.character(n_proteins),
      n_samples       = as.character(n_samples %||% "unknown"),
      precursor_mz_range = as.character(mz_range %||% "unknown"),
      charge_range    = as.character(charge_range %||% "unknown"),
      fasta_file      = as.character(fasta_file %||% "unknown"),
      fasta_entries   = as.character(fasta_entries %||% "unknown"),
      mass_acc_ms2    = as.character(sp$mass_acc %||% "auto"),
      mass_acc_ms1    = as.character(sp$mass_acc_ms1 %||% "auto"),
      scan_window     = as.character(sp$scan_window %||% "auto"),
      enzyme          = as.character(sp$enzyme %||% "Trypsin/P"),
      missed_cleavages = as.character(sp$missed_cleavages %||% "unknown"),
      search_mz_min   = as.character(sp$min_pr_mz %||% "300"),
      search_mz_max   = as.character(sp$max_pr_mz %||% "1800"),
      min_fr_mz       = as.character(sp$min_fr_mz %||% "200"),
      max_fr_mz       = as.character(sp$max_fr_mz %||% "1800"),
      min_pr_charge   = as.character(sp$min_pr_charge %||% "1"),
      max_pr_charge   = as.character(sp$max_pr_charge %||% "4"),
      min_pep_len     = as.character(sp$min_pep_len %||% "7"),
      max_pep_len     = as.character(sp$max_pep_len %||% "30"),
      var_mods        = paste(Filter(nzchar, c(
        if (isTRUE(sp$mod_met_ox)) "UniMod:35 (Met oxidation)",
        if (isTRUE(sp$mod_nterm_acetyl)) "UniMod:1 (N-term acetylation)",
        if (nzchar(sp$extra_var_mods %||% "")) sp$extra_var_mods
      )), collapse = "; "),
      fixed_mods      = if (isTRUE(sp$unimod4)) "UniMod:4 (Carbamidomethylation)" else "none",
      max_var_mods    = as.character(sp$max_var_mods %||% "1"),
      mbr             = as.character(sp$mbr %||% "unknown"),
      search_mode     = as.character(ss$search_mode %||% "libfree"),
      library         = as.character(sp$library %||% ss$library_path %||% "FASTA-predicted"),
      # Dynamic range (from raw precursor-level intensities, log2 scale)
      dynamic_range_orders    = if (!is.null(dynamic_range)) as.character(dynamic_range$global_orders) else NA,
      dynamic_range_log2_min  = if (!is.null(dynamic_range)) as.character(round(dynamic_range$global_min, 1)) else NA,
      dynamic_range_log2_max  = if (!is.null(dynamic_range)) as.character(round(dynamic_range$global_max, 1)) else NA,
      dynamic_range_median    = if (!is.null(dynamic_range)) as.character(round(dynamic_range$global_median, 1)) else NA,
      dynamic_range_iqr       = if (!is.null(dynamic_range)) as.character(round(dynamic_range$global_iqr, 1)) else NA,
      dynamic_range_pct_missing = if (!is.null(dynamic_range)) as.character(dynamic_range$pct_missing) else NA
    ),
    contrasts   = contrasts
  )
}

#' Parse Spectronaut Candidates.tsv into per-comparison DE results
#' @param de_path Path to Candidates.tsv
#' @return Named list of data.frames (one per comparison), or NULL
parse_spectronaut_candidates <- function(de_path) {
  de_df <- data.table::fread(de_path, data.table = TRUE)
  de_prot_col  <- grep("ProteinGroup|ProteinAccession|UniProtIds|^Group$", names(de_df), value = TRUE)[1]
  de_gene_col  <- grep("^PG\\.Genes$|^Gene$|Genes", names(de_df), value = TRUE)[1]
  de_logfc_col <- grep("Log2Ratio|log2FC|AVG\\.Log2\\.Ratio|AVG Log2 Ratio", names(de_df), value = TRUE)[1]
  de_qval_col  <- grep("Qvalue|Q\\.Value|q\\.value|adj", names(de_df), ignore.case = TRUE, value = TRUE)[1]
  de_pval_col  <- grep("^.*Pvalue$|^.*p\\.value$|PValue", names(de_df), ignore.case = TRUE, value = TRUE)[1]
  nratios_col  <- grep("# of Ratios|NrOfRatios|NumberOfRatios", names(de_df), value = TRUE)[1]

  if (is.na(de_prot_col)) stop("Candidates file missing protein column")
  if (is.na(de_logfc_col)) stop("Candidates file missing Log2Ratio / logFC column")

  # Handle multiple comparisons
  comp_col <- grep("^Comparison|^R\\.Condition$|^Contrast$", names(de_df), value = TRUE)[1]
  build_de_frame <- function(sub_df) {
    sub_ids <- normalize_protein_id(sub_df[[de_prot_col]])
    out <- data.frame(
      protein_id = sub_ids,
      gene       = if (!is.na(de_gene_col)) as.character(sub_df[[de_gene_col]]) else sub_ids,
      logFC      = as.numeric(sub_df[[de_logfc_col]]),
      P.Value    = as.numeric(if (!is.na(de_pval_col)) sub_df[[de_pval_col]] else NA),
      adj.P.Val  = as.numeric(if (!is.na(de_qval_col)) sub_df[[de_qval_col]] else NA),
      stringsAsFactors = FALSE
    )
    if (!is.na(nratios_col)) out$n_ratios <- as.integer(sub_df[[nratios_col]])
    out
  }

  if (!is.na(comp_col)) {
    comparisons <- unique(de_df[[comp_col]])
    de_stats_list <- lapply(comparisons, function(comp) {
      build_de_frame(de_df[de_df[[comp_col]] == comp, ])
    })
    names(de_stats_list) <- comparisons
    de_stats_list
  } else {
    list(comparison = build_de_frame(de_df))
  }
}

#' Parse Spectronaut RunSummaries directory into per-sample QC table
#' @param run_summary_paths Character vector of paths to RunSummary TSV files
#' @return data.frame with per-sample QC metrics, or NULL
parse_spectronaut_run_summaries <- function(run_summary_paths) {
  rows <- lapply(run_summary_paths, function(p) {
    tryCatch({
      df <- data.table::fread(p, data.table = FALSE)
      # RunSummary files are typically key-value pairs or single-row tables
      # Try single-row table first (most common Spectronaut export)
      if (nrow(df) >= 1 && ncol(df) > 3) {
        find_col <- function(patterns) {
          for (pat in patterns) {
            m <- grep(pat, names(df), ignore.case = TRUE, value = TRUE)
            if (length(m) > 0) return(m[1])
          }
          NA_character_
        }
        data.frame(
          file_name       = basename(tools::file_path_sans_ext(p)),
          precursors      = as.integer(df[[find_col(c("Precursors$", "^EG\\.Count"))]] %||% NA),
          peptides        = as.integer(df[[find_col(c("Peptides$", "StrippedSequences"))]] %||% NA),
          protein_groups  = as.integer(df[[find_col(c("Protein.?Groups?$", "^PG\\.Count"))]] %||% NA),
          avg_pep_per_pg  = as.numeric(df[[find_col(c("AVG Peptides", "AvgPeptides"))]] %||% NA),
          ms1_ppm         = as.numeric(df[[find_col(c("MS1.*Delta.*ppm", "MS1.*Mass.*Acc"))]] %||% NA),
          ms2_ppm         = as.numeric(df[[find_col(c("MS2.*Delta.*ppm", "MS2.*Mass.*Acc"))]] %||% NA),
          instrument      = as.character(df[[find_col(c("Instrument.?Model", "^R\\.Instrument"))]] %||% NA),
          raw_file        = as.character(df[[find_col(c("Raw.?File.?Name", "^R\\.RawFileName", "^R\\.FileName"))]] %||% NA),
          gradient_min    = as.numeric(df[[find_col(c("Gradient.?Length", "^R\\.GradientLength"))]] %||% NA),
          stringsAsFactors = FALSE
        )
      } else NULL
    }, error = function(e) NULL)
  })
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

#' Parse Spectronaut ConditionSetup.tsv into sample map
#' @param path Path to ConditionSetup.tsv
#' @return data.frame with columns: run_label, condition, file_name
parse_spectronaut_condition_setup <- function(path) {
  df <- data.table::fread(path, data.table = FALSE)
  # Expected columns: Run Label, Condition, Replicate, File Name (Spectronaut standard)
  label_col <- grep("Run.?Label", names(df), ignore.case = TRUE, value = TRUE)[1]
  cond_col  <- grep("^Condition$", names(df), ignore.case = TRUE, value = TRUE)[1]
  file_col  <- grep("File.?Name", names(df), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(file_col)) stop("ConditionSetup.tsv missing 'File Name' column")
  if (is.na(cond_col)) stop("ConditionSetup.tsv missing 'Condition' column")

  data.frame(
    run_label  = if (!is.na(label_col)) as.character(df[[label_col]]) else as.character(df[[file_col]]),
    condition  = as.character(df[[cond_col]]),
    file_name  = as.character(df[[file_col]]),
    stringsAsFactors = FALSE
  )
}

#' Parse Spectronaut ExperimentSetupOverview text file into settings list
#' @param path Path to ExperimentSetupOverview.txt
#' @return Named list of settings key-value pairs
parse_spectronaut_setup_overview <- function(path) {
  lines <- readLines(path, warn = FALSE)
  # Strip tree drawing characters (│ ├─ └─ ─) and leading whitespace
  lines <- gsub("[\u2502\u251c\u2514\u2500\u2510\u250c\u2518\u2524\u252c\u2534\u253c]", "", lines)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  settings <- list()
  # First line often has version: "Spectronaut 20.5.260227.92449"
  raw_first <- trimws(readLines(path, n = 1, warn = FALSE))
  ver_match <- regmatches(raw_first, regexpr("Spectronaut [0-9][0-9.]+", raw_first))
  if (length(ver_match) > 0 && nzchar(ver_match)) settings[["spectronaut_version"]] <- ver_match

  for (line in lines) {
    # Try "Key: Value" or "Key\tValue" patterns
    if (grepl(":\t|:\\s{2,}", line)) {
      parts <- strsplit(line, ":\t|:\\s{2,}")[[1]]
      if (length(parts) == 2) settings[[trimws(parts[1])]] <- trimws(parts[2])
    } else if (grepl("^[^:]+:\\s+.+", line)) {
      key <- sub(":.*", "", line)
      val <- sub("^[^:]+:\\s+", "", line)
      settings[[trimws(key)]] <- trimws(val)
    }
  }
  # Extract TopN Max values contextually (flat parser can't distinguish "Max" under different parents)
  raw_lines_topn <- readLines(path, warn = FALSE)
  raw_topn <- gsub("[\u2502\u251c\u2514\u2500\u2510\u250c\u2518\u2524\u252c\u2534\u253c]", "", raw_lines_topn)
  raw_topn <- trimws(raw_topn)
  in_major_topn <- FALSE
  in_minor_topn <- FALSE
  for (line in raw_topn) {
    if (grepl("^Major Group Top N", line)) { in_major_topn <- TRUE; in_minor_topn <- FALSE; next }
    if (grepl("^Minor Group Top N", line)) { in_minor_topn <- TRUE; in_major_topn <- FALSE; next }
    if (grepl("^Max", line)) {
      val <- sub("^Max[:\\s\t]+", "", line)
      val <- trimws(val)
      if (in_major_topn) { settings$topn_protein_max <- val; in_major_topn <- FALSE }
      else if (in_minor_topn) { settings$topn_peptide_max <- val; in_minor_topn <- FALSE }
    }
    # Reset context if we hit a different key at same level
    if (!grepl("^(Max|Min)$", sub("[:\\s\t].*", "", line))) {
      in_major_topn <- FALSE
      in_minor_topn <- FALSE
    }
  }

  # Collect all FASTA databases (Original File appears multiple times under Protein Databases Used)
  raw_lines <- readLines(path, warn = FALSE)
  raw_clean <- gsub("[\u2502\u251c\u2514\u2500\u2510\u250c\u2518\u2524\u252c\u2534\u253c]", "", raw_lines)
  raw_clean <- trimws(raw_clean)
  fasta_files <- character(0)
  fasta_entries <- integer(0)
  for (i in seq_along(raw_clean)) {
    if (grepl("^Original File", raw_clean[i])) {
      val <- sub("^Original File[:\\s\t]+", "", raw_clean[i])
      val <- trimws(val)
      if (nzchar(val) && grepl("\\.(fasta|fa)$", val, ignore.case = TRUE)) {
        fasta_files <- c(fasta_files, val)
      }
    }
    if (grepl("^Protein Entries", raw_clean[i])) {
      val <- sub("^Protein Entries[:\\s\t]+", "", raw_clean[i])
      val <- trimws(val)
      n <- suppressWarnings(as.integer(gsub(",", "", val)))
      if (!is.na(n)) fasta_entries <- c(fasta_entries, n)
    }
  }
  # Deduplicate — same FASTA listed per sample
  if (length(fasta_files) > 0 && length(fasta_files) == length(fasta_entries)) {
    fasta_df <- unique(data.frame(file = fasta_files, entries = fasta_entries, stringsAsFactors = FALSE))
    fasta_files <- fasta_df$file
    fasta_entries <- fasta_df$entries
  } else if (length(fasta_files) > 0) {
    fasta_files <- unique(fasta_files)
  }
  settings$fasta_databases <- if (length(fasta_files) > 0) {
    paste(fasta_files, collapse = " + ")
  } else NULL
  settings$fasta_total_entries <- if (length(fasta_entries) > 0) {
    format(sum(fasta_entries), big.mark = ",")
  } else NULL

  message("[Comparator] Parsed ", length(settings), " settings from ExperimentSetupOverview")
  if (length(fasta_files) > 0) {
    message("[Comparator] FASTA databases: ", paste(fasta_files, "(", format(fasta_entries, big.mark = ","), "entries)", collapse = " + "))
  }
  settings
}

#' Extract typed search settings from ExperimentSetupOverview key-value pairs
#' @param overview Named list from parse_spectronaut_setup_overview()
#' @return Named list of search settings with standardized keys
parse_spectronaut_search_settings <- function(overview) {
  get_val <- function(...) {
    keys <- c(...)
    for (k in keys) {
      matches <- grep(k, names(overview), ignore.case = TRUE, value = TRUE)
      if (length(matches) > 0) return(overview[[matches[1]]])
    }
    "Not available"
  }
  list(
    missed_cleavages   = get_val("^Missed Cleavages$"),
    max_peptide_length = get_val("^Max Peptide Length$"),
    min_peptide_length = get_val("^Min Peptide Length$"),
    variable_mods      = get_val("^Variable Modifications"),
    normalization      = get_val("^Normalization Strategy$"),
    cross_run_norm     = get_val("^Cross-Run Normalization$"),
    topn_protein       = get_val("^Major Group Top N$"),
    topn_peptide       = get_val("^Minor Group Top N$"),
    de_test            = get_val("^Differential Abundance Testing$"),
    lfc_filter         = get_val("^Log2 Ratio Candidate Filter$"),
    confidence_filter  = get_val("^Confidence$"),
    interference_corr  = get_val("^Interference Correction$"),
    imputation         = get_val("^Imputation Strategy$"),
    use_all_ms_levels  = get_val("^Use All MS-Level Quantities$"),
    quantity_ms_level  = get_val("^Quantity MS Level$"),
    quantity_type      = get_val("^Quantity Type$"),
    protein_lfq_method = get_val("^Protein LFQ Method$"),
    # FDR thresholds
    protein_qvalue     = get_val("^Protein Qvalue Cutoff \\(Experiment\\)$", "^Protein Group FDR$"),
    precursor_qvalue   = get_val("^Precursor Qvalue Cutoff$"),
    # Enzyme / search
    enzyme             = get_val("^Enzymes / Cleavage Rules$", "^Enzyme$", "^Protease$"),
    digest_type        = get_val("^Digest Type$")
  )
}

#' Parse Spectronaut AnalysisLog.txt for library stats and version
#' @param log_lines Character vector of log file lines
#' @return Named list with version and library composition
parse_spectronaut_log <- function(log_lines) {
  extract_int <- function(pattern) {
    m <- grep(pattern, log_lines, value = TRUE, ignore.case = TRUE)
    if (length(m) == 0) return(NA_integer_)
    as.integer(gsub("[^0-9]", "", m[length(m)]))
  }
  extract_val <- function(pattern) {
    m <- grep(pattern, log_lines, value = TRUE, ignore.case = TRUE)
    if (length(m) == 0) return(NULL)
    # Extract value after the last colon (handles tree-style "├─ Key:    Value")
    val <- sub("^.*:\\s+", "", trimws(m[1]))
    if (nzchar(val)) val else NULL
  }
  list(
    spectronaut_version = {
      m <- grep("Spectronaut [0-9]", log_lines, value = TRUE)
      if (length(m) > 0) sub(".*(Spectronaut [0-9.]+).*", "\\1", m[1]) else "Unknown"
    },
    library_precursors  = {
      v <- extract_int("precursors.*library|library.*precursors")
      # Fallback: summary line "Precursors: 94202" (no "library" keyword)
      if (is.na(v)) {
        m <- grep("Precursors:\\s*[0-9,]+", log_lines, value = TRUE)
        if (length(m) > 0) {
          nums <- regmatches(m[1], gregexpr("[0-9,]+", m[1]))[[1]]
          # Take the number right after "Precursors:"
          prec_idx <- grep("Precursors", strsplit(m[1], ",")[[1]])
          if (length(prec_idx) > 0) {
            chunk <- strsplit(m[1], ",")[[1]][prec_idx[1]]
            v <- as.integer(gsub("[^0-9]", "", sub(".*Precursors:?\\s*", "", chunk)))
          }
        }
      }
      v
    },
    library_peptides    = extract_int("peptides.*library|library.*peptides"),
    library_proteins    = extract_int("protein groups.*library|library.*protein groups"),
    unique_precursors   = extract_int("unique precursor"),
    unique_proteins     = extract_int("unique protein"),
    pulsar_candidates   = extract_int("candidate peptides|Pulsar"),
    # Settings from log tree
    normalization_strategy = extract_val("Normalization Strategy"),
    cross_run_normalization = extract_val("Cross-Run Normalization"),
    normalization_filter   = extract_val("Normalization Filter Type"),
    use_all_ms_levels      = extract_val("Use All MS-Level Quantities"),
    topn_protein           = extract_val("Major Group Top N"),
    topn_peptide           = extract_val("Minor Group Top N"),
    de_test                = extract_val("Differential Abundance.*Test|Statistical Testing"),
    imputation             = extract_val("Imputation")
  )
}

#' Parse Spectronaut ZIP export (primary Mode B input)
#' Detects and parses: ConditionSetup.tsv, Pivot quantities, Candidates.tsv,
#' RunSummaries/*.tsv, ExperimentSetupOverview.txt
#' @param zip_path Path to the ZIP file
#' @return Named list compatible with comparator pipeline (same structure as parse_spectronaut)
parse_spectronaut_zip <- function(zip_path) {
  tmp <- tempdir()
  zip_dir <- file.path(tmp, paste0("spectronaut_zip_", format(Sys.time(), "%H%M%S")))
  dir.create(zip_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(zip_dir, recursive = TRUE), add = TRUE)

  utils::unzip(zip_path, exdir = zip_dir)
  all_files <- list.files(zip_dir, recursive = TRUE, full.names = TRUE)
  basenames <- tolower(basename(all_files))
  message("[Comparator] ZIP contains ", length(all_files), " files: ",
          paste(basename(all_files), collapse = ", "))

  # --- Detect files ---
  detect <- function(patterns, files = all_files, bases = basenames) {
    for (pat in patterns) {
      idx <- grep(pat, bases)
      if (length(idx) > 0) return(files[idx[1]])
    }
    NULL
  }

  condition_file <- detect(c("conditionsetup\\.tsv$", "conditionsetup\\.csv$"))
  pivot_file     <- detect(c("de-limp.*pivot.*\\.tsv$", "pivot.*\\.tsv$",
                             "bgs factory report.*\\.tsv$",
                             "de-limp.*pivot.*\\.csv$", "pivot.*\\.csv$",
                             "bgs factory report.*\\.csv$"))
  candidates_file <- detect(c("^candidates\\.tsv$", "^candidates\\.csv$",
                              "candidates.*\\.tsv$"))
  setup_file     <- detect(c("experimentsetupoverview.*\\.txt$",
                             "experimentsetupoverview.*\\.tsv$"))
  analysis_log   <- detect(c("^analysislog\\.txt$", "^analysislog\\.log$",
                             "spectronaut.*log\\.txt$", ".*\\.log\\.txt$",
                             "analysis.*log.*\\.txt$"))
  # AnalysisOverview (note: Spectronaut typos it as "AnalyisOverview")
  analysis_overview_file <- detect(c("analy.?is.?overview.*\\.txt$",
                                     "analysisoverview.*\\.txt$"))

  # RunSummaries: look for a RunSummaries directory or matching files
  run_summary_files <- all_files[grep("runsummar", basenames)]
  # If there's a RunSummaries subdirectory, get all .tsv files inside it
  rs_dirs <- list.dirs(zip_dir, recursive = TRUE, full.names = TRUE)
  rs_dir <- rs_dirs[grep("RunSummar", basename(rs_dirs), ignore.case = TRUE)]
  if (length(rs_dir) > 0) {
    rs_candidates <- list.files(rs_dir[1], pattern = "\\.(tsv|csv)$",
                                full.names = TRUE, recursive = FALSE)
    if (length(rs_candidates) > 0) run_summary_files <- rs_candidates
  }

  # Build detection manifest
  manifest <- list(
    condition_setup = !is.null(condition_file),
    quantities      = !is.null(pivot_file),
    candidates      = !is.null(candidates_file),
    run_summaries   = length(run_summary_files) > 0,
    n_run_summaries = length(run_summary_files),
    settings           = !is.null(setup_file),
    analysis_log       = !is.null(analysis_log),
    analysis_overview  = !is.null(analysis_overview_file)
  )

  if (!manifest$quantities) stop("ZIP must contain a Pivot quantities file (*Pivot*.tsv or *BGS Factory Report*.tsv)")

  # --- Parse each component ---

  # 1. ConditionSetup (sample map)
  sample_map <- if (manifest$condition_setup) {
    tryCatch(parse_spectronaut_condition_setup(condition_file),
             error = function(e) { message("Warning: Could not parse ConditionSetup: ", e$message); NULL })
  } else NULL

  # 2. Pivot quantities (reuse existing column-detection logic)
  df <- data.table::fread(pivot_file, data.table = TRUE)
  protein_col <- grep("ProteinGroup|ProteinAccession|UniProtIds", names(df), value = TRUE)[1]
  gene_col    <- grep("^PG\\.Genes$|^Gene$|Genes", names(df), value = TRUE)[1]
  quant_cols  <- grep("PG\\.Quantity$", names(df), value = TRUE)
  npep_col    <- grep("NrOfStrippedSequences", names(df), value = TRUE)[1]

  if (is.na(protein_col)) stop("Pivot file missing protein column (PG.ProteinGroups)")
  if (length(quant_cols) == 0) stop("Pivot file has no quantity columns (*.PG.Quantity)")

  protein_ids <- normalize_protein_id(df[[protein_col]])
  quant_mat <- as.data.frame(df[, quant_cols, with = FALSE])

  # If we have a sample_map, rename quant columns from Spectronaut labels to File Name keys
  # This makes match_samples() work against DE-LIMP colnames which are file-based
  if (!is.null(sample_map)) {
    # Spectronaut pivot columns are like "[N] RunLabel.PG.Quantity" — extract the prefix
    col_prefixes <- sub("\\.PG\\.Quantity$", "", quant_cols)
    # Strip Spectronaut's "[N] " column numbering prefix (e.g., "[1] AD-1" → "AD-1")
    col_prefixes <- sub("^\\[\\d+\\]\\s*", "", col_prefixes)
    # Match to sample_map run_label to get file_name
    map_idx <- match(col_prefixes, sample_map$run_label)
    new_names <- ifelse(is.na(map_idx), col_prefixes, sample_map$file_name[map_idx])
    colnames(quant_mat) <- new_names
  } else {
    # Strip .PG.Quantity suffix and [N] prefix for cleaner column names
    clean_names <- sub("\\.PG\\.Quantity$", "", quant_cols)
    clean_names <- sub("^\\[\\d+\\]\\s*", "", clean_names)
    colnames(quant_mat) <- clean_names
  }

  # Inline DE columns in pivot (single-file mode fallback)
  logfc_col <- grep("Log2Ratio|log2FC", names(df), value = TRUE)[1]
  qval_col  <- grep("Qvalue|Q\\.Value|q\\.value|adj", names(df), ignore.case = TRUE, value = TRUE)[1]
  pval_col  <- grep("^.*Pvalue$|^.*p\\.value$", names(df), ignore.case = TRUE, value = TRUE)[1]
  has_inline_de <- !is.na(logfc_col) && !is.na(pval_col)

  # 3. Candidates (DE results) — prefer separate file over inline
  de_stats_df <- NULL
  has_de <- FALSE
  if (manifest$candidates) {
    tryCatch({
      de_stats_df <- parse_spectronaut_candidates(candidates_file)
      has_de <- TRUE
    }, error = function(e) message("Warning: Could not parse Candidates: ", e$message))
  }
  if (!has_de && has_inline_de) {
    de_stats_df <- list(comparison = data.frame(
      protein_id = protein_ids,
      gene       = if (!is.na(gene_col)) as.character(df[[gene_col]]) else protein_ids,
      logFC      = as.numeric(df[[logfc_col]]),
      P.Value    = as.numeric(if (!is.na(pval_col)) df[[pval_col]] else NA),
      adj.P.Val  = as.numeric(if (!is.na(qval_col)) df[[qval_col]] else NA),
      stringsAsFactors = FALSE
    ))
    has_de <- TRUE
  }

  # 4. RunSummaries (per-sample QC)
  run_qc <- if (manifest$run_summaries) {
    tryCatch(parse_spectronaut_run_summaries(run_summary_files),
             error = function(e) { message("Warning: Could not parse RunSummaries: ", e$message); NULL })
  } else NULL

  message("[Comparator] Manifest: quantities=", manifest$quantities,
          ", candidates=", manifest$candidates,
          ", condition_setup=", manifest$condition_setup,
          ", settings=", manifest$settings,
          ", analysis_log=", manifest$analysis_log,
          ", analysis_overview=", manifest$analysis_overview,
          ", run_summaries=", manifest$run_summaries, " (", manifest$n_run_summaries, ")")

  # 5. ExperimentSetupOverview (settings)
  setup_overview <- if (manifest$settings) {
    tryCatch(parse_spectronaut_setup_overview(setup_file),
             error = function(e) { message("Warning: Could not parse SetupOverview: ", e$message); NULL })
  } else NULL
  if (!is.null(setup_overview)) {
    message("[Comparator] ExperimentSetupOverview keys: ", paste(names(setup_overview), collapse = ", "))
  }

  # 6. AnalysisLog.txt (Spectronaut version + library composition)
  library_info <- if (manifest$analysis_log) {
    tryCatch({
      log_lines <- readLines(analysis_log, warn = FALSE)
      parse_spectronaut_log(log_lines)
    }, error = function(e) { message("Warning: Could not parse AnalysisLog: ", e$message); NULL })
  } else NULL

  # 6b. AnalysisOverview (precursor/protein counts, data completeness)
  analysis_overview <- if (manifest$analysis_overview) {
    tryCatch({
      ao_lines <- readLines(analysis_overview_file, warn = FALSE)
      ao <- list()
      for (line in ao_lines) {
        if (grepl("^\t|^$", line)) next
        parts <- strsplit(line, "\t")[[1]]
        if (length(parts) == 2) ao[[trimws(parts[1])]] <- trimws(parts[2])
      }
      ao
    }, error = function(e) { message("Warning: Could not parse AnalysisOverview: ", e$message); NULL })
  } else NULL

  # Enrich library_info with AnalysisOverview precursor counts (more reliable than log grep)
  if (!is.null(analysis_overview)) {
    ao_precursors <- analysis_overview[["Avg Profiles - Precursor"]]
    if (!is.null(ao_precursors)) {
      # Format: "66,071 of 94,202" — total (after "of") is library precursor count
      total_match <- regmatches(ao_precursors, regexpr("of ([0-9,]+)", ao_precursors))
      if (length(total_match) > 0) {
        total_val <- as.integer(gsub("[^0-9]", "", sub("^of ", "", total_match)))
        if (!is.na(total_val)) {
          if (is.null(library_info)) library_info <- list()
          if (is.na(library_info$library_precursors) || is.null(library_info$library_precursors))
            library_info$library_precursors <- total_val
        }
      }
      # Avg (before "of") is identified precursors
      avg_match <- regmatches(ao_precursors, regexpr("^[0-9,]+", ao_precursors))
      if (length(avg_match) > 0) {
        avg_val <- as.integer(gsub("[^0-9]", "", avg_match))
        if (!is.na(avg_val)) {
          if (is.null(library_info)) library_info <- list()
          library_info$avg_precursors_identified <- avg_val
        }
      }
    }
    # Also get protein group count
    ao_pg <- analysis_overview[["Avg Profiles - Protein Groups"]]
    if (!is.null(ao_pg)) {
      pg_val <- as.integer(gsub("[^0-9]", "", ao_pg))
      if (!is.na(pg_val) && !is.null(library_info) &&
          (is.na(library_info$library_proteins) || is.null(library_info$library_proteins)))
        library_info$library_proteins <- pg_val
    }
  }

  # 7. Typed search settings from ExperimentSetupOverview
  search_settings <- if (!is.null(setup_overview)) {
    tryCatch(parse_spectronaut_search_settings(setup_overview),
             error = function(e) { message("Warning: Could not parse search settings: ", e$message); NULL })
  } else NULL

  # 8. Extract n_ratios from Candidates (if parsed)
  n_ratios <- NULL
  if (has_de && manifest$candidates) {
    # Check the first comparison's data for n_ratios column
    first_comp <- de_stats_df[[1]]
    if ("n_ratios" %in% names(first_comp)) {
      n_ratios <- data.frame(
        protein_id = first_comp$protein_id,
        n_ratios   = first_comp$n_ratios,
        stringsAsFactors = FALSE
      )
    }
  }

  # --- Peptide counts ---
  # Spectronaut pivot has per-sample peptide counts: [N] Sample.PG.NrOfStrippedSequencesIdentified
  # Aggregate by taking the mean across samples (round to integer)
  n_pep <- NULL
  all_npep_cols <- grep("NrOfStrippedSequences", names(df), value = TRUE)
  if (length(all_npep_cols) > 1) {
    npep_mat <- as.data.frame(df[, all_npep_cols, with = FALSE])
    n_pep <- data.frame(
      protein_id = protein_ids,
      n_peptides = as.integer(round(rowMeans(npep_mat, na.rm = TRUE))),
      stringsAsFactors = FALSE
    )
  } else if (!is.na(npep_col)) {
    n_pep <- data.frame(
      protein_id = protein_ids,
      n_peptides = as.integer(df[[npep_col]]),
      stringsAsFactors = FALSE
    )
  }

  # --- Build settings from parsed components ---
  spec_version <- if (!is.null(library_info) && !is.null(library_info$spectronaut_version) &&
                      library_info$spectronaut_version != "Unknown") {
    library_info$spectronaut_version
  } else if (!is.null(setup_overview) && !is.null(setup_overview[["spectronaut_version"]])) {
    setup_overview[["spectronaut_version"]]
  } else if (!is.null(setup_overview) && !is.null(setup_overview[["Software Version"]])) {
    setup_overview[["Software Version"]]
  } else detect_spectronaut_version(df)

  settings <- list(
    software           = "Spectronaut",
    version            = spec_version,
    library_type       = "directDIA+ (auto-generated)",
    library_precursors = if (!is.null(library_info)) library_info$library_precursors else NA_integer_,
    library_proteins   = if (!is.null(library_info)) library_info$library_proteins else NA_integer_,
    normalization      = {
      # Try: search_settings (ExperimentSetupOverview) → library_info (AnalysisLog) → setup_overview → fallback
      norm_val <- coalesce_setting(
        if (!is.null(search_settings)) search_settings$normalization else NULL,
        coalesce_setting(
          if (!is.null(library_info)) library_info$normalization_strategy else NULL,
          if (!is.null(setup_overview) && !is.null(setup_overview[["Normalization"]])) setup_overview[["Normalization"]]
          else NULL
        )
      )
      # Build descriptive string from log fields
      if (!is.null(norm_val) && !identical(norm_val, "Not available")) {
        cross_run <- if (!is.null(library_info) && !is.null(library_info$cross_run_normalization)) {
          paste0(" (Cross-run: ", library_info$cross_run_normalization, ")")
        } else ""
        paste0(norm_val, cross_run)
      } else {
        "Local regression normalization (Spectronaut default)"
      }
    },
    topn_protein       = coalesce_setting(
      if (!is.null(setup_overview) && !is.null(setup_overview$topn_protein_max))
        paste0(setup_overview$topn_protein_max, " (TopN enabled)")
      else if (!is.null(search_settings) && !identical(search_settings$topn_protein, "Not available"))
        search_settings$topn_protein
      else NULL,
      "3 (default)"
    ),
    topn_peptide       = coalesce_setting(
      if (!is.null(setup_overview) && !is.null(setup_overview$topn_peptide_max))
        paste0(setup_overview$topn_peptide_max, " (TopN enabled)")
      else if (!is.null(search_settings) && !identical(search_settings$topn_peptide, "Not available"))
        search_settings$topn_peptide
      else NULL,
      "3 (default)"
    ),
    rollup_method      = "MaxLFQ on TopN peptides",
    de_engine          = if (has_de) {
      coalesce_setting(
        if (!is.null(search_settings)) search_settings$de_test else NULL,
        "Spectronaut internal t-test")
    } else "None (quantities only)",
    use_all_ms_levels  = coalesce_setting(
      if (!is.null(search_settings)) search_settings$use_all_ms_levels else NULL, "Not available"),
    missed_cleavages   = coalesce_setting(
      if (!is.null(search_settings)) search_settings$missed_cleavages else NULL, "Not available"),
    interference_corr  = coalesce_setting(
      if (!is.null(search_settings)) search_settings$interference_corr else NULL, "Not available"),
    de_significance    = "0.05 (applied post-hoc to Spectronaut q-values)",
    identification_fdr = coalesce_setting(
      if (!is.null(search_settings)) search_settings$protein_qvalue else NULL,
      "Not available in export"
    ),
    lfc_threshold      = coalesce_setting(
      if (!is.null(search_settings)) search_settings$lfc_filter else NULL,
      "N/A"
    ),
    imputation         = coalesce_setting(
      if (!is.null(search_settings)) search_settings$imputation else NULL, "Not available"),
    min_peptides       = "N/A (TopN controls peptide selection)",
    enzyme             = coalesce_setting(
      if (!is.null(search_settings)) search_settings$enzyme else NULL, "Not available"),
    fasta_file         = coalesce_setting(
      if (!is.null(setup_overview) && !is.null(setup_overview$fasta_databases)) setup_overview$fasta_databases
      else NULL, "Not available"),
    fasta_entries      = coalesce_setting(
      if (!is.null(setup_overview) && !is.null(setup_overview$fasta_total_entries)) setup_overview$fasta_total_entries
      else NULL, "Not available"),
    n_precursors       = coalesce_setting(
      if (!is.null(library_info) && !is.na(library_info$unique_precursors))
        format(library_info$unique_precursors, big.mark = ",") else NULL,
      "Not available"
    ),
    n_proteins_total   = as.character(length(protein_ids)),
    n_samples          = as.character(ncol(quant_mat)),
    unique_precursors  = if (!is.null(library_info) && !is.na(library_info$unique_precursors)) {
      format(library_info$unique_precursors, big.mark = ",")
    } else NULL,
    unique_proteins    = if (!is.null(library_info) && !is.na(library_info$unique_proteins)) {
      format(library_info$unique_proteins, big.mark = ",")
    } else NULL
  )

  contrasts <- if (!is.null(de_stats_df)) names(de_stats_df) else character(0)

  list(
    source          = "spectronaut",
    de_stats        = de_stats_df,
    intensities     = quant_mat,
    protein_ids     = protein_ids,
    n_peptides      = n_pep,
    n_ratios        = n_ratios,
    missing_pct     = setNames(apply(quant_mat, 1, function(x) mean(is.na(x) | x == 0)), protein_ids),
    settings        = settings,
    contrasts       = contrasts,
    sample_map      = sample_map,
    run_qc          = run_qc,
    library_info    = library_info,
    search_settings = search_settings,
    manifest        = manifest
  )
}

#' Parse Spectronaut protein group report (single TSV fallback for Mode B)
#' Accepts quantities TSV + optional DE candidates TSV (legacy 2-file mode)
parse_spectronaut <- function(file_path, de_file_path = NULL) {
  df <- data.table::fread(file_path, data.table = TRUE)

  protein_col <- grep("ProteinGroup|ProteinAccession|UniProtIds", names(df), value = TRUE)[1]
  gene_col    <- grep("^PG\\.Genes$|^Gene$|Genes", names(df), value = TRUE)[1]
  quant_cols  <- grep("PG\\.Quantity$", names(df), value = TRUE)
  npep_col    <- grep("NrOfStrippedSequences", names(df), value = TRUE)[1]

  # The quantities file may also contain DE columns (single-file export)
  logfc_col   <- grep("Log2Ratio|log2FC", names(df), value = TRUE)[1]
  qval_col    <- grep("Qvalue|Q\\.Value|q\\.value|adj", names(df), ignore.case = TRUE, value = TRUE)[1]
  pval_col    <- grep("^.*Pvalue$|^.*p\\.value$", names(df), ignore.case = TRUE, value = TRUE)[1]

  # Validate required columns
  if (is.na(protein_col)) stop("Spectronaut export missing protein column (PG.ProteinGroups)")
  if (length(quant_cols) == 0 && is.na(logfc_col)) stop("Spectronaut export has no quantities or DE stats")

  protein_ids <- normalize_protein_id(df[[protein_col]])

  # Build DE stats from quantities file (single-file mode) or separate candidates file
  de_stats_df <- NULL
  has_de <- !is.na(logfc_col) && !is.na(pval_col)

  if (!is.null(de_file_path)) {
    de_stats_df <- parse_spectronaut_candidates(de_file_path)
    has_de <- TRUE
  } else if (has_de) {
    # Single-file mode: DE columns in the same file
    de_stats_df <- list(comparison = data.frame(
      protein_id = protein_ids,
      gene       = if (!is.na(gene_col)) as.character(df[[gene_col]]) else protein_ids,
      logFC      = as.numeric(df[[logfc_col]]),
      P.Value    = as.numeric(if (!is.na(pval_col)) df[[pval_col]] else NA),
      adj.P.Val  = as.numeric(if (!is.na(qval_col)) df[[qval_col]] else NA),
      stringsAsFactors = FALSE
    ))
  }

  n_pep <- NULL
  if (!is.na(npep_col)) {
    n_pep <- data.frame(
      protein_id = protein_ids,
      n_peptides = as.integer(df[[npep_col]]),
      stringsAsFactors = FALSE
    )
  }

  quant_mat <- if (length(quant_cols) > 0) as.data.frame(df[, quant_cols, with = FALSE]) else NULL
  if (!is.null(quant_mat)) {
    # Strip .PG.Quantity suffix and [N] prefix (Spectronaut column numbering)
    clean_names <- sub("\\.PG\\.Quantity$", "", quant_cols)
    clean_names <- sub("^\\[\\d+\\]\\s*", "", clean_names)
    colnames(quant_mat) <- clean_names
  }

  contrasts <- if (!is.null(de_stats_df)) names(de_stats_df) else character(0)

  list(
    source      = "spectronaut",
    de_stats    = de_stats_df,
    intensities = quant_mat,
    protein_ids = protein_ids,
    n_peptides  = n_pep,
    missing_pct = if (!is.null(quant_mat)) {
      setNames(apply(quant_mat, 1, function(x) mean(is.na(x) | x == 0)), protein_ids)
    } else setNames(rep(NA_real_, length(protein_ids)), protein_ids),
    settings    = list(
      software         = "Spectronaut",
      version          = detect_spectronaut_version(df),
      normalization    = "Local regression (Spectronaut default)",
      rollup_method    = "Protein group quantity (PG.Quantity)",
      de_engine        = if (has_de) "Spectronaut internal statistics" else "None (quantities only)",
      de_significance  = "0.05 (applied post-hoc)",
      identification_fdr = "Not available in export",
      lfc_threshold    = "Not available in export",
      min_peptides     = "Not available in export",
      n_proteins_total = as.character(length(protein_ids)),
      n_samples        = if (!is.null(quant_mat)) as.character(ncol(quant_mat)) else "unknown"
    ),
    contrasts = contrasts
  )
}

#' Parse FragPipe-Analyst DE export
parse_fragpipe_analyst <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  sep <- if (ext == "tsv") "\t" else ","
  df <- data.table::fread(file_path, sep = sep, data.table = TRUE)

  # Validate required columns
  required <- c("Protein", "log2FC", "pvalue", "adj.pvalue")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) stop(paste("FragPipe-Analyst export missing columns:", paste(missing, collapse = ", ")))

  quant_cols <- detect_fp_intensity_cols(df)
  protein_ids <- normalize_protein_id(df$Protein)

  de_stats_df <- data.frame(
    protein_id = protein_ids,
    gene       = if ("Gene" %in% names(df)) as.character(df$Gene) else protein_ids,
    logFC      = as.numeric(df$log2FC),
    P.Value    = as.numeric(df$pvalue),
    adj.P.Val  = as.numeric(df$adj.pvalue),
    stringsAsFactors = FALSE
  )

  quant_mat <- if (length(quant_cols) > 0) as.data.frame(df[, quant_cols, with = FALSE]) else NULL

  list(
    source      = "fragpipe_analyst",
    de_stats    = list(comparison = de_stats_df),
    intensities = quant_mat,
    protein_ids = protein_ids,
    n_peptides  = NULL,
    missing_pct = if (!is.null(quant_mat)) {
      setNames(apply(quant_mat, 1, function(x) mean(is.na(x) | x == 0)), protein_ids)
    } else setNames(rep(NA_real_, length(protein_ids)), protein_ids),
    settings    = list(
      software         = "FragPipe + FragPipe-Analyst",
      fp_version       = detect_fp_version(df),
      normalization    = "IonQuant MaxLFQ (optional VSN in FragPipe-Analyst)",
      rollup_method    = "MaxLFQ consensus ratio protein rollup",
      de_engine        = "limma (via FragPipe-Analyst)",
      imputation       = "Perseus-style (FragPipe-Analyst default)",
      fdr_method       = "BH (FragPipe-Analyst default)",
      de_significance  = "Not available in export",
      identification_fdr = "0.01 (FragPipe default)",
      lfc_threshold    = "Not available in export",
      n_proteins_total = as.character(length(protein_ids)),
      n_samples        = if (!is.null(quant_mat)) as.character(ncol(quant_mat)) else "unknown"
    ),
    contrasts = "comparison"
  )
}

#' Parse FragPipe combined_protein.tsv (raw output, no DE stats)
parse_fragpipe_combined_protein <- function(file_path) {
  df <- data.table::fread(file_path, data.table = TRUE)

  if (!"Protein" %in% names(df)) stop("combined_protein.tsv missing 'Protein' column")

  maxlfq_cols <- grep("MaxLFQ Intensity$", names(df), value = TRUE)
  int_cols    <- grep(" Intensity$", names(df), value = TRUE)
  quant_cols  <- if (length(maxlfq_cols) > 0) maxlfq_cols else int_cols
  quant_cols  <- quant_cols[!grepl("^(Total|Combined)", quant_cols)]

  if (length(quant_cols) == 0) stop("combined_protein.tsv has no intensity columns")

  protein_ids <- normalize_protein_id(df$Protein)
  quant_mat   <- as.data.frame(df[, quant_cols, with = FALSE])

  n_pep <- NULL
  if ("Unique Peptides" %in% names(df)) {
    n_pep <- data.frame(
      protein_id = protein_ids,
      n_peptides = as.integer(df[["Unique Peptides"]]),
      stringsAsFactors = FALSE
    )
  }

  list(
    source      = "fragpipe_raw",
    de_stats    = NULL,
    intensities = quant_mat,
    protein_ids = protein_ids,
    n_peptides  = n_pep,
    missing_pct = setNames(
      apply(quant_mat, 1, function(x) mean(is.na(x) | x == 0)),
      protein_ids
    ),
    settings    = list(
      software        = "FragPipe (raw output)",
      normalization   = "IonQuant MaxLFQ",
      rollup_method   = "MaxLFQ consensus ratio protein rollup",
      de_engine       = "None (run FragPipe-Analyst for DE stats)",
      de_significance = "N/A",
      identification_fdr = "0.01 (FragPipe default)",
      n_proteins_total = length(unique(protein_ids)),
      n_samples        = length(quant_cols)
    ),
    contrasts = NULL
  )
}

# --- Sample Matching ---

#' Match samples between two runs by normalized filename
match_samples <- function(names_a, names_b, source_b = "delimp",
                          meta_a = NULL, meta_b = NULL) {
  strip_ext <- function(x) {
    x <- tools::file_path_sans_ext(basename(x))
    # Strip common suffixes
    x <- gsub("\\.(raw|d|mzML|wiff)$", "", x, ignore.case = TRUE)
    # Strip trailing dots (Spectronaut appends dot when label ends in digit)
    x <- gsub("\\.$", "", x)
    tolower(trimws(x))
  }

  strip_fp <- function(x) {
    # "experiment_biorep MaxLFQ Intensity"  ->  "biorep"
    x <- sub(" (MaxLFQ )?Intensity$", "", x)
    x <- sub("^[^_]+_", "", x)
    tolower(trimws(x))
  }

  stripped_a <- strip_ext(names_a)
  stripped_b <- if (source_b %in% c("fragpipe_analyst", "fragpipe_raw")) {
    strip_fp(names_b)
  } else {
    strip_ext(names_b)
  }

  matched_idx <- match(stripped_a, stripped_b)

  # Fallback: suffix-based matching for Spectronaut short labels
  # E.g., DE-LIMP "exp12102021_2_ad-1" should match Spectronaut "ad-1"
  if (any(is.na(matched_idx)) && source_b == "spectronaut") {
    for (i in which(is.na(matched_idx))) {
      # Check if any name_b is a suffix of name_a (after _ separator)
      suffix_matches <- which(endsWith(stripped_a[i], stripped_b) &
                                (nchar(stripped_a[i]) == nchar(stripped_b) |
                                 substr(stripped_a[i],
                                        nchar(stripped_a[i]) - nchar(stripped_b),
                                        nchar(stripped_a[i]) - nchar(stripped_b)) == "_"))
      if (length(suffix_matches) == 1) {
        matched_idx[i] <- suffix_matches
      }
    }
  }

  # Fallback: group-based matching when filenames differ (e.g., same samples on different instruments)
  # Match by group assignment + within-group replicate order
  if (any(is.na(matched_idx)) && !is.null(meta_a) && !is.null(meta_b)) {
    tryCatch({
      # Build group lookup for each run (File.Name -> Group, case-insensitive group matching)
      grp_a <- setNames(tolower(trimws(meta_a$Group)), meta_a$File.Name)
      grp_b <- setNames(tolower(trimws(meta_b$Group)), meta_b$File.Name)
      # Only proceed if group structures match
      groups_a <- sort(unique(grp_a[grp_a != ""]))
      groups_b <- sort(unique(grp_b[grp_b != ""]))
      if (length(groups_a) == length(groups_b) && all(groups_a == groups_b)) {
        # Within each group, match by replicate order (1st rep -> 1st rep, etc.)
        for (grp in groups_a) {
          samples_a_in_grp <- names_a[which(grp_a[names_a] == grp)]
          samples_b_in_grp <- names_b[which(grp_b[names_b] == grp)]
          if (length(samples_a_in_grp) == length(samples_b_in_grp)) {
            for (j in seq_along(samples_a_in_grp)) {
              idx_in_a <- which(names_a == samples_a_in_grp[j])
              if (is.na(matched_idx[idx_in_a])) {
                idx_in_b <- which(names_b == samples_b_in_grp[j])
                matched_idx[idx_in_a] <- idx_in_b
              }
            }
          }
        }
        message("[DE-LIMP] match_samples: group-based fallback matched ",
                sum(!is.na(matched_idx)), "/", length(names_a), " samples")
      }
    }, error = function(e) {
      message("[DE-LIMP] match_samples: group-based fallback error: ", e$message)
    })
  }

  data.frame(
    run_a    = names_a,
    run_b    = ifelse(is.na(matched_idx), NA_character_, names_b[matched_idx]),
    status   = ifelse(is.na(matched_idx), "unresolved", "matched"),
    match_method = {
      filename_idx <- match(stripped_a, stripped_b)
      ifelse(is.na(matched_idx), NA_character_,
        ifelse(!is.na(filename_idx) & filename_idx == matched_idx, "filename", "group"))
    },
    stringsAsFactors = FALSE
  )
}

# --- Analysis Functions ---

#' Merge DIA-NN log-derived fields into a run's settings list
#' Log values take precedence over "unknown" / NULL session values
merge_log_into_settings <- function(settings, log_parsed) {
  if (is.null(log_parsed) || !isTRUE(log_parsed$success)) return(settings)
  p <- log_parsed$params

  # Helper: use log value if setting is missing/unknown
  prefer <- function(current, new_val) {
    if (is.null(new_val)) return(current)
    if (is.null(current) || current == "unknown" || current == "auto") {
      as.character(new_val)
    } else current
  }

  settings$diann_version      <- prefer(settings$diann_version, log_parsed$version)
  settings$mass_acc_ms2       <- prefer(settings$mass_acc_ms2, p$mass_acc)
  settings$mass_acc_ms1       <- prefer(settings$mass_acc_ms1, p$mass_acc_ms1)
  settings$scan_window        <- prefer(settings$scan_window, p$scan_window)
  settings$enzyme             <- prefer(settings$enzyme, p$enzyme)
  settings$search_mz_min      <- prefer(settings$search_mz_min, p$min_pr_mz)
  settings$search_mz_max      <- prefer(settings$search_mz_max, p$max_pr_mz)
  settings$fasta_file         <- prefer(settings$fasta_file,
    if (length(log_parsed$fasta_files) > 0) paste(basename(log_parsed$fasta_files), collapse = ", ") else NULL)

  # Log-only fields (not in session)
  settings$pg_level           <- prefer(settings$pg_level, p$pg_level)
  settings$proteoforms        <- prefer(settings$proteoforms,
    if (isTRUE(p$proteoforms)) "yes" else "no")
  settings$reanalyse          <- prefer(settings$reanalyse,
    if (isTRUE(p$mbr)) "yes" else "no")
  settings$n_precursors_lib   <- prefer(settings$n_precursors_lib,
    if (!is.null(log_parsed$n_precursors_library)) format(log_parsed$n_precursors_library, big.mark = ",") else NULL)
  settings$library_source     <- prefer(settings$library_source,
    if (!is.null(log_parsed$lib_path)) basename(log_parsed$lib_path)
    else if (isTRUE(log_parsed$search_mode == "libfree")) "(in silico from FASTA)" else NULL)
  settings$pipeline_step      <- prefer(settings$pipeline_step, log_parsed$pipeline_step)
  settings$min_fr_mz          <- prefer(settings$min_fr_mz, p$min_fr_mz)
  settings$max_fr_mz          <- prefer(settings$max_fr_mz, p$max_fr_mz)
  settings$min_pep_len        <- prefer(settings$min_pep_len, p$min_pep_len)
  settings$max_pep_len        <- prefer(settings$max_pep_len, p$max_pep_len)
  settings$missed_cleavages   <- prefer(settings$missed_cleavages, p$missed_cleavages)

  settings
}

#' Build settings comparison table
build_settings_diff <- function(run_a, run_b) {
  # Unified parameter names — grouped by category
  params <- c(
    # Pipeline
    "software", "version", "delimp_version", "fp_version",
    "diann_version", "limpa_version", "limma_version",
    # DE settings
    "de_engine", "normalization", "rollup_method",
    "de_significance", "identification_fdr", "lfc_threshold", "min_peptides",
    "imputation", "fdr_method", "covariates", "contrast_string",
    # Search parameters
    "enzyme", "mass_acc_ms2", "mass_acc_ms1", "scan_window",
    "search_mz_min", "search_mz_max", "library",
    "fasta_file", "fasta_entries",
    # DIA-NN log-derived
    "pg_level", "proteoforms", "reanalyse",
    "library_source", "n_precursors_lib", "pipeline_step",
    "min_fr_mz", "max_fr_mz", "min_pep_len", "max_pep_len",
    "min_pr_charge", "max_pr_charge",
    "missed_cleavages", "var_mods", "fixed_mods", "max_var_mods",
    "mbr", "search_mode",
    # Data stats
    "n_precursors", "n_proteins_total", "n_samples",
    "precursor_mz_range", "charge_range",
    # Dynamic range
    "dynamic_range_orders", "dynamic_range_log2_min", "dynamic_range_log2_max",
    "dynamic_range_median", "dynamic_range_iqr", "dynamic_range_pct_missing"
  )

  # Nice labels
  labels <- c(
    # Pipeline
    "Software", "Version", "DE-LIMP Version", "FragPipe Version",
    "DIA-NN Version", "limpa Version", "limma Version",
    # DE settings
    "DE Engine", "Normalization", "Protein Rollup",
    "DE Significance (adj.P.Val cutoff)", "Identification FDR (protein group)", "logFC Threshold", "Min Peptides",
    "Imputation", "FDR Method", "Covariates", "Contrast Formula",
    # Search parameters
    "Enzyme", "MS2 Mass Accuracy", "MS1 Mass Accuracy", "Scan Window",
    "Search Min m/z", "Search Max m/z", "Spectral Library",
    "FASTA File", "FASTA Entries",
    # DIA-NN log-derived
    "Protein Grouping (pg-level)", "--proteoforms", "--reanalyse",
    "Library Source", "Precursors in Library", "Pipeline Step",
    "Min Fragment m/z", "Max Fragment m/z", "Min Peptide Length", "Max Peptide Length",
    "Min Precursor Charge", "Max Precursor Charge",
    "Missed Cleavages", "Variable Mods", "Fixed Mods", "Max Variable Mods",
    "Match Between Runs", "Search Mode",
    # Data stats
    "Precursors", "Protein Groups", "Samples",
    "Precursor m/z Range", "Charge Range",
    # Dynamic range
    "Dynamic Range (orders of magnitude)", "Log2 Intensity Min", "Log2 Intensity Max",
    "Median Log2 Intensity", "Log2 Intensity IQR", "% Missing Values"
  )

  # Section boundaries (insert header before these indices)
  section_starts <- list(
    c(1, "--- Pipeline ---"),
    c(which(params == "de_engine"), "--- DE Analysis ---"),
    c(which(params == "enzyme"), "--- Search Parameters ---"),
    c(which(params == "pg_level"), "--- DIA-NN Search (from log) ---"),
    c(which(params == "n_precursors"), "--- Data Statistics ---"),
    c(which(params == "dynamic_range_orders"), "--- Dynamic Range (precursor-level) ---")
  )

  rows <- lapply(seq_along(params), function(i) {
    val_a <- run_a$settings[[params[i]]] %||% NA_character_
    val_b <- run_b$settings[[params[i]]] %||% NA_character_
    if (is.na(val_a) && is.na(val_b)) return(NULL)
    data.frame(
      Parameter = labels[i],
      Run_A = as.character(val_a %||% "N/A"),
      Run_B = as.character(val_b %||% "N/A"),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(result) || nrow(result) == 0) {
    return(data.frame(Parameter = "No settings to compare",
                      Run_A = "", Run_B = "", match = "match",
                      stringsAsFactors = FALSE))
  }

  result$match <- ifelse(
    result$Run_A == "N/A" | result$Run_B == "N/A", "unknown",
    ifelse(result$Run_A == result$Run_B, "match", "differs")
  )
  result$note <- ""

  # --- Spectronaut-specific rows (Mode B) ---
  if (identical(run_b$source, "spectronaut")) {
    ss <- run_b$search_settings  # may be NULL
    li <- run_b$library_info     # may be NULL
    use_all_ms <- if (!is.null(ss)) trimws(coalesce_setting(ss$use_all_ms_levels, "Not available")) else "Not available"

    spec_rows <- data.frame(
      Parameter = c(
        "--- Spectronaut-Specific ---",
        "Library Type",
        "Library Precursors",
        "Library Protein Groups",
        "Protein TopN (max peptides for quant)",
        "Peptide TopN (max precursors per peptide)",
        "Protein Quantity Calculation",
        "Statistical Test",
        "Interference Correction",
        "Use All MS-Level Quantities (Quant3)"
      ),
      Run_A = c(
        "",
        coalesce_setting(run_a$settings$library_source, "DIA-NN empirical library"),
        format_or_na(run_a$settings$n_precursors_lib),
        format_or_na(run_a$settings$n_proteins_total),
        "ALL (DPC-Quant uses all detected precursors)",
        "ALL",
        "DPC-Quant: empirical Bayes precursor aggregation",
        "limma moderated t-test (empirical Bayes variance shrinkage)",
        "DIA-NN: integrated signal correction",
        "N/A -- DE-LIMP uses single protein intensity per sample"
      ),
      Run_B = c(
        "",
        coalesce_setting(run_b$settings$library_type, "directDIA+ (auto-generated)"),
        format_or_na(if (!is.null(li)) li$library_precursors else NA),
        format_or_na(if (!is.null(li)) li$library_proteins else NA),
        coalesce_setting(run_b$settings$topn_protein, "3 (default)"),
        coalesce_setting(run_b$settings$topn_peptide, "3 (default)"),
        "MaxLFQ on TopN peptides",
        coalesce_setting(run_b$settings$de_engine, "Unpaired t-test (Welch)"),
        coalesce_setting(if (!is.null(ss)) ss$interference_corr else NULL, "Not available"),
        use_all_ms
      ),
      match = c(
        "match",
        "structural_difference", "differs", "differs",
        "structural_difference", "structural_difference",
        "structural_difference", "structural_difference",
        "unknown",
        if (identical(use_all_ms, "True")) "severe" else "structural_difference"
      ),
      note = c(
        "",
        "Spectronaut builds its own in-experiment library from the data",
        "Larger library = broader protein coverage",
        "",
        "Critical: Spectronaut TopN limits quant to N peptides; DPC-Quant is unbounded",
        "",
        "These algorithms weight evidence differently",
        "Moderated t-test is substantially more conservative",
        "",
        if (identical(use_all_ms, "True"))
          "ENABLED -- Quant3 doubles effective observation count; p-values not comparable to limma"
        else "Disabled -- standard single-level quantification"
      ),
      stringsAsFactors = FALSE
    )

    # Ensure result has a note column before rbind
    result <- rbind(result, spec_rows)
  }

  result
}

#' Classify proteins into universe tiers
classify_protein_universe <- function(run_a, run_b) {
  ids_a <- run_a$protein_ids
  ids_b <- run_b$protein_ids
  all_ids <- union(ids_a, ids_b)

  # Helper: check if a string looks like a gene symbol (not an accession)
  is_gene_symbol <- function(x) {
    # Gene symbols: short (1-15 chars), mostly letters, may have digits at end (e.g. TP53, BRCA1)
    # Accessions: have semicolons, pipes, long digit strings, underscores (e.g. A0A075B6K5;P80748)
    !is.na(x) & nzchar(x) &
      nchar(x) <= 15 &
      !grepl(";", x) &
      !grepl("\\|", x) &
      !grepl("^[A-Z][0-9][A-Z0-9]{3}[0-9]$", x) &  # UniProt accession pattern
      !grepl("^[A-Z0-9]{6,}", x) &  # Long accession-like strings
      grepl("^[A-Za-z]", x)  # Starts with letter
  }

  # Try to extract gene symbols from DE stats
  gene_map <- character(0)
  for (run in list(run_a, run_b)) {
    if (!is.null(run$de_stats)) {
      for (ct in names(run$de_stats)) {
        tbl <- run$de_stats[[ct]]
        # Prefer "Gene" (from bitr mapping) over "Genes" (raw DIA-NN, often accessions)
        gene_col <- NULL
        if ("Gene" %in% names(tbl)) gene_col <- "Gene"
        else if ("gene" %in% names(tbl)) gene_col <- "gene"
        else if ("Genes" %in% names(tbl)) gene_col <- "Genes"
        if (is.null(gene_col)) next

        candidates <- setNames(as.character(tbl[[gene_col]]), tbl$protein_id)
        candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
        # Filter to actual gene symbols, not accessions
        valid <- is_gene_symbol(candidates)
        candidates <- candidates[valid]
        gene_map <- c(gene_map, candidates[!names(candidates) %in% names(gene_map)])
      }
    }
  }

  # Fallback: extract gene from sp|ACC|GENE_SPECIES format
  if (length(gene_map) == 0) {
    pipe_ids <- grep("\\|", all_ids, value = TRUE)
    if (length(pipe_ids) > 0) {
      parts <- strsplit(pipe_ids, "\\|")
      genes <- sapply(parts, function(p) {
        if (length(p) >= 3) sub("_[A-Z]+$", "", p[3]) else NA_character_
      })
      gene_map <- setNames(genes, normalize_protein_id(pipe_ids))
    }
  }

  df <- data.frame(
    protein_id = all_ids,
    in_a       = all_ids %in% ids_a,
    in_b       = all_ids %in% ids_b,
    tier       = ifelse(
      all_ids %in% ids_a & all_ids %in% ids_b, "shared",
      ifelse(all_ids %in% ids_a, "a_only", "b_only")
    ),
    stringsAsFactors = FALSE
  )

  # Add gene symbols where available
  df$gene_symbol <- gene_map[df$protein_id]
  df$gene_symbol[is.na(df$gene_symbol)] <- ""
  df
}

#' Compute quantification comparison for shared proteins
compute_quant_comparison <- function(run_a, run_b, universe, sample_map) {
  shared_ids <- universe$protein_id[universe$tier == "shared"]

  # Get indices in each run
  idx_a <- match(shared_ids, run_a$protein_ids)
  idx_b <- match(shared_ids, run_b$protein_ids)

  # Remove any with NA indices (shouldn't happen but safety)
  valid <- !is.na(idx_a) & !is.na(idx_b)
  shared_ids <- shared_ids[valid]
  idx_a <- idx_a[valid]
  idx_b <- idx_b[valid]

  if (length(shared_ids) == 0) return(NULL)

  # Mean log2 intensity per protein
  mat_a <- as.matrix(run_a$intensities[idx_a, , drop = FALSE])
  mat_b <- as.matrix(run_b$intensities[idx_b, , drop = FALSE])

  # Log2 transform if not already (check if max > 30 suggests non-log scale)
  if (max(mat_a, na.rm = TRUE) > 30) mat_a <- log2_safe(mat_a)
  if (max(mat_b, na.rm = TRUE) > 30) mat_b <- log2_safe(mat_b)

  mean_a <- rowMeans(mat_a, na.rm = TRUE)
  mean_b <- rowMeans(mat_b, na.rm = TRUE)

  # Per-protein offset
  offsets <- mean_a - mean_b
  median_offset <- median(offsets, na.rm = TRUE)

  # Per-sample correlation matrix (matched samples)
  matched <- sample_map[sample_map$status == "matched", ]
  if (nrow(matched) > 0) {
    cor_a_cols <- match(matched$run_a, colnames(mat_a))
    cor_b_cols <- match(matched$run_b, colnames(mat_b))
    valid_pairs <- !is.na(cor_a_cols) & !is.na(cor_b_cols)

    if (sum(valid_pairs) > 0) {
      cor_a_cols <- cor_a_cols[valid_pairs]
      cor_b_cols <- cor_b_cols[valid_pairs]
      sample_labels <- matched$run_a[valid_pairs]

      cor_vals <- sapply(seq_along(cor_a_cols), function(i) {
        a_vec <- mat_a[, cor_a_cols[i]]
        b_vec <- mat_b[, cor_b_cols[i]]
        complete <- complete.cases(a_vec, b_vec)
        if (sum(complete) < 3) return(NA_real_)
        cor(a_vec[complete], b_vec[complete])
      })
      names(cor_vals) <- sample_labels
    } else {
      cor_vals <- numeric(0)
    }
  } else {
    cor_vals <- numeric(0)
  }

  list(
    scatter_data = data.frame(
      protein_id = shared_ids,
      mean_a     = mean_a,
      mean_b     = mean_b,
      stringsAsFactors = FALSE
    ),
    offsets       = offsets,
    median_offset = median_offset,
    sd_offset     = sd(offsets, na.rm = TRUE),
    cor_per_sample = cor_vals,
    min_cor       = if (length(cor_vals) > 0) min(cor_vals, na.rm = TRUE) else NA,
    max_cor       = if (length(cor_vals) > 0) max(cor_vals, na.rm = TRUE) else NA,
    mat_a         = mat_a,
    mat_b         = mat_b
  )
}

#' Compute DE concordance (3x3: Up/Down/NS for each run)
compute_de_concordance <- function(run_a, run_b, universe, contrast_a, contrast_b,
                                   source_b, global_offset = 0,
                                   imputation_b = NULL) {
  shared_ids <- universe$protein_id[universe$tier == "shared"]

  # Get DE stats for selected contrast
  de_a <- if (is.data.frame(run_a$de_stats)) {
    run_a$de_stats
  } else if (is.list(run_a$de_stats)) {
    run_a$de_stats[[contrast_a]]
  } else return(NULL)

  de_b <- if (is.data.frame(run_b$de_stats)) {
    run_b$de_stats
  } else if (is.list(run_b$de_stats)) {
    run_b$de_stats[[contrast_b]]
  } else return(NULL)

  if (is.null(de_a) || is.null(de_b)) return(NULL)

  # Filter to shared proteins
  de_a <- de_a[de_a$protein_id %in% shared_ids, ]
  de_b <- de_b[de_b$protein_id %in% shared_ids, ]

  # Merge on protein_id
  merged <- merge(de_a, de_b, by = "protein_id", suffixes = c("_A", "_B"))

  if (nrow(merged) == 0) return(NULL)

  # Use correct column names based on source
  logfc_a <- if ("logFC_A" %in% names(merged)) merged$logFC_A else merged$logFC
  logfc_b <- if ("logFC_B" %in% names(merged)) merged$logFC_B else merged$logFC
  adjp_a  <- if ("adj.P.Val_A" %in% names(merged)) merged$adj.P.Val_A else merged$adj.P.Val
  adjp_b  <- if ("adj.P.Val_B" %in% names(merged)) merged$adj.P.Val_B else merged$adj.P.Val

  # Classify each protein: Up / Down / NS (NaN-safe)
  classify_de <- function(logfc, adjp, threshold = 0.05) {
    ifelse(!is.finite(adjp) | adjp >= threshold, "NS",
           ifelse(is.finite(logfc) & logfc > 0, "Up",
                  ifelse(is.finite(logfc) & logfc < 0, "Down", "NS")))
  }

  merged$status_a <- classify_de(logfc_a, adjp_a)
  merged$status_b <- classify_de(logfc_b, adjp_b)
  merged$logFC_A  <- logfc_a
  merged$logFC_B  <- logfc_b
  merged$adjP_A   <- adjp_a
  merged$adjP_B   <- adjp_b

  # 3x3 concordance matrix
  levels_3 <- c("Up", "Down", "NS")
  conc_3x3 <- table(
    factor(merged$status_a, levels = levels_3),
    factor(merged$status_b, levels = levels_3)
  )

  # Add n_ratios to merged data for rescue stats (Mode B)
  merged$n_ratios_B <- NA_integer_
  if (source_b == "spectronaut" && !is.null(run_b$n_ratios)) {
    nr_match <- match(merged$protein_id, run_b$n_ratios$protein_id)
    merged$n_ratios_B <- run_b$n_ratios$n_ratios[nr_match]
  }

  # Identify discordant proteins (different DE status)
  merged$concordant <- merged$status_a == merged$status_b
  discordant <- merged[!merged$concordant, ]

  # Compute "rescue" stats: proteins with 0 ratios in Spectronaut but tested in DE-LIMP
  rescue_stats <- NULL
  if (source_b == "spectronaut" && any(!is.na(merged$n_ratios_B))) {
    zero_ratio <- merged$n_ratios_B == 0 & !is.na(merged$n_ratios_B)
    n_zero_ratio <- sum(zero_ratio)
    # Of those, how many are significant in DE-LIMP?
    sig_in_a <- zero_ratio & merged$status_a != "NS"
    n_rescued_sig <- sum(sig_in_a)
    # Low ratio (<25% of max possible)
    n_group_1 <- n_group_2 <- 1L
    if (!is.null(run_b$sample_map)) {
      conds <- unique(run_b$sample_map$condition)
      if (length(conds) == 2) {
        n_group_1 <- sum(run_b$sample_map$condition == conds[1])
        n_group_2 <- sum(run_b$sample_map$condition == conds[2])
      }
    }
    max_ratios <- n_group_1 * n_group_2
    low_ratio <- !is.na(merged$n_ratios_B) & merged$n_ratios_B > 0 &
                 merged$n_ratios_B < (max_ratios * 0.25)
    n_low_ratio <- sum(low_ratio)
    # Imputation context
    imp_setting <- if (!is.null(imputation_b)) trimws(imputation_b) else "Unknown"
    imp_none <- grepl("^none$", imp_setting, ignore.case = TRUE)
    rescue_stats <- list(
      n_zero_ratio   = n_zero_ratio,
      n_rescued_sig  = n_rescued_sig,
      n_low_ratio    = n_low_ratio,
      max_ratios     = max_ratios,
      imputation     = imp_setting,
      imputation_off = imp_none
    )
  }

  # Get peptide counts and missing percentages
  if (nrow(discordant) > 0 && !is.null(run_a$n_peptides)) {
    pep_a <- run_a$n_peptides[, c("protein_id", "n_peptides"), drop = FALSE]
    names(pep_a)[2] <- "n_peptides_A"
    discordant <- merge(discordant, pep_a, by = "protein_id", all.x = TRUE)
  } else {
    discordant$n_peptides_A <- NA_integer_
  }
  if (nrow(discordant) > 0 && !is.null(run_b$n_peptides)) {
    pep_b <- run_b$n_peptides[, c("protein_id", "n_peptides"), drop = FALSE]
    names(pep_b)[2] <- "n_peptides_B"
    discordant <- merge(discordant, pep_b, by = "protein_id", all.x = TRUE)
  } else {
    discordant$n_peptides_B <- NA_integer_
  }

  # Add missing percentages
  discordant$missing_pct_A <- run_a$missing_pct[match(discordant$protein_id, names(run_a$missing_pct))]
  discordant$missing_pct_B <- run_b$missing_pct[match(discordant$protein_id, names(run_b$missing_pct))]

  # n_ratios_B already inherited from merged (added above)

  # Compute group sizes for n_ratios context
  n_group_a <- n_group_b <- 1L
  if (!is.null(run_b$sample_map)) {
    conds <- unique(run_b$sample_map$condition)
    if (length(conds) == 2) {
      n_group_a <- sum(run_b$sample_map$condition == conds[1])
      n_group_b <- sum(run_b$sample_map$condition == conds[2])
    }
  }
  use_all_ms_true <- if (!is.null(run_b$search_settings)) {
    identical(trimws(coalesce_setting(run_b$search_settings$use_all_ms_levels, "")), "True")
  } else if (!is.null(run_b$settings$use_all_ms_levels)) {
    identical(trimws(run_b$settings$use_all_ms_levels), "True")
  } else FALSE

  # Assign hypotheses (pass context via hidden columns)
  if (nrow(discordant) > 0) {
    discordant$.use_all_ms_true <- use_all_ms_true
    discordant$.n_group_a <- n_group_a
    discordant$.n_group_b <- n_group_b
    discordant$.imputation_b <- imputation_b %||% "Unknown"
    hyp_results <- lapply(seq_len(nrow(discordant)), function(i) {
      assign_hypothesis(discordant[i, ], source_b, global_offset)
    })
    discordant$hypothesis  <- sapply(hyp_results, `[[`, "hypothesis")
    discordant$confidence  <- sapply(hyp_results, `[[`, "confidence")
    discordant$hypothesis_category <- sapply(hyp_results, `[[`, "category")
  } else {
    discordant$hypothesis  <- character(0)
    discordant$confidence  <- character(0)
    discordant$hypothesis_category <- character(0)
  }

  # Sort discordant: High confidence first, then by FC difference
  if (nrow(discordant) > 0) {
    # Remove context columns
    discordant$.use_all_ms_true <- NULL
    discordant$.n_group_a <- NULL
    discordant$.n_group_b <- NULL
    discordant$.imputation_b <- NULL
    conf_order <- c("High" = 1, "Medium" = 2, "Low" = 3)
    discordant$conf_rank <- conf_order[discordant$confidence]
    discordant <- discordant[order(discordant$conf_rank,
                                   -abs(discordant$logFC_A - discordant$logFC_B)), ]
    discordant$conf_rank <- NULL
  }

  list(
    concordance_3x3  = conc_3x3,
    merged           = merged,
    discordant_table = discordant,
    n_concordant     = sum(merged$concordant),
    n_discordant     = sum(!merged$concordant),
    n_total          = nrow(merged),
    rescue_stats     = rescue_stats
  )
}

#' Assign diagnostic hypothesis for a discordant protein (8 rules + Rule 0, tool-aware)
assign_hypothesis <- function(row, source_b, global_offset = 0) {
  lfc_a <- if (is.finite(row$logFC_A)) row$logFC_A else 0
  lfc_b <- if (is.finite(row$logFC_B)) row$logFC_B else 0
  adjp_a <- if (is.finite(row$adjP_A)) row$adjP_A else 1
  adjp_b <- if (is.finite(row$adjP_B)) row$adjP_B else 1

  # Rule 0: Spectronaut had 0 computable ratios (highest priority)
  if (source_b == "spectronaut" && !is.na(row$n_ratios_B) && row$n_ratios_B == 0) {
    imp <- row$.imputation_b %||% "Unknown"
    imp_none <- grepl("^none$", imp, ignore.case = TRUE)
    if (imp_none) {
      msg <- paste0(
        "Spectronaut had 0 computable ratios for this protein (complete missingness in one condition). ",
        "Imputation was disabled in Spectronaut, so no statistical test was possible. ",
        "DE-LIMP tested this protein via limpa empirical Bayes modelling with DIA-NN's match-between-runs quantification.")
    } else if (grepl("unknown", imp, ignore.case = TRUE)) {
      msg <- paste0(
        "Spectronaut had 0 computable ratios for this protein (complete missingness in one condition). ",
        "No statistical test was possible in Spectronaut. ",
        "DE-LIMP tested this protein via limpa empirical Bayes modelling.")
    } else {
      msg <- paste0(
        "Spectronaut had 0 computable ratios despite imputation being enabled (", imp, "). ",
        "This protein had complete missingness in one condition that imputation could not rescue. ",
        "DE-LIMP tested this protein via limpa empirical Bayes modelling with DIA-NN's match-between-runs quantification.")
    }
    return(list(hypothesis = msg, confidence = "High", category = "Untestable in Spectronaut"))
  }

  logfc_diff   <- abs(lfc_a - lfc_b)
  same_dir     <- sign(lfc_a) == sign(lfc_b)
  p_diff_large <- abs(log10(adjp_a + 1e-10) - log10(adjp_b + 1e-10)) > 1
  missing_diff <- abs((row$missing_pct_A %||% 0) - (row$missing_pct_B %||% 0)) > 0.2
  peptide_diff <- !is.na(row$n_peptides_A) && !is.na(row$n_peptides_B) &&
                  abs(row$n_peptides_A - row$n_peptides_B) > 2
  sys_offset   <- abs(global_offset) > 0.2

  # Rule 1: Direction reversal
  if (!same_dir && logfc_diff > 0.5) return(list(
    hypothesis = "logFC direction reverses between runs - likely a quantification or normalization discrepancy at the precursor level",
    confidence = "High",
    category   = "Direction reversal"
  ))

  # Rule 2: Systematic normalization offset
  if (sys_offset && same_dir && logfc_diff < 0.5) {
    tool_note <- switch(source_b,
      "spectronaut"      = " (Spectronaut local regression vs DIA-NN RT-dependent + DPC-CN normalization)",
      "fragpipe_analyst" = " (IonQuant normalization vs DIA-NN RT-dependent + DPC-CN normalization)",
      "fragpipe_raw"     = " (MaxLFQ vs DPC-Quant)",
      ""
    )
    return(list(
      hypothesis = paste0("Global intensity offset of ", round(global_offset, 2),
                         " log2 units detected", tool_note,
                         ". Borderline significance may reflect normalization scale difference rather than biology."),
      confidence = "High",
      category   = "Normalization offset"
    ))
  }

  # Rule 3: Same FC, divergent p-value
  if (same_dir && logfc_diff < 0.3 && p_diff_large) {
    engine_note <- if (source_b == "fragpipe_analyst") {
      " Both use limma, so this likely reflects missing-value handling differences: FragPipe-Analyst uses Perseus-style imputation; DE-LIMP uses DPC-Quant probabilistic modelling (not imputation)."
    } else if (source_b == "spectronaut") {
      if (isTRUE(row$.use_all_ms_true)) {
        " Spectronaut used 'Use All MS-Level Quantities' (Quant3 method), which combines MS1 and MS2 observations to double effective sample size in its t-test. This is the most likely cause of the p-value divergence."
      } else {
        " Spectronaut uses an unpaired t-test; DE-LIMP uses limma empirical Bayes variance shrinkage. Borderline proteins near the FDR threshold will differ between methods."
      }
    } else ""
    return(list(
      hypothesis = paste0("Similar fold-change but divergent p-value - likely variance estimation or effective sample size difference.", engine_note),
      confidence = "High",
      category   = "Variance estimation"
    ))
  }

  # Rule 4: Missing value pattern differs
  if (missing_diff) {
    mbr_note <- if (source_b %in% c("fragpipe_analyst", "fragpipe_raw")) {
      " FragPipe uses IonQuant MBR; DE-LIMP uses DIA-NN MBR. Different implementations produce different dropout patterns."
    } else ""
    return(list(
      hypothesis = paste0("Missing value pattern differs (",
                         round((row$missing_pct_A %||% 0) * 100), "% vs ",
                         round((row$missing_pct_B %||% 0) * 100), "% missing).", mbr_note),
      confidence = "High",
      category   = "Missing values"
    ))
  }

  # Rule 8: Low ratio count in Spectronaut (Mode B only)
  if (source_b == "spectronaut" && !is.na(row$n_ratios_B %||% NA)) {
    max_ratios <- (row$.n_group_a %||% 1L) * (row$.n_group_b %||% 1L)
    if (max_ratios > 0) {
      pct_ratios <- round(100 * row$n_ratios_B / max_ratios)
      if (pct_ratios < 25) {
        return(list(
          hypothesis = sprintf(
            "Very low ratio count in Spectronaut: only %d of %d possible pairwise ratios were computable (%d%% completeness). This protein has high missingness in Spectronaut's quantification.",
            row$n_ratios_B, max_ratios, pct_ratios),
          confidence = "High",
          category   = "Low ratio count"
        ))
      } else if (pct_ratios < 50) {
        return(list(
          hypothesis = sprintf(
            "Moderate ratio count in Spectronaut (%d of %d possible, %d%%). Protein has partial missingness; Spectronaut's t-test may be underpowered relative to limma.",
            row$n_ratios_B, max_ratios, pct_ratios),
          confidence = "Medium",
          category   = "Low ratio count"
        ))
      }
    }
  }

  # Rule 5: Peptide count differs
  if (peptide_diff) {
    rollup_note <- if (source_b %in% c("fragpipe_analyst", "fragpipe_raw")) {
      paste0(" DIA-NN and FragPipe use different rollup strategies (DPC-Quant vs MaxLFQ).",
             " Proteins with few peptides are most sensitive to this.")
    } else ""
    return(list(
      hypothesis = paste0("Peptide count differs (", row$n_peptides_A, " vs ",
                         row$n_peptides_B, " peptides).", rollup_note),
      confidence = "Medium",
      category   = "Peptide count"
    ))
  }

  # Rule 6: Same direction, large FC magnitude difference
  if (same_dir && logfc_diff > 0.5) return(list(
    hypothesis = "Same direction but different magnitude - normalization scale difference between runs",
    confidence = "Medium",
    category   = "FC magnitude"
  ))

  # Rule 7: Default - borderline significance
  list(
    hypothesis = "Borderline significance in both runs - stochastic FDR threshold effect (protein near cutoff in both)",
    confidence = "Low",
    category   = "Borderline"
  )
}

# --- AI Prompt Builders ---

#' Build Gemini prompt for run comparator (tool-aware)
build_gemini_comparator_prompt <- function(comp_results, mofa_obj = NULL, instrument_meta = NULL) {
  stats    <- comp_results$summary_stats
  top_disc <- head(comp_results$de_concordance$discordant_table, 10)
  source_b <- stats$source_b

  tool_context <- switch(source_b,
    "delimp" =
      "Both runs used DE-LIMP (DIA-NN -> limpa/limma pipeline). The comparison isolates the effect of search or analysis parameter differences on the same raw data.",
    "spectronaut" = {
      run_b <- comp_results$run_b
      ss <- run_b$search_settings
      li <- run_b$library_info
      topn_val <- coalesce_setting(run_b$settings$topn_protein, "3")
      de_test <- coalesce_setting(run_b$settings$de_engine, "unpaired t-test")
      use_all_ms <- if (!is.null(ss)) trimws(coalesce_setting(ss$use_all_ms_levels, "")) else ""

      # Compute stats
      pct_limited <- if (!is.null(run_b$n_peptides)) {
        topn_int <- suppressWarnings(as.integer(gsub("[^0-9]", "", topn_val)))
        if (is.na(topn_int)) topn_int <- 3L
        round(100 * mean(run_b$n_peptides$n_peptides > topn_int, na.rm = TRUE))
      } else "unknown"
      mean_pep <- if (!is.null(run_b$n_peptides)) {
        round(mean(run_b$n_peptides$n_peptides, na.rm = TRUE), 1)
      } else "unknown"

      # Sig counts from first comparison
      first_de <- if (!is.null(run_b$de_stats)) run_b$de_stats[[1]] else NULL
      n_sig_spec <- if (!is.null(first_de)) sum(first_de$adj.P.Val <= 0.05, na.rm = TRUE) else "unknown"
      n_total_b <- if (!is.null(first_de)) nrow(first_de) else 0
      pct_sig_spec <- if (n_total_b > 0) round(100 * n_sig_spec / n_total_b) else "?"
      median_lfc_sig <- if (!is.null(first_de)) {
        sig_rows <- first_de[!is.na(first_de$adj.P.Val) & first_de$adj.P.Val <= 0.05, ]
        if (nrow(sig_rows) > 0) round(median(abs(sig_rows$logFC), na.rm = TRUE), 2) else "N/A"
      } else "unknown"

      quant3_note <- if (identical(use_all_ms, "True")) {
        paste0("\n\nQUANT3 SETTING: ",
               "Spectronaut's 'Use All MS-Level Quantities' was ON. ",
               "This combines MS1 and MS2 measurements, increasing the observation count ",
               "per protein in the t-test. This can increase statistical power but also ",
               "raises questions about observation independence (MS1 and MS2 from the same ",
               "peptide are correlated). Consider whether this affects the p-value comparison.")
      } else ""

      paste0(
        "TOOL COMPARISON CONTEXT:\n",
        "Run A: DE-LIMP -- DIA-NN search -> DPC-Quant protein rollup (all precursors, empirical Bayes weighting) -> limma moderated t-test (variance stabilization across proteins).\n",
        "Run B: Spectronaut ", coalesce_setting(run_b$settings$version, ""), " -- directDIA+ search -> MaxLFQ (TopN=", topn_val, ", uses ", topn_val, " most intense peptides per protein) -> ", de_test, ".\n",
        "\nQUANTIFICATION APPROACH DIFFERENCE:\n",
        "Spectronaut: TopN=", topn_val, " peptide selection for protein quantification. ",
        "The TopN cap applied to ~", pct_limited, "% of proteins (average ", mean_pep, " peptides detected). ",
        "DE-LIMP: DPC-Quant aggregates all detected precursors with empirical Bayes weighting. ",
        "Both approaches have trade-offs: TopN may be more robust to noisy peptides; all-precursor uses more data but may include low-quality signals.\n",
        "\nSTATISTICAL MODEL DIFFERENCE:\n",
        "Run A (limma): moderated t-test with empirical Bayes variance shrinkage — called ", stats$n_sig_a, " proteins significant. ",
        "Run B (", de_test, "): called ", n_sig_spec, " significant at q<=0.05 (", pct_sig_spec, "% of ", n_total_b, " total). ",
        "Median |logFC| among Run B significant proteins: ", median_lfc_sig, " log2. ",
        "These methods differ in how they estimate per-protein variance, which can produce different significance calls even from similar fold-changes.\n",
        "\nPRE-FILTERING DIFFERENCES:\n",
        "Run B applies a Log2 Ratio Candidate Filter (", coalesce_setting(run_b$settings$lfc_threshold, "unknown"), " log2) ",
        "BEFORE statistical testing — proteins with insufficient fold-change evidence are excluded from DE testing entirely (NaN ratios). ",
        "Run B also uses ", coalesce_setting(run_b$settings$imputation, "unknown"), " imputation. ",
        "Run A (limpa) uses DPC-Quant, which models missing values probabilistically via a detection probability curve (not imputation), and tests all quantified proteins.\n",
        "\nNORMALIZATION DIFFERENCE:\n",
        "Run A: DIA-NN RT-dependent normalization + DPC-CN cyclic loess. ",
        "Run B: ", coalesce_setting(run_b$settings$normalization, "local regression"), ". ",
        "Different normalization strategies can introduce a global intensity offset.",
        quant3_note
      )
    },
    "fragpipe_analyst" =
      paste("Run A: DE-LIMP (DIA-NN search -> DPC-Quant rollup -> limma DE).",
            "Run B: FragPipe + FragPipe-Analyst (MSFragger -> IonQuant MaxLFQ rollup -> limma DE).",
            "Key structural differences:",
            "- Search engine: DIA-NN (DIA-optimized) vs MSFragger (DDA/DIA)",
            "- Protein rollup: DPC-Quant empirical Bayes vs MaxLFQ consensus ratios",
            "- Normalization: DIA-NN RT-dependent vs IonQuant (optional VSN)",
            "- Missingness: DIA-NN MBR vs IonQuant MBR (different implementations)",
            "- Missing values: DPC-Quant probabilistic modelling (not imputation) vs Perseus-style imputation (FP-Analyst default)",
            "Both use limma for DE, so p-value calibration should be broadly comparable."),
    "fragpipe_raw" =
      "Run A: DE-LIMP (DIA-NN). Run B: raw FragPipe output (no DE stats). Only quantification differences can be assessed.",
    ""
  )

  # Format discordant table
  disc_text <- if (nrow(top_disc) > 0) {
    paste(apply(top_disc[, c("protein_id", "logFC_A", "logFC_B", "adjP_A", "adjP_B",
                             "hypothesis_category"), drop = FALSE], 1,
      function(r) paste(r, collapse = " | ")), collapse = "\n")
  } else "No discordant proteins"

  # Settings diff summary
  settings_diff <- comp_results$settings_diff
  diff_rows <- settings_diff[settings_diff$match %in% c("differs", "structural_difference", "severe"), ]
  settings_summary <- if (nrow(diff_rows) > 0) {
    paste(apply(diff_rows[, c("Parameter", "Run_A", "Run_B")], 1,
      function(r) paste0("- ", r[1], ": ", r[2], " vs ", r[3])), collapse = "\n")
  } else "All settings identical"

  # MOFA2 section (if available)
  mofa_section <- ""
  if (!is.null(mofa_obj) && requireNamespace("MOFA2", quietly = TRUE)) {
    tryCatch({
      var_exp    <- MOFA2::get_variance_explained(mofa_obj)
      r2_a       <- var_exp$r2_per_factor[[1]][, "run_a"]
      r2_b       <- var_exp$r2_per_factor[[1]][, "run_b"]
      specificity <- r2_a - r2_b
      key_factor  <- names(which.max(specificity))

      top_weights <- MOFA2::get_weights(mofa_obj, views = "run_a",
                                         factors = key_factor,
                                         as.data.frame = TRUE) |>
        dplyr::arrange(dplyr::desc(abs(value))) |>
        dplyr::slice_head(n = 5) |>
        dplyr::pull(feature) |>
        paste(collapse = ", ")

      mofa_section <- paste0(
        "\n\nMOFA2 FACTOR DECOMPOSITION:\n",
        "- ", key_factor, " is the most run-A-specific factor\n",
        "  (explains ", round(r2_a[key_factor] * 100, 1), "% variance in Run A, ",
        round(r2_b[key_factor] * 100, 1), "% in Run B)\n",
        "- Top proteins driving this factor: ", top_weights, "\n",
        "- This suggests Run B's pipeline is suppressing the signal captured by this factor."
      )
    }, error = function(e) "")
  }

  # Instrument context (brief)
  instrument_section <- ""
  if (!is.null(instrument_meta)) {
    inst_parts <- character(0)
    if (!is.null(instrument_meta$instrument_model) && nzchar(instrument_meta$instrument_model))
      inst_parts <- c(inst_parts, instrument_meta$instrument_model)
    if (!is.null(instrument_meta$lc_system) && nzchar(instrument_meta$lc_system))
      inst_parts <- c(inst_parts, instrument_meta$lc_system)
    if (!is.null(instrument_meta$lc_method) && nzchar(instrument_meta$lc_method))
      inst_parts <- c(inst_parts, instrument_meta$lc_method)
    if (!is.null(instrument_meta$evosep_spd) && !is.na(instrument_meta$evosep_spd))
      inst_parts <- c(inst_parts, paste0(instrument_meta$evosep_spd, " SPD"))
    if (!is.null(instrument_meta$evosep_gradient_min) && !is.na(instrument_meta$evosep_gradient_min))
      inst_parts <- c(inst_parts, paste0(instrument_meta$evosep_gradient_min, " min gradient"))
    if (length(inst_parts) > 0)
      instrument_section <- paste0("INSTRUMENT: ", paste(inst_parts, collapse = ", "), "\n")
  }

  paste0(
    "You are analyzing a proteomics cross-tool comparison. Be objective — neither tool is inherently superior. Evaluate based on the data.\n\n",
    instrument_section,
    "TOOL COMPARISON CONTEXT:\n", tool_context, "\n\n",
    "COMPARISON OVERVIEW:\n",
    "- Run A: ", stats$run_a_label, " (", stats$n_samples, " samples)\n",
    "- Run B: ", stats$run_b_label, "\n",
    "- Contrast: ", stats$contrast, "\n",
    "- Proteins shared: ", stats$n_shared, "\n",
    "- Concordant DE: ", stats$n_concordant, "\n",
    "- Discordant: ", stats$n_discordant, " (",
      stats$n_a_only_de, " DE in A only, ", stats$n_b_only_de, " DE in B only)\n",
    "- Global intensity offset (median): ", round(stats$global_offset, 3), " log2 units\n",
    "- Per-sample correlation range: ",
      if (!is.na(stats$min_cor)) paste0(round(stats$min_cor, 2), "-", round(stats$max_cor, 2)) else "N/A",
    "\n",
    # Library size section from DIA-NN logs
    if (!is.null(comp_results$diann_log_a) && !is.null(comp_results$diann_log_a$n_precursors_library)) {
      log_a <- comp_results$diann_log_a
      log_b <- comp_results$diann_log_b
      paste0("\nDIA-NN LIBRARY SIZES:\n",
        "Run A: ", format(log_a$n_precursors_library, big.mark = ","), " precursors",
        if (!is.null(log_a$params$pg_level)) paste0(" (pg-level ", log_a$params$pg_level, ")") else "",
        if (isTRUE(log_a$params$proteoforms)) ", --proteoforms" else "", "\n",
        if (!is.null(log_b) && !is.null(log_b$n_precursors_library)) {
          paste0("Run B: ", format(log_b$n_precursors_library, big.mark = ","), " precursors",
            if (!is.null(log_b$params$pg_level)) paste0(" (pg-level ", log_b$params$pg_level, ")") else "",
            if (isTRUE(log_b$params$proteoforms)) ", --proteoforms" else "", "\n")
        } else "Run B: log not provided\n"
      )
    } else "",
    "\n",
    "SETTINGS DIFFERENCES:\n", settings_summary, "\n\n",
    "TOP DISCORDANT PROTEINS (protein_id | logFC_A | logFC_B | adjP_A | adjP_B | category):\n",
    disc_text, "\n",
    mofa_section, "\n\n",
    "Structure your response EXACTLY as follows:\n\n",
    "## 1. Factual Observations\n",
    "List only what the data shows — concordance rates, fold-change distributions, correlation ",
    "ranges, number of significant proteins per tool, global offset. No interpretation yet. ",
    "Cite specific numbers from the comparison overview above.\n\n",
    "## 2. Sources of Disagreement\n",
    "For each major source of disagreement, state:\n",
    "- What the difference is (quantification method, statistical model, normalization, etc.)\n",
    "- How many proteins it likely affects (systematic vs protein-specific)\n",
    "- The expected direction of its effect (more/fewer DE calls, shifted fold-changes, etc.)\n",
    "Ground claims in proteomics literature where possible (e.g., limma empirical Bayes: Smyth 2004; ",
    "MaxLFQ: Cox et al. 2014; DIA-NN: Demichev et al. 2020; Spectronaut: Bruderer et al. 2015).\n\n",
    "## 3. Case for Run A\n",
    "Argue the strongest case for trusting Run A's results for this specific experiment. ",
    "What does Run A do well here? What in the data supports this?\n\n",
    "## 4. Case for Run B\n",
    "Argue the strongest case for trusting Run B's results for this specific experiment. ",
    "What does Run B do well here? What in the data supports this?\n\n",
    "## 5. Settings Audit\n",
    "Review the SETTINGS DIFFERENCES table above. For each difference:\n",
    "- Could this specific setting explain any of the observed discrepancies?\n",
    "- Are any settings misconfigured, unusual, or potentially erroneous for this type of experiment?\n",
    "- Flag anything that looks like a mistake (e.g., wrong enzyme, unexpected mass accuracy, ",
    "mismatched FASTA, unusual normalization choice for the sample size).\n",
    "Be specific about which discrepancies each setting could cause.\n\n",
    "## 6. Biology of Concordant Proteins\n",
    "The concordant DE proteins (significant in both pipelines with the same direction) represent ",
    "the highest-confidence biological signal. Briefly characterize:\n",
    "- Are they enriched in specific pathways, compartments, or functional categories?\n",
    "- Do the concordant proteins suggest a coherent biological narrative for this comparison?\n",
    "- Are there notable proteins in the concordant set that are well-known markers or ",
    "functionally relevant?\n",
    "If the protein IDs are UniProt accessions rather than gene symbols, note any you recognize.\n\n",
    "## 7. Synthesis\n",
    "Now weigh the cases from sections 3-4. Where do the pipelines agree? Where is the evidence ",
    "genuinely ambiguous? If one pipeline is clearly more appropriate for this experiment design ",
    "(sample size, data type), explain why with reference to the literature. If the evidence does ",
    "not clearly favor one, say so explicitly — do not force a recommendation.\n\n",
    "## 8. Recommended Follow-ups\n",
    "Two concrete, actionable follow-up analyses — one investigating a potential issue in each ",
    "pipeline (not both focused on the same tool). Be specific about which proteins to examine, ",
    "which settings to change, or which validation to run.\n\n",
    "GUIDELINES:\n",
    "- Both pipelines are established and peer-reviewed. Do not assume one is inherently superior.\n",
    "- Every claim must be supported by a specific number from the data or a literature reference.\n",
    "- If a setting (e.g., Quant3) is unusual, discuss its potential impact but do not treat it ",
    "as automatically invalidating results — the actual effect depends on the experiment.\n",
    "- Label speculation clearly: 'This suggests...' or 'One possible explanation is...' ",
    "vs stating as fact.\n",
    "- If discordant proteins show a clear pattern (e.g., all low-abundance, all membrane proteins), ",
    "note it. If no pattern is evident, say so."
  )
}

#' Build Claude export prompt
build_claude_comparator_prompt <- function(comp_results, gemini_narrative = NULL, instrument_meta = NULL) {
  stats    <- comp_results$summary_stats
  source_b <- stats$source_b

  tool_label <- switch(source_b,
    "delimp"           = "two DE-LIMP sessions (same pipeline, different settings)",
    "spectronaut"      = "DE-LIMP (DIA-NN/DPC-Quant/limma) vs Spectronaut",
    "fragpipe_analyst" = "DE-LIMP (DIA-NN/DPC-Quant/limma) vs FragPipe+FragPipe-Analyst (MSFragger/MaxLFQ/limma)",
    "fragpipe_raw"     = "DE-LIMP (DIA-NN) vs FragPipe raw output (no DE stats)",
    "unknown comparison"
  )

  # Dominant hypothesis category
  dom_hyp <- if (!is.null(comp_results$de_concordance) &&
                 nrow(comp_results$de_concordance$discordant_table) > 0) {
    cats <- comp_results$de_concordance$discordant_table$hypothesis_category
    names(sort(table(cats), decreasing = TRUE))[1]
  } else "N/A"

  gemini_section <- if (!is.null(gemini_narrative)) {
    paste0("\n\nGEMINI PRE-ANALYSIS:\n", gemini_narrative)
  } else ""

  # Instrument context (brief)
  instrument_line <- ""
  if (!is.null(instrument_meta)) {
    inst_parts <- character(0)
    if (!is.null(instrument_meta$instrument_model) && nzchar(instrument_meta$instrument_model))
      inst_parts <- c(inst_parts, instrument_meta$instrument_model)
    if (!is.null(instrument_meta$lc_system) && nzchar(instrument_meta$lc_system))
      inst_parts <- c(inst_parts, instrument_meta$lc_system)
    if (!is.null(instrument_meta$lc_method) && nzchar(instrument_meta$lc_method))
      inst_parts <- c(inst_parts, instrument_meta$lc_method)
    if (!is.null(instrument_meta$evosep_spd) && !is.na(instrument_meta$evosep_spd))
      inst_parts <- c(inst_parts, paste0(instrument_meta$evosep_spd, " SPD"))
    if (!is.null(instrument_meta$evosep_gradient_min) && !is.na(instrument_meta$evosep_gradient_min))
      inst_parts <- c(inst_parts, paste0(instrument_meta$evosep_gradient_min, " min gradient"))
    if (length(inst_parts) > 0)
      instrument_line <- paste0("INSTRUMENT: ", paste(inst_parts, collapse = ", "), "\n")
  }

  paste0(
    "I am sharing a proteomics run comparison from DE-LIMP.\n\n",
    instrument_line,
    "COMPARISON TYPE: ", tool_label, "\n",
    "EXPERIMENT: ", stats$contrast, " contrast, ", stats$n_samples, " samples.\n\n",
    "KEY FINDING: Of ", stats$n_shared, " shared proteins, ",
    stats$n_discordant, " are discordant",
    " (", stats$n_a_only_de, " DE in Run A only, ", stats$n_b_only_de, " DE in Run B only).\n",
    "Most common discordance pattern: ", dom_hyp, ".\n",
    "Global intensity offset: ", round(stats$global_offset, 3), " log2 units",
    " (", if (abs(stats$global_offset) > 0.2) "SYSTEMATIC BIAS DETECTED" else "no systematic bias", ").\n",
    gemini_section, "\n\n",
    "FILES ATTACHED:\n",
    "- discordant_proteins.csv: All disagreements with per-protein diagnostic flags",
    if (source_b == "spectronaut") " (includes n_ratios_B column)" else "", "\n",
    "- de_results_combined.csv: Full DE stats from both runs side-by-side\n",
    "- settings_diff.csv: Parameter comparison table\n",
    "- protein_universe.csv: All proteins with tier classification\n",
    "- diann_search_params.txt: DIA-NN search parameters for Run A (if available)\n",
    "- precursor_summary_discordant.csv: Precursor-level data for discordant proteins (if available)\n",
    if (source_b == "spectronaut") {
      paste0(
        "- spectronaut_run_qc.csv: Per-sample QC metrics from Spectronaut (if available)\n",
        "- spectronaut_library_info.csv: Spectronaut library composition and version\n")
    } else "",
    "\nQUESTIONS:\n",
    "1. Based on the discordant protein patterns and tool differences, what is the most likely root cause?\n",
    "2. Are there proteins where the two runs disagree biologically (not just statistically)?\n",
    if (source_b == "spectronaut") {
      paste0(
        "3. Do low n_ratios_B proteins in discordant_proteins.csv cluster biologically, or are they random?\n",
        "4. If spectronaut_run_qc.csv is present, are outlier samples driving specific discordant proteins?\n",
        "5. Among proteins significant only in DE-LIMP: do they have high or low n_ratios_B in Spectronaut?\n")
    } else {
      paste0(
        "3. If precursor_summary_discordant.csv is included, do any discordant proteins have unusual ",
        "precursor characteristics that might explain the disagreement?\n",
        "4. What additional information would resolve the ambiguity?\n")
    }
  )
}


# ==============================================================================
#  SERVER MODULE
# ==============================================================================

server_comparator <- function(input, output, session, values, add_to_log) {

  # --- Info modal ---
  # --- Spectronaut Setup Guide download ---
  output$download_spectronaut_schema <- downloadHandler(
    filename = function() "DE-LIMP_Spectronaut_Export_Guide.txt",
    content = function(file) {
      writeLines(c(
        "==========================================================",
        "  DE-LIMP Run Comparator: Spectronaut Export Guide",
        "  For use with Mode B: DE-LIMP vs Spectronaut comparison",
        "  Spectronaut 20 / 20.5 | Updated March 2026",
        "==========================================================",
        "",
        "OVERVIEW",
        "--------",
        "The DE-LIMP Run Comparator (Mode B) accepts a single ZIP file",
        "exported directly from Spectronaut's Report Exporter. This ZIP",
        "contains everything the comparator needs: protein quantities,",
        "DE results, per-run QC summaries, settings, and the sample name",
        "map that bridges Spectronaut's short labels to your original",
        "raw file names.",
        "",
        "  One ZIP, no manual column selection. The de-limp report schema",
        "  handles all required columns automatically.",
        "",
        "",
        "QUICK START (3 STEPS)",
        "---------------------",
        "",
        "  Step  Action                                Where in Spectronaut",
        "  ----  ------------------------------------  ------------------------------------",
        "  1     Install the de-limp report schema     Report Perspective > Import Schema",
        "  2     Run Report Exporter with checkboxes   Pipeline > Report Exporter",
        "  3     Upload the ZIP to DE-LIMP             Run Comparator > Mode B",
        "",
        "",
        "STEP 1 -- Install the de-limp Report Schema",
        "============================================",
        "",
        "The de-limp report schema is a pre-configured Spectronaut report",
        "layout that exports exactly the columns DE-LIMP needs. Install it",
        "once and it will be available for all future experiments.",
        "",
        "Download the schema file:",
        "  https://github.com/bsphinney/DE-LIMP/releases/download/v3.5.0/de-limp.rs",
        "  (or find it in the DE-LIMP GitHub Releases page under current release assets)",
        "",
        "Import into Spectronaut:",
        "  1. Open Spectronaut, go to Report Perspective (toolbar icon or View > Report)",
        "  2. Click the Import button (folder icon) in the Report Schema panel",
        "  3. Browse to de-limp.rs and click Open",
        "  4. The schema named 'de-limp (Pivot)' will appear in the Run Reports list",
        "",
        "  NOTE: The schema only needs to be installed once per Spectronaut installation.",
        "",
        "",
        "STEP 2 -- Export from Report Exporter",
        "=====================================",
        "",
        "Spectronaut's Report Exporter (new in Spectronaut 18+) exports",
        "multiple report types and ancillary files into a single ZIP in",
        "one operation. This is the recommended export method for DE-LIMP.",
        "",
        "Opening Report Exporter:",
        "  - From the Pipeline Perspective: click Report Exporter in the toolbar",
        "  - Or right-click experiment > Export > Report Exporter",
        "  - Or from Post Analysis Perspective: Actions > Export > Report Exporter",
        "",
        "What to check in Report Exporter:",
        "",
        "  Run Reports",
        "    [x] de-limp (Pivot)           REQUIRED   Pivot report with PG.Quantity and",
        "                                             PG.NrOfStrippedSequencesIdentified",
        "    [ ] BGS Factory Report         optional   Not needed for DE-LIMP comparator",
        "",
        "  Log and Setup Reports",
        "    [x] Spectronaut Log            REQUIRED   Version info, unique precursor/protein counts",
        "    [x] Experiment Settings Report  REQUIRED   FASTA, enzyme, modifications, FDR, normalization",
        "    [x] Run Meta Report            REQUIRED   Per-sample QC: precursors, peptides, PGs,",
        "                                             mass accuracy, instrument, gradient length",
        "    [x] Post Analysis Reports      REQUIRED   Includes Candidates.tsv (DE results: logFC,",
        "                                             p-value, q-value per protein group)",
        "    [x] Condition Setup Report     REQUIRED   Maps Spectronaut short labels (e.g. AD-1) to",
        "                                             raw file names -- critical for sample matching",
        "    [ ] Memory Consumption Report  skip       Not used by DE-LIMP",
        "",
        "  Other Reports",
        "    [ ] Spectronaut Experiment     skip       Very large file -- not needed",
        "        Save File (.sne)",
        "    [ ] Normalization Report       skip       PDF only -- not machine-readable",
        "",
        "  Leave all Calibration Reports unchecked -- they are not used by DE-LIMP.",
        "",
        "Click Export:",
        "  5. Set the Destination folder at the top of the dialog",
        "  6. Check the boxes as shown above",
        "  7. Click the Export button at the bottom",
        "  8. Spectronaut creates a ZIP named after your experiment",
        "",
        "  NOTE: The Run Meta Report generates one TSV per sample in a RunSummaries/",
        "  subfolder. This is normal -- they are all included in the ZIP automatically.",
        "",
        "",
        "STEP 3 -- Upload to DE-LIMP Run Comparator",
        "===========================================",
        "",
        "  9.  Open DE-LIMP and go to the Run Comparator tab",
        "  10. Select Mode B: DE-LIMP vs Spectronaut",
        "  11. Under 'Spectronaut Export (ZIP)', click Browse and select the ZIP",
        "  12. DE-LIMP will auto-detect and parse all included files",
        "  13. Select your DE-LIMP session or upload a session .rds file as Run A",
        "  14. Click Run Comparison",
        "",
        "  The comparator uses ConditionSetup.tsv to automatically map Spectronaut",
        "  sample labels to DE-LIMP column names -- no manual matching needed.",
        "",
        "",
        "WHAT EACH FILE ENABLES",
        "----------------------",
        "",
        "  File in ZIP                     Enables in Comparator",
        "  ----------------------------    -------------------------------------------",
        "  de-limp (Pivot).tsv             Protein Universe, Quantification, Peptides",
        "  ConditionSetup.tsv              Automatic sample name matching (CRITICAL)",
        "  Candidates.tsv                  DE Concordance: logFC/q-value scatter",
        "  ExperimentSetupOverview_*.txt   Settings Diff: FASTA, enzyme, mods, FDR",
        "  AnalysisLog.txt                 Software version, unique precursor/protein counts",
        "  RunSummaries/*.tsv              Per-sample QC: precursors, peptides, PGs",
        "",
        "",
        "FALLBACK: Manual Export (without Report Exporter)",
        "=================================================",
        "",
        "If using Spectronaut 17 or earlier, or prefer manual export, you can",
        "export files individually. Upload via 'Or upload individual files...'",
        "",
        "File 1: Protein Quantities (required)",
        "  Export from: Report Perspective > de-limp (Pivot) schema > Export",
        "  Columns: PG.ProteinGroups, PG.Genes, [N] Sample.PG.Quantity,",
        "           [N] Sample.PG.NrOfStrippedSequencesIdentified",
        "",
        "  NOTE: Spectronaut 20 changed the peptide count column name from",
        "  PG.NrOfStrippedSequences to PG.NrOfStrippedSequencesIdentified.",
        "  DE-LIMP detects both names automatically.",
        "",
        "File 2: DE Results (required for DE Concordance)",
        "  Export from: Post Analysis > Differential Abundance > Candidates > Export Table",
        "  Set the q-value filter to 1.0 before exporting to get ALL proteins.",
        "  Required columns: ProteinGroups, AVG Log2 Ratio, Pvalue, Qvalue, # of Ratios",
        "",
        "  NOTE: In Spectronaut 20 the Candidates export moved from Report Perspective",
        "  to Post Analysis Perspective > Differential Abundance section.",
        "",
        "File 3: Condition Setup (required for sample matching)",
        "  Export from: right-click experiment > Export > Condition Setup",
        "  Maps Spectronaut short labels (e.g. 'AD-1') to raw file names.",
        "  Without this file, matching falls back to fuzzy suffix matching.",
        "",
        "File 4: Experiment Settings (recommended)",
        "  Export from: right-click experiment tab > Export Experiment Settings",
        "  Extracted fields: FASTA, enzyme, missed cleavages, modifications,",
        "  FDR thresholds, normalization strategy, DE test type",
        "",
        "",
        "TROUBLESHOOTING",
        "---------------",
        "",
        "  Problem                          Likely cause / Fix",
        "  -------------------------------- -----------------------------------------",
        "  'No samples matched' error       ConditionSetup.tsv missing from ZIP.",
        "                                   Re-export with Condition Setup Report checked.",
        "",
        "  Sample names all unmatched       Labels differ from raw file stems.",
        "                                   Check the sample map table and correct manually.",
        "",
        "  Candidates.tsv has 0 rows        q-value filter set to 0.05 before export.",
        "                                   Set Confidence Candidate Filter to 1.0, re-export.",
        "",
        "  Settings diff shows 'not         Experiment Settings Report not in ZIP.",
        "  detected'                        Re-export with that checkbox checked.",
        "",
        "  de-limp schema not visible       Schema not yet imported.",
        "  in Report Exporter               Import de-limp.rs via Report Perspective.",
        "",
        "  ZIP upload fails                 Check file size (<500 MB). Spaces in",
        "                                   experiment names are handled automatically.",
        "",
        "  AD12. / AD14. trailing dot       Spectronaut appends dot when label ends in",
        "  in sample names                  digit. DE-LIMP strips trailing dots during matching.",
        "",
        "",
        "COLUMN NAME REFERENCE",
        "---------------------",
        "",
        "DE-LIMP auto-detects column names across Spectronaut versions:",
        "",
        "  Data Element       Spectronaut 20 column              Fallback / alternate",
        "  -----------------  -----------------------------------  --------------------------",
        "  Protein group ID   PG.ProteinGroups                     PG.UniProtIds",
        "  Gene symbol        PG.Genes                             --",
        "  Per-sample qty     [N] Sample.PG.Quantity               Sample.PG.Quantity",
        "  Peptide count      [N] Sample.PG.NrOfStripped-          PG.NrOfStrippedSequences",
        "                     SequencesIdentified                  (Spectronaut <=19)",
        "  DE logFC           AVG Log2 Ratio                       AVG.Log2.Ratio",
        "  DE p-value         Pvalue                               PValue",
        "  DE q-value (FDR)   Qvalue                               Q.Value",
        "  Number of ratios   # of Ratios                          NrOfRatios",
        "  Comparison label   Comparison (group1/group2)           Comparison",
        "",
        "Generated by DE-LIMP v3.5.0 | UC Davis Proteomics Core",
        "https://github.com/bsphinney/DE-LIMP"
      ), file)
    }
  )

  observeEvent(input$comparator_info_btn, {
    showModal(modalDialog(
      title = "Run Comparator",
      tags$p("Compare two analyses of the same samples to understand why DE results differ."),
      tags$h6("Modes"),
      tags$ul(
        tags$li(tags$b("DE-LIMP vs DE-LIMP"), " - Compare two sessions with different settings"),
        tags$li(tags$b("DE-LIMP vs Spectronaut"), " - Cross-platform DE comparison"),
        tags$li(tags$b("DE-LIMP vs FragPipe"), " - Cross-pipeline comparison (with or without DE stats)")
      ),
      tags$h6("Diagnostic Layers"),
      tags$ol(
        tags$li(tags$b("Settings Diff"), " \u2014 Side-by-side parameter comparison"),
        tags$li(tags$b("Protein Universe"), " \u2014 Overlap of detected proteins (shared vs run-specific)"),
        tags$li(tags$b("Quantification"), " \u2014 Intensity correlation scatter plot, systematic bias detection"),
        tags$li(tags$b("DE Concordance"), " \u2014 3\u00d73 matrix classifying every shared protein as ",
          tags$b("Up"), " (significant, logFC > 0), ",
          tags$b("Down"), " (significant, logFC < 0), or ",
          tags$b("NS"), " (Not Significant). ",
          "Concordant = same classification in both runs. Discordant proteins are analyzed by a ",
          "hypothesis engine that suggests likely causes (normalization offset, borderline significance, ",
          "missing values, direction reversal, etc.).")
      ),
      tags$h6("AI Analysis"),
      tags$p("Generate a Gemini summary of the comparison, optionally run MOFA2 factor decomposition, ",
             "or export a ZIP for Claude analysis."),
      tags$hr(),
      tags$p(class = "text-muted small",
        icon("lightbulb"), " Tip: Click the ", icon("question-circle"),
        " button next to the Hypothesis Distribution chart in the DE Concordance tab ",
        "for detailed explanations of each hypothesis category."),
      easyClose = TRUE, size = "l",
      footer = modalButton("Close")
    ))
  })

  # --- Hypothesis category info modal ---
  observeEvent(input$comparator_hypothesis_info_btn, {
    showModal(modalDialog(
      title = "Hypothesis Categories Explained",
      tags$p("When two runs disagree on a protein, the hypothesis engine examines the data to suggest ",
             tags$b("why"), " they disagree. Each discordant protein is assigned one of these categories:"),
      tags$dl(style = "margin-top: 10px;",
        tags$dt(style = "color: #c0392b;", "Direction reversal"),
        tags$dd("The protein's fold-change points in ", tags$b("opposite directions"),
                " between runs (e.g., up in Run A, down in Run B). Usually indicates a quantification ",
                "or normalization discrepancy at the precursor level. Requires careful investigation."),

        tags$dt(style = "color: #e67e22;", "Normalization offset"),
        tags$dd("A ", tags$b("global intensity shift"), " was detected between runs (e.g., all logFC values ",
                "shifted by ~0.2 log2 units). The protein has the same fold-change direction in both runs, ",
                "but one tool's normalization pushed it across the significance threshold. ",
                "Common when comparing tools with different normalization strategies (e.g., Spectronaut's ",
                "local regression vs DIA-NN's median centering)."),

        tags$dt(style = "color: #8e44ad;", "Variance estimation"),
        tags$dd("The protein has a ", tags$b("similar fold-change in both runs"), " but very different ",
                "p-values. This happens when the two tools use different statistical models ",
                "(e.g., limma's empirical Bayes variance shrinkage vs Spectronaut's unpaired t-test). ",
                "Proteins near the significance threshold are most affected. ",
                "If Spectronaut used 'Use All MS-Level Quantities' (Quant3), this effectively doubles ",
                "the observation count, making its p-values artificially small."),

        tags$dt(style = "color: #2980b9;", "Missing values"),
        tags$dd("One tool detected the protein in significantly ", tags$b("more samples"), " than the other. ",
                "Different match-between-runs (MBR) implementations and imputation strategies ",
                "produce different missingness patterns, which directly affects statistical power."),

        tags$dt(style = "color: #16a085;", "Low ratio count"),
        tags$dd(tags$b("Spectronaut-specific."), " Spectronaut reports how many valid pairwise ratios ",
                "it could form for each protein. A low ratio count (e.g., 5 out of 420 possible) means ",
                "the protein had high missingness in Spectronaut's quantification. Its significance call ",
                "may be underpowered relative to limma, which borrows strength across proteins."),

        tags$dt(style = "color: #27ae60;", "Peptide count"),
        tags$dd("The two tools quantified the protein using a ", tags$b("different number of peptides"), ". ",
                "This matters because protein-level intensity is aggregated from peptide/precursor intensities. ",
                "More peptides generally means more stable quantification. Spectronaut's TopN (default 3) ",
                "caps the number of peptides used, while DIA-NN's DPC-Quant uses all detected precursors."),

        tags$dt(style = "color: #d35400;", "FC magnitude"),
        tags$dd("Same direction, but the ", tags$b("magnitude differs substantially"), " (>0.5 log2). ",
                "Usually a normalization scale difference between tools."),

        tags$dt(style = "color: #7f8c8d;", "Borderline"),
        tags$dd("The protein sits ", tags$b("near the significance threshold in both runs"), ". ",
                "Small stochastic differences (noise, imputation, FDR correction) push it above the ",
                "cutoff in one run and below in the other. This is expected and not actionable ",
                "\u2014 neither result is wrong.")
      ),
      easyClose = TRUE, size = "l",
      footer = modalButton("Close")
    ))
  })

  # --- Sub-tab info modals ---

  observeEvent(input$comparator_settings_info_btn, {
    showModal(modalDialog(
      title = "Settings Diff",
      tags$p("Side-by-side comparison of analysis parameters between the two runs."),
      tags$h6("Color Coding"),
      tags$dl(
        tags$dt(style = "color: #856404;", "Amber (differs)"),
        tags$dd("The parameter has a different value in each run. Review whether this could ",
                "contribute to discordant DE results."),
        tags$dt(style = "color: #004085;", "Blue (structural difference)"),
        tags$dd("An inherent difference between the two tools (e.g., TopN vs DPC-Quant rollup). ",
                "Expected and not necessarily problematic."),
        tags$dt(style = "color: #721c24;", "Red (severe)"),
        tags$dd("A setting that is known to cause major comparability problems. For example, ",
                "Spectronaut's 'Use All MS-Level Quantities' (Quant3) doubles effective observations ",
                "in its built-in t-test, inflating significance counts."),
        tags$dt(style = "color: #6c757d;", "Grey (unknown)"),
        tags$dd("The parameter could not be extracted from one or both runs.")
      ),
      tags$h6("Per-Sample QC (Spectronaut Mode)"),
      tags$p("When Spectronaut RunSummaries are available, a per-sample QC bar chart and table ",
             "appear below the settings diff. Samples with fewer than 60% of the median protein count ",
             "are flagged as potential outliers."),
      easyClose = TRUE, size = "l",
      footer = modalButton("Close")
    ))
  })

  observeEvent(input$comparator_universe_info_btn, {
    showModal(modalDialog(
      title = "Protein Universe",
      tags$p("Shows the overlap of quantified proteins between the two runs."),
      tags$h6("Venn Diagram"),
      tags$p("Circle areas are proportional to the number of proteins in each run. ",
             "The overlap region shows proteins quantified in both runs (shared). ",
             "Non-overlapping regions show proteins unique to one run."),
      tags$h6("Why Protein Universes Differ"),
      tags$ul(
        tags$li(tags$b("Search space:"), " Different FASTA databases, missed cleavages, or variable modifications ",
                "can cause one tool to identify more proteins."),
        tags$li(tags$b("Library strategy:"), " Spectronaut directDIA+ and DIA-NN empirical libraries use different ",
                "algorithms to build spectral libraries, leading to different peptide/protein coverage."),
        tags$li(tags$b("FDR control:"), " Different FDR methods may let slightly different protein sets pass the cutoff."),
        tags$li(tags$b("Match-between-runs:"), " Different MBR implementations transfer IDs across runs differently.")
      ),
      tags$h6("Protein Details Table"),
      tags$p("Use the filter buttons (All / Shared / Run A only / Run B only) to explore which proteins ",
             "are detected by each tool. The Gene column shows gene symbols where available."),
      easyClose = TRUE, size = "l",
      footer = modalButton("Close")
    ))
  })

  observeEvent(input$comparator_quant_info_btn, {
    showModal(modalDialog(
      title = "Quantification",
      tags$p("Compares protein-level intensities between the two runs for all shared proteins."),
      tags$h6("Scatter Plot"),
      tags$p("Each dot is a protein, plotted by its mean log2 intensity in Run A (x-axis) vs Run B (y-axis). ",
             "Points near the diagonal indicate agreement. A systematic shift off the diagonal suggests ",
             "a normalization difference."),
      tags$h6("Per-Sample Correlation"),
      tags$p("Bar chart showing Pearson correlation for each matched sample pair between the two runs. ",
             "High correlation (r > 0.95) indicates the tools largely agree on relative protein abundances."),
      tags$h6("Bias Density"),
      tags$p("Distribution of per-protein intensity offsets (Run A \u2212 Run B). ",
             "A distribution centered near zero means no systematic bias. ",
             "A shifted distribution indicates one tool consistently reports higher intensities."),
      tags$h6("TopN Effect (Spectronaut Mode)"),
      tags$p("When Spectronaut uses TopN quantification (e.g., Top 3 peptides), proteins with more ",
             "identified peptides than the TopN cap have extra signal discarded. This scatter plot shows ",
             "whether quantification divergence increases for proteins above the TopN threshold."),
      tags$ul(
        tags$li(tags$b("X-axis:"), " Number of peptides identified by Spectronaut (log scale)"),
        tags$li(tags$b("Y-axis:"), " Absolute log2 quantification difference between the two tools"),
        tags$li(tags$b("Vertical dashed line:"), " The TopN cap (e.g., 3). Proteins to the right ",
                "have more peptides than Spectronaut uses for quantification."),
        tags$li(tags$b("LOESS curve:"), " If the curve rises to the right of the dashed line, ",
                "TopN capping is a driver of quantification divergence.")
      ),
      easyClose = TRUE, size = "l",
      footer = modalButton("Close")
    ))
  })

  observeEvent(input$comparator_concordance_info_btn, {
    showModal(modalDialog(
      title = "DE Concordance",
      tags$p("Classifies every shared protein by its differential expression status in each run, ",
             "then identifies why discordant proteins disagree."),
      tags$h6("3\u00d73 Concordance Matrix"),
      tags$p("Each protein is classified as:"),
      tags$ul(
        tags$li(tags$b("Up"), " \u2014 significant (adj. P < 0.05) with positive logFC"),
        tags$li(tags$b("Down"), " \u2014 significant with negative logFC"),
        tags$li(tags$b("NS"), " \u2014 not significant (adj. P \u2265 0.05)")
      ),
      tags$p("Diagonal cells (Up/Up, Down/Down, NS/NS) are ", tags$b("concordant"),
             " \u2014 both runs agree. Off-diagonal cells are ", tags$b("discordant"), "."),
      tags$h6("Hypothesis Engine"),
      tags$p("Each discordant protein is analyzed by an 8-rule hypothesis engine that examines ",
             "logFC values, p-values, peptide counts, and other evidence to suggest the most likely ",
             "cause of disagreement. Click the ", icon("question-circle"),
             " next to 'Hypothesis Distribution' for detailed descriptions of each category."),
      tags$h6("Volcano Overlay"),
      tags$p("Superimposes both runs' volcano plots. Blue dots = Run A, orange = Run B. ",
             "This reveals systematic differences in fold-change distributions and p-value scales."),
      tags$h6("Discordant Proteins Table"),
      tags$p("Lists each discordant protein with its logFC and adj. P-value in both runs, ",
             "the assigned hypothesis, and confidence level. Sort or filter to find proteins of interest."),
      easyClose = TRUE, size = "l",
      footer = modalButton("Close")
    ))
  })

  observeEvent(input$comparator_ai_info_btn, {
    showModal(modalDialog(
      title = "AI Analysis",
      tags$h6("Gemini Summary"),
      tags$p("Generates a narrative interpretation of the comparison using Google Gemini. ",
             "The prompt includes the settings diff, concordance statistics, systematic bias metrics, ",
             "and top discordant proteins with hypotheses. Requires a Gemini API key in AI Chat settings."),
      tags$h6("Export ZIP for Claude"),
      tags$p("Downloads a ZIP file containing structured CSVs and a pre-written prompt ",
             "for analysis in Claude or another LLM. Includes settings diff, protein universe, ",
             "DE results, discordant proteins with hypotheses, and comparison context."),
      tags$h6("MOFA2 Decomposition"),
      tags$p("Treats the two runs as two 'views' of the same samples and uses MOFA2 to decompose ",
             "joint variance into shared and run-specific factors. Helps identify whether differences ",
             "are driven by a global shift (one factor) or multiple independent sources."),
      tags$p(class = "text-muted small",
        "Requires \u2265 4 matched sample pairs. Runtime: ~1\u20132 minutes."),
      easyClose = TRUE, size = "l",
      footer = modalButton("Close")
    ))
  })

  # --- Internal reactive values ---
  comp_run_a <- reactiveVal(NULL)
  comp_run_b <- reactiveVal(NULL)
  comp_results <- reactiveVal(NULL)
  comp_diann_log_a <- reactiveVal(NULL)
  comp_diann_log_b <- reactiveVal(NULL)

  # --- DIA-NN log upload observers ---
  observeEvent(input$comparator_diann_log_a, {
    req(input$comparator_diann_log_a)
    parsed <- tryCatch(parse_diann_log(input$comparator_diann_log_a$datapath),
                       error = function(e) list(success = FALSE))
    comp_diann_log_a(if (isTRUE(parsed$success)) parsed else NULL)
    update_diann_log_status()
  })

  observeEvent(input$comparator_diann_log_b, {
    req(input$comparator_diann_log_b)
    parsed <- tryCatch(parse_diann_log(input$comparator_diann_log_b$datapath),
                       error = function(e) list(success = FALSE))
    comp_diann_log_b(if (isTRUE(parsed$success)) parsed else NULL)
    update_diann_log_status()
  })

  # Log status display
  update_diann_log_status <- function() {
    log_a <- comp_diann_log_a()
    log_b <- comp_diann_log_b()
    items <- character()
    for (run in c("A", "B")) {
      p <- if (run == "A") log_a else log_b
      if (!is.null(p)) {
        ver <- p$version %||% "unknown"
        nprec <- p$n_precursors_library
        nprec_str <- if (!is.null(nprec)) format(nprec, big.mark = ",") else "?"
        pg <- p$params$pg_level
        pg_str <- if (!is.null(pg)) paste0(" | pg-level ", pg) else ""
        proto <- if (isTRUE(p$params$proteoforms)) " | proteoforms" else ""
        rean <- if (isTRUE(p$params$mbr)) " | reanalyse" else ""
        items <- c(items, sprintf(
          "<b>Run %s:</b> DIA-NN %s | %s precursors%s%s%s",
          run, ver, nprec_str, pg_str, proto, rean
        ))
      }
    }
    html <- if (length(items) > 0) {
      paste0('<div class="small text-muted">',
             paste(items, collapse = "<br>"), '</div>')
    } else ""
    shinyjs::html("comparator_diann_log_status", html)
  }

  # --- Parse Run A on source change or file upload ---
  observe({
    if (isTRUE(input$comparator_run_a_source == "current")) {
      req(values$fit, values$y_protein)
      tryCatch({
        parsed <- parse_delimp_session(values = values)
        comp_run_a(parsed)
      }, error = function(e) {
        showNotification(paste("Error parsing current session:", e$message), type = "error")
        comp_run_a(NULL)
      })
    }
  })

  observeEvent(input$comparator_run_a_file, {
    req(input$comparator_run_a_source == "file")
    tryCatch({
      parsed <- parse_delimp_session(rds_path = input$comparator_run_a_file$datapath)
      comp_run_a(parsed)
      showNotification("Run A loaded successfully", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading Run A:", e$message), type = "error")
      comp_run_a(NULL)
    })
  })

  # --- Parse Run B based on mode ---
  observeEvent(input$comparator_run_b_rds, {
    req(input$comparator_mode == "delimp_delimp")
    tryCatch({
      parsed <- parse_delimp_session(rds_path = input$comparator_run_b_rds$datapath)
      comp_run_b(parsed)
      showNotification("Run B loaded successfully", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading Run B:", e$message), type = "error")
      comp_run_b(NULL)
    })
  })

  # --- Load from History (compare 2 analyses) ---
  observeEvent(values$comparator_compare_from_history, {
    hist <- values$comparator_compare_from_history
    req(hist)
    tryCatch({
      if (!is.null(values$comparator_run_a)) {
        comp_run_a(values$comparator_run_a)
        values$comparator_run_a <- NULL
      }
      if (!is.null(values$comparator_run_b)) {
        comp_run_b(values$comparator_run_b)
        values$comparator_run_b <- NULL
      }
      if (!is.null(values$comparator_diann_log_a)) {
        comp_diann_log_a(values$comparator_diann_log_a)
        values$comparator_diann_log_a <- NULL
      }
      if (!is.null(values$comparator_diann_log_b)) {
        comp_diann_log_b(values$comparator_diann_log_b)
        values$comparator_diann_log_b <- NULL
      }
      update_diann_log_status()
      values$comparator_compare_from_history <- NULL
    }, error = function(e) {
      showNotification(paste("Error loading comparison:", e$message), type = "error")
    })
  }, ignoreInit = TRUE)

  # --- Spectronaut ZIP upload (primary Mode B input) ---
  spec_manifest <- reactiveVal(NULL)

  observeEvent(input$comparator_run_b_spec_zip, {
    req(input$comparator_mode == "delimp_spectronaut")
    tryCatch({
      parsed <- parse_spectronaut_zip(input$comparator_run_b_spec_zip$datapath)
      comp_run_b(parsed)
      spec_manifest(parsed$manifest)
      n_files <- sum(c(parsed$manifest$quantities, parsed$manifest$candidates,
                       parsed$manifest$run_summaries, parsed$manifest$settings,
                       parsed$manifest$condition_setup, isTRUE(parsed$manifest$analysis_log)))
      showNotification(
        sprintf("Spectronaut ZIP loaded: %d components detected", n_files),
        type = "message"
      )
    }, error = function(e) {
      showNotification(paste("Error parsing Spectronaut ZIP:", e$message), type = "error")
      comp_run_b(NULL)
      spec_manifest(NULL)
    })
  })

  # Manifest display (green checkmarks for detected files)
  output$spectronaut_zip_manifest <- renderUI({
    m <- spec_manifest()
    if (is.null(m)) return(NULL)
    check <- function(ok, label) {
      icon_tag <- if (ok) tags$span("\u2705", style = "margin-right: 4px;")
                  else     tags$span("\u274C", style = "margin-right: 4px;")
      tags$div(icon_tag, tags$small(label))
    }
    div(style = "background: #f0fdf4; border: 1px solid #bbf7d0; border-radius: 6px; padding: 8px 10px; margin-bottom: 8px;",
      tags$small(tags$b("Found:")),
      check(m$quantities, "Quantities"),
      check(m$candidates, "DE Results"),
      check(m$run_summaries, paste0("Run Summaries (", m$n_run_summaries, ")")),
      check(m$settings, "Settings"),
      check(m$condition_setup, "Sample Map"),
      check(isTRUE(m$analysis_log), "Analysis Log")
    )
  })

  # --- Spectronaut individual TSV uploads (fallback) ---
  observeEvent(input$comparator_run_b_spectronaut, {
    req(input$comparator_mode == "delimp_spectronaut")
    tryCatch({
      de_path <- input$comparator_run_b_spectronaut_de$datapath
      parsed <- parse_spectronaut(input$comparator_run_b_spectronaut$datapath,
                                  de_file_path = de_path)
      comp_run_b(parsed)
      spec_manifest(NULL)  # Clear ZIP manifest when using individual files
      msg <- if (!is.null(de_path)) "Spectronaut quantities + DE loaded" else "Spectronaut quantities loaded (no DE file)"
      showNotification(msg, type = "message")
    }, error = function(e) {
      showNotification(paste("Error parsing Spectronaut:", e$message), type = "error")
      comp_run_b(NULL)
    })
  })

  # If DE file uploaded after quantities, re-parse with both
  observeEvent(input$comparator_run_b_spectronaut_de, {
    req(input$comparator_mode == "delimp_spectronaut",
        input$comparator_run_b_spectronaut)
    tryCatch({
      parsed <- parse_spectronaut(input$comparator_run_b_spectronaut$datapath,
                                  de_file_path = input$comparator_run_b_spectronaut_de$datapath)
      comp_run_b(parsed)
      showNotification("Spectronaut DE stats added", type = "message")
    }, error = function(e) {
      showNotification(paste("Error parsing Spectronaut DE:", e$message), type = "error")
    })
  })

  observeEvent(input$comparator_fp_analyst_file, {
    req(input$comparator_mode == "delimp_fragpipe",
        input$comparator_fragpipe_type == "fp_analyst")
    tryCatch({
      parsed <- parse_fragpipe_analyst(input$comparator_fp_analyst_file$datapath)
      comp_run_b(parsed)
      showNotification("FragPipe-Analyst data loaded successfully", type = "message")
    }, error = function(e) {
      showNotification(paste("Error parsing FragPipe-Analyst:", e$message), type = "error")
      comp_run_b(NULL)
    })
  })

  observeEvent(input$comparator_fp_combined_protein, {
    req(input$comparator_mode == "delimp_fragpipe",
        input$comparator_fragpipe_type == "fp_raw")
    tryCatch({
      parsed <- parse_fragpipe_combined_protein(input$comparator_fp_combined_protein$datapath)
      comp_run_b(parsed)
      showNotification("FragPipe data loaded successfully", type = "message")
    }, error = function(e) {
      showNotification(paste("Error parsing combined_protein.tsv:", e$message), type = "error")
      comp_run_b(NULL)
    })
  })

  # --- Contrast selector UI (dynamic based on loaded runs) ---
  output$comparator_contrast_selectors <- renderUI({
    run_a <- comp_run_a()
    run_b <- comp_run_b()
    if (is.null(run_a) || is.null(run_b)) return(NULL)

    contrasts_a <- run_a$contrasts
    contrasts_b <- run_b$contrasts

    # For external tools with a single comparison, auto-fill
    if (is.null(contrasts_b) || length(contrasts_b) == 1) {
      label_b <- if (is.null(contrasts_b)) "N/A (no DE stats)" else contrasts_b
      tagList(
        div(style = "display: flex; gap: 15px; margin-top: 10px;",
          div(style = "flex: 1;",
            selectInput("comparator_contrast_a", "Run A Contrast",
                        choices = contrasts_a, selected = contrasts_a[1])
          ),
          div(style = "flex: 1;",
            tags$label("Run B Contrast", class = "control-label"),
            tags$p(label_b, class = "form-control-static text-muted",
                   style = "background: #f5f5f5; padding: 6px 12px; border-radius: 4px;")
          )
        )
      )
    } else {
      # Both runs have multiple contrasts — try matching by name
      # Normalize for fuzzy matching: lowercase, strip spaces/punctuation
      norm <- function(x) tolower(gsub("[^a-z0-9]", "", x))
      exact_match <- contrasts_a[1] %in% contrasts_b
      fuzzy_match <- if (!exact_match) {
        m <- match(norm(contrasts_a[1]), norm(contrasts_b))
        if (!is.na(m)) contrasts_b[m] else NULL
      } else NULL
      auto_match <- if (exact_match) contrasts_a[1]
                    else if (!is.null(fuzzy_match)) fuzzy_match
                    else contrasts_b[1]
      # Check if ANY contrast overlaps (exact or fuzzy)
      any_overlap <- any(norm(contrasts_a) %in% norm(contrasts_b))
      mismatch_warn <- if (!any_overlap) {
        div(class = "alert alert-warning py-2 px-3 mt-2", style = "font-size: 0.88em;",
          icon("triangle-exclamation"), " ",
          tags$b("Contrast names don't match between runs. "),
          "The ", tags$b(run_b$source %||% "external tool"),
          " export used different condition groupings than DE-LIMP. ",
          "Select the equivalent contrasts manually to compare the same biological question, ",
          "or re-export from ", tags$b(run_b$source %||% "the external tool"),
          " with matching conditions."
        )
      }
      tagList(
        div(style = "display: flex; gap: 15px; margin-top: 10px;",
          div(style = "flex: 1;",
            selectInput("comparator_contrast_a", "Run A Contrast",
                        choices = contrasts_a, selected = contrasts_a[1])
          ),
          div(style = "flex: 1;",
            selectInput("comparator_contrast_b", "Run B Contrast",
                        choices = contrasts_b, selected = auto_match)
          )
        ),
        mismatch_warn
      )
    }
  })

  # --- Sample matching status preview ---
  output$comparator_sample_status <- renderUI({
    run_a <- comp_run_a()
    run_b <- comp_run_b()
    if (is.null(run_a) || is.null(run_b)) return(NULL)
    if (is.null(run_a$intensities) || is.null(run_b$intensities)) return(NULL)

    sample_map <- match_samples(colnames(run_a$intensities),
                                colnames(run_b$intensities),
                                run_b$source,
                                meta_a = run_a$metadata,
                                meta_b = run_b$metadata)
    n_matched    <- sum(sample_map$status == "matched")
    n_unresolved <- sum(sample_map$status == "unresolved")
    n_group_matched <- sum(sample_map$match_method == "group", na.rm = TRUE)

    if (n_unresolved > 0) {
      unmatched <- sample_map$run_a[sample_map$status == "unresolved"]
      div(class = "alert alert-warning mt-2",
        icon("triangle-exclamation"), " ",
        tags$b(paste0(n_unresolved, " sample(s) could not be matched: ")),
        paste(unmatched, collapse = ", "),
        tags$br(),
        "Check that both runs analyze the same samples."
      )
    } else if (n_group_matched > 0) {
      div(class = "alert alert-info mt-2 py-2",
        icon("link"), " ",
        paste0("All ", n_matched, " samples matched by group assignment ",
               "(filenames differ — matched ", n_group_matched,
               " sample(s) by group + replicate order).")
      )
    } else {
      div(class = "alert alert-success mt-2 py-2",
        icon("check-circle"), " ",
        paste0("All ", n_matched, " samples matched between runs.")
      )
    }
  })

  # --- Main comparison pipeline ---
  observeEvent(input$run_comparison, {
    run_a <- comp_run_a()
    run_b <- comp_run_b()
    req(run_a, run_b)

    # Validate sample matching
    sample_map <- match_samples(colnames(run_a$intensities),
                                colnames(run_b$intensities),
                                run_b$source,
                                meta_a = run_a$metadata,
                                meta_b = run_b$metadata)
    if (any(sample_map$status == "unresolved")) {
      showNotification("Cannot run comparison: unresolved sample matches. Check sample names.",
                       type = "error")
      return()
    }

    # Get selected contrasts
    contrast_a <- input$comparator_contrast_a %||% run_a$contrasts[1]
    contrast_b <- input$comparator_contrast_b %||%
                  (if (!is.null(run_b$contrasts)) run_b$contrasts[1] else NULL)

    withProgress(message = "Running comparison...", value = 0, {

      # Layer 1: Settings Diff — merge log-derived fields if available
      setProgress(0.1, detail = "Comparing settings...")
      log_a <- comp_diann_log_a()
      log_b <- comp_diann_log_b()
      if (!is.null(log_a)) run_a$settings <- merge_log_into_settings(run_a$settings, log_a)
      if (!is.null(log_b)) run_b$settings <- merge_log_into_settings(run_b$settings, log_b)
      settings_diff <- build_settings_diff(run_a, run_b)

      # Layer 2: Protein Universe
      setProgress(0.25, detail = "Classifying protein universe...")
      universe <- classify_protein_universe(run_a, run_b)

      # Layer 3: Quantification Comparison
      setProgress(0.4, detail = "Comparing quantification...")
      quant <- tryCatch(
        compute_quant_comparison(run_a, run_b, universe, sample_map),
        error = function(e) {
          showNotification(paste("Quant comparison error:", e$message), type = "warning")
          NULL
        }
      )

      global_offset <- if (!is.null(quant)) quant$median_offset else 0

      # Layer 4: DE Concordance (skip if Run B has no DE stats)
      setProgress(0.6, detail = "Analyzing DE concordance...")
      de_conc <- NULL
      # Extract Spectronaut imputation setting for rescue stats
      imputation_b <- coalesce_setting(
        run_b$search_settings$imputation,
        coalesce_setting(run_b$settings$imputation, NULL)
      )

      if (!is.null(run_b$de_stats)) {
        de_conc <- tryCatch(
          compute_de_concordance(run_a, run_b, universe,
                                contrast_a, contrast_b %||% "comparison",
                                run_b$source, global_offset, imputation_b),
          error = function(e) {
            showNotification(paste("DE concordance error:", e$message), type = "warning")
            NULL
          }
        )
      }

      # Summary stats
      setProgress(0.85, detail = "Building summary...")
      n_shared <- sum(universe$tier == "shared")
      summary_stats <- list(
        source_b      = run_b$source,
        run_a_label   = paste0(run_a$settings$software, " (", run_a$settings$delimp_version %||% "", ")"),
        run_b_label   = run_b$settings$software %||% "Unknown",
        contrast      = contrast_a,
        n_samples     = ncol(run_a$intensities),
        n_proteins_a  = length(run_a$protein_ids),
        n_proteins_b  = length(run_b$protein_ids),
        n_shared      = n_shared,
        n_a_only      = sum(universe$tier == "a_only"),
        n_b_only      = sum(universe$tier == "b_only"),
        n_concordant  = if (!is.null(de_conc)) de_conc$n_concordant else NA,
        n_discordant  = if (!is.null(de_conc)) de_conc$n_discordant else NA,
        n_sig_a       = if (!is.null(de_conc)) {
          sum(de_conc$merged$status_a != "NS")
        } else NA,
        n_a_only_de   = if (!is.null(de_conc)) {
          sum(de_conc$discordant_table$status_a != "NS" & de_conc$discordant_table$status_b == "NS")
        } else NA,
        n_b_only_de   = if (!is.null(de_conc)) {
          sum(de_conc$discordant_table$status_a == "NS" & de_conc$discordant_table$status_b != "NS")
        } else NA,
        global_offset  = global_offset,
        min_cor        = if (!is.null(quant)) quant$min_cor else NA,
        max_cor        = if (!is.null(quant)) quant$max_cor else NA,
        rescue_stats   = if (!is.null(de_conc)) de_conc$rescue_stats else NULL
      )

      results <- list(
        settings_diff    = settings_diff,
        protein_universe = universe,
        quant_comparison = quant,
        de_concordance   = de_conc,
        summary_stats    = summary_stats,
        sample_map       = sample_map,
        diann_log_a      = log_a,
        diann_log_b      = log_b,
        run_b            = run_b
      )

      setProgress(1.0, detail = "Done")
      comp_results(results)
      values$comparator_results <- results
      values$comparator_run_a <- run_a
      values$comparator_run_b <- run_b
      values$comparator_mode <- input$comparator_mode
    })

    add_to_log("Run Comparator", c(
      sprintf("# Mode: %s", input$comparator_mode),
      sprintf("# Run A: %s | Run B: %s", run_a$settings$software, run_b$settings$software),
      sprintf("# Shared proteins: %d", sum(universe$tier == "shared"))
    ))
  })

  # --- Restore from session load ---
  observe({
    if (!is.null(values$comparator_results) && is.null(comp_results())) {
      comp_results(values$comparator_results)
      comp_run_a(values$comparator_run_a)
      comp_run_b(values$comparator_run_b)
    }
  })

  # ==========================================================================
  #  SUMMARY BANNER (top of results)
  # ==========================================================================

  output$comparator_summary_banner <- renderUI({
    res <- comp_results()
    req(res)
    stats <- res$summary_stats

    # Concordance rate
    conc_rate <- if (!is.na(stats$n_concordant) && !is.na(stats$n_discordant) &&
                     (stats$n_concordant + stats$n_discordant) > 0) {
      round(stats$n_concordant / (stats$n_concordant + stats$n_discordant) * 100, 1)
    } else NA

    # Dominant hypothesis
    dom_hyp <- if (!is.null(res$de_concordance) &&
                   nrow(res$de_concordance$discordant_table) > 0) {
      cats <- res$de_concordance$discordant_table$hypothesis_category
      names(sort(table(cats), decreasing = TRUE))[1]
    } else NULL

    bias_tag <- if (abs(stats$global_offset) > 0.2) {
      tags$span(class = "badge bg-warning text-dark", "Systematic bias detected")
    } else {
      tags$span(class = "badge bg-success", "No systematic bias")
    }

    # Rescue stats line (Spectronaut Mode B)
    rescue_line <- NULL
    rs <- stats$rescue_stats
    if (!is.null(rs) && rs$n_zero_ratio > 0) {
      imp_note <- if (rs$imputation_off) " (imputation disabled)" else ""
      rescue_line <- if (rs$n_rescued_sig > 0) {
        div(style = "width: 100%; font-size: 0.9em; margin-top: 4px; padding-top: 4px; border-top: 1px solid rgba(0,0,0,0.1);",
          icon("flask"), " ",
          tags$b(rs$n_rescued_sig), " proteins rescued by DE-LIMP — ",
          "untestable in Spectronaut (0 of ", rs$max_ratios, " ratios)", imp_note,
          " but significant in DE-LIMP via limpa modelling. ",
          tags$span(style = "opacity: 0.75;",
            sprintf("(%d total with 0 ratios, %d with low ratios)", rs$n_zero_ratio, rs$n_low_ratio))
        )
      } else {
        div(style = "width: 100%; font-size: 0.9em; margin-top: 4px; padding-top: 4px; border-top: 1px solid rgba(0,0,0,0.1);",
          sprintf("%d proteins untestable in Spectronaut (0 ratios)%s", rs$n_zero_ratio, imp_note)
        )
      }
    }

    div(class = "alert alert-info py-2 mb-3",
      style = "display: flex; flex-wrap: wrap; gap: 15px; align-items: center;",
      div(
        tags$b(paste0(stats$n_shared, " shared proteins")),
        " | ",
        if (!is.na(conc_rate)) paste0(conc_rate, "% concordance") else "No DE comparison",
        if (!is.na(stats$n_discordant)) paste0(" | ", stats$n_discordant, " discordant") else ""
      ),
      bias_tag,
      if (!is.null(dom_hyp)) {
        tags$span(class = "badge bg-secondary",
                  paste("Dominant cause:", dom_hyp))
      },
      rescue_line
    )
  })

  # ==========================================================================
  #  LAYER 1: Settings Diff
  # ==========================================================================

  output$comparator_settings_diff <- DT::renderDT({
    res <- comp_results()
    req(res)
    sd <- res$settings_diff
    # Add note column if missing (backwards compat)
    if (!"note" %in% names(sd)) sd$note <- ""
    DT::datatable(sd,
      rownames = FALSE,
      colnames = c("Parameter", "Run A", "Run B", "Status", "Note"),
      options  = list(
        pageLength = 50, dom = 't', ordering = FALSE,
        columnDefs = list(
          list(visible = FALSE, targets = which(names(sd) == "match") - 1),
          list(width = '300px', targets = which(names(sd) == "note") - 1)
        )
      ),
      class = "compact stripe"
    ) |>
      DT::formatStyle("match",
        target          = "row",
        backgroundColor = DT::styleEqual(
          c("differs", "unknown", "match", "structural_difference", "severe"),
          c("#fff3cd", "#f8f9fa", "transparent", "#cce5ff", "#f8d7da")
        )
      )
  })

  # Pipeline step / library size warnings above settings diff
  observe({
    res <- comp_results()
    req(res)
    log_a <- res$diann_log_a
    log_b <- res$diann_log_b
    warnings <- character()

    # Pipeline step mismatch warning
    if (!is.null(log_a) && !is.null(log_b)) {
      is_libpred <- function(s) {
        !is.null(s) && grepl("libpred|lib_pred|step.*1.*lib|Library Prediction",
                             s, ignore.case = TRUE)
      }
      if (is_libpred(log_a$pipeline_step) || is_libpred(log_b$pipeline_step)) {
        warnings <- c(warnings, paste0(
          '<div class="alert alert-warning py-2 mb-2">',
          '<i class="fas fa-triangle-exclamation"></i> <b>Log file mismatch:</b> ',
          'One of the uploaded logs appears to be a library prediction step, not a full search. ',
          'The settings comparison may not reflect the actual search parameters. ',
          'Upload the final search step log (e.g., <code>report_log.txt</code>) for an accurate comparison.',
          '</div>'
        ))
      }
    }

    # Library size mismatch note
    if (!is.null(log_a) && !is.null(log_b) &&
        !is.null(log_a$n_precursors_library) && !is.null(log_b$n_precursors_library)) {
      ratio <- max(log_a$n_precursors_library, log_b$n_precursors_library) /
               max(min(log_a$n_precursors_library, log_b$n_precursors_library), 1)
      if (ratio > 1.2) {
        warnings <- c(warnings, paste0(
          '<div class="alert alert-info py-2 mb-2">',
          '<i class="fas fa-info-circle"></i> <b>Library size difference:</b> ',
          format(log_a$n_precursors_library, big.mark = ","), ' vs ',
          format(log_b$n_precursors_library, big.mark = ","),
          ' precursors (', round(ratio, 1), 'x). ',
          'Different search spaces affect peptide-level evidence and protein rollup.',
          '</div>'
        ))
      }
    }

    # Spectronaut-specific banners (Mode B)
    run_b <- comp_run_b()
    if (!is.null(run_b) && identical(run_b$source, "spectronaut")) {
      # Quant3 danger banner
      use_all_ms <- if (!is.null(run_b$search_settings)) {
        trimws(coalesce_setting(run_b$search_settings$use_all_ms_levels, ""))
      } else if (!is.null(run_b$settings$use_all_ms_levels)) {
        trimws(run_b$settings$use_all_ms_levels)
      } else ""
      if (identical(use_all_ms, "True")) {
        warnings <- c(warnings, paste0(
          '<div class="alert alert-danger py-2 mb-2">',
          '<i class="fas fa-circle-exclamation"></i> ',
          '<b>Quant3 (Use All MS-Level Quantities) was ENABLED in Spectronaut.</b> ',
          'Spectronaut combined MS1 and MS2 measurements to double the effective ',
          'sample size in its t-test. This inflates statistical power significantly ',
          'beyond the actual number of biological replicates. p-values and q-values ',
          'from this Spectronaut analysis are NOT directly comparable to limma ',
          'moderated t-test results.',
          '</div>'
        ))
      }

      # TopN banner
      if (!is.null(run_b$n_peptides) && nrow(run_b$n_peptides) > 0) {
        topn_val <- suppressWarnings(as.integer(
          gsub("[^0-9]", "", coalesce_setting(run_b$settings$topn_protein, "3"))))
        if (is.na(topn_val)) topn_val <- 3L
        pep_counts <- run_b$n_peptides$n_peptides
        pct_limited <- round(100 * mean(pep_counts > topn_val, na.rm = TRUE))
        mean_pep <- round(mean(pep_counts, na.rm = TRUE), 1)
        if (pct_limited > 30) {
          warnings <- c(warnings, paste0(
            '<div class="alert alert-warning py-2 mb-2">',
            '<i class="fas fa-triangle-exclamation"></i> ',
            '<b>TopN=', topn_val, ' quantification is limiting for ', pct_limited,
            '% of proteins.</b> ',
            'Spectronaut used at most ', topn_val, ' peptides per protein for quantification, ',
            'but the average protein had ', mean_pep, ' peptides detected. ',
            "DE-LIMP's DPC-Quant used all detected precursors.",
            '</div>'
          ))
        }
      }
    }

    shinyjs::html("comparator_pipeline_warning", paste(warnings, collapse = ""))
  })

  # ==========================================================================
  #  LAYER 2: Protein Universe
  # ==========================================================================

  output$comparator_universe_plot <- plotly::renderPlotly({
    res <- comp_results()
    req(res)

    universe <- res$protein_universe
    n_shared <- sum(universe$tier == "shared")
    n_a_only <- sum(universe$tier == "a_only")
    n_b_only <- sum(universe$tier == "b_only")
    n_a <- n_shared + n_a_only
    n_b <- n_shared + n_b_only
    n_total <- n_a + n_b - n_shared
    overlap_pct <- round(n_shared / n_total * 100, 0)

    # Fixed-position Venn with proportional radii
    # Use fixed centers to avoid label collision at high overlap
    r_a <- 1.4
    r_b <- 1.4 * sqrt(n_b / max(n_a, 1))
    r_b <- max(r_b, 0.8)  # minimum visible size

    # Separation: high overlap = closer circles, but with minimum gap for labels
    jaccard <- n_shared / max(n_total, 1)
    sep <- max((r_a + r_b) * (1 - jaccard * 0.7), r_a * 0.6)

    cx_a <- -sep / 2
    cx_b <- sep / 2
    mid_x <- (cx_a + cx_b) / 2

    plotly::plot_ly() |>
      plotly::layout(
        shapes = list(
          list(type = "circle", xref = "x", yref = "y",
               x0 = cx_a - r_a, x1 = cx_a + r_a, y0 = -r_a, y1 = r_a,
               line = list(color = "#3498db", width = 2),
               fillcolor = "rgba(52, 152, 219, 0.15)"),
          list(type = "circle", xref = "x", yref = "y",
               x0 = cx_b - r_b, x1 = cx_b + r_b, y0 = -r_b, y1 = r_b,
               line = list(color = "#e67e22", width = 2),
               fillcolor = "rgba(230, 126, 34, 0.15)")
        ),
        annotations = list(
          # A only count — far left of circle A
          list(x = cx_a - r_a * 0.5, y = 0,
               text = paste0("<b>", n_a_only, "</b><br>A only"),
               showarrow = FALSE, font = list(size = 13, color = "#3498db")),
          # Shared count — center of overlap
          list(x = mid_x, y = 0,
               text = paste0("<b>", n_shared, "</b><br>Shared"),
               showarrow = FALSE, font = list(size = 14, color = "#27ae60")),
          # B only count — far right of circle B
          list(x = cx_b + r_b * 0.5, y = 0,
               text = paste0("<b>", n_b_only, "</b><br>B only"),
               showarrow = FALSE, font = list(size = 13, color = "#e67e22")),
          # Run labels well above circles (no overlap)
          list(x = cx_a - r_a * 0.3, y = r_a + 0.35,
               text = paste0("<b>Run A</b> (", format(n_a, big.mark = ","), ")"),
               showarrow = FALSE, font = list(size = 11, color = "#2c3e50")),
          list(x = cx_b + r_b * 0.3, y = r_b + 0.35,
               text = paste0("<b>Run B</b> (", format(n_b, big.mark = ","), ")"),
               showarrow = FALSE, font = list(size = 11, color = "#2c3e50"))
        ),
        xaxis = list(visible = FALSE, range = list(-3.5, 3.5), fixedrange = TRUE),
        yaxis = list(visible = FALSE, scaleanchor = "x", range = list(-2.2, 2.4), fixedrange = TRUE),
        margin = list(l = 5, r = 5, t = 35, b = 5),
        title = list(
          text = paste0(format(n_total, big.mark = ","), " proteins total \u2014 ",
                        overlap_pct, "% shared"),
          font = list(size = 13))
      ) |>
      plotly::config(displayModeBar = FALSE)
  })

  output$comparator_universe_summary <- renderUI({
    res <- comp_results()
    req(res)
    stats <- res$summary_stats

    div(style = "display: flex; gap: 12px; flex-wrap: wrap; margin-top: 10px;",
      div(class = "card", style = "flex: 1; min-width: 120px; text-align: center; padding: 10px;",
        tags$h4(stats$n_proteins_a, style = "margin: 0; color: #3498db;"),
        tags$small("Run A proteins")
      ),
      div(class = "card", style = "flex: 1; min-width: 120px; text-align: center; padding: 10px;",
        tags$h4(stats$n_proteins_b, style = "margin: 0; color: #e67e22;"),
        tags$small("Run B proteins")
      ),
      div(class = "card", style = "flex: 1; min-width: 120px; text-align: center; padding: 10px;",
        tags$h4(stats$n_shared, style = "margin: 0; color: #2ecc71;"),
        tags$small("Shared")
      ),
      div(class = "card", style = "flex: 1; min-width: 120px; text-align: center; padding: 10px;",
        tags$h4(paste0(round(stats$n_shared / max(stats$n_proteins_a, 1) * 100, 0), "%"),
                style = "margin: 0; color: #6c757d;"),
        tags$small("Overlap (of A)")
      )
    )
  })

  # Protein Universe table with tier filter
  universe_filter <- reactiveVal("All")

  observeEvent(input$universe_filter_all,    universe_filter("All"))
  observeEvent(input$universe_filter_shared, universe_filter("Shared"))
  observeEvent(input$universe_filter_a_only, universe_filter("Run A only"))
  observeEvent(input$universe_filter_b_only, universe_filter("Run B only"))

  output$comparator_universe_table <- DT::renderDT({
    res <- comp_results()
    req(res, res$protein_universe)

    universe <- res$protein_universe
    # Add display-friendly tier label
    universe$Category <- ifelse(universe$tier == "shared", "Shared",
                         ifelse(universe$tier == "a_only", "Run A only", "Run B only"))

    # Apply filter
    filt <- universe_filter()
    if (filt != "All") universe <- universe[universe$Category == filt, ]

    # Build display columns
    display_df <- data.frame(
      Protein = universe$protein_id,
      Category = universe$Category,
      stringsAsFactors = FALSE
    )

    # Add gene symbol column only if we have actual gene names (not blank/accessions)
    if ("gene_symbol" %in% names(universe) && any(nzchar(universe$gene_symbol))) {
      display_df$Gene <- universe$gene_symbol
      display_df <- display_df[, c("Protein", "Gene", "Category")]
    }

    DT::datatable(display_df,
      rownames = FALSE,
      options = list(
        pageLength = 15,
        dom = 'frtip',
        scrollX = TRUE,
        language = list(
          search = "Search proteins:",
          info = "Showing _START_ to _END_ of _TOTAL_ proteins"
        )
      ),
      class = "compact stripe hover"
    ) |>
      DT::formatStyle("Category",
        color = DT::styleEqual(
          c("Shared", "Run A only", "Run B only"),
          c("#2ecc71", "#3498db", "#e67e22")
        ),
        fontWeight = "bold"
      )
  })

  # Export protein universe CSV (respects current filter)
  output$download_universe_csv <- downloadHandler(
    filename = function() {
      filt <- universe_filter()
      suffix <- if (filt == "All") "all" else gsub(" ", "_", tolower(filt))
      paste0("protein_universe_", suffix, ".csv")
    },
    content = function(file) {
      res <- comp_results()
      req(res, res$protein_universe)
      universe <- res$protein_universe
      universe$Category <- ifelse(universe$tier == "shared", "Shared",
                           ifelse(universe$tier == "a_only", "Run A only", "Run B only"))
      filt <- universe_filter()
      if (filt != "All") universe <- universe[universe$Category == filt, ]

      export_df <- data.frame(
        Protein = universe$protein_id,
        Gene = if ("gene_symbol" %in% names(universe)) universe$gene_symbol else "",
        Category = universe$Category,
        In_Run_A = universe$in_a,
        In_Run_B = universe$in_b,
        stringsAsFactors = FALSE
      )
      write.csv(export_df, file, row.names = FALSE)
    }
  )

  # ==========================================================================
  #  LAYER 3: Quantification Comparison
  # ==========================================================================

  # 3a. Global scatter
  output$comparator_quant_scatter <- plotly::renderPlotly({
    res <- comp_results()
    req(res, res$quant_comparison)

    sdata <- res$quant_comparison$scatter_data

    # Color by DE status if concordance available
    if (!is.null(res$de_concordance)) {
      merged <- res$de_concordance$merged
      sdata <- merge(sdata, merged[, c("protein_id", "concordant", "status_a", "status_b")],
                     by = "protein_id", all.x = TRUE)
      sdata$de_label <- ifelse(is.na(sdata$concordant), "No DE data",
                        ifelse(sdata$concordant, "Concordant",
                        ifelse(sdata$status_a != "NS" & sdata$status_b == "NS", "DE in A only",
                        ifelse(sdata$status_a == "NS" & sdata$status_b != "NS", "DE in B only",
                               "Direction change"))))
      color_map <- c("Concordant" = "#95a5a6", "DE in A only" = "#3498db",
                     "DE in B only" = "#e67e22", "Direction change" = "#e74c3c",
                     "No DE data" = "#bdc3c7")
    } else {
      sdata$de_label <- "Shared"
      color_map <- c("Shared" = "#3498db")
    }

    plotly::plot_ly(sdata,
      x = ~mean_a, y = ~mean_b,
      color = ~de_label, colors = color_map,
      text = ~paste0(protein_id, "<br>Run A: ", round(mean_a, 2),
                     "<br>Run B: ", round(mean_b, 2)),
      hoverinfo = "text",
      type = "scatter", mode = "markers",
      marker = list(size = 4, opacity = 0.6)
    ) |>
      plotly::add_trace(
        x = c(min(sdata$mean_a, na.rm = TRUE), max(sdata$mean_a, na.rm = TRUE)),
        y = c(min(sdata$mean_a, na.rm = TRUE), max(sdata$mean_a, na.rm = TRUE)),
        type = "scatter", mode = "lines",
        line = list(color = "grey", dash = "dash", width = 1),
        showlegend = FALSE, inherit = FALSE
      ) |>
      plotly::layout(
        xaxis = list(title = "Mean log2 Intensity (Run A)"),
        yaxis = list(title = "Mean log2 Intensity (Run B)"),
        title = list(text = paste0("Quantification: R = ",
                     round(cor(sdata$mean_a, sdata$mean_b, use = "complete.obs"), 3),
                     " | Offset = ", round(res$quant_comparison$median_offset, 2), " log2"),
                     font = list(size = 13)),
        legend = list(orientation = "h", y = -0.15)
      )
  })

  # 3b. Per-sample correlation (plotly heatmap)
  output$comparator_correlation_heatmap <- plotly::renderPlotly({
    res <- comp_results()
    req(res, res$quant_comparison)

    cor_vals <- res$quant_comparison$cor_per_sample
    if (length(cor_vals) == 0) {
      return(plotly::plot_ly() |>
        plotly::layout(title = "No per-sample correlations available"))
    }

    # Build correlation matrix (samples as rows/cols)
    n <- length(cor_vals)
    labels <- names(cor_vals)
    # Simple display: bar chart of per-sample correlations
    df <- data.frame(
      sample = factor(labels, levels = labels[order(cor_vals)]),
      correlation = cor_vals,
      stringsAsFactors = FALSE
    )

    color_scale <- ifelse(cor_vals < 0.95, "#e74c3c",
                   ifelse(cor_vals < 0.99, "#f39c12", "#2ecc71"))

    plotly::plot_ly(df,
      y = ~sample, x = ~correlation,
      type = "bar", orientation = "h",
      marker = list(color = color_scale),
      text = ~paste0(sample, ": r = ", round(correlation, 4)),
      hoverinfo = "text"
    ) |>
      plotly::layout(
        xaxis = list(title = "Pearson r (Run A vs Run B)", range = c(
          max(0.8, min(cor_vals, na.rm = TRUE) - 0.02), 1.0)),
        yaxis = list(title = ""),
        title = list(text = paste0("Per-Sample Correlations (range: ",
                     round(min(cor_vals, na.rm = TRUE), 3), " - ",
                     round(max(cor_vals, na.rm = TRUE), 3), ")"),
                     font = list(size = 13)),
        margin = list(l = 120)
      )
  })

  # 3c. Bias density
  output$comparator_bias_density <- plotly::renderPlotly({
    res <- comp_results()
    req(res, res$quant_comparison)

    offsets <- res$quant_comparison$offsets
    offsets <- offsets[is.finite(offsets)]

    plotly::plot_ly(x = offsets, type = "histogram",
      nbinsx = 50,
      marker = list(color = "#3498db", line = list(color = "white", width = 0.5)),
      hoverinfo = "x+y"
    ) |>
      plotly::add_trace(
        x = c(res$quant_comparison$median_offset, res$quant_comparison$median_offset),
        y = c(0, length(offsets) / 5),
        type = "scatter", mode = "lines",
        line = list(color = "#e74c3c", width = 2, dash = "dash"),
        showlegend = FALSE, inherit = FALSE
      ) |>
      plotly::layout(
        xaxis = list(title = "log2(Run A) - log2(Run B) per protein"),
        yaxis = list(title = "Count"),
        title = list(
          text = paste0("Intensity Bias Distribution (median offset: ",
                       round(res$quant_comparison$median_offset, 3), " log2)"),
          font = list(size = 13)),
        showlegend = FALSE
      )
  })

  # ==========================================================================
  #  TopN Effect Scatter (Mode B only)
  # ==========================================================================

  output$comparator_topn_effect_section <- renderUI({
    res <- comp_results()
    run_b <- comp_run_b()
    if (is.null(res) || is.null(run_b) || !identical(run_b$source, "spectronaut") ||
        is.null(run_b$n_peptides) || is.null(res$quant_comparison)) return(NULL)

    tagList(
      tags$hr(),
      tags$h6("TopN Effect on Quantification", class = "mt-3"),
      plotly::plotlyOutput("comparator_topn_effect", height = "380px"),
      uiOutput("comparator_topn_interpretation")
    )
  })

  output$comparator_topn_effect <- plotly::renderPlotly({
    res <- comp_results()
    run_b <- comp_run_b()
    req(res, run_b, identical(run_b$source, "spectronaut"),
        run_b$n_peptides, res$quant_comparison)

    qc_data <- res$quant_comparison$scatter_data
    req(nrow(qc_data) > 0)

    topn_val <- suppressWarnings(as.integer(
      gsub("[^0-9]", "", coalesce_setting(run_b$settings$topn_protein, "3"))))
    if (is.na(topn_val)) topn_val <- 3L

    joined <- merge(qc_data, run_b$n_peptides, by = "protein_id", all.x = TRUE)
    joined <- joined[!is.na(joined$n_peptides) & !is.na(joined$mean_a) & !is.na(joined$mean_b), ]
    req(nrow(joined) > 0)

    joined$quant_diff_abs <- abs(joined$mean_a - joined$mean_b)
    joined$topn_limited <- joined$n_peptides > topn_val
    joined$color_label <- ifelse(joined$topn_limited,
      paste0(">", topn_val, " peptides (TopN limiting)"),
      paste0("<=", topn_val, " peptides"))
    joined$hover <- paste0(
      joined$gene %||% joined$protein_id, "\n",
      joined$n_peptides, " peptides\n",
      "Quant diff: ", round(joined$quant_diff_abs, 2), " log2")

    p <- ggplot2::ggplot(joined, ggplot2::aes(
      x = n_peptides, y = quant_diff_abs, color = color_label, text = hover)) +
      ggplot2::geom_point(alpha = 0.4, size = 1.5) +
      ggplot2::geom_vline(xintercept = topn_val + 0.5, linetype = "dashed", color = "black", alpha = 0.7) +
      ggplot2::geom_smooth(ggplot2::aes(x = n_peptides, y = quant_diff_abs),
                           method = "loess", se = TRUE, color = "black",
                           inherit.aes = FALSE, linewidth = 0.8, data = joined) +
      ggplot2::scale_color_manual(
        values = c(setNames("#27AE60", paste0("<=", topn_val, " peptides")),
                   setNames("#E67E22", paste0(">", topn_val, " peptides (TopN limiting)"))),
        name = NULL) +
      ggplot2::scale_x_continuous(trans = "log1p",
                                  breaks = c(1, 2, 3, 5, 10, 20, 50),
                                  labels = c(1, 2, 3, 5, 10, 20, 50)) +
      ggplot2::labs(
        x = "Peptides identified by Spectronaut (log scale)",
        y = "|log2 quant difference| (DE-LIMP vs Spectronaut)",
        title = sprintf("TopN=%d Effect on Quantification", topn_val)) +
      ggplot2::theme_minimal(base_size = 12)

    plotly::ggplotly(p, tooltip = "text")
  })

  output$comparator_topn_interpretation <- renderUI({
    res <- comp_results()
    run_b <- comp_run_b()
    req(res, run_b, identical(run_b$source, "spectronaut"),
        run_b$n_peptides, res$quant_comparison)

    qc_data <- res$quant_comparison$scatter_data
    req(nrow(qc_data) > 0)

    topn_val <- suppressWarnings(as.integer(
      gsub("[^0-9]", "", coalesce_setting(run_b$settings$topn_protein, "3"))))
    if (is.na(topn_val)) topn_val <- 3L

    joined <- merge(qc_data, run_b$n_peptides, by = "protein_id", all.x = TRUE)
    joined <- joined[!is.na(joined$n_peptides) & !is.na(joined$mean_a) & !is.na(joined$mean_b), ]
    req(nrow(joined) > 0)
    joined$quant_diff_abs <- abs(joined$mean_a - joined$mean_b)

    summary <- compute_topn_effect_summary(joined, topn_val)
    if (is.null(summary)) return(NULL)

    alert_class <- if (summary$severity == "warning") "alert-warning" else "alert-info"
    div(class = paste("alert", alert_class, "mt-2"),
      icon(if (summary$severity == "warning") "triangle-exclamation" else "info-circle"),
      " ", summary$text
    )
  })

  # ==========================================================================
  #  Per-Sample QC Section (Mode B — Settings Diff)
  # ==========================================================================

  output$comparator_sample_qc_section <- renderUI({
    run_b <- comp_run_b()
    if (is.null(run_b) || !identical(run_b$source, "spectronaut") || is.null(run_b$run_qc)) return(NULL)

    tagList(
      tags$hr(),
      tags$h6("Per-Sample QC (from Spectronaut RunSummaries)", class = "mt-3"),
      plotly::plotlyOutput("comparator_sample_qc_plot", height = "400px"),
      uiOutput("comparator_sample_qc_outlier"),
      DT::DTOutput("comparator_sample_qc_table")
    )
  })

  output$comparator_sample_qc_plot <- plotly::renderPlotly({
    run_b <- comp_run_b()
    req(run_b, identical(run_b$source, "spectronaut"), run_b$run_qc)

    qc <- run_b$run_qc
    # Use the right column names (flexible naming from parse_spectronaut_run_summaries)
    prec_col <- intersect(c("precursors", "n_precursors"), names(qc))[1]
    pg_col <- intersect(c("protein_groups", "n_protein_groups"), names(qc))[1]
    name_col <- intersect(c("file_name", "run_label"), names(qc))[1]
    req(!is.na(prec_col), !is.na(name_col))

    prec_vals <- as.numeric(qc[[prec_col]])
    qc$run_name <- qc[[name_col]]
    med <- median(prec_vals, na.rm = TRUE)
    qc$outlier <- prec_vals < 0.6 * med
    qc$prec <- prec_vals
    qc$pg <- if (!is.na(pg_col)) as.numeric(qc[[pg_col]]) else NA

    qc$hover <- paste0(
      qc$run_name, "\n",
      format(round(qc$prec), big.mark = ","), " precursors",
      if (!is.na(pg_col)) paste0("\n", format(round(qc$pg), big.mark = ","), " PGs") else "")

    p <- ggplot2::ggplot(qc, ggplot2::aes(
      x = stats::reorder(run_name, prec), y = prec, fill = outlier, text = hover)) +
      ggplot2::geom_col() +
      ggplot2::geom_hline(yintercept = med, linetype = "dashed", color = "black") +
      ggplot2::geom_hline(yintercept = 0.6 * med, linetype = "dotted", color = "red", alpha = 0.5) +
      ggplot2::scale_fill_manual(values = c("FALSE" = "#4A90D9", "TRUE" = "#E74C3C"), guide = "none") +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(labels = scales::comma) +
      ggplot2::labs(x = NULL, y = "Precursors identified",
                    caption = "Dashed = median; dotted red = 60% threshold") +
      ggplot2::theme_minimal(base_size = 12)

    plotly::ggplotly(p, tooltip = "text")
  })

  output$comparator_sample_qc_outlier <- renderUI({
    run_b <- comp_run_b()
    req(run_b, identical(run_b$source, "spectronaut"), run_b$run_qc)

    qc <- run_b$run_qc
    prec_col <- intersect(c("precursors", "n_precursors"), names(qc))[1]
    name_col <- intersect(c("file_name", "run_label"), names(qc))[1]
    req(!is.na(prec_col), !is.na(name_col))

    prec_vals <- as.numeric(qc[[prec_col]])
    med <- median(prec_vals, na.rm = TRUE)
    outliers <- qc[[name_col]][prec_vals < 0.6 * med]

    if (length(outliers) == 0) return(NULL)

    div(class = "alert alert-danger mt-2",
      icon("circle-exclamation"), " ",
      tags$b(sprintf("%d sample(s) appear low-quality: %s.",
                     length(outliers), paste(outliers, collapse = ", "))),
      " These have <60% of the median precursor count and may drive discordant DE results."
    )
  })

  output$comparator_sample_qc_table <- DT::renderDT({
    run_b <- comp_run_b()
    req(run_b, identical(run_b$source, "spectronaut"), run_b$run_qc)
    DT::datatable(run_b$run_qc, rownames = FALSE,
      options = list(pageLength = 25, dom = 'tp', scrollX = TRUE),
      class = "compact stripe"
    )
  })

  # ==========================================================================
  #  LAYER 4: DE Concordance
  # ==========================================================================

  output$comparator_layer4_content <- renderUI({
    res <- comp_results()
    req(res)

    # No DE stats in Run B -> show banner (mode-aware message)
    if (is.null(res$de_concordance)) {
      mode <- input$comparator_mode
      hint <- if (identical(mode, "delimp_spectronaut")) {
        tagList("Ensure the Spectronaut export includes a ",
                tags$b("Candidates.tsv"), " file with DE statistics (Log2Ratio, Pvalue, Qvalue columns).")
      } else if (identical(mode, "delimp_fragpipe")) {
        tagList("Run ", tags$a("FragPipe-Analyst", href = "https://fragpipe-analyst.org/",
                               target = "_blank"),
                " and upload its DE results export to enable this layer.")
      } else {
        "The uploaded file contains intensities only."
      }
      return(div(class = "alert alert-warning mt-3",
        icon("triangle-exclamation"), " ",
        tags$b("DE Concordance requires DE statistics. "), hint
      ))
    }

    conc <- res$de_concordance
    conc_pct <- if (conc$n_total > 0) round(100 * conc$n_concordant / conc$n_total, 1) else NA
    n_disc <- conc$n_discordant

    tagList(
      # Explanation
      div(class = "alert alert-info py-2 px-3 mb-3", style = "font-size: 0.88em;",
        tags$b("How to read this:"), " Each protein is classified as ",
        tags$b("Up"), " (significant & positive logFC), ",
        tags$b("Down"), " (significant & negative logFC), or ",
        tags$b("NS"), " (Not Significant, adj. P \u2265 0.05) in each run. ",
        "The matrix shows how many proteins fall into each combination. ",
        "Diagonal cells (Up/Up, Down/Down, NS/NS) are ", tags$b("concordant"), ". ",
        "Off-diagonal cells are ", tags$b("discordant"), " \u2014 proteins where the two runs disagree.",
        if (!is.na(conc_pct)) tagList(
          tags$br(),
          sprintf("Overall concordance: %s%% (%d discordant proteins).", conc_pct, n_disc)
        )
      ),

      # 3x3 concordance matrix + hypothesis chart
      div(style = "display: flex; gap: 20px; flex-wrap: wrap; margin-bottom: 15px;",
        div(style = "flex: 1; min-width: 300px;",
          tags$h6("Concordance Matrix (Run A rows \u00d7 Run B columns)"),
          DT::DTOutput("comparator_concordance_3x3", height = "auto")
        ),
        div(style = "flex: 1; min-width: 300px;",
          div(style = "display: flex; align-items: center; gap: 6px;",
            tags$h6("Hypothesis Distribution", style = "margin: 0;"),
            actionButton("comparator_hypothesis_info_btn", icon("question-circle"),
                         class = "btn-outline-info btn-sm", title = "What do these categories mean?",
                         style = "padding: 1px 5px; font-size: 0.75em;")
          ),
          tags$p(style = "font-size: 0.82em; color: #666; margin-bottom: 4px;",
            "Why do discordant proteins disagree? Each is assigned a likely cause."),
          plotly::plotlyOutput("comparator_hypothesis_dist", height = "200px")
        )
      ),

      # Volcano overlay
      plotly::plotlyOutput("comparator_volcano_overlay", height = "450px"),

      # Discordant protein table
      tags$h6("Discordant Proteins", class = "mt-3"),
      tags$p(style = "font-size: 0.82em; color: #666; margin-top: -4px;",
        "Proteins classified differently between runs. The hypothesis column suggests why they disagree."),
      DT::DTOutput("comparator_discordant_table")
    )
  })

  # 3x3 concordance matrix
  output$comparator_concordance_3x3 <- DT::renderDT({
    res <- comp_results()
    req(res, res$de_concordance)

    mat <- res$de_concordance$concordance_3x3
    df <- as.data.frame.matrix(mat)
    df$`Run A` <- rownames(df)
    df <- df[, c("Run A", "Up", "Down", "NS")]

    DT::datatable(df,
      rownames = FALSE,
      options = list(dom = 't', ordering = FALSE, pageLength = 3),
      class = "compact"
    ) |>
      DT::formatStyle(c("Up", "Down", "NS"),
        backgroundColor = DT::styleInterval(
          c(0, 5, 20),
          c("transparent", "#fff3cd", "#ffeeba", "#f8d7da")
        )
      )
  })

  # Hypothesis distribution chart
  output$comparator_hypothesis_dist <- plotly::renderPlotly({
    res <- comp_results()
    req(res, res$de_concordance)
    disc <- res$de_concordance$discordant_table
    if (nrow(disc) == 0) {
      return(plotly::plot_ly() |>
        plotly::layout(title = "No discordant proteins"))
    }

    cats <- table(disc$hypothesis_category)
    df <- data.frame(
      category = names(cats),
      count = as.integer(cats),
      stringsAsFactors = FALSE
    )
    df <- df[order(-df$count), ]
    df$category <- factor(df$category, levels = df$category)

    cat_colors <- c("Direction reversal" = "#e74c3c", "Normalization offset" = "#f39c12",
                    "Variance estimation" = "#9b59b6", "Missing values" = "#3498db",
                    "Peptide count" = "#1abc9c", "FC magnitude" = "#e67e22",
                    "Borderline" = "#95a5a6")
    colors <- cat_colors[as.character(df$category)]
    colors[is.na(colors)] <- "#95a5a6"

    plotly::plot_ly(df,
      x = ~category, y = ~count,
      type = "bar",
      marker = list(color = colors),
      text = ~count, textposition = "outside",
      hoverinfo = "x+y"
    ) |>
      plotly::layout(
        xaxis = list(title = ""),
        yaxis = list(title = "Count"),
        margin = list(b = 80),
        title = list(text = paste0(sum(df$count), " discordant proteins"),
                     font = list(size = 12))
      ) |>
      plotly::config(displayModeBar = FALSE)
  })

  # Volcano overlay
  output$comparator_volcano_overlay <- plotly::renderPlotly({
    res <- comp_results()
    req(res, res$de_concordance)
    merged <- res$de_concordance$merged

    # Build volcano data
    vol_a <- data.frame(
      protein_id = merged$protein_id,
      logFC      = merged$logFC_A,
      neg_log10p = -log10(merged$adjP_A + 1e-300),
      run        = "Run A",
      concordant = merged$concordant,
      stringsAsFactors = FALSE
    )
    vol_b <- data.frame(
      protein_id = merged$protein_id,
      logFC      = merged$logFC_B,
      neg_log10p = -log10(merged$adjP_B + 1e-300),
      run        = "Run B",
      concordant = merged$concordant,
      stringsAsFactors = FALSE
    )

    plotly::plot_ly() |>
      plotly::add_trace(data = vol_a,
        x = ~logFC, y = ~neg_log10p,
        color = I(ifelse(vol_a$concordant, "#3498db40", "#3498dbCC")),
        text = ~paste0(protein_id, " (A)<br>logFC: ", round(logFC, 3),
                       "<br>-log10(adj.P): ", round(neg_log10p, 2)),
        hoverinfo = "text", type = "scatter", mode = "markers",
        marker = list(size = 4), name = "Run A"
      ) |>
      plotly::add_trace(data = vol_b,
        x = ~logFC, y = ~neg_log10p,
        color = I(ifelse(vol_b$concordant, "#e67e2240", "#e67e22CC")),
        text = ~paste0(protein_id, " (B)<br>logFC: ", round(logFC, 3),
                       "<br>-log10(adj.P): ", round(neg_log10p, 2)),
        hoverinfo = "text", type = "scatter", mode = "markers",
        marker = list(size = 4), name = "Run B"
      ) |>
      plotly::layout(
        xaxis = list(title = "logFC"),
        yaxis = list(title = "-log10(adj.P.Val)"),
        title = list(text = "Volcano Overlay (Run A blue, Run B orange)",
                     font = list(size = 13)),
        shapes = list(
          list(type = "line", x0 = -10, x1 = 10,
               y0 = -log10(0.05), y1 = -log10(0.05),
               line = list(color = "grey", dash = "dash", width = 1))
        ),
        legend = list(orientation = "h", y = -0.12)
      )
  })

  # Discordant protein table
  output$comparator_discordant_table <- DT::renderDT({
    res <- comp_results()
    req(res, res$de_concordance)
    disc <- res$de_concordance$discordant_table

    if (nrow(disc) == 0) {
      return(DT::datatable(data.frame(Message = "No discordant proteins found")))
    }

    display_cols <- c("protein_id", "status_a", "status_b",
                      "logFC_A", "logFC_B", "adjP_A", "adjP_B",
                      "hypothesis_category", "confidence")
    available <- intersect(display_cols, names(disc))
    show_df <- disc[, available, drop = FALSE]

    # Round numeric columns
    for (col in c("logFC_A", "logFC_B")) {
      if (col %in% names(show_df)) show_df[[col]] <- round(show_df[[col]], 3)
    }
    for (col in c("adjP_A", "adjP_B")) {
      if (col %in% names(show_df)) show_df[[col]] <- signif(show_df[[col]], 3)
    }

    DT::datatable(show_df,
      rownames = FALSE,
      colnames = c("Protein", "Status A", "Status B", "logFC A", "logFC B",
                    "adj.P A", "adj.P B", "Diagnosis", "Confidence")[seq_along(available)],
      options = list(pageLength = 20, scrollX = TRUE,
                     order = list(list(which(available == "confidence") - 1, "asc"))),
      class = "compact stripe hover",
      selection = "single"
    ) |>
      DT::formatStyle("confidence",
        backgroundColor = DT::styleEqual(
          c("High", "Medium", "Low"),
          c("#f8d7da", "#fff3cd", "#f8f9fa")
        )
      )
  })

  # ==========================================================================
  #  AI ANALYSIS TAB
  # ==========================================================================

  # Gemini summary
  observeEvent(input$comparator_gemini_btn, {
    res <- comp_results()
    req(res, input$user_api_key)

    mofa_obj <- values$comparator_mofa
    prompt <- build_gemini_comparator_prompt(res, mofa_obj, values$instrument_metadata)

    withProgress(message = "Generating Gemini summary...", value = 0.3, {
      tryCatch({
        response <- ask_gemini_text_chat(prompt, input$user_api_key,
                                         input$model_name %||% "gemini-2.0-flash")
        setProgress(1.0, detail = "Done")
        values$comparator_gemini_narrative <- response
      }, error = function(e) {
        showNotification(paste("Gemini error:", e$message), type = "error")
      })
    })
  })

  # --- View Prompt (transparency) ---
  observeEvent(input$comparator_view_prompt_btn, {
    res <- comp_results()
    req(res)

    mofa_obj <- values$comparator_mofa
    prompt <- build_gemini_comparator_prompt(res, mofa_obj, values$instrument_metadata)

    # JS to copy prompt text to clipboard
    copy_js <- paste0(
      "navigator.clipboard.writeText(document.getElementById('comparator_prompt_text').innerText)",
      ".then(function(){ Shiny.setInputValue('comparator_prompt_copied', Math.random()); })"
    )

    showModal(modalDialog(
      title = "Gemini Comparator Prompt",
      size = "l",
      easyClose = TRUE,
      tags$p(class = "text-muted small",
        "This is the exact prompt sent to Gemini. Review for bias or missing context. ",
        "You can copy it to use with any AI tool."),
      tags$div(style = "margin-bottom: 10px;",
        actionButton("comparator_copy_prompt_btn", "Copy to Clipboard",
                      icon = icon("clipboard"), class = "btn-outline-secondary btn-sm",
                      onclick = copy_js)
      ),
      tags$pre(id = "comparator_prompt_text",
               style = "max-height: 500px; overflow-y: auto; white-space: pre-wrap; word-wrap: break-word; font-size: 0.85em; background: #f8f9fa; padding: 12px; border-radius: 6px;",
               prompt),
      footer = modalButton("Close")
    ))
  })

  # Clipboard copy confirmation
  observeEvent(input$comparator_prompt_copied, {
    showNotification("Prompt copied to clipboard", type = "message", duration = 2)
  })

  # Inject Gemini narrative into static div (avoids bslib uiOutput disappearing bug)
  observe({
    narrative <- values$comparator_gemini_narrative
    if (is.null(narrative)) {
      shinyjs::html("comparator_gemini_container", "")
    } else {
      html_content <- paste0(
        '<div class="card mt-3"><div class="card-header"><b>Gemini Analysis</b></div>',
        '<div class="card-body" style="max-height:500px;overflow-y:auto;">',
        markdown::markdownToHTML(text = narrative, fragment.only = TRUE),
        '</div></div>'
      )
      shinyjs::html("comparator_gemini_container", html_content)
    }
  })

  # --- MOFA2 Decomposition (optional) ---

  # MOFA2 run — uses callr::r() subprocess to isolate basilisk/Python from Shiny
  observeEvent(input$comparator_mofa_btn, {
    res <- comp_results()
    req(res, res$quant_comparison)

    quant <- res$quant_comparison
    n_proteins <- nrow(quant$mat_a)

    if (n_proteins < 10) {
      showNotification("Too few shared proteins for MOFA2 (need >= 10).", type = "warning")
      return()
    }
    if (is.null(res$sample_map)) {
      showNotification("No sample mapping available for MOFA2.", type = "warning")
      return()
    }

    # Prepare data before launching subprocess
    matched <- res$sample_map[res$sample_map$status == "matched", ]
    if (nrow(matched) < 4) {
      showNotification("Need >= 4 matched sample pairs for MOFA2.", type = "warning")
      return()
    }

    # Subset to matched columns and assign common names
    cols_a <- match(matched$run_a, colnames(quant$mat_a))
    cols_b <- match(matched$run_b, colnames(quant$mat_b))
    valid <- !is.na(cols_a) & !is.na(cols_b)
    cols_a <- cols_a[valid]
    cols_b <- cols_b[valid]
    common_names <- paste0("Sample_", seq_along(cols_a))

    view_a <- quant$mat_a[, cols_a, drop = FALSE]
    view_b <- quant$mat_b[, cols_b, drop = FALSE]
    colnames(view_a) <- common_names
    colnames(view_b) <- common_names

    # Subset to top 500 most variable proteins to keep MOFA2 fast
    max_features <- 500
    if (nrow(view_a) > max_features) {
      # Combined variance across both views
      var_a <- apply(view_a, 1, var, na.rm = TRUE)
      var_b <- apply(view_b, 1, var, na.rm = TRUE)
      combined_var <- var_a + var_b
      top_idx <- order(combined_var, decreasing = TRUE)[seq_len(max_features)]
      view_a <- view_a[top_idx, , drop = FALSE]
      view_b <- view_b[top_idx, , drop = FALSE]
    }

    # Make feature names unique across views to avoid MOFA2 warning
    rownames(view_a) <- paste0("A_", rownames(view_a) %||% seq_len(nrow(view_a)))
    rownames(view_b) <- paste0("B_", rownames(view_b) %||% seq_len(nrow(view_b)))

    n_factors <- min(5, length(common_names) - 1)
    mofa_data <- list(run_a = view_a, run_b = view_b)

    withProgress(message = "Running MOFA2 in subprocess...", value = 0.2, {
      tryCatch({
        outfile <- tempfile(fileext = ".hdf5")

        setProgress(0.3, detail = "Training model (this may take 1-2 minutes)...")

        # Run in isolated subprocess — if basilisk/Python crashes, only the
        # subprocess dies, not the Shiny app (same pattern as server_mofa.R)
        callr::r(
          function(mofa_data, outfile, n_factors) {
            library(MOFA2)

            mofa_obj <- create_mofa(mofa_data)

            data_opts <- get_default_data_options(mofa_obj)
            data_opts$scale_views <- TRUE

            model_opts <- get_default_model_options(mofa_obj)
            model_opts$num_factors <- n_factors

            train_opts <- get_default_training_options(mofa_obj)
            train_opts$convergence_mode <- "fast"
            train_opts$seed <- 42
            train_opts$verbose <- FALSE

            mofa_obj <- prepare_mofa(mofa_obj,
              data_options = data_opts,
              model_options = model_opts,
              training_options = train_opts
            )

            run_mofa(mofa_obj, outfile = outfile, use_basilisk = TRUE)
          },
          args = list(
            mofa_data = mofa_data,
            outfile = outfile,
            n_factors = n_factors
          ),
          show = FALSE,
          timeout = 300  # 5 minute timeout
        )

        setProgress(0.8, detail = "Loading trained model...")

        # Load the trained model back into the main process
        mofa_obj <- MOFA2::load_model(outfile)
        unlink(outfile)

        values$comparator_mofa <- mofa_obj
        shinyjs::hide("comparator_mofa_btn")
        showNotification("MOFA2 decomposition complete", type = "message")
        setProgress(1.0, detail = "Done")
      }, error = function(e) {
        msg <- conditionMessage(e)
        if (!is.null(e$parent)) msg <- conditionMessage(e$parent)
        showNotification(paste("MOFA2 failed:", msg),
                         type = "error", duration = 10)
      })
    })
  })

  # MOFA2 variance explained (plotly heatmap)
  output$comparator_mofa_variance <- plotly::renderPlotly({
    req(values$comparator_mofa)
    mofa_obj <- values$comparator_mofa

    var_exp <- MOFA2::get_variance_explained(mofa_obj)
    r2 <- var_exp$r2_per_factor[[1]]  # matrix: factors x views

    plotly::plot_ly(
      z = t(r2) * 100,
      x = rownames(r2),
      y = colnames(r2),
      type = "heatmap",
      colorscale = list(c(0, "#f7fbff"), c(0.5, "#6baed6"), c(1, "#08306b")),
      text = round(t(r2) * 100, 1),
      texttemplate = "%{text}%",
      hoverinfo = "x+y+z"
    ) |>
      plotly::layout(
        xaxis = list(title = "Factor"),
        yaxis = list(title = ""),
        title = list(
          text = "Variance Explained (%) per Run and Factor",
          font = list(size = 13)
        )
      )
  })

  # MOFA2 factor weights scatter
  output$comparator_mofa_weights <- plotly::renderPlotly({
    req(values$comparator_mofa, comp_results())
    mofa_obj <- values$comparator_mofa

    disc_ids <- if (!is.null(comp_results()$de_concordance)) {
      comp_results()$de_concordance$discordant_table$protein_id
    } else character(0)

    var_exp <- MOFA2::get_variance_explained(mofa_obj)
    r2 <- var_exp$r2_per_factor[[1]]
    top_factors <- rownames(r2)[order(-rowSums(r2))][1:min(2, nrow(r2))]

    weights <- MOFA2::get_weights(mofa_obj,
                                   views = "run_a",
                                   factors = top_factors,
                                   as.data.frame = TRUE)
    names(weights)[names(weights) == "feature"] <- "protein_id"
    weights_wide <- tidyr::pivot_wider(weights, names_from = factor, values_from = value)

    weights_wide$is_discordant <- weights_wide$protein_id %in% disc_ids
    weights_wide$color_group <- ifelse(weights_wide$is_discordant, "Discordant", "Concordant")

    f1 <- top_factors[1]
    f2 <- if (length(top_factors) >= 2) top_factors[2] else top_factors[1]

    plotly::plot_ly(weights_wide,
      x = as.formula(paste0("~`", f1, "`")),
      y = as.formula(paste0("~`", f2, "`")),
      color = ~color_group,
      colors = c("Concordant" = "grey80", "Discordant" = "#e74c3c"),
      text = ~protein_id,
      hoverinfo = "text",
      type = "scatter", mode = "markers",
      marker = list(size = 5, opacity = 0.7)
    ) |>
      plotly::layout(
        xaxis = list(title = f1),
        yaxis = list(title = f2),
        title = list(text = "Factor weights - discordant proteins highlighted",
                     font = list(size = 13))
      )
  })

  # MOFA2 top weights table
  output$comparator_mofa_top_weights <- DT::renderDT({
    req(values$comparator_mofa, comp_results())
    mofa_obj <- values$comparator_mofa

    var_exp <- MOFA2::get_variance_explained(mofa_obj)
    r2 <- var_exp$r2_per_factor[[1]]
    r2_a <- r2[, "run_a"]
    r2_b <- r2[, "run_b"]
    specificity <- r2_a - r2_b
    key_factor <- names(which.max(specificity))

    disc_ids <- if (!is.null(comp_results()$de_concordance)) {
      comp_results()$de_concordance$discordant_table$protein_id
    } else character(0)

    weights <- MOFA2::get_weights(mofa_obj,
                                   views = "run_a",
                                   factors = key_factor,
                                   as.data.frame = TRUE)
    names(weights)[names(weights) == "feature"] <- "protein_id"
    names(weights)[names(weights) == "value"] <- "weight"

    weights$abs_weight <- abs(weights$weight)
    weights$is_discordant <- weights$protein_id %in% disc_ids
    weights <- weights[order(-weights$abs_weight), ]
    weights <- head(weights, 15)
    weights$weight <- round(weights$weight, 4)
    weights$abs_weight <- round(weights$abs_weight, 4)

    DT::datatable(
      weights[, c("protein_id", "weight", "abs_weight", "is_discordant")],
      caption = paste0("Top proteins driving ", key_factor,
                       " (most run-A-specific factor; explains ",
                       round(r2_a[key_factor] * 100, 1), "% in A vs ",
                       round(r2_b[key_factor] * 100, 1), "% in B)"),
      rownames = FALSE,
      colnames = c("Protein", "Weight", "|Weight|", "Discordant?"),
      options = list(pageLength = 15, dom = 't'),
      class = "compact"
    ) |>
      DT::formatStyle("is_discordant",
        target = "row",
        backgroundColor = DT::styleEqual(TRUE, "#fde8e8")
      )
  })

  # ==========================================================================
  #  CLAUDE ZIP EXPORT
  # ==========================================================================

  output$comparator_claude_export <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_Run_Comparison_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip")
    },
    content = function(file) {
      res <- comp_results()
      req(res)

      tryCatch({
        withProgress(message = "Building export...", value = 0.1, {
          tmp_dir <- file.path(tempdir(), paste0("comparator_export_", format(Sys.time(), "%Y%m%d_%H%M%S")))
          dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
          files_to_zip <- character(0)

          # 1. Settings diff CSV
          setProgress(0.2, detail = "Writing settings...")
          settings_file <- file.path(tmp_dir, "settings_diff.csv")
          write.csv(res$settings_diff, settings_file, row.names = FALSE)
          files_to_zip <- c(files_to_zip, settings_file)

          # 2. Protein universe CSV
          setProgress(0.3, detail = "Writing protein universe...")
          universe_file <- file.path(tmp_dir, "protein_universe.csv")
          write.csv(res$protein_universe, universe_file, row.names = FALSE)
          files_to_zip <- c(files_to_zip, universe_file)

          # 3. DE results combined (if available)
          if (!is.null(res$de_concordance)) {
            setProgress(0.4, detail = "Writing DE results...")

            merged_file <- file.path(tmp_dir, "de_results_combined.csv")
            merged_export <- res$de_concordance$merged[, c("protein_id", "logFC_A", "logFC_B",
                                                            "adjP_A", "adjP_B",
                                                            "status_a", "status_b", "concordant")]
            write.csv(merged_export, merged_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, merged_file)

            # 4. Discordant proteins CSV
            disc_file <- file.path(tmp_dir, "discordant_proteins.csv")
            disc_export <- res$de_concordance$discordant_table
            export_cols <- intersect(
              c("protein_id", "logFC_A", "logFC_B", "adjP_A", "adjP_B",
                "status_a", "status_b", "n_peptides_A", "n_peptides_B",
                "n_ratios_B",
                "missing_pct_A", "missing_pct_B", "hypothesis", "hypothesis_category", "confidence"),
              names(disc_export)
            )
            write.csv(disc_export[, export_cols], disc_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, disc_file)
          }

          # 4b. Spectronaut-specific files (Mode B)
          run_b_data <- comp_run_b()
          if (!is.null(run_b_data) && identical(run_b_data$source, "spectronaut")) {
            if (!is.null(run_b_data$run_qc)) {
              qc_file <- file.path(tmp_dir, "spectronaut_run_qc.csv")
              write.csv(run_b_data$run_qc, qc_file, row.names = FALSE)
              files_to_zip <- c(files_to_zip, qc_file)
            }
            if (!is.null(run_b_data$library_info)) {
              li_df <- data.frame(
                field = names(run_b_data$library_info),
                value = as.character(unlist(run_b_data$library_info)),
                stringsAsFactors = FALSE
              )
              li_file <- file.path(tmp_dir, "spectronaut_library_info.csv")
              write.csv(li_df, li_file, row.names = FALSE)
              files_to_zip <- c(files_to_zip, li_file)
            }
          }

          # 5. Comparison context
          setProgress(0.6, detail = "Writing context...")
          context_lines <- c(
            "# Run Comparison Context",
            "",
            paste("Run A:", res$summary_stats$run_a_label),
            paste("Run B:", res$summary_stats$run_b_label),
            paste("Contrast:", res$summary_stats$contrast),
            paste("Samples:", res$summary_stats$n_samples),
            paste("Shared proteins:", res$summary_stats$n_shared),
            paste("Global offset:", round(res$summary_stats$global_offset, 3), "log2 units"),
            ""
          )

          # Include Run A settings from session
          run_a_data <- comp_run_a()
          if (!is.null(run_a_data$settings)) {
            context_lines <- c(context_lines, "## Run A Settings", "")
            for (nm in names(run_a_data$settings)) {
              val <- run_a_data$settings[[nm]]
              val_str <- tryCatch(as.character(val), error = function(e) NULL)
              if (!is.null(val_str) && length(val_str) == 1 && !is.na(val_str) &&
                  nzchar(val_str) && val_str != "Unknown")
                context_lines <- c(context_lines, paste0("- **", nm, "**: ", val_str))
            }
            context_lines <- c(context_lines, "")
          }

          # Include Run B settings
          if (!is.null(run_b_data) && !is.null(run_b_data$settings)) {
            context_lines <- c(context_lines, "## Run B Settings", "")
            for (nm in names(run_b_data$settings)) {
              val <- run_b_data$settings[[nm]]
              val_str <- tryCatch(as.character(val), error = function(e) NULL)
              if (!is.null(val_str) && length(val_str) == 1 && !is.na(val_str) &&
                  nzchar(val_str) && val_str != "Unknown")
                context_lines <- c(context_lines, paste0("- **", nm, "**: ", val_str))
            }
            context_lines <- c(context_lines, "")
          }

          # Include concordance summary
          if (!is.null(res$de_concordance)) {
            conc <- res$de_concordance
            conc_pct <- if (conc$n_total > 0) round(100 * conc$n_concordant / conc$n_total, 1) else NA
            context_lines <- c(context_lines,
              "## DE Concordance Summary", "",
              paste("Total shared proteins tested:", conc$n_total),
              paste("Concordant:", conc$n_concordant, sprintf("(%.1f%%)", conc_pct)),
              paste("Discordant:", conc$n_discordant),
              paste("Run A significant:", sum(conc$merged$status_a != "NS")),
              paste("Run B significant:", sum(conc$merged$status_b != "NS")),
              ""
            )
          }

          if (!is.null(values$comparator_gemini_narrative)) {
            context_lines <- c(context_lines,
              "## Gemini Analysis", "",
              values$comparator_gemini_narrative
            )
          }
          context_file <- file.path(tmp_dir, "comparison_context.md")
          writeLines(context_lines, context_file)
          files_to_zip <- c(files_to_zip, context_file)

          # 5b. Excluded files (if any)
          if (!is.null(values$excluded_files) && nrow(values$excluded_files) > 0) {
            excl_file <- file.path(tmp_dir, "excluded_files.csv")
            write.csv(values$excluded_files, excl_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, excl_file)
          }

          # 6. Claude prompt
          setProgress(0.8, detail = "Building prompt...")
          prompt <- build_claude_comparator_prompt(res, values$comparator_gemini_narrative, values$instrument_metadata)
          prompt_file <- file.path(tmp_dir, "claude_prompt.md")
          writeLines(prompt, prompt_file)
          files_to_zip <- c(files_to_zip, prompt_file)

          # 7. DIA-NN search provenance (graceful — skip if unavailable)
          setProgress(0.85, detail = "Extracting search provenance...")
          provenance_notes <- character(0)  # track what was/wasn't available

          # 7a. Search params from session
          tryCatch({
            ss <- values$diann_search_settings
            sp <- ss$search_params
            param_lines <- c("# DIA-NN Search Parameters (Run A — current session)", "")

            if (!is.null(ss)) {
              if (!is.null(ss$search_name))   param_lines <- c(param_lines, paste("Search name:", ss$search_name))
              if (!is.null(ss$output_dir))     param_lines <- c(param_lines, paste("Output dir:", ss$output_dir))
              if (!is.null(ss$search_mode))    param_lines <- c(param_lines, paste("Search mode:", ss$search_mode))
              if (length(ss$fasta_files) > 0)  param_lines <- c(param_lines, paste("FASTA:", paste(basename(ss$fasta_files), collapse = ", ")))
              if (!is.null(ss$fasta_seq_count)) param_lines <- c(param_lines, paste("FASTA sequences:", ss$fasta_seq_count))
            }

            if (!is.null(sp)) {
              if (!is.null(sp$mass_acc_mode))  param_lines <- c(param_lines, paste("Mass accuracy mode:", sp$mass_acc_mode))
              if (!is.null(sp$mass_acc))       param_lines <- c(param_lines, paste("MS2 mass accuracy:", sp$mass_acc))
              if (!is.null(sp$mass_acc_ms1))   param_lines <- c(param_lines, paste("MS1 mass accuracy:", sp$mass_acc_ms1))
              if (!is.null(sp$scan_window))    param_lines <- c(param_lines, paste("Scan window:", sp$scan_window))
              if (!is.null(sp$enzyme))         param_lines <- c(param_lines, paste("Enzyme:", sp$enzyme))
              if (!is.null(sp$mbr))            param_lines <- c(param_lines, paste("MBR:", sp$mbr))
              if (!is.null(sp$min_pr_mz))      param_lines <- c(param_lines, paste("Min precursor m/z:", sp$min_pr_mz))
              if (!is.null(sp$max_pr_mz))      param_lines <- c(param_lines, paste("Max precursor m/z:", sp$max_pr_mz))
              if (!is.null(sp$min_pr_charge))   param_lines <- c(param_lines, paste("Min precursor charge:", sp$min_pr_charge))
              if (!is.null(sp$max_pr_charge))   param_lines <- c(param_lines, paste("Max precursor charge:", sp$max_pr_charge))
              if (nzchar(sp$extra_cli_flags %||% "")) param_lines <- c(param_lines, paste("Extra CLI flags:", sp$extra_cli_flags))
            }

            # Library info — check multiple possible locations
            lib_path <- sp$library %||% ss$library_path %||% values$library_path %||% NULL
            if (!is.null(lib_path)) {
              param_lines <- c(param_lines, "", paste("Spectral library:", basename(lib_path)))
            }

            if (length(param_lines) > 2) {
              params_file <- file.path(tmp_dir, "diann_search_params.txt")
              writeLines(param_lines, params_file)
              files_to_zip <- c(files_to_zip, params_file)
              provenance_notes <- c(provenance_notes, "- `diann_search_params.txt` - DIA-NN search parameters for Run A")
            } else {
              provenance_notes <- c(provenance_notes, "- DIA-NN search parameters: unavailable (no search settings in session)")
            }
          }, error = function(e) {
            provenance_notes <<- c(provenance_notes, paste0("- DIA-NN search parameters: error extracting (", e$message, ")"))
          })

          # 7b. Precursor summary for discordant proteins
          tryCatch({
            rd <- values$raw_data
            disc_ids <- if (!is.null(res$de_concordance)) res$de_concordance$discordant_table$protein_id else character(0)

            if (!is.null(rd) && !is.null(rd$genes) && length(disc_ids) > 0 &&
                "Precursor.Mz" %in% colnames(rd$genes) && "Precursor.Charge" %in% colnames(rd$genes)) {

              genes_df <- rd$genes
              # Normalize protein IDs for matching
              if ("Protein.Group" %in% colnames(genes_df)) {
                genes_df$norm_id <- normalize_protein_id(genes_df$Protein.Group)
              } else {
                genes_df$norm_id <- normalize_protein_id(rownames(rd$E))
              }

              # Filter to discordant proteins
              disc_genes <- genes_df[genes_df$norm_id %in% disc_ids, ]

              if (nrow(disc_genes) > 0) {
                # Summarize per protein: precursor count, charge states, m/z range
                prec_summary <- do.call(rbind, lapply(split(disc_genes, disc_genes$norm_id), function(pg) {
                  data.frame(
                    Protein.Group = pg$norm_id[1],
                    n_precursors  = nrow(pg),
                    charge_states = paste(sort(unique(pg$Precursor.Charge)), collapse = ","),
                    mz_min        = round(min(pg$Precursor.Mz, na.rm = TRUE), 2),
                    mz_max        = round(max(pg$Precursor.Mz, na.rm = TRUE), 2),
                    mz_mean       = round(mean(pg$Precursor.Mz, na.rm = TRUE), 2),
                    stringsAsFactors = FALSE
                  )
                }))

                prec_file <- file.path(tmp_dir, "precursor_summary_discordant.csv")
                write.csv(prec_summary, prec_file, row.names = FALSE)
                files_to_zip <- c(files_to_zip, prec_file)
                provenance_notes <- c(provenance_notes,
                  paste0("- `precursor_summary_discordant.csv` - Precursor-level summary for ",
                         nrow(prec_summary), " discordant proteins"))
              } else {
                provenance_notes <- c(provenance_notes, "- Precursor summary: no discordant proteins matched in raw data")
              }
            } else {
              reason <- if (is.null(rd)) "raw_data not loaded"
                else if (is.null(rd$genes)) "no genes annotation"
                else if (length(disc_ids) == 0) "no discordant proteins"
                else "Precursor.Mz/Charge columns not in raw data"
              provenance_notes <- c(provenance_notes, paste0("- Precursor summary: unavailable (", reason, ")"))
            }
          }, error = function(e) {
            provenance_notes <<- c(provenance_notes, paste0("- Precursor summary: error (", e$message, ")"))
          })

          # 7c. DIA-NN log-derived parameter files
          for (run_id in c("a", "b")) {
            log_obj <- res[[paste0("diann_log_", run_id)]]
            if (!is.null(log_obj) && isTRUE(log_obj$success)) {
              run_label <- paste("Run", toupper(run_id))
              log_file <- file.path(tmp_dir, paste0("diann_params_run_", run_id, ".txt"))
              p <- log_obj$params
              log_lines <- c(
                paste0("DIA-NN Search Parameters - ", run_label),
                paste0(rep("=", 50), collapse = ""),
                paste0("DIA-NN Version:        ", log_obj$version %||% "unknown"),
                paste0("Pipeline Step:         ", log_obj$pipeline_step %||% "full search"),
                "",
                "--- Library ---",
                paste0("Library Source:        ", if (!is.null(log_obj$lib_path)) basename(log_obj$lib_path) else "(in silico from FASTA)"),
                paste0("FASTA File(s):         ", paste(basename(log_obj$fasta_files %||% "unknown"), collapse = "; ")),
                paste0("Precursors in Library: ", if (!is.null(log_obj$n_precursors_library)) format(log_obj$n_precursors_library, big.mark = ",") else "unknown"),
                paste0("Output Library:        ", if (!is.null(log_obj$out_lib_path)) basename(log_obj$out_lib_path) else "none"),
                "",
                "--- Protein Grouping ---",
                paste0("--pg-level:            ", p$pg_level %||% "unknown"),
                paste0("--proteoforms:         ", if (isTRUE(p$proteoforms)) "yes" else "no"),
                "",
                "--- Search Parameters ---",
                paste0("--window:              ", p$scan_window %||% "unknown"),
                paste0("--mass-acc (MS2):      ", p$mass_acc %||% "unknown"),
                paste0("--mass-acc-ms1:        ", p$mass_acc_ms1 %||% "unknown"),
                paste0("--min-pr-mz:           ", p$min_pr_mz %||% "unknown"),
                paste0("--max-pr-mz:           ", p$max_pr_mz %||% "unknown"),
                paste0("--min-pr-charge:       ", p$min_pr_charge %||% "unknown"),
                paste0("--max-pr-charge:       ", p$max_pr_charge %||% "unknown"),
                paste0("--min-fr-mz:           ", p$min_fr_mz %||% "unknown"),
                paste0("--max-fr-mz:           ", p$max_fr_mz %||% "unknown"),
                paste0("--min-pep-len:         ", p$min_pep_len %||% "unknown"),
                paste0("--max-pep-len:         ", p$max_pep_len %||% "unknown"),
                paste0("--missed-cleavages:    ", p$missed_cleavages %||% "unknown"),
                paste0("--enzyme (--cut):      ", p$enzyme %||% "unknown"),
                paste0("--qvalue:              ", p$qvalue %||% "unknown"),
                paste0("--reanalyse (MBR):     ", if (isTRUE(p$mbr)) "yes" else "no"),
                paste0("--fasta-search:        ", if (isTRUE(log_obj$search_mode == "libfree")) "yes" else "no")
              )
              writeLines(log_lines, log_file)
              files_to_zip <- c(files_to_zip, log_file)
              provenance_notes <- c(provenance_notes,
                paste0("- `diann_params_run_", run_id, ".txt` - Parsed DIA-NN parameters from ", run_label, " log"))
            }
          }

          # 7. README
          readme_lines <- c(
            "# DE-LIMP Run Comparison Export",
            "",
            "## How to use with Claude",
            "1. Upload all CSV files and claude_prompt.md to a Claude conversation",
            "2. Paste the contents of claude_prompt.md as your first message",
            "3. Claude will analyze the comparison and answer the questions",
            "",
            "## Files",
            "- `claude_prompt.md` - Opening prompt with context and questions",
            "- `comparison_context.md` - Tool context and Gemini analysis (if generated)",
            "- `settings_diff.csv` - Parameter comparison table",
            "- `protein_universe.csv` - All proteins with tier classification",
            "- `de_results_combined.csv` - Full DE stats from both runs (if available)",
            "- `discordant_proteins.csv` - Disagreements with per-protein diagnostics (if available)",
            "- `diann_params_run_a.txt` / `diann_params_run_b.txt` - Parsed DIA-NN log parameters (if uploaded)",
            "",
            "## Search Provenance",
            provenance_notes
          )
          readme_file <- file.path(tmp_dir, "README.md")
          writeLines(readme_lines, readme_file)
          files_to_zip <- c(files_to_zip, readme_file)

          setProgress(0.9, detail = "Creating zip...")
          old_wd <- setwd(tmp_dir)
          on.exit(setwd(old_wd), add = TRUE)
          zip(file, basename(files_to_zip))

          message(sprintf("[DE-LIMP] Comparator export: %d files", length(files_to_zip)))
        })
      }, error = function(e) {
        message("[DE-LIMP] Comparator export FAILED: ", e$message)
        showNotification(paste("Export error:", e$message), type = "error", duration = 15)
      })
    },
    contentType = "application/zip"
  )

}
