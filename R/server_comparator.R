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
  ss <- session_obj$diann_search_settings
  sp <- ss$search_params

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
    settings    = list(
      software        = "DE-LIMP",
      delimp_version  = session_obj$app_version %||% "unknown",
      limpa_version   = tryCatch(as.character(packageVersion("limpa")), error = function(e) "unknown"),
      limma_version   = tryCatch(as.character(packageVersion("limma")), error = function(e) "unknown"),
      diann_version   = session_obj$diann_version %||% "unknown",
      fdr_threshold   = session_obj$fdr_threshold %||% "0.05",
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
      search_mz_min   = as.character(sp$min_pr_mz %||% "300"),
      search_mz_max   = as.character(sp$max_pr_mz %||% "1800"),
      library         = as.character(sp$library %||% ss$library_path %||% "FASTA-predicted")
    ),
    contrasts   = contrasts
  )
}

#' Parse Spectronaut protein group report (quantities file + optional DE candidates file)
parse_spectronaut <- function(file_path, de_file_path = NULL) {
  df <- data.table::fread(file_path, data.table = TRUE)

  protein_col <- grep("ProteinGroup|ProteinAccession", names(df), value = TRUE)[1]
  gene_col    <- grep("^PG\\.Genes$|^Gene$|Genes", names(df), value = TRUE)[1]
  quant_cols  <- grep("PG\\.Quantity$", names(df), value = TRUE)
  npep_col    <- grep("NrOfStrippedSequences", names(df), value = TRUE)[1]

  # The quantities file may also contain DE columns (single-file export)
  logfc_col   <- grep("Log2Ratio|log2FC", names(df), value = TRUE)[1]
  qval_col    <- grep("Qvalue|q\\.value|adj", names(df), ignore.case = TRUE, value = TRUE)[1]
  pval_col    <- grep("^.*Pvalue$|^.*p\\.value$", names(df), ignore.case = TRUE, value = TRUE)[1]

  # Validate required columns
  if (is.na(protein_col)) stop("Spectronaut export missing protein column (PG.ProteinGroups)")
  if (length(quant_cols) == 0 && is.na(logfc_col)) stop("Spectronaut export has no quantities or DE stats")

  protein_ids <- normalize_protein_id(df[[protein_col]])

  # Build DE stats from quantities file (single-file mode) or separate candidates file
  de_stats_df <- NULL
  has_de <- !is.na(logfc_col) && !is.na(pval_col)

  if (!is.null(de_file_path)) {
    # Two-file mode: separate Candidates export
    de_df <- data.table::fread(de_file_path, data.table = TRUE)
    de_prot_col  <- grep("ProteinGroup|ProteinAccession", names(de_df), value = TRUE)[1]
    de_gene_col  <- grep("^PG\\.Genes$|^Gene$|Genes", names(de_df), value = TRUE)[1]
    de_logfc_col <- grep("Log2Ratio|log2FC|AVG\\.Log2\\.Ratio", names(de_df), value = TRUE)[1]
    de_qval_col  <- grep("Qvalue|q\\.value|adj", names(de_df), ignore.case = TRUE, value = TRUE)[1]
    de_pval_col  <- grep("^.*Pvalue$|^.*p\\.value$|PValue", names(de_df), ignore.case = TRUE, value = TRUE)[1]

    if (is.na(de_prot_col)) stop("Candidates file missing protein column")
    if (is.na(de_logfc_col)) stop("Candidates file missing Log2Ratio / logFC column")

    de_protein_ids <- normalize_protein_id(de_df[[de_prot_col]])

    # Handle multiple comparisons: detect comparison columns
    comp_col <- grep("^Comparison$|^R\\.Condition$|^Contrast$", names(de_df), value = TRUE)[1]
    if (!is.na(comp_col)) {
      comparisons <- unique(de_df[[comp_col]])
      de_stats_list <- lapply(comparisons, function(comp) {
        sub_df <- de_df[de_df[[comp_col]] == comp, ]
        sub_ids <- normalize_protein_id(sub_df[[de_prot_col]])
        data.frame(
          protein_id = sub_ids,
          gene       = if (!is.na(de_gene_col)) as.character(sub_df[[de_gene_col]]) else sub_ids,
          logFC      = as.numeric(sub_df[[de_logfc_col]]),
          P.Value    = as.numeric(if (!is.na(de_pval_col)) sub_df[[de_pval_col]] else NA),
          adj.P.Val  = as.numeric(if (!is.na(de_qval_col)) sub_df[[de_qval_col]] else NA),
          stringsAsFactors = FALSE
        )
      })
      names(de_stats_list) <- comparisons
      de_stats_df <- de_stats_list
    } else {
      de_stats_df <- list(comparison = data.frame(
        protein_id = de_protein_ids,
        gene       = if (!is.na(de_gene_col)) as.character(de_df[[de_gene_col]]) else de_protein_ids,
        logFC      = as.numeric(de_df[[de_logfc_col]]),
        P.Value    = as.numeric(if (!is.na(de_pval_col)) de_df[[de_pval_col]] else NA),
        adj.P.Val  = as.numeric(if (!is.na(de_qval_col)) de_df[[de_qval_col]] else NA),
        stringsAsFactors = FALSE
      ))
    }
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
      fdr_threshold    = "Not available in export",
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
      fdr_threshold    = "Not available in export",
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
      fdr_threshold   = "N/A",
      n_proteins_total = length(unique(protein_ids)),
      n_samples        = length(quant_cols)
    ),
    contrasts = NULL
  )
}

# --- Sample Matching ---

#' Match samples between two runs by normalized filename
match_samples <- function(names_a, names_b, source_b = "delimp") {
  strip_ext <- function(x) {
    x <- tools::file_path_sans_ext(basename(x))
    # Strip common suffixes
    x <- gsub("\\.(raw|d|mzML|wiff)$", "", x, ignore.case = TRUE)
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

  data.frame(
    run_a    = names_a,
    run_b    = ifelse(is.na(matched_idx), NA_character_, names_b[matched_idx]),
    status   = ifelse(is.na(matched_idx), "unresolved", "matched"),
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
    "fdr_threshold", "lfc_threshold", "min_peptides",
    "imputation", "fdr_method", "covariates", "contrast_string",
    # Search parameters
    "enzyme", "mass_acc_ms2", "mass_acc_ms1", "scan_window",
    "search_mz_min", "search_mz_max", "library",
    "fasta_file", "fasta_entries",
    # DIA-NN log-derived
    "pg_level", "proteoforms", "reanalyse",
    "library_source", "n_precursors_lib", "pipeline_step",
    "min_fr_mz", "max_fr_mz", "min_pep_len", "max_pep_len",
    "missed_cleavages",
    # Data stats
    "n_precursors", "n_proteins_total", "n_samples",
    "precursor_mz_range", "charge_range"
  )

  # Nice labels
  labels <- c(
    # Pipeline
    "Software", "Version", "DE-LIMP Version", "FragPipe Version",
    "DIA-NN Version", "limpa Version", "limma Version",
    # DE settings
    "DE Engine", "Normalization", "Protein Rollup",
    "FDR Threshold", "logFC Threshold", "Min Peptides",
    "Imputation", "FDR Method", "Covariates", "Contrast Formula",
    # Search parameters
    "Enzyme", "MS2 Mass Accuracy", "MS1 Mass Accuracy", "Scan Window",
    "Search Min m/z", "Search Max m/z", "Spectral Library",
    "FASTA File", "FASTA Entries",
    # DIA-NN log-derived
    "Protein Grouping (pg-level)", "--proteoforms", "--reanalyse",
    "Library Source", "Precursors in Library", "Pipeline Step",
    "Min Fragment m/z", "Max Fragment m/z", "Min Peptide Length", "Max Peptide Length",
    "Missed Cleavages",
    # Data stats
    "Precursors", "Protein Groups", "Samples",
    "Precursor m/z Range", "Charge Range"
  )

  # Section boundaries (insert header before these indices)
  section_starts <- list(
    c(1, "--- Pipeline ---"),
    c(which(params == "de_engine"), "--- DE Analysis ---"),
    c(which(params == "enzyme"), "--- Search Parameters ---"),
    c(which(params == "pg_level"), "--- DIA-NN Search (from log) ---"),
    c(which(params == "n_precursors"), "--- Data Statistics ---")
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
                                   source_b, global_offset = 0) {
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

  # Classify each protein: Up / Down / NS
  classify_de <- function(logfc, adjp, threshold = 0.05) {
    ifelse(is.na(adjp) | adjp >= threshold, "NS",
           ifelse(logfc > 0, "Up", "Down"))
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

  # Identify discordant proteins (different DE status)
  merged$concordant <- merged$status_a == merged$status_b
  discordant <- merged[!merged$concordant, ]

  # Get peptide counts and missing percentages
  if (!is.null(run_a$n_peptides)) {
    discordant <- merge(discordant,
      run_a$n_peptides[, c("protein_id", "n_peptides")],
      by = "protein_id", all.x = TRUE, suffixes = c("", "_npA"))
    names(discordant)[names(discordant) == "n_peptides"] <- "n_peptides_A"
  } else {
    discordant$n_peptides_A <- NA_integer_
  }
  if (!is.null(run_b$n_peptides)) {
    discordant <- merge(discordant,
      run_b$n_peptides[, c("protein_id", "n_peptides")],
      by = "protein_id", all.x = TRUE)
    names(discordant)[names(discordant) == "n_peptides"] <- "n_peptides_B"
  } else {
    discordant$n_peptides_B <- NA_integer_
  }

  # Add missing percentages
  discordant$missing_pct_A <- run_a$missing_pct[match(discordant$protein_id, names(run_a$missing_pct))]
  discordant$missing_pct_B <- run_b$missing_pct[match(discordant$protein_id, names(run_b$missing_pct))]

  # Assign hypotheses
  if (nrow(discordant) > 0) {
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
    n_total          = nrow(merged)
  )
}

#' Assign diagnostic hypothesis for a discordant protein (7 rules, tool-aware)
assign_hypothesis <- function(row, source_b, global_offset = 0) {
  logfc_diff   <- abs(row$logFC_A - row$logFC_B)
  same_dir     <- sign(row$logFC_A) == sign(row$logFC_B)
  p_diff_large <- abs(log10(row$adjP_A + 1e-10) - log10(row$adjP_B + 1e-10)) > 1
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
      "spectronaut"      = " (Spectronaut local regression vs DIA-NN RT-dependent normalization)",
      "fragpipe_analyst" = " (MaxLFQ rollup vs DPC-Quant - a structural difference between these protein quantification algorithms)",
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
      " Both use limma, so this likely reflects imputation differences: FragPipe-Analyst uses Perseus-style; DE-LIMP uses DPC-CN modelling."
    } else if (source_b == "spectronaut") {
      " Spectronaut uses a different variance model than limma/limpa."
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
build_gemini_comparator_prompt <- function(comp_results, mofa_obj = NULL) {
  stats    <- comp_results$summary_stats
  top_disc <- head(comp_results$de_concordance$discordant_table, 10)
  source_b <- stats$source_b

  tool_context <- switch(source_b,
    "delimp" =
      "Both runs used DE-LIMP (DIA-NN -> limpa/limma pipeline). The comparison isolates the effect of search or analysis parameter differences on the same raw data.",
    "spectronaut" =
      paste("Run A: DE-LIMP (DIA-NN search -> DPC-Quant rollup -> limma DE).",
            "Run B: Spectronaut (proprietary DIA search -> PG.Quantity rollup -> internal statistics).",
            "Key structural differences:",
            "- Normalization: DIA-NN RT-dependent vs Spectronaut local regression",
            "- Protein rollup: DPC-Quant empirical Bayes vs PG.Quantity",
            "- Statistical model: limma moderated t-test vs Spectronaut paired t-test",
            "A non-zero global intensity offset is expected."),
    "fragpipe_analyst" =
      paste("Run A: DE-LIMP (DIA-NN search -> DPC-Quant rollup -> limma DE).",
            "Run B: FragPipe + FragPipe-Analyst (MSFragger -> IonQuant MaxLFQ rollup -> limma DE).",
            "Key structural differences:",
            "- Search engine: DIA-NN (DIA-optimized) vs MSFragger (DDA/DIA)",
            "- Protein rollup: DPC-Quant empirical Bayes vs MaxLFQ consensus ratios",
            "- Normalization: DIA-NN RT-dependent vs IonQuant (optional VSN)",
            "- Missingness: DIA-NN MBR vs IonQuant MBR (different implementations)",
            "- Imputation: DPC-CN modelling vs Perseus-style (FP-Analyst default)",
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
  diff_rows <- settings_diff[settings_diff$match == "differs", ]
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

  paste0(
    "You are analyzing a proteomics run comparison in DE-LIMP.\n\n",
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
    "Please provide:\n",
    "1. A plain-language summary of the dominant cause of disagreement, accounting for the specific tools compared\n",
    "2. Whether the difference appears systematic (many proteins) or idiosyncratic (specific proteins)\n",
    "3. Which run's results you would trust more for this experiment type, and why\n",
    "4. Two specific follow-up actions that would resolve the ambiguity"
  )
}

#' Build Claude export prompt
build_claude_comparator_prompt <- function(comp_results, gemini_narrative = NULL) {
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

  paste0(
    "I am sharing a proteomics run comparison from DE-LIMP.\n\n",
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
    "- discordant_proteins.csv: All disagreements with per-protein diagnostic flags\n",
    "- de_results_combined.csv: Full DE stats from both runs side-by-side\n",
    "- settings_diff.csv: Parameter comparison table\n",
    "- protein_universe.csv: All proteins with tier classification\n",
    "- diann_search_params.txt: DIA-NN search parameters for Run A (if available)\n",
    "- precursor_summary_discordant.csv: Precursor-level data for discordant proteins (if available)\n\n",
    "QUESTIONS:\n",
    "1. Based on the discordant protein patterns and tool differences, what is the most likely root cause?\n",
    "2. Are there proteins where the two runs disagree biologically (not just statistically)?\n",
    "3. If precursor_summary_discordant.csv is included, do any discordant proteins have unusual ",
    "precursor characteristics (few precursors, extreme m/z, single charge state) that might explain the disagreement?\n",
    "4. What additional information would resolve the ambiguity?"
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
        "==========================================================",
        "",
        "Export TWO files from Spectronaut for the Run Comparator:",
        "",
        "FILE 1: Protein Quantities (required)",
        "---------------------------------------",
        "  Report > Report Schema > Protein Group Pivot Report",
        "",
        "  Required columns:",
        "    - PG.ProteinGroups      (protein accessions)",
        "    - [Sample].PG.Quantity  (per-sample quantities, auto-included in pivot)",
        "",
        "  Recommended columns:",
        "    - PG.Genes                    (gene symbols)",
        "    - PG.NrOfStrippedSequences    (peptide count per protein)",
        "",
        "  How to set up:",
        "    1. Go to Report tab in Spectronaut",
        "    2. Select 'Protein Group Pivot Report' schema (or create custom)",
        "    3. Make sure the report level is 'Protein Group'",
        "    4. Ensure columns above are checked",
        "    5. Export as .tsv (tab-separated)",
        "",
        "",
        "FILE 2: Candidates / DE Results (optional, enables DE Concordance)",
        "------------------------------------------------------------------",
        "  Report > Report Schema > Candidates Report",
        "",
        "  Required columns:",
        "    - PG.ProteinGroups          (protein accessions)",
        "    - AVG.Log2.Ratio            (log2 fold change)",
        "    - Qvalue                    (adjusted p-value / FDR)",
        "",
        "  Recommended columns:",
        "    - PG.Genes                  (gene symbols)",
        "    - PValue                    (raw p-value)",
        "    - Comparison                (comparison label, for multi-comparison exports)",
        "",
        "  How to set up:",
        "    1. Go to Report tab in Spectronaut",
        "    2. Run your statistical analysis (Post Analysis > Differential Expression)",
        "    3. Select 'Candidate' report schema",
        "    4. Ensure columns above are checked",
        "    5. Export as .tsv",
        "",
        "",
        "WHAT EACH FILE ENABLES",
        "----------------------",
        "  Quantities only:     Settings Diff, Protein Universe, Quantification",
        "  Quantities + DE:     All above + DE Concordance, Hypothesis Engine, AI Analysis",
        "",
        "",
        "TIPS",
        "----",
        "  - Column names are auto-detected, so exact schema isn't critical",
        "  - If you have a single combined file with both quantities and DE stats,",
        "    you can upload it as File 1 and skip File 2",
        "  - The comparator normalizes protein IDs across formats (sp|ACC|NAME -> ACC)",
        "  - For best results, use the same FASTA for both DIA-NN and Spectronaut",
        "",
        "Generated by DE-LIMP Run Comparator"
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

  # Re-parse Spectronaut when either file changes
  observeEvent(input$comparator_run_b_spectronaut, {
    req(input$comparator_mode == "delimp_spectronaut")
    tryCatch({
      de_path <- input$comparator_run_b_spectronaut_de$datapath
      parsed <- parse_spectronaut(input$comparator_run_b_spectronaut$datapath,
                                  de_file_path = de_path)
      comp_run_b(parsed)
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
      auto_match <- if (contrasts_a[1] %in% contrasts_b) contrasts_a[1] else contrasts_b[1]
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
        )
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
                                run_b$source)
    n_matched    <- sum(sample_map$status == "matched")
    n_unresolved <- sum(sample_map$status == "unresolved")

    if (n_unresolved > 0) {
      unmatched <- sample_map$run_a[sample_map$status == "unresolved"]
      div(class = "alert alert-warning mt-2",
        icon("triangle-exclamation"), " ",
        tags$b(paste0(n_unresolved, " sample(s) could not be matched: ")),
        paste(unmatched, collapse = ", "),
        tags$br(),
        "Check that both runs analyze the same samples."
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
                                run_b$source)
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
      if (!is.null(run_b$de_stats)) {
        de_conc <- tryCatch(
          compute_de_concordance(run_a, run_b, universe,
                                contrast_a, contrast_b %||% "comparison",
                                run_b$source, global_offset),
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
        n_a_only_de   = if (!is.null(de_conc)) {
          sum(de_conc$discordant_table$status_a != "NS" & de_conc$discordant_table$status_b == "NS")
        } else NA,
        n_b_only_de   = if (!is.null(de_conc)) {
          sum(de_conc$discordant_table$status_a == "NS" & de_conc$discordant_table$status_b != "NS")
        } else NA,
        global_offset = global_offset,
        min_cor       = if (!is.null(quant)) quant$min_cor else NA,
        max_cor       = if (!is.null(quant)) quant$max_cor else NA
      )

      results <- list(
        settings_diff    = settings_diff,
        protein_universe = universe,
        quant_comparison = quant,
        de_concordance   = de_conc,
        summary_stats    = summary_stats,
        sample_map       = sample_map,
        diann_log_a      = log_a,
        diann_log_b      = log_b
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
      }
    )
  })

  # ==========================================================================
  #  LAYER 1: Settings Diff
  # ==========================================================================

  output$comparator_settings_diff <- DT::renderDT({
    res <- comp_results()
    req(res)
    DT::datatable(res$settings_diff,
      rownames = FALSE,
      colnames = c("Parameter", "Run A", "Run B", "Status"),
      options  = list(pageLength = 25, dom = 't', ordering = FALSE),
      class    = "compact stripe"
    ) |>
      DT::formatStyle("match",
        target          = "row",
        backgroundColor = DT::styleEqual(
          c("differs", "unknown", "match"),
          c("#fff3cd", "#f8f9fa", "transparent")
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
  #  LAYER 4: DE Concordance
  # ==========================================================================

  output$comparator_layer4_content <- renderUI({
    res <- comp_results()
    req(res)

    # No DE stats in Run B -> show banner
    if (is.null(res$de_concordance)) {
      return(div(class = "alert alert-warning mt-3",
        icon("triangle-exclamation"), " ",
        tags$b("DE Concordance requires DE statistics. "),
        "The uploaded file contains intensities only. ",
        "Run ", tags$a("FragPipe-Analyst", href = "https://fragpipe-analyst.org/",
                       target = "_blank"),
        " and upload its DE results export to enable this layer."
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
          tags$h6("Hypothesis Distribution"),
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
    prompt <- build_gemini_comparator_prompt(res, mofa_obj)

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
                "missing_pct_A", "missing_pct_B", "hypothesis", "hypothesis_category", "confidence"),
              names(disc_export)
            )
            write.csv(disc_export[, export_cols], disc_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, disc_file)
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
          if (!is.null(values$comparator_gemini_narrative)) {
            context_lines <- c(context_lines,
              "## Gemini Analysis", "",
              values$comparator_gemini_narrative
            )
          }
          context_file <- file.path(tmp_dir, "comparison_context.md")
          writeLines(context_lines, context_file)
          files_to_zip <- c(files_to_zip, context_file)

          # 6. Claude prompt
          setProgress(0.8, detail = "Building prompt...")
          prompt <- build_claude_comparator_prompt(res, values$comparator_gemini_narrative)
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
