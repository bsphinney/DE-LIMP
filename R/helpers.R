# ==============================================================================
#  HELPER FUNCTIONS — General utilities
# ==============================================================================

# --- Covariate coercion for the DE design matrix ---------------------------
#
# Decide whether a metadata column should enter the design matrix as a
# numeric (continuous) covariate or a factor (categorical). The point is to
# stop users from accidentally turning per-sample identifiers (e.g. Run order
# 707, 708, 813, …) into a factor with ~N levels, which makes the design
# matrix rank-deficient and breaks limma's lmFit with "NA/NaN/Inf in 'y'".
#
# Heuristics, in order:
#   1. If every non-empty value parses cleanly as a finite number AND there
#      are at least, say, 5 distinct values, treat it as numeric.
#   2. Otherwise treat it as a factor.
#
# Returns: list(values = <numeric or factor vector>, kind = "numeric"|"factor",
#                n_levels = <int>, has_singletons = <logical>, singleton_levels = <chr>)
coerce_covariate_column <- function(x) {
  raw <- as.character(x)
  raw[is.na(raw)] <- ""
  nonempty <- raw[nzchar(raw)]
  if (length(nonempty) == 0) {
    return(list(values = factor(raw), kind = "factor", n_levels = 0,
                has_singletons = FALSE, singleton_levels = character(0)))
  }
  # Numeric heuristic: every non-empty value parses; ≥ 5 distinct numeric values
  numeric_try <- suppressWarnings(as.numeric(nonempty))
  all_numeric <- all(is.finite(numeric_try))
  enough_unique_numeric <- length(unique(numeric_try)) >= 5
  if (all_numeric && enough_unique_numeric) {
    out <- suppressWarnings(as.numeric(raw))
    # leave NAs in place; lmFit + limma handle row-level NAs OK as long as
    # the design column itself isn't all NA
    return(list(values = out, kind = "numeric",
                n_levels = length(unique(numeric_try)),
                has_singletons = FALSE, singleton_levels = character(0)))
  }
  # Factor path
  fac <- factor(raw)
  tab <- table(fac[nzchar(as.character(fac))])
  singletons <- names(tab)[tab == 1]
  list(values = fac, kind = "factor", n_levels = length(levels(fac)),
       has_singletons = length(singletons) > 0,
       singleton_levels = singletons)
}

# --- Diagnose a rank-deficient design matrix --------------------------------
#
# Run before lmFit / dpcDE. Returns NULL if the design is full-rank;
# otherwise returns a single string naming the columns that are not
# estimable so the caller can build a helpful user-facing error.
diagnose_design_rank <- function(design) {
  qr_d <- qr(design)
  if (qr_d$rank == ncol(design)) return(NULL)
  estimable <- qr_d$pivot[seq_len(qr_d$rank)]
  not_estimable <- setdiff(seq_len(ncol(design)), estimable)
  bad_cols <- colnames(design)[not_estimable]
  if (length(bad_cols) == 0) bad_cols <- as.character(not_estimable)
  preview <- if (length(bad_cols) > 6) {
    paste0(paste(head(bad_cols, 6), collapse = ", "),
           " (+", length(bad_cols) - 6, " more)")
  } else {
    paste(bad_cols, collapse = ", ")
  }
  sprintf("Design matrix is rank-deficient (rank %d of %d). %d coefficient(s) are not estimable: %s",
          qr_d$rank, ncol(design), length(not_estimable), preview)
}

# --- QuantUMS Score Pre-Filter ----------------------------------------------
#
# Filter a DIA-NN report.parquet by QuantUMS quality scores BEFORE handing
# the file to limpa::readDIANN(). limpa drops the QuantUMS columns during
# read, so any filtering on Empirical.Quality / PG.MaxLFQ.Quality must
# happen at the parquet stage.
#
# Reference: Moschem et al., J. Proteome Res. 2025, 24, 3860–3873
# (DOI: 10.1021/acs.jproteome.5c00009). Recommended cutoffs ≥ 0.75 for
# eQ (Empirical.Quality) and pgQ (PG.MaxLFQ.Quality); the qQ
# (Quantity.Quality) score is intentionally not exposed because the paper
# shows it has negligible impact.
#
# Behaviour:
#  - If both cutoffs are <= 0, returns `parquet_path` unchanged (no work).
#  - Otherwise reads the parquet via arrow, drops rows where the named
#    column is below the cutoff, writes the survivors to a temp parquet,
#    and returns the temp path.
#  - If a column is missing (older DIA-NN that predates QuantUMS), the
#    corresponding cutoff is silently skipped and a message is emitted.
#
# Returns: list(path = <parquet path to use>, n_in = <input rows>,
#               n_out = <surviving rows>, applied = <character vec of
#               filters that ran>)
filter_quantums_parquet <- function(parquet_path, eq_cutoff = 0, pgq_cutoff = 0) {
  if ((is.null(eq_cutoff)  || is.na(eq_cutoff)  || eq_cutoff  <= 0) &&
      (is.null(pgq_cutoff) || is.na(pgq_cutoff) || pgq_cutoff <= 0)) {
    return(list(path = parquet_path, n_in = NA_integer_, n_out = NA_integer_,
                applied = character(0)))
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    message("[QuantUMS filter] arrow package missing — skipping filter.")
    return(list(path = parquet_path, n_in = NA_integer_, n_out = NA_integer_,
                applied = character(0)))
  }

  ds <- arrow::open_dataset(parquet_path, format = "parquet")
  cols <- names(ds$schema)
  applied <- character(0)
  flt <- ds

  if (!is.null(eq_cutoff) && !is.na(eq_cutoff) && eq_cutoff > 0) {
    if ("Empirical.Quality" %in% cols) {
      flt <- dplyr::filter(flt, Empirical.Quality >= !!eq_cutoff)
      applied <- c(applied, sprintf("Empirical.Quality >= %.2f", eq_cutoff))
    } else {
      message("[QuantUMS filter] Empirical.Quality column absent — eQ filter skipped (DIA-NN < 1.8.2 β39?)")
    }
  }
  if (!is.null(pgq_cutoff) && !is.na(pgq_cutoff) && pgq_cutoff > 0) {
    if ("PG.MaxLFQ.Quality" %in% cols) {
      flt <- dplyr::filter(flt, PG.MaxLFQ.Quality >= !!pgq_cutoff)
      applied <- c(applied, sprintf("PG.MaxLFQ.Quality >= %.2f", pgq_cutoff))
    } else {
      message("[QuantUMS filter] PG.MaxLFQ.Quality column absent — pgQ filter skipped.")
    }
  }

  if (length(applied) == 0) {
    return(list(path = parquet_path, n_in = NA_integer_, n_out = NA_integer_,
                applied = character(0)))
  }

  n_in <- tryCatch(as.integer(ds %>% dplyr::summarise(n = dplyr::n()) %>% dplyr::collect() %>% .$n),
                   error = function(e) NA_integer_)

  out_path <- tempfile(pattern = "quantums_filtered_", fileext = ".parquet")
  arrow::write_parquet(dplyr::collect(flt), out_path)

  n_out <- tryCatch(as.integer(arrow::open_dataset(out_path, format = "parquet") %>%
                               dplyr::summarise(n = dplyr::n()) %>%
                               dplyr::collect() %>% .$n),
                    error = function(e) NA_integer_)

  message(sprintf("[QuantUMS filter] %s — kept %s / %s precursors (%.1f%%)",
                  paste(applied, collapse = " AND "),
                  format(n_out, big.mark = ","),
                  format(n_in,  big.mark = ","),
                  100 * n_out / max(n_in, 1)))

  list(path = out_path, n_in = n_in, n_out = n_out, applied = applied)
}

# --- QC Stats Calculation ---
# Memory-optimized: reads only needed columns via Arrow col_select,
# then aggregates before collecting into R memory.
get_diann_stats_r <- function(file_path) {
  tryCatch({
    # Check which columns are available without reading data
    available_cols <- names(arrow::read_parquet(file_path, as_data_frame = FALSE))

    needed_cols <- c("Run", "Protein.Group", "Q.Value")
    has_pg_q <- "PG.Q.Value" %in% available_cols
    has_ms1  <- "Ms1.Apex.Area" %in% available_cols
    if (has_pg_q) needed_cols <- c(needed_cols, "PG.Q.Value")
    if (has_ms1)  needed_cols <- c(needed_cols, "Ms1.Apex.Area")

    # Read only the needed columns (saves 70-80% memory for large files)
    df <- arrow::read_parquet(file_path, col_select = dplyr::all_of(needed_cols))
    if ("Q.Value" %in% names(df)) df <- df %>% dplyr::filter(Q.Value <= 0.01)

    stats_df <- df %>%
      dplyr::group_by(Run) %>%
      dplyr::summarise(
        Precursors = dplyr::n(),
        Proteins = if(has_pg_q) {
          dplyr::n_distinct(Protein.Group[PG.Q.Value <= 0.01])
        } else {
          dplyr::n_distinct(Protein.Group)
        },
        MS1_Signal = if(has_ms1) sum(Ms1.Apex.Area, na.rm = TRUE) else NA_real_,
        .groups = 'drop'
      ) %>% dplyr::arrange(Run)

    # Free the large intermediate df immediately
    rm(df); gc(verbose = FALSE)

    return(stats_df)
  }, error = function(e) { data.frame(Run = "Error", Precursors = 0, Proteins = 0, MS1_Signal = 0) })
}

# --- Z-Score Utility ---
cal_z_score <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }

# --- DPC-Quant Detection Class ---
# Classifies each protein based on n.observations across all samples.
# Returns a character vector: "Detected_All", "Detected_Partial", or "Inferred_All".
# n_obs_mat: matrix with proteins as rows, samples as columns (from y_protein$other$n.observations)
# protein_ids: character vector of protein IDs to classify (must match rownames of n_obs_mat)
compute_detection_class <- function(n_obs_mat, protein_ids) {
  if (is.null(n_obs_mat)) return(rep(NA_character_, length(protein_ids)))
  rn <- rownames(n_obs_mat)
  vapply(protein_ids, function(pid) {
    idx <- match(pid, rn)
    if (is.na(idx)) return(NA_character_)
    obs <- n_obs_mat[idx, ]
    if (all(obs > 0)) {
      "Detected_All"
    } else if (all(obs == 0)) {
      "Inferred_All"
    } else {
      "Detected_Partial"
    }
  }, character(1))
}

# --- Auto-detect Organism ---
detect_organism_db <- function(protein_ids) {
  ORGANISM_DB_MAP <- list(
    "_HUMAN" = "org.Hs.eg.db", "_MOUSE" = "org.Mm.eg.db", "_RAT"   = "org.Rn.eg.db",
    "_BOVIN" = "org.Bt.eg.db", "_CANLF" = "org.Cf.eg.db", "_CHICK" = "org.Gg.eg.db",
    "_DROME" = "org.Dm.eg.db", "_CAEEL" = "org.Ce.eg.db", "_DANRE" = "org.Dr.eg.db",
    "_YEAST" = "org.Sc.sgd.db", "_ARATH" = "org.At.tair.db", "_PIG"   = "org.Ss.eg.db"
  )
  for (suffix in names(ORGANISM_DB_MAP)) {
    if (any(grepl(suffix, protein_ids, ignore.case = TRUE))) {
      return(ORGANISM_DB_MAP[[suffix]])
    }
  }
  return("org.Hs.eg.db")
}
