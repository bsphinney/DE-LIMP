# ==============================================================================
#  HELPER FUNCTIONS — General utilities
# ==============================================================================

# --- Pipeline descriptors — single source of truth for "what pipeline ran"
#
# Built to retire the `values$pipeline_mode_used` reactiveVal and the
# scattered `if (isTRUE(values$pipeline_mode_used == "maxlfq"))` guards that
# proliferated through 7 files in v3.9.x. Each pipeline puts a descriptor on
# its returned y_protein object; downstream code reads from there. Adding a
# third pipeline (e.g. DDA) requires zero edits to downstream files.
#
# Schema (every pipeline must populate all fields):
#   pipeline_id        — short slug, e.g. "dpc", "maxlfq"
#   display_label      — human label, e.g. "MaxLFQ + limma (Moschem 2025)"
#   rollup_method      — short string for tables, e.g. "DIA-NN PG.MaxLFQ"
#   normalization      — short string, e.g. "quantile (limma::normalizeBetweenArrays)"
#   de_engine          — short string, e.g. "limma::lmFit + eBayes"
#   missing_value_policy — "imputed via DPC detection model" / "kept as NA, dropped per row by limma"
#   citation           — for methods text, e.g. "Moschem et al., J. Proteome Res. 2025; 24:3860"

# Default descriptor for the historical limpa DPC-Quant pipeline. Attached to
# y_protein after dpcQuant() runs (server_data.R does this).
dpc_pipeline_descriptor <- function() {
  list(
    pipeline_id        = "dpc",
    display_label      = "DPC-Quant + limma (limpa)",
    rollup_method      = "DPC-Quant (Detection Probability Curve Quantification)",
    normalization      = "DPC-CN (Data Point Correspondence cyclic loess)",
    de_engine          = "limpa::dpcDE → contrasts.fit → eBayes",
    missing_value_policy = "Missing precursors contribute via the detection probability model — not imputed, not filled in",
    citation           = "Law CW, Smyth GK (limpa, Bioconductor)"
  )
}

# Descriptor for the MaxLFQ + limma pipeline. Attached inside
# build_maxlfq_pipeline() at quantification time.
maxlfq_pipeline_descriptor <- function() {
  list(
    pipeline_id        = "maxlfq",
    display_label      = "MaxLFQ + limma (Moschem 2025)",
    rollup_method      = "DIA-NN PG.MaxLFQ",
    normalization      = "quantile normalization (limma::normalizeBetweenArrays)",
    de_engine          = "limma::lmFit → contrasts.fit → eBayes (NA-tolerant per-row)",
    missing_value_policy = "NAs left in place; limma drops them per row at fit time. Proteins entirely missing in one condition surface in the On/Off Proteins panel.",
    citation           = "Moschem et al., J. Proteome Res. 2025; 24:3860 (DOI: 10.1021/acs.jproteome.5c00009)"
  )
}

# Read the descriptor off a y_protein. Falls back to the DPC descriptor for
# legacy y_protein objects (saved in older session.rds files) that didn't
# carry one. NULL when y_protein itself is NULL.
pipeline_descriptor <- function(y_protein) {
  if (is.null(y_protein)) return(NULL)
  d <- y_protein$other$descriptor
  if (!is.null(d)) return(d)
  # Legacy: read $other$pipeline if present (set by build_maxlfq_pipeline pre-v3.9.17)
  legacy_id <- y_protein$other$pipeline %||% NA_character_
  if (!is.na(legacy_id) && legacy_id == "maxlfq") return(maxlfq_pipeline_descriptor())
  # Default: legacy session predates descriptors → assume DPC-Quant
  dpc_pipeline_descriptor()
}

# Convenience wrappers — the canonical accessors used by every downstream file.
# Replace `if (isTRUE(values$pipeline_mode_used == "maxlfq"))` with `is_maxlfq(values$y_protein)`.
is_maxlfq <- function(y_protein) {
  d <- pipeline_descriptor(y_protein)
  !is.null(d) && identical(d$pipeline_id, "maxlfq")
}

pipeline_label <- function(y_protein) {
  d <- pipeline_descriptor(y_protein)
  if (is.null(d)) "(no pipeline ran)" else d$display_label
}

# --- safe_section: wraps an export sub-step so silent failures become visible
# Replaces the dozens of `tryCatch(error = function(e) NULL)` blocks scattered
# through claude_export_content / Complete export. Each safe_section call
# either records a successful section or appends a "SKIPPED" line to a
# manifest the export bundles into the ZIP root, so reviewers can see what's
# missing and why instead of being silently shipped a 9-file ZIP that should
# have had 12 files.
#
# Usage:
#   manifest <- new.env(parent = emptyenv()); manifest$lines <- character(0)
#   safe_section(manifest, "QC stats CSV", {
#     write.csv(qc_df, file = qc_file, row.names = FALSE)
#     files_to_zip <- c(files_to_zip, qc_file)   # caller maintains this
#   })
#   # ...then later, write the manifest:
#   writeLines(c("DE-LIMP export manifest", manifest$lines),
#              file.path(tmp_dir, "MANIFEST.txt"))
#
# Returns: TRUE on success, FALSE on caught failure (after appending to manifest).
safe_section <- function(manifest, name, expr) {
  start <- Sys.time()
  ok <- tryCatch({
    force(expr)
    elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    manifest$lines <- c(manifest$lines,
      sprintf("[OK]      %-50s (%.1fs)", name, elapsed))
    TRUE
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (nchar(msg) > 200) msg <- paste0(substr(msg, 1, 197), "...")
    manifest$lines <- c(manifest$lines,
      sprintf("[SKIPPED] %-50s -- %s", name, msg))
    message(sprintf("[DE-LIMP] safe_section '%s' failed: %s", name, msg))
    FALSE
  })
  invisible(ok)
}

# --- req_nzchar: like req() but treats "" as missing
# Standard Shiny req() considers "" truthy for character inputs, which means
# downstream code happily runs with empty strings and produces nonsense
# (limma::topTable(fit, coef = "") errors deep in the call stack). Use this
# whenever an input must be a non-empty character.
req_nzchar <- function(...) {
  args <- list(...)
  for (a in args) shiny::req(!is.null(a) && length(a) > 0 && nzchar(as.character(a)))
}

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

# --- MaxLFQ + limma pipeline (paper-faithful Moschem 2025) -----------------
#
# Read a DIA-NN report.parquet, apply the same Q-value + QuantUMS filters
# the user specified, and produce a Protein.Group x Run matrix from
# DIA-NN's already-computed PG.MaxLFQ values. NAs are preserved (not
# imputed) so limma's per-row NA handling can do its thing — that matches
# the paper's pipeline.
#
# Returns an EList-shaped list compatible with downstream DE-LIMP code:
#   $E       — log2(MaxLFQ) matrix, proteins x samples, NAs preserved
#   $genes   — data.frame(Protein.Group, Genes, Protein.Names)
#   $targets — data.frame(File.Name) — sample sheet placeholder
#   $other$n.observations — 1 / 0 matrix marking detection in each cell
#   $other$pipeline — "maxlfq"  (used by downstream code to branch)
#   $other$filters_applied — character vector of filters used (for Methods)
build_maxlfq_pipeline <- function(parquet_path, q_cutoff = 0.01,
                                   eq_cutoff = 0, pgq_cutoff = 0,
                                   keep_runs = NULL) {
  if (!requireNamespace("arrow", quietly = TRUE))
    stop("arrow package required for the MaxLFQ pipeline.")

  ds <- arrow::open_dataset(parquet_path, format = "parquet")
  cols <- names(ds$schema)

  needed <- c("Run", "Protein.Group", "PG.MaxLFQ",
              "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value")
  optional <- c("Empirical.Quality", "PG.MaxLFQ.Quality",
                "Genes", "Protein.Names")
  missing_needed <- setdiff(needed, cols)
  if (length(missing_needed) > 0)
    stop("MaxLFQ pipeline: missing required columns in parquet: ",
         paste(missing_needed, collapse = ", "))

  select_cols <- c(needed, intersect(optional, cols))
  flt <- ds %>% dplyr::select(dplyr::all_of(select_cols))

  filters_applied <- character(0)
  filter_counts <- list()

  count_rows <- function(plan) {
    tryCatch(as.integer(plan %>% dplyr::summarise(n = dplyr::n()) %>%
                          dplyr::collect() %>% .$n),
             error = function(e) NA_integer_)
  }
  filter_counts$input <- count_rows(flt)

  # Identification FDR (matches paper's Methods, page 3861)
  if (!is.null(q_cutoff) && !is.na(q_cutoff) && q_cutoff > 0) {
    flt <- flt %>%
      dplyr::filter(Q.Value <= !!q_cutoff,
                    Lib.Q.Value <= !!q_cutoff,
                    Lib.PG.Q.Value <= !!q_cutoff)
    filters_applied <- c(filters_applied,
      sprintf("Q.Value, Lib.Q.Value, Lib.PG.Q.Value <= %.3f", q_cutoff))
    filter_counts$after_fdr <- count_rows(flt)
  }
  # QuantUMS — eQ
  if (!is.null(eq_cutoff) && !is.na(eq_cutoff) && eq_cutoff > 0 &&
      "Empirical.Quality" %in% cols) {
    flt <- flt %>% dplyr::filter(Empirical.Quality >= !!eq_cutoff)
    filters_applied <- c(filters_applied,
      sprintf("Empirical.Quality >= %.2f", eq_cutoff))
    filter_counts$after_eq <- count_rows(flt)
  }
  # QuantUMS — pgQ
  if (!is.null(pgq_cutoff) && !is.na(pgq_cutoff) && pgq_cutoff > 0 &&
      "PG.MaxLFQ.Quality" %in% cols) {
    flt <- flt %>% dplyr::filter(PG.MaxLFQ.Quality >= !!pgq_cutoff)
    filters_applied <- c(filters_applied,
      sprintf("PG.MaxLFQ.Quality >= %.2f", pgq_cutoff))
    filter_counts$after_pgq <- count_rows(flt)
  }

  # Restrict to the user's kept runs (excluded_files honoured) before pivot
  if (!is.null(keep_runs) && length(keep_runs) > 0) {
    flt <- flt %>% dplyr::filter(Run %in% !!keep_runs)
    filter_counts$after_excluded_files <- count_rows(flt)
  }

  rows <- flt %>% dplyr::collect()
  if (nrow(rows) == 0)
    stop("MaxLFQ pipeline: no precursor rows survived the filters. ",
         "Loosen the QuantUMS cutoffs and try again.")

  # One PG.MaxLFQ value per (Protein.Group, Run). Take max in case multiple
  # precursor rows duplicate the protein-group MaxLFQ value (DIA-NN does
  # broadcast it across rows of a PG within a run).
  pg_run <- rows %>%
    dplyr::group_by(Protein.Group, Run) %>%
    dplyr::summarise(PG.MaxLFQ = max(PG.MaxLFQ, na.rm = TRUE),
                     .groups = "drop") %>%
    dplyr::mutate(PG.MaxLFQ = ifelse(is.finite(PG.MaxLFQ), PG.MaxLFQ, NA_real_))

  # Pivot wide: rows = proteins, cols = runs
  wide <- pg_run %>%
    tidyr::pivot_wider(id_cols = Protein.Group,
                       names_from = Run, values_from = PG.MaxLFQ)
  prot_ids <- wide$Protein.Group
  E <- as.matrix(wide[, -1, drop = FALSE])
  rownames(E) <- prot_ids

  # log2 transform (NaN/Inf -> NA)
  E[E <= 0 | !is.finite(E)] <- NA_real_
  E_pre <- log2(E)

  # Quantile-normalize across samples — standard practice for DIA matrices
  # before limma. Median-centering alone leaves between-sample variance
  # differences leaking into eBayes and crushes statistical power.
  # Use limma::normalizeBetweenArrays(method = "quantile"), the same
  # default used by FragPipe-Analyst, Spectronaut/MSstats, and DIA-NN's
  # own analyzer for cross-sample alignment.
  if (requireNamespace("limma", quietly = TRUE)) {
    E <- limma::normalizeBetweenArrays(E_pre, method = "quantile")
  } else {
    # Fallback: median-center
    col_med <- apply(E_pre, 2, function(x) stats::median(x, na.rm = TRUE))
    global_med <- stats::median(col_med, na.rm = TRUE)
    E <- sweep(E_pre, 2, col_med - global_med, FUN = "-")
  }
  rownames(E) <- prot_ids

  # Detection (1 if not NA, 0 otherwise) — used by downstream nObs logic
  n_obs <- ifelse(is.na(E), 0L, 1L)

  # Best-effort gene annotation per Protein.Group: take the most common Genes /
  # Protein.Names string per PG (in case multiple precursor rows disagree).
  ann_cols <- intersect(c("Genes", "Protein.Names"), names(rows))
  ann <- if (length(ann_cols) > 0) {
    rows %>%
      dplyr::group_by(Protein.Group) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(ann_cols),
                                     ~ names(sort(table(.x), decreasing = TRUE))[1] %||% NA_character_),
                       .groups = "drop")
  } else {
    data.frame(Protein.Group = unique(rows$Protein.Group), stringsAsFactors = FALSE)
  }
  genes <- merge(data.frame(Protein.Group = prot_ids, stringsAsFactors = FALSE),
                 ann, by = "Protein.Group", all.x = TRUE, sort = FALSE)
  rownames(genes) <- genes$Protein.Group

  list(
    E = E,
    genes = genes,
    targets = data.frame(File.Name = colnames(E), stringsAsFactors = FALSE),
    other = list(
      n.observations = n_obs,
      pipeline = "maxlfq",                       # legacy slug — retained for compat
      descriptor = maxlfq_pipeline_descriptor(), # canonical (v3.9.17+)
      filters_applied = filters_applied,
      n_proteins_in_matrix = nrow(E),
      n_runs = ncol(E),
      n_cells_total = length(E),
      n_cells_missing = sum(is.na(E)),
      # Per-filter precursor row counts so the UI can show "kept X / Y at eQ ≥ 0.75"
      filter_counts = filter_counts,
      # Pre-normalization log2 matrix retained for Norm QC visualization
      E_log2_raw = E_pre
    )
  )
}

# --- Compute "On/Off" proteins per contrast --------------------------------
#
# Surface proteins detected in ≥ n_min samples of one condition AND in 0
# samples of the other — these get NA logFC from limma so they're invisible
# in the volcano. Returns a data.frame: Protein.Group, Gene, Contrast,
# Direction (one of "Group1_only" / "Group2_only"), n_in_group1, n_in_group2,
# total_in_group1, total_in_group2.
compute_onoff_proteins <- function(E, group_factor, contrasts_list = NULL,
                                    n_min = 2, gene_lookup = NULL) {
  # contrasts_list can be a list of c(g1, g2) character pairs OR a 2-row matrix
  # where each column is a contrast. Normalize to list-of-pairs at function entry.
  if (is.matrix(contrasts_list)) {
    contrasts_list <- lapply(seq_len(ncol(contrasts_list)),
                              function(i) as.character(contrasts_list[, i]))
  }
  stopifnot(ncol(E) == length(group_factor))
  groups <- as.character(group_factor)
  unique_groups <- unique(groups[!is.na(groups) & nzchar(groups)])

  # If no contrasts given, generate all pairs
  if (is.null(contrasts_list) || length(contrasts_list) == 0) {
    if (length(unique_groups) < 2) return(NULL)
    pairs <- utils::combn(unique_groups, 2, simplify = FALSE)
    contrasts_list <- lapply(pairs, function(p) c(p[2], p[1]))
  }

  # Detection per cell
  detected <- !is.na(E)

  out <- list()
  for (con in contrasts_list) {
    g1 <- con[1]; g2 <- con[2]
    cols_g1 <- which(groups == g1)
    cols_g2 <- which(groups == g2)
    if (length(cols_g1) == 0 || length(cols_g2) == 0) next

    n1 <- rowSums(detected[, cols_g1, drop = FALSE])
    n2 <- rowSums(detected[, cols_g2, drop = FALSE])

    g1_only <- (n1 >= n_min) & (n2 == 0)
    g2_only <- (n2 >= n_min) & (n1 == 0)

    if (sum(g1_only) > 0) {
      df <- data.frame(
        Protein.Group = rownames(E)[g1_only],
        Contrast = paste0(g1, " - ", g2),
        Direction = paste0(g1, "_only"),
        n_in_group1 = n1[g1_only],
        n_in_group2 = n2[g1_only],
        total_in_group1 = length(cols_g1),
        total_in_group2 = length(cols_g2),
        stringsAsFactors = FALSE
      )
      out[[length(out) + 1]] <- df
    }
    if (sum(g2_only) > 0) {
      df <- data.frame(
        Protein.Group = rownames(E)[g2_only],
        Contrast = paste0(g1, " - ", g2),
        Direction = paste0(g2, "_only"),
        n_in_group1 = n1[g2_only],
        n_in_group2 = n2[g2_only],
        total_in_group1 = length(cols_g1),
        total_in_group2 = length(cols_g2),
        stringsAsFactors = FALSE
      )
      out[[length(out) + 1]] <- df
    }
  }
  if (length(out) == 0) return(NULL)
  result <- do.call(rbind, out)
  if (!is.null(gene_lookup)) {
    result$Gene <- gene_lookup[result$Protein.Group]
  }
  rownames(result) <- NULL
  result
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
