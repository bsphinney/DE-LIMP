#!/usr/bin/env Rscript
# build_maxlfq.R  --  Port of DE-LIMP's build_maxlfq_pipeline() (R/helpers.R).
# Reads a DIA-NN report, applies ID-FDR + optional QuantUMS filters, pivots
# PG.MaxLFQ to a protein x run matrix, log2-transforms, and quantile-normalizes
# with limma::normalizeBetweenArrays (the DE-LIMP default for the MaxLFQ path).
#
# Returns: list(E, genes, descriptor, n_obs, filters_applied).

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

build_maxlfq <- function(report_path, format = "parquet", q_cutoff = 0.01,
                         eq_cutoff = 0, pgq_cutoff = 0, keep_runs = NULL) {
  stopifnot(requireNamespace("dplyr", quietly = TRUE),
            requireNamespace("tidyr", quietly = TRUE))

  if (identical(format, "parquet")) {
    if (!requireNamespace("arrow", quietly = TRUE)) stop("arrow required for parquet input.")
    ds   <- arrow::open_dataset(report_path, format = "parquet")
    cols <- names(ds$schema)
  } else {
    ds   <- arrow::read_delim_arrow(report_path, delim = "\t")  # arrow handles tsv too
    cols <- names(ds)
  }

  needed   <- c("Run", "Protein.Group", "PG.MaxLFQ", "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value")
  # PG.Q.Value / Global.* are optional because older reports lack them, but they are
  # what actually controls FDR at the protein and experiment level — see below.
  optional <- c("Empirical.Quality", "PG.MaxLFQ.Quality", "Genes", "Protein.Names",
                "PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value")
  miss <- setdiff(needed, cols)
  if (length(miss)) stop("MaxLFQ: missing required columns: ", paste(miss, collapse = ", "))

  sel <- c(needed, intersect(optional, cols))
  flt <- if (identical(format, "parquet")) dplyr::select(ds, dplyr::all_of(sel)) else ds[, sel]
  filters_applied <- character(0)

  # Run-level + library q-values do NOT control FDR across an experiment: a union of
  # run-level-passing IDs over many runs sits well above the nominal cutoff, and the
  # inflated protein list also inflates the family size m that BH corrects over.
  # Add DIA-NN's protein-level and experiment-wide q-values when present.
  if (!is.na(q_cutoff) && q_cutoff > 0) {
    flt <- dplyr::filter(flt, Q.Value <= !!q_cutoff,
                              Lib.Q.Value <= !!q_cutoff,
                              Lib.PG.Q.Value <= !!q_cutoff)
    .fdr <- c("Q", "Lib.Q", "Lib.PG.Q")
    for (.qc in c("PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value")) {
      if (.qc %in% cols) {
        flt <- dplyr::filter(flt, .data[[.qc]] <= !!q_cutoff)
        .fdr <- c(.fdr, .qc)
      }
    }
    filters_applied <- c(filters_applied,
                         sprintf("%s <= %.3f", paste(.fdr, collapse = "/"), q_cutoff))
  }
  if (!is.na(eq_cutoff) && eq_cutoff > 0 && "Empirical.Quality" %in% cols) {
    flt <- dplyr::filter(flt, Empirical.Quality >= !!eq_cutoff)
    filters_applied <- c(filters_applied, sprintf("Empirical.Quality >= %.2f", eq_cutoff))
  }
  if (!is.na(pgq_cutoff) && pgq_cutoff > 0 && "PG.MaxLFQ.Quality" %in% cols) {
    flt <- dplyr::filter(flt, PG.MaxLFQ.Quality >= !!pgq_cutoff)
    filters_applied <- c(filters_applied, sprintf("PG.MaxLFQ.Quality >= %.2f", pgq_cutoff))
  }
  if (!is.null(keep_runs) && length(keep_runs))
    flt <- dplyr::filter(flt, Run %in% !!keep_runs)

  rows <- if (identical(format, "parquet")) dplyr::collect(flt) else as.data.frame(flt)
  if (!nrow(rows)) stop("MaxLFQ: no rows survived the filters. Loosen QuantUMS cutoffs.")

  # one PG.MaxLFQ per (Protein.Group, Run): DIA-NN broadcasts it across precursor rows.
  # Assert the broadcast rather than trust it — if it fails, max() silently biases the
  # affected cells upward. It is also the reason the QuantUMS cutoffs above only censor
  # cells and never improve a retained value: see docs/QUANTUMS_MAXLFQ_NOTES.md.
  pg_run <- rows |>
    dplyr::group_by(Protein.Group, Run) |>
    dplyr::summarise(.nd = dplyr::n_distinct(PG.MaxLFQ),
                     PG.MaxLFQ = max(PG.MaxLFQ, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(PG.MaxLFQ = ifelse(is.finite(PG.MaxLFQ), PG.MaxLFQ, NA_real_))
  .nbad <- sum(pg_run$.nd > 1, na.rm = TRUE)
  if (.nbad > 0)
    warning(sprintf("MaxLFQ: PG.MaxLFQ not broadcast in %d (Protein.Group, Run) cell(s); max() may bias upward.", .nbad))
  pg_run$.nd <- NULL

  wide <- tidyr::pivot_wider(pg_run, id_cols = Protein.Group,
                             names_from = Run, values_from = PG.MaxLFQ)
  prot_ids <- wide$Protein.Group
  E <- as.matrix(wide[, -1, drop = FALSE]); rownames(E) <- prot_ids
  E[E <= 0 | !is.finite(E)] <- NA_real_
  E_pre <- log2(E)

  # limma::normalizeQuantiles() calls approx() per column and dies with
  # "need at least two non-NA values to interpolate" if a run has <2 quantified
  # proteins. An aggressive QuantUMS cutoff does exactly that: at eq/pgq >= 0.75 one
  # run lost every row (so it never became a column and the job survived), but at
  # >= 0.50 a run survived with a single value and the analysis crashed.
  .nf <- colSums(is.finite(E_pre))
  .degen <- colnames(E_pre)[.nf < 2]
  if (length(.degen)) {
    warning(sprintf("MaxLFQ: dropping %d run(s) with <2 quantified proteins after filtering: %s",
                    length(.degen), paste(.degen, collapse = ", ")))
    E_pre <- E_pre[, .nf >= 2, drop = FALSE]
    filters_applied <- c(filters_applied,
                         sprintf("dropped %d degenerate run(s) (<2 proteins)", length(.degen)))
  }
  if (ncol(E_pre) < 2) stop("MaxLFQ: fewer than 2 usable runs survived the filters.")

  if (requireNamespace("limma", quietly = TRUE)) {
    E <- limma::normalizeBetweenArrays(E_pre, method = "quantile")
  } else {
    cm <- apply(E_pre, 2, stats::median, na.rm = TRUE)
    E  <- sweep(E_pre, 2, cm - stats::median(cm, na.rm = TRUE))
  }
  rownames(E) <- prot_ids
  n_obs <- ifelse(is.na(E), 0L, 1L)

  ann_cols <- intersect(c("Genes", "Protein.Names"), names(rows))
  ann <- if (length(ann_cols)) {
    rows |>
      dplyr::group_by(Protein.Group) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(ann_cols),
        ~ names(sort(table(.x), decreasing = TRUE))[1] %||% NA_character_), .groups = "drop")
  } else data.frame(Protein.Group = unique(rows$Protein.Group))
  genes <- merge(data.frame(Protein.Group = prot_ids), ann, by = "Protein.Group",
                 all.x = TRUE, sort = FALSE)
  rownames(genes) <- genes$Protein.Group

  list(
    E = E, genes = genes, n_obs = n_obs, filters_applied = filters_applied,
    descriptor = list(
      pipeline_id    = "maxlfq",
      display_label  = "MaxLFQ + limma",
      rollup_method  = "DIA-NN PG.MaxLFQ",
      de_engine      = "limma::lmFit -> contrasts.fit -> eBayes (NA-tolerant per row)",
      missing_policy = "NAs left in place; limma drops them per row. All-missing-in-one-condition proteins are on/off calls.",
      citation       = "Quantification: DIA-NN MaxLFQ (Demichev et al. 2020, Nat Methods 17:41). DE: limma (Ritchie et al. 2015, NAR 43:e47)."
    )
  )
}
