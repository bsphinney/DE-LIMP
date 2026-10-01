#!/usr/bin/env Rscript
# build_maxlfq.R  --  Port of DE-LIMP's build_maxlfq_pipeline() (R/helpers.R).
# Reads a DIA-NN report, applies ID-FDR + optional QuantUMS filters, pivots
# PG.MaxLFQ to a protein x run matrix, log2-transforms, and quantile-normalizes
# with limma::normalizeBetweenArrays (the DE-LIMP default for the MaxLFQ path).
# An adapted report that declares its own quantity (parquet metadata, declared_quantity()) is
# described by what it declares -- a Sage report's rows are peptides, not DIA-NN's PG.MaxLFQ.
#
# Returns: list(E, genes, descriptor, n_obs, filters_applied, contaminants).

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || (length(a) == 1 && is.na(a))) b else a

# ONE definition of the q-value column set (mirrored in diann_q_columns.py).
# run_de.R sources this file, so the constants are usually already present;
# source them here too so build_maxlfq.R also works when loaded on its own.
if (!exists("DIANN_FDR_REQUIRED")) local({
  d <- local({
    f <- grep("^--file=", commandArgs(), value = TRUE)[1]
    if (is.na(f)) getwd() else dirname(normalizePath(sub("^--file=", "", f), mustWork = FALSE))
  })
  for (p in c(file.path(d, "diann_q_columns.R"), "diann_q_columns.R"))
    if (file.exists(p)) { source(p); return(invisible()) }
  stop("diann_q_columns.R not found next to build_maxlfq.R -- it defines the ",
       "identification-FDR columns and there is no safe default to guess.")
})
# ONE definition of the contaminant filter (shared with run_de.R's dpc path).
if (!exists("CONTAMINANT_TAG")) local({
  d <- local({
    f <- grep("^--file=", commandArgs(), value = TRUE)[1]
    if (is.na(f)) getwd() else dirname(normalizePath(sub("^--file=", "", f), mustWork = FALSE))
  })
  for (p in c(file.path(d, "contaminants.R"), "contaminants.R"))
    if (file.exists(p)) { source(p); return(invisible()) }
  stop("contaminants.R not found next to build_maxlfq.R -- it defines the contaminant filter.")
})

# What the report says its quantity is: parquet schema metadata "delimp.quantity.<key>", written by
# the adapter that made the report (run_search.py QUANTITY_META: level, label, value,
# identification, q_columns = placeholder | real, citation, caveat). Sage declares level "peptide"; FragPipe-DDA, AlphaDIA and
# Radiant declare "protein". DIA-NN's own report declares nothing -- its PG.MaxLFQ IS DIA-NN's
# MaxLFQ, one value per protein and run. NULL when nothing is declared.
declared_quantity <- function(ds) {
  m <- tryCatch(ds$schema$metadata, error = function(e) NULL)
  k <- grep("^delimp\\.quantity\\.", names(m), value = TRUE)
  if (!length(k)) return(NULL)
  stats::setNames(lapply(k, function(x) m[[x]]), sub("^delimp\\.quantity\\.", "", k))
}

# THE description of this path's quantification (architectural rule 1): methods.txt,
# de_provenance.json, make_methods.py, the report, the AI brief and the reproducibility log all
# read it; none of them names a quantification of its own. A report that declares PEPTIDE-level
# rows (Sage) is rolled up below as the highest row per protein and run, and says exactly that.
# `cq`: diann_cont_quant_exclude() for the report (contaminants.R), read only for an undeclared
# (DIA-NN) report -- what its PG.MaxLFQ left out is said from how DIA-NN actually ran.
maxlfq_descriptor <- function(declared = NULL, cq = NULL) {
  de_engine <- "limma::lmFit -> contrasts.fit -> eBayes (NA-tolerant per row)"
  missing_policy <- paste("NAs left in place; limma drops them per row.",
                          "All-missing-in-one-condition proteins are on/off calls.")
  if (!is.null(declared) && identical(declared$level, "peptide")) {
    return(list(
      pipeline_id    = "peptide_max",
      display_label  = "Highest-peptide intensity + limma",
      # each protein's value is computed here from the peptide rows the contaminant filter
      # kept, so a kept row (a keratin sample's keratin) can carry it: the filter (contaminant
      # step, above the rollup in build_maxlfq) runs before max()
      requantified_from_precursors = TRUE,
      kept_contaminant_quant = sprintf(paste(
        "Each protein's value is taken here from the peptide rows the filter kept, so those",
        "keratin peptides count in it; %s was given no contaminant exclusion."),
        declared$label %||% "the search engine"),
      kept_keratin_under_quantified = FALSE,
      rollup_method  = sprintf(paste("protein intensity = the highest %s per protein and run",
                                     "(one peptide carries each value; no protein rollup model)"),
                               declared$value %||% "peptide-level intensity"),
      identification_fdr = declared$identification,
      q_columns_role = declared$q_columns,
      de_engine      = de_engine,
      missing_policy = missing_policy,
      caveat = paste(
        "Protein quantities are cruder than a protein rollup model: each protein's value in a",
        "run is its single most intense peptide, so different peptides can carry it in",
        "different runs, one interfered peptide can set it, and the other peptides' evidence",
        "is not used. Treat fold changes as indicative and check key proteins on their",
        "peptides."),
      plain_language = paste(
        "**Highest-peptide intensity + limma**: each protein's amount in a sample is the",
        "intensity of its most intense peptide in the search engine's label-free",
        "quantification -- a simple stand-in, not a model combining all of its peptides; limma",
        "then tests for differences with empirical-Bayes-moderated statistics. Missing values",
        "are left in place and handled per protein."),
      citation = sprintf(paste("Quantification: %s; each protein's value is its most intense",
                               "peptide per run. DE: limma (Ritchie et al. 2015, NAR 43:e47)."),
                         declared$citation %||% "the search engine's peptide LFQ")
    ))
  }
  if (!is.null(declared) && identical(declared$level, "protein")) {
    # The engine's own protein quantity, one value per protein and run: passed through as the
    # engine reported it (max() over one value is that value), so no rollup is described.
    label <- declared$label %||% "Search-engine protein quantities"
    value <- declared$value %||% "the search engine's protein quantity"
    return(list(
      pipeline_id    = "engine_protein",
      display_label  = sprintf("%s + limma", label),
      # the engine's own protein quantity as reported: the contaminant filter decides which
      # rows reach the matrix, never what the engine put in a quantity
      requantified_from_precursors = FALSE,
      # Not FALSE: the engine was given no contaminant exclusion, but its protein inference (e.g.
      # FragPipe/Philosopher's razor peptides) can still have moved peptides shared with those
      # entries to another protein -- nothing here records whether it did (sage-review N4)
      kept_contaminant_quant = sprintf(paste(
        "Their quantities are %s's own protein quantities, used as reported (not re-derived",
        "here); it was given no contaminant exclusion, and whether its protein inference moved",
        "peptides shared with those entries to another protein is not recorded."), label),
      kept_keratin_under_quantified = NA,
      rollup_method  = sprintf(paste("protein intensity = %s, as the search engine reported it",
                                     "(one value per protein and run; no re-rollup here)"), value),
      identification_fdr = declared$identification,
      q_columns_role = declared$q_columns,
      de_engine      = de_engine,
      missing_policy = missing_policy,
      caveat         = declared$caveat,
      plain_language = sprintf(paste(
        "**%s + limma**: the search engine's own protein quantities (%s) are compared between",
        "groups with limma's empirical-Bayes-moderated statistics. Missing values are left in",
        "place and handled per protein."), label, value),
      citation = sprintf("Quantification: %s. DE: limma (Ritchie et al. 2015, NAR 43:e47).",
                         declared$citation %||% "the search engine's protein quantification")
    ))
  }
  .dk <- diann_kept_quant(cq %||% list(recorded = FALSE, tag = NA_character_,
                                       source = "no report was named to read it from"),
                          requantified = FALSE)
  list(
    pipeline_id    = "maxlfq",
    display_label  = "MaxLFQ + limma",
    # the engine's protein quantity (PG.MaxLFQ) as reported: the contaminant filter decides
    # which rows reach the matrix, never what the engine put in a quantity
    requantified_from_precursors = FALSE,
    kept_contaminant_quant = .dk$text,
    kept_keratin_under_quantified = .dk$under_quantified,
    rollup_method  = "DIA-NN PG.MaxLFQ",
    de_engine      = de_engine,
    missing_policy = missing_policy,
    citation       = "Quantification: DIA-NN MaxLFQ (Demichev et al. 2020, Nat Methods 17:41). DE: limma (Ritchie et al. 2015, NAR 43:e47)."
  )
}

build_maxlfq <- function(report_path, format = "parquet", q_cutoff = 0.01,
                         eq_cutoff = 0, pgq_cutoff = 0, keep_runs = NULL,
                         drop_contaminants = TRUE, contaminant_exempt = NULL) {
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
  declared <- if (identical(format, "parquet")) declared_quantity(ds) else NULL
  per_peptide <- identical(declared$level, "peptide")

  needed   <- c("Run", "Protein.Group", "PG.MaxLFQ", DIANN_FDR_REQUIRED)
  # DIANN_FDR_OPTIONAL is optional only because older reports lack those columns,
  # not because applying them is discretionary. They MUST be listed here as well
  # as filtered on: filtering an arrow dataset on a column select() dropped
  # returns ZERO ROWS silently -- the defect that broke MaxLFQ in the DE-LIMP app.
  # Precursor.Id, the accession columns and the share intensity feed the contaminant
  # census and filter (contaminants.R); an adapted protein-level report lacks most of them.
  optional <- c("Empirical.Quality", "PG.MaxLFQ.Quality", "Genes", "Protein.Names",
                DIANN_FDR_OPTIONAL, CONTAMINANT_ID_COLUMNS, CONTAMINANT_FEATURE_COLUMNS,
                CONTAMINANT_SHARE_COLUMNS)
  miss <- setdiff(needed, cols)
  if (length(miss)) stop("MaxLFQ: missing required columns: ", paste(miss, collapse = ", "))

  sel <- c(needed, intersect(optional, cols))
  flt <- if (identical(format, "parquet")) dplyr::select(ds, dplyr::all_of(sel)) else ds[, sel]
  filters_applied <- character(0)

  # Run-level + library q-values do NOT control FDR across an experiment: a union of
  # run-level-passing IDs over many runs sits well above the nominal cutoff, and the
  # inflated protein list also inflates the family size m that BH corrects over.
  # Add DIA-NN's protein-level and experiment-wide q-values when present.
  q_columns <- character(0)
  if (!is.na(q_cutoff) && q_cutoff > 0) {
    flt <- dplyr::filter(flt, Q.Value <= !!q_cutoff,
                              Lib.Q.Value <= !!q_cutoff,
                              Lib.PG.Q.Value <= !!q_cutoff)
    q_columns <- DIANN_FDR_REQUIRED
    # Derived, not restated: this label was a third hand-written copy of the
    # required set and would have gone stale the moment the set changed.
    .fdr <- sub("\\.Value$", "", DIANN_FDR_REQUIRED)
    for (.qc in DIANN_FDR_OPTIONAL) {
      if (.qc %in% cols) {
        .cut <- diann_cutoff_for(.qc, q_cutoff)
        flt <- dplyr::filter(flt, .data[[.qc]] <= !!.cut)
        # Label the column with its own cutoff when it differs, so the recorded
        # provenance says what actually ran rather than implying one uniform value.
        .fdr <- c(.fdr, if (isTRUE(all.equal(.cut, q_cutoff))) .qc
                        else sprintf("%s@%.3f", .qc, .cut))
        q_columns <- c(q_columns, .qc)
      }
    }
    .qtxt <- sprintf("%s <= %.3f", paste(.fdr, collapse = "/"), q_cutoff)
    # an adapted report whose q-columns are 0.0 placeholders (its engine filtered upstream, and
    # it says so: declared_quantity() q_columns) -- the filter ran and kept every row; say that
    if (identical(declared$q_columns, "placeholder"))
      .qtxt <- paste(.qtxt, "(placeholder q-value columns, all 0.0, in this adapted report: it",
                     "keeps every row; the engine's own FDR is stated in the methods)")
    filters_applied <- c(filters_applied, .qtxt)
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

  # Contaminants: counted on the rows that passed the filters above, then removed before
  # the matrix is built -- so they never reach quantile normalisation, lmFit or BH.
  # `contaminant_exempt`: a keratin sample's keratin-family contaminant accessions, which stay in
  # (contaminants.R keratin_exemption: list(shared, alone)); n_exempt counts what that kept. tag counts which tags
  # (Cont_, FragPipe's contam_) the removed rows carried; flagged_ids are their accession lists,
  # one per precursor, for keratin_unchecked().
  # What one counted item is (contaminants.R contaminant_unit, from the declared level): DIA-NN
  # precursors, Sage peptides (Peptide.Id), a protein-level adapter's own groups. Counted
  # distinct, never item x run rows.
  .fcol <- contaminant_feature_column(names(rows))
  .unit <- contaminant_unit(declared$level, has_feature = !is.na(.fcol))
  cont <- list(census = NULL, share = NULL, id_column = contaminant_id_column(names(rows)),
               intensity_column = NA_character_, n_exempt = 0L, tag = CONTAMINANT_TAG,
               flagged_ids = character(0), unit = .unit)
  if (!is.na(cont$id_column)) {
    .flag <- is_contaminant(rows[[cont$id_column]], exempt = contaminant_exempt)
    .feat <- if (!is.na(.fcol)) rows[[.fcol]]
             else if (identical(.unit$level, "protein")) rows$Protein.Group else NULL
    if (length(contaminant_exempt$shared)) {
      .kept <- is_contaminant(rows[[cont$id_column]]) & !.flag
      cont$n_exempt <- if (is.null(.feat)) sum(.kept) else length(unique(.feat[.kept]))
    }
    cont$census <- contaminant_census(group = rows$Protein.Group, is_cont = .flag,
                                      feature = .feat, genes = rows$Genes,
                                      exempt = contaminant_exempt, unit = .unit)
    cont$tag <- contaminant_tag_label(contaminant_tag_counts(rows[[cont$id_column]], .flag, .feat))
    cont$flagged_ids <- .flagged_ids(rows[[cont$id_column]], .flag, .feat)
    .sc <- intersect(CONTAMINANT_SHARE_COLUMNS, names(rows))[1]
    if (!is.na(.sc)) {
      cont$intensity_column <- .sc
      cont$share <- contaminant_share(rows[[.sc]], .flag, run = rows$Run)
    }
    if (any(.flag)) {
      if (drop_contaminants) {
        rows <- rows[!.flag, , drop = FALSE]
        filters_applied <- c(filters_applied, sprintf(
          "contaminants removed: %s mapping to a %s entry (%s); %d %s protein groups",
          count_of(cont$census$n_precursors, .unit),
          cont$tag, cont$id_column, cont$census$n_groups_contaminant, cont$tag))
        if (!nrow(rows)) stop("MaxLFQ: every row maps to a contaminant entry.")
      } else {
        filters_applied <- c(filters_applied, sprintf(
          "contaminants kept (--keep-contaminants): %d %s protein groups",
          cont$census$n_groups_contaminant, cont$tag))
      }
    }
  }

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
  # A peptide-level report has many rows per cell BY DECLARATION: max() is its rollup, and the
  # descriptor says so. Only an undeclared report's extra rows are a surprise worth a warning.
  if (.nbad > 0 && !per_peptide)
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
    contaminants = cont,
    # The q-value columns actually filtered on. Recorded rather than assumed so the
    # emitted reproducibility script names the same columns this run used -- older
    # reports lack PG.Q.Value / Global.*, and hard-coding them would emit a script
    # that filters on a column the run never had.
    q_columns = q_columns,
    descriptor = maxlfq_descriptor(declared,
                                   if (is.null(declared)) diann_cont_quant_exclude(report_path)),
    # what the report declared about itself, verbatim, for the record
    declared_quantity = declared
  )
}
