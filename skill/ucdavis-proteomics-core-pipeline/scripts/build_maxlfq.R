#!/usr/bin/env Rscript
# build_maxlfq.R  --  Port of DE-LIMP's build_maxlfq_pipeline() (R/helpers.R).
# Reads a DIA-NN report, applies ID-FDR + optional QuantUMS filters, pivots
# PG.MaxLFQ to a protein x run matrix and log2-transforms it. A between-run step on that matrix
# is decided by ONE rule, maxlfq_between_run() (run_de.R --between-run, default auto): none for
# a DIA-NN PG.MaxLFQ DIA-NN already normalised (since 2.14.2: a quantile step on top of it made
# 61% of the significant calls wrong on the PXD063416 4 m/z benchmark) or an engine's quantities
# the engine normalised; a location-only median-ratio step for quantities no engine normalised
# (a --no-norm search's PG.MaxLFQ, Sage's LFQ; 2.15.1); nothing with quantities = "raw";
# quantile, median, median_ratio or none on request.
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
# `rollup`: how a peptide-level report's rows became one value per protein and run (maxlfq_matrix
# says which ran: "median_polish", or "max" for a report with no peptide column). `diann_shaped`:
# an undeclared report that has DIA-NN's precursor columns (Precursor.Id, Precursor.Quantity) --
# only then is its PG.MaxLFQ described as DIA-NN's (lane3 #16: a hand-made file was).
# The missing-value policy here says nothing about on/off proteins: run_de.R adds that sentence
# only when OnOff_Proteins.csv was actually written (onoff record), never from this text.
maxlfq_descriptor <- function(declared = NULL, cq = NULL, rollup = "median_polish",
                              diann_shaped = TRUE) {
  de_engine <- "limma::lmFit -> contrasts.fit -> eBayes (NA-tolerant per row)"
  missing_policy <- "NAs left in place; limma drops them per row."
  if (!is.null(declared) && identical(declared$level, "peptide") && identical(rollup, "median_polish")) {
    return(list(
      pipeline_id    = "peptide_medpolish",
      display_label  = "Peptide median polish + limma",
      requantified_from_precursors = TRUE,
      kept_contaminant_quant = sprintf(paste(
        "Each protein's value is computed here from the peptide rows the filter kept, so those",
        "keratin peptides count in it; %s was given no contaminant exclusion."),
        declared$label %||% "the search engine"),
      kept_keratin_under_quantified = FALSE,
      contaminant_rule_origin = contaminant_rule_origin(declared$label %||% "search engine's"),
      rollup_method  = sprintf(paste("protein intensity = Tukey median polish of the log2 %s",
                                     "(peptides x runs, per protein): the run effect plus the",
                                     "overall level, in runs with at least one peptide value"),
                               declared$value %||% "peptide-level intensity"),
      identification_fdr = declared$identification,
      q_columns_role = declared$q_columns,
      de_engine      = de_engine,
      missing_policy = missing_policy,
      caveat = paste(
        "Protein quantities are a summary of the search engine's peptide quantities (median",
        "polish, robust to one interfered peptide), not a protein rollup model with",
        "uncertainty: a protein carried by one peptide rests on that peptide alone. Check key",
        "proteins on their peptides."),
      plain_language = paste(
        "**Peptide median polish + limma**: each protein's amount in a sample is a robust",
        "summary (Tukey's median polish) of its peptides' intensities from the search engine's",
        "label-free quantification; limma then tests for differences with",
        "empirical-Bayes-moderated statistics. Missing values are left in place and handled",
        "per protein."),
      citation = sprintf(paste("Quantification: %s; each protein's value is a Tukey median polish",
                               "of its peptides (as in MSstats, Choi et al. 2014, Bioinformatics",
                               "30:2524). DE: limma (Ritchie et al. 2015, NAR 43:e47)."),
                         declared$citation %||% "the search engine's peptide LFQ")
    ))
  }
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
      contaminant_rule_origin = contaminant_rule_origin(declared$label %||% "search engine's"),
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
      contaminant_rule_origin = contaminant_rule_origin(declared$label %||% "search engine's"),
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
  if (!isTRUE(diann_shaped)) {
    # Neither declared by an adapter nor shaped like a DIA-NN report: nothing says what its
    # PG.MaxLFQ column is, so it is not called DIA-NN's MaxLFQ (lane3 #16)
    return(list(
      pipeline_id    = "undeclared_pg_maxlfq",
      display_label  = "PG.MaxLFQ column of an undeclared report + limma",
      requantified_from_precursors = FALSE,
      kept_contaminant_quant = paste("The report's own protein quantities, used as given; how they",
                                     "were made is NOT RECORDED [confirm]."),
      kept_keratin_under_quantified = NA,
      contaminant_rule_origin = contaminant_rule_origin(),
      rollup_method  = paste("the report's PG.MaxLFQ column as given; the report declares no",
                             "quantity and has no DIA-NN precursor columns, so which tool made it",
                             "and how is NOT RECORDED [confirm]"),
      de_engine      = de_engine,
      missing_policy = missing_policy,
      caveat = paste("This report declares no quantity (no delimp.quantity metadata) and is not a",
                     "DIA-NN precursor report: the tool and method behind its PG.MaxLFQ values",
                     "are not recorded. Say which before using these results."),
      citation = "Quantification: NOT RECORDED [confirm]. DE: limma (Ritchie et al. 2015, NAR 43:e47)."
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
    contaminant_rule_origin = contaminant_rule_origin(),
    rollup_method  = "DIA-NN PG.MaxLFQ",
    de_engine      = de_engine,
    missing_policy = missing_policy,
    citation       = paste("Quantification: DIA-NN PG.MaxLFQ (Demichev et al. 2020, Nat Methods 17:41),",
                           "DIA-NN's implementation of MaxLFQ (Cox et al. 2014, Mol Cell Proteomics",
                           "13:2513), on precursor quantities from QuantUMS in DIA-NN 1.9 and later",
                           "(Grossmann et al. 2026, Nat Biotechnol, doi:10.1038/s41587-026-03131-2).",
                           "DE: limma (Ritchie et al. 2015, NAR 43:e47).")
  )
}

# Is this DIA-NN report non-normalised -- searched with --no-norm, so its PG.MaxLFQ is too
# (DIA-NN README: "if normalisation is disabled, these quantities are actually non-normalised")?
# Read from the data, never assumed: every row's Precursor.Normalised equals its
# Precursor.Quantity. -> TRUE / FALSE, or NA when the report has neither column.
report_is_unnormalised <- function(ds, cols) {
  if (!all(c("Precursor.Quantity", "Precursor.Normalised") %in% cols)) return(NA)
  d <- dplyr::collect(dplyr::summarise(
    dplyr::filter(dplyr::select(ds, Precursor.Quantity, Precursor.Normalised),
                  Precursor.Quantity > 0, Precursor.Normalised > 0),
    worst = max(abs(Precursor.Normalised / Precursor.Quantity - 1), na.rm = TRUE)))
  isTRUE(d$worst[1] < 1e-6)
}

# Did the ENGINE normalise its declared quantities across runs? -> list(value = TRUE / FALSE / NA,
# source). Read from the declaration (run_search.py, key engine_normalised "true" / "false"), never
# from its prose. A peptide-level report declaring nothing is Sage's (the only peptide-level adapter):
# Sage's LFQ has no cross-run normalisation. Anything else that does not say: NA, not known.
engine_normalised <- function(declared) {
  if (is.null(declared)) return(list(value = NA, source = "not an adapted report"))
  v <- tolower(trimws(as.character(declared$engine_normalised %||% "")))
  if (v %in% c("true", "false"))
    return(list(value = identical(v, "true"),
                source = sprintf("the report's declaration (engine_normalised %s)", v)))
  # declared as not known (run_search.ENGINE_NORMALISED_UNKNOWN; Radiant, 2.15.1): never read as
  # Sage's, whatever its level
  if (identical(v, "unknown"))
    return(list(value = NA, source = paste("the report declares that whether its engine normalised",
                                           "is not known (engine_normalised unknown)")))
  if (identical(declared$level, "peptide"))
    return(list(value = FALSE, source = paste("a peptide-level report (Sage's LFQ, which has no",
                                              "cross-run normalisation); declared before 2.15.1")))
  list(value = NA, source = "the report does not declare whether its engine normalised (engine_normalised)")
}

# THE between-run rule of the MaxLFQ route (architectural rule 3): which step, if any, runs on
# the log2 protein x run matrix, and why. The descriptor's `between_run` and its `normalisation`
# text (maxlfq_normalisation_text), run_de.R's --between-run and the reproducibility log all read
# this answer; none of them decides it again.
#   between_run  "auto" (default) | "none" | "quantile" | "median" | "median_ratio"
#   quantities   "normalised" | "raw"                                (run_de.R --quantities)
#   declared     declared_quantity() of the report; NULL for a DIA-NN report
#   dia_nn_norm  a DIA-NN report's list(no_norm = TRUE / FALSE / NA, source); NULL otherwise
# auto asks ONE question of every report: did the engine already normalise these quantities
# across runs (TRUE / FALSE / not known)? For DIA-NN the answer is "not --no-norm"
# (dia_nn_norm); for an adapted report its declaration (engine_normalised()).
#   yes (DIA-NN's PG.MaxLFQ; IonQuant      none -- no second step (2.14.2, I20: a quantile step on
#   normalization=1, directLFQ)            top of DIA-NN's gave 61% wrong calls, PXD063416 4 m/z)
#   no (a --no-norm search's PG.MaxLFQ;    median_ratio, a location-only step (2.15.1; up to 2.15.0
#   Sage's LFQ)                            quantile -- on real --no-norm PXD063416 data quantile
#                                          gave 5.5/46/72/69% wrong calls at 2/4/6/10 m/z,
#                                          median_ratio 3.9/3.7/11.6/3.3%: measured on HIVE, 2026-11-17)
#   not known, DIA-NN                      none, tagged [confirm]: DIA-NN normalises unless it is
#                                          run with --no-norm, and the report's own columns are
#                                          read first (report_is_unnormalised)
#   not known, an adapted report           median_ratio, tagged [confirm]: on quantities an engine
#                                          already normalised it is close to a no-op
# A record made before 2.15.1 keeps the step it ran: provenance.py's replay passes the applied
# method (between_run.method), and the normalisation check's gate compares applied steps.
# raw: none; on an adapted report only when its engine did NOT normalise (Sage: its quantities as
# Sage measured them -- the non-normalised route the normalisation check offers, audit L2-2).
# Any step asked for on raw quantities is refused.
# median = subtract each run's median log2, then add back the median of the run medians.
# median_ratio = maxlfq_median_ratio_factors(): each run's median log2 ratio to the protein's
# median across runs, over the proteins quantified in every run. Chosen over quantile and median
# for engine quantities that need a step: a location-only step keeps each run's distribution
# shape, and fitted on proteins quantified in every run it is not moved by missing values. It IS
# moved by proteins that change in one direction, in proportion to how many: the group bias of
# its factors (simulated, 4 vs 4, 6 seeds; 2.15.1 review exp1.R) is 0.00 log2 for a null, +0.06
# with 10% of proteins up 1.5 log2, +0.19 with 25% up, +0.47 with 40% up 2 log2, -0.18 with 30%
# down -- never zero once proteins change one way (each protein's median across runs sits
# between the groups); hence the lopsided-volcano guard. On simulated data (8 seeds each, 4 vs 4,
# 2,000 proteins, intensity-dependent missing values) its wrong calls among the significant ones were 5.1-12.5% where the groups were
# loaded differently (quantile 9.7%, no step 76-93%), 7.2% with 10% of proteins up (quantile 8.5%)
# and 25% with 25% up 2.8-fold (quantile 49%); with 40% up 4-fold it is as wrong as quantile
# (53%), which the normalisation check's lopsided-volcano test exists to stop. Measured on the
# PXD063416 benchmark's DIA-NN --no-norm quantities (2.15.1, above; references/de-analysis.md):
# better than quantile and no step at 4/6/10 m/z, but 11.6% wrong at 6 m/z, its weak window.
# -> list(method = "none" | "quantile" | "median" | "median_ratio", why, requested = between_run,
#         engine_normalised = engine_normalised() for an adapted report)
BETWEEN_RUN_CHOICES <- c("auto", "none", "quantile", "median", "median_ratio")
maxlfq_between_run <- function(between_run = "auto", quantities = "normalised", declared = NULL,
                               dia_nn_norm = NULL) {
  if (!is.character(between_run) || length(between_run) != 1 ||
      !between_run %in% BETWEEN_RUN_CHOICES)
    stop("--between-run must be one of ",
         paste(paste(head(BETWEEN_RUN_CHOICES, -1), collapse = ", "), "or",
               tail(BETWEEN_RUN_CHOICES, 1)),
         " (got ", paste(format(between_run), collapse = " "), ")", call. = FALSE)
  en <- if (!is.null(declared)) engine_normalised(declared)
  if (identical(quantities, "raw")) {
    if (!between_run %in% c("auto", "none"))
      stop(sprintf(paste(
        "--between-run %s with --quantities raw: raw quantities are the non-normalised ones, and",
        "a between-run step would normalise them. Either --quantities normalised --between-run",
        "%s, or --quantities raw with no --between-run."), between_run, between_run), call. = FALSE)
    return(list(method = "none", requested = between_run, engine_normalised = en,
                why = if (is.null(declared))
                  "--quantities raw: a --no-norm search's PG.MaxLFQ, with no between-run step"
                else paste0("--quantities raw: ", declared$label %||% "the engine",
                            "'s quantities as it measured them (it applies no cross-run ",
                            "normalisation), with no between-run step")))
  }
  if (!identical(between_run, "auto"))
    return(list(method = between_run, requested = between_run, engine_normalised = en,
                why = sprintf("--between-run %s (chosen for this DE)", between_run)))
  # THE question (above): did the engine already normalise these quantities across runs?
  normed <- if (!is.null(declared)) en$value else !as.logical(dia_nn_norm$no_norm %||% NA)
  dia <- is.null(declared)
  # a DIA-NN report's answer carries no engine_normalised (its source is dia_nn_normalisation)
  if (isTRUE(normed))
    return(c(list(method = "none", requested = "auto"), list(engine_normalised = en)[!dia],
             list(why = if (dia) paste(
      "PG.MaxLFQ already carries DIA-NN's cross-run normalisation; no second step (a quantile",
      "step on top gave 61% wrong calls on the PXD063416 4 m/z benchmark)") else paste(
      "the engine already normalised its quantities across runs; no second step (a quantile",
      "step on top of already-normalised quantities gave 61% wrong calls on the PXD063416",
      "4 m/z benchmark)"))))
  if (isFALSE(normed))
    return(c(list(method = "median_ratio", requested = "auto"), list(engine_normalised = en)[!dia],
             list(why = if (dia)
      paste("the search ran with --no-norm, so PG.MaxLFQ carries no cross-run normalisation: a",
            "location-only step, each run centred on its median log2 ratio over the proteins",
            "quantified in every run (2.15.1; quantile up to 2.15.0)")
      else paste("the engine applies no cross-run normalisation, so a location-only step: each run",
                 "centred on its median log2 ratio over the proteins quantified in every run"))))
  if (dia)
    return(list(method = "none", requested = "auto", why = paste(
      "DIA-NN normalises by default; whether it did for this report is NOT RECORDED [confirm]")))
  list(method = "median_ratio", requested = "auto", engine_normalised = en, why = paste(
    "a location-only step (median log2 ratio over the proteins quantified in every run);",
    "whether the engine already normalised is NOT RECORDED [confirm] -- if it did, this step",
    "is close to a no-op"))
}

# The proteins the median-ratio step is computed on: those quantified in every run (no missing-
# value pattern can move them), or -- with fewer than MEDIAN_RATIO_MIN_PROTEINS of those -- the
# proteins quantified in at least half the runs, ratios taken over the runs that have them.
MEDIAN_RATIO_MIN_PROTEINS <- 50L
# E: log2 protein x run matrix -> list(factors (one per run, subtracted from it), n_proteins, basis)
maxlfq_median_ratio_factors <- function(E) {
  nobs <- rowSums(is.finite(E))
  full <- nobs == ncol(E)
  basis <- "proteins quantified in every run"
  use <- full
  if (sum(full) < MEDIAN_RATIO_MIN_PROTEINS) {
    use <- nobs >= max(2, ceiling(ncol(E) / 2))
    basis <- sprintf("proteins quantified in at least half the runs (only %d in every run)", sum(full))
  }
  if (sum(use) < 10)
    stop("MaxLFQ: the median-ratio step needs at least 10 proteins quantified in most runs; ",
         sum(use), " are. Use --between-run none or quantile, and say why.", call. = FALSE)
  X <- E[use, , drop = FALSE]
  R <- X - apply(X, 1, stats::median, na.rm = TRUE)
  f <- apply(R, 2, stats::median, na.rm = TRUE)
  list(factors = f, n_proteins = sum(use), basis = basis)
}

# The Methods' words for the between-run normalisation the route applied (descriptor
# `normalisation`): what the quantities carried before this step, then the step -- built from
# maxlfq_between_run()'s answer `br`, never decided here.
maxlfq_normalisation_text <- function(br, quantities, declared = NULL, dia_nn_norm = NULL) {
  if (identical(quantities, "raw"))
    return(if (is.null(declared))
      "none: PG.MaxLFQ of a DIA-NN report searched with --no-norm, no between-run step"
    else paste0("none: ", declared$label %||% "the engine", "'s quantities as it measured them ",
                "(no cross-run normalisation by the engine, ", br$engine_normalised$source %||%
                  "source not recorded", "), no between-run step here"))
  step <- switch(br$method,
                 quantile = "quantile normalisation (limma::normalizeBetweenArrays)",
                 median = paste("median centring (each run's log2 values shifted so that its",
                                "median is the median of the run medians)"),
                 median_ratio = paste0("median-ratio centring (each run shifted by its median log2 ",
                                       "ratio to the protein's median across runs, over the ",
                                       br$basis %||% "proteins quantified in every run",
                                       if (!is.null(br$n_proteins)) sprintf(", %d proteins", br$n_proteins) else "",
                                       "; a location-only step)"),
                 none = NULL)
  nn <- dia_nn_norm$no_norm
  # what the engine itself did across runs, as its config records it (run_search.py declares
  # `engine_normalisation` for FragPipe IonQuant and AlphaDIA): said before this DE's step, so
  # the Methods never leave out the first of two normalisations (C2 statistician M2)
  eng <- if (!is.null(declared)) declared$engine_normalisation else NULL
  eng <- if (is.character(eng) && length(eng) == 1 && nzchar(eng)) eng else NULL
  txt <- if (!is.null(declared) && !is.null(eng)) {
    if (is.null(step)) paste0(eng, "; no between-run step here")
    else paste0(eng, ", then ", step, " of the engine's quantities")
  } else if (!is.null(declared)) {
    if (is.null(step)) "none: the engine's quantities as it reported them, with no between-run step here"
    else paste(step, "of the engine's quantities")
  } else if (isTRUE(nn)) {
    if (is.null(step)) paste("none: PG.MaxLFQ of a DIA-NN search run with --no-norm and no",
                             "between-run step here, so no between-run normalisation at all")
    else paste(step, "only: the DIA-NN search ran with --no-norm, so PG.MaxLFQ carries no",
               "DIA-NN cross-run normalisation")
  } else if (isFALSE(nn)) {
    if (is.null(step)) "DIA-NN's cross-run normalisation (PG.MaxLFQ already carries it); no second step"
    else paste0("DIA-NN's cross-run normalisation (PG.MaxLFQ), then ", step)
  } else {
    if (is.null(step)) paste("PG.MaxLFQ as DIA-NN reported it, no second step; whether DIA-NN",
                             "normalised it is NOT RECORDED [confirm] (DIA-NN normalises by",
                             "default)")
    else paste(step, "of PG.MaxLFQ; whether DIA-NN also normalised it is NOT RECORDED [confirm]")
  }
  if (!identical(br$requested, "auto"))
    txt <- sprintf("%s (--between-run %s, chosen for this DE)", txt, br$requested)
  else if (!is.null(declared) && is.null(eng) && !is.null(br$engine_normalised) &&
           is.na(br$engine_normalised$value) && !is.null(step))
    txt <- paste(txt, "(whether the engine normalised is NOT RECORDED [confirm])")
  txt
}

# A peptide-level report's rows -> one log2 value per (protein, run): Tukey's median polish of each
# protein's peptide x run log2 matrix (stats::medpolish), the overall level plus the run effect, in
# the runs where at least one of its peptides has a value. Robust to one interfered peptide, and
# every peptide counts -- unlike the single most intense peptide per run 2.15.0 used, which could
# be a different peptide in each run (audit L2-4). `feature`: the peptide column.
# -> list(E = log2 protein x run matrix, prot_ids)
peptide_medpolish <- function(rows, feature) {
  d <- data.frame(p = as.character(rows$Protein.Group), f = as.character(rows[[feature]]),
                  r = as.character(rows$Run), v = suppressWarnings(log2(rows$PG.MaxLFQ)),
                  stringsAsFactors = FALSE)
  d <- d[is.finite(d$v), , drop = FALSE]
  runs <- sort(unique(as.character(rows$Run)))
  prots <- sort(unique(d$p))
  E <- matrix(NA_real_, length(prots), length(runs), dimnames = list(prots, runs))
  for (pr in split(d, d$p)) {
    m <- tapply(pr$v, list(pr$f, factor(pr$r, levels = runs)), max)
    if (nrow(m) == 1) { E[pr$p[1], ] <- m[1, ]; next }
    mp <- stats::medpolish(m, na.rm = TRUE, trace.iter = FALSE)
    val <- mp$overall + mp$col
    val[colSums(is.finite(m)) == 0] <- NA_real_
    E[pr$p[1], ] <- val
  }
  list(E = E, prot_ids = prots)
}

# rows (one per precursor / peptide, after the identification, quality and contaminant filters) ->
# the log2 protein x run matrix, with the between-run step `br` (maxlfq_between_run) applied.
# ONE definition, used by the run itself and by the QuantUMS trade-off (every cutoff's matrix is
# built by this function, so a trade-off row at the cutoffs actually used equals the DE's matrix).
# per_peptide: a peptide-level report -- median polish (peptide_medpolish) when it has a peptide
# column, else the highest row per cell (what 2.15.0 did, said by the descriptor). Returns also
# `rollup` (which ran), `br` (with the median-ratio step's basis filled in) and `identified`: for a
# report carrying Identified.In.Run (Sage, 2.15.1), a protein x run logical -- TRUE where a peptide
# of the protein was identified by MS2 in that run AND has MS1 signal there, FALSE where its value
# comes only from peptides measured by matching between runs (Sage integrated their MS1 signal at
# the aligned retention time; only their identity comes from other runs), NA where there is no
# value; NULL when the report does not say. A row with no MS1 signal (PG.MaxLFQ 0 or missing) is
# never a value: the rollups above drop it (log2 of 0 is not finite), so it is NA, never present.
maxlfq_matrix <- function(rows, br, quantities, per_peptide = FALSE) {
  filters_applied <- character(0)
  .fcol <- if (per_peptide) contaminant_feature_column(names(rows)) else NA_character_
  rollup <- if (per_peptide && !is.na(.fcol)) "median_polish" else if (per_peptide) "max" else "engine"
  if (identical(rollup, "median_polish")) {
    pm <- peptide_medpolish(rows, .fcol)
    E_pre <- pm$E; prot_ids <- pm$prot_ids
  } else {
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
  }
  E_pre[!is.finite(E_pre)] <- NA_real_

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

  # the step maxlfq_between_run() decided, and only that one: never a quiet substitute (a quantile
  # step without limma used to fall back to median centring while the record said quantile)
  if (identical(br$method, "quantile")) {
    if (!requireNamespace("limma", quietly = TRUE))
      stop("MaxLFQ: the quantile step (", br$why, ") needs limma: BiocManager::install('limma').",
           call. = FALSE)
    E <- limma::normalizeBetweenArrays(E_pre, method = "quantile")
  } else if (identical(br$method, "median")) {
    cm <- apply(E_pre, 2, stats::median, na.rm = TRUE)
    E  <- sweep(E_pre, 2, cm - stats::median(cm, na.rm = TRUE))
  } else if (identical(br$method, "median_ratio")) {
    mr <- maxlfq_median_ratio_factors(E_pre)
    E  <- sweep(E_pre, 2, mr$factors)
    br$basis <- mr$basis; br$n_proteins <- mr$n_proteins
    br$factors <- as.list(round(mr$factors, 4))
  } else {
    E <- E_pre                        # no between-run step here
  }
  if (identical(quantities, "raw"))
    filters_applied <- c(filters_applied, paste(
      "no between-run normalisation (--quantities raw: the engine's non-normalised quantities,",
      "no between-run step)"))
  rownames(E) <- prot_ids
  identified <- NULL
  if ("Identified.In.Run" %in% names(rows)) {
    # identified AND measured: an identified peptide with no MS1 signal in that run gives the
    # cell no value, so it cannot make the cell "identified" (sage_lfq_check.classify_rows)
    .sig <- suppressWarnings(as.numeric(rows$PG.MaxLFQ))
    idr <- rows[!is.na(rows$Identified.In.Run) & as.integer(rows$Identified.In.Run) == 1L &
                  is.finite(.sig) & .sig > 0, , drop = FALSE]
    identified <- matrix(FALSE, nrow(E), ncol(E), dimnames = dimnames(E))
    ii <- match(as.character(idr$Protein.Group), rownames(E)); jj <- match(as.character(idr$Run), colnames(E))
    ok <- !is.na(ii) & !is.na(jj)
    identified[cbind(ii[ok], jj[ok])] <- TRUE
    identified[is.na(E)] <- NA
  }
  list(E = E, prot_ids = prot_ids, filters_applied = filters_applied, rollup = rollup, br = br,
       identified = identified, feature = if (identical(rollup, "median_polish")) .fcol else NULL)
}

# THE coverage rule of the MaxLFQ route (run_de.R's filter, the QuantUMS trade-off's refits and the
# reproducibility log all use it): a protein is TESTED when it is quantified in at least
# max(2, ceiling(frac x all samples)) samples -- 2.15.0's rule, kept so nothing it tested is lost --
# OR in at least max(2, ceiling(frac x n)) of the n runs of ANY ONE group. The second half is 2.15.1's
# (audit L2-1): counted over all samples only, a protein seen in every run of a small group and
# nowhere else (a bait-specific interactor) was deleted from every output, e.g. 3 of 3 vs 0 of 5.
# Such a protein now reaches limma: a contrast where it has values on both sides is tested, one
# where a side is empty gives no fold change (Evidence "presence call"; OnOff_Proteins.csv lists the
# simple A-B ones). `groups`: the group of each column of E, in order.
# Such a protein stays in the ONE lmFit/eBayes fit, so its within-group variance feeds limma's
# variance prior (2.15.1, HIVE F3: kept on purpose). It is a real replicate variance, as for any
# protein missing from a group, which limma has always pooled; with 3+ groups the protein IS
# tested between the groups that have it, and those contrasts need the prior every other protein
# gets -- a separate fit would give them another. Measured: on PXD063416 these 28-88 rows moved
# adj.P by at most 0.0018 and changed no call in the default routes (measured on HIVE, 2026-11-17).
# -> list(keep (logical per row), min_all, min_group (named per group), n_by_rule, text)
coverage_rule <- function(E, groups, frac) {
  g <- as.character(groups)
  stopifnot(length(g) == ncol(E))
  obs <- !is.na(E)
  min_all <- max(2L, as.integer(ceiling(frac * ncol(E))))
  lv <- sort(unique(g))
  n_g <- vapply(lv, function(x) sum(g == x), 0L)
  min_g <- vapply(n_g, function(n) max(2L, as.integer(ceiling(frac * n))), 0L)
  by_all <- rowSums(obs) >= min_all
  by_grp <- if (length(lv)) Reduce(`|`, lapply(lv, function(x)
    rowSums(obs[, g == x, drop = FALSE]) >= min_g[[x]])) else rep(FALSE, nrow(E))
  keep <- by_all | by_grp
  list(keep = keep, min_all = min_all, min_group = as.list(min_g), n_group = as.list(n_g),
       n_kept = sum(keep), n_dropped = sum(!keep),
       n_kept_by_group_only = sum(keep & !by_all),
       text = sprintf(paste0("coverage: quantified in >= %d/%d samples (%.0f%%), or in >= %.0f%% ",
                             "(at least 2) of the runs of any one group (%s); %d proteins dropped, ",
                             "%d kept by the per-group rule only"),
                      min_all, ncol(E), 100 * frac, 100 * frac,
                      paste(sprintf("%s >= %d/%d", lv, min_g, n_g), collapse = ", "),
                      sum(!keep), sum(keep & !by_all)))
}

# `no_norm`: list(no_norm = TRUE / FALSE / NULL, source) -- whether the search ran with --no-norm,
# from its provenance (run_de.R asks normalization_check.py searched_no_norm). NULL no_norm: not
# recorded, and a normalised MaxLFQ then reads it from the data (report_is_unnormalised).
# `between_run`: run_de.R --between-run, decided by maxlfq_between_run() (above).
build_maxlfq <- function(report_path, format = "parquet", q_cutoff = 0.01,
                         eq_cutoff = 0, pgq_cutoff = 0, keep_runs = NULL,
                         drop_contaminants = TRUE, contaminant_exempt = NULL,
                         quantities = "normalised", no_norm = NULL, between_run = "auto",
                         tradeoff = NULL) {
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
  # The between-run step (maxlfq_between_run, THE rule): a raw request is decided -- or refused --
  # before the report is read; a normalised one needs what PG.MaxLFQ carries (dia_nn_norm, below).
  br <- if (identical(quantities, "raw")) maxlfq_between_run(between_run, quantities, declared)
  # --quantities raw: no between-run normalisation anywhere. PG.MaxLFQ is DIA-NN's NORMALISED
  # protein quantity unless the search ran with --no-norm, and an adapted report's engine
  # quantities carry whatever their engine applied -- so only a non-normalised DIA-NN report
  # qualifies, checked from its data. Otherwise refuse: never "raw" in name only.
  if (identical(quantities, "raw") && !is.null(declared)) {
    # an adapted report: its quantities ARE the non-normalised ones when its engine applied no
    # cross-run normalisation (Sage's LFQ) -- the normalisation check's raw DE (audit L2-2).
    # An engine that normalised (IonQuant normalization=1, directLFQ), or one that does not say,
    # has no non-normalised quantity in this report: refused, with the way to get one.
    en <- engine_normalised(declared)
    if (!isFALSE(en$value))
      stop("--quantities raw on --method maxlfq: this report's quantities (",
           declared$label %||% "an adapted report", ") ",
           if (isTRUE(en$value)) "were already normalised across runs by the engine"
           else "may have been normalised across runs by the engine (not declared)",
           " (", en$source, "), so they are not raw. For FragPipe, search again with ",
           "ionquant.normalization=0 (make_presets.py) and run the raw DE on that report; ",
           "otherwise run the final DE on normalised quantities and say why.", call. = FALSE)
  } else if (identical(quantities, "raw")) {
    un <- report_is_unnormalised(ds, cols)
    if (!isTRUE(un))
      stop("--quantities raw on --method maxlfq: PG.MaxLFQ in this report is DIA-NN's ",
           "NORMALISED protein quantity (",
           if (is.na(un)) "the report has no Precursor.Quantity / Precursor.Normalised to check"
           else "its Precursor.Normalised differs from Precursor.Quantity",
           "). Next step, either: (1) --method dpc, which reads the non-normalised ",
           "Precursor.Quantity from this same report, no new search -- for step 8's ",
           "normalisation check, run its two DEs and the final DE with --method dpc; or (2) to ",
           "keep MaxLFQ, search again with --no-norm (diann_parallel.py --no-norm writes ",
           "no_norm_report.parquet), run this DE on that report, and give it to the check as ",
           "normalization_check.py run ... --report-raw <no_norm_report.parquet> ",
           "(SKILL.md step 8).", call. = FALSE)
  }

  # --quantities normalised on a DIA-NN report: what PG.MaxLFQ carries before any step here
  # -- DIA-NN's cross-run normalisation, or none when the search ran with --no-norm. The search's
  # provenance says which; only when it does not are the data read (a pass over two columns).
  # Said in the descriptor, never assumed (the DE once described a --no-norm report's MaxLFQ as
  # DIA-NN-normalised).
  dia_nn_norm <- NULL
  if (!identical(quantities, "raw") && is.null(declared)) {
    nn <- no_norm$no_norm
    src <- no_norm$source %||% "not recorded"
    if (is.null(nn)) {
      un <- report_is_unnormalised(ds, cols)
      if (!is.na(un)) {
        src <- paste0("read from the data: Precursor.Normalised ",
                      if (un) "equals" else "differs from", " Precursor.Quantity (", src, ")")
        nn <- un
      }
    }
    dia_nn_norm <- list(no_norm = if (is.null(nn)) NA else nn, source = src)
  }
  if (is.null(br)) br <- maxlfq_between_run(between_run, quantities, declared, dia_nn_norm)

  needed   <- c("Run", "Protein.Group", "PG.MaxLFQ", DIANN_FDR_REQUIRED)
  # DIANN_FDR_OPTIONAL is optional only because older reports lack those columns,
  # not because applying them is discretionary. They MUST be listed here as well
  # as filtered on: filtering an arrow dataset on a column select() dropped
  # returns ZERO ROWS silently -- the defect that broke MaxLFQ in the DE-LIMP app.
  # Precursor.Id, the accession columns and the share intensity feed the contaminant
  # census and filter (contaminants.R); an adapted protein-level report lacks most of them.
  optional <- c("Empirical.Quality", "PG.MaxLFQ.Quality", "Genes", "Protein.Names",
                "Identified.In.Run",
                DIANN_FDR_OPTIONAL, CONTAMINANT_ID_COLUMNS, CONTAMINANT_FEATURE_COLUMNS,
                CONTAMINANT_SHARE_COLUMNS)
  miss <- setdiff(needed, cols)
  if (length(miss)) stop("MaxLFQ: missing required columns: ", paste(miss, collapse = ", "))
  # Engine protein quantities (Ruling 53): MaxLFQ reads PG.MaxLFQ, and a report whose PG.MaxLFQ is
  # 0 or missing in every row -- a DIA-NN --dda first-pass report -- has none. It used to drop
  # every run and stop with "fewer than 2 usable runs survived the filters", which points at the
  # filters. Asked of dda_first_pass.py (pg_maxlfq_state: the ONE definition, which the search's
  # record uses too); its refusal text is that file's MAXLFQ_REFUSAL. Not answered: said, and run.
  .pq <- python_json("dda_first_pass.py", "pg-maxlfq",
                     c("--report", shQuote(normalizePath(report_path, mustWork = FALSE))),
                     "check whether the report has engine protein quantities (PG.MaxLFQ)")
  if (identical(.pq$res$pg_maxlfq, "absent")) stop(.pq$res$refusal, call. = FALSE)
  if (is.null(.pq$res$pg_maxlfq))
    message("[build_maxlfq] note: whether this report has engine protein quantities was not ",
            "checked (", .pq$why %||% .pq$res$why, ")")

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
  # QuantUMS trade-off (run_de.R): the rows are read ONCE, here, before the quality cutoffs, and every
  # cutoff set is applied in memory below -- the report is never re-read per cutoff. Needs both
  # QuantUMS columns; a report without them (Sage, FragPipe, adapted engines) gets a recorded reason.
  tradeoff_pre <- NULL; tradeoff_skip <- NULL
  if (!is.null(tradeoff)) {
    .miss_q <- setdiff(c("Empirical.Quality", "PG.MaxLFQ.Quality"), cols)
    if (length(.miss_q)) {
      tradeoff_skip <- "this report has no QuantUMS quality columns (Empirical.Quality and the protein-group quality score), so there is nothing to cut on"
    } else {
      flt <- tradeoff_pre <- dplyr::collect(flt)
    }
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
  # the matrix is built -- so they never reach the between-run step, lmFit or BH.
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

  mm <- maxlfq_matrix(rows, br, quantities, per_peptide)
  E <- mm$E; prot_ids <- mm$prot_ids; br <- mm$br
  filters_applied <- c(filters_applied, mm$filters_applied)
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

  # One matrix per cutoff set: the same filters as above (identification q-values already in
  # tradeoff_pre; keep_runs; the cutoffs; contaminants) and the same maxlfq_matrix().
  tradeoff_out <- NULL
  if (!is.null(tradeoff_pre)) {
    rp <- tradeoff_pre
    if (!is.null(keep_runs) && length(keep_runs)) rp <- rp[rp$Run %in% keep_runs, , drop = FALSE]
    tradeoff_out <- lapply(tradeoff, function(co) {
      r <- rp
      if (co$eq > 0) r <- r[!is.na(r$Empirical.Quality) & r$Empirical.Quality >= co$eq, , drop = FALSE]
      if (co$pgq > 0) r <- r[!is.na(r$PG.MaxLFQ.Quality) & r$PG.MaxLFQ.Quality >= co$pgq, , drop = FALSE]
      if (!is.na(cont$id_column) && drop_contaminants)
        r <- r[!is_contaminant(r[[cont$id_column]], exempt = contaminant_exempt), , drop = FALSE]
      if (!nrow(r)) return(list(label = co$label, pgq = co$pgq, eq = co$eq, E = NULL,
                                why = "no rows survive these cutoffs"))
      m <- tryCatch(suppressWarnings(maxlfq_matrix(r, br, quantities, per_peptide)),
                    error = function(e) conditionMessage(e))
      if (is.character(m)) return(list(label = co$label, pgq = co$pgq, eq = co$eq, E = NULL, why = m))
      list(label = co$label, pgq = co$pgq, eq = co$eq, E = m$E, why = NULL)
    })
    rm(rp)
  }
  list(
    tradeoff = tradeoff_out, tradeoff_skip = tradeoff_skip,
    E = E, genes = genes, n_obs = n_obs, filters_applied = filters_applied,
    contaminants = cont,
    # The q-value columns actually filtered on. Recorded rather than assumed so the
    # emitted reproducibility script names the same columns this run used -- older
    # reports lack PG.Q.Value / Global.*, and hard-coding them would emit a script
    # that filters on a column the run never had.
    q_columns = q_columns,
    # an undeclared report is described as DIA-NN's only when it has DIA-NN's precursor columns
    descriptor = c(maxlfq_descriptor(declared,
                                     if (is.null(declared)) diann_cont_quant_exclude(report_path),
                                     rollup = mm$rollup,
                                     diann_shaped = !is.null(declared) ||
                                       all(c("Precursor.Id", "Precursor.Quantity") %in% cols)),
                   list(quantities = quantities,
                        # the peptide column a median polish ran over (NULL otherwise): the
                        # reproducibility log reads it
                        rollup_feature = mm$feature,
                        # the step this run applied and why (maxlfq_between_run), and the
                        # Methods' words for it, built from that answer
                        between_run = br,
                        normalisation = maxlfq_normalisation_text(br, quantities, declared,
                                                                  dia_nn_norm),
                        dia_nn_normalisation = dia_nn_norm)),
    # what the report declared about itself, verbatim, for the record
    declared_quantity = declared,
    # Sage (2.15.1): which protein x run values rest on an MS2 identification in that run (TRUE)
    # and which only on peptides measured by matching between runs (FALSE); NULL when not
    # declared
    identified = mm$identified
  )
}
