# =============================================================================
# blocking.R -- a blocking factor for run_de.R (--block <column>).
#
# Samples that share a source -- the IPs cut from one mouse brain, the biopsies from one
# patient -- are correlated. Modelling them as independent throws that pairing away.
# Why (PROT_0756, Dickson lab, 2026-09): 30 IPs from 6 mouse brains, one JPH3/JPH4/
# Kv2.1/RyR/IgG IP per mouse. Every bait-vs-IgG contrast compares IPs from the SAME
# mice, yet run_de.R fitted all 30 as independent samples.
#
# The model is limma's for multi-level experiments: the block is a random effect, one
# consensus within-block correlation is estimated across all proteins by
# duplicateCorrelation(), and lmFit(block =, correlation =) fits generalised least
# squares with it. The block may be NESTED in a fixed factor (mouse within age) -- the
# design this exists for; contrasts between mice (Old vs Young) and within mice (bait vs
# IgG) are both handled. It must not ALSO be a fixed effect: its correlation is then not
# estimable, and duplicateCorrelation() returns zero with only a warning -- an unblocked
# analysis recorded as a blocked one. block_check() stops on that instead.
#
# How each --method fits it (the record's `fit` says which ran):
#   dpc     limpa::dpcDE(y, design, block = b). dpcDE(y, design, plot, ...) passes `...`
#           to voomaLmFitWithImputation(), whose `block` argument (limpa 1.4.0 source):
#           fits, builds the vooma precision weights, estimates the correlation with
#           duplicateCorrelation(y, design, block, weights) ("First intra-block
#           correlation"), refits lmFit(block, correlation), recomputes the weights,
#           re-estimates it ("Final intra-block correlation"), and fits the final
#           lmFit(y, design, block, correlation, weights). Nothing to re-implement.
#   maxlfq  duplicateCorrelation(E, design, block = b) -> lmFit(E, design, block =,
#           correlation =). One pass: this path fits without precision weights.
# =============================================================================

# Warning thresholds -- they only decide whether the run says CAUTION; the fit never
# depends on them. The consensus is a trimmed mean over proteins, so it needs many.
BLOCK_MIN_PROTEINS   <- 50
# limpa estimates the correlation twice (before and after the weights are recomputed);
# a large move between the two means the estimate is not settled.
BLOCK_MAX_PASS_SHIFT <- 0.1
BLOCK_ESTIMATOR <- paste0("limma::duplicateCorrelation -- per-protein correlations from a ",
                          "two-level mixed model, consensus = tanh of their 15%-trimmed mean ",
                          "on the atanh scale")
BLOCK_FIT <- list(
  dpc = paste0("limpa::dpcDE(block =) -> voomaLmFitWithImputation: duplicateCorrelation with ",
               "the vooma precision weights, re-estimated after the weights are recomputed, ",
               "then lmFit(block =, correlation =, weights =)"),
  maxlfq = "limma::duplicateCorrelation(E, design, block =) -> lmFit(E, design, block =, correlation =)")

# The checks that need only conditions.csv -- run before quantification, which can take
# an hour, so a wrong column fails in seconds.
block_validate_column <- function(meta, block_col, covariates) {
  if (!block_col %in% names(meta))
    stop(sprintf("--block %s: no such column in the metadata (columns: %s).",
                 block_col, paste(names(meta), collapse = ", ")), call. = FALSE)
  if (block_col %in% covariates)
    stop(sprintf(paste0(
      "--block %s: '%s' is also a fixed-effect covariate in the design (~ 0 + groups + %s; ",
      "run_de.R fits any column named %s as one). A column can be a fixed effect or a random ",
      "blocking factor, not both -- as a fixed effect its within-block correlation cannot be ",
      "estimated. Keep it once: rename the column (e.g. 'Mouse') and pass --block with the ",
      "new name, or drop --block to keep it fixed."),
      block_col, block_col, paste(covariates, collapse = " + "),
      paste(covariates, collapse = " / ")), call. = FALSE)
  if (block_col %in% c("File.Name", "Group"))
    stop(sprintf("--block %s: that column is the %s; name the column holding the blocking unit (e.g. Mouse).",
                 block_col, if (block_col == "Group") "fixed grouping itself" else "sample identifier"),
         call. = FALSE)
  b <- trimws(as.character(meta[[block_col]]))
  if (any(is.na(b) | !nzchar(b)))
    stop(sprintf("--block %s: %d sample(s) have no value in that column, e.g. %s. Every sample needs one.",
                 block_col, sum(is.na(b) | !nzchar(b)),
                 meta$File.Name[which(is.na(b) | !nzchar(b))[1]]), call. = FALSE)
  invisible(b)
}

# The checks that need the design as it will be fitted. Returns the block sizes.
block_check <- function(block, design, block_col) {
  sizes <- table(block)
  if (length(sizes) < 2)
    stop(sprintf("--block %s: every sample is in the same block ('%s'); there is nothing to estimate.",
                 block_col, names(sizes)[1]), call. = FALSE)
  small <- names(sizes)[sizes < 2]
  if (length(small))
    stop(sprintf(paste0(
      "--block %s: %d block(s) hold a single analysed sample (%s). A block needs at least 2 ",
      "samples for its within-block correlation to mean anything -- check the %s column ",
      "(a sample removed from the analysis can leave its block with one)."),
      block_col, length(small), paste(utils::head(small, 5), collapse = ", "), block_col),
      call. = FALSE)
  # limma's own test (duplicateCorrelation): the block is "already encoded in the design"
  # when its indicators lie in the design's column space -- then limma sets the correlation
  # to zero and carries on. Stop instead.
  qrd <- qr(design)
  z <- stats::model.matrix(~ factor(block))[, -1, drop = FALSE]
  r <- qr.qty(qrd, z)
  if (qrd$rank < nrow(design) && max(abs(r[-seq_len(qrd$rank), , drop = FALSE])) < 1e-8)
    stop(sprintf(paste0(
      "--block %s is already a fixed effect in the design (%s): its levels coincide with the ",
      "groups or a covariate, so a within-block correlation cannot be estimated. Use it as ",
      "one or the other, not both."), block_col, paste(colnames(design), collapse = ", ")),
      call. = FALSE)
  invisible(sizes)
}

# Which blocks each contrast compares: "within" (the same blocks on both sides -- paired,
# where blocking gains power), "between" (disjoint blocks -- e.g. Old vs Young mice, where
# the evidence is the number of blocks, not samples), "partial", or "other" (a contrast
# that does not compare groups).
block_contrast_structure <- function(block, groups, cmat) {
  g <- as.character(groups)
  out <- list()
  for (cn in colnames(cmat)) {
    w <- cmat[, cn]
    bp <- unique(block[g %in% names(w)[w > 0]]); bn <- unique(block[g %in% names(w)[w < 0]])
    out[[cn]] <- if (!length(bp) || !length(bn)) "other"
                 else if (setequal(bp, bn)) "within"
                 else if (!length(intersect(bp, bn))) "between"
                 else "partial"
  }
  out
}

# The ONE record of the blocking step: de_provenance.json's `block`, methods.txt's
# Blocking line, the de_engine label and make_methods.py's sentence all read it.
block_record <- function(block_col, block, method, consensus, atanh_per_protein, n_proteins,
                         first_pass = NULL, structure = NULL) {
  arho <- atanh_per_protein[is.finite(atanh_per_protein)]
  q <- if (length(arho)) tanh(stats::quantile(arho, c(0.25, 0.5, 0.75), names = FALSE))
       else rep(NA_real_, 3)
  nb <- length(unique(block))
  warn <- character(0)
  if (!is.finite(consensus) || consensus <= 0)
    warn <- c(warn, sprintf(paste0(
      "the consensus within-%s correlation is %s (<= 0): samples sharing a %s are no more ",
      "alike than samples from different ones, so blocking gains nothing here. Check the %s ",
      "assignments; if they are right, the unblocked analysis is the simpler equivalent."),
      block_col, format(round(consensus, 3)), block_col, block_col))
  if (length(arho) < BLOCK_MIN_PROTEINS)
    warn <- c(warn, sprintf(paste0(
      "the within-%s correlation is unstable: only %d protein(s) gave an estimate (fewer than ",
      "%d), and the consensus is an average over proteins."),
      block_col, length(arho), BLOCK_MIN_PROTEINS))
  if (!is.null(first_pass) && is.finite(first_pass) && is.finite(consensus) &&
      abs(consensus - first_pass) > BLOCK_MAX_PASS_SHIFT)
    warn <- c(warn, sprintf(paste0(
      "the within-%s correlation is unstable: it moved from %.3f to %.3f between limpa's two ",
      "estimation passes (more than %.2f)."),
      block_col, first_pass, consensus, BLOCK_MAX_PASS_SHIFT))
  if (nb < 3)
    warn <- c(warn, sprintf(paste0(
      "only %d %s levels: each protein's between-%s variance rests on %d degree(s) of freedom, ",
      "so the per-protein estimates are very noisy."), nb, block_col, block_col, nb - 1))
  list(
    applied = TRUE, column = block_col, n_blocks = nb,
    block_sizes = as.list(table(block)),
    consensus_correlation = consensus,
    first_pass_correlation = first_pass,
    n_proteins = n_proteins, n_proteins_estimated = length(arho),
    per_protein_correlation = list(q25 = q[1], median = q[2], q75 = q[3]),
    estimator = BLOCK_ESTIMATOR, fit = BLOCK_FIT[[method]],
    contrast_structure = structure,
    warnings = as.list(warn))
}

block_none_record <- function()
  list(applied = FALSE, note = "no --block: samples modelled as independent")

# Appended to the descriptor's de_engine so every reader of that label (methods.txt,
# de_provenance.json, the AI brief) sees the blocking too.
block_engine_suffix <- function(rec) {
  if (!isTRUE(rec$applied)) return("")
  sprintf("; block = %s (random effect, consensus correlation %.3f)",
          rec$column, rec$consensus_correlation)
}

block_methods_lines <- function(rec) {
  pad <- "                "
  if (!isTRUE(rec$applied))
    return("Blocking      : none -- samples modelled as independent (no --block)")
  sz <- unlist(rec$block_sizes)
  sz_txt <- if (length(unique(sz)) == 1) sprintf("%d samples each", sz[1])
            else sprintf("%d-%d samples each", min(sz), max(sz))
  st <- unlist(rec$contrast_structure)
  pc <- rec$per_protein_correlation
  out <- c(sprintf("Blocking      : %s as a random effect (%d levels, %s)", rec$column, rec$n_blocks, sz_txt),
    sprintf("%sconsensus within-%s correlation %.3f, from %d of %d proteins", pad, rec$column,
            rec$consensus_correlation, rec$n_proteins_estimated, rec$n_proteins),
    if (is.finite(pc$median))
      sprintf("%s(per-protein median %.2f, IQR %.2f to %.2f)", pad, pc$median, pc$q25, pc$q75),
    if (!is.null(rec$first_pass_correlation) && is.finite(rec$first_pass_correlation))
      sprintf("%sfirst-pass estimate %.3f", pad, rec$first_pass_correlation),
    sprintf("%sEstimator: %s", pad, rec$estimator),
    sprintf("%sFit: %s", pad, rec$fit),
    if (length(st)) vapply(c("within", "between", "partial"), function(k)
      if (any(st == k)) sprintf("%s%s-%s contrasts%s: %s", pad, tools::toTitleCase(k), rec$column,
                                if (k == "between") sprintf(" (judged on the %s levels, not the samples)",
                                                            rec$column) else "",
                                paste(names(st)[st == k], collapse = ", ")) else NA_character_,
      character(1)) else NULL,
    if (length(rec$warnings)) sprintf("%sCAUTION: %s", pad, unlist(rec$warnings)))
  out[!is.na(out)]
}

# Columns that look like a blocking unit the user did not pass: a level that recurs
# across groups (the same mouse in several conditions). Used for a hint only.
block_candidates <- function(meta, covariates) {
  hits <- character(0)
  for (cc in setdiff(names(meta), c("File.Name", "Group", covariates))) {
    v <- as.character(meta[[cc]])
    if (anyNA(v) || !all(nzchar(v))) next
    tab <- table(v)
    if (length(tab) < 2 || length(tab) == length(v) || any(tab < 2)) next
    if (any(tapply(meta$Group, v, function(g) length(unique(g))) > 1)) hits <- c(hits, cc)
  }
  hits
}
