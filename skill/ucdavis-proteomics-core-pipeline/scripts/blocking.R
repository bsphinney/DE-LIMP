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
#
# WHICH contrasts the blocked fit reports (--block-scope, default "within"):
#   within  a contrast BETWEEN blocks that uses at most one sample per block (Old vs Young
#           mice for one bait: 3 vs 3 mice, one IP each) comes from the same quantification
#           fitted with samples as independent; every other contrast -- within, partial, or
#           between but pooling several samples per block (Old vs Young over all baits,
#           which WOULD be pseudo-replicated unblocked) -- from the blocked fit. If any
#           block holds two samples of one group, the independent fit's variance is itself
#           pseudo-replicated, so every contrast comes from the blocked fit.
#   all     every contrast from the blocked fit.
# Why "within" is the default (PROT_0756, 6 mice x 5 IPs, consensus correlation 0.17): a
# between-block contrast that already uses one sample per block (each bait, 3 vs 3 mice)
# compares independent samples, and in a balanced design the independent fit's variance is
# unbiased for every protein. It is not exact: its pooled residuals come from the same mice
# across groups, so they share the mouse effect and their degrees of freedom are overstated
# -- a mild anti-conservatism (blocksim/between.R: type I 0.047-0.071 across per-protein
# correlation 0.05-0.85; a fit on the two groups' samples alone is exact, follow-up P1).
# The blocked fit is worse there. It applies ONE consensus
# correlation to all proteins, so for proteins with strong block-to-block variation it
# understates the between-block variance: blocked/independent SE ratio 0.86 at per-protein
# correlation > 0.6 (1.04 at <= 0), matching the value the design predicts to 0.015 -- and
# the 77 Old-vs-Young calls only the blocked fit made were exactly those proteins (median
# correlation 0.44 vs 0.14 overall: red-cell, complement, tRNA-synthetase proteins). For
# within-block contrasts the blocked fit is the right model: +14-58% calls, none lost.
#
# FIXED or RANDOM (--block-effect, default "auto"; review C2, blocksim/paired.R).
# When the block is CROSSED with the groups (it stays full rank as a fixed term) and every
# contrast compares samples within one block -- before/after in the same patient -- the
# fixed subject effect IS the exact paired analysis: type I 0.050 and power 0.95 at every
# per-protein correlation. The random effect drifts (type I 0.078 -> 0.010 from
# correlation 0.05 to 0.85, power 0.96 -> 0.76) because one consensus correlation is
# applied to every protein. So auto = fixed there; random (duplicateCorrelation) stays for
# nested / multi-level designs such as PROT_0756 (mice within age), where a fixed mouse
# term is aliased with the age groups and between-mouse contrasts need the random effect.
#
# POOLED between-block contrasts (review C3, blocksim/between.R): a contrast between blocks
# that averages several samples per block (Old vs Young over all baits) can only come from
# the blocked fit -- unblocked it is pseudo-replicated (type I 0.07-0.36) -- but the
# blocked fit is itself anti-conservative for proteins with strong block effects (type I
# 0.10 / 0.22 at per-protein correlation 0.6 / 0.85 for a nominal 0.05). Such a contrast is
# reported from the blocked fit WITH a CAUTION; define between-block contrasts one sample
# per block (per bait), which scope "within" reports from the independent fit.
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
BLOCK_FIT_FIXED <- list(
  dpc = "limpa::dpcDE(y, design + fixed block columns) -> voomaLmFitWithImputation (no block argument)",
  maxlfq = "limma::lmFit(E, design + fixed block columns)")

# A block was asked for, so the fit that reports it must actually CARRY it -- never an
# unblocked fit recorded as a blocked one. Random: the within-block correlation lmFit(block =,
# correlation =) used comes back as fit$correlation; a fit without it treated the samples as
# independent (a limpa whose dpcDE accepted block= and dropped it would do exactly that --
# 1.2.5 and 1.4.x pass it through). Fixed: the block columns are in fit$design.
block_fit_check <- function(fit, effect, block_col, block_cols = NULL) {
  if (identical(effect, "random")) {
    rho <- fit$correlation
    if (is.null(rho) || length(rho) != 1L || !is.finite(rho))
      stop(sprintf(paste0(
        "--block %s: the fit came back WITHOUT a within-block correlation -- it treated the ",
        "samples as independent. Refusing to report it as blocked. limpa %s: its dpcDE() must ",
        "pass block= to voomaLmFitWithImputation() (limpa 1.2.5 and 1.4.x do)."),
        block_col, tryCatch(as.character(utils::packageVersion("limpa")), error = function(e) "?")),
        call. = FALSE)
  } else if (identical(effect, "fixed")) {
    have <- colnames(fit$design)
    if (!length(block_cols) || is.null(have) || !all(block_cols %in% have))
      stop(sprintf(paste0(
        "--block %s: the fixed-effect fit came back WITHOUT the %s columns in its design -- it ",
        "treated the samples as independent. Refusing to report it as blocked."),
        block_col, block_col), call. = FALSE)
  }
  invisible(TRUE)
}

# The block as fixed design columns (treatment coding; the first level is absorbed by the
# group means). One definition: run_de.R fits them, repro_script.R emits the same code.
block_fixed_columns <- function(block, block_col) {
  f <- factor(block)
  z <- stats::model.matrix(~ f)[, -1, drop = FALSE]
  colnames(z) <- make.names(paste0(block_col, "_", levels(f)[-1]))
  z
}

# Fixed or random (header, review C2). auto: fixed when the block is crossed with the groups
# (full rank beside the group columns) AND every contrast is within one block; random
# otherwise. A covariate the block is nested in (each patient's pre/post pair run in one
# Batch) is ABSORBED by a fixed block -- its columns are sums of the block's -- so the
# fixed fit drops it, and says so (`absorbed`). What the block is aliased with is named:
# the groups, or covariate X. A requested "fixed" that the groups alias stops.
# `covariates` are the design's covariate terms, in model.matrix order (attr "assign").
block_choose_effect <- function(requested, block, design, structure, block_col,
                                covariates = character(0)) {
  bcols <- block_fixed_columns(block, block_col)
  full_rank <- function(X) qr(X)$rank == ncol(X)
  asg <- attr(design, "assign")
  if (is.null(asg)) asg <- rep(1L, ncol(design))           # no term map: all group columns
  term_cols <- function(k) design[, asg == k, drop = FALSE]  # 1 = groups, 1 + i = covariates[i]
  groups_ok <- full_rank(cbind(term_cols(1L), bcols))
  # Absorbed = ALL of the covariate's columns lie in the span of groups + block (adding them
  # adds no rank). Rank-deficient is not enough: a 3-level Batch can be PARTLY aliased --
  # nested for some subjects, confounded with treatment within others (S1-2 both runs in
  # b1; S3-6 Ctrl in b2, Trt in b3) -- and dropping it would move logFC by up to 4.1
  # (blocksim/absorb.R). Such a covariate is kept; the fixed design is then rank-deficient,
  # so auto goes random and a forced fixed stops.
  rank_gb <- qr(cbind(term_cols(1L), bcols))$rank
  nests_in <- if (groups_ok) Filter(function(cv)
    qr(cbind(term_cols(1L), bcols, term_cols(1L + match(cv, covariates))))$rank == rank_gb,
    covariates) else character(0)
  keep <- setdiff(covariates, nests_in)
  dfix <- cbind(design[, asg %in% c(1L, 1L + match(keep, covariates)), drop = FALSE], bcols)
  crossed <- groups_ok && full_rank(dfix)
  aliased_with <- if (!groups_ok) "the groups"
                  else if (!crossed) "the groups and covariates together"
                  else NULL
  st <- unlist(structure)
  all_within <- length(st) > 0 && all(st == "within")
  auto <- if (crossed && all_within) "fixed" else "random"
  effect <- if (identical(requested, "auto")) auto else requested
  if (identical(effect, "fixed") && !crossed)
    stop(sprintf(paste0("--block-effect fixed: %s is nested in %s (as fixed columns it is aliased ",
                        "with them%s), so it cannot be a fixed effect here. Use --block-effect ",
                        "random (the default for this design)."),
                 block_col, aliased_with, if (!groups_ok) " -- mice within age" else ""),
         call. = FALSE)
  nest_txt <- if (length(nests_in))
    sprintf("%s is nested in covariate %s (each %s sits in one %s)", block_col,
            paste(nests_in, collapse = " / "), block_col, paste(nests_in, collapse = " / ")) else NULL
  why <- if (crossed && all_within)
    paste0(sprintf("%s is crossed with the groups and every contrast is within one %s: the fixed effect is the exact paired analysis",
                   block_col, block_col),
           if (length(nests_in)) sprintf("; %s, so its fixed effect absorbs %s, which is dropped from the design",
                                         nest_txt, paste(nests_in, collapse = " / ")) else "")
  else if (!crossed)
    sprintf("%s is nested in %s (aliased as a fixed term): random effect", block_col, aliased_with)
  else paste0(sprintf("not every contrast is within one %s: random effect", block_col),
              if (length(nests_in)) sprintf(" (%s)", nest_txt) else "")
  fixed <- identical(effect, "fixed")
  list(effect = effect,
       choice = sprintf("%s -- auto would be %s: %s", if (identical(requested, "auto")) "auto"
                        else paste("requested", requested), auto, why),
       absorbed = if (fixed) nests_in else character(0),
       design = if (fixed) dfix else design)
}

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
# technical = TRUE: the block is the biological sample of technical replicates (run_de.R's
# Sample column), where a sample injected once is a block of one -- duplicateCorrelation takes
# it -- so only "no sample has two runs left" stops.
block_check <- function(block, design, block_col, technical = FALSE) {
  sizes <- table(block)
  if (length(sizes) < 2)
    stop(sprintf("--block %s: every sample is in the same block ('%s'); there is nothing to estimate.",
                 block_col, names(sizes)[1]), call. = FALSE)
  if (technical && !any(sizes >= 2))
    stop(sprintf(paste0(
      "%s: no sample has two analysed runs left, so there are no technical replicates to block ",
      "on (a run removed from the analysis can leave its sample with one). Drop the %s column, ",
      "or check which runs were removed."), block_col, block_col), call. = FALSE)
  small <- if (technical) character(0) else names(sizes)[sizes < 2]
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
                         first_pass = NULL, groups = NULL, cmat = NULL, scope = "within",
                         effect = "random", effect_choice = NULL, absorbed = character(0)) {
  structure <- if (!is.null(cmat)) block_contrast_structure(block, groups, cmat) else NULL
  model <- if (!is.null(cmat)) block_contrast_model(block, groups, cmat, scope) else NULL
  fixed <- identical(effect, "fixed")
  if (fixed && length(model)) model <- lapply(model, function(x) "blocked")
  arho <- atanh_per_protein[is.finite(atanh_per_protein)]
  q <- if (length(arho)) tanh(stats::quantile(arho, c(0.25, 0.5, 0.75), names = FALSE))
       else rep(NA_real_, 3)
  nb <- length(unique(block))
  warn <- character(0)
  if (!fixed) {
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
  # A between-block contrast on the blocked fit: its single consensus correlation understates
  # the between-block variance of proteins with strong block effects (review C3).
  bb <- names(structure)[unlist(structure) == "between" & unlist(model)[names(structure)] == "blocked"]
  if (length(bb))
    warn <- c(warn, sprintf(paste0(
      "%s compare%s different %s levels but %s reported from the blocked fit (%s). Its single ",
      "consensus correlation understates the between-%s variance of proteins with strong ",
      "%s-to-%s variation, so these p-values are anti-conservative for exactly those proteins ",
      "(simulated type I error 0.10 / 0.22 at per-protein correlation 0.6 / 0.85, nominal 0.05). ",
      "Prefer between-%s contrasts that use one sample per %s (e.g. per bait), which ",
      "--block-scope within reports from the fit with samples independent."),
      paste(bb, collapse = ", "), if (length(bb) == 1) "s" else "", block_col,
      if (length(bb) == 1) "is" else "are",
      if (identical(scope, "all")) "--block-scope all"
      else "the fit with samples independent would pseudo-replicate it: several samples per level",
      block_col, block_col, block_col, block_col, block_col))
  }
  # applied = the blocked fit reports at least one contrast -- what readers of the record
  # ask ("does the design carry the block?"). All-between under "within" = it reports none.
  used <- !length(model) || any(unlist(model) == "blocked")
  list(
    applied = used, column = block_col, scope = scope,
    # "fixed": a block term in the design (the exact paired analysis); "random":
    # duplicateCorrelation. effect_choice says why (block_choose_effect).
    effect = effect, effect_choice = effect_choice,
    # covariates a fixed block absorbed and the fit therefore dropped (block_choose_effect)
    absorbed_covariates = as.list(absorbed),
    note = if (!used) sprintf(paste0("--block %s given, but every contrast compares different %s ",
                                     "levels: all were reported from the independent fit ",
                                     "(--block-scope within)"), block_col, block_col) else NULL,
    n_blocks = nb,
    block_sizes = as.list(table(block)),
    consensus_correlation = if (fixed) NULL else consensus,
    first_pass_correlation = if (fixed) NULL else first_pass,
    n_proteins = n_proteins, n_proteins_estimated = if (fixed) NULL else length(arho),
    per_protein_correlation = if (fixed) NULL else list(q25 = q[1], median = q[2], q75 = q[3]),
    estimator = if (fixed) sprintf("none: %s is a fixed effect (one coefficient per level)", block_col)
                else BLOCK_ESTIMATOR,
    fit = if (fixed) BLOCK_FIT_FIXED[[method]] else BLOCK_FIT[[method]],
    contrast_structure = structure,
    contrast_model = model,
    contrast_model_rule = if (fixed) sprintf("fixed %s effect: every contrast from that fit", block_col)
      else if (identical(scope, "all")) "all: every contrast from the blocked fit"
      else paste0("within: a contrast comparing different ", block_col, " levels with at most one ",
                  "sample per ", block_col, " is reported from the fit with samples independent; ",
                  "every other contrast from the blocked fit",
                  if (block_reps_within_group(block, groups))
                    sprintf(paste0(" (here every contrast: some %s holds two samples of one group, ",
                                   "which the independent fit would pseudo-replicate)"), block_col)
                  else ""),
    warnings = as.list(warn))
}

# Does any block hold two samples of the same group (technical replicates of one mouse)?
# Then even the independent fit's pooled variance is pseudo-replicated.
block_reps_within_group <- function(block, groups)
  !is.null(groups) && anyDuplicated(paste(block, as.character(groups), sep = "\r")) > 0

# The fit that reports each contrast: "blocked" or "independent" (see the header). Under
# "within", independent only for a between-block contrast that takes at most one sample
# from each block, in a design with no block holding two samples of one group.
block_contrast_model <- function(block, groups, cmat, scope) {
  st <- block_contrast_structure(block, groups, cmat)
  g <- as.character(groups)
  reps <- block_reps_within_group(block, groups)
  out <- list()
  for (cn in colnames(cmat)) {
    used <- g %in% rownames(cmat)[cmat[, cn] != 0]
    one_each <- !anyDuplicated(block[used])
    out[[cn]] <- if (identical(scope, "within") && identical(st[[cn]], "between") &&
                     one_each && !reps) "independent" else "blocked"
  }
  out
}

# One fit object whose contrast columns come from each contrast's reporting fit, so
# topTable() on it returns exactly that fit's table (it reads coefficients,
# stdev.unscaled, t, p.value, lods and Amean -- identical Amean in both fits). The
# per-protein variance fields (sigma, s2.post, df.total) stay the blocked fit's, so the
# moderated F -- which would mix the two -- is dropped; the session keeps the independent
# fit whole beside it.
block_merge_fits <- function(fit_blocked, fit_independent, model) {
  # every contrast independent: nothing of the blocked fit (s2.post, correlation) may ride along
  if (!any(model == "blocked")) { fit_independent$contrast_model <- model; return(fit_independent) }
  ind <- names(model)[model == "independent"]
  if (!length(ind)) { fit_blocked$contrast_model <- model; return(fit_blocked) }
  stopifnot(identical(rownames(fit_blocked$coefficients), rownames(fit_independent$coefficients)))
  out <- fit_blocked
  for (el in c("coefficients", "stdev.unscaled", "t", "p.value", "lods"))
    out[[el]][, ind] <- fit_independent[[el]][, ind]
  out$F <- NULL; out$F.p.value <- NULL
  out$contrast_model <- model
  out
}

block_none_record <- function()
  list(applied = FALSE, note = "no --block: samples modelled as independent")

# How the contrasts split between the fits, in one phrase for labels and methods.
block_scope_phrase <- function(rec) {
  if (identical(rec$effect, "fixed"))
    return(sprintf("all contrasts from the fit with %s as a fixed effect", rec$column))
  m <- unlist(rec$contrast_model)
  if (!any(m == "independent")) return(sprintf("all contrasts from the blocked fit (scope %s)", rec$scope))
  sprintf(paste0("between-%s contrasts using at most one sample per %s from the fit with ",
                 "samples independent, all others from the blocked fit (scope %s)"),
          rec$column, rec$column, rec$scope)
}

# Appended to the descriptor's de_engine so every reader of that label (methods.txt,
# de_provenance.json, the AI brief) sees the blocking too.
block_engine_suffix <- function(rec) {
  if (!isTRUE(rec$applied)) return("")
  if (identical(rec$effect, "fixed"))
    return(sprintf("; block = %s (fixed effect: crossed with the groups, every contrast within one %s)",
                   rec$column, rec$column))
  sprintf("; block = %s (random effect, consensus correlation %.3f; %s)",
          rec$column, rec$consensus_correlation, block_scope_phrase(rec))
}

block_methods_lines <- function(rec) {
  pad <- "                "
  tr <- rec$technical_replicates
  tech <- if (!is.null(tr))
    sprintf(paste0("Replicates    : %d sample(s) were injected more than once (%d runs; %s column): ",
                   "technical replicates, blocked on %s below, never counted as independent samples"),
            tr$n_samples, tr$n_runs, tr$column, tr$column)
  if (!isTRUE(rec$applied) && is.null(rec$column))
    return("Blocking      : none -- samples modelled as independent (no --block)")
  if (!isTRUE(rec$applied))
    return(sprintf("Blocking      : none used -- %s", rec$note))
  sz <- unlist(rec$block_sizes)
  sz_txt <- if (length(unique(sz)) == 1) sprintf("%d samples each", sz[1])
            else sprintf("%d-%d samples each", min(sz), max(sz))
  st <- unlist(rec$contrast_structure)
  pc <- rec$per_protein_correlation
  fixed <- identical(rec$effect, "fixed")
  out <- c(sprintf("Blocking      : %s as a %s effect (%d levels, %s)", rec$column,
                   if (fixed) "FIXED" else "random", rec$n_blocks, sz_txt),
    if (!is.null(rec$effect_choice)) sprintf("%sEffect: %s", pad, rec$effect_choice),
    if (length(rec$absorbed_covariates))
      sprintf("%sDropped from the design: %s (absorbed by the fixed %s effect)", pad,
              paste(unlist(rec$absorbed_covariates), collapse = ", "), rec$column),
    if (!fixed) sprintf("%sconsensus within-%s correlation %.3f, from %d of %d proteins", pad, rec$column,
            rec$consensus_correlation, rec$n_proteins_estimated, rec$n_proteins),
    if (!fixed && is.finite(pc$median))
      sprintf("%s(per-protein median %.2f, IQR %.2f to %.2f)", pad, pc$median, pc$q25, pc$q75),
    if (!is.null(rec$first_pass_correlation) && is.finite(rec$first_pass_correlation))
      sprintf("%sfirst-pass estimate %.3f", pad, rec$first_pass_correlation),
    sprintf("%sEstimator: %s", pad, rec$estimator),
    sprintf("%sFit: %s", pad, rec$fit),
    sprintf("%sScope: %s", pad, block_scope_phrase(rec)),
    if (length(st)) {
      md <- unlist(rec$contrast_model)[names(st)]
      unlist(lapply(c("within", "partial", "between", "other"), function(k)
        vapply(c("blocked", "independent"), function(m) {
          hit <- names(st)[st == k & md == m]
          if (!length(hit)) return(NA_character_)
          sprintf("%s%s contrasts, %s fit%s: %s", pad,
                  if (k == "other") "Other" else paste0(tools::toTitleCase(k), "-", rec$column), m,
                  if (k == "between") sprintf(" (judged on the %s levels, not the samples)", rec$column) else "",
                  paste(hit, collapse = ", "))
        }, character(1))))
    } else NULL,
    if (length(rec$warnings)) sprintf("%sCAUTION: %s", pad, unlist(rec$warnings)))
  c(tech, out[!is.na(out)])
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

# ...and the NESTED case block_candidates() cannot see: a column that splits a group into
# units holding several runs each -- technical replicates, several injections or fractions
# of one animal (Mouse = Ctrl_M1 x3, Ctrl_M2 x3 inside Ctrl). Counted as independent they
# inflate n. A column constant within every group is a coarser label, not a unit, and is
# not flagged. Used for a hint only.
block_candidates_within <- function(meta, covariates) {
  hits <- character(0)
  for (cc in setdiff(names(meta), c("File.Name", "Group", covariates))) {
    v <- as.character(meta[[cc]])
    if (anyNA(v) || !all(nzchar(v))) next
    splits <- vapply(split(v, meta$Group), function(x)
      length(unique(x)) >= 2 && anyDuplicated(x) > 0, logical(1))
    if (any(splits)) hits <- c(hits, cc)
  }
  hits
}
