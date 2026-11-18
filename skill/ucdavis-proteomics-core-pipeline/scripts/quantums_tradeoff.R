# quantums_tradeoff.R -- the QuantUMS cutoff trade-off table (MaxLFQ route only).
#
# What a QuantUMS cutoff costs and buys, shown, never applied silently: the DE is re-fitted at a
# fixed set of cutoffs and, per contrast, the table reports proteins tested, significant calls,
# calls on incomplete data (incomplete_data.R) and proteins lost relative to no filter. The cutoffs
# actually used are marked. The choice stays the user's: the skill's defaults are NOT changed by
# this table (pgQ / eQ default to off), and a recommendation is only a label.
#
# Source: filter-tuning prototype, Phinney lab, 2026-11-13 (FILTER_TUNING_PROTOTYPE.md), PXD063416
# known-answer benchmark: pgQ 0.75 + eQ 0.75 was close to the best cutoff in every DIA window and
# clearly better than none; no rule that works without the answer beat it; above pgQ 0.85 the
# true calls fall apart (pgQ 0.9 lost 27-57% of them).
#
# ONE DEFINITION (rule 3): the five cutoffs and the recommended row are constants HERE. run_de.R
# builds each cutoff's matrix with build_maxlfq.R maxlfq_matrix() (the run's own function), fits it
# with refit_fn (the run's own design, contrasts, block, coverage filter), and calls it with the
# run's own rule (adj.P.Val < the threshold), through the same incomplete_flag_table().

QUANTUMS_TRADEOFF_CUTOFFS <- list(
  list(label = "none",               pgq = 0,    eq = 0),
  list(label = "pgQ 0.5",            pgq = 0.5,  eq = 0),
  list(label = "pgQ 0.75",           pgq = 0.75, eq = 0),
  list(label = "pgQ 0.75 + eQ 0.75", pgq = 0.75, eq = 0.75),
  list(label = "pgQ 0.85",           pgq = 0.85, eq = 0))

QUANTUMS_RECOMMENDED <- list(
  pgq = 0.75, eq = 0.75, tag = "benchmark-recommended (PXD063416)",
  citation = paste0("Phinney lab filter-tuning prototype, 2026-11-13 (FILTER_TUNING_PROTOTYPE.md): ",
                    "on the PXD063416 known-answer benchmark pgQ 0.75 + eQ 0.75 was within about 1-3 wrong ",
                    "calls of the best cutoff in every DIA window and clearly better than no filter"))

QUANTUMS_DEFAULT_TAG <- "DEFAULT — not user-confirmed"

QUANTUMS_COLUMNS <- c("Cutoff", "pgQ_min", "eQ_min", "Contrast", "proteins_tested", "significant",
                      "called_on_incomplete_data", "called_single_sample", "proteins_lost_vs_none",
                      "used_in_this_run", "cutoff_source", "note")

.q_same <- function(a, b) isTRUE(all.equal(as.numeric(a), as.numeric(b)))

# The cutoff sets to build: the fixed five, plus the run's own when it is not one of them.
quantums_cutoffs_for_run <- function(pgq_used, eq_used) {
  cs <- QUANTUMS_TRADEOFF_CUTOFFS
  if (!any(vapply(cs, function(c) .q_same(c$pgq, pgq_used) && .q_same(c$eq, eq_used), TRUE)))
    cs[[length(cs) + 1L]] <- list(label = sprintf("as run (pgQ %g, eQ %g)", pgq_used, eq_used),
                                  pgq = pgq_used, eq = eq_used)
  cs
}

# Where the used cutoffs came from (rule 2): a cutoff the user did not set is DEFAULT, never silent.
quantums_source_text <- function(pgq_set, eq_set) {
  if (!pgq_set && !eq_set) return(QUANTUMS_DEFAULT_TAG)
  sprintf("pgQ: %s; eQ: %s", if (pgq_set) "user-chosen" else QUANTUMS_DEFAULT_TAG,
          if (eq_set) "user-chosen" else QUANTUMS_DEFAULT_TAG)
}

# mats: build_maxlfq()$tradeoff (list of list(label, pgq, eq, E, why)); run_names: the DE's runs
# (columns, in order); keep_fn(Ec) -> logical per row: the DE's own coverage rule (build_maxlfq.R
# coverage_rule, 2.15.1: all samples or any one group); cmat: the contrast
# matrix; groups: factor of run groups; adjp: the call threshold; used: c(pgq, eq) of the run;
# source: quantums_source_text(); refit_fn(Ec) -> list(<contrast> = data.frame(P.Value) with the
# proteins as row names), the run's own fit; untested_fn is not needed: a protein measured in no
# run of the groups a contrast compares is untested here exactly as in the DE tables.
# -> data.frame with QUANTUMS_COLUMNS.
quantums_tradeoff_table <- function(mats, run_names, keep_fn, cmat, groups, adjp, used, source,
                                    refit_fn) {
  cn_all <- colnames(cmat)
  cg_of <- function(cn) intersect(c(rownames(cmat)[cmat[, cn] > 0], rownames(cmat)[cmat[, cn] < 0]),
                                  levels(groups))
  tested_none <- list()
  rows <- list()
  note_for <- function(m) {
    n <- character(0)
    if (.q_same(m$pgq, QUANTUMS_RECOMMENDED$pgq) && .q_same(m$eq, QUANTUMS_RECOMMENDED$eq))
      n <- c(n, QUANTUMS_RECOMMENDED$tag)
    paste(n, collapse = "; ")
  }
  for (m in mats) {
    is_used <- .q_same(m$pgq, used[1]) && .q_same(m$eq, used[2])
    base <- list(Cutoff = m$label, pgQ_min = m$pgq, eQ_min = m$eq,
                 used_in_this_run = is_used, cutoff_source = if (is_used) source else "")
    why <- m$why
    Ec <- m$E
    if (is.null(why) && !all(run_names %in% colnames(Ec)))
      why <- sprintf("run(s) %s have too few quantified proteins at these cutoffs and were dropped, so the design cannot be fitted",
                     paste(setdiff(run_names, colnames(Ec)), collapse = ", "))
    if (is.null(why)) {
      Ec <- Ec[, run_names, drop = FALSE]
      Ec <- Ec[keep_fn(Ec), , drop = FALSE]
      if (nrow(Ec) < 10) why <- sprintf("only %d proteins pass the coverage filter", nrow(Ec))
    }
    res <- NULL
    if (is.null(why)) res <- tryCatch(refit_fn(Ec), error = function(e) { why <<- paste("fit failed:", conditionMessage(e)); NULL })
    for (cn in cn_all) {
      row <- base; row$Contrast <- cn
      if (is.null(res)) {
        row <- c(row, list(proteins_tested = NA_integer_, significant = NA_integer_,
                           called_on_incomplete_data = NA_integer_, called_single_sample = NA_integer_,
                           proteins_lost_vs_none = NA_integer_,
                           note = paste(c(note_for(m), paste("not fitted:", why)), collapse = "; ")))
      } else {
        cg <- cg_of(cn)
        k <- vapply(cg, function(g) rowSums(!is.na(Ec[, as.character(groups) == g, drop = FALSE])),
                    numeric(nrow(Ec)))
        k <- matrix(k, nrow = nrow(Ec))
        n <- vapply(cg, function(g) sum(as.character(groups) == g), 0L)
        tt <- res[[cn]]
        p <- tt$P.Value[match(rownames(Ec), rownames(tt))]
        p[rowSums(k) == 0] <- NA_real_                    # measured in no compared run: not tested
        adj <- stats::p.adjust(p, method = "BH")
        sig <- !is.na(adj) & adj < adjp
        fl <- incomplete_flag_table(k, n, cg, sig)
        tested <- rownames(Ec)[!is.na(p)]
        if (identical(m$label, QUANTUMS_TRADEOFF_CUTOFFS[[1]]$label)) tested_none[[cn]] <- tested
        row <- c(row, list(proteins_tested = length(tested), significant = sum(sig),
                           called_on_incomplete_data = sum(fl$called_on_incomplete_data %in% TRUE),
                           called_single_sample = sum(fl$incomplete_tier %in% INCOMPLETE_TIERS[["single"]]),
                           proteins_lost_vs_none = if (is.null(tested_none[[cn]])) NA_integer_
                                                   else sum(!tested_none[[cn]] %in% tested),
                           note = note_for(m)))
      }
      rows[[length(rows) + 1L]] <- as.data.frame(row[QUANTUMS_COLUMNS], stringsAsFactors = FALSE,
                                                 check.names = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
