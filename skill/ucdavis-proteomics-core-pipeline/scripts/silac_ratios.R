# =============================================================================
# silac_ratios.R -- the SILAC ratio statistics, as plain functions (data.table + limma).
#
# run_silac.R sources this file for its numbers, and copies it VERBATIM into the analysis's
# reproducibility_log.R, so the replay runs the code that produced the tables -- not a second
# description of it (DE-LIMP rule 3: one definition).
#
# The method, as validated on a triple-SILAC experiment (DIA-NN 2.7.0 plexDIA, 4 replicates,
# 2026-10-14; references/silac.md) and reviewed by a statistician (2026-10-15):
#   * precursors passing DIA-NN's q-value columns (`cutoffs`: the run-specific Q.Value and the
#     channel-specific PG.Q.Value, and the experiment-wide Global.* / Lib.* columns whenever the
#     report has them -- diann_q_columns.R, the label-free path's set), contaminants removed;
#   * a precursor ratio needs both channels at Channel.Q.Value <= channel_q, with a non-zero raw
#     quantity (Precursor.Quantity: DIA-NN's channel normalisation is not used);
#   * per protein group and RUN, the median of its precursor log2 ratios, median-centred within
#     the run (the zero is the median detected protein; the uncentred ratios are kept: they carry
#     the mixing offset);
#   * per REPLICATE, the mean of its runs' ratios -- re-injections of one mix are one replicate,
#     never several;
#   * proteins with a ratio in >= min_reps replicates: moderated one-sample t-test against 0
#     (limma lmFit, intercept only; eBayes(robust = TRUE)), Benjamini-Hochberg per comparison;
#   * proteins detected in >= min_reps replicates of some channels and in none of another have
#     no ratio: they are listed as a detection pattern (no test), never dropped silently.
# =============================================================================

# The density peak of x: where most protein ratios sit. Compared with the median that sets the
# zero -- the median moves when one channel's strongest losses fall below detection.
silac_mode <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 10) return(NA_real_)
  d <- stats::density(x, n = 1024)
  d$x[which.max(d$y)]
}

# Rows of the report the analysis reads. `design`: one row per (Run, Channel) that holds a
# sample, with Replicate and Group. `cutoffs`: named q-value cutoffs, applied in order (a column
# the report lacks is never asked for). `cont`: TRUE for a contaminant row (any accession
# tagged). -> list(d = the kept rows with Replicate, Group and `ok` (the channel passed), counts =
# how many rows each filter removed, in order).
silac_prepare <- function(d, design, cutoffs, channel_q, cont) {
  d <- data.table::as.data.table(d)
  counts <- list(rows_read = nrow(d))
  keep <- d$Run %in% design$Run
  counts$rows_not_in_design <- sum(!keep)
  d <- d[keep]; cont <- cont[keep]
  counts$contaminant_rows <- sum(cont)
  counts$contaminant_protein_groups <- data.table::uniqueN(d$Protein.Group[cont])
  d <- d[!cont]
  for (col in names(cutoffs)) {
    n <- nrow(d)
    v <- d[[col]]
    d <- d[!is.na(v) & v <= cutoffs[[col]]]
    counts[[paste0("rows_failing_", col)]] <- n - nrow(d)
  }
  n <- nrow(d)
  d <- merge(d, data.table::as.data.table(design)[, .(Run, Channel, Replicate, Group)],
             by = c("Run", "Channel"))
  counts$rows_in_unassigned_channels <- n - nrow(d)
  d[, ok := !is.na(Channel.Q.Value) & Channel.Q.Value <= channel_q &
        !is.na(Precursor.Quantity) & Precursor.Quantity > 0]
  counts$rows_kept <- nrow(d)
  counts$rows_channel_passing <- sum(d$ok)
  list(d = d, counts = counts)
}

# Precursor quantities side by side: one row per (Run, Replicate, Protein.Group, Precursor.Id),
# one column per group, from the channel-passing rows only.
silac_wide <- function(d) {
  data.table::dcast(d[ok == TRUE], Run + Replicate + Protein.Group + Precursor.Id ~ Group,
                    value.var = "Precursor.Quantity", fun.aggregate = sum, fill = NA)
}

# One comparison num/den over the replicates `reps` (those holding both groups).
# -> list(table, offsets (per run: median and mode), n_tested, fit, mode_check).
silac_compare <- function(w, num, den, reps, min_reps, adjp = 0.05, mode_tol = 0.25) {
  r <- data.table::copy(w[Replicate %in% reps])
  if (!all(c(num, den) %in% names(r))) r <- r[0]
  r[, lr := if (nrow(r)) log2(get(num) / get(den)) else numeric(0)]
  r <- r[is.finite(lr)]
  # per run: the protein's median precursor ratio, centred on the run's median protein
  pr <- r[, .(ratio = stats::median(lr), n_prec = .N), by = .(Run, Replicate, Protein.Group)]
  off <- pr[, .(offset = stats::median(ratio), mode = silac_mode(ratio), n_proteins = .N),
            by = .(Run, Replicate)]
  pr <- merge(pr, off[, .(Run, offset, mode)], by = "Run")
  pr[, `:=`(centred = ratio - offset, centred_mode = ratio - mode)]
  # per replicate: the mean of its runs (one run: that run), and its distinct precursors
  rp <- pr[, .(centred = mean(centred), ratio = mean(ratio), centred_mode = mean(centred_mode),
               n_runs = .N), by = .(Replicate, Protein.Group)]
  rp <- merge(rp, r[, .(n_prec = data.table::uniqueN(Precursor.Id)), by = .(Replicate, Protein.Group)],
              by = c("Replicate", "Protein.Group"))
  wide <- function(v, fill = NA) {
    x <- data.table::dcast(rp, Protein.Group ~ Replicate, value.var = v, fill = fill)
    for (rpn in reps) if (!rpn %in% names(x)) x[, (rpn) := fill]
    x[, c("Protein.Group", reps), with = FALSE]
  }
  cw <- wide("centred"); uw <- wide("ratio"); nw <- wide("n_prec", 0L)
  m <- as.matrix(cw[, reps, with = FALSE]); rownames(m) <- cw$Protein.Group
  n_rep <- rowSums(!is.na(m))
  keep <- n_rep >= min_reps
  tab <- data.table::data.table(Protein.Group = cw$Protein.Group, n_replicates = n_rep,
                                logFC = NA_real_, t = NA_real_, P.Value = NA_real_,
                                adj.P.Val = NA_real_)
  fitrec <- list(df_prior = NA_real_, s2_prior = NA_real_, p_cutoff = NA_real_,
                 p_cutoff_kind = NA_character_, p_cutoff_basis = NA_character_,
                 critical_log2_raw_p05 = NA_real_, critical_log2_significant = NA_real_)
  fit_bh <- function(mm) {
    k <- rowSums(!is.na(mm)) >= min_reps
    if (sum(k) < 2 || length(reps) < min_reps) return(NULL)
    f <- limma::eBayes(limma::lmFit(mm[k, , drop = FALSE]), robust = TRUE)
    list(fit = f, tt = limma::topTable(f, coef = 1, number = Inf, sort.by = "none"))
  }
  res <- fit_bh(m)
  if (!is.null(res)) {
    tt <- res$tt; fit <- res$fit
    i <- match(rownames(tt), tab$Protein.Group)
    tab$logFC[i] <- tt$logFC; tab$t[i] <- tt$t
    tab$P.Value[i] <- tt$P.Value; tab$adj.P.Val[i] <- tt$adj.P.Val
    # The precision of a typical estimate, from the median posterior SD, df and replicates (a
    # rough reading, never a rule), and the cutoff the correction applied: Benjamini-Hochberg's
    # realised step-up cutoff adjp * k / m with k proteins significant of m tested -- or, with
    # none, adjp / m, BH's first step, its STRICTEST: what one protein changing alone needs. That
    # is a property of the cutoff, not of the replicates or the variance, and is never comparable
    # with a realised cutoff (a zero-call comparison once read as four-fold less sensitive than
    # one with ~800 calls, though its variance was the lower; statistics review 2026-10-16). No
    # "power" value: a true change found 80% of the time is not a detection limit, and rules
    # nothing out.
    se <- sqrt(stats::median(fit$s2.post)) / sqrt(stats::median(n_rep[keep]))
    dfm <- stats::median(fit$df.total)
    k_sig <- sum(tt$adj.P.Val < adjp)
    pst <- adjp * max(1, k_sig) / nrow(tt)
    fitrec <- list(
      df_prior = stats::median(fit$df.prior), s2_prior = stats::median(fit$s2.prior),
      p_cutoff = pst, p_cutoff_kind = if (k_sig) "realised" else "first_step",
      p_cutoff_basis = if (k_sig) sprintf("BH's realised cutoff adjp x k / m (%d significant of %d tested)",
                                          k_sig, nrow(tt))
                       else sprintf("adjp / m, BH's first step (none of %d tested significant)", nrow(tt)),
      # the half-width of a typical estimate's 95% interval -- also the OBSERVED |log2| that gives
      # p < 0.05 before the adjustment -- and the observed |log2| that passes p_cutoff
      critical_log2_raw_p05 = stats::qt(0.975, dfm) * se,
      # (upper tails with lower.tail = FALSE: a cutoff near 1e-50 is lost in 1 - p)
      critical_log2_significant = stats::qt(pst / 2, dfm, lower.tail = FALSE) * se)
  }
  # the zero point: per run, the median against the density peak (mode) of the uncentred ratios
  dm <- suppressWarnings(max(abs(off$offset - off$mode), na.rm = TRUE))
  mode_check <- list(max_abs_difference = if (is.finite(dm)) dm else NA_real_,
                     tolerance = mode_tol, n_significant_mode_centred = NA_integer_)
  if (is.finite(mode_check$max_abs_difference) && mode_check$max_abs_difference > mode_tol) {
    mm <- as.matrix(wide("centred_mode")[, reps, with = FALSE])
    alt <- fit_bh(mm)
    if (!is.null(alt)) {
      s <- alt$tt$adj.P.Val < adjp
      mode_check$n_significant_mode_centred <- sum(s)
      mode_check$n_up_mode_centred <- sum(s & alt$tt$logFC > 0)
      mode_check$n_down_mode_centred <- sum(s & alt$tt$logFC < 0)
    }
  }
  setnames <- data.table::setnames
  cc <- setnames(data.table::copy(cw), reps, paste0("centred_", reps))
  uu <- setnames(data.table::copy(uw), reps, paste0("uncentred_", reps))
  nn <- setnames(data.table::copy(nw), reps, paste0("n_precursors_", reps))
  tab <- merge(merge(merge(tab, cc, by = "Protein.Group"), uu, by = "Protein.Group"),
               nn, by = "Protein.Group")
  tab[, log2_ratio_uncentred_mean := rowMeans(as.matrix(.SD), na.rm = TRUE),
      .SDcols = paste0("uncentred_", reps)]
  # the fewest precursors behind any of the protein's replicate ratios
  pm <- as.matrix(tab[, paste0("n_precursors_", reps), with = FALSE])
  pm[pm == 0] <- NA
  tab[, min_precursors := suppressWarnings(apply(pm, 1, min, na.rm = TRUE))]
  tab[!is.finite(min_precursors), min_precursors := NA_real_]
  tab[n_replicates < min_reps, Evidence := sprintf(
    "ratio in %d of %d replicates: not tested (minimum %d)", n_replicates, length(reps), min_reps)]
  tab[n_replicates >= min_reps & !is.na(P.Value), Evidence := "tested"]
  tab[is.na(Evidence), Evidence := "not tested: too few protein groups had enough replicates to fit"]
  data.table::setorderv(tab, c("P.Value", "Protein.Group"), na.last = TRUE)
  list(table = tab, offsets = off, n_tested = sum(!is.na(tab$P.Value)), fit = fitrec,
       mode_check = mode_check)
}

# The precursor's C-terminal residue (K or R for a tryptic peptide), from Precursor.Id: the
# sequence with its modifications -- (...) and [...] -- and the charge removed.
silac_cterm <- function(id) {
  s <- sub("[0-9]+$", "", gsub("\\([^)]*\\)|\\[[^]]*\\]", "", id))
  substr(s, nchar(s), nchar(s))
}

# A label that reads lysine and arginine differently: per run and protein group, the median log2
# ratio num/den of its K-ending precursors minus that of its R-ending ones, over the replicates
# `reps`. A deuterated lysine (Lys4, 2H(4)) elutes slightly earlier than its light partner, so
# the two channels are not integrated over the same peak; a protein's ratio then depends on its
# mix of K- and R-ending peptides (a triple-SILAC experiment, 2026-10-16: M/L K-ending 0.29 log2
# below R-ending in every run). -> list(median (over protein-run pairs), per_run (named medians),
# n_pairs, n_proteins).
silac_residue_bias <- function(w, num, den, reps) {
  r <- data.table::copy(w[Replicate %in% reps])
  if (!all(c(num, den) %in% names(r)) || !nrow(r))
    return(list(median = NA_real_, per_run = list(), n_pairs = 0L, n_proteins = 0L))
  r[, lr := log2(get(num) / get(den))]
  r <- r[is.finite(lr)]
  r[, res := silac_cterm(Precursor.Id)]
  pr <- r[res %in% c("K", "R"), .(k = stats::median(lr[res == "K"]), r = stats::median(lr[res == "R"])),
          by = .(Run, Protein.Group)]
  pr <- pr[is.finite(k) & is.finite(r)]
  pr[, d := k - r]
  per_run <- pr[, .(d = stats::median(d)), by = Run]
  list(median = if (nrow(pr)) stats::median(pr$d) else NA_real_,
       per_run = stats::setNames(as.list(per_run$d), per_run$Run),
       n_pairs = nrow(pr), n_proteins = data.table::uniqueN(pr$Protein.Group))
}

# Detection of each protein group in each group's channel: the replicates where >= 1 of its
# precursors passed the channel q-value. `loaded`: list(group = replicates that hold it).
# -> one row per protein group: <group>_replicates, and the detection pattern -- every group
# either detected (>= min_reps replicates, or all that hold it) or not detected (0), at least one
# of each -- or NA. A pattern, not a test: not detected is not absent.
silac_detection <- function(d, loaded, min_reps) {
  det <- d[ok == TRUE, .(n = data.table::uniqueN(Replicate)), by = .(Protein.Group, Group)]
  wide <- data.table::dcast(det, Protein.Group ~ Group, value.var = "n", fill = 0L)
  groups <- names(loaded)
  for (g in groups) if (!g %in% names(wide)) wide[, (g) := 0L]
  wide <- wide[, c("Protein.Group", groups), with = FALSE]
  data.table::setnames(wide, groups, paste0(groups, "_replicates"))
  m <- as.matrix(wide[, paste0(groups, "_replicates"), with = FALSE])
  nl <- vapply(loaded, length, integer(1))
  need <- pmin(min_reps, nl)
  pres <- sweep(m, 2, need, ">=")
  abs_ <- m == 0
  ok <- rowSums(pres | abs_) == length(groups) & rowSums(pres) > 0 & rowSums(abs_) > 0
  wide[, pattern := NA_character_]
  if (any(ok)) wide$pattern[ok] <- vapply(which(ok), function(i) {
    p <- pres[i, ]
    sprintf("detected in %s; not detected in %s",
            paste(sprintf("%s (%d/%d)", groups[p], m[i, p], nl[p]), collapse = ", "),
            paste(sprintf("%s (0/%d)", groups[!p], nl[!p]), collapse = ", "))
  }, "")
  wide
}

# Each protein group's abundance in each group's channel: the median over replicates of its log2
# summed channel-passing precursor quantity, and its percentile among the protein groups detected
# in that channel -- so an abundant protein missing elsewhere reads differently from one near the
# detection limit.
silac_intensity <- function(d) {
  x <- d[ok == TRUE, .(q = sum(Precursor.Quantity)), by = .(Protein.Group, Group, Replicate)]
  x <- x[, .(log2_intensity = stats::median(log2(q))), by = .(Protein.Group, Group)]
  x[, percentile := round(100 * data.table::frank(log2_intensity, ties.method = "average") / .N),
    by = Group]
  data.table::dcast(x, Protein.Group ~ Group, value.var = c("log2_intensity", "percentile"))
}

# Each group's share of the channel-passing precursor signal (raw Precursor.Quantity), and the
# precursors and protein groups passing, per run.
silac_channel_summary <- function(d) {
  s <- d[ok == TRUE, .(signal = sum(Precursor.Quantity),
                       precursors = data.table::uniqueN(Precursor.Id),
                       protein_groups = data.table::uniqueN(Protein.Group)),
         by = .(Replicate, Run, Channel, Group)]
  s[, share_percent := 100 * signal / sum(signal), by = Run]
  data.table::setorderv(s, c("Replicate", "Run", "Channel"))
  s
}

# Label incorporation from label-check runs: a run holding ONE labelled channel only (the cells
# grown in labelled medium, injected alone). Per run, over precursors whose labelled channel passes
# `cutoffs` and channel_q with a non-zero quantity: the labelled share H / (H + L) of the raw
# quantity, as the median over precursors and summed over signal. The light row is taken as
# reported, whatever its own q-values (in a labelled-only sample it should fail them); a precursor
# with no light row (DIA-NN wrote none: the light channel did not pass --qvalue) is counted, not
# estimated from. Light signal is unlabelled protein -- or noise and interference -- and the
# precursors estimated from are the ones with light signal, so this is a LOWER bound.
# `d`: the check runs' rows, contaminants removed, no q-value filter. `checks`: list(run =
# labelled channel).
silac_incorporation <- function(d, checks, light, cutoffs, channel_q) {
  out <- list()
  for (run in names(checks)) {
    ch <- checks[[run]]
    x <- d[Run == run & Channel %in% c(ch, light)]
    lab <- x[Channel == ch & Channel.Q.Value <= channel_q & Precursor.Quantity > 0]
    for (col in names(cutoffs)) lab <- lab[!is.na(get(col)) & get(col) <= cutoffs[[col]]]
    okp <- unique(lab$Precursor.Id)
    w <- data.table::dcast(x, Precursor.Id ~ Channel, value.var = "Precursor.Quantity",
                           fun.aggregate = sum, fill = NA)
    for (cc in c(ch, light)) if (!cc %in% names(w)) w[, (cc) := NA_real_]
    w <- w[Precursor.Id %in% okp & !is.na(get(ch))]
    # The report holds a channel's row only when it passed --qvalue: no light row usually means
    # no confident light signal. Those precursors are counted, not estimated from.
    no_light <- sum(is.na(w[[light]]))
    w <- w[!is.na(get(light))]
    hv <- w[[ch]]; lv <- pmax(w[[light]], 0)
    out[[run]] <- list(run = run, channel = ch, precursors = nrow(w),
                       precursors_without_light_row = no_light,
                       median_percent = if (nrow(w)) 100 * stats::median(hv / (hv + lv)) else NA_real_,
                       signal_percent = if (nrow(w)) 100 * sum(hv) / sum(hv + lv) else NA_real_)
  }
  out
}
