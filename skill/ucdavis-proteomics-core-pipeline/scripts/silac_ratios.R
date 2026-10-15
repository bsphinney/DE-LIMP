# =============================================================================
# silac_ratios.R -- the SILAC ratio statistics, as plain functions (data.table + limma).
#
# run_silac.R sources this file for its numbers, and copies it VERBATIM into the analysis's
# reproducibility_log.R, so the replay runs the code that produced the tables -- not a second
# description of it (DE-LIMP rule 3: one definition).
#
# The method, as validated on a triple-SILAC experiment (DIA-NN 2.7.0 plexDIA, 4 replicates,
# 2026-10-14; references/silac.md):
#   * precursors at Q.Value <= q and channel-specific PG.Q.Value <= pg_q, contaminants removed;
#   * a precursor ratio needs both channels at Channel.Q.Value <= channel_q, with a non-zero raw
#     quantity (Precursor.Quantity: DIA-NN's channel normalisation is not used);
#   * per protein group and replicate, the median of its precursor log2 ratios;
#   * per replicate, median-centring (the uncentred ratios are kept: they carry the mixing
#     offset, e.g. one sample that put 3-5 times more protein into a 1:1:1-by-volume mix);
#   * proteins with a ratio in >= min_reps replicates: moderated one-sample t-test against 0
#     (limma lmFit, intercept only; eBayes(robust = TRUE)), Benjamini-Hochberg per comparison;
#   * proteins detected in >= min_reps replicates of some channels and in none of another have
#     no ratio: they are listed as channel-specific detections, never dropped silently.
# =============================================================================

# Rows of the report the analysis reads. `design`: one row per (Run, Channel) that holds a
# sample, with Replicate and Group. `cont`: TRUE for a contaminant row (any accession tagged).
# -> list(d = the kept rows with Replicate, Group and `ok` (the channel passed), counts = how
# many rows each filter removed, in order).
silac_prepare <- function(d, design, q, pg_q, channel_q, cont) {
  d <- data.table::as.data.table(d)
  counts <- list(rows_read = nrow(d))
  keep <- d$Run %in% design$Run
  counts$rows_not_in_design <- sum(!keep)
  d <- d[keep]; cont <- cont[keep]
  counts$contaminant_rows <- sum(cont)
  counts$contaminant_protein_groups <- data.table::uniqueN(d$Protein.Group[cont])
  d <- d[!cont]
  n <- nrow(d); d <- d[!is.na(Q.Value) & Q.Value <= q]; counts$rows_failing_q <- n - nrow(d)
  n <- nrow(d); d <- d[!is.na(PG.Q.Value) & PG.Q.Value <= pg_q]; counts$rows_failing_pg_q <- n - nrow(d)
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

# Precursor quantities side by side: one row per (Replicate, Protein.Group, Precursor.Id), one
# column per group, from the channel-passing rows only.
silac_wide <- function(d) {
  data.table::dcast(d[ok == TRUE], Replicate + Protein.Group + Precursor.Id ~ Group,
                    value.var = "Precursor.Quantity", fun.aggregate = sum, fill = NA)
}

# One comparison num/den over the replicates `reps` (those holding both groups).
# -> list(table, offsets, n_tested, fit = list(df_prior, s2_prior, detectable)).
silac_compare <- function(w, num, den, reps, min_reps) {
  r <- data.table::copy(w[Replicate %in% reps])
  if (!all(c(num, den) %in% names(r))) r <- r[0]
  r[, lr := if (nrow(r)) log2(get(num) / get(den)) else numeric(0)]
  r <- r[is.finite(lr)]
  pr <- r[, .(ratio = stats::median(lr), n_prec = .N), by = .(Replicate, Protein.Group)]
  off <- pr[, .(offset = stats::median(ratio)), by = Replicate]
  pr <- merge(pr, off, by = "Replicate")
  pr[, centred := ratio - offset]
  wide <- function(v, fill = NA) {
    x <- data.table::dcast(pr, Protein.Group ~ Replicate, value.var = v, fill = fill)
    for (rp in reps) if (!rp %in% names(x)) x[, (rp) := fill]
    x[, c("Protein.Group", reps), with = FALSE]
  }
  cw <- wide("centred"); uw <- wide("ratio"); nw <- wide("n_prec", 0L)
  m <- as.matrix(cw[, reps, with = FALSE]); rownames(m) <- cw$Protein.Group
  n_rep <- rowSums(!is.na(m))
  keep <- n_rep >= min_reps
  tab <- data.table::data.table(Protein.Group = cw$Protein.Group, n_replicates = n_rep,
                                logFC = NA_real_, t = NA_real_, P.Value = NA_real_,
                                adj.P.Val = NA_real_)
  fitrec <- list(df_prior = NA_real_, s2_prior = NA_real_, detectable = NA_real_,
                 detectable_alone = NA_real_)
  if (sum(keep) >= 2 && length(reps) >= min_reps) {
    fit <- limma::eBayes(limma::lmFit(m[keep, , drop = FALSE]), robust = TRUE)
    tt <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
    i <- match(rownames(tt), tab$Protein.Group)
    tab$logFC[i] <- tt$logFC; tab$t[i] <- tt$t
    tab$P.Value[i] <- tt$P.Value; tab$adj.P.Val[i] <- tt$adj.P.Val
    # How large a centred log2 ratio a typical protein needs (median posterior SD, median df and
    # replicates): to reach p < 0.05 before the multiple-testing adjustment, and to reach
    # BH-adjusted p < 0.05 when it is the only change among the m tested (p < 0.05 / m, BH's
    # first step). A rough reading of the power, for "not significant" -- never a rule.
    s_med <- sqrt(stats::median(fit$s2.post)) / sqrt(stats::median(n_rep[keep]))
    df_med <- stats::median(fit$df.total)
    fitrec <- list(df_prior = stats::median(fit$df.prior), s2_prior = stats::median(fit$s2.prior),
                   detectable = stats::qt(0.975, df_med) * s_med,
                   detectable_alone = stats::qt(1 - 0.025 / sum(keep), df_med) * s_med)
  }
  setnames <- data.table::setnames
  cc <- setnames(data.table::copy(cw), reps, paste0("centred_", reps))
  uu <- setnames(data.table::copy(uw), reps, paste0("uncentred_", reps))
  nn <- setnames(data.table::copy(nw), reps, paste0("n_precursors_", reps))
  tab <- merge(merge(merge(tab, cc, by = "Protein.Group"), uu, by = "Protein.Group"),
               nn, by = "Protein.Group")
  tab[, log2_ratio_uncentred_mean := rowMeans(as.matrix(.SD), na.rm = TRUE),
      .SDcols = paste0("uncentred_", reps)]
  tab[n_replicates < min_reps, Evidence := sprintf(
    "ratio in %d of %d replicates: not tested (minimum %d)", n_replicates, length(reps), min_reps)]
  tab[n_replicates >= min_reps & !is.na(P.Value), Evidence := "tested"]
  tab[is.na(Evidence), Evidence := "not tested: too few protein groups had enough replicates to fit"]
  data.table::setorderv(tab, c("P.Value", "Protein.Group"), na.last = TRUE)
  list(table = tab, offsets = off, n_tested = sum(!is.na(tab$P.Value)), fit = fitrec)
}

# Detection of each protein group in each group's channel: the replicates where >= 1 of its
# precursors passed the channel q-value. `loaded`: list(group = replicates that hold it).
# -> one row per protein group: <group>_replicates, and the presence/absence pattern -- every
# group either present (>= min_reps replicates) or absent (0), at least one of each -- or NA.
silac_detection <- function(d, loaded, min_reps) {
  det <- d[ok == TRUE, .(n = data.table::uniqueN(Replicate)), by = .(Protein.Group, Group)]
  wide <- data.table::dcast(det, Protein.Group ~ Group, value.var = "n", fill = 0L)
  groups <- names(loaded)
  for (g in groups) if (!g %in% names(wide)) wide[, (g) := 0L]
  wide <- wide[, c("Protein.Group", groups), with = FALSE]
  data.table::setnames(wide, groups, paste0(groups, "_replicates"))
  m <- as.matrix(wide[, paste0(groups, "_replicates"), with = FALSE])
  need <- pmin(min_reps, vapply(loaded, length, integer(1)))
  pres <- sweep(m, 2, need, ">=")
  abs_ <- m == 0
  ok <- rowSums(pres | abs_) == length(groups) & rowSums(pres) > 0 & rowSums(abs_) > 0
  wide[, pattern := NA_character_]
  if (any(ok)) wide$pattern[ok] <- apply(pres[ok, , drop = FALSE], 1, function(p)
    sprintf("in %s; not in %s", paste(groups[p], collapse = ", "),
            paste(groups[!p], collapse = ", ")))
  wide
}

# Each group's share of the channel-passing precursor signal (raw Precursor.Quantity), and the
# precursors and protein groups passing, per replicate.
silac_channel_summary <- function(d) {
  s <- d[ok == TRUE, .(signal = sum(Precursor.Quantity),
                       precursors = data.table::uniqueN(Precursor.Id),
                       protein_groups = data.table::uniqueN(Protein.Group)),
         by = .(Replicate, Run, Channel, Group)]
  s[, share_percent := 100 * signal / sum(signal), by = Replicate]
  data.table::setorderv(s, c("Replicate", "Channel"))
  s
}

# Label incorporation from label-check runs: a run holding ONE labelled channel only (the cells
# grown in labelled medium, injected alone). Per run, over precursors whose labelled channel passes
# q, pg_q and channel_q with a non-zero quantity: the labelled share H / (H + L) of the raw
# quantity, as the median over precursors and summed over signal. The light row is taken as
# reported, whatever its own q-values (in a labelled-only sample it should fail them); a precursor
# with no light row (DIA-NN wrote none: the light channel did not pass --qvalue) is counted, not
# estimated from. Light signal is unlabelled protein -- or noise and interference -- and the
# precursors estimated from are the ones with light signal, so this is a LOWER bound.
# `d`: the check runs' rows, contaminants removed, no q-value filter. `checks`: list(run =
# labelled channel).
silac_incorporation <- function(d, checks, light, q, pg_q, channel_q) {
  out <- list()
  for (run in names(checks)) {
    ch <- checks[[run]]
    x <- d[Run == run & Channel %in% c(ch, light)]
    okp <- unique(x[Channel == ch & Q.Value <= q & PG.Q.Value <= pg_q &
                      Channel.Q.Value <= channel_q & Precursor.Quantity > 0]$Precursor.Id)
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
