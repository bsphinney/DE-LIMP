#!/usr/bin/env Rscript
# =============================================================================
# sim_set_tests.R  --  the simulations behind run_sets.R's choices
# (references/set-tests.md has the numbers and what they mean).
#
#   Rscript sim_set_tests.R --what camera [--reps 300] [--seed 20260928]
#   Rscript sim_set_tests.R --what ip     [--reps 200] [--seed 20260928] [--out sim.csv]
#
# --what camera: camera's inter-protein correlation, estimated per set (run_sets.R's default)
#   or limma's preset 0.01, and fry. (a) Null sets whose proteins share a latent sample factor
#   (complexes, co-regulated proteins) against independent null sets, and power for a shifted
#   set, at 3 and 6 runs per group. (b) A GLOBAL run factor that every protein loads on (a
#   pulldown's recovery, a sample's loading), with a set-specific factor on one set and a shift
#   on another: random sets then correlate across runs, and camera's power falls with it.
#
# --what ip: a pulldown, simulated the way one is measured -- a protein's intensity in an IP is
#   its nonspecific background (in every IP and in the control) PLUS its specific signal (bait
#   interactors only), on the linear scale; specific enrichment runs from 0.3 to 3 log2 over the
#   background, measurement noise 0.4 log2. Old IPs recover the bait 1 log2 unit less; the Old
#   runs' background is 0 or +0.6 log2 higher (a condition-dependent background). 20 interactors
#   fall 0.8 relative to the complex in Old, 20 rise 0.8.
#   (1) interactome selection -- enriched at EVERY condition, pooled over conditions, pooled with
#       at least 1 or 2 log2 (2- or 4-fold; 2 is run_sets.R's default): how many members of the planted sets each
#       keeps, how often that is run_sets.R's minimum set size (P(tested); an untested set counts
#       as not significant), the between-condition tests after the interactome-median offset on
#       the planted sets and on null sets of random and of the 20 LEAST-enriched interactors, and
#       how often run_sets.R's "mostly weak interactors" flag (set_tests.R weak_interactors) fires
#       -- next to the old rule (within 1 log2 of the minimum) it replaced;
#   (2) the reference -- interactome median or bait protein -- with the pooled >= 2 log2
#       selection, in three scenarios: the planted sets only; "broad" (half the interactome +0.8
#       in Old); "bait_epitope" (the bait's own read-out -0.8 in Old, not the complex).
# =============================================================================
suppressMessages(library(limma))
.here <- local({
  f <- grep("^--file=", commandArgs(), value = TRUE)[1]
  if (is.na(f)) getwd() else dirname(normalizePath(sub("^--file=", "", f)))
})
source(file.path(.here, "set_tests.R"))
arg <- function(flag, default) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i) || i == length(a)) default else a[i + 1]
}
what <- arg("--what", "ip")
if (!what %in% c("camera", "ip")) stop("--what camera | ip")
reps <- as.integer(arg("--reps", if (what == "camera") "300" else "200"))
seed <- as.integer(arg("--seed", "20260928"))
out <- arg("--out", NA)
rate <- function(p) mean(p < 0.05, na.rm = TRUE)

# ---- --what camera ------------------------------------------------------------------
tests_of <- function(E, d, con, idx) {
  cs <- contrast_stats(lmFit(E, d), con); z <- z_from_t(cs$t, cs$df_total)
  eff <- lm_effects(E, d, con); U <- residual_U(eff)
  list(estimated = camera_fit(z, U, idx, ncol(U), NA)$PValue,
       preset_0.01 = camera_fit(z, U, idx, ncol(U), 0.01)$PValue,
       fry = fry_effects(eff, idx)$PValue, global = global_correlation(U))
}
camera_once <- function(n, shift, latent, G = 1000, m = 20) {
  grp <- factor(rep(c("A", "B"), each = n)); d <- model.matrix(~0 + grp); colnames(d) <- c("A", "B")
  E <- matrix(rnorm(G * 2 * n), G, 2 * n) + rnorm(G, 10, 2)
  E[1:m, ] <- E[1:m, ] + latent * matrix(rnorm(2 * n), m, 2 * n, byrow = TRUE)
  E[1:m, grp == "B"] <- E[1:m, grp == "B"] + shift
  r <- tests_of(E, d, c(A = -1, B = 1), list(1:m))
  unlist(r[c("estimated", "preset_0.01", "fry")])
}
global_once <- function(n, G = 1000) {
  grp <- factor(rep(c("A", "B"), each = n)); d <- model.matrix(~0 + grp); colnames(d) <- c("A", "B")
  E <- matrix(rnorm(G * 2 * n), G) + 0.6 * matrix(rnorm(2 * n), G, 2 * n, byrow = TRUE) * runif(G, 0.5, 1.5) +
       rnorm(G, 20, 2)                                    # every protein loads on one run factor
  E[1:20, ] <- E[1:20, ] + 0.7 * matrix(rnorm(2 * n), 20, 2 * n, byrow = TRUE)   # and set 1 on its own
  E[21:40, grp == "B"] <- E[21:40, grp == "B"] + 0.5                              # set 2 truly shifted
  r <- tests_of(E, d, c(A = -1, B = 1), list(specific = 1:20, shifted = 21:40, random = sample(41:G, 20)))
  data.frame(set = c("specific", "shifted", "random"), estimated = r$estimated,
             preset_0.01 = r$preset_0.01, fry = r$fry, global = r$global)
}
if (what == "camera") {
  set.seed(seed)
  grid <- expand.grid(latent = c(0, 0.7), shift = c(0, 0.4), n = c(3, 6))
  res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    p <- t(replicate(reps, camera_once(grid$n[i], grid$shift[i], grid$latent[i])))
    data.frame(grid[i, ], t(colMeans(p < 0.05)))
  }))
  cat(sprintf("sim_set_tests --what camera: %d replicates per row, seed %d, 1,000 proteins, one set of 20\n", reps, seed))
  cat("(a) rate of p < 0.05 (shift 0: false-positive rate, 0.05 = calibrated; shift 0.4: power)\n")
  cat("latent 0.7: the set's proteins share a sample factor (within-set correlation about 0.33)\n")
  print(res, row.names = FALSE, digits = 3)
  cat("\n(b) a global run factor on every protein (loading 0.6 x U(0.5, 1.5)); set 'specific' shares a further\n")
  cat("    factor, set 'shifted' moves 0.5 between groups, set 'random' is random; rate of p < 0.05\n")
  for (n in c(3, 6)) {
    g <- do.call(rbind, lapply(seq_len(reps), function(r) global_once(n)))
    for (st in c("specific", "random", "shifted")) {
      o <- g[g$set == st, ]
      cat(sprintf("  n = %d per group  %-9s camera estimated %.3f | camera preset 0.01 %.3f | fry %.3f   (random-set correlation %.2f)\n",
                  n, st, rate(o$estimated), rate(o$preset_0.01), rate(o$fry), mean(o$global)))
    }
  }
  quit(save = "no")
}

# ---- --what ip ----------------------------------------------------------------------
# One simulated pulldown. Returns, per selection rule and reference, the tests of the planted and
# null sets after the offset.
ip_once <- function(bg_old, scenario = "planted", n = 3, n_int = 400, n_bg = 2000, beta = 0.8,
                    rules = c("every_condition", "pooled", "pooled_min1", "pooled_min2"), refs = "interactome") {
  G <- n_int + n_bg; int <- seq_len(n_int); bait <- 1L
  grp <- factor(rep(c("Old_IP", "Young_IP", "Old_IgG", "Young_IgG"), each = n),
                levels = c("Old_IP", "Young_IP", "Old_IgG", "Young_IgG"))
  ip <- grp %in% c("Old_IP", "Young_IP"); old <- grp %in% c("Old_IP", "Old_IgG")
  a <- rnorm(G, 20, 2)
  e <- c(3.5, runif(n_int - 1, 0.3, 3), rep(-Inf, n_bg))     # the bait (row 1) strongly enriched
  rec <- ifelse(grp == "Old_IP", -1, 0) + rnorm(length(grp), 0, 0.3)
  bg <- ifelse(old, bg_old, 0) + rnorm(length(grp), 0, 0.1)
  comp <- rep(0, G); down <- 2:21; up <- 22:41
  comp[down] <- -beta; comp[up] <- beta
  if (scenario == "broad") comp[setdiff(2:(n_int / 2 + 1), c(down, up))] <- beta
  if (scenario == "bait_epitope") comp[bait] <- -beta
  lin <- 2^(outer(a, bg, "+"))
  spec <- 2^(outer(a + e, rec, "+") + outer(comp, as.numeric(grp == "Old_IP")))
  lin[, ip] <- lin[, ip] + spec[, ip]
  E <- log2(lin) + matrix(rnorm(G * length(grp), 0, 0.4), G)
  design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
  cm <- makeContrasts(Old_IP - Old_IgG, Young_IP - Young_IgG, levels = design)
  fit <- lmFit(E, design)
  fc <- eBayes(contrasts.fit(fit, cm))
  adj <- apply(fc$p.value, 2, p.adjust, "BH")
  pooled <- contrast_stats(fit, pooled_enrichment_contrast(cm, 1:2))
  padj <- p.adjust(2 * pt(-abs(pooled$t), pooled$df_total), "BH")
  sel <- list(every_condition = which(fc$coefficients[, 1] > 0 & adj[, 1] < 0.05 &
                                      fc$coefficients[, 2] > 0 & adj[, 2] < 0.05),
              pooled = which(pooled$logFC > 0 & padj < 0.05),
              pooled_min1 = which(pooled$logFC >= 1 & padj < 0.05),
              pooled_min2 = which(pooled$logFC >= 2 & padj < 0.05))[rules]
  runs <- which(ip)
  d2 <- model.matrix(~0 + droplevels(grp[runs])); colnames(d2) <- c("Old_IP", "Young_IP")
  con <- c(Old_IP = 1, Young_IP = -1)
  out <- list()
  for (rule in rules) for (ref in refs) {
    s <- sel[[rule]]
    if (ref == "bait" && !(bait %in% s)) next
    uni <- if (ref == "bait") setdiff(s, bait) else s
    unpl <- setdiff(intersect(uni, int), c(bait, down, up, which(comp != 0)))
    weak <- unpl[order(e[unpl])][seq_len(min(20, length(unpl)))]
    sets <- list(planted_down = intersect(down, uni), planted_up = intersect(up, uni),
                 null_random = if (length(unpl) >= 20) sample(unpl, 20) else unpl, null_weakest = weak)
    off <- ip_offsets(E, runs, s, ref, bait)
    En <- sweep(E[uni, runs, drop = FALSE], 2, off)
    cs <- contrast_stats(lmFit(En, d2), con)
    m <- list(E = En, design = d2, contrast = con, weights = NULL, block = NULL, correlation = NULL,
              z = z_from_t(cs$t, cs$df_total))
    idx <- reindex(sets, uni)
    ok <- lengths(idx) >= SET_DEFAULTS$min_size           # run_sets.R tests no smaller set
    tt <- if (any(ok)) test_sets(m, idx[ok]) else NULL
    minr <- c(every_condition = 0, pooled = 0, pooled_min1 = 1, pooled_min2 = 2)[[rule]]
    weak <- logical(G); weak[s] <- weak_interactors(pooled$logFC[s])
    frac <- function(flag) vapply(idx, function(i) if (length(i)) mean(flag[uni[i]]) else NA_real_, 0)
    o <- data.frame(scenario = scenario, bg_old = bg_old, rule = rule, reference = ref, set = names(sets),
                    kept = lengths(idx), tested = ok, camera_p = NA_real_, fry_p = NA_real_,
                    weak_flag = frac(weak) > 0.5, weak_flag_old = frac(pooled$logFC < minr + 1) > 0.5,
                    n_interactome = length(s), n_true = length(intersect(s, int)))
    if (!is.null(tt)) o[ok, c("camera_p", "fry_p")] <- tt[, c("camera_PValue", "fry_PValue")]
    out[[paste(rule, ref)]] <- o
  }
  do.call(rbind, out)
}

set.seed(seed)
cat(sprintf("sim_set_tests --what ip: %d replicates per condition, seed %d; 400 interactors (0.3-3 log2 over background), 2,000 background proteins, 3 IPs per condition\n", reps, seed))
sel <- do.call(rbind, lapply(c(0, 0.6), function(bg) do.call(rbind, lapply(seq_len(reps), function(r) ip_once(bg)))))
sig <- function(p) mean(!is.na(p) & p < 0.05)            # an untested set is not significant
cat(sprintf(paste0("\n(1) interactome selection (interactome-median reference). kept: members of the planted set of 20 in the\n",
                   "    interactome; tested: kept >= %d (run_sets.R's minimum); fry / camera: rate of p < 0.05 over all\n",
                   "    replicates (planted: power of the whole pipeline; null: false-positive rate); weak: flagged 'mostly\n",
                   "    weak interactors' (over half in the interactome's bottom quarter of enrichment), old: by the old\n",
                   "    rule (within 1 log2 of the minimum); fry unflagged: significant and not flagged\n"), SET_DEFAULTS$min_size))
for (bg in c(0, 0.6)) {
  o <- sel[sel$bg_old == bg, ]
  cat(sprintf("  Old background %+.1f log2; interactome size: every_condition %.0f, pooled %.0f, pooled >= 1 log2 %.0f, >= 2 log2 %.0f\n", bg,
              mean(o$n_interactome[o$rule == "every_condition"]), mean(o$n_interactome[o$rule == "pooled"]),
              mean(o$n_interactome[o$rule == "pooled_min1"]), mean(o$n_interactome[o$rule == "pooled_min2"])))
  for (st in c("planted_down", "planted_up", "null_random", "null_weakest")) for (rl in c("every_condition", "pooled", "pooled_min1", "pooled_min2")) {
    x <- o[o$set == st & o$rule == rl, ]
    cat(sprintf("    %-13s %-16s kept %5.1f  tested %.3f   fry %.3f   camera %.3f   weak %.3f  old %.3f   fry unflagged %.3f\n",
                st, rl, mean(x$kept), mean(x$tested), sig(x$fry_p), sig(x$camera_p), mean(x$weak_flag, na.rm = TRUE),
                mean(x$weak_flag_old, na.rm = TRUE), mean(!is.na(x$fry_p) & x$fry_p < 0.05 & !x$weak_flag)))
  }
}
nr <- sel[sel$set == "null_random" & sel$rule == "pooled_min2", ]
cat(sprintf("  null_random fry false-positive rate at the default (pooled >= 2 log2), both backgrounds: %.4f (SE %.4f, n %d)\n",
            sig(nr$fry_p), sqrt(0.05 * 0.95 / nrow(nr)), nrow(nr)))
refr <- do.call(rbind, lapply(c("planted", "broad", "bait_epitope"), function(sc) do.call(rbind, lapply(seq_len(reps), function(r)
  ip_once(0.6, sc, rules = "pooled_min2", refs = c("interactome", "bait"))))))
cat("\n(2) the reference, with the pooled >= 2 log2 interactome and +0.6 Old background; rate of p < 0.05 (untested = not significant)\n")
for (sc in c("planted", "broad", "bait_epitope")) for (st in c("planted_down", "planted_up", "null_random", "null_weakest")) {
  x <- refr[refr$scenario == sc & refr$set == st, ]
  cat(sprintf("  %-12s %-13s interactome: fry %.3f camera %.3f | bait: fry %.3f camera %.3f\n", sc, st,
              sig(x$fry_p[x$reference == "interactome"]), sig(x$camera_p[x$reference == "interactome"]),
              sig(x$fry_p[x$reference == "bait"]), sig(x$camera_p[x$reference == "bait"])))
}
if (!is.na(out)) utils::write.csv(rbind(sel, refr), out, row.names = FALSE)
