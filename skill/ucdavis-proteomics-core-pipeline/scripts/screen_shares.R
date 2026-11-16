# =============================================================================
# screen_shares.R -- the protein screen's shares: how much each sample's protein value rests on each
# found precursor value, as the derivative of limpa's fit of that protein. One definition (repo rule
# 3), sourced by make_screen_data.R and by its tests (test_screen_data.py).
#
# It mirrors limpa 1.4.0's objective in peptides2ProteinBFGS: the found values enter only through
# sum(((y - X beta) / sigma)^2), with limpa's design X (shares_design); missing values enter through
# the detection curve, which does not involve y. A limpa that changes either would make these
# derivatives wrong while limpa's own fit still reproduces, so make_screen_data.R checks them on every
# run against nudges of fully converged fits (shares_nudge) and withholds all shares when they differ.
#
# The detection chances (design section 5.7) come from the same fit: the chance limpa's fit gives
# each precursor of being found in each sample (detection_at_fit). They mirror the missing-value term
# of that objective, so each protein's are checked against limpa's own objective value at its fit
# (detection_objective) and withheld when it differs.
#
# h, everywhere: list(dpc, prior.mean, prior.sd, prior.logFC), as dpcQuant used them.
# =============================================================================

# limpa's start values for one protein (peptides2Proteins, limpa 1.4.0), from its rows yi of the
# whole-matrix imputation. A protein quantified apart from the others gets other start values, and
# BFGS then stops up to ~1e-5 from the delivered value.
shares_start <- function(E, yi) {
  ns <- ncol(E); np <- nrow(E)
  if (np > 1L) {
    beta <- numeric(ns + np - 1L)
    beta[seq_len(ns)] <- colMeans(yi)
    b <- rowMeans(yi); b <- b - mean(b)
    beta[ns + seq_len(np - 1L)] <- b[-np]
  } else beta <- drop(yi)
  imp <- which(colSums(!is.na(E)) == 0L)
  if (length(imp) > 1L) beta[imp] <- mean(beta[imp])
  beta
}

# limpa's fit of one protein with dpcQuant's settings, plus newton.polish (start NULL: limpa's own
# per-protein start). limpa 1.4.0 applies the polish, and returns the Hessian, only inside its
# standard-errors branch, hence standard.errors = TRUE (dpcQuant uses it too). $par.bgfs is where BFGS
# stopped (dpcQuant's value: it runs newton.polish = FALSE); $par is one Newton step on; $hessian is
# at $par.bgfs.
shares_fit <- function(E, sg, h, start = NULL)
  limpa::peptides2ProteinBFGS(E, sigma = sg, dpc = h$dpc, prior.mean = h$prior.mean,
    prior.sd = h$prior.sd, prior.logFC = h$prior.logFC, standard.errors = TRUE,
    newton.polish = TRUE, start = start)

# limpa's design for one protein's np x ns matrix (peptides2ProteinBFGS): sample effects, then
# sum-to-zero precursor effects; row i is as.vector(E)[i].
shares_design <- function(np, ns) {
  Samples <- factor(rep.int(seq_len(ns), rep.int(np, ns)))
  if (np < 2L) return(stats::model.matrix(~0 + Samples))
  Precursors <- factor(rep.int(seq_len(np), ns))
  stats::contrasts(Precursors) <- stats::contr.sum(np)
  stats::model.matrix(~0 + Samples + Precursors)
}

# The derivatives of limpa's protein values at its optimum, given fit f1 of E (shares_fit). A second
# fit, started from f1's polished optimum, gives limpa's Hessian H at the optimum (its BFGS stop, which
# must be within settle_tol of its own polish). Moving one found value y_i moves the optimum by
# d beta / d y_i = H^-1 (2 / sigma^2) X[i, ] (implicit-function theorem). Returns list(Wm), Wm an
# ns x (np * ns) matrix whose column p + (k - 1) * np is d values / d E[p, k] (0 at a gap), or
# list(withheld = <reason>).
shares_derivative <- function(E, sg, h, f1, settle_tol = 1e-5) {
  ns <- ncol(E); np <- nrow(E)
  if (sum(rowSums(!is.na(E)) == 0L) > 1L)    # limpa merges those rows, so its design is not shares_design()
    return(list(withheld = "two or more of its precursors have no value in any sample; shares are not computed for that"))
  f2 <- shares_fit(E, sg, h, drop(f1$par))
  moved <- max(abs(f2$par - f2$par.bgfs))
  if (!is.finite(moved) || moved > settle_tol)
    return(list(withheld = sprintf("limpa's fit did not settle at its optimum (a Newton step still moves it by %s)",
                                   format(signif(moved, 3)))))
  X <- shares_design(np, ns)
  if (!identical(dim(f2$hessian), c(ncol(X), ncol(X)))) stop("limpa's Hessian does not match its design")
  obs <- which(!is.na(as.vector(E)))
  Wm <- matrix(0, ns, np * ns)
  Wm[, obs] <- solve(f2$hessian, (2 / sg^2) * t(X[obs, , drop = FALSE]))[seq_len(ns), , drop = FALSE]
  list(Wm = Wm)
}

# own[p, j] = d value_j / d E[p, j] (NA at a gap), total[j] = sum of all of value_j's derivatives,
# other[j] = total[j] - sum over p of own[p, j].
shares_split <- function(Wm, E) {
  ns <- ncol(E); np <- nrow(E)
  own <- matrix(Wm[cbind(rep(seq_len(ns), each = np), seq_len(np * ns))], np, ns)
  own[is.na(E)] <- NA
  total <- rowSums(Wm)
  list(own = own, other = total - colSums(own, na.rm = TRUE), total = total)
}

# limpa's protein values at its optimum, converged fully: limpa's fit from its own start, then refit
# from the polished point until a Newton step moves it less than 1e-12 (at most 4 times).
shares_converged <- function(E, sg, h) {
  f <- shares_fit(E, sg, h)
  for (r in 1:4) {
    f <- shares_fit(E, sg, h, drop(f$par))
    if (max(abs(f$par - f$par.bgfs)) < 1e-12) break
  }
  drop(f$par[seq_len(ncol(E))])
}

# The same derivatives by nudging: a central difference of fully converged fits, one column per
# found value in `cells` (indices into as.vector(E)), the columns of shares_derivative's Wm.
shares_nudge <- function(E, sg, h, cells, step = 1e-3)
  matrix(vapply(cells, function(i) {
    up <- E; dn <- E; up[i] <- up[i] + step; dn[i] <- dn[i] - step
    (shares_converged(up, sg, h) - shares_converged(dn, sg, h)) / (2 * step)
  }, numeric(ncol(E))), ncol(E))

# The detection chances of E at the point beta of limpa's fit (sample effects, then the sum-to-zero
# precursor effects of shares_design). mu[p, j] is the fitted level of precursor p in sample j; G is
# the chance of missing it that limpa's objective uses at a gap, over limpa's 16-point normal Gauss
# quadrature: G = sum_k plogis(b0 + b1 (mu + sg z_k), lower = FALSE) w_k. pf_fit = 1 - G (the weights
# sum to 1: = sum_k plogis(b0 + b1 (mu + sg z_k)) w_k), for every cell; pf_obs = plogis(b0 + b1 E)
# at a found value, NA at a gap. All np x ns.
detection_chances <- function(E, beta, sg, dpc) {
  np <- nrow(E); ns <- ncol(E)
  mu <- drop(shares_design(np, ns) %*% beta)
  gq <- statmod::gauss.quad.prob(16, dist = "normal")
  eta <- dpc[1] + dpc[2] * (mu + sg * matrix(gq$nodes, length(mu), 16L, byrow = TRUE))
  G <- matrix(drop(stats::plogis(eta, lower.tail = FALSE) %*% gq$weights), np, ns)
  list(pf_fit = 1 - G, pf_obs = stats::plogis(dpc[1] + dpc[2] * E), G = G)
}

# limpa's objective at beta (peptides2ProteinBFGS's minusLogPosterior, limpa 1.4.0), rebuilt with the
# chances of missing G at the gaps: found residuals + the detection term + the prior on the values.
detection_objective <- function(E, beta, sg, h, G) {
  ns <- ncol(E); y <- as.vector(E); gap <- is.na(y)
  mu <- drop(shares_design(nrow(E), ns) %*% beta)
  v <- beta[seq_len(ns)]; m <- mean(v)
  sum(((y[!gap] - mu[!gap]) / sg)^2) - 2 * sum(log(G[gap])) +
    (m - h$prior.mean)^2 / h$prior.sd^2 + sum((v - m)^2) / h$prior.logFC^2
}

# The detection chances of limpa's fit f of E (shares_fit): at f$par.bgfs, where BFGS stopped -- the
# point dpcQuant delivers (it runs newton.polish = FALSE), not the polished f$par. Returns
# detection_chances' list plus objective_diff, the relative difference between the rebuilt objective
# and limpa's own value there (f$value.bfgs); or list(withheld = <reason>) when that exceeds tol.
detection_at_fit <- function(E, f, sg, h, tol = 1e-8) {
  beta <- drop(f$par.bgfs); v <- f$value.bfgs
  if (!is.numeric(v) || length(v) != 1L)
    return(list(withheld = "limpa's fit did not return its objective value where BFGS stopped (value.bfgs)"))
  d <- detection_chances(E, beta, sg, h$dpc)
  diff <- abs(detection_objective(E, beta, sg, h, d$G) - v) / max(1, abs(v))
  if (!is.finite(diff) || diff > tol)
    return(list(withheld = sprintf(paste("limpa's objective at its fit differs by %s (relative) from the one",
                                         "rebuilt from the detection chances, so they may not be the ones limpa used"),
                                   format(signif(diff, 3)))))
  d$objective_diff <- diff
  d
}
