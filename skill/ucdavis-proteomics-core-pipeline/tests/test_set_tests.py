#!/usr/bin/env python3
"""Protein-set tests: set_tests.R's functions and run_sets.R end to end.

  * camera_fit() is limma::camera() exactly with no block and no weights (estimated or preset
    inter-gene correlation), and limma::cameraPR() on the z-scores with a preset correlation --
    so the blocked route is camera's own algorithm on the blocked fit; fry_fit() is limma::fry()
    up to the basis its robust variance sees; both work on one canonical residual basis, so
    constant weights change nothing and neither does how the design is parametrised;
  * a planted shift is found by camera and fry; null sets give valid p-values (fry uniform,
    camera's estimated-correlation mode conservative), and a correlated null set is kept near
    5% by the estimated mode but not by the 0.01 preset -- why run_sets.R estimates it;
  * a set that tracks run depth is significant alone and "lost" with depth in the model;
  * a random block: fry takes block + correlation, and whitening lowers the residual
    correlation camera sees;
  * presence-call fractions and the majority flag; GMT reading; size limits; the organism;
  * run_de.R writes set_test_inputs.rds and suggests set tests below --sets-trigger; run_sets.R
    on a synthetic pulldown (random Mouse block, --ip-map, a GMT) writes the tables, the
    bait-normalised table, the record and the methods, and refuses a GMT without an origin;
  * setup.sh installs the GO annotation packages (conda first, then Bioconductor) and says so in
    setup.json -- never fatal, since DE does not need them.

Skips without R (and, for the end-to-end class, without limpa / arrow / dplyr / tidyr).
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, HERE)

from test_run_de_contaminants import r_has, rscript   # noqa: E402
import test_setup_limpa as setup_limpa                 # noqa: E402  (module: its tests are not re-collected)

SET_TESTS_R = os.path.join(SCRIPTS, "set_tests.R")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")
RUN_SETS = os.path.join(SCRIPTS, "run_sets.R")
SIM = os.path.join(SCRIPTS, "sim_set_tests.R")
PRE = f'suppressMessages(library(limma)); source("{SET_TESTS_R}")\n'


def r_json(body):
    """Run R code that ends by printing one JSON value; return it parsed."""
    out = rscript(PRE + body)
    return json.loads(out.strip().splitlines()[-1])


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@unittest.skipUnless(r_has("limma", "statmod", "jsonlite"), "needs R with limma, statmod and jsonlite")
class CameraIsCamera(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = r_json(r'''
set.seed(1); G <- 600; n <- 12
grp <- factor(rep(c("A", "B", "C"), each = 4)); design <- model.matrix(~0 + grp); colnames(design) <- levels(grp)
E <- matrix(rnorm(G * n), G, n, dimnames = list(paste0("p", 1:G), NULL)); E[1:30, 5:8] <- E[1:30, 5:8] + 0.8
W <- matrix(runif(G * n, 0.3, 2), G, n)
idx <- list(s1 = 1:30, s2 = 31:70, s3 = sample(G, 50))
con <- makeContrasts(B - A, levels = design)[, 1]
z_of <- function(fit, c) { cs <- contrast_stats(fit, c); z_from_t(cs$t, cs$df_total) }
# no weights: camera() and fry() exactly
d_cam <- c()
for (cor in list(NA, 0.01)) {
  ref <- camera(E, idx, design, contrast = con, inter.gene.cor = cor, sort = FALSE)
  eff <- lm_effects(E, design, con)
  mine <- camera_fit(z_of(lmFit(E, design), con), residual_U(eff), idx, ncol(eff) - 1, cor)
  d_cam <- c(d_cam, max(abs(ref$PValue - mine$PValue)), sum(ref$Direction != mine$Direction))
}
d_fry <- max(abs(fry(E, idx, design, con, sort = "none")$PValue - fry_fit(E, idx, design, con)$PValue))
blk <- rep(1:4, 3); rho <- duplicateCorrelation(E, design, block = blk)$consensus
d_fry_blk <- max(abs(fry(E, idx, design, con, block = blk, correlation = rho, sort = "none")$PValue -
                     fry_fit(E, idx, design, con, NULL, blk, rho)$PValue))
# without weights the ONLY basis-dependent step is the robust variance's largest coordinate:
# taken from limma's own effects, ours is limma::fry exactly (block or not)
d_fry_u2 <- max(abs(fry(E, idx, design, con, sort = "none")$PValue -
                    fry_effects(lm_effects(E, design, con), idx,
                                u2_eff = limma:::.lmEffects(E, design = design, contrast = con))$PValue))
d_fry_u2_blk <- max(abs(fry(E, idx, design, con, block = blk, correlation = rho, sort = "none")$PValue -
                        fry_effects(lm_effects(E, design, con, NULL, blk, rho), idx,
                                    u2_eff = limma:::.lmEffects(E, design = design, contrast = con,
                                                                block = blk, correlation = rho))$PValue))
z <- z_of(lmFit(E, design, block = blk, correlation = rho), con)
U <- residual_U(lm_effects(E, design, con, NULL, blk, rho))
pr <- max(abs(camera_fit(z, U, idx, ncol(U), 0.01)$PValue - cameraPR(z, idx, inter.gene.cor = 0.01, sort = FALSE)$PValue))
# constant weights: the common basis is the unweighted one
Wc <- matrix(2, G, n)
m0 <- list(E = E, design = design, contrast = con, weights = NULL, block = NULL, correlation = NULL, z = z_of(lmFit(E, design), con))
m1 <- modifyList(m0, list(weights = Wc)); m1$z <- z_of(lmFit(E, design, weights = Wc), con)
a <- test_sets(m0, idx); b <- test_sets(m1, idx)
d_const <- max(abs(c(a$camera_PValue - b$camera_PValue, a$fry_PValue - b$fry_PValue)))
# varying weights: the same answer whatever the design's parametrisation (column order and
# sign are arbitrary, as limma's per-protein QR is), and close to limma's own weighted camera
mw <- modifyList(m0, list(weights = W)); mw$z <- z_of(lmFit(E, design, weights = W), con)
design2 <- -design[, c(3, 1, 2)]; con2 <- -con[c(3, 1, 2)]
mw2 <- modifyList(mw, list(design = design2, contrast = con2))
x <- test_sets(mw, idx); y <- test_sets(mw2, idx)
d_param <- c(max(abs(x$camera_PValue - y$camera_PValue)), max(abs(x$fry_PValue - y$fry_PValue)))
lc <- camera(E, idx, design, contrast = con, weights = W, inter.gene.cor = NA, sort = FALSE)
cat(jsonlite::toJSON(list(camera = d_cam, fry = d_fry, fry_block = d_fry_blk, camerapr = pr,
                          fry_u2 = d_fry_u2, fry_u2_block = d_fry_u2_blk,
                          const = d_const, param = d_param,
                          limma_weighted = cor(log(lc$PValue), log(x$camera_PValue))),
                     auto_unbox = TRUE, digits = NA))
''')

    def test_no_weights_is_limma_camera_and_fry(self):
        self.assertLess(max(self.r["camera"]), 1e-12, self.r)
        # fry: with the largest-coordinate term from limma's basis, limma::fry exactly
        self.assertLess(self.r["fry_u2"], 1e-12, self.r)
        self.assertLess(self.r["fry_u2_block"], 1e-12, self.r)
        # and on the canonical basis alone, close (that term only)
        self.assertLess(self.r["fry"], 0.05, self.r)
        self.assertLess(self.r["fry_block"], 0.05, self.r)

    def test_blocked_preset_is_camerapr(self):
        self.assertLess(self.r["camerapr"], 1e-12)

    def test_weights_on_one_basis(self):
        self.assertLess(self.r["const"], 1e-10)         # constant weights change nothing
        self.assertLess(max(self.r["param"]), 1e-10)    # nor does how the design is written
        self.assertGreater(self.r["limma_weighted"], 0.95)


@unittest.skipUnless(r_has("limma", "jsonlite"), "needs R with limma and jsonlite")
class PlantedNullAndCorrelation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = r_json(r'''
set.seed(20260928); G <- 1000
one <- function(n, shift, latent) {
  grp <- factor(rep(c("A", "B"), each = n)); d <- model.matrix(~0 + grp); colnames(d) <- c("A", "B")
  con <- c(A = -1, B = 1)
  E <- matrix(rnorm(G * 2 * n), G, 2 * n) + rnorm(G, 10, 2)
  E[1:20, ] <- E[1:20, ] + latent * matrix(rnorm(2 * n), 20, 2 * n, byrow = TRUE)
  E[1:20, grp == "B"] <- E[1:20, grp == "B"] + shift
  m <- list(E = E, design = d, contrast = con, weights = NULL, block = NULL, correlation = NULL)
  cs <- contrast_stats(lmFit(E, d), con); m$z <- z_from_t(cs$t, cs$df_total)
  a <- test_sets(m, list(s = 1:20)); b <- test_sets(m, list(s = 1:20), 0.01)
  c(a$camera_PValue, b$camera_PValue, a$fry_PValue, a$camera_Direction == "Up", a$fry_Direction == "Up")
}
null <- t(replicate(300, one(4, 0, 0))); cor <- t(replicate(300, one(4, 0, 0.7)))
planted <- one(4, 1.5, 0)
ks <- function(p) suppressWarnings(ks.test(p, "punif")$p.value)
cat(jsonlite::toJSON(list(
  null_rate = colMeans(null[, 1:3] < 0.05), null_ks = c(ks(null[, 2]), ks(null[, 3])),
  cor_rate = colMeans(cor[, 1:3] < 0.05), planted = planted), digits = NA))
''')

    def test_planted_shift_found_by_both(self):
        cam, _, fry, cam_up, fry_up = self.r["planted"]
        self.assertLess(cam, 1e-3)
        self.assertLess(fry, 1e-3)
        self.assertEqual((cam_up, fry_up), (1, 1))

    def test_null_p_values_valid(self):
        est, preset, fry = self.r["null_rate"]
        # fry and camera with a preset correlation: uniform under independence
        self.assertGreater(self.r["null_ks"][0], 0.001, "camera (preset) null p not uniform")
        self.assertGreater(self.r["null_ks"][1], 0.001, "fry null p not uniform")
        self.assertTrue(0.01 < fry < 0.10, fry)
        # camera estimating the correlation: valid, and conservative at small n
        self.assertLessEqual(est, 0.07)

    def test_correlated_null_is_why_the_correlation_is_estimated(self):
        est, preset, fry = self.r["cor_rate"]
        self.assertLess(est, 0.15)
        self.assertGreater(preset, 0.3)
        self.assertLess(fry, 0.12)


@unittest.skipUnless(r_has("limma", "jsonlite"), "needs R with limma and jsonlite")
class DepthBlockAndPresence(unittest.TestCase):
    def test_depth_tracking_set_is_lost_with_depth(self):
        r = r_json(r'''
set.seed(5); G <- 800; n <- 5
grp <- factor(rep(c("A", "B"), each = n)); d <- model.matrix(~0 + grp); colnames(d) <- c("A", "B")
con <- c(A = -1, B = 1)
depth <- c(rnorm(n, 0, 0.3), rnorm(n, 0.9, 0.3))          # B runs deeper
E <- matrix(rnorm(G * 2 * n, 0, 0.3), G, 2 * n) + rnorm(G, 10, 2)
E[1:25, ] <- E[1:25, ] + 1.2 * matrix(depth, 25, 2 * n, byrow = TRUE)   # tracks depth, not group
idx <- list(depth_set = 1:25, other = 26:60)
f0 <- fit_like_run_de("maxlfq", E, d); c0 <- contrast_stats(f0$fit, con)
m0 <- list(E = E, design = d, contrast = con, weights = NULL, block = NULL, correlation = NULL,
           z = z_from_t(c0$t, c0$df_total))
a <- add_depth(d, con, depth - mean(depth))
f1 <- fit_like_run_de("maxlfq", E, a$design); c1 <- contrast_stats(f1$fit, a$contrast)
m1 <- list(E = E, design = a$design, contrast = a$contrast, weights = NULL, block = NULL, correlation = NULL,
           z = z_from_t(c1$t, c1$df_total))
t0 <- test_sets(m0, idx); t1 <- test_sets(m1, idx)
cat(jsonlite::toJSON(list(
  fry0 = t0$fry_FDR[1], fry1 = t1$fry_FDR[1],
  holds = set_holds(t0$fry_FDR, t0$fry_Direction, t1$fry_FDR, t1$fry_Direction, 0.05),
  cam_holds = set_holds(t0$camera_FDR, t0$camera_Direction, t1$camera_FDR, t1$camera_Direction, 0.05)),
  auto_unbox = TRUE, digits = NA))
''')
        self.assertLess(r["fry0"], 0.05)
        self.assertGreater(r["fry1"], 0.05)
        self.assertEqual(r["holds"][0], "lost")
        self.assertIn(r["cam_holds"][0], ("lost", ""))

    def test_random_block_whitening_and_fry(self):
        r = r_json(r'''
set.seed(9); G <- 600
subj <- rep(paste0("S", 1:6), 2); cond <- factor(rep(c("A", "B"), each = 6))
d <- model.matrix(~0 + cond); colnames(d) <- c("A", "B"); con <- c(A = -1, B = 1)
se <- rnorm(6, 0, 1); names(se) <- paste0("S", 1:6)
E <- matrix(rnorm(G * 12, 0, 0.3), G, 12) + rnorm(G, 10, 2) +
     matrix(se[subj], G, 12, byrow = TRUE)                   # a subject effect every protein shares
E[1:20, cond == "B"] <- E[1:20, cond == "B"] + 0.5
rho <- duplicateCorrelation(E, d, block = subj)$consensus
f <- lmFit(E, d, block = subj, correlation = rho); cs <- contrast_stats(f, con)
cor_of <- function(U) { m <- 20; v <- m * mean(colMeans(U[21:40, ])^2); (v - 1) / (m - 1) }
Ub <- residual_U(lm_effects(E, d, con, NULL, subj, rho)); Uu <- residual_U(lm_effects(E, d, con))
fb <- fry_fit(E, list(s = 1:20), d, con, NULL, subj, rho); fu <- fry_fit(E, list(s = 1:20), d, con)
cat(jsonlite::toJSON(list(rho = rho, cor_blocked = cor_of(Ub), cor_unblocked = cor_of(Uu),
                          fry_blocked = fb$PValue, fry_unblocked = fu$PValue), auto_unbox = TRUE, digits = NA))
''')
        self.assertGreater(r["rho"], 0.5)
        self.assertLess(r["cor_blocked"], r["cor_unblocked"])
        self.assertLess(r["fry_blocked"], r["fry_unblocked"])
        self.assertLess(r["fry_blocked"], 0.01)

    def test_global_correlation_depth_by_role_and_pooled_enrichment(self):
        r = r_json(r'''
set.seed(3); G <- 800; n <- 6
grp <- factor(rep(c("A", "B"), each = 3)); d <- model.matrix(~0 + grp); colnames(d) <- c("A", "B"); con <- c(A = -1, B = 1)
E0 <- matrix(rnorm(G * n), G) + rnorm(G, 20, 2)
E1 <- E0 + 0.8 * matrix(rnorm(n), G, n, byrow = TRUE)          # every protein on one run factor
g0 <- global_correlation(residual_U(lm_effects(E0, d, con)))
g1 <- global_correlation(residual_U(lm_effects(E1, d, con)))
dep <- c(1, 2, 3, 10, 11, 12); role <- c("bait", "control", "bait", "control", "bait", "control")
a <- add_depth(d, con, dep, role); b <- add_depth(d, con, dep)
cm <- cbind(x = c(1, 0, -1, 0), y = c(0, 1, 0, -1))
cat(jsonlite::toJSON(list(g0 = jsonlite::unbox(g0), g1 = jsonlite::unbox(g1), cols = a$columns, con = unname(a$contrast),
  bait_col = unname(a$design[, "run_depth_bait"]), ctrl_col = unname(a$design[, "run_depth_control"]),
  one = I(b$columns), pooled = unname(pooled_enrichment_contrast(cm, 1:2))), digits = NA))
''')
        self.assertLess(r["g0"], 0.02)
        self.assertGreater(r["g1"], 0.1)                  # camera's power problem, measured
        self.assertEqual(r["cols"], ["run_depth_bait", "run_depth_control"])
        self.assertEqual(r["con"], [-1, 1, 0, 0])
        self.assertEqual(r["bait_col"], [-4, 0, -2, 0, 6, 0])    # centred within bait runs (1, 3, 11)
        self.assertEqual(r["ctrl_col"], [0, -6, 0, 2, 0, 4])
        self.assertEqual(r["one"], ["run_depth"])
        self.assertEqual(r["pooled"], [0.5, 0.5, -0.5, -0.5])

    def test_weak_interactors_are_the_interactomes_bottom_quarter(self):
        # review M3: relative to the interactome, so a random set is rarely "mostly weak" -- the old
        # "within 1 log2 of the minimum" flagged nearly every set
        r = r_json(r'''
set.seed(5); e <- runif(400, 2, 5)
w <- weak_interactors(e)
fr <- replicate(2000, mean(w[sample(400, 20)]) > 0.5)
cat(jsonlite::toJSON(list(share = jsonlite::unbox(mean(w)), cut = jsonlite::unbox(max(e[w]) <= min(e[!w])),
  random_flagged = jsonlite::unbox(mean(fr)), weakest_flagged = jsonlite::unbox(mean(w[order(e)[1:20]]) > 0.5),
  old_rule = jsonlite::unbox(mean(replicate(2000, mean(e[sample(400, 20)] < 2 + 1) > 0.5))),
  bias2 = jsonlite::unbox(ip_direction_bias(2)), bias1 = jsonlite::unbox(ip_direction_bias(1)),
  bias0 = jsonlite::unbox(is.na(ip_direction_bias(0)))), digits = NA))
''')
        self.assertAlmostEqual(r["share"], 0.25, places=2)
        self.assertTrue(r["cut"])
        self.assertLess(r["random_flagged"], 0.02)
        self.assertTrue(r["weakest_flagged"])
        self.assertGreater(r["old_rule"], r["random_flagged"])
        # review M2/W2: what the enrichment minimum costs -- stated, never as an absolute
        self.assertTrue(r["bias2"].startswith("Losses from a complex are harder to detect than gains"))
        self.assertIn("4-fold", r["bias2"])
        self.assertNotIn("only gains", r["bias2"])
        self.assertIn("2-fold", r["bias1"])
        self.assertTrue(r["bias0"])

    def test_losses_per_bait(self):
        # review W2: each bait's own numbers -- the share of its interactome a fall would push below
        # the minimum, and the share of tested sets that would keep min_size members
        r = r_json(r'''
e <- c(rep(2.1, 4), rep(2.5, 4), rep(3.5, 12))      # 20 interactome proteins, minimum 2 log2
idx <- list(a = 1:20, b = 9:20, c = 1:10)
ls <- ip_loss_sensitivity(e, idx, 2, min_size = 10)
cat(jsonlite::toJSON(list(below = ls$below_minimum, kept = ls$sets_still_tested,
  txt = jsonlite::unbox(ip_losses_reading(ls, 2)), none = jsonlite::unbox(is.na(ip_losses_reading(ls, 0)))), digits = NA))
''')
        self.assertEqual(r["below"], [0.2, 0.4])          # 0.8 fall: 2.1s drop; 1.6 fall: 2.1s and 2.5s
        self.assertEqual(r["kept"], [0.667, 0.667])       # a and b keep >= 10; c does not
        self.assertEqual(r["txt"], "Losses are harder to detect than gains here: a 0.8 log2 fall in one condition "
                         "would put about 20% of this interactome below the 4-fold minimum, and about 67% of the "
                         "sets tested here would keep enough members to be tested (about 67% for a 1.6 log2 fall).")
        self.assertTrue(r["none"])

    def test_broad_shift_needs_the_proteins(self):
        # review W1: "most proteins moved together" only when >= 60% of proteins moved that way and
        # their median did; otherwise the sets are reported without it, with the bait's recovery
        r = r_json(r'''
base <- list(n_pos = 3, n_neg = 3, n_words = "samples, one per Mouse", df_residual = 20, camera_low_power = TRUE,
  camera_global_correlation = 0.13, camera_vif_50 = 7.4, between_block_caution = FALSE, camera = 0, both = 0,
  fry = 45, fry_up = 0, fry_down = 45, pos_groups = list("Old_JPH3"), depth_tested = TRUE, fry_depth_holds = 45,
  depth_difference_log2 = 0.1)
no <- modifyList(base, list(ip_role = "bait_between_conditions", protein_share_up = 0.47, protein_share_down = 0.53,
  protein_median_logFC = 0, recovery = list(bait = "JPH3", reference = "interactome", difference_log2 = -0.94)))
yes <- modifyList(base, list(ip_role = "control_vs_control", protein_share_up = 0.2, protein_share_down = 0.8,
  protein_median_logFC = -0.3))
cat(jsonlite::toJSON(list(no = jsonlite::unbox(contrast_reading(no, 0.05)),
  yes = jsonlite::unbox(contrast_reading(yes, 0.05))), digits = NA))
''')
        self.assertTrue(r["no"].startswith("45 of fry's 45 significant sets are lower in Old_JPH3 (FDR < 0.05; "), r["no"])
        self.assertIn("did not move that way (53% of the proteins lower, median log2 fold change 0.000)", r["no"])
        self.assertNotIn("broad shift", r["no"].lower())
        self.assertNotIn("most proteins moved together", r["no"])
        self.assertIn("Not corrected for bait recovery (the Old_JPH3 IPs brought down 0.94 log2 less of the JPH3 "
                      "complex, by its interactome reference)", r["no"])
        self.assertIn("camera has little power here", r["no"])
        self.assertTrue(r["yes"].startswith("A broad shift: 45 of fry's 45 significant sets are lower in Old_JPH3, as "
                                            "are 80% of the proteins (median log2 fold change -0.300) -- most proteins "
                                            "moved together; whether any category"), r["yes"])
        self.assertEqual(r["yes"].count("little power"), 1)
        self.assertEqual(r["yes"].count("most proteins"), 1)

    def test_a_value_that_rounds_to_zero_has_no_sign(self):
        r = r_json(r'''cat(jsonlite::toJSON(signed(c(0.0004, -0.0004, 0, -0.00049, 0.0006, -0.1234, 0.41), 3)))''')
        self.assertEqual(r, ["0.000", "0.000", "0.000", "0.000", "+0.001", "-0.123", "+0.410"])

    def test_reference_dependent_calls_give_direction_effect_and_recovery(self):
        # review: JPH3's 6 sets, -0.12 to -0.09 log2, none holding with the bait protein, whose recovery
        # estimate differs by 0.19 log2 -- within the reference uncertainty
        r = r_json(r'''
x <- list(bait = "JPH3", reference = "interactome", cross_check = "bait", n_pos = 3, n_neg = 3, df_residual = 4,
  camera = 0, fry = 6, both = 0, pos_group = "Old_JPH3", n_interactome = 212, n_prior_proteins = 212,
  camera_low_power = FALSE, recovery_difference = 0.19, losses = "LOSSES.",
  sig = data.frame(direction = "Down", effect = c(-0.116, -0.116, -0.115, -0.105, -0.094, -0.09),
                   depends = TRUE, weak = FALSE, flagged = TRUE))
y <- modifyList(x, list(recovery_difference = 0.05))
z <- x; z$sig <- data.frame(direction = c("Up", "Down"), effect = c(0.3, -0.2), depends = FALSE,
                            weak = c(TRUE, FALSE), flagged = c(TRUE, FALSE))
cat(jsonlite::toJSON(list(x = jsonlite::unbox(bait_reading(x)), y = jsonlite::unbox(bait_reading(y)),
  z = jsonlite::unbox(bait_reading(z))), digits = NA))
''')
        self.assertEqual(r["x"], "Relative to the JPH3 complex (3 vs 3 IPs, 4 residual df): camera 0, fry 6, both 0 "
                         "sets -- 6 lower in Old_JPH3 (-0.12 to -0.09 log2), relative to the interactome median. None "
                         "holds with the bait-protein reference, whose recovery estimate differs by 0.19 log2, so the "
                         "effect is within the reference uncertainty -- unresolved. LOSSES.")
        self.assertIn("so the call rests on the choice of reference -- unresolved", r["y"])
        self.assertIn("1 lower in Old_JPH3 (-0.20 log2) and 1 higher (+0.30 log2), relative to", r["z"])
        self.assertIn("All hold with the bait-protein reference. 1 is mostly weak interactors", r["z"])
        self.assertNotIn("Every one is flagged", r["z"])

    def test_presence_fraction_and_reindex(self):
        r = r_json(r'''
ev <- c("measured in both", "presence call", "presence call", "partly inferred", NA, "presence call")
idx <- list(a = 1:3, b = c(1, 4, 5), c = c(2, 3, 6))
m <- which(!(ev %in% PRESENCE_CALL_LABEL))
cat(jsonlite::toJSON(list(f = unname(presence_fraction(idx, ev)), re = unname(reindex(idx, m))), digits = NA))
''')
        self.assertAlmostEqual(r["f"][0], 2 / 3, places=6)
        self.assertEqual(r["f"][1], 0)
        self.assertEqual(r["f"][2], 1)
        self.assertEqual(r["re"], [[1], [1, 2, 3], []])


@unittest.skipUnless(r_has("limma", "jsonlite"), "needs R with limma and jsonlite")
class SetsInOut(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def gmt(self, name, lines):
        p = os.path.join(self.tmp, name)
        with open(p, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        return p

    def test_read_gmt_symbols_entrez_and_errors(self):
        sym = self.gmt("s.gmt", ["setA\tdesc A\tKcnb1\tVapa\tvapa", "setB\tdesc B\tRyr2"])
        ent = self.gmt("e.gmt", ["x\tdesc\t16535\t22099"])
        dup = self.gmt("d.gmt", ["x\td\tA", "x\td\tB"])
        short = self.gmt("t.gmt", ["x\tonly-two"])
        r = r_json(f'''
s <- read_gmt("{sym}"); e <- read_gmt("{ent}")
err <- function(p) tryCatch({{ read_gmt(p); "" }}, error = function(e) conditionMessage(e))
cat(jsonlite::toJSON(list(sets = s$sets, type = c(s$id_type, e$id_type), dup = err("{dup}"), short = err("{short}")),
                     auto_unbox = TRUE))
''')
        self.assertEqual(r["sets"]["setA"], ["KCNB1", "VAPA"])
        self.assertEqual(r["type"], ["gene symbol (case-insensitive)", "Entrez Gene ID"])
        self.assertIn("duplicated set names", r["dup"])
        self.assertIn("fewer than 3", r["short"])

    def test_size_limits_count_tested_proteins(self):
        r = r_json(r'''
coll <- list(sets = list(a = c("1", "2", "3"), b = c("1", "9"), c = as.character(1:6)), id = c("a", "b", "c"),
             name = c("A", "B", "C"), source = rep("x", 3))
ids <- list("1", c("2", "3"), "3", "4", "5", "6")
x <- index_sets(coll, ids, 2, 4)
cat(jsonlite::toJSON(list(kept = I(names(x$index)), a = x$index$a, counts = x$counts), auto_unbox = TRUE))
''')
        self.assertEqual(r["kept"], ["a"])            # b has 1 tested protein, c has 6
        self.assertEqual(r["a"], [1, 2, 3])
        self.assertEqual(r["counts"]["below_min"], 1)
        self.assertEqual(r["counts"]["above_max"], 1)

    def test_organism_from_the_fasta_sidecar_or_stop(self):
        meta = os.path.join(self.tmp, "search.fasta.meta.json")
        with open(meta, "w") as fh:
            json.dump({"organism": "Mus musculus", "taxid": 10090}, fh)
        r = r_json(f'''
o <- resolve_organism(NULL, "{meta}"); a <- resolve_organism("human")
err <- tryCatch({{ resolve_organism(NULL, NULL); "" }}, error = function(e) conditionMessage(e))
cat(jsonlite::toJSON(list(o = o$code, src = o$source, a = a$code, err = err), auto_unbox = TRUE))
''')
        self.assertEqual(r["o"], "Mm")
        self.assertIn("taxid 10090", r["src"])
        self.assertEqual(r["a"], "Hs")
        self.assertIn("--organism", r["err"])


# A pulldown-shaped DIA-NN report: bait (B) and control (C) IPs from 6 mice, 3 Old + 3 Young;
# each mouse gives one B and one C run. 80 interactors enriched in B; 15 of them (G1..G15) +1.2
# in Old B relative to the rest of the complex; Old recovers the bait 0.7 log2 units less.
SYNTH_IP_R = r'''
set.seed(7)
mice <- paste0("M", 1:6); age <- rep(c("Old", "Young"), each = 3)
runs <- data.frame(Run = c(paste0("b_", mice), paste0("c_", mice)),
                   Group = c(paste0(age, "_B"), paste0(age, "_C")), Mouse = c(mice, mice), stringsAsFactors = FALSE)
mouse_eff <- setNames(rnorm(6, 0, 0.3), mice)
rec <- setNames(ifelse(age == "Old", -0.7, 0) + rnorm(6, 0, 0.25), mice)
rows <- list(); n <- 0
for (i in 1:320) {
  pg <- sprintf("P%05d", i); gene <- sprintf("G%d", i); base <- runif(1, 12, 20); interactor <- i <= 80
  for (k in 1:3) {
    pb <- base + rnorm(1, 0, 0.8)
    for (r in seq_len(nrow(runs))) {
      m <- runs$Mouse[r]; isb <- startsWith(runs$Run[r], "b_"); old <- startsWith(runs$Group[r], "Old")
      lv <- pb + mouse_eff[[m]] + rnorm(1, 0, 0.25)
      if (interactor) lv <- lv + if (isb) 3 + rec[[m]] + (i <= 15 && old) * 1.2 else -1
      if (runif(1) < plogis(-(lv - 12.5) * 2)) next
      n <- n + 1
      rows[[n]] <- data.frame(Run = runs$Run[r], Precursor.Id = sprintf("PEP%dK%d2", i, k),
        Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0(gene, "_X"), Genes = gene,
        Proteotypic = 1L, Precursor.Normalised = 2^lv, Precursor.Quantity = 2^lv,
        Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001,
        Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
    }
  }
}
d <- do.call(rbind, rows)
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum); names(pgm)[3] <- "PG.MaxLFQ"
d <- merge(d, pgm, by = c("Protein.Group", "Run"))
arrow::write_parquet(d, file.path(OUT, "report.parquet"))
write.csv(data.frame(File.Name = runs$Run, Group = runs$Group, Mouse = runs$Mouse),
          file.path(OUT, "conditions.csv"), row.names = FALSE)
write.csv(data.frame(Group = c("Old_B", "Young_B", "Old_C", "Young_C"), Role = c("bait", "bait", "control", "control"),
                     Bait = c("B", "B", "", ""), Condition = c("Old", "Young", "Old", "Young"),
                     Bait_gene = c("G80", "G80", "", "")), file.path(OUT, "ip_map.csv"), row.names = FALSE)
sets <- c(list(planted = paste0("G", 1:15)),
          setNames(lapply(1:6, function(j) paste0("G", sample(16:79, 15))), paste0("int_null", 1:6)),
          setNames(lapply(1:20, function(j) paste0("G", sample(81:320, 20))), paste0("bg_null", 1:20)))
writeLines(vapply(names(sets), function(s) paste(c(s, "synthetic", sets[[s]]), collapse = "\t"), ""),
           file.path(OUT, "sets.gmt"))
'''
E2E_NEEDS = ("limpa", "limma", "arrow", "dplyr", "tidyr", "jsonlite")
CONTRASTS = "Old_B-Old_C,Young_B-Young_C,Old_B-Young_B,Old_C-Young_C,(Old_B+Old_C)/2-(Young_B+Young_C)/2"
POOLED = "(Old_B+Old_C)/2-(Young_B+Young_C)/2"
TABLE = {"Old_B-Old_C": "Old_B.Old_C", "Young_B-Young_C": "Young_B.Young_C", "Old_B-Young_B": "Old_B.Young_B",
         "Old_C-Young_C": "Old_C.Young_C", POOLED: "X.Old_B.Old_C..2..Young_B.Young_C..2"}


@unittest.skipUnless(r_has(*E2E_NEEDS), "needs R with " + "/".join(E2E_NEEDS))
class RunSetsEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_IP_R)      # the GMT is written BEFORE run_de
        # a session's layout: run_de into output/tables, run_sets' defaults beside it
        cls.session = os.path.join(cls.tmp, "sess")
        cls.de = os.path.join(cls.session, "output", "tables")
        cls.fig = os.path.join(cls.session, "output", "figures")
        os.makedirs(cls.fig)
        cls.de_proc = subprocess.run(
            ["Rscript", RUN_DE, "--input", "report.parquet", "--metadata", "conditions.csv",
             "--method", "dpc", "--contrasts", CONTRASTS, "--block", "Mouse", "--outdir", cls.de],
            capture_output=True, text=True, cwd=cls.tmp)
        cls.sets = cls.de                                  # run_sets' default --outdir
        cls.sets_proc = subprocess.run(
            ["Rscript", RUN_SETS, "--de-dir", cls.de, "--sets", "none", "--gmt", "sets.gmt",
             "--gmt-origin", "synthetic sets written before run_de", "--ip-map", "ip_map.csv",
             "--min-size", "5"],
            capture_output=True, text=True, cwd=cls.tmp)
        cls.prov = None
        if cls.sets_proc.returncode == 0:
            with open(os.path.join(cls.sets, "sets_provenance.json")) as fh:
                cls.prov = json.load(fh)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def setUp(self):
        self.assertEqual(self.de_proc.returncode, 0, self.de_proc.stderr[-2000:])
        self.assertEqual(self.sets_proc.returncode, 0, self.sets_proc.stderr[-2000:])

    def test_run_de_writes_inputs_and_suggests(self):
        self.assertTrue(os.path.exists(os.path.join(self.de, "set_test_inputs.rds")))
        with open(os.path.join(self.de, "de_provenance.json")) as fh:
            st = json.load(fh)["set_tests"]
        self.assertEqual(st["inputs"], "set_test_inputs.rds")
        self.assertEqual(st["suggest_below"], 10)
        self.assertIn("DEFAULT", st["suggest_below_source"])
        self.assertTrue({"Old_B-Young_B", "Old_C-Young_C"} <= set(st["suggested_for"]))
        with open(os.path.join(self.de, "de_provenance.json")) as fh:
            dp = json.load(fh)
        self.assertRegex(dp["first_written"], r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$")
        self.assertIn("protein-set test", self.de_proc.stderr)

    def test_the_model_is_run_des(self):
        self.assertLess(self.prov["model"]["refit_check"]["max_abs_t_difference"], 1e-6)
        self.assertEqual(self.prov["model"]["block"]["column"], "Mouse")
        self.assertEqual(self.prov["model"]["block"]["effect"], "random")
        roles = self.prov["pulldown"]["contrast_roles"]
        self.assertEqual(roles, {"Old_B-Old_C": "bait_vs_control", "Young_B-Young_C": "bait_vs_control",
                                 "Old_B-Young_B": "bait_between_conditions",
                                 "Old_C-Young_C": "control_vs_control", POOLED: "other"})
        c = self.prov["contrasts"]
        self.assertEqual(c["Old_B-Old_C"]["model"], "blocked")
        self.assertEqual(c["Old_B-Young_B"]["model"], "independent")
        self.assertFalse(c["Old_B-Old_C"]["depth_tested"])      # the depth difference IS the enrichment
        self.assertTrue(c["Old_C-Young_C"]["depth_tested"])
        # one depth slope per run type in a pulldown (review S3)
        self.assertEqual(self.prov["depth"]["columns"], ["run_depth_bait", "run_depth_control"])
        # a pooled between-mouse contrast comes from the blocked fit: it carries the CAUTION (S5)
        self.assertTrue(c[POOLED]["between_block_caution"])
        self.assertFalse(c["Old_B-Young_B"]["between_block_caution"])
        self.assertIn("anti-conservative", c[POOLED]["reading"])
        # "none" always carries n and residual df, and is never "unchanged" (B2)
        r = c["Old_C-Young_C"]["reading"]
        self.assertTrue(r.startswith("No set change detectable in this design (3 vs 3 samples, one per Mouse, "), r)
        self.assertIn("residual df", r)
        self.assertIn("not evidence that nothing changed", r)
        # camera's power is measured and said (S1)
        self.assertTrue(c["Old_B-Old_C"]["camera_low_power"])
        self.assertIn("camera has little power here", c["Old_B-Old_C"]["reading"])
        self.assertIn("What co-purifies with the bait", c["Old_B-Old_C"]["reading"])

    def test_tables_and_columns(self):
        for cn, fn in TABLE.items():
            rows = read_csv(os.path.join(self.sets, "Sets_%s.csv" % fn))
            self.assertEqual(len(rows), 27)
            for col in ("NGenes", "Mean_logFC", "Call", "camera_Direction", "camera_PValue", "camera_FDR",
                        "fry_Direction", "fry_PValue", "fry_FDR", "fry_Reference_Stable", "camera_Depth",
                        "fry_Depth", "depth_Mean_logFC", "depth_fry_FDR", "Presence_Call_Fraction",
                        "Mostly_Presence_Calls", "measured_camera_FDR", "measured_fry_FDR", "Flags"):
                self.assertIn(col, rows[0])
            for r in rows:     # every fry-significant set says whether its call held on every basis
                self.assertEqual(bool(r["fry_Reference_Stable"]), float(r["fry_FDR"]) < 0.05)
        # bait vs control: every interactome set is enriched (fry), none of the background ones
        rows = {r["Set_ID"]: r for r in read_csv(os.path.join(self.sets, "Sets_Old_B.Old_C.csv"))}
        self.assertLess(float(rows["planted"]["fry_FDR"]), 0.05)
        self.assertEqual(rows["planted"]["fry_Direction"], "Up")
        self.assertTrue(all(float(rows[f"bg_null{i}"]["fry_FDR"]) > 0.05 for i in range(1, 21)))

    def test_bait_normalised_finds_the_planted_set(self):
        b = self.prov["pulldown"]["baits"]["B"]
        self.assertEqual(b["n_interactome"], 80)
        self.assertIn("pooled over conditions", self.prov["pulldown"]["interactome_rule"])
        self.assertEqual(self.prov["pulldown"]["min_enrichment"]["value"], 2)
        self.assertIn("DEFAULT", self.prov["pulldown"]["min_enrichment"]["source"])
        self.assertEqual(b["bait_protein"], "P00080")
        rows = {r["Set_ID"]: r for r in read_csv(os.path.join(self.sets, "Sets_baitnorm_Old_B.Young_B.csv"))}
        self.assertEqual(set(rows), {"planted"} | {f"int_null{i}" for i in range(1, 7)})
        p = rows["planted"]
        self.assertEqual((p["Call"], p["camera_Direction"], p["fry_Direction"]), ("camera + fry", "Up", "Up"))
        self.assertIn("Weak_Fraction", p)
        self.assertIn("Relative to the B complex", b["contrasts"]["Old_B-Young_B"]["reading"])
        self.assertIn("variance prior rests on only 80 proteins", b["contrasts"]["Old_B-Young_B"]["reading"])
        # what a fall in one condition would cost, with this bait's own numbers (M2, W2)
        self.assertIn("Losses from a complex are harder to detect than gains", self.prov["pulldown"]["direction_bias"])
        self.assertEqual(b["losses"]["fall_log2"], [0.8, 1.6])
        self.assertTrue(b["losses_reading"].startswith("Losses are harder to detect than gains here: a 0.8 log2 fall"))
        self.assertIn(b["losses_reading"], b["contrasts"]["Old_B-Young_B"]["reading"])
        self.assertIn("relative to the interactome median", b["contrasts"]["Old_B-Young_B"]["reading"])
        self.assertIn("bottom 25%", self.prov["pulldown"]["weak_rule"])
        self.assertTrue(0 <= float(p["Weak_Fraction"]) <= 1)
        self.assertEqual(p["Reference_Robust"], "holds with the other reference")
        # the planted 15 of 80 pull the interactome median up: fry calls the rest "down" relative to
        # it, the bait-protein reference does not -- and the table says so
        for i in range(1, 7):
            r = rows[f"int_null{i}"]
            if r["fry_FDR"] and float(r["fry_FDR"]) < 0.05:
                self.assertEqual(r["Reference_Robust"], "depends on the reference")
                self.assertIn("depends on the reference", r["Flags"])

    def test_record_methods_and_defaults(self):
        s = self.prov["settings"]
        self.assertEqual(s["min_size"], {"value": 5, "source": "user"})
        self.assertIn("DEFAULT", s["max_size"]["source"])
        self.assertIn("DEFAULT", s["camera_inter_gene_cor"]["source"])
        self.assertEqual(s["fdr"]["source"], "run_de --adjp")
        self.assertTrue(self.prov["sources"]["gmt"]["predates_results"])
        self.assertEqual(len(self.prov["sources"]["gmt"]["sha256"]), 64)
        m = self.prov["methods_paragraph"]
        for words in ("camera", "fry", "Mouse as a random effect", "block-whitened", "Benjamini-Hochberg",
                      "presence calls", "interactome", "synthetic sets written before run_de",
                      "no correction was made across the 5 comparisons", "fewer than 10 significant proteins",
                      "pooled over conditions", "one slope for each run type", "re-tested on 4 other reference bases",
                      "losses relative to the complex are harder to detect than gains: a 0.8 log2 fall",
                      "slightly optimistic",
                      "least-enriched quarter",
                      "a user-made list (DEFAULT"):
            self.assertIn(words, m)
        with open(os.path.join(self.sets, "sets_methods.txt")) as fh:
            self.assertIn("Wu D, Smyth GK (2012)", fh.read())
        self.assertEqual(self.prov["figure_dir"], "../figures")
        self.assertIn("Members", self.prov["columns"])
        members = {r["Set_ID"]: r for r in read_csv(os.path.join(self.sets, self.prov["members"]))}
        self.assertEqual(members["planted"]["Members"], ";".join(sorted(f"G{i}" for i in range(1, 16))))
        self.assertTrue(any("Nothing significant" in r for r in self.prov["reading"]))

    def test_the_report_page_has_the_section(self):
        out = os.path.join(self.session, "output", "Analysis_Report.html")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--tables", self.de, "--figures", self.fig, "--out", out, "--no-pdf"],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.splitext(out)[0] + ".md", encoding="utf-8") as fh:
            md = fh.read()
        sec = md[md.index("## Protein-set tests"):]
        sec = sec[:sec.index("\n## ", 5)] if "\n## " in sec[5:] else sec
        self.assertIn("| Comparison | Significant proteins | Sets tested | camera | fry | What it shows |", sec)
        self.assertIn("(little power)", sec)
        self.assertIn(self.prov["contrasts"]["Old_C-Young_C"]["reading"], sec)
        self.assertNotIn("no set was significant in both", sec)
        self.assertIn("control vs control: the lysate background", sec)
        self.assertIn("(DEFAULT — not user-confirmed)", sec)
        self.assertIn("| B | Old B vs Young B | 80 |", sec)
        self.assertIn("holds with the other reference", sec)      # the planted set's detail row
        bias = self.prov["pulldown"]["direction_bias"]           # once, above the table (M2)
        self.assertEqual(sec.count(bias), 1)
        self.assertIn("**What these tests see less well.** " + bias, sec)
        self.assertIn(self.prov["pulldown"]["baits"]["B"]["contrasts"]["Old_B-Young_B"]["reading"], sec)
        self.assertIn("Nothing significant means", sec)            # the record's reading rules
        self.assertIn("*Methods:* Protein-set tests used", sec)
        # the figures of the comparisons run_de suggested set tests for, embedded; none "stale"
        for cn in ("Old_B.Young_B", "Old_C.Young_C"):
            self.assertIn(f"sets_{cn}.png", sec)
        self.assertNotIn("not in this run's figures.json", md)
        self.assertNotIn("sets_", p.stderr)

    def test_brief_agents_readme_and_methods_quote_the_record(self):
        brief = os.path.join(self.tmp, "brief.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                            "--de-dir", self.de, "--out", brief], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(brief, encoding="utf-8") as fh:
            b = fh.read()
        self.assertIn("### Protein-set tests", b)
        self.assertIn("start from the record's own reading of each comparison", b)
        for cn, c in self.prov["contrasts"].items():
            self.assertIn(c["reading"], b)
        self.assertIn("`Reference_Robust`:", b)
        self.assertIn("Never write that the complex lost nothing.", b)
        for rule in self.prov["reading"]:
            self.assertIn(rule, b)
        sys.path.insert(0, SCRIPTS)
        import session_docs
        f = session_docs.gather(self.session)
        agents = session_docs.agents_md(f)
        self.assertIn("| Protein-set tests (camera + fry) per contrast, and how they were run | "
                      "`output/tables/Sets_<contrast>.csv` (6 files) + "
                      "`output/tables/sets_provenance.json` |", agents)
        self.assertIn("- `Presence_Call_Fraction` — ", agents)
        self.assertIn("read them by run_sets.R's rules", agents)
        self.assertIn("protein-set tests `Sets_*.csv`", session_docs.readme_md(f))
        raw = os.path.join(self.tmp, "x.raw")
        with open(raw, "wb") as fh:
            fh.write(b"\0" * 64)
        out = os.path.join(self.tmp, "methods.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"), "--raw", raw,
                            "--out", out, "--de-dir", self.de], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(out, encoding="utf-8") as fh:
            m = fh.read()
        self.assertIn("## Protein-set tests\n\n" + self.prov["methods_paragraph"], m)

    def test_gmt_needs_an_origin_and_a_late_gmt_is_flagged(self):
        p = subprocess.run(["Rscript", RUN_SETS, "--de-dir", self.de, "--sets", "none", "--gmt", "sets.gmt",
                            "--outdir", os.path.join(self.tmp, "x")], capture_output=True, text=True, cwd=self.tmp)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("--gmt-origin", p.stderr)
        late = os.path.join(self.tmp, "late.gmt")
        shutil.copy(os.path.join(self.tmp, "sets.gmt"), late)
        future = time.time() + 5
        os.utime(late, (future, future))
        out = os.path.join(self.tmp, "late")
        p = subprocess.run(["Rscript", RUN_SETS, "--de-dir", self.de, "--sets", "none", "--gmt", late,
                            "--gmt-origin", "made after the results", "--contrasts", "Old_B-Young_B",
                            "--outdir", out, "--figdir", out], capture_output=True, text=True, cwd=self.tmp)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(out, "sets_provenance.json")) as fh:
            g = json.load(fh)["sources"]["gmt"]
        self.assertFalse(g["predates_results"])
        self.assertIn("AFTER the first DE results", p.stderr)
        # a public collection, named with its version: recorded as such, not warned about
        out2 = os.path.join(self.tmp, "public")
        p = subprocess.run(["Rscript", RUN_SETS, "--de-dir", self.de, "--sets", "none", "--gmt", late,
                            "--gmt-origin", "MSigDB Hallmark v2024.1", "--gmt-kind", "public",
                            "--contrasts", "Old_B-Young_B", "--outdir", out2, "--figdir", out2],
                           capture_output=True, text=True, cwd=self.tmp)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertNotIn("AFTER the first DE results", p.stderr)
        with open(os.path.join(out2, "sets_provenance.json")) as fh:
            pr = json.load(fh)
        self.assertEqual(pr["sources"]["gmt"]["kind"], {"value": "public", "source": "user"})
        self.assertIn("a public collection", pr["methods_paragraph"])

    def test_without_inputs_it_says_rerun_run_de(self):
        empty = os.path.join(self.tmp, "empty")
        os.makedirs(empty, exist_ok=True)
        p = subprocess.run(["Rscript", RUN_SETS, "--de-dir", empty, "--gmt", "sets.gmt", "--gmt-origin", "x"],
                           capture_output=True, text=True, cwd=self.tmp)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("re-run run_de.R", p.stderr)


@unittest.skipUnless(r_has(*E2E_NEEDS), "needs R with " + "/".join(E2E_NEEDS))
class RunSetsPairedAndMaxlfq(unittest.TestCase):
    """The other two model routes: a crossed Subject block (auto -> FIXED effect: block columns
    in the design, no correlation) on dpc, and maxlfq (a matrix with missing values: the tests
    run on complete rows). test_run_de_contaminants' report: 80 proteins, A x3 / B x3, G1-G5
    up 1.5 in B."""

    @classmethod
    def setUpClass(cls):
        from test_run_de_contaminants import SYNTH_R
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        with open(os.path.join(cls.tmp, "paired.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["File.Name", "Group", "Subject"])
            for i, g in enumerate(["A"] * 3 + ["B"] * 3):
                w.writerow([f"run{i + 1:02d}", g, f"S{i % 3 + 1}"])
        with open(os.path.join(cls.tmp, "sets.gmt"), "w") as fh:
            fh.write("up\tplanted\t" + "\t".join(f"G{i}" for i in range(1, 6)) + "\n")
            for j in range(6):
                fh.write(f"null{j}\tnull\t" + "\t".join(f"G{i}" for i in range(10 + 10 * j, 18 + 10 * j)) + "\n")
        cls.runs = {}
        for key, args in (("fixed", ["--method", "dpc", "--metadata", "paired.csv", "--block", "Subject"]),
                          ("maxlfq", ["--method", "maxlfq", "--metadata", "conditions.csv"])):
            out = os.path.join(cls.tmp, key)
            de = subprocess.run(["Rscript", RUN_DE, "--input", "report.parquet", "--outdir", out, *args],
                                capture_output=True, text=True, cwd=cls.tmp)
            rs = subprocess.run(["Rscript", RUN_SETS, "--de-dir", out, "--sets", "none", "--gmt", "sets.gmt",
                                 "--gmt-origin", "synthetic", "--min-size", "5"],
                                capture_output=True, text=True, cwd=cls.tmp)
            cls.runs[key] = (de, rs, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def check(self, key):
        de, rs, out = self.runs[key]
        self.assertEqual(de.returncode, 0, de.stderr[-1500:])
        self.assertEqual(rs.returncode, 0, rs.stderr[-1500:])
        with open(os.path.join(out, "sets_provenance.json")) as fh:
            prov = json.load(fh)
        rows = {r["Set_ID"]: r for r in read_csv(os.path.join(out, "Sets_B.A.csv"))}
        return prov, rows

    def test_fixed_block(self):
        prov, rows = self.check("fixed")
        self.assertEqual(prov["model"]["block"]["effect"], "fixed")
        self.assertLess(prov["model"]["refit_check"]["max_abs_t_difference"], 1e-6)
        self.assertIn("Subject as a fixed effect", prov["methods_paragraph"])
        self.assertNotIn("random effect", prov["methods_paragraph"])
        self.assertEqual((rows["up"]["fry_Direction"], rows["up"]["camera_Direction"]), ("Up", "Up"))
        self.assertLess(float(rows["up"]["fry_FDR"]), 0.05)

    def test_first_written_survives_a_rerun(self):
        de, rs, out = self.runs["maxlfq"]
        self.assertEqual(de.returncode, 0, de.stderr[-1500:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            first = json.load(fh)
        time.sleep(1.1)
        p = subprocess.run(["Rscript", RUN_DE, "--input", "report.parquet", "--outdir", out,
                            "--method", "maxlfq", "--metadata", "conditions.csv"],
                           capture_output=True, text=True, cwd=self.tmp)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            again = json.load(fh)
        self.assertEqual(again["first_written"], first["first_written"])
        self.assertNotEqual(again["written"], first["written"])

    def test_maxlfq(self):
        prov, rows = self.check("maxlfq")
        self.assertEqual(prov["method"], "maxlfq")
        self.assertLess(prov["model"]["refit_check"]["max_abs_t_difference"], 1e-6)
        t = prov["tested_proteins"]
        self.assertLessEqual(t["n"], t["of"])
        self.assertNotIn("precision weights", prov["methods_paragraph"])
        self.assertEqual(rows["up"]["fry_Direction"], "Up")


# setup.sh's GO-annotation step, with test_setup_limpa's stubs (no conda, no network): this
# Rscript stub reports GO.db and org.Mm.eg.db missing until a Bioconductor annotation-repo
# install.packages() "installs" them; the stub micromamba's install changes nothing, so the
# conda-first attempt is followed by the Bioconductor fallback.
FAKE_RSCRIPT_ANN = r'''#!/bin/bash
echo "$*" >> "$FAKE_R_LOG"
case "$*" in
  *"cat(p[!vapply"*) [ -f "$FAKE_ANN_STATE" ] || printf 'GO.db org.Mm.eg.db'; exit 0 ;;
  *BioCann*) [ -n "${FAKE_ANN_INSTALL_FAIL:-}" ] && { echo "cannot open URL" >&2; exit 1; }
             : > "$FAKE_ANN_STATE"; exit 0 ;;     # no touch on the stub PATH
  *"cat(tryCatch"*) printf '1.4.0'; exit 0 ;;
esac
exit 0
'''


class SetupInstallsGoAnnotation(setup_limpa.setup_t.SetupCheckHarness):
    run_setup = setup_limpa.SetupLimpaGate.run_setup

    def setUp(self):            # SetupLimpaGate.setUp, whose super() needs that class
        super().setUp()
        for tool in ("chmod", "ln"):
            os.symlink(shutil.which(tool), os.path.join(self.sys, tool))
        setup_limpa._exe(os.path.join(self.sys, "micromamba"), setup_limpa.FAKE_MICROMAMBA)
        self.ensure = setup_limpa._exe(os.path.join(self.d, "ensure_dotnet8_stub.sh"), setup_limpa.FAKE_ENSURE)
        self.dotnet_root = setup_limpa.fake_dotnet_root(os.path.join(self.d, ".proteomics-pipeline", "dotnet8"))
        self.mm_log = os.path.join(self.d, "mm.log")
        self.r_log = os.path.join(self.d, "r.log")
        self.state = os.path.join(self.d, "limpa_version")
        self.prefix = os.path.join(self.home, "micromamba", "envs", "proteomics-pipeline")

    def prebuilt_env(self):
        setup_limpa.SetupLimpaGate.prebuilt_env(self, "1.4.0")
        rs = os.path.join(self.prefix, "bin", "Rscript")
        os.remove(rs)
        setup_limpa._exe(rs, FAKE_RSCRIPT_ANN)
        self.ann_state = os.path.join(self.d, "ann_installed")

    def test_missing_packages_conda_first_then_bioconductor(self):
        self.prebuilt_env()
        r, s = self.run_setup(FAKE_ANN_STATE=self.ann_state)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        mm = [c for c in setup_limpa._lines(self.mm_log) if "bioconductor-go.db" in c]
        self.assertEqual(len(mm), 1, setup_limpa._lines(self.mm_log))
        self.assertIn("--override-channels -c conda-forge -c bioconda", mm[0])
        for pkg in ("bioconductor-annotationdbi", "bioconductor-org.hs.eg.db", "bioconductor-org.mm.eg.db"):
            self.assertIn(pkg, mm[0])
        ann = [c for c in setup_limpa._lines(self.r_log) if "BioCann" in c]
        self.assertEqual(len(ann), 1)
        self.assertIn("strsplit('GO.db org.Mm.eg.db', ' ')", ann[0])      # only what is missing
        self.assertIn("https://bioconductor.org/packages/3.23/data/annotation", ann[0])
        self.assertEqual(s["set_annotation"]["ok"], True)
        self.assertEqual(s["set_annotation"]["missing"], "")

    def test_a_failed_install_is_a_note_not_a_failure(self):
        self.prebuilt_env()
        r, s = self.run_setup(FAKE_ANN_STATE=self.ann_state, FAKE_ANN_INSTALL_FAIL="1")
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertEqual(s["set_annotation"]["ok"], False)
        self.assertEqual(s["set_annotation"]["missing"], "GO.db org.Mm.eg.db")
        note = [n for n in s["notes"] if "run_sets.R" in n]
        self.assertEqual(len(note), 1, s["notes"])
        self.assertIn("c('GO.db','org.Mm.eg.db')", note[0])
        self.assertIn("DE does not need them", note[0])

    def test_check_only_installs_nothing(self):
        self.prebuilt_env()
        r, s = self.run_setup("--check", FAKE_ANN_STATE=self.ann_state)
        self.assertFalse(any("BioCann" in c for c in setup_limpa._lines(self.r_log)))
        self.assertFalse(any("bioconductor-go.db" in c for c in setup_limpa._lines(self.mm_log)))
        self.assertEqual(s["set_annotation"]["ok"], False)


@unittest.skipUnless(r_has("limma"), "needs R with limma")
class SimulationRuns(unittest.TestCase):
    def test_both_simulations_run(self):
        for what in ("camera", "ip"):
            p = subprocess.run(["Rscript", SIM, "--what", what, "--reps", "2"], capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr[-1500:])
            self.assertIn("rate of p < 0.05", p.stdout)
            if what == "ip":         # run_sets.R's minimum size and weak flag, applied (review M3)
                self.assertIn("tested: kept >= 10", p.stdout)
                self.assertIn("null_random fry false-positive rate at the default", p.stdout)


if __name__ == "__main__":
    unittest.main()
