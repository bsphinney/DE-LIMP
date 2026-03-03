# Tests for R/helpers.R — pure utility functions

# =============================================================================
# cal_z_score
# =============================================================================

test_that("cal_z_score returns z-scores", {
  x <- c(1, 2, 3, 4, 5)
  z <- cal_z_score(x)
  expect_equal(mean(z), 0, tolerance = 1e-10)
  expect_equal(sd(z), 1, tolerance = 1e-10)
  expect_length(z, 5)
})

test_that("cal_z_score handles NAs", {
  x <- c(1, 2, NA, 4, 5)
  z <- cal_z_score(x)
  expect_true(is.na(z[3]))
  expect_equal(mean(z, na.rm = TRUE), 0, tolerance = 1e-10)
})

test_that("cal_z_score handles constant vector", {
  x <- c(5, 5, 5)
  z <- cal_z_score(x)
  # sd = 0 → division by zero → NaN

expect_true(all(is.nan(z)))
})

test_that("cal_z_score handles single value", {
  z <- cal_z_score(42)
  # sd(42) is NA → division produces NaN or NA
  expect_true(is.na(z))
})

# =============================================================================
# detect_organism_db
# =============================================================================

test_that("detect_organism_db identifies human proteins", {
  ids <- c("sp|P04406|G3P_HUMAN", "sp|P68871|HBB_HUMAN")
  expect_equal(detect_organism_db(ids), "org.Hs.eg.db")
})

test_that("detect_organism_db identifies mouse proteins", {
  ids <- c("sp|P10126|EF1A1_MOUSE", "sp|P16858|G3P_MOUSE")
  expect_equal(detect_organism_db(ids), "org.Mm.eg.db")
})

test_that("detect_organism_db identifies yeast proteins", {
  ids <- c("sp|P00560|PGK_YEAST")
  expect_equal(detect_organism_db(ids), "org.Sc.sgd.db")
})

test_that("detect_organism_db defaults to human when unknown", {
  ids <- c("PROT001", "PROT002")
  expect_equal(detect_organism_db(ids), "org.Hs.eg.db")
})

test_that("detect_organism_db handles mixed organisms (first wins)", {
  ids <- c("sp|P04406|G3P_HUMAN", "sp|P16858|G3P_MOUSE")
  # Should match HUMAN first (it's checked first in the map)
  expect_equal(detect_organism_db(ids), "org.Hs.eg.db")
})

test_that("detect_organism_db is case-insensitive", {
  ids <- c("sp|P04406|G3P_human")
  expect_equal(detect_organism_db(ids), "org.Hs.eg.db")
})
