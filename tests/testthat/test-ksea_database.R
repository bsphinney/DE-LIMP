# The KSEA kinase-substrate database (issue #3).
#
# KSEAapp bundles KSData: 9,757 PhosphoSitePlus records plus a 10,000-row slice
# of NetworKIN. The slice is NOT a quality filter — it is neither the top 10,000
# by score nor the first 10,000 in file order, and its score distribution is
# indistinguishable from the full set. So which substrates a kinase is credited
# with was arbitrary.
#
# It matters more than the row count suggests: KSEA.Scores() runs with
# NetworKIN.cutoff = 5, and only 759 of the sampled rows clear it against 17,896
# in the complete table.

test_that("every annotation the analysis can use is present", {
  skip_if_not_installed("KSEAapp")
  db <- load_ksea_database()
  expect_equal(sum(db$Source == "PhosphoSitePlus"), 9757)   # all of them
  expect_equal(sum(db$Source == "NetworKIN"), 17896)        # all that clear cutoff = 5
  expect_equal(nrow(db), 27653)
})

test_that("nothing usable was dropped when the table was trimmed for size", {
  # The shipped NetworKIN file is the >= 5 subset of the 230,992-row table,
  # because KSEA.Scores(NetworKIN.cutoff = 5) discards the rest before scoring.
  # This asserts the trim boundary: no row below the cutoff, none above missing.
  f <- system.file_or_local("data/networkin_kinase_substrate_ge5_July2016.csv")
  skip_if(!nzchar(f), "NetworKIN table not present")
  nk <- read.csv(f, stringsAsFactors = FALSE)
  sc <- suppressWarnings(as.numeric(nk$networkin_score))
  expect_true(all(sc >= 5, na.rm = TRUE))
  expect_equal(nrow(nk), 17896)
})

test_that("it is the annotations that survive the cutoff that matter", {
  skip_if_not_installed("KSEAapp")
  db <- load_ksea_database()
  usable <- sum(db$Source == "NetworKIN" &
                  suppressWarnings(as.numeric(db$networkin_score)) >= 5, na.rm = TRUE)
  expect_equal(usable, 17896)                          # was 759 with the sampled slice

  data("KSData", package = "KSEAapp", envir = environment())
  ks <- get("KSData", envir = environment())
  before <- sum(ks$Source == "NetworKIN" &
                  suppressWarnings(as.numeric(ks$networkin_score)) >= 5, na.rm = TRUE)
  expect_gt(usable / before, 20)
})

test_that("DE-LIMP ships no PhosphoSitePlus data of its own", {
  # PSP is licensed by Cell Signaling Technology under a formal licence
  # agreement (phosphosite.org/staticLicensing, checked 2026-08-27). The curated
  # half must keep arriving via the KSEAapp dependency, never from this repo.
  f <- system.file_or_local("data/networkin_kinase_substrate_ge5_July2016.csv")
  skip_if(!nzchar(f), "NetworKIN table not present")
  nk <- read.csv(f, stringsAsFactors = FALSE)
  expect_true(all(nk$Source == "NetworKIN"))
  expect_false(any(nk$Source == "PhosphoSitePlus"))
})

test_that("the schema still matches what KSEA.Scores expects", {
  skip_if_not_installed("KSEAapp")
  data("KSData", package = "KSEAapp", envir = environment())
  ks <- get("KSData", envir = environment())
  expect_equal(names(load_ksea_database()), names(ks))
})

test_that("what was actually used is recorded, never silently substituted", {
  skip_if_not_installed("KSEAapp")
  # Architectural rule #2: the run log quotes this, so it has to name the real
  # table — including when the full one is missing and the sample was used.
  expect_match(attr(load_ksea_database(), "ksea_source"), "NetworKIN 17896")
  expect_match(attr(load_ksea_database(path = "/no/such/file.csv.gz"), "ksea_source"),
               "not found")
})

test_that("the file is found from any working directory the app runs in", {
  # runApp() from the repo root, the container from /app, testthat from
  # tests/testthat. A getwd()-only lookup silently degrades in two of the three.
  expect_true(nzchar(system.file_or_local("data/networkin_kinase_substrate_ge5_July2016.csv")))
})
