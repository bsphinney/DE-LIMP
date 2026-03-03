# Data integrity test for the example dataset pipeline
#
# These tests verify that the example HeLa dataset produces consistent
# outputs through the LIMPA pipeline. They require the limpa package
# and a local copy of the parquet file.
#
# Run manually: testthat::test_file("tests/testthat/test-data_integrity.R")

test_that("example dataset loads and produces consistent dimensions", {
  skip_if_not_installed("limpa")
  skip_if_not_installed("arrow")

  local_path <- file.path(project_root, "Affinisep_vs_evosep_noNorm.parquet")
  skip_if(!file.exists(local_path), "Example parquet not found locally")

  dat <- limpa::readDIANN(local_path, format = "parquet", q.cutoffs = 0.01)

  # EList from limpa (S4 class)
  expect_true(is(dat, "EList"))
  expect_true("E" %in% names(dat))

  # 12 samples in the HeLa dataset
  expect_equal(ncol(dat$E), 12)
  # Precursor-level rows (40k-50k range at q=0.01)
  expect_true(nrow(dat$E) >= 30000)
  expect_true(nrow(dat$E) <= 60000)
})

test_that("example dataset sample names are consistent", {
  skip_if_not_installed("limpa")
  skip_if_not_installed("arrow")

  local_path <- file.path(project_root, "Affinisep_vs_evosep_noNorm.parquet")
  skip_if(!file.exists(local_path), "Example parquet not found locally")

  dat <- limpa::readDIANN(local_path, format = "parquet", q.cutoffs = 0.01)
  sample_names <- sort(colnames(dat$E))

  # Verify both prep methods are present
  expect_true(any(grepl("affinisep", sample_names, ignore.case = TRUE)))
  expect_equal(length(sample_names), 12)
})

test_that("example dataset protein IDs are detected as human", {
  skip_if_not_installed("limpa")
  skip_if_not_installed("arrow")

  local_path <- file.path(project_root, "Affinisep_vs_evosep_noNorm.parquet")
  skip_if(!file.exists(local_path), "Example parquet not found locally")

  dat <- limpa::readDIANN(local_path, format = "parquet", q.cutoffs = 0.01)
  protein_ids <- rownames(dat$E)

  org_db <- detect_organism_db(protein_ids)
  expect_equal(org_db, "org.Hs.eg.db")
})

test_that("QC stats computation produces expected structure", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("dplyr")

  local_path <- file.path(project_root, "Affinisep_vs_evosep_noNorm.parquet")
  skip_if(!file.exists(local_path), "Example parquet not found locally")

  # Need dplyr loaded for get_diann_stats_r
  library(dplyr)
  qc <- get_diann_stats_r(local_path)

  expect_s3_class(qc, "data.frame")
  expect_true(nrow(qc) > 0)
  # Expected columns
  expect_true("Run" %in% names(qc))
  expect_true("Precursors" %in% names(qc))
  expect_true("Proteins" %in% names(qc))
  # Should have 12 runs
  expect_equal(nrow(qc), 12)
  # Protein counts should be reasonable (5000-10000 per run)
  expect_true(all(qc$Proteins > 3000))
})
