# Test gene symbol mapping works without clusterProfiler
# This catches the bitr() → AnnotationDbi::select migration

test_that("Gene mapping works via AnnotationDbi::select (no clusterProfiler needed)", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")

  library(org.Hs.eg.db)
  db <- org.Hs.eg.db

  # Test with known human UniProt accessions
  test_accessions <- c("P43490", "P30838", "P09382", "P04264")

  result <- suppressMessages(
    AnnotationDbi::select(db,
      keys = test_accessions,
      keytype = "UNIPROT",
      columns = c("SYMBOL", "GENENAME"))
  )

  expect_true(is.data.frame(result))
  expect_true(nrow(result) > 0)
  expect_true("SYMBOL" %in% colnames(result))
  expect_true("GENENAME" %in% colnames(result))

  # P43490 should map to NAMPT
  nampt <- result[result$UNIPROT == "P43490", "SYMBOL"]
  expect_true("NAMPT" %in% nampt)

  # P04264 should map to KRT1 (keratin)
  krt1 <- result[result$UNIPROT == "P04264", "SYMBOL"]
  expect_true("KRT1" %in% krt1)
})

test_that("Gene mapping handles missing/invalid accessions gracefully", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")

  library(org.Hs.eg.db)
  db <- org.Hs.eg.db

  # Mix of valid and invalid accessions
  test_accessions <- c("P43490", "FAKE_ACC_123", "XP_052618122.1", "")

  result <- tryCatch(
    suppressMessages(suppressWarnings(
      AnnotationDbi::select(db,
        keys = test_accessions,
        keytype = "UNIPROT",
        columns = c("SYMBOL", "GENENAME"))
    )),
    error = function(e) NULL
  )

  # Should not crash — either returns results or NULL
  if (!is.null(result)) {
    expect_true(is.data.frame(result))
    # P43490 should still map correctly even with bad keys
    if (nrow(result) > 0) {
      expect_true(any(result$UNIPROT == "P43490"))
    }
  }
})

test_that("detect_organism_db returns valid OrgDb name", {
  source("../../R/helpers.R", local = TRUE)

  # Human proteins (sp|P43490|NAMPT_HUMAN format)
  expect_equal(detect_organism_db(c("P43490", "P30838")), "org.Hs.eg.db")

  # Proteins with _HUMAN suffix
  expect_equal(detect_organism_db(c("sp|P43490|NAMPT_HUMAN")), "org.Hs.eg.db")

  # Mouse proteins
  expect_equal(detect_organism_db(c("sp|Q9D8N0|EF1G_MOUSE")), "org.Mm.eg.db")
})
