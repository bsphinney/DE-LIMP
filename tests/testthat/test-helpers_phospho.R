# Tests for R/helpers_phospho.R — phosphoproteomics helper functions

# =============================================================================
# parse_phospho_positions
# =============================================================================

test_that("parse_phospho_positions finds single phosphosite", {
  # S at position 5, modified with UniMod:21
  result <- parse_phospho_positions(
    mod_seq = "AAGAS(UniMod:21)PEPTK",
    stripped_seq = "AAGASPEPTK"
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$site_residue, "S")
  expect_equal(result$site_peptide_pos, 5L)
})

test_that("parse_phospho_positions finds multiple phosphosites", {
  result <- parse_phospho_positions(
    mod_seq = "AAS(UniMod:21)GT(UniMod:21)PK",
    stripped_seq = "AASGTPK"
  )
  expect_equal(nrow(result), 2)
  expect_equal(result$site_residue, c("S", "T"))
  expect_equal(result$site_peptide_pos, c(3L, 5L))
})

test_that("parse_phospho_positions returns empty for no phospho", {
  result <- parse_phospho_positions(
    mod_seq = "AAGASPEPTK",
    stripped_seq = "AAGASPEPTK"
  )
  expect_equal(nrow(result), 0)
  expect_equal(names(result), c("site_residue", "site_peptide_pos"))
})

test_that("parse_phospho_positions ignores non-phospho modifications", {
  # Met oxidation (UniMod:35) should be ignored
  result <- parse_phospho_positions(
    mod_seq = "AAM(UniMod:35)S(UniMod:21)PK",
    stripped_seq = "AAMSPK"
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$site_residue, "S")
  expect_equal(result$site_peptide_pos, 4L)
})

test_that("parse_phospho_positions handles Y phosphorylation", {
  result <- parse_phospho_positions(
    mod_seq = "AAGAY(UniMod:21)PEPTK",
    stripped_seq = "AAGAYPEPTK"
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$site_residue, "Y")
  expect_equal(result$site_peptide_pos, 5L)
})

test_that("parse_phospho_positions handles T phosphorylation", {
  result <- parse_phospho_positions(
    mod_seq = "AT(UniMod:21)PK",
    stripped_seq = "ATPK"
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$site_residue, "T")
  expect_equal(result$site_peptide_pos, 2L)
})

test_that("parse_phospho_positions handles phospho at first position", {
  result <- parse_phospho_positions(
    mod_seq = "S(UniMod:21)PEPTK",
    stripped_seq = "SPEPTK"
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$site_residue, "S")
  expect_equal(result$site_peptide_pos, 1L)
})

test_that("parse_phospho_positions handles phospho at last position", {
  result <- parse_phospho_positions(
    mod_seq = "PEPTS(UniMod:21)",
    stripped_seq = "PEPTS"
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$site_residue, "S")
  expect_equal(result$site_peptide_pos, 5L)
})

test_that("parse_phospho_positions handles mixed mods with phospho", {
  # N-term acetyl + Met ox + phospho
  result <- parse_phospho_positions(
    mod_seq = "(UniMod:1)AM(UniMod:35)GS(UniMod:21)PK",
    stripped_seq = "AMGSPK"
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$site_residue, "S")
  expect_equal(result$site_peptide_pos, 4L)
})
