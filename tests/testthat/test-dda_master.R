# Tests for the de novo Master Table join + species attribution.
# These lock in the v3.11.2x fixes:
#   - species must come from the taxonomy/LCA, never the mangled accession suffix
#   - the Sage / Casanovo / LCA join must match on the canonical (I/L) key
#   - negative-confidence (mass-mismatch) peptides must NOT be dropped by the
#     join — filtering is the caller's job (they still BLAST usefully)

test_that("build_dda_canonical_peptide strips mods + non-letters and uppercases", {
  expect_equal(build_dda_canonical_peptide("[Acetyl]-SAGGPPPGDPPK"), "SAGGPPPGDPPK")
  expect_equal(build_dda_canonical_peptide("PEPM+15.995TIDEK"), "PEPMTIDEK")
  expect_equal(build_dda_canonical_peptide("C[+57.021]VLK"), "CVLK")
  expect_equal(build_dda_canonical_peptide("peptidek"), "PEPTIDEK")
})

test_that("dda_blast_species: NCBI nr accession is NOT mangled when LCA is present", {
  lca <- data.frame(peptide = "ELVISLLVESK", lca_name = "Felidae",
                    stringsAsFactors = FALSE)
  sp <- dda_blast_species("ELVISLLVESK", "XP_025773238.1", lca)
  expect_equal(sp, "Felidae")          # regression: must NOT be "025773238.1"
})

test_that("dda_blast_species: UniProt mnemonic fallback when no LCA table", {
  expect_equal(
    dda_blast_species("PEPTIDEK", "sp|P12345|KRT5_FELCA", lca_tbl = NULL),
    "FELCA")
})

test_that("dda_blast_species: LCA join is I/L-insensitive and order-safe", {
  lca <- data.frame(peptide = c("AAACVLR", "ELVISLLVESK"),
                    lca_name = c("Bos taurus", "Felidae"),
                    stringsAsFactors = FALSE)
  sp <- dda_blast_species(c("ELVLSLLVESK", "AAACVLR"),   # note I->L in query
                          c("XP_1.1", "XP_2.1"), lca)
  expect_equal(sp, c("Felidae", "Bos taurus"))
})

test_that("build_denovo_master joins Casanovo + Sage + LCA on the canonical key", {
  classified <- data.frame(
    seq_stripped = c("PEPTIDEK", "ELVISLLVESK", "ELVISLLVESK", "MASSMISMATCHK"),
    score        = c(0.95, 0.60, 0.30, -0.50),
    match_type   = c("confirmed", "novel", "novel", "novel"),
    stringsAsFactors = FALSE)
  sage <- data.frame(peptide = "PEPTIDEK", proteins = "sp|P1|PROT_HUMAN",
                     stringsAsFactors = FALSE)
  lca <- data.frame(
    peptide    = c("ELVISLLVESK", "MASSMISMATCHK"),
    lca_name   = c("Felidae", "Panthera onca"),
    lca_rank   = c("family", "species"),
    category   = c("host", "host"),
    top_pident = c(95, 88),
    diagnostic = c(0L, 1L),
    stringsAsFactors = FALSE)

  m <- build_denovo_master(classified, sage, lca)

  expect_equal(nrow(m), 3)                       # 3 unique peptides (ELVIS deduped)

  elvis <- m[m$Peptide == "ELVISLLVESK", ]
  expect_equal(elvis$Confidence, 0.6)            # max score per peptide
  expect_equal(elvis$Species_or_clade, "Felidae")  # regression: LCA join populated
  expect_equal(elvis$Type, "host")

  pep <- m[m$Peptide == "PEPTIDEK", ]
  expect_true(pep$Found_by_Sage)
  expect_equal(pep$Sage_protein, "sp|P1|PROT_HUMAN")
  expect_true(is.na(pep$Species_or_clade))       # no nr hit -> NA, not error

  mm <- m[m$Peptide == "MASSMISMATCHK", ]
  expect_true(mm$Diagnostic)                     # diagnostic flag -> logical TRUE
  expect_true(mm$Confidence < 0)                 # negative-score peptide KEPT
  expect_equal(mm$Species_or_clade, "Panthera onca")
})

test_that("build_denovo_master tolerates missing Sage + LCA (de-novo-only)", {
  classified <- data.frame(seq_stripped = "NOVELONLYK", score = 0.8,
                           match_type = "novel", stringsAsFactors = FALSE)
  m <- build_denovo_master(classified, sage_psms = NULL, lca = NULL)
  expect_equal(nrow(m), 1)
  expect_false(m$Found_by_Sage)
  expect_true(is.na(m$Sage_protein))
  expect_false("Species_or_clade" %in% names(m))  # no LCA -> no species columns
})

test_that("build_denovo_master returns empty frame on empty input", {
  expect_equal(nrow(build_denovo_master(NULL)), 0)
  expect_equal(nrow(build_denovo_master(data.frame())), 0)
})

test_that("build_denovo_score_calibration: target hit-rate + decoy FDR", {
  classified <- data.frame(
    seq_stripped = c("HIGHAAAAK", "HIGHBBBBR", "MIDCCCCK", "MIDDDDDR",
                     "LOWEEEEK", "LOWFFFFR"),
    score = c(0.9, 0.8, 0.4, 0.3, -0.5, -0.6),
    match_type = "novel", stringsAsFactors = FALSE)
  blast <- data.frame(peptide = c("HIGHAAAAK", "HIGHBBBBR", "MIDCCCCK", "MIDDDDDR"),
                      pident = c(98, 95, 92, 90), stringsAsFactors = FALSE)
  decoy <- data.frame(peptide = "LOWEEEEK", pident = 82, stringsAsFactors = FALSE)

  cal <- build_denovo_score_calibration(classified, blast, decoy, bin = 0.5)
  expect_true(all(c("hit_rate", "mean_pident", "decoy_hit_rate", "fdr", "cum_fdr")
                  %in% names(cal)))
  # peptides at/above 0 are all real (no decoy hits) -> cum FDR 0
  hi <- cal[cal$bin >= 0, ]
  expect_true(all(hi$cum_fdr == 0 | is.na(hi$cum_fdr)))
  # including the negative bins pulls in a decoy hit -> cum FDR rises
  expect_true(max(cal$cum_fdr, na.rm = TRUE) > 0)
})

test_that("build_denovo_score_calibration: target-only (no decoy) omits FDR", {
  classified <- data.frame(seq_stripped = c("AAAAAAAK", "BBBBBBBR"),
                           score = c(0.7, 0.2), match_type = "novel",
                           stringsAsFactors = FALSE)
  blast <- data.frame(peptide = "AAAAAAAK", pident = 95, stringsAsFactors = FALSE)
  cal <- build_denovo_score_calibration(classified, blast, NULL, bin = 0.5)
  expect_true("hit_rate" %in% names(cal))
  expect_false("fdr" %in% names(cal))
  expect_equal(nrow(build_denovo_score_calibration(NULL, NULL)), 0)
})
