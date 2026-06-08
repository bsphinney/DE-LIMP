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
  expect_equal(elvis$Casanovo_score, 0.6)        # max score per peptide
  expect_equal(elvis$Species_or_clade, "Felidae")  # regression: LCA join populated
  expect_equal(elvis$Type, "host")

  pep <- m[m$Peptide == "PEPTIDEK", ]
  expect_true(pep$Found_by_Sage)
  expect_equal(pep$Sage_protein, "sp|P1|PROT_HUMAN")
  expect_true(is.na(pep$Species_or_clade))       # no nr hit -> NA, not error

  mm <- m[m$Peptide == "MASSMISMATCHK", ]
  expect_true(mm$Diagnostic)                     # diagnostic flag -> logical TRUE
  expect_true(mm$Casanovo_score < 0)             # negative-score peptide KEPT
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

test_that("build_denovo_master surfaces query coverage + e-value from blast", {
  classified <- data.frame(seq_stripped = "FNGGYDNFNFNGLTGTGVLTESNK",  # 24 aa
                           score = 0.036, match_type = "novel", stringsAsFactors = FALSE)
  blast <- data.frame(peptide = "FNGGYDNFNFNGLTGTGVLTESNK", subject = "UJW42468.1",
                      pident = 100, length = 18, evalue = 0.04, bitscore = 39.7,
                      stringsAsFactors = FALSE)
  m <- build_denovo_master(classified, NULL, NULL, blast)
  expect_equal(m$Query_coverage, 75)        # 18/24 — partial, not a full match
  expect_equal(m$E_value, 0.04)
  expect_equal(m$Bitscore, 39.7)
})

test_that("build_denovo_master returns empty frame on empty input", {
  expect_equal(nrow(build_denovo_master(NULL)), 0)
  expect_equal(nrow(build_denovo_master(data.frame())), 0)
})

test_that("best_blast_hit_per_peptide carries the best hit's alignment when present", {
  blast <- data.frame(
    peptide  = c("PEPTIDEK", "PEPTIDEK"),
    length   = c(8, 8),
    evalue   = c(0.1, 1e-5),
    bitscore = c(20, 40),                 # second row is the better hit
    qseq     = c("PEPTIDEK", "PEPTIDEK"),
    sseq     = c("PEPTLDEK", "PEPTVDEK"),
    btop     = c("4LL3", "4VV3"),
    qstart   = c(1L, 1L),
    stringsAsFactors = FALSE)
  bb <- best_blast_hit_per_peptide(blast)
  expect_true(all(c("blast_qseq", "blast_sseq", "blast_btop", "blast_qstart") %in% names(bb)))
  expect_equal(bb$blast_btop, "4VV3")     # the higher-bitscore hit's alignment
})

test_that("build_denovo_master adds CWI/HCA/Call from btop + aa_lookup", {
  classified <- data.frame(seq_stripped = "AAAAAAQAAAAAAAAAAAAA", score = 0.7,
                           match_type = "novel", stringsAsFactors = FALSE)
  blast <- data.frame(peptide = "AAAAAAQAAAAAAAAAAAAA", subject = "XP_1.1",
                      pident = 95, length = 20, evalue = 1e-4, bitscore = 50,
                      qseq = "AAAAAAQAAAAAAAAAAAAA", sseq = "AAAAAAHAAAAAAAAAAAAA",
                      btop = "6QH13", qstart = 1L, stringsAsFactors = FALSE)
  k <- gsub("I", "L", build_dda_canonical_peptide("AAAAAAQAAAAAAAAAAAAA"))
  aa <- stats::setNames(paste(c(rep(0.99, 6), 0.30, rep(0.99, 13)), collapse = ","), k)
  m <- build_denovo_master(classified, NULL, NULL, blast, aa_lookup = aa)
  expect_true(all(c("CWI", "HCA", "Call") %in% names(m)))
  expect_true(m$CWI >= 95)        # only mismatch is low-confidence -> high CWI
  expect_equal(m$HCA, 100)        # every high-confidence residue agrees
  expect_equal(m$Call, "Confident")

  # No aa_lookup -> CWI still computes (unweighted), HCA NA, Call NA
  m2 <- build_denovo_master(classified, NULL, NULL, blast, aa_lookup = NULL)
  expect_true("CWI" %in% names(m2) && is.finite(m2$CWI))
  expect_true(is.na(m2$HCA))

  # blast without btop -> no CWI columns at all (never fabricated)
  m3 <- build_denovo_master(classified, NULL, NULL,
                            blast[, c("peptide","subject","pident","length","evalue","bitscore")])
  expect_false("CWI" %in% names(m3))
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

test_that("parse_btop reconstructs the gapped alignment from btop + qseq", {
  expect_equal(parse_btop("20", "ACDEFGHIKLMNPQRSTVWY")$qaln, "ACDEFGHIKLMNPQRSTVWY")
  m <- parse_btop("6QH13", "AAAAAAQAAAAAAAAAAAAA")        # 6 match, Q->H, 13 match
  expect_equal(nchar(m$qaln), 20); expect_equal(nchar(m$saln), 20)
  expect_equal(substr(m$qaln, 7, 7), "Q"); expect_equal(substr(m$saln, 7, 7), "H")
  g <- parse_btop("3-S4", "ACDEFGH")                      # gap in query
  expect_true(grepl("-", g$qaln)); expect_equal(nchar(g$qaln), nchar(g$saln))
})

test_that("compute_cwi: low-confidence mismatches barely penalize (the 'good hit' case)", {
  aln <- parse_btop("6QH13", "AAAAAAQAAAAAAAAAAAAA")
  aa <- rep(0.99, 20); aa[7] <- 0.30                      # the only mismatch is low-confidence
  cw <- compute_cwi(aln$qaln, aln$saln, aa, qstart = 1)
  expect_true(cw$cwi >= 95)                               # high CWI despite a mismatch
  expect_equal(cw$hca, 100)                               # every high-confidence residue agrees
  # opposite: a HIGH-confidence mismatch should drag CWI + HCA down
  aa2 <- rep(0.99, 20)                                    # mismatch position now high-confidence
  cw2 <- compute_cwi(aln$qaln, aln$saln, aa2, qstart = 1)
  expect_true(cw2$cwi < cw$cwi); expect_true(cw2$hca < 100)
})

test_that("build_denovo_score_calibration: min_length drops short peptides from BOTH sides", {
  classified <- data.frame(
    seq_stripped = c("SHORTK", "LONGPEPTIDEAAAK", "LONGPEPTIDEBBBR"),
    score = c(-0.5, 0.7, 0.6), match_type = "novel", stringsAsFactors = FALSE)
  blast <- data.frame(peptide = "LONGPEPTIDEAAAK", pident = 95, stringsAsFactors = FALSE)
  decoy <- data.frame(peptide = "SHORTK", pident = 80, stringsAsFactors = FALSE)  # spurious short decoy
  c7 <- build_denovo_score_calibration(classified, blast, decoy, bin = 0.5, min_length = 7)
  # the short peptide's score bin (negative) is gone -> its decoy hit can't count
  expect_equal(nrow(c7[c7$bin < 0, ]), 0)
  expect_true(all(nchar(gsub("[^0-9.-]", "", as.character(c7$bin))) >= 0))  # sane bins
})

test_that("build_denovo_protein_groups: parsimony + razor counts (FragPipe/IDPicker model)", {
  # protA explains {p1,p2,shared}; protB explains {shared,p4}. Greedy picks protA
  # first (covers 3), razor-assigns shared to it; protB keeps only p4.
  blast <- data.frame(
    peptide  = c("PEPONEK", "PEPTWOK", "SHAREDDK", "SHAREDDK", "PEPFOURK"),
    subject  = c("protA", "protA", "protA", "protB", "protB"),
    pident   = c(95, 95, 95, 95, 95),
    bitscore = c(50, 50, 50, 40, 50),
    sstart   = c(1, 10, 20, 5, 15), send = c(8, 17, 27, 12, 23),
    stringsAsFactors = FALSE)
  g <- build_denovo_protein_groups(blast, master = NULL, min_pident = 50)
  expect_equal(nrow(g), 2)
  A <- g[g$Subject == "protA", ]; B <- g[g$Subject == "protB", ]
  expect_equal(A$Razor_peptides, 3); expect_equal(A$Unique_peptides, 2); expect_equal(A$Total_peptides, 3)
  expect_equal(B$Razor_peptides, 1); expect_equal(B$Unique_peptides, 1); expect_equal(B$Total_peptides, 2)
  expect_equal(g$Subject[1], "protA")                      # ordered by razor desc
})

test_that("build_denovo_protein_groups: min_pident drops weak edges; empty in -> empty out", {
  blast <- data.frame(peptide = c("AAAK", "AAAK"), subject = c("hi", "lo"),
                      pident = c(95, 40), bitscore = c(50, 20),
                      sstart = c(1, 1), send = c(4, 4), stringsAsFactors = FALSE)
  g <- build_denovo_protein_groups(blast, min_pident = 50)
  expect_equal(nrow(g), 1); expect_equal(g$Subject, "hi")   # the 40% edge is excluded
  expect_equal(nrow(build_denovo_protein_groups(NULL)), 0)
})

test_that("build_protein_coverage_track: whole protein from residue 1 + confidence colouring", {
  blast <- data.frame(peptide = "ACDEFG", subject = "X", pident = 83, bitscore = 20,
                      qseq = "ACDEFG", sseq = "ACDQFG", btop = "3EQ2",
                      qstart = 1, sstart = 10, send = 15, stringsAsFactors = FALSE)
  tr <- build_protein_coverage_track("X", blast, aa_lookup = NULL)
  expect_equal(tr$pos[1], 1L)                                # whole protein: starts at residue 1
  expect_equal(max(tr$pos), 15L)
  expect_equal(tr$state[tr$pos < 10], rep("uncovered", 9))   # N-terminus before the peptide = dots
  at <- function(p, col) tr[[col]][tr$pos == p]
  expect_equal(at(10, "ref"), "A"); expect_equal(at(13, "ref"), "Q")  # ref = subject residues
  expect_equal(at(13, "obs"), "E")                           # de novo substitution E vs ref Q
  expect_equal(at(10, "state"), "match")
  expect_equal(at(13, "state"), "amber")                     # difference, no score
  k <- gsub("I", "L", build_dda_canonical_peptide("ACDEFG"))
  hi <- stats::setNames(paste(rep(0.99, 6), collapse = ","), k)
  th <- build_protein_coverage_track("X", blast, aa_lookup = hi)
  expect_equal(th$state[th$pos == 13], "variant")            # confident difference
  lo <- stats::setNames(paste(c(0.99,0.99,0.99,0.30,0.99,0.99), collapse = ","), k)
  tl <- build_protein_coverage_track("X", blast, aa_lookup = lo)
  expect_equal(tl$state[tl$pos == 13], "error")              # low-confidence difference
})

test_that("build_trusted_substitutions: confident mismatches only, bracketed flag", {
  aln <- parse_btop("6QH13", "AAAAAAQAAAAAAAAAAAAA")    # Q->H mismatch at position 7
  aa <- rep(0.99, 20)
  ts <- build_trusted_substitutions(aln$qaln, aln$saln, aa, qstart = 1)
  expect_equal(nrow(ts), 1)
  expect_equal(ts$pos, 7); expect_equal(ts$ref, "H"); expect_equal(ts$denovo, "Q")
  expect_true(ts$bracketed)                              # flanked by confident A matches
  # a low-confidence mismatch is NOT a trusted substitution (likely de novo error)
  aa2 <- rep(0.99, 20); aa2[7] <- 0.30
  expect_equal(nrow(build_trusted_substitutions(aln$qaln, aln$saln, aa2, qstart = 1)), 0)
})

test_that("compute_cwi: DIAMOND-masked 'X' residues are excluded (not scored as differences)", {
  # query has a masked low-complexity stretch (X); without exclusion these would
  # be counted as high-confidence mismatches and tank CWI/HCA.
  qa <- "ACDXXXFGHK"; sa <- "ACDEEEFGHK"   # positions 4-6 masked in the query
  aa <- rep(0.99, 10)
  cw <- compute_cwi(qa, sa, aa, qstart = 1)
  expect_equal(cw$n_pos, 7)                # 10 aligned - 3 masked
  expect_equal(cw$cwi, 100)                # every NON-masked residue matches
  expect_equal(cw$hca, 100)
})

test_that("build_protein_coverage_track: masked 'X' positions get state 'masked', not a difference", {
  blast <- data.frame(peptide = "ACXFG", subject = "X", pident = 80, bitscore = 20,
                      qseq = "ACXFG", sseq = "ACDFG", btop = "2XD2",
                      qstart = 1, sstart = 10, send = 14, stringsAsFactors = FALSE)
  tr <- build_protein_coverage_track("X", blast, aa_lookup = NULL)
  expect_equal(tr$state[tr$pos == 12], "masked")   # the X position is masked, not a substitution
  expect_equal(tr$n_mismatch[tr$pos == 12], 0)     # masked excluded from the mismatch tally
})

test_that("build_protein_coverage_track: a position is a difference only when the CONSENSUS differs", {
  # 3 peptides match the reference at every position; 1 peptide misreads position 4.
  # The consensus still matches, so position 4 must be 'match' (regression: it was
  # wrongly coloured as a difference whenever ANY peptide mismatched).
  blast <- data.frame(
    peptide  = c("ACDEFG", "ACDEFG", "ACDEFG", "ACDQFG"),
    subject  = "X", pident = c(100,100,100,83), bitscore = c(40,40,40,20),
    qseq = c("ACDEFG","ACDEFG","ACDEFG","ACDQFG"),
    sseq = c("ACDEFG","ACDEFG","ACDEFG","ACDEFG"),
    btop = c("6","6","6","3QE2"), qstart = 1, sstart = 10, send = 15,
    stringsAsFactors = FALSE)
  tr <- build_protein_coverage_track("X", blast, aa_lookup = NULL)
  expect_equal(tr$state[tr$pos == 13], "match")   # consensus E == ref E -> match, not a difference
  expect_equal(tr$obs[tr$pos == 13], "E")
  expect_true(tr$n_mismatch[tr$pos == 13] >= 1)   # the lone misread is still counted
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
