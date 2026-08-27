# Contaminant detection: ONE definition, several conventions.
#
# The failure this guards is silent and total. The app used to test only
# `grepl("^Cont_", x)` — the tag its own FASTA builder stamps — so a user who
# searched against a contaminant database they built themselves had every
# "Exclude contaminants" checkbox quietly do nothing and the Contaminant
# Analysis tab report none present (issue #67, DIA-NN 2.6.1 + cRAP_).
#
# The opposite error is worse: a false positive DELETES a real protein from the
# analysis. Both directions are asserted.

test_that("every shipped tagging convention is recognised", {
  expect_true(all(is_contaminant_accession(c(
    "Cont_P00761",              # DE-LIMP / pipeline skill
    "cRAP_P00761", "cRAP-P00761",  # GPM cRAP as several builders prefix it
    "CON__P00761",              # MaxQuant
    "contam_sp|P00761|TRYP_PIG" # FragPipe / Philosopher
  ))))
})

test_that("the tag is found inside a UniProt-style header, not just at the start", {
  # DIA-NN may report the whole header, so an anchor of "^" alone would miss it.
  expect_true(is_contaminant_accession("sp|Cont_P00761|TRYP_PIG"))
  expect_true(is_contaminant_accession("tr|cRAP_Q9XYZ1|SOMETHING"))
})

test_that("real proteins whose names merely contain the letters are NOT flagged", {
  # A false positive removes real data. CONTACTIN and Q9CONT both contain
  # "CONT"; neither is a contaminant.
  expect_false(any(is_contaminant_accession(c(
    "P12345", "sp|P12345|ALBU_HUMAN", "CONTACTIN_HUMAN", "Q9CONT", "CONTRA_MOUSE"
  ))))
})

test_that("stripping the tag leaves an accession that still resolves", {
  # This feeds the UniProt links, so the result has to be the bare accession.
  expect_equal(strip_contaminant_prefix("Cont_P04264"), "P04264")
  expect_equal(strip_contaminant_prefix("cRAP_P00761"), "P00761")
  expect_equal(strip_contaminant_prefix("sp|Cont_P00761|TRYP_PIG"), "sp|P00761|TRYP_PIG")
  expect_equal(strip_contaminant_prefix("P12345"), "P12345")   # untouched
})

test_that("a user-supplied tag is honoured and taken literally", {
  old <- getOption("delimp.contaminant_tags")
  on.exit(options(delimp.contaminant_tags = old), add = TRUE)

  expect_false(is_contaminant_accession("myCONTAM_P1"))
  options(delimp.contaminant_tags = "myCONTAM_")
  expect_true(is_contaminant_accession("myCONTAM_P1"))
  expect_false(is_contaminant_accession("P12345"))       # still no false positives

  # Typed into a text box, so it may contain regex metacharacters. It must match
  # literally rather than erroring or matching everything.
  options(delimp.contaminant_tags = "MY(CON[")
  expect_true(is_contaminant_accession("MY(CON[P1"))
  expect_false(is_contaminant_accession("MYXCONYP1"))
})

test_that("a protein group counts as contaminant only if EVERY accession is one", {
  # A group with a real protein in it is real — dropping it would lose that protein.
  expect_true(is_contaminant_protein_group("Cont_P00761;cRAP_P02769"))
  expect_false(is_contaminant_protein_group("Cont_P00761;P12345"))
  expect_false(is_contaminant_protein_group("P12345"))
})

test_that("detected_contaminant_tags reports what is actually present", {
  expect_equal(detected_contaminant_tags(c("Cont_P1", "P2")), "Cont_")
  expect_setequal(detected_contaminant_tags(c("cRAP_P1", "CON__P2")),
                  c("cRAP", "MaxQuant"))
  expect_length(detected_contaminant_tags(c("P1", "P2")), 0)
})

test_that("empty and NA input do not error", {
  expect_length(is_contaminant_accession(character(0)), 0)
  expect_false(is_contaminant_accession(NA))
})
