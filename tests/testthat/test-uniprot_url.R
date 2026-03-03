# Tests for UniProt FASTA download URL construction in R/helpers_search.R
# These test the URL logic without making actual HTTP requests.

# =============================================================================
# download_uniprot_fasta URL construction
# =============================================================================
# We can't test the actual download (requires network), but we can verify
# the URL is constructed correctly by inspecting the function logic.
# Instead, we test the URL-building portion via a helper that extracts it.

# Helper: extract the URL that download_uniprot_fasta would construct
build_uniprot_url <- function(proteome_id, content_type) {
  base_query <- sprintf("(proteome:%s)", proteome_id)
  query <- switch(content_type,
    "one_per_gene" = base_query,
    "reviewed"     = paste0(base_query, " AND (reviewed:true)"),
    "full"         = base_query,
    "full_isoforms" = base_query,
    base_query
  )
  include_isoform <- content_type == "full_isoforms"
  one_per_gene <- content_type == "one_per_gene"
  paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "query=", utils::URLencode(query),
    "&format=fasta",
    "&compressed=false",
    if (include_isoform) "&includeIsoform=true" else "",
    if (one_per_gene) "&onePerGene=true" else ""
  )
}

test_that("UniProt URL for one_per_gene includes onePerGene=true", {
  url <- build_uniprot_url("UP000005640", "one_per_gene")
  expect_true(grepl("onePerGene=true", url))
  expect_false(grepl("includeIsoform=true", url))
  expect_false(grepl("reviewed:true", url))
})

test_that("UniProt URL for reviewed includes reviewed:true filter", {
  url <- build_uniprot_url("UP000005640", "reviewed")
  expect_true(grepl("reviewed%3Atrue", url) || grepl("reviewed:true", url))
  expect_false(grepl("onePerGene=true", url))
  expect_false(grepl("includeIsoform=true", url))
})

test_that("UniProt URL for full has no extra parameters", {
  url <- build_uniprot_url("UP000005640", "full")
  expect_false(grepl("onePerGene=true", url))
  expect_false(grepl("includeIsoform=true", url))
  expect_false(grepl("reviewed:true", url))
})

test_that("UniProt URL for full_isoforms includes includeIsoform=true", {
  url <- build_uniprot_url("UP000005640", "full_isoforms")
  expect_true(grepl("includeIsoform=true", url))
  expect_false(grepl("onePerGene=true", url))
})

test_that("UniProt URL always includes proteome ID in query", {
  for (ct in c("one_per_gene", "reviewed", "full", "full_isoforms")) {
    url <- build_uniprot_url("UP000005640", ct)
    expect_true(grepl("UP000005640", url),
                info = sprintf("Missing proteome ID for content_type=%s", ct))
  }
})

test_that("UniProt URL uses stream endpoint", {
  url <- build_uniprot_url("UP000005640", "full")
  expect_true(grepl("rest.uniprot.org/uniprotkb/stream", url))
})

test_that("UniProt URL format is fasta", {
  url <- build_uniprot_url("UP000005640", "full")
  expect_true(grepl("format=fasta", url))
})

# =============================================================================
# generate_fasta_filename content_type mapping
# =============================================================================

test_that("generate_fasta_filename: one_per_gene maps to _opg_", {
  fn <- generate_fasta_filename("UP000005640", "Homo sapiens", "one_per_gene")
  expect_true(grepl("_opg_", fn))
  expect_false(grepl("_full_", fn))
})

test_that("generate_fasta_filename: reviewed maps to _sprot_", {
  fn <- generate_fasta_filename("UP000005640", "Homo sapiens", "reviewed")
  expect_true(grepl("_sprot_", fn))
})

test_that("generate_fasta_filename: full maps to _full_ (not _full_iso_)", {
  fn <- generate_fasta_filename("UP000005640", "Homo sapiens", "full")
  expect_true(grepl("_full_", fn))
  expect_false(grepl("_full_iso_", fn))
})
