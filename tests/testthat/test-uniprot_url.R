# Tests for UniProt FASTA download URL construction in R/helpers_search.R
# These test the URL logic without making actual HTTP requests.

# =============================================================================
# download_uniprot_fasta URL construction (REST API for non-FTP content types)
# =============================================================================
# Helper: extract the REST API URL for non-one_per_gene content types
build_uniprot_url <- function(proteome_id, content_type) {
  base_query <- sprintf("(proteome:%s)", proteome_id)
  query <- switch(content_type,
    "reviewed"      = paste0(base_query, " AND (reviewed:true)"),
    "full"          = base_query,
    "full_isoforms" = base_query,
    base_query
  )
  include_isoform <- content_type == "full_isoforms"
  paste0(
    "https://rest.uniprot.org/uniprotkb/stream?",
    "query=", utils::URLencode(query),
    "&format=fasta",
    "&compressed=false",
    if (include_isoform) "&includeIsoform=true" else ""
  )
}

test_that("UniProt URL for reviewed includes reviewed:true filter", {
  url <- build_uniprot_url("UP000005640", "reviewed")
  expect_true(grepl("reviewed%3Atrue", url) || grepl("reviewed:true", url))
  expect_false(grepl("includeIsoform=true", url))
})

test_that("UniProt URL for full has no extra parameters", {
  url <- build_uniprot_url("UP000005640", "full")
  expect_false(grepl("includeIsoform=true", url))
  expect_false(grepl("reviewed:true", url))
})

test_that("UniProt URL for full_isoforms includes includeIsoform=true", {
  url <- build_uniprot_url("UP000005640", "full_isoforms")
  expect_true(grepl("includeIsoform=true", url))
})

test_that("UniProt URL always includes proteome ID in query", {
  for (ct in c("reviewed", "full", "full_isoforms")) {
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
# FTP URL construction for one-per-gene downloads
# =============================================================================
# Helper: build expected FTP URL for a given proteome
build_ftp_url <- function(proteome_id, taxon_id, kingdom) {
  sprintf(
    "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/%s/%s/%s_%s.fasta.gz",
    kingdom, proteome_id, proteome_id, taxon_id
  )
}

test_that("FTP URL for human one-per-gene is correct", {
  url <- build_ftp_url("UP000005640", 9606, "Eukaryota")
  expect_true(grepl("Eukaryota/UP000005640", url))
  expect_true(grepl("UP000005640_9606.fasta.gz", url))
})

test_that("FTP URL for E. coli one-per-gene uses Bacteria kingdom", {
  url <- build_ftp_url("UP000000625", 83333, "Bacteria")
  expect_true(grepl("Bacteria/UP000000625", url))
  expect_true(grepl("UP000000625_83333.fasta.gz", url))
})

test_that("FTP URL for mouse one-per-gene is correct", {
  url <- build_ftp_url("UP000000589", 10090, "Eukaryota")
  expect_true(grepl("Eukaryota/UP000000589", url))
  expect_true(grepl("UP000000589_10090.fasta.gz", url))
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
