# Claude Export Consistency Test
#
# Tests that Claude's analysis of the DE-LIMP export zip file remains
# structurally consistent over time. Uses an actual zip file exported
# from the app (tests/fixtures/example_claude_export.zip).
#
# The export contains:
#   PROMPT.md, DE_Results_Full.csv, Expression_Matrix.csv, QC_Metrics.csv,
#   Group_Assignments.csv, Analysis_Parameters.txt, Methods_and_References.txt,
#   Reproducibility_Code.R, Session.rds
#
# What we test:
#   1. Zip file structure hasn't changed (all expected files present)
#   2. Data files are internally consistent (row counts, contrasts, etc.)
#   3. Claude's API response mentions key cross-comparison proteins
#   4. Claude's response has all requested sections
#   5. Over time: gene mention overlap doesn't drift below threshold
#
# Requires: ANTHROPIC_API_KEY env var for the API test
#
# Run: Rscript -e 'testthat::test_file("tests/testthat/test-claude_export.R")'

# Known ground truth from the example dataset export
EXPECTED_FILES <- c(
  "DE_Results_Full.csv", "QC_Metrics.csv", "Expression_Matrix.csv",
  "Group_Assignments.csv", "Analysis_Parameters.txt", "PROMPT.md"
)
EXPECTED_CONTRASTS <- c(
  "affinisepACN - Affinisep", "affinisepIPA - Affinisep",
  "affinisepIPA - affinisepACN",
  "Evosep - Affinisep", "Evosep - affinisepACN", "Evosep - affinisepIPA"
)
# Cross-comparison markers (significant in >= 3 comparisons)
KEY_PROTEINS <- c("NAMPT", "ALDH3A1", "LGALS1", "C3", "ITGA2",
                   "PRSS23", "TIMP1", "GSTM3", "CGA")
# Stable biomarkers (lowest CV)
STABLE_MARKERS <- c("MED14", "TUBA1B", "TBC1D17")

EXPECTED_N_PROTEINS <- 4854
EXPECTED_N_SAMPLES <- 12

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
get_fixture_zip <- function() {
  zip_path <- file.path(project_root, "tests", "fixtures", "example_claude_export.zip")
  if (!file.exists(zip_path)) return(NULL)
  zip_path
}

unzip_to_temp <- function(zip_path) {
  tmp <- tempfile("claude_export_")
  dir.create(tmp)
  unzip(zip_path, exdir = tmp)
  tmp
}

call_claude_api <- function(prompt_text, de_csv_head, api_key) {
  user_content <- paste0(
    prompt_text, "\n\n",
    "--- DE_Results_Full.csv (first 300 rows) ---\n",
    de_csv_head
  )
  body <- list(
    model = "claude-sonnet-4-20250514",
    max_tokens = 4096,
    messages = list(list(role = "user", content = user_content))
  )
  resp <- httr2::request("https://api.anthropic.com/v1/messages") |>
    httr2::req_headers(
      `x-api-key` = api_key,
      `anthropic-version` = "2023-06-01",
      `content-type` = "application/json"
    ) |>
    httr2::req_body_json(body) |>
    httr2::req_timeout(120) |>
    httr2::req_perform()
  result <- httr2::resp_body_json(resp)
  # Return both text and model metadata
  list(
    text = result$content[[1]]$text,
    model = result$model %||% body$model,
    input_tokens = result$usage$input_tokens %||% NA,
    output_tokens = result$usage$output_tokens %||% NA
  )
}

extract_findings <- function(response_text) {
  text_lower <- tolower(response_text)

  # Sections present
  sections <- c("overview", "qc assessment", "key findings",
                 "cross-comparison", "biomarker", "biological interpretation",
                 "how this analysis works")
  sections_found <- vapply(sections, function(s) grepl(s, text_lower), logical(1))

  # Gene names mentioned (uppercase 2-10 char words)
  gene_pattern <- "\\b[A-Z][A-Z0-9]{1,9}\\b"
  genes_raw <- unique(regmatches(response_text,
    gregexpr(gene_pattern, response_text))[[1]])
  noise <- c("THE", "AND", "FOR", "ARE", "NOT", "BUT", "HAS", "WAS", "WITH",
             "THIS", "THAT", "FROM", "WERE", "BEEN", "HAVE", "EACH", "ALL",
             "CAN", "MORE", "BOTH", "ALSO", "THAN", "MAY", "DIA", "CSV",
             "QC", "DE", "MS1", "MS2", "FDR", "PCA", "NES", "RNA", "DNA",
             "TOP", "LOG", "PER", "USE", "RUN", "TWO", "ONE", "NEW", "KEY",
             "LOW", "HIGH", "UP", "DOWN", "LC", "DDA", "ANY", "ITS",
             "RAW", "FASTA", "SLURM", "HPC", "SSH")
  genes_mentioned <- sort(setdiff(genes_raw, noise))

  # Key protein coverage
  key_found <- vapply(KEY_PROTEINS, function(g)
    grepl(g, response_text, fixed = TRUE), logical(1))
  stable_found <- vapply(STABLE_MARKERS, function(g)
    grepl(g, response_text, fixed = TRUE), logical(1))

  # Directional language
  n_up <- max(0L, length(gregexpr("upregulated|up-regulated|increased|elevated",
    text_lower)[[1]]))
  n_down <- max(0L, length(gregexpr("downregulated|down-regulated|decreased|reduced",
    text_lower)[[1]]))

  # Statistical rigor: does Claude cite specific numbers?
  n_pvalues <- length(gregexpr("p[- ]?value|adj\\.?p|fdr|q[- ]?value", text_lower)[[1]])
  n_foldchanges <- length(gregexpr("fold[- ]?change|log.?fc|logfc", text_lower)[[1]])
  n_numbers <- length(gregexpr("\\b\\d+\\.\\d+\\b", response_text)[[1]])  # decimal numbers

  # Hedging language (scientific caution)
  n_hedges <- length(gregexpr("suggest|may|could|potential|appears|likely|possible",
    text_lower)[[1]])

  # Confidence language (definitive statements)
  n_confident <- length(gregexpr("clearly|strongly|significantly|definitively|robust",
    text_lower)[[1]])

  list(
    sections_found = sections_found,
    genes_mentioned = genes_mentioned,
    n_genes = length(genes_mentioned),
    key_proteins_found = key_found,
    stable_markers_found = stable_found,
    n_key_found = sum(key_found),
    n_stable_found = sum(stable_found),
    n_up_mentions = n_up,
    n_down_mentions = n_down,
    n_pvalue_refs = n_pvalues,
    n_fc_refs = n_foldchanges,
    n_decimal_numbers = n_numbers,
    n_hedges = n_hedges,
    n_confident = n_confident,
    word_count = length(strsplit(response_text, "\\s+")[[1]]),
    timestamp = Sys.time()
  )
}

# ===========================================================================
# TEST 1: Zip file structure
# ===========================================================================
test_that("export zip contains all expected files", {
  zip_path <- get_fixture_zip()
  skip_if(is.null(zip_path), "No fixture zip at tests/fixtures/example_claude_export.zip")

  contents <- unzip(zip_path, list = TRUE)$Name
  for (f in EXPECTED_FILES) {
    expect_true(f %in% contents, info = paste("Missing file:", f))
  }
  expect_true("PROMPT.md" %in% contents)
})

# ===========================================================================
# TEST 2: DE results CSV internal consistency
# ===========================================================================
test_that("DE_Results_Full.csv has expected contrasts and structure", {
  zip_path <- get_fixture_zip()
  skip_if(is.null(zip_path), "No fixture zip")

  tmp <- unzip_to_temp(zip_path)
  on.exit(unlink(tmp, recursive = TRUE))

  de <- read.csv(file.path(tmp, "DE_Results_Full.csv"))

  # Expected columns
  expect_true(all(c("Protein.Group", "logFC", "adj.P.Val", "Contrast") %in% names(de)))

  # Expected contrasts
  contrasts_found <- unique(de$Contrast)
  for (ct in EXPECTED_CONTRASTS) {
    expect_true(ct %in% contrasts_found,
      info = paste("Missing contrast:", ct))
  }

  # Expected protein count per contrast
  per_contrast <- table(de$Contrast)
  expect_true(all(per_contrast == EXPECTED_N_PROTEINS),
    info = paste("Not all contrasts have", EXPECTED_N_PROTEINS, "proteins"))
})

# ===========================================================================
# TEST 3: Expression matrix dimensions
# ===========================================================================
test_that("Expression_Matrix.csv has expected dimensions", {
  zip_path <- get_fixture_zip()
  skip_if(is.null(zip_path), "No fixture zip")

  tmp <- unzip_to_temp(zip_path)
  on.exit(unlink(tmp, recursive = TRUE))

  expr <- read.csv(file.path(tmp, "Expression_Matrix.csv"))

  # Should have Protein.Group + 12 sample columns
  expect_equal(ncol(expr), EXPECTED_N_SAMPLES + 1)
  expect_true("Protein.Group" %in% names(expr))
  expect_equal(nrow(expr), EXPECTED_N_PROTEINS)
})

# ===========================================================================
# TEST 4: QC metrics match sample count
# ===========================================================================
test_that("QC_Metrics.csv has one row per sample", {
  zip_path <- get_fixture_zip()
  skip_if(is.null(zip_path), "No fixture zip")

  tmp <- unzip_to_temp(zip_path)
  on.exit(unlink(tmp, recursive = TRUE))

  qc <- read.csv(file.path(tmp, "QC_Metrics.csv"))
  expect_equal(nrow(qc), EXPECTED_N_SAMPLES)
  expect_true("Run" %in% names(qc))
  expect_true("Proteins" %in% names(qc))
})

# ===========================================================================
# TEST 5: PROMPT.md contains data summaries
# ===========================================================================
test_that("PROMPT.md references expected contrasts and proteins", {
  zip_path <- get_fixture_zip()
  skip_if(is.null(zip_path), "No fixture zip")

  tmp <- unzip_to_temp(zip_path)
  on.exit(unlink(tmp, recursive = TRUE))

  prompt <- paste(readLines(file.path(tmp, "PROMPT.md")), collapse = "\n")

  # Should reference the cross-comparison proteins
  for (gene in KEY_PROTEINS) {
    expect_true(grepl(gene, prompt, fixed = TRUE),
      info = paste("PROMPT.md missing key protein:", gene))
  }

  # Should mention the number of comparisons
  expect_true(grepl("6 comparison", prompt))

  # Should have section headers
  expect_true(grepl("## Overview", prompt))
  expect_true(grepl("## Key Findings", prompt))
  expect_true(grepl("## Cross-Comparison", prompt))
})

# ===========================================================================
# TEST 6: Claude API analysis — structural consistency
# ===========================================================================
test_that("Claude analysis mentions key proteins and has all sections", {
  skip_if_not_installed("httr2")

  api_key <- Sys.getenv("ANTHROPIC_API_KEY")
  skip_if(!nzchar(api_key), "ANTHROPIC_API_KEY not set")

  zip_path <- get_fixture_zip()
  skip_if(is.null(zip_path), "No fixture zip")

  tmp <- unzip_to_temp(zip_path)
  on.exit(unlink(tmp, recursive = TRUE))

  # Read prompt and DE data
  prompt_text <- paste(readLines(file.path(tmp, "PROMPT.md")), collapse = "\n")
  de <- read.csv(file.path(tmp, "DE_Results_Full.csv"))
  de_head <- paste(capture.output(write.csv(head(de, 300), row.names = FALSE)),
    collapse = "\n")

  # Call Claude API
  api_result <- call_claude_api(prompt_text, de_head, api_key)
  response <- api_result$text
  findings <- extract_findings(response)

  # Store model metadata in findings for drift tracking
  findings$model <- api_result$model
  findings$input_tokens <- api_result$input_tokens
  findings$output_tokens <- api_result$output_tokens

  # --- Save for golden comparison ---
  golden_dir <- file.path(project_root, "tests", "golden")
  if (!dir.exists(golden_dir)) dir.create(golden_dir, recursive = TRUE)

  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  writeLines(response, file.path(golden_dir, paste0("claude_response_", ts, ".md")))
  saveRDS(findings, file.path(golden_dir, paste0("claude_findings_", ts, ".rds")))

  # --- Structural assertions ---

  # All major sections present
  expect_true(findings$sections_found["overview"],
    info = "Missing 'Overview' section")
  expect_true(findings$sections_found["key findings"],
    info = "Missing 'Key Findings' section")
  expect_true(findings$sections_found["cross-comparison"],
    info = "Missing 'Cross-Comparison' section")
  expect_true(findings$sections_found["biological interpretation"],
    info = "Missing 'Biological Interpretation' section")

  # Must mention at least 5 of 9 key cross-comparison proteins
  expect_true(findings$n_key_found >= 5,
    info = sprintf("Only %d/%d key proteins mentioned: %s",
      findings$n_key_found, length(KEY_PROTEINS),
      paste(KEY_PROTEINS[findings$key_proteins_found], collapse = ", ")))

  # Must mention at least 1 stable biomarker
  expect_true(findings$n_stable_found >= 1,
    info = sprintf("No stable biomarkers mentioned (expected: %s)",
      paste(STABLE_MARKERS, collapse = ", ")))

  # Must discuss both directions

  expect_true(findings$n_up_mentions > 0, info = "No mention of upregulation")
  expect_true(findings$n_down_mentions > 0, info = "No mention of downregulation")

  # Substantive response (>= 500 words)
  expect_true(findings$word_count >= 500,
    info = paste("Response only", findings$word_count, "words"))

  # Must mention at least 10 gene names total
  expect_true(findings$n_genes >= 10,
    info = paste("Only", findings$n_genes, "genes mentioned"))

  # --- Compare against previous golden run ---
  golden_files <- sort(list.files(golden_dir, "^claude_findings_.*\\.rds$",
    full.names = TRUE), decreasing = TRUE)
  current_file <- file.path(golden_dir, paste0("claude_findings_", ts, ".rds"))
  golden_files <- setdiff(golden_files, current_file)

  if (length(golden_files) > 0) {
    prev <- readRDS(golden_files[1])

    # Gene overlap: at least 30% of previously mentioned genes still appear
    overlap <- intersect(findings$genes_mentioned, prev$genes_mentioned)
    overlap_pct <- length(overlap) / max(length(prev$genes_mentioned), 1)
    expect_true(overlap_pct >= 0.3,
      info = sprintf(
        "Gene overlap dropped to %.0f%% (threshold: 30%%)\nOverlap: %s\nNew: %s\nMissing: %s",
        overlap_pct * 100,
        paste(head(overlap, 15), collapse = ", "),
        paste(setdiff(findings$genes_mentioned, prev$genes_mentioned)[1:10], collapse = ", "),
        paste(setdiff(prev$genes_mentioned, findings$genes_mentioned)[1:10], collapse = ", ")))

    # Key protein coverage shouldn't drop by more than 2
    expect_true(findings$n_key_found >= prev$n_key_found - 2,
      info = sprintf("Key protein coverage dropped: %d -> %d",
        prev$n_key_found, findings$n_key_found))

    # Section coverage shouldn't regress
    expect_true(sum(findings$sections_found) >= sum(prev$sections_found) - 1,
      info = sprintf("Section coverage dropped: %d -> %d",
        sum(prev$sections_found), sum(findings$sections_found)))

    message(sprintf(paste0(
      "\n=== Claude Export Drift Report ===\n",
      "  Compared against: %s\n",
      "  Gene overlap: %.0f%% (%d shared / %d previous / %d current)\n",
      "  Key proteins: %d/%d (was %d/%d)\n",
      "  Stable markers: %d/%d (was %d/%d)\n",
      "  Sections: %d/%d (was %d/%d)\n",
      "  Word count: %d (was %d)\n"),
      basename(golden_files[1]),
      overlap_pct * 100, length(overlap),
      length(prev$genes_mentioned), length(findings$genes_mentioned),
      findings$n_key_found, length(KEY_PROTEINS),
      prev$n_key_found, length(KEY_PROTEINS),
      findings$n_stable_found, length(STABLE_MARKERS),
      prev$n_stable_found, length(STABLE_MARKERS),
      sum(findings$sections_found), length(findings$sections_found),
      sum(prev$sections_found), length(prev$sections_found),
      findings$word_count, prev$word_count
    ))
  } else {
    message(sprintf(paste0(
      "\n=== Claude Export Test — First Run (Baseline) ===\n",
      "  Key proteins found: %d/%d (%s)\n",
      "  Stable markers found: %d/%d (%s)\n",
      "  Sections: %d/%d\n",
      "  Genes mentioned: %d (%s)\n",
      "  Word count: %d\n",
      "  Golden baseline saved to: tests/golden/\n"),
      findings$n_key_found, length(KEY_PROTEINS),
      paste(KEY_PROTEINS[findings$key_proteins_found], collapse = ", "),
      findings$n_stable_found, length(STABLE_MARKERS),
      paste(STABLE_MARKERS[findings$stable_markers_found], collapse = ", "),
      sum(findings$sections_found), length(findings$sections_found),
      findings$n_genes, paste(head(findings$genes_mentioned, 20), collapse = ", "),
      findings$word_count
    ))
  }
})
