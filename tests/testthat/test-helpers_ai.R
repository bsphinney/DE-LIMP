# Tests for the AI provider abstraction and payload formatter in R/helpers_ai.R
# Pure functions only — no HTTP requests are made.

# =============================================================================
# ai_providers() — the single provider registry (CLAUDE.md rule #3)
# =============================================================================

test_that("ai_providers returns both providers with the required fields", {
  p <- ai_providers()
  expect_true(all(c("gemini", "openai_compat") %in% names(p)))
  for (nm in names(p)) {
    expect_true(all(c("label", "destination", "key_label", "key_placeholder",
                      "default_model", "supports_file_api") %in% names(p[[nm]])),
                info = nm)
    expect_true(nzchar(p[[nm]]$destination), info = nm)
  }
})

test_that("only Gemini claims File API support", {
  p <- ai_providers()
  expect_true(p$gemini$supports_file_api)
  expect_false(p$openai_compat$supports_file_api)
})

test_that("ai_provider_field falls back to gemini for an unknown provider", {
  expect_equal(ai_provider_field("gemini", "default_model"),
               ai_providers()$gemini$default_model)
  expect_equal(ai_provider_field("nonsense-provider", "default_model"),
               ai_providers()$gemini$default_model)
  expect_equal(ai_provider_field(NULL, "default_model"),
               ai_providers()$gemini$default_model)
})

test_that("every provider names a concrete destination for the privacy modal", {
  # Rule #1: user-facing text must describe where data actually goes.
  for (nm in names(ai_providers())) {
    d <- ai_provider_field(nm, "destination")
    expect_type(d, "character")
    expect_gt(nchar(d), 3)
  }
})

# =============================================================================
# format_ai_table() — the shared payload formatter
# =============================================================================

make_de_df <- function(n = 5) {
  data.frame(
    Protein   = sprintf("sp|P%05d|GENE%d_HUMAN", seq_len(n), seq_len(n)),
    logFC     = seq(-2.123456, 2.123456, length.out = n),
    AveExpr   = seq(10.987654, 18.987654, length.out = n),
    t         = seq(-5.5555555, 5.5555555, length.out = n),
    P.Value   = seq(1e-9, 0.04, length.out = n),
    adj.P.Val = seq(1e-6, 0.049, length.out = n),
    B         = seq(-1.111111, 9.111111, length.out = n),
    stringsAsFactors = FALSE
  )
}

test_that("format_ai_table emits tab-separated text with a header", {
  out <- format_ai_table(make_de_df(3))
  lines <- strsplit(out, "\n", fixed = TRUE)[[1]]
  expect_equal(length(lines), 4L)             # header + 3 rows
  expect_true(grepl("\t", lines[1], fixed = TRUE))
  expect_false(grepl(",", lines[1], fixed = TRUE))
})

test_that("format_ai_table strips sp|...| accession decoration", {
  out <- format_ai_table(make_de_df(2))
  expect_false(grepl("sp|", out, fixed = TRUE))
  expect_true(grepl("P00001", out, fixed = TRUE))
})

test_that("format_ai_table leaves bare accessions and gene names untouched", {
  df <- make_de_df(2)
  df$Protein <- c("P12345", "ALBU_HUMAN")
  out <- format_ai_table(df)
  expect_true(grepl("P12345", out, fixed = TRUE))
  expect_true(grepl("ALBU_HUMAN", out, fixed = TRUE))
})

test_that("format_ai_table rounds to 3 significant figures", {
  out <- format_ai_table(make_de_df(3))
  # 6-decimal precision must not survive
  expect_false(grepl("2.123456", out, fixed = TRUE))
  expect_false(grepl("10.987654", out, fixed = TRUE))
})

test_that("format_ai_table renders ordinary magnitudes as plain decimals", {
  # Scientific notation on a logFC is both harder for the model to read and
  # MORE tokens than the plain decimal it replaces.
  df <- make_de_df(2)
  df$logFC <- c(-2.123456, 3.987654)
  df$AveExpr <- c(18.987654, 12.345678)
  out <- format_ai_table(df)
  expect_true(grepl("-2.12", out, fixed = TRUE))
  expect_true(grepl("3.99", out, fixed = TRUE))
  expect_true(grepl("19", out, fixed = TRUE))
  expect_false(grepl("e+00", out, fixed = TRUE))
  expect_false(grepl("e+01", out, fixed = TRUE))
})

test_that("format_ai_table keeps scientific notation for small p-values", {
  df <- make_de_df(2)
  df$adj.P.Val <- c(4.5e-7, 1.1e-3)
  out <- format_ai_table(df)
  expect_true(grepl("4.5e-07", out, fixed = TRUE) || grepl("4.5e-7", out, fixed = TRUE))
})

test_that("format_ai_table drops per-sample expression columns by default", {
  df <- make_de_df(3)
  df$Sample_1.d <- runif(3, 8, 24)
  df$Sample_2.d <- runif(3, 8, 24)
  out <- format_ai_table(df)
  expect_false(grepl("Sample_1.d", out, fixed = TRUE))
  expect_false(grepl("Sample_2.d", out, fixed = TRUE))
  expect_true(grepl("logFC", out, fixed = TRUE))
  expect_true(grepl("adj.P.Val", out, fixed = TRUE))
})

test_that("format_ai_table keeps group-level Mean_/SD_ summary columns", {
  df <- make_de_df(3)
  df$Mean_Control <- runif(3, 8, 24)
  df$SD_Control   <- runif(3, 0, 2)
  df$Sample_9.d   <- runif(3, 8, 24)
  out <- format_ai_table(df)
  expect_true(grepl("Mean_Control", out, fixed = TRUE))
  expect_true(grepl("SD_Control", out, fixed = TRUE))
  expect_false(grepl("Sample_9.d", out, fixed = TRUE))
})

test_that("format_ai_table is substantially smaller than the raw CSV it replaces", {
  df <- make_de_df(200)
  for (i in 1:12) df[[paste0("Sample_", i, ".d")]] <- runif(200, 8, 24)
  raw <- paste(capture.output(write.csv(df, row.names = FALSE)), collapse = "\n")
  out <- format_ai_table(df)
  expect_lt(nchar(out), nchar(raw) * 0.5)
})

# =============================================================================
# Token budget — MEASURED against the gateway: numeric DE tables tokenise at
# ~1.01 chars/token (English prose is ~6.0). A character budget is therefore a
# token budget for this payload, within a few percent.
# =============================================================================

test_that("ai_chars_per_token reflects the measured numeric-table ratio", {
  expect_lt(ai_chars_per_token(), 1.2)
  expect_gt(ai_chars_per_token(), 0.9)
})

test_that("format_ai_table trims rows to fit max_chars", {
  df <- make_de_df(800)
  out <- format_ai_table(df, max_chars = 5000)
  expect_lte(nchar(out), 5000)
  expect_gt(nchar(out), 3000)   # should use most of the budget, not bail out
})

test_that("format_ai_table keeps the header when trimming", {
  df <- make_de_df(800)
  lines <- strsplit(format_ai_table(df, max_chars = 3000), "\n", fixed = TRUE)[[1]]
  expect_true(grepl("Protein", lines[1], fixed = TRUE))
  expect_true(grepl("adj.P.Val", lines[1], fixed = TRUE))
  expect_gt(length(lines), 2L)
})

test_that("format_ai_table trims from the bottom, keeping most-significant rows", {
  df <- make_de_df(100)
  df$Protein <- sprintf("PROT%03d", seq_len(100))
  out <- format_ai_table(df, max_chars = 1500)
  expect_true(grepl("PROT001", out, fixed = TRUE))   # top of the table survives
  expect_false(grepl("PROT100", out, fixed = TRUE))  # bottom is dropped
})

test_that("format_ai_table does not trim when the table already fits", {
  df <- make_de_df(5)
  untrimmed <- format_ai_table(df)
  expect_equal(format_ai_table(df, max_chars = 100000), untrimmed)
})

test_that("format_ai_table with an impossibly small budget still returns a header", {
  df <- make_de_df(50)
  out <- format_ai_table(df, max_chars = 10)
  expect_true(grepl("Protein", out, fixed = TRUE))
})

test_that("providers declare a payload budget that fits their context window", {
  p <- ai_providers()
  # openai_compat must fit the smallest gateway model (65,536 tokens) with room
  # for the QC table, system prompt, reasoning and the answer.
  expect_lt(p$openai_compat$max_payload_chars, 65536 * 0.75)
  # Gemini's 1M window can take the full table
  expect_gt(p$gemini$max_payload_chars, 100000)
})

# =============================================================================
# Request timeout — configurable, because gateway latency varies wildly
# =============================================================================

test_that("providers declare a request timeout", {
  for (nm in names(ai_providers())) {
    t <- ai_provider_field(nm, "timeout_s")
    expect_type(t, "double")
    expect_gt(t, 0)
  }
})

test_that("the OpenAI-compatible timeout is generous enough for a slow gateway", {
  # Measured 2026-09-10: qwen3.8-27b exceeded 300s on a Data Chat request.
  expect_gte(ai_provider_field("openai_compat", "timeout_s"), 600)
})

test_that("ask_openai_compat accepts an explicit timeout override", {
  expect_true("timeout_s" %in% names(formals(ask_openai_compat)))
  expect_true("timeout_s" %in% names(formals(ask_ai_text)))
  expect_true("timeout_s" %in% names(formals(ask_ai_data)))
})

# =============================================================================
# Annotation columns — the fix for fabricated protein identities.
# DIA-NN already carries Genes / Protein.Names / NPrec / PropObs in
# y_protein$genes. Sending them means the model never has to recall an
# identity, which is where every measured hallucination came from.
# =============================================================================

test_that("format_ai_table keeps Genes and Protein.Names when supplied", {
  df <- make_de_df(3)
  df$Genes <- c("ALDH3A1", "DERPC", "TMSB4X")
  df$Protein.Names <- c("AL3A1_HUMAN", "DERPC_HUMAN", "TYB4_HUMAN")
  out <- format_ai_table(df)
  expect_true(grepl("Genes", out, fixed = TRUE))
  expect_true(grepl("ALDH3A1", out, fixed = TRUE))
  expect_true(grepl("AL3A1_HUMAN", out, fixed = TRUE))
})

test_that("format_ai_table keeps the evidence-strength columns", {
  # PropObs turns "this hit is a single-run artifact" from an inference into a
  # number the model can read directly.
  df <- make_de_df(2)
  df$NPrec   <- c(2L, 10L)
  df$PropObs <- c(0.0833, 0.6417)
  out <- format_ai_table(df)
  expect_true(grepl("NPrec", out, fixed = TRUE))
  expect_true(grepl("PropObs", out, fixed = TRUE))
  expect_true(grepl("0.0833", out, fixed = TRUE))
})

test_that("format_ai_table still drops per-sample columns when annotations present", {
  df <- make_de_df(2)
  df$Genes <- c("A", "B")
  df$Sample_1.d <- c(1, 2)
  out <- format_ai_table(df)
  expect_true(grepl("Genes", out, fixed = TRUE))
  expect_false(grepl("Sample_1.d", out, fixed = TRUE))
})

test_that("format_ai_table tolerates missing annotations without emitting NA", {
  df <- make_de_df(3)
  df$Genes <- c("ALDH3A1", NA_character_, "")
  out <- format_ai_table(df)
  expect_false(grepl("\tNA\t", out, fixed = TRUE))
  expect_true(grepl("ALDH3A1", out, fixed = TRUE))
})

# =============================================================================
# Conversation memory — prior turns threaded back into the request.
# No model has memory between calls; what looks like memory is the client
# resending history. DE-LIMP kept chat_history for display only.
# =============================================================================

mk_hist <- function(n) {
  h <- list()
  for (i in seq_len(n)) {
    h <- append(h, list(list(role = "user", content = paste0("question ", i))))
    h <- append(h, list(list(role = "ai",   content = paste0("answer ", i))))
  }
  h
}

test_that("ai_history_messages maps DE-LIMP roles to OpenAI roles", {
  m <- ai_history_messages(mk_hist(1))
  expect_equal(vapply(m, function(x) x$role, character(1)), c("user", "assistant"))
})

test_that("ai_history_messages keeps the MOST RECENT turns when capped", {
  m <- ai_history_messages(mk_hist(10), max_turns = 4)
  expect_equal(length(m), 4L)
  txt <- paste(vapply(m, function(x) x$content, character(1)), collapse = " ")
  expect_true(grepl("answer 10", txt, fixed = TRUE))   # newest survives
  expect_false(grepl("question 1 ", txt, fixed = TRUE)) # oldest dropped
})

test_that("ai_history_messages respects a character budget", {
  h <- list(list(role = "user", content = strrep("x", 5000)),
            list(role = "ai",   content = strrep("y", 5000)),
            list(role = "user", content = "short recent question"))
  m <- ai_history_messages(h, max_turns = 10, max_chars = 1000)
  total <- sum(vapply(m, function(x) nchar(x$content), numeric(1)))
  expect_lte(total, 1000)
  expect_true(grepl("short recent question",
                    paste(vapply(m, function(x) x$content, character(1)), collapse = " ")))
})

test_that("ai_history_messages returns an empty list for no history", {
  expect_equal(ai_history_messages(NULL), list())
  expect_equal(ai_history_messages(list()), list())
})

test_that("ai_history_messages drops blank and malformed entries", {
  h <- list(list(role = "user", content = ""), list(role = "ai", content = NA_character_),
            list(role = "user", content = "real question"))
  m <- ai_history_messages(h)
  expect_equal(length(m), 1L)
  expect_equal(m[[1]]$content, "real question")
})

test_that("ai_history_text renders history for the single-prompt Gemini path", {
  out <- ai_history_text(mk_hist(2))
  expect_true(grepl("question 2", out, fixed = TRUE))
  expect_true(grepl("answer 1", out, fixed = TRUE))
  expect_equal(ai_history_text(NULL), "")
})

# =============================================================================
# Project memory — user-authored notes only. NEVER model output (rule #2).
# =============================================================================

test_that("format_project_context renders project and notes when present", {
  out <- format_project_context(project = "Evosep pilot", notes = "Batch 2 rerun after column swap")
  expect_true(grepl("Evosep pilot", out, fixed = TRUE))
  expect_true(grepl("column swap", out, fixed = TRUE))
})

test_that("format_project_context emits nothing when there is nothing to say", {
  # Rule #2: never fabricate a default that reaches user-facing text.
  expect_equal(format_project_context(NULL, NULL), "")
  expect_equal(format_project_context("", ""), "")
  expect_equal(format_project_context(NA, NA), "")
})

test_that("format_project_context omits only the missing half", {
  a <- format_project_context(project = "Evosep pilot", notes = NULL)
  expect_true(grepl("Evosep pilot", a, fixed = TRUE))
  expect_false(grepl("Notes", a, fixed = TRUE))
  b <- format_project_context(project = NA, notes = "sample 7 degraded")
  expect_true(grepl("sample 7 degraded", b, fixed = TRUE))
  expect_false(grepl("Project", b, fixed = TRUE))
})

test_that("format_project_context marks the content as user-authored", {
  # The model must not confuse a human note with its own earlier conclusion.
  out <- format_project_context("P1", "some note")
  expect_match(out, "user|analyst|recorded by", ignore.case = TRUE)
})

test_that("format_project_context truncates an over-long note", {
  out <- format_project_context("P1", strrep("z", 5000), max_chars = 500)
  expect_lte(nchar(out), 700)
})

# =============================================================================
# Prompt safety rules. Each pins a fix for a failure measured against the
# UC Davis gateway on 2026-09-10/11 — see docs and CHANGELOG for the benchmark.
# =============================================================================

ai_data_prompt <- function() {
  # Capture the system prompt ask_ai_data builds, without making a request.
  # The gemini branch folds everything into one string, so intercept there.
  captured <- NULL
  orig <- get("ask_gemini_text_chat", envir = globalenv())
  on.exit(assign("ask_gemini_text_chat", orig, envir = globalenv()), add = TRUE)
  assign("ask_gemini_text_chat",
         function(prompt, api_key, model_name) { captured <<- prompt; "" },
         envir = globalenv())
  ask_ai_data("Q", "Protein\tlogFC\nP1\t1", NULL, "k", "m", provider = "gemini")
  captured
}

test_that("the data prompt forbids naming a protein from memory", {
  p <- ai_data_prompt()
  expect_match(p, "RULE ON PROTEIN IDENTITY")
  expect_match(p, "[Nn]ever supply a protein name.*from memory")
})

test_that("the data prompt defines most increased/decreased by fold-change", {
  # Measured: two models answered the most *significant* protein when asked for
  # the most decreased one. Both cited correct statistics for the wrong protein.
  p <- ai_data_prompt()
  expect_match(p, "LARGEST|largest")
  expect_match(p, "effect size|not about significance|NOT about significance")
})

test_that("the data prompt explains B correctly", {
  # Measured: one model described B as a 'batch effect estimate' and inverted it.
  p <- ai_data_prompt()
  expect_match(p, "log-odds")
  expect_match(p, "NOT a batch effect|not a batch effect")
})

test_that("the data prompt requires evidence strength to be cited", {
  p <- ai_data_prompt()
  expect_match(p, "NPrec")
  expect_match(p, "PropObs")
})

test_that("format_ai_table handles NA without emitting the string 'NA' as a value", {
  df <- make_de_df(3)
  df$logFC[2] <- NA_real_
  out <- format_ai_table(df)
  expect_false(grepl("\tNA\t", out, fixed = TRUE))
})

test_that("format_ai_table on an empty frame returns a header-only string", {
  df <- make_de_df(0)
  out <- format_ai_table(df)
  expect_type(out, "character")
  expect_equal(length(strsplit(out, "\n", fixed = TRUE)[[1]]), 1L)
})

# =============================================================================
# strip_reasoning() — thinking models must not leak traces into summaries
# =============================================================================

test_that("strip_reasoning removes a complete think block", {
  expect_equal(strip_reasoning("<think>weighing options</think>The answer."), "The answer.")
})

test_that("strip_reasoning removes a multi-line think block", {
  txt <- "<think>\nline one\nline two\n</think>\n\nFinal answer here."
  expect_equal(strip_reasoning(txt), "Final answer here.")
})

test_that("strip_reasoning truncates an unterminated think block", {
  # Model hit its token budget mid-thought — emit nothing rather than raw trace
  expect_equal(strip_reasoning("Preamble. <think>still thinking and then cut"), "Preamble.")
})

test_that("strip_reasoning leaves ordinary text untouched", {
  expect_equal(strip_reasoning("## Overview\nNo reasoning tags here."),
               "## Overview\nNo reasoning tags here.")
})

test_that("strip_reasoning handles NULL and empty input", {
  expect_equal(strip_reasoning(NULL), "")
  expect_equal(strip_reasoning(""), "")
})

test_that("strip_reasoning does not eat markdown that merely mentions think", {
  txt <- "The protein is thought to bind DNA."
  expect_equal(strip_reasoning(txt), txt)
})

# =============================================================================
# extract_openai_reply() — what the user actually sees
# =============================================================================

test_that("extract_openai_reply returns clean content", {
  expect_equal(extract_openai_reply(list(content = "## Overview\nAll good.")),
               "## Overview\nAll good.")
})

test_that("extract_openai_reply strips inline think tags from content", {
  expect_equal(extract_openai_reply(list(content = "<think>hmm</think>Answer.")), "Answer.")
})

test_that("extract_openai_reply ignores reasoning_content when content is present", {
  # The gateway's models return reasoning in a separate field; it is internal
  # chain-of-thought and must never be mixed into a scientific summary.
  out <- extract_openai_reply(list(content = "The answer.",
                                   reasoning_content = "Let me think step by step..."))
  expect_equal(out, "The answer.")
  expect_false(grepl("step by step", out))
})

test_that("extract_openai_reply does NOT pass off reasoning as the answer", {
  # Model burned its whole budget thinking. Returning the raw trace would
  # present internal deliberation to a scientist as if it were analysis.
  out <- extract_openai_reply(list(content = "",
                                   reasoning_content = "First I should check the logFC column..."))
  expect_false(grepl("First I should check", out))
  expect_match(out, "output budget|max_tokens", ignore.case = TRUE)
})

test_that("extract_openai_reply reports a wholly empty message", {
  out <- extract_openai_reply(list(content = "", reasoning_content = ""))
  expect_match(out, "empty", ignore.case = TRUE)
})

test_that("extract_openai_reply handles missing fields", {
  expect_match(extract_openai_reply(list()), "empty", ignore.case = TRUE)
  expect_match(extract_openai_reply(NULL), "empty", ignore.case = TRUE)
})

test_that("format_ai_table preserves row order (rank by significance)", {
  df <- make_de_df(4)
  df$Protein <- c("AAA", "BBB", "CCC", "DDD")
  lines <- strsplit(format_ai_table(df), "\n", fixed = TRUE)[[1]]
  expect_true(startsWith(lines[2], "AAA"))
  expect_true(startsWith(lines[5], "DDD"))
})
