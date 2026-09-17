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

test_that("format_ai_table drops per-group Mean_/SD_ columns (documented as not sent)", {
  # v4.1.0 review finding I/L: the CHANGELOG and privacy modal say per-group
  # Mean/SD are not sent. The formatter is an allowlist, so they cannot leak in
  # via a caller that happens to carry them.
  df <- make_de_df(3)
  df$Mean_Control <- runif(3, 8, 24)
  df$SD_Control   <- runif(3, 0, 2)
  df$Sample_9.d   <- runif(3, 8, 24)
  out <- format_ai_table(df)
  expect_false(grepl("Mean_Control", out, fixed = TRUE))
  expect_false(grepl("SD_Control", out, fixed = TRUE))
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

test_that("format_ai_table keeps Genes but drops Protein.Names", {
  # Finding I: topTable() carries fit$genes columns, so Protein.Names reached the
  # model although the comments and CHANGELOG said it was dropped.
  df <- make_de_df(3)
  df$Genes <- c("ALDH3A1", "DERPC", "TMSB4X")
  df$Protein.Names <- c("AL3A1_HUMAN", "DERPC_HUMAN", "TYB4_HUMAN")
  out <- format_ai_table(df)
  expect_true(grepl("Genes", out, fixed = TRUE))
  expect_true(grepl("ALDH3A1", out, fixed = TRUE))
  expect_false(grepl("Protein.Names", out, fixed = TRUE))
  expect_false(grepl("AL3A1_HUMAN", out, fixed = TRUE))
})

test_that("format_ai_table keeps evidence columns only when the pipeline declares them", {
  # PropObs turns "this hit rests on few observations" from an inference into a
  # number the model can read directly — but which columns exist is the
  # pipeline's business (rule #1), passed in via extra_cols.
  df <- make_de_df(2)
  df$NPrec   <- c(2L, 10L)
  df$PropObs <- c(0.0833, 0.6417)
  out <- format_ai_table(df, extra_cols = c("NPrec", "PropObs"))
  expect_true(grepl("NPrec", out, fixed = TRUE))
  expect_true(grepl("PropObs", out, fixed = TRUE))
  expect_true(grepl("0.0833", out, fixed = TRUE))
  expect_false(grepl("NPrec", format_ai_table(df), fixed = TRUE))
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
            list(role = "user", content = "short recent question"),
            list(role = "ai",   content = "short recent answer"))
  m <- ai_history_messages(h, max_turns = 10, max_chars = 1000)
  total <- sum(vapply(m, function(x) nchar(x$content), numeric(1)))
  expect_lte(total, 1000)
  expect_true(grepl("short recent question",
                    paste(vapply(m, function(x) x$content, character(1)), collapse = " ")))
})

test_that("ai_history_messages shortens a single over-budget exchange instead of dropping it", {
  h <- list(list(role = "user", content = strrep("x", 5000)),
            list(role = "ai",   content = strrep("y", 5000)))
  m <- ai_history_messages(h, max_turns = 10, max_chars = 1000)
  expect_equal(vapply(m, function(x) x$role, character(1)), c("user", "assistant"))
  expect_lte(sum(vapply(m, function(x) nchar(x$content), numeric(1))), 1000)
})

test_that("ai_history_messages returns an empty list for no history", {
  expect_equal(ai_history_messages(NULL), list())
  expect_equal(ai_history_messages(list()), list())
})

test_that("ai_history_messages drops blank and malformed entries", {
  h <- list(list(role = "user", content = ""), list(role = "ai", content = NA_character_),
            list(role = "user", content = "real question"),
            list(role = "ai", content = "real answer"))
  m <- ai_history_messages(h)
  expect_equal(length(m), 2L)
  expect_equal(m[[1]]$content, "real question")
  expect_equal(m[[2]]$content, "real answer")
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

ai_data_prompt <- function(data_table = "Protein\tlogFC\nP1\t1", ...) {
  # Capture the system prompt ask_ai_data builds, without making a request.
  # The gemini branch folds everything into one string, so intercept there.
  captured <- NULL
  orig <- get("ask_gemini_text_chat", envir = globalenv())
  on.exit(assign("ask_gemini_text_chat", orig, envir = globalenv()), add = TRUE)
  assign("ask_gemini_text_chat",
         function(prompt, api_key, model_name, ...) { captured <<- prompt; "" },
         envir = globalenv())
  ask_ai_data("Q", data_table, NULL, "k", "m", provider = "gemini", ...)
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

test_that("the data prompt describes evidence columns when the pipeline supplies them", {
  ev <- dpc_pipeline_descriptor()$evidence_columns
  p <- ai_data_prompt("Protein\tlogFC\tNPrec\tPropObs\nP1\t1\t2\t0.5",
                      evidence = ev,
                      evidence_guidance = dpc_pipeline_descriptor()$evidence_guidance)
  expect_match(p, "NPrec is the number of precursors")
  expect_match(p, "PropObs is the proportion of this protein's precursor-by-run measurements")
  expect_match(p, "PropObs < 0.5")
})

test_that("the data prompt never describes NPrec/PropObs when they are not sent (MaxLFQ)", {
  # Finding G: the prompt described NPrec/PropObs even under MaxLFQ, where
  # the pipeline does not produce them.
  p <- ai_data_prompt("Protein\tlogFC\nP1\t1",
                      evidence = maxlfq_pipeline_descriptor()$evidence_columns)
  expect_false(grepl("NPrec", p, fixed = TRUE))
  expect_false(grepl("PropObs", p, fixed = TRUE))
  expect_match(p, "No per-protein measurement-depth statistics")
  # declared by the pipeline but absent from the table: still not described
  p2 <- ai_data_prompt("Protein\tlogFC\nP1\t1", evidence = dpc_pipeline_descriptor()$evidence_columns,
                       evidence_guidance = dpc_pipeline_descriptor()$evidence_guidance)
  expect_false(grepl("PropObs", p2, fixed = TRUE))
})

# Source-level guard: no DE-LIMP prompt may invite the model to identify a
# protein from memory. Measured 2026-09-10: that phrasing produced ALDH3A1
# reported as a haemoglobin, with invented interpretation built on top.
# This scans the real files, so a new prompt anywhere inherits the rule.
test_that("no prompt asks the model to recognise proteins from memory", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  banned <- c("note any you recognize", "where you recognize the gene name",
              "any you recognise")
  hits <- character(0)
  for (f in files) {
    src <- readLines(f, warn = FALSE)
    # ignore comment lines — they may legitimately quote the removed phrasing
    src <- src[!grepl("^\\s*#", src)]
    for (b in banned) {
      w <- grep(b, src, fixed = TRUE)
      if (length(w)) hits <- c(hits, paste0(basename(f), ": ", b))
    }
  }
  expect_equal(hits, character(0))
})

test_that("prompts do not ask for biology the model merely knows", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
  hits <- character(0)
  for (f in files) {
    src <- readLines(f, warn = FALSE)
    src <- src[!grepl("^\\s*#", src)]
    w <- grep("discuss their known biological functions", src, fixed = TRUE)
    if (length(w)) hits <- c(hits, basename(f))
  }
  expect_equal(hits, character(0))
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

test_that("strip_reasoning removes a trace whose opening <think> came from the template", {
  # Finding M: some chat templates inject <think> into the prompt, so the reply
  # holds only the reasoning and a lone closing tag.
  expect_equal(strip_reasoning("Let me check the logFC column first.\n</think>\n\n## Overview\nAnswer."),
               "## Overview\nAnswer.")
  expect_equal(strip_reasoning("trace</think>Answer <think>cut off"), "Answer")
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

# =============================================================================
# v4.1.0 review fixes. Each block names the finding it pins.
# No test here makes a network request: every path below returns before
# req_perform(), or uses an injected resolver / a hand-built response object.
# =============================================================================

fake_resolver <- function(ips) function(host) ips
pub  <- ai_deployment_policy(TRUE)
priv <- ai_deployment_policy(FALSE)

# --- A: credentials never cross providers ------------------------------------

gemini_like_key <- paste0("AIza", strrep("B", 35))

test_that("A: a Gemini-format key is recognised as belonging to another provider", {
  expect_true(ai_key_belongs_to_other_provider("openai_compat", gemini_like_key))
  expect_false(ai_key_belongs_to_other_provider("gemini", gemini_like_key))
  expect_false(ai_key_belongs_to_other_provider("openai_compat", "sk-local-gateway-key"))
  expect_false(ai_key_belongs_to_other_provider("openai_compat", ""))
  expect_false(ai_key_belongs_to_other_provider("openai_compat", NULL))
})

test_that("A: the OpenAI-compatible path refuses to send a Gemini key", {
  out <- ask_openai_compat("sys", "q", gemini_like_key, "m",
                           base_url = "https://llm.metabolomics.us/v1", policy = priv)
  expect_true(is_ai_error(out))
  expect_match(out, "Gemini key")
  expect_false(grepl(gemini_like_key, out, fixed = TRUE))
  expect_true(is_ai_error(list_openai_compat_models(gemini_like_key, policy = priv)))
})

test_that("A: the provider-switch observer clears the key (source guard)", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  src <- paste(readLines(file.path(root, "R", "server_ai.R"), warn = FALSE), collapse = "\n")
  block <- regmatches(src, regexpr("(?s)observeEvent\\(input\\$ai_provider,.*?ignoreInit = TRUE\\)", src, perl = TRUE))
  expect_length(block, 1)
  expect_match(block, "\"user_api_key\"")
  expect_match(block, "value = \"\"", fixed = TRUE)
})

# --- B: endpoint validation (SSRF) -------------------------------------------

test_that("B: IP classification covers loopback, private, link-local and metadata ranges", {
  expect_equal(ai_ip_class("8.8.8.8"), "public")
  expect_equal(ai_ip_class("2606:4700::1111"), "public")
  for (ip in c("127.0.0.1", "10.0.0.5", "172.16.0.1", "172.31.255.255", "192.168.1.1",
               "100.64.0.1", "100.100.100.200", "::1", "fd00:ec2::254", "::ffff:10.0.0.1"))
    expect_equal(ai_ip_class(ip), "private", info = ip)
  for (ip in c("169.254.169.254", "0.0.0.0", "224.0.0.1", "255.255.255.255", "192.0.0.192",
               "fe80::1", "::", "ff02::1", "::ffff:169.254.169.254", "64:ff9b::a9fe:a9fe",
               "not-an-ip", "999.1.1.1"))
    expect_equal(ai_ip_class(ip), "blocked", info = ip)
  expect_equal(ai_ip_class("172.32.0.1"), "public")
  expect_equal(ai_ip_class("::ffff:8.8.8.8"), "public")
})

test_that("B: a public deployment refuses internal endpoints, whatever the URL says", {
  internal <- list(
    list(url = "https://internal.example.org/v1", ips = "10.0.0.8"),
    list(url = "https://metadata.example.org/v1", ips = "169.254.169.254"),
    list(url = "https://localhost/v1", ips = "127.0.0.1"),
    list(url = "https://mixed.example.org/v1", ips = c("8.8.8.8", "192.168.0.10")),  # any bad address fails
    list(url = "https://127.0.0.1/v1", ips = character(0)),                          # literal IP
    list(url = "https://[::1]:8443/v1", ips = character(0)),
    list(url = "https://2130706433/v1", ips = "127.0.0.1")                            # decimal-encoded loopback
  )
  for (case in internal) {
    v <- ai_validate_endpoint(case$url, policy = pub, resolver = fake_resolver(case$ips))
    expect_false(v$ok, info = case$url)
  }
})

test_that("B: an unresolvable host and an internal host give the same public message", {
  a <- ai_validate_endpoint("https://nope.example.org/v1", policy = pub, resolver = fake_resolver(NULL))
  b <- ai_validate_endpoint("https://nope.example.org/v1", policy = pub, resolver = fake_resolver("10.1.1.1"))
  expect_false(a$ok); expect_false(b$ok)
  expect_identical(a$error, b$error)
})

test_that("B: https to a public host is accepted and empty input means the provider default", {
  v <- ai_validate_endpoint("https://api.example.org/v1/", policy = pub, resolver = fake_resolver("93.184.216.34"))
  expect_true(v$ok)
  expect_equal(v$url, "https://api.example.org/v1")
  expect_equal(v$port, 443L)
  d <- ai_validate_endpoint("", policy = pub, resolver = fake_resolver("93.184.216.34"))
  expect_true(d$ok)
  expect_equal(d$url, ai_provider_field("openai_compat", "base_url"))
  expect_true(ai_validate_endpoint(NULL, policy = pub, resolver = fake_resolver("93.184.216.34"))$ok)
})

test_that("B: plain http is refused except for localhost on a local install", {
  expect_false(ai_validate_endpoint("http://api.example.org/v1", policy = pub,
                                    resolver = fake_resolver("93.184.216.34"))$ok)
  expect_false(ai_validate_endpoint("http://api.example.org/v1", policy = priv,
                                    resolver = fake_resolver("93.184.216.34"))$ok)
  expect_false(ai_validate_endpoint("http://localhost:8000/v1", policy = pub,
                                    resolver = fake_resolver("127.0.0.1"))$ok)
  expect_true(ai_validate_endpoint("http://localhost:8000/v1", policy = priv,
                                   resolver = fake_resolver("127.0.0.1"))$ok)
  expect_false(ai_validate_endpoint("http://192.168.1.20:8000/v1", policy = priv,
                                    resolver = fake_resolver("192.168.1.20"))$ok)
})

test_that("B: a local install may use a private https endpoint but never metadata", {
  expect_true(ai_validate_endpoint("https://vllm.lab.local/v1", policy = priv,
                                   resolver = fake_resolver("192.168.1.20"))$ok)
  expect_false(ai_validate_endpoint("https://meta.example/v1", policy = priv,
                                    resolver = fake_resolver("169.254.169.254"))$ok)
})

test_that("B: URLs with credentials, queries, other schemes or junk are refused", {
  r <- fake_resolver("93.184.216.34")
  for (u in c("https://user:pw@api.example.org/v1", "https://api.example.org/v1?x=1",
              "https://api.example.org/v1#frag", "ftp://api.example.org/v1", "file:///etc/passwd",
              "gopher://api.example.org/", "api.example.org/v1", "https://api.example.org/v 1",
              "https://api.example.org:99999/v1"))
    expect_false(ai_validate_endpoint(u, policy = priv, resolver = r)$ok, info = u)
})

test_that("B: the request is pinned to the validated addresses and never follows redirects", {
  v <- ai_validate_endpoint("https://api.example.org/v1", policy = pub,
                            resolver = fake_resolver(c("93.184.216.34", "2606:2800:220:1::1")))
  req <- ai_harden_request(httr2::request("https://api.example.org/v1/chat/completions"), v)
  expect_equal(req$options$followlocation, 0L)
  expect_equal(req$options$resolve, "api.example.org:443:93.184.216.34,[2606:2800:220:1::1]")
  # no endpoint (Gemini's fixed host): still no redirects
  expect_equal(ai_harden_request(httr2::request("https://x.example"))$options$followlocation, 0L)
})

test_that("B: ask_openai_compat returns an error for a refused endpoint without a request", {
  out <- ask_openai_compat("sys", "q", "sk-test-key", "m", base_url = "https://127.0.0.1/v1", policy = pub)
  expect_true(is_ai_error(out))
  expect_match(out, "not allowed")
  out2 <- list_openai_compat_models("sk-test-key", base_url = "http://10.0.0.1/v1", policy = pub)
  expect_true(is_ai_error(out2))
})

test_that("B: upstream error bodies are reduced to a short sanitised status line", {
  key <- "sk-very-secret-key-123"
  body <- paste0('{"error":{"message":"Invalid key ', key, ' <script>alert(1)</script>', strrep("z", 500), '"}}')
  resp <- httr2::response(status_code = 401, headers = list(`Content-Type` = "application/json"),
                          body = charToRaw(body))
  cnd <- structure(class = c("httr2_http_401", "httr2_http", "httr2_error", "error", "condition"),
                   list(message = "HTTP 401 Unauthorized.", call = NULL, resp = resp))
  out <- ai_error_from_condition(cnd, secrets = key)
  expect_true(is_ai_error(out))
  expect_match(out, "HTTP 401 Unauthorized", fixed = TRUE)
  expect_false(grepl(key, out, fixed = TRUE))
  expect_false(grepl("<script>", out, fixed = TRUE))
  expect_lt(nchar(out), 300)

  html <- httr2::response(status_code = 500, headers = list(`Content-Type` = "text/html"),
                          body = charToRaw(paste0("<html>internal admin panel ", strrep("x", 5000), "</html>")))
  cnd2 <- structure(class = c("httr2_http_500", "httr2_http", "httr2_error", "error", "condition"),
                    list(message = "HTTP 500.", call = NULL, resp = html))
  out2 <- ai_error_from_condition(cnd2)
  expect_false(grepl("admin panel", out2, fixed = TRUE))
  expect_lt(nchar(out2), 80)
})

# --- C: project notes only for the loaded dataset and user -------------------

mk_log <- function() data.frame(
  id = 1:6,
  event_type = c("search_submitted", "analysis_completed", "data_loaded",
                 "search_submitted", "analysis_completed", "search_submitted"),
  user = c("alice", "alice", "alice", "bob", NA, "alice"),
  output_dir = c("/q/alice/run1", "/q/alice/run1/", "/q/alice/run1",
                 "/q/alice/run1", "/q/alice/run1", "/q/alice/run2"),
  project = c("Liver", "", NA, "BOB PROJECT", "", "Kidney"),
  notes = c("column swapped before batch 2", "", activity_note_job_loaded("run1", "123"),
            "bob's private note", "", "someone else's newest row"),
  stringsAsFactors = FALSE
)

test_that("C: the LAST row of the log is never used just because it is last", {
  ctx <- activity_project_context(mk_log(), "/q/alice/run1", current_users = "alice")
  expect_equal(ctx$project, "Liver")
  expect_equal(ctx$notes, "column swapped before batch 2")
  expect_false(identical(ctx$notes, "someone else's newest row"))
})

test_that("C: another user's row for the same folder is excluded", {
  ctx <- activity_project_context(mk_log(), "/q/alice/run1", current_users = "carol")
  expect_null(ctx$project)      # alice's rows excluded; the NA-user row has no project
  expect_null(ctx$notes)
  ctx_bob <- activity_project_context(mk_log(), "/q/alice/run1", current_users = "bob")
  expect_equal(ctx_bob$notes, "bob's private note")
})

test_that("C: nothing is sent when no dataset folder is known or nothing matches", {
  none <- list(project = NULL, notes = NULL)
  expect_equal(activity_project_context(mk_log(), NULL, "alice"), none)
  expect_equal(activity_project_context(mk_log(), "", "alice"), none)
  expect_equal(activity_project_context(mk_log(), NA_character_, "alice"), none)
  expect_equal(activity_project_context(mk_log(), "/q/other", "alice"), none)
  expect_equal(activity_project_context(data.frame(), "/q/alice/run1", "alice"), none)
  expect_equal(activity_project_context(NULL, "/q/alice/run1", "alice"), none)
})

test_that("C: app-written placeholder notes are never forwarded as the user's", {
  expect_true(all(activity_note_is_auto(c(ACTIVITY_NOTE_SESSION_RESTORED, ACTIVITY_NOTE_BACKFILLED,
                                          activity_note_job_loaded("My search", "4411")))))
  expect_false(any(activity_note_is_auto(c("Job: rerun with MBR on", "sample 7 degraded", NA))))
  log <- mk_log()[3, , drop = FALSE]
  expect_null(activity_project_context(log, "/q/alice/run1", "alice")$notes)
})

test_that("C: a public deployment never sends activity-log notes", {
  expect_false(ai_deployment_policy(TRUE)$send_activity_notes)
  expect_true(ai_deployment_policy(FALSE)$send_activity_notes)
  root <- normalizePath(file.path(getwd(), "..", ".."))
  src <- paste(readLines(file.path(root, "R", "server_ai.R"), warn = FALSE), collapse = "\n")
  expect_match(src, "if (!isTRUE(ai_policy$send_activity_notes)) return(none)", fixed = TRUE)
})

test_that("C: the server matches on the loaded dataset, not the newest row (source guard)", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  src <- paste(readLines(file.path(root, "R", "server_ai.R"), warn = FALSE), collapse = "\n")
  expect_false(grepl("log[nrow(log), , drop = FALSE]", src, fixed = TRUE))
  expect_match(src, "activity_project_context(log, od)", fixed = TRUE)
})

# --- D: timeouts applied everywhere, capped on a public deployment ------------

test_that("D: a public deployment caps every request timeout", {
  expect_equal(ai_request_timeout("openai_compat", NULL, pub), AI_PUBLIC_MAX_TIMEOUT_S)
  expect_equal(ai_request_timeout("gemini", NULL, pub), min(300, AI_PUBLIC_MAX_TIMEOUT_S))
  expect_equal(ai_request_timeout("openai_compat", 5000, pub), AI_PUBLIC_MAX_TIMEOUT_S)
  expect_equal(ai_request_timeout("openai_compat", 30, pub), 30)
  expect_lte(AI_PUBLIC_MAX_TIMEOUT_S, 300)
})

test_that("D: a local install keeps the provider's own ceiling", {
  expect_equal(ai_request_timeout("openai_compat", NULL, priv), 900)
  expect_equal(ai_request_timeout("gemini", NULL, priv), 300)
  expect_equal(ai_request_timeout("gemini", "junk", priv), 300)
})

test_that("D: the default policy is the strict (public) one", {
  expect_true(ai_deployment_policy()$public)
  expect_true(ai_deployment_policy(NULL)$public)
  expect_false(ai_deployment_policy(FALSE)$public)
})

test_that("D: every helper that performs a request sets a timeout", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  exprs <- parse(file.path(root, "R", "helpers_ai.R"))
  offenders <- character(0)
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
        identical(e[[3]][[1]], as.name("function"))) {
      body_txt <- paste(deparse(e[[3]]), collapse = "\n")
      if (grepl("req_perform(", body_txt, fixed = TRUE) && !grepl("req_timeout(", body_txt, fixed = TRUE))
        offenders <- c(offenders, as.character(e[[2]]))
    }
  }
  expect_equal(offenders, character(0))
})

# --- E: errors are detectable and never results ------------------------------

test_that("E: is_ai_error recognises every failure shape", {
  expect_true(is_ai_error(ai_error("boom")))
  expect_true(is_ai_error("  API Error: x"))
  expect_true(is_ai_error(NULL))
  expect_true(is_ai_error(NA_character_))
  expect_true(is_ai_error(""))
  expect_true(is_ai_error(c("a", "b")))
  expect_false(is_ai_error("## Overview\nAll good."))
})

test_that("E: extract_gemini_reply returns an error for blocked or empty answers", {
  expect_true(is_ai_error(extract_gemini_reply(list(candidates = list(),
                                                    promptFeedback = list(blockReason = "SAFETY")))))
  expect_true(is_ai_error(extract_gemini_reply(list(candidates = list(list(
    content = list(parts = list()), finishReason = "MAX_TOKENS"))))))
  ok <- extract_gemini_reply(list(candidates = list(list(content = list(parts = list(
    list(text = "internal plan", thought = TRUE), list(text = "The answer."), list(text = " More.")))))))
  expect_equal(ok, "The answer. More.")
})

test_that("E: no call site stores or renders an AI reply without checking it (source guard)", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  ai <- paste(readLines(file.path(root, "R", "server_ai.R"), warn = FALSE), collapse = "\n")
  expect_match(ai, "if (is_ai_error(ai_summary))", fixed = TRUE)
  comp <- paste(readLines(file.path(root, "R", "server_comparator.R"), warn = FALSE), collapse = "\n")
  expect_match(comp, "if (is_ai_error(response))", fixed = TRUE)
})

# --- F / G / J: AI Summary payload -------------------------------------------

mk_tt <- function() data.frame(
  Protein.Group = c("P1", "P2", "P3", "Cont_P00761", "P5"),
  Gene          = c("AAA", "BBB", "", "TRYP", "EEE"),
  logFC         = c(0.4, -0.3, 5.2, 1.1, -4.8),
  adj.P.Val     = c(1e-9, 1e-8, 0.03, 0.001, 0.04),
  NPrec         = c(12, 8, 1, 3, 2),
  PropObs       = c(0.9, 0.8, 0.1, 0.7, 0.4),
  stringsAsFactors = FALSE
)

test_that("F: the summary lists largest effect sizes separately from significance", {
  txt <- ai_contrast_summary_text(mk_tt(), "A - B", top_n = 2, top_fc_n = 1)
  sig_block <- sub("Largest INCREASES.*$", "", txt)
  up_block <- sub("^.*Largest INCREASES", "", sub("Largest DECREASES.*$", "", txt))
  down_block <- sub("^.*Largest DECREASES", "", txt)
  expect_match(sig_block, "P1")                 # most significant
  expect_false(grepl("P3", sig_block))          # biggest increase is NOT among the top-2 by p
  expect_match(up_block, "P3")                  # ...but it is the largest increase
  expect_match(down_block, "P5")
  expect_match(txt, "Significant proteins (adj.P.Val < 0.05): 5 (3 up, 2 down)", fixed = TRUE)
  expect_match(txt, "EFFECT SIZE")
})

test_that("G: evidence columns appear in the summary only when passed", {
  with_ev <- ai_contrast_summary_text(mk_tt(), "A - B", evidence_cols = c("NPrec", "PropObs"))
  expect_match(with_ev, "NPrec")
  without <- ai_contrast_summary_text(mk_tt(), "A - B")
  expect_false(grepl("NPrec", without, fixed = TRUE))
})

test_that("J: contaminants are flagged from the single definition, not described in prose", {
  txt <- ai_contrast_summary_text(mk_tt(), "A - B")
  expect_match(txt, "Contaminant")
  df <- ai_flag_contaminants(data.frame(Protein = c("P1", "Cont_P2", "sp|CON__P3|X", "contam_P4",
                                                    "cRAP-P5", "P6"), stringsAsFactors = FALSE))
  expect_equal(df$Contaminant, c("", "yes", "yes", "yes", "yes", ""))
  clean <- ai_flag_contaminants(data.frame(Protein = c("P1", "P2"), stringsAsFactors = FALSE))
  expect_false("Contaminant" %in% names(clean))
})

test_that("J: no prompt defines contaminants by an accession pattern (source guard)", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  hits <- character(0)
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) {
    src <- readLines(f, warn = FALSE)
    src <- src[!grepl("^\\s*#", src)]
    if (any(grepl("beginning Cont_", src, fixed = TRUE))) hits <- c(hits, basename(f))
  }
  expect_equal(hits, character(0))
})

test_that("G: the AI Summary prompt describes evidence only for pipelines that have it", {
  ctx <- list(n_contrasts = 1, contrast_text = "x", cross_text = "y", stable_prots_text = "z")
  dpc <- build_ai_summary_prompt(ctx, dpc_pipeline_descriptor()$evidence_columns,
                                 dpc_pipeline_descriptor()$evidence_guidance, "DPC-Quant + limma (limpa)")
  expect_match(dpc, "NPrec is the number of precursors")
  expect_match(dpc, "PropObs < 0.5")
  mx <- build_ai_summary_prompt(ctx, maxlfq_pipeline_descriptor()$evidence_columns, "",
                                "MaxLFQ + limma (Moschem 2025)")
  expect_false(grepl("NPrec", mx, fixed = TRUE))
  expect_false(grepl("PropObs", mx, fixed = TRUE))
  expect_match(mx, "No per-protein measurement-depth statistics")
  expect_match(mx, "MaxLFQ + limma", fixed = TRUE)
  # the rules that were already pinned still hold
  expect_match(dpc, "Never from memory", fixed = TRUE)
  expect_false(grepl("where you recognize the gene name", dpc, fixed = TRUE))
  expect_match(dpc, "largest and smallest logFC")
})

test_that("G: pipeline_evidence_columns reads the descriptor and the columns actually present", {
  dpc_y <- list(genes = data.frame(Genes = "A", NPrec = 2, PropObs = 0.5),
                other = list(descriptor = dpc_pipeline_descriptor()))
  expect_equal(names(pipeline_evidence_columns(dpc_y)), c("NPrec", "PropObs"))
  expect_match(pipeline_evidence_guidance(dpc_y), "PropObs")

  mx_y <- list(genes = data.frame(Protein.Group = "P1", Genes = "A"),
               other = list(descriptor = maxlfq_pipeline_descriptor()))
  expect_length(pipeline_evidence_columns(mx_y), 0)
  expect_equal(pipeline_evidence_guidance(mx_y), "")

  legacy <- dpc_pipeline_descriptor(); legacy$evidence_columns <- NULL; legacy$evidence_guidance <- NULL
  legacy_y <- list(genes = data.frame(NPrec = 1, PropObs = 1), other = list(descriptor = legacy))
  expect_equal(names(pipeline_evidence_columns(legacy_y)), c("NPrec", "PropObs"))

  single <- list(genes = data.frame(PropObs = 1), other = list(descriptor = dpc_pipeline_descriptor()))
  expect_equal(names(pipeline_evidence_columns(single)), "PropObs")   # dpcQuantByRow drops NPrec
  expect_equal(pipeline_evidence_guidance(single), "")                 # guidance names NPrec

  unknown <- list(genes = data.frame(NPrec = 1), other = list(descriptor = list(pipeline_id = "third")))
  expect_length(pipeline_evidence_columns(unknown), 0)
  expect_length(pipeline_evidence_columns(NULL), 0)
})

# --- H: selected proteins are never trimmed -----------------------------------

test_that("H: selected rows are listed first and survive trimming", {
  df <- make_de_df(400)
  df$Protein <- sprintf("PROT%03d", seq_len(400))
  out <- format_ai_table(df, max_chars = 2000, pin = c("PROT400", "PROT399"))
  lines <- strsplit(out, "\n", fixed = TRUE)[[1]]
  expect_true(startsWith(lines[2], "PROT399"))
  expect_true(startsWith(lines[3], "PROT400"))
  expect_true(grepl("PROT001", out, fixed = TRUE))
  expect_false(grepl("PROT200", out, fixed = TRUE))
  expect_lte(nchar(out), 2000)
})

test_that("H: selected rows are kept even when they alone exceed the budget", {
  df <- make_de_df(50)
  df$Protein <- sprintf("PROT%03d", seq_len(50))
  out <- format_ai_table(df, max_chars = 10, pin = "PROT050")
  expect_true(grepl("PROT050", out, fixed = TRUE))
})

test_that("H: the selection list in the prompt uses the same IDs as the table", {
  p <- ai_data_prompt("Protein\tlogFC\nP12345\t1", selected_ids = "sp|P12345|ALBU_HUMAN")
  expect_match(p, "P12345")
  expect_false(grepl("sp|P12345|", p, fixed = TRUE))
  expect_match(p, "listed first")
})

# --- K: exports name the provider and model that ran -------------------------

test_that("K: ai_source_label names provider, model and (for gateways) the host", {
  g <- ai_source_label("gemini", "gemini-2.5-flash")
  expect_match(g, "Google Gemini", fixed = TRUE)
  expect_match(g, "gemini-2.5-flash", fixed = TRUE)
  o <- ai_source_label("openai_compat", "qwen3.8-flash-next", "https://llm.example.org/v1")
  expect_match(o, "OpenAI-compatible endpoint at llm.example.org", fixed = TRUE)
  expect_match(o, "qwen3.8-flash-next", fixed = TRUE)
  expect_match(ai_source_label("gemini", ""), "not recorded")
})

test_that("K: the chat transcript attributes each turn correctly", {
  h <- list(list(role = "user", content = "q1"),
            list(role = "ai", content = "a1", source = "OpenAI-compatible endpoint at h, model m"),
            list(role = "ai", content = ai_error("boom"), source = "Google Gemini (cloud), model g", error = TRUE),
            list(role = "ai", content = "Please load data and run analysis first.", notice = TRUE),
            list(role = "ai", content = "legacy turn"))
  out <- ai_chat_transcript(h)
  expect_match(out[1], "^YOU: q1")
  expect_match(out[2], "^AI \\(OpenAI-compatible endpoint at h, model m\\): a1")
  expect_match(out[3], "^ERROR \\(Google Gemini")
  expect_match(out[4], "^DE-LIMP: ")
  expect_match(out[5], "not recorded")
  expect_false(any(grepl("GEMINI:", out, fixed = TRUE)))
})

test_that("K: no export hardcodes Gemini as the author of AI text (source guard)", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  banned <- c("\"GEMINI: \"", "## Gemini Analysis", "GEMINI PRE-ANALYSIS", "<b>Gemini Analysis</b>")
  hits <- character(0)
  for (f in list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)) {
    src <- readLines(f, warn = FALSE)
    src <- src[!grepl("^\\s*#", src)]
    for (b in banned) if (any(grepl(b, src, fixed = TRUE))) hits <- c(hits, paste0(basename(f), ": ", b))
  }
  expect_equal(hits, character(0))
})

# --- M: see strip_reasoning tests above --------------------------------------

# --- N: [[SELECT: ...]] maps back to real row names ---------------------------

test_that("N: ai_short_protein_id is the single shortening used for the payload", {
  expect_equal(ai_short_protein_id(c("sp|P12345|ALBU_HUMAN", "P99999", "tr|A0A0|X_Y;sp|Q1|Z", NA)),
               c("P12345", "P99999", "A0A0", ""))
})

test_that("N: IDs the model returns map back to the full row names", {
  full <- c("sp|P12345|ALBU_HUMAN", "P67890;Q11111", "sp|O00001|AAA_HUMAN;sp|O00002|BBB_HUMAN", "P22222")
  m <- ai_match_protein_ids(c("P12345", "Q11111", "O00002", "P22222", "NOPE1"), full)
  expect_setequal(m$matched, full)
  expect_equal(m$unmatched, "NOPE1")
  none <- ai_match_protein_ids("ALDH3A1", full)
  expect_length(none$matched, 0)
  expect_equal(none$unmatched, "ALDH3A1")
})

test_that("N: the SELECT directive is parsed out of the reply", {
  r <- ai_parse_select_directive("Look at these.\n[[SELECT: P12345; P67890, P12345]]")
  expect_equal(r$text, "Look at these.")
  expect_equal(r$ids, c("P12345", "P67890"))
  expect_equal(ai_parse_select_directive("No directive.")$ids, character(0))
})

test_that("N: the chat only claims a plot update when something matched (source guard)", {
  root <- normalizePath(file.path(getwd(), "..", ".."))
  src <- paste(readLines(file.path(root, "R", "server_ai.R"), warn = FALSE), collapse = "\n")
  expect_match(src, "if (length(m$matched) > 0) {", fixed = TRUE)
  expect_false(grepl("values$plot_selected_proteins <- trimws(id_vec)", src, fixed = TRUE))
})

# --- O: history is valid for strict-alternation chat templates ----------------

roles_of <- function(m) vapply(m, function(x) x$role, character(1))
alternates <- function(r) length(r) == 0 ||
  (r[1] == "user" && r[length(r)] == "assistant" && all(r[-1] != r[-length(r)]))

test_that("O: error turns and notices are never re-sent as history", {
  h <- list(list(role = "user", content = "q1"),
            list(role = "ai", content = ai_error("HTTP 500")),
            list(role = "user", content = "q2"),
            list(role = "ai", content = "a2"),
            list(role = "user", content = "q3"),
            list(role = "ai", content = "some text", error = TRUE),
            list(role = "user", content = "q4"),
            list(role = "ai", content = "Please load data and run analysis first.", notice = TRUE))
  m <- ai_history_messages(h, max_turns = 20)
  txt <- vapply(m, function(x) x$content, character(1))
  expect_false(any(grepl("API Error", txt, fixed = TRUE)))
  expect_false(any(txt %in% c("some text", "Please load data and run analysis first.")))
  expect_equal(txt, c("q2", "a2"))
  expect_true(alternates(roles_of(m)))
})

test_that("O: history never starts with the assistant or ends with the user", {
  h <- list(list(role = "ai", content = "orphan answer"),
            list(role = "user", content = "q1"),
            list(role = "user", content = "q1 again"),
            list(role = "ai", content = "a1"),
            list(role = "user", content = "unanswered"))
  m <- ai_history_messages(h)
  expect_equal(vapply(m, function(x) x$content, character(1)), c("q1 again", "a1"))
  expect_true(alternates(roles_of(m)))
})

test_that("O: alternation holds for arbitrary histories, caps and budgets", {
  set.seed(42)
  for (i in 1:200) {
    n <- sample(0:12, 1)
    h <- lapply(seq_len(n), function(j) {
      kind <- sample(c("user", "ai", "err", "notice", "blank"), 1, prob = c(.4, .35, .1, .05, .1))
      switch(kind,
             user = list(role = "user", content = strrep("u", sample(1:3000, 1))),
             ai = list(role = "ai", content = strrep("a", sample(1:3000, 1))),
             err = list(role = "ai", content = ai_error("x")),
             notice = list(role = "ai", content = "n", notice = TRUE),
             blank = list(role = sample(c("user", "ai"), 1), content = ""))
    })
    m <- ai_history_messages(h, max_turns = sample(1:8, 1), max_chars = sample(c(500, 6000), 1))
    expect_true(alternates(roles_of(m)), info = paste("iteration", i))
  }
})
