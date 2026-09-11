# ==============================================================================
#  HELPER FUNCTIONS — AI Integration (Google Gemini + OpenAI-compatible)
# ==============================================================================
#
#  Two providers are supported:
#    gemini        — Google's generativeLanguage API (File API for data upload)
#    openai_compat — any OpenAI-compatible /v1/chat/completions endpoint
#                    (vLLM, llama.cpp, LiteLLM, the UC Davis metabolomics
#                    gateway, ...). No File API: the data table is inlined.
#
#  Call sites MUST use the ask_ai_* dispatchers, never ask_gemini_* directly,
#  so that adding a third provider requires no edits outside this file
#  (CLAUDE.md architectural rule #3: concepts have one definition).
# ==============================================================================

# --- PROVIDER REGISTRY (the single definition) ---
ai_providers <- function() {
  list(
    gemini = list(
      label             = "Google Gemini (cloud)",
      destination       = "Google's Gemini API",
      key_label         = "Gemini API Key",
      key_placeholder   = "AIzaSy...",
      default_model     = "gemini-2.5-flash",
      base_url          = "https://generativelanguage.googleapis.com/v1beta",
      supports_file_api = TRUE,
      # Google responds quickly; no need for a long ceiling
      timeout_s         = 300,
      # ~1M token window: the full 800-protein table fits comfortably
      max_payload_chars = 400000
    ),
    openai_compat = list(
      label             = "OpenAI-compatible endpoint",
      destination       = "the OpenAI-compatible endpoint you configured",
      key_label         = "API Key",
      key_placeholder   = "sk-...",
      # Measured on the UC Davis gateway (2026-09-10): flash-next answered a
      # trivial prompt in 3s vs 35s / 67s / >90s for the other three, and
      # produced AI Summary output of equal quality. Fastest wins here because
      # every one of these is a thinking model and latency is queue-dominated.
      default_model     = "qwen3.8-flash-next",
      base_url          = "https://llm.metabolomics.us/v1",
      supports_file_api = FALSE,
      # Measured 2026-09-10: qwen3.8-27b exceeded 300s on a Data Chat request
      # and 90s on a two-token prompt. Latency here is queue-dominated, so the
      # ceiling is deliberately generous while the gateway is being tuned.
      timeout_s         = 900,
      # The smallest gateway models expose 65,536 tokens. DE tables tokenise at
      # ~1.01 chars/token (measured — prose is ~6.0), so chars ~= tokens here.
      # 40,000 leaves ~25k for the QC table, system prompt, reasoning and answer.
      max_payload_chars = 40000
    )
  )
}

# Measured chars-per-token for the numeric DE tables this app sends.
# Calibrated against the gateway's own prompt_tokens accounting on 2026-09-10:
# 50/200/800-protein tables gave 1.00 / 1.01 / 1.01. Do NOT substitute the
# usual ~4 chars/token rule of thumb — that is a PROSE figure (measured 6.0
# here) and under-counts a numeric table by roughly 6x.
ai_chars_per_token <- function() 1.01

# Read one field for a provider, falling back to gemini for unknown/NULL input
# so a stale saved session can never blank out a user-facing label.
ai_provider_field <- function(provider, field) {
  p <- ai_providers()
  key <- if (is.null(provider) || !nzchar(provider) || !(provider %in% names(p))) "gemini" else provider
  p[[key]][[field]]
}

# --- SHARED PAYLOAD FORMATTER ---
# Formats a DE result table for an LLM prompt. Used by BOTH providers so the
# model sees identical data regardless of where the request is sent.
#
# Three deliberate reductions (measured ~55% token saving vs. the raw CSV):
#   1. 3 significant figures  — 6dp of a log2 intensity is precision the
#                               measurement does not have, and costs tokens.
#   2. bare accessions        — "sp|P12345|ALBU_HUMAN" -> "P12345"
#   3. no per-sample columns  — limpa/limma already did the statistics; the
#                               model interprets logFC/adj.P.Val, it does not
#                               recompute from intensities. Group-level
#                               Mean_/SD_ summaries ARE kept.
#
# `max_chars` caps the result. Rows are dropped from the BOTTOM, which is the
# least-significant end because callers pass a topTable() ordered by p-value —
# so a trim costs the weakest evidence, never the strongest.
format_ai_table <- function(df, sig = 3, max_chars = Inf) {
  if (is.null(df) || nrow(df) == 0) {
    keep <- setdiff(names(df), character(0))
    return(paste(keep, collapse = "\t"))
  }

  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # Promote rownames to a Protein column if the caller did not supply one
  if (!"Protein" %in% names(df) && !is.null(rownames(df))) {
    df <- cbind(Protein = rownames(df), df, stringsAsFactors = FALSE)
  }

  # Drop per-sample expression columns; keep DE statistics, group summaries,
  # and the annotation/evidence columns DIA-NN already provides.
  #
  # Genes / Protein.Names come from the search FASTA via y_protein$genes. Sending
  # them is what stops a model inventing a protein identity — measured 2026-09-10,
  # where two models fabricated names for accessions whose real names were sitting
  # one column away. A join cannot be wrong; recall can.
  #
  # NPrec / PropObs are evidence strength: PropObs 0.083 says "seen in 8% of runs",
  # which makes a single-run artifact readable rather than something the model has
  # to infer from a zero standard deviation.
  stat_cols <- c("Protein", "Gene", "Genes", "Protein.Names",
                 "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
                 "NPrec", "PropObs")
  keep <- names(df)[names(df) %in% stat_cols | grepl("^(Mean|SD)_", names(df))]
  # Preserve a readable order rather than whatever the caller cbind()ed
  ord <- c(intersect(stat_cols, keep), setdiff(keep, stat_cols))
  df <- df[, ord, drop = FALSE]

  # Bare accessions: sp|P12345|ALBU_HUMAN -> P12345 (leave anything else alone)
  if ("Protein" %in% names(df)) {
    df$Protein <- vapply(as.character(df$Protein), function(x) {
      if (is.na(x)) return("")
      if (grepl("^[a-z]{2}\\|[^|]+\\|", x)) sub("^[a-z]{2}\\|([^|]+)\\|.*$", "\\1", x) else x
    }, character(1), USE.NAMES = FALSE)
  }

  # Round numerics; render NA as empty so the model does not read "NA" as data
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      # "g" picks the shorter of fixed/scientific per value: logFC stays
      # "-2.12", a p-value stays "1.2e-09". Forcing either wastes tokens.
      out <- formatC(signif(col, sig), digits = sig, format = "g")
      out <- trimws(out)
      out[is.na(col)] <- ""
      out
    } else {
      out <- as.character(col)
      out[is.na(out)] <- ""
      out
    }
  })

  header <- paste(names(df), collapse = "\t")
  rows   <- apply(df, 1, function(r) paste(r, collapse = "\t"))

  if (is.finite(max_chars)) {
    # Cumulative length including the newline that joins each row
    running <- nchar(header) + cumsum(nchar(rows) + 1L)
    n_keep  <- sum(running <= max_chars)
    if (n_keep < length(rows)) {
      message(sprintf("[DE-LIMP] AI payload trimmed: %d of %d rows kept to fit %d chars",
                      n_keep, length(rows), max_chars))
      rows <- if (n_keep > 0) rows[seq_len(n_keep)] else character(0)
    }
  }

  paste(c(header, rows), collapse = "\n")
}

# ==============================================================================
#  MEMORY — what the client resends, because the model remembers nothing
# ==============================================================================
#
#  No LLM has memory between API calls: every request is stateless and the model
#  sees only the `messages` array. What looks like memory in any chat product is
#  the client resending history. DE-LIMP kept `values$chat_history` for the UI
#  and the download button but never sent it, so a follow-up like "why is that?"
#  reached the model with no referent.
#
#  ONLY two sources feed memory here, both of which are incapable of fabricating:
#    - prior turns of THIS conversation (what was actually said)
#    - the user's own project name and notes from the activity log
#
#  Model conclusions are deliberately NOT persisted across sessions. A confident
#  fabrication ("P30838 is a hemoglobin" - measured, and wrong) would become
#  durable project knowledge and compound. That is CLAUDE.md rule #2 at scale.
# ==============================================================================

# Convert DE-LIMP's chat_history into OpenAI messages, newest-first-preserved.
ai_history_messages <- function(chat_history, max_turns = 6, max_chars = 6000) {
  if (is.null(chat_history) || length(chat_history) == 0) return(list())

  clean <- Filter(Negate(is.null), lapply(chat_history, function(m) {
    txt <- m$content %||% ""
    if (length(txt) != 1 || is.na(txt) || !nzchar(trimws(txt))) return(NULL)
    list(role = if (identical(m$role, "user")) "user" else "assistant",
         content = as.character(txt))
  }))
  if (length(clean) == 0) return(list())

  # Keep the most recent turns: drop from the front
  if (length(clean) > max_turns) clean <- clean[(length(clean) - max_turns + 1):length(clean)]

  # Then trim to the character budget, again dropping oldest first
  while (length(clean) > 1 &&
         sum(vapply(clean, function(x) nchar(x$content), numeric(1))) > max_chars) {
    clean <- clean[-1]
  }
  # A single surviving turn may still exceed the budget on its own
  if (length(clean) == 1 && nchar(clean[[1]]$content) > max_chars) {
    clean[[1]]$content <- substr(clean[[1]]$content, 1, max_chars)
  }
  clean
}

# The Gemini path takes one text blob, not a messages array.
ai_history_text <- function(chat_history, max_turns = 6, max_chars = 6000) {
  m <- ai_history_messages(chat_history, max_turns, max_chars)
  if (length(m) == 0) return("")
  paste0("\n--- EARLIER IN THIS CONVERSATION ---\n",
         paste(vapply(m, function(x)
           paste0(if (x$role == "user") "User: " else "Assistant: ", x$content),
           character(1)), collapse = "\n"),
         "\n--- END ---\n")
}

# Project context from the user's OWN notes. Emits nothing when there is nothing
# to say - never a placeholder, never an invented default (rule #2).
format_project_context <- function(project = NULL, notes = NULL, max_chars = 1500) {
  ok <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) && nzchar(trimws(as.character(x)))
  parts <- character(0)
  if (ok(project)) parts <- c(parts, paste0("Project: ", trimws(as.character(project))))
  if (ok(notes)) {
    n <- trimws(as.character(notes))
    if (nchar(n) > max_chars) n <- paste0(substr(n, 1, max_chars), " [truncated]")
    parts <- c(parts, paste0("Notes: ", n))
  }
  if (length(parts) == 0) return("")
  paste0("\n--- PROJECT CONTEXT (recorded by the user, not generated by AI) ---\n",
         paste(parts, collapse = "\n"), "\n--- END ---\n")
}

# ==============================================================================
#  OPENAI-COMPATIBLE PROVIDER (/v1/chat/completions)
# ==============================================================================

# Strip reasoning traces. Qwen3 and other thinking models served through vLLM
# emit <think>...</think> inline; some gateways return a separate
# reasoning_content field. Neither belongs in a user-facing summary.
strip_reasoning <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return("")
  txt <- gsub("(?s)<think>.*?</think>", "", txt, perl = TRUE)
  # An unterminated <think> means the model ran out of budget mid-thought
  txt <- gsub("(?s)<think>.*$", "", txt, perl = TRUE)
  trimws(txt)
}

# Turn one /chat/completions message object into the text a user should see.
#
# Reasoning models expose their chain-of-thought either inline as <think>...
# </think> or, on this gateway, as a separate `reasoning_content` field. Both
# are internal deliberation. Neither may reach a scientific summary — a reader
# would take "let me check the logFC column" as a finding. If the model spent
# its whole budget thinking and produced no answer, say so plainly rather than
# dressing the trace up as analysis (CLAUDE.md rule #2: never silently
# substitute a fabricated value for a missing one).
extract_openai_reply <- function(msg) {
  if (is.null(msg)) return("API Error: endpoint returned an empty response.")
  content <- strip_reasoning(msg$content %||% "")
  if (nzchar(content)) return(content)

  reasoning <- msg$reasoning_content %||% ""
  if (nzchar(reasoning)) {
    return(paste("API Error: the model used its entire output budget on internal",
                 "reasoning and returned no answer. Raise max_tokens, or choose a",
                 "model with a smaller reasoning overhead."))
  }
  "API Error: endpoint returned an empty response."
}

# max_tokens covers reasoning AND the answer on this gateway. Measured: 4096
# was entirely consumed by deepseek-v4-flash's internal reasoning on a Data Chat
# prompt, leaving an empty answer. 8192 leaves room for both.
ask_openai_compat <- function(system_prompt, user_prompt, api_key, model_name,
                              base_url = NULL, max_tokens = 8192,
                              timeout_s = NULL, history = list()) {
  timeout_s <- timeout_s %||% ai_provider_field("openai_compat", "timeout_s")
  base_url <- base_url %||% ai_provider_field("openai_compat", "base_url")
  base_url <- sub("/+$", "", base_url)

  body <- list(
    model = model_name,
    messages = c(
      list(list(role = "system", content = system_prompt)),
      history,                                    # prior turns, oldest first
      list(list(role = "user", content = user_prompt))
    ),
    temperature = 0.3,
    max_tokens  = max_tokens
  )

  req <- request(paste0(base_url, "/chat/completions")) %>%
    req_headers(
      "Authorization" = paste("Bearer", api_key),
      "Content-Type"  = "application/json"
    ) %>%
    req_body_json(body) %>%
    req_timeout(timeout_s)

  tryCatch({
    resp <- req_perform(req)
    payload <- resp_body_json(resp)
    if (length(payload$choices) == 0) return("API Error: endpoint returned no choices.")
    extract_openai_reply(payload$choices[[1]]$message)
  }, error = function(e) {
    err_msg <- if (!is.null(e$resp)) {
      tryCatch(resp_body_string(e$resp), error = function(z) e$message)
    } else e$message
    paste("API Error:", err_msg)
  })
}

list_openai_compat_models <- function(api_key, base_url = NULL) {
  base_url <- base_url %||% ai_provider_field("openai_compat", "base_url")
  base_url <- sub("/+$", "", base_url)
  req <- request(paste0(base_url, "/models")) %>%
    req_headers("Authorization" = paste("Bearer", api_key)) %>%
    req_timeout(60)
  tryCatch({
    data <- resp_body_json(req_perform(req))
    ids <- vapply(data$data, function(x) x$id %||% "", character(1))
    ids[nzchar(ids)]
  }, error = function(e) paste("Error listing models:", e$message))
}

# ==============================================================================
#  DISPATCHERS — the ONLY entry points call sites may use
# ==============================================================================

# Plain text prompt, no data table (AI Summary, Comparator hypothesis engine)
ask_ai_text <- function(user_query, api_key, model_name, provider = "gemini",
                        base_url = NULL, timeout_s = NULL, chat_history = NULL) {
  if (identical(provider, "openai_compat")) {
    ask_openai_compat(
      system_prompt = "You are a PhD-level expert in proteomics and systems biology.",
      user_prompt   = user_query,
      api_key = api_key, model_name = model_name, base_url = base_url,
      timeout_s = timeout_s, history = ai_history_messages(chat_history)
    )
  } else {
    ask_gemini_text_chat(paste0(ai_history_text(chat_history), user_query),
                         api_key, model_name)
  }
}

# Data-grounded chat. `data_table` is pre-formatted text from format_ai_table().
# Both providers inline it: at ~10k tokens there is nothing for Gemini's File
# API to solve, and one code path is easier to reason about than two.
ask_ai_data <- function(user_query, data_table, qc_df, api_key, model_name,
                        provider = "gemini", selected_ids = NULL, base_url = NULL,
                        timeout_s = NULL, chat_history = NULL,
                        project = NULL, notes = NULL) {

  qc_text <- if (is.null(qc_df)) "No QC data available." else
    paste(capture.output(write.table(qc_df, sep = "\t", row.names = FALSE, quote = FALSE)), collapse = "\n")

  selection_context <- ""
  if (!is.null(selected_ids) && length(selected_ids) > 0) {
    selection_context <- paste0(
      "\n!!! USER SELECTION ACTIVE !!!\n",
      "Focus analysis on these specific proteins:\n",
      paste(selected_ids, collapse = ", "), "\n"
    )
  }

  system_prompt <- paste0(
    "You are a PhD-level expert in proteomics. You have two data sources.\n\n",
    "SOURCE 1 - QC METRICS (tab-separated):\n",
    "Use the 'Group' column to compare technical quality (Precursors, MS1) between groups.\n",
    "--- START QC DATA ---\n", qc_text, "\n--- END QC DATA ---\n\n",
    "SOURCE 2 - DIFFERENTIAL EXPRESSION RESULTS (tab-separated):\n",
    "These are the output of limpa/limma. The statistics are already computed - ",
    "interpret them, do not recompute. logFC is log2 fold-change; adj.P.Val is ",
    "BH-adjusted; B is the log-odds of differential expression (higher = stronger ",
    "evidence, NOT a batch effect). Mean_/SD_ columns, when present, are per-group ",
    "summaries. Genes/Protein.Names are the identities from the search FASTA. ",
    "NPrec is the number of precursors and PropObs the proportion of runs in which ",
    "the protein was observed - a low PropObs means the result rests on very few ",
    "measurements and should be reported with that caveat.\n\n",
    "RULE ON EFFECT SIZE vs SIGNIFICANCE: 'most increased' and 'most decreased' ",
    "mean the LARGEST and SMALLEST logFC. That is a question about effect size, ",
    "NOT about significance - the protein with the best p-value or the highest B ",
    "is often a different one. If they differ, report both and label which is ",
    "which.\n",
    "RULE ON PROTEIN IDENTITY: name a protein ONLY using the Genes or ",
    "Protein.Names column supplied above. If those columns are absent or blank for ",
    "a row, refer to it by accession alone and say the identity was not supplied. ",
    "Never supply a protein name, family, or function from memory - a wrong ",
    "identity invalidates the interpretation built on it.\n",
    "--- START DE DATA ---\n", data_table, "\n--- END DE DATA ---\n\n",
    "BI-DIRECTIONAL CONTROL:\n",
    "1. If the user asks about 'selected proteins', see the USER SELECTION section.\n",
    "2. If you find interesting proteins, output their IDs at the end like this:\n",
    "   [[SELECT: P12345; P67890]]\n"
  )

  proj_context <- format_project_context(project, notes)

  if (identical(provider, "openai_compat")) {
    ask_openai_compat(
      system_prompt = paste0(system_prompt, proj_context, selection_context),
      user_prompt   = user_query,
      api_key = api_key, model_name = model_name, base_url = base_url,
      timeout_s = timeout_s, history = ai_history_messages(chat_history)
    )
  } else {
    ask_gemini_text_chat(
      paste0(system_prompt, proj_context, selection_context,
             ai_history_text(chat_history), "\n\nUser Question: ", user_query),
      api_key, model_name
    )
  }
}

# Model listing, dispatched by provider
list_ai_models <- function(api_key, provider = "gemini", base_url = NULL) {
  if (identical(provider, "openai_compat")) {
    list_openai_compat_models(api_key, base_url)
  } else {
    list_google_models(api_key)
  }
}

# --- CHECK AVAILABLE MODELS ---
list_google_models <- function(api_key) {
  req <- request("https://generativelanguage.googleapis.com/v1beta/models") %>%
    req_url_query(key = api_key)
  tryCatch({
    resp <- req_perform(req)
    data <- resp_body_json(resp)
    models <- sapply(data$models, function(x) x$name)
    models <- gsub("^models/", "", models)
    return(models)
  }, error = function(e) { return(paste("Error listing models:", e$message)) })
}

# --- FILE API UPLOADER (LEGACY — no longer on the Data Chat path) ---
# Kept for reference and any out-of-band use. Since v4.1.0 the DE table is
# trimmed by format_ai_table() to ~10k tokens and inlined for both providers,
# so there is nothing left for the File API to solve. Note the File API never
# reduced context cost — uploaded files are tokenised into the request exactly
# like inline text; it only avoided re-uploading bytes.
upload_csv_to_gemini <- function(df, api_key) {
  temp_file <- tempfile(fileext = ".csv")
  write.csv(df, temp_file, row.names = FALSE)
  file_size <- file.size(temp_file)

  req <- request("https://generativelanguage.googleapis.com/upload/v1beta/files") %>%
    req_url_query(key = api_key) %>%
    req_headers(
      "X-Goog-Upload-Protocol" = "raw",
      "X-Goog-Upload-Command" = "start, upload, finalize",
      "X-Goog-Upload-Header-Content-Length" = as.character(file_size),
      "X-Goog-Upload-Header-Content-Type" = "text/csv",
      "Content-Type" = "text/csv"
    ) %>%
    req_body_file(temp_file)

  resp <- req_perform(req)
  file_info <- resp_body_json(resp)
  return(file_info$file$uri)
}

# --- AI CHAT FUNCTION ---
ask_gemini_file_chat <- function(user_query, file_uri, qc_df, api_key, model_name, selected_ids = NULL) {

  qc_text <- "No QC Data Available"
  if(!is.null(qc_df)) {
    qc_text <- paste(capture.output(write.csv(qc_df, row.names=FALSE)), collapse="\n")
  }

  selection_context <- ""
  if (!is.null(selected_ids) && length(selected_ids) > 0) {
    selection_context <- paste0(
      "\n!!! URGENT: USER SELECTION ACTIVE !!!\n",
      "Focus analysis on these specific proteins:\n",
      paste(selected_ids, collapse=", "), "\n"
    )
  }

  system_instruction <- paste0(
    "You are a PhD-level expert in proteomics. ",
    "You have access to two data sources:\n",
    "SOURCE 1: QC METRICS (In Text Below)\n",
    "This table includes 'Group' columns. Use it to compare technical quality (Precursors, MS1) between experimental groups.\n",
    "--- START QC DATA ---\n",
    qc_text,
    "\n--- END QC DATA ---\n\n",
    "SOURCE 2: EXPRESSION DATA (In Uploaded File)\n",
    "Use this file to answer biological questions. It contains the Top 800 proteins.\n\n",
    "IMPORTANT: BI-DIRECTIONAL CONTROL.\n",
    "1. If the user asks about 'selected proteins', refer to the 'URGENT' section below.\n",
    "2. If you find interesting proteins, OUTPUT their IDs at the end like this:\n",
    "   [[SELECT: P12345; P67890]]\n"
  )

  base_url <- "https://generativelanguage.googleapis.com/v1beta/models/"
  clean_model <- gsub("^models/", "", model_name)
  full_url <- paste0(base_url, clean_model, ":generateContent")

  body <- list(contents = list(list(parts = list(
    list(text = paste0(system_instruction, selection_context, "\n\nUser Question: ", user_query)),
    list(file_data = list(file_uri = file_uri, mime_type = "text/csv"))
  ))))

  req <- request(full_url) %>%
    req_url_query(key = api_key) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_body_json(body)

  tryCatch({
    resp <- req_perform(req)
    return(resp_body_json(resp)$candidates[[1]]$content$parts[[1]]$text)
  }, error = function(e) {
    err_msg <- "Unknown Error"
    if (!is.null(e$resp)) { err_msg <- tryCatch(resp_body_string(e$resp), error = function(z) e$message) } else { err_msg <- e$message }
    return(paste("API Error:", err_msg))
  })
}

# --- AI TEXT CHAT FUNCTION ---
ask_gemini_text_chat <- function(user_query, api_key, model_name) {
  base_url <- "https://generativelanguage.googleapis.com/v1beta/models/"
  clean_model <- gsub("^models/", "", model_name)
  full_url <- paste0(base_url, clean_model, ":generateContent")

  body <- list(contents = list(list(parts = list(list(text = user_query)))))

  req <- request(full_url) %>%
    req_url_query(key = api_key) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_body_json(body)

  tryCatch({
    resp <- req_perform(req)
    return(resp_body_json(resp)$candidates[[1]]$content$parts[[1]]$text)
  }, error = function(e) {
    err_msg <- "Unknown Error"
    if (!is.null(e$resp)) { err_msg <- tryCatch(resp_body_string(e$resp), error = function(z) e$message) } else { err_msg <- e$message }
    return(paste("API Error:", err_msg))
  })
}
