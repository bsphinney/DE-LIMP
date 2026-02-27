server_ai <- function(input, output, session, values) {

  # --- AI SUMMARY (Data Overview Tab) — Analyzes ALL contrasts ---
  observeEvent(input$generate_ai_summary_overview, {
    req(values$fit, values$y_protein, input$user_api_key)

    withProgress(message = "Generating AI Summary...", value = 0, {
      incProgress(0.1, detail = "Gathering DE data across all comparisons...")

      all_contrasts <- colnames(values$fit$contrasts)
      n_contrasts <- length(all_contrasts)

      # Scale top-N per contrast to manage token budget
      top_n <- if (n_contrasts <= 3) 30 else if (n_contrasts <= 6) 20 else 10

      # --- Gene mapping (done once, shared across contrasts) ---
      first_tt <- topTable(values$fit, coef = all_contrasts[1], number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(first_tt)) first_tt <- first_tt %>% rownames_to_column("Protein.Group")
      first_tt$Accession <- str_split_fixed(first_tt$Protein.Group, "[; ]", 2)[,1]
      org_db_name <- detect_organism_db(first_tt$Protein.Group)

      id_map <- tryCatch({
        if (!requireNamespace(org_db_name, quietly = TRUE)) BiocManager::install(org_db_name, ask = FALSE)
        library(org_db_name, character.only = TRUE)
        db_obj <- get(org_db_name)
        AnnotationDbi::select(db_obj, keys = first_tt$Accession, columns = c("SYMBOL"), keytype = "UNIPROT") %>%
          dplyr::rename(Accession = UNIPROT, Gene = SYMBOL) %>% distinct(Accession, .keep_all = TRUE)
      }, error = function(e) data.frame(Accession = first_tt$Accession, Gene = first_tt$Accession))

      # --- Per-contrast DE summaries ---
      contrast_texts <- list()
      all_sig_proteins <- list()  # track which proteins are sig in which contrasts

      for (i in seq_along(all_contrasts)) {
        cname <- all_contrasts[i]
        incProgress(0.1 + 0.3 * (i / n_contrasts), detail = paste0("Analyzing: ", cname, "..."))

        tt <- topTable(values$fit, coef = cname, number = Inf) %>% as.data.frame()
        if (!"Protein.Group" %in% colnames(tt)) tt <- tt %>% rownames_to_column("Protein.Group")
        tt$Accession <- str_split_fixed(tt$Protein.Group, "[; ]", 2)[,1]
        tt <- left_join(tt, id_map, by = "Accession")
        tt$Gene[is.na(tt$Gene)] <- tt$Accession[is.na(tt$Gene)]

        sig <- tt %>% filter(adj.P.Val < 0.05)
        n_up <- sum(sig$logFC > 0)
        n_down <- sum(sig$logFC < 0)

        # Track per-protein significance across contrasts
        if (nrow(sig) > 0) {
          for (pid in sig$Protein.Group) {
            if (is.null(all_sig_proteins[[pid]])) all_sig_proteins[[pid]] <- list()
            row <- sig[sig$Protein.Group == pid, ]
            all_sig_proteins[[pid]][[cname]] <- list(
              gene = row$Gene[1], logFC = round(row$logFC[1], 3), pval = round(row$adj.P.Val[1], 4)
            )
          }
        }

        top_hits <- sig %>% arrange(adj.P.Val) %>% head(top_n) %>%
          dplyr::select(Gene, logFC, adj.P.Val) %>%
          mutate(across(where(is.numeric), ~round(.x, 3)))

        top_text <- paste(capture.output(print(as.data.frame(top_hits))), collapse = "\n")

        contrast_texts[[cname]] <- paste0(
          "### ", cname, "\n",
          "Significant proteins: ", nrow(sig), " (", n_up, " up, ", n_down, " down)\n\n",
          "Top ", min(top_n, nrow(top_hits)), " by significance:\n", top_text
        )
      }

      incProgress(0.5, detail = "Identifying cross-comparison biomarkers...")

      # --- Cross-contrast proteins (significant in >= 2 comparisons) ---
      multi_contrast <- names(all_sig_proteins)[sapply(all_sig_proteins, length) >= 2]
      cross_text <- if (length(multi_contrast) > 0) {
        cross_df <- do.call(rbind, lapply(head(multi_contrast, 10), function(pid) {
          info <- all_sig_proteins[[pid]]
          gene <- info[[1]]$gene
          contrasts_str <- paste(names(info), collapse = ", ")
          fc_str <- paste(sapply(names(info), function(cn) {
            paste0(cn, ": ", sprintf("%+.2f", info[[cn]]$logFC))
          }), collapse = "; ")
          data.frame(Gene = gene, N_Comparisons = length(info), Contrasts = contrasts_str, LogFC = fc_str)
        }))
        cross_df <- cross_df[order(-cross_df$N_Comparisons), ]
        paste(capture.output(print(as.data.frame(cross_df), row.names = FALSE)), collapse = "\n")
      } else {
        "No proteins were significant in more than one comparison."
      }

      incProgress(0.6, detail = "Computing stable biomarkers...")

      # --- Stable biomarkers (lowest CV across all significant proteins) ---
      stable_prots_text <- tryCatch({
        all_sig_pids <- names(all_sig_proteins)
        valid_pids <- intersect(all_sig_pids, rownames(values$y_protein$E))
        if (length(valid_pids) == 0) return("No significant proteins to assess for stability.")

        raw_exprs <- values$y_protein$E[valid_pids, , drop = FALSE]
        linear_exprs <- 2^raw_exprs
        cv_list <- list()

        for (g in unique(values$metadata$Group)) {
          if (g == "") next
          files_in_group <- values$metadata$File.Name[values$metadata$Group == g]
          group_cols <- intersect(colnames(linear_exprs), files_in_group)
          if (length(group_cols) > 1) {
            group_data <- linear_exprs[, group_cols, drop = FALSE]
            cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
          }
        }

        if (length(cv_list) == 0) return("Could not calculate CVs (not enough replicates).")

        cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
        cv_df$Avg_CV <- rowMeans(cv_df[, grep("^CV_", colnames(cv_df)), drop = FALSE], na.rm = TRUE)

        # Get gene names and which contrasts each is in
        stable_df <- cv_df %>% arrange(Avg_CV) %>% head(5)
        stable_df$Gene <- sapply(stable_df$Protein.Group, function(pid) {
          info <- all_sig_proteins[[pid]]
          if (!is.null(info)) info[[1]]$gene else pid
        })
        stable_df$Significant_In <- sapply(stable_df$Protein.Group, function(pid) {
          info <- all_sig_proteins[[pid]]
          if (!is.null(info)) paste(names(info), collapse = "; ") else ""
        })

        out <- stable_df %>% dplyr::select(Gene, Avg_CV, Significant_In) %>%
          mutate(Avg_CV = round(Avg_CV, 2))
        paste(capture.output(print(as.data.frame(out), row.names = FALSE)), collapse = "\n")
      }, error = function(e) "Could not calculate stable proteins.")

      incProgress(0.7, detail = "Constructing prompt...")

      # --- Build the full prompt ---
      system_prompt <- paste0(
        "You are a senior proteomics and systems biology consultant. Write a comprehensive ",
        "analysis of the differential expression results across ALL comparisons below.\n\n",
        "Structure your response with these markdown sections:\n\n",
        "## Overview\n",
        "Number of comparisons analyzed, total significant proteins per comparison (up/down split). ",
        "Overall assessment of the experiment's quality and scope.\n\n",
        "## Key Findings Per Comparison\n",
        "For each comparison: highlight the top upregulated and downregulated proteins by fold-change ",
        "(use gene names). Note any comparison with unusually few or many significant hits.\n\n",
        "## Cross-Comparison Biomarkers\n",
        "Proteins significant in multiple comparisons are highest-confidence candidates. ",
        "Discuss consistency of direction (always up, always down, or mixed across comparisons).\n\n",
        "## High-Confidence Biomarker Insights\n",
        "For the most stable proteins (lowest coefficient of variation): discuss their known biological functions, ",
        "pathway involvement, and disease associations where you recognize the gene name. ",
        "Assess their potential as reliable biomarkers based on the combination of low CV, ",
        "significant p-value, and meaningful fold-change.\n\n",
        "## Biological Interpretation\n",
        "Suggest what biological processes or pathways may be affected based on the protein lists. ",
        "Note any well-known protein families, complexes, or signaling cascades represented. ",
        "If the data suggests a clear biological narrative, describe it.\n\n",
        "Use markdown formatting with headers. Be scientific but accessible."
      )

      all_contrast_text <- paste(contrast_texts, collapse = "\n\n")

      final_prompt <- paste0(
        system_prompt,
        "\n\n--- DATA FOR ANALYSIS ---\n\n",
        "Number of comparisons: ", n_contrasts, "\n\n",
        all_contrast_text, "\n\n",
        "--- CROSS-COMPARISON PROTEINS (significant in >= 2 comparisons) ---\n",
        cross_text, "\n\n",
        "--- MOST STABLE SIGNIFICANT PROTEINS (lowest CV across replicates) ---\n",
        stable_prots_text
      )

      message(sprintf("[DE-LIMP] AI Summary prompt: %d characters, %d contrasts", nchar(final_prompt), n_contrasts))

      incProgress(0.8, detail = "Asking AI...")
      ai_summary <- ask_gemini_text_chat(final_prompt, input$user_api_key, input$model_name)

      # Store for export
      values$ai_summary_text <- ai_summary

      # Render the summary to the output area
      output$ai_summary_output <- renderUI({
        div(style = "background-color: #ffffff; padding: 20px; border: 1px solid #dee2e6; border-radius: 8px;",
          tags$h5(icon("check-circle"), " Analysis Complete", style = "color: #28a745; margin-bottom: 15px;"),
          HTML(markdown::markdownToHTML(text = ai_summary, fragment.only = TRUE))
        )
      })
    })
  })

  # --- AI Summary HTML Export ---
  output$download_ai_summary_html <- downloadHandler(
    filename = function() {
      paste0("AI_Analysis_Report_", format(Sys.time(), "%Y%m%d_%H%M"), ".html")
    },
    content = function(file) {
      req(values$ai_summary_text)
      summary_html <- markdown::markdownToHTML(text = values$ai_summary_text, fragment.only = TRUE)

      # Build a styled standalone HTML report
      html_doc <- paste0(
        '<!DOCTYPE html>\n<html lang="en">\n<head>\n',
        '<meta charset="UTF-8">\n',
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">\n',
        '<title>DE-LIMP AI Analysis Report</title>\n',
        '<style>\n',
        '  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif; ',
        'max-width: 800px; margin: 40px auto; padding: 0 20px; color: #333; line-height: 1.7; }\n',
        '  .header { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); ',
        'color: white; padding: 25px 30px; border-radius: 10px; margin-bottom: 30px; }\n',
        '  .header h1 { margin: 0 0 5px 0; font-size: 1.5em; }\n',
        '  .header p { margin: 0; opacity: 0.9; font-size: 0.9em; }\n',
        '  .content { background: #fff; padding: 25px; border: 1px solid #dee2e6; border-radius: 8px; }\n',
        '  h2 { color: #2c3e50; border-bottom: 2px solid #667eea; padding-bottom: 8px; margin-top: 30px; }\n',
        '  h3 { color: #34495e; }\n',
        '  table { border-collapse: collapse; width: 100%; margin: 15px 0; }\n',
        '  th, td { border: 1px solid #dee2e6; padding: 8px 12px; text-align: left; }\n',
        '  th { background-color: #f1f3f5; font-weight: 600; }\n',
        '  tr:nth-child(even) { background-color: #f8f9fa; }\n',
        '  code { background: #f1f3f5; padding: 2px 6px; border-radius: 3px; font-size: 0.9em; }\n',
        '  strong { color: #2c3e50; }\n',
        '  .footer { text-align: center; margin-top: 30px; padding: 15px; color: #6c757d; font-size: 0.85em; border-top: 1px solid #dee2e6; }\n',
        '  @media print { body { margin: 20px; } .header { -webkit-print-color-adjust: exact; print-color-adjust: exact; } }\n',
        '</style>\n</head>\n<body>\n',
        '<div class="header">\n',
        '  <h1>DE-LIMP AI Analysis Report</h1>\n',
        '  <p>Generated: ', format(Sys.time(), "%B %d, %Y at %I:%M %p"), '</p>\n',
        '</div>\n',
        '<div class="content">\n',
        summary_html, '\n',
        '</div>\n',
        '<div class="footer">\n',
        '  Generated by <strong>DE-LIMP Proteomics</strong> &mdash; ',
        '  <a href="https://github.com/bsphinney/DE-LIMP">github.com/bsphinney/DE-LIMP</a>\n',
        '</div>\n',
        '</body>\n</html>'
      )

      writeLines(html_doc, file)
    }
  )

  # --- AI Summary Info Modal ---
  observeEvent(input$ai_summary_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About AI Summary"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("How it works"),
        p("The AI Summary analyzes ", strong("all comparisons"), " in your experiment at once, not just the currently selected contrast. ",
          "It identifies the top differentially expressed proteins per comparison, finds proteins that are significant across multiple ",
          "comparisons (cross-comparison biomarkers), and highlights the most reproducibly measured proteins (lowest CV) as high-confidence candidates."),
        p("The AI then provides biological interpretation, discussing known functions, pathway involvement, and disease associations for the top biomarkers."),
        tags$h6("What data is sent to Google Gemini"),
        tags$ul(
          tags$li("Top significant proteins per comparison (gene names, log2 fold-changes, adjusted p-values)"),
          tags$li("Proteins significant across multiple comparisons with their fold-changes"),
          tags$li("Most stable significant proteins (lowest coefficient of variation across replicates)"),
          tags$li("Number of comparisons and significance counts (up/down)")
        ),
        tags$h6("What is NOT sent"),
        tags$ul(
          tags$li("Raw expression values or intensity data"),
          tags$li("Individual sample names or file paths"),
          tags$li("Metadata details (groups, batches, covariates)"),
          tags$li("QC statistics or run-level information")
        ),
        tags$h6("Privacy"),
        p("Data is sent to Google's Gemini API and processed according to Google's API terms of service. ",
          "No data is stored permanently by this app \u2014 uploaded files are deleted when your session ends."),
        tags$h6("API key"),
        p("You need a Google Gemini API key (enter in the sidebar). Get one free at ",
          tags$a(href = "https://aistudio.google.com/apikey", target = "_blank", "Google AI Studio"), ".")
      )
    ))
  })

  # --- Data Chat Info Modal ---
  observeEvent(input$data_chat_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About AI Analysis"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("How it works"),
        p("Data Chat uses the Google Gemini API to provide AI-powered analysis of your proteomics data. ",
          "Your QC statistics and the top 800 proteins are uploaded to Gemini for context-aware responses."),
        tags$h6("What data is sent"),
        tags$ul(
          tags$li("QC statistics (precursor counts, protein counts, MS1 signal per sample)"),
          tags$li("Top 800 DE proteins with fold-changes and p-values"),
          tags$li("Your chat messages")
        ),
        tags$h6("Privacy"),
        p("Data is sent to Google's Gemini API. It is processed according to Google's API terms of service. ",
          "No data is stored permanently by this app \u2014 uploaded files are deleted when your session ends."),
        tags$h6("Plot selection integration"),
        p("If you select proteins in the volcano plot or results table, the chat knows about your selection. ",
          "The AI can also suggest proteins to highlight \u2014 look for the ",
          tags$em("'I have updated your plots'"), " message after AI responses."),
        tags$h6("API key"),
        p("You need a Google Gemini API key (enter in the sidebar). Get one free at ",
          tags$a(href = "https://aistudio.google.com/apikey", target = "_blank", "Google AI Studio"), ".")
      )
    ))
  })

  observeEvent(input$check_models, { if (nchar(input$user_api_key) < 10) { showNotification("Please enter a valid API Key first.", type="error"); return() }; withProgress(message = "Checking Google Models...", { models <- list_google_models(input$user_api_key); if (length(models) > 0 && !grepl("Error", models[1])) { showModal(modalDialog(title = "Available Models for Your Key", p("Copy one of these into the Model Name box:"), tags$textarea(paste(models, collapse="\n"), rows=10, style="width:100%;"), easyClose = TRUE)) } else { showNotification(paste("Failed to list models:", models), type="error") } }) })
  output$chat_selection_indicator <- renderText({ if (!is.null(values$plot_selected_proteins)) { paste("\u2705 Current Selection:", length(values$plot_selected_proteins), "Proteins from Plots.") } else { "\u2139\ufe0f No proteins selected in plots." } })

  observeEvent(input$summarize_data, {
    req(input$user_api_key)
    auto_prompt <- "Analyze this dataset. Identify key quality control issues (if any) by looking at the Group QC stats. Then, summarize the main biological findings from the expression data, focusing on the most significantly differentially expressed proteins."
    values$chat_history <- append(values$chat_history, list(list(role = "user", content = "(Auto-Query: Summarize & Analyze)")))
    withProgress(message = "Auto-Analyzing Dataset...", {
      if (!is.null(values$fit) && !is.null(values$y_protein)) {
        # Scale protein count to stay within Gemini's token limit (~1M tokens)
        n_samples <- ncol(values$y_protein$E)
        n_max <- if (n_samples > 200) 100 else if (n_samples > 100) 200 else if (n_samples > 50) 400 else 800
        message(sprintf("[DE-LIMP] AI data: %d proteins x %d samples (scaled from 800)", n_max, n_samples))

        df_de <- topTable(values$fit, coef=input$contrast_selector, number=n_max)

        # For large datasets (>100 samples), send group-level summary stats
        # instead of per-sample expression to stay within token limits
        if (n_samples > 100 && !is.null(values$metadata)) {
          exprs_mat <- values$y_protein$E[rownames(df_de), , drop = FALSE]
          meta <- values$metadata[values$metadata$Group != "", ]
          group_stats <- do.call(cbind, lapply(unique(meta$Group), function(g) {
            cols <- intersect(meta$File.Name[meta$Group == g], colnames(exprs_mat))
            if (length(cols) == 0) return(NULL)
            data.frame(
              setNames(list(
                rowMeans(exprs_mat[, cols, drop = FALSE], na.rm = TRUE),
                apply(exprs_mat[, cols, drop = FALSE], 1, sd, na.rm = TRUE)
              ), c(paste0("Mean_", g), paste0("SD_", g)))
            )
          }))
          df_full <- cbind(Protein = rownames(df_de), df_de, group_stats)
        } else {
          df_exprs <- as.data.frame(values$y_protein$E[rownames(df_de), ])
          df_full <- cbind(Protein = rownames(df_de), df_de, df_exprs)
        }

        incProgress(0.3, detail = "Sending data file..."); current_file_uri <- upload_csv_to_gemini(df_full, input$user_api_key)
        qc_final <- NULL; if(!is.null(values$qc_stats) && !is.null(values$metadata)) { qc_final <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>% dplyr::select(Run, Group, Precursors, Proteins, MS1_Signal) }
        # Append phospho context if phospho analysis is active
        auto_msg <- auto_prompt
        if (!is.null(values$phospho_fit) && !is.null(input$phospho_contrast_selector)) {
          phospho_ctx <- tryCatch(
            phospho_ai_context(values$phospho_fit, input$phospho_contrast_selector, values$ksea_results),
            error = function(e) ""
          )
          if (nzchar(phospho_ctx)) auto_msg <- paste0(auto_msg, phospho_ctx)
        }
        incProgress(0.7, detail = "Thinking..."); ai_reply <- ask_gemini_file_chat(auto_msg, current_file_uri, qc_final, input$user_api_key, input$model_name, values$plot_selected_proteins)
      } else { ai_reply <- "Please load data and run analysis first." }
      values$chat_history <- append(values$chat_history, list(list(role = "ai", content = ai_reply)))
    })
  })

  observeEvent(input$send_chat, {
    req(input$chat_input, input$user_api_key)
    values$chat_history <- append(values$chat_history, list(list(role = "user", content = input$chat_input)))
    withProgress(message = "Processing...", {
      if (!is.null(values$fit) && !is.null(values$y_protein)) {
        # Scale protein count to stay within Gemini's token limit (~1M tokens)
        n_samples <- ncol(values$y_protein$E)
        n_max <- if (n_samples > 200) 100 else if (n_samples > 100) 200 else if (n_samples > 50) 400 else 800
        df_de <- topTable(values$fit, coef=input$contrast_selector, number=n_max)
        if (!is.null(values$plot_selected_proteins)) { missing_ids <- setdiff(values$plot_selected_proteins, rownames(df_de)); if (length(missing_ids) > 0) { valid_missing <- intersect(missing_ids, rownames(values$fit$coefficients)); if(length(valid_missing) > 0) { df_extra <- topTable(values$fit, coef=input$contrast_selector, number=Inf)[valid_missing, ]; df_de <- rbind(df_de, df_extra) } } }
        # For large datasets (>100 samples), send group-level summary stats
        if (n_samples > 100 && !is.null(values$metadata)) {
          exprs_mat <- values$y_protein$E[rownames(df_de), , drop = FALSE]
          meta <- values$metadata[values$metadata$Group != "", ]
          group_stats <- do.call(cbind, lapply(unique(meta$Group), function(g) {
            cols <- intersect(meta$File.Name[meta$Group == g], colnames(exprs_mat))
            if (length(cols) == 0) return(NULL)
            data.frame(
              setNames(list(
                rowMeans(exprs_mat[, cols, drop = FALSE], na.rm = TRUE),
                apply(exprs_mat[, cols, drop = FALSE], 1, sd, na.rm = TRUE)
              ), c(paste0("Mean_", g), paste0("SD_", g)))
            )
          }))
          df_full <- cbind(Protein = rownames(df_de), df_de, group_stats)
        } else {
          df_exprs <- as.data.frame(values$y_protein$E[rownames(df_de), ]); df_full <- cbind(Protein = rownames(df_de), df_de, df_exprs)
        }
        incProgress(0.3, detail = "Sending data file..."); current_file_uri <- upload_csv_to_gemini(df_full, input$user_api_key)
        qc_final <- NULL; if(!is.null(values$qc_stats) && !is.null(values$metadata)) { qc_final <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>% dplyr::select(Run, Group, Precursors, Proteins, MS1_Signal) }
        # Append phospho context if phospho analysis is active
        chat_msg <- input$chat_input
        if (!is.null(values$phospho_fit) && !is.null(input$phospho_contrast_selector)) {
          phospho_ctx <- tryCatch(
            phospho_ai_context(values$phospho_fit, input$phospho_contrast_selector, values$ksea_results),
            error = function(e) ""
          )
          if (nzchar(phospho_ctx)) chat_msg <- paste0(chat_msg, phospho_ctx)
        }
        incProgress(0.7, detail = "Thinking..."); ai_reply <- ask_gemini_file_chat(chat_msg, current_file_uri, qc_final, input$user_api_key, input$model_name, values$plot_selected_proteins)
      } else { ai_reply <- "Please load data and run analysis first." }
    })

    ai_selected <- str_extract(ai_reply, "\\[\\[SELECT:.*?\\]\\]")
    if (!is.na(ai_selected)) { raw_ids <- gsub("\\[\\[SELECT:|\\]\\]", "", ai_selected); id_vec <- unlist(strsplit(raw_ids, "[,;]\\s*")); values$plot_selected_proteins <- trimws(id_vec); ai_reply <- gsub("\\[\\[SELECT:.*?\\]\\]", "", ai_reply); ai_reply <- paste0(ai_reply, "\n\n*(I have updated your plots with these highlighted proteins.)*") }
    values$chat_history <- append(values$chat_history, list(list(role = "ai", content = ai_reply))); updateTextAreaInput(session, "chat_input", value = "")
  })

  output$chat_window <- renderUI({ chat_content <- lapply(values$chat_history, function(msg) { if (msg$role == "user") { div(class = "user-msg", span(msg$content)) } else { div(class = "ai-msg", span(markdown(msg$content))) } }); div(class = "chat-container", chat_content) })

  output$download_chat_txt <- downloadHandler(
    filename = function() { req(values$chat_history); paste0("Limpa_Chat_History_", Sys.Date(), ".txt") },
    content = function(file) { req(values$chat_history); text_out <- sapply(values$chat_history, function(msg) { paste0(if(msg$role == "user") "YOU: " else "GEMINI: ", msg$content, "\n---\n") }); writeLines(unlist(text_out), file) }
  )

}
