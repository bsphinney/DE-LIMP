server_ai <- function(input, output, session, values) {

  # --- AI SUMMARY (Data Overview Tab) ---
  observeEvent(input$generate_ai_summary_overview, {
    req(values$fit, input$contrast_selector, values$y_protein, input$user_api_key)

    withProgress(message = "Generating AI Summary...", value = 0, {
      incProgress(0.2, detail = "Gathering DE data...")

      # Get volcano data for the current contrast
      de_results_full <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(de_results_full)) {
        de_results_full <- de_results_full %>% rownames_to_column("Protein.Group")
      }

      # Add gene symbols
      de_results_full$Accession <- str_split_fixed(de_results_full$Protein.Group, "[; ]", 2)[,1]
      org_db_name <- detect_organism_db(de_results_full$Protein.Group)

      id_map <- tryCatch({
        if (!requireNamespace(org_db_name, quietly = TRUE)) { BiocManager::install(org_db_name, ask = FALSE) }
        library(org_db_name, character.only = TRUE)
        db_obj <- get(org_db_name)
        AnnotationDbi::select(db_obj, keys = de_results_full$Accession, columns = c("SYMBOL"), keytype = "UNIPROT") %>%
          rename(Accession = UNIPROT, Gene = SYMBOL) %>% distinct(Accession, .keep_all = TRUE)
      }, error = function(e) data.frame(Accession = de_results_full$Accession, Gene = de_results_full$Accession))

      de_results_with_genes <- left_join(de_results_full, id_map, by = "Accession")
      de_results_with_genes$Gene[is.na(de_results_with_genes$Gene)] <- de_results_with_genes$Accession[is.na(de_results_with_genes$Gene)]

      # Filter for significant and get top 50
      top_de_data <- de_results_with_genes %>%
        filter(adj.P.Val < 0.05) %>%
        arrange(adj.P.Val) %>%
        head(50) %>%
        dplyr::select(Gene, logFC, adj.P.Val) %>%
        mutate(across(where(is.numeric), ~round(.x, 3)))

      top_de_text <- paste(capture.output(print(as.data.frame(top_de_data))), collapse = "\n")

      incProgress(0.5, detail = "Gathering stable proteins...")

      stable_prots_df <- tryCatch({
        df_res <- de_results_with_genes %>% filter(adj.P.Val < 0.05)
        if(nrow(df_res) == 0) return(data.frame(Info="No significant proteins to assess for stability."))

        protein_ids_for_cv <- df_res$Protein.Group
        raw_exprs <- values$y_protein$E[protein_ids_for_cv, , drop = FALSE]
        linear_exprs <- 2^raw_exprs
        cv_list <- list()

        for(g in unique(values$metadata$Group)) {
          if (g == "") next
          files_in_group <- values$metadata$File.Name[values$metadata$Group == g]
          group_cols <- intersect(colnames(linear_exprs), files_in_group)
          if (length(group_cols) > 1) {
            group_data <- linear_exprs[, group_cols, drop = FALSE]
            cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
          } else {
            cv_list[[paste0("CV_", g)]] <- NA
          }
        }

        cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
        left_join(df_res, cv_df, by = "Protein.Group") %>%
          rowwise() %>%
          mutate(Avg_CV = mean(c_across(starts_with("CV_")), na.rm = TRUE)) %>%
          ungroup() %>%
          arrange(Avg_CV) %>%
          head(3) %>%
          dplyr::select(Gene, Avg_CV, logFC, adj.P.Val) %>%
          mutate(across(where(is.numeric), ~round(.x, 2)))
      }, error = function(e) data.frame(Error = "Could not calculate stable proteins."))

      stable_prots_text <- paste(capture.output(print(as.data.frame(stable_prots_df))), collapse = "\n")

      incProgress(0.7, detail = "Constructing prompt...")

      system_prompt <- paste0(
        "You are a senior proteomics consultant at a core facility. Write a 3-paragraph summary of the differential expression results below.\n\n",
        "Paragraph 1: Overview of the comparison being analyzed and total number of significant proteins (FDR < 0.05).\n",
        "Paragraph 2: Highlight the top upregulated and downregulated proteins by fold change, using gene names when available.\n",
        "Paragraph 3: Mention the most stable differentially expressed proteins (lowest CV) as potential high-confidence biomarkers.\n\n",
        "Use markdown formatting. Be concise and scientific."
      )

      final_prompt <- paste0(
        system_prompt,
        "\n\n--- DATA FOR SUMMARY ---\n\n",
        "Comparison: ", input$contrast_selector, "\n\n",
        "Top Significant Proteins:\n", top_de_text, "\n\n",
        "Most Stable DE Proteins:\n", stable_prots_text
      )

      incProgress(0.8, detail = "Asking AI...")
      ai_summary <- ask_gemini_text_chat(final_prompt, input$user_api_key, input$model_name)

      # Render the summary to the output area
      output$ai_summary_output <- renderUI({
        div(style = "background-color: #ffffff; padding: 20px; border: 1px solid #dee2e6; border-radius: 8px;",
          tags$h5(icon("check-circle"), " Analysis Complete", style = "color: #28a745; margin-bottom: 15px;"),
          HTML(markdown::markdownToHTML(text = ai_summary, fragment.only = TRUE))
        )
      })
    })
  })

  # --- Data Chat Info Modal ---
  observeEvent(input$data_chat_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About Data Chat"),
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
