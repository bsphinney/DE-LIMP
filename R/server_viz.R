# ==============================================================================
#  SERVER MODULE — Grid View, Signal Distribution, Dataset Summary
#  Called from app.R as: server_viz(input, output, session, values, add_to_log, is_hf_space)
# ==============================================================================

server_viz <- function(input, output, session, values, add_to_log, is_hf_space) {

  # --- STATUS MESSAGE (line 835) ---
  output$run_status_msg <- renderText({ values$status })

  # --- DATASET SUMMARY (lines 838-956) ---
  # Dataset summary as tab content (instead of modal)
  output$dataset_summary_content <- renderUI({
    req(values$metadata)

    summary_elements <- list()

    # File summary
    summary_elements[[length(summary_elements) + 1]] <- div(
      style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
      tags$h4(icon("file"), " File Summary"),
      tags$hr(),
      tags$p(style = "font-size: 1.1em;",
        icon("folder-open"), " ",
        strong("Total Files: "), nrow(values$metadata)
      ),
      tags$p(style = "font-size: 1.1em;",
        icon("users"), " ",
        strong("Assigned Groups: "),
        length(unique(values$metadata$Group[values$metadata$Group != ""]))
      )
    )

    # Dataset metrics (if pipeline has run)
    if (!is.null(values$y_protein)) {
      avg_signal <- rowMeans(values$y_protein$E, na.rm = TRUE)
      min_linear <- 2^min(avg_signal, na.rm = TRUE)
      max_linear <- 2^max(avg_signal, na.rm = TRUE)

      dynamic_range_text <- if (min_linear > 1e-10) {
        orders_of_magnitude <- log10(max_linear / min_linear)
        paste(round(orders_of_magnitude, 1), "orders of magnitude")
      } else {
        "N/A (Min signal is zero)"
      }

      summary_elements[[length(summary_elements) + 1]] <- div(
        style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
        tags$h4(icon("chart-bar"), " Dataset Metrics"),
        tags$hr(),
        tags$p(style = "font-size: 1.1em;",
          icon("signal"), " ",
          strong("Signal Dynamic Range: "), dynamic_range_text
        ),
        tags$p(style = "font-size: 1.1em;",
          icon("dna"), " ",
          strong("Total Proteins Quantified: "), nrow(values$y_protein$E)
        )
      )
    }

    # Differential expression summary (if DE analysis has run)
    if (!is.null(values$fit)) {
      # Calculate DE proteins per comparison
      all_comparisons <- colnames(values$fit$contrasts)
      de_summary_list <- list()

      for (comp in all_comparisons) {
        de_results <- topTable(values$fit, coef = comp, number = Inf)
        n_sig <- sum(de_results$adj.P.Val < 0.05, na.rm = TRUE)
        n_up <- sum(de_results$adj.P.Val < 0.05 & de_results$logFC > 0, na.rm = TRUE)
        n_down <- sum(de_results$adj.P.Val < 0.05 & de_results$logFC < 0, na.rm = TRUE)

        # Parse comparison name (format: "GroupA - GroupB")
        comp_parts <- strsplit(comp, " - ")[[1]]
        if (length(comp_parts) == 2) {
          group_a <- trimws(comp_parts[1])
          group_b <- trimws(comp_parts[2])
        } else {
          group_a <- "Group 1"
          group_b <- "Group 2"
        }

        # Create explicit, easy-to-read summary
        de_summary_list[[length(de_summary_list) + 1]] <- div(
          style = "margin-bottom: 15px; padding: 12px; background-color: white; border-left: 4px solid #667eea; border-radius: 4px;",
          # Comparison header
          tags$p(style = "font-size: 1.05em; margin-bottom: 8px; font-weight: 500;",
            icon("microscope"), " ",
            strong(comp),
            span(style = "margin-left: 10px; color: #6c757d; font-weight: normal; font-size: 0.9em;",
              paste0("(", n_sig, " significant proteins)")
            )
          ),
          # Detailed breakdown
          if (n_sig > 0) {
            tagList(
              tags$div(style = "margin-left: 25px; font-size: 0.95em;",
                tags$div(style = "margin-bottom: 4px;",
                  span(style = "color: #e41a1c; font-weight: 500;", "\u2191 ", n_up),
                  " proteins higher in ",
                  strong(style = "color: #e41a1c;", group_a)
                ),
                tags$div(
                  span(style = "color: #377eb8; font-weight: 500;", "\u2193 ", n_down),
                  " proteins higher in ",
                  strong(style = "color: #377eb8;", group_b)
                )
              )
            )
          } else {
            tags$div(style = "margin-left: 25px; font-size: 0.9em; color: #6c757d; font-style: italic;",
              "No significant differences detected"
            )
          }
        )
      }

      summary_elements[[length(summary_elements) + 1]] <- div(
        style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
        tags$h4(icon("flask"), " Differential Expression Summary"),
        tags$hr(),
        tags$p(style = "font-size: 0.85em; color: #6c757d; margin-bottom: 15px;",
          "Proteins with FDR-adjusted p-value < 0.05. Arrows indicate direction of change."
        ),
        do.call(tagList, de_summary_list)
      )

      # Complete dataset export button
      summary_elements[[length(summary_elements) + 1]] <- div(
        style = "background-color: #eef2ff; padding: 20px; border-radius: 8px; border: 1px solid #c7d2fe;",
        div(style = "display: flex; justify-content: space-between; align-items: center;",
          div(
            tags$h4(icon("download"), " Export Complete Dataset", style = "margin: 0 0 5px 0;"),
            tags$p(style = "font-size: 0.85em; color: #6c757d; margin: 0;",
              "DE statistics for all comparisons, expression values, gene symbols, and sample metadata — ready for downstream tools."
            )
          ),
          downloadButton("download_complete_dataset", "Download Complete Dataset",
            class = "btn-primary", style = "white-space: nowrap;")
        )
      )
    }

    tagList(summary_elements)
  })

  # --- COMPLETE DATASET EXPORT ---
  output$download_complete_dataset <- downloadHandler(
    filename = function() {
      paste0("Limpa_Complete_Dataset_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(values$fit, values$y_protein, values$metadata)

      withProgress(message = "Building complete dataset...", value = 0, {

        # --- 1. Gene symbol mapping ---
        incProgress(0.1, detail = "Mapping gene symbols...")
        protein_ids <- rownames(values$y_protein$E)
        accessions <- str_split_fixed(protein_ids, "[; ]", 2)[,1]
        org_db_name <- detect_organism_db(protein_ids)

        id_map <- tryCatch({
          if (!requireNamespace(org_db_name, quietly = TRUE)) BiocManager::install(org_db_name, ask = FALSE)
          library(org_db_name, character.only = TRUE)
          db_obj <- get(org_db_name)
          AnnotationDbi::select(db_obj, keys = accessions, columns = c("SYMBOL"), keytype = "UNIPROT") %>%
            dplyr::rename(Accession = UNIPROT, Gene = SYMBOL) %>% distinct(Accession, .keep_all = TRUE)
        }, error = function(e) data.frame(Accession = accessions, Gene = accessions))

        gene_df <- data.frame(Protein.Group = protein_ids, Accession = accessions, stringsAsFactors = FALSE)
        gene_df <- left_join(gene_df, id_map, by = "Accession")
        gene_df$Gene[is.na(gene_df$Gene)] <- gene_df$Accession[is.na(gene_df$Gene)]

        # --- 2. DE stats for ALL contrasts ---
        incProgress(0.3, detail = "Gathering DE statistics for all comparisons...")
        all_contrasts <- colnames(values$fit$contrasts)
        de_combined <- data.frame(Protein.Group = protein_ids, stringsAsFactors = FALSE)

        for (cname in all_contrasts) {
          tt <- topTable(values$fit, coef = cname, number = Inf) %>% as.data.frame()
          if (!"Protein.Group" %in% colnames(tt)) tt <- tt %>% rownames_to_column("Protein.Group")
          safe_name <- make.names(cname)
          tt_subset <- tt %>% dplyr::select(Protein.Group, logFC, P.Value, adj.P.Val)
          colnames(tt_subset) <- c("Protein.Group",
            paste0("logFC_", safe_name),
            paste0("P.Value_", safe_name),
            paste0("adj.P.Val_", safe_name))
          de_combined <- left_join(de_combined, tt_subset, by = "Protein.Group")
        }

        # --- 3. Expression matrix ---
        incProgress(0.6, detail = "Adding expression values...")
        exprs_df <- as.data.frame(values$y_protein$E) %>% rownames_to_column("Protein.Group")

        # --- 4. Combine everything ---
        incProgress(0.8, detail = "Assembling and writing file...")
        full_export <- gene_df %>%
          dplyr::select(Protein.Group, Accession, Gene) %>%
          left_join(de_combined, by = "Protein.Group") %>%
          left_join(exprs_df, by = "Protein.Group")

        # --- 5. Add sample metadata as header rows ---
        # Build metadata annotation lines that map to sample columns
        sample_cols <- colnames(values$y_protein$E)
        meta <- values$metadata
        non_sample_cols <- setdiff(colnames(full_export), sample_cols)

        # Get custom covariate names
        cov1_name <- if (!is.null(values$cov1_name) && nzchar(values$cov1_name)) values$cov1_name else "Covariate1"
        cov2_name <- if (!is.null(values$cov2_name) && nzchar(values$cov2_name)) values$cov2_name else "Covariate2"

        # Build annotation rows
        annot_rows <- list()
        annot_fields <- list(
          Group = "Group", Batch = "Batch",
          Covariate1 = cov1_name, Covariate2 = cov2_name
        )
        for (col_name in names(annot_fields)) {
          if (col_name %in% colnames(meta) && any(nzchar(meta[[col_name]]))) {
            label <- annot_fields[[col_name]]
            vals <- meta[[col_name]][match(sample_cols, meta$File.Name)]
            vals[is.na(vals)] <- ""
            row <- c(paste0("#", label), rep("", length(non_sample_cols) - 1), vals)
            annot_rows[[length(annot_rows) + 1]] <- row
          }
        }

        # Write: annotation header rows, then column names, then data
        con <- file(file, "w")
        for (arow in annot_rows) {
          writeLines(paste(arow, collapse = ","), con)
        }
        close(con)
        write.table(full_export, file, sep = ",", row.names = FALSE, quote = TRUE, append = TRUE)

        incProgress(1.0, detail = "Done!")
      })
    }
  )

  # --- GRID VIEW & PLOT LOGIC (lines 1064-1191) ---
  grid_react_df <- reactive({
    req(values$fit, values$y_protein, values$metadata, input$contrast_selector)

    # Build volcano-style DE data directly (volcano_data() is local to server_de)
    df_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>% as.data.frame()
    if (!"Protein.Group" %in% colnames(df_raw)) {
      df_raw <- df_raw %>% rownames_to_column("Protein.Group")
    }
    df_raw$Accession <- str_split_fixed(df_raw$Protein.Group, "[; ]", 2)[, 1]

    org_db_name <- detect_organism_db(df_raw$Protein.Group)
    id_map <- tryCatch({
      if (!requireNamespace(org_db_name, quietly = TRUE)) BiocManager::install(org_db_name, ask = FALSE)
      library(org_db_name, character.only = TRUE)
      bitr(df_raw$Accession, fromType = "UNIPROT", toType = c("SYMBOL", "GENENAME"), OrgDb = get(org_db_name))
    }, error = function(e) NULL)

    if (!is.null(id_map)) {
      id_map <- id_map %>% distinct(UNIPROT, .keep_all = TRUE)
      df_raw <- df_raw %>%
        left_join(id_map, by = c("Accession" = "UNIPROT")) %>%
        mutate(Gene = ifelse(is.na(SYMBOL), Accession, SYMBOL),
               Protein.Name = ifelse(is.na(GENENAME), Protein.Group, GENENAME))
    } else {
      df_raw$Gene <- df_raw$Accession
      df_raw$Protein.Name <- df_raw$Protein.Group
    }

    df_raw$Significance <- "Not Sig"
    df_raw$Significance[df_raw$adj.P.Val < 0.05 & abs(df_raw$logFC) > input$logfc_cutoff] <- "Significant"

    df_volc <- df_raw
    df_exprs <- as.data.frame(values$y_protein$E) %>% rownames_to_column("Protein.Group")
    df_merged <- left_join(df_volc, df_exprs, by="Protein.Group")
    if (!is.null(values$plot_selected_proteins)) { df_merged <- df_merged %>% filter(Protein.Group %in% values$plot_selected_proteins) }

    meta_sorted <- values$metadata %>% arrange(Group, File.Name)
    ordered_files <- meta_sorted$File.Name
    valid_cols <- intersect(ordered_files, colnames(df_exprs))
    run_ids <- values$metadata$ID[match(valid_cols, values$metadata$File.Name)]
    new_headers <- as.character(run_ids)

    df_final <- df_merged %>%
      dplyr::select(Protein.Group, Gene, Protein.Name, Significance, logFC, P.Value, adj.P.Val, all_of(valid_cols)) %>%
      mutate(across(where(is.numeric), ~round(., 2)))

    df_final$Original.ID <- df_final$Protein.Group
    fixed_cols <- c("Protein.Group", "Gene", "Protein.Name", "Significance", "logFC", "P.Value", "adj.P.Val")
    colnames(df_final) <- c(fixed_cols, new_headers, "Original.ID")

    unique_groups <- sort(unique(meta_sorted$Group))
    group_colors <- setNames(rainbow(length(unique_groups), v=0.85, s=0.8), unique_groups)

    list(data = df_final, fixed_cols = fixed_cols, expr_cols = new_headers, valid_cols_map = valid_cols, meta_sorted = meta_sorted, group_colors = group_colors)
  })

  # Grid View legend UI
  output$grid_legend_ui <- renderUI({
    gdata <- grid_react_df()
    tags$div(
      style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
      tags$strong(icon("palette"), " Condition Legend:"), tags$br(),
      lapply(names(gdata$group_colors), function(grp) {
        tags$span(
          style = paste0("background-color:", gdata$group_colors[[grp]], "; color:white; padding:4px 10px; margin-right:8px; border-radius:4px; display:inline-block; margin-top:5px;"),
          grp
        )
      })
    )
  })

  # Grid View file mapping UI
  output$grid_file_map_ui <- renderUI({
    gdata <- grid_react_df()
    tags$details(
      style = "margin-bottom: 10px;",
      tags$summary(
        style = "cursor: pointer; color: #0d6efd;",
        icon("list"), " Click to view File ID Mapping (Run # \u2192 Filename)"
      ),
      tags$div(
        style = "max-height: 200px; overflow-y: auto; background: #f9f9f9; padding: 10px; border: 1px solid #dee2e6; border-radius: 4px; margin-top: 8px;",
        lapply(1:nrow(gdata$meta_sorted), function(i) {
          row <- gdata$meta_sorted[i, ]
          tags$div(
            style = "padding: 2px 0;",
            tags$span(style = "font-weight:bold; color:#007bff;", paste0("[", row$ID, "] ")),
            tags$span(row$File.Name),
            tags$span(style = "color:#6c757d; font-size:0.9em;", paste0(" (", row$Group, ")"))
          )
        })
      )
    )
  })

  observeEvent(input$grid_reset_selection, { values$plot_selected_proteins <- NULL })

  output$grid_view_table <- renderDT({
    gdata <- grid_react_df(); df_display <- gdata$data
    acc <- str_split_fixed(df_display$Original.ID, "[; ]", 2)[,1]
    df_display$Protein.Group <- paste0("<a href='https://www.uniprot.org/uniprotkb/", acc, "/entry' target='_blank'>", df_display$Protein.Group, "</a>")
    df_display <- df_display %>% dplyr::select(-Original.ID)
    fixed_cols <- gdata$fixed_cols; expr_cols <- gdata$expr_cols
    valid_cols_map <- gdata$valid_cols_map; meta_sorted <- gdata$meta_sorted; group_colors <- gdata$group_colors

    header_html <- tags$table(class = "display", tags$thead(tags$tr(lapply(seq_along(colnames(df_display)), function(i) {
      col_name <- colnames(df_display)[i]
      if (i > length(fixed_cols)) {
        original_name <- valid_cols_map[i - length(fixed_cols)]; grp <- meta_sorted$Group[meta_sorted$File.Name == original_name]; bg_color <- group_colors[grp]
        tags$th(title = paste("File:", original_name, "\nGroup:", grp), col_name, style = paste0("background-color: ", bg_color, "; color: white; text-align: center;"))
      } else { tags$th(col_name) }
    }))))

    expression_matrix <- as.matrix(df_display[, expr_cols]); brks <- quantile(expression_matrix, probs = seq(.05, .95, .05), na.rm = TRUE); clrs <- colorRampPalette(c("#4575b4", "white", "#d73027"))(length(brks) + 1)
    datatable(df_display, container = header_html, selection = 'single', escape = FALSE, options = list(dom = 'frtip', pageLength = 15, scrollX = TRUE, columnDefs = list(list(className = 'dt-center', targets = (length(fixed_cols)):(ncol(df_display)-1)))), rownames = FALSE) %>% formatStyle(expr_cols, backgroundColor = styleInterval(brks, clrs))
  })

  output$download_grid_data <- downloadHandler(
    filename = function() { paste0("DE_LIMP_Grid_Export_", Sys.Date(), ".csv") },
    content = function(file) {
      gdata <- grid_react_df(); df_export <- gdata$data %>% dplyr::select(-Original.ID)
      fixed_cols <- gdata$fixed_cols; real_names <- gdata$valid_cols_map
      colnames(df_export) <- c(fixed_cols, real_names)
      write.csv(df_export, file, row.names = FALSE)

      # Log export
      add_to_log("Export Grid View to CSV", c(
        sprintf("# Exported full expression matrix with DE stats"),
        sprintf("# File: DE_LIMP_Grid_Export_%s.csv", Sys.Date()),
        "# (Combined DE results + expression values for all samples)"
      ))
    }
  )

  observeEvent(input$grid_view_table_rows_selected, {
    req(grid_react_df()); selected_idx <- input$grid_view_table_rows_selected
    if (length(selected_idx) > 0) {
      gdata <- grid_react_df(); selected_id <- gdata$data$Original.ID[selected_idx]; values$grid_selected_protein <- selected_id
      xic_btn <- if (!is_hf_space) actionButton("show_xic_from_grid", "\U0001F4C8 XICs", class="btn-info") else NULL
      showModal(modalDialog(title = paste("Expression Plot:", selected_id), size = "xl", plotOutput("violin_plot_grid", height = "600px"), footer = tagList(xic_btn, actionButton("back_to_grid", "Back to Grid", class="btn-info"), modalButton("Close")), easyClose = TRUE))
    }
  })

  observeEvent(input$back_to_grid, { click("show_grid_view") })

  output$violin_plot_grid <- renderPlot({
    req(values$y_protein, values$grid_selected_protein, values$metadata)
    prot_id <- values$grid_selected_protein
    exprs_mat <- values$y_protein$E[prot_id, , drop=FALSE]
    long_df <- as.data.frame(exprs_mat) %>% rownames_to_column("Protein") %>% pivot_longer(-Protein, names_to = "File.Name", values_to = "LogIntensity")
    long_df <- left_join(long_df, values$metadata, by="File.Name")

    ggplot(long_df, aes(x = Group, y = LogIntensity, fill = Group)) +
      geom_violin(alpha = 0.5, trim = FALSE) + geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      facet_wrap(~Protein, scales = "free_y") + theme_bw() +
      labs(title = paste("Protein:", prot_id), y = "Log2 Intensity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, height = 600) # FIXED HEIGHT

  # --- SIGNAL DISTRIBUTION (lines 1193-1242) ---
  output$protein_signal_plot <- renderPlot({
    req(values$y_protein)
    avg_signal <- rowMeans(values$y_protein$E, na.rm = TRUE)
    plot_df <- data.frame(Protein.Group = names(avg_signal), Average_Signal_Log2 = avg_signal) %>%
      mutate(Average_Signal_Log10 = Average_Signal_Log2 / log2(10))

    # Always show DE coloring when results are available
    if (!is.null(values$fit) && !is.null(input$contrast_selector_signal) && nchar(input$contrast_selector_signal) > 0) {
      de_data_raw <- topTable(values$fit, coef = input$contrast_selector_signal, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(de_data_raw)) {
        de_data_intermediate <- de_data_raw %>% rownames_to_column("Protein.Group")
      } else {
        de_data_intermediate <- de_data_raw
      }
      de_data <- de_data_intermediate %>%
        mutate(DE_Status = case_when(
          adj.P.Val < 0.05 & logFC > input$logfc_cutoff ~ "Up-regulated",
          adj.P.Val < 0.05 & logFC < -input$logfc_cutoff ~ "Down-regulated",
          TRUE ~ "Not Significant"
        )) %>%
        dplyr::select(Protein.Group, DE_Status)
      plot_df <- left_join(plot_df, de_data, by = "Protein.Group")
      plot_df$DE_Status[is.na(plot_df$DE_Status)] <- "Not Significant"
    } else {
      plot_df$DE_Status <- "Not Significant"
    }

    if (!is.null(values$plot_selected_proteins)) {
      plot_df$Is_Selected <- plot_df$Protein.Group %in% values$plot_selected_proteins
    } else {
      plot_df$Is_Selected <- FALSE
    }
    selected_df <- filter(plot_df, Is_Selected)

    # Build plot - always use DE coloring when available
    p <- ggplot(plot_df, aes(x = reorder(Protein.Group, -Average_Signal_Log10), y = Average_Signal_Log10))
    if (!is.null(values$fit)) {
      p <- p + geom_point(aes(color = DE_Status), size = 1.5) +
        scale_color_manual(name = "DE Status",
          values = c("Up-regulated" = "#e41a1c", "Down-regulated" = "#377eb8", "Not Significant" = "grey70"))
    } else {
      p <- p + geom_point(color = "cornflowerblue", size = 1.5)
    }

    p + labs(title = "Signal Distribution Across All Protein Groups", x = NULL, y = "Average Signal (Log10 Intensity)") +
      theme_minimal() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_discrete(expand = expansion(add = 1)) +
      geom_point(data = selected_df, color = "black", shape = 1, size = 4, stroke = 1) +
      geom_text_repel(data = selected_df, aes(label = Protein.Group), size = 4, max.overlaps = 20)
  }, height = 500) # FIXED HEIGHT

  # --- GROUP QC SUMMARY TABLE (lines 1244-1248) ---
  output$group_summary_table <- renderDT({
    req(values$qc_stats, values$metadata)
    summary_df <- values$qc_stats %>% left_join(values$metadata, by = c("Run" = "File.Name")) %>% filter(Group != "") %>% group_by(Group) %>% summarise(`Avg. Precursors` = mean(Precursors, na.rm = TRUE), `Avg. MS1 Signal` = mean(MS1_Signal, na.rm = TRUE), `Avg. Proteins` = mean(Proteins, na.rm = TRUE)) %>% mutate(across(where(is.numeric), ~round(.x, 0)))
    datatable(summary_df, options = list(dom = 't', pageLength = 10), rownames = FALSE)
  })

  # --- FULLSCREEN: Signal Distribution (lines 1734-1790) ---
  observeEvent(input$fullscreen_signal, {
    showModal(modalDialog(
      title = "Signal Distribution - Fullscreen View",
      plotOutput("protein_signal_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$protein_signal_plot_fs <- renderPlot({
    req(values$y_protein)
    avg_signal <- rowMeans(values$y_protein$E, na.rm = TRUE)
    plot_df <- data.frame(Protein.Group = names(avg_signal), Average_Signal_Log2 = avg_signal) %>%
      mutate(Average_Signal_Log10 = Average_Signal_Log2 / log2(10))

    # Always show DE coloring when results are available
    if (!is.null(values$fit) && !is.null(input$contrast_selector_signal) && nchar(input$contrast_selector_signal) > 0) {
      de_data_raw <- topTable(values$fit, coef = input$contrast_selector_signal, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(de_data_raw)) {
        de_data_intermediate <- de_data_raw %>% rownames_to_column("Protein.Group")
      } else {
        de_data_intermediate <- de_data_raw
      }
      de_data <- de_data_intermediate %>%
        mutate(DE_Status = case_when(
          adj.P.Val < 0.05 & logFC > input$logfc_cutoff ~ "Up-regulated",
          adj.P.Val < 0.05 & logFC < -input$logfc_cutoff ~ "Down-regulated",
          TRUE ~ "Not Significant"
        )) %>%
        dplyr::select(Protein.Group, DE_Status)
      plot_df <- left_join(plot_df, de_data, by = "Protein.Group")
      plot_df$DE_Status[is.na(plot_df$DE_Status)] <- "Not Significant"
    } else {
      plot_df$DE_Status <- "Not Significant"
    }

    if (!is.null(values$plot_selected_proteins)) {
      plot_df$Is_Selected <- plot_df$Protein.Group %in% values$plot_selected_proteins
    } else {
      plot_df$Is_Selected <- FALSE
    }
    selected_df <- filter(plot_df, Is_Selected)

    # Build plot - always use DE coloring when available
    p <- ggplot(plot_df, aes(x = reorder(Protein.Group, -Average_Signal_Log10), y = Average_Signal_Log10))
    if (!is.null(values$fit)) {
      p <- p + geom_point(aes(color = DE_Status), size = 1.5) +
        scale_color_manual(name = "DE Status",
          values = c("Up-regulated" = "#e41a1c", "Down-regulated" = "#377eb8", "Not Significant" = "grey70"))
    } else {
      p <- p + geom_point(color = "cornflowerblue", size = 1.5)
    }

    p + labs(title = "Signal Distribution Across All Protein Groups", x = NULL, y = "Average Signal (Log10 Intensity)") +
      theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_discrete(expand = expansion(add = 1)) +
      geom_point(data = selected_df, color = "black", shape = 1, size = 4, stroke = 1) +
      geom_text_repel(data = selected_df, aes(label = Protein.Group), size = 4, max.overlaps = 20)
  }, height = 700)

}
