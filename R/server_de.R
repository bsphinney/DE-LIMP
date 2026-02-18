# ==============================================================================
#  SERVER MODULE — DE Dashboard, Volcano, Consistent DE, Selection Sync
#  Called from app.R as: server_de(input, output, session, values, add_to_log)
# ==============================================================================

server_de <- function(input, output, session, values, add_to_log) {

  # --- volcano_data() reactive (app.R lines 805-830) ---
  volcano_data <- reactive({
    req(values$fit, input$contrast_selector)
    df_raw <- topTable(values$fit, coef=input$contrast_selector, number=Inf) %>% as.data.frame()
    if (!"Protein.Group" %in% colnames(df_raw)) { df <- df_raw %>% rownames_to_column("Protein.Group") } else { df <- df_raw }

    org_db_name <- detect_organism_db(df$Protein.Group)
    df$Accession <- str_split_fixed(df$Protein.Group, "[; ]", 2)[,1]

    id_map <- tryCatch({
      if (!requireNamespace(org_db_name, quietly = TRUE)) { BiocManager::install(org_db_name, ask = FALSE) }
      library(org_db_name, character.only = TRUE)
      bitr(df$Accession, fromType = "UNIPROT", toType = c("SYMBOL", "GENENAME"), OrgDb = get(org_db_name))
    }, error = function(e) { NULL })

    if (!is.null(id_map)) {
      id_map <- id_map %>% distinct(UNIPROT, .keep_all = TRUE)
      df <- df %>% left_join(id_map, by = c("Accession" = "UNIPROT")) %>%
        mutate(Gene = ifelse(is.na(SYMBOL), Accession, SYMBOL), Protein.Name = ifelse(is.na(GENENAME), Protein.Group, GENENAME))
    } else {
      df$Gene <- df$Accession; df$Protein.Name <- df$Protein.Group
    }

    df$Significance <- "Not Sig"; df$Significance[df$adj.P.Val < 0.05 & abs(df$logFC) > input$logfc_cutoff] <- "Significant"
    df$Selected <- "No"; if (!is.null(values$plot_selected_proteins)) { df$Selected[df$Protein.Group %in% values$plot_selected_proteins] <- "Yes" }
    df
  })

  # --- Fullscreen Volcano (app.R lines 2357-2405) ---
  observeEvent(input$fullscreen_volcano, {
    showModal(modalDialog(
      title = "Volcano Plot - Fullscreen View",
      plotlyOutput("volcano_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$volcano_plot_fs <- renderPlotly({
    df <- volcano_data()
    cols <- c("Not Sig" = "grey", "Significant" = "red")

    p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), text = paste("Protein:", Protein.Group), key = Protein.Group, color = Significance)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = cols) +

      # Threshold lines with color
      geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff),
                 linetype = "dashed", color = "#FFA500", size = 0.8) +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed", color = "#4169E1", size = 0.8) +

      theme_minimal() +
      labs(y = "-log10(P-Value)", title = paste0("Volcano Plot: ", input$contrast_selector))

    df_sel <- df %>% filter(Selected == "Yes")
    if (nrow(df_sel) > 0) {
      p <- p + geom_point(data = df_sel, aes(x = logFC, y = -log10(P.Value)),
                         shape = 21, size = 4, fill = NA, color = "blue", stroke = 2)
    }

    # Convert to plotly and add annotations using plotly's native system
    ggplotly(p, tooltip = "text") %>%
      layout(
        annotations = list(
          # Significance criteria box text (using plotly annotations for better positioning)
          list(x = 0.02, y = 0.98, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top",
               text = "<b>Significant if:</b>", showarrow = FALSE, font = list(size = 14)),
          list(x = 0.02, y = 0.93, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top",
               text = paste0("• FDR-adj. p < 0.05<br>• |log2FC| > ", round(input$logfc_cutoff, 2)),
               showarrow = FALSE, font = list(size = 12, color = "#555555"))
        ),
        shapes = list(
          # White background box for legend
          list(type = "rect", x0 = 0.01, x1 = 0.28, y0 = 0.88, y1 = 0.99,
               xref = "paper", yref = "paper", fillcolor = "white", opacity = 0.85,
               line = list(color = "#333333", width = 1))
        )
      )
  })

  # --- Fullscreen Heatmap (app.R lines 2408-2430) ---
  observeEvent(input$fullscreen_heatmap, {
    showModal(modalDialog(
      title = "Heatmap - Fullscreen View",
      plotOutput("heatmap_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$heatmap_plot_fs <- renderPlot({
    req(values$fit, values$y_protein, input$contrast_selector)
    df_volc <- volcano_data(); prot_ids <- NULL
    if (!is.null(input$de_table_rows_selected)) {
      current_table_data <- df_volc
      if (!is.null(values$plot_selected_proteins)) current_table_data <- current_table_data %>% filter(Protein.Group %in% values$plot_selected_proteins)
      prot_ids <- current_table_data$Protein.Group[input$de_table_rows_selected]
    } else if (!is.null(values$plot_selected_proteins)) {
      prot_ids <- values$plot_selected_proteins; if (length(prot_ids) > 50) prot_ids <- head(prot_ids, 50)
    } else { top_prots <- topTable(values$fit, coef = input$contrast_selector, number = 20); prot_ids <- rownames(top_prots) }
    valid_ids <- intersect(prot_ids, rownames(values$y_protein$E)); if (length(valid_ids) == 0) return(NULL)
    mat <- values$y_protein$E[valid_ids, , drop = FALSE]; mat_z <- t(apply(mat, 1, cal_z_score))
    meta <- values$metadata[match(colnames(mat), values$metadata$File.Name), ]; groups <- factor(meta$Group)
    ha <- HeatmapAnnotation(Group = groups, col = list(Group = setNames(rainbow(length(levels(groups))), levels(groups))))
    Heatmap(mat_z, name = "Z-score", top_annotation = ha, cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = FALSE)
  }, height = 700)

  # --- DE Table (app.R lines 2474-2522) ---
  output$de_table <- renderDT({
    req(values$fit, input$contrast_selector)

    # Build table data independently (not using volcano_data() to avoid reactive loops)
    df_raw <- topTable(values$fit, coef=input$contrast_selector, number=Inf) %>% as.data.frame()
    if (!"Protein.Group" %in% colnames(df_raw)) {
      df_full <- df_raw %>% rownames_to_column("Protein.Group")
    } else {
      df_full <- df_raw
    }

    # Add gene symbol and protein name
    org_db_name <- detect_organism_db(df_full$Protein.Group)
    df_full$Accession <- str_split_fixed(df_full$Protein.Group, "[; ]", 2)[,1]

    id_map <- tryCatch({
      if (!requireNamespace(org_db_name, quietly = TRUE)) { BiocManager::install(org_db_name, ask = FALSE) }
      library(org_db_name, character.only = TRUE)
      bitr(df_full$Accession, fromType = "UNIPROT", toType = c("SYMBOL", "GENENAME"), OrgDb = get(org_db_name))
    }, error = function(e) { NULL })

    if (!is.null(id_map)) {
      id_map <- id_map %>% distinct(UNIPROT, .keep_all = TRUE)
      df_full <- df_full %>% left_join(id_map, by = c("Accession" = "UNIPROT")) %>%
        mutate(Gene = ifelse(is.na(SYMBOL), Accession, SYMBOL),
               Protein.Name = ifelse(is.na(GENENAME), Protein.Group, GENENAME))
    } else {
      df_full$Gene <- df_full$Accession
      df_full$Protein.Name <- df_full$Protein.Group
    }

    df_full$Significance <- "Not Sig"
    df_full$Significance[df_full$adj.P.Val < 0.05 & abs(df_full$logFC) > input$logfc_cutoff] <- "Significant"

    # Filter to selected proteins from volcano/AI selection
    if (!is.null(values$plot_selected_proteins) && length(values$plot_selected_proteins) > 0) {
      df_full <- df_full %>% filter(Protein.Group %in% values$plot_selected_proteins)
    }

    df_display <- df_full %>% mutate(across(where(is.numeric), function(x) round(x,4))) %>%
      mutate(Protein.Name_Link = ifelse(!is.na(Accession) & str_detect(Accession, "^[A-Z0-9]{6,}$"),
                                       paste0("<a href='https://www.uniprot.org/uniprotkb/", Accession,
                                              "/entry' target='_blank' onclick='window.open(this.href, \"_blank\"); return false;'>",
                                              Protein.Name, "</a>"),
                                       Protein.Name)) %>%
      dplyr::select(Gene, `Protein Name` = Protein.Name_Link, logFC, P.Value, adj.P.Val, Significance)

    datatable(df_display, selection = "multiple", options = list(pageLength = 10, scrollX = TRUE), escape = FALSE, rownames = FALSE)
  })

  # --- Heatmap Plot (app.R lines 2524-2533) ---
  output$heatmap_plot <- renderPlot({
    req(values$fit, values$y_protein, input$contrast_selector)
    df_volc <- volcano_data(); prot_ids <- NULL
    if (!is.null(input$de_table_rows_selected)) { current_table_data <- df_volc; if (!is.null(values$plot_selected_proteins)) current_table_data <- current_table_data %>% filter(Protein.Group %in% values$plot_selected_proteins); prot_ids <- current_table_data$Protein.Group[input$de_table_rows_selected]
    } else if (!is.null(values$plot_selected_proteins)) { prot_ids <- values$plot_selected_proteins; if(length(prot_ids) > 50) prot_ids <- head(prot_ids, 50)
    } else { top_prots <- topTable(values$fit, coef=input$contrast_selector, number=20); prot_ids <- rownames(top_prots) }
    valid_ids <- intersect(prot_ids, rownames(values$y_protein$E)); if (length(valid_ids) == 0) return(NULL)
    mat <- values$y_protein$E[valid_ids, , drop=FALSE]; mat_z <- t(apply(mat, 1, cal_z_score)); meta <- values$metadata[match(colnames(mat), values$metadata$File.Name), ]; groups <- factor(meta$Group)
    ha <- HeatmapAnnotation(Group = groups, col = list(Group = setNames(rainbow(length(levels(groups))), levels(groups)))); Heatmap(mat_z, name="Z-score", top_annotation = ha, cluster_rows=TRUE, cluster_columns=TRUE, show_column_names=FALSE)
  }, height = 400) # FIXED HEIGHT

  # --- Consistent DE Table (app.R lines 2535-2550) ---
  output$consistent_table <- renderDT({
    req(values$fit, values$y_protein, input$contrast_selector, values$metadata)
    df_res_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>% as.data.frame() %>% filter(adj.P.Val < 0.05)
    if (!"Protein.Group" %in% colnames(df_res_raw)) { df_res <- df_res_raw %>% rownames_to_column("Protein.Group") } else { df_res <- df_res_raw }
    if(nrow(df_res) == 0) return(datatable(data.frame(Status="No significant proteins found.")))
    protein_ids_for_cv <- df_res$Protein.Group; raw_exprs <- values$y_protein$E[protein_ids_for_cv, , drop = FALSE]; linear_exprs <- 2^raw_exprs; cv_list <- list()
    for(g in unique(values$metadata$Group)) {
      if (g == "") next
      files_in_group <- values$metadata$File.Name[values$metadata$Group == g]; group_cols <- intersect(colnames(linear_exprs), files_in_group)
      if (length(group_cols) > 1) { group_data <- linear_exprs[, group_cols, drop = FALSE]; cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
      } else { cv_list[[paste0("CV_", g)]] <- NA }
    }
    cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
    df_final <- left_join(df_res, cv_df, by = "Protein.Group") %>% rowwise() %>% mutate(Avg_CV = mean(c_across(starts_with("CV_")), na.rm = TRUE)) %>% ungroup() %>% arrange(Avg_CV) %>% mutate(Stability = ifelse(Avg_CV < 20, "High", "Low")) %>% dplyr::select(Protein.Group, Stability, Avg_CV, logFC, adj.P.Val, starts_with("CV_")) %>% mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(df_final, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })

  # --- CV Histogram (app.R lines 2553-2630) ---
  output$cv_histogram <- renderPlot({
    req(values$fit, values$y_protein, input$contrast_selector, values$metadata)

    # Get significant proteins
    df_res_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>%
      as.data.frame() %>%
      filter(adj.P.Val < 0.05)
    if (!"Protein.Group" %in% colnames(df_res_raw)) {
      df_res <- df_res_raw %>% rownames_to_column("Protein.Group")
    } else {
      df_res <- df_res_raw
    }
    if(nrow(df_res) == 0) {
      # Show message if no significant proteins
      plot.new()
      text(0.5, 0.5, "No significant proteins found.\nAdjust significance threshold or check data.",
           cex = 1.2, col = "gray50")
      return()
    }

    # Calculate CV for each group
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
        cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) {
          (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
        })
      } else {
        cv_list[[paste0("CV_", g)]] <- NA
      }
    }

    # Convert to long format for ggplot
    cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
    cv_long <- cv_df %>%
      pivot_longer(cols = starts_with("CV_"),
                   names_to = "Group",
                   values_to = "CV") %>%
      mutate(Group = gsub("CV_", "", Group)) %>%
      filter(!is.na(CV))

    # Calculate averages for each group
    cv_averages <- cv_long %>%
      group_by(Group) %>%
      summarise(Avg_CV = mean(CV, na.rm = TRUE), .groups = 'drop')

    # Create histogram with facets
    ggplot(cv_long, aes(x = CV)) +
      geom_histogram(aes(fill = Group), bins = 30, alpha = 0.7, color = "white") +
      geom_vline(data = cv_averages, aes(xintercept = Avg_CV, color = Group),
                 linetype = "dashed", size = 1.2) +
      geom_text(data = cv_averages,
                aes(x = Avg_CV, y = Inf, label = paste0("Avg: ", round(Avg_CV, 1), "%")),
                vjust = 1.5, hjust = -0.1, size = 3.5, fontface = "bold") +
      facet_wrap(~ Group, ncol = 2, scales = "free_y") +
      labs(title = paste0("CV Distribution by Group (", nrow(df_res), " significant proteins)"),
           subtitle = "Dashed line shows average CV for each group",
           x = "Coefficient of Variation (%)",
           y = "Number of Proteins") +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "none",
        strip.background = element_rect(fill = "#667eea", color = NA),
        strip.text = element_text(color = "white", face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray40", size = 11),
        panel.grid.minor = element_blank()
      )
  })

  # --- Fullscreen CV Histogram (app.R lines 2633-2715) ---
  observeEvent(input$fullscreen_cv_hist, {
    req(values$fit, values$y_protein, input$contrast_selector, values$metadata)

    # Get significant proteins
    df_res_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>%
      as.data.frame() %>%
      filter(adj.P.Val < 0.05)
    if (!"Protein.Group" %in% colnames(df_res_raw)) {
      df_res <- df_res_raw %>% rownames_to_column("Protein.Group")
    } else {
      df_res <- df_res_raw
    }
    if(nrow(df_res) == 0) {
      showNotification("No significant proteins found for CV analysis.", type = "warning")
      return()
    }

    # Calculate CV for each group
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
        cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) {
          (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
        })
      } else {
        cv_list[[paste0("CV_", g)]] <- NA
      }
    }

    # Convert to long format for ggplot
    cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
    cv_long <- cv_df %>%
      pivot_longer(cols = starts_with("CV_"),
                   names_to = "Group",
                   values_to = "CV") %>%
      mutate(Group = gsub("CV_", "", Group)) %>%
      filter(!is.na(CV))

    # Calculate averages for each group
    cv_averages <- cv_long %>%
      group_by(Group) %>%
      summarise(Avg_CV = mean(CV, na.rm = TRUE), .groups = 'drop')

    # Create fullscreen histogram
    p <- ggplot(cv_long, aes(x = CV)) +
      geom_histogram(aes(fill = Group), bins = 40, alpha = 0.7, color = "white") +
      geom_vline(data = cv_averages, aes(xintercept = Avg_CV, color = Group),
                 linetype = "dashed", size = 1.5) +
      geom_text(data = cv_averages,
                aes(x = Avg_CV, y = Inf, label = paste0("Avg: ", round(Avg_CV, 1), "%")),
                vjust = 1.5, hjust = -0.1, size = 4, fontface = "bold") +
      facet_wrap(~ Group, ncol = 2, scales = "free_y") +
      labs(title = paste0("CV Distribution by Group (", nrow(df_res), " significant proteins)"),
           subtitle = "Dashed line shows average CV for each group. Lower CV = more stable/reproducible biomarker",
           x = "Coefficient of Variation (%)",
           y = "Number of Proteins") +
      theme_bw(base_size = 16) +
      theme(
        legend.position = "none",
        strip.background = element_rect(fill = "#667eea", color = NA),
        strip.text = element_text(color = "white", face = "bold", size = 14),
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(color = "gray40", size = 12),
        panel.grid.minor = element_blank()
      )

    showModal(modalDialog(
      title = "CV Distribution - Fullscreen View",
      renderPlot({ p }, height = 700, width = 1000),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # --- Download Results CSV (app.R lines 3075-3093) ---
  output$download_result_csv <- downloadHandler(
    filename = function() { paste0("Limpa_Results_", make.names(input$contrast_selector), ".csv") },
    content = function(file) {
      req(values$fit, values$y_protein)
      de_stats <- topTable(values$fit, coef=input$contrast_selector, number=Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(de_stats)) de_stats <- de_stats %>% rownames_to_column("Protein.Group")
      exprs_data <- as.data.frame(values$y_protein$E) %>% rownames_to_column("Protein.Group")
      full_data <- left_join(de_stats, exprs_data, by="Protein.Group")
      write.csv(full_data, file, row.names=FALSE)

      # Log export
      add_to_log("Export Results to CSV", c(
        sprintf("# Exported: %s", basename(file)),
        sprintf("de_stats <- topTable(fit, coef='%s', number=Inf) %%>%% rownames_to_column('Protein.Group')", input$contrast_selector),
        "exprs_data <- as.data.frame(y_protein$E) %>% rownames_to_column('Protein.Group')",
        "full_data <- left_join(de_stats, exprs_data, by='Protein.Group')",
        sprintf("write.csv(full_data, 'Limpa_Results_%s.csv', row.names=FALSE)", make.names(input$contrast_selector))
      ))
    }
  )

  # --- Volcano Plot Interactive (app.R lines 3162-3205) ---
  output$volcano_plot_interactive <- renderPlotly({
    df <- volcano_data()
    cols <- c("Not Sig" = "grey", "Significant" = "red")

    # Use non-adjusted P.Value for y-axis, but color by adjusted p-value significance
    p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), text = paste("Protein:", Protein.Group), key = Protein.Group, color = Significance)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = cols) +

      # Threshold lines with color
      geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff),
                 linetype = "dashed", color = "#FFA500", size = 0.7) +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed", color = "#4169E1", size = 0.7) +

      theme_minimal() +
      labs(y = "-log10(P-Value)", title = paste0("Volcano Plot: ", input$contrast_selector))

    df_sel <- df %>% filter(Selected == "Yes")
    if (nrow(df_sel) > 0) {
      p <- p + geom_point(data = df_sel, aes(x = logFC, y = -log10(P.Value)),
                         shape = 21, size = 4, fill = NA, color = "blue", stroke = 2)
    }

    # Convert to plotly and add annotations using plotly's native system
    ggplotly(p, tooltip = "text", source = "volcano_source") %>%
      layout(
        dragmode = "select",
        annotations = list(
          # Significance criteria box text (using plotly annotations for better positioning)
          list(x = 0.02, y = 0.98, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top",
               text = "<b>Significant if:</b>", showarrow = FALSE, font = list(size = 12)),
          list(x = 0.02, y = 0.93, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top",
               text = paste0("• FDR-adj. p < 0.05<br>• |log2FC| > ", round(input$logfc_cutoff, 2)),
               showarrow = FALSE, font = list(size = 11, color = "#555555"))
        ),
        shapes = list(
          # White background box for legend
          list(type = "rect", x0 = 0.01, x1 = 0.28, y0 = 0.88, y1 = 0.99,
               xref = "paper", yref = "paper", fillcolor = "white", opacity = 0.85,
               line = list(color = "#333333", width = 1))
        )
      )
  })

  # --- Selection Sync: plotly_selected, plotly_click, clear (app.R lines 3207-3209) ---
  observeEvent(event_data("plotly_selected", source = "volcano_source"), { select_data <- event_data("plotly_selected", source = "volcano_source"); if (!is.null(select_data)) values$plot_selected_proteins <- select_data$key })
  observeEvent(event_data("plotly_click", source = "volcano_source"), { click_data <- event_data("plotly_click", source = "volcano_source"); if (!is.null(click_data)) values$plot_selected_proteins <- click_data$key })
  observeEvent(input$clear_plot_selection, { values$plot_selected_proteins <- NULL })

  # --- Table Row Selection Sync (app.R lines 3212-3227) ---
  observeEvent(input$de_table_rows_selected, {
    req(input$de_table_rows_selected, length(input$de_table_rows_selected) > 0)
    df_full <- volcano_data()

    # If table is filtered (volcano/AI selection active), row indices refer to filtered data
    current_selection <- isolate(values$plot_selected_proteins)
    if (!is.null(current_selection) && length(current_selection) > 0) {
      df_full <- df_full %>% filter(Protein.Group %in% current_selection)
    }

    selected_proteins <- df_full$Protein.Group[input$de_table_rows_selected]

    if (length(selected_proteins) > 0) {
      values$plot_selected_proteins <- selected_proteins
    }
  })

  # --- Violin Plot Popup (app.R lines 3230-3277) ---
  observeEvent(input$show_violin, {
    if (is.null(values$plot_selected_proteins) || length(values$plot_selected_proteins) == 0) {
      showNotification("\u26a0\ufe0f Please select a protein in the Volcano Plot or Table first!", type = "warning")
      return()
    }
    # Store all selected proteins for plotting
    values$temp_violin_target <- values$plot_selected_proteins

    n_proteins <- length(values$plot_selected_proteins)
    title_text <- if (n_proteins == 1) {
      paste("Expression Profile:", values$plot_selected_proteins[1])
    } else {
      paste("Expression Profiles for", n_proteins, "Selected Proteins")
    }

    showModal(modalDialog(
      title = title_text,
      size = "xl",
      plotOutput("violin_plot_de_popup", height = paste0(max(400, 200 * ceiling(n_proteins / 2)), "px")),
      footer = modalButton("Close"),
      easyClose = TRUE
    ))
  })

  output$violin_plot_de_popup <- renderPlot({
    req(values$y_protein, values$temp_violin_target, values$metadata)
    prot_ids <- values$temp_violin_target

    # Get expression data for all selected proteins
    exprs_mat <- values$y_protein$E[prot_ids, , drop=FALSE]
    long_df <- as.data.frame(exprs_mat) %>%
      rownames_to_column("Protein") %>%
      pivot_longer(-Protein, names_to = "File.Name", values_to = "LogIntensity")
    long_df <- left_join(long_df, values$metadata, by="File.Name")

    # Create violin plots with faceting for multiple proteins
    ggplot(long_df, aes(x = Group, y = LogIntensity, fill = Group)) +
      geom_violin(alpha = 0.5, trim = FALSE) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      facet_wrap(~Protein, scales = "free_y", ncol = 2) +
      theme_bw() +
      labs(y = "Log2 Intensity", x = "Group") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(face = "bold")
      )
  })

}
