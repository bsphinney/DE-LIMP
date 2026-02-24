# ==============================================================================
#  SERVER MODULE — XIC Viewer (chromatogram inspection)
#  Called from app.R as: server_xic(input, output, session, values, is_hf_space)
# ==============================================================================

server_xic <- function(input, output, session, values, is_hf_space) {

  # ============================================================================
  #      XIC Viewer: Server Logic
  # ============================================================================

  # --- XIC Directory Loading ---
  observeEvent(input$xic_load_dir, {
    req(input$xic_dir_input)
    xic_path <- trimws(input$xic_dir_input)

    # Smart path resolution: if user enters a .parquet file path, derive the _xic sibling
    if (grepl("\\.parquet$", xic_path, ignore.case = TRUE) && file.exists(xic_path)) {
      report_stem <- tools::file_path_sans_ext(basename(xic_path))
      derived_dir <- file.path(dirname(xic_path), paste0(report_stem, "_xic"))
      if (dir.exists(derived_dir)) {
        xic_path <- derived_dir
        updateTextInput(session, "xic_dir_input", value = xic_path)
        showNotification(paste0("Derived XIC directory: ", basename(xic_path)),
          type = "message", duration = 3)
      }
    }

    # If path doesn't end with _xic, try appending _xic (user may have pasted the parent dir)
    if (!dir.exists(xic_path) && !grepl("_xic$", basename(xic_path))) {
      candidate <- paste0(xic_path, "_xic")
      if (dir.exists(candidate)) {
        xic_path <- candidate
        updateTextInput(session, "xic_dir_input", value = xic_path)
      }
    }

    if (!dir.exists(xic_path)) {
      showNotification("Directory not found. Check the path and try again.",
        type = "error", duration = 5)
      values$xic_available <- FALSE
      return()
    }

    xic_files <- list.files(xic_path, pattern = "\\.xic\\.parquet$",
                            full.names = TRUE, recursive = TRUE)

    if (length(xic_files) > 0) {
      values$xic_dir <- xic_path
      values$xic_available <- TRUE

      # Detect XIC format version (v1 = DIA-NN 1.x wide, v2 = DIA-NN 2.x long)
      values$xic_format <- detect_xic_format(xic_path)
      message(paste("Detected XIC format:", values$xic_format))

      # Detect mobilogram files and check if they contain non-zero data
      # DIA-NN names these: .ms1_mobilogram.parquet, .ms2_mobilogram.parquet
      mob_files <- list.files(xic_path, pattern = "_mobilogram\\.parquet$",
                              full.names = TRUE, recursive = TRUE)
      values$mobilogram_files_found <- length(mob_files)
      if (length(mob_files) > 0) {
        # Sample one file to check for non-zero data (IM instruments only)
        mob_has_data <- tryCatch({
          mob_sample <- arrow::read_parquet(mob_files[1], as_data_frame = TRUE)
          any(mob_sample$data != 0, na.rm = TRUE)
        }, error = function(e) FALSE)
        values$mobilogram_available <- mob_has_data
        if (mob_has_data) values$mobilogram_dir <- xic_path
        message(paste("Mobilogram files:", length(mob_files),
                      "| Has IM data:", mob_has_data))
      } else {
        values$mobilogram_available <- FALSE
      }

      # Build precursor map directly from loaded data (no file re-reading needed)
      if (!is.null(values$raw_data)) {
        tryCatch({
          precursor_ids <- rownames(values$raw_data$E)
          protein_groups <- values$raw_data$genes$Protein.Group
          values$xic_report_map <- data.frame(
            Precursor.Id = precursor_ids,
            Protein.Group = protein_groups,
            stringsAsFactors = FALSE
          )
          message(paste("XIC precursor map built from loaded data:",
                        length(unique(protein_groups)), "proteins,",
                        length(unique(precursor_ids)), "precursors"))
        }, error = function(e) {
          values$xic_report_map <- NULL
          showNotification(paste("Error building precursor map:", e$message),
            type = "error", duration = 8)
        })
      } else {
        values$xic_report_map <- NULL
        showNotification(
          "Please load your DIA-NN report data first (Step 1), then load XICs.",
          type = "warning", duration = 8)
      }

      status_msg <- paste("Found", length(xic_files), "XIC files")
      if (values$mobilogram_available) {
        status_msg <- paste0(status_msg, " + ion mobility data")
      } else if (values$mobilogram_files_found > 0) {
        status_msg <- paste0(status_msg, " (mobilogram files present but contain no IM data)")
      }
      showNotification(paste0(status_msg, ". Select a protein to view chromatograms."),
        type = "message", duration = 5)
    } else {
      values$xic_available <- FALSE
      showNotification(
        "No .xic.parquet files found in selected directory.",
        type = "warning", duration = 5)
    }
  })

  # --- XIC Status Badge ---
  output$xic_status_badge <- renderUI({
    if (values$xic_available) {
      xic_files <- list.files(values$xic_dir, pattern = "\\.xic\\.parquet$")
      badge_text <- paste(length(xic_files), "XIC files ready")
      if (values$mobilogram_available) badge_text <- paste0(badge_text, " + IM")
      format_label <- if (values$xic_format == "v2") "DIA-NN 2.x" else "DIA-NN 1.x"
      div(class = "alert alert-success py-1 px-2 mt-2 mb-0",
        style = "font-size: 0.85em;",
        icon("check-circle"),
        badge_text, br(),
        span(class = "text-muted", style = "font-size: 0.85em;", format_label))
    } else {
      NULL
    }
  })

  # --- XIC Modal Function ---
  show_xic_modal <- function(session, values) {
    showModal(modalDialog(
      title = div(style = "display: flex; align-items: center; gap: 10px;",
        icon("wave-square"),
        span("XIC Chromatograms:"),
        span(values$xic_protein, style = "color: #667eea; font-weight: bold;"),
        span(textOutput("xic_precursor_count", inline = TRUE),
          style = "font-size: 0.8em; color: #6c757d; margin-left: 10px;")
      ),

      # Controls row
      div(style = "display: flex; flex-wrap: wrap; gap: 12px; align-items: flex-end; margin-bottom: 12px;",
        selectInput("xic_display_mode", "Display:",
          choices = c("Facet by sample" = "overlay",
                      "Facet by fragment" = "facet",
                      "Intensity alignment" = "alignment"),
          selected = "overlay", width = "220px"),

        selectInput("xic_precursor_select", "Precursor:",
          choices = NULL, width = "280px"),

        selectInput("xic_group_filter", "Filter Group:",
          choices = NULL, width = "180px"),

        # MS1 and IM controls hidden in alignment mode (irrelevant for bar charts)
        conditionalPanel(
          condition = "input.xic_display_mode != 'alignment'",
          div(style = "display: flex; flex-wrap: wrap; gap: 12px; align-items: center;",
            checkboxInput("xic_show_ms1", "Show MS1 (split axis)", value = FALSE),
            # Ion mobility toggle — only shown when mobilogram data is available
            if (values$mobilogram_available) {
              div(style = "display: flex; align-items: center; gap: 6px; padding: 4px 10px; border-radius: 6px; background: linear-gradient(135deg, #e0f2fe, #bae6fd); border: 1px solid #7dd3fc;",
                checkboxInput("xic_show_mobilogram", "Ion Mobility", value = FALSE),
                icon("bolt", style = "color: #0284c7; font-size: 1.1em;")
              )
            }
          )
        )
      ),

      # Mobilogram mode banner — visible when IM is active
      uiOutput("xic_mobilogram_banner"),

      # Alignment guidance banner — visible in alignment mode
      uiOutput("xic_alignment_banner"),

      # Main plot
      plotlyOutput("xic_plot", height = "calc(100vh - 380px)"),

      # Info panel
      div(class = "mt-2 p-2 bg-light rounded",
        style = "font-size: 0.85em;",
        uiOutput("xic_info_panel")
      ),

      size = "xl",
      easyClose = TRUE,
      footer = div(style = "display: flex; gap: 8px; align-items: center;",
        actionButton("xic_prev_protein", label = "Prev",
          icon = icon("arrow-left"), class = "btn-outline-secondary btn-sm"),
        actionButton("xic_next_protein", label = "Next",
          icon = icon("arrow-right"), class = "btn-outline-secondary btn-sm"),
        downloadButton("xic_download_plot", "Download",
          class = "btn-outline-success btn-sm"),
        modalButton("Close")
      )
    ))
  }

  # --- XIC from DE Dashboard ---
  observeEvent(input$show_xic, {
    if (!values$xic_available) {
      showNotification("Load XIC files first (sidebar > '5. XIC Viewer')", type = "warning")
      return()
    }
    if (is.null(values$plot_selected_proteins) || length(values$plot_selected_proteins) == 0) {
      showNotification("Select a protein in the Volcano Plot or Table first!", type = "warning")
      return()
    }
    values$xic_protein <- values$plot_selected_proteins[1]
    values$xic_load_trigger <- Sys.time()
    show_xic_modal(session, values)
  })

  # --- XIC from Grid View ---
  observeEvent(input$show_xic_from_grid, {
    if (!values$xic_available) {
      showNotification("Load XIC files first (sidebar > '5. XIC Viewer')", type = "warning")
      return()
    }
    req(values$grid_selected_protein)
    values$xic_protein <- values$grid_selected_protein
    values$xic_load_trigger <- Sys.time()
    removeModal()
    show_xic_modal(session, values)
  })

  # --- Reactive XIC Data Loading ---
  observe({
    # Depend on trigger so reopening same protein still reloads
    values$xic_load_trigger
    req(values$xic_protein, values$xic_available, values$xic_dir,
        values$xic_report_map)

    withProgress(message = "Loading chromatograms...", {
      incProgress(0.3, detail = "Reading XIC files...")

      tryCatch({
        xic_raw <- load_xic_for_protein(
          xic_dir = values$xic_dir,
          protein_id = values$xic_protein,
          report_map = values$xic_report_map,
          xic_format = values$xic_format
        )

        if (!is.null(xic_raw) && nrow(xic_raw) > 0) {
          incProgress(0.6, detail = "Reshaping data...")

          values$xic_data <- reshape_xic_for_plotting(
            xic_raw, values$metadata, xic_format = values$xic_format)

          # Update precursor selector
          precursors <- unique(values$xic_data$Precursor.Id)
          updateSelectInput(session, "xic_precursor_select",
            choices = c("All Precursors" = "all",
                        setNames(precursors, precursors)),
            selected = "all")

          # Update group filter
          groups <- unique(na.omit(values$xic_data$Group))
          if (length(groups) == 0) groups <- unique(values$metadata$Group)
          updateSelectInput(session, "xic_group_filter",
            choices = c("All Groups" = "all", setNames(groups, groups)),
            selected = "all")
        } else {
          values$xic_data <- NULL
          showNotification(
            paste("No XIC data found for", values$xic_protein,
                  "-- this protein may have been identified via MBR."),
            type = "warning", duration = 5)
        }
      }, error = function(e) {
        values$xic_data <- NULL
        showNotification(paste("Error loading XICs:", e$message),
          type = "error", duration = 8)
      })
    })
  })

  # --- XIC Precursor Count ---
  output$xic_precursor_count <- renderText({
    req(values$xic_data)
    n_prec <- n_distinct(values$xic_data$Precursor.Id)
    n_frag <- values$xic_data %>%
      filter(MS.Level == 2) %>%
      distinct(Fragment.Label) %>%
      nrow()
    paste0(n_prec, " precursor(s), ", n_frag, " fragments")
  })

  # --- Mobilogram Mode Banner ---
  output$xic_mobilogram_banner <- renderUI({
    if (isTRUE(input$xic_show_mobilogram)) {
      div(style = "background: linear-gradient(135deg, #0284c7, #0369a1); color: white; padding: 8px 16px; border-radius: 8px; margin-bottom: 8px; display: flex; align-items: center; gap: 10px; font-weight: 500;",
        icon("bolt", style = "font-size: 1.3em;"),
        span("Ion Mobility Mode", style = "font-size: 1.05em;"),
        span("-- Showing mobilogram data (1/K0 vs intensity). Requires timsTOF or PASEF instrument data.",
          style = "font-weight: 400; font-size: 0.9em; opacity: 0.9;")
      )
    } else {
      NULL
    }
  })

  # --- XIC Alignment Data (shared reactive for stacked bar chart + banner) ---
  xic_alignment_data <- reactive({
    req(values$xic_data, input$xic_display_mode == "alignment")

    xic <- values$xic_data

    # Filter by selected precursor
    if (!is.null(input$xic_precursor_select) &&
        input$xic_precursor_select != "all") {
      xic <- xic %>% filter(Precursor.Id == input$xic_precursor_select)
    }

    # Filter by group
    if (!is.null(input$xic_group_filter) &&
        input$xic_group_filter != "all") {
      xic <- xic %>% filter(Group == input$xic_group_filter)
    }

    # Keep only MS2 fragments
    xic_ms2 <- xic %>% filter(MS.Level == 2)
    if (nrow(xic_ms2) == 0) return(NULL)

    # Compute AUC per (Sample x Fragment) by summing intensity across RT points
    auc_data <- xic_ms2 %>%
      group_by(File.Name, ID, Group, Precursor.Id, Fragment.Label) %>%
      summarise(AUC = sum(Intensity, na.rm = TRUE), .groups = "drop")

    # Compute per-sample total AUC and proportions
    auc_data <- auc_data %>%
      group_by(File.Name, ID, Group) %>%
      mutate(
        Total_AUC = sum(AUC, na.rm = TRUE),
        Proportion = ifelse(Total_AUC > 0, AUC / Total_AUC, 0)
      ) %>%
      ungroup()

    # Compute median proportion per fragment across all samples (reference pattern)
    ref_pattern <- auc_data %>%
      group_by(Fragment.Label) %>%
      summarise(Median_Proportion = median(Proportion, na.rm = TRUE),
                .groups = "drop")

    auc_data <- auc_data %>%
      left_join(ref_pattern, by = "Fragment.Label")

    # Score each sample: sum of absolute deviations from reference
    sample_scores <- auc_data %>%
      group_by(File.Name, ID, Group) %>%
      summarise(
        Deviation_Score = sum(abs(Proportion - Median_Proportion), na.rm = TRUE),
        # Cosine similarity for tooltip
        Cosine_Sim = {
          p <- Proportion
          m <- Median_Proportion
          denom <- sqrt(sum(p^2)) * sqrt(sum(m^2))
          if (denom > 0) sum(p * m) / denom else NA_real_
        },
        .groups = "drop"
      )

    # Flag outliers: deviation > mean + 2*SD
    mean_dev <- mean(sample_scores$Deviation_Score, na.rm = TRUE)
    sd_dev <- sd(sample_scores$Deviation_Score, na.rm = TRUE)
    threshold <- mean_dev + 2 * sd_dev
    # If only 1-2 samples, don't flag (SD is unreliable)
    if (is.na(sd_dev) || nrow(sample_scores) < 3) threshold <- Inf

    sample_scores <- sample_scores %>%
      mutate(Flagged = Deviation_Score > threshold)

    # Add flag info back to auc_data for plotting
    auc_data <- auc_data %>%
      left_join(sample_scores %>% dplyr::select(File.Name, Deviation_Score, Cosine_Sim, Flagged),
                by = "File.Name")

    # Create ordered sample label (by Group then ID)
    sample_order <- auc_data %>%
      distinct(File.Name, ID, Group) %>%
      arrange(Group, ID)
    auc_data$Sample_Label <- factor(
      paste0(auc_data$ID, " (", auc_data$Group, ")"),
      levels = paste0(sample_order$ID, " (", sample_order$Group, ")")
    )

    list(
      auc_data = auc_data,
      sample_scores = sample_scores,
      threshold = threshold,
      ref_pattern = ref_pattern
    )
  })

  # --- XIC Alignment Banner ---
  output$xic_alignment_banner <- renderUI({
    if (!isTRUE(input$xic_display_mode == "alignment")) return(NULL)

    align_data <- xic_alignment_data()
    if (is.null(align_data)) return(NULL)

    scores <- align_data$sample_scores
    flagged <- scores %>% filter(Flagged)
    n_flagged <- nrow(flagged)

    if (n_flagged == 0) {
      div(style = "background: linear-gradient(135deg, #059669, #047857); color: white; padding: 10px 16px; border-radius: 8px; margin-bottom: 8px; display: flex; align-items: center; gap: 10px;",
        icon("circle-check", style = "font-size: 1.3em;"),
        span("All samples consistent", style = "font-weight: 600; font-size: 1.05em;"),
        span(paste0("-- Fragment ion ratios are consistent across all ",
                     nrow(scores), " samples. No interference detected."),
          style = "font-weight: 400; font-size: 0.9em; opacity: 0.9;")
      )
    } else {
      flagged_ids <- paste(flagged$ID, collapse = ", ")
      div(style = "background: linear-gradient(135deg, #d97706, #b45309); color: white; padding: 10px 16px; border-radius: 8px; margin-bottom: 8px;",
        div(style = "display: flex; align-items: center; gap: 10px;",
          icon("triangle-exclamation", style = "font-size: 1.3em;"),
          span(paste0(n_flagged, " sample(s) flagged"),
            style = "font-weight: 600; font-size: 1.05em;"),
          span(paste0("-- Inconsistent fragment ratios in: ", flagged_ids),
            style = "font-weight: 400; font-size: 0.9em; opacity: 0.9;")
        ),
        div(style = "margin-top: 6px; font-size: 0.85em; opacity: 0.9; padding-left: 30px;",
          "Possible causes: co-elution interference, ion suppression, chimeric spectra, or injection issues.",
          " Check chromatogram view for irregular peak shapes in flagged samples."
        )
      )
    }
  })

  # --- XIC Plot ---
  output$xic_plot <- renderPlotly({
    req(values$xic_data)

    xic <- values$xic_data
    display_mode <- input$xic_display_mode

    # --- Alignment mode: short-circuit before chromatogram filtering ---
    if (isTRUE(display_mode == "alignment")) {
      align_data <- xic_alignment_data()
      if (is.null(align_data)) {
        p <- ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "No MS2 fragment data available",
                   size = 5, color = "gray60") +
          theme_void()
        return(ggplotly(p) %>% config(displayModeBar = FALSE))
      }

      auc <- align_data$auc_data
      scores <- align_data$sample_scores

      # Create tooltip text
      auc <- auc %>%
        mutate(
          tooltip_text = paste0(
            "<b>Sample:</b> ", ID, " (", Group, ")",
            "<br><b>Fragment:</b> ", Fragment.Label,
            "<br><b>AUC:</b> ", format(round(AUC), big.mark = ","),
            "<br><b>Proportion:</b> ", round(Proportion * 100, 1), "%",
            "<br><b>Median:</b> ", round(Median_Proportion * 100, 1), "%",
            "<br><b>Deviation score:</b> ", round(Deviation_Score, 3),
            "<br><b>Cosine sim:</b> ", round(Cosine_Sim, 3),
            ifelse(Flagged, "<br><b style='color:#ef4444;'>FLAGGED</b>", "")
          )
        )

      # Identify group boundaries for vertical lines
      sample_groups <- auc %>%
        distinct(Sample_Label, Group) %>%
        arrange(Sample_Label) %>%
        mutate(x_pos = as.numeric(Sample_Label))

      group_boundaries <- sample_groups %>%
        group_by(Group) %>%
        summarise(max_x = max(x_pos), .groups = "drop") %>%
        filter(max_x < max(sample_groups$x_pos))

      # Build plot
      p <- ggplot(auc,
          aes(x = Sample_Label, y = AUC, fill = Fragment.Label,
              text = tooltip_text)) +
        geom_col(position = position_stack(), width = 0.85, color = "white", linewidth = 0.2) +
        scale_y_continuous(labels = scales::label_comma(), expand = c(0, 0)) +
        theme_minimal() +
        labs(
          title = paste("MS2 Intensity Alignment --", values$xic_protein),
          subtitle = "Each bar = one sample. Bar height = total intensity. Consistent proportions = reliable quantification.",
          x = NULL, y = "Summed Fragment Intensity",
          fill = "Fragment"
        ) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          legend.position = "bottom",
          legend.text = element_text(size = 7),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank()
        )

      # Add group boundary lines
      if (nrow(group_boundaries) > 0) {
        p <- p + geom_vline(xintercept = group_boundaries$max_x + 0.5,
                            linetype = "dashed", color = "gray40", linewidth = 0.5)
      }

      # Convert to plotly
      pl <- ggplotly(p, tooltip = "text")

      # Add warning markers for flagged samples
      flagged_samples <- scores %>% filter(Flagged)
      if (nrow(flagged_samples) > 0) {
        flagged_labels <- paste0(flagged_samples$ID, " (", flagged_samples$Group, ")")
        flagged_x <- match(flagged_labels, levels(auc$Sample_Label))
        flagged_x <- flagged_x[!is.na(flagged_x)]

        if (length(flagged_x) > 0) {
          # Position markers just above each bar's total intensity
          sample_totals <- auc %>%
            group_by(Sample_Label) %>%
            summarise(Total = sum(AUC), .groups = "drop")
          annotations <- lapply(flagged_x, function(x) {
            sample_total <- sample_totals$Total[x]
            list(
              x = x - 1,  # plotly uses 0-indexed for categorical
              y = sample_total * 1.03,
              text = "\u26A0",
              showarrow = FALSE,
              font = list(size = 16, color = "#ef4444"),
              xref = "x", yref = "y"
            )
          })
          pl <- pl %>% layout(annotations = annotations)
        }
      }

      return(
        pl %>%
          layout(legend = list(orientation = "h", y = -0.25),
                 margin = list(b = 120)) %>%
          config(displayModeBar = TRUE)
      )
    }

    # --- Chromatogram modes (overlay / facet) ---

    # Filter by selected precursor
    if (!is.null(input$xic_precursor_select) &&
        input$xic_precursor_select != "all") {
      xic <- xic %>% filter(Precursor.Id == input$xic_precursor_select)
    }

    # Filter by group
    if (!is.null(input$xic_group_filter) &&
        input$xic_group_filter != "all") {
      xic <- xic %>% filter(Group == input$xic_group_filter)
    }

    # MS level filter — when MS1 shown, add panel label for split axis
    show_ms1 <- isTRUE(input$xic_show_ms1)
    if (!show_ms1) {
      xic_plot <- xic %>% filter(MS.Level == 2)
      if (nrow(xic_plot) == 0) xic_plot <- xic  # fallback if only MS1 exists
    } else {
      xic_plot <- xic %>%
        mutate(MS_Panel = factor(
          ifelse(MS.Level == 1, "MS1 Precursor", "MS2 Fragments"),
          levels = c("MS1 Precursor", "MS2 Fragments")))
    }

    # Cap at top 6 precursors for very large proteins
    if (n_distinct(xic_plot$Precursor.Id) > 6 &&
        (is.null(input$xic_precursor_select) || input$xic_precursor_select == "all")) {
      top_prec <- xic_plot %>%
        filter(MS.Level == 2) %>%
        count(Precursor.Id) %>%
        slice_max(n, n = 6) %>%
        pull(Precursor.Id)
      # Keep MS1 for selected precursors too
      xic_plot <- xic_plot %>% filter(Precursor.Id %in% top_prec)
      showNotification(
        paste("Showing top 6 of", n_distinct(xic$Precursor.Id),
              "precursors. Use the selector to view others."),
        type = "message", duration = 4)
    }

    # Guard: ensure data exists and has faceting variables before plotting
    if (nrow(xic_plot) == 0 || n_distinct(xic_plot$Fragment.Label) == 0) {
      p <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No chromatogram data available for this selection",
                 size = 5, color = "gray60") +
        theme_void()
      return(ggplotly(p) %>% config(displayModeBar = FALSE))
    }

    # Sample label for faceting
    xic_plot <- xic_plot %>%
      mutate(
        Sample.Label = paste0(ID, " (", Group, ")"),
        tooltip_text = paste0(
          "<b>Sample:</b> ", ID, " (", Group, ")",
          "<br><b>Fragment:</b> ", Fragment.Label,
          "<br><b>RT:</b> ", round(RT, 3), " min",
          "<br><b>Intensity:</b> ", format(round(Intensity), big.mark = ",")
        )
      )

    if (display_mode == "overlay") {
      if (show_ms1) {
        # Split axis: MS1 on top, MS2 fragments below, per sample
        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Fragment.Label,
                group = interaction(File.Name, Fragment.Label),
                text = tooltip_text)) +
          geom_line(alpha = 0.7, linewidth = 0.5) +
          facet_grid(MS_Panel ~ Sample.Label, scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            subtitle = "MS1 precursor (top) vs MS2 fragments (bottom) per sample",
            x = "Retention Time (min)", y = "Intensity", color = "Fragment"
          ) +
          theme(legend.position = "bottom", legend.text = element_text(size = 7),
                strip.text = element_text(size = 8),
                strip.text.y = element_text(angle = 0, face = "bold"))
      } else {
        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Fragment.Label,
                group = interaction(File.Name, Fragment.Label),
                text = tooltip_text)) +
          geom_line(alpha = 0.7, linewidth = 0.5) +
          facet_wrap(~ Sample.Label, scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            subtitle = "Each panel = one sample, colors = fragment ions",
            x = "Retention Time (min)", y = "Intensity", color = "Fragment"
          ) +
          theme(legend.position = "bottom", legend.text = element_text(size = 7),
                strip.text = element_text(size = 8))
      }

    } else if (display_mode == "facet") {
      if (show_ms1) {
        # Split axis: MS1 on top row, each fragment in its own panel below
        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Group,
                group = File.Name, text = tooltip_text)) +
          geom_line(alpha = 0.6, linewidth = 0.4) +
          facet_grid(MS_Panel ~ Fragment.Label, scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            subtitle = "MS1 precursor (top) vs MS2 fragments (bottom), colors = groups",
            x = "Retention Time (min)", y = "Intensity", color = "Group"
          ) +
          theme(strip.text = element_text(size = 7),
                strip.text.y = element_text(angle = 0, face = "bold"))
      } else {
        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Group,
                group = File.Name, text = tooltip_text)) +
          geom_line(alpha = 0.6, linewidth = 0.4) +
          facet_wrap(~Fragment.Label, scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            subtitle = "Each panel = one fragment ion, colors = groups",
            x = "Retention Time (min)", y = "Intensity", color = "Group"
          ) +
          theme(strip.text = element_text(size = 8))
      }

    }

    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "h", y = -0.15),
             margin = list(b = 80)) %>%
      config(displayModeBar = TRUE)
  })

  # --- XIC Info Panel ---
  output$xic_info_panel <- renderUI({
    req(values$xic_data, values$xic_protein)

    xic <- values$xic_data
    n_precursors <- n_distinct(xic$Precursor.Id)
    n_fragments <- xic %>% filter(MS.Level == 2) %>%
      distinct(Fragment.Label) %>% nrow()
    n_samples <- n_distinct(xic$File.Name)
    rt_range <- range(xic$RT, na.rm = TRUE)

    # Get DE stats for this protein from current contrast
    de_info <- tryCatch({
      tt <- topTable(values$fit, coef = input$contrast_selector, number = Inf)
      tt[values$xic_protein, c("logFC", "adj.P.Val")]
    }, error = function(e) NULL)

    tagList(
      div(class = "row",
        div(class = "col-md-4",
          strong("Protein: "), values$xic_protein, br(),
          strong("Precursors: "), n_precursors, br(),
          strong("Fragment ions: "), n_fragments
        ),
        div(class = "col-md-4",
          strong("Samples: "), n_samples, br(),
          strong("RT range: "),
          paste0(round(rt_range[1], 2), " -- ", round(rt_range[2], 2), " min")
        ),
        div(class = "col-md-4",
          if (!is.null(de_info)) {
            tagList(
              strong("log2 FC: "),
              span(round(de_info$logFC, 3),
                style = paste0("color: ",
                  ifelse(de_info$logFC > 0, "#d32f2f", "#1976d2"), ";")),
              br(),
              strong("adj. p-value: "),
              formatC(de_info$adj.P.Val, format = "e", digits = 2)
            )
          } else {
            em("DE stats unavailable for current contrast")
          }
        )
      ),
      hr(style = "margin: 5px 0;"),
      if (isTRUE(input$xic_display_mode == "alignment")) {
        div(class = "text-muted", style = "font-size: 0.8em;",
          icon("chart-bar"),
          " Each bar shows relative fragment ion proportions for one sample.",
          " Consistent bars across samples indicate reliable quantification.",
          " Bars with unusual proportions (\u26A0) may have co-elution interference, ion suppression, or chimeric spectra.",
          " Group separators (dashed lines) help compare within and between conditions."
        )
      } else {
        div(class = "text-muted", style = "font-size: 0.8em;",
          icon("info-circle"),
          " Co-eluting fragment ions with similar peak shapes indicate reliable identification.",
          " Consistent peak areas across replicates within a group support accurate quantification.",
          " Irregular peaks or missing fragments may indicate interference or low-confidence IDs."
        )
      }
    )
  })

  # --- XIC Protein Navigation (Prev/Next) ---
  observeEvent(input$xic_prev_protein, {
    req(values$xic_protein, values$fit, input$contrast_selector)
    tt <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>%
      filter(adj.P.Val < 0.05) %>% arrange(adj.P.Val)
    prot_list <- rownames(tt)
    current_idx <- which(prot_list == values$xic_protein)
    if (length(current_idx) > 0 && current_idx > 1) {
      values$xic_protein <- prot_list[current_idx - 1]
      values$plot_selected_proteins <- values$xic_protein
    }
  })

  observeEvent(input$xic_next_protein, {
    req(values$xic_protein, values$fit, input$contrast_selector)
    tt <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>%
      filter(adj.P.Val < 0.05) %>% arrange(adj.P.Val)
    prot_list <- rownames(tt)
    current_idx <- which(prot_list == values$xic_protein)
    if (length(current_idx) > 0 && current_idx < length(prot_list)) {
      values$xic_protein <- prot_list[current_idx + 1]
      values$plot_selected_proteins <- values$xic_protein
    }
  })

  # --- XIC Download ---
  output$xic_download_plot <- downloadHandler(
    filename = function() {
      mode_tag <- if (isTRUE(input$xic_display_mode == "alignment")) "Alignment" else "XIC"
      paste0(mode_tag, "_", gsub("[^A-Za-z0-9]", "_", values$xic_protein), "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(values$xic_data)

      if (isTRUE(input$xic_display_mode == "alignment")) {
        # --- Alignment mode: stacked bar chart ---
        align_data <- xic_alignment_data()
        req(align_data)
        auc <- align_data$auc_data
        scores <- align_data$sample_scores
        group_boundaries_df <- auc %>%
          distinct(Sample_Label, Group) %>%
          arrange(Sample_Label) %>%
          mutate(x_pos = as.numeric(Sample_Label)) %>%
          group_by(Group) %>%
          summarise(max_x = max(x_pos), .groups = "drop") %>%
          filter(max_x < max(as.numeric(auc$Sample_Label)))

        p <- ggplot(auc,
            aes(x = Sample_Label, y = AUC, fill = Fragment.Label)) +
          geom_col(position = position_stack(), width = 0.85, color = "white", linewidth = 0.2) +
          scale_y_continuous(labels = scales::label_comma(), expand = c(0, 0)) +
          theme_minimal() +
          labs(
            title = paste("MS2 Intensity Alignment --", values$xic_protein),
            subtitle = paste0("Flagged: ", sum(scores$Flagged), " of ", nrow(scores), " samples"),
            x = NULL, y = "Summed Fragment Intensity", fill = "Fragment"
          ) +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "bottom",
            panel.grid.major.x = element_blank(),
            panel.grid.minor = element_blank()
          )

        if (nrow(group_boundaries_df) > 0) {
          p <- p + geom_vline(xintercept = group_boundaries_df$max_x + 0.5,
                              linetype = "dashed", color = "gray40", linewidth = 0.5)
        }

        # Mark flagged samples with red triangles
        flagged <- scores %>% filter(Flagged)
        if (nrow(flagged) > 0) {
          flagged_labels <- paste0(flagged$ID, " (", flagged$Group, ")")
          sample_totals <- auc %>%
            group_by(Sample_Label) %>%
            summarise(Total = sum(AUC), .groups = "drop")
          flagged_y <- sample_totals$Total[match(flagged_labels, sample_totals$Sample_Label)]
          flagged_y[is.na(flagged_y)] <- max(sample_totals$Total)
          p <- p + annotate("text", x = flagged_labels, y = flagged_y * 1.03,
                            label = "\u26A0", color = "#ef4444", size = 5)
        }

        ggsave(file, plot = p, width = 14, height = 10, dpi = 150)

      } else {
        # --- Chromatogram mode ---
        xic <- values$xic_data

        if (!isTRUE(input$xic_show_ms1)) {
          xic_plot <- xic %>% filter(MS.Level == 2)
          if (nrow(xic_plot) == 0) xic_plot <- xic
        } else {
          xic_plot <- xic
        }

        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Fragment.Label,
                group = interaction(File.Name, Fragment.Label))) +
          geom_line(alpha = 0.7, linewidth = 0.5) +
          facet_wrap(~ paste0(ID, " (", Group, ")"), scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            x = "Retention Time (min)", y = "Intensity", color = "Fragment"
          ) +
          theme(legend.position = "bottom", legend.text = element_text(size = 7),
                strip.text = element_text(size = 8))

        ggsave(file, plot = p, width = 14, height = 10, dpi = 150)
      }
    }
  )
}
