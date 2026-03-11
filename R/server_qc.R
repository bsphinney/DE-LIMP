# ==============================================================================
#  SERVER MODULE — QC Trends, QC Diagnostic Plots
#  Called from app.R as: server_qc(input, output, session, values)
# ==============================================================================

server_qc <- function(input, output, session, values) {

  # ============================================================================
  #  1. QC Trend Plots (Precursors, Proteins, MS1 Signal)
  # ============================================================================

  # Shared reactive: faceted QC metrics data (Precursors, Proteins, MS1 Signal, Data Completeness)
  qc_metrics_data <- reactive({
    req(values$qc_stats, values$metadata)
    df <- left_join(values$qc_stats, values$metadata, by = c("Run" = "File.Name")) %>%
      mutate(Run_Number = as.numeric(str_extract(Run, "\\d+$")))

    # Compute per-sample data completeness from precursor matrix
    if (!is.null(values$raw_data)) {
      raw_mat <- values$raw_data$E
      completeness <- colMeans(!is.na(raw_mat)) * 100  # % non-NA per sample
      df$Completeness <- completeness[match(df$Run, names(completeness))]
    } else {
      df$Completeness <- NA_real_
    }

    if (input$qc_sort_order == "Group") {
      df <- df %>% arrange(Group, Run_Number)
    } else {
      df <- df %>% arrange(Run_Number)
    }
    df$Sort_Index <- 1:nrow(df)

    # Pivot to long format for faceting
    pivot_cols <- c("Precursors", "Proteins", "MS1_Signal")
    metric_levels <- c("Precursors", "Proteins", "MS1_Signal")
    metric_labels <- c("Precursors", "Proteins", "MS1 Signal")
    if (!all(is.na(df$Completeness))) {
      pivot_cols <- c(pivot_cols, "Completeness")
      metric_levels <- c(metric_levels, "Completeness")
      metric_labels <- c(metric_labels, "Data Completeness (%)")
    }

    long <- df %>%
      pivot_longer(
        cols = all_of(pivot_cols),
        names_to = "Metric",
        values_to = "Value"
      ) %>%
      mutate(
        Metric = factor(Metric, levels = metric_levels, labels = metric_labels),
        Tooltip = paste0("<b>File:</b> ", Run, "<br><b>Group:</b> ", Group,
                         "<br><b>", Metric, ":</b> ", round(Value, 2))
      )

    # Group averages per metric per group
    group_stats <- long %>%
      group_by(Group, Metric) %>%
      summarise(
        mean_value = mean(Value, na.rm = TRUE),
        x_min = min(Sort_Index),
        x_max = max(Sort_Index),
        .groups = "drop"
      )

    list(long = long, group_stats = group_stats)
  })

  # Build faceted trend plot from prepared data
  build_qc_metrics_plot <- function(data) {
    long <- data$long
    group_stats <- data$group_stats

    # Split data: bars for counts/intensity, dots for completeness
    is_completeness <- "Data Completeness (%)"
    bar_data <- long %>% filter(Metric != is_completeness)
    dot_data <- long %>% filter(Metric == is_completeness)
    bar_stats <- group_stats %>% filter(Metric != is_completeness)
    dot_stats <- group_stats %>% filter(Metric == is_completeness)

    p <- ggplot(long, aes(x = Sort_Index, y = Value, text = Tooltip)) +
      geom_bar(data = bar_data, aes(fill = Group),
               stat = "identity", width = 0.8) +
      geom_point(data = dot_data, aes(color = Group),
                 size = 3, alpha = 0.8) +
      geom_segment(data = bar_stats,
                   aes(x = x_min - 0.5, xend = x_max + 0.5,
                       y = mean_value, yend = mean_value, color = Group),
                   linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE,
                   show.legend = FALSE) +
      geom_segment(data = dot_stats,
                   aes(x = x_min - 0.5, xend = x_max + 0.5,
                       y = mean_value, yend = mean_value, color = Group),
                   linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE,
                   show.legend = FALSE) +
      geom_smooth(aes(group = 1), method = "loess", se = FALSE,
                  color = "black", linewidth = 1, span = 0.75) +
      facet_wrap(~ Metric, ncol = 1, scales = "free_y",
                 strip.position = "top") +
      scale_color_discrete(guide = "none") +
      theme_minimal() +
      labs(x = "Sample Index (Sorted)", y = NULL) +
      theme(
        panel.grid.major.x = element_blank(),
        axis.text.x = element_text(size = 8),
        strip.background = element_rect(fill = "#2c3e50", color = NA),
        strip.text = element_text(color = "white", face = "bold", size = 11)
      )

    ggplotly(p, tooltip = "text") %>% config(displayModeBar = TRUE)
  }

  # Main faceted trend plot
  output$qc_metrics_trend <- renderPlotly({
    build_qc_metrics_plot(qc_metrics_data())
  })

  # Fullscreen modal
  observeEvent(input$fullscreen_qc_metrics, {
    showModal(modalDialog(
      title = "Sample Metrics - Fullscreen View",
      plotlyOutput("qc_metrics_trend_fs", height = "800px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$qc_metrics_trend_fs <- renderPlotly({
    build_qc_metrics_plot(qc_metrics_data())
  })

  # ============================================================================
  #  2. QC Stats Table
  # ============================================================================

  output$r_qc_table <- renderDT({ req(values$qc_stats); df_display <- values$qc_stats %>% arrange(Run) %>% mutate(ID = 1:n()) %>% dplyr::select(ID, Run, everything()); datatable(df_display, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) })

  # QC Stats CSV export
  output$download_qc_stats_csv <- downloadHandler(
    filename = function() {
      paste0("QC_Stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$qc_stats)
      write.csv(values$qc_stats %>% arrange(Run), file, row.names = FALSE)
    }
  )

  # QC Stats info modal
  observeEvent(input$qc_stats_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " QC Statistics Table"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("Per-run QC statistics extracted from your DIA-NN report: precursor counts, protein counts, and MS1 signal intensity."),
        p("Use this table to identify outlier runs with unusually low precursor/protein counts or signal intensity. ",
          "Export to CSV for external QC tracking or reporting.")
      )
    ))
  })

  # Sample Metrics info modal (combined)
  observeEvent(input$qc_metrics_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " Sample Metrics"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("This faceted plot shows key quality metrics for every run in your experiment, ",
          "stacked vertically so you can spot correlated trends at a glance."),
        tags$h6("The metrics"),
        tags$ul(
          tags$li(strong("Precursors: "), "Number of peptide precursors identified at your Q-value cutoff. ",
            "Higher counts indicate better instrument sensitivity and sample quality."),
          tags$li(strong("Proteins: "), "Number of protein groups quantified per run. ",
            "Should be relatively stable across runs; lower counts may indicate sample quality issues."),
          tags$li(strong("MS1 Signal: "), "Overall MS1 intensity per run. ",
            "Consistent signal indicates stable instrument performance and uniform sample loading."),
          tags$li(strong("Data Completeness (%): "), "Percentage of precursors detected (non-missing) per sample ",
            "in the raw expression matrix. Samples with low completeness had many missed detections, ",
            "which may indicate injection failures or low sample loading.")
        ),
        tags$h6("Reading the plot"),
        tags$ul(
          tags$li(strong("Bars: "), "Per-run values, colored by experimental group."),
          tags$li(strong("Dashed lines: "), "Group averages \u2014 compare groups at a glance."),
          tags$li(strong("Black trend line (LOESS): "), "Smoothed trend across all runs. ",
            "A flat line means stable performance; a downward slope suggests instrument drift ",
            "(e.g., column degradation, source contamination).")
        ),
        tags$h6("What to look for"),
        tags$ul(
          tags$li("Sudden drops in a single sample flag potential injection failures or outliers"),
          tags$li("Gradual downward trends across all metrics suggest instrument degradation"),
          tags$li("Use ", strong("Sort Order: Run Order"), " to see acquisition-time drift; ",
            strong("Group"), " to compare conditions side by side")
        )
      )
    ))
  })

  # ============================================================================
  #  3. Group QC Violin Plot
  # ============================================================================

  output$qc_group_violin <- renderPlotly({
    req(values$qc_stats, values$metadata, input$qc_violin_metric)
    df <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")); metric <- input$qc_violin_metric
    df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Val:</b> ", round(df[[metric]], 2))
    p <- ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) + geom_violin(alpha = 0.5, trim = FALSE) + geom_jitter(aes(text = Tooltip), width = 0.2, size = 2, alpha = 0.8, color = "black") + theme_bw() + labs(title = paste("Distribution of", metric), x = "Group", y = metric) + theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })

  # ============================================================================
  #  4. DPC Plot
  # ============================================================================

  output$dpc_plot <- renderPlot({ req(values$dpc_fit); limpa::plotDPC(values$dpc_fit) }) # Height controlled by UI (70vh)

  # ============================================================================
  #  5. MDS Plot
  # ============================================================================

  # Update MDS "Color by" dropdown when metadata changes (includes custom covariate names)
  observeEvent(values$metadata, {
    req(values$metadata)
    meta <- values$metadata
    color_choices <- "Group"
    if ("Batch" %in% colnames(meta) && any(nzchar(meta$Batch))) color_choices <- c(color_choices, "Batch")
    # Add custom covariates if they have data
    cov1_name <- if (!is.null(values$cov1_name) && nzchar(values$cov1_name)) values$cov1_name else "Covariate1"
    cov2_name <- if (!is.null(values$cov2_name) && nzchar(values$cov2_name)) values$cov2_name else "Covariate2"
    if ("Covariate1" %in% colnames(meta) && any(nzchar(meta$Covariate1))) {
      color_choices <- c(color_choices, setNames("Covariate1", cov1_name))
    }
    if ("Covariate2" %in% colnames(meta) && any(nzchar(meta$Covariate2))) {
      color_choices <- c(color_choices, setNames("Covariate2", cov2_name))
    }
    updateSelectInput(session, "mds_color_by", choices = color_choices, selected = "Group")
  })

  # Helper: get MDS color variable from metadata
  mds_color_data <- function(meta) {
    color_by <- input$mds_color_by %||% "Group"
    col_name <- if (color_by %in% colnames(meta)) color_by else "Group"
    vals <- meta[[col_name]]
    vals[is.na(vals) | vals == ""] <- "(unassigned)"
    grps <- factor(vals)
    # Use a colorblind-friendly palette (up to 12 levels, then fall back to rainbow)
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7", "#999999", "#000000", "#66A61E", "#E6AB02", "#A6761D")
    n_lvl <- length(levels(grps))
    cols <- if (n_lvl <= length(palette)) palette[1:n_lvl] else rainbow(n_lvl)
    # Build label for legend header
    label <- color_by
    if (color_by == "Covariate1" && !is.null(values$cov1_name) && nzchar(values$cov1_name)) label <- values$cov1_name
    if (color_by == "Covariate2" && !is.null(values$cov2_name) && nzchar(values$cov2_name)) label <- values$cov2_name
    list(grps = grps, cols = cols, label = label)
  }

  output$mds_plot <- renderPlot({
    req(values$y_protein, values$metadata)
    meta <- values$metadata[match(colnames(values$y_protein$E), values$metadata$File.Name), ]
    cd <- mds_color_data(meta)
    limpa::plotMDSUsingSEs(values$y_protein, pch = 16,
      main = paste0("MDS Plot (colored by ", cd$label, ")"), col = cd$cols[cd$grps])
    legend("bottomright", legend = levels(cd$grps), col = cd$cols[1:length(levels(cd$grps))],
           pch = 16, bg = "white", box.col = "gray80", cex = 0.9, title = cd$label)
  }) # Height controlled by UI (70vh)

  # ============================================================================
  #  6. Normalization Diagnostic
  # ============================================================================

  # DIA-NN normalization status badge
  output$diann_norm_status_badge <- renderUI({
    status <- values$diann_norm_detected
    if (status == "on") {
      span(class = "badge bg-info", style = "margin-right: 10px;",
        icon("check-circle"), " DIA-NN normalization: ON (RT-dependent)")
    } else if (status == "off") {
      span(class = "badge bg-warning", style = "margin-right: 10px;",
        icon("exclamation-triangle"), " DIA-NN normalization: OFF")
    } else {
      span(class = "badge bg-secondary", style = "margin-right: 10px;",
        icon("question-circle"), " DIA-NN normalization: unknown")
    }
  })

  # Health assessment reactive
  assess_distribution_health <- reactive({
    req(values$raw_data, values$y_protein)
    pre_mat <- values$raw_data$E
    post_mat <- values$y_protein$E
    pre_medians <- apply(pre_mat, 2, median, na.rm = TRUE)
    post_medians <- apply(post_mat, 2, median, na.rm = TRUE)
    pre_cv <- sd(pre_medians) / abs(mean(pre_medians))
    post_cv <- sd(post_medians) / abs(mean(post_medians))
    pre_outliers <- which(abs(pre_medians - mean(pre_medians)) > 2 * sd(pre_medians))
    post_outliers <- which(abs(post_medians - mean(post_medians)) > 2 * sd(post_medians))
    list(
      pre_cv = pre_cv, post_cv = post_cv,
      pre_outlier_samples = names(pre_outliers),
      post_outlier_samples = names(post_outliers),
      status = if (pre_cv > 0.05) "warning" else if (length(pre_outliers) > 0) "caution" else "good"
    )
  })

  # Contextual guidance banners
  output$norm_diag_guidance <- renderUI({
    req(values$raw_data, values$y_protein)
    health <- assess_distribution_health()
    diann_status <- values$diann_norm_detected
    warnings <- list()

    # SCENARIO 1: DIA-NN normalization OFF + bad distributions
    if (diann_status == "off" && health$status == "warning") {
      warnings <- c(warnings, list(
        div(class = "alert alert-warning", role = "alert",
          icon("exclamation-triangle"),
          strong(" Unnormalized data detected. "),
          "Your DIA-NN output does not appear to be normalized (the ",
          tags$code("Precursor.Normalised"), " and ", tags$code("Precursor.Quantity"),
          " columns are identical). The sample distributions look uneven, which can lead ",
          "to unreliable differential expression results.",
          br(), br(),
          strong("What to do: "),
          "For most experiments, re-process your data in DIA-NN with ",
          tags$b("RT-dependent normalization"), " enabled (this is the default setting). ",
          "This corrects for differences in sample loading and LC-MS run variability.",
          br(), br(),
          em("Exception: "), "If you are analyzing AP-MS/Co-IP, fractionated samples, or ",
          "isotope labeling time-courses, unnormalized data may be appropriate \u2014 ",
          "but you should apply your own normalization before using DE-LIMP."
        )
      ))
    }

    # SCENARIO 2: DIA-NN normalization OFF but distributions look OK
    if (diann_status == "off" && health$status == "good") {
      warnings <- c(warnings, list(
        div(class = "alert alert-info", role = "alert",
          icon("info-circle"),
          strong(" DIA-NN normalization was off, "),
          "but your sample distributions look reasonably aligned. This can happen if your ",
          "samples had very consistent loading and LC-MS performance. Results may still be ",
          "valid, but consider whether normalization would improve your analysis."
        )
      ))
    }

    # SCENARIO 3: DIA-NN normalization ON but distributions still bad
    if (diann_status == "on" && health$status == "warning") {
      warnings <- c(warnings, list(
        div(class = "alert alert-danger", role = "alert",
          icon("times-circle"),
          strong(" Sample distributions are uneven despite normalization. "),
          "DIA-NN normalization was applied but the per-sample distributions still show ",
          "substantial differences. This could indicate:",
          tags$ul(
            tags$li("A failed or low-quality injection (check the QC Trends tab)"),
            tags$li("Very different sample types being compared (e.g., tissue vs plasma)"),
            tags$li("Severe batch effects that normalization couldn't fully correct")
          ),
          strong("What to do: "),
          "Check the QC Trends tab for outlier samples. Consider whether any samples ",
          "should be excluded. If batch effects are suspected, make sure you've assigned ",
          "batch information in the Assign Groups modal."
        )
      ))
    }

    # SCENARIO 4: Outlier sample(s) detected
    if (length(health$pre_outlier_samples) > 0) {
      outlier_names <- paste(health$pre_outlier_samples, collapse = ", ")
      warnings <- c(warnings, list(
        div(class = "alert alert-warning", role = "alert",
          icon("user-times"),
          strong(" Possible outlier sample(s): "),
          tags$code(outlier_names),
          br(),
          "These samples have median intensities substantially different from the rest. ",
          "This could indicate a failed injection, sample preparation issue, or ",
          "biological outlier. Check the QC Trends tab and MDS plot for confirmation. ",
          "If the sample is clearly problematic, consider re-running the analysis ",
          "without it."
        )
      ))
    }

    # SCENARIO 5: Everything looks good
    if (health$status == "good" && diann_status %in% c("on", "unknown") &&
        length(health$pre_outlier_samples) == 0) {
      warnings <- c(warnings, list(
        div(class = "alert alert-success", role = "alert",
          icon("check-circle"),
          strong(" Distributions look good. "),
          "Per-sample intensity distributions are well-aligned. ",
          "No outlier samples detected."
        )
      ))
    }

    do.call(tagList, warnings)
  })

  # Shared reactive for the diagnostic plot (regular + fullscreen)
  generate_norm_diagnostic_plot <- reactive({
    req(values$raw_data, values$y_protein, values$metadata)

    pre_mat <- values$raw_data$E   # precursor-level, log2, NAs present
    post_mat <- values$y_protein$E # protein-level, log2, no NAs
    meta <- values$metadata

    if (input$norm_diag_type == "boxplot") {
      # === BOX PLOT VIEW ===
      # Subsample precursor matrix if very large (performance)
      if (nrow(pre_mat) > 10000) {
        sample_idx <- sample(nrow(pre_mat), 10000)
        pre_mat_plot <- pre_mat[sample_idx, ]
      } else {
        pre_mat_plot <- pre_mat
      }

      pre_long <- as.data.frame(pre_mat_plot) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
        mutate(Stage = "Precursor Input\n(DIA-NN normalized)") %>%
        filter(!is.na(Log2Intensity))

      post_long <- as.data.frame(post_mat) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
        mutate(Stage = "Protein Output\n(DPC-Quant)")

      pre_long$Group <- meta$Group[match(pre_long$Sample, meta$File.Name)]
      post_long$Group <- meta$Group[match(post_long$Sample, meta$File.Name)]

      combined <- bind_rows(pre_long, post_long)
      combined$Stage <- factor(combined$Stage,
        levels = c("Precursor Input\n(DIA-NN normalized)", "Protein Output\n(DPC-Quant)"))

      sample_order <- meta %>% arrange(Group, File.Name) %>% pull(File.Name)
      combined$Sample <- factor(combined$Sample, levels = sample_order)
      combined$SampleID <- meta$ID[match(combined$Sample, meta$File.Name)]

      p <- ggplot(combined, aes(x = factor(SampleID), y = Log2Intensity, fill = Group)) +
        geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
        facet_wrap(~Stage, scales = "free_y", ncol = 2) +
        theme_minimal() +
        labs(
          title = "Pipeline Diagnostic: Precursor Input \u2192 Protein Output",
          subtitle = "Left: DIA-NN normalized precursors | Right: DPC-Quant protein estimates",
          x = "Sample ID", y = "Log2 Intensity"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

      ggplotly(p, tooltip = c("x", "y")) %>% layout(boxmode = "group")

    } else {
      # === DENSITY OVERLAY VIEW ===
      pre_long <- as.data.frame(pre_mat) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
        mutate(Stage = "Precursor Input (DIA-NN)") %>%
        filter(!is.na(Log2Intensity))

      post_long <- as.data.frame(post_mat) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
        mutate(Stage = "Protein Output (DPC-Quant)")

      pre_long$Group <- meta$Group[match(pre_long$Sample, meta$File.Name)]
      post_long$Group <- meta$Group[match(post_long$Sample, meta$File.Name)]

      combined <- bind_rows(pre_long, post_long)
      combined$Stage <- factor(combined$Stage,
        levels = c("Precursor Input (DIA-NN)", "Protein Output (DPC-Quant)"))

      p <- ggplot(combined, aes(x = Log2Intensity, color = Group, group = Sample)) +
        geom_density(alpha = 0.3, linewidth = 0.4) +
        facet_wrap(~Stage, ncol = 2) +
        theme_minimal() +
        labs(
          title = "Pipeline Diagnostic: Per-Sample Density Curves",
          subtitle = "Well-aligned input distributions should remain aligned after quantification",
          x = "Log2 Intensity", y = "Density"
        )

      ggplotly(p)
    }
  })

  output$norm_diagnostic_plot <- renderPlotly({ generate_norm_diagnostic_plot() })

  # Fullscreen modal for pipeline diagnostic
  observeEvent(input$fullscreen_norm_diag, {
    showModal(modalDialog(
      title = "Pipeline Diagnostic - Fullscreen View",
      plotlyOutput("norm_diagnostic_plot_fullscreen", height = "700px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$norm_diagnostic_plot_fullscreen <- renderPlotly({ generate_norm_diagnostic_plot() })

  observeEvent(input$norm_diag_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What am I looking at?"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Reading this plot"),
        p("Each box (or density curve) represents one sample's intensity distribution \u2014 ",
          "essentially, how bright all the detected peptides/proteins are in that sample."),
        p(strong("Left panel: "), "What DIA-NN gave us. These are the peptide-level intensities ",
          "after DIA-NN's normalization (if it was enabled). ",
          strong("Right panel: "), "What our pipeline produced. These are the final protein-level ",
          "estimates after aggregating peptides and handling missing values."),
        tags$h6("What 'good' looks like"),
        p("The boxes (or curves) should sit at roughly the same height across all samples. ",
          "Small differences are normal. If one sample is dramatically higher or lower than ",
          "the rest, that sample may be problematic."),
        tags$h6("What 'bad' looks like"),
        p("If all the boxes are at very different heights, your samples aren't comparable ",
          "and the statistical results may not be reliable. The most common cause is that ",
          "DIA-NN normalization was turned off when the data was processed."),
        tags$h6("Why doesn't the right panel 'fix' bad data?"),
        p("Unlike some other tools, this pipeline does not apply its own normalization. ",
          "The protein quantification step (DPC-Quant) aggregates peptides into proteins and ",
          "handles missing values, but it ",
          strong("trusts the input intensities as-is"), ". ",
          "If the input is unnormalized, the output will be too. ",
          "Normalization happens in DIA-NN, before the data reaches this tool."),
        tags$h6("The DIA-NN normalization badge"),
        p("The badge next to the plot controls tells you whether DIA-NN applied normalization to your data. ",
          tags$span(class = "badge bg-info", "ON"), " = DIA-NN's RT-dependent normalization was active (recommended). ",
          tags$span(class = "badge bg-warning", "OFF"), " = Data was exported without normalization. ",
          tags$span(class = "badge bg-secondary", "Unknown"), " = Couldn't determine (older DIA-NN version or non-standard export).")
      )
    ))
  })

  # --- DPC Fit Info Modal ---
  observeEvent(input$dpc_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What is DPC Fit?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Data Point Correspondence (DPC)"),
        p("DPC is the normalization and quantification method used by the LIMPA pipeline. ",
          "It models the relationship between peptide-level measurements and protein-level estimates, ",
          "accounting for missing values and variable peptide behavior."),
        tags$h6("What this plot shows"),
        p("The DPC fit plot visualizes how well the model fits your data. Each point represents a peptide-protein ",
          "relationship, and the fitted curve shows the expected correspondence."),
        tags$h6("What 'good' looks like"),
        tags$ul(
          tags$li("Points should cluster tightly around the fitted line"),
          tags$li("No strong systematic deviations or outlier clusters"),
          tags$li("The fit should be smooth without sharp jumps")
        ),
        tags$h6("What 'bad' looks like"),
        tags$ul(
          tags$li("Large scatter around the fitted line suggests noisy data or poor peptide-to-protein mapping"),
          tags$li("Systematic curvature away from the fit may indicate batch effects or normalization issues"),
          tags$li("Distinct outlier clusters could indicate contaminated samples or misassigned peptides")
        )
      )
    ))
  })

  # --- MDS Plot Info Modal ---
  observeEvent(input$mds_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What is the MDS Plot?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Multidimensional Scaling (MDS)"),
        p("MDS reduces the high-dimensional protein expression data into two dimensions so you can ",
          "visualize how similar or different your samples are. Think of it as a map where samples ",
          "that are biologically similar appear close together."),
        tags$h6("What 'good' looks like"),
        tags$ul(
          tags$li("Samples from the same experimental group cluster together"),
          tags$li("Different groups are clearly separated"),
          tags$li("Replicates within a group are tightly clustered")
        ),
        tags$h6("What 'bad' looks like"),
        tags$ul(
          tags$li(strong("One sample far from its group: "), "Possible outlier \u2014 check sample quality, injection issues, or mislabeling"),
          tags$li(strong("Groups overlap completely: "), "Little biological difference between conditions, or high technical variability masking real signal"),
          tags$li(strong("Samples cluster by batch, not group: "), "Batch effect \u2014 consider adding batch as a covariate in the model")
        ),
        tags$h6("Reading the axes"),
        p("The axes show 'leading z-statistic dimensions' with the percentage of variance explained in parentheses. ",
          "Dimension 1 (x-axis) captures the largest source of variation, dimension 2 (y-axis) the second largest. ",
          "High percentage on dim 1 means most variation is along that axis.")
      )
    ))
  })

  # --- Group Distribution Info Modal ---
  observeEvent(input$group_dist_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What is the Group Distribution?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Group-level QC Violin Plots"),
        p("These violin plots show the distribution of a QC metric across your experimental groups. ",
          "The width of the violin indicates how many samples have that value \u2014 wider means more samples."),
        tags$h6("Available metrics"),
        tags$ul(
          tags$li(strong("Precursors: "), "Number of peptide precursors identified per sample. More = better sensitivity."),
          tags$li(strong("Proteins: "), "Number of proteins quantified per sample. Should be consistent across groups."),
          tags$li(strong("MS1 Signal: "), "Overall MS1 intensity. Large differences may indicate loading or injection issues.")
        ),
        tags$h6("What to look for"),
        tags$ul(
          tags$li("Groups should have similar distributions (overlapping violins)"),
          tags$li("A group with consistently lower values may have systematic quality issues"),
          tags$li("Individual outlier dots indicate samples worth investigating in more detail")
        )
      )
    ))
  })

  # ============================================================================
  #  Fullscreen Modals for QC Plot Panels
  # ============================================================================

  # --- DPC Fit (QC Plots) ---
  observeEvent(input$fullscreen_dpc, {
    showModal(modalDialog(
      title = "DPC Fit - Fullscreen View",
      plotOutput("dpc_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$dpc_plot_fs <- renderPlot({ req(values$dpc_fit); limpa::plotDPC(values$dpc_fit) }, height = 700)

  # --- MDS Plot (QC Plots) ---
  observeEvent(input$fullscreen_mds, {
    showModal(modalDialog(
      title = "MDS Plot - Fullscreen View",
      plotOutput("mds_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$mds_plot_fs <- renderPlot({
    req(values$y_protein, values$metadata)
    meta <- values$metadata[match(colnames(values$y_protein$E), values$metadata$File.Name), ]
    cd <- mds_color_data(meta)
    limpa::plotMDSUsingSEs(values$y_protein, pch = 16,
      main = paste0("MDS Plot (colored by ", cd$label, ")"), col = cd$cols[cd$grps])
    legend("bottomright", legend = levels(cd$grps), col = cd$cols[1:length(levels(cd$grps))],
           pch = 16, bg = "white", box.col = "gray80", cex = 0.9, title = cd$label)
  }, height = 700)

  # --- Group QC Distribution Violin (QC Plots) ---
  observeEvent(input$fullscreen_qc_violin, {
    showModal(modalDialog(
      title = "Group QC Distribution - Fullscreen View",
      plotlyOutput("qc_group_violin_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$qc_group_violin_fs <- renderPlotly({
    req(values$qc_stats, values$metadata, input$qc_violin_metric)
    df <- left_join(values$qc_stats, values$metadata, by = c("Run" = "File.Name"))
    metric <- input$qc_violin_metric
    df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Val:</b> ", round(df[[metric]], 2))
    p <- ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) +
      geom_violin(alpha = 0.5, trim = FALSE) +
      geom_jitter(aes(text = Tooltip), width = 0.2, size = 2, alpha = 0.8, color = "black") +
      theme_bw() + labs(title = paste("Distribution of", metric), x = "Group", y = metric) +
      theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })

  # ============================================================================
  #  7. P-value Distribution
  # ============================================================================

  # Reactive to assess p-value distribution health
  assess_pvalue_health <- reactive({
    req(values$fit, input$contrast_selector_pvalue)

    # Get p-values
    de_results <- topTable(values$fit, coef = input$contrast_selector_pvalue, number = Inf)
    pvalues <- de_results$P.Value

    n_proteins <- length(pvalues)

    # Bin the p-values into 10 bins
    breaks <- seq(0, 1, by = 0.1)
    hist_counts <- hist(pvalues, breaks = breaks, plot = FALSE)$counts

    # Expected count per bin if uniform
    expected_per_bin <- n_proteins / length(hist_counts)

    # Calculate ratios for different regions
    low_pval_ratio <- sum(pvalues < 0.05) / (n_proteins * 0.05)  # Ratio vs expected 5%
    mid_pval_ratio <- sum(pvalues >= 0.3 & pvalues <= 0.7) / (n_proteins * 0.4)  # Ratio vs expected 40%

    # Detect patterns
    has_spike <- low_pval_ratio > 2  # Spike at zero if > 2x expected
    has_inflation <- mid_pval_ratio > 1.3  # Inflation if mid-range > 1.3x expected
    has_depletion <- low_pval_ratio < 0.5  # Depletion if < 0.5x expected

    # U-shaped: high at both ends
    first_bin_ratio <- hist_counts[1] / expected_per_bin
    last_bin_ratio <- hist_counts[length(hist_counts)] / expected_per_bin
    is_u_shaped <- (first_bin_ratio > 1.5 && last_bin_ratio > 1.5)

    # Completely uniform: no spike, no inflation
    is_uniform <- !has_spike && !has_inflation && (low_pval_ratio > 0.8 && low_pval_ratio < 1.2)

    # Determine overall status
    if (is_u_shaped) {
      status <- "u_shaped"
    } else if (has_inflation) {
      status <- "inflation"
    } else if (has_depletion && !has_spike) {
      status <- "low_power"
    } else if (is_uniform) {
      status <- "uniform"
    } else if (has_spike) {
      status <- "healthy"
    } else {
      status <- "unknown"
    }

    list(
      status = status,
      n_proteins = n_proteins,
      n_significant = sum(de_results$adj.P.Val < 0.05),
      low_pval_ratio = low_pval_ratio,
      mid_pval_ratio = mid_pval_ratio,
      has_spike = has_spike,
      has_inflation = has_inflation,
      has_depletion = has_depletion,
      is_u_shaped = is_u_shaped
    )
  })

  # Render contextual guidance banner
  output$pvalue_guidance <- renderUI({
    health <- assess_pvalue_health()

    if (health$status == "healthy") {
      # Green success banner
      div(class = "alert alert-success", role = "alert",
        icon("check-circle"),
        strong(" P-value distribution looks healthy. "),
        sprintf("Good spike near p=0 (%d significant proteins after FDR correction). ", health$n_significant),
        "This indicates genuine differential expression with proper statistical power."
      )

    } else if (health$status == "inflation") {
      # Yellow warning - p-value inflation
      div(class = "alert alert-warning", role = "alert",
        icon("exclamation-triangle"),
        strong(" Possible p-value inflation detected. "),
        "Too many intermediate p-values (0.3-0.7) relative to expectation. This may indicate:",
        tags$ul(
          tags$li("Unmodeled batch effects → Add batch covariate in Assign Groups tab"),
          tags$li("Variance heterogeneity → Check MDS Plot for outliers"),
          tags$li("Small sample size → Consider adding biological replicates")
        ),
        "If this pattern persists, consider checking the Normalization Diagnostic tab."
      )

    } else if (health$status == "low_power") {
      # Yellow warning - low power
      div(class = "alert alert-warning", role = "alert",
        icon("battery-quarter"),
        strong(" Low statistical power detected. "),
        sprintf("Fewer small p-values than expected (%d significant proteins). ", health$n_significant),
        "Possible causes:",
        tags$ul(
          tags$li("Small sample size → Increase biological replicates if possible"),
          tags$li("High biological variability → Check CV Distribution tab"),
          tags$li("Effect sizes too small to detect with current sample size"),
          tags$li("Over-conservative FDR correction → Consider less stringent threshold")
        )
      )

    } else if (health$status == "u_shaped") {
      # Red danger - U-shaped distribution
      div(class = "alert alert-danger", role = "alert",
        icon("times-circle"),
        strong(" Statistical model issue detected. "),
        "U-shaped p-value distribution (enrichment at both p~0 and p~1) suggests problems with the statistical model or data quality. ",
        "Recommended actions:",
        tags$ul(
          tags$li("Check Normalization Diagnostic - samples may not be properly normalized"),
          tags$li("Review MDS Plot for outlier samples or batch structure"),
          tags$li("Verify group assignments are correct"),
          tags$li("Consider whether this comparison is biologically appropriate")
        )
      )

    } else if (health$status == "uniform") {
      # Blue info - no signal
      div(class = "alert alert-info", role = "alert",
        icon("info-circle"),
        strong(" No differential expression signal detected. "),
        sprintf("P-values are uniformly distributed (only %d proteins pass FDR < 0.05). ", health$n_significant),
        "This could mean:",
        tags$ul(
          tags$li("Groups are truly similar (no biological difference)"),
          tags$li("Test lacks power to detect existing differences"),
          tags$li("Technical variation masks biological signal")
        ),
        "Consider checking QC plots to rule out technical issues."
      )

    } else {
      # Default blue info banner
      div(style = "background-color: #e7f3ff; padding: 12px; border-radius: 5px;",
        icon("info-circle"),
        strong(" P-value Diagnostic: "),
        "This histogram shows the distribution of raw p-values from your differential expression test. ",
        "A healthy analysis shows mostly uniform distribution (flat histogram) with enrichment near p=0 for true positives."
      )
    }
  })

  output$pvalue_histogram <- renderPlot({
    req(values$fit, input$contrast_selector_pvalue)

    # Get all p-values for the current contrast
    de_results <- topTable(values$fit, coef = input$contrast_selector_pvalue, number = Inf)
    pvalues <- de_results$P.Value

    # Calculate expected uniform distribution and bin counts
    n_proteins <- length(pvalues)
    n_bins <- 30
    expected_per_bin <- n_proteins / n_bins
    h <- hist(pvalues, breaks = seq(0, 1, length.out = n_bins + 1), plot = FALSE)
    first_bin_count <- h$counts[1]
    other_max <- max(h$counts[-1])

    # Cap y-axis so the distribution shape is visible; annotate the clipped first bin
    y_max <- max(other_max * 1.5, expected_per_bin * 3)
    hist_data <- data.frame(PValue = pvalues)

    ggplot(hist_data, aes(x = PValue)) +
      geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                     fill = "#4A90E2", color = "white", alpha = 0.7) +
      geom_hline(yintercept = expected_per_bin, linetype = "dashed", color = "red", size = 1) +
      annotate("text", x = 0.75, y = expected_per_bin * 1.15,
               label = "Expected under null (uniform)",
               color = "red", size = 3.5, fontface = "italic") +
      {if (first_bin_count > y_max)
        annotate("text", x = h$mids[1], y = y_max * 0.92,
                 label = paste0("n = ", format(first_bin_count, big.mark = ",")),
                 size = 3.5, fontface = "bold", color = "#2c3e50")
      } +
      coord_cartesian(ylim = c(0, y_max)) +
      labs(
        title = paste0("P-value Distribution (", nrow(de_results), " proteins tested)"),
        subtitle = paste0("Comparison: ", input$contrast_selector_pvalue),
        x = "P-value",
        y = "Number of Proteins"
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray40", size = 11),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  })

  # P-value histogram info modal
  observeEvent(input$pvalue_hist_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " How do I interpret this?"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What this plot shows"),
        p("This histogram displays the distribution of raw (unadjusted) p-values from your differential expression analysis. ",
          "Each bar represents how many proteins have p-values falling in that range."),
        tags$h6("What 'good' looks like"),
        tags$ul(
          tags$li(strong("Flat with a spike at zero: "), "Most p-values uniformly distributed (flat histogram) with a peak near p=0. ",
            "This indicates a mix of non-changing proteins (uniform) and true positives (spike at zero)."),
          tags$li(strong("Expected under the null: "), "For proteins that are truly not changing, p-values should be uniformly distributed between 0 and 1. ",
            "The dashed red line shows this expected uniform distribution.")
        ),
        tags$h6("Warning signs"),
        tags$ul(
          tags$li(strong("Too many intermediate p-values (0.3-0.7): "), "May indicate p-value inflation due to unmodeled variance, batch effects, or outliers."),
          tags$li(strong("Depletion near zero: "), "Too few small p-values suggests the test is overly conservative or lacks statistical power."),
          tags$li(strong("U-shaped distribution: "), "Enrichment at both ends (near 0 and 1) can indicate problems with the statistical model or data quality."),
          tags$li(strong("Completely uniform: "), "No enrichment at p=0 means no differential expression detected, or the test has no power.")
        ),
        tags$h6("What to do if it looks wrong"),
        tags$ul(
          tags$li("Check the Normalization Diagnostic tab to ensure samples are properly normalized"),
          tags$li("Review the MDS plot for outlier samples or unwanted variation"),
          tags$li("Consider adding batch or other covariates to the model if appropriate"),
          tags$li("Verify that sample sizes are adequate for the comparison")
        )
      )
    ))
  })

  # Fullscreen modal for p-value histogram
  observeEvent(input$fullscreen_pvalue_hist, {
    req(values$fit, input$contrast_selector_pvalue)

    # Get all p-values for the current contrast
    de_results <- topTable(values$fit, coef = input$contrast_selector_pvalue, number = Inf)
    pvalues <- de_results$P.Value

    # Calculate expected uniform distribution and bin counts
    n_proteins <- length(pvalues)
    n_bins <- 40  # More bins for fullscreen
    expected_per_bin <- n_proteins / n_bins
    h <- hist(pvalues, breaks = seq(0, 1, length.out = n_bins + 1), plot = FALSE)
    first_bin_count <- h$counts[1]
    other_max <- max(h$counts[-1])

    # Cap y-axis so the distribution shape is visible
    y_max <- max(other_max * 1.5, expected_per_bin * 3)
    hist_data <- data.frame(PValue = pvalues)

    # Create enhanced plot for fullscreen
    p <- ggplot(hist_data, aes(x = PValue)) +
      geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                     fill = "#4A90E2", color = "white", alpha = 0.7) +
      geom_hline(yintercept = expected_per_bin, linetype = "dashed", color = "red", size = 1.2) +
      annotate("text", x = 0.75, y = expected_per_bin * 1.15,
               label = "Expected uniform distribution",
               color = "red", size = 4, fontface = "italic") +
      {if (first_bin_count > y_max)
        annotate("text", x = h$mids[1], y = y_max * 0.92,
                 label = paste0("n = ", format(first_bin_count, big.mark = ",")),
                 size = 4, fontface = "bold", color = "#2c3e50")
      } +
      coord_cartesian(ylim = c(0, y_max)) +
      labs(
        title = paste0("P-value Distribution: ", input$contrast_selector_pvalue),
        subtitle = paste0(nrow(de_results), " proteins tested | ",
                         sum(de_results$adj.P.Val < 0.05), " significant after FDR correction"),
        x = "Raw P-value",
        y = "Number of Proteins"
      ) +
      theme_bw(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(color = "gray40", size = 13),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(breaks = seq(0, 1, 0.1)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

    showModal(modalDialog(
      title = "P-value Distribution - Fullscreen View",
      renderPlot({ p }, height = 700, width = 1000),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # ============================================================================
  #    Section 8: Chromatography QC — TIC traces, diagnostics
  # ============================================================================

  # Flag for conditionalPanel
  output$tic_qc_has_data <- reactive({
    !is.null(values$tic_traces) && length(values$tic_traces) > 0
  })
  outputOptions(output, "tic_qc_has_data", suspendWhenHidden = FALSE)

  # Status badges
  output$tic_qc_status_badges <- renderUI({
    req(values$tic_metrics)
    m <- values$tic_metrics

    n_pass <- sum(m$status == "pass")
    n_warn <- sum(m$status == "warn")
    n_fail <- sum(m$status == "fail")
    n_total <- nrow(m)

    div(style = "display: flex; gap: 8px; margin-bottom: 8px; flex-wrap: wrap;",
      tags$span(class = "badge bg-success", style = "font-size: 0.85em; padding: 5px 10px;",
        sprintf("%d Pass", n_pass)),
      if (n_warn > 0) tags$span(class = "badge bg-warning text-dark",
        style = "font-size: 0.85em; padding: 5px 10px;",
        sprintf("%d Warn", n_warn)),
      if (n_fail > 0) tags$span(class = "badge bg-danger",
        style = "font-size: 0.85em; padding: 5px 10px;",
        sprintf("%d Fail", n_fail)),
      tags$span(class = "badge bg-secondary", style = "font-size: 0.85em; padding: 5px 10px;",
        sprintf("%d Total", n_total))
    )
  })

  # Helper: build TIC plot
  build_tic_plot <- function(view_mode, traces, metrics) {
    status_colors <- c(pass = "#28a745", warn = "#ffc107", fail = "#dc3545")

    if (view_mode == "faceted") {
      # Compute median trace on common grid for overlay
      grid_rt <- seq(
        max(sapply(traces, function(x) min(x$rt_min))),
        min(sapply(traces, function(x) max(x$rt_min))),
        length.out = 300
      )
      interp_matrix <- sapply(traces, function(x) {
        stats::approx(x$rt_min, x$tic, xout = grid_rt, rule = 2)$y
      })
      median_trace <- if (is.matrix(interp_matrix)) {
        apply(interp_matrix, 1, median)
      } else {
        interp_matrix  # single file case
      }

      # Build native plotly subplots (avoids ggplotly facet_wrap first-panel bug)
      run_names <- names(traces)
      ncol <- min(4L, length(run_names))
      nrow <- ceiling(length(run_names) / ncol)

      subplots <- lapply(seq_along(run_names), function(i) {
        nm <- run_names[i]
        df <- traces[[nm]]
        run_status <- metrics$status[metrics$run == nm]
        line_col <- status_colors[run_status] %||% "#666"
        run_label <- sub("\\.d$", "", nm)

        plotly::plot_ly() %>%
          plotly::add_lines(x = df$rt_min, y = df$tic / 1e6,
            line = list(color = line_col, width = 1.5),
            hoverinfo = "text",
            text = sprintf("%s<br>RT: %.1f min<br>TIC: %.1fM", run_label, df$rt_min, df$tic / 1e6),
            showlegend = FALSE) %>%
          plotly::add_lines(x = grid_rt, y = median_trace / 1e6,
            line = list(color = "#007bff", width = 1, dash = "dash"),
            opacity = 0.7, hoverinfo = "skip", showlegend = FALSE) %>%
          plotly::layout(
            annotations = list(list(
              text = run_label, x = 0.5, y = 1.05, xref = "paper", yref = "paper",
              showarrow = FALSE, font = list(size = 10, color = "white"),
              bgcolor = "#2c3e50", borderpad = 3, xanchor = "center"
            )),
            xaxis = list(title = if (i > length(run_names) - ncol) "RT (min)" else "",
                          showgrid = TRUE, gridcolor = "#eee"),
            yaxis = list(title = if ((i - 1) %% ncol == 0) "TIC (M)" else "",
                          showgrid = TRUE, gridcolor = "#eee")
          )
      })

      plotly::subplot(subplots, nrows = nrow, shareX = TRUE, titleX = TRUE, titleY = TRUE) %>%
        plotly::layout(
          title = list(text = "TIC Traces by Run<br><sup style='color:gray'>Blue dashed = median trace across all runs</sup>",
                        font = list(size = 14)),
          showlegend = FALSE
        )

    } else if (view_mode == "overlay") {
      # All runs normalized, overlaid
      plot_data <- do.call(rbind, lapply(names(traces), function(nm) {
        df <- traces[[nm]]
        df$run <- sub("\\.d$", "", nm)
        df$status <- metrics$status[metrics$run == nm]
        df
      }))

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = rt_min, y = tic_norm,
                                                     group = run, color = status,
                                                     text = run)) +
        ggplot2::geom_line(alpha = 0.5, linewidth = 0.4) +
        ggplot2::scale_color_manual(values = status_colors) +
        ggplot2::labs(x = "Retention Time (min)", y = "Normalized TIC",
          title = "TIC Overlay (Normalized 0-1)",
          subtitle = sprintf("%d runs colored by QC status", length(traces))) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(color = "gray50", size = 11),
          legend.position = "bottom"
        )

      plotly::ggplotly(p, tooltip = c("text", "x", "y")) %>%
        plotly::layout(legend = list(orientation = "h", x = 0.3, y = -0.15))

    } else {
      # Metrics view: bar chart of AUC by run
      bar_data <- metrics[metrics$valid, ]
      bar_data$run_label <- sub("\\.d$", "", bar_data$run)
      bar_data$auc_m <- bar_data$total_auc / 1e6
      med_auc <- median(bar_data$auc_m, na.rm = TRUE)

      # Sort by AUC
      bar_data <- bar_data[order(bar_data$auc_m), ]
      bar_data$run_label <- factor(bar_data$run_label, levels = bar_data$run_label)

      p <- ggplot2::ggplot(bar_data, ggplot2::aes(
          x = run_label, y = auc_m, fill = status,
          text = sprintf("Run: %s\nAUC: %.1fM\nPeak RT: %.1f min\nShape r: %.3f\nStatus: %s",
                          run_label, auc_m, peak_rt_min, shape_r, status))) +
        ggplot2::geom_col(width = 0.7) +
        ggplot2::geom_hline(yintercept = med_auc, linetype = "dashed", color = "#666", linewidth = 0.5) +
        ggplot2::scale_fill_manual(values = status_colors) +
        ggplot2::labs(x = NULL, y = "Total AUC (M)",
          title = "TIC Area Under Curve by Run",
          subtitle = sprintf("Dashed line = median AUC (%.1fM)", med_auc)) +
        ggplot2::coord_flip() +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(face = "bold", size = 14),
          plot.subtitle = ggplot2::element_text(color = "gray50", size = 11),
          legend.position = "none"
        )

      plotly::ggplotly(p, tooltip = "text")
    }
  }

  # Main TIC plot
  output$tic_qc_main_plot <- plotly::renderPlotly({
    req(values$tic_traces, values$tic_metrics)
    view_mode <- input$tic_view_mode %||% "faceted"
    build_tic_plot(view_mode, values$tic_traces, values$tic_metrics)
  })

  # Metrics DT table
  output$tic_metrics_table <- DT::renderDT({
    req(values$tic_metrics)
    m <- values$tic_metrics[values$tic_metrics$valid, ]

    size_col <- if ("size_mb" %in% names(m) && !all(is.na(m$size_mb))) {
      data.frame(`Size (MB)` = round(m$size_mb), check.names = FALSE)
    } else {
      NULL
    }

    display <- data.frame(
      Run = sub("\\.d$", "", m$run),
      Status = ifelse(m$status == "pass",
        '<span class="badge bg-success">Pass</span>',
        ifelse(m$status == "warn",
          '<span class="badge bg-warning text-dark">Warn</span>',
          '<span class="badge bg-danger">Fail</span>')),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    if (!is.null(size_col)) display <- cbind(display, size_col)
    display <- cbind(display, data.frame(
      `AUC (M)` = round(m$total_auc / 1e6, 1),
      `Peak RT` = m$peak_rt_min,
      `Gradient Width` = m$gradient_width_min,
      `Baseline Ratio` = sprintf("%.1f%%", m$baseline_ratio * 100),
      `Late Signal` = sprintf("%.1f%%", m$late_signal_ratio * 100),
      `Shape r` = m$shape_r,
      Flags = ifelse(nzchar(m$flags), m$flags, "\u2014"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))

    DT::datatable(display, escape = FALSE, rownames = FALSE,
      options = list(
        pageLength = 50, dom = "t",
        columnDefs = list(list(className = "dt-center", targets = 1))
      )
    )
  })

  # Diagnostics panel
  output$tic_qc_diagnostics <- renderUI({
    req(values$tic_metrics)
    m <- values$tic_metrics
    flagged <- m[m$status != "pass" & nzchar(m$flags), ]

    if (nrow(flagged) == 0) {
      div(class = "alert alert-success", style = "margin-top: 12px;",
        icon("check-circle"), " All runs pass chromatography QC checks.")
    } else {
      alerts <- lapply(seq_len(nrow(flagged)), function(i) {
        row <- flagged[i, ]
        alert_class <- if (row$status == "fail") "alert-danger" else "alert-warning"
        size_info <- if (!is.null(row$size_mb) && !is.na(row$size_mb)) {
          sprintf(" [%d MB]", round(row$size_mb))
        } else ""
        div(class = paste("alert", alert_class), style = "padding: 8px; margin-bottom: 4px; font-size: 0.88em;",
          tags$b(sub("\\.d$", "", row$run)), size_info, ": ",
          row$flags
        )
      })
      div(style = "margin-top: 12px;", alerts)
    }
  })

  # Info modal
  observeEvent(input$tic_qc_info_btn, {
    showModal(modalDialog(
      title = "Chromatography QC",
      size = "l",
      easyClose = TRUE,
      tags$h5("Overview"),
      tags$p("TIC (Total Ion Chromatogram) traces show the summed MS1 signal across the gradient.",
        "Comparing TIC shapes across runs reveals injection failures, loading anomalies, RT drift, and carryover."),
      tags$h5("Views"),
      tags$ul(
        tags$li(tags$b("Faceted:"), " Each run in its own panel with median trace overlay (blue dashed). Color indicates QC status."),
        tags$li(tags$b("Overlay:"), " All runs on one axis, normalized 0-1. Quickly spot outlier shapes."),
        tags$li(tags$b("Metrics:"), " Bar chart of total AUC per run with metrics table below.")
      ),
      tags$h5("Automated Diagnostics"),
      tags$ul(
        tags$li(tags$b("Shape deviation:"), " Pearson r of each trace vs median. r < 0.90 = fail, r < 0.95 = warn."),
        tags$li(tags$b("RT shift:"), " Peak RT > 3 MAD from median."),
        tags$li(tags$b("Loading anomaly:"), " AUC > 3x or < 0.3x median = fail; 2x/0.5x = warn."),
        tags$li(tags$b("Late elution:"), " > 15% signal in last 20% of gradient (carryover)."),
        tags$li(tags$b("Elevated baseline:"), " Baseline > 10% of peak intensity."),
        tags$li(tags$b("Narrow gradient:"), " Effective width < 70% of median.")
      ),
      footer = modalButton("Close")
    ))
  })

  # Fullscreen modal
  observeEvent(input$tic_qc_fullscreen_btn, {
    req(values$tic_traces, values$tic_metrics)
    view_mode <- input$tic_view_mode %||% "faceted"
    p <- build_tic_plot(view_mode, values$tic_traces, values$tic_metrics)

    showModal(modalDialog(
      title = "Chromatography QC \u2014 Fullscreen",
      size = "xl",
      easyClose = TRUE,
      plotly::plotlyOutput("tic_qc_fullscreen_plot", height = "800px"),
      footer = modalButton("Close")
    ))
    output$tic_qc_fullscreen_plot <- plotly::renderPlotly({ p })
  })

  # Extract TIC from QC tab — uses raw file paths from search settings or scanned files
  observeEvent(input$tic_extract_from_qc_btn, {
    # Try to find .d file paths from various sources
    d_files <- NULL

    # Source 1: Already scanned raw files
    if (!is.null(values$diann_raw_files) && nrow(values$diann_raw_files) > 0) {
      df <- values$diann_raw_files
      d_idx <- grepl("\\.d$", df$filename, ignore.case = TRUE)
      if (any(d_idx)) d_files <- df[d_idx, ]
    }

    # Source 2: Search settings output_dir — scan for .d files in parent dir
    if (is.null(d_files)) {
      ss <- values$diann_search_settings
      if (!is.null(ss) && !is.null(ss$output_dir) && nzchar(ss$output_dir %||% "")) {
        # Raw files are usually in the parent of output_dir or a sibling
        raw_dir <- dirname(ss$output_dir)
        cfg <- if (isTRUE(values$ssh_connected) && !is.null(input$ssh_host) && nzchar(input$ssh_host %||% ""))
          list(host = input$ssh_host, user = input$ssh_user,
               port = input$ssh_port %||% 22, key_path = input$ssh_key_path) else NULL

        if (!is.null(cfg)) {
          res <- tryCatch(
            ssh_exec(cfg, sprintf("ls -d %s/*.d 2>/dev/null", shQuote(raw_dir)), timeout = 15),
            error = function(e) list(status = 1, stdout = character()))
          if (res$status == 0 && length(res$stdout) > 0) {
            paths <- trimws(res$stdout)
            paths <- paths[nzchar(paths) & grepl("\\.d$", paths)]
            if (length(paths) > 0) {
              d_files <- data.frame(
                filename = basename(paths),
                full_path = paths,
                size_mb = NA_real_,
                stringsAsFactors = FALSE)
              # Store so the extraction code can find them
              values$diann_raw_files <- d_files
              updateTextInput(session, "ssh_raw_data_dir", value = raw_dir)
            }
          }
        } else if (dir.exists(raw_dir)) {
          d_paths <- list.dirs(raw_dir, recursive = FALSE, full.names = TRUE)
          d_paths <- d_paths[grepl("\\.d$", d_paths)]
          if (length(d_paths) > 0) {
            d_files <- data.frame(
              filename = basename(d_paths),
              full_path = d_paths,
              size_mb = file.size(file.path(d_paths, "analysis.tdf")) / 1e6,
              stringsAsFactors = FALSE)
            values$diann_raw_files <- d_files
          }
        }
      }
    }

    if (is.null(d_files) || nrow(d_files) == 0) {
      showNotification(
        "No .d files found. Go to New Search tab and scan a raw file directory first.",
        type = "warning", duration = 8)
      return()
    }

    # Extract TIC traces
    n_files <- nrow(d_files)
    cfg <- if (isTRUE(values$ssh_connected) && !is.null(input$ssh_host) && nzchar(input$ssh_host %||% ""))
      list(host = input$ssh_host, user = input$ssh_user,
           port = input$ssh_port %||% 22, key_path = input$ssh_key_path) else NULL

    traces <- list()
    withProgress(message = sprintf("Extracting TIC from %d files...", n_files), value = 0, {
      for (i in seq_len(n_files)) {
        fname <- d_files$filename[i]
        incProgress(1 / n_files, detail = fname)
        tic_df <- NULL
        if (!is.null(cfg)) {
          tryCatch({
            raw_dir <- dirname(d_files$full_path[i])
            remote_tdf <- file.path(raw_dir, fname, "analysis.tdf")
            local_tdf <- file.path(tempdir(), paste0("tic_qc_", i, "_analysis.tdf"))
            dl <- scp_download(cfg, remote_tdf, local_tdf)
            if (dl$status == 0 && file.exists(local_tdf)) {
              tic_df <- extract_tic_timstof(local_tdf)
              unlink(local_tdf)
            }
          }, error = function(e) NULL)
        } else {
          tryCatch({
            tdf_path <- file.path(d_files$full_path[i], "analysis.tdf")
            tic_df <- extract_tic_timstof(tdf_path)
          }, error = function(e) NULL)
        }
        if (!is.null(tic_df)) traces[[fname]] <- tic_df
      }
    })

    if (length(traces) == 0) {
      showNotification("TIC extraction failed for all files", type = "error")
      return()
    }

    # Compute metrics (same as server_search.R)
    metrics_list <- lapply(names(traces), function(nm) compute_tic_metrics(traces[[nm]], nm))
    metrics_df <- do.call(rbind, lapply(metrics_list, function(m) {
      data.frame(
        run = m$run, valid = m$valid,
        total_auc = m$total_auc %||% NA_real_,
        peak_rt_min = m$peak_rt_min %||% NA_real_,
        peak_tic = m$peak_tic %||% NA_real_,
        ramp_rt_min = m$ramp_rt_min %||% NA_real_,
        tail_rt_min = m$tail_rt_min %||% NA_real_,
        gradient_width_min = m$gradient_width_min %||% NA_real_,
        baseline_ratio = m$baseline_ratio %||% NA_real_,
        late_signal_ratio = m$late_signal_ratio %||% NA_real_,
        asymmetry = m$asymmetry %||% NA_real_,
        stringsAsFactors = FALSE)
    }))
    metrics_df$size_mb <- d_files$size_mb[match(metrics_df$run, d_files$filename)]
    shape_df <- compute_shape_similarity(traces)
    if (!is.null(shape_df)) {
      metrics_df <- merge(metrics_df, shape_df, by = "run", all.x = TRUE)
    } else {
      metrics_df$shape_r <- 1.0
    }
    diag_results <- lapply(seq_len(nrow(metrics_df)), function(i) {
      diagnose_run(as.list(metrics_df[i, ]), metrics_df, metrics_df$shape_r[i])
    })
    metrics_df$status <- sapply(diag_results, function(d) d$status)
    metrics_df$flags <- sapply(diag_results, function(d) paste(d$flags, collapse = "; "))
    traces <- lapply(traces, normalize_tic)

    values$tic_traces <- traces
    values$tic_metrics <- metrics_df

    n_pass <- sum(metrics_df$status == "pass")
    n_warn <- sum(metrics_df$status == "warn")
    n_fail <- sum(metrics_df$status == "fail")
    showNotification(
      sprintf("TIC extracted: %d pass, %d warn, %d fail", n_pass, n_warn, n_fail),
      type = if (n_fail > 0) "warning" else "message", duration = 6)
  })

}
