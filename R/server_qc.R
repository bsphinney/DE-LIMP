# ==============================================================================
#  SERVER MODULE — QC Trends, QC Diagnostic Plots
#  Called from app.R as: server_qc(input, output, session, values)
# ==============================================================================

server_qc <- function(input, output, session, values) {

  # ============================================================================
  #  1. QC Trend Plots (Precursors, Proteins, MS1 Signal)
  # ============================================================================

  # Helper function to generate QC trend plot (parameterized by metric)
  generate_qc_trend_plot <- function(metric) {
    reactive({
      req(values$qc_stats, values$metadata)
      df <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>%
        mutate(Run_Number = as.numeric(str_extract(Run, "\\d+$")))

      if (input$qc_sort_order == "Group") {
        df <- df %>% arrange(Group, Run_Number)
      } else {
        df <- df %>% arrange(Run_Number)
      }

      df$Sort_Index <- 1:nrow(df)
      df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Group:</b> ", df$Group,
                           "<br><b>", metric, ":</b> ", round(df[[metric]], 2))

      # Calculate group means and ranges for average lines
      group_stats <- df %>%
        group_by(Group) %>%
        summarise(
          mean_value = mean(.data[[metric]], na.rm = TRUE),
          x_min = min(Sort_Index),
          x_max = max(Sort_Index),
          .groups = 'drop'
        )

      # Create plot with bars and group average lines
      p <- ggplot(df, aes(x = Sort_Index, y = .data[[metric]], fill = Group, text = Tooltip)) +
        geom_bar(stat = "identity", width = 0.8) +
        geom_segment(data = group_stats,
                     aes(x = x_min - 0.5, xend = x_max + 0.5,
                         y = mean_value, yend = mean_value,
                         color = Group),
                     linewidth = 1, linetype = "dashed", inherit.aes = FALSE,
                     show.legend = FALSE) +
        scale_color_discrete(guide = "none") +
        theme_minimal() +
        labs(title = paste(metric, "per Run (dashed lines = group averages)"),
             x = "Sample Index (Sorted)", y = metric) +
        theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size=8))

      ggplotly(p, tooltip = "text") %>% config(displayModeBar = TRUE)
    })
  }

  # Three separate plot outputs for the three metrics
  output$qc_trend_plot_precursors <- renderPlotly({
    generate_qc_trend_plot("Precursors")()
  })

  output$qc_trend_plot_proteins <- renderPlotly({
    generate_qc_trend_plot("Proteins")()
  })

  output$qc_trend_plot_ms1 <- renderPlotly({
    generate_qc_trend_plot("MS1_Signal")()
  })

  # Fullscreen modals for each metric
  observeEvent(input$fullscreen_trend_precursors, {
    showModal(modalDialog(
      title = "Precursors Trend - Fullscreen View",
      plotlyOutput("qc_trend_plot_precursors_fs", height = "700px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$qc_trend_plot_precursors_fs <- renderPlotly({
    generate_qc_trend_plot("Precursors")()
  })

  observeEvent(input$fullscreen_trend_proteins, {
    showModal(modalDialog(
      title = "Proteins Trend - Fullscreen View",
      plotlyOutput("qc_trend_plot_proteins_fs", height = "700px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$qc_trend_plot_proteins_fs <- renderPlotly({
    generate_qc_trend_plot("Proteins")()
  })

  observeEvent(input$fullscreen_trend_ms1, {
    showModal(modalDialog(
      title = "MS1 Signal Trend - Fullscreen View",
      plotlyOutput("qc_trend_plot_ms1_fs", height = "700px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$qc_trend_plot_ms1_fs <- renderPlotly({
    generate_qc_trend_plot("MS1_Signal")()
  })

  # ============================================================================
  #  2. QC Stats Table
  # ============================================================================

  output$r_qc_table <- renderDT({ req(values$qc_stats); df_display <- values$qc_stats %>% arrange(Run) %>% mutate(ID = 1:n()) %>% dplyr::select(ID, Run, everything()); datatable(df_display, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) })

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

}
