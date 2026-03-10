# ==============================================================================
#  server_session.R
#  Info modals (contextual help), reproducibility log, session save/load,
#  methodology text, and template import/export.
#
#  NOTE: Fullscreen modals and plot renderers live in their respective modules:
#    - QC trend/diagnostic fullscreens → server_qc.R
#    - Volcano/heatmap/CV fullscreens  → server_de.R
#    - Signal distribution fullscreen  → server_viz.R
#    - P-value diagnostic/histogram    → server_qc.R
# ==============================================================================

server_session <- function(input, output, session, values, add_to_log) {

  # ============================================================================
  #      Load Shared Analysis from URL Query Parameter (?analysis=UUID)
  #      Enables live links from core facility reports on HF Spaces
  # ============================================================================

  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (is.null(query$analysis)) return()

    analysis_id <- query$analysis

    # Validate UUID format
    if (!grepl("^[a-f0-9-]{36}$", analysis_id)) {
      showNotification("Invalid analysis link.", type = "error")
      return()
    }

    # Download .rds from HF dataset repo
    rds_url <- paste0(
      "https://huggingface.co/datasets/brettsp/delimp-shared-analyses/",
      "resolve/main/", analysis_id, ".rds"
    )

    withProgress(message = "Loading shared analysis...", {
      tryCatch({
        temp_rds <- tempfile(fileext = ".rds")
        download.file(rds_url, temp_rds, mode = "wb", quiet = TRUE)
        state <- readRDS(temp_rds)

        # Restore all reactive values (same pattern as session load)
        values$metadata   <- state$metadata
        values$fit        <- state$fit
        values$y_protein  <- state$y_protein
        values$dpc_fit    <- state$dpc_fit
        values$design     <- state$design
        values$qc_stats   <- state$qc_stats
        values$gsea_results <- state$gsea_results
        values$diann_norm_detected <- state$diann_norm_detected %||% "unknown"

        if (!is.null(state$gsea_results_cache)) {
          values$gsea_results_cache <- state$gsea_results_cache
        }
        values$gsea_last_contrast <- state$gsea_last_contrast
        values$gsea_last_org_db   <- state$gsea_last_org_db

        # Phospho
        values$phospho_detected <- state$phospho_detected
        values$phospho_site_matrix <- state$phospho_site_matrix
        values$phospho_site_info <- state$phospho_site_info
        values$phospho_fit <- state$phospho_fit
        values$phospho_site_matrix_filtered <- state$phospho_site_matrix_filtered
        values$phospho_input_mode <- state$phospho_input_mode
        values$ksea_results <- state$ksea_results
        values$ksea_last_contrast <- state$ksea_last_contrast
        values$phospho_annotations <- state$phospho_annotations

        # MOFA (if saved)
        if (!is.null(state$mofa_object)) {
          values$mofa_object <- state$mofa_object
          values$mofa_factors <- state$mofa_factors
          values$mofa_variance_explained <- state$mofa_variance_explained
        }

        # UI state
        values$color_plot_by_de <- state$color_plot_by_de %||% FALSE
        values$cov1_name <- state$cov1_name %||% "Covariate1"
        values$cov2_name <- state$cov2_name %||% "Covariate2"

        # Update contrast selectors
        if (!is.null(values$fit)) {
          cn <- colnames(values$fit$contrasts)
          sel <- state$contrast %||% cn[1]
          updateSelectInput(session, "contrast_selector",
                            choices = cn, selected = sel)
          updateSelectInput(session, "contrast_selector_signal",
                            choices = cn, selected = sel)
          updateSelectInput(session, "contrast_selector_grid",
                            choices = cn, selected = sel)
          updateSelectInput(session, "contrast_selector_pvalue",
                            choices = cn, selected = sel)
        }

        values$status <- "Loaded from shared link"

        title <- if (!is.null(state$metadata_extra))
          state$metadata_extra$title %||% "Shared Analysis"
        else "Shared Analysis"

        showNotification(
          sprintf("Loaded: %s", title),
          type = "message", duration = 8
        )
        nav_select("main_tabs", "DE Dashboard")

      }, error = function(e) {
        showNotification(
          paste("Could not load analysis:", e$message),
          type = "error", duration = 10
        )
      })
    })
  }) |> bindEvent(session$clientData$url_search, once = TRUE)


  # ============================================================================
  #      Info Modals (Contextual Help)
  # ============================================================================

  # --- DE Dashboard Info Modal ---
  observeEvent(input$de_dashboard_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " How to Use the DE Dashboard"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("The Volcano Plot"),
        p("The volcano plot displays all proteins, with log2 fold-change (effect size) on the x-axis and ",
          "-log10 adjusted p-value (statistical significance) on the y-axis. Proteins further from the center ",
          "and higher up are the most interesting."),
        tags$ul(
          tags$li(strong("Click"), " a point to select that protein"),
          tags$li(strong("Box-select"), " (drag) to select multiple proteins"),
          tags$li("Selected proteins are highlighted and the results table filters to show only those proteins")
        ),
        tags$h6("The Results Table"),
        tags$ul(
          tags$li(strong("Gene: "), "Gene symbol (auto-mapped from UniProt accession)"),
          tags$li(strong("Protein Name: "), "Clickable link to UniProt entry"),
          tags$li(strong("logFC: "), "Log2 fold-change between groups. Positive = higher in first group, negative = higher in second group"),
          tags$li(strong("adj.P.Val: "), "FDR-adjusted p-value (Benjamini-Hochberg). Below 0.05 = statistically significant")
        ),
        p("Click rows to select proteins for violin plots, XICs, or heatmap visualization."),
        tags$h6("Threshold Legend"),
        p("The colored lines on the volcano plot show your current significance thresholds: ",
          "blue horizontal line = FDR 0.05, orange vertical lines = fold-change cutoff (adjustable via sidebar slider).")
      )
    ))
  })

  # --- Signal Distribution Info Modal ---
  observeEvent(input$signal_dist_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What is Signal Distribution?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What this shows"),
        p("This boxplot displays the overall protein expression distribution for each sample. ",
          "Each box represents one sample's intensity values across all quantified proteins."),
        tags$h6("DE Coloring"),
        p("When differential expression results are available, samples are colored by the currently selected ",
          "comparison. This helps you visually verify that the groups being compared have distinguishable expression patterns."),
        tags$h6("What to look for"),
        tags$ul(
          tags$li("Boxes should be at roughly the same height \u2014 large shifts indicate normalization issues"),
          tags$li("Samples in the same group should have similar medians"),
          tags$li("Outlier samples (dramatically different from their group) may need investigation")
        )
      )
    ))
  })

  # --- Expression Grid Info Modal ---
  observeEvent(input$expression_grid_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What is the Expression Grid?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What this shows"),
        p("The expression grid is a table of all differentially expressed proteins with their expression values ",
          "across all samples. It provides a complete view of the data behind the statistical results."),
        tags$h6("How to interact"),
        tags$ul(
          tags$li(strong("Click a row"), " to open a violin plot showing that protein's expression across groups"),
          tags$li(strong("Color coding: "), "Expression values are colored from blue (low) through white to red (high)"),
          tags$li(strong("File map: "), "Sample column headers are abbreviated \u2014 the file map above the table shows the full sample names"),
          tags$li(strong("Export: "), "Download the full table as CSV for external analysis")
        ),
        tags$h6("Comparison selector"),
        p("The purple banner controls which comparison's DE results are shown. Changing it updates which proteins ",
          "appear in the grid (only those significant in the selected comparison).")
      )
    ))
  })

  # --- Consistent DE Info Modal ---
  observeEvent(input$consistent_de_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About CV Analysis"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What this tab shows"),
        p("CV Analysis visualizes the relationship between differential expression (logFC) and ",
          "measurement reproducibility (CV) for significant proteins. Lower CV means more reproducible ",
          "quantification \u2014 these are your strongest biomarker candidates."),
        tags$h6("Scatter plot: logFC vs CV"),
        p("Each point is a significant protein. Points in the ", strong("lower-left or lower-right"),
          " corners have both a large effect size and low variability \u2014 the best candidates for ",
          "follow-up validation."),
        tags$h6("How CV is calculated"),
        tags$ol(
          tags$li("Start with all ", strong("significant proteins"), " (FDR-adjusted p-value < 0.05) for the current contrast"),
          tags$li("Convert log2 expression values back to linear scale (2\u207F)"),
          tags$li("Calculate the ", strong("%CV"), " for each protein within each experimental group: ",
                   tags$code("CV = (SD / Mean) \u00d7 100")),
          tags$li("Average the per-group CVs into a single ", strong("Avg CV"), " score")
        ),
        tags$h6("Color coding"),
        tags$ul(
          tags$li(tags$span(style = "background: #d4edda; padding: 2px 6px;", "Green"), " \u2014 CV < 20% (excellent reproducibility)"),
          tags$li(tags$span(style = "background: #fff3cd; padding: 2px 6px;", "Yellow"), " \u2014 CV 20-35% (acceptable)"),
          tags$li(tags$span(style = "background: #f8d7da; padding: 2px 6px;", "Red"), " \u2014 CV > 35% (high variability)")
        ),
        tags$h6("Summary stats"),
        p("The cards below the scatter plot show per-group median CV and the percentage of proteins ",
          "with CV < 20%. Groups with high median CV may have more technical variability or biological heterogeneity."),
        tags$h6("Results Table"),
        p("The Avg CV (%) column also appears in the ", strong("Results Table"), " sub-tab, making it ",
          "easy to sort and filter proteins by reproducibility alongside fold-change and p-value."),
        tags$h6("CV Distribution"),
        p("The histogram below shows the full CV distribution for all significant proteins, ",
          "broken down by group. Tight distributions centered below 20% indicate good data quality.")
      )
    ))
  })

  # --- CV Distribution Info Modal ---
  observeEvent(input$cv_dist_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What is the CV Distribution?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What this shows"),
        p("This histogram displays the distribution of CV values for all significant proteins, ",
          "faceted by experimental group. The dashed vertical line marks the average CV for each group."),
        tags$h6("What to look for"),
        tags$ul(
          tags$li(strong("Peak position: "), "Where most proteins fall. Lower is better \u2014 most proteins should have CV < 30%."),
          tags$li(strong("Long right tail: "), "A few proteins with very high CV. These may be unreliable or biologically variable."),
          tags$li(strong("Group comparison: "), "If one group has consistently higher CVs, it may have more technical variability or biological heterogeneity.")
        ),
        tags$h6("Typical values"),
        p("In well-controlled proteomics experiments, most proteins have CV < 20-30%. ",
          "CVs above 50% suggest the measurement is highly variable for that protein.")
      )
    ))
  })

  # --- Dataset Summary Info Modal ---
  observeEvent(input$dataset_summary_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About Dataset Summary"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("Overview of your uploaded dataset including the number of samples, proteins quantified, ",
          "experimental groups, and data completeness metrics."),
        p("Use this to verify that your data loaded correctly and that the expected number of ",
          "samples and proteins are present before running the analysis pipeline.")
      )
    ))
  })

  # --- Replicate Consistency Info Modal ---
  observeEvent(input$replicate_consistency_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About Replicate Consistency"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Correlation Heatmap"),
        p("Pairwise Pearson correlation of protein intensities across all samples. ",
          "Values near 1.0 indicate high reproducibility. Hierarchical clustering groups ",
          "similar samples together, making it easy to spot outlier replicates that appear ",
          "as rows/columns with lower correlation values."),
        tags$h6("Per-Group Replicate Statistics"),
        tags$ul(
          tags$li(strong("Median CV (%) "), "measures intensity variability across replicates. ",
            "Values below 20% are excellent for proteomics data."),
          tags$li(strong("Mean Correlation "), "is the average pairwise Pearson r within a group. ",
            "Values above 0.95 indicate strong replicate agreement."),
          tags$li(strong("Completeness (%) "), "shows the percentage of proteins quantified in ",
            "every replicate of the group. Higher is better.")
        ),
        p("Export to CSV to include in QC reports or supplementary materials.")
      )
    ))
  })

  # --- Replicate Consistency CSV Export (handler in server_viz.R, co-located with replicate_stats_data) ---

  # --- QC Trends Info Modal ---
  observeEvent(input$qc_trends_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What are QC Trends?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What these plots show"),
        p("The Sample Metrics plot tracks key quality metrics (Precursors, Proteins, MS1 Signal, Data Completeness) ",
          "across all runs in a single faceted view. A LOESS trendline on each facet makes injection drift ",
          "immediately visible."),
        tags$h6("The metrics"),
        tags$ul(
          tags$li(strong("Precursors: "), "Number of peptide precursors identified. A sudden drop may indicate instrument issues or sample problems."),
          tags$li(strong("Proteins: "), "Number of proteins quantified. Should be relatively stable across runs."),
          tags$li(strong("MS1 Signal: "), "Overall MS1 intensity. Gradual decline may indicate column degradation or source contamination."),
          tags$li(strong("Data Completeness (%): "), "Percentage of precursors detected per sample. Low completeness indicates many missed detections.")
        ),
        tags$h6("Sort order"),
        p(strong("Run Order"), " shows samples in the order they were acquired \u2014 useful for spotting instrument drift over time. ",
          strong("Group"), " sorts by experimental condition \u2014 useful for comparing groups side by side."),
        tags$h6("What to look for"),
        tags$ul(
          tags$li("A flat LOESS trendline means stable instrument performance"),
          tags$li("A downward-sloping trendline suggests drift (column degradation, source contamination)"),
          tags$li("Sudden drops or spikes in individual samples flag potential outliers"),
          tags$li("Dashed lines show group averages for quick between-group comparison")
        )
      )
    ))
  })

  # --- Methodology Info Modal ---
  observeEvent(input$methodology_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About the Statistical Methods"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("The LIMPA Pipeline"),
        p("This app uses the LIMPA (LIMma for Proteomics Analysis) package, which adapts the widely-used ",
          "limma framework from genomics for proteomics data."),
        tags$h6("Key steps"),
        tags$ul(
          tags$li(strong("DPC Normalization: "), "Data Point Correspondence handles peptide-to-protein aggregation and missing value estimation"),
          tags$li(strong("Linear model: "), "limma fits a linear model to each protein's expression across groups"),
          tags$li(strong("Empirical Bayes (eBayes): "), "Borrows information across proteins to stabilize variance estimates, especially helpful with small sample sizes"),
          tags$li(strong("FDR correction: "), "Benjamini-Hochberg adjustment controls the false discovery rate at 5%")
        ),
        tags$h6("Why limma?"),
        p("limma has been the gold standard for differential expression analysis for over 20 years. ",
          "Its empirical Bayes approach is particularly powerful for proteomics where sample sizes are often small ",
          "and variance estimates for individual proteins are noisy."),
        tags$h6("Covariates"),
        p("If batch, instrument, or other covariates were specified during group assignment, they are included in the linear model ",
          "to remove their effects before testing for group differences.")
      )
    ))
  })

  # ============================================================================
  #      Settings Modal (Navbar Gear Icon)
  # ============================================================================

  observeEvent(input$open_settings, {
    showModal(modalDialog(
      title = tagList(icon("gear"), " Settings"),
      size = "m",
      passwordInput("settings_api_key", "Gemini API Key",
        value = input$user_api_key %||% "", placeholder = "AIzaSy..."),
      p(class = "text-muted small",
        "Required for AI Analysis. Get a key at ",
        a(href = "https://ai.google.dev", target = "_blank", "ai.google.dev")),
      textInput("settings_model_name", "Model Name",
        value = input$model_name %||% "gemini-3-flash-preview",
        placeholder = "gemini-3-flash-preview"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("save_settings", "Save", class = "btn-success", icon = icon("check"))
      )
    ))
  })

  observeEvent(input$save_settings, {
    updateTextInput(session, "user_api_key", value = input$settings_api_key)
    updateTextInput(session, "model_name", value = input$settings_model_name)
    removeModal()
    showNotification("Settings saved", type = "message", duration = 3)
  })

  # ============================================================================
  #      Reproducibility
  # ============================================================================

  output$reproducible_code <- renderText({
    req(values$repro_log)
    log_content <- paste(values$repro_log, collapse = "\n")
    session_info_text <- paste(capture.output(sessionInfo()), collapse = "\n")

    # Add helpful footer
    footer <- c(
      "",
      "# ==============================================================================",
      "# How to Use This Log:",
      "# 1. Replace 'path/to/your/report.parquet' with your actual file path",
      "# 2. Ensure all required packages are installed (see Session Info below)",
      "# 3. Run sections sequentially from top to bottom",
      "# 4. Each section is timestamped showing when you performed that action",
      "# ==============================================================================",
      "",
      "# --- Session Info (Package Versions) ---",
      session_info_text
    )

    paste(c(log_content, footer), collapse = "\n")
  })

  # Download handler for reproducibility log
  output$download_repro_log <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_reproducibility_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R")
    },
    content = function(file) {
      req(values$repro_log)
      log_content <- paste(values$repro_log, collapse = "\n")
      session_info_text <- paste(capture.output(sessionInfo()), collapse = "\n")

      footer <- c(
        "",
        "# ==============================================================================",
        "# How to Use This Log:",
        "# 1. Replace 'path/to/your/report.parquet' with your actual file path",
        "# 2. Ensure all required packages are installed (see Session Info below)",
        "# 3. Run sections sequentially from top to bottom",
        "# 4. Each section is timestamped showing when you performed that action",
        "# ==============================================================================",
        "",
        "# --- Session Info (Package Versions) ---",
        session_info_text
      )

      writeLines(paste(c(log_content, footer), collapse = "\n"), file)
    }
  )

  # ============================================================================
  #      Save / Load Session (RDS)
  # ============================================================================

  output$save_session <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      session_data <- list(
        raw_data   = values$raw_data,
        metadata   = values$metadata,
        fit        = values$fit,
        y_protein  = values$y_protein,
        dpc_fit    = values$dpc_fit,
        design     = values$design,
        qc_stats   = values$qc_stats,
        gsea_results = values$gsea_results,
        gsea_results_cache = values$gsea_results_cache,
        gsea_last_contrast = values$gsea_last_contrast,
        gsea_last_org_db = values$gsea_last_org_db,
        repro_log  = values$repro_log,
        color_plot_by_de = values$color_plot_by_de,
        # Store UI state so it can be restored
        contrast   = input$contrast_selector,
        logfc_cutoff = input$logfc_cutoff,
        q_cutoff   = input$q_cutoff,
        # Covariate settings
        cov1_name  = values$cov1_name,
        cov2_name  = values$cov2_name,
        # Phosphoproteomics state
        phospho_detected = values$phospho_detected,
        phospho_site_matrix = values$phospho_site_matrix,
        phospho_site_info = values$phospho_site_info,
        phospho_fit = values$phospho_fit,
        phospho_site_matrix_filtered = values$phospho_site_matrix_filtered,
        phospho_input_mode = values$phospho_input_mode,
        # Phospho Phase 2/3
        ksea_results = values$ksea_results,
        ksea_last_contrast = values$ksea_last_contrast,
        phospho_fasta_sequences = values$phospho_fasta_sequences,
        phospho_annotations = values$phospho_annotations,
        # DIA-NN Search state
        diann_jobs = values$diann_jobs,
        diann_fasta_files = values$diann_fasta_files,
        fasta_info = values$fasta_info,
        # MOFA2 Multi-View Integration
        mofa_view_configs = values$mofa_view_configs,
        mofa_views = values$mofa_views,
        mofa_view_fits = values$mofa_view_fits,
        mofa_sample_metadata = values$mofa_sample_metadata,
        mofa_object = values$mofa_object,
        mofa_factors = values$mofa_factors,
        mofa_weights = values$mofa_weights,
        mofa_variance_explained = values$mofa_variance_explained,
        mofa_last_run_params = values$mofa_last_run_params,
        # Instrument metadata
        instrument_metadata = values$instrument_metadata,
        # TIC Chromatography QC
        tic_traces = values$tic_traces,
        tic_metrics = values$tic_metrics,
        # Run Comparator state
        comparator_results          = values$comparator_results,
        comparator_run_a            = values$comparator_run_a,
        comparator_run_b            = values$comparator_run_b,
        comparator_mode             = values$comparator_mode,
        comparator_gemini_narrative = values$comparator_gemini_narrative,
        comparator_mofa             = values$comparator_mofa,
        # Save timestamp & version
        saved_at   = Sys.time(),
        app_version = paste0("DE-LIMP v", values$app_version)
      )
      saveRDS(session_data, file)
      showNotification("Session saved successfully!", type = "message")
    }
  )

  observeEvent(input$load_session_btn, {
    showModal(modalDialog(
      title = "Load Saved Session",
      fileInput("session_file", "Choose .rds session file", accept = ".rds"),
      div(style = "background-color: #fff3cd; padding: 10px; border-radius: 5px;",
        icon("exclamation-triangle"),
        " Loading a session will replace all current data and results."
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_load_session", "Load Session", class = "btn-primary", icon = icon("upload"))
      )
    ))
  })

  observeEvent(input$confirm_load_session, {
    req(input$session_file)
    tryCatch({
      session_data <- readRDS(input$session_file$datapath)

      # Validate that this is a DE-LIMP session file
      required_fields <- c("raw_data", "metadata", "fit")
      if (!all(required_fields %in% names(session_data))) {
        showNotification("Invalid session file: missing required data fields.", type = "error")
        return()
      }

      # Restore reactive values
      values$raw_data   <- session_data$raw_data
      values$metadata   <- session_data$metadata
      values$fit        <- session_data$fit
      values$y_protein  <- session_data$y_protein
      values$dpc_fit    <- session_data$dpc_fit
      values$design     <- session_data$design
      values$qc_stats   <- session_data$qc_stats
      values$gsea_results <- session_data$gsea_results
      if (!is.null(session_data$gsea_results_cache)) values$gsea_results_cache <- session_data$gsea_results_cache
      values$gsea_last_contrast <- session_data$gsea_last_contrast
      values$gsea_last_org_db <- session_data$gsea_last_org_db
      values$color_plot_by_de <- session_data$color_plot_by_de %||% FALSE
      values$cov1_name  <- session_data$cov1_name %||% "Covariate1"
      values$cov2_name  <- session_data$cov2_name %||% "Covariate2"

      # Restore phospho state
      values$phospho_detected <- session_data$phospho_detected
      values$phospho_site_matrix <- session_data$phospho_site_matrix
      values$phospho_site_info <- session_data$phospho_site_info
      values$phospho_fit <- session_data$phospho_fit
      values$phospho_site_matrix_filtered <- session_data$phospho_site_matrix_filtered
      values$phospho_input_mode <- session_data$phospho_input_mode
      # Phospho Phase 2/3
      values$ksea_results <- session_data$ksea_results
      values$ksea_last_contrast <- session_data$ksea_last_contrast
      values$phospho_fasta_sequences <- session_data$phospho_fasta_sequences
      values$phospho_annotations <- session_data$phospho_annotations
      # DIA-NN Search state
      values$diann_jobs <- session_data$diann_jobs %||% list()
      values$diann_fasta_files <- session_data$diann_fasta_files %||% character()
      values$fasta_info <- session_data$fasta_info

      # MOFA2 Multi-View Integration
      if (!is.null(session_data$mofa_object)) {
        values$mofa_view_configs <- session_data$mofa_view_configs %||% list()
        values$mofa_views <- session_data$mofa_views %||% list()
        values$mofa_view_fits <- session_data$mofa_view_fits %||% list()
        values$mofa_sample_metadata <- session_data$mofa_sample_metadata
        values$mofa_object <- session_data$mofa_object
        values$mofa_factors <- session_data$mofa_factors
        values$mofa_weights <- session_data$mofa_weights
        values$mofa_variance_explained <- session_data$mofa_variance_explained
        values$mofa_last_run_params <- session_data$mofa_last_run_params
      }

      # Restore instrument metadata
      if (!is.null(session_data$instrument_metadata)) {
        values$instrument_metadata <- session_data$instrument_metadata
      }

      # Restore TIC Chromatography QC
      if (!is.null(session_data$tic_traces)) {
        values$tic_traces <- session_data$tic_traces
        values$tic_metrics <- session_data$tic_metrics
      }

      # Restore Run Comparator state
      if (!is.null(session_data$comparator_results)) {
        values$comparator_results          <- session_data$comparator_results
        values$comparator_run_a            <- session_data$comparator_run_a
        values$comparator_run_b            <- session_data$comparator_run_b
        values$comparator_mode             <- session_data$comparator_mode
        values$comparator_gemini_narrative <- session_data$comparator_gemini_narrative
        values$comparator_mofa             <- session_data$comparator_mofa
      }

      # Restore repro log and append load event
      values$repro_log  <- session_data$repro_log %||% values$repro_log
      add_to_log("Session Loaded", c(
        sprintf("# Loaded session saved at: %s", session_data$saved_at),
        sprintf("# App version: %s", session_data$app_version %||% "unknown")
      ))

      # Restore UI state: update contrast choices from the fit object
      if (!is.null(values$fit)) {
        contrast_names <- colnames(values$fit$contrasts)
        selected_contrast <- session_data$contrast %||% contrast_names[1]

        # Update all four comparison selectors
        updateSelectInput(session, "contrast_selector", choices = contrast_names, selected = selected_contrast)
        updateSelectInput(session, "contrast_selector_signal", choices = contrast_names, selected = selected_contrast)
        updateSelectInput(session, "contrast_selector_grid", choices = contrast_names, selected = selected_contrast)
        updateSelectInput(session, "contrast_selector_pvalue", choices = contrast_names, selected = selected_contrast)
      }
      # Restore phospho contrast selector
      if (!is.null(values$phospho_fit)) {
        phospho_contrasts <- colnames(values$phospho_fit$contrasts)
        updateSelectInput(session, "phospho_contrast_selector",
                         choices = phospho_contrasts, selected = phospho_contrasts[1])
      }

      if (!is.null(session_data$logfc_cutoff)) {
        updateSliderInput(session, "logfc_cutoff", value = session_data$logfc_cutoff)
      }
      if (!is.null(session_data$q_cutoff)) {
        updateNumericInput(session, "q_cutoff", value = session_data$q_cutoff)
      }

      values$status <- "Session loaded successfully"
      removeModal()
      showNotification(
        paste0("Session loaded! (saved ", format(session_data$saved_at, "%Y-%m-%d %H:%M"), ")"),
        type = "message", duration = 5
      )

      # Record to activity log
      tryCatch({
        record_activity(list(
          event_type = "session_restored",
          timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          user = Sys.getenv("USER", "unknown"),
          search_name = basename(input$session_file$name),
          fasta_files = if (!is.null(values$diann_search_settings))
            paste(basename(values$diann_search_settings$fasta_files), collapse = ", ") else NA,
          fasta_seq_count = if (!is.null(values$diann_search_settings)) values$diann_search_settings$fasta_seq_count else NA,
          n_proteins = if (!is.null(values$y_protein)) nrow(values$y_protein$E) else NA,
          n_samples = if (!is.null(values$y_protein)) ncol(values$y_protein$E) else NA,
          n_contrasts = if (!is.null(values$fit)) length(colnames(values$fit$contrasts)) else NA,
          n_de_proteins = if (!is.null(values$fit)) count_de_proteins(values$fit) else NA,
          app_version = values$app_version %||% "unknown",
          source_type = "session-load",
          notes = "Loaded from session file"
        ))
      }, error = function(e) message("[DE-LIMP] Activity log record failed: ", e$message))

    }, error = function(e) {
      showNotification(paste("Error loading session:", e$message), type = "error")
    })
  })

  # ============================================================================
  #      Methodology Text
  # ============================================================================

  # Helper function to build methodology text (used by render + Claude export)
  build_methodology_text <- function() {
    req(values$fit %||% values$phospho_fit)

    ss <- values$diann_search_settings

    # Instrument & acquisition section (conditional)
    # Also check if stored in search_settings (loaded from job queue)
    instrument_section <- ""
    inst_meta <- values$instrument_metadata
    if (is.null(inst_meta) && !is.null(ss) && is.list(ss))
      inst_meta <- ss$instrument_metadata
    if (!is.null(inst_meta)) {
      inst_text <- format_instrument_methods_text(inst_meta)
      if (nzchar(inst_text)) {
        instrument_section <- paste0("0a. INSTRUMENT & DATA ACQUISITION\n", inst_text, "\n\n")
      }
    }

    # DIA-NN search parameters section (conditional)
    diann_section <- ""
    if (!is.null(ss) && is.list(ss)) {
      sp <- ss$search_params
      tryCatch({
        if (!is.null(sp) && is.list(sp)) {
          # Build enzyme name (null-safe)
          enzyme_name <- if (!is.null(sp$enzyme) && !is.na(sp$enzyme)) {
            switch(sp$enzyme,
              "K*,R*" = "Trypsin/P (cleaves after K, R)",
              "K*,R*,!*P" = "Trypsin (not before P)",
              "R*" = "Arg-C",
              "K*" = "Lys-C",
              sp$enzyme)
          } else "not specified"

          # Build modifications list (null-safe)
          mods <- c()
          if (isTRUE(sp$mod_met_ox)) mods <- c(mods, "Oxidation (M) [UniMod:35]")
          if (isTRUE(sp$mod_nterm_acetyl)) mods <- c(mods, "N-terminal acetylation [UniMod:1]")
          extra_vm <- sp$extra_var_mods
          if (!is.null(extra_vm) && !is.na(extra_vm) && nzchar(extra_vm)) {
            extra <- trimws(strsplit(extra_vm, "\n")[[1]])
            mods <- c(mods, extra[nzchar(extra)])
          }
          if (isTRUE(sp$unimod4)) mods <- c(mods, "Carbamidomethylation (C) [UniMod:4, fixed]")

          # FASTA description
          fasta_desc <- if (length(ss$fasta_files) > 0)
            paste(basename(ss$fasta_files), collapse = ", ") else "not specified"
          contam_desc <- if (!is.null(ss$contaminant_library) && !is.na(ss$contaminant_library) &&
                            ss$contaminant_library != "none")
            sprintf(" with %s contaminant library", gsub("_", " ", ss$contaminant_library))
          else ""

          # Mass accuracy description (null-safe)
          mass_acc_mode <- sp$mass_acc_mode %||% "auto"
          mass_acc_desc <- if (identical(mass_acc_mode, "auto")) {
            "automatic (DIA-NN optimized)"
          } else {
            sprintf("MS2 %s ppm, MS1 %s ppm",
                    sp$mass_acc %||% "?", sp$mass_acc_ms1 %||% "?")
          }

          # DIA-NN engine description (Local vs Docker vs HPC/Singularity vs imported log)
          engine_desc <- if (isTRUE(ss$imported_from_log) && !is.null(ss$diann_version)) {
            sprintf("DIA-NN %s (settings imported from log file)", ss$diann_version)
          } else if (isTRUE(ss$imported_from_log)) {
            "DIA-NN (settings imported from log file)"
          } else if (!is.null(ss$local)) {
            "DIA-NN (local binary)"
          } else if (!is.null(ss$docker_image)) {
            sprintf("DIA-NN (%s, Docker)", ss$docker_image)
          } else if (!is.null(ss$diann_sif)) {
            sprintf("DIA-NN (%s)", ss$diann_sif)
          } else {
            "DIA-NN"
          }

          # Compute resource description
          resource_desc <- if (!is.null(ss$local)) {
            sprintf("Local execution: %d threads.", ss$local$threads %||% 1)
          } else if (!is.null(ss$docker)) {
            sprintf("Docker resources: %d CPUs, %d GB RAM (image: %s).",
                    ss$docker$cpus %||% 1, ss$docker$mem_gb %||% 8, ss$docker$image %||% "?")
          } else if (!is.null(ss$slurm)) {
            sprintf("SLURM resources: %s CPUs, %s GB RAM, %s hour(s), partition: %s.",
                    ss$slurm$cpus %||% "?", ss$slurm$mem_gb %||% "?",
                    ss$slurm$time_hours %||% "?", ss$slurm$partition %||% "?")
          } else if (!is.null(ss$parallel)) {
            sprintf("SLURM parallel 5-step pipeline: %s CPUs/file, %s GB/file, max %s simultaneous.",
                    ss$parallel$cpus_per_file %||% "?", ss$parallel$mem_per_file %||% "?",
                    ss$parallel$max_simultaneous %||% "?")
          } else {
            ""
          }

          # Build search section lines
          search_lines <- c(
            "0b. DIA-NN DATABASE SEARCH",
            sprintf("Raw data files (%s %s files) were searched using %s in %s mode%s.",
                    ss$n_raw_files %||% "?", toupper(ss$raw_file_type %||% "?"), engine_desc,
                    ss$search_mode %||% "library-free",
                    if (!is.null(ss$speclib)) paste0(" with spectral library (", ss$speclib, ")") else ""),
            sprintf("Sequence database: %s%s%s.",
                    fasta_desc, contam_desc,
                    if (!is.null(ss$fasta_seq_count) && !is.na(ss$fasta_seq_count))
                      sprintf(" (%s sequences)", format(ss$fasta_seq_count, big.mark = ","))
                    else ""),
            sprintf("Enzyme: %s with up to %s missed cleavage(s).",
                    enzyme_name, sp$missed_cleavages %||% "?"),
            sprintf("Peptide length: %s-%s amino acids. Precursor m/z: %s-%s. Precursor charge: %s-%s.",
                    sp$min_pep_len %||% "?", sp$max_pep_len %||% "?",
                    sp$min_pr_mz %||% "?", sp$max_pr_mz %||% "?",
                    sp$min_pr_charge %||% "?", sp$max_pr_charge %||% "?"),
            sprintf("Mass accuracy: %s.", mass_acc_desc)
          )

          if (length(mods) > 0) {
            search_lines <- c(search_lines,
              sprintf("Variable modifications (max %s): %s.",
                      sp$max_var_mods %||% "?", paste(mods, collapse = "; ")))
          }

          search_lines <- c(search_lines,
            sprintf("FDR threshold: %s at precursor and protein levels.",
                    sp$qvalue %||% "0.01"),
            sprintf("Match-between-runs: %s. RT-dependent normalization: %s.",
                    if (isTRUE(sp$mbr)) "enabled" else "disabled",
                    if (isTRUE(sp$rt_profiling)) "enabled" else "disabled")
          )

          # Scan window (if manual)
          if (!is.null(sp$scan_window) && !is.na(sp$scan_window) && sp$scan_window > 0) {
            search_lines <- c(search_lines,
              sprintf("Scan window: %s.", sp$scan_window))
          }

          # DIA-NN normalization
          diann_norm <- ss$normalization %||% sp$normalization
          if (!is.null(diann_norm) && !is.na(diann_norm) && nzchar(diann_norm)) {
            norm_label <- switch(diann_norm,
              "on" = "RT-dependent (default)", "off" = "disabled",
              "global" = "global", diann_norm)
            search_lines <- c(search_lines,
              sprintf("DIA-NN normalization: %s.", norm_label))
          }

          # Phospho mode
          if (identical(ss$search_mode, "phospho")) {
            search_lines <- c(search_lines,
              "Phosphoproteomics mode: STY phosphorylation (UniMod:21, +79.966 Da), phospho-specific output and library info reporting enabled.")
          }

          # Extra CLI flags
          extra_flags <- sp$extra_cli_flags
          if (!is.null(extra_flags) && !is.na(extra_flags) && nzchar(extra_flags)) {
            search_lines <- c(search_lines,
              sprintf("Additional flags: %s", extra_flags))
          }

          # Resource description
          if (nzchar(resource_desc)) {
            search_lines <- c(search_lines, resource_desc)
          }

          diann_section <- paste0(paste(search_lines, collapse = "\n"), "\n\n")
        }
      }, error = function(e) {
        message("[DE-LIMP] Methods text DIA-NN section error: ", e$message)
      })
    }

    methodology <- paste(
      "METHODOLOGY\n",
      "===========\n",
      "Data Processing and Statistical Analysis Pipeline\n",
      "---------------------------------------------------\n\n",

      instrument_section,
      diann_section,

      "1. DATA INPUT\n",
      "Raw DIA-NN output files (parquet format) containing precursor-level quantification were\n",
      "imported using the limpa package. Input data includes precursor intensities, protein\n",
      "grouping information, and quality metrics (Q-values).\n\n",

      "2. NORMALIZATION\n",
      "Data Point Correspondence - Cyclic Normalization (DPC-CN) was applied using the dpcCN()\n",
      "function. This method normalizes signal intensities across runs by identifying invariant\n",
      "data points and applying cyclic loess normalization to correct for systematic technical\n",
      "variation while preserving biological differences. DPC-CN is specifically designed for\n",
      "DIA-NN data and performs robust normalization without requiring reference proteins or\n",
      "assuming equal protein abundances across samples.\n\n",

      "3. PROTEIN QUANTIFICATION\n",
      "Normalized precursor-level data were aggregated to protein-level quantification using\n",
      "dpcQuant(). This function employs a modified version of the maxLFQ algorithm, which:\n",
      "  \u2022 Identifies peptides/precursors unique to each protein group\n",
      "  \u2022 Uses pairwise ratios to estimate relative protein abundance\n",
      "  \u2022 Maximizes information from all available peptides while handling missing values\n",
      "  \u2022 Produces log2-transformed protein intensities for downstream analysis\n\n",

      "4. DIFFERENTIAL EXPRESSION ANALYSIS\n",
      "Statistical analysis was performed using the limma framework (Linear Models for\n",
      "Microarray Data), adapted for proteomics data through the dpcDE() function.\n",
      "The analysis workflow includes:\n\n",

      "  a) Linear Model Fitting:\n",
      "     A linear model was fit to the log2-transformed protein intensities with experimental\n",
      "     groups as factors. This model accounts for the mean-variance relationship in the data.\n\n",

      "  b) Empirical Bayes Moderation:\n",
      "     Variance estimates were moderated across proteins using empirical Bayes methods\n",
      "     (eBayes()). This 'borrows information' across proteins to stabilize variance estimates,\n",
      "     particularly beneficial for experiments with limited replicates.\n\n",

      "  c) Contrast Analysis:\n",
      "     Pairwise comparisons between experimental groups were performed using contrast matrices.\n",
      "     Each contrast produces:\n",
      "       \u2022 Log2 fold change (logFC): Effect size of differential expression\n",
      "       \u2022 Average expression (AveExpr): Mean log2 intensity across all samples\n",
      "       \u2022 t-statistic: Test statistic for differential expression\n",
      "       \u2022 P-value: Statistical significance of the change\n",
      "       \u2022 Adjusted P-value (adj.P.Val): FDR-corrected p-value using Benjamini-Hochberg method\n\n",

      "5. MULTIPLE TESTING CORRECTION\n",
      "False Discovery Rate (FDR) control was applied using the Benjamini-Hochberg procedure.\n",
      "This method controls the expected proportion of false positives among rejected hypotheses,\n",
      "providing adjusted p-values (adj.P.Val) that account for testing thousands of proteins\n",
      "simultaneously. Proteins with adj.P.Val < 0.05 are considered statistically significant\n",
      "at 5% FDR.\n\n",

      "6. GENE SET ENRICHMENT ANALYSIS (Optional)\n",
      "When performed, Gene Set Enrichment Analysis (GSEA) was conducted using clusterProfiler.\n",
      "Proteins were ranked by log2 fold-change and mapped from UniProt accessions to Entrez\n",
      "Gene IDs using organism-specific annotation databases (org.Hs.eg.db for human,\n",
      "org.Mm.eg.db for mouse, etc.).\n\n",

      "  Available enrichment databases:\n",
      "    \u2022 GO Biological Process (BP): Biological programs and pathways\n",
      "    \u2022 GO Molecular Function (MF): Molecular activities of gene products\n",
      "    \u2022 GO Cellular Component (CC): Subcellular localization of proteins\n",
      "    \u2022 KEGG Pathways: Metabolic and signaling pathway maps\n\n",

      "  GSEA parameters:\n",
      "    \u2022 Minimum gene set size: 10\n",
      "    \u2022 Maximum gene set size: 500\n",
      "    \u2022 P-value cutoff: 0.05\n",
      "    \u2022 Multiple testing correction: Benjamini-Hochberg FDR\n\n",

      "  The Normalized Enrichment Score (NES) indicates direction and magnitude:\n",
      "    \u2022 Positive NES: Gene set enriched in up-regulated proteins\n",
      "    \u2022 Negative NES: Gene set enriched in down-regulated proteins\n\n\n",

      if (!is.null(values$mofa_object)) paste(
        "7. MULTI-OMICS FACTOR ANALYSIS (MOFA2)\n",
        "Unsupervised multi-view integration was performed using MOFA2 (Multi-Omics Factor\n",
        "Analysis v2). MOFA2 identifies latent factors that capture the major sources of variation\n",
        "across data views, enabling deconvolution of shared vs view-specific biology.\n\n",
        sprintf("  Views: %d\n", values$mofa_last_run_params$n_views %||% 0),
        sprintf("  Common samples: %d\n", values$mofa_last_run_params$n_samples %||% 0),
        sprintf("  Active factors: %d\n", values$mofa_last_run_params$n_factors %||% 0),
        sprintf("  Convergence mode: %s\n", values$mofa_last_run_params$convergence %||% "medium"),
        sprintf("  Scale views: %s\n", values$mofa_last_run_params$scale_views %||% TRUE),
        sprintf("  Random seed: %d\n\n", values$mofa_last_run_params$seed %||% 42),
        "  Reference: Argelaguet R, et al. (2020) MOFA+: a statistical framework for\n",
        "  comprehensive integration of multi-modal single-cell data. Genome Biology 21:111.\n\n\n",
        sep = ""
      ) else "",

      if (!is.null(values$phospho_fit)) paste(
        "8. PHOSPHOSITE-LEVEL DIFFERENTIAL EXPRESSION\n",
        "Unlike standard proteomics where peptides are aggregated to proteins, phosphoproteomics\n",
        "requires site-level analysis because a single protein can harbour dozens of independently\n",
        "regulated phosphorylation sites. DE-LIMP implements site-level DE in a dedicated pipeline\n",
        "that runs in parallel with (and independently of) the protein-level analysis.\n\n",

        "  a) Site Matrix Input:\n",
        "     Two input paths are supported. Path A (recommended): Upload the DIA-NN 1.9+\n",
        "     site-localization matrix (report.phosphosites_*.tsv) which provides protein-relative\n",
        "     residue positions and localization confidence scores. Path B: Parse phosphosites\n",
        "     directly from the DIA-NN report.parquet using UniMod:21 annotations, with peptide-\n",
        "     relative positions.\n\n",

        "  b) Data Filtering:\n",
        "     Sites are filtered for quantification completeness: each site must have at least 2\n",
        "     non-missing values per experimental group. Sites failing this threshold are removed\n",
        "     before imputation. This prevents false discoveries driven by imputation of sites\n",
        "     detected in only one group.\n\n",

        "  c) Missing Value Imputation:\n",
        "     Remaining missing values are imputed using a tail-based (Perseus-style) approach.\n",
        "     Random values are drawn from a downshifted normal distribution:\n",
        "       mean = global_mean - 1.8 * global_sd\n",
        "       sd   = 0.3 * global_sd\n",
        "     This assumes missing values are predominantly 'missing not at random' (MNAR), i.e.,\n",
        "     below the limit of detection. The shift of 1.8 standard deviations places imputed\n",
        "     values in the lower tail of the observed distribution. A fixed random seed (42) is\n",
        "     used for reproducibility.\n\n",

        "     Reference: Tyanova S et al. (2016) The Perseus computational platform for\n",
        "     comprehensive analysis of (prote)omics data. Nature Methods 13:731-740.\n\n",

        "  d) Normalization (Optional):\n",
        "     Three normalization modes are available for the site-level matrix:\n",
        "       - None (default): Use DIA-NN's built-in normalization as-is\n",
        "       - Median centering: Subtract each sample's median, then add back the grand mean.\n",
        "         This aligns sample distributions without changing relative differences.\n",
        "       - Quantile normalization: Forces all sample distributions to be identical using\n",
        "         limma::normalizeBetweenArrays(method='quantile'). This is the most aggressive\n",
        "         option and assumes most sites are unchanged between conditions.\n\n",
        "     CAUTION for phospho-enriched samples: DIA-NN's global normalization assumes most\n",
        "     peptides are unchanged. In phospho-enriched datasets (>30% phosphopeptides), this\n",
        "     assumption may be partially violated. Median centering is recommended if systematic\n",
        "     biases are observed in the QC completeness plot.\n\n",

        "  e) Statistical Testing:\n",
        "     Site-level differential expression is performed using the limma framework:\n",
        "       1. A linear model is fit to each site: lmFit(site_matrix, design)\n",
        "       2. Pairwise contrasts are computed: contrasts.fit(fit, contrasts)\n",
        "       3. Empirical Bayes moderation stabilises variance: eBayes(fit)\n",
        "     This is identical to the protein-level approach but applied to log2-transformed\n",
        "     site intensities rather than aggregated protein quantities. Each site receives:\n",
        "       - logFC: log2 fold change between conditions\n",
        "       - P.Value: raw p-value from the moderated t-test\n",
        "       - adj.P.Val: Benjamini-Hochberg FDR-corrected p-value\n\n",

        "  f) Phosphosite Volcano Plot:\n",
        "     The volcano plot visualises all tested phosphosites with logFC on the x-axis and\n",
        "     -log10(adj.P.Val) on the y-axis. Significance thresholds are |logFC| > 1 (2-fold\n",
        "     change) AND FDR < 0.05. Sites are coloured:\n",
        "       - Red: Significant (both thresholds)\n",
        "       - Blue: FDR < 0.05 only\n",
        "       - Grey: Not significant\n",
        "     The top 15 significant sites are labelled with Gene + Residue + Position.\n",
        "     When protein-level abundance correction is enabled, logFC values represent\n",
        "     phosphorylation stoichiometry changes (site logFC minus protein logFC).\n\n",

        "  g) Residue Distribution:\n",
        "     A bar chart comparing residue composition (Ser/Thr/Tyr) across all quantified\n",
        "     sites versus significant sites. Expected distribution from TiO2 or IMAC enrichment:\n",
        "     ~85% pSer, ~14% pThr, ~1% pTyr. Enrichment of pTyr among significant hits suggests\n",
        "     active tyrosine kinase signalling.\n\n",

        "     Reference: Olsen JV et al. (2006) Global, in vivo, and site-specific\n",
        "     phosphorylation dynamics in signaling networks. Cell 127:635-648.\n\n\n",

        "9. KINASE-SUBSTRATE ENRICHMENT ANALYSIS (KSEA)\n",
        "KSEA infers upstream kinase activity from the phosphosite fold-changes in your data,\n",
        "analogous to gene set enrichment but for kinase-substrate relationships.\n\n",

        "  a) Kinase-Substrate Database:\n",
        "     The KSEAapp R package bundles known kinase-substrate relationships from\n",
        "     PhosphoSitePlus (curated literature evidence) and NetworKIN (computationally\n",
        "     predicted based on motifs and network context). When NetworKIN mode is enabled\n",
        "     (default), predicted substrates with a confidence score >= 5 are included,\n",
        "     substantially increasing coverage.\n\n",

        "  b) Scoring Algorithm:\n",
        "     For each kinase with >= 1 matched substrate in the dataset:\n",
        "       1. Substrate fold-changes are collected from the DE results\n",
        "       2. A mean enrichment score (mS) is computed across matched substrates\n",
        "       3. The z-score normalises the enrichment against the background:\n",
        "            z = (mS * sqrt(m)) / delta\n",
        "          where m = number of substrates and delta = standard deviation of all log2FC\n",
        "       4. P-values are computed from the z-score (normal distribution)\n",
        "       5. FDR correction (Benjamini-Hochberg) is applied across all kinases\n\n",

        "  c) Interpretation:\n",
        "     - Positive z-score (red bars): Substrates collectively up-phosphorylated,\n",
        "       suggesting the kinase is MORE ACTIVE in the treatment condition\n",
        "     - Negative z-score (blue bars): Substrates collectively down-phosphorylated,\n",
        "       suggesting the kinase is LESS ACTIVE in the treatment condition\n",
        "     - The 'n=' label indicates how many substrates contributed to each score;\n",
        "       higher n provides more statistical power and confidence\n\n",

        "     References:\n",
        "     Casado P et al. (2013) Kinase-substrate enrichment analysis provides insights\n",
        "     into the heterogeneity of signaling pathway activation in leukemia cells.\n",
        "     Sci Signaling 6:rs6.\n",
        "     Wiredja DD et al. (2017) The KSEA App: a web-based tool for kinase activity\n",
        "     inference from quantitative phosphoproteomics. Bioinformatics 33:3489-3491.\n\n\n",

        "10. SEQUENCE MOTIF ANALYSIS (Optional)\n",
        "Sequence logos reveal enriched amino acid motifs around regulated phosphosites,\n",
        "which can suggest the kinase families responsible for observed phosphorylation changes.\n\n",

        "  a) Flanking Sequence Extraction:\n",
        "     For each significant phosphosite (FDR < 0.05 and |logFC| > 1), the +/- 7 amino\n",
        "     acid flanking sequence is extracted from the uploaded FASTA protein sequences.\n",
        "     Sites near protein termini are padded with '_' characters. A minimum of 10 sites\n",
        "     is required per direction (up/down) to generate a logo.\n\n",

        "  b) Sequence Logo Visualisation:\n",
        "     Logos are rendered using Shannon information content (bits). Amino acids enriched\n",
        "     above background frequency appear taller. Common motifs:\n",
        "       - Proline at +1: CDK, MAPK, GSK3 (proline-directed kinases)\n",
        "       - Acidic (D/E) at +1 to +3: CK2 family\n",
        "       - Arginine at -3: PKA, PKC, CaMKII (basophilic kinases)\n\n",

        "     Reference: Schwartz D & Gygi SP (2005) An iterative statistical approach to\n",
        "     the identification of protein phosphorylation motifs from large-scale data sets.\n",
        "     Nature Biotechnology 23:1391-1398.\n\n\n",

        "11. PHOSPHOSITE ANNOTATION (Optional)\n",
        "Each phosphosite is annotated as Known (previously reported) or Novel by querying:\n",
        "  - UniProt REST API: Curated 'Modified residue' features with evidence codes\n",
        "    (ECO:0000269 = published experiment, ECO:0007744 = large-scale study)\n",
        "  - PhosphoSitePlus: Kinase-substrate database bundled with KSEAapp\n\n",
        "Matching uses protein accession + residue type + position number. Path A (site matrix)\n",
        "provides protein-relative positions that match database numbering more accurately than\n",
        "Path B (peptide-relative positions).\n\n",

        "  Reference: Hornbeck PV et al. (2015) PhosphoSitePlus, 2014: mutations, PTMs and\n",
        "  recalibrations. Nucleic Acids Res 43:D512-D520.\n\n\n",

        sep = ""
      ) else "",

      "SOFTWARE AND PACKAGES\n",
      "---------------------\n",
      sprintf("Primary analysis: limpa v%s (Bioconductor)\n",
        tryCatch(as.character(packageVersion("limpa")), error = function(e) "unknown")),
      sprintf("Statistical framework: limma v%s (Linear Models for Microarray and RNA-Seq Data)\n",
        tryCatch(as.character(packageVersion("limma")), error = function(e) "unknown")),
      "Data manipulation: dplyr, tidyr\n",
      "Visualization: ggplot2, ComplexHeatmap, plotly\n",
      "Enrichment: clusterProfiler (gseGO, gseKEGG), enrichplot\n",
      "Annotation DBs: org.Hs.eg.db, org.Mm.eg.db (auto-detected)\n",
      if (!is.null(values$mofa_object)) paste(
        "Multi-view integration: MOFA2 (Bioconductor), basilisk (Python env management)\n",
        sep = ""
      ) else "",
      if (!is.null(values$phospho_fit)) paste(
        "Phosphoproteomics: KSEAapp (kinase activity inference), ggseqlogo (motif analysis)\n",
        "Phosphosite databases: PhosphoSitePlus (via KSEAapp), UniProt REST API\n",
        sep = ""
      ) else "",
      sprintf("R version: %s\n", R.version.string),
      sprintf("DE-LIMP version: %s\n\n\n", values$app_version %||% "unknown"),

      "REFERENCES\n",
      "----------\n",
      "\u2022 limpa package: Bioconductor (https://bioconductor.org/packages/limpa/)\n",
      "\u2022 DPC normalization: Designed for DIA-NN proteomics data\n",
      "\u2022 limma: Ritchie ME, et al. (2015) Nucleic Acids Research 43(7):e47\n",
      "\u2022 Empirical Bayes: Smyth GK (2004) Statistical Applications in Genetics and\n",
      "  Molecular Biology 3:Article3\n",
      "\u2022 FDR control: Benjamini Y, Hochberg Y (1995) Journal of the Royal Statistical\n",
      "  Society 57(1):289-300\n",
      "\u2022 clusterProfiler: Yu G, et al. (2012) OMICS 16(5):284-287\n",
      "\u2022 GSEA method: Subramanian A, et al. (2005) PNAS 102(43):15545-15550\n",
      "\u2022 KEGG database: Kanehisa M, Goto S (2000) Nucleic Acids Research 28(1):27-30\n",
      if (!is.null(values$mofa_object)) paste(
        "\u2022 MOFA2: Argelaguet R, et al. (2020) Genome Biology 21:111\n",
        sep = ""
      ) else "",
      if (!is.null(values$phospho_fit)) paste(
        "\u2022 Perseus imputation: Tyanova S, et al. (2016) Nature Methods 13:731-740\n",
        "\u2022 Phosphoproteome dynamics: Olsen JV, et al. (2006) Cell 127:635-648\n",
        "\u2022 KSEA: Casado P, et al. (2013) Sci Signaling 6:rs6\n",
        "\u2022 KSEAapp: Wiredja DD, et al. (2017) Bioinformatics 33:3489-3491\n",
        "\u2022 PhosphoSitePlus: Hornbeck PV, et al. (2015) Nucleic Acids Res 43:D512-D520\n",
        "\u2022 Phospho motifs: Schwartz D & Gygi SP (2005) Nature Biotechnology 23:1391-1398\n",
        "\u2022 ggseqlogo: Wagih O (2017) Bioinformatics 33:3645-3647\n",
        sep = ""
      ) else "",
      "\n\n",

      "CITATION\n",
      "--------\n",
      "If you use this analysis in your research, please cite:\n",
      "\u2022 The limpa package (Bioconductor)\n",
      "\u2022 The limma package: Ritchie ME, et al. (2015) Nucleic Acids Research\n",
      "\u2022 DIA-NN: Demichev V, et al. (2020) Nature Methods 17:41-44",

      sep = ""
    )

    methodology
  }

  output$methodology_text <- renderText({ build_methodology_text() })

  # Store methodology text in values so other modules (e.g., Claude export) can access it
  observe({
    tryCatch({
      values$methodology_text <- build_methodology_text()
    }, error = function(e) NULL)
  })

  # ============================================================================
  #      Template Import / Export
  # ============================================================================

  # Export group assignment template
  output$export_template <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_group_template_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      # Get current table data (including any edits)
      template_data <- if (!is.null(input$hot_metadata)) {
        hot_to_r(input$hot_metadata)
      } else {
        values$metadata
      }

      # Export with current custom covariate names
      write.csv(template_data, file, row.names = FALSE)
      showNotification("Template exported successfully!", type = "message", duration = 3)
    }
  )

  # Import group assignment template
  observeEvent(input$import_template, {
    showModal(modalDialog(
      title = "Import Group Assignment Template",
      fileInput("template_file", "Choose CSV File",
                accept = c("text/csv", "text/comma-separated-values", ".csv")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_import", "Import", class = "btn-primary")
      )
    ))
  })

  observeEvent(input$confirm_import, {
    req(input$template_file)

    tryCatch({
      imported_data <- read.csv(input$template_file$datapath, stringsAsFactors = FALSE)

      # Validate columns
      required_cols <- c("ID", "File.Name", "Group", "Batch", "Covariate1", "Covariate2")
      if (!all(required_cols %in% colnames(imported_data))) {
        showNotification("Error: Template must have columns: ID, File.Name, Group, Batch, Covariate1, Covariate2",
                        type = "error", duration = 10)
        return()
      }

      # Validate File.Name matches
      if (!all(imported_data$File.Name %in% values$metadata$File.Name)) {
        showNotification("Warning: Some file names in template don't match current data. Using matching rows only.",
                        type = "warning", duration = 8)
      }

      # Update metadata with imported data (match by File.Name)
      values$metadata <- imported_data

      showNotification("Template imported successfully!", type = "message", duration = 3)
      removeModal()  # Close import dialog

      # Navigate to Assign Groups sub-tab to show imported data
      nav_select("main_tabs", "Data Overview")
      nav_select("data_overview_tabs", "Assign Groups & Run")

    }, error = function(e) {
      showNotification(paste("Error importing template:", e$message), type = "error", duration = 10)
    })
  })

  # ============================================================================
  #      About Tab — Version Info, Community Stats, Trend Sparklines
  # ============================================================================

  output$about_version_text <- renderText({
    paste0("v", values$app_version)
  })

  output$community_stats_cards <- renderUI({
    stats <- values$community_stats
    if (is.null(stats)) {
      return(div(style = "text-align: center; padding: 20px; color: #a0aec0;",
        icon("chart-bar"), " Community stats not yet available \u2014 they'll appear after the first GitHub Actions run."
      ))
    }

    gh <- stats$github
    make_card <- function(value, label, icon_name, color) {
      div(style = paste0(
        "flex: 1; min-width: 150px; text-align: center; padding: 15px; ",
        "background: white; border-radius: 8px; border: 1px solid #e2e8f0; ",
        "box-shadow: 0 1px 3px rgba(0,0,0,0.05);"
      ),
        icon(icon_name, style = paste0("font-size: 1.5em; color: ", color, "; margin-bottom: 8px;")),
        tags$div(style = "font-size: 1.8em; font-weight: 600;", format(value, big.mark = ",")),
        tags$div(style = "color: #718096; font-size: 0.85em;", label)
      )
    }

    div(style = "display: flex; gap: 15px; flex-wrap: wrap; justify-content: center; margin-bottom: 20px;",
      make_card(gh$stars %||% 0, "Stars", "star", "#f6ad55"),
      make_card(gh$forks %||% 0, "Forks", "code-branch", "#68d391"),
      make_card(gh$unique_visitors_14d %||% 0, "Visitors (14d)", "eye", "#63b3ed"),
      make_card(gh$unique_users_14d %||% 0, "Clones (14d)", "download", "#b794f4")
    )
  })

  output$community_trend_plots <- renderUI({
    stats <- values$community_stats
    if (is.null(stats) || is.null(stats$trends)) return(NULL)

    has_clones <- !is.null(stats$trends$clones) && length(stats$trends$clones$date) > 1
    has_views <- !is.null(stats$trends$views) && length(stats$trends$views$date) > 1

    if (!has_clones && !has_views) return(NULL)

    div(style = "display: flex; gap: 15px; flex-wrap: wrap; justify-content: center; margin-bottom: 15px;",
      if (has_views) div(style = "flex: 1; min-width: 300px;",
        tags$p(style = "text-align: center; font-size: 0.85em; color: #718096; margin-bottom: 4px;",
          "Unique Visitors (14-day)"),
        plotlyOutput("about_views_sparkline", height = "120px")
      ),
      if (has_clones) div(style = "flex: 1; min-width: 300px;",
        tags$p(style = "text-align: center; font-size: 0.85em; color: #718096; margin-bottom: 4px;",
          "Unique Clones (14-day)"),
        plotlyOutput("about_clones_sparkline", height = "120px")
      )
    )
  })

  output$about_views_sparkline <- renderPlotly({
    stats <- values$community_stats
    req(stats, stats$trends, stats$trends$views)
    trend <- stats$trends$views
    req(length(trend$date) > 1)

    plot_ly(x = as.Date(trend$date), y = trend$count, type = "scatter", mode = "lines+markers",
            line = list(color = "#63b3ed", width = 2),
            marker = list(color = "#63b3ed", size = 5),
            hovertemplate = "%{x|%b %d}: %{y} visitors<extra></extra>") %>%
      layout(
        xaxis = list(title = "", showgrid = FALSE, tickformat = "%b %d"),
        yaxis = list(title = "", showgrid = TRUE, gridcolor = "#f0f0f0", rangemode = "tozero"),
        margin = list(l = 30, r = 10, t = 5, b = 30),
        plot_bgcolor = "white", paper_bgcolor = "white"
      ) %>%
      config(displayModeBar = FALSE)
  })

  output$about_clones_sparkline <- renderPlotly({
    stats <- values$community_stats
    req(stats, stats$trends, stats$trends$clones)
    trend <- stats$trends$clones
    req(length(trend$date) > 1)

    plot_ly(x = as.Date(trend$date), y = trend$count, type = "scatter", mode = "lines+markers",
            line = list(color = "#b794f4", width = 2),
            marker = list(color = "#b794f4", size = 5),
            hovertemplate = "%{x|%b %d}: %{y} clones<extra></extra>") %>%
      layout(
        xaxis = list(title = "", showgrid = FALSE, tickformat = "%b %d"),
        yaxis = list(title = "", showgrid = TRUE, gridcolor = "#f0f0f0", rangemode = "tozero"),
        margin = list(l = 30, r = 10, t = 5, b = 30),
        plot_bgcolor = "white", paper_bgcolor = "white"
      ) %>%
      config(displayModeBar = FALSE)
  })

  output$community_discussions <- renderUI({
    stats <- values$community_stats
    if (is.null(stats) || is.null(stats$discussions)) return(NULL)

    disc <- stats$discussions
    total <- disc$total %||% 0
    recent <- disc$recent

    if (is.null(recent) || length(recent$title) == 0) return(NULL)

    # Category icon mapping
    cat_icons <- list(
      Announcements = "bullhorn", General = "comments",
      Ideas = "lightbulb", "Q&A" = "circle-question",
      Polls = "square-poll-vertical", "Show and tell" = "hands-clapping"
    )

    discussion_items <- lapply(seq_along(recent$title), function(i) {
      cat_name <- recent$category[i]
      ic <- cat_icons[[cat_name]] %||% "comment"
      comment_count <- recent$comments[i] %||% 0
      comment_badge <- if (comment_count > 0) {
        tags$span(class = "badge bg-secondary ms-1", style = "font-size: 0.75em;",
          icon("comment", style = "font-size: 0.8em;"), " ", comment_count)
      }

      div(style = "padding: 8px 12px; border-bottom: 1px solid #f0f0f0;",
        div(style = "display: flex; align-items: center; gap: 8px;",
          icon(ic, style = "color: #718096; width: 16px;"),
          tags$a(href = recent$url[i], target = "_blank",
            style = "color: #2d3748; text-decoration: none; font-size: 0.9em; font-weight: 500;",
            recent$title[i]),
          comment_badge
        ),
        div(style = "margin-left: 24px; color: #a0aec0; font-size: 0.78em;",
          cat_name, " \u00b7 ", recent$author[i], " \u00b7 ", recent$created[i]
        )
      )
    })

    div(style = "margin-top: 25px; max-width: 700px; margin-left: auto; margin-right: auto;",
      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;",
        tags$h6(style = "margin: 0;", icon("comments"), " Recent Discussions"),
        tags$a(href = "https://github.com/bsphinney/DE-LIMP/discussions", target = "_blank",
          style = "font-size: 0.8em; color: #63b3ed;",
          sprintf("View all (%d)", total))
      ),
      div(style = "background: white; border: 1px solid #e2e8f0; border-radius: 8px; overflow: hidden;",
        discussion_items
      )
    )
  })

  output$stats_updated_at <- renderUI({
    stats <- values$community_stats
    if (is.null(stats) || is.null(stats$updated_at)) return(NULL)

    updated <- tryCatch(
      format(as.POSIXct(stats$updated_at, format = "%Y-%m-%dT%H:%M:%S"), "%B %d, %Y at %H:%M UTC"),
      error = function(e) stats$updated_at
    )
    div(style = "text-align: center; margin-top: 15px; color: #a0aec0; font-size: 0.8em;",
      paste("Last updated:", updated)
    )
  })

  # ==========================================================================
  #   Unified History — activity log, projects, notes
  # ==========================================================================

  history_refresh <- reactiveVal(0)

  # Migrate old CSVs + backfill from job queue on startup
  session$onFlushed(once = TRUE, function() {
    tryCatch({
      migrated <- migrate_to_activity_log()
      if (isTRUE(migrated)) message("[DE-LIMP] Migration completed")
      jobs <- isolate(values$diann_jobs)
      if (length(jobs) > 0) {
        n <- backfill_activity_log(jobs, app_version = isolate(values$app_version) %||% "unknown")
        if (n > 0) history_refresh(isolate(history_refresh()) + 1)
      }
    }, error = function(e) message("[DE-LIMP] History startup error: ", e$message))
  })

  # Populate project filter dropdown
  observe({
    history_refresh()
    choices <- get_projects()
    updateSelectizeInput(session, "project_filter",
      choices = c("All projects" = "", choices),
      selected = input$project_filter %||% "",
      server = FALSE)
  })

  # Populate user filter dropdown
  observe({
    history_refresh()
    log <- activity_log_read()
    users <- if (nrow(log) > 0 && "user" %in% names(log)) {
      sort(unique(na.omit(log$user[nzchar(log$user)])))
    } else character(0)
    updateSelectizeInput(session, "history_user_filter",
      choices = c("All users" = "", users),
      selected = input$history_user_filter %||% "",
      server = FALSE)
  })

  # Project summary cards (shown when a project is selected)
  output$project_summary_cards <- renderUI({
    proj_name <- input$project_filter
    if (is.null(proj_name) || !nzchar(proj_name)) return(NULL)

    history_refresh()
    log <- activity_log_read()
    if (nrow(log) == 0 || !"project" %in% names(log)) return(NULL)

    proj_log <- log[!is.na(log$project) & log$project == proj_name, , drop = FALSE]
    if (nrow(proj_log) == 0) return(NULL)

    n_entries <- nrow(proj_log)
    total_samples <- sum(as.numeric(proj_log$n_samples), na.rm = TRUE)
    total_de <- sum(as.numeric(proj_log$n_de_proteins), na.rm = TRUE)
    dates <- sort(proj_log$timestamp)
    date_range <- if (length(dates) >= 2) {
      paste(substr(dates[1], 1, 10), "\u2013", substr(dates[length(dates)], 1, 10))
    } else if (length(dates) == 1) substr(dates[1], 1, 10) else ""

    div(style = "margin-bottom: 15px;",
      div(style = "display: flex; gap: 10px; flex-wrap: wrap; margin-bottom: 8px;",
        div(style = "flex: 1; min-width: 120px; background: #e8f5e9; padding: 10px; border-radius: 6px; text-align: center;",
          tags$strong(n_entries), tags$br(), tags$small("Entries")),
        div(style = "flex: 1; min-width: 120px; background: #e3f2fd; padding: 10px; border-radius: 6px; text-align: center;",
          tags$strong(total_samples), tags$br(), tags$small("Total Samples")),
        div(style = "flex: 1; min-width: 120px; background: #fff3e0; padding: 10px; border-radius: 6px; text-align: center;",
          tags$strong(total_de), tags$br(), tags$small("Total DE Proteins")),
        div(style = "flex: 1; min-width: 120px; background: #f3e5f5; padding: 10px; border-radius: 6px; text-align: center;",
          tags$strong(date_range), tags$br(), tags$small("Date Range"))
      )
    )
  })

  # ==========================================================================
  #   Unified History Table
  # ==========================================================================

  observeEvent(input$history_refresh_btn, {
    history_refresh(isolate(history_refresh()) + 1)
  }, ignoreInit = TRUE)

  output$history_table <- renderDT({
    history_refresh()
    message("[DE-LIMP] Rendering unified history table...")
    log <- activity_log_read()
    if (nrow(log) == 0) return(NULL)

    log <- log[order(log$timestamp, decreasing = TRUE), ]

    # Join live status from job queue
    tryCatch({
      jobs <- values$diann_jobs
      if (length(jobs) > 0) {
        job_map <- list()
        for (j in jobs) {
          if (!is.null(j$output_dir) && nzchar(j$output_dir))
            job_map[[j$output_dir]] <- j$status
        }
        for (i in seq_len(nrow(log))) {
          od <- log$output_dir[i]
          if (!is.na(od) && nzchar(od) && !is.null(job_map[[od]])) {
            log$status[i] <- job_map[[od]]
          }
        }
      }
    }, error = function(e) NULL)

    # Apply filters
    proj_filter <- input$project_filter
    if (!is.null(proj_filter) && nzchar(proj_filter)) {
      log <- log[!is.na(log$project) & log$project == proj_filter, , drop = FALSE]
      if (nrow(log) == 0) return(NULL)
    }
    status_filter <- input$history_status_filter
    if (!is.null(status_filter) && nzchar(status_filter)) {
      log <- log[!is.na(log$status) & log$status == status_filter, , drop = FALSE]
      if (nrow(log) == 0) return(NULL)
    }
    user_filter <- input$history_user_filter
    if (!is.null(user_filter) && nzchar(user_filter)) {
      log <- log[!is.na(log$user) & log$user == user_filter, , drop = FALSE]
      if (nrow(log) == 0) return(NULL)
    }

    # Event type badge
    log$type_badge <- vapply(log$event_type, function(et) {
      switch(et %||% "",
        "search_submitted" = '<span class="badge bg-info">Search</span>',
        "search_completed" = '<span class="badge bg-info">Search</span>',
        "search_failed"    = '<span class="badge bg-info">Search</span>',
        "analysis_completed" = '<span class="badge bg-success">Analysis</span>',
        "data_loaded"      = '<span class="badge bg-secondary">Load</span>',
        "session_restored" = '<span class="badge bg-warning">Restore</span>',
        '<span class="badge bg-light text-dark">Other</span>')
    }, character(1))

    # Status badge
    log$status_badge <- vapply(log$status %||% "", function(s) {
      switch(s %||% "",
        "completed" = '<span class="badge bg-success">Completed</span>',
        "failed"    = '<span class="badge bg-danger">Failed</span>',
        "running"   = '<span class="badge bg-primary">Running</span>',
        "submitted" = '<span class="badge bg-info">Submitted</span>',
        "cancelled" = '<span class="badge bg-warning">Cancelled</span>',
        "")
    }, character(1))

    # Duration
    log$duration_fmt <- vapply(log$duration_min, function(d) {
      if (is.na(d)) return("")
      d <- as.numeric(d)
      if (d < 60) sprintf("%.0f min", d) else sprintf("%.1f hr", d / 60)
    }, character(1))

    # Truncated name
    trunc <- function(x, n = 30) {
      x <- as.character(x %||% "")
      ifelse(nchar(x) > n, paste0(substr(x, 1, n), "\u2026"), x)
    }
    log$name_short <- trunc(log$search_name, 30)

    # Action buttons
    btn_style <- "font-size:0.75em;padding:2px 6px;margin:1px;"
    log$actions <- vapply(seq_len(nrow(log)), function(i) {
      od <- log$output_dir[i]
      id <- log$id[i]
      sf <- if ("session_file" %in% names(log)) log$session_file[i] else NA
      btns <- ""
      if (!is.na(od) && nzchar(od)) {
        od_esc <- gsub("'", "\\\\'", od)
        id_esc <- gsub("'", "\\\\'", id %||% "")

        # Log button (view search_info.md)
        btns <- sprintf(
          '<button class="btn btn-outline-info btn-xs" style="%s" onclick="event.stopPropagation();Shiny.setInputValue(\'history_log_click\', {od: \'%s\', ts: Date.now()})"><i class="fa fa-circle-info"></i></button>',
          btn_style, od_esc)

        # Load button — green for full session (RDS), outline for raw data (parquet)
        has_session <- !is.na(sf) && nzchar(sf %||% "")
        is_analysis <- !is.na(log$event_type[i]) && log$event_type[i] == "analysis_completed"
        if (has_session || is_analysis) {
          # Full session restore (RDS exists)
          sf_target <- if (has_session) sf else file.path(od, "session.rds")
          sf_esc <- gsub("'", "\\\\'", sf_target)
          btns <- paste0(btns, sprintf(
            '<button class="btn btn-success btn-xs" style="%s" onclick="event.stopPropagation();Shiny.setInputValue(\'history_load_click\', {sf: \'%s\', od: \'%s\', ts: Date.now()})"><i class="fa fa-download"></i> Load</button>',
            btn_style, sf_esc, od_esc))
        } else if (!is.na(log$status[i]) && log$status[i] %in% c("completed", "")) {
          # Search-only entry — load raw report.parquet (no pipeline results)
          btns <- paste0(btns, sprintf(
            '<button class="btn btn-outline-success btn-xs" style="%s" onclick="event.stopPropagation();Shiny.setInputValue(\'history_load_click\', {sf: \'\', od: \'%s\', ts: Date.now()})"><i class="fa fa-file-arrow-down"></i> Raw</button>',
            btn_style, od_esc))
        }

        # Settings button (for search events)
        if (!is.na(log$event_type[i]) && grepl("search", log$event_type[i])) {
          btns <- paste0(btns, sprintf(
            '<button class="btn btn-outline-warning btn-xs" style="%s" onclick="event.stopPropagation();Shiny.setInputValue(\'history_settings_click\', {od: \'%s\', ts: Date.now()})"><i class="fa fa-file-import"></i></button>',
            btn_style, od_esc))
        }

        # Notes button
        btns <- paste0(btns, sprintf(
          '<button class="btn btn-outline-secondary btn-xs" style="%s" onclick="event.stopPropagation();Shiny.setInputValue(\'history_notes_click\', {id: \'%s\', od: \'%s\', notes: \'%s\', name: \'%s\', ts: Date.now()})"><i class="fa fa-pen"></i></button>',
          btn_style, id_esc, od_esc,
          gsub("'", "\\\\'", log$notes[i] %||% ""),
          gsub("'", "\\\\'", log$search_name[i] %||% "")))

        # Project button
        btns <- paste0(btns, sprintf(
          '<button class="btn btn-outline-secondary btn-xs" style="%s" onclick="event.stopPropagation();Shiny.setInputValue(\'history_project_click\', {od: \'%s\', ts: Date.now()})"><i class="fa fa-folder-open"></i></button>',
          btn_style, od_esc))
      }
      btns
    }, character(1))

    # Build hidden details for child row expansion
    log$details <- vapply(seq_len(nrow(log)), function(i) {
      items <- c()
      add <- function(label, val) {
        if (!is.null(val) && !is.na(val) && nzchar(as.character(val)))
          items <<- c(items, paste0("<b>", label, ":</b> ", htmltools::htmlEscape(as.character(val))))
      }
      add("FASTA", log$fasta_files[i])
      add("Enzyme", log$enzyme[i])
      add("Mass Acc", log$mass_acc_mode[i])
      if (!is.na(log$mass_acc[i])) add("MS2 Acc", log$mass_acc[i])
      if (!is.na(log$mass_acc_ms1[i])) add("MS1 Acc", log$mass_acc_ms1[i])
      if (!is.na(log$scan_window[i]) && log$scan_window[i] != 0) add("Scan Window", log$scan_window[i])
      add("Normalization", log$normalization[i])
      if (!is.na(log$extra_cli_flags[i]) && nzchar(log$extra_cli_flags[i])) add("Extra Flags", log$extra_cli_flags[i])
      add("Job ID", log$job_id[i])
      if (!is.na(log$output_dir[i]))
        items <- c(items, paste0("<b>Output:</b> <code style='font-size:0.85em;'>",
          htmltools::htmlEscape(log$output_dir[i]), "</code>"))
      if (!is.na(log$n_proteins[i])) add("Proteins", log$n_proteins[i])
      if (!is.na(log$n_de_proteins[i])) add("DE Proteins", log$n_de_proteins[i])
      if (!is.na(log$notes[i]) && nzchar(log$notes[i]))
        items <- c(items, paste0("<b>Notes:</b> <em>", htmltools::htmlEscape(log$notes[i]), "</em>"))
      if (length(items) == 0) return("")
      paste0('<div style="padding:8px 12px;background:#f8f9fa;color:#2d3748;font-size:0.9em;line-height:1.8;">',
        paste(items, collapse = "<br>"), '</div>')
    }, character(1))

    # Compare checkboxes — only for entries with a session.rds (post-pipeline)
    log$compare_cb <- vapply(seq_len(nrow(log)), function(i) {
      has_od <- !is.na(log$output_dir[i]) && nzchar(log$output_dir[i])
      has_sf <- "session_file" %in% names(log) && !is.na(log$session_file[i]) && nzchar(log$session_file[i] %||% "")
      is_analysis <- !is.na(log$event_type[i]) && log$event_type[i] == "analysis_completed"
      # Show checkbox if: explicit session_file exists, OR analysis_completed with output_dir
      # (new auto-save puts session.rds in output_dir)
      if (has_od && (has_sf || is_analysis)) {
        od_esc <- gsub("'", "\\\\'", log$output_dir[i])
        nm_esc <- gsub("'", "\\\\'", log$search_name[i] %||% "")
        sf_val <- if (has_sf) log$session_file[i] else file.path(log$output_dir[i], "session.rds")
        sf_esc <- gsub("'", "\\\\'", sf_val)
        sprintf('<input type="checkbox" class="history-compare-cb" data-od="%s" data-name="%s" data-sf="%s" onclick="event.stopPropagation();window._delimp_updateCompare();">', od_esc, nm_esc, sf_esc)
      } else ""
    }, character(1))

    display_cols <- c("details", "timestamp", "user", "name_short", "type_badge", "n_files",
      "n_proteins", "n_de_proteins", "status_badge", "duration_fmt", "project", "actions", "compare_cb")
    display_cols <- intersect(display_cols, names(log))

    no_sort_cols <- which(display_cols %in% c("details", "type_badge", "status_badge", "actions", "compare_cb")) - 1

    col_names <- c(details = "", timestamp = "Time", user = "User", name_short = "Name",
      type_badge = "Type", n_files = "Files", n_proteins = "Proteins",
      n_de_proteins = "DE", status_badge = "Status", duration_fmt = "Duration",
      project = "Project", actions = "", compare_cb = "\u2194")

    child_row_js <- DT::JS("
      // Compare checkbox logic
      window._delimp_updateCompare = function() {
        var checked = $('.history-compare-cb:checked');
        var btn = $('#history_compare_btn');
        if (checked.length === 2) {
          btn.show();
        } else {
          btn.hide();
        }
        // Enforce max 2: uncheck others when 2 already selected
        if (checked.length > 2) {
          checked.last().prop('checked', false);
          btn.show();
        }
      };
      // Child row expansion (skip checkbox and actions columns)
      table.on('click', 'tbody tr td:not(:last-child):not(:nth-last-child(2))', function() {
        // Ignore clicks on checkboxes
        if ($(event.target).is('input[type=checkbox]')) return;
        var tr = $(this).closest('tr');
        var row = table.row(tr);
        if (row.child.isShown()) {
          row.child.hide();
          tr.removeClass('shown');
        } else {
          var d = row.data()[0];
          if (d && d.length > 0) {
            row.child(d).show();
            tr.addClass('shown');
          }
        }
      });
    ")

    datatable(log[, display_cols, drop = FALSE],
      options = list(
        pageLength = 15, dom = "ftp",
        columnDefs = list(
          list(visible = FALSE, targets = 0),
          list(orderable = FALSE, targets = no_sort_cols),
          list(width = "30px", targets = which(display_cols == "compare_cb") - 1)
        ),
        order = list(list(1, "desc"))
      ),
      callback = child_row_js,
      rownames = FALSE, escape = FALSE,
      colnames = unname(col_names[display_cols]))
  }, server = FALSE)

  # CSV export
  output$history_export_csv <- downloadHandler(
    filename = function() paste0("delimp_history_", format(Sys.time(), "%Y%m%d"), ".csv"),
    content = function(file) {
      log <- activity_log_read()
      write.csv(log, file, row.names = FALSE)
    }
  )

  # ==========================================================================
  #   History action handlers
  # ==========================================================================

  # View search_info.md
  observeEvent(input$history_log_click, {
    out_dir <- input$history_log_click$od
    if (is.null(out_dir) || !nzchar(out_dir)) return()

    log <- activity_log_read()
    search_name <- ""
    if (nrow(log) > 0) {
      match_row <- which(log$output_dir == out_dir)
      if (length(match_row) > 0) search_name <- log$search_name[match_row[1]] %||% ""
    }

    info_content <- ""
    info_path <- file.path(out_dir, "search_info.md")

    cfg <- if (isTRUE(values$ssh_connected) && nzchar(input$ssh_host %||% ""))
      list(host = input$ssh_host, user = input$ssh_user,
           port = input$ssh_port %||% 22, key_path = input$ssh_key_path) else NULL
    if (!is.null(cfg)) {
      result <- tryCatch(
        ssh_exec(cfg, sprintf("cat %s 2>/dev/null", shQuote(info_path)), timeout = 15),
        error = function(e) list(status = 1, stdout = character()))
      if (result$status == 0 && length(result$stdout) > 0)
        info_content <- paste(result$stdout, collapse = "\n")
    }
    if (!nzchar(info_content) && file.exists(info_path)) {
      info_content <- tryCatch(paste(readLines(info_path, warn = FALSE), collapse = "\n"),
        error = function(e) "")
    }

    if (nzchar(info_content)) {
      showModal(modalDialog(
        title = paste("Search Info:", search_name),
        size = "l", easyClose = TRUE,
        tags$pre(style = "max-height:500px;overflow-y:auto;white-space:pre-wrap;font-size:0.85em;",
          info_content),
        footer = modalButton("Close")))
    } else {
      showNotification("No search_info.md found for this entry.", type = "warning")
    }
  }, ignoreInit = TRUE)

  # Import settings from a history row
  observeEvent(input$history_settings_click, {
    out_dir <- input$history_settings_click$od
    if (is.null(out_dir) || !nzchar(out_dir)) return()

    log <- activity_log_read()
    if (nrow(log) == 0) return()
    match_row <- which(!is.na(log$output_dir) & log$output_dir == out_dir)
    if (length(match_row) == 0) return()
    row <- log[match_row[length(match_row)], ]

    tryCatch({
      if (!is.na(row$mass_acc_mode)) updateRadioButtons(session, "mass_acc_mode", selected = row$mass_acc_mode)
      if (!is.na(row$mass_acc)) updateNumericInput(session, "diann_mass_acc", value = as.numeric(row$mass_acc))
      if (!is.na(row$mass_acc_ms1)) updateNumericInput(session, "diann_mass_acc_ms1", value = as.numeric(row$mass_acc_ms1))
      if (!is.na(row$scan_window)) updateNumericInput(session, "diann_scan_window", value = as.integer(row$scan_window))
      if (!is.na(row$enzyme) && nzchar(row$enzyme)) updateSelectInput(session, "diann_enzyme", selected = row$enzyme)
      if (!is.na(row$search_mode)) updateRadioButtons(session, "search_mode", selected = row$search_mode)
      if (!is.na(row$normalization)) updateRadioButtons(session, "diann_normalization", selected = row$normalization)
      if (!is.na(row$extra_cli_flags) && nzchar(row$extra_cli_flags))
        updateTextInput(session, "extra_cli_flags", value = row$extra_cli_flags)

      showNotification(sprintf("Imported settings from '%s'", row$search_name %||% "search"),
        type = "message", duration = 5)
    }, error = function(e) showNotification(sprintf("Failed: %s", e$message), type = "error"))
  }, ignoreInit = TRUE)

  # Unified Load handler — tries session.rds first, falls back to report.parquet
  observeEvent(input$history_load_click, {
    sf <- input$history_load_click$sf
    od <- input$history_load_click$od
    if (is.null(od) || !nzchar(od)) return()

    cfg <- if (isTRUE(values$ssh_connected) && nzchar(input$ssh_host %||% ""))
      list(host = input$ssh_host, user = input$ssh_user,
           port = input$ssh_port %||% 22, key_path = input$ssh_key_path) else NULL

    tryCatch({
      # Try to load session.rds first (full state restore)
      rds_loaded <- FALSE
      rds_path <- NULL

      if (!is.null(sf) && nzchar(sf %||% "")) {
        if (!is.null(cfg) && isTRUE(values$ssh_connected)) {
          find_result <- tryCatch(ssh_exec(cfg, paste("ls", shQuote(sf), "2>/dev/null")),
            error = function(e) list(status = 1))
          if (find_result$status == 0) {
            showNotification("Downloading session...", type = "message", duration = 30, id = "hist_load")
            local_rds <- file.path(tempdir(), basename(sf))
            dl <- scp_download(cfg, sf, local_rds)
            if (dl$status == 0) rds_path <- local_rds
          }
        } else if (file.exists(sf)) {
          rds_path <- sf
        }
      }

      if (!is.null(rds_path)) {
        showNotification("Restoring session...", type = "message", duration = 30, id = "hist_load")
        session_data <- readRDS(rds_path)

        if (all(c("raw_data", "metadata", "fit") %in% names(session_data))) {
          # Full restore
          values$raw_data   <- session_data$raw_data
          values$metadata   <- session_data$metadata
          values$fit        <- session_data$fit
          values$y_protein  <- session_data$y_protein
          values$dpc_fit    <- session_data$dpc_fit
          values$design     <- session_data$design
          values$qc_stats   <- session_data$qc_stats
          values$gsea_results <- session_data$gsea_results
          if (!is.null(session_data$gsea_results_cache)) values$gsea_results_cache <- session_data$gsea_results_cache
          values$gsea_last_contrast <- session_data$gsea_last_contrast
          values$gsea_last_org_db <- session_data$gsea_last_org_db
          values$color_plot_by_de <- session_data$color_plot_by_de %||% FALSE
          values$cov1_name  <- session_data$cov1_name %||% "Covariate1"
          values$cov2_name  <- session_data$cov2_name %||% "Covariate2"
          values$phospho_detected <- session_data$phospho_detected
          values$phospho_site_matrix <- session_data$phospho_site_matrix
          values$phospho_site_info <- session_data$phospho_site_info
          values$phospho_fit <- session_data$phospho_fit
          values$phospho_site_matrix_filtered <- session_data$phospho_site_matrix_filtered
          values$phospho_input_mode <- session_data$phospho_input_mode
          values$ksea_results <- session_data$ksea_results
          values$ksea_last_contrast <- session_data$ksea_last_contrast
          values$phospho_fasta_sequences <- session_data$phospho_fasta_sequences
          values$phospho_annotations <- session_data$phospho_annotations
          if (!is.null(session_data$mofa_object)) {
            values$mofa_view_configs <- session_data$mofa_view_configs %||% list()
            values$mofa_views <- session_data$mofa_views %||% list()
            values$mofa_view_fits <- session_data$mofa_view_fits %||% list()
            values$mofa_sample_metadata <- session_data$mofa_sample_metadata
            values$mofa_object <- session_data$mofa_object
            values$mofa_factors <- session_data$mofa_factors
            values$mofa_weights <- session_data$mofa_weights
            values$mofa_variance_explained <- session_data$mofa_variance_explained
            values$mofa_last_run_params <- session_data$mofa_last_run_params
          }
          if (!is.null(session_data$instrument_metadata)) values$instrument_metadata <- session_data$instrument_metadata
          if (!is.null(session_data$tic_traces)) {
            values$tic_traces <- session_data$tic_traces
            values$tic_metrics <- session_data$tic_metrics
          }
          if (!is.null(session_data$comparator_results)) {
            values$comparator_results <- session_data$comparator_results
            values$comparator_run_a <- session_data$comparator_run_a
            values$comparator_run_b <- session_data$comparator_run_b
            values$comparator_mode <- session_data$comparator_mode
            values$comparator_gemini_narrative <- session_data$comparator_gemini_narrative
            values$comparator_mofa <- session_data$comparator_mofa
          }
          if (!is.null(session_data$diann_search_settings)) {
            values$diann_search_settings <- session_data$diann_search_settings
            values$diann_search_settings$output_dir <- od
          } else {
            values$diann_search_settings <- list(output_dir = od)
          }
          values$repro_log <- session_data$repro_log %||% values$repro_log
          if (!is.null(values$fit)) {
            cn <- colnames(values$fit$contrasts)
            sel <- session_data$contrast %||% cn[1]
            updateSelectInput(session, "contrast_selector", choices = cn, selected = sel)
            updateSelectInput(session, "contrast_selector_signal", choices = cn, selected = sel)
            updateSelectInput(session, "contrast_selector_grid", choices = cn, selected = sel)
            updateSelectInput(session, "contrast_selector_pvalue", choices = cn, selected = sel)
          }
          if (!is.null(session_data$logfc_cutoff)) updateSliderInput(session, "logfc_cutoff", value = session_data$logfc_cutoff)
          if (!is.null(session_data$q_cutoff)) updateNumericInput(session, "q_cutoff", value = session_data$q_cutoff)

          rds_loaded <- TRUE
          removeNotification("hist_load")
          showNotification(
            sprintf("Session restored! (saved %s)", format(session_data$saved_at %||% Sys.time(), "%Y-%m-%d %H:%M")),
            type = "message", duration = 5)
          nav_select("main_tabs", "Data Overview", session = session)
        }
      }

      # Fallback: load report.parquet only
      if (!rds_loaded) {
        report_path <- NULL
        if (!is.null(cfg) && isTRUE(values$ssh_connected)) {
          for (pattern in c("report.parquet", "no_norm_report.parquet")) {
            remote_path <- file.path(od, pattern)
            check <- tryCatch(ssh_exec(cfg, paste("ls", shQuote(remote_path), "2>/dev/null")),
              error = function(e) list(status = 1))
            if (check$status == 0) {
              showNotification("Downloading report...", type = "message", duration = 30, id = "hist_load")
              local_path <- file.path(tempdir(), paste0("hist_", basename(od), "_", pattern))
              dl <- scp_download(cfg, remote_path, local_path)
              if (dl$status == 0) { report_path <- local_path; break }
            }
          }
        } else {
          for (f in c(file.path(od, "report.parquet"), file.path(od, "no_norm_report.parquet"))) {
            if (file.exists(f)) { report_path <- f; break }
          }
        }

        if (is.null(report_path)) {
          removeNotification("hist_load")
          showNotification("No session.rds or report.parquet found.", type = "error", duration = 8)
          return()
        }

        showNotification("Reading DIA-NN report...", type = "message", duration = 30, id = "hist_load")
        raw_data <- suppressMessages(suppressWarnings(limpa::readDIANN(report_path, format = "parquet")))
        values$raw_data <- raw_data
        values$qc_stats <- get_diann_stats_r(report_path)
        values$uploaded_report_path <- report_path
        values$original_report_name <- basename(report_path)
        values$metadata <- data.frame(
          ID = seq_along(colnames(raw_data$E)),
          File.Name = colnames(raw_data$E),
          Group = "", Batch = "", Covariate1 = "", Covariate2 = "",
          stringsAsFactors = FALSE)

        # Reconstruct search settings from activity log
        tryCatch({
          log <- activity_log_read()
          match_rows <- log[!is.na(log$output_dir) & log$output_dir == od, , drop = FALSE]
          if (nrow(match_rows) > 0) {
            row <- match_rows[nrow(match_rows), ]
            values$diann_search_settings <- list(
              output_dir = od, search_name = row$search_name %||% NA,
              fasta_files = if (!is.na(row$fasta_files)) strsplit(row$fasta_files, ",\\s*")[[1]] else character(),
              fasta_seq_count = if (!is.na(row$fasta_seq_count)) as.integer(row$fasta_seq_count) else NULL,
              search_mode = row$search_mode %||% NA,
              n_raw_files = row$n_files %||% NA, raw_file_type = NA,
              normalization = row$normalization %||% NA,
              search_params = list(
                enzyme = row$enzyme %||% "K*,R*", missed_cleavages = 1,
                mass_acc_mode = row$mass_acc_mode %||% "auto",
                mass_acc = row$mass_acc %||% NA, mass_acc_ms1 = row$mass_acc_ms1 %||% NA,
                scan_window = row$scan_window %||% NA, mbr = row$mbr %||% NA,
                normalization = row$normalization %||% NA,
                extra_cli_flags = row$extra_cli_flags %||% "",
                qvalue = 0.01, max_var_mods = 1,
                min_pep_len = 7, max_pep_len = 30,
                min_pr_mz = 300, max_pr_mz = 1800,
                min_pr_charge = 1, max_pr_charge = 4,
                min_fr_mz = 200, max_fr_mz = 1800,
                unimod4 = TRUE, met_excision = TRUE,
                mod_met_ox = TRUE, mod_nterm_acetyl = FALSE,
                rt_profiling = TRUE
              ), imported_from_log = TRUE)
          } else {
            values$diann_search_settings <- list(output_dir = od)
          }
        }, error = function(e) values$diann_search_settings <- list(output_dir = od))

        removeNotification("hist_load")
        showNotification(sprintf("Loaded %s: %d proteins, %d samples",
          basename(report_path), nrow(raw_data$E), ncol(raw_data$E)), type = "message", duration = 10)
        nav_select("main_tabs", "Data Overview", session = session)
      }
    }, error = function(e) {
      removeNotification("hist_load")
      showNotification(sprintf("Load failed: %s", e$message), type = "error", duration = 10)
    })
  }, ignoreInit = TRUE)

  # ==========================================================================
  #   Notes modal — triggered on search completion or manual edit
  # ==========================================================================

  observeEvent(values$pending_notes_od, {
    od <- values$pending_notes_od
    if (is.null(od) || !nzchar(od)) return()

    search_name <- values$pending_notes_name %||% basename(od)

    showModal(modalDialog(
      title = sprintf("Add Notes: %s", search_name),
      size = "m", easyClose = TRUE,
      textAreaInput("notes_text", "Notes (optional)",
        value = "", rows = 4,
        placeholder = "Describe what you were testing, key parameters, or observations..."),
      footer = tagList(
        actionButton("notes_skip_btn", "Skip", class = "btn-default"),
        actionButton("notes_save_btn", "Save", class = "btn-primary")
      )
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$notes_save_btn, {
    od <- values$pending_notes_od
    notes_text <- input$notes_text %||% ""
    if (!is.null(od) && nzchar(od) && nzchar(notes_text)) {
      tryCatch({
        update_activity(od, list(notes = notes_text))
        history_refresh(isolate(history_refresh()) + 1)
      }, error = function(e) NULL)
    }
    values$pending_notes_od <- NULL
    values$pending_notes_name <- NULL
    removeModal()
  }, ignoreInit = TRUE)

  observeEvent(input$notes_skip_btn, {
    values$pending_notes_od <- NULL
    values$pending_notes_name <- NULL
    removeModal()
  }, ignoreInit = TRUE)

  # Manual notes edit from history table
  observeEvent(input$history_notes_click, {
    id <- input$history_notes_click$id
    od <- input$history_notes_click$od
    existing_notes <- input$history_notes_click$notes %||% ""
    name <- input$history_notes_click$name %||% ""

    showModal(modalDialog(
      title = sprintf("Edit Notes: %s", name),
      size = "m", easyClose = TRUE,
      tags$input(type = "hidden", id = "edit_notes_id", value = id %||% ""),
      tags$input(type = "hidden", id = "edit_notes_od", value = od %||% ""),
      textAreaInput("edit_notes_text", "Notes", value = existing_notes, rows = 4),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("edit_notes_save_btn", "Save", class = "btn-primary")
      )
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$edit_notes_save_btn, {
    id <- input$edit_notes_id
    od <- input$edit_notes_od
    notes_text <- input$edit_notes_text %||% ""

    tryCatch({
      if (!is.null(id) && nzchar(id)) {
        update_activity_by_id(id, list(notes = notes_text))
      } else if (!is.null(od) && nzchar(od)) {
        update_activity(od, list(notes = notes_text))
      }
      history_refresh(isolate(history_refresh()) + 1)
      removeModal()
      showNotification("Notes saved.", type = "message", duration = 3)
    }, error = function(e) showNotification(sprintf("Failed: %s", e$message), type = "error"))
  }, ignoreInit = TRUE)

  # ==========================================================================
  #   Project assignment
  # ==========================================================================

  assign_od <- reactiveVal(NULL)

  observeEvent(input$history_project_click, {
    od <- input$history_project_click$od
    if (is.null(od) || !nzchar(od)) return()

    assign_od(od)
    existing_projects <- get_projects()

    showModal(modalDialog(
      title = "Assign to Project",
      size = "s", easyClose = TRUE,
      selectizeInput("assign_project_name", "Project",
        choices = existing_projects,
        options = list(create = TRUE, placeholder = "Select or type new project name...")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("assign_project_confirm", "Assign", class = "btn-primary")
      )
    ))
  }, ignoreInit = TRUE)

  observeEvent(input$assign_project_confirm, {
    proj_name <- input$assign_project_name
    od <- assign_od()

    if (is.null(proj_name) || !nzchar(proj_name)) {
      showNotification("Please enter a project name.", type = "warning")
      return()
    }
    if (is.null(od) || !nzchar(od)) return()

    tryCatch({
      set_project(od, proj_name)
      history_refresh(isolate(history_refresh()) + 1)
      removeModal()
      showNotification(sprintf("Assigned to project: %s", proj_name), type = "message")
    }, error = function(e) showNotification(sprintf("Failed: %s", e$message), type = "error"))
  }, ignoreInit = TRUE)

  # ==========================================================================
  #   Compare two analyses — load RDS files and send to Run Comparator
  # ==========================================================================
  observeEvent(input$history_compare_click, {
    od_a <- input$history_compare_click$od_a
    od_b <- input$history_compare_click$od_b
    sf_a <- input$history_compare_click$sf_a
    sf_b <- input$history_compare_click$sf_b
    name_a <- input$history_compare_click$name_a
    name_b <- input$history_compare_click$name_b

    if (is.null(od_a) || is.null(od_b) || !nzchar(od_a) || !nzchar(od_b)) {
      showNotification("Please select exactly 2 analyses to compare.", type = "warning")
      return()
    }

    cfg <- if (isTRUE(values$ssh_connected) && nzchar(input$ssh_host %||% ""))
      list(host = input$ssh_host, user = input$ssh_user,
           port = input$ssh_port %||% 22, key_path = input$ssh_key_path) else NULL

    showNotification("Loading sessions for comparison...", type = "message", duration = 30, id = "compare_load")

    tryCatch({
      # Helper to load an RDS from local or remote — tries explicit path, then output_dir/session.rds
      load_rds <- function(od, sf = NULL) {
        candidates <- unique(c(
          if (!is.null(sf) && nzchar(sf %||% "")) sf,
          file.path(od, "session.rds")
        ))
        for (rds_path in candidates) {
          if (!is.null(cfg) && isTRUE(values$ssh_connected)) {
            check <- tryCatch(ssh_exec(cfg, paste("ls", shQuote(rds_path), "2>/dev/null")),
              error = function(e) list(status = 1))
            if (check$status == 0) {
              local_rds <- file.path(tempdir(), paste0("compare_", basename(od), "_", basename(rds_path)))
              dl <- scp_download(cfg, rds_path, local_rds)
              if (dl$status == 0) return(local_rds)
            }
          } else if (file.exists(rds_path)) {
            return(rds_path)
          }
        }
        stop(sprintf("No session.rds found for %s. Run the LIMPA pipeline on this data first, then the RDS will be auto-saved.", basename(od)))
      }

      # Helper to find and parse DIA-NN log from output_dir
      find_diann_log <- function(od) {
        tryCatch({
          log_dir <- file.path(od, "logs")
          log_files <- NULL
          if (!is.null(cfg) && isTRUE(values$ssh_connected)) {
            result <- tryCatch(ssh_exec(cfg, paste("ls", shQuote(log_dir), "2>/dev/null")),
              error = function(e) list(status = 1, stdout = ""))
            if (result$status == 0 && nzchar(result$stdout)) {
              files <- trimws(strsplit(result$stdout, "\n")[[1]])
              log_files <- files[grepl("\\.(log|out)$", files)]
            }
          } else if (dir.exists(log_dir)) {
            log_files <- list.files(log_dir, pattern = "\\.(log|out)$")
          }
          if (is.null(log_files) || length(log_files) == 0) return(NULL)

          # Try each log file until one parses successfully
          for (lf in log_files) {
            remote_path <- file.path(log_dir, lf)
            local_path <- NULL
            if (!is.null(cfg) && isTRUE(values$ssh_connected)) {
              local_path <- file.path(tempdir(), paste0("compare_log_", lf))
              dl <- scp_download(cfg, remote_path, local_path)
              if (dl$status != 0) next
            } else if (file.exists(remote_path)) {
              local_path <- remote_path
            } else next

            parsed <- tryCatch(parse_diann_log(local_path), error = function(e) list(success = FALSE))
            if (isTRUE(parsed$success)) return(parsed)
          }
          NULL
        }, error = function(e) NULL)
      }

      # Load both RDS files
      rds_a <- load_rds(od_a, sf_a)
      rds_b <- load_rds(od_b, sf_b)

      # Parse both sessions
      parsed_a <- parse_delimp_session(rds_path = rds_a)
      parsed_b <- parse_delimp_session(rds_path = rds_b)

      # Set comparator mode and data
      updateRadioButtons(session, "comparator_mode", selected = "delimp_delimp")
      updateRadioButtons(session, "comparator_run_a_source", selected = "file")

      # Set Run A and Run B via internal reactiveVals
      values$comparator_run_a <- parsed_a
      values$comparator_run_b <- parsed_b
      values$comparator_compare_from_history <- list(
        od_a = od_a, od_b = od_b,
        name_a = name_a %||% basename(od_a),
        name_b = name_b %||% basename(od_b)
      )

      # Try to find DIA-NN logs for enriched comparison
      log_a <- find_diann_log(od_a)
      log_b <- find_diann_log(od_b)
      values$comparator_diann_log_a <- log_a
      values$comparator_diann_log_b <- log_b

      removeNotification("compare_load")
      showNotification(
        sprintf("Loaded: %s vs %s. Navigating to Comparator...",
          name_a %||% basename(od_a), name_b %||% basename(od_b)),
        type = "message", duration = 5)
      nav_select("main_tabs", "comparator_tab", session = session)

    }, error = function(e) {
      removeNotification("compare_load")
      showNotification(sprintf("Compare failed: %s", e$message), type = "error", duration = 10)
    })
  }, ignoreInit = TRUE)

}
