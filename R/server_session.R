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
      title = tagList(icon("question-circle"), " What is the Consistency Table?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What is %CV?"),
        p("The Coefficient of Variation (CV) measures how consistent a protein's expression is within each group. ",
          "It's calculated as (standard deviation / mean) \u00d7 100. Lower CV = more reproducible measurement."),
        tags$h6("Why this matters"),
        p("A protein that is differentially expressed AND has low CV across replicates is a stronger biomarker candidate. ",
          "High CV means the measurement is noisy and the DE result may not be reproducible."),
        tags$h6("Reading the table"),
        tags$ul(
          tags$li("Proteins are ranked by average %CV across all groups (lowest first)"),
          tags$li("Each group's %CV is shown in its own column"),
          tags$li(strong("Top candidates: "), "Significant proteins (adj.P.Val < 0.05) with %CV < 20% in all groups")
        )
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

  # --- QC Trends Info Modal ---
  observeEvent(input$qc_trends_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " What are QC Trends?"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What these plots show"),
        p("QC trend plots track key quality metrics across all samples in your experiment. ",
          "They help identify systematic issues like instrument drift, sample degradation, or injection problems."),
        tags$h6("The metrics"),
        tags$ul(
          tags$li(strong("Precursors: "), "Number of peptide precursors identified. A sudden drop may indicate instrument issues or sample problems."),
          tags$li(strong("Proteins: "), "Number of proteins quantified. Should be relatively stable across runs."),
          tags$li(strong("MS1 Signal: "), "Overall MS1 intensity. Gradual decline may indicate column degradation or source contamination.")
        ),
        tags$h6("Sort order"),
        p(strong("Run Order"), " shows samples in the order they were acquired \u2014 useful for spotting instrument drift over time. ",
          strong("Group"), " sorts by experimental condition \u2014 useful for comparing groups side by side."),
        tags$h6("What to look for"),
        tags$ul(
          tags$li("Stable, flat trends are ideal"),
          tags$li("Sudden drops or spikes in individual samples flag potential outliers"),
          tags$li("Gradual downward trends may indicate column or instrument degradation")
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
        # Save timestamp & version
        saved_at   = Sys.time(),
        app_version = "DE-LIMP v2.4"
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
    }, error = function(e) {
      showNotification(paste("Error loading session:", e$message), type = "error")
    })
  })

  # ============================================================================
  #      Methodology Text
  # ============================================================================

  output$methodology_text <- renderText({
    req(values$fit)

    methodology <- paste(
      "METHODOLOGY\n",
      "===========\n",
      "Data Processing and Statistical Analysis Pipeline\n",
      "---------------------------------------------------\n\n",

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
      "When performed, over-representation analysis (ORA) was conducted using clusterProfiler.\n",
      "Significant proteins (adj.P.Val < 0.05) were mapped to Gene Ontology (GO) terms, and\n",
      "enrichment was tested using hypergeometric distribution with FDR correction.\n\n\n",

      "SOFTWARE AND PACKAGES\n",
      "---------------------\n",
      "Primary analysis: limpa R package (Bioconductor 3.22+)\n",
      "Statistical framework: limma (Linear Models for Microarray and RNA-Seq Data)\n",
      "Data manipulation: dplyr, tidyr\n",
      "Visualization: ggplot2, ComplexHeatmap, plotly\n",
      "Enrichment: clusterProfiler, enrichplot\n",
      sprintf("R version: %s\n\n\n", R.version.string),

      "REFERENCES\n",
      "----------\n",
      "\u2022 limpa package: Bioconductor (https://bioconductor.org/packages/limpa/)\n",
      "\u2022 DPC normalization: Designed for DIA-NN proteomics data\n",
      "\u2022 limma: Ritchie ME, et al. (2015) Nucleic Acids Research 43(7):e47\n",
      "\u2022 Empirical Bayes: Smyth GK (2004) Statistical Applications in Genetics and\n",
      "  Molecular Biology 3:Article3\n",
      "\u2022 FDR control: Benjamini Y, Hochberg Y (1995) Journal of the Royal Statistical\n",
      "  Society 57(1):289-300\n",
      "\u2022 clusterProfiler: Yu G, et al. (2012) OMICS 16(5):284-287\n\n\n",

      "CITATION\n",
      "--------\n",
      "If you use this analysis in your research, please cite:\n",
      "\u2022 The limpa package (Bioconductor)\n",
      "\u2022 The limma package: Ritchie ME, et al. (2015) Nucleic Acids Research\n",
      "\u2022 DIA-NN: Demichev V, et al. (2020) Nature Methods 17:41-44",

      sep = ""
    )

    methodology
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

}
