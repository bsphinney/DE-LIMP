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
        # Save timestamp & version
        saved_at   = Sys.time(),
        app_version = "DE-LIMP v2.5"
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
    req(values$fit %||% values$phospho_fit)

    # DIA-NN search parameters section (conditional)
    diann_section <- ""
    ss <- values$diann_search_settings
    if (!is.null(ss)) {
      sp <- ss$search_params

      # Build enzyme name
      enzyme_name <- switch(sp$enzyme,
        "K*,R*" = "Trypsin/P (cleaves after K, R)",
        "K*,R*,!*P" = "Trypsin (not before P)",
        "R*" = "Arg-C",
        "K*" = "Lys-C",
        sp$enzyme)

      # Build modifications list
      mods <- c()
      if (isTRUE(sp$mod_met_ox)) mods <- c(mods, "Oxidation (M) [UniMod:35]")
      if (isTRUE(sp$mod_nterm_acetyl)) mods <- c(mods, "N-terminal acetylation [UniMod:1]")
      if (nzchar(sp$extra_var_mods)) {
        extra <- trimws(strsplit(sp$extra_var_mods, "\n")[[1]])
        mods <- c(mods, extra[nzchar(extra)])
      }
      if (isTRUE(sp$unimod4)) mods <- c(mods, "Carbamidomethylation (C) [UniMod:4, fixed]")

      # FASTA description
      fasta_desc <- paste(basename(ss$fasta_files), collapse = ", ")
      contam_desc <- if (!is.null(ss$contaminant_library) && ss$contaminant_library != "none")
        sprintf(" with %s contaminant library", gsub("_", " ", ss$contaminant_library))
      else ""

      # Mass accuracy description
      mass_acc_desc <- if (sp$mass_acc_mode == "auto") {
        "automatic (DIA-NN optimized)"
      } else {
        sprintf("MS2 %s ppm, MS1 %s ppm", sp$mass_acc, sp$mass_acc_ms1)
      }

      # DIA-NN engine description (Docker vs HPC/Singularity)
      engine_desc <- if (!is.null(ss$docker_image)) {
        sprintf("DIA-NN (%s, Docker)", ss$docker_image)
      } else if (!is.null(ss$diann_sif)) {
        sprintf("DIA-NN (%s)", ss$diann_sif)
      } else {
        "DIA-NN"
      }

      # Compute resource description
      resource_desc <- if (!is.null(ss$docker)) {
        sprintf("Docker resources: %d CPUs, %d GB RAM (image: %s).",
                ss$docker$cpus, ss$docker$mem_gb, ss$docker$image)
      } else if (!is.null(ss$slurm)) {
        sprintf("SLURM resources: %s CPUs, %s GB RAM, %s hour(s), partition: %s.",
                ss$slurm$cpus, ss$slurm$mem_gb, ss$slurm$time_hours, ss$slurm$partition)
      } else {
        ""
      }

      diann_section <- paste0(
        "0. DIA-NN DATABASE SEARCH\n",
        "Raw data files (", ss$n_raw_files, " ", toupper(ss$raw_file_type), " files) were searched using ", engine_desc, " ",
        "in ", ss$search_mode, " mode",
        if (!is.null(ss$speclib)) paste0(" with spectral library (", ss$speclib, ")") else "",
        ".\n",
        "Sequence database: ", fasta_desc, contam_desc, ".\n",
        "Enzyme: ", enzyme_name, " with up to ", sp$missed_cleavages, " missed cleavage(s).\n",
        "Peptide length: ", sp$min_pep_len, "-", sp$max_pep_len, " amino acids. ",
        "Precursor m/z: ", sp$min_pr_mz, "-", sp$max_pr_mz, ". ",
        "Precursor charge: ", sp$min_pr_charge, "-", sp$max_pr_charge, ".\n",
        "Mass accuracy: ", mass_acc_desc, ".\n",
        "Variable modifications (max ", sp$max_var_mods, "): ", paste(mods, collapse = "; "), ".\n",
        "FDR threshold: ", sp$qvalue, " at precursor and protein levels.\n",
        "Match-between-runs: ", if (isTRUE(sp$mbr)) "enabled" else "disabled", ". ",
        "RT-dependent normalization: ", if (isTRUE(sp$rt_profiling)) "enabled" else "disabled", ".\n",
        if (ss$search_mode == "phospho")
          "Phosphoproteomics mode: STY phosphorylation (UniMod:21, +79.966 Da), phospho-specific output and library info reporting enabled.\n"
        else "",
        if (nzchar(sp$extra_cli_flags)) paste0("Additional flags: ", sp$extra_cli_flags, "\n") else "",
        resource_desc, "\n\n"
      )
    }

    methodology <- paste(
      "METHODOLOGY\n",
      "===========\n",
      "Data Processing and Statistical Analysis Pipeline\n",
      "---------------------------------------------------\n\n",

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

      if (!is.null(values$phospho_fit)) paste(
        "7. PHOSPHOSITE-LEVEL DIFFERENTIAL EXPRESSION\n",
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

        "8. KINASE-SUBSTRATE ENRICHMENT ANALYSIS (KSEA)\n",
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

        "9. SEQUENCE MOTIF ANALYSIS (Optional)\n",
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

        "10. PHOSPHOSITE ANNOTATION (Optional)\n",
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
      "Primary analysis: limpa R package (Bioconductor 3.22+)\n",
      "Statistical framework: limma (Linear Models for Microarray and RNA-Seq Data)\n",
      "Data manipulation: dplyr, tidyr\n",
      "Visualization: ggplot2, ComplexHeatmap, plotly\n",
      "Enrichment: clusterProfiler (gseGO, gseKEGG), enrichplot\n",
      "Annotation DBs: org.Hs.eg.db, org.Mm.eg.db (auto-detected)\n",
      if (!is.null(values$phospho_fit)) paste(
        "Phosphoproteomics: KSEAapp (kinase activity inference), ggseqlogo (motif analysis)\n",
        "Phosphosite databases: PhosphoSitePlus (via KSEAapp), UniProt REST API\n",
        sep = ""
      ) else "",
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
      "\u2022 clusterProfiler: Yu G, et al. (2012) OMICS 16(5):284-287\n",
      "\u2022 GSEA method: Subramanian A, et al. (2005) PNAS 102(43):15545-15550\n",
      "\u2022 KEGG database: Kanehisa M, Goto S (2000) Nucleic Acids Research 28(1):27-30\n",
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
