# ==============================================================================
#  SERVER MODULE -- DDA Search (Sage pipeline on Hive)
#  Called from app.R as: server_dda(input, output, session, values, add_to_log)
# ==============================================================================

server_dda <- function(input, output, session, values, add_to_log) {

  # Load sage_bin path from config
  config <- tryCatch(yaml::read_yaml("config.yml"), error = function(e) list())
  sage_bin <- config$tools$sage_bin %||%
    "/quobyte/proteomics-grp/de-limp/cascadia/sage-v0.14.7-x86_64-unknown-linux-gnu/sage"
  slurm_account   <- config$slurm$account   %||% "genome-center-grp"
  slurm_partition <- config$slurm$partition  %||% "high"

  # --- Mode observer: sync input to reactive values ---
  observeEvent(input$acquisition_mode, {
    values$acquisition_mode <- input$acquisition_mode
  })

  # --- Contextual label for the mode switcher ---
  output$mode_context_label <- renderUI({
    mode <- input$acquisition_mode %||% "dia"
    switch(mode,
      "dia" = tags$span(
        style = "font-size: 12px; color: #0d6efd; font-style: italic;",
        icon("circle-check", style = "color: #198754;"),
        " DIA-NN + limpa pipeline"
      ),
      "dda" = tags$span(
        style = "font-size: 12px; color: #6c757d; font-style: italic;",
        icon("circle-info", style = "color: #0d6efd;"),
        " Sage + Casanovo pipeline"
      ),
      "xlms" = tags$span(
        style = "font-size: 12px; color: #6c757d; font-style: italic;",
        icon("diagram-project", style = "color: #6f42c1;"),
        " MeroX + xiSearch + network"
      )
    )
  })

  # ============================================================================
  #    SSH config helper (reuses DIA-NN search SSH settings)
  # ============================================================================
  dda_ssh_config <- reactive({
    req(values$ssh_connected)
    list(
      host     = isolate(input$ssh_host),
      user     = isolate(input$ssh_user),
      port     = isolate(input$ssh_port) %||% 22,
      key_path = isolate(input$ssh_key_path),
      modules  = isolate(input$ssh_modules) %||% ""
    )
  })

  # ============================================================================
  #    File scan: list .d files in remote directory
  # ============================================================================
  observeEvent(input$dda_scan_files, {
    req(values$ssh_connected, input$dda_raw_dir)
    raw_dir <- trimws(input$dda_raw_dir)
    if (!nzchar(raw_dir)) {
      showNotification("Please enter a raw file directory path.", type = "warning")
      return()
    }

    ssh_cfg <- dda_ssh_config()
    result <- ssh_exec(ssh_cfg,
      paste0("ls -1d ", shQuote(raw_dir), "/*.d 2>/dev/null | head -200"),
      timeout = 15)

    if (result$status != 0 || length(result$stdout) == 0 ||
        all(!nzchar(trimws(result$stdout)))) {
      showNotification("No .d directories found in the specified path.", type = "warning")
      values$dda_raw_files <- character(0)
      return()
    }

    files <- trimws(result$stdout)
    files <- files[nzchar(files)]
    values$dda_raw_files <- files
    showNotification(paste("Found", length(files), ".d files"), type = "message")
  })

  # File list preview
  output$dda_file_list_preview <- renderUI({
    files <- values$dda_raw_files
    if (is.null(files) || length(files) == 0) {
      return(tags$p(style = "color: #6c757d; font-style: italic;",
        "No files scanned yet. Enter a directory path and click Scan."))
    }
    basenames <- basename(files)
    tags$div(
      style = "max-height: 200px; overflow-y: auto; background: #f8f9fa; padding: 8px; border-radius: 6px; font-size: 12px;",
      tags$strong(paste(length(basenames), "files:")),
      tags$ul(style = "margin: 4px 0; padding-left: 16px;",
        lapply(basenames, function(f) tags$li(f))
      )
    )
  })

  # ============================================================================
  #    Submit Sage search
  # ============================================================================
  observeEvent(input$run_dda_search, {
    req(values$ssh_connected)

    raw_dir    <- trimws(input$dda_raw_dir %||% "")
    fasta_path <- trimws(input$dda_fasta_path %||% "")
    exp_name   <- trimws(input$dda_experiment_name %||% "dda_search")

    # Validation
    if (!nzchar(raw_dir)) {
      showNotification("Please enter a raw file directory path.", type = "error")
      return()
    }
    if (!nzchar(fasta_path)) {
      showNotification("Please enter a FASTA file path.", type = "error")
      return()
    }
    if (is.null(values$dda_raw_files) || length(values$dda_raw_files) == 0) {
      showNotification("Please scan for files first.", type = "error")
      return()
    }

    ssh_cfg <- dda_ssh_config()

    # Build output directory on HPC
    output_dir <- file.path(
      "/quobyte/proteomics-grp/de-limp",
      ssh_cfg$user,
      "dda_output",
      gsub("[^a-zA-Z0-9_.-]", "_", exp_name)
    )

    withProgress(message = "Submitting Sage search...", value = 0.1, {
      # Create output dir + logs dir on HPC
      mkdir_result <- ssh_exec(ssh_cfg,
        paste0("mkdir -p ", shQuote(file.path(output_dir, "logs"))),
        timeout = 15)
      if (mkdir_result$status != 0) {
        showNotification("Failed to create output directory on HPC.", type = "error")
        return()
      }
      setProgress(0.2, detail = "Generating Sage config...")

      # Generate sage.json locally, then upload
      local_tmp <- tempdir()
      raw_paths <- values$dda_raw_files

      config_path_local <- generate_sage_config(
        fasta_path       = fasta_path,
        raw_paths        = raw_paths,
        output_dir       = output_dir,
        preset           = input$dda_preset %||% "standard",
        missed_cleavages = input$dda_missed_cleavages %||% 2,
        precursor_tol_ppm = input$dda_precursor_tol %||% 20,
        fragment_tol_da   = input$dda_fragment_tol %||% 0.05,
        min_peaks         = 6
      )

      # Upload sage.json to HPC
      remote_config <- file.path(output_dir, "sage.json")
      scp_result <- scp_upload(ssh_cfg, config_path_local, remote_config)
      if (scp_result$status != 0) {
        showNotification("Failed to upload Sage config to HPC.", type = "error")
        return()
      }
      setProgress(0.4, detail = "Generating sbatch script...")

      # Generate sbatch script
      sbatch_content <- generate_sage_sbatch(
        sage_bin        = sage_bin,
        config_path     = remote_config,
        raw_dir         = raw_dir,
        output_dir      = output_dir,
        experiment_name = exp_name,
        cpus            = input$dda_cpus %||% 32,
        mem_gb          = input$dda_mem %||% 64,
        time_limit      = input$dda_time_limit %||% "02:00:00",
        account         = slurm_account,
        partition       = slurm_partition
      )

      # Write sbatch script locally, upload
      local_sbatch <- file.path(local_tmp, "sage_search.sbatch")
      writeLines(sbatch_content, local_sbatch)
      remote_sbatch <- file.path(output_dir, "sage_search.sbatch")
      scp_upload(ssh_cfg, local_sbatch, remote_sbatch)

      setProgress(0.6, detail = "Submitting to SLURM...")

      # Submit via sbatch
      sbatch_path <- values$ssh_sbatch_path %||% "sbatch"
      submit_result <- ssh_exec(ssh_cfg,
        paste(sbatch_path, shQuote(remote_sbatch)),
        timeout = 30)

      if (submit_result$status != 0) {
        showNotification(
          paste("sbatch submission failed:", paste(submit_result$stdout, collapse = " ")),
          type = "error")
        return()
      }

      # Parse job ID from "Submitted batch job 12345"
      job_line <- grep("Submitted batch job", submit_result$stdout, value = TRUE)
      if (length(job_line) == 0) {
        showNotification("Could not parse job ID from sbatch output.", type = "error")
        return()
      }
      job_id <- trimws(sub(".*Submitted batch job\\s+", "", job_line[1]))
      message("[DDA] Sage job submitted: ", job_id)

      # Store state
      values$dda_job_id     <- job_id
      values$dda_output_dir <- output_dir
      values$dda_status     <- "running"
      values$dda_search_params <- list(
        preset           = input$dda_preset %||% "standard",
        fasta_path       = fasta_path,
        raw_dir          = raw_dir,
        n_files          = length(raw_paths),
        missed_cleavages = input$dda_missed_cleavages %||% 2,
        precursor_tol    = input$dda_precursor_tol %||% 20,
        fragment_tol     = input$dda_fragment_tol %||% 0.05,
        normalization    = input$dda_norm_method %||% "cyclicloess",
        imputation       = input$dda_impute_method %||% "perseus",
        min_valid        = input$dda_min_valid %||% 0.5,
        submitted_at     = Sys.time(),
        sage_bin         = sage_bin
      )

      setProgress(1.0, detail = "Job submitted!")
      showNotification(
        paste("Sage search submitted! Job ID:", job_id),
        type = "message", duration = 10)
    })
  })

  # ============================================================================
  #    Job polling (every 15 seconds when a job is running)
  # ============================================================================
  observe({
    req(values$dda_status == "running", values$dda_job_id, values$ssh_connected)
    invalidateLater(15000)

    ssh_cfg <- isolate(dda_ssh_config())
    job_id  <- isolate(values$dda_job_id)

    result <- tryCatch(
      ssh_exec(ssh_cfg,
        paste0("sacct -j ", job_id, " --format=JobID,State --noheader --parsable2"),
        timeout = 15),
      error = function(e) list(status = 1, stdout = character(0))
    )

    if (result$status != 0 || length(result$stdout) == 0) return()

    # Parse SLURM state -- filter out .extern/.batch substeps
    lines <- trimws(result$stdout)
    lines <- lines[nzchar(lines)]
    main_lines <- lines[!grepl("\\.", lines)]  # exclude substeps

    if (length(main_lines) == 0) return()

    # Get state from the main job line
    parts <- strsplit(main_lines[1], "\\|")[[1]]
    if (length(parts) < 2) return()
    state <- trimws(parts[2])

    if (state %in% c("COMPLETED")) {
      message("[DDA] Sage job completed: ", job_id)
      values$dda_status <- "loading"
      showNotification("Sage search completed! Loading results...",
        type = "message", duration = 8)
      # Trigger result loading
      load_sage_results_from_hpc()
    } else if (state %in% c("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED", "NODE_FAIL")) {
      message("[DDA] Sage job failed: ", job_id, " (", state, ")")
      values$dda_status <- "error"
      showNotification(
        paste("Sage search failed:", state),
        type = "error", duration = 15)
    }
    # PENDING, RUNNING, COMPLETING -> keep polling
  })

  # ============================================================================
  #    Load results from HPC
  # ============================================================================
  load_sage_results_from_hpc <- function() {
    ssh_cfg    <- dda_ssh_config()
    output_dir <- values$dda_output_dir

    withProgress(message = "Loading Sage results...", value = 0.1, {
      local_tmp <- file.path(tempdir(), "sage_results")
      dir.create(local_tmp, recursive = TRUE, showWarnings = FALSE)

      # Download results files
      files_to_get <- c("results.sage.tsv", "lfq.tsv")
      for (f in files_to_get) {
        remote_path <- file.path(output_dir, f)
        local_path  <- file.path(local_tmp, f)
        dl <- scp_download(ssh_cfg, remote_path, local_path)
        if (dl$status != 0) {
          showNotification(paste("Failed to download", f, "from HPC."), type = "error")
          values$dda_status <- "error"
          return()
        }
      }
      setProgress(0.4, detail = "Parsing Sage output...")

      # Also try to download sage report JSON (may have various names)
      for (rpt in c("results.json", "sage_report.json")) {
        remote_rpt <- file.path(output_dir, rpt)
        local_rpt  <- file.path(local_tmp, rpt)
        tryCatch(scp_download(ssh_cfg, remote_rpt, local_rpt), error = function(e) NULL)
      }

      # Parse results
      parsed <- tryCatch(
        parse_sage_results(
          results_path = file.path(local_tmp, "results.sage.tsv"),
          lfq_path     = file.path(local_tmp, "lfq.tsv")
        ),
        error = function(e) {
          showNotification(paste("Error parsing Sage output:", e$message), type = "error")
          NULL
        }
      )

      if (is.null(parsed)) {
        values$dda_status <- "error"
        return()
      }

      setProgress(0.6, detail = "Storing results...")
      values$dda_sage_psms    <- parsed$psms
      values$dda_lfq_wide     <- parsed$lfq_wide
      values$dda_protein_meta <- parsed$protein_meta

      # Parse report JSON if available
      for (rpt in c("results.json", "sage_report.json")) {
        local_rpt <- file.path(local_tmp, rpt)
        if (file.exists(local_rpt)) {
          values$dda_sage_report <- parse_sage_report(local_rpt)
          break
        }
      }

      # Compute QC metrics
      values$dda_qc_metrics <- compute_dda_qc_metrics(parsed$psms, parsed$lfq_wide)

      setProgress(0.8, detail = "Done!")
      values$dda_status <- "done"
      showNotification(
        paste("Sage results loaded:",
              nrow(parsed$lfq_wide), "proteins,",
              nrow(parsed$psms), "PSMs"),
        type = "message", duration = 10)
    })
  }

  # ============================================================================
  #    Run DDA pipeline (normalize + filter + impute + build EList)
  # ============================================================================
  observeEvent(input$run_dda_pipeline, {
    req(values$dda_lfq_wide, values$dda_protein_meta)

    # Check that group assignment exists
    if (is.null(values$metadata) || !"Group" %in% colnames(values$metadata)) {
      showNotification(
        "Please assign sample groups before running the pipeline.",
        type = "warning")
      return()
    }

    withProgress(message = "Running DDA pipeline...", value = 0.1, {
      result <- tryCatch(
        run_dda_pipeline(
          lfq_wide           = values$dda_lfq_wide,
          protein_meta       = values$dda_protein_meta,
          metadata_df        = values$metadata,
          norm_method        = input$dda_norm_method %||% "cyclicloess",
          min_valid_fraction = input$dda_min_valid %||% 0.5,
          impute_method      = input$dda_impute_method %||% "perseus",
          perseus_width      = input$dda_perseus_width %||% 0.3,
          perseus_shift      = input$dda_perseus_shift %||% 1.8
        ),
        error = function(e) {
          showNotification(paste("Pipeline error:", e$message), type = "error")
          NULL
        }
      )

      if (is.null(result)) return()

      setProgress(0.6, detail = "Building EList...")
      values$dda_elist                <- result$elist
      values$dda_n_proteins_prefilter <- result$n_prefilter
      values$dda_n_proteins_postfilter <- result$n_postfilter

      # Store as y_protein for downstream DE/QC/viz modules
      values$y_protein <- result$elist

      setProgress(0.7, detail = "Running limma DE...")

      # Run limma DE pipeline (same as DIA path in server_data.R)
      tryCatch({
        groups <- factor(values$metadata$Group)
        group_sizes <- table(groups)
        has_replicates <- all(group_sizes >= 2)

        if (has_replicates) {
          design <- model.matrix(~ 0 + groups)
          colnames(design) <- levels(groups)

          # Standard limma pipeline (not dpcDE -- MaxLFQ is already protein-level)
          fit <- limma::lmFit(result$elist, design)

          # Generate all pairwise contrasts
          combs <- combn(levels(groups), 2)
          forms <- apply(combs, 2, function(x) paste(x[2], "-", x[1]))
          contrast_matrix <- limma::makeContrasts(contrasts = forms, levels = design)
          fit <- limma::contrasts.fit(fit, contrast_matrix)
          fit <- limma::eBayes(fit)

          values$fit <- fit
          values$design <- design

          # Clear stale GSEA cache
          values$gsea_results_cache <- list()
          values$gsea_last_contrast <- NULL

          # Update all four comparison selectors
          updateSelectInput(session, "contrast_selector", choices = forms)
          updateSelectInput(session, "contrast_selector_signal", choices = forms, selected = forms[1])
          updateSelectInput(session, "contrast_selector_grid", choices = forms, selected = forms[1])
          updateSelectInput(session, "contrast_selector_pvalue", choices = forms, selected = forms[1])

          # Show DE Dashboard
          nav_show("main_tabs", "DE Dashboard")
          nav_show("main_tabs", "Gene Set Enrichment")
          nav_show("main_tabs", "AI Analysis")
          nav_show("main_tabs", "Output")
          nav_select("main_tabs", "DE Dashboard")

          n_de <- sum(limma::topTable(fit, coef = forms[1], number = Inf)$adj.P.Val < 0.05)
          showNotification(
            paste("DE analysis complete!", n_de, "significant proteins in", forms[1]),
            type = "message", duration = 10)
        } else {
          showNotification(
            "Some groups have <2 replicates. Skipping DE -- quantification-only mode.",
            type = "warning", duration = NULL)
        }
      }, error = function(e) {
        message("[DDA] limma DE failed: ", e$message)
        showNotification(paste("DE analysis failed:", e$message), type = "warning")
      })

      setProgress(1.0, detail = "Done!")
      showNotification(
        paste("DDA pipeline complete:",
              result$n_postfilter, "proteins after filtering"),
        type = "message", duration = 8)

      add_to_log("DDA Pipeline", c(
        sprintf("# Sage results: %d PSMs, %d proteins", nrow(values$dda_sage_psms), result$n_prefilter),
        sprintf("# Normalization: %s", input$dda_norm_method %||% "cyclicloess"),
        sprintf("# Imputation: %s", input$dda_impute_method %||% "perseus"),
        sprintf("# Valid value filter: %.0f%%", (input$dda_min_valid %||% 0.5) * 100),
        sprintf("# After filter: %d proteins", result$n_postfilter)
      ))
    })
  })

  # ============================================================================
  #    Group assignment for DDA samples
  # ============================================================================
  output$dda_group_assignment_ui <- renderUI({
    req(values$dda_lfq_wide)
    samples <- colnames(values$dda_lfq_wide)
    if (length(samples) == 0) return(NULL)

    tags$div(
      style = "background: #f8f9fa; padding: 12px; border-radius: 8px; margin-bottom: 12px;",
      tags$h6(icon("users"), " Assign Sample Groups"),
      tags$p(style = "font-size: 12px; color: #6c757d;",
        "Assign each sample to a group for differential expression analysis."),
      lapply(seq_along(samples), function(i) {
        div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 4px;",
          tags$span(style = "font-size: 12px; min-width: 200px; font-family: monospace;",
            samples[i]),
          textInput(
            paste0("dda_group_", i),
            label = NULL,
            value = "",
            width = "150px",
            placeholder = "Group name"
          )
        )
      }),
      actionButton("dda_apply_groups", "Apply Groups",
        icon = icon("check"), class = "btn-primary btn-sm mt-2"),
      actionButton("run_dda_pipeline", "Run DE Pipeline",
        icon = icon("play"), class = "btn-success btn-sm mt-2 ms-2")
    )
  })

  # Apply group assignments
 observeEvent(input$dda_apply_groups, {
    req(values$dda_lfq_wide)
    samples <- colnames(values$dda_lfq_wide)
    groups <- vapply(seq_along(samples), function(i) {
      input[[paste0("dda_group_", i)]] %||% ""
    }, character(1))

    if (any(!nzchar(groups))) {
      showNotification("Please assign all samples to a group.", type = "warning")
      return()
    }

    values$metadata <- data.frame(
      SampleID = samples,
      File.Name = samples,
      Group = groups,
      stringsAsFactors = FALSE
    )
    showNotification("Groups assigned!", type = "message")
  })

  # ============================================================================
  #    Status UI
  # ============================================================================
  output$dda_job_status_ui <- renderUI({
    status <- values$dda_status %||% "idle"

    switch(status,
      "idle" = NULL,
      "running" = {
        job_id <- values$dda_job_id %||% "?"
        div(
          class = "alert alert-info",
          style = "margin-top: 12px;",
          icon("spinner", class = "fa-spin"),
          paste(" Sage search running... Job ID:", job_id),
          tags$br(),
          tags$small("Polling every 15 seconds. You can navigate to other tabs.")
        )
      },
      "loading" = {
        div(
          class = "alert alert-info",
          style = "margin-top: 12px;",
          icon("spinner", class = "fa-spin"),
          " Loading Sage results from HPC..."
        )
      },
      "done" = {
        qc <- values$dda_qc_metrics
        div(
          class = "alert alert-success",
          style = "margin-top: 12px;",
          icon("check-circle"),
          " Sage search complete!",
          if (!is.null(qc)) {
            tags$div(
              style = "margin-top: 8px; font-size: 13px;",
              tags$strong(format(qc$n_psms, big.mark = ",")), " PSMs | ",
              tags$strong(format(qc$n_peptides, big.mark = ",")), " peptides | ",
              tags$strong(format(qc$n_proteins, big.mark = ",")), " proteins"
            )
          }
        )
      },
      "error" = {
        div(
          class = "alert alert-danger",
          style = "margin-top: 12px;",
          icon("exclamation-triangle"),
          " Sage search failed. Check SLURM logs for details.",
          if (!is.null(values$dda_output_dir)) {
            tags$small(style = "display: block; margin-top: 4px;",
              paste("Log dir:", file.path(values$dda_output_dir, "logs/")))
          }
        )
      }
    )
  })

  # ============================================================================
  #    Results summary UI (after pipeline)
  # ============================================================================
  output$dda_results_summary_ui <- renderUI({
    req(values$dda_elist)
    elist <- values$dda_elist
    n_pre  <- values$dda_n_proteins_prefilter %||% nrow(values$dda_lfq_wide)
    n_post <- values$dda_n_proteins_postfilter %||% nrow(elist$E)

    div(
      class = "alert alert-success",
      style = "margin-top: 12px;",
      icon("chart-bar"),
      tags$strong(" DDA pipeline complete"),
      tags$div(
        style = "margin-top: 8px; font-size: 13px;",
        tags$strong(format(n_post, big.mark = ",")), " proteins (from ",
        format(n_pre, big.mark = ","), " before filtering) | ",
        tags$strong(ncol(elist$E)), " samples",
        tags$br(),
        tags$small(
          "Normalization: ", values$dda_search_params$normalization %||% "cyclicloess",
          " | Imputation: ", values$dda_search_params$imputation %||% "perseus"
        )
      )
    )
  })

  # ============================================================================
  #    DDA QC Summary Card (rendered in QC tab or DDA panel)
  # ============================================================================
  output$dda_qc_summary_card <- renderUI({
    req(values$acquisition_mode == "dda")
    qc <- values$dda_qc_metrics
    if (is.null(qc)) return(NULL)

    div(
      style = paste(
        "background: white; border: 1px solid #dee2e6; border-radius: 8px;",
        "padding: 16px; margin-bottom: 16px;"
      ),
      tags$h6(icon("magnifying-glass"), " Sage DDA Search Summary",
        style = "margin-bottom: 12px; color: #2c3e50;"),
      div(
        class = "row text-center",
        div(class = "col-2",
          tags$h5(format(qc$n_psms, big.mark = ","), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "PSMs")
        ),
        div(class = "col-2",
          tags$h5(format(qc$n_peptides, big.mark = ","), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "Peptides")
        ),
        div(class = "col-2",
          tags$h5(format(qc$n_proteins, big.mark = ","), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "Proteins")
        ),
        div(class = "col-2",
          tags$h5(sprintf("%.1f", qc$med_pep_per_prot), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "Med. pep/prot")
        ),
        div(class = "col-2",
          tags$h5(paste0(qc$pct_missed_cleavage, "%"), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "Missed cleav.")
        ),
        div(class = "col-2",
          tags$h5(
            if (!is.na(qc$mass_error_ppm)) paste0(qc$mass_error_ppm, " ppm") else "N/A",
            style = "margin: 0; color: #2c3e50;"
          ),
          tags$small(style = "color: #6c757d;", "Mass error")
        )
      )
    )
  })

}
