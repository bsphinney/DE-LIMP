# ==============================================================================
#  SERVER MODULE — De Novo Sequencing (Cascadia SSL Integration, DIAMOND BLAST)
#  Called from app.R as: server_denovo(input, output, session, values, add_to_log)
# ==============================================================================

server_denovo <- function(input, output, session, values, add_to_log) {

  # --- Local reactive: SSH config mirror (same pattern as server_search.R) ---
  ssh_config <- reactive({
    if (!isTRUE(values$ssh_connected)) return(NULL)
    list(
      host  = values$ssh_host  %||% "",
      user  = values$ssh_user  %||% "",
      port  = values$ssh_port  %||% 22,
      key_path = values$ssh_key_path %||% "",
      modules  = values$ssh_modules %||% ""
    )
  })

  # ============================================================================
  #  1. SSL File Upload Handler
  # ============================================================================


  observeEvent(input$ssl_files, {
    req(input$ssl_files)

    add_to_log("De novo: uploading SSL files", "denovo")

    tryCatch({
      ssl_paths <- input$ssl_files$datapath
      ssl_names <- input$ssl_files$name

      # Validate file extensions
      exts <- tolower(tools::file_ext(ssl_names))
      if (!all(exts == "ssl")) {
        showNotification("Only .ssl files are accepted.", type = "error")
        return()
      }

      threshold <- input$denovo_score_threshold %||% 0.8

      # Parse all SSL files
      ssl_data <- parse_cascadia_ssl(ssl_paths, score_threshold = threshold)

      if (is.null(ssl_data) || nrow(ssl_data) == 0) {
        showNotification(
          "No de novo predictions above the score threshold.",
          type = "warning"
        )
        return()
      }

      # Store parsed data
      values$denovo_ssl_paths <- ssl_names
      values$denovo_data <- ssl_data
      values$denovo_score_threshold <- threshold

      # Classify against DIA-NN results if available
      classify_and_store(ssl_data)

      n_total <- nrow(ssl_data)
      showNotification(
        sprintf("Loaded %s de novo predictions from %d SSL file(s).",
                format(n_total, big.mark = ","), length(ssl_names)),
        type = "message"
      )
      add_to_log(
        sprintf("De novo: parsed %d predictions from %d files (threshold=%.2f)",
                n_total, length(ssl_names), threshold),
        "denovo"
      )

    }, error = function(e) {
      showNotification(
        paste("Error parsing SSL files:", conditionMessage(e)),
        type = "error"
      )
      add_to_log(paste("De novo SSL parse error:", conditionMessage(e)), "error")
    })
  })

  # ============================================================================
  #  2. Score Threshold Slider — Re-filter and Reclassify
  # ============================================================================

  observeEvent(input$denovo_score_threshold, {
    req(values$denovo_data)

    new_threshold <- input$denovo_score_threshold
    if (identical(new_threshold, values$denovo_score_threshold)) return()

    values$denovo_score_threshold <- new_threshold

    # Re-parse from stored full data (re-read if we have paths, or filter in memory)
    # SSL files store the full parsed data; we filter in memory for speed
    full_data <- values$denovo_data

    # If the stored data was already filtered at a higher threshold,
    # we need to re-parse from files. But since fileInput paths are transient,
    # we keep the full unfiltered data in a separate slot on first load.
    if (!is.null(values$denovo_data_full)) {
      filtered <- values$denovo_data_full[values$denovo_data_full$score >= new_threshold, ]
    } else {
      # Fallback: can only filter down from current data
      filtered <- full_data[full_data$score >= new_threshold, ]
    }

    values$denovo_data <- filtered
    classify_and_store(filtered)

    add_to_log(
      sprintf("De novo: threshold adjusted to %.2f (%d predictions retained)",
              new_threshold, nrow(filtered)),
      "denovo"
    )
  }, ignoreInit = TRUE)

  # ============================================================================
  #  Helper: Classify de novo peptides and store results
  # ============================================================================

  classify_and_store <- function(ssl_data) {
    if (is.null(values$raw_data)) {
      # No DIA-NN data loaded yet — store as all unclassified
      values$denovo_classified <- list(
        confirmed = ssl_data[0, ],
        novel     = ssl_data,
        protein_summary = data.frame(
          Protein.Group = character(0),
          n_denovo_confirmed = integer(0),
          denovo_max_score = numeric(0),
          stringsAsFactors = FALSE
        )
      )
      values$denovo_protein_summary <- values$denovo_classified$protein_summary
      return()
    }

    # Build DIA-NN peptide reference from raw_data
    diann_report <- tryCatch({
      if (!is.null(values$raw_data$genes) &&
          "Stripped.Sequence" %in% names(values$raw_data$genes)) {
        values$raw_data$genes
      } else {
        NULL
      }
    }, error = function(e) NULL)

    if (is.null(diann_report)) {
      # Cannot cross-reference — treat all as novel
      values$denovo_classified <- list(
        confirmed = ssl_data[0, ],
        novel     = ssl_data,
        protein_summary = data.frame(
          Protein.Group = character(0),
          n_denovo_confirmed = integer(0),
          denovo_max_score = numeric(0),
          stringsAsFactors = FALSE
        )
      )
      values$denovo_protein_summary <- values$denovo_classified$protein_summary
      return()
    }

    # Cross-reference using classify_denovo_peptides from helpers_denovo.R
    classified <- classify_denovo_peptides(ssl_data, diann_report)
    values$denovo_classified <- classified
    values$denovo_protein_summary <- classified$protein_summary
  }

  # ============================================================================
  #  3. Summary Cards
  # ============================================================================

  output$denovo_summary_cards <- renderUI({
    req(values$denovo_data)

    n_total <- nrow(values$denovo_data)
    classified <- values$denovo_classified

    n_confirmed <- if (!is.null(classified$confirmed)) nrow(classified$confirmed) else 0
    n_novel <- if (!is.null(classified$novel)) nrow(classified$novel) else 0

    confirmation_rate <- if (n_total > 0) {
      round(100 * n_confirmed / n_total, 1)
    } else {
      0
    }

    n_proteins_confirmed <- if (!is.null(values$denovo_protein_summary)) {
      nrow(values$denovo_protein_summary)
    } else {
      0
    }

    tags$div(
      class = "row",
      style = "margin-bottom: 15px;",

      tags$div(
        class = "col-md-3",
        tags$div(
          class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #3498db; padding: 15px;",
          tags$h4(format(n_total, big.mark = ","), style = "margin: 0; color: #3498db;"),
          tags$small("Total Predictions")
        )
      ),

      tags$div(
        class = "col-md-3",
        tags$div(
          class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #2ecc71; padding: 15px;",
          tags$h4(format(n_confirmed, big.mark = ","), style = "margin: 0; color: #2ecc71;"),
          tags$small("Confirmed Peptides")
        )
      ),

      tags$div(
        class = "col-md-3",
        tags$div(
          class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #e67e22; padding: 15px;",
          tags$h4(format(n_novel, big.mark = ","), style = "margin: 0; color: #e67e22;"),
          tags$small("Novel Peptides")
        )
      ),

      tags$div(
        class = "col-md-3",
        tags$div(
          class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #9b59b6; padding: 15px;",
          tags$h4(paste0(confirmation_rate, "%"), style = "margin: 0; color: #9b59b6;"),
          tags$small(paste0("Confirmation Rate (", n_proteins_confirmed, " proteins)"))
        )
      )
    )
  })

  # ============================================================================
  #  4. Confirmed Peptides Table
  # ============================================================================

  output$denovo_confirmed_table <- DT::renderDT({
    req(values$denovo_classified)
    confirmed <- values$denovo_classified$confirmed
    req(nrow(confirmed) > 0)

    # Build display table
    display_df <- data.frame(
      Sequence       = confirmed$sequence,
      Stripped       = confirmed$seq_stripped,
      Score          = round(confirmed$score, 3),
      Charge         = confirmed$charge,
      RT_min         = if ("retention_time" %in% names(confirmed)) {
        round(confirmed$retention_time, 2)
      } else {
        NA_real_
      },
      Protein.Group  = if ("Protein.Group" %in% names(confirmed)) {
        confirmed$Protein.Group
      } else {
        NA_character_
      },
      Source_File     = confirmed$source_file,
      stringsAsFactors = FALSE
    )

    DT::datatable(
      display_df,
      rownames = FALSE,
      filter   = "top",
      selection = "multiple",
      options  = list(
        pageLength = 25,
        scrollX    = TRUE,
        order      = list(list(2, "desc")),  # Sort by Score descending
        dom        = "Bfrtip",
        buttons    = list("csv", "excel")
      ),
      extensions = "Buttons",
      caption = "Confirmed: de novo peptides matching DIA-NN database search results (I/L normalized)"
    )
  })

  # ============================================================================
  #  5. Novel Peptides Table
  # ============================================================================

  output$denovo_novel_table <- DT::renderDT({
    req(values$denovo_classified)
    novel <- values$denovo_classified$novel
    req(nrow(novel) > 0)

    # Base columns
    display_df <- data.frame(
      Sequence    = novel$sequence,
      Stripped    = novel$seq_stripped,
      Score       = round(novel$score, 3),
      Charge      = novel$charge,
      RT_min      = if ("retention_time" %in% names(novel)) {
        round(novel$retention_time, 2)
      } else {
        NA_real_
      },
      Source_File  = novel$source_file,
      stringsAsFactors = FALSE
    )

    # Append DIAMOND BLAST results if available
    blast <- values$denovo_novel_blast
    if (!is.null(blast) && nrow(blast) > 0) {
      # Join on stripped sequence (I/L normalized)
      blast_lookup <- blast[!duplicated(blast$peptide_sequence), ]
      blast_map <- stats::setNames(blast_lookup$subject, blast_lookup$peptide_sequence)
      identity_map <- stats::setNames(blast_lookup$identity, blast_lookup$peptide_sequence)
      evalue_map <- stats::setNames(blast_lookup$evalue, blast_lookup$peptide_sequence)

      display_df$BLAST_Hit     <- blast_map[novel$seq_stripped]
      display_df$Identity_Pct  <- round(identity_map[novel$seq_stripped], 1)
      display_df$E_Value       <- evalue_map[novel$seq_stripped]
    }

    DT::datatable(
      display_df,
      rownames = FALSE,
      filter   = "top",
      selection = "multiple",
      options  = list(
        pageLength = 25,
        scrollX    = TRUE,
        order      = list(list(2, "desc")),
        dom        = "Bfrtip",
        buttons    = list("csv", "excel")
      ),
      extensions = "Buttons",
      caption = htmltools::tags$caption(
        style = "caption-side: top; color: #e67e22; font-weight: bold;",
        "Novel: de novo peptides NOT found in DIA-NN results.",
        tags$br(),
        tags$small(
          style = "color: #666; font-weight: normal;",
          "These may represent sequence variants, unexpected organisms, or proteins absent from your reference FASTA."
        )
      )
    )
  })

  # ============================================================================
  #  6. Run DIAMOND BLAST Button
  # ============================================================================

  observeEvent(input$run_diamond_blast, {
    req(values$denovo_classified)
    req(nrow(values$denovo_classified$novel) > 0)

    novel_peptides <- unique(values$denovo_classified$novel$seq_stripped)

    if (length(novel_peptides) == 0) {
      showNotification("No novel peptides to BLAST.", type = "warning")
      return()
    }

    # Determine FASTA path
    fasta_path <- values$diann_fasta_files[1] %||% values$diann_fasta_path %||% input$diamond_fasta_path
    if (is.null(fasta_path) || !nzchar(fasta_path)) {
      showNotification(
        "No FASTA file specified. Set the reference FASTA for DIAMOND BLAST.",
        type = "error"
      )
      return()
    }

    cfg <- ssh_config()

    # --- SSH/HPC path ---
    if (!is.null(cfg) && isTRUE(values$ssh_connected)) {
      tryCatch({
        withProgress(message = "Running DIAMOND BLAST on HPC...", value = 0.1, {

          output_dir <- values$diann_output_dir %||% paste0("/tmp/delimp_denovo_", Sys.getpid())
          denovo_dir <- file.path(output_dir, "denovo")

          # Create remote directory
          ssh_exec(cfg, paste("mkdir -p", shQuote(denovo_dir)), timeout = 15)
          setProgress(0.2, detail = "Created remote directory")

          # Write query FASTA locally, then SCP upload
          query_fasta_local <- tempfile(fileext = ".fasta")
          query_lines <- paste0(">denovo_", seq_along(novel_peptides), "\n", novel_peptides)
          writeLines(query_lines, query_fasta_local)

          query_fasta_remote <- file.path(denovo_dir, "novel_denovo_queries.fasta")
          scp_upload(cfg, query_fasta_local, query_fasta_remote)
          setProgress(0.3, detail = "Uploaded query FASTA")

          # Build DIAMOND DB if not cached
          diamond_bin <- "/cvmfs/hpc.ucdavis.edu/sw/spack/environments/main/view/generic/diamond-2.1.7/bin/diamond"
          diamond_db_remote <- file.path(denovo_dir, "ref_diamond.dmnd")

          # Check if DB already exists (skip rebuild)
          db_check <- ssh_exec(cfg, paste("test -f", shQuote(diamond_db_remote), "&& echo EXISTS"), timeout = 10)
          if (!any(grepl("EXISTS", db_check$stdout))) {
            db_build_cmd <- paste(
              diamond_bin, "makedb",
              "--in", shQuote(fasta_path),
              "--db", shQuote(diamond_db_remote),
              "--threads 4 --quiet"
            )
            ssh_exec(cfg, db_build_cmd, timeout = 300)
          }
          setProgress(0.5, detail = "DIAMOND database ready")

          # Run DIAMOND blastp
          blast_out_remote <- file.path(denovo_dir, "novel_denovo_blast.tsv")
          blast_cmd <- paste(
            diamond_bin, "blastp",
            "--query", shQuote(query_fasta_remote),
            "--db", shQuote(diamond_db_remote),
            "--out", shQuote(blast_out_remote),
            "--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore",
            "--id", input$diamond_min_identity %||% 90,
            "--threads 4 --sensitive --quiet"
          )
          ssh_exec(cfg, blast_cmd, timeout = 600)
          setProgress(0.8, detail = "BLAST complete, downloading results")

          # Download results
          blast_out_local <- tempfile(fileext = ".tsv")
          scp_download(cfg, blast_out_remote, blast_out_local)

          if (file.exists(blast_out_local) && file.size(blast_out_local) > 0) {
            hits <- data.table::fread(blast_out_local, header = FALSE)
            names(hits) <- c("query_idx", "subject", "identity", "length",
                             "qlen", "slen", "evalue", "bitscore")

            # Extract numeric index from query name (denovo_1, denovo_2, ...)
            hit_idx <- as.integer(gsub("denovo_", "", hits$query_idx))
            hits$peptide_sequence <- novel_peptides[hit_idx]

            # Extract protein accession
            hits$protein <- stringr::str_extract(hits$subject, "(?<=\\|)[^|]+(?=\\|)")
            if (all(is.na(hits$protein))) {
              hits$protein <- hits$subject
            }

            values$denovo_novel_blast <- as.data.frame(hits)
            n_hits <- length(unique(hits$peptide_sequence))
            showNotification(
              sprintf("DIAMOND BLAST: %d novel peptides mapped to %d protein hits.",
                      n_hits, nrow(hits)),
              type = "message"
            )
            add_to_log(
              sprintf("De novo DIAMOND BLAST: %d peptides -> %d hits (SSH)", n_hits, nrow(hits)),
              "denovo"
            )
          } else {
            values$denovo_novel_blast <- data.frame()
            showNotification("No DIAMOND BLAST hits found.", type = "warning")
            add_to_log("De novo DIAMOND BLAST: no hits found", "denovo")
          }

          setProgress(1.0, detail = "Done")
        })

      }, error = function(e) {
        showNotification(
          paste("DIAMOND BLAST error:", conditionMessage(e)),
          type = "error"
        )
        add_to_log(paste("De novo DIAMOND BLAST error:", conditionMessage(e)), "error")
      })

    } else {
      # --- Local path ---
      tryCatch({
        # Check if diamond is available locally
        diamond_check <- tryCatch(
          system2("diamond", args = "version", stdout = TRUE, stderr = TRUE),
          error = function(e) NULL
        )

        if (is.null(diamond_check)) {
          showNotification(
            "DIAMOND not found locally. Connect via SSH for HPC execution, or install DIAMOND.",
            type = "error"
          )
          return()
        }

        withProgress(message = "Running DIAMOND BLAST locally...", value = 0.1, {

          output_dir <- values$diann_output_dir %||% tempdir()
          denovo_dir <- file.path(output_dir, "denovo")
          dir.create(denovo_dir, recursive = TRUE, showWarnings = FALSE)

          blast_results <- run_diamond_blast(
            novel_peptides = novel_peptides,
            fasta_path     = fasta_path,
            diamond_db     = NULL,
            output_dir     = denovo_dir,
            min_identity   = input$diamond_min_identity %||% 90,
            threads        = 4
          )

          setProgress(0.9, detail = "Processing results")

          values$denovo_novel_blast <- blast_results

          if (nrow(blast_results) > 0) {
            n_hits <- length(unique(blast_results$peptide_sequence))
            showNotification(
              sprintf("DIAMOND BLAST: %d novel peptides mapped.", n_hits),
              type = "message"
            )
            add_to_log(
              sprintf("De novo DIAMOND BLAST: %d peptides -> %d hits (local)",
                      n_hits, nrow(blast_results)),
              "denovo"
            )
          } else {
            showNotification("No DIAMOND BLAST hits found.", type = "warning")
            add_to_log("De novo DIAMOND BLAST: no hits found (local)", "denovo")
          }

          setProgress(1.0, detail = "Done")
        })

      }, error = function(e) {
        showNotification(
          paste("Local DIAMOND error:", conditionMessage(e)),
          type = "error"
        )
        add_to_log(paste("De novo DIAMOND error (local):", conditionMessage(e)), "error")
      })
    }
  })

  # ============================================================================
  #  7. Cascadia Job Submission (conditional on SSH + checkbox)
  # ============================================================================

  observeEvent(input$submit_cascadia_job, {
    req(isTRUE(values$ssh_connected))
    req(values$diann_raw_files)
    req(input$cascadia_model_path)

    cfg <- ssh_config()
    req(cfg)

    tryCatch({
      withProgress(message = "Submitting Cascadia de novo job...", value = 0.1, {

        analysis_name <- values$diann_analysis_name %||% "cascadia_run"
        # Sanitize for SLURM job name (same pattern as DIA-NN)
        safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", analysis_name)

        raw_files <- values$diann_raw_files$full_path
        output_dir <- values$diann_output_dir
        req(output_dir)

        denovo_dir <- file.path(output_dir, "denovo")
        model_ckpt <- input$cascadia_model_path
        conda_env  <- input$cascadia_conda_env %||% "cascadia"
        min_score  <- input$cascadia_min_score %||% 0.8
        partition   <- input$cascadia_partition %||% "gpu"
        account     <- input$diann_account %||% "genome-center-grp"
        gpu_type    <- input$cascadia_gpu_type %||% "1"
        time_hours  <- input$cascadia_time_hours %||% 4

        # Generate sbatch script using helper from helpers_denovo.R
        sbatch_content <- generate_cascadia_sbatch(
          analysis_name = safe_name,
          raw_files     = raw_files,
          output_dir    = output_dir,
          model_ckpt    = model_ckpt,
          conda_env     = conda_env,
          min_score     = min_score,
          partition     = partition,
          account       = account,
          gpu_type      = gpu_type,
          time_hours    = time_hours
        )

        setProgress(0.3, detail = "Uploading sbatch script")

        # Write locally and SCP upload
        local_script <- tempfile(fileext = ".sbatch")
        writeLines(sbatch_content, local_script)

        remote_script <- file.path(output_dir, "run_cascadia.sbatch")

        # Create remote dirs
        ssh_exec(cfg, paste("mkdir -p", shQuote(denovo_dir), shQuote(file.path(output_dir, "logs"))),
                 timeout = 15)
        scp_upload(cfg, local_script, remote_script)

        setProgress(0.5, detail = "Submitting to SLURM")

        # Submit via sbatch
        sbatch_bin <- values$ssh_sbatch_path %||% "sbatch"
        submit_result <- ssh_exec(cfg, paste(sbatch_bin, shQuote(remote_script)),
                                  timeout = 30)

        # Parse job ID
        job_id <- trimws(stringr::str_extract(submit_result$stdout, "\\d+"))

        if (is.na(job_id) || !nzchar(job_id)) {
          showNotification(
            paste("Cascadia submission failed:", submit_result$stderr),
            type = "error"
          )
          add_to_log(
            paste("Cascadia submit failed:", submit_result$stderr),
            "error"
          )
          return()
        }

        setProgress(0.8, detail = paste("Job ID:", job_id))

        # Store job info
        values$denovo_job_id <- job_id
        values$denovo_job_status <- "queued"
        values$cascadia_model_ckpt <- model_ckpt

        showNotification(
          sprintf("Cascadia job submitted: %s (partition=%s)", job_id, partition),
          type = "message"
        )
        add_to_log(
          sprintf("Cascadia job submitted: ID=%s, partition=%s, files=%d, model=%s",
                  job_id, partition, length(raw_files), basename(model_ckpt)),
          "denovo"
        )

        setProgress(1.0, detail = "Submitted")
      })

    }, error = function(e) {
      showNotification(
        paste("Cascadia submission error:", conditionMessage(e)),
        type = "error"
      )
      add_to_log(paste("Cascadia submit error:", conditionMessage(e)), "error")
    })
  })

  # --- Cascadia job status monitor ---
  observe({
    req(values$denovo_job_id)
    req(values$denovo_job_status %in% c("queued", "running"))
    req(isTRUE(values$ssh_connected))

    cfg <- isolate(ssh_config())
    req(cfg)

    # Poll every 30 seconds
    invalidateLater(30000, session)

    tryCatch({
      sacct_bin <- values$ssh_sacct_path %||% "sacct"
      result <- ssh_exec(
        cfg,
        paste(sacct_bin, "-j", values$denovo_job_id,
              "--format=JobID,State --noheader --parsable2"),
        timeout = 15
      )

      if (!is.null(result$stdout) && nzchar(result$stdout)) {
        lines <- strsplit(trimws(result$stdout), "\n")[[1]]
        # Filter out .extern/.batch substeps
        main_lines <- lines[!grepl("\\.", lines)]

        if (length(main_lines) > 0) {
          state <- trimws(strsplit(main_lines[1], "\\|")[[1]][2])

          if (state %in% c("COMPLETED")) {
            values$denovo_job_status <- "completed"
            showNotification("Cascadia job completed! Upload SSL files to continue.",
                             type = "message", duration = 10)
            add_to_log(
              sprintf("Cascadia job %s completed", values$denovo_job_id),
              "denovo"
            )
          } else if (state %in% c("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY")) {
            values$denovo_job_status <- "failed"
            showNotification(
              sprintf("Cascadia job %s: %s", values$denovo_job_id, state),
              type = "error", duration = 10
            )
            add_to_log(
              sprintf("Cascadia job %s failed: %s", values$denovo_job_id, state),
              "error"
            )
          } else if (state %in% c("RUNNING")) {
            values$denovo_job_status <- "running"
          }
          # PENDING stays as "queued"
        }
      }
    }, error = function(e) {
      # Silently ignore polling errors (network glitch, etc.)
    })
  })

  # ============================================================================
  #  Score Distribution Plot (plotly for bslib safety)
  # ============================================================================

  output$denovo_score_plot <- plotly::renderPlotly({
    req(values$denovo_data)
    req(nrow(values$denovo_data) > 0)

    df <- values$denovo_data
    threshold <- values$denovo_score_threshold %||% 0.8

    # Assign classification for color
    match_type <- if (!is.null(values$denovo_classified)) {
      confirmed_seqs <- values$denovo_classified$confirmed$seq_norm
      ifelse(df$seq_norm %in% confirmed_seqs, "Confirmed", "Novel")
    } else {
      rep("Unclassified", nrow(df))
    }

    plot_df <- data.frame(
      score = df$score,
      type  = match_type,
      stringsAsFactors = FALSE
    )

    colors <- c("Confirmed" = "#2ecc71", "Novel" = "#e67e22", "Unclassified" = "#95a5a6")

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = score, fill = type)) +
      ggplot2::geom_histogram(bins = 50, alpha = 0.8, position = "stack") +
      ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "#e74c3c") +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(
        x = "Cascadia Confidence Score",
        y = "Count",
        fill = "Classification"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")

    plotly::ggplotly(p) %>%
      plotly::layout(
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.05)
      )
  })

  # ============================================================================
  #  Cascadia job status badge output (for UI conditionalPanel)
  # ============================================================================

  output$denovo_job_status_text <- renderUI({
    status <- values$denovo_job_status %||% "none"
    job_id <- values$denovo_job_id

    if (status == "none") return(NULL)

    badge_color <- switch(status,
      queued    = "#f39c12",
      running   = "#3498db",
      completed = "#2ecc71",
      failed    = "#e74c3c",
      "#95a5a6"
    )

    badge_label <- switch(status,
      queued    = "QUEUED",
      running   = "RUNNING",
      completed = "COMPLETED",
      failed    = "FAILED",
      toupper(status)
    )

    tags$div(
      style = "margin-top: 8px;",
      tags$span(
        class = "badge",
        style = sprintf("background: %s; color: white; padding: 4px 10px; font-size: 12px;",
                        badge_color),
        paste("Cascadia:", badge_label)
      ),
      if (!is.null(job_id)) {
        tags$small(style = "margin-left: 6px; color: #888;", paste("Job:", job_id))
      }
    )
  })

  # ============================================================================
  #  8. AI Context Injection Helper
  # ============================================================================

  # Returns a character string summarizing de novo results for Gemini/Claude prompts.
  # Called from server_ai.R when building the AI context block.
  build_denovo_ai_context <- function() {
    if (is.null(values$denovo_classified)) return("")

    classified <- values$denovo_classified
    n_total    <- nrow(values$denovo_data %||% data.frame())
    n_confirmed <- nrow(classified$confirmed)
    n_novel     <- nrow(classified$novel)
    n_proteins  <- nrow(values$denovo_protein_summary %||% data.frame())
    threshold   <- values$denovo_score_threshold %||% 0.8

    context <- paste0(
      "DE NOVO SEQUENCING CONTEXT (Cascadia):\n",
      sprintf("- %s total de novo peptide predictions (score >= %.2f)\n",
              format(n_total, big.mark = ","), threshold),
      sprintf("- %s confirmed (matched DIA-NN database search, I/L normalized)\n",
              format(n_confirmed, big.mark = ",")),
      sprintf("- %s novel (not in DIA-NN results)\n",
              format(n_novel, big.mark = ",")),
      sprintf("- %d proteins have de novo-confirmed peptides\n", n_proteins)
    )

    # Add DIAMOND BLAST summary if available
    blast <- values$denovo_novel_blast
    if (!is.null(blast) && nrow(blast) > 0) {
      n_mapped <- length(unique(blast$peptide_sequence))
      n_unique_proteins <- length(unique(blast$protein))
      median_identity <- round(median(blast$identity, na.rm = TRUE), 1)
      context <- paste0(
        context,
        sprintf("- DIAMOND BLAST: %d novel peptides mapped to %d proteins (median identity: %.1f%%)\n",
                n_mapped, n_unique_proteins, median_identity)
      )
    }

    # Per-protein detail (top confirmed proteins)
    if (n_proteins > 0) {
      top_proteins <- head(
        values$denovo_protein_summary[
          order(-values$denovo_protein_summary$n_denovo_confirmed), ],
        10
      )
      context <- paste0(
        context,
        "\nTop de novo-confirmed proteins:\n",
        paste0(
          sprintf("  %s: %d confirmed peptides (best score: %.2f)",
                  top_proteins$Protein.Group,
                  top_proteins$n_denovo_confirmed,
                  top_proteins$denovo_max_score),
          collapse = "\n"
        ),
        "\n"
      )
    }

    context <- paste0(
      context,
      "\nProteins marked [De Novo Confirmed] have orthogonal sequence-level validation ",
      "from de novo sequencing, increasing confidence in their identification. ",
      "Novel peptides not in the FASTA may represent sequence variants, unexpected organisms, ",
      "or proteoforms missed in the original database.\n"
    )

    return(context)
  }

  # Expose the helper to other modules via values
  values$build_denovo_ai_context <- build_denovo_ai_context

  # ============================================================================
  #  Re-classify when DIA-NN data becomes available
  # ============================================================================

  observeEvent(values$raw_data, {
    if (!is.null(values$denovo_data) && nrow(values$denovo_data) > 0) {
      classify_and_store(values$denovo_data)
      add_to_log("De novo: reclassified against updated DIA-NN data", "denovo")
    }
  }, ignoreInit = TRUE)

  # ============================================================================
  #  Download handlers
  # ============================================================================

  output$download_confirmed_csv <- downloadHandler(
    filename = function() {
      paste0("denovo_confirmed_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$denovo_classified$confirmed)
      utils::write.csv(values$denovo_classified$confirmed, file, row.names = FALSE)
    }
  )

  output$download_novel_csv <- downloadHandler(
    filename = function() {
      paste0("denovo_novel_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$denovo_classified$novel)
      novel <- values$denovo_classified$novel

      # Append BLAST results if available
      blast <- values$denovo_novel_blast
      if (!is.null(blast) && nrow(blast) > 0) {
        blast_dedup <- blast[!duplicated(blast$peptide_sequence), ]
        novel <- merge(novel, blast_dedup[, c("peptide_sequence", "subject", "identity", "evalue")],
                       by.x = "seq_stripped", by.y = "peptide_sequence",
                       all.x = TRUE)
      }

      utils::write.csv(novel, file, row.names = FALSE)
    }
  )

  output$download_protein_summary_csv <- downloadHandler(
    filename = function() {
      paste0("denovo_protein_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$denovo_protein_summary)
      utils::write.csv(values$denovo_protein_summary, file, row.names = FALSE)
    }
  )

}
