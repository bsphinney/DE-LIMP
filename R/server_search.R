# ==============================================================================
#  server_search.R
#  DIA-NN HPC Search Integration — New Search tab server logic
#  Handles: file browsing, UniProt FASTA download, sbatch generation,
#  job submission, monitoring, auto-load, and job queue management.
# ==============================================================================

server_search <- function(input, output, session, values, add_to_log, hpc_mode) {

  # Early return if not HPC mode
  if (!hpc_mode) return(invisible())

  # ============================================================================
  #    SSH Config Reactive
  # ============================================================================

  ssh_config <- reactive({
    if (is.null(input$search_connection_mode) ||
        input$search_connection_mode != "ssh") return(NULL)
    list(
      host = input$ssh_host,
      user = input$ssh_user,
      port = input$ssh_port %||% 22,
      key_path = input$ssh_key_path,
      modules = input$ssh_modules %||% ""
    )
  })

  # ============================================================================
  #    Job Queue Persistence (survives app restarts)
  # ============================================================================

  job_queue_path <- file.path(Sys.getenv("HOME"), ".delimp_job_queue.rds")
  job_queue_loaded <- reactiveVal(FALSE)

  # Load saved jobs on startup
  observe({
    if (file.exists(job_queue_path)) {
      tryCatch({
        saved_jobs <- readRDS(job_queue_path)
        if (is.list(saved_jobs) && length(saved_jobs) > 0) {
          values$diann_jobs <- saved_jobs
          n_active <- sum(vapply(saved_jobs, function(j)
            j$status %in% c("queued", "running"), logical(1)))
          if (n_active > 0) {
            showNotification(
              sprintf("Restored %d job(s) from previous session (%d active).",
                      length(saved_jobs), n_active),
              type = "message", duration = 5)
          }
        }
      }, error = function(e) {
        message("[DE-LIMP] Failed to load saved job queue: ", e$message)
      })
    }
    job_queue_loaded(TRUE)
  }) |> bindEvent(TRUE)  # Run once on startup

  # Save jobs to disk whenever the queue changes (after initial load)
  observe({
    jobs <- values$diann_jobs
    req(job_queue_loaded())
    tryCatch({
      saveRDS(jobs, job_queue_path)
    }, error = function(e) {
      message("[DE-LIMP] Failed to save job queue: ", e$message)
    })
  })

  # ============================================================================
  #    SSH Connection Test
  # ============================================================================

  observeEvent(input$test_ssh_btn, {
    cfg <- ssh_config()
    if (is.null(cfg)) return()

    withProgress(message = "Testing SSH connection...", {
      result <- test_ssh_connection(cfg)
    })

    output$ssh_status_ui <- renderUI({
      if (result$success) {
        div(class = "alert alert-success py-1 px-2 mt-2",
          style = "font-size: 0.82em;",
          icon("check-circle"), " ", result$message)
      } else {
        div(class = "alert alert-danger py-1 px-2 mt-2",
          style = "font-size: 0.82em;",
          icon("times-circle"), " ", result$message)
      }
    })

    values$ssh_connected <- result$success
    values$ssh_sbatch_path <- result$sbatch_path
  })

  # ============================================================================
  #    SSH Remote File Scanning
  # ============================================================================

  observeEvent(input$ssh_scan_raw_btn, {
    cfg <- ssh_config()
    req(cfg, input$ssh_raw_data_dir, nzchar(input$ssh_raw_data_dir))

    withProgress(message = "Scanning remote directory...", {
      raw_files <- ssh_scan_raw_files(cfg, input$ssh_raw_data_dir)
    })

    if (nrow(raw_files) > 0) {
      # Add full_path column for sbatch script generation
      raw_files$full_path <- file.path(input$ssh_raw_data_dir, raw_files$filename)
    }
    values$diann_raw_files <- raw_files
  })

  observeEvent(input$ssh_scan_fasta_btn, {
    cfg <- ssh_config()
    req(cfg, input$ssh_fasta_browse_dir, nzchar(input$ssh_fasta_browse_dir))

    withProgress(message = "Scanning remote FASTA files...", {
      fasta_files <- ssh_scan_fasta_files(cfg, input$ssh_fasta_browse_dir)
    })

    values$diann_fasta_files <- as.character(fasta_files)

    output$browsed_fasta_info <- renderUI({
      if (length(fasta_files) == 0) {
        div(class = "text-muted small mt-2", "No FASTA files found in remote directory.")
      } else {
        div(class = "alert alert-success py-1 px-2 mt-2",
          style = "font-size: 0.82em;",
          sprintf("%d FASTA file(s): %s", length(fasta_files),
                  paste(names(fasta_files), collapse = ", ")))
      }
    })
  })

  # ============================================================================
  #    shinyFiles Initialization (local mode only)
  # ============================================================================

  volumes <- c(Home = Sys.getenv("HOME"), Root = "/")

  shinyFiles::shinyDirChoose(input, "raw_data_dir", roots = volumes, session = session)
  shinyFiles::shinyDirChoose(input, "fasta_browse_dir", roots = volumes, session = session)
  shinyFiles::shinyDirChoose(input, "output_base_dir", roots = volumes, session = session)
  shinyFiles::shinyFileChoose(input, "lib_file", roots = volumes, session = session,
    filetypes = c("speclib", "tsv", "csv"))

  # ============================================================================
  #    File Selection Observers
  # ============================================================================

  # Raw data directory selection
  observeEvent(input$raw_data_dir, {
    if (is.integer(input$raw_data_dir)) return()  # Initial NULL state

    dir_path <- shinyFiles::parseDirPath(volumes, input$raw_data_dir)
    if (length(dir_path) == 0 || !nzchar(dir_path)) return()

    raw_files <- scan_raw_files(as.character(dir_path))
    values$diann_raw_files <- raw_files
  })

  output$raw_file_summary <- renderUI({
    req(values$diann_raw_files)
    df <- values$diann_raw_files

    if (nrow(df) == 0) {
      return(div(class = "alert alert-warning",
        style = "margin-top: 8px; padding: 8px; font-size: 0.85em;",
        icon("exclamation-triangle"),
        " No .d / .raw / .mzML files found in selected directory."
      ))
    }

    n_files <- nrow(df)
    total_size <- sum(df$size_mb)
    types <- paste(unique(df$type), collapse = ", ")

    div(class = "alert alert-success",
      style = "margin-top: 8px; padding: 8px; font-size: 0.85em;",
      icon("check-circle"),
      sprintf(" %d files found (%s) — %.1f GB total", n_files, types, total_size / 1024)
    )
  })

  # Spectral library selection
  observeEvent(input$lib_file, {
    if (is.integer(input$lib_file)) return()

    file_info <- shinyFiles::parseFilePaths(volumes, input$lib_file)
    if (nrow(file_info) == 0) return()

    values$diann_speclib <- as.character(file_info$datapath)

    # Auto-switch to library mode if speclib selected
    updateRadioButtons(session, "search_mode", selected = "library")
  })

  # SSH mode: spectral library path from text input
  observeEvent(input$ssh_lib_file, {
    if (nzchar(input$ssh_lib_file %||% "")) {
      values$diann_speclib <- input$ssh_lib_file
      updateRadioButtons(session, "search_mode", selected = "library")
    }
  })

  output$lib_file_info <- renderUI({
    if (is.null(values$diann_speclib)) return(NULL)
    div(class = "alert alert-info",
      style = "margin-top: 8px; padding: 8px; font-size: 0.85em;",
      icon("book"), " ", basename(values$diann_speclib)
    )
  })

  # ============================================================================
  #    Phosphoproteomics Search Mode — auto-configure settings
  # ============================================================================

  observeEvent(input$search_mode, {
    if (input$search_mode == "phospho") {
      # Phospho-optimized DIA-NN settings
      updateNumericInput(session, "diann_max_var_mods", value = 3)
      updateCheckboxInput(session, "mod_met_ox", value = TRUE)
      updateTextAreaInput(session, "extra_var_mods",
        value = "UniMod:21,79.966331,STY")
      updateNumericInput(session, "diann_missed_cleavages", value = 2)
      showNotification(
        paste("Phospho mode: STY phosphorylation (UniMod:21) added,",
              "max var mods = 3, missed cleavages = 2"),
        type = "message", duration = 8)
    }
  }, ignoreInit = TRUE)

  # ============================================================================
  #    UniProt FASTA Download
  # ============================================================================

  observeEvent(input$search_uniprot, {
    req(nzchar(input$uniprot_search_query))

    withProgress(message = "Searching UniProt...", {
      results <- search_uniprot_proteomes(input$uniprot_search_query)
      values$uniprot_results <- results
    })

    if (nrow(values$uniprot_results) == 0) {
      showNotification("No proteomes found. Try a different search term.", type = "warning")
    }
  })

  output$uniprot_results_table <- DT::renderDT({
    req(values$uniprot_results, nrow(values$uniprot_results) > 0)

    display_df <- values$uniprot_results[, c("upid", "organism", "common_name", "protein_count")]
    colnames(display_df) <- c("ID", "Organism", "Common Name", "Proteins")

    DT::datatable(display_df,
      selection = "single",
      options = list(
        pageLength = 5, dom = "t", scrollY = "150px",
        columnDefs = list(list(width = "80px", targets = 0))
      ),
      rownames = FALSE,
      class = "compact"
    )
  })

  # FASTA filename preview
  output$fasta_filename_preview <- renderUI({
    req(values$uniprot_results, nrow(values$uniprot_results) > 0)
    sel <- input$uniprot_results_table_rows_selected
    req(length(sel) > 0)

    row <- values$uniprot_results[sel, ]
    fname <- generate_fasta_filename(row$upid, row$organism, input$fasta_content_type)

    div(style = "font-size: 0.8em; color: #6c757d; margin-top: 5px;",
      icon("file"), " ", fname
    )
  })

  # Download FASTA button handler
  observeEvent(input$download_fasta_btn, {
    req(values$uniprot_results, nrow(values$uniprot_results) > 0)
    sel <- input$uniprot_results_table_rows_selected

    if (length(sel) == 0) {
      showNotification("Please select a proteome from the table first.", type = "warning")
      return()
    }

    row <- values$uniprot_results[sel, ]
    fname <- generate_fasta_filename(row$upid, row$organism, input$fasta_content_type)

    # Determine output directory for FASTA
    fasta_dir <- file.path(Sys.getenv("HOME"), "proteomics_databases")
    output_path <- file.path(fasta_dir, fname)

    withProgress(message = sprintf("Downloading %s from UniProt...", row$upid), {
      result <- download_uniprot_fasta(
        proteome_id = row$upid,
        content_type = input$fasta_content_type,
        output_path = output_path
      )
    })

    if (result$success) {
      cfg <- ssh_config()
      if (!is.null(cfg)) {
        # SSH mode: check if FASTA already exists on remote
        remote_fasta_dir <- file.path(output_base(), "databases")
        remote_path <- file.path(remote_fasta_dir, fname)

        exists_check <- ssh_exec(cfg,
          paste("test -f", shQuote(remote_path), "&& echo EXISTS"))
        if (any(grepl("EXISTS", exists_check$stdout))) {
          # Already exists — just use it
          values$diann_fasta_files <- remote_path
          values$fasta_info <- result
          showNotification(
            sprintf("FASTA already exists on HPC, using: %s", remote_path),
            type = "message", duration = 8)
        } else {
          # Upload to remote
          ssh_exec(cfg, paste("mkdir -p", shQuote(remote_fasta_dir)))

          withProgress(message = "Uploading FASTA to remote HPC...", {
            up_result <- scp_upload(cfg, output_path, remote_path)
          })

          if (up_result$status != 0) {
            showNotification(
              paste("FASTA downloaded locally but upload to HPC failed:",
                    paste(up_result$stdout, collapse = " ")),
              type = "error", duration = 10)
            return()
          }

          values$diann_fasta_files <- remote_path
          values$fasta_info <- result
          showNotification(
            sprintf("FASTA uploaded to HPC: %d proteins (%.1f MB)\n%s",
              result$n_sequences,
              result$file_size / 1e6,
              remote_path),
            type = "message", duration = 10)
        }
      } else {
        # Local mode: use local path directly
        values$diann_fasta_files <- output_path
        values$fasta_info <- result
        showNotification(
          sprintf("FASTA downloaded: %d proteins (%.1f MB)",
            result$n_sequences,
            result$file_size / 1e6),
          type = "message", duration = 8
        )
      }
    } else {
      showNotification(paste("Download failed:", result$error), type = "error")
    }
  })

  # ============================================================================
  #    Pre-staged FASTA Selection
  # ============================================================================

  # Scan for pre-staged databases on startup
  observe({
    fasta_dir <- getOption("delimp.fasta_dir",
      default = "/share/proteomics/databases/fasta")
    databases <- scan_prestaged_databases(fasta_dir)
    if (length(databases) > 0) {
      updateSelectInput(session, "prestaged_fasta", choices = databases)
    }
  }) |> bindEvent(TRUE)  # Run once on startup

  output$prestaged_fasta_info <- renderUI({
    req(nzchar(input$prestaged_fasta))
    if (!file.exists(input$prestaged_fasta)) return(NULL)

    size_mb <- round(file.size(input$prestaged_fasta) / 1e6, 1)
    n_seqs <- tryCatch({
      sum(grepl("^>", readLines(input$prestaged_fasta, n = 200000, warn = FALSE)))
    }, error = function(e) NA)

    div(class = "alert alert-info",
      style = "margin-top: 8px; padding: 6px 10px; font-size: 0.82em;",
      icon("info-circle"),
      sprintf(" %s MB", size_mb),
      if (!is.na(n_seqs)) sprintf(", ~%d sequences", n_seqs) else ""
    )
  })

  observeEvent(input$prestaged_fasta, {
    req(nzchar(input$prestaged_fasta))
    values$diann_fasta_files <- input$prestaged_fasta
  })

  # ============================================================================
  #    Browsed FASTA Selection
  # ============================================================================

  observeEvent(input$fasta_browse_dir, {
    if (is.integer(input$fasta_browse_dir)) return()

    dir_path <- shinyFiles::parseDirPath(volumes, input$fasta_browse_dir)
    if (length(dir_path) == 0 || !nzchar(dir_path)) return()

    fasta_files <- list.files(as.character(dir_path),
      pattern = "\\.(fasta|fa)$", ignore.case = TRUE, full.names = TRUE)

    if (length(fasta_files) == 0) {
      showNotification("No FASTA files found in selected directory.", type = "warning")
      return()
    }

    # Use all FASTA files found (DIA-NN can take multiple)
    values$diann_fasta_files <- fasta_files
  })

  output$browsed_fasta_info <- renderUI({
    req(length(values$diann_fasta_files) > 0)

    n_files <- length(values$diann_fasta_files)
    fnames <- paste(basename(values$diann_fasta_files), collapse = ", ")

    div(class = "alert alert-success",
      style = "margin-top: 8px; padding: 8px; font-size: 0.85em;",
      icon("check-circle"),
      sprintf(" %d FASTA file%s: %s", n_files, if (n_files > 1) "s" else "", fnames)
    )
  })

  # ============================================================================
  #    Normalization Guidance
  # ============================================================================

  output$norm_guidance_search <- renderUI({
    if (input$diann_normalization == "on") {
      div(class = "alert alert-info",
        style = "padding: 6px 10px; font-size: 0.8em; margin-top: 5px;",
        icon("info-circle"),
        " RT-dependent normalization is recommended for standard proteomics experiments."
      )
    } else {
      div(class = "alert alert-warning",
        style = "padding: 6px 10px; font-size: 0.8em; margin-top: 5px;",
        icon("exclamation-triangle"),
        " Normalization OFF is recommended for AP-MS, Co-IP, or proximity labeling ",
        "where protein abundance differences are expected."
      )
    }
  })

  # ============================================================================
  #    Output Path Display
  # ============================================================================

  # Track selected output base directory
  output_base <- reactiveVal(file.path(Sys.getenv("HOME"), "diann_output"))

  observeEvent(input$output_base_dir, {
    if (is.integer(input$output_base_dir)) return()
    dir_path <- shinyFiles::parseDirPath(volumes, input$output_base_dir)
    if (length(dir_path) > 0 && nzchar(dir_path)) {
      output_base(as.character(dir_path))
    }
  })

  # SSH mode: update output base from text input
  observeEvent(input$ssh_output_base_dir, {
    if (nzchar(input$ssh_output_base_dir %||% "")) {
      output_base(input$ssh_output_base_dir)
    }
  })

  output$full_output_path <- renderText({
    base <- output_base()
    name <- gsub("[^A-Za-z0-9._-]", "_", input$analysis_name)
    file.path(base, name)
  })

  # ============================================================================
  #    Time Estimate
  # ============================================================================

  output$time_estimate_ui <- renderUI({
    req(values$diann_raw_files, nrow(values$diann_raw_files) > 0)

    est <- estimate_search_time(
      n_files = nrow(values$diann_raw_files),
      search_mode = input$search_mode,
      cpus = input$diann_cpus
    )

    div(class = "alert alert-info",
      style = "padding: 8px; font-size: 0.85em;",
      icon("clock"), " Estimated time: ", strong(est)
    )
  })

  # ============================================================================
  #    Job Submission
  # ============================================================================

  observeEvent(input$submit_diann, {
    tryCatch({

    # --- Validation ---
    errors <- character()

    if (is.null(values$diann_raw_files) || nrow(values$diann_raw_files) == 0) {
      errors <- c(errors, "No raw data files selected.")
    }
    has_fasta <- length(values$diann_fasta_files) > 0 &&
      all(nzchar(values$diann_fasta_files))
    has_speclib <- !is.null(values$diann_speclib) && nzchar(values$diann_speclib)
    if (!has_fasta && !has_speclib) {
      errors <- c(errors, "No FASTA database or spectral library selected.")
    }

    sif_path <- input$diann_sif_path
    cfg <- ssh_config()
    if (is.null(cfg)) {
      if (!file.exists(sif_path)) {
        errors <- c(errors, sprintf("DIA-NN container not found: %s", sif_path))
      }
    } else {
      # SSH mode: check remote file existence
      sif_check <- ssh_exec(cfg, paste("test -f", shQuote(sif_path), "&& echo EXISTS"))
      if (!any(grepl("EXISTS", sif_check$stdout))) {
        errors <- c(errors, sprintf("DIA-NN container not found on remote: %s", sif_path))
      }
    }

    if (!nzchar(input$analysis_name)) {
      errors <- c(errors, "Analysis name is required.")
    }

    if (length(errors) > 0) {
      showNotification(
        HTML(paste("<b>Cannot submit:</b><br>",
          paste("&bull;", errors, collapse = "<br>"))),
        type = "error", duration = 10
      )
      return()
    }

    # --- Prepare submission ---
    analysis_name <- gsub("[^A-Za-z0-9._-]", "_", input$analysis_name)
    output_dir <- file.path(output_base(), analysis_name)

    if (!is.null(cfg)) {
      ssh_exec(cfg, paste("mkdir -p", shQuote(output_dir)))
    } else {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    # Collect search params
    search_params <- list(
      qvalue = input$diann_fdr %||% 0.01,
      max_var_mods = input$diann_max_var_mods,
      scan_window = input$diann_scan_window %||% 6,
      mass_acc_mode = input$mass_acc_mode,
      mass_acc = input$diann_mass_acc %||% 14,
      mass_acc_ms1 = input$diann_mass_acc_ms1 %||% 14,
      unimod4 = input$diann_unimod4 %||% TRUE,
      met_excision = input$diann_met_excision %||% TRUE,
      min_pep_len = input$min_pep_len %||% 7,
      max_pep_len = input$max_pep_len %||% 30,
      min_pr_mz = input$min_pr_mz %||% 300,
      max_pr_mz = input$max_pr_mz %||% 1200,
      min_pr_charge = 1, max_pr_charge = 4,
      min_fr_mz = 200, max_fr_mz = 1200,
      enzyme = input$diann_enzyme,
      missed_cleavages = input$diann_missed_cleavages,
      mbr = input$diann_mbr %||% TRUE,
      rt_profiling = input$diann_rt_profiling %||% TRUE,
      xic = input$diann_xic %||% TRUE,
      mod_met_ox = input$mod_met_ox,
      mod_nterm_acetyl = input$mod_nterm_acetyl,
      extra_var_mods = input$extra_var_mods %||% "",
      extra_cli_flags = input$extra_cli_flags %||% ""
    )

    # Handle contaminant library — add as separate FASTA file
    fasta_files <- values$diann_fasta_files
    contam_lib <- input$contaminant_library
    if (!is.null(contam_lib) && contam_lib != "none") {
      contam_result <- get_contaminant_fasta(contam_lib)

      if (contam_result$success) {
        if (!is.null(cfg)) {
          # SSH mode: upload contaminant FASTA to same remote dir as proteome
          remote_contam_dir <- file.path(output_base(), "databases")
          remote_contam_path <- file.path(remote_contam_dir, basename(contam_result$path))

          exists_check <- ssh_exec(cfg,
            paste("test -f", shQuote(remote_contam_path), "&& echo EXISTS"))
          if (!any(grepl("EXISTS", exists_check$stdout))) {
            ssh_exec(cfg, paste("mkdir -p", shQuote(remote_contam_dir)))
            scp_upload(cfg, contam_result$path, remote_contam_path)
          }
          fasta_files <- c(fasta_files, remote_contam_path)
        } else {
          fasta_files <- c(fasta_files, contam_result$path)
        }
        showNotification(
          sprintf("Added %s contaminant library (%d proteins)",
                  gsub("_", " ", contam_lib), contam_result$n_sequences),
          type = "message", duration = 5)
      } else {
        showNotification(
          paste("Warning: Contaminant library not found:", contam_result$error),
          type = "warning", duration = 8)
      }
    }

    # Generate sbatch script
    script_content <- generate_sbatch_script(
      analysis_name = analysis_name,
      raw_files = values$diann_raw_files$full_path,
      fasta_files = fasta_files,
      speclib_path = values$diann_speclib,
      output_dir = output_dir,
      diann_sif = sif_path,
      normalization = input$diann_normalization,
      search_mode = input$search_mode,
      cpus = input$diann_cpus,
      mem_gb = input$diann_mem_gb,
      time_hours = input$diann_time_hours,
      partition = input$diann_partition,
      account = input$diann_account,
      search_params = search_params
    )

    # Write sbatch script and submit
    script_path <- file.path(output_dir, "diann_search.sbatch")

    if (!is.null(cfg)) {
      # SSH mode: write script locally, SCP to remote, then submit
      local_tmp <- tempfile(fileext = ".sbatch")
      writeLines(script_content, local_tmp)
      on.exit(unlink(local_tmp), add = TRUE)

      scp_result <- scp_upload(cfg, local_tmp, script_path)
      if (scp_result$status != 0) {
        showNotification(
          paste("Failed to write sbatch script to remote host:",
                paste(scp_result$stdout, collapse = " ")),
          type = "error")
        return()
      }

      # Use stored full sbatch path to avoid slow login shell initialization
      sbatch_bin <- values$ssh_sbatch_path %||% "sbatch"
      sbatch_cmd <- paste(sbatch_bin, shQuote(script_path))
      submit_result <- tryCatch({
        result <- ssh_exec(cfg, sbatch_cmd,
                           login_shell = is.null(values$ssh_sbatch_path))
        list(success = result$status == 0, stdout = result$stdout,
             error = if (result$status != 0) paste(result$stdout, collapse = " ") else NULL)
      }, error = function(e) {
        list(success = FALSE, error = e$message)
      })
    } else {
      # Local mode: write and submit locally
      writeLines(script_content, script_path)

      submit_result <- tryCatch({
        stdout <- system2("sbatch", args = script_path, stdout = TRUE, stderr = TRUE)
        list(success = TRUE, stdout = stdout)
      }, error = function(e) {
        list(success = FALSE, error = e$message)
      })
    }

    if (!submit_result$success) {
      showNotification(paste("sbatch submission failed:", submit_result$error), type = "error")
      return()
    }

    job_id <- parse_sbatch_output(submit_result$stdout)
    if (is.null(job_id)) {
      showNotification(
        paste("Could not parse job ID from sbatch output:",
          paste(submit_result$stdout, collapse = " ")),
        type = "error"
      )
      return()
    }

    # Add to job queue
    job_entry <- list(
      job_id = job_id,
      name = analysis_name,
      status = "queued",
      output_dir = output_dir,
      script_path = script_path,
      submitted_at = Sys.time(),
      n_files = nrow(values$diann_raw_files),
      search_mode = input$search_mode,
      search_settings = list(
        search_params = search_params,
        fasta_files = fasta_files,
        contaminant_library = contam_lib,
        n_raw_files = nrow(values$diann_raw_files),
        raw_file_type = if (nrow(values$diann_raw_files) > 0)
          tools::file_ext(values$diann_raw_files$name[1]) else "unknown",
        search_mode = input$search_mode,
        normalization = input$diann_normalization,
        diann_sif = basename(sif_path),
        speclib = if (!is.null(values$diann_speclib) && nzchar(values$diann_speclib))
          basename(values$diann_speclib) else NULL,
        slurm = list(
          cpus = input$diann_cpus,
          mem_gb = input$diann_mem_gb,
          time_hours = input$diann_time_hours,
          partition = input$diann_partition
        )
      ),
      auto_load = input$auto_load_results,
      log_content = "",
      completed_at = NULL,
      loaded = FALSE,
      is_ssh = !is.null(cfg)
    )

    values$diann_jobs <- c(values$diann_jobs, list(job_entry))

    # Log the submission
    add_to_log("DIA-NN Search Submitted", c(
      sprintf("# Job ID: %s", job_id),
      sprintf("# Analysis: %s", analysis_name),
      sprintf("# Files: %d raw data files", nrow(values$diann_raw_files)),
      sprintf("# Mode: %s", input$search_mode),
      sprintf("# Output: %s", output_dir),
      sprintf("# sbatch script: %s", script_path)
    ))

    showNotification(
      sprintf("Job %s submitted successfully! Monitoring in background.", job_id),
      type = "message", duration = 8
    )

    }, error = function(e) {
      showNotification(
        paste("Submission error:", e$message),
        type = "error", duration = 15
      )
      message("[DE-LIMP] Submit error: ", e$message)
    })
  })

  # ============================================================================
  #    Job Monitoring (polls every 15 seconds)
  # ============================================================================

  observe({
    req(length(values$diann_jobs) > 0)

    # Only poll if there are active jobs
    active_jobs <- vapply(values$diann_jobs, function(j) {
      j$status %in% c("queued", "running")
    }, logical(1))

    if (!any(active_jobs)) return()

    invalidateLater(15000)  # Poll every 15 seconds

    jobs <- values$diann_jobs
    changed <- FALSE

    # Get SSH config once for this polling cycle
    cfg <- isolate(ssh_config())

    for (i in seq_along(jobs)) {
      if (!jobs[[i]]$status %in% c("queued", "running")) next

      # Use SSH for jobs submitted via SSH, local for local jobs
      job_cfg <- if (isTRUE(jobs[[i]]$is_ssh)) cfg else NULL
      new_status <- check_slurm_status(jobs[[i]]$job_id, ssh_config = job_cfg,
                                        sbatch_path = values$ssh_sbatch_path)

      if (new_status != jobs[[i]]$status) {
        jobs[[i]]$status <- new_status
        changed <- TRUE

        if (new_status == "completed") {
          jobs[[i]]$completed_at <- Sys.time()
          showNotification(
            sprintf("DIA-NN search '%s' completed!", jobs[[i]]$name),
            type = "message", duration = 15
          )
        } else if (new_status == "failed") {
          showNotification(
            sprintf("DIA-NN search '%s' failed. Check log for details.", jobs[[i]]$name),
            type = "error", duration = 15
          )
        }
      }

      # Tail the log file (local or remote)
      if (isTRUE(jobs[[i]]$is_ssh) && !is.null(cfg)) {
        log_result <- ssh_exec(cfg, sprintf(
          "ls -t %s/diann_*.{out,err} 2>/dev/null | head -1 | xargs tail -50 2>/dev/null",
          shQuote(jobs[[i]]$output_dir)))
        if (log_result$status == 0 && length(log_result$stdout) > 0) {
          jobs[[i]]$log_content <- paste(log_result$stdout, collapse = "\n")
          changed <- TRUE
        }
      } else {
        log_files <- list.files(jobs[[i]]$output_dir,
          pattern = "^diann_.*\\.(out|err)$", full.names = TRUE)
        if (length(log_files) > 0) {
          log_file <- log_files[which.max(file.mtime(log_files))]
          log_lines <- tryCatch(
            tail(readLines(log_file, warn = FALSE), 50),
            error = function(e) character(0)
          )
          jobs[[i]]$log_content <- paste(log_lines, collapse = "\n")
          changed <- TRUE
        }
      }
    }

    if (changed) {
      values$diann_jobs <- jobs
    }
  })

  # ============================================================================
  #    Auto-Load Results
  # ============================================================================

  observe({
    req(length(values$diann_jobs) > 0)

    cfg <- isolate(ssh_config())

    for (i in seq_along(values$diann_jobs)) {
      job <- values$diann_jobs[[i]]

      if (job$status != "completed" || !isTRUE(job$auto_load) || isTRUE(job$loaded)) next

      # Look for report.parquet in output directory
      report_name <- if (grepl("no_norm", job$name, ignore.case = TRUE)) {
        "no_norm_report.parquet"
      } else {
        "report.parquet"
      }

      remote_report <- file.path(job$output_dir, report_name)

      if (isTRUE(job$is_ssh) && !is.null(cfg)) {
        # SSH mode: check remote, then SCP download
        find_result <- ssh_exec(cfg, paste("ls", shQuote(remote_report), "2>/dev/null"))
        if (find_result$status != 0) {
          # Try finding any report parquet
          find_result <- ssh_exec(cfg, sprintf(
            "ls %s/report*.parquet 2>/dev/null | head -1", shQuote(job$output_dir)))
          if (find_result$status != 0 || length(find_result$stdout) == 0 ||
              !nzchar(trimws(find_result$stdout[1]))) next
          remote_report <- trimws(find_result$stdout[1])
        }

        local_report <- file.path(tempdir(), paste0(job$name, "_", basename(remote_report)))
        dl_result <- scp_download(cfg, remote_report, local_report)
        if (dl_result$status != 0) {
          showNotification(sprintf("SCP download failed for '%s'.", job$name),
            type = "error", duration = 10)
          next
        }
        report_path <- local_report
      } else {
        # Local mode: direct file access
        report_path <- remote_report
        if (!file.exists(report_path)) {
          parquet_files <- list.files(job$output_dir, pattern = "report.*\\.parquet$",
            full.names = TRUE)
          if (length(parquet_files) > 0) {
            report_path <- parquet_files[1]
          } else {
            next
          }
        }
      }

      # Load the results into DE-LIMP pipeline
      tryCatch({
        withProgress(message = sprintf("Loading results from %s...", job$name), {
          raw_data <- limpa::readDIANN(report_path)

          values$raw_data <- raw_data
          values$uploaded_report_path <- report_path
          values$original_report_name <- basename(report_path)

          # Initialize metadata from raw_data
          sample_names <- colnames(raw_data$E)
          values$metadata <- data.frame(
            ID = seq_along(sample_names),
            File.Name = sample_names,
            Group = "",
            Batch = "",
            Covariate1 = "",
            Covariate2 = "",
            stringsAsFactors = FALSE
          )

          # Run phospho detection
          tryCatch({
            report_df <- arrow::read_parquet(report_path,
              col_select = c("Modified.Sequence"))
            values$phospho_detected <- detect_phospho(report_df)
          }, error = function(e) NULL)

          # Check for XIC files (local mode only)
          if (!isTRUE(job$is_ssh)) {
            xic_dir <- paste0(tools::file_path_sans_ext(report_path), "_xic")
            if (dir.exists(xic_dir)) {
              values$xic_dir <- xic_dir
              values$xic_available <- TRUE
            }
          }

          # Save search settings for methodology tab
          if (!is.null(job$search_settings)) {
            values$diann_search_settings <- job$search_settings
          }

          # Mark job as loaded
          jobs <- values$diann_jobs
          jobs[[i]]$loaded <- TRUE
          values$diann_jobs <- jobs

          # Build log with key search parameters
          log_lines <- c(
            sprintf("# Loaded from: %s", report_path),
            sprintf("# Job ID: %s, Analysis: %s", job$job_id, job$name),
            sprintf("# Mode: %s", if (isTRUE(job$is_ssh)) "SSH (SCP download)" else "Local")
          )
          if (!is.null(job$search_settings)) {
            ss <- job$search_settings
            sp <- ss$search_params
            log_lines <- c(log_lines,
              sprintf("# Search mode: %s", ss$search_mode),
              sprintf("# FASTA: %s", paste(basename(ss$fasta_files), collapse = ", ")),
              sprintf("# Enzyme: %s, Missed cleavages: %d", sp$enzyme, sp$missed_cleavages),
              sprintf("# FDR: %s, MBR: %s", sp$qvalue, sp$mbr)
            )
          }
          log_lines <- c(log_lines,
            sprintf("raw_data <- limpa::readDIANN('%s')", report_path)
          )
          add_to_log("Auto-Load DIA-NN Results", log_lines)

          # Navigate to Assign Groups tab
          nav_select("main_tabs", "Data Overview")
          nav_select("data_overview_tabs", "Assign Groups & Run")

          showNotification(
            sprintf("Results loaded from '%s'! Assign groups and run the pipeline.", job$name),
            type = "message", duration = 10
          )
        })
      }, error = function(e) {
        showNotification(
          sprintf("Failed to auto-load results from '%s': %s", job$name, e$message),
          type = "error", duration = 10
        )
      })
    }
  })

  # ============================================================================
  #    Job Queue UI
  # ============================================================================

  output$search_queue_ui <- renderUI({
    jobs <- values$diann_jobs
    if (length(jobs) == 0) {
      return(div(style = "color: #999; font-size: 0.85em; text-align: center; padding: 10px;",
        "No jobs submitted yet."
      ))
    }

    # Refresh all button at top
    has_unknown <- any(vapply(jobs, function(j) j$status == "unknown", logical(1)))

    job_rows <- lapply(seq_along(jobs), function(i) {
      job <- jobs[[i]]

      status_badge <- switch(job$status,
        "queued"    = span(class = "badge bg-secondary", "Queued"),
        "running"   = span(class = "badge bg-primary", "Running"),
        "completed" = span(class = "badge bg-success", "Completed"),
        "failed"    = span(class = "badge bg-danger", "Failed"),
        "cancelled" = span(class = "badge bg-warning", "Cancelled"),
        "unknown"   = span(class = "badge bg-light text-dark", "Unknown"),
        span(class = "badge bg-light text-dark", job$status)
      )

      elapsed <- if (!is.null(job$completed_at)) {
        difftime(job$completed_at, job$submitted_at, units = "mins")
      } else {
        difftime(Sys.time(), job$submitted_at, units = "mins")
      }
      elapsed_str <- if (as.numeric(elapsed) < 60) {
        sprintf("%.0f min", as.numeric(elapsed))
      } else {
        sprintf("%.1f hrs", as.numeric(elapsed) / 60)
      }

      div(style = "border: 1px solid #dee2e6; border-radius: 5px; padding: 8px; margin-bottom: 8px; font-size: 0.82em;",
        div(style = "display: flex; justify-content: space-between; align-items: center;",
          div(
            strong(job$name), " ",
            span(style = "color: #999;", sprintf("(#%s)", job$job_id))
          ),
          status_badge
        ),
        div(style = "display: flex; justify-content: space-between; align-items: center; margin-top: 4px;",
          span(style = "color: #666;",
            sprintf("%d files | %s", job$n_files, elapsed_str)
          ),
          div(style = "display: flex; gap: 4px;",
            actionButton(sprintf("view_log_%d", i), "Log",
              class = "btn-outline-secondary btn-xs",
              style = "font-size: 0.75em; padding: 2px 6px;"),
            if (job$status == "unknown") {
              actionButton(sprintf("refresh_job_%d", i), "Refresh",
                class = "btn-outline-info btn-xs",
                style = "font-size: 0.75em; padding: 2px 6px;")
            },
            if (job$status %in% c("queued", "running")) {
              actionButton(sprintf("cancel_job_%d", i), "Cancel",
                class = "btn-outline-danger btn-xs",
                style = "font-size: 0.75em; padding: 2px 6px;")
            },
            if (job$status == "completed" && !isTRUE(job$loaded)) {
              actionButton(sprintf("load_results_%d", i), "Load",
                class = "btn-outline-success btn-xs",
                style = "font-size: 0.75em; padding: 2px 6px;")
            }
          )
        )
      )
    })

    tagList(
      if (has_unknown) div(style = "text-align: right; margin-bottom: 6px;",
        actionButton("refresh_all_jobs", "Refresh All",
          class = "btn-outline-info btn-xs",
          style = "font-size: 0.75em; padding: 2px 8px;",
          icon = icon("arrows-rotate"))
      ),
      job_rows
    )
  })

  # ============================================================================
  #    Dynamic Observers for Job Queue Buttons
  # ============================================================================

  # Track which observers have been registered to avoid duplicates
  registered_observers <- reactiveVal(character())

  observe({
    jobs <- values$diann_jobs
    existing <- registered_observers()

    for (i in seq_along(jobs)) {
      job_key <- as.character(i)
      if (job_key %in% existing) next

      local({
        idx <- i

        # View log modal
        observeEvent(input[[sprintf("view_log_%d", idx)]], {
          job <- values$diann_jobs[[idx]]
          showModal(modalDialog(
            title = sprintf("Log: %s (#%s)", job$name, job$job_id),
            size = "l", easyClose = TRUE, footer = modalButton("Close"),
            pre(style = "max-height: 500px; overflow-y: auto; font-size: 0.8em;",
              job$log_content
            )
          ))
        }, ignoreInit = TRUE)

        # Refresh job status
        observeEvent(input[[sprintf("refresh_job_%d", idx)]], {
          job <- values$diann_jobs[[idx]]
          tryCatch({
            job_cfg <- if (isTRUE(job$is_ssh)) isolate(ssh_config()) else NULL
            new_status <- check_slurm_status(job$job_id, ssh_config = job_cfg,
                                              sbatch_path = values$ssh_sbatch_path)
            jobs <- values$diann_jobs
            jobs[[idx]]$status <- new_status
            if (new_status == "completed" && is.null(jobs[[idx]]$completed_at)) {
              jobs[[idx]]$completed_at <- Sys.time()
            }
            values$diann_jobs <- jobs
            showNotification(sprintf("Job %s: %s", job$job_id, new_status), type = "message")
          }, error = function(e) {
            showNotification(sprintf("Refresh failed: %s", e$message), type = "error")
          })
        }, ignoreInit = TRUE)

        # Cancel job
        observeEvent(input[[sprintf("cancel_job_%d", idx)]], {
          job <- values$diann_jobs[[idx]]
          tryCatch({
            if (isTRUE(job$is_ssh)) {
              cfg <- ssh_config()
              if (!is.null(cfg)) {
                scancel_cmd <- if (!is.null(values$ssh_sbatch_path)) {
                  file.path(dirname(values$ssh_sbatch_path), "scancel")
                } else "scancel"
                ssh_exec(cfg, paste(scancel_cmd, job$job_id))
              }
            } else {
              system2("scancel", args = job$job_id, stdout = TRUE, stderr = TRUE)
            }
            jobs <- values$diann_jobs
            jobs[[idx]]$status <- "cancelled"
            values$diann_jobs <- jobs
            showNotification(sprintf("Job %s cancelled.", job$job_id), type = "message")
          }, error = function(e) {
            showNotification(sprintf("Failed to cancel job: %s", e$message), type = "error")
          })
        }, ignoreInit = TRUE)

        # Manual load results
        observeEvent(input[[sprintf("load_results_%d", idx)]], {
          job <- values$diann_jobs[[idx]]

          tryCatch({
            if (isTRUE(job$is_ssh)) {
              # SSH mode: SCP download first
              cfg <- ssh_config()
              if (is.null(cfg)) {
                showNotification("SSH not configured. Re-enter connection details.", type = "error")
                return()
              }

              remote_report <- file.path(job$output_dir, "report.parquet")
              find_result <- ssh_exec(cfg, paste("ls", shQuote(remote_report), "2>/dev/null"))
              if (find_result$status != 0) {
                find_result <- ssh_exec(cfg, sprintf(
                  "ls %s/report*.parquet 2>/dev/null | head -1", shQuote(job$output_dir)))
                if (find_result$status != 0 || length(find_result$stdout) == 0 ||
                    !nzchar(trimws(find_result$stdout[1]))) {
                  showNotification("No report.parquet found on remote.", type = "error")
                  return()
                }
                remote_report <- trimws(find_result$stdout[1])
              }

              withProgress(message = "Downloading results via SCP...", {
                local_report <- file.path(tempdir(), paste0(job$name, "_", basename(remote_report)))
                dl_result <- scp_download(cfg, remote_report, local_report)
                if (dl_result$status != 0) {
                  showNotification("SCP download failed.", type = "error")
                  return()
                }
                report_path <- local_report
              })
            } else {
              # Local mode: direct access
              report_path <- file.path(job$output_dir, "report.parquet")
              if (!file.exists(report_path)) {
                parquet_files <- list.files(job$output_dir, pattern = "report.*\\.parquet$",
                  full.names = TRUE)
                if (length(parquet_files) > 0) {
                  report_path <- parquet_files[1]
                } else {
                  showNotification("No report.parquet found in output directory.", type = "error")
                  return()
                }
              }
            }

            raw_data <- limpa::readDIANN(report_path)
            values$raw_data <- raw_data
            values$uploaded_report_path <- report_path
            values$original_report_name <- basename(report_path)

            sample_names <- colnames(raw_data$E)
            values$metadata <- data.frame(
              ID = seq_along(sample_names),
              File.Name = sample_names,
              Group = "", Batch = "",
              Covariate1 = "", Covariate2 = "",
              stringsAsFactors = FALSE
            )

            jobs <- values$diann_jobs
            jobs[[idx]]$loaded <- TRUE
            values$diann_jobs <- jobs

            nav_select("main_tabs", "Data Overview")
            nav_select("data_overview_tabs", "Assign Groups & Run")

            showNotification("Results loaded! Assign groups and run pipeline.",
              type = "message")
          }, error = function(e) {
            showNotification(paste("Failed to load:", e$message), type = "error")
          })
        }, ignoreInit = TRUE)
      })

      existing <- c(existing, job_key)
    }

    registered_observers(existing)
  })

  # Refresh all jobs with unknown status
  observeEvent(input$refresh_all_jobs, {
    jobs <- values$diann_jobs
    cfg <- isolate(ssh_config())
    changed <- FALSE

    for (i in seq_along(jobs)) {
      if (jobs[[i]]$status != "unknown") next
      tryCatch({
        job_cfg <- if (isTRUE(jobs[[i]]$is_ssh)) cfg else NULL
        new_status <- check_slurm_status(jobs[[i]]$job_id, ssh_config = job_cfg,
                                          sbatch_path = values$ssh_sbatch_path)
        jobs[[i]]$status <- new_status
        if (new_status == "completed" && is.null(jobs[[i]]$completed_at)) {
          jobs[[i]]$completed_at <- Sys.time()
        }
        changed <- TRUE
      }, error = function(e) NULL)
    }

    if (changed) values$diann_jobs <- jobs
    showNotification("Job statuses refreshed.", type = "message", duration = 3)
  })

}
