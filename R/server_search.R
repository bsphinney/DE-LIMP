# ==============================================================================
#  server_search.R
#  DIA-NN Search Integration — New Search tab server logic
#  Supports three backends: Local embedded, Local Docker, and HPC (SSH/SLURM).
#  Handles: file browsing, UniProt FASTA download, sbatch generation,
#  Docker execution, local execution, job submission, monitoring, auto-load, and job queue.
# ==============================================================================

server_search <- function(input, output, session, values, add_to_log,
                          search_enabled, docker_available, docker_config,
                          hpc_available, local_sbatch,
                          local_diann = FALSE, delimp_data_dir = "",
                          is_core_facility = FALSE, cf_config = NULL) {

  # Early return if no search backend available
  if (!search_enabled) return(invisible())

  # ============================================================================
  #    SSH Config Reactive (HPC backend only)
  # ============================================================================

  ssh_config <- reactive({
    if (is.null(input$search_backend) || input$search_backend != "hpc") return(NULL)
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
  #    Docker Backend UI (image status, resource controls, output path)
  # ============================================================================

  # Docker image status
  output$docker_image_status <- renderUI({
    if (!docker_available) return(NULL)
    img <- input$docker_image_name %||% docker_config$diann_image %||% "diann:2.0"
    result <- check_diann_image(img)

    if (result$exists) {
      # Image found — check for ARM/Rosetta
      arch <- Sys.info()[["machine"]]
      arm_warning <- if (arch %in% c("arm64", "aarch64")) {
        tags$div(class = "alert alert-warning py-1 px-2 mt-1",
          style = "font-size: 0.82em;",
          icon("triangle-exclamation"),
          tags$strong(" Apple Silicon detected."),
          " DIA-NN runs under Rosetta 2 emulation (~3-5x slower). ",
          "Fine for small datasets; use HPC for large experiments.")
      }
      tagList(
        tags$div(class = "alert alert-success py-1 px-2",
          style = "font-size: 0.85em;",
          icon("check-circle"),
          sprintf(" DIA-NN Docker image ready: %s", img)),
        arm_warning
      )
    } else {
      tags$div(class = "alert alert-warning py-2 px-3",
        icon("docker"),
        tags$strong(" DIA-NN Docker image not found."),
        tags$p("Image ", tags$code(img), " is not available locally. ",
          "DIA-NN must be built locally due to licensing restrictions."),
        tags$p("Run the build script included with DE-LIMP:"),
        tags$pre(style = "font-size: 0.8em; margin-bottom: 4px;",
          "bash build_diann_docker.sh"),
        tags$small(class = "text-muted",
          "See ", tags$a(href = "https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md",
                         "DIA-NN license", target = "_blank"), " for terms.")
      )
    }
  })

  # Docker resource controls (CPU/memory sliders)
  output$docker_resources_ui <- renderUI({
    res <- get_host_resources()
    max_cpus <- res$cpus
    max_mem <- res$memory_gb
    tagList(
      div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
        div(style = "flex: 1; min-width: 150px;",
          sliderInput("docker_cpus", "CPUs:",
            min = 1, max = max_cpus,
            value = min(max_cpus, 16), step = 1)
        ),
        div(style = "flex: 1; min-width: 150px;",
          sliderInput("docker_mem_gb", "Memory (GB):",
            min = 4, max = max_mem,
            value = min(max_mem, 64), step = 4)
        )
      ),
      tags$p(class = "text-muted", style = "font-size: 0.8em;",
        sprintf("System: %d CPUs, %d GB RAM. Leave headroom for OS + DE-LIMP.", max_cpus, max_mem))
    )
  })

  # Docker output path display
  output$docker_output_path <- renderText({
    dir_chosen <- shinyFiles::parseDirPath(volumes, input$docker_output_dir)
    if (length(dir_chosen) > 0) as.character(dir_chosen) else "(not selected)"
  })

  # ============================================================================
  #    Local (Embedded) Backend UI
  # ============================================================================

  # Local resource controls (threads slider)
  output$local_resources_ui <- renderUI({
    n_cores <- parallel::detectCores(logical = TRUE)
    sliderInput("local_diann_threads", "Threads:",
      min = 1, max = n_cores,
      value = min(n_cores, 16), step = 1)
  })

  # Local output path display (native mode — container mode uses fixed textInput)
  output$local_output_path <- renderText({
    dir_chosen <- shinyFiles::parseDirPath(volumes, input$local_output_dir_browse)
    if (length(dir_chosen) > 0) as.character(dir_chosen) else "(not selected)"
  })

  # Local output dir observer (native mode — update output_base when user picks a folder)
  observeEvent(input$local_output_dir_browse, {
    if (is.integer(input$local_output_dir_browse)) return()
    dir_path <- shinyFiles::parseDirPath(volumes, input$local_output_dir_browse)
    if (length(dir_path) > 0 && nzchar(dir_path)) {
      output_base(as.character(dir_path))
    }
  })

  # Local output dir observer (container mode — update output_base from text input)
  observeEvent(input$local_output_dir, {
    if (nzchar(input$local_output_dir %||% "")) {
      output_base(input$local_output_dir)
    }
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
          # Backward compat: add backend field to old jobs
          for (i in seq_along(saved_jobs)) {
            if (is.null(saved_jobs[[i]]$backend)) saved_jobs[[i]]$backend <- "hpc"
          }
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
  # CRITICAL: ignoreInit = TRUE prevents overwriting saved jobs with the empty
  # initial value of values$diann_jobs before the load observer restores them.
  observeEvent(values$diann_jobs, {
    req(job_queue_loaded())
    tryCatch({
      saveRDS(values$diann_jobs, job_queue_path)
    }, error = function(e) {
      message("[DE-LIMP] Failed to save job queue: ", e$message)
    })
  }, ignoreInit = TRUE)

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

  volumes <- if (nzchar(delimp_data_dir)) {
    c(Data = delimp_data_dir)
  } else {
    c(Home = Sys.getenv("HOME"), Root = "/")
  }

  shinyFiles::shinyDirChoose(input, "raw_data_dir", roots = volumes, session = session)
  shinyFiles::shinyDirChoose(input, "fasta_browse_dir", roots = volumes, session = session)
  shinyFiles::shinyDirChoose(input, "output_base_dir", roots = volumes, session = session)
  shinyFiles::shinyFileChoose(input, "lib_file", roots = volumes, session = session,
    filetypes = c("speclib", "tsv", "csv"))
  shinyFiles::shinyDirChoose(input, "docker_output_dir", roots = volumes, session = session)
  if (local_diann && !nzchar(delimp_data_dir)) {
    shinyFiles::shinyDirChoose(input, "local_output_dir_browse", roots = volumes, session = session)
  }

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

  # Docker mode: update output base from directory chooser
  observeEvent(input$docker_output_dir, {
    if (is.integer(input$docker_output_dir)) return()
    dir_path <- shinyFiles::parseDirPath(volumes, input$docker_output_dir)
    if (length(dir_path) > 0 && nzchar(dir_path)) {
      output_base(as.character(dir_path))
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

    backend <- input$search_backend %||% "hpc"

    # --- Validation (shared) ---
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
    if (!nzchar(input$analysis_name)) {
      errors <- c(errors, "Analysis name is required.")
    }

    # Backend-specific validation
    if (backend == "local") {
      diann_bin <- Sys.which("diann")
      if (!nzchar(diann_bin)) diann_bin <- Sys.which("diann-linux")
      if (!nzchar(diann_bin)) {
        errors <- c(errors, "DIA-NN binary not found on PATH.")
      }
    } else if (backend == "docker") {
      img <- input$docker_image_name %||% docker_config$diann_image %||% "diann:2.0"
      img_check <- check_diann_image(img)
      if (!img_check$exists) {
        errors <- c(errors, sprintf(
          "DIA-NN Docker image '%s' not found. Run build_diann_docker.sh first.", img))
      }
    } else {
      sif_path <- input$diann_sif_path
      cfg <- ssh_config()
      if (is.null(cfg)) {
        if (!file.exists(sif_path)) {
          errors <- c(errors, sprintf("DIA-NN container not found: %s", sif_path))
        }
      } else {
        sif_check <- ssh_exec(cfg, paste("test -f", shQuote(sif_path), "&& echo EXISTS"))
        if (!any(grepl("EXISTS", sif_check$stdout))) {
          errors <- c(errors, sprintf("DIA-NN container not found on remote: %s", sif_path))
        }
      }
    }

    if (length(errors) > 0) {
      showNotification(
        HTML(paste("<b>Cannot submit:</b><br>",
          paste("&bull;", errors, collapse = "<br>"))),
        type = "error", duration = 10
      )
      return()
    }

    # --- Prepare submission (shared) ---
    analysis_name <- gsub("[^A-Za-z0-9._-]", "_", input$analysis_name)
    output_dir <- file.path(output_base(), analysis_name)

    if (backend == "local") {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      cfg <- NULL
    } else if (backend == "docker") {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      cfg <- NULL
    } else {
      cfg <- ssh_config()
      if (!is.null(cfg)) {
        ssh_exec(cfg, paste("mkdir -p", shQuote(output_dir)))
      } else {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      }
    }

    # Collect search params (shared between backends)
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
        if (backend == "hpc" && !is.null(cfg)) {
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
          # Docker or local HPC: use local path directly
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

    # ====================================================================
    #  Backend-specific submission
    # ====================================================================

    if (backend == "local") {
      # --- Local (embedded) submission via processx ---
      threads <- input$local_diann_threads %||% 4

      speclib_path <- if (!is.null(values$diann_speclib) && nzchar(values$diann_speclib)) {
        values$diann_speclib
      } else NULL
      diann_flags <- build_diann_flags(search_params, input$search_mode,
                                        input$diann_normalization, speclib_path)

      log_file <- file.path(output_dir, paste0("diann_", analysis_name, ".log"))

      submit_result <- tryCatch({
        result <- run_local_diann(
          raw_files = values$diann_raw_files$full_path,
          fasta_files = fasta_files,
          output_dir = output_dir,
          diann_flags = diann_flags,
          threads = threads,
          log_file = log_file,
          speclib_path = speclib_path
        )
        list(success = TRUE, process = result$process, pid = result$pid, log_file = result$log_file)
      }, error = function(e) {
        list(success = FALSE, error = e$message)
      })

      if (!submit_result$success) {
        showNotification(paste("Local DIA-NN launch failed:", submit_result$error),
          type = "error", duration = 15)
        return()
      }

      job_id <- sprintf("local_%s_%s", analysis_name, format(Sys.time(), "%Y%m%d_%H%M%S"))

      # Create local job entry
      job_entry <- list(
        job_id = job_id,
        backend = "local",
        name = analysis_name,
        status = "running",
        output_dir = output_dir,
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
          speclib = if (!is.null(values$diann_speclib) && nzchar(values$diann_speclib))
            basename(values$diann_speclib) else NULL,
          local = list(threads = threads)
        ),
        auto_load = input$auto_load_results,
        log_content = "",
        log_file = log_file,
        pid = submit_result$pid,
        process = submit_result$process,
        completed_at = NULL,
        loaded = FALSE,
        is_ssh = FALSE
      )

    } else if (backend == "docker") {
      # --- Docker submission ---
      img <- input$docker_image_name %||% docker_config$diann_image %||% "diann:2.0"
      cpus <- input$docker_cpus %||% 8
      mem_gb <- input$docker_mem_gb %||% 32

      # Build DIA-NN flags (shared with HPC via build_diann_flags)
      speclib_mount <- if (!is.null(values$diann_speclib) && nzchar(values$diann_speclib)) {
        sprintf("/work/lib/%s", basename(values$diann_speclib))
      } else NULL
      diann_flags <- build_diann_flags(search_params, input$search_mode,
                                        input$diann_normalization, speclib_mount)

      # Generate unique container name
      container_name <- sprintf("delimp_%s_%s", analysis_name,
                                 format(Sys.time(), "%Y%m%d_%H%M%S"))

      # Build docker run command
      docker_args <- build_docker_command(
        raw_files = values$diann_raw_files$full_path,
        fasta_files = fasta_files,
        output_dir = output_dir,
        image_name = img,
        diann_flags = diann_flags,
        cpus = cpus,
        mem_gb = mem_gb,
        container_name = container_name,
        speclib_path = values$diann_speclib
      )

      # Launch Docker container (detached mode — returns container ID)
      submit_result <- tryCatch({
        stdout <- suppressWarnings(
          system2("docker", args = docker_args, stdout = TRUE, stderr = TRUE)
        )
        exit_status <- attr(stdout, "status")
        if (!is.null(exit_status) && exit_status != 0) {
          list(success = FALSE, error = paste(stdout, collapse = "\n"))
        } else {
          container_id <- trimws(stdout[length(stdout)])
          list(success = TRUE, container_id = container_id)
        }
      }, error = function(e) {
        list(success = FALSE, error = e$message)
      })

      if (!submit_result$success) {
        showNotification(paste("Docker launch failed:", submit_result$error),
          type = "error", duration = 15)
        return()
      }

      job_id <- container_name

      # Create Docker job entry
      job_entry <- list(
        job_id = job_id,
        container_id = submit_result$container_id,
        backend = "docker",
        name = analysis_name,
        status = "running",
        output_dir = output_dir,
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
          docker_image = img,
          speclib = if (!is.null(values$diann_speclib) && nzchar(values$diann_speclib))
            basename(values$diann_speclib) else NULL,
          docker = list(cpus = cpus, mem_gb = mem_gb, image = img)
        ),
        auto_load = input$auto_load_results,
        log_content = "",
        completed_at = NULL,
        loaded = FALSE,
        is_ssh = FALSE
      )

    } else {
      # --- HPC (SLURM) submission ---
      sif_path <- input$diann_sif_path

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

      # Create HPC job entry
      job_entry <- list(
        job_id = job_id,
        backend = "hpc",
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
    }

    # --- Shared: add to queue & notify ---
    values$diann_jobs <- c(values$diann_jobs, list(job_entry))

    # Record in SQLite if core facility mode is active
    if (is_core_facility && !is.null(cf_config)) {
      tryCatch({
        cf_record_search(cf_config$db_path, list(
          analysis_name = analysis_name,
          submitted_by  = input$staff_selector %||% "unknown",
          lab           = input$search_lab %||% "",
          instrument    = input$search_instrument %||% "",
          organism      = input$diann_organism %||% "",
          fasta_file    = if (length(values$diann_fasta_files) > 0)
                            basename(values$diann_fasta_files[1]) else "",
          n_raw_files   = if (!is.null(values$diann_raw_files))
                            nrow(values$diann_raw_files) else 0L,
          search_mode   = input$search_mode %||% "libfree",
          slurm_job_id  = if (backend == "hpc") job_id else NA,
          container_id  = if (backend == "docker") job_id else NA,
          backend       = backend,
          output_dir    = output_dir
        ))
      }, error = function(e) {
        message("Core facility DB recording failed: ", e$message)
      })
    }

    add_to_log("DIA-NN Search Submitted", c(
      sprintf("# Job ID: %s", job_id),
      sprintf("# Backend: %s", backend),
      sprintf("# Analysis: %s", analysis_name),
      sprintf("# Files: %d raw data files", nrow(values$diann_raw_files)),
      sprintf("# Mode: %s", input$search_mode),
      sprintf("# Output: %s", output_dir)
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

      if (isTRUE(jobs[[i]]$backend == "local")) {
        # --- Local (embedded) monitoring via processx ---
        proc <- jobs[[i]]$process
        log_path <- jobs[[i]]$log_file

        if (!is.null(proc) && inherits(proc, "process")) {
          result <- check_local_diann_status(proc, log_path)
          new_status <- result$status
          if (nzchar(result$log_tail)) {
            jobs[[i]]$log_content <- result$log_tail
            changed <- TRUE
          }
        } else {
          # Process handle lost (e.g., app restart) — check log file for completion markers
          new_status <- "unknown"
          if (!is.null(log_path) && file.exists(log_path)) {
            log_lines <- tryCatch(readLines(log_path, warn = FALSE), error = function(e) character(0))
            if (any(grepl("Processing finished|report.*saved", log_lines, ignore.case = TRUE))) {
              new_status <- "completed"
            }
            jobs[[i]]$log_content <- paste(tail(log_lines, 30), collapse = "\n")
            changed <- TRUE
          }
        }

      } else if (isTRUE(jobs[[i]]$backend == "docker")) {
        # --- Docker monitoring ---
        cid <- jobs[[i]]$container_id %||% jobs[[i]]$job_id
        result <- check_docker_container_status(cid)
        new_status <- result$status

        if (nzchar(result$log_tail)) {
          jobs[[i]]$log_content <- result$log_tail
          changed <- TRUE
        }
      } else {
        # --- HPC (SLURM) monitoring ---
        job_cfg <- if (isTRUE(jobs[[i]]$is_ssh)) cfg else NULL
        new_status <- check_slurm_status(jobs[[i]]$job_id, ssh_config = job_cfg,
                                          sbatch_path = values$ssh_sbatch_path)

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

      if (new_status != jobs[[i]]$status) {
        jobs[[i]]$status <- new_status
        changed <- TRUE

        if (new_status == "completed") {
          jobs[[i]]$completed_at <- Sys.time()
          showNotification(
            sprintf("DIA-NN search '%s' completed!", jobs[[i]]$name),
            type = "message", duration = 15
          )
          # Docker cleanup: remove stopped container
          if (isTRUE(jobs[[i]]$backend == "docker")) {
            cid <- jobs[[i]]$container_id %||% jobs[[i]]$job_id
            tryCatch(system2("docker", c("rm", cid),
              stdout = FALSE, stderr = FALSE), error = function(e) NULL)
          }
        } else if (new_status == "failed") {
          showNotification(
            sprintf("DIA-NN search '%s' failed. Check log for details.", jobs[[i]]$name),
            type = "error", duration = 15
          )
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
          raw_data <- suppressMessages(suppressWarnings(
            limpa::readDIANN(report_path, format = "parquet")))

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

      backend_icon <- if (isTRUE(job$backend == "local")) {
        span(class = "badge bg-success text-white", style = "font-size: 0.7em; margin-right: 4px;",
          icon("microchip"), " Local")
      } else if (isTRUE(job$backend == "docker")) {
        span(class = "badge bg-info text-white", style = "font-size: 0.7em; margin-right: 4px;",
          icon("docker", lib = "font-awesome"), " Docker")
      } else {
        span(class = "badge bg-secondary", style = "font-size: 0.7em; margin-right: 4px;",
          icon("server"), " HPC")
      }

      div(style = "border: 1px solid #dee2e6; border-radius: 5px; padding: 8px; margin-bottom: 8px; font-size: 0.82em;",
        div(style = "display: flex; justify-content: space-between; align-items: center;",
          div(
            backend_icon,
            strong(job$name), " ",
            span(style = "color: #999;", sprintf("(#%s)", substr(job$job_id, 1, 16)))
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
          safe_log <- iconv(job$log_content %||% "", from = "", to = "UTF-8", sub = "")
          showModal(modalDialog(
            title = sprintf("Log: %s (#%s)", job$name, job$job_id),
            size = "l", easyClose = TRUE, footer = modalButton("Close"),
            pre(style = "max-height: 500px; overflow-y: auto; font-size: 0.8em;",
              safe_log
            )
          ))
        }, ignoreInit = TRUE)

        # Refresh job status
        observeEvent(input[[sprintf("refresh_job_%d", idx)]], {
          job <- values$diann_jobs[[idx]]
          tryCatch({
            if (isTRUE(job$backend == "local")) {
              proc <- job$process
              log_path <- job$log_file
              if (!is.null(proc) && inherits(proc, "process")) {
                result <- check_local_diann_status(proc, log_path)
                new_status <- result$status
              } else {
                # Process handle lost — check log
                new_status <- "unknown"
                if (!is.null(log_path) && file.exists(log_path)) {
                  log_lines <- tryCatch(readLines(log_path, warn = FALSE), error = function(e) character(0))
                  if (any(grepl("Processing finished|report.*saved", log_lines, ignore.case = TRUE))) {
                    new_status <- "completed"
                  }
                }
              }
            } else if (isTRUE(job$backend == "docker")) {
              cid <- job$container_id %||% job$job_id
              result <- check_docker_container_status(cid)
              new_status <- result$status
            } else {
              job_cfg <- if (isTRUE(job$is_ssh)) isolate(ssh_config()) else NULL
              new_status <- check_slurm_status(job$job_id, ssh_config = job_cfg,
                                                sbatch_path = values$ssh_sbatch_path)
            }
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
            if (isTRUE(job$backend == "local")) {
              # Local: kill processx process
              proc <- job$process
              if (!is.null(proc) && inherits(proc, "process") && proc$is_alive()) {
                proc$kill()
              }
            } else if (isTRUE(job$backend == "docker")) {
              # Docker: stop + remove container
              cid <- job$container_id %||% job$job_id
              system2("docker", c("stop", cid), stdout = TRUE, stderr = TRUE)
              tryCatch(system2("docker", c("rm", cid),
                stdout = FALSE, stderr = FALSE), error = function(e) NULL)
            } else if (isTRUE(job$is_ssh)) {
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
          report_path <- NULL

          tryCatch({
            if (isTRUE(job$is_ssh)) {
              # SSH mode: SCP download first
              cfg <- ssh_config()
              if (is.null(cfg)) {
                showNotification("SSH not configured. Test connection first.", type = "error", duration = 8)
                return()
              }

              if (!nzchar(job$output_dir %||% "")) {
                showNotification("No output directory known for this job. Try Recover first.", type = "error", duration = 8)
                return()
              }

              showNotification("Locating report on remote...", type = "message", duration = 3, id = "load_progress")

              # Find the report.parquet file on remote
              remote_report <- file.path(job$output_dir, "report.parquet")
              find_result <- ssh_exec(cfg, paste("ls", shQuote(remote_report), "2>/dev/null"))
              if (find_result$status != 0) {
                find_result <- ssh_exec(cfg, sprintf(
                  "ls %s/report*.parquet 2>/dev/null | head -1", shQuote(job$output_dir)))
                if (find_result$status != 0 || length(find_result$stdout) == 0 ||
                    !nzchar(trimws(find_result$stdout[1]))) {
                  showNotification("No report.parquet found on remote.", type = "error", duration = 8)
                  return()
                }
                remote_report <- trimws(find_result$stdout[1])
              }

              # Download via SCP
              showNotification("Downloading report via SCP...", type = "message", duration = 30, id = "load_progress")
              local_report <- file.path(tempdir(), paste0(job$name, "_", basename(remote_report)))
              dl_result <- scp_download(cfg, remote_report, local_report)
              if (dl_result$status != 0) {
                showNotification("SCP download failed.", type = "error", duration = 8)
                return()
              }
              report_path <- local_report

            } else {
              # Local mode: direct access
              report_path <- file.path(job$output_dir, "report.parquet")
              if (!file.exists(report_path)) {
                parquet_files <- list.files(job$output_dir, pattern = "report.*\\.parquet$",
                  full.names = TRUE)
                if (length(parquet_files) > 0) {
                  report_path <- parquet_files[1]
                } else {
                  showNotification("No report.parquet found in output directory.", type = "error", duration = 8)
                  return()
                }
              }
            }

            if (is.null(report_path) || !file.exists(report_path)) {
              showNotification("Report file not available.", type = "error", duration = 8)
              return()
            }

            showNotification("Reading DIA-NN report...", type = "message", duration = 30, id = "load_progress")
            raw_data <- suppressMessages(suppressWarnings(
              limpa::readDIANN(report_path, format = "parquet")))
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

            removeNotification(id = "load_progress")
            nav_select("main_tabs", "Data Overview")
            nav_select("data_overview_tabs", "Assign Groups & Run")

            showNotification("Results loaded! Assign groups and run pipeline.",
              type = "message", duration = 8)
          }, error = function(e) {
            removeNotification(id = "load_progress")
            err_msg <- tryCatch(
              iconv(conditionMessage(e), from = "", to = "UTF-8", sub = ""),
              error = function(e2) "Unknown error (possible encoding issue)"
            )
            showNotification(paste("Failed to load:", err_msg), type = "error", duration = 10)
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
        if (isTRUE(jobs[[i]]$backend == "local")) {
          proc <- jobs[[i]]$process
          log_path <- jobs[[i]]$log_file
          if (!is.null(proc) && inherits(proc, "process")) {
            result <- check_local_diann_status(proc, log_path)
            new_status <- result$status
          } else {
            new_status <- "unknown"
            if (!is.null(log_path) && file.exists(log_path)) {
              log_lines <- tryCatch(readLines(log_path, warn = FALSE), error = function(e) character(0))
              if (any(grepl("Processing finished|report.*saved", log_lines, ignore.case = TRUE))) {
                new_status <- "completed"
              }
            }
          }
        } else if (isTRUE(jobs[[i]]$backend == "docker")) {
          cid <- jobs[[i]]$container_id %||% jobs[[i]]$job_id
          result <- check_docker_container_status(cid)
          new_status <- result$status
        } else {
          job_cfg <- if (isTRUE(jobs[[i]]$is_ssh)) cfg else NULL
          new_status <- check_slurm_status(jobs[[i]]$job_id, ssh_config = job_cfg,
                                            sbatch_path = values$ssh_sbatch_path)
        }
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

  # ============================================================================
  #    Recover Jobs from SLURM / Docker
  # ============================================================================

  observeEvent(input$recover_jobs_btn, {
    recovered <- 0
    updated <- 0

    # --- Recover HPC jobs via sacct ---
    if (hpc_available) {
      cfg <- isolate(ssh_config())
      withProgress(message = "Scanning SLURM for previous DIA-NN jobs...", {
        slurm_jobs <- recover_slurm_jobs(
          ssh_config = cfg,
          sbatch_path = values$ssh_sbatch_path,
          days_back = 14
        )
      })

      if (nrow(slurm_jobs) > 0) {
        existing_ids <- vapply(values$diann_jobs, function(j) j$job_id, character(1))

        for (i in seq_len(nrow(slurm_jobs))) {
          row <- slurm_jobs[i, ]

          # Map SLURM state to DE-LIMP status
          status <- switch(toupper(row$state),
            "COMPLETED" = "completed",
            "RUNNING"   = "running",
            "PENDING"   = "queued",
            "FAILED"    = "failed",
            "CANCELLED" = "cancelled",
            "TIMEOUT"   = "failed",
            "unknown"
          )

          # Find the actual log file and output directory.
          # Strategy 1: scontrol show job → StdOut path (most reliable for recent jobs)
          # Strategy 2: sacct SubmitLine → script path → derive output dir
          # Strategy 3: find in common HPC paths
          output_dir <- ""
          log_content <- ""
          n_files <- 0
          log_file <- ""

          run_ssh <- function(cmd) {
            if (!is.null(cfg)) ssh_exec(cfg, cmd, timeout = 30)
            else {
              out <- tryCatch(system2("bash", c("-c", cmd), stdout = TRUE, stderr = TRUE),
                error = function(e) character())
              list(status = 0, stdout = out)
            }
          }

          # Derive SLURM tool paths from sbatch path
          slurm_bin_dir <- if (!is.null(values$ssh_sbatch_path) && nzchar(values$ssh_sbatch_path)) {
            dirname(values$ssh_sbatch_path)
          } else ""

          scontrol_bin <- if (nzchar(slurm_bin_dir)) file.path(slurm_bin_dir, "scontrol") else "scontrol"
          sacct_bin <- if (nzchar(slurm_bin_dir)) file.path(slurm_bin_dir, "sacct") else "sacct"

          # Strategy 1: scontrol show job → extract StdOut path
          # Use sed instead of grep -oP for portability (not all systems have PCRE grep)
          scontrol_result <- tryCatch({
            run_ssh(sprintf(
              "%s show job %s 2>/dev/null | sed -n 's/.*StdOut=//p' | tr -d ' '",
              scontrol_bin, row$job_id))
          }, error = function(e) list(status = 1, stdout = character()))

          if (scontrol_result$status == 0 && length(scontrol_result$stdout) > 0 &&
              nzchar(trimws(scontrol_result$stdout[1]))) {
            log_file <- trimws(scontrol_result$stdout[1])
            output_dir <- dirname(log_file)
          }

          # Strategy 2: sacct SubmitLine → derive from script path
          if (!nzchar(log_file)) {
            submit_result <- tryCatch({
              run_ssh(sprintf(
                "%s -j %s --format=SubmitLine%%300 --parsable2 --noheader 2>/dev/null | head -1",
                sacct_bin, row$job_id))
            }, error = function(e) list(status = 1, stdout = character()))

            if (submit_result$status == 0 && length(submit_result$stdout) > 0 &&
                nzchar(trimws(submit_result$stdout[1]))) {
              submit_line <- trimws(submit_result$stdout[1])
              parts <- strsplit(submit_line, "[[:space:]]+")[[1]]
              script_path <- parts[grepl("/.*\\.sbatch$", parts)]
              if (length(script_path) > 0) {
                output_dir <- dirname(script_path[1])
                log_file <- file.path(output_dir, sprintf("diann_%s.out", row$job_id))
              }
            }
          }

          # Strategy 3: search configured output base + common HPC paths
          # Use timeout to avoid long waits on large shared filesystems
          if (!nzchar(log_file)) {
            search_base <- isolate(output_base())
            find_result <- tryCatch({
              find_cmd <- sprintf(paste0(
                "timeout 10 find %s -maxdepth 4 -name 'diann_%s.out' 2>/dev/null | head -1"),
                shQuote(search_base), row$job_id)
              if (!is.null(cfg)) ssh_exec(cfg, find_cmd, timeout = 15)
              else {
                out <- system2("bash", c("-c", find_cmd), stdout = TRUE, stderr = TRUE)
                list(status = 0, stdout = out)
              }
            }, error = function(e) list(status = 1, stdout = character()))

            if (find_result$status == 0 && length(find_result$stdout) > 0 &&
                nzchar(trimws(find_result$stdout[1]))) {
              log_file <- trimws(find_result$stdout[1])
              output_dir <- dirname(log_file)
            }
          }

          # Fetch actual log content and file count
          if (nzchar(log_file)) {
            # Get file count from the "N files will be processed" line near the top
            count_result <- tryCatch(
              run_ssh(sprintf(
                "grep -m1 'files will be processed' %s 2>/dev/null", shQuote(log_file))),
              error = function(e) list(status = 1, stdout = character()))

            if (count_result$status == 0 && length(count_result$stdout) > 0 &&
                nzchar(count_result$stdout[1])) {
              # Line format: "[HH:MM] N files will be processed"
              m <- regexpr("[0-9]+(?=\\s+files will be processed)",
                count_result$stdout[1], perl = TRUE)
              if (m > 0) n_files <- as.integer(regmatches(count_result$stdout[1], m))
            }

            # Tail the log for display
            tail_result <- tryCatch(
              run_ssh(sprintf("tail -150 %s 2>/dev/null", shQuote(log_file))),
              error = function(e) list(status = 1, stdout = character()))

            if (tail_result$status == 0 && length(tail_result$stdout) > 0) {
              log_content <- iconv(paste(tail_result$stdout, collapse = "\n"),
                from = "", to = "UTF-8", sub = "")
            }
          }

          if (!nzchar(log_content)) {
            log_content <- sprintf(paste0(
              "Recovered from SLURM sacct.\nState: %s, Elapsed: %s\n",
              "Output dir: %s\n\n",
              "Could not locate log file diann_%s.out on the cluster.\n",
              "Tried: scontrol show job, sacct SubmitLine, find in common paths."),
              row$state, row$elapsed,
              if (nzchar(output_dir)) output_dir else "(unknown)", row$job_id)
          }

          # Check if job already exists in queue
          existing_idx <- match(row$job_id, existing_ids)

          if (!is.na(existing_idx)) {
            # Update existing entry with fresh data from cluster
            jobs <- values$diann_jobs
            jobs[[existing_idx]]$status <- status
            jobs[[existing_idx]]$log_content <- log_content
            if (nzchar(output_dir)) jobs[[existing_idx]]$output_dir <- output_dir
            if (n_files > 0) jobs[[existing_idx]]$n_files <- n_files
            if (status %in% c("completed", "failed", "cancelled") &&
                is.null(jobs[[existing_idx]]$completed_at)) {
              jobs[[existing_idx]]$completed_at <- Sys.time()
            }
            values$diann_jobs <- jobs
            updated <- updated + 1
          } else {
            # Add new entry
            job_entry <- list(
              job_id = row$job_id,
              backend = "hpc",
              name = row$name,
              status = status,
              output_dir = output_dir,
              submitted_at = Sys.time(),
              n_files = n_files,
              search_mode = "unknown",
              search_settings = NULL,
              auto_load = FALSE,
              log_content = log_content,
              completed_at = if (status %in% c("completed", "failed", "cancelled")) Sys.time() else NULL,
              loaded = FALSE,
              is_ssh = !is.null(cfg)
            )
            values$diann_jobs <- c(values$diann_jobs, list(job_entry))
            recovered <- recovered + 1
          }
        }
      }
    }

    # --- Recover Docker jobs ---
    if (docker_available) {
      withProgress(message = "Scanning Docker for previous DIA-NN containers...", {
        docker_jobs <- recover_docker_jobs()
      })

      if (nrow(docker_jobs) > 0) {
        existing_ids <- vapply(values$diann_jobs, function(j) j$job_id, character(1))
        for (i in seq_len(nrow(docker_jobs))) {
          row <- docker_jobs[i, ]

          # Check actual container status
          result <- check_docker_container_status(row$container_id)

          existing_idx <- match(row$name, existing_ids)

          if (!is.na(existing_idx)) {
            jobs <- values$diann_jobs
            jobs[[existing_idx]]$status <- result$status
            jobs[[existing_idx]]$log_content <- result$log_tail
            values$diann_jobs <- jobs
            updated <- updated + 1
          } else {
            job_entry <- list(
              job_id = row$name,
              container_id = row$container_id,
              backend = "docker",
              name = sub("^delimp_", "", row$name),
              status = result$status,
              output_dir = "",
              submitted_at = Sys.time(),
              n_files = 0,
              search_mode = "unknown",
              search_settings = NULL,
              auto_load = FALSE,
              log_content = result$log_tail,
              completed_at = if (result$status %in% c("completed", "failed")) Sys.time() else NULL,
              loaded = FALSE,
              is_ssh = FALSE
            )
            values$diann_jobs <- c(values$diann_jobs, list(job_entry))
            recovered <- recovered + 1
          }
        }
      }
    }

    if (recovered > 0 || updated > 0) {
      parts <- c()
      if (recovered > 0) parts <- c(parts, sprintf("%d new job(s) recovered", recovered))
      if (updated > 0) parts <- c(parts, sprintf("%d existing job(s) updated", updated))
      showNotification(paste(parts, collapse = ", "), type = "message", duration = 8)
    } else {
      showNotification("No DIA-NN jobs found on cluster.", type = "message", duration = 5)
    }
  })

}
