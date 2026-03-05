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
                          is_core_facility = FALSE, cf_config = NULL,
                          local_sbatch_path = "") {

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
  #    Parallel Search Mode UI (HPC backend, >= 8 files)
  # ============================================================================

  output$parallel_mode_ui <- renderUI({
    # Only show for HPC backend with >= 8 files
    req(input$search_backend == "hpc")
    n_files <- if (!is.null(values$diann_raw_files)) nrow(values$diann_raw_files) else 0
    if (n_files < 8) return(NULL)

    # Recommendation badge based on file count
    rec_badge <- if (n_files >= 50) {
      span(class = "badge bg-danger", "Strongly recommended")
    } else if (n_files >= 20) {
      span(class = "badge bg-warning text-dark", "Recommended")
    } else {
      span(class = "badge bg-info", "Optional")
    }

    tagList(
      hr(),
      tags$h6(icon("layer-group"), " Parallel Search Mode ", rec_badge),
      checkboxInput("parallel_search", "Enable Parallel Search (split across nodes)",
        value = FALSE),
      tags$p(class = "text-muted", style = "font-size: 0.8em;",
        sprintf("With %d files, parallel search splits processing into 5 steps: ", n_files),
        "library prediction, per-file first-pass, library assembly, ",
        "per-file final-pass, and cross-run report. Each file runs as a ",
        "separate SLURM array task."),

      conditionalPanel("input.parallel_search",
        tags$div(class = "alert alert-info py-1 px-2",
          style = "font-size: 0.82em;",
          icon("circle-info"),
          " Mass accuracy is forced to Manual mode in parallel search. ",
          "MBR is disabled (replaced by the 5-step workflow)."),

        tags$h6("Per-File Resources", style = "margin-top: 10px;"),
        div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
          div(style = "flex: 1; min-width: 100px;",
            numericInput("parallel_cpus", "CPUs/file:", value = 16, min = 4, max = 32, step = 4)
          ),
          div(style = "flex: 1; min-width: 100px;",
            numericInput("parallel_mem_gb", "Memory/file (GB):", value = 32, min = 8, max = 128, step = 8)
          )
        ),
        div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
          div(style = "flex: 1; min-width: 100px;",
            numericInput("parallel_time_hours", "Time/file (hrs):", value = 2, min = 1, max = 8, step = 1)
          ),
          div(style = "flex: 1; min-width: 100px;",
            numericInput("max_simultaneous", "Max concurrent:", value = 20, min = 4, max = 100, step = 4)
          )
        ),
        tags$p(class = "text-muted", style = "font-size: 0.78em;",
          "Assembly/report steps use the main SLURM resource settings above. ",
          "Per-file resources apply to the array jobs (steps 2 & 4).")
      )
    )
  })

  # Instrument-aware mass accuracy hint
  output$mass_acc_hint <- renderUI({
    ext <- NULL
    if (!is.null(values$diann_raw_files) && nrow(values$diann_raw_files) > 0) {
      ext <- tolower(tools::file_ext(values$diann_raw_files$filename[1]))
    }
    hint <- if (identical(ext, "d")) {
      "timsTOF detected \u2014 recommended: MS2 15, MS1 15"
    } else if (identical(ext, "raw")) {
      "Orbitrap detected \u2014 recommended: MS2 10, MS1 5"
    } else if (identical(ext, "mzml")) {
      "Instrument unknown (.mzML) \u2014 typical: MS2 10\u201315, MS1 5\u201315"
    } else {
      "These values are passed directly to DIA-NN"
    }
    tags$p(class = "text-muted", style = "font-size: 0.78em; margin-top: -4px;", hint)
  })

  # Force mass accuracy to manual when parallel mode is enabled
  observeEvent(input$parallel_search, {
    if (isTRUE(input$parallel_search)) {
      updateRadioButtons(session, "mass_acc_mode", selected = "manual")

      # Set instrument-aware defaults based on file extensions
      if (!is.null(values$diann_raw_files) && nrow(values$diann_raw_files) > 0) {
        ext <- tolower(tools::file_ext(values$diann_raw_files$filename[1]))
        if (ext == "d") {
          # timsTOF
          updateNumericInput(session, "diann_mass_acc", value = 15)
          updateNumericInput(session, "diann_mass_acc_ms1", value = 15)
        } else if (ext == "raw") {
          # Orbitrap
          updateNumericInput(session, "diann_mass_acc", value = 10)
          updateNumericInput(session, "diann_mass_acc_ms1", value = 5)
        }
        # .mzML: keep current values (can't detect instrument)
      }

      showNotification(
        "Parallel mode: mass accuracy set to Manual with instrument-aware defaults. MBR disabled.",
        type = "message", duration = 6)
    }
  }, ignoreInit = TRUE)

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
      # Exclude removed jobs from persistence to avoid unbounded growth
      active_jobs <- Filter(function(j) !isTRUE(j$removed), values$diann_jobs)
      saveRDS(active_jobs, job_queue_path)
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

    # Trigger immediate cluster resource check + auto-select on successful connect
    if (result$success) {
      tryCatch({
        res <- check_cluster_resources(
          ssh_config = cfg, account = "genome-center-grp",
          partition = "high", sbatch_path = result$sbatch_path
        )
        values$cluster_resources <- res
      }, error = function(e) {
        values$cluster_resources <- list(success = FALSE, error = e$message)
      })
      tryCatch({
        pub_res <- check_cluster_resources(
          ssh_config = cfg, account = "publicgrp",
          partition = "low", sbatch_path = result$sbatch_path
        )
        values$public_resources <- pub_res
      }, error = function(e) NULL)

      # Auto-select best partition immediately
      best <- select_best_partition(values$cluster_resources, values$public_resources, 64)
      values$auto_partition <- best
      if (!isTRUE(isolate(input$partition_override))) {
        updateTextInput(session, "diann_account", value = best$account)
        updateTextInput(session, "diann_partition", value = best$partition)
      }
    } else {
      values$cluster_resources <- NULL
      values$public_resources <- NULL
    }
  })

  # ============================================================================
  #    Auto-connect SSH on startup (if credentials pre-filled)
  # ============================================================================

  observe({
    # Wait for inputs to initialize
    req(input$ssh_host, input$ssh_user, input$ssh_key_path)
    req(!isTRUE(values$ssh_connected))
    cfg <- ssh_config()
    req(cfg)

    message("[DE-LIMP] Auto-connecting SSH to ", cfg$host, "...")
    result <- test_ssh_connection(cfg)

    output$ssh_status_ui <- renderUI({
      if (result$success) {
        div(class = "alert alert-success py-1 px-2 mt-2",
          style = "font-size: 0.82em;",
          icon("check-circle"), " ", result$message)
      } else {
        div(class = "alert alert-warning py-1 px-2 mt-2",
          style = "font-size: 0.82em;",
          icon("info-circle"), " Auto-connect failed. Click Test Connection to retry.")
      }
    })

    values$ssh_connected <- result$success
    values$ssh_sbatch_path <- result$sbatch_path

    if (result$success) {
      tryCatch({
        res <- check_cluster_resources(
          ssh_config = cfg, account = "genome-center-grp",
          partition = "high", sbatch_path = result$sbatch_path)
        values$cluster_resources <- res
      }, error = function(e) {
        values$cluster_resources <- list(success = FALSE, error = e$message)
      })
      tryCatch({
        pub_res <- check_cluster_resources(
          ssh_config = cfg, account = "publicgrp",
          partition = "low", sbatch_path = result$sbatch_path)
        values$public_resources <- pub_res
      }, error = function(e) NULL)

      best <- select_best_partition(values$cluster_resources, values$public_resources, 64)
      values$auto_partition <- best
      if (!isTRUE(isolate(input$partition_override))) {
        updateTextInput(session, "diann_account", value = best$account)
        updateTextInput(session, "diann_partition", value = best$partition)
      }
    }
  }) |> bindEvent(input$ssh_host, once = TRUE)

  # ============================================================================
  #    Re-verify stale job statuses on SSH connect
  # ============================================================================
  #
  # Jobs saved as "completed" or "running" in RDS may have stale statuses
  # (e.g., a FAILED job showing as completed due to the .extern sacct bug,
  # or a running job that finished while the app was closed). Re-check once
  # when SSH connects.

  jobs_reverified <- reactiveVal(FALSE)

  observe({
    req(isTRUE(values$ssh_connected))
    req(!jobs_reverified())
    req(length(values$diann_jobs) > 0)
    message("[DE-LIMP] Re-verifying ", length(values$diann_jobs), " jobs after SSH connect...")

    jobs <- values$diann_jobs
    cfg <- isolate(ssh_config())
    slurm_path <- isolate(values$ssh_sbatch_path)
    changed <- FALSE
    n_updated <- 0

    for (i in seq_along(jobs)) {
      if (isTRUE(jobs[[i]]$removed)) next
      if (!isTRUE(jobs[[i]]$is_ssh)) next
      if (!jobs[[i]]$status %in% c("completed", "running", "queued")) next

      tryCatch({
        new_status <- check_slurm_status(
          jobs[[i]]$job_id, ssh_config = cfg, sbatch_path = slurm_path)

        # For "completed" jobs, also verify report.parquet exists on remote.
        # DIA-NN can hit internal errors (e.g. library mismatch) and exit 0
        # without producing output — SLURM says COMPLETED but there's no report.
        if (identical(new_status, "completed") && !isTRUE(jobs[[i]]$loaded)) {
          out_dir <- jobs[[i]]$output_dir
          if (!is.null(out_dir) && nzchar(out_dir)) {
            report_check <- ssh_exec(cfg,
              sprintf("test -f %s && echo EXISTS || echo MISSING",
                shQuote(file.path(out_dir, "report.parquet"))))
            if (report_check$status == 0 &&
                any(grepl("MISSING", report_check$stdout))) {
              message(sprintf("[DE-LIMP] Job %s: SLURM says completed but no report.parquet — marking failed",
                jobs[[i]]$job_id))
              new_status <- "failed"
              jobs[[i]]$failure_reason <- "DIA-NN completed without producing report.parquet (check logs)"
            }
          }
        }

        if (!is.null(new_status) && new_status != jobs[[i]]$status) {
          message(sprintf("[DE-LIMP] Job %s status corrected: %s -> %s",
            jobs[[i]]$job_id, jobs[[i]]$status, new_status))
          jobs[[i]]$status <- new_status
          if (new_status == "completed" && is.null(jobs[[i]]$completed_at)) {
            jobs[[i]]$completed_at <- Sys.time()
          }
          changed <- TRUE
          n_updated <- n_updated + 1
        }
      }, error = function(e) {
        message("[DE-LIMP] Re-verify failed for job ", jobs[[i]]$job_id, ": ", e$message)
      })
    }

    if (changed) {
      values$diann_jobs <- jobs
      showNotification(
        sprintf("Re-verified job statuses: %d updated", n_updated),
        type = "message", duration = 5)
    }
    jobs_reverified(TRUE)
  })

  # ============================================================================
  #    Cluster Resource Indicator (auto-refresh every 60s)
  # ============================================================================

  observe({
    invalidateLater(60000)
    req(isTRUE(values$ssh_connected))
    cfg <- isolate(ssh_config())
    req(cfg)
    sbatch_path <- isolate(values$ssh_sbatch_path)

    # Always check both accounts
    tryCatch({
      res <- check_cluster_resources(cfg, "genome-center-grp", "high", sbatch_path)
      values$cluster_resources <- res
    }, error = function(e) NULL)

    tryCatch({
      pub_res <- check_cluster_resources(cfg, "publicgrp", "low", sbatch_path)
      values$public_resources <- pub_res
    }, error = function(e) NULL)

    # Auto-select best partition
    peak_cpus <- if (isTRUE(isolate(input$parallel_search))) {
      cpus_per <- isolate(input$parallel_cpus) %||% 16
      max_sim <- isolate(input$max_simultaneous) %||% 20
      max(32, cpus_per * max_sim)
    } else {
      isolate(input$diann_cpus) %||% 64
    }

    best <- select_best_partition(values$cluster_resources, values$public_resources, peak_cpus)
    values$auto_partition <- best

    # Update hidden inputs unless user is overriding
    if (!isTRUE(isolate(input$partition_override))) {
      updateTextInput(session, "diann_account", value = best$account)
      updateTextInput(session, "diann_partition", value = best$partition)
    }
  })

  output$cluster_status_ui <- renderUI({
    res <- values$cluster_resources
    if (is.null(res) || !isTRUE(res$success)) return(NULL)

    has_group <- !is.na(res$group_limit)
    has_partition <- !is.na(res$partition_idle)

    if (!has_group && !has_partition) return(NULL)

    # Determine traffic light color from group utilization
    if (has_group && res$group_limit > 0) {
      pct_used <- res$group_used / res$group_limit
      if (pct_used > 0.8) {
        color <- "#dc3545"; bg <- "#f8d7da"; border <- "#f5c2c7"
      } else if (pct_used > 0.5) {
        color <- "#ffc107"; bg <- "#fff3cd"; border <- "#ffecb5"
      } else {
        color <- "#198754"; bg <- "#d1e7dd"; border <- "#badbcc"
      }
    } else if (has_partition && res$partition_total > 0) {
      pct_idle <- res$partition_idle / res$partition_total
      if (pct_idle < 0.2) {
        color <- "#dc3545"; bg <- "#f8d7da"; border <- "#f5c2c7"
      } else if (pct_idle < 0.5) {
        color <- "#ffc107"; bg <- "#fff3cd"; border <- "#ffecb5"
      } else {
        color <- "#198754"; bg <- "#d1e7dd"; border <- "#badbcc"
      }
    } else {
      color <- "#6c757d"; bg <- "#e9ecef"; border <- "#dee2e6"
    }

    # Build text lines
    group_line <- if (has_group) {
      sprintf("genome-center-grp: %s/%s CPUs used (%s available)",
        format(res$group_used, big.mark = ","),
        format(res$group_limit, big.mark = ","),
        format(res$group_available, big.mark = ","))
    } else if (!is.na(res$group_used)) {
      sprintf("genome-center-grp: %s CPUs in use",
        format(res$group_used, big.mark = ","))
    } else NULL

    pub_res <- values$public_resources
    pub_line <- if (!is.null(pub_res) && isTRUE(pub_res$success) && !is.na(pub_res$partition_idle)) {
      sprintf("publicgrp/low: %s idle of %s total",
        format(pub_res$partition_idle, big.mark = ","),
        format(pub_res$partition_total, big.mark = ","))
    } else NULL

    indicator <- span(style = sprintf("color: %s; font-size: 1.1em;", color),
      HTML("&#9679;"))

    div(
      style = sprintf(
        "background: %s; border: 1px solid %s; border-radius: 6px; padding: 6px 10px; margin-top: 8px; font-size: 0.82em;",
        bg, border),
      div(indicator, " ", if (!is.null(group_line)) group_line),
      if (!is.null(pub_line)) div(
        style = "margin-left: 20px; color: #555;", pub_line)
    )
  })

  # ============================================================================
  #    Auto-Select Partition UI + Override
  # ============================================================================

  output$partition_selector_ui <- renderUI({
    selected <- values$auto_partition
    override <- isTRUE(input$partition_override)

    tagList(
      # Auto-selected display (when not overriding)
      if (!override && !is.null(selected)) {
        div(style = "background: #d1e7dd; border: 1px solid #badbcc; border-radius: 4px; padding: 6px 10px; margin-bottom: 6px; font-size: 0.85em;",
          icon("magic"), " ",
          tags$strong(sprintf("%s / %s", selected$account, selected$partition)),
          span(style = "color: #555; margin-left: 6px;", sprintf("(%s)", selected$reason))
        )
      },
      # Override toggle
      checkboxInput("partition_override", "Override account/partition", value = override),
      # Manual inputs (shown only when override checked)
      conditionalPanel("input.partition_override",
        div(style = "display: flex; gap: 8px;",
          div(style = "flex: 1;", textInput("diann_account_override", "Account:", value = "genome-center-grp")),
          div(style = "flex: 1;", textInput("diann_partition_override", "Partition:", value = "high"))
        )
      )
    )
  })

  # Override toggle: sync hidden inputs from override inputs or restore auto-selected
  observeEvent(input$partition_override, {
    if (isTRUE(input$partition_override)) {
      updateTextInput(session, "diann_account", value = input$diann_account_override %||% "genome-center-grp")
      updateTextInput(session, "diann_partition", value = input$diann_partition_override %||% "high")
    } else {
      best <- values$auto_partition
      if (!is.null(best)) {
        updateTextInput(session, "diann_account", value = best$account)
        updateTextInput(session, "diann_partition", value = best$partition)
      }
    }
  })

  # Sync hidden inputs when override inputs change
  observeEvent(c(input$diann_account_override, input$diann_partition_override), {
    if (isTRUE(input$partition_override)) {
      updateTextInput(session, "diann_account", value = input$diann_account_override)
      updateTextInput(session, "diann_partition", value = input$diann_partition_override)
    }
  }, ignoreInit = TRUE)

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
    vols <- c(Home = Sys.getenv("HOME"), Root = "/")
    # Add proteomics shared storage if available (HPC default)
    if (dir.exists("/quobyte/proteomics-grp")) {
      vols <- c(Proteomics = "/quobyte/proteomics-grp", vols)
    }
    vols
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

  # Open UniProt search modal
  observeEvent(input$open_uniprot_modal, {
    showModal(modalDialog(
      title = tagList(icon("dna"), " UniProt FASTA Database Search"),
      size = "l",
      easyClose = TRUE,
      div(style = "display: flex; gap: 8px; margin-bottom: 12px;",
        div(style = "flex: 1;",
          textInput("uniprot_search_query", NULL,
            placeholder = "e.g., human, mouse, E. coli", width = "100%")
        ),
        actionButton("search_uniprot", "Search",
          class = "btn-info", style = "margin-top: 0;")
      ),
      DTOutput("uniprot_results_table"),
      hr(),
      div(style = "display: flex; gap: 12px; align-items: flex-end;",
        div(style = "flex: 1;",
          selectInput("fasta_content_type", "Content:",
            choices = c(
              "One per gene (recommended)" = "one_per_gene",
              "Swiss-Prot reviewed" = "reviewed",
              "Swiss-Prot + isoforms" = "reviewed_isoforms",
              "Full proteome" = "full",
              "Full + isoforms" = "full_isoforms"
            ), selected = "one_per_gene", width = "100%")
        ),
        div(style = "flex: 1;",
          uiOutput("fasta_filename_preview_modal")
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("download_fasta_btn", "Download FASTA",
          class = "btn-success", icon = icon("download"))
      )
    ))
  })

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
        pageLength = 10, dom = "tip", scrollY = "300px",
        columnDefs = list(list(width = "90px", targets = 0))
      ),
      rownames = FALSE,
      class = "compact stripe"
    )
  })

  # Filename preview inside modal
  output$fasta_filename_preview_modal <- renderUI({
    req(values$uniprot_results, nrow(values$uniprot_results) > 0)
    sel <- input$uniprot_results_table_rows_selected
    req(length(sel) > 0)

    row <- values$uniprot_results[sel, ]
    fname <- generate_fasta_filename(row$upid, row$organism, input$fasta_content_type)

    div(style = "font-size: 0.85em; color: #6c757d; padding-top: 28px;",
      icon("file"), " ", fname
    )
  })

  # Filename preview in sidebar (after download)
  output$fasta_filename_preview <- renderUI({
    req(length(values$diann_fasta_files) > 0, all(nzchar(values$diann_fasta_files)))
    div(style = "font-size: 0.8em; color: #6c757d; margin-top: 5px;",
      icon("check-circle", style = "color: #28a745;"), " ",
      basename(values$diann_fasta_files[1])
    )
  })

  # Summary of selected proteome in sidebar
  output$fasta_selected_summary <- renderUI({
    req(values$uniprot_results, nrow(values$uniprot_results) > 0)
    sel <- input$uniprot_results_table_rows_selected
    req(length(sel) > 0)
    row <- values$uniprot_results[sel, ]

    if (!is.null(values$fasta_info) && !is.null(values$fasta_info$n_sequences)) {
      # FASTA has been downloaded — show actual sequence count
      div(style = "font-size: 0.8em; color: #495057; margin-top: 5px;",
        tags$strong(row$common_name), " \u2014 ",
        format(values$fasta_info$n_sequences, big.mark = ","), " sequences downloaded"
      )
    } else {
      # Not yet downloaded — show organism name only
      div(style = "font-size: 0.8em; color: #495057; margin-top: 5px;",
        tags$strong(row$common_name), " \u2014 proteome ",
        tags$code(row$upid)
      )
    }
  })

  # Clear stale download info when user changes selection or content type
  observeEvent(input$uniprot_results_table_rows_selected, { values$fasta_info <- NULL })
  observeEvent(input$fasta_content_type, { values$fasta_info <- NULL })

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
      # Warn if FTP one-per-gene wasn't available and we fell back to full proteome
      if (!is.null(result$warning)) {
        showNotification(result$warning, type = "warning", duration = 12)
      }
      removeModal()
      cfg <- ssh_config()
      if (!is.null(cfg)) {
        # SSH mode: check if FASTA already exists on remote
        remote_fasta_dir <- file.path(output_base(), "databases")
        remote_path <- file.path(remote_fasta_dir, fname)

        # Check if FASTA already exists on remote AND has matching sequence count
        needs_upload <- TRUE
        exists_check <- ssh_exec(cfg,
          paste("test -f", shQuote(remote_path), "&& grep -c '^>' ", shQuote(remote_path)))
        remote_count <- suppressWarnings(
          as.integer(trimws(paste(exists_check$stdout, collapse = ""))))
        if (!is.na(remote_count) && remote_count == result$n_sequences) {
          needs_upload <- FALSE
          values$diann_fasta_files <- remote_path
          values$fasta_info <- result
          showNotification(
            sprintf("FASTA already exists on HPC (%d sequences): %s",
              remote_count, remote_path),
            type = "message", duration = 8)
        } else if (!is.na(remote_count)) {
          showNotification(
            sprintf("Existing FASTA has %d sequences but download has %d — re-uploading.",
              remote_count, result$n_sequences),
            type = "warning", duration = 8)
        }
        if (needs_upload) {
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
  #    Shared FASTA Database Library
  # ============================================================================

  # Reactive: FASTA library catalog (reloaded when triggered)
  fasta_library_catalog <- reactiveVal(list())

  # Load catalog on startup and when refresh is triggered
  observe({
    fasta_library_catalog(fasta_library_load())
  }) |> bindEvent(TRUE)

  # Open the library modal
  observeEvent(input$open_fasta_library_modal, {
    # Refresh catalog each time modal opens
    fasta_library_catalog(fasta_library_load())

    showModal(modalDialog(
      title = tagList(icon("book"), " Shared FASTA Database Library"),
      size = "xl",
      easyClose = TRUE,
      div(
        # Status banner: shared vs local
        if (fasta_library_is_shared()) {
          div(class = "alert alert-success py-1 px-3 mb-2",
            style = "font-size: 0.85em;",
            icon("network-wired"), " Connected to shared proteomics volume"
          )
        } else {
          div(class = "alert alert-warning py-1 px-3 mb-2",
            style = "font-size: 0.85em;",
            icon("user"), " Using local library (shared volume not mounted)"
          )
        },
        # Catalog table
        DTOutput("fasta_library_table"),
        hr(),
        # Detail panel (shown on row select)
        uiOutput("fasta_library_detail_panel")
      ),
      footer = tagList(
        actionButton("fasta_library_refresh_btn", "Refresh",
          class = "btn-outline-secondary btn-sm", icon = icon("sync")),
        modalButton("Cancel"),
        actionButton("fasta_library_use_btn", "Use This Database",
          class = "btn-success", icon = icon("check"))
      )
    ))
  })

  # Render the library catalog table
  output$fasta_library_table <- DT::renderDT({
    catalog <- fasta_library_catalog()
    display_df <- fasta_library_display_df(catalog)

    if (nrow(display_df) == 0) {
      # Return empty table with proper columns
      return(DT::datatable(
        data.frame(
          Name = character(), Organism = character(), Proteins = character(),
          Age = character(), Status = character(), `Created By` = character(),
          stringsAsFactors = FALSE, check.names = FALSE
        ),
        selection = "single",
        options = list(dom = "t", language = list(
          emptyTable = "No databases in library. Download from UniProt and add to library."
        )),
        rownames = FALSE,
        class = "compact stripe"
      ))
    }

    # Format proteins with comma separator
    display_df$Proteins <- format(display_df$Proteins, big.mark = ",")

    # Color-code Status column
    display_df$Status <- vapply(display_df$Status, function(s) {
      switch(s,
        "fresh"    = '<span class="badge bg-success">Fresh</span>',
        "expiring" = '<span class="badge bg-warning text-dark">Expiring soon</span>',
        "expired"  = '<span class="badge bg-danger">Expired</span>',
        s
      )
    }, character(1))

    # Hide the id column (used for lookup)
    show_df <- display_df[, !names(display_df) %in% "id", drop = FALSE]

    DT::datatable(show_df,
      selection = "single",
      escape = FALSE,  # Allow HTML in Status column
      options = list(
        pageLength = 10,
        dom = "ftip",
        scrollY = "300px",
        columnDefs = list(
          list(width = "180px", targets = 0),  # Name
          list(width = "120px", targets = 1),  # Organism
          list(width = "80px", targets = 2),   # Proteins
          list(width = "80px", targets = 3),   # Age
          list(width = "100px", targets = 4),  # Status
          list(width = "80px", targets = 5)    # Created By
        )
      ),
      rownames = FALSE,
      class = "compact stripe"
    )
  })

  # Detail panel on row selection
  output$fasta_library_detail_panel <- renderUI({
    sel <- input$fasta_library_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      return(div(class = "text-muted text-center py-3",
        icon("hand-pointer"), " Select a database from the table above to see details"
      ))
    }

    catalog <- fasta_library_catalog()
    if (sel > length(catalog)) return(NULL)
    entry <- catalog[[sel]]

    # Check file existence
    files_ok <- fasta_library_verify_files(entry)
    age_status <- fasta_library_check_age(entry)

    # Format file size
    size_mb <- if (!is.null(entry$file_size_bytes) && entry$file_size_bytes > 0) {
      sprintf("%.1f MB", entry$file_size_bytes / 1e6)
    } else "Unknown"

    # Speclib info
    speclib_info <- if (!is.null(entry$speclib_path) && nzchar(entry$speclib_path %||% "")) {
      tags$span(class = "badge bg-info", icon("bolt"), " Predicted speclib available")
    } else {
      tags$span(class = "text-muted", "None")
    }

    # Search settings
    ss <- entry$search_settings %||% list()

    div(class = "card",
      div(class = "card-body", style = "font-size: 0.88em; padding: 12px;",
        div(class = "row",
          div(class = "col-md-6",
            tags$dl(class = "row mb-0",
              tags$dt(class = "col-sm-5", "Organism:"),
              tags$dd(class = "col-sm-7",
                tags$strong(entry$organism %||% ""),
                if (nzchar(entry$organism_common %||% ""))
                  sprintf(" (%s)", entry$organism_common)
              ),
              tags$dt(class = "col-sm-5", "UniProt proteome:"),
              tags$dd(class = "col-sm-7", tags$code(entry$proteome_id %||% "N/A")),
              tags$dt(class = "col-sm-5", "Content type:"),
              tags$dd(class = "col-sm-7", entry$content_type %||% ""),
              tags$dt(class = "col-sm-5", "Protein count:"),
              tags$dd(class = "col-sm-7",
                format(entry$protein_count %||% 0L, big.mark = ","), " sequences"),
              tags$dt(class = "col-sm-5", "File size:"),
              tags$dd(class = "col-sm-7", size_mb),
              tags$dt(class = "col-sm-5", "Contaminants:"),
              tags$dd(class = "col-sm-7",
                if (!is.null(entry$contaminant_library))
                  sprintf("%s (%s proteins)",
                    entry$contaminant_library,
                    format(entry$contaminant_count %||% 0L, big.mark = ","))
                else "None"
              ),
              tags$dt(class = "col-sm-5", "Custom sequences:"),
              tags$dd(class = "col-sm-7",
                if ((entry$custom_sequence_count %||% 0L) > 0)
                  sprintf("%d sequences", entry$custom_sequence_count)
                else "None"
              )
            )
          ),
          div(class = "col-md-6",
            tags$dl(class = "row mb-0",
              tags$dt(class = "col-sm-5", "Enzyme:"),
              tags$dd(class = "col-sm-7", ss$enzyme %||% "N/A"),
              tags$dt(class = "col-sm-5", "Missed cleavages:"),
              tags$dd(class = "col-sm-7", as.character(ss$missed_cleavages %||% "")),
              tags$dt(class = "col-sm-5", "Variable mods:"),
              tags$dd(class = "col-sm-7", ss$var_mods %||% "None"),
              tags$dt(class = "col-sm-5", "Fixed mods:"),
              tags$dd(class = "col-sm-7", ss$fixed_mods %||% "None"),
              tags$dt(class = "col-sm-5", "Peptide length:"),
              tags$dd(class = "col-sm-7",
                sprintf("%d-%d aa",
                  ss$min_pep_len %||% 7L, ss$max_pep_len %||% 30L)),
              tags$dt(class = "col-sm-5", "Precursor m/z:"),
              tags$dd(class = "col-sm-7",
                sprintf("%d-%d",
                  as.integer(ss$min_pr_mz %||% 300),
                  as.integer(ss$max_pr_mz %||% 1200))),
              tags$dt(class = "col-sm-5", "Predicted speclib:"),
              tags$dd(class = "col-sm-7", speclib_info)
            )
          )
        ),
        # FASTA files
        div(style = "margin-top: 8px;",
          tags$strong("FASTA files: "),
          tags$ul(style = "margin-bottom: 4px;",
            lapply(entry$fasta_files %||% character(), function(f) {
              tags$li(tags$code(f))
            })
          )
        ),
        # Metadata footer
        div(class = "d-flex justify-content-between align-items-center",
          style = "margin-top: 8px; padding-top: 8px; border-top: 1px solid #dee2e6;",
          div(
            tags$small(class = "text-muted",
              sprintf("Created %s by %s",
                entry$created_at %||% "Unknown",
                entry$created_by %||% "Unknown")),
            if (nzchar(entry$notes %||% ""))
              div(tags$small(class = "text-muted fst-italic",
                icon("comment"), " ", entry$notes))
          ),
          div(
            # File status indicator
            if (!files_ok) {
              tags$span(class = "badge bg-danger",
                icon("exclamation-triangle"), " Files missing")
            } else if (age_status == "expired") {
              tags$span(class = "badge bg-danger",
                icon("clock"), " Expired")
            } else if (age_status == "expiring") {
              tags$span(class = "badge bg-warning text-dark",
                icon("clock"), " Expiring soon")
            } else {
              tags$span(class = "badge bg-success",
                icon("check-circle"), " Ready")
            },
            # Delete button
            actionButton("fasta_library_delete_btn", "Delete",
              class = "btn-outline-danger btn-sm ms-2", icon = icon("trash"))
          )
        )
      )
    )
  })

  # Refresh catalog button
  observeEvent(input$fasta_library_refresh_btn, {
    fasta_library_catalog(fasta_library_load())
    showNotification("Library catalog refreshed", type = "message", duration = 3)
  })

  # "Use This Database" button
  observeEvent(input$fasta_library_use_btn, {
    sel <- input$fasta_library_table_rows_selected

    if (is.null(sel) || length(sel) == 0) {
      showNotification("Please select a database from the table first.", type = "warning")
      return()
    }

    catalog <- fasta_library_catalog()
    if (sel > length(catalog)) return()
    entry <- catalog[[sel]]

    # Verify files exist
    if (!fasta_library_verify_files(entry)) {
      showNotification(
        "FASTA files for this entry are missing from disk. The entry may need to be removed.",
        type = "error", duration = 8)
      return()
    }

    # Check age/expiration
    age_status <- fasta_library_check_age(entry)
    if (age_status == "expired") {
      showNotification(
        paste("This database is over 6 months old and cannot be used for new searches.",
              "Please download a fresh version from UniProt."),
        type = "error", duration = 10)
      return()
    }

    # Get file paths — use remote paths if HPC/SSH mode
    cfg <- ssh_config()
    use_remote <- !is.null(cfg)
    fasta_paths <- fasta_library_file_paths(entry, use_remote = use_remote)

    # Set the FASTA files
    values$diann_fasta_files <- fasta_paths

    # Store the selected library entry for reference
    values$fasta_info <- list(
      n_sequences = entry$protein_count %||% 0L,
      file_size = entry$file_size_bytes %||% 0L,
      library_entry_id = entry$id,
      library_entry_name = entry$name
    )

    # Check for linked speclib
    if (!is.null(entry$speclib_path) && nzchar(entry$speclib_path %||% "")) {
      # Verify speclib still exists
      speclib_exists <- if (use_remote && !is.null(cfg)) {
        check_result <- ssh_exec(cfg,
          paste("test -f", shQuote(entry$speclib_path), "&& echo EXISTS"))
        any(grepl("EXISTS", check_result$stdout))
      } else {
        file.exists(entry$speclib_path)
      }

      if (speclib_exists && age_status != "expired") {
        values$diann_speclib <- entry$speclib_path
        showNotification(
          sprintf("Database loaded: %s\nPredicted speclib available — Step 1 will be skipped.",
            entry$name),
          type = "message", duration = 8)
      } else {
        showNotification(
          sprintf("Database loaded: %s (%s proteins)",
            entry$name,
            format(entry$protein_count %||% 0L, big.mark = ",")),
          type = "message", duration = 6)
      }
    } else {
      # Show expiring warning if applicable
      notify_msg <- sprintf("Database loaded: %s (%s proteins)",
        entry$name,
        format(entry$protein_count %||% 0L, big.mark = ","))
      if (age_status == "expiring") {
        notify_msg <- paste0(notify_msg,
          "\nNote: This database is nearing expiration. Consider refreshing soon.")
      }
      showNotification(notify_msg, type = "message", duration = 6)
    }

    removeModal()
  })

  # Delete library entry
  observeEvent(input$fasta_library_delete_btn, {
    sel <- input$fasta_library_table_rows_selected
    if (is.null(sel) || length(sel) == 0) return()

    catalog <- fasta_library_catalog()
    if (sel > length(catalog)) return()
    entry <- catalog[[sel]]

    showModal(modalDialog(
      title = "Confirm Delete",
      div(
        tags$p(sprintf("Are you sure you want to remove '%s' from the library?",
          entry$name)),
        checkboxInput("fasta_library_delete_files", "Also delete FASTA files from disk",
          value = FALSE)
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("fasta_library_confirm_delete", "Delete",
          class = "btn-danger", icon = icon("trash"))
      )
    ))
  })

  observeEvent(input$fasta_library_confirm_delete, {
    sel <- input$fasta_library_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      removeModal()
      return()
    }

    catalog <- fasta_library_catalog()
    if (sel > length(catalog)) {
      removeModal()
      return()
    }
    entry <- catalog[[sel]]

    # Remove entry
    success <- fasta_library_remove(entry$id,
      delete_files = isTRUE(input$fasta_library_delete_files))

    if (success) {
      showNotification(sprintf("Removed '%s' from library", entry$name),
        type = "message", duration = 5)
      # Refresh catalog
      fasta_library_catalog(fasta_library_load())
    } else {
      showNotification("Failed to remove entry from library", type = "error")
    }

    # Re-open the library modal
    removeModal()
    # Slight delay to let modal close before reopening
    shinyjs::delay(300, {
      shinyjs::click("open_fasta_library_modal")
    })
  })

  # Summary display in sidebar when library DB is selected
  output$fasta_library_selected_summary <- renderUI({
    req(length(values$diann_fasta_files) > 0)

    # Check if the current selection came from the library
    finfo <- values$fasta_info
    if (!is.null(finfo$library_entry_name)) {
      div(style = "font-size: 0.8em; color: #495057; margin-top: 5px;",
        icon("check-circle", style = "color: #28a745;"), " ",
        tags$strong(finfo$library_entry_name), " \u2014 ",
        format(finfo$n_sequences %||% 0L, big.mark = ","), " sequences"
      )
    }
  })

  # ============================================================================
  #    "Add to Library" after UniProt download
  # ============================================================================

  # Show "Add to Library" button after a successful UniProt download
  output$fasta_add_to_library_btn_ui <- renderUI({
    req(values$fasta_info)
    req(values$fasta_info$success %||% !is.null(values$fasta_info$n_sequences))
    # Only show if this wasn't already from the library
    if (!is.null(values$fasta_info$library_entry_id)) return(NULL)
    # Only show if we have downloaded fasta files
    req(length(values$diann_fasta_files) > 0)

    div(style = "margin-top: 8px;",
      actionButton("add_fasta_to_library", "Add to Library",
        class = "btn-outline-primary btn-sm w-100",
        icon = icon("book-medical")),
      tags$small(class = "text-muted d-block mt-1",
        "Save this database to the shared library for reuse")
    )
  })

  observeEvent(input$add_fasta_to_library, {
    req(values$fasta_info, values$uniprot_results)

    # Get the selected UniProt row
    sel <- input$uniprot_results_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      # Try to reconstruct from fasta_info if no selection (modal was closed)
      showNotification(
        "Could not determine the UniProt source. Please re-download the FASTA first.",
        type = "warning")
      return()
    }

    uniprot_row <- values$uniprot_results[sel, ]
    download_result <- values$fasta_info

    # Get contaminant info if used
    contam_name <- input$contaminant_library %||% "none"
    contam_info <- NULL
    if (contam_name != "none") {
      contam_info <- get_contaminant_fasta(contam_name)
    }

    # Get custom sequences
    custom_seq <- input$custom_fasta_sequences
    if (!is.null(custom_seq) && !nzchar(trimws(custom_seq))) custom_seq <- NULL

    # Collect current search params
    search_params <- list(
      enzyme = input$diann_enzyme %||% "K*,R*",
      missed_cleavages = input$diann_missed_cleavages %||% 1L,
      mod_met_ox = isTRUE(input$mod_met_ox),
      unimod4 = isTRUE(input$diann_unimod4),
      min_pep_len = input$min_pep_len %||% 7L,
      max_pep_len = input$max_pep_len %||% 30L,
      min_pr_mz = input$min_pr_mz %||% 300,
      max_pr_mz = input$max_pr_mz %||% 1200
    )

    # Show a dialog to add notes before saving
    showModal(modalDialog(
      title = tagList(icon("book-medical"), " Add to FASTA Library"),
      div(
        tags$p(sprintf("Adding %s to the shared FASTA library.",
          uniprot_row$common_name %||% uniprot_row$organism)),
        textInput("fasta_library_add_notes", "Notes (optional):",
          placeholder = "e.g., Standard human database for routine DIA searches",
          width = "100%"),
        textInput("fasta_library_add_created_by", "Your name:",
          value = Sys.info()[["user"]], width = "100%")
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("fasta_library_confirm_add", "Add to Library",
          class = "btn-success", icon = icon("plus"))
      )
    ))
  })

  observeEvent(input$fasta_library_confirm_add, {
    req(values$fasta_info, values$uniprot_results)

    sel <- input$uniprot_results_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      removeModal()
      showNotification("UniProt selection lost. Please try again.", type = "warning")
      return()
    }

    uniprot_row <- values$uniprot_results[sel, ]
    download_result <- values$fasta_info

    # Get contaminant info
    contam_name <- input$contaminant_library %||% "none"
    contam_info <- NULL
    if (contam_name != "none") {
      contam_info <- get_contaminant_fasta(contam_name)
    }

    # Custom sequences
    custom_seq <- input$custom_fasta_sequences
    if (!is.null(custom_seq) && !nzchar(trimws(custom_seq))) custom_seq <- NULL

    # Search params
    search_params <- list(
      enzyme = input$diann_enzyme %||% "K*,R*",
      missed_cleavages = input$diann_missed_cleavages %||% 1L,
      mod_met_ox = isTRUE(input$mod_met_ox),
      unimod4 = isTRUE(input$diann_unimod4),
      min_pep_len = input$min_pep_len %||% 7L,
      max_pep_len = input$max_pep_len %||% 30L,
      min_pr_mz = input$min_pr_mz %||% 300,
      max_pr_mz = input$max_pr_mz %||% 1200
    )

    # Build the catalog entry
    entry <- fasta_library_build_entry(
      download_result = download_result,
      uniprot_row = uniprot_row,
      content_type = input$fasta_content_type %||% "one_per_gene",
      contam_info = contam_info,
      contam_name = contam_name,
      custom_sequences = custom_seq,
      search_params = search_params,
      created_by = input$fasta_library_add_created_by %||% Sys.info()[["user"]],
      notes = input$fasta_library_add_notes %||% ""
    )

    # Copy FASTA files to library directory
    lib_path <- fasta_library_path()
    entry_dir <- file.path(lib_path, entry$fasta_dir)

    tryCatch({
      dir.create(entry_dir, recursive = TRUE, showWarnings = FALSE)

      # Copy the main FASTA file
      main_fasta_path <- values$diann_fasta_files[1]
      if (file.exists(main_fasta_path)) {
        file.copy(main_fasta_path,
          file.path(entry_dir, basename(main_fasta_path)),
          overwrite = TRUE)
      }

      # Copy contaminant FASTA if used
      if (!is.null(contam_info) && isTRUE(contam_info$success)) {
        file.copy(contam_info$path,
          file.path(entry_dir, basename(contam_info$path)),
          overwrite = TRUE)
      }

      # Write custom sequences if provided
      if (!is.null(custom_seq) && nzchar(trimws(custom_seq))) {
        writeLines(custom_seq,
          file.path(entry_dir, "custom_proteins.fasta"))
        entry$fasta_files <- c(entry$fasta_files, "custom_proteins.fasta")
      }

      # Write metadata.json for non-R tools
      tryCatch({
        jsonlite::write_json(
          entry[!names(entry) %in% "custom_sequences"],
          file.path(entry_dir, "metadata.json"),
          pretty = TRUE, auto_unbox = TRUE)
      }, error = function(e) NULL)  # Non-critical

      # Add to catalog
      success <- fasta_library_add(entry)

      if (success) {
        removeModal()
        showNotification(
          sprintf("Added '%s' to FASTA library", entry$name),
          type = "message", duration = 8)
        # Update fasta_info to reflect library membership
        values$fasta_info$library_entry_id <- entry$id
        values$fasta_info$library_entry_name <- entry$name
      } else {
        showNotification("Failed to save catalog entry", type = "error")
      }
    }, error = function(e) {
      showNotification(
        sprintf("Failed to copy files to library: %s", e$message),
        type = "error", duration = 10)
    })
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

  # SSH mode: derive output base from raw data directory
  observeEvent(input$ssh_output_base_dir, {
    if (nzchar(input$ssh_output_base_dir %||% "")) {
      output_base(input$ssh_output_base_dir)
    }
  })
  observeEvent(input$ssh_raw_data_dir, {
    if (nzchar(input$ssh_raw_data_dir %||% "")) {
      output_base(input$ssh_raw_data_dir)
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
    # Preview shows output as subfolder of input directory
    input_dir <- tryCatch({
      backend <- input$search_backend %||% "hpc"
      if (backend == "hpc" && nzchar(input$ssh_raw_data_dir %||% "")) {
        input$ssh_raw_data_dir
      } else if (!is.null(values$diann_raw_files) && nrow(values$diann_raw_files) > 0) {
        dirname(values$diann_raw_files$full_path[1])
      } else {
        output_base()
      }
    }, error = function(e) output_base())
    file.path(input_dir, paste0("output_", format(Sys.time(), "%Y%m%d_%H%M")))
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

    # Auto-set output directory as <input_dir>/output_YYYYMMDD_HHMM
    timestamp_suffix <- format(Sys.time(), "%Y%m%d_%H%M")
    input_dir <- tryCatch({
      if (backend == "hpc" && nzchar(input$ssh_raw_data_dir %||% "")) {
        input$ssh_raw_data_dir
      } else if (!is.null(values$diann_raw_files) && nrow(values$diann_raw_files) > 0) {
        dirname(values$diann_raw_files$full_path[1])
      } else {
        output_base()
      }
    }, error = function(e) output_base())
    output_dir <- gsub("//+", "/", file.path(input_dir, paste0(analysis_name, "_", timestamp_suffix)))

    if (backend == "local") {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      cfg <- NULL
    } else if (backend == "docker") {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
      cfg <- NULL
    } else {
      cfg <- ssh_config()
      if (!is.null(cfg)) {
        mkdir_res <- ssh_exec(cfg, sprintf("mkdir -p %s %s/logs", shQuote(output_dir), shQuote(output_dir)))
        if (mkdir_res$status != 0) {
          showNotification(paste("Failed to create remote directory:",
            paste(mkdir_res$stdout, collapse = " ")), type = "error")
          return()
        }
      } else {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
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
    #  Custom FASTA sequences — write temp file, append to fasta_files
    # ====================================================================
    custom_fasta_text <- NULL
    if (!is.null(input$custom_fasta_sequences) &&
        nzchar(trimws(input$custom_fasta_sequences))) {
      custom_seq <- trimws(input$custom_fasta_sequences)
      # Basic validation: must start with >
      if (!grepl("^>", custom_seq)) {
        showNotification("Custom sequences must be in FASTA format (start with '>')",
          type = "warning", duration = 8)
      } else {
        custom_fasta_text <- custom_seq
        custom_fasta_local <- file.path(tempdir(), "custom_proteins.fasta")
        writeLines(custom_seq, custom_fasta_local)

        if (backend == "hpc" && !is.null(cfg)) {
          # SSH mode: upload to output dir
          remote_custom_path <- file.path(output_dir, "custom_proteins.fasta")
          scp_upload(cfg, custom_fasta_local, remote_custom_path)
          fasta_files <- c(fasta_files, remote_custom_path)
        } else {
          # Docker/local: write to output dir
          local_custom_path <- file.path(output_dir, "custom_proteins.fasta")
          file.copy(custom_fasta_local, local_custom_path, overwrite = TRUE)
          fasta_files <- c(fasta_files, local_custom_path)
        }
        showNotification("Added custom protein sequences to search",
          type = "message", duration = 5)
      }
    }

    # ====================================================================
    #  Predicted library cache lookup — reuse if same FASTA + params
    # ====================================================================
    cached_entry <- NULL
    if (is.null(values$diann_speclib) || !nzchar(values$diann_speclib)) {
      fasta_seq_count <- values$fasta_info$n_sequences
      cached_entry <- speclib_cache_lookup(fasta_files, search_params, input$search_mode,
                                           custom_fasta_text, fasta_seq_count)
      if (!is.null(cached_entry)) {
        # Verify the cached speclib file still exists
        speclib_exists <- FALSE
        if (backend == "hpc" && !is.null(cfg)) {
          res <- tryCatch(
            ssh_exec(cfg, paste("test -f", shQuote(cached_entry$speclib_path),
                                "&& echo EXISTS")),
            error = function(e) list(stdout = ""))
          speclib_exists <- any(grepl("EXISTS", res$stdout))
        } else {
          speclib_exists <- file.exists(cached_entry$speclib_path)
        }

        if (speclib_exists) {
          values$diann_speclib <- cached_entry$speclib_path
          showNotification(
            sprintf("Reusing predicted library from '%s' \u2014 skipping Step 1",
                    cached_entry$analysis_name),
            type = "message", duration = 10)
        } else {
          cached_entry <- NULL  # File gone, can't reuse
        }
      }

      # Notify when no cache hit and Step 1 will run
      if (is.null(cached_entry) && (is.null(values$diann_speclib) || !nzchar(values$diann_speclib))) {
        showNotification(
          "No cached library found \u2014 Step 1 (library prediction) will run (~30\u201360 min for human proteome)",
          type = "message", duration = 8)
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

      log_file <- file.path(output_dir, "logs", paste0("diann_", analysis_name, ".log"))

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
          fasta_seq_count = values$fasta_info$n_sequences,
          contaminant_library = contam_lib,
          n_raw_files = nrow(values$diann_raw_files),
          raw_file_type = if (nrow(values$diann_raw_files) > 0)
            tools::file_ext(values$diann_raw_files$filename[1]) else "unknown",
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
          fasta_seq_count = values$fasta_info$n_sequences,
          contaminant_library = contam_lib,
          n_raw_files = nrow(values$diann_raw_files),
          raw_file_type = if (nrow(values$diann_raw_files) > 0)
            tools::file_ext(values$diann_raw_files$filename[1]) else "unknown",
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

    } else if (isTRUE(input$parallel_search)) {
      # --- HPC Parallel (5-step SLURM array) submission ---
      sif_path <- input$diann_sif_path
      sbatch_bin <- values$ssh_sbatch_path %||% "sbatch"
      use_login_shell <- is.null(values$ssh_sbatch_path)

      # Generate all 5 scripts — all steps use the same partition/account
      scripts <- generate_parallel_scripts(
        analysis_name = analysis_name,
        raw_files = values$diann_raw_files$full_path,
        fasta_files = fasta_files,
        speclib_path = values$diann_speclib,
        output_dir = output_dir,
        diann_sif = sif_path,
        normalization = input$diann_normalization,
        search_mode = input$search_mode,
        cpus_per_file = input$parallel_cpus %||% 16,
        mem_per_file = input$parallel_mem_gb %||% 32,
        time_per_file = input$parallel_time_hours %||% 2,
        assembly_cpus = 32,
        assembly_mem = 256,
        assembly_time = input$diann_time_hours,
        partition = input$diann_partition,
        account = input$diann_account,
        search_params = search_params,
        max_simultaneous = input$max_simultaneous %||% 20
      )

      # --- Upload all files + launcher, then submit via one SSH call ---
      # Minimizes SSH connections to avoid HPC MaxStartups throttling.
      # Total: 1 SSH (mkdir) + 1 SCP (all files) + 1 SSH (launcher) = 3 connections.
      has_step1 <- !is.null(scripts$step1_library)
      script_names <- c("step1_libpred.sbatch", "step2_firstpass.sbatch",
                         "step3_assembly.sbatch", "step4_finalpass.sbatch",
                         "step5_report.sbatch")
      script_contents <- list(scripts$step1_library, scripts$step2_firstpass,
                               scripts$step3_assembly, scripts$step4_finalpass,
                               scripts$step5_report)
      step_script_paths <- file.path(output_dir, script_names)

      # Build launcher script that chains sbatch submissions with dependencies
      launcher_lines <- c("#!/bin/bash", "set -e", "")
      if (has_step1) {
        launcher_lines <- c(launcher_lines,
          sprintf('JOB1=$(%s %s 2>&1)', sbatch_bin, step_script_paths[1]),
          'JOB1_ID=$(echo "$JOB1" | grep -oP "[0-9]+$")',
          'echo "STEP1:$JOB1_ID"',
          sprintf('JOB2=$(%s --dependency=afterok:$JOB1_ID %s 2>&1)',
                  sbatch_bin, step_script_paths[2])
        )
      } else {
        launcher_lines <- c(launcher_lines,
          sprintf('JOB2=$(%s %s 2>&1)', sbatch_bin, step_script_paths[2]),
          'echo "STEP1:skipped"'
        )
      }
      launcher_lines <- c(launcher_lines,
        'JOB2_ID=$(echo "$JOB2" | grep -oP "[0-9]+$")',
        'echo "STEP2:$JOB2_ID"', "",
        sprintf('JOB3=$(%s --dependency=afterok:$JOB2_ID %s 2>&1)',
                sbatch_bin, step_script_paths[3]),
        'JOB3_ID=$(echo "$JOB3" | grep -oP "[0-9]+$")',
        'echo "STEP3:$JOB3_ID"', "",
        sprintf('JOB4=$(%s --dependency=afterok:$JOB3_ID %s 2>&1)',
                sbatch_bin, step_script_paths[4]),
        'JOB4_ID=$(echo "$JOB4" | grep -oP "[0-9]+$")',
        'echo "STEP4:$JOB4_ID"', "",
        sprintf('JOB5=$(%s --dependency=afterok:$JOB4_ID %s 2>&1)',
                sbatch_bin, step_script_paths[5]),
        'JOB5_ID=$(echo "$JOB5" | grep -oP "[0-9]+$")',
        'echo "STEP5:$JOB5_ID"'
      )

      # Write file_list.txt locally
      file_list_local <- write_file_list(values$diann_raw_files$full_path, tempdir())

      # Write everything to a temp dir for upload
      upload_dir <- tempfile("delimp_scripts_")
      dir.create(upload_dir)
      on.exit(unlink(upload_dir, recursive = TRUE), add = TRUE)
      for (i in seq_along(script_names)) {
        if (!is.null(script_contents[[i]])) {
          writeLines(script_contents[[i]], file.path(upload_dir, script_names[i]))
        }
      }
      file.copy(file_list_local, file.path(upload_dir, "file_list.txt"))
      writeLines(paste(launcher_lines, collapse = "\n"),
                 file.path(upload_dir, "submit_all.sh"))

      # Write search_info.md — archives all metadata so recovery works
      # even after SLURM purges job records
      search_info <- generate_search_info(
        analysis_name = analysis_name,
        output_dir = output_dir,
        raw_files = values$diann_raw_files$full_path,
        fasta_files = fasta_files,
        search_params = search_params,
        search_mode = input$search_mode,
        normalization = input$diann_normalization,
        sif_path = sif_path,
        job_ids = NULL,  # Updated after submission with actual IDs
        parallel = TRUE,
        resources = list(
          "Step 1 (Library Prediction)" = list(cpus = 16, mem = 64, time = 4),
          "Steps 2/4 (Per-file Quant)" = list(
            cpus = input$parallel_cpus %||% 16,
            mem = input$parallel_mem_gb %||% 32,
            time = input$parallel_time_hours %||% 2),
          "Steps 3/5 (Assembly/Report)" = list(
            cpus = 32, mem = 256,
            time = input$diann_time_hours)
        ),
        partition = input$diann_partition,
        account = input$diann_account,
        cached_speclib = cached_entry,
        custom_fasta_sequences = custom_fasta_text
      )
      writeLines(search_info, file.path(upload_dir, "search_info.md"))

      step_ids <- new.env(parent = emptyenv())
      pstate <- new.env(parent = emptyenv())
      pstate$failed <- FALSE

      if (!is.null(cfg)) {
        # SSH mode: 1 mkdir + 1 SCP + 1 bash = 3 SSH connections total
        mkdir_cmd <- sprintf("mkdir -p %s %s/logs %s/quant_step2 %s/quant_step4",
                              shQuote(output_dir), shQuote(output_dir),
                              shQuote(output_dir), shQuote(output_dir))
        mkdir_result <- ssh_exec(cfg, mkdir_cmd)
        if (mkdir_result$status != 0) {
          showNotification(paste("Failed to create remote directories:",
            paste(mkdir_result$stdout, collapse = " ")), type = "error")
          return()
        }

        # Single SCP: upload all sbatch scripts + file_list + launcher
        local_files <- list.files(upload_dir, full.names = TRUE)
        scp_args <- c(
          "-i", cfg$key_path,
          "-P", as.character(cfg$port %||% 22),
          "-o", "StrictHostKeyChecking=accept-new",
          "-o", "ConnectTimeout=10",
          "-o", "BatchMode=yes",
          ssh_mux_args(cfg),
          local_files,
          paste0(cfg$user, "@", cfg$host, ":", output_dir, "/")
        )
        message("[DE-LIMP] Uploading ", length(local_files), " files to ", output_dir)
        scp_result <- tryCatch({
          if (requireNamespace("processx", quietly = TRUE)) {
            res <- processx::run("scp", args = scp_args, timeout = 120,
                                 error_on_status = FALSE,
                                 env = c("current", MallocStackLogging = ""))
            list(status = res$status, stdout = iconv(paste0(res$stdout, res$stderr),
                                                      to = "UTF-8", sub = ""))
          } else {
            out <- system2("scp", args = scp_args, stdout = TRUE, stderr = TRUE)
            list(status = attr(out, "status") %||% 0L,
                 stdout = iconv(out, to = "UTF-8", sub = ""))
          }
        }, error = function(e) list(status = 1L, stdout = e$message))

        if (scp_result$status != 0) {
          showNotification(paste("Failed to upload scripts:", scp_result$stdout),
            type = "error")
          return()
        }
        message("[DE-LIMP] All files uploaded. Executing launcher...")

        # Execute launcher: one SSH call submits all 5 sbatch jobs
        launcher_remote <- file.path(output_dir, "submit_all.sh")
        result <- ssh_exec(cfg, paste("bash", launcher_remote),
                           login_shell = use_login_shell, timeout = 120)
        message("[DE-LIMP] Launcher status=", result$status,
                " stdout=", paste(result$stdout, collapse = " | "))

        if (result$status != 0) {
          showNotification(paste("Parallel submission failed:",
            paste(result$stdout, collapse = " ")), type = "error", duration = 15)
          return()
        }

        # Parse STEP1:id through STEP5:id from output
        for (line in result$stdout) {
          m <- regmatches(line, regexec("^STEP([1-5]):(.+)$", line))[[1]]
          if (length(m) == 3) {
            step_name <- paste0("step", m[2])
            job_id_val <- trimws(m[3])
            if (job_id_val != "skipped" && nzchar(job_id_val)) {
              step_ids[[step_name]] <- job_id_val
            }
            message("[DE-LIMP] ", step_name, " = ", job_id_val)
          }
        }

        if (is.null(step_ids$step5)) {
          showNotification("Could not parse all step job IDs from sbatch output.",
            type = "error", duration = 15)
          return()
        }

      } else {
        # Local mode: copy scripts to output_dir (SSH mode uploads via SCP)
        dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
        dir.create(file.path(output_dir, "quant_step2"), showWarnings = FALSE)
        dir.create(file.path(output_dir, "quant_step4"), showWarnings = FALSE)
        local_files <- list.files(upload_dir, full.names = TRUE)
        file.copy(local_files, output_dir, overwrite = TRUE)

        # Local mode: submit sequentially
        local_sbatch_bin <- if (nzchar(local_sbatch_path)) local_sbatch_path else "sbatch"
        has_step1 <- !is.null(scripts$step1_library)

        withProgress(message = "Submitting 5-step parallel search...", value = 0, {
          # Step 1
          if (!pstate$failed && has_step1) {
            incProgress(0.1, detail = "Step 1: Library prediction")
            tryCatch({
              out <- system2(local_sbatch_bin,
                args = file.path(output_dir, "step1_libpred.sbatch"),
                stdout = TRUE, stderr = TRUE)
              step_ids$step1 <- parse_sbatch_output(out)
              if (is.null(step_ids$step1)) pstate$failed <- TRUE
            }, error = function(e) { pstate$failed <- TRUE })
          }
          # Step 2
          if (!pstate$failed) {
            incProgress(0.2, detail = "Step 2: First-pass array")
            dep <- if (!is.null(step_ids$step1))
              sprintf("--dependency=afterok:%s", step_ids$step1)
            tryCatch({
              out <- system2(local_sbatch_bin,
                args = c(dep, file.path(output_dir, "step2_firstpass.sbatch")),
                stdout = TRUE, stderr = TRUE)
              step_ids$step2 <- parse_sbatch_output(out)
              if (is.null(step_ids$step2)) pstate$failed <- TRUE
            }, error = function(e) { pstate$failed <- TRUE })
          }
          # Step 3
          if (!pstate$failed) {
            incProgress(0.2, detail = "Step 3: Library assembly")
            tryCatch({
              out <- system2(local_sbatch_bin,
                args = c(sprintf("--dependency=afterok:%s", step_ids$step2),
                         file.path(output_dir, "step3_assembly.sbatch")),
                stdout = TRUE, stderr = TRUE)
              step_ids$step3 <- parse_sbatch_output(out)
              if (is.null(step_ids$step3)) pstate$failed <- TRUE
            }, error = function(e) { pstate$failed <- TRUE })
          }
          # Step 4
          if (!pstate$failed) {
            incProgress(0.2, detail = "Step 4: Final-pass array")
            tryCatch({
              out <- system2(local_sbatch_bin,
                args = c(sprintf("--dependency=afterok:%s", step_ids$step3),
                         file.path(output_dir, "step4_finalpass.sbatch")),
                stdout = TRUE, stderr = TRUE)
              step_ids$step4 <- parse_sbatch_output(out)
              if (is.null(step_ids$step4)) pstate$failed <- TRUE
            }, error = function(e) { pstate$failed <- TRUE })
          }
          # Step 5
          if (!pstate$failed) {
            incProgress(0.2, detail = "Step 5: Cross-run report")
            tryCatch({
              out <- system2(local_sbatch_bin,
                args = c(sprintf("--dependency=afterok:%s", step_ids$step4),
                         file.path(output_dir, "step5_report.sbatch")),
                stdout = TRUE, stderr = TRUE)
              step_ids$step5 <- parse_sbatch_output(out)
              if (is.null(step_ids$step5)) pstate$failed <- TRUE
            }, error = function(e) { pstate$failed <- TRUE })
          }
        })

        if (pstate$failed) {
          showNotification("Parallel submission failed (local mode).", type = "error")
          return()
        }
      }

      # Abort if any step failed to submit
      if (pstate$failed) return()

      # The main job_id is the final step (its completion = workflow done)
      job_id <- step_ids$step5

      # Create parallel HPC job entry
      job_entry <- list(
        job_id = job_id,
        backend = "hpc",
        name = analysis_name,
        status = "queued",
        output_dir = output_dir,
        script_path = file.path(output_dir, "step5_report.sbatch"),
        submitted_at = Sys.time(),
        n_files = nrow(values$diann_raw_files),
        search_mode = input$search_mode,
        parallel = TRUE,
        parallel_steps = as.list(step_ids),
        parallel_n_files = nrow(values$diann_raw_files),
        parallel_current_step = if (is.null(step_ids$step1)) 2L else 1L,
        parallel_step_status = list(
          step1 = if (is.null(step_ids$step1)) "skipped" else "queued",
          step2 = "queued", step3 = "queued",
          step4 = "queued", step5 = "queued"
        ),
        search_settings = list(
          search_params = search_params,
          fasta_files = fasta_files,
          fasta_seq_count = values$fasta_info$n_sequences,
          contaminant_library = contam_lib,
          custom_fasta_text = custom_fasta_text,
          n_raw_files = nrow(values$diann_raw_files),
          raw_file_type = if (nrow(values$diann_raw_files) > 0)
            tools::file_ext(values$diann_raw_files$filename[1]) else "unknown",
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
          ),
          parallel = list(
            cpus_per_file = input$parallel_cpus %||% 16,
            mem_per_file = input$parallel_mem_gb %||% 32,
            time_per_file = input$parallel_time_hours %||% 2,
            max_simultaneous = input$max_simultaneous %||% 20
          )
        ),
        auto_load = input$auto_load_results,
        log_content = "",
        completed_at = NULL,
        loaded = FALSE,
        is_ssh = !is.null(cfg),
        speclib_cached = !is.null(cached_entry)
      )

      # Update search_info.md with actual job IDs
      tryCatch({
        job_id_list <- as.list(step_ids)
        updated_info <- generate_search_info(
          analysis_name = analysis_name,
          output_dir = output_dir,
          raw_files = values$diann_raw_files$full_path,
          fasta_files = fasta_files,
          search_params = search_params,
          search_mode = input$search_mode,
          normalization = input$diann_normalization,
          sif_path = sif_path,
          job_ids = job_id_list,
          parallel = TRUE,
          resources = list(
            "Step 1 (Library Prediction)" = list(cpus = 16, mem = 64, time = 4),
            "Steps 2/4 (Per-file Quant)" = list(
              cpus = input$parallel_cpus %||% 16,
              mem = input$parallel_mem_gb %||% 32,
              time = input$parallel_time_hours %||% 2),
            "Steps 3/5 (Assembly/Report)" = list(
              cpus = 32, mem = 256,
              time = input$diann_time_hours)
          ),
          partition = input$diann_partition,
          account = input$diann_account,
          cached_speclib = cached_entry,
          custom_fasta_sequences = custom_fasta_text
        )
        local_info <- tempfile(fileext = ".md")
        writeLines(updated_info, local_info)
        if (!is.null(cfg)) {
          scp_upload(cfg, local_info, file.path(output_dir, "search_info.md"))
        } else {
          file.copy(local_info, file.path(output_dir, "search_info.md"), overwrite = TRUE)
        }
        unlink(local_info)
      }, error = function(e) message("[DE-LIMP] Could not update search_info.md: ", e$message))

    } else {
      # --- HPC (SLURM) standard single-job submission ---
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
        search_params = search_params,
        requeue = (tolower(input$diann_partition) == "low")
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

        # Use full sbatch path (may be on CVMFS inside Apptainer container)
        sbatch_local <- if (nzchar(local_sbatch_path)) local_sbatch_path else "sbatch"
        submit_result <- tryCatch({
          stdout <- system2(sbatch_local, args = script_path, stdout = TRUE, stderr = TRUE)
          exit_code <- attr(stdout, "status")
          list(success = is.null(exit_code) || exit_code == 0L, stdout = stdout)
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
          fasta_seq_count = values$fasta_info$n_sequences,
          contaminant_library = contam_lib,
          n_raw_files = nrow(values$diann_raw_files),
          raw_file_type = if (nrow(values$diann_raw_files) > 0)
            tools::file_ext(values$diann_raw_files$filename[1]) else "unknown",
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

      # Write search_info.md for single-job submission
      tryCatch({
        search_info <- generate_search_info(
          analysis_name = analysis_name,
          output_dir = output_dir,
          raw_files = values$diann_raw_files$full_path,
          fasta_files = fasta_files,
          search_params = search_params,
          search_mode = input$search_mode,
          normalization = input$diann_normalization,
          sif_path = sif_path,
          job_ids = job_id,
          parallel = FALSE,
          resources = list(
            "Single Job" = list(
              cpus = input$diann_cpus,
              mem = input$diann_mem_gb,
              time = input$diann_time_hours)
          ),
          partition = input$diann_partition,
          account = input$diann_account,
          cached_speclib = cached_entry,
          custom_fasta_sequences = custom_fasta_text
        )
        local_info <- tempfile(fileext = ".md")
        writeLines(search_info, local_info)
        if (!is.null(cfg)) {
          scp_upload(cfg, local_info, file.path(output_dir, "search_info.md"))
        } else {
          file.copy(local_info, file.path(output_dir, "search_info.md"), overwrite = TRUE)
        }
        unlink(local_info)
      }, error = function(e) message("[DE-LIMP] Could not write search_info.md: ", e$message))
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
          lc_method     = input$search_lc_method %||% "",
          project       = input$search_project %||% "",
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
      if (isTRUE(jobs[[i]]$removed)) next
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
      } else if (isTRUE(jobs[[i]]$parallel)) {
        # --- HPC Parallel (5-step) monitoring ---
        job_cfg <- if (isTRUE(jobs[[i]]$is_ssh)) cfg else NULL
        slurm_path <- if (isTRUE(jobs[[i]]$is_ssh)) {
          values$ssh_sbatch_path
        } else if (nzchar(local_sbatch_path)) {
          local_sbatch_path
        } else NULL

        steps <- jobs[[i]]$parallel_steps
        step_status <- jobs[[i]]$parallel_step_status %||% list()
        step_names <- c("step1", "step2", "step3", "step4", "step5")
        current_step <- jobs[[i]]$parallel_current_step %||% 1L
        n_files <- jobs[[i]]$parallel_n_files %||% 0

        # Poll each non-terminal step
        for (sn in step_names) {
          sid <- steps[[sn]]
          if (is.null(sid)) next
          prev <- step_status[[sn]] %||% "queued"
          if (prev %in% c("completed", "skipped", "failed", "cancelled")) next

          s_status <- check_slurm_status(sid, ssh_config = job_cfg, sbatch_path = slurm_path)
          step_status[[sn]] <- s_status
        }
        jobs[[i]]$parallel_step_status <- step_status

        # Determine current_step (first non-completed step)
        for (si in seq_along(step_names)) {
          ss <- step_status[[step_names[si]]] %||% "queued"
          if (!ss %in% c("completed", "skipped")) {
            current_step <- si
            break
          }
        }
        jobs[[i]]$parallel_current_step <- current_step

        # For array steps (2 & 4), count completed tasks via sacct
        for (array_step in c("step2", "step4")) {
          arr_id <- steps[[array_step]]
          arr_status <- step_status[[array_step]] %||% "queued"
          if (is.null(arr_id) || !arr_status %in% c("running", "queued")) next

          slurm_cmd_fn <- function(cmd) {
            if (!is.null(slurm_path)) file.path(dirname(slurm_path), cmd)
            else cmd
          }
          sacct_cmd <- sprintf(
            "%s -j %s --format=State --noheader --parsable2 2>/dev/null",
            slurm_cmd_fn("sacct"), arr_id)

          sacct_result <- if (!is.null(job_cfg)) {
            ssh_exec(job_cfg, sacct_cmd, login_shell = is.null(slurm_path))
          } else {
            tryCatch({
              out <- system2(slurm_cmd_fn("sacct"),
                args = c("-j", arr_id, "--format=State", "--noheader", "--parsable2"),
                stdout = TRUE, stderr = TRUE)
              list(status = 0, stdout = out)
            }, error = function(e) list(status = 1, stdout = character()))
          }

          if (sacct_result$status == 0 && length(sacct_result$stdout) > 0) {
            states <- toupper(trimws(sacct_result$stdout))
            states <- states[nzchar(states)]
            n_done <- sum(grepl("COMPLETED", states))
            n_running <- sum(grepl("RUNNING", states))
            n_pending <- sum(grepl("PENDING", states))
            n_failed <- sum(grepl("FAILED|TIMEOUT|OUT_OF_ME", states))
            jobs[[i]][[paste0(array_step, "_progress")]] <- list(
              completed = n_done, running = n_running,
              pending = n_pending, failed = n_failed)
          }
        }

        # Overall status: check step5 final state
        step5_status <- step_status[["step5"]] %||% "queued"
        new_status <- if (step5_status == "completed") "completed"
          else if (any(vapply(step_status, function(s) s %in% c("failed"), logical(1)))) "failed"
          else if (any(vapply(step_status, function(s) s %in% c("cancelled"), logical(1)))) "cancelled"
          else if (any(vapply(step_status, function(s) s %in% c("running"), logical(1)))) "running"
          else "queued"

        # Tail the most recent log for display
        if (isTRUE(jobs[[i]]$is_ssh) && !is.null(cfg)) {
          log_result <- ssh_exec(cfg, sprintf(
            "{ ls -t %1$s/logs/diann_*.out %1$s/logs/diann_*.err %1$s/diann_*.out %1$s/diann_*.err 2>/dev/null; } | head -1 | xargs tail -50 2>/dev/null",
            shQuote(jobs[[i]]$output_dir)))
          if (log_result$status == 0 && length(log_result$stdout) > 0) {
            jobs[[i]]$log_content <- paste(log_result$stdout, collapse = "\n")
          }
        } else {
          log_dirs <- c(file.path(jobs[[i]]$output_dir, "logs"), jobs[[i]]$output_dir)
          log_files <- unlist(lapply(log_dirs, list.files,
            pattern = "^diann_.*\\.(out|err)$", full.names = TRUE))
          if (length(log_files) > 0) {
            log_file <- log_files[which.max(file.mtime(log_files))]
            log_lines <- tryCatch(
              tail(readLines(log_file, warn = FALSE), 50),
              error = function(e) character(0))
            jobs[[i]]$log_content <- paste(log_lines, collapse = "\n")
          }
        }
        changed <- TRUE

      } else {
        # --- HPC (SLURM) standard single-job monitoring ---
        job_cfg <- if (isTRUE(jobs[[i]]$is_ssh)) cfg else NULL
        # Use SSH sbatch path for remote jobs, local path for local jobs
        slurm_path <- if (isTRUE(jobs[[i]]$is_ssh)) {
          values$ssh_sbatch_path
        } else if (nzchar(local_sbatch_path)) {
          local_sbatch_path
        } else {
          NULL
        }
        new_status <- check_slurm_status(jobs[[i]]$job_id, ssh_config = job_cfg,
                                          sbatch_path = slurm_path)

        # Tail the log file (local or remote)
        if (isTRUE(jobs[[i]]$is_ssh) && !is.null(cfg)) {
          log_result <- ssh_exec(cfg, sprintf(
            "{ ls -t %1$s/logs/diann_*.out %1$s/logs/diann_*.err %1$s/diann_*.out %1$s/diann_*.err 2>/dev/null; } | head -1 | xargs tail -50 2>/dev/null",
            shQuote(jobs[[i]]$output_dir)))
          if (log_result$status == 0 && length(log_result$stdout) > 0) {
            jobs[[i]]$log_content <- paste(log_result$stdout, collapse = "\n")
            changed <- TRUE
          }
        } else {
          log_dirs <- c(file.path(jobs[[i]]$output_dir, "logs"), jobs[[i]]$output_dir)
          log_files <- unlist(lapply(log_dirs, list.files,
            pattern = "^diann_.*\\.(out|err)$", full.names = TRUE))
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

        # Sync status to Core Facility SQLite database
        if (is_core_facility && !is.null(cf_config)) {
          tryCatch({
            job_id_key <- jobs[[i]]$container_id %||% jobs[[i]]$job_id
            cf_update_search_status(cf_config$db_path, job_id_key, new_status)
          }, error = function(e) message("CF DB update failed: ", e$message))
        }

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
          # QC auto-ingest: download report.parquet and ingest metrics
          if (isTRUE(jobs[[i]]$qc_run) && !isTRUE(jobs[[i]]$qc_ingested) &&
              is_core_facility && !is.null(cf_config)) {
            tryCatch({
              remote_report <- file.path(jobs[[i]]$output_dir, "report.parquet")
              local_report <- file.path(tempdir(),
                paste0(jobs[[i]]$name, "_report.parquet"))

              dl_ok <- FALSE
              if (isTRUE(jobs[[i]]$is_ssh) && !is.null(cfg)) {
                dl_result <- scp_download(cfg, remote_report, local_report)
                dl_ok <- dl_result$status == 0 && file.exists(local_report)
              } else if (file.exists(remote_report)) {
                file.copy(remote_report, local_report)
                dl_ok <- TRUE
              }

              if (dl_ok) {
                cf_ingest_qc_metrics(
                  db_path = cf_config$db_path,
                  report_path = local_report,
                  instrument = jobs[[i]]$qc_instrument %||% "Unknown",
                  run_name = jobs[[i]]$name,
                  search_id = NULL,
                  ng_loaded = jobs[[i]]$qc_ng_loaded,
                  gradient = jobs[[i]]$qc_gradient
                )
                jobs[[i]]$qc_ingested <- TRUE
                changed <- TRUE
                showNotification(
                  sprintf("QC metrics auto-ingested for '%s'!", jobs[[i]]$name),
                  type = "message", duration = 10)
              } else {
                message("[DE-LIMP] QC auto-ingest: report.parquet not found for ",
                  jobs[[i]]$name)
              }
            }, error = function(e) {
              message("[DE-LIMP] QC auto-ingest error: ", e$message)
              showNotification(
                sprintf("QC auto-ingest failed for '%s': %s",
                  jobs[[i]]$name, e$message),
                type = "warning", duration = 10)
            })
          }

          # Register predicted spectral library in cache
          if (isTRUE(jobs[[i]]$parallel) && !isTRUE(jobs[[i]]$speclib_cached)) {
            tryCatch({
              ss <- jobs[[i]]$search_settings
              speclib_path <- file.path(jobs[[i]]$output_dir, "step1.predicted.speclib")
              speclib_exists <- if (isTRUE(jobs[[i]]$is_ssh) && !is.null(cfg)) {
                res <- ssh_exec(cfg, paste("test -f", shQuote(speclib_path),
                                            "&& echo EXISTS"))
                any(grepl("EXISTS", res$stdout))
              } else {
                file.exists(speclib_path)
              }
              if (speclib_exists) {
                speclib_cache_register(
                  fasta_files = ss$fasta_files,
                  search_params = ss$search_params,
                  search_mode = ss$search_mode,
                  speclib_path = speclib_path,
                  analysis_name = jobs[[i]]$name,
                  output_dir = jobs[[i]]$output_dir,
                  custom_fasta_text = ss$custom_fasta_text,
                  fasta_seq_count = ss$fasta_seq_count
                )
                jobs[[i]]$speclib_cached <- TRUE
                changed <- TRUE
              }
            }, error = function(e) {
              message("[DE-LIMP] speclib cache registration failed: ", e$message)
            })
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

          # Update Core Facility DB with protein/peptide counts
          if (is_core_facility && !is.null(cf_config)) {
            tryCatch({
              job_id_key <- job$container_id %||% job$job_id
              n_prot <- if (!is.null(values$raw_data) && !is.null(values$raw_data$genes$Protein.Group))
                length(unique(values$raw_data$genes$Protein.Group))
              else if (!is.null(values$raw_data)) nrow(values$raw_data$E) else NA
              n_pep <- if (!is.null(values$raw_data) && !is.null(values$raw_data$genes$Stripped.Sequence)) {
                length(unique(values$raw_data$genes$Stripped.Sequence))
              } else NA
              cf_update_search_status(cf_config$db_path, job_id_key, "completed",
                                      n_proteins = n_prot, n_peptides = n_pep)
            }, error = function(e) message("CF DB stats update failed: ", e$message))
          }

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

          # Record to analysis history log
          tryCatch({
            record_analysis(list(
              timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              user = Sys.getenv("USER", "unknown"),
              source_type = "auto-load",
              source_file = basename(report_path),
              source_path = report_path,
              remote_path = if (isTRUE(job$is_ssh)) remote_report else NA,
              fasta_file = if (!is.null(job$search_settings))
                paste(basename(job$search_settings$fasta_files), collapse = ", ") else NA,
              fasta_seq_count = if (!is.null(job$search_settings)) job$search_settings$fasta_seq_count else NA,
              n_proteins = if (!is.null(values$raw_data) && !is.null(values$raw_data$genes$Protein.Group))
                length(unique(values$raw_data$genes$Protein.Group))
              else if (!is.null(values$raw_data)) nrow(values$raw_data$E) else NA,
              n_samples = if (!is.null(values$raw_data)) ncol(values$raw_data$E) else NA,
              n_contrasts = NA,
              n_de_proteins = NA,
              output_dir = job$output_dir,
              app_version = values$app_version %||% "unknown",
              notes = sprintf("Job: %s (%s)", job$name, job$job_id)
            ))
          }, error = function(e) message("[DE-LIMP] History record failed: ", e$message))

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
    active_jobs <- Filter(function(j) !isTRUE(j$removed), jobs)
    if (length(active_jobs) == 0) {
      return(div(style = "color: #999; font-size: 0.85em; text-align: center; padding: 10px;",
        "No jobs submitted yet."
      ))
    }

    # Refresh all button at top
    has_unknown <- any(vapply(jobs, function(j) !isTRUE(j$removed) && j$status == "unknown", logical(1)))

    job_rows <- lapply(seq_along(jobs), function(i) {
      job <- jobs[[i]]
      if (isTRUE(job$removed)) return(NULL)

      status_badge <- switch(job$status,
        "queued"    = span(class = "badge bg-secondary", "Queued"),
        "running"   = span(class = "badge bg-primary", "Running"),
        "completed" = span(class = "badge bg-success", "Completed"),
        "failed"    = span(class = "badge bg-danger", title = job$failure_reason %||% "", "Failed"),
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
      } else if (isTRUE(job$parallel)) {
        span(class = "badge bg-primary text-white", style = "font-size: 0.7em; margin-right: 4px;",
          icon("layer-group"), " HPC Parallel")
      } else {
        span(class = "badge bg-secondary", style = "font-size: 0.7em; margin-right: 4px;",
          icon("server"), " HPC")
      }

      # Build parallel step progress display
      parallel_progress_ui <- if (isTRUE(job$parallel)) {
        step_labels <- c("1: Lib Predict", "2: First Pass", "3: Assembly",
                          "4: Final Pass", "5: Report")
        step_names <- c("step1", "step2", "step3", "step4", "step5")
        step_status <- job$parallel_step_status %||% list()
        n_files <- job$parallel_n_files %||% 0

        step_icons <- lapply(seq_along(step_names), function(si) {
          sn <- step_names[si]
          ss <- step_status[[sn]] %||% "queued"
          step_icon <- switch(ss,
            "completed" = icon("check", style = "color: #28a745;"),
            "running"   = icon("spinner", class = "fa-spin", style = "color: #007bff;"),
            "failed"    = icon("xmark", style = "color: #dc3545;"),
            "skipped"   = icon("forward", style = "color: #999;"),
            icon("clock", style = "color: #999;"))

          # Array step progress for steps 2 & 4
          progress_text <- ""
          if (sn %in% c("step2", "step4") && ss == "running") {
            prog <- job[[paste0(sn, "_progress")]]
            if (!is.null(prog)) {
              pct <- if (n_files > 0) round(prog$completed / n_files * 100) else 0
              progress_text <- sprintf(" (%d/%d, %d%%)", prog$completed, n_files, pct)
            }
          }

          div(style = "display: flex; align-items: center; gap: 4px; font-size: 0.78em; padding: 1px 0;",
            step_icon, span(paste0(step_labels[si], progress_text)))
        })

        div(style = "margin-top: 4px; border-top: 1px dashed #dee2e6; padding-top: 4px;",
          step_icons)
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
            actionButton(sprintf("view_info_%d", i), "Info",
              class = "btn-outline-info btn-xs",
              style = "font-size: 0.75em; padding: 2px 6px;",
              icon = icon("circle-info")),
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
            } else if (job$status == "completed" && isTRUE(job$loaded)) {
              actionButton(sprintf("load_results_%d", i), "Reload",
                class = "btn-outline-secondary btn-xs",
                style = "font-size: 0.75em; padding: 2px 6px;",
                icon = icon("rotate-right"))
            },
            if (job$status %in% c("failed", "cancelled") &&
                isTRUE(job$backend == "hpc")) {
              actionButton(sprintf("resubmit_job_%d", i), "Resubmit",
                class = "btn-outline-warning btn-xs",
                style = "font-size: 0.75em; padding: 2px 6px;",
                icon = icon("rotate-right"))
            },
            if (job$status %in% c("completed", "failed", "cancelled")) {
              actionButton(sprintf("remove_job_%d", i), NULL,
                class = "btn-outline-secondary btn-xs",
                style = "font-size: 0.75em; padding: 2px 6px;",
                icon = icon("xmark"))
            }
          )
        ),
        parallel_progress_ui
      )
    })

    # Count terminal and failed jobs for action buttons
    n_terminal <- sum(vapply(jobs, function(j)
      !isTRUE(j$removed) && j$status %in% c("completed", "failed", "cancelled"), logical(1)))
    n_failed <- sum(vapply(jobs, function(j)
      !isTRUE(j$removed) && j$status %in% c("failed", "cancelled"), logical(1)))

    tagList(
      div(style = "display: flex; justify-content: flex-end; gap: 6px; margin-bottom: 6px;",
        if (n_failed >= 1) actionButton("clear_failed_jobs", "Clear Failed",
          class = "btn-outline-danger btn-xs",
          style = "font-size: 0.75em; padding: 2px 8px;",
          icon = icon("trash-can")),
        if (n_terminal >= 2) actionButton("clear_finished_jobs", "Clear Finished",
          class = "btn-outline-secondary btn-xs",
          style = "font-size: 0.75em; padding: 2px 8px;",
          icon = icon("broom")),
        if (has_unknown) actionButton("refresh_all_jobs", "Refresh All",
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

          # On-demand log fetch: if cached log says "Could not locate" but we have
          # an output_dir, try fetching the log from the cluster now
          if (grepl("Could not locate log file", job$log_content %||% "") &&
              nzchar(job$output_dir %||% "") && job$output_dir != "(unknown)" &&
              isTRUE(job$is_ssh)) {
            cfg <- isolate(ssh_config())
            if (!is.null(cfg)) {
              log_path <- file.path(job$output_dir, "logs", sprintf("diann_%s.out", job$job_id))
              # Fallback to old location (pre-logs-subdir) if not found
              tail_cmd <- sprintf(
                "if [ -f %s ]; then tail -150 %s; else tail -150 %s 2>/dev/null; fi",
                shQuote(log_path), shQuote(log_path),
                shQuote(file.path(job$output_dir, sprintf("diann_%s.out", job$job_id))))
              tail_result <- tryCatch(
                ssh_exec(cfg, tail_cmd, timeout = 15),
                error = function(e) list(status = 1, stdout = character()))
              if (tail_result$status == 0 && length(tail_result$stdout) > 0) {
                fetched <- iconv(paste(tail_result$stdout, collapse = "\n"),
                  from = "", to = "UTF-8", sub = "")
                if (nzchar(fetched)) {
                  jobs <- values$diann_jobs
                  jobs[[idx]]$log_content <- fetched
                  values$diann_jobs <- jobs
                  job <- jobs[[idx]]
                }
              }
            }
          }

          safe_log <- iconv(job$log_content %||% "", from = "", to = "UTF-8", sub = "")
          showModal(modalDialog(
            title = sprintf("Log: %s (#%s)", job$name, job$job_id),
            size = "l", easyClose = TRUE, footer = modalButton("Close"),
            pre(style = "max-height: 500px; overflow-y: auto; font-size: 0.8em;",
              safe_log
            )
          ))
        }, ignoreInit = TRUE)

        # View search_info.md
        observeEvent(input[[sprintf("view_info_%d", idx)]], {
          job <- values$diann_jobs[[idx]]
          info_content <- ""

          if (nzchar(job$output_dir %||% "") && job$output_dir != "(unknown)") {
            info_path <- file.path(job$output_dir, "search_info.md")
            if (isTRUE(job$is_ssh)) {
              cfg <- isolate(ssh_config())
              if (!is.null(cfg)) {
                result <- tryCatch(
                  ssh_exec(cfg, sprintf("cat %s 2>/dev/null", shQuote(info_path)), timeout = 15),
                  error = function(e) list(status = 1, stdout = character()))
                if (result$status == 0 && length(result$stdout) > 0) {
                  info_content <- paste(result$stdout, collapse = "\n")
                }
              }
            } else if (file.exists(info_path)) {
              info_content <- paste(readLines(info_path, warn = FALSE), collapse = "\n")
            }
          }

          if (!nzchar(info_content)) {
            info_content <- "No search_info.md found in output directory."
          }

          showModal(modalDialog(
            title = sprintf("Search Info: %s", job$name),
            size = "l", easyClose = TRUE, footer = modalButton("Close"),
            pre(style = "max-height: 500px; overflow-y: auto; font-size: 0.8em; white-space: pre-wrap;",
              info_content)
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
                # Cancel all step IDs for parallel jobs
                if (isTRUE(job$parallel) && !is.null(job$parallel_steps)) {
                  all_ids <- paste(Filter(Negate(is.null), job$parallel_steps), collapse = " ")
                  ssh_exec(cfg, paste(scancel_cmd, all_ids))
                } else {
                  ssh_exec(cfg, paste(scancel_cmd, job$job_id))
                }
              }
            } else {
              if (isTRUE(job$parallel) && !is.null(job$parallel_steps)) {
                all_ids <- Filter(Negate(is.null), job$parallel_steps)
                for (sid in all_ids) {
                  system2("scancel", args = sid, stdout = TRUE, stderr = TRUE)
                }
              } else {
                system2("scancel", args = job$job_id, stdout = TRUE, stderr = TRUE)
              }
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

        # Remove job from queue (mark as removed to preserve indices)
        observeEvent(input[[sprintf("remove_job_%d", idx)]], {
          job <- values$diann_jobs[[idx]]
          jobs <- values$diann_jobs
          jobs[[idx]]$removed <- TRUE
          values$diann_jobs <- jobs
          showNotification(sprintf("Removed job '%s' from queue.", job$name), type = "message")
        }, ignoreInit = TRUE)

        # Resubmit failed/cancelled HPC job
        observeEvent(input[[sprintf("resubmit_job_%d", idx)]], {
          job <- values$diann_jobs[[idx]]

          # --- Smart resume for parallel jobs ---
          if (isTRUE(job$parallel)) {
            cfg <- if (isTRUE(job$is_ssh)) isolate(ssh_config()) else NULL
            if (is.null(cfg)) {
              showNotification("Parallel resubmit requires SSH connection.", type = "error")
              return()
            }

            # Find first failed/cancelled step
            step_status <- job$parallel_step_status %||% list()
            resume_from <- 1L
            for (s in 1:5) {
              st <- step_status[[paste0("step", s)]] %||% "unknown"
              if (st %in% c("completed", "skipped")) next
              resume_from <- s
              break
            }

            # Verify prerequisites exist on remote
            output_dir <- job$output_dir
            # Step 1 was skipped if user provided a speclib — check the original speclib path instead
            step1_skipped <- identical(step_status[["step1"]], "skipped")
            speclib_check <- if (step1_skipped) {
              # User-provided speclib is at its original path (stored in search_settings)
              speclib_path <- job$search_settings$speclib
              if (!is.null(speclib_path) && nzchar(speclib_path %||% "")) {
                sprintf("test -f %s", shQuote(speclib_path))
              } else {
                "true"  # Can't verify, assume OK
              }
            } else {
              sprintf("test -f %s/step1.predicted.speclib", shQuote(output_dir))
            }
            prereq_checks <- list(
              step2 = speclib_check,
              # Step 3 resume: check for backup quant files first, restore if needed
              step3 = sprintf(paste0(
                "if [ -d %s/quant_step2_orig ]; then ",
                "echo 'Restoring Step 2 quant files from backup...'; ",
                "rm -rf %s/quant_step2; ",
                "cp -r %s/quant_step2_orig %s/quant_step2; ",
                "fi && ls %s/quant_step2/*.quant 2>/dev/null | head -1"),
                output_dir, output_dir, output_dir, output_dir, output_dir),
              step4 = sprintf("test -f %s/empirical.parquet", shQuote(output_dir)),
              step5 = sprintf("test -f %s/empirical.parquet && ls %s/quant_step4/*.quant 2>/dev/null | head -1",
                              shQuote(output_dir), output_dir)
            )

            if (resume_from > 1) {
              check_key <- paste0("step", resume_from)
              if (!is.null(prereq_checks[[check_key]])) {
                check <- ssh_exec(cfg, prereq_checks[[check_key]])
                if (check$status != 0 || length(check$stdout) == 0 || !nzchar(check$stdout[1])) {
                  # Fall back: restart from Step 1 (or Step 2 if Step 1 was originally skipped)
                  fallback <- if (step1_skipped) 2L else 1L
                  showNotification(
                    sprintf("Step %d prerequisites not found on remote. Restarting from Step %d.",
                            resume_from, fallback),
                    type = "warning", duration = 8)
                  resume_from <- fallback
                }
              }
            }

            # Build script paths from output_dir
            step_names <- c("step1_libpred.sbatch", "step2_firstpass.sbatch",
                             "step3_assembly.sbatch", "step4_finalpass.sbatch",
                             "step5_report.sbatch")
            step_script_paths <- file.path(output_dir, step_names)

            # Verify scripts exist on remote
            check_cmd <- paste("ls", paste(shQuote(step_script_paths[resume_from:5]), collapse = " "),
                                "2>/dev/null | wc -l")
            check <- ssh_exec(cfg, check_cmd)
            expected <- 5 - resume_from + 1
            if (check$status != 0 || as.integer(trimws(check$stdout[1])) < expected) {
              showNotification("Some sbatch scripts are missing on the remote. Cannot resume.",
                               type = "error", duration = 8)
              return()
            }

            # Generate resume launcher
            sbatch_bin <- values$ssh_sbatch_path %||% "sbatch"
            resume_launcher <- generate_resume_launcher(resume_from, sbatch_bin, step_script_paths)

            # Upload + execute
            resume_file <- tempfile("resume_", fileext = ".sh")
            writeLines(resume_launcher, resume_file)
            on.exit(unlink(resume_file), add = TRUE)

            remote_launcher <- file.path(output_dir, "resume_submit.sh")
            scp_upload(cfg, resume_file, remote_launcher)
            result <- ssh_exec(cfg, paste("bash", shQuote(remote_launcher)))

            if (result$status != 0) {
              showNotification(paste("Resume submission failed:",
                paste(result$stdout, collapse = " ")), type = "error")
              return()
            }

            # Parse step IDs from launcher output
            new_step_ids <- list()
            new_step_status <- list()
            for (line in result$stdout) {
              m <- regexec("^STEP([1-5]):(.+)$", trimws(line))
              if (m[[1]][1] != -1) {
                parts <- regmatches(trimws(line), m)[[1]]
                step_num <- as.integer(parts[2])
                step_val <- trimws(parts[3])
                step_key <- paste0("step", step_num)
                if (step_val == "skipped") {
                  new_step_ids[[step_key]] <- job$parallel_steps[[step_key]]
                  new_step_status[[step_key]] <- "skipped"
                } else {
                  new_step_ids[[step_key]] <- step_val
                  new_step_status[[step_key]] <- "queued"
                }
              }
            }

            if (length(new_step_ids) == 0) {
              showNotification("Could not parse job IDs from resume output.", type = "error")
              return()
            }

            # Create new job entry
            new_entry <- job
            # Use the last submitted step's ID as the main job_id
            last_step_key <- paste0("step", 5)
            new_entry$job_id <- new_step_ids[[last_step_key]] %||%
              new_step_ids[[tail(names(Filter(function(x) x != "skipped",
                new_step_status)), 1)]]
            new_entry$status <- "queued"
            new_entry$parallel_steps <- new_step_ids
            new_entry$parallel_step_status <- new_step_status
            new_entry$parallel_current_step <- resume_from
            new_entry$submitted_at <- Sys.time()
            new_entry$completed_at <- NULL
            new_entry$log_content <- ""
            new_entry$loaded <- FALSE

            values$diann_jobs <- c(values$diann_jobs, list(new_entry))

            skipped_msg <- if (resume_from > 1) {
              sprintf(" (skipped Steps 1-%d, reusing existing results)", resume_from - 1)
            } else ""
            showNotification(
              sprintf("Resumed from Step %d%s", resume_from, skipped_msg),
              type = "message", duration = 10)
            return()
          }

          # --- Single-job resubmit ---
          script_path <- job$script_path

          tryCatch({
            cfg <- if (isTRUE(job$is_ssh)) isolate(ssh_config()) else NULL

            # If script_path is missing, try to recover from output_dir or scontrol
            if (is.null(script_path) || !nzchar(script_path %||% "")) {
              # Try inferring from output_dir
              if (!is.null(job$output_dir) && nzchar(job$output_dir %||% "")) {
                script_path <- file.path(job$output_dir, "diann_search.sbatch")
              }

              # Still missing — try scontrol show job to get Command field
              if (is.null(script_path) || !nzchar(script_path %||% "")) {
                if (!is.null(cfg)) {
                  scontrol_bin <- if (!is.null(values$ssh_sbatch_path)) {
                    file.path(dirname(values$ssh_sbatch_path), "scontrol")
                  } else "scontrol"
                  sctl <- ssh_exec(cfg, paste(scontrol_bin, "show job", job$job_id, "2>/dev/null"),
                                   login_shell = is.null(values$ssh_sbatch_path))
                  if (sctl$status == 0) {
                    cmd_line <- grep("Command=", sctl$stdout, value = TRUE)
                    if (length(cmd_line) > 0) {
                      script_path <- trimws(sub(".*Command=", "", cmd_line[1]))
                    }
                  }
                }
              }

              if (is.null(script_path) || !nzchar(script_path %||% "")) {
                showNotification(
                  "Cannot resubmit: sbatch script path unknown. This job was recovered without full metadata.",
                  type = "error", duration = 10)
                return()
              }
            }

            # Verify script still exists
            if (!is.null(cfg)) {
              check <- ssh_exec(cfg, paste("test -f", shQuote(script_path), "&& echo OK"))
              if (check$status != 0 || !any(grepl("OK", check$stdout))) {
                showNotification("Sbatch script no longer exists on remote.", type = "error")
                return()
              }
            } else {
              if (!file.exists(script_path)) {
                showNotification("Sbatch script no longer exists locally.", type = "error")
                return()
              }
            }

            # Submit via sbatch
            if (!is.null(cfg)) {
              sbatch_bin <- values$ssh_sbatch_path %||% "sbatch"
              sbatch_cmd <- paste(sbatch_bin, shQuote(script_path))
              result <- ssh_exec(cfg, sbatch_cmd,
                                 login_shell = is.null(values$ssh_sbatch_path))
              if (result$status != 0) {
                showNotification(paste("sbatch failed:",
                  paste(result$stdout, collapse = " ")), type = "error")
                return()
              }
              new_job_id <- parse_sbatch_output(result$stdout)
            } else {
              local_sbatch <- Sys.which("sbatch")
              if (!nzchar(local_sbatch)) local_sbatch <- "sbatch"
              stdout <- system2(local_sbatch, args = script_path,
                                stdout = TRUE, stderr = TRUE)
              new_job_id <- parse_sbatch_output(stdout)
            }

            if (is.null(new_job_id)) {
              showNotification("Could not parse new job ID from sbatch output.", type = "error")
              return()
            }

            # Clone job entry with new ID and reset status
            new_entry <- job
            new_entry$job_id <- new_job_id
            new_entry$status <- "queued"
            new_entry$script_path <- script_path
            new_entry$submitted_at <- Sys.time()
            new_entry$completed_at <- NULL
            new_entry$log_content <- ""
            new_entry$loaded <- FALSE

            values$diann_jobs <- c(values$diann_jobs, list(new_entry))
            showNotification(sprintf("Resubmitted as job %s", new_job_id),
              type = "message", duration = 8)
          }, error = function(e) {
            showNotification(sprintf("Resubmit failed: %s", e$message), type = "error")
          })
        }, ignoreInit = TRUE)
      })

      existing <- c(existing, job_key)
    }

    registered_observers(existing)
  })

  # Clear failed/cancelled jobs (mark as removed to preserve observer indices)
  observeEvent(input$clear_failed_jobs, {
    jobs <- values$diann_jobs
    failed_statuses <- c("failed", "cancelled")
    n_removed <- 0L
    for (j in seq_along(jobs)) {
      if (!isTRUE(jobs[[j]]$removed) && jobs[[j]]$status %in% failed_statuses) {
        jobs[[j]]$removed <- TRUE
        n_removed <- n_removed + 1L
      }
    }
    values$diann_jobs <- jobs
    showNotification(sprintf("Cleared %d failed/cancelled job(s).", n_removed), type = "message")
  }, ignoreInit = TRUE)

  # Clear all finished jobs (mark as removed to preserve observer indices)
  observeEvent(input$clear_finished_jobs, {
    jobs <- values$diann_jobs
    terminal <- c("completed", "failed", "cancelled")
    n_removed <- 0L
    for (j in seq_along(jobs)) {
      if (!isTRUE(jobs[[j]]$removed) && jobs[[j]]$status %in% terminal) {
        jobs[[j]]$removed <- TRUE
        n_removed <- n_removed + 1L
      }
    }
    values$diann_jobs <- jobs
    showNotification(sprintf("Cleared %d finished job(s).", n_removed), type = "message")
  }, ignoreInit = TRUE)

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
          sbatch_path = values$ssh_sbatch_path %||%
            (if (nzchar(local_sbatch_path)) local_sbatch_path else NULL),
          days_back = 14
        )
      })

      if (nrow(slurm_jobs) > 0) {
        existing_ids <- if (length(values$diann_jobs) > 0) {
          vapply(values$diann_jobs, function(j) j$job_id %||% "", character(1))
        } else {
          character(0)
        }

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
          # Strategy 0: StdOut from bulk sacct query (most reliable — works for old jobs)
          # Strategy 1: scontrol show job → StdOut path (fallback for recent jobs)
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

          # Strategy 0: Use StdOut from bulk sacct query (most reliable — works for old jobs)
          # StdOut contains the log path template with %j/%A placeholders
          if (nzchar(row$std_out %||% "")) {
            expanded <- row$std_out
            expanded <- gsub("%j", row$job_id, expanded, fixed = TRUE)
            expanded <- gsub("%A", row$job_id, expanded, fixed = TRUE)
            # Verify file exists on cluster
            check_result <- tryCatch(
              run_ssh(sprintf("ls %s 2>/dev/null", shQuote(expanded))),
              error = function(e) list(status = 1, stdout = character()))
            if (check_result$status == 0 && length(check_result$stdout) > 0 &&
                nzchar(trimws(check_result$stdout[1]))) {
              log_file <- trimws(check_result$stdout[1])
              output_dir <- dirname(log_file)
            }
          }

          # Strategy 1: scontrol show job → extract StdOut path (skip if Strategy 0 worked)
          # Use sed instead of grep -oP for portability (not all systems have PCRE grep)
          if (!nzchar(log_file)) {
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
                log_file <- file.path(output_dir, "logs", sprintf("diann_%s.out", row$job_id))
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
              log_parent <- dirname(log_file)
              # If found in logs/ subdir, output_dir is the parent
              output_dir <- if (basename(log_parent) == "logs") dirname(log_parent) else log_parent
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
        existing_ids <- if (length(values$diann_jobs) > 0) {
          vapply(values$diann_jobs, function(j) j$job_id %||% "", character(1))
        } else {
          character(0)
        }
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
