# ==============================================================================
#  server_facility.R
#  Core Facility Mode — Server module
#  Features: Report generation, staff SSH selector, searchable job queue,
#  instrument QC dashboard, template library.
#  Only active when is_core_facility == TRUE.
# ==============================================================================

server_facility <- function(input, output, session, values, add_to_log,
                            is_core_facility, cf_config,
                            search_enabled = FALSE) {

  # Early return if not in core facility mode
  if (!is_core_facility) return(invisible())

  # ============================================================================
  #    Staff SSH Selector (P1)
  # ============================================================================

  # Auto-connect when staff member is selected
  observeEvent(input$staff_selector, {
    req(input$staff_selector, nzchar(input$staff_selector))

    staff <- Filter(
      function(s) s$name == input$staff_selector,
      cf_config$staff$staff
    )
    if (length(staff) == 0) return()
    staff <- staff[[1]]

    # Update SSH fields (hidden but still used by server_search.R)
    updateTextInput(session, "ssh_host", value = staff$hpc_host %||% "")
    updateTextInput(session, "ssh_user", value = staff$hpc_username %||% "")
    updateTextInput(session, "ssh_key_path", value = staff$ssh_key %||% "")
    updateTextInput(session, "diann_account", value = staff$slurm_account %||% "")
    updateTextInput(session, "diann_partition", value = staff$slurm_partition %||% "")

    # Force SSH connection mode (only if HPC backend selected)
    if (!is.null(input$search_backend) && input$search_backend == "hpc") {
      updateRadioButtons(session, "search_connection_mode", selected = "ssh")
    }

    # Auto-test connection
    cfg <- list(
      host = staff$hpc_host %||% "",
      user = staff$hpc_username %||% "",
      port = 22,
      key_path = staff$ssh_key %||% ""
    )

    # Only test if we have the function and required fields
    if (nzchar(cfg$host) && nzchar(cfg$user) && nzchar(cfg$key_path) &&
        exists("test_ssh_connection", mode = "function")) {
      result <- tryCatch(
        test_ssh_connection(cfg),
        error = function(e) list(success = FALSE, error = e$message)
      )

      output$staff_connection_status <- renderUI({
        if (isTRUE(result$success)) {
          tags$div(class = "alert alert-success py-2 px-3",
            style = "font-size: 0.85em;",
            icon("check-circle"),
            sprintf(" Connected to %s as %s", cfg$host, cfg$user)
          )
        } else {
          tags$div(class = "alert alert-danger py-2 px-3",
            style = "font-size: 0.85em;",
            icon("times-circle"),
            sprintf(" Connection failed: %s", result$error %||% "Unknown error"),
            tags$br(),
            tags$small("Check SSH key or ask admin to run ssh-copy-id.")
          )
        }
      })

      values$ssh_connected <- isTRUE(result$success)
    } else {
      output$staff_connection_status <- renderUI({
        tags$div(class = "alert alert-info py-2 px-3",
          style = "font-size: 0.85em;",
          icon("info-circle"),
          sprintf(" Staff: %s (SSH fields updated, use Test Connection to verify)",
                  staff$name)
        )
      })
    }
  })


  # ============================================================================
  #    Project Auto-Populate (staff + lab → previous projects)
  # ============================================================================

  # When staff or lab changes, query DB for existing projects from that combo
  observe({
    staff <- input$staff_selector %||% ""
    lab   <- input$search_lab %||% ""

    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    on.exit(DBI::dbDisconnect(db))

    # Build query: filter by staff+lab if set, otherwise show all projects
    query <- "SELECT DISTINCT project FROM searches WHERE project IS NOT NULL AND project != ''"
    params <- list()
    i <- 1
    if (nzchar(staff)) {
      query <- paste(query, sprintf("AND submitted_by = ?%d", i))
      params[[i]] <- staff
      i <- i + 1
    }
    if (nzchar(lab)) {
      query <- paste(query, sprintf("AND lab = ?%d", i))
      params[[i]] <- lab
      i <- i + 1
    }
    query <- paste(query, "ORDER BY project")

    projects <- if (length(params) > 0) {
      DBI::dbGetQuery(db, query, params = params)$project
    } else {
      DBI::dbGetQuery(db, query)$project
    }

    updateSelectizeInput(session, "search_project",
      choices = projects, selected = character(0),
      server = FALSE)
  }) |> bindEvent(input$staff_selector, input$search_lab, ignoreInit = FALSE)


  # ============================================================================
  #    Shareable HTML Reports (P2)
  # ============================================================================

  # Report generation implementation (called from Job History or modal)
  generate_report_impl <- function(title, lab, instrument, lc_method = "",
                                    project = "", search_id = NULL) {
    if (is.null(values$fit)) {
      showNotification(
        "No DE results loaded. Run the analysis pipeline first.",
        type = "warning", duration = 8)
      return(invisible(NULL))
    }

    analysis_id <- uuid::UUIDgenerate()

    # -- Build state object --
    state <- list(
      raw_data   = NULL,
      metadata   = values$metadata,
      fit        = values$fit,
      y_protein  = values$y_protein,
      dpc_fit    = values$dpc_fit,
      design     = values$design,
      qc_stats   = values$qc_stats,
      gsea_results       = values$gsea_results,
      gsea_results_cache = values$gsea_results_cache,
      gsea_last_contrast = values$gsea_last_contrast,
      gsea_last_org_db   = values$gsea_last_org_db,
      phospho_detected           = values$phospho_detected,
      phospho_site_matrix        = values$phospho_site_matrix,
      phospho_site_info          = values$phospho_site_info,
      phospho_fit                = values$phospho_fit,
      phospho_site_matrix_filtered = values$phospho_site_matrix_filtered,
      phospho_input_mode         = values$phospho_input_mode,
      ksea_results               = values$ksea_results,
      ksea_last_contrast         = values$ksea_last_contrast,
      phospho_annotations        = values$phospho_annotations,
      mofa_object            = values$mofa_object,
      mofa_factors           = values$mofa_factors,
      mofa_variance_explained = values$mofa_variance_explained,
      contrast        = input$contrast_selector,
      logfc_cutoff    = input$logfc_cutoff,
      q_cutoff        = input$q_cutoff,
      color_plot_by_de = values$color_plot_by_de,
      cov1_name       = values$cov1_name,
      cov2_name       = values$cov2_name,
      diann_norm_detected = values$diann_norm_detected,
      metadata_extra = list(
        analysis_id    = analysis_id,
        title          = title,
        lab            = lab,
        instrument     = instrument,
        lc_method      = lc_method,
        project        = project,
        created        = Sys.time(),
        created_by     = input$staff_selector %||% "unknown",
        app_version    = paste0("DE-LIMP v", values$app_version),
        n_proteins     = if (!is.null(values$y_protein)) nrow(values$y_protein) else NA,
        n_contrasts    = if (!is.null(values$fit)) ncol(values$fit$contrasts) else NA,
        contrast_names = if (!is.null(values$fit)) colnames(values$fit$contrasts) else NULL
      ),
      saved_at = Sys.time()
    )

    # -- Save state to disk --
    state_path <- file.path(cf_config$state_dir, paste0(analysis_id, ".rds"))
    saveRDS(state, state_path)

    # -- Render Quarto report --
    html_path <- file.path(cf_config$reports_dir, paste0(analysis_id, ".html"))
    render_success <- FALSE

    withProgress(message = "Generating report...", {
      incProgress(0.2, detail = "Rendering HTML...")

      if (file.exists(cf_config$template_qmd) &&
          requireNamespace("quarto", quietly = TRUE)) {
        render_success <- tryCatch({
          # quarto output_file must be just a filename, not a full path
          quarto::quarto_render(
            input = cf_config$template_qmd,
            output_file = paste0(analysis_id, ".html"),
            execute_params = list(
              analysis_id = analysis_id,
              state_dir = cf_config$state_dir,
              db_path = cf_config$db_path
            )
          )
          # quarto renders next to the input file; move to reports dir
          rendered_at <- file.path(dirname(cf_config$template_qmd),
                                   paste0(analysis_id, ".html"))
          if (file.exists(rendered_at) && rendered_at != html_path) {
            file.rename(rendered_at, html_path)
          }
          TRUE
        }, error = function(e) {
          showNotification(paste("Report render error:", e$message),
                           type = "warning", duration = 10)
          FALSE
        })
      } else {
        showNotification(
          "Quarto template not found. State saved, but HTML report not generated.",
          type = "warning", duration = 8)
      }

      # -- Upload to HF --
      incProgress(0.6, detail = "Uploading for live link...")
      hf_url <- tryCatch({
        upload_state_to_hf(state_path, analysis_id)
      }, error = function(e) {
        showNotification("HF upload failed - state file still saved locally.",
                         type = "warning", duration = 8)
        NULL
      })

      # -- Record in SQLite --
      incProgress(0.8, detail = "Recording...")
      db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
      on.exit(DBI::dbDisconnect(db))

      qc_bracket <- cf_get_bracketing_qc(db, instrument, as.character(Sys.time()))

      DBI::dbExecute(db, "
        INSERT INTO reports (id, title, created_by, lab, instrument,
                             n_proteins, n_contrasts, contrast_names,
                             html_path, state_path, hf_state_url,
                             qc_before_id, qc_after_id)
        VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13)",
        params = list(
          analysis_id, title,
          input$staff_selector %||% "unknown",
          lab, instrument,
          state$metadata_extra$n_proteins,
          state$metadata_extra$n_contrasts,
          paste(state$metadata_extra$contrast_names, collapse = "; "),
          if (render_success) html_path else NA_character_,
          state_path, if (is.null(hf_url)) NA_character_ else hf_url,
          if (nrow(qc_bracket$before) > 0) qc_bracket$before$id[1] else NA_character_,
          if (nrow(qc_bracket$after) > 0) qc_bracket$after$id[1] else NA_character_
        )
      )

      incProgress(1, detail = "Done!")
    })

    # -- Build live link --
    live_link <- if (!is.null(hf_url)) {
      paste0("https://huggingface.co/spaces/brettsp/de-limp-proteomics?analysis=",
             analysis_id)
    } else NULL

    # -- Show result in sidebar --
    output$report_link_ui <- renderUI({
      tags$div(
        style = "margin-top: 10px; padding: 10px; background: #f0f8f0; border-radius: 6px;",
        tags$p(style = "font-weight: 600; margin-bottom: 4px;",
          icon("check-circle"), " Report generated!"),
        if (render_success) {
          tags$p(style = "font-size: 0.85em;",
            "HTML report: ", tags$code(basename(html_path)))
        },
        tags$p(style = "font-size: 0.85em;",
          "State saved: ", tags$code(basename(state_path))),
        if (!is.null(live_link)) {
          tags$p(style = "font-size: 0.85em;",
            "Live link: ",
            tags$a(href = live_link, target = "_blank", live_link))
        },
        if (render_success) {
          downloadButton("download_report_html", "Download HTML Report",
                         class = "btn-outline-success btn-sm w-100 mt-2")
        }
      )
    })

    if (render_success) {
      output$download_report_html <- downloadHandler(
        filename = function() {
          paste0("DE-LIMP_report_", make.names(title), ".html")
        },
        content = function(file) {
          file.copy(html_path, file)
        },
        contentType = "text/html"
      )
    }

    add_to_log("Generate Core Facility Report", c(
      sprintf("# Analysis ID: %s", analysis_id),
      sprintf("# Title: %s", title),
      sprintf("# Lab: %s, Instrument: %s", lab, instrument),
      sprintf("# LC Method: %s, Project: %s", lc_method, project),
      sprintf("# State: %s", state_path),
      if (render_success) sprintf("# HTML: %s", html_path) else "# HTML: not rendered",
      if (!is.null(live_link)) sprintf("# Live link: %s", live_link)
    ))
  }

  # Generate report from Job History — if a job is selected, pull metadata from DB;
  # otherwise show a modal for manual entry
  observeEvent(input$job_generate_report, {
    if (is.null(values$fit)) {
      showNotification(
        "No DE results loaded. Load example data or upload a report, assign groups, and run the pipeline first.",
        type = "warning", duration = 8)
      return()
    }

    sel <- input$job_history_table_rows_selected
    df <- job_history_data()

    if (length(sel) > 0 && !is.null(df) && nrow(df) >= sel) {
      # Job selected — pull metadata from the SQLite record
      job <- df[sel, ]
      generate_report_impl(
        title      = job$analysis_name %||% "Untitled",
        lab        = job$lab %||% "",
        instrument = job$instrument %||% "",
        lc_method  = job$lc_method %||% "",
        project    = job$project %||% "",
        search_id  = job$id
      )
    } else {
      # No job selected — show modal for manual metadata entry
      showModal(modalDialog(
        title = "Generate Report",
        textInput("report_modal_title", "Analysis Title",
          placeholder = "e.g., Smith Lab AP-MS Feb 2026"),
        selectInput("report_modal_lab", "Lab",
          choices = c("(select)" = "", cf_lab_names(cf_config)),
          selected = ""),
        selectInput("report_modal_instrument", "Instrument",
          choices = c("(select)" = "", cf_instrument_names(cf_config)),
          selected = ""),
        selectInput("report_modal_lc_method", "LC / Gradient Method",
          choices = c("(select)" = "", cf_lc_method_names(cf_config)),
          selected = ""),
        selectizeInput("report_modal_project", "Project",
          choices = NULL, selected = NULL,
          options = list(create = TRUE,
                         placeholder = "Select or type new project...")),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirm_report_modal", "Generate", class = "btn-success")
        )
      ))
    }
  })

  # Confirm modal report generation
  observeEvent(input$confirm_report_modal, {
    removeModal()
    generate_report_impl(
      title      = input$report_modal_title %||% "Untitled",
      lab        = input$report_modal_lab %||% "",
      instrument = input$report_modal_instrument %||% "",
      lc_method  = input$report_modal_lc_method %||% "",
      project    = input$report_modal_project %||% ""
    )
  })


  # ============================================================================
  #    Searchable Job Queue (P1)
  # ============================================================================

  # Reactive: query SQLite with filters
  job_history_data <- reactive({
    # Re-query every 30 seconds or when filters change
    invalidateLater(30000, session)
    input$job_search_text
    input$job_filter_lab
    input$job_filter_status
    input$job_filter_staff
    input$job_filter_instrument
    input$job_filter_lc_method

    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    on.exit(DBI::dbDisconnect(db))

    query <- "SELECT * FROM searches WHERE 1=1"
    params <- list()
    i <- 1

    if (nzchar(input$job_search_text %||% "")) {
      query <- paste(query, sprintf("AND analysis_name LIKE ?%d", i))
      params[[i]] <- paste0("%", input$job_search_text, "%")
      i <- i + 1
    }
    if (nzchar(input$job_filter_lab %||% "")) {
      query <- paste(query, sprintf("AND lab = ?%d", i))
      params[[i]] <- input$job_filter_lab
      i <- i + 1
    }
    if (nzchar(input$job_filter_status %||% "")) {
      query <- paste(query, sprintf("AND status = ?%d", i))
      params[[i]] <- input$job_filter_status
      i <- i + 1
    }
    if (nzchar(input$job_filter_staff %||% "")) {
      query <- paste(query, sprintf("AND submitted_by = ?%d", i))
      params[[i]] <- input$job_filter_staff
      i <- i + 1
    }
    if (nzchar(input$job_filter_instrument %||% "")) {
      query <- paste(query, sprintf("AND instrument = ?%d", i))
      params[[i]] <- input$job_filter_instrument
      i <- i + 1
    }
    if (nzchar(input$job_filter_lc_method %||% "")) {
      query <- paste(query, sprintf("AND lc_method = ?%d", i))
      params[[i]] <- input$job_filter_lc_method
      i <- i + 1
    }

    query <- paste(query, "ORDER BY submitted_at DESC LIMIT 500")
    if (length(params) > 0) {
      DBI::dbGetQuery(db, query, params = params)
    } else {
      DBI::dbGetQuery(db, query)
    }
  })

  output$job_history_table <- DT::renderDT({
    df <- job_history_data()
    if (is.null(df) || nrow(df) == 0) {
      return(DT::datatable(data.frame(Message = "No jobs found"),
                           options = list(dom = "t"), rownames = FALSE))
    }

    # Build display columns with graceful fallback for missing columns
    col_names <- c("analysis_name", "lab", "submitted_by",
                    "submitted_at", "status", "backend",
                    "instrument", "lc_method", "project",
                    "n_proteins", "organism")
    # Only include columns that exist in the data
    col_names <- intersect(col_names, names(df))
    display <- df[, col_names, drop = FALSE]

    # Rename for display
    name_map <- c(
      analysis_name = "Name", lab = "Lab", submitted_by = "By",
      submitted_at = "Submitted", status = "Status", backend = "Backend",
      instrument = "Instrument", lc_method = "LC Method", project = "Project",
      n_proteins = "Proteins", organism = "Organism"
    )
    names(display) <- name_map[names(display)]

    # Format timestamps
    if ("Submitted" %in% names(display)) {
      display$Submitted <- tryCatch(
        format(as.POSIXct(display$Submitted), "%b %d %H:%M"),
        error = function(e) display$Submitted
      )
    }

    DT::datatable(display,
      selection = "single",
      options = list(pageLength = 15, scrollX = TRUE, dom = "tip"),
      rownames = FALSE,
      class = "compact stripe hover"
    ) |>
      DT::formatStyle("Status",
        backgroundColor = DT::styleEqual(
          c("completed", "running", "failed", "queued"),
          c("#d4edda", "#cce5ff", "#f8d7da", "#fff3cd")
        )
      )
  })

  # Job count text
  output$job_count_text <- renderText({
    df <- job_history_data()
    if (is.null(df)) return("0 jobs")
    sprintf("%d job%s", nrow(df), if (nrow(df) != 1) "s" else "")
  })

  # Populate filter dropdowns from DB
  observe({
    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    on.exit(DBI::dbDisconnect(db))

    labs <- DBI::dbGetQuery(db,
      "SELECT DISTINCT lab FROM searches WHERE lab IS NOT NULL AND lab != ''")$lab
    staff <- DBI::dbGetQuery(db,
      "SELECT DISTINCT submitted_by FROM searches WHERE submitted_by IS NOT NULL")$submitted_by
    instruments <- DBI::dbGetQuery(db,
      "SELECT DISTINCT instrument FROM searches WHERE instrument IS NOT NULL AND instrument != ''")$instrument
    lc_methods <- DBI::dbGetQuery(db,
      "SELECT DISTINCT lc_method FROM searches WHERE lc_method IS NOT NULL AND lc_method != ''")$lc_method

    updateSelectInput(session, "job_filter_lab",
      choices = c("All labs" = "", labs))
    updateSelectInput(session, "job_filter_staff",
      choices = c("All staff" = "", staff))
    updateSelectInput(session, "job_filter_instrument",
      choices = c("All instruments" = "", instruments))
    updateSelectInput(session, "job_filter_lc_method",
      choices = c("All LC methods" = "", lc_methods))
  })

  # Load results from selected job
  observeEvent(input$job_load_results, {
    sel <- input$job_history_table_rows_selected
    req(sel)
    df <- job_history_data()
    req(nrow(df) >= sel)

    job <- df[sel, ]
    if (job$status != "completed" || is.na(job$output_dir) || !nzchar(job$output_dir)) {
      showNotification("Selected job has no results to load.", type = "warning")
      return()
    }

    # Find report.parquet in the output directory
    report_path <- file.path(job$output_dir, "report.parquet")
    if (!file.exists(report_path)) {
      # Try any report parquet
      parquet_files <- list.files(job$output_dir, pattern = "report.*\\.parquet$",
                                  full.names = TRUE)
      if (length(parquet_files) > 0) {
        report_path <- parquet_files[1]
      } else {
        showNotification(
          sprintf("report.parquet not found in %s", job$output_dir),
          type = "error", duration = 8)
        return()
      }
    }

    # Load the results into the pipeline
    tryCatch({
      withProgress(message = sprintf("Loading results from %s...", job$analysis_name), {
        incProgress(0.1, detail = "Reading parquet...")
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

        incProgress(0.5, detail = "Checking phospho...")
        # Run phospho detection
        tryCatch({
          report_df <- arrow::read_parquet(report_path,
            col_select = c("Modified.Sequence"))
          values$phospho_detected <- detect_phospho(report_df)
        }, error = function(e) NULL)

        incProgress(0.8, detail = "Done!")

        # Log the action
        add_to_log("Load from Job History", c(
          sprintf("# Job: %s (ID: %s)", job$analysis_name, job$id),
          sprintf("# Output dir: %s", job$output_dir),
          sprintf("raw_data <- limpa::readDIANN('%s', format = 'parquet')", report_path)
        ))

        # Navigate to Assign Groups tab
        nav_select("main_tabs", "Data Overview")
        nav_select("data_overview_tabs", "Assign Groups & Run")

        showNotification(
          sprintf("Results loaded from '%s'! Assign groups and run the pipeline.",
                  job$analysis_name),
          type = "message", duration = 10)
      })
    }, error = function(e) {
      showNotification(
        sprintf("Failed to load results: %s", e$message),
        type = "error", duration = 10)
    })
  })

  # (job_generate_report observer is defined above in the Reports section)


  # ============================================================================
  #    Instrument QC Dashboard (P3)
  # ============================================================================

  # Reactive: fetch QC data from SQLite
  qc_dashboard_data <- reactive({
    invalidateLater(60000, session)  # refresh every minute
    input$qc_instrument_filter
    input$qc_date_range

    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    on.exit(DBI::dbDisconnect(db))

    days <- as.integer(input$qc_date_range %||% 90)
    cutoff <- as.character(Sys.time() - days * 86400)

    query <- "SELECT * FROM qc_runs WHERE run_date >= ?1"
    params <- list(cutoff)

    if (nzchar(input$qc_instrument_filter %||% "")) {
      query <- paste(query, "AND instrument = ?2")
      params[[2]] <- input$qc_instrument_filter
    }

    query <- paste(query, "ORDER BY run_date ASC")
    DBI::dbGetQuery(db, query, params = params)
  })

  # Populate instrument filter from DB
  observe({
    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    on.exit(DBI::dbDisconnect(db))
    instruments <- DBI::dbGetQuery(db,
      "SELECT DISTINCT instrument FROM qc_runs")$instrument
    updateSelectInput(session, "qc_instrument_filter",
      choices = c("All Instruments" = "", instruments))
  })

  # -- Protein ID Trend --
  output$qc_protein_trend <- plotly::renderPlotly({
    df <- qc_dashboard_data()
    req(nrow(df) > 0)

    df$run_date <- as.POSIXct(df$run_date)
    df$is_excluded <- !is.na(df$excluded) & df$excluded == 1

    # Compute thresholds from non-excluded runs only
    included_df <- df[!df$is_excluded, ]
    thresholds <- if (nrow(included_df) > 0) {
      included_df |>
        dplyr::group_by(instrument) |>
        dplyr::summarise(
          mean_val = mean(n_proteins, na.rm = TRUE),
          sd_val = stats::sd(n_proteins, na.rm = TRUE),
          lower = mean_val - 2 * sd_val,
          .groups = "drop"
        )
    } else NULL

    # Plot included runs
    p <- plotly::plot_ly()
    if (nrow(included_df) > 0) {
      p <- p |> plotly::add_trace(
        data = included_df,
        x = ~run_date, y = ~n_proteins, color = ~instrument,
        text = ~paste0(run_name, "<br>", n_proteins, " proteins<br>",
                       format(run_date, "%b %d %H:%M")),
        hoverinfo = "text",
        type = "scatter", mode = "markers",
        marker = list(size = 8, opacity = 0.7))
    }

    # Plot excluded runs as faded X markers
    excluded_df <- df[df$is_excluded, ]
    if (nrow(excluded_df) > 0) {
      p <- p |> plotly::add_trace(
        data = excluded_df,
        x = ~run_date, y = ~n_proteins,
        text = ~paste0(run_name, "<br>", n_proteins, " proteins<br>",
                       format(run_date, "%b %d %H:%M"),
                       "<br><b>EXCLUDED</b>",
                       ifelse(!is.na(exclude_reason), paste0("<br>", exclude_reason), "")),
        hoverinfo = "text",
        type = "scatter", mode = "markers",
        marker = list(size = 8, opacity = 0.25, symbol = "x",
                       color = "gray"),
        name = "Excluded",
        showlegend = TRUE)
    }

    p <- p |> plotly::layout(
      title = list(text = "Protein Identifications", font = list(size = 14)),
      xaxis = list(title = ""),
      yaxis = list(title = "Proteins"),
      legend = list(orientation = "h", y = -0.15),
      hovermode = "closest"
    )

    # Add threshold lines per instrument
    if (!is.null(thresholds)) {
      for (i in seq_len(nrow(thresholds))) {
        if (!is.na(thresholds$lower[i])) {
          p <- p |> plotly::add_segments(
            x = min(df$run_date), xend = max(df$run_date),
            y = thresholds$lower[i], yend = thresholds$lower[i],
            line = list(color = "red", dash = "dash", width = 1),
            showlegend = FALSE,
            inherit = FALSE
          )
        }
      }
    }

    p
  })

  # -- Precursor ID Trend --
  output$qc_precursor_trend <- plotly::renderPlotly({
    df <- qc_dashboard_data()
    req(nrow(df) > 0)
    df$run_date <- as.POSIXct(df$run_date)
    df$is_excluded <- !is.na(df$excluded) & df$excluded == 1

    included_df <- df[!df$is_excluded, ]
    p <- plotly::plot_ly()
    if (nrow(included_df) > 0) {
      p <- p |> plotly::add_trace(
        data = included_df,
        x = ~run_date, y = ~n_precursors, color = ~instrument,
        text = ~paste0(run_name, "<br>", n_precursors, " precursors"),
        hoverinfo = "text",
        type = "scatter", mode = "markers",
        marker = list(size = 8, opacity = 0.7))
    }
    excluded_df <- df[df$is_excluded, ]
    if (nrow(excluded_df) > 0) {
      p <- p |> plotly::add_trace(
        data = excluded_df,
        x = ~run_date, y = ~n_precursors,
        text = ~paste0(run_name, "<br>", n_precursors, " precursors<br><b>EXCLUDED</b>"),
        hoverinfo = "text",
        type = "scatter", mode = "markers",
        marker = list(size = 8, opacity = 0.25, symbol = "x", color = "gray"),
        name = "Excluded", showlegend = FALSE)
    }
    p |> plotly::layout(
      title = list(text = "Precursor Identifications", font = list(size = 14)),
      xaxis = list(title = ""),
      yaxis = list(title = "Precursors"),
      legend = list(orientation = "h", y = -0.15)
    )
  })

  # -- MS1 TIC Trend --
  output$qc_tic_trend <- plotly::renderPlotly({
    df <- qc_dashboard_data()
    req(nrow(df) > 0)
    df$run_date <- as.POSIXct(df$run_date)
    df <- df[!is.na(df$median_ms1_tic), ]
    req(nrow(df) > 0)
    df$is_excluded <- !is.na(df$excluded) & df$excluded == 1

    included_df <- df[!df$is_excluded, ]
    p <- plotly::plot_ly()
    if (nrow(included_df) > 0) {
      p <- p |> plotly::add_trace(
        data = included_df,
        x = ~run_date, y = ~median_ms1_tic, color = ~instrument,
        text = ~paste0(run_name, "<br>TIC: ", signif(median_ms1_tic, 3)),
        hoverinfo = "text",
        type = "scatter", mode = "markers",
        marker = list(size = 8, opacity = 0.7))
    }
    excluded_df <- df[df$is_excluded, ]
    if (nrow(excluded_df) > 0) {
      p <- p |> plotly::add_trace(
        data = excluded_df,
        x = ~run_date, y = ~median_ms1_tic,
        text = ~paste0(run_name, "<br>TIC: ", signif(median_ms1_tic, 3),
                        "<br><b>EXCLUDED</b>"),
        hoverinfo = "text",
        type = "scatter", mode = "markers",
        marker = list(size = 8, opacity = 0.25, symbol = "x", color = "gray"),
        name = "Excluded", showlegend = FALSE)
    }
    p |> plotly::layout(
      title = list(text = "Median MS1 TIC", font = list(size = 14)),
      xaxis = list(title = ""),
      yaxis = list(title = "TIC", type = "log"),
      legend = list(orientation = "h", y = -0.15)
    )
  })

  # -- QC Runs Table --
  output$qc_runs_table <- DT::renderDT({
    df <- qc_dashboard_data()
    req(nrow(df) > 0)

    display <- df |>
      dplyr::arrange(dplyr::desc(run_date)) |>
      dplyr::transmute(
        Instrument = instrument,
        File = run_name,
        Date = format(as.POSIXct(run_date), "%b %d %H:%M"),
        ng = ng_loaded,
        Gradient = gradient,
        Proteins = n_proteins,
        Precursors = n_precursors,
        `MS1 TIC` = signif(median_ms1_tic, 3),
        Status = ifelse(!is.na(excluded) & excluded == 1,
          paste0("Excluded", ifelse(!is.na(exclude_reason) & nzchar(exclude_reason),
                                     paste0(": ", exclude_reason), "")),
          "")
      )

    dt <- DT::datatable(display,
      selection = "multiple",
      options = list(pageLength = 10, scrollX = TRUE, dom = "tip"),
      rownames = FALSE,
      class = "compact stripe hover"
    )
    # Highlight rows with non-empty Status (excluded runs)
    excluded_rows <- which(nzchar(display$Status))
    if (length(excluded_rows) > 0) {
      dt <- dt |> DT::formatStyle("Status",
        backgroundColor = DT::styleEqual(
          display$Status[excluded_rows],
          rep("#fff3cd", length(excluded_rows))
        ),
        color = DT::styleEqual(
          display$Status[excluded_rows],
          rep("#856404", length(excluded_rows))
        )
      )
    }
    dt
  })


  # ============================================================================
  #    QC Run Ingestion (two-mode modal: Scan & Submit / Ingest Existing)
  # ============================================================================

  # Reactive to hold scan results
  qc_scan_results <- reactiveVal(NULL)

  # QC modal — two-mode: Scan & Submit vs Ingest Existing
  observeEvent(input$qc_ingest_btn, {
    instruments <- cf_instrument_names(cf_config)
    if (length(instruments) == 0) instruments <- c("timsTOF HT", "Exploris 480")

    qc_search_cfg <- cf_config$qc$qc_search
    has_qc_search <- !is.null(qc_search_cfg)

    showModal(modalDialog(
      title = "QC Run Management",
      size = "l",
      # Mode toggle
      radioButtons("qc_modal_mode", NULL,
        choices = c("Scan & Submit" = "scan", "Ingest Existing" = "ingest"),
        selected = if (has_qc_search) "scan" else "ingest",
        inline = TRUE),

      # --- Scan & Submit mode ---
      conditionalPanel(
        condition = "input.qc_modal_mode == 'scan'",
        if (has_qc_search) {
          tagList(
            div(style = "display: flex; gap: 10px; align-items: flex-end; margin-bottom: 10px;",
              div(style = "flex: 0 0 200px;",
                selectizeInput("qc_scan_instrument", "Instrument",
                  choices = instruments, options = list(create = TRUE))
              ),
              div(style = "flex: 1;",
                textInput("qc_scan_dir", "Remote Directory",
                  placeholder = "/path/to/instrument/data/")
              ),
              div(style = "flex: 0 0 auto; padding-bottom: 0px;",
                actionButton("qc_scan_btn", "Scan", icon = icon("search"),
                  class = "btn-primary btn-sm")
              )
            ),
            tags$small(class = "text-muted",
              "Scans for *hela* and *qc* files (.d, .raw) in the specified directory."),
            hr(),
            uiOutput("qc_scan_results_ui"),
            div(style = "margin-top: 10px;",
              actionButton("qc_submit_searches", "Submit DIA-NN Searches",
                icon = icon("rocket"), class = "btn-success",
                disabled = TRUE)
            )
          )
        } else {
          tags$div(class = "alert alert-info",
            icon("info-circle"),
            " QC search configuration not found in qc_config.yml. ",
            "Add a ", tags$code("qc_search"), " section to enable scan & submit.")
        }
      ),

      # --- Ingest Existing mode ---
      conditionalPanel(
        condition = "input.qc_modal_mode == 'ingest'",
        tags$p("Record QC metrics from an existing DIA-NN report.parquet (e.g., HeLa digest standard)."),
        textInput("qc_ingest_path", "Path to report.parquet",
          placeholder = "/path/to/qc_output/report.parquet"),
        selectizeInput("qc_ingest_instrument", "Instrument",
          choices = instruments, options = list(create = TRUE)),
        div(style = "display: flex; gap: 10px;",
          div(style = "flex: 1;",
            textInput("qc_ingest_run_name", "Run Name (optional)",
              placeholder = "Auto-detected from directory name")
          ),
          div(style = "flex: 0 0 120px;",
            numericInput("qc_ingest_ng", "ng Loaded", value = NA, min = 1)
          ),
          div(style = "flex: 0 0 150px;",
            textInput("qc_ingest_gradient", "Gradient",
              placeholder = "e.g., 60SPD")
          )
        )
      ),

      footer = tagList(
        modalButton("Cancel"),
        conditionalPanel(
          condition = "input.qc_modal_mode == 'ingest'",
          style = "display: inline;",
          actionButton("confirm_qc_ingest", "Record QC Metrics",
            icon = icon("check"), class = "btn-success")
        )
      )
    ))
  })

  # Scan for QC files via SSH
  observeEvent(input$qc_scan_btn, {
    req(input$qc_scan_dir, nzchar(input$qc_scan_dir))

    cfg <- isolate(ssh_config())
    if (is.null(cfg) || !isTRUE(values$ssh_connected)) {
      showNotification("SSH not connected. Select a staff member first.",
        type = "warning")
      return()
    }

    scan_dir <- trimws(input$qc_scan_dir)

    withProgress(message = "Scanning for QC files...", {
      # Find .d directories and .raw files matching HeLa/QC patterns
      scan_cmd <- sprintf(
        "find %s -maxdepth 1 \\( -name '*.d' -o -name '*.raw' \\) | grep -i -E '(hela|qc)' | sort",
        shQuote(scan_dir))
      result <- ssh_exec(cfg, scan_cmd)

      if (result$status != 0 || length(result$stdout) == 0 ||
          !any(nzchar(trimws(result$stdout)))) {
        showNotification("No QC files found in that directory.", type = "warning")
        qc_scan_results(NULL)
        return()
      }

      files <- trimws(result$stdout)
      files <- files[nzchar(files)]

      # Check which files are already ingested
      db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
      existing <- DBI::dbGetQuery(db,
        "SELECT DISTINCT run_name FROM qc_runs")$run_name
      DBI::dbDisconnect(db)

      # Build scan results data frame
      scan_df <- data.frame(
        path = files,
        filename = basename(files),
        stringsAsFactors = FALSE
      )

      # Parse metadata from filenames
      parsed <- lapply(scan_df$filename, parse_qc_filename)
      scan_df$ng_loaded <- vapply(parsed, function(p) {
        if (is.na(p$ng_loaded)) NA_integer_ else p$ng_loaded
      }, integer(1))
      scan_df$gradient <- vapply(parsed, function(p) {
        if (is.na(p$gradient)) NA_character_ else p$gradient
      }, character(1))

      scan_df$status <- ifelse(scan_df$filename %in% existing, "Ingested", "New")
      scan_df$selected <- scan_df$status == "New"

      qc_scan_results(scan_df)
    })
  })

  # Render scan results table
  output$qc_scan_results_ui <- renderUI({
    df <- qc_scan_results()
    if (is.null(df) || nrow(df) == 0) {
      return(tags$p(class = "text-muted", "No scan results. Click 'Scan' to search for QC files."))
    }

    n_new <- sum(df$status == "New")

    # Enable/disable submit button
    shinyjs_available <- requireNamespace("shinyjs", quietly = TRUE)

    tagList(
      tags$p(sprintf("Found %d file%s (%d new, %d already ingested)",
        nrow(df), if (nrow(df) != 1) "s" else "",
        n_new, sum(df$status == "Ingested"))),

      # Editable table of scan results
      tags$div(style = "max-height: 300px; overflow-y: auto;",
        tags$table(class = "table table-sm table-striped",
          tags$thead(
            tags$tr(
              tags$th(style = "width: 40px;", ""),
              tags$th("File"),
              tags$th(style = "width: 60px;", "Status"),
              tags$th(style = "width: 100px;", "ng Loaded"),
              tags$th(style = "width: 120px;", "Gradient")
            )
          ),
          tags$tbody(
            lapply(seq_len(nrow(df)), function(i) {
              is_new <- df$status[i] == "New"
              tags$tr(
                style = if (!is_new) "opacity: 0.5;" else NULL,
                tags$td(
                  checkboxInput(
                    paste0("qc_sel_", i), NULL,
                    value = is_new, width = "30px")
                ),
                tags$td(df$filename[i]),
                tags$td(
                  if (is_new) {
                    tags$span(class = "badge bg-success", "New")
                  } else {
                    tags$span(class = "badge bg-secondary", "Ingested")
                  }
                ),
                tags$td(
                  numericInput(paste0("qc_ng_", i), NULL,
                    value = df$ng_loaded[i], min = 1, width = "90px")
                ),
                tags$td(
                  textInput(paste0("qc_grad_", i), NULL,
                    value = df$gradient[i] %||% "", width = "110px")
                )
              )
            })
          )
        )
      )
    )
  })

  # Submit QC DIA-NN searches
  observeEvent(input$qc_submit_searches, {
    df <- qc_scan_results()
    req(df, nrow(df) > 0)

    qc_search_cfg <- cf_config$qc$qc_search
    req(qc_search_cfg)

    cfg <- isolate(ssh_config())
    req(cfg, values$ssh_connected)

    instrument <- input$qc_scan_instrument %||% ""

    # Collect selected files with their metadata
    selected_files <- list()
    for (i in seq_len(nrow(df))) {
      if (isTRUE(input[[paste0("qc_sel_", i)]])) {
        selected_files[[length(selected_files) + 1]] <- list(
          path = df$path[i],
          filename = df$filename[i],
          ng_loaded = input[[paste0("qc_ng_", i)]],
          gradient = input[[paste0("qc_grad_", i)]] %||% ""
        )
      }
    }

    if (length(selected_files) == 0) {
      showNotification("No files selected for submission.", type = "warning")
      return()
    }

    # SLURM settings from config
    slurm_cfg <- qc_search_cfg$slurm %||% list()
    cpus <- slurm_cfg$cpus %||% 32
    mem_gb <- slurm_cfg$mem_gb %||% 128
    time_hours <- slurm_cfg$time_hours %||% 4
    account <- input$diann_account %||% "genome-center-grp"
    partition <- input$diann_partition %||% "high"

    # Search params from config
    sp <- qc_search_cfg$default_search_params %||% list()
    search_params <- list(
      qvalue = sp$qvalue %||% 0.01,
      max_var_mods = sp$max_var_mods %||% 1,
      mbr = sp$mbr %||% TRUE,
      rt_profiling = sp$rt_profiling %||% FALSE,
      xic = sp$xic %||% FALSE,
      enzyme = sp$enzyme %||% "K*,R*",
      missed_cleavages = sp$missed_cleavages %||% 1
    )

    withProgress(message = "Submitting QC searches...",
      value = 0, max = length(selected_files), {

      jobs <- values$diann_jobs %||% list()

      for (idx in seq_along(selected_files)) {
        f <- selected_files[[idx]]
        analysis_name <- paste0("QC_", tools::file_path_sans_ext(f$filename))
        date_str <- format(Sys.Date(), "%Y-%m-%d")
        output_dir <- file.path(qc_search_cfg$output_base, instrument,
                                 date_str, analysis_name)

        incProgress(1, detail = sprintf("Submitting %s...", f$filename))

        tryCatch({
          # Create output directory
          ssh_exec(cfg, paste("mkdir -p", shQuote(output_dir)))

          # Generate sbatch script (library-only, no FASTA)
          script <- generate_sbatch_script(
            analysis_name = analysis_name,
            raw_files = f$path,
            fasta_files = character(0),
            speclib_path = qc_search_cfg$hela_library,
            output_dir = output_dir,
            diann_sif = qc_search_cfg$diann_sif,
            search_mode = "library",
            cpus = cpus,
            mem_gb = mem_gb,
            time_hours = time_hours,
            partition = partition,
            account = account,
            search_params = search_params
          )

          # Upload sbatch script via SCP
          script_local <- tempfile(fileext = ".sh")
          writeLines(script, script_local)
          script_remote <- file.path(output_dir,
            paste0("qc_search_", analysis_name, ".sh"))
          scp_upload(cfg, script_local, script_remote)
          unlink(script_local)

          # Submit via sbatch
          sbatch_path <- values$ssh_sbatch_path %||% "sbatch"
          submit_result <- ssh_exec(cfg, paste(sbatch_path, shQuote(script_remote)))

          if (submit_result$status != 0) {
            showNotification(
              sprintf("Failed to submit %s: %s", f$filename,
                paste(submit_result$stderr, collapse = " ")),
              type = "error", duration = 10)
            next
          }

          # Parse SLURM job ID
          slurm_id <- regmatches(
            paste(submit_result$stdout, collapse = " "),
            regexpr("\\d+", paste(submit_result$stdout, collapse = " "))
          )
          if (length(slurm_id) == 0) slurm_id <- paste0("qc_", idx)

          # Add to job queue with QC metadata
          jobs[[length(jobs) + 1]] <- list(
            job_id = slurm_id,
            name = analysis_name,
            status = "queued",
            output_dir = output_dir,
            n_files = 1,
            auto_load = FALSE,
            backend = "hpc",
            is_ssh = TRUE,
            qc_run = TRUE,
            qc_instrument = instrument,
            qc_ng_loaded = f$ng_loaded,
            qc_gradient = if (nzchar(f$gradient)) f$gradient else NA_character_,
            qc_raw_file = f$path,
            submitted_at = Sys.time()
          )

          # Record in facility DB
          if (!is.null(cf_config$db_path)) {
            tryCatch(
              cf_record_search(cf_config$db_path, list(
                analysis_name = analysis_name,
                submitted_by = input$staff_selector %||% "unknown",
                lab = "",
                instrument = instrument,
                lc_method = f$gradient,
                project = "QC",
                organism = "Human",
                fasta_file = "",
                n_raw_files = 1L,
                search_mode = "library",
                slurm_job_id = slurm_id,
                backend = "hpc",
                output_dir = output_dir
              )),
              error = function(e) message("CF DB record failed: ", e$message)
            )
          }
        }, error = function(e) {
          showNotification(
            sprintf("Error submitting %s: %s", f$filename, e$message),
            type = "error", duration = 10)
        })
      }

      values$diann_jobs <- jobs
    })

    removeModal()
    showNotification(
      sprintf("Submitted %d QC search%s!", length(selected_files),
        if (length(selected_files) != 1) "es" else ""),
      type = "message", duration = 10)

    add_to_log("Submit QC Searches", c(
      sprintf("# Instrument: %s", instrument),
      sprintf("# Files: %d", length(selected_files)),
      sprintf("# Library: %s", qc_search_cfg$hela_library),
      vapply(selected_files, function(f) sprintf("#   %s (%s ng, %s)",
        f$filename, f$ng_loaded %||% "?", f$gradient %||% "?"), character(1))
    ))
  })

  # Ingest Existing — confirm handler
  observeEvent(input$confirm_qc_ingest, {
    req(input$qc_ingest_path)
    req(input$qc_ingest_instrument)

    report_path <- trimws(input$qc_ingest_path)
    instrument <- input$qc_ingest_instrument
    run_name <- if (nzchar(input$qc_ingest_run_name %||% "")) {
      input$qc_ingest_run_name
    } else NULL
    ng_loaded <- input$qc_ingest_ng
    gradient <- if (nzchar(input$qc_ingest_gradient %||% "")) {
      input$qc_ingest_gradient
    } else NULL

    if (!file.exists(report_path)) {
      showNotification(
        sprintf("File not found: %s", report_path),
        type = "error", duration = 8)
      return()
    }

    # Check if already ingested
    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    existing <- DBI::dbGetQuery(db, "SELECT id FROM qc_runs WHERE file_path = ?1",
      params = list(report_path))
    DBI::dbDisconnect(db)

    if (nrow(existing) > 0) {
      showNotification("This file has already been ingested.", type = "warning")
      return()
    }

    tryCatch({
      withProgress(message = "Ingesting QC metrics...", {
        row_id <- cf_ingest_qc_metrics(
          db_path = cf_config$db_path,
          report_path = report_path,
          instrument = instrument,
          run_name = run_name,
          ng_loaded = ng_loaded,
          gradient = gradient
        )

        add_to_log("Ingest QC Run", c(
          sprintf("# Instrument: %s", instrument),
          sprintf("# File: %s", report_path),
          sprintf("# ng: %s, Gradient: %s", ng_loaded %||% "NA", gradient %||% "NA"),
          sprintf("cf_ingest_qc_metrics('%s', '%s', '%s')",
                  cf_config$db_path, report_path, instrument)
        ))
      })

      removeModal()
      showNotification(
        sprintf("QC run ingested for %s!", instrument),
        type = "message", duration = 8)
    }, error = function(e) {
      showNotification(
        sprintf("QC ingestion failed: %s", e$message),
        type = "error", duration = 10)
    })
  })

  # ============================================================================
  #    QC Exclude / Re-include
  # ============================================================================

  # Exclude selected QC runs
  observeEvent(input$qc_exclude_btn, {
    sel <- input$qc_runs_table_rows_selected
    req(length(sel) > 0)

    showModal(modalDialog(
      title = "Exclude QC Runs",
      tags$p(sprintf("Exclude %d selected run%s from trend analysis?",
        length(sel), if (length(sel) != 1) "s" else "")),
      textInput("qc_exclude_reason", "Reason (optional)",
        placeholder = "e.g., Loading error, column degraded"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_qc_exclude", "Exclude",
          icon = icon("ban"), class = "btn-warning")
      )
    ))
  })

  observeEvent(input$confirm_qc_exclude, {
    sel <- input$qc_runs_table_rows_selected
    req(length(sel) > 0)

    df <- qc_dashboard_data()
    req(nrow(df) >= max(sel))

    # Get IDs in display order (descending by run_date)
    display_df <- df |> dplyr::arrange(dplyr::desc(run_date))
    run_ids <- display_df$id[sel]

    tryCatch({
      cf_update_qc_excluded(cf_config$db_path, run_ids,
        excluded = TRUE, reason = input$qc_exclude_reason)
      removeModal()
      showNotification(
        sprintf("Excluded %d QC run%s.", length(sel),
          if (length(sel) != 1) "s" else ""),
        type = "message", duration = 5)
    }, error = function(e) {
      showNotification(sprintf("Exclude failed: %s", e$message),
        type = "error", duration = 8)
    })
  })

  # Re-include selected QC runs
  observeEvent(input$qc_include_btn, {
    sel <- input$qc_runs_table_rows_selected
    req(length(sel) > 0)

    df <- qc_dashboard_data()
    req(nrow(df) >= max(sel))

    display_df <- df |> dplyr::arrange(dplyr::desc(run_date))
    run_ids <- display_df$id[sel]

    tryCatch({
      cf_update_qc_excluded(cf_config$db_path, run_ids, excluded = FALSE)
      showNotification(
        sprintf("Re-included %d QC run%s.", length(sel),
          if (length(sel) != 1) "s" else ""),
        type = "message", duration = 5)
    }, error = function(e) {
      showNotification(sprintf("Re-include failed: %s", e$message),
        type = "error", duration = 8)
    })
  })


  # ============================================================================
  #    Template Library (P4)
  # ============================================================================

  # Populate template dropdown
  observe({
    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    on.exit(DBI::dbDisconnect(db))
    templates <- DBI::dbGetQuery(db,
      "SELECT id, name, type FROM templates ORDER BY name")
    choices <- c("(none)" = "")
    if (nrow(templates) > 0) {
      choices <- c(choices,
        setNames(templates$id,
                 paste0(templates$name, " [", templates$type, "]")))
    }
    updateSelectInput(session, "template_selector", choices = choices)
  })

  # Save current settings as template
  observeEvent(input$save_template, {
    showModal(modalDialog(
      title = "Save as Template",
      textInput("template_name", "Template Name",
        placeholder = "e.g., Standard AP-MS"),
      selectInput("template_type", "Type",
        choices = c("DIA-NN Search Preset" = "search_preset",
                    "Group Assignment" = "group_template",
                    "GSEA Configuration" = "gsea_config")),
      textAreaInput("template_notes", "Notes", rows = 2),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_save_template", "Save", class = "btn-primary")
      )
    ))
  })

  observeEvent(input$confirm_save_template, {
    req(input$template_name, nzchar(input$template_name))

    config <- switch(input$template_type,
      "search_preset" = jsonlite::toJSON(list(
        cpus = input$diann_cpus,
        mem_gb = input$diann_mem_gb,
        time_hours = input$diann_time_hours,
        partition = input$diann_partition,
        account = input$diann_account,
        mass_acc_mode = input$mass_acc_mode,
        mass_acc = input$diann_mass_acc,
        ms1_acc = input$diann_mass_acc_ms1,
        search_mode = input$search_mode,
        enzyme = input$diann_enzyme,
        missed_cleavages = input$diann_missed_cleavages,
        max_var_mods = input$diann_max_var_mods,
        normalization = input$diann_normalization
      ), auto_unbox = TRUE),
      "group_template" = jsonlite::toJSON(list(
        metadata = values$metadata
      ), auto_unbox = TRUE),
      "gsea_config" = jsonlite::toJSON(list(
        ontology = input$gsea_ontology,
        logfc_cutoff = input$logfc_cutoff
      ), auto_unbox = TRUE)
    )

    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    on.exit(DBI::dbDisconnect(db))
    DBI::dbExecute(db, "
      INSERT INTO templates (name, type, created_by, config_json, notes)
      VALUES (?1, ?2, ?3, ?4, ?5)",
      params = list(
        input$template_name,
        input$template_type,
        input$staff_selector %||% "unknown",
        config,
        input$template_notes %||% ""
      )
    )

    removeModal()
    showNotification("Template saved!", type = "message")
  })

  # Load (apply) a template
  observeEvent(input$load_template, {
    req(input$template_selector, nzchar(input$template_selector))

    db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
    on.exit(DBI::dbDisconnect(db))

    tpl <- DBI::dbGetQuery(db,
      "SELECT * FROM templates WHERE id = ?1",
      params = list(as.integer(input$template_selector))
    )
    req(nrow(tpl) > 0)

    config <- tryCatch(
      jsonlite::fromJSON(tpl$config_json[1]),
      error = function(e) NULL
    )
    req(!is.null(config))

    if (tpl$type[1] == "search_preset") {
      # Apply search settings
      if (!is.null(config$cpus))
        updateNumericInput(session, "diann_cpus", value = config$cpus)
      if (!is.null(config$mem_gb))
        updateNumericInput(session, "diann_mem_gb", value = config$mem_gb)
      if (!is.null(config$time_hours))
        updateNumericInput(session, "diann_time_hours", value = config$time_hours)
      if (!is.null(config$partition))
        updateTextInput(session, "diann_partition", value = config$partition)
      if (!is.null(config$account))
        updateTextInput(session, "diann_account", value = config$account)
      if (!is.null(config$mass_acc_mode))
        updateSelectInput(session, "mass_acc_mode", selected = config$mass_acc_mode)
      if (!is.null(config$mass_acc))
        updateNumericInput(session, "diann_mass_acc", value = config$mass_acc)
      if (!is.null(config$ms1_acc))
        updateNumericInput(session, "diann_mass_acc_ms1", value = config$ms1_acc)
      if (!is.null(config$search_mode))
        updateRadioButtons(session, "search_mode", selected = config$search_mode)
      if (!is.null(config$enzyme))
        updateSelectInput(session, "diann_enzyme", selected = config$enzyme)
      if (!is.null(config$missed_cleavages))
        updateNumericInput(session, "diann_missed_cleavages", value = config$missed_cleavages)
      if (!is.null(config$max_var_mods))
        updateNumericInput(session, "diann_max_var_mods", value = config$max_var_mods)
      if (!is.null(config$normalization))
        updateRadioButtons(session, "diann_normalization", selected = config$normalization)

      showNotification(sprintf("Applied search preset: %s", tpl$name[1]),
                       type = "message")

    } else if (tpl$type[1] == "group_template") {
      # Restore metadata
      if (!is.null(config$metadata)) {
        values$metadata <- config$metadata
        showNotification(sprintf("Applied group template: %s", tpl$name[1]),
                         type = "message")
      }

    } else if (tpl$type[1] == "gsea_config") {
      if (!is.null(config$ontology))
        updateSelectInput(session, "gsea_ontology", selected = config$ontology)
      if (!is.null(config$logfc_cutoff))
        updateSliderInput(session, "logfc_cutoff", value = config$logfc_cutoff)

      showNotification(sprintf("Applied GSEA config: %s", tpl$name[1]),
                       type = "message")
    }
  })

}
