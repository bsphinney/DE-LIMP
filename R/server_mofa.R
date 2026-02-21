# ==============================================================================
#  server_mofa.R
#  Multi-View Integration (MOFA2) — view management, training, results.
#  Called from app.R as: server_mofa(input, output, session, values, add_to_log)
# ==============================================================================

server_mofa <- function(input, output, session, values, add_to_log) {

  # Detect HF Spaces for resource limits
  is_hf <- nzchar(Sys.getenv("SPACE_ID", ""))

  # ============================================================================
  #      Tab Content (rendered dynamically)
  # ============================================================================

  output$mofa_tab_content <- renderUI({
    tagList(
      # Header
      div(style = "background: linear-gradient(135deg, #2c3e50, #3498db); color: white; padding: 15px 20px; border-radius: 8px; margin-bottom: 15px;",
        div(style = "display: flex; justify-content: space-between; align-items: center;",
          div(
            tags$h4("Multi-Omics MOFA2", style = "margin: 0;"),
            tags$p("Unsupervised discovery of shared and view-specific variation patterns",
                   style = "margin: 5px 0 0 0; opacity: 0.9; font-size: 0.9em;")
          ),
          div(style = "display: flex; gap: 8px; align-items: center;",
            actionButton("load_mofa_example", "Mouse Brain (2-view)",
                         class = "btn-info btn-sm", icon = icon("download")),
            actionButton("load_mofa_tcga", "TCGA Breast (3-view)",
                         class = "btn-info btn-sm", icon = icon("download")),
            actionButton("mofa_info_btn", icon("question-circle"),
                         class = "btn-outline-light btn-sm", title = "About MOFA")
          )
        )
      ),

      # HF warning
      if (is_hf) {
        div(class = "alert alert-warning py-2 mb-3",
          icon("cloud"),
          "Running on Hugging Face Spaces with limited resources. ",
          "For large datasets or thorough training, consider local installation."
        )
      },

      # Two-column layout: setup + results
      fluidRow(
        # Left column: View configuration + parameters
        column(5,
          # View cards
          div(style = "border: 1px solid #dee2e6; border-radius: 8px; padding: 15px; margin-bottom: 15px; background: #f8f9fa;",
            tags$h5("Data Views (2-6 supported)", style = "margin-top: 0;"),
            uiOutput("mofa_view_cards")
          ),

          # Sample matching
          uiOutput("mofa_sample_matching"),

          # Parameters
          div(style = "border: 1px solid #dee2e6; border-radius: 8px; padding: 15px; margin-bottom: 15px;",
            div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
              tags$h5("MOFA Parameters", style = "margin: 0;"),
              actionButton("mofa_params_info", icon("question-circle"),
                           class = "btn-outline-info btn-sm", title = "Parameter help")
            ),
            fluidRow(
              column(6,
                numericInput("mofa_num_factors", "Number of factors:",
                             value = if (is_hf) 10 else 15,
                             min = 2, max = if (is_hf) 10 else 50, step = 1)
              ),
              column(6,
                div(style = "margin-top: 30px;",
                  checkboxInput("mofa_auto_factors", "Auto-select (data-driven)", value = TRUE)
                )
              )
            ),
            radioButtons("mofa_convergence", "Convergence mode:",
              choices = c(
                "Fast (~2 min)" = "fast",
                "Medium (~10 min)" = "medium",
                "Thorough (~30+ min)" = "slow"
              ),
              selected = if (is_hf) "fast" else "medium"
            ),
            tags$small(class = "text-muted", "Faster = fewer iterations, may miss subtle factors"),

            tags$details(
              tags$summary("Advanced options"),
              class = "mt-3",
              sliderInput("mofa_min_variance", "Min variance explained:",
                          min = 0, max = 0.05, value = 0.01, step = 0.005,
                          post = " (drop factors below this)"),
              checkboxInput("mofa_scale_views", "Scale views (unit variance)", value = TRUE),
              tags$small(class = "text-muted d-block mb-2",
                         "Recommended ON for multi-omics (equalizes view contributions)"),
              checkboxInput("mofa_center_features", "Center features (zero mean)", value = TRUE),
              numericInput("mofa_seed", "Random seed:", value = 42, min = 1, max = 99999)
            ),

            tags$hr(),
            actionButton("run_mofa", "Train MOFA Model",
                         class = "btn-warning btn-lg w-100",
                         icon = icon("play")),

            # Template download
            div(style = "margin-top: 10px; text-align: center;",
              downloadButton("download_mofa_template", "Download template CSV",
                             class = "btn-sm btn-outline-secondary")
            )
          )
        ),

        # Right column: Results
        column(7,
          uiOutput("mofa_results_panel")
        )
      )
    )
  })

  # ============================================================================
  #      View Card UI Generation
  # ============================================================================

  # Auto-populate View 1 from current analysis
  observe({
    req(values$fit)
    if (length(values$mofa_view_configs) == 0) {
      values$mofa_view_configs <- list(
        list(
          id = "v1",
          num = 1,
          name = "Global Proteomics",
          type = "proteomics_other",
          source = "current",
          matrix = values$raw_data,
          fit = values$fit,
          n_features = nrow(values$raw_data),
          n_samples = ncol(values$raw_data),
          status = "ready"
        )
      )
      # Also register as a loaded view
      values$mofa_views[["Global Proteomics"]] <- values$raw_data
      values$mofa_view_fits[["v1"]] <- values$fit
    }
  })

  # Render all view cards
  output$mofa_view_cards <- renderUI({
    configs <- values$mofa_view_configs
    if (length(configs) < 1) {
      return(div(class = "text-muted", "Load data and run the DE pipeline first to populate View 1."))
    }

    tagList(
      lapply(seq_along(configs), function(i) {
        config <- configs[[i]]
        view_id <- config$id
        is_first <- (i == 1)

        div(style = "border: 1px solid #ccc; border-radius: 6px; padding: 12px; margin-bottom: 10px; background: white;",
          # Header row
          div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 8px;",
            tags$span(
              tags$strong(paste("VIEW", config$num, ":")),
              tags$span(config$name, style = "margin-left: 5px;")
            ),
            if (!is_first) {
              actionButton(paste0("mofa_remove_view_", view_id),
                label = NULL, icon = icon("times"),
                class = "btn-sm btn-outline-danger p-1")
            } else {
              tags$span(class = "badge bg-secondary", "Required")
            }
          ),

          # View name
          textInput(paste0("mofa_view_name_", view_id), label = NULL,
                    value = config$name,
                    placeholder = "e.g., Global Proteomics, RNA-seq"),

          # Source and upload (non-View 1 only)
          if (!is_first) {
            tagList(
              selectInput(paste0("mofa_view_type_", view_id), "Data type:",
                choices = c(
                  "Proteomics (DIA-NN)" = "proteomics_diann",
                  "Proteomics (other)" = "proteomics_other",
                  "Phosphoproteomics" = "phospho",
                  "Transcriptomics (RNA-seq)" = "transcriptomics",
                  "Metabolomics" = "metabolomics",
                  "Lipidomics" = "lipidomics",
                  "Other quantitative matrix" = "other"
                ),
                selected = config$type
              ),
              uiOutput(paste0("mofa_view_source_ui_", view_id)),
              uiOutput(paste0("mofa_view_upload_ui_", view_id))
            )
          } else {
            # View 1 status
            if (config$status == "ready") {
              div(class = "alert alert-success py-1 mb-0 small",
                icon("check"),
                sprintf("Current analysis: %d features x %d samples",
                        config$n_features, config$n_samples)
              )
            }
          },

          # Status indicator
          uiOutput(paste0("mofa_view_status_", view_id))
        )
      }),

      # Add view button
      if (length(configs) < 6) {
        div(class = "d-flex align-items-center gap-3 mt-3",
          actionButton("mofa_add_view", "+ Add View",
                       class = "btn-outline-primary", icon = icon("plus")),
          tags$span(class = "text-muted",
                    sprintf("(%d/6 views configured)", length(configs)))
        )
      } else {
        div(class = "alert alert-info py-2 mt-3",
          icon("info-circle"), "Maximum 6 views reached."
        )
      }
    )
  })

  # ============================================================================
  #      Dynamic View Management
  # ============================================================================

  # Add View handler
  observeEvent(input$mofa_add_view, {
    configs <- values$mofa_view_configs
    n <- length(configs)
    if (n >= 6) {
      showNotification("Maximum 6 views supported", type = "warning")
      return()
    }
    new_id <- paste0("v", n + 1)
    new_config <- default_view_config(new_id, n + 1)
    configs[[n + 1]] <- new_config
    values$mofa_view_configs <- configs
    showNotification(sprintf("Added View %d", n + 1), type = "message", duration = 3)
  })

  # Remove View handlers (Views 2-6)
  lapply(2:6, function(i) {
    view_id <- paste0("v", i)
    observeEvent(input[[paste0("mofa_remove_view_", view_id)]], {
      configs <- values$mofa_view_configs
      # Find and remove
      configs <- configs[sapply(configs, function(x) x$id != view_id)]
      # Renumber
      for (j in seq_along(configs)) {
        configs[[j]]$num <- j
      }
      # Clean up loaded data
      old_name <- NULL
      for (vc in values$mofa_view_configs) {
        if (vc$id == view_id) old_name <- vc$name
      }
      if (!is.null(old_name) && old_name %in% names(values$mofa_views)) {
        values$mofa_views[[old_name]] <- NULL
      }
      values$mofa_view_fits[[view_id]] <- NULL
      values$mofa_view_configs <- configs
      showNotification("View removed", type = "message", duration = 3)
    }, ignoreInit = TRUE)
  })

  # ============================================================================
  #      Source / Upload UI for Views 2-6
  # ============================================================================

  # Create source + upload UI for each possible view
  lapply(2:6, function(i) {
    view_id <- paste0("v", i)

    # Source radio buttons
    output[[paste0("mofa_view_source_ui_", view_id)]] <- renderUI({
      vtype <- input[[paste0("mofa_view_type_", view_id)]]
      req(vtype)

      if (vtype == "proteomics_diann") {
        radioButtons(paste0("mofa_view_source_", view_id), "Source:",
          choices = c(
            "Upload pre-processed matrix (.csv, .tsv)" = "upload_matrix",
            "Load DE-LIMP session (.rds)" = "upload_rds"
          ),
          selected = "upload_rds"
        )
      } else if (vtype == "phospho") {
        has_phospho <- !is.null(values$phospho_site_matrix_filtered)
        choices <- c(
          "Upload phosphosite matrix (.csv, .tsv)" = "upload_matrix",
          "Load DE-LIMP session (.rds)" = "upload_rds"
        )
        if (has_phospho) {
          choices <- c("From Phospho tab" = "phospho_tab", choices)
        }
        radioButtons(paste0("mofa_view_source_", view_id), "Source:",
          choices = choices,
          selected = if (has_phospho) "phospho_tab" else "upload_matrix"
        )
      } else {
        radioButtons(paste0("mofa_view_source_", view_id), "Source:",
          choices = c(
            "Upload matrix (.csv, .tsv, .parquet)" = "upload_matrix",
            "Upload R object (.rds)" = "upload_rds"
          ),
          selected = "upload_matrix"
        )
      }
    })

    # File upload UI
    output[[paste0("mofa_view_upload_ui_", view_id)]] <- renderUI({
      src <- input[[paste0("mofa_view_source_", view_id)]]
      req(src)

      if (src == "phospho_tab") {
        if (!is.null(values$phospho_site_matrix_filtered)) {
          div(class = "alert alert-success py-2",
            icon("check"),
            sprintf("Phospho data available: %d sites x %d samples",
                    nrow(values$phospho_site_matrix_filtered),
                    ncol(values$phospho_site_matrix_filtered))
          )
        } else {
          div(class = "alert alert-warning py-2",
            icon("exclamation-triangle"),
            "No phospho data available. Run phosphosite analysis first."
          )
        }
      } else if (src == "upload_rds") {
        tagList(
          fileInput(paste0("mofa_view_file_", view_id), "Upload RDS file:",
                    accept = ".rds"),
          tags$small(class = "text-muted",
            icon("info-circle"),
            "Accepts: DE-LIMP session, limma object, or matrix."
          )
        )
      } else {
        tagList(
          fileInput(paste0("mofa_view_file_", view_id), "Upload matrix:",
                    accept = c(".csv", ".tsv", ".txt", ".parquet")),
          tags$small(class = "text-muted",
            "Format: Features (rows) x Samples (columns). First col = feature IDs."
          )
        )
      }
    })

    # Status indicator
    output[[paste0("mofa_view_status_", view_id)]] <- renderUI({
      configs <- values$mofa_view_configs
      this_config <- NULL
      for (vc in configs) {
        if (vc$id == view_id) this_config <- vc
      }
      if (is.null(this_config) || this_config$status == "pending") return(NULL)

      if (this_config$status == "ready") {
        div(class = "alert alert-success py-1 mb-0 small mt-2",
          icon("check"),
          sprintf("Ready: %d features x %d samples",
                  this_config$n_features, this_config$n_samples)
        )
      } else if (this_config$status == "error") {
        div(class = "alert alert-danger py-1 mb-0 small mt-2",
          icon("times-circle"), "Error loading data"
        )
      }
    })

    # File upload handler
    observeEvent(input[[paste0("mofa_view_file_", view_id)]], {
      file_info <- input[[paste0("mofa_view_file_", view_id)]]
      req(file_info)

      src <- input[[paste0("mofa_view_source_", view_id)]]
      view_name <- input[[paste0("mofa_view_name_", view_id)]] %||% paste("View", i)

      tryCatch({
        result <- if (src == "upload_rds") {
          parse_rds_for_mofa(file_info$datapath)
        } else {
          parse_matrix_file(file_info$datapath, file_info$name)
        }

        mat <- result$matrix

        # Store in views
        values$mofa_views[[view_name]] <- mat

        # Store fit if available (from RDS)
        if (!is.null(result$fit)) {
          values$mofa_view_fits[[view_id]] <- result$fit
        }

        # Update config
        configs <- values$mofa_view_configs
        for (j in seq_along(configs)) {
          if (configs[[j]]$id == view_id) {
            configs[[j]]$name <- view_name
            configs[[j]]$matrix <- mat
            configs[[j]]$fit <- result$fit
            configs[[j]]$n_features <- nrow(mat)
            configs[[j]]$n_samples <- ncol(mat)
            configs[[j]]$status <- "ready"
            configs[[j]]$file_name <- file_info$name
            break
          }
        }
        values$mofa_view_configs <- configs

        showNotification(result$message, type = "message", duration = 5)

      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error")
        configs <- values$mofa_view_configs
        for (j in seq_along(configs)) {
          if (configs[[j]]$id == view_id) {
            configs[[j]]$status <- "error"
            break
          }
        }
        values$mofa_view_configs <- configs
      })
    })

    # Phospho tab source handler
    observeEvent(input[[paste0("mofa_view_source_", view_id)]], {
      src <- input[[paste0("mofa_view_source_", view_id)]]
      if (is.null(src) || src != "phospho_tab") return()
      req(values$phospho_site_matrix_filtered)

      view_name <- input[[paste0("mofa_view_name_", view_id)]] %||% "Phosphoproteomics"
      mat <- values$phospho_site_matrix_filtered

      values$mofa_views[[view_name]] <- mat

      configs <- values$mofa_view_configs
      for (j in seq_along(configs)) {
        if (configs[[j]]$id == view_id) {
          configs[[j]]$name <- view_name
          configs[[j]]$matrix <- mat
          configs[[j]]$n_features <- nrow(mat)
          configs[[j]]$n_samples <- ncol(mat)
          configs[[j]]$status <- "ready"
          break
        }
      }
      values$mofa_view_configs <- configs

      showNotification(sprintf("Loaded phospho data: %d sites x %d samples",
                               nrow(mat), ncol(mat)), type = "message")
    }, ignoreInit = TRUE)

    # View name change -> update the views list key
    observeEvent(input[[paste0("mofa_view_name_", view_id)]], {
      new_name <- input[[paste0("mofa_view_name_", view_id)]]
      req(nzchar(new_name))

      configs <- values$mofa_view_configs
      for (j in seq_along(configs)) {
        if (configs[[j]]$id == view_id) {
          old_name <- configs[[j]]$name
          if (old_name != new_name && old_name %in% names(values$mofa_views)) {
            values$mofa_views[[new_name]] <- values$mofa_views[[old_name]]
            values$mofa_views[[old_name]] <- NULL
          }
          configs[[j]]$name <- new_name
          break
        }
      }
      values$mofa_view_configs <- configs
    }, ignoreInit = TRUE)
  })

  # View 1 status (always present)
  output$mofa_view_status_v1 <- renderUI(NULL)

  # ============================================================================
  #      Sample Matching Panel
  # ============================================================================

  output$mofa_sample_matching <- renderUI({
    views <- values$mofa_views
    if (length(views) < 2) return(NULL)

    # Filter out NULL views
    valid_views <- views[!sapply(views, is.null)]
    if (length(valid_views) < 2) return(NULL)

    overlap <- compute_sample_overlap(valid_views)

    # Status styling
    if (overlap$n_common == overlap$n_total) {
      status_class <- "alert-success"
      status_msg <- sprintf("All %d samples present in all views", overlap$n_common)
    } else if (overlap$pct_overlap >= 80) {
      status_class <- "alert-info"
      status_msg <- sprintf("%d/%d samples matched (%.0f%%)",
                            overlap$n_common, overlap$n_total, overlap$pct_overlap)
    } else if (overlap$n_common >= 3) {
      status_class <- "alert-warning"
      status_msg <- sprintf("Only %d/%d samples matched -- check sample naming",
                            overlap$n_common, overlap$n_total)
    } else {
      status_class <- "alert-danger"
      status_msg <- "Insufficient sample overlap -- cannot train MOFA"
    }

    view_breakdown <- lapply(overlap$per_view, function(pv) {
      tags$li(sprintf("%s: %d samples (%d shared, %d unique)",
                      pv$name, pv$total, pv$shared, pv$unique))
    })

    div(style = "border: 1px solid #dee2e6; border-radius: 8px; padding: 15px; margin-bottom: 15px;",
      tags$h5("Sample Matching", style = "margin-top: 0;"),
      div(class = paste("alert", status_class, "py-2"), status_msg),
      tags$details(
        tags$summary("View sample breakdown"),
        tags$ul(view_breakdown)
      )
    )
  })

  # ============================================================================
  #      MOFA Training
  # ============================================================================

  observeEvent(input$run_mofa, {
    views <- values$mofa_views
    valid_views <- views[!sapply(views, is.null)]

    if (length(valid_views) < 2) {
      showNotification("At least 2 data views are required to train MOFA", type = "error")
      return()
    }

    # Find common samples
    overlap <- compute_sample_overlap(valid_views)
    if (overlap$n_common < 3) {
      showNotification("Need at least 3 common samples across views", type = "error")
      return()
    }

    # Check MOFA2 available
    if (!requireNamespace("MOFA2", quietly = TRUE)) {
      showNotification("MOFA2 package not installed. Run: BiocManager::install('MOFA2')",
                       type = "error", duration = 10)
      return()
    }

    # Capture parameters before entering withProgress (input$ not always available inside)
    scale_views <- isTRUE(input$mofa_scale_views)
    auto_factors <- isTRUE(input$mofa_auto_factors)
    num_factors <- if (auto_factors) { if (is_hf) 10 else 15 } else { input$mofa_num_factors }
    convergence_mode <- input$mofa_convergence %||% "medium"
    seed_val <- input$mofa_seed %||% 42
    min_var <- input$mofa_min_variance %||% 0.01

    withProgress(message = "Training MOFA model...", {

      tryCatch({

        # Step 1: Prepare data
        incProgress(0.1, detail = "Preparing data matrices...")
        common_samples <- overlap$common_samples

        mofa_data <- lapply(valid_views, function(mat) {
          mat[, common_samples, drop = FALSE]
        })

        # Step 2: Train in isolated subprocess
        # CRITICAL: basilisk's Python process management crashes the R session
        # when run inside Shiny's httpuv event loop. Using callr::r() runs the
        # entire MOFA pipeline (create → configure → train) in a separate R
        # process, completely isolating basilisk/Python from Shiny. The trained
        # model is saved to HDF5, then loaded back in the main process.
        incProgress(0.2, detail = "Training model in subprocess (this may take several minutes)...")

        outfile <- tempfile(fileext = ".hdf5")

        tryCatch({
          callr::r(
            function(mofa_data, outfile, scale_views, num_factors,
                     convergence_mode, seed_val) {
              library(MOFA2)

              mofa_obj <- create_mofa(mofa_data)

              data_opts <- get_default_data_options(mofa_obj)
              data_opts$scale_views <- scale_views

              model_opts <- get_default_model_options(mofa_obj)
              model_opts$num_factors <- num_factors

              train_opts <- get_default_training_options(mofa_obj)
              train_opts$convergence_mode <- convergence_mode
              train_opts$seed <- seed_val
              train_opts$verbose <- FALSE

              mofa_obj <- prepare_mofa(mofa_obj,
                data_options = data_opts,
                model_options = model_opts,
                training_options = train_opts
              )

              run_mofa(mofa_obj, outfile = outfile, use_basilisk = TRUE)
            },
            args = list(
              mofa_data = mofa_data,
              outfile = outfile,
              scale_views = scale_views,
              num_factors = num_factors,
              convergence_mode = convergence_mode,
              seed_val = seed_val
            ),
            show = FALSE,
            timeout = if (is_hf) 600 else 1800  # 10 min HF, 30 min local
          )
        }, error = function(ce) {
          # Extract the real error from callr's wrapper
          msg <- conditionMessage(ce)
          # callr errors have $parent with the original error
          if (!is.null(ce$parent)) msg <- conditionMessage(ce$parent)
          # Also check for stderr output
          if (!is.null(ce$stderr) && nzchar(ce$stderr)) {
            msg <- paste0(msg, "\nSubprocess stderr: ",
                          substr(ce$stderr, max(1, nchar(ce$stderr) - 500), nchar(ce$stderr)))
          }
          stop(msg, call. = FALSE)
        })

        # Step 3: Load trained model back from HDF5 (no basilisk needed)
        incProgress(0.85, detail = "Loading trained model...")
        mofa_trained <- MOFA2::load_model(outfile)

        # Step 4: Post-process
        incProgress(0.9, detail = "Extracting results...")

        # Drop factors below variance threshold
        if (min_var > 0) {
          r2 <- MOFA2::get_variance_explained(mofa_trained)$r2_total[[1]]
          keep_factors <- names(r2)[r2 >= min_var * 100]
          if (length(keep_factors) > 0 && length(keep_factors) < length(r2)) {
            mofa_trained <- MOFA2::subset_factors(mofa_trained, factors = keep_factors)
            showNotification(sprintf("Dropped %d factors below %.1f%% variance threshold",
                                     length(r2) - length(keep_factors),
                                     min_var * 100), type = "message")
          }
        }

        # Store results
        values$mofa_object <- mofa_trained
        values$mofa_factors <- MOFA2::get_factors(mofa_trained)[[1]]
        values$mofa_weights <- MOFA2::get_weights(mofa_trained)
        values$mofa_variance_explained <- MOFA2::get_variance_explained(mofa_trained)

        # Clean up temp file
        unlink(outfile)

        n_factors <- ncol(values$mofa_factors)
        view_names <- names(valid_views)

        values$mofa_last_run_params <- list(
          n_views = length(valid_views),
          n_samples = length(common_samples),
          n_factors = n_factors,
          view_names = view_names,
          convergence = convergence_mode,
          scale_views = scale_views,
          seed = seed_val,
          timestamp = Sys.time()
        )

        showNotification(
          sprintf("MOFA training complete: %d factors across %d views",
                  n_factors, length(valid_views)),
          type = "message", duration = 10
        )

        add_to_log("MOFA Training", generate_mofa_code(values$mofa_last_run_params))

      }, error = function(e) {
        showNotification(paste("MOFA training failed:", e$message),
                         type = "error", duration = 15)
        values$mofa_object <- NULL
      })
    })
  })

  # ============================================================================
  #      Results Panel
  # ============================================================================

  output$mofa_results_panel <- renderUI({
    if (is.null(values$mofa_object)) {
      return(div(
        style = "border: 1px solid #dee2e6; border-radius: 8px; padding: 40px; text-align: center; color: #6c757d;",
        icon("layer-group", style = "font-size: 3em; margin-bottom: 15px;"),
        tags$h5("No MOFA results yet"),
        tags$p("Configure at least 2 views, then click 'Train MOFA Model' to discover latent factors.")
      ))
    }

    params <- values$mofa_last_run_params
    div(
      # Results header
      div(class = "alert alert-success py-2 mb-3",
        sprintf("Model trained: %d views, %d samples, %d active factors",
                params$n_views, params$n_samples, params$n_factors)
      ),

      # Results tabs
      navset_card_tab(
        id = "mofa_results_tabs",

        nav_panel("Variance Explained",
          plotOutput("mofa_variance_heatmap", height = "400px"),
          tags$details(class = "mt-2",
            tags$summary("How to interpret"),
            tags$ul(
              tags$li(tags$b("Factor loads heavily on one view:"), " View-specific variation (e.g., technical batch, unique biology)"),
              tags$li(tags$b("Factor loads similarly across views:"), " Shared variation (e.g., treatment effect)"),
              tags$li(tags$b("Phospho-heavy factor:"), " Signaling/kinase activity independent of abundance"),
              tags$li(tags$b("Global-heavy factor:"), " Expression-level changes driving phospho (abundance effect)")
            )
          )
        ),

        nav_panel("Factor Weights",
          fluidRow(
            column(3, uiOutput("mofa_weights_controls")),
            column(9, plotly::plotlyOutput("mofa_weights_plot", height = "500px"))
          )
        ),

        nav_panel("Sample Scores",
          fluidRow(
            column(3, uiOutput("mofa_scores_controls")),
            column(9, plotly::plotlyOutput("mofa_scores_plot", height = "500px"))
          )
        ),

        nav_panel("Top Features",
          fluidRow(
            column(3,
              uiOutput("mofa_table_controls")
            ),
            column(9,
              DT::DTOutput("mofa_top_features_table")
            )
          )
        ),

        nav_panel("Factor-DE Correlation",
          uiOutput("mofa_de_controls"),
          plotOutput("mofa_de_correlation", height = "500px")
        )
      )
    )
  })

  # ============================================================================
  #      Tab 1: Variance Explained Heatmap
  # ============================================================================

  output$mofa_variance_heatmap <- renderPlot({
    req(values$mofa_variance_explained)

    r2 <- values$mofa_variance_explained$r2_per_factor

    # r2_per_factor is a list of groups, each a matrix (factors x views)
    var_list <- list()
    for (group_name in names(r2)) {
      mat <- r2[[group_name]]
      for (i in seq_len(nrow(mat))) {
        for (j in seq_len(ncol(mat))) {
          var_list[[length(var_list) + 1]] <- data.frame(
            View = colnames(mat)[j],
            Factor = rownames(mat)[i],
            Variance = mat[i, j],
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(var_list) == 0) return(NULL)
    var_df <- do.call(rbind, var_list)

    # Preserve natural factor order (Factor1, Factor2, ...)
    factor_levels <- unique(var_df$Factor)
    factor_levels <- factor_levels[order(as.integer(gsub("\\D+", "", factor_levels)))]
    var_df$Factor <- factor(var_df$Factor, levels = factor_levels)

    ggplot(var_df, aes(x = Factor, y = View, fill = Variance)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.1f%%", Variance)), size = 3.5) +
      scale_fill_gradient2(low = "white", mid = "#457B9D", high = "#1D3557",
                           midpoint = max(var_df$Variance, na.rm = TRUE) / 2,
                           name = "Variance\nExplained (%)") +
      labs(
        title = "Variance Explained per View per Factor",
        subtitle = "Factors capturing view-specific vs shared variation",
        x = NULL, y = NULL
      ) +
      theme_minimal(base_size = 14) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold")
      )
  }, height = 400)

  # ============================================================================
  #      Tab 2: Factor Weights Browser
  # ============================================================================

  output$mofa_weights_controls <- renderUI({
    req(values$mofa_weights)

    view_names <- names(values$mofa_weights)
    factor_names <- colnames(values$mofa_weights[[1]])

    tagList(
      selectInput("mofa_weight_view", "View:", choices = view_names),
      selectInput("mofa_weight_factor", "Factor:", choices = factor_names),
      sliderInput("mofa_weight_n_top", "Top N features:",
                  min = 10, max = 50, value = 20, step = 5),
      tags$hr(),
      downloadButton("download_mofa_weights", "Download all weights (.csv)")
    )
  })

  output$mofa_weights_plot <- plotly::renderPlotly({
    req(values$mofa_weights, input$mofa_weight_view, input$mofa_weight_factor)

    view <- input$mofa_weight_view
    factor <- input$mofa_weight_factor

    req(view %in% names(values$mofa_weights))
    w_mat <- values$mofa_weights[[view]]
    req(factor %in% colnames(w_mat))

    weights <- w_mat[, factor]
    n_top <- input$mofa_weight_n_top %||% 20
    top_idx <- order(abs(weights), decreasing = TRUE)[1:min(n_top, length(weights))]

    df <- data.frame(
      Feature = names(weights)[top_idx],
      Weight = weights[top_idx],
      Direction = ifelse(weights[top_idx] > 0, "Positive", "Negative"),
      stringsAsFactors = FALSE
    )
    df$Feature <- factor(df$Feature, levels = df$Feature[order(df$Weight)])

    p <- ggplot(df, aes(x = Weight, y = Feature, fill = Direction)) +
      geom_col() +
      scale_fill_manual(values = c("Positive" = "#E63946", "Negative" = "#457B9D")) +
      labs(
        title = sprintf("Top %d Features for %s", n_top, factor),
        subtitle = paste("View:", view),
        x = "Weight", y = NULL
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")

    plotly::ggplotly(p, tooltip = c("x", "y")) %>%
      plotly::layout(hoverlabel = list(bgcolor = "white"))
  })

  # Download all weights
  output$download_mofa_weights <- downloadHandler(
    filename = function() {
      paste0("mofa_weights_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$mofa_weights)
      all_weights <- lapply(names(values$mofa_weights), function(view) {
        w <- as.data.frame(values$mofa_weights[[view]])
        w$Feature <- rownames(w)
        w$View <- view
        w
      })
      combined <- do.call(rbind, all_weights)
      write.csv(combined, file, row.names = FALSE)
    }
  )

  # ============================================================================
  #      Tab 3: Sample Factor Scores
  # ============================================================================

  output$mofa_scores_controls <- renderUI({
    req(values$mofa_factors)
    factor_names <- colnames(values$mofa_factors)

    tagList(
      selectInput("mofa_score_x", "X-axis factor:",
                  choices = factor_names, selected = factor_names[1]),
      selectInput("mofa_score_y", "Y-axis factor:",
                  choices = factor_names,
                  selected = if (length(factor_names) >= 2) factor_names[2] else factor_names[1]),
      tags$hr(),
      tags$small(class = "text-muted",
        "Samples cluster by shared biology captured in these factors. ",
        "Points colored by experimental group from your metadata."),
      tags$hr(),
      downloadButton("download_mofa_scores", "Download factor scores (.csv)")
    )
  })

  output$mofa_scores_plot <- plotly::renderPlotly({
    req(values$mofa_factors, input$mofa_score_x, input$mofa_score_y)

    factors_df <- as.data.frame(values$mofa_factors)
    factors_df$Sample <- rownames(factors_df)

    # Add group coloring from metadata
    # Try MOFA-specific metadata first (from test data), then main DE pipeline metadata
    if (!is.null(values$mofa_sample_metadata) && "Sample" %in% names(values$mofa_sample_metadata) &&
        "Group" %in% names(values$mofa_sample_metadata)) {
      factors_df <- merge(factors_df, values$mofa_sample_metadata[, c("Sample", "Group")],
                          by = "Sample", all.x = TRUE)
    } else if (!is.null(values$metadata) && "Run" %in% names(values$metadata) &&
               "Group" %in% names(values$metadata)) {
      factors_df <- merge(factors_df, values$metadata[, c("Run", "Group")],
                          by.x = "Sample", by.y = "Run", all.x = TRUE)
    } else {
      factors_df$Group <- "Unknown"
    }

    x_col <- input$mofa_score_x
    y_col <- input$mofa_score_y
    req(x_col %in% names(factors_df), y_col %in% names(factors_df))

    p <- ggplot(factors_df, aes(x = .data[[x_col]], y = .data[[y_col]],
                                 color = Group, label = Sample)) +
      geom_point(size = 3, alpha = 0.8) +
      labs(
        title = sprintf("%s vs %s", x_col, y_col),
        subtitle = "Sample positions in latent factor space"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")

    plotly::ggplotly(p, tooltip = c("label", "color", "x", "y"))
  })

  # Download factor scores
  output$download_mofa_scores <- downloadHandler(
    filename = function() {
      paste0("mofa_factor_scores_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$mofa_factors)
      df <- as.data.frame(values$mofa_factors)
      df$Sample <- rownames(df)
      write.csv(df, file, row.names = FALSE)
    }
  )

  # ============================================================================
  #      Tab 4: Top Features per Factor (Table)
  # ============================================================================

  output$mofa_table_controls <- renderUI({
    req(values$mofa_weights)
    factor_names <- colnames(values$mofa_weights[[1]])

    tagList(
      selectInput("mofa_table_factor", "Factor:", choices = factor_names),
      sliderInput("mofa_table_n_top", "Top N per view:",
                  min = 10, max = 100, value = 50, step = 10)
    )
  })

  output$mofa_top_features_table <- DT::renderDataTable({
    req(values$mofa_weights, input$mofa_table_factor)

    factor <- input$mofa_table_factor
    n_top <- input$mofa_table_n_top %||% 50

    all_top <- lapply(names(values$mofa_weights), function(view) {
      w <- values$mofa_weights[[view]][, factor]
      top_idx <- order(abs(w), decreasing = TRUE)[1:min(n_top, length(w))]
      data.frame(
        View = view,
        Feature = names(w)[top_idx],
        Weight = round(w[top_idx], 4),
        AbsWeight = round(abs(w[top_idx]), 4),
        stringsAsFactors = FALSE
      )
    })

    combined <- do.call(rbind, all_top)
    combined <- combined[order(-combined$AbsWeight), ]
    combined$Rank <- seq_len(nrow(combined))

    DT::datatable(
      combined[, c("Rank", "View", "Feature", "Weight")],
      options = list(
        pageLength = 25,
        dom = "Bfrtip",
        buttons = c("csv", "excel")
      ),
      extensions = "Buttons",
      rownames = FALSE,
      caption = sprintf("Top features driving %s across all views", factor)
    )
  })

  # ============================================================================
  #      Tab 5: Factor-DE Correlation
  # ============================================================================

  output$mofa_de_controls <- renderUI({
    req(values$mofa_weights)

    # Check if any fit object is available (main pipeline or imported from RDS views)
    has_main_fit <- !is.null(values$fit) && !is.null(values$fit$contrasts)
    has_view_fits <- any(sapply(values$mofa_view_fits, function(f) !is.null(f)))

    if (!has_main_fit && !has_view_fits) {
      return(div(class = "alert alert-info mt-3",
        icon("info-circle"),
        "Factor-DE correlation requires differential expression results. ",
        "Run the DE pipeline on your data first (Data Overview tab), or ",
        "upload a DE-LIMP session (.rds) as one of your MOFA views to include its fit object."
      ))
    }

    # Collect all available fit objects
    fits_available <- character()
    if (has_main_fit) {
      fits_available <- c("Main DE pipeline" = "current")
    }
    for (vid in names(values$mofa_view_fits)) {
      if (!is.null(values$mofa_view_fits[[vid]])) {
        for (vc in values$mofa_view_configs) {
          if (vc$id == vid) {
            fits_available[vc$name] <- vid
            break
          }
        }
      }
    }

    # Get contrast names from first available fit
    first_fit <- if (has_main_fit) values$fit else values$mofa_view_fits[[names(which(sapply(values$mofa_view_fits, function(f) !is.null(f))))[1]]]
    contrast_names <- colnames(first_fit$contrasts)

    div(class = "mb-3",
      fluidRow(
        column(4,
          selectInput("mofa_de_contrast", "DE contrast:", choices = contrast_names)
        ),
        column(4,
          selectInput("mofa_de_view", "Weights from view:",
                      choices = names(values$mofa_weights))
        ),
        column(4,
          selectInput("mofa_de_fit", "Fit object:", choices = fits_available)
        )
      )
    )
  })

  output$mofa_de_correlation <- renderPlot({
    req(values$mofa_weights, input$mofa_de_contrast, input$mofa_de_view)

    # Get the fit object
    fit_source <- input$mofa_de_fit %||% "current"
    fit_obj <- if (fit_source == "current") {
      values$fit
    } else {
      values$mofa_view_fits[[fit_source]]
    }
    req(fit_obj)

    # Get DE results
    contrast <- input$mofa_de_contrast
    req(contrast %in% colnames(fit_obj$contrasts))

    de_res <- limma::topTable(fit_obj, coef = contrast, number = Inf)

    # Get weights for selected view
    view <- input$mofa_de_view
    req(view %in% names(values$mofa_weights))
    view_weights <- values$mofa_weights[[view]]

    # Match features
    common_features <- intersect(rownames(de_res), rownames(view_weights))
    if (length(common_features) < 10) {
      plot(NULL, xlim = c(0, 1), ylim = c(0, 1),
           main = "Insufficient feature overlap for correlation",
           xlab = "", ylab = "")
      text(0.5, 0.5, sprintf("Only %d features in common\n(need at least 10)",
                              length(common_features)), cex = 1.2)
      return()
    }

    # Calculate correlations
    correlations <- sapply(colnames(view_weights), function(f) {
      cor(de_res[common_features, "logFC"],
          view_weights[common_features, f],
          use = "complete.obs")
    })

    df <- data.frame(
      Factor = names(correlations),
      Correlation = correlations,
      stringsAsFactors = FALSE
    )
    df$Factor <- factor(df$Factor,
                        levels = df$Factor[order(abs(df$Correlation), decreasing = TRUE)])

    ggplot(df, aes(x = Factor, y = Correlation, fill = Correlation > 0)) +
      geom_col() +
      scale_fill_manual(values = c("TRUE" = "#E63946", "FALSE" = "#457B9D"), guide = "none") +
      geom_hline(yintercept = 0) +
      coord_flip() +
      labs(
        title = sprintf("Factor-DE Correlation: %s", contrast),
        subtitle = "Which MOFA factors align with your differential expression results?",
        x = NULL,
        y = "Pearson correlation (logFC vs factor weights)"
      ) +
      theme_minimal(base_size = 12)
  }, height = 500)

  # ============================================================================
  #      Template Download
  # ============================================================================

  output$download_mofa_template <- downloadHandler(
    filename = function() {
      paste0("mofa_view_template_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$raw_data)
      sample_names <- colnames(values$raw_data)
      template <- data.frame(
        FeatureID = c("GENE1_or_PROTEIN1", "GENE2_or_PROTEIN2", "GENE3_or_PROTEIN3"),
        matrix("", nrow = 3, ncol = length(sample_names),
               dimnames = list(NULL, sample_names)),
        check.names = FALSE
      )
      write.csv(template, file, row.names = FALSE)
    }
  )

  # ============================================================================
  #      Info Modal
  # ============================================================================

  observeEvent(input$mofa_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("layer-group"), " About Multi-Omics Factor Analysis (MOFA)"),
      size = "l", easyClose = TRUE,

      tags$h5("What does MOFA do?"),
      tags$p(
        "MOFA is like PCA, but for multiple datasets measured on the same samples. ",
        "It finds 'factors' -- hidden patterns that explain variation in your data."
      ),

      tags$h5("Why use MOFA with proteomics?"),
      tags$ul(
        tags$li(tags$b("Global + Phospho:"), " Separate abundance effects from true phospho-regulation"),
        tags$li(tags$b("Two experiments:"), " Find what's shared vs unique between measurements"),
        tags$li(tags$b("QC:"), " Discover hidden batch effects or sample outliers"),
        tags$li(tags$b("Multi-omics:"), " Integrate proteomics with RNA-seq, metabolomics, etc.")
      ),

      tags$h5("How to interpret results?"),
      tags$ul(
        tags$li(tags$b("Variance heatmap:"), " Shows which views each factor explains. ",
                "A factor explaining 30% of phospho but 2% of global = phospho-specific signal."),
        tags$li(tags$b("Factor weights:"), " The proteins/sites driving each factor. ",
                "High |weight| = strong contributor."),
        tags$li(tags$b("Sample scores:"), " Where samples fall in factor space. ",
                "Clustering = shared biology.")
      ),

      tags$h5("Key terminology"),
      tags$ul(
        tags$li(tags$b("View:"), " A data matrix (e.g., 'proteomics', 'phosphoproteomics')"),
        tags$li(tags$b("Factor:"), " A latent variable capturing coordinated variation"),
        tags$li(tags$b("Weight:"), " How much each feature contributes to a factor"),
        tags$li(tags$b("Score:"), " Each sample's position along a factor")
      ),

      tags$hr(),
      tags$h5("Example Datasets"),

      tags$h6("Mouse Brain (2-view)"),
      tags$ul(
        tags$li(tags$b("Source:"), " Bhattacharya et al. ",
                tags$em("Open-source multi-omics dataset of mouse somatosensory cortex")),
        tags$li(tags$b("GitHub:"), " ",
                tags$a(href = "https://github.com/gggrecco/Open-Source-Multi-Omics-Dataset-of-the-Somatosensory-Cortex-Mice",
                       target = "_blank", "gggrecco/Open-Source-Multi-Omics-Dataset")),
        tags$li(tags$b("Views:"), " Global proteomics (10,333 proteins) + Phosphoproteomics (89 sites)"),
        tags$li(tags$b("Design:"), " 16 samples, prenatal methadone exposure study (PME vs PSE, Male vs Female)"),
        tags$li(tags$b("Best for:"), " Testing abundance vs phospho-regulation deconvolution")
      ),

      tags$h6("TCGA Breast Cancer (3-view)"),
      tags$ul(
        tags$li(tags$b("Source:"), " The Cancer Genome Atlas (TCGA) via mixOmics R package (breast.TCGA dataset)"),
        tags$li(tags$b("Publication:"), " ",
                tags$a(href = "https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02015-1",
                       target = "_blank", "Rohart et al. (2017) PLOS Comp. Bio.")),
        tags$li(tags$b("Views:"), " mRNA (~200 genes) + miRNA (184 features) + Protein (142 features)"),
        tags$li(tags$b("Design:"), " 150 samples, 3 breast cancer subtypes (Basal, Her2, LumA)"),
        tags$li(tags$b("Best for:"), " Testing 3-view multi-omics integration")
      ),

      tags$hr(),
      tags$p(class = "text-muted",
        "MOFA reference: Argelaguet et al. (2020) Genome Biology. ",
        tags$a(href = "https://doi.org/10.1186/s13059-020-02015-1",
               target = "_blank", "DOI:10.1186/s13059-020-02015-1")
      ),

      footer = modalButton("Close")
    ))
  })

  # Parameters info modal
  observeEvent(input$mofa_params_info, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " MOFA Training Parameters"),
      size = "l", easyClose = TRUE,

      tags$h5("Number of Factors"),
      tags$p("How many latent factors to discover. Auto-select lets MOFA determine ",
             "the optimal number based on variance explained. Typical range: 5-15."),

      tags$h5("Convergence Mode"),
      tags$ul(
        tags$li(tags$b("Fast:"), " ~500 iterations, good for quick exploration"),
        tags$li(tags$b("Medium:"), " ~1000 iterations, recommended for final results"),
        tags$li(tags$b("Thorough:"), " ~5000 iterations, for publication-quality analysis")
      ),

      tags$h5("Advanced Options"),
      tags$ul(
        tags$li(tags$b("Scale views:"), " Equalizes the contribution of each view. ",
                "Recommended ON when views have very different numbers of features."),
        tags$li(tags$b("Center features:"), " Subtracts the mean of each feature. ",
                "Should almost always be ON."),
        tags$li(tags$b("Min variance:"), " Drops factors explaining less than this threshold. ",
                "Removes noise factors.")
      ),

      footer = modalButton("Close")
    ))
  })

  # ============================================================================
  #      Download variance explained
  # ============================================================================

  output$download_mofa_variance <- downloadHandler(
    filename = function() {
      paste0("mofa_variance_explained_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(values$mofa_variance_explained)
      r2 <- values$mofa_variance_explained$r2_per_factor
      var_list <- list()
      for (group_name in names(r2)) {
        mat <- r2[[group_name]]
        for (i in seq_len(nrow(mat))) {
          for (j in seq_len(ncol(mat))) {
            var_list[[length(var_list) + 1]] <- data.frame(
              View = colnames(mat)[j],
              Factor = rownames(mat)[i],
              Variance = mat[i, j],
              stringsAsFactors = FALSE
            )
          }
        }
      }
      combined <- do.call(rbind, var_list)
      write.csv(combined, file, row.names = FALSE)
    }
  )

  # ============================================================================
  #      Load MOFA Test Data (Mouse Brain Proteomics + Phosphoproteomics)
  # ============================================================================

  observeEvent(input$load_mofa_example, {
    withProgress(message = "Loading MOFA test data...", {

      base_url <- "https://github.com/bsphinney/DE-LIMP/releases/download/v1.0"
      prot_url   <- paste0(base_url, "/proteomics.csv")
      phospho_url <- paste0(base_url, "/phospho.csv")
      meta_url   <- paste0(base_url, "/metadata.csv")

      tryCatch({
        # Download files
        incProgress(0.1, detail = "Downloading proteomics...")
        prot_tmp <- tempfile(fileext = ".csv")
        curl::curl_download(prot_url, prot_tmp, quiet = TRUE)

        incProgress(0.3, detail = "Downloading phosphoproteomics...")
        phospho_tmp <- tempfile(fileext = ".csv")
        curl::curl_download(phospho_url, phospho_tmp, quiet = TRUE)

        incProgress(0.5, detail = "Downloading metadata...")
        meta_tmp <- tempfile(fileext = ".csv")
        curl::curl_download(meta_url, meta_tmp, quiet = TRUE)

        # Parse matrices
        incProgress(0.6, detail = "Loading matrices...")
        prot_df <- read.csv(prot_tmp, row.names = 1, check.names = FALSE)
        prot_mat <- as.matrix(prot_df)

        phospho_df <- read.csv(phospho_tmp, row.names = 1, check.names = FALSE)
        phospho_mat <- as.matrix(phospho_df)

        meta_df <- read.csv(meta_tmp, stringsAsFactors = FALSE)

        # Build group labels from metadata (Sex_Treatment)
        if (all(c("Sex", "Treatment") %in% names(meta_df))) {
          meta_df$Group <- paste(meta_df$Sex, meta_df$Treatment, sep = "_")
        }

        incProgress(0.8, detail = "Configuring views...")

        # Store as MOFA views
        values$mofa_views[["Global Proteomics"]] <- prot_mat
        values$mofa_views[["Phosphoproteomics"]] <- phospho_mat

        # Store metadata for sample score coloring
        values$mofa_sample_metadata <- meta_df

        # Build view configs
        values$mofa_view_configs <- list(
          list(
            id = "v1",
            num = 1,
            name = "Global Proteomics",
            type = "proteomics_other",
            source = "example",
            matrix = prot_mat,
            fit = NULL,
            n_features = nrow(prot_mat),
            n_samples = ncol(prot_mat),
            status = "ready"
          ),
          list(
            id = "v2",
            num = 2,
            name = "Phosphoproteomics",
            type = "phospho",
            source = "example",
            matrix = phospho_mat,
            fit = NULL,
            n_features = nrow(phospho_mat),
            n_samples = ncol(phospho_mat),
            status = "ready"
          )
        )

        # Clear any previous MOFA results
        values$mofa_object <- NULL
        values$mofa_factors <- NULL
        values$mofa_weights <- list()
        values$mofa_variance_explained <- NULL
        values$mofa_last_run_params <- NULL
        values$mofa_view_fits <- list()

        showNotification(
          sprintf(paste0(
            "MOFA test data loaded!\n",
            "View 1: Global Proteomics (%d proteins x %d samples)\n",
            "View 2: Phosphoproteomics (%d sites x %d samples)\n",
            "Design: %s"),
            nrow(prot_mat), ncol(prot_mat),
            nrow(phospho_mat), ncol(phospho_mat),
            paste(names(table(meta_df$Group)), collapse = ", ")),
          type = "message", duration = 8
        )

        add_to_log("Load MOFA Test Data", c(
          "# Mouse brain proteomics + phosphoproteomics test dataset",
          sprintf("# View 1: Global Proteomics (%d x %d)", nrow(prot_mat), ncol(prot_mat)),
          sprintf("# View 2: Phosphoproteomics (%d x %d)", nrow(phospho_mat), ncol(phospho_mat)),
          sprintf("# Groups: %s", paste(names(table(meta_df$Group)), collapse = ", "))
        ))

      }, error = function(e) {
        showNotification(paste("Error loading test data:", e$message),
                         type = "error", duration = 10)
      })
    })
  })

  # ============================================================================
  #      Load TCGA Breast Cancer Multi-Omics (3-view: mRNA + miRNA + Protein)
  # ============================================================================

  observeEvent(input$load_mofa_tcga, {
    withProgress(message = "Loading TCGA breast cancer data...", {

      base_url <- "https://github.com/bsphinney/DE-LIMP/releases/download/v1.0"
      mrna_url    <- paste0(base_url, "/tcga_mRNA.csv")
      mirna_url   <- paste0(base_url, "/tcga_miRNA.csv")
      protein_url <- paste0(base_url, "/tcga_protein.csv")
      meta_url    <- paste0(base_url, "/tcga_metadata.csv")

      tryCatch({
        # Download files
        incProgress(0.1, detail = "Downloading mRNA...")
        mrna_tmp <- tempfile(fileext = ".csv")
        curl::curl_download(mrna_url, mrna_tmp, quiet = TRUE)

        incProgress(0.25, detail = "Downloading miRNA...")
        mirna_tmp <- tempfile(fileext = ".csv")
        curl::curl_download(mirna_url, mirna_tmp, quiet = TRUE)

        incProgress(0.4, detail = "Downloading protein...")
        protein_tmp <- tempfile(fileext = ".csv")
        curl::curl_download(protein_url, protein_tmp, quiet = TRUE)

        incProgress(0.55, detail = "Downloading metadata...")
        meta_tmp <- tempfile(fileext = ".csv")
        curl::curl_download(meta_url, meta_tmp, quiet = TRUE)

        # Parse matrices
        incProgress(0.7, detail = "Loading matrices...")
        mrna_df <- read.csv(mrna_tmp, row.names = 1, check.names = FALSE)
        mrna_mat <- as.matrix(mrna_df)

        mirna_df <- read.csv(mirna_tmp, row.names = 1, check.names = FALSE)
        mirna_mat <- as.matrix(mirna_df)

        protein_df <- read.csv(protein_tmp, row.names = 1, check.names = FALSE)
        protein_mat <- as.matrix(protein_df)

        meta_df <- read.csv(meta_tmp, stringsAsFactors = FALSE)

        # Use Subtype as Group
        if ("Subtype" %in% names(meta_df)) {
          meta_df$Group <- meta_df$Subtype
        }

        incProgress(0.85, detail = "Configuring views...")

        # Store as MOFA views
        values$mofa_views[["mRNA"]] <- mrna_mat
        values$mofa_views[["miRNA"]] <- mirna_mat
        values$mofa_views[["Protein"]] <- protein_mat

        # Store metadata for sample score coloring
        values$mofa_sample_metadata <- meta_df

        # Build view configs
        values$mofa_view_configs <- list(
          list(
            id = "v1", num = 1, name = "mRNA",
            type = "proteomics_other", source = "example",
            matrix = mrna_mat, fit = NULL,
            n_features = nrow(mrna_mat), n_samples = ncol(mrna_mat),
            status = "ready"
          ),
          list(
            id = "v2", num = 2, name = "miRNA",
            type = "proteomics_other", source = "example",
            matrix = mirna_mat, fit = NULL,
            n_features = nrow(mirna_mat), n_samples = ncol(mirna_mat),
            status = "ready"
          ),
          list(
            id = "v3", num = 3, name = "Protein",
            type = "proteomics_other", source = "example",
            matrix = protein_mat, fit = NULL,
            n_features = nrow(protein_mat), n_samples = ncol(protein_mat),
            status = "ready"
          )
        )

        # Clear any previous MOFA results
        values$mofa_object <- NULL
        values$mofa_factors <- NULL
        values$mofa_weights <- list()
        values$mofa_variance_explained <- NULL
        values$mofa_last_run_params <- NULL
        values$mofa_view_fits <- list()

        showNotification(
          sprintf(paste0(
            "TCGA breast cancer data loaded!\n",
            "View 1: mRNA (%d genes x %d samples)\n",
            "View 2: miRNA (%d features x %d samples)\n",
            "View 3: Protein (%d features x %d samples)\n",
            "Subtypes: %s"),
            nrow(mrna_mat), ncol(mrna_mat),
            nrow(mirna_mat), ncol(mirna_mat),
            nrow(protein_mat), ncol(protein_mat),
            paste(names(table(meta_df$Group)), collapse = ", ")),
          type = "message", duration = 8
        )

        add_to_log("Load TCGA MOFA Test Data", c(
          "# TCGA breast cancer multi-omics (mRNA + miRNA + Protein)",
          sprintf("# View 1: mRNA (%d x %d)", nrow(mrna_mat), ncol(mrna_mat)),
          sprintf("# View 2: miRNA (%d x %d)", nrow(mirna_mat), ncol(mirna_mat)),
          sprintf("# View 3: Protein (%d x %d)", nrow(protein_mat), ncol(protein_mat)),
          sprintf("# Subtypes: %s", paste(names(table(meta_df$Group)), collapse = ", "))
        ))

      }, error = function(e) {
        showNotification(paste("Error loading TCGA data:", e$message),
                         type = "error", duration = 10)
      })
    })
  })
}
