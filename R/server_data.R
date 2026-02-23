# ==============================================================================
#  SERVER MODULE — Data Loading, Pipeline, Metadata, Contrast Sync
#  Called from app.R as: server_data(input, output, session, values, add_to_log, is_hf_space)
# ==============================================================================

server_data <- function(input, output, session, values, add_to_log, is_hf_space) {

  # ============================================================================
  #      2. Main Data Loading & Processing Pipeline
  # ============================================================================

  # Load example data from GitHub releases
  observeEvent(input$load_example, {
    withProgress(message = "Downloading example data...", {
      example_url <- "https://github.com/bsphinney/DE-LIMP/releases/download/v1.0/Affinisep_vs_evosep_noNorm.parquet"
      temp_file <- tempfile(fileext = ".parquet")

      tryCatch({
        incProgress(0.3, detail = "Downloading from GitHub...")
        download.file(example_url, temp_file, mode = "wb", quiet = TRUE)

        incProgress(0.5, detail = "Calculating Trends...")
        values$qc_stats <- get_diann_stats_r(temp_file)

        incProgress(0.7, detail = "Reading Matrix...")
        values$raw_data <- limpa::readDIANN(temp_file, format="parquet", q.cutoffs=input$q_cutoff)
        fnames <- sort(colnames(values$raw_data$E))
        values$metadata <- data.frame(
          ID = 1:length(fnames),
          File.Name = fnames,
          Group = rep("", length(fnames)),
          Batch = rep("", length(fnames)),
          Covariate1 = rep("", length(fnames)),
          Covariate2 = rep("", length(fnames)),
          stringsAsFactors=FALSE
        )

        # Initialize custom covariate names
        if(is.null(values$cov1_name)) values$cov1_name <- "Covariate1"
        if(is.null(values$cov2_name)) values$cov2_name <- "Covariate2"

        # Flag this as example data for auto-guess logic
        values$is_example_data <- TRUE

        # Detect DIA-NN normalization status
        values$diann_norm_detected <- tryCatch({
          raw_parquet <- arrow::read_parquet(temp_file,
            col_select = c("Precursor.Quantity", "Precursor.Normalised"))
          has_both_cols <- all(c("Precursor.Quantity", "Precursor.Normalised") %in% names(raw_parquet))
          if (has_both_cols) {
            sample_rows <- head(raw_parquet, 1000)
            ratio <- sample_rows$Precursor.Normalised / sample_rows$Precursor.Quantity
            ratios_vary <- sd(ratio, na.rm = TRUE) > 0.001
            if (ratios_vary) "on" else "off"
          } else { "unknown" }
        }, error = function(e) "unknown")

        # Store report path for XIC viewer precursor mapping (copy to session dir)
        session_report <- file.path(tempdir(), "de_limp_report.parquet")
        file.copy(temp_file, session_report, overwrite = TRUE)
        values$uploaded_report_path <- session_report
        values$original_report_name <- "Affinisep_vs_evosep_noNorm.parquet"

        # Auto-detect XIC directory in working directory for example data (local/HPC only)
        if (!is_hf_space) {
          tryCatch({
            cand <- file.path(getwd(), "Affinisep_vs_evosep_noNorm_xic")
            if (dir.exists(cand) && length(list.files(cand, pattern = "\\.xic\\.parquet$")) > 0) {
              updateTextInput(session, "xic_dir_input", value = cand)
              # Auto-load XICs after a short delay (let updateTextInput propagate)
              shinyjs::delay(500, shinyjs::click("xic_load_dir"))
            }
          }, error = function(e) NULL)
        }

        # Auto-detect phospho data
        values$phospho_detected <- detect_phospho(session_report)

        incProgress(0.9, detail = "Opening setup...")

        # Log to reproducibility
        add_to_log("Example Data Loaded", c(
          "# Example data: Affinisep vs Evosep (50ng Thermo Hela digest)",
          sprintf("# Downloaded from: %s", example_url),
          sprintf("dat <- readDIANN('Affinisep_vs_evosep_noNorm.parquet', format='parquet', q.cutoffs=%s)", input$q_cutoff)
        ))

        showNotification("Example data loaded successfully!", type = "message", duration = 3)
        # Navigate to Assign Groups sub-tab
        nav_select("main_tabs", "Data Overview")
        nav_select("data_overview_tabs", "Assign Groups & Run")

      }, error = function(e) {
        showNotification(paste("Error loading example data:", e$message), type = "error", duration = 10)
      })
    })
  })

  observeEvent(input$report_file, {
    req(input$report_file)
    withProgress(message = "Loading...", {
      file_mb <- round(file.size(input$report_file$datapath) / 1e6, 1)
      message(sprintf("[DE-LIMP] Loading report: %s (%.1f MB)", input$report_file$name, file_mb))

      incProgress(0.15, detail = sprintf("Calculating QC Trends (%.0f MB file)...", file_mb))
      values$qc_stats <- get_diann_stats_r(input$report_file$datapath)
      gc(verbose = FALSE)  # free QC stats intermediate memory

      incProgress(0.4, detail = "Reading expression matrix (this may take a while for large files)...")
      message(sprintf("[DE-LIMP] Memory before readDIANN: %.0f MB used", sum(gc()[,2])))
      tryCatch({
        values$raw_data <- limpa::readDIANN(input$report_file$datapath, format="parquet", q.cutoffs=input$q_cutoff)
        gc(verbose = FALSE)  # free readDIANN intermediates
        message(sprintf("[DE-LIMP] Memory after readDIANN: %.0f MB used", sum(gc()[,2])))
        fnames <- sort(colnames(values$raw_data$E))
        values$metadata <- data.frame(
          ID = 1:length(fnames),
          File.Name = fnames,
          Group = rep("", length(fnames)),
          Batch = rep("", length(fnames)),
          Covariate1 = rep("", length(fnames)),
          Covariate2 = rep("", length(fnames)),
          stringsAsFactors=FALSE
        )
        # Initialize custom covariate names (user can change these)
        if(is.null(values$cov1_name)) values$cov1_name <- "Covariate1"
        if(is.null(values$cov2_name)) values$cov2_name <- "Covariate2"
        # Clear example data flag for user uploads
        values$is_example_data <- FALSE

        # Detect DIA-NN normalization status
        values$diann_norm_detected <- tryCatch({
          raw_parquet <- arrow::read_parquet(input$report_file$datapath,
            col_select = c("Precursor.Quantity", "Precursor.Normalised"))
          has_both_cols <- all(c("Precursor.Quantity", "Precursor.Normalised") %in% names(raw_parquet))
          if (has_both_cols) {
            sample_rows <- head(raw_parquet, 1000)
            ratio <- sample_rows$Precursor.Normalised / sample_rows$Precursor.Quantity
            ratios_vary <- sd(ratio, na.rm = TRUE) > 0.001
            if (ratios_vary) "on" else "off"
          } else { "unknown" }
        }, error = function(e) "unknown")

        # Store report path for XIC viewer precursor mapping
        # For large files (>1 GB), avoid copying — Shiny keeps the upload alive for the session.
        # For smaller files, copy to a stable path (upload tempfiles can have unpredictable names).
        if (file_mb > 1000) {
          values$uploaded_report_path <- input$report_file$datapath
          message("[DE-LIMP] Large file — using upload path directly (skipping copy)")
        } else {
          session_report <- file.path(tempdir(), "de_limp_report.parquet")
          file.copy(input$report_file$datapath, session_report, overwrite = TRUE)
          values$uploaded_report_path <- session_report
        }
        values$original_report_name <- input$report_file$name

        # Auto-detect phospho data
        values$phospho_detected <- detect_phospho(values$uploaded_report_path)

        # Auto-detect XIC directory next to the uploaded report (local/HPC only)
        if (!is_hf_space) {
          tryCatch({
            report_name <- tools::file_path_sans_ext(input$report_file$name)
            # Check common locations: working directory, or if user uploaded from a known path
            candidate_dirs <- c(
              file.path(getwd(), paste0(report_name, "_xic")),
              file.path(dirname(input$report_file$datapath), paste0(report_name, "_xic"))
            )
            for (cand in candidate_dirs) {
              if (dir.exists(cand) && length(list.files(cand, pattern = "\\.xic\\.parquet$")) > 0) {
                updateTextInput(session, "xic_dir_input", value = cand)
                # Auto-load XICs after a short delay (let updateTextInput propagate)
                shinyjs::delay(500, shinyjs::click("xic_load_dir"))
                break
              }
            }
          }, error = function(e) NULL)
        }

        # Navigate to Assign Groups sub-tab
        nav_select("main_tabs", "Data Overview")
        nav_select("data_overview_tabs", "Assign Groups & Run")
      }, error=function(e) { showNotification(paste("Error:", e$message), type="error") })
    })
  })

  observeEvent(input$report_file, {
    req(input$report_file)
    add_to_log("Data Upload", c(
      sprintf("# File: %s", input$report_file$name),
      sprintf("dat <- readDIANN('%s', format='parquet', q.cutoffs=%s)",
              "path/to/your/report.parquet", input$q_cutoff)
    ))
  })

  # Old standalone "Run Pipeline" button observer - REMOVED
  # Pipeline now runs from the "Assign Groups & Run" sub-tab via run_pipeline observer

  # ============================================================================
  #      3. Metadata Handling (Assign Groups sub-tab)
  # ============================================================================

  output$hot_metadata <- renderRHandsontable({
    req(values$metadata)

    # Get custom covariate names
    cov1_display <- if(!is.null(input$cov1_label) && input$cov1_label != "") input$cov1_label else "Covariate1"
    cov2_display <- if(!is.null(input$cov2_label) && input$cov2_label != "") input$cov2_label else "Covariate2"

    # Store for later use
    values$cov1_name <- cov1_display
    values$cov2_name <- cov2_display

    # Create display dataframe with custom column names
    display_df <- values$metadata
    colnames(display_df) <- c("ID", "File.Name", "Group", "Batch", cov1_display, cov2_display)

    rhandsontable(display_df, rowHeaders=NULL, stretchH="all", height=500, width="100%") %>%
      hot_col("ID", readOnly=TRUE, width=50) %>%
      hot_col("File.Name", readOnly=TRUE) %>%
      hot_col("Group", type="text") %>%
      hot_col("Batch", type="text", width=100) %>%
      hot_col(cov1_display, type="text", width=100) %>%
      hot_col(cov2_display, type="text", width=100)
  })

  observeEvent(input$guess_groups, {
    req(values$metadata)
    meta <- if(!is.null(input$hot_metadata)) hot_to_r(input$hot_metadata) else values$metadata
    n <- nrow(meta)
    if (n < 2) return()

    # Use basenames, strip common extensions and leading date stamps
    fnames <- basename(meta$File.Name)
    fnames <- sub("\\.(d|raw|mzML|parquet)$", "", fnames, ignore.case = TRUE)
    fnames <- sub("^\\d{6,10}_", "", fnames)

    # --- Special case: example data (filenames don't distinguish groups) ---
    if (isTRUE(values$is_example_data) && n == 6) {
      meta$Group <- c(rep("Affinisep", 3), rep("Evosep", 3))
      values$metadata <- meta
      return()
    }

    # --- Strategy 1: Try known keywords first ---
    keywords <- c("affinisepACN", "affinisepIPA", "Control", "Treatment",
                   "Evosep", "Affinisep", "EGF", "untreat", "untreated",
                   "treated", "KO", "WT", "wildtype", "mutant", "vehicle",
                   "drug", "stim", "unstim", "inhibitor", "DMSO")

    find_keyword_match <- function(fname) {
      matches <- keywords[stringr::str_detect(fname, stringr::regex(keywords, ignore_case = TRUE))]
      if (length(matches) == 0) return("")
      matches[which.max(nchar(matches))]
    }

    guessed <- vapply(fnames, find_keyword_match, character(1), USE.NAMES = FALSE)

    # Check if keywords produced at least 2 groups
    keyword_groups <- unique(guessed[guessed != ""])
    if (length(keyword_groups) >= 2) {
      # Keywords worked — assign remaining as Sample_X
      sample_counter <- 0
      for (i in which(guessed == "")) {
        sample_counter <- sample_counter + 1
        guessed[i] <- paste0("Sample_", sample_counter)
      }
      meta$Group <- guessed
      values$metadata <- meta
      return()
    }

    # --- Strategy 2: Token-based auto-detection ---
    # Split filenames into tokens, find tokens that partition samples into groups
    tokens_per_file <- strsplit(fnames, "[_\\-\\.]+")

    # Collect all unique tokens (excluding pure numbers and very short ones)
    all_tokens <- unique(unlist(tokens_per_file))
    all_tokens <- all_tokens[nchar(all_tokens) >= 3]
    all_tokens <- all_tokens[!grepl("^[0-9]+$", all_tokens)]

    # For each token, check how many files contain it
    token_presence <- vapply(all_tokens, function(tok) {
      sum(vapply(tokens_per_file, function(toks) tok %in% toks, logical(1)))
    }, integer(1))

    # Good discriminating tokens appear in SOME but not ALL files (2+ groups)
    discriminating <- all_tokens[token_presence > 0 & token_presence < n]

    if (length(discriminating) > 0) {
      # Score tokens: prefer those that create balanced groups (close to n/2)
      token_scores <- vapply(discriminating, function(tok) {
        count <- token_presence[tok]
        # Penalize tokens in too many or too few files; prefer near n/2
        balance <- 1 - abs(count / n - 0.5) * 2
        # Prefer longer tokens (more specific)
        specificity <- min(nchar(tok) / 10, 1)
        balance * 0.7 + specificity * 0.3
      }, numeric(1))

      best_token <- discriminating[which.max(token_scores)]

      # Assign groups based on presence/absence of best token
      has_token <- vapply(tokens_per_file, function(toks) best_token %in% toks, logical(1))

      # Try to find a second discriminating token for the "other" group
      other_indices <- which(!has_token)
      other_tokens <- unique(unlist(tokens_per_file[other_indices]))
      other_tokens <- setdiff(other_tokens, unlist(tokens_per_file[has_token]))
      other_tokens <- other_tokens[nchar(other_tokens) >= 3 & !grepl("^[0-9]+$", other_tokens)]

      other_label <- if (length(other_tokens) > 0) {
        # Pick the most common non-shared token among the "other" group
        other_counts <- vapply(other_tokens, function(tok) {
          sum(vapply(tokens_per_file[other_indices], function(toks) tok %in% toks, logical(1)))
        }, integer(1))
        other_tokens[which.max(other_counts)]
      } else {
        paste0("non_", best_token)
      }

      guessed <- ifelse(has_token, best_token, other_label)
      meta$Group <- guessed
      values$metadata <- meta
      return()
    }

    # --- Fallback: number all samples ---
    meta$Group <- paste0("Sample_", seq_len(n))
    values$metadata <- meta
  })

  # Run pipeline from "Assign Groups & Run" sub-tab - saves groups and runs pipeline
  observeEvent(input$run_pipeline, {
    req(input$hot_metadata, values$metadata, values$raw_data)

    # First, save the groups
    old_meta <- values$metadata
    new_meta <- hot_to_r(input$hot_metadata)
    changed_indices <- which(old_meta$Group != new_meta$Group)
    if (length(changed_indices) > 0) {
      code_lines <- sprintf("metadata$Group[%d] <- '%s'  # %s",
                           changed_indices,
                           new_meta$Group[changed_indices],
                           new_meta$File.Name[changed_indices])
      add_to_log("Manual Group Assignment", code_lines)
    }
    values$metadata <- new_meta

    # Validate groups
    meta <- values$metadata
    meta$Group <- trimws(meta$Group)
    if(length(unique(meta$Group)) < 2) {
      showNotification("Error: Need at least 2 groups to run pipeline.", type="error", duration = 10)
      return()
    }

    showNotification("Groups saved! Running pipeline...", type="message")

    # Build covariates list for logging
    covariates_to_log <- character(0)
    cov_display_names <- character(0)

    if (isTRUE(input$include_batch) && length(unique(meta$Batch[meta$Batch != ""])) > 1) {
      covariates_to_log <- c(covariates_to_log, "Batch")
      cov_display_names <- c(cov_display_names, "Batch")
    }
    if (isTRUE(input$include_cov1) && length(unique(meta$Covariate1[meta$Covariate1 != ""])) > 1) {
      covariates_to_log <- c(covariates_to_log, "Covariate1")
      cov_display_names <- c(cov_display_names, values$cov1_name %||% "Covariate1")
    }
    if (isTRUE(input$include_cov2) && length(unique(meta$Covariate2[meta$Covariate2 != ""])) > 1) {
      covariates_to_log <- c(covariates_to_log, "Covariate2")
      cov_display_names <- c(cov_display_names, values$cov2_name %||% "Covariate2")
    }

    # Generate pipeline code for reproducibility log
    pipeline_code <- c(
      "# Normalization & Quantification",
      "dpcfit <- dpcCN(dat)",
      "y_protein <- dpcQuant(dat, 'Protein.Group', dpc=dpcfit)",
      ""
    )

    if (length(covariates_to_log) > 0) {
      # Log with covariates
      pipeline_code <- c(pipeline_code,
        sprintf("# Experimental Design (with covariates: %s)", paste(covariates_to_log, collapse = ", ")),
        "group_map <- c(",
        paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Group), collapse=",\n"),
        ")"
      )

      # Add covariate maps
      if ("Batch" %in% covariates_to_log) {
        pipeline_code <- c(pipeline_code,
          "batch_map <- c(",
          paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Batch), collapse=",\n"),
          ")"
        )
      }
      if ("Covariate1" %in% covariates_to_log) {
        cov1_name <- values$cov1_name %||% "Covariate1"
        pipeline_code <- c(pipeline_code,
          sprintf("%s_map <- c(", tolower(gsub(" ", "_", cov1_name))),
          paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Covariate1), collapse=",\n"),
          ")"
        )
      }
      if ("Covariate2" %in% covariates_to_log) {
        cov2_name <- values$cov2_name %||% "Covariate2"
        pipeline_code <- c(pipeline_code,
          sprintf("%s_map <- c(", tolower(gsub(" ", "_", cov2_name))),
          paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Covariate2), collapse=",\n"),
          ")"
        )
      }

      # Build metadata dataframe
      df_cols <- "Group = group_map"
      if ("Batch" %in% covariates_to_log) df_cols <- paste0(df_cols, ", Batch = batch_map")
      if ("Covariate1" %in% covariates_to_log) {
        cov1_name <- values$cov1_name %||% "Covariate1"
        cov1_var <- tolower(gsub(" ", "_", cov1_name))
        df_cols <- paste0(df_cols, sprintf(", %s = %s_map", cov1_name, cov1_var))
      }
      if ("Covariate2" %in% covariates_to_log) {
        cov2_name <- values$cov2_name %||% "Covariate2"
        cov2_var <- tolower(gsub(" ", "_", cov2_name))
        df_cols <- paste0(df_cols, sprintf(", %s = %s_map", cov2_name, cov2_var))
      }

      pipeline_code <- c(pipeline_code,
        sprintf("metadata <- data.frame(File.Name = names(group_map), %s)", df_cols),
        "metadata <- metadata[match(colnames(dat$E), metadata$File.Name), ]",
        "groups <- factor(metadata$Group)"
      )

      # Add factor creation for each covariate
      if ("Batch" %in% covariates_to_log) {
        pipeline_code <- c(pipeline_code, "batch <- factor(metadata$Batch)")
      }
      if ("Covariate1" %in% covariates_to_log) {
        cov1_name <- values$cov1_name %||% "Covariate1"
        cov1_var <- tolower(gsub(" ", "_", cov1_name))
        pipeline_code <- c(pipeline_code, sprintf("%s <- factor(metadata$%s)", cov1_var, cov1_name))
      }
      if ("Covariate2" %in% covariates_to_log) {
        cov2_name <- values$cov2_name %||% "Covariate2"
        cov2_var <- tolower(gsub(" ", "_", cov2_name))
        pipeline_code <- c(pipeline_code, sprintf("%s <- factor(metadata$%s)", cov2_var, cov2_name))
      }

      # Build design formula with custom names
      formula_parts <- c("groups")
      if ("Batch" %in% covariates_to_log) formula_parts <- c(formula_parts, "batch")
      if ("Covariate1" %in% covariates_to_log) formula_parts <- c(formula_parts, tolower(gsub(" ", "_", values$cov1_name %||% "covariate1")))
      if ("Covariate2" %in% covariates_to_log) formula_parts <- c(formula_parts, tolower(gsub(" ", "_", values$cov2_name %||% "covariate2")))
      formula_str <- paste0("~ 0 + ", paste(formula_parts, collapse = " + "))
      pipeline_code <- c(pipeline_code,
        sprintf("design <- model.matrix(%s)", formula_str),
        "colnames(design) <- gsub('groups', '', colnames(design))"
      )

      pipeline_code <- c(pipeline_code, "",
        "# Differential Expression Model (with covariates)",
        "fit <- dpcDE(y_protein, design, plot=FALSE)"
      )
    } else {
      # Log without covariates
      pipeline_code <- c(pipeline_code,
        "# Experimental Design",
        "group_map <- c(",
        paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Group), collapse=",\n"),
        ")",
        "metadata <- data.frame(File.Name = names(group_map), Group = group_map)",
        "metadata <- metadata[match(colnames(dat$E), metadata$File.Name), ]",
        "groups <- factor(metadata$Group)",
        "design <- model.matrix(~ 0 + groups)",
        "colnames(design) <- levels(groups)",
        "",
        "# Differential Expression Model",
        "fit <- dpcDE(y_protein, design, plot=FALSE)"
      )
    }

    add_to_log("Run Pipeline (Main Analysis)", pipeline_code)

    withProgress(message='Running Pipeline...', {
      tryCatch({
        dat <- values$raw_data
        message(sprintf("[DE-LIMP] Pipeline start — %d samples, memory: %.0f MB",
                        ncol(dat$E), sum(gc(verbose = FALSE)[,2])))

        incProgress(0.2, detail = "Normalizing (DPC-CN)...")
        dpcfit <- limpa::dpcCN(dat)
        values$dpc_fit <- dpcfit
        gc(verbose = FALSE)

        incProgress(0.5, detail = "Protein quantification (maxLFQ)...")
        values$y_protein <- tryCatch({
          result <- limpa::dpcQuant(dat, "Protein.Group", dpc=dpcfit)
          gc(verbose = FALSE)
          message(sprintf("[DE-LIMP] Quantification done — %d proteins, memory: %.0f MB",
                          nrow(result$E), sum(gc(verbose = FALSE)[,2])))
          result
        }, error = function(e) {
          showNotification(paste("Protein quantification failed:", e$message), type = "error", duration = NULL)
          return(NULL)
        })

        req(values$y_protein)

        rownames(meta) <- meta$File.Name
        meta <- meta[colnames(dat$E), ]
        meta$Group <- make.names(meta$Group)
        groups <- factor(meta$Group)

        # Build design formula with selected covariates
        covariates_to_include <- character(0)
        covariate_messages <- character(0)

        # Check each covariate
        if (isTRUE(input$include_batch) && length(unique(meta$Batch[meta$Batch != ""])) > 1) {
          meta$Batch <- make.names(meta$Batch)
          covariates_to_include <- c(covariates_to_include, "batch")
          covariate_messages <- c(covariate_messages, "Batch")
        }

        if (isTRUE(input$include_cov1) && length(unique(meta$Covariate1[meta$Covariate1 != ""])) > 1) {
          meta$Covariate1 <- make.names(meta$Covariate1)
          cov1_var <- tolower(gsub(" ", "_", values$cov1_name %||% "covariate1"))
          covariates_to_include <- c(covariates_to_include, cov1_var)
          covariate_messages <- c(covariate_messages, values$cov1_name %||% "Covariate1")
        }

        if (isTRUE(input$include_cov2) && length(unique(meta$Covariate2[meta$Covariate2 != ""])) > 1) {
          meta$Covariate2 <- make.names(meta$Covariate2)
          cov2_var <- tolower(gsub(" ", "_", values$cov2_name %||% "covariate2"))
          covariates_to_include <- c(covariates_to_include, cov2_var)
          covariate_messages <- c(covariate_messages, values$cov2_name %||% "Covariate2")
        }

        # Build design matrix
        if (length(covariates_to_include) > 0) {
          # Create formula dynamically
          formula_str <- paste0("~ 0 + groups + ", paste(covariates_to_include, collapse = " + "))

          # Build data frame for design matrix
          design_df <- data.frame(groups = groups)
          if("batch" %in% covariates_to_include) design_df$batch <- factor(meta$Batch)
          if(isTRUE(input$include_cov1) && length(unique(meta$Covariate1[meta$Covariate1 != ""])) > 1) {
            cov1_var <- tolower(gsub(" ", "_", values$cov1_name %||% "covariate1"))
            design_df[[cov1_var]] <- factor(meta$Covariate1)
          }
          if(isTRUE(input$include_cov2) && length(unique(meta$Covariate2[meta$Covariate2 != ""])) > 1) {
            cov2_var <- tolower(gsub(" ", "_", values$cov2_name %||% "covariate2"))
            design_df[[cov2_var]] <- factor(meta$Covariate2)
          }

          design <- model.matrix(as.formula(formula_str), data = design_df)
          colnames(design) <- gsub("groups", "", colnames(design))

          showNotification(
            paste0("Including covariates: ", paste(covariate_messages, collapse = ", ")),
            type = "message", duration = 5
          )
        } else {
          # Standard design without covariates
          design <- model.matrix(~ 0 + groups)
          colnames(design) <- levels(groups)
        }

        combs <- combn(levels(groups), 2)
        forms <- apply(combs, 2, function(x) paste(x[2], "-", x[1]))

        fit <- limpa::dpcDE(values$y_protein, design, plot=FALSE)
        fit <- contrasts.fit(fit, makeContrasts(contrasts=forms, levels=design))
        fit <- eBayes(fit)
        values$fit <- fit

        # Update all four comparison selectors
        updateSelectInput(session, "contrast_selector", choices=forms)
        updateSelectInput(session, "contrast_selector_signal", choices=forms, selected=forms[1])
        updateSelectInput(session, "contrast_selector_grid", choices=forms, selected=forms[1])
        updateSelectInput(session, "contrast_selector_pvalue", choices=forms, selected=forms[1])
        values$status <- "\u2705 Complete!"

        # Log contrasts
        contrast_code <- c(
          sprintf("# Available contrasts: %s", paste(forms, collapse=", ")),
          "combs <- combn(levels(groups), 2)",
          "forms <- apply(combs, 2, function(x) paste(x[2], '-', x[1]))",
          "fit <- contrasts.fit(fit, makeContrasts(contrasts=forms, levels=design))",
          "fit <- eBayes(fit)"
        )
        add_to_log("Contrast Fitting", contrast_code)

        nav_select("main_tabs", "QC")
        showNotification("\u2713 Pipeline complete! View results in tabs below.", type="message", duration=10)
      }, error = function(e) {
        showNotification(paste("Pipeline error:", e$message), type = "error", duration = NULL)
      })
    })
  })

  # ============================================================================
  #      4. Parameter Change Logging
  # ============================================================================

  # Log contrast changes
  observeEvent(input$contrast_selector, {
    req(input$contrast_selector)
    if (!is.null(values$fit)) {
      add_to_log("Select Contrast for Visualization", c(
        sprintf("# Viewing contrast: %s", input$contrast_selector),
        sprintf("results <- topTable(fit, coef='%s', number=Inf)", input$contrast_selector)
      ))
    }

    # Sync with Signal Distribution, Expression Grid, and P-value Distribution selectors
    if (!is.null(input$contrast_selector_signal) && input$contrast_selector_signal != input$contrast_selector) {
      updateSelectInput(session, "contrast_selector_signal", selected = input$contrast_selector)
    }
    if (!is.null(input$contrast_selector_grid) && input$contrast_selector_grid != input$contrast_selector) {
      updateSelectInput(session, "contrast_selector_grid", selected = input$contrast_selector)
    }
    if (!is.null(input$contrast_selector_pvalue) && input$contrast_selector_pvalue != input$contrast_selector) {
      updateSelectInput(session, "contrast_selector_pvalue", selected = input$contrast_selector)
    }
  })

  # Sync Signal Distribution selector with main selector
  observeEvent(input$contrast_selector_signal, {
    req(input$contrast_selector_signal)
    if (!is.null(input$contrast_selector) && input$contrast_selector != input$contrast_selector_signal) {
      updateSelectInput(session, "contrast_selector", selected = input$contrast_selector_signal)
    }
  })

  # Sync Expression Grid selector with main selector
  observeEvent(input$contrast_selector_grid, {
    req(input$contrast_selector_grid)
    if (!is.null(input$contrast_selector) && input$contrast_selector != input$contrast_selector_grid) {
      updateSelectInput(session, "contrast_selector", selected = input$contrast_selector_grid)
    }
  })

  # Sync P-value Distribution selector with main selector
  observeEvent(input$contrast_selector_pvalue, {
    req(input$contrast_selector_pvalue)
    if (!is.null(input$contrast_selector) && input$contrast_selector != input$contrast_selector_pvalue) {
      updateSelectInput(session, "contrast_selector", selected = input$contrast_selector_pvalue)
    }
  })

  # Log logFC threshold changes
  observeEvent(input$logfc_cutoff, {
    req(input$logfc_cutoff)
    if (!is.null(values$fit)) {
      add_to_log("Change LogFC Threshold", c(
        sprintf("# LogFC cutoff set to: %.2f", input$logfc_cutoff),
        sprintf("# Filter: adj.P.Val < 0.05 & abs(logFC) > %.2f", input$logfc_cutoff)
      ))
    }
  })

}
