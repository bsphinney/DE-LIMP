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
        # NOTE (v3.9.7): QuantUMS filtering moved to pipeline run-time inside
        # build_maxlfq_pipeline(). Load handlers always read the unfiltered
        # parquet so DPC-Quant gets paper-faithful input regardless of slider
        # values. The filter sliders only take effect when MaxLFQ + limma runs.
        values$raw_data <- limpa::readDIANN(temp_file, format="parquet", q.cutoffs=input$q_cutoff)
        values$quantums_filter_applied <- character(0)
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
        # NOTE (v3.9.7): QuantUMS filtering happens at pipeline run-time only
        # (build_maxlfq_pipeline). Load always reads the unfiltered parquet.
        values$raw_data <- limpa::readDIANN(input$report_file$datapath, format="parquet", q.cutoffs=input$q_cutoff)
        values$quantums_filter_applied <- character(0)
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
    message(sprintf("[DE-LIMP] Rendering metadata table: %d rows x %d cols",
                    nrow(values$metadata), ncol(values$metadata)))

    # Get custom covariate names (Batch column also renamable)
    batch_display <- if (!is.null(input$batch_label) && input$batch_label != "") input$batch_label else "Batch"
    cov1_display  <- if (!is.null(input$cov1_label)  && input$cov1_label  != "") input$cov1_label  else "Covariate1"
    cov2_display  <- if (!is.null(input$cov2_label)  && input$cov2_label  != "") input$cov2_label  else "Covariate2"

    # Store for later use
    values$batch_name <- batch_display
    values$cov1_name  <- cov1_display
    values$cov2_name  <- cov2_display

    # Append excluded files as read-only rows with "Excluded" status
    display_df <- values$metadata
    n_active <- nrow(display_df)
    excluded_rows <- integer(0)

    if (!is.null(values$excluded_files) && nrow(values$excluded_files) > 0) {
      ef <- values$excluded_files
      excl_df <- data.frame(
        ID = seq(n_active + 1, n_active + nrow(ef)),
        File.Name = ef$filename,
        Group = ifelse(nzchar(ef$group), ef$group, "[Excluded]"),
        Batch = "",
        Covariate1 = "",
        Covariate2 = "",
        stringsAsFactors = FALSE
      )
      display_df <- rbind(display_df, excl_df)
      excluded_rows <- seq(n_active + 1, nrow(display_df))
    }

    colnames(display_df) <- c("ID", "File.Name", "Group", batch_display, cov1_display, cov2_display)

    # Add Status column to mark excluded rows
    display_df$Status <- ""
    if (length(excluded_rows) > 0) {
      display_df$Status[excluded_rows] <- "Excluded"
    }

    # Custom renderer JS for red-background excluded rows
    excluded_renderer <- if (length(excluded_rows) > 0) {
      sprintf("function(instance, td, row, col, prop, value, cellProperties) {
        Handsontable.renderers.TextRenderer.apply(this, arguments);
        var excludedRows = [%s];
        if (excludedRows.indexOf(row) > -1) {
          td.style.backgroundColor = '#f8d7da';
          td.style.color = '#721c24';
          td.style.fontStyle = 'italic';
        }
      }", paste(excluded_rows - 1, collapse = ","))
    } else NULL

    hot <- rhandsontable(display_df, rowHeaders=NULL, stretchH="all", height=500, width="100%")

    if (!is.null(excluded_renderer)) {
      hot <- hot %>%
        hot_col("ID", readOnly=TRUE, width=50, renderer = excluded_renderer) %>%
        hot_col("File.Name", readOnly=TRUE, renderer = excluded_renderer) %>%
        hot_col("Group", type="text", renderer = excluded_renderer) %>%
        hot_col(batch_display, type="text", width=100, renderer = excluded_renderer) %>%
        hot_col(cov1_display, type="text", width=100, renderer = excluded_renderer) %>%
        hot_col(cov2_display, type="text", width=100, renderer = excluded_renderer) %>%
        hot_col("Status", readOnly=TRUE, width=70, renderer = excluded_renderer)
    } else {
      hot <- hot %>%
        hot_col("ID", readOnly=TRUE, width=50) %>%
        hot_col("File.Name", readOnly=TRUE) %>%
        hot_col("Group", type="text") %>%
        hot_col(batch_display, type="text", width=100) %>%
        hot_col(cov1_display, type="text", width=100) %>%
        hot_col(cov2_display, type="text", width=100) %>%
        hot_col("Status", readOnly=TRUE, width=70)
    }

    hot
  })

  # Helper: sync excluded file groups from the rhandsontable back to values$excluded_files
  sync_excluded_groups <- function() {
    if (is.null(values$excluded_files) || nrow(values$excluded_files) == 0) return()
    if (is.null(values$metadata) || is.null(input$hot_metadata)) return()

    tbl <- hot_to_r(input$hot_metadata)
    n_active <- nrow(values$metadata)
    n_total <- nrow(tbl)
    if (n_total <= n_active) return()

    excl_rows <- tbl[(n_active + 1):n_total, ]
    # Column 3 is Group regardless of display name
    for (i in seq_len(min(nrow(excl_rows), nrow(values$excluded_files)))) {
      grp <- as.character(excl_rows[i, 3])
      if (!identical(grp, "[Excluded]") && nzchar(grp)) {
        values$excluded_files$group[i] <- grp
      }
    }
  }

  observeEvent(input$guess_groups, {
    req(values$metadata)
    meta <- if(!is.null(input$hot_metadata)) {
      tbl <- hot_to_r(input$hot_metadata)
      colnames(tbl) <- c("ID", "File.Name", "Group", "Batch", "Covariate1", "Covariate2", "Status")
      # Only keep active rows (exclude excluded files) and drop Status column
      tbl[seq_len(nrow(values$metadata)), c("ID", "File.Name", "Group", "Batch", "Covariate1", "Covariate2")]
    } else values$metadata
    n <- nrow(meta)
    if (n < 2) return()

    # Use basenames, strip common extensions and leading date stamps
    fnames <- basename(meta$File.Name)
    fnames <- sub("\\.(d|raw|mzML|parquet)$", "", fnames, ignore.case = TRUE)
    fnames <- sub("^\\d{6,10}_", "", fnames)

    # --- Special case: example data (filenames don't auto-guess cleanly) ---
    if (isTRUE(values$is_example_data) && !isTRUE(values$is_example_phospho)) {
      meta$Group <- ifelse(grepl("affinisepIPA", meta$File.Name), "affinisepIPA",
                    ifelse(grepl("affinisepACN", meta$File.Name), "affinisepACN",
                    ifelse(grepl("affinisep", meta$File.Name, ignore.case = TRUE), "Affinisep",
                    "Evosep")))
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

    # First, save the groups — separate active rows from excluded rows
    old_meta <- values$metadata
    full_table <- hot_to_r(input$hot_metadata)
    # Table has Status column (7th) — name all columns including Status
    colnames(full_table) <- c("ID", "File.Name", "Group", "Batch", "Covariate1", "Covariate2", "Status")
    n_active <- nrow(old_meta)
    new_meta <- full_table[seq_len(n_active), c("ID", "File.Name", "Group", "Batch", "Covariate1", "Covariate2")]

    # Sync excluded file groups from table (rows after active)
    sync_excluded_groups()

    changed_indices <- which(old_meta$Group != new_meta$Group)
    if (length(changed_indices) > 0) {
      code_lines <- sprintf("metadata$Group[%d] <- '%s'  # %s",
                           changed_indices,
                           new_meta$Group[changed_indices],
                           new_meta$File.Name[changed_indices])
      add_to_log("Manual Group Assignment", code_lines)
    }
    values$metadata <- new_meta

    # Validate groups (active files only)
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

    # Generate pipeline code for reproducibility log.
    # Branch on the actual pipeline that's about to run so the log accurately
    # reflects what R did (v3.9.12).
    is_maxlfq_log <- (input$pipeline_mode %||% "dpc") == "maxlfq" &&
                      !isTRUE(input$use_limpa_with_filter)
    if (is_maxlfq_log) {
      qc <- input$q_cutoff %||% 0.01
      eqc <- input$eq_cutoff %||% 0
      pgqc <- input$pgq_cutoff %||% 0
      cov_pct <- 100 * (input$coverage_min_frac %||% 0.5)
      pipeline_code <- c(
        "# ─── MaxLFQ + limma pipeline (Moschem et al. 2025) ───",
        "# DOI: 10.1021/acs.jproteome.5c00009",
        "library(arrow); library(dplyr); library(tidyr); library(limma)",
        "",
        "# 1. Read parquet, apply FDR + QuantUMS filters lazily",
        "rep <- arrow::open_dataset('report.parquet') %>%",
        "  dplyr::select(Run, Protein.Group, PG.MaxLFQ,",
        "                Q.Value, Lib.Q.Value, Lib.PG.Q.Value,",
        "                Empirical.Quality, PG.MaxLFQ.Quality, Genes) %>%",
        sprintf("  dplyr::filter(Q.Value <= %.3f, Lib.Q.Value <= %.3f, Lib.PG.Q.Value <= %.3f%s%s) %%>%%",
                qc, qc, qc,
                if (eqc > 0) sprintf(",\n                Empirical.Quality >= %.2f", eqc) else "",
                if (pgqc > 0) sprintf(",\n                PG.MaxLFQ.Quality >= %.2f", pgqc) else ""),
        "  dplyr::group_by(Protein.Group, Run) %>%",
        "  dplyr::summarise(PG.MaxLFQ = max(PG.MaxLFQ, na.rm = TRUE), .groups = 'drop') %>%",
        "  dplyr::collect()",
        "",
        "# 2. Pivot wide → log2 → quantile-normalize",
        "wide <- rep %>% tidyr::pivot_wider(",
        "  id_cols = Protein.Group, names_from = Run, values_from = PG.MaxLFQ)",
        "E <- as.matrix(wide[, -1, drop = FALSE]); rownames(E) <- wide$Protein.Group",
        "E[E <= 0 | !is.finite(E)] <- NA_real_",
        "E <- log2(E)",
        "E <- limma::normalizeBetweenArrays(E, method = 'quantile')",
        "",
        sprintf("# 3. Coverage filter — drop proteins with < %.0f%% non-NA samples", cov_pct),
        sprintf("keep <- rowSums(!is.na(E)) >= ceiling(%.2f * ncol(E))", cov_pct / 100),
        "E_fit <- E[keep, , drop = FALSE]",
        ""
      )
    } else {
      qc_dpc <- input$q_cutoff %||% 0.01
      report_path_dpc <- values$uploaded_report_path %||% "report.parquet"
      pipeline_code <- c(
        "# ─── DPC-Quant pipeline (limpa, default DE-LIMP path) ───",
        "library(limpa); library(limma)",
        "",
        "# 1. Read DIA-NN report.parquet, applying identification FDR cutoffs",
        sprintf("dat <- limpa::readDIANN('%s', format = 'parquet', q.cutoffs = %.3f)",
                report_path_dpc, qc_dpc),
        "",
        "# 2. Normalize + quantify proteins via DPC-CN + DPC-Quant",
        "dpcfit <- limpa::dpcCN(dat)",
        "y_protein <- limpa::dpcQuant(dat, 'Protein.Group', dpc = dpcfit)",
        ""
      )
    }

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
        if (is_maxlfq_log) "fit <- limma::lmFit(E_fit, design)" else "fit <- dpcDE(y_protein, design, plot=FALSE)"
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
        if (is_maxlfq_log) "fit <- limma::lmFit(E_fit, design)" else "fit <- dpcDE(y_protein, design, plot=FALSE)"
      )
    }

    add_to_log("Run Pipeline (Main Analysis)", pipeline_code)

    withProgress(message='Running Pipeline...', {
      tryCatch({
        dat <- values$raw_data
        message(sprintf("[DE-LIMP] Pipeline start — %d samples, memory: %.0f MB",
                        ncol(dat$E), sum(gc(verbose = FALSE)[,2])))

        # v3.9 — choose between DPC-Quant (limpa) and MaxLFQ + limma (Moschem 2025)
        pipeline_mode <- input$pipeline_mode %||% "dpc"
        use_limpa_override <- isTRUE(input$use_limpa_with_filter)
        use_maxlfq <- (pipeline_mode == "maxlfq") && !use_limpa_override

        if (use_maxlfq) {
          # ---- Paper-faithful MaxLFQ + limma branch ----
          incProgress(0.5, detail = "Building MaxLFQ matrix (paper-faithful)...")
          parquet_path <- values$uploaded_report_path
          if (is.null(parquet_path) || !file.exists(parquet_path)) {
            showNotification("MaxLFQ pipeline needs the loaded report.parquet on disk. Reload the file and try again.",
                             type = "error", duration = NULL)
            return(invisible(NULL))
          }
          # Honour the user's excluded_files set: pass meta$File.Name as the
          # keep-list so the MaxLFQ matrix matches the metadata table exactly.
          keep_runs_maxlfq <- meta$File.Name
          values$y_protein <- tryCatch({
            res <- build_maxlfq_pipeline(parquet_path,
                     q_cutoff   = input$q_cutoff   %||% 0.01,
                     eq_cutoff  = input$eq_cutoff  %||% 0,
                     pgq_cutoff = input$pgq_cutoff %||% 0,
                     keep_runs  = keep_runs_maxlfq)
            gc(verbose = FALSE)
            # Per-filter precursor counts (visible feedback for the user)
            fc <- res$other$filter_counts
            if (!is.null(fc$input) && !is.na(fc$input)) {
              message(sprintf("[DE-LIMP] MaxLFQ filters — input precursor rows: %s",
                              format(fc$input, big.mark = ",")))
              if (!is.null(fc$after_fdr))
                message(sprintf("[DE-LIMP]   after FDR (Q ≤ %.3f): %s (%.1f%% kept)",
                                input$q_cutoff %||% 0.01,
                                format(fc$after_fdr, big.mark = ","),
                                100 * fc$after_fdr / fc$input))
              if (!is.null(fc$after_eq))
                message(sprintf("[DE-LIMP]   after eQ ≥ %.2f: %s (%.1f%% kept of FDR pool)",
                                input$eq_cutoff %||% 0,
                                format(fc$after_eq, big.mark = ","),
                                100 * fc$after_eq / max(fc$after_fdr %||% fc$input, 1)))
              if (!is.null(fc$after_pgq))
                message(sprintf("[DE-LIMP]   after pgQ ≥ %.2f: %s (%.1f%% kept of eQ pool)",
                                input$pgq_cutoff %||% 0,
                                format(fc$after_pgq, big.mark = ","),
                                100 * fc$after_pgq / max(fc$after_eq %||% fc$after_fdr %||% fc$input, 1)))
              if (!is.null(fc$after_excluded_files))
                message(sprintf("[DE-LIMP]   after excluded-runs filter: %s",
                                format(fc$after_excluded_files, big.mark = ",")))
            }
            message(sprintf("[DE-LIMP] MaxLFQ pipeline: %d proteins x %d runs, %d cells missing (%.1f%%).",
                            res$other$n_proteins_in_matrix, res$other$n_runs,
                            res$other$n_cells_missing,
                            100 * res$other$n_cells_missing / res$other$n_cells_total))
            res
          }, error = function(e) {
            showNotification(paste("MaxLFQ pipeline failed:", e$message),
                             type = "error", duration = NULL)
            return(NULL)
          })
          values$dpc_fit <- NULL
          values$pipeline_mode_used <- "maxlfq"
        } else {
          # ---- DPC-Quant (limpa) branch ----
          if (pipeline_mode == "maxlfq" && use_limpa_override) {
            showNotification(
              paste0("Experimental: running limpa DPC-Quant on QuantUMS-filtered precursors. ",
                     "This combination is not tested in either paper — DPC-Quant assumes no pre-filtering."),
              type = "warning", duration = 15)
            message("[DE-LIMP] Running experimental combo: QuantUMS filter + limpa DPC-Quant.")
          }
          incProgress(0.2, detail = "Normalizing (DPC-CN)...")
          dpcfit <- limpa::dpcCN(dat)
          values$dpc_fit <- dpcfit
          gc(verbose = FALSE)

          incProgress(0.5, detail = "Protein quantification (DPC-Quant)...")
          values$y_protein <- tryCatch({
            result <- limpa::dpcQuant(dat, "Protein.Group", dpc=dpcfit)
            gc(verbose = FALSE)
            message(sprintf("[DE-LIMP] Quantification done — %d proteins, memory: %.0f MB",
                            nrow(result$E), sum(gc(verbose = FALSE)[,2])))
            # Attach canonical pipeline descriptor (v3.9.17+) so downstream
            # code reads `is_maxlfq(y_protein)` / `pipeline_label(y_protein)`
            # instead of the volatile `values$pipeline_mode_used` reactiveVal.
            result$other$descriptor <- dpc_pipeline_descriptor()
            result$other$pipeline   <- "dpc"
            result
          }, error = function(e) {
            showNotification(paste("Protein quantification failed:", e$message), type = "error", duration = NULL)
            return(NULL)
          })
          values$pipeline_mode_used <- if (use_limpa_override) "dpc_with_filter_experimental" else "dpc"
        }

        req(values$y_protein)

        rownames(meta) <- meta$File.Name
        # Align metadata to the SAMPLE matrix actually produced by the chosen
        # pipeline. Under DPC-Quant the matrix is `dat$E`; under MaxLFQ it's
        # `values$y_protein$E`. limma::lmFit matches by column position, so
        # this row order MUST match the matrix column order exactly.
        sample_cols_used <- if (isTRUE(values$pipeline_mode_used == "maxlfq")) {
          colnames(values$y_protein$E)
        } else {
          colnames(dat$E)
        }
        missing_meta <- setdiff(sample_cols_used, meta$File.Name)
        if (length(missing_meta) > 0) {
          showNotification(paste0("Pipeline aborted: ", length(missing_meta),
            " sample(s) in the matrix have no metadata row (",
            paste(head(missing_meta, 5), collapse = ", "),
            if (length(missing_meta) > 5) ", ..." else "",
            "). This usually means the parquet contains runs that were excluded ",
            "in the metadata table — re-load the report after fixing exclusions."),
            type = "error", duration = NULL)
          return(invisible(NULL))
        }
        meta <- meta[sample_cols_used, , drop = FALSE]
        meta$Group <- make.names(meta$Group)
        groups <- factor(meta$Group)

        # Build design formula with selected covariates.
        # Each covariate is auto-coerced to numeric (continuous) when the
        # column looks like a numeric identifier (Run order, age, etc.) so
        # users can't accidentally turn a per-sample number into a 200-level
        # factor and blow up the design matrix rank.
        covariates_to_include <- character(0)
        covariate_messages <- character(0)
        design_df <- data.frame(groups = groups)
        warnings_out <- character(0)

        add_covariate <- function(slot_name, raw_values, label) {
          info <- coerce_covariate_column(raw_values)
          if (info$kind == "numeric") {
            design_df[[slot_name]] <<- info$values
            covariates_to_include <<- c(covariates_to_include, slot_name)
            covariate_messages <<- c(covariate_messages,
                                     sprintf("%s (numeric)", label))
            message(sprintf("[DE-LIMP] Covariate '%s' treated as numeric (%d distinct values).",
                            label, info$n_levels))
          } else {
            if (info$n_levels < 2) {
              warnings_out <<- c(warnings_out,
                sprintf("Covariate '%s' has fewer than 2 levels — skipped.", label))
              return(invisible(NULL))
            }
            if (info$has_singletons) {
              n_sing <- length(info$singleton_levels)
              warnings_out <<- c(warnings_out, sprintf(
                "Covariate '%s' has %d level(s) that occur in only one sample (%s) — these break the model. Either drop those rows or merge them into another level.",
                label, n_sing,
                if (n_sing > 5) paste0(paste(head(info$singleton_levels, 5), collapse = ", "), ", …")
                else paste(info$singleton_levels, collapse = ", ")))
              return(invisible(NULL))
            }
            design_df[[slot_name]] <<- info$values
            covariates_to_include <<- c(covariates_to_include, slot_name)
            covariate_messages <<- c(covariate_messages,
                                     sprintf("%s (factor, %d levels)", label, info$n_levels))
            message(sprintf("[DE-LIMP] Covariate '%s' treated as factor (%d levels).",
                            label, info$n_levels))
          }
        }

        if (isTRUE(input$include_batch)) {
          add_covariate("batch", meta$Batch, values$batch_name %||% "Batch")
        }
        if (isTRUE(input$include_cov1)) {
          slot1 <- tolower(gsub(" ", "_", values$cov1_name %||% "covariate1"))
          add_covariate(slot1, meta$Covariate1, values$cov1_name %||% "Covariate1")
        }
        if (isTRUE(input$include_cov2)) {
          slot2 <- tolower(gsub(" ", "_", values$cov2_name %||% "covariate2"))
          add_covariate(slot2, meta$Covariate2, values$cov2_name %||% "Covariate2")
        }

        if (length(warnings_out) > 0) {
          for (w in warnings_out) message("[DE-LIMP] Covariate warning: ", w)
          showNotification(paste(warnings_out, collapse = "  "),
                           type = "warning", duration = 15)
        }

        # Build design matrix
        if (length(covariates_to_include) > 0) {
          formula_str <- paste0("~ 0 + groups + ", paste(covariates_to_include, collapse = " + "))
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

        # Pre-flight: refuse a rank-deficient design before limma blows up
        # with the cryptic "NA/NaN/Inf in 'y'" error.
        rank_problem <- diagnose_design_rank(design)
        if (!is.null(rank_problem)) {
          message("[DE-LIMP] ", rank_problem)
          showNotification(
            paste0("Differential expression skipped — ", rank_problem,
                   ". Untick the offending covariate(s) in the sidebar (often Run order, ",
                   "or a covariate with a level that appears in only one sample) and click Run Pipeline again. ",
                   "QC, Expression Grid, and PCA are still available."),
            type = "error", duration = NULL)
          values$status <- "⚠ DE skipped (rank-deficient design)"
          return(invisible(NULL))
        }

        combs <- combn(levels(groups), 2)
        forms <- apply(combs, 2, function(x) paste(x[2], "-", x[1]))

        # Check if DE analysis is possible (need >= 2 replicates in at least one group)
        group_sizes <- table(groups)
        has_replicates <- any(group_sizes >= 2)

        if (!has_replicates) {
          showNotification(
            paste0("No replicates detected (", length(group_sizes), " groups, 1 sample each). ",
                   "Skipping differential expression — QC, Expression Grid, Signal Distribution, ",
                   "and PCA are still available."),
            type = "warning", duration = NULL)
          message("[DE-LIMP] Skipping DE: no replicates (", paste(names(group_sizes), "=",
                  group_sizes, collapse = ", "), ")")
        } else {
          tryCatch({
            if (isTRUE(values$pipeline_mode_used == "maxlfq")) {
              # Paper-faithful MaxLFQ + limma path. Apply coverage filter
              # (UC Davis Bioinformatics Core's recommendation, also reviewer
              # request HIGH #4) so eBayes isn't moderating against rows with
              # only 1-2 finite values. Default 50% of samples non-NA.
              cov_frac <- input$coverage_min_frac %||% 0.5
              n_samples <- ncol(values$y_protein$E)
              min_obs <- max(2, ceiling(cov_frac * n_samples))
              n_obs_per_row <- rowSums(!is.na(values$y_protein$E))
              keep <- n_obs_per_row >= min_obs
              n_dropped <- sum(!keep)
              message(sprintf("[DE-LIMP] MaxLFQ coverage filter: keep proteins with ≥ %d / %d non-NA (%.0f%%). Kept %d, dropped %d to On/Off panel only.",
                              min_obs, n_samples, 100 * cov_frac,
                              sum(keep), n_dropped))
              if (sum(keep) < 10) {
                stop(sprintf("Coverage filter left only %d testable proteins (threshold: ≥ %d non-NA). Loosen the QuantUMS cutoffs or the coverage filter.",
                             sum(keep), min_obs))
              }
              E_for_fit <- values$y_protein$E[keep, , drop = FALSE]
              values$maxlfq_dropped_for_coverage <- n_dropped
              message("[DE-LIMP] Running plain limma::lmFit on MaxLFQ matrix (paper-faithful, no DPC-Quant).")
              fit <- limma::lmFit(E_for_fit, design)
            } else {
              fit <- limpa::dpcDE(values$y_protein, design, plot=FALSE)
            }
            fit <- contrasts.fit(fit, makeContrasts(contrasts=forms, levels=design))
            fit <- eBayes(fit)
            values$fit <- fit

            # On/Off proteins (only meaningful when matrix has missing values,
            # i.e. when the MaxLFQ pipeline is in use).
            values$onoff_proteins <- tryCatch({
              if (any(is.na(values$y_protein$E))) {
                gene_lookup <- if (!is.null(values$y_protein$genes$Genes)) {
                  setNames(values$y_protein$genes$Genes,
                           values$y_protein$genes$Protein.Group)
                } else NULL
                # Pass the contrast matrix directly so we don't depend on
                # parsing "X - Y" strings — Blocker #3 in v3.9.1 review.
                # `forms` was built as paste(x[2], "-", x[1]), so the limma
                # contrast is g2 - g1. Flip rows of `combs` so the on/off
                # function sees (g2, g1) and emits Contrast = "g2 - g1"
                # matching limma's convention.
                combs_flipped <- combs[c(2L, 1L), , drop = FALSE]
                compute_onoff_proteins(values$y_protein$E, groups,
                                       contrasts_list = combs_flipped,
                                       n_min = input$onoff_min_n %||% 2,
                                       gene_lookup = gene_lookup)
              } else NULL
            }, error = function(e) {
              message("[DE-LIMP] On/Off computation skipped: ", e$message); NULL
            })

            # Clear stale GSEA cache from previous pipeline run
            values$gsea_results_cache <- list()
            values$gsea_last_contrast <- NULL

            # Update all four comparison selectors
            updateSelectInput(session, "contrast_selector", choices=forms)
            updateSelectInput(session, "contrast_selector_signal", choices=forms, selected=forms[1])
            updateSelectInput(session, "contrast_selector_grid", choices=forms, selected=forms[1])
            updateSelectInput(session, "contrast_selector_pvalue", choices=forms, selected=forms[1])

            # Log contrasts
            contrast_code <- c(
              sprintf("# Available contrasts: %s", paste(forms, collapse=", ")),
              "combs <- combn(levels(groups), 2)",
              "forms <- apply(combs, 2, function(x) paste(x[2], '-', x[1]))",
              "fit <- contrasts.fit(fit, makeContrasts(contrasts = forms, levels = design))",
              "fit <- eBayes(fit)",
              "",
              "# Per-contrast DE tables",
              "for (cname in colnames(fit$contrasts)) {",
              "  tt <- limma::topTable(fit, coef = cname, number = Inf, adjust.method = 'BH')",
              "  write.csv(tt, paste0('DE_', make.names(cname), '.csv'), row.names = TRUE)",
              "}"
            )
            add_to_log("Contrast Fitting", contrast_code)

            # === Provenance block (v3.9.13) ===
            # End-of-pipeline record: parquet MD5, R sessionInfo, DE-LIMP version,
            # timestamp. Anyone re-running the script can confirm they have the
            # exact same input and R environment as the original analysis.
            tryCatch({
              report_path_prov <- values$uploaded_report_path %||% NA_character_
              parquet_md5 <- if (!is.na(report_path_prov) && file.exists(report_path_prov)) {
                tryCatch(unname(tools::md5sum(report_path_prov)),
                         error = function(e) NA_character_)
              } else NA_character_
              parquet_size <- if (!is.na(report_path_prov) && file.exists(report_path_prov))
                file.info(report_path_prov)$size else NA_integer_
              sess_lines <- tryCatch(
                capture.output(print(utils::sessionInfo())),
                error = function(e) c("(sessionInfo capture failed: ", e$message, ")")
              )
              prov_lines <- c(
                "# ─── Provenance ───────────────────────────────────────────────",
                sprintf("# DE-LIMP version    : %s", values$app_version %||% "unknown"),
                sprintf("# Pipeline used      : %s",
                        values$pipeline_mode_used %||% "dpc"),
                sprintf("# Timestamp          : %s",
                        format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
                sprintf("# Input parquet path : %s",
                        if (is.na(report_path_prov)) "(not recorded)" else report_path_prov),
                sprintf("# Input parquet MD5  : %s",
                        if (is.na(parquet_md5)) "(not computed)" else parquet_md5),
                sprintf("# Input parquet size : %s bytes",
                        if (is.na(parquet_size)) "(unknown)" else format(parquet_size, big.mark = ",")),
                "# To verify the input file matches:",
                "#   tools::md5sum('report.parquet')",
                "",
                "# R sessionInfo() at the time of the run:",
                paste("#  ", sess_lines)
              )
              add_to_log("Provenance", prov_lines)
            }, error = function(e) {
              message("[DE-LIMP] Provenance block could not be written: ", e$message)
            })
          }, error = function(e) {
            message("[DE-LIMP] DE fitting failed: ", e$message)
            showNotification(
              paste0("Differential expression failed: ", e$message,
                     ". QC, Expression Grid, and PCA are still available."),
              type = "warning", duration = NULL)
          })
        }  # end has_replicates

        values$status <- "\u2705 Complete!"

        if (!is.null(values$fit)) {
          if (isTRUE(values$phospho_detected$detected)) {
            nav_select("main_tabs", "Phosphoproteomics")
          } else {
            nav_select("main_tabs", "DE Dashboard")
          }
          showNotification("\u2713 Pipeline complete! View results in tabs below.", type="message", duration=10)
        } else {
          nav_select("main_tabs", "Data Overview")
          showNotification("\u2713 Quantification complete! View Expression Grid and QC tabs.", type="message", duration=10)
        }

        # Auto-save session .rds + record to activity log
        tryCatch({
          ss <- values$diann_search_settings
          out_dir <- if (!is.null(ss)) ss$output_dir else NA

          # Auto-save session .rds — deterministic path: {output_dir}/session.rds
          rds_path <- NA
          tryCatch({
            session_data <- list(
              raw_data = values$raw_data, metadata = values$metadata,
              fit = values$fit, y_protein = values$y_protein,
              dpc_fit = values$dpc_fit, design = values$design,
              qc_stats = values$qc_stats,
              gsea_results = values$gsea_results,
              gsea_results_cache = values$gsea_results_cache,
              repro_log = values$repro_log,
              color_plot_by_de = values$color_plot_by_de,
              contrast = input$contrast_selector,
              logfc_cutoff = input$logfc_cutoff, q_cutoff = input$q_cutoff,
              phospho_detected = values$phospho_detected,
              phospho_site_matrix = values$phospho_site_matrix,
              phospho_site_info = values$phospho_site_info,
              phospho_fit = values$phospho_fit,
              diann_search_settings = ss,
              instrument_metadata = values$instrument_metadata,
              tic_traces = values$tic_traces, tic_metrics = values$tic_metrics,
              saved_at = Sys.time(),
              app_version = paste0("DE-LIMP v", values$app_version),
              auto_saved = TRUE
            )

            # Save to temp file first (staging for SCP upload)
            local_rds <- tempfile(pattern = "delimp_session_", fileext = ".rds")
            saveRDS(session_data, local_rds)
            on.exit(unlink(local_rds), add = TRUE)

            # Save to {output_dir}/session.rds — remote via SCP or local copy
            cfg <- if (nzchar(input$ssh_host %||% "") && nzchar(input$ssh_user %||% ""))
              list(host = input$ssh_host, user = input$ssh_user,
                   port = input$ssh_port %||% 22, key_path = input$ssh_key_path)
            else NULL
            if (!is.null(cfg) && isTRUE(values$ssh_connected) &&
                !is.na(out_dir) && nzchar(out_dir %||% "")) {
              # Upload to remote {output_dir}/session.rds via SCP
              remote_rds <- file.path(out_dir, "session.rds")
              ul <- scp_upload(cfg, local_rds, remote_rds)
              if (ul$status == 0) {
                rds_path <- remote_rds
                message("[DE-LIMP] Auto-saved session to remote: ", remote_rds)
              } else {
                rds_path <- NA
                message("[DE-LIMP] SCP upload failed for session auto-save")
              }
            } else if (!is.na(out_dir) && nzchar(out_dir %||% "") && dir.exists(out_dir)) {
              # Local output_dir exists (Docker/local backend)
              rds_path <- file.path(out_dir, "session.rds")
              file.copy(local_rds, rds_path, overwrite = TRUE)
              message("[DE-LIMP] Auto-saved session: ", rds_path)
            } else {
              rds_path <- NA
              message("[DE-LIMP] No output_dir available for session auto-save")
            }
          }, error = function(e) {
            message("[DE-LIMP] Auto-save failed: ", e$message)
            rds_path <- NA
          })

          record_activity(list(
            event_type = "analysis_completed",
            timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            user = Sys.getenv("USER", "unknown"),
            search_name = if (!is.null(ss)) ss$search_name %||% values$original_report_name else values$original_report_name %||% NA,
            fasta_files = if (!is.null(ss)) paste(basename(ss$fasta_files), collapse = ", ") else NA,
            fasta_seq_count = if (!is.null(ss)) ss$fasta_seq_count else NA,
            n_proteins = nrow(values$y_protein$E),
            n_samples = ncol(values$y_protein$E),
            n_contrasts = if (!is.null(values$fit)) length(colnames(values$fit$contrasts)) else 0L,
            n_de_proteins = if (!is.null(values$fit)) count_de_proteins(values$fit) else 0L,
            output_dir = out_dir,
            session_file = if (!is.na(rds_path)) rds_path else NA,
            app_version = values$app_version %||% "unknown",
            source_type = if (!is.null(ss)) "search" else "upload"
          ))
        }, error = function(e) message("[DE-LIMP] Activity log record failed: ", e$message))

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
