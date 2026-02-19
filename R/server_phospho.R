# ==============================================================================
#  SERVER MODULE — Phosphoproteomics (Phase 1)
#  Site-level DE, volcano, site table, residue distribution, completeness QC
#  Called from app.R as: server_phospho(input, output, session, values, add_to_log)
# ==============================================================================

server_phospho <- function(input, output, session, values, add_to_log) {

  # ============================================================================
  #  Flag for conditionalPanel (sidebar controls)
  # ============================================================================
  output$phospho_detected_flag <- reactive({
    !is.null(values$phospho_detected) && isTRUE(values$phospho_detected$detected)
  })
  outputOptions(output, "phospho_detected_flag", suspendWhenHidden = FALSE)

  # ============================================================================
  #  Detection banner in Data Overview
  # ============================================================================
  output$phospho_detection_banner <- renderUI({
    req(values$phospho_detected)
    pd <- values$phospho_detected
    if (!pd$detected) return(NULL)

    enrichment_msg <- if (pd$is_enriched) {
      sprintf("Phospho-enriched dataset (%s%% phosphopeptides).", pd$pct_phospho)
    } else {
      sprintf("Phosphopeptides detected (%s%% of precursors).", pd$pct_phospho)
    }

    # Show different guidance depending on whether the pipeline has been run
    pipeline_done <- !is.null(values$phospho_fit)

    tags$div(
      class = if (pipeline_done) "alert alert-success py-2 px-3 mb-3" else "alert alert-info py-2 px-3 mb-3",
      style = if (pipeline_done) "border-left: 4px solid #198754;" else "border-left: 4px solid #0d6efd;",
      tags$div(
        style = "display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px;",
        tags$span(
          icon("flask"), tags$strong(" Phosphoproteomics Data Detected — "),
          enrichment_msg,
          sprintf(" %s phosphoprecursors across %s total.",
                  format(pd$n_phospho, big.mark = ","),
                  format(pd$n_total, big.mark = ","))
        ),
        actionButton("go_to_phospho",
                      if (pipeline_done) "View Phospho Results \u2192" else "Open Phospho Tab \u2192",
                      class = if (pipeline_done) "btn btn-sm btn-success" else "btn btn-sm btn-outline-primary")
      ),
      if (pipeline_done) {
        tags$div(
          style = "margin-top: 8px; font-size: 0.85em;",
          icon("check-circle"),
          " Pipeline complete! Visit the ", tags$strong("Phosphoproteomics tab"),
          " to explore site-level volcano plots, kinase activity (KSEA), ",
          "motif analysis, and site annotations."
        )
      } else {
        tags$div(
          style = "margin-top: 8px; font-size: 0.85em; color: #6c757d;",
          icon("info-circle"),
          " Run the pipeline below, then visit the ", tags$strong("Phosphoproteomics tab"),
          " for site-level differential expression, kinase activity inference, and motif analysis."
        )
      }
    )
  })

  # Navigate to phospho tab on banner button click
  observeEvent(input$go_to_phospho, {
    nav_select("main_tabs", "Phosphoproteomics")
  })

  # ============================================================================
  #  Load Example Phospho Data (downloads from GitHub release)
  # ============================================================================
  observeEvent(input$load_example_phospho, {
    withProgress(message = "Loading example phospho data...", {
      tryCatch({
        # Download the phospho report (for detection + parsed_report path)
        incProgress(0.1, detail = "Downloading phospho report...")
        report_url <- "https://github.com/bsphinney/DE-LIMP/releases/download/v1.0/phospho.parquet"
        report_tmp <- tempfile(fileext = ".parquet")
        download.file(report_url, report_tmp, mode = "wb", quiet = TRUE)

        incProgress(0.4, detail = "Downloading site matrix...")
        # Download the site localization TSV (90% confidence)
        site_url <- "https://github.com/bsphinney/DE-LIMP/releases/download/v1.0/report.phosphosites_90.tsv"
        site_tmp <- tempfile(fileext = ".tsv")
        download.file(site_url, site_tmp, mode = "wb", quiet = TRUE)

        incProgress(0.6, detail = "Loading report...")
        # Load main report into pipeline
        values$raw_data <- limpa::readDIANN(report_tmp, format = "parquet",
                                             q.cutoffs = input$q_cutoff)
        fnames <- sort(colnames(values$raw_data$E))
        values$metadata <- data.frame(
          ID = seq_along(fnames),
          File.Name = fnames,
          Group = rep("", length(fnames)),
          Batch = rep("", length(fnames)),
          Covariate1 = rep("", length(fnames)),
          Covariate2 = rep("", length(fnames)),
          stringsAsFactors = FALSE
        )

        # Store report path
        session_report <- file.path(tempdir(), "de_limp_phospho_report.parquet")
        file.copy(report_tmp, session_report, overwrite = TRUE)
        values$uploaded_report_path <- session_report
        values$original_report_name <- "phospho.parquet"

        # Detect phospho data
        values$phospho_detected <- detect_phospho(session_report)

        incProgress(0.7, detail = "Downloading FASTA...")
        # Download FASTA for motif analysis
        fasta_url <- "https://github.com/bsphinney/DE-LIMP/releases/download/v1.0/UP000005640_9606-2.fasta"
        fasta_tmp <- tempfile(fileext = ".fasta")
        tryCatch({
          download.file(fasta_url, fasta_tmp, mode = "wb", quiet = TRUE)
          values$phospho_fasta_sequences <- read_fasta_sequences(fasta_tmp)
        }, error = function(e) {
          showNotification(
            paste("FASTA download failed (motif analysis unavailable):", e$message),
            type = "warning", duration = 8)
        })

        incProgress(0.8, detail = "Loading site matrix...")
        # Parse the site matrix TSV
        mat_df <- utils::read.delim(site_tmp, check.names = FALSE, stringsAsFactors = FALSE)

        # Identify numeric vs info columns by checking actual content
        is_numeric_col <- vapply(mat_df, function(col) {
          is.numeric(col) || all(grepl("^[0-9eE.+\\-]*$", col[!is.na(col) & col != ""]))
        }, logical(1))
        numeric_cols <- names(mat_df)[is_numeric_col]
        info_col_names <- names(mat_df)[!is_numeric_col]

        mat <- as.matrix(mat_df[, numeric_cols, drop = FALSE])
        storage.mode(mat) <- "double"
        mat[mat == 0] <- NA
        mat <- log2(mat)

        # Build site info from non-numeric columns
        site_info <- if (length(info_col_names) > 0) {
          mat_df[, info_col_names, drop = FALSE]
        } else {
          data.frame(row.names = seq_len(nrow(mat)))
        }

        # Map column names to standard names for downstream compatibility
        col_map <- c(
          "Protein" = "Protein.Group", "Protein.Group" = "Protein.Group",
          "Gene.Names" = "Genes", "Genes" = "Genes",
          "Protein.Names" = "Protein.Names",
          "Modified.Sequence" = "Modified.Sequence",
          "Sequence" = "Stripped.Sequence", "Stripped.Sequence" = "Stripped.Sequence"
        )
        for (old_name in names(col_map)) {
          if (old_name %in% names(site_info) && old_name != col_map[[old_name]]) {
            names(site_info)[names(site_info) == old_name] <- col_map[[old_name]]
          }
        }

        # Build SiteIDs from Residue + Site columns (DIA-NN format) or parse from sequences
        if (all(c("Residue", "Site") %in% names(site_info))) {
          protein_col <- if ("Protein.Group" %in% names(site_info)) site_info$Protein.Group else "Unknown"
          site_info$Position <- site_info$Site
          site_info$SiteID <- paste0(protein_col, "_", site_info$Residue, site_info$Site)
        } else if ("Modified.Sequence" %in% names(site_info)) {
          parsed <- parse_phospho_positions(site_info$Modified.Sequence)
          site_info$Residue <- parsed$residue
          site_info$Position <- parsed$position
          site_info$SiteID <- paste0(
            if ("Protein.Group" %in% names(site_info)) site_info$Protein.Group else "Unknown",
            "_", parsed$residue, parsed$position
          )
        } else {
          site_info$SiteID <- paste0("Site_", seq_len(nrow(mat)))
          site_info$Residue <- NA
          site_info$Position <- NA
        }

        rownames(mat) <- site_info$SiteID
        values$phospho_site_matrix <- mat
        values$phospho_site_info <- site_info
        values$phospho_input_mode <- "site_matrix"

        incProgress(1.0, detail = "Complete!")

        add_to_log("Example Phospho Data Loaded", c(
          "# Example phospho dataset downloaded from GitHub releases v1.0",
          "# Report: phospho.parquet",
          "# Site matrix: report.phosphosites_90.tsv (90% localization confidence)",
          "# FASTA: UP000005640_9606-2.fasta (human, for motif analysis)",
          sprintf("# %d sites x %d samples", nrow(mat), ncol(mat))
        ))

        showNotification(
          sprintf("Example phospho data loaded: %d sites x %d samples. Assign groups and run the pipeline!",
                  nrow(mat), ncol(mat)),
          type = "message", duration = 8
        )

        nav_select("main_tabs", "Data Overview")
        nav_select("data_overview_tabs", "Assign Groups & Run")

      }, error = function(e) {
        showNotification(paste("Error loading example phospho data:", e$message),
          type = "error", duration = 10)
      })
    })
  })

  # ============================================================================
  #  Site matrix upload (Path A: DIA-NN 1.9+ site matrix)
  # ============================================================================
  observeEvent(input$phospho_site_matrix_file, {
    req(input$phospho_site_matrix_file)
    tryCatch({
      # Detect file format from original filename extension
      orig_name <- tolower(input$phospho_site_matrix_file$name)
      if (grepl("\\.(tsv|txt)$", orig_name)) {
        mat_df <- utils::read.delim(input$phospho_site_matrix_file$datapath,
                                     check.names = FALSE, stringsAsFactors = FALSE)
      } else {
        mat_df <- arrow::read_parquet(input$phospho_site_matrix_file$datapath,
                                       as_data_frame = TRUE)
      }

      # First column is SiteID / row index; rest are sample intensities
      site_ids <- mat_df[[1]]
      mat <- as.matrix(mat_df[, -1, drop = FALSE])
      rownames(mat) <- site_ids

      # Log2 transform if values appear to be on linear scale
      if (median(mat, na.rm = TRUE) > 100) {
        mat <- log2(mat)
        mat[is.infinite(mat)] <- NA
      }

      values$phospho_site_matrix <- mat
      values$phospho_input_mode  <- "site_matrix"

      # Parse site info from SiteID format (ProteinGroup_Residue_Position)
      site_info <- data.frame(SiteID = site_ids, stringsAsFactors = FALSE)
      # Attempt to parse components from SiteID
      parts <- strsplit(site_ids, "_")
      site_info$Protein.Group <- vapply(parts, function(p) {
        if (length(p) >= 2) paste(p[1:(length(p)-1)], collapse = "_") else p[1]
      }, character(1))
      site_info$Residue <- gsub("[0-9]", "", vapply(parts, function(p) p[length(p)], character(1)))
      site_info$Position <- as.integer(gsub("[^0-9]", "", vapply(parts, function(p) p[length(p)], character(1))))
      site_info$Genes <- NA_character_  # Not available from matrix alone
      values$phospho_site_info <- site_info

      showNotification(
        sprintf("Site matrix loaded: %d sites x %d samples",
                nrow(mat), ncol(mat)),
        type = "message", duration = 5
      )
    }, error = function(e) {
      showNotification(paste("Error reading site matrix:", e$message),
                       type = "error", duration = 10)
    })
  })

  # ============================================================================
  #  Run Phosphosite Analysis pipeline
  # ============================================================================
  observeEvent(input$run_phospho_pipeline, {
    req(values$metadata)

    # Validate groups
    meta <- values$metadata
    meta$Group <- trimws(meta$Group)
    if (length(unique(meta$Group[meta$Group != ""])) < 2) {
      showNotification("Need at least 2 groups. Run the main pipeline first.", type = "error")
      return()
    }

    withProgress(message = "Running phosphosite analysis...", {
      incProgress(0.1, detail = "Preparing site matrix...")

      # --- Get or build site matrix ---
      if (input$phospho_input_mode == "parsed_report") {
        req(values$uploaded_report_path)
        loc_thresh <- input$phospho_loc_threshold
        result <- extract_phosphosites(
          values$uploaded_report_path,
          loc_threshold = loc_thresh,
          q_cutoff = 0.01
        )

        if (is.null(result$matrix)) {
          showNotification(
            paste("Phosphosite extraction failed:", result$message),
            type = "error", duration = 10)
          return()
        }

        values$phospho_site_matrix <- result$matrix
        values$phospho_site_info   <- result$info
        values$phospho_input_mode  <- "parsed_report"

        showNotification(
          sprintf("Extracted %d phosphosites from report",
                  nrow(result$matrix)),
          type = "message", duration = 4)
      }

      req(values$phospho_site_matrix)
      mat <- values$phospho_site_matrix

      # --- Match matrix columns to metadata ---
      # Sample names in matrix may not exactly match metadata$File.Name
      # Try to match by finding common samples
      meta_samples <- meta$File.Name[meta$Group != ""]
      mat_samples  <- colnames(mat)

      common <- intersect(mat_samples, meta_samples)
      if (length(common) == 0) {
        # Try matching by basename (strip path)
        mat_base  <- basename(mat_samples)
        meta_base <- basename(meta_samples)
        common_base <- intersect(mat_base, meta_base)
        if (length(common_base) > 0) {
          mat_idx  <- match(common_base, mat_base)
          meta_idx <- match(common_base, meta_base)
          mat <- mat[, mat_idx, drop = FALSE]
          colnames(mat) <- meta_samples[meta_idx]
          common <- meta_samples[meta_idx]
        }
      }
      if (length(common) == 0) {
        # Try matching by basename with extensions stripped (.d, .raw, .mzML)
        mat_clean  <- sub("\\.(d|raw|mzML)$", "", basename(mat_samples), ignore.case = TRUE)
        meta_clean <- sub("\\.(d|raw|mzML)$", "", basename(meta_samples), ignore.case = TRUE)
        common_clean <- intersect(mat_clean, meta_clean)
        if (length(common_clean) > 0) {
          mat_idx  <- match(common_clean, mat_clean)
          meta_idx <- match(common_clean, meta_clean)
          mat <- mat[, mat_idx, drop = FALSE]
          colnames(mat) <- meta_samples[meta_idx]
          common <- meta_samples[meta_idx]
        }
      }

      if (length(common) < 3) {
        showNotification(
          sprintf("Only %d samples match between site matrix and metadata. Need at least 3.",
                  length(common)),
          type = "error", duration = 10)
        return()
      }

      # Subset and align
      mat  <- mat[, common, drop = FALSE]
      meta_sub <- meta[match(common, meta$File.Name), ]
      groups <- factor(make.names(meta_sub$Group))

      incProgress(0.3, detail = "Filtering sites...")

      # --- Filter: require >=2 non-NA per group ---
      keep <- apply(mat, 1, function(row) {
        all(tapply(!is.na(row), groups, sum) >= 2)
      })
      n_removed <- sum(!keep)
      mat_f <- mat[keep, , drop = FALSE]

      if (nrow(mat_f) < 10) {
        showNotification(
          sprintf("Only %d sites passed filtering. Need at least 10.", nrow(mat_f)),
          type = "error", duration = 10)
        return()
      }

      incProgress(0.4, detail = "Imputing missing values...")

      # --- Tail-based imputation (Perseus-style) ---
      if (any(is.na(mat_f))) {
        gm <- mean(mat_f, na.rm = TRUE)
        gs <- sd(mat_f, na.rm = TRUE)
        set.seed(42)
        na_idx <- which(is.na(mat_f), arr.ind = TRUE)
        mat_f[na_idx] <- rnorm(nrow(na_idx), mean = gm - 1.8 * gs, sd = 0.3 * gs)
      }

      incProgress(0.5, detail = "Normalizing...")

      # --- Optional normalization ---
      norm_method <- input$phospho_norm
      if (!is.null(norm_method) && norm_method == "median") {
        col_meds <- apply(mat_f, 2, median, na.rm = TRUE)
        mat_f <- sweep(mat_f, 2, col_meds - mean(col_meds))
      } else if (!is.null(norm_method) && norm_method == "quantile") {
        mat_f <- limma::normalizeBetweenArrays(mat_f, method = "quantile")
      }

      incProgress(0.6, detail = "Fitting linear model...")

      # --- limma DE ---
      design <- model.matrix(~ 0 + groups)
      colnames(design) <- levels(groups)

      fit <- limma::lmFit(mat_f, design)
      combs <- combn(levels(groups), 2)
      forms <- apply(combs, 2, function(x) paste(x[2], "-", x[1]))

      contrasts_mat <- limma::makeContrasts(contrasts = forms, levels = design)
      fit <- limma::contrasts.fit(fit, contrasts_mat)
      fit <- limma::eBayes(fit)

      values$phospho_fit <- fit
      values$phospho_site_matrix_filtered <- mat_f

      # Update contrast selector
      updateSelectInput(session, "phospho_contrast_selector", choices = forms)

      incProgress(1.0, detail = "Complete!")

      showNotification(
        sprintf("\u2713 Phosphosite DE complete: %d sites tested, %d contrasts",
                nrow(mat_f), length(forms)),
        type = "message", duration = 8
      )

      # Log to reproducibility
      add_to_log("Phosphosite DE Analysis", c(
        sprintf("# Input mode: %s", values$phospho_input_mode),
        sprintf("# Site-level matrix: %d sites x %d samples", nrow(mat_f), ncol(mat_f)),
        sprintf("# Sites removed (missing data): %d", n_removed),
        sprintf("# Normalization: %s", norm_method %||% "none"),
        "# Imputation: tail-based (mean - 1.8 SD, width 0.3 SD)",
        "fit_phospho <- lmFit(site_matrix_filtered, design)",
        sprintf("fit_phospho <- contrasts.fit(fit_phospho, makeContrasts(contrasts=c('%s'), levels=design))",
                paste(forms, collapse = "', '")),
        "fit_phospho <- eBayes(fit_phospho)"
      ))

      # Navigate to phospho tab
      nav_select("main_tabs", "Phosphoproteomics")
    })
  })

  # ============================================================================
  #  Phospho Tab Content (dynamic UI)
  # ============================================================================
  output$phospho_tab_content <- renderUI({
    if (is.null(values$phospho_detected) || !isTRUE(values$phospho_detected$detected)) {
      return(
        tags$div(class = "text-center text-muted py-5",
          icon("flask", class = "fa-3x mb-3"),
          tags$h5("No Phosphoproteomics Data Detected"),
          tags$p("This tab activates when your DIA-NN report contains phosphopeptides ",
                 "(searched with STY phosphorylation as a variable modification)."),
          tags$p("To enable, search in DIA-NN with:"),
          tags$code("--var-mod UniMod:21,79.966331,STY --peptidoforms"),
          tags$p(class = "mt-3 text-muted small",
            "Upload your report.parquet and the app will auto-detect phospho data.")
        )
      )
    }

    # --- Phospho detected: full UI ---
    tagList(
      # Status & controls card
      tags$div(
        style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 15px;",

        uiOutput("phospho_status_banner"),
        uiOutput("phospho_norm_warning"),

        # Contrast selector
        tags$div(
          style = "display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-top: 10px;",
          tags$strong("Phosphosite Contrast:"),
          tags$div(style = "flex-grow: 1; max-width: 400px;",
            selectInput("phospho_contrast_selector", NULL, choices = NULL, width = "100%")
          )
        ),

        # Educational expandable
        tags$details(
          style = "margin-top: 10px;",
          tags$summary(
            style = "cursor: pointer; color: #6c757d; font-size: 0.85em;",
            icon("info-circle"), " About phosphosite-level analysis"
          ),
          tags$div(
            style = "padding: 8px 0; font-size: 0.85em; color: #6c757d;",
            tags$p(
              "Unlike standard proteomics where peptides are aggregated into ",
              "protein-level measurements, phosphoproteomics requires ",
              tags$strong("site-level analysis"), ". A single protein can have ",
              "dozens of phosphorylation sites, each independently regulated."
            ),
            tags$p(
              tags$strong("Site localization confidence"), " indicates certainty that ",
              "the phosphate is on the correct residue. Class I sites (\u2265 0.75) are ",
              "reliably localized."
            ),
            tags$p(
              tags$strong("Imputation"), ": Missing values are filled using a ",
              "downshifted normal distribution (tail-based, Perseus default). ",
              "This assumes missing = below detection limit."
            )
          )
        )
      ),

      # Results tabs
      navset_card_tab(
        id = "phospho_results_tabs",

        nav_panel("Phospho Volcano",
          tags$div(style = "display: flex; justify-content: flex-end; gap: 5px; margin-bottom: 10px;",
            actionButton("phospho_volcano_info_btn", NULL, icon = icon("question-circle"),
                         class = "btn-outline-info btn-sm", title = "About this plot"),
            actionButton("phospho_volcano_fullscreen_btn", NULL, icon = icon("expand"),
                         class = "btn-outline-secondary btn-sm", title = "Fullscreen"),
            downloadButton("download_phospho_volcano", "Save Plot",
                           class = "btn-outline-secondary btn-sm")
          ),
          plotOutput("phospho_volcano", height = "550px")
        ),

        nav_panel("Site Table",
          tags$div(style = "display: flex; justify-content: flex-end; gap: 5px; margin-bottom: 10px;",
            actionButton("phospho_table_info_btn", NULL, icon = icon("question-circle"),
                         class = "btn-outline-info btn-sm", title = "About this table"),
            downloadButton("download_phospho_table", "Export CSV",
                           class = "btn-success btn-sm")
          ),
          DT::DTOutput("phospho_site_table")
        ),

        nav_panel("Residue Distribution",
          tags$div(style = "display: flex; justify-content: flex-end; gap: 5px; margin-bottom: 10px;",
            actionButton("phospho_residue_info_btn", NULL, icon = icon("question-circle"),
                         class = "btn-outline-info btn-sm", title = "About this plot"),
            actionButton("phospho_residue_fullscreen_btn", NULL, icon = icon("expand"),
                         class = "btn-outline-secondary btn-sm", title = "Fullscreen")
          ),
          plotOutput("phospho_residue_dist", height = "450px")
        ),

        nav_panel("QC: Completeness",
          tags$div(style = "display: flex; justify-content: flex-end; gap: 5px; margin-bottom: 10px;",
            actionButton("phospho_completeness_info_btn", NULL, icon = icon("question-circle"),
                         class = "btn-outline-info btn-sm", title = "About this plot"),
            actionButton("phospho_completeness_fullscreen_btn", NULL, icon = icon("expand"),
                         class = "btn-outline-secondary btn-sm", title = "Fullscreen")
          ),
          plotOutput("phospho_completeness", height = "450px")
        ),

        nav_panel("Kinase Activity",
          tags$div(
            style = "display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-bottom: 15px;",
            actionButton("run_ksea", "Run Kinase Enrichment (KSEA)",
                         class = "btn-info", icon = icon("project-diagram")),
            tags$span(id = "ksea_status_text", class = "text-muted small"),
            tags$div(style = "margin-left: auto; display: flex; gap: 5px;",
              actionButton("ksea_info_btn", NULL, icon = icon("question-circle"),
                           class = "btn-outline-info btn-sm", title = "About KSEA"),
              actionButton("ksea_fullscreen_btn", NULL, icon = icon("expand"),
                           class = "btn-outline-secondary btn-sm", title = "Fullscreen")
            )
          ),
          hr(),
          plotOutput("ksea_barplot", height = "500px"),
          hr(),
          tags$div(style = "text-align: right; margin-bottom: 10px;",
            downloadButton("download_ksea_table", "Export CSV",
                           class = "btn-success btn-sm")
          ),
          DT::DTOutput("ksea_results_table")
        ),

        nav_panel("Motif Analysis",
          tags$div(style = "display: flex; justify-content: flex-end; gap: 5px; margin-bottom: 10px;",
            actionButton("phospho_motif_info_btn", NULL, icon = icon("question-circle"),
                         class = "btn-outline-info btn-sm", title = "About motif analysis"),
            actionButton("phospho_motif_fullscreen_btn", NULL, icon = icon("expand"),
                         class = "btn-outline-secondary btn-sm", title = "Fullscreen")
          ),
          plotOutput("phospho_motif_up", height = "220px"),
          plotOutput("phospho_motif_down", height = "220px")
        ),

        nav_panel("Site Annotation",
          tags$div(
            style = "display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-bottom: 15px;",
            actionButton("run_site_annotation", "Query UniProt for Known Sites",
                         class = "btn-info", icon = icon("database")),
            tags$span(class = "text-muted small",
              "Queries UniProt & PhosphoSitePlus to annotate each site as Known or Novel."
            ),
            tags$div(style = "margin-left: auto;",
              actionButton("site_annotation_info_btn", NULL, icon = icon("question-circle"),
                           class = "btn-outline-info btn-sm", title = "About site annotation")
            )
          ),
          uiOutput("site_annotation_summary"),
          tags$p(class = "text-muted small",
            "Known sites have been experimentally confirmed in UniProt (curated literature evidence) ",
            "or appear as kinase substrates in PhosphoSitePlus. Novel sites are not yet in these databases. ",
            "Note: matching uses protein-relative positions (best with Path A / DIA-NN site matrix)."
          ),
          hr(),
          tags$div(style = "text-align: right; margin-bottom: 10px;",
            downloadButton("download_annotation_table", "Export CSV",
                           class = "btn-success btn-sm")
          ),
          DT::DTOutput("site_annotation_table")
        )
      )
    )
  })

  # ============================================================================
  #  Status banner (inside phospho tab)
  # ============================================================================
  output$phospho_status_banner <- renderUI({
    req(values$phospho_detected)
    pd <- values$phospho_detected

    n_sites_tested <- if (!is.null(values$phospho_site_matrix_filtered)) {
      nrow(values$phospho_site_matrix_filtered)
    } else { NULL }

    tags$div(
      class = "alert alert-info py-2 px-3 mb-2",
      style = "font-size: 0.9em; border-left: 4px solid #0d6efd;",
      icon("flask"),
      sprintf(" %s phosphoprecursors detected (%s%% of total). ",
              format(pd$n_phospho, big.mark = ","), pd$pct_phospho),
      if (!is.null(n_sites_tested)) {
        sprintf("%s unique sites tested.", format(n_sites_tested, big.mark = ","))
      } else {
        "Run the phosphosite pipeline to begin analysis."
      }
    )
  })

  # ============================================================================
  #  Normalization warning (phospho-enriched)
  # ============================================================================
  output$phospho_norm_warning <- renderUI({
    req(values$phospho_detected)
    if (!isTRUE(values$phospho_detected$is_enriched)) return(NULL)

    tags$div(
      class = "alert alert-warning py-2 px-3 mb-2",
      style = "font-size: 0.85em;",
      icon("exclamation-triangle"),
      tags$strong(" Phospho-enriched data: "),
      "DIA-NN's normalization assumes most peptides are unchanged. ",
      "For phospho-enriched samples, this may be partially violated. ",
      "Consider 'median centering' if you see systematic biases in QC."
    )
  })

  # ============================================================================
  #  Phospho Volcano Plot
  # ============================================================================
  output$phospho_volcano <- renderPlot({
    req(values$phospho_fit, input$phospho_contrast_selector)

    de <- limma::topTable(values$phospho_fit,
                          coef = input$phospho_contrast_selector,
                          number = Inf, sort.by = "none")
    de$SiteID <- rownames(de)

    # Merge with site info for gene names
    if (!is.null(values$phospho_site_info)) {
      de <- merge(de, values$phospho_site_info, by = "SiteID", all.x = TRUE)
    }

    # Apply protein correction if active
    correction_active <- FALSE
    if (isTRUE(input$phospho_protein_correction) && !is.null(values$fit)) {
      corrected <- phospho_corrected_de()
      if (!is.null(corrected)) {
        # Replace logFC with corrected values
        corr_map <- stats::setNames(corrected$corrected_logFC, corrected$SiteID)
        matched <- corr_map[de$SiteID]
        de$logFC[!is.na(matched)] <- matched[!is.na(matched)]
        correction_active <- TRUE
      }
    }

    # Significance categories
    de$sig <- dplyr::case_when(
      de$adj.P.Val < 0.05 & abs(de$logFC) > 1 ~ "Significant",
      de$adj.P.Val < 0.05 ~ "FDR < 0.05",
      TRUE ~ "NS"
    )

    n_sig  <- sum(de$sig == "Significant")
    n_up   <- sum(de$sig == "Significant" & de$logFC > 0)
    n_down <- sum(de$sig == "Significant" & de$logFC < 0)

    # Label: Gene Residue+Position
    de$label <- ifelse(
      !is.na(de$Genes) & !is.na(de$Residue) & de$Genes != "",
      paste0(de$Genes, " ", de$Residue, de$Position),
      de$SiteID
    )

    # Top labeled hits (safe subset)
    sig_de <- de[de$sig == "Significant", ]
    if (nrow(sig_de) > 0) {
      sig_de <- sig_de[order(sig_de$adj.P.Val), ]
      label_de <- head(sig_de, 15)
    } else {
      label_de <- sig_de[0, ]
    }

    p <- ggplot2::ggplot(de, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
      ggplot2::geom_point(alpha = 0.6, size = 1.5) +
      ggplot2::scale_color_manual(values = c(
        "Significant" = "#E63946",
        "FDR < 0.05"  = "#457B9D",
        "NS"          = "gray70"
      )) +
      ggrepel::geom_text_repel(
        data = label_de,
        ggplot2::aes(label = label),
        size = 3, max.overlaps = 20, color = "black"
      ) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
      ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = paste("Phosphosite Volcano:", input$phospho_contrast_selector,
                      if (correction_active) "(protein-corrected)" else ""),
        subtitle = sprintf("%d sites | %d significant (\u2191%d \u2193%d) at |FC|>2 & FDR<0.05",
                           nrow(de), n_sig, n_up, n_down),
        x = if (correction_active) "Corrected log2 FC (phospho - protein)" else "log2 Fold Change (phosphosite)",
        y = "-log10(adjusted p-value)"
      ) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(legend.position = "bottom",
                     plot.subtitle = ggplot2::element_text(color = "gray40"))

    p
  })

  # Download volcano
  output$download_phospho_volcano <- downloadHandler(
    filename = function() {
      paste0("phospho_volcano_", gsub(" ", "_", input$phospho_contrast_selector), ".pdf")
    },
    content = function(file) {
      grDevices::pdf(file, width = 10, height = 8)
      print(
        # Re-render the volcano (same code as above)
        {
          req(values$phospho_fit, input$phospho_contrast_selector)
          de <- limma::topTable(values$phospho_fit, coef = input$phospho_contrast_selector,
                                number = Inf, sort.by = "none")
          de$SiteID <- rownames(de)
          if (!is.null(values$phospho_site_info)) {
            de <- merge(de, values$phospho_site_info, by = "SiteID", all.x = TRUE)
          }
          de$sig <- dplyr::case_when(
            de$adj.P.Val < 0.05 & abs(de$logFC) > 1 ~ "Significant",
            de$adj.P.Val < 0.05 ~ "FDR < 0.05",
            TRUE ~ "NS"
          )
          de$label <- ifelse(!is.na(de$Genes) & !is.na(de$Residue) & de$Genes != "",
                            paste0(de$Genes, " ", de$Residue, de$Position), de$SiteID)
          sig_de <- de[de$sig == "Significant", ]
          label_de <- if (nrow(sig_de) > 0) head(sig_de[order(sig_de$adj.P.Val), ], 15) else sig_de[0, ]
          ggplot2::ggplot(de, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
            ggplot2::geom_point(alpha = 0.6, size = 1.5) +
            ggplot2::scale_color_manual(values = c("Significant" = "#E63946", "FDR < 0.05" = "#457B9D", "NS" = "gray70")) +
            ggrepel::geom_text_repel(
              data = label_de,
              ggplot2::aes(label = label), size = 3, max.overlaps = 20, color = "black") +
            ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
            ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
            ggplot2::labs(title = paste("Phosphosite Volcano:", input$phospho_contrast_selector),
                         x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
            ggplot2::theme_bw(base_size = 14) +
            ggplot2::theme(legend.position = "bottom")
        }
      )
      grDevices::dev.off()
    }
  )

  # ============================================================================
  #  Site Table
  # ============================================================================
  output$phospho_site_table <- DT::renderDT({
    req(values$phospho_fit, input$phospho_contrast_selector)

    de <- limma::topTable(values$phospho_fit,
                          coef = input$phospho_contrast_selector,
                          number = Inf)
    de$SiteID <- rownames(de)

    if (!is.null(values$phospho_site_info)) {
      de <- merge(de, values$phospho_site_info, by = "SiteID", all.x = TRUE)
    }

    # Select display columns
    display_cols <- c("SiteID", "Genes", "Residue", "Position",
                      "logFC", "AveExpr", "P.Value", "adj.P.Val")
    if ("Best.Loc.Conf" %in% names(de)) {
      display_cols <- c(display_cols, "Best.Loc.Conf")
    }
    display_cols <- intersect(display_cols, names(de))
    de <- de[, display_cols, drop = FALSE]

    # Round numeric columns
    num_cols <- c("logFC", "AveExpr", "P.Value", "adj.P.Val", "Best.Loc.Conf")
    for (col in intersect(num_cols, names(de))) {
      if (col %in% c("P.Value", "adj.P.Val")) {
        de[[col]] <- signif(de[[col]], 3)
      } else {
        de[[col]] <- round(de[[col]], 3)
      }
    }

    DT::datatable(de,
      rownames = FALSE,
      filter = "top",
      selection = "multiple",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(which(names(de) == "adj.P.Val") - 1, "asc"))
      )
    )
  })

  # Download site table
  output$download_phospho_table <- downloadHandler(
    filename = function() {
      paste0("phosphosites_", gsub(" ", "_", input$phospho_contrast_selector), ".csv")
    },
    content = function(file) {
      de <- limma::topTable(values$phospho_fit,
                            coef = input$phospho_contrast_selector,
                            number = Inf)
      de$SiteID <- rownames(de)
      if (!is.null(values$phospho_site_info)) {
        de <- merge(de, values$phospho_site_info, by = "SiteID", all.x = TRUE)
      }
      write.csv(de, file, row.names = FALSE)
    }
  )

  # ============================================================================
  #  Residue Distribution
  # ============================================================================
  output$phospho_residue_dist <- renderPlot({
    req(values$phospho_fit, input$phospho_contrast_selector, values$phospho_site_info)

    de <- limma::topTable(values$phospho_fit,
                          coef = input$phospho_contrast_selector,
                          number = Inf)
    de$SiteID <- rownames(de)
    de <- merge(de, values$phospho_site_info, by = "SiteID")

    # All sites
    all_counts <- table(factor(de$Residue, levels = c("S", "T", "Y")))
    # Significant sites
    sig_de <- de[de$adj.P.Val < 0.05, ]
    sig_counts <- table(factor(sig_de$Residue, levels = c("S", "T", "Y")))

    plot_data <- data.frame(
      Residue  = rep(c("S", "T", "Y"), 2),
      Count    = c(as.numeric(all_counts), as.numeric(sig_counts)),
      Category = rep(c("All quantified", "Significant (FDR < 0.05)"), each = 3)
    )
    plot_data$Count[is.na(plot_data$Count)] <- 0

    ggplot2::ggplot(plot_data, ggplot2::aes(x = Residue, y = Count, fill = Category)) +
      ggplot2::geom_col(position = "dodge", alpha = 0.8) +
      ggplot2::scale_fill_manual(values = c(
        "All quantified"            = "#457B9D",
        "Significant (FDR < 0.05)"  = "#E63946"
      )) +
      ggplot2::labs(
        title = "Phosphosite Residue Distribution",
        subtitle = "Expected: ~85% Ser / ~14% Thr / ~1% Tyr (typical TiO2/IMAC enrichment)",
        x = "Phosphorylated Residue", y = "Number of Sites"
      ) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(legend.position = "bottom")
  })

  # ============================================================================
  #  Completeness QC
  # ============================================================================
  output$phospho_completeness <- renderPlot({
    req(values$phospho_site_matrix)

    mat <- values$phospho_site_matrix
    n_sites   <- nrow(mat)
    n_samples <- ncol(mat)

    site_completeness <- rowSums(!is.na(mat)) / n_samples * 100
    hist_data <- data.frame(Completeness = site_completeness)

    ggplot2::ggplot(hist_data, ggplot2::aes(x = Completeness)) +
      ggplot2::geom_histogram(bins = 20, fill = "#457B9D", color = "white", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
      ggplot2::annotate("text", x = 52, y = Inf, vjust = 2, hjust = 0,
                        label = "50% threshold", color = "red", size = 3.5) +
      ggplot2::labs(
        title = "Phosphosite Quantification Completeness",
        subtitle = sprintf("%s sites across %d samples | Median completeness: %.0f%%",
                           format(n_sites, big.mark = ","), n_samples,
                           median(site_completeness)),
        x = "% of Samples with Quantification", y = "Number of Phosphosites"
      ) +
      ggplot2::theme_bw(base_size = 14)
  })

  # ============================================================================
  #  Phase 2: FASTA Upload
  # ============================================================================
  observeEvent(input$phospho_fasta_file, {
    req(input$phospho_fasta_file)
    tryCatch({
      seqs <- read_fasta_sequences(input$phospho_fasta_file$datapath)
      values$phospho_fasta_sequences <- seqs
      showNotification(
        sprintf("FASTA loaded: %d protein sequences", length(seqs)),
        type = "message", duration = 5
      )
    }, error = function(e) {
      showNotification(paste("Error reading FASTA:", e$message),
                       type = "error", duration = 10)
    })
  })

  # ============================================================================
  #  Phase 2: KSEA Kinase Activity Inference
  # ============================================================================
  observeEvent(input$run_ksea, {
    req(values$phospho_fit, input$phospho_contrast_selector)

    if (!requireNamespace("KSEAapp", quietly = TRUE)) {
      showNotification(
        "KSEAapp package not installed. Run: install.packages('KSEAapp')",
        type = "error", duration = 10)
      return()
    }

    withProgress(message = "Running kinase activity inference...", {
      incProgress(0.2, detail = "Preparing input...")

      ksea_input <- prepare_ksea_input(
        values$phospho_fit,
        input$phospho_contrast_selector,
        values$phospho_site_info
      )

      if (nrow(ksea_input) < 5) {
        showNotification(
          sprintf("Only %d sites with gene annotations. Need more for KSEA.", nrow(ksea_input)),
          type = "error", duration = 8)
        return()
      }

      incProgress(0.5, detail = "Computing KSEA scores...")

      tryCatch({
        # Load PhosphoSitePlus + NetworKIN kinase-substrate database
        ks_data <- NULL
        data("KSData", package = "KSEAapp", envir = environment())
        ks_data <- get("KSData", envir = environment())

        # Suppress KSEAapp's built-in plot output
        grDevices::png(tempfile())
        ksea_scores <- KSEAapp::KSEA.Scores(
          ks_data,
          ksea_input,
          NetworKIN    = TRUE,
          NetworKIN.cutoff = 5
        )
        grDevices::dev.off()

        values$ksea_results <- ksea_scores
        values$ksea_last_contrast <- input$phospho_contrast_selector

        incProgress(1.0, detail = "Complete!")

        n_sig <- sum(ksea_scores$FDR < 0.05, na.rm = TRUE)
        showNotification(
          sprintf("\u2713 KSEA complete: %d kinases scored, %d significant (FDR < 0.05).",
                  nrow(ksea_scores), n_sig),
          type = "message", duration = 8
        )

        add_to_log("KSEA Kinase Activity", c(
          sprintf("# Contrast: %s", input$phospho_contrast_selector),
          sprintf("# Input: %d phosphosites with gene annotations", nrow(ksea_input)),
          "data('KSData', package = 'KSEAapp')",
          "ksea_scores <- KSEAapp::KSEA.Scores(KSData, ksea_input, NetworKIN=TRUE, NetworKIN.cutoff=5)"
        ))
      }, error = function(e) {
        try(grDevices::dev.off(), silent = TRUE)
        showNotification(paste("KSEA error:", e$message),
                         type = "error", duration = 10)
      })
    })
  })

  # KSEA Bar Plot
  output$ksea_barplot <- renderPlot({
    req(values$ksea_results)

    ks <- values$ksea_results
    ks <- ks[order(ks$z.score, decreasing = TRUE), ]

    # Show significant kinases, or top 20 if none significant
    sig_ks <- ks[ks$FDR < 0.05, ]
    if (nrow(sig_ks) == 0) sig_ks <- utils::head(ks, 20)

    # Take top 15 activated + top 15 inhibited
    top_up   <- utils::head(sig_ks[sig_ks$z.score > 0, ], 15)
    top_down <- utils::tail(sig_ks[sig_ks$z.score < 0, ], 15)
    plot_ks  <- rbind(top_up, top_down)

    if (nrow(plot_ks) == 0) {
      plot.new()
      text(0.5, 0.5, "No kinases scored (insufficient substrate matches)")
      return()
    }

    plot_ks <- plot_ks[order(plot_ks$z.score), ]
    plot_ks$Kinase.Gene <- factor(plot_ks$Kinase.Gene, levels = plot_ks$Kinase.Gene)
    plot_ks$Direction <- ifelse(plot_ks$z.score > 0, "Activated", "Inhibited")

    ggplot2::ggplot(plot_ks, ggplot2::aes(x = z.score, y = Kinase.Gene, fill = Direction)) +
      ggplot2::geom_col(alpha = 0.85) +
      ggplot2::scale_fill_manual(values = c("Activated" = "#E63946", "Inhibited" = "#457B9D")) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("n=%d", m)),
        hjust = ifelse(plot_ks$z.score > 0, -0.1, 1.1), size = 3
      ) +
      ggplot2::labs(
        title = paste("Kinase Activity:", values$ksea_last_contrast),
        subtitle = "KSEA z-scores (PhosphoSitePlus + NetworKIN substrates)",
        x = "KSEA z-score", y = NULL
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(legend.position = "bottom")
  })

  # KSEA Results Table
  output$ksea_results_table <- DT::renderDT({
    req(values$ksea_results)
    ks <- values$ksea_results
    # Round numeric columns
    num_cols <- intersect(c("z.score", "p.value", "FDR", "mS"), names(ks))
    for (col in num_cols) {
      ks[[col]] <- signif(ks[[col]], 3)
    }
    DT::datatable(ks, rownames = FALSE, filter = "top",
      options = list(pageLength = 20, scrollX = TRUE,
        order = list(list(which(names(ks) == "FDR") - 1, "asc")))
    )
  })

  # Download KSEA table
  output$download_ksea_table <- downloadHandler(
    filename = function() {
      paste0("ksea_results_", gsub(" ", "_", values$ksea_last_contrast %||% ""), ".csv")
    },
    content = function(file) {
      req(values$ksea_results)
      write.csv(values$ksea_results, file, row.names = FALSE)
    }
  )

  # ============================================================================
  #  Phase 2: Motif Analysis (ggseqlogo)
  # ============================================================================

  # Compute flanking sequences reactively when FASTA + DE are available
  phospho_flanking_seqs <- reactive({
    req(values$phospho_fit, values$phospho_fasta_sequences,
        values$phospho_site_info, input$phospho_contrast_selector)

    de <- limma::topTable(values$phospho_fit,
                          coef = input$phospho_contrast_selector,
                          number = Inf)
    extract_flanking_sequences(
      values$phospho_site_info, de,
      values$phospho_fasta_sequences, window = 7L
    )
  })

  output$phospho_motif_up <- renderPlot({
    seqs <- phospho_flanking_seqs()

    if (is.null(seqs) || length(seqs$up) < 10) {
      plot.new()
      if (is.null(values$phospho_fasta_sequences)) {
        text(0.5, 0.5, "Upload a FASTA file (sidebar) to enable motif analysis",
             cex = 1.1, col = "gray50")
      } else {
        text(0.5, 0.5,
             sprintf("Insufficient up-regulated sites for motif (need >= 10, have %d)",
                     length(seqs$up)),
             cex = 1.1, col = "gray50")
      }
      return()
    }

    if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
      plot.new()
      text(0.5, 0.5, "ggseqlogo package not installed.\nRun: install.packages('ggseqlogo')",
           cex = 1.1, col = "gray50")
      return()
    }

    ggseqlogo::ggseqlogo(seqs$up, method = "bits", seq_type = "aa") +
      ggplot2::ggtitle(sprintf("Up-regulated phosphosites (n=%d)", length(seqs$up))) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12))
  })

  output$phospho_motif_down <- renderPlot({
    seqs <- phospho_flanking_seqs()

    if (is.null(seqs) || length(seqs$down) < 10) {
      plot.new()
      if (is.null(values$phospho_fasta_sequences)) {
        text(0.5, 0.5, "Upload a FASTA file (sidebar) to enable motif analysis",
             cex = 1.1, col = "gray50")
      } else {
        text(0.5, 0.5,
             sprintf("Insufficient down-regulated sites for motif (need >= 10, have %d)",
                     length(seqs$down)),
             cex = 1.1, col = "gray50")
      }
      return()
    }

    if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
      plot.new()
      text(0.5, 0.5, "ggseqlogo package not installed.\nRun: install.packages('ggseqlogo')",
           cex = 1.1, col = "gray50")
      return()
    }

    ggseqlogo::ggseqlogo(seqs$down, method = "bits", seq_type = "aa") +
      ggplot2::ggtitle(sprintf("Down-regulated phosphosites (n=%d)", length(seqs$down))) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 12))
  })

  # ============================================================================
  #  Phase 3: Protein-Level Abundance Correction
  # ============================================================================

  # Reactive: corrected DE table (only when checkbox is checked + both fits exist)
  phospho_corrected_de <- reactive({
    req(input$phospho_protein_correction,
        values$phospho_fit, values$fit, values$phospho_site_info,
        input$phospho_contrast_selector)

    # Check if the same contrast exists in the protein-level fit
    protein_contrasts <- colnames(values$fit$contrasts)
    phospho_contrast  <- input$phospho_contrast_selector

    if (!phospho_contrast %in% protein_contrasts) return(NULL)

    correct_phospho_for_protein(
      values$phospho_fit, values$fit,
      phospho_contrast, values$phospho_site_info
    )
  })

  # Override volcano with corrected values when correction is active
  observe({
    if (isTRUE(input$phospho_protein_correction) && !is.null(phospho_corrected_de())) {
      values$phospho_corrected_active <- TRUE
    } else {
      values$phospho_corrected_active <- FALSE
    }
  })

  # ============================================================================
  #  Site Annotation — Known vs Novel phosphosites
  # ============================================================================
  observeEvent(input$run_site_annotation, {
    req(values$phospho_site_info)

    site_info <- values$phospho_site_info

    # Extract unique accessions (first accession from Protein.Group)
    accessions <- unique(vapply(
      strsplit(site_info$Protein.Group, "[;]"),
      function(x) trimws(x[1]),
      character(1)
    ))
    accessions <- accessions[!is.na(accessions) & nzchar(accessions)]

    withProgress(message = "Annotating phosphosites...", {
      incProgress(0.05, detail = sprintf("Querying UniProt for %d proteins...", length(accessions)))

      # --- Query UniProt ---
      uniprot_sites <- tryCatch(
        query_uniprot_phosphosites(accessions, progress_fn = function(frac) {
          incProgress(frac * 0.7,
            detail = sprintf("UniProt: %d / %d proteins...",
                             round(frac * length(accessions)), length(accessions)))
        }),
        error = function(e) {
          showNotification(paste("UniProt query error:", e$message),
                           type = "warning", duration = 8)
          data.frame(Accession = character(0), Position = integer(0),
                     Residue = character(0), Description = character(0),
                     Evidence.Code = character(0), Evidence.Source = character(0),
                     Evidence.ID = character(0), stringsAsFactors = FALSE)
        }
      )

      incProgress(0.75, detail = "Checking PhosphoSitePlus...")

      # --- Get PhosphoSitePlus sites ---
      psp_sites <- get_psp_known_sites()

      incProgress(0.85, detail = "Building annotation table...")

      # --- Get DE results if available ---
      de_results <- NULL
      if (!is.null(values$phospho_fit) && !is.null(input$phospho_contrast_selector)) {
        de_results <- tryCatch(
          limma::topTable(values$phospho_fit,
                          coef = input$phospho_contrast_selector,
                          number = Inf),
          error = function(e) NULL
        )
      }

      # --- Build annotated table ---
      ann_table <- build_site_annotation_table(
        site_info, de_results, uniprot_sites, psp_sites
      )

      values$phospho_annotations <- ann_table

      incProgress(1.0, detail = "Complete!")

      n_known <- sum(ann_table$Status == "Known")
      n_total <- nrow(ann_table)
      showNotification(
        sprintf("\u2713 Annotation complete: %d / %d sites are known (%.0f%%).",
                n_known, n_total, 100 * n_known / max(n_total, 1)),
        type = "message", duration = 8
      )

      add_to_log("Phosphosite Annotation", c(
        sprintf("# Queried %d unique accessions against UniProt REST API", length(accessions)),
        sprintf("# UniProt: %d known phosphosites found", nrow(uniprot_sites)),
        sprintf("# PhosphoSitePlus (KSEAapp): %d substrate sites in database", nrow(psp_sites)),
        sprintf("# Result: %d known / %d total sites", n_known, n_total)
      ))
    })
  })

  # Summary banner
  output$site_annotation_summary <- renderUI({
    req(values$phospho_annotations)
    ann <- values$phospho_annotations

    n_total <- nrow(ann)
    n_known <- sum(ann$Status == "Known")
    n_novel <- n_total - n_known
    pct_known <- round(100 * n_known / max(n_total, 1))

    tags$div(
      class = "alert alert-success py-2 px-3 mb-3",
      style = "font-size: 0.9em;",
      tags$div(
        style = "display: flex; gap: 20px; align-items: center; flex-wrap: wrap;",
        tags$span(
          icon("check-circle"),
          sprintf(" %d of %d sites annotated (%d%% known)", n_known, n_total, pct_known)
        ),
        tags$span(
          style = "color: #2e7d32; font-weight: bold;",
          sprintf("\u2705 Known: %d", n_known)
        ),
        tags$span(
          style = "color: #e65100; font-weight: bold;",
          sprintf("\u2728 Novel: %d", n_novel)
        )
      )
    )
  })

  # Annotation table
  output$site_annotation_table <- DT::renderDT({
    req(values$phospho_annotations)
    ann <- values$phospho_annotations

    # Round numeric columns
    if ("logFC" %in% names(ann))     ann$logFC     <- round(ann$logFC, 3)
    if ("adj.P.Val" %in% names(ann)) ann$adj.P.Val <- signif(ann$adj.P.Val, 3)

    DT::datatable(
      ann,
      rownames = FALSE,
      escape = FALSE,  # Allow HTML in Evidence column
      filter = "top",
      selection = "multiple",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(which(names(ann) == "Status") - 1, "asc")),
        columnDefs = list(
          list(
            targets = which(names(ann) == "Status") - 1,
            render = DT::JS(
              "function(data, type, row) {",
              "  if (type === 'display') {",
              "    if (data === 'Known') {",
              "      return '<span style=\"background:#c8e6c9; color:#2e7d32; padding:2px 8px; border-radius:4px; font-weight:600;\">' + data + '</span>';",
              "    } else {",
              "      return '<span style=\"background:#ffe0b2; color:#e65100; padding:2px 8px; border-radius:4px; font-weight:600;\">' + data + '</span>';",
              "    }",
              "  }",
              "  return data;",
              "}"
            )
          )
        )
      )
    )
  })

  # Download annotation table
  output$download_annotation_table <- downloadHandler(
    filename = function() {
      paste0("phosphosite_annotations_", format(Sys.time(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(values$phospho_annotations)
      # Strip HTML from Evidence column for CSV export
      ann <- values$phospho_annotations
      ann$Evidence <- gsub("<[^>]+>", "", ann$Evidence)
      write.csv(ann, file, row.names = FALSE)
    }
  )

  # ============================================================================
  #  Info Modal Popups (? buttons)
  # ============================================================================

  observeEvent(input$phospho_volcano_info_btn, {
    showModal(modalDialog(
      title = "Phosphosite Volcano Plot",
      tags$p(
        "The volcano plot displays differential phosphosite abundance between experimental ",
        "conditions. Each point represents a single phosphorylation site. The x-axis shows ",
        "log2 fold change and the y-axis shows -log10 adjusted p-value (BH-corrected)."
      ),
      tags$p(
        tags$strong("Significance thresholds: "),
        "Sites passing both |log2FC| > 1 (2-fold change) and FDR < 0.05 are colored red. ",
        "Sites passing only FDR < 0.05 are colored blue. Top 15 significant sites are labeled."
      ),
      tags$p(
        "When protein-level abundance correction is enabled, log2FC values reflect ",
        "phosphorylation stoichiometry changes (site logFC minus protein logFC), ",
        "isolating regulation from protein expression changes."
      ),
      tags$hr(),
      tags$p(class = "text-muted small",
        tags$strong("References:"), tags$br(),
        "Ritchie ME et al. (2015) limma powers differential expression analyses for RNA-sequencing ",
        "and microarray studies. ", tags$em("Nucleic Acids Res"), " 43(7):e47.", tags$br(),
        "Demichev V et al. (2020) DIA-NN: neural networks and interference correction enable deep ",
        "proteome coverage in high throughput. ", tags$em("Nature Methods"), " 17:41-44.", tags$br(),
        "Benjamini Y & Hochberg Y (1995) Controlling the false discovery rate. ",
        tags$em("J Royal Stat Soc B"), " 57:289-300."
      ),
      easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  observeEvent(input$phospho_table_info_btn, {
    showModal(modalDialog(
      title = "Phosphosite Results Table",
      tags$p(
        "This table contains site-level differential expression statistics for every ",
        "quantified phosphorylation site. Each row represents a unique phosphosite ",
        "(Protein + Residue + Position)."
      ),
      tags$p(tags$strong("Key columns:")),
      tags$ul(
        tags$li(tags$strong("logFC"), " \u2014 log2 fold change between conditions"),
        tags$li(tags$strong("adj.P.Val"), " \u2014 Benjamini-Hochberg adjusted p-value (FDR)"),
        tags$li(tags$strong("Residue"), " \u2014 phosphorylated amino acid (S, T, or Y)"),
        tags$li(tags$strong("Position"), " \u2014 residue position (peptide-relative for Path B, ",
                "protein-relative for Path A / site matrix)")
      ),
      tags$p(
        "The table is filterable and sortable. Use the column filters to search for ",
        "specific genes or sites of interest."
      ),
      tags$hr(),
      tags$p(class = "text-muted small",
        tags$strong("References:"), tags$br(),
        "Ritchie ME et al. (2015) limma powers differential expression analyses. ",
        tags$em("Nucleic Acids Res"), " 43(7):e47.", tags$br(),
        "Tyanova S et al. (2016) The Perseus computational platform for comprehensive ",
        "analysis of (prote)omics data. ", tags$em("Nature Methods"), " 13:731-740."
      ),
      easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  observeEvent(input$phospho_residue_info_btn, {
    showModal(modalDialog(
      title = "Residue Distribution",
      tags$p(
        "This bar chart shows the distribution of phosphorylated residues (Serine, Threonine, ",
        "Tyrosine) across all quantified sites and among statistically significant sites."
      ),
      tags$p(
        tags$strong("Expected distribution "), "(from TiO2 or IMAC enrichment): ",
        "~85% phosphoserine (pS), ~14% phosphothreonine (pT), ~1% phosphotyrosine (pY). ",
        "Substantial deviations may indicate enrichment bias or biology-driven shifts."
      ),
      tags$p(
        "A higher-than-expected pY proportion among significant hits could suggest ",
        "active tyrosine kinase signaling in your experimental condition."
      ),
      tags$hr(),
      tags$p(class = "text-muted small",
        tags$strong("References:"), tags$br(),
        "Olsen JV et al. (2006) Global, in vivo, and site-specific phosphorylation dynamics ",
        "in signaling networks. ", tags$em("Cell"), " 127:635-648.", tags$br(),
        "Sharma K et al. (2014) Ultradeep human phosphoproteome reveals a distinct regulatory ",
        "nature of Tyr and Ser/Thr-based signaling. ", tags$em("Cell Reports"), " 8:1583-1594.", tags$br(),
        "Humphrey SJ et al. (2015) High-throughput phosphoproteomics reveals in vivo insulin ",
        "signaling dynamics. ", tags$em("Nature Biotechnology"), " 33:990-995."
      ),
      easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  observeEvent(input$phospho_completeness_info_btn, {
    showModal(modalDialog(
      title = "Phosphosite Completeness QC",
      tags$p(
        "This histogram shows the data completeness for each phosphosite across all samples. ",
        "A site with 100% completeness was quantified in every sample; a site with 50% ",
        "completeness was quantified in half of the samples."
      ),
      tags$p(
        tags$strong("Why this matters: "),
        "Phosphoproteomics data is inherently more missing than protein-level data because ",
        "phosphopeptides are less abundant. Sites with very low completeness may represent ",
        "stochastic detection events rather than reproducible biology."
      ),
      tags$p(
        "The red dashed line at 50% marks the default filtering threshold. Sites below this ",
        "threshold are removed before imputation and differential expression testing to reduce ",
        "false discoveries from imputation-driven signals."
      ),
      tags$hr(),
      tags$p(class = "text-muted small",
        tags$strong("References:"), tags$br(),
        "Webb-Robertson BJ et al. (2015) Review, evaluation, and discussion of the challenges ",
        "of missing value imputation for mass spectrometry-based label-free global proteomics. ",
        tags$em("J Proteome Res"), " 14:1993-2001.", tags$br(),
        "Lazar C et al. (2016) Accounting for the multiple natures of missing values in label-free ",
        "quantitative proteomics data sets to compare imputation strategies. ",
        tags$em("J Proteome Res"), " 15:1116-1125."
      ),
      easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  observeEvent(input$ksea_info_btn, {
    showModal(modalDialog(
      title = "Kinase-Substrate Enrichment Analysis (KSEA)",
      tags$p(
        "KSEA infers upstream kinase activity from the phosphosite fold-changes in your data. ",
        "It uses known kinase-substrate relationships from the ",
        tags$strong("PhosphoSitePlus"), " and ", tags$strong("NetworKIN"), " databases to ",
        "aggregate substrate fold-changes for each kinase into a z-score."
      ),
      tags$p(tags$strong("How to interpret:")),
      tags$ul(
        tags$li(tags$strong("Positive z-score (red)"), " \u2014 substrates of this kinase are ",
                "collectively up-regulated, suggesting the kinase is ", tags$em("more active"), "."),
        tags$li(tags$strong("Negative z-score (blue)"), " \u2014 substrates are collectively ",
                "down-regulated, suggesting the kinase is ", tags$em("less active"), "."),
        tags$li(tags$strong("n="), " \u2014 number of substrates from your data matched to that kinase. ",
                "Higher n gives more statistical power.")
      ),
      tags$p(
        "FDR < 0.05 kinases are shown by default. The bar labels indicate the number of ",
        "substrate sites contributing to each kinase score."
      ),
      tags$hr(),
      tags$p(class = "text-muted small",
        tags$strong("References:"), tags$br(),
        "Casado P et al. (2013) Kinase-substrate enrichment analysis provides insights into ",
        "the heterogeneity of signaling pathway activation in leukemia cells. ",
        tags$em("Sci Signaling"), " 6:rs6.", tags$br(),
        "Wiredja DD et al. (2017) The KSEA App: a web-based tool for kinase activity inference ",
        "from quantitative phosphoproteomics. ", tags$em("Bioinformatics"), " 33:3489-3491.", tags$br(),
        "Hornbeck PV et al. (2015) PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. ",
        tags$em("Nucleic Acids Res"), " 43:D512-D520.", tags$br(),
        "Linding R et al. (2007) Systematic discovery of in vivo phosphorylation networks. ",
        tags$em("Cell"), " 129:1415-1426."
      ),
      easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  observeEvent(input$phospho_motif_info_btn, {
    showModal(modalDialog(
      title = "Phosphosite Motif Analysis",
      tags$p(
        "Sequence logos display the amino acid composition around regulated phosphosites ",
        "(\u00b17 residues flanking the phosphorylated position). Logos are shown separately ",
        "for up-regulated and down-regulated sites (FDR < 0.05, |log2FC| > 1)."
      ),
      tags$p(tags$strong("How to interpret motifs:")),
      tags$ul(
        tags$li(tags$strong("Proline at +1"), " \u2014 proline-directed kinases (CDKs, MAPKs, GSK3)"),
        tags$li(tags$strong("Acidic residues (D/E) at +1 to +3"), " \u2014 CK2 family"),
        tags$li(tags$strong("Arginine at -3"), " \u2014 basophilic kinases (PKA, PKC, CaMKII)"),
        tags$li(tags$strong("No clear motif"), " \u2014 mixed kinase regulation or insufficient data")
      ),
      tags$p(
        "A minimum of 10 sites is required per direction. Requires FASTA protein sequences ",
        "(uploaded via the sidebar) to extract flanking residues."
      ),
      tags$hr(),
      tags$p(class = "text-muted small",
        tags$strong("References:"), tags$br(),
        "Wagih O (2017) ggseqlogo: a versatile R package for drawing sequence logos. ",
        tags$em("Bioinformatics"), " 33:3645-3647.", tags$br(),
        "Schwartz D & Gygi SP (2005) An iterative statistical approach to the identification of ",
        "protein phosphorylation motifs from large-scale data sets. ",
        tags$em("Nature Biotechnology"), " 23:1391-1398.", tags$br(),
        "Miller ML et al. (2008) Linear motif atlas for phosphorylation-dependent signaling. ",
        tags$em("Sci Signaling"), " 1:ra2."
      ),
      easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  observeEvent(input$site_annotation_info_btn, {
    showModal(modalDialog(
      title = "Phosphosite Annotation: Known vs Novel",
      tags$p(
        "This tool queries two major phosphosite databases to classify each site in your ",
        "dataset as ", tags$strong("Known"), " (previously reported) or ",
        tags$strong("Novel"), " (not yet in databases)."
      ),
      tags$p(tags$strong("Data sources:")),
      tags$ul(
        tags$li(tags$strong("UniProt"), " \u2014 curated phosphorylation sites with experimental ",
                "evidence codes (ECO:0000269 = published experiment, ECO:0007744 = large-scale study)"),
        tags$li(tags$strong("PhosphoSitePlus"), " \u2014 kinase-substrate database bundled with the ",
                "KSEAapp R package (Hornbeck et al. 2015)")
      ),
      tags$p(
        "A site is marked 'Known' if it matches either database at the same protein + residue + position. ",
        "Novel sites may represent genuinely undiscovered regulation or positional mismatches ",
        "(e.g., peptide-relative vs protein-relative numbering)."
      ),
      tags$p(
        tags$em("Tip: "), "Path A (DIA-NN site matrix upload) provides protein-relative positions ",
        "that match database numbering more accurately than Path B (parsed from report)."
      ),
      tags$hr(),
      tags$p(class = "text-muted small",
        tags$strong("References:"), tags$br(),
        "UniProt Consortium (2023) UniProt: the Universal Protein Knowledgebase in 2023. ",
        tags$em("Nucleic Acids Res"), " 51:D523-D531.", tags$br(),
        "Hornbeck PV et al. (2015) PhosphoSitePlus, 2014: mutations, PTMs and recalibrations. ",
        tags$em("Nucleic Acids Res"), " 43:D512-D520.", tags$br(),
        "Ochoa D et al. (2020) The functional landscape of the human phosphoproteome. ",
        tags$em("Nature Biotechnology"), " 38:365-373."
      ),
      easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  # ============================================================================
  #  Fullscreen Modal Plots
  # ============================================================================

  # --- Phospho Volcano Fullscreen ---
  observeEvent(input$phospho_volcano_fullscreen_btn, {
    showModal(modalDialog(
      title = "Phosphosite Volcano Plot",
      plotOutput("phospho_volcano_fullscreen", height = "85vh"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  output$phospho_volcano_fullscreen <- renderPlot({
    req(values$phospho_fit, input$phospho_contrast_selector)
    de <- limma::topTable(values$phospho_fit,
                          coef = input$phospho_contrast_selector,
                          number = Inf, sort.by = "none")
    de$SiteID <- rownames(de)
    if (!is.null(values$phospho_site_info)) {
      de <- merge(de, values$phospho_site_info, by = "SiteID", all.x = TRUE)
    }
    correction_active <- FALSE
    if (isTRUE(input$phospho_protein_correction) && !is.null(values$fit)) {
      corrected <- phospho_corrected_de()
      if (!is.null(corrected)) {
        corr_map <- stats::setNames(corrected$corrected_logFC, corrected$SiteID)
        matched <- corr_map[de$SiteID]
        de$logFC[!is.na(matched)] <- matched[!is.na(matched)]
        correction_active <- TRUE
      }
    }
    de$sig <- dplyr::case_when(
      de$adj.P.Val < 0.05 & abs(de$logFC) > 1 ~ "Significant",
      de$adj.P.Val < 0.05 ~ "FDR < 0.05",
      TRUE ~ "NS"
    )
    n_sig  <- sum(de$sig == "Significant")
    n_up   <- sum(de$sig == "Significant" & de$logFC > 0)
    n_down <- sum(de$sig == "Significant" & de$logFC < 0)
    de$label <- ifelse(
      !is.na(de$Genes) & !is.na(de$Residue) & de$Genes != "",
      paste0(de$Genes, " ", de$Residue, de$Position), de$SiteID
    )
    sig_de <- de[de$sig == "Significant", ]
    if (nrow(sig_de) > 0) {
      sig_de <- sig_de[order(sig_de$adj.P.Val), ]
      label_de <- head(sig_de, 20)
    } else {
      label_de <- sig_de[0, ]
    }
    ggplot2::ggplot(de, ggplot2::aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
      ggplot2::geom_point(alpha = 0.6, size = 2) +
      ggplot2::scale_color_manual(values = c(
        "Significant" = "#E63946", "FDR < 0.05" = "#457B9D", "NS" = "gray70"
      )) +
      ggrepel::geom_text_repel(
        data = label_de, ggplot2::aes(label = label),
        size = 3.5, max.overlaps = 25, color = "black"
      ) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
      ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = paste("Phosphosite Volcano:", input$phospho_contrast_selector,
                      if (correction_active) "(protein-corrected)" else ""),
        subtitle = sprintf("%d sites | %d significant (\u2191%d \u2193%d) at |FC|>2 & FDR<0.05",
                           nrow(de), n_sig, n_up, n_down),
        x = if (correction_active) "Corrected log2 FC (phospho - protein)" else "log2 Fold Change (phosphosite)",
        y = "-log10(adjusted p-value)"
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(legend.position = "bottom",
                     plot.subtitle = ggplot2::element_text(color = "gray40"))
  })

  # --- Residue Distribution Fullscreen ---
  observeEvent(input$phospho_residue_fullscreen_btn, {
    showModal(modalDialog(
      title = "Residue Distribution",
      plotOutput("phospho_residue_fullscreen", height = "85vh"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  output$phospho_residue_fullscreen <- renderPlot({
    req(values$phospho_fit, input$phospho_contrast_selector, values$phospho_site_info)
    de <- limma::topTable(values$phospho_fit,
                          coef = input$phospho_contrast_selector,
                          number = Inf)
    de$SiteID <- rownames(de)
    de <- merge(de, values$phospho_site_info, by = "SiteID")
    all_counts <- table(factor(de$Residue, levels = c("S", "T", "Y")))
    sig_de <- de[de$adj.P.Val < 0.05, ]
    sig_counts <- table(factor(sig_de$Residue, levels = c("S", "T", "Y")))
    plot_data <- data.frame(
      Residue  = rep(c("S", "T", "Y"), 2),
      Count    = c(as.numeric(all_counts), as.numeric(sig_counts)),
      Category = rep(c("All quantified", "Significant (FDR < 0.05)"), each = 3)
    )
    plot_data$Count[is.na(plot_data$Count)] <- 0
    ggplot2::ggplot(plot_data, ggplot2::aes(x = Residue, y = Count, fill = Category)) +
      ggplot2::geom_col(position = "dodge", alpha = 0.8) +
      ggplot2::scale_fill_manual(values = c(
        "All quantified" = "#457B9D", "Significant (FDR < 0.05)" = "#E63946"
      )) +
      ggplot2::labs(
        title = "Phosphosite Residue Distribution",
        subtitle = "Expected: ~85% Ser / ~14% Thr / ~1% Tyr (typical TiO2/IMAC enrichment)",
        x = "Phosphorylated Residue", y = "Number of Sites"
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(legend.position = "bottom")
  })

  # --- Completeness QC Fullscreen ---
  observeEvent(input$phospho_completeness_fullscreen_btn, {
    showModal(modalDialog(
      title = "Phosphosite Completeness QC",
      plotOutput("phospho_completeness_fullscreen", height = "85vh"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  output$phospho_completeness_fullscreen <- renderPlot({
    req(values$phospho_site_matrix)
    mat <- values$phospho_site_matrix
    n_sites   <- nrow(mat)
    n_samples <- ncol(mat)
    site_completeness <- rowSums(!is.na(mat)) / n_samples * 100
    hist_data <- data.frame(Completeness = site_completeness)
    ggplot2::ggplot(hist_data, ggplot2::aes(x = Completeness)) +
      ggplot2::geom_histogram(bins = 20, fill = "#457B9D", color = "white", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = 50, linetype = "dashed", color = "red") +
      ggplot2::annotate("text", x = 52, y = Inf, vjust = 2, hjust = 0,
                        label = "50% threshold", color = "red", size = 4.5) +
      ggplot2::labs(
        title = "Phosphosite Quantification Completeness",
        subtitle = sprintf("%s sites across %d samples | Median completeness: %.0f%%",
                           format(n_sites, big.mark = ","), n_samples,
                           median(site_completeness)),
        x = "% of Samples with Quantification", y = "Number of Phosphosites"
      ) +
      ggplot2::theme_bw(base_size = 16)
  })

  # --- KSEA Barplot Fullscreen ---
  observeEvent(input$ksea_fullscreen_btn, {
    showModal(modalDialog(
      title = "Kinase Activity (KSEA)",
      plotOutput("ksea_barplot_fullscreen", height = "85vh"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  output$ksea_barplot_fullscreen <- renderPlot({
    req(values$ksea_results)
    ks <- values$ksea_results
    ks <- ks[order(ks$z.score, decreasing = TRUE), ]
    sig_ks <- ks[ks$FDR < 0.05, ]
    if (nrow(sig_ks) == 0) sig_ks <- utils::head(ks, 20)
    top_up   <- utils::head(sig_ks[sig_ks$z.score > 0, ], 15)
    top_down <- utils::tail(sig_ks[sig_ks$z.score < 0, ], 15)
    plot_ks  <- rbind(top_up, top_down)
    if (nrow(plot_ks) == 0) {
      plot.new()
      text(0.5, 0.5, "No kinases scored (insufficient substrate matches)")
      return()
    }
    plot_ks <- plot_ks[order(plot_ks$z.score), ]
    plot_ks$Kinase.Gene <- factor(plot_ks$Kinase.Gene, levels = plot_ks$Kinase.Gene)
    plot_ks$Direction <- ifelse(plot_ks$z.score > 0, "Activated", "Inhibited")
    ggplot2::ggplot(plot_ks, ggplot2::aes(x = z.score, y = Kinase.Gene, fill = Direction)) +
      ggplot2::geom_col(alpha = 0.85) +
      ggplot2::scale_fill_manual(values = c("Activated" = "#E63946", "Inhibited" = "#457B9D")) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("n=%d", m)),
        hjust = ifelse(plot_ks$z.score > 0, -0.1, 1.1), size = 3.5
      ) +
      ggplot2::labs(
        title = paste("Kinase Activity:", values$ksea_last_contrast),
        subtitle = "KSEA z-scores (PhosphoSitePlus + NetworKIN substrates)",
        x = "KSEA z-score", y = NULL
      ) +
      ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(legend.position = "bottom")
  })

  # --- Motif Analysis Fullscreen ---
  observeEvent(input$phospho_motif_fullscreen_btn, {
    showModal(modalDialog(
      title = "Motif Analysis",
      plotOutput("phospho_motif_up_fullscreen", height = "40vh"),
      plotOutput("phospho_motif_down_fullscreen", height = "40vh"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  output$phospho_motif_up_fullscreen <- renderPlot({
    seqs <- phospho_flanking_seqs()
    if (is.null(seqs) || length(seqs$up) < 10) {
      plot.new()
      text(0.5, 0.5, "Insufficient up-regulated sites for motif", cex = 1.2, col = "gray50")
      return()
    }
    if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
      plot.new()
      text(0.5, 0.5, "ggseqlogo package not installed", cex = 1.2, col = "gray50")
      return()
    }
    ggseqlogo::ggseqlogo(seqs$up, method = "bits", seq_type = "aa") +
      ggplot2::ggtitle(sprintf("Up-regulated phosphosites (n=%d)", length(seqs$up))) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 14))
  })

  output$phospho_motif_down_fullscreen <- renderPlot({
    seqs <- phospho_flanking_seqs()
    if (is.null(seqs) || length(seqs$down) < 10) {
      plot.new()
      text(0.5, 0.5, "Insufficient down-regulated sites for motif", cex = 1.2, col = "gray50")
      return()
    }
    if (!requireNamespace("ggseqlogo", quietly = TRUE)) {
      plot.new()
      text(0.5, 0.5, "ggseqlogo package not installed", cex = 1.2, col = "gray50")
      return()
    }
    ggseqlogo::ggseqlogo(seqs$down, method = "bits", seq_type = "aa") +
      ggplot2::ggtitle(sprintf("Down-regulated phosphosites (n=%d)", length(seqs$down))) +
      ggplot2::theme(plot.title = ggplot2::element_text(size = 14))
  })

}
