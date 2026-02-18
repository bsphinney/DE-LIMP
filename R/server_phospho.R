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

    tags$div(
      class = "alert alert-info py-2 px-3 mb-3",
      style = "border-left: 4px solid #0d6efd;",
      tags$div(
        style = "display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px;",
        tags$span(
          icon("flask"), tags$strong(" Phosphoproteomics Data Detected — "),
          enrichment_msg,
          sprintf(" %s phosphoprecursors across %s total.",
                  format(pd$n_phospho, big.mark = ","),
                  format(pd$n_total, big.mark = ","))
        ),
        actionButton("go_to_phospho", "Open Phospho Tab \u2192",
                      class = "btn btn-sm btn-outline-primary")
      )
    )
  })

  # Navigate to phospho tab on banner button click
  observeEvent(input$go_to_phospho, {
    nav_select("main_tabs", "Phosphoproteomics")
  })

  # ============================================================================
  #  Site matrix upload (Path A: DIA-NN 1.9+ site matrix)
  # ============================================================================
  observeEvent(input$phospho_site_matrix_file, {
    req(input$phospho_site_matrix_file)
    tryCatch({
      mat_df <- arrow::read_parquet(input$phospho_site_matrix_file$datapath,
                                     as_data_frame = TRUE)

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
          # Use base-matched samples
          mat_idx  <- match(common_base, mat_base)
          meta_idx <- match(common_base, meta_base)
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
          tags$div(style = "text-align: right; margin-bottom: 10px;",
            downloadButton("download_phospho_volcano", "Save Plot",
                           class = "btn-outline-secondary btn-sm")
          ),
          plotOutput("phospho_volcano", height = "550px")
        ),

        nav_panel("Site Table",
          tags$div(style = "text-align: right; margin-bottom: 10px;",
            downloadButton("download_phospho_table", "Export CSV",
                           class = "btn-success btn-sm")
          ),
          DT::DTOutput("phospho_site_table")
        ),

        nav_panel("Residue Distribution",
          plotOutput("phospho_residue_dist", height = "450px")
        ),

        nav_panel("QC: Completeness",
          plotOutput("phospho_completeness", height = "450px")
        ),

        nav_panel("Kinase Activity",
          tags$div(
            style = "display: flex; align-items: center; gap: 15px; flex-wrap: wrap; margin-bottom: 15px;",
            actionButton("run_ksea", "Run Kinase Enrichment (KSEA)",
                         class = "btn-info", icon = icon("project-diagram")),
            tags$span(id = "ksea_status_text", class = "text-muted small")
          ),
          tags$p(class = "text-muted small",
            "Infers upstream kinase activity from phosphosite fold-changes using the ",
            "PhosphoSitePlus + NetworKIN kinase-substrate database. ",
            "Requires the ", tags$code("KSEAapp"), " R package. ",
            "Reference: Casado et al. 2013, Sci Signal 6:rs6; Wiredja et al. 2017."
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
          tags$p(class = "text-muted small",
            "Sequence logos show the amino acid enrichment around regulated phosphosites. ",
            "Distinct motifs suggest specific kinase families are active. For example, ",
            "proline at +1 indicates proline-directed kinases (CDKs, MAPKs); ",
            "acidic residues (D/E) at +1 to +3 suggest CK2. ",
            "Requires FASTA upload (sidebar) and the ", tags$code("ggseqlogo"), " R package."
          ),
          tags$p(class = "text-muted small",
            "Reference: Wagih O (2017) ggseqlogo, Bioinformatics 33:3645-3647."
          ),
          hr(),
          plotOutput("phospho_motif_up", height = "220px"),
          plotOutput("phospho_motif_down", height = "220px")
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
          NetworKIN.cutoff = 5,
          m.cutoff     = 2,
          p.cutoff     = 0.05
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
          "ksea_scores <- KSEAapp::KSEA.Scores(KSData, ksea_input, NetworKIN=TRUE, NetworKIN.cutoff=5, m.cutoff=2)"
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

}
