# ==============================================================================
#  SERVER MODULE -- De Novo Controls (Confidence Filtering, Manuscript Stats,
#                   GO Annotation, Disagreement Analysis)
#  Called from app.R as: server_denovo_controls(input, output, session, values)
#
#  IMPORTANT: This module uses its own filtered reactive to avoid conflicts with
#  server_dda.R (which owns the raw data and classification logic).
# ==============================================================================

server_denovo_controls <- function(input, output, session, values) {

  # ============================================================================
  #  Feature 1: Interactive Confidence Filtering
  # ============================================================================


  # Core filtered reactive — ALL downstream renders should use this

  filtered_casanovo_psms <- reactive({
    # Use unified reactive (works for both Casanovo and Cascadia)
    psms <- values$denovo_psms %||% values$dda_casanovo_psms
    req(psms)
    req(nrow(psms) > 0)

    threshold <- input$dda_denovo_score_threshold %||% 0.9
    psms[psms$score >= threshold, ]
  })


  # Filtered classification: re-classifies using the filtered PSMs
  filtered_classification <- reactive({
    psms <- filtered_casanovo_psms()
    req(nrow(psms) > 0)

    # For Cascadia mode, use the pre-computed classification from the adapter
    if (isTRUE(values$denovo_engine == "cascadia") && !is.null(values$denovo_classification)) {
      return(values$denovo_classification)
    }

    # For Casanovo mode, re-classify with the threshold-filtered PSMs
    sage_psms <- values$dda_sage_psms
    if (is.null(sage_psms) || nrow(sage_psms) == 0) {
      # No Sage data — everything is unclassified
      return(list(
        classified      = psms,
        confirmed       = psms[0, ],
        novel           = psms,
        protein_summary = data.frame(
          proteins = character(0),
          n_casanovo_confirmed = integer(0),
          casanovo_max_score = numeric(0),
          casanovo_mean_aa_score = numeric(0),
          stringsAsFactors = FALSE
        ),
        summary_stats = list(
          n_total = nrow(psms), n_confirmed = 0L, n_novel = nrow(psms),
          pct_confirmed = 0, pct_novel = 100
        )
      ))
    }

    tryCatch(
      classify_dda_denovo(psms, sage_psms),
      error = function(e) {
        message("[denovo_controls] Classification error: ", e$message)
        NULL
      }
    )
  })


  # Sync filtered classification to values so click handlers in server_dda.R
  # can reference the correct data for per-residue visualization
  observe({
    cls <- filtered_classification()
    values$dda_filtered_classification <- cls
  })


  # --- Threshold count display ---
  output$dda_denovo_threshold_count <- renderUI({
    req(values$dda_casanovo_psms)
    total_all <- nrow(values$dda_casanovo_psms)
    threshold <- input$dda_denovo_score_threshold %||% 0.9

    n_above <- sum(values$dda_casanovo_psms$score >= threshold, na.rm = TRUE)
    pct <- round(100 * n_above / max(total_all, 1), 1)

    cls <- filtered_classification()
    n_conf <- if (!is.null(cls)) cls$summary_stats$n_confirmed else 0
    n_novel <- if (!is.null(cls)) cls$summary_stats$n_novel else 0

    tags$div(
      style = "font-size: 14px; line-height: 1.6;",
      tags$span(
        style = "font-weight: 600; color: #2c3e50;",
        format(n_above, big.mark = ","), " / ", format(total_all, big.mark = ","),
        " PSMs above threshold (", pct, "%)"
      ),
      tags$br(),
      tags$span(style = "color: #2ecc71; font-weight: 500;",
        icon("check-circle"), " ", format(n_conf, big.mark = ","), " confirmed"),
      tags$span(style = "margin-left: 12px; color: #e67e22; font-weight: 500;",
        icon("question-circle"), " ", format(n_novel, big.mark = ","), " novel")
    )
  })


  # --- Override summary cards to use filtered data ---
  output$dda_denovo_summary_cards <- renderUI({
    cls <- filtered_classification()
    req(cls)

    n_total     <- cls$summary_stats$n_total
    n_confirmed <- cls$summary_stats$n_confirmed
    n_novel     <- cls$summary_stats$n_novel
    pct_conf    <- cls$summary_stats$pct_confirmed
    n_proteins  <- if (!is.null(cls$protein_summary)) nrow(cls$protein_summary) else 0

    # BLAST stats (use full blast data, filtered to peptides above threshold)
    blast <- values$denovo_blast %||% values$dda_casanovo_blast
    n_blast <- 0L
    if (!is.null(blast) && nrow(blast) > 0) {
      novel_seqs <- cls$novel$seq_stripped
      n_blast <- length(unique(blast$peptide[blast$peptide %in% novel_seqs]))
    }

    tags$div(
      class = "row",
      style = "margin-bottom: 15px;",
      tags$div(class = "col-md-2",
        tags$div(class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #3498db; padding: 15px;",
          tags$h4(format(n_total, big.mark = ","), style = "margin: 0; color: #3498db;"),
          tags$small("Above Threshold")
        )
      ),
      tags$div(class = "col-md-2",
        tags$div(class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #2ecc71; padding: 15px;",
          tags$h4(format(n_confirmed, big.mark = ","), style = "margin: 0; color: #2ecc71;"),
          tags$small("Confirmed (in Sage)")
        )
      ),
      tags$div(class = "col-md-2",
        tags$div(class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #e67e22; padding: 15px;",
          tags$h4(format(n_novel, big.mark = ","), style = "margin: 0; color: #e67e22;"),
          tags$small("Novel Peptides")
        )
      ),
      tags$div(class = "col-md-2",
        tags$div(class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #9b59b6; padding: 15px;",
          tags$h4(paste0(pct_conf, "%"), style = "margin: 0; color: #9b59b6;"),
          tags$small(paste0("Confirm Rate (", n_proteins, " prot)"))
        )
      ),
      tags$div(class = "col-md-2",
        tags$div(class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #1abc9c; padding: 15px;",
          tags$h4(format(n_blast, big.mark = ","), style = "margin: 0; color: #1abc9c;"),
          tags$small("BLAST Hits")
        )
      ),
      tags$div(class = "col-md-2",
        tags$div(class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #e74c3c; padding: 15px;",
          tags$h4(
            sprintf("%.2f", input$dda_denovo_score_threshold %||% 0.9),
            style = "margin: 0; color: #e74c3c;"
          ),
          tags$small("Score Cutoff")
        )
      )
    )
  })


  # --- Override confirmed table to use filtered data ---
  output$dda_denovo_confirmed_table <- DT::renderDT({
    cls <- filtered_classification()
    req(cls)
    confirmed <- cls$confirmed
    req(nrow(confirmed) > 0)

    display_df <- data.frame(
      Sequence    = confirmed$sequence,
      Stripped    = confirmed$seq_stripped,
      Score       = round(confirmed$score, 3),
      Charge      = confirmed$charge,
      AA_Scores   = if ("mean_aa_score" %in% names(confirmed)) {
        round(confirmed$mean_aa_score, 3)
      } else {
        NA_real_
      },
      Protein     = if ("proteins" %in% names(confirmed)) {
        confirmed$proteins
      } else {
        NA_character_
      },
      Source_File  = confirmed$source_file,
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
        order      = list(list(2, "desc")),
        dom        = "Bfrtip",
        buttons    = list("csv", "excel")
      ),
      extensions = "Buttons",
      caption = htmltools::tags$caption(
        style = "caption-side: top; font-weight: bold; color: #2ecc71;",
        paste0("Confirmed peptides (score >= ",
               input$dda_denovo_score_threshold %||% 0.9, ")")
      )
    )
  })


  # --- Override novel table to use filtered data ---
  output$dda_denovo_novel_table <- DT::renderDT({
    cls <- filtered_classification()
    req(cls)
    novel <- cls$novel
    req(nrow(novel) > 0)

    display_df <- data.frame(
      Sequence    = novel$sequence,
      Stripped    = novel$seq_stripped,
      Score       = round(novel$score, 3),
      Charge      = novel$charge,
      AA_Scores   = if ("mean_aa_score" %in% names(novel)) {
        round(novel$mean_aa_score, 3)
      } else {
        NA_real_
      },
      Source_File  = novel$source_file,
      stringsAsFactors = FALSE
    )

    # Append DIAMOND BLAST results if available
    blast <- values$denovo_blast %||% values$dda_casanovo_blast
    if (!is.null(blast) && nrow(blast) > 0) {
      blast_dedup <- blast[!duplicated(blast$peptide), ]
      blast_map    <- stats::setNames(blast_dedup$subject, blast_dedup$peptide)
      identity_map <- stats::setNames(blast_dedup$pident, blast_dedup$peptide)
      evalue_map   <- stats::setNames(blast_dedup$evalue, blast_dedup$peptide)

      display_df$BLAST_Hit    <- blast_map[novel$seq_stripped]
      display_df$Identity_Pct <- round(identity_map[novel$seq_stripped], 1)
      display_df$E_Value      <- evalue_map[novel$seq_stripped]
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
        paste0("Novel peptides (score >= ",
               input$dda_denovo_score_threshold %||% 0.9, ")")
      )
    )
  })


  # --- Override score distribution to use filtered threshold line ---
  output$dda_denovo_score_dist <- plotly::renderPlotly({
    req(values$dda_casanovo_psms)
    req(nrow(values$dda_casanovo_psms) > 0)

    # Use ALL PSMs for histogram, but show threshold line
    df <- values$dda_casanovo_psms
    threshold <- input$dda_denovo_score_threshold %||% 0.9

    cls <- filtered_classification()
    confirmed_seqs <- if (!is.null(cls)) cls$confirmed$seq_norm else character(0)

    match_type <- ifelse(df$seq_norm %in% confirmed_seqs & df$score >= threshold,
      "Confirmed",
      ifelse(df$score >= threshold, "Novel (above)", "Below threshold")
    )

    plot_df <- data.frame(
      score = df$score,
      type  = match_type,
      stringsAsFactors = FALSE
    )

    colors <- c("Confirmed" = "#2ecc71", "Novel (above)" = "#e67e22",
                "Below threshold" = "#bdc3c7")

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = score, fill = type)) +
      ggplot2::geom_histogram(bins = 50, alpha = 0.8, position = "stack") +
      ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                          color = "#e74c3c", linewidth = 1) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(
        x = "Casanovo Confidence Score",
        y = "Count",
        fill = "Classification",
        subtitle = paste0("Threshold: ", threshold,
                          " | ", sum(df$score >= threshold), " / ",
                          nrow(df), " PSMs above")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")

    plotly::ggplotly(p) %>%
      plotly::layout(
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.05)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })


  # ============================================================================
  #  Feature 2: Manuscript Summary Statistics (Table 1)
  # ============================================================================

  manuscript_data <- reactive({
    psms <- filtered_casanovo_psms()
    req(nrow(psms) > 0)

    cls <- filtered_classification()
    sage_psms <- values$dda_sage_psms
    blast <- values$denovo_blast %||% values$dda_casanovo_blast

    # Get per-source_file breakdown
    source_files <- unique(psms$source_file)

    rows <- lapply(source_files, function(sf) {
      sf_psms <- psms[psms$source_file == sf, ]
      n_above <- nrow(sf_psms)

      # Total Sage spectra for this file
      n_sage_spectra <- 0L
      if (!is.null(sage_psms)) {
        # Sage filename column is 'filename'
        sage_fn_norm <- gsub("\\.(d|mzML|mgf)$", "", basename(sage_psms$filename))
        n_sage_spectra <- sum(sage_fn_norm == sf, na.rm = TRUE)
      }

      # Classification counts
      n_confirmed <- 0L
      n_novel <- 0L
      if (!is.null(cls)) {
        conf_seqs <- if (nrow(cls$confirmed) > 0) {
          cls$confirmed$seq_norm[cls$confirmed$source_file == sf]
        } else character(0)
        n_confirmed <- length(conf_seqs)
        n_novel <- sum(cls$novel$source_file == sf, na.rm = TRUE)
      }

      # BLAST stats for this file's novel peptides
      n_blast_hits <- 0L
      n_unique_proteins <- 0L
      pct_contaminant <- 0
      if (!is.null(blast) && nrow(blast) > 0 && !is.null(cls) && nrow(cls$novel) > 0) {
        novel_seqs <- cls$novel$seq_stripped[cls$novel$source_file == sf]
        blast_hits <- blast[blast$peptide %in% novel_seqs, ]
        n_blast_hits <- length(unique(blast_hits$peptide))

        if (nrow(blast_hits) > 0) {
          # Extract accessions
          accessions <- unique(blast_hits$subject)
          n_unique_proteins <- length(accessions)

          # Contaminant check
          if ("is_contaminant" %in% names(blast_hits)) {
            n_contam <- length(unique(blast_hits$subject[blast_hits$is_contaminant]))
            pct_contaminant <- round(100 * n_contam / max(n_unique_proteins, 1), 1)
          }
        }
      }

      # Median confidence
      median_score <- round(median(sf_psms$score, na.rm = TRUE), 3)

      data.frame(
        Sample             = sf,
        Sage_Spectra       = n_sage_spectra,
        Casanovo_PSMs      = n_above,
        Confirmed          = n_confirmed,
        Novel              = n_novel,
        BLAST_Hits         = n_blast_hits,
        Unique_Proteins    = n_unique_proteins,
        Contaminant_Pct    = pct_contaminant,
        Median_Score       = median_score,
        stringsAsFactors   = FALSE
      )
    })

    result <- do.call(rbind, rows)

    # Totals row
    totals <- data.frame(
      Sample             = "TOTAL",
      Sage_Spectra       = sum(result$Sage_Spectra),
      Casanovo_PSMs      = sum(result$Casanovo_PSMs),
      Confirmed          = sum(result$Confirmed),
      Novel              = sum(result$Novel),
      BLAST_Hits         = sum(result$BLAST_Hits),
      Unique_Proteins    = length(unique(
        if (!is.null(blast) && nrow(blast) > 0) blast$subject else character(0)
      )),
      Contaminant_Pct    = {
        if (!is.null(blast) && nrow(blast) > 0 && "is_contaminant" %in% names(blast)) {
          round(100 * length(unique(blast$subject[blast$is_contaminant])) /
                  max(length(unique(blast$subject)), 1), 1)
        } else 0
      },
      Median_Score       = round(median(psms$score, na.rm = TRUE), 3),
      stringsAsFactors   = FALSE
    )

    rbind(result, totals)
  })

  output$dda_manuscript_summary <- DT::renderDT({
    df <- manuscript_data()
    req(nrow(df) > 0)

    DT::datatable(
      df,
      rownames = FALSE,
      selection = "none",
      options = list(
        pageLength = 50,
        paging = FALSE,
        searching = FALSE,
        scrollX = TRUE,
        dom = "t",
        columnDefs = list(
          list(className = "dt-right", targets = 1:8)
        )
      ),
      caption = htmltools::tags$caption(
        style = "caption-side: top; font-weight: bold; color: #2d6a2d;",
        paste0("Per-sample summary (score >= ",
               input$dda_denovo_score_threshold %||% 0.9, ")")
      )
    ) %>%
      DT::formatStyle(
        "Sample",
        target = "row",
        fontWeight = DT::styleEqual("TOTAL", "bold"),
        backgroundColor = DT::styleEqual("TOTAL", "#f0f0f0")
      )
  })

  output$dda_denovo_manuscript_csv <- downloadHandler(
    filename = function() {
      paste0("denovo_manuscript_summary_", Sys.Date(), ".csv")
    },
    content = function(file) {
      df <- manuscript_data()
      utils::write.csv(df, file, row.names = FALSE)
    }
  )


  # ============================================================================
  #  Feature 3: GO/Functional Annotation (protein name-based classification)
  # ============================================================================

  # Classify proteins into functional categories from SwissProt names
  functional_categories <- reactive({
    blast <- values$denovo_blast %||% values$dda_casanovo_blast
    req(blast, nrow(blast) > 0)

    cls <- filtered_classification()
    req(cls)

    # Get novel peptides at current threshold
    novel_seqs <- cls$novel$seq_stripped

    # Filter BLAST to novel peptides from filtered set
    blast_filt <- blast[blast$peptide %in% novel_seqs, ]
    if (nrow(blast_filt) == 0) return(NULL)

    # Parse protein name from SwissProt format: sp|ACC|PROTNAME_SPECIES -> PROTNAME
    blast_filt$protein_name <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", blast_filt$subject))

    # Classification rules (case-insensitive on raw subject for broader matching)
    classify_protein <- function(prot_name, subject) {
      pn <- toupper(prot_name)
      subj <- toupper(subject)

      # Structural / keratinization
      if (grepl("^KRT|^K1C|^K2C|^K22|^KR[0-9]|KERA|COLL|^CO[0-9]|^COL[0-9]|ELAS|FIBR|^FBN|LAMI|^DSP|^PKP|PLAK|DESM|CORN|LORI|^SPRR|^IVL|^FLG|FILGG|^DSC|^DSG|DESMC|PLEC", pn))
        return("Structural/Keratinization")

      # Histones / nuclear
      if (grepl("^H[1234]|^HIST|HISTON|^H2A|^H2B|^H3|^H4|NUCLE|LAMIN|^NUP", pn))
        return("Nuclear/Chromatin")

      # Common lab contaminants
      if (grepl("TRYP|ALBU.*BOVIN|^CAS[12A]|TRFE.*BOVIN|LACB|OVAL.*CHICK", pn))
        return("Lab Contaminant")

      # Metabolic enzymes
      if (grepl("DEHYDR|OXIDA|REDUC|KINASE|PHOSPH|SYNTH|LYASE|LIGASE|CARBOX|ENOLA|ALDOL|GAPDH|^PGK|^PKM|^ENO|ISOMER", pn))
        return("Metabolic/Enzymatic")

      # Immune / defense
      if (grepl("IMMUN|ANTI|COMPL|INTERL|^IG[HKLM]|DEFEN|LYSO|CATHL", pn))
        return("Immune/Defense")

      # Cytoskeletal
      if (grepl("ACTI|TUBU|MYOS|VIMEN|TROPO|GELS|COFI|PROFI|^MYH|^MYL|^ACT[ABCG]|^TBB|^TBA", pn))
        return("Cytoskeletal")

      # Ribosomal / translation
      if (grepl("^RS[0-9]|^RL[0-9]|RIBOS|ELONGA|^EF1|^EIF|^RPS|^RPL", pn))
        return("Ribosomal/Translation")

      # Signaling
      if (grepl("RECEPT|SIGNAL|NOTCH|WNT|HEDGEHOG|MAPK|^RAS|GTPASE|GPCR", pn))
        return("Signaling")

      # Heat shock / chaperone
      if (grepl("HSP|HEAT|CHAPER|^GRP|CALRET|CALNEX|^CCT|^TCP", pn))
        return("Chaperone/Stress")

      # Hemoglobin / blood
      if (grepl("^HB[ABGDEZ]|HEMOG|GLOBIN|^HBA|^HBB|FERRI|TRANSF", pn))
        return("Hemoglobin/Blood")

      "Other"
    }

    blast_filt$category <- mapply(classify_protein,
      blast_filt$protein_name, blast_filt$subject,
      USE.NAMES = FALSE
    )

    blast_filt
  })


  output$dda_denovo_go_summary <- renderUI({
    fc <- functional_categories()
    if (is.null(fc) || nrow(fc) == 0) {
      return(tags$div(
        style = "padding: 20px; text-align: center; color: #6c757d;",
        icon("info-circle"),
        " Run DIAMOND BLAST first to see functional annotation."
      ))
    }

    # Summary
    cats <- table(fc$category)
    n_cats <- length(cats)
    top_cat <- names(sort(cats, decreasing = TRUE))[1]
    n_top <- max(cats)

    # Paleoproteomics highlight
    n_kerat <- sum(grepl("Keratinization", names(cats)))
    paleo_note <- if (n_kerat > 0) {
      tags$div(
        style = "background: #fff8e1; border-left: 4px solid #ffa000; padding: 10px; border-radius: 4px; margin-top: 8px;",
        icon("feather", style = "color: #ffa000;"),
        tags$span(style = "font-weight: 600;",
          " Paleoproteomics note: "),
        paste0("Structural/Keratinization proteins detected (",
               cats["Structural/Keratinization"], " BLAST hits). ",
               "Alpha-keratins and feather beta-keratins are key markers in ancient specimens.")
      )
    } else NULL

    tags$div(
      tags$div(
        style = "background: #f0f7ff; border-left: 4px solid #1565c0; padding: 12px; border-radius: 4px; margin-bottom: 12px;",
        tags$span(style = "font-weight: 600; color: #1565c0;",
          icon("layer-group"), " Functional Classification"),
        tags$br(),
        sprintf("BLAST hits classified into %d functional categories. ", n_cats),
        sprintf("Most abundant: %s (%d hits).", top_cat, n_top),
        tags$br(),
        tags$small(style = "color: #666;",
          "Categories assigned from SwissProt protein names using keyword rules. ",
          "Not a formal GO annotation — for exploratory/manuscript context.")
      ),
      paleo_note
    )
  })


  output$dda_denovo_go_bar <- plotly::renderPlotly({
    fc <- functional_categories()
    req(fc, nrow(fc) > 0)

    cat_counts <- as.data.frame(table(fc$category), stringsAsFactors = FALSE)
    names(cat_counts) <- c("Category", "Count")
    cat_counts <- cat_counts[order(cat_counts$Count, decreasing = TRUE), ]

    # Color palette (keratinization highlighted)
    cat_colors <- c(
      "Structural/Keratinization" = "#e65100",
      "Nuclear/Chromatin"         = "#1565c0",
      "Lab Contaminant"           = "#c62828",
      "Metabolic/Enzymatic"       = "#2e7d32",
      "Immune/Defense"            = "#6a1b9a",
      "Cytoskeletal"              = "#00838f",
      "Ribosomal/Translation"     = "#ef6c00",
      "Signaling"                 = "#ad1457",
      "Chaperone/Stress"          = "#4e342e",
      "Hemoglobin/Blood"          = "#b71c1c",
      "Other"                     = "#78909c"
    )

    # Map colors to actual categories
    bar_colors <- vapply(cat_counts$Category, function(cat) {
      if (cat %in% names(cat_colors)) cat_colors[[cat]] else "#78909c"
    }, character(1))

    # Also compute unique proteins per category
    prot_counts <- vapply(cat_counts$Category, function(cat) {
      length(unique(fc$subject[fc$category == cat]))
    }, integer(1))

    cat_counts$Unique_Proteins <- prot_counts

    p <- plotly::plot_ly(
      data = cat_counts,
      x = ~reorder(Category, Count),
      y = ~Count,
      type = "bar",
      marker = list(color = bar_colors),
      text = ~paste0(Category, "\n", Count, " BLAST hits\n",
                     Unique_Proteins, " unique proteins"),
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "BLAST Hit Count"),
        title = list(
          text = "Functional Category Distribution (BLAST Hits)",
          font = list(size = 14)
        ),
        margin = list(b = 120)
      )

    p %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })


  output$dda_denovo_go_table <- DT::renderDT({
    fc <- functional_categories()
    req(fc, nrow(fc) > 0)

    # Build per-category summary table
    cats <- unique(fc$category)
    summary_rows <- lapply(cats, function(cat) {
      sub <- fc[fc$category == cat, ]
      data.frame(
        Category        = cat,
        BLAST_Hits      = nrow(sub),
        Unique_Proteins = length(unique(sub$subject)),
        Unique_Peptides = length(unique(sub$peptide)),
        Median_Identity = round(median(sub$pident, na.rm = TRUE), 1),
        Top_Protein     = {
          prot_tab <- sort(table(sub$subject), decreasing = TRUE)
          if (length(prot_tab) > 0) {
            nm <- names(prot_tab)[1]
            sub("^[a-z]+\\|[^|]+\\|", "", nm)  # strip sp|ACC| prefix
          } else ""
        },
        stringsAsFactors = FALSE
      )
    })
    summary_df <- do.call(rbind, summary_rows)
    summary_df <- summary_df[order(summary_df$BLAST_Hits, decreasing = TRUE), ]

    DT::datatable(
      summary_df,
      rownames = FALSE,
      selection = "none",
      options = list(
        pageLength = 20,
        paging = FALSE,
        scrollX = TRUE,
        dom = "t",
        order = list(list(1, "desc"))
      ),
      caption = htmltools::tags$caption(
        style = "caption-side: top; font-weight: bold; color: #1565c0;",
        "Functional category summary"
      )
    )
  })


  # ============================================================================
  #  Feature 4: Sage vs Casanovo Disagreement Analysis
  # ============================================================================

  disagreement_data <- reactive({
    cas_psms <- filtered_casanovo_psms()
    sage_psms <- values$dda_sage_psms
    req(cas_psms, sage_psms)
    req(nrow(cas_psms) > 0, nrow(sage_psms) > 0)

    # --- Build matching keys ---
    # Casanovo: psm_id (scan number) + source_file
    # Sage: scannr + filename (basename without extension)
    #
    # Sage column names from results.sage.tsv:
    #   scannr, filename, peptide, charge, hyperscore, spectrum_q, protein_q, proteins, ...
    # Casanovo mzTab columns (parsed): psm_id, source_file, sequence, score, charge, ...

    # Normalize sage filenames for matching
    sage_fn <- gsub("\\.(d|mzML|mgf)$", "", basename(sage_psms$filename))
    sage_key <- paste0(sage_fn, ":", sage_psms$scannr)

    cas_key <- paste0(cas_psms$source_file, ":", cas_psms$psm_id)

    # Find spectra present in BOTH tools
    shared_keys <- intersect(cas_key, sage_key)
    if (length(shared_keys) == 0) {
      # Try matching by psm_id only (single-file case)
      sage_key <- as.character(sage_psms$scannr)
      cas_key <- as.character(cas_psms$psm_id)
      shared_keys <- intersect(cas_key, sage_key)
    }

    if (length(shared_keys) == 0) return(NULL)

    # Merge
    cas_idx <- match(shared_keys, cas_key)
    sage_idx <- match(shared_keys, sage_key)

    merged <- data.frame(
      scan_key       = shared_keys,
      casanovo_seq   = cas_psms$seq_stripped[cas_idx],
      sage_seq       = sage_psms$peptide[sage_idx],
      casanovo_score = round(cas_psms$score[cas_idx], 3),
      sage_score     = round(sage_psms$hyperscore[sage_idx], 2),
      sage_q         = formatC(sage_psms$spectrum_q[sage_idx], format = "e", digits = 2),
      charge         = cas_psms$charge[cas_idx],
      sage_protein   = sage_psms$proteins[sage_idx],
      source_file    = cas_psms$source_file[cas_idx],
      stringsAsFactors = FALSE
    )

    # I/L normalize both for comparison
    cas_norm <- gsub("I", "L", toupper(merged$casanovo_seq))
    sage_norm <- gsub("I", "L", toupper(merged$sage_seq))

    # Classify disagreement type
    merged$agreement <- mapply(function(c_seq, s_seq) {
      if (c_seq == s_seq) return("Exact match")

      # Check I/L swap only
      c_il <- gsub("I", "L", c_seq)
      s_il <- gsub("I", "L", s_seq)
      if (c_il == s_il) return("I/L swap (benign)")

      # Same length?
      if (nchar(c_il) == nchar(s_il)) {
        # Count mismatches
        c_chars <- strsplit(c_il, "")[[1]]
        s_chars <- strsplit(s_il, "")[[1]]
        n_diff <- sum(c_chars != s_chars)
        if (n_diff == 1) return("Single AA substitution")
        if (n_diff == 2) return("Two AA substitutions")
        return("Multiple substitutions")
      }

      # Different lengths
      "Completely different"
    }, toupper(merged$casanovo_seq), toupper(merged$sage_seq),
    USE.NAMES = FALSE)

    # Only keep disagreements (exclude exact + I/L swap)
    disagreements <- merged[!merged$agreement %in% c("Exact match", "I/L swap (benign)"), ]

    # Add the matched spectra for context
    attr(disagreements, "all_matched") <- merged
    attr(disagreements, "n_total_matched") <- nrow(merged)
    attr(disagreements, "n_exact") <- sum(merged$agreement == "Exact match")
    attr(disagreements, "n_il_swap") <- sum(merged$agreement == "I/L swap (benign)")

    disagreements
  })


  output$dda_denovo_disagree_summary <- renderUI({
    # Disagreement analysis only works in DDA mode (Sage vs Casanovo)
    if (isTRUE(values$denovo_engine == "cascadia")) {
      return(tags$div(
        style = "padding: 20px; text-align: center; color: #6c757d;",
        icon("info-circle"),
        " Disagreement analysis compares database search vs de novo sequencing on the same spectra. ",
        "This requires DDA mode (Sage + Casanovo). Not available for DIA/Cascadia data."
      ))
    }

    cas_psms <- filtered_casanovo_psms()
    sage_psms <- values$dda_sage_psms

    if (is.null(sage_psms) || nrow(sage_psms) == 0) {
      return(tags$div(
        style = "padding: 20px; text-align: center; color: #6c757d;",
        icon("info-circle"),
        " Requires both Sage and Casanovo results to compare."
      ))
    }

    disagree <- disagreement_data()
    if (is.null(disagree)) {
      return(tags$div(
        style = "padding: 20px; text-align: center; color: #6c757d;",
        icon("info-circle"),
        " No matched spectra found between Sage and Casanovo. ",
        "Check that scan numbers and filenames align."
      ))
    }

    n_total    <- attr(disagree, "n_total_matched")
    n_exact    <- attr(disagree, "n_exact")
    n_il       <- attr(disagree, "n_il_swap")
    n_disagree <- nrow(disagree)
    pct_agree  <- round(100 * (n_exact + n_il) / max(n_total, 1), 1)

    # Breakdown of disagreement types
    type_counts <- if (n_disagree > 0) {
      table(disagree$agreement)
    } else {
      integer(0)
    }

    tags$div(
      tags$div(
        style = "background: #f0f7ff; border-left: 4px solid #1565c0; padding: 12px; border-radius: 4px; margin-bottom: 12px;",
        tags$span(style = "font-weight: 600; color: #1565c0;",
          icon("code-compare"), " Sage vs Casanovo Spectrum-Level Comparison"),
        tags$br(),
        sprintf("%s spectra matched between tools. ", format(n_total, big.mark = ",")),
        tags$br(),
        tags$span(style = "color: #2ecc71; font-weight: 500;",
          icon("check"), sprintf(" %s exact matches", format(n_exact, big.mark = ","))),
        tags$span(style = "margin-left: 12px; color: #3498db;",
          sprintf("+ %s I/L swaps (benign)", format(n_il, big.mark = ","))),
        tags$span(style = "margin-left: 12px; color: #e74c3c; font-weight: 500;",
          sprintf("= %s disagreements", format(n_disagree, big.mark = ","))),
        tags$br(),
        tags$span(style = "font-weight: 600;",
          sprintf("Agreement rate: %s%%", pct_agree))
      ),

      if (n_disagree > 0) {
        tags$div(
          style = "display: flex; gap: 12px; margin-bottom: 12px; flex-wrap: wrap;",
          lapply(names(type_counts), function(typ) {
            color <- switch(typ,
              "Single AA substitution" = "#e67e22",
              "Two AA substitutions"   = "#d35400",
              "Multiple substitutions" = "#c0392b",
              "Completely different"   = "#e74c3c",
              "#95a5a6"
            )
            tags$div(
              style = paste0("background: #f8f9fa; border-left: 4px solid ", color,
                             "; padding: 8px 12px; border-radius: 4px; min-width: 160px;"),
              tags$span(style = paste0("font-weight: 600; color: ", color, ";"),
                type_counts[typ]),
              tags$br(),
              tags$small(typ)
            )
          })
        )
      } else {
        tags$div(
          style = "background: #e8f5e9; padding: 12px; border-radius: 4px; text-align: center;",
          icon("check-circle", style = "color: #2e7d32; font-size: 18px;"),
          tags$span(style = "color: #2e7d32; font-weight: 600;",
            " Perfect agreement: all matched spectra have consistent sequences.")
        )
      }
    )
  })


  output$dda_denovo_disagree_table <- DT::renderDT({
    disagree <- disagreement_data()
    req(disagree, nrow(disagree) > 0)

    display_df <- data.frame(
      Scan            = disagree$scan_key,
      Casanovo_Seq    = disagree$casanovo_seq,
      Sage_Seq        = disagree$sage_seq,
      Type            = disagree$agreement,
      Casanovo_Score  = disagree$casanovo_score,
      Sage_Hyperscore = disagree$sage_score,
      Sage_Q          = disagree$sage_q,
      Charge          = disagree$charge,
      Sage_Protein    = disagree$sage_protein,
      File            = disagree$source_file,
      stringsAsFactors = FALSE
    )

    DT::datatable(
      display_df,
      rownames = FALSE,
      filter = "top",
      selection = "none",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(4, "desc")),  # Sort by Casanovo score desc
        dom = "Bfrtip",
        buttons = list("csv", "excel")
      ),
      extensions = "Buttons",
      caption = htmltools::tags$caption(
        style = "caption-side: top; font-weight: bold; color: #e74c3c;",
        paste0("Spectra where Sage and Casanovo assigned different sequences ",
               "(", nrow(display_df), " disagreements)")
      )
    ) %>%
      DT::formatStyle(
        "Type",
        backgroundColor = DT::styleEqual(
          c("Single AA substitution", "Two AA substitutions",
            "Multiple substitutions", "Completely different"),
          c("#fff3e0", "#ffe0b2", "#ffccbc", "#ffcdd2")
        )
      )
  })

  # --- Source engine badge ---
  output$denovo_source_badge <- renderUI({
    engine <- values$denovo_engine
    if (is.null(engine)) return(NULL)

    label <- if (engine == "casanovo") "Casanovo (DDA)" else "Cascadia (DIA)"
    ref <- if (engine == "casanovo") "Sage" else "DIA-NN"
    color <- if (engine == "casanovo") "#6a1b9a" else "#1565c0"
    ico <- if (engine == "casanovo") "wand-magic-sparkles" else "dna"

    div(style = "margin-bottom: 8px;",
      tags$span(
        style = paste0("display: inline-block; padding: 4px 12px; border-radius: 12px; ",
                       "background: ", color, "; color: white; font-size: 12px; font-weight: 600;"),
        icon(ico), paste0(" De novo: ", label, " vs ", ref, " database search")
      )
    )
  })

  # ============================================================================
  #  INFO MODALS — De Novo Controls Sub-tabs
  # ============================================================================

  observeEvent(input$denovo_score_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " Score Distribution & QC"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("Quality metrics for de novo peptide predictions from Casanovo."),
        tags$h6("Score Distribution"),
        tags$ul(
          tags$li(strong("Confidence Score (0-1): "), "Casanovo's overall peptide-level score. ",
            "Products of per-residue softmax probabilities. Higher = more confident."),
          tags$li(strong("Score Threshold: "), "Use the slider above to filter peptides by minimum confidence. ",
            "Default 0.9 retains high-quality predictions; lower to 0.7-0.8 for exploratory analysis.")
        ),
        tags$h6("Peptide Length Distribution"),
        p("Most tryptic peptides are 7-25 amino acids. Very short (<6 AA) or very long (>30 AA) peptides ",
          "may indicate non-specific cleavage or sequencing artifacts."),
        tags$h6("Charge State Distribution"),
        p("Expected: mostly 2+ and 3+ for tryptic peptides. A high proportion of 4+ or 5+ may indicate ",
          "incomplete digestion or unusual peptide properties."),
        tags$hr(),
        tags$p(style = "color: #666; font-size: 0.85em;",
          icon("camera"), " Click the camera icon on any plot to download as SVG for publication figures.")
      )
    ))
  })

  observeEvent(input$denovo_go_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " GO/Functional Annotation"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("Functional annotation of BLAST-matched proteins based on protein name patterns."),
        tags$h6("Categories"),
        tags$ul(
          tags$li(strong("Keratin: "), "Hair, feather, and skin structural proteins (keratins, corneous proteins)."),
          tags$li(strong("Collagen: "), "Connective tissue structural proteins."),
          tags$li(strong("Histone: "), "Chromatin packaging proteins — highly conserved across species."),
          tags$li(strong("Hemoglobin: "), "Oxygen transport proteins — indicates blood contamination in tissue samples."),
          tags$li(strong("Ribosomal: "), "Translation machinery — ubiquitous, often contaminants."),
          tags$li(strong("Metabolic: "), "Enzymes in metabolic pathways."),
          tags$li(strong("Cytoskeletal: "), "Actin, tubulin, intermediate filaments.")
        ),
        p("This is a name-based heuristic, not a formal Gene Ontology analysis. ",
          "For full GO enrichment, use the GSEA tab after database search.")
      )
    ))
  })

  observeEvent(input$denovo_disagree_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " Sage vs Casanovo Disagreements"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("Spectra where Sage (database search) and Casanovo (de novo) assigned different sequences."),
        tags$h6("Disagreement Types"),
        tags$ul(
          tags$li(strong("Single AA substitution: "), "One amino acid differs. Common for I/L ambiguity ",
            "(isoleucine and leucine are isobaric — indistinguishable by mass)."),
          tags$li(strong("Two AA substitutions: "), "Two positions differ. May indicate a genuine variant."),
          tags$li(strong("Multiple substitutions: "), "Three or more differences. Review per-residue scores."),
          tags$li(strong("Completely different: "), "Entirely different sequences. One tool is likely wrong.")
        ),
        tags$h6("Interpretation"),
        p("Disagreements are expected — they highlight where de novo and database approaches diverge. ",
          "Use per-residue confidence scores to judge which assignment is more reliable. ",
          "I/L substitutions (highlighted amber) are mass-equivalent and should be ignored."),
        tags$p(style = "color: #888; font-size: 0.85em;",
          "Note: Disagreement analysis requires both Sage AND Casanovo results for the same raw files.")
      )
    ))
  })

}
