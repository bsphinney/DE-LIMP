# ==============================================================================
#  SERVER MODULE — Grid View, Signal Distribution, Dataset Summary
#  Called from app.R as: server_viz(input, output, session, values, add_to_log, is_hf_space)
# ==============================================================================

server_viz <- function(input, output, session, values, add_to_log, is_hf_space) {

  # --- STATUS MESSAGE (line 835) ---
  output$run_status_msg <- renderText({ values$status })

  # --- DATASET SUMMARY (lines 838-956) ---
  # Dataset summary as tab content (instead of modal)
  output$dataset_summary_content <- renderUI({
    req(values$metadata)

    summary_elements <- list()

    # File summary
    summary_elements[[length(summary_elements) + 1]] <- div(
      style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
      tags$h4(icon("file"), " File Summary"),
      tags$hr(),
      tags$p(style = "font-size: 1.1em;",
        icon("folder-open"), " ",
        strong("Total Files: "), nrow(values$metadata)
      ),
      tags$p(style = "font-size: 1.1em;",
        icon("users"), " ",
        strong("Assigned Groups: "),
        length(unique(values$metadata$Group[values$metadata$Group != ""]))
      )
    )

    # Dataset metrics (if pipeline has run)
    if (!is.null(values$y_protein)) {
      avg_signal <- rowMeans(values$y_protein$E, na.rm = TRUE)
      min_linear <- 2^min(avg_signal, na.rm = TRUE)
      max_linear <- 2^max(avg_signal, na.rm = TRUE)

      dynamic_range_text <- if (min_linear > 1e-10) {
        orders_of_magnitude <- log10(max_linear / min_linear)
        paste(round(orders_of_magnitude, 1), "orders of magnitude")
      } else {
        "N/A (Min signal is zero)"
      }

      summary_elements[[length(summary_elements) + 1]] <- div(
        style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
        tags$h4(icon("chart-bar"), " Dataset Metrics"),
        tags$hr(),
        tags$p(style = "font-size: 1.1em;",
          icon("signal"), " ",
          strong("Signal Dynamic Range: "), dynamic_range_text
        ),
        tags$p(style = "font-size: 1.1em;",
          icon("dna"), " ",
          strong("Total Proteins Quantified: "), nrow(values$y_protein$E)
        )
      )
    }

    # Differential expression summary (if DE analysis has run)
    if (!is.null(values$fit)) {
      # Calculate DE proteins per comparison
      all_comparisons <- colnames(values$fit$contrasts)
      de_summary_list <- list()

      for (comp in all_comparisons) {
        de_results <- topTable(values$fit, coef = comp, number = Inf)
        n_sig <- sum(de_results$adj.P.Val < 0.05, na.rm = TRUE)
        n_up <- sum(de_results$adj.P.Val < 0.05 & de_results$logFC > 0, na.rm = TRUE)
        n_down <- sum(de_results$adj.P.Val < 0.05 & de_results$logFC < 0, na.rm = TRUE)

        # Parse comparison name (format: "GroupA - GroupB")
        comp_parts <- strsplit(comp, " - ")[[1]]
        if (length(comp_parts) == 2) {
          group_a <- trimws(comp_parts[1])
          group_b <- trimws(comp_parts[2])
        } else {
          group_a <- "Group 1"
          group_b <- "Group 2"
        }

        # Create explicit, easy-to-read summary
        de_summary_list[[length(de_summary_list) + 1]] <- div(
          style = "margin-bottom: 15px; padding: 12px; background-color: white; border-left: 4px solid #667eea; border-radius: 4px;",
          # Comparison header
          tags$p(style = "font-size: 1.05em; margin-bottom: 8px; font-weight: 500;",
            icon("microscope"), " ",
            strong(comp),
            span(style = "margin-left: 10px; color: #6c757d; font-weight: normal; font-size: 0.9em;",
              paste0("(", n_sig, " significant proteins)")
            )
          ),
          # Detailed breakdown
          if (n_sig > 0) {
            tagList(
              tags$div(style = "margin-left: 25px; font-size: 0.95em;",
                tags$div(style = "margin-bottom: 4px;",
                  span(style = "color: #e41a1c; font-weight: 500;", "\u2191 ", n_up),
                  " proteins higher in ",
                  strong(style = "color: #e41a1c;", group_a)
                ),
                tags$div(
                  span(style = "color: #377eb8; font-weight: 500;", "\u2193 ", n_down),
                  " proteins higher in ",
                  strong(style = "color: #377eb8;", group_b)
                )
              )
            )
          } else {
            tags$div(style = "margin-left: 25px; font-size: 0.9em; color: #6c757d; font-style: italic;",
              "No significant differences detected"
            )
          }
        )
      }

      summary_elements[[length(summary_elements) + 1]] <- div(
        style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;",
        tags$h4(icon("flask"), " Differential Expression Summary"),
        tags$hr(),
        tags$p(style = "font-size: 0.85em; color: #6c757d; margin-bottom: 15px;",
          "Proteins with FDR-adjusted p-value < 0.05. Arrows indicate direction of change."
        ),
        do.call(tagList, de_summary_list)
      )

      # Complete dataset export button
      summary_elements[[length(summary_elements) + 1]] <- div(
        style = "background-color: #eef2ff; padding: 20px; border-radius: 8px; border: 1px solid #c7d2fe;",
        div(style = "display: flex; justify-content: space-between; align-items: center;",
          div(
            tags$h4(icon("download"), " Export Complete Dataset", style = "margin: 0 0 5px 0;"),
            tags$p(style = "font-size: 0.85em; color: #6c757d; margin: 0;",
              "DE statistics for all comparisons, expression values, gene symbols, and sample metadata — ready for downstream tools."
            )
          ),
          downloadButton("download_complete_dataset", "Download Complete Dataset",
            class = "btn-primary", style = "white-space: nowrap;")
        )
      )
    }

    tagList(summary_elements)
  })

  # --- COMPLETE DATASET EXPORT ---
  output$download_complete_dataset <- downloadHandler(
    filename = function() {
      paste0("Limpa_Complete_Dataset_", format(Sys.Date(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(values$fit, values$y_protein, values$metadata)

      withProgress(message = "Building complete dataset...", value = 0, {

        # --- 1. Gene symbol mapping ---
        incProgress(0.1, detail = "Mapping gene symbols...")
        protein_ids <- rownames(values$y_protein$E)
        accessions <- str_split_fixed(protein_ids, "[; ]", 2)[,1]
        org_db_name <- detect_organism_db(protein_ids)

        id_map <- tryCatch({
          if (!requireNamespace(org_db_name, quietly = TRUE)) BiocManager::install(org_db_name, ask = FALSE)
          library(org_db_name, character.only = TRUE)
          db_obj <- get(org_db_name)
          AnnotationDbi::select(db_obj, keys = accessions, columns = c("SYMBOL"), keytype = "UNIPROT") %>%
            dplyr::rename(Accession = UNIPROT, Gene = SYMBOL) %>% distinct(Accession, .keep_all = TRUE)
        }, error = function(e) data.frame(Accession = accessions, Gene = accessions))

        gene_df <- data.frame(Protein.Group = protein_ids, Accession = accessions, stringsAsFactors = FALSE)
        gene_df <- left_join(gene_df, id_map, by = "Accession")
        gene_df$Gene[is.na(gene_df$Gene)] <- gene_df$Accession[is.na(gene_df$Gene)]

        # --- 2. DE stats for ALL contrasts ---
        incProgress(0.3, detail = "Gathering DE statistics for all comparisons...")
        all_contrasts <- colnames(values$fit$contrasts)
        de_combined <- data.frame(Protein.Group = protein_ids, stringsAsFactors = FALSE)

        for (cname in all_contrasts) {
          tt <- topTable(values$fit, coef = cname, number = Inf) %>% as.data.frame()
          if (!"Protein.Group" %in% colnames(tt)) tt <- tt %>% rownames_to_column("Protein.Group")
          safe_name <- make.names(cname)
          tt_subset <- tt %>% dplyr::select(Protein.Group, logFC, P.Value, adj.P.Val)
          colnames(tt_subset) <- c("Protein.Group",
            paste0("logFC_", safe_name),
            paste0("P.Value_", safe_name),
            paste0("adj.P.Val_", safe_name))
          de_combined <- left_join(de_combined, tt_subset, by = "Protein.Group")
        }

        # --- 3. Expression matrix ---
        incProgress(0.6, detail = "Adding expression values...")
        exprs_df <- as.data.frame(values$y_protein$E) %>% rownames_to_column("Protein.Group")

        # --- 4. Combine everything ---
        incProgress(0.8, detail = "Assembling and writing file...")
        full_export <- gene_df %>%
          dplyr::select(Protein.Group, Accession, Gene) %>%
          left_join(de_combined, by = "Protein.Group") %>%
          left_join(exprs_df, by = "Protein.Group")

        # --- 5. Add sample metadata as header rows ---
        # Build metadata annotation lines that map to sample columns
        sample_cols <- colnames(values$y_protein$E)
        meta <- values$metadata
        non_sample_cols <- setdiff(colnames(full_export), sample_cols)

        # Get custom covariate names
        cov1_name <- if (!is.null(values$cov1_name) && nzchar(values$cov1_name)) values$cov1_name else "Covariate1"
        cov2_name <- if (!is.null(values$cov2_name) && nzchar(values$cov2_name)) values$cov2_name else "Covariate2"

        # Build annotation rows
        annot_rows <- list()
        annot_fields <- list(
          Group = "Group", Batch = "Batch",
          Covariate1 = cov1_name, Covariate2 = cov2_name
        )
        for (col_name in names(annot_fields)) {
          if (col_name %in% colnames(meta) && any(nzchar(meta[[col_name]]))) {
            label <- annot_fields[[col_name]]
            vals <- meta[[col_name]][match(sample_cols, meta$File.Name)]
            vals[is.na(vals)] <- ""
            row <- c(paste0("#", label), rep("", length(non_sample_cols) - 1), vals)
            annot_rows[[length(annot_rows) + 1]] <- row
          }
        }

        # Write: annotation header rows, then column names, then data
        con <- file(file, "w")
        for (arow in annot_rows) {
          writeLines(paste(arow, collapse = ","), con)
        }
        close(con)
        write.table(full_export, file, sep = ",", row.names = FALSE, quote = TRUE, append = TRUE)

        incProgress(1.0, detail = "Done!")
      })
    }
  )

  # --- GRID VIEW & PLOT LOGIC (lines 1064-1191) ---
  grid_react_df <- reactive({
    req(values$y_protein, values$metadata)

    if (!is.null(values$fit) && !is.null(input$contrast_selector) && nzchar(input$contrast_selector)) {
      # Full DE mode — use topTable for gene names and significance
      df_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(df_raw)) {
        df_raw <- df_raw %>% rownames_to_column("Protein.Group")
      }
    } else {
      # No-DE mode (no replicates) — build from y_protein directly
      df_raw <- data.frame(
        Protein.Group = rownames(values$y_protein$E),
        logFC = 0, P.Value = 1, adj.P.Val = 1,
        stringsAsFactors = FALSE
      )
      # Add gene info from y_protein$genes if available
      if (!is.null(values$y_protein$genes)) {
        gene_cols <- intersect(c("Genes", "Protein.Names"), colnames(values$y_protein$genes))
        if (length(gene_cols) > 0) {
          df_raw <- cbind(df_raw, values$y_protein$genes[, gene_cols, drop = FALSE])
        }
      }
    }

    df_raw$Accession <- str_split_fixed(df_raw$Protein.Group, "[; ]", 2)[, 1]

    # Clean contaminant prefixes (Cont_P04264 → P04264 for UniProt lookup)
    df_raw$Accession_clean <- gsub("^Cont_", "", df_raw$Accession)
    df_raw$is_contaminant <- grepl("^Cont_", df_raw$Accession)

    # Detect if accessions are NCBI RefSeq (XP_, NP_, WP_) vs UniProt
    non_contam <- df_raw$Accession_clean[!df_raw$is_contaminant]
    is_ncbi <- length(non_contam) > 0 && any(grepl("^[XNW]P_", head(non_contam, 50)))

    if (is_ncbi) {
      # NCBI accessions — try gene_map.tsv locally, then via SSH
      ncbi_gene_map <- NULL

      # Search common local locations
      search_dirs <- c(tempdir(), "/data/fasta", "/quobyte/proteomics-grp/de-limp/fasta")
      if (!is.null(values$diann_fasta_files)) {
        search_dirs <- c(dirname(values$diann_fasta_files), search_dirs)
      }
      for (d in unique(search_dirs)) {
        if (!dir.exists(d)) next
        gmaps <- list.files(d, pattern = "gene_map\\.tsv$", full.names = TRUE)
        if (length(gmaps) > 0) {
          ncbi_gene_map <- tryCatch(read.delim(gmaps[1], stringsAsFactors = FALSE), error = function(e) NULL)
          if (!is.null(ncbi_gene_map) && nrow(ncbi_gene_map) > 0) {
            message("[Grid] Loaded gene map: ", gmaps[1], " (", nrow(ncbi_gene_map), " entries)")
            break
          }
        }
      }

      # SSH fallback: download gene map from HIVE if not found locally
      if (is.null(ncbi_gene_map) && isTRUE(values$ssh_connected)) {
        tryCatch({
          cfg <- list(host = isolate(input$ssh_host), user = isolate(input$ssh_user),
                      port = isolate(input$ssh_port) %||% 22L,
                      key_path = isolate(input$ssh_key_path))
          # Find gene_map.tsv on remote
          remote_result <- ssh_exec(cfg,
            "ls /quobyte/proteomics-grp/de-limp/fasta/*gene_map.tsv 2>/dev/null | head -1",
            timeout = 10)
          if (remote_result$status == 0 && length(remote_result$stdout) > 0 &&
              nzchar(trimws(remote_result$stdout[1]))) {
            remote_path <- trimws(remote_result$stdout[1])
            local_path <- file.path(tempdir(), basename(remote_path))
            dl <- scp_download(cfg, remote_path, local_path)
            if (dl$status == 0 && file.exists(local_path)) {
              ncbi_gene_map <- tryCatch(read.delim(local_path, stringsAsFactors = FALSE), error = function(e) NULL)
              message("[Grid] Downloaded gene map via SSH: ", nrow(ncbi_gene_map), " entries")
            }
          }
        }, error = function(e) message("[Grid] SSH gene map download failed: ", e$message))
      }

      if (!is.null(ncbi_gene_map) && nrow(ncbi_gene_map) > 0 && "gene_symbol" %in% colnames(ncbi_gene_map)) {
        # Use gene map from batch E-utilities lookup
        ncbi_gene_map <- ncbi_gene_map[!duplicated(ncbi_gene_map$accession), ]
        df_raw <- df_raw %>%
          left_join(ncbi_gene_map, by = c("Accession" = "accession")) %>%
          mutate(
            Gene = ifelse(!is.na(gene_symbol) & nzchar(gene_symbol), gene_symbol, Accession),
            Protein.Name = ifelse(!is.na(protein_name) & nzchar(protein_name), protein_name, Protein.Group)
          )
      } else {
        # Fallback: parse from y_protein$genes (DIA-NN FASTA header parsing)
        id_map <- NULL
        genes_df <- values$y_protein$genes
        if (!is.null(genes_df)) {
          pn_col <- intersect(c("Protein.Names", "Protein.Name"), colnames(genes_df))
          gn_col <- intersect(c("Genes", "Gene.Names", "Gene"), colnames(genes_df))

        if (length(pn_col) > 0 || length(gn_col) > 0) {
          pg_col <- if ("Protein.Group" %in% colnames(genes_df)) "Protein.Group" else NULL
          fasta_map <- data.frame(
            Protein.Group = if (!is.null(pg_col)) genes_df[[pg_col]] else rownames(genes_df),
            stringsAsFactors = FALSE
          )
          if (length(gn_col) > 0) {
            fasta_map$SYMBOL <- genes_df[[gn_col[1]]]
          }
          if (length(pn_col) > 0) {
            # Clean NCBI protein names: remove "[Organism name]" suffix
            raw_names <- genes_df[[pn_col[1]]]
            fasta_map$GENENAME <- gsub("\\s*\\[.*\\]\\s*$", "", raw_names)
          }
          fasta_map <- fasta_map %>% distinct(Protein.Group, .keep_all = TRUE)

          df_raw <- df_raw %>%
            left_join(fasta_map, by = "Protein.Group") %>%
            mutate(
              Gene = ifelse(!is.null(SYMBOL) & !is.na(SYMBOL) & nzchar(SYMBOL),
                            SYMBOL, Accession),
              Protein.Name = ifelse(!is.null(GENENAME) & !is.na(GENENAME) & nzchar(GENENAME),
                                    GENENAME, Protein.Group)
            )
        } else {
          df_raw$Gene <- df_raw$Accession
          df_raw$Protein.Name <- df_raw$Protein.Group
        }
      } else {
        df_raw$Gene <- df_raw$Accession
        df_raw$Protein.Name <- df_raw$Protein.Group
      }
      }  # end fallback else (no gene map TSV)
    } else {
      # UniProt accessions — use AnnotationDbi::select for gene symbol mapping
      org_db_name <- detect_organism_db(df_raw$Protein.Group)
      id_map <- tryCatch({
        if (!requireNamespace(org_db_name, quietly = TRUE)) {
          tryCatch(BiocManager::install(org_db_name, ask = FALSE, update = FALSE), error = function(e) NULL)
        }
        if (requireNamespace(org_db_name, quietly = TRUE)) {
          library(org_db_name, character.only = TRUE)
          db <- get(org_db_name)
          suppressMessages(AnnotationDbi::select(db, keys = unique(df_raw$Accession),
            keytype = "UNIPROT", columns = c("SYMBOL", "GENENAME")))
        } else NULL
      }, error = function(e) NULL)

      if (!is.null(id_map) && nrow(id_map) > 0) {
        colnames(id_map)[colnames(id_map) == "UNIPROT"] <- "Accession"
        id_map <- id_map %>% distinct(Accession, .keep_all = TRUE)
        df_raw <- df_raw %>%
          left_join(id_map, by = "Accession") %>%
          mutate(Gene = ifelse(is.na(SYMBOL), Accession, SYMBOL),
                 Protein.Name = ifelse(is.na(GENENAME), Protein.Group, GENENAME))
      } else {
        df_raw$Gene <- df_raw$Accession
        df_raw$Protein.Name <- df_raw$Protein.Group
      }
    }

    df_raw$Significance <- "Not Sig"
    if ("adj.P.Val" %in% colnames(df_raw)) {
      df_raw$Significance[df_raw$adj.P.Val < 0.05 & abs(df_raw$logFC) > (input$logfc_cutoff %||% 0.6)] <- "Significant"
    }

    df_volc <- df_raw
    df_exprs <- as.data.frame(values$y_protein$E) %>% rownames_to_column("Protein.Group")
    df_merged <- left_join(df_volc, df_exprs, by="Protein.Group")
    if (!is.null(values$plot_selected_proteins)) { df_merged <- df_merged %>% filter(Protein.Group %in% values$plot_selected_proteins) }

    meta_sorted <- values$metadata %>% arrange(Group, File.Name)
    ordered_files <- meta_sorted$File.Name
    valid_cols <- intersect(ordered_files, colnames(df_exprs))
    run_ids <- values$metadata$ID[match(valid_cols, values$metadata$File.Name)]
    new_headers <- as.character(run_ids)

    # Add Type column before final select
    df_merged$Type <- ifelse(grepl("^Cont_", df_merged$Protein.Group), "Contaminant", "Sample")

    df_final <- df_merged %>%
      dplyr::select(Protein.Group, Gene, Protein.Name, Significance, logFC, P.Value, adj.P.Val, Type, all_of(valid_cols)) %>%
      mutate(across(where(is.numeric), ~round(., 2)))

    df_final$Original.ID <- df_final$Protein.Group
    fixed_cols <- c("Protein.Group", "Gene", "Protein.Name", "Significance", "logFC", "P.Value", "adj.P.Val", "Type")
    colnames(df_final) <- c(fixed_cols, new_headers, "Original.ID")

    unique_groups <- sort(unique(meta_sorted$Group))
    group_colors <- setNames(rainbow(length(unique_groups), v=0.85, s=0.8), unique_groups)

    list(data = df_final, fixed_cols = fixed_cols, expr_cols = new_headers, valid_cols_map = valid_cols, meta_sorted = meta_sorted, group_colors = group_colors)
  })

  # Grid View legend UI
  output$grid_legend_ui <- renderUI({
    gdata <- grid_react_df()
    tags$div(
      style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
      tags$strong(icon("palette"), " Condition Legend:"), tags$br(),
      lapply(names(gdata$group_colors), function(grp) {
        tags$span(
          style = paste0("background-color:", gdata$group_colors[[grp]], "; color:white; padding:4px 10px; margin-right:8px; border-radius:4px; display:inline-block; margin-top:5px;"),
          grp
        )
      })
    )
  })

  # Grid View file mapping UI
  output$grid_file_map_ui <- renderUI({
    gdata <- grid_react_df()
    tags$details(
      style = "margin-bottom: 10px;",
      tags$summary(
        style = "cursor: pointer; color: #0d6efd;",
        icon("list"), " Click to view File ID Mapping (Run # \u2192 Filename)"
      ),
      tags$div(
        style = "max-height: 200px; overflow-y: auto; background: #f9f9f9; padding: 10px; border: 1px solid #dee2e6; border-radius: 4px; margin-top: 8px;",
        lapply(1:nrow(gdata$meta_sorted), function(i) {
          row <- gdata$meta_sorted[i, ]
          tags$div(
            style = "padding: 2px 0;",
            tags$span(style = "font-weight:bold; color:#007bff;", paste0("[", row$ID, "] ")),
            tags$span(row$File.Name),
            tags$span(style = "color:#6c757d; font-size:0.9em;", paste0(" (", row$Group, ")"))
          )
        })
      )
    )
  })

  observeEvent(input$grid_reset_selection, { values$plot_selected_proteins <- NULL })

  output$grid_view_table <- renderDT({
    gdata <- grid_react_df(); df_display <- gdata$data
    acc <- str_split_fixed(df_display$Original.ID, "[; ]", 2)[,1]
    # Link to correct database: NCBI for XP_/NP_/WP_, UniProt for others (strip Cont_ prefix)
    link_urls <- ifelse(grepl("^[XNW]P_", acc),
      paste0("https://www.ncbi.nlm.nih.gov/protein/", acc),
      paste0("https://www.uniprot.org/uniprotkb/", gsub("^Cont_", "", acc), "/entry")
    )
    df_display$Protein.Group <- paste0("<a href='", link_urls, "' target='_blank'>", df_display$Protein.Group, "</a>")
    df_display <- df_display %>% dplyr::select(-Original.ID)
    fixed_cols <- gdata$fixed_cols; expr_cols <- gdata$expr_cols
    valid_cols_map <- gdata$valid_cols_map; meta_sorted <- gdata$meta_sorted; group_colors <- gdata$group_colors

    header_html <- tags$table(class = "display", tags$thead(tags$tr(lapply(seq_along(colnames(df_display)), function(i) {
      col_name <- colnames(df_display)[i]
      if (i > length(fixed_cols)) {
        original_name <- valid_cols_map[i - length(fixed_cols)]; grp <- meta_sorted$Group[meta_sorted$File.Name == original_name]; bg_color <- group_colors[grp]
        tags$th(title = paste("File:", original_name, "\nGroup:", grp), col_name, style = paste0("background-color: ", bg_color, "; color: white; text-align: center;"))
      } else { tags$th(col_name) }
    }))))

    expression_matrix <- as.matrix(df_display[, expr_cols]); brks <- quantile(expression_matrix, probs = seq(.05, .95, .05), na.rm = TRUE); clrs <- colorRampPalette(c("#4575b4", "white", "#d73027"))(length(brks) + 1)
    # Find Type column index (0-based for JS)
    type_col_idx <- which(colnames(df_display) == "Type") - 1

    datatable(df_display, container = header_html, selection = 'single', escape = FALSE, options = list(dom = 'frtip', pageLength = 15, scrollX = TRUE, columnDefs = list(
      list(className = 'dt-center', targets = (length(fixed_cols)):(ncol(df_display)-1)),
      list(visible = FALSE, targets = type_col_idx)  # hide Type column but keep for styling
    )), rownames = FALSE) %>%
      formatStyle(expr_cols, backgroundColor = styleInterval(brks, clrs)) %>%
      formatStyle("Protein.Group", "Type",
        backgroundColor = styleEqual("Contaminant", "#fff0f0"),
        color = styleEqual("Contaminant", "#cc3333"))
  })

  output$download_grid_data <- downloadHandler(
    filename = function() { paste0("DE_LIMP_Grid_Export_", Sys.Date(), ".csv") },
    content = function(file) {
      gdata <- grid_react_df(); df_export <- gdata$data %>% dplyr::select(-Original.ID)
      fixed_cols <- gdata$fixed_cols; real_names <- gdata$valid_cols_map
      colnames(df_export) <- c(fixed_cols, real_names)
      write.csv(df_export, file, row.names = FALSE)

      # Log export
      add_to_log("Export Grid View to CSV", c(
        sprintf("# Exported full expression matrix with DE stats"),
        sprintf("# File: DE_LIMP_Grid_Export_%s.csv", Sys.Date()),
        "# (Combined DE results + expression values for all samples)"
      ))
    }
  )

  observeEvent(input$grid_view_table_rows_selected, {
    req(grid_react_df()); selected_idx <- input$grid_view_table_rows_selected
    if (length(selected_idx) > 0) {
      gdata <- grid_react_df(); selected_id <- gdata$data$Original.ID[selected_idx]; values$grid_selected_protein <- selected_id
      xic_btn <- if (!is_hf_space) actionButton("show_xic_from_grid", "\U0001F4C8 XICs", class="btn-info") else NULL
      showModal(modalDialog(title = paste("Expression Plot:", selected_id), size = "xl", plotOutput("violin_plot_grid", height = "600px"), footer = tagList(xic_btn, actionButton("back_to_grid", "Back to Grid", class="btn-info"), modalButton("Close")), easyClose = TRUE))
    }
  })

  observeEvent(input$back_to_grid, { click("show_grid_view") })

  output$violin_plot_grid <- renderPlot({
    req(values$y_protein, values$grid_selected_protein, values$metadata)
    prot_id <- values$grid_selected_protein
    exprs_mat <- values$y_protein$E[prot_id, , drop=FALSE]
    long_df <- as.data.frame(exprs_mat) %>% rownames_to_column("Protein") %>% pivot_longer(-Protein, names_to = "File.Name", values_to = "LogIntensity")
    long_df <- left_join(long_df, values$metadata, by="File.Name")

    ggplot(long_df, aes(x = Group, y = LogIntensity, fill = Group)) +
      geom_violin(alpha = 0.5, trim = FALSE) + geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      facet_wrap(~Protein, scales = "free_y") + theme_bw() +
      labs(title = paste("Protein:", prot_id), y = "Log2 Intensity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, height = 600) # FIXED HEIGHT

  # --- SIGNAL DISTRIBUTION (lines 1193-1242) ---
  output$protein_signal_plot <- renderPlot({
    req(values$y_protein)
    avg_signal <- rowMeans(values$y_protein$E, na.rm = TRUE)
    plot_df <- data.frame(Protein.Group = names(avg_signal), Average_Signal_Log2 = avg_signal) %>%
      mutate(Average_Signal_Log10 = Average_Signal_Log2 / log2(10))

    # Tag contaminants
    plot_df$Is_Contaminant <- grepl("^Cont_", plot_df$Protein.Group)

    # Base view: always show sample proteins
    sample_df <- plot_df[!plot_df$Is_Contaminant, ]
    contam_df <- plot_df[plot_df$Is_Contaminant, ]
    show_contam <- isTRUE(input$signal_overlay_contam)

    # DE coloring for sample proteins
    if (!is.null(values$fit) && !is.null(input$contrast_selector_signal) && nchar(input$contrast_selector_signal) > 0) {
      de_data_raw <- topTable(values$fit, coef = input$contrast_selector_signal, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(de_data_raw)) {
        de_data_intermediate <- de_data_raw %>% rownames_to_column("Protein.Group")
      } else {
        de_data_intermediate <- de_data_raw
      }
      de_data <- de_data_intermediate %>%
        mutate(DE_Status = case_when(
          adj.P.Val < 0.05 & logFC > input$logfc_cutoff ~ "Up-regulated",
          adj.P.Val < 0.05 & logFC < -input$logfc_cutoff ~ "Down-regulated",
          TRUE ~ "Not Significant"
        )) %>%
        dplyr::select(Protein.Group, DE_Status)
      sample_df <- left_join(sample_df, de_data, by = "Protein.Group")
      sample_df$DE_Status[is.na(sample_df$DE_Status)] <- "Not Significant"
    } else {
      sample_df$DE_Status <- "Not Significant"
    }
    contam_df$DE_Status <- "Contaminant"

    # Combine for plotting
    plot_df <- if (show_contam) rbind(sample_df, contam_df) else sample_df

    if (!is.null(values$plot_selected_proteins)) {
      plot_df$Is_Selected <- plot_df$Protein.Group %in% values$plot_selected_proteins
    } else {
      plot_df$Is_Selected <- FALSE
    }
    selected_df <- filter(plot_df, Is_Selected)

    # Color palette — includes contaminant color
    color_values <- c(
      "Up-regulated" = "#e41a1c", "Down-regulated" = "#377eb8",
      "Not Significant" = "grey70", "Contaminant" = "#ff8c00"
    )

    # Build plot
    p <- ggplot(plot_df, aes(x = reorder(Protein.Group, -Average_Signal_Log10), y = Average_Signal_Log10))
    if (!is.null(values$fit) || show_contam) {
      p <- p + geom_point(aes(color = DE_Status), size = 1.5) +
        scale_color_manual(name = "Status", values = color_values)
    } else {
      p <- p + geom_point(color = "cornflowerblue", size = 1.5)
    }

    n_label <- sprintf("(%s proteins)", format(nrow(plot_df), big.mark = ","))
    p + labs(title = paste("Signal Distribution", n_label), x = NULL, y = "Average Signal (Log10 Intensity)") +
      theme_minimal() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_discrete(expand = expansion(add = 1)) +
      geom_point(data = selected_df, color = "black", shape = 1, size = 4, stroke = 1) +
      geom_text_repel(data = selected_df, aes(label = Protein.Group), size = 4, max.overlaps = 20)
  }, height = 500) # FIXED HEIGHT

  # --- SAMPLE CORRELATION HEATMAP ---
  # Helper: build correlation heatmap object (returns Heatmap, caller draws it)
  build_correlation_heatmap <- function(font_size = 9, cell_font_size = NULL) {
    mat <- values$y_protein$E
    mat <- mat[complete.cases(mat), ]
    if (nrow(mat) < 2 || ncol(mat) < 2) return(NULL)

    cor_mat <- cor(mat, use = "pairwise.complete.obs")

    # Build group annotation
    meta <- values$metadata
    col_names <- colnames(cor_mat)
    group_vec <- meta$Group[match(col_names, meta$File.Name)]
    group_vec[is.na(group_vec) | group_vec == ""] <- "Unassigned"
    groups <- factor(group_vec)
    group_colors <- setNames(rainbow(length(levels(groups))), levels(groups))

    # Short sample labels
    sample_labels <- paste0("S", seq_along(col_names))
    rownames(cor_mat) <- sample_labels
    colnames(cor_mat) <- sample_labels

    ha <- HeatmapAnnotation(
      Group = groups,
      col = list(Group = group_colors),
      show_annotation_name = TRUE
    )

    col_fun <- circlize::colorRamp2(
      c(min(cor_mat, na.rm = TRUE), mean(c(min(cor_mat, na.rm = TRUE), 1)), 1),
      c("#2166AC", "white", "#B2182B")
    )

    n_samples <- ncol(cor_mat)
    default_cell_size <- if (is.null(cell_font_size)) {
      if (n_samples <= 8) 9 else 7
    } else cell_font_size
    cell_fn <- if (n_samples <= 12) {
      function(j, i, x, y, width, height, fill) {
        grid::grid.text(sprintf("%.2f", cor_mat[i, j]), x, y,
          gp = grid::gpar(fontsize = default_cell_size))
      }
    } else NULL

    Heatmap(
      cor_mat,
      name = "Pearson r",
      col = col_fun,
      top_annotation = ha,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      cell_fun = cell_fn,
      row_names_gp = grid::gpar(fontsize = font_size),
      column_names_gp = grid::gpar(fontsize = font_size),
      column_title = "Sample Correlation Heatmap (Pearson r)"
    )
  }

  # Render to temp PNG with explicit dimensions (avoids zero-dimension container crash)
  render_heatmap_png <- function(font_size = 9, cell_font_size = NULL,
                                  width = 700, height = 500) {
    ht <- build_correlation_heatmap(font_size, cell_font_size)
    if (is.null(ht)) return(NULL)
    tmp <- tempfile(fileext = ".png")
    png(tmp, width = width, height = height, res = 96)
    ComplexHeatmap::draw(ht)
    dev.off()
    tmp
  }

  output$correlation_heatmap <- renderImage({
    req(values$y_protein, values$metadata)
    tmp <- render_heatmap_png(font_size = 9, width = 700, height = 500)
    req(tmp)
    list(src = tmp, width = 700, height = 500, alt = "Sample Correlation Heatmap")
  }, deleteFile = TRUE)

  # --- Fullscreen Correlation Heatmap ---
  observeEvent(input$fullscreen_corr_heatmap, {
    req(values$y_protein, values$metadata)
    showModal(modalDialog(
      title = "Sample Correlation Heatmap - Fullscreen View",
      renderImage({
        tmp <- render_heatmap_png(font_size = 11, cell_font_size = 10,
                                   width = 1000, height = 700)
        req(tmp)
        list(src = tmp, width = 1000, height = 700, alt = "Sample Correlation Heatmap")
      }, deleteFile = TRUE),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })

  # --- PER-GROUP REPLICATE STATISTICS ---
  replicate_stats_data <- reactive({
    req(values$y_protein, values$metadata)
    meta <- values$metadata
    mat <- values$y_protein$E
    groups <- unique(meta$Group[meta$Group != ""])
    total_proteins <- nrow(mat)

    stats_list <- lapply(groups, function(g) {
      files <- meta$File.Name[meta$Group == g]
      group_cols <- intersect(colnames(mat), files)
      n_samples <- length(group_cols)
      if (n_samples == 0) return(NULL)

      group_mat <- mat[, group_cols, drop = FALSE]

      # Median CV: log2 -> linear -> SD/mean * 100
      linear_mat <- 2^group_mat
      cvs <- apply(linear_mat, 1, function(x) {
        x <- x[!is.na(x)]
        if (length(x) > 1) (sd(x) / mean(x)) * 100 else NA
      })
      median_cv <- round(median(cvs, na.rm = TRUE), 2)

      # Mean within-group pairwise correlation
      if (n_samples > 1) {
        group_cor <- cor(group_mat, use = "pairwise.complete.obs")
        # Extract upper triangle (exclude diagonal)
        upper_vals <- group_cor[upper.tri(group_cor)]
        mean_cor <- round(mean(upper_vals, na.rm = TRUE), 4)
      } else {
        mean_cor <- NA
      }

      # Proteins in all reps (no NA in any replicate)
      proteins_all_reps <- sum(complete.cases(group_mat))
      completeness <- round(100 * proteins_all_reps / total_proteins, 2)

      data.frame(
        Group = g,
        `N Samples` = n_samples,
        `Median CV (%)` = median_cv,
        `Mean Correlation` = mean_cor,
        `Proteins in All Reps` = proteins_all_reps,
        `Completeness (%)` = completeness,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, stats_list)
  })

  output$replicate_stats_table <- renderDT({
    df <- replicate_stats_data()
    req(df)
    dt <- datatable(df, options = list(dom = 't', pageLength = 20), rownames = FALSE)
    dt <- formatStyle(dt, "Median CV (%)",
      backgroundColor = styleInterval(c(20, 35), c("#d4edda", "#fff3cd", "#f8d7da")))
    dt <- formatStyle(dt, "Mean Correlation",
      backgroundColor = styleInterval(c(0.90, 0.95), c("#f8d7da", "#fff3cd", "#d4edda")))
    dt <- formatStyle(dt, "Completeness (%)",
      backgroundColor = styleInterval(c(70, 90), c("#f8d7da", "#fff3cd", "#d4edda")))
    dt
  })

  # --- Replicate Consistency CSV Export ---
  output$download_replicate_csv <- downloadHandler(
    filename = function() {
      paste0("Replicate_Consistency_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      df <- replicate_stats_data()
      req(df)
      write.csv(df, file, row.names = FALSE)
    }
  )

  # ============================================================================
  #  PCA Plot (Data Overview > PCA sub-tab)
  # ============================================================================

  # Reactive: PCA computation
  pca_result <- reactive({
    req(values$y_protein)
    mat <- values$y_protein$E
    mat_complete <- mat[complete.cases(mat), ]
    req(nrow(mat_complete) > 2, ncol(mat_complete) > 2)
    prcomp(t(mat_complete), center = TRUE, scale. = TRUE)
  })

  # Observer: dynamically populate PCA color selector (same pattern as MDS in server_qc.R)
  observeEvent(values$metadata, {
    req(values$metadata)
    meta <- values$metadata
    color_choices <- "Group"
    if ("Batch" %in% colnames(meta) && any(nzchar(meta$Batch)))
      color_choices <- c(color_choices, "Batch")
    cov1_name <- if (!is.null(values$cov1_name) && nzchar(values$cov1_name)) values$cov1_name else "Covariate1"
    cov2_name <- if (!is.null(values$cov2_name) && nzchar(values$cov2_name)) values$cov2_name else "Covariate2"
    if ("Covariate1" %in% colnames(meta) && any(nzchar(meta$Covariate1)))
      color_choices <- c(color_choices, setNames("Covariate1", cov1_name))
    if ("Covariate2" %in% colnames(meta) && any(nzchar(meta$Covariate2)))
      color_choices <- c(color_choices, setNames("Covariate2", cov2_name))
    updateSelectInput(session, "pca_color_by", choices = color_choices, selected = "Group")
  })

  # Helper: build PCA ggplot object (shared between main and fullscreen render)
  build_pca_plot <- function(height_mode = "main") {
    pca <- pca_result()
    meta <- values$metadata

    # Parse axis selection
    axes <- strsplit(input$pca_axes %||% "1_2", "_")[[1]]
    pc_x <- as.integer(axes[1])
    pc_y <- as.integer(axes[2])

    # Variance explained
    var_pct <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

    # Build plot data
    scores <- as.data.frame(pca$x)
    scores$Sample <- rownames(scores)
    scores <- merge(scores, meta, by.x = "Sample", by.y = "File.Name", all.x = TRUE)

    # Color variable
    color_by <- input$pca_color_by %||% "Group"
    col_name <- if (color_by %in% colnames(scores)) color_by else "Group"
    scores$ColorVar <- scores[[col_name]]
    scores$ColorVar[is.na(scores$ColorVar) | scores$ColorVar == ""] <- "(unassigned)"

    # Resolve display label for legend
    color_label <- color_by
    if (color_by == "Covariate1" && !is.null(values$cov1_name) && nzchar(values$cov1_name))
      color_label <- values$cov1_name
    if (color_by == "Covariate2" && !is.null(values$cov2_name) && nzchar(values$cov2_name))
      color_label <- values$cov2_name

    pc_x_col <- paste0("PC", pc_x)
    pc_y_col <- paste0("PC", pc_y)

    # Colorblind-friendly palette (same as MDS)
    cb_pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999",
                "#000000", "#E41A1C", "#377EB8", "#4DAF4A")
    n_groups <- length(unique(scores$ColorVar))
    pal <- if (n_groups <= length(cb_pal)) cb_pal[1:n_groups] else rainbow(n_groups, v = 0.85, s = 0.8)

    p <- ggplot(scores, aes(x = .data[[pc_x_col]], y = .data[[pc_y_col]],
                             color = ColorVar, text = paste0("Sample: ", Sample, "\n", color_label, ": ", ColorVar))) +
      stat_ellipse(aes(group = ColorVar), level = 0.95, linetype = "dashed", show.legend = FALSE) +
      geom_point(size = 3, alpha = 0.85) +
      scale_color_manual(values = pal, name = color_label) +
      labs(
        x = sprintf("PC%d (%.1f%%)", pc_x, var_pct[pc_x]),
        y = sprintf("PC%d (%.1f%%)", pc_y, var_pct[pc_y]),
        title = "Principal Component Analysis"
      ) +
      theme_bw(base_size = 13) +
      theme(legend.position = "right")

    p
  }

  # Render: PCA plotly scatter
  output$pca_plot <- renderPlotly({
    req(pca_result())
    p <- build_pca_plot("main")
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "v"))
  })

  # Fullscreen handler
  observeEvent(input$fullscreen_pca, {
    showModal(modalDialog(
      title = "PCA - Fullscreen View",
      plotlyOutput("pca_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$pca_plot_fs <- renderPlotly({
    req(pca_result())
    p <- build_pca_plot("fullscreen")
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "v"))
  })

  # PNG export
  output$download_pca_png <- downloadHandler(
    filename = function() {
      paste0("PCA_", format(Sys.Date(), "%Y%m%d"), ".png")
    },
    content = function(file) {
      req(pca_result())
      p <- build_pca_plot("export")
      ggsave(file, plot = p, width = 10, height = 7, dpi = 150, bg = "white")
    }
  )

  # PCA Info Modal
  observeEvent(input$pca_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About PCA"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Principal Component Analysis"),
        p("PCA reduces your high-dimensional protein expression data into a small number of ",
          "principal components that capture the most variation. Each point is one sample, ",
          "and the axes show the directions of greatest variance in the data."),
        tags$h6("What 'good' looks like"),
        tags$ul(
          tags$li("Samples from the same group cluster together"),
          tags$li("Groups are well-separated along PC1 or PC2"),
          tags$li("The first two PCs explain a large fraction of total variance (shown in axis labels)")
        ),
        tags$h6("What 'bad' looks like"),
        tags$ul(
          tags$li("A single outlier sample dominates PC1 — consider removing it"),
          tags$li("Groups overlap completely — the biological effect may be subtle or absent"),
          tags$li("Samples cluster by batch instead of group — add batch as a covariate")
        ),
        tags$h6("PCA vs MDS"),
        p("Both are dimensionality reduction techniques. PCA is based on variance decomposition ",
          "of the expression matrix; MDS is based on pairwise sample distances. They usually give ",
          "very similar results, but PCA additionally provides variance explained percentages ",
          "and loading information for each component."),
        tags$h6("Controls"),
        tags$ul(
          tags$li(strong("Color by: "), "Switch between Group, Batch, or covariates to check for confounding"),
          tags$li(strong("Axes: "), "Switch between PC1/2, PC1/3, PC2/3 to explore additional dimensions"),
          tags$li("Dashed ellipses show the 95% confidence region for each group")
        )
      )
    ))
  })

  # ============================================================================
  #  CONTAMINANT ANALYSIS
  # ============================================================================

  # Reactive: contaminant analysis data
  contaminant_data <- reactive({
    req(values$y_protein)
    mat <- values$y_protein$E
    protein_ids <- rownames(mat)
    is_contam <- grepl("^Cont_", protein_ids)

    if (sum(is_contam) == 0) return(NULL)

    # Expression matrix is log2 — convert to linear for intensity sums
    linear_mat <- 2^mat

    # Per-sample intensity breakdown
    contam_intensity <- colSums(linear_mat[is_contam, , drop = FALSE], na.rm = TRUE)
    sample_intensity <- colSums(linear_mat[!is_contam, , drop = FALSE], na.rm = TRUE)
    total_intensity <- contam_intensity + sample_intensity
    contam_pct <- round(100 * contam_intensity / total_intensity, 2)

    sample_df <- data.frame(
      Sample = colnames(mat),
      Contaminant = contam_intensity,
      Sample_Protein = sample_intensity,
      Total = total_intensity,
      Contaminant_Pct = contam_pct,
      stringsAsFactors = FALSE
    )

    # Add group info from metadata
    if (!is.null(values$metadata)) {
      sample_df$Group <- values$metadata$Group[match(sample_df$Sample, values$metadata$File.Name)]
    } else {
      sample_df$Group <- "Unknown"
    }

    # Top contaminants by total intensity
    contam_mat <- linear_mat[is_contam, , drop = FALSE]
    avg_intensity <- rowMeans(contam_mat, na.rm = TRUE)
    total_int <- rowSums(contam_mat, na.rm = TRUE)
    presence <- rowSums(!is.na(mat[is_contam, , drop = FALSE]))
    overall_total <- sum(linear_mat, na.rm = TRUE)

    # Get gene info from y_protein$genes if available
    genes_df <- values$y_protein$genes
    contam_ids <- protein_ids[is_contam]
    gene_names <- rep("", length(contam_ids))
    protein_names <- rep("", length(contam_ids))

    if (!is.null(genes_df)) {
      gn_col <- intersect(c("Genes", "Gene.Names", "Gene"), colnames(genes_df))
      pn_col <- intersect(c("Protein.Names", "Protein.Name"), colnames(genes_df))

      if (length(gn_col) > 0) {
        # Match by row position — y_protein$genes rows correspond to y_protein$E rows
        contam_idx <- which(is_contam)
        gene_names <- as.character(genes_df[[gn_col[1]]][contam_idx])
        gene_names[is.na(gene_names)] <- ""
      }
      if (length(pn_col) > 0) {
        contam_idx <- which(is_contam)
        protein_names <- as.character(genes_df[[pn_col[1]]][contam_idx])
        protein_names[is.na(protein_names)] <- ""
      }
    }

    # Identify keratins
    is_keratin <- grepl("^(KRT|K1C|K2C|K22)", gene_names, ignore.case = TRUE) |
                  grepl("keratin", protein_names, ignore.case = TRUE)

    contam_df <- data.frame(
      Protein.Group = contam_ids,
      Gene = gene_names,
      Protein.Name = protein_names,
      Avg_Intensity = round(avg_intensity, 0),
      Pct_of_Total = round(100 * total_int / overall_total, 3),
      Present_in = paste0(presence, "/", ncol(mat)),
      Is_Keratin = is_keratin,
      stringsAsFactors = FALSE
    ) %>% arrange(desc(Avg_Intensity))

    list(
      n_contam = sum(is_contam),
      n_sample = sum(!is_contam),
      pct_contam = round(100 * sum(is_contam) / length(protein_ids), 1),
      sample_df = sample_df,
      contam_df = contam_df,
      contam_mat = mat[is_contam, , drop = FALSE]  # keep log2 for heatmap
    )
  })

  # Summary cards
  output$contaminant_summary_cards <- renderUI({
    cdata <- contaminant_data()

    if (is.null(cdata)) {
      return(div(
        style = "background-color: #d4edda; padding: 20px; border-radius: 8px; text-align: center;",
        icon("check-circle", style = "color: #28a745; font-size: 1.5em;"),
        tags$h5("No contaminant proteins detected", style = "color: #28a745; margin-top: 10px;"),
        tags$p("No proteins with 'Cont_' prefix found in the dataset.",
          style = "color: #6c757d;")
      ))
    }

    median_pct <- round(median(cdata$sample_df$Contaminant_Pct), 2)
    max_pct <- round(max(cdata$sample_df$Contaminant_Pct), 2)
    n_keratins <- sum(cdata$contam_df$Is_Keratin)

    # Color code the contamination level
    pct_color <- if (median_pct < 1) "#28a745" else if (median_pct < 5) "#ffc107" else "#dc3545"
    pct_bg <- if (median_pct < 1) "#d4edda" else if (median_pct < 5) "#fff3cd" else "#f8d7da"

    div(style = "display: flex; flex-wrap: wrap; gap: 12px;",
      # Contaminant count
      div(style = "flex: 1; min-width: 150px; background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #6c757d;",
        tags$p(style = "font-size: 0.85em; color: #6c757d; margin-bottom: 4px;", "Contaminant Proteins"),
        tags$h4(style = "margin: 0; font-weight: 600;", cdata$n_contam),
        tags$p(style = "font-size: 0.8em; color: #6c757d; margin: 0;",
          paste0(cdata$pct_contam, "% of ", cdata$n_contam + cdata$n_sample, " total"))
      ),
      # Sample proteins
      div(style = "flex: 1; min-width: 150px; background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #0d6efd;",
        tags$p(style = "font-size: 0.85em; color: #6c757d; margin-bottom: 4px;", "Sample Proteins"),
        tags$h4(style = "margin: 0; font-weight: 600;", cdata$n_sample)
      ),
      # Median contaminant intensity %
      div(style = paste0("flex: 1; min-width: 150px; background-color: ", pct_bg, "; padding: 15px; border-radius: 8px; border-left: 4px solid ", pct_color, ";"),
        tags$p(style = "font-size: 0.85em; color: #6c757d; margin-bottom: 4px;", "Median Contam. Intensity"),
        tags$h4(style = paste0("margin: 0; font-weight: 600; color: ", pct_color, ";"),
          paste0(median_pct, "%")),
        tags$p(style = "font-size: 0.8em; color: #6c757d; margin: 0;",
          paste0("Max: ", max_pct, "%"))
      ),
      # Keratins
      div(style = paste0("flex: 1; min-width: 150px; background-color: ", if (n_keratins > 0) "#fff3cd" else "#f8f9fa",
        "; padding: 15px; border-radius: 8px; border-left: 4px solid ", if (n_keratins > 0) "#ffc107" else "#6c757d", ";"),
        tags$p(style = "font-size: 0.85em; color: #6c757d; margin-bottom: 4px;", "Keratin Contaminants"),
        tags$h4(style = "margin: 0; font-weight: 600;", n_keratins),
        tags$p(style = "font-size: 0.8em; color: #6c757d; margin: 0;",
          if (n_keratins > 0) "KRT/K1C/K2C detected" else "None detected")
      )
    )
  })

  # Per-sample contaminant bar chart
  output$contaminant_bar_chart <- renderPlotly({
    cdata <- contaminant_data()
    req(cdata)

    df <- cdata$sample_df %>% arrange(desc(Contaminant_Pct))

    # Use short sample IDs if metadata available
    if (!is.null(values$metadata)) {
      df$Label <- paste0("S", values$metadata$ID[match(df$Sample, values$metadata$File.Name)])
    } else {
      df$Label <- df$Sample
    }
    df$Label <- factor(df$Label, levels = df$Label)

    plot_ly(df, y = ~Label, x = ~Contaminant, type = "bar", orientation = "h",
            name = "Contaminant", marker = list(color = "#dc3545"),
            hoverinfo = "text",
            text = ~paste0(Sample, "\nContaminant: ", format(round(Contaminant), big.mark = ","),
                          "\n", Contaminant_Pct, "% of total")) %>%
      add_trace(x = ~Sample_Protein, name = "Sample", marker = list(color = "#0d6efd"),
                hoverinfo = "text",
                text = ~paste0(Sample, "\nSample: ", format(round(Sample_Protein), big.mark = ","),
                              "\n", round(100 - Contaminant_Pct, 2), "% of total")) %>%
      layout(
        barmode = "stack",
        xaxis = list(title = "Total Intensity (linear)"),
        yaxis = list(title = "", categoryorder = "trace"),
        legend = list(orientation = "h", xanchor = "center", x = 0.5, y = 1.05),
        margin = list(l = 60)
      )
  })

  # Top contaminants table
  output$contaminant_top_table <- renderDT({
    cdata <- contaminant_data()
    req(cdata)

    df <- cdata$contam_df %>%
      dplyr::select(Protein.Group, Gene, Protein.Name, Avg_Intensity, Pct_of_Total, Present_in, Is_Keratin)

    colnames(df) <- c("Protein Group", "Gene", "Protein Name", "Avg Intensity", "% of Total", "Present In", "Keratin")
    df$Keratin <- ifelse(df$Keratin, "Yes", "")

    datatable(df,
      options = list(dom = 'frtip', pageLength = 20, scrollX = TRUE,
        order = list(list(3, "desc"))),
      rownames = FALSE
    ) %>%
      formatStyle("Keratin",
        backgroundColor = styleEqual("Yes", "#fff3cd"),
        fontWeight = styleEqual("Yes", "bold")
      ) %>%
      formatRound("Avg Intensity", digits = 0) %>%
      formatRound("% of Total", digits = 3)
  })

  # Contaminant heatmap
  output$contaminant_heatmap <- renderPlotly({
    cdata <- contaminant_data()
    req(cdata)

    # Top 20 contaminants by average intensity
    top_n <- min(20, nrow(cdata$contam_df))
    top_ids <- cdata$contam_df$Protein.Group[1:top_n]

    hm_mat <- cdata$contam_mat[top_ids, , drop = FALSE]

    # Use gene names as labels where available
    row_labels <- cdata$contam_df$Gene[1:top_n]
    row_labels <- ifelse(nzchar(row_labels), row_labels, cdata$contam_df$Protein.Group[1:top_n])
    # Mark keratins
    is_ker <- cdata$contam_df$Is_Keratin[1:top_n]
    row_labels <- ifelse(is_ker, paste0(row_labels, " *"), row_labels)

    # Use short sample IDs
    if (!is.null(values$metadata)) {
      col_labels <- paste0("S", values$metadata$ID[match(colnames(hm_mat), values$metadata$File.Name)])
    } else {
      col_labels <- colnames(hm_mat)
    }

    plot_ly(
      x = col_labels,
      y = row_labels,
      z = hm_mat,
      type = "heatmap",
      colors = colorRamp(c("#2166AC", "#F7F7F7", "#B2182B")),
      hoverinfo = "text",
      text = outer(1:nrow(hm_mat), 1:ncol(hm_mat), Vectorize(function(i, j) {
        paste0(row_labels[i], "\n", col_labels[j], "\nLog2 Intensity: ", round(hm_mat[i, j], 2))
      }))
    ) %>%
      layout(
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "", autorange = "reversed"),
        margin = list(l = 120, b = 80),
        annotations = list(
          list(text = "* = Keratin contaminant", x = 1, y = -0.15,
               xref = "paper", yref = "paper", showarrow = FALSE,
               font = list(size = 11, color = "#6c757d"))
        )
      )
  })

  # Contaminant info modal
  observeEvent(input$contaminant_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About Contaminant Analysis"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What are contaminant proteins?"),
        p("Contaminant proteins (prefixed with 'Cont_') come from the contaminant FASTA library ",
          "included in DIA-NN searches. These are common lab contaminants like keratins (skin), ",
          "trypsin (digestion enzyme), albumin (serum), and other environmental proteins."),
        tags$h6("Why monitor them?"),
        tags$ul(
          tags$li("High contaminant levels indicate sample preparation issues"),
          tags$li("Keratins suggest skin contact during prep — check glove usage"),
          tags$li("Sample-specific spikes may indicate individual prep failures"),
          tags$li("Consistent high contaminants across all samples may indicate reagent contamination")
        ),
        tags$h6("Interpreting the results"),
        tags$ul(
          tags$li(strong("< 1% contaminant intensity: "), "Excellent — typical for well-prepared samples"),
          tags$li(strong("1-5% contaminant intensity: "), "Acceptable — some contamination present"),
          tags$li(strong("> 5% contaminant intensity: "), "Concerning — investigate sample prep workflow")
        ),
        tags$h6("Heatmap"),
        p("The heatmap shows the top 20 contaminant proteins across all samples. ",
          "Look for samples with unusually high intensity (red) in specific contaminants, ",
          "which may indicate sample-specific issues. Keratin contaminants are marked with *."),
        tags$h6("Expression Grid"),
        p("Contaminant proteins are highlighted with a pink background and red text in the Expression Grid tab.")
      )
    ))
  })

  # --- FULLSCREEN: Signal Distribution (lines 1734-1790) ---
  observeEvent(input$fullscreen_signal, {
    showModal(modalDialog(
      title = "Signal Distribution - Fullscreen View",
      plotOutput("protein_signal_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$protein_signal_plot_fs <- renderPlot({
    req(values$y_protein)
    avg_signal <- rowMeans(values$y_protein$E, na.rm = TRUE)
    plot_df <- data.frame(Protein.Group = names(avg_signal), Average_Signal_Log2 = avg_signal) %>%
      mutate(Average_Signal_Log10 = Average_Signal_Log2 / log2(10))

    # Always show DE coloring when results are available
    if (!is.null(values$fit) && !is.null(input$contrast_selector_signal) && nchar(input$contrast_selector_signal) > 0) {
      de_data_raw <- topTable(values$fit, coef = input$contrast_selector_signal, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(de_data_raw)) {
        de_data_intermediate <- de_data_raw %>% rownames_to_column("Protein.Group")
      } else {
        de_data_intermediate <- de_data_raw
      }
      de_data <- de_data_intermediate %>%
        mutate(DE_Status = case_when(
          adj.P.Val < 0.05 & logFC > input$logfc_cutoff ~ "Up-regulated",
          adj.P.Val < 0.05 & logFC < -input$logfc_cutoff ~ "Down-regulated",
          TRUE ~ "Not Significant"
        )) %>%
        dplyr::select(Protein.Group, DE_Status)
      plot_df <- left_join(plot_df, de_data, by = "Protein.Group")
      plot_df$DE_Status[is.na(plot_df$DE_Status)] <- "Not Significant"
    } else {
      plot_df$DE_Status <- "Not Significant"
    }

    if (!is.null(values$plot_selected_proteins)) {
      plot_df$Is_Selected <- plot_df$Protein.Group %in% values$plot_selected_proteins
    } else {
      plot_df$Is_Selected <- FALSE
    }
    selected_df <- filter(plot_df, Is_Selected)

    # Build plot - always use DE coloring when available
    p <- ggplot(plot_df, aes(x = reorder(Protein.Group, -Average_Signal_Log10), y = Average_Signal_Log10))
    if (!is.null(values$fit)) {
      p <- p + geom_point(aes(color = DE_Status), size = 1.5) +
        scale_color_manual(name = "DE Status",
          values = c("Up-regulated" = "#e41a1c", "Down-regulated" = "#377eb8", "Not Significant" = "grey70"))
    } else {
      p <- p + geom_point(color = "cornflowerblue", size = 1.5)
    }

    p + labs(title = "Signal Distribution Across All Protein Groups", x = NULL, y = "Average Signal (Log10 Intensity)") +
      theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_discrete(expand = expansion(add = 1)) +
      geom_point(data = selected_df, color = "black", shape = 1, size = 4, stroke = 1) +
      geom_text_repel(data = selected_df, aes(label = Protein.Group), size = 4, max.overlaps = 20)
  }, height = 700)

  # ==============================================================================
  #  DATA EXPLORER — Abundance Profiles & Sample Scatter
  # ==============================================================================

  # --- Info modal ---
  observeEvent(input$data_explorer_info_btn, {
    showModal(modalDialog(
      title = "Data Explorer",
      size = "l",
      easyClose = TRUE,
      tags$h5("Abundance Profiles (Quartile Analysis)"),
      tags$p("Proteins are split into four quartiles (Q1 = highest, Q4 = lowest) by their average intensity across all samples. ",
        "The heatmap shows the top 10 proteins in each quartile (40 total). Each cell is colored by which quartile that protein falls in ",
        "FOR THAT SPECIFIC SAMPLE — not the average. This reveals proteins whose relative abundance shifts across samples."),
      tags$p("The 'Variable Proteins' table below lists proteins whose per-sample quartile assignment varies by 2 or more quartiles, ",
        "indicating potential biological regulation or technical variability."),
      tags$hr(),
      tags$h5("Sample-Sample Scatter"),
      tags$p("Select two samples to compare their protein intensities directly. Each point is a protein. ",
        "Points far from the identity line (y = x) represent proteins with large intensity differences between the two samples."),
      tags$p("Outliers (>4-fold difference) are labeled with gene names when available. ",
        "Contaminant proteins are shown as orange triangles when not excluded.")
    ))
  })

  # --- Helper: get gene label for a protein ID ---
  explorer_gene_label <- function(protein_id, genes_df) {
    if (is.null(genes_df)) return(protein_id)
    # Try Genes column first
    gene_col <- intersect(c("Genes", "Gene.Names", "Gene"), colnames(genes_df))
    if (length(gene_col) > 0) {
      idx <- match(protein_id, rownames(genes_df))
      if (!is.na(idx)) {
        g <- genes_df[idx, gene_col[1]]
        if (!is.na(g) && nzchar(g) && nchar(g) < 20 && !grepl("[;|]", g)) return(g)
      }
    }
    # Try parsing from sp|ACC|GENE format
    parsed <- sub("^sp\\|[^|]+\\|([^_]+)_.*$", "\\1", protein_id)
    if (parsed != protein_id && nchar(parsed) < 20) return(parsed)
    # Truncate long accessions
    if (nchar(protein_id) > 15) return(substr(protein_id, 1, 15))
    protein_id
  }

  # --- Quartile heatmap reactive ---
  explorer_quartile_data <- reactive({
    req(values$y_protein)
    mat <- values$y_protein$E
    genes_df <- values$y_protein$genes

    # Exclude contaminants if requested
    if (isTRUE(input$explorer_exclude_contam_profile)) {
      keep <- !grepl("^Cont_", rownames(mat))
      mat <- mat[keep, , drop = FALSE]
      if (!is.null(genes_df)) genes_df <- genes_df[keep, , drop = FALSE]
    }

    if (nrow(mat) < 4) return(NULL)

    # Step 1: Average intensity per protein
    avg_intensity <- rowMeans(mat, na.rm = TRUE)

    # Step 2: Assign quartiles based on average (Q1 = highest intensity)
    # Use rank-based assignment for guaranteed equal groups
    ranks <- rank(-avg_intensity, ties.method = "first")  # highest = rank 1
    n <- length(ranks)
    avg_quartile <- factor(
      ifelse(ranks <= n/4, "Q1",
        ifelse(ranks <= n/2, "Q2",
          ifelse(ranks <= 3*n/4, "Q3", "Q4"))),
      levels = c("Q1", "Q2", "Q3", "Q4")
    )
    names(avg_quartile) <- names(avg_intensity)

    # Step 3: Top 10 per quartile
    top_proteins <- character(0)
    for (q in c("Q1", "Q2", "Q3", "Q4")) {
      in_q <- names(avg_quartile)[avg_quartile == q]
      if (length(in_q) == 0) next
      ordered <- in_q[order(-avg_intensity[in_q])]
      top_proteins <- c(top_proteins, head(ordered, 10))
    }

    # Step 4: Per-sample quartile assignment for selected proteins
    heatmap_mat <- mat[top_proteins, , drop = FALSE]
    sample_quartiles <- matrix(NA_integer_, nrow = length(top_proteins), ncol = ncol(mat))
    rownames(sample_quartiles) <- top_proteins
    colnames(sample_quartiles) <- colnames(mat)

    for (j in seq_len(ncol(mat))) {
      col_vals <- mat[, j]
      col_ranks <- rank(-col_vals, ties.method = "first")
      col_n <- length(col_ranks)
      col_q <- ifelse(col_ranks <= col_n/4, 1L,
                ifelse(col_ranks <= col_n/2, 2L,
                  ifelse(col_ranks <= 3*col_n/4, 3L, 4L)))
      sample_quartiles[, j] <- col_q[match(top_proteins, names(col_vals))]
    }

    # Gene labels for rows
    row_labels <- vapply(top_proteins, function(pid) {
      explorer_gene_label(pid, genes_df)
    }, character(1))

    # Short sample labels
    if (!is.null(values$metadata)) {
      col_labels <- paste0("S", values$metadata$ID[match(colnames(mat), values$metadata$File.Name)])
      col_labels[is.na(col_labels)] <- colnames(mat)[is.na(col_labels)]
    } else {
      col_labels <- colnames(mat)
    }

    # Step 5: Variable proteins (quartile range >= 2 across ALL proteins)
    all_sample_q <- matrix(NA_integer_, nrow = nrow(mat), ncol = ncol(mat))
    rownames(all_sample_q) <- rownames(mat)
    for (j in seq_len(ncol(mat))) {
      col_vals <- mat[, j]
      col_ranks <- rank(-col_vals, ties.method = "first")
      col_n <- length(col_ranks)
      all_sample_q[, j] <- ifelse(col_ranks <= col_n/4, 1L,
                             ifelse(col_ranks <= col_n/2, 2L,
                               ifelse(col_ranks <= 3*col_n/4, 3L, 4L)))
    }

    min_q <- apply(all_sample_q, 1, min, na.rm = TRUE)
    max_q <- apply(all_sample_q, 1, max, na.rm = TRUE)
    q_range <- max_q - min_q
    variable_idx <- which(q_range >= 2)

    variable_df <- NULL
    if (length(variable_idx) > 0) {
      variable_df <- data.frame(
        Protein.Group = rownames(mat)[variable_idx],
        Gene = vapply(rownames(mat)[variable_idx], function(pid) explorer_gene_label(pid, genes_df), character(1)),
        Avg_Intensity = round(avg_intensity[variable_idx], 2),
        Min_Quartile = min_q[variable_idx],
        Max_Quartile = max_q[variable_idx],
        Quartile_Range = q_range[variable_idx],
        Samples = ncol(mat),
        stringsAsFactors = FALSE
      )
      variable_df <- variable_df[order(-variable_df$Quartile_Range, -variable_df$Avg_Intensity), ]
    }

    # Quartile group boundaries for divider lines
    group_sizes <- vapply(c("Q1", "Q2", "Q3", "Q4"), function(q) sum(avg_quartile[top_proteins] == q), integer(1))

    list(
      sample_quartiles = sample_quartiles,
      row_labels = row_labels,
      col_labels = col_labels,
      avg_quartile = avg_quartile[top_proteins],
      group_sizes = group_sizes,
      variable_df = variable_df
    )
  })

  # --- Quartile heatmap plot ---
  output$explorer_quartile_heatmap <- renderPlotly({
    qdata <- explorer_quartile_data()
    req(qdata)

    sq <- qdata$sample_quartiles
    row_labels <- qdata$row_labels
    avg_q <- qdata$avg_quartile
    # Don't reverse — use yaxis autorange="reversed" to put Q1 at top

    # Build hover text
    hover_text <- matrix("", nrow = nrow(sq), ncol = ncol(sq))
    for (i in seq_len(nrow(sq))) {
      for (j in seq_len(ncol(sq))) {
        hover_text[i, j] <- paste0(
          "Protein: ", row_labels[i],
          "<br>Sample: ", qdata$col_labels[j],
          "<br>Sample Quartile: Q", sq[i, j],
          "<br>Average Quartile: ", avg_q[i]
        )
      }
    }

    # Color scale: Q1 (highest, value 1) = dark blue, Q4 (lowest, value 4) = red
    colorscale <- list(
      list(0, "#08306b"),     # Q1 - darkest blue
      list(0.333, "#2171b5"), # Q2 - medium blue
      list(0.667, "#6baed6"), # Q3 - light blue
      list(1, "#ef3b2c")     # Q4 - red
    )

    p <- plot_ly(
      z = sq,
      x = qdata$col_labels,
      y = row_labels,
      type = "heatmap",
      colorscale = colorscale,
      zmin = 1, zmax = 4,
      text = hover_text,
      hoverinfo = "text",
      colorbar = list(
        title = "Quartile",
        tickvals = c(1, 2, 3, 4),
        ticktext = c("Q1 (High)", "Q2", "Q3", "Q4 (Low)"),
        len = 0.5
      )
    )

    # Add divider lines and quartile labels between groups
    # Row order (top to bottom with reversed y-axis): Q1, Q2, Q3, Q4
    gs <- qdata$group_sizes[c("Q1", "Q2", "Q3", "Q4")]
    shapes <- list()
    annotations <- list()
    cum <- 0
    q_labels <- c("Q1 (High)", "Q2", "Q3", "Q4 (Low)")
    for (i in seq_along(gs)) {
      # Section label at midpoint
      mid_y <- cum + gs[i] / 2 - 0.5
      annotations[[length(annotations) + 1]] <- list(
        text = q_labels[i], x = -0.15, y = mid_y,
        xref = "paper", yref = "y",
        showarrow = FALSE, font = list(size = 11, color = "#555", family = "Arial Black"),
        xanchor = "right"
      )

      cum <- cum + gs[i]
      # Divider line
      if (i < length(gs) && cum > 0 && cum < nrow(sq)) {
        shapes[[length(shapes) + 1]] <- list(
          type = "line",
          x0 = -0.5, x1 = ncol(sq) - 0.5,
          y0 = cum - 0.5, y1 = cum - 0.5,
          line = list(color = "white", width = 3)
        )
      }
    }

    p %>% layout(
      xaxis = list(title = "", tickangle = -45, tickfont = list(size = 10)),
      yaxis = list(title = "", tickfont = list(size = 9), dtick = 1, autorange = "reversed"),
      shapes = shapes,
      annotations = annotations,
      margin = list(l = 150, b = 80)
    )
  })

  # --- Variable proteins table ---
  output$explorer_variable_proteins_table <- renderDT({
    qdata <- explorer_quartile_data()
    req(qdata, qdata$variable_df)

    df <- qdata$variable_df
    colnames(df) <- c("Protein Group", "Gene", "Avg Intensity", "Min Quartile", "Max Quartile", "Quartile Range", "Samples")

    datatable(df,
      options = list(dom = 'frtip', pageLength = 15, scrollX = TRUE,
        order = list(list(5, "desc"))),
      rownames = FALSE
    ) %>%
      formatRound("Avg Intensity", digits = 1) %>%
      formatStyle("Quartile Range",
        backgroundColor = styleInterval(c(2, 3), c("#fff3cd", "#f8d7da", "#f8d7da")),
        fontWeight = "bold"
      )
  })

  # --- Update sample selectors when y_protein changes ---
  observeEvent(values$y_protein, {
    req(values$y_protein)
    sample_names <- colnames(values$y_protein$E)

    if (!is.null(values$metadata)) {
      short_ids <- paste0("S", values$metadata$ID[match(sample_names, values$metadata$File.Name)])
      short_ids[is.na(short_ids)] <- sample_names[is.na(short_ids)]
      choices <- setNames(sample_names, short_ids)
    } else {
      choices <- setNames(sample_names, sample_names)
    }

    updateSelectInput(session, "explorer_sample_a", choices = choices,
      selected = if (length(choices) >= 1) choices[1] else NULL)
    updateSelectInput(session, "explorer_sample_b", choices = choices,
      selected = if (length(choices) >= 2) choices[2] else NULL)
  })

  # --- Sample-Sample Scatter ---
  output$explorer_sample_scatter <- renderPlotly({
    req(values$y_protein, input$explorer_sample_a, input$explorer_sample_b)
    req(input$explorer_sample_a != input$explorer_sample_b)

    mat <- values$y_protein$E
    genes_df <- values$y_protein$genes
    sa <- input$explorer_sample_a
    sb <- input$explorer_sample_b
    req(sa %in% colnames(mat), sb %in% colnames(mat))

    # Build data frame
    df <- data.frame(
      Protein.Group = rownames(mat),
      A = mat[, sa],
      B = mat[, sb],
      stringsAsFactors = FALSE
    )
    df$is_contaminant <- grepl("^Cont_", df$Protein.Group)

    # Gene labels
    df$Gene <- vapply(df$Protein.Group, function(pid) explorer_gene_label(pid, genes_df), character(1))

    # Protein name for hover
    if (!is.null(genes_df) && "Protein.Names" %in% colnames(genes_df)) {
      name_idx <- match(df$Protein.Group, rownames(genes_df))
      df$Protein.Name <- ifelse(is.na(name_idx), "", genes_df$Protein.Names[name_idx])
    } else {
      df$Protein.Name <- ""
    }

    # Remove rows with NA in either sample
    df <- df[!is.na(df$A) & !is.na(df$B), ]

    # Exclude contaminants if requested
    exclude_contam <- isTRUE(input$explorer_exclude_contam_scatter)
    if (exclude_contam) {
      df <- df[!df$is_contaminant, ]
    }

    if (nrow(df) < 2) return(NULL)

    # Compute diff and stats
    df$diff <- df$A - df$B
    df$abs_diff <- abs(df$diff)
    df$fold_diff <- round(2^df$abs_diff, 1)
    df$is_outlier <- df$abs_diff > 2  # 4-fold

    pearson_r <- round(cor(df$A, df$B, use = "complete.obs"), 3)
    n_proteins <- nrow(df)
    n_outliers <- sum(df$is_outlier)

    # Short sample labels for axis
    sa_label <- sa
    sb_label <- sb
    if (!is.null(values$metadata)) {
      sa_label <- paste0("S", values$metadata$ID[match(sa, values$metadata$File.Name)])
      sb_label <- paste0("S", values$metadata$ID[match(sb, values$metadata$File.Name)])
      if (is.na(sa_label)) sa_label <- sa
      if (is.na(sb_label)) sb_label <- sb
    }

    # Color by distance from identity
    df$color_val <- pmin(df$abs_diff, 4)  # cap for color scaling

    # Hover text
    df$hover <- paste0(
      "Gene: ", df$Gene,
      "<br>Protein: ", df$Protein.Name,
      "<br>", sa_label, ": ", round(df$A, 2),
      "<br>", sb_label, ": ", round(df$B, 2),
      "<br>Fold diff: ", df$fold_diff, "x"
    )

    # Subtitle
    subtitle_text <- paste0(
      "Pearson r = ", pearson_r,
      " | N = ", format(n_proteins, big.mark = ","), " proteins",
      " | ", n_outliers, " with >4-fold difference"
    )

    # Separate contaminants and sample proteins
    contam_df <- df[df$is_contaminant, ]
    sample_df <- df[!df$is_contaminant, ]
    outlier_df <- df[df$is_outlier & !df$is_contaminant, ]

    # Build ggplot
    axis_range <- range(c(df$A, df$B), na.rm = TRUE)
    axis_pad <- diff(axis_range) * 0.05

    p <- ggplot() +
      # Identity line
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey60", linewidth = 0.5) +
      # Sample proteins colored by distance
      geom_point(data = sample_df, aes(x = A, y = B, color = color_val, text = hover),
        size = 1.5, alpha = 0.7) +
      scale_color_gradient(low = "grey70", high = "#e41a1c", name = "|Diff| (log2)",
        limits = c(0, 4), guide = "colorbar")

    # Add contaminants as orange triangles if not excluded
    if (!exclude_contam && nrow(contam_df) > 0) {
      p <- p + geom_point(data = contam_df, aes(x = A, y = B, text = hover),
        shape = 17, color = "#ff8c00", size = 2.5, alpha = 0.8)
    }

    # Label outliers
    if (isTRUE(input$explorer_label_outliers) && nrow(outlier_df) > 0) {
      # Show top 30 outliers max to avoid clutter
      outlier_label_df <- head(outlier_df[order(-outlier_df$abs_diff), ], 30)
      p <- p + geom_text(data = outlier_label_df, aes(x = A, y = B, label = Gene),
        size = 3, hjust = -0.15, vjust = -0.3, check_overlap = TRUE, color = "#333333")
    }

    p <- p +
      labs(
        title = paste0("Sample Scatter: ", sa_label, " vs ", sb_label),
        subtitle = subtitle_text,
        x = paste0(sa_label, " (log2 intensity)"),
        y = paste0(sb_label, " (log2 intensity)")
      ) +
      coord_fixed(ratio = 1,
        xlim = c(axis_range[1] - axis_pad, axis_range[2] + axis_pad),
        ylim = c(axis_range[1] - axis_pad, axis_range[2] + axis_pad)) +
      theme_bw() +
      theme(plot.subtitle = element_text(size = 10, color = "#555555"))

    ggplotly(p, tooltip = "text") %>%
      layout(
        margin = list(t = 60),
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.15)
      )
  })

  # ==============================================================================
  #  DATA EXPLORER — Export for Claude
  # ==============================================================================

  output$export_explorer_claude <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_Explorer_Claude_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip")
    },
    content = function(file) {
      req(values$y_protein)

      tryCatch({
      withProgress(message = "Building explorer export...", value = 0, {
        tmp_dir <- file.path(tempdir(), paste0("explorer_claude_", format(Sys.time(), "%Y%m%d%H%M%S")))
        dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
        files_to_zip <- character(0)

        mat <- values$y_protein$E
        genes_df <- values$y_protein$genes
        n_proteins <- nrow(mat)
        n_samples <- ncol(mat)

        # --- 0. Load NCBI gene map if applicable ---
        export_gene_map <- NULL
        n_gene_mapped <- 0L
        non_contam_ids <- rownames(mat)[!grepl("^Cont_", rownames(mat))]
        first_accessions <- sub(";.*", "", head(non_contam_ids, 50))
        is_ncbi_export <- length(first_accessions) > 0 && any(grepl("^[XNW]P_", first_accessions))

        if (is_ncbi_export) {
          # Search common local locations (same logic as grid_react_df)
          search_dirs <- c(tempdir(), "/data/fasta", "/quobyte/proteomics-grp/de-limp/fasta")
          if (!is.null(values$diann_fasta_files)) {
            search_dirs <- c(dirname(values$diann_fasta_files), search_dirs)
          }
          for (d in unique(search_dirs)) {
            if (!dir.exists(d)) next
            gmaps <- list.files(d, pattern = "gene_map\\.tsv$", full.names = TRUE)
            if (length(gmaps) > 0) {
              export_gene_map <- tryCatch(read.delim(gmaps[1], stringsAsFactors = FALSE), error = function(e) NULL)
              if (!is.null(export_gene_map) && nrow(export_gene_map) > 0) {
                message("[Explorer Export] Loaded gene map: ", gmaps[1], " (", nrow(export_gene_map), " entries)")
                break
              }
            }
          }

          # SSH fallback
          if (is.null(export_gene_map) && isTRUE(values$ssh_connected)) {
            tryCatch({
              cfg <- list(host = isolate(input$ssh_host), user = isolate(input$ssh_user),
                          port = isolate(input$ssh_port) %||% 22L,
                          key_path = isolate(input$ssh_key_path))
              remote_result <- ssh_exec(cfg,
                "ls /quobyte/proteomics-grp/de-limp/fasta/*gene_map.tsv 2>/dev/null | head -1",
                timeout = 10)
              if (remote_result$status == 0 && length(remote_result$stdout) > 0 &&
                  nzchar(trimws(remote_result$stdout[1]))) {
                remote_path <- trimws(remote_result$stdout[1])
                local_path <- file.path(tempdir(), basename(remote_path))
                dl <- scp_download(cfg, remote_path, local_path)
                if (dl$status == 0 && file.exists(local_path)) {
                  export_gene_map <- tryCatch(read.delim(local_path, stringsAsFactors = FALSE), error = function(e) NULL)
                  message("[Explorer Export] Downloaded gene map via SSH: ", nrow(export_gene_map), " entries")
                }
              }
            }, error = function(e) message("[Explorer Export] SSH gene map download failed: ", e$message))
          }

          # Deduplicate
          if (!is.null(export_gene_map) && nrow(export_gene_map) > 0 && "gene_symbol" %in% colnames(export_gene_map)) {
            export_gene_map <- export_gene_map[!duplicated(export_gene_map$accession), ]
          } else {
            export_gene_map <- NULL
          }
        }

        # Helper: resolve gene label using gene_map.tsv (NCBI) or explorer_gene_label (UniProt)
        resolve_gene <- function(protein_id) {
          if (!is.null(export_gene_map)) {
            acc <- sub(";.*", "", protein_id)  # first accession from semicolon group
            acc <- sub("^Cont_", "", acc)
            idx <- match(acc, export_gene_map$accession)
            if (!is.na(idx)) {
              gs <- export_gene_map$gene_symbol[idx]
              if (!is.na(gs) && nzchar(gs)) return(gs)
            }
          }
          explorer_gene_label(protein_id, genes_df)
        }

        # --- 1. Expression matrix CSV ---
        incProgress(0.1, detail = "Expression matrix...")
        expr_df <- as.data.frame(mat)
        # Add gene symbols
        expr_df$Protein.Group <- rownames(mat)
        expr_df$Gene <- vapply(rownames(mat), resolve_gene, character(1))
        # Count how many NCBI proteins were successfully mapped via gene_map.tsv
        if (!is.null(export_gene_map)) {
          non_contam_genes <- expr_df$Gene[!grepl("^Cont_", expr_df$Protein.Group)]
          non_contam_pids <- expr_df$Protein.Group[!grepl("^Cont_", expr_df$Protein.Group)]
          n_gene_mapped <- sum(non_contam_genes != non_contam_pids &
            non_contam_genes != substr(non_contam_pids, 1, 15))
        }
        if (!is.null(genes_df) && "Protein.Names" %in% colnames(genes_df)) {
          name_idx <- match(rownames(mat), rownames(genes_df))
          expr_df$Protein.Name <- ifelse(is.na(name_idx), "", genes_df$Protein.Names[name_idx])
        } else {
          expr_df$Protein.Name <- ""
        }
        # Reorder: ID columns first, then samples
        id_cols <- c("Protein.Group", "Gene", "Protein.Name")
        expr_df <- expr_df[, c(id_cols, setdiff(colnames(expr_df), id_cols))]
        expr_file <- file.path(tmp_dir, "expression_matrix.csv")
        write.csv(expr_df, expr_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, expr_file)

        # --- 1b. Detection matrix from raw precursor data (shows real missing values) ---
        if (!is.null(values$raw_data) && !is.null(values$raw_data$E)) {
          tryCatch({
            raw_mat <- values$raw_data$E
            raw_genes <- values$raw_data$genes
            if (!is.null(raw_genes) && "Protein.Group" %in% colnames(raw_genes)) {
              # Per protein group: count detected precursors per sample
              pg <- raw_genes$Protein.Group
              det_counts <- do.call(rbind, lapply(unique(pg), function(p) {
                rows <- which(pg == p)
                sub_mat <- raw_mat[rows, , drop = FALSE]
                detected <- colSums(!is.na(sub_mat) & is.finite(sub_mat))
                total_prec <- nrow(sub_mat)
                c(Protein.Group = p, Total_Precursors = total_prec,
                  setNames(as.list(detected), paste0("Detected_", colnames(raw_mat))))
              }))
              det_df <- as.data.frame(det_counts, stringsAsFactors = FALSE)
              # Add gene symbols
              det_df$Gene <- vapply(det_df$Protein.Group, function(pid) resolve_gene(pid), character(1))
              det_df <- det_df[, c("Protein.Group", "Gene", "Total_Precursors",
                grep("^Detected_", colnames(det_df), value = TRUE))]
              det_file <- file.path(tmp_dir, "detection_matrix.csv")
              write.csv(det_df, det_file, row.names = FALSE)
              files_to_zip <- c(files_to_zip, det_file)
            }
          }, error = function(e) message("[Export] Detection matrix failed: ", e$message))
        }

        # --- 1c. DIA-NN protein group matrix (with real missing values) ---
        # Try to include pg_matrix.tsv from the search output directory
        tryCatch({
          ss <- values$diann_search_settings
          if (!is.null(ss) && !is.null(ss$output_dir)) {
            od <- translate_storage_path(ss$output_dir, to = "hpc")
            pg_matrix_remote <- file.path(od, "report.pg_matrix.tsv")

            # Check locally first, then via SSH
            pg_local <- file.path(ss$output_dir, "report.pg_matrix.tsv")
            if (file.exists(pg_local)) {
              pg_file <- file.path(tmp_dir, "diann_pg_matrix.tsv")
              file.copy(pg_local, pg_file)
              files_to_zip <- c(files_to_zip, pg_file)
              message("[Export] Included local pg_matrix.tsv")
            } else if (isTRUE(values$ssh_connected)) {
              cfg <- list(host = isolate(input$ssh_host), user = isolate(input$ssh_user),
                          port = isolate(input$ssh_port) %||% 22L,
                          key_path = isolate(input$ssh_key_path))
              pg_file <- file.path(tmp_dir, "diann_pg_matrix.tsv")
              dl <- scp_download(cfg, pg_matrix_remote, pg_file)
              if (dl$status == 0 && file.exists(pg_file)) {
                files_to_zip <- c(files_to_zip, pg_file)
                message("[Export] Downloaded pg_matrix.tsv via SSH")
              }
            }
          }
        }, error = function(e) message("[Export] pg_matrix.tsv not available: ", e$message))

        # --- 1d. Data quality summary (per-sample protein counts + missingness) ---
        tryCatch({
          # Use pg_matrix if we downloaded it, otherwise derive from raw_data
          pg_file_path <- file.path(tmp_dir, "diann_pg_matrix.tsv")
          if (file.exists(pg_file_path)) {
            pg <- read.delim(pg_file_path, stringsAsFactors = FALSE, check.names = FALSE)
            # Intensity columns are after the annotation columns
            annot_cols <- c("Protein.Group", "Protein.Names", "Genes",
                            "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences")
            int_cols <- setdiff(colnames(pg), annot_cols)

            if (length(int_cols) > 0) {
              pg_mat <- as.matrix(pg[, int_cols])
              # 0 = not detected in DIA-NN pg_matrix
              detected <- colSums(pg_mat > 0, na.rm = TRUE)
              total_pg <- nrow(pg_mat)
              contam_count <- sum(grepl("^Cont_", pg$Protein.Group))

              quality_df <- data.frame(
                Sample = int_cols,
                Proteins_Detected = detected,
                Total_Protein_Groups = total_pg,
                Pct_Detected = round(100 * detected / total_pg, 1),
                Missing = total_pg - detected,
                Pct_Missing = round(100 * (total_pg - detected) / total_pg, 1),
                Contaminant_Proteins = contam_count,
                stringsAsFactors = FALSE
              )
              # Add short sample labels from metadata
              if (!is.null(values$metadata)) {
                quality_df$Group <- values$metadata$Group[match(quality_df$Sample,
                  values$metadata$File.Name)]
                quality_df$Sample_ID <- paste0("S", values$metadata$ID[match(quality_df$Sample,
                  values$metadata$File.Name)])
              }
              quality_file <- file.path(tmp_dir, "data_quality_summary.csv")
              write.csv(quality_df, quality_file, row.names = FALSE)
              files_to_zip <- c(files_to_zip, quality_file)
            }
          }
        }, error = function(e) message("[Export] data quality summary failed: ", e$message))

        # --- 2. Quartile profiles CSV (recompute without contaminant exclusion) ---
        incProgress(0.2, detail = "Quartile profiles...")
        avg_intensity <- rowMeans(mat, na.rm = TRUE)
        ranks <- rank(-avg_intensity, ties.method = "first")
        n <- length(ranks)
        avg_quartile <- ifelse(ranks <= n/4, "Q1",
          ifelse(ranks <= n/2, "Q2",
            ifelse(ranks <= 3*n/4, "Q3", "Q4")))
        names(avg_quartile) <- names(avg_intensity)

        # Per-sample quartile assignment for ALL proteins
        all_sample_q <- matrix(NA_integer_, nrow = nrow(mat), ncol = ncol(mat))
        rownames(all_sample_q) <- rownames(mat)
        colnames(all_sample_q) <- colnames(mat)
        for (j in seq_len(ncol(mat))) {
          col_vals <- mat[, j]
          col_ranks <- rank(-col_vals, ties.method = "first")
          col_n <- length(col_ranks)
          all_sample_q[, j] <- ifelse(col_ranks <= col_n/4, 1L,
            ifelse(col_ranks <= col_n/2, 2L,
              ifelse(col_ranks <= 3*col_n/4, 3L, 4L)))
        }

        min_q <- apply(all_sample_q, 1, min, na.rm = TRUE)
        max_q <- apply(all_sample_q, 1, max, na.rm = TRUE)
        q_range <- max_q - min_q

        quartile_df <- data.frame(
          Protein.Group = rownames(mat),
          Gene = vapply(rownames(mat), resolve_gene, character(1)),
          stringsAsFactors = FALSE
        )
        if (!is.null(genes_df) && "Protein.Names" %in% colnames(genes_df)) {
          name_idx <- match(rownames(mat), rownames(genes_df))
          quartile_df$Protein.Name <- ifelse(is.na(name_idx), "", genes_df$Protein.Names[name_idx])
        } else {
          quartile_df$Protein.Name <- ""
        }
        quartile_df$Avg_Intensity <- round(avg_intensity, 2)
        quartile_df$Avg_Quartile <- avg_quartile

        # Per-sample quartile columns
        for (j in seq_len(ncol(mat))) {
          col_name <- paste0("Q_", colnames(mat)[j])
          quartile_df[[col_name]] <- paste0("Q", all_sample_q[, j])
        }
        quartile_df$Quartile_Range <- q_range
        quartile_df$Is_Contaminant <- grepl("^Cont_", rownames(mat))

        quartile_df <- quartile_df[order(-quartile_df$Avg_Intensity), ]
        quartile_file <- file.path(tmp_dir, "quartile_profiles.csv")
        write.csv(quartile_df, quartile_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, quartile_file)

        # --- 3. Variable proteins CSV ---
        incProgress(0.3, detail = "Variable proteins...")
        variable_idx <- which(q_range >= 2)
        n_variable <- length(variable_idx)
        if (n_variable > 0) {
          var_df <- quartile_df[quartile_df$Protein.Group %in% rownames(mat)[variable_idx], ]
          var_df <- var_df[order(-var_df$Quartile_Range, -var_df$Avg_Intensity), ]
          var_file <- file.path(tmp_dir, "variable_proteins.csv")
          write.csv(var_df, var_file, row.names = FALSE)
          files_to_zip <- c(files_to_zip, var_file)
        } else {
          n_variable <- 0
          # Write empty CSV with header
          var_file <- file.path(tmp_dir, "variable_proteins.csv")
          write.csv(quartile_df[0, ], var_file, row.names = FALSE)
          files_to_zip <- c(files_to_zip, var_file)
        }

        # --- 4. Sample metadata CSV ---
        incProgress(0.4, detail = "Sample metadata...")
        if (!is.null(values$metadata)) {
          meta_file <- file.path(tmp_dir, "sample_metadata.csv")
          write.csv(values$metadata, meta_file, row.names = FALSE)
          files_to_zip <- c(files_to_zip, meta_file)
        }

        # --- 5. Contaminant summary CSV ---
        incProgress(0.45, detail = "Contaminant summary...")
        contam_mask <- grepl("^Cont_", rownames(mat))
        n_contam <- sum(contam_mask)
        contam_note <- "No contaminant proteins detected in this dataset."
        if (n_contam > 0) {
          contam_mat <- mat[contam_mask, , drop = FALSE]
          contam_summary <- data.frame(
            Protein.Group = rownames(contam_mat),
            Gene = vapply(rownames(contam_mat), resolve_gene, character(1)),
            Avg_Intensity = round(rowMeans(contam_mat, na.rm = TRUE), 2),
            stringsAsFactors = FALSE
          )
          # Per-sample intensities
          for (j in seq_len(ncol(contam_mat))) {
            contam_summary[[colnames(contam_mat)[j]]] <- round(contam_mat[, j], 2)
          }
          contam_summary <- contam_summary[order(-contam_summary$Avg_Intensity), ]
          contam_file <- file.path(tmp_dir, "contaminant_summary.csv")
          write.csv(contam_summary, contam_file, row.names = FALSE)
          files_to_zip <- c(files_to_zip, contam_file)
          contam_note <- paste0(n_contam, " contaminant proteins detected (",
            round(n_contam / n_proteins * 100, 1), "% of total).")
        }

        # --- 6. Session RDS ---
        incProgress(0.5, detail = "Saving session state...")
        tryCatch({
          session_data <- list(
            raw_data = values$raw_data, metadata = values$metadata,
            y_protein = values$y_protein, design = values$design,
            qc_stats = values$qc_stats,
            repro_log = values$repro_log,
            instrument_metadata = values$instrument_metadata,
            diann_search_settings = values$diann_search_settings,
            saved_at = Sys.time(),
            app_version = paste0("DE-LIMP v", values$app_version)
          )
          rds_file <- file.path(tmp_dir, "session.rds")
          saveRDS(session_data, rds_file)
          files_to_zip <- c(files_to_zip, rds_file)
        }, error = function(e) message("[DE-LIMP] Explorer export: RDS save error: ", e$message))

        # --- 7. Methods text ---
        incProgress(0.55, detail = "Methods...")
        tryCatch({
          params <- c(
            "DE-LIMP Data Explorer Export",
            paste0("Export date: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
            paste0("App version: DE-LIMP v", values$app_version),
            paste0("R version: ", R.version.string),
            "",
            "MODE: Exploratory analysis (no differential expression)",
            paste0("Total proteins: ", n_proteins),
            paste0("Total samples: ", n_samples),
            paste0("Variable proteins (quartile range >= 2): ", n_variable),
            paste0("Contaminant proteins: ", n_contam),
            ""
          )
          # Groups
          if (!is.null(values$metadata)) {
            grp_counts <- table(values$metadata$Group[values$metadata$Group != ""])
            if (length(grp_counts) > 0) {
              params <- c(params, "GROUPS:",
                paste0("  ", names(grp_counts), ": n=", grp_counts), "")
            }
          }
          # DIA-NN search settings
          ss <- values$diann_search_settings
          if (!is.null(ss) && is.list(ss)) {
            sp <- ss$search_params
            params <- c(params, "DIA-NN SEARCH SETTINGS:",
              if (!is.null(ss$diann_version) && nzchar(ss$diann_version))
                paste0("  DIA-NN version: ", ss$diann_version) else NULL,
              paste0("  FASTA: ", paste(basename(ss$fasta_files), collapse = ", ")),
              paste0("  Enzyme: ", sp$enzyme),
              paste0("  MBR: ", if (isTRUE(sp$mbr)) "enabled" else "disabled"),
              if (!is.null(sp$mass_acc) && !is.na(sp$mass_acc))
                paste0("  Mass accuracy (MS2): ", sp$mass_acc, " ppm") else NULL,
              if (!is.null(ss$normalization) && nzchar(ss$normalization))
                paste0("  Normalization: ", ss$normalization) else NULL,
              "")
          }
          # Instrument metadata
          if (!is.null(values$instrument_metadata)) {
            meta <- values$instrument_metadata
            params <- c(params, "INSTRUMENT:",
              if (!is.null(meta$instrument_model)) paste0("  Model: ", meta$instrument_model) else NULL,
              if (!is.null(meta$lc_method_name)) paste0("  LC method: ", meta$lc_method_name) else NULL,
              if (!is.null(meta$gradient_length_min)) paste0("  Gradient: ", meta$gradient_length_min, " min") else NULL,
              "")
          }
          # Package versions
          params <- c(params, "PACKAGE VERSIONS:",
            paste0("  limpa: ", tryCatch(as.character(packageVersion("limpa")), error = function(e) "unknown")),
            paste0("  DE-LIMP: v", values$app_version %||% "unknown"))

          methods_file <- file.path(tmp_dir, "methods.txt")
          writeLines(params, methods_file)
          files_to_zip <- c(files_to_zip, methods_file)
        }, error = function(e) message("[DE-LIMP] Explorer export: methods error: ", e$message))

        # --- 7b. search_info.md (DIA-NN search parameters and metadata) ---
        tryCatch({
          ss <- values$diann_search_settings
          if (!is.null(ss) && !is.null(ss$output_dir)) {
            od <- translate_storage_path(ss$output_dir, to = "hpc")
            si_remote <- file.path(od, "search_info.md")
            si_local <- file.path(ss$output_dir, "search_info.md")
            si_file <- file.path(tmp_dir, "search_info.md")

            if (file.exists(si_local)) {
              file.copy(si_local, si_file)
              files_to_zip <- c(files_to_zip, si_file)
            } else if (isTRUE(values$ssh_connected)) {
              cfg <- list(host = isolate(input$ssh_host), user = isolate(input$ssh_user),
                          port = isolate(input$ssh_port) %||% 22L,
                          key_path = isolate(input$ssh_key_path))
              dl <- scp_download(cfg, si_remote, si_file)
              if (dl$status == 0 && file.exists(si_file)) {
                files_to_zip <- c(files_to_zip, si_file)
              }
            }
          }
        }, error = function(e) message("[Export] search_info.md not available: ", e$message))

        # --- 8. Reproducibility log ---
        incProgress(0.6, detail = "Reproducibility log...")
        if (!is.null(values$repro_log) && length(values$repro_log) > 0) {
          repro_file <- file.path(tmp_dir, "reproducibility_log.R")
          log_content <- paste(values$repro_log, collapse = "\n")
          writeLines(log_content, repro_file)
          files_to_zip <- c(files_to_zip, repro_file)
        }

        # --- 9. PROMPT.md ---
        incProgress(0.7, detail = "Building prompt...")

        # Detect organism
        organism_info <- tryCatch({
          org_db <- detect_organism_db(rownames(mat))
          org_map <- c("org.Hs.eg.db" = "Human (Homo sapiens)",
            "org.Mm.eg.db" = "Mouse (Mus musculus)",
            "org.Rn.eg.db" = "Rat (Rattus norvegicus)",
            "org.Bt.eg.db" = "Bovine (Bos taurus)",
            "org.Cf.eg.db" = "Dog (Canis lupus familiaris)",
            "org.Gg.eg.db" = "Chicken (Gallus gallus)",
            "org.Dm.eg.db" = "Fruit fly (Drosophila melanogaster)",
            "org.Ce.eg.db" = "C. elegans",
            "org.Dr.eg.db" = "Zebrafish (Danio rerio)",
            "org.Sc.sgd.db" = "Yeast (Saccharomyces cerevisiae)",
            "org.At.tair.db" = "Arabidopsis thaliana",
            "org.Ss.eg.db" = "Pig (Sus scrofa)")
          paste0("Organism: ", org_map[org_db] %||% "Unknown", " (OrgDb: ", org_db, ")")
        }, error = function(e) "Organism: Unknown (detection failed)")

        # Group info
        group_info <- "No group assignments available."
        if (!is.null(values$metadata)) {
          grp_counts <- table(values$metadata$Group[values$metadata$Group != ""])
          if (length(grp_counts) > 0) {
            group_lines <- paste0("- ", names(grp_counts), ": n=", grp_counts)
            group_info <- paste(group_lines, collapse = "\n")
          }
        }

        prompt_text <- paste0(
'# DE-LIMP Data Exploration Analysis

## Context
This is a proteomics dataset analyzed with DE-LIMP v', values$app_version, '. The dataset has ', format(n_proteins, big.mark = ","),
' proteins quantified across ', n_samples, ' samples. **No differential expression analysis was performed** (either because there are no replicates, or the user chose exploratory analysis only).

## Your Task
You are a proteomics bioinformatics expert. Analyze this dataset and provide biological insights WITHOUT relying on statistical significance (there are no p-values or fold changes).

## Data Files

**IMPORTANT**: All CSV files use the **Gene** column for gene symbols (not Protein.Group, which contains raw database accessions like XP_ or UniProt IDs). Always use the Gene column for biological interpretation.

| File | Description | Key Columns |
|------|-------------|-------------|
| `expression_matrix.csv` | Log2 protein intensities (', format(n_proteins, big.mark = ","), ' proteins x ', n_samples, ' samples) | **Gene** (symbol), Protein.Name (description), then sample intensity columns |
| `diann_pg_matrix.tsv` | DIA-NN protein group matrix with **real missing values** (0 = not detected). This is BEFORE maxLFQ — use this to see which proteins were directly quantified vs inferred. Compare with expression_matrix.csv to understand data completeness. | Protein.Group, Genes, sample intensity columns (0 = missing) |
| `data_quality_summary.csv` | Per-sample data quality: proteins detected, % missing, contaminant counts. Use this to assess sample quality and compare completeness across runs | Sample, Proteins_Detected, Pct_Missing, Group |
| `detection_matrix.csv` | Per-protein precursor detection counts from DIA-NN (BEFORE maxLFQ). Shows how many precursors were directly detected per sample — use this to assess data completeness | **Gene**, Total_Precursors, Detected_SampleN columns |
| `quartile_profiles.csv` | Per-sample quartile assignments (Q1=top 25%, Q4=bottom 25%) | **Gene**, Avg_Quartile, per-sample Q columns, Quartile_Range |
| `variable_proteins.csv` | ', n_variable, ' proteins shifting 2+ quartiles across samples | **Gene**, Avg_Intensity, Quartile_Range |
| `sample_metadata.csv` | Sample groups and identifiers | |
| `contaminant_summary.csv` | Contaminant protein statistics | |
| `session.rds` | Full DE-LIMP session state (reload via DE-LIMP > Load Session) | |
| `methods.txt` | Pipeline parameters, normalization, app version | |
| `reproducibility_log.R` | R code log recording every analysis step | |

## Analysis Requested

### 1. Quartile-Based Biological Insights
For each intensity quartile (Q1 through Q4), analyze the top proteins:
- **Q1 (Most Abundant)**: What biological processes dominate? Are these expected housekeeping/structural proteins? Any surprises?
- **Q2-Q3 (Mid-Range)**: What functional categories are enriched? Are signaling/regulatory proteins concentrated here?
- **Q4 (Low Abundance)**: Are there transcription factors, kinases, or other low-abundance regulators? What biological signals might be hiding in this quartile?

### 2. Variable Protein Analysis
The `variable_proteins.csv` contains ', n_variable, ' proteins whose intensity rank shifts dramatically across samples (e.g., top 25% in one sample, bottom 50% in another). For these proteins:
- What biological processes do they represent?
- Are any known biomarkers or disease-associated proteins?
- Do the variable proteins cluster into functional groups (e.g., immune response, metabolism, stress response)?
- Which variable proteins would you prioritize for follow-up experiments?

### 3. Sample Comparison
Compare the protein profiles across samples:
- Which samples are most similar/different based on protein abundance patterns?
- Are there sample-specific protein signatures?
- Do any samples show signs of technical issues (unusual contaminant levels, missing proteins)?

### 4. Contaminant Assessment
Review the contaminant data:
- ', contam_note, '
- Are there sample-specific contamination patterns suggesting prep issues?
- Are keratin levels consistent or variable across samples?

### 5. Biological Hypotheses
Based on all the data, propose 3-5 testable biological hypotheses that could be validated with follow-up experiments (e.g., with replicates for proper statistical testing).

## Organism & Database
', organism_info, '

## Sample Groups
', group_info, '

## Important Notes — Read Before Analysis
- All intensity values are **log2-transformed**
- The expression matrix is **complete** (no missing values). This is NOT from imputation — it is a mathematical property of the **maxLFQ algorithm** (via limpa DPC-Quant), which computes protein-level quantities from pairwise precursor ratios. Do NOT describe this as "MBR imputation" or "gap-filling" — maxLFQ produces valid quantitative estimates for all proteins across all samples.
- **DIA-NN MBR** (Match Between Runs) transfers peptide IDs across runs at the precursor level, but the complete protein matrix comes from maxLFQ aggregation, not MBR.
- Quartile assignments are computed independently per sample (a protein can be Q1 in one sample and Q3 in another)
- Variable proteins are candidates, not statistically validated — they need replicated experiments for confirmation
- Contaminant proteins are prefixed with "Cont_"
',
        if (!is.null(export_gene_map) && n_gene_mapped > 0) {
          paste0('- **Gene symbols mapped from NCBI RefSeq accessions** using gene_map.tsv (', n_gene_mapped, ' of ',
            length(non_contam_ids), ' non-contaminant proteins mapped to gene symbols via NCBI E-utilities lookup). ',
            'Unmapped proteins retain their accession IDs in the Gene column.\n')
        } else if (is_ncbi_export && is.null(export_gene_map)) {
          '- **WARNING: NCBI RefSeq accessions detected but gene_map.tsv not found.** Gene column contains truncated protein descriptions, not proper gene symbols. For accurate gene symbols, ensure gene_map.tsv is available in the FASTA directory.\n'
        } else { '' },
'
')

        prompt_file <- file.path(tmp_dir, "PROMPT.md")
        writeLines(prompt_text, prompt_file)
        files_to_zip <- c(files_to_zip, prompt_file)

        # --- Create ZIP ---
        incProgress(0.9, detail = "Creating ZIP...")
        # Use basename so zip doesn't include full tmp path
        old_wd <- setwd(tmp_dir)
        on.exit(setwd(old_wd), add = TRUE)
        zip(file, basename(files_to_zip))

        message("[DE-LIMP] Explorer Claude export complete: ", length(files_to_zip), " files")
      })
      }, error = function(e) {
        message("[DE-LIMP] Explorer Claude export FAILED: ", e$message)
        showNotification(paste("Export error:", e$message), type = "error", duration = 15)
      })
    },
    contentType = "application/zip"
  )

}
