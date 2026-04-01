server_ai <- function(input, output, session, values) {

  # --- Helper: Build DE data context for AI prompts ---
  build_ai_data_context <- function() {
    req(values$fit, values$y_protein)

    all_contrasts <- colnames(values$fit$contrasts)
    n_contrasts <- length(all_contrasts)
    top_n <- if (n_contrasts <= 3) 30 else if (n_contrasts <= 6) 20 else 10

    # Gene mapping
    first_tt <- topTable(values$fit, coef = all_contrasts[1], number = Inf) %>% as.data.frame()
    if (!"Protein.Group" %in% colnames(first_tt)) first_tt <- first_tt %>% rownames_to_column("Protein.Group")
    first_tt$Accession <- str_split_fixed(first_tt$Protein.Group, "[; ]", 2)[,1]
    org_db_name <- detect_organism_db(first_tt$Protein.Group)

    id_map <- tryCatch({
      if (!requireNamespace(org_db_name, quietly = TRUE)) BiocManager::install(org_db_name, ask = FALSE)
      library(org_db_name, character.only = TRUE)
      db_obj <- get(org_db_name)
      AnnotationDbi::select(db_obj, keys = first_tt$Accession, columns = c("SYMBOL"), keytype = "UNIPROT") %>%
        dplyr::rename(Accession = UNIPROT, Gene = SYMBOL) %>% distinct(Accession, .keep_all = TRUE)
    }, error = function(e) data.frame(Accession = first_tt$Accession, Gene = first_tt$Accession))

    # Per-contrast DE summaries
    contrast_texts <- list()
    all_sig_proteins <- list()

    for (i in seq_along(all_contrasts)) {
      cname <- all_contrasts[i]
      tt <- topTable(values$fit, coef = cname, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(tt)) tt <- tt %>% rownames_to_column("Protein.Group")
      tt$Accession <- str_split_fixed(tt$Protein.Group, "[; ]", 2)[,1]
      tt <- left_join(tt, id_map, by = "Accession")
      tt$Gene[is.na(tt$Gene)] <- tt$Accession[is.na(tt$Gene)]

      sig <- tt %>% filter(adj.P.Val < 0.05)
      n_up <- sum(sig$logFC > 0)
      n_down <- sum(sig$logFC < 0)

      if (nrow(sig) > 0) {
        for (pid in sig$Protein.Group) {
          if (is.null(all_sig_proteins[[pid]])) all_sig_proteins[[pid]] <- list()
          row <- sig[sig$Protein.Group == pid, ]
          all_sig_proteins[[pid]][[cname]] <- list(
            gene = row$Gene[1], logFC = round(row$logFC[1], 3), pval = round(row$adj.P.Val[1], 4)
          )
        }
      }

      top_hits <- sig %>% arrange(adj.P.Val) %>% head(top_n) %>%
        dplyr::select(Gene, logFC, adj.P.Val) %>%
        mutate(across(where(is.numeric), ~round(.x, 3)))

      top_text <- paste(capture.output(print(as.data.frame(top_hits))), collapse = "\n")

      contrast_texts[[cname]] <- paste0(
        "### ", cname, "\n",
        "Significant proteins: ", nrow(sig), " (", n_up, " up, ", n_down, " down)\n\n",
        "Top ", min(top_n, nrow(top_hits)), " by significance:\n", top_text
      )
    }

    # Cross-contrast proteins
    multi_contrast <- names(all_sig_proteins)[sapply(all_sig_proteins, length) >= 2]
    cross_text <- if (length(multi_contrast) > 0) {
      cross_df <- do.call(rbind, lapply(head(multi_contrast, 10), function(pid) {
        info <- all_sig_proteins[[pid]]
        gene <- info[[1]]$gene
        contrasts_str <- paste(names(info), collapse = ", ")
        fc_str <- paste(sapply(names(info), function(cn) {
          paste0(cn, ": ", sprintf("%+.2f", info[[cn]]$logFC))
        }), collapse = "; ")
        data.frame(Gene = gene, N_Comparisons = length(info), Contrasts = contrasts_str, LogFC = fc_str)
      }))
      cross_df <- cross_df[order(-cross_df$N_Comparisons), ]
      paste(capture.output(print(as.data.frame(cross_df), row.names = FALSE)), collapse = "\n")
    } else {
      "No proteins were significant in more than one comparison."
    }

    # Stable biomarkers (lowest CV)
    # NOTE: Do NOT use return() inside tryCatch — it exits the enclosing function,
    # not just the tryCatch block. Use if/else instead. (Same gotcha as withProgress.)
    stable_prots_text <- tryCatch({
      all_sig_pids <- names(all_sig_proteins)
      valid_pids <- intersect(all_sig_pids, rownames(values$y_protein$E))
      if (length(valid_pids) == 0) {
        "No significant proteins to assess for stability."
      } else {
        raw_exprs <- values$y_protein$E[valid_pids, , drop = FALSE]
        linear_exprs <- 2^raw_exprs
        cv_list <- list()

        for (g in unique(values$metadata$Group)) {
          if (g == "") next
          files_in_group <- values$metadata$File.Name[values$metadata$Group == g]
          group_cols <- intersect(colnames(linear_exprs), files_in_group)
          if (length(group_cols) > 1) {
            group_data <- linear_exprs[, group_cols, drop = FALSE]
            cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
          }
        }

        if (length(cv_list) == 0) {
          "Could not calculate CVs (not enough replicates)."
        } else {
          cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
          cv_df$Avg_CV <- rowMeans(cv_df[, grep("^CV_", colnames(cv_df)), drop = FALSE], na.rm = TRUE)

          stable_df <- cv_df %>% arrange(Avg_CV) %>% head(5)
          stable_df$Gene <- sapply(stable_df$Protein.Group, function(pid) {
            info <- all_sig_proteins[[pid]]
            if (!is.null(info)) info[[1]]$gene else pid
          })
          stable_df$Significant_In <- sapply(stable_df$Protein.Group, function(pid) {
            info <- all_sig_proteins[[pid]]
            if (!is.null(info)) paste(names(info), collapse = "; ") else ""
          })

          out <- stable_df %>% dplyr::select(Gene, Avg_CV, Significant_In) %>%
            mutate(Avg_CV = round(Avg_CV, 2))
          paste(capture.output(print(as.data.frame(out), row.names = FALSE)), collapse = "\n")
        }
      }
    }, error = function(e) "Could not calculate stable proteins.")

    all_contrast_text <- paste(contrast_texts, collapse = "\n\n")

    list(
      n_contrasts = n_contrasts,
      contrast_text = all_contrast_text,
      cross_text = cross_text,
      stable_prots_text = stable_prots_text
    )
  }

  # --- AI SUMMARY (Data Overview Tab) — Analyzes ALL contrasts ---
  observeEvent(input$generate_ai_summary_overview, {
    req(values$fit, values$y_protein, input$user_api_key)

    withProgress(message = "Generating AI Summary...", value = 0, {
      incProgress(0.1, detail = "Gathering DE data across all comparisons...")

      ctx <- build_ai_data_context()

      incProgress(0.7, detail = "Constructing prompt...")

      # --- Build the full prompt ---
      system_prompt <- paste0(
        "You are a senior proteomics and systems biology consultant. Write a comprehensive ",
        "analysis of the differential expression results across ALL comparisons below.\n\n",
        "Structure your response with these markdown sections:\n\n",
        "## Overview\n",
        "Number of comparisons analyzed, total significant proteins per comparison (up/down split). ",
        "Overall assessment of the experiment's quality and scope.\n\n",
        "## Key Findings Per Comparison\n",
        "For each comparison: highlight the top upregulated and downregulated proteins by fold-change ",
        "(use gene names). Note any comparison with unusually few or many significant hits.\n\n",
        "## Cross-Comparison Biomarkers\n",
        "Proteins significant in multiple comparisons are highest-confidence candidates. ",
        "Discuss consistency of direction (always up, always down, or mixed across comparisons).\n\n",
        "## High-Confidence Biomarker Insights\n",
        "For the most stable proteins (lowest coefficient of variation): discuss their known biological functions, ",
        "pathway involvement, and disease associations where you recognize the gene name. ",
        "Assess their potential as reliable biomarkers based on the combination of low CV, ",
        "significant p-value, and meaningful fold-change.\n\n",
        "## Biological Interpretation\n",
        "Suggest what biological processes or pathways may be affected based on the protein lists. ",
        "Note any well-known protein families, complexes, or signaling cascades represented. ",
        "If the data suggests a clear biological narrative, describe it.\n\n",
        "Use markdown formatting with headers. Be scientific but accessible."
      )

      final_prompt <- paste0(
        system_prompt,
        "\n\n--- DATA FOR ANALYSIS ---\n\n",
        "Number of comparisons: ", ctx$n_contrasts, "\n\n",
        ctx$contrast_text, "\n\n",
        "--- CROSS-COMPARISON PROTEINS (significant in >= 2 comparisons) ---\n",
        ctx$cross_text, "\n\n",
        "--- MOST STABLE SIGNIFICANT PROTEINS (lowest CV across replicates) ---\n",
        ctx$stable_prots_text
      )

      message(sprintf("[DE-LIMP] AI Summary prompt: %d characters, %d contrasts", nchar(final_prompt), ctx$n_contrasts))

      incProgress(0.8, detail = "Asking AI...")
      ai_summary <- ask_gemini_text_chat(final_prompt, input$user_api_key, input$model_name)

      # Store for export and show download buttons
      values$ai_summary_text <- ai_summary
      shinyjs::show("download_ai_summary_html")

      # Render the summary to the output area
      output$ai_summary_output <- renderUI({
        div(style = "background-color: #ffffff; padding: 20px; border: 1px solid #dee2e6; border-radius: 8px;",
          tags$h5(icon("check-circle"), " Analysis Complete", style = "color: #28a745; margin-bottom: 15px;"),
          HTML(markdown::markdownToHTML(text = ai_summary, fragment.only = TRUE))
        )
      })
    })
  })

  # --- AI Summary Markdown Export ---
  output$download_ai_summary_html <- downloadHandler(
    filename = function() {
      paste0("AI_Analysis_Report_", format(Sys.time(), "%Y%m%d_%H%M"), ".md")
    },
    content = function(file) {
      req(values$ai_summary_text)

      md_doc <- paste0(
        "# DE-LIMP AI Analysis Report\n\n",
        "*Generated: ", format(Sys.time(), "%B %d, %Y at %I:%M %p"), "*\n\n",
        "---\n\n",
        values$ai_summary_text, "\n\n",
        "---\n\n",
        "*Generated by [DE-LIMP Proteomics](https://github.com/bsphinney/DE-LIMP)*\n"
      )

      writeLines(md_doc, file)
    }
  )

  # --- Export Prompt for Claude (downloads .zip with prompt + full data CSVs) ---
  # Shared content function for both download buttons (AI Summary tab + AI Chat tab)
  claude_export_content <- function(file) {
      req(values$fit, values$y_protein)

      tryCatch({
      withProgress(message = "Building export...", value = 0, {
        tmp_dir <- tempdir()
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
        files_to_zip <- character(0)

        message("[DE-LIMP] Claude export: starting...")
        incProgress(0.1, detail = "Gathering DE data...")
        ctx <- tryCatch(build_ai_data_context(), error = function(e) {
          message("[DE-LIMP] Claude export: build_ai_data_context FAILED: ", e$message)
          list(n_contrasts = length(colnames(values$fit$contrasts)),
               contrast_text = "(Error building DE summary)",
               cross_text = "(Error)", stable_prots_text = "(Error)")
        })
        message("[DE-LIMP] Claude export: ctx OK")

        # --- 1. Full DE results CSV (all proteins, all contrasts) ---
        incProgress(0.3, detail = "Exporting full DE results...")
        all_contrasts <- colnames(values$fit$contrasts)
        full_results <- do.call(rbind, lapply(all_contrasts, function(cname) {
          tt <- topTable(values$fit, coef = cname, number = Inf) %>% as.data.frame()
          if (!"Protein.Group" %in% colnames(tt)) tt <- tt %>% rownames_to_column("Protein.Group")
          tt$Contrast <- cname
          tt
        }))
        results_file <- file.path(tmp_dir, "DE_Results_Full.csv")
        write.csv(full_results, results_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, results_file)
        message("[DE-LIMP] Claude export: DE results OK (", nrow(full_results), " rows)")

        # --- 2. QC stats CSV ---
        if (!is.null(values$qc_stats) && is.data.frame(values$qc_stats) && !is.null(values$metadata)) {
          incProgress(0.4, detail = "Exporting QC metrics...")
          qc_df <- left_join(values$qc_stats, values$metadata,
            by = c("Run" = "File.Name")) %>%
            arrange(Group, Run)
          qc_file <- file.path(tmp_dir, "QC_Metrics.csv")
          write.csv(qc_df, qc_file, row.names = FALSE)
          files_to_zip <- c(files_to_zip, qc_file)
          message("[DE-LIMP] Claude export: QC stats OK")
        } else {
          message("[DE-LIMP] Claude export: QC stats skipped (NULL or not data.frame)")
        }

        # --- 3. Expression matrix CSV ---
        incProgress(0.5, detail = "Exporting expression matrix...")
        expr_mat <- values$y_protein$E
        expr_df <- as.data.frame(expr_mat) %>% rownames_to_column("Protein.Group")
        # Add Detection_Class column (DPC-Quant transparency)
        n_obs_export <- values$y_protein$other$n.observations
        expr_df$Detection_Class <- compute_detection_class(n_obs_export, rownames(expr_mat))
        # Reorder: Protein.Group, Detection_Class, then samples
        if ("Detection_Class" %in% colnames(expr_df)) {
          id_cols_ai <- c("Protein.Group", "Detection_Class")
          expr_df <- expr_df[, c(id_cols_ai, setdiff(colnames(expr_df), id_cols_ai))]
        }
        expr_file <- file.path(tmp_dir, "Expression_Matrix.csv")
        write.csv(expr_df, expr_file, row.names = FALSE)
        files_to_zip <- c(files_to_zip, expr_file)
        message("[DE-LIMP] Claude export: expression matrix OK")

        # --- 4. Phospho results CSV (if available) ---
        phospho_note <- ""
        if (!is.null(values$phospho_fit)) {
          incProgress(0.6, detail = "Exporting phospho data...")
          phospho_note <- "\n\n**Phosphoproteomics data included** — see `Phospho_DE_Results.csv` for site-level results."
          tryCatch({
            phospho_contrasts <- colnames(values$phospho_fit$contrasts)
            phospho_results <- do.call(rbind, lapply(phospho_contrasts, function(cname) {
              tt <- limma::topTable(values$phospho_fit, coef = cname, number = Inf) %>% as.data.frame()
              tt$SiteID <- rownames(tt)
              tt$Contrast <- cname
              tt
            }))
            phospho_file <- file.path(tmp_dir, "Phospho_DE_Results.csv")
            write.csv(phospho_results, phospho_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, phospho_file)
          }, error = function(e) NULL)
        }

        # --- 5. Session RDS (full app state — reload into DE-LIMP) ---
        incProgress(0.55, detail = "Saving session state...")
        rds_note <- ""
        tryCatch({
          session_data <- list(
            raw_data = values$raw_data, metadata = values$metadata,
            fit = values$fit, y_protein = values$y_protein,
            dpc_fit = values$dpc_fit, design = values$design,
            qc_stats = values$qc_stats,
            gsea_results = values$gsea_results,
            gsea_results_cache = values$gsea_results_cache,
            gsea_last_contrast = values$gsea_last_contrast,
            gsea_last_org_db = values$gsea_last_org_db,
            repro_log = values$repro_log,
            phospho_detected = values$phospho_detected,
            phospho_site_matrix = values$phospho_site_matrix,
            phospho_fit = values$phospho_fit,
            ksea_results = values$ksea_results,
            mofa_object = values$mofa_object,
            mofa_variance_explained = values$mofa_variance_explained,
            mofa_last_run_params = values$mofa_last_run_params,
            diann_search_settings = values$diann_search_settings,
            saved_at = Sys.time(),
            app_version = paste0("DE-LIMP v", values$app_version)
          )
          rds_file <- file.path(tmp_dir, "Session.rds")
          saveRDS(session_data, rds_file)
          files_to_zip <- c(files_to_zip, rds_file)
          rds_size <- round(file.size(rds_file) / 1024 / 1024, 1)
          rds_note <- paste0("\n- **`Session.rds`** — Full DE-LIMP session state (", rds_size,
                             " MB). Reload via DE-LIMP > Load Session to restore all results\n")
        }, error = function(e) message("[DE-LIMP] Could not save RDS: ", e$message))

        # --- 5b. Group assignments CSV ---
        groups_note <- ""
        if (!is.null(values$metadata)) {
          groups_df <- values$metadata %>%
            dplyr::select(File.Name, Group) %>%
            filter(Group != "")
          if (nrow(groups_df) > 0) {
            groups_file <- file.path(tmp_dir, "Group_Assignments.csv")
            write.csv(groups_df, groups_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, groups_file)
            groups_note <- "\n- **`Group_Assignments.csv`** — Sample-to-group mapping used in this analysis\n"
          }
        }

        # --- 5c. Pipeline parameters summary ---
        params_note <- ""
        tryCatch({
          params <- c(
            "DE-LIMP Analysis Parameters",
            paste0("Export date: ", format(Sys.time(), "%Y-%m-%d %H:%M")),
            paste0("App version: DE-LIMP v", values$app_version),
            paste0("R version: ", R.version.string),
            ""
          )
          # Contrasts
          params <- c(params, "CONTRASTS:", paste0("  ", all_contrasts), "")
          # Groups
          if (!is.null(values$metadata)) {
            grp_counts <- table(values$metadata$Group[values$metadata$Group != ""])
            params <- c(params, "GROUPS:",
              paste0("  ", names(grp_counts), ": n=", grp_counts), "")
          }
          # Covariates
          if (!is.null(values$cov1_name) && nzchar(values$cov1_name))
            params <- c(params, paste0("Covariate 1: ", values$cov1_name))
          if (!is.null(values$cov2_name) && nzchar(values$cov2_name))
            params <- c(params, paste0("Covariate 2: ", values$cov2_name))
          # DIA-NN search settings
          ss <- values$diann_search_settings
          if (!is.null(ss) && is.list(ss)) {
            sp <- ss$search_params
            params <- c(params, "", "DIA-NN SEARCH SETTINGS:",
              if (!is.null(ss$diann_version) && nzchar(ss$diann_version))
                paste0("  DIA-NN version: ", ss$diann_version) else NULL,
              paste0("  Search mode: ", ss$search_mode %||% "unknown"),
              paste0("  FASTA: ", paste(basename(ss$fasta_files), collapse = ", ")),
              if (!is.null(ss$fasta_seq_count) && !is.na(ss$fasta_seq_count))
                paste0("  FASTA sequences: ", format(ss$fasta_seq_count, big.mark = ","))
              else NULL,
              paste0("  Enzyme: ", sp$enzyme),
              paste0("  Missed cleavages: ", sp$missed_cleavages),
              paste0("  FDR: ", sp$qvalue),
              paste0("  MBR: ", if (isTRUE(sp$mbr)) "enabled" else "disabled"),
              if (!is.null(sp$mass_acc_mode) && nzchar(sp$mass_acc_mode))
                paste0("  Mass accuracy mode: ", sp$mass_acc_mode) else NULL,
              if (!is.null(sp$mass_acc) && !is.na(sp$mass_acc))
                paste0("  Mass accuracy (MS2): ", sp$mass_acc, " ppm") else NULL,
              if (!is.null(sp$mass_acc_ms1) && !is.na(sp$mass_acc_ms1))
                paste0("  Mass accuracy (MS1): ", sp$mass_acc_ms1, " ppm") else NULL,
              if (!is.null(sp$scan_window) && !is.na(sp$scan_window) && sp$scan_window > 0)
                paste0("  Scan window: ", sp$scan_window) else NULL,
              if (isTRUE(sp$mod_met_ox))
                "  Variable mod: Methionine oxidation (UniMod:35)" else NULL,
              if (isTRUE(sp$mod_nterm_acetyl))
                "  Variable mod: N-terminal acetylation (UniMod:1)" else NULL,
              if (!is.null(sp$extra_var_mods) && nzchar(sp$extra_var_mods))
                paste0("  Extra variable mods: ", sp$extra_var_mods) else NULL,
              if (!is.null(sp$min_fr_mz) && !is.na(sp$min_fr_mz))
                paste0("  Fragment m/z range: ", sp$min_fr_mz, "-", sp$max_fr_mz %||% 1800) else NULL,
              if (!is.null(sp$min_pr_charge) && !is.na(sp$min_pr_charge))
                paste0("  Precursor charge range: ", sp$min_pr_charge, "-", sp$max_pr_charge %||% 4) else NULL,
              if (isTRUE(sp$rt_profiling))
                "  RT profiling: enabled" else NULL,
              if (!is.null(ss$normalization) && nzchar(ss$normalization))
                paste0("  DIA-NN normalization: ", ss$normalization) else NULL,
              if (!is.null(sp$extra_cli_flags) && nzchar(sp$extra_cli_flags))
                paste0("  Extra CLI flags: ", sp$extra_cli_flags) else NULL,
              if (!is.null(ss$n_raw_files) && !is.na(ss$n_raw_files))
                paste0("  Raw files searched: ", ss$n_raw_files) else NULL,
              if (isTRUE(ss$imported_from_log))
                "  (Settings imported from DIA-NN log file)" else NULL)
          }
          # Package versions
          params <- c(params, "", "PACKAGE VERSIONS:",
            paste0("  limpa: ", tryCatch(as.character(packageVersion("limpa")), error = function(e) "unknown")),
            paste0("  limma: ", tryCatch(as.character(packageVersion("limma")), error = function(e) "unknown")),
            paste0("  R: ", R.version.string),
            paste0("  DE-LIMP: v", values$app_version %||% "unknown"))
          # Input file info
          n_samples <- tryCatch(ncol(values$y_protein$E), error = function(e) 0)
          n_proteins <- tryCatch(nrow(values$y_protein$E), error = function(e) 0)
          params <- c(params, "", "INPUT DATA:",
            paste0("  Total samples: ", n_samples),
            paste0("  Total proteins: ", n_proteins),
            paste0("  Source: ", if (!is.null(values$original_report_name)) values$original_report_name else "unknown"))
          params_file <- file.path(tmp_dir, "Analysis_Parameters.txt")
          writeLines(params, params_file)
          files_to_zip <- c(files_to_zip, params_file)
          params_note <- "\n- **`Analysis_Parameters.txt`** — Pipeline settings, contrasts, group sizes, DIA-NN search parameters\n"
        }, error = function(e) message("[DE-LIMP] Claude export: params section error: ", e$message))

        # --- 6. GSEA results CSV (if any ontologies have been run) ---
        gsea_note <- ""
        if (!is.null(values$gsea_results_cache) && length(values$gsea_results_cache) > 0) {
          incProgress(0.65, detail = "Exporting GSEA results...")
          tryCatch({
            gsea_all <- do.call(rbind, lapply(names(values$gsea_results_cache), function(ont) {
              res <- values$gsea_results_cache[[ont]]
              if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
                df <- as.data.frame(res)
                df$Ontology <- ont
                df
              }
            }))
            if (!is.null(gsea_all) && nrow(gsea_all) > 0) {
              gsea_file <- file.path(tmp_dir, "GSEA_Results.csv")
              write.csv(gsea_all, gsea_file, row.names = FALSE)
              files_to_zip <- c(files_to_zip, gsea_file)
              n_terms <- nrow(gsea_all)
              ontologies <- paste(names(values$gsea_results_cache), collapse = ", ")
              gsea_note <- paste0(
                "\n- **`GSEA_Results.csv`** — Gene Set Enrichment Analysis results (",
                n_terms, " terms across ", ontologies,
                "). Columns: ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, Ontology\n"
              )
            }
          }, error = function(e) NULL)
        }

        # --- 6. Methodology text file ---
        methods_note <- ""
        if (!is.null(values$methodology_text) && nzchar(values$methodology_text)) {
          methods_file <- file.path(tmp_dir, "Methods_and_References.txt")
          writeLines(values$methodology_text, methods_file)
          files_to_zip <- c(files_to_zip, methods_file)
          methods_note <- "\n- **`Methods_and_References.txt`** — Full statistical methodology, software versions, and literature references\n"
        }

        # --- 6b. Instrument metadata CSV ---
        instrument_note <- ""
        if (!is.null(values$instrument_metadata)) {
          tryCatch({
            meta <- values$instrument_metadata
            meta_clean <- meta[!sapply(meta, is.null)]
            inst_df <- data.frame(
              Parameter = names(meta_clean),
              Value = vapply(meta_clean, function(x) as.character(x), character(1)),
              stringsAsFactors = FALSE
            )
            inst_file <- file.path(tmp_dir, "Instrument_Metadata.csv")
            write.csv(inst_df, inst_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, inst_file)
            instrument_note <- "\n- **`Instrument_Metadata.csv`** — Instrument model, m/z range, gradient length, and acquisition parameters extracted from raw files\n"
          }, error = function(e) NULL)
        }

        # --- 6c. TIC chromatography QC CSV ---
        tic_note <- ""
        if (!is.null(values$tic_metrics) && nrow(values$tic_metrics) > 0) {
          tryCatch({
            tic_file <- file.path(tmp_dir, "TIC_QC_Metrics.csv")
            write.csv(values$tic_metrics, tic_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, tic_file)
            n_pass <- sum(values$tic_metrics$status == "pass", na.rm = TRUE)
            n_warn <- sum(values$tic_metrics$status == "warn", na.rm = TRUE)
            n_fail <- sum(values$tic_metrics$status == "fail", na.rm = TRUE)
            tic_note <- sprintf(
              "\n- **`TIC_QC_Metrics.csv`** — Per-run chromatography QC: AUC, peak RT, gradient width, baseline ratio, late signal, shape correlation (%d pass, %d warn, %d fail)\n",
              n_pass, n_warn, n_fail)
          }, error = function(e) NULL)
        }

        # --- 6d. Excluded files CSV ---
        excluded_note <- ""
        if (!is.null(values$excluded_files) && nrow(values$excluded_files) > 0) {
          tryCatch({
            excl_file <- file.path(tmp_dir, "Excluded_Files.csv")
            write.csv(values$excluded_files, excl_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, excl_file)
            excluded_note <- sprintf(
              "\n- **`Excluded_Files.csv`** \u2014 %d file(s) excluded from analysis with reasons, timestamps, group assignments, and user notes\n",
              nrow(values$excluded_files))
          }, error = function(e) NULL)
        }

        # --- 7. Reproducibility R code log ---
        repro_note <- ""
        if (!is.null(values$repro_log) && length(values$repro_log) > 0) {
          repro_file <- file.path(tmp_dir, "Reproducibility_Code.R")
          log_content <- paste(values$repro_log, collapse = "\n")
          session_info_text <- paste(capture.output(sessionInfo()), collapse = "\n")
          writeLines(paste(log_content, "\n\n# --- Session Info ---\n", session_info_text), repro_file)
          files_to_zip <- c(files_to_zip, repro_file)
          repro_note <- "\n- **`Reproducibility_Code.R`** — R code log recording every analysis step with timestamps (can reproduce the full analysis)\n"
        }

        # --- 7b. search_info.md and DIA-NN pg_matrix.tsv ---
        tryCatch({
          ss <- values$diann_search_settings
          if (!is.null(ss) && !is.null(ss$output_dir)) {
            od <- translate_storage_path(ss$output_dir, to = "hpc")
            cfg <- if (isTRUE(values$ssh_connected))
              list(host = isolate(input$ssh_host), user = isolate(input$ssh_user),
                   port = isolate(input$ssh_port) %||% 22L, key_path = isolate(input$ssh_key_path))
            else NULL

            for (fname in c("search_info.md", "report.pg_matrix.tsv")) {
              local_path <- file.path(ss$output_dir, fname)
              dest <- file.path(tmp_dir, fname)
              if (file.exists(local_path)) {
                file.copy(local_path, dest)
                files_to_zip <- c(files_to_zip, dest)
              } else if (!is.null(cfg)) {
                dl <- scp_download(cfg, file.path(od, fname), dest)
                if (dl$status == 0 && file.exists(dest))
                  files_to_zip <- c(files_to_zip, dest)
              }
            }
          }
        }, error = function(e) message("[Export] search_info/pg_matrix: ", e$message))

        # --- 7b2. Protein confidence (DPC-Quant n.observations + standard.error) ---
        tryCatch({
          n_obs <- values$y_protein$other$n.observations
          se_mat <- values$y_protein$other$standard.error
          if (!is.null(n_obs) && !is.null(se_mat)) {
            conf_df <- data.frame(Protein.Group = rownames(n_obs), stringsAsFactors = FALSE)
            for (j in seq_len(ncol(n_obs))) {
              conf_df[[paste0("nObs_", colnames(n_obs)[j])]] <- n_obs[, j]
            }
            for (j in seq_len(ncol(se_mat))) {
              conf_df[[paste0("SE_", colnames(se_mat)[j])]] <- round(se_mat[, j], 4)
            }
            conf_file <- file.path(tmp_dir, "protein_confidence.csv")
            write.csv(conf_df, conf_file, row.names = FALSE)
            files_to_zip <- c(files_to_zip, conf_file)
          }
        }, error = function(e) message("[Export] protein_confidence: ", e$message))

        # --- 7c. Data quality summary (per-sample protein counts + missingness) ---
        tryCatch({
          pg_file_path <- file.path(tmp_dir, "report.pg_matrix.tsv")
          if (file.exists(pg_file_path)) {
            pg <- read.delim(pg_file_path, stringsAsFactors = FALSE, check.names = FALSE)
            annot_cols <- c("Protein.Group", "Protein.Names", "Genes",
                            "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences")
            int_cols <- setdiff(colnames(pg), annot_cols)
            if (length(int_cols) > 0) {
              pg_mat <- as.matrix(pg[, int_cols])
              detected <- colSums(pg_mat > 0, na.rm = TRUE)
              total_pg <- nrow(pg_mat)
              quality_df <- data.frame(
                Sample = int_cols,
                Proteins_Detected = detected,
                Total_Protein_Groups = total_pg,
                Pct_Detected = round(100 * detected / total_pg, 1),
                Missing = total_pg - detected,
                Pct_Missing = round(100 * (total_pg - detected) / total_pg, 1),
                Contaminant_Proteins = sum(grepl("^Cont_", pg$Protein.Group)),
                stringsAsFactors = FALSE
              )
              if (!is.null(values$metadata)) {
                quality_df$Group <- values$metadata$Group[match(quality_df$Sample, values$metadata$File.Name)]
              }
              quality_file <- file.path(tmp_dir, "data_quality_summary.csv")
              write.csv(quality_df, quality_file, row.names = FALSE)
              files_to_zip <- c(files_to_zip, quality_file)
            }
          }
        }, error = function(e) message("[Export] data quality summary: ", e$message))

        # --- 7d. Detection matrix (per-protein precursor counts) ---
        tryCatch({
          if (!is.null(values$raw_data) && !is.null(values$raw_data$E)) {
            raw_mat <- values$raw_data$E
            raw_genes <- values$raw_data$genes
            if (!is.null(raw_genes) && "Protein.Group" %in% colnames(raw_genes)) {
              pg <- raw_genes$Protein.Group
              det_counts <- do.call(rbind, lapply(unique(pg), function(p) {
                rows <- which(pg == p)
                sub_mat <- raw_mat[rows, , drop = FALSE]
                detected <- colSums(!is.na(sub_mat) & is.finite(sub_mat))
                det_list <- as.list(detected)
                names(det_list) <- paste0("Detected_", colnames(raw_mat))
                data.frame(Protein.Group = p, Total_Precursors = nrow(sub_mat),
                  det_list, stringsAsFactors = FALSE, check.names = FALSE)
              }))
              det_df <- as.data.frame(det_counts, stringsAsFactors = FALSE)
              det_file <- file.path(tmp_dir, "detection_matrix.csv")
              write.csv(det_df, det_file, row.names = FALSE)
              files_to_zip <- c(files_to_zip, det_file)
            }
          }
        }, error = function(e) message("[Export] detection matrix: ", e$message))

        # --- 7e. SVG Figures ---
        figures_note <- ""
        tryCatch({
          fig_dir <- file.path(tmp_dir, "figures")
          dir.create(fig_dir, showWarnings = FALSE)
          fig_files <- character(0)

          # Volcano plot SVG (first contrast)
          tryCatch({
            df_volc <- topTable(values$fit, coef = all_contrasts[1], number = Inf) %>% as.data.frame()
            if (!"Protein.Group" %in% colnames(df_volc)) df_volc <- df_volc %>% rownames_to_column("Protein.Group")
            df_volc$Significance <- ifelse(df_volc$adj.P.Val < 0.05, "Significant", "Not Sig")

            p_volcano <- ggplot2::ggplot(df_volc, ggplot2::aes(x = logFC, y = -log10(P.Value), color = Significance)) +
              ggplot2::geom_point(alpha = 0.5, size = 1) +
              ggplot2::scale_color_manual(values = c("Significant" = "#E74C3C", "Not Sig" = "#BDC3C7")) +
              ggplot2::geom_vline(xintercept = c(-0.6, 0.6), linetype = "dashed", color = "orange") +
              ggplot2::labs(title = paste("Volcano:", all_contrasts[1]), x = "logFC", y = "-log10(P-Value)") +
              ggplot2::theme_bw()
            svg(file.path(fig_dir, "volcano.svg"), width = 8, height = 6)
            print(p_volcano)
            dev.off()
            fig_files <- c(fig_files, file.path(fig_dir, "volcano.svg"))
            message("[DE-LIMP] Claude export: volcano SVG OK")
          }, error = function(e) { try(dev.off(), silent = TRUE); message("[Export] Volcano SVG: ", e$message) })

          # Heatmap SVG (top 20 DE proteins)
          tryCatch({
            df_volc_ht <- topTable(values$fit, coef = all_contrasts[1], number = Inf) %>% as.data.frame()
            top20 <- head(df_volc_ht[order(df_volc_ht$adj.P.Val), ], 20)
            mat <- values$y_protein$E[rownames(top20), , drop = FALSE]
            mat_z <- t(apply(mat, 1, function(x) (x - mean(x)) / sd(x)))

            svg(file.path(fig_dir, "heatmap_top20.svg"), width = 10, height = 8)
            meta_ht <- values$metadata[match(colnames(mat), values$metadata$File.Name), ]
            ha <- ComplexHeatmap::HeatmapAnnotation(Group = factor(meta_ht$Group))
            ComplexHeatmap::draw(ComplexHeatmap::Heatmap(mat_z, name = "Z-score", top_annotation = ha,
              cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = FALSE))
            dev.off()
            fig_files <- c(fig_files, file.path(fig_dir, "heatmap_top20.svg"))
            message("[DE-LIMP] Claude export: heatmap SVG OK")
          }, error = function(e) { try(dev.off(), silent = TRUE); message("[Export] Heatmap SVG: ", e$message) })

          # Violin plots SVG (top 10 up and top 10 down)
          tryCatch({
            df_volc_vln <- topTable(values$fit, coef = all_contrasts[1], number = Inf) %>% as.data.frame()
            if (!"Protein.Group" %in% colnames(df_volc_vln)) df_volc_vln <- df_volc_vln %>% rownames_to_column("Protein.Group")
            sig_prots <- df_volc_vln[df_volc_vln$adj.P.Val < 0.05, ]

            for (direction in c("up", "down")) {
              prot_ids <- if (direction == "up") {
                head(sig_prots[order(-sig_prots$logFC), "Protein.Group"], 10)
              } else {
                head(sig_prots[order(sig_prots$logFC), "Protein.Group"], 10)
              }
              if (length(prot_ids) == 0) next
              mat_vln <- values$y_protein$E[prot_ids, , drop = FALSE]
              long <- as.data.frame(mat_vln) %>%
                tibble::rownames_to_column("Protein") %>%
                tidyr::pivot_longer(-Protein, names_to = "File.Name", values_to = "LogIntensity")
              long <- dplyr::left_join(long, values$metadata, by = "File.Name")

              p <- ggplot(long, aes(x = Group, y = LogIntensity, fill = Group)) +
                geom_violin(alpha = 0.5, trim = FALSE) +
                geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
                facet_wrap(~Protein, scales = "free_y", ncol = 5) +
                theme_bw() +
                labs(title = paste("Top 10", if (direction == "up") "Upregulated" else "Downregulated", "Proteins"),
                     y = "Log2 Intensity") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))

              svg(file.path(fig_dir, paste0("violin_top10_", direction, ".svg")), width = 14, height = 6)
              print(p)
              dev.off()
              fig_files <- c(fig_files, file.path(fig_dir, paste0("violin_top10_", direction, ".svg")))
            }
            message("[DE-LIMP] Claude export: violin SVGs OK")
          }, error = function(e) message("[Export] Violin SVGs: ", e$message))

          # PCA plot SVG
          tryCatch({
            pca_res <- prcomp(t(values$y_protein$E), scale. = TRUE)
            pca_df <- data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                                  File.Name = rownames(pca_res$x))
            pca_df <- dplyr::left_join(pca_df, values$metadata, by = "File.Name")
            var_exp <- round(100 * summary(pca_res)$importance[2, 1:2], 1)

            p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
              geom_point(size = 3) +
              stat_ellipse(level = 0.95, linetype = "dashed") +
              labs(title = "PCA", x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
              theme_bw()
            svg(file.path(fig_dir, "pca.svg"), width = 8, height = 6)
            print(p_pca)
            dev.off()
            fig_files <- c(fig_files, file.path(fig_dir, "pca.svg"))
            message("[DE-LIMP] Claude export: PCA SVG OK")
          }, error = function(e) { try(dev.off(), silent = TRUE); message("[Export] PCA SVG: ", e$message) })

          if (length(fig_files) > 0) {
            files_to_zip <- c(files_to_zip, fig_files)
            figures_note <- paste0(
              "\n- **`figures/`** — Publication-quality SVG figures: volcano plot, heatmap (top 20 DE), PCA, violin plots (top 10 up/down)\n")
            message("[DE-LIMP] Claude export: ", length(fig_files), " SVG figures generated")
          }
        }, error = function(e) message("[Export] SVG figures section: ", e$message))

        # --- 8. Build the prompt .md ---
        incProgress(0.7, detail = "Assembling prompt...")

        message("[DE-LIMP] Claude export: sections 1-7 OK, building prompt...")
        # QC summary for inline prompt
        qc_inline <- ""
        if (!is.null(values$qc_stats) && is.data.frame(values$qc_stats) && !is.null(values$metadata)) {
          qc_summary <- left_join(values$qc_stats, values$metadata,
            by = c("Run" = "File.Name")) %>%
            dplyr::select(Run, Group, Precursors, Proteins, MS1_Signal) %>%
            arrange(Group, Run)
          qc_text <- paste(capture.output(print(as.data.frame(qc_summary), row.names = FALSE)), collapse = "\n")
          qc_inline <- paste0(
            "\n\n--- QC METRICS SUMMARY ---\n",
            "Full data in `QC_Metrics.csv`. Key columns:\n",
            "- Precursors: identified precursor ions per run\n",
            "- Proteins: identified protein groups per run\n",
            "- MS1_Signal: median MS1 intensity (total signal proxy)\n\n",
            qc_text
          )
        }

        # Experimental design (with covariates if present)
        design_section <- ""
        if (!is.null(values$metadata)) {
          group_summary <- values$metadata %>%
            filter(Group != "") %>%
            group_by(Group) %>%
            summarise(N_Replicates = n(), .groups = "drop")
          design_text <- paste(capture.output(print(as.data.frame(group_summary), row.names = FALSE)), collapse = "\n")
          design_section <- paste0("\n\n--- EXPERIMENTAL DESIGN ---\n", design_text)

          # Add covariate info if present
          cov1 <- values$cov1_name
          cov2 <- values$cov2_name
          has_cov1 <- !is.null(cov1) && nzchar(cov1) && cov1 != "Covariate1" &&
            any(nzchar(values$metadata$Covariate1))
          has_cov2 <- !is.null(cov2) && nzchar(cov2) && cov2 != "Covariate2" &&
            any(nzchar(values$metadata$Covariate2))
          if (has_cov1 || has_cov2) {
            cov_lines <- "\n\nCOVARIATES IN MODEL:"
            if (has_cov1) {
              cov1_vals <- values$metadata %>%
                filter(Group != "", nzchar(Covariate1)) %>%
                dplyr::select(File.Name, Group, Covariate1)
              colnames(cov1_vals)[3] <- cov1
              cov1_text <- paste(capture.output(print(as.data.frame(cov1_vals), row.names = FALSE)), collapse = "\n")
              cov_lines <- paste0(cov_lines, "\n\n", cov1, ":\n", cov1_text)
            }
            if (has_cov2) {
              cov2_vals <- values$metadata %>%
                filter(Group != "", nzchar(Covariate2)) %>%
                dplyr::select(File.Name, Group, Covariate2)
              colnames(cov2_vals)[3] <- cov2
              cov2_text <- paste(capture.output(print(as.data.frame(cov2_vals), row.names = FALSE)), collapse = "\n")
              cov_lines <- paste0(cov_lines, "\n\n", cov2, ":\n", cov2_text)
            }
            design_section <- paste0(design_section, cov_lines)
          }
        }

        # DIA-NN search settings inline block for prompt
        search_settings_inline <- ""
        ss <- values$diann_search_settings
        if (!is.null(ss) && is.list(ss)) {
          sp <- ss$search_params
          ss_lines <- c(
            "--- DIA-NN SEARCH SETTINGS ---",
            if (!is.null(ss$diann_version) && nzchar(ss$diann_version))
              paste0("DIA-NN version: ", ss$diann_version) else NULL,
            paste0("Search mode: ", ss$search_mode %||% "unknown"),
            paste0("FASTA database: ", paste(basename(ss$fasta_files), collapse = ", ")),
            if (!is.null(ss$fasta_seq_count) && !is.na(ss$fasta_seq_count))
              paste0("FASTA sequences: ", format(ss$fasta_seq_count, big.mark = ",")) else NULL,
            paste0("Enzyme: ", sp$enzyme),
            paste0("Missed cleavages: ", sp$missed_cleavages),
            paste0("Precursor FDR: ", sp$qvalue),
            paste0("Match between runs (MBR): ", if (isTRUE(sp$mbr)) "enabled" else "disabled"),
            if (!is.null(sp$mass_acc_mode) && nzchar(sp$mass_acc_mode))
              paste0("Mass accuracy mode: ", sp$mass_acc_mode) else NULL,
            if (!is.null(sp$mass_acc) && !is.na(sp$mass_acc))
              paste0("Mass accuracy MS2: ", sp$mass_acc, " ppm") else NULL,
            if (!is.null(sp$mass_acc_ms1) && !is.na(sp$mass_acc_ms1))
              paste0("Mass accuracy MS1: ", sp$mass_acc_ms1, " ppm") else NULL,
            if (!is.null(sp$scan_window) && !is.na(sp$scan_window) && sp$scan_window > 0)
              paste0("Scan window: ", sp$scan_window) else NULL,
            if (isTRUE(sp$mod_met_ox)) "Variable mod: Methionine oxidation" else NULL,
            if (isTRUE(sp$mod_nterm_acetyl)) "Variable mod: N-terminal acetylation" else NULL,
            if (!is.null(sp$extra_var_mods) && nzchar(sp$extra_var_mods))
              paste0("Extra variable mods: ", sp$extra_var_mods) else NULL,
            if (!is.null(sp$min_fr_mz) && !is.na(sp$min_fr_mz))
              paste0("Fragment m/z range: ", sp$min_fr_mz, "-", sp$max_fr_mz %||% 1800) else NULL,
            if (!is.null(sp$min_pr_charge) && !is.na(sp$min_pr_charge))
              paste0("Precursor charge range: ", sp$min_pr_charge, "-", sp$max_pr_charge %||% 4) else NULL,
            if (isTRUE(sp$rt_profiling)) "RT profiling: enabled" else NULL,
            if (!is.null(ss$normalization) && nzchar(ss$normalization))
              paste0("DIA-NN normalization: ", ss$normalization) else NULL,
            if (!is.null(sp$extra_cli_flags) && nzchar(sp$extra_cli_flags))
              paste0("Extra CLI flags: ", sp$extra_cli_flags) else NULL,
            if (!is.null(ss$n_raw_files) && !is.na(ss$n_raw_files))
              paste0("Raw files searched: ", ss$n_raw_files) else NULL
          )
          search_settings_inline <- paste0("\n\n", paste(ss_lines, collapse = "\n"))
        }

        # Instrument metadata inline block for prompt
        instrument_inline <- ""
        if (!is.null(values$instrument_metadata)) {
          inst_text <- format_instrument_for_prompt(values$instrument_metadata)
          if (nzchar(inst_text)) {
            instrument_inline <- paste0("\n\n--- INSTRUMENT & ACQUISITION ---\n", inst_text)
          }
        }

        # TIC chromatography QC inline summary
        tic_inline <- ""
        if (!is.null(values$tic_metrics) && nrow(values$tic_metrics) > 0) {
          tryCatch({
            tm <- values$tic_metrics
            n_pass <- sum(tm$status == "pass", na.rm = TRUE)
            n_warn <- sum(tm$status == "warn", na.rm = TRUE)
            n_fail <- sum(tm$status == "fail", na.rm = TRUE)
            tic_lines <- c(
              "--- CHROMATOGRAPHY QC (TIC) ---",
              sprintf("Runs analyzed: %d (%d pass, %d warn, %d fail)", nrow(tm), n_pass, n_warn, n_fail))
            if (any(tm$valid, na.rm = TRUE)) {
              valid <- tm[tm$valid, ]
              tic_lines <- c(tic_lines,
                sprintf("Median AUC: %.1fM (range: %.1f-%.1fM)",
                  median(valid$total_auc, na.rm = TRUE) / 1e6,
                  min(valid$total_auc, na.rm = TRUE) / 1e6,
                  max(valid$total_auc, na.rm = TRUE) / 1e6),
                sprintf("Median gradient width: %.1f min", median(valid$gradient_width_min, na.rm = TRUE)),
                sprintf("Shape correlation range: %.3f-%.3f",
                  min(valid$shape_r, na.rm = TRUE), max(valid$shape_r, na.rm = TRUE)))
            }
            # List flagged runs
            flagged <- tm[nzchar(tm$flags) & tm$flags != "", ]
            if (nrow(flagged) > 0) {
              tic_lines <- c(tic_lines, "", "Flagged runs:")
              for (r in seq_len(nrow(flagged))) {
                tic_lines <- c(tic_lines, sprintf("  %s [%s]: %s",
                  sub("\\.d$", "", flagged$run[r]), flagged$status[r], flagged$flags[r]))
              }
            }
            tic_inline <- paste0("\n\n", paste(tic_lines, collapse = "\n"))
          }, error = function(e) NULL)
        }

        # Excluded files inline summary
        excluded_inline <- ""
        if (!is.null(values$excluded_files) && nrow(values$excluded_files) > 0) {
          ef <- values$excluded_files
          exc_lines <- c(
            "--- EXCLUDED FILES ---",
            sprintf("%d file(s) were excluded before analysis:", nrow(ef)))
          for (i in seq_len(nrow(ef))) {
            note_part <- if (nzchar(ef$user_note[i])) paste0(" | Note: ", ef$user_note[i]) else ""
            grp_part <- if (nzchar(ef$group[i]) && ef$group[i] != "[Excluded]")
              paste0(" | Group: ", ef$group[i]) else ""
            exc_lines <- c(exc_lines, sprintf("  %s: %s (source: %s)%s%s",
              ef$filename[i], ef$reason[i], ef$source[i], grp_part, note_part))
          }
          # Check for group bias
          if (any(nzchar(ef$group) & ef$group != "[Excluded]")) {
            grp_counts <- table(ef$group[nzchar(ef$group) & ef$group != "[Excluded]"])
            exc_lines <- c(exc_lines, "",
              "Excluded files by group: ",
              paste(paste0("  ", names(grp_counts), ": ", grp_counts), collapse = "\n"),
              "IMPORTANT: Assess whether exclusions are biased toward specific experimental groups.")
          }
          excluded_inline <- paste0("\n\n", paste(exc_lines, collapse = "\n"))
        }

        # Missingness summary (pre vs post DPC-Quant protein quantification)
        missingness_inline <- ""
        tryCatch({
          if (!is.null(values$raw_data) && !is.null(values$y_protein) && !is.null(values$metadata)) {
            raw_mat <- values$raw_data$E
            post_mat <- values$y_protein$E
            meta <- values$metadata %>% filter(Group != "")

            # Per-group missingness in raw data (precursor-level)
            miss_lines <- c("--- DATA COMPLETENESS ---",
              paste0("Precursor-level (raw): ", nrow(raw_mat), " precursors x ", ncol(raw_mat), " samples"),
              paste0("Protein-level (after pipeline): ", nrow(post_mat), " proteins x ", ncol(post_mat), " samples"),
              paste0("Overall raw missingness: ", round(100 * sum(is.na(raw_mat)) / length(raw_mat), 1), "%"),
              paste0("Overall post-pipeline missingness: ", round(100 * sum(is.na(post_mat)) / length(post_mat), 1), "%"),
              "", "Per-group raw missingness (precursor-level):")
            for (grp in unique(meta$Group)) {
              grp_samples <- meta$File.Name[meta$Group == grp]
              grp_cols <- colnames(raw_mat) %in% grp_samples
              if (any(grp_cols)) {
                grp_mat <- raw_mat[, grp_cols, drop = FALSE]
                pct <- round(100 * sum(is.na(grp_mat)) / length(grp_mat), 1)
                miss_lines <- c(miss_lines, paste0("  ", grp, ": ", pct, "% missing (n=", sum(grp_cols), " samples)"))
              }
            }
            missingness_inline <- paste0("\n\n", paste(miss_lines, collapse = "\n"))
          }
        }, error = function(e) NULL)

        # Dynamic range from precursor-level raw intensities
        dynamic_range_inline <- ""
        tryCatch({
          if (!is.null(values$raw_data) && !is.null(values$raw_data$E)) {
            raw_mat <- values$raw_data$E
            dr_lines <- c("--- DYNAMIC RANGE (precursor-level, log2 intensities) ---")
            all_vals <- raw_mat[!is.na(raw_mat)]
            dr_lines <- c(dr_lines,
              sprintf("Global: log2 %.1f - %.1f (%.1f orders of magnitude)",
                min(all_vals), max(all_vals), (max(all_vals) - min(all_vals)) / log2(10)),
              sprintf("Median log2 intensity: %.1f, IQR: %.1f", median(all_vals), IQR(all_vals)),
              sprintf("Precursors: %d, Samples: %d", nrow(raw_mat), ncol(raw_mat)),
              "", "Per-sample dynamic range:")
            for (j in seq_len(ncol(raw_mat))) {
              v <- raw_mat[, j]
              v <- v[!is.na(v)]
              if (length(v) > 0) {
                dr_lines <- c(dr_lines, sprintf("  %s: %.1f - %.1f (%.1f orders, %d precursors, %.0f%% missing)",
                  colnames(raw_mat)[j], min(v), max(v),
                  (max(v) - min(v)) / log2(10), length(v),
                  100 * mean(is.na(raw_mat[, j]))))
              }
            }
            dynamic_range_inline <- paste0("\n\n", paste(dr_lines, collapse = "\n"))
          }
        }, error = function(e) NULL)

        # MOFA2 variance explained summary
        mofa_inline <- ""
        tryCatch({
          if (!is.null(values$mofa_variance_explained)) {
            r2 <- values$mofa_variance_explained$r2_per_factor
            mofa_lines <- c("--- MOFA2 MULTI-OMICS INTEGRATION ---")
            for (group_name in names(r2)) {
              mat <- r2[[group_name]]
              if (length(names(r2)) > 1)
                mofa_lines <- c(mofa_lines, paste0("Group: ", group_name))
              for (i in seq_len(nrow(mat))) {
                vals <- paste(paste0(colnames(mat), ": ", round(mat[i, ], 1), "%"), collapse = ", ")
                mofa_lines <- c(mofa_lines, paste0("  ", rownames(mat)[i], " — ", vals))
              }
            }
            # Total variance per view
            r2_total <- values$mofa_variance_explained$r2_total
            if (!is.null(r2_total)) {
              for (group_name in names(r2_total)) {
                totals <- r2_total[[group_name]]
                mofa_lines <- c(mofa_lines, "",
                  paste0("Total variance explained per view: ",
                    paste(paste0(names(totals), ": ", round(totals, 1), "%"), collapse = ", ")))
              }
            }
            if (!is.null(values$mofa_last_run_params)) {
              mp <- values$mofa_last_run_params
              mofa_lines <- c(mofa_lines, "",
                paste0("MOFA2 params: ", mp$n_factors %||% "?", " factors, ",
                  length(mp$views %||% list()), " views (",
                  paste(mp$views %||% "?", collapse = ", "), ")"))
            }
            mofa_inline <- paste0("\n\n", paste(mofa_lines, collapse = "\n"))
          }
        }, error = function(e) NULL)

        # Phosphoproteomics summary
        phospho_inline <- ""
        tryCatch({
          if (!is.null(values$phospho_fit)) {
            phospho_contrasts <- colnames(values$phospho_fit$contrasts)
            phospho_lines <- c("--- PHOSPHOPROTEOMICS SUMMARY ---",
              paste0("Total phosphosites quantified: ", nrow(values$phospho_fit$coefficients)),
              paste0("Contrasts: ", paste(phospho_contrasts, collapse = ", ")))
            for (pc in phospho_contrasts) {
              tt <- limma::topTable(values$phospho_fit, coef = pc, number = Inf)
              n_sig <- sum(tt$adj.P.Val < 0.05, na.rm = TRUE)
              n_up <- sum(tt$adj.P.Val < 0.05 & tt$logFC > 0, na.rm = TRUE)
              n_down <- sum(tt$adj.P.Val < 0.05 & tt$logFC < 0, na.rm = TRUE)
              phospho_lines <- c(phospho_lines,
                paste0("  ", pc, ": ", n_sig, " significant sites (", n_up, " up, ", n_down, " down)"))
              # Top 5 sites by absolute logFC
              top5 <- head(tt[order(-abs(tt$logFC)), ], 5)
              if (nrow(top5) > 0) {
                site_ids <- if (!is.null(top5$SiteID)) top5$SiteID else rownames(top5)
                for (j in seq_len(nrow(top5))) {
                  phospho_lines <- c(phospho_lines,
                    paste0("    ", site_ids[j], ": logFC=", round(top5$logFC[j], 2),
                      ", adj.P=", formatC(top5$adj.P.Val[j], format = "e", digits = 2)))
                }
              }
            }
            phospho_inline <- paste0("\n\n", paste(phospho_lines, collapse = "\n"))
          }
        }, error = function(e) NULL)

        n_total <- nrow(topTable(values$fit, coef = all_contrasts[1], number = Inf))

        # Ensure all inline/note variables are character strings (not NULL/list)
        safe_str <- function(x) if (is.null(x) || length(x) == 0) "" else paste(as.character(x), collapse = "")
        design_section <- safe_str(design_section)
        qc_inline <- safe_str(qc_inline)
        search_settings_inline <- safe_str(search_settings_inline)
        instrument_inline <- safe_str(instrument_inline)
        tic_inline <- safe_str(tic_inline)
        excluded_inline <- safe_str(excluded_inline)
        missingness_inline <- safe_str(missingness_inline)
        dynamic_range_inline <- safe_str(dynamic_range_inline)
        mofa_inline <- safe_str(mofa_inline)
        phospho_inline <- safe_str(phospho_inline)
        gsea_note <- safe_str(gsea_note)
        phospho_note <- safe_str(phospho_note)
        instrument_note <- safe_str(instrument_note)
        tic_note <- safe_str(tic_note)
        excluded_note <- safe_str(excluded_note)
        methods_note <- safe_str(methods_note)
        repro_note <- safe_str(repro_note)
        rds_note <- safe_str(rds_note)
        groups_note <- safe_str(groups_note)
        params_note <- safe_str(params_note)
        figures_note <- safe_str(figures_note)

        prompt <- paste0(
          "# Proteomics Differential Expression Analysis\n\n",
          "You are a senior proteomics and systems biology consultant. ",
          "Analyze the following differential expression results from a DIA-NN / limma pipeline. ",
          "The data was processed by DE-LIMP (Differential Expression - LIMPA Pipeline).\n\n",
          "## Attached Data Files\n\n",
          "- **`DE_Results_Full.csv`** — Complete DE statistics for all ", n_total, " proteins across ",
          ctx$n_contrasts, " comparison(s). Columns: Protein.Group, logFC, AveExpr, t, P.Value, adj.P.Val, B, Contrast\n",
          "- **`Expression_Matrix.csv`** — Log2 expression values for all proteins across all samples\n",
          if (!is.null(values$qc_stats)) "- **`QC_Metrics.csv`** — Per-sample QC metrics (precursor/protein counts, MS1 signal)\n" else "",
          gsea_note,
          phospho_note,
          instrument_note,
          tic_note,
          excluded_note,
          methods_note,
          repro_note,
          rds_note,
          groups_note,
          params_note,
          figures_note, "\n\n",
          "## OUTPUT FORMAT\n\n",
          "Generate a publication-quality report suitable for sharing with collaborators. ",
          "The report should be professional, well-structured, and written in a scientific but accessible tone. ",
          "Use clear section headers, bullet points for key findings, and data tables where appropriate. ",
          "Include specific protein names, fold changes, and p-values to support all claims. ",
          "Write it as if preparing a formal analysis report for a PI — not a casual summary.\n\n",
          "Format the output as a structured document with numbered sections. ",
          "For tables, use markdown table syntax. For key warnings or critical findings, use bold text. ",
          "The report should stand alone — a reader should understand the experiment, methods, findings, and recommendations without additional context.\n\n",
          "Include references to the attached figures in your analysis (e.g., 'As shown in the volcano plot (figures/volcano.svg)...'). ",
          "When discussing top DE proteins, reference the violin plots showing their expression distributions across groups.\n\n",
          "Please provide a comprehensive analysis with these sections:\n\n",
          "## Overview\n",
          "Number of comparisons analyzed, total significant proteins per comparison (up/down split). ",
          "Overall assessment of the experiment's quality and scope.\n\n",
          "## QC Assessment\n",
          "Evaluate the technical quality of the experiment based on the QC metrics. ",
          "Comment on consistency of precursor/protein identifications across replicates and groups, ",
          "and flag any outlier samples or systematic biases.\n\n",
          "## Key Findings Per Comparison\n",
          "For each comparison: highlight the top upregulated and downregulated proteins by fold-change ",
          "(use gene names). Note any comparison with unusually few or many significant hits.\n\n",
          "## Cross-Comparison Biomarkers\n",
          "Proteins significant in multiple comparisons are highest-confidence candidates. ",
          "Discuss consistency of direction (always up, always down, or mixed across comparisons).\n\n",
          "## High-Confidence Biomarker Insights\n",
          "For the most stable proteins (lowest coefficient of variation): discuss their known biological functions, ",
          "pathway involvement, and disease associations. ",
          "Assess their potential as reliable biomarkers based on the combination of low CV, ",
          "significant p-value, and meaningful fold-change.\n\n",
          if (nzchar(gsea_note)) paste0(
            "## Pathway & Gene Set Enrichment Analysis\n",
            "GSEA results are in `GSEA_Results.csv`. Summarize the top enriched pathways by ontology. ",
            "Highlight pathways with the highest normalized enrichment scores (NES). ",
            "Connect enriched pathways to the DE protein findings above.\n\n"
          ) else "",
          if (nzchar(mofa_inline)) paste0(
            "## Multi-Omics Integration (MOFA2)\n",
            "MOFA2 variance explained data is provided below. Discuss:\n",
            "- Which factors capture the most variance and in which views\n",
            "- Whether the multi-omics integration reveals structure not visible in single-view analysis\n",
            "- How the MOFA factors relate to the experimental groups and DE findings\n\n"
          ) else "",
          if (nzchar(phospho_inline)) paste0(
            "## Phosphoproteomics\n",
            "Phosphosite-level DE results are provided below and in `Phospho_DE_Results.csv`. Discuss:\n",
            "- Number and direction of significant phosphosites per comparison\n",
            "- Top regulated phosphosites and their known regulatory roles\n",
            "- Whether phospho changes are concordant or discordant with protein-level changes\n",
            "- Any kinase activity implications from the regulated sites\n\n"
          ) else "",
          if (nzchar(excluded_inline)) paste0(
            "## Excluded Files Assessment\n",
            "Files were excluded from analysis (details below and in `Excluded_Files.csv`). Assess:\n",
            "- Whether exclusions are biased toward specific experimental groups (could introduce bias)\n",
            "- Whether the QC reasons (TIC diagnostics) suggest technical failures vs potential biological effects\n",
            "- Impact on statistical power given the remaining sample sizes per group\n",
            "- Any recommendations for re-inclusion or follow-up experiments\n\n"
          ) else "",
          if (nzchar(missingness_inline)) paste0(
            "## Data Completeness\n",
            "Missingness data is provided below. Comment on:\n",
            "- Overall data completeness and whether missingness varies across groups\n",
            "- Whether group-specific missingness could introduce bias\n",
            "- The effectiveness of DPC-Quant quantification (compare pre vs post-pipeline missingness — note: limpa uses probabilistic modelling, not imputation)\n\n"
          ) else "",
          "## Biological Interpretation\n",
          "Suggest what biological processes or pathways may be affected based on the protein lists. ",
          "Note any well-known protein families, complexes, or signaling cascades represented. ",
          "If the data suggests a clear biological narrative, describe it.\n\n",
          "## How This Analysis Works\n",
          "Write an educational background section that a PhD student or biologist with no mass spectrometry ",
          "or bioinformatics background can understand. Cover each stage of the pipeline in plain language, ",
          "using analogies where helpful. Include these topics:\n\n",
          "### Liquid Chromatography–Mass Spectrometry (LC-MS/MS)\n",
          "Explain what LC-MS/MS does at a high level: proteins are digested into peptides, separated by ",
          "liquid chromatography (like sorting by stickiness), then ionized and measured by mass. ",
          "Explain that the mass spectrometer measures both the mass of intact peptides (MS1) and breaks them ",
          "into fragments to identify the sequence (MS2). Keep it intuitive.\n\n",
          "### Data-Independent Acquisition (DIA)\n",
          "Explain the difference between DDA (picks the loudest signals one at a time) and DIA (systematically ",
          "scans all peptides in windows across the full mass range). Explain why DIA gives more complete, ",
          "reproducible quantification — every peptide gets measured every time, not just the most abundant ones. ",
          "Mention that DIA produces more complex data that requires specialized software to deconvolve.\n\n",
          "### DIA-NN Software\n",
          "Explain that DIA-NN is the software that takes the raw mass spectrometry data and figures out which ",
          "peptides (and therefore which proteins) are present and how abundant they are. Mention that it uses ",
          "neural networks to score peptide identifications and that it performs library-free search (predicting ",
          "what peptides should look like rather than requiring a pre-built library). Explain that it outputs a ",
          "report with protein quantities per sample.\n\n",
          "### LIMPA / limma Statistical Framework\n",
          "Explain that once we have protein quantities, we need statistics to determine which proteins are truly ",
          "different between groups vs. random noise. Explain limma's key innovation in plain terms: it borrows ",
          "information across all proteins to get better variance estimates, which is especially powerful when you ",
          "have few replicates (common in proteomics). Mention empirical Bayes moderation — the idea that a protein's ",
          "variance estimate is improved by considering how variable all the other proteins are. Explain that LIMPA ",
          "is an R package that wraps limma with proteomics-specific preprocessing (normalization, filtering).\n\n",
          "### Key Statistical Concepts\n",
          "Define these terms in plain language with brief examples from this dataset:\n",
          "- **log2 Fold Change (logFC)**: How much a protein goes up or down between groups (logFC of 1 = doubled, -1 = halved)\n",
          "- **P-value**: The probability of seeing this difference by chance alone\n",
          "- **Adjusted P-value (FDR)**: P-values corrected for testing thousands of proteins at once (Benjamini-Hochberg). ",
          "Explain the multiple testing problem with an intuitive example (e.g., flipping coins)\n",
          "- **Volcano plot**: Why it's shaped like a volcano and how to read it (x = effect size, y = significance)\n",
          "- **Coefficient of Variation (CV)**: A measure of measurement reproducibility — lower is more reliable\n",
          "- **Normalization**: Why raw intensities need correction (loading differences between samples) and how DPC-CN (Data Point Correspondence - Cyclic Normalization) works conceptually\n\n",
          "Keep the tone approachable and encouraging. Avoid jargon where possible, and define it when unavoidable.\n\n",
          if (nzchar(search_settings_inline)) paste0(
            "## DIA-NN Search Methods\n",
            "DIA-NN search parameters are provided below under 'DIA-NN SEARCH SETTINGS'. ",
            "Write a publication-ready Methods subsection describing the database search. Include:\n",
            "- The FASTA database used (name, number of sequences if available)\n",
            "- DIA-NN version and search mode (library-free vs spectral library)\n",
            "- Enzyme specificity and missed cleavages\n",
            "- Mass accuracy settings (MS1 and MS2) and whether they were auto-optimized or manually set\n",
            "- Variable modifications (e.g., methionine oxidation, N-terminal acetylation)\n",
            "- FDR threshold and match-between-runs status\n",
            "- Any other relevant search parameters (scan window, normalization, extra flags)\n",
            "Write this in third person past tense, suitable for a journal Methods section. ",
            "Cite DIA-NN (Demichev et al., Nature Methods, 2020) appropriately.\n\n"
          ) else "",
          if (nzchar(instrument_inline)) paste0(
            "## Instrument & Acquisition\n",
            "Instrument metadata was extracted from the raw files. Use this to write a complete ",
            "Sample Preparation & Data Acquisition section for the Methods. Include the instrument model, ",
            "m/z scan range, gradient length, and ion mobility range (if timsTOF). ",
            "Write in third person past tense.\n\n"
          ) else "",
          if (nzchar(methods_note)) paste0(
            "## Methodology & Reproducibility\n",
            "A detailed methodology is in `Methods_and_References.txt`. ",
            "Summarize the key steps of the analysis pipeline (normalization, statistical testing, ",
            "multiple testing correction) in a concise Methods section suitable for a publication. ",
            "Include the R code log (`Reproducibility_Code.R`) as context — cite the specific ",
            "parameters and software versions used. Include proper literature citations from the ",
            "methods file.\n\n"
          ) else "",
          "Use markdown formatting with headers. Be scientific but accessible.\n",
          "Reference specific proteins from the CSV files to support your analysis.\n",
          design_section,
          qc_inline,
          search_settings_inline,
          instrument_inline,
          tic_inline,
          excluded_inline,
          missingness_inline,
          dynamic_range_inline,
          mofa_inline,
          phospho_inline,
          "\n\n--- TOP DE PROTEINS (summary — full data in CSV) ---\n\n",
          "Number of comparisons: ", ctx$n_contrasts, "\n\n",
          ctx$contrast_text, "\n\n",
          "--- CROSS-COMPARISON PROTEINS (significant in >= 2 comparisons) ---\n",
          ctx$cross_text, "\n\n",
          "--- MOST STABLE SIGNIFICANT PROTEINS (lowest CV across replicates) ---\n",
          ctx$stable_prots_text
        )

        prompt_file <- file.path(tmp_dir, "PROMPT.md")
        # Ensure prompt is a single character string (not a list)
        if (is.list(prompt)) prompt <- paste(unlist(prompt), collapse = "")
        if (!is.character(prompt)) prompt <- as.character(prompt)
        writeLines(prompt, prompt_file)
        files_to_zip <- c(files_to_zip, prompt_file)

        incProgress(0.9, detail = "Creating zip...")
        # zip from tmp_dir so figures/ subdirectory is preserved in ZIP
        old_wd <- setwd(tmp_dir)
        on.exit(setwd(old_wd), add = TRUE)
        # Convert absolute paths to relative (strip tmp_dir prefix, handle trailing /)
        norm_tmp <- normalizePath(tmp_dir, mustWork = FALSE)
        rel_paths <- vapply(files_to_zip, function(fp) {
          nfp <- normalizePath(fp, mustWork = FALSE)
          sub(paste0(norm_tmp, "/"), "", nfp, fixed = TRUE)
        }, character(1), USE.NAMES = FALSE)
        zip(file, rel_paths)

        message(sprintf("[DE-LIMP] Claude export: %d files, prompt %d chars", length(files_to_zip), nchar(prompt)))
      })
      }, error = function(e) {
        message("[DE-LIMP] Claude export FAILED: ", e$message)
        message("[DE-LIMP] Error call: ", deparse(e$call))
        showNotification(
          paste("Export error:", e$message, "\nCheck R console for details."),
          type = "error", duration = 15)
      })
  }

  output$download_claude_prompt <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_Claude_Export_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip")
    },
    content = claude_export_content,
    contentType = "application/zip"
  )

  # --- Export for Claude (AI Chat tab) — same content as download_claude_prompt ---
  output$download_claude_prompt_chat <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_Claude_Export_", format(Sys.time(), "%Y%m%d_%H%M"), ".zip")
    },
    content = claude_export_content,
    contentType = "application/zip"
  )

  # --- AI Summary Info Modal ---
  observeEvent(input$ai_summary_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About AI Summary"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("How it works"),
        p("The AI Summary analyzes ", strong("all comparisons"), " in your experiment at once, not just the currently selected contrast. ",
          "It identifies the top differentially expressed proteins per comparison, finds proteins that are significant across multiple ",
          "comparisons (cross-comparison biomarkers), and highlights the most reproducibly measured proteins (lowest CV) as high-confidence candidates."),
        p("The AI then provides biological interpretation, discussing known functions, pathway involvement, and disease associations for the top biomarkers."),
        tags$h6("What data is sent to Google Gemini"),
        tags$ul(
          tags$li("Top significant proteins per comparison (gene names, log2 fold-changes, adjusted p-values)"),
          tags$li("Proteins significant across multiple comparisons with their fold-changes"),
          tags$li("Most stable significant proteins (lowest coefficient of variation across replicates)"),
          tags$li("Number of comparisons and significance counts (up/down)")
        ),
        tags$h6("What is NOT sent"),
        tags$ul(
          tags$li("Raw expression values or intensity data"),
          tags$li("Individual sample names or file paths"),
          tags$li("Metadata details (groups, batches, covariates)"),
          tags$li("QC statistics or run-level information")
        ),
        tags$h6("Privacy"),
        p("Data is sent to Google's Gemini API and processed according to Google's API terms of service. ",
          "No data is stored permanently by this app \u2014 uploaded files are deleted when your session ends."),
        tags$h6("API key"),
        p("You need a Google Gemini API key (enter in the sidebar). Get one free at ",
          tags$a(href = "https://aistudio.google.com/apikey", target = "_blank", "Google AI Studio"), ".")
      )
    ))
  })

  # --- Data Chat Info Modal ---
  observeEvent(input$data_chat_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " About AI Analysis"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("How it works"),
        p("Data Chat uses the Google Gemini API to provide AI-powered analysis of your proteomics data. ",
          "Your QC statistics and the top 100-800 differentially expressed proteins (scaled by dataset size) are uploaded to Gemini for context-aware responses."),
        tags$h6("What data is sent"),
        tags$ul(
          tags$li("QC statistics (precursor counts, protein counts, MS1 signal per sample)"),
          tags$li("Top 800 DE proteins with fold-changes and p-values"),
          tags$li("Your chat messages")
        ),
        tags$h6("Privacy"),
        p("Data is sent to Google's Gemini API. It is processed according to Google's API terms of service. ",
          "No data is stored permanently by this app \u2014 uploaded files are deleted when your session ends."),
        tags$h6("Plot selection integration"),
        p("If you select proteins in the volcano plot or results table, the chat knows about your selection. ",
          "The AI can also suggest proteins to highlight \u2014 look for the ",
          tags$em("'I have updated your plots'"), " message after AI responses."),
        tags$h6("API key"),
        p("You need a Google Gemini API key (enter in the sidebar). Get one free at ",
          tags$a(href = "https://aistudio.google.com/apikey", target = "_blank", "Google AI Studio"), ".")
      )
    ))
  })

  observeEvent(input$check_models, { if (nchar(input$user_api_key) < 10) { showNotification("Please enter a valid API Key first.", type="error"); return() }; withProgress(message = "Checking Google Models...", { models <- list_google_models(input$user_api_key); if (length(models) > 0 && !grepl("Error", models[1])) { showModal(modalDialog(title = "Available Models for Your Key", p("Copy one of these into the Model Name box:"), tags$textarea(paste(models, collapse="\n"), rows=10, style="width:100%;"), easyClose = TRUE)) } else { showNotification(paste("Failed to list models:", models), type="error") } }) })
  output$chat_selection_indicator <- renderText({ if (!is.null(values$plot_selected_proteins)) { paste("\u2705 Current Selection:", length(values$plot_selected_proteins), "Proteins from Plots.") } else { "\u2139\ufe0f No proteins selected in plots." } })

  observeEvent(input$summarize_data, {
    req(input$user_api_key)
    auto_prompt <- "Analyze this dataset. Identify key quality control issues (if any) by looking at the Group QC stats. Then, summarize the main biological findings from the expression data, focusing on the most significantly differentially expressed proteins."
    values$chat_history <- append(values$chat_history, list(list(role = "user", content = "(Auto-Query: Summarize & Analyze)")))
    withProgress(message = "Auto-Analyzing Dataset...", {
      if (!is.null(values$fit) && !is.null(values$y_protein)) {
        # Scale protein count to stay within Gemini's token limit (~1M tokens)
        n_samples <- ncol(values$y_protein$E)
        n_max <- if (n_samples > 200) 100 else if (n_samples > 100) 200 else if (n_samples > 50) 400 else 800
        message(sprintf("[DE-LIMP] AI data: %d proteins x %d samples (scaled from 800)", n_max, n_samples))

        df_de <- topTable(values$fit, coef=input$contrast_selector, number=n_max)

        # For large datasets (>100 samples), send group-level summary stats
        # instead of per-sample expression to stay within token limits
        if (n_samples > 100 && !is.null(values$metadata)) {
          row_idx <- match(rownames(df_de), rownames(values$y_protein$E))
          row_idx <- row_idx[!is.na(row_idx)]
          exprs_mat <- values$y_protein$E[row_idx, , drop = FALSE]
          meta <- values$metadata[values$metadata$Group != "", ]
          group_stats <- do.call(cbind, lapply(unique(meta$Group), function(g) {
            cols <- intersect(meta$File.Name[meta$Group == g], colnames(exprs_mat))
            if (length(cols) == 0) return(NULL)
            data.frame(
              setNames(list(
                rowMeans(exprs_mat[, cols, drop = FALSE], na.rm = TRUE),
                apply(exprs_mat[, cols, drop = FALSE], 1, sd, na.rm = TRUE)
              ), c(paste0("Mean_", g), paste0("SD_", g)))
            )
          }))
          df_full <- cbind(Protein = rownames(df_de), df_de, group_stats)
        } else {
          row_idx <- match(rownames(df_de), rownames(values$y_protein$E))
          row_idx <- row_idx[!is.na(row_idx)]
          df_exprs <- as.data.frame(values$y_protein$E[row_idx, ])
          df_full <- cbind(Protein = rownames(df_de), df_de, df_exprs)
        }

        incProgress(0.3, detail = "Sending data file..."); current_file_uri <- upload_csv_to_gemini(df_full, input$user_api_key)
        qc_final <- NULL; if(!is.null(values$qc_stats) && !is.null(values$metadata)) { qc_final <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>% dplyr::select(Run, Group, Precursors, Proteins, MS1_Signal) }
        # Append phospho context if phospho analysis is active
        auto_msg <- auto_prompt
        if (!is.null(values$phospho_fit) && !is.null(input$phospho_contrast_selector)) {
          phospho_ctx <- tryCatch(
            phospho_ai_context(values$phospho_fit, input$phospho_contrast_selector, values$ksea_results),
            error = function(e) ""
          )
          if (nzchar(phospho_ctx)) auto_msg <- paste0(auto_msg, phospho_ctx)
        }
        incProgress(0.7, detail = "Thinking..."); ai_reply <- ask_gemini_file_chat(auto_msg, current_file_uri, qc_final, input$user_api_key, input$model_name, values$plot_selected_proteins)
      } else { ai_reply <- "Please load data and run analysis first." }
      values$chat_history <- append(values$chat_history, list(list(role = "ai", content = ai_reply)))
    })
  })

  observeEvent(input$send_chat, {
    req(input$chat_input, input$user_api_key)
    values$chat_history <- append(values$chat_history, list(list(role = "user", content = input$chat_input)))
    withProgress(message = "Processing...", {
      if (!is.null(values$fit) && !is.null(values$y_protein)) {
        # Scale protein count to stay within Gemini's token limit (~1M tokens)
        n_samples <- ncol(values$y_protein$E)
        n_max <- if (n_samples > 200) 100 else if (n_samples > 100) 200 else if (n_samples > 50) 400 else 800
        df_de <- topTable(values$fit, coef=input$contrast_selector, number=n_max)
        if (!is.null(values$plot_selected_proteins)) { missing_ids <- setdiff(values$plot_selected_proteins, rownames(df_de)); if (length(missing_ids) > 0) { valid_missing <- intersect(missing_ids, rownames(values$fit$coefficients)); if(length(valid_missing) > 0) { df_extra <- topTable(values$fit, coef=input$contrast_selector, number=Inf)[valid_missing, ]; df_de <- rbind(df_de, df_extra) } } }
        # For large datasets (>100 samples), send group-level summary stats
        if (n_samples > 100 && !is.null(values$metadata)) {
          row_idx <- match(rownames(df_de), rownames(values$y_protein$E))
          row_idx <- row_idx[!is.na(row_idx)]
          exprs_mat <- values$y_protein$E[row_idx, , drop = FALSE]
          meta <- values$metadata[values$metadata$Group != "", ]
          group_stats <- do.call(cbind, lapply(unique(meta$Group), function(g) {
            cols <- intersect(meta$File.Name[meta$Group == g], colnames(exprs_mat))
            if (length(cols) == 0) return(NULL)
            data.frame(
              setNames(list(
                rowMeans(exprs_mat[, cols, drop = FALSE], na.rm = TRUE),
                apply(exprs_mat[, cols, drop = FALSE], 1, sd, na.rm = TRUE)
              ), c(paste0("Mean_", g), paste0("SD_", g)))
            )
          }))
          df_full <- cbind(Protein = rownames(df_de), df_de, group_stats)
        } else {
          row_idx <- match(rownames(df_de), rownames(values$y_protein$E))
          row_idx <- row_idx[!is.na(row_idx)]
          df_exprs <- as.data.frame(values$y_protein$E[row_idx, ]); df_full <- cbind(Protein = rownames(df_de), df_de, df_exprs)
        }
        incProgress(0.3, detail = "Sending data file..."); current_file_uri <- upload_csv_to_gemini(df_full, input$user_api_key)
        qc_final <- NULL; if(!is.null(values$qc_stats) && !is.null(values$metadata)) { qc_final <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>% dplyr::select(Run, Group, Precursors, Proteins, MS1_Signal) }
        # Append phospho context if phospho analysis is active
        chat_msg <- input$chat_input
        if (!is.null(values$phospho_fit) && !is.null(input$phospho_contrast_selector)) {
          phospho_ctx <- tryCatch(
            phospho_ai_context(values$phospho_fit, input$phospho_contrast_selector, values$ksea_results),
            error = function(e) ""
          )
          if (nzchar(phospho_ctx)) chat_msg <- paste0(chat_msg, phospho_ctx)
        }
        incProgress(0.7, detail = "Thinking..."); ai_reply <- ask_gemini_file_chat(chat_msg, current_file_uri, qc_final, input$user_api_key, input$model_name, values$plot_selected_proteins)
      } else { ai_reply <- "Please load data and run analysis first." }
    })

    ai_selected <- str_extract(ai_reply, "\\[\\[SELECT:.*?\\]\\]")
    if (!is.na(ai_selected)) { raw_ids <- gsub("\\[\\[SELECT:|\\]\\]", "", ai_selected); id_vec <- unlist(strsplit(raw_ids, "[,;]\\s*")); values$plot_selected_proteins <- trimws(id_vec); ai_reply <- gsub("\\[\\[SELECT:.*?\\]\\]", "", ai_reply); ai_reply <- paste0(ai_reply, "\n\n*(I have updated your plots with these highlighted proteins.)*") }
    values$chat_history <- append(values$chat_history, list(list(role = "ai", content = ai_reply))); updateTextAreaInput(session, "chat_input", value = "")
  })

  output$chat_window <- renderUI({ chat_content <- lapply(values$chat_history, function(msg) { if (msg$role == "user") { div(class = "user-msg", span(msg$content)) } else { div(class = "ai-msg", span(markdown(msg$content))) } }); div(class = "chat-container", chat_content) })

  output$download_chat_txt <- downloadHandler(
    filename = function() { req(values$chat_history); paste0("Limpa_Chat_History_", Sys.Date(), ".txt") },
    content = function(file) { req(values$chat_history); text_out <- sapply(values$chat_history, function(msg) { paste0(if(msg$role == "user") "YOU: " else "GEMINI: ", msg$content, "\n---\n") }); writeLines(unlist(text_out), file) }
  )

}
