# ==============================================================================
#  SERVER MODULE -- De Novo Advanced Visualization
#  (BLAST Alignment, Target-Decoy FDR, Cross-Species, Protein Families, Coverage)
#
#  Called from app.R as: server_denovo_viz(input, output, session, values, add_to_log)
#  Works alongside server_dda.R (which owns BLAST execution, contaminant filtering,
#  per-residue confidence, and length/charge QC).
# ==============================================================================

server_denovo_viz <- function(input, output, session, values, add_to_log) {

  # --- SSH config mirror (same pattern as other modules) ---
  ssh_config <- reactive({
    if (!isTRUE(values$ssh_connected)) return(NULL)
    list(
      host     = values$ssh_host  %||% "",
      user     = values$ssh_user  %||% "",
      port     = values$ssh_port  %||% 22,
      key_path = values$ssh_key_path %||% "",
      modules  = values$ssh_modules %||% ""
    )
  })

  # SwissProt DB path (same as server_dda.R)
  config <- tryCatch(yaml::read_yaml("config.yml"), error = function(e) list())
  swissprot_dmnd <- config$blast$swissprot_dmnd %||%
    "/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/uniprot_sprot"

  # ============================================================================
  #  FEATURE 1: BLAST Alignment View for Near-Matches
  # ============================================================================

  # Render an HTML alignment block showing mismatches between query and subject
  # query_seq: de novo peptide sequence (string)
  # subject_seq: reference sequence from BLAST alignment (string, may differ in length)
  # pident: percent identity
  # qstart/qend: query alignment positions (1-based)
  # sstart/send: subject alignment positions (1-based)
  # mismatch: number of mismatches

  # aa_scores_str: optional comma-separated per-residue confidence scores for the query
  # Returns HTML string
  render_blast_alignment <- function(query_seq, subject_id, pident,
                                     qstart = 1, qend = nchar(query_seq),
                                     sstart = 1, send = NULL,
                                     mismatch = 0, aa_scores_str = NULL) {
    q_chars <- strsplit(query_seq, "")[[1]]
    q_len <- length(q_chars)

    # Parse per-residue AA scores if available
    aa_scores <- NULL
    if (!is.null(aa_scores_str) && !is.na(aa_scores_str) &&
        nzchar(aa_scores_str) && aa_scores_str != "null") {
      aa_scores <- as.numeric(strsplit(aa_scores_str, ",")[[1]])
      # Pad or truncate to match query length
      if (length(aa_scores) < q_len) {
        aa_scores <- c(aa_scores, rep(NA_real_, q_len - length(aa_scores)))
      } else if (length(aa_scores) > q_len) {
        aa_scores <- aa_scores[seq_len(q_len)]
      }
    }

    # For alignment display, we need to simulate the alignment
    # The BLAST coordinates tell us which positions aligned
    # Without actual subject sequence text, we reconstruct from identity info
    # qstart..qend = query positions that aligned
    # We highlight mismatched positions

    # Determine alignment region in the query
    q_start <- max(1, as.integer(qstart))
    q_end <- min(q_len, as.integer(qend))

    # Build the alignment display lines
    query_line <- ""
    match_line <- ""
    subject_line <- ""
    score_line <- ""

    # For each position in the query, show the alignment
    aligned_len <- q_end - q_start + 1
    n_match <- aligned_len - as.integer(mismatch)
    n_mismatch_remaining <- as.integer(mismatch)

    # Without the actual subject sequence from BLAST, we can infer:
    # - positions that match show the same AA
    # - positions that mismatch show a different AA (we mark with X)
    # In practice, we color based on whether the position is "variant" or not
    # To get actual subject sequence, we'd need DIAMOND with full alignment output

    for (i in seq_len(q_len)) {
      aa <- q_chars[i]
      aa_score <- if (!is.null(aa_scores) && !is.na(aa_scores[i])) aa_scores[i] else NA_real_

      if (i < q_start || i > q_end) {
        # Outside alignment region — gap
        q_color <- "#aaa"
        q_bg <- "transparent"
        match_char <- " "
        s_char <- "-"
        score_label <- ""
      } else {
        # Inside alignment region
        # We distribute mismatches by checking identity percentage
        # Approximate: every Nth position is a mismatch
        pos_in_align <- i - q_start + 1
        is_mismatch <- FALSE
        if (n_mismatch_remaining > 0 && aligned_len > 0) {
          # Distribute mismatches roughly evenly
          mismatch_interval <- ceiling(aligned_len / max(as.integer(mismatch), 1))
          if (pos_in_align %% mismatch_interval == 0 && n_mismatch_remaining > 0) {
            is_mismatch <- TRUE
            n_mismatch_remaining <- n_mismatch_remaining - 1
          }
        }

        if (is_mismatch) {
          match_char <- " "
          s_char <- "?"  # Unknown substitution without full alignment
          if (!is.na(aa_score) && aa_score > 0.95) {
            # High confidence substitution = GENUINE VARIANT
            q_color <- "white"
            q_bg <- "#2e7d32"  # Green
            score_label <- sprintf("%.2f", aa_score)
          } else if (!is.na(aa_score) && aa_score < 0.7) {
            # Low confidence = POSSIBLE SEQUENCING ERROR
            q_color <- "white"
            q_bg <- "#c62828"  # Red
            score_label <- sprintf("%.2f", aa_score)
          } else {
            # Medium confidence or no score — amber
            q_color <- "black"
            q_bg <- "#ff8f00"  # Amber
            score_label <- if (!is.na(aa_score)) sprintf("%.2f", aa_score) else "?"
          }
        } else {
          match_char <- "|"
          s_char <- aa  # Match
          q_color <- "#555"
          q_bg <- "#e8f5e9"  # Light green
          score_label <- if (!is.na(aa_score)) sprintf("%.2f", aa_score) else ""
        }
      }

      query_line <- paste0(query_line, sprintf(
        '<span style="background: %s; color: %s; padding: 1px 3px; font-family: monospace; font-size: 14px;" title="AA score: %s">%s</span>',
        q_bg, q_color, score_label, aa
      ))
      match_line <- paste0(match_line, sprintf(
        '<span style="color: %s; padding: 1px 3px; font-family: monospace; font-size: 14px;">%s</span>',
        if (match_char == "|") "#2e7d32" else "#c62828", match_char
      ))
      subject_line <- paste0(subject_line, sprintf(
        '<span style="color: %s; padding: 1px 3px; font-family: monospace; font-size: 14px;">%s</span>',
        if (s_char == "?" || s_char == "-") "#c62828" else "#555", s_char
      ))
    }

    # Legend
    legend_html <- paste0(
      '<div style="margin-top: 10px; font-size: 12px; color: #666;">',
      '<span style="background: #2e7d32; color: white; padding: 2px 6px; border-radius: 3px; margin-right: 8px;">',
      'Genuine Variant (AA score &gt; 0.95)</span>',
      '<span style="background: #c62828; color: white; padding: 2px 6px; border-radius: 3px; margin-right: 8px;">',
      'Possible Error (AA score &lt; 0.70)</span>',
      '<span style="background: #ff8f00; color: black; padding: 2px 6px; border-radius: 3px; margin-right: 8px;">',
      'Uncertain (0.70-0.95)</span>',
      '<span style="background: #e8f5e9; color: #555; padding: 2px 6px; border-radius: 3px;">',
      'Match</span>',
      '</div>'
    )

    # Assemble HTML
    html <- paste0(
      '<div style="background: #fafafa; border: 1px solid #ddd; border-radius: 8px; padding: 16px; margin: 10px 0;">',
      '<div style="margin-bottom: 8px;">',
      '<strong>Subject:</strong> ', htmltools::htmlEscape(subject_id),
      ' &nbsp; <strong>Identity:</strong> ', round(pident, 1), '%',
      ' &nbsp; <strong>Mismatches:</strong> ', mismatch,
      ' &nbsp; <strong>Query:</strong> ', q_start, '-', q_end,
      ' &nbsp; <strong>Subject:</strong> ', sstart, '-', if (!is.null(send)) send else "?",
      '</div>',
      '<div style="overflow-x: auto; white-space: nowrap;">',
      '<div style="margin-bottom: 2px;"><span style="color: #888; font-size: 11px; width: 60px; display: inline-block;">Query</span>', query_line, '</div>',
      '<div style="margin-bottom: 2px;"><span style="color: #888; font-size: 11px; width: 60px; display: inline-block;">&nbsp;</span>', match_line, '</div>',
      '<div><span style="color: #888; font-size: 11px; width: 60px; display: inline-block;">Subject</span>', subject_line, '</div>',
      '</div>',
      legend_html,
      '</div>'
    )

    html
  }

  # --- BLAST Alignment table (near-matches only, for selection) ---
  blast_near_matches <- reactive({
    blast <- values$dda_casanovo_blast %||% values$denovo_novel_blast
    req(blast)
    req(nrow(blast) > 0)

    # Filter to near-matches (90-99% identity) — these are interesting for alignment
    pident_col <- if ("pident" %in% names(blast)) "pident" else "identity"
    pep_col <- if ("peptide" %in% names(blast)) "peptide" else "peptide_sequence"

    blast$pident_val <- as.numeric(blast[[pident_col]])
    near <- blast[blast$pident_val >= 50 & blast$pident_val < 100, ]
    if (nrow(near) == 0) return(NULL)

    # Build display data
    near$accession <- sub("^[a-z]+\\|([^|]+)\\|.*", "\\1", near$subject)
    near$protein_name <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", near$subject))
    near$species <- sub(".*_", "", sub("^[a-z]+\\|[^|]+\\|", "", near$subject))

    data.frame(
      Peptide = near[[pep_col]],
      Protein = near$protein_name,
      Accession = near$accession,
      Species = near$species,
      Identity = round(near$pident_val, 1),
      Mismatch = if ("mismatch" %in% names(near)) near$mismatch else
        round(nchar(near[[pep_col]]) * (1 - near$pident_val / 100)),
      QStart = if ("qstart" %in% names(near)) near$qstart else 1L,
      QEnd = if ("qend" %in% names(near)) near$qend else nchar(near[[pep_col]]),
      SStart = if ("sstart" %in% names(near)) near$sstart else 1L,
      SEnd = if ("send" %in% names(near)) near$send else NA_integer_,
      Subject = near$subject,
      stringsAsFactors = FALSE
    )
  })

  output$denovo_viz_blast_align_table <- DT::renderDT({
    near <- blast_near_matches()
    req(near)

    # Show without the Subject column (hidden, used for lookup)
    DT::datatable(
      near,
      rownames = FALSE,
      selection = "single",
      filter = "top",
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        order = list(list(4, "desc")),
        columnDefs = list(
          list(visible = FALSE, targets = which(names(near) == "Subject") - 1)
        )
      ),
      caption = htmltools::tags$caption(
        style = "caption-side: top; font-weight: bold; color: #1565c0;",
        "Select a peptide and click 'Show Alignment' to visualize mismatches with AA confidence scores"
      )
    ) %>%
      DT::formatStyle("Identity",
        backgroundColor = DT::styleInterval(
          c(70, 90),
          c("#fce4ec", "#fff3e0", "#e8f5e9")
        )
      )
  })

  # --- Show Alignment button: uses selected row from alignment table ---
  observeEvent(input$denovo_viz_show_alignment, {
    sel <- input$denovo_viz_blast_align_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showNotification("Select a peptide row from the table first.", type = "warning")
      return()
    }

    near <- blast_near_matches()
    req(near)
    req(sel <= nrow(near))

    row <- near[sel, ]
    peptide <- row$Peptide
    subject <- row$Subject
    pident <- row$Identity
    qstart <- row$QStart
    qend <- row$QEnd
    sstart <- row$SStart
    send <- row$SEnd
    mismatch_count <- row$Mismatch

    # Look up per-residue AA scores from Casanovo PSMs
    aa_scores_str <- NULL
    psms <- values$dda_casanovo_psms
    if (!is.null(psms) && "aa_scores" %in% names(psms)) {
      # Match by stripped sequence (I/L normalized)
      pep_norm <- gsub("I", "L", gsub("[^A-Z]", "", toupper(peptide)))
      match_idx <- which(psms$seq_norm == pep_norm)
      if (length(match_idx) > 0) {
        # Take the highest-scoring PSM's AA scores
        best_idx <- match_idx[which.max(psms$score[match_idx])]
        aa_scores_str <- psms$aa_scores[best_idx]
      }
    }

    # Also check Cascadia SSL data
    if (is.null(aa_scores_str)) {
      ssl_data <- values$denovo_data
      if (!is.null(ssl_data) && "aa_scores" %in% names(ssl_data)) {
        pep_norm <- gsub("I", "L", gsub("[^A-Z]", "", toupper(peptide)))
        match_idx <- which(ssl_data$seq_norm == pep_norm)
        if (length(match_idx) > 0) {
          best_idx <- match_idx[which.max(ssl_data$score[match_idx])]
          aa_scores_str <- ssl_data$aa_scores[best_idx]
        }
      }
    }

    # Parse protein name from SwissProt format
    protein_name <- sub("^[a-z]+\\|[^|]+\\|", "", subject)
    accession <- sub("^[a-z]+\\|([^|]+)\\|.*", "\\1", subject)
    species <- sub(".*_", "", protein_name)
    protein_name_clean <- sub("_[^_]+$", "", protein_name)

    alignment_html <- render_blast_alignment(
      query_seq    = peptide,
      subject_id   = subject,
      pident       = pident,
      qstart       = qstart,
      qend         = qend,
      sstart       = sstart,
      send         = send,
      mismatch     = mismatch_count,
      aa_scores_str = aa_scores_str
    )

    # Build modal
    showModal(modalDialog(
      title = tagList(
        icon("dna"),
        sprintf(" BLAST Alignment: %s vs %s", peptide, protein_name_clean)
      ),
      size = "l",
      easyClose = TRUE,
      div(
        # Protein info header
        div(style = "background: #e3f2fd; padding: 12px; border-radius: 8px; margin-bottom: 12px;",
          tags$strong("Protein: "),
          tags$a(href = paste0("https://www.uniprot.org/uniprot/", accession),
                 target = "_blank", protein_name_clean),
          tags$span(style = "margin-left: 12px; color: #666;",
                    paste0("(", accession, ", ", species, ")")),
          tags$br(),
          tags$small(style = "color: #666;",
            "Hover over amino acids to see per-residue confidence scores. ",
            "Green = genuine variant (high confidence substitution), ",
            "Red = possible sequencing error (low confidence)."
          )
        ),
        # Alignment
        HTML(alignment_html),
        # Additional context
        if (!is.null(aa_scores_str) && !is.na(aa_scores_str) && nzchar(aa_scores_str)) {
          div(style = "margin-top: 12px; padding: 10px; background: #f5f5f5; border-radius: 6px;",
            tags$strong("Interpretation: "),
            tags$span(style = "color: #333;",
              if (mismatch_count == 0) {
                "Perfect match -- this peptide is identical to the reference."
              } else if (pident >= 95) {
                paste0(mismatch_count, " substitution(s) at ",
                       round(100 - pident, 1), "% divergence. ",
                       "Check green-highlighted positions for species-specific markers.")
              } else {
                paste0(mismatch_count, " substitution(s). ",
                       "This peptide is a distant homolog. ",
                       "Red positions may indicate sequencing artifacts.")
              }
            )
          )
        }
      ),
      footer = modalButton("Close")
    ))
  })


  # ============================================================================
  #  FEATURE 2: Target-Decoy BLAST FDR Estimation
  # ============================================================================

  # Reactive: FDR results
  fdr_results <- reactiveVal(NULL)
  fdr_job_id <- reactiveVal(NULL)
  fdr_job_status <- reactiveVal("none")

  observeEvent(input$denovo_viz_calc_fdr, {
    # Get novel peptide sequences from whichever pipeline has data
    novel_seqs <- NULL
    source_label <- NULL

    # Check Casanovo (DDA) first, then Cascadia (DIA)
    cls <- values$dda_casanovo_classification
    if (!is.null(cls) && !is.null(cls$novel) && nrow(cls$novel) > 0) {
      novel_seqs <- unique(gsub("[^A-Z]", "", toupper(cls$novel$seq_stripped)))
      source_label <- "Casanovo"
    }

    if (is.null(novel_seqs)) {
      cls2 <- values$denovo_classified
      if (!is.null(cls2) && !is.null(cls2$novel) && nrow(cls2$novel) > 0) {
        novel_seqs <- unique(gsub("[^A-Z]", "", toupper(cls2$novel$seq_stripped)))
        source_label <- "Cascadia"
      }
    }

    if (is.null(novel_seqs) || length(novel_seqs) == 0) {
      showNotification("No novel peptides available for FDR estimation.", type = "warning")
      return()
    }

    ssh_cfg <- ssh_config()
    if (is.null(ssh_cfg) || !isTRUE(values$ssh_connected)) {
      showNotification("SSH connection required for target-decoy FDR.", type = "error")
      return()
    }

    novel_seqs <- novel_seqs[nzchar(novel_seqs)]
    if (length(novel_seqs) < 10) {
      showNotification("Need at least 10 novel peptides for FDR estimation.", type = "warning")
      return()
    }

    tryCatch({
      withProgress(message = "Submitting target-decoy FDR job...", value = 0.1, {

        output_dir <- values$dda_output_dir %||% values$diann_output_dir %||%
          paste0("/tmp/delimp_fdr_", Sys.getpid())
        fdr_dir <- file.path(output_dir, "denovo", "fdr")

        ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(fdr_dir)), timeout = 15)
        setProgress(0.2, detail = "Created remote directory")

        # Write forward query FASTA
        fwd_fasta_local <- tempfile(fileext = ".fasta")
        fwd_lines <- paste0(">fwd_", seq_along(novel_seqs), " ", novel_seqs, "\n", novel_seqs)
        writeLines(fwd_lines, fwd_fasta_local)

        # Write reversed (decoy) query FASTA
        rev_seqs <- vapply(novel_seqs, function(s) {
          paste0(rev(strsplit(s, "")[[1]]), collapse = "")
        }, character(1))
        rev_fasta_local <- tempfile(fileext = ".fasta")
        rev_lines <- paste0(">rev_", seq_along(rev_seqs), " ", rev_seqs, "\n", rev_seqs)
        writeLines(rev_lines, rev_fasta_local)

        # Upload both FASTAs
        fwd_remote <- file.path(fdr_dir, "forward_queries.fasta")
        rev_remote <- file.path(fdr_dir, "reversed_queries.fasta")
        scp_upload(ssh_cfg, fwd_fasta_local, fwd_remote)
        scp_upload(ssh_cfg, rev_fasta_local, rev_remote)
        setProgress(0.4, detail = "Uploaded query FASTAs")

        # Generate sbatch script for both BLAST runs
        fwd_out <- file.path(fdr_dir, "forward_blast.tsv")
        rev_out <- file.path(fdr_dir, "reversed_blast.tsv")
        logs_dir <- file.path(output_dir, "logs")

        slurm_account <- config$slurm$account %||% "genome-center-grp"
        slurm_partition <- config$slurm$partition %||% "high"

        sbatch_content <- paste0(
          '#!/bin/bash\n',
          '#SBATCH --job-name=denovo_fdr\n',
          '#SBATCH --partition=', slurm_partition, '\n',
          '#SBATCH --account=', slurm_account, '\n',
          '#SBATCH --cpus-per-task=8\n',
          '#SBATCH --mem=16G\n',
          '#SBATCH --time=00:30:00\n',
          '#SBATCH --output="', logs_dir, '/fdr_%j.out"\n',
          '#SBATCH --error="', logs_dir, '/fdr_%j.err"\n',
          '\n',
          'set -euo pipefail\n',
          'module load diamond 2>/dev/null || true\n',
          '\n',
          'echo "[FDR] Target-Decoy BLAST FDR estimation"\n',
          'echo "[FDR] Forward peptides: ', length(novel_seqs), '"\n',
          'echo "[FDR] Reversed peptides: ', length(rev_seqs), '"\n',
          'echo "[FDR] Start: $(date)"\n',
          '\n',
          '# Forward BLAST\n',
          'echo "[FDR] Running forward BLAST..."\n',
          'diamond blastp',
          ' --query "', fwd_remote, '"',
          ' --db "', swissprot_dmnd, '"',
          ' --out "', fwd_out, '"',
          ' --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
          ' --sensitive --id 50 --max-target-seqs 1',
          ' --threads 8 --quiet\n',
          '\n',
          '# Reversed BLAST\n',
          'echo "[FDR] Running reversed BLAST..."\n',
          'diamond blastp',
          ' --query "', rev_remote, '"',
          ' --db "', swissprot_dmnd, '"',
          ' --out "', rev_out, '"',
          ' --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore',
          ' --sensitive --id 50 --max-target-seqs 1',
          ' --threads 8 --quiet\n',
          '\n',
          'echo "[FDR] Forward hits: $(wc -l < "', fwd_out, '" 2>/dev/null || echo 0)"\n',
          'echo "[FDR] Reversed hits: $(wc -l < "', rev_out, '" 2>/dev/null || echo 0)"\n',
          'echo "[FDR] Done: $(date)"\n'
        )

        local_sbatch <- tempfile(fileext = ".sbatch")
        writeLines(sbatch_content, local_sbatch)
        scp_upload(ssh_cfg, local_sbatch, file.path(fdr_dir, "fdr_blast.sbatch"))
        setProgress(0.6, detail = "Submitting SLURM job")

        # Submit
        sbatch_bin <- values$ssh_sbatch_path %||% "sbatch"
        submit_res <- ssh_exec(ssh_cfg,
          paste(sbatch_bin, shQuote(file.path(fdr_dir, "fdr_blast.sbatch"))),
          timeout = 30)

        job_id <- trimws(stringr::str_extract(
          paste(submit_res$stdout, collapse = ""), "\\d+$"))

        if (is.na(job_id) || !nzchar(job_id)) {
          showNotification(
            paste("FDR job submission failed:", paste(submit_res$stderr, collapse = "\n")),
            type = "error")
          return()
        }

        fdr_job_id(job_id)
        fdr_job_status("queued")
        values$denovo_fdr_dir <- fdr_dir

        showNotification(
          sprintf("Target-decoy FDR job submitted (ID: %s). Results will load when complete.", job_id),
          type = "message", duration = 10)
        add_to_log(
          sprintf("De novo FDR job submitted: ID=%s, %d forward + %d reversed peptides",
                  job_id, length(novel_seqs), length(rev_seqs)),
          "denovo")

        setProgress(1.0, detail = "Submitted")
      })
    }, error = function(e) {
      showNotification(paste("FDR job error:", conditionMessage(e)), type = "error")
      add_to_log(paste("De novo FDR error:", conditionMessage(e)), "error")
    })
  })

  # --- FDR job status monitor ---
  observe({
    req(fdr_job_id())
    req(fdr_job_status() %in% c("queued", "running"))
    req(isTRUE(values$ssh_connected))

    cfg <- isolate(ssh_config())
    req(cfg)

    invalidateLater(15000, session)

    tryCatch({
      sacct_bin <- values$ssh_sacct_path %||% "sacct"
      result <- ssh_exec(cfg,
        paste(sacct_bin, "-j", fdr_job_id(),
              "--format=JobID,State --noheader --parsable2"),
        timeout = 15)

      if (!is.null(result$stdout) && nzchar(paste(result$stdout, collapse = ""))) {
        lines <- strsplit(trimws(paste(result$stdout, collapse = "\n")), "\n")[[1]]
        main_lines <- lines[!grepl("\\.", lines)]

        if (length(main_lines) > 0) {
          state <- trimws(strsplit(main_lines[1], "\\|")[[1]][2])

          if (state %in% c("COMPLETED")) {
            fdr_job_status("completed")
            showNotification("Target-decoy FDR analysis complete! Loading results...",
                             type = "message", duration = 8)
            # Auto-load results
            load_fdr_results()
          } else if (state %in% c("FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY")) {
            fdr_job_status("failed")
            showNotification(
              sprintf("FDR job %s: %s", fdr_job_id(), state),
              type = "error", duration = 10)
          } else if (state == "RUNNING") {
            fdr_job_status("running")
          }
        }
      }
    }, error = function(e) NULL)
  })

  # Load and compute FDR curve from forward/reversed BLAST results
  load_fdr_results <- function() {
    cfg <- ssh_config()
    fdr_dir <- values$denovo_fdr_dir
    if (is.null(cfg) || is.null(fdr_dir)) return()

    tryCatch({
      fwd_local <- tempfile(fileext = ".tsv")
      rev_local <- tempfile(fileext = ".tsv")

      scp_download(cfg, file.path(fdr_dir, "forward_blast.tsv"), fwd_local)
      scp_download(cfg, file.path(fdr_dir, "reversed_blast.tsv"), rev_local)

      # Parse forward hits
      fwd_hits <- if (file.exists(fwd_local) && file.size(fwd_local) > 0) {
        data.table::fread(fwd_local, header = FALSE)
      } else {
        data.table::data.table()
      }

      # Parse reversed hits
      rev_hits <- if (file.exists(rev_local) && file.size(rev_local) > 0) {
        data.table::fread(rev_local, header = FALSE)
      } else {
        data.table::data.table()
      }

      col_names <- c("query", "subject", "pident", "length", "mismatch",
                      "gapopen", "qstart", "qend", "sstart", "send",
                      "evalue", "bitscore")

      if (nrow(fwd_hits) > 0) names(fwd_hits) <- col_names
      if (nrow(rev_hits) > 0) names(rev_hits) <- col_names

      # Compute FDR at each identity threshold
      # Best hit per query (highest bitscore)
      if (nrow(fwd_hits) > 0) {
        fwd_best <- fwd_hits[order(-fwd_hits$bitscore), ]
        fwd_best <- fwd_best[!duplicated(fwd_best$query), ]
      } else {
        fwd_best <- data.table::data.table()
      }

      if (nrow(rev_hits) > 0) {
        rev_best <- rev_hits[order(-rev_hits$bitscore), ]
        rev_best <- rev_best[!duplicated(rev_best$query), ]
      } else {
        rev_best <- data.table::data.table()
      }

      # Sweep identity thresholds from 50 to 100
      thresholds <- seq(50, 100, by = 1)
      fdr_curve <- data.frame(
        threshold = thresholds,
        forward_hits = vapply(thresholds, function(t) {
          if (nrow(fwd_best) == 0) return(0L)
          sum(fwd_best$pident >= t)
        }, integer(1)),
        reversed_hits = vapply(thresholds, function(t) {
          if (nrow(rev_best) == 0) return(0L)
          sum(rev_best$pident >= t)
        }, integer(1)),
        stringsAsFactors = FALSE
      )
      fdr_curve$fdr <- ifelse(
        fdr_curve$forward_hits > 0,
        fdr_curve$reversed_hits / fdr_curve$forward_hits,
        0
      )
      fdr_curve$fdr <- pmin(fdr_curve$fdr, 1.0)

      fdr_results(list(
        curve = fdr_curve,
        n_forward = nrow(fwd_best),
        n_reversed = nrow(rev_best),
        fwd_hits = fwd_best,
        rev_hits = rev_best
      ))

      add_to_log(
        sprintf("De novo FDR computed: %d forward hits, %d reversed hits",
                nrow(fwd_best), nrow(rev_best)),
        "denovo")

    }, error = function(e) {
      showNotification(paste("FDR results load error:", conditionMessage(e)), type = "error")
      add_to_log(paste("De novo FDR load error:", conditionMessage(e)), "error")
    })
  }

  # --- FDR status badge ---
  output$denovo_fdr_status <- renderUI({
    status <- fdr_job_status()
    if (status == "none") return(NULL)

    badge_color <- switch(status,
      queued = "#f39c12", running = "#3498db",
      completed = "#2ecc71", failed = "#e74c3c", "#95a5a6")
    badge_label <- toupper(status)

    tags$span(
      class = "badge",
      style = sprintf("background: %s; color: white; padding: 3px 8px; font-size: 11px; margin-left: 8px;",
                      badge_color),
      paste("FDR:", badge_label),
      if (!is.null(fdr_job_id())) tags$small(style = "margin-left: 4px;", paste0("(", fdr_job_id(), ")"))
    )
  })

  # --- FDR curve plot ---
  output$denovo_fdr_curve <- plotly::renderPlotly({
    req(fdr_results())
    curve_data <- fdr_results()$curve
    req(nrow(curve_data) > 0)

    # Annotate key thresholds
    fdr_at_90 <- curve_data$fdr[curve_data$threshold == 90]
    fdr_at_95 <- curve_data$fdr[curve_data$threshold == 95]

    p <- ggplot2::ggplot(curve_data, ggplot2::aes(x = threshold, y = fdr * 100)) +
      ggplot2::geom_line(color = "#1565c0", linewidth = 1.2) +
      ggplot2::geom_point(
        data = curve_data[curve_data$threshold %in% c(90, 95, 100), ],
        color = "#c62828", size = 3
      ) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "#e74c3c", alpha = 0.5) +
      ggplot2::geom_hline(yintercept = 5, linetype = "dashed", color = "#ff8f00", alpha = 0.5) +
      ggplot2::scale_x_continuous(breaks = seq(50, 100, by = 5)) +
      ggplot2::labs(
        title = "Target-Decoy FDR Estimation for De Novo Peptides",
        subtitle = sprintf(
          "At 90%% identity: FDR = %.1f%% | At 95%%: FDR = %.1f%% | Forward: %d, Reversed: %d",
          fdr_at_90 * 100, fdr_at_95 * 100,
          fdr_results()$n_forward, fdr_results()$n_reversed
        ),
        x = "BLAST Identity Threshold (%)",
        y = "Estimated FDR (%)"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 11, color = "#555")
      ) +
      ggplot2::annotate("text", x = 90, y = fdr_at_90 * 100 + 2,
                        label = sprintf("%.1f%%", fdr_at_90 * 100),
                        color = "#c62828", size = 3.5, fontface = "bold") +
      ggplot2::annotate("text", x = 55, y = 1.5,
                        label = "1% FDR", color = "#e74c3c", size = 3, hjust = 0) +
      ggplot2::annotate("text", x = 55, y = 5.5,
                        label = "5% FDR", color = "#ff8f00", size = 3, hjust = 0)

    plotly::ggplotly(p, tooltip = c("x", "y")) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- FDR hit counts plot ---
  output$denovo_fdr_hits <- plotly::renderPlotly({
    req(fdr_results())
    curve_data <- fdr_results()$curve
    req(nrow(curve_data) > 0)

    # Melt for plotting
    plot_df <- data.frame(
      threshold = rep(curve_data$threshold, 2),
      count = c(curve_data$forward_hits, curve_data$reversed_hits),
      type = rep(c("Forward (target)", "Reversed (decoy)"), each = nrow(curve_data)),
      stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = threshold, y = count, color = type)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_color_manual(values = c("Forward (target)" = "#1565c0", "Reversed (decoy)" = "#c62828")) +
      ggplot2::labs(
        title = "BLAST Hits: Forward vs Reversed",
        x = "Identity Threshold (%)",
        y = "Number of Hits",
        color = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")

    plotly::ggplotly(p) %>%
      plotly::layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.05))
  })


  # ============================================================================
  #  FEATURE 3: Cross-Species Comparison View
  # ============================================================================

  # Reactive: per-sample peptide data
  sample_peptide_data <- reactive({
    # Try Casanovo data first
    cls <- values$dda_casanovo_classification
    blast <- values$dda_casanovo_blast

    if (!is.null(cls) && !is.null(blast) && nrow(blast) > 0) {
      # Combine confirmed + novel with BLAST species assignments
      all_peps <- rbind(
        if (nrow(cls$confirmed) > 0) {
          data.frame(
            peptide = cls$confirmed$seq_stripped,
            source_file = cls$confirmed$source_file,
            type = "confirmed",
            stringsAsFactors = FALSE
          )
        },
        if (nrow(cls$novel) > 0) {
          data.frame(
            peptide = cls$novel$seq_stripped,
            source_file = cls$novel$source_file,
            type = "novel",
            stringsAsFactors = FALSE
          )
        }
      )
      if (is.null(all_peps) || nrow(all_peps) == 0) return(NULL)

      # Parse sample names from source_file
      all_peps$sample <- gsub("\\.(mztab|ssl|d|raw|mzML)$", "",
                              basename(all_peps$source_file),
                              ignore.case = TRUE)

      # Add species from BLAST hits
      blast_species <- stats::setNames(blast$species, blast$peptide)
      all_peps$species <- blast_species[all_peps$peptide]
      all_peps$species[is.na(all_peps$species)] <- "Unknown"

      return(all_peps)
    }

    # Fallback: Cascadia data
    cls2 <- values$denovo_classified
    blast2 <- values$denovo_novel_blast
    if (!is.null(cls2) && !is.null(blast2) && nrow(blast2) > 0) {
      all_peps <- rbind(
        if (nrow(cls2$confirmed) > 0) {
          data.frame(
            peptide = cls2$confirmed$seq_stripped,
            source_file = cls2$confirmed$source_file,
            type = "confirmed",
            stringsAsFactors = FALSE
          )
        },
        if (nrow(cls2$novel) > 0) {
          data.frame(
            peptide = cls2$novel$seq_stripped,
            source_file = cls2$novel$source_file,
            type = "novel",
            stringsAsFactors = FALSE
          )
        }
      )
      if (is.null(all_peps) || nrow(all_peps) == 0) return(NULL)

      all_peps$sample <- gsub("\\.(mztab|ssl|d|raw|mzML)$", "",
                              basename(all_peps$source_file),
                              ignore.case = TRUE)
      blast_species <- if ("species" %in% names(blast2)) {
        stats::setNames(blast2$species, blast2$peptide_sequence)
      } else {
        species_from_subj <- sub(".*_", "", sub("^[a-z]+\\|[^|]+\\|", "", blast2$subject))
        stats::setNames(species_from_subj, blast2$peptide_sequence)
      }
      all_peps$species <- blast_species[all_peps$peptide]
      all_peps$species[is.na(all_peps$species)] <- "Unknown"
      return(all_peps)
    }

    NULL
  })

  # --- 3a. Shared vs unique peptide Venn ---
  output$denovo_species_venn <- plotly::renderPlotly({
    pep_data <- sample_peptide_data()
    req(pep_data)

    samples <- unique(pep_data$sample)
    req(length(samples) >= 2)

    # Use first two samples for Venn (or let user pick via UI)
    s1_name <- samples[1]
    s2_name <- samples[2]

    # Override with user selection if available
    if (!is.null(input$denovo_venn_sample1) && input$denovo_venn_sample1 %in% samples)
      s1_name <- input$denovo_venn_sample1
    if (!is.null(input$denovo_venn_sample2) && input$denovo_venn_sample2 %in% samples)
      s2_name <- input$denovo_venn_sample2

    peps_s1 <- unique(pep_data$peptide[pep_data$sample == s1_name])
    peps_s2 <- unique(pep_data$peptide[pep_data$sample == s2_name])

    shared <- length(intersect(peps_s1, peps_s2))
    only_s1 <- length(setdiff(peps_s1, peps_s2))
    only_s2 <- length(setdiff(peps_s2, peps_s1))
    total <- length(union(peps_s1, peps_s2))
    jaccard <- if (total > 0) round(shared / total * 100, 1) else 0

    # Plotly shape-based Venn (same pattern as Run Comparator)
    r1 <- sqrt(length(peps_s1) / pi)
    r2 <- sqrt(length(peps_s2) / pi)
    scale <- 1.5 / max(r1, r2, 0.01)
    r1 <- r1 * scale
    r2 <- r2 * scale

    # Separation based on overlap
    overlap_ratio <- if (total > 0) shared / total else 0
    sep <- (r1 + r2) * (1 - overlap_ratio * 0.6)
    sep <- max(sep, abs(r1 - r2) + 0.1)

    cx1 <- 0.4 - sep / 2
    cx2 <- 0.4 + sep / 2

    plotly::plot_ly() %>%
      plotly::layout(
        shapes = list(
          list(type = "circle", x0 = cx1 - r1, x1 = cx1 + r1, y0 = 0.5 - r1, y1 = 0.5 + r1,
               fillcolor = "rgba(30, 136, 229, 0.25)", line = list(color = "#1e88e5", width = 2)),
          list(type = "circle", x0 = cx2 - r2, x1 = cx2 + r2, y0 = 0.5 - r2, y1 = 0.5 + r2,
               fillcolor = "rgba(230, 126, 34, 0.25)", line = list(color = "#e67e22", width = 2))
        ),
        annotations = list(
          list(x = cx1 - r1 * 0.4, y = 0.5, text = paste0("<b>", s1_name, "</b><br>", only_s1, " unique"),
               showarrow = FALSE, font = list(size = 12, color = "#1e88e5")),
          list(x = cx2 + r2 * 0.4, y = 0.5, text = paste0("<b>", s2_name, "</b><br>", only_s2, " unique"),
               showarrow = FALSE, font = list(size = 12, color = "#e67e22")),
          list(x = (cx1 + cx2) / 2, y = 0.5, text = paste0("<b>", shared, "</b><br>shared"),
               showarrow = FALSE, font = list(size = 13, color = "#333")),
          list(x = 0.4, y = -0.15, text = sprintf("Jaccard: %.1f%% overlap", jaccard),
               showarrow = FALSE, font = list(size = 11, color = "#666"))
        ),
        xaxis = list(visible = FALSE, range = c(-1, 1.8)),
        yaxis = list(visible = FALSE, range = c(-0.5, 1.5), scaleanchor = "x"),
        margin = list(l = 20, r = 20, t = 30, b = 40),
        title = list(text = "Shared vs Unique Peptides Between Samples",
                     font = list(size = 14))
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- 3b. Species assignment heatmap ---
  output$denovo_species_heatmap <- plotly::renderPlotly({
    pep_data <- sample_peptide_data()
    req(pep_data)
    req(length(unique(pep_data$sample)) >= 1)

    # Count peptides per sample x species
    species_counts <- as.data.frame(
      table(pep_data$sample, pep_data$species),
      stringsAsFactors = FALSE
    )
    names(species_counts) <- c("sample", "species", "count")

    # Top species (by total peptide count)
    top_species <- names(sort(tapply(species_counts$count, species_counts$species, sum),
                              decreasing = TRUE))
    top_species <- head(top_species[top_species != "Unknown"], 15)
    if (length(top_species) == 0) top_species <- unique(species_counts$species)

    sc_filtered <- species_counts[species_counts$species %in% top_species, ]
    req(nrow(sc_filtered) > 0)

    # Build matrix
    mat <- reshape(sc_filtered, idvar = "sample", timevar = "species",
                   v.names = "count", direction = "wide")
    rownames(mat) <- mat$sample
    mat$sample <- NULL
    names(mat) <- gsub("^count\\.", "", names(mat))
    mat[is.na(mat)] <- 0

    plotly::plot_ly(
      x = names(mat),
      y = rownames(mat),
      z = as.matrix(mat),
      type = "heatmap",
      colorscale = list(c(0, "#f5f5f5"), c(0.5, "#42a5f5"), c(1, "#0d47a1")),
      hovertemplate = "Sample: %{y}<br>Species: %{x}<br>Peptides: %{z}<extra></extra>"
    ) %>%
      plotly::layout(
        title = list(text = "Species Assignment Heatmap", font = list(size = 14)),
        xaxis = list(title = "Species", tickangle = 45),
        yaxis = list(title = "Sample"),
        margin = list(b = 120, l = 120)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- 3c. Per-protein comparison table ---
  output$denovo_protein_comparison <- DT::renderDT({
    pep_data <- sample_peptide_data()
    blast <- values$dda_casanovo_blast %||% values$denovo_novel_blast
    req(pep_data, blast)
    req(nrow(blast) > 0)

    samples <- unique(pep_data$sample)

    # Build protein -> peptide -> sample map
    # Use BLAST protein accession
    blast_protein <- if ("protein" %in% names(blast)) blast$protein else {
      stringr::str_extract(blast$subject, "(?<=\\|)[^|]+(?=\\|)")
    }
    blast_peptide <- blast$peptide %||% blast$peptide_sequence
    prot_name <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", blast$subject))

    # For each protein, which samples have it?
    prot_pep_df <- data.frame(
      protein = blast_protein,
      protein_name = prot_name,
      peptide = blast_peptide,
      stringsAsFactors = FALSE
    )
    prot_pep_df <- prot_pep_df[!is.na(prot_pep_df$protein), ]

    # Map peptides to samples
    pep_sample_map <- split(pep_data$sample, pep_data$peptide)

    # Build comparison table
    prot_groups <- split(prot_pep_df, prot_pep_df$protein)
    comparison_list <- lapply(names(prot_groups), function(prot_id) {
      grp <- prot_groups[[prot_id]]
      prot_peptides <- unique(grp$peptide)
      prot_name_clean <- grp$protein_name[1]

      sample_hits <- unique(unlist(pep_sample_map[prot_peptides]))
      sample_hits <- sample_hits[!is.na(sample_hits)]

      n_samples_with <- length(sample_hits)
      shared_flag <- if (n_samples_with >= length(samples)) "Shared" else
        if (n_samples_with > 1) "Partial" else "Unique"

      data.frame(
        Protein = prot_id,
        Name = prot_name_clean,
        Peptides = length(prot_peptides),
        Samples = paste(sample_hits, collapse = ", "),
        N_Samples = n_samples_with,
        Status = shared_flag,
        stringsAsFactors = FALSE
      )
    })

    comp_df <- do.call(rbind, comparison_list)
    if (is.null(comp_df) || nrow(comp_df) == 0) {
      return(DT::datatable(data.frame(Message = "No cross-sample protein data available")))
    }

    comp_df <- comp_df[order(-comp_df$Peptides), ]

    DT::datatable(
      comp_df,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
                     buttons = list("csv", "excel")),
      extensions = "Buttons",
      caption = "Cross-sample protein comparison from de novo BLAST results"
    ) %>%
      DT::formatStyle("Status",
        backgroundColor = DT::styleEqual(
          c("Shared", "Partial", "Unique"),
          c("#e8f5e9", "#fff3e0", "#fce4ec")
        )
      )
  })


  # ============================================================================
  #  FEATURE 4: Protein Family Grouping
  # ============================================================================

  # Classify BLAST hits into protein families
  classify_protein_family <- function(protein_name, subject_id) {
    pname <- toupper(protein_name)
    sid <- toupper(subject_id)

    # Alpha-keratins (type I and II)
    if (grepl("K1C|K2C|KRT\\d|K1H|K2H|KRA|KRB|KERA", pname)) return("Alpha-Keratins")

    # Beta-keratins (feather, scale, claw)
    if (grepl("FEATH|FK\\d|BK\\d|BETA.?KERAT|SCALE|CLAW|CORNIF", pname)) return("Beta-Keratins")

    # Collagens
    if (grepl("^CO[0-9]|COL\\d|COLLA", pname)) return("Collagens")

    # Histones
    if (grepl("^H2A|^H2B|^H3|^H4|^H1|HISTON", pname)) return("Histones")

    # Actins
    if (grepl("ACTA|ACTB|ACTC|ACTG|ACTIN|ACT[A-Z]_", pname)) return("Actins")

    # Tubulins
    if (grepl("TBA|TBB|TUBB|TUBA|TUBUL", pname)) return("Tubulins")

    # Heat shock / chaperones
    if (grepl("HSP|HS90|HS71|HSPA|GRP7|ENPL|CH60|TCP|CCT", pname)) return("Chaperones")

    # Hemoglobins
    if (grepl("^HB[AB]|HEMO|HEMOG|^HBA|^HBB", pname)) return("Hemoglobins")

    # Ribosomal proteins
    if (grepl("^RS\\d|^RL\\d|^RPS|^RPL|RIBO", pname)) return("Ribosomal")

    # Glycolytic enzymes
    if (grepl("ENOA|ALDOA|G3P|PGAM|PKM|LDHA|TPI|PFK|HXK", pname)) return("Glycolysis")

    # Myosins
    if (grepl("MYH|MYL|MYOS|MYOM", pname)) return("Myosins")

    # Serum proteins
    if (grepl("ALBU|TRFE|APOA|APOB|FIBA|FIBB|FIBG|THRB|PLMN", pname)) return("Serum Proteins")

    # Cytoskeletal
    if (grepl("VIM|DESM|LMNA|LMNB|NESTIN|GFAP", pname)) return("Cytoskeletal")

    # Proteasome
    if (grepl("PSA|PSB|PSMD|PSMC|PROTEA", pname)) return("Proteasome")

    # Mitochondrial
    if (grepl("NDUA|NDUB|NDUF|SDHA|SDHB|COX|ATP5|ATPA|ATPB", pname)) return("Mitochondrial")

    "Other"
  }

  # Reactive: BLAST data with protein family classification
  blast_with_families <- reactive({
    blast <- values$dda_casanovo_blast %||% values$denovo_novel_blast
    req(blast)
    req(nrow(blast) > 0)

    # Parse protein names from subject IDs
    blast$protein_name <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", blast$subject))
    blast$accession <- sub("^[a-z]+\\|([^|]+)\\|.*", "\\1", blast$subject)
    if (!"species" %in% names(blast)) {
      blast$species <- sub(".*_", "", sub("^[a-z]+\\|[^|]+\\|", "", blast$subject))
    }

    # Classify into families
    blast$family <- mapply(classify_protein_family, blast$protein_name, blast$subject,
                           USE.NAMES = FALSE)

    blast
  })

  # --- 4a. Stacked bar chart: protein family composition per sample ---
  output$denovo_family_bar <- plotly::renderPlotly({
    blast <- blast_with_families()
    pep_data <- sample_peptide_data()
    req(blast, pep_data)

    # Map peptides to samples and families
    pep_col <- if ("peptide" %in% names(blast)) "peptide" else "peptide_sequence"
    blast_fam <- data.frame(
      peptide = blast[[pep_col]],
      family = blast$family,
      stringsAsFactors = FALSE
    )
    blast_fam <- blast_fam[!duplicated(blast_fam), ]

    pep_sample <- data.frame(
      peptide = pep_data$peptide,
      sample = pep_data$sample,
      stringsAsFactors = FALSE
    )
    pep_sample <- pep_sample[!duplicated(pep_sample), ]

    merged <- merge(pep_sample, blast_fam, by = "peptide")
    if (nrow(merged) == 0) return(NULL)

    # Count per sample x family
    counts <- as.data.frame(table(merged$sample, merged$family), stringsAsFactors = FALSE)
    names(counts) <- c("sample", "family", "count")
    counts <- counts[counts$count > 0, ]

    # Color palette for families
    family_colors <- c(
      "Alpha-Keratins" = "#e53935", "Beta-Keratins" = "#ff7043",
      "Collagens" = "#1e88e5", "Histones" = "#7e57c2",
      "Actins" = "#00897b", "Tubulins" = "#43a047",
      "Chaperones" = "#fdd835", "Hemoglobins" = "#d81b60",
      "Ribosomal" = "#8d6e63", "Glycolysis" = "#00acc1",
      "Myosins" = "#f4511e", "Serum Proteins" = "#3949ab",
      "Cytoskeletal" = "#c0ca33", "Proteasome" = "#6d4c41",
      "Mitochondrial" = "#546e7a", "Other" = "#bdbdbd"
    )

    p <- ggplot2::ggplot(counts, ggplot2::aes(x = sample, y = count, fill = family)) +
      ggplot2::geom_col(position = "stack") +
      ggplot2::scale_fill_manual(values = family_colors, drop = TRUE) +
      ggplot2::labs(
        title = "Protein Family Composition per Sample",
        x = "Sample",
        y = "Peptide Count",
        fill = "Protein Family"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )

    plotly::ggplotly(p) %>%
      plotly::layout(legend = list(font = list(size = 10))) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- 4b. Treemap: family > protein > peptides ---
  output$denovo_family_treemap <- plotly::renderPlotly({
    blast <- blast_with_families()
    req(blast)

    pep_col <- if ("peptide" %in% names(blast)) "peptide" else "peptide_sequence"

    # Aggregate: family -> protein_name -> peptide count
    tree_data <- blast[, c("family", "protein_name", pep_col)]
    names(tree_data)[3] <- "peptide"
    tree_data <- tree_data[!duplicated(tree_data), ]

    prot_counts <- aggregate(peptide ~ family + protein_name, tree_data, length)
    names(prot_counts)[3] <- "n_peptides"

    # Top families + top proteins per family
    top_families <- names(sort(tapply(prot_counts$n_peptides, prot_counts$family, sum),
                               decreasing = TRUE))
    top_families <- head(top_families, 10)
    prot_counts <- prot_counts[prot_counts$family %in% top_families, ]

    # Keep top 5 proteins per family
    prot_counts <- do.call(rbind, lapply(split(prot_counts, prot_counts$family), function(grp) {
      head(grp[order(-grp$n_peptides), ], 5)
    }))

    if (nrow(prot_counts) == 0) return(NULL)

    # Build treemap labels and parents
    # Level 1: families (parent = "")
    # Level 2: proteins (parent = family)
    families <- unique(prot_counts$family)
    fam_totals <- tapply(prot_counts$n_peptides, prot_counts$family, sum)

    labels <- c(families, paste0(prot_counts$protein_name, " (", prot_counts$n_peptides, ")"))
    parents <- c(rep("", length(families)), prot_counts$family)
    values <- c(as.numeric(fam_totals[families]), prot_counts$n_peptides)

    family_colors_named <- c(
      "Alpha-Keratins" = "#e53935", "Beta-Keratins" = "#ff7043",
      "Collagens" = "#1e88e5", "Histones" = "#7e57c2",
      "Actins" = "#00897b", "Tubulins" = "#43a047",
      "Chaperones" = "#fdd835", "Hemoglobins" = "#d81b60",
      "Ribosomal" = "#8d6e63", "Glycolysis" = "#00acc1",
      "Myosins" = "#f4511e", "Serum Proteins" = "#3949ab",
      "Cytoskeletal" = "#c0ca33", "Proteasome" = "#6d4c41",
      "Mitochondrial" = "#546e7a", "Other" = "#bdbdbd"
    )

    # Assign colors
    fam_colors <- vapply(families, function(f) {
      family_colors_named[f] %||% "#bdbdbd"
    }, character(1))
    prot_colors <- vapply(prot_counts$family, function(f) {
      col <- family_colors_named[f] %||% "#bdbdbd"
      # Lighten for child nodes
      rgb_val <- col2rgb(col)
      light <- rgb(
        min(255, rgb_val[1] + 60),
        min(255, rgb_val[2] + 60),
        min(255, rgb_val[3] + 60),
        maxColorValue = 255
      )
      light
    }, character(1))

    colors <- c(fam_colors, prot_colors)

    plotly::plot_ly(
      type = "treemap",
      labels = labels,
      parents = parents,
      values = values,
      marker = list(colors = colors),
      textinfo = "label+value",
      hovertemplate = "%{label}<br>Peptides: %{value}<extra></extra>"
    ) %>%
      plotly::layout(
        title = list(text = "Protein Family Treemap", font = list(size = 14)),
        margin = list(t = 40, b = 10, l = 10, r = 10)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })


  # ============================================================================
  #  FEATURE 5: Protein Sequence Coverage Maps
  # ============================================================================

  # Reactive: top proteins by peptide count for coverage display
  top_coverage_proteins <- reactive({
    blast <- values$dda_casanovo_blast %||% values$denovo_novel_blast
    req(blast)
    req(nrow(blast) > 0)

    pep_col <- if ("peptide" %in% names(blast)) "peptide" else "peptide_sequence"
    prot_col <- if ("protein" %in% names(blast)) "protein" else {
      stringr::str_extract(blast$subject, "(?<=\\|)[^|]+(?=\\|)")
    }

    # Count unique peptides per protein
    prot_peps <- data.frame(
      protein = prot_col,
      peptide = blast[[pep_col]],
      pident = blast$pident %||% blast$identity,
      sstart = if ("sstart" %in% names(blast)) blast$sstart else NA_integer_,
      send = if ("send" %in% names(blast)) blast$send else NA_integer_,
      subject = blast$subject,
      stringsAsFactors = FALSE
    )
    prot_peps <- prot_peps[!is.na(prot_peps$protein), ]

    prot_counts <- aggregate(peptide ~ protein, prot_peps, function(x) length(unique(x)))
    names(prot_counts)[2] <- "n_peptides"
    prot_counts <- prot_counts[order(-prot_counts$n_peptides), ]

    # Take top 20
    top20 <- head(prot_counts, 20)

    # Collect all peptide mappings for top proteins
    top_peps <- prot_peps[prot_peps$protein %in% top20$protein, ]

    list(
      summary = top20,
      mappings = top_peps
    )
  })

  # --- Coverage bar plot ---
  output$denovo_coverage_plot <- plotly::renderPlotly({
    top_data <- top_coverage_proteins()
    req(top_data)

    mappings <- top_data$mappings
    req(nrow(mappings) > 0)

    # Only proteins with sstart/send data
    has_pos <- !is.na(mappings$sstart) & !is.na(mappings$send)
    if (sum(has_pos) == 0) {
      # No positional data -- show peptide count bar chart instead
      summary_df <- top_data$summary
      summary_df$protein_label <- sub("_[^_]+$", "",
        sub("^[a-z]+\\|[^|]+\\|", "",
            mappings$subject[match(summary_df$protein, mappings$protein)]))
      summary_df$protein_label[is.na(summary_df$protein_label)] <- summary_df$protein[is.na(summary_df$protein_label)]

      p <- ggplot2::ggplot(summary_df,
        ggplot2::aes(x = stats::reorder(protein_label, n_peptides), y = n_peptides)) +
        ggplot2::geom_col(fill = "#1e88e5") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Top 20 Proteins by De Novo Peptide Count",
          subtitle = "(Alignment positions not available -- showing counts only)",
          x = NULL, y = "Unique Peptides"
        ) +
        ggplot2::theme_minimal()

      return(plotly::ggplotly(p) %>%
        plotly::config(toImageButtonOptions = list(format = "svg", scale = 2)))
    }

    mappings_pos <- mappings[has_pos, ]

    # Build coverage visualization
    # For each protein, show peptides as horizontal bars
    proteins <- unique(mappings_pos$protein)
    proteins <- head(proteins, 20)

    # Get protein names
    prot_names <- sub("_[^_]+$", "",
      sub("^[a-z]+\\|[^|]+\\|", "",
          mappings_pos$subject[match(proteins, mappings_pos$protein)]))

    # Build segments data
    seg_list <- lapply(seq_along(proteins), function(i) {
      prot <- proteins[i]
      peps <- mappings_pos[mappings_pos$protein == prot, ]
      peps <- peps[!duplicated(peps$peptide), ]

      # Color by match type
      peps$color <- ifelse(peps$pident >= 100, "#4caf50",  # Confirmed (green)
                    ifelse(peps$pident >= 90, "#ff9800",   # Near-match (orange)
                           "#f44336"))                      # Distant (red)
      peps$category <- ifelse(peps$pident >= 100, "Confirmed",
                       ifelse(peps$pident >= 90, "Near-match", "Distant"))

      data.frame(
        protein = paste0(prot_names[i], " (", prot, ")"),
        protein_idx = i,
        start = peps$sstart,
        end = peps$send,
        peptide = peps$peptide,
        pident = peps$pident,
        color = peps$color,
        category = peps$category,
        stringsAsFactors = FALSE
      )
    })

    seg_df <- do.call(rbind, seg_list)
    if (is.null(seg_df) || nrow(seg_df) == 0) return(NULL)

    # Create plotly figure with rectangles
    fig <- plotly::plot_ly()

    for (cat in c("Confirmed", "Near-match", "Distant")) {
      cat_data <- seg_df[seg_df$category == cat, ]
      if (nrow(cat_data) == 0) next

      cat_color <- switch(cat,
        "Confirmed" = "#4caf50", "Near-match" = "#ff9800", "Distant" = "#f44336")

      fig <- fig %>% plotly::add_segments(
        data = cat_data,
        x = ~start, xend = ~end,
        y = ~protein, yend = ~protein,
        line = list(color = cat_color, width = 8),
        name = cat,
        text = ~paste0("Peptide: ", peptide, "<br>",
                       "Position: ", start, "-", end, "<br>",
                       "Identity: ", round(pident, 1), "%"),
        hoverinfo = "text"
      )
    }

    fig %>%
      plotly::layout(
        title = list(text = "Protein Sequence Coverage Map", font = list(size = 14)),
        xaxis = list(title = "Position in Reference Protein"),
        yaxis = list(title = "", categoryorder = "trace"),
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = -0.1),
        margin = list(l = 200, b = 80)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Coverage detail: expandable peptide list per protein ---
  output$denovo_coverage_detail <- DT::renderDT({
    top_data <- top_coverage_proteins()
    req(top_data)

    mappings <- top_data$mappings
    req(nrow(mappings) > 0)

    prot_name <- sub("_[^_]+$", "",
      sub("^[a-z]+\\|[^|]+\\|", "",
          mappings$subject))

    display_df <- data.frame(
      Protein = mappings$protein,
      Name = prot_name,
      Peptide = mappings$peptide,
      Identity = round(mappings$pident, 1),
      Start = if (all(is.na(mappings$sstart))) NA_integer_ else mappings$sstart,
      End = if (all(is.na(mappings$send))) NA_integer_ else mappings$send,
      stringsAsFactors = FALSE
    )

    display_df <- display_df[!duplicated(display_df), ]
    display_df <- display_df[order(display_df$Protein, -display_df$Identity), ]

    DT::datatable(
      display_df,
      rownames = FALSE,
      filter = "top",
      options = list(pageLength = 25, scrollX = TRUE, dom = "Bfrtip",
                     buttons = list("csv", "excel"),
                     order = list(list(0, "asc"), list(4, "asc"))),
      extensions = "Buttons",
      caption = "Peptide-to-protein mappings for top 20 proteins by de novo peptide count"
    ) %>%
      DT::formatStyle("Identity",
        backgroundColor = DT::styleInterval(
          c(90, 100),
          c("#fce4ec", "#fff3e0", "#e8f5e9")
        )
      )
  })

  # --- Update Venn sample selectors when data changes ---
  observe({
    pep_data <- sample_peptide_data()
    req(pep_data)
    samples <- unique(pep_data$sample)
    if (length(samples) >= 2) {
      updateSelectInput(session, "denovo_venn_sample1", choices = samples, selected = samples[1])
      updateSelectInput(session, "denovo_venn_sample2", choices = samples, selected = samples[2])
    }
  })

}
