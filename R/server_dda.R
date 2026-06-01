# ==============================================================================
#  SERVER MODULE -- DDA Search (Sage pipeline on Hive)
#  Called from app.R as: server_dda(input, output, session, values, add_to_log)
# ==============================================================================

server_dda <- function(input, output, session, values, add_to_log) {

  # Load sage_bin path from config
  config <- tryCatch(yaml::read_yaml("config.yml"), error = function(e) list())
  sage_bin <- config$tools$sage_bin %||%
    "/quobyte/proteomics-grp/de-limp/cascadia/sage-v0.14.7-x86_64-unknown-linux-gnu/sage"
  slurm_account   <- config$slurm$account   %||% "genome-center-grp"
  slurm_partition <- config$slurm$partition  %||% "high"

  # ─── DDA + Casanovo job queue persistence ───────────────────────────────
  # Mirrors values$dda_jobs (list of submitted runs) to disk so the queue
  # survives Shiny restarts. Same per-user-local pattern as DIA-NN job_queue
  # and proteog_builds. Active job state (values$dda_job_id / dda_status /
  # dda_output_dir) continues to reflect the most-recently-active run.
  dda_queue_path <- function() file.path(path.expand("~"), ".delimp_dda_queue.rds")
  dda_queue_save <- function(jobs) {
    tryCatch(saveRDS(jobs, dda_queue_path()),
             error = function(e) message("[DDA] queue save failed: ",
                                          conditionMessage(e)))
  }
  dda_queue_load <- function() {
    p <- dda_queue_path()
    if (!file.exists(p)) return(list())
    tryCatch(readRDS(p), error = function(e) {
      message("[DDA] queue load failed: ", conditionMessage(e)); list()
    })
  }
  # Restore on module init — wrap in isolate() because reactiveValues access
  # outside a reactive consumer throws.
  isolate({
    if (length(values$dda_jobs %||% list()) == 0) {
      restored <- dda_queue_load()
      if (length(restored) > 0) {
        values$dda_jobs <- restored
        message(sprintf("[DDA] restored %d job(s) from disk", length(restored)))
      }
    }
  })
  # Persist on every change — never clobber a non-empty file with an empty
  # list (avoids the startup-race bug we hit on proteog).
  observe({
    jobs <- values$dda_jobs %||% list()
    if (length(jobs) == 0 && file.exists(dda_queue_path())) {
      existing <- tryCatch(readRDS(dda_queue_path()), error = function(e) NULL)
      if (length(existing) > 0) return()
    }
    dda_queue_save(jobs)
  })

  # Sync `values$dda_status` transitions into the matching queue entry so the
  # persisted queue reflects "running" → "done" / "error" without touching
  # every set-status call site.
  observe({
    st <- values$dda_status
    jid <- isolate(values$dda_job_id)
    jobs <- isolate(values$dda_jobs %||% list())
    if (is.null(st) || is.null(jid) || length(jobs) == 0) return()
    idx <- which(vapply(jobs, function(j) identical(as.character(j$job_id),
                                                     as.character(jid)),
                        logical(1)))
    if (length(idx) == 0) return()
    if (!identical(jobs[[idx[1]]]$status, st)) {
      jobs[[idx[1]]]$status <- st
      if (st %in% c("done", "error", "loaded"))
        jobs[[idx[1]]]$finished_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
      values$dda_jobs <- jobs
    }
  })

  # Casanovo config defaults
  casanovo_conda_env   <- config$tools$casanovo_conda_env %||%
    "/quobyte/proteomics-grp/conda_envs/cassonovo_env"
  casanovo_model_ckpt  <- config$tools$casanovo_model_ckpt %||%
    "/quobyte/proteomics-grp/bioinformatics_programs/casanovo_modles/casanovo_v4_2_0.ckpt"
  casanovo_converter   <- config$tools$casanovo_converter %||%
    "/quobyte/proteomics-grp/de-limp/python/bruker_to_mgf.py"
  casanovo_gpu_partition <- config$slurm$gpu_partition %||% "gpu-a100"
  casanovo_gpu_qos       <- config$slurm$gpu_qos %||% "genome-center-grp-gpu-a100-qos"

  # DIAMOND BLAST database paths (pre-built SwissProt/TrEMBL on shared storage)
  swissprot_dmnd <- config$blast$swissprot_dmnd %||%
    "/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/uniprot_sprot"
  trembl_dmnd    <- config$blast$trembl_dmnd %||%
    "/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/uniprot_trembl"

  # ---- DDA workflow shortcuts → jump to Run Search with mode pre-set ----
  # The DDA workflows dropdown has 3 friendly landing pages (De Novo,
  # Peptidomics, HLA). Each has a button that flips the user to the actual
  # Run Search tab with acquisition_mode=DDA and the right preset selected.
  dda_workflow_jump <- function(preset_value) {
    updateRadioButtons(session, "acquisition_mode", selected = "dda")
    updateSelectInput(session, "dda_preset", selected = preset_value)
    nav_select(id = "main_tabs", selected = "search_tab", session = session)
  }
  observeEvent(input$dda_workflow_open_denovo,      dda_workflow_jump("standard"))
  observeEvent(input$dda_workflow_open_peptidomics, dda_workflow_jump("peptidomics"))
  observeEvent(input$dda_workflow_open_hla, {
    # User picked Class I vs Class II in the radio on the HLA landing panel
    chosen <- input$dda_workflow_hla_class %||% "hla_class_i"
    dda_workflow_jump(chosen)
  })

  # ---- Preset / mode hint (visible under the preset selector) ----
  # Lets the user see at a glance what enzyme / length / mass / charge
  # window each mode uses, so picking HLA Class I doesn't feel magical.
  output$dda_preset_hint <- renderUI({
    p <- input$dda_preset %||% "standard"
    info <- switch(p,
      "standard"     = "Trypsin/P, 7–50 AA, ox(M) + acetyl-Nt. Default proteomics.",
      "phospho"      = "Trypsin/P, 7–50 AA, +Phospho(STY) variable mods. PhosphoSTY discovery.",
      "tmt"          = "Trypsin/P, 7–50 AA, TMT static (Nt + K). 11-plex iTRAQ-style runs.",
      "peptidomics"  = "Nonspecific, 5–25 AA, 400–5000 Da. Mods: ox(M), pyro-Glu, C-term amidation, Nt acetyl. Endogenous peptides.",
      "hla_class_i"  = "Nonspecific, 8–12 AA, 700–1500 Da. Mods: ox(M), deamid(N/Q). Charge 1–3 (z=1 dominant on TOF).",
      "hla_class_ii" = "Nonspecific, 13–25 AA, 1300–3000 Da. Mods: ox(M), deamid(N/Q).",
      "(unknown preset)"
    )
    tags$small(style = "color: #6c757d; display: block; margin-top: 4px; line-height: 1.3;",
               icon("circle-info"), " ", info)
  })

  # ---- Casanovo GPU vs CPU compute toggle --------------------------------
  # Live probe of gpu-a100 queue depth + node state. Drives the hint text
  # under the GPU/CPU radio so the user can pick the faster path.
  gpu_queue_state <- reactiveVal(NULL)
  gpu_queue_last_check <- reactiveVal(0)

  probe_gpu_queue <- function() {
    if (!isTRUE(values$ssh_connected)) return(invisible(NULL))
    ssh_cfg <- tryCatch(dda_ssh_config(), error = function(e) NULL)
    if (is.null(ssh_cfg)) return(invisible(NULL))
    # Show a brief "probing" placeholder so the user knows it's working
    gpu_queue_state(list(success = FALSE,
                         reason  = "Probing gpu-a100 queue...",
                         recommend = "gpu"))
    res <- tryCatch(
      check_casanovo_gpu_queue(
        ssh_cfg, "gpu-a100",
        sbatch_path = values$ssh_sbatch_path   # full sbatch path → skip login_shell (fast)
      ),
      error = function(e) {
        list(success = FALSE, reason = paste("Probe failed:", e$message),
             recommend = "gpu")
      }
    )
    gpu_queue_state(res)
    gpu_queue_last_check(as.numeric(Sys.time()))
    # Auto-flip the radio to whatever the probe recommends (user can override).
    if (isTRUE(res$success) && !is.null(res$recommend) &&
        !identical(input$dda_casanovo_compute, res$recommend)) {
      updateRadioButtons(session, "dda_casanovo_compute",
                        selected = res$recommend)
    }
  }

  # Probe when the user opens the Casanovo accordion section
  observeEvent(input$dda_run_casanovo, {
    if (isTRUE(input$dda_run_casanovo)) probe_gpu_queue()
  }, ignoreInit = TRUE)

  # Manual refresh button
  observeEvent(input$dda_refresh_gpu_queue, {
    probe_gpu_queue()
  })

  output$dda_gpu_queue_hint <- renderUI({
    s <- gpu_queue_state()
    if (is.null(s)) {
      return(tags$small(
        style = "color: #6c757d;",
        icon("circle-info"),
        " GPU queue not checked yet. Toggle Casanovo on, or click Refresh."
      ))
    }
    if (!isTRUE(s$success)) {
      return(tags$small(style = "color: #b16e1f;",
                        icon("triangle-exclamation"), " ", s$reason))
    }
    badge_color <- if (identical(s$recommend, "gpu")) "#28a745" else "#dc3545"
    rec_label   <- if (identical(s$recommend, "gpu"))
                     "GPU recommended" else "CPU recommended"
    last_chk <- gpu_queue_last_check()
    age_sec  <- if (last_chk > 0) round(as.numeric(Sys.time()) - last_chk) else NA
    tags$div(
      style = "font-size: 12px; margin: 4px 0;",
      tags$span(
        style = sprintf(
          "display: inline-block; padding: 2px 8px; border-radius: 10px; background: %s; color: white; font-weight: 600;",
          badge_color),
        rec_label
      ),
      tags$span(style = "margin-left: 8px; color: #495057;",
        sprintf("%d pending · %d running · %d/%d GPUs free on usable nodes",
                s$pending %||% 0L, s$running %||% 0L,
                s$gpus_free %||% 0L, s$gpus_total %||% 0L)
      ),
      tags$br(),
      tags$small(style = "color: #6c757d;",
        s$reason,
        if (!is.na(age_sec)) sprintf("  (refreshed %ds ago)", age_sec) else ""
      )
    )
  })

  # --- Mode observer: sync input to reactive values ---
  observeEvent(input$acquisition_mode, {
    values$acquisition_mode <- input$acquisition_mode
  })

  # --- Contextual label for the mode switcher ---
  output$mode_context_label <- renderUI({
    mode <- input$acquisition_mode %||% "dia"
    switch(mode,
      "dia" = tags$span(
        style = "font-size: 12px; color: #0d6efd; font-style: italic;",
        icon("circle-check", style = "color: #198754;"),
        " DIA-NN + limpa pipeline"
      ),
      "dda" = tags$span(
        style = "font-size: 12px; color: #6c757d; font-style: italic;",
        icon("circle-info", style = "color: #0d6efd;"),
        " Sage + Casanovo pipeline"
      ),
      "xlms" = tags$span(
        style = "font-size: 12px; color: #6c757d; font-style: italic;",
        icon("diagram-project", style = "color: #6f42c1;"),
        " MeroX + xiSearch + network"
      )
    )
  })

  # ============================================================================
  #    SSH config helper (reuses DIA-NN search SSH settings)
  # ============================================================================
  dda_ssh_config <- reactive({
    req(values$ssh_connected)
    list(
      host     = isolate(input$ssh_host),
      user     = isolate(input$ssh_user),
      port     = isolate(input$ssh_port) %||% 22,
      key_path = isolate(input$ssh_key_path),
      modules  = isolate(input$ssh_modules) %||% ""
    )
  })

  # ============================================================================
  #    FASTA Database — UniProt download, SSH browse, contaminant append
  # ============================================================================

  # Track the resolved DDA FASTA path (from any source)
  dda_fasta_resolved <- reactiveVal(NULL)

  # --- UniProt modal ---
  observeEvent(input$dda_open_uniprot_modal, {
    showModal(modalDialog(
      title = tagList(icon("dna"), " UniProt FASTA Database Search (DDA)"),
      size = "l",
      easyClose = TRUE,
      div(style = "display: flex; gap: 8px; margin-bottom: 12px;",
        div(style = "flex: 1;",
          textInput("dda_uniprot_search_query", NULL,
            placeholder = "e.g., human, mouse, E. coli", width = "100%")
        ),
        actionButton("dda_search_uniprot", "Search",
          class = "btn-info", style = "margin-top: 0;")
      ),
      DTOutput("dda_uniprot_results_table"),
      hr(),
      div(style = "display: flex; gap: 12px; align-items: flex-end;",
        div(style = "flex: 1;",
          selectInput("dda_fasta_content_type", "Content:",
            choices = c(
              "One per gene (recommended)" = "one_per_gene",
              "Swiss-Prot reviewed" = "reviewed",
              "Swiss-Prot + isoforms" = "reviewed_isoforms",
              "Full proteome" = "full",
              "Full + isoforms" = "full_isoforms"
            ), selected = "one_per_gene", width = "100%")
        ),
        div(style = "flex: 1;",
          uiOutput("dda_fasta_filename_preview_modal")
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("dda_download_fasta_btn", "Download FASTA",
          class = "btn-success", icon = icon("download"))
      )
    ))
  })

  # UniProt search
  observeEvent(input$dda_search_uniprot, {
    req(nzchar(input$dda_uniprot_search_query))
    withProgress(message = "Searching UniProt...", {
      results <- search_uniprot_proteomes(input$dda_uniprot_search_query)
      values$dda_uniprot_results <- results
    })
    if (nrow(values$dda_uniprot_results) == 0) {
      showNotification("No proteomes found. Try a different search term.", type = "warning")
    }
  })

  # UniProt results table
  output$dda_uniprot_results_table <- DT::renderDT({
    req(values$dda_uniprot_results, nrow(values$dda_uniprot_results) > 0)
    display_df <- values$dda_uniprot_results[, c("upid", "organism", "common_name", "protein_count")]
    colnames(display_df) <- c("ID", "Organism", "Common Name", "Proteins")
    DT::datatable(display_df,
      selection = "single",
      options = list(pageLength = 10, dom = "tip", scrollY = "300px",
        columnDefs = list(list(width = "90px", targets = 0))),
      rownames = FALSE, class = "compact stripe")
  })

  # Filename preview in modal
  output$dda_fasta_filename_preview_modal <- renderUI({
    req(values$dda_uniprot_results, nrow(values$dda_uniprot_results) > 0)
    sel <- input$dda_uniprot_results_table_rows_selected
    req(length(sel) > 0)
    row <- values$dda_uniprot_results[sel, ]
    fname <- generate_fasta_filename(row$upid, row$organism, input$dda_fasta_content_type)
    div(style = "font-size: 0.85em; color: #6c757d; padding-top: 28px;",
      icon("file"), " ", fname)
  })

  # Download FASTA from UniProt
  observeEvent(input$dda_download_fasta_btn, {
    req(values$dda_uniprot_results, nrow(values$dda_uniprot_results) > 0)
    sel <- input$dda_uniprot_results_table_rows_selected
    if (length(sel) == 0) {
      showNotification("Please select a proteome from the table first.", type = "warning")
      return()
    }

    row <- values$dda_uniprot_results[sel, ]
    fname <- generate_fasta_filename(row$upid, row$organism, input$dda_fasta_content_type)

    # Download locally first
    fasta_dir <- getOption("delimp.fasta_dir",
      default = "/quobyte/proteomics-grp/de-limp/fasta")
    if (!dir.exists(fasta_dir)) dir.create(fasta_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(fasta_dir)) fasta_dir <- tempdir()
    output_path <- file.path(fasta_dir, fname)

    withProgress(message = sprintf("Downloading %s from UniProt...", row$upid), {
      result <- download_uniprot_fasta(
        proteome_id  = row$upid,
        content_type = input$dda_fasta_content_type,
        output_path  = output_path
      )
    })

    if (!result$success) {
      showNotification(paste("Download failed:", result$error), type = "error")
      return()
    }
    if (!is.null(result$warning)) {
      showNotification(result$warning, type = "warning", duration = 12)
    }

    removeModal()

    # Upload to HPC if SSH connected
    ssh_cfg <- tryCatch(dda_ssh_config(), error = function(e) NULL)
    if (!is.null(ssh_cfg)) {
      remote_fasta_dir <- file.path(
        "/quobyte/proteomics-grp/de-limp", ssh_cfg$user, "databases")
      remote_path <- file.path(remote_fasta_dir, fname)

      # Check if already exists with same sequence count
      needs_upload <- TRUE
      exists_check <- ssh_exec(ssh_cfg,
        paste("test -f", shQuote(remote_path), "&& grep -c '^>'", shQuote(remote_path)))
      remote_count <- suppressWarnings(
        as.integer(trimws(paste(exists_check$stdout, collapse = ""))))
      if (!is.na(remote_count) && remote_count == result$n_sequences) {
        needs_upload <- FALSE
      }
      if (needs_upload) {
        ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(remote_fasta_dir)))
        withProgress(message = "Uploading FASTA to HPC...", {
          scp_upload(ssh_cfg, output_path, remote_path)
        })
      }
      dda_fasta_resolved(remote_path)
      showNotification(
        sprintf("FASTA ready on HPC: %s (%s sequences)",
          basename(remote_path), format(result$n_sequences, big.mark = ",")),
        type = "message", duration = 8)
    } else {
      # No SSH — use local path (won't work for HPC submit but shown for reference)
      dda_fasta_resolved(output_path)
      showNotification(
        sprintf("FASTA downloaded: %s (%s sequences). Connect SSH to upload to HPC.",
          basename(output_path), format(result$n_sequences, big.mark = ",")),
        type = "warning", duration = 10)
    }

    values$dda_fasta_info <- list(
      organism = row$common_name,
      n_sequences = result$n_sequences,
      filename = fname
    )
  })

  # Show selected FASTA info
  output$dda_fasta_selected_info <- renderUI({
    fpath <- dda_fasta_resolved()
    info  <- values$dda_fasta_info
    if (!is.null(fpath) && nzchar(fpath)) {
      div(style = "font-size: 0.85em; margin-top: 8px; padding: 8px; background: #e8f5e9; border-radius: 6px;",
        icon("check-circle", style = "color: #28a745;"), " ",
        tags$strong(basename(fpath)),
        if (!is.null(info$n_sequences))
          paste0(" (", format(info$n_sequences, big.mark = ","), " sequences)"),
        if (!is.null(info$organism))
          paste0(" -- ", info$organism)
      )
    }
  })

  # Keep dda_fasta_path synced when user types in the browse/path textInput
  observeEvent(input$dda_fasta_path, {
    p <- trimws(input$dda_fasta_path %||% "")
    if (nzchar(p)) dda_fasta_resolved(p)
  }, ignoreInit = TRUE)

  # Sync resolved FASTA path to values for cross-module access (DIAMOND BLAST)
  observe({
    p <- dda_fasta_resolved()
    if (!is.null(p) && nzchar(p)) values$dda_fasta_path <- p
  })

  # --- NCBI modal (reuse DIA pattern) ---
  observeEvent(input$dda_open_ncbi_modal, {
    showModal(modalDialog(
      title = "Download FASTA from NCBI",
      textInput("dda_ncbi_organism", "Organism name", placeholder = "e.g., Gallus gallus"),
      actionButton("dda_ncbi_search_btn", "Search NCBI", class = "btn-success btn-sm", icon = icon("search")),
      uiOutput("dda_ncbi_results_ui"),
      footer = modalButton("Close"),
      size = "l"
    ))
  })

  observeEvent(input$dda_ncbi_search_btn, {
    req(nzchar(input$dda_ncbi_organism))
    tryCatch({
      out_dir <- file.path(tempdir(), "ncbi_fasta")
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      fasta_path <- ncbi_download_proteome(input$dda_ncbi_organism, output_dir = out_dir)
      if (!is.null(fasta_path) && file.exists(fasta_path)) {
        n_seq <- length(grep("^>", readLines(fasta_path, warn = FALSE)))
        # Upload to HPC
        ssh_cfg <- dda_ssh_config()
        remote_dir <- "/quobyte/proteomics-grp/de-limp/fasta"
        ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(remote_dir)), timeout = 10)
        remote_path <- file.path(remote_dir, basename(fasta_path))
        scp_upload(ssh_cfg, fasta_path, remote_path)
        dda_fasta_resolved(remote_path)
        output$dda_ncbi_fasta_selected_info <- renderUI({
          div(class = "alert alert-success small py-1 mt-1",
            icon("check"), sprintf(" %s (%d sequences)", basename(remote_path), n_seq))
        })
        removeModal()
        showNotification(sprintf("NCBI FASTA uploaded: %s", basename(remote_path)), type = "message")
      } else {
        showNotification("NCBI download returned no FASTA. Try a different organism name or accession.", type = "warning")
      }
    }, error = function(e) {
      showNotification(paste("NCBI download failed:", e$message), type = "error")
    })
  })

  # --- Database Library ---
  output$dda_fasta_library_ui <- renderUI({
    req(values$ssh_connected)
    ssh_cfg <- dda_ssh_config()
    lib_path <- "/quobyte/proteomics-grp/dia-nn/fasta_library"
    result <- ssh_exec(ssh_cfg, paste("ls", shQuote(lib_path), "2>/dev/null"), timeout = 10)
    fastas <- grep("\\.(fasta|fa|faa)$", trimws(result$stdout), value = TRUE)
    if (length(fastas) == 0) {
      div(class = "text-muted small", "No FASTA files found in library")
    } else {
      tagList(
        selectInput("dda_library_fasta", "Select database",
          choices = setNames(file.path(lib_path, fastas), fastas)),
        actionButton("dda_select_library_fasta", "Use selected", class = "btn-sm btn-outline-primary")
      )
    }
  })

  observeEvent(input$dda_select_library_fasta, {
    req(input$dda_library_fasta)
    dda_fasta_resolved(input$dda_library_fasta)
    showNotification(sprintf("Using: %s", basename(input$dda_library_fasta)), type = "message")
  })

  # --- SSH file browser for FASTA ---
  observeEvent(input$dda_ssh_browse_fasta_btn, {
    req(values$ssh_connected)
    # Reuse the SSH file browser infrastructure from server_search.R
    # Open file browser modal with FASTA filter
    ssh_cfg <- dda_ssh_config()
    start_dir <- "/quobyte/proteomics-grp/dia-nn/fasta_library"

    # Check if the fasta library dir exists, fall back to user home
    dir_check <- ssh_exec(ssh_cfg,
      paste("test -d", shQuote(start_dir), "&& echo EXISTS"), timeout = 10)
    if (!any(grepl("EXISTS", dir_check$stdout))) {
      start_dir <- paste0("/quobyte/proteomics-grp/de-limp/", ssh_cfg$user)
    }

    # List .fasta files in the directory
    ls_result <- ssh_exec(ssh_cfg,
      paste0("ls -1 ", shQuote(start_dir), "/*.fasta ", shQuote(start_dir), "/*.fa 2>/dev/null | head -100"),
      timeout = 15)

    fasta_files <- character(0)
    if (ls_result$status == 0 && length(ls_result$stdout) > 0) {
      fasta_files <- trimws(ls_result$stdout)
      fasta_files <- fasta_files[nzchar(fasta_files)]
    }

    # Also list subdirectories
    ls_dirs <- ssh_exec(ssh_cfg,
      paste0("ls -1d ", shQuote(start_dir), "/*/ 2>/dev/null | head -50"),
      timeout = 15)
    subdirs <- character(0)
    if (ls_dirs$status == 0 && length(ls_dirs$stdout) > 0) {
      subdirs <- trimws(ls_dirs$stdout)
      subdirs <- subdirs[nzchar(subdirs)]
    }

    values$dda_fasta_browser_dir <- start_dir

    showModal(modalDialog(
      title = tagList(icon("folder-open"), " Browse FASTA Files on HPC"),
      size = "l", easyClose = TRUE,
      div(style = "margin-bottom: 12px;",
        div(style = "display: flex; gap: 8px; align-items: flex-end;",
          div(style = "flex: 1;",
            textInput("dda_fasta_browse_dir", "Directory:",
              value = start_dir, width = "100%")
          ),
          actionButton("dda_fasta_browse_go", "Go",
            class = "btn-outline-primary btn-sm", style = "margin-bottom: 15px;")
        )
      ),
      uiOutput("dda_fasta_browser_content"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("dda_fasta_browser_select", "Select",
          class = "btn-success", icon = icon("check"))
      )
    ))
  })

  # Navigate within FASTA browser
  observeEvent(input$dda_fasta_browse_go, {
    req(values$ssh_connected)
    browse_dir <- trimws(input$dda_fasta_browse_dir %||% "")
    req(nzchar(browse_dir))
    values$dda_fasta_browser_dir <- browse_dir
  })

  # Click on a directory in the browser
  observeEvent(input$dda_fasta_browser_click_dir, {
    req(nzchar(input$dda_fasta_browser_click_dir))
    values$dda_fasta_browser_dir <- input$dda_fasta_browser_click_dir
    updateTextInput(session, "dda_fasta_browse_dir", value = input$dda_fasta_browser_click_dir)
  })

  # Click on a file in the browser to select it
  observeEvent(input$dda_fasta_browser_click_file, {
    req(nzchar(input$dda_fasta_browser_click_file))
    values$dda_fasta_browser_selected <- input$dda_fasta_browser_click_file
  })

  # Render the file browser content
  output$dda_fasta_browser_content <- renderUI({
    browse_dir <- values$dda_fasta_browser_dir
    req(nzchar(browse_dir))
    ssh_cfg <- dda_ssh_config()

    # List directory contents
    ls_result <- ssh_exec(ssh_cfg,
      paste0("ls -1ap ", shQuote(browse_dir), " 2>/dev/null | head -200"),
      timeout = 15)

    if (ls_result$status != 0) {
      return(div(class = "alert alert-warning", "Could not read directory: ", browse_dir))
    }

    items <- trimws(ls_result$stdout)
    items <- items[nzchar(items) & items != "./" & items != "../"]

    # Separate dirs and files
    dirs  <- items[grepl("/$", items)]
    files <- items[!grepl("/$", items)]
    # Filter to FASTA files only
    files <- files[grepl("\\.(fasta|fa|faa)$", files, ignore.case = TRUE)]

    selected <- values$dda_fasta_browser_selected

    dir_items <- lapply(dirs, function(d) {
      full_path <- file.path(browse_dir, sub("/$", "", d))
      tags$div(
        style = "padding: 4px 8px; cursor: pointer; border-bottom: 1px solid #eee;",
        onclick = sprintf("Shiny.setInputValue('dda_fasta_browser_click_dir', '%s', {priority: 'event'})", full_path),
        icon("folder", style = "color: #0d6efd; margin-right: 8px;"),
        tags$span(d, style = "font-weight: 500;")
      )
    })

    file_items <- lapply(files, function(f) {
      full_path <- file.path(browse_dir, f)
      is_selected <- identical(full_path, selected)
      bg <- if (is_selected) "background: #d4edda;" else ""
      tags$div(
        style = paste0("padding: 4px 8px; cursor: pointer; border-bottom: 1px solid #eee; ", bg),
        onclick = sprintf("Shiny.setInputValue('dda_fasta_browser_click_file', '%s', {priority: 'event'})", full_path),
        icon("file-alt", style = "color: #28a745; margin-right: 8px;"),
        tags$span(f)
      )
    })

    # Parent directory link
    parent <- dirname(browse_dir)
    parent_link <- if (parent != browse_dir) {
      tags$div(
        style = "padding: 4px 8px; cursor: pointer; border-bottom: 1px solid #eee;",
        onclick = sprintf("Shiny.setInputValue('dda_fasta_browser_click_dir', '%s', {priority: 'event'})", parent),
        icon("level-up-alt", style = "color: #6c757d; margin-right: 8px;"),
        tags$span(".. (parent directory)", style = "color: #6c757d;")
      )
    }

    div(style = "max-height: 400px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 6px;",
      parent_link,
      dir_items,
      if (length(file_items) == 0 && length(dir_items) == 0)
        div(style = "padding: 16px; color: #6c757d; text-align: center;",
          "No FASTA files found in this directory.")
      else
        file_items
    )
  })

  # Select button in browser modal
  observeEvent(input$dda_fasta_browser_select, {
    selected <- values$dda_fasta_browser_selected
    if (is.null(selected) || !nzchar(selected)) {
      showNotification("Click a FASTA file to select it first.", type = "warning")
      return()
    }
    dda_fasta_resolved(selected)
    updateTextInput(session, "dda_fasta_path", value = selected)
    values$dda_fasta_info <- list(filename = basename(selected))
    removeModal()
    showNotification(paste("Selected:", basename(selected)), type = "message")
  })

  # ============================================================================
  #    File scan: list .d files in remote directory
  # ============================================================================
  observeEvent(input$dda_scan_files, {
    req(values$ssh_connected, input$dda_raw_dir)
    raw_dir <- trimws(input$dda_raw_dir)
    if (!nzchar(raw_dir)) {
      showNotification("Please enter a raw file directory path.", type = "warning")
      return()
    }

    ssh_cfg <- dda_ssh_config()
    result <- ssh_exec(ssh_cfg,
      paste0("{ ls -1d ", shQuote(raw_dir), "/*.d 2>/dev/null; ls -1 ", shQuote(raw_dir), "/*.raw 2>/dev/null; } | head -200"),
      timeout = 15)

    if (result$status != 0 || length(result$stdout) == 0 ||
        all(!nzchar(trimws(result$stdout)))) {
      showNotification("No .d or .raw files found in the specified path.", type = "warning")
      values$dda_raw_files <- character(0)
      return()
    }

    files <- trimws(result$stdout)
    files <- files[nzchar(files)]
    values$dda_raw_files <- files
    showNotification(paste("Found", length(files), ".d files"), type = "message")
  })

  # File list preview
  output$dda_file_list_preview <- renderUI({
    files <- values$dda_raw_files
    if (is.null(files) || length(files) == 0) {
      return(tags$p(style = "color: #6c757d; font-style: italic;",
        "No files scanned yet. Enter a directory path and click Scan."))
    }
    basenames <- basename(files)
    tags$div(
      style = "max-height: 200px; overflow-y: auto; background: #f8f9fa; padding: 8px; border-radius: 6px; font-size: 12px;",
      tags$strong(paste(length(basenames), "files:")),
      tags$ul(style = "margin: 4px 0; padding-left: 16px;",
        lapply(basenames, function(f) tags$li(f))
      )
    )
  })

  # ============================================================================
  #    Load existing DDA results from HPC
  # ============================================================================

  # Load Results modal — triggered from DDA search panel OR De Novo tab top button
  load_results_modal <- function() {
    is_hf <- nzchar(Sys.getenv("SPACE_ID", ""))
    # On HF: ZIP-only (no SSH possible). Elsewhere: both options.
    hpc_panel <- if (!is_hf) {
      tabPanel("From HPC (SSH)",
        tags$br(),
        tags$p(style = "color: #6c757d; font-size: 0.9em;",
          "Requires an SSH connection to Hive. Downloads results from the given output directory."),
        textInput("dda_load_path", "Output directory on HPC",
          placeholder = "/quobyte/proteomics-grp/de-limp/brettsp/dda_output/dda_search",
          width = "100%")
      )
    } else NULL
    showModal(modalDialog(
      title = if (is_hf) "Load DDA / De Novo Results (Upload ZIP)" else "Load DDA / De Novo Results",
      size = "l",
      tabsetPanel(
        hpc_panel,
        tabPanel("Upload ZIP",
          tags$br(),
          tags$p(style = "color: #6c757d; font-size: 0.9em;",
            "Upload a ZIP containing your DDA / de novo results. Useful on Hugging Face ",
            "or when you've already pulled the files locally. Any subset of the files below works:"),
          tags$pre(style = "font-size: 0.82em; background: #f8f9fa; padding: 10px; border-radius: 4px;",
"results.zip
├── results.sage.tsv          (Sage PSMs — required for DB-search results)
├── lfq.tsv                    (optional: Sage label-free quant matrix)
├── report.parquet             (alternative to Sage — DIA-NN DDA)
└── casanovo/
    └── mztab/
        └── *.mztab            (Casanovo de novo results)"
          ),
          fileInput("dda_load_zip", NULL,
            accept = c(".zip", "application/zip"), width = "100%")
        )
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("dda_load_confirm", "Load", class = "btn-primary", icon = icon("download"))
      )
    ))
  }

  observeEvent(input$load_dda_results, load_results_modal())
  observeEvent(input$load_dda_results_top, load_results_modal())
  observeEvent(input$load_dda_results_top2, load_results_modal())

  # ── Info modal: explains both SSH-load and ZIP-upload formats ────────────
  # Aimed at users who got a results ZIP from a colleague and need to load it
  # into DE-LIMP on Hugging Face (no HPC access). Also useful for HPC users.
  observeEvent(input$load_dda_results_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " Loading DDA / De Novo Results"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.93em; line-height: 1.6;",

        tags$h6("What this loads"),
        tags$p("This panel ingests results from a DDA proteomics + de novo sequencing run that ",
               "happened ", tags$em("elsewhere"), " (e.g. on the UC Davis Hive HPC, where Sage / DIA-NN / ",
               "Casanovo ran as SLURM jobs). It does not run any search itself — it loads the outputs."),

        tags$h6("Two ways to load"),
        tags$ol(
          tags$li(tags$strong("From HPC (SSH): "),
                  "Paste the absolute path of the run's output directory on Hive ",
                  "(e.g. ", tags$code("/quobyte/proteomics-grp/de-limp/brettsp/dda_output/dda_search"),
                  "). DE-LIMP scp's the result files down. Requires the SSH connection in the sidebar."),
          tags$li(tags$strong("Upload ZIP: "),
                  "Drag-and-drop a ZIP file containing the results. ", tags$strong("This is how Hugging Face users load shared results"),
                  " — DE-LIMP unpacks the ZIP locally and parses the contents. No HPC connection needed.")
        ),

        tags$h6("ZIP layout"),
        tags$p("Include any subset of these files. The flat layout below is what HPC users see after a run ",
               "— if you have nested folders, DE-LIMP will find files by basename anywhere in the ZIP."),
        tags$pre(style = "font-size: 0.85em; background: #f8f9fa; padding: 12px; border-radius: 6px; border: 1px solid #e0e0e0;",
"results.zip
├── results.sage.tsv          ← Sage PSMs (DB-search). Required if you want to see Confirmed peptides.
├── lfq.tsv                    ← Optional: Sage label-free quant matrix → enables Quantification panel.
├── report.parquet             ← Alternative DB-search engine: DIA-NN DDA. Use instead of results.sage.tsv.
└── casanovo/
    └── mztab/
        └── *.mztab            ← Optional: Casanovo de novo PSMs → enables Novel Peptides panel + BLAST."
        ),

        tags$h6("What each file unlocks"),
        tags$ul(
          tags$li(tags$strong("Sage TSV (or DIA-NN parquet): "),
                  "Confirmed peptides table, Score Distribution, FDR analysis"),
          tags$li(tags$strong("lfq.tsv: "), "Per-sample quantification panel"),
          tags$li(tags$strong("Casanovo mztab: "),
                  "Novel peptides (no DB match) → classification → de novo confidence plot"),
          tags$li(tags$strong("blast_results.tsv at "), tags$code("denovo/blast_results.tsv"), ": ",
                  "Pre-computed DIAMOND BLAST hits — skips re-running BLAST on upload")
        ),

        div(class = "alert alert-info py-2 px-3 mt-2",
            style = "font-size: 0.88em;",
            icon("info-circle"),
            tags$strong(" For collaborators on HF: "),
            "you only need the ZIP. Brett (or whoever shared the link) will pre-package one ",
            "for you — drag it into the Upload ZIP tab and click Load. ",
            "Everything works without an HPC account."),

        div(class = "alert alert-warning py-2 px-3 mt-2",
            style = "font-size: 0.88em;",
            icon("exclamation-triangle"),
            tags$strong(" Note (ZIP mode): "),
            "DE-LIMP cannot ", tags$em("submit"), " new BLAST jobs when loading from a ZIP — ",
            "include ", tags$code("denovo/blast_results.tsv"), " in the ZIP if you want BLAST results.")
      )
    ))
  })

  # Flag for conditional panel — hide Load button when data exists
  output$denovo_has_data <- reactive({
    !is.null(values$denovo_classification) || !is.null(values$dda_casanovo_classification)
  })
  outputOptions(output, "denovo_has_data", suspendWhenHidden = FALSE)

  # ─── Per-file filter ──────────────────────────────────────────────────
  # Every downstream viz (length hist, anchor logos, cleavage motifs,
  # confirmed/novel tables, etc.) reads from `dda_filtered_psms()` instead
  # of `values$dda_sage_psms` directly. Empty selection = show all files.
  observe({
    req(values$dda_sage_psms)
    files <- sort(unique(as.character(values$dda_sage_psms$filename %||% character(0))))
    if (length(files) == 0) return()
    # Strip path + extension for display
    short <- basename(files)
    short <- sub("\\.(mzML|raw|d)$", "", short, ignore.case = TRUE)
    names(files) <- short
    updateSelectizeInput(session, "dda_file_filter",
      choices = files, selected = character(0), server = TRUE)
  })

  observeEvent(input$dda_file_filter_all, {
    updateSelectizeInput(session, "dda_file_filter", selected = character(0))
  })

  # TRUE for each PSM whose ENTIRE protein group is contaminant — i.e. every
  # accession carries the `Cont_` tag stamped by the contaminant FASTA we
  # appended at search time. Peptides shared with a real protein are kept.
  .dda_is_contaminant <- function(proteins) {
    vapply(strsplit(as.character(proteins), ";", fixed = TRUE), function(accs) {
      accs <- accs[nzchar(accs)]
      length(accs) > 0 && all(grepl("Cont_", accs, fixed = TRUE))
    }, logical(1))
  }

  dda_filtered_psms <- reactive({
    req(values$dda_sage_psms)
    dt <- values$dda_sage_psms
    # File filter (empty selection = all files)
    sel <- input$dda_file_filter
    if (!is.null(sel) && length(sel) > 0) {
      dt <- dt[as.character(dt$filename) %in% as.character(sel), ]
    }
    # Contaminant filter — drop Cont_-tagged groups (from the searched contaminant DB)
    if (isTRUE(input$dda_results_exclude_contaminants %||% TRUE) &&
        "proteins" %in% names(dt)) {
      dt <- dt[!.dda_is_contaminant(dt$proteins), ]
    }
    dt
  })

  output$dda_file_filter_summary <- renderUI({
    req(values$dda_sage_psms)
    n_all_psms   <- nrow(values$dda_sage_psms)
    n_filt_psms  <- nrow(dda_filtered_psms())
    n_all_files  <- length(unique(values$dda_sage_psms$filename))
    n_sel_files  <- length(input$dda_file_filter %||% character(0))
    contam_note  <- ""
    if (isTRUE(input$dda_results_exclude_contaminants %||% TRUE) &&
        "proteins" %in% names(values$dda_sage_psms)) {
      n_contam <- sum(.dda_is_contaminant(values$dda_sage_psms$proteins))
      if (n_contam > 0)
        contam_note <- sprintf(" — %s contaminant PSMs hidden",
                               format(n_contam, big.mark = ","))
    }
    if (n_sel_files == 0 || n_sel_files == n_all_files) {
      tags$small(style = "color: #6c757d; display: block; margin-top: 6px;",
        sprintf("Combined view — %s PSMs across %d files%s.",
                format(n_filt_psms, big.mark = ","), n_all_files, contam_note))
    } else {
      tags$small(style = "color: #1565c0; display: block; margin-top: 6px;",
        sprintf("Filtered — %s PSMs from %d of %d files%s.",
                format(n_filt_psms, big.mark = ","), n_sel_files, n_all_files, contam_note))
    }
  })

  # ─── Current loaded mode (drives mode-specific viz + ZIP export) ──────
  # Looks up the mode tag from the queue row whose output_dir matches the
  # currently loaded search. Falls back to "standard" if no match found.
  current_dda_mode <- reactive({
    od <- values$dda_loaded$output_dir %||% values$dda_output_dir
    if (is.null(od)) return("standard")
    jobs <- values$dda_jobs %||% list()
    for (j in jobs) {
      if (identical(j$output_dir, od)) return(j$mode %||% "standard")
    }
    "standard"
  })

  output$dda_mode_is_hla <- reactive({
    current_dda_mode() %in% c("hla_class_i", "hla_class_ii")
  })
  outputOptions(output, "dda_mode_is_hla", suspendWhenHidden = FALSE)

  output$dda_mode_is_peptidomics <- reactive({
    identical(current_dda_mode(), "peptidomics")
  })
  outputOptions(output, "dda_mode_is_peptidomics", suspendWhenHidden = FALSE)

  # TRUE once Casanovo de novo results are loaded — used to hide HLA/peptidomics
  # anchor panels, which are meaningless for a de novo run.
  output$dda_is_denovo <- reactive({
    !is.null(values$dda_casanovo_classification) || !is.null(values$dda_casanovo_psms)
  })
  outputOptions(output, "dda_is_denovo", suspendWhenHidden = FALSE)

  # ─── Length distribution (universal for every DDA mode) ─────────────────
  .dda_clean_peptide <- function(x) gsub("\\[.*?\\]|[^A-Z]", "", x)

  output$dda_length_hist <- plotly::renderPlotly({
    psms <- dda_filtered_psms()
    req(psms, nrow(psms) > 0)
    compare <- input$dda_compare_mode %||% "combined"
    if (compare == "per_file") {
      # Multi-trace overlay: one line per file
      pep <- .dda_clean_peptide(psms$peptide)
      df <- data.frame(
        filename = basename(as.character(psms$filename)),
        length   = nchar(pep),
        stringsAsFactors = FALSE
      )
      df <- df[df$length > 0 & df$length < 60, ]
      if (nrow(df) == 0) return(plotly::plot_ly() |> plotly::layout(title = "No peptides"))
      df$filename <- sub("\\.(mzML|raw|d)$", "", df$filename, ignore.case = TRUE)
      counts <- as.data.frame(table(filename = df$filename, length = df$length),
                              stringsAsFactors = FALSE)
      counts$length <- as.integer(counts$length)
      counts <- counts[counts$Freq > 0, ]
      plotly::plot_ly(counts, x = ~length, y = ~Freq, color = ~filename,
                      type = "scatter", mode = "lines+markers",
                      hovertemplate = "%{fullData.name}<br>len=%{x}<br>n=%{y}<extra></extra>") |>
        plotly::layout(
          xaxis = list(title = "Peptide length (AA)", dtick = 1),
          yaxis = list(title = "N PSMs"),
          legend = list(orientation = "h", y = -0.2),
          margin = list(t = 20, b = 80, l = 60, r = 20)
        ) |>
        plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
    } else {
      # Combined: single bar series with mode-signature highlighting
      pep <- .dda_clean_peptide(psms$peptide)
      lens <- nchar(pep)
      lens <- lens[lens > 0 & lens < 60]
      if (length(lens) == 0) return(plotly::plot_ly() |> plotly::layout(title = "No peptides"))
      df <- as.data.frame(table(length = lens), stringsAsFactors = FALSE)
      df$length <- as.integer(df$length)
      mode <- current_dda_mode()
      df$is_signature <- switch(mode,
        "hla_class_i"  = df$length %in% 8:11,
        "hla_class_ii" = df$length %in% 13:18,
        "peptidomics"  = df$length %in% 5:25,
        rep(FALSE, nrow(df))
      )
      bar_color <- ifelse(df$is_signature, "#1565c0", "#9aa6b2")
      plotly::plot_ly(df, x = ~length, y = ~Freq, type = "bar",
                      marker = list(color = bar_color),
                      hovertemplate = "len=%{x}<br>n=%{y}<extra></extra>") |>
        plotly::layout(
          xaxis = list(title = "Peptide length (AA)", dtick = 1),
          yaxis = list(title = "N PSMs"),
          margin = list(t = 20, b = 50, l = 60, r = 20),
          showlegend = FALSE
        ) |>
        plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
    }
  })

  # Helper: bar-frequency renderer shared between HLA anchor + peptidomics cleavage.
  # In per-file compare mode, each file gets its own facet (subplot row).
  .render_pos_freq <- function(psms, positions, names_lbl, colors,
                               min_len = 8, compare_mode = "combined") {
    pep <- .dda_clean_peptide(psms$peptide)
    keep <- nchar(pep) >= min_len
    pep <- pep[keep]
    fn  <- if (compare_mode == "per_file") basename(as.character(psms$filename[keep])) else NULL
    if (!is.null(fn)) fn <- sub("\\.(mzML|raw|d)$", "", fn, ignore.case = TRUE)
    if (length(pep) == 0) return(plotly::plot_ly() |> plotly::layout(title = "No peptides"))
    extract <- function(i) {
      if (i > 0) substr(pep, i, i) else substr(pep, nchar(pep), nchar(pep))
    }
    p_aas <- lapply(positions, extract)
    all_aa <- sort(unique(unlist(p_aas)))
    if (compare_mode == "per_file" && !is.null(fn)) {
      # Long format: file × position × residue → freq
      rows <- do.call(rbind, lapply(seq_along(positions), function(j) {
        data.frame(
          file = fn, position = names_lbl[j], residue = p_aas[[j]],
          stringsAsFactors = FALSE
        )
      }))
      counts <- aggregate(rep(1L, nrow(rows)), list(file = rows$file,
                          position = rows$position, residue = rows$residue),
                          sum)
      names(counts)[4] <- "n"
      totals <- aggregate(counts$n,
                          list(file = counts$file, position = counts$position), sum)
      names(totals)[3] <- "tot"
      counts <- merge(counts, totals, by = c("file", "position"))
      counts$pct <- 100 * counts$n / counts$tot
      counts$facet <- paste(counts$file, counts$position, sep = " · ")
      plotly::plot_ly(counts, x = ~residue, y = ~pct, color = ~position,
                      colors = colors, type = "bar",
                      hovertemplate = "%{x}: %{y:.1f}%<extra></extra>") |>
        plotly::layout(
          barmode = "group", xaxis = list(title = "Residue"),
          yaxis = list(title = "Frequency (%)"),
          legend = list(orientation = "h", y = 1.1),
          margin = list(t = 20, b = 50, l = 60, r = 20)
        ) |>
        plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
    } else {
      plt <- plotly::plot_ly()
      for (j in seq_along(positions)) {
        pcts <- as.integer(table(factor(p_aas[[j]], levels = all_aa)))
        pcts <- 100 * pcts / sum(pcts)
        plt <- plotly::add_bars(plt, x = all_aa, y = pcts,
                                name = names_lbl[j],
                                marker = list(color = colors[j]),
                                hovertemplate = sprintf("%s %%{x}: %%{y:.1f}%%<extra></extra>",
                                                        names_lbl[j]))
      }
      plotly::layout(plt,
        barmode = "group",
        xaxis = list(title = "Residue", categoryorder = "category ascending"),
        yaxis = list(title = "Frequency (%)"),
        margin = list(t = 20, b = 50, l = 60, r = 20),
        legend = list(orientation = "h", y = 1.1)
      ) |>
        plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
    }
  }

  # ─── HLA anchor residues (P2 + PΩ) ──────────────────────────────────────
  output$dda_anchor_freq <- plotly::renderPlotly({
    psms <- dda_filtered_psms(); req(psms, nrow(psms) > 0)
    .render_pos_freq(psms, positions = c(2, -1),
                     names_lbl = c("P2", "PΩ"),
                     colors = c("#1565c0", "#b16e1f"),
                     min_len = 8,
                     compare_mode = input$dda_compare_mode %||% "combined")
  })

  # ─── Peptidomics N-term + C-term flanks ─────────────────────────────────
  output$dda_cleavage_freq <- plotly::renderPlotly({
    psms <- dda_filtered_psms(); req(psms, nrow(psms) > 0)
    .render_pos_freq(psms, positions = c(1, -1),
                     names_lbl = c("N-term", "C-term"),
                     colors = c("#198754", "#dc3545"),
                     min_len = 5,
                     compare_mode = input$dda_compare_mode %||% "combined")
  })

  # ─── Peptide × File matrix (presence + PSM counts) ──────────────────────
  # Rows = unique stripped peptides; columns = files; cells = PSM count.
  # Sortable, downloadable. Lets you spot peptides unique to a condition.
  output$dda_peptide_file_matrix <- DT::renderDT({
    psms <- dda_filtered_psms(); req(psms, nrow(psms) > 0)
    pep_stripped <- .dda_clean_peptide(psms$peptide)
    fn <- sub("\\.(mzML|raw|d)$", "", basename(as.character(psms$filename)),
              ignore.case = TRUE)
    tbl <- table(peptide = pep_stripped, file = fn)
    mat <- as.matrix(unclass(tbl))
    df <- data.frame(
      peptide  = rownames(mat),
      length   = nchar(rownames(mat)),
      n_files  = rowSums(mat > 0),
      total    = rowSums(mat),
      check.names = FALSE
    )
    df <- cbind(df, as.data.frame(mat, check.names = FALSE))
    df <- df[order(-df$total), ]   # most-abundant peptides first
    DT::datatable(df,
      rownames = FALSE, filter = "top",
      extensions = c("Buttons", "FixedColumns", "Scroller"),
      options = list(
        scrollX = TRUE, scrollY = "55vh", scroller = TRUE,
        deferRender = TRUE,
        fixedColumns = list(leftColumns = 4),
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel"),
        pageLength = 200
      )) |>
      DT::formatStyle(columns = colnames(mat),
                      background = DT::styleInterval(c(0),
                        c("transparent", "#e3f2fd")))
  })

  # ─── HF-viewable ZIP export ─────────────────────────────────────────────
  # Bundles the current loaded DDA search into a flat .zip that renders well
  # on a HuggingFace Space: methods.md, settings.json, sage_results.tsv,
  # casanovo/*.mztab, diamond_hits.tsv, mode-specific summary CSVs, PROMPT.md,
  # MANIFEST.txt. Mode comes from values$dda_loaded$mode or the queue row
  # for the currently loaded output_dir.
  output$dda_export_zip <- downloadHandler(
    filename = function() {
      mode <- isolate(values$dda_loaded$mode %||% "standard")
      sprintf("delimp_dda_%s_%s.zip", gsub("[^a-zA-Z0-9]", "_", mode),
              format(Sys.time(), "%Y%m%d_%H%M%S"))
    },
    content = function(file) {
      output_dir <- isolate(values$dda_loaded$output_dir %||% values$dda_output_dir)
      mode       <- isolate(values$dda_loaded$mode %||% "standard")
      app_ver    <- isolate(values$app_version %||% "unknown")
      if (is.null(output_dir) || !dir.exists(output_dir)) {
        # In SSH/HPC mode the canonical artifacts live on Hive — they get
        # downloaded into a local tempdir during the Load step. Refuse rather
        # than build a half-empty bundle.
        showNotification(
          "Export requires a local output directory. Load the search first.",
          type = "error", duration = 10)
        stop("dda_export_zip: no local output_dir available")
      }
      zip_path <- generate_dda_export_zip(
        output_dir  = output_dir,
        mode        = mode,
        app_version = app_ver
      )
      file.copy(zip_path, file, overwrite = TRUE)
    }
  )

  # Pre-flight check — verify the requested output dir actually has the
  # expected files before kicking off scp_download + parse. Refuses to load
  # search dirs where the underlying SLURM jobs failed (so users don't load
  # zero/garbage results).
  dda_preflight_check <- function(remote_dir, ssh_cfg) {
    if (is.null(ssh_cfg) || !nzchar(remote_dir %||% "")) {
      return(list(ok = TRUE, msg = ""))  # ZIP mode skips this
    }
    # Use one ssh_exec with stat checks for all files we care about
    cmd <- sprintf(paste(
      "echo SAGE=$(test -s %s/results.sage.tsv && stat -c%%s %s/results.sage.tsv || echo 0);",
      "echo DIANN=$(test -s %s/report.parquet && stat -c%%s %s/report.parquet || echo 0);",
      "echo MZTAB=$(ls %s/casanovo/mztab/*.mztab 2>/dev/null | wc -l);",
      "echo BLAST=$(test -s %s/denovo/blast_results.tsv && echo 1 || echo 0)"),
      shQuote(remote_dir), shQuote(remote_dir),
      shQuote(remote_dir), shQuote(remote_dir),
      shQuote(remote_dir), shQuote(remote_dir))
    res <- tryCatch(ssh_exec(ssh_cfg, cmd, login_shell = FALSE, timeout = 15),
                    error = function(e) NULL)
    if (is.null(res) || !identical(res$status, 0L)) {
      return(list(ok = FALSE, msg = "Could not reach the output directory on Hive."))
    }
    txt <- paste(res$stdout %||% character(), collapse = "\n")
    sage_size  <- as.numeric(sub(".*SAGE=([0-9]+).*",  "\\1", txt))
    diann_size <- as.numeric(sub(".*DIANN=([0-9]+).*", "\\1", txt))
    n_mztab    <- as.numeric(sub(".*MZTAB=([0-9]+).*", "\\1", txt))
    has_blast  <- grepl("BLAST=1", txt)
    has_sage_or_diann <- (!is.na(sage_size) && sage_size > 100) ||
                         (!is.na(diann_size) && diann_size > 100)
    has_mztab <- !is.na(n_mztab) && n_mztab > 0
    if (!has_sage_or_diann && !has_mztab) {
      return(list(ok = FALSE, msg = paste0(
        "No usable results found in ", remote_dir, ".\n",
        "Expected: results.sage.tsv (Sage), report.parquet (DIA-NN), or casanovo/mztab/*.mztab.\n",
        "All searches in this directory may have failed — check sacct for the job IDs."
      )))
    }
    # Build a status summary so the user knows what's about to load
    parts <- c()
    if (!is.na(sage_size) && sage_size > 100)
      parts <- c(parts, sprintf("Sage (%s)", format(sage_size, big.mark = ",")))
    if (!is.na(diann_size) && diann_size > 100)
      parts <- c(parts, sprintf("DIA-NN (%s)", format(diann_size, big.mark = ",")))
    if (has_mztab)
      parts <- c(parts, sprintf("Casanovo (%d .mztab)", n_mztab))
    if (has_blast)
      parts <- c(parts, "BLAST")
    list(ok = TRUE, msg = paste(parts, collapse = " + "))
  }

  observeEvent(input$dda_load_confirm, {
    # Either a ZIP upload (HF/local mode) or an HPC path must be provided
    zip_uploaded <- !is.null(input$dda_load_zip) &&
                     is.data.frame(input$dda_load_zip) &&
                     nrow(input$dda_load_zip) > 0 &&
                     file.exists(input$dda_load_zip$datapath[1])

    if (!zip_uploaded && !nzchar(input$dda_load_path %||% "")) {
      showNotification(
        "Provide a ZIP file or an HPC output directory.",
        type = "warning", duration = 6)
      return()
    }

    # Pre-flight: refuse to load search dirs whose jobs all failed.
    if (!zip_uploaded) {
      pre <- dda_preflight_check(trimws(input$dda_load_path),
                                  dda_ssh_config())
      if (!isTRUE(pre$ok)) {
        showNotification(
          tags$div(tags$strong("Pre-flight check failed."),
                   tags$br(), tags$pre(style="margin-top:6px; font-size:0.85em;", pre$msg)),
          type = "error", duration = 15)
        return()
      } else if (nzchar(pre$msg)) {
        showNotification(sprintf("Loading: %s", pre$msg),
                         type = "default", duration = 5)
      }
    }

    # Set up local working dir + populate it (from ZIP or SSH)
    local_tmp <- file.path(tempdir(), sprintf("dda_load_%s",
                                              format(Sys.time(), "%Y%m%d_%H%M%S")))
    if (dir.exists(local_tmp)) unlink(local_tmp, recursive = TRUE)
    dir.create(local_tmp, showWarnings = FALSE, recursive = TRUE)

    ssh_cfg <- NULL
    remote_dir <- ""

    if (zip_uploaded) {
      # ── ZIP-upload path: unpack into local_tmp, then flatten layout so the
      # downstream code (which expects results.sage.tsv at local_tmp/results.sage.tsv,
      # .mztab at local_tmp/mztab/, etc.) finds files at canonical paths regardless
      # of how the user structured their ZIP.
      withProgress(message = "Unpacking ZIP...", value = 0.05, {
        tryCatch({
          utils::unzip(input$dda_load_zip$datapath[1], exdir = local_tmp)
        }, error = function(e) {
          showNotification(sprintf("ZIP unpack failed: %s",
                                    conditionMessage(e)),
                           type = "error", duration = 10)
          stop(e)
        })
      })
      # Walk the unzipped tree and copy files of interest into local_tmp/
      # at the same flat layout the SSH path produces. We support any
      # directory structure inside the ZIP (with or without a wrapper dir).
      all_files <- list.files(local_tmp, recursive = TRUE, full.names = TRUE)
      copy_first <- function(pattern, dest) {
        hits <- all_files[grepl(pattern, basename(all_files), ignore.case = TRUE)]
        if (length(hits) > 0 && !file.exists(dest)) {
          file.copy(hits[1], dest, overwrite = FALSE)
        }
      }
      copy_first("^results\\.sage\\.tsv$",   file.path(local_tmp, "results.sage.tsv"))
      copy_first("^lfq\\.tsv$",              file.path(local_tmp, "lfq.tsv"))
      copy_first("^report\\.parquet$",       file.path(local_tmp, "report.parquet"))
      copy_first("^blast_results\\.tsv$",    file.path(local_tmp, "blast_results.tsv"))
      # mztab files: copy all into a single mztab/ subdir
      mztab_hits <- all_files[grepl("\\.mztab$", basename(all_files), ignore.case = TRUE)]
      if (length(mztab_hits) > 0) {
        mztab_dir <- file.path(local_tmp, "mztab")
        dir.create(mztab_dir, showWarnings = FALSE)
        for (mt in mztab_hits) {
          dest <- file.path(mztab_dir, basename(mt))
          if (!file.exists(dest)) file.copy(mt, dest)
        }
      }
      remote_dir <- local_tmp  # for display; populates values$dda_output_dir
      message(sprintf("[DDA Load] Unpacked ZIP to %s (%d files, %d mztabs)",
                       local_tmp, length(all_files), length(mztab_hits)))
    } else {
      # ── HPC-SSH path: original behavior ──────────────────────────────────
      ssh_cfg <- dda_ssh_config()
      remote_dir <- trimws(input$dda_load_path)
    }

    withProgress(message = "Loading DDA results...", value = 0.1, {
      tryCatch({

        # Try to load database search results: Sage first, then DIA-NN
        # Track which engine was used for source badge
        db_engine <- "Sage"  # default

        # --- Sage: results.sage.tsv ---
        results_remote <- file.path(remote_dir, "results.sage.tsv")
        results_local <- file.path(local_tmp, "results.sage.tsv")
        tryCatch(scp_download(ssh_cfg, results_remote, results_local), error = function(e) NULL)

        parsed <- NULL
        sage_found <- file.exists(results_local) && file.info(results_local)$size > 100

        # --- DIA-NN: report.parquet ---
        diann_remote <- file.path(remote_dir, "report.parquet")
        diann_local <- file.path(local_tmp, "report.parquet")
        tryCatch(scp_download(ssh_cfg, diann_remote, diann_local), error = function(e) NULL)
        diann_found <- file.exists(diann_local) && file.info(diann_local)$size > 100

        if (sage_found) {
          setProgress(0.3, detail = "Parsing Sage results...")

          # Check for lfq.tsv
          lfq_remote <- file.path(remote_dir, "lfq.tsv")
          lfq_local <- file.path(local_tmp, "lfq.tsv")
          tryCatch(scp_download(ssh_cfg, lfq_remote, lfq_local), error = function(e) NULL)

          # Parse
          parsed <- parse_sage_results(
            results_local,
            if (file.exists(lfq_local)) lfq_local else NULL
          )
          values$dda_sage_psms    <- parsed$psms
          values$dda_lfq_wide     <- parsed$lfq_wide
          values$dda_protein_meta <- parsed$protein_meta
          db_engine <- "Sage"
          message("[DDA Load] Loaded ", nrow(parsed$psms), " Sage PSMs")
        } else if (diann_found) {
          setProgress(0.3, detail = "Parsing DIA-NN DDA results...")

          diann_parsed <- tryCatch(
            parse_diann_dda_results(diann_local, fdr_threshold = 0.01),
            error = function(e) {
              message("[DDA Load] DIA-NN parse error: ", e$message)
              NULL
            }
          )

          if (!is.null(diann_parsed) && nrow(diann_parsed$psms) > 0) {
            # Store DIA-NN PSMs in the same slot as Sage for classification
            values$dda_sage_psms <- diann_parsed$psms
            db_engine <- "DIA-NN"
            # Create a minimal parsed list so downstream code works
            parsed <- list(psms = diann_parsed$psms)
            message("[DDA Load] Loaded ", nrow(diann_parsed$psms), " DIA-NN PSMs")
          } else {
            message("[DDA Load] DIA-NN parquet found but no passing PSMs")
          }
        } else {
          message("[DDA Load] No Sage or DIA-NN results found — loading Casanovo/BLAST only")
        }
        values$dda_output_dir    <- remote_dir
        values$dda_status        <- "loaded"
        values$dda_db_engine     <- db_engine
        values$dda_loaded_at     <- Sys.time()

        setProgress(0.5, detail = "Checking for Casanovo results...")

        # Check for Casanovo mztab files in multiple locations:
        # 1. {output_dir}/casanovo/mztab/*.mztab (DE-LIMP generated)
        # 2. {raw_data_dir}/*.mztab (pre-existing, e.g. Glendon's feather data)
        # ZIP-mode short-circuit: the unzip step already copied .mztab into local_tmp/mztab/.
        mztab_remote <- character(0)
        if (is.null(ssh_cfg)) {
          local_mztab_dir <- file.path(local_tmp, "mztab")
          if (dir.exists(local_mztab_dir)) {
            mztab_remote <- list.files(local_mztab_dir, pattern = "\\.mztab$",
                                        full.names = TRUE)
            if (length(mztab_remote) > 0)
              message("[DDA Load] Found ", length(mztab_remote),
                      " mztab files in ZIP")
          }
        } else {
          for (mztab_search_dir in c(
            file.path(remote_dir, "casanovo", "mztab"),
            file.path(remote_dir, "casanovo"),
            file.path(remote_dir, "denovo")
          )) {
            check <- ssh_exec(ssh_cfg,
              paste0("ls ", shQuote(mztab_search_dir), "/*.mztab 2>/dev/null"),
              timeout = 10)
            if (check$status == 0 && length(check$stdout) > 0) {
              found <- trimws(check$stdout)
              found <- found[nzchar(found)]
              if (length(found) > 0) {
                mztab_remote <- found
                message("[DDA Load] Found ", length(found), " mztab files in ", mztab_search_dir)
                break
              }
            }
          }
        }

        # Also check if raw data dir has mztabs (common for pre-existing Casanovo runs)
        if (length(mztab_remote) == 0) {
          # Try to find the raw data directory from search_info.md or sage config
          raw_dir_check <- ssh_exec(ssh_cfg,
            paste0("grep -r 'mztab\\|raw_data' ", shQuote(file.path(remote_dir, "search_info.md")),
                   " 2>/dev/null | head -1"),
            timeout = 10)

          # Also try parent directory and common locations
          for (search_path in c(
            dirname(remote_dir),
            "/quobyte/proteomics-grp/brett/glendon/glendon_feathers"
          )) {
            check <- ssh_exec(ssh_cfg,
              paste0("ls ", shQuote(search_path), "/*.mztab 2>/dev/null"),
              timeout = 10)
            if (check$status == 0 && length(check$stdout) > 0) {
              found <- trimws(check$stdout)
              found <- found[nzchar(found)]
              if (length(found) > 0) {
                mztab_remote <- found
                message("[DDA Load] Found ", length(found), " mztab files in ", search_path)
                break
              }
            }
          }
        }

        if (length(mztab_remote) > 0) {
          mztab_local_dir <- file.path(local_tmp, "mztab")
          dir.create(mztab_local_dir, showWarnings = FALSE)
          if (is.null(ssh_cfg)) {
            # ZIP mode: mztab_remote already contains LOCAL paths from list.files
            mztab_local <- mztab_remote
          } else {
            for (mt in mztab_remote) {
              tryCatch(scp_download(ssh_cfg, mt, file.path(mztab_local_dir, basename(mt))),
                error = function(e) NULL)
            }
            mztab_local <- list.files(mztab_local_dir, pattern = "\\.mztab$", full.names = TRUE)
          }
          if (length(mztab_local) > 0) {
            # Load ALL Casanovo PSMs (no parse-time score cutoff) so the
            # confidence slider filters the full range AND every de novo peptide
            # can be cross-referenced against the DIAMOND BLAST hits. A hard
            # 0.9 cutoff here previously dropped ~98% of PSMs (440k -> 7k),
            # leaving only short high-confidence peptides that don't BLAST — so
            # the BLAST tab read 0. The slider (default 0.9) still gates the view.
            casanovo_psms <- parse_casanovo_mztab(mztab_local, score_threshold = -Inf)
            if (nrow(casanovo_psms) > 0) {
              db_psms <- if (!is.null(parsed)) parsed$psms else NULL
              classified <- classify_dda_denovo(casanovo_psms, db_psms,
                db_engine = db_engine)
              values$dda_casanovo_psms <- casanovo_psms
              values$dda_casanovo_classification <- classified
              values$dda_casanovo_status <- "done"
              message("[DDA Load] Loaded ", nrow(casanovo_psms), " Casanovo PSMs (vs ",
                      db_engine, "), ",
                      nrow(classified$confirmed), " confirmed, ",
                      nrow(classified$novel), " novel")
            }
          }
        }

        setProgress(0.6, detail = "Checking for BLAST results...")

        # Load existing BLAST results if available (skip re-running BLAST)
        blast_remote <- file.path(remote_dir, "denovo", "blast_results.tsv")
        blast_local <- file.path(local_tmp, "blast_results.tsv")
        blast_loaded <- FALSE
        tryCatch({
          if (is.null(ssh_cfg)) {
            # ZIP mode: blast_results.tsv was already copied into local_tmp/
            # by the normalize step (if present in the ZIP). Skip SSH probe.
            blast_present <- file.exists(blast_local) && file.info(blast_local)$size > 100
          } else {
            check <- ssh_exec(ssh_cfg,
              paste("test -s", shQuote(blast_remote), "&& echo YES || echo NO"),
              timeout = 10)
            blast_present <- grepl("YES", paste(check$stdout, collapse = ""))
            if (blast_present) {
              scp_download(ssh_cfg, blast_remote, blast_local)
            }
          }
          if (blast_present) {
            if (file.exists(blast_local) && file.info(blast_local)$size > 100) {
              blast_df <- data.table::fread(blast_local, header = FALSE)
              if (nrow(blast_df) > 0) {
                col_names <- c("peptide", "subject", "pident", "length", "mismatch",
                               "gapopen", "qstart", "qend", "sstart", "send",
                               "evalue", "bitscore")
                if (ncol(blast_df) >= length(col_names)) {
                  names(blast_df)[seq_along(col_names)] <- col_names
                }
                # Add species + category columns
                blast_df$species <- sub(".*_", "", blast_df$subject)
                blast_df$category <- ifelse(blast_df$pident >= 100, "Conserved",
                  ifelse(blast_df$pident >= 90, "Near-match", "Distant"))
                # Best hit per peptide
                blast_df <- blast_df[order(-blast_df$bitscore), ]
                blast_df <- blast_df[!duplicated(blast_df$peptide), ]
                values$dda_casanovo_blast <- as.data.frame(blast_df)
                blast_loaded <- TRUE
                message("[DDA Load] Loaded ", nrow(blast_df), " BLAST hits")
              }
            }
          } else {
            message("[DDA Load] No BLAST results file at ", blast_remote)
          }
        }, error = function(e) {
          message("[DDA Load] BLAST load error: ", e$message)
        })

        # Load the LCA species-attribution table if present (nr + LCA pipeline).
        # Prefers the relaxed-evalue (*_e1) table when both exist.
        tryCatch({
          lca_local <- file.path(local_tmp, "peptide_lca.tsv")
          lca_present <- FALSE
          if (is.null(ssh_cfg)) {
            cand <- list.files(local_tmp, pattern = "_peptide_lca\\.tsv$",
                               full.names = TRUE, recursive = TRUE)
            if (length(cand) > 0) { lca_local <- tail(sort(cand), 1); lca_present <- TRUE }
          } else {
            found <- ssh_exec(ssh_cfg, paste0(
              "ls -1 ", shQuote(file.path(remote_dir, "denovo")),
              "/*_peptide_lca.tsv 2>/dev/null | sort | tail -1"), timeout = 10)
            remote_lca <- trimws(paste(found$stdout, collapse = ""))
            if (nzchar(remote_lca)) {
              scp_download(ssh_cfg, remote_lca, lca_local)
              lca_present <- file.exists(lca_local) && file.info(lca_local)$size > 50
            }
          }
          if (lca_present) {
            lca_df <- data.table::fread(lca_local, sep = "\t", header = TRUE)
            values$dda_lca <- as.data.frame(lca_df)
            message("[DDA Load] Loaded LCA species table: ", nrow(lca_df), " peptides")
          }
        }, error = function(e) message("[DDA Load] LCA load error: ", e$message))

        # Only submit new BLAST if no results loaded AND we have novel peptides
        # AND we have an SSH connection (ZIP-only mode can't submit jobs).
        if (!blast_loaded &&
            !is.null(ssh_cfg) &&
            !is.null(values$dda_casanovo_classification) &&
            nrow(values$dda_casanovo_classification$novel) > 0) {
          tryCatch({
            setProgress(0.6, detail = "Submitting DIAMOND BLAST job...")
            novel <- values$dda_casanovo_classification$novel
            # Strip modification masses from sequences — keep only amino acid letters
            clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(novel$seq_stripped))
            unique_seqs <- unique(clean_seqs[nzchar(clean_seqs)])
            message("[DDA] Submitting BLAST: ", length(unique_seqs), " novel peptides")

            # Write FASTA of novel peptides
            local_fasta <- file.path(local_tmp, "novel_peptides.fasta")
            fasta_lines <- unlist(lapply(unique_seqs, function(s) {
              c(paste0(">", s), s)
            }))
            writeLines(fasta_lines, local_fasta)

            # Upload to HPC
            remote_denovo_dir <- file.path(remote_dir, "denovo")
            ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(remote_denovo_dir)), timeout = 10)
            scp_upload(ssh_cfg, local_fasta, file.path(remote_denovo_dir, "novel_peptides.fasta"))

            # Write sbatch for DIAMOND — uses pre-built SwissProt DB
            blast_out <- file.path(remote_denovo_dir, "blast_results.tsv")
            logs_dir <- file.path(remote_dir, "logs")

            blast_sbatch <- paste0(
'#!/bin/bash -l
#SBATCH --job-name=delimp_diamond_blast
#SBATCH --partition=high
#SBATCH --account=', slurm_account, '
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output="', logs_dir, '/diamond_%j.out"
#SBATCH --error="', logs_dir, '/diamond_%j.err"

set -euo pipefail
module load diamond

echo "[DIAMOND] Start: $(date)"
echo "[DIAMOND] Novel peptides: ', length(unique_seqs), '"
echo "[DIAMOND] Database: UniProt SwissProt"

# Run BLAST against pre-built SwissProt DB
echo "[DIAMOND] Running blastp against SwissProt..."
diamond blastp \\
  --query "', file.path(remote_denovo_dir, "novel_peptides.fasta"), '" \\
  --db "', swissprot_dmnd, '" \\
  --out "', blast_out, '" \\
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
  --sensitive --id 50 --max-target-seqs 5 \\
  --threads 8 --ignore-warnings

echo "[DIAMOND] Results: $(wc -l < "', blast_out, '") hits"
echo "[DIAMOND] Done: $(date)"
')
            local_sbatch <- file.path(local_tmp, "diamond_blast.sbatch")
            writeLines(blast_sbatch, local_sbatch)
            scp_upload(ssh_cfg, local_sbatch, file.path(remote_denovo_dir, "diamond_blast.sbatch"))

            # Submit
            sbatch_path <- values$ssh_sbatch_path %||% "sbatch"
            submit_result <- ssh_exec(ssh_cfg,
              paste(sbatch_path, shQuote(file.path(remote_denovo_dir, "diamond_blast.sbatch"))),
              timeout = 15)

            if (submit_result$status == 0) {
              blast_jid <- trimws(sub(".*Submitted batch job\\s+", "",
                grep("Submitted batch job", submit_result$stdout, value = TRUE)[1]))
              values$dda_blast_job_id <- blast_jid
              message("[DDA] DIAMOND BLAST submitted: job ", blast_jid)
              showNotification(
                paste("DIAMOND BLAST submitted (job", blast_jid, ")— results will load when complete."),
                type = "message", duration = 8)

              # Push to the DDA queue so the user can see it in the table.
              # Add as a new entry rather than mutating the parent search,
              # since BLAST is a distinct compute step with its own job_id.
              values$dda_jobs <- c(values$dda_jobs %||% list(), list(list(
                job_id           = blast_jid,
                sage_job_id      = NA_character_,
                casanovo_job_id  = NA_character_,
                blast_job_id     = blast_jid,
                sage_status      = "none",
                casanovo_status  = "none",
                blast_status     = "running",
                name             = sprintf("%s — DIAMOND BLAST",
                                            basename(remote_dir)),
                output_dir       = remote_dir,
                status           = "running",
                submitted_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                n_files          = NA_integer_,
                engine           = "DIAMOND BLAST",
                casanovo         = FALSE
              )))
            }
          }, error = function(e) message("[DDA] Auto-BLAST submission failed: ", e$message))
        }

        # (mztab and BLAST loading handled above)

        setProgress(0.95, detail = "Done!")

        n_psms <- if (!is.null(parsed)) nrow(parsed$psms) else 0
        n_casanovo <- if (!is.null(values$dda_casanovo_psms)) nrow(values$dda_casanovo_psms) else 0
        n_blast <- if (!is.null(values$dda_casanovo_blast)) nrow(values$dda_casanovo_blast) else 0
        removeModal()
        parts <- character(0)
        if (n_psms > 0) {
          parts <- c(parts, sprintf("%s %s PSMs",
            format(n_psms, big.mark = ","), db_engine))
        }
        if (n_casanovo > 0) parts <- c(parts, sprintf("%s Casanovo PSMs", format(n_casanovo, big.mark = ",")))
        if (n_blast > 0) parts <- c(parts, sprintf("%s BLAST hits", format(n_blast, big.mark = ",")))
        if (length(parts) == 0) parts <- "No results found"
        showNotification(
          paste("Loaded:", paste(parts, collapse = ", ")),
          type = "message", duration = 10)

      }, error = function(e) {
        showNotification(paste("Load failed:", e$message), type = "error", duration = 15)
      })
    })
  })

  # ============================================================================
  #    Submit Sage search
  # ============================================================================
  observeEvent(input$run_dda_search, {
    req(values$ssh_connected)

    raw_dir    <- trimws(input$dda_raw_dir %||% "")
    exp_name   <- trimws(input$dda_experiment_name %||% "")

    # Resolve FASTA path from whichever source was used
    fasta_source <- input$dda_fasta_source %||% "browse"
    if (fasta_source == "uniprot") {
      fasta_path <- dda_fasta_resolved() %||% ""
    } else {
      fasta_path <- trimws(input$dda_fasta_path %||% "")
    }

    # Validation
    if (!nzchar(exp_name)) {
      showNotification("Please enter an experiment name (e.g. ocelot_dda_2026_05).",
                       type = "error", duration = 8)
      return()
    }
    if (!nzchar(raw_dir)) {
      showNotification("Please enter a raw file directory path.", type = "error")
      return()
    }
    if (!nzchar(fasta_path)) {
      showNotification("Please select or enter a FASTA file path.", type = "error")
      return()
    }
    if (is.null(values$dda_raw_files) || length(values$dda_raw_files) == 0) {
      showNotification("Please scan for files first.", type = "error")
      return()
    }

    ssh_cfg <- dda_ssh_config()

    # Build output directory on HPC
    output_dir <- file.path(
      "/quobyte/proteomics-grp/de-limp",
      ssh_cfg$user,
      "dda_output",
      gsub("[^a-zA-Z0-9_.-]", "_", exp_name)
    )

    withProgress(message = "Submitting Sage search...", value = 0.1, {
      # Create output dir + logs dir on HPC
      mkdir_result <- ssh_exec(ssh_cfg,
        paste0("mkdir -p ", shQuote(file.path(output_dir, "logs"))),
        timeout = 15)
      if (mkdir_result$status != 0) {
        showNotification("Failed to create output directory on HPC.", type = "error")
        return()
      }
      setProgress(0.2, detail = "Preparing FASTA database...")

      # Handle contaminant library — append to FASTA on HPC
      contam_lib <- input$dda_contaminant_library %||% "none"
      if (contam_lib != "none") {
        contam_result <- get_contaminant_fasta(contam_lib)
        if (contam_result$success) {
          # Upload contaminant FASTA to HPC
          remote_contam_dir <- file.path(output_dir, "databases")
          remote_contam_path <- file.path(remote_contam_dir, basename(contam_result$path))

          exists_check <- ssh_exec(ssh_cfg,
            paste("test -f", shQuote(remote_contam_path), "&& echo EXISTS"))
          if (!any(grepl("EXISTS", exists_check$stdout))) {
            ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(remote_contam_dir)))
            scp_upload(ssh_cfg, contam_result$path, remote_contam_path)
          }

          # Concatenate proteome + contaminant into combined FASTA on HPC
          # Sage takes a single FASTA path, unlike DIA-NN which accepts multiple --fasta args
          combined_fasta <- file.path(output_dir, "databases",
            paste0("combined_", basename(fasta_path)))
          ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(dirname(combined_fasta))))
          cat_result <- ssh_exec(ssh_cfg,
            paste("cat", shQuote(fasta_path), shQuote(remote_contam_path),
              ">", shQuote(combined_fasta)),
            timeout = 30)
          if (cat_result$status == 0) {
            message("[DDA] Combined FASTA: ", combined_fasta,
              " (proteome + ", contam_lib, " contaminants)")
            fasta_path <- combined_fasta
          } else {
            showNotification("Warning: Could not append contaminant library. Using proteome only.",
              type = "warning")
          }
        } else {
          showNotification(paste("Warning: Contaminant library not found:", contam_result$error),
            type = "warning")
        }
      }

      setProgress(0.3, detail = "Generating Sage config...")

      # Generate sage.json locally, then upload
      local_tmp <- tempdir()
      raw_paths <- values$dda_raw_files

      config_path_local <- generate_sage_config(
        fasta_path       = fasta_path,
        raw_paths        = raw_paths,
        output_dir       = output_dir,
        preset           = input$dda_preset %||% "standard",
        missed_cleavages = input$dda_missed_cleavages %||% 2,
        precursor_tol_ppm = input$dda_precursor_tol %||% 20,
        fragment_tol_ppm  = input$dda_fragment_tol %||% 20,
        min_peaks         = 6
      )

      # Upload sage.json to HPC
      remote_config <- file.path(output_dir, "sage.json")
      scp_result <- scp_upload(ssh_cfg, config_path_local, remote_config)
      if (scp_result$status != 0) {
        showNotification("Failed to upload Sage config to HPC.", type = "error")
        return()
      }
      setProgress(0.4, detail = "Generating sbatch script...")

      # Generate sbatch script
      sbatch_content <- generate_sage_sbatch(
        sage_bin        = sage_bin,
        config_path     = remote_config,
        raw_dir         = raw_dir,
        output_dir      = output_dir,
        experiment_name = exp_name,
        cpus            = input$dda_cpus %||% 32,
        mem_gb          = input$dda_mem %||% 64,
        time_limit      = input$dda_time_limit %||% "02:00:00",
        account         = slurm_account,
        partition       = slurm_partition
      )

      # Write sbatch script locally, upload
      local_sbatch <- file.path(local_tmp, "sage_search.sbatch")
      writeLines(sbatch_content, local_sbatch)
      remote_sbatch <- file.path(output_dir, "sage_search.sbatch")
      scp_upload(ssh_cfg, local_sbatch, remote_sbatch)

      setProgress(0.6, detail = "Submitting to SLURM...")

      # Submit via sbatch
      sbatch_path <- values$ssh_sbatch_path %||% "sbatch"
      submit_result <- ssh_exec(ssh_cfg,
        paste(sbatch_path, shQuote(remote_sbatch)),
        timeout = 30)

      if (submit_result$status != 0) {
        showNotification(
          paste("sbatch submission failed:", paste(submit_result$stdout, collapse = " ")),
          type = "error")
        return()
      }

      # Parse job ID from "Submitted batch job 12345"
      job_line <- grep("Submitted batch job", submit_result$stdout, value = TRUE)
      if (length(job_line) == 0) {
        showNotification("Could not parse job ID from sbatch output.", type = "error")
        return()
      }
      job_id <- trimws(sub(".*Submitted batch job\\s+", "", job_line[1]))
      message("[DDA] Sage job submitted: ", job_id)

      # ─── Write search_info.md on Hive ─────────────────────────────────────
      # Mirrors what generate_search_info() does for DIA. This is the file
      # Discover-from-Hive parses to populate the queue table (job ID,
      # submitted, n_files, etc.) AND the View Info modal reads later.
      info_md <- sprintf("# DDA Search Info

**Pipeline**: Sage + Casanovo (DE-LIMP DDA mode)
**Submitted**: %s
**App version**: DE-LIMP %s

## Search Parameters
- **Search engine**: Sage (binary: %s)
- **Preset**: %s
- **FASTA**: `%s`
- **Contaminant library**: %s
- **Enzyme**: Trypsin/P, %d missed cleavages
- **Precursor tolerance**: ±%s ppm
- **Fragment tolerance**: ±%s ppm
- **Normalization**: %s
- **Imputation**: %s
- **Casanovo enabled**: %s
- **Casanovo model**: %s

## Files
- **Raw files**: %d files in `%s`
- **Output**: `%s`
- **Sage job ID**: %s
- **Casanovo job ID**: (pending — chained submission)
",
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        as.character(values$app_version %||% "v3.11.0"),
        basename(sage_bin),
        input$dda_preset %||% "standard",
        basename(fasta_path %||% "(unspecified)"),
        contam_lib %||% "(none)",
        as.integer(input$dda_missed_cleavages %||% 2),
        as.character(input$dda_precursor_tol %||% 20),
        as.character(input$dda_fragment_tol %||% 20),
        as.character(input$dda_norm_method %||% "cyclicloess"),
        as.character(input$dda_impute_method %||% "perseus"),
        if (isTRUE(input$dda_run_casanovo)) "TRUE" else "FALSE",
        if (isTRUE(input$dda_run_casanovo))
          as.character(input$dda_casanovo_model %||% "casanovo_v4_2_0")
        else "(disabled)",
        length(raw_paths), raw_dir, output_dir, job_id)
      tryCatch({
        local_info <- file.path(tempdir(), "search_info.md")
        writeLines(info_md, local_info)
        scp_upload(ssh_cfg, local_info, file.path(output_dir, "search_info.md"))
        unlink(local_info)
      }, error = function(e) message("[DDA] could not write search_info.md: ",
                                     conditionMessage(e)))

      # Store state — both the "active job" scalars (used by polling + status
      # UI + result loaders) AND a queue entry (persisted to disk via the
      # observer at module init, so the queue survives Shiny restarts).
      values$dda_job_id     <- job_id
      values$dda_output_dir <- output_dir
      values$dda_status     <- "running"

      # Use whichever name the user actually typed in the UI.
      # (input ID is dda_experiment_name — was previously misnamed dda_analysis_name here.)
      dda_name <- if (nzchar(exp_name)) exp_name else
                    sprintf("DDA_%s", format(Sys.time(), "%Y%m%d_%H%M%S"))
      values$dda_jobs <- c(values$dda_jobs %||% list(), list(list(
        job_id           = job_id,
        sage_job_id      = job_id,
        casanovo_job_id  = NA_character_,  # populated when casanovo sbatch lands
        sage_status      = "running",
        casanovo_status  = if (isTRUE(input$dda_run_casanovo)) "pending" else "none",
        name             = dda_name,
        output_dir       = output_dir,
        status           = "running",
        submitted_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        n_files          = length(raw_paths),
        engine           = "Sage",
        casanovo         = isTRUE(input$dda_run_casanovo),
        casanovo_model = if (isTRUE(input$dda_run_casanovo))
                           input$dda_casanovo_model %||% NA_character_
                         else NA_character_,
        # Analysis mode drives results-page routing + ZIP export contents.
        # Comes straight from the preset selector. New value → no edits anywhere else.
        mode             = input$dda_preset %||% "standard"
      )))

      run_casanovo <- isTRUE(input$dda_run_casanovo)
      values$dda_search_params <- list(
        preset              = input$dda_preset %||% "standard",
        fasta_path          = fasta_path,
        raw_dir             = raw_dir,
        n_files             = length(raw_paths),
        missed_cleavages    = input$dda_missed_cleavages %||% 2,
        precursor_tol       = input$dda_precursor_tol %||% 20,
        fragment_tol        = input$dda_fragment_tol %||% 20,
        normalization       = input$dda_norm_method %||% "cyclicloess",
        imputation          = input$dda_impute_method %||% "perseus",
        min_valid           = input$dda_min_valid %||% 0.5,
        contaminant_library = contam_lib,
        submitted_at        = Sys.time(),
        sage_bin            = sage_bin,
        casanovo_enabled    = run_casanovo
      )

      # --- Write search_info.md to output directory ---
      tryCatch({
        sp <- values$dda_search_params
        search_info <- paste0(
          "# DDA Search Info\n\n",
          "**Pipeline**: Sage + Casanovo (DE-LIMP DDA mode)\n",
          "**Submitted**: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
          "**App version**: DE-LIMP v", values$app_version %||% "unknown", "\n\n",
          "## Search Parameters\n\n",
          "- **Search engine**: Sage v0.14.6\n",
          "- **Preset**: ", sp$preset, "\n",
          "- **FASTA**: `", basename(sp$fasta_path), "`\n",
          "- **Contaminant library**: ", sp$contaminant_library, "\n",
          "- **Enzyme**: Trypsin/P, ", sp$missed_cleavages, " missed cleavages\n",
          "- **Precursor tolerance**: ±", sp$precursor_tol, " ppm\n",
          "- **Fragment tolerance**: ±", sp$fragment_tol, " ppm\n",
          "- **Casanovo enabled**: ", sp$casanovo_enabled, "\n\n",
          "## Files\n\n",
          "- **Raw files**: ", sp$n_files, " files in `", sp$raw_dir, "`\n",
          "- **Output**: `", output_dir, "`\n",
          "- **Sage job ID**: ", job_id, "\n",
          if (run_casanovo) paste0("- **Casanovo job ID**: (pending)\n") else "",
          "\n## File List\n\n",
          paste0("- `", basename(raw_paths), "`\n", collapse = ""),
          "\n## Quantification Settings\n\n",
          "- **Normalization**: ", sp$normalization, "\n",
          "- **Imputation**: ", sp$imputation, "\n",
          "- **Min valid fraction**: ", sp$min_valid, "\n",
          "\n**Log files**: `", output_dir, "/logs/`\n"
        )
        local_info <- file.path(tempdir(), "search_info.md")
        writeLines(search_info, local_info)
        scp_upload(ssh_cfg, local_info, file.path(output_dir, "search_info.md"))
        message("[DDA] search_info.md written to ", output_dir)
      }, error = function(e) message("[DDA] search_info.md write failed: ", e$message))

      # --- Casanovo submission (optional, GPU) ---
      if (run_casanovo) {
        setProgress(0.7, detail = "Submitting Casanovo de novo...")

        tryCatch({
          # Resolve env + ckpt + CLI version from the user's model pick.
          # Casanovo v5 ships a newer depthcharge with `depthcharge.tokenizers`
          # that the v4 env (cassonovo_env) cannot import — they MUST be paired.
          model_id <- input$dda_casanovo_model %||% "casanovo_v4_2_0"
          is_v5 <- grepl("^casanovo_v5", model_id, ignore.case = TRUE)
          casanovo_env_use <- if (is_v5) {
            "/quobyte/proteomics-grp/conda_envs/casanovo5"
          } else {
            casanovo_conda_env
          }
          casanovo_ckpt_use <- file.path(
            "/quobyte/proteomics-grp/bioinformatics_programs/casanovo_modles",
            paste0(model_id, ".ckpt")
          )
          casanovo_version_use <- if (is_v5) "v5" else "v4"
          message("[DDA] Casanovo route: model=", model_id,
                  " version=", casanovo_version_use,
                  " env=", casanovo_env_use,
                  " ckpt=", casanovo_ckpt_use)

          compute_mode_use <- input$dda_casanovo_compute %||% "gpu"
          message("[DDA] Casanovo compute_mode = ", compute_mode_use)

          casanovo_scripts <- generate_casanovo_sbatch(
            raw_dir          = raw_dir,
            output_dir       = output_dir,
            experiment_name  = exp_name,
            conda_env_path   = casanovo_env_use,
            model_ckpt       = casanovo_ckpt_use,
            casanovo_version = casanovo_version_use,
            compute_mode     = compute_mode_use,
            converter_script = casanovo_converter,
            n_files          = length(raw_paths),
            account          = slurm_account,
            gpu_partition    = casanovo_gpu_partition,
            gpu_qos          = casanovo_gpu_qos
          )

          # Create casanovo subdirs on HPC
          ssh_exec(ssh_cfg,
            paste0("mkdir -p ",
              shQuote(casanovo_scripts$mgf_dir), " ",
              shQuote(casanovo_scripts$mztab_dir)),
            timeout = 15)

          # Upload bruker_to_mgf.py converter to HPC
          local_converter <- file.path(getwd(), "python", "bruker_to_mgf.py")
          if (file.exists(local_converter)) {
            scp_upload(ssh_cfg, local_converter, casanovo_converter)
          }

          # Write and upload sbatch scripts
          local_convert_sbatch <- file.path(local_tmp, "casanovo_convert.sbatch")
          writeLines(casanovo_scripts$convert_script, local_convert_sbatch)
          remote_convert_sbatch <- file.path(output_dir, "casanovo_convert.sbatch")
          scp_upload(ssh_cfg, local_convert_sbatch, remote_convert_sbatch)

          local_casanovo_sbatch <- file.path(local_tmp, "casanovo_sequence.sbatch")
          writeLines(casanovo_scripts$casanovo_script, local_casanovo_sbatch)
          remote_casanovo_sbatch <- file.path(output_dir, "casanovo_sequence.sbatch")
          scp_upload(ssh_cfg, local_casanovo_sbatch, remote_casanovo_sbatch)

          # Write and upload launcher script.
          # Pass Sage's job_id so the convert step waits — Sage now pre-converts
          # .raw → mzML and Casanovo reads those mzML files.
          launcher_content <- generate_casanovo_launcher(
            remote_convert_sbatch, remote_casanovo_sbatch,
            sage_job_id = job_id)
          local_launcher <- file.path(local_tmp, "casanovo_submit.sh")
          writeLines(launcher_content, local_launcher)
          remote_launcher <- file.path(output_dir, "casanovo_submit.sh")
          scp_upload(ssh_cfg, local_launcher, remote_launcher)

          # Submit Casanovo pipeline
          setProgress(0.85, detail = "Submitting Casanovo to GPU queue...")
          casanovo_submit <- ssh_exec(ssh_cfg,
            paste("bash -l", shQuote(remote_launcher)),
            timeout = 30)

          if (casanovo_submit$status == 0) {
            # Parse job IDs from launcher output
            convert_line <- grep("^CONVERT:", casanovo_submit$stdout, value = TRUE)
            casanovo_line <- grep("^CASANOVO:", casanovo_submit$stdout, value = TRUE)

            convert_jid <- if (length(convert_line) > 0)
              trimws(sub("^CONVERT:", "", convert_line[1])) else NULL
            casanovo_jid <- if (length(casanovo_line) > 0)
              trimws(sub("^CASANOVO:", "", casanovo_line[1])) else NULL

            values$dda_casanovo_convert_job_id <- convert_jid
            values$dda_casanovo_job_id  <- casanovo_jid
            values$dda_casanovo_status  <- "running"
            values$dda_casanovo_mztab_dir <- casanovo_scripts$mztab_dir

            message("[DDA] Casanovo MGF convert job: ", convert_jid,
                    ", Casanovo sequence job: ", casanovo_jid)
            showNotification(
              paste("Casanovo submitted! Convert:", convert_jid,
                    "| Sequence:", casanovo_jid),
              type = "message", duration = 10)

            # ── Auto-chain DIAMOND BLAST after Casanovo (afterok dependency) ──
            # Runs against ALL Casanovo peptides (over-broad — DE-LIMP filters
            # to "novel only" at Load time). Tiny extra compute, but means BLAST
            # is finished by the time you click Load Results, no second wait.
            blast_jid_chain <- tryCatch({
              denovo_dir   <- file.path(output_dir, "denovo")
              mztab_dir    <- file.path(output_dir, "casanovo", "mztab")
              chain_blast_sbatch <- paste0('#!/bin/bash -l
#SBATCH --job-name=delimp_diamond_chained
#SBATCH --partition=high
#SBATCH --account=', slurm_account, '
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output="', output_dir, '/logs/diamond_chain_%j.out"
#SBATCH --error="', output_dir, '/logs/diamond_chain_%j.err"

set -euo pipefail
shopt -s nullglob
module load diamond

MZTAB_DIR="', mztab_dir, '"
OUT_DIR="', denovo_dir, '"
SWISSPROT_DB="', swissprot_dmnd, '"
TREMBL_DB="', trembl_dmnd, '"

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# Extract unique peptide sequences from all Casanovo .mztab files.
# mztab PSM rows start with "PSM"; sequence is column 2 (TSV).
echo "[DIAMOND-chain] extracting peptides from $MZTAB_DIR"
awk -F"\t" \'$1 == "PSM" { gsub(/\\[[^]]*\\]/, "", $2); gsub(/[^A-Z]/, "", $2); if (length($2) >= 5) print $2 }\' "$MZTAB_DIR"/*.mztab \
  | sort -u > all_casanovo_peptides.txt
N_PEP=$(wc -l < all_casanovo_peptides.txt)
echo "[DIAMOND-chain] $N_PEP unique peptides"
awk \'{ print ">" $1; print $1 }\' all_casanovo_peptides.txt > all_casanovo_peptides.fasta

# ── STEP 1: SwissProt (fast — most peptides hit here) ─────────────────
echo "[DIAMOND-chain] step 1/2: SwissProt"
diamond blastp \\
  -d "$SWISSPROT_DB" \\
  -q all_casanovo_peptides.fasta \\
  -o blast_results_swissprot.tsv \\
  --outfmt 6 \\
  --threads 8 \\
  --max-target-seqs 25 \\
  --evalue 1e-3 \\
  --ignore-warnings
N_SP_HITS=$(wc -l < blast_results_swissprot.tsv)
N_SP_PEPS=$(cut -f1 blast_results_swissprot.tsv | sort -u | wc -l)
echo "[DIAMOND-chain] SwissProt: $N_SP_HITS hits across $N_SP_PEPS peptides"

# ── STEP 2: TrEMBL on ALL peptides (not just SwissProt-misses) ─────────
# Search every peptide against TrEMBL too. Combined with SwissProt + an LCA
# over the top hits, this avoids the model-organism skew you get when a
# "good enough" SwissProt hit (mouse/human/bovine) freezes the call before the
# real species-specific TrEMBL match (e.g. cat/felid) is ever considered.
# Keep up to 25 hits/peptide so the LCA sees the species spread, not one best hit.
echo "[DIAMOND-chain] step 2/2: TrEMBL on ALL peptides"
diamond blastp \\
  -d "$TREMBL_DB" \\
  -q all_casanovo_peptides.fasta \\
  -o blast_results_trembl.tsv \\
  --outfmt 6 \\
  --threads 8 \\
  --max-target-seqs 25 \\
  --evalue 1e-3 \\
  --ignore-warnings
N_TR_HITS=$(wc -l < blast_results_trembl.tsv)
N_TR_PEPS=$(cut -f1 blast_results_trembl.tsv | sort -u | wc -l)
echo "[DIAMOND-chain] TrEMBL: $N_TR_HITS hits across $N_TR_PEPS peptides"

# ── Combine into a single blast_results.tsv with a db source column ───
# Source-tag each hit so downstream code can tell DBs apart.
awk -v db=SwissProt \'BEGIN{OFS="\t"} { print $0, db }\' blast_results_swissprot.tsv  > blast_results.tsv
awk -v db=TrEMBL    \'BEGIN{OFS="\t"} { print $0, db }\' blast_results_trembl.tsv    >> blast_results.tsv

echo "[DIAMOND-chain] done. Combined blast_results.tsv: $(wc -l < blast_results.tsv) hits"
echo "[DIAMOND-chain] summary: $N_SP_HITS SwissProt + $N_TR_HITS TrEMBL hits"
')
              local_chain <- file.path(local_tmp, "diamond_chain.sbatch")
              writeLines(chain_blast_sbatch, local_chain)
              remote_chain <- file.path(output_dir, "diamond_chain.sbatch")
              scp_upload(ssh_cfg, local_chain, remote_chain)
              chain_res <- ssh_exec(ssh_cfg,
                paste(sbatch_path, "--dependency=afterok:", shQuote(casanovo_jid), " ",
                       shQuote(remote_chain)),
                timeout = 30)
              chain_jid_line <- grep("Submitted batch job", chain_res$stdout, value = TRUE)
              if (length(chain_jid_line)) {
                cj <- trimws(sub(".*Submitted batch job\\s+", "", chain_jid_line[1]))
                values$dda_blast_job_id <- cj
                message("[DDA] DIAMOND BLAST chained: job ", cj,
                        " (afterok:", casanovo_jid, ")")
                cj
              } else NA_character_
            }, error = function(e) {
              message("[DDA] Auto-chain BLAST failed: ", conditionMessage(e))
              NA_character_
            })
          } else {
            message("[DDA] Casanovo submission failed: ",
                    paste(casanovo_submit$stdout, collapse = " "))
            showNotification(
              "Casanovo submission failed. Sage search continues.",
              type = "warning", duration = 10)
            values$dda_casanovo_status <- "error"
          }
        }, error = function(e) {
          message("[DDA] Casanovo submission error: ", e$message)
          showNotification(
            paste("Casanovo error:", e$message, "- Sage search continues."),
            type = "warning", duration = 10)
          values$dda_casanovo_status <- "error"
        })
      } else {
        values$dda_casanovo_status <- "disabled"
      }

      setProgress(1.0, detail = "Job(s) submitted!")
      msg <- paste("Sage search submitted! Job ID:", job_id)
      if (run_casanovo && !is.null(values$dda_casanovo_job_id)) {
        msg <- paste(msg, "| Casanovo:", values$dda_casanovo_job_id)
      }
      showNotification(msg, type = "message", duration = 10)
    })
  })

  # ============================================================================
  #    Job polling (every 15 seconds when a job is running)
  # ============================================================================
  observe({
    req(values$dda_status == "running", values$dda_job_id, values$ssh_connected)
    invalidateLater(15000)

    ssh_cfg <- isolate(dda_ssh_config())
    job_id  <- isolate(values$dda_job_id)

    result <- tryCatch(
      ssh_exec(ssh_cfg,
        paste0("sacct -j ", job_id, " --format=JobID,State --noheader --parsable2"),
        timeout = 15),
      error = function(e) list(status = 1, stdout = character(0))
    )

    if (result$status != 0 || length(result$stdout) == 0) return()

    # Parse SLURM state -- filter out .extern/.batch substeps
    lines <- trimws(result$stdout)
    lines <- lines[nzchar(lines)]
    main_lines <- lines[!grepl("\\.", lines)]  # exclude substeps

    if (length(main_lines) == 0) return()

    # Get state from the main job line
    parts <- strsplit(main_lines[1], "\\|")[[1]]
    if (length(parts) < 2) return()
    state <- trimws(parts[2])

    # Helper: write back per-engine status to the matching queue row so the
    # Sage column flips "running" → "done"/"failed" (without it the queue
    # display stays stale even after Sage finishes — observed v3.11.2/3).
    update_queue_sage_status <- function(new_status) {
      jobs <- isolate(values$dda_jobs %||% list())
      idx <- which(vapply(jobs, function(j)
                          identical(as.character(j$job_id),
                                    as.character(job_id)),
                          logical(1)))
      if (length(idx) == 0) return()
      jobs[[idx[1]]]$sage_status <- new_status
      values$dda_jobs <- jobs
    }

    if (state %in% c("COMPLETED")) {
      message("[DDA] Sage job completed: ", job_id)
      update_queue_sage_status("done")
      values$dda_status <- "loading"
      showNotification("Sage search completed! Loading results...",
        type = "message", duration = 8)
      # Trigger result loading
      load_sage_results_from_hpc()
    } else if (state %in% c("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED", "NODE_FAIL")) {
      message("[DDA] Sage job failed: ", job_id, " (", state, ")")
      update_queue_sage_status("failed")
      values$dda_status <- "error"
      showNotification(
        paste("Sage search failed:", state),
        type = "error", duration = 15)
    }
    # PENDING, RUNNING, COMPLETING -> keep polling
  })

  # ============================================================================
  #    Casanovo job polling (every 15 seconds when running)
  # ============================================================================
  observe({
    req(values$dda_casanovo_status == "running",
        values$dda_casanovo_job_id,
        values$ssh_connected)
    invalidateLater(15000)

    ssh_cfg      <- isolate(dda_ssh_config())
    casanovo_jid <- isolate(values$dda_casanovo_job_id)

    # Check the array job status
    result <- tryCatch(
      ssh_exec(ssh_cfg,
        paste0("sacct -j ", casanovo_jid,
               " --format=JobID,State --noheader --parsable2"),
        timeout = 15),
      error = function(e) list(status = 1, stdout = character(0))
    )

    if (result$status != 0 || length(result$stdout) == 0) return()

    lines <- trimws(result$stdout)
    lines <- lines[nzchar(lines)]
    # For array jobs: filter to task lines (contain _) but not substeps (contain .)
    task_lines <- lines[grepl("_", lines) & !grepl("\\.", lines)]
    if (length(task_lines) == 0) {
      # Not an array yet or single job — check main line
      main_lines <- lines[!grepl("[_.]", lines)]
      if (length(main_lines) == 0) return()
      parts <- strsplit(main_lines[1], "\\|")[[1]]
      if (length(parts) < 2) return()
      state <- trimws(parts[2])

      if (state %in% c("PENDING")) return()  # still queued
      if (state %in% c("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED", "NODE_FAIL")) {
        message("[DDA] Casanovo job failed: ", casanovo_jid, " (", state, ")")
        # Sync queue row before we req() ourselves out of further polling
        jobs <- isolate(values$dda_jobs %||% list())
        sage_id <- isolate(values$dda_job_id)
        idx <- which(vapply(jobs, function(j)
                            identical(as.character(j$job_id),
                                      as.character(sage_id)),
                            logical(1)))
        if (length(idx) > 0) {
          jobs[[idx[1]]]$casanovo_status <- "failed"
          values$dda_jobs <- jobs
        }
        values$dda_casanovo_status <- "error"
        showNotification(
          paste("Casanovo failed:", state, "- Sage results still available."),
          type = "warning", duration = 10)
        return()
      }
      return()  # RUNNING
    }

    # Parse array task states
    task_states <- vapply(task_lines, function(l) {
      parts <- strsplit(l, "\\|")[[1]]
      if (length(parts) >= 2) trimws(parts[2]) else "UNKNOWN"
    }, character(1))

    n_completed <- sum(task_states == "COMPLETED")
    n_failed    <- sum(task_states %in% c("FAILED", "TIMEOUT", "OUT_OF_MEMORY"))
    n_total     <- length(task_states)
    n_pending   <- sum(task_states %in% c("PENDING", "RUNNING", "COMPLETING"))

    # Update progress message
    message(sprintf("[DDA] Casanovo progress: %d/%d completed, %d failed, %d pending",
      n_completed, n_total, n_failed, n_pending))

    # Helper: write back per-engine status to the matching queue row, keyed
    # on the Sage job_id (the queue row's primary key, not casanovo_jid).
    update_queue_casanovo_status <- function(new_status) {
      jobs <- isolate(values$dda_jobs %||% list())
      sage_id <- isolate(values$dda_job_id)
      idx <- which(vapply(jobs, function(j)
                          identical(as.character(j$job_id),
                                    as.character(sage_id)),
                          logical(1)))
      if (length(idx) == 0) return()
      jobs[[idx[1]]]$casanovo_status <- new_status
      values$dda_jobs <- jobs
    }

    if (n_pending == 0) {
      # All tasks finished
      if (n_completed > 0) {
        message("[DDA] Casanovo completed: ", n_completed, "/", n_total, " tasks")
        update_queue_casanovo_status(if (n_failed > 0) "partial" else "done")
        values$dda_casanovo_status <- "loading"
        showNotification(
          paste("Casanovo completed!", n_completed, "/", n_total, "files"),
          type = "message", duration = 8)
        # Trigger Casanovo result loading
        load_casanovo_results_from_hpc()
      } else {
        update_queue_casanovo_status("failed")
        values$dda_casanovo_status <- "error"
        showNotification("All Casanovo tasks failed.", type = "warning")
      }
    }
  })

  # ============================================================================
  #    Load Casanovo results from HPC
  # ============================================================================
  load_casanovo_results_from_hpc <- function() {
    ssh_cfg   <- dda_ssh_config()
    mztab_dir <- values$dda_casanovo_mztab_dir

    if (is.null(mztab_dir)) {
      values$dda_casanovo_status <- "error"
      return()
    }

    withProgress(message = "Loading Casanovo results...", value = 0.1, {
      # List mztab files on HPC
      list_result <- ssh_exec(ssh_cfg,
        paste0("ls -1 ", shQuote(mztab_dir), "/*.mztab 2>/dev/null"),
        timeout = 15)

      if (list_result$status != 0 || length(list_result$stdout) == 0) {
        showNotification("No Casanovo .mztab files found.", type = "warning")
        values$dda_casanovo_status <- "error"
        return()
      }

      remote_mztabs <- trimws(list_result$stdout)
      remote_mztabs <- remote_mztabs[nzchar(remote_mztabs)]
      message("[DDA] Found ", length(remote_mztabs), " Casanovo .mztab files")

      setProgress(0.3, detail = paste("Downloading", length(remote_mztabs), "files..."))

      # Download all mztab files
      local_mztab_dir <- file.path(tempdir(), "casanovo_mztab")
      dir.create(local_mztab_dir, recursive = TRUE, showWarnings = FALSE)

      local_paths <- character(0)
      for (remote_path in remote_mztabs) {
        local_path <- file.path(local_mztab_dir, basename(remote_path))
        dl <- tryCatch(
          scp_download(ssh_cfg, remote_path, local_path),
          error = function(e) list(status = 1)
        )
        if (dl$status == 0) {
          local_paths <- c(local_paths, local_path)
        }
      }

      if (length(local_paths) == 0) {
        showNotification("Failed to download Casanovo results.", type = "error")
        values$dda_casanovo_status <- "error"
        return()
      }

      setProgress(0.6, detail = "Parsing mzTab files...")

      # Parse mzTab files
      casanovo_psms <- tryCatch(
        parse_casanovo_mztab(local_paths),
        error = function(e) {
          message("[DDA] Casanovo parse error: ", e$message)
          showNotification(paste("Casanovo parse error:", e$message), type = "error")
          NULL
        }
      )

      if (is.null(casanovo_psms) || nrow(casanovo_psms) == 0) {
        values$dda_casanovo_status <- "error"
        return()
      }

      # Determine which database search engine is loaded
      db_engine <- values$dda_db_engine %||% "Sage"
      setProgress(0.8, detail = paste0("Cross-referencing with ", db_engine, "..."))

      # Store raw Casanovo results
      values$dda_casanovo_psms <- casanovo_psms

      # Cross-reference with database search if available
      if (!is.null(values$dda_sage_psms)) {
        classification <- tryCatch(
          classify_dda_denovo(casanovo_psms, values$dda_sage_psms,
            db_engine = db_engine),
          error = function(e) {
            message("[DDA] Classification error: ", e$message)
            NULL
          }
        )

        if (!is.null(classification)) {
          values$dda_casanovo_classification <- classification
          message(sprintf(
            "[DDA] Casanovo classification (vs %s): %d confirmed, %d novel",
            db_engine,
            classification$summary_stats$n_confirmed,
            classification$summary_stats$n_novel
          ))
        }
      }

      # Auto-submit DIAMOND BLAST on novel peptides (score >= 0.8)
      if (!is.null(values$dda_casanovo_classification) &&
          nrow(values$dda_casanovo_classification$novel) > 0) {
        tryCatch({
          setProgress(0.9, detail = "Submitting DIAMOND BLAST...")
          novel <- values$dda_casanovo_classification$novel
          # Filter to score >= 0.8 for BLAST, strip to pure amino acid letters
          novel_hc <- novel[novel$score >= 0.8, ]
          clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(novel_hc$seq_stripped))
          unique_seqs <- unique(clean_seqs[nzchar(clean_seqs)])
          message("[DDA] Submitting BLAST: ", length(unique_seqs), " novel peptides (score >= 0.8)")

          local_fasta <- file.path(tempdir(), "novel_peptides.fasta")
          fasta_lines <- unlist(lapply(unique_seqs, function(s) {
            c(paste0(">", s), s)
          }))
          writeLines(fasta_lines, local_fasta)

          remote_denovo_dir <- file.path(values$dda_output_dir, "denovo")
          ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(remote_denovo_dir)), timeout = 10)
          scp_upload(ssh_cfg, local_fasta, file.path(remote_denovo_dir, "novel_peptides.fasta"))

          blast_out <- file.path(remote_denovo_dir, "blast_results.tsv")
          logs_dir <- file.path(values$dda_output_dir, "logs")

          blast_sbatch <- paste0(
'#!/bin/bash -l
#SBATCH --job-name=delimp_diamond_blast
#SBATCH --partition=high
#SBATCH --account=', slurm_account, '
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output="', logs_dir, '/diamond_%j.out"
#SBATCH --error="', logs_dir, '/diamond_%j.err"

set -euo pipefail
module load diamond
echo "[DIAMOND] Start: $(date)"
echo "[DIAMOND] Novel peptides (score>=0.8): ', length(unique_seqs), '"
echo "[DIAMOND] Database: UniProt SwissProt"

# Run BLAST against pre-built SwissProt DB
diamond blastp \\
  --query "', file.path(remote_denovo_dir, "novel_peptides.fasta"), '" \\
  --db "', swissprot_dmnd, '" \\
  --out "', blast_out, '" \\
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \\
  --sensitive --id 50 --max-target-seqs 5 --threads 8 --ignore-warnings

echo "[DIAMOND] Hits: $(wc -l < "', blast_out, '")"
echo "[DIAMOND] Done: $(date)"
')
          local_sbatch <- file.path(tempdir(), "diamond_blast.sbatch")
          writeLines(blast_sbatch, local_sbatch)
          scp_upload(ssh_cfg, local_sbatch, file.path(remote_denovo_dir, "diamond_blast.sbatch"))

          sbatch_path <- values$ssh_sbatch_path %||% "sbatch"
          submit_result <- ssh_exec(ssh_cfg,
            paste(sbatch_path, shQuote(file.path(remote_denovo_dir, "diamond_blast.sbatch"))),
            timeout = 15)
          if (submit_result$status == 0) {
            blast_jid <- trimws(sub(".*Submitted batch job\\s+", "",
              grep("Submitted batch job", submit_result$stdout, value = TRUE)[1]))
            values$dda_blast_job_id <- blast_jid
            message("[DDA] DIAMOND BLAST submitted: ", blast_jid)
          }
        }, error = function(e) message("[DDA] Auto-BLAST failed: ", e$message))
      }

      setProgress(1.0, detail = "Done!")
      values$dda_casanovo_status <- "done"
      showNotification(
        paste("Casanovo loaded:", nrow(casanovo_psms), "de novo sequences"),
        type = "message", duration = 10)
    })
  }

  # ============================================================================
  #    Load results from HPC
  # ============================================================================
  load_sage_results_from_hpc <- function() {
    ssh_cfg    <- dda_ssh_config()
    output_dir <- values$dda_output_dir

    withProgress(message = "Loading Sage results...", value = 0.1, {
      local_tmp <- file.path(tempdir(), "sage_results")
      dir.create(local_tmp, recursive = TRUE, showWarnings = FALSE)

      # Download results files
      files_to_get <- c("results.sage.tsv", "lfq.tsv")
      for (f in files_to_get) {
        remote_path <- file.path(output_dir, f)
        local_path  <- file.path(local_tmp, f)
        dl <- scp_download(ssh_cfg, remote_path, local_path)
        if (dl$status != 0) {
          showNotification(paste("Failed to download", f, "from HPC."), type = "error")
          values$dda_status <- "error"
          return()
        }
      }
      setProgress(0.4, detail = "Parsing Sage output...")

      # Also try to download sage report JSON (may have various names)
      for (rpt in c("results.json", "sage_report.json")) {
        remote_rpt <- file.path(output_dir, rpt)
        local_rpt  <- file.path(local_tmp, rpt)
        tryCatch(scp_download(ssh_cfg, remote_rpt, local_rpt), error = function(e) NULL)
      }

      # Parse results
      parsed <- tryCatch(
        parse_sage_results(
          results_path = file.path(local_tmp, "results.sage.tsv"),
          lfq_path     = file.path(local_tmp, "lfq.tsv")
        ),
        error = function(e) {
          showNotification(paste("Error parsing Sage output:", e$message), type = "error")
          NULL
        }
      )

      if (is.null(parsed)) {
        values$dda_status <- "error"
        return()
      }

      setProgress(0.6, detail = "Storing results...")
      values$dda_sage_psms    <- parsed$psms
      values$dda_lfq_wide     <- parsed$lfq_wide
      values$dda_protein_meta <- parsed$protein_meta
      values$dda_db_engine    <- "Sage"

      # Parse report JSON if available
      for (rpt in c("results.json", "sage_report.json")) {
        local_rpt <- file.path(local_tmp, rpt)
        if (file.exists(local_rpt)) {
          values$dda_sage_report <- parse_sage_report(local_rpt)
          break
        }
      }

      # Compute QC metrics
      values$dda_qc_metrics <- compute_dda_qc_metrics(parsed$psms, parsed$lfq_wide)

      setProgress(0.8, detail = "Done!")
      values$dda_status <- "done"
      showNotification(
        paste("Sage results loaded:",
              nrow(parsed$lfq_wide), "proteins,",
              nrow(parsed$psms), "PSMs"),
        type = "message", duration = 10)
    })
  }

  # ============================================================================
  #    Run DDA pipeline (normalize + filter + impute + build EList)
  # ============================================================================
  observeEvent(input$run_dda_pipeline, {
    req(values$dda_lfq_wide, values$dda_protein_meta)

    # Check that group assignment exists
    if (is.null(values$metadata) || !"Group" %in% colnames(values$metadata)) {
      showNotification(
        "Please assign sample groups before running the pipeline.",
        type = "warning")
      return()
    }

    withProgress(message = "Running DDA pipeline...", value = 0.1, {
      result <- tryCatch(
        run_dda_pipeline(
          lfq_wide           = values$dda_lfq_wide,
          protein_meta       = values$dda_protein_meta,
          metadata_df        = values$metadata,
          norm_method        = input$dda_norm_method %||% "cyclicloess",
          min_valid_fraction = input$dda_min_valid %||% 0.5,
          impute_method      = input$dda_impute_method %||% "perseus",
          perseus_width      = input$dda_perseus_width %||% 0.3,
          perseus_shift      = input$dda_perseus_shift %||% 1.8
        ),
        error = function(e) {
          showNotification(paste("Pipeline error:", e$message), type = "error")
          NULL
        }
      )

      if (is.null(result)) return()

      setProgress(0.6, detail = "Building EList...")
      values$dda_elist                <- result$elist
      values$dda_n_proteins_prefilter <- result$n_prefilter
      values$dda_n_proteins_postfilter <- result$n_postfilter

      # Store as y_protein for downstream DE/QC/viz modules
      values$y_protein <- result$elist

      setProgress(0.7, detail = "Running limma DE...")

      # Run limma DE pipeline (same as DIA path in server_data.R)
      tryCatch({
        groups <- factor(values$metadata$Group)
        group_sizes <- table(groups)
        has_replicates <- all(group_sizes >= 2)

        if (has_replicates) {
          design <- model.matrix(~ 0 + groups)
          colnames(design) <- levels(groups)

          # Standard limma pipeline (not dpcDE -- MaxLFQ is already protein-level)
          fit <- limma::lmFit(result$elist, design)

          # Generate all pairwise contrasts
          combs <- combn(levels(groups), 2)
          forms <- apply(combs, 2, function(x) paste(x[2], "-", x[1]))
          contrast_matrix <- limma::makeContrasts(contrasts = forms, levels = design)
          fit <- limma::contrasts.fit(fit, contrast_matrix)
          fit <- limma::eBayes(fit)

          values$fit <- fit
          values$design <- design

          # Clear stale GSEA cache
          values$gsea_results_cache <- list()
          values$gsea_last_contrast <- NULL

          # Update all four comparison selectors
          updateSelectInput(session, "contrast_selector", choices = forms)
          updateSelectInput(session, "contrast_selector_signal", choices = forms, selected = forms[1])
          updateSelectInput(session, "contrast_selector_grid", choices = forms, selected = forms[1])
          updateSelectInput(session, "contrast_selector_pvalue", choices = forms, selected = forms[1])

          # Show DE Dashboard
          nav_show("main_tabs", "DE Dashboard")
          nav_show("main_tabs", "Gene Set Enrichment")
          nav_show("main_tabs", "AI Analysis")
          nav_show("main_tabs", "Output")
          nav_select("main_tabs", "DE Dashboard")

          n_de <- sum(limma::topTable(fit, coef = forms[1], number = Inf)$adj.P.Val < 0.05)
          showNotification(
            paste("DE analysis complete!", n_de, "significant proteins in", forms[1]),
            type = "message", duration = 10)
        } else {
          showNotification(
            "Some groups have <2 replicates. Skipping DE -- quantification-only mode.",
            type = "warning", duration = NULL)
        }
      }, error = function(e) {
        message("[DDA] limma DE failed: ", e$message)
        showNotification(paste("DE analysis failed:", e$message), type = "warning")
      })

      setProgress(1.0, detail = "Done!")
      showNotification(
        paste("DDA pipeline complete:",
              result$n_postfilter, "proteins after filtering"),
        type = "message", duration = 8)

      add_to_log("DDA Pipeline", c(
        sprintf("# Sage results: %d PSMs, %d proteins", nrow(values$dda_sage_psms), result$n_prefilter),
        sprintf("# Normalization: %s", input$dda_norm_method %||% "cyclicloess"),
        sprintf("# Imputation: %s", input$dda_impute_method %||% "perseus"),
        sprintf("# Valid value filter: %.0f%%", (input$dda_min_valid %||% 0.5) * 100),
        sprintf("# After filter: %d proteins", result$n_postfilter)
      ))
    })
  })

  # ============================================================================
  #    Annotate y_protein with Casanovo de novo confirmation columns
  #    Fires when both DDA pipeline and Casanovo classification are available
  # ============================================================================
  observe({
    req(values$y_protein, values$dda_casanovo_classification)

    cls <- values$dda_casanovo_classification
    prot_summary <- cls$protein_summary

    if (is.null(prot_summary) || nrow(prot_summary) == 0) return()

    genes_df <- values$y_protein$genes
    if ("DeNovo_Confirmed" %in% colnames(genes_df)) return()  # already annotated

    # Match protein IDs (Sage uses semicolon-separated protein groups)
    protein_ids <- genes_df$Protein.Group

    genes_df$DeNovo_Confirmed <- vapply(protein_ids, function(pid) {
      # Check if any protein in a semicolon-separated group has Casanovo confirmation
      ids <- trimws(strsplit(pid, ";")[[1]])
      match_idx <- which(prot_summary$proteins %in% ids)
      if (length(match_idx) > 0) sum(prot_summary$n_casanovo_confirmed[match_idx]) else 0L
    }, integer(1))

    genes_df$DeNovo_MaxScore <- vapply(protein_ids, function(pid) {
      ids <- trimws(strsplit(pid, ";")[[1]])
      match_idx <- which(prot_summary$proteins %in% ids)
      if (length(match_idx) > 0) max(prot_summary$casanovo_max_score[match_idx], na.rm = TRUE) else NA_real_
    }, numeric(1))

    genes_df$DeNovo_AvgAAScore <- vapply(protein_ids, function(pid) {
      ids <- trimws(strsplit(pid, ";")[[1]])
      match_idx <- which(prot_summary$proteins %in% ids)
      if (length(match_idx) > 0) {
        scores <- prot_summary$casanovo_mean_aa_score[match_idx]
        scores <- scores[!is.na(scores)]
        if (length(scores) > 0) mean(scores) else NA_real_
      } else NA_real_
    }, numeric(1))

    values$y_protein$genes <- genes_df
    message("[DDA] Added Casanovo annotation columns to y_protein$genes: ",
            sum(genes_df$DeNovo_Confirmed > 0), " proteins with de novo confirmation")
  })

  # ============================================================================
  #    Group assignment for DDA samples
  # ============================================================================
  output$dda_group_assignment_ui <- renderUI({
    req(values$dda_lfq_wide)
    samples <- colnames(values$dda_lfq_wide)
    if (length(samples) == 0) return(NULL)

    tags$div(
      style = "background: #f8f9fa; padding: 12px; border-radius: 8px; margin-bottom: 12px;",
      tags$h6(icon("users"), " Assign Sample Groups"),
      tags$p(style = "font-size: 12px; color: #6c757d;",
        "Assign each sample to a group for differential expression analysis."),
      lapply(seq_along(samples), function(i) {
        div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 4px;",
          tags$span(style = "font-size: 12px; min-width: 200px; font-family: monospace;",
            samples[i]),
          textInput(
            paste0("dda_group_", i),
            label = NULL,
            value = "",
            width = "150px",
            placeholder = "Group name"
          )
        )
      }),
      actionButton("dda_apply_groups", "Apply Groups",
        icon = icon("check"), class = "btn-primary btn-sm mt-2"),
      actionButton("run_dda_pipeline", "Run DE Pipeline",
        icon = icon("play"), class = "btn-success btn-sm mt-2 ms-2")
    )
  })

  # Apply group assignments
 observeEvent(input$dda_apply_groups, {
    req(values$dda_lfq_wide)
    samples <- colnames(values$dda_lfq_wide)
    groups <- vapply(seq_along(samples), function(i) {
      input[[paste0("dda_group_", i)]] %||% ""
    }, character(1))

    if (any(!nzchar(groups))) {
      showNotification("Please assign all samples to a group.", type = "warning")
      return()
    }

    values$metadata <- data.frame(
      SampleID = samples,
      File.Name = samples,
      Group = groups,
      stringsAsFactors = FALSE
    )
    showNotification("Groups assigned!", type = "message")
  })

  # ============================================================================
  #    Status UI
  # ============================================================================
  # Note: output$denovo_source_badge lives in server_denovo_controls.R —
  # see "Source engine + loaded-dataset banner" there. server_denovo_controls
  # is initialized after server_dda in app.R, so any renderer registered
  # here would be overwritten anyway.

  output$dda_job_status_ui <- renderUI({
    status <- values$dda_status %||% "idle"

    active_block <- switch(status,
      "idle" = NULL,
      "running" = {
        job_id <- values$dda_job_id %||% "?"
        div(
          class = "alert alert-info",
          style = "margin-top: 12px;",
          icon("spinner", class = "fa-spin"),
          paste(" Sage search running... Job ID:", job_id),
          tags$br(),
          tags$small("Polling every 15 seconds. You can navigate to other tabs.")
        )
      },
      "loading" = {
        div(
          class = "alert alert-info",
          style = "margin-top: 12px;",
          icon("spinner", class = "fa-spin"),
          " Loading Sage results from HPC..."
        )
      },
      "done" = {
        qc <- values$dda_qc_metrics
        div(
          class = "alert alert-success",
          style = "margin-top: 12px;",
          icon("check-circle"),
          " Sage search complete!",
          if (!is.null(qc)) {
            tags$div(
              style = "margin-top: 8px; font-size: 13px;",
              tags$strong(format(qc$n_psms, big.mark = ",")), " PSMs | ",
              tags$strong(format(qc$n_peptides, big.mark = ",")), " peptides | ",
              tags$strong(format(qc$n_proteins, big.mark = ",")), " proteins"
            )
          }
        )
      },
      "error" = {
        div(
          class = "alert alert-danger",
          style = "margin-top: 12px;",
          icon("exclamation-triangle"),
          " Sage search failed. Check SLURM logs for details.",
          if (!is.null(values$dda_output_dir)) {
            tags$small(style = "display: block; margin-top: 4px;",
              paste("Log dir:", file.path(values$dda_output_dir, "logs/")))
          }
        )
      }
    )

    # ── Persistent queue table — survives Shiny restarts via ~/.delimp_dda_queue.rds ──
    jobs <- values$dda_jobs %||% list()
    discover_btn <- tagList(
      actionButton("dda_discover_btn",
        tagList(icon("magnifying-glass"), " Discover from Hive"),
        class = "btn-outline-primary btn-sm",
        title = "Scan Hive for past DDA / Sage searches and add them to the queue"),
      if (length(jobs) > 0)
        actionButton("dda_clear_queue_btn",
          tagList(icon("broom"), " Clear queue"),
          class = "btn-outline-danger btn-sm",
          style = "margin-left: 4px;",
          title = "Empty the local queue (does not touch Hive results)")
      else NULL
    )

    queue_block <- if (length(jobs) > 0) {
      badge_color <- function(s) switch(s,
        "running" = "#f39c12",
        "loading" = "#3498db",
        "done"    = "#27ae60",
        "loaded"  = "#27ae60",
        "error"   = "#c0392b",
        "#6c757d")
      rows <- lapply(rev(jobs), function(j) {
        out_dir <- as.character(j$output_dir %||% "")
        action_btns <- if (nzchar(out_dir)) {
          dir_esc <- htmltools::htmlEscape(out_dir, attribute = TRUE)
          HTML(sprintf(
            '<div class="btn-group btn-group-sm" role="group">
               <button class="btn btn-outline-primary" title="Load this run" onclick="Shiny.setInputValue(\'dda_load_from_queue\', {dir: \'%s\', _ts: Date.now()}, {priority:\'event\'})"><i class="fa fa-download"></i></button>
               <button class="btn btn-outline-info" title="View search_info.md" onclick="Shiny.setInputValue(\'dda_view_info_from_queue\', {dir: \'%s\', _ts: Date.now()}, {priority:\'event\'})"><i class="fa fa-info-circle"></i></button>
               <button class="btn btn-outline-danger" title="Remove from local queue (does not touch Hive)" onclick="Shiny.setInputValue(\'dda_delete_from_queue\', {dir: \'%s\', _ts: Date.now()}, {priority:\'event\'})"><i class="fa fa-times"></i></button>
             </div>',
            dir_esc, dir_esc, dir_esc))
        } else NULL
        # Per-engine status badge: P/running/✓/✗/—
        engine_badge <- function(s) {
          s <- as.character(s %||% "none")
          spec <- switch(s,
            "done"     = list("✓",  "#27ae60", "completed"),
            "running"  = list("…",  "#3498db", "running"),
            "pending"  = list("P",  "#f39c12", "pending"),
            "failed"   = list("✗",  "#c0392b", "failed"),
            "cancelled"= list("⊘",  "#6c757d", "cancelled"),
            "none"     = list("—",  "#dee2e6", "not run"),
                         list(s,    "#6c757d", s))
          tags$span(
            style = sprintf("background:%s; color:white; padding:2px 8px; border-radius:4px; font-size: 0.78em; display:inline-block; min-width: 24px; text-align: center;",
                            spec[[2]]),
            title = spec[[3]],
            spec[[1]])
        }
        tags$tr(
          tags$td(j$name %||% "?"),
          tags$td(engine_badge(j$sage_status %||% j$status %||% "?")),
          tags$td(engine_badge(j$casanovo_status %||%
                                (if (isTRUE(j$casanovo)) "done" else "none"))),
          tags$td(engine_badge(j$blast_status %||% "none")),
          tags$td(tags$small(as.character(j$sage_job_id %||% j$job_id %||% "?"))),
          tags$td(tags$small(as.character(j$casanovo_job_id %||% "—"))),
          tags$td(tags$small(as.character(j$blast_job_id %||% "—"))),
          tags$td(j$n_files %||% "?"),
          tags$td(tags$small(j$submitted_at %||% "?")),
          tags$td(tags$small(style = "font-family: monospace;",
                              basename(out_dir))),
          tags$td(action_btns)
        )
      })
      div(style = "margin-top: 16px;",
        div(style = "display: flex; justify-content: space-between; align-items: center;",
          tags$h6(icon("list-check"), " DDA + Casanovo job queue ",
                  tags$small(style = "color: #6c757d; font-weight: normal;",
                             sprintf("(%d %s)", length(jobs),
                                     if (length(jobs) == 1) "entry" else "entries"))),
          discover_btn
        ),
        tags$table(class = "table table-sm",
          style = "font-size: 0.88em;",
          tags$thead(tags$tr(
            tags$th("Name"),
            tags$th(title = "Sage database search status", "Sage"),
            tags$th(title = "Casanovo de novo sequencing status", "Casanovo"),
            tags$th(title = "DIAMOND BLAST status (runs on Casanovo peptides)", "BLAST"),
            tags$th("Sage Job ID"),
            tags$th("Casanovo Job ID"),
            tags$th("BLAST Job ID"),
            tags$th("# Files"),
            tags$th("Submitted"),
            tags$th("Output dir"),
            tags$th("Actions"))),
          tags$tbody(rows)),
        tags$small(style = "color: #6c757d;",
          "Queue persists across app restarts at ", tags$code("~/.delimp_dda_queue.rds"))
      )
    } else {
      # Empty queue: still show the Discover button so users can pull in
      # historical searches submitted before persistence existed.
      div(style = "margin-top: 16px; padding: 12px; background: #f0f7ff; border: 1px solid #b8d4f0; border-radius: 6px;",
        div(style = "display: flex; justify-content: space-between; align-items: center;",
          tags$div(
            tags$strong("No DDA jobs in queue yet."),
            tags$br(),
            tags$small(style = "color: #6c757d;",
              "Submit a search above, or pull in previous runs from Hive.")
          ),
          discover_btn
        )
      )
    }

    tagList(active_block, queue_block)
  })

  # ── Discover-from-Hive handler ──────────────────────────────────────────
  # Scans common DDA output roots on Hive for past Sage searches and adds
  # them to the queue. Matches the proteog "Discover from Hive" pattern.
  # ── View search_info.md from a queue row ─────────────────────────────
  observeEvent(input$dda_view_info_from_queue, {
    payload <- input$dda_view_info_from_queue
    if (is.null(payload) || is.null(payload$dir)) return()
    sc <- dda_ssh_config()
    if (is.null(sc)) {
      showNotification("Connect to Hive first to view search info.",
                       type = "warning", duration = 6)
      return()
    }
    res <- tryCatch(
      ssh_exec(sc, sprintf("cat %s/search_info.md 2>/dev/null",
                            shQuote(payload$dir)),
               login_shell = FALSE, timeout = 15),
      error = function(e) NULL)
    txt <- if (!is.null(res) && identical(res$status, 0L))
             paste(res$stdout %||% character(), collapse = "\n")
           else ""
    if (!nzchar(txt)) {
      txt <- sprintf("(No search_info.md found at %s. Runs submitted before DE-LIMP v3.11.0+'s search_info.md was wired up won't have one.)",
                     payload$dir)
    }
    showModal(modalDialog(
      title = tagList(icon("info-circle"), " Search Info"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      tags$pre(style = "font-size: 0.85em; max-height: 60vh; overflow-y: auto; background: #f8f9fa; padding: 12px; border-radius: 6px;",
               txt)
    ))
  })

  # ── Delete a single queue entry ─────────────────────────────────────
  observeEvent(input$dda_delete_from_queue, {
    payload <- input$dda_delete_from_queue
    if (is.null(payload) || is.null(payload$dir)) return()
    jobs <- values$dda_jobs %||% list()
    target <- as.character(payload$dir)
    keep <- vapply(jobs, function(j) !identical(as.character(j$output_dir %||% ""), target),
                    logical(1))
    values$dda_jobs <- jobs[keep]
  })

  # ── Clear entire local queue (Hive results untouched) ─────────────────
  observeEvent(input$dda_clear_queue_btn, {
    showModal(modalDialog(
      title = "Clear DDA queue?",
      "This removes all entries from the local queue. The actual search results on Hive are NOT deleted.",
      footer = tagList(modalButton("Cancel"),
        actionButton("dda_clear_queue_confirm", "Clear queue",
                     class = "btn-danger", icon = icon("trash")))
    ))
  })
  observeEvent(input$dda_clear_queue_confirm, {
    values$dda_jobs <- list()
    removeModal()
    showNotification("DDA queue cleared.", type = "default", duration = 5)
  })

  # ── Browse for raw data directory on Hive (mirrors DIA Browse) ──────
  observeEvent(input$dda_browse_raw_btn, {
    # Reuse the existing SSH file browser by trigger-clicking it, with
    # a callback path set so the selection is written back to dda_raw_dir.
    sc <- dda_ssh_config()
    if (is.null(sc)) {
      showNotification("Connect to Hive first (in the DIA Search tab).",
                       type = "warning", duration = 6)
      return()
    }
    values$ssh_browse_target_input <- "dda_raw_dir"
    shinyjs::click("ssh_browse_raw_btn")
  })

  # Per-row "Load" button on queued/discovered jobs: open the Load Results
  # modal pre-populated with this run's output_dir, so the user can review
  # + click Confirm to run the existing scp_download + parse pipeline
  # (Sage TSV + LFQ + Casanovo .mztab + BLAST + search_info.md ingestion).
  observeEvent(input$dda_load_from_queue, {
    payload <- input$dda_load_from_queue
    if (is.null(payload) || is.null(payload$dir)) return()
    load_results_modal()
    # The modal builds its inputs on-render; update the path field once it's
    # in the DOM. session$onFlushed fires after the next UI flush.
    session$onFlushed(once = TRUE, function() {
      updateTextInput(session, "dda_load_path",
                      value = as.character(payload$dir))
    })
  })

  observeEvent(input$dda_discover_btn, {
    sc <- dda_ssh_config()
    if (is.null(sc)) {
      showNotification("Connect to Hive first to discover past searches.",
                       type = "warning", duration = 8)
      return()
    }
    # Look for results.sage.tsv under typical DDA roots. Use sort -u in case
    # multiple search roots resolve to the same tree (e.g. /brettsp/ ==
    # /$USER/), de-dupe BOTH the input queue and the find output. Filter
    # by mtime within the last 60 days so legacy test runs from many months
    # ago don't auto-pollute the queue (user can manually load those via
    # the Load Results modal if needed).
    user_root <- sprintf("/quobyte/proteomics-grp/de-limp/%s", sc$user %||% "")
    roots <- unique(c("/quobyte/proteomics-grp/de-limp/brettsp", user_root))
    roots <- roots[nzchar(roots)]
    cmd <- sprintf(
      "find %s -maxdepth 5 -name results.sage.tsv -type f -mtime -60 2>/dev/null | sort -u",
      paste(shQuote(roots), collapse = " ")
    )
    res <- tryCatch(
      ssh_exec(sc, cmd, login_shell = FALSE, timeout = 30),
      error = function(e) list(status = 1L, stdout = character(),
                                stderr = conditionMessage(e)))
    paths <- trimws(unlist(strsplit(paste(res$stdout %||% character(),
                                          collapse = "\n"), "\n")))
    paths <- unique(paths[!is.na(paths) & nzchar(paths)])

    # De-dupe the existing queue by output_dir (fix old duplicates from
    # an earlier double-scan).
    queue_in <- values$dda_jobs %||% list()
    if (length(queue_in) > 0) {
      dirs_seen <- character(0)
      deduped <- list()
      for (j in queue_in) {
        d <- as.character(j$output_dir %||% "")
        if (!nzchar(d) || d %in% dirs_seen) next
        dirs_seen <- c(dirs_seen, d)
        deduped[[length(deduped) + 1L]] <- j
      }
      if (length(deduped) != length(queue_in)) {
        message(sprintf("[DDA Discover] dedupe: %d → %d queue entries",
                         length(queue_in), length(deduped)))
      }
      values$dda_jobs <- deduped
    }
    if (length(paths) == 0) {
      showNotification("No past DDA searches found on Hive.",
                       type = "default", duration = 6)
      return()
    }
    current <- values$dda_jobs %||% list()
    existing_dirs <- vapply(current, function(j)
                            as.character(j$output_dir %||% ""),
                            character(1))
    added <- 0L
    for (p in paths) {
      output_dir <- dirname(p)
      if (output_dir %in% existing_dirs) next
      # Compute per-engine status from disk artifacts. The results.sage.tsv
      # presence drives Sage status (since `find` matched it, we know it
      # exists), .mztab presence drives Casanovo status. If search_info.md
      # mentions a Casanovo job ID but no mztabs are on disk, Casanovo is
      # in-flight (queued/running).
      cas <- tryCatch(ssh_exec(sc,
        sprintf("ls %s/casanovo/mztab/*.mztab 2>/dev/null | head -1",
                shQuote(output_dir)),
        login_shell = FALSE, timeout = 10),
        error = function(e) list(status = 1L, stdout = character()))
      has_mztab <- !is.null(cas) && identical(cas$status, 0L) &&
                    nzchar(trimws(paste(cas$stdout %||% character(),
                                         collapse = "")))
      casanovo_job_id <- NA_character_
      info_peek <- tryCatch(ssh_exec(sc,
        sprintf("grep -iE 'Casanovo (enabled|job ID)' %s/search_info.md 2>/dev/null",
                shQuote(output_dir)),
        login_shell = FALSE, timeout = 10),
        error = function(e) NULL)
      casanovo_announced <- FALSE
      if (!is.null(info_peek) && identical(info_peek$status, 0L)) {
        peek_txt <- paste(info_peek$stdout %||% character(), collapse = "\n")
        casanovo_announced <- grepl("enabled[^\\n]*TRUE|job ID[^\\n]*[0-9]+",
                                    peek_txt, ignore.case = TRUE)
        cas_id <- regmatches(peek_txt, regexpr("Casanovo job ID[^0-9]+[0-9]+",
                                                peek_txt, ignore.case = TRUE))
        if (length(cas_id)) {
          jd <- gsub("[^0-9]", "", cas_id)
          if (nzchar(jd)) casanovo_job_id <- jd
        }
      }
      sage_status <- "done"  # find matched results.sage.tsv → done
      casanovo_status <- if (has_mztab) "done"
                         else if (casanovo_announced) "pending"
                         else "none"
      has_casanovo <- has_mztab || casanovo_announced
      # Parse job_id / submitted_at / n_files from search_info.md OR (when
      # search_info.md is absent) fall back to file mtime + simple counts.
      # Verified search_info.md field formats on Hive:
      #   **Sage job ID**: 14532676
      #   **Submitted**: 2026-05-27 16:41:17
      #   **Raw files**: 3 files in /quobyte/...
      # We use a line-based parser (more reliable than regex on the full
      # blob — \\s inside an R character class doesn't mean whitespace).
      job_id <- NA_character_; submitted_at <- NA_character_
      n_files <- NA_integer_
      info <- tryCatch(ssh_exec(sc,
        sprintf("cat %s/search_info.md 2>/dev/null", shQuote(output_dir)),
        login_shell = FALSE, timeout = 10),
        error = function(e) NULL)
      info_lines <- if (!is.null(info) && identical(info$status, 0L))
                      unlist(strsplit(paste(info$stdout %||% character(),
                                            collapse = "\n"), "\n"))
                    else character()
      if (length(info_lines) > 0) {
        jid_line <- grep("Sage job ID", info_lines, value = TRUE)
        if (length(jid_line)) {
          jid_digits <- gsub("[^0-9]", "", jid_line[1])
          if (nzchar(jid_digits)) job_id <- jid_digits
        }
        ts_line <- grep("^\\s*\\*\\*Submitted\\*\\*", info_lines, value = TRUE)
        if (length(ts_line)) {
          ts_m <- regmatches(ts_line[1],
            regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}[ T][0-9]{2}:[0-9]{2}:[0-9]{2}",
                    ts_line[1]))
          if (length(ts_m) && nzchar(ts_m)) submitted_at <- ts_m
        }
        nf_line <- grep("Raw files", info_lines, value = TRUE)
        if (length(nf_line)) {
          nf_m <- regmatches(nf_line[1], regexpr("[0-9]+", nf_line[1]))
          if (length(nf_m) && nzchar(nf_m))
            n_files <- as.integer(nf_m)
        }
      }
      # Fallbacks when search_info.md is missing or didn't yield values —
      # use file system metadata on the results.sage.tsv itself.
      if (is.na(submitted_at) || !nzchar(submitted_at)) {
        st <- tryCatch(ssh_exec(sc,
          sprintf("stat -c %%y %s 2>/dev/null", shQuote(p)),
          login_shell = FALSE, timeout = 10),
          error = function(e) NULL)
        if (!is.null(st) && identical(st$status, 0L)) {
          stxt <- trimws(paste(st$stdout %||% character(), collapse = "\n"))
          # stat -c %y format: 2026-05-28 10:12:34.123456789 -0700
          ts_m <- regmatches(stxt, regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}", stxt))
          if (length(ts_m) && nzchar(ts_m)) submitted_at <- ts_m
        }
      }
      if (is.na(n_files) && length(info_lines) > 0) {
        # Count lines containing .raw or .mzML in the File List section
        n_guess <- sum(grepl("\\.raw|\\.mzML|\\.d\\b|\\.mgf", info_lines, ignore.case = TRUE))
        if (n_guess > 0) n_files <- as.integer(n_guess)
      }
      current[[length(current) + 1L]] <- list(
        job_id           = job_id,         # Sage job id (back-compat)
        sage_job_id      = job_id,
        casanovo_job_id  = casanovo_job_id,
        sage_status      = sage_status,
        casanovo_status  = casanovo_status,
        name             = basename(output_dir),
        output_dir       = output_dir,
        status           = "discovered",
        submitted_at     = submitted_at,
        n_files          = n_files,
        engine           = "Sage",
        casanovo         = has_casanovo
      )
      added <- added + 1L
    }
    values$dda_jobs <- current
    showNotification(sprintf("Discovered %d past DDA search(es). %d already in queue.",
                             added, length(paths) - added),
                     type = "default", duration = 8)
  })

  # ============================================================================
  #    Casanovo Status UI
  # ============================================================================
  output$dda_casanovo_status_ui <- renderUI({
    status <- values$dda_casanovo_status %||% "disabled"

    switch(status,
      "disabled" = NULL,
      "running" = {
        jid <- values$dda_casanovo_job_id %||% "?"
        convert_jid <- values$dda_casanovo_convert_job_id %||% "?"
        div(
          class = "alert alert-info",
          style = "margin-top: 8px; border-left: 4px solid #6f42c1;",
          icon("wand-magic-sparkles"),
          paste(" Casanovo de novo running... Array job:", jid),
          tags$br(),
          tags$small(paste("MGF convert:", convert_jid, "| GPU array:", jid))
        )
      },
      "loading" = {
        div(
          class = "alert alert-info",
          style = "margin-top: 8px; border-left: 4px solid #6f42c1;",
          icon("spinner", class = "fa-spin"),
          " Loading Casanovo results..."
        )
      },
      "done" = {
        cls <- values$dda_casanovo_classification
        n_psms <- if (!is.null(values$dda_casanovo_psms)) nrow(values$dda_casanovo_psms) else 0
        db_eng <- cls$db_engine %||% values$dda_db_engine %||% "Sage"
        db_badge_color <- if (db_eng == "DIA-NN") "#e74c3c" else "#2980b9"
        div(
          class = "alert alert-success",
          style = "margin-top: 8px; border-left: 4px solid #6f42c1;",
          icon("wand-magic-sparkles"),
          paste(" Casanovo complete!", format(n_psms, big.mark = ","), "de novo sequences"),
          tags$span(
            style = paste0("margin-left: 6px; padding: 2px 6px; border-radius: 3px; ",
              "font-size: 11px; color: white; background-color: ", db_badge_color, ";"),
            paste("vs", db_eng)
          ),
          if (!is.null(cls)) {
            tags$div(
              style = "margin-top: 4px; font-size: 12px;",
              tags$strong(cls$summary_stats$n_confirmed), " confirmed (",
              cls$summary_stats$pct_confirmed, "%) | ",
              tags$strong(cls$summary_stats$n_novel), " novel (",
              cls$summary_stats$pct_novel, "%)"
            )
          }
        )
      },
      "error" = {
        div(
          class = "alert alert-warning",
          style = "margin-top: 8px; border-left: 4px solid #6f42c1;",
          icon("exclamation-triangle"),
          " Casanovo failed. Sage results are unaffected."
        )
      }
    )
  })

  # ============================================================================
  #    Results summary UI (after pipeline)
  # ============================================================================
  output$dda_results_summary_ui <- renderUI({
    req(values$dda_elist)
    elist <- values$dda_elist
    n_pre  <- values$dda_n_proteins_prefilter %||% nrow(values$dda_lfq_wide)
    n_post <- values$dda_n_proteins_postfilter %||% nrow(elist$E)

    div(
      class = "alert alert-success",
      style = "margin-top: 12px;",
      icon("chart-bar"),
      tags$strong(" DDA pipeline complete"),
      tags$div(
        style = "margin-top: 8px; font-size: 13px;",
        tags$strong(format(n_post, big.mark = ",")), " proteins (from ",
        format(n_pre, big.mark = ","), " before filtering) | ",
        tags$strong(ncol(elist$E)), " samples",
        tags$br(),
        tags$small(
          "Normalization: ", values$dda_search_params$normalization %||% "cyclicloess",
          " | Imputation: ", values$dda_search_params$imputation %||% "perseus"
        )
      )
    )
  })

  # ============================================================================
  #    DDA QC Summary Card (rendered in QC tab or DDA panel)
  # ============================================================================
  output$dda_qc_summary_card <- renderUI({
    req(values$acquisition_mode == "dda")
    qc <- values$dda_qc_metrics
    if (is.null(qc)) return(NULL)

    div(
      style = paste(
        "background: white; border: 1px solid #dee2e6; border-radius: 8px;",
        "padding: 16px; margin-bottom: 16px;"
      ),
      tags$h6(icon("magnifying-glass"), " Sage DDA Search Summary",
        style = "margin-bottom: 12px; color: #2c3e50;"),
      div(
        class = "row text-center",
        div(class = "col-2",
          tags$h5(format(qc$n_psms, big.mark = ","), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "PSMs")
        ),
        div(class = "col-2",
          tags$h5(format(qc$n_peptides, big.mark = ","), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "Peptides")
        ),
        div(class = "col-2",
          tags$h5(format(qc$n_proteins, big.mark = ","), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "Proteins")
        ),
        div(class = "col-2",
          tags$h5(sprintf("%.1f", qc$med_pep_per_prot), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "Med. pep/prot")
        ),
        div(class = "col-2",
          tags$h5(paste0(qc$pct_missed_cleavage, "%"), style = "margin: 0; color: #2c3e50;"),
          tags$small(style = "color: #6c757d;", "Missed cleav.")
        ),
        div(class = "col-2",
          tags$h5(
            if (!is.na(qc$mass_error_ppm)) paste0(qc$mass_error_ppm, " ppm") else "N/A",
            style = "margin: 0; color: #2c3e50;"
          ),
          tags$small(style = "color: #6c757d;", "Mass error")
        )
      )
    )
  })

  # ============================================================================
  #    DDA De Novo Tab — Summary Cards, Tables, DIAMOND BLAST
  #    Renders in the De Novo > Casanovo nav_panel (ui.R)
  # ============================================================================

  # --- Summary cards: MOVED to server_denovo_controls.R (confidence-filtered) ---

  # --- Confirmed peptides table: MOVED to server_denovo_controls.R (confidence-filtered) ---

  # --- Novel peptides table: MOVED to server_denovo_controls.R (confidence-filtered) ---

  # ==========================================================================
  #    BLAST Results Visualization (comprehensive SwissProt BLAST view)
  # ==========================================================================

  # Helper: ensure species + category + contaminant_type columns exist on blast data
  blast_with_species <- reactive({
    # Depend on session trigger to force re-evaluation after session restore
    values$denovo_session_trigger
    blast_data <- values$denovo_blast %||% values$dda_casanovo_blast
    req(blast_data)
    blast <- blast_data
    req(nrow(blast) > 0)

    # Normalize identity column name (HPC load path uses "pident", DIAMOND path uses "identity")
    if ("identity" %in% names(blast) && !"pident" %in% names(blast)) {
      blast$pident <- blast$identity
    }

    # Fix legacy denovo_/casanovo_ IDs in peptide column — map back to sequences
    if (any(grepl("^(denovo|casanovo)_", blast$peptide))) {
      if (!is.null(values$dda_casanovo_classification)) {
        novel <- values$dda_casanovo_classification$novel
        if (!is.null(novel) && nrow(novel) > 0) {
          # Build map from denovo_N -> sequence (matching original FASTA order)
          novel_seqs <- unique(gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(novel$seq_stripped)))
          id_to_seq <- setNames(novel_seqs, paste0("denovo_", seq_along(novel_seqs)))
          # Also map casanovo_N format
          casanovo_map <- setNames(novel_seqs, paste0("casanovo_", seq_along(novel_seqs)))
          id_to_seq <- c(id_to_seq, casanovo_map)
          mapped <- id_to_seq[blast$peptide]
          blast$peptide[!is.na(mapped)] <- mapped[!is.na(mapped)]
        }
      }
      # Strip any remaining prefixes; fall back to query column
      blast$peptide <- sub("^(denovo|casanovo)_\\d+\\s*", "", blast$peptide)
      if ("query" %in% names(blast)) {
        empty <- !nzchar(blast$peptide)
        blast$peptide[empty] <- sub("^(denovo|casanovo)_\\d+\\s*", "", blast$query[empty])
      }
    }

    # --- Species attribution (single source of truth) ---
    # NCBI nr accessions (XP_/NP_/WP_/YP_…) carry NO `_SPECIES` mnemonic, so the
    # UniProt-style "text after the last underscore" parse produces garbage
    # (XP_025773238.1 -> "025773238.1"). When the taxonomy-aware per-peptide LCA
    # table is loaded it is the authoritative species/clade source — join it on
    # the canonical (mods-stripped, I/L-normalized) peptide key. UniProt searches
    # (no LCA table) fall back to the mnemonic parse. One definition of "species
    # for a de novo peptide" — CLAUDE.md architectural rule #3.
    lca_tbl <- values$dda_lca
    if (!is.null(lca_tbl) && nrow(lca_tbl) > 0 &&
        all(c("peptide", "lca_name") %in% names(lca_tbl))) {
      bkey <- gsub("I", "L", build_dda_canonical_peptide(blast$peptide))
      lkey <- gsub("I", "L", build_dda_canonical_peptide(lca_tbl$peptide))
      mi <- match(bkey, lkey)
      blast$species      <- lca_tbl$lca_name[mi]
      blast$lca_rank     <- if ("lca_rank"   %in% names(lca_tbl)) lca_tbl$lca_rank[mi]   else NA_character_
      blast$lca_category <- if ("category"   %in% names(lca_tbl)) lca_tbl$category[mi]   else NA_character_
      blast$diagnostic   <- if ("diagnostic" %in% names(lca_tbl)) lca_tbl$diagnostic[mi] else NA_integer_
      blast$species[is.na(blast$species)] <- "Unresolved"
    } else if (!"species" %in% names(blast)) {
      # UniProt mnemonic fallback: sp|ACC|PROT_SPECIES -> SPECIES
      blast$species <- sub(".*_", "", sub("^[a-z]+\\|[^|]+\\|", "", blast$subject))
    }
    if (!"category" %in% names(blast)) {
      blast$category <- ifelse(
        blast$pident >= 100, "Conserved",
        ifelse(blast$pident >= 90, "Near-match", "Distant")
      )
    }

    # --- Contaminant classification (paleoproteomics-aware) ---
    # Parse protein name from SwissProt ID: sp|ACC|PROT_SPECIES -> PROT
    blast$protein_name_raw <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", blast$subject))

    # Known human keratins (definite contaminant at 100% identity)
    human_keratin_pattern <- "^(KRT[0-9]|K1C[0-9]|K2C[0-9]|K22[EO]|KR[0-9])"
    # Common lab contaminant proteins
    lab_contam_pattern <- "TRYP_|TRYL_|ALBU_BOVIN|CAS[12]_BOVIN|CASA[12]_BOVIN|TRFE_BOVIN|ACTB_|TBB5_|ACTG_"
    # Common contaminant species
    contam_species <- c("HUMAN", "MOUSE", "BOVIN", "SHEEP", "RAT", "PIG")

    blast$contaminant_type <- vapply(seq_len(nrow(blast)), function(i) {
      sp <- blast$species[i]
      prot <- blast$protein_name_raw[i]
      pid <- blast$pident[i]

      is_keratin <- grepl(human_keratin_pattern, prot, ignore.case = TRUE)
      is_lab_contam <- grepl(lab_contam_pattern, blast$subject[i], ignore.case = TRUE)

      if (pid >= 99.5 && (is_keratin || is_lab_contam)) {
        # 100% identity to known contaminant proteins = definite contaminant
        return("Definite")
      }
      if (pid >= 85 && pid < 99.5 && sp == "HUMAN" && is_keratin) {
        # 85-95% to human keratin = likely avian keratin with no closer reference (NOT contaminant)
        return("Sample")
      }
      if (pid >= 99.5 && sp %in% contam_species) {
        # 100% identity to common contaminant species (non-keratin) = possible contaminant
        return("Possible")
      }
      "Sample"
    }, character(1))

    # Legacy flag for backward compat
    blast$is_contaminant <- blast$contaminant_type %in% c("Definite", "Possible")
    blast
  })

  # Filtered blast reactive — the DIAMOND BLAST tab shows NOVEL peptides only
  # (de novo BLAST is meaningful for peptides Sage's DB search missed) and
  # applies the user-controlled "Protein filter" dropdown.
  #
  # NOTE: we deliberately do NOT apply a keratin/name-based contaminant
  # heuristic here — for hair/feather/skin projects keratins are the SIGNAL.
  # The user opts into dropping them via the "Exclude skin & hair" mode.
  # Novel peptides above the current confidence threshold. The stored
  # classification is the FULL ~440k-PSM set; without gating here the DIAMOND
  # tab reports ~311k "novel" instead of the threshold-consistent count shown
  # everywhere else. Lowering the slider reveals more peptides (and more hits).
  novel_thresholded <- reactive({
    n <- values$dda_casanovo_classification$novel
    if (is.null(n) || nrow(n) == 0) return(n)
    thr <- input$dda_denovo_score_threshold %||% 0
    if ("score" %in% names(n)) n[n$score >= thr, , drop = FALSE] else n
  })

  blast_filtered <- reactive({
    blast <- blast_with_species()
    novel <- novel_thresholded()
    if (!is.null(novel) && nrow(novel) > 0 &&
        !is.null(blast) && nrow(blast) > 0 && "peptide" %in% names(blast)) {
      novel_norm <- unique(gsub("I", "L", novel$seq_stripped))
      blast_norm <- gsub("I", "L", build_dda_canonical_peptide(blast$peptide))
      blast <- blast[blast_norm %in% novel_norm, ]
    }
    dda_apply_protein_filter(blast, input$dda_protein_family_filter %||% "all")
  })

  # --- Species (LCA) panel: host/microbiome/conserved from the nr LCA table ---
  output$denovo_lca_category <- plotly::renderPlotly({
    lca <- lca_filtered()
    req(lca, nrow(lca) > 0, "category" %in% names(lca))
    tab <- as.data.frame(table(lca$category), stringsAsFactors = FALSE)
    names(tab) <- c("category", "n")
    tab <- tab[order(tab$n), ]
    plotly::plot_ly(tab, x = ~n, y = ~factor(category, levels = category),
                    type = "bar", orientation = "h",
                    marker = list(color = "#43a047")) %>%
      plotly::layout(
        title = list(text = paste0("Peptides by category (n=", nrow(lca), ")"), font = list(size = 13)),
        xaxis = list(title = "peptides"), yaxis = list(title = ""), margin = list(l = 130))
  })

  output$denovo_lca_top_species <- plotly::renderPlotly({
    lca <- lca_filtered()
    req(lca, nrow(lca) > 0, "lca_name" %in% names(lca))
    d <- lca[lca$category %in% "host", , drop = FALSE]
    req(nrow(d) > 0)
    cnt <- sort(table(d$lca_name), decreasing = TRUE)
    cnt <- cnt[seq_len(min(15, length(cnt)))]
    tab <- data.frame(taxon = names(cnt), n = as.integer(cnt), stringsAsFactors = FALSE)
    tab <- tab[order(tab$n), ]
    plotly::plot_ly(tab, x = ~n, y = ~factor(taxon, levels = taxon),
                    type = "bar", orientation = "h",
                    marker = list(color = "#1565c0")) %>%
      plotly::layout(
        title = list(text = "Top host taxa (LCA)", font = list(size = 13)),
        xaxis = list(title = "peptides"), yaxis = list(title = ""), margin = list(l = 170))
  })

  output$denovo_lca_table <- DT::renderDT({
    lca <- lca_filtered()
    req(lca, nrow(lca) > 0)
    ord <- if ("top_pident" %in% names(lca)) order(-lca$top_pident) else seq_len(nrow(lca))
    d <- lca[ord, , drop = FALSE]
    if (all(c("lca_name", "lca_taxid") %in% names(d)))
      d$lca_name <- ncbi_tax_link(d$lca_name, d$lca_taxid)
    DT::datatable(d, rownames = FALSE, filter = "top", escape = FALSE,
                  options = list(pageLength = 25, scrollX = TRUE))
  })

  observeEvent(input$denovo_lca_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " LCA species attribution"),
      tags$p("Each de novo peptide's nr BLAST hits are mapped to NCBI taxa; the lowest common ancestor (LCA) of the top hits is the most specific taxon consistent with all of them."),
      tags$ul(
        tags$li(tags$b("Diagnostic"), " — LCA resolves to species/genus (e.g. Panthera onca, Leopardus): real species signal."),
        tags$li(tags$b("Conserved"), " — LCA is family or higher (e.g. Felidae, Eutheria): peptide shared across many taxa, not species-diagnostic."),
        tags$li(tags$b("Microbiome"), " — LCA falls in Bacteria / Archaea / Viruses.")),
      tags$p("For a non-model organism not in nr (e.g. ocelot), peptides correctly converge on the family (Felidae) and nearest relatives rather than over-committing to one species."),
      easyClose = TRUE, size = "l", footer = modalButton("Close")))
  })

  # ════════════════════════════════════════════════════════════════════════
  # MASTER PEPTIDE TABLE  (v3.11.22)
  # One row per de novo peptide, joining the three evidence streams on a single
  # canonical key (mods stripped, I/L-normalized): Casanovo confidence +
  # confirmed/novel, Sage protein hit, and the nr LCA taxon attribution. The
  # biologist-facing "what is this sample" view and the single source the
  # confidence slider, the LCA panel, and the export all read from. Built data-
  # first — nothing here is dataset-specific.
  # ════════════════════════════════════════════════════════════════════════

  # Per-canonical-peptide max Casanovo confidence. The nr BLAST/LCA ran on ALL
  # de novo peptides with no score cutoff, so the slider is applied here.
  casanovo_pep_score <- reactive({
    cls <- values$dda_casanovo_classification
    if (is.null(cls) || is.null(cls$classified) || nrow(cls$classified) == 0)
      return(NULL)
    cas <- data.table::as.data.table(cls$classified)
    if (!"seq_norm" %in% names(cas)) {
      base_seq <- if ("seq_stripped" %in% names(cas)) cas$seq_stripped else cas$sequence
      cas$seq_norm <- gsub("I", "L", build_dda_canonical_peptide(base_seq))
    }
    if (!"score" %in% names(cas)) cas$score <- NA_real_
    cas[, list(best_score = suppressWarnings(max(score, na.rm = TRUE))),
        by = "seq_norm"]
  })

  # LCA gated by the confidence slider; peptides below threshold (by best
  # Casanovo PSM) are hidden, no-score peptides keep showing. Internal var
  # `lcat` so the global `lca <- ...` rename doesn't recurse.
  lca_filtered <- reactive({
    lcat <- values$dda_lca
    req(lcat, nrow(lcat) > 0)
    sc <- casanovo_pep_score()
    thr <- input$dda_denovo_score_threshold %||% 0
    if (!is.null(sc) && "peptide" %in% names(lcat)) {
      lkey <- gsub("I", "L", build_dda_canonical_peptide(lcat$peptide))
      sco  <- sc$best_score[match(lkey, sc$seq_norm)]
      lcat <- lcat[is.na(sco) | sco >= thr, , drop = FALSE]
    }
    lcat
  })

  # The master join (unfiltered — confidence filter applied downstream).
  denovo_master <- reactive({
    values$denovo_session_trigger
    cls <- values$dda_casanovo_classification
    req(cls, !is.null(cls$classified), nrow(cls$classified) > 0)
    build_denovo_master(cls$classified, values$dda_sage_psms, values$dda_lca)
  })

  # Master gated by the confidence slider (hide low-confidence by default).
  denovo_master_filtered <- reactive({
    m <- denovo_master()
    req(nrow(m) > 0)
    thr <- input$dda_denovo_score_threshold %||% 0
    if ("Casanovo_score" %in% names(m) && any(is.finite(m$Casanovo_score)))
      m <- m[is.na(m$Casanovo_score) | m$Casanovo_score >= thr, , drop = FALSE]
    m
  })

  output$denovo_master_table <- DT::renderDT({
    m <- denovo_master_filtered()
    req(nrow(m) > 0)
    # Linkify the taxon to its NCBI Taxonomy page (new tab), then drop the raw taxid.
    if (all(c("Species_or_clade", "lca_taxid") %in% names(m)))
      m$Species_or_clade <- ncbi_tax_link(m$Species_or_clade, m$lca_taxid)
    show <- m[, setdiff(names(m), c("seq_norm", "lca_taxid")), drop = FALSE]
    pref <- intersect(c("Peptide", "Casanovo_score", "Found_by_Sage",
                        "Species_or_clade", "Rank", "Type", "Best_pct_ID",
                        "Diagnostic", "Sage_protein", "n_PSMs"), names(show))
    show <- show[, c(pref, setdiff(names(show), pref)), drop = FALSE]
    if ("Best_pct_ID" %in% names(show))
      show <- show[order(-replace(show$Best_pct_ID, is.na(show$Best_pct_ID), -1)), ]
    dt <- DT::datatable(show, rownames = FALSE, filter = "top", escape = FALSE,
                        options = list(pageLength = 25, scrollX = TRUE))
    if ("Type" %in% names(show))
      dt <- DT::formatStyle(dt, "Type", backgroundColor = DT::styleEqual(
        c("host", "microbiome", "other/conserved", "plant/fungal"),
        c("#e8f5e9", "#fff3e0", "#eceff1", "#f3e5f5")))
    dt
  })

  # Plain-language species verdict — computed from the data, never hardcoded.
  # Injected as static HTML (shinyjs::html) into a plain div because renderUI is
  # unreliable inside navset_card_tab sub-tabs (CLAUDE.md UI gotcha).
  observe({
    m <- denovo_master_filtered()
    req(nrow(m) > 0)
    if (!all(c("Type", "Species_or_clade") %in% names(m))) {
      shinyjs::html("denovo_master_verdict_box", as.character(
        div(class = "alert alert-info",
            "Run the nr BLAST + LCA to populate a species attribution here.")))
      return(invisible(NULL))
    }
    host    <- m[!is.na(m$Type) & m$Type == "host", , drop = FALSE]
    n_micro <- sum(m$Type == "microbiome", na.rm = TRUE)
    n_cons  <- sum(m$Type %in% "other/conserved", na.rm = TRUE)
    diag    <- if ("Diagnostic" %in% names(host))
                 host[!is.na(host$Diagnostic) & host$Diagnostic, , drop = FALSE]
               else host[0, , drop = FALSE]
    top_diag <- head(sort(table(diag$Species_or_clade), decreasing = TRUE), 5)
    top_host <- head(sort(table(host$Species_or_clade), decreasing = TRUE), 5)
    fmt <- function(tb) if (length(tb))
      paste(sprintf("%s (%d)", names(tb), as.integer(tb)), collapse = ", ") else "—"
    shinyjs::html("denovo_master_verdict_box", as.character(
      div(class = "alert alert-success",
        tags$h5(icon("dna"), " Species attribution (confidence-filtered)"),
        tags$p(sprintf("%s peptides placed — %d host, %d conserved, %d microbiome.",
                       format(nrow(m), big.mark = ","), nrow(host), n_cons, n_micro)),
        tags$p(tags$b("Diagnostic (species/genus) taxa: "), fmt(top_diag)),
        tags$p(tags$b("All host taxa (incl. conserved clades): "), fmt(top_host)),
        tags$small(class = "text-muted",
          "Conserved peptides resolve only to family or higher and are NOT pinned to one species. For a non-model organism absent from NCBI, peptides correctly converge on the nearest represented clade."))))
  })

  output$dda_master_download <- downloadHandler(
    filename = function() paste0("denovo_master_table_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      m <- denovo_master_filtered()
      utils::write.csv(m[, setdiff(names(m), "seq_norm"), drop = FALSE],
                       file, row.names = FALSE)
    }
  )

  # Venn-style overlap stats (Sage vs Casanovo), confidence-filtered. Computed
  # from the master join + the full Sage peptide set on the shared canonical key.
  observe({
    m <- denovo_master_filtered()
    req(nrow(m) > 0, "seq_norm" %in% names(m))
    cas_keys <- unique(m$seq_norm)
    cas_n    <- length(cas_keys)
    sage <- values$dda_sage_psms
    sage_keys <- if (!is.null(sage) && "peptide" %in% names(sage))
      unique(gsub("I", "L", build_dda_canonical_peptide(sage$peptide))) else character(0)
    sage_n      <- length(sage_keys)
    overlap     <- length(intersect(cas_keys, sage_keys))
    denovo_only <- cas_n - overlap
    sage_only   <- max(sage_n - overlap, 0L)
    with_sp <- if ("Species_or_clade" %in% names(m)) sum(!is.na(m$Species_or_clade)) else 0L
    diagn   <- if ("Diagnostic" %in% names(m)) sum(m$Diagnostic %in% TRUE, na.rm = TRUE) else 0L
    nf <- function(x) format(as.integer(x), big.mark = ",")
    card <- function(val, lab, col, bg) sprintf(
      '<div style="flex:1;min-width:120px;background:%s;padding:12px;border-radius:8px;text-align:center;"><div style="font-size:23px;font-weight:700;color:%s;">%s</div><div style="font-size:12px;color:#555;">%s</div></div>',
      bg, col, val, lab)
    cards <- paste0(
      card(nf(sage_n),      "Sage peptides",        "#1565c0", "#e3f2fd"),
      card(nf(cas_n),       "Casanovo peptides",    "#6a1b9a", "#f3e5f5"),
      card(nf(overlap),     "Found by both",        "#2e7d32", "#e8f5e9"),
      card(nf(denovo_only), "De-novo only",         "#e65100", "#fff3e0"),
      card(nf(with_sp),     "With nr species",      "#00695c", "#e0f2f1"),
      card(nf(diagn),       "Diagnostic (sp/genus)","#ad1457", "#fce4ec"))
    shinyjs::html("denovo_master_stats_box", paste0(
      '<div style="display:flex;gap:8px;flex-wrap:wrap;">', cards, '</div>',
      '<div style="font-size:11px;color:#888;margin-top:4px;">',
      'Sage-only peptides: ', nf(sage_only),
      ' &middot; counts reflect the confidence slider (low-confidence de novo hidden).</div>'))
  })

  output$dda_lca_download <- downloadHandler(
    filename = function() paste0("peptide_lca_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv"),
    content = function(file) {
      lcat <- lca_filtered()
      utils::write.csv(as.data.frame(lcat), file, row.names = FALSE)
    }
  )

  # --- Summary cards ---
  output$dda_blast_summary_cards <- renderUI({
    blast_all <- blast_with_species()
    blast <- blast_filtered()
    req(nrow(blast_all) > 0)
    novel <- novel_thresholded()
    n_novel <- length(unique(novel$seq_stripped))
    n_with_hits <- length(unique(blast$peptide))
    n_no_hits <- n_novel - n_with_hits
    pct_hits <- round(100 * n_with_hits / max(n_novel, 1), 1)

    # Top species by best-hit count (deduplicate to best hit per peptide)
    best_hits <- blast[!duplicated(blast$peptide), ]
    top_sp <- names(sort(table(best_hits$species), decreasing = TRUE))[1]
    mean_id <- round(mean(blast$pident, na.rm = TRUE), 1)

    # Contaminant stats
    n_definite <- sum(blast_all$contaminant_type == "Definite")
    n_possible <- sum(blast_all$contaminant_type == "Possible")
    contam_text <- paste0(n_definite, " definite")
    if (n_possible > 0) contam_text <- paste0(contam_text, " + ", n_possible, " possible")

    tagList(
      div(class = "row", style = "margin-bottom: 15px;",
        div(class = "col-md-2",
          div(style = "background: #e8f5e9; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #2e7d32;", n_novel),
            tags$small("Novel peptides")
          )
        ),
        div(class = "col-md-2",
          div(style = "background: #e3f2fd; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #1565c0;", paste0(n_with_hits, " (", pct_hits, "%)")),
            tags$small("With hits")
          )
        ),
        div(class = "col-md-2",
          div(style = "background: #fff3e0; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #e65100;", n_no_hits),
            tags$small("No hits (truly novel)")
          )
        ),
        div(class = "col-md-2",
          div(style = "background: #f3e5f5; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #7b1fa2;", top_sp %||% "N/A"),
            tags$small("Top species")
          )
        ),
        div(class = "col-md-2",
          div(style = "background: #fce4ec; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #c62828;", paste0(mean_id, "%")),
            tags$small("Mean identity")
          )
        ),
        div(class = "col-md-2",
          div(style = "background: #fff8e1; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #f57f17;", contam_text),
            tags$small("Contaminants")
          )
        )
      )
    )
  })

  # --- Taxonomic breakdown: donut chart ---
  output$dda_blast_species_donut <- plotly::renderPlotly({
    blast <- blast_filtered()
    req(nrow(blast) > 0)
    # Best hit per peptide for species assignment
    best_hits <- blast[order(blast$pident, decreasing = TRUE), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide), ]

    sp_counts <- sort(table(best_hits$species), decreasing = TRUE)
    top_n <- min(10, length(sp_counts))
    top_sp <- sp_counts[seq_len(top_n)]
    if (length(sp_counts) > top_n) {
      top_sp <- c(top_sp, "Other" = sum(sp_counts[(top_n + 1):length(sp_counts)]))
    }

    plotly::plot_ly(
      labels = names(top_sp),
      values = as.numeric(top_sp),
      type = "pie",
      hole = 0.4,
      textinfo = "label+percent",
      textposition = "auto",
      marker = list(colors = grDevices::rainbow(length(top_sp), s = 0.6, v = 0.9))
    ) %>%
      plotly::layout(
        title = list(text = "Species Distribution (best hit per peptide)", font = list(size = 14)),
        showlegend = TRUE,
        legend = list(orientation = "v", x = 1.02, y = 0.5)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Taxonomic breakdown: bar chart ---
  output$dda_blast_species_bar <- plotly::renderPlotly({
    blast <- blast_filtered()
    req(nrow(blast) > 0)
    best_hits <- blast[order(blast$pident, decreasing = TRUE), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide), ]

    sp_counts <- sort(table(best_hits$species), decreasing = TRUE)
    top_n <- min(15, length(sp_counts))
    top_sp <- sp_counts[seq_len(top_n)]

    sp_df <- data.frame(
      species = factor(names(top_sp), levels = rev(names(top_sp))),
      count   = as.numeric(top_sp),
      stringsAsFactors = FALSE
    )

    plotly::plot_ly(
      data = sp_df,
      x = ~count, y = ~species,
      type = "bar", orientation = "h",
      marker = list(color = "#5c6bc0")
    ) %>%
      plotly::layout(
        title = list(text = "Species Hit Counts", font = list(size = 14)),
        xaxis = list(title = "Number of novel peptides"),
        yaxis = list(title = ""),
        margin = list(l = 120)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Species summary text ---
  output$dda_blast_species_summary <- renderUI({
    blast <- blast_filtered()
    req(nrow(blast) > 0)
    best_hits <- blast[order(blast$pident, decreasing = TRUE), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide), ]

    sp_counts <- sort(table(best_hits$species), decreasing = TRUE)
    total <- sum(sp_counts)
    top_lines <- vapply(seq_len(min(5, length(sp_counts))), function(i) {
      pct <- round(100 * sp_counts[i] / total, 1)
      paste0(pct, "% ", names(sp_counts)[i])
    }, character(1))

    tags$p(style = "color: #555; font-size: 0.95em; margin-top: 10px;",
      paste("Species breakdown:", paste(top_lines, collapse = ", "),
        if (length(sp_counts) > 5) paste0(", + ", length(sp_counts) - 5, " more") else "")
    )
  })

  # --- Identity distribution histogram (colored by species) ---
  output$dda_blast_identity_hist <- plotly::renderPlotly({
    blast <- blast_filtered()
    req(nrow(blast) > 0)
    # Best hit per peptide
    best_hits <- blast[order(blast$pident, decreasing = TRUE), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide), ]

    # Top 5 species, rest as "Other"
    sp_counts <- sort(table(best_hits$species), decreasing = TRUE)
    top5 <- names(sp_counts)[seq_len(min(5, length(sp_counts)))]
    best_hits$sp_group <- ifelse(best_hits$species %in% top5, best_hits$species, "Other")
    best_hits$sp_group <- factor(best_hits$sp_group, levels = c(top5, "Other"))

    colors <- c(
      grDevices::rainbow(length(top5), s = 0.5, v = 0.85),
      "#cccccc"
    )
    names(colors) <- c(top5, "Other")

    p <- ggplot2::ggplot(best_hits, ggplot2::aes(x = pident, fill = sp_group)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.85, position = "stack") +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::geom_vline(xintercept = 90, linetype = "dashed", color = "#e65100", alpha = 0.7) +
      ggplot2::geom_vline(xintercept = 100, linetype = "dashed", color = "#2e7d32", alpha = 0.7) +
      ggplot2::annotate("text", x = 91, y = Inf, label = "Near-match", vjust = 1.5,
        hjust = 0, color = "#e65100", size = 3) +
      ggplot2::annotate("text", x = 99, y = Inf, label = "Identical", vjust = 1.5,
        hjust = 1, color = "#2e7d32", size = 3) +
      ggplot2::labs(
        x = "% Identity to SwissProt",
        y = "Number of peptides",
        fill = "Species",
        subtitle = "100% = identical to reference; 90-99% = likely variants; <80% = distant homologs"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")

    plotly::ggplotly(p) %>%
      plotly::layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.05)) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Top proteins by peptide count ---
  output$dda_blast_top_proteins <- plotly::renderPlotly({
    blast <- tryCatch(blast_filtered(), error = function(e) NULL)
    req(!is.null(blast), nrow(blast) > 1)

    # Protein label: UniProt mnemonic, or FULL accession for nr (single helper).
    blast$protein_name <- dda_protein_label(blast$subject)
    blast$protein_label <- paste0(blast$protein_name, " (", blast$species, ")")

    # Count unique peptides per protein (best hit per peptide)
    best_hits <- blast[order(-blast$bitscore), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide), ]

    prot_counts <- as.data.frame(table(best_hits$protein_label), stringsAsFactors = FALSE)
    colnames(prot_counts) <- c("Protein", "Peptides")
    prot_counts <- prot_counts[order(-prot_counts$Peptides), ]
    top_n <- min(50, nrow(prot_counts))
    prot_counts <- prot_counts[seq_len(top_n), ]
    prot_counts$Protein <- factor(prot_counts$Protein, levels = rev(prot_counts$Protein))

    # Color by contaminant status
    contam_prots <- grepl("KRT\\d|K1C\\d|K2C\\d|TRYP_|ALBU_BOVIN|CASA|_HUMAN|_MOUSE|_BOVIN|_SHEEP|_RAT",
      prot_counts$Protein)
    prot_counts$Type <- ifelse(contam_prots, "Contaminant", "Sample")

    plotly::plot_ly(prot_counts,
      y = ~Protein, x = ~Peptides, color = ~Type,
      colors = c("Sample" = "#2196F3", "Contaminant" = "#FF9800"),
      type = "bar", orientation = "h",
      hoverinfo = "text",
      text = ~paste(Protein, "<br>", Peptides, "peptides")
    ) %>%
      plotly::layout(
        title = list(text = "Top Proteins by De Novo Peptide Count", font = list(size = 14)),
        xaxis = list(title = "Unique Peptides"),
        yaxis = list(title = "", tickfont = list(size = 9)),
        showlegend = TRUE,
        legend = list(x = 0.7, y = 0.1),
        margin = list(l = 200),
        height = max(400, top_n * 18)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Species common name mapping ---
  species_common_names <- c(
    CHICK = "Chicken", HUMAN = "Human", MOUSE = "Mouse", BOVIN = "Bovine",
    COLLI = "Columba livia (Pigeon)", ANAPL = "Anas platyrhynchos (Mallard)",
    CATAU = "Cathartes aura (Vulture)", MYCAM = "Mycteria americana (Stork)",
    SHEEP = "Sheep", RAT = "Rat", PIG = "Pig", HORSE = "Horse",
    RABIT = "Rabbit", CANFA = "Dog", FELCA = "Cat", DANRE = "Zebrafish",
    DROME = "Fruit fly", CAEEL = "C. elegans", YEAST = "Yeast",
    ECOLI = "E. coli", ARATH = "Arabidopsis"
  )

  get_common_name <- function(sp) {
    cn <- species_common_names[sp]
    ifelse(is.na(cn), sp, cn)
  }

  # --- Compute species resolution data (shared by bar chart + diagnostic table) ---
  blast_species_resolution <- reactive({
    blast <- tryCatch(blast_filtered(), error = function(e) NULL)
    req(!is.null(blast), nrow(blast) > 1)

    # For each peptide, rank hits by bitscore
    blast <- blast[order(blast$peptide, -blast$bitscore), ]

    peptides <- unique(blast$peptide)
    resolution <- do.call(rbind, lapply(peptides, function(pep) {
      hits <- blast[blast$peptide == pep, ]
      if (nrow(hits) == 0) return(NULL)

      best <- hits[1, ]
      best_sp <- best$species
      best_id <- best$pident
      best_prot <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", best$subject))
      best_acc <- sub("^[a-z]+\\|([^|]+)\\|.*", "\\1", best$subject)

      # Find second-best hit from a DIFFERENT species
      other_sp <- hits[hits$species != best_sp, ]
      if (nrow(other_sp) > 0) {
        second <- other_sp[1, ]
        second_sp <- second$species
        second_id <- second$pident
      } else {
        second_sp <- NA_character_
        second_id <- NA_real_
      }

      delta <- if (!is.na(second_id)) best_id - second_id else best_id

      data.frame(
        peptide = pep,
        best_species = best_sp,
        best_identity = best_id,
        second_species = second_sp,
        second_identity = second_id,
        delta_identity = delta,
        protein_name = best_prot,
        accession = best_acc,
        stringsAsFactors = FALSE
      )
    }))
    resolution
  })

  # --- Feature 2: Top Diagnostic Peptides Summary Card (shown first) ---
  output$dda_blast_diagnostic_card <- renderUI({
    res <- tryCatch(blast_species_resolution(), error = function(e) NULL)
    req(!is.null(res), nrow(res) > 0)

    # Top peptides by delta_identity
    res <- res[order(-res$delta_identity), ]
    top_n <- min(10, nrow(res))
    top <- res[seq_len(top_n), ]

    # Count diagnostic peptides (delta > 15%)
    n_diagnostic <- sum(res$delta_identity > 15, na.rm = TRUE)
    diagnostic_species <- if (n_diagnostic > 0) {
      diag_rows <- res[res$delta_identity > 15, ]
      sp_tab <- sort(table(diag_rows$best_species), decreasing = TRUE)
      paste(paste0(get_common_name(names(sp_tab)), " (", sp_tab, ")"), collapse = ", ")
    } else {
      "none"
    }

    # Identify dominant protein families among diagnostic peptides
    protein_summary <- ""
    if (n_diagnostic > 0) {
      diag_rows <- res[res$delta_identity > 15, ]
      prot_tab <- sort(table(diag_rows$protein_name), decreasing = TRUE)
      top_prots <- names(prot_tab)[seq_len(min(3, length(prot_tab)))]
      protein_summary <- paste(top_prots, collapse = ", ")
    }

    # Build table rows
    table_rows <- lapply(seq_len(nrow(top)), function(i) {
      r <- top[i, ]
      pep_display <- if (nchar(r$peptide) > 20) paste0(substr(r$peptide, 1, 17), "...") else r$peptide
      second_sp <- if (is.na(r$second_species)) "-" else get_common_name(r$second_species)
      second_id <- if (is.na(r$second_identity)) "-" else paste0(round(r$second_identity, 1), "%")
      delta_color <- if (r$delta_identity > 15) "#2e7d32" else if (r$delta_identity > 5) "#f57c00" else "#757575"
      uniprot_link <- paste0('<a href="https://www.uniprot.org/uniprot/', r$accession,
        '" target="_blank">', r$protein_name, '</a>')

      tags$tr(
        tags$td(style = "font-family: monospace; font-size: 0.85em;", pep_display),
        tags$td(get_common_name(r$best_species)),
        tags$td(paste0(round(r$best_identity, 1), "%")),
        tags$td(second_sp),
        tags$td(second_id),
        tags$td(style = paste0("font-weight: bold; color: ", delta_color, ";"),
          paste0(round(r$delta_identity, 1), "%")),
        tags$td(HTML(uniprot_link))
      )
    })

    # Summary sentence
    summary_text <- if (n_diagnostic > 0) {
      base <- paste0(n_diagnostic, " peptide", if (n_diagnostic > 1) "s" else "",
        " are species-diagnostic (delta > 15%): ", diagnostic_species, ".")
      if (nzchar(protein_summary)) {
        paste0(base, " Primarily from ", protein_summary, " proteins.")
      } else {
        base
      }
    } else {
      "No peptides exceed the 15% delta threshold for species-diagnostic classification."
    }

    tagList(
      div(style = "background: linear-gradient(135deg, #e3f2fd, #f3e5f5); border-radius: 8px; padding: 16px; margin-bottom: 16px; border: 1px solid #90caf9;",
        tags$h5(icon("crosshairs"), " Top Diagnostic Peptides",
          style = "margin-top: 0; margin-bottom: 8px; color: #1565c0;"),
        tags$p(style = "color: #37474f; font-size: 0.92em; margin-bottom: 12px;",
          "Peptides with the largest identity gap between their best and second-best species hit. ",
          "High delta values indicate species-specific sequences useful for taxonomic identification."),
        div(style = "overflow-x: auto;",
          tags$table(class = "table table-sm table-hover", style = "font-size: 0.88em; margin-bottom: 8px;",
            tags$thead(tags$tr(
              tags$th("Sequence"), tags$th("Best Species"), tags$th("% Identity"),
              tags$th("2nd Species"), tags$th("2nd %"), tags$th("Delta"),
              tags$th("Source Protein")
            )),
            tags$tbody(table_rows)
          )
        ),
        tags$p(style = "margin-bottom: 0; font-size: 0.9em; color: #1b5e20; font-style: italic;",
          icon("info-circle"), " ", summary_text)
      )
    )
  })

  # --- Feature 1: Species Resolution Bar Chart ---
  output$dda_blast_species_resolution <- plotly::renderPlotly({
    res <- tryCatch(blast_species_resolution(), error = function(e) NULL)
    req(!is.null(res), nrow(res) > 0)

    # Top 30 by delta_identity
    res <- res[order(-res$delta_identity), ]
    top_n <- min(30, nrow(res))
    res <- res[seq_len(top_n), ]

    # Truncate peptide labels
    res$pep_label <- ifelse(nchar(res$peptide) > 18,
      paste0(substr(res$peptide, 1, 15), "..."), res$peptide)
    # Make labels unique for factor ordering
    res$pep_label <- make.unique(res$pep_label, sep = " ")

    # Sort ascending (bottom to top in horizontal bar)
    res <- res[order(res$delta_identity), ]
    res$pep_label <- factor(res$pep_label, levels = res$pep_label)

    # Color by best species
    res$species_label <- get_common_name(res$best_species)

    # Build hover text
    res$hover <- vapply(seq_len(nrow(res)), function(i) {
      r <- res[i, ]
      second_sp <- if (is.na(r$second_species)) "No other species" else get_common_name(r$second_species)
      second_id <- if (is.na(r$second_identity)) "-" else paste0(round(r$second_identity, 1), "%")
      paste0(
        "<b>", r$peptide, "</b>",
        "<br>Best: ", get_common_name(r$best_species), " (", round(r$best_identity, 1), "%)",
        "<br>2nd: ", second_sp, " (", second_id, ")",
        "<br>Delta: ", round(r$delta_identity, 1), "%",
        "<br>Protein: ", r$protein_name
      )
    }, character(1))

    # Diagnostic threshold line
    diagnostic_threshold <- 15

    p <- plotly::plot_ly(res,
      y = ~pep_label, x = ~delta_identity, color = ~species_label,
      type = "bar", orientation = "h",
      hoverinfo = "text", text = ~hover
    ) %>%
      plotly::layout(
        title = list(text = "Species Resolution: Identity Gap Between Best and Second-Best Species",
          font = list(size = 14)),
        xaxis = list(title = "Delta Identity (%)", zeroline = FALSE),
        yaxis = list(title = "", tickfont = list(size = 9)),
        showlegend = TRUE,
        legend = list(title = list(text = "Best Species"), x = 0.7, y = 0.15),
        margin = list(l = 160),
        height = max(400, top_n * 22),
        shapes = list(
          list(type = "line", x0 = diagnostic_threshold, x1 = diagnostic_threshold,
            y0 = -0.5, y1 = top_n - 0.5,
            line = list(color = "#2e7d32", width = 2, dash = "dash"))
        ),
        annotations = list(
          list(x = diagnostic_threshold, y = top_n - 0.5,
            text = "Species-diagnostic", showarrow = FALSE,
            font = list(color = "#2e7d32", size = 11),
            xanchor = "left", xshift = 5)
        )
      )
    p %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Feature 3: Taxonomic Coverage Dot Plot ---
  output$dda_blast_taxonomic_coverage <- plotly::renderPlotly({
    blast <- tryCatch(blast_filtered(), error = function(e) NULL)
    req(!is.null(blast), nrow(blast) > 1)

    # Best hit per peptide-species combination
    blast_best <- blast[order(blast$pident, decreasing = TRUE), ]
    blast_best <- blast_best[!duplicated(paste(blast_best$peptide, blast_best$species)), ]

    # Parse protein name for grouping
    blast_best$prot_group <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", blast_best$subject))

    # Top 5 species by frequency of best hits
    best_per_pep <- blast_best[!duplicated(blast_best$peptide), ]
    sp_counts <- sort(table(best_per_pep$species), decreasing = TRUE)
    top_species <- names(sp_counts)[seq_len(min(5, length(sp_counts)))]

    # Subset to top species and get peptides that appear in at least 2 species
    sub <- blast_best[blast_best$species %in% top_species, ]
    pep_sp_count <- tapply(sub$species, sub$peptide, function(x) length(unique(x)))
    multi_sp_peps <- names(pep_sp_count[pep_sp_count >= 2])

    if (length(multi_sp_peps) == 0) {
      # Fallback: use all peptides with top species
      multi_sp_peps <- unique(sub$peptide)
    }

    # Limit to top 40 peptides by average identity for readability
    sub <- sub[sub$peptide %in% multi_sp_peps, ]
    pep_avg <- tapply(sub$pident, sub$peptide, mean, na.rm = TRUE)
    pep_avg <- sort(pep_avg, decreasing = TRUE)
    top_peps <- names(pep_avg)[seq_len(min(40, length(pep_avg)))]
    sub <- sub[sub$peptide %in% top_peps, ]

    # Assign protein group per peptide (from best hit overall)
    pep_prot <- blast_best[!duplicated(blast_best$peptide), c("peptide", "prot_group")]
    sub <- merge(sub, pep_prot[, c("peptide", "prot_group")],
      by = "peptide", all.x = TRUE, suffixes = c("", ".best"))
    sub$prot_group <- sub$prot_group.best
    sub$prot_group.best <- NULL

    # Sort peptides by protein group then identity for visual grouping
    pep_order <- unique(sub[order(sub$prot_group, -sub$pident), "peptide"])
    sub$peptide_label <- ifelse(nchar(sub$peptide) > 12,
      paste0(substr(sub$peptide, 1, 9), "..."), sub$peptide)
    # Maintain order
    label_map <- setNames(sub$peptide_label, sub$peptide)
    label_map <- label_map[!duplicated(names(label_map))]
    ordered_labels <- label_map[pep_order]
    # Make unique for factoring
    ordered_labels <- make.unique(ordered_labels, sep = " ")
    sub$peptide_label <- make.unique(sub$peptide_label, sep = " ")
    sub$peptide_label <- factor(sub$peptide_label, levels = ordered_labels)

    # Species label
    sub$species_label <- get_common_name(sub$species)

    # Hover text
    sub$hover <- paste0(
      "<b>", sub$peptide, "</b>",
      "<br>Species: ", sub$species_label,
      "<br>Identity: ", round(sub$pident, 1), "%",
      "<br>Protein: ", sub$prot_group
    )

    # Color by protein group (up to 10 colors)
    unique_prots <- unique(sub$prot_group)
    prot_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
    prot_color_map <- setNames(
      prot_colors[seq_len(min(length(unique_prots), length(prot_colors)))],
      unique_prots[seq_len(min(length(unique_prots), length(prot_colors)))]
    )

    n_species <- length(top_species)
    n_peps <- length(top_peps)
    plot_height <- max(400, n_peps * 16 + 100)

    p <- plotly::plot_ly()

    for (sp in top_species) {
      sp_data <- sub[sub$species == sp, ]
      if (nrow(sp_data) == 0) next

      sp_label <- get_common_name(sp)

      # Add lines connecting dots within same protein
      for (prot in unique(sp_data$prot_group)) {
        prot_data <- sp_data[sp_data$prot_group == prot, ]
        if (nrow(prot_data) > 1) {
          prot_data <- prot_data[order(prot_data$peptide_label), ]
          p <- p %>% plotly::add_trace(
            data = prot_data,
            x = ~pident, y = ~peptide_label,
            type = "scatter", mode = "lines",
            line = list(color = prot_color_map[prot] %||% "#999999", width = 1, dash = "dot"),
            showlegend = FALSE, hoverinfo = "skip",
            legendgroup = sp_label
          )
        }
      }

      # Add dots
      sp_data$dot_color <- vapply(sp_data$prot_group, function(pg) {
        prot_color_map[pg] %||% "#999999"
      }, character(1))

      p <- p %>% plotly::add_trace(
        data = sp_data,
        x = ~pident, y = ~peptide_label,
        type = "scatter", mode = "markers",
        marker = list(size = 8, color = sp_data$dot_color, opacity = 0.85,
          line = list(width = 1, color = "#333")),
        name = sp_label,
        hoverinfo = "text", text = ~hover,
        legendgroup = sp_label
      )
    }

    # Add facet-like annotations for species
    p %>% plotly::layout(
      title = list(
        text = paste0("Taxonomic Coverage: % Identity Across Top ", n_species, " Species"),
        font = list(size = 14)
      ),
      xaxis = list(title = "% Identity", range = c(
        max(40, min(sub$pident, na.rm = TRUE) - 5), 105)),
      yaxis = list(title = "", tickfont = list(size = 8)),
      showlegend = TRUE,
      legend = list(title = list(text = "Species"), x = 1.02, y = 1),
      margin = list(l = 120, r = 120),
      height = plot_height
    ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Peptide-Species heatmap (legacy, collapsible) ---
  output$dda_blast_heatmap <- plotly::renderPlotly({
    blast <- tryCatch(blast_filtered(), error = function(e) NULL)
    req(!is.null(blast), nrow(blast) > 1)

    # Best hit per peptide-species combination
    blast_best <- blast[order(blast$pident, decreasing = TRUE), ]
    blast_best <- blast_best[!duplicated(paste(blast_best$peptide, blast_best$species)), ]

    # Top 10 species by frequency
    sp_counts <- sort(table(blast_best$species), decreasing = TRUE)
    top_species <- names(sp_counts)[seq_len(min(10, length(sp_counts)))]

    # Top 50 peptides by average identity
    pep_avg <- tapply(blast_best$pident, blast_best$peptide, mean, na.rm = TRUE)
    pep_avg <- sort(pep_avg, decreasing = TRUE)
    top_peptides <- names(pep_avg)[seq_len(min(50, length(pep_avg)))]

    # Build matrix
    mat <- matrix(NA_real_, nrow = length(top_peptides), ncol = length(top_species),
      dimnames = list(top_peptides, top_species))

    sub <- blast_best[blast_best$peptide %in% top_peptides &
                      blast_best$species %in% top_species, ]
    for (i in seq_len(nrow(sub))) {
      mat[sub$peptide[i], sub$species[i]] <- sub$pident[i]
    }

    # Truncate long peptide labels
    row_labels <- ifelse(nchar(top_peptides) > 15,
      paste0(substr(top_peptides, 1, 12), "..."), top_peptides)

    plotly::plot_ly(
      z = mat,
      x = colnames(mat),
      y = row_labels,
      type = "heatmap",
      colorscale = list(
        list(0, "#ffffff"),
        list(0.5, "#ffcc80"),
        list(0.9, "#ef6c00"),
        list(1, "#b71c1c")
      ),
      zmin = 50, zmax = 100,
      hovertemplate = "Peptide: %{y}<br>Species: %{x}<br>Identity: %{z:.1f}%<extra></extra>",
      colorbar = list(title = "% Identity")
    ) %>%
      plotly::layout(
        title = list(text = "Peptide-Species Identity Matrix (top 50 x top 10)", font = list(size = 14)),
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "", tickfont = list(size = 9)),
        margin = list(l = 120, b = 80)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Enhanced BLAST results table ---
  output$dda_denovo_blast_table <- DT::renderDT({
    blast <- blast_with_species()

    # Parse protein name and accession from SwissProt subject: sp|ACC|PROT_SPECIES
    blast$accession <- sub("^[a-z]+\\|([^|]+)\\|.*", "\\1", blast$subject)
    blast$protein_name <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", blast$subject))
    # Create clickable UniProt link
    blast$protein_link <- paste0(
      '<a href="https://www.uniprot.org/uniprot/', blast$accession,
      '" target="_blank">', blast$protein_name, '</a>')

    display_df <- data.frame(
      Peptide     = blast$peptide,
      Accession   = blast$accession,
      Protein     = blast$protein_link,
      Species     = blast$species,
      Category    = blast$category,
      Identity    = round(blast$pident, 1),
      Length      = blast$length,
      E_Value     = formatC(blast$evalue, format = "e", digits = 2),
      Bitscore    = round(blast$bitscore, 1),
      Contaminant = ifelse(blast$contaminant_type == "Sample", "",
                     blast$contaminant_type),
      stringsAsFactors = FALSE
    )

    # Apply filter if set
    filt <- input$dda_blast_filter %||% "All"
    if (filt == "Conserved") display_df <- display_df[display_df$Category == "Conserved", ]
    if (filt == "Near-match") display_df <- display_df[display_df$Category == "Near-match", ]
    if (filt == "Distant") display_df <- display_df[display_df$Category == "Distant", ]

    DT::datatable(
      display_df,
      rownames = FALSE,
      filter   = "top",
      selection = "none",
      escape   = FALSE,
      options  = list(
        pageLength = 25,
        scrollX    = TRUE,
        order      = list(list(5, "desc")),
        dom        = "Bfrtip",
        buttons    = list("csv", "excel")
      ),
      extensions = "Buttons",
      caption = htmltools::tags$caption(
        style = "caption-side: top; font-weight: bold; color: #1565c0;",
        "DIAMOND BLAST hits for novel peptides against UniProt SwissProt + TrEMBL (SwissProt ~572k reviewed first, then TrEMBL on the misses)"
      )
    ) %>%
      DT::formatStyle("Category",
        backgroundColor = DT::styleEqual(
          c("Conserved", "Near-match", "Distant"),
          c("#e8f5e9", "#fff3e0", "#fce4ec")
        )
      ) %>%
      DT::formatStyle("Contaminant",
        backgroundColor = DT::styleEqual(
          c("Definite", "Possible"),
          c("#ffcdd2", "#fff9c4")
        ),
        fontWeight = DT::styleEqual(
          c("Definite", "Possible"),
          c("bold", "normal")
        )
      )
  })

  # --- Score distribution: MOVED to server_denovo_controls.R (with threshold line) ---

  # ==========================================================================
  #    PRIORITY 2: Per-Residue Confidence Visualization
  # ==========================================================================

  # Click handler for confirmed peptide table (uses filtered classification from server_denovo_controls.R)
  observeEvent(input$dda_denovo_confirmed_table_rows_selected, {
    sel <- input$dda_denovo_confirmed_table_rows_selected
    req(length(sel) > 0)
    sel_row <- sel[length(sel)]  # Use last selected row

    cls <- values$dda_filtered_classification %||% values$dda_casanovo_classification
    req(cls)
    confirmed <- cls$confirmed
    req(nrow(confirmed) >= sel_row)

    row <- confirmed[sel_row, ]
    html <- build_residue_confidence_html(row, values$denovo_blast %||% values$dda_casanovo_blast)
    shinyjs::html("dda_confirmed_residue_viz", html)
  })

  # Click handler for novel peptide table (uses filtered classification from server_denovo_controls.R)
  observeEvent(input$dda_denovo_novel_table_rows_selected, {
    sel <- input$dda_denovo_novel_table_rows_selected
    req(length(sel) > 0)
    sel_row <- sel[length(sel)]

    cls <- values$denovo_classification %||% values$dda_filtered_classification %||% values$dda_casanovo_classification
    req(cls)
    novel <- cls$novel
    req(nrow(novel) >= sel_row)

    row <- novel[sel_row, ]
    html <- build_residue_confidence_html(row, values$denovo_blast %||% values$dda_casanovo_blast)
    shinyjs::html("dda_novel_residue_viz", html)
  })

  # Helper: build per-residue colored HTML for a PSM row
  build_residue_confidence_html <- function(row, blast_data = NULL) {
    seq <- row$seq_stripped
    aa_str <- row$aa_scores
    score <- round(row$score, 3)
    charge <- row$charge
    mean_aa <- if (!is.null(row$mean_aa_score) && !is.na(row$mean_aa_score)) {
      round(row$mean_aa_score, 3)
    } else {
      NA
    }

    # Parse per-residue scores
    residues <- strsplit(seq, "")[[1]]
    aa_vals <- if (!is.na(aa_str) && nzchar(aa_str) && aa_str != "null") {
      as.numeric(strsplit(aa_str, ",")[[1]])
    } else {
      NULL
    }

    # Build colored sequence HTML
    if (!is.null(aa_vals) && length(aa_vals) == length(residues)) {
      colored_spans <- vapply(seq_along(residues), function(i) {
        v <- aa_vals[i]
        color <- if (is.na(v)) {
          "#999999"
        } else if (v >= 0.95) {
          "#2e7d32"  # green
        } else if (v >= 0.7) {
          "#f9a825"  # yellow/amber
        } else {
          "#c62828"  # red
        }
        bg <- if (!is.na(v) && v < 0.7) {
          "background: #fce4ec;"
        } else {
          ""
        }
        sprintf(
          '<span style="color: %s; font-weight: bold; font-family: monospace; font-size: 1.3em; %s" title="%.3f">%s</span>',
          color, bg, if (is.na(v)) 0 else v, residues[i]
        )
      }, character(1))
      seq_html <- paste(colored_spans, collapse = "")
    } else {
      seq_html <- paste0(
        '<span style="font-family: monospace; font-size: 1.3em; color: #555;">',
        seq, '</span>')
    }

    # BLAST hit info if available
    blast_html <- ""
    if (!is.null(blast_data) && nrow(blast_data) > 0) {
      clean_seq <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(seq))
      hit <- blast_data[blast_data$peptide == clean_seq, ]
      if (nrow(hit) > 0) {
        hit <- hit[which.max(hit$bitscore), ]
        prot_name <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", hit$subject))
        acc <- sub("^[a-z]+\\|([^|]+)\\|.*", "\\1", hit$subject)
        blast_html <- sprintf(
          '<div style="margin-top: 8px; padding: 8px; background: #e3f2fd; border-radius: 6px; font-size: 0.9em;">
            <strong>BLAST Hit:</strong>
            <a href="https://www.uniprot.org/uniprot/%s" target="_blank">%s</a>
            (%s) | Identity: %.1f%% | E-value: %s
          </div>',
          acc, prot_name, hit$species,
          hit$pident %||% hit$identity,
          formatC(hit$evalue, format = "e", digits = 2)
        )
      }
    }

    # Legend
    legend_html <- paste0(
      '<div style="margin-top: 6px; font-size: 0.8em; color: #666;">',
      '<span style="color: #2e7d32; font-weight: bold;">Green</span> >= 0.95 | ',
      '<span style="color: #f9a825; font-weight: bold;">Yellow</span> 0.70-0.95 | ',
      '<span style="color: #c62828; font-weight: bold;">Red</span> < 0.70 (potential error)',
      '</div>'
    )

    # Stats line
    stats_parts <- c(paste0("Score: ", score), paste0("Charge: ", charge, "+"))
    if (!is.na(mean_aa)) stats_parts <- c(stats_parts, paste0("Mean AA: ", mean_aa))
    if (!is.null(aa_vals)) {
      n_high <- sum(aa_vals >= 0.95, na.rm = TRUE)
      n_low <- sum(aa_vals < 0.7, na.rm = TRUE)
      stats_parts <- c(stats_parts,
        paste0(n_high, "/", length(aa_vals), " high-conf"),
        if (n_low > 0) paste0(n_low, " low-conf") else NULL
      )
    }
    stats_html <- paste0(
      '<div style="margin-top: 4px; font-size: 0.85em; color: #444;">',
      paste(stats_parts, collapse = " | "),
      '</div>'
    )

    paste0(
      '<div style="padding: 12px; background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 8px; margin-top: 10px;">',
      '<div style="margin-bottom: 6px; font-weight: 600; color: #333;">Per-Residue Confidence</div>',
      seq_html, stats_html, legend_html, blast_html,
      '</div>'
    )
  }

  # ==========================================================================
  #    PRIORITY 3: Length and Charge Distribution QC
  # ==========================================================================

  output$dda_denovo_length_charge_qc <- plotly::renderPlotly({
    req(values$dda_casanovo_classification)
    cls <- values$dda_casanovo_classification

    confirmed <- cls$confirmed
    novel <- cls$novel
    req(nrow(confirmed) > 0 || nrow(novel) > 0)

    # Compute peptide lengths
    conf_lengths <- nchar(confirmed$seq_stripped)
    novel_lengths <- nchar(novel$seq_stripped)

    plot_df <- data.frame(
      length = c(conf_lengths, novel_lengths),
      type   = c(rep("Confirmed", length(conf_lengths)), rep("Novel", length(novel_lengths))),
      stringsAsFactors = FALSE
    )

    colors <- c("Confirmed" = "#2ecc71", "Novel" = "#e67e22")

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = length, fill = type)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.75, position = "dodge") +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::geom_vline(xintercept = c(7, 25), linetype = "dashed", color = "#c62828", alpha = 0.5) +
      ggplot2::annotate("rect", xmin = 7, xmax = 25, ymin = -Inf, ymax = Inf,
        fill = "#e8f5e9", alpha = 0.15) +
      ggplot2::annotate("text", x = 16, y = Inf, label = "Expected tryptic range (7-25)",
        vjust = 1.5, color = "#2e7d32", size = 3.5) +
      ggplot2::labs(
        x = "Peptide Length (amino acids)",
        y = "Count",
        fill = "Type",
        subtitle = sprintf("Confirmed: median %d aa | Novel: median %d aa",
          as.integer(median(conf_lengths)), as.integer(median(novel_lengths)))
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")

    plotly::ggplotly(p) %>%
      plotly::layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.05)) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  output$dda_denovo_charge_dist <- plotly::renderPlotly({
    req(values$dda_casanovo_classification)
    cls <- values$dda_casanovo_classification
    confirmed <- cls$confirmed
    novel <- cls$novel

    plot_df <- data.frame(
      charge = c(confirmed$charge, novel$charge),
      type   = c(rep("Confirmed", nrow(confirmed)), rep("Novel", nrow(novel))),
      stringsAsFactors = FALSE
    )
    plot_df$charge <- factor(plot_df$charge)

    colors <- c("Confirmed" = "#2ecc71", "Novel" = "#e67e22")

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = charge, fill = type)) +
      ggplot2::geom_bar(position = "dodge", alpha = 0.8) +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(
        x = "Charge State",
        y = "Count",
        fill = "Type",
        subtitle = "Expected: 2+ and 3+ dominant for tryptic peptides"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")

    plotly::ggplotly(p) %>%
      plotly::layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.05)) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  output$dda_denovo_qc_summary <- renderUI({
    req(values$dda_casanovo_classification)
    cls <- values$dda_casanovo_classification
    novel <- cls$novel
    confirmed <- cls$confirmed

    novel_lengths <- nchar(novel$seq_stripped)
    novel_charges <- novel$charge
    conf_lengths <- nchar(confirmed$seq_stripped)
    conf_charges <- confirmed$charge

    # Tryptic quality: length 7-25 and charge 2-3
    novel_tryptic <- sum(novel_lengths >= 7 & novel_lengths <= 25 &
                         novel_charges >= 2 & novel_charges <= 3)
    novel_pct <- round(100 * novel_tryptic / max(nrow(novel), 1), 1)

    conf_tryptic <- sum(conf_lengths >= 7 & conf_lengths <= 25 &
                        conf_charges >= 2 & conf_charges <= 3)
    conf_pct <- round(100 * conf_tryptic / max(nrow(confirmed), 1), 1)

    # Flags
    n_short <- sum(novel_lengths < 7)
    n_long <- sum(novel_lengths > 25)
    n_charge1 <- sum(novel_charges == 1)

    flags <- character(0)
    if (n_short > 0) flags <- c(flags, sprintf("%d novel peptides < 7 aa (unreliable)", n_short))
    if (n_long > 0) flags <- c(flags, sprintf("%d novel peptides > 25 aa (unusual)", n_long))
    if (n_charge1 > 0) flags <- c(flags, sprintf("%d novel peptides at 1+ charge (suspicious)", n_charge1))

    tagList(
      div(style = "padding: 12px; background: #f0f7ff; border-radius: 8px; margin-top: 10px;",
        tags$p(style = "margin: 0; font-size: 0.95em;",
          sprintf("%.1f%% of novel peptides have tryptic characteristics (length 7-25, charge 2-3+).", novel_pct),
          sprintf(" Confirmed peptides: %.1f%%.", conf_pct)
        ),
        if (length(flags) > 0) {
          tags$div(style = "margin-top: 8px;",
            lapply(flags, function(f) {
              tags$p(style = "margin: 2px 0; color: #c62828; font-size: 0.9em;",
                icon("triangle-exclamation"), " ", f)
            })
          )
        }
      )
    )
  })

  # ==========================================================================
  #    PRIORITY 4: Modification Tracking (Deamidation as Authenticity Marker)
  # ==========================================================================

  output$dda_denovo_modifications <- renderUI({
    req(values$dda_casanovo_psms)
    psms <- values$dda_casanovo_psms
    req(nrow(psms) > 0)

    seqs <- psms$sequence

    # Parse modification masses from sequences
    # Patterns: N+0.984 (deamidation N), Q+0.984 (deamidation Q), M+15.995 (oxidation M)
    # Also: [+mass] format

    # Count residues and modifications
    n_total_psms <- nrow(psms)

    # Count N residues and N-deamidation (+0.984 after N)
    n_N_residues <- sum(vapply(seqs, function(s) {
      stripped <- gsub("\\+[0-9.]+", "", gsub("\\[|\\]", "", s))
      nchar(gsub("[^N]", "", stripped))
    }, integer(1)))

    n_N_deamid <- sum(grepl("N[+]0\\.98[0-9]|N\\[\\+0\\.98[0-9]", seqs))

    # Count Q residues and Q-deamidation
    n_Q_residues <- sum(vapply(seqs, function(s) {
      stripped <- gsub("\\+[0-9.]+", "", gsub("\\[|\\]", "", s))
      nchar(gsub("[^Q]", "", stripped))
    }, integer(1)))

    n_Q_deamid <- sum(grepl("Q[+]0\\.98[0-9]|Q\\[\\+0\\.98[0-9]", seqs))

    # Count M residues and oxidation (+15.995)
    n_M_residues <- sum(vapply(seqs, function(s) {
      stripped <- gsub("\\+[0-9.]+", "", gsub("\\[|\\]", "", s))
      nchar(gsub("[^M]", "", stripped))
    }, integer(1)))

    n_M_oxidized <- sum(grepl("M[+]15\\.99[0-9]|M\\[\\+15\\.99[0-9]", seqs))

    # Any PSM with any modification
    n_modified_psms <- sum(grepl("[+][0-9]", seqs))

    # Rates
    n_deamid_rate <- if (n_N_residues > 0) round(100 * n_N_deamid / n_N_residues, 2) else 0
    q_deamid_rate <- if (n_Q_residues > 0) round(100 * n_Q_deamid / n_Q_residues, 2) else 0
    m_oxid_rate <- if (n_M_residues > 0) round(100 * n_M_oxidized / n_M_residues, 2) else 0
    nq_ratio <- if (n_Q_deamid > 0) round(n_N_deamid / n_Q_deamid, 1) else
      if (n_N_deamid > 0) "Inf" else "N/A"

    # Authenticity assessment
    authenticity_html <- ""
    if (n_N_deamid > 0 || n_Q_deamid > 0) {
      if (is.numeric(nq_ratio) && nq_ratio > 3) {
        authenticity_html <- paste0(
          '<div style="margin-top: 10px; padding: 10px; background: #e8f5e9; border-radius: 6px; border-left: 4px solid #2e7d32;">',
          '<strong style="color: #2e7d32;">Authenticity signal detected:</strong> ',
          'N/Q deamidation ratio = ', nq_ratio, '. ',
          'High N-deamidation with low Q-deamidation is characteristic of genuine ancient proteins ',
          '(spontaneous asparagine deamidation accumulates over time, while glutamine deamidation is more random).',
          '</div>'
        )
      } else if (is.numeric(nq_ratio) && nq_ratio < 1.5) {
        authenticity_html <- paste0(
          '<div style="margin-top: 10px; padding: 10px; background: #fff3e0; border-radius: 6px; border-left: 4px solid #e65100;">',
          '<strong style="color: #e65100;">Low authenticity signal:</strong> ',
          'N/Q deamidation ratio = ', nq_ratio, '. ',
          'Similar N and Q deamidation rates may indicate sample preparation artifacts ',
          'rather than time-dependent degradation.',
          '</div>'
        )
      }
    }

    tagList(
      div(class = "row", style = "margin-bottom: 12px;",
        div(class = "col-md-3",
          div(style = "background: #f3e5f5; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #7b1fa2;",
              paste0(n_modified_psms, "/", n_total_psms)),
            tags$small("Modified PSMs")
          )
        ),
        div(class = "col-md-3",
          div(style = "background: #e8f5e9; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #2e7d32;",
              paste0(n_deamid_rate, "%")),
            tags$small(paste0("N-Deamidation (", n_N_deamid, "/", n_N_residues, " N)"))
          )
        ),
        div(class = "col-md-3",
          div(style = "background: #fff3e0; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #e65100;",
              paste0(q_deamid_rate, "%")),
            tags$small(paste0("Q-Deamidation (", n_Q_deamid, "/", n_Q_residues, " Q)"))
          )
        ),
        div(class = "col-md-3",
          div(style = "background: #e3f2fd; padding: 12px; border-radius: 8px; text-align: center;",
            tags$h4(style = "margin: 0; color: #1565c0;",
              paste0(m_oxid_rate, "%")),
            tags$small(paste0("M-Oxidation (", n_M_oxidized, "/", n_M_residues, " M)"))
          )
        )
      ),
      HTML(authenticity_html)
    )
  })

  # Modification types bar chart
  output$dda_denovo_mod_bar <- plotly::renderPlotly({
    req(values$dda_casanovo_psms)
    psms <- values$dda_casanovo_psms
    req(nrow(psms) > 0)

    seqs <- psms$sequence

    # Extract all modification masses
    mod_masses <- unlist(regmatches(seqs, gregexpr("[A-Z][+][0-9.]+", seqs)))

    if (length(mod_masses) == 0) {
      return(plotly::plot_ly() %>% plotly::layout(
        title = list(text = "No modifications detected", font = list(size = 14))))
    }

    # Classify modifications
    mod_type <- vapply(mod_masses, function(m) {
      mass <- as.numeric(sub("^[A-Z][+]", "", m))
      aa <- substr(m, 1, 1)
      if (abs(mass - 0.984) < 0.01 && aa == "N") return("N-Deamidation")
      if (abs(mass - 0.984) < 0.01 && aa == "Q") return("Q-Deamidation")
      if (abs(mass - 15.995) < 0.01) return("Oxidation (M)")
      if (abs(mass - 57.021) < 0.01) return("Carbamidomethyl (C)")
      if (abs(mass - 42.011) < 0.01) return("Acetylation")
      paste0(aa, "+", round(mass, 3))
    }, character(1))

    mod_counts <- sort(table(mod_type), decreasing = TRUE)
    top_n <- min(15, length(mod_counts))
    mod_df <- data.frame(
      Modification = factor(names(mod_counts)[seq_len(top_n)],
        levels = rev(names(mod_counts)[seq_len(top_n)])),
      Count = as.numeric(mod_counts[seq_len(top_n)]),
      stringsAsFactors = FALSE
    )

    # Color deamidation types specially
    mod_df$color <- ifelse(grepl("Deamid", mod_df$Modification), "#2e7d32",
      ifelse(grepl("Oxid", mod_df$Modification), "#1565c0", "#7b1fa2"))

    plotly::plot_ly(mod_df,
      y = ~Modification, x = ~Count,
      type = "bar", orientation = "h",
      marker = list(color = mod_df$color)
    ) %>%
      plotly::layout(
        title = list(text = "Modification Types", font = list(size = 14)),
        xaxis = list(title = "PSM Count"),
        yaxis = list(title = ""),
        margin = list(l = 160)
      ) %>%
      plotly::config(toImageButtonOptions = list(format = "svg", scale = 2))
  })

  # --- Manuscript summary: MOVED to server_denovo_controls.R (confidence-filtered) ---

  # ============================================================================
  #    DIAMOND BLAST for DDA Novel Peptides
  # ============================================================================

  observeEvent(input$dda_run_diamond_blast, {
    cls <- values$denovo_classification %||% values$dda_casanovo_classification
    req(cls)
    novel <- cls$novel
    req(nrow(novel) > 0)

    # Strip modification masses — keep only amino acid letters
    clean_seqs <- gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", toupper(novel$seq_stripped))
    novel_peptides <- unique(clean_seqs[nzchar(clean_seqs)])

    if (length(novel_peptides) == 0) {
      showNotification("No novel peptides to BLAST.", type = "warning")
      return()
    }

    ssh_cfg <- tryCatch(dda_ssh_config(), error = function(e) NULL)
    if (is.null(ssh_cfg) || !isTRUE(values$ssh_connected)) {
      showNotification(
        "SSH connection required. DIAMOND BLAST runs on HPC via SSH.",
        type = "error"
      )
      return()
    }

    tryCatch({
      withProgress(message = "Running DIAMOND BLAST on HPC...", value = 0.1, {

        output_dir <- values$dda_output_dir %||% paste0("/tmp/delimp_dda_denovo_", Sys.getpid())
        denovo_dir <- file.path(output_dir, "denovo")

        # Create remote directory
        ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(denovo_dir)), timeout = 15)
        setProgress(0.2, detail = "Created remote directory")

        # Write query FASTA locally, then SCP upload
        query_fasta_local <- tempfile(fileext = ".fasta")
        query_lines <- paste0(">", novel_peptides, "\n", novel_peptides)
        writeLines(query_lines, query_fasta_local)

        query_fasta_remote <- file.path(denovo_dir, "novel_casanovo_queries.fasta")
        scp_upload(ssh_cfg, query_fasta_local, query_fasta_remote)
        setProgress(0.3, detail = "Uploaded query FASTA")

        # Run DIAMOND blastp against pre-built SwissProt DB
        diamond_bin <- "diamond"
        setProgress(0.5, detail = "Running BLAST against SwissProt...")

        blast_out_remote <- file.path(denovo_dir, "novel_casanovo_blast.tsv")
        blast_cmd <- paste0(
          "module load diamond 2>/dev/null && ",
          diamond_bin, " blastp",
          " --query ", shQuote(query_fasta_remote),
          " --db ", shQuote(swissprot_dmnd),
          " --out ", shQuote(blast_out_remote),
          " --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore",
          " --sensitive --id 50 --max-target-seqs 5",
          " --threads 4 --quiet"
        )
        blast_res <- ssh_exec(ssh_cfg, blast_cmd, login_shell = TRUE, timeout = 600)
        if ((blast_res$status %||% 0L) != 0L) {
          showNotification(
            paste("DIAMOND blastp failed:", paste(blast_res$stderr, collapse = "\n")),
            type = "error"
          )
          return()
        }
        setProgress(0.8, detail = "BLAST complete, downloading results")

        # Download results
        blast_out_local <- tempfile(fileext = ".tsv")
        scp_download(ssh_cfg, blast_out_remote, blast_out_local)

        if (file.exists(blast_out_local) && file.size(blast_out_local) > 0) {
          hits <- data.table::fread(blast_out_local, header = FALSE)
          names(hits) <- c("query", "subject", "identity", "length",
                           "qlen", "slen", "evalue", "bitscore")

          # Query column IS the peptide sequence (FASTA header = sequence)
          hits$peptide <- hits$query

          # Extract protein accession (handles sp|ACC|NAME format)
          hits$protein <- stringr::str_extract(hits$subject, "(?<=\\|)[^|]+(?=\\|)")
          no_match <- is.na(hits$protein)
          if (any(no_match)) {
            hits$protein[no_match] <- hits$subject[no_match]
          }

          # Species via the single taxonomy-aware definition (LCA join when an
          # nr LCA table is loaded; UniProt mnemonic otherwise).
          hits$species <- dda_blast_species(hits$peptide, hits$subject, values$dda_lca)
          # Classify by identity
          hits$category <- ifelse(
            hits$identity >= 100, "Conserved",
            ifelse(hits$identity >= 90, "Near-match", "Distant")
          )

          values$dda_casanovo_blast <- as.data.frame(hits)
          n_hits <- length(unique(hits$peptide))
          n_proteins <- length(unique(hits$protein))
          showNotification(
            sprintf("DIAMOND BLAST: %d novel peptides mapped to %d protein hits.",
                    n_hits, n_proteins),
            type = "message"
          )
          add_to_log(
            sprintf("DDA DIAMOND BLAST: %d peptides -> %d hits", n_hits, nrow(hits)),
            "denovo"
          )
        } else {
          values$dda_casanovo_blast <- data.frame()
          showNotification("No DIAMOND BLAST hits found.", type = "warning")
          add_to_log("DDA DIAMOND BLAST: no hits found", "denovo")
        }

        setProgress(1.0, detail = "Done")
      })

    }, error = function(e) {
      showNotification(
        paste("DIAMOND BLAST error:", conditionMessage(e)),
        type = "error"
      )
      add_to_log(paste("DDA DIAMOND BLAST error:", conditionMessage(e)), "error")
    })
  })

  # ============================================================================
  #  ADAPTER: Populate unified de novo reactives from Casanovo data
  # ============================================================================

  observe({
    req(values$dda_casanovo_classification)
    values$denovo_classification <- values$dda_casanovo_classification
    values$denovo_psms <- values$dda_casanovo_psms
    values$denovo_engine <- "casanovo"
    values$denovo_reference <- "Sage"
  })

  observe({
    req(values$dda_casanovo_blast)
    values$denovo_blast <- values$dda_casanovo_blast
  })

  # ============================================================================
  #  INFO MODALS — De Novo Sub-tabs
  # ============================================================================

  observeEvent(input$denovo_confirmed_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " Sage DB hits"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("Peptides sequenced de novo by Casanovo that were also identified by database search (Sage)."),
        tags$ul(
          tags$li(strong("Sequence: "), "The amino acid sequence predicted by Casanovo."),
          tags$li(strong("Score: "), "Casanovo's confidence score (0-1). Higher = more reliable."),
          tags$li(strong("AA Scores: "), "Per-residue confidence. Each position gets its own probability."),
          tags$li(strong("Charge: "), "Precursor charge state from the spectrum."),
          tags$li(strong("Matched Protein: "), "Sage's database match for the same spectrum.")
        ),
        p("Confirmed peptides validate the de novo algorithm — these sequences match known proteins. ",
          "The confirmation rate (confirmed / total) indicates de novo sequencing accuracy."),
        tags$hr(),
        tags$p(style = "color: #666; font-size: 0.85em;",
          icon("camera"), " Click the camera icon on any plot to download as SVG for publication figures.")
      )
    ))
  })

  observeEvent(input$denovo_novel_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " De novo only peptides"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("Peptides sequenced de novo by Casanovo that were ", strong("not"), " found by database search."),
        p("These are the most scientifically interesting sequences — they may represent:"),
        tags$ul(
          tags$li("Species-specific proteins not in the search database"),
          tags$li("Post-translationally modified peptides missed by database search"),
          tags$li("Degraded/ancient peptides (in paleoproteomics)"),
          tags$li("Sequencing errors (check per-residue AA scores)")
        ),
        p("Use the ", strong("DIAMOND BLAST"), " tab to search these against UniProt SwissProt ",
          "and identify the closest known proteins."),
        tags$hr(),
        tags$p(style = "color: #666; font-size: 0.85em;",
          "Click a row to see per-residue confidence coloring below the table. ",
          "Green = high confidence, Red = potential error.")
      )
    ))
  })

  observeEvent(input$denovo_blast_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " DIAMOND BLAST Analysis"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("DIAMOND BLAST searches novel de novo peptides against UniProt SwissProt + TrEMBL (SwissProt ~572k reviewed first, then TrEMBL on the misses) ",
          "to identify the closest known homologs."),
        tags$h6("Key Visualizations"),
        tags$ul(
          tags$li(strong("Species Donut: "), "Distribution of best-hit species per peptide. Reveals sample composition."),
          tags$li(strong("Identity Histogram: "), "Distribution of BLAST identity scores. ",
            "100% = exact match, 90-99% = near-match (potential variant), <90% = distant homolog."),
          tags$li(strong("Top Proteins: "), "Proteins ranked by number of matching de novo peptides."),
          tags$li(strong("Species Resolution: "), "For each peptide, computes the identity gap (delta) between the ",
            "best-matching species and the second-best species. ",
            "Example: 95% identity to chicken, 70% to pigeon = delta of 25%. ",
            "The vertical dashed line at delta = 15% separates: ",
            tags$ul(
              tags$li("Right of line (delta > 15%): ", strong("Species-diagnostic"), " — this peptide is specific to one species. Strong evidence for species ID."),
              tags$li("Left of line (delta < 15%): ", strong("Conserved"), " — similar identity to multiple species. Less useful for distinguishing species.")
            )),
          tags$li(strong("Taxonomic Coverage: "), "Dot plot showing peptide identity across species, grouped by protein.")
        ),
        tags$h6("Interpretation"),
        tags$ul(
          tags$li(strong("Conserved (100%): "), "Identical to known protein — high confidence."),
          tags$li(strong("Near-match (90-99%): "), "1-2 AA differences — potential species variant or PTM artifact."),
          tags$li(strong("Distant (<90%): "), "Low homology — novel protein family or sequencing error.")
        ),
        tags$h6("Contaminant Filtering"),
        p("When 'Exclude contaminant proteins' is checked, common lab contaminants (keratins, trypsin, BSA) ",
          "are removed. Feather keratins are NOT contaminants — they're real sample proteins."),
        tags$hr(),
        tags$p(style = "color: #666; font-size: 0.85em;",
          icon("camera"), " Click the camera icon on any plot to download as SVG.")
      )
    ))
  })

  observeEvent(input$denovo_mods_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("question-circle"), " Modification Analysis"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        p("Post-translational modifications detected in de novo peptide sequences."),
        tags$h6("Key Modifications"),
        tags$ul(
          tags$li(strong("Oxidation (M, +15.995): "), "Methionine oxidation — common artifact from sample handling."),
          tags$li(strong("Deamidation (N, +0.984): "), "Asparagine deamidation — in paleoproteomics, indicates authentic ancient protein. ",
            "Time-dependent process; higher rates = older protein."),
          tags$li(strong("Deamidation (Q, +0.984): "), "Glutamine deamidation — less time-dependent than N-deamidation. ",
            "High Q-deamidation suggests contamination or sample handling artifact."),
          tags$li(strong("Carbamidomethyl (C, +57.021): "), "Cysteine alkylation — expected from sample preparation.")
        ),
        tags$h6("Paleoproteomics Authenticity"),
        p("High N-deamidation + Low Q-deamidation = authentic ancient protein (genuine endogenous). ",
          "High Q-deamidation or very low N-deamidation = possible contamination or modern protein.")
      )
    ))
  })

  # ==========================================================================
  #  Force key outputs to evaluate even when De Novo tab is hidden.
  #  Required for session restore — values are set before tab becomes visible.
  # ==========================================================================
  outputOptions(output, "dda_blast_summary_cards",   suspendWhenHidden = FALSE)
  outputOptions(output, "dda_blast_diagnostic_card", suspendWhenHidden = FALSE)

}
