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

  observeEvent(input$load_dda_results, {
    showModal(modalDialog(
      title = "Load DDA Results from HPC",
      textInput("dda_load_path", "Output directory on HPC",
        placeholder = "/quobyte/proteomics-grp/de-limp/brettsp/dda_output/dda_search"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("dda_load_confirm", "Load", class = "btn-primary", icon = icon("download"))
      )
    ))
  })

  observeEvent(input$dda_load_confirm, {
    req(nzchar(input$dda_load_path))
    ssh_cfg <- dda_ssh_config()
    remote_dir <- trimws(input$dda_load_path)

    withProgress(message = "Loading DDA results...", value = 0.1, {
      tryCatch({
        # Download results.sage.tsv
        local_tmp <- file.path(tempdir(), "dda_load")
        dir.create(local_tmp, showWarnings = FALSE, recursive = TRUE)

        results_remote <- file.path(remote_dir, "results.sage.tsv")
        results_local <- file.path(local_tmp, "results.sage.tsv")
        scp_download(ssh_cfg, results_remote, results_local)

        if (!file.exists(results_local) || file.info(results_local)$size < 100) {
          showNotification("results.sage.tsv not found or empty.", type = "error")
          return()
        }
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
        values$dda_output_dir   <- remote_dir
        values$dda_status       <- "loaded"

        setProgress(0.5, detail = "Checking for Casanovo results...")

        # Check for Casanovo mztab files
        mztab_check <- ssh_exec(ssh_cfg,
          paste0("ls ", shQuote(file.path(remote_dir, "casanovo", "mztab")), "/*.mztab 2>/dev/null"),
          timeout = 10)
        if (mztab_check$status == 0 && length(mztab_check$stdout) > 0) {
          mztab_remote <- trimws(mztab_check$stdout)
          mztab_remote <- mztab_remote[nzchar(mztab_remote)]
          if (length(mztab_remote) > 0) {
            mztab_local_dir <- file.path(local_tmp, "mztab")
            dir.create(mztab_local_dir, showWarnings = FALSE)
            for (mt in mztab_remote) {
              tryCatch(scp_download(ssh_cfg, mt, file.path(mztab_local_dir, basename(mt))),
                error = function(e) NULL)
            }
            mztab_local <- list.files(mztab_local_dir, pattern = "\\.mztab$", full.names = TRUE)
            if (length(mztab_local) > 0) {
              casanovo_psms <- parse_casanovo_mztab(mztab_local, score_threshold = 0.9)
              if (nrow(casanovo_psms) > 0) {
                classified <- classify_dda_denovo(casanovo_psms, parsed$psms)
                values$dda_casanovo_psms <- casanovo_psms
                values$dda_casanovo_classification <- classified
                values$dda_casanovo_status <- "done"
              }
            }
          }
        }

        # Auto-submit DIAMOND BLAST SLURM job on novel peptides (SwissProt DB)
        if (!is.null(values$dda_casanovo_classification) &&
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
            fasta_lines <- unlist(lapply(seq_along(unique_seqs), function(i) {
              c(paste0(">denovo_", i, " ", unique_seqs[i]), unique_seqs[i])
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
'#!/bin/bash
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
  --threads 8

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
            }
          }, error = function(e) message("[DDA] Auto-BLAST submission failed: ", e$message))
        }

        # Also load existing MGF-based Casanovo results from raw directory
        setProgress(0.7, detail = "Checking for existing mztab files...")
        raw_dir <- trimws(input$dda_raw_dir %||% "")
        if (nzchar(raw_dir)) {
          mztab_raw_check <- ssh_exec(ssh_cfg,
            paste0("ls ", shQuote(raw_dir), "/*.mztab 2>/dev/null"),
            timeout = 10)
          if (mztab_raw_check$status == 0 && length(mztab_raw_check$stdout) > 0 &&
              is.null(values$dda_casanovo_psms)) {
            mztab_raw <- trimws(mztab_raw_check$stdout)
            mztab_raw <- mztab_raw[nzchar(mztab_raw)]
            if (length(mztab_raw) > 0) {
              mztab_dl_dir <- file.path(local_tmp, "mztab_raw")
              dir.create(mztab_dl_dir, showWarnings = FALSE)
              for (mt in mztab_raw) {
                tryCatch(scp_download(ssh_cfg, mt, file.path(mztab_dl_dir, basename(mt))),
                  error = function(e) NULL)
              }
              mztab_dl <- list.files(mztab_dl_dir, pattern = "\\.mztab$", full.names = TRUE)
              if (length(mztab_dl) > 0) {
                casanovo_psms <- parse_casanovo_mztab(mztab_dl, score_threshold = 0.9)
                if (nrow(casanovo_psms) > 0) {
                  classified <- classify_dda_denovo(casanovo_psms, parsed$psms)
                  values$dda_casanovo_psms <- casanovo_psms
                  values$dda_casanovo_classification <- classified
                  values$dda_casanovo_status <- "done"
                }
              }
            }
          }
        }

        # Load existing DIAMOND BLAST results if present
        setProgress(0.85, detail = "Checking for BLAST results...")
        tryCatch({
          blast_remote <- file.path(remote_dir, "denovo", "blast_results.tsv")
          blast_check <- ssh_exec(ssh_cfg, paste("test -s", shQuote(blast_remote), "&& echo YES || echo NO"), timeout = 10)
          if (grepl("YES", paste(blast_check$stdout, collapse = ""))) {
            blast_local <- file.path(local_tmp, "blast_results.tsv")
            scp_download(ssh_cfg, blast_remote, blast_local)
            if (file.exists(blast_local) && file.info(blast_local)$size > 0) {
              blast_df <- read.delim(blast_local, header = FALSE, stringsAsFactors = FALSE)
              colnames(blast_df) <- c("query", "subject", "pident", "length", "mismatch",
                "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
              # Parse species from SwissProt IDs: sp|ACC|PROT_SPECIES
              blast_df$species <- sub(".*_", "", blast_df$subject)
              # Map denovo IDs back to peptide sequences from FASTA headers
              fasta_remote <- file.path(remote_dir, "denovo", "novel_peptides.fasta")
              fasta_local <- file.path(local_tmp, "novel_peptides.fasta")
              tryCatch(scp_download(ssh_cfg, fasta_remote, fasta_local), error = function(e) NULL)
              if (file.exists(fasta_local)) {
                fasta_lines <- readLines(fasta_local, warn = FALSE)
                header_lines <- grep("^>", fasta_lines, value = TRUE)
                id_to_seq <- setNames(
                  sub("^>\\S+\\s+", "", header_lines),
                  sub("^>(\\S+).*", "\\1", header_lines)
                )
                blast_df$peptide <- unname(id_to_seq[blast_df$query])
                blast_df$peptide[is.na(blast_df$peptide)] <- blast_df$query[is.na(blast_df$peptide)]
              } else {
                blast_df$peptide <- blast_df$query
              }
              # Category
              blast_df$category <- ifelse(blast_df$pident == 100, "Conserved",
                ifelse(blast_df$pident >= 90, "Near-match", "Distant"))
              # Best hit per peptide
              blast_best <- blast_df[order(-blast_df$bitscore), ]
              blast_best <- blast_best[!duplicated(blast_best$peptide), ]
              values$dda_casanovo_blast <- blast_best
              message("[DDA] Loaded BLAST results: ", nrow(blast_best), " peptides with hits")
            }
          }
        }, error = function(e) message("[DDA] BLAST load: ", e$message))

        setProgress(0.95, detail = "Done!")

        n_psms <- nrow(parsed$psms)
        n_casanovo <- if (!is.null(values$dda_casanovo_psms)) nrow(values$dda_casanovo_psms) else 0
        removeModal()
        showNotification(
          sprintf("Loaded: %s Sage PSMs%s",
            format(n_psms, big.mark = ","),
            if (n_casanovo > 0) sprintf(", %s Casanovo PSMs", format(n_casanovo, big.mark = ",")) else ""),
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
    exp_name   <- trimws(input$dda_experiment_name %||% "dda_search")

    # Resolve FASTA path from whichever source was used
    fasta_source <- input$dda_fasta_source %||% "browse"
    if (fasta_source == "uniprot") {
      fasta_path <- dda_fasta_resolved() %||% ""
    } else {
      fasta_path <- trimws(input$dda_fasta_path %||% "")
    }

    # Validation
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
        fragment_tol_da   = input$dda_fragment_tol %||% 0.05,
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

      # Store state
      values$dda_job_id     <- job_id
      values$dda_output_dir <- output_dir
      values$dda_status     <- "running"

      run_casanovo <- isTRUE(input$dda_run_casanovo)
      values$dda_search_params <- list(
        preset              = input$dda_preset %||% "standard",
        fasta_path          = fasta_path,
        raw_dir             = raw_dir,
        n_files             = length(raw_paths),
        missed_cleavages    = input$dda_missed_cleavages %||% 2,
        precursor_tol       = input$dda_precursor_tol %||% 20,
        fragment_tol        = input$dda_fragment_tol %||% 0.05,
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
          "- **Fragment tolerance**: ±", sp$fragment_tol, " Da\n",
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
          casanovo_scripts <- generate_casanovo_sbatch(
            raw_dir          = raw_dir,
            output_dir       = output_dir,
            experiment_name  = exp_name,
            conda_env_path   = casanovo_conda_env,
            model_ckpt       = casanovo_model_ckpt,
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

          # Write and upload launcher script
          launcher_content <- generate_casanovo_launcher(
            remote_convert_sbatch, remote_casanovo_sbatch)
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

    if (state %in% c("COMPLETED")) {
      message("[DDA] Sage job completed: ", job_id)
      values$dda_status <- "loading"
      showNotification("Sage search completed! Loading results...",
        type = "message", duration = 8)
      # Trigger result loading
      load_sage_results_from_hpc()
    } else if (state %in% c("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED", "NODE_FAIL")) {
      message("[DDA] Sage job failed: ", job_id, " (", state, ")")
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

    if (n_pending == 0) {
      # All tasks finished
      if (n_completed > 0) {
        message("[DDA] Casanovo completed: ", n_completed, "/", n_total, " tasks")
        values$dda_casanovo_status <- "loading"
        showNotification(
          paste("Casanovo completed!", n_completed, "/", n_total, "files"),
          type = "message", duration = 8)
        # Trigger Casanovo result loading
        load_casanovo_results_from_hpc()
      } else {
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

      setProgress(0.8, detail = "Cross-referencing with Sage...")

      # Store raw Casanovo results
      values$dda_casanovo_psms <- casanovo_psms

      # Cross-reference with Sage if available
      if (!is.null(values$dda_sage_psms)) {
        classification <- tryCatch(
          classify_dda_denovo(casanovo_psms, values$dda_sage_psms),
          error = function(e) {
            message("[DDA] Classification error: ", e$message)
            NULL
          }
        )

        if (!is.null(classification)) {
          values$dda_casanovo_classification <- classification
          message(sprintf(
            "[DDA] Casanovo classification: %d confirmed, %d novel",
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
          fasta_lines <- unlist(lapply(seq_along(unique_seqs), function(i) {
            c(paste0(">denovo_", i, " ", unique_seqs[i]), unique_seqs[i])
          }))
          writeLines(fasta_lines, local_fasta)

          remote_denovo_dir <- file.path(values$dda_output_dir, "denovo")
          ssh_exec(ssh_cfg, paste("mkdir -p", shQuote(remote_denovo_dir)), timeout = 10)
          scp_upload(ssh_cfg, local_fasta, file.path(remote_denovo_dir, "novel_peptides.fasta"))

          blast_out <- file.path(remote_denovo_dir, "blast_results.tsv")
          logs_dir <- file.path(values$dda_output_dir, "logs")

          blast_sbatch <- paste0(
'#!/bin/bash
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
  --sensitive --id 50 --max-target-seqs 5 --threads 8

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
  output$dda_job_status_ui <- renderUI({
    status <- values$dda_status %||% "idle"

    switch(status,
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
        div(
          class = "alert alert-success",
          style = "margin-top: 8px; border-left: 4px solid #6f42c1;",
          icon("wand-magic-sparkles"),
          paste(" Casanovo complete!", format(n_psms, big.mark = ","), "de novo sequences"),
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

  # --- Summary cards for DDA Casanovo de novo results ---
  output$dda_denovo_summary_cards <- renderUI({
    req(values$dda_casanovo_classification)
    cls <- values$dda_casanovo_classification

    n_total     <- cls$summary_stats$n_total
    n_confirmed <- cls$summary_stats$n_confirmed
    n_novel     <- cls$summary_stats$n_novel
    pct_conf    <- cls$summary_stats$pct_confirmed

    n_proteins  <- if (!is.null(cls$protein_summary)) nrow(cls$protein_summary) else 0

    tags$div(
      class = "row",
      style = "margin-bottom: 15px;",

      tags$div(
        class = "col-md-3",
        tags$div(
          class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #3498db; padding: 15px;",
          tags$h4(format(n_total, big.mark = ","), style = "margin: 0; color: #3498db;"),
          tags$small("Total De Novo PSMs")
        )
      ),

      tags$div(
        class = "col-md-3",
        tags$div(
          class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #2ecc71; padding: 15px;",
          tags$h4(format(n_confirmed, big.mark = ","), style = "margin: 0; color: #2ecc71;"),
          tags$small("Confirmed (in Sage)")
        )
      ),

      tags$div(
        class = "col-md-3",
        tags$div(
          class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #e67e22; padding: 15px;",
          tags$h4(format(n_novel, big.mark = ","), style = "margin: 0; color: #e67e22;"),
          tags$small("Novel Peptides")
        )
      ),

      tags$div(
        class = "col-md-3",
        tags$div(
          class = "card text-center",
          style = "background: #f8f9fa; border-left: 4px solid #9b59b6; padding: 15px;",
          tags$h4(paste0(pct_conf, "%"), style = "margin: 0; color: #9b59b6;"),
          tags$small(paste0("Confirmation Rate (", n_proteins, " proteins)"))
        )
      )
    )
  })

  # --- Confirmed peptides table (DDA Casanovo) ---
  output$dda_denovo_confirmed_table <- DT::renderDT({
    req(values$dda_casanovo_classification)
    confirmed <- values$dda_casanovo_classification$confirmed
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
      caption = "Confirmed: Casanovo de novo peptides matching Sage database search results (I/L normalized)"
    )
  })

  # --- Novel peptides table (DDA Casanovo) ---
  output$dda_denovo_novel_table <- DT::renderDT({
    req(values$dda_casanovo_classification)
    novel <- values$dda_casanovo_classification$novel
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
    blast <- values$dda_casanovo_blast
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
        "Novel: Casanovo de novo peptides NOT found in Sage results.",
        tags$br(),
        tags$small(
          style = "color: #666; font-weight: normal;",
          "These may represent sequence variants, unexpected organisms, or proteins absent from your reference FASTA. ",
          "Use DIAMOND BLAST to map them to known proteins."
        )
      )
    )
  })

  # ==========================================================================
  #    BLAST Results Visualization (comprehensive SwissProt BLAST view)
  # ==========================================================================

  # Helper: ensure species + category columns exist on blast data
  blast_with_species <- reactive({
    req(values$dda_casanovo_blast)
    blast <- values$dda_casanovo_blast
    req(nrow(blast) > 0)
    # Parse species from SwissProt IDs if not already done
    if (!"species" %in% names(blast)) {
      blast$species <- sub(".*_", "", sub("^[a-z]+\\|[^|]+\\|", "", blast$subject))
    }
    if (!"category" %in% names(blast)) {
      blast$category <- ifelse(
        blast$pident >= 100, "Conserved",
        ifelse(blast$pident >= 90, "Near-match", "Distant")
      )
    }
    blast
  })

  # --- Summary cards ---
  output$dda_blast_summary_cards <- renderUI({
    blast <- blast_with_species()
    req(nrow(blast) > 0)
    novel <- values$dda_casanovo_classification$novel
    n_novel <- length(unique(novel$seq_stripped))
    n_with_hits <- length(unique(blast$peptide))
    n_no_hits <- n_novel - n_with_hits
    pct_hits <- round(100 * n_with_hits / max(n_novel, 1), 1)

    # Top species by best-hit count (deduplicate to best hit per peptide)
    best_hits <- blast[!duplicated(blast$peptide), ]
    top_sp <- names(sort(table(best_hits$species), decreasing = TRUE))[1]
    mean_id <- round(mean(blast$pident, na.rm = TRUE), 1)

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
      div(class = "col-md-3",
        div(style = "background: #f3e5f5; padding: 12px; border-radius: 8px; text-align: center;",
          tags$h4(style = "margin: 0; color: #7b1fa2;", top_sp %||% "N/A"),
          tags$small("Top species")
        )
      ),
      div(class = "col-md-3",
        div(style = "background: #fce4ec; padding: 12px; border-radius: 8px; text-align: center;",
          tags$h4(style = "margin: 0; color: #c62828;", paste0(mean_id, "%")),
          tags$small("Mean identity")
        )
      )
    )
  })

  # --- Taxonomic breakdown: donut chart ---
  output$dda_blast_species_donut <- plotly::renderPlotly({
    blast <- blast_with_species()
    # Best hit per peptide for species assignment
    best_hits <- blast[order(blast$pident, decreasing = TRUE), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide_sequence), ]

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
      )
  })

  # --- Taxonomic breakdown: bar chart ---
  output$dda_blast_species_bar <- plotly::renderPlotly({
    blast <- blast_with_species()
    best_hits <- blast[order(blast$pident, decreasing = TRUE), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide_sequence), ]

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
      )
  })

  # --- Species summary text ---
  output$dda_blast_species_summary <- renderUI({
    blast <- blast_with_species()
    best_hits <- blast[order(blast$pident, decreasing = TRUE), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide_sequence), ]

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
    blast <- blast_with_species()
    # Best hit per peptide
    best_hits <- blast[order(blast$pident, decreasing = TRUE), ]
    best_hits <- best_hits[!duplicated(best_hits$peptide_sequence), ]

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

    p <- ggplot2::ggplot(best_hits, ggplot2::aes(x = identity, fill = sp_group)) +
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
      plotly::layout(legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.05))
  })

  # --- Peptide-Species heatmap ---
  output$dda_blast_heatmap <- plotly::renderPlotly({
    blast <- blast_with_species()

    # Best hit per peptide-species combination
    blast_best <- blast[order(blast$pident, decreasing = TRUE), ]
    blast_best <- blast_best[!duplicated(paste(blast_best$peptide_sequence, blast_best$species)), ]

    # Top 10 species by frequency
    sp_counts <- sort(table(blast_best$species), decreasing = TRUE)
    top_species <- names(sp_counts)[seq_len(min(10, length(sp_counts)))]

    # Top 50 peptides by average identity
    pep_avg <- tapply(blast_best$identity, blast_best$peptide_sequence, mean, na.rm = TRUE)
    pep_avg <- sort(pep_avg, decreasing = TRUE)
    top_peptides <- names(pep_avg)[seq_len(min(50, length(pep_avg)))]

    # Build matrix
    mat <- matrix(NA_real_, nrow = length(top_peptides), ncol = length(top_species),
      dimnames = list(top_peptides, top_species))

    sub <- blast_best[blast_best$peptide_sequence %in% top_peptides &
                      blast_best$species %in% top_species, ]
    for (i in seq_len(nrow(sub))) {
      mat[sub$peptide_sequence[i], sub$species[i]] <- sub$identity[i]
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
      )
  })

  # --- Enhanced BLAST results table ---
  output$dda_denovo_blast_table <- DT::renderDT({
    blast <- blast_with_species()

    display_df <- data.frame(
      Peptide     = blast$peptide,
      Hit         = blast$subject,
      Protein     = blast$protein,
      Species     = blast$species,
      Category    = blast$category,
      Identity    = round(blast$pident, 1),
      Length      = blast$length,
      E_Value     = formatC(blast$evalue, format = "e", digits = 2),
      Bitscore    = round(blast$bitscore, 1),
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
        "DIAMOND BLAST hits for novel peptides against UniProt SwissProt (572k reviewed proteins)"
      )
    ) %>%
      DT::formatStyle("Category",
        backgroundColor = DT::styleEqual(
          c("Conserved", "Near-match", "Distant"),
          c("#e8f5e9", "#fff3e0", "#fce4ec")
        )
      )
  })

  # --- Score distribution plot (DDA Casanovo) ---
  output$dda_denovo_score_dist <- plotly::renderPlotly({
    req(values$dda_casanovo_psms)
    req(nrow(values$dda_casanovo_psms) > 0)

    df <- values$dda_casanovo_psms
    cls <- values$dda_casanovo_classification

    match_type <- if (!is.null(cls)) {
      confirmed_seqs <- cls$confirmed$seq_norm
      ifelse(df$seq_norm %in% confirmed_seqs, "Confirmed", "Novel")
    } else {
      rep("Unclassified", nrow(df))
    }

    plot_df <- data.frame(
      score = df$score,
      type  = match_type,
      stringsAsFactors = FALSE
    )

    colors <- c("Confirmed" = "#2ecc71", "Novel" = "#e67e22", "Unclassified" = "#95a5a6")

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = score, fill = type)) +
      ggplot2::geom_histogram(bins = 50, alpha = 0.8, position = "stack") +
      ggplot2::scale_fill_manual(values = colors) +
      ggplot2::labs(
        x = "Casanovo Confidence Score",
        y = "Count",
        fill = "Classification"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")

    plotly::ggplotly(p) %>%
      plotly::layout(
        legend = list(orientation = "h", x = 0.5, xanchor = "center", y = 1.05)
      )
  })

  # ============================================================================
  #    DIAMOND BLAST for DDA Novel Peptides
  # ============================================================================

  observeEvent(input$dda_run_diamond_blast, {
    req(values$dda_casanovo_classification)
    novel <- values$dda_casanovo_classification$novel
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
        query_lines <- paste0(">casanovo_", seq_along(novel_peptides), "\n", novel_peptides)
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
          names(hits) <- c("query_idx", "subject", "identity", "length",
                           "qlen", "slen", "evalue", "bitscore")

          # Extract numeric index from query name (casanovo_1, casanovo_2, ...)
          hit_idx <- as.integer(gsub("casanovo_", "", hits$query_idx))
          hits$peptide_sequence <- novel_peptides[hit_idx]

          # Extract protein accession (handles sp|ACC|NAME format)
          hits$protein <- stringr::str_extract(hits$subject, "(?<=\\|)[^|]+(?=\\|)")
          no_match <- is.na(hits$protein)
          if (any(no_match)) {
            hits$protein[no_match] <- hits$subject[no_match]
          }

          # Extract species from SwissProt ID: sp|P12345|PROT_SPECIES -> SPECIES
          hits$species <- sub(".*_", "", sub("^[a-z]+\\|[^|]+\\|", "", hits$subject))
          # Classify by identity
          hits$category <- ifelse(
            hits$identity >= 100, "Conserved",
            ifelse(hits$identity >= 90, "Near-match", "Distant")
          )

          values$dda_casanovo_blast <- as.data.frame(hits)
          n_hits <- length(unique(hits$peptide_sequence))
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

}
