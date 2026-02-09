# ==============================================================================
#  DE-LIMP: Differential Expression & Limpa Proteomics App
#  (Formerly LIMP-D)
#  Status: Production Ready (Hugging Face Compatible v1.2)
# ==============================================================================

# --- 1. AUTO-INSTALLATION & SETUP ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("BiocManager not found. Installing...")
  install.packages("BiocManager")
}

required_pkgs <- c("limma", "limpa", "ComplexHeatmap", "shinyjs", "plotly", "DT", "tidyr", "tibble", "stringr", "curl", "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db", "enrichplot", "ggridges", "ggrepel", "httr2", "markdown")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Package '", pkg, "' is missing. Auto-installing..."))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# --- 2. SERVER CONFIGURATION ---
options(repos = c(BiocManager::repositories(), CRAN = "https://cloud.r-project.org"))

library(shiny)
library(bslib)
library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(httr2)      # CRITICAL for AI Chat
library(rhandsontable)
library(DT)     
library(arrow)  
library(ComplexHeatmap)
library(shinyjs)
library(plotly)
library(stringr)
library(clusterProfiler)
library(AnnotationDbi)
library(enrichplot)
library(ggrepel)
library(markdown) # Needed for AI formatting

options(shiny.maxRequestSize = 500 * 1024^2)

# Verify Limpa
if (!requireNamespace("limpa", quietly = TRUE)) {
  stop("CRITICAL ERROR: 'limpa' failed to install. Please run: BiocManager::install('limpa') manually.")
}
library(limpa) 

# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

# --- QC Stats Calculation ---
get_diann_stats_r <- function(file_path) {
  tryCatch({
    df <- arrow::read_parquet(file_path)
    has_pg_q <- "PG.Q.Value" %in% names(df)
    has_ms1  <- "Ms1.Apex.Area" %in% names(df)
    if ("Q.Value" %in% names(df)) df <- df %>% filter(Q.Value <= 0.01)
    
    stats_df <- df %>%
      group_by(Run) %>%
      summarise(
        Precursors = n(),
        Proteins = if(has_pg_q) n_distinct(Protein.Group[PG.Q.Value <= 0.01]) else n_distinct(Protein.Group),
        MS1_Signal = if(has_ms1) sum(Ms1.Apex.Area, na.rm = TRUE) else NA_real_
      ) %>% ungroup() %>% arrange(Run) 
    
    return(stats_df)
  }, error = function(e) { data.frame(Run = "Error", Precursors = 0, Proteins = 0, MS1_Signal = 0) })
}

# --- CHECK AVAILABLE MODELS ---
list_google_models <- function(api_key) {
  req <- request("https://generativelanguage.googleapis.com/v1beta/models") %>%
    req_url_query(key = api_key)
  tryCatch({
    resp <- req_perform(req)
    data <- resp_body_json(resp)
    models <- sapply(data$models, function(x) x$name)
    models <- gsub("^models/", "", models)
    return(models)
  }, error = function(e) { return(paste("Error listing models:", e$message)) })
}

# --- FILE API UPLOADER ---
upload_csv_to_gemini <- function(df, api_key) {
  temp_file <- tempfile(fileext = ".csv")
  write.csv(df, temp_file, row.names = FALSE)
  file_size <- file.size(temp_file)
  
  req <- request("https://generativelanguage.googleapis.com/upload/v1beta/files") %>%
    req_url_query(key = api_key) %>%
    req_headers(
      "X-Goog-Upload-Protocol" = "raw",
      "X-Goog-Upload-Command" = "start, upload, finalize",
      "X-Goog-Upload-Header-Content-Length" = as.character(file_size),
      "X-Goog-Upload-Header-Content-Type" = "text/csv",
      "Content-Type" = "text/csv"
    ) %>%
    req_body_file(temp_file)
  
  resp <- req_perform(req)
  file_info <- resp_body_json(resp)
  return(file_info$file$uri)
}

# --- AI CHAT FUNCTION ---
ask_gemini_file_chat <- function(user_query, file_uri, qc_df, api_key, model_name, selected_ids = NULL) {
  
  qc_text <- "No QC Data Available"
  if(!is.null(qc_df)) {
    qc_text <- paste(capture.output(write.csv(qc_df, row.names=FALSE)), collapse="\n")
  }
  
  selection_context <- ""
  if (!is.null(selected_ids) && length(selected_ids) > 0) {
    selection_context <- paste0(
      "\n!!! URGENT: USER SELECTION ACTIVE !!!\n",
      "Focus analysis on these specific proteins:\n",
      paste(selected_ids, collapse=", "), "\n"
    )
  }
  
  system_instruction <- paste0(
    "You are a PhD-level expert in proteomics. ",
    "You have access to two data sources:\n",
    "SOURCE 1: QC METRICS (In Text Below)\n",
    "This table includes 'Group' columns. Use it to compare technical quality (Precursors, MS1) between experimental groups.\n",
    "--- START QC DATA ---\n",
    qc_text,
    "\n--- END QC DATA ---\n\n",
    "SOURCE 2: EXPRESSION DATA (In Uploaded File)\n",
    "Use this file to answer biological questions. It contains the Top 800 proteins.\n\n",
    "IMPORTANT: BI-DIRECTIONAL CONTROL.\n",
    "1. If the user asks about 'selected proteins', refer to the 'URGENT' section below.\n",
    "2. If you find interesting proteins, OUTPUT their IDs at the end like this:\n",
    "   [[SELECT: P12345; P67890]]\n"
  )
  
  base_url <- "https://generativelanguage.googleapis.com/v1beta/models/"
  clean_model <- gsub("^models/", "", model_name)
  full_url <- paste0(base_url, clean_model, ":generateContent")
  
  body <- list(contents = list(list(parts = list(
    list(text = paste0(system_instruction, selection_context, "\n\nUser Question: ", user_query)),
    list(file_data = list(file_uri = file_uri, mime_type = "text/csv"))
  ))))
  
  req <- request(full_url) %>%
    req_url_query(key = api_key) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_body_json(body)
  
  tryCatch({
    resp <- req_perform(req)
    return(resp_body_json(resp)$candidates[[1]]$content$parts[[1]]$text)
  }, error = function(e) { 
    err_msg <- "Unknown Error"
    if (!is.null(e$resp)) { err_msg <- tryCatch(resp_body_string(e$resp), error = function(z) e$message) } else { err_msg <- e$message }
    return(paste("API Error:", err_msg)) 
  })
}

# --- AI TEXT CHAT FUNCTION ---
ask_gemini_text_chat <- function(user_query, api_key, model_name) {
  base_url <- "https://generativelanguage.googleapis.com/v1beta/models/"
  clean_model <- gsub("^models/", "", model_name)
  full_url <- paste0(base_url, clean_model, ":generateContent")
  
  body <- list(contents = list(list(parts = list(list(text = user_query)))))
  
  req <- request(full_url) %>%
    req_url_query(key = api_key) %>%
    req_headers("Content-Type" = "application/json") %>%
    req_body_json(body)
  
  tryCatch({
    resp <- req_perform(req)
    return(resp_body_json(resp)$candidates[[1]]$content$parts[[1]]$text)
  }, error = function(e) { 
    err_msg <- "Unknown Error"
    if (!is.null(e$resp)) { err_msg <- tryCatch(resp_body_string(e$resp), error = function(z) e$message) } else { err_msg <- e$message }
    return(paste("API Error:", err_msg)) 
  })
}

cal_z_score <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }

# --- Auto-detect Organism ---
detect_organism_db <- function(protein_ids) {
  ORGANISM_DB_MAP <- list(
    "_HUMAN" = "org.Hs.eg.db", "_MOUSE" = "org.Mm.eg.db", "_RAT"   = "org.Rn.eg.db",
    "_BOVIN" = "org.Bt.eg.db", "_CANLF" = "org.Cf.eg.db", "_CHICK" = "org.Gg.eg.db",
    "_DROME" = "org.Dm.eg.db", "_CAEEL" = "org.Ce.eg.db", "_DANRE" = "org.Dr.eg.db",
    "_YEAST" = "org.Sc.sgd.db", "_ARATH" = "org.At.tair.db", "_PIG"   = "org.Ss.eg.db"
  )
  for (suffix in names(ORGANISM_DB_MAP)) {
    if (any(grepl(suffix, protein_ids, ignore.case = TRUE))) {
      return(ORGANISM_DB_MAP[[suffix]])
    }
  }
  return("org.Hs.eg.db")
}

# ==============================================================================
#  USER INTERFACE (UI)
# ==============================================================================

ui <- page_sidebar(
  title = "DE-LIMP Proteomics",
  theme = bs_theme(bootswatch = "flatly"),
  useShinyjs(), 
  
  tags$head(tags$style(HTML(" 
    .chat-container { height: 500px; overflow-y: auto; border: 1px solid #ddd; padding: 15px; background-color: #f8f9fa; border-radius: 5px; margin-bottom: 15px; }
    .user-msg { text-align: right; margin: 10px 0; }
    .user-msg span { background-color: #007bff; color: white; padding: 8px 12px; border-radius: 15px 15px 0 15px; display: inline-block; max-width: 80%; }
    .ai-msg { text-align: left; margin: 10px 0; }
    .ai-msg span { background-color: #e9ecef; color: #333; padding: 8px 12px; border-radius: 15px 15px 15px 0; display: inline-block; max-width: 80%; }
    .selection-banner { background-color: #d4edda; color: #155724; padding: 10px; border-radius: 5px; margin-bottom: 10px; font-weight: bold; border: 1px solid #c3e6cb; }
  "))),
  
  tags$head(tags$script(HTML(" 
    $(document).on('shown.bs.modal', function() { setTimeout(function() { $(window).trigger('resize'); }, 200); });
  "))),
  
  sidebar = sidebar(
    width = 320,
    title = "Controls",
    
    # --- MOVED TO TOP: DOCS & REPO LINKS ---
    div(style="display: flex; gap: 5px; margin-bottom: 15px;",
      tags$a(href="https://github.com/bsphinney/DE-LIMP/blob/main/USER_GUIDE.md", target="_blank", class="btn btn-info w-50", icon("book"), "Guide", style="color:white; font-weight:bold;"),
      tags$a(href="https://github.com/bsphinney/DE-LIMP", target="_blank", class="btn btn-secondary w-50", icon("github"), "Code", style="color:white; font-weight:bold;")
    ),
    hr(),
    
    h5("1. Upload"),
    fileInput("report_file", "DIA-NN Report (.parquet)", accept = c(".parquet")),
    numericInput("q_cutoff", "Q-Value Cutoff", value = 0.01, min = 0, max = 0.1, step = 0.01),
    hr(),
    h5("2. Setup"),
    actionButton("open_setup", "Assign Groups", class = "btn-info w-100", icon = icon("table")),
    hr(),
    h5("3. Analyze"),
    actionButton("run", "Run Pipeline", class = "btn-primary w-100", icon = icon("play")),
    textOutput("run_status_msg"),
    hr(),
    h5("4. Explore Results"),
    selectInput("contrast_selector", "Comparison:", choices=NULL, width="100%"),
    sliderInput("logfc_cutoff", "Min Log2 Fold Change:", min=0, max=5, value=1, step=0.1),
    hr(),
    h5("5. AI Chat"),
    passwordInput("user_api_key", "Gemini API Key", value = "", placeholder = "AIzaSy..."),
    actionButton("check_models", "Check Models", class="btn-warning btn-xs w-100"),
    br(), br(),
    textInput("model_name", "Model Name", value = "gemini-3-flash-preview", placeholder = "gemini-3-flash-preview")
  ),
  
  navset_card_tab(
    id = "main_tabs", 
    
    nav_panel("Data Overview", icon = icon("database"),
              div(style="display: flex; gap: 10px; margin-bottom: 10px;",
                  actionButton("show_summary_modal", "View Full Summary", icon = icon("list-alt"), class = "btn-primary w-50"),
                  actionButton("show_grid_view", "Open Grid View", icon = icon("th"), class = "btn-success w-50")
              ),
              card(
                card_header("Signal Distribution Across All Protein Groups"),
                card_body(
                  div(
                    actionButton("color_de", "Color by DE Status", icon = icon("paint-brush"), class = "btn-info btn-sm"),
                    actionButton("reset_color", "Reset Colors", icon = icon("undo"), class = "btn-secondary btn-sm")
                  ),
                  hr(),
                  plotOutput("protein_signal_plot", height = "500px")
                )
              ),
              card(
                card_header("Group QC Summary"),
                card_body(DTOutput("group_summary_table"))
              )
    ),
    
    nav_panel("QC Trends", icon = icon("chart-bar"),
              layout_columns(col_widths=c(12,12),
                             card(
                               card_header("Trend Analysis"),
                               card_body(
                                 selectInput("qc_metric_select", "Metric:", choices = c("Precursors", "Proteins", "MS1_Signal")),
                                 radioButtons("qc_sort_order", "Order By:", choices = c("Run Order", "Group"), inline = TRUE),
                                 plotlyOutput("qc_trend_plot", height = "500px")
                               )
                             ),
                             card(card_header("Stats Table"), DTOutput("r_qc_table"))
              )
    ),
    
    nav_panel("QC Plots", icon = icon("chart-scatter"),
              layout_columns(col_widths=c(6,6), 
                             card(card_header("DPC Fit"), plotOutput("dpc_plot", height="400px")),
                             card(card_header("MDS Plot"), plotOutput("mds_plot", height="400px"))),
              layout_columns(col_widths=c(12),
                             card(
                               card_header("Group QC Distribution (Hover for Info)"),
                               card_body(
                                 selectInput("qc_violin_metric", "Metric:", choices = c("Precursors", "Proteins", "MS1_Signal"), width = "200px"),
                                 plotlyOutput("qc_group_violin", height = "400px")
                               )
                             )
              )
    ),
    
    nav_panel("DE Dashboard", icon = icon("table-columns"),
              layout_columns(col_widths = c(6, 6),
                             card(card_header(div(style="display: flex; justify-content: space-between; align-items: center;", span("Results Table"), div(actionButton("generate_ai_summary", "🤖 Generate AI Summary", class="btn-info btn-sm"), actionButton("clear_plot_selection", "Reset", class="btn-warning btn-xs"), actionButton("show_violin", "📊 Violin Plot", class="btn-primary btn-xs"), downloadButton("download_result_csv", "💾 Export Results", class="btn-success btn-xs")))), DTOutput("de_table")),
                             card(card_header("Volcano Plot (Click/Box Select to Filter Table)"), plotlyOutput("volcano_plot_interactive", height = "600px"))),
              card(card_header("Heatmap"), plotOutput("heatmap_plot", height="400px"))),
    
    nav_panel("Consistent DE", icon = icon("check-double"),
              card(
                card_header("High-Consistency Significant Proteins (Ranked by %CV)"),
                card_body(
                  p("Ranking by %CV (Coefficient of Variation) to find stable markers."),
                  DTOutput("consistent_table")
                )
              )
    ),
    
    nav_panel("Reproducibility", icon = icon("code"),
              navset_card_tab(
                nav_panel("Code Log",
                          card_body(
                            p("Copy the R code below to reproduce your analysis steps. It's timestamped and includes comments."),
                            verbatimTextOutput("reproducible_code")
                          )
                ),
                nav_panel("Methodology",
                          card_body(
                            textOutput("methodology_text")
                          )
                )
              )
    ),
    
    nav_panel("Gene Set Enrichment", icon = icon("sitemap"),
              card(
                card_header("Gene Ontology (GO) Analysis"),
                card_body(
                  actionButton("run_gsea", "Run GSEA", class = "btn-success w-100", icon = icon("play")),
                  verbatimTextOutput("gsea_status"),
                  hr(),
                  p("This panel performs Gene Set Enrichment Analysis on the ranked list of proteins from the DE results. It automatically detects the organism (Human or Mouse) to use the correct annotation database.", class = "text-muted small")
                )
              ),
              card(
                card_header("GSEA Results"),
                navset_card_tab(
                  nav_panel("Dot Plot", plotOutput("gsea_dot_plot", height = "500px")),
                  nav_panel("Enrichment Map", plotOutput("gsea_emapplot", height = "500px")),
                  nav_panel("Ridgeplot", plotOutput("gsea_ridgeplot", height = "600px")),
                  nav_panel("Results Table", DTOutput("gsea_results_table"))
                )
              )
    ),
    
    nav_panel("Data Chat", icon = icon("comments"),
              card(
                card_header(div(style="display: flex; justify-content: space-between; align-items: center;", span("Chat with Full Data (QC + Expression)"), downloadButton("download_chat_txt", "💾 Save Chat", class="btn-secondary btn-sm"))),
                card_body(
                  verbatimTextOutput("chat_selection_indicator"),
                  uiOutput("chat_window"),
                  tags$div(style="margin-top: 15px; display: flex; gap: 10px;",
                           textAreaInput("chat_input", label=NULL, placeholder="Ask e.g. 'Which group has higher precursor counts?'", width="100%", rows=2),
                           actionButton("summarize_data", "🤖 Auto-Analyze", class="btn-info", style="height: 54px; margin-top: 2px;"),
                           actionButton("send_chat", "Send", icon=icon("paper-plane"), class="btn-primary", style="height: 54px; margin-top: 2px;")
                  ),
                  p("Note: QC Stats (with Groups) + Top 800 Expression Data are sent to AI.", style="font-size: 0.8em; color: green; font-weight: bold; margin-top: 5px;")
                )
              )
    )
  )
)

# ==============================================================================
#  SERVER LOGIC
# ==============================================================================
server <- function(input, output, session) { 
  
  # ============================================================================
  #      1. Reactive Values & Global Helpers
  # ============================================================================
  values <- reactiveValues(
    raw_data = NULL, metadata = NULL, fit = NULL, y_protein = NULL, 
    dpc_fit = NULL, status = "Waiting...", design = NULL, qc_stats = NULL,
    plot_selected_proteins = NULL, chat_history = list(),
    current_file_uri = NULL, gsea_results = NULL,
    repro_log = character(0),
    color_plot_by_de = FALSE,
    grid_selected_protein = NULL,
    temp_violin_target = NULL # Added for violin popup
  )
  
  # ============================================================================
  #      2. Main Data Loading & Processing Pipeline
  # ============================================================================
  observeEvent(input$report_file, {
    req(input$report_file)
    withProgress(message = "Loading...", {
      incProgress(0.2, detail = "Calculating Trends...")
      values$qc_stats <- get_diann_stats_r(input$report_file$datapath)
      incProgress(0.5, detail = "Reading Matrix...")
      tryCatch({
        values$raw_data <- limpa::readDIANN(input$report_file$datapath, format="parquet", q.cutoffs=input$q_cutoff)
        fnames <- sort(colnames(values$raw_data$E))
        values$metadata <- data.frame(ID = 1:length(fnames), File.Name = fnames, Group = rep("", length(fnames)), stringsAsFactors=FALSE)
        click("open_setup") 
      }, error=function(e) { showNotification(paste("Error:", e$message), type="error") })
    })
  })

  observeEvent(input$run, {
    req(values$raw_data, values$metadata); meta <- values$metadata; meta$Group <- trimws(meta$Group)
    if(length(unique(meta$Group)) < 2) { showNotification("Error: Need at least 2 groups.", type="error"); return() }

    # --- Start Reproducibility Log ---
    base_script <- c(
      "# === DE-LIMP Reproducibility Log ===",
      sprintf("# Log generated on %s", Sys.time()),
      "\n# --- 1. Load Libraries ---",
      "library(limpa); library(limma); library(dplyr); library(stringr); library(ggrepel);",
      "\n# --- 2. Data Loading ---",
      sprintf("dat <- readDIANN('path/to/your/report.parquet', format='parquet', q.cutoffs=%s)", input$q_cutoff),
      "\n# --- 3. Normalization & Quantification ---",
      "dpcfit <- dpcCN(dat)",
      "y_protein <- dpcQuant(dat, 'Protein.Group', dpc=dpcfit)",
      "\n# --- 4. Experimental Design ---",
      "group_map <- c(", paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Group), collapse=",\n"), ")",
      "metadata <- data.frame(File.Name = names(group_map), Group = group_map)",
      "metadata <- metadata[match(colnames(dat$E), metadata$File.Name), ]",
      "groups <- factor(metadata$Group)",
      "design <- model.matrix(~ 0 + groups); colnames(design) <- levels(groups)",
      "\n# --- 5. Initial Model Fit ---",
      "fit <- dpcDE(y_protein, design, plot=FALSE)"
    )
    
    values$repro_log <- base_script 
    
    withProgress(message='Running Pipeline...', {
      tryCatch({
        dat <- values$raw_data; 
        dpcfit <- limpa::dpcCN(dat);
        values$dpc_fit <- dpcfit;
        
        values$y_protein <- tryCatch({
          limpa::dpcQuant(dat, "Protein.Group", dpc=dpcfit)
        }, error = function(e) {
          showNotification(paste("Protein quantification failed:", e$message), type = "error", duration = NULL)
          return(NULL)
        })
        
        req(values$y_protein)
        
        rownames(meta) <- meta$File.Name; meta <- meta[colnames(dat$E), ]; meta$Group <- make.names(meta$Group); groups <- factor(meta$Group)
        design <- model.matrix(~ 0 + groups); colnames(design) <- levels(groups)
        combs <- combn(levels(groups), 2); forms <- apply(combs, 2, function(x) paste(x[2], "-", x[1]))
        fit <- limpa::dpcDE(values$y_protein, design, plot=FALSE); fit <- contrasts.fit(fit, makeContrasts(contrasts=forms, levels=design)); fit <- eBayes(fit)
        values$fit <- fit; updateSelectInput(session, "contrast_selector", choices=forms); values$status <- "✅ Complete!"; nav_select("main_tabs", "QC Plots")
      }, error=function(e) showNotification(paste("Error during DE analysis:", e$message), type="error", duration = NULL))
    })
  })
  
  # ============================================================================
  #      3. Setup Modal & Metadata Handling
  # ============================================================================
  observeEvent(input$open_setup, { 
    req(values$metadata)
    showModal(modalDialog(
      title = "Assign Groups", size = "xl",
      actionButton("guess_groups", "🪄 Auto-Guess", class="btn-info btn-sm"), br(), br(),
      rHandsontableOutput("hot_metadata_modal"),
      footer = tagList(modalButton("Cancel"), actionButton("save_groups", "Save & Close", class="btn-success"))
    ))
  })

  output$hot_metadata_modal <- renderRHandsontable({ 
    req(values$metadata)
    rhandsontable(values$metadata, rowHeaders=NULL, stretchH="all", height=500, width="100%") %>%
      hot_col("ID", readOnly=TRUE, width=50) %>% hot_col("File.Name", readOnly=TRUE) %>% hot_col("Group", type="text") 
  })

  observeEvent(input$guess_groups, { 
    req(values$metadata)
    meta <- if(!is.null(input$hot_metadata_modal)) hot_to_r(input$hot_metadata_modal) else values$metadata
    
    cleaned_filenames <- str_remove(meta$File.Name, "^\\d{8}_")
    cleaned_filenames <- str_remove_all(cleaned_filenames, "_deepfrac|\\.parquet")
    keywords <- c("affinisepACN", "affinisepIPA", "Control", "Treatment", "Evosep", "Affinisep")
    
    find_best_match <- function(filename_clean) {
      matches <- keywords[str_detect(filename_clean, regex(keywords, ignore_case = TRUE))]
      if (length(matches) == 0) return("") else return(matches[which.max(nchar(matches))])
    }
    
    guessed_groups <- sapply(cleaned_filenames, find_best_match)
    sample_counter <- 0
    for (i in seq_along(guessed_groups)) {
      if (guessed_groups[i] == "") {
        sample_counter <- sample_counter + 1
        guessed_groups[i] <- paste0("Sample_", sample_counter)
      }
    }
    meta$Group <- guessed_groups
    values$metadata <- meta
  })

  observeEvent(input$save_groups, { 
    req(input$hot_metadata_modal, values$metadata)
    old_meta <- values$metadata
    new_meta <- hot_to_r(input$hot_metadata_modal)
    changed_indices <- which(old_meta$Group != new_meta$Group)
    if (length(changed_indices) > 0) {
      code_lines <- sprintf("metadata$Group[%d] <- '%s'", changed_indices, new_meta$Group[changed_indices])
      values$repro_log <- c(values$repro_log, "", c(sprintf("# --- Manual Group Edits [%s] ---", Sys.time()), code_lines))
    }
    values$metadata <- new_meta
    removeModal()
    showNotification("Groups saved!", type="message")
  })

  # ============================================================================
  #      4. Reactive Data Expressions
  # ============================================================================
  volcano_data <- reactive({
    req(values$fit, input$contrast_selector)
    df_raw <- topTable(values$fit, coef=input$contrast_selector, number=Inf) %>% as.data.frame()
    if (!"Protein.Group" %in% colnames(df_raw)) { df <- df_raw %>% rownames_to_column("Protein.Group") } else { df <- df_raw }
    
    org_db_name <- detect_organism_db(df$Protein.Group)
    df$Accession <- str_split_fixed(df$Protein.Group, "[; ]", 2)[,1]
    
    id_map <- tryCatch({
      if (!requireNamespace(org_db_name, quietly = TRUE)) { BiocManager::install(org_db_name, ask = FALSE) }
      library(org_db_name, character.only = TRUE)
      bitr(df$Accession, fromType = "UNIPROT", toType = c("SYMBOL", "GENENAME"), OrgDb = get(org_db_name))
    }, error = function(e) { NULL })
    
    if (!is.null(id_map)) {
      id_map <- id_map %>% distinct(UNIPROT, .keep_all = TRUE)
      df <- df %>% left_join(id_map, by = c("Accession" = "UNIPROT")) %>%
        mutate(Gene = ifelse(is.na(SYMBOL), Accession, SYMBOL), Protein.Name = ifelse(is.na(GENENAME), Protein.Group, GENENAME))
    } else {
      df$Gene <- df$Accession; df$Protein.Name <- df$Protein.Group
    }
    
    df$Significance <- "Not Sig"; df$Significance[df$adj.P.Val < 0.05 & abs(df$logFC) > input$logfc_cutoff] <- "Significant"
    df$Selected <- "No"; if (!is.null(values$plot_selected_proteins)) { df$Selected[df$Protein.Group %in% values$plot_selected_proteins] <- "Yes" }
    df
  })

  # ============================================================================
  #      5. UI Outputs & Renderers
  # ============================================================================
  output$run_status_msg <- renderText({ values$status })
  
  observeEvent(input$show_summary_modal, {
    req(values$metadata)
    summary_elements <- list()
    summary_elements[[length(summary_elements) + 1]] <- tags$h5("File Summary")
    summary_elements[[length(summary_elements) + 1]] <- tags$p(paste("Total Files:", nrow(values$metadata)))
    summary_elements[[length(summary_elements) + 1]] <- tags$p(paste("Assigned Groups:", length(unique(values$metadata$Group[values$metadata$Group != ""])) ))
    
    if (!is.null(values$y_protein)) {
      summary_elements[[length(summary_elements) + 1]] <- tags$hr()
      summary_elements[[length(summary_elements) + 1]] <- tags$h5("Dataset Metrics")
      avg_signal <- rowMeans(values$y_protein$E, na.rm = TRUE)
      min_linear <- 2^min(avg_signal, na.rm = TRUE)
      max_linear <- 2^max(avg_signal, na.rm = TRUE)
      if (min_linear > 1e-10) {
        orders_of_magnitude <- log10(max_linear / min_linear)
        summary_elements[[length(summary_elements) + 1]] <- tags$p(paste("Signal Dynamic Range:", round(orders_of_magnitude, 1), "orders of magnitude"))
      } else {
        summary_elements[[length(summary_elements) + 1]] <- tags$p("Dynamic Range: N/A (Min signal is zero)")
      }
    }
    showModal(modalDialog(title = "Dataset Summary", tagList(summary_elements), easyClose = TRUE, footer = modalButton("Close")))
  })

  # --- GRID VIEW & PLOT LOGIC ---
  grid_react_df <- reactive({
    req(values$fit, values$y_protein, values$metadata, input$contrast_selector)
    df_volc <- volcano_data()
    df_exprs <- as.data.frame(values$y_protein$E) %>% rownames_to_column("Protein.Group")
    df_merged <- left_join(df_volc, df_exprs, by="Protein.Group")
    if (!is.null(values$plot_selected_proteins)) { df_merged <- df_merged %>% filter(Protein.Group %in% values$plot_selected_proteins) }
    
    meta_sorted <- values$metadata %>% arrange(Group, File.Name)
    ordered_files <- meta_sorted$File.Name
    valid_cols <- intersect(ordered_files, colnames(df_exprs))
    run_ids <- values$metadata$ID[match(valid_cols, values$metadata$File.Name)]
    new_headers <- as.character(run_ids)
    
    df_final <- df_merged %>%
      dplyr::select(Protein.Group, Gene, Protein.Name, Significance, logFC, adj.P.Val, all_of(valid_cols)) %>%
      mutate(across(where(is.numeric), ~round(., 2)))
    
    df_final$Original.ID <- df_final$Protein.Group
    fixed_cols <- c("Protein.Group", "Gene", "Protein.Name", "Significance", "logFC", "adj.P.Val")
    colnames(df_final) <- c(fixed_cols, new_headers, "Original.ID")
    
    unique_groups <- sort(unique(meta_sorted$Group))
    group_colors <- setNames(rainbow(length(unique_groups), v=0.85, s=0.8), unique_groups)
    
    list(data = df_final, fixed_cols = fixed_cols, expr_cols = new_headers, valid_cols_map = valid_cols, meta_sorted = meta_sorted, group_colors = group_colors)
  })

  observeEvent(input$show_grid_view, {
    gdata <- grid_react_df()
    legend_ui <- tags$div(style = "margin-bottom: 10px;", tags$strong("Condition Legend:"), tags$br(),
      lapply(names(gdata$group_colors), function(grp) { tags$span(style = paste0("background-color:", gdata$group_colors[[grp]], "; color:white; padding:2px 6px; margin-right:5px; border-radius:3px;"), grp) })
    )
    file_map_ui <- tags$details(
      tags$summary(tags$strong("Click to view File ID Mapping (Run # -> Filename)")),
      tags$div(style = "max-height: 150px; overflow-y: auto; background: #f9f9f9; padding: 10px; border: 1px solid #eee; margin-top: 5px;",
        lapply(1:nrow(gdata$meta_sorted), function(i) {
          row <- gdata$meta_sorted[i, ]; tags$div(tags$span(style="font-weight:bold; color:#007bff;", paste0("[", row$ID, "] ")), tags$span(row$File.Name), tags$span(style="color:gray; font-size:0.9em;", paste0(" (", row$Group, ")")))
        })
      )
    )

    showModal(modalDialog(
      title = "Expression Grid View", size = "xl", legend_ui, file_map_ui, hr(),
      DTOutput("grid_view_table"),
      footer = tagList(actionButton("grid_reset_selection", "Show All / Clear Selection", class="btn-warning"), downloadButton("download_grid_data", "Export Full Table (with Filenames)", class="btn-success"), modalButton("Close")), easyClose = TRUE
    ))
  })
  
  observeEvent(input$grid_reset_selection, { values$plot_selected_proteins <- NULL })
  
  output$grid_view_table <- renderDT({
    gdata <- grid_react_df(); df_display <- gdata$data
    acc <- str_split_fixed(df_display$Original.ID, "[; ]", 2)[,1]
    df_display$Protein.Group <- paste0("<a href='https://www.uniprot.org/uniprotkb/", acc, "/entry' target='_blank'>", df_display$Protein.Group, "</a>")
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
    datatable(df_display, container = header_html, selection = 'single', escape = FALSE, options = list(dom = 'frtip', pageLength = 15, scrollX = TRUE, columnDefs = list(list(className = 'dt-center', targets = (length(fixed_cols)):(ncol(df_display)-1)))), rownames = FALSE) %>% formatStyle(expr_cols, backgroundColor = styleInterval(brks, clrs))
  })
  
  output$download_grid_data <- downloadHandler(
    filename = function() { paste0("DE_LIMP_Grid_Export_", Sys.Date(), ".csv") },
    content = function(file) {
      gdata <- grid_react_df(); df_export <- gdata$data %>% dplyr::select(-Original.ID)
      fixed_cols <- gdata$fixed_cols; real_names <- gdata$valid_cols_map
      colnames(df_export) <- c(fixed_cols, real_names)
      write.csv(df_export, file, row.names = FALSE)
    }
  )
  
  observeEvent(input$grid_view_table_rows_selected, {
    req(grid_react_df()); selected_idx <- input$grid_view_table_rows_selected
    if (length(selected_idx) > 0) {
      gdata <- grid_react_df(); selected_id <- gdata$data$Original.ID[selected_idx]; values$grid_selected_protein <- selected_id
      showModal(modalDialog(title = paste("Expression Plot:", selected_id), size = "xl", plotOutput("violin_plot_grid", height = "600px"), footer = tagList(actionButton("back_to_grid", "Back to Grid", class="btn-info"), modalButton("Close")), easyClose = TRUE))
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

  observeEvent(input$color_de, { values$color_plot_by_de <- TRUE })
  observeEvent(input$reset_color, { values$color_plot_by_de <- FALSE })

  output$protein_signal_plot <- renderPlot({
    req(values$y_protein)
    avg_signal <- rowMeans(values$y_protein$E, na.rm = TRUE)
    plot_df <- data.frame(Protein.Group = names(avg_signal), Average_Signal_Log2 = avg_signal) %>% mutate(Average_Signal_Log10 = Average_Signal_Log2 / log2(10))
    
    if (values$color_plot_by_de && !is.null(values$fit) && !is.null(input$contrast_selector) && nchar(input$contrast_selector) > 0) {
      de_data_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(de_data_raw)) { de_data_intermediate <- de_data_raw %>% rownames_to_column("Protein.Group") } else { de_data_intermediate <- de_data_raw }
      de_data <- de_data_intermediate %>% mutate(DE_Status = case_when(adj.P.Val < 0.05 & logFC > input$logfc_cutoff ~ "Up-regulated", adj.P.Val < 0.05 & logFC < -input$logfc_cutoff ~ "Down-regulated", TRUE ~ "Not Significant")) %>% dplyr::select(Protein.Group, DE_Status)
      plot_df <- left_join(plot_df, de_data, by = "Protein.Group"); plot_df$DE_Status[is.na(plot_df$DE_Status)] <- "Not Significant"
    } else { plot_df$DE_Status <- "Not Significant" }
    
    if (!is.null(values$plot_selected_proteins)) { plot_df$Is_Selected <- plot_df$Protein.Group %in% values$plot_selected_proteins } else { plot_df$Is_Selected <- FALSE }
    selected_df <- filter(plot_df, Is_Selected)
    
    p <- ggplot(plot_df, aes(x = reorder(Protein.Group, -Average_Signal_Log10), y = Average_Signal_Log10))
    if (values$color_plot_by_de && !is.null(values$fit)) {
      p <- p + geom_point(aes(color = DE_Status), size = 1.5) + scale_color_manual(name = "DE Status", values = c("Up-regulated" = "#e41a1c", "Down-regulated" = "#377eb8", "Not Significant" = "grey70"))
    } else { p <- p + geom_point(color = "cornflowerblue", size = 1.5) }
    
    p + labs(title = "Signal Distribution Across All Protein Groups", x = NULL, y = "Average Signal (Log10 Intensity)") +
      theme_minimal() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      scale_x_discrete(expand = expansion(add = 1)) +
      geom_point(data = selected_df, color = "black", shape = 1, size = 4, stroke = 1) +
      geom_text_repel(data = selected_df, aes(label = Protein.Group), size = 4, max.overlaps = 20)
  }, height = 500) # FIXED HEIGHT

  output$group_summary_table <- renderDT({
    req(values$qc_stats, values$metadata)
    summary_df <- values$qc_stats %>% left_join(values$metadata, by = c("Run" = "File.Name")) %>% filter(Group != "") %>% group_by(Group) %>% summarise(`Avg. Precursors` = mean(Precursors, na.rm = TRUE), `Avg. MS1 Signal` = mean(MS1_Signal, na.rm = TRUE), `Avg. Proteins` = mean(Proteins, na.rm = TRUE)) %>% mutate(across(where(is.numeric), ~round(.x, 0)))
    datatable(summary_df, options = list(dom = 't', pageLength = 10), rownames = FALSE)
  })
  
  output$qc_trend_plot <- renderPlotly({
    req(values$qc_stats, input$qc_metric_select, values$metadata)
    df <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>% mutate(Run_Number = as.numeric(str_extract(Run, "\\d+$")))
    if (input$qc_sort_order == "Group") { df <- df %>% arrange(Group, Run_Number) } else { df <- df %>% arrange(Run_Number) }
    df$Sort_Index <- factor(1:nrow(df), levels = 1:nrow(df)); metric <- input$qc_metric_select
    df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Group:</b> ", df$Group, "<br><b>", metric, ":</b> ", round(df[[metric]], 2))
    p <- ggplot(df, aes(x = Sort_Index, y = .data[[metric]], fill = Group, text = Tooltip)) + geom_bar(stat = "identity", width = 0.8) + theme_minimal() + labs(title = paste(metric, "per Run"), x = "Sample Index (Sorted)", y = metric) + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size=8))
    ggplotly(p, tooltip = "text") %>% config(displayModeBar = FALSE)
  })
  
  output$r_qc_table <- renderDT({ req(values$qc_stats); df_display <- values$qc_stats %>% arrange(Run) %>% mutate(ID = 1:n()) %>% dplyr::select(ID, Run, everything()); datatable(df_display, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) })
  
  output$qc_group_violin <- renderPlotly({
    req(values$qc_stats, values$metadata, input$qc_violin_metric)
    df <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")); metric <- input$qc_violin_metric
    df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Val:</b> ", round(df[[metric]], 2))
    p <- ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) + geom_violin(alpha = 0.5, trim = FALSE) + geom_jitter(aes(text = Tooltip), width = 0.2, size = 2, alpha = 0.8, color = "black") + theme_bw() + labs(title = paste("Distribution of", metric), x = "Group", y = metric) + theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })
  
  output$dpc_plot <- renderPlot({ req(values$dpc_fit); limpa::plotDPC(values$dpc_fit) }, height = 400) # FIXED HEIGHT
  
  output$mds_plot <- renderPlot({
    req(values$y_protein, values$metadata)
    meta <- values$metadata[match(colnames(values$y_protein$E), values$metadata$File.Name), ]; grps <- factor(meta$Group); cols <- rainbow(length(levels(grps)))
    par(xpd = TRUE); limpa::plotMDSUsingSEs(values$y_protein, pch=16, main="MDS Plot", col=cols[grps]); legend(x = "right", inset = c(-0.2, 0), legend=levels(grps), col=cols, pch=16, bty = "n")
  }, height = 400) # FIXED HEIGHT

  output$de_table <- renderDT({
    df_full <- volcano_data(); df_filtered <- if (!is.null(values$plot_selected_proteins)) { df_full %>% filter(Protein.Group %in% values$plot_selected_proteins) } else { df_full }
    df_display <- df_filtered %>% mutate(across(where(is.numeric), function(x) round(x,4))) %>%
      mutate(Protein.Name_Link = ifelse(!is.na(Accession) & str_detect(Accession, "^[A-Z0-9]{6,}$"), paste0("<a href='https://www.uniprot.org/uniprotkb/", Accession, "/entry' target='_blank' onclick='window.open(this.href, \"_blank\"); return false;'>", Protein.Name, "</a>"), Protein.Name)) %>%
      dplyr::select(Gene, `Protein Name` = Protein.Name_Link, logFC, adj.P.Val, Significance)
    datatable(df_display, selection = "multiple", options = list(pageLength = 10, scrollX = TRUE), escape = FALSE, rownames = FALSE)
  })

  observeEvent(input$generate_ai_summary, {
    req(values$fit, input$contrast_selector, values$y_protein, input$user_api_key)
    withProgress(message = "Generating AI Summary...", value = 0, {
      incProgress(0.2, detail = "Gathering DE data...")
      top_de_data <- volcano_data() %>% filter(Significance != "Not Significant") %>% arrange(adj.P.Val) %>% head(50) %>% dplyr::select(Gene, logFC, adj.P.Val) %>% mutate(across(where(is.numeric), ~round(.x, 3)))
      top_de_text <- paste(capture.output(print(as.data.frame(top_de_data))), collapse = "\n")
      
      incProgress(0.5, detail = "Gathering stable proteins...")
      stable_prots_df <- tryCatch({
        df_res_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>% as.data.frame() %>% filter(adj.P.Val < 0.05)
        if (!"Protein.Group" %in% colnames(df_res_raw)) { df_res <- df_res_raw %>% rownames_to_column("Protein.Group") } else { df_res <- df_res_raw }
        if(nrow(df_res) == 0) return(data.frame(Info="No significant proteins to assess for stability."))
        protein_ids_for_cv <- df_res$Protein.Group; raw_exprs <- values$y_protein$E[protein_ids_for_cv, , drop = FALSE]; linear_exprs <- 2^raw_exprs; cv_list <- list()
        for(g in unique(values$metadata$Group)) {
          if (g == "") next
          files_in_group <- values$metadata$File.Name[values$metadata$Group == g]; group_cols <- intersect(colnames(linear_exprs), files_in_group)
          if (length(group_cols) > 1) { group_data <- linear_exprs[, group_cols, drop = FALSE]; cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
          } else { cv_list[[paste0("CV_", g)]] <- NA }
        }
        cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
        left_join(df_res, cv_df, by = "Protein.Group") %>% rowwise() %>% mutate(Avg_CV = mean(c_across(starts_with("CV_")), na.rm = TRUE)) %>% ungroup() %>% arrange(Avg_CV) %>% head(3) %>% dplyr::select(Protein.Group, Avg_CV, logFC, adj.P.Val) %>% mutate(across(where(is.numeric), ~round(.x, 2)))
      }, error = function(e) data.frame(Error = "Could not calculate stable proteins."))
      stable_prots_text <- paste(capture.output(print(as.data.frame(stable_prots_df))), collapse = "\n")

      incProgress(0.7, detail = "Constructing prompt...")
      system_prompt <- "You are a senior proteomics consultant at a core facility. Write a 3-paragraph summary..."
      final_prompt <- paste0(system_prompt, "\n\n--- DATA FOR SUMMARY ---\n\n...", top_de_text, "...", stable_prots_text, "...")
      incProgress(0.8, detail = "Asking AI...")
      ai_summary <- ask_gemini_text_chat(final_prompt, input$user_api_key, input$model_name)
      showModal(modalDialog(title = "AI-Generated Summary", HTML(markdown::markdownToHTML(text = ai_summary, fragment.only = TRUE)), easyClose = TRUE, footer = modalButton("Close")))
    })
  })
  
  output$heatmap_plot <- renderPlot({
    req(values$fit, values$y_protein, input$contrast_selector)
    df_volc <- volcano_data(); prot_ids <- NULL
    if (!is.null(input$de_table_rows_selected)) { current_table_data <- df_volc; if (!is.null(values$plot_selected_proteins)) current_table_data <- current_table_data %>% filter(Protein.Group %in% values$plot_selected_proteins); prot_ids <- current_table_data$Protein.Group[input$de_table_rows_selected]
    } else if (!is.null(values$plot_selected_proteins)) { prot_ids <- values$plot_selected_proteins; if(length(prot_ids) > 50) prot_ids <- head(prot_ids, 50)
    } else { top_prots <- topTable(values$fit, coef=input$contrast_selector, number=20); prot_ids <- rownames(top_prots) }
    valid_ids <- intersect(prot_ids, rownames(values$y_protein$E)); if (length(valid_ids) == 0) return(NULL)
    mat <- values$y_protein$E[valid_ids, , drop=FALSE]; mat_z <- t(apply(mat, 1, cal_z_score)); meta <- values$metadata[match(colnames(mat), values$metadata$File.Name), ]; groups <- factor(meta$Group)
    ha <- HeatmapAnnotation(Group = groups, col = list(Group = setNames(rainbow(length(levels(groups))), levels(groups)))); Heatmap(mat_z, name="Z-score", top_annotation = ha, cluster_rows=TRUE, cluster_columns=TRUE, show_column_names=FALSE)
  }, height = 400) # FIXED HEIGHT
  
  output$consistent_table <- renderDT({
    req(values$fit, values$y_protein, input$contrast_selector, values$metadata)
    df_res_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>% as.data.frame() %>% filter(adj.P.Val < 0.05)
    if (!"Protein.Group" %in% colnames(df_res_raw)) { df_res <- df_res_raw %>% rownames_to_column("Protein.Group") } else { df_res <- df_res_raw }
    if(nrow(df_res) == 0) return(datatable(data.frame(Status="No significant proteins found.")))
    protein_ids_for_cv <- df_res$Protein.Group; raw_exprs <- values$y_protein$E[protein_ids_for_cv, , drop = FALSE]; linear_exprs <- 2^raw_exprs; cv_list <- list()
    for(g in unique(values$metadata$Group)) {
      if (g == "") next
      files_in_group <- values$metadata$File.Name[values$metadata$Group == g]; group_cols <- intersect(colnames(linear_exprs), files_in_group)
      if (length(group_cols) > 1) { group_data <- linear_exprs[, group_cols, drop = FALSE]; cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
      } else { cv_list[[paste0("CV_", g)]] <- NA }
    }
    cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
    df_final <- left_join(df_res, cv_df, by = "Protein.Group") %>% rowwise() %>% mutate(Avg_CV = mean(c_across(starts_with("CV_")), na.rm = TRUE)) %>% ungroup() %>% arrange(Avg_CV) %>% mutate(Stability = ifelse(Avg_CV < 20, "High", "Low")) %>% dplyr::select(Protein.Group, Stability, Avg_CV, logFC, adj.P.Val, starts_with("CV_")) %>% mutate(across(where(is.numeric), ~round(.x, 2)))
    datatable(df_final, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
  })
  
  output$reproducible_code <- renderText({ req(values$repro_log); log_content <- paste(values$repro_log, collapse = "\n"); session_info_text <- paste(capture.output(sessionInfo()), collapse = "\n"); paste(log_content, "\n\n# --- Session Info ---", session_info_text, sep = "\n") })
  output$methodology_text <- renderText({ req(values$fit); paste0("Data was processed using the limpa R package (v1.0+)...") })
  
  observeEvent(input$run_gsea, {
    req(values$fit, input$contrast_selector)
    output$gsea_status <- renderText("Running GSEA... This may take a few minutes.")
    withProgress(message = "Running GSEA", {
      tryCatch({
        incProgress(0.1, detail = "Preparing gene list...")
        de_results <- topTable(values$fit, coef = input$contrast_selector, number = Inf)
        org_db_name <- detect_organism_db(rownames(de_results))
        if (!requireNamespace(org_db_name, quietly = TRUE)) { BiocManager::install(org_db_name, ask = FALSE) }
        library(org_db_name, character.only = TRUE)
        output$gsea_status <- renderText(paste("Detected Organism DB:", org_db_name))
        incProgress(0.3, detail = "Converting Gene IDs...")
        protein_ids <- str_extract(rownames(de_results), "[A-Z0-9]+")
        id_map <- bitr(protein_ids, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = get(org_db_name))
        de_results_mapped <- de_results[rownames(de_results) %in% id_map$UNIPROT, ]; de_results_mapped$ENTREZID <- id_map$ENTREZID[match(str_extract(rownames(de_results_mapped), "[A-Z0-9]+"), id_map$UNIPROT)]
        gene_list <- de_results_mapped$logFC; names(gene_list) <- de_results_mapped$ENTREZID
        gene_list <- sort(gene_list, decreasing = TRUE); gene_list <- gene_list[!duplicated(names(gene_list))]
        incProgress(0.6, detail = "Running GO enrichment...")
        gsea_res <- gseGO(geneList=gene_list, OrgDb=get(org_db_name), keyType="ENTREZID", ont="BP", minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, verbose=FALSE)
        values$gsea_results <- gsea_res; incProgress(1, detail = "Complete.")
        output$gsea_status <- renderText("GSEA Complete.")
      }, error = function(e) { output$gsea_status <- renderText(paste("GSEA Error:", e$message)) })
    })
  })
  
  output$gsea_dot_plot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) dotplot(values$gsea_results, showCategory = 20) + ggtitle("GSEA GO Biological Process") else plot(NULL, xlim=c(0,1), ylim=c(0,1), main="No significant enrichment found.", xaxt='n', yaxt='n') }, height = 500) # FIXED HEIGHT
  output$gsea_emapplot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) emapplot(pairwise_termsim(values$gsea_results), showCategory=20) }, height = 500) # FIXED HEIGHT
  output$gsea_ridgeplot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) ridgeplot(values$gsea_results) }, height = 600) # FIXED HEIGHT
  output$gsea_results_table <- renderDT({ req(values$gsea_results); datatable(as.data.frame(values$gsea_results), options=list(pageLength=10, scrollX=TRUE)) })
  
  observeEvent(input$check_models, { if (nchar(input$user_api_key) < 10) { showNotification("Please enter a valid API Key first.", type="error"); return() }; withProgress(message = "Checking Google Models...", { models <- list_google_models(input$user_api_key); if (length(models) > 0 && !grepl("Error", models[1])) { showModal(modalDialog(title = "Available Models for Your Key", p("Copy one of these into the Model Name box:"), tags$textarea(paste(models, collapse="\n"), rows=10, style="width:100%;"), easyClose = TRUE)) } else { showNotification(paste("Failed to list models:", models), type="error") } }) })
  output$chat_selection_indicator <- renderText({ if (!is.null(values$plot_selected_proteins)) { paste("✅ Current Selection:", length(values$plot_selected_proteins), "Proteins from Plots.") } else { "ℹ️ No proteins selected in plots." } })
  
  observeEvent(input$summarize_data, {
    req(input$user_api_key)
    auto_prompt <- "Analyze this dataset. Identify key quality control issues (if any) by looking at the Group QC stats. Then, summarize the main biological findings from the expression data, focusing on the most significantly differentially expressed proteins."
    values$chat_history <- append(values$chat_history, list(list(role = "user", content = "(Auto-Query: Summarize & Analyze)")))
    withProgress(message = "Auto-Analyzing Dataset...", {
      if (!is.null(values$fit) && !is.null(values$y_protein)) {
        n_max <- 800; df_de <- topTable(values$fit, coef=input$contrast_selector, number=n_max); df_exprs <- as.data.frame(values$y_protein$E[rownames(df_de), ]); df_full <- cbind(Protein = rownames(df_de), df_de, df_exprs)
        incProgress(0.3, detail = "Sending data file..."); current_file_uri <- upload_csv_to_gemini(df_full, input$user_api_key)
        qc_final <- NULL; if(!is.null(values$qc_stats) && !is.null(values$metadata)) { qc_final <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>% dplyr::select(Run, Group, Precursors, Proteins, MS1_Signal) }
        incProgress(0.7, detail = "Thinking..."); ai_reply <- ask_gemini_file_chat(auto_prompt, current_file_uri, qc_final, input$user_api_key, input$model_name, values$plot_selected_proteins)
      } else { ai_reply <- "Please load data and run analysis first." }
      values$chat_history <- append(values$chat_history, list(list(role = "ai", content = ai_reply)))
    })
  })
  
  observeEvent(input$send_chat, {
    req(input$chat_input, input$user_api_key)
    values$chat_history <- append(values$chat_history, list(list(role = "user", content = input$chat_input)))
    withProgress(message = "Processing...", {
      if (!is.null(values$fit) && !is.null(values$y_protein)) {
        n_max <- 800; df_de <- topTable(values$fit, coef=input$contrast_selector, number=n_max)
        if (!is.null(values$plot_selected_proteins)) { missing_ids <- setdiff(values$plot_selected_proteins, rownames(df_de)); if (length(missing_ids) > 0) { valid_missing <- intersect(missing_ids, rownames(values$fit$coefficients)); if(length(valid_missing) > 0) { df_extra <- topTable(values$fit, coef=input$contrast_selector, number=Inf)[valid_missing, ]; df_de <- rbind(df_de, df_extra) } } }
        df_exprs <- as.data.frame(values$y_protein$E[rownames(df_de), ]); df_full <- cbind(Protein = rownames(df_de), df_de, df_exprs)
        incProgress(0.3, detail = "Sending data file..."); current_file_uri <- upload_csv_to_gemini(df_full, input$user_api_key)
        qc_final <- NULL; if(!is.null(values$qc_stats) && !is.null(values$metadata)) { qc_final <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>% dplyr::select(Run, Group, Precursors, Proteins, MS1_Signal) }
        incProgress(0.7, detail = "Thinking..."); ai_reply <- ask_gemini_file_chat(input$chat_input, current_file_uri, qc_final, input$user_api_key, input$model_name, values$plot_selected_proteins)
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
  
  output$download_result_csv <- downloadHandler(
    filename = function() { paste0("Limpa_Results_", make.names(input$contrast_selector), ".csv") },
    content = function(file) { req(values$fit, values$y_protein); de_stats <- topTable(values$fit, coef=input$contrast_selector, number=Inf) %>% rownames_to_column("Protein.Group"); exprs_data <- as.data.frame(values$y_protein$E) %>% rownames_to_column("Protein.Group"); full_data <- left_join(de_stats, exprs_data, by="Protein.Group"); write.csv(full_data, file, row.names=FALSE) }
  )
  
  output$volcano_plot_interactive <- renderPlotly({
    df <- volcano_data(); cols <- c("Not Sig" = "grey", "Significant" = "red")
    p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), text = paste("Protein:", Protein.Group), key = Protein.Group, color = Significance)) +
      geom_point(alpha = 0.6) + scale_color_manual(values = cols) +
      geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff), linetype="dashed") + geom_hline(yintercept = -log10(0.05), linetype="dashed") + theme_minimal()
    df_sel <- df %>% filter(Selected == "Yes")
    if (nrow(df_sel) > 0) p <- p + geom_point(data = df_sel, aes(x=logFC, y=-log10(adj.P.Val)), shape=21, size=4, fill=NA, color="blue", stroke=2)
    ggplotly(p, tooltip = "text", source = "volcano_source") %>% layout(dragmode = "select")
  })
  
  observeEvent(event_data("plotly_selected", source = "volcano_source"), { select_data <- event_data("plotly_selected", source = "volcano_source"); if (!is.null(select_data)) values$plot_selected_proteins <- select_data$key })
  observeEvent(event_data("plotly_click", source = "volcano_source"), { click_data <- event_data("plotly_click", source = "volcano_source"); if (!is.null(click_data)) values$plot_selected_proteins <- click_data$key })
  observeEvent(input$clear_plot_selection, { values$plot_selected_proteins <- NULL })
  
  # --- ADDED: Violin Plot Button Logic ---
  observeEvent(input$show_violin, {
    if (is.null(values$plot_selected_proteins) || length(values$plot_selected_proteins) == 0) {
      showNotification("⚠️ Please select a protein in the Volcano Plot or Table first!", type = "warning")
      return()
    }
    target_protein <- values$plot_selected_proteins[1]; values$temp_violin_target <- target_protein
    showModal(modalDialog(title = paste("Expression Profile:", target_protein), size = "xl", plotOutput("violin_plot_de_popup", height = "500px"), footer = modalButton("Close"), easyClose = TRUE))
  })
  
  output$violin_plot_de_popup <- renderPlot({
    req(values$y_protein, values$temp_violin_target, values$metadata)
    prot_id <- values$temp_violin_target
    exprs_mat <- values$y_protein$E[prot_id, , drop=FALSE]
    long_df <- as.data.frame(exprs_mat) %>% rownames_to_column("Protein") %>% pivot_longer(-Protein, names_to = "File.Name", values_to = "LogIntensity")
    long_df <- left_join(long_df, values$metadata, by="File.Name")
    ggplot(long_df, aes(x = Group, y = LogIntensity, fill = Group)) + geom_violin(alpha = 0.5, trim = FALSE) + geom_jitter(width = 0.2, size = 2, alpha = 0.8) + theme_bw() + labs(title = paste("Protein:", prot_id), y = "Log2 Intensity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }, height = 500)
}

shinyApp(ui, server)
