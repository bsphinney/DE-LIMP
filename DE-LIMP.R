# ==============================================================================
#  LIMP-D: Limpa Proteomics Dashboard (v1.2 - Pathways & Methodology)
#  filename: DIA-NN_LIMPA.R
# ==============================================================================

# --- 1. AUTO-INSTALLATION & SETUP ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("BiocManager not found. Installing...")
  install.packages("BiocManager")
}

# Added gprofiler2 for lightweight Pathway Analysis
required_pkgs <- c("limma", "limpa", "ComplexHeatmap", "shiny", "shinyjs", "plotly", "DT", "tidyr", "tibble", "stringr", "curl", "gprofiler2")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Package '", pkg, "' is missing. Auto-installing..."))
    if (pkg == "gprofiler2") {
      install.packages("gprofiler2")
    } else {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
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
library(httr2)
library(rhandsontable)
library(DT)     
library(arrow)  
library(ComplexHeatmap)
library(shinyjs)
library(plotly)
library(stringr)
library(gprofiler2) # New for Pathways

options(shiny.maxRequestSize = 500 * 1024^2)

if (!requireNamespace("limpa", quietly = TRUE)) {
  stop("CRITICAL ERROR: 'limpa' failed to install. Please run: BiocManager::install('limpa') manually.")
}
library(limpa) 

# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

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

list_google_models <- function(api_key) {
  req <- request("https://generativelanguage.googleapis.com/v1beta/models") %>% req_url_query(key = api_key)
  tryCatch({
    resp <- req_perform(req)
    data <- resp_body_json(resp)
    models <- sapply(data$models, function(x) x$name)
    models <- gsub("^models/", "", models)
    return(models)
  }, error = function(e) { return(paste("Error listing models:", e$message)) })
}

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

ask_gemini_file_chat <- function(user_query, file_uri, qc_df, api_key, model_name, selected_ids = NULL) {
  qc_text <- "No QC Data Available"
  if(!is.null(qc_df)) { qc_text <- paste(capture.output(write.csv(qc_df, row.names=FALSE)), collapse="\n") }

  selection_context <- ""
  if (!is.null(selected_ids) && length(selected_ids) > 0) {
    selection_context <- paste0("\n!!! URGENT: USER SELECTION ACTIVE !!!\nFocus on these proteins:\n", paste(selected_ids, collapse=", "), "\n")
  }

  system_instruction <- paste0(
    "You are a PhD-level expert in proteomics. Sources:\n",
    "1. QC METRICS (Text below): Use for technical quality/groups.\n--- QC START ---\n", qc_text, "\n--- QC END ---\n",
    "2. EXPRESSION DATA (File): Top 800 proteins + User Selection.\n",
    "IMPORTANT: If user asks about 'selected proteins', refer to URGENT section.\n",
    "Output interesting protein IDs like this: [[SELECT: P12345; P67890]]"
  )
  
  base_url <- "https://generativelanguage.googleapis.com/v1beta/models/"
  clean_model <- gsub("^models/", "", model_name)
  full_url <- paste0(base_url, clean_model, ":generateContent")
  
  body <- list(contents = list(list(parts = list(
        list(text = paste0(system_instruction, selection_context, "\n\nUser Question: ", user_query)),
        list(file_data = list(file_uri = file_uri, mime_type = "text/csv"))
  ))))
  
  req <- request(full_url) %>% req_url_query(key = api_key) %>% req_headers("Content-Type" = "application/json") %>% req_body_json(body)
  
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

# ==============================================================================
#  USER INTERFACE (UI)
# ==============================================================================

ui <- page_sidebar(
  title = "LIMP-D Proteomics",
  theme = bs_theme(bootswatch = "flatly"),
  useShinyjs(), 
  
  tags$head(tags$style(HTML("
    .chat-container { height: 500px; overflow-y: auto; border: 1px solid #ddd; padding: 15px; background-color: #f8f9fa; border-radius: 5px; margin-bottom: 15px; }
    .user-msg { text-align: right; margin: 10px 0; }
    .user-msg span { background-color: #007bff; color: white; padding: 8px 12px; border-radius: 15px 15px 0 15px; display: inline-block; max-width: 80%; }
    .ai-msg { text-align: left; margin: 10px 0; }
    .ai-msg span { background-color: #e9ecef; color: #333; padding: 8px 12px; border-radius: 15px 15px 15px 0; display: inline-block; max-width: 80%; }
  "))),
  
  sidebar = sidebar(
    width = 320,
    title = "Controls",
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
    radioButtons("pval_mode", "Volcano Y-Axis:", choices = c("Adjusted P-Value" = "adj", "Raw P-Value" = "raw"), inline = TRUE),
    hr(),
    h5("5. AI Chat"),
    passwordInput("user_api_key", "Gemini API Key", value = "", placeholder = "AIzaSy..."),
    actionButton("check_models", "Check Models", class="btn-warning btn-xs w-100"),
    br(), br(),
    textInput("model_name", "Model Name", value = "gemini-1.5-flash", placeholder = "gemini-1.5-flash")
  ),
  
  navset_card_tab(
    id = "main_tabs", 
    nav_panel("Data Overview", icon = icon("database"), card(card_header("Summary"), card_body(verbatimTextOutput("file_summary")))),
    nav_panel("QC Trends", icon = icon("chart-bar"), layout_columns(col_widths=c(12,12), card(card_header("Trend Analysis"), card_body(layout_columns(col_widths=c(4, 8), selectInput("qc_metric_select", "Metric:", choices = c("Precursors", "Proteins", "MS1_Signal")), radioButtons("qc_sort_order", "Order By:", choices = c("Run Order", "Group"), inline = TRUE)), plotlyOutput("qc_trend_plot", height = "500px"))), card(card_header("Stats Table"), DTOutput("r_qc_table")))),
    nav_panel("QC Plots", icon = icon("chart-scatter"), layout_columns(col_widths=c(6,6), card(card_header("DPC Fit"), plotOutput("dpc_plot", height="400px")), card(card_header("MDS Plot"), plotOutput("mds_plot", height="400px"))), layout_columns(col_widths=c(12), card(card_header("Group QC Distribution (Hover for Info)"), card_body(selectInput("qc_violin_metric", "Metric:", choices = c("Precursors", "Proteins", "MS1_Signal"), width = "200px"), plotlyOutput("qc_group_violin", height = "400px"))))),
    
    # --- DE DASHBOARD ---
    nav_panel("DE Dashboard", icon = icon("table-columns"), 
              # Added Current Contrast Header
              card(card_body(h4(textOutput("current_contrast_label"), style="text-align:center; color: #2c3e50; font-weight:bold;"))),
              layout_columns(col_widths = c(6, 6), 
                card(card_header(div(style="display: flex; justify-content: space-between; align-items: center;", span("Results Table"), div(actionButton("clear_plot_selection", "Reset", class="btn-warning btn-xs"), actionButton("show_violin", "📊 Violin Plot", class="btn-primary btn-xs"), downloadButton("download_results", "💾 Export Results", class="btn-success btn-xs")))), DTOutput("de_table")), 
                card(card_header("Volcano Plot (Click/Box Select)"), plotlyOutput("volcano_plot_interactive", height = "600px"))
              ), 
              card(card_header("Heatmap"), plotOutput("heatmap_plot", height="400px"))
    ),
    
    # --- NEW: PATHWAY ANALYSIS TAB ---
    nav_panel("Pathway Analysis", icon = icon("network-wired"),
              card(
                card_header("Enrichment Analysis (KEGG, Reactome, GO)"),
                card_body(
                  layout_columns(col_widths = c(3, 9),
                    div(
                      h5("Settings"),
                      textInput("gsea_org", "Organism (gProfiler code)", value = "hsapiens", placeholder = "hsapiens, mmusculus..."),
                      p("Common: hsapiens (Human), mmusculus (Mouse), rnorvegicus (Rat)"),
                      hr(),
                      p("This analysis uses the significant proteins from the currently selected contrast in the DE Dashboard."),
                      actionButton("run_gsea", "Run Enrichment", class = "btn-primary w-100")
                    ),
                    div(
                      h5("Enrichment Results"),
                      DTOutput("gsea_table")
                    )
                  )
                )
              )
    ),
    
    nav_panel("Consistent DE", icon = icon("check-double"), card(card_header("High-Consistency Significant Proteins (Ranked by %CV)"), card_body(p("Ranking by %CV (Coefficient of Variation) to find stable markers."), DTOutput("consistent_table")))),
    
    # --- UPDATED: METHODOLOGY TAB ---
    nav_panel("Methodology", icon = icon("file-lines"),
              card(
                card_header("Analysis Methodology (Copy for Publication)"),
                card_body(
                  verbatimTextOutput("methodology_text")
                )
              )
    ),
    
    nav_panel("Reproducibility", icon = icon("code"), card(card_header(div(style="display:flex; justify-content:space-between;", span("Session Code Log"), actionButton("clear_log", "Clear Log", class="btn-danger btn-sm"))), card_body(p("This log records your analysis steps in real-time."), verbatimTextOutput("reproducible_code")))),
    nav_panel("Data Chat", icon = icon("comments"), card(card_header(div(style="display: flex; justify-content: space-between; align-items: center;", span("Chat with Full Data (QC + Expression)"), downloadButton("download_chat", "💾 Save Chat", class="btn-secondary btn-sm"))), card_body(verbatimTextOutput("chat_selection_indicator"), uiOutput("chat_window"), tags$div(style="margin-top: 15px; display: flex; gap: 10px;", textAreaInput("chat_input", label=NULL, placeholder="Ask e.g. 'Which group has higher precursor counts?'", width="100%", rows=2), actionButton("send_chat", "Send", icon=icon("paper-plane"), class="btn-primary", style="height: 54px; margin-top: 2px;")), p("Note: QC Stats (with Groups) + Top 800 Expression Data are sent to AI.", style="font-size: 0.8em; color: green; font-weight: bold; margin-top: 5px;"))))
  )
)

# ==============================================================================
#  SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {
  
  values <- reactiveValues(
    raw_data=NULL, metadata=NULL, fit=NULL, y_protein=NULL, dpc_fit=NULL, 
    status="Waiting...", design=NULL, qc_stats=NULL,
    plot_selected_proteins = NULL, chat_history = list(),
    current_file_uri = NULL,
    gsea_results = NULL,
    code_log = c("# --- LIMP-D Log ---", "library(limpa)", "library(limma)")
  )
  
  add_to_log <- function(comment, code) { values$code_log <- c(values$code_log, c(paste0("\n# ", format(Sys.time(), "[%H:%M:%S] "), comment), code)) }
  observeEvent(input$clear_log, { values$code_log <- c("# --- Log Cleared ---", "library(limpa)", "library(limma)") })

  observeEvent(input$report_file, {
    req(input$report_file)
    withProgress(message = "Loading...", {
      incProgress(0.2, detail = "QC..."); values$qc_stats <- get_diann_stats_r(input$report_file$datapath)
      incProgress(0.5, detail = "Reading...")
      tryCatch({
        values$raw_data <- limpa::readDIANN(input$report_file$datapath, format="parquet", q.cutoffs=input$q_cutoff)
        fnames <- sort(colnames(values$raw_data$E))
        values$metadata <- data.frame(ID = 1:length(fnames), File.Name = fnames, Group = rep("", length(fnames)), stringsAsFactors=FALSE)
        click("open_setup")
        add_to_log("USER LOADED DATA", paste0("dat <- readDIANN('", input$report_file$name, "', format='parquet', q.cutoffs=", input$q_cutoff, ")"))
      }, error=function(e) { showNotification(paste("Error:", e$message), type="error") })
    })
  })
  
  output$file_summary <- renderText({ if(is.null(values$metadata)) return("No data loaded."); paste("Total Files:", nrow(values$metadata), "\nAssigned Groups:", length(unique(values$metadata$Group[values$metadata$Group!=""]))) })
  output$qc_trend_plot <- renderPlotly({ req(values$qc_stats, input$qc_metric_select, values$metadata); df <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")); if (input$qc_sort_order == "Group") { df <- df %>% arrange(Group, Run) } else { df <- df %>% arrange(Run) }; df$Sort_Index <- factor(1:nrow(df), levels = 1:nrow(df)); metric <- input$qc_metric_select; df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Group:</b> ", df$Group, "<br><b>", metric, ":</b> ", round(df[[metric]], 2)); p <- ggplot(df, aes(x = Sort_Index, y = .data[[metric]], fill = Group, text = Tooltip)) + geom_bar(stat = "identity", width = 0.8) + theme_minimal() + labs(title = paste(metric, "per Run"), x = "Sample Index (Sorted)", y = metric) + theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size=8)); ggplotly(p, tooltip = "text") %>% config(displayModeBar = FALSE) })
  output$r_qc_table <- renderDT({ req(values$qc_stats); df_display <- values$qc_stats %>% arrange(Run) %>% mutate(ID = 1:n()) %>% select(ID, Run, everything()); datatable(df_display, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) })
  output$qc_group_violin <- renderPlotly({ req(values$qc_stats, values$metadata, input$qc_violin_metric); df <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")); metric <- input$qc_violin_metric; df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Val:</b> ", round(df[[metric]], 2)); p <- ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) + geom_violin(alpha = 0.5, trim = FALSE) + geom_jitter(aes(text = Tooltip), width = 0.2, size = 2, alpha = 0.8, color = "black") + theme_bw() + labs(title = paste("Distribution of", metric), x = "Group", y = metric) + theme(legend.position = "none"); ggplotly(p, tooltip = "text") })

  observeEvent(input$open_setup, { req(values$metadata); showModal(modalDialog(title = "Assign Groups", size = "xl", actionButton("guess_groups", "🪄 Auto-Guess", class="btn-info btn-sm"), br(), br(), rHandsontableOutput("hot_metadata_modal"), footer = tagList(modalButton("Cancel"), actionButton("save_groups", "Save & Close", class="btn-success")))) })
  output$hot_metadata_modal <- renderRHandsontable({ req(values$metadata); rhandsontable(values$metadata, rowHeaders=NULL, stretchH="all", height=500, width="100%") %>% hot_col("ID", readOnly=TRUE, width=50) %>% hot_col("File.Name", readOnly=TRUE) %>% hot_col("Group", type="text") })
  observeEvent(input$guess_groups, { req(values$metadata); meta <- if(!is.null(input$hot_metadata_modal)) hot_to_r(input$hot_metadata_modal) else values$metadata; meta$Group <- ifelse(grepl("affinisep", meta$File.Name, ignore.case=T), "Affinisep", ifelse(grepl("evosep", meta$File.Name, ignore.case=T), "Evosep", ifelse(grepl("control", meta$File.Name, ignore.case=T), "Control", ifelse(grepl("treat", meta$File.Name, ignore.case=T), "Treatment", "")))); values$metadata <- meta; output$hot_metadata_modal <- renderRHandsontable({ rhandsontable(values$metadata, rowHeaders=NULL, stretchH="all", height=500, width="100%") %>% hot_col("ID", readOnly=TRUE, width=50) %>% hot_col("File.Name", readOnly=TRUE) %>% hot_col("Group", type="text") }) })
  observeEvent(input$save_groups, { if(!is.null(input$hot_metadata_modal)) values$metadata <- hot_to_r(input$hot_metadata_modal); removeModal(); showNotification("Groups saved!", type="message"); add_to_log("GROUPS SAVED", "metadata <- ... (user defined groups)") })
  
  observeEvent(input$run, {
    req(values$raw_data, values$metadata); meta <- values$metadata; meta$Group <- trimws(meta$Group)
    if(length(unique(meta$Group)) < 2) { showNotification("Error: Need 2+ groups.", type="error"); return() }
    withProgress(message='Running Pipeline...', {
      tryCatch({
        dat <- values$raw_data; dpcfit <- limpa::dpc(dat); values$dpc_fit <- dpcfit; values$y_protein <- limpa::dpcQuant(dat, "Protein.Group", dpc=dpcfit)
        rownames(meta) <- meta$File.Name; meta <- meta[colnames(dat$E), ]; meta$Group <- make.names(meta$Group); groups <- factor(meta$Group)
        design <- model.matrix(~ 0 + groups); colnames(design) <- levels(groups)
        combs <- combn(levels(groups), 2); forms <- apply(combs, 2, function(x) paste(x[2], "-", x[1]))
        fit <- limpa::dpcDE(values$y_protein, design, plot=FALSE); fit <- contrasts.fit(fit, makeContrasts(contrasts=forms, levels=design)); fit <- eBayes(fit)
        values$fit <- fit; updateSelectInput(session, "contrast_selector", choices=forms); values$status <- "✅ Complete!"; nav_select("main_tabs", "QC Plots")
        add_to_log("PIPELINE EXECUTION", paste("Contrasts:", paste(forms, collapse=", ")))
      }, error=function(e) showNotification(paste("Error:", e$message), type="error"))
    })
  })
  output$run_status_msg <- renderText({ values$status })

  output$dpc_plot <- renderPlot({ req(values$dpc_fit); limpa::plotDPC(values$dpc_fit) })
  output$mds_plot <- renderPlot({ req(values$y_protein, values$metadata); meta <- values$metadata[match(colnames(values$y_protein$E), values$metadata$File.Name), ]; grps <- factor(meta$Group); cols <- rainbow(length(levels(grps))); limpa::plotMDSUsingSEs(values$y_protein, pch=16, main="MDS Plot", col=cols[grps]); legend("topright", legend=levels(grps), col=cols, pch=16) })
  observeEvent(input$contrast_selector, { req(values$fit); add_to_log("CONTRAST CHANGE", paste0("coef='", input$contrast_selector, "'")) })

  # --- CURRENT CONTRAST LABEL ---
  output$current_contrast_label <- renderText({ 
    req(input$contrast_selector)
    paste("Currently Viewing Contrast:", input$contrast_selector)
  })

  # --- UPDATED VOLCANO with SUBTITLE ---
  volcano_data <- reactive({
    req(values$fit, input$contrast_selector)
    df <- topTable(values$fit, coef=input$contrast_selector, number=Inf)
    
    if (input$pval_mode == "raw") { df$Plot_PVal <- df$P.Value; y_label <- "-Log10(Raw P-Value)" } else { df$Plot_PVal <- df$adj.P.Val; y_label <- "-Log10(Adj. P-Value)" }
    
    df$Significance <- "Not Sig"
    df$Significance[df$adj.P.Val < 0.05 & abs(df$logFC) > input$logfc_cutoff] <- "Significant"
    df$Selected <- "No"; if (!is.null(values$plot_selected_proteins)) { df$Selected[df$Protein.Group %in% values$plot_selected_proteins] <- "Yes" }
    
    df$Tooltip <- paste0("<b>Protein:</b> ", df$Protein.Group, "<br><b>LogFC:</b> ", round(df$logFC, 3), "<br><b>Raw P:</b> ", formatC(df$P.Value, format="e", digits=2), "<br><b>Adj P:</b> ", formatC(df$adj.P.Val, format="e", digits=2))
    list(df = df, y_label = y_label)
  })
  
  output$volcano_plot_interactive <- renderPlotly({
    data_list <- volcano_data(); df <- data_list$df; y_lbl <- data_list$y_label
    cols <- c("Not Sig" = "grey", "Significant" = "red")
    
    # NEW: Subtitle Logic
    sig_label <- paste0("Significance Defined By: Adj.P < 0.05 & |LogFC| > ", input$logfc_cutoff)
    
    p <- ggplot(df, aes(x = logFC, y = -log10(Plot_PVal), text = Tooltip, key = Protein.Group, color = Significance)) +
      geom_point(alpha = 0.6) + scale_color_manual(values = cols) + 
      geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff), linetype="dashed") + 
      geom_hline(yintercept = -log10(0.05), linetype="dashed") + 
      theme_minimal() + 
      labs(y = y_lbl, subtitle = sig_label) +
      theme(plot.subtitle = element_text(size = 9, color = "blue", face = "italic"))
    
    df_sel <- df %>% filter(Selected == "Yes"); if (nrow(df_sel) > 0) p <- p + geom_point(data = df_sel, aes(x=logFC, y=-log10(Plot_PVal)), shape=21, size=4, fill=NA, color="blue", stroke=2)
    ggplotly(p, tooltip = "text", source = "volcano_source") %>% layout(dragmode = "select")
  })
  
  observeEvent(event_data("plotly_selected", source = "volcano_source"), { select_data <- event_data("plotly_selected", source = "volcano_source"); if (!is.null(select_data)) values$plot_selected_proteins <- select_data$key })
  observeEvent(event_data("plotly_click", source = "volcano_source"), { click_data <- event_data("plotly_click", source = "volcano_source"); if (!is.null(click_data)) values$plot_selected_proteins <- click_data$key })
  observeEvent(input$clear_plot_selection, { values$plot_selected_proteins <- NULL })
  
  output$de_table <- renderDT({
    df <- volcano_data()$df %>% mutate(across(c(logFC, P.Value, adj.P.Val), is.numeric)) %>% select(Protein.Group, logFC, P.Value, adj.P.Val, Significance)
    if (!is.null(values$plot_selected_proteins)) df <- df %>% filter(Protein.Group %in% values$plot_selected_proteins)
    datatable(df, selection = "multiple", options = list(pageLength = 10, scrollX = TRUE)) %>% formatRound(columns = "logFC", digits = 3) %>% formatSignif(columns = c("P.Value", "adj.P.Val"), digits = 3)
  })
  
  output$heatmap_plot <- renderPlot({ req(values$fit, values$y_protein, input$contrast_selector); df_volc <- volcano_data()$df; prot_ids <- NULL; if (!is.null(input$de_table_rows_selected)) { current_table_data <- df_volc; if (!is.null(values$plot_selected_proteins)) current_table_data <- current_table_data %>% filter(Protein.Group %in% values$plot_selected_proteins); prot_ids <- current_table_data$Protein.Group[input$de_table_rows_selected] } else if (!is.null(values$plot_selected_proteins)) { prot_ids <- values$plot_selected_proteins; if(length(prot_ids) > 50) prot_ids <- head(prot_ids, 50) } else { top_prots <- topTable(values$fit, coef=input$contrast_selector, number=20); prot_ids <- rownames(top_prots) }; valid_ids <- intersect(prot_ids, rownames(values$y_protein$E)); if (length(valid_ids) == 0) return(NULL); mat <- values$y_protein$E[valid_ids, , drop=FALSE]; mat_z <- t(apply(mat, 1, cal_z_score)); meta <- values$metadata[match(colnames(mat), values$metadata$File.Name), ]; groups <- factor(meta$Group); ha <- HeatmapAnnotation(Group = groups, col = list(Group = setNames(rainbow(length(levels(groups))), levels(groups)))); Heatmap(mat_z, name="Z-score", top_annotation = ha, cluster_rows=TRUE, cluster_columns=TRUE, show_column_names=FALSE) })
  
  output$consistent_table <- renderDT({ req(values$fit, values$y_protein, input$contrast_selector, values$metadata); df_res <- topTable(values$fit, coef=input$contrast_selector, number=Inf) %>% filter(adj.P.Val < 0.05) %>% rownames_to_column("Protein_ID_Key"); if(nrow(df_res) == 0) return(NULL); raw_exprs <- values$y_protein$E[df_res$Protein_ID_Key, , drop=FALSE]; meta <- values$metadata; linear_exprs <- 2^raw_exprs; unique_groups <- unique(meta$Group); cv_list <- list(); for(g in unique_groups) { files_in_group <- meta$File.Name[meta$Group == g]; group_data <- linear_exprs[, intersect(colnames(linear_exprs), files_in_group), drop=FALSE]; group_cv <- apply(group_data, 1, function(x) { if(all(is.na(x)) || mean(x, na.rm=T) == 0) return(NA); (sd(x, na.rm=T) / mean(x, na.rm=T)) * 100 }); cv_list[[paste0("CV_", g)]] <- group_cv }; cv_df <- as.data.frame(cv_list); cv_df$Protein_ID_Key <- rownames(cv_df); cv_df$Avg_CV <- rowMeans(cv_df[, -ncol(cv_df)], na.rm=TRUE); df_final <- left_join(df_res, cv_df, by="Protein_ID_Key") %>% arrange(Avg_CV) %>% select(Protein_ID_Key, Avg_CV, starts_with("CV_"), logFC, adj.P.Val) %>% mutate(across(where(is.numeric), \(x) round(x, 2))) %>% mutate(Stability = ifelse(Avg_CV < 20, "High", "Low")) %>% rename(Protein.Group = Protein_ID_Key) %>% select(Protein.Group, Stability, Avg_CV, everything()); datatable(df_final, options = list(pageLength = 15, scrollX = TRUE)) })
  
  # --- NEW: PATHWAY ANALYSIS LOGIC ---
  observeEvent(input$run_gsea, {
    req(values$fit, input$contrast_selector)
    withProgress(message = "Running Pathway Enrichment (ORA)...", {
      
      # 1. Get Significant Proteins
      df_all <- topTable(values$fit, coef=input$contrast_selector, number=Inf)
      sig_prots <- rownames(df_all)[df_all$adj.P.Val < 0.05 & abs(df_all$logFC) > input$logfc_cutoff]
      
      if (length(sig_prots) == 0) {
        showNotification("No significant proteins found for enrichment.", type="warning")
        values$gsea_results <- NULL
        return()
      }
      
      incProgress(0.2, detail = paste("Querying gProfiler with", length(sig_prots), "proteins..."))
      
      # 2. Run gProfiler2 (KEGG, REAC, GO)
      tryCatch({
        gost_res <- gprofiler2::gost(
          query = sig_proteins = sig_prots, 
          organism = input$gsea_org, 
          ordered_query = FALSE, 
          multi_query = FALSE, 
          significant = TRUE, 
          exclude_iea = TRUE, 
          measure_underrepresentation = FALSE, 
          evcodes = FALSE, 
          user_threshold = 0.05, 
          correction_method = "g_SCS", 
          domain_scope = "annotated", 
          custom_bg = NULL, 
          numeric_ns = "", 
          sources = c("GO:BP", "KEGG", "REAC")
        )
        
        if (is.null(gost_res) || is.null(gost_res$result)) {
          showNotification("No enriched pathways found.", type="warning")
          values$gsea_results <- NULL
        } else {
          values$gsea_results <- gost_res$result %>%
            select(source, term_name, term_size, intersection_size, p_value) %>%
            arrange(p_value)
        }
        
      }, error = function(e) {
        showNotification(paste("GSEA Error:", e$message), type="error")
      })
    })
  })
  
  output$gsea_table <- renderDT({
    req(values$gsea_results)
    datatable(values$gsea_results, options = list(pageLength = 10, scrollX = TRUE)) %>%
      formatSignif(columns = "p_value", digits = 3)
  })

  # --- NEW: DYNAMIC METHODOLOGY TEXT ---
  output$methodology_text <- renderText({
    req(values$raw_data)
    
    # Gather Versions
    r_ver <- R.version.string
    limma_ver <- packageVersion("limma")
    limpa_ver <- tryCatch(packageVersion("limpa"), error = function(e) "Unknown")
    
    # Build Text
    paste0(
      "Data analysis was performed using the LIMP-D pipeline running in ", r_ver, ".\n\n",
      "1. Data Processing:\n",
      "   Raw DIA-NN output was processed using the 'limpa' package (v", limpa_ver, "). ",
      "Data Processing & Cleaning (DPC) normalization was applied to correct for systematic bias.\n\n",
      "2. Statistical Analysis:\n",
      "   Differential expression analysis was conducted using 'limma' (v", limma_ver, "). ",
      "Linear models were fitted for each protein, and empirical Bayes moderation was used to calculate statistics.\n\n",
      "3. Significance Thresholds:\n",
      "   Proteins were considered significantly differentially expressed if they met the following criteria:\n",
      "   - Adjusted P-value (Benjamini-Hochberg) < 0.05\n",
      "   - Absolute Log2 Fold Change > ", input$logfc_cutoff, "\n\n",
      "4. Quality Control:\n",
      "   Q-value cutoff for precursor filtering was set to ", input$q_cutoff, "."
    )
  })

  output$download_results <- downloadHandler(filename = function() { paste0("Limpa_Results_", make.names(input$contrast_selector), ".csv") }, content = function(file) { req(values$fit, values$y_protein); de_stats <- topTable(values$fit, coef=input$contrast_selector, number=Inf); de_stats$Row_ID_Temp <- rownames(de_stats); exprs_data <- as.data.frame(values$y_protein$E); exprs_data$Row_ID_Temp <- rownames(exprs_data); full_data <- left_join(de_stats, exprs_data, by="Row_ID_Temp"); if(!"Protein.Group" %in% names(full_data)) { full_data <- full_data %>% rename(Protein.Group = Row_ID_Temp) } else { full_data <- full_data %>% select(-Row_ID_Temp) }; write.csv(full_data, file, row.names=FALSE); add_to_log("USER EXPORTED CSV", paste0("write.csv(results, 'Limpa_Results_", make.names(input$contrast_selector), ".csv')")) })
  output$reproducible_code <- renderText({ paste(values$code_log, collapse="\n") })
  observeEvent(input$check_models, { if (nchar(input$user_api_key) < 10) { showNotification("Please enter a valid API Key first.", type="error"); return() }; withProgress(message = "Checking Google Models...", { models <- list_google_models(input$user_api_key); if (length(models) > 0 && !grepl("Error", models[1])) { showModal(modalDialog(title = "Available Models for Your Key", p("Copy one of these into the Model Name box:"), tags$textarea(paste(models, collapse="\n"), rows=10, style="width:100%;"), easyClose = TRUE)) } else { showNotification(paste("Failed to list models:", models), type="error") } }) })
  output$chat_selection_indicator <- renderText({ if (!is.null(values$plot_selected_proteins)) { n_sel <- length(values$plot_selected_proteins); paste("✅ Current Selection:", n_sel, "Proteins from Plots.") } else { "ℹ️ No proteins selected in plots." } })
  observeEvent(input$send_chat, { req(input$chat_input, input$user_api_key); values$chat_history <- append(values$chat_history, list(list(role = "user", content = input$chat_input))); withProgress(message = "Processing (Top 800 + Selection + QC)...", { if (!is.null(values$fit) && !is.null(values$y_protein)) { n_max <- 800; df_de <- topTable(values$fit, coef=input$contrast_selector, number=n_max); if (!is.null(values$plot_selected_proteins)) { missing_ids <- setdiff(values$plot_selected_proteins, rownames(df_de)); if (length(missing_ids) > 0) { valid_missing <- intersect(missing_ids, rownames(values$fit$coefficients)); if(length(valid_missing) > 0) { df_extra <- topTable(values$fit, coef=input$contrast_selector, number=Inf)[valid_missing, ]; df_de <- rbind(df_de, df_extra) } } }; df_exprs <- as.data.frame(values$y_protein$E[rownames(df_de), ]); df_full <- cbind(Protein = rownames(df_de), df_de, df_exprs); incProgress(0.3, detail = "Sending data file..."); current_file_uri <- upload_csv_to_gemini(df_full, input$user_api_key); qc_final <- NULL; if(!is.null(values$qc_stats) && !is.null(values$metadata)) { qc_final <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>% select(Run, Group, Precursors, Proteins, MS1_Signal) }; incProgress(0.7, detail = "Thinking..."); ai_reply <- ask_gemini_file_chat(input$chat_input, current_file_uri, qc_final, input$user_api_key, input$model_name, values$plot_selected_proteins) } else { ai_reply <- "Please load data and run analysis first." } }); ai_selected <- str_extract(ai_reply, "\\[\\[SELECT:.*?\\]\\]"); if (!is.na(ai_selected)) { raw_ids <- gsub("\\[\\[SELECT:|\\]\\]", "", ai_selected); id_vec <- unlist(strsplit(raw_ids, "[,;]\\s*")); values$plot_selected_proteins <- trimws(id_vec); ai_reply <- gsub("\\[\\[SELECT:.*?\\]\\]", "", ai_reply); ai_reply <- paste0(ai_reply, "\n\n*(I have updated your plots with these highlighted proteins.)*") }; values$chat_history <- append(values$chat_history, list(list(role = "ai", content = ai_reply))); updateTextAreaInput(session, "chat_input", value = "") })
  output$chat_window <- renderUI({ chat_content <- lapply(values$chat_history, function(msg) { if (msg$role == "user") { div(class = "user-msg", span(msg$content)) } else { div(class = "ai-msg", span(markdown(msg$content))) } }); div(class = "chat-container", chat_content) })
  output$download_chat <- downloadHandler(filename = function() { paste0("Limpa_Chat_History_", Sys.Date(), ".txt") }, content = function(file) { text_out <- sapply(values$chat_history, function(msg) { role <- if(msg$role == "user") "YOU: " else "GEMINI: "; paste0(role, msg$content, "\n--------------------------------------------------\n") }); writeLines(unlist(text_out), file) })
}

shinyApp(ui, server)
