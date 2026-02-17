# ==============================================================================
#  DE-LIMP: Differential Expression & Limpa Proteomics App
#  (Formerly LIMP-D)
#  Status: Production Ready (Hugging Face Compatible v1.2)
# ==============================================================================

# Set CRAN mirror to avoid interactive popup (especially in VS Code)
options(repos = c(CRAN = "https://cloud.r-project.org"))

# --- 1. AUTO-INSTALLATION & SETUP ---
# IMPORTANT: Install packages BEFORE loading libraries to avoid conflicts

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("BiocManager not found. Installing...")
  install.packages("BiocManager", quiet = TRUE)
}

# Check for limpa and install if needed
if (!requireNamespace("limpa", quietly = TRUE)) {
  message("Package 'limpa' is missing. Attempting installation...")

  r_version <- getRversion()
  bioc_version <- as.character(BiocManager::version())
  message(paste0("R version: ", r_version, ", Bioconductor version: ", bioc_version))

  # Try installing from Bioconductor (may require devel version)
  limpa_installed <- tryCatch({
    suppressWarnings({
      BiocManager::install("limpa", ask = FALSE, update = FALSE, quiet = TRUE)
    })
    requireNamespace("limpa", quietly = TRUE)
  }, error = function(e) FALSE)

  # If standard Bioconductor failed, try development version
  if (!limpa_installed) {
    message("limpa not found in release version. Trying Bioconductor devel...")
    limpa_installed <- tryCatch({
      suppressWarnings({
        BiocManager::install(version = "devel", ask = FALSE, update = FALSE)
        BiocManager::install("limpa", ask = FALSE, update = FALSE, quiet = TRUE)
      })
      requireNamespace("limpa", quietly = TRUE)
    }, error = function(e) FALSE)
  }

  # Final check
  if (!limpa_installed) {

    # Detect platform for specific instructions
    os_type <- Sys.info()["sysname"]
    download_url <- if (os_type == "Darwin") {
      "https://cloud.r-project.org/bin/macosx/"
    } else if (os_type == "Windows") {
      "https://cloud.r-project.org/bin/windows/base/"
    } else {
      "https://cloud.r-project.org/bin/linux/"
    }

    stop(paste0(
      "\n\n╔════════════════════════════════════════════════════════════════╗\n",
      "║          LIMPA INSTALLATION FAILED - R UPGRADE NEEDED          ║\n",
      "╚════════════════════════════════════════════════════════════════╝\n\n",
      "Current setup:\n",
      "  • R version: ", r_version, " (NEED: 4.5+)\n",
      "  • Bioconductor: ", bioc_version, " (NEED: 3.22+)\n",
      "  • Platform: ", os_type, "\n\n",
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n",
      "UPGRADE INSTRUCTIONS:\n",
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
      "1. Download R 4.5+ for your platform:\n",
      "   ", download_url, "\n\n",
      if (os_type == "Darwin") {
        "2. Install the .pkg file (R will be upgraded in-place)\n"
      } else if (os_type == "Windows") {
        "2. Run the installer .exe file\n"
      } else {
        "2. Follow platform-specific installation instructions\n"
      },
      "\n3. Restart VSCode/RStudio completely\n",
      "\n4. Verify upgrade by running: R.version.string\n",
      "\n5. Rerun this script - limpa will install automatically\n\n",
      "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n",
      "More info: https://bioconductor.org/packages/limpa/\n\n"
    ))
  } else {
    message("✓ limpa installed successfully!")
  }
}

# Required packages (excluding limpa which was handled above)
required_pkgs <- c("shiny", "bslib", "readr", "tibble", "dplyr", "tidyr",
                   "ggplot2", "httr2", "rhandsontable", "DT", "arrow",
                   "ComplexHeatmap", "shinyjs", "plotly", "stringr", "limma",
                   "clusterProfiler", "AnnotationDbi", "org.Hs.eg.db", "org.Mm.eg.db",
                   "enrichplot", "ggridges", "ggrepel", "markdown", "curl")

# Only install truly missing packages (don't update already-loaded packages)
missing_pkgs <- character(0)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_pkgs <- c(missing_pkgs, pkg)
  }
}

if (length(missing_pkgs) > 0) {
  message(paste0("Installing missing packages: ", paste(missing_pkgs, collapse = ", ")))
  BiocManager::install(missing_pkgs, ask = FALSE, update = FALSE, quiet = TRUE)
}

# --- 2. SERVER CONFIGURATION ---
options(repos = c(BiocManager::repositories(), CRAN = "https://cloud.r-project.org"))

library(shiny)
library(bslib)

# Verify bslib version supports responsive UI components
if (packageVersion("bslib") < "0.5.0") {
  stop(paste0(
    "bslib >= 0.5.0 required for responsive UI components.\n",
    "Current version: ", packageVersion("bslib"), "\n",
    "Please upgrade: install.packages('bslib')"
  ))
}

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

# Verify Limpa installation
if (!requireNamespace("limpa", quietly = TRUE)) {
  os_type <- Sys.info()["sysname"]
  download_url <- if (os_type == "Darwin") {
    "https://cloud.r-project.org/bin/macosx/"
  } else if (os_type == "Windows") {
    "https://cloud.r-project.org/bin/windows/base/"
  } else {
    "https://cloud.r-project.org/bin/linux/"
  }

  stop(paste0(
    "\n\n╔══════════════════════════════════════════════════════════╗\n",
    "║     CRITICAL: limpa package not found                    ║\n",
    "╚══════════════════════════════════════════════════════════╝\n\n",
    "Your R version: ", getRversion(), " (NEED: 4.5+)\n\n",
    "Upgrade R from: ", download_url, "\n",
    "Then run: BiocManager::install('limpa')\n\n"
  ))
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
        Proteins = if(has_pg_q) {
          n_distinct(Protein.Group[PG.Q.Value <= 0.01])
        } else {
          n_distinct(Protein.Group)
        },
        MS1_Signal = if(has_ms1) sum(Ms1.Apex.Area, na.rm = TRUE) else NA_real_,
        .groups = 'drop'
      ) %>% arrange(Run) 
    
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

    /* === RESPONSIVE UI ADDITIONS === */

    /* DE Dashboard responsive grid */
    .de-dashboard-grid {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 1rem;
      margin-bottom: 1rem;
    }

    @media (max-width: 1200px) {
      .de-dashboard-grid {
        grid-template-columns: 1fr;
      }
    }

    /* Card body with internal scroll for tables */
    .card-body-scroll {
      overflow-y: auto;
      overflow-x: auto;
      max-height: calc(100vh - 380px);
    }

    /* Accordion compact styling */
    .accordion {
      margin-top: 1rem;
    }

    .accordion-button {
      padding: 0.75rem 1rem;
      font-size: 0.95rem;
    }

    /* Viewport-relative plot containers */
    .plot-container-vh {
      min-height: 400px;
      max-height: 85vh;
    }

    /* Compact inline controls */
    .controls-inline {
      display: flex;
      align-items: center;
      gap: 10px;
      flex-wrap: wrap;
    }

    /* Sub-tab navigation styling */
    .nav-tabs .nav-link {
      padding: 0.5rem 1rem;
      font-size: 0.9rem;
    }
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
    actionButton("load_example", "📊 Load Example Data", class = "btn-info btn-sm w-100",
                 style = "margin-bottom: 10px;"),
    numericInput("q_cutoff", "Q-Value Cutoff", value = 0.01, min = 0, max = 0.1, step = 0.01),
    hr(),
    h5("2. Session"),
    div(style="display: flex; gap: 5px; margin-bottom: 5px;",
      downloadButton("save_session", "Save", class = "btn-primary w-50", icon = icon("download")),
      actionButton("load_session_btn", "Load", class = "btn-outline-primary w-50", icon = icon("upload"))
    ),
    hr(),
    h5("3. Explore Results"),
    sliderInput("logfc_cutoff", "Min Log2 Fold Change:", min=0, max=5, value=1, step=0.1),
    hr(),
    h5("4. AI Chat"),
    passwordInput("user_api_key", "Gemini API Key", value = "", placeholder = "AIzaSy..."),
    actionButton("check_models", "Check Models", class="btn-warning btn-xs w-100"),
    br(), br(),
    textInput("model_name", "Model Name", value = "gemini-3-flash-preview", placeholder = "gemini-3-flash-preview"),
    hr(),
    h5("5. XIC Viewer"),
    p(class = "text-muted small",
      "Load .xic.parquet files from DIA-NN to inspect chromatograms."),
    textInput("xic_dir_input", "XIC Directory Path:",
      placeholder = "Auto-detected or paste path here"),
    actionButton("xic_load_dir", "Load XICs", class = "btn-outline-info btn-sm w-100",
      icon = icon("wave-square")),
    uiOutput("xic_status_badge")
  ),
  
  navset_card_tab(
    id = "main_tabs", 
    
    nav_panel("Data Overview", icon = icon("database"),
              # Data views as tabs
              navset_card_tab(
                id = "data_overview_tabs",

                nav_panel("Assign Groups & Run",
                  icon = icon("table"),
                  # Tip banner
                  div(style="background-color: #e7f3ff; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
                    icon("info-circle"),
                    strong(" Tip: "),
                    "Assign experimental groups (required). Covariate columns are optional - customize names and include in model as needed."
                  ),

                  # Top row: Auto-Guess + Covariates + Run Pipeline (responsive)
                  div(style="display: flex; flex-wrap: wrap; gap: 12px; margin-bottom: 12px; align-items: flex-start;",
                    # Auto-Guess + Template buttons
                    div(style="min-width: 160px;",
                      actionButton("guess_groups", "Auto-Guess Groups", class="btn-info btn-sm w-100",
                        icon = icon("wand-magic-sparkles")),
                      div(style="display: flex; gap: 5px; margin-top: 8px;",
                        downloadButton("export_template", "Export", class="btn-outline-secondary btn-sm"),
                        actionButton("import_template", "Import", class="btn-outline-secondary btn-sm")
                      )
                    ),

                    # Covariates (compact inline)
                    div(style="flex: 1; min-width: 250px;",
                      strong("Covariates:", style="font-size: 0.9em;"),
                      div(style="display: flex; flex-wrap: wrap; gap: 8px; margin-top: 5px;",
                        div(style="min-width: 110px;",
                          checkboxInput("include_batch", "Batch", value = FALSE),
                          textInput("batch_label", NULL, value = "Batch", placeholder = "e.g., Batch")
                        ),
                        div(style="min-width: 130px;",
                          checkboxInput("include_cov1", NULL, value = FALSE),
                          textInput("cov1_label", "Name:", value = "Covariate1",
                                   placeholder = "e.g., Sex, Diet")
                        ),
                        div(style="min-width: 130px;",
                          checkboxInput("include_cov2", NULL, value = FALSE),
                          textInput("cov2_label", "Name:", value = "Covariate2",
                                   placeholder = "e.g., Age, Time")
                        )
                      )
                    ),

                    # Run Pipeline button
                    div(style="min-width: 150px; display: flex; align-items: center;",
                      actionButton("run_pipeline", "Run Pipeline",
                        class="btn-success btn-lg w-100", icon = icon("play"),
                        style="padding: 12px; font-size: 1.05em; white-space: nowrap;")
                    )
                  ),

                  # Metadata table (with overflow scroll)
                  div(style="overflow-y: auto; max-height: calc(100vh - 420px);",
                    rHandsontableOutput("hot_metadata")
                  )
                ),

                nav_panel("Signal Distribution",
                  icon = icon("chart-area"),
                  # Comparison selector banner
                  div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 12px 15px; border-radius: 8px; margin-bottom: 15px; display: flex; align-items: center; justify-content: space-between; gap: 10px; flex-wrap: wrap;",
                    div(style = "display: flex; align-items: center; gap: 10px;",
                      icon("microscope"),
                      span("Viewing Comparison:", style = "font-weight: 500;")
                    ),
                    div(style = "flex-grow: 1; max-width: 400px;",
                      selectInput("contrast_selector_signal", NULL,
                        choices = NULL,
                        width = "100%"
                      )
                    )
                  ),
                  # Control buttons
                  div(style = "display: flex; justify-content: flex-end; align-items: center; margin-bottom: 10px;",
                    actionButton("fullscreen_signal", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
                  ),
                  plotOutput("protein_signal_plot", height = "calc(100vh - 370px)")
                ),

                nav_panel("Dataset Summary",
                  icon = icon("info-circle"),
                  uiOutput("dataset_summary_content")
                ),

                nav_panel("Group QC Summary",
                  icon = icon("table"),
                  DTOutput("group_summary_table")
                ),

                nav_panel("Expression Grid",
                  icon = icon("th"),
                  # Comparison selector banner
                  div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 12px 15px; border-radius: 8px; margin-bottom: 15px; display: flex; align-items: center; justify-content: space-between; gap: 10px; flex-wrap: wrap;",
                    div(style = "display: flex; align-items: center; gap: 10px;",
                      icon("microscope"),
                      span("Viewing Comparison:", style = "font-weight: 500;")
                    ),
                    div(style = "flex-grow: 1; max-width: 400px;",
                      selectInput("contrast_selector_grid", NULL,
                        choices = NULL,
                        width = "100%"
                      )
                    )
                  ),
                  # Legend and file mapping
                  div(style = "margin-bottom: 15px;",
                    uiOutput("grid_legend_ui"),
                    uiOutput("grid_file_map_ui")
                  ),
                  # Control buttons
                  div(style = "margin-bottom: 10px;",
                    actionButton("grid_reset_selection", "Show All / Clear Selection", class = "btn-warning btn-sm"),
                    downloadButton("download_grid_data", "💾 Export Full Table", class = "btn-success btn-sm")
                  ),
                  # Grid table
                  DTOutput("grid_view_table")
                ),

                nav_panel("AI Summary",
                  icon = icon("robot"),
                  div(style = "padding: 20px;",
                    # Header section
                    div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
                      tags$h4(icon("robot"), " AI-Powered Analysis Summary", style = "margin: 0; font-weight: 500;")
                    ),

                    # Instructions
                    div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
                      tags$p(class = "mb-2",
                        icon("info-circle"), " ",
                        strong("How it works:"),
                        " Click the button below to generate an AI-powered summary of your differential expression results."
                      ),
                      tags$p(class = "mb-0", style = "font-size: 0.9em; color: #6c757d;",
                        "The AI will analyze the top differentially expressed proteins and most stable proteins from the current comparison."
                      )
                    ),

                    # Generate button
                    div(style = "text-align: center; margin-bottom: 20px;",
                      actionButton("generate_ai_summary_overview",
                        "🤖 Generate AI Summary",
                        class = "btn-info btn-lg",
                        style = "padding: 12px 30px; font-size: 1.1em;"
                      )
                    ),

                    # Output area
                    uiOutput("ai_summary_output")
                  )
                )
              )
    ),

    nav_panel("QC Trends", icon = icon("chart-bar"),
              # Global sort order control
              div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
                div(style = "display: flex; align-items: center; gap: 15px;",
                  icon("sort", style = "color: #6c757d;"),
                  strong("Sort Order:"),
                  radioButtons("qc_sort_order", NULL,
                    choices = c("Run Order", "Group"),
                    inline = TRUE,
                    selected = "Run Order"
                  ),
                  span(style = "color: #6c757d; font-size: 0.85em;",
                    "(Applies to all metric tabs)")
                )
              ),

              # Metric tabs
              navset_card_tab(
                id = "qc_trends_tabs",

                nav_panel("Precursors",
                  icon = icon("dna"),
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_trend_precursors", "🔍 Fullscreen",
                      class = "btn-outline-secondary btn-sm")
                  ),
                  plotlyOutput("qc_trend_plot_precursors", height = "calc(100vh - 380px)")
                ),

                nav_panel("Proteins",
                  icon = icon("shapes"),
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_trend_proteins", "🔍 Fullscreen",
                      class = "btn-outline-secondary btn-sm")
                  ),
                  plotlyOutput("qc_trend_plot_proteins", height = "calc(100vh - 380px)")
                ),

                nav_panel("MS1 Signal",
                  icon = icon("signal"),
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_trend_ms1", "🔍 Fullscreen",
                      class = "btn-outline-secondary btn-sm")
                  ),
                  plotlyOutput("qc_trend_plot_ms1", height = "calc(100vh - 380px)")
                ),

                nav_panel("Stats Table",
                  icon = icon("table"),
                  DTOutput("r_qc_table")
                )
              )
    ),
    
    nav_panel("QC Plots", icon = icon("chart-line"),
              navset_card_tab(
                id = "qc_subtabs",

                # TAB 1: Normalization Diagnostic (FIRST - MOST IMPORTANT)
                nav_panel("Normalization Diagnostic",
                  icon = icon("stethoscope"),
                  card_body(
                    # Guidance banner (keep dynamic uiOutput)
                    uiOutput("norm_diag_guidance"),

                    # Control row
                    div(style = "display: flex; justify-content: space-between; align-items: center; margin: 10px 0;",
                      div(style = "display: flex; gap: 15px; align-items: center;",
                        uiOutput("diann_norm_status_badge", inline = TRUE),
                        radioButtons("norm_diag_type", NULL,
                          choices = c("Box Plots" = "boxplot", "Density Overlay" = "density"),
                          inline = TRUE
                        )
                      ),
                      actionButton("fullscreen_norm_diag", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),

                    # Plot with viewport height
                    plotlyOutput("norm_diagnostic_plot", height = "calc(100vh - 380px)"),

                    # Expandable help
                    tags$details(
                      tags$summary(style = "cursor: pointer; color: #0d6efd; font-size: 0.9em; margin-top: 10px;",
                        icon("question-circle"), " What am I looking at?"
                      ),
                      div(style = "background-color: #f8f9fa; padding: 12px; border-radius: 5px;
                                   margin-top: 8px; font-size: 0.85em; line-height: 1.6;",
                        tags$h6("Reading this plot"),
                        p("Each box (or density curve) represents one sample's intensity distribution \u2014 ",
                          "essentially, how bright all the detected peptides/proteins are in that sample."),
                        p(strong("Left panel: "), "What DIA-NN gave us. These are the peptide-level intensities ",
                          "after DIA-NN's normalization (if it was enabled). ",
                          strong("Right panel: "), "What our pipeline produced. These are the final protein-level ",
                          "estimates after aggregating peptides and handling missing values."),
                        tags$h6("What 'good' looks like"),
                        p("The boxes (or curves) should sit at roughly the same height across all samples. ",
                          "Small differences are normal. If one sample is dramatically higher or lower than ",
                          "the rest, that sample may be problematic."),
                        tags$h6("What 'bad' looks like"),
                        p("If all the boxes are at very different heights, your samples aren't comparable ",
                          "and the statistical results may not be reliable. The most common cause is that ",
                          "DIA-NN normalization was turned off when the data was processed."),
                        tags$h6("Why doesn't the right panel 'fix' bad data?"),
                        p("Unlike some other tools, this pipeline does not apply its own normalization. ",
                          "The protein quantification step (DPC-Quant) aggregates peptides into proteins and ",
                          "handles missing values, but it ",
                          strong("trusts the input intensities as-is"), ". ",
                          "If the input is unnormalized, the output will be too. ",
                          "Normalization happens in DIA-NN, before the data reaches this tool."),
                        tags$h6("The DIA-NN normalization badge"),
                        p("The badge at the top of this plot tells you whether DIA-NN applied normalization ",
                          "to your data. ",
                          tags$span(class = "badge bg-info", "ON"), " = DIA-NN's RT-dependent normalization was active (recommended). ",
                          tags$span(class = "badge bg-warning", "OFF"), " = Data was exported without normalization. ",
                          tags$span(class = "badge bg-secondary", "Unknown"), " = Couldn't determine (older DIA-NN version or non-standard export).")
                      )
                    )
                  )
                ),

                # TAB 2: DPC Fit
                nav_panel("DPC Fit",
                  icon = icon("chart-scatter"),
                  card_body(
                    div(style = "text-align: right; margin-bottom: 10px;",
                      actionButton("fullscreen_dpc", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),
                    plotOutput("dpc_plot", height = "70vh")
                  )
                ),

                # TAB 3: MDS Plot
                nav_panel("MDS Plot",
                  icon = icon("project-diagram"),
                  card_body(
                    div(style = "text-align: right; margin-bottom: 10px;",
                      actionButton("fullscreen_mds", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),
                    plotOutput("mds_plot", height = "70vh")
                  )
                ),

                # TAB 4: Group Distribution
                nav_panel("Group Distribution",
                  icon = icon("chart-area"),
                  card_body(
                    div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                      selectInput("qc_violin_metric", "Metric:",
                        choices = c("Precursors", "Proteins", "MS1_Signal"),
                        width = "200px"
                      ),
                      actionButton("fullscreen_qc_violin", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),
                    plotlyOutput("qc_group_violin", height = "calc(100vh - 320px)")
                  )
                ),

                # TAB 5: P-value Distribution Diagnostic
                nav_panel("P-value Distribution",
                  icon = icon("chart-column"),
                  # Comparison selector banner
                  div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 12px 15px; border-radius: 8px; margin-bottom: 15px; display: flex; align-items: center; justify-content: space-between; gap: 10px; flex-wrap: wrap;",
                    div(style = "display: flex; align-items: center; gap: 10px;",
                      icon("microscope"),
                      span("Viewing Comparison:", style = "font-weight: 500;")
                    ),
                    div(style = "flex-grow: 1; max-width: 400px;",
                      selectInput("contrast_selector_pvalue", NULL,
                        choices = NULL,
                        width = "100%"
                      )
                    )
                  ),
                  card_body(
                    # Automated contextual guidance (changes based on detected pattern)
                    uiOutput("pvalue_guidance"),

                    # Control row
                    div(style = "display: flex; justify-content: flex-end; align-items: center; margin-bottom: 10px;",
                      actionButton("fullscreen_pvalue_hist", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),

                    # Plot
                    plotOutput("pvalue_histogram", height = "calc(100vh - 450px)"),

                    # Expandable interpretation guide
                    tags$details(
                      tags$summary(style = "cursor: pointer; color: #0d6efd; font-size: 0.9em; margin-top: 10px;",
                        icon("question-circle"), " How do I interpret this?"
                      ),
                      div(style = "background-color: #f8f9fa; padding: 12px; border-radius: 5px;
                                   margin-top: 8px; font-size: 0.85em; line-height: 1.6;",
                        tags$h6("What this plot shows"),
                        p("This histogram displays the distribution of raw (unadjusted) p-values from your differential expression analysis. ",
                          "Each bar represents how many proteins have p-values falling in that range."),

                        tags$h6("What 'good' looks like"),
                        tags$ul(
                          tags$li(strong("Flat with a spike at zero: "), "Most p-values uniformly distributed (flat histogram) with a peak near p=0. ",
                            "This indicates a mix of non-changing proteins (uniform) and true positives (spike at zero)."),
                          tags$li(strong("Expected under the null: "), "For proteins that are truly not changing, p-values should be uniformly distributed between 0 and 1. ",
                            "The dashed red line shows this expected uniform distribution.")
                        ),

                        tags$h6("Warning signs"),
                        tags$ul(
                          tags$li(strong("Too many intermediate p-values (0.3-0.7): "), "May indicate p-value inflation due to unmodeled variance, batch effects, or outliers."),
                          tags$li(strong("Depletion near zero: "), "Too few small p-values suggests the test is overly conservative or lacks statistical power."),
                          tags$li(strong("U-shaped distribution: "), "Enrichment at both ends (near 0 and 1) can indicate problems with the statistical model or data quality."),
                          tags$li(strong("Completely uniform: "), "No enrichment at p=0 means no differential expression detected, or the test has no power.")
                        ),

                        tags$h6("What to do if it looks wrong"),
                        tags$ul(
                          tags$li("Check the Normalization Diagnostic tab to ensure samples are properly normalized"),
                          tags$li("Review the MDS plot for outlier samples or unwanted variation"),
                          tags$li("Consider adding batch or other covariates to the model if appropriate"),
                          tags$li("Verify that sample sizes are adequate for the comparison")
                        )
                      )
                    )
                  )
                )
              )
    ),
    
    nav_panel("DE Dashboard", icon = icon("table-columns"),
              # Interactive comparison selector
              div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 8px; margin-bottom: 15px;",
                div(style = "display: flex; align-items: center; gap: 15px;",
                  icon("microscope"),
                  span("Viewing Comparison:", style = "font-weight: 500;"),
                  selectInput("contrast_selector",
                    label = NULL,
                    choices = NULL,
                    width = "300px"
                  )
                )
              ),

              # Responsive two-column grid (stacks on small screens)
              div(class = "de-dashboard-grid",
                # Results table card (left column)
                card(
                  card_header(
                    div(style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px;",
                      span("Results Table"),
                      div(
                        actionButton("clear_plot_selection", "Reset", class="btn-warning btn-xs"),
                        actionButton("show_violin", "📊 Violin", class="btn-primary btn-xs"),
                        actionButton("show_xic", "📈 XICs", class="btn-info btn-xs"),
                        downloadButton("download_result_csv", "💾 Export", class="btn-success btn-xs")
                      )
                    )
                  ),
                  card_body(
                    class = "card-body-scroll",
                    DTOutput("de_table")
                  )
                ),

                # Volcano plot card (right column)
                card(
                  card_header(
                    div(style="display: flex; justify-content: space-between; align-items: center;",
                      span("Volcano Plot (Click/Box Select to Filter Table)"),
                      actionButton("fullscreen_volcano", "🔍 Fullscreen", class="btn-outline-secondary btn-sm")
                    )
                  ),
                  card_body(
                    plotlyOutput("volcano_plot_interactive", height = "calc(100vh - 380px)")
                  )
                )
              ),

              # Heatmap as accordion (secondary visualization, collapsed by default)
              accordion(
                id = "de_heatmap_accordion",
                open = FALSE,
                accordion_panel(
                  "Heatmap of Selected/Top Proteins",
                  icon = icon("grip"),
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_heatmap", "🔍 Fullscreen", class="btn-outline-secondary btn-sm")
                  ),
                  plotOutput("heatmap_plot", height = "450px")
                )
              )
    ),
    
    nav_panel("Consistent DE", icon = icon("check-double"),
              navset_card_tab(
                id = "consistent_de_tabs",

                # TAB 1: High-Consistency Table
                nav_panel("High-Consistency Table",
                  icon = icon("table"),
                  card_body(
                    p("Ranking by %CV (Coefficient of Variation) to find stable markers across all experimental groups.",
                      class = "text-muted small mb-3"),
                    DTOutput("consistent_table")
                  )
                ),

                # TAB 2: CV Distribution Histogram
                nav_panel("CV Distribution",
                  icon = icon("chart-bar"),
                  card_body(
                    # Control row
                    div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                      p("Distribution of Coefficient of Variation (CV) for significant proteins, broken down by experimental group.",
                        class = "text-muted small mb-0"),
                      actionButton("fullscreen_cv_hist", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),
                    plotOutput("cv_histogram", height = "calc(100vh - 320px)")
                  )
                )
              )
    ),
    
    nav_panel("Reproducibility", icon = icon("code"),
              navset_card_tab(
                nav_panel("Code Log",
                          card_body(
                            div(style="background-color: #d1ecf1; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
                              icon("info-circle"),
                              strong(" Action Log:"),
                              "This code recreates your analysis step-by-step. Each section shows:",
                              tags$ul(
                                tags$li(strong("Action name"), " - what you did (e.g., 'Run Pipeline')"),
                                tags$li(strong("Timestamp"), " - when you did it"),
                                tags$li(strong("R code"), " - how to reproduce it")
                              ),
                              p(style="margin-bottom: 0;", "Copy this entire code block to reproduce your analysis in a fresh R session.")
                            ),
                            downloadButton("download_repro_log", "💾 Download Reproducibility Log", class="btn-success mb-3"),
                            verbatimTextOutput("reproducible_code")
                          )
                ),
                nav_panel("Methodology",
                          card_body(
                            verbatimTextOutput("methodology_text")
                          )
                )
              )
    ),
    
    nav_panel("Gene Set Enrichment", icon = icon("sitemap"),
              # Compact control bar
              card(
                card_body(
                  div(style = "display: flex; align-items: center; gap: 15px; flex-wrap: wrap;",
                    actionButton("run_gsea", "▶ Run GSEA", class = "btn-success", icon = icon("play")),
                    div(style = "flex-grow: 1;",
                      verbatimTextOutput("gsea_status", placeholder = TRUE) |>
                        tagAppendAttributes(style = "margin: 0; padding: 5px 10px; min-height: 38px;")
                    )
                  ),
                  p("Performs Gene Ontology enrichment analysis on DE results. Auto-detects organism (Human/Mouse).",
                    class = "text-muted small", style = "margin: 10px 0 0 0;")
                )
              ),

              # Results tabs with full-height plots
              navset_card_tab(
                id = "gsea_results_tabs",

                nav_panel("Dot Plot",
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_gsea_dot", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
                  ),
                  plotOutput("gsea_dot_plot", height = "calc(100vh - 340px)")
                ),

                nav_panel("Enrichment Map",
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_gsea_emap", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
                  ),
                  plotOutput("gsea_emapplot", height = "calc(100vh - 340px)")
                ),

                nav_panel("Ridgeplot",
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_gsea_ridge", "🔍 Fullscreen", class = "btn-outline-secondary btn-sm")
                  ),
                  plotOutput("gsea_ridgeplot", height = "calc(100vh - 340px)")
                ),

                nav_panel("Results Table",
                  DTOutput("gsea_results_table")
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
    ),

    nav_panel("Education", icon = icon("graduation-cap"),
              card(
                card_header("Proteomics Resources & Training"),
                card_body(
                  tags$iframe(
                    src = "https://bsphinney.github.io/DE-LIMP/",
                    style = "width: 100%; height: 700px; border: 1px solid #e2e8f0; border-radius: 8px;",
                    frameborder = "0",
                    scrolling = "yes"
                  ),
                  p("Explore video tutorials, training courses, and methodology citations.",
                    style = "margin-top: 10px; color: #718096; font-size: 0.9em; text-align: center;")
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
    repro_log = c(
      "# ==============================================================================",
      "# DE-LIMP Reproducibility Log",
      sprintf("# Session started: %s", Sys.time()),
      "# ==============================================================================",
      "",
      "# --- Load Required Libraries ---",
      "library(limpa); library(limma); library(dplyr); library(stringr); library(ggrepel);"
    ),
    color_plot_by_de = FALSE,
    grid_selected_protein = NULL,
    temp_violin_target = NULL, # Added for violin popup
    diann_norm_detected = "unknown", # "on", "off", or "unknown" — DIA-NN normalization status

    # === XIC Viewer ===
    xic_dir = NULL,              # Path to directory containing .xic.parquet files
    xic_available = FALSE,       # Whether XIC files were detected
    xic_format = "v2",           # "v1" (DIA-NN 1.x wide) or "v2" (DIA-NN 2.x long)
    xic_protein = NULL,          # Currently selected protein for XIC viewing
    xic_data = NULL,             # Loaded XIC data (tibble, only for selected protein)
    xic_report_map = NULL,       # Protein → Precursor mapping from report
    uploaded_report_path = NULL,  # Path to uploaded report.parquet for re-reading
    original_report_name = NULL,  # Original filename of uploaded report (e.g., "report.parquet")
    mobilogram_available = FALSE, # Whether mobilogram files with non-zero data exist
    mobilogram_dir = NULL         # Path to mobilogram directory (same as xic_dir)
  )

  # Helper function to append to reproducibility log
  add_to_log <- function(action_name, code_lines) {
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    header <- c(
      "",
      paste0("# --- ", action_name, " [", timestamp, "] ---")
    )
    values$repro_log <- c(values$repro_log, header, code_lines)
  }

  # Helper: Detect XIC format version from column names
  detect_xic_format <- function(xic_dir) {
    xic_files <- list.files(xic_dir, pattern = "\\.xic\\.parquet$",
                            full.names = TRUE, recursive = TRUE)
    if (length(xic_files) == 0) return("unknown")
    cols <- tryCatch(names(arrow::read_parquet(xic_files[1], as_data_frame = FALSE)),
                     error = function(e) character(0))
    if (all(c("pr", "feature", "rt", "value") %in% cols)) return("v2")
    if ("Precursor.Id" %in% cols) return("v1")
    return("unknown")
  }

  # Helper: Load XIC data for a single protein using Arrow predicate pushdown
  # Supports both DIA-NN 1.x (wide format) and 2.x (long format: pr/feature/rt/value)
  load_xic_for_protein <- function(xic_dir, protein_id, report_map, xic_format = "v2") {
    # In DIA-NN 2.x, Precursor.Id in the report already matches the XIC "pr" column
    # (both use StrippedSequence+Charge format, e.g., "PEPTIDEK2")
    # In DIA-NN 1.x, Precursor.Id matches the XIC "Precursor.Id" column directly
    target_prs <- report_map %>%
      filter(Protein.Group == protein_id) %>%
      distinct(Precursor.Id) %>%
      pull(Precursor.Id)

    if (length(target_prs) == 0) return(NULL)

    # Read only .xic.parquet files (exclude mobilograms)
    xic_files <- list.files(xic_dir, pattern = "\\.xic\\.parquet$",
                            full.names = TRUE, recursive = TRUE)

    # Per-file read with source file tagging (needed for v2 which has no File.Name)
    xic_list <- lapply(xic_files, function(f) {
      tryCatch({
        df <- arrow::read_parquet(f, as_data_frame = TRUE)
        # Filter by precursor — use base R to avoid rlang issues
        if (xic_format == "v2") {
          df <- df[df$pr %in% target_prs, , drop = FALSE]
        } else {
          df <- df[df$Precursor.Id %in% target_prs, , drop = FALSE]
        }
        if (nrow(df) > 0) {
          sample_name <- sub("\\.xic\\.parquet$", "", basename(f))
          df$Source.File <- sample_name
        }
        df
      }, error = function(e2) NULL)
    })
    result <- bind_rows(Filter(function(x) !is.null(x) && nrow(x) > 0, xic_list))
    if (nrow(result) == 0) return(NULL)
    result
  }

  # Helper: Reshape XIC data for plotting — handles both v1 and v2 formats
  reshape_xic_for_plotting <- function(xic_raw, metadata, xic_format = "v2") {

    if (xic_format == "v2") {
      # ── DIA-NN 2.x: already long format (pr, feature, info, rt, value) ──
      # Rename columns using base R to avoid tidy evaluation issues with "pr"
      names(xic_raw)[names(xic_raw) == "pr"] <- "Precursor.Id"
      names(xic_raw)[names(xic_raw) == "feature"] <- "Fragment.Label"
      names(xic_raw)[names(xic_raw) == "rt"] <- "RT"
      names(xic_raw)[names(xic_raw) == "value"] <- "Intensity"

      xic_plot <- xic_raw %>%
        mutate(
          MS.Level = ifelse(Fragment.Label == "ms1", 1L, 2L),
          RT = as.numeric(RT),
          Intensity = as.numeric(Intensity)
        ) %>%
        filter(!(Intensity == 0 & RT == 0))

      # Match Source.File to metadata File.Name via fuzzy basename matching
      if ("Source.File" %in% names(xic_plot)) {
        meta_lookup <- metadata %>%
          mutate(File.Name.Base = basename(tools::file_path_sans_ext(File.Name)))
        xic_plot <- xic_plot %>%
          mutate(File.Name = Source.File,
                 File.Name.Base = Source.File) %>%
          left_join(
            meta_lookup %>% dplyr::select(File.Name.Base, Group, ID),
            by = "File.Name.Base"
          )
      }

      # Keep only the columns we need — use base R to avoid arrow::select conflict
      keep_cols <- intersect(c("File.Name", "ID", "Group", "Precursor.Id",
                               "MS.Level", "Fragment.Label", "RT", "Intensity"),
                             names(xic_plot))
      xic_plot <- xic_plot[, keep_cols, drop = FALSE]

      return(xic_plot)

    } else {
      # ── DIA-NN 1.x: wide format with numbered columns ──
      num_cols <- names(xic_raw)[grepl("^\\d+$", names(xic_raw))]
      if (length(num_cols) == 0) {
        warning("No numbered columns found in XIC data")
        return(NULL)
      }

      make_key <- function(df) {
        paste(df$File.Name, df$Precursor.Id, df$MS.Level,
              df$Theoretical.Mz, df$FragmentType, df$FragmentCharge,
              df$FragmentSeriesNumber, df$FragmentLossType, sep = "|")
      }

      rt_rows <- xic_raw %>% filter(Retention.Times == 1)
      int_rows <- xic_raw %>% filter(Intensities == 1)

      rt_long <- rt_rows %>%
        mutate(.key = make_key(rt_rows)) %>%
        dplyr::select(.key, all_of(num_cols)) %>%
        pivot_longer(cols = all_of(num_cols), names_to = "point_idx", values_to = "RT")

      int_long <- int_rows %>%
        mutate(.key = make_key(int_rows)) %>%
        dplyr::select(.key, File.Name, Precursor.Id, Modified.Sequence, MS.Level,
               Theoretical.Mz, Reference.Intensity, FragmentType, FragmentCharge,
               FragmentSeriesNumber, FragmentLossType, all_of(num_cols)) %>%
        pivot_longer(cols = all_of(num_cols), names_to = "point_idx", values_to = "Intensity")

      xic_plot <- inner_join(
        rt_long %>% dplyr::select(.key, point_idx, RT),
        int_long,
        by = c(".key", "point_idx")
      ) %>%
        mutate(
          RT = as.numeric(RT),
          Intensity = as.numeric(Intensity),
          Fragment.Label = case_when(
            MS.Level == 1 ~ paste0("MS1 (", round(as.numeric(Theoretical.Mz), 2), ")"),
            TRUE ~ paste0(FragmentType, FragmentSeriesNumber,
                          ifelse(as.integer(FragmentCharge) > 1,
                                 paste0("+", FragmentCharge), ""),
                          ifelse(FragmentLossType != "noloss",
                                 paste0("-", FragmentLossType), ""))
          )
        ) %>%
        filter(!(Intensity == 0 & RT == 0)) %>%
        mutate(File.Name.Base = basename(tools::file_path_sans_ext(
          tools::file_path_sans_ext(File.Name)))) %>%
        left_join(
          metadata %>%
            mutate(File.Name.Base = basename(tools::file_path_sans_ext(File.Name))) %>%
            dplyr::select(File.Name.Base, Group, ID),
          by = "File.Name.Base"
        ) %>%
        dplyr::select(File.Name, ID, Group, Precursor.Id, Modified.Sequence,
               MS.Level, Fragment.Label, Theoretical.Mz, Reference.Intensity,
               RT, Intensity)

      return(xic_plot)
    }
  }

  # ============================================================================
  #      2. Main Data Loading & Processing Pipeline
  # ============================================================================

  # Load example data from GitHub releases
  observeEvent(input$load_example, {
    withProgress(message = "Downloading example data...", {
      example_url <- "https://github.com/bsphinney/DE-LIMP/releases/download/v1.0/Affinisep_vs_evosep_noNorm.parquet"
      temp_file <- tempfile(fileext = ".parquet")

      tryCatch({
        incProgress(0.3, detail = "Downloading from GitHub...")
        download.file(example_url, temp_file, mode = "wb", quiet = TRUE)

        incProgress(0.5, detail = "Calculating Trends...")
        values$qc_stats <- get_diann_stats_r(temp_file)

        incProgress(0.7, detail = "Reading Matrix...")
        values$raw_data <- limpa::readDIANN(temp_file, format="parquet", q.cutoffs=input$q_cutoff)
        fnames <- sort(colnames(values$raw_data$E))
        values$metadata <- data.frame(
          ID = 1:length(fnames),
          File.Name = fnames,
          Group = rep("", length(fnames)),
          Batch = rep("", length(fnames)),
          Covariate1 = rep("", length(fnames)),
          Covariate2 = rep("", length(fnames)),
          stringsAsFactors=FALSE
        )

        # Initialize custom covariate names
        if(is.null(values$cov1_name)) values$cov1_name <- "Covariate1"
        if(is.null(values$cov2_name)) values$cov2_name <- "Covariate2"

        # Flag this as example data for auto-guess logic
        values$is_example_data <- TRUE

        # Detect DIA-NN normalization status
        values$diann_norm_detected <- tryCatch({
          raw_parquet <- arrow::read_parquet(temp_file,
            col_select = c("Precursor.Quantity", "Precursor.Normalised"))
          has_both_cols <- all(c("Precursor.Quantity", "Precursor.Normalised") %in% names(raw_parquet))
          if (has_both_cols) {
            sample_rows <- head(raw_parquet, 1000)
            ratio <- sample_rows$Precursor.Normalised / sample_rows$Precursor.Quantity
            ratios_vary <- sd(ratio, na.rm = TRUE) > 0.001
            if (ratios_vary) "on" else "off"
          } else { "unknown" }
        }, error = function(e) "unknown")

        # Store report path for XIC viewer precursor mapping (copy to session dir)
        session_report <- file.path(tempdir(), "de_limp_report.parquet")
        file.copy(temp_file, session_report, overwrite = TRUE)
        values$uploaded_report_path <- session_report
        values$original_report_name <- "Affinisep_vs_evosep_noNorm.parquet"

        # Auto-detect XIC directory in working directory for example data
        tryCatch({
          cand <- file.path(getwd(), "Affinisep_vs_evosep_noNorm_xic")
          if (dir.exists(cand) && length(list.files(cand, pattern = "\\.xic\\.parquet$")) > 0) {
            updateTextInput(session, "xic_dir_input", value = cand)
            showNotification("Auto-detected XIC directory for example data.",
              type = "message", duration = 4)
          }
        }, error = function(e) NULL)

        incProgress(0.9, detail = "Opening setup...")

        # Log to reproducibility
        add_to_log("Example Data Loaded", c(
          "# Example data: Affinisep vs Evosep (50ng Thermo Hela digest)",
          sprintf("# Downloaded from: %s", example_url),
          sprintf("dat <- readDIANN('Affinisep_vs_evosep_noNorm.parquet', format='parquet', q.cutoffs=%s)", input$q_cutoff)
        ))

        showNotification("Example data loaded successfully!", type = "message", duration = 3)
        # Navigate to Assign Groups sub-tab
        nav_select("main_tabs", "Data Overview")
        nav_select("data_overview_tabs", "Assign Groups & Run")

      }, error = function(e) {
        showNotification(paste("Error loading example data:", e$message), type = "error", duration = 10)
      })
    })
  })

  observeEvent(input$report_file, {
    req(input$report_file)
    withProgress(message = "Loading...", {
      incProgress(0.2, detail = "Calculating Trends...")
      values$qc_stats <- get_diann_stats_r(input$report_file$datapath)
      incProgress(0.5, detail = "Reading Matrix...")
      tryCatch({
        values$raw_data <- limpa::readDIANN(input$report_file$datapath, format="parquet", q.cutoffs=input$q_cutoff)
        fnames <- sort(colnames(values$raw_data$E))
        values$metadata <- data.frame(
          ID = 1:length(fnames),
          File.Name = fnames,
          Group = rep("", length(fnames)),
          Batch = rep("", length(fnames)),
          Covariate1 = rep("", length(fnames)),
          Covariate2 = rep("", length(fnames)),
          stringsAsFactors=FALSE
        )
        # Initialize custom covariate names (user can change these)
        if(is.null(values$cov1_name)) values$cov1_name <- "Covariate1"
        if(is.null(values$cov2_name)) values$cov2_name <- "Covariate2"
        # Clear example data flag for user uploads
        values$is_example_data <- FALSE

        # Detect DIA-NN normalization status
        values$diann_norm_detected <- tryCatch({
          raw_parquet <- arrow::read_parquet(input$report_file$datapath,
            col_select = c("Precursor.Quantity", "Precursor.Normalised"))
          has_both_cols <- all(c("Precursor.Quantity", "Precursor.Normalised") %in% names(raw_parquet))
          if (has_both_cols) {
            sample_rows <- head(raw_parquet, 1000)
            ratio <- sample_rows$Precursor.Normalised / sample_rows$Precursor.Quantity
            ratios_vary <- sd(ratio, na.rm = TRUE) > 0.001
            if (ratios_vary) "on" else "off"
          } else { "unknown" }
        }, error = function(e) "unknown")

        # Store report path for XIC viewer precursor mapping (copy to session dir)
        session_report <- file.path(tempdir(), "de_limp_report.parquet")
        file.copy(input$report_file$datapath, session_report, overwrite = TRUE)
        values$uploaded_report_path <- session_report
        values$original_report_name <- input$report_file$name

        # Auto-detect XIC directory next to the uploaded report
        # For local Shiny, check if _xic sibling directory exists
        tryCatch({
          report_name <- tools::file_path_sans_ext(input$report_file$name)
          # Check common locations: working directory, or if user uploaded from a known path
          candidate_dirs <- c(
            file.path(getwd(), paste0(report_name, "_xic")),
            file.path(dirname(input$report_file$datapath), paste0(report_name, "_xic"))
          )
          for (cand in candidate_dirs) {
            if (dir.exists(cand) && length(list.files(cand, pattern = "\\.xic\\.parquet$")) > 0) {
              updateTextInput(session, "xic_dir_input", value = cand)
              showNotification(paste0("Auto-detected XIC directory: ", basename(cand)),
                type = "message", duration = 4)
              break
            }
          }
        }, error = function(e) NULL)

        # Navigate to Assign Groups sub-tab
        nav_select("main_tabs", "Data Overview")
        nav_select("data_overview_tabs", "Assign Groups & Run")
      }, error=function(e) { showNotification(paste("Error:", e$message), type="error") })
    })
  })

  observeEvent(input$report_file, {
    req(input$report_file)
    add_to_log("Data Upload", c(
      sprintf("# File: %s", input$report_file$name),
      sprintf("dat <- readDIANN('%s', format='parquet', q.cutoffs=%s)",
              "path/to/your/report.parquet", input$q_cutoff)
    ))
  })

  # Old standalone "Run Pipeline" button observer - REMOVED
  # Pipeline now runs from the "Assign Groups & Run" sub-tab via run_pipeline observer
  
  # ============================================================================
  #      3. Metadata Handling (Assign Groups sub-tab)
  # ============================================================================

  output$hot_metadata <- renderRHandsontable({
    req(values$metadata)

    # Get custom covariate names
    cov1_display <- if(!is.null(input$cov1_label) && input$cov1_label != "") input$cov1_label else "Covariate1"
    cov2_display <- if(!is.null(input$cov2_label) && input$cov2_label != "") input$cov2_label else "Covariate2"

    # Store for later use
    values$cov1_name <- cov1_display
    values$cov2_name <- cov2_display

    # Create display dataframe with custom column names
    display_df <- values$metadata
    colnames(display_df) <- c("ID", "File.Name", "Group", "Batch", cov1_display, cov2_display)

    rhandsontable(display_df, rowHeaders=NULL, stretchH="all", height=500, width="100%") %>%
      hot_col("ID", readOnly=TRUE, width=50) %>%
      hot_col("File.Name", readOnly=TRUE) %>%
      hot_col("Group", type="text") %>%
      hot_col("Batch", type="text", width=100) %>%
      hot_col(cov1_display, type="text", width=100) %>%
      hot_col(cov2_display, type="text", width=100)
  })

  observeEvent(input$guess_groups, {
    req(values$metadata)
    meta <- if(!is.null(input$hot_metadata)) hot_to_r(input$hot_metadata) else values$metadata

    cleaned_filenames <- str_remove(meta$File.Name, "^\\d{8}_")
    cleaned_filenames <- str_remove_all(cleaned_filenames, "_deepfrac|\\.parquet")
    keywords <- c("affinisepACN", "affinisepIPA", "Control", "Treatment", "Evosep", "Affinisep")

    find_best_match <- function(filename_clean) {
      matches <- keywords[str_detect(filename_clean, regex(keywords, ignore_case = TRUE))]
      if (length(matches) == 0) return("") else return(matches[which.max(nchar(matches))])
    }

    guessed_groups <- sapply(cleaned_filenames, find_best_match)

    # Find indices of unmatched samples
    unmatched_indices <- which(guessed_groups == "")

    # Use flag to detect if this is the example data
    if (!is.null(values$is_example_data) && values$is_example_data && length(unmatched_indices) >= 3) {
      # For example data: set last 3 unmatched samples to "Evosep"
      last_three <- tail(unmatched_indices, 3)
      guessed_groups[last_three] <- "Evosep"
      # Set any remaining unmatched to generic Sample_X
      remaining_unmatched <- setdiff(unmatched_indices, last_three)
      sample_counter <- 0
      for (i in remaining_unmatched) {
        sample_counter <- sample_counter + 1
        guessed_groups[i] <- paste0("Sample_", sample_counter)
      }
    } else {
      # Standard behavior: all unmatched become Sample_X
      sample_counter <- 0
      for (i in seq_along(guessed_groups)) {
        if (guessed_groups[i] == "") {
          sample_counter <- sample_counter + 1
          guessed_groups[i] <- paste0("Sample_", sample_counter)
        }
      }
    }

    meta$Group <- guessed_groups
    values$metadata <- meta
  })

  # Run pipeline from "Assign Groups & Run" sub-tab - saves groups and runs pipeline
  observeEvent(input$run_pipeline, {
    req(input$hot_metadata, values$metadata, values$raw_data)

    # First, save the groups
    old_meta <- values$metadata
    new_meta <- hot_to_r(input$hot_metadata)
    changed_indices <- which(old_meta$Group != new_meta$Group)
    if (length(changed_indices) > 0) {
      code_lines <- sprintf("metadata$Group[%d] <- '%s'  # %s",
                           changed_indices,
                           new_meta$Group[changed_indices],
                           new_meta$File.Name[changed_indices])
      add_to_log("Manual Group Assignment", code_lines)
    }
    values$metadata <- new_meta

    # Validate groups
    meta <- values$metadata
    meta$Group <- trimws(meta$Group)
    if(length(unique(meta$Group)) < 2) {
      showNotification("Error: Need at least 2 groups to run pipeline.", type="error", duration = 10)
      return()
    }

    showNotification("Groups saved! Running pipeline...", type="message")

    # Build covariates list for logging
    covariates_to_log <- character(0)
    cov_display_names <- character(0)

    if (isTRUE(input$include_batch) && length(unique(meta$Batch[meta$Batch != ""])) > 1) {
      covariates_to_log <- c(covariates_to_log, "Batch")
      cov_display_names <- c(cov_display_names, "Batch")
    }
    if (isTRUE(input$include_cov1) && length(unique(meta$Covariate1[meta$Covariate1 != ""])) > 1) {
      covariates_to_log <- c(covariates_to_log, "Covariate1")
      cov_display_names <- c(cov_display_names, values$cov1_name %||% "Covariate1")
    }
    if (isTRUE(input$include_cov2) && length(unique(meta$Covariate2[meta$Covariate2 != ""])) > 1) {
      covariates_to_log <- c(covariates_to_log, "Covariate2")
      cov_display_names <- c(cov_display_names, values$cov2_name %||% "Covariate2")
    }

    # Generate pipeline code for reproducibility log
    pipeline_code <- c(
      "# Normalization & Quantification",
      "dpcfit <- dpcCN(dat)",
      "y_protein <- dpcQuant(dat, 'Protein.Group', dpc=dpcfit)",
      ""
    )

    if (length(covariates_to_log) > 0) {
      # Log with covariates
      pipeline_code <- c(pipeline_code,
        sprintf("# Experimental Design (with covariates: %s)", paste(covariates_to_log, collapse = ", ")),
        "group_map <- c(",
        paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Group), collapse=",\n"),
        ")"
      )

      # Add covariate maps
      if ("Batch" %in% covariates_to_log) {
        pipeline_code <- c(pipeline_code,
          "batch_map <- c(",
          paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Batch), collapse=",\n"),
          ")"
        )
      }
      if ("Covariate1" %in% covariates_to_log) {
        cov1_name <- values$cov1_name %||% "Covariate1"
        pipeline_code <- c(pipeline_code,
          sprintf("%s_map <- c(", tolower(gsub(" ", "_", cov1_name))),
          paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Covariate1), collapse=",\n"),
          ")"
        )
      }
      if ("Covariate2" %in% covariates_to_log) {
        cov2_name <- values$cov2_name %||% "Covariate2"
        pipeline_code <- c(pipeline_code,
          sprintf("%s_map <- c(", tolower(gsub(" ", "_", cov2_name))),
          paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Covariate2), collapse=",\n"),
          ")"
        )
      }

      # Build metadata dataframe
      df_cols <- "Group = group_map"
      if ("Batch" %in% covariates_to_log) df_cols <- paste0(df_cols, ", Batch = batch_map")
      if ("Covariate1" %in% covariates_to_log) {
        cov1_name <- values$cov1_name %||% "Covariate1"
        cov1_var <- tolower(gsub(" ", "_", cov1_name))
        df_cols <- paste0(df_cols, sprintf(", %s = %s_map", cov1_name, cov1_var))
      }
      if ("Covariate2" %in% covariates_to_log) {
        cov2_name <- values$cov2_name %||% "Covariate2"
        cov2_var <- tolower(gsub(" ", "_", cov2_name))
        df_cols <- paste0(df_cols, sprintf(", %s = %s_map", cov2_name, cov2_var))
      }

      pipeline_code <- c(pipeline_code,
        sprintf("metadata <- data.frame(File.Name = names(group_map), %s)", df_cols),
        "metadata <- metadata[match(colnames(dat$E), metadata$File.Name), ]",
        "groups <- factor(metadata$Group)"
      )

      # Add factor creation for each covariate
      if ("Batch" %in% covariates_to_log) {
        pipeline_code <- c(pipeline_code, "batch <- factor(metadata$Batch)")
      }
      if ("Covariate1" %in% covariates_to_log) {
        cov1_name <- values$cov1_name %||% "Covariate1"
        cov1_var <- tolower(gsub(" ", "_", cov1_name))
        pipeline_code <- c(pipeline_code, sprintf("%s <- factor(metadata$%s)", cov1_var, cov1_name))
      }
      if ("Covariate2" %in% covariates_to_log) {
        cov2_name <- values$cov2_name %||% "Covariate2"
        cov2_var <- tolower(gsub(" ", "_", cov2_name))
        pipeline_code <- c(pipeline_code, sprintf("%s <- factor(metadata$%s)", cov2_var, cov2_name))
      }

      # Build design formula with custom names
      formula_parts <- c("groups")
      if ("Batch" %in% covariates_to_log) formula_parts <- c(formula_parts, "batch")
      if ("Covariate1" %in% covariates_to_log) formula_parts <- c(formula_parts, tolower(gsub(" ", "_", values$cov1_name %||% "covariate1")))
      if ("Covariate2" %in% covariates_to_log) formula_parts <- c(formula_parts, tolower(gsub(" ", "_", values$cov2_name %||% "covariate2")))
      formula_str <- paste0("~ 0 + ", paste(formula_parts, collapse = " + "))
      pipeline_code <- c(pipeline_code,
        sprintf("design <- model.matrix(%s)", formula_str),
        "colnames(design) <- gsub('groups', '', colnames(design))"
      )

      pipeline_code <- c(pipeline_code, "",
        "# Differential Expression Model (with covariates)",
        "fit <- dpcDE(y_protein, design, plot=FALSE)"
      )
    } else {
      # Log without covariates
      pipeline_code <- c(pipeline_code,
        "# Experimental Design",
        "group_map <- c(",
        paste(sprintf("  '%s' = '%s'", meta$File.Name, meta$Group), collapse=",\n"),
        ")",
        "metadata <- data.frame(File.Name = names(group_map), Group = group_map)",
        "metadata <- metadata[match(colnames(dat$E), metadata$File.Name), ]",
        "groups <- factor(metadata$Group)",
        "design <- model.matrix(~ 0 + groups)",
        "colnames(design) <- levels(groups)",
        "",
        "# Differential Expression Model",
        "fit <- dpcDE(y_protein, design, plot=FALSE)"
      )
    }

    add_to_log("Run Pipeline (Main Analysis)", pipeline_code)

    withProgress(message='Running Pipeline...', {
      tryCatch({
        dat <- values$raw_data
        dpcfit <- limpa::dpcCN(dat)
        values$dpc_fit <- dpcfit

        values$y_protein <- tryCatch({
          limpa::dpcQuant(dat, "Protein.Group", dpc=dpcfit)
        }, error = function(e) {
          showNotification(paste("Protein quantification failed:", e$message), type = "error", duration = NULL)
          return(NULL)
        })

        req(values$y_protein)

        rownames(meta) <- meta$File.Name
        meta <- meta[colnames(dat$E), ]
        meta$Group <- make.names(meta$Group)
        groups <- factor(meta$Group)

        # Build design formula with selected covariates
        covariates_to_include <- character(0)
        covariate_messages <- character(0)

        # Check each covariate
        if (isTRUE(input$include_batch) && length(unique(meta$Batch[meta$Batch != ""])) > 1) {
          meta$Batch <- make.names(meta$Batch)
          covariates_to_include <- c(covariates_to_include, "batch")
          covariate_messages <- c(covariate_messages, "Batch")
        }

        if (isTRUE(input$include_cov1) && length(unique(meta$Covariate1[meta$Covariate1 != ""])) > 1) {
          meta$Covariate1 <- make.names(meta$Covariate1)
          cov1_var <- tolower(gsub(" ", "_", values$cov1_name %||% "covariate1"))
          covariates_to_include <- c(covariates_to_include, cov1_var)
          covariate_messages <- c(covariate_messages, values$cov1_name %||% "Covariate1")
        }

        if (isTRUE(input$include_cov2) && length(unique(meta$Covariate2[meta$Covariate2 != ""])) > 1) {
          meta$Covariate2 <- make.names(meta$Covariate2)
          cov2_var <- tolower(gsub(" ", "_", values$cov2_name %||% "covariate2"))
          covariates_to_include <- c(covariates_to_include, cov2_var)
          covariate_messages <- c(covariate_messages, values$cov2_name %||% "Covariate2")
        }

        # Build design matrix
        if (length(covariates_to_include) > 0) {
          # Create formula dynamically
          formula_str <- paste0("~ 0 + groups + ", paste(covariates_to_include, collapse = " + "))

          # Build data frame for design matrix
          design_df <- data.frame(groups = groups)
          if("batch" %in% covariates_to_include) design_df$batch <- factor(meta$Batch)
          if(isTRUE(input$include_cov1) && length(unique(meta$Covariate1[meta$Covariate1 != ""])) > 1) {
            cov1_var <- tolower(gsub(" ", "_", values$cov1_name %||% "covariate1"))
            design_df[[cov1_var]] <- factor(meta$Covariate1)
          }
          if(isTRUE(input$include_cov2) && length(unique(meta$Covariate2[meta$Covariate2 != ""])) > 1) {
            cov2_var <- tolower(gsub(" ", "_", values$cov2_name %||% "covariate2"))
            design_df[[cov2_var]] <- factor(meta$Covariate2)
          }

          design <- model.matrix(as.formula(formula_str), data = design_df)
          colnames(design) <- gsub("groups", "", colnames(design))

          showNotification(
            paste0("Including covariates: ", paste(covariate_messages, collapse = ", ")),
            type = "message", duration = 5
          )
        } else {
          # Standard design without covariates
          design <- model.matrix(~ 0 + groups)
          colnames(design) <- levels(groups)
        }

        combs <- combn(levels(groups), 2)
        forms <- apply(combs, 2, function(x) paste(x[2], "-", x[1]))

        fit <- limpa::dpcDE(values$y_protein, design, plot=FALSE)
        fit <- contrasts.fit(fit, makeContrasts(contrasts=forms, levels=design))
        fit <- eBayes(fit)
        values$fit <- fit

        # Update all four comparison selectors
        updateSelectInput(session, "contrast_selector", choices=forms)
        updateSelectInput(session, "contrast_selector_signal", choices=forms, selected=forms[1])
        updateSelectInput(session, "contrast_selector_grid", choices=forms, selected=forms[1])
        updateSelectInput(session, "contrast_selector_pvalue", choices=forms, selected=forms[1])
        values$status <- "✅ Complete!"

        # Log contrasts
        contrast_code <- c(
          sprintf("# Available contrasts: %s", paste(forms, collapse=", ")),
          "combs <- combn(levels(groups), 2)",
          "forms <- apply(combs, 2, function(x) paste(x[2], '-', x[1]))",
          "fit <- contrasts.fit(fit, makeContrasts(contrasts=forms, levels=design))",
          "fit <- eBayes(fit)"
        )
        add_to_log("Contrast Fitting", contrast_code)

        nav_select("main_tabs", "QC Plots")
        showNotification("✓ Pipeline complete! View results in tabs below.", type="message", duration=10)
      }, error = function(e) {
        showNotification(paste("Pipeline error:", e$message), type = "error", duration = NULL)
      })
    })
  })

  # ============================================================================
  #      4. Parameter Change Logging
  # ============================================================================

  # Log contrast changes
  observeEvent(input$contrast_selector, {
    req(input$contrast_selector)
    if (!is.null(values$fit)) {
      add_to_log("Select Contrast for Visualization", c(
        sprintf("# Viewing contrast: %s", input$contrast_selector),
        sprintf("results <- topTable(fit, coef='%s', number=Inf)", input$contrast_selector)
      ))
    }

    # Sync with Signal Distribution, Expression Grid, and P-value Distribution selectors
    if (!is.null(input$contrast_selector_signal) && input$contrast_selector_signal != input$contrast_selector) {
      updateSelectInput(session, "contrast_selector_signal", selected = input$contrast_selector)
    }
    if (!is.null(input$contrast_selector_grid) && input$contrast_selector_grid != input$contrast_selector) {
      updateSelectInput(session, "contrast_selector_grid", selected = input$contrast_selector)
    }
    if (!is.null(input$contrast_selector_pvalue) && input$contrast_selector_pvalue != input$contrast_selector) {
      updateSelectInput(session, "contrast_selector_pvalue", selected = input$contrast_selector)
    }
  })

  # Sync Signal Distribution selector with main selector
  observeEvent(input$contrast_selector_signal, {
    req(input$contrast_selector_signal)
    if (!is.null(input$contrast_selector) && input$contrast_selector != input$contrast_selector_signal) {
      updateSelectInput(session, "contrast_selector", selected = input$contrast_selector_signal)
    }
  })

  # Sync Expression Grid selector with main selector
  observeEvent(input$contrast_selector_grid, {
    req(input$contrast_selector_grid)
    if (!is.null(input$contrast_selector) && input$contrast_selector != input$contrast_selector_grid) {
      updateSelectInput(session, "contrast_selector", selected = input$contrast_selector_grid)
    }
  })

  # Sync P-value Distribution selector with main selector
  observeEvent(input$contrast_selector_pvalue, {
    req(input$contrast_selector_pvalue)
    if (!is.null(input$contrast_selector) && input$contrast_selector != input$contrast_selector_pvalue) {
      updateSelectInput(session, "contrast_selector", selected = input$contrast_selector_pvalue)
    }
  })

  # Log logFC threshold changes
  observeEvent(input$logfc_cutoff, {
    req(input$logfc_cutoff)
    if (!is.null(values$fit)) {
      add_to_log("Change LogFC Threshold", c(
        sprintf("# LogFC cutoff set to: %.2f", input$logfc_cutoff),
        sprintf("# Filter: adj.P.Val < 0.05 & abs(logFC) > %.2f", input$logfc_cutoff)
      ))
    }
  })

  # ============================================================================
  #      5. Reactive Data Expressions
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
  #      6. UI Outputs & Renderers
  # ============================================================================
  output$run_status_msg <- renderText({ values$status })
  
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
        style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px;",
        tags$h4(icon("flask"), " Differential Expression Summary"),
        tags$hr(),
        tags$p(style = "font-size: 0.85em; color: #6c757d; margin-bottom: 15px;",
          "Proteins with FDR-adjusted p-value < 0.05. Arrows indicate direction of change."
        ),
        do.call(tagList, de_summary_list)
      )
    }

    tagList(summary_elements)
  })

  # --- AI SUMMARY (Data Overview Tab) ---
  observeEvent(input$generate_ai_summary_overview, {
    req(values$fit, input$contrast_selector, values$y_protein, input$user_api_key)

    withProgress(message = "Generating AI Summary...", value = 0, {
      incProgress(0.2, detail = "Gathering DE data...")

      # Get volcano data for the current contrast
      de_results_full <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>% as.data.frame()
      if (!"Protein.Group" %in% colnames(de_results_full)) {
        de_results_full <- de_results_full %>% rownames_to_column("Protein.Group")
      }

      # Add gene symbols
      de_results_full$Accession <- str_split_fixed(de_results_full$Protein.Group, "[; ]", 2)[,1]
      org_db_name <- detect_organism_db(de_results_full$Protein.Group)

      id_map <- tryCatch({
        if (!requireNamespace(org_db_name, quietly = TRUE)) { BiocManager::install(org_db_name, ask = FALSE) }
        library(org_db_name, character.only = TRUE)
        db_obj <- get(org_db_name)
        AnnotationDbi::select(db_obj, keys = de_results_full$Accession, columns = c("SYMBOL"), keytype = "UNIPROT") %>%
          rename(Accession = UNIPROT, Gene = SYMBOL) %>% distinct(Accession, .keep_all = TRUE)
      }, error = function(e) data.frame(Accession = de_results_full$Accession, Gene = de_results_full$Accession))

      de_results_with_genes <- left_join(de_results_full, id_map, by = "Accession")
      de_results_with_genes$Gene[is.na(de_results_with_genes$Gene)] <- de_results_with_genes$Accession[is.na(de_results_with_genes$Gene)]

      # Filter for significant and get top 50
      top_de_data <- de_results_with_genes %>%
        filter(adj.P.Val < 0.05) %>%
        arrange(adj.P.Val) %>%
        head(50) %>%
        dplyr::select(Gene, logFC, adj.P.Val) %>%
        mutate(across(where(is.numeric), ~round(.x, 3)))

      top_de_text <- paste(capture.output(print(as.data.frame(top_de_data))), collapse = "\n")

      incProgress(0.5, detail = "Gathering stable proteins...")

      stable_prots_df <- tryCatch({
        df_res <- de_results_with_genes %>% filter(adj.P.Val < 0.05)
        if(nrow(df_res) == 0) return(data.frame(Info="No significant proteins to assess for stability."))

        protein_ids_for_cv <- df_res$Protein.Group
        raw_exprs <- values$y_protein$E[protein_ids_for_cv, , drop = FALSE]
        linear_exprs <- 2^raw_exprs
        cv_list <- list()

        for(g in unique(values$metadata$Group)) {
          if (g == "") next
          files_in_group <- values$metadata$File.Name[values$metadata$Group == g]
          group_cols <- intersect(colnames(linear_exprs), files_in_group)
          if (length(group_cols) > 1) {
            group_data <- linear_exprs[, group_cols, drop = FALSE]
            cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100)
          } else {
            cv_list[[paste0("CV_", g)]] <- NA
          }
        }

        cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
        left_join(df_res, cv_df, by = "Protein.Group") %>%
          rowwise() %>%
          mutate(Avg_CV = mean(c_across(starts_with("CV_")), na.rm = TRUE)) %>%
          ungroup() %>%
          arrange(Avg_CV) %>%
          head(3) %>%
          dplyr::select(Gene, Avg_CV, logFC, adj.P.Val) %>%
          mutate(across(where(is.numeric), ~round(.x, 2)))
      }, error = function(e) data.frame(Error = "Could not calculate stable proteins."))

      stable_prots_text <- paste(capture.output(print(as.data.frame(stable_prots_df))), collapse = "\n")

      incProgress(0.7, detail = "Constructing prompt...")

      system_prompt <- paste0(
        "You are a senior proteomics consultant at a core facility. Write a 3-paragraph summary of the differential expression results below.\n\n",
        "Paragraph 1: Overview of the comparison being analyzed and total number of significant proteins (FDR < 0.05).\n",
        "Paragraph 2: Highlight the top upregulated and downregulated proteins by fold change, using gene names when available.\n",
        "Paragraph 3: Mention the most stable differentially expressed proteins (lowest CV) as potential high-confidence biomarkers.\n\n",
        "Use markdown formatting. Be concise and scientific."
      )

      final_prompt <- paste0(
        system_prompt,
        "\n\n--- DATA FOR SUMMARY ---\n\n",
        "Comparison: ", input$contrast_selector, "\n\n",
        "Top Significant Proteins:\n", top_de_text, "\n\n",
        "Most Stable DE Proteins:\n", stable_prots_text
      )

      incProgress(0.8, detail = "Asking AI...")
      ai_summary <- ask_gemini_text_chat(final_prompt, input$user_api_key, input$model_name)

      # Render the summary to the output area
      output$ai_summary_output <- renderUI({
        div(style = "background-color: #ffffff; padding: 20px; border: 1px solid #dee2e6; border-radius: 8px;",
          tags$h5(icon("check-circle"), " Analysis Complete", style = "color: #28a745; margin-bottom: 15px;"),
          HTML(markdown::markdownToHTML(text = ai_summary, fragment.only = TRUE))
        )
      })
    })
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
      dplyr::select(Protein.Group, Gene, Protein.Name, Significance, logFC, P.Value, adj.P.Val, all_of(valid_cols)) %>%
      mutate(across(where(is.numeric), ~round(., 2)))

    df_final$Original.ID <- df_final$Protein.Group
    fixed_cols <- c("Protein.Group", "Gene", "Protein.Name", "Significance", "logFC", "P.Value", "adj.P.Val")
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
        icon("list"), " Click to view File ID Mapping (Run # → Filename)"
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
      showModal(modalDialog(title = paste("Expression Plot:", selected_id), size = "xl", plotOutput("violin_plot_grid", height = "600px"), footer = tagList(actionButton("show_xic_from_grid", "📈 XICs", class="btn-info"), actionButton("back_to_grid", "Back to Grid", class="btn-info"), modalButton("Close")), easyClose = TRUE))
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

  output$protein_signal_plot <- renderPlot({
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
  
  # Helper function to generate QC trend plot (parameterized by metric)
  generate_qc_trend_plot <- function(metric) {
    reactive({
      req(values$qc_stats, values$metadata)
      df <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")) %>%
        mutate(Run_Number = as.numeric(str_extract(Run, "\\d+$")))

      if (input$qc_sort_order == "Group") {
        df <- df %>% arrange(Group, Run_Number)
      } else {
        df <- df %>% arrange(Run_Number)
      }

      df$Sort_Index <- 1:nrow(df)
      df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Group:</b> ", df$Group,
                           "<br><b>", metric, ":</b> ", round(df[[metric]], 2))

      # Calculate group means and ranges for average lines
      group_stats <- df %>%
        group_by(Group) %>%
        summarise(
          mean_value = mean(.data[[metric]], na.rm = TRUE),
          x_min = min(Sort_Index),
          x_max = max(Sort_Index),
          .groups = 'drop'
        )

      # Create plot with bars and group average lines
      p <- ggplot(df, aes(x = Sort_Index, y = .data[[metric]], fill = Group, text = Tooltip)) +
        geom_bar(stat = "identity", width = 0.8) +
        geom_segment(data = group_stats,
                     aes(x = x_min - 0.5, xend = x_max + 0.5,
                         y = mean_value, yend = mean_value,
                         color = Group),
                     linewidth = 1, linetype = "dashed", inherit.aes = FALSE,
                     show.legend = FALSE) +
        scale_color_discrete(guide = "none") +
        theme_minimal() +
        labs(title = paste(metric, "per Run (dashed lines = group averages)"),
             x = "Sample Index (Sorted)", y = metric) +
        theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size=8))

      ggplotly(p, tooltip = "text") %>% config(displayModeBar = TRUE)
    })
  }

  # Three separate plot outputs for the three metrics
  output$qc_trend_plot_precursors <- renderPlotly({
    generate_qc_trend_plot("Precursors")()
  })

  output$qc_trend_plot_proteins <- renderPlotly({
    generate_qc_trend_plot("Proteins")()
  })

  output$qc_trend_plot_ms1 <- renderPlotly({
    generate_qc_trend_plot("MS1_Signal")()
  })

  # Fullscreen modals for each metric
  observeEvent(input$fullscreen_trend_precursors, {
    showModal(modalDialog(
      title = "Precursors Trend - Fullscreen View",
      plotlyOutput("qc_trend_plot_precursors_fs", height = "700px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$qc_trend_plot_precursors_fs <- renderPlotly({
    generate_qc_trend_plot("Precursors")()
  })

  observeEvent(input$fullscreen_trend_proteins, {
    showModal(modalDialog(
      title = "Proteins Trend - Fullscreen View",
      plotlyOutput("qc_trend_plot_proteins_fs", height = "700px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$qc_trend_plot_proteins_fs <- renderPlotly({
    generate_qc_trend_plot("Proteins")()
  })

  observeEvent(input$fullscreen_trend_ms1, {
    showModal(modalDialog(
      title = "MS1 Signal Trend - Fullscreen View",
      plotlyOutput("qc_trend_plot_ms1_fs", height = "700px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$qc_trend_plot_ms1_fs <- renderPlotly({
    generate_qc_trend_plot("MS1_Signal")()
  })

  output$r_qc_table <- renderDT({ req(values$qc_stats); df_display <- values$qc_stats %>% arrange(Run) %>% mutate(ID = 1:n()) %>% dplyr::select(ID, Run, everything()); datatable(df_display, options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE) })
  
  output$qc_group_violin <- renderPlotly({
    req(values$qc_stats, values$metadata, input$qc_violin_metric)
    df <- left_join(values$qc_stats, values$metadata, by=c("Run"="File.Name")); metric <- input$qc_violin_metric
    df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Val:</b> ", round(df[[metric]], 2))
    p <- ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) + geom_violin(alpha = 0.5, trim = FALSE) + geom_jitter(aes(text = Tooltip), width = 0.2, size = 2, alpha = 0.8, color = "black") + theme_bw() + labs(title = paste("Distribution of", metric), x = "Group", y = metric) + theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })
  
  output$dpc_plot <- renderPlot({ req(values$dpc_fit); limpa::plotDPC(values$dpc_fit) }) # Height controlled by UI (70vh)
  
  output$mds_plot <- renderPlot({
    req(values$y_protein, values$metadata)
    meta <- values$metadata[match(colnames(values$y_protein$E), values$metadata$File.Name), ]; grps <- factor(meta$Group); cols <- rainbow(length(levels(grps)))
    par(xpd = TRUE); limpa::plotMDSUsingSEs(values$y_protein, pch=16, main="MDS Plot", col=cols[grps]); legend(x = "right", inset = c(-0.2, 0), legend=levels(grps), col=cols, pch=16, bty = "n")
  }) # Height controlled by UI (70vh)

  # ============================================================================
  #      Pipeline Diagnostic: Input → Output Distributions
  # ============================================================================

  # DIA-NN normalization status badge
  output$diann_norm_status_badge <- renderUI({
    status <- values$diann_norm_detected
    if (status == "on") {
      span(class = "badge bg-info", style = "margin-right: 10px;",
        icon("check-circle"), " DIA-NN normalization: ON (RT-dependent)")
    } else if (status == "off") {
      span(class = "badge bg-warning", style = "margin-right: 10px;",
        icon("exclamation-triangle"), " DIA-NN normalization: OFF")
    } else {
      span(class = "badge bg-secondary", style = "margin-right: 10px;",
        icon("question-circle"), " DIA-NN normalization: unknown")
    }
  })

  # Health assessment reactive
  assess_distribution_health <- reactive({
    req(values$raw_data, values$y_protein)
    pre_mat <- values$raw_data$E
    post_mat <- values$y_protein$E
    pre_medians <- apply(pre_mat, 2, median, na.rm = TRUE)
    post_medians <- apply(post_mat, 2, median, na.rm = TRUE)
    pre_cv <- sd(pre_medians) / abs(mean(pre_medians))
    post_cv <- sd(post_medians) / abs(mean(post_medians))
    pre_outliers <- which(abs(pre_medians - mean(pre_medians)) > 2 * sd(pre_medians))
    post_outliers <- which(abs(post_medians - mean(post_medians)) > 2 * sd(post_medians))
    list(
      pre_cv = pre_cv, post_cv = post_cv,
      pre_outlier_samples = names(pre_outliers),
      post_outlier_samples = names(post_outliers),
      status = if (pre_cv > 0.05) "warning" else if (length(pre_outliers) > 0) "caution" else "good"
    )
  })

  # Contextual guidance banners
  output$norm_diag_guidance <- renderUI({
    req(values$raw_data, values$y_protein)
    health <- assess_distribution_health()
    diann_status <- values$diann_norm_detected
    warnings <- list()

    # SCENARIO 1: DIA-NN normalization OFF + bad distributions
    if (diann_status == "off" && health$status == "warning") {
      warnings <- c(warnings, list(
        div(class = "alert alert-warning", role = "alert",
          icon("exclamation-triangle"),
          strong(" Unnormalized data detected. "),
          "Your DIA-NN output does not appear to be normalized (the ",
          tags$code("Precursor.Normalised"), " and ", tags$code("Precursor.Quantity"),
          " columns are identical). The sample distributions look uneven, which can lead ",
          "to unreliable differential expression results.",
          br(), br(),
          strong("What to do: "),
          "For most experiments, re-process your data in DIA-NN with ",
          tags$b("RT-dependent normalization"), " enabled (this is the default setting). ",
          "This corrects for differences in sample loading and LC-MS run variability.",
          br(), br(),
          em("Exception: "), "If you are analyzing AP-MS/Co-IP, fractionated samples, or ",
          "isotope labeling time-courses, unnormalized data may be appropriate \u2014 ",
          "but you should apply your own normalization before using DE-LIMP."
        )
      ))
    }

    # SCENARIO 2: DIA-NN normalization OFF but distributions look OK
    if (diann_status == "off" && health$status == "good") {
      warnings <- c(warnings, list(
        div(class = "alert alert-info", role = "alert",
          icon("info-circle"),
          strong(" DIA-NN normalization was off, "),
          "but your sample distributions look reasonably aligned. This can happen if your ",
          "samples had very consistent loading and LC-MS performance. Results may still be ",
          "valid, but consider whether normalization would improve your analysis."
        )
      ))
    }

    # SCENARIO 3: DIA-NN normalization ON but distributions still bad
    if (diann_status == "on" && health$status == "warning") {
      warnings <- c(warnings, list(
        div(class = "alert alert-danger", role = "alert",
          icon("times-circle"),
          strong(" Sample distributions are uneven despite normalization. "),
          "DIA-NN normalization was applied but the per-sample distributions still show ",
          "substantial differences. This could indicate:",
          tags$ul(
            tags$li("A failed or low-quality injection (check the QC Trends tab)"),
            tags$li("Very different sample types being compared (e.g., tissue vs plasma)"),
            tags$li("Severe batch effects that normalization couldn't fully correct")
          ),
          strong("What to do: "),
          "Check the QC Trends tab for outlier samples. Consider whether any samples ",
          "should be excluded. If batch effects are suspected, make sure you've assigned ",
          "batch information in the Assign Groups modal."
        )
      ))
    }

    # SCENARIO 4: Outlier sample(s) detected
    if (length(health$pre_outlier_samples) > 0) {
      outlier_names <- paste(health$pre_outlier_samples, collapse = ", ")
      warnings <- c(warnings, list(
        div(class = "alert alert-warning", role = "alert",
          icon("user-times"),
          strong(" Possible outlier sample(s): "),
          tags$code(outlier_names),
          br(),
          "These samples have median intensities substantially different from the rest. ",
          "This could indicate a failed injection, sample preparation issue, or ",
          "biological outlier. Check the QC Trends tab and MDS plot for confirmation. ",
          "If the sample is clearly problematic, consider re-running the analysis ",
          "without it."
        )
      ))
    }

    # SCENARIO 5: Everything looks good
    if (health$status == "good" && diann_status %in% c("on", "unknown") &&
        length(health$pre_outlier_samples) == 0) {
      warnings <- c(warnings, list(
        div(class = "alert alert-success", role = "alert",
          icon("check-circle"),
          strong(" Distributions look good. "),
          "Per-sample intensity distributions are well-aligned. ",
          "No outlier samples detected."
        )
      ))
    }

    do.call(tagList, warnings)
  })

  # Shared reactive for the diagnostic plot (regular + fullscreen)
  generate_norm_diagnostic_plot <- reactive({
    req(values$raw_data, values$y_protein, values$metadata)

    pre_mat <- values$raw_data$E   # precursor-level, log2, NAs present
    post_mat <- values$y_protein$E # protein-level, log2, no NAs
    meta <- values$metadata

    if (input$norm_diag_type == "boxplot") {
      # === BOX PLOT VIEW ===
      # Subsample precursor matrix if very large (performance)
      if (nrow(pre_mat) > 10000) {
        sample_idx <- sample(nrow(pre_mat), 10000)
        pre_mat_plot <- pre_mat[sample_idx, ]
      } else {
        pre_mat_plot <- pre_mat
      }

      pre_long <- as.data.frame(pre_mat_plot) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
        mutate(Stage = "Precursor Input\n(DIA-NN normalized)") %>%
        filter(!is.na(Log2Intensity))

      post_long <- as.data.frame(post_mat) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
        mutate(Stage = "Protein Output\n(DPC-Quant)")

      pre_long$Group <- meta$Group[match(pre_long$Sample, meta$File.Name)]
      post_long$Group <- meta$Group[match(post_long$Sample, meta$File.Name)]

      combined <- bind_rows(pre_long, post_long)
      combined$Stage <- factor(combined$Stage,
        levels = c("Precursor Input\n(DIA-NN normalized)", "Protein Output\n(DPC-Quant)"))

      sample_order <- meta %>% arrange(Group, File.Name) %>% pull(File.Name)
      combined$Sample <- factor(combined$Sample, levels = sample_order)
      combined$SampleID <- meta$ID[match(combined$Sample, meta$File.Name)]

      p <- ggplot(combined, aes(x = factor(SampleID), y = Log2Intensity, fill = Group)) +
        geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
        facet_wrap(~Stage, scales = "free_y", ncol = 2) +
        theme_minimal() +
        labs(
          title = "Pipeline Diagnostic: Precursor Input \u2192 Protein Output",
          subtitle = "Left: DIA-NN normalized precursors | Right: DPC-Quant protein estimates",
          x = "Sample ID", y = "Log2 Intensity"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

      ggplotly(p, tooltip = c("x", "y")) %>% layout(boxmode = "group")

    } else {
      # === DENSITY OVERLAY VIEW ===
      pre_long <- as.data.frame(pre_mat) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
        mutate(Stage = "Precursor Input (DIA-NN)") %>%
        filter(!is.na(Log2Intensity))

      post_long <- as.data.frame(post_mat) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
        mutate(Stage = "Protein Output (DPC-Quant)")

      pre_long$Group <- meta$Group[match(pre_long$Sample, meta$File.Name)]
      post_long$Group <- meta$Group[match(post_long$Sample, meta$File.Name)]

      combined <- bind_rows(pre_long, post_long)
      combined$Stage <- factor(combined$Stage,
        levels = c("Precursor Input (DIA-NN)", "Protein Output (DPC-Quant)"))

      p <- ggplot(combined, aes(x = Log2Intensity, color = Group, group = Sample)) +
        geom_density(alpha = 0.3, linewidth = 0.4) +
        facet_wrap(~Stage, ncol = 2) +
        theme_minimal() +
        labs(
          title = "Pipeline Diagnostic: Per-Sample Density Curves",
          subtitle = "Well-aligned input distributions should remain aligned after quantification",
          x = "Log2 Intensity", y = "Density"
        )

      ggplotly(p)
    }
  })

  output$norm_diagnostic_plot <- renderPlotly({ generate_norm_diagnostic_plot() })

  # Fullscreen modal for pipeline diagnostic
  observeEvent(input$fullscreen_norm_diag, {
    showModal(modalDialog(
      title = "Pipeline Diagnostic - Fullscreen View",
      plotlyOutput("norm_diagnostic_plot_fullscreen", height = "700px"),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  output$norm_diagnostic_plot_fullscreen <- renderPlotly({ generate_norm_diagnostic_plot() })

  # ============================================================================
  #      Fullscreen Modals for All Plot Panels
  # ============================================================================

  # --- 1. Signal Distribution (Data Overview) ---
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

  # --- 2. DPC Fit (QC Plots) ---
  observeEvent(input$fullscreen_dpc, {
    showModal(modalDialog(
      title = "DPC Fit - Fullscreen View",
      plotOutput("dpc_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$dpc_plot_fs <- renderPlot({ req(values$dpc_fit); limpa::plotDPC(values$dpc_fit) }, height = 700)

  # --- 3. MDS Plot (QC Plots) ---
  observeEvent(input$fullscreen_mds, {
    showModal(modalDialog(
      title = "MDS Plot - Fullscreen View",
      plotOutput("mds_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$mds_plot_fs <- renderPlot({
    req(values$y_protein, values$metadata)
    meta <- values$metadata[match(colnames(values$y_protein$E), values$metadata$File.Name), ]
    grps <- factor(meta$Group); cols <- rainbow(length(levels(grps)))
    par(xpd = TRUE)
    limpa::plotMDSUsingSEs(values$y_protein, pch = 16, main = "MDS Plot", col = cols[grps])
    legend(x = "right", inset = c(-0.15, 0), legend = levels(grps), col = cols, pch = 16, bty = "n")
  }, height = 700)

  # --- 4. Group QC Distribution Violin (QC Plots) ---
  observeEvent(input$fullscreen_qc_violin, {
    showModal(modalDialog(
      title = "Group QC Distribution - Fullscreen View",
      plotlyOutput("qc_group_violin_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$qc_group_violin_fs <- renderPlotly({
    req(values$qc_stats, values$metadata, input$qc_violin_metric)
    df <- left_join(values$qc_stats, values$metadata, by = c("Run" = "File.Name"))
    metric <- input$qc_violin_metric
    df$Tooltip <- paste0("<b>File:</b> ", df$Run, "<br><b>Val:</b> ", round(df[[metric]], 2))
    p <- ggplot(df, aes(x = Group, y = .data[[metric]], fill = Group)) +
      geom_violin(alpha = 0.5, trim = FALSE) +
      geom_jitter(aes(text = Tooltip), width = 0.2, size = 2, alpha = 0.8, color = "black") +
      theme_bw() + labs(title = paste("Distribution of", metric), x = "Group", y = metric) +
      theme(legend.position = "none")
    ggplotly(p, tooltip = "text")
  })

  # --- 5. P-value Distribution Histogram ---
  # Reactive to assess p-value distribution health
  assess_pvalue_health <- reactive({
    req(values$fit, input$contrast_selector_pvalue)

    # Get p-values
    de_results <- topTable(values$fit, coef = input$contrast_selector_pvalue, number = Inf)
    pvalues <- de_results$P.Value

    n_proteins <- length(pvalues)

    # Bin the p-values into 10 bins
    breaks <- seq(0, 1, by = 0.1)
    hist_counts <- hist(pvalues, breaks = breaks, plot = FALSE)$counts

    # Expected count per bin if uniform
    expected_per_bin <- n_proteins / length(hist_counts)

    # Calculate ratios for different regions
    low_pval_ratio <- sum(pvalues < 0.05) / (n_proteins * 0.05)  # Ratio vs expected 5%
    mid_pval_ratio <- sum(pvalues >= 0.3 & pvalues <= 0.7) / (n_proteins * 0.4)  # Ratio vs expected 40%

    # Detect patterns
    has_spike <- low_pval_ratio > 2  # Spike at zero if > 2x expected
    has_inflation <- mid_pval_ratio > 1.3  # Inflation if mid-range > 1.3x expected
    has_depletion <- low_pval_ratio < 0.5  # Depletion if < 0.5x expected

    # U-shaped: high at both ends
    first_bin_ratio <- hist_counts[1] / expected_per_bin
    last_bin_ratio <- hist_counts[length(hist_counts)] / expected_per_bin
    is_u_shaped <- (first_bin_ratio > 1.5 && last_bin_ratio > 1.5)

    # Completely uniform: no spike, no inflation
    is_uniform <- !has_spike && !has_inflation && (low_pval_ratio > 0.8 && low_pval_ratio < 1.2)

    # Determine overall status
    if (is_u_shaped) {
      status <- "u_shaped"
    } else if (has_inflation) {
      status <- "inflation"
    } else if (has_depletion && !has_spike) {
      status <- "low_power"
    } else if (is_uniform) {
      status <- "uniform"
    } else if (has_spike) {
      status <- "healthy"
    } else {
      status <- "unknown"
    }

    list(
      status = status,
      n_proteins = n_proteins,
      n_significant = sum(de_results$adj.P.Val < 0.05),
      low_pval_ratio = low_pval_ratio,
      mid_pval_ratio = mid_pval_ratio,
      has_spike = has_spike,
      has_inflation = has_inflation,
      has_depletion = has_depletion,
      is_u_shaped = is_u_shaped
    )
  })

  # Render contextual guidance banner
  output$pvalue_guidance <- renderUI({
    health <- assess_pvalue_health()

    if (health$status == "healthy") {
      # Green success banner
      div(class = "alert alert-success", role = "alert",
        icon("check-circle"),
        strong(" P-value distribution looks healthy. "),
        sprintf("Good spike near p=0 (%d significant proteins after FDR correction). ", health$n_significant),
        "This indicates genuine differential expression with proper statistical power."
      )

    } else if (health$status == "inflation") {
      # Yellow warning - p-value inflation
      div(class = "alert alert-warning", role = "alert",
        icon("exclamation-triangle"),
        strong(" Possible p-value inflation detected. "),
        "Too many intermediate p-values (0.3-0.7) relative to expectation. This may indicate:",
        tags$ul(
          tags$li("Unmodeled batch effects → Add batch covariate in Assign Groups tab"),
          tags$li("Variance heterogeneity → Check MDS Plot for outliers"),
          tags$li("Small sample size → Consider adding biological replicates")
        ),
        "If this pattern persists, consider checking the Normalization Diagnostic tab."
      )

    } else if (health$status == "low_power") {
      # Yellow warning - low power
      div(class = "alert alert-warning", role = "alert",
        icon("battery-quarter"),
        strong(" Low statistical power detected. "),
        sprintf("Fewer small p-values than expected (%d significant proteins). ", health$n_significant),
        "Possible causes:",
        tags$ul(
          tags$li("Small sample size → Increase biological replicates if possible"),
          tags$li("High biological variability → Check CV Distribution tab"),
          tags$li("Effect sizes too small to detect with current sample size"),
          tags$li("Over-conservative FDR correction → Consider less stringent threshold")
        )
      )

    } else if (health$status == "u_shaped") {
      # Red danger - U-shaped distribution
      div(class = "alert alert-danger", role = "alert",
        icon("times-circle"),
        strong(" Statistical model issue detected. "),
        "U-shaped p-value distribution (enrichment at both p~0 and p~1) suggests problems with the statistical model or data quality. ",
        "Recommended actions:",
        tags$ul(
          tags$li("Check Normalization Diagnostic - samples may not be properly normalized"),
          tags$li("Review MDS Plot for outlier samples or batch structure"),
          tags$li("Verify group assignments are correct"),
          tags$li("Consider whether this comparison is biologically appropriate")
        )
      )

    } else if (health$status == "uniform") {
      # Blue info - no signal
      div(class = "alert alert-info", role = "alert",
        icon("info-circle"),
        strong(" No differential expression signal detected. "),
        sprintf("P-values are uniformly distributed (only %d proteins pass FDR < 0.05). ", health$n_significant),
        "This could mean:",
        tags$ul(
          tags$li("Groups are truly similar (no biological difference)"),
          tags$li("Test lacks power to detect existing differences"),
          tags$li("Technical variation masks biological signal")
        ),
        "Consider checking QC plots to rule out technical issues."
      )

    } else {
      # Default blue info banner
      div(style = "background-color: #e7f3ff; padding: 12px; border-radius: 5px;",
        icon("info-circle"),
        strong(" P-value Diagnostic: "),
        "This histogram shows the distribution of raw p-values from your differential expression test. ",
        "A healthy analysis shows mostly uniform distribution (flat histogram) with enrichment near p=0 for true positives."
      )
    }
  })

  output$pvalue_histogram <- renderPlot({
    req(values$fit, input$contrast_selector_pvalue)

    # Get all p-values for the current contrast
    de_results <- topTable(values$fit, coef = input$contrast_selector_pvalue, number = Inf)
    pvalues <- de_results$P.Value

    # Calculate expected uniform distribution
    n_proteins <- length(pvalues)
    n_bins <- 30
    expected_per_bin <- n_proteins / n_bins

    # Create histogram data
    hist_data <- data.frame(PValue = pvalues)

    # Create the plot
    ggplot(hist_data, aes(x = PValue)) +
      geom_histogram(bins = n_bins, fill = "#4A90E2", color = "white", alpha = 0.7) +
      geom_hline(yintercept = expected_per_bin, linetype = "dashed", color = "red", size = 1) +
      annotate("text", x = 0.75, y = expected_per_bin * 1.1,
               label = "Expected under null (uniform)",
               color = "red", size = 3.5, fontface = "italic") +
      labs(
        title = paste0("P-value Distribution (", nrow(de_results), " proteins tested)"),
        subtitle = paste0("Comparison: ", input$contrast_selector_pvalue),
        x = "P-value",
        y = "Number of Proteins"
      ) +
      theme_bw(base_size = 14) +
      theme(
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray40", size = 11),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))
  })

  # Fullscreen modal for p-value histogram
  observeEvent(input$fullscreen_pvalue_hist, {
    req(values$fit, input$contrast_selector_pvalue)

    # Get all p-values for the current contrast
    de_results <- topTable(values$fit, coef = input$contrast_selector_pvalue, number = Inf)
    pvalues <- de_results$P.Value

    # Calculate expected uniform distribution
    n_proteins <- length(pvalues)
    n_bins <- 40  # More bins for fullscreen
    expected_per_bin <- n_proteins / n_bins

    # Create histogram data
    hist_data <- data.frame(PValue = pvalues)

    # Create enhanced plot for fullscreen
    p <- ggplot(hist_data, aes(x = PValue)) +
      geom_histogram(bins = n_bins, fill = "#4A90E2", color = "white", alpha = 0.7) +
      geom_hline(yintercept = expected_per_bin, linetype = "dashed", color = "red", size = 1.2) +
      annotate("text", x = 0.75, y = expected_per_bin * 1.15,
               label = "Expected uniform distribution",
               color = "red", size = 4, fontface = "italic") +
      labs(
        title = paste0("P-value Distribution: ", input$contrast_selector_pvalue),
        subtitle = paste0(nrow(de_results), " proteins tested | ",
                         sum(de_results$adj.P.Val < 0.05), " significant after FDR correction"),
        x = "Raw P-value",
        y = "Number of Proteins"
      ) +
      theme_bw(base_size = 16) +
      theme(
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(color = "gray40", size = 13),
        panel.grid.minor = element_blank()
      ) +
      scale_x_continuous(breaks = seq(0, 1, 0.1), limits = c(0, 1))

    showModal(modalDialog(
      title = "P-value Distribution - Fullscreen View",
      renderPlot({ p }, height = 700, width = 1000),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  # --- 6. Volcano Plot (DE Dashboard) ---
  observeEvent(input$fullscreen_volcano, {
    showModal(modalDialog(
      title = "Volcano Plot - Fullscreen View",
      plotlyOutput("volcano_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$volcano_plot_fs <- renderPlotly({
    df <- volcano_data()
    cols <- c("Not Sig" = "grey", "Significant" = "red")

    p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), text = paste("Protein:", Protein.Group), key = Protein.Group, color = Significance)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = cols) +

      # Threshold lines with color
      geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff),
                 linetype = "dashed", color = "#FFA500", size = 0.8) +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed", color = "#4169E1", size = 0.8) +

      theme_minimal() +
      labs(y = "-log10(P-Value)", title = paste0("Volcano Plot: ", input$contrast_selector))

    df_sel <- df %>% filter(Selected == "Yes")
    if (nrow(df_sel) > 0) {
      p <- p + geom_point(data = df_sel, aes(x = logFC, y = -log10(P.Value)),
                         shape = 21, size = 4, fill = NA, color = "blue", stroke = 2)
    }

    # Convert to plotly and add annotations using plotly's native system
    ggplotly(p, tooltip = "text") %>%
      layout(
        annotations = list(
          # Significance criteria box text (using plotly annotations for better positioning)
          list(x = 0.02, y = 0.98, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top",
               text = "<b>Significant if:</b>", showarrow = FALSE, font = list(size = 14)),
          list(x = 0.02, y = 0.93, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top",
               text = paste0("• FDR-adj. p < 0.05<br>• |log2FC| > ", round(input$logfc_cutoff, 2)),
               showarrow = FALSE, font = list(size = 12, color = "#555555"))
        ),
        shapes = list(
          # White background box for legend
          list(type = "rect", x0 = 0.01, x1 = 0.28, y0 = 0.88, y1 = 0.99,
               xref = "paper", yref = "paper", fillcolor = "white", opacity = 0.85,
               line = list(color = "#333333", width = 1))
        )
      )
  })

  # --- 6. Heatmap (DE Dashboard) ---
  observeEvent(input$fullscreen_heatmap, {
    showModal(modalDialog(
      title = "Heatmap - Fullscreen View",
      plotOutput("heatmap_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$heatmap_plot_fs <- renderPlot({
    req(values$fit, values$y_protein, input$contrast_selector)
    df_volc <- volcano_data(); prot_ids <- NULL
    if (!is.null(input$de_table_rows_selected)) {
      current_table_data <- df_volc
      if (!is.null(values$plot_selected_proteins)) current_table_data <- current_table_data %>% filter(Protein.Group %in% values$plot_selected_proteins)
      prot_ids <- current_table_data$Protein.Group[input$de_table_rows_selected]
    } else if (!is.null(values$plot_selected_proteins)) {
      prot_ids <- values$plot_selected_proteins; if (length(prot_ids) > 50) prot_ids <- head(prot_ids, 50)
    } else { top_prots <- topTable(values$fit, coef = input$contrast_selector, number = 20); prot_ids <- rownames(top_prots) }
    valid_ids <- intersect(prot_ids, rownames(values$y_protein$E)); if (length(valid_ids) == 0) return(NULL)
    mat <- values$y_protein$E[valid_ids, , drop = FALSE]; mat_z <- t(apply(mat, 1, cal_z_score))
    meta <- values$metadata[match(colnames(mat), values$metadata$File.Name), ]; groups <- factor(meta$Group)
    ha <- HeatmapAnnotation(Group = groups, col = list(Group = setNames(rainbow(length(levels(groups))), levels(groups))))
    Heatmap(mat_z, name = "Z-score", top_annotation = ha, cluster_rows = TRUE, cluster_columns = TRUE, show_column_names = FALSE)
  }, height = 700)

  # --- 7. GSEA Dot Plot ---
  observeEvent(input$fullscreen_gsea_dot, {
    showModal(modalDialog(
      title = "GSEA Dot Plot - Fullscreen View",
      plotOutput("gsea_dot_plot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$gsea_dot_plot_fs <- renderPlot({
    req(values$gsea_results)
    if (nrow(values$gsea_results) > 0) dotplot(values$gsea_results, showCategory = 20) + ggtitle("GSEA GO Biological Process")
    else plot(NULL, xlim = c(0, 1), ylim = c(0, 1), main = "No significant enrichment found.", xaxt = 'n', yaxt = 'n')
  }, height = 700)

  # --- 8. GSEA Enrichment Map ---
  observeEvent(input$fullscreen_gsea_emap, {
    showModal(modalDialog(
      title = "GSEA Enrichment Map - Fullscreen View",
      plotOutput("gsea_emapplot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$gsea_emapplot_fs <- renderPlot({
    req(values$gsea_results)
    if (nrow(values$gsea_results) > 0) emapplot(pairwise_termsim(values$gsea_results), showCategory = 20)
  }, height = 700)

  # --- 9. GSEA Ridgeplot ---
  observeEvent(input$fullscreen_gsea_ridge, {
    showModal(modalDialog(
      title = "GSEA Ridgeplot - Fullscreen View",
      plotOutput("gsea_ridgeplot_fs", height = "700px"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close")
    ))
  })
  output$gsea_ridgeplot_fs <- renderPlot({
    req(values$gsea_results)
    if (nrow(values$gsea_results) > 0) ridgeplot(values$gsea_results)
  }, height = 700)

  # ============================================================================

  output$de_table <- renderDT({
    req(values$fit, input$contrast_selector)

    # Build table data independently (not using volcano_data() to avoid reactive loops)
    df_raw <- topTable(values$fit, coef=input$contrast_selector, number=Inf) %>% as.data.frame()
    if (!"Protein.Group" %in% colnames(df_raw)) {
      df_full <- df_raw %>% rownames_to_column("Protein.Group")
    } else {
      df_full <- df_raw
    }

    # Add gene symbol and protein name
    org_db_name <- detect_organism_db(df_full$Protein.Group)
    df_full$Accession <- str_split_fixed(df_full$Protein.Group, "[; ]", 2)[,1]

    id_map <- tryCatch({
      if (!requireNamespace(org_db_name, quietly = TRUE)) { BiocManager::install(org_db_name, ask = FALSE) }
      library(org_db_name, character.only = TRUE)
      bitr(df_full$Accession, fromType = "UNIPROT", toType = c("SYMBOL", "GENENAME"), OrgDb = get(org_db_name))
    }, error = function(e) { NULL })

    if (!is.null(id_map)) {
      id_map <- id_map %>% distinct(UNIPROT, .keep_all = TRUE)
      df_full <- df_full %>% left_join(id_map, by = c("Accession" = "UNIPROT")) %>%
        mutate(Gene = ifelse(is.na(SYMBOL), Accession, SYMBOL),
               Protein.Name = ifelse(is.na(GENENAME), Protein.Group, GENENAME))
    } else {
      df_full$Gene <- df_full$Accession
      df_full$Protein.Name <- df_full$Protein.Group
    }

    df_full$Significance <- "Not Sig"
    df_full$Significance[df_full$adj.P.Val < 0.05 & abs(df_full$logFC) > input$logfc_cutoff] <- "Significant"

    # Prepare display (don't filter, keep all rows for multi-select)
    df_display <- df_full %>% mutate(across(where(is.numeric), function(x) round(x,4))) %>%
      mutate(Protein.Name_Link = ifelse(!is.na(Accession) & str_detect(Accession, "^[A-Z0-9]{6,}$"),
                                       paste0("<a href='https://www.uniprot.org/uniprotkb/", Accession,
                                              "/entry' target='_blank' onclick='window.open(this.href, \"_blank\"); return false;'>",
                                              Protein.Name, "</a>"),
                                       Protein.Name)) %>%
      dplyr::select(Gene, `Protein Name` = Protein.Name_Link, logFC, P.Value, adj.P.Val, Significance)

    datatable(df_display, selection = "multiple", options = list(pageLength = 10, scrollX = TRUE), escape = FALSE, rownames = FALSE)
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

  # CV Histogram: Distribution of CV values by group
  output$cv_histogram <- renderPlot({
    req(values$fit, values$y_protein, input$contrast_selector, values$metadata)

    # Get significant proteins
    df_res_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>%
      as.data.frame() %>%
      filter(adj.P.Val < 0.05)
    if (!"Protein.Group" %in% colnames(df_res_raw)) {
      df_res <- df_res_raw %>% rownames_to_column("Protein.Group")
    } else {
      df_res <- df_res_raw
    }
    if(nrow(df_res) == 0) {
      # Show message if no significant proteins
      plot.new()
      text(0.5, 0.5, "No significant proteins found.\nAdjust significance threshold or check data.",
           cex = 1.2, col = "gray50")
      return()
    }

    # Calculate CV for each group
    protein_ids_for_cv <- df_res$Protein.Group
    raw_exprs <- values$y_protein$E[protein_ids_for_cv, , drop = FALSE]
    linear_exprs <- 2^raw_exprs
    cv_list <- list()

    for(g in unique(values$metadata$Group)) {
      if (g == "") next
      files_in_group <- values$metadata$File.Name[values$metadata$Group == g]
      group_cols <- intersect(colnames(linear_exprs), files_in_group)

      if (length(group_cols) > 1) {
        group_data <- linear_exprs[, group_cols, drop = FALSE]
        cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) {
          (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
        })
      } else {
        cv_list[[paste0("CV_", g)]] <- NA
      }
    }

    # Convert to long format for ggplot
    cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
    cv_long <- cv_df %>%
      pivot_longer(cols = starts_with("CV_"),
                   names_to = "Group",
                   values_to = "CV") %>%
      mutate(Group = gsub("CV_", "", Group)) %>%
      filter(!is.na(CV))

    # Calculate averages for each group
    cv_averages <- cv_long %>%
      group_by(Group) %>%
      summarise(Avg_CV = mean(CV, na.rm = TRUE), .groups = 'drop')

    # Create histogram with facets
    ggplot(cv_long, aes(x = CV)) +
      geom_histogram(aes(fill = Group), bins = 30, alpha = 0.7, color = "white") +
      geom_vline(data = cv_averages, aes(xintercept = Avg_CV, color = Group),
                 linetype = "dashed", size = 1.2) +
      geom_text(data = cv_averages,
                aes(x = Avg_CV, y = Inf, label = paste0("Avg: ", round(Avg_CV, 1), "%")),
                vjust = 1.5, hjust = -0.1, size = 3.5, fontface = "bold") +
      facet_wrap(~ Group, ncol = 2, scales = "free_y") +
      labs(title = paste0("CV Distribution by Group (", nrow(df_res), " significant proteins)"),
           subtitle = "Dashed line shows average CV for each group",
           x = "Coefficient of Variation (%)",
           y = "Number of Proteins") +
      theme_bw(base_size = 14) +
      theme(
        legend.position = "none",
        strip.background = element_rect(fill = "#667eea", color = NA),
        strip.text = element_text(color = "white", face = "bold", size = 12),
        plot.title = element_text(face = "bold", size = 16),
        plot.subtitle = element_text(color = "gray40", size = 11),
        panel.grid.minor = element_blank()
      )
  })

  # Fullscreen modal for CV histogram
  observeEvent(input$fullscreen_cv_hist, {
    req(values$fit, values$y_protein, input$contrast_selector, values$metadata)

    # Get significant proteins
    df_res_raw <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>%
      as.data.frame() %>%
      filter(adj.P.Val < 0.05)
    if (!"Protein.Group" %in% colnames(df_res_raw)) {
      df_res <- df_res_raw %>% rownames_to_column("Protein.Group")
    } else {
      df_res <- df_res_raw
    }
    if(nrow(df_res) == 0) {
      showNotification("No significant proteins found for CV analysis.", type = "warning")
      return()
    }

    # Calculate CV for each group
    protein_ids_for_cv <- df_res$Protein.Group
    raw_exprs <- values$y_protein$E[protein_ids_for_cv, , drop = FALSE]
    linear_exprs <- 2^raw_exprs
    cv_list <- list()

    for(g in unique(values$metadata$Group)) {
      if (g == "") next
      files_in_group <- values$metadata$File.Name[values$metadata$Group == g]
      group_cols <- intersect(colnames(linear_exprs), files_in_group)

      if (length(group_cols) > 1) {
        group_data <- linear_exprs[, group_cols, drop = FALSE]
        cv_list[[paste0("CV_", g)]] <- apply(group_data, 1, function(x) {
          (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)) * 100
        })
      } else {
        cv_list[[paste0("CV_", g)]] <- NA
      }
    }

    # Convert to long format for ggplot
    cv_df <- as.data.frame(cv_list) %>% rownames_to_column("Protein.Group")
    cv_long <- cv_df %>%
      pivot_longer(cols = starts_with("CV_"),
                   names_to = "Group",
                   values_to = "CV") %>%
      mutate(Group = gsub("CV_", "", Group)) %>%
      filter(!is.na(CV))

    # Calculate averages for each group
    cv_averages <- cv_long %>%
      group_by(Group) %>%
      summarise(Avg_CV = mean(CV, na.rm = TRUE), .groups = 'drop')

    # Create fullscreen histogram
    p <- ggplot(cv_long, aes(x = CV)) +
      geom_histogram(aes(fill = Group), bins = 40, alpha = 0.7, color = "white") +
      geom_vline(data = cv_averages, aes(xintercept = Avg_CV, color = Group),
                 linetype = "dashed", size = 1.5) +
      geom_text(data = cv_averages,
                aes(x = Avg_CV, y = Inf, label = paste0("Avg: ", round(Avg_CV, 1), "%")),
                vjust = 1.5, hjust = -0.1, size = 4, fontface = "bold") +
      facet_wrap(~ Group, ncol = 2, scales = "free_y") +
      labs(title = paste0("CV Distribution by Group (", nrow(df_res), " significant proteins)"),
           subtitle = "Dashed line shows average CV for each group. Lower CV = more stable/reproducible biomarker",
           x = "Coefficient of Variation (%)",
           y = "Number of Proteins") +
      theme_bw(base_size = 16) +
      theme(
        legend.position = "none",
        strip.background = element_rect(fill = "#667eea", color = NA),
        strip.text = element_text(color = "white", face = "bold", size = 14),
        plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(color = "gray40", size = 12),
        panel.grid.minor = element_blank()
      )

    showModal(modalDialog(
      title = "CV Distribution - Fullscreen View",
      renderPlot({ p }, height = 700, width = 1000),
      size = "xl",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })
  
  output$reproducible_code <- renderText({
    req(values$repro_log)
    log_content <- paste(values$repro_log, collapse = "\n")
    session_info_text <- paste(capture.output(sessionInfo()), collapse = "\n")

    # Add helpful footer
    footer <- c(
      "",
      "# ==============================================================================",
      "# How to Use This Log:",
      "# 1. Replace 'path/to/your/report.parquet' with your actual file path",
      "# 2. Ensure all required packages are installed (see Session Info below)",
      "# 3. Run sections sequentially from top to bottom",
      "# 4. Each section is timestamped showing when you performed that action",
      "# ==============================================================================",
      "",
      "# --- Session Info (Package Versions) ---",
      session_info_text
    )

    paste(c(log_content, footer), collapse = "\n")
  })

  # Download handler for reproducibility log
  output$download_repro_log <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_reproducibility_log_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R")
    },
    content = function(file) {
      req(values$repro_log)
      log_content <- paste(values$repro_log, collapse = "\n")
      session_info_text <- paste(capture.output(sessionInfo()), collapse = "\n")

      footer <- c(
        "",
        "# ==============================================================================",
        "# How to Use This Log:",
        "# 1. Replace 'path/to/your/report.parquet' with your actual file path",
        "# 2. Ensure all required packages are installed (see Session Info below)",
        "# 3. Run sections sequentially from top to bottom",
        "# 4. Each section is timestamped showing when you performed that action",
        "# ==============================================================================",
        "",
        "# --- Session Info (Package Versions) ---",
        session_info_text
      )

      writeLines(paste(c(log_content, footer), collapse = "\n"), file)
    }
  )

  # ============================================================================
  #      Save / Load Session (RDS)
  # ============================================================================

  output$save_session <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_session_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      session_data <- list(
        raw_data   = values$raw_data,
        metadata   = values$metadata,
        fit        = values$fit,
        y_protein  = values$y_protein,
        dpc_fit    = values$dpc_fit,
        design     = values$design,
        qc_stats   = values$qc_stats,
        gsea_results = values$gsea_results,
        repro_log  = values$repro_log,
        color_plot_by_de = values$color_plot_by_de,
        # Store UI state so it can be restored
        contrast   = input$contrast_selector,
        logfc_cutoff = input$logfc_cutoff,
        q_cutoff   = input$q_cutoff,
        # Covariate settings
        cov1_name  = values$cov1_name,
        cov2_name  = values$cov2_name,
        # Save timestamp & version
        saved_at   = Sys.time(),
        app_version = "DE-LIMP v1.2"
      )
      saveRDS(session_data, file)
      showNotification("Session saved successfully!", type = "message")
    }
  )

  observeEvent(input$load_session_btn, {
    showModal(modalDialog(
      title = "Load Saved Session",
      fileInput("session_file", "Choose .rds session file", accept = ".rds"),
      div(style = "background-color: #fff3cd; padding: 10px; border-radius: 5px;",
        icon("exclamation-triangle"),
        " Loading a session will replace all current data and results."
      ),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_load_session", "Load Session", class = "btn-primary", icon = icon("upload"))
      )
    ))
  })

  observeEvent(input$confirm_load_session, {
    req(input$session_file)
    tryCatch({
      session_data <- readRDS(input$session_file$datapath)

      # Validate that this is a DE-LIMP session file
      required_fields <- c("raw_data", "metadata", "fit")
      if (!all(required_fields %in% names(session_data))) {
        showNotification("Invalid session file: missing required data fields.", type = "error")
        return()
      }

      # Restore reactive values
      values$raw_data   <- session_data$raw_data
      values$metadata   <- session_data$metadata
      values$fit        <- session_data$fit
      values$y_protein  <- session_data$y_protein
      values$dpc_fit    <- session_data$dpc_fit
      values$design     <- session_data$design
      values$qc_stats   <- session_data$qc_stats
      values$gsea_results <- session_data$gsea_results
      values$color_plot_by_de <- session_data$color_plot_by_de %||% FALSE
      values$cov1_name  <- session_data$cov1_name %||% "Covariate1"
      values$cov2_name  <- session_data$cov2_name %||% "Covariate2"

      # Restore repro log and append load event
      values$repro_log  <- session_data$repro_log %||% values$repro_log
      add_to_log("Session Loaded", c(
        sprintf("# Loaded session saved at: %s", session_data$saved_at),
        sprintf("# App version: %s", session_data$app_version %||% "unknown")
      ))

      # Restore UI state: update contrast choices from the fit object
      if (!is.null(values$fit)) {
        contrast_names <- colnames(values$fit$contrasts)
        selected_contrast <- session_data$contrast %||% contrast_names[1]

        # Update all four comparison selectors
        updateSelectInput(session, "contrast_selector", choices = contrast_names, selected = selected_contrast)
        updateSelectInput(session, "contrast_selector_signal", choices = contrast_names, selected = selected_contrast)
        updateSelectInput(session, "contrast_selector_grid", choices = contrast_names, selected = selected_contrast)
        updateSelectInput(session, "contrast_selector_pvalue", choices = contrast_names, selected = selected_contrast)
      }
      if (!is.null(session_data$logfc_cutoff)) {
        updateSliderInput(session, "logfc_cutoff", value = session_data$logfc_cutoff)
      }
      if (!is.null(session_data$q_cutoff)) {
        updateNumericInput(session, "q_cutoff", value = session_data$q_cutoff)
      }

      values$status <- "Session loaded successfully"
      removeModal()
      showNotification(
        paste0("Session loaded! (saved ", format(session_data$saved_at, "%Y-%m-%d %H:%M"), ")"),
        type = "message", duration = 5
      )
    }, error = function(e) {
      showNotification(paste("Error loading session:", e$message), type = "error")
    })
  })

  output$methodology_text <- renderText({
    req(values$fit)

    methodology <- paste(
      "METHODOLOGY\n",
      "===========\n",
      "Data Processing and Statistical Analysis Pipeline\n",
      "---------------------------------------------------\n\n",

      "1. DATA INPUT\n",
      "Raw DIA-NN output files (parquet format) containing precursor-level quantification were\n",
      "imported using the limpa package. Input data includes precursor intensities, protein\n",
      "grouping information, and quality metrics (Q-values).\n\n",

      "2. NORMALIZATION\n",
      "Data Point Correspondence - Cyclic Normalization (DPC-CN) was applied using the dpcCN()\n",
      "function. This method normalizes signal intensities across runs by identifying invariant\n",
      "data points and applying cyclic loess normalization to correct for systematic technical\n",
      "variation while preserving biological differences. DPC-CN is specifically designed for\n",
      "DIA-NN data and performs robust normalization without requiring reference proteins or\n",
      "assuming equal protein abundances across samples.\n\n",

      "3. PROTEIN QUANTIFICATION\n",
      "Normalized precursor-level data were aggregated to protein-level quantification using\n",
      "dpcQuant(). This function employs a modified version of the maxLFQ algorithm, which:\n",
      "  • Identifies peptides/precursors unique to each protein group\n",
      "  • Uses pairwise ratios to estimate relative protein abundance\n",
      "  • Maximizes information from all available peptides while handling missing values\n",
      "  • Produces log2-transformed protein intensities for downstream analysis\n\n",

      "4. DIFFERENTIAL EXPRESSION ANALYSIS\n",
      "Statistical analysis was performed using the limma framework (Linear Models for\n",
      "Microarray Data), adapted for proteomics data through the dpcDE() function.\n",
      "The analysis workflow includes:\n\n",

      "  a) Linear Model Fitting:\n",
      "     A linear model was fit to the log2-transformed protein intensities with experimental\n",
      "     groups as factors. This model accounts for the mean-variance relationship in the data.\n\n",

      "  b) Empirical Bayes Moderation:\n",
      "     Variance estimates were moderated across proteins using empirical Bayes methods\n",
      "     (eBayes()). This 'borrows information' across proteins to stabilize variance estimates,\n",
      "     particularly beneficial for experiments with limited replicates.\n\n",

      "  c) Contrast Analysis:\n",
      "     Pairwise comparisons between experimental groups were performed using contrast matrices.\n",
      "     Each contrast produces:\n",
      "       • Log2 fold change (logFC): Effect size of differential expression\n",
      "       • Average expression (AveExpr): Mean log2 intensity across all samples\n",
      "       • t-statistic: Test statistic for differential expression\n",
      "       • P-value: Statistical significance of the change\n",
      "       • Adjusted P-value (adj.P.Val): FDR-corrected p-value using Benjamini-Hochberg method\n\n",

      "5. MULTIPLE TESTING CORRECTION\n",
      "False Discovery Rate (FDR) control was applied using the Benjamini-Hochberg procedure.\n",
      "This method controls the expected proportion of false positives among rejected hypotheses,\n",
      "providing adjusted p-values (adj.P.Val) that account for testing thousands of proteins\n",
      "simultaneously. Proteins with adj.P.Val < 0.05 are considered statistically significant\n",
      "at 5% FDR.\n\n",

      "6. GENE SET ENRICHMENT ANALYSIS (Optional)\n",
      "When performed, over-representation analysis (ORA) was conducted using clusterProfiler.\n",
      "Significant proteins (adj.P.Val < 0.05) were mapped to Gene Ontology (GO) terms, and\n",
      "enrichment was tested using hypergeometric distribution with FDR correction.\n\n\n",

      "SOFTWARE AND PACKAGES\n",
      "---------------------\n",
      "Primary analysis: limpa R package (Bioconductor 3.22+)\n",
      "Statistical framework: limma (Linear Models for Microarray and RNA-Seq Data)\n",
      "Data manipulation: dplyr, tidyr\n",
      "Visualization: ggplot2, ComplexHeatmap, plotly\n",
      "Enrichment: clusterProfiler, enrichplot\n",
      sprintf("R version: %s\n\n\n", R.version.string),

      "REFERENCES\n",
      "----------\n",
      "• limpa package: Bioconductor (https://bioconductor.org/packages/limpa/)\n",
      "• DPC normalization: Designed for DIA-NN proteomics data\n",
      "• limma: Ritchie ME, et al. (2015) Nucleic Acids Research 43(7):e47\n",
      "• Empirical Bayes: Smyth GK (2004) Statistical Applications in Genetics and\n",
      "  Molecular Biology 3:Article3\n",
      "• FDR control: Benjamini Y, Hochberg Y (1995) Journal of the Royal Statistical\n",
      "  Society 57(1):289-300\n",
      "• clusterProfiler: Yu G, et al. (2012) OMICS 16(5):284-287\n\n\n",

      "CITATION\n",
      "--------\n",
      "If you use this analysis in your research, please cite:\n",
      "• The limpa package (Bioconductor)\n",
      "• The limma package: Ritchie ME, et al. (2015) Nucleic Acids Research\n",
      "• DIA-NN: Demichev V, et al. (2020) Nature Methods 17:41-44",

      sep = ""
    )

    methodology
  })
  
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

        # Log GSEA analysis
        gsea_code <- c(
          sprintf("# Contrast: %s", input$contrast_selector),
          sprintf("# Organism DB: %s", org_db_name),
          "library(clusterProfiler)",
          sprintf("library(%s)", org_db_name),
          "",
          sprintf("de_results <- topTable(fit, coef='%s', number=Inf)", input$contrast_selector),
          "protein_ids <- str_extract(rownames(de_results), '[A-Z0-9]+')",
          sprintf("id_map <- bitr(protein_ids, fromType='UNIPROT', toType='ENTREZID', OrgDb=%s)", org_db_name),
          "de_results_mapped <- de_results[rownames(de_results) %in% id_map$UNIPROT, ]",
          "de_results_mapped$ENTREZID <- id_map$ENTREZID[match(str_extract(rownames(de_results_mapped), '[A-Z0-9]+'), id_map$UNIPROT)]",
          "gene_list <- de_results_mapped$logFC",
          "names(gene_list) <- de_results_mapped$ENTREZID",
          "gene_list <- sort(gene_list, decreasing=TRUE)",
          "gene_list <- gene_list[!duplicated(names(gene_list))]",
          "",
          sprintf("gsea_res <- gseGO(geneList=gene_list, OrgDb=%s, keyType='ENTREZID',", org_db_name),
          "                  ont='BP', minGSSize=10, maxGSSize=500, pvalueCutoff=0.05)"
        )
        add_to_log("GSEA Analysis", gsea_code)

      }, error = function(e) { output$gsea_status <- renderText(paste("GSEA Error:", e$message)) })
    })
  })
  
  output$gsea_dot_plot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) dotplot(values$gsea_results, showCategory = 20) + ggtitle("GSEA GO Biological Process") else plot(NULL, xlim=c(0,1), ylim=c(0,1), main="No significant enrichment found.", xaxt='n', yaxt='n') }) # Height controlled by UI (calc(100vh - 340px))
  output$gsea_emapplot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) emapplot(pairwise_termsim(values$gsea_results), showCategory=20) }) # Height controlled by UI (calc(100vh - 340px))
  output$gsea_ridgeplot <- renderPlot({ req(values$gsea_results); if(nrow(values$gsea_results) > 0) ridgeplot(values$gsea_results) }) # Height controlled by UI (calc(100vh - 340px))
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
    content = function(file) {
      req(values$fit, values$y_protein)
      de_stats <- topTable(values$fit, coef=input$contrast_selector, number=Inf) %>% rownames_to_column("Protein.Group")
      exprs_data <- as.data.frame(values$y_protein$E) %>% rownames_to_column("Protein.Group")
      full_data <- left_join(de_stats, exprs_data, by="Protein.Group")
      write.csv(full_data, file, row.names=FALSE)

      # Log export
      add_to_log("Export Results to CSV", c(
        sprintf("# Exported: %s", basename(file)),
        sprintf("de_stats <- topTable(fit, coef='%s', number=Inf) %%>%% rownames_to_column('Protein.Group')", input$contrast_selector),
        "exprs_data <- as.data.frame(y_protein$E) %>% rownames_to_column('Protein.Group')",
        "full_data <- left_join(de_stats, exprs_data, by='Protein.Group')",
        sprintf("write.csv(full_data, 'Limpa_Results_%s.csv', row.names=FALSE)", make.names(input$contrast_selector))
      ))
    }
  )

  # Export group assignment template
  output$export_template <- downloadHandler(
    filename = function() {
      paste0("DE-LIMP_group_template_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      # Get current table data (including any edits)
      template_data <- if (!is.null(input$hot_metadata)) {
        hot_to_r(input$hot_metadata)
      } else {
        values$metadata
      }

      # Export with current custom covariate names
      write.csv(template_data, file, row.names = FALSE)
      showNotification("Template exported successfully!", type = "message", duration = 3)
    }
  )

  # Import group assignment template
  observeEvent(input$import_template, {
    showModal(modalDialog(
      title = "Import Group Assignment Template",
      fileInput("template_file", "Choose CSV File",
                accept = c("text/csv", "text/comma-separated-values", ".csv")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_import", "Import", class = "btn-primary")
      )
    ))
  })

  observeEvent(input$confirm_import, {
    req(input$template_file)

    tryCatch({
      imported_data <- read.csv(input$template_file$datapath, stringsAsFactors = FALSE)

      # Validate columns
      required_cols <- c("ID", "File.Name", "Group", "Batch", "Covariate1", "Covariate2")
      if (!all(required_cols %in% colnames(imported_data))) {
        showNotification("Error: Template must have columns: ID, File.Name, Group, Batch, Covariate1, Covariate2",
                        type = "error", duration = 10)
        return()
      }

      # Validate File.Name matches
      if (!all(imported_data$File.Name %in% values$metadata$File.Name)) {
        showNotification("Warning: Some file names in template don't match current data. Using matching rows only.",
                        type = "warning", duration = 8)
      }

      # Update metadata with imported data (match by File.Name)
      values$metadata <- imported_data

      showNotification("Template imported successfully!", type = "message", duration = 3)
      removeModal()  # Close import dialog

      # Navigate to Assign Groups sub-tab to show imported data
      nav_select("main_tabs", "Data Overview")
      nav_select("data_overview_tabs", "Assign Groups & Run")

    }, error = function(e) {
      showNotification(paste("Error importing template:", e$message), type = "error", duration = 10)
    })
  })

  output$volcano_plot_interactive <- renderPlotly({
    df <- volcano_data()
    cols <- c("Not Sig" = "grey", "Significant" = "red")

    # Use non-adjusted P.Value for y-axis, but color by adjusted p-value significance
    p <- ggplot(df, aes(x = logFC, y = -log10(P.Value), text = paste("Protein:", Protein.Group), key = Protein.Group, color = Significance)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = cols) +

      # Threshold lines with color
      geom_vline(xintercept = c(-input$logfc_cutoff, input$logfc_cutoff),
                 linetype = "dashed", color = "#FFA500", size = 0.7) +
      geom_hline(yintercept = -log10(0.05),
                 linetype = "dashed", color = "#4169E1", size = 0.7) +

      theme_minimal() +
      labs(y = "-log10(P-Value)", title = paste0("Volcano Plot: ", input$contrast_selector))

    df_sel <- df %>% filter(Selected == "Yes")
    if (nrow(df_sel) > 0) {
      p <- p + geom_point(data = df_sel, aes(x = logFC, y = -log10(P.Value)),
                         shape = 21, size = 4, fill = NA, color = "blue", stroke = 2)
    }

    # Convert to plotly and add annotations using plotly's native system
    ggplotly(p, tooltip = "text", source = "volcano_source") %>%
      layout(
        dragmode = "select",
        annotations = list(
          # Significance criteria box text (using plotly annotations for better positioning)
          list(x = 0.02, y = 0.98, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top",
               text = "<b>Significant if:</b>", showarrow = FALSE, font = list(size = 12)),
          list(x = 0.02, y = 0.93, xref = "paper", yref = "paper", xanchor = "left", yanchor = "top",
               text = paste0("• FDR-adj. p < 0.05<br>• |log2FC| > ", round(input$logfc_cutoff, 2)),
               showarrow = FALSE, font = list(size = 11, color = "#555555"))
        ),
        shapes = list(
          # White background box for legend
          list(type = "rect", x0 = 0.01, x1 = 0.28, y0 = 0.88, y1 = 0.99,
               xref = "paper", yref = "paper", fillcolor = "white", opacity = 0.85,
               line = list(color = "#333333", width = 1))
        )
      )
  })
  
  observeEvent(event_data("plotly_selected", source = "volcano_source"), { select_data <- event_data("plotly_selected", source = "volcano_source"); if (!is.null(select_data)) values$plot_selected_proteins <- select_data$key })
  observeEvent(event_data("plotly_click", source = "volcano_source"), { click_data <- event_data("plotly_click", source = "volcano_source"); if (!is.null(click_data)) values$plot_selected_proteins <- click_data$key })
  observeEvent(input$clear_plot_selection, { values$plot_selected_proteins <- NULL })

  # Sync table row selection with plot selection system
  observeEvent(input$de_table_rows_selected, {
    req(input$de_table_rows_selected, length(input$de_table_rows_selected) > 0)
    df_full <- volcano_data()

    # Table always shows full data (not filtered), so row indices map directly
    selected_proteins <- df_full$Protein.Group[input$de_table_rows_selected]

    if (length(selected_proteins) > 0) {
      values$plot_selected_proteins <- selected_proteins
    }
  })
  
  # --- ADDED: Violin Plot Button Logic ---
  observeEvent(input$show_violin, {
    if (is.null(values$plot_selected_proteins) || length(values$plot_selected_proteins) == 0) {
      showNotification("⚠️ Please select a protein in the Volcano Plot or Table first!", type = "warning")
      return()
    }
    # Store all selected proteins for plotting
    values$temp_violin_target <- values$plot_selected_proteins

    n_proteins <- length(values$plot_selected_proteins)
    title_text <- if (n_proteins == 1) {
      paste("Expression Profile:", values$plot_selected_proteins[1])
    } else {
      paste("Expression Profiles for", n_proteins, "Selected Proteins")
    }

    showModal(modalDialog(
      title = title_text,
      size = "xl",
      plotOutput("violin_plot_de_popup", height = paste0(max(400, 200 * ceiling(n_proteins / 2)), "px")),
      footer = modalButton("Close"),
      easyClose = TRUE
    ))
  })

  output$violin_plot_de_popup <- renderPlot({
    req(values$y_protein, values$temp_violin_target, values$metadata)
    prot_ids <- values$temp_violin_target

    # Get expression data for all selected proteins
    exprs_mat <- values$y_protein$E[prot_ids, , drop=FALSE]
    long_df <- as.data.frame(exprs_mat) %>%
      rownames_to_column("Protein") %>%
      pivot_longer(-Protein, names_to = "File.Name", values_to = "LogIntensity")
    long_df <- left_join(long_df, values$metadata, by="File.Name")

    # Create violin plots with faceting for multiple proteins
    ggplot(long_df, aes(x = Group, y = LogIntensity, fill = Group)) +
      geom_violin(alpha = 0.5, trim = FALSE) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      facet_wrap(~Protein, scales = "free_y", ncol = 2) +
      theme_bw() +
      labs(y = "Log2 Intensity", x = "Group") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_rect(fill = "lightblue"),
        strip.text = element_text(face = "bold")
      )
  })

  # ============================================================================
  #      XIC Viewer: Server Logic
  # ============================================================================

  # --- XIC Directory Loading ---
  observeEvent(input$xic_load_dir, {
    req(input$xic_dir_input)
    xic_path <- trimws(input$xic_dir_input)

    # Smart path resolution: if user enters a .parquet file path, derive the _xic sibling
    if (grepl("\\.parquet$", xic_path, ignore.case = TRUE) && file.exists(xic_path)) {
      report_stem <- tools::file_path_sans_ext(basename(xic_path))
      derived_dir <- file.path(dirname(xic_path), paste0(report_stem, "_xic"))
      if (dir.exists(derived_dir)) {
        xic_path <- derived_dir
        updateTextInput(session, "xic_dir_input", value = xic_path)
        showNotification(paste0("Derived XIC directory: ", basename(xic_path)),
          type = "message", duration = 3)
      }
    }

    # If path doesn't end with _xic, try appending _xic (user may have pasted the parent dir)
    if (!dir.exists(xic_path) && !grepl("_xic/?$", xic_path)) {
      candidate <- paste0(xic_path, "_xic")
      if (dir.exists(candidate)) {
        xic_path <- candidate
        updateTextInput(session, "xic_dir_input", value = xic_path)
      }
    }

    if (!dir.exists(xic_path)) {
      showNotification("Directory not found. Check the path and try again.",
        type = "error", duration = 5)
      values$xic_available <- FALSE
      return()
    }

    xic_files <- list.files(xic_path, pattern = "\\.xic\\.parquet$",
                            full.names = TRUE, recursive = TRUE)

    if (length(xic_files) > 0) {
      values$xic_dir <- xic_path
      values$xic_available <- TRUE

      # Detect XIC format version (v1 = DIA-NN 1.x wide, v2 = DIA-NN 2.x long)
      values$xic_format <- detect_xic_format(xic_path)
      message(paste("Detected XIC format:", values$xic_format))

      # Detect mobilogram files and check if they contain non-zero data
      mob_files <- list.files(xic_path, pattern = "\\.mobilogram\\.parquet$",
                              full.names = TRUE, recursive = TRUE)
      if (length(mob_files) > 0) {
        # Sample one file to check for non-zero data (IM instruments only)
        mob_has_data <- tryCatch({
          mob_sample <- arrow::read_parquet(mob_files[1], as_data_frame = TRUE)
          any(mob_sample$data != 0, na.rm = TRUE)
        }, error = function(e) FALSE)
        values$mobilogram_available <- mob_has_data
        if (mob_has_data) values$mobilogram_dir <- xic_path
        message(paste("Mobilogram files:", length(mob_files),
                      "| Has IM data:", mob_has_data))
      } else {
        values$mobilogram_available <- FALSE
      }

      # Build precursor map directly from loaded data (no file re-reading needed)
      if (!is.null(values$raw_data)) {
        tryCatch({
          precursor_ids <- rownames(values$raw_data$E)
          protein_groups <- values$raw_data$genes$Protein.Group
          values$xic_report_map <- data.frame(
            Precursor.Id = precursor_ids,
            Protein.Group = protein_groups,
            stringsAsFactors = FALSE
          )
          message(paste("XIC precursor map built from loaded data:",
                        length(unique(protein_groups)), "proteins,",
                        length(unique(precursor_ids)), "precursors"))
        }, error = function(e) {
          values$xic_report_map <- NULL
          showNotification(paste("Error building precursor map:", e$message),
            type = "error", duration = 8)
        })
      } else {
        values$xic_report_map <- NULL
        showNotification(
          "Please load your DIA-NN report data first (Step 1), then load XICs.",
          type = "warning", duration = 8)
      }

      status_msg <- paste("Found", length(xic_files), "XIC files")
      if (values$mobilogram_available) {
        status_msg <- paste0(status_msg, " + ion mobility data")
      }
      showNotification(paste0(status_msg, ". Select a protein to view chromatograms."),
        type = "message", duration = 5)
    } else {
      values$xic_available <- FALSE
      showNotification(
        "No .xic.parquet files found in selected directory.",
        type = "warning", duration = 5)
    }
  })

  # --- XIC Status Badge ---
  output$xic_status_badge <- renderUI({
    if (values$xic_available) {
      xic_files <- list.files(values$xic_dir, pattern = "\\.xic\\.parquet$")
      badge_text <- paste(length(xic_files), "XIC files ready")
      if (values$mobilogram_available) badge_text <- paste0(badge_text, " + IM")
      format_label <- if (values$xic_format == "v2") "DIA-NN 2.x" else "DIA-NN 1.x"
      div(class = "alert alert-success py-1 px-2 mt-2 mb-0",
        style = "font-size: 0.85em;",
        icon("check-circle"),
        badge_text, br(),
        span(class = "text-muted", style = "font-size: 0.85em;", format_label))
    } else {
      NULL
    }
  })

  # --- XIC Modal Function ---
  show_xic_modal <- function(session, values) {
    showModal(modalDialog(
      title = div(style = "display: flex; align-items: center; gap: 10px;",
        icon("wave-square"),
        span("XIC Chromatograms:"),
        span(values$xic_protein, style = "color: #667eea; font-weight: bold;"),
        span(textOutput("xic_precursor_count", inline = TRUE),
          style = "font-size: 0.8em; color: #6c757d; margin-left: 10px;")
      ),

      # Controls row
      div(style = "display: flex; flex-wrap: wrap; gap: 12px; align-items: flex-end; margin-bottom: 12px;",
        selectInput("xic_display_mode", "Display:",
          choices = c("Facet by sample" = "overlay",
                      "Facet by fragment" = "facet"),
          selected = "overlay", width = "220px"),

        selectInput("xic_precursor_select", "Precursor:",
          choices = NULL, width = "280px"),

        selectInput("xic_group_filter", "Filter Group:",
          choices = NULL, width = "180px"),

        checkboxInput("xic_show_ms1", "Show MS1 (split axis)", value = FALSE),

        # Ion mobility toggle — only shown when mobilogram data is available
        if (values$mobilogram_available) {
          div(style = "display: flex; align-items: center; gap: 6px; padding: 4px 10px; border-radius: 6px; background: linear-gradient(135deg, #e0f2fe, #bae6fd); border: 1px solid #7dd3fc;",
            checkboxInput("xic_show_mobilogram", "Ion Mobility", value = FALSE),
            icon("bolt", style = "color: #0284c7; font-size: 1.1em;")
          )
        }
      ),

      # Mobilogram mode banner — visible when IM is active
      uiOutput("xic_mobilogram_banner"),

      # Main plot
      plotlyOutput("xic_plot", height = "calc(100vh - 380px)"),

      # Info panel
      div(class = "mt-2 p-2 bg-light rounded",
        style = "font-size: 0.85em;",
        uiOutput("xic_info_panel")
      ),

      size = "xl",
      easyClose = TRUE,
      footer = div(style = "display: flex; gap: 8px; align-items: center;",
        actionButton("xic_prev_protein", label = "Prev",
          icon = icon("arrow-left"), class = "btn-outline-secondary btn-sm"),
        actionButton("xic_next_protein", label = "Next",
          icon = icon("arrow-right"), class = "btn-outline-secondary btn-sm"),
        downloadButton("xic_download_plot", "Download",
          class = "btn-outline-success btn-sm"),
        modalButton("Close")
      )
    ))
  }

  # --- XIC from DE Dashboard ---
  observeEvent(input$show_xic, {
    if (!values$xic_available) {
      showNotification("Load XIC files first (sidebar > '5. XIC Viewer')", type = "warning")
      return()
    }
    if (is.null(values$plot_selected_proteins) || length(values$plot_selected_proteins) == 0) {
      showNotification("Select a protein in the Volcano Plot or Table first!", type = "warning")
      return()
    }
    values$xic_protein <- values$plot_selected_proteins[1]
    show_xic_modal(session, values)
  })

  # --- XIC from Grid View ---
  observeEvent(input$show_xic_from_grid, {
    if (!values$xic_available) {
      showNotification("Load XIC files first (sidebar > '5. XIC Viewer')", type = "warning")
      return()
    }
    req(values$grid_selected_protein)
    values$xic_protein <- values$grid_selected_protein
    removeModal()
    show_xic_modal(session, values)
  })

  # --- Reactive XIC Data Loading ---
  observe({
    req(values$xic_protein, values$xic_available, values$xic_dir,
        values$xic_report_map)

    withProgress(message = "Loading chromatograms...", {
      incProgress(0.3, detail = "Reading XIC files...")

      tryCatch({
        xic_raw <- load_xic_for_protein(
          xic_dir = values$xic_dir,
          protein_id = values$xic_protein,
          report_map = values$xic_report_map,
          xic_format = values$xic_format
        )

        if (!is.null(xic_raw) && nrow(xic_raw) > 0) {
          incProgress(0.6, detail = "Reshaping data...")

          values$xic_data <- reshape_xic_for_plotting(
            xic_raw, values$metadata, xic_format = values$xic_format)

          # Update precursor selector
          precursors <- unique(values$xic_data$Precursor.Id)
          updateSelectInput(session, "xic_precursor_select",
            choices = c("All Precursors" = "all",
                        setNames(precursors, precursors)),
            selected = "all")

          # Update group filter
          groups <- unique(na.omit(values$xic_data$Group))
          if (length(groups) == 0) groups <- unique(values$metadata$Group)
          updateSelectInput(session, "xic_group_filter",
            choices = c("All Groups" = "all", setNames(groups, groups)),
            selected = "all")
        } else {
          values$xic_data <- NULL
          showNotification(
            paste("No XIC data found for", values$xic_protein,
                  "-- this protein may have been identified via MBR."),
            type = "warning", duration = 5)
        }
      }, error = function(e) {
        values$xic_data <- NULL
        showNotification(paste("Error loading XICs:", e$message),
          type = "error", duration = 8)
      })
    })
  })

  # --- XIC Precursor Count ---
  output$xic_precursor_count <- renderText({
    req(values$xic_data)
    n_prec <- n_distinct(values$xic_data$Precursor.Id)
    n_frag <- values$xic_data %>%
      filter(MS.Level == 2) %>%
      distinct(Fragment.Label) %>%
      nrow()
    paste0(n_prec, " precursor(s), ", n_frag, " fragments")
  })

  # --- Mobilogram Mode Banner ---
  output$xic_mobilogram_banner <- renderUI({
    if (isTRUE(input$xic_show_mobilogram)) {
      div(style = "background: linear-gradient(135deg, #0284c7, #0369a1); color: white; padding: 8px 16px; border-radius: 8px; margin-bottom: 8px; display: flex; align-items: center; gap: 10px; font-weight: 500;",
        icon("bolt", style = "font-size: 1.3em;"),
        span("Ion Mobility Mode", style = "font-size: 1.05em;"),
        span("— Showing mobilogram data (1/K0 vs intensity). Requires timsTOF or PASEF instrument data.",
          style = "font-weight: 400; font-size: 0.9em; opacity: 0.9;")
      )
    } else {
      NULL
    }
  })

  # --- XIC Plot ---
  output$xic_plot <- renderPlotly({
    req(values$xic_data)

    xic <- values$xic_data
    display_mode <- input$xic_display_mode

    # Filter by selected precursor
    if (!is.null(input$xic_precursor_select) &&
        input$xic_precursor_select != "all") {
      xic <- xic %>% filter(Precursor.Id == input$xic_precursor_select)
    }

    # Filter by group
    if (!is.null(input$xic_group_filter) &&
        input$xic_group_filter != "all") {
      xic <- xic %>% filter(Group == input$xic_group_filter)
    }

    # MS level filter — when MS1 shown, add panel label for split axis
    show_ms1 <- isTRUE(input$xic_show_ms1)
    if (!show_ms1) {
      xic_plot <- xic %>% filter(MS.Level == 2)
      if (nrow(xic_plot) == 0) xic_plot <- xic  # fallback if only MS1 exists
    } else {
      xic_plot <- xic %>%
        mutate(MS_Panel = factor(
          ifelse(MS.Level == 1, "MS1 Precursor", "MS2 Fragments"),
          levels = c("MS1 Precursor", "MS2 Fragments")))
    }

    # Cap at top 6 precursors for very large proteins
    if (n_distinct(xic_plot$Precursor.Id) > 6 &&
        (is.null(input$xic_precursor_select) || input$xic_precursor_select == "all")) {
      top_prec <- xic_plot %>%
        filter(MS.Level == 2) %>%
        count(Precursor.Id) %>%
        slice_max(n, n = 6) %>%
        pull(Precursor.Id)
      # Keep MS1 for selected precursors too
      xic_plot <- xic_plot %>% filter(Precursor.Id %in% top_prec)
      showNotification(
        paste("Showing top 6 of", n_distinct(xic$Precursor.Id),
              "precursors. Use the selector to view others."),
        type = "message", duration = 4)
    }

    # Sample label for faceting
    xic_plot <- xic_plot %>%
      mutate(
        Sample.Label = paste0(ID, " (", Group, ")"),
        tooltip_text = paste0(
          "<b>Sample:</b> ", ID, " (", Group, ")",
          "<br><b>Fragment:</b> ", Fragment.Label,
          "<br><b>RT:</b> ", round(RT, 3), " min",
          "<br><b>Intensity:</b> ", format(round(Intensity), big.mark = ",")
        )
      )

    if (display_mode == "overlay") {
      if (show_ms1) {
        # Split axis: MS1 on top, MS2 fragments below, per sample
        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Fragment.Label,
                group = interaction(File.Name, Fragment.Label),
                text = tooltip_text)) +
          geom_line(alpha = 0.7, linewidth = 0.5) +
          facet_grid(MS_Panel ~ Sample.Label, scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            subtitle = "MS1 precursor (top) vs MS2 fragments (bottom) per sample",
            x = "Retention Time (min)", y = "Intensity", color = "Fragment"
          ) +
          theme(legend.position = "bottom", legend.text = element_text(size = 7),
                strip.text = element_text(size = 8),
                strip.text.y = element_text(angle = 0, face = "bold"))
      } else {
        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Fragment.Label,
                group = interaction(File.Name, Fragment.Label),
                text = tooltip_text)) +
          geom_line(alpha = 0.7, linewidth = 0.5) +
          facet_wrap(~ Sample.Label, scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            subtitle = "Each panel = one sample, colors = fragment ions",
            x = "Retention Time (min)", y = "Intensity", color = "Fragment"
          ) +
          theme(legend.position = "bottom", legend.text = element_text(size = 7),
                strip.text = element_text(size = 8))
      }

    } else if (display_mode == "facet") {
      if (show_ms1) {
        # Split axis: MS1 on top row, each fragment in its own panel below
        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Group,
                group = File.Name, text = tooltip_text)) +
          geom_line(alpha = 0.6, linewidth = 0.4) +
          facet_grid(MS_Panel ~ Fragment.Label, scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            subtitle = "MS1 precursor (top) vs MS2 fragments (bottom), colors = groups",
            x = "Retention Time (min)", y = "Intensity", color = "Group"
          ) +
          theme(strip.text = element_text(size = 7),
                strip.text.y = element_text(angle = 0, face = "bold"))
      } else {
        p <- ggplot(xic_plot,
            aes(x = RT, y = Intensity, color = Group,
                group = File.Name, text = tooltip_text)) +
          geom_line(alpha = 0.6, linewidth = 0.4) +
          facet_wrap(~Fragment.Label, scales = "free_y") +
          theme_minimal() +
          labs(
            title = paste("Fragment XICs --", values$xic_protein),
            subtitle = "Each panel = one fragment ion, colors = groups",
            x = "Retention Time (min)", y = "Intensity", color = "Group"
          ) +
          theme(strip.text = element_text(size = 8))
      }

    }

    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(orientation = "h", y = -0.15),
             margin = list(b = 80)) %>%
      config(displayModeBar = TRUE)
  })

  # --- XIC Info Panel ---
  output$xic_info_panel <- renderUI({
    req(values$xic_data, values$xic_protein)

    xic <- values$xic_data
    n_precursors <- n_distinct(xic$Precursor.Id)
    n_fragments <- xic %>% filter(MS.Level == 2) %>%
      distinct(Fragment.Label) %>% nrow()
    n_samples <- n_distinct(xic$File.Name)
    rt_range <- range(xic$RT, na.rm = TRUE)

    # Get DE stats for this protein from current contrast
    de_info <- tryCatch({
      tt <- topTable(values$fit, coef = input$contrast_selector, number = Inf)
      tt[values$xic_protein, c("logFC", "adj.P.Val")]
    }, error = function(e) NULL)

    tagList(
      div(class = "row",
        div(class = "col-md-4",
          strong("Protein: "), values$xic_protein, br(),
          strong("Precursors: "), n_precursors, br(),
          strong("Fragment ions: "), n_fragments
        ),
        div(class = "col-md-4",
          strong("Samples: "), n_samples, br(),
          strong("RT range: "),
          paste0(round(rt_range[1], 2), " -- ", round(rt_range[2], 2), " min")
        ),
        div(class = "col-md-4",
          if (!is.null(de_info)) {
            tagList(
              strong("log2 FC: "),
              span(round(de_info$logFC, 3),
                style = paste0("color: ",
                  ifelse(de_info$logFC > 0, "#d32f2f", "#1976d2"), ";")),
              br(),
              strong("adj. p-value: "),
              formatC(de_info$adj.P.Val, format = "e", digits = 2)
            )
          } else {
            em("DE stats unavailable for current contrast")
          }
        )
      ),
      hr(style = "margin: 5px 0;"),
      div(class = "text-muted", style = "font-size: 0.8em;",
        icon("info-circle"),
        " Co-eluting fragment ions with similar peak shapes indicate reliable identification.",
        " Consistent peak areas across replicates within a group support accurate quantification.",
        " Irregular peaks or missing fragments may indicate interference or low-confidence IDs."
      )
    )
  })

  # --- XIC Protein Navigation (Prev/Next) ---
  observeEvent(input$xic_prev_protein, {
    req(values$xic_protein, values$fit, input$contrast_selector)
    tt <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>%
      filter(adj.P.Val < 0.05) %>% arrange(adj.P.Val)
    prot_list <- rownames(tt)
    current_idx <- which(prot_list == values$xic_protein)
    if (length(current_idx) > 0 && current_idx > 1) {
      values$xic_protein <- prot_list[current_idx - 1]
      values$plot_selected_proteins <- values$xic_protein
    }
  })

  observeEvent(input$xic_next_protein, {
    req(values$xic_protein, values$fit, input$contrast_selector)
    tt <- topTable(values$fit, coef = input$contrast_selector, number = Inf) %>%
      filter(adj.P.Val < 0.05) %>% arrange(adj.P.Val)
    prot_list <- rownames(tt)
    current_idx <- which(prot_list == values$xic_protein)
    if (length(current_idx) > 0 && current_idx < length(prot_list)) {
      values$xic_protein <- prot_list[current_idx + 1]
      values$plot_selected_proteins <- values$xic_protein
    }
  })

  # --- XIC Download ---
  output$xic_download_plot <- downloadHandler(
    filename = function() {
      paste0("XIC_", gsub("[^A-Za-z0-9]", "_", values$xic_protein), "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      req(values$xic_data)
      xic <- values$xic_data

      if (!isTRUE(input$xic_show_ms1)) {
        xic_plot <- xic %>% filter(MS.Level == 2)
        if (nrow(xic_plot) == 0) xic_plot <- xic
      } else {
        xic_plot <- xic
      }

      p <- ggplot(xic_plot,
          aes(x = RT, y = Intensity, color = Fragment.Label,
              group = interaction(File.Name, Fragment.Label))) +
        geom_line(alpha = 0.7, linewidth = 0.5) +
        facet_wrap(~ paste0(ID, " (", Group, ")"), scales = "free_y") +
        theme_minimal() +
        labs(
          title = paste("Fragment XICs --", values$xic_protein),
          x = "Retention Time (min)", y = "Intensity", color = "Fragment"
        ) +
        theme(legend.position = "bottom", legend.text = element_text(size = 7),
              strip.text = element_text(size = 8))

      ggsave(file, plot = p, width = 14, height = 10, dpi = 150)
    }
  )
}

shinyApp(ui, server)