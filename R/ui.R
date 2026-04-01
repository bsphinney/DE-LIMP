# ==============================================================================
#  USER INTERFACE (UI) — Build the complete app UI
#  Called from app.R as: ui <- build_ui(is_hf_space, search_enabled, docker_available, hpc_available, local_sbatch)
# ==============================================================================

build_ui <- function(is_hf_space, search_enabled = FALSE,
                     docker_available = FALSE, hpc_available = FALSE,
                     local_sbatch = FALSE, local_diann = FALSE,
                     delimp_data_dir = "",
                     is_core_facility = FALSE, cf_config = NULL,
                     deploy_env = "Local",
                     config = list(), is_hive = FALSE) {

  # Environment badge colors
  env_colors <- list(
    Docker = "#e74c3c",       # red
    HPC    = "#27ae60",       # green
    Local  = "#3498db",       # blue
    `Hugging Face` = "#f39c12" # orange
  )
  env_color <- env_colors[[deploy_env]] %||% "#6c757d"

  page_navbar(
  title = tags$span(
    "DE-LIMP Proteomics",
    tags$span(
      deploy_env,
      style = sprintf(
        "font-size: 0.65em; background: %s; color: white; padding: 2px 8px; border-radius: 10px; margin-left: 8px; vertical-align: middle; font-weight: 500;",
        env_color
      )
    )
  ),
  window_title = paste0("DE-LIMP (", deploy_env, ")"),
  id = "main_tabs",
  theme = bs_theme(bootswatch = "flatly"),
  navbar_options = navbar_options(bg = "#2c3e50"),
  header = tagList(
    useShinyjs(),
    tags$head(tags$style(HTML("
    /* Custom properties matching mockup */
    :root {
      --flatly-primary: #2c3e50;
      --flatly-info: #3498db;
      --flatly-success: #18bc9c;
      --flatly-body: #f5f7f9;
      --flatly-border: #dee2e6;
      --flatly-muted: #6c757d;
    }
    /* Methodology text: wrap long lines instead of horizontal overflow */
    #methodology_text {
      white-space: pre-wrap;
      word-wrap: break-word;
    }

    /* Navbar dropdown styling (open/close handled by JS) */
    .navbar .dropdown-menu {
      border-radius: 0 0 6px 6px;
      box-shadow: 0 6px 20px rgba(0,0,0,0.12);
      min-width: 230px;
    }

    /* Force white text on dark navbar */
    .navbar .nav-link,
    .navbar .navbar-nav .nav-link {
      color: rgba(255,255,255,0.75) !important;
    }
    .navbar .nav-link:hover,
    .navbar .navbar-nav .nav-link:hover {
      color: #ffffff !important;
    }
    .navbar .nav-link.active,
    .navbar .navbar-nav .nav-link.active {
      color: #ffffff !important;
      border-bottom: 3px solid var(--flatly-success);
    }
    .navbar .navbar-brand {
      color: #ffffff !important;
    }

    /* Ensure hidden nav items are truly invisible (progressive reveal) */
    .navbar .nav-item[style*='display: none'],
    .navbar .nav-item.d-none {
      width: 0 !important;
      overflow: hidden !important;
    }

    /* Card consistency */
    .card { border-radius: 8px; }

    /* Sidebar accordion refinement */
    .sidebar .accordion-button {
      padding: 10px 12px;
      font-size: 0.82rem;
      font-weight: 600;
    }

    .chat-container { height: 500px; overflow-y: auto; border: 1px solid #ddd; padding: 15px; background-color: #f8f9fa; border-radius: 5px; margin-bottom: 15px; }
    .user-msg { text-align: right; margin: 10px 0; }
    .user-msg span { background-color: #007bff; color: white; padding: 8px 12px; border-radius: 15px 15px 0 15px; display: inline-block; max-width: 80%; }
    .ai-msg { text-align: left; margin: 10px 0; }
    .ai-msg span { background-color: #e9ecef; color: #333; padding: 8px 12px; border-radius: 15px 15px 15px 0; display: inline-block; max-width: 80%; }
    .selection-banner { background-color: #d4edda; color: #155724; padding: 10px; border-radius: 5px; margin-bottom: 10px; font-weight: bold; border: 1px solid #c3e6cb; }

    /* === RESPONSIVE UI ADDITIONS === */

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

    /* Remove default bottom margin from selectInput inside gradient banners */
    div[style*='linear-gradient'] .form-group,
    div[style*='linear-gradient'] .shiny-input-container {
      margin-bottom: 0 !important;
    }

    /* Ensure selectize dropdowns render above plots and card bodies */
    .selectize-dropdown {
      z-index: 10000 !important;
    }
    .tab-pane, .tab-content {
      overflow: visible !important;
    }
    /* Compact covariate checkboxes and inputs */
    .covariate-row .checkbox { margin-top: 0; margin-bottom: 0; }
    .covariate-row .form-group { margin-bottom: 0; }
    .covariate-row .form-control { padding: 4px 8px; height: auto; font-size: 0.85em; }

    /* Re-enable horizontal scrolling for DataTables inside tab panes */
    .dataTables_wrapper {
      overflow-x: auto !important;
      overflow-y: visible;
      max-width: 100%;
      box-sizing: border-box;
    }
  "))),

  tags$head(tags$script(HTML("
    // Resize plotly charts when tabs or modals become visible
    function resizePlotlyAll() {
      var plots = document.querySelectorAll('.js-plotly-plot');
      plots.forEach(function(gd) {
        if (gd.offsetParent !== null && gd.data) {
          Plotly.Plots.resize(gd);
        }
      });
    }
    $(document).on('shown.bs.modal', function() {
      setTimeout(resizePlotlyAll, 200);
    });
    $(document).on('shown.bs.tab', function() {
      setTimeout(resizePlotlyAll, 150);
      setTimeout(resizePlotlyAll, 500);
    });

    // Navbar dropdown: hover-to-open using Bootstrap 5 Dropdown API
    $(document).ready(function() {
      var closeTimer = null;

      // Open dropdown on hover (event delegation for DOM safety)
      $(document).on('mouseenter', '.navbar .nav-item.dropdown', function() {
        clearTimeout(closeTimer);
        var toggle = $(this).children('.dropdown-toggle')[0];
        if (!toggle) return;
        // Close all OTHER dropdowns via Bootstrap API
        $('.navbar .nav-item.dropdown').not(this).each(function() {
          var t = $(this).children('.dropdown-toggle')[0];
          if (t) { var i = bootstrap.Dropdown.getInstance(t); if (i) i.hide(); }
        });
        // Open THIS dropdown via Bootstrap API
        bootstrap.Dropdown.getOrCreateInstance(toggle).show();
      });

      // Close dropdown on mouse leave (with small delay for UX)
      $(document).on('mouseleave', '.navbar .nav-item.dropdown', function() {
        var toggle = $(this).children('.dropdown-toggle')[0];
        if (!toggle) return;
        closeTimer = setTimeout(function() {
          var inst = bootstrap.Dropdown.getInstance(toggle);
          if (inst) inst.hide();
        }, 150);
      });

      // Close ALL dropdowns when a menu item is clicked
      $(document).on('click', '.navbar .dropdown-menu .dropdown-item', function() {
        clearTimeout(closeTimer);
        $('.navbar .nav-item.dropdown .dropdown-toggle').each(function() {
          var inst = bootstrap.Dropdown.getInstance(this);
          if (inst) inst.hide();
        });
      });

      // Prevent default click-toggle on dropdown buttons (hover handles open/close)
      $(document).on('click', '.navbar .nav-item.dropdown > .dropdown-toggle', function(e) {
        e.stopPropagation();
        e.preventDefault();
      });
    });

    // Inject section labels into Analysis dropdown
    $(document).ready(function() {
      setTimeout(function() {
        var items = $('a.dropdown-item');
        items.each(function() {
          var text = $(this).text().trim();
          if (text === 'Data Overview') {
            $('<h6 class=\"dropdown-header\" style=\"font-size:0.72rem;font-weight:700;text-transform:uppercase;letter-spacing:0.06em;\">Setup</h6>').insertBefore($(this));
          }
          if (text === 'DE Dashboard') {
            $('<div class=\"dropdown-divider\"></div><h6 class=\"dropdown-header\" style=\"font-size:0.72rem;font-weight:700;text-transform:uppercase;letter-spacing:0.06em;\">Results</h6>').insertBefore($(this));
          }
          if (text === 'AI Analysis') {
            $('<div class=\"dropdown-divider\"></div><h6 class=\"dropdown-header\" style=\"font-size:0.72rem;font-weight:700;text-transform:uppercase;letter-spacing:0.06em;\">AI</h6>').insertBefore($(this));
          }
        });
      }, 300);
    });
  "))),
  ),  # end header

  sidebar = sidebar(
    width = 300,

    # --- Top links ---
    div(style="display: flex; gap: 5px; margin-bottom: 10px;",
      tags$a(href="https://github.com/bsphinney/DE-LIMP/blob/main/USER_GUIDE.md", target="_blank", class="btn btn-info w-50", icon("book"), "Guide", style="color:white; font-weight:bold;"),
      tags$a(href="https://github.com/bsphinney/DE-LIMP", target="_blank", class="btn btn-secondary w-50", icon("github"), "Code", style="color:white; font-weight:bold;")
    ),

    # --- Acquisition Mode Switcher (DIA / DDA / XL-MS) ---
    # Gated by config.yml feature flags + HPC availability
    {
      mode_choices <- c("DIA" = "dia")
      if (is_hive && isTRUE(config$features$enable_dda)) {
        mode_choices <- c(mode_choices, "DDA" = "dda")
      }
      if (is_hive && isTRUE(config$features$enable_xlms)) {
        mode_choices <- c(mode_choices, "XL-MS" = "xlms")
      }
      if (length(mode_choices) > 1) {
        div(
          class = "acquisition-mode-switcher",
          style = paste(
            "display: flex; align-items: center; gap: 12px;",
            "padding: 10px 16px; border-radius: 8px;",
            "background: linear-gradient(135deg, #f8fafc 0%, #e8f4fd 100%);",
            "border: 1px solid #c8dff0; margin-bottom: 16px;"
          ),
          div(
            style = paste(
              "font-size: 11px; font-weight: 600; color: #6c757d;",
              "letter-spacing: 0.05em; text-transform: uppercase;"
            ),
            "Mode"
          ),
          radioButtons(
            inputId  = "acquisition_mode",
            label    = NULL,
            choices  = mode_choices,
            selected = "dia",
            inline   = TRUE
          ),
          uiOutput("mode_context_label")
        )
      } else {
        # Single mode (DIA only) — hidden input so conditionalPanel works
        div(style = "display:none;",
          radioButtons("acquisition_mode", NULL, choices = c("DIA" = "dia"), selected = "dia")
        )
      }
    },

    # --- Main accordion ---
    accordion(
      id = "sidebar_sections", multiple = TRUE, open = "Upload Data",

      accordion_panel("Upload Data", icon = icon("file-arrow-up"),
        fileInput("report_file", "DIA-NN Report (.parquet)", accept = c(".parquet")),
        if (!is_hf_space) conditionalPanel("output.ssh_connected_flag == true",
          actionButton("load_from_hpc_btn", "Load from HPC",
            icon = icon("server"), class = "btn-outline-info btn-sm w-100",
            style = "margin-bottom: 5px;")
        ),
        actionButton("load_example", "\U0001F4CA Load Example Data", class = "btn-info btn-sm w-100",
                     style = "margin-bottom: 5px;"),
        actionButton("load_example_phospho", "Load Example Phospho Data",
          class = "btn-outline-info btn-sm w-100", icon = icon("flask"),
          style = "margin-bottom: 10px;"),
        tags$details(
          tags$summary(icon("dna"), " De Novo Sequencing (Cascadia)"),
          fileInput("ssl_files", "Cascadia SSL Files (.ssl)",
            multiple = TRUE, accept = ".ssl"),
          sliderInput("denovo_score_threshold", "Min. Confidence",
            min = 0.5, max = 1.0, value = 0.8, step = 0.05),
          checkboxInput("denovo_enable", "Enable Cascadia integration", value = FALSE)
        ),
        numericInput("q_cutoff", "Q-Value Cutoff", value = 0.01, min = 0, max = 0.1, step = 0.01)
      ),

      accordion_panel("Pipeline Settings", icon = icon("sliders"),
        sliderInput("logfc_cutoff", "Min Log2 Fold Change:", min=0, max=5, value=0.6, step=0.1)
      ),

      accordion_panel("AI Chat", icon = icon("robot"),
        passwordInput("user_api_key", "Gemini API Key", value = "", placeholder = "AIzaSy..."),
        tags$details(style = "margin: 5px 0 10px 0; font-size: 0.85em; color: #6c757d;",
          tags$summary(style = "cursor: pointer; color: #17a2b8;", "How to get a free API key"),
          tags$ol(style = "margin-top: 5px; padding-left: 20px;",
            tags$li("Go to ", tags$a("Google AI Studio", href = "https://aistudio.google.com/apikey",
              target = "_blank", rel = "noopener noreferrer")),
            tags$li("Sign in with your Google account"),
            tags$li('Click "Create API Key"'),
            tags$li("Copy and paste the key above")
          ),
          tags$p(style = "margin-bottom: 0;", "The free tier is sufficient for DE-LIMP.")
        ),
        actionButton("check_models", "Check Models", class="btn-warning btn-xs w-100"),
        br(), br(),
        textInput("model_name", "Model Name", value = "gemini-2.5-flash", placeholder = "gemini-2.5-flash")
      )
    ),

    # --- Session buttons (outside accordion) ---
    div(style="margin-top: 10px;",
      div(style="display: flex; gap: 5px; margin-bottom: 5px;",
        downloadButton("save_session", "Save", class = "btn-primary w-50", icon = icon("download")),
        actionButton("load_session_btn", "Load", class = "btn-outline-primary w-50", icon = icon("upload"))
      ),
      actionButton("prepare_next_btn", "Prepare Next Analysis",
        class = "btn-outline-secondary btn-sm w-100",
        icon = icon("broom"),
        style = "margin-top: 4px;")
    ),

    # Core facility report link
    if (is_core_facility) tagList(
      hr(),
      uiOutput("report_link_ui")
    ),

    # Core Facility: Templates
    if (is_core_facility) tagList(
      hr(),
      h5(icon("bookmark"), " Templates"),
      selectInput("template_selector", NULL,
        choices = c("(none)" = ""), selected = ""),
      div(style = "display: flex; gap: 5px;",
        actionButton("load_template", "Apply", class = "btn-sm btn-outline-primary w-50"),
        actionButton("save_template", "Save Current", class = "btn-sm btn-outline-success w-50")
      )
    ),

    # Phospho controls (conditional on detection)
    conditionalPanel(
      condition = "output.phospho_detected_flag",
      accordion(
        id = "sidebar_phospho", open = "Phosphoproteomics",
        accordion_panel("Phosphoproteomics", icon = icon("flask"),
          radioButtons("phospho_input_mode", "Site Quantification Source",
            choices = c(
              "DIA-NN site matrix (recommended)" = "site_matrix",
              "Parse from report.parquet" = "parsed_report"
            ),
            selected = "site_matrix"
          ),
          conditionalPanel(
            condition = "input.phospho_input_mode == 'site_matrix'",
            fileInput("phospho_site_matrix_file", "Upload site matrix (.tsv or .parquet)",
                      accept = c(".tsv", ".parquet", ".txt")),
            tags$p(class = "text-muted small",
              "Upload the site localization matrix from DIA-NN (TSV or parquet format)."
            )
          ),
          conditionalPanel(
            condition = "input.phospho_input_mode == 'parsed_report'",
            sliderInput("phospho_loc_threshold", "Site Localization Confidence",
              min = 0.5, max = 1.0, value = 0.75, step = 0.05),
            tags$p(class = "text-muted small",
              "Recommended: 0.75 for exploratory, 0.9 for high-confidence sites."
            )
          ),
          radioButtons("phospho_norm", "Site-Level Normalization",
            choices = c(
              "None (DIA-NN normalized)" = "none",
              "Median centering" = "median",
              "Quantile normalization" = "quantile"
            ),
            selected = "none"
          ),
          actionButton("run_phospho_pipeline", "Run Phosphosite Analysis",
                       class = "btn-warning w-100", icon = icon("bolt")),
          hr(),
          tags$p(class = "text-muted small", style = "margin-top: 8px;",
            tags$strong("Advanced (Phase 2/3):")
          ),
          fileInput("phospho_fasta_file", "Upload FASTA (for motifs)",
                    accept = c(".fasta", ".fa", ".faa")),
          tags$p(class = "text-muted small",
            "Protein FASTA enables accurate motif extraction around phosphosites."
          ),
          checkboxInput("phospho_protein_correction",
            "Normalize to protein abundance", value = FALSE),
          tags$p(class = "text-muted small",
            "Subtracts protein-level logFC from site logFC (requires total proteome pipeline to be run first)."
          )
        )
      )
    ),

    # XIC Viewer (conditional on !is_hf_space)
    if (!is_hf_space) accordion(
      id = "sidebar_xic",
      accordion_panel("XIC Viewer", icon = icon("wave-square"),
        p(class = "text-muted small",
          "Load .xic.parquet files from DIA-NN to inspect chromatograms."),
        textInput("xic_dir_input", "XIC Directory Path:",
          placeholder = "Auto-detected or paste path here"),
        actionButton("xic_load_dir", "Load XICs", class = "btn-outline-info btn-sm w-100",
          icon = icon("wave-square")),
        uiOutput("xic_status_badge")
      )
    ),
    if (is_hf_space) div(
      style = "padding: 8px; margin-top: 4px; background: linear-gradient(135deg, #e0f2fe, #f0f9ff); border: 1px solid #bae6fd; border-radius: 8px; font-size: 0.82em;",
      icon("chart-line", style = "color: #0284c7;"),
      span(style = "font-weight: 600; color: #0c4a6e;", " XIC Viewer"),
      p(style = "margin: 4px 0 0 0; color: #475569;",
        "Fragment-level chromatogram inspection is available when running DE-LIMP locally or on HPC.",
        tags$a(href = "https://github.com/bsphinney/DE-LIMP", target = "_blank", " Download here."))
    )
  ),

  # ============================================================================
  #  MAIN CONTENT — nav items directly inside page_navbar
  #  Layout: Search | QC | Analysis v | Output v | Education | Facility v
  # ============================================================================

    # ==========================================================================
    # SEARCH (standalone, conditional) — visible when Docker or HPC backend available
    # ==========================================================================
    if (search_enabled) nav_panel("New Search", icon = icon("rocket"),

      # --- DIA panel (existing — shown when acquisition_mode === 'dia') ---
      conditionalPanel(
        condition = "input.acquisition_mode === 'dia'",

      # Three-panel wizard layout
      layout_column_wrap(
        width = 1/3,

        # === PANEL 1: FILES ===
        card(
          card_header(tagList(icon("folder-open"), " 1. Files")),
          card_body(
            style = "overflow-y: auto; max-height: calc(100vh - 200px);",

            textInput("analysis_name", "Analysis Name",
              value = "",
              placeholder = "e.g., HeLa_DIA_2026"),
            textInput("search_notes", NULL,
              value = "",
              placeholder = "Notes (optional) — e.g., testing sprot+iso FASTA..."),

            # Core facility: lab, instrument, project for search tracking
            if (is_core_facility) tagList(
              div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
                div(style = "flex: 1; min-width: 130px;",
                  selectInput("search_lab", "Lab",
                    choices = c("(select)" = "", cf_lab_names(cf_config)),
                    selected = "")
                ),
                div(style = "flex: 1; min-width: 130px;",
                  selectInput("search_instrument", "Instrument",
                    choices = c("(select)" = "", cf_instrument_names(cf_config)),
                    selected = "")
                )
              ),
              selectizeInput("search_project", "Project",
                choices = NULL, selected = NULL,
                options = list(create = TRUE,
                               placeholder = "Select or type new project..."))
            ),

            hr(),
            tags$h6(icon("hard-drive"), " Raw Data"),
            # Local file browser: shown for Docker backend OR local-sbatch HPC
            conditionalPanel(
              "input.search_backend == 'local' || input.search_backend == 'docker' || (input.search_backend == 'hpc' && input.search_connection_mode != 'ssh')",
              shinyFiles::shinyDirButton("raw_data_dir", "Select Raw Data Folder",
                title = "Choose directory with .d / .raw / .mzML files",
                class = "btn-outline-primary btn-sm w-100")
            ),
            # SSH remote path: only for HPC + SSH mode
            conditionalPanel(
              "input.search_backend == 'hpc' && input.search_connection_mode == 'ssh'",
              div(style = "display: flex; gap: 5px;",
                div(style = "flex: 1;",
                  textInput("ssh_raw_data_dir", NULL,
                    value = "",
                    placeholder = "/quobyte/proteomics-grp/raw/experiment")
                ),
                actionButton("ssh_browse_raw_btn", NULL, icon = icon("folder-open"),
                  class = "btn-outline-primary btn-sm",
                  style = "margin-top: 0;", title = "Browse remote directories"),
                actionButton("ssh_scan_raw_btn", "Scan", icon = icon("magnifying-glass"),
                  class = "btn-outline-secondary btn-sm", style = "margin-top: 0;")
              )
            ),
            uiOutput("raw_file_summary"),
            uiOutput("tic_extract_ui"),

            hr(),
            tags$h6(icon("dna"), " FASTA Database"),
            selectInput("fasta_source", NULL,
              choices = c("Download from UniProt" = "uniprot",
                          "Download from NCBI" = "ncbi",
                          "Database Library" = "library",
                          "Pre-staged on server" = "prestaged",
                          "Browse / enter path" = "browse"),
              width = "100%"),

            # --- UniProt source ---
            conditionalPanel("input.fasta_source == 'uniprot'",
              actionButton("open_uniprot_modal", "Search UniProt",
                class = "btn-info btn-sm w-100", icon = icon("search")),
              uiOutput("fasta_filename_preview"),
              uiOutput("fasta_selected_summary"),
              uiOutput("fasta_add_to_library_btn_ui")
            ),

            # --- NCBI source ---
            conditionalPanel("input.fasta_source == 'ncbi'",
              actionButton("open_ncbi_modal", "Search NCBI",
                class = "btn-success btn-sm w-100", icon = icon("search")),
              uiOutput("ncbi_fasta_selected_summary")
            ),

            # --- Database Library source ---
            conditionalPanel("input.fasta_source == 'library'",
              actionButton("open_fasta_library_modal", "Browse Speclib Library",
                class = "btn-primary btn-sm w-100", icon = icon("bolt")),
              uiOutput("fasta_library_selected_summary")
            ),

            # --- Pre-staged source ---
            conditionalPanel("input.fasta_source == 'prestaged'",
              selectInput("prestaged_fasta", "Available Databases:",
                choices = NULL, width = "100%"),
              uiOutput("prestaged_fasta_info")
            ),

            # --- Browse / path source ---
            conditionalPanel("input.fasta_source == 'browse'",
              conditionalPanel(
                "input.search_backend == 'local' || input.search_backend == 'docker' || (input.search_backend == 'hpc' && input.search_connection_mode != 'ssh')",
                shinyFiles::shinyDirButton("fasta_browse_dir", "Browse for FASTA Folder",
                  title = "Navigate to FASTA directory",
                  class = "btn-outline-primary btn-sm w-100")
              ),
              conditionalPanel(
                "input.search_backend == 'hpc' && input.search_connection_mode == 'ssh'",
                div(style = "display: flex; gap: 5px;",
                  div(style = "flex: 1;",
                    textInput("ssh_fasta_browse_dir", NULL,
                      placeholder = "/share/proteomics/fasta/")
                  ),
                  actionButton("ssh_browse_fasta_btn", NULL, icon = icon("folder-open"),
                    class = "btn-outline-primary btn-sm",
                    style = "margin-top: 0;", title = "Browse remote directories"),
                  actionButton("ssh_scan_fasta_btn", "Scan", icon = icon("magnifying-glass"),
                    class = "btn-outline-secondary btn-sm", style = "margin-top: 0;")
                )
              ),
              uiOutput("browsed_fasta_info")
            ),

            div(style = "margin-top: 10px;",
              selectInput("contaminant_library", "Add Contaminant Library:",
                choices = c(
                  "None" = "none",
                  "Universal (Recommended)" = "universal",
                  "Cell Culture" = "cell_culture",
                  "Mouse Tissue" = "mouse_tissue",
                  "Rat Tissue" = "rat_tissue",
                  "Neuron Culture" = "neuron_culture",
                  "Stem Cell Culture" = "stem_cell_culture"
                ),
                selected = "universal", width = "100%"),
              tags$small(class = "text-muted",
                "Contaminant libraries from ",
                tags$a(href = "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics",
                       "HaoGroup-ProtContLib", target = "_blank"))
            ),

            div(style = "margin-top: 10px;",
              textAreaInput("custom_fasta_sequences",
                "Custom Protein Sequences (FASTA format):",
                placeholder = ">sp|CUSTOM1|My_Protein\nMSEQUENCE...",
                rows = 3, width = "100%"),
              tags$small(class = "text-muted",
                "Paste FASTA-formatted sequences for custom proteins, tagged constructs, etc. ",
                "These are added to your search alongside the main database.")
            ),

            hr(),
            tags$h6(icon("book"), " Spectral Library (optional)"),
            conditionalPanel(
              "input.search_backend == 'local' || input.search_backend == 'docker' || (input.search_backend == 'hpc' && input.search_connection_mode != 'ssh')",
              shinyFiles::shinyFilesButton("lib_file", "Select .speclib File",
                title = "Choose spectral library",
                class = "btn-outline-secondary btn-sm w-100",
                multiple = FALSE)
            ),
            conditionalPanel(
              "input.search_backend == 'hpc' && input.search_connection_mode == 'ssh'",
              textInput("ssh_lib_file", NULL,
                placeholder = "/share/proteomics/libraries/my.speclib")
            ),
            uiOutput("lib_file_info")
          )
        ),

        # === PANEL 2: SEARCH SETTINGS (shared across backends) ===
        card(
          card_header(tagList(icon("sliders"), " 2. Search Settings")),
          card_body(
            style = "overflow-y: auto; max-height: calc(100vh - 200px);",

            tags$div(
              style = "margin-bottom: 12px; padding: 8px 10px; background: #f0f4f8; border-radius: 6px; border: 1px solid #dee2e6;",
              div(style = "display: flex; align-items: center; gap: 8px;",
                div(style = "flex: 1;",
                  fileInput("diann_log_file", NULL,
                    placeholder = "Import from DIA-NN log...",
                    width = "100%")
                ),
                actionButton("import_log_info_btn", icon("question-circle"),
                  class = "btn-outline-info btn-sm",
                  style = "margin-top: -18px;",
                  title = "How to use log import")
              ),
              uiOutput("log_import_feedback")
            ),

            radioButtons("search_mode", "Search Mode:",
              choices = c(
                "Library-free (default)" = "libfree",
                "Phosphoproteomics" = "phospho",
                "Use spectral library" = "library"
              ), selected = "libfree"
            ),
            conditionalPanel("input.search_mode == 'phospho'",
              tags$div(class = "alert alert-info py-1 px-2 mb-2",
                style = "font-size: 0.82em;",
                icon("flask"),
                " Phospho mode: STY phosphorylation (UniMod:21), max 3 var mods,",
                " 2 missed cleavages, --phospho-output enabled"
              )
            ),

            radioButtons("diann_normalization", "DIA-NN Normalization:",
              choices = c(
                "RT-dependent (default)" = "on",
                "Off (for AP-MS / Co-IP)" = "off"
              ), selected = "on"
            ),
            uiOutput("norm_guidance_search"),

            # Core facility: LC / Gradient method
            if (is_core_facility) tagList(
              hr(),
              selectInput("search_lc_method", "LC / Gradient Method",
                choices = c("(select)" = "", cf_lc_method_names(cf_config)),
                selected = "")
            ),

            hr(),
            tags$h6("Basic Parameters"),
            div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
              div(style = "flex: 1; min-width: 120px;",
                selectInput("diann_enzyme", "Enzyme:",
                  choices = c("Trypsin/P (K*,R*)" = "K*,R*",
                              "Trypsin strict" = "K,R",
                              "Lys-C" = "K",
                              "Chymotrypsin" = "F,W,Y,L",
                              "None" = "-"),
                  selected = "K*,R*")
              ),
              div(style = "flex: 1; min-width: 100px;",
                numericInput("diann_missed_cleavages", "Missed Cleavages:",
                  value = 1, min = 0, max = 5)
              )
            ),

            div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
              div(style = "flex: 1; min-width: 120px;",
                numericInput("diann_max_var_mods", "Max Variable Mods:",
                  value = 1, min = 0, max = 5)
              ),
              div(style = "flex: 1; min-width: 120px;",
                selectInput("mass_acc_mode", "Mass Accuracy:",
                  choices = c("Automatic" = "auto", "Manual" = "manual"),
                  selected = "auto")
              )
            ),

            conditionalPanel(
              condition = "input.mass_acc_mode == 'auto'",
              tags$p(class = "text-muted", style = "font-size: 0.78em; margin-top: -4px;",
                "DIA-NN will determine optimal MS1/MS2 mass accuracy from the data.")
            ),
            conditionalPanel(
              condition = "input.mass_acc_mode == 'manual'",
              div(style = "display: flex; gap: 8px;",
                div(style = "flex: 1;",
                  numericInput("diann_mass_acc", "MS2 (ppm):", value = 14, min = 1, max = 50)
                ),
                div(style = "flex: 1;",
                  numericInput("diann_mass_acc_ms1", "MS1 (ppm):", value = 14, min = 1, max = 50)
                )
              ),
              uiOutput("mass_acc_hint")
            ),

            tags$h6("Variable Modifications"),
            checkboxInput("mod_met_ox", "Methionine oxidation (UniMod:35)", TRUE),
            checkboxInput("mod_nterm_acetyl", "N-term acetylation (UniMod:1)", FALSE),
            textAreaInput("extra_var_mods", "Additional Mods (one per line):",
              placeholder = "e.g., UniMod:21,79.966331,STY", rows = 2),

            # Advanced accordion
            accordion(
              id = "search_advanced_accordion",
              open = FALSE,
              accordion_panel(
                "Advanced Options",
                icon = icon("gear"),

                tags$h6("Peptide/Precursor Ranges"),
                div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
                  div(style = "flex: 1; min-width: 100px;",
                    numericInput("min_pep_len", "Min Peptide:", value = 7, min = 4, max = 15)
                  ),
                  div(style = "flex: 1; min-width: 100px;",
                    numericInput("max_pep_len", "Max Peptide:", value = 30, min = 15, max = 52)
                  )
                ),
                div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
                  div(style = "flex: 1; min-width: 100px;",
                    numericInput("min_pr_mz", "Min m/z:", value = 300, min = 100, max = 500)
                  ),
                  div(style = "flex: 1; min-width: 100px;",
                    numericInput("max_pr_mz", "Max m/z:", value = 1800, min = 800, max = 2000)
                  )
                ),

                tags$h6("Processing Options"),
                checkboxInput("diann_mbr", "Match Between Runs (MBR)", TRUE),
                checkboxInput("diann_rt_profiling", "RT profiling", TRUE),
                checkboxInput("diann_xic", "Generate XICs", TRUE),
                checkboxInput("diann_unimod4", "UniMod4 (carbamidomethylation)", TRUE),
                checkboxInput("diann_met_excision", "N-term methionine excision", TRUE),

                numericInput("diann_scan_window", "Scan Window:", value = 6, min = 0, max = 20),

                textAreaInput("extra_cli_flags", "Extra DIA-NN Flags:",
                  placeholder = "e.g., --peptidoforms --min-pr-charge 2", rows = 2),

                numericInput("diann_fdr", "FDR Threshold:", value = 0.01, min = 0.001, max = 0.1, step = 0.005)
              )
            )
          )
        ),

        # === PANEL 3: RESOURCES & SUBMIT ===
        card(
          card_header(tagList(icon("server"), " 3. Resources & Submit")),
          card_body(
            style = "overflow-y: auto; max-height: calc(100vh - 200px);",

            # Backend selector
            tags$h6(icon("microchip"), " Compute Backend"),
            {
              backends <- c()
              if (local_diann)      backends <- c(backends, "Local (Embedded)" = "local")
              if (docker_available) backends <- c(backends, "Local (Docker)" = "docker")
              if (hpc_available)    backends <- c(backends, "HPC (SSH/SLURM)" = "hpc")
              # Default to HPC when SSH key is available (Docker + HPC is the common setup)
              ssh_key_detected <- nzchar(Sys.getenv("DELIMP_SSH_KEY", "")) ||
                file.exists("/tmp/.ssh/id_ed25519") ||
                (nzchar(delimp_data_dir) && length(list.files(
                  file.path(delimp_data_dir, "ssh"), pattern = "^[^.]", full.names = FALSE)) > 0)
              default_val <- if (hpc_available && ssh_key_detected) "hpc"
                            else unname(backends[1])

              if (length(backends) > 1) {
                radioButtons("search_backend", NULL,
                  choices = backends, selected = default_val, inline = TRUE)
              } else {
                # Single backend — use real Shiny input (hidden) so conditionalPanel works
                div(style = "display: none;",
                  radioButtons("search_backend", NULL,
                    choices = backends, selected = default_val))
              }
            },

            # ---------- Local (embedded) backend controls ----------
            conditionalPanel("input.search_backend == 'local'",
              tags$div(class = "alert alert-success py-1 px-2",
                style = "font-size: 0.85em;",
                icon("check-circle"),
                " DIA-NN binary detected. Ready for local searches."),
              hr(),
              tags$h6("Resources"),
              uiOutput("local_resources_ui"),
              hr(),
              tags$h6("Output Directory"),
              if (nzchar(delimp_data_dir)) {
                # Container mode: fixed output path
                textInput("local_output_dir", NULL,
                  value = file.path(delimp_data_dir, "output"))
              } else {
                # Native mode: file browser
                tagList(
                  shinyFiles::shinyDirButton("local_output_dir_browse",
                    "Select Output Folder",
                    title = "Choose output directory for DIA-NN results",
                    class = "btn-outline-primary btn-sm w-100"),
                  verbatimTextOutput("local_output_path")
                )
              }
            ),

            # ---------- Docker backend controls ----------
            conditionalPanel("input.search_backend == 'docker'",
              uiOutput("docker_image_status"),
              hr(),
              tags$h6("Docker Resources"),
              uiOutput("docker_resources_ui"),
              hr(),
              tags$h6("Output Directory"),
              shinyFiles::shinyDirButton("docker_output_dir", "Select Output Folder",
                title = "Choose output directory for DIA-NN results",
                class = "btn-outline-primary btn-sm w-100"),
              verbatimTextOutput("docker_output_path"),
              textInput("docker_image_name", "DIA-NN Docker Image:",
                value = "diann:2.0")
            ),

            # ---------- HPC backend controls ----------
            conditionalPanel("input.search_backend == 'hpc'",

              # Core facility: staff selector replaces manual SSH fields
              if (is_core_facility) tagList(
                tags$h6(icon("user"), " Staff Identity"),
                selectInput("staff_selector", NULL,
                  choices = c("(select)" = "", cf_staff_names(cf_config))
                ),
                uiOutput("staff_connection_status"),
                hr()
              ),

              tags$h6(icon("plug"), " Connection Mode"),
              radioButtons("search_connection_mode", NULL,
                choices = c("Local (on HPC)" = "local", "Remote (SSH)" = "ssh"),
                selected = if (local_sbatch || nzchar(Sys.getenv("APPTAINER_CONTAINER", "")) ||
                               nzchar(Sys.getenv("SINGULARITY_CONTAINER", ""))) "local"
                           else "ssh",
                inline = TRUE),
              conditionalPanel("input.search_connection_mode == 'ssh'",
                # SSH fields: shown always (core facility auto-fills them from staff selector)
                if (!is_core_facility) tagList(
                  textInput("ssh_host", "HPC Hostname",
                    value = "hive.hpc.ucdavis.edu"),
                  div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
                    div(style = "flex: 1; min-width: 100px;",
                      textInput("ssh_user", "Username",
                        value = Sys.getenv("DELIMP_SSH_USER", ""))
                    ),
                    div(style = "flex: 1; min-width: 80px;",
                      numericInput("ssh_port", "Port", value = 22, min = 1, max = 65535)
                    )
                  ),
                  div(style = "display: flex; gap: 5px; align-items: flex-end;",
                    div(style = "flex: 1;",
                      textInput("ssh_key_path", "SSH Key Path",
                        value = {
                          # Auto-detect SSH key: env var > Docker mount > data/ssh > ~/.ssh
                          env_key <- Sys.getenv("DELIMP_SSH_KEY", "")
                          if (nzchar(env_key) && file.exists(env_key)) env_key
                          else if (file.exists("/home/shiny/.ssh/id_ed25519")) "/home/shiny/.ssh/id_ed25519"
                          else if (nzchar(delimp_data_dir) && file.exists(paste0(delimp_data_dir, "/ssh/id_ed25519")))
                            paste0(delimp_data_dir, "/ssh/id_ed25519")
                          else paste0(Sys.getenv("HOME"), "/.ssh/id_ed25519")
                        })
                    ),
                    shinyFiles::shinyFilesButton("ssh_key_browse", "Browse",
                      title = "Select SSH private key", multiple = FALSE,
                      class = "btn-outline-secondary btn-sm", style = "margin-bottom: 15px;")
                  ),
                  textInput("ssh_modules", "Modules to Load (optional)",
                    value = "",
                    placeholder = "e.g., slurm apptainer"),
                  actionButton("test_ssh_btn", "Test Connection",
                    icon = icon("plug"), class = "btn-outline-info btn-sm"),
                  uiOutput("ssh_status_ui")
                ),
                # Hidden SSH inputs for core facility mode (auto-filled by staff selector)
                if (is_core_facility) div(style = "display: none;",
                  textInput("ssh_host", "HPC Hostname", value = ""),
                  textInput("ssh_user", "Username", value = ""),
                  numericInput("ssh_port", "Port", value = 22, min = 1, max = 65535),
                  textInput("ssh_key_path", "SSH Key Path", value = ""),
                  textInput("ssh_modules", "Modules to Load", value = "")
                )
              ),
              # Cluster status indicator — outside SSH conditionalPanel so it shows
              # in both "Local (on HPC)" mode (via SLURM proxy) and "Remote (SSH)" mode
              uiOutput("cluster_status_ui"),

              hr(),
              tags$h6("SLURM Resources"),
              div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
                div(style = "flex: 1; min-width: 100px;",
                  numericInput("diann_cpus", "CPUs:", value = 64, min = 4, max = 128, step = 4)
                ),
                div(style = "flex: 1; min-width: 100px;",
                  numericInput("diann_mem_gb", "Memory (GB):", value = 128, min = 16, max = 1024, step = 16)
                )
              ),
              div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
                div(style = "flex: 1; min-width: 100px;",
                  numericInput("diann_time_hours", "Time (hrs):", value = 12, min = 1, max = 48)
                )
              ),
              # Auto-select partition/account display + override
              uiOutput("partition_selector_ui"),
              # Hidden textInputs — actual values used by job submission
              div(style = "display: none;",
                textInput("diann_partition", NULL, value = "high"),
                textInput("diann_account", NULL, value = "genome-center-grp")
              ),

              # Cluster Monitor — usage history + grant reporting
              accordion(
                id = "cluster_monitor_accordion",
                open = FALSE,
                accordion_panel(
                  "Cluster Monitor",
                  icon = icon("chart-line"),
                  uiOutput("cluster_capacity_alert"),
                  div(style = "display: flex; align-items: center; gap: 8px; margin-bottom: 6px;",
                    radioButtons("cluster_history_range", NULL,
                      choices = c("24h" = "24", "7d" = "168", "30d" = "720", "All" = "0"),
                      selected = "168", inline = TRUE),
                    actionButton("cluster_monitor_expand_btn", icon("expand"),
                      class = "btn-outline-primary btn-xs", style = "padding: 1px 5px;",
                      title = "Open in full window"),
                    actionButton("cluster_monitor_info_btn", icon("question-circle"),
                      class = "btn-outline-info btn-xs", style = "padding: 1px 5px;")
                  ),
                  plotlyOutput("cluster_usage_chart", height = "260px"),
                  tags$h6("Group Members", style = "margin-top: 12px; margin-bottom: 4px;"),
                  plotlyOutput("per_user_chart", height = "200px"),
                  hr(style = "margin: 10px 0;"),
                  div(style = "display: flex; align-items: center; gap: 8px;",
                    checkboxInput("auto_queue_switch", "Auto-switch pending jobs to publicgrp/low",
                      value = TRUE, width = "100%")
                  ),
                  conditionalPanel("input.auto_queue_switch",
                    div(style = "display: flex; align-items: center; gap: 6px; margin-bottom: 8px;",
                      tags$small("After waiting"),
                      numericInput("queue_wait_minutes", NULL, value = 5, min = 1, max = 60,
                        step = 1, width = "65px"),
                      tags$small("min")
                    )
                  ),
                  div(style = "margin-top: 8px; text-align: right;",
                    downloadButton("export_cluster_csv", "Export for Grant",
                      class = "btn-outline-primary btn-sm", icon = icon("file-csv"))
                  )
                )
              ),

              # Parallel search mode (rendered server-side based on file count)
              uiOutput("parallel_mode_ui"),

              hr(),
              tags$h6("DIA-NN Container"),
              textInput("diann_sif_path", "Apptainer SIF Path:",
                value = "/quobyte/proteomics-grp/dia-nn/diann_2.3.0.sif",
                placeholder = "/path/to/diann_2.3.0.sif"),

              hr(),
              conditionalPanel("input.search_connection_mode != 'ssh'",
                tags$h6("Output Directory"),
                shinyFiles::shinyDirButton("output_base_dir", "Select Output Folder",
                  title = "Choose base output directory",
                  class = "btn-outline-primary btn-sm w-100"),
                verbatimTextOutput("full_output_path")
              ),
              # SSH mode: output dir auto-generated from raw data dir, keep input hidden
              div(style = "display: none;",
                textInput("ssh_output_base_dir", NULL,
                  value = "",
                  placeholder = "/share/proteomics/results/")
              )
            ),

            hr(),
            uiOutput("time_estimate_ui"),

            conditionalPanel("output.ssh_connected_flag == true",
              checkboxInput("add_cascadia_denovo",
                tagList(icon("dna"), " Add de novo sequencing (Cascadia)"),
                value = FALSE),
              conditionalPanel("input.add_cascadia_denovo",
                tags$small(class = "text-muted", style = "display: block; margin: -8px 0 8px 25px;",
                  "Runs Cascadia on GPU in parallel with DIA-NN. Results appear in the De Novo dropdown.")
              )
            ),

            actionButton("submit_diann", "Submit DIA-NN Search",
              class = "btn-success btn-lg w-100",
              icon = icon("rocket"),
              style = "margin-top: 10px;"),

            checkboxInput("auto_load_results", "Auto-load results when complete", FALSE),

            hr(),
            div(style = "display: flex; justify-content: space-between; align-items: center;",
              tags$h6(icon("list-check"), " Job Queue", style = "margin-bottom: 0;"),
              actionButton("recover_jobs_btn", "Recover",
                class = "btn-outline-info btn-sm",
                style = "font-size: 0.75em; padding: 2px 8px;",
                icon = icon("magnifying-glass"))
            ),
            uiOutput("search_queue_ui"),

            # License attribution
            tags$div(class = "text-muted", style = "font-size: 0.78em; margin-top: 12px; border-top: 1px solid #dee2e6; padding-top: 8px;",
              "DIA-NN by Vadim Demichev. ",
              tags$a(href = "https://github.com/vdemichev/DiaNN/blob/master/LICENSE.md",
                     "License", target = "_blank"), " | ",
              "Demichev V et al. (2020) ", tags$em("Nature Methods"), " 17:41-44")
          )
        )
      )
      ), # end conditionalPanel DIA

      # --- DDA panel (shown when acquisition_mode === 'dda') ---
      conditionalPanel(
        condition = "input.acquisition_mode === 'dda'",
        div(
          style = "padding: 20px; max-width: 900px; margin: 0 auto;",

          # Info banner
          div(class = "alert alert-info", role = "alert",
            icon("flask"), " ",
            tags$strong("DDA Workflow"),
            " powered by ",
            tags$strong("Sage"), " (database search) + ",
            tags$strong("Casanovo"), " (de novo, optional). ",
            "Requires Hive HPC connection. timsTOF ddaPASEF .d files only."
          ),

          # SSH required warning
          conditionalPanel(
            condition = "!output.ssh_connected_flag",
            div(class = "alert alert-warning", role = "alert",
              icon("exclamation-triangle"),
              " Connect to Hive via SSH (in the DIA Search tab) before submitting a DDA search."
            )
          ),

          # --- Raw Files ---
          div(
            style = "background: white; border: 1px solid #dee2e6; border-radius: 8px; padding: 16px; margin-bottom: 16px;",
            tags$h6(icon("folder-open"), " Raw Files", style = "margin-bottom: 12px;"),
            div(style = "display: flex; gap: 8px; align-items: flex-end;",
              div(style = "flex: 1;",
                textInput("dda_raw_dir", "Raw file directory (Hive path)",
                  placeholder = "/quobyte/proteomics-grp/to-hive/mass-spec-archive/...",
                  width = "100%")
              ),
              actionButton("dda_scan_files", "Scan",
                icon = icon("search"), class = "btn-outline-primary btn-sm",
                style = "margin-bottom: 15px;")
            ),
            uiOutput("dda_file_list_preview")
          ),

          # --- Database ---
          div(
            style = "background: white; border: 1px solid #dee2e6; border-radius: 8px; padding: 16px; margin-bottom: 16px;",
            tags$h6(icon("database"), " Database", style = "margin-bottom: 12px;"),
            textInput("dda_fasta_path", "FASTA path (Hive)",
              placeholder = "/quobyte/proteomics-grp/brett/data/uniprotkb_proteome_UP000005640_2026_01_14.fasta",
              width = "100%"),
            tags$small(style = "color: #6c757d;",
              "Recommended: one-protein-per-gene (OPG) FASTA for cleaner protein inference.")
          ),

          # --- Search Parameters ---
          div(
            style = "background: white; border: 1px solid #dee2e6; border-radius: 8px; padding: 16px; margin-bottom: 16px;",
            tags$h6(icon("sliders"), " Search Parameters", style = "margin-bottom: 12px;"),
            div(style = "display: flex; gap: 16px; flex-wrap: wrap;",
              div(style = "min-width: 180px;",
                textInput("dda_experiment_name", "Experiment name",
                  value = "dda_search", width = "100%")
              ),
              div(style = "min-width: 150px;",
                selectInput("dda_preset", "Preset",
                  choices = c("Standard" = "standard",
                              "Phosphoproteomics" = "phospho",
                              "TMT Labeling" = "tmt"),
                  selected = "standard", width = "100%")
              )
            ),
            div(style = "display: flex; gap: 16px; flex-wrap: wrap; margin-top: 8px;",
              div(style = "min-width: 130px;",
                numericInput("dda_missed_cleavages", "Missed cleavages",
                  value = 2, min = 0, max = 4, step = 1, width = "100%")
              ),
              div(style = "min-width: 130px;",
                numericInput("dda_precursor_tol", "Precursor tol. (ppm)",
                  value = 20, min = 1, max = 100, step = 1, width = "100%")
              ),
              div(style = "min-width: 130px;",
                numericInput("dda_fragment_tol", "Fragment tol. (Da)",
                  value = 0.05, min = 0.01, max = 0.5, step = 0.01, width = "100%")
              )
            ),
            # Advanced SLURM controls (collapsed)
            tags$details(
              style = "margin-top: 12px;",
              tags$summary(style = "cursor: pointer; color: #6c757d; font-size: 13px;",
                icon("gear"), " Advanced SLURM settings"),
              div(style = "display: flex; gap: 16px; flex-wrap: wrap; margin-top: 8px;",
                div(style = "min-width: 100px;",
                  numericInput("dda_cpus", "CPUs", value = 32, min = 4, max = 128, step = 4, width = "100%")
                ),
                div(style = "min-width: 100px;",
                  numericInput("dda_mem", "Memory (GB)", value = 64, min = 8, max = 512, step = 8, width = "100%")
                ),
                div(style = "min-width: 120px;",
                  textInput("dda_time_limit", "Time limit", value = "02:00:00", width = "100%")
                )
              )
            )
          ),

          # --- Normalization & Imputation ---
          div(
            style = "background: white; border: 1px solid #dee2e6; border-radius: 8px; padding: 16px; margin-bottom: 16px;",
            tags$h6(icon("chart-line"), " Normalization & Imputation", style = "margin-bottom: 12px;"),
            div(style = "display: flex; gap: 16px; flex-wrap: wrap;",
              div(style = "min-width: 180px;",
                selectInput("dda_norm_method", "Normalization",
                  choices = c("Cyclic Loess" = "cyclicloess",
                              "Median centering" = "median",
                              "Quantile" = "quantile",
                              "None" = "none"),
                  selected = "cyclicloess", width = "100%")
              ),
              div(style = "min-width: 180px;",
                selectInput("dda_impute_method", "Imputation",
                  choices = c("Perseus (MNAR)" = "perseus",
                              "MinProb (MNAR)" = "minprob",
                              "MinDet (deterministic)" = "mindet",
                              "None" = "none"),
                  selected = "perseus", width = "100%")
              ),
              div(style = "min-width: 130px;",
                sliderInput("dda_min_valid", "Min. valid fraction",
                  min = 0.3, max = 1.0, value = 0.5, step = 0.1, width = "100%")
              )
            ),
            # Perseus parameters (conditional)
            conditionalPanel(
              condition = "input.dda_impute_method === 'perseus' || input.dda_impute_method === 'minprob'",
              div(style = "display: flex; gap: 16px; flex-wrap: wrap; margin-top: 8px;",
                div(style = "min-width: 130px;",
                  numericInput("dda_perseus_width", "Width", value = 0.3, min = 0.1, max = 1.0, step = 0.1, width = "100%")
                ),
                div(style = "min-width: 130px;",
                  numericInput("dda_perseus_shift", "Shift (SD)", value = 1.8, min = 0.5, max = 3.0, step = 0.1, width = "100%")
                )
              )
            )
          ),

          # --- Casanovo (disabled) ---
          div(
            style = "background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 8px; padding: 16px; margin-bottom: 16px; opacity: 0.7;",
            tags$h6(icon("wand-magic-sparkles"), " De Novo Sequencing (Casanovo)", style = "margin-bottom: 12px; color: #6c757d;"),
            div(
              style = "display: flex; align-items: center; gap: 8px;",
              tags$input(type = "checkbox", disabled = "disabled",
                style = "width: 16px; height: 16px;"),
              tags$span(style = "color: #6c757d;", "Run Casanovo de novo sequencing"),
              tags$span(
                class = "badge bg-secondary",
                style = "font-size: 10px; margin-left: 8px;",
                "Coming soon"
              )
            ),
            tags$small(style = "color: #adb5bd; display: block; margin-top: 8px;",
              "Casanovo GPU-accelerated de novo sequencing will be available in a future update.")
          ),

          # --- Submit ---
          conditionalPanel(
            condition = "output.ssh_connected_flag",
            div(
              style = "padding: 8px 0;",
              actionButton("run_dda_search", "Submit Sage Search",
                icon  = icon("rocket"),
                class = "btn-primary btn-lg",
                width = "100%"
              )
            )
          ),

          # --- Status + Results ---
          uiOutput("dda_job_status_ui"),

          # --- Group Assignment + Pipeline (shown after results loaded) ---
          uiOutput("dda_group_assignment_ui"),

          uiOutput("dda_results_summary_ui"),

          # --- QC Summary Card ---
          uiOutput("dda_qc_summary_card")
        )
      ),

      # --- XL-MS panel (stub — shown when acquisition_mode === 'xlms') ---
      conditionalPanel(
        condition = "input.acquisition_mode === 'xlms'",
        div(
          style = "padding: 40px; max-width: 800px; margin: 0 auto;",
          div(class = "alert alert-info", role = "alert",
            icon("diagram-project"), " ",
            tags$strong("Cross-Linking Mass Spectrometry"),
            " powered by ",
            tags$strong("MeroX"), " + ",
            tags$strong("xiSearch"), " for protein-protein interaction mapping."
          ),
          div(
            style = paste(
              "text-align: center; padding: 60px 20px;",
              "background: #f8f9fa; border-radius: 12px;",
              "border: 2px dashed #dee2e6;"
            ),
            icon("link", style = "font-size: 48px; color: #6f42c1; margin-bottom: 16px;"),
            tags$h4("XL-MS Pipeline Coming Soon", style = "color: #495057;"),
            tags$p(
              style = "color: #6c757d; max-width: 500px; margin: 0 auto;",
              "The MeroX + xiSearch cross-linking pipeline will appear here. ",
              "Crosslink identification, FDR control, and protein-protein ",
              "interaction network visualization."
            )
          )
        )
      )

    ),

    # ==========================================================================
    # QC (standalone, merged from QC Trends + QC Plots)
    # ==========================================================================
    nav_panel("QC", icon = icon("chart-bar"),
              # DDA QC Summary Card (shown when in DDA mode)
              uiOutput("dda_qc_summary_card"),
              # Global sort order control (from QC Trends)
              div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
                div(style = "display: flex; align-items: center; justify-content: space-between;",
                  div(style = "display: flex; align-items: center; gap: 15px;",
                    icon("sort", style = "color: #6c757d;"),
                    strong("Sort Order:"),
                    radioButtons("qc_sort_order", NULL,
                      choices = c("Run Order", "Group"),
                      inline = TRUE,
                      selected = "Run Order"
                    ),
                    span(style = "color: #6c757d; font-size: 0.85em;",
                      "(Applies to Sample Metrics)")
                  ),
                  actionButton("qc_trends_info_btn", icon("question-circle"), title = "What are QC Trends?",
                    class = "btn-outline-info btn-sm")
                )
              ),

              # Merged QC sub-tabs: Sample Metrics + Diagnostics
              navset_card_tab(
                id = "qc_merged_tabs",

                # ── Data Completeness (Detected vs Inferred proteins) ──
                nav_panel("Data Completeness",
                  icon = icon("eye"),
                  div(style = "overflow-y: auto; max-height: calc(100vh - 250px);",
                    # Warning banner (conditional)
                    uiOutput("completeness_warning_banner"),
                    # Summary cards
                    uiOutput("completeness_summary_cards"),
                    # Info button row
                    div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                      actionButton("completeness_info_btn", icon("question-circle"),
                        title = "About Data Completeness", class = "btn-outline-info btn-sm")
                    ),
                    # 1. Detected vs Inferred stacked bar
                    tags$h5("Detected vs Inferred Proteins per Sample", style = "margin-top: 8px;"),
                    plotlyOutput("completeness_stacked_bar", height = "400px"),
                    tags$hr(),
                    # 2. Precursor Evidence Heatmap
                    tags$h5("Precursor Evidence Heatmap (Top 50 Most Variable)", style = "margin-top: 12px;"),
                    plotlyOutput("completeness_evidence_heatmap", height = "500px"),
                    tags$hr(),
                    # 3. Cumulative Detection Curve
                    tags$h5("Cumulative Detection Curve", style = "margin-top: 12px;"),
                    plotlyOutput("completeness_cumulative_curve", height = "350px"),
                    tags$hr(),
                    # 4. Sample Clustering by Detection Pattern
                    tags$h5("Sample Clustering by Detection Pattern (Jaccard Distance)", style = "margin-top: 12px;"),
                    plotlyOutput("completeness_dendrogram", height = "400px"),
                    tags$hr(),
                    # 5. Precursor Evidence Distribution
                    tags$h5("Precursor Count per Protein (per Sample)", style = "margin-top: 12px;"),
                    plotlyOutput("completeness_precursor_violin", height = "400px")
                  )
                ),

                # ── Sample Metrics (faceted: Precursors, Proteins, MS1 Signal) ──
                nav_panel("Sample Metrics",
                  icon = icon("chart-line"),
                  div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                    actionButton("qc_metrics_info_btn", icon("question-circle"),
                      title = "About Sample Metrics", class = "btn-outline-info btn-sm"),
                    actionButton("fullscreen_qc_metrics", "\U0001F50D Fullscreen",
                      class = "btn-outline-secondary btn-sm")
                  ),
                  plotlyOutput("qc_metrics_trend", height = "calc(100vh - 380px)")
                ),

                # ── Chromatography QC (TIC traces, run diagnostics) ──
                nav_panel("Chromatography QC",
                  icon = icon("chart-area"),
                  div(style = "overflow-y: auto; max-height: calc(100vh - 250px);",
                    conditionalPanel(
                      condition = "typeof output.tic_qc_has_data !== 'undefined' && output.tic_qc_has_data",
                      uiOutput("tic_qc_status_badges"),
                      div(style = "display: flex; justify-content: space-between; align-items: center; margin: 8px 0;",
                        div(style = "display: flex; gap: 15px; align-items: center;",
                          radioButtons("tic_view_mode", NULL,
                            choices = c("Faceted" = "faceted", "Overlay" = "overlay", "Metrics" = "metrics"),
                            inline = TRUE
                          ),
                          conditionalPanel(
                            condition = "input.tic_view_mode == 'faceted'",
                            radioButtons("tic_facet_mode", NULL,
                              choices = c("By Run" = "run", "By Group" = "group"),
                              selected = "run", inline = TRUE)
                          )
                        ),
                        div(style = "display: flex; gap: 5px;",
                          actionButton("tic_qc_info_btn", icon("question-circle"),
                            title = "About Chromatography QC", class = "btn-outline-info btn-sm"),
                          actionButton("tic_qc_fullscreen_btn", "\U0001F50D Fullscreen",
                            class = "btn-outline-secondary btn-sm")
                        )
                      ),
                      div(
                        uiOutput("tic_qc_plot_container")
                      ),
                      conditionalPanel(
                        condition = "input.tic_view_mode == 'metrics'",
                        div(style = "margin-top: 12px;",
                          DTOutput("tic_metrics_table")
                        )
                      ),
                      conditionalPanel(
                        condition = "input.tic_view_mode != 'faceted'",
                        uiOutput("tic_qc_diagnostics")
                      )
                    ),
                    conditionalPanel(
                      condition = "typeof output.tic_qc_has_data === 'undefined' || !output.tic_qc_has_data",
                      div(class = "alert alert-info", style = "margin-top: 20px;",
                        icon("info-circle"),
                        " No TIC data available. Extract from the Search tab before searching, or use the button below to extract from an existing output directory.",
                        div(style = "margin-top: 8px;",
                          actionButton("tic_extract_from_qc_btn", "Extract TIC from Raw Files",
                            class = "btn-outline-primary btn-sm", icon = icon("chart-area"))
                        )
                      )
                    )
                  )
                ),

                nav_panel("Stats Table",
                  icon = icon("table"),
                  div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                    actionButton("qc_stats_info_btn", icon("question-circle"), title = "QC Statistics",
                      class = "btn-outline-info btn-sm"),
                    downloadButton("download_qc_stats_csv", tagList(icon("download"), " CSV"),
                      class = "btn-success btn-sm")
                  ),
                  DTOutput("r_qc_table")
                ),

                # ── Diagnostics (from QC Plots) ──
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
                      div(style = "display: flex; gap: 8px;",
                        actionButton("norm_diag_info_btn", icon("question-circle"), title = "What am I looking at?",
                          class = "btn-outline-info btn-sm"),
                        actionButton("fullscreen_norm_diag", "\U0001F50D Fullscreen",
                          class = "btn-outline-secondary btn-sm")
                      )
                    ),

                    # Plot with viewport height
                    plotlyOutput("norm_diagnostic_plot", height = "calc(100vh - 340px)")
                  )
                ),

                nav_panel("DPC Fit",
                  icon = icon("chart-line"),
                  card_body(
                    div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                      actionButton("dpc_info_btn", icon("question-circle"), title = "What is DPC Fit?",
                        class = "btn-outline-info btn-sm"),
                      actionButton("fullscreen_dpc", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),
                    plotOutput("dpc_plot", height = "70vh")
                  )
                ),

                nav_panel("MDS Plot",
                  icon = icon("project-diagram"),
                  card_body(
                    div(style = "display: flex; justify-content: flex-end; gap: 8px; align-items: center; margin-bottom: 10px;",
                      span("Color by:", style = "font-weight: 500; font-size: 0.85em; color: #555;"),
                      div(style = "width: 160px;",
                        selectInput("mds_color_by", label = NULL,
                          choices = c("Group", "Batch"), selected = "Group", width = "100%")
                      ),
                      actionButton("mds_info_btn", icon("question-circle"), title = "What is MDS?",
                        class = "btn-outline-info btn-sm"),
                      actionButton("fullscreen_mds", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),
                    plotOutput("mds_plot", height = "70vh")
                  )
                ),

                nav_panel("Group Distribution",
                  icon = icon("chart-area"),
                  card_body(
                    div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                      selectInput("qc_violin_metric", "Metric:",
                        choices = c("Precursors", "Proteins", "MS1_Signal"),
                        width = "200px"
                      ),
                      div(style = "display: flex; gap: 8px;",
                        actionButton("group_dist_info_btn", icon("question-circle"), title = "What is this?",
                          class = "btn-outline-info btn-sm"),
                        actionButton("fullscreen_qc_violin", "\U0001F50D Fullscreen",
                          class = "btn-outline-secondary btn-sm")
                      )
                    ),
                    plotlyOutput("qc_group_violin", height = "calc(100vh - 320px)")
                  )
                ),

                nav_panel("P-value Distribution",
                  icon = icon("chart-column"),
                  # Comparison selector banner — two-row layout for full-width dropdown
                  div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 12px 15px; border-radius: 8px; margin-bottom: 15px; position: relative; z-index: 10;",
                    div(style = "display: flex; align-items: center; justify-content: space-between; margin-bottom: 8px;",
                      div(style = "display: flex; align-items: center; gap: 10px;",
                        icon("microscope"),
                        span("Viewing Comparison:", style = "font-weight: 500;")
                      ),
                      div(style = "display: flex; gap: 6px;",
                        actionButton("pvalue_hist_info_btn", icon("question-circle"), title = "How do I interpret this?",
                          class = "btn-outline-light btn-sm"),
                        actionButton("fullscreen_pvalue_hist", "\U0001F50D Fullscreen",
                          class = "btn-outline-light btn-sm")
                      )
                    ),
                    selectInput("contrast_selector_pvalue", NULL,
                      choices = NULL,
                      width = "100%"
                    )
                  ),
                  # Plot (plain div — card_body creates stacking context that clips dropdown)
                  plotOutput("pvalue_histogram", height = "calc(100vh - 400px)"),

                  # Automated contextual guidance (below plot)
                  uiOutput("pvalue_guidance")
                )
              )
    ),

    # ==========================================================================
    # ANALYSIS dropdown — Data Overview, DE, Phospho, GSEA, MOFA2, AI Analysis
    # ==========================================================================
    nav_menu("Analysis", icon = icon("microscope"),

      # -- Setup section --
      nav_panel("Data Overview", icon = icon("database"),
                # Data views as tabs
                navset_card_tab(
                  id = "data_overview_tabs",

                  nav_panel("Assign Groups & Run",
                    icon = icon("table"),
                    # Phospho detection banner (shown when phospho data detected)
                    uiOutput("phospho_detection_banner"),
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

                      # Covariates (compact single-row)
                      div(style="flex: 1; min-width: 250px;",
                        strong("Covariates:", style="font-size: 0.85em; display: block; margin-bottom: 4px;"),
                        div(class="covariate-row", style="display: flex; flex-wrap: wrap; gap: 6px; align-items: center;",
                          div(style="display: flex; align-items: center; gap: 4px;",
                            checkboxInput("include_batch", NULL, value = FALSE),
                            div(style="margin-top: -15px; width: 80px;",
                              textInput("batch_label", NULL, value = "Batch", placeholder = "Batch"))
                          ),
                          div(style="display: flex; align-items: center; gap: 4px;",
                            checkboxInput("include_cov1", NULL, value = FALSE),
                            div(style="margin-top: -15px; width: 110px;",
                              textInput("cov1_label", NULL, value = "Covariate1", placeholder = "e.g., Sex"))
                          ),
                          div(style="display: flex; align-items: center; gap: 4px;",
                            checkboxInput("include_cov2", NULL, value = FALSE),
                            div(style="margin-top: -15px; width: 110px;",
                              textInput("cov2_label", NULL, value = "Covariate2", placeholder = "e.g., Age"))
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
                    # min-height prevents collapse in HF iframe where 100vh can be tiny
                    div(style="overflow-y: auto; max-height: calc(100vh - 420px); min-height: 300px;",
                      rHandsontableOutput("hot_metadata")
                    )
                  ),

                  nav_panel("Signal Distribution",
                    icon = icon("chart-area"),
                    # Comparison selector banner
                    div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 12px 15px; border-radius: 8px; margin-bottom: 15px; display: flex; align-items: center; gap: 10px; flex-wrap: nowrap; position: relative; z-index: 10;",
                      div(style = "display: flex; align-items: center; gap: 10px; white-space: nowrap;",
                        icon("microscope"),
                        span("Viewing Comparison:", style = "font-weight: 500;")
                      ),
                      div(style = "flex: 1 1 auto; min-width: 200px;",
                        selectInput("contrast_selector_signal", NULL,
                          choices = NULL,
                          width = "100%"
                        )
                      )
                    ),
                    # Control buttons
                    div(style = "display: flex; justify-content: space-between; align-items: center; gap: 8px; margin-bottom: 10px;",
                      div(style = "display: flex; align-items: center; gap: 10px;",
                        checkboxInput("signal_overlay_contam", "Overlay Contaminants", value = FALSE)
                      ),
                      div(
                        actionButton("signal_dist_info_btn", icon("question-circle"), title = "What is this?",
                          class = "btn-outline-info btn-sm"),
                        actionButton("fullscreen_signal", "\U0001F50D Fullscreen", class = "btn-outline-secondary btn-sm")
                      )
                    ),
                    plotOutput("protein_signal_plot", height = "calc(100vh - 400px)")
                  ),

                  nav_panel("Dataset Summary",
                    icon = icon("info-circle"),
                    div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                      actionButton("dataset_summary_info_btn", icon("question-circle"),
                        title = "About Dataset Summary", class = "btn-outline-info btn-sm")
                    ),
                    uiOutput("dataset_summary_content")
                  ),

                  nav_panel("Replicate Consistency",
                    icon = icon("arrows-to-circle"),
                    div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                      actionButton("replicate_consistency_info_btn", icon("question-circle"),
                        title = "About Replicate Consistency", class = "btn-outline-info btn-sm"),
                      actionButton("fullscreen_corr_heatmap", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm"),
                      downloadButton("download_replicate_csv", tagList(icon("download"), " CSV"),
                        class = "btn-success btn-sm")
                    ),
                    imageOutput("correlation_heatmap", height = "500px"),
                    div(style = "margin-top: 16px;",
                      tags$h6(icon("table"), " Per-Group Replicate Statistics",
                        style = "font-weight: 600; margin-bottom: 8px;"),
                      DTOutput("replicate_stats_table")
                    )
                  ),

                  nav_panel("Contaminant Analysis",
                    icon = icon("shield-virus"),
                    div(style = "overflow-y: auto; max-height: calc(100vh - 200px);",
                      # Header with info button
                      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                        tags$h5(icon("shield-virus"), " Contaminant Protein Analysis",
                          style = "margin: 0; font-weight: 600;"),
                        actionButton("contaminant_info_btn", icon("question-circle"),
                          title = "About Contaminant Analysis", class = "btn-outline-info btn-sm")
                      ),

                      # Summary statistics cards
                      div(style = "min-height: 80px;",
                        uiOutput("contaminant_summary_cards")
                      ),

                      # Per-sample contaminant breakdown bar chart
                      div(style = "min-height: 400px; margin-top: 15px;",
                        tags$h6(icon("chart-bar"), " Per-Sample Contaminant Intensity",
                          style = "font-weight: 600; margin-bottom: 8px;"),
                        plotlyOutput("contaminant_bar_chart", height = "380px")
                      ),

                      # Top contaminants table
                      div(style = "margin-top: 20px;",
                        tags$h6(icon("table"), " Top Contaminant Proteins",
                          style = "font-weight: 600; margin-bottom: 8px;"),
                        DTOutput("contaminant_top_table")
                      ),

                      # Contaminant heatmap
                      div(style = "min-height: 450px; margin-top: 20px;",
                        tags$h6(icon("fire"), " Contaminant Expression Heatmap (Top 20)",
                          style = "font-weight: 600; margin-bottom: 8px;"),
                        plotlyOutput("contaminant_heatmap", height = "420px")
                      )
                    )
                  ),

                  nav_panel("Expression Grid",
                    icon = icon("th"),
                    # Comparison selector banner
                    div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 12px 15px; border-radius: 8px; margin-bottom: 15px; display: flex; align-items: center; gap: 10px; flex-wrap: nowrap; position: relative; z-index: 10;",
                      div(style = "display: flex; align-items: center; gap: 10px; white-space: nowrap;",
                        icon("microscope"),
                        span("Viewing Comparison:", style = "font-weight: 500;")
                      ),
                      div(style = "flex: 1 1 auto; min-width: 200px;",
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
                    div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                      div(
                        actionButton("grid_reset_selection", "Show All / Clear Selection", class = "btn-warning btn-sm"),
                        downloadButton("download_grid_data", "\U0001F4BE Export Full Table", class = "btn-success btn-sm")
                      ),
                      actionButton("expression_grid_info_btn", icon("question-circle"), title = "What is this?",
                        class = "btn-outline-info btn-sm")
                    ),
                    # Grid table
                    div(style = "overflow-x: auto; width: 100%;",
                      DTOutput("grid_view_table")
                    )
                  ),

                  nav_panel("Data Explorer",
                    icon = icon("search-plus"),
                    div(style = "overflow-y: auto; max-height: calc(100vh - 200px); padding: 15px;",
                      # Header with info button
                      div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
                        tags$h4(icon("search-plus"), " Data Explorer", style = "margin: 0; font-weight: 600;"),
                        actionButton("data_explorer_info_btn", icon("question-circle"),
                          class = "btn-outline-info btn-sm", title = "About Data Explorer"),
                        downloadButton("export_explorer_claude", "Export for Claude",
                          class = "btn-outline-primary btn-sm", icon = icon("download"))
                      ),

                      # --- Panel 1: Abundance Profiles ---
                      tags$h5(icon("layer-group"), " Abundance Profiles (Quartile Analysis)",
                        style = "font-weight: 600; margin-bottom: 10px; padding-top: 5px; border-top: 2px solid #dee2e6;"),
                      div(style = "background-color: #e8f4f8; padding: 12px; border-radius: 8px; margin-bottom: 15px; font-size: 0.9em;",
                        icon("info-circle"), " ",
                        "Proteins split into quartiles by average intensity. Colors show per-sample quartile assignment. ",
                        "Proteins that change quartile across samples may be biologically interesting."
                      ),
                      div(style = "margin-bottom: 15px;",
                        checkboxInput("explorer_exclude_contam_profile", "Exclude contaminants", value = TRUE)
                      ),
                      div(style = "min-height: 500px;",
                        plotlyOutput("explorer_quartile_heatmap", height = "500px")
                      ),
                      tags$h6(icon("exchange-alt"), " Variable Proteins (Quartile Range >= 2)",
                        style = "font-weight: 600; margin-top: 20px; margin-bottom: 10px;"),
                      DTOutput("explorer_variable_proteins_table"),

                      # --- Panel 2: Sample-Sample Scatter ---
                      tags$h5(icon("braille"), " Sample-Sample Scatter",
                        style = "font-weight: 600; margin-top: 30px; margin-bottom: 10px; padding-top: 15px; border-top: 2px solid #dee2e6;"),
                      div(style = "display: flex; gap: 15px; flex-wrap: wrap; align-items: flex-end; margin-bottom: 15px;",
                        div(style = "flex: 1; min-width: 150px;",
                          selectInput("explorer_sample_a", "Sample A", choices = NULL, width = "100%")
                        ),
                        div(style = "flex: 1; min-width: 150px;",
                          selectInput("explorer_sample_b", "Sample B", choices = NULL, width = "100%")
                        ),
                        div(style = "flex: 0 0 auto;",
                          checkboxInput("explorer_label_outliers", "Label outliers (>4-fold)", value = TRUE),
                          checkboxInput("explorer_exclude_contam_scatter", "Exclude contaminants", value = TRUE)
                        )
                      ),
                      div(style = "min-height: 550px;",
                        plotlyOutput("explorer_sample_scatter", height = "550px")
                      )
                    )
                  ),

                  nav_panel("AI Summary",
                    icon = icon("robot"),
                    div(style = "padding: 20px;",
                      # Header section
                      div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 8px; margin-bottom: 20px; display: flex; justify-content: space-between; align-items: center;",
                        tags$h4(icon("robot"), " AI-Powered Analysis Summary", style = "margin: 0; font-weight: 500;"),
                        actionButton("ai_summary_info_btn", icon("question-circle"),
                          class = "btn-outline-light btn-sm", title = "About AI Summary")
                      ),

                      # Instructions
                      div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px;",
                        tags$p(class = "mb-2",
                          icon("info-circle"), " ",
                          strong("How it works:"),
                          " Click the button below to generate a comprehensive AI-powered analysis across all comparisons in your experiment."
                        ),
                        tags$p(class = "mb-0", style = "font-size: 0.9em; color: #6c757d;",
                          "The AI will identify key DE proteins per comparison, cross-comparison biomarkers, ",
                          "and provide biological insights on high-confidence candidates."
                        )
                      ),

                      # Generate button + export + Claude prompt
                      div(style = "text-align: center; margin-bottom: 20px; display: flex; justify-content: center; gap: 12px; flex-wrap: wrap;",
                        actionButton("generate_ai_summary_overview",
                          "\U0001F916 Generate AI Summary",
                          class = "btn-info btn-lg",
                          style = "padding: 12px 30px; font-size: 1.1em;"
                        ),
                        shinyjs::hidden(
                          downloadButton("download_ai_summary_html",
                            tagList(icon("download"), " Download as Markdown"),
                            class = "btn-success btn-lg",
                            style = "padding: 12px 30px; font-size: 1.1em;"
                          )
                        ),
                        downloadButton("download_claude_prompt",
                          tagList(icon("download"), " Export for Claude"),
                          class = "btn-outline-secondary btn-lg",
                          style = "padding: 12px 30px; font-size: 1.1em;"
                        )
                      ),

                      # Output area
                      uiOutput("ai_summary_output")
                    )
                  )
                )
      ),

      nav_spacer(),  # visual divider between Setup and Results

      # -- Results section --
      nav_panel("DE Dashboard", icon = icon("table-columns"),
                # Interactive comparison selector banner
                div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; padding: 15px; border-radius: 8px; margin-bottom: 15px; position: relative; z-index: 10;",
                  div(style = "display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;",
                    div(style = "display: flex; align-items: center; gap: 15px;",
                      icon("microscope"),
                      span("Viewing Comparison:", style = "font-weight: 500;"),
                      selectInput("contrast_selector",
                        label = NULL,
                        choices = NULL,
                        width = "300px"
                      )
                    ),
                    actionButton("de_dashboard_info_btn", icon("question-circle"),
                      title = "How to use this dashboard",
                      class = "btn-outline-light btn-sm")
                  )
                ),

                # Sub-tabs
                navset_card_tab(
                  id = "de_dashboard_subtabs",

                  nav_panel("Volcano", icon = icon("chart-simple"),
                    div(style = "display: grid; grid-template-columns: 1fr 1fr; gap: 16px; align-items: start;",
                      # Left: Volcano plot
                      div(
                        div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                          actionButton("clear_plot_selection_volcano", "Reset Selection", class="btn-warning btn-sm"),
                          actionButton("fullscreen_volcano", "\U0001F50D Fullscreen", class="btn-outline-secondary btn-sm")
                        ),
                        plotlyOutput("volcano_plot_interactive", height = "calc(100vh - 340px)")
                      ),
                      # Right: Heatmap
                      div(
                        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                          span("Heatmap of Selected/Top Proteins", style = "font-weight: 600; font-size: 0.9rem;"),
                          div(style = "display: flex; gap: 8px;",
                            downloadButton("download_heatmap_png", tagList(icon("image"), " PNG"),
                              class = "btn-outline-secondary btn-sm"),
                            downloadButton("download_heatmap_svg", tagList(icon("download"), " SVG"),
                              class = "btn-outline-secondary btn-sm"),
                            actionButton("fullscreen_heatmap", "\U0001F50D Fullscreen", class="btn-outline-secondary btn-sm")
                          )
                        ),
                        plotOutput("heatmap_plot", height = "calc(100vh - 340px)")
                      )
                    )
                  ),

                  nav_panel("Results Table", icon = icon("table"),
                    div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px; flex-wrap: wrap;",
                      actionButton("clear_plot_selection", "Reset", class="btn-warning btn-xs"),
                      actionButton("show_violin", "\U0001F4CA Violin", class="btn-primary btn-xs"),
                      if (!is_hf_space) actionButton("show_xic", "\U0001F4C8 XICs", class="btn-info btn-xs"),
                      downloadButton("download_result_csv", "\U0001F4BE Export", class="btn-success btn-xs")
                    ),
                    DTOutput("de_table")
                  ),

                  nav_panel("PCA", icon = icon("compass"),
                    div(style = "display: flex; justify-content: flex-end; align-items: center; gap: 8px; flex-wrap: wrap; margin-bottom: 10px;",
                      span("Color by:", style = "font-weight: 500; font-size: 0.85em; color: #555;"),
                      div(style = "width: 160px;",
                        selectInput("pca_color_by", label = NULL,
                          choices = c("Group"), selected = "Group", width = "100%")
                      ),
                      span("Axes:", style = "font-weight: 500; font-size: 0.85em; color: #555;"),
                      div(style = "width: 160px;",
                        selectInput("pca_axes", label = NULL,
                          choices = c("PC1 vs PC2" = "1_2", "PC1 vs PC3" = "1_3", "PC2 vs PC3" = "2_3"),
                          width = "100%")
                      ),
                      actionButton("pca_info_btn", icon("question-circle"),
                        title = "About PCA", class = "btn-outline-info btn-sm"),
                      downloadButton("download_pca_png", tagList(icon("image"), " PNG"),
                        class = "btn-outline-secondary btn-sm"),
                      actionButton("fullscreen_pca", "\U0001F50D Fullscreen",
                        class = "btn-outline-secondary btn-sm")
                    ),
                    plotlyOutput("pca_plot", height = "calc(100vh - 370px)")
                  ),

                  nav_panel("CV Analysis", icon = icon("check-double"),
                    div(style = "overflow-y: auto; max-height: calc(100vh - 200px); padding-right: 5px;",
                      # Controls row: info + CSV download
                      div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                        actionButton("consistent_de_info_btn", icon("question-circle"),
                          title = "About CV Analysis", class = "btn-outline-info btn-sm"),
                        downloadButton("download_consistent_csv", tagList(icon("download"), " CSV"),
                          class = "btn-success btn-sm"),
                        actionButton("fullscreen_cv_scatter", "\U0001F50D Fullscreen",
                          class = "btn-outline-secondary btn-sm")
                      ),
                      # logFC vs Avg CV scatter plot
                      div(style = "min-height: 500px;",
                        plotlyOutput("cv_scatter_plot", height = "500px")
                      ),
                      hr(),
                      # CV Distribution histogram
                      div(
                        div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                          p("Distribution of Coefficient of Variation (CV) for significant proteins, broken down by experimental group.",
                            class = "text-muted small mb-0"),
                          div(style = "display: flex; gap: 8px;",
                            actionButton("cv_dist_info_btn", icon("question-circle"), title = "What is this?",
                              class = "btn-outline-info btn-sm"),
                            downloadButton("download_cv_hist_png", tagList(icon("image"), " PNG"),
                              class = "btn-outline-secondary btn-sm"),
                            actionButton("fullscreen_cv_hist", "\U0001F50D Fullscreen",
                              class = "btn-outline-secondary btn-sm")
                          )
                        ),
                        plotOutput("cv_histogram", height = "450px")
                      )
                    )
                  )

                )
      ),

      nav_panel("Phosphoproteomics", icon = icon("flask"),
        uiOutput("phospho_tab_content")
      ),

      nav_panel("Gene Set Enrichment", icon = icon("sitemap"),
                # Contrast indicator
                uiOutput("gsea_contrast_indicator"),
                # Compact control bar
                card(
                  card_body(
                    div(style = "display: flex; align-items: center; gap: 15px; flex-wrap: wrap;",
                      div(style = "min-width: 220px;",
                        selectInput("gsea_ontology", NULL,
                          choices = c(
                            "GO: Biological Process (BP)" = "BP",
                            "GO: Molecular Function (MF)" = "MF",
                            "GO: Cellular Component (CC)" = "CC",
                            "KEGG Pathways" = "KEGG"
                          ),
                          selected = "BP", width = "100%"
                        )
                      ),
                      actionButton("run_gsea", "\u25B6 Run GSEA", class = "btn-success", icon = icon("play")),
                      div(style = "flex-grow: 1;",
                        verbatimTextOutput("gsea_status", placeholder = TRUE) |>
                          tagAppendAttributes(style = "margin: 0; padding: 5px 10px; min-height: 38px;")
                      ),
                      actionButton("gsea_info_btn", icon("question-circle"), title = "What is GSEA?",
                        class = "btn-outline-info btn-sm")
                    ),
                    p("Enrichment analysis on DE results. Auto-detects organism. Results cached per ontology.",
                      class = "text-muted small", style = "margin: 10px 0 0 0;")
                  )
                ),

                # Results tabs with full-height plots
                navset_card_tab(
                  id = "gsea_results_tabs",

                  nav_panel("Dot Plot",
                    div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                      downloadButton("download_gsea_dot_png", tagList(icon("image"), " PNG"),
                        class = "btn-outline-secondary btn-sm"),
                      actionButton("fullscreen_gsea_dot", "\U0001F50D Fullscreen", class = "btn-outline-secondary btn-sm")
                    ),
                    plotOutput("gsea_dot_plot", height = "calc(100vh - 340px)")
                  ),

                  nav_panel("Enrichment Map",
                    div(style = "text-align: right; margin-bottom: 10px;",
                      actionButton("fullscreen_gsea_emap", "\U0001F50D Fullscreen", class = "btn-outline-secondary btn-sm")
                    ),
                    plotOutput("gsea_emapplot", height = "calc(100vh - 340px)")
                  ),

                  nav_panel("Ridgeplot",
                    div(style = "text-align: right; margin-bottom: 10px;",
                      actionButton("fullscreen_gsea_ridge", "\U0001F50D Fullscreen", class = "btn-outline-secondary btn-sm")
                    ),
                    plotOutput("gsea_ridgeplot", height = "calc(100vh - 340px)")
                  ),

                  nav_panel("Results Table",
                    div(style = "display: flex; justify-content: flex-end; gap: 8px; margin-bottom: 10px;",
                      actionButton("gsea_table_info_btn", icon("question-circle"), title = "Column definitions",
                        class = "btn-outline-info btn-sm"),
                      downloadButton("download_gsea_csv", tagList(icon("download"), " CSV"),
                        class = "btn-success btn-sm")
                    ),
                    DTOutput("gsea_results_table")
                  )
                )
      ),

      nav_panel("Multi-Omics MOFA2", icon = icon("layer-group"),
                value = "mofa_tab",
                uiOutput("mofa_tab_content")
      ),

      # -- Run Comparator --
      nav_panel("Run Comparator", icon = icon("code-compare"),
                value = "comparator_tab",
        div(style = "overflow-y: auto; max-height: calc(100vh - 120px); padding: 15px;",
          # Configure Comparison card
          div(style = "border: 1px solid #dee2e6; border-radius: 8px; padding: 15px; margin-bottom: 15px; background: #f8f9fa;",
            div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
              tags$h5("Configure Comparison", style = "margin: 0;"),
              actionButton("comparator_info_btn", icon("question-circle"),
                           class = "btn-outline-info btn-sm", title = "About Run Comparator")
            ),
            div(style = "display: flex; gap: 15px; flex-wrap: wrap;",
              # Column 1: Mode
              div(style = "flex: 1; min-width: 200px;",
                radioButtons("comparator_mode", "Comparison Type",
                  choices = c(
                    "DE-LIMP vs DE-LIMP"     = "delimp_delimp",
                    "DE-LIMP vs Spectronaut" = "delimp_spectronaut",
                    "DE-LIMP vs FragPipe"    = "delimp_fragpipe"
                  )
                )
              ),
              # Column 2: Run A
              div(style = "flex: 1; min-width: 200px;",
                radioButtons("comparator_run_a_source", "Run A Source",
                  choices = c("Current session" = "current", "Load from file" = "file")),
                conditionalPanel("input.comparator_run_a_source == 'file'",
                  fileInput("comparator_run_a_file", "Load Run A (.rds)", accept = ".rds")
                )
              ),
              # Column 3: Run B (conditional on mode)
              div(style = "flex: 1; min-width: 200px;",
                conditionalPanel("input.comparator_mode == 'delimp_delimp'",
                  fileInput("comparator_run_b_rds", "Load Run B (.rds)", accept = ".rds")
                ),
                conditionalPanel("input.comparator_mode == 'delimp_spectronaut'",
                  fileInput("comparator_run_b_spec_zip", "Spectronaut Export (.zip)",
                            accept = ".zip"),
                  uiOutput("spectronaut_zip_manifest"),
                  tags$details(style = "margin-top: 4px; margin-bottom: 8px;",
                    tags$summary(tags$small(class = "text-muted", "Or upload individual files...")),
                    div(style = "padding-top: 6px;",
                      fileInput("comparator_run_b_spectronaut", "Protein Quantities (.tsv)",
                                accept = c(".tsv", ".csv", ".txt")),
                      fileInput("comparator_run_b_spectronaut_de", "Candidates / DE Stats (.tsv, optional)",
                                accept = c(".tsv", ".csv", ".txt"))
                    )
                  ),
                  div(style = "display: flex; gap: 6px; align-items: center; margin-top: -8px;",
                    downloadButton("download_spectronaut_schema",
                      "Spectronaut Setup Guide", class = "btn-outline-info btn-xs"),
                    tags$small(class = "text-muted", "How to export from Spectronaut")
                  )
                ),
                conditionalPanel("input.comparator_mode == 'delimp_fragpipe'",
                  radioButtons("comparator_fragpipe_type", "FragPipe output type",
                    choices = c(
                      "FragPipe-Analyst DE export" = "fp_analyst",
                      "combined_protein.tsv (intensities only)" = "fp_raw"
                    )
                  ),
                  conditionalPanel("input.comparator_fragpipe_type == 'fp_analyst'",
                    fileInput("comparator_fp_analyst_file",
                              "FragPipe-Analyst DE results (.csv/.tsv)",
                              accept = c(".csv", ".tsv"))
                  ),
                  conditionalPanel("input.comparator_fragpipe_type == 'fp_raw'",
                    fileInput("comparator_fp_combined_protein", "combined_protein.tsv",
                              accept = ".tsv"),
                    div(class = "alert alert-info py-1 mt-1", style = "font-size: 0.8em;",
                      icon("info-circle"), " Layers 1-3 only. ",
                      "Run FragPipe-Analyst for full DE comparison."
                    )
                  )
                )
              )
            ),
            # Optional DIA-NN log upload (Mode A only)
            conditionalPanel("input.comparator_mode == 'delimp_delimp'",
              div(style = "margin-top: 8px; border-top: 1px solid #dee2e6; padding-top: 8px;",
                actionLink("toggle_diann_logs",
                  tagList(icon("chevron-right", id = "diann_log_chevron"),
                          " Attach DIA-NN log files (optional \u2014 fills in search parameters)")),
                conditionalPanel("input.toggle_diann_logs % 2 == 1",
                  div(class = "mt-2",
                    div(class = "alert alert-info py-1 px-2", style = "font-size: 0.82em;",
                      icon("info-circle"), " Upload ",
                      tags$code("report_log.txt"), " or the SLURM ",
                      tags$code(".out"), " file from each DIA-NN search. ",
                      "Only the command line and summary stats are read."
                    ),
                    div(style = "display: flex; gap: 12px; flex-wrap: wrap;",
                      div(style = "flex: 1; min-width: 200px;",
                        fileInput("comparator_diann_log_a", "Run A \u2014 DIA-NN log",
                                  accept = c(".txt", ".log", ".out"))
                      ),
                      div(style = "flex: 1; min-width: 200px;",
                        fileInput("comparator_diann_log_b", "Run B \u2014 DIA-NN log",
                                  accept = c(".txt", ".log", ".out"))
                      )
                    ),
                    div(id = "comparator_diann_log_status")
                  )
                )
              )
            ),
            # Contrast selectors + sample status (dynamic)
            uiOutput("comparator_contrast_selectors"),
            uiOutput("comparator_sample_status"),
            div(style = "margin-top: 10px;",
              actionButton("run_comparison", "Run Comparison",
                           class = "btn-primary", icon = icon("play"))
            )
          ),

          # Results (shown after comparison runs)
          conditionalPanel("input.run_comparison > 0",
            uiOutput("comparator_summary_banner"),
            navset_card_tab(id = "comparator_subtabs",
              nav_panel("Settings Diff", icon = icon("sliders"),
                div(style = "overflow-y: auto; max-height: calc(100vh - 200px); min-height: 300px;",
                  div(style = "display: flex; justify-content: flex-end; margin-bottom: 4px;",
                    actionButton("comparator_settings_info_btn", icon("question-circle"),
                      class = "btn-outline-info btn-sm", title = "About Settings Diff")
                  ),
                  div(id = "comparator_pipeline_warning"),
                  DT::DTOutput("comparator_settings_diff"),
                  # Per-Sample QC (Mode B only — shown when RunSummaries present)
                  uiOutput("comparator_sample_qc_section")
                )
              ),
              nav_panel("Protein Universe", icon = icon("circle-nodes"),
                div(style = "overflow-y: auto; max-height: calc(100vh - 200px);",
                  div(style = "display: flex; justify-content: flex-end; margin-bottom: 4px;",
                    actionButton("comparator_universe_info_btn", icon("question-circle"),
                      class = "btn-outline-info btn-sm", title = "About Protein Universe")
                  ),
                  div(style = "display: flex; gap: 15px; flex-wrap: wrap; align-items: start;",
                    div(style = "flex: 1; min-width: 350px;",
                      plotly::plotlyOutput("comparator_universe_plot", height = "300px")
                    ),
                    div(style = "flex: 1; min-width: 350px;",
                      uiOutput("comparator_universe_summary")
                    )
                  ),
                  div(style = "display: flex; align-items: center; gap: 12px; margin-top: 12px;",
                    tags$h6("Protein Details", style = "margin: 0;"),
                    div(style = "display: flex; gap: 6px;",
                      actionButton("universe_filter_all", "All",
                        class = "btn-outline-secondary btn-sm active"),
                      actionButton("universe_filter_shared", "Shared",
                        class = "btn-outline-success btn-sm"),
                      actionButton("universe_filter_a_only", "Run A only",
                        class = "btn-outline-primary btn-sm"),
                      actionButton("universe_filter_b_only", "Run B only",
                        class = "btn-outline-warning btn-sm")
                    ),
                    div(style = "margin-left: auto;",
                      downloadButton("download_universe_csv", "Export CSV",
                        class = "btn-outline-secondary btn-sm")
                    )
                  ),
                  DT::DTOutput("comparator_universe_table")
                )
              ),
              nav_panel("Quantification", icon = icon("chart-line"),
                div(style = "overflow-y: auto; max-height: calc(100vh - 200px);",
                  div(style = "display: flex; justify-content: flex-end; margin-bottom: 4px;",
                    actionButton("comparator_quant_info_btn", icon("question-circle"),
                      class = "btn-outline-info btn-sm", title = "About Quantification")
                  ),
                  div(style = "display: flex; gap: 15px; flex-wrap: wrap; min-height: 350px;",
                    div(style = "flex: 1; min-width: 400px;",
                      plotly::plotlyOutput("comparator_quant_scatter", height = "380px")
                    ),
                    div(style = "flex: 1; min-width: 300px;",
                      plotly::plotlyOutput("comparator_correlation_heatmap", height = "380px")
                    )
                  ),
                  plotly::plotlyOutput("comparator_bias_density", height = "250px"),
                  # TopN Effect scatter (Mode B only)
                  uiOutput("comparator_topn_effect_section")
                )
              ),
              nav_panel("DE Concordance", icon = icon("code-compare"),
                div(
                  div(style = "display: flex; justify-content: flex-end; margin-bottom: 4px;",
                    actionButton("comparator_concordance_info_btn", icon("question-circle"),
                      class = "btn-outline-info btn-sm", title = "About DE Concordance")
                  ),
                  uiOutput("comparator_layer4_content")
                )
              ),
              nav_panel("AI Analysis", icon = icon("robot"), value = "comparator_ai_tab",
                tags$div(style = "padding: 12px 4px; min-height: 400px;",
                  div(style = "display: flex; justify-content: space-between; align-items: center;",
                    tags$h6("AI-Powered Comparison Analysis"),
                    actionButton("comparator_ai_info_btn", icon("question-circle"),
                      class = "btn-outline-info btn-sm", title = "About AI Analysis")
                  ),
                  tags$p(class = "text-muted small",
                    "Generate an AI narrative summary or export data for external analysis."),
                  tags$div(style = "display: flex; gap: 10px; align-items: center; flex-wrap: wrap; margin-bottom: 16px;",
                    actionButton("comparator_gemini_btn", "Generate Gemini Summary",
                                 icon = icon("wand-magic-sparkles"),
                                 class = "btn-outline-primary"),
                    actionButton("comparator_view_prompt_btn", "View Prompt",
                                 icon = icon("eye"),
                                 class = "btn-outline-info btn-sm"),
                    downloadButton("comparator_claude_export",
                                   "Export ZIP for Claude Analysis",
                                   icon = icon("file-zipper"),
                                   class = "btn-outline-secondary")
                  ),
                  tags$div(id = "comparator_gemini_container"),
                  tags$hr(),
                  tags$h6("MOFA2 Factor Decomposition", class = "text-muted"),
                  tags$p(class = "text-muted small",
                    "Treats Run A and Run B as two views of the same samples and decomposes ",
                    "joint variance into shared and run-specific factors."),
                  actionButton("comparator_mofa_btn",
                    "Run MOFA2 Decomposition (~1-2 min)",
                    icon = icon("circle-nodes"),
                    class = "btn-outline-secondary"),
                  tags$div(style = "margin-top: 12px;",
                    plotly::plotlyOutput("comparator_mofa_variance", height = "320px"),
                    plotly::plotlyOutput("comparator_mofa_weights", height = "380px"),
                    DT::DTOutput("comparator_mofa_top_weights")
                  )
                )
              )
            )
          )
        )
      ),

      nav_spacer(),  # visual divider before AI section

      # -- AI section --
      nav_panel("AI Analysis", icon = icon("robot"),
                card(
                  card_header(div(style="display: flex; justify-content: space-between; align-items: center;",
                    span("Chat with Full Data (QC + Expression)"),
                    div(style = "display: flex; gap: 8px;",
                      actionButton("data_chat_info_btn", icon("question-circle"), title = "About Data Chat",
                        class = "btn-outline-info btn-sm"),
                      downloadButton("download_chat_txt", "\U0001F4BE Save Chat", class="btn-secondary btn-sm"),
                      downloadButton("download_claude_prompt_chat", tagList(icon("download"), " Export for Claude"),
                        class = "btn-outline-secondary btn-sm")
                    )
                  )),
                  card_body(
                    verbatimTextOutput("chat_selection_indicator"),
                    uiOutput("chat_window"),
                    tags$div(style="margin-top: 15px; display: flex; gap: 10px;",
                             textAreaInput("chat_input", label=NULL, placeholder="Ask e.g. 'Which group has higher precursor counts?'", width="100%", rows=2),
                             actionButton("summarize_data", "\U0001F916 Auto-Analyze", class="btn-info", style="height: 54px; margin-top: 2px;"),
                             actionButton("send_chat", "Send", icon=icon("paper-plane"), class="btn-primary", style="height: 54px; margin-top: 2px;")
                    ),
                    p("Note: QC Stats (with Groups) + Top 800 Expression Data are sent to AI.", style="font-size: 0.8em; color: green; font-weight: bold; margin-top: 5px;")
                  )
                )
      )
    ),

    # ==========================================================================
    # DE NOVO dropdown — Cascadia de novo sequencing
    # ==========================================================================
    nav_menu("De Novo", icon = icon("dna"),
      nav_panel("Cascadia", value = "cascadia_tab", icon = icon("dna"),
        div(style = "overflow-y: auto; max-height: calc(100vh - 200px);",
          uiOutput("denovo_summary_cards"),
          navset_card_tab(
            nav_panel("Confirmed Peptides",
              DT::DTOutput("denovo_confirmed_table")),
            nav_panel("Novel Peptides",
              DT::DTOutput("denovo_novel_table")),
            nav_panel("DIAMOND BLAST",
              div(style = "margin-bottom: 15px;",
                actionButton("run_diamond_blast", "Run DIAMOND BLAST",
                  icon = icon("search"), class = "btn-info btn-sm")),
              DT::DTOutput("denovo_blast_table")),
            nav_panel("Score Distribution",
              plotlyOutput("denovo_score_dist", height = "400px"))
          )
        )
      ),
      nav_panel("Submit Cascadia Job", value = "cascadia_submit_tab", icon = icon("rocket"),
        div(style = "max-width: 800px; margin: 0 auto; padding: 20px;",
          tags$h4("Submit Cascadia De Novo Sequencing"),
          tags$p("Run Cascadia on the same raw files alongside DIA-NN for orthogonal validation."),
          uiOutput("cascadia_submit_ui")
        )
      )
    ),

    # ==========================================================================
    # OUTPUT dropdown — Methods & Code
    # ==========================================================================
    nav_menu("Output", icon = icon("file-export"),
      nav_panel("Export Data", icon = icon("download"),
        div(style = "max-width: 800px; margin: 30px auto; padding: 20px;",

          # --- Complete Analysis ZIP ---
          div(style = "background-color: #f0f7ff; padding: 20px; border-radius: 8px; margin-bottom: 20px; border-left: 4px solid #6f42c1;",
            tags$h4(icon("file-archive"), " Export Complete Analysis"),
            tags$p(class = "text-muted",
              "Download everything needed to reproduce and share this analysis. ",
              "Includes all data files, DIA-NN search parameters, and session state."
            ),
            tags$details(
              tags$summary(style = "cursor: pointer; color: #6f42c1; font-weight: 500;",
                "What's included (click to expand)"),
              tags$ul(style = "font-size: 0.88em; color: #555; margin-top: 8px;",
                tags$li(tags$strong("expression_matrix.csv"), " -- Normalized protein intensities (DPC-Quant, complete, no missing values)"),
                tags$li(tags$strong("diann_pg_matrix.tsv"), " -- DIA-NN protein-level matrix with real missing values (0 = not detected, ~200 KB)"),
                tags$li(tags$strong("data_quality_summary.csv"), " -- Per-sample protein counts, % detected, contaminant counts"),
                tags$li(tags$strong("detection_matrix.csv"), " -- Per-protein precursor detection counts per sample"),
                tags$li(tags$strong("quartile_profiles.csv"), " -- Intensity quartile assignments per sample"),
                tags$li(tags$strong("variable_proteins.csv"), " -- Proteins with inconsistent abundance across samples"),
                tags$li(tags$strong("sample_metadata.csv"), " -- Sample groups and identifiers"),
                tags$li(tags$strong("contaminant_summary.csv"), " -- Contaminant protein statistics"),
                tags$li(tags$strong("search_info.md"), " -- Full DIA-NN search parameters and job metadata"),
                tags$li(tags$strong("session.rds"), " -- Complete session state (reload in DE-LIMP)"),
                tags$li(tags$strong("methods.txt"), " -- Pipeline parameters, normalization, app version"),
                tags$li(tags$strong("reproducibility_log.R"), " -- R code log to reproduce every step"),
                tags$li(tags$strong("PROMPT.md"), " -- AI analysis prompt with biological questions")
              )
            ),
            downloadButton("export_complete_analysis", tagList(icon("download"), " Export Complete Analysis ZIP"),
              class = "btn-primary btn-lg mt-2")
          ),

          # --- DE Results CSV ---
          div(style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px; border-left: 4px solid #28a745;",
            tags$h4(icon("file-csv"), " DE Results Table"),
            tags$p(class = "text-muted",
              "Quick export of the DE results for the selected comparison. ",
              "Includes gene symbols, logFC, P.Value, adj.P.Val, and per-sample expression values. ",
              "One CSV file — no search parameters or session data."
            ),
            downloadButton("download_result_csv_output", tagList(icon("download"), " Export Results CSV"),
              class = "btn-success mt-2")
          ),

          # --- CV Analysis CSV ---
          div(style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px; border-left: 4px solid #17a2b8;",
            tags$h4(icon("chart-bar"), " CV Analysis"),
            tags$p(class = "text-muted",
              "Coefficient of variation for significant proteins. ",
              "Includes per-group CV and average CV values. One CSV file."
            ),
            downloadButton("download_consistent_csv_output", tagList(icon("download"), " Export CV Analysis CSV"),
              class = "btn-info mt-2")
          ),

          # --- DIA-NN Output Location ---
          div(style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px; border-left: 4px solid #6c757d;",
            tags$h4(icon("server"), " Full DIA-NN Output"),
            tags$p(class = "text-muted",
              "The complete DIA-NN search output (report.parquet, precursor matrices, spectral libraries, logs) ",
              "is stored on the HPC cluster. These files can be large (100 MB+) and are not included in the analysis export."
            ),
            uiOutput("diann_output_path_display")
          )
        )
      ),
      nav_panel("Methods & Code", icon = icon("scroll"),
                navset_card_tab(
                  nav_panel("R Code Log",
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
                              downloadButton("download_repro_log", "\U0001F4BE Download Reproducibility Log", class="btn-success mb-3"),
                              verbatimTextOutput("reproducible_code")
                            )
                  ),
                  nav_panel("Methods Summary",
                            card_body(
                              div(style = "display: flex; justify-content: flex-end; margin-bottom: 10px;",
                                actionButton("methodology_info_btn", icon("question-circle"), title = "About the methods",
                                  class = "btn-outline-info btn-sm")
                              ),
                              div(style = "overflow: auto; max-height: calc(100vh - 200px);",
                                verbatimTextOutput("methodology_text")
                              )
                            )
                  )
                )
      )
    ),

    # ==========================================================================
    # ABOUT dropdown — Community stats + Analysis History as separate tabs
    # ==========================================================================
    nav_menu("About", icon = icon("circle-info"),
      nav_panel("Community", value = "about_tab", icon = icon("chart-line"),
        div(style = "padding: 20px; max-width: 900px; margin: 0 auto;",
          # Header with version
          div(style = "text-align: center; margin-bottom: 25px;",
            tags$h3("DE-LIMP"),
            tags$p(class = "text-muted", textOutput("about_version_text", inline = TRUE)),
            tags$p(style = "color: #718096; font-size: 0.9em;",
              "Differential Expression \u2014 LIMPA Pipeline")
          ),

          # Stats cards row
          uiOutput("community_stats_cards"),

          # Trend sparklines row
          uiOutput("community_trend_plots"),

          # Recent discussions
          uiOutput("community_discussions"),

          # Links section
          div(style = "text-align: center; margin-top: 25px;",
            tags$h6("Links"),
            tags$a(href = "https://github.com/bsphinney/DE-LIMP", target = "_blank",
              icon("github"), " GitHub"), " | ",
            tags$a(href = "https://huggingface.co/spaces/brettsp/de-limp-proteomics", target = "_blank",
              icon("rocket"), " Hugging Face"), " | ",
            tags$a(href = "https://bsphinney.github.io/DE-LIMP/", target = "_blank",
              icon("book"), " Documentation"), " | ",
            tags$a(href = "https://github.com/bsphinney/DE-LIMP/discussions", target = "_blank",
              icon("comments"), " Discussions")
          ),

          # Stats freshness note
          uiOutput("stats_updated_at")
        )
      ),
      if (!is_hf_space) nav_panel("History", value = "history_tab", icon = icon("clock-rotate-left"),
        div(style = "padding: 20px;",
          div(style = "display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 8px; margin-bottom: 10px;",
            h4("History", style = "margin: 0;"),
            div(style = "display: flex; align-items: center; gap: 6px;",
              div(style = "width: 180px;",
                selectizeInput("project_filter", NULL, choices = NULL,
                  options = list(placeholder = "Filter by project...", allowEmptyOption = TRUE))
              ),
              div(style = "width: 140px;",
                selectInput("history_status_filter", NULL,
                  choices = c("All statuses" = "", "completed", "submitted", "running", "failed"),
                  selected = "")
              ),
              div(style = "width: 140px;",
                selectizeInput("history_user_filter", NULL, choices = NULL,
                  options = list(placeholder = "Filter by user...", allowEmptyOption = TRUE))
              ),
              actionButton("history_refresh_btn", "Refresh",
                icon = icon("arrows-rotate"), class = "btn-outline-primary btn-sm"),
              downloadButton("history_export_csv", "CSV",
                class = "btn-outline-secondary btn-sm"),
              tags$button(id = "history_compare_btn", type = "button",
                class = "btn btn-outline-warning btn-sm action-button shiny-bound-input",
                style = "display:none;",
                onclick = "var cbs=$('.history-compare-cb:checked');if(cbs.length===2){Shiny.setInputValue('history_compare_click',{od_a:$(cbs[0]).data('od'),od_b:$(cbs[1]).data('od'),sf_a:$(cbs[0]).data('sf'),sf_b:$(cbs[1]).data('sf'),name_a:$(cbs[0]).data('name'),name_b:$(cbs[1]).data('name'),ts:Date.now()});}",
                icon("code-compare"), " Compare")
            )
          ),
          uiOutput("history_source_badge"),
          p(class = "text-muted", "Searches and analyses from this machine/volume. Click a row to expand details. Green ", tags$b("Load"), " = full session restore (post-pipeline). Outline ", tags$b("Raw"), " = load report.parquet only. Check two analyses to compare."),
          uiOutput("project_summary_cards"),
          DTOutput("history_table")
        )
      )
    ),

    # ==========================================================================
    # EDUCATION (standalone)
    # ==========================================================================
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
    ),

    # ==========================================================================
    # FACILITY dropdown (conditional — core facility mode only)
    # ==========================================================================
    if (is_core_facility) nav_menu("Facility", icon = icon("building"),

      # ---------- Search DB ----------
      nav_panel("Search DB", icon = icon("database"),
        div(style = "padding: 15px;",
          # Header with count
          div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
            tags$h5(icon("database"), " Search Database", style = "margin: 0;"),
            tags$span(class = "text-muted", textOutput("job_count_text", inline = TRUE))
          ),

          # Filter row
          div(style = "display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 12px;",
            div(style = "flex: 2; min-width: 180px;",
              textInput("job_search_text", NULL, placeholder = "Search by name...")
            ),
            div(style = "flex: 1; min-width: 120px;",
              selectInput("job_filter_lab", NULL,
                choices = c("All labs" = ""), selected = "")
            ),
            div(style = "flex: 1; min-width: 120px;",
              selectInput("job_filter_status", NULL,
                choices = c("All" = "", "Queued" = "queued", "Running" = "running",
                            "Completed" = "completed", "Failed" = "failed"),
                selected = "")
            ),
            div(style = "flex: 1; min-width: 120px;",
              selectInput("job_filter_staff", NULL,
                choices = c("All staff" = ""), selected = "")
            ),
            div(style = "flex: 1; min-width: 120px;",
              selectInput("job_filter_instrument", NULL,
                choices = c("All instruments" = ""), selected = "")
            ),
            div(style = "flex: 1; min-width: 120px;",
              selectInput("job_filter_lc_method", NULL,
                choices = c("All LC methods" = ""), selected = "")
            )
          ),

          # Job history table
          DTOutput("job_history_table"),

          # Action buttons
          div(style = "margin-top: 10px; display: flex; gap: 8px; align-items: center;",
            actionButton("job_load_results", "Load Selected Results",
              icon = icon("folder-open"),
              class = "btn-outline-primary btn-sm"),
            actionButton("job_generate_report", "Generate Report",
              icon = icon("file-export"),
              class = "btn-outline-success btn-sm"),
            div(style = "flex: 1;"),
            uiOutput("report_link_ui")
          )
        )
      ),

      # ---------- Instrument QC ----------
      nav_panel("Instrument QC", icon = icon("heartbeat"),
        div(style = "padding: 15px;",
          div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;",
            tags$h5("Instrument Performance Dashboard", style = "margin: 0;"),
            div(style = "display: flex; gap: 10px; align-items: center;",
              actionButton("qc_ingest_btn", "Ingest QC Run",
                icon = icon("plus-circle"),
                class = "btn-outline-success btn-sm"),
              selectInput("qc_instrument_filter", NULL,
                choices = c("All Instruments" = ""),
                selected = "", width = "200px"),
              selectInput("qc_date_range", NULL,
                choices = c("Last 30 days" = "30", "Last 90 days" = "90",
                            "Last 180 days" = "180", "All time" = "9999"),
                selected = "90", width = "150px")
            )
          ),
          plotly::plotlyOutput("qc_protein_trend", height = "280px"),
          plotly::plotlyOutput("qc_precursor_trend", height = "280px"),
          plotly::plotlyOutput("qc_tic_trend", height = "280px"),
          hr(),
          div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
            tags$h6("QC Runs", style = "margin: 0;"),
            div(style = "display: flex; gap: 8px;",
              actionButton("qc_exclude_btn", "Exclude Selected",
                icon = icon("ban"), class = "btn-outline-warning btn-sm"),
              actionButton("qc_include_btn", "Re-include Selected",
                icon = icon("undo"), class = "btn-outline-info btn-sm")
            )
          ),
          DTOutput("qc_runs_table")
        )
      )
    ),

    # Gear icon pushed to far-right of navbar
    nav_spacer(),
    nav_item(actionLink("open_settings", label = NULL, icon = icon("gear"), title = "Settings"))
  )
}
