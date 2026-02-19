# ==============================================================================
#  USER INTERFACE (UI) — Build the complete app UI
#  Called from app.R as: ui <- build_ui(is_hf_space, hpc_mode, local_sbatch)
# ==============================================================================

build_ui <- function(is_hf_space, hpc_mode = FALSE, local_sbatch = FALSE) {
  page_sidebar(
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
    actionButton("load_example", "\U0001F4CA Load Example Data", class = "btn-info btn-sm w-100",
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
    # Phospho controls (conditional on detection)
    conditionalPanel(
      condition = "output.phospho_detected_flag",
      hr(),
      h5("5. Phosphoproteomics", style = "color: #6c757d;"),
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
      actionButton("load_example_phospho", "Load Example Phospho Data",
        class = "btn-outline-info btn-sm w-100", icon = icon("flask"),
        style = "margin-bottom: 10px;"),
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
    ),

    if (!is_hf_space) tagList(
      h5("6. XIC Viewer"),
      p(class = "text-muted small",
        "Load .xic.parquet files from DIA-NN to inspect chromatograms."),
      textInput("xic_dir_input", "XIC Directory Path:",
        placeholder = "Auto-detected or paste path here"),
      actionButton("xic_load_dir", "Load XICs", class = "btn-outline-info btn-sm w-100",
        icon = icon("wave-square")),
      uiOutput("xic_status_badge")
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

  navset_card_tab(
    id = "main_tabs",

    # ==========================================================================
    # New Search tab (HPC only — visible when sbatch is on PATH)
    # ==========================================================================
    if (hpc_mode) nav_panel("New Search", icon = icon("rocket"),
      # Three-panel wizard layout
      layout_column_wrap(
        width = 1/3,

        # === PANEL 1: FILES ===
        card(
          card_header(tagList(icon("folder-open"), " 1. Files")),
          card_body(
            style = "overflow-y: auto; max-height: calc(100vh - 200px);",

            textInput("analysis_name", "Analysis Name",
              value = paste0("search_", format(Sys.Date(), "%Y%m%d")),
              placeholder = "e.g., HeLa_DIA_2026"),

            hr(),
            tags$h6(icon("hard-drive"), " Raw Data"),
            conditionalPanel("input.search_connection_mode != 'ssh'",
              shinyFiles::shinyDirButton("raw_data_dir", "Select Raw Data Folder",
                title = "Choose directory with .d / .raw / .mzML files",
                class = "btn-outline-primary btn-sm w-100")
            ),
            conditionalPanel("input.search_connection_mode == 'ssh'",
              div(style = "display: flex; gap: 5px;",
                div(style = "flex: 1;",
                  textInput("ssh_raw_data_dir", NULL,
                    placeholder = "/share/proteomics/raw/experiment_001")
                ),
                actionButton("ssh_scan_raw_btn", "Scan", icon = icon("magnifying-glass"),
                  class = "btn-outline-secondary btn-sm", style = "margin-top: 0;")
              )
            ),
            uiOutput("raw_file_summary"),

            hr(),
            tags$h6(icon("dna"), " FASTA Database"),
            selectInput("fasta_source", NULL,
              choices = c("Download from UniProt" = "uniprot",
                          "Pre-staged on server" = "prestaged",
                          "Browse / enter path" = "browse"),
              width = "100%"),

            # --- UniProt source ---
            conditionalPanel("input.fasta_source == 'uniprot'",
              div(style = "display: flex; gap: 5px; margin-bottom: 8px;",
                div(style = "flex: 1;",
                  textInput("uniprot_search_query", NULL,
                    placeholder = "e.g., human, mouse, E. coli")
                ),
                actionButton("search_uniprot", "Search",
                  class = "btn-info btn-sm", style = "margin-top: 0;")
              ),
              DTOutput("uniprot_results_table", height = "150px"),
              selectInput("fasta_content_type", "Content:",
                choices = c(
                  "One per gene (recommended)" = "one_per_gene",
                  "Swiss-Prot reviewed" = "reviewed",
                  "Full proteome" = "full",
                  "Full + isoforms" = "full_isoforms"
                ), selected = "one_per_gene", width = "100%"
              ),
              checkboxInput("append_contaminants", "Append contaminants", TRUE),
              uiOutput("fasta_filename_preview"),
              actionButton("download_fasta_btn", "Download FASTA",
                class = "btn-success btn-sm w-100", icon = icon("download"))
            ),

            # --- Pre-staged source ---
            conditionalPanel("input.fasta_source == 'prestaged'",
              selectInput("prestaged_fasta", "Available Databases:",
                choices = NULL, width = "100%"),
              uiOutput("prestaged_fasta_info")
            ),

            # --- Browse / path source ---
            conditionalPanel("input.fasta_source == 'browse'",
              conditionalPanel("input.search_connection_mode != 'ssh'",
                shinyFiles::shinyDirButton("fasta_browse_dir", "Browse for FASTA Folder",
                  title = "Navigate to FASTA directory",
                  class = "btn-outline-primary btn-sm w-100")
              ),
              conditionalPanel("input.search_connection_mode == 'ssh'",
                div(style = "display: flex; gap: 5px;",
                  div(style = "flex: 1;",
                    textInput("ssh_fasta_browse_dir", NULL,
                      placeholder = "/share/proteomics/fasta/")
                  ),
                  actionButton("ssh_scan_fasta_btn", "Scan", icon = icon("magnifying-glass"),
                    class = "btn-outline-secondary btn-sm", style = "margin-top: 0;")
                )
              ),
              uiOutput("browsed_fasta_info")
            ),

            hr(),
            tags$h6(icon("book"), " Spectral Library (optional)"),
            conditionalPanel("input.search_connection_mode != 'ssh'",
              shinyFiles::shinyFilesButton("lib_file", "Select .speclib File",
                title = "Choose spectral library",
                class = "btn-outline-secondary btn-sm w-100",
                multiple = FALSE)
            ),
            conditionalPanel("input.search_connection_mode == 'ssh'",
              textInput("ssh_lib_file", NULL,
                placeholder = "/share/proteomics/libraries/my.speclib")
            ),
            uiOutput("lib_file_info")
          )
        ),

        # === PANEL 2: SEARCH SETTINGS ===
        card(
          card_header(tagList(icon("sliders"), " 2. Search Settings")),
          card_body(
            style = "overflow-y: auto; max-height: calc(100vh - 200px);",

            radioButtons("search_mode", "Search Mode:",
              choices = c(
                "Library-free (default)" = "libfree",
                "Use spectral library" = "library"
              ), selected = "libfree"
            ),

            radioButtons("diann_normalization", "DIA-NN Normalization:",
              choices = c(
                "RT-dependent (default)" = "on",
                "Off (for AP-MS / Co-IP)" = "off"
              ), selected = "on"
            ),
            uiOutput("norm_guidance_search"),

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
              condition = "input.mass_acc_mode == 'manual'",
              div(style = "display: flex; gap: 8px;",
                div(style = "flex: 1;",
                  numericInput("diann_mass_acc", "MS2 (ppm):", value = 14, min = 1, max = 50)
                ),
                div(style = "flex: 1;",
                  numericInput("diann_mass_acc_ms1", "MS1 (ppm):", value = 14, min = 1, max = 50)
                )
              )
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
                    numericInput("max_pr_mz", "Max m/z:", value = 1200, min = 800, max = 2000)
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

            tags$h6(icon("plug"), " Connection Mode"),
            radioButtons("search_connection_mode", NULL,
              choices = c("Local (on HPC)" = "local", "Remote (SSH)" = "ssh"),
              selected = if (local_sbatch) "local" else "ssh", inline = TRUE),
            conditionalPanel("input.search_connection_mode == 'ssh'",
              textInput("ssh_host", "HPC Hostname",
                placeholder = "hive.genomecenter.ucdavis.edu"),
              div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
                div(style = "flex: 1; min-width: 100px;",
                  textInput("ssh_user", "Username", value = Sys.getenv("USER"))
                ),
                div(style = "flex: 1; min-width: 80px;",
                  numericInput("ssh_port", "Port", value = 22, min = 1, max = 65535)
                )
              ),
              textInput("ssh_key_path", "SSH Key Path",
                value = paste0(Sys.getenv("HOME"), "/.ssh/id_rsa")),
              actionButton("test_ssh_btn", "Test Connection",
                icon = icon("plug"), class = "btn-outline-info btn-sm"),
              uiOutput("ssh_status_ui")
            ),

            hr(),
            tags$h6("SLURM Resources"),
            div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
              div(style = "flex: 1; min-width: 100px;",
                numericInput("diann_cpus", "CPUs:", value = 64, min = 4, max = 128, step = 4)
              ),
              div(style = "flex: 1; min-width: 100px;",
                numericInput("diann_mem_gb", "Memory (GB):", value = 512, min = 16, max = 1024, step = 16)
              )
            ),
            div(style = "display: flex; gap: 8px; flex-wrap: wrap;",
              div(style = "flex: 1; min-width: 100px;",
                numericInput("diann_time_hours", "Time (hrs):", value = 12, min = 1, max = 48)
              ),
              div(style = "flex: 1; min-width: 100px;",
                textInput("diann_partition", "Partition:", value = "high")
              )
            ),
            textInput("diann_account", "Account:", value = "genome-center-grp"),

            hr(),
            tags$h6("DIA-NN Container"),
            textInput("diann_sif_path", "Apptainer SIF Path:",
              value = "/share/proteomics/sw/diann_2.3.0.sif",
              placeholder = "/path/to/diann_2.3.0.sif"),

            hr(),
            tags$h6("Output Directory"),
            conditionalPanel("input.search_connection_mode != 'ssh'",
              shinyFiles::shinyDirButton("output_base_dir", "Select Output Folder",
                title = "Choose base output directory",
                class = "btn-outline-primary btn-sm w-100")
            ),
            conditionalPanel("input.search_connection_mode == 'ssh'",
              textInput("ssh_output_base_dir", NULL, value = "",
                placeholder = "/share/proteomics/results/")
            ),
            verbatimTextOutput("full_output_path"),

            hr(),
            uiOutput("time_estimate_ui"),

            actionButton("submit_diann", "Submit DIA-NN Search",
              class = "btn-success btn-lg w-100",
              icon = icon("rocket"),
              style = "margin-top: 10px;"),

            checkboxInput("auto_load_results", "Auto-load results when complete", TRUE),

            hr(),
            tags$h6(icon("list-check"), " Job Queue"),
            uiOutput("search_queue_ui")
          )
        )
      )
    ),

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
                  div(style = "display: flex; justify-content: flex-end; align-items: center; gap: 8px; margin-bottom: 10px;",
                    actionButton("signal_dist_info_btn", icon("question-circle"), title = "What is this?",
                      class = "btn-outline-info btn-sm"),
                    actionButton("fullscreen_signal", "\U0001F50D Fullscreen", class = "btn-outline-secondary btn-sm")
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
                  DTOutput("grid_view_table")
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

                    # Generate button
                    div(style = "text-align: center; margin-bottom: 20px;",
                      actionButton("generate_ai_summary_overview",
                        "\U0001F916 Generate AI Summary",
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
                      "(Applies to all metric tabs)")
                  ),
                  actionButton("qc_trends_info_btn", icon("question-circle"), title = "What are QC Trends?",
                    class = "btn-outline-info btn-sm")
                )
              ),

              # Metric tabs
              navset_card_tab(
                id = "qc_trends_tabs",

                nav_panel("Precursors",
                  icon = icon("dna"),
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_trend_precursors", "\U0001F50D Fullscreen",
                      class = "btn-outline-secondary btn-sm")
                  ),
                  plotlyOutput("qc_trend_plot_precursors", height = "calc(100vh - 380px)")
                ),

                nav_panel("Proteins",
                  icon = icon("shapes"),
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_trend_proteins", "\U0001F50D Fullscreen",
                      class = "btn-outline-secondary btn-sm")
                  ),
                  plotlyOutput("qc_trend_plot_proteins", height = "calc(100vh - 380px)")
                ),

                nav_panel("MS1 Signal",
                  icon = icon("signal"),
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_trend_ms1", "\U0001F50D Fullscreen",
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

                # TAB 2: DPC Fit
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

                # TAB 3: MDS Plot
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

                # TAB 4: Group Distribution
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

                # TAB 5: P-value Distribution Diagnostic
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

    nav_panel("DE Dashboard", icon = icon("table-columns"),
              # Interactive comparison selector
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

              # Responsive two-column grid (stacks on small screens)
              div(class = "de-dashboard-grid",
                # Results table card (left column)
                card(
                  card_header(
                    div(style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px;",
                      span("Results Table"),
                      div(
                        actionButton("clear_plot_selection", "Reset", class="btn-warning btn-xs"),
                        actionButton("show_violin", "\U0001F4CA Violin", class="btn-primary btn-xs"),
                        if (!is_hf_space) actionButton("show_xic", "\U0001F4C8 XICs", class="btn-info btn-xs"),
                        downloadButton("download_result_csv", "\U0001F4BE Export", class="btn-success btn-xs")
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
                      actionButton("fullscreen_volcano", "\U0001F50D Fullscreen", class="btn-outline-secondary btn-sm")
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
                open = TRUE,
                accordion_panel(
                  "Heatmap of Selected/Top Proteins",
                  icon = icon("grip"),
                  div(style = "text-align: right; margin-bottom: 10px;",
                    actionButton("fullscreen_heatmap", "\U0001F50D Fullscreen", class="btn-outline-secondary btn-sm")
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
                    div(style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px;",
                      p("Ranking by %CV (Coefficient of Variation) to find stable markers across all experimental groups.",
                        class = "text-muted small mb-0"),
                      actionButton("consistent_de_info_btn", icon("question-circle"), title = "What is this?",
                        class = "btn-outline-info btn-sm")
                    ),
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
                      div(style = "display: flex; gap: 8px;",
                        actionButton("cv_dist_info_btn", icon("question-circle"), title = "What is this?",
                          class = "btn-outline-info btn-sm"),
                        actionButton("fullscreen_cv_hist", "\U0001F50D Fullscreen",
                          class = "btn-outline-secondary btn-sm")
                      )
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
                            downloadButton("download_repro_log", "\U0001F4BE Download Reproducibility Log", class="btn-success mb-3"),
                            verbatimTextOutput("reproducible_code")
                          )
                ),
                nav_panel("Methodology",
                          card_body(
                            div(style = "display: flex; justify-content: flex-end; margin-bottom: 10px;",
                              actionButton("methodology_info_btn", icon("question-circle"), title = "About the methods",
                                class = "btn-outline-info btn-sm")
                            ),
                            verbatimTextOutput("methodology_text")
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
                  div(style = "text-align: right; margin-bottom: 10px;",
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
                  div(style = "display: flex; justify-content: flex-end; margin-bottom: 10px;",
                    actionButton("gsea_table_info_btn", icon("question-circle"), title = "Column definitions",
                      class = "btn-outline-info btn-sm")
                  ),
                  DTOutput("gsea_results_table")
                )
              )
    ),

    nav_panel("Data Chat", icon = icon("comments"),
              card(
                card_header(div(style="display: flex; justify-content: space-between; align-items: center;",
                  span("Chat with Full Data (QC + Expression)"),
                  div(style = "display: flex; gap: 8px;",
                    actionButton("data_chat_info_btn", icon("question-circle"), title = "About Data Chat",
                      class = "btn-outline-info btn-sm"),
                    downloadButton("download_chat_txt", "\U0001F4BE Save Chat", class="btn-secondary btn-sm")
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
}
