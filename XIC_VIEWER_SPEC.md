# Feature Spec: XIC Viewer for Differentially Expressed Proteins

## Overview

Add an on-demand XIC (Extracted Ion Chromatogram) viewer to DE-LIMP that lets users visualize
fragment-level chromatograms for differentially expressed proteins. When a user selects a protein
from the volcano plot, results table, or grid view, they can launch the XIC viewer to inspect
the underlying chromatographic evidence — peak shapes, co-elution of fragments, and cross-sample
consistency — directly within DE-LIMP.

This is a **validation feature**: after statistical analysis identifies DE proteins, the XIC viewer
lets bench scientists visually confirm that the quantification is backed by clean, well-shaped peaks
rather than noisy or interfered signals.

## Background: DIA-NN XIC Output

### What `--xic` produces

When DIA-NN is run with the `--xic` flag (already enabled in DE-LIMP's sbatch templates), it
generates `.xic.parquet` files — one per raw data file — in the same directory as the raw data.
These files contain extracted ion chromatograms for every identified precursor.

### XIC Parquet File Structure

Each `.xic.parquet` file is a table with the following columns:

| Column | Description |
|--------|-------------|
| `File.Name` | The raw data file this XIC came from |
| `Precursor.Id` | Precursor identifier (e.g., `PEPTIDEK2`) — the key for lookup |
| `Modified.Sequence` | Peptide sequence with modifications |
| `Stripped.Sequence` | Unmodified peptide sequence |
| `Q.Value` | Precursor-level q-value |
| `MS.Level` | 1 = MS1 (precursor), 2 = MS2 (fragment) |
| `Intensities` | Row type flag: 0 = retention times row, 1 = intensities row |
| `Retention.Times` | Row type flag: 0 = intensities row, 1 = retention times row |
| `Theoretical.Mz` | m/z of this ion |
| `Reference.Intensity` | Library/predicted intensity (for MS2 fragments) |
| `FragmentType` | `y` or `b` (for MS2 rows) |
| `FragmentCharge` | Fragment ion charge state |
| `FragmentSeriesNumber` | Position in the y/b series |
| `FragmentLossType` | `noloss`, `H2O`, `NH3`, etc. |
| `0`, `1`, `2`, ... `N` | **The actual chromatogram data points** |

### Critical: How the numbered columns work

The numbered columns (`0`, `1`, `2`, ...) contain the chromatogram data. The number of columns
is set by the `--xic N` parameter (default varies; we recommend `--xic 15` minimum for
meaningful chromatograms).

For each precursor + ion combination, there are **two rows**:
- **Retention times row** (`Retention.Times = 1`, `Intensities = 0`): The numbered columns
  contain RT values (in minutes) for each data point
- **Intensities row** (`Retention.Times = 0`, `Intensities = 1`): The numbered columns contain
  intensity values at each corresponding RT point

So to plot a chromatogram for one fragment, you read the RT row and the intensity row for the
same `Precursor.Id` + `Theoretical.Mz` combination, then plot RT (x-axis) vs intensity (y-axis).

### File Size Considerations

XIC parquet files can be very large:
- **Typical size**: 100 MB – 2 GB per raw file, depending on library size and `--xic N` setting
- **Total for an experiment**: 10–50+ GB for a 20-sample experiment
- **Why so large**: Every precursor × every fragment × 2 rows (RT + intensity) × N data points

This is the primary engineering challenge. DE-LIMP **must not** attempt to load entire XIC files
into memory. The solution is selective, on-demand reading using Apache Arrow's predicate pushdown.

## Deployment Scope

### Local / HPC: Full support
XIC files reside on the filesystem alongside the raw data. Users point DE-LIMP to the XIC
directory and everything works.

### Hugging Face: Not supported
XIC files are far too large for Hugging Face deployment. The XIC viewer UI elements should be
**conditionally hidden** when XIC files are not available, with an informational message if
the user tries to access the feature.

## Implementation Spec

### Architecture: On-Demand Lazy Loading

```
User selects protein → Look up precursors in report.parquet
                     → For each sample's .xic.parquet:
                        Arrow::open_dataset() with predicate pushdown
                        Filter: Precursor.Id %in% target_precursors
                        Collect only matching rows
                     → Reshape for plotting
                     → Render interactive plotly chromatograms
```

The key insight is that `arrow::open_dataset()` with `filter()` + `collect()` uses predicate
pushdown to read only the relevant row groups from the parquet files, avoiding loading the
entire file into memory. This is effectively instantaneous even for very large XIC files.

### New Reactive Values

Add to the existing `reactiveValues` block (around line 1040 in v2.1):

```r
values <- reactiveValues(
  # ... existing values ...
  raw_data = NULL, metadata = NULL, fit = NULL, y_protein = NULL,
  dpc_fit = NULL, status = "Waiting...", design = NULL, qc_stats = NULL,
  plot_selected_proteins = NULL, chat_history = list(),
  current_file_uri = NULL, gsea_results = NULL,
  repro_log = c(...),
  color_plot_by_de = FALSE,
  grid_selected_protein = NULL,
  temp_violin_target = NULL,
  diann_norm_detected = "unknown",

  # === NEW: XIC Viewer ===
  xic_dir = NULL,              # Path to directory containing .xic.parquet files
  xic_available = FALSE,       # Whether XIC files were detected
  xic_protein = NULL,          # Currently selected protein for XIC viewing
  xic_data = NULL,             # Loaded XIC data (tibble, only for selected protein)
  xic_report_map = NULL,       # Protein → Precursor mapping from report
  uploaded_report_path = NULL   # Path to uploaded report.parquet for re-reading
)
```

### Store Report Path on Upload

The existing upload handler (line ~1140 in v2.1) reads the report via `limpa::readDIANN()`
but doesn't retain the file path. Add path storage:

```r
# In observeEvent(input$report_file, { ... })  (around line 1140)
# Add after: values$raw_data <- limpa::readDIANN(...)
values$uploaded_report_path <- input$report_file$datapath

# Similarly in the example data loader (around line 1090):
values$uploaded_report_path <- temp_file
```

### UI: XIC Directory Input

Add to the **sidebar** (around line 436), after the existing "4. AI Chat" section:

```r
# Add after the AI Chat section, before sidebar closes
hr(),
h5("5. XIC Viewer"),
p(class = "text-muted small",
  "Load .xic.parquet files from DIA-NN to inspect chromatograms."),
textInput("xic_dir_input", "XIC Directory Path:",
  placeholder = "/path/to/diann/output/"),
actionButton("xic_load_dir", "Load XICs", class = "btn-outline-info btn-sm w-100",
  icon = icon("wave-square")),
uiOutput("xic_status_badge")
```

**Why `textInput` instead of `shinyDirButton`**: Avoids adding `shinyFiles` as a new
dependency. The directory chooser widget doesn't work well on headless/HPC environments anyway.
Users on HPC/local deployments are comfortable pasting paths. On Hugging Face, this entire
section is hidden since XIC files won't be available.

### Server: XIC Directory Detection

```r
# Handle XIC directory loading
observeEvent(input$xic_load_dir, {
  req(input$xic_dir_input)
  xic_path <- trimws(input$xic_dir_input)

  if (!dir.exists(xic_path)) {
    showNotification("Directory not found. Check the path and try again.",
      type = "error", duration = 5)
    values$xic_available <- FALSE
    return()
  }

  # Look for .xic.parquet files
  xic_files <- list.files(xic_path, pattern = "\\.xic\\.parquet$",
                          full.names = TRUE, recursive = TRUE)

  if (length(xic_files) > 0) {
    values$xic_dir <- xic_path
    values$xic_available <- TRUE

    # Pre-load precursor mapping if report is available
    if (!is.null(values$uploaded_report_path)) {
      tryCatch({
        values$xic_report_map <- arrow::read_parquet(
          values$uploaded_report_path,
          col_select = c("Protein.Group", "Precursor.Id",
                         "Modified.Sequence", "Stripped.Sequence", "File.Name")
        ) %>% distinct()

        message(paste("XIC precursor map loaded:",
                      n_distinct(values$xic_report_map$Protein.Group), "proteins,",
                      n_distinct(values$xic_report_map$Precursor.Id), "precursors"))
      }, error = function(e) {
        warning(paste("Could not load precursor mapping:", e$message))
        values$xic_report_map <- NULL
      })
    }

    showNotification(
      paste("Found", length(xic_files), "XIC files. Select a protein to view chromatograms."),
      type = "message", duration = 5)
  } else {
    values$xic_available <- FALSE
    showNotification(
      "No .xic.parquet files found in selected directory.",
      type = "warning", duration = 5)
  }
})

output$xic_status_badge <- renderUI({
  if (values$xic_available) {
    xic_files <- list.files(values$xic_dir, pattern = "\\.xic\\.parquet$")
    div(class = "alert alert-success py-1 px-2 mt-2 mb-0",
      style = "font-size: 0.85em;",
      icon("check-circle"),
      paste(length(xic_files), "XIC files ready"))
  } else {
    NULL
  }
})
```

### UI: XIC Viewer Trigger Buttons

Add "View XICs" buttons to two interaction points, following the existing v2.1 patterns:

#### 1. DE Dashboard Results Table Header (line ~850)

Add alongside the existing `Violin`, `Export`, and `Reset` buttons:

```r
# In the DE Dashboard card_header, alongside existing buttons:
# Current v2.1 (line 850-854):
div(
  actionButton("clear_plot_selection", "Reset", class="btn-warning btn-xs"),
  actionButton("show_violin", "📊 Violin", class="btn-primary btn-xs"),
  actionButton("show_xic", "📈 XICs", class="btn-info btn-xs"),  # <-- NEW
  downloadButton("download_result_csv", "💾 Export", class="btn-success btn-xs")
)
```

#### 2. Grid View Protein Detail Modal (line ~1965)

Add alongside the existing "Back to Grid" button in the grid view modal footer:

```r
# Current v2.1 (line 1965):
showModal(modalDialog(
  title = paste("Expression Plot:", selected_id),
  size = "xl",
  plotOutput("violin_plot_grid", height = "600px"),
  footer = tagList(
    actionButton("show_xic_from_grid", "📈 XICs", class="btn-info"),  # <-- NEW
    actionButton("back_to_grid", "Back to Grid", class="btn-info"),
    modalButton("Close")
  ),
  easyClose = TRUE
))
```

### Server: XIC Viewer Modal

The XIC viewer opens as a full-width modal, following the same pattern as the fullscreen
views and violin plot popup already in v2.1:

```r
# === XIC from DE Dashboard ===
observeEvent(input$show_xic, {
  # Same guard pattern as show_violin (line 3601-3604 in v2.1)
  if (!values$xic_available) {
    showNotification("⚠️ Load XIC files first (sidebar → '5. XIC Viewer')", type = "warning")
    return()
  }
  if (is.null(values$plot_selected_proteins) || length(values$plot_selected_proteins) == 0) {
    showNotification("⚠️ Select a protein in the Volcano Plot or Table first!", type = "warning")
    return()
  }

  # Use first selected protein (one XIC at a time)
  values$xic_protein <- values$plot_selected_proteins[1]
  show_xic_modal(session, values)
})

# === XIC from Grid View ===
observeEvent(input$show_xic_from_grid, {
  if (!values$xic_available) {
    showNotification("⚠️ Load XIC files first (sidebar → '5. XIC Viewer')", type = "warning")
    return()
  }
  req(values$grid_selected_protein)
  values$xic_protein <- values$grid_selected_protein

  # Close the current grid modal, open XIC modal
  removeModal()
  show_xic_modal(session, values)
})

# === Shared modal function ===
show_xic_modal <- function(session, values) {
  showModal(modalDialog(
    title = div(style = "display: flex; align-items: center; gap: 10px;",
      icon("wave-square"),
      span("XIC Chromatograms:"),
      span(values$xic_protein, style = "color: #667eea; font-weight: bold;"),
      span(textOutput("xic_precursor_count", inline = TRUE),
        style = "font-size: 0.8em; color: #6c757d; margin-left: 10px;")
    ),

    # Controls row (compact inline, matches v2.1 .controls-inline pattern)
    div(class = "controls-inline mb-3",
      selectInput("xic_display_mode", "Display:",
        choices = c("Overlay (all fragments)" = "overlay",
                    "Facet by fragment" = "facet",
                    "Facet by sample" = "sample_facet"),
        selected = "overlay", width = "220px"),

      selectInput("xic_precursor_select", "Precursor:",
        choices = NULL, width = "280px"),  # populated reactively

      selectInput("xic_group_filter", "Filter Group:",
        choices = NULL, width = "180px"),  # populated reactively

      checkboxInput("xic_show_ms1", "Show MS1", value = FALSE)
    ),

    # Main plot area with viewport-relative height (matches v2.1 pattern)
    plotlyOutput("xic_plot", height = "calc(100vh - 380px)"),

    # Info panel
    div(class = "mt-2 p-2 bg-light rounded",
      style = "font-size: 0.85em;",
      uiOutput("xic_info_panel")
    ),

    size = "xl",
    easyClose = TRUE,
    footer = div(style = "display: flex; gap: 8px; align-items: center;",
      actionButton("xic_prev_protein", icon("arrow-left"),
        label = "Prev", class = "btn-outline-secondary btn-sm"),
      actionButton("xic_next_protein", icon("arrow-right"),
        label = "Next", class = "btn-outline-secondary btn-sm"),
      downloadButton("xic_download_plot", "💾 Download",
        class = "btn-outline-success btn-sm"),
      modalButton("Close")
    )
  ))
}
```

### Server: Core XIC Data Loading

Uses Arrow predicate pushdown to read only target rows from potentially massive files:

```r
# Core XIC loading function
load_xic_for_protein <- function(xic_dir, protein_id, report_map) {
  # Step 1: Look up precursors for this protein
  target_precursors <- report_map %>%
    filter(Protein.Group == protein_id) %>%
    distinct(Precursor.Id) %>%
    pull(Precursor.Id)

  if (length(target_precursors) == 0) return(NULL)

  # Step 2: Use Arrow dataset for predicate pushdown across all XIC files
  tryCatch({
    ds <- arrow::open_dataset(xic_dir, format = "parquet")

    # Arrow pushes this filter down — only matching row groups are read from disk
    xic_data <- ds %>%
      filter(Precursor.Id %in% target_precursors) %>%
      collect()

    if (nrow(xic_data) == 0) return(NULL)
    xic_data
  }, error = function(e) {
    # Fallback: read files individually (handles non-uniform schemas)
    xic_files <- list.files(xic_dir, pattern = "\\.xic\\.parquet$",
                            full.names = TRUE, recursive = TRUE)
    xic_list <- lapply(xic_files, function(f) {
      tryCatch({
        arrow::read_parquet(f) %>%
          filter(Precursor.Id %in% target_precursors)
      }, error = function(e2) NULL)
    })
    bind_rows(Filter(Negate(is.null), xic_list))
  })
}
```

### Server: XIC Data Reshaping

Pairs the RT and intensity rows and creates a tidy long-format dataframe for ggplot:

```r
reshape_xic_for_plotting <- function(xic_raw, metadata) {
  # Identify numbered columns (the chromatogram data points)
  num_cols <- names(xic_raw)[grepl("^\\d+$", names(xic_raw))]

  if (length(num_cols) == 0) {
    warning("No numbered columns found in XIC data — check --xic parameter")
    return(NULL)
  }

  # Create a join key for pairing RT and intensity rows
  make_key <- function(df) {
    paste(df$File.Name, df$Precursor.Id, df$MS.Level,
          df$Theoretical.Mz, df$FragmentType, df$FragmentCharge,
          df$FragmentSeriesNumber, df$FragmentLossType, sep = "|")
  }

  # Separate RT rows and intensity rows
  rt_rows <- xic_raw %>% filter(Retention.Times == 1)
  int_rows <- xic_raw %>% filter(Intensities == 1)

  rt_keys <- make_key(rt_rows)
  int_keys <- make_key(int_rows)

  # Pivot RT rows to long format
  rt_long <- rt_rows %>%
    mutate(.key = rt_keys) %>%
    select(.key, all_of(num_cols)) %>%
    pivot_longer(cols = all_of(num_cols), names_to = "point_idx", values_to = "RT")

  # Pivot intensity rows to long format
  int_long <- int_rows %>%
    mutate(.key = int_keys) %>%
    select(.key, File.Name, Precursor.Id, Modified.Sequence, MS.Level,
           Theoretical.Mz, Reference.Intensity, FragmentType, FragmentCharge,
           FragmentSeriesNumber, FragmentLossType, all_of(num_cols)) %>%
    pivot_longer(cols = all_of(num_cols), names_to = "point_idx", values_to = "Intensity")

  # Join RT and intensity by key + point index
  xic_plot <- inner_join(
    rt_long %>% select(.key, point_idx, RT),
    int_long,
    by = c(".key", "point_idx")
  ) %>%
    mutate(
      RT = as.numeric(RT),
      Intensity = as.numeric(Intensity),
      # Human-readable fragment label
      Fragment.Label = case_when(
        MS.Level == 1 ~ paste0("MS1 (", round(as.numeric(Theoretical.Mz), 2), ")"),
        TRUE ~ paste0(FragmentType, FragmentSeriesNumber,
                      ifelse(as.integer(FragmentCharge) > 1,
                             paste0("+", FragmentCharge), ""),
                      ifelse(FragmentLossType != "noloss",
                             paste0("-", FragmentLossType), ""))
      )
    ) %>%
    # Remove zero-padded points
    filter(!(Intensity == 0 & RT == 0)) %>%
    # Add sample metadata (with fuzzy File.Name matching)
    mutate(File.Name.Base = basename(tools::file_path_sans_ext(
      tools::file_path_sans_ext(File.Name)))) %>%
    left_join(
      metadata %>%
        mutate(File.Name.Base = basename(tools::file_path_sans_ext(File.Name))) %>%
        select(File.Name.Base, Group, ID),
      by = "File.Name.Base"
    ) %>%
    select(.key, File.Name, ID, Group, Precursor.Id, Modified.Sequence,
           MS.Level, Fragment.Label, Theoretical.Mz, Reference.Intensity,
           RT, Intensity)

  return(xic_plot)
}
```

### Server: Reactive XIC Data Loading on Protein Change

```r
# Load XIC data when protein changes
observe({
  req(values$xic_protein, values$xic_available, values$xic_dir,
      values$xic_report_map)

  withProgress(message = "Loading chromatograms...", {
    incProgress(0.3, detail = "Reading XIC files...")

    tryCatch({
      xic_raw <- load_xic_for_protein(
        xic_dir = values$xic_dir,
        protein_id = values$xic_protein,
        report_map = values$xic_report_map
      )

      if (!is.null(xic_raw) && nrow(xic_raw) > 0) {
        incProgress(0.6, detail = "Reshaping data...")

        values$xic_data <- reshape_xic_for_plotting(xic_raw, values$metadata)

        # Update precursor selector
        precursors <- unique(values$xic_data$Precursor.Id)
        updateSelectInput(session, "xic_precursor_select",
          choices = c("All Precursors" = "all",
                      setNames(precursors, precursors)),
          selected = "all")

        # Update group filter
        groups <- unique(values$metadata$Group)
        updateSelectInput(session, "xic_group_filter",
          choices = c("All Groups" = "all", setNames(groups, groups)),
          selected = "all")

      } else {
        values$xic_data <- NULL
        showNotification(
          paste("No XIC data found for", values$xic_protein,
                "— this protein may have been identified via MBR."),
          type = "warning", duration = 5)
      }
    }, error = function(e) {
      values$xic_data <- NULL
      showNotification(paste("Error loading XICs:", e$message),
        type = "error", duration = 8)
    })
  })
})

# Precursor count display
output$xic_precursor_count <- renderText({
  req(values$xic_data)
  n_prec <- n_distinct(values$xic_data$Precursor.Id)
  n_frag <- values$xic_data %>%
    filter(MS.Level == 2) %>%
    distinct(Fragment.Label) %>%
    nrow()
  paste0(n_prec, " precursor(s), ", n_frag, " fragments")
})
```

### Server: XIC Plot Rendering

```r
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

  # MS level filter
  if (!isTRUE(input$xic_show_ms1)) {
    xic_plot <- xic %>% filter(MS.Level == 2)
    if (nrow(xic_plot) == 0) xic_plot <- xic  # fallback to MS1
  } else {
    xic_plot <- xic
  }

  # Cap very large proteins (e.g. Titin) at top 6 precursors by fragment count
  if (n_distinct(xic_plot$Precursor.Id) > 6 &&
      input$xic_precursor_select == "all") {
    top_prec <- xic_plot %>%
      filter(MS.Level == 2) %>%
      count(Precursor.Id) %>%
      slice_max(n, n = 6) %>%
      pull(Precursor.Id)
    xic_plot <- xic_plot %>% filter(Precursor.Id %in% top_prec)
    showNotification(
      paste("Showing top 6 of", n_distinct(xic$Precursor.Id),
            "precursors. Use the selector to view others."),
      type = "message", duration = 4)
  }

  # Tooltip text
  xic_plot <- xic_plot %>%
    mutate(tooltip_text = paste0(
      "<b>Sample:</b> ", ID, " (", Group, ")",
      "<br><b>Fragment:</b> ", Fragment.Label,
      "<br><b>RT:</b> ", round(RT, 3), " min",
      "<br><b>Intensity:</b> ", format(round(Intensity), big.mark = ",")
    ))

  if (display_mode == "overlay") {
    # === OVERLAY: All fragments per sample, faceted by sample ===
    p <- ggplot(xic_plot,
        aes(x = RT, y = Intensity, color = Fragment.Label,
            group = interaction(File.Name, Fragment.Label),
            text = tooltip_text)) +
      geom_line(alpha = 0.7, linewidth = 0.5) +
      facet_wrap(~ paste0(ID, " (", Group, ")"), scales = "free_y") +
      theme_minimal() +
      labs(
        title = paste("Fragment XICs —", values$xic_protein),
        subtitle = "Each panel = one sample, colors = fragment ions",
        x = "Retention Time (min)", y = "Intensity", color = "Fragment"
      ) +
      theme(legend.position = "bottom", legend.text = element_text(size = 7),
            strip.text = element_text(size = 8))

  } else if (display_mode == "facet") {
    # === FACET BY FRAGMENT: Compare groups per fragment ion ===
    p <- ggplot(xic_plot,
        aes(x = RT, y = Intensity, color = Group,
            group = File.Name, text = tooltip_text)) +
      geom_line(alpha = 0.6, linewidth = 0.4) +
      facet_wrap(~Fragment.Label, scales = "free_y") +
      theme_minimal() +
      labs(
        title = paste("Fragment XICs —", values$xic_protein),
        subtitle = "Each panel = one fragment ion, colors = groups",
        x = "Retention Time (min)", y = "Intensity", color = "Group"
      ) +
      theme(strip.text = element_text(size = 8))

  } else if (display_mode == "sample_facet") {
    # === FACET BY SAMPLE: All fragments overlaid per sample ===
    p <- ggplot(xic_plot,
        aes(x = RT, y = Intensity, color = Fragment.Label,
            group = Fragment.Label, text = tooltip_text)) +
      geom_line(alpha = 0.7, linewidth = 0.5) +
      facet_wrap(~ paste0(ID, " (", Group, ")"), scales = "free_y") +
      theme_minimal() +
      labs(
        title = paste("Fragment XICs —", values$xic_protein),
        subtitle = "Each panel = one sample, all fragments overlaid",
        x = "Retention Time (min)", y = "Intensity", color = "Fragment"
      ) +
      theme(legend.position = "bottom", legend.text = element_text(size = 7),
            strip.text = element_text(size = 8))
  }

  ggplotly(p, tooltip = "text") %>%
    layout(legend = list(orientation = "h", y = -0.15),
           margin = list(b = 80)) %>%
    config(displayModeBar = TRUE)
})
```

### Server: XIC Info Panel

```r
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
        paste0(round(rt_range[1], 2), " – ", round(rt_range[2], 2), " min")
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
```

### Server: Protein Navigation (Prev/Next)

Step through DE protein list from within the XIC viewer, matching the v2.1 contrast
selector approach:

```r
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
```

## Integration Points in v2.1 Codebase

| What | Where in v2.1 | Action |
|------|---------------|--------|
| Reactive values | Line ~1040 | Add 5 new values |
| Store report path | Lines ~1147, ~1109 | Add `values$uploaded_report_path` |
| Sidebar XIC section | After line ~436 | Add section 5 UI |
| DE Dashboard buttons | Line ~852 | Add `show_xic` button |
| Grid view modal | Line ~1965 | Add `show_xic_from_grid` button |
| CSS `.controls-inline` | Already exists (line ~388) | Reuse for XIC controls |
| Height pattern | `calc(100vh - 380px)` | Reuse for XIC plot |
| Modal pattern | Lines ~3616-3622 | Follow violin popup pattern |
| `useShinyjs()` | Line 338 | Already present |
| Gradient styling | Lines 832, 537 | Reuse `#667eea → #764ba2` palette |
| `navset_card_tab` | Used throughout v2.1 | Follow same sub-tab pattern if XIC grows |

## Dependencies

### New R Packages: **None**

| Package | Source | Purpose | Already in DE-LIMP? |
|---------|--------|---------|---------------------|
| `arrow` | CRAN | Read `.xic.parquet` with predicate pushdown | **Yes** (via `limpa`) |
| `shinyjs` | CRAN | Show/hide elements | **Yes** (line 338) |
| `plotly` | CRAN | Interactive chromatogram plots | **Yes** |
| `ggplot2` | CRAN | Plot construction | **Yes** |
| `dplyr` / `tidyr` | CRAN | Data reshaping | **Yes** |

No new dependencies needed. Everything required is already in the v2.1 stack.

## Edge Cases & Error Handling

1. **File.Name mismatch**: XIC files may contain full paths while metadata has basenames.
   Use `basename()` fuzzy matching (implemented in `reshape_xic_for_plotting`).

2. **Large protein (Titin, etc.)**: Cap at top 6 precursors by fragment count. Show
   notification and allow manual precursor selection.

3. **XIC from different DIA-NN run**: If precursor IDs don't match, show clear error.

4. **Memory pressure**: If filtered XIC data exceeds ~50 MB, warn and offer subset.

5. **MBR-only proteins**: Some DE proteins identified via Match Between Runs may lack
   XIC data. Show explanatory message.

6. **XIC window too short**: If `--xic` was set to a very small N (e.g., 2-3), the
   chromatograms will have too few points to be useful. Detect and warn.

7. **Non-uniform XIC schemas**: Different DIA-NN versions or runs may have different
   numbers of data point columns. The fallback in `load_xic_for_protein` handles this
   by reading files individually when `open_dataset()` fails.

## Future Enhancements

- **Peak boundary markers**: Vertical lines showing DIA-NN's integration boundaries
  (requires `RT.Start` / `RT.Stop` from main report)
- **Reference spectrum overlay**: Show library/predicted spectrum alongside observed
- **Batch XIC export**: Multi-page PDF of XICs for all significant proteins
  (supplementary materials for publications)
- **Co-elution quality score**: Pearson correlation between fragment traces as a
  quantitative confidence metric
- **Auto-detection from HPC integration**: When DIA-NN HPC submission is active,
  auto-populate XIC directory from the job output path
- **Gemini integration**: Let AI chat reference XIC quality when discussing proteins

## Testing Checklist

- [ ] XIC directory path input accepts valid paths
- [ ] .xic.parquet files are correctly detected and counted
- [ ] Status badge shows correct file count
- [ ] Precursor lookup maps protein → precursors correctly
- [ ] Arrow predicate pushdown reads only target rows (profile memory usage)
- [ ] RT/intensity row pairing produces correct chromatograms
- [ ] All three display modes render correctly (overlay, facet, sample_facet)
- [ ] Precursor selector filters correctly
- [ ] Group filter works
- [ ] MS1 checkbox toggles MS1 traces
- [ ] Protein navigation (prev/next) steps through DE proteins
- [ ] Loading progress bar shows during data fetch
- [ ] Graceful handling when no XIC data found for a protein
- [ ] Info panel shows correct DE stats for current contrast
- [ ] Download plot works
- [ ] Memory stays reasonable with large XIC files
- [ ] Feature buttons disabled/warning when XICs not loaded
- [ ] File.Name fuzzy matching handles path differences
- [ ] Large protein cap at 6 precursors works
- [ ] Works with both library-mode and library-free DIA-NN output
