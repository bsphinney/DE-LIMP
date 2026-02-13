# Feature Spec: Normalization Before/After Diagnostic Plot

## Overview

Add a side-by-side visualization to the **QC Plots** tab that shows how normalization transforms
the data distributions across samples. This helps users confirm normalization is working correctly
and understand the magnitude of correction being applied.

## Background: The Normalization Chain

Data arriving in DE-LIMP passes through **up to three normalization layers**:

### Layer 1: DIA-NN Normalization (Pre-DE-LIMP)
- **Default**: RT-dependent normalization (recommended, always-on unless user disabled it)
- **How it works**: Selects top N precursors with lowest CVs, normalizes their summed intensities
  to be equal across runs, then applies RT-bin-specific correction factors
- **Applied at**: Precursor level, before DE-LIMP ever sees the data
- **Key columns in parquet file**:
  - `Precursor.Quantity` = raw, unnormalized intensity
  - `Precursor.Normalised` = after DIA-NN normalization
  - Normalization factor = `Precursor.Normalised / Precursor.Quantity`
- **Critical**: `limpa::readDIANN()` reads `Precursor.Normalised` by default (its `qty.column`
  parameter defaults to `"Precursor.Normalised"`)

### Layer 2: limpa DPC-CN (`dpcCN()`)
- Cyclic Normalization using the Detection Probability Curve
- Corrects for systematic technical variation across runs
- Accounts for intensity-dependent missing values
- **Input**: `values$raw_data` (EList from `readDIANN()` — already DIA-NN normalized)
- **Output**: DPC fit object (`values$dpc_fit`)

### Layer 3: limpa DPC-Quant (`dpcQuant()`)
- Aggregates precursors → proteins using modified maxLFQ with DPC-based missing value handling
- Uses the DPC fit from Layer 2
- **Output**: `values$y_protein` (protein-level EList with complete data, no NAs)

### What this means for the diagnostic

`values$raw_data$E` is **precursor-level, log2, already DIA-NN normalized** data (with NAs).
`values$y_protein$E` is **protein-level, log2, DPC-CN normalized + DPC-Quant aggregated** data
(complete, no NAs).

Comparing these two shows the combined effect of limpa's normalization + protein aggregation.
They are at different granularities (precursor vs protein), but per-sample distribution comparison
(box plots, density plots) is still valid and informative — users can see whether medians align and
distributions become more uniform.

## Implementation Spec

### Location in UI

Add to the **QC Plots** tab (`nav_panel("QC Plots")`), as a new `layout_columns` row **below**
the existing DPC Fit + MDS row and above the Group QC Distribution violin.

### UI Components

```r
layout_columns(col_widths = c(12),
  card(
    card_header(
      div(style = "display: flex; justify-content: space-between; align-items: center;",
        span("Normalization Diagnostic"),
        div(
          # Info badge about DIA-NN normalization status
          uiOutput("diann_norm_status_badge", inline = TRUE),
          actionButton("fullscreen_norm_diag", "🔍 View Fullscreen", class = "btn-info btn-sm")
        )
      )
    ),
    card_body(
      radioButtons("norm_diag_type", "View:",
        choices = c("Box Plots" = "boxplot", "Density Overlay" = "density"),
        inline = TRUE
      ),
      plotlyOutput("norm_diagnostic_plot", height = "450px")
    )
  )
)
```

### DIA-NN Normalization Detection

At data load time (in both the `input$report_file` and `input$load_example` observers), after
calling `readDIANN()`, read the parquet file independently to check whether DIA-NN normalization
was applied:

```r
# After readDIANN() call, detect DIA-NN normalization status
values$diann_norm_detected <- tryCatch({
  raw_parquet <- arrow::read_parquet(file_path,
    col_select = c("Precursor.Quantity", "Precursor.Normalised"))

  has_both_cols <- all(c("Precursor.Quantity", "Precursor.Normalised") %in% names(raw_parquet))

  if (has_both_cols) {
    # Sample a subset to check if values differ
    sample_rows <- head(raw_parquet, 1000)
    ratio <- sample_rows$Precursor.Normalised / sample_rows$Precursor.Quantity
    # If all ratios ≈ 1, normalization was OFF (or --no-norm was used)
    ratios_vary <- sd(ratio, na.rm = TRUE) > 0.001
    if (ratios_vary) "on" else "off"
  } else {
    "unknown"  # Columns missing, can't determine
  }
}, error = function(e) "unknown")
```

**Add to reactive values initialization:**
```r
values <- reactiveValues(
  # ... existing values ...
  diann_norm_detected = "unknown"  # "on", "off", or "unknown"
)
```

### DIA-NN Status Badge (UI Output)

```r
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
```

### Main Plot Logic

```r
output$norm_diagnostic_plot <- renderPlotly({
  req(values$raw_data, values$y_protein, values$metadata)

  # --- Pre-normalization: per-sample medians from precursor data ---
  pre_mat <- values$raw_data$E  # precursor-level, log2, NAs present
  pre_medians <- data.frame(
    Sample = colnames(pre_mat),
    Median = apply(pre_mat, 2, median, na.rm = TRUE),
    Stage = "Pre-Normalization\n(Precursor-level, DIA-NN normalized)"
  )

  # --- Post-normalization: per-sample medians from protein data ---
  post_mat <- values$y_protein$E  # protein-level, log2, no NAs
  post_medians <- data.frame(
    Sample = colnames(post_mat),
    Median = apply(post_mat, 2, median, na.rm = TRUE),
    Stage = "Post-Normalization\n(Protein-level, DPC-CN + DPC-Quant)"
  )

  # Add group info
  meta <- values$metadata
  pre_medians$Group <- meta$Group[match(pre_medians$Sample, meta$File.Name)]
  post_medians$Group <- meta$Group[match(post_medians$Sample, meta$File.Name)]

  if (input$norm_diag_type == "boxplot") {
    # === BOX PLOT VIEW ===
    # Build long-form data for side-by-side box plots
    pre_long <- as.data.frame(pre_mat) %>%
      pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
      mutate(Stage = "Before\n(Precursor)") %>%
      filter(!is.na(Log2Intensity))

    post_long <- as.data.frame(post_mat) %>%
      pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
      mutate(Stage = "After\n(Protein)")

    # Add group info to both
    pre_long$Group <- meta$Group[match(pre_long$Sample, meta$File.Name)]
    post_long$Group <- meta$Group[match(post_long$Sample, meta$File.Name)]

    combined <- bind_rows(pre_long, post_long)
    combined$Stage <- factor(combined$Stage, levels = c("Before\n(Precursor)", "After\n(Protein)"))

    # Order samples by group then name
    sample_order <- meta %>% arrange(Group, File.Name) %>% pull(File.Name)
    combined$Sample <- factor(combined$Sample, levels = sample_order)

    # Short sample labels (use ID numbers)
    combined$SampleID <- meta$ID[match(combined$Sample, meta$File.Name)]

    p <- ggplot(combined, aes(x = factor(SampleID), y = Log2Intensity, fill = Group)) +
      geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
      facet_wrap(~Stage, scales = "free_y", ncol = 2) +
      theme_minimal() +
      labs(
        title = "Normalization Effect: Per-Sample Intensity Distributions",
        subtitle = "Medians should align after normalization",
        x = "Sample ID",
        y = "Log2 Intensity"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

    ggplotly(p, tooltip = c("x", "y")) %>% layout(boxmode = "group")

  } else {
    # === DENSITY OVERLAY VIEW ===
    # One density curve per sample, colored by group, before vs after

    pre_long <- as.data.frame(pre_mat) %>%
      pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
      mutate(Stage = "Before (Precursor)") %>%
      filter(!is.na(Log2Intensity))

    post_long <- as.data.frame(post_mat) %>%
      pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
      mutate(Stage = "After (Protein)")

    pre_long$Group <- meta$Group[match(pre_long$Sample, meta$File.Name)]
    post_long$Group <- meta$Group[match(post_long$Sample, meta$File.Name)]

    combined <- bind_rows(pre_long, post_long)
    combined$Stage <- factor(combined$Stage,
      levels = c("Before (Precursor)", "After (Protein)"))

    p <- ggplot(combined, aes(x = Log2Intensity, color = Group, group = Sample)) +
      geom_density(alpha = 0.3, linewidth = 0.4) +
      facet_wrap(~Stage, ncol = 2) +
      theme_minimal() +
      labs(
        title = "Normalization Effect: Per-Sample Density Curves",
        subtitle = "Curves should overlap more tightly after normalization",
        x = "Log2 Intensity",
        y = "Density"
      )

    ggplotly(p)
  }
})
```

### Fullscreen Modal (same pattern as QC Trends)

```r
observeEvent(input$fullscreen_norm_diag, {
  showModal(modalDialog(
    title = "Normalization Diagnostic - Fullscreen View",
    plotlyOutput("norm_diagnostic_plot_fullscreen", height = "700px"),
    size = "xl",
    easyClose = TRUE,
    footer = modalButton("Close")
  ))
})

output$norm_diagnostic_plot_fullscreen <- renderPlotly({
  # Same logic as norm_diagnostic_plot — extract to a shared reactive
  # (see Refactoring Notes below)
})
```

### Refactoring Notes

To avoid code duplication between regular and fullscreen views, extract the plot logic into a
reactive expression (same pattern used for `generate_qc_trend_plot`):

```r
generate_norm_diagnostic_plot <- reactive({
  req(values$raw_data, values$y_protein, values$metadata)
  # ... all the plot logic from above ...
})

output$norm_diagnostic_plot <- renderPlotly({ generate_norm_diagnostic_plot() })
output$norm_diagnostic_plot_fullscreen <- renderPlotly({ generate_norm_diagnostic_plot() })
```

## What Users Should See

### Healthy Data (Normalization Working Well)
- **Before**: Box plot medians scattered at different heights across samples
- **After**: Box plot medians well-aligned at similar heights
- **Density**: Curves converge and overlap more tightly after normalization
- DIA-NN badge shows "ON" (green)

### Already-Normalized Data (e.g., `noNorm.parquet` example file)
- **Before**: Medians may already be somewhat scattered (raw DIA-NN output)
- **After**: limpa's DPC-CN provides additional alignment
- DIA-NN badge shows "OFF" (yellow)
- This is actually the ideal scenario for showing DPC-CN's effect in isolation

### Problematic Data
- **After** distributions still highly scattered → possible batch effects not corrected by normalization
- One sample's distribution dramatically different → possible outlier/failed sample
- Very wide distributions post-normalization → high missing value rate, DPC imputation active

## Important Caveats to Display

Add an info box below the plot (or as a collapsible details section):

```r
div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;
             font-size: 0.85em; color: #6c757d;",
  icon("info-circle"),
  strong(" Note: "),
  "'Before' shows precursor-level data (with NAs) as read by limpa from DIA-NN's ",
  tags$code("Precursor.Normalised"), " column. ",
  "'After' shows protein-level data after DPC-CN normalization and DPC-Quant aggregation. ",
  "The two panels differ in both normalization state AND granularity (precursor vs protein), ",
  "so distributions will naturally differ in shape. Focus on whether ",
  strong("per-sample medians align"), " after normalization."
)
```

## Files to Modify

1. **DE-LIMP.R** (and then copy to app.R):
   - Add `diann_norm_detected` to `reactiveValues()` initialization (~line 555)
   - Add DIA-NN norm detection code to both data loading observers (~lines 599 and 640)
   - Add UI elements to `nav_panel("QC Plots")` (~line 425)
   - Add server-side reactive + render functions (after existing QC renderers ~line 1000)

2. **Dockerfile**: No changes needed (no new packages required — uses existing ggplot2, plotly, tidyr)

3. **CLAUDE.md**: Add to Recent Changes section documenting the new feature

4. **USER_GUIDE.md**: Add section about interpreting the normalization diagnostic

## Testing Checklist

- [ ] Plot renders correctly with example data (which has DIA-NN norm OFF)
- [ ] Plot renders correctly with user-uploaded data (DIA-NN norm likely ON)
- [ ] DIA-NN badge correctly detects ON vs OFF vs unknown
- [ ] Box plot view works and medians are visible
- [ ] Density overlay view works and curves are distinguishable
- [ ] Fullscreen modal works
- [ ] Plot handles edge cases: single group, many samples (>20), very few samples (2-3)
- [ ] No errors when pipeline hasn't been run yet (plot should just not render — `req()` guards)
- [ ] Tooltip/hover info works in plotly
- [ ] Performance OK with large datasets (precursor matrix can be huge — may need to subsample)

## Performance Consideration

The precursor matrix (`values$raw_data$E`) can have 50,000+ rows. For the density plot, this is
fine (ggplot2 handles it). For box plots, the `pivot_longer` creates a very large data frame.

**Mitigation**: For box plots, compute summary statistics (quartiles, whiskers) per sample instead
of passing all data points. Or subsample to max 10,000 precursors per sample for the box plot view:

```r
# If precursor matrix is very large, subsample for box plot
if (nrow(pre_mat) > 10000) {
  sample_idx <- sample(nrow(pre_mat), 10000)
  pre_mat_plot <- pre_mat[sample_idx, ]
} else {
  pre_mat_plot <- pre_mat
}
```
