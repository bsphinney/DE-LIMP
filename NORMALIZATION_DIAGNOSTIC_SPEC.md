# Feature Spec: Normalization Diagnostic Plot

## Overview

Add a diagnostic visualization to the **QC Plots** tab that shows per-sample intensity distributions
at two stages of the pipeline: (1) precursor-level input from DIA-NN, and (2) protein-level output
from DPC-Quant. This helps users confirm that:

- DIA-NN normalization produced well-aligned precursor distributions across samples
- Protein-level quantification preserved that alignment
- No individual samples are behaving as outliers

**Important**: DIA-NN's RT-dependent normalization is the **only** cross-run normalization in this
pipeline. Limpa's `dpcCN()` estimates the missingness model and `dpcQuant()` performs protein
aggregation — neither applies normalization. This plot is therefore primarily a check on DIA-NN's
normalization quality, with a secondary check that protein quantification didn't introduce artifacts.

## Background: The Pipeline — Normalization, Missingness Modeling, and Quantification

Understanding what each step does (and doesn't do) is critical for this feature.

### Step 1: DIA-NN Normalization (Pre-DE-LIMP) — THE ONLY NORMALIZATION
- **Default**: RT-dependent normalization (recommended, always-on unless user disabled it)
- **How it works**: Selects top N precursors with lowest CVs, normalizes their summed intensities
  to be equal across runs, then applies RT-bin-specific correction factors (running median of
  fold changes across retention time bins)
- **Applied at**: Precursor level, before DE-LIMP ever sees the data
- **Key columns in parquet file**:
  - `Precursor.Quantity` = raw, unnormalized intensity
  - `Precursor.Normalised` = after DIA-NN normalization
  - Normalization factor = `Precursor.Normalised / Precursor.Quantity`
- **Critical**: `limpa::readDIANN()` reads `Precursor.Normalised` by default (its `qty.column`
  parameter defaults to `"Precursor.Normalised"`)
- **This is the only cross-run normalization in the entire pipeline.** Limpa does NOT add its
  own normalization layer. If DIA-NN normalization was off, the data flows through the rest of
  the pipeline unnormalized.

### Step 2: limpa `dpcCN()` — MISSINGNESS MODEL (NOT normalization)
- **What it is**: Estimates the Detection Probability Curve (DPC) under the Complete Normal model
- **What "CN" stands for**: "Complete Normal" — the statistical model assumption that the
  *complete* (observed + unobserved) log-intensities are normally distributed. **NOT** "Cyclic
  Normalization" as previously mislabeled.
- **What it does**: Learns the relationship between a peptide's intensity and its probability of
  being detected. Low-intensity peptides are more likely to be missing — the DPC quantifies
  this precisely with an intercept and slope parameter.
- **What it does NOT do**: It does NOT adjust, shift, scale, or transform the intensity values.
  The precursor intensities go in and come out unchanged. It only estimates the DPC parameters
  (stored in `values$dpc_fit`).
- **Input**: `values$raw_data` (EList from `readDIANN()`)
- **Output**: DPC fit object (`values$dpc_fit`) containing `$dpc` (intercept + slope)

### Step 3: limpa `dpcQuant()` — PROTEIN QUANTIFICATION (NOT normalization)
- **What it does**: Aggregates precursor-level data → protein-level expression estimates
- Uses the DPC from Step 2 to probabilistically recover information from missing values
- Fits an additive model (protein expression + peptide baseline) for each protein
- Missing values are represented by the DPC rather than imputed with point estimates
- Empirical Bayes prior borrows information across all proteins
- **Output**: `values$y_protein` — protein-level EList with COMPLETE data (no NAs), plus
  standard errors for each protein-in-sample estimate
- **What it does NOT do**: No cross-run normalization. The protein expression estimates
  reflect whatever normalization state the input precursor data was in.

### Step 4: limpa `dpcDE()` — DIFFERENTIAL EXPRESSION
- Wraps limma's linear model framework with precision weights from DPC-Quant standard errors
- Proteins with more missing precursors → larger standard errors → downweighted in DE analysis
- This is where the uncertainty from missingness is properly propagated

### What this means for the diagnostic

`values$raw_data$E` is **precursor-level, log2, DIA-NN normalized** data (with NAs).
`values$y_protein$E` is **protein-level, log2, DPC-Quant aggregated** data (complete, no NAs).

The difference between these two is driven by **protein aggregation and missing value recovery**,
NOT by a second normalization step. Per-sample distribution changes (median shifts, shape changes)
between the two panels reflect the effect of rolling up precursors to proteins and filling in missing
values — not normalization correction.

This means the diagnostic is really showing two things:
1. Whether DIA-NN's normalization produced well-aligned precursor distributions (left panel)
2. Whether the protein-level quantification preserved that alignment (right panel)

If the left panel already looks bad (scattered medians), the right panel won't magically fix it —
limpa trusts the input normalization.

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
        span("Pipeline Diagnostic: Input → Output Distributions"),
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
    Stage = "Precursor Input\n(DIA-NN normalized, with NAs)"
  )

  # --- Post-quantification: per-sample medians from protein data ---
  post_mat <- values$y_protein$E  # protein-level, log2, no NAs
  post_medians <- data.frame(
    Sample = colnames(post_mat),
    Median = apply(post_mat, 2, median, na.rm = TRUE),
    Stage = "Protein Output\n(DPC-Quant aggregated, complete)"
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
      mutate(Stage = "Precursor Input\n(DIA-NN normalized)") %>%
      filter(!is.na(Log2Intensity))

    post_long <- as.data.frame(post_mat) %>%
      pivot_longer(everything(), names_to = "Sample", values_to = "Log2Intensity") %>%
      mutate(Stage = "Protein Output\n(DPC-Quant)")

    # Add group info to both
    pre_long$Group <- meta$Group[match(pre_long$Sample, meta$File.Name)]
    post_long$Group <- meta$Group[match(post_long$Sample, meta$File.Name)]

    combined <- bind_rows(pre_long, post_long)
    combined$Stage <- factor(combined$Stage, levels = c("Precursor Input\n(DIA-NN normalized)", "Protein Output\n(DPC-Quant)"))

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
        title = "Pipeline Diagnostic: Precursor Input → Protein Output",
        subtitle = "Left: DIA-NN normalized precursors | Right: DPC-Quant protein estimates",
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
    title = "Pipeline Diagnostic - Fullscreen View",
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

## User-Facing Guidance System

The diagnostic plot is only useful if users can interpret it and know what to do when something
looks wrong. This section specifies the automated guidance that should appear contextually.

### 1. Automatic Health Assessment

After the pipeline runs, compute a simple diagnostic score and display a status indicator on the
plot card header (next to the DIA-NN badge):

```r
# Compute after pipeline completes
assess_distribution_health <- reactive({
  req(values$raw_data, values$y_protein)

  pre_mat <- values$raw_data$E
  post_mat <- values$y_protein$E

  # Per-sample medians
  pre_medians <- apply(pre_mat, 2, median, na.rm = TRUE)
  post_medians <- apply(post_mat, 2, median, na.rm = TRUE)

  # Metric: coefficient of variation of medians (lower = better aligned)
  pre_cv <- sd(pre_medians) / abs(mean(pre_medians))
  post_cv <- sd(post_medians) / abs(mean(post_medians))

  # Metric: any single sample > 2 SD from mean of medians?
  pre_outliers <- which(abs(pre_medians - mean(pre_medians)) > 2 * sd(pre_medians))
  post_outliers <- which(abs(post_medians - mean(post_medians)) > 2 * sd(post_medians))

  list(
    pre_cv = pre_cv,
    post_cv = post_cv,
    pre_outlier_samples = names(pre_outliers),
    post_outlier_samples = names(post_outliers),
    status = if (pre_cv > 0.05) "warning" else if (length(pre_outliers) > 0) "caution" else "good"
  )
})
```

### 2. Contextual Warning Banners

Display above the plot, color-coded by severity. These fire automatically based on the health
assessment and DIA-NN normalization detection.

```r
output$norm_diag_guidance <- renderUI({
  req(values$raw_data, values$y_protein)

  health <- assess_distribution_health()
  diann_status <- values$diann_norm_detected

  warnings <- list()

  # --- SCENARIO 1: DIA-NN normalization OFF + bad distributions ---
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
        "isotope labeling time-courses, unnormalized data may be appropriate — ",
        "but you should apply your own normalization before using DE-LIMP."
      )
    ))
  }

  # --- SCENARIO 2: DIA-NN normalization OFF but distributions look OK ---
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

  # --- SCENARIO 3: DIA-NN normalization ON but distributions still bad ---
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

  # --- SCENARIO 4: Outlier sample(s) detected ---
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

  # --- SCENARIO 5: Everything looks good ---
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
```

Place `uiOutput("norm_diag_guidance")` **above** the plot in the card body.

### 3. "What does this mean?" Expandable Section

Below the plot, add a collapsible explanation written for bench scientists:

```r
tags$details(
  tags$summary(style = "cursor: pointer; color: #0d6efd; font-size: 0.9em;",
    icon("question-circle"), " What am I looking at?"
  ),
  div(style = "background-color: #f8f9fa; padding: 12px; border-radius: 5px;
               margin-top: 8px; font-size: 0.85em; line-height: 1.6;",

    tags$h6("Reading this plot"),
    p("Each box (or density curve) represents one sample's intensity distribution — ",
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
```

### 4. Updated Card Body Layout

The full card body should be structured as:

```r
card_body(
  # Contextual warnings (auto-generated)
  uiOutput("norm_diag_guidance"),

  # View toggle
  radioButtons("norm_diag_type", "View:",
    choices = c("Box Plots" = "boxplot", "Density Overlay" = "density"),
    inline = TRUE
  ),

  # The plot
  plotlyOutput("norm_diagnostic_plot", height = "450px"),

  # Expandable explanation
  tags$details(
    tags$summary(style = "cursor: pointer; color: #0d6efd; font-size: 0.9em; margin-top: 10px;",
      icon("question-circle"), " What am I looking at?"
    ),
    # ... (explanation content from section 3 above)
  )
)
```

## What Developers Should Know (Internal Reference)

These are the expected outcomes for each data scenario. Use for testing and QA.

### Healthy Data (DIA-NN Normalization ON — most common case)
- **Left panel**: Box plot medians well-aligned across samples
- **Right panel**: Medians remain aligned, distributions slightly narrower (fewer features, no NAs)
- **DIA-NN badge**: Blue "ON"
- **Automated guidance**: Green "Distributions look good" banner

### Data with DIA-NN Normalization OFF
- **Left panel**: Medians may be scattered depending on sample loading variation
- **Right panel**: Medians will still be scattered — limpa does not normalize
- **DIA-NN badge**: Yellow "OFF"
- **Automated guidance**: If medians scattered → yellow warning with re-processing instructions.
  If medians happen to be aligned → blue info note.

### Problematic Data (even with normalization ON)
- One sample diverges in both panels → red banner naming the outlier sample
- All medians scattered despite norm ON → red banner about possible batch effects
- Right panel diverges but left looks fine → unusual, may indicate extreme missingness in some samples

## Files to Modify

1. **DE-LIMP.R** (and then copy to app.R):
   - Add `diann_norm_detected` to `reactiveValues()` initialization (~line 555)
   - Add DIA-NN norm detection code to both data loading observers (~lines 599 and 640)
   - Add UI elements to `nav_panel("QC Plots")` (~line 425)
   - Add `assess_distribution_health` reactive (~after pipeline runs)
   - Add `output$norm_diag_guidance` renderUI
   - Add plot reactive + render functions (after existing QC renderers ~line 1000)

2. **Dockerfile**: No changes needed (no new packages required — uses existing ggplot2, plotly, tidyr)

3. **CLAUDE.md**: Add to Recent Changes section documenting the new feature

4. **USER_GUIDE.md**: Add section about interpreting the diagnostic plot

## Testing Checklist

- [ ] Plot renders correctly with example data (which has DIA-NN norm OFF)
- [ ] Plot renders correctly with user-uploaded data (DIA-NN norm likely ON)
- [ ] DIA-NN badge correctly detects ON vs OFF vs unknown
- [ ] Box plot view works and medians are visible
- [ ] Density overlay view works and curves are distinguishable
- [ ] Fullscreen modal works
- [ ] **Green banner appears when distributions are healthy**
- [ ] **Yellow warning appears when DIA-NN norm is OFF + distributions scattered**
- [ ] **Red banner appears when distributions bad despite normalization ON**
- [ ] **Outlier sample banner names the correct sample(s)**
- [ ] "What am I looking at?" section expands/collapses correctly
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
