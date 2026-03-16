# RUN COMPARATOR — MOFA2 ADDENDUM

> **Extends:** RUN_COMPARATOR_SPEC_V2.md
> **Placement:** Optional step in the AI Analysis tab of the Run Comparator module
> **Prerequisite:** MOFA2_INTEGRATION_SPEC.md infrastructure (server_mofa.R)
> **Availability:** Modes A, B, C1 only (requires quantification in both runs)

---

## Rationale

The Run Comparator's rule-based hypothesis engine diagnoses discordant proteins
individually — it explains *why a specific protein disagrees* (missing values,
peptide count, normalization offset, etc.). What it cannot do is answer:

> *"Do the discordant proteins share a hidden pattern that a single rule can't see?"*

Examples of patterns MOFA2 can detect that the hypothesis engine cannot:
- All discordant proteins load heavily on a single latent factor that is strong
  in Run A but suppressed in Run B → Run B's normalization is collapsing a
  real biological signal
- Discordant proteins split into two groups with opposite factor loadings →
  two different mechanisms are causing disagreement simultaneously
- A batch-like factor appears only in Run B → a normalization artifact in
  Spectronaut or FragPipe that DIA-NN didn't introduce
- Discordant proteins are enriched in a specific GO term or pathway →
  the disagreement is not random, it's biologically coherent

MOFA2 treats Run A and Run B as two "views" of the same samples and
decomposes the joint variance into shared and view-specific factors. This
gives a global decomposition of *why the runs differ* that complements the
per-protein rule-based diagnoses.

---

## Where It Lives

The MOFA2 decomposition is an **optional step** in the **AI Analysis tab**
of the Run Comparator, positioned between the Gemini summary and the Claude
ZIP export:

```
AI Analysis tab layout:
  ┌─────────────────────────────────────────────┐
  │ [Generate Gemini Summary]                   │
  │ <Gemini response card>                      │
  │                                             │
  │ ─── Optional: MOFA2 Decomposition ──────── │
  │ [Run MOFA2 on Discordant Proteins]          │  ← NEW
  │ <MOFA2 results (collapsible card)>          │  ← NEW
  │                                             │
  │ ─────────────────────────────────────────── │
  │ [Export ZIP for Claude Analysis]            │
  └─────────────────────────────────────────────┘
```

It is gated behind a button so users who don't need it aren't slowed down.
MOFA2 adds ~10–30 seconds of compute depending on protein count.

---

## Data Preparation

MOFA2 runs on the **Shared & quantified** proteins only — those in Layer 2's
top tier. Only two views are used: Run A intensities and Run B intensities,
sample-matched via `comparator_sample_map`.

```r
prepare_mofa_views <- function(comparator_data) {
  shared <- comparator_data$protein_universe |>
    dplyr::filter(tier == "shared_quantified") |>
    dplyr::pull(protein_id)

  # Run A: log2-transformed intensity matrix, proteins × samples
  mat_a <- comparator_data$run_a$intensities |>
    dplyr::filter(rownames(.data) %in% shared) |>
    as.matrix() |>
    log2_safe()   # log2(x + 1), preserves NAs

  # Run B: same, using matched sample order from comparator_sample_map
  sample_order <- comparator_data$comparator_sample_map$run_b
  mat_b <- comparator_data$run_b$intensities[, sample_order, drop = FALSE] |>
    dplyr::filter(rownames(.data) %in% shared) |>
    as.matrix() |>
    log2_safe()

  # Align row names
  common_proteins <- intersect(rownames(mat_a), rownames(mat_b))
  list(
    run_a = mat_a[common_proteins, ],
    run_b = mat_b[common_proteins, ],
    n_proteins = length(common_proteins),
    n_samples  = ncol(mat_a)
  )
}

log2_safe <- function(mat) {
  mat[mat == 0] <- NA
  log2(mat)
}
```

---

## MOFA2 Model

A lightweight MOFA2 run — fewer factors than the full MOFA tab, focused on
diagnostic decomposition rather than exploratory analysis:

```r
run_comparator_mofa <- function(views, n_factors = 5) {
  mofa_obj <- MOFA2::create_mofa(list(
    run_a = views$run_a,
    run_b = views$run_b
  ))

  # Conservative settings for diagnostic use — fast convergence
  data_opts           <- MOFA2::get_default_data_options(mofa_obj)
  data_opts$scale_views <- TRUE   # important: views may have different scales

  model_opts           <- MOFA2::get_default_model_options(mofa_obj)
  model_opts$num_factors <- n_factors

  train_opts           <- MOFA2::get_default_training_options(mofa_obj)
  train_opts$maxiter   <- 500     # sufficient for diagnostic use
  train_opts$seed      <- 42
  train_opts$verbose   <- FALSE

  mofa_obj <- MOFA2::prepare_mofa(mofa_obj,
    data_options     = data_opts,
    model_options    = model_opts,
    training_options = train_opts
  )

  MOFA2::run_mofa(mofa_obj, use_basilisk = TRUE)
}
```

---

## Outputs — Three Focused Visualizations

### Output 1: Variance Explained Heatmap

The core diagnostic. Shows which factors are shared between runs vs
run-specific.

```r
output$comparator_mofa_variance <- renderPlot({
  req(values$comparator_mofa)
  MOFA2::plot_variance_explained(
    values$comparator_mofa,
    x = "view", y = "factor",
    plot_total = TRUE
  ) +
    labs(title = "Variance explained per run and factor",
         subtitle = paste0(
           "Shared factors (high in both) = biology both runs agree on. ",
           "Run-specific factors = where the runs diverge."
         )) +
    theme_minimal(base_size = 12)
})
```

**How to read it:**
- Factor with high % in **both** views → shared biology, runs agree
- Factor with high % in **Run A only** → signal Run B is missing or suppressing
- Factor with high % in **Run B only** → artifact or signal unique to Run B's pipeline

### Output 2: Discordant Protein Factor Weights

Shows where the discordant proteins sit in factor space. Do they cluster
on a single factor, or are they scattered?

```r
output$comparator_mofa_weights <- renderPlotly({
  req(values$comparator_mofa, comparator_data())

  discordant_ids <- comparator_data()$discordant_table$protein_id

  # Get top 2 factors by variance explained in Run A
  var_exp    <- MOFA2::get_variance_explained(values$comparator_mofa)
  top_factors <- names(sort(
    var_exp$r2_per_factor[[1]][, "run_a"],
    decreasing = TRUE
  ))[1:2]

  weights <- MOFA2::get_weights(values$comparator_mofa,
                                 views = "run_a",
                                 factors = top_factors,
                                 as.data.frame = TRUE) |>
    dplyr::rename(protein_id = feature) |>
    tidyr::pivot_wider(names_from = factor, values_from = value)

  weights$is_discordant <- weights$protein_id %in% discordant_ids
  weights$label <- ifelse(weights$is_discordant,
                          weights$protein_id, "")

  plotly::plot_ly(weights,
    x = ~get(top_factors[1]),
    y = ~get(top_factors[2]),
    color = ~is_discordant,
    colors = c("grey80", "#e74c3c"),
    text  = ~label,
    type  = "scatter", mode = "markers",
    marker = list(size = 6, opacity = 0.7)
  ) |>
    plotly::layout(
      xaxis = list(title = top_factors[1]),
      yaxis = list(title = top_factors[2]),
      title = "Factor weights — discordant proteins highlighted"
    )
})
```

### Output 3: Top Factor Interpretation Summary

A small table showing the top 10 proteins by absolute weight on the most
run-A-specific factor (the factor most likely responsible for discordance):

```r
output$comparator_mofa_top_weights <- renderDT({
  req(values$comparator_mofa, comparator_data())

  # Identify most run-A-specific factor
  var_exp    <- MOFA2::get_variance_explained(values$comparator_mofa)
  r2_a       <- var_exp$r2_per_factor[[1]][, "run_a"]
  r2_b       <- var_exp$r2_per_factor[[1]][, "run_b"]
  specificity <- r2_a - r2_b   # positive = more explained in run_a than run_b
  key_factor  <- names(which.max(specificity))

  discordant_ids <- comparator_data()$discordant_table$protein_id

  MOFA2::get_weights(values$comparator_mofa,
                     views = "run_a",
                     factors = key_factor,
                     as.data.frame = TRUE) |>
    dplyr::rename(protein_id = feature, weight = value) |>
    dplyr::mutate(
      abs_weight    = abs(weight),
      is_discordant = protein_id %in% discordant_ids
    ) |>
    dplyr::arrange(dplyr::desc(abs_weight)) |>
    dplyr::slice_head(n = 15) |>
    datatable(
      caption = paste0("Top proteins driving ", key_factor,
                       " (most run-A-specific factor)"),
      rownames = FALSE,
      options  = list(pageLength = 15, dom = 't'),
      class    = "compact"
    ) |>
    formatStyle("is_discordant",
      target          = "row",
      backgroundColor = styleEqual(TRUE, "#fde8e8")
    )
})
```

Discordant proteins highlighted in pink. If many of the top weight proteins
are discordant, the factor is likely *the* cause of disagreement.

---

## UI — Collapsible Card in AI Analysis Tab

```r
# Inside the AI Analysis nav_panel, after Gemini card:
hr(),
h5("Optional: MOFA2 Factor Decomposition",
   class = "text-muted mt-3"),
p(class = "text-muted small",
  "Treats Run A and Run B as two views of the same samples and decomposes ",
  "joint variance into shared and run-specific factors. Helps identify whether ",
  "discordant proteins share a hidden pattern (e.g., all driven by the same ",
  "latent factor that one run's pipeline is suppressing)."),

conditionalPanel(
  "input.comparator_mofa_btn == 0",
  actionButton("comparator_mofa_btn",
               "Run MOFA2 Decomposition (~15–30 sec)",
               icon = icon("circle-nodes"),
               class = "btn-outline-secondary")
),

conditionalPanel(
  "input.comparator_mofa_btn > 0",
  uiOutput("comparator_mofa_status"),   # "Running..." spinner or "Complete"

  bslib::accordion(
    open = TRUE,
    bslib::accordion_panel("Variance Explained by Run and Factor",
      plotOutput("comparator_mofa_variance", height = "320px")
    ),
    bslib::accordion_panel("Discordant Proteins in Factor Space",
      plotlyOutput("comparator_mofa_weights", height = "380px")
    ),
    bslib::accordion_panel("Top Proteins on Most Run-Specific Factor",
      DTOutput("comparator_mofa_top_weights")
    )
  )
)
```

---

## Reactive Logic

```r
# Add to app.R reactiveValues():
comparator_mofa = NULL   # MOFA2 model object for comparator

# In server_comparator.R:
observeEvent(input$comparator_mofa_btn, {
  req(comparator_data())

  views <- prepare_mofa_views(comparator_data())

  # Guard: need at least 10 shared proteins and 4 samples
  if (views$n_proteins < 10) {
    showNotification("Too few shared proteins for MOFA2 (need ≥ 10).",
                     type = "warning")
    return()
  }
  if (views$n_samples < 4) {
    showNotification("Too few samples for MOFA2 (need ≥ 4).",
                     type = "warning")
    return()
  }

  withProgress(message = "Running MOFA2 decomposition...", value = 0.3, {
    tryCatch({
      mofa_obj <- run_comparator_mofa(views, n_factors = min(5, views$n_samples - 1))
      setProgress(0.9, message = "Extracting results...")
      values$comparator_mofa <- mofa_obj
      setProgress(1.0, message = "Done")
    }, error = function(e) {
      showNotification(paste("MOFA2 failed:", conditionMessage(e)),
                       type = "error", duration = 10)
    })
  })
})
```

---

## Gemini Prompt Enhancement

When MOFA2 has been run, the Gemini prompt is automatically enriched with
the factor structure. The `build_gemini_comparator_prompt()` function gains
an optional `mofa_summary` argument:

```r
build_gemini_comparator_prompt <- function(comparator_data,
                                            mofa_obj = NULL) {
  # ... existing prompt content ...

  mofa_section <- if (!is.null(mofa_obj)) {
    var_exp    <- MOFA2::get_variance_explained(mofa_obj)
    r2_a       <- var_exp$r2_per_factor[[1]][, "run_a"]
    r2_b       <- var_exp$r2_per_factor[[1]][, "run_b"]
    specificity <- r2_a - r2_b
    key_factor  <- names(which.max(specificity))

    top_weights <- MOFA2::get_weights(mofa_obj, views = "run_a",
                                       factors = key_factor,
                                       as.data.frame = TRUE) |>
      dplyr::arrange(dplyr::desc(abs(value))) |>
      dplyr::slice_head(n = 5) |>
      dplyr::pull(feature) |>
      paste(collapse = ", ")

    glue::glue("
MOFA2 FACTOR DECOMPOSITION:
- {key_factor} is the most run-A-specific factor
  (explains {round(r2_a[key_factor]*100, 1)}% variance in Run A,
   {round(r2_b[key_factor]*100, 1)}% in Run B)
- Top proteins driving this factor: {top_weights}
- This suggests Run B's pipeline is suppressing the signal captured by
  this factor. Check whether the discordant proteins overlap with these
  top-weight proteins.
")
  } else ""

  glue::glue("
{existing_prompt_content}
{mofa_section}
")
}
```

This means Gemini sees not just per-protein patterns but the global factor
structure — giving it much better context to identify the root cause.

---

## Claude ZIP Enhancement

When MOFA2 has been run, two additional files are added to the ZIP:

```
plots/
  ├── mofa_variance_explained.png   # NEW
  └── mofa_factor_weights.png       # NEW
```

And `claude_prompt.md` gains a MOFA section with the same factor summary
text included in the Gemini prompt.

---

## Availability Guard

MOFA2 decomposition is only available when:
1. Mode is A, B, or C1 (Run B has quantification data)
2. ≥ 10 shared & quantified proteins
3. ≥ 4 matched samples
4. MOFA2 package is installed (check via `requireNamespace("MOFA2")`)

If MOFA2 is not installed, replace the button with:
```r
div(class = "alert alert-secondary",
  icon("circle-info"), " ",
  "MOFA2 is not installed. To enable this feature, run: ",
  tags$code('BiocManager::install("MOFA2")')
)
```

---

## Additional Dependencies

| Package | Use | Status |
|---------|-----|--------|
| `MOFA2` | Factor decomposition | Already in MOFA2_INTEGRATION_SPEC — check if installed |
| `basilisk` | MOFA2 Python backend | Installed with MOFA2 |

No new dependencies beyond what MOFA2_INTEGRATION_SPEC already requires.

---

## Testing Checklist

- [ ] Button appears after comparison runs (not before)
- [ ] Button hidden in Mode C2 (no Run B quantification)
- [ ] Spinner shows during MOFA2 run, disappears on completion
- [ ] Variance heatmap renders — 2 rows (run_a, run_b), N factor columns
- [ ] Factor weights scatter: discordant proteins highlighted in red
- [ ] Top weights table: discordant proteins highlighted in pink
- [ ] Guard: < 10 shared proteins → warning notification, no crash
- [ ] Guard: < 4 samples → warning notification, no crash
- [ ] Guard: MOFA2 not installed → informational message with install command
- [ ] Gemini prompt enriched with MOFA section when run before Gemini
- [ ] Gemini prompt NOT enriched (no MOFA section) when run before MOFA2
- [ ] ZIP includes `mofa_variance_explained.png` when MOFA2 was run
- [ ] ZIP does NOT include MOFA plots when MOFA2 was not run
- [ ] Session save: `comparator_mofa` object persisted correctly
- [ ] Session load: MOFA plots restore from saved object

---

## CLAUDE.md Updates

After implementation, add:
- [ ] Run Comparator: optional MOFA2 factor decomposition in AI Analysis tab ✅
- [ ] `prepare_mofa_views()` — prepares shared-protein two-view matrix ✅
- [ ] `run_comparator_mofa()` — lightweight 5-factor diagnostic MOFA run ✅
- [ ] Gemini prompt auto-enriched with factor structure when MOFA2 is run ✅
