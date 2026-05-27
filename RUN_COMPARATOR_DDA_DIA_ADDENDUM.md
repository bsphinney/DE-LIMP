# Run Comparator Addendum: Mode D — DIA vs DDA (Same Samples)

> **Addendum to**: `RUN_COMPARATOR_SPEC_V2.md`
> **Version**: 1.0 — March 2026
> **Author**: Brett Phinney / UC Davis Proteomics Core
> **Prereqs**: `DDA_SAGE_CASANOVO_SPEC.md` implemented (§§1–16);
>   `RUN_COMPARATOR_SPEC_V2.md` Mode A–C implemented

---

## Overview

Mode D compares a **DE-LIMP DIA run** (DIA-NN + limpa + dpcQuant) against a
**DE-LIMP DDA run** (Sage + MaxLFQ + cyclic loess + Perseus imputation) on the
same biological samples. Unlike Modes A–C, both runs live *inside* DE-LIMP —
no external file import needed. The session already holds both ELists.

This is the most structurally informative comparison DE-LIMP can offer: same
samples, same biology, two fundamentally different acquisition strategies and
quantification pipelines. The diagnostic output explains *why* results differ
in terms users can act on.

---

## When Mode D is Available

Mode D appears in the comparator mode selector **only when both are present
in the current session:**

```r
# Availability check — in server_comparator.R
mode_d_available <- reactive({
  !is.null(values$y_protein) &&          # DIA pipeline has run
  !is.null(values$dda$elist) &&           # DDA pipeline has run
  !is.null(values$fit) &&                 # DIA DE fit exists
  !is.null(values$dda$de_fit)             # DDA DE fit exists
})
```

The mode selector in the comparator UI gains a fourth option, conditionally shown:

```r
observe({
  choices <- c(
    "DE-LIMP vs DE-LIMP"     = "delimp_delimp",
    "DE-LIMP vs Spectronaut" = "delimp_spectronaut",
    "DE-LIMP vs FragPipe"    = "delimp_fragpipe"
  )
  if (isTRUE(mode_d_available())) {
    choices <- c(choices,
      "DIA vs DDA (same samples)" = "dia_vs_dda"
    )
  }
  updateRadioButtons(session, "comparator_mode", choices = choices)
})
```

---

## Architecture

Mode D bypasses file parsing entirely — both runs are already in reactive memory.
The comparator reads directly from `values$*` (DIA) and `values$dda$*` (DDA).

```
DIA pipeline (values$*)          DDA pipeline (values$dda$*)
  y_protein$E                      dda$elist$E
  y_protein$genes$NPeptides        dda$protein_meta$NPeptides
  fit (limma MArrayLM)             dda$de_fit (limma MArrayLM)
  metadata                         metadata (shared)
  norm: DIA-NN RT-dependent        norm: cyclic loess (or user choice)
  rollup: dpcQuant                 rollup: Sage MaxLFQ
  imputation: dpcCN (model)        imputation: Perseus (user choice)
        │                                    │
        └──────── comparator_data() ─────────┘
                  (Mode D parser — no I/O)
```

---

## Mode D Parser — `parse_mode_d()`

No file I/O. Pulls directly from reactive values.

```r
#' Build comparator-compatible data structures from live DIA + DDA sessions
#'
#' @param values Shiny reactiveValues — must contain both DIA and DDA results
#' @param dia_contrast Contrast name to use from DIA fit (string)
#' @param dda_contrast Contrast name to use from DDA fit (string)
#' @return Named list matching the structure returned by Modes A–C parsers,
#'   with source = "dia_delimp" (Run A) and "dda_delimp" (Run B)
parse_mode_d <- function(values, dia_contrast, dda_contrast) {

  # ── DIA (Run A) ────────────────────────────────────────────────────────
  dia_tt <- limma::topTable(values$fit, coef = dia_contrast, number = Inf,
                            sort.by = "none")
  dia_tt$protein_id <- normalize_protein_id(rownames(dia_tt))

  run_a <- list(
    source      = "dia_delimp",
    de_stats    = dia_tt,
    intensities = as.data.frame(values$y_protein$E),
    n_peptides  = setNames(
      values$y_protein$genes$NPeptides,
      normalize_protein_id(rownames(values$y_protein$genes))
    ),
    missing_pct = apply(values$y_protein$E, 1,
                        function(x) mean(is.na(x))),
    prop_obs    = values$y_protein$genes$PropObs,   # DIA-only: DPC detection rate
    settings    = list(
      software        = "DE-LIMP (DIA)",
      acquisition     = "DIA (diaPASEF)",
      search_engine   = glue::glue("DIA-NN {values$diann_version %||% 'unknown'}"),
      normalization   = "DIA-NN RT-dependent (Precursor.Normalised)",
      rollup_method   = "DPC-Quant (empirical Bayes precursor aggregation)",
      imputation      = "dpcCN probabilistic model (not Perseus-style)",
      de_engine       = "limma moderated t-test",
      fdr_threshold   = values$fdr_threshold %||% "unknown",
      lfc_threshold   = values$lfc_threshold %||% "unknown",
      fasta           = values$fasta_path    %||% "unknown",
      contrast        = dia_contrast
    )
  )

  # ── DDA (Run B) ────────────────────────────────────────────────────────
  dda_tt <- limma::topTable(values$dda$de_fit, coef = dda_contrast,
                            number = Inf, sort.by = "none")
  dda_tt$protein_id <- normalize_protein_id(rownames(dda_tt))

  run_b <- list(
    source      = "dda_delimp",
    de_stats    = dda_tt,
    intensities = as.data.frame(values$dda$elist$E),
    n_peptides  = setNames(
      values$dda$protein_meta$NPeptides,
      normalize_protein_id(values$dda$protein_meta$ProteinID)
    ),
    missing_pct = values$dda$imputation_frac,    # pre-imputation NA fraction (§16.10)
    impute_frac = values$dda$imputation_frac,    # same field, DDA-specific label
    settings    = list(
      software        = "DE-LIMP (DDA)",
      acquisition     = "DDA (ddaPASEF)",
      search_engine   = glue::glue("Sage {values$dda$search_params$sage_version %||% 'unknown'}"),
      normalization   = values$dda$norm_method   %||% "cyclicloess",
      rollup_method   = "Sage MaxLFQ (consensus ratio protein rollup)",
      imputation      = glue::glue(
        "Perseus-style (width={values$dda$search_params$perseus_width %||% 0.3}, ",
        "shift={values$dda$search_params$perseus_shift %||% 1.8})"
      ),
      de_engine       = "limma moderated t-test",
      fdr_threshold   = values$fdr_threshold %||% "unknown",   # shared threshold
      lfc_threshold   = values$lfc_threshold %||% "unknown",
      valid_filter    = glue::glue(
        "{values$dda$search_params$min_valid * 100 %||% 50}% per group"
      ),
      fasta           = values$dda$search_params$fasta_path %||% "unknown",
      contrast        = dda_contrast
    )
  )

  list(run_a = run_a, run_b = run_b, mode = "dia_vs_dda")
}
```

---

## Contrast Matching for Mode D

Both DIA and DDA fits can have multiple contrasts. The selector shows
dropdowns for each independently and warns when names don't match:

```r
output$comparator_contrast_selectors <- renderUI({
  req(input$comparator_mode == "dia_vs_dda")
  req(values$fit, values$dda$de_fit)

  dia_contrasts <- colnames(values$fit$contrasts)
  dda_contrasts <- colnames(values$dda$de_fit$contrasts)

  tagList(
    div(class = "row",
      div(class = "col-6",
        selectInput("comp_dia_contrast", "DIA contrast",
          choices = dia_contrasts, selected = dia_contrasts[1])
      ),
      div(class = "col-6",
        selectInput("comp_dda_contrast", "DDA contrast",
          choices = dda_contrasts,
          selected = {
            # Auto-match if names are identical
            if (dia_contrasts[1] %in% dda_contrasts) dia_contrasts[1]
            else dda_contrasts[1]
          }
        )
      )
    ),
    # Warn if auto-match failed
    if (!identical(input$comp_dia_contrast, input$comp_dda_contrast)) {
      div(class = "alert alert-warning mt-1",
        icon("triangle-exclamation"), " ",
        "DIA and DDA contrast names differ. Confirm you are comparing ",
        "the same biological question before proceeding."
      )
    }
  )
})
```

---

## Layer 1 — Settings Diff (Mode D)

The settings diff in Mode D is more structured than external-tool modes because
we have complete metadata for both runs. Key rows added specifically for DIA↔DDA:

```r
build_settings_diff_mode_d <- function(run_a, run_b) {
  # All rows from existing build_settings_diff(), plus DDA-specific additions:
  dda_specific <- data.frame(
    parameter = c(
      "Acquisition mode",
      "Search engine",
      "Normalization",
      "Protein rollup",
      "Missingness model",
      "Imputation method",
      "Valid value filter",
      "FASTA used"
    ),
    run_a_value = c(
      run_a$settings$acquisition,
      run_a$settings$search_engine,
      run_a$settings$normalization,
      run_a$settings$rollup_method,
      run_a$settings$imputation,
      "dpcCN (Detection Probability Curve model)",
      "N/A — dpcQuant handles missingness probabilistically",
      run_a$settings$fasta
    ),
    run_b_value = c(
      run_b$settings$acquisition,
      run_b$settings$search_engine,
      run_b$settings$normalization,
      run_b$settings$rollup_method,
      run_b$settings$imputation,
      run_b$settings$imputation,
      run_b$settings$valid_filter,
      run_b$settings$fasta
    ),
    match = c(
      "differs",   # always — DIA vs DDA
      "differs",   # always — DIA-NN vs Sage
      "differs",   # always — different norm strategies
      "differs",   # always — dpcQuant vs MaxLFQ
      "differs",   # always — model vs Perseus
      "differs",
      "differs",
      # FASTA: match only if same path
      ifelse(run_a$settings$fasta == run_b$settings$fasta, "match", "differs")
    )
  )

  # Note: FDR and logFC thresholds should match (shared from DE-LIMP UI)
  # Only flag as "differs" if user changed them between runs
  dda_specific
}
```

---

## Layer 2 — Protein Universe (Mode D)

Same UpSet structure as Modes A–C. For DIA↔DDA the tiers have specific meaning:

| Tier | Biological interpretation |
|------|--------------------------|
| **Shared & quantified** | Robustly detected by both acquisition modes — highest confidence |
| **DIA only** | Low-abundance proteins rescued by DIA's non-stochastic sampling; OR proteins where DDA precursor competition/co-elution caused dropout |
| **DDA only** | Proteins where DIA spectral library had poor coverage; OR proteins in a variable window gap |
| **Missingness differs** | Detected in both but quantified in different subsets of samples |

Add a Mode D-specific interpretation card below the UpSet plot:

```r
# Conditional interpretation card for Mode D
output$comparator_mode_d_interpretation <- renderUI({
  req(input$comparator_mode == "dia_vs_dda")
  comp <- comparator_data()
  req(comp)

  n_dia_only <- sum(comp$protein_universe$in_a & !comp$protein_universe$in_b)
  n_dda_only <- sum(!comp$protein_universe$in_a & comp$protein_universe$in_b)
  n_shared   <- sum(comp$protein_universe$in_a & comp$protein_universe$in_b)
  pct_shared <- round(100 * n_shared / (n_dia_only + n_dda_only + n_shared))

  bslib::card(
    bslib::card_header("Acquisition mode comparison"),
    bslib::card_body(
      p(glue::glue(
        "{pct_shared}% of detected proteins are shared between DIA and DDA. ",
        "{n_dia_only} proteins were detected only by DIA; ",
        "{n_dda_only} only by DDA."
      )),
      p(tags$strong("DIA-only proteins "), "are typically low-abundance
        candidates rescued by DIA's exhaustive fragmentation. They are not
        missing from DDA due to absence — they were below the ddaPASEF
        precursor selection threshold."),
      p(tags$strong("DDA-only proteins "), "may reflect gaps in the DIA
        spectral library, proteins excluded by the DIA isolation window scheme,
        or proteins with unusual fragmentation patterns that the DIA-NN model
        scores poorly.")
    )
  )
})
```

---

## Layer 3 — Quantification (Mode D)

Standard Layer 3 outputs apply. Mode D adds one extra panel:

**3d. Peptide evidence comparison (DIA vs DDA)**

For shared proteins, scatter of:
- X: `NPeptides` in DIA (precursor count from dpcQuant)
- Y: `NPeptides` in DDA (unique peptides from Sage PSMs)
- Color: concordance tier from Layer 4

Interpretation: proteins far above the diagonal have more DDA peptide evidence
than DIA; below the diagonal have richer DIA coverage. Proteins with very few
DDA peptides but many DIA precursors are likely low-abundance candidates where
DDA sampled fewer MS2 events.

```r
output$comparator_peptide_scatter <- renderPlotly({
  req(input$comparator_mode == "dia_vs_dda")
  comp <- comparator_data()
  req(comp$shared_proteins)

  shared <- comp$shared_proteins
  plot_df <- data.frame(
    protein    = shared$protein_id,
    n_pep_dia  = comp$run_a$n_peptides[shared$protein_id],
    n_pep_dda  = comp$run_b$n_peptides[shared$protein_id],
    concordant = shared$concordant_de
  )
  plot_df <- plot_df[!is.na(plot_df$n_pep_dia) & !is.na(plot_df$n_pep_dda), ]

  plotly::plot_ly(plot_df,
    x    = ~n_pep_dia,
    y    = ~n_pep_dda,
    text = ~protein,
    color= ~concordant,
    colors = c("TRUE" = "#1d9e75", "FALSE" = "#d85a30"),
    type = "scatter", mode = "markers",
    marker = list(size = 7, opacity = 0.7)
  ) |>
    plotly::layout(
      xaxis  = list(title = "DIA precursor count (dpcQuant NPeptides)"),
      yaxis  = list(title = "DDA unique peptides (Sage PSMs)"),
      shapes = list(list(
        type = "line", x0 = 0, x1 = max(plot_df$n_pep_dia, na.rm=TRUE),
        y0   = 0, y1 = max(plot_df$n_pep_dia, na.rm=TRUE),
        line = list(dash = "dot", color = "#aaa", width = 1)
      ))
    )
})
```

---

## Layer 4 — DE Concordance (Mode D)

Full Layer 4 is available — both runs have limma fits. Mode D adds DDA-specific
columns to the discordant protein table and extends the hypothesis engine.

### Extended discordant table columns

| Column | DIA value | DDA value | Mode D addition |
|--------|-----------|-----------|-----------------|
| `logFC` | DIA logFC | DDA logFC | Both present as in Modes A–C |
| `adj.P.Val` | DIA FDR | DDA FDR | Both present |
| `n_peptides` | DIA precursor count | DDA unique peptide count | **Labeled by mode** |
| `missing_pct` | DIA missing rate (post-dpcQuant, usually ~0) | DDA pre-imputation NA fraction | **Labeled by mode** |
| `impute_flag` | N/A | ⚠ if >50% imputed (from §16.10) | **DDA only** |
| `prop_obs` | DIA PropObs (from dpcQuant) | N/A | **DIA only** |

---

## Mode D Hypothesis Engine Extensions

The existing 7-rule engine applies unchanged. Mode D adds 4 new rules
inserted at appropriate priority positions:

```r
assign_hypothesis_mode_d <- function(row, global_offset = 0) {

  # Run base engine first (rules 1–7)
  base <- assign_hypothesis(row, source_b = "dda_delimp", global_offset)

  # Rules 8–11 — Mode D specific, check after base rules
  # (only reached if base returned a low-confidence or generic result)

  imputed_heavily  <- isTRUE(row$impute_frac_b > 0.5)
  low_prop_obs_dia <- isTRUE(row$prop_obs_a < 0.4)
  few_dda_peptides <- isTRUE(row$n_peptides_b <= 2)
  many_dia_prec    <- isTRUE(row$n_peptides_a >= 6)

  # Rule 8: Heavy DDA imputation driving the result
  if (imputed_heavily && base$confidence != "High") return(list(
    hypothesis = glue::glue(
      ">{round(row$impute_frac_b * 100)}% of DDA values for this protein ",
      "were imputed (Perseus-style). The DDA DE call is largely driven by ",
      "synthetic below-detection values — treat as hypothesis-generating."
    ),
    confidence = "High"
  ))

  # Rule 9: DIA low detection rate (PropObs) but DDA has good coverage
  if (low_prop_obs_dia && !few_dda_peptides && base$confidence != "High") return(list(
    hypothesis = glue::glue(
      "DIA PropObs = {round(row$prop_obs_a, 2)} — this protein is detected ",
      "in only {round(row$prop_obs_a * 100)}% of DIA precursor slots, ",
      "meaning its DIA quantification relies heavily on dpcCN modelling. ",
      "The DDA result ({row$n_peptides_b} unique peptides) may be more reliable."
    ),
    confidence = "Medium"
  ))

  # Rule 10: Rich DIA precursor evidence but sparse DDA coverage
  if (many_dia_prec && few_dda_peptides && base$confidence != "High") return(list(
    hypothesis = glue::glue(
      "DIA has {row$n_peptides_a} precursors for this protein; DDA identified ",
      "only {row$n_peptides_b} unique peptide(s). Poor DDA sampling is the ",
      "likely cause of the quantification discrepancy — this protein may be ",
      "near the ddaPASEF detection threshold."
    ),
    confidence = "Medium"
  ))

  # Rule 11: Structural rollup difference (dpcQuant vs MaxLFQ)
  # Fires when fold-changes are in the same direction, similar magnitude,
  # but one is significant and the other isn't — and no other rule fired
  same_dir    <- sign(row$logFC_A) == sign(row$logFC_B)
  fc_similar  <- abs(row$logFC_A - row$logFC_B) < 0.5
  p_diverges  <- abs(log10(row$adjP_A + 1e-10) -
                     log10(row$adjP_B + 1e-10)) > 1
  if (same_dir && fc_similar && p_diverges) return(list(
    hypothesis = paste0(
      "Similar logFC but divergent p-values. The most likely cause is the ",
      "difference in protein rollup algorithm: DIA uses DPC-Quant (empirical ",
      "Bayes, precision-weighted) while DDA uses Sage MaxLFQ (consensus ratios). ",
      "DPC-Quant's precision weights give high-evidence proteins more statistical ",
      "power; MaxLFQ applies equal weighting per precursor. Check whether this ",
      "protein has high NPeptides in DIA — if so, DPC-Quant is likely giving it ",
      "a power advantage the DDA pipeline doesn't have."
    ),
    confidence = "Medium"
  ))

  # Fall through to base result
  base
}
```

---

## Gemini Prompt — Mode D

```r
build_gemini_comparator_prompt_mode_d <- function(stats, top_disc) {

  tool_context <- glue::glue(
    "Run A: DE-LIMP DIA pipeline (DIA-NN {stats$diann_version} → ",
    "dpcQuant rollup → limma moderated t-test). ",
    "Run B: DE-LIMP DDA pipeline (Sage {stats$sage_version} → ",
    "Sage MaxLFQ rollup → cyclic loess normalization → Perseus imputation ",
    "(width={stats$perseus_width}, shift={stats$perseus_shift}) → limma). ",
    "Both runs processed the same biological samples but with different ",
    "acquisition strategies: DIA provides exhaustive, non-stochastic ",
    "fragmentation; DDA uses data-dependent precursor selection (ddaPASEF). ",
    "\n\nKey structural differences between pipelines:\n",
    "- Normalization: DIA-NN RT-dependent (built-in) vs cyclic loess (post-MaxLFQ)\n",
    "- Missingness: DPC probabilistic model (DIA) vs filter + Perseus imputation (DDA)\n",
    "- Protein rollup: DPC-Quant empirical Bayes (DIA) vs MaxLFQ consensus ratios (DDA)\n",
    "- Coverage: DIA typically identifies more proteins at low abundance; ",
    "DDA has stochastic sampling and may miss low-abundance proteins.\n",
    "- A non-zero global intensity offset between pipelines is expected because ",
    "dpcQuant and MaxLFQ use fundamentally different scaling approaches.\n"
  )

  glue::glue("
You are analyzing a DIA vs DDA proteomics comparison in DE-LIMP.

{tool_context}

COMPARISON OVERVIEW:
- DIA proteins quantified: {stats$n_dia_proteins}
- DDA proteins quantified: {stats$n_dda_proteins}
- Shared proteins (quantified in both): {stats$n_shared}
- DIA-only: {stats$n_dia_only} | DDA-only: {stats$n_dda_only}
- Concordant DE: {stats$n_concordant_de} of {stats$n_any_de} proteins with DE in either run
- Global intensity offset (median DIA − DDA): {round(stats$global_offset, 3)} log2 units
- Per-sample correlation range: {round(stats$min_cor, 2)}–{round(stats$max_cor, 2)}

DDA PIPELINE SETTINGS:
- Normalization: {stats$dda_norm}
- Valid value filter: {stats$valid_filter}
- Imputation: Perseus (width={stats$perseus_width}, shift={stats$perseus_shift})
- DDA proteins with >50% imputed values: {stats$n_heavily_imputed}

TOP DISCORDANT PROTEINS:
{format_discord_table(top_disc)}

Please provide:
1. What fraction of the DE signal is shared between DIA and DDA, and how
   consistent are the fold-changes for concordant proteins?
2. Whether the DIA-only or DDA-only proteins suggest a specific biological
   or technical pattern (e.g., low-abundance bias, imputation artifacts)
3. For the top discordant proteins shown, the most likely cause of
   disagreement given what you know about both pipelines
4. Which pipeline's results you would prioritize for the primary
   biological conclusions, and why
5. Whether the two pipelines are sufficiently concordant to draw
   confident conclusions from either alone
")
}
```

---

## Layer Availability — Updated Matrix

Extend the existing table in `RUN_COMPARATOR_SPEC_V2.md`:

| Layer | Mode A | Mode B | Mode C1 | Mode C2 | **Mode D** |
|-------|--------|--------|---------|---------|------------|
| 1 — Settings Diff | ✅ | ✅ | ✅ | ✅ | ✅ (extended) |
| 2 — Protein Universe | ✅ | ✅ | ✅ | ✅ | ✅ (+ interpretation card) |
| 3 — Quantification | ✅ | ✅ | ✅ | ✅ | ✅ (+ peptide scatter 3d) |
| 3d — Peptide evidence scatter | ❌ | ❌ | ❌ | ❌ | ✅ |
| 4 — DE Concordance | ✅ | ✅ | ✅ | ❌ | ✅ (+ extended columns) |
| AI Gemini analysis | ✅ | ✅ | ✅ | ✅ | ✅ (Mode D prompt) |

---

## Mode D-Specific UI Additions to `ui_comparator.R`

### Mode selector extension (conditional)

```r
# Add to comparator mode radio group — shown only when mode_d_available()
conditionalPanel(
  condition = "output.mode_d_available",
  div(
    class = "alert alert-success mt-2 mb-0",
    style = "padding: 8px 12px; font-size: 12px;",
    icon("circle-check"),
    " Both DIA and DDA pipelines are complete — ",
    tags$strong("DIA vs DDA comparison available.")
  )
)
```

### Mode D contrast selectors

```r
conditionalPanel(
  condition = "input.comparator_mode === 'dia_vs_dda'",
  uiOutput("comparator_contrast_selectors")   # dynamic — see Contrast Matching §above
)
```

### Mode D UpSet interpretation card (Layer 2)

```r
conditionalPanel(
  condition = "input.comparator_mode === 'dia_vs_dda'",
  uiOutput("comparator_mode_d_interpretation")
)
```

### Mode D peptide evidence scatter (Layer 3, new tab)

Add `nav_panel` inside the existing `navset_card_tab` in the comparator results:

```r
# In the comparator results navset_card_tab, after "Quantification":
conditionalPanel(
  condition = "input.comparator_mode === 'dia_vs_dda'",
  nav_panel("Peptide Evidence",
    plotlyOutput("comparator_peptide_scatter", height = "500px"),
    div(class = "text-muted mt-2", style = "font-size: 12px;",
      "Comparing DIA precursor count (NPeptides from dpcQuant) vs ",
      "DDA unique peptide count (Sage PSMs). Proteins above the diagonal ",
      "have richer DDA evidence; below the diagonal have richer DIA coverage."
    )
  )
)
```

---

## Reactive Values — Additions to `app.R`

No new top-level reactiveValues needed. Mode D reads directly from existing
`values$*` (DIA) and `values$dda$*` (DDA) namespaces. Add only:

```r
values <- reactiveValues(
  # ... existing ...

  # Comparator: add mode_d state
  comparator_mode_d_available = FALSE    # set by observer when both pipelines complete
)

# Observer to update availability
observe({
  values$comparator_mode_d_available <- (
    !is.null(values$y_protein) &&
    !is.null(values$dda$elist) &&
    !is.null(values$fit) &&
    !is.null(values$dda$de_fit)
  )
})
```

---

## Session Save/Load — Additions

Mode D results are derived from `values$y_protein` and `values$dda$elist`,
both of which are already saved/loaded. No new session fields needed — the
comparator result is recomputed on demand when the session is reloaded and
the user runs the comparison again.

---

## Implementation Plan

Mode D builds on the Modes A–C infrastructure. Implement after Mode A is stable.

### Step 1: `parse_mode_d()` + availability observer (~30 min)
> "Add `parse_mode_d()` to `server_comparator.R`. Add the `mode_d_available`
> reactive and the observer that updates it. Add the conditional fourth radio
> option to the mode selector UI."

Checkpoint: Mode D option appears when both DIA and DDA pipelines have run.

### Step 2: Layer 1 + 2 for Mode D (~1 hr)
> "Add `build_settings_diff_mode_d()`. Wire `parse_mode_d()` into the
> `comparator_data()` reactive when mode = 'dia_vs_dda'. Add the Mode D
> interpretation card to the Layer 2 UpSet output. Confirm Layer 1 settings
> diff renders the extended DDA rows."

Checkpoint: Settings diff and UpSet plot render for a real DIA+DDA session.

### Step 3: Layer 3 extensions + peptide scatter (~1 hr)
> "Add the peptide evidence scatter (3d) as a new nav_panel in the
> comparator results navset. Add the conditionalPanel wrapper so it
> only appears in Mode D. Wire the plotly output."

Checkpoint: Peptide scatter renders for shared proteins with both peptide counts.

### Step 4: Layer 4 extensions + hypothesis engine (~1 hr)
> "Add the Mode D columns (impute_flag, prop_obs) to the discordant table.
> Add `assign_hypothesis_mode_d()` — call this instead of `assign_hypothesis()`
> when mode = 'dia_vs_dda'. Verify all 11 rules fire correctly on test data."

Checkpoint: Discordant table shows Mode D columns; hypothesis text is accurate.

### Step 5: Gemini prompt + polish (~30 min)
> "Add `build_gemini_comparator_prompt_mode_d()`. Wire it into the Gemini
> button handler when mode = 'dia_vs_dda'. Confirm the prompt includes
> correct DDA pipeline settings from the reactive values."

Checkpoint: Gemini generates a Mode D-aware narrative.

---

## Testing Checklist (Mode D)

### Availability
- [ ] Mode D option hidden when only DIA pipeline complete
- [ ] Mode D option hidden when only DDA pipeline complete
- [ ] Mode D option appears when both pipelines complete (both de_fit present)
- [ ] Mode D option disappears if session reset clears DDA results

### Parse + Data Flow
- [ ] `parse_mode_d()` extracts correct topTable for selected contrasts
- [ ] Contrast mismatch warning appears when names differ
- [ ] `protein_id` normalization consistent between DIA and DDA outputs
- [ ] Settings list populates correctly from reactive values (no "unknown" for known fields)

### Layer 1
- [ ] "Acquisition mode" row shows "DIA (diaPASEF)" vs "DDA (ddaPASEF)"
- [ ] FASTA row shows "match" if same path was used for both runs
- [ ] All expected rows present; no missing parameter rows

### Layer 2
- [ ] UpSet plot renders with correct intersection counts
- [ ] Mode D interpretation card appears below UpSet
- [ ] % shared, n_dia_only, n_dda_only are accurate
- [ ] Biological interpretation text is contextually correct

### Layer 3
- [ ] Global scatter renders for shared + quantified proteins
- [ ] Peptide evidence scatter (3d) appears only in Mode D
- [ ] Peptide scatter axes labelled correctly (DIA precursors vs DDA peptides)
- [ ] Identity diagonal line renders correctly

### Layer 4
- [ ] 2×2 concordance table renders with correct counts
- [ ] Discordant table shows `impute_flag` column for DDA-heavy proteins
- [ ] Discordant table shows `prop_obs` column for DIA results
- [ ] Rule 8 (heavy DDA imputation) fires for proteins with impute_frac > 0.5
- [ ] Rule 9 (low DIA PropObs) fires for appropriate proteins
- [ ] Rule 10 (rich DIA, sparse DDA) fires when n_pep_dia >> n_pep_dda
- [ ] Rule 11 (rollup structural difference) fires for same-dir, similar FC, divergent p

### Gemini
- [ ] Mode D prompt includes DIA-NN version and Sage version
- [ ] Prompt includes Perseus width/shift parameters
- [ ] Prompt includes n_heavily_imputed count
- [ ] Response addresses DIA vs DDA acquisition difference specifically

### Regression (Modes A–C unaffected)
- [ ] Mode A–C still parse and render correctly after Mode D additions
- [ ] `assign_hypothesis()` (Modes A–C) unchanged; Mode D uses `assign_hypothesis_mode_d()`
- [ ] Session save/load for Mode A–C unaffected

---

*Addendum version 1.0 — Brett Phinney / UC Davis Proteomics Core — March 2026*
*Extends RUN_COMPARATOR_SPEC_V2.md with Mode D (DIA vs DDA)*
