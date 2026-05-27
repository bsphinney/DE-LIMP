# Three-Mode Acquisition Switcher — Addendum

> **Extends:** `DDA_SAGE_CASANOVO_SPEC.md` (§2 Mode Switcher) and
> `XLMS_INTEGRATION_SPEC.md` (§6 UI)
>
> **Purpose:** Locks in the decision to make XL-MS a **sibling acquisition mode**
> alongside DIA and DDA — not nested under DDA. Both parent specs remain authoritative
> for their own pipelines. This addendum defines only what changes at the nav/switcher
> level and what is shared across all three modes.
>
> **When implementing:** Read this addendum FIRST, then read the relevant parent spec.
> This addendum takes precedence over §2 of DDA_SAGE_CASANOVO_SPEC.md wherever they conflict.

---

## Decision Summary

Three sibling modes. One radio switcher. All share the same Search tab entry point.

```
○ DIA   → DIA-NN + limpa pipeline       (existing)
○ DDA   → Sage + Casanovo pipeline      (DDA spec)
○ XL-MS → MeroX + xiSearch + network   (XLMS spec)
```

XL-MS is DDA acquisition in hardware, but answers a completely different biological
question (protein contacts vs protein abundance). Grouping as sibling modes keeps
each workflow's UI clean and prevents users from seeing irrelevant options.

---

## 1. Mode Switcher UI — Updated Implementation

Replaces the two-choice switcher in `DDA_SAGE_CASANOVO_SPEC.md` §2.
Drop this into `ui.R` — top of sidebar, below the DE-LIMP logo, above the accordion.

```r
# ui.R — acquisition mode switcher (three modes)
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
    "Acquisition Mode"
  ),
  shinyWidgets::radioGroupButtons(
    inputId  = "acquisition_mode",
    label    = NULL,
    choices  = c(
      "<i class='fa fa-infinity'></i> DIA"    = "dia",
      "<i class='fa fa-list-ul'></i> DDA"     = "dda",
      "<i class='fa fa-link'></i> XL-MS"      = "xlms"
    ),
    selected  = "dia",   # DIA always default — zero friction for existing users
    justified = FALSE,
    size      = "sm",
    status    = "outline-primary"
  ),
  uiOutput("mode_context_label")
)
```

### Contextual label (server-side)

```r
# server — mode_context_label reactive
# Replace the two-branch version in DDA spec §2 with this three-branch version
output$mode_context_label <- renderUI({
  mode <- input$acquisition_mode %||% "dia"
  switch(mode,
    "dia" = tags$span(
      style = "font-size: 12px; color: #0d6efd; font-style: italic;",
      icon("circle-check", style = "color: #198754;"),
      " DIA-NN + limpa pipeline"
    ),
    "dda" = tags$span(
      style = "font-size: 12px; color: #6c757d; font-style: italic;",
      icon("circle-info", style = "color: #0d6efd;"),
      " Sage + Casanovo pipeline"
    ),
    "xlms" = tags$span(
      style = "font-size: 12px; color: #6c757d; font-style: italic;",
      icon("diagram-project", style = "color: #6f42c1;"),
      " MeroX + xiSearch + network"
    )
  )
})
```

---

## 2. Search Tab — Three-Way `conditionalPanel` Split

Replace the two-panel split in `DDA_SAGE_CASANOVO_SPEC.md` §2 with this:

```r
# Search tab body — three mutually exclusive panels

# DIA panel (existing — zero changes)
conditionalPanel(
  condition = "input.acquisition_mode === 'dia'",
  # ... all existing DIA-NN search UI unchanged ...
)

# DDA panel (from DDA spec §10.2 — zero changes)
conditionalPanel(
  condition = "input.acquisition_mode === 'dda'",
  # ... Sage + Casanovo search UI ...
)

# XL-MS panel (from XLMS spec §6.1 — zero changes)
conditionalPanel(
  condition = "input.acquisition_mode === 'xlms'",
  # ... MeroX + xiSearch UI from xlms_setup_ui() ...
)
```

---

## 3. Reactive Values — Mode Observer

The existing DDA mode observer in `server_dda.R` watches `input$acquisition_mode`
and syncs to `values$acquisition_mode`. Extend it to handle `"xlms"`:

```r
# In server_dda.R (or a shared server_mode.R if preferred)
observeEvent(input$acquisition_mode, {
  values$acquisition_mode <- input$acquisition_mode

  # XL-MS mode: ensure xlms reactive namespace is initialized
  if (input$acquisition_mode == "xlms" && is.null(xlms$experiment_name)) {
    # xlms$ namespace already initialized in app.R reactiveValues
    # Nothing extra needed — server_xlms.R handles its own state
  }
})
```

---

## 4. `values$acquisition_mode` — Valid States

```r
# Valid values — add "xlms" to the existing "dia" | "dda" enum
values$acquisition_mode  # "dia" | "dda" | "xlms"
```

Downstream panels that check `values$acquisition_mode` (QC, DE, Export) only need
to handle `"dia"` and `"dda"` — XL-MS has its own completely separate Results &
Network panel and does not flow through the limma DE pipeline. No changes needed
to `server_de.R`, `server_qc.R`, or `server_gsea.R` for XL-MS.

---

## 5. Feature Flags in `config.yml`

All three modes are independently gate-able. XL-MS and DDA are off by default on
HuggingFace Spaces (no Hive access). DIA is always on.

```yaml
# config.yml — add xlms flag alongside existing dda flag
features:
  enable_dda:   true    # Hive deployment only
  enable_xlms:  true    # Hive deployment only

# Tool paths — all three pipelines in one place
tools:
  # DIA-NN (existing)
  diann_sif: "/share/proteomics/tools/diann_2.3.0.sif"

  # DDA — Sage + Casanovo
  sage_bin:  "/share/proteomics/tools/sage"

  # XL-MS — MeroX + xiSearch + timsrust
  merox_jar:      "/share/proteomics/tools/xlms/MeroX_2025.jar"
  xi_jar:         "/share/proteomics/tools/xlms/xiSearch.jar"
  xifdr_jar:      "/share/proteomics/tools/xlms/xiFDR.jar"
  timsrust_bin:   "/share/proteomics/tools/timsrust_mgf"  # compiled Rust binary
```

UI gating (add to `ui.R`):

```r
# Hide DDA and XL-MS switcher options on HF Spaces
# (same pattern as existing DIA-NN search tab hiding)
is_hive <- nzchar(Sys.which("sbatch"))

mode_choices <- c("<i class='fa fa-infinity'></i> DIA" = "dia")
if (is_hive && isTRUE(config$features$enable_dda)) {
  mode_choices <- c(mode_choices, "<i class='fa fa-list-ul'></i> DDA" = "dda")
}
if (is_hive && isTRUE(config$features$enable_xlms)) {
  mode_choices <- c(mode_choices, "<i class='fa fa-link'></i> XL-MS" = "xlms")
}

# Pass mode_choices into radioGroupButtons() choices argument above
```

---

## 6. Shared Infrastructure Across All Three Modes

These components are shared — implement once, used by all three:

| Component | Where defined | Used by |
|---|---|---|
| `check_array_job()` | `server_xlms.R` (or extract to `helpers_hpc.R`) | DDA + XL-MS |
| `submit_conversion_array_job()` | `server_xlms.R` | XL-MS (DDA uses Sage direct `.d` read — no conversion needed) |
| SLURM account/partition config | `config.yml` | DDA + XL-MS |
| Claude export ZIP infrastructure | `server_ai.R` | All three |
| Session save/load | `server_session.R` | All three |
| `mode_context_label` reactive | `server_dda.R` or `server_mode.R` | All three |

**Note on conversion:** DDA with Sage does NOT need timsrust MGF conversion — Sage
reads `.d` natively via timsrust. Only XL-MS needs the conversion step (MeroX and
xiSearch require MGF input). The `submit_conversion_array_job()` function is
XL-MS-specific.

---

## 7. Nav Tab Visibility — Mode-Aware Progressive Reveal

The existing progressive reveal logic (tabs appear as data becomes available)
extends naturally:

```r
# Existing DIA reveal: QC appears after data loads, DE after pipeline runs
# DDA reveal: same pattern, driven by values$dda$elist being non-null
# XL-MS reveal: Results & Network tab appears after xlms$consensus is non-null

# In server.R nav visibility observer — add xlms condition:
observe({
  # XL-MS Results tab — only shown when search is complete
  shinyjs::toggleElement(
    "xlms_results_tab",
    condition = !is.null(xlms$consensus)
  )
})
```

---

## 8. Session Save/Load — Mode Persistence

The DDA spec already saves/loads `values$acquisition_mode`. Ensure the load
handler also covers `"xlms"`:

```r
# server_session.R — load handler
# The existing code already does:
#   values$acquisition_mode <- session_data$acquisition_mode
# This already handles "xlms" because it's just a string.
# No change needed as long as "xlms" is a valid value.

# However, add XL-MS state to the save block:
session_data$xlms <- list(
  consensus       = xlms$consensus,
  protein_pairs   = xlms$protein_pairs,
  qc_stats        = xlms$qc_stats,
  experiment_name = xlms$experiment_name,
  crosslinker     = xlms$crosslinker
  # Note: job IDs not saved — HPC jobs don't survive session reload
)

# And the load block:
if (!is.null(session_data$xlms) && !is.null(session_data$xlms$consensus)) {
  xlms$consensus       <- session_data$xlms$consensus
  xlms$protein_pairs   <- session_data$xlms$protein_pairs
  xlms$qc_stats        <- session_data$xlms$qc_stats
  xlms$experiment_name <- session_data$xlms$experiment_name
  xlms$crosslinker     <- session_data$xlms$crosslinker
}
```

---

## 9. Optional: Sage Linear Search Alongside XL-MS

The XL-MS mode can optionally run Sage on the same `.d` files as a linear peptide
search to confirm crosslinker efficiency and protein coverage. This is a QC feature,
not required for MVP.

When implemented, add a checkbox to the XL-MS Setup panel:

```r
checkboxInput("xlms_run_sage_qc",
  "Also run Sage linear search for protein coverage QC (optional)",
  value = FALSE)
```

If checked, submit a Sage SLURM job alongside the MeroX/xiSearch arrays, using the
same FASTA. The Sage `lfq.tsv` output feeds a "Crosslinker Efficiency" QC card showing:
- % of FASTA proteins detected as linear peptides
- Crosslinker modification rate (dead-end + crosslinked vs unmodified)
- Estimated crosslinking efficiency per protein

This reuses `submit_sage_job()` from `server_dda.R` directly — zero new code needed
beyond wiring the checkbox.

---

## 10. Implementation Order

When Claude Code implements this, do it in this order:

1. **First: extend the mode switcher** in `ui.R` from 2 choices to 3 (30 min)
2. **First: extend the mode observer** in `server_dda.R` to handle `"xlms"` (10 min)
3. **First: add feature flags** to `config.yml` and mode_choices gating (20 min)
4. **Then: implement DDA spec** as written — no changes needed to DDA internals
5. **Then: implement XLMS spec** as written — no changes needed to XLMS internals
6. **Last: optional** Sage linear search QC checkbox in XL-MS mode

The first three steps are the only actual cross-cutting changes. Everything else
is the two parent specs implemented independently.

---

## 11. What Does NOT Change

To be explicit — these things are unchanged by this addendum:

- All existing DIA-NN pipeline code (`server_data.R`, `server_de.R`, `server_qc.R`, etc.)
- The DDA spec's reactive namespace (`dda$`), job submission, parsing, DE pipeline
- The XLMS spec's reactive namespace (`xlms$`), job submission, parsing, network viz
- HF Spaces deployment — DDA and XL-MS tabs simply don't appear (no sbatch available)
- The Claude export system — each mode adds its own template; shared ZIP infrastructure unchanged
- The `values$fit` / `values$y_protein` / `values$metadata` DIA pipeline — XL-MS never touches these
