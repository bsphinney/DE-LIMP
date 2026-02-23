# Core Facility Mode — Phase 2 Specification

> **Version**: Targets DE-LIMP v3.1+
> **Prerequisite**: Phase 1 is implemented (branch `claude/competent-shtern` merged from `claude/heuristic-wu`)
> **Date**: February 2026

---

## What's Already Implemented (Phase 1)

Before tackling Phase 2, here's what Phase 1 delivered and where the code lives:

| Feature | Files | Status |
|---------|-------|--------|
| Core facility detection (`DELIMP_CORE_DIR` env var) | `app.R` lines 150-195 | Done |
| SQLite database (4 tables: searches, qc_runs, reports, templates) | `R/helpers_facility.R` `cf_init_db()` | Done |
| Staff SSH selector + auto-connect | `R/server_facility.R` lines 20-88 | Done |
| Project auto-populate (staff+lab filtering) | `R/server_facility.R` lines 92-128 | Done |
| Report generation (Quarto HTML + state .rds) | `R/server_facility.R` `generate_report_impl()` lines 136-329 | Done |
| HF state upload | `R/helpers_facility.R` `upload_state_to_hf()` | Done |
| Search DB tab (6-filter job history table) | `R/ui.R` lines 1288-1341, `R/server_facility.R` lines 396-560 | Done |
| Instrument QC dashboard (protein/precursor/TIC trends) | `R/ui.R` lines 1343-1364, `R/server_facility.R` lines 564-710 | Done |
| Template library (save/load search presets, group templates) | `R/server_facility.R` lines 714-864 | Done |
| Job recording on search submit | `R/server_search.R` line 1062 → `cf_record_search()` | Done |
| Job status updates | `R/helpers_facility.R` `cf_update_search_status()` | Done |
| QC bracket in report template | `report_template.qmd` lines 65-127 | Done |
| Synthetic test data generator | `seed_test_db.R` | Done |

### Current SQLite Schema

```sql
-- searches: DIA-NN job tracking
searches (id, analysis_name, submitted_by, lab, submitted_at, completed_at,
          status, slurm_job_id, container_id, backend, instrument,
          lc_method, project, organism, fasta_file, n_raw_files,
          search_mode, n_proteins, n_peptides, n_precursors,
          output_dir, report_id, notes)

-- qc_runs: HeLa QC metrics per instrument
qc_runs (id, run_name, instrument, run_date, ingested_at, file_path,
         n_proteins, n_peptides, n_precursors, median_ms1_tic,
         median_cv, protein_intensities, search_id, report_path)

-- reports: Generated report metadata
reports (id, title, created_at, created_by, lab, instrument,
         n_proteins, n_contrasts, contrast_names, search_id,
         qc_before_id, qc_after_id, html_path, state_path, hf_state_url)

-- templates: Saved search presets / group templates
templates (id, name, type, created_by, created_at, config_json, notes)
```

---

## Phase 2 Features — Table of Contents

1. [QC Run Ingestion](#1-qc-run-ingestion)
2. [Report Template Polish](#2-report-template-polish)
3. [Report Comparison](#3-report-comparison)
4. [HF State Upload/Download](#4-hf-state-uploaddownload)
5. [Template Application on Search Submit](#5-template-application-on-search-submit)
6. [Audit Log](#6-audit-log)
7. [Multi-Instrument QC Alerts](#7-multi-instrument-qc-alerts)
8. [End-to-End Testing](#8-end-to-end-testing)

### Implementation Priority

| Priority | Feature | Effort | Dependencies |
|----------|---------|--------|-------------|
| **P1** | QC Run Ingestion | 4 hrs | None (uses existing `qc_runs` table) |
| **P1** | End-to-End Testing | 3 hrs | None (validates Phase 1) |
| **P2** | Report Template Polish | 4 hrs | None |
| **P2** | Multi-Instrument QC Alerts | 3 hrs | QC Run Ingestion (needs real data) |
| **P3** | HF State Upload/Download | 4 hrs | None (extends existing `upload_state_to_hf`) |
| **P3** | Template Application on Search Submit | 2 hrs | None (templates already saved in SQLite) |
| **P3** | Audit Log | 3 hrs | None (new SQLite table) |
| **P4** | Report Comparison | 5 hrs | Report Template Polish |

**Suggested order**: QC Ingestion → E2E Testing → Report Polish → QC Alerts → HF State → Template Auto-Apply → Audit Log → Report Comparison

---

## 1. QC Run Ingestion

### Concept

When a user loads a `report.parquet` file, auto-detect whether it's a HeLa QC digest run. If so, extract metrics (protein count, precursor count, MS1 TIC) and record them in the `qc_runs` SQLite table. This populates the Instrument QC dashboard with real data without requiring an external cron job.

### Detection Strategy

Two detection methods (user chooses or auto-detect):

1. **Filename pattern**: `(?i)hela|qc[_-]?digest|qc[_-]?standard` in the uploaded filename
2. **Checkbox**: Add a "This is a QC run" checkbox in the Data Overview sidebar that appears in core facility mode

Auto-detection shows a confirmation banner: "This looks like a HeLa QC run. Record QC metrics?" with Yes/No buttons.

### Files to Modify

| File | Changes |
|------|---------|
| `R/ui.R` | Add QC ingestion checkbox + confirmation banner UI (core facility sidebar, ~15 lines) |
| `R/server_data.R` | Add auto-detection logic after file upload (~30 lines) |
| `R/server_facility.R` | Add `ingest_qc_run()` observer (~60 lines) |
| `R/helpers_facility.R` | Add `cf_ingest_qc_metrics()` helper function (~50 lines) |

### UI Changes — `R/ui.R`

In the sidebar, after the file upload section, add (wrapped in `is_core_facility`):

```r
if (is_core_facility) tagList(
  checkboxInput("is_qc_run", "Record as QC run", value = FALSE),
  conditionalPanel("input.is_qc_run",
    selectInput("qc_run_instrument", "Instrument",
      choices = c("(select)" = "", cf_instrument_names(cf_config))),
    tags$small(class = "text-muted",
      "QC metrics will be recorded after pipeline completes.")
  ),
  uiOutput("qc_detection_banner")
)
```

### Server Logic — `R/server_data.R`

After successful file upload (in the existing upload observer), add auto-detection:

```r
# After values$uploaded_report_path is set and values$original_report_name is known:
if (is_core_facility) {
  fname <- tolower(values$original_report_name %||% "")
  if (grepl("hela|qc[_-]?digest|qc[_-]?standard", fname)) {
    output$qc_detection_banner <- renderUI({
      tags$div(class = "alert alert-info py-2 px-3 mt-2",
        style = "font-size: 0.85em;",
        icon("flask"), " This looks like a QC run.",
        actionButton("confirm_qc_ingest", "Record QC Metrics",
          class = "btn-sm btn-outline-primary ms-2"),
        actionButton("dismiss_qc_banner", "Dismiss",
          class = "btn-sm btn-outline-secondary ms-1")
      )
    })
  }
}
```

### Server Logic — `R/server_facility.R`

Add QC ingestion observer:

```r
# Trigger after pipeline completes (values$fit becomes non-NULL) if QC checkbox is checked
observeEvent(input$confirm_qc_ingest, {
  req(values$y_protein)
  req(input$qc_run_instrument, nzchar(input$qc_run_instrument))

  cf_ingest_qc_metrics(
    db_path     = cf_config$db_path,
    run_name    = values$original_report_name %||% "QC_run",
    instrument  = input$qc_run_instrument,
    y_protein   = values$y_protein,
    qc_stats    = values$qc_stats,
    file_path   = values$uploaded_report_path
  )

  output$qc_detection_banner <- renderUI({
    tags$div(class = "alert alert-success py-2 px-3 mt-2",
      style = "font-size: 0.85em;",
      icon("check-circle"), " QC metrics recorded for ",
      tags$strong(input$qc_run_instrument))
  })

  showNotification("QC run recorded!", type = "message")
})

# Also support the checkbox-based flow (record after pipeline runs)
observe({
  req(input$is_qc_run)
  req(values$fit)  # pipeline has completed
  req(input$qc_run_instrument, nzchar(input$qc_run_instrument))

  cf_ingest_qc_metrics(
    db_path     = cf_config$db_path,
    run_name    = values$original_report_name %||% "QC_run",
    instrument  = input$qc_run_instrument,
    y_protein   = values$y_protein,
    qc_stats    = values$qc_stats,
    file_path   = values$uploaded_report_path
  )

  updateCheckboxInput(session, "is_qc_run", value = FALSE)
  showNotification("QC metrics recorded!", type = "message")
}) |> bindEvent(values$fit, ignoreInit = TRUE)
```

### Helper — `R/helpers_facility.R`

```r
#' Extract and record QC metrics from a completed analysis
#' @param db_path Path to SQLite database
#' @param run_name Filename of the QC run
#' @param instrument Instrument name
#' @param y_protein EList or matrix of protein-level quantification
#' @param qc_stats Optional QC stats from DIA-NN
#' @param file_path Optional path to the source file
cf_ingest_qc_metrics <- function(db_path, run_name, instrument,
                                  y_protein, qc_stats = NULL,
                                  file_path = NULL) {
  # Extract expression matrix
  y <- if (is.list(y_protein) && "E" %in% names(y_protein)) {
    y_protein$E
  } else {
    as.matrix(y_protein)
  }

  n_proteins <- nrow(y)

  # Precursor count from qc_stats if available
  n_precursors <- if (!is.null(qc_stats) && "n_precursors" %in% names(qc_stats)) {
    qc_stats$n_precursors
  } else NA

  n_peptides <- if (!is.null(qc_stats) && "n_peptides" %in% names(qc_stats)) {
    qc_stats$n_peptides
  } else NA

  # Median MS1 TIC from qc_stats
  median_tic <- if (!is.null(qc_stats) && "median_ms1_tic" %in% names(qc_stats)) {
    qc_stats$median_ms1_tic
  } else NA

  # Median CV across proteins (use column medians)
  med_cv <- tryCatch({
    protein_means <- rowMeans(y, na.rm = TRUE)
    protein_sds   <- apply(y, 1, sd, na.rm = TRUE)
    cvs <- protein_sds / protein_means * 100
    median(cvs, na.rm = TRUE)
  }, error = function(e) NA)

  # Top protein intensities as JSON (for rolling CV computation)
  protein_intensities <- tryCatch({
    protein_medians <- apply(y, 1, median, na.rm = TRUE)
    top_n <- min(500, length(protein_medians))
    top <- sort(protein_medians, decreasing = TRUE)[seq_len(top_n)]
    jsonlite::toJSON(as.list(top), auto_unbox = TRUE)
  }, error = function(e) "{}")

  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(db))

  DBI::dbExecute(db, "
    INSERT INTO qc_runs (run_name, instrument, run_date, file_path,
                         n_proteins, n_peptides, n_precursors,
                         median_ms1_tic, median_cv, protein_intensities)
    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10)",
    params = list(
      run_name, instrument, as.character(Sys.time()),
      file_path %||% NA_character_,
      n_proteins, n_peptides %||% NA, n_precursors %||% NA,
      median_tic %||% NA, round(med_cv, 1) %||% NA,
      as.character(protein_intensities)
    )
  )
}
```

### Gotchas

- `y_protein` is an EList, not a matrix — always extract `$E`
- `qc_stats` structure depends on what DIA-NN reports; check field names defensively
- Duplicate run prevention: check `run_name + instrument` before INSERT to avoid re-ingesting the same run
- The ingestion should happen AFTER the pipeline completes (not on upload), because we need protein counts from the quantified data

---

## 2. Report Template Polish

### Concept

Enhance `report_template.qmd` with additional analysis sections and configurable branding.

### Changes

#### A. Add GSEA Section

After the Phosphoproteomics section in `report_template.qmd`:

```r
## Gene Set Enrichment Analysis

```{r gsea-section, results='asis', fig.height=5}
has_gsea <- !is.null(state$gsea_results_cache) &&
            length(state$gsea_results_cache) > 0

if (!has_gsea) {
  cat("*No GSEA results included in this analysis.*\n")
}

if (has_gsea) {
  library(enrichplot)
  library(clusterProfiler)

  for (ont_name in names(state$gsea_results_cache)) {
    res <- state$gsea_results_cache[[ont_name]]
    if (is.null(res) || nrow(res@result) == 0) next

    cat(sprintf("\n### %s\n\n", ont_name))

    n_sig <- sum(res@result$p.adjust < 0.05, na.rm = TRUE)
    cat(sprintf("**%d significant terms** (FDR < 0.05)\n\n", n_sig))

    # Dot plot (top 15)
    if (n_sig > 0) {
      tryCatch({
        print(dotplot(res, showCategory = 15, font.size = 9) +
          ggplot2::ggtitle(paste("GSEA:", ont_name)))
      }, error = function(e) {
        cat(sprintf("*Could not render dot plot: %s*\n", e$message))
      })
    }
    cat("\n\n")
  }
}
```

**Dependency**: `enrichplot` and `clusterProfiler` must be available at render time. If not installed, skip the section gracefully with `requireNamespace()` guards.

#### B. Add MOFA Variance Explained Section

```r
## Multi-View Integration (MOFA2)

```{r mofa-section, results='asis', fig.height=4}
has_mofa <- !is.null(state$mofa_variance_explained)

if (!has_mofa) {
  cat("*No MOFA2 integration results included in this analysis.*\n")
}

if (has_mofa) {
  ve <- state$mofa_variance_explained

  # Heatmap of variance explained per view per factor
  if (is.matrix(ve) || is.data.frame(ve)) {
    ve_mat <- as.matrix(ve)

    # Simple heatmap using base R
    par(mar = c(5, 8, 3, 2))
    image(t(ve_mat[nrow(ve_mat):1, , drop = FALSE]),
          axes = FALSE,
          col = colorRampPalette(c("white", "#2166AC"))(50),
          main = "Variance Explained per View per Factor")
    axis(1, at = seq(0, 1, length.out = ncol(ve_mat)),
         labels = colnames(ve_mat), las = 2, cex.axis = 0.8)
    axis(2, at = seq(0, 1, length.out = nrow(ve_mat)),
         labels = rev(rownames(ve_mat)), las = 1, cex.axis = 0.8)

    # Add text values
    for (i in seq_len(nrow(ve_mat))) {
      for (j in seq_len(ncol(ve_mat))) {
        text((j - 1) / (ncol(ve_mat) - 1),
             1 - (i - 1) / (nrow(ve_mat) - 1),
             sprintf("%.1f%%", ve_mat[i, j]),
             cex = 0.7)
      }
    }
  }

  cat(sprintf("\n\n**%d views** integrated across **%d factors**.\n",
              nrow(ve_mat), ncol(ve_mat)))
}
```

#### C. Configurable Logo/Header

Add an optional `logo_path` parameter to the Quarto template:

```yaml
params:
  analysis_id: "unknown"
  state_dir: "."
  db_path: "delimp.db"
  logo_path: ""           # NEW: path to facility logo (PNG/SVG)
  facility_name: ""       # NEW: e.g., "UC Davis Proteomics Core"
```

In the report header section:

```r
```{r header, results='asis'}
if (nzchar(params$facility_name)) {
  cat(sprintf('<div style="text-align: center; margin-bottom: 20px;">\n'))
  if (nzchar(params$logo_path) && file.exists(params$logo_path)) {
    cat(sprintf('<img src="%s" alt="Logo" style="max-height: 80px; margin-bottom: 10px;">\n',
                params$logo_path))
  }
  cat(sprintf('<h2>%s</h2>\n', params$facility_name))
  cat('</div>\n')
}
```

#### D. Configuration

Add to `qc_config.yml`:

```yaml
report:
  facility_name: "UC Davis Proteomics Core"
  logo_path: "/srv/delimp/logo.png"
```

Update `generate_report_impl()` in `server_facility.R` to pass these as `execute_params`:

```r
quarto::quarto_render(
  input = cf_config$template_qmd,
  output_file = paste0(analysis_id, ".html"),
  execute_params = list(
    analysis_id = analysis_id,
    state_dir   = cf_config$state_dir,
    db_path     = cf_config$db_path,
    logo_path   = cf_config$qc$report$logo_path %||% "",
    facility_name = cf_config$qc$report$facility_name %||% ""
  )
)
```

### Files to Modify

| File | Changes |
|------|---------|
| `report_template.qmd` | Add GSEA section (~40 lines), MOFA section (~35 lines), logo/header (~15 lines), new params |
| `R/server_facility.R` | Pass new params to `quarto_render()` (~5 lines) |
| `R/helpers_facility.R` | Read report config from `qc_config.yml` (already loaded in `cf_config$qc`) |

---

## 3. Report Comparison

### Concept

Side-by-side comparison of two reports on the same instrument: QC bracket metrics, DE summary, and protein overlap. Useful for longitudinal tracking (e.g., "Did the January batch perform as well as December?").

### UI — New Sub-Tab in Search DB

Add a "Compare Reports" button next to "Generate Report" in the Search DB tab:

```r
actionButton("job_compare_reports", "Compare Two Reports",
  icon = icon("columns"),
  class = "btn-outline-info btn-sm")
```

Clicking it opens a modal:

```r
observeEvent(input$job_compare_reports, {
  db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
  on.exit(DBI::dbDisconnect(db))

  reports <- DBI::dbGetQuery(db,
    "SELECT id, title, created_at, instrument, n_proteins, lab
     FROM reports ORDER BY created_at DESC LIMIT 50")

  choices <- setNames(reports$id,
    paste0(reports$title, " (", reports$instrument, ", ",
           format(as.POSIXct(reports$created_at), "%b %d"), ")"))

  showModal(modalDialog(
    title = "Compare Two Reports",
    selectInput("compare_report_a", "Report A (baseline)",
      choices = c("(select)" = "", choices)),
    selectInput("compare_report_b", "Report B (comparison)",
      choices = c("(select)" = "", choices)),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("run_comparison", "Compare", class = "btn-primary")
    )
  ))
})
```

### Comparison Output

Renders a comparison card (inline or as a modal with tabs):

**Tab 1: QC Bracket Comparison**

| Metric | Report A | Report B | Delta |
|--------|----------|----------|-------|
| QC Proteins (before) | 6,234 | 6,189 | -45 |
| QC Proteins (after) | 6,301 | 6,250 | -51 |
| QC Median CV (before) | 12.3% | 13.1% | +0.8% |

**Tab 2: DE Summary Comparison**

| Contrast | Up (A) | Up (B) | Down (A) | Down (B) |
|----------|--------|--------|----------|----------|
| Treatment_vs_Control | 245 | 312 | 189 | 201 |

**Tab 3: Protein Overlap**

- Venn-style counts: unique to A, shared, unique to B
- Jaccard similarity index

### Server Logic

```r
observeEvent(input$run_comparison, {
  removeModal()
  req(input$compare_report_a, input$compare_report_b)
  req(input$compare_report_a != input$compare_report_b)

  # Load both state files
  db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
  on.exit(DBI::dbDisconnect(db))

  report_a <- DBI::dbGetQuery(db, "SELECT * FROM reports WHERE id = ?1",
    params = list(input$compare_report_a))
  report_b <- DBI::dbGetQuery(db, "SELECT * FROM reports WHERE id = ?1",
    params = list(input$compare_report_b))

  state_a <- readRDS(report_a$state_path)
  state_b <- readRDS(report_b$state_path)

  # Build comparison data...
  # (QC bracket from reports table qc_before_id/qc_after_id)
  # (DE counts from state$fit)
  # (Protein overlap from rownames of state$y_protein$E or state$y_protein)

  # Render comparison UI
  output$comparison_output <- renderUI({ ... })
})
```

### Files to Modify

| File | Changes |
|------|---------|
| `R/ui.R` | Add "Compare" button (~5 lines), add `uiOutput("comparison_output")` |
| `R/server_facility.R` | Add comparison modal + observer (~100 lines) |

### Gotchas

- Both `.rds` state files must still exist on disk
- Handle cases where instruments differ between reports (warn user)
- Memory: loading two full state objects simultaneously — consider extracting only needed fields

---

## 4. HF State Upload/Download

### Concept

Complete the HF Spaces live link feature. Phase 1 implemented `upload_state_to_hf()` but the **download side on HF Spaces** is not yet wired up. When a user visits `?analysis=UUID`, the HF Space should download the `.rds` and restore the full session.

### Prerequisites

1. Create HF dataset repo: `brettsp/delimp-shared-analyses` (private or gated)
2. Set `HF_TOKEN` secret in HF Space settings

### Server Logic — `R/server_session.R`

Add URL query parameter handler (runs only on HF Spaces):

```r
# Load shared analysis from URL query parameter
# Only on HF Spaces — lab workstation loads from local disk
if (is_hf_space) {
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (is.null(query$analysis)) return()

    analysis_id <- query$analysis

    # Validate UUID format (prevent path traversal)
    if (!grepl("^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$",
               analysis_id)) {
      showNotification("Invalid analysis link.", type = "error")
      return()
    }

    rds_url <- paste0(
      "https://huggingface.co/datasets/brettsp/delimp-shared-analyses/",
      "resolve/main/", analysis_id, ".rds")

    withProgress(message = "Loading shared analysis...", {
      tryCatch({
        temp_rds <- tempfile(fileext = ".rds")
        download.file(rds_url, temp_rds, mode = "wb", quiet = TRUE)
        state <- readRDS(temp_rds)

        # Restore reactive values (same pattern as session load)
        values$metadata   <- state$metadata
        values$fit        <- state$fit
        values$y_protein  <- state$y_protein
        values$dpc_fit    <- state$dpc_fit
        values$design     <- state$design
        values$qc_stats   <- state$qc_stats
        values$gsea_results_cache <- state$gsea_results_cache
        values$gsea_last_contrast <- state$gsea_last_contrast
        values$gsea_last_org_db   <- state$gsea_last_org_db
        values$diann_norm_detected <- state$diann_norm_detected %||% "unknown"

        # Phospho
        values$phospho_detected <- state$phospho_detected
        values$phospho_site_matrix <- state$phospho_site_matrix
        values$phospho_site_info <- state$phospho_site_info
        values$phospho_fit <- state$phospho_fit
        values$phospho_site_matrix_filtered <- state$phospho_site_matrix_filtered
        values$phospho_input_mode <- state$phospho_input_mode

        # MOFA
        if (!is.null(state$mofa_object)) {
          values$mofa_object <- state$mofa_object
          values$mofa_factors <- state$mofa_factors
          values$mofa_variance_explained <- state$mofa_variance_explained
        }

        # UI state
        values$color_plot_by_de <- state$color_plot_by_de %||% FALSE

        # Update contrast selectors
        if (!is.null(values$fit)) {
          cn <- colnames(values$fit$contrasts)
          sel <- state$contrast %||% cn[1]
          for (sel_id in c("contrast_selector", "contrast_selector_signal",
                           "contrast_selector_grid", "contrast_selector_pvalue")) {
            updateSelectInput(session, sel_id, choices = cn, selected = sel)
          }
        }

        values$status <- "Loaded from shared link"
        showNotification(
          sprintf("Loaded: %s", state$metadata_extra$title %||% "Shared Analysis"),
          type = "message", duration = 8)
        bslib::nav_select("main_tabs", "DE Dashboard")

      }, error = function(e) {
        showNotification(
          paste("Could not load analysis:", e$message),
          type = "error", duration = 10)
      })
    })
  }) |> bindEvent(session$clientData$url_search, once = TRUE)
}
```

### Files to Modify

| File | Changes |
|------|---------|
| `R/server_session.R` | Add URL query parameter handler (~70 lines) |
| `app.R` | Pass `is_hf_space` to `server_session()` if not already available |

### Gotchas

- UUID validation is critical — prevents path traversal attacks
- HF dataset repo must be public or the HF Space must have a token with read access
- `.rds` files can be large (50-200 MB for MOFA objects) — show progress
- Consider a max file size check before downloading
- The `bindEvent(..., once = TRUE)` ensures it only fires once per session

---

## 5. Template Application on Search Submit

### Concept

When submitting a new DIA-NN search, optionally auto-apply a saved template to pre-fill search settings. This saves time for repetitive workflows (e.g., "every AP-MS uses the same DIA-NN settings").

### UI Changes — `R/ui.R`

In the New Search tab's settings section, add a template selector (core facility mode only):

```r
if (is_core_facility) tagList(
  selectInput("search_template", "Apply Template",
    choices = c("(none)" = ""), selected = ""),
  tags$small(class = "text-muted",
    "Template pre-fills search settings below.")
)
```

### Server Logic — `R/server_facility.R`

```r
# Populate search template dropdown (only search_preset templates)
observe({
  db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
  on.exit(DBI::dbDisconnect(db))

  templates <- DBI::dbGetQuery(db,
    "SELECT id, name FROM templates WHERE type = 'search_preset' ORDER BY name")
  choices <- c("(none)" = "")
  if (nrow(templates) > 0) {
    choices <- c(choices, setNames(templates$id, templates$name))
  }
  updateSelectInput(session, "search_template", choices = choices)
})

# Auto-apply when template is selected
observeEvent(input$search_template, {
  req(input$search_template, nzchar(input$search_template))

  db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
  on.exit(DBI::dbDisconnect(db))

  tpl <- DBI::dbGetQuery(db,
    "SELECT config_json FROM templates WHERE id = ?1",
    params = list(as.integer(input$search_template)))
  req(nrow(tpl) > 0)

  config <- jsonlite::fromJSON(tpl$config_json[1])

  # Apply each field if present (reuse existing load_template logic)
  field_map <- list(
    cpus = "diann_cpus", mem_gb = "diann_mem_gb",
    time_hours = "diann_time_hours", partition = "diann_partition",
    account = "diann_account", mass_acc_mode = "mass_acc_mode",
    mass_acc = "diann_mass_acc", ms1_acc = "diann_mass_acc_ms1",
    search_mode = "search_mode", enzyme = "diann_enzyme",
    missed_cleavages = "diann_missed_cleavages",
    max_var_mods = "diann_max_var_mods",
    normalization = "diann_normalization"
  )

  for (cfg_key in names(field_map)) {
    if (!is.null(config[[cfg_key]])) {
      input_id <- field_map[[cfg_key]]
      val <- config[[cfg_key]]
      # Determine update function based on input type
      if (cfg_key %in% c("search_mode", "normalization")) {
        updateRadioButtons(session, input_id, selected = val)
      } else if (cfg_key %in% c("mass_acc_mode", "enzyme")) {
        updateSelectInput(session, input_id, selected = val)
      } else if (cfg_key %in% c("partition", "account")) {
        updateTextInput(session, input_id, value = val)
      } else {
        updateNumericInput(session, input_id, value = val)
      }
    }
  }

  showNotification(sprintf("Applied template: %s", tpl$name[1] %||% ""),
                   type = "message", duration = 3)
})
```

### Files to Modify

| File | Changes |
|------|---------|
| `R/ui.R` | Add `selectInput("search_template", ...)` in search settings (~8 lines) |
| `R/server_facility.R` | Add template populate + apply observers (~40 lines) |

### Note

This largely reuses the existing `load_template` observer logic from Phase 1. The difference is placement (in the search tab rather than sidebar) and filtering to `search_preset` type only.

---

## 6. Audit Log

### Concept

Track who did what and when. Every significant action in core facility mode gets a timestamped audit entry: report generation, job submission, template save, QC ingestion, session load.

### Schema — New SQLite Table

Add to `cf_init_db()` in `R/helpers_facility.R`:

```sql
CREATE TABLE IF NOT EXISTS audit_log (
  id          INTEGER PRIMARY KEY AUTOINCREMENT,
  timestamp   DATETIME DEFAULT CURRENT_TIMESTAMP,
  user_name   TEXT NOT NULL,
  action      TEXT NOT NULL,
  entity_type TEXT,
  entity_id   TEXT,
  details     TEXT,
  ip_address  TEXT
)

CREATE INDEX IF NOT EXISTS idx_audit_timestamp ON audit_log(timestamp)
CREATE INDEX IF NOT EXISTS idx_audit_user ON audit_log(user_name)
CREATE INDEX IF NOT EXISTS idx_audit_action ON audit_log(action)
```

### Helper — `R/helpers_facility.R`

```r
#' Record an audit log entry
#' @param db_path Path to SQLite database
#' @param user_name Staff name
#' @param action Action type (e.g., "generate_report", "submit_search")
#' @param entity_type Type of entity (e.g., "report", "search", "template")
#' @param entity_id ID of the entity
#' @param details Free-text description
#' @param session Optional Shiny session for IP extraction
cf_audit_log <- function(db_path, user_name, action,
                          entity_type = NULL, entity_id = NULL,
                          details = NULL, session = NULL) {
  ip <- tryCatch({
    if (!is.null(session)) session$request$REMOTE_ADDR else NA_character_
  }, error = function(e) NA_character_)

  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(db))

  DBI::dbExecute(db, "
    INSERT INTO audit_log (user_name, action, entity_type, entity_id,
                           details, ip_address)
    VALUES (?1, ?2, ?3, ?4, ?5, ?6)",
    params = list(
      user_name %||% "unknown",
      action,
      entity_type %||% NA_character_,
      entity_id %||% NA_character_,
      details %||% NA_character_,
      ip %||% NA_character_
    )
  )
}
```

### Integration Points

Add `cf_audit_log()` calls at these locations:

| File | Location | Action |
|------|----------|--------|
| `R/server_facility.R` | `generate_report_impl()` after successful render | `"generate_report"` |
| `R/server_search.R` | After `cf_record_search()` | `"submit_search"` |
| `R/server_facility.R` | `confirm_save_template` observer | `"save_template"` |
| `R/server_facility.R` | `load_template` observer | `"apply_template"` |
| `R/server_facility.R` | `cf_ingest_qc_metrics()` (new) | `"ingest_qc_run"` |
| `R/server_facility.R` | `job_load_results` observer | `"load_results"` |
| `R/server_session.R` | Session save (in core facility mode) | `"save_session"` |

### Viewing Audit Log

Add a collapsible section at the bottom of the Search DB tab:

```r
tags$details(
  tags$summary(style = "cursor: pointer; color: #6c757d; font-size: 0.9em;",
    icon("history"), " Audit Log"),
  DTOutput("audit_log_table")
)
```

Server:

```r
output$audit_log_table <- DT::renderDT({
  invalidateLater(60000, session)

  db <- DBI::dbConnect(RSQLite::SQLite(), cf_config$db_path)
  on.exit(DBI::dbDisconnect(db))

  log <- DBI::dbGetQuery(db,
    "SELECT timestamp, user_name as 'User', action as 'Action',
            entity_type as 'Type', details as 'Details'
     FROM audit_log ORDER BY timestamp DESC LIMIT 200")

  DT::datatable(log,
    options = list(pageLength = 10, dom = "tip", scrollX = TRUE),
    rownames = FALSE, class = "compact stripe hover")
})
```

### Files to Modify

| File | Changes |
|------|---------|
| `R/helpers_facility.R` | Add `audit_log` table to `cf_init_db()` (~10 lines), add `cf_audit_log()` helper (~25 lines) |
| `R/server_facility.R` | Add `cf_audit_log()` calls at 5 locations (~10 lines total), add audit log table render (~20 lines) |
| `R/server_search.R` | Add `cf_audit_log()` call on submit (~3 lines) |
| `R/ui.R` | Add collapsible audit log section in Search DB tab (~8 lines) |

---

## 7. Multi-Instrument QC Alerts

### Concept

Automatically flag instruments where the most recent QC run's protein count drops below the rolling mean - 2*SD. Show alerts in the Instrument QC dashboard and optionally in the sidebar.

### Helper — `R/helpers_facility.R`

```r
#' Check all instruments for QC alert conditions
#' @param db_path Path to SQLite database
#' @return Data frame with columns: instrument, latest_proteins, mean_prot, sd_prot, alert
cf_check_qc_alerts <- function(db_path) {
  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(db))

  instruments <- DBI::dbGetQuery(db,
    "SELECT DISTINCT instrument FROM qc_runs WHERE n_proteins IS NOT NULL")$instrument

  alerts <- lapply(instruments, function(inst) {
    # Latest run
    latest <- DBI::dbGetQuery(db, "
      SELECT n_proteins, run_date, run_name FROM qc_runs
      WHERE instrument = ?1 AND n_proteins IS NOT NULL
      ORDER BY run_date DESC LIMIT 1",
      params = list(inst))

    if (nrow(latest) == 0) return(NULL)

    # Rolling stats (90 days, excluding latest)
    stats <- DBI::dbGetQuery(db, "
      SELECT AVG(n_proteins) as mean_prot,
             COUNT(*) as n_runs
      FROM qc_runs
      WHERE instrument = ?1
        AND n_proteins IS NOT NULL
        AND run_date >= datetime('now', '-90 days')
        AND id != (SELECT id FROM qc_runs
                   WHERE instrument = ?1
                   ORDER BY run_date DESC LIMIT 1)",
      params = list(inst))

    # Need SD separately (SQLite has no STDEV)
    vals <- DBI::dbGetQuery(db, "
      SELECT n_proteins FROM qc_runs
      WHERE instrument = ?1
        AND n_proteins IS NOT NULL
        AND run_date >= datetime('now', '-90 days')",
      params = list(inst))$n_proteins

    sd_prot <- if (length(vals) >= 3) sd(vals) else NA

    threshold <- if (!is.na(stats$mean_prot) && !is.na(sd_prot)) {
      stats$mean_prot - 2 * sd_prot
    } else NA

    data.frame(
      instrument = inst,
      latest_proteins = latest$n_proteins,
      latest_run = latest$run_name,
      latest_date = latest$run_date,
      mean_prot = round(stats$mean_prot, 0),
      sd_prot = round(sd_prot, 0),
      threshold = round(threshold, 0),
      n_runs = stats$n_runs,
      alert = !is.na(threshold) && latest$n_proteins < threshold,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), alerts))
}
```

### UI — Alert Banner in Instrument QC Tab

Add at the top of the Instrument QC dashboard, before the trend plots:

```r
uiOutput("qc_alert_banner")
```

### Server Logic — `R/server_facility.R`

```r
# Check QC alerts periodically
output$qc_alert_banner <- renderUI({
  invalidateLater(300000, session)  # check every 5 minutes

  alerts <- tryCatch(
    cf_check_qc_alerts(cf_config$db_path),
    error = function(e) NULL
  )

  if (is.null(alerts) || nrow(alerts) == 0) return(NULL)

  flagged <- alerts[alerts$alert == TRUE, , drop = FALSE]
  if (nrow(flagged) == 0) return(NULL)

  tags$div(class = "alert alert-warning",
    style = "margin-bottom: 15px;",
    icon("exclamation-triangle"),
    tags$strong(" QC Alert: "),
    paste(
      vapply(seq_len(nrow(flagged)), function(i) {
        sprintf("%s at %d proteins (threshold: %d, 90-day mean: %d)",
                flagged$instrument[i],
                flagged$latest_proteins[i],
                flagged$threshold[i],
                flagged$mean_prot[i])
      }, ""),
      collapse = "; "
    )
  )
})
```

### Optional: Sidebar Alert Badge

In `R/ui.R`, add a small alert indicator next to the Instrument QC tab label:

```r
# This requires dynamic tab labels which bslib doesn't natively support.
# Alternative: show alert in sidebar when in core facility mode.
if (is_core_facility) tagList(
  uiOutput("sidebar_qc_alert")
)
```

In `server_facility.R`:

```r
output$sidebar_qc_alert <- renderUI({
  invalidateLater(300000, session)
  alerts <- tryCatch(cf_check_qc_alerts(cf_config$db_path), error = function(e) NULL)
  if (is.null(alerts)) return(NULL)
  n_flagged <- sum(alerts$alert, na.rm = TRUE)
  if (n_flagged == 0) return(NULL)

  tags$div(class = "alert alert-warning py-1 px-2 mt-2",
    style = "font-size: 0.8em;",
    icon("exclamation-triangle"),
    sprintf(" %d instrument%s below QC threshold",
            n_flagged, if (n_flagged > 1) "s" else ""))
})
```

### Files to Modify

| File | Changes |
|------|---------|
| `R/helpers_facility.R` | Add `cf_check_qc_alerts()` (~50 lines) |
| `R/server_facility.R` | Add alert banner render + sidebar alert (~30 lines) |
| `R/ui.R` | Add `uiOutput("qc_alert_banner")` in Instrument QC tab (~3 lines), sidebar alert (~3 lines) |

### Gotchas

- SQLite has no built-in `STDEV()` — compute in R after fetching values
- Need at least 3 QC runs for meaningful SD calculation
- Exclude the latest run from the mean/SD calculation to avoid self-comparison
- `invalidateLater(300000)` = every 5 minutes (avoids excessive DB queries)

---

## 8. End-to-End Testing

### Concept

Validate the complete core facility workflow with real or synthetic data. This is a test plan, not code to ship.

### Test Environment Setup

```bash
# 1. Create test directory
mkdir -p /tmp/delimp_test/.core_facility_test
export DELIMP_CORE_DIR=/tmp/delimp_test/.core_facility_test

# 2. Create minimal staff.yml
cat > /tmp/delimp_test/.core_facility_test/staff.yml << 'EOF'
staff:
  - name: "Test User"
    hpc_username: "testuser"
    ssh_key: "/tmp/fake_key"
    hpc_host: "localhost"
    slurm_account: "test-grp"
    slurm_partition: "low"
    lab: "Test Lab"
EOF

# 3. Create qc_config.yml
cat > /tmp/delimp_test/.core_facility_test/qc_config.yml << 'EOF'
instruments:
  - name: "timsTOF HT"
  - name: "timsTOF Ultra"
lc_methods:
  - name: "Evosep 100SPD"
  - name: "Evosep 60SPD"
report:
  facility_name: "Test Proteomics Core"
  logo_path: ""
EOF

# 4. Copy report template
cp report_template.qmd /tmp/delimp_test/.core_facility_test/

# 5. Seed test database
Rscript seed_test_db.R

# 6. Start app
Rscript -e "shiny::runApp('.', port=3838)"
```

### Test Checklist

#### A. Core Facility Detection & UI
- [ ] App starts with `is_core_facility = TRUE`
- [ ] "Search DB" tab appears in navigation
- [ ] "Instrument QC" tab appears in navigation
- [ ] Staff selector appears in New Search sidebar
- [ ] Template selector appears in sidebar
- [ ] Non-facility UI elements (manual SSH fields) are hidden

#### B. Job History (Search DB Tab)
- [ ] Seeded jobs appear in the table (20 rows from `seed_test_db.R`)
- [ ] Text search filters by analysis name
- [ ] Lab filter dropdown populated from DB
- [ ] Status filter works (completed/running/failed/queued)
- [ ] Staff filter dropdown populated from DB
- [ ] Instrument filter works
- [ ] LC method filter works
- [ ] Combined filters work together
- [ ] Table auto-refreshes every 30 seconds
- [ ] Row selection highlights correctly
- [ ] "Load Results" shows appropriate message for completed/non-completed jobs

#### C. Instrument QC Dashboard
- [ ] Seeded QC runs appear in trend plots (30 runs from `seed_test_db.R`)
- [ ] Protein trend plot shows both instruments
- [ ] Precursor trend plot renders
- [ ] TIC trend plot renders
- [ ] Instrument filter narrows to single instrument
- [ ] Date range filter works (30/90/180 days)
- [ ] Red dashed threshold lines (mean - 2SD) appear on protein trend
- [ ] QC runs table is sortable
- [ ] Hover tooltips show run details

#### D. Report Generation (requires DE results loaded)
- [ ] Load example data → Assign groups → Run pipeline
- [ ] Click "Generate Report" in Search DB tab (no job selected)
- [ ] Modal appears with metadata fields
- [ ] Fill in title, lab, instrument
- [ ] Click Generate → progress bar appears
- [ ] State `.rds` saved to state directory
- [ ] HTML report rendered (if Quarto available)
- [ ] HTML report opens in browser with:
  - [ ] Metadata table
  - [ ] QC bracket section (with seeded QC data)
  - [ ] Volcano plots (one per contrast)
  - [ ] DE summary table
  - [ ] Top DE proteins table
  - [ ] Normalization boxplot
  - [ ] Missing data barplot
  - [ ] Session info
- [ ] Download button provides the HTML file
- [ ] Report recorded in SQLite `reports` table

#### E. Report Generation from Job Selection
- [ ] Select a completed job row in Search DB
- [ ] Click "Generate Report" → uses job metadata (no modal)
- [ ] Report title/lab/instrument pre-filled from job record

#### F. QC Run Ingestion (Phase 2)
- [ ] Upload a file named "HeLa_QC_test.parquet" → auto-detection banner appears
- [ ] Click "Record QC Metrics" → success notification
- [ ] Check Instrument QC tab → new data point appears
- [ ] Alternatively: check "Record as QC run" checkbox → select instrument → run pipeline → QC recorded

#### G. Templates
- [ ] Save current search settings as template
- [ ] Template appears in template selector dropdown
- [ ] Loading template applies settings correctly
- [ ] Multiple templates coexist

#### H. Staff Selector
- [ ] Selecting "Test User" auto-fills SSH fields
- [ ] Connection status shown (will fail without real SSH — expected)

#### I. QC Alerts (Phase 2)
- [ ] Insert a low-protein QC run via SQLite
- [ ] Refresh Instrument QC tab → alert banner appears
- [ ] Sidebar alert indicator shows

#### J. Audit Log (Phase 2)
- [ ] Generate a report → audit log entry created
- [ ] Save a template → audit log entry created
- [ ] View audit log in Search DB tab (collapsible section)

### Automated Test Script

Consider creating `tests/test_facility.R` for non-interactive validation:

```r
# tests/test_facility.R — Run with: Rscript tests/test_facility.R
library(testthat)
library(DBI)
library(RSQLite)

source("R/helpers_facility.R")

test_that("cf_init_db creates all tables", {
  db <- dbConnect(SQLite(), ":memory:")
  cf_init_db(db)

  tables <- dbListTables(db)
  expect_true("searches" %in% tables)
  expect_true("qc_runs" %in% tables)
  expect_true("reports" %in% tables)
  expect_true("templates" %in% tables)
  expect_true("audit_log" %in% tables)  # Phase 2

  dbDisconnect(db)
})

test_that("cf_record_search inserts correctly", {
  db_path <- tempfile(fileext = ".db")
  db <- dbConnect(SQLite(), db_path)
  cf_init_db(db)
  dbDisconnect(db)

  cf_record_search(db_path, list(
    analysis_name = "test_search",
    submitted_by = "Test User",
    lab = "Test Lab",
    instrument = "timsTOF HT"
  ))

  db <- dbConnect(SQLite(), db_path)
  rows <- dbGetQuery(db, "SELECT * FROM searches")
  expect_equal(nrow(rows), 1)
  expect_equal(rows$analysis_name, "test_search")
  dbDisconnect(db)
  unlink(db_path)
})

test_that("cf_check_qc_alerts flags low protein count", {
  db_path <- tempfile(fileext = ".db")
  db <- dbConnect(SQLite(), db_path)
  cf_init_db(db)

  # Insert 10 "normal" runs
  for (i in 1:10) {
    dbExecute(db, "INSERT INTO qc_runs (run_name, instrument, run_date, n_proteins)
      VALUES (?1, ?2, ?3, ?4)",
      params = list(paste0("QC_", i), "TestInst",
                    as.character(Sys.time() - (11 - i) * 86400),
                    6000 + sample(-200:200, 1)))
  }

  # Insert one very low run (latest)
  dbExecute(db, "INSERT INTO qc_runs (run_name, instrument, run_date, n_proteins)
    VALUES ('QC_LOW', 'TestInst', ?1, 4000)",
    params = list(as.character(Sys.time())))

  dbDisconnect(db)

  alerts <- cf_check_qc_alerts(db_path)
  expect_true(any(alerts$alert))
  unlink(db_path)
})

test_that("cf_ingest_qc_metrics records run", {
  db_path <- tempfile(fileext = ".db")
  db <- dbConnect(SQLite(), db_path)
  cf_init_db(db)
  dbDisconnect(db)

  # Create a fake y_protein matrix
  y <- matrix(rnorm(500 * 6, mean = 20, sd = 3), nrow = 500, ncol = 6)
  rownames(y) <- paste0("P", 1:500)
  colnames(y) <- paste0("Sample", 1:6)

  cf_ingest_qc_metrics(db_path, "test_hela.parquet", "TestInst", y)

  db <- dbConnect(SQLite(), db_path)
  runs <- dbGetQuery(db, "SELECT * FROM qc_runs")
  expect_equal(nrow(runs), 1)
  expect_equal(runs$n_proteins, 500)
  expect_equal(runs$instrument, "TestInst")
  dbDisconnect(db)
  unlink(db_path)
})
```

### Files to Create

| File | Purpose |
|------|---------|
| `tests/test_facility.R` | Automated tests for facility helpers (~100 lines) |

---

## Summary of All File Changes

### New Files

| File | Lines (est.) | Purpose |
|------|-------------|---------|
| `tests/test_facility.R` | ~100 | Automated tests for facility helpers |

### Modified Files

| File | Feature(s) | Lines Added (est.) |
|------|-----------|-------------------|
| `R/helpers_facility.R` | QC ingestion, audit log, QC alerts | ~140 |
| `R/server_facility.R` | QC ingestion UI, template auto-apply, audit log, QC alerts, report comparison | ~250 |
| `R/server_session.R` | HF state download | ~70 |
| `R/ui.R` | QC ingestion checkbox, template selector in search, audit log, QC alert banners, compare button | ~40 |
| `report_template.qmd` | GSEA section, MOFA section, logo/header, new params | ~100 |
| `seed_test_db.R` | Add audit_log test data | ~15 |

**Estimated total**: ~715 lines across 7 files, ~28 hours of implementation

### New Packages Required

None. All features use packages already in the base image (DBI, RSQLite, jsonlite, DT, plotly, quarto). GSEA section in report template uses `enrichplot`/`clusterProfiler` which are already installed for the GSEA tab — but wrapped in `requireNamespace()` guards for graceful degradation.
