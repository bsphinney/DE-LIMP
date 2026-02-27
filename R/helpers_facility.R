# ==============================================================================
#  helpers_facility.R
#  Core Facility Mode — Pure utility functions (no Shiny reactivity)
#  Database initialization, HF upload, CV computation, report helpers.
# ==============================================================================

#' Initialize SQLite database with core facility tables
#' @param db DBI connection to SQLite database
cf_init_db <- function(db) {
  # Enable WAL mode for concurrent reads across multiple staff members
  DBI::dbExecute(db, "PRAGMA journal_mode=WAL")

  # -- Searches / Job Queue --
  DBI::dbExecute(db, "
    CREATE TABLE IF NOT EXISTS searches (
      id              INTEGER PRIMARY KEY AUTOINCREMENT,
      analysis_name   TEXT NOT NULL,
      submitted_by    TEXT NOT NULL,
      lab             TEXT,
      submitted_at    DATETIME DEFAULT CURRENT_TIMESTAMP,
      completed_at    DATETIME,
      status          TEXT DEFAULT 'queued',
      slurm_job_id    TEXT,
      container_id    TEXT,
      backend         TEXT,
      instrument      TEXT,
      lc_method       TEXT,
      project         TEXT,
      organism        TEXT,
      fasta_file      TEXT,
      n_raw_files     INTEGER,
      search_mode     TEXT,
      n_proteins      INTEGER,
      n_peptides      INTEGER,
      n_precursors    INTEGER,
      output_dir      TEXT,
      report_id       TEXT,
      notes           TEXT
    )
  ")

  # -- Migrations for existing databases --
  cols <- DBI::dbListFields(db, "searches")
  if (!"lc_method" %in% cols) {
    DBI::dbExecute(db, "ALTER TABLE searches ADD COLUMN lc_method TEXT")
  }
  if (!"project" %in% cols) {
    DBI::dbExecute(db, "ALTER TABLE searches ADD COLUMN project TEXT")
  }

  # -- QC Runs --
  DBI::dbExecute(db, "
    CREATE TABLE IF NOT EXISTS qc_runs (
      id                  INTEGER PRIMARY KEY AUTOINCREMENT,
      run_name            TEXT NOT NULL,
      instrument          TEXT NOT NULL,
      run_date            DATETIME,
      ingested_at         DATETIME DEFAULT CURRENT_TIMESTAMP,
      file_path           TEXT,
      n_proteins          INTEGER,
      n_peptides          INTEGER,
      n_precursors        INTEGER,
      median_ms1_tic      REAL,
      median_cv           REAL,
      protein_intensities TEXT,
      search_id           INTEGER REFERENCES searches(id),
      report_path         TEXT
    )
  ")

  # -- Shared Reports --
  DBI::dbExecute(db, "
    CREATE TABLE IF NOT EXISTS reports (
      id            TEXT PRIMARY KEY,
      title         TEXT,
      created_at    DATETIME DEFAULT CURRENT_TIMESTAMP,
      created_by    TEXT,
      lab           TEXT,
      instrument    TEXT,
      n_proteins    INTEGER,
      n_contrasts   INTEGER,
      contrast_names TEXT,
      search_id     INTEGER REFERENCES searches(id),
      qc_before_id  INTEGER REFERENCES qc_runs(id),
      qc_after_id   INTEGER REFERENCES qc_runs(id),
      html_path     TEXT,
      state_path    TEXT,
      hf_state_url  TEXT
    )
  ")

  # -- Templates --
  DBI::dbExecute(db, "
    CREATE TABLE IF NOT EXISTS templates (
      id            INTEGER PRIMARY KEY AUTOINCREMENT,
      name          TEXT NOT NULL,
      type          TEXT NOT NULL,
      created_by    TEXT,
      created_at    DATETIME DEFAULT CURRENT_TIMESTAMP,
      config_json   TEXT NOT NULL,
      notes         TEXT
    )
  ")

  # -- Indexes --
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_searches_status ON searches(status)")
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_searches_lab ON searches(lab)")
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_searches_submitted_by ON searches(submitted_by)")
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_searches_submitted_at ON searches(submitted_at)")
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_searches_project ON searches(project)")
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_qc_instrument ON qc_runs(instrument)")
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_qc_date ON qc_runs(run_date)")
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_reports_lab ON reports(lab)")
  DBI::dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_templates_type ON templates(type)")
}


#' Upload analysis state to HF dataset repo
#' @param state_path Path to local .rds file
#' @param analysis_id UUID for this analysis
#' @return URL to the uploaded file, or NULL on failure
upload_state_to_hf <- function(state_path, analysis_id) {
  hf_token <- Sys.getenv("HF_TOKEN", "")
  if (!nzchar(hf_token)) return(NULL)

  repo_id <- "brettsp/delimp-shared-analyses"
  filename <- paste0(analysis_id, ".rds")

  # Upload via HF API
  url <- paste0("https://huggingface.co/api/datasets/", repo_id,
                "/upload/main/", filename)

  resp <- tryCatch({
    httr2::request(url) |>
      httr2::req_headers(
        Authorization = paste("Bearer", hf_token),
        `Content-Type` = "application/octet-stream"
      ) |>
      httr2::req_body_file(state_path) |>
      httr2::req_method("PUT") |>
      httr2::req_perform()
  }, error = function(e) {
    message("HF upload error: ", e$message)
    NULL
  })

  if (!is.null(resp) && httr2::resp_status(resp) < 300) {
    paste0("https://huggingface.co/datasets/", repo_id,
           "/resolve/main/", filename)
  } else {
    NULL
  }
}


#' Compute median CV of top N proteins across recent QC runs for an instrument
#' @param db SQLite connection
#' @param instrument Character
#' @param n_runs Number of recent runs to include (default 10)
#' @param top_n Number of most abundant proteins to track (default 500)
#' @return Median CV percentage, or NA if insufficient data
compute_rolling_cv <- function(db, instrument, n_runs = 10, top_n = 500) {
  runs <- DBI::dbGetQuery(db, "
    SELECT id, run_date, protein_intensities FROM qc_runs
    WHERE instrument = ?1 AND protein_intensities IS NOT NULL
    ORDER BY run_date DESC LIMIT ?2",
    params = list(instrument, n_runs)
  )

  if (nrow(runs) < 3) return(NA)

  # Parse JSON intensities for each run
  all_quant <- lapply(seq_len(nrow(runs)), function(i) {
    vals <- tryCatch(
      jsonlite::fromJSON(runs$protein_intensities[i]),
      error = function(e) NULL
    )
    if (is.null(vals) || length(vals) == 0) return(NULL)
    data.frame(protein = names(vals), intensity = as.numeric(vals),
               run_id = runs$id[i], stringsAsFactors = FALSE)
  })
  combined <- do.call(rbind, Filter(Negate(is.null), all_quant))
  if (is.null(combined) || nrow(combined) == 0) return(NA)

  # Find proteins present in all runs
  protein_counts <- table(combined$protein)
  common_proteins <- names(protein_counts[protein_counts == nrow(runs)])

  if (length(common_proteins) < 50) return(NA)

  # Take top N by mean intensity
  subset_df <- combined[combined$protein %in% common_proteins, ]
  mean_int <- tapply(subset_df$intensity, subset_df$protein, mean, na.rm = TRUE)
  top_proteins <- names(sort(mean_int, decreasing = TRUE))[seq_len(min(top_n, length(mean_int)))]

  # Compute CV per protein across runs
  top_subset <- subset_df[subset_df$protein %in% top_proteins, ]
  cvs <- tapply(top_subset$intensity, top_subset$protein, function(x) {
    sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
  })

  median(cvs, na.rm = TRUE)
}


#' Record a search job in the SQLite database
#' @param db_path Path to SQLite database
#' @param params Named list of job parameters
cf_record_search <- function(db_path, params) {
  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(db))

  DBI::dbExecute(db, "
    INSERT INTO searches (analysis_name, submitted_by, lab, instrument,
                          lc_method, project,
                          organism, fasta_file, n_raw_files, search_mode,
                          slurm_job_id, container_id, backend, status, output_dir)
    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14, ?15)",
    params = list(
      params$analysis_name %||% "",
      params$submitted_by %||% "unknown",
      params$lab %||% "",
      params$instrument %||% "",
      params$lc_method %||% "",
      params$project %||% "",
      params$organism %||% "",
      params$fasta_file %||% "",
      params$n_raw_files %||% 0L,
      params$search_mode %||% "library-free",
      params$slurm_job_id %||% NA,
      params$container_id %||% NA,
      params$backend %||% "hpc",
      "running",
      params$output_dir %||% ""
    )
  )
}


#' Update search status in SQLite
#' @param db_path Path to SQLite database
#' @param job_id SLURM job ID or Docker container ID
#' @param status New status
#' @param n_proteins Optional protein count
#' @param n_peptides Optional peptide count
#' @param n_precursors Optional precursor count
cf_update_search_status <- function(db_path, job_id, status,
                                     n_proteins = NULL, n_peptides = NULL,
                                     n_precursors = NULL) {
  db <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  on.exit(DBI::dbDisconnect(db))

  if (!is.null(n_proteins)) {
    DBI::dbExecute(db, "
      UPDATE searches
      SET status = ?1, completed_at = CURRENT_TIMESTAMP,
          n_proteins = ?2, n_peptides = ?3, n_precursors = ?4
      WHERE slurm_job_id = ?5 OR container_id = ?5",
      params = list(status, n_proteins, n_peptides %||% NA,
                    n_precursors %||% NA, job_id)
    )
  } else {
    DBI::dbExecute(db, "
      UPDATE searches SET status = ?1,
      completed_at = CASE WHEN ?1 IN ('completed', 'failed') THEN CURRENT_TIMESTAMP ELSE completed_at END
      WHERE slurm_job_id = ?2 OR container_id = ?2",
      params = list(status, job_id)
    )
  }
}


#' Get bracketing QC runs for a given instrument and timestamp
#' @param db DBI connection
#' @param instrument Instrument name
#' @param timestamp Character timestamp
#' @return List with qc_before and qc_after data frames
cf_get_bracketing_qc <- function(db, instrument, timestamp) {
  qc_before <- DBI::dbGetQuery(db, "
    SELECT * FROM qc_runs
    WHERE instrument = ?1 AND run_date <= ?2
    ORDER BY run_date DESC LIMIT 1",
    params = list(instrument, timestamp)
  )

  qc_after <- DBI::dbGetQuery(db, "
    SELECT * FROM qc_runs
    WHERE instrument = ?1 AND run_date >= ?2
    ORDER BY run_date ASC LIMIT 1",
    params = list(instrument, timestamp)
  )

  list(before = qc_before, after = qc_after)
}


#' Get 90-day rolling statistics for an instrument
#' @param db DBI connection
#' @param instrument Instrument name
#' @return Data frame with mean_prot, sd_prot, n_runs
cf_get_rolling_stats <- function(db, instrument) {
  DBI::dbGetQuery(db, "
    SELECT AVG(n_proteins) as mean_prot,
           COUNT(*) as n_runs
    FROM qc_runs
    WHERE instrument = ?1
      AND run_date >= datetime('now', '-90 days')
      AND n_proteins IS NOT NULL",
    params = list(instrument)
  )
}


#' Get staff list from config
#' @param cf_config Core facility config list
#' @return Named character vector (name = name, suitable for selectInput)
cf_staff_names <- function(cf_config) {
  if (is.null(cf_config) || is.null(cf_config$staff) || is.null(cf_config$staff$staff)) {
    return(c())
  }
  nms <- vapply(cf_config$staff$staff, `[[`, "", "name")
  setNames(nms, nms)
}


#' Get instrument names from QC config
#' @param cf_config Core facility config list
#' @return Character vector of instrument names
cf_instrument_names <- function(cf_config) {
  if (is.null(cf_config) || is.null(cf_config$qc) || is.null(cf_config$qc$instruments)) {
    return(c())
  }
  vapply(cf_config$qc$instruments, `[[`, "", "name")
}


#' Get lab names from staff config
#' @param cf_config Core facility config list
#' @return Unique lab names
cf_lab_names <- function(cf_config) {
  if (is.null(cf_config) || is.null(cf_config$staff) || is.null(cf_config$staff$staff)) {
    return(c())
  }
  unique(vapply(cf_config$staff$staff, function(s) s$lab %||% "", ""))
}


#' Get LC method names from QC config
#' @param cf_config Core facility config list
#' @return Character vector of LC method names
cf_lc_method_names <- function(cf_config) {
  if (is.null(cf_config) || is.null(cf_config$qc) || is.null(cf_config$qc$lc_methods)) {
    return(c())
  }
  vapply(cf_config$qc$lc_methods, `[[`, "", "name")
}
