#!/usr/bin/env Rscript
# ==============================================================================
#  seed_test_db.R — Populate core facility test DB with synthetic data
#  Run: Rscript seed_test_db.R
# ==============================================================================

library(DBI)
library(RSQLite)

db_path <- file.path(Sys.getenv("DELIMP_CORE_DIR",
  unset = file.path(dirname(getwd()), ".core_facility_test")),
  "delimp.db")

# Also check current-dir-relative path
if (!dir.exists(dirname(db_path))) {
  db_path <- file.path(".core_facility_test", "delimp.db")
}

cat("DB path:", db_path, "\n")

# Source helpers to get cf_init_db
source("R/helpers_facility.R")

db <- dbConnect(SQLite(), db_path)
cf_init_db(db)

set.seed(42)

# ==============================================================================
#  Staff & Labs
# ==============================================================================
staff <- c("Brett Phinney", "Test Analyst")
labs <- c(
  "Bhatt Lab", "Bhatt Lab", "Bhatt Lab",
  "Tagkopoulos Lab", "Tagkopoulos Lab",
  "Bhatt Lab", "Smith Lab", "Smith Lab",
  "Chen Lab", "Chen Lab",
  "Tagkopoulos Lab", "Smith Lab",
  "Bhatt Lab", "Chen Lab",
  "Smith Lab", "Smith Lab",
  "Chen Lab", "Bhatt Lab",
  "Tagkopoulos Lab", "Chen Lab"
)
instruments <- c("timsTOF HT", "timsTOF Ultra")
lc_methods <- c("Evosep 100SPD", "Evosep 60SPD", "Evosep 30SPD",
                 "nanoElute 120min", "nanoElute 60min")
search_modes <- c("libfree", "libfree", "libfree", "phospho", "library")
organisms <- c("Homo sapiens", "Homo sapiens", "Mus musculus",
               "Homo sapiens", "Rattus norvegicus")
projects <- c(
  "AP-MS Screen 2026Q1", "AP-MS Screen 2026Q1",
  "Phospho Timecourse", "Phospho Timecourse",
  "Liver Proteome Atlas", "Liver Proteome Atlas",
  "Drug Response Panel", "Drug Response Panel",
  "AP-MS Screen 2026Q1", "Cardiac Remodeling",
  "Cardiac Remodeling", "Liver Proteome Atlas",
  "Drug Response Panel", "AP-MS Screen 2026Q1",
  "Phospho Timecourse", "Cardiac Remodeling",
  "AP-MS Screen 2026Q1", "Drug Response Panel",
  "Liver Proteome Atlas", "Phospho Timecourse"
)

statuses <- c("completed", "completed", "completed", "completed",
              "completed", "completed", "completed", "running",
              "completed", "completed",
              "completed", "completed", "failed", "completed",
              "queued", "completed", "completed", "completed",
              "completed", "running")

# ==============================================================================
#  Insert 20 synthetic searches (spanning ~30 days)
# ==============================================================================
cat("Inserting searches...\n")

base_time <- as.POSIXct("2026-01-20 09:00:00")

for (i in seq_along(labs)) {
  submit_time <- base_time + (i - 1) * 86400 + sample(0:28800, 1)  # ~1 per day, random offset
  complete_time <- if (statuses[i] %in% c("completed", "failed")) {
    submit_time + sample(3600:14400, 1)  # 1-4 hours later
  } else NA

  n_raw <- sample(4:24, 1)
  n_prot <- if (statuses[i] == "completed") sample(4500:8500, 1) else NA
  n_pep  <- if (statuses[i] == "completed") as.integer(n_prot * runif(1, 4.5, 6.0)) else NA
  n_prec <- if (statuses[i] == "completed") as.integer(n_pep * runif(1, 1.3, 1.8)) else NA

  inst <- sample(instruments, 1)
  lc   <- sample(lc_methods, 1)
  mode <- sample(search_modes, 1)
  org  <- sample(organisms, 1)
  staff_member <- sample(staff, 1)

  fasta <- switch(org,
    "Homo sapiens" = "UP000005640_human_one_per_gene.fasta",
    "Mus musculus" = "UP000000589_mouse_one_per_gene.fasta",
    "Rattus norvegicus" = "UP000002494_rat_one_per_gene.fasta"
  )

  analysis_name <- sprintf("%s_%s_%s",
    gsub(" ", "_", labs[i]),
    gsub(" ", "", mode),
    format(as.Date(submit_time), "%Y%m%d"))

  dbExecute(db, "
    INSERT INTO searches (analysis_name, submitted_by, lab, instrument,
                          lc_method, project,
                          organism, fasta_file, n_raw_files, search_mode,
                          submitted_at, completed_at, status, backend,
                          n_proteins, n_peptides, n_precursors,
                          slurm_job_id, output_dir)
    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9, ?10, ?11, ?12, ?13, ?14,
            ?15, ?16, ?17, ?18, ?19)",
    params = list(
      analysis_name, staff_member, labs[i], inst,
      lc, projects[i],
      org, fasta, n_raw, mode,
      format(submit_time, "%Y-%m-%d %H:%M:%S"),
      if (!is.na(complete_time)) format(complete_time, "%Y-%m-%d %H:%M:%S") else NA,
      statuses[i], sample(c("hpc", "docker"), 1, prob = c(0.7, 0.3)),
      n_prot, n_pep, n_prec,
      if (statuses[i] != "queued") sprintf("%d", sample(100000:999999, 1)) else NA,
      sprintf("/quobyte/proteomics-grp/de-limp/output/%s", analysis_name)
    )
  )
}

cat(sprintf("  Inserted %d searches\n",
    dbGetQuery(db, "SELECT COUNT(*) FROM searches")[[1]]))

# ==============================================================================
#  Insert QC runs — 3 per instrument, spanning 90 days
# ==============================================================================
cat("Inserting QC runs...\n")

qc_proteins <- list(
  "timsTOF HT"   = c(mean = 5200, sd = 300),
  "timsTOF Ultra" = c(mean = 6800, sd = 350)
)

for (inst in instruments) {
  params <- qc_proteins[[inst]]
  n_qc <- 15  # 15 QC runs per instrument over ~90 days

  for (j in seq_len(n_qc)) {
    run_date <- base_time - 90 * 86400 + j * (90 * 86400 / n_qc) + sample(-3600:3600, 1)

    n_prot  <- as.integer(rnorm(1, params["mean"], params["sd"]))
    n_pep   <- as.integer(n_prot * runif(1, 4.5, 5.5))
    n_prec  <- as.integer(n_pep * runif(1, 1.3, 1.7))
    ms1_tic <- runif(1, 1e9, 5e9)
    med_cv  <- runif(1, 8, 22)

    # Generate fake protein intensities JSON (top 20 for demo)
    prot_names <- paste0("P", sprintf("%05d", sample(10000:99999, 20)))
    prot_vals  <- round(rlnorm(20, meanlog = 25, sdlog = 2))
    intensities_json <- jsonlite::toJSON(setNames(as.list(prot_vals), prot_names),
                                          auto_unbox = TRUE)

    run_name <- sprintf("QC_HeLa_%s_%s",
      gsub(" ", "", inst),
      format(as.Date(run_date), "%Y%m%d"))

    dbExecute(db, "
      INSERT INTO qc_runs (run_name, instrument, run_date, n_proteins,
                            n_peptides, n_precursors, median_ms1_tic,
                            median_cv, protein_intensities)
      VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
      params = list(
        run_name, inst,
        format(run_date, "%Y-%m-%d %H:%M:%S"),
        n_prot, n_pep, n_prec, ms1_tic, round(med_cv, 1),
        as.character(intensities_json)
      )
    )
  }
}

cat(sprintf("  Inserted %d QC runs\n",
    dbGetQuery(db, "SELECT COUNT(*) FROM qc_runs")[[1]]))

# ==============================================================================
#  Insert 3 sample reports
# ==============================================================================
cat("Inserting reports...\n")

for (k in 1:3) {
  report_id <- sprintf("test-report-%04d", k)
  title <- c("Bhatt Lab AP-MS Feb 2026",
             "Smith Lab Phospho Analysis",
             "Chen Lab Liver Atlas")[k]
  lab <- c("Bhatt Lab", "Smith Lab", "Chen Lab")[k]
  inst <- sample(instruments, 1)

  dbExecute(db, "
    INSERT INTO reports (id, title, created_by, lab, instrument,
                         n_proteins, n_contrasts, contrast_names,
                         state_path)
    VALUES (?1, ?2, ?3, ?4, ?5, ?6, ?7, ?8, ?9)",
    params = list(
      report_id, title, "Brett Phinney", lab, inst,
      sample(5000:7000, 1), sample(2:6, 1),
      paste(paste0("Group", 1:sample(2:4, 1), "_vs_Control"), collapse = "; "),
      sprintf("/tmp/delimp_states/%s.rds", report_id)
    )
  )
}

cat(sprintf("  Inserted %d reports\n",
    dbGetQuery(db, "SELECT COUNT(*) FROM reports")[[1]]))

# ==============================================================================
#  Insert 2 sample templates
# ==============================================================================
cat("Inserting templates...\n")

dbExecute(db, "
  INSERT INTO templates (name, type, created_by, config_json, notes)
  VALUES (?1, ?2, ?3, ?4, ?5)",
  params = list(
    "Standard DIA Search", "search_preset", "Brett Phinney",
    jsonlite::toJSON(list(
      cpus = 64, mem_gb = 512, time_hours = 12,
      partition = "high", account = "genome-center-grp",
      search_mode = "libfree", enzyme = "K*,R*",
      missed_cleavages = 1, max_var_mods = 1,
      normalization = "on"
    ), auto_unbox = TRUE),
    "Standard library-free search for most experiments"
  )
)

dbExecute(db, "
  INSERT INTO templates (name, type, created_by, config_json, notes)
  VALUES (?1, ?2, ?3, ?4, ?5)",
  params = list(
    "Phospho Search", "search_preset", "Brett Phinney",
    jsonlite::toJSON(list(
      cpus = 64, mem_gb = 512, time_hours = 18,
      partition = "high", account = "genome-center-grp",
      search_mode = "phospho", enzyme = "K*,R*",
      missed_cleavages = 2, max_var_mods = 3,
      normalization = "on"
    ), auto_unbox = TRUE),
    "Phosphoproteomics with STY mod, extra missed cleavages"
  )
)

cat(sprintf("  Inserted %d templates\n",
    dbGetQuery(db, "SELECT COUNT(*) FROM templates")[[1]]))

# ==============================================================================
#  Summary
# ==============================================================================
cat("\n=== Database Summary ===\n")
for (tbl in c("searches", "qc_runs", "reports", "templates")) {
  n <- dbGetQuery(db, sprintf("SELECT COUNT(*) FROM %s", tbl))[[1]]
  cat(sprintf("  %-15s %d rows\n", tbl, n))
}

# Show project distribution
cat("\n=== Projects by Staff + Lab ===\n")
proj_summary <- dbGetQuery(db, "
  SELECT submitted_by, lab, project, COUNT(*) as n
  FROM searches
  WHERE project IS NOT NULL AND project != ''
  GROUP BY submitted_by, lab, project
  ORDER BY submitted_by, lab, project")
print(proj_summary)

# Show instrument distribution
cat("\n=== Instrument Distribution ===\n")
inst_summary <- dbGetQuery(db, "
  SELECT instrument, lc_method, COUNT(*) as n
  FROM searches
  GROUP BY instrument, lc_method
  ORDER BY instrument, lc_method")
print(inst_summary)

dbDisconnect(db)
cat("\nDone! DB seeded at:", db_path, "\n")
