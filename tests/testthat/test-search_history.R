# Tests for search history helper functions (helpers_search.R)

test_that("search_history_path returns a valid path", {
  path <- search_history_path()
  expect_true(is.character(path))
  expect_true(nzchar(path))
  expect_true(grepl("search_history\\.csv$", path))
})

test_that("search_history_read returns empty data.frame for missing file", {
  tmp <- tempfile(fileext = ".csv")
  result <- search_history_read(tmp)
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 0)
})

test_that("record_search creates CSV with correct headers on first write", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp, paste0(tmp, ".lock"))), add = TRUE)

  entry <- list(
    timestamp = "2024-06-15 10:30:00",
    completed_at = NA,
    user = "testuser",
    search_name = "test_search",
    backend = "hpc",
    search_mode = "libfree",
    parallel = TRUE,
    n_files = 10,
    fasta_files = "human.fasta",
    fasta_seq_count = 20000,
    normalization = "on",
    enzyme = "K*,R*",
    mass_acc_mode = "auto",
    mass_acc = NA,
    mass_acc_ms1 = NA,
    scan_window = 6,
    mbr = TRUE,
    extra_cli_flags = "",
    output_dir = "/tmp/test_output",
    job_id = "12345",
    status = "submitted",
    duration_min = NA,
    speclib_cached = FALSE,
    imported_from_log = FALSE,
    app_version = "3.2.0",
    notes = ""
  )

  record_search(entry, path = tmp)
  result <- search_history_read(tmp)

  expect_equal(nrow(result), 1)
  expect_equal(names(result), search_history_headers)
  expect_equal(result$search_name, "test_search")
  expect_equal(result$backend, "hpc")
  expect_equal(result$output_dir, "/tmp/test_output")
  expect_equal(result$status, "submitted")
})

test_that("record_search appends rows on subsequent writes", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp, paste0(tmp, ".lock"))), add = TRUE)

  base_entry <- list(
    timestamp = "2024-06-15 10:30:00", completed_at = NA,
    user = "testuser", search_name = "search_1", backend = "hpc",
    search_mode = "libfree", parallel = FALSE, n_files = 5,
    fasta_files = "human.fasta", fasta_seq_count = 20000,
    normalization = "on", enzyme = "K*,R*", mass_acc_mode = "auto",
    mass_acc = NA, mass_acc_ms1 = NA, scan_window = 6, mbr = TRUE,
    extra_cli_flags = "", output_dir = "/tmp/out1", job_id = "111",
    status = "submitted", duration_min = NA, speclib_cached = FALSE,
    imported_from_log = FALSE, app_version = "3.2.0", notes = ""
  )

  record_search(base_entry, path = tmp)

  entry2 <- base_entry
  entry2$search_name <- "search_2"
  entry2$output_dir <- "/tmp/out2"
  entry2$job_id <- "222"
  entry2$backend <- "docker"
  record_search(entry2, path = tmp)

  result <- search_history_read(tmp)
  expect_equal(nrow(result), 2)
  expect_equal(result$search_name, c("search_1", "search_2"))
  expect_equal(result$backend, c("hpc", "docker"))
})

test_that("update_search_status modifies correct row by output_dir", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp, paste0(tmp, ".lock"))), add = TRUE)

  entry <- list(
    timestamp = "2024-06-15 10:00:00", completed_at = NA,
    user = "testuser", search_name = "my_search", backend = "hpc",
    search_mode = "libfree", parallel = FALSE, n_files = 8,
    fasta_files = "human.fasta", fasta_seq_count = 20000,
    normalization = "on", enzyme = "K*,R*", mass_acc_mode = "auto",
    mass_acc = NA, mass_acc_ms1 = NA, scan_window = 6, mbr = TRUE,
    extra_cli_flags = "", output_dir = "/data/search_001", job_id = "555",
    status = "submitted", duration_min = NA, speclib_cached = FALSE,
    imported_from_log = FALSE, app_version = "3.2.0", notes = ""
  )
  record_search(entry, path = tmp)

  update_search_status(
    output_dir = "/data/search_001",
    status = "completed",
    completed_at = "2024-06-15 12:30:00",
    duration_min = 150.5,
    path = tmp
  )

  result <- search_history_read(tmp)
  expect_equal(nrow(result), 1)
  expect_equal(result$status, "completed")
  expect_equal(result$completed_at, "2024-06-15 12:30:00")
  expect_equal(result$duration_min, 150.5)
})

test_that("update_search_status doesn't modify other rows", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp, paste0(tmp, ".lock"))), add = TRUE)

  base <- list(
    timestamp = "2024-06-15 10:00:00", completed_at = NA,
    user = "testuser", search_name = "search_A", backend = "hpc",
    search_mode = "libfree", parallel = FALSE, n_files = 5,
    fasta_files = "human.fasta", fasta_seq_count = 20000,
    normalization = "on", enzyme = "K*,R*", mass_acc_mode = "auto",
    mass_acc = NA, mass_acc_ms1 = NA, scan_window = 6, mbr = TRUE,
    extra_cli_flags = "", output_dir = "/data/A", job_id = "100",
    status = "submitted", duration_min = NA, speclib_cached = FALSE,
    imported_from_log = FALSE, app_version = "3.2.0", notes = ""
  )
  record_search(base, path = tmp)

  entry_b <- base
  entry_b$search_name <- "search_B"
  entry_b$output_dir <- "/data/B"
  entry_b$job_id <- "200"
  record_search(entry_b, path = tmp)

  update_search_status("/data/B", "completed", "2024-06-15 14:00:00", 240, path = tmp)

  result <- search_history_read(tmp)
  expect_equal(result$status[1], "submitted")  # Row A unchanged
  expect_equal(result$status[2], "completed")  # Row B updated
  expect_equal(result$search_name[1], "search_A")
  expect_equal(result$search_name[2], "search_B")
})

test_that("update_search_status handles missing file gracefully", {
  tmp <- tempfile(fileext = ".csv")
  # Should not error
  expect_silent(
    update_search_status("/nonexistent", "completed", path = tmp)
  )
})

test_that("round-trip: record, read, verify all fields", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp, paste0(tmp, ".lock"))), add = TRUE)

  entry <- list(
    timestamp = "2024-06-15 09:00:00",
    completed_at = NA,
    user = "bphinney",
    search_name = "hela_digest",
    backend = "docker",
    search_mode = "phospho",
    parallel = FALSE,
    n_files = 3,
    fasta_files = "human.fasta, contam.fasta",
    fasta_seq_count = 21500,
    normalization = "off",
    enzyme = "K*,R*",
    mass_acc_mode = "manual",
    mass_acc = 14,
    mass_acc_ms1 = 12,
    scan_window = 8,
    mbr = FALSE,
    extra_cli_flags = "--min-corr 2.0",
    output_dir = "/results/hela_001",
    job_id = "delimp_hela_20240615",
    status = "submitted",
    duration_min = NA,
    speclib_cached = TRUE,
    imported_from_log = TRUE,
    app_version = "3.2.0",
    notes = "Test run with phospho"
  )

  record_search(entry, path = tmp)
  result <- search_history_read(tmp)

  expect_equal(nrow(result), 1)
  expect_equal(result$user, "bphinney")
  expect_equal(result$search_name, "hela_digest")
  expect_equal(result$backend, "docker")
  expect_equal(result$search_mode, "phospho")
  expect_equal(result$n_files, 3)
  expect_equal(result$fasta_files, "human.fasta, contam.fasta")
  expect_equal(result$fasta_seq_count, 21500)
  expect_equal(result$normalization, "off")
  expect_equal(result$mass_acc_mode, "manual")
  expect_equal(result$mass_acc, 14)
  expect_equal(result$mass_acc_ms1, 12)
  expect_equal(result$scan_window, 8)
  expect_equal(result$extra_cli_flags, "--min-corr 2.0")
  expect_equal(result$output_dir, "/results/hela_001")
  expect_equal(result$job_id, "delimp_hela_20240615")
  expect_equal(result$notes, "Test run with phospho")
})

test_that("backfill_search_history populates from job queue entries", {
  tmp <- tempfile(fileext = ".csv")
  on.exit(unlink(c(tmp, paste0(tmp, ".lock"))), add = TRUE)

  jobs <- list(
    list(
      name = "test_job_1",
      backend = "hpc",
      status = "completed",
      search_mode = "libfree",
      parallel = TRUE,
      n_files = 10,
      output_dir = "/data/job1",
      job_id = "12345",
      submitted_at = as.POSIXct("2024-06-15 10:00:00"),
      completed_at = as.POSIXct("2024-06-15 12:30:00"),
      speclib_cached = TRUE,
      search_settings = list(
        search_params = list(
          enzyme = "K*,R*", mass_acc_mode = "auto", mass_acc = NA,
          mass_acc_ms1 = NA, scan_window = 6, mbr = TRUE,
          extra_cli_flags = ""
        ),
        fasta_files = c("/data/human.fasta", "/data/contam.fasta"),
        fasta_seq_count = 20000,
        search_mode = "libfree",
        normalization = "on",
        n_raw_files = 10
      )
    ),
    list(
      name = "test_job_2",
      backend = "docker",
      status = "failed",
      search_mode = "library",
      parallel = FALSE,
      n_files = 5,
      output_dir = "/data/job2",
      job_id = "docker_abc",
      submitted_at = as.POSIXct("2024-06-16 08:00:00"),
      completed_at = NULL,
      speclib_cached = FALSE,
      search_settings = list(
        search_params = list(
          enzyme = "K*,R*", mass_acc_mode = "manual", mass_acc = 14,
          mass_acc_ms1 = 12, scan_window = 8, mbr = FALSE,
          extra_cli_flags = "--min-corr 2.0"
        ),
        fasta_files = "/data/mouse.fasta",
        fasta_seq_count = 17000,
        search_mode = "library",
        normalization = "off",
        n_raw_files = 5
      )
    )
  )

  n <- backfill_search_history(jobs, path = tmp, app_version = "3.2.0")
  expect_equal(n, 2)

  result <- search_history_read(tmp)
  expect_equal(nrow(result), 2)
  expect_equal(result$search_name, c("test_job_1", "test_job_2"))
  expect_equal(result$backend, c("hpc", "docker"))
  expect_equal(result$status, c("completed", "failed"))
  expect_equal(result$output_dir, c("/data/job1", "/data/job2"))
  expect_equal(result$duration_min[1], 150)
  expect_true(is.na(result$duration_min[2]))

  # Backfill again — should skip existing entries
  n2 <- backfill_search_history(jobs, path = tmp, app_version = "3.2.0")
  expect_equal(n2, 0)
  expect_equal(nrow(search_history_read(tmp)), 2)
})

test_that("backfill_search_history handles empty job list", {
  tmp <- tempfile(fileext = ".csv")
  n <- backfill_search_history(list(), path = tmp)
  expect_null(n)
})
