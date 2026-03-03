# Tests for generate_resume_launcher() in R/helpers_search.R

# =============================================================================
# Input validation
# =============================================================================

test_that("generate_resume_launcher rejects resume_from < 1", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  expect_error(generate_resume_launcher(0, "sbatch", paths))
})

test_that("generate_resume_launcher rejects resume_from > 5", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  expect_error(generate_resume_launcher(6, "sbatch", paths))
})

test_that("generate_resume_launcher rejects wrong number of script paths", {
  expect_error(generate_resume_launcher(1, "sbatch", paste0("/out/step", 1:3, ".sbatch")))
  expect_error(generate_resume_launcher(1, "sbatch", paste0("/out/step", 1:6, ".sbatch")))
})

test_that("generate_resume_launcher coerces numeric resume_from to integer", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  # Should not error — 3.0 coerces cleanly to 3L
  script <- generate_resume_launcher(3.0, "sbatch", paths)
  expect_type(script, "character")
})

# =============================================================================
# Resume from Step 1 (full restart)
# =============================================================================

test_that("generate_resume_launcher from Step 1 submits all 5 steps", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  script <- generate_resume_launcher(1, "/usr/bin/sbatch", paths)

  # No skipped steps
  expect_false(grepl("skipped", script))

  # All 5 steps submitted
  for (s in 1:5) {
    expect_true(grepl(sprintf("STEP%d:", s), script))
  }

  # First step has no --dependency
  expect_true(grepl("/usr/bin/sbatch /out/step1.sbatch", script))

  # Steps 2-5 have --dependency
  expect_true(grepl("--dependency=afterok:", script))
})

test_that("generate_resume_launcher from Step 1 uses correct sbatch path", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  script <- generate_resume_launcher(1, "/opt/slurm/bin/sbatch", paths)
  expect_true(grepl("/opt/slurm/bin/sbatch", script))
})

# =============================================================================
# Resume from middle step
# =============================================================================

test_that("generate_resume_launcher from Step 3 skips Steps 1-2", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  script <- generate_resume_launcher(3, "sbatch", paths)

  # Steps 1-2 marked as skipped
  expect_true(grepl('echo "STEP1:skipped"', script))
  expect_true(grepl('echo "STEP2:skipped"', script))

  # Step 3 submitted without dependency (it's the first actual submission)
  expect_true(grepl("JOB3=\\$\\(sbatch /out/step3.sbatch", script))

  # Steps 4-5 have dependencies
  expect_true(grepl("--dependency=afterok:\\$JOB3_ID /out/step4.sbatch", script))
  expect_true(grepl("--dependency=afterok:\\$JOB4_ID /out/step5.sbatch", script))
})

test_that("generate_resume_launcher from Step 2 skips Step 1 only", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  script <- generate_resume_launcher(2, "sbatch", paths)

  expect_true(grepl('echo "STEP1:skipped"', script))
  expect_false(grepl('echo "STEP2:skipped"', script))

  # Step 2 submitted without dependency
  expect_true(grepl("JOB2=\\$\\(sbatch /out/step2.sbatch", script))
})

# =============================================================================
# Resume from Step 5 (only final report)
# =============================================================================

test_that("generate_resume_launcher from Step 5 skips Steps 1-4", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  script <- generate_resume_launcher(5, "sbatch", paths)

  for (s in 1:4) {
    expect_true(grepl(sprintf('echo "STEP%d:skipped"', s), script))
  }

  # Step 5 submitted without dependency
  expect_true(grepl("JOB5=\\$\\(sbatch /out/step5.sbatch", script))
  # No dependency flags at all
  expect_false(grepl("--dependency", script))
})

# =============================================================================
# Script structure
# =============================================================================

test_that("generate_resume_launcher produces valid bash script", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  script <- generate_resume_launcher(3, "sbatch", paths)

  # Starts with shebang
  expect_true(grepl("^#!/bin/bash", script))
  # Has set -e for fail-fast
  expect_true(grepl("set -e", script))
})

test_that("generate_resume_launcher dependency chain is correct", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  script <- generate_resume_launcher(2, "sbatch", paths)

  # Step 2 → JOB2, Step 3 depends on JOB2, Step 4 depends on JOB3, etc.
  lines <- strsplit(script, "\n")[[1]]

  # Find submitted steps and verify chain
  submit_lines <- grep("--dependency=afterok:", lines, value = TRUE)
  expect_length(submit_lines, 3)  # Steps 3, 4, 5 each depend on previous

  expect_true(grepl("afterok:\\$JOB2_ID", submit_lines[1]))
  expect_true(grepl("afterok:\\$JOB3_ID", submit_lines[2]))
  expect_true(grepl("afterok:\\$JOB4_ID", submit_lines[3]))
})

test_that("generate_resume_launcher outputs STEP lines for all 5 steps", {
  paths <- paste0("/out/step", 1:5, ".sbatch")
  script <- generate_resume_launcher(3, "sbatch", paths)
  lines <- strsplit(script, "\n")[[1]]

  # Every step gets a STEP output line
  for (s in 1:5) {
    step_echo <- grep(sprintf('echo "STEP%d:', s), lines)
    expect_true(length(step_echo) >= 1,
                info = sprintf("Missing STEP%d output line", s))
  }
})
