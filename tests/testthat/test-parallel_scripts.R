# Tests for generate_parallel_scripts() in R/helpers_search.R
# Focus: step_flags filtering, multi-directory bind mounts, structural correctness
#
# Return structure is a flat list:
#   step1_library, step2_firstpass, step3_assembly, step4_finalpass, step5_report

# =============================================================================
# Step flag filtering (prevents FASTA re-digest bug)
# =============================================================================

test_that("parallel scripts: Step 1 has --fasta-search and --predictor", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/human.fasta",
    output_dir = "/out/test",
    diann_sif = "/sif/diann.sif"
  )
  step1 <- result$step1_library
  expect_true(grepl("--fasta-search", step1))
  expect_true(grepl("--predictor", step1))
})

test_that("parallel scripts: Steps 2-5 do NOT have --fasta-search or --predictor", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/human.fasta",
    output_dir = "/out/test",
    diann_sif = "/sif/diann.sif"
  )

  for (step_name in c("step2_firstpass", "step3_assembly",
                       "step4_finalpass", "step5_report")) {
    step_script <- result[[step_name]]
    expect_false(grepl("--fasta-search", step_script),
                 info = sprintf("%s should NOT have --fasta-search", step_name))
    expect_false(grepl("--predictor", step_script),
                 info = sprintf("%s should NOT have --predictor", step_name))
  }
})

test_that("parallel scripts: step_flags exclude MBR", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif",
    search_params = list(mbr = TRUE)
  )
  # --reanalyse (MBR) should not appear in any step (parallel replaces it)
  for (step_name in names(result)) {
    if (!is.null(result[[step_name]])) {
      expect_false(grepl("--reanalyse", result[[step_name]]),
                   info = sprintf("%s should NOT have --reanalyse", step_name))
    }
  }
})

# =============================================================================
# Multi-directory bind mounts
# =============================================================================

test_that("parallel scripts: single data dir uses /work/data mount", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  # Step 3 (assembly) uses full bind mount with all data dirs
  step3 <- result$step3_assembly
  expect_true(grepl("/data:/work/data", step3))
  expect_false(grepl("/work/data1", step3))
})

test_that("parallel scripts: multiple data dirs use numbered /work/dataN mounts", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data1/s1.raw", "/data2/s2.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  # Assembly steps should have both data dirs
  step3 <- result$step3_assembly
  expect_true(grepl("/data1:/work/data1", step3))
  expect_true(grepl("/data2:/work/data2", step3))
})

test_that("parallel scripts: --f flags map to correct mount points", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data1/s1.raw", "/data2/s2.raw", "/data1/s3.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  step3 <- result$step3_assembly
  # s1 and s3 are in /data1 -> /work/data1
  expect_true(grepl("--f /work/data1/s1.raw", step3))
  expect_true(grepl("--f /work/data1/s3.raw", step3))
  # s2 is in /data2 -> /work/data2
  expect_true(grepl("--f /work/data2/s2.raw", step3))
})

# =============================================================================
# Speclib-provided mode (skip Step 1)
# =============================================================================

test_that("parallel scripts: Step 1 is NULL when speclib provided", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    speclib_path = "/lib/my.speclib",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  expect_null(result$step1_library)
})

test_that("parallel scripts: Steps 2-5 use user speclib when provided", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    speclib_path = "/lib/my.speclib",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  # Step 2 should reference the user-provided lib path (mounted)
  step2 <- result$step2_firstpass
  expect_true(grepl("--lib /work/lib/my.speclib", step2))
})

# =============================================================================
# Return structure
# =============================================================================

test_that("parallel scripts: returns expected structure", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  expect_type(result, "list")

  # 5 script entries (flat list)
  script_names <- c("step1_library", "step2_firstpass", "step3_assembly",
                     "step4_finalpass", "step5_report")
  for (nm in script_names) {
    expect_true(nm %in% names(result),
                info = sprintf("Missing script: %s", nm))
  }
})

# =============================================================================
# SBATCH headers
# =============================================================================

test_that("parallel scripts: Step 1 uses libpred resource settings", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif",
    libpred_cpus = 8, libpred_mem = 32
  )
  step1 <- result$step1_library
  expect_true(grepl("--cpus-per-task=8", step1))
  expect_true(grepl("--mem=32G", step1))
})

test_that("parallel scripts: assembly steps use assembly resources", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif",
    assembly_cpus = 128, assembly_mem = 1024
  )
  step3 <- result$step3_assembly
  expect_true(grepl("--cpus-per-task=128", step3))
  expect_true(grepl("--mem=1024G", step3))
})

test_that("parallel scripts: normalization off uses no_norm_report.parquet", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif",
    normalization = "off"
  )
  step5 <- result$step5_report
  expect_true(grepl("no_norm_report.parquet", step5))
})

# =============================================================================
# Step 3 does not re-digest FASTA (regression test for d2b9bc6)
# =============================================================================

test_that("parallel scripts: Step 3 uses --use-quant, not --fasta-search", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/human.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  step3 <- result$step3_assembly
  expect_true(grepl("--use-quant", step3))
  expect_false(grepl("--fasta-search", step3))
})

# =============================================================================
# Step 5 uses --matrices for final report
# =============================================================================

test_that("parallel scripts: Step 5 has --matrices for output", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  step5 <- result$step5_report
  expect_true(grepl("--matrices", step5))
})

# =============================================================================
# Array job specs
# =============================================================================

# =============================================================================
# Mass accuracy forced to manual (prevents DIA-NN warning with --use-quant)
# =============================================================================

test_that("parallel scripts: mass accuracy flags present even when mode is auto", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif",
    search_params = list(mass_acc_mode = "auto", mass_acc = 14, mass_acc_ms1 = 14,
                          scan_window = 6)
  )
  # Even with auto mode, parallel scripts should have fixed mass accuracy
  step3 <- result$step3_assembly
  expect_true(grepl("--mass-acc 14", step3),
              info = "Step 3 should have --mass-acc even in auto mode")
  expect_true(grepl("--mass-acc-ms1 14", step3),
              info = "Step 3 should have --mass-acc-ms1 even in auto mode")
  expect_true(grepl("--window 6", step3),
              info = "Step 3 should have --window even in auto mode")

  step5 <- result$step5_report
  expect_true(grepl("--mass-acc 14", step5),
              info = "Step 5 should have --mass-acc even in auto mode")
})

test_that("parallel scripts: mass accuracy flags use manual values when set", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif",
    search_params = list(mass_acc_mode = "manual", mass_acc = 10, mass_acc_ms1 = 5,
                          scan_window = 8)
  )
  step3 <- result$step3_assembly
  expect_true(grepl("--mass-acc 10", step3))
  expect_true(grepl("--mass-acc-ms1 5", step3))
  expect_true(grepl("--window 8", step3))
})

# =============================================================================
# Quant file verification blocks
# =============================================================================

test_that("parallel scripts: Step 3 has quant file verification block", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out/test_dir",
    diann_sif = "/sif/d.sif"
  )
  step3 <- result$step3_assembly
  expect_true(grepl("Verifying Step 2 quant files", step3),
              info = "Step 3 should verify Step 2 quant files")
  expect_true(grepl("quant_step2", step3),
              info = "Step 3 should check quant_step2 directory")
  expect_true(grepl("exit 1", step3),
              info = "Step 3 should abort if quant files missing")
  expect_true(grepl("file_list.txt", step3),
              info = "Step 3 should read from file_list.txt")
})

test_that("parallel scripts: Step 5 has quant file verification block", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out/test_dir",
    diann_sif = "/sif/d.sif"
  )
  step5 <- result$step5_report
  expect_true(grepl("Verifying Step 4 quant files", step5),
              info = "Step 5 should verify Step 4 quant files")
  expect_true(grepl("quant_step4", step5),
              info = "Step 5 should check quant_step4 directory")
})

test_that("parallel scripts: Steps 2/4 do NOT have quant verification", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  step2 <- result$step2_firstpass
  step4 <- result$step4_finalpass
  expect_false(grepl("Verifying.*quant files", step2),
               info = "Step 2 creates quant files, should not verify them")
  expect_false(grepl("Verifying.*quant files", step4),
               info = "Step 4 creates quant files, should not verify them")
})

# =============================================================================
# Array job specs
# =============================================================================

test_that("parallel scripts: Step 2 array spec matches file count", {
  result <- generate_parallel_scripts(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw", "/data/s3.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif",
    max_simultaneous = 10
  )
  step2 <- result$step2_firstpass
  # 3 files → array 0-2%10
  expect_true(grepl("--array=0-2%10", step2))
})
