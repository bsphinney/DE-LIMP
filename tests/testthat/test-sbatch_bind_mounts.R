# Tests for bind mount construction in generate_sbatch_script()
# Focuses on multi-directory data/fasta handling

# =============================================================================
# Single directory (standard case)
# =============================================================================

test_that("sbatch: single data dir maps to /work/data", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = c("/data/s1.raw", "/data/s2.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  expect_true(grepl("/data:/work/data", script))
  expect_true(grepl("--f /work/data/s1.raw", script))
  expect_true(grepl("--f /work/data/s2.raw", script))
})

# =============================================================================
# Multiple data directories
# =============================================================================

test_that("sbatch: multiple data dirs get numbered mounts", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = c("/data1/s1.raw", "/data2/s2.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  # Each dir gets its own numbered mount
  expect_true(grepl("/data1:/work/data1", script))
  expect_true(grepl("/data2:/work/data2", script))
  # Files map to correct mount points
  expect_true(grepl("--f /work/data1/s1.raw", script))
  expect_true(grepl("--f /work/data2/s2.raw", script))
})

test_that("sbatch: three data dirs all mounted correctly", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = c("/a/s1.raw", "/b/s2.raw", "/c/s3.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  expect_true(grepl("/a:/work/data1", script))
  expect_true(grepl("/b:/work/data2", script))
  expect_true(grepl("/c:/work/data3", script))
  expect_true(grepl("--f /work/data1/s1.raw", script))
  expect_true(grepl("--f /work/data2/s2.raw", script))
  expect_true(grepl("--f /work/data3/s3.raw", script))
})

test_that("sbatch: mixed dirs — files from same dir share mount", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = c("/data1/s1.raw", "/data2/s2.raw", "/data1/s3.raw"),
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  # Only 2 unique dirs
  expect_true(grepl("/data1:/work/data1", script))
  expect_true(grepl("/data2:/work/data2", script))
  # s1 and s3 both from /data1
  expect_true(grepl("--f /work/data1/s1.raw", script))
  expect_true(grepl("--f /work/data1/s3.raw", script))
  expect_true(grepl("--f /work/data2/s2.raw", script))
})

# =============================================================================
# Multiple FASTA directories
# =============================================================================

test_that("sbatch: multiple FASTA dirs get numbered mounts", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = "/data/s1.raw",
    fasta_files = c("/fasta1/human.fasta", "/fasta2/contaminants.fasta"),
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  expect_true(grepl("/fasta1:/work/fasta1", script))
  expect_true(grepl("/fasta2:/work/fasta2", script))
  expect_true(grepl("--fasta /work/fasta1/human.fasta", script))
  expect_true(grepl("--fasta /work/fasta2/contaminants.fasta", script))
})

# =============================================================================
# Speclib bind mount
# =============================================================================

test_that("sbatch: speclib dir mounted as /work/lib", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = "/data/s1.raw",
    fasta_files = "/fasta/h.fasta",
    speclib_path = "/lib/my_library.speclib",
    output_dir = "/out",
    diann_sif = "/sif/d.sif",
    search_mode = "library"
  )
  expect_true(grepl("/lib:/work/lib", script))
  expect_true(grepl("--lib /work/lib/my_library.speclib", script))
})

# =============================================================================
# Output dir mount
# =============================================================================

test_that("sbatch: output dir mounted as /work/out", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = "/data/s1.raw",
    fasta_files = "/fasta/h.fasta",
    output_dir = "/results/my_search",
    diann_sif = "/sif/d.sif"
  )
  expect_true(grepl("/results/my_search:/work/out", script))
  expect_true(grepl("--out /work/out/report.parquet", script))
})

# =============================================================================
# Log directory
# =============================================================================

test_that("sbatch: logs go to logs/ subdirectory", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = "/data/s1.raw",
    fasta_files = "/fasta/h.fasta",
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  expect_true(grepl("/out/logs/diann_", script))
})
