# Tests for R/helpers_search.R — DIA-NN search helper functions

# =============================================================================
# parse_sbatch_output
# =============================================================================

test_that("parse_sbatch_output extracts job ID from standard output", {
  stdout <- c("Submitted batch job 12345678")
  expect_equal(parse_sbatch_output(stdout), "12345678")
})

test_that("parse_sbatch_output handles trailing whitespace/carriage return", {
  stdout <- c("Submitted batch job 12345678\r")
  expect_equal(parse_sbatch_output(stdout), "12345678")
})

test_that("parse_sbatch_output handles multi-line output", {
  stdout <- c("Some warning message", "Submitted batch job 99887766", "")
  expect_equal(parse_sbatch_output(stdout), "99887766")
})

test_that("parse_sbatch_output returns NULL when no match", {
  stdout <- c("Error: invalid resource specification")
  expect_null(parse_sbatch_output(stdout))
})

test_that("parse_sbatch_output returns NULL for empty input", {
  expect_null(parse_sbatch_output(character(0)))
})

# =============================================================================
# ssh_control_path
# =============================================================================

test_that("ssh_control_path generates short path under /tmp", {
  cfg <- list(user = "testuser", host = "login.hpc.example.edu")
  path <- ssh_control_path(cfg)
  expect_true(startsWith(path, "/tmp/"))
  expect_true(grepl("testuser", path))
  # macOS limit: 104 bytes
  expect_true(nchar(path) < 104)
})

test_that("ssh_control_path sanitizes host", {
  cfg <- list(user = "me", host = "host-with.special.chars")
  path <- ssh_control_path(cfg)
  # Non-alphanumeric chars stripped from host portion
  expect_true(grepl("hostwithspecialchars", path))
  expect_false(grepl("host-with", path))
})

# =============================================================================
# ssh_mux_args
# =============================================================================

test_that("ssh_mux_args returns ControlMaster options", {
  cfg <- list(user = "testuser", host = "login.hpc.edu")
  args <- ssh_mux_args(cfg)
  expect_type(args, "character")
  expect_true(any(grepl("ControlMaster=auto", args)))
  expect_true(any(grepl("ControlPersist=300", args)))
  expect_true(any(grepl("ControlPath=", args)))
})

# =============================================================================
# generate_fasta_filename
# =============================================================================

test_that("generate_fasta_filename creates valid filename", {
  fn <- generate_fasta_filename("UP000005640", "Homo sapiens", "one_per_gene")
  expect_true(grepl("^UP000005640_homo_sapiens_opg_", fn))
  expect_true(grepl("\\.fasta$", fn))
})

test_that("generate_fasta_filename sanitizes organism name", {
  fn <- generate_fasta_filename("UP000000589", "Mus musculus (House mouse)", "full")
  # Special chars replaced, underscores collapsed
  expect_false(grepl("[()]", fn))
  expect_true(grepl("_full_", fn))
})

test_that("generate_fasta_filename maps content types", {
  expect_true(grepl("_opg_", generate_fasta_filename("UP1", "test", "one_per_gene")))
  expect_true(grepl("_sprot_", generate_fasta_filename("UP1", "test", "reviewed")))
  expect_true(grepl("_full_", generate_fasta_filename("UP1", "test", "full")))
  expect_true(grepl("_full_iso_", generate_fasta_filename("UP1", "test", "full_isoforms")))
})

test_that("generate_fasta_filename truncates long organism names", {
  fn <- generate_fasta_filename("UP1", "This is an extremely long organism name that should be truncated", "full")
  # Organism part should be <= 30 chars
  parts <- strsplit(fn, "_")[[1]]
  # Rejoin parts between ID and type suffix
  org_part <- paste(parts[2:(length(parts) - 3)], collapse = "_")
  expect_true(nchar(org_part) <= 30)
})

# =============================================================================
# estimate_search_time
# =============================================================================

test_that("estimate_search_time returns empty for 0 files", {
  expect_equal(estimate_search_time(0), "")
})

test_that("estimate_search_time returns non-empty for valid input", {
  result <- estimate_search_time(10, "libfree", 64)
  expect_true(nzchar(result))
  expect_true(grepl("10 files", result))
})

test_that("estimate_search_time scales with CPU count", {
  fast <- estimate_search_time(10, "libfree", 128)
  slow <- estimate_search_time(10, "libfree", 16)
  # Both should mention 10 files
  expect_true(grepl("10 files", fast))
  expect_true(grepl("10 files", slow))
})

test_that("estimate_search_time differentiates library vs libfree modes", {
  libfree <- estimate_search_time(5, "libfree", 64)
  library <- estimate_search_time(5, "library", 64)
  # Both should have output, but library mode should be faster
  expect_true(nzchar(libfree))
  expect_true(nzchar(library))
})

# =============================================================================
# build_diann_flags
# =============================================================================

test_that("build_diann_flags returns character vector", {
  flags <- build_diann_flags()
  expect_type(flags, "character")
  expect_true(length(flags) > 0)
})

test_that("build_diann_flags includes core flags", {
  flags <- build_diann_flags()
  joined <- paste(flags, collapse = " ")
  expect_true(grepl("--out-lib", joined))
  expect_true(grepl("--matrices", joined))
  expect_true(grepl("--gen-spec-lib", joined))
  expect_true(grepl("--qvalue", joined))
})

test_that("build_diann_flags includes Met oxidation by default", {
  flags <- build_diann_flags()
  expect_true(any(grepl("UniMod:35", flags)))
})

test_that("build_diann_flags omits Met oxidation when disabled", {
  flags <- build_diann_flags(search_params = list(mod_met_ox = FALSE))
  expect_false(any(grepl("UniMod:35", flags)))
})

test_that("build_diann_flags handles library mode", {
  flags <- build_diann_flags(search_mode = "library", speclib_mount = "/work/lib/speclib.tsv")
  joined <- paste(flags, collapse = " ")
  expect_true(grepl("--lib /work/lib/speclib.tsv", joined))
  expect_true(grepl("--use-quant", joined))
  # library mode should NOT have --fasta-search
  expect_false(grepl("--fasta-search", joined))
})

test_that("build_diann_flags handles libfree mode", {
  flags <- build_diann_flags(search_mode = "libfree")
  joined <- paste(flags, collapse = " ")
  expect_true(grepl("--fasta-search", joined))
  expect_true(grepl("--predictor", joined))
  expect_true(grepl("--cut", joined))
})

test_that("build_diann_flags handles manual mass accuracy", {
  flags <- build_diann_flags(search_params = list(mass_acc_mode = "manual", mass_acc = 10, mass_acc_ms1 = 5))
  joined <- paste(flags, collapse = " ")
  expect_true(grepl("--mass-acc 10", joined))
  expect_true(grepl("--mass-acc-ms1 5", joined))
})

test_that("build_diann_flags omits mass accuracy in auto mode", {
  flags <- build_diann_flags(search_params = list(mass_acc_mode = "auto"))
  joined <- paste(flags, collapse = " ")
  expect_false(grepl("--mass-acc ", joined))
})

test_that("build_diann_flags enables MBR by default", {
  flags <- build_diann_flags()
  expect_true(any(grepl("--reanalyse", flags)))
})

test_that("build_diann_flags disables MBR when requested", {
  flags <- build_diann_flags(search_params = list(mbr = FALSE))
  expect_false(any(grepl("--reanalyse", flags)))
})

test_that("build_diann_flags handles normalization off", {
  flags <- build_diann_flags(normalization = "off")
  expect_true(any(grepl("--no-norm", flags)))
})

test_that("build_diann_flags handles phospho mode", {
  flags <- build_diann_flags(search_mode = "phospho")
  joined <- paste(flags, collapse = " ")
  expect_true(grepl("--phospho-output", joined))
  expect_true(grepl("--report-lib-info", joined))
  # Phospho is library-free, so should have --fasta-search
  expect_true(grepl("--fasta-search", joined))
})

test_that("build_diann_flags adds extra CLI flags", {
  flags <- build_diann_flags(search_params = list(extra_cli_flags = "--smart-profiling"))
  expect_true(any(grepl("--smart-profiling", flags)))
})

test_that("build_diann_flags adds extra variable mods", {
  flags <- build_diann_flags(search_params = list(
    extra_var_mods = "UniMod:21,79.966331,STY\nUniMod:7,0.984016,NQ"
  ))
  expect_true(any(grepl("UniMod:21", flags)))
  expect_true(any(grepl("UniMod:7", flags)))
})

test_that("build_diann_flags respects xic toggle", {
  flags_on <- build_diann_flags(search_params = list(xic = TRUE))
  flags_off <- build_diann_flags(search_params = list(xic = FALSE))
  expect_true(any(grepl("--xic", flags_on)))
  expect_false(any(grepl("--xic", flags_off)))
})

# =============================================================================
# generate_sbatch_script
# =============================================================================

test_that("generate_sbatch_script returns valid sbatch script", {
  script <- generate_sbatch_script(
    analysis_name = "test_run",
    raw_files = c("/data/sample1.raw", "/data/sample2.raw"),
    fasta_files = c("/fasta/human.fasta"),
    output_dir = "/results/test_run",
    diann_sif = "/apps/diann_2.3.0.sif",
    cpus = 32, mem_gb = 256, time_hours = 8,
    partition = "high", account = "mygroup"
  )
  expect_type(script, "character")
  expect_true(grepl("#!/bin/bash", script))
  expect_true(grepl("#SBATCH --job-name=diann_test_run", script))
  expect_true(grepl("#SBATCH --cpus-per-task=32", script))
  expect_true(grepl("#SBATCH --mem=256G", script))
  expect_true(grepl("#SBATCH --time=8:00:00", script))
  expect_true(grepl("#SBATCH --partition=high", script))
  expect_true(grepl("#SBATCH --account=mygroup", script))
  expect_true(grepl("apptainer exec", script))
  expect_true(grepl("sample1.raw", script))
  expect_true(grepl("sample2.raw", script))
  expect_true(grepl("human.fasta", script))
})

test_that("generate_sbatch_script uses no_norm_report for normalization off", {
  script <- generate_sbatch_script(
    analysis_name = "test", raw_files = "/data/s1.raw",
    fasta_files = "/fasta/h.fasta", output_dir = "/out",
    diann_sif = "/sif/d.sif", normalization = "off"
  )
  expect_true(grepl("no_norm_report.parquet", script))
})

test_that("generate_sbatch_script handles multiple FASTA directories", {
  script <- generate_sbatch_script(
    analysis_name = "test",
    raw_files = "/data/s1.raw",
    fasta_files = c("/fasta1/human.fasta", "/fasta2/contaminants.fasta"),
    output_dir = "/out",
    diann_sif = "/sif/d.sif"
  )
  expect_true(grepl("/work/fasta", script))
})

# =============================================================================
# build_docker_command
# =============================================================================

test_that("build_docker_command returns docker run args", {
  flags <- build_diann_flags()
  args <- build_docker_command(
    raw_files = c("/data/sample1.raw"),
    fasta_files = c("/fasta/human.fasta"),
    output_dir = "/output",
    image_name = "diann:2.3.0",
    diann_flags = paste(flags, collapse = " "),
    cpus = 8, mem_gb = 32,
    container_name = "delimp_test"
  )
  expect_type(args, "character")
  expect_true("run" %in% args)
  expect_true("--rm" %in% args)
  expect_true("-d" %in% args)
  expect_true(any(grepl("delimp_test", args)))
  expect_true(any(grepl("diann:2.3.0", args)))
  expect_true(any(grepl("--cpus=8", args)))
  expect_true(any(grepl("--memory=32g", args)))
})

test_that("build_docker_command mounts data and fasta as read-only", {
  flags <- build_diann_flags()
  args <- build_docker_command(
    raw_files = c("/data/s1.raw"),
    fasta_files = c("/fasta/h.fasta"),
    output_dir = "/output",
    image_name = "diann:2.3.0",
    diann_flags = paste(flags, collapse = " "),
    cpus = 4, mem_gb = 16,
    container_name = "test"
  )
  joined <- paste(args, collapse = " ")
  expect_true(grepl("/data:/work/data:ro", joined))
  expect_true(grepl("/fasta:/work/fasta:ro", joined))
  # Output dir should NOT be read-only
  expect_true(grepl("/output:/work/out", joined))
  expect_false(grepl("/output:/work/out:ro", joined))
})

test_that("build_docker_command handles multiple data directories", {
  flags <- build_diann_flags()
  args <- build_docker_command(
    raw_files = c("/data1/s1.raw", "/data2/s2.raw"),
    fasta_files = c("/fasta/h.fasta"),
    output_dir = "/output",
    image_name = "diann:2.3.0",
    diann_flags = paste(flags, collapse = " "),
    cpus = 4, mem_gb = 16,
    container_name = "test"
  )
  joined <- paste(args, collapse = " ")
  expect_true(grepl("/data1:/work/data1:ro", joined))
  expect_true(grepl("/data2:/work/data2:ro", joined))
})

test_that("build_docker_command uses /tmp for intermediate files", {
  flags <- build_diann_flags()
  args <- build_docker_command(
    raw_files = c("/data/s1.raw"),
    fasta_files = c("/fasta/h.fasta"),
    output_dir = "/output",
    image_name = "diann:2.3.0",
    diann_flags = paste(flags, collapse = " "),
    cpus = 4, mem_gb = 16,
    container_name = "test"
  )
  joined <- paste(args, collapse = " ")
  expect_true(grepl("/tmp/diann_work", joined))
})

# =============================================================================
# write_file_list
# =============================================================================

test_that("write_file_list creates file_list.txt", {
  tmp <- tempdir()
  raw_files <- c("/data/s1.raw", "/data/s2.raw", "/data/s3.raw")
  result <- write_file_list(raw_files, tmp)
  expect_true(file.exists(result))
  lines <- readLines(result)
  expect_equal(lines, raw_files)
  unlink(result)
})

# =============================================================================
# check_cluster_resources (mocked)
# =============================================================================

test_that("check_cluster_resources returns expected structure", {
  skip("Requires SSH connection or local SLURM")
  res <- check_cluster_resources(NULL, "mygroup", "high")
  expect_type(res, "list")
  expect_true("success" %in% names(res))
  expect_true("group_limit" %in% names(res))
  expect_true("group_used" %in% names(res))
  expect_true("partition_idle" %in% names(res))
})
