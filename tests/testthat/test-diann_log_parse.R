# Tests for parse_diann_log() — DIA-NN log file parser

# Helper: write a temp log file with given command line
write_log <- function(cmd_line, version_line = "DIA-NN 2.0.1") {
  tmp <- tempfile(fileext = ".log")
  writeLines(c(
    version_line,
    "Thread 1",
    cmd_line
  ), tmp)
  tmp
}

# =============================================================================
# Basic functionality
# =============================================================================

test_that("parse_diann_log extracts numeric parameters", {
  log <- write_log("diann --mass-acc 24 --mass-acc-ms1 7 --window 13 --qvalue 0.05 --var-mods 2")
  result <- parse_diann_log(log)
  expect_true(result$success)
  expect_equal(result$params$mass_acc, 24)
  expect_equal(result$params$mass_acc_ms1, 7)
  expect_equal(result$params$scan_window, 13L)
  expect_equal(result$params$qvalue, 0.05)
  expect_equal(result$params$max_var_mods, 2L)
  expect_equal(result$params$mass_acc_mode, "manual")
})

test_that("parse_diann_log extracts peptide/precursor range parameters", {
  log <- write_log("diann --min-pep-len 8 --max-pep-len 25 --min-pr-mz 350 --max-pr-mz 1800")
  result <- parse_diann_log(log)
  expect_true(result$success)
  expect_equal(result$params$min_pep_len, 8L)
  expect_equal(result$params$max_pep_len, 25L)
  expect_equal(result$params$min_pr_mz, 350L)
  expect_equal(result$params$max_pr_mz, 1800L)
})

test_that("parse_diann_log extracts enzyme and missed cleavages", {
  log <- write_log("diann --cut K*,R* --missed-cleavages 2")
  result <- parse_diann_log(log)
  expect_true(result$success)
  expect_equal(result$params$enzyme, "K*,R*")
  expect_equal(result$params$missed_cleavages, 2L)
})

test_that("parse_diann_log extracts DIA-NN version", {
  log <- write_log("diann --qvalue 0.01", version_line = "DIA-NN 2.3.0")
  result <- parse_diann_log(log)
  expect_equal(result$version, "2.3.0")
})

test_that("parse_diann_log returns raw command line", {
  cmd <- "diann --qvalue 0.01 --reanalyse"
  log <- write_log(cmd)
  result <- parse_diann_log(log)
  expect_equal(result$command_line, cmd)
})

# =============================================================================
# Boolean flags
# =============================================================================

test_that("parse_diann_log detects boolean flags when present", {
  log <- write_log("diann --reanalyse --rt-profiling --xic --unimod4 --met-excision")
  result <- parse_diann_log(log)
  expect_true(result$params$mbr)
  expect_true(result$params$rt_profiling)
  expect_true(result$params$xic)
  expect_true(result$params$unimod4)
  expect_true(result$params$met_excision)
})

test_that("parse_diann_log sets boolean flags FALSE when absent", {
  log <- write_log("diann --qvalue 0.01")
  result <- parse_diann_log(log)
  expect_false(result$params$mbr)
  expect_false(result$params$rt_profiling)
  expect_false(result$params$xic)
  expect_false(result$params$unimod4)
  expect_false(result$params$met_excision)
})

test_that("parse_diann_log handles mixed boolean flags", {
  log <- write_log("diann --reanalyse --unimod4")
  result <- parse_diann_log(log)
  expect_true(result$params$mbr)
  expect_true(result$params$unimod4)
  expect_false(result$params$rt_profiling)
  expect_false(result$params$xic)
  expect_false(result$params$met_excision)
})

# =============================================================================
# Variable modifications
# =============================================================================

test_that("parse_diann_log detects methionine oxidation", {
  log <- write_log("diann --var-mod UniMod:35,15.994915,M")
  result <- parse_diann_log(log)
  expect_true(result$params$mod_met_ox)
  expect_false(result$params$mod_nterm_acetyl)
})

test_that("parse_diann_log detects N-term acetylation", {
  log <- write_log("diann --var-mod UniMod:1,42.010565,*n")
  result <- parse_diann_log(log)
  expect_false(result$params$mod_met_ox)
  expect_true(result$params$mod_nterm_acetyl)
})

test_that("parse_diann_log collects custom var-mods", {
  log <- write_log("diann --var-mod UniMod:35,15.994915,M --var-mod UniMod:21,79.966331,STY")
  result <- parse_diann_log(log)
  expect_true(result$params$mod_met_ox)
  expect_equal(result$params$extra_var_mods, "UniMod:21,79.966331,STY")
})

test_that("parse_diann_log handles multiple custom var-mods", {
  log <- write_log("diann --var-mod UniMod:21,79.966331,STY --var-mod UniMod:7,0.984016,NQ")
  result <- parse_diann_log(log)
  expect_equal(result$params$extra_var_mods,
    "UniMod:21,79.966331,STY\nUniMod:7,0.984016,NQ")
})

# =============================================================================
# Search mode detection
# =============================================================================

test_that("parse_diann_log detects library-free mode", {
  log <- write_log("diann --fasta-search --predictor --cut K*,R*")
  result <- parse_diann_log(log)
  expect_equal(result$search_mode, "libfree")
})

test_that("parse_diann_log detects library mode", {
  log <- write_log("diann --lib /path/to/lib.speclib --qvalue 0.01")
  result <- parse_diann_log(log)
  expect_equal(result$search_mode, "library")
})

test_that("parse_diann_log detects phospho mode", {
  log <- write_log("diann --fasta-search --phospho-output --var-mod UniMod:21,79.966331,STY")
  result <- parse_diann_log(log)
  expect_equal(result$search_mode, "phospho")
})

test_that("parse_diann_log treats --lib + --fasta + --cut as libfree (second-pass log)", {
  log <- write_log("diann --lib /path/to/pred.speclib --fasta /db/human.fasta --cut K*,R* --mass-acc 10")
  result <- parse_diann_log(log)
  expect_equal(result$search_mode, "libfree")
})

# =============================================================================
# Normalization detection
# =============================================================================

test_that("parse_diann_log detects --no-norm", {
  log <- write_log("diann --no-norm --qvalue 0.01")
  result <- parse_diann_log(log)
  expect_equal(result$normalization, "off")
})

test_that("parse_diann_log defaults normalization to on", {
  log <- write_log("diann --qvalue 0.01")
  result <- parse_diann_log(log)
  expect_equal(result$normalization, "on")
})

# =============================================================================
# Mass accuracy mode
# =============================================================================

test_that("parse_diann_log detects manual mass accuracy", {
  log <- write_log("diann --mass-acc 10 --mass-acc-ms1 5")
  result <- parse_diann_log(log)
  expect_equal(result$params$mass_acc_mode, "manual")
})

test_that("parse_diann_log defaults to auto mass accuracy when no flags", {
  log <- write_log("diann --qvalue 0.01")
  result <- parse_diann_log(log)
  expect_equal(result$params$mass_acc_mode, "auto")
})

# =============================================================================
# Extra flags and hardcoded params
# =============================================================================

test_that("parse_diann_log parses fr_mz and pr_charge into params (not extra_cli_flags)", {
  log <- write_log("diann --max-fr-mz 1800 --min-pr-charge 2")
  result <- parse_diann_log(log)
  expect_equal(result$params$max_fr_mz, 1800)
  expect_equal(result$params$min_pr_charge, 2)
  # These should NOT appear in extra_cli_flags
  if (!is.null(result$params$extra_cli_flags)) {
    expect_false(grepl("--max-fr-mz", result$params$extra_cli_flags))
    expect_false(grepl("--min-pr-charge", result$params$extra_cli_flags))
  }
})

test_that("parse_diann_log parses default fr_mz and pr_charge values", {
  log <- write_log("diann --min-pr-charge 1 --max-pr-charge 4 --min-fr-mz 200 --max-fr-mz 1200")
  result <- parse_diann_log(log)
  expect_equal(result$params$min_pr_charge, 1)
  expect_equal(result$params$max_pr_charge, 4)
  expect_equal(result$params$min_fr_mz, 200)
  expect_equal(result$params$max_fr_mz, 1200)
})

test_that("parse_diann_log puts unrecognized flags in extra_cli_flags", {
  log <- write_log("diann --proteoforms --pg-level 0 --ids-to-names")
  result <- parse_diann_log(log)
  # proteoforms and pg-level are now recognized params (not extra_cli_flags)
  expect_true(result$params$proteoforms)
  expect_equal(result$params$pg_level, 0)
  # ids-to-names is still unrecognized
  expect_true(grepl("--ids-to-names", result$params$extra_cli_flags))
})

test_that("parse_diann_log puts invalid enzyme in extra_cli_flags", {
  log <- write_log("diann --cut X,Y,Z")
  result <- parse_diann_log(log)
  expect_null(result$params$enzyme)
  expect_true(grepl("--cut X,Y,Z", result$params$extra_cli_flags))
})

# =============================================================================
# FASTA and raw file counting
# =============================================================================

test_that("parse_diann_log counts raw files", {
  log <- write_log("diann --f /data/s1.raw --f /data/s2.raw --f /data/s3.raw --qvalue 0.01")
  result <- parse_diann_log(log)
  expect_equal(result$n_raw_files, 3L)
})

test_that("parse_diann_log collects FASTA paths", {
  log <- write_log("diann --fasta /db/human.fasta --fasta /db/contam.fasta --fasta-search")
  result <- parse_diann_log(log)
  expect_equal(length(result$fasta_files), 2)
  expect_true("/db/human.fasta" %in% result$fasta_files)
  expect_true("/db/contam.fasta" %in% result$fasta_files)
})

# =============================================================================
# Binary path formats
# =============================================================================

test_that("parse_diann_log handles Linux binary path", {
  log <- write_log("/usr/local/bin/diann-linux --qvalue 0.01 --reanalyse")
  result <- parse_diann_log(log)
  expect_true(result$success)
  expect_equal(result$params$qvalue, 0.01)
  expect_true(result$params$mbr)
})

test_that("parse_diann_log handles Windows .exe binary", {
  log <- write_log("diann.exe --qvalue 0.01 --mass-acc 10")
  result <- parse_diann_log(log)
  expect_true(result$success)
  expect_equal(result$params$mass_acc, 10)
})

test_that("parse_diann_log handles versioned binary name", {
  log <- write_log("diann-2.0 --qvalue 0.01")
  result <- parse_diann_log(log)
  expect_true(result$success)
  expect_equal(result$params$qvalue, 0.01)
})

# =============================================================================
# Error handling
# =============================================================================

test_that("parse_diann_log fails on missing file", {
  result <- parse_diann_log("/nonexistent/path/diann.log")
  expect_false(result$success)
  expect_true(grepl("not found", result$message))
})

test_that("parse_diann_log fails on empty file", {
  tmp <- tempfile(fileext = ".log")
  writeLines(character(0), tmp)
  result <- parse_diann_log(tmp)
  expect_false(result$success)
})

test_that("parse_diann_log fails when no command line found", {
  tmp <- tempfile(fileext = ".log")
  writeLines(c("DIA-NN 2.0.1", "Thread 1", "Loading library..."), tmp)
  result <- parse_diann_log(tmp)
  expect_false(result$success)
  expect_true(grepl("No DIA-NN command line", result$message))
})

# =============================================================================
# Windows CRLF handling
# =============================================================================

test_that("parse_diann_log handles Windows CRLF line endings", {
  tmp <- tempfile(fileext = ".log")
  writeBin(charToRaw("DIA-NN 2.0.1\r\nThread 1\r\ndiann --mass-acc 10 --reanalyse\r\n"), tmp)
  result <- parse_diann_log(tmp)
  expect_true(result$success)
  expect_equal(result$params$mass_acc, 10)
  expect_true(result$params$mbr)
})

# =============================================================================
# Comprehensive realistic command line
# =============================================================================

test_that("parse_diann_log handles realistic full command line", {
  cmd <- paste(
    "diann",
    "--f /data/sample1.raw --f /data/sample2.raw --f /data/sample3.raw",
    "--fasta /db/UP000005640.fasta",
    "--fasta-search --predictor",
    "--out /results/report.parquet --out-lib /results/report-lib.parquet",
    "--gen-spec-lib --matrices --verbose 1",
    "--var-mod UniMod:35,15.994915,M --var-mod UniMod:1,42.010565,*n",
    "--qvalue 0.01 --var-mods 2",
    "--cut K*,R* --missed-cleavages 1",
    "--min-pep-len 7 --max-pep-len 30",
    "--min-pr-mz 300 --max-pr-mz 1800",
    "--min-pr-charge 1 --max-pr-charge 4",
    "--min-fr-mz 200 --max-fr-mz 1200",
    "--mass-acc 24 --mass-acc-ms1 7 --window 13",
    "--unimod4 --met-excision --xic",
    "--reanalyse --rt-profiling",
    "--proteoforms --pg-level 0 --ids-to-names"
  )
  log <- write_log(cmd)
  result <- parse_diann_log(log)

  expect_true(result$success)
  expect_equal(result$version, "2.0.1")
  expect_equal(result$search_mode, "libfree")
  expect_equal(result$normalization, "on")
  expect_equal(result$n_raw_files, 3L)
  expect_equal(length(result$fasta_files), 1)
  expect_equal(basename(result$fasta_files), "UP000005640.fasta")

  p <- result$params
  expect_equal(p$qvalue, 0.01)
  expect_equal(p$max_var_mods, 2L)
  expect_equal(p$mass_acc, 24)
  expect_equal(p$mass_acc_ms1, 7)
  expect_equal(p$scan_window, 13L)
  expect_equal(p$mass_acc_mode, "manual")
  expect_equal(p$enzyme, "K*,R*")
  expect_equal(p$missed_cleavages, 1L)
  expect_equal(p$min_pep_len, 7L)
  expect_equal(p$max_pep_len, 30L)
  expect_equal(p$min_pr_mz, 300L)
  expect_equal(p$max_pr_mz, 1800L)
  expect_true(p$mod_met_ox)
  expect_true(p$mod_nterm_acetyl)
  expect_true(p$mbr)
  expect_true(p$rt_profiling)
  expect_true(p$xic)
  expect_true(p$unimod4)
  expect_true(p$met_excision)

  # max-pr-mz 1800 != default 1200, so it should be in extra
  expect_true(grepl("--max-pr-mz", result$command_line))
  # proteoforms and pg-level are now recognized params
  expect_true(p$proteoforms)
  expect_equal(p$pg_level, 0)
  # ids-to-names is still unrecognized
  expect_true(grepl("--ids-to-names", p$extra_cli_flags))
})
