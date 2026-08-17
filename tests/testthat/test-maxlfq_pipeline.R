# =============================================================================
# build_maxlfq_pipeline() + diann_q_columns()
#
# Written after v4.0.2-4.0.3 shipped with MaxLFQ COMPLETELY BROKEN: the function
# select()ed a column subset that omitted PG.Q.Value / Global.*, then filtered on
# exactly those columns. Filtering an arrow dataset on a select()-dropped column
# via .data[[name]] returns ZERO ROWS instead of erroring, so every modern DIA-NN
# report produced nothing and the user was told to "loosen the QuantUMS cutoffs"
# -- cutoffs that were 0 and never applied.
#
# There was no test on this path at all, which is why it shipped twice. The first
# test below is the one that would have caught it: a perfectly ordinary report in
# which every q-value passes must yield a non-empty matrix.
# =============================================================================

skip_if_no_arrow <- function() {
  skip_if_not_installed("arrow")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("tidyr")
}

# A minimal DIA-NN-shaped report. `extra_q` controls whether the experiment-wide
# columns are present -- the whole point is that BOTH shapes must work.
make_report <- function(path, n_prot = 6, runs = c("runA", "runB", "runC"),
                        extra_q = TRUE, q_override = NULL) {
  grid <- expand.grid(prot = seq_len(n_prot), Run = runs, stringsAsFactors = FALSE)
  df <- data.frame(
    Run            = grid$Run,
    Protein.Group  = sprintf("P%03d", grid$prot),
    Genes          = sprintf("GENE%03d", grid$prot),
    Protein.Names  = sprintf("Prot %d", grid$prot),
    # PG.MaxLFQ is broadcast across a protein's precursor rows within a run --
    # constant per (Protein.Group, Run), which the pipeline asserts.
    PG.MaxLFQ      = 1000 * grid$prot + match(grid$Run, runs),
    Q.Value        = 0.001,
    Lib.Q.Value    = 0.001,
    Lib.PG.Q.Value = 0.001,
    Empirical.Quality = 0.9,
    PG.MaxLFQ.Quality = 0.9,
    stringsAsFactors = FALSE
  )
  if (extra_q) {
    df$PG.Q.Value        <- 0.001
    df$Global.Q.Value    <- 0.001
    df$Global.PG.Q.Value <- 0.001
  }
  if (!is.null(q_override)) df <- q_override(df)
  arrow::write_parquet(df, path)
  path
}

# ---------------------------------------------------------------- regression --

test_that("a report WITH experiment-wide q-columns yields a non-empty matrix", {
  skip_if_no_arrow()
  # THE REGRESSION TEST. Before the fix this stopped with
  # "no precursor rows survived the filters" on exactly this input.
  p <- make_report(tempfile(fileext = ".parquet"), extra_q = TRUE)
  res <- build_maxlfq_pipeline(p, q_cutoff = 0.01)

  expect_true(nrow(res$E) > 0)
  expect_equal(nrow(res$E), 6)
  expect_equal(ncol(res$E), 3)
  expect_false(all(is.na(res$E)))
})

test_that("the extra q-columns are actually APPLIED, not merely tolerated", {
  skip_if_no_arrow()
  # The mirror image of the regression: proving the filter still bites. If a
  # future refactor drops these columns from the filter to "fix" the zero-row
  # bug, this fails. One protein fails Global.PG.Q.Value in every run.
  p <- make_report(tempfile(fileext = ".parquet"), extra_q = TRUE,
                   q_override = function(df) {
                     df$Global.PG.Q.Value[df$Protein.Group == "P001"] <- 0.5
                     df
                   })
  res <- build_maxlfq_pipeline(p, q_cutoff = 0.01)

  expect_equal(nrow(res$E), 5)
  expect_false("P001" %in% rownames(res$E))
  expect_match(paste(res$other$filters_applied, collapse = " "), "Global.PG.Q.Value")
})

test_that("a legacy report WITHOUT the experiment-wide columns still works", {
  skip_if_no_arrow()
  # Older reports lack them. Naming an absent column must not error, and the
  # run-level three must still be applied.
  p <- make_report(tempfile(fileext = ".parquet"), extra_q = FALSE)
  res <- build_maxlfq_pipeline(p, q_cutoff = 0.01)

  expect_equal(nrow(res$E), 6)
  applied <- paste(res$other$filters_applied, collapse = " ")
  expect_match(applied, "Q.Value")
  expect_false(grepl("Global.PG.Q.Value", applied, fixed = TRUE))
})

test_that("run-level q-columns are applied on a report that has them", {
  skip_if_no_arrow()
  p <- make_report(tempfile(fileext = ".parquet"), extra_q = TRUE,
                   q_override = function(df) {
                     df$Q.Value[df$Protein.Group == "P002"] <- 0.5
                     df
                   })
  res <- build_maxlfq_pipeline(p, q_cutoff = 0.01)
  expect_false("P002" %in% rownames(res$E))
})

test_that("filter_counts records a real reduction, not a silent wipe", {
  skip_if_no_arrow()
  # after_fdr == 0 with a non-zero input is the exact signature of the shipped
  # bug, so assert against it directly.
  p <- make_report(tempfile(fileext = ".parquet"), extra_q = TRUE)
  res <- build_maxlfq_pipeline(p, q_cutoff = 0.01)

  expect_true(res$other$filter_counts$input > 0)
  expect_true(res$other$filter_counts$after_fdr > 0)
  expect_equal(res$other$filter_counts$after_fdr, res$other$filter_counts$input)
})

# ------------------------------------------------------- diann_q_columns() ----

test_that("diann_q_columns returns all six when the report has them", {
  skip_if_no_arrow()
  p <- make_report(tempfile(fileext = ".parquet"), extra_q = TRUE)
  expect_setequal(diann_q_columns(p),
                  c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
                    "PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value"))
})

test_that("diann_q_columns returns only what is present", {
  skip_if_no_arrow()
  p <- make_report(tempfile(fileext = ".parquet"), extra_q = FALSE)
  expect_setequal(diann_q_columns(p),
                  c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"))
})

test_that("diann_q_columns falls back to limpa's default when it cannot read", {
  # An unreadable path must not take a data load down -- fall back to exactly the
  # behaviour that existed before the helper.
  expect_equal(diann_q_columns("/nonexistent/nope.parquet"), DIANN_Q_COLUMNS_BASE)
  expect_equal(diann_q_columns(cols = character(0)), DIANN_Q_COLUMNS_BASE)
  expect_equal(diann_q_columns(), DIANN_Q_COLUMNS_BASE)
})

test_that("diann_q_columns never returns a column the report lacks", {
  # Naming an absent column makes limpa::readDIANN() error, so this is the
  # invariant the two load handlers depend on.
  expect_equal(diann_q_columns(cols = c("Run", "Q.Value")), "Q.Value")
  expect_false("Global.Q.Value" %in%
                 diann_q_columns(cols = c("Q.Value", "Lib.Q.Value", "PG.Q.Value")))
})

test_that("diann_q_columns_code emits parseable R naming the same columns", {
  qc <- c("Q.Value", "Global.PG.Q.Value")
  code <- diann_q_columns_code(qc)
  expect_equal(eval(parse(text = code)), qc)
})

# --------------------------------------------------------- degenerate runs ----

test_that("a run with <2 quantified proteins is dropped, not left to crash limma", {
  skip_if_no_arrow()
  # limma::normalizeQuantiles() calls approx() per column and dies with
  # "need at least two non-NA values to interpolate" when a run keeps <2
  # proteins. It must be dropped and named instead.
  p <- make_report(tempfile(fileext = ".parquet"), n_prot = 6,
                   runs = c("runA", "runB", "runC"), extra_q = TRUE,
                   q_override = function(df) {
                     # runC keeps a single protein
                     kill <- df$Run == "runC" & df$Protein.Group != "P001"
                     df$Q.Value[kill] <- 0.5
                     df
                   })
  expect_warning(res <- build_maxlfq_pipeline(p, q_cutoff = 0.01),
                 "dropping 1 run")
  expect_false("runC" %in% colnames(res$E))
  expect_equal(ncol(res$E), 2)
})

# ------------------------------------------------------ the underlying trap ---

test_that("arrow silently returns 0 rows when filtering a dropped column", {
  skip_if_no_arrow()
  # Not testing our code -- pinning the ARROW behaviour that made the bug
  # invisible, so that if a future arrow release turns this into an error (or
  # this assumption stops holding) we find out from the suite rather than from
  # a user with an empty result.
  p <- make_report(tempfile(fileext = ".parquet"), extra_q = TRUE)
  ds <- arrow::open_dataset(p, format = "parquet")
  dropped <- dplyr::select(ds, dplyr::all_of(c("Run", "Protein.Group")))

  n <- tryCatch(
    as.integer(dplyr::collect(dplyr::summarise(
      dplyr::filter(dropped, .data[["Q.Value"]] <= 0.01), n = dplyr::n()))$n),
    error = function(e) NA_integer_)

  # Either arrow errors (NA) or it silently yields 0 -- both mean "never filter
  # on a column you did not select". A non-zero count would mean arrow started
  # reading unselected columns, which would also be worth knowing about.
  expect_true(is.na(n) || n == 0)
})
