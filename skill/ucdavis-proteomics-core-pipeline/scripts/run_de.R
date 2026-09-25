#!/usr/bin/env Rscript
# =============================================================================
# run_de.R  --  Differential expression for the proteomics-pipeline skill.
#
# Two pipelines, selectable with --method. Both are faithful to the recipes
# DE-LIMP logs in R/server_data.R and R/helpers.R:
#
#   dpc     DPC-Quant + limma (limpa). Default DE-LIMP path.
#           readDIANN -> dpcCN -> dpcQuant -> dpcDE -> contrasts.fit -> eBayes.
#           Missing precursors are handled by the detection-probability model;
#           nothing is imputed or dropped. dpcDE wraps voomaLmFitWithImputation,
#           so the imputation uncertainty is propagated into the limma fit.
#           Cite: Li M, Cobbold SA, Smyth GK (2025) bioRxiv 10.1101/2025.04.28.651125;
#                 Li M, Smyth GK (2023) Bioinformatics 39(5):btad200.
#
#   maxlfq  MaxLFQ + limma. DE-LIMP alternative path.
#           DIA-NN PG.MaxLFQ -> log2 -> quantile normalize
#           (limma::normalizeBetweenArrays) -> lmFit -> contrasts.fit -> eBayes.
#           NAs are left in place; limma drops them per row at fit time.
#           Proteins entirely missing in one condition are reported separately
#           as qualitative on/off calls.
#
# Input contract: a DIA-NN report (report.parquet preferred, or report.tsv).
# For Sage / FragPipe DDA results, first convert to a DIA-NN-shaped precursor
# report or a peptide matrix (see references/de-analysis.md, "Feeding non
# DIA-NN engines into the DE step"); then point --input at that.
#
# Usage:
#   Rscript run_de.R --input report.parquet --metadata conditions.csv \
#                    --method dpc --outdir de_results
#
# metadata CSV columns: File.Name,Group[,Batch,Covariate1,Covariate2][,<block column>]
#   File.Name must match the Run / column names in the report.
#
# --block <column>: samples sharing a value of that column (the mouse several IPs came
#   from, the patient several biopsies came from) are fitted as correlated -- a random
#   effect via limma's duplicateCorrelation. See blocking.R.
# --block-scope within|all (default within): which contrasts the blocked fit reports.
#   within: contrasts inside a block (bait vs IgG in the same mice) come from the blocked
#   fit, contrasts BETWEEN blocks (Old vs Young mice) from the same data fitted without it;
#   all: every contrast from the blocked fit. blocking.R says why.
# --block-effect auto|fixed|random (default auto): auto fits the block as a FIXED effect
#   when it is crossed with the groups and every contrast is within one block (a paired
#   design: the exact analysis), and as a random effect otherwise (nested / multi-level).
# =============================================================================

suppressWarnings(suppressMessages({
  ok <- requireNamespace("limma", quietly = TRUE)
}))
if (!ok) stop("limma is required. Install with BiocManager::install('limma').")

# ---- minimal dependency-free arg parser -------------------------------------
args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) return(TRUE)            # bare flag
  val <- args[i + 1]
  if (startsWith(val, "--")) return(TRUE)        # bare flag followed by another flag
  val
}
input     <- getarg("--input")
format    <- getarg("--format", "parquet")
meta_path <- getarg("--metadata")
method    <- getarg("--method", "dpc")
contrasts <- getarg("--contrasts", NULL)
q_cutoff  <- as.numeric(getarg("--q-cutoff", "0.01"))
eq_cutoff <- as.numeric(getarg("--eq-cutoff", "0"))
pgq_cutoff<- as.numeric(getarg("--pgq-cutoff", "0"))
# Fraction of samples a protein must be quantified in to be TESTED (maxlfq path).
# 0.5 matches DE-LIMP's default and the UC Davis Bioinformatics Core recommendation.
cov_min_frac <- as.numeric(getarg("--coverage-min", "0.5"))
outdir    <- getarg("--outdir", "de_results")
# --logfc is a REFERENCE LINE ONLY -- it is drawn and labelled on the volcano and never
# filters anything. Significance is the BH-adjusted p-value alone, which is the exact
# hypothesis eBayes/topTable tested (true logFC != 0). Filtering that list on the
# observed |logFC| afterwards would report a stronger claim ("changed by at least Nx")
# than the error rate covers, since the observed fold change is a point estimate carrying
# no uncertainty. Testing a fold-change threshold properly needs limma::treat()'s interval
# null, not a post-hoc cut (McCarthy & Smyth 2009).
logfc_ref <- as.numeric(getarg("--logfc", "1.0"))
adjp_thr  <- as.numeric(getarg("--adjp", "0.05"))
# Contaminants (Cont_-tagged entries) are REMOVED before quantification by default --
# see contaminants.R for the rule and why. --keep-contaminants quantifies and tests them
# with the sample proteins instead; either way de_provenance.json records which.
keep_contaminants <- isTRUE(getarg("--keep-contaminants", FALSE))
# fetch_fasta.py's sidecar: says whether real proteins sat in the database only as Cont_
# entries, which the filter would then remove too. Same default as the auditors.
fasta_meta <- getarg("--fasta-meta", NULL)
# A conditions.csv column naming the unit samples are nested in (e.g. Mouse). Fitted as a
# random effect -- blocking.R has the model and why. Unset = samples are independent.
block_col <- getarg("--block", NULL)
if (isTRUE(block_col)) stop("--block needs a metadata column name (e.g. --block Mouse)")
block_scope <- getarg("--block-scope", NULL)
if (!is.null(block_scope) && is.null(block_col))
  stop("--block-scope needs --block <column>")
if (is.null(block_scope)) block_scope <- "within"
if (!block_scope %in% c("within", "all"))
  stop("--block-scope must be 'within' (default) or 'all'")
block_effect_req <- getarg("--block-effect", NULL)
if (!is.null(block_effect_req) && is.null(block_col))
  stop("--block-effect needs --block <column>")
if (is.null(block_effect_req)) block_effect_req <- "auto"
if (!block_effect_req %in% c("auto", "fixed", "random"))
  stop("--block-effect must be 'auto' (default), 'fixed' or 'random'")
if (isTRUE(fasta_meta)) stop("--fasta-meta needs a path (<fasta>.meta.json)")
if (is.null(fasta_meta) && file.exists("search.fasta.meta.json"))
  fasta_meta <- "search.fasta.meta.json"
if (!is.null(fasta_meta) && !file.exists(fasta_meta))
  stop("--fasta-meta ", fasta_meta, " does not exist")

if (is.null(input) || is.null(meta_path))
  stop("Required: --input <report> and --metadata <conditions.csv>")

# Where this script lives — used to find its sibling helpers (build_maxlfq.R,
# repro_script.R) whether it is run from the skill directory or a copy of it.
.script_dir <- local({
  f <- grep("^--file=", commandArgs(), value = TRUE)[1]
  if (is.na(f)) getwd() else dirname(normalizePath(sub("^--file=", "", f), mustWork = FALSE))
})
.sibling <- function(nm) {
  for (p in c(file.path(.script_dir, nm), nm)) if (file.exists(p)) return(p)
  NULL
}
# ONE definition of the q-value column set, shared with build_maxlfq.R (and
# mirrored in diann_q_columns.py for the Python scripts). Hand-copied duplicates
# are what let this path drift away from build_maxlfq.R and apply a different
# identification FDR to the same report -- see SKILL_OPEN_DEFECTS #2.
local({
  f <- .sibling("diann_q_columns.R")
  if (is.null(f))
    stop("diann_q_columns.R not found next to run_de.R -- it defines the ",
         "identification-FDR columns and there is no safe default to guess.")
  source(f)
})
# ONE definition of the contaminant filter, shared with build_maxlfq.R.
local({
  f <- .sibling("contaminants.R")
  if (is.null(f))
    stop("contaminants.R not found next to run_de.R -- it defines the contaminant ",
         "filter, and skipping it would leave contaminants in the DE unrecorded.")
  source(f)
})
# ONE definition of the blocking step (--block): validation, record, methods lines.
local({
  f <- .sibling("blocking.R")
  if (is.null(f))
    stop("blocking.R not found next to run_de.R -- it defines the --block model and its record.")
  source(f)
})
if (!method %in% c("dpc", "maxlfq"))
  stop("--method must be 'dpc' or 'maxlfq'")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- read experimental design -----------------------------------------------
meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("File.Name", "Group") %in% names(meta)))
  stop("metadata CSV must have at least File.Name and Group columns")
covariates <- intersect(c("Batch", "Covariate1", "Covariate2"), names(meta))

# ~ 0 + groups [+ covariates]. A function because the design is built twice: from the
# whole metadata up front, so a --block that cannot be fitted fails before quantification
# (which can take an hour), and again from the runs that survived it.
build_design <- function(meta, covariates) {
  groups <- factor(meta$Group)
  ft <- list(groups = groups)
  formula_parts <- c("groups")
  for (cv in covariates) {
    ft[[cv]] <- factor(meta[[cv]])
    formula_parts <- c(formula_parts, cv)
  }
  design <- stats::model.matrix(
    stats::as.formula(paste0("~ 0 + ", paste(formula_parts, collapse = " + "))),
    data = ft)
  colnames(design) <- sub("^groups", "", colnames(design))
  list(design = design, groups = groups, formula_parts = formula_parts)
}
# rank check (DE-LIMP helpers.R guards against this before fitting). A covariate whose
# removal restores the rank and that recurs across groups is usually the animal/subject
# the samples came from, nested in the groups -- a random effect, not a fixed one.
check_design_rank <- function(meta, covariates) {
  full_rank <- function(cv) { d <- build_design(meta, cv)$design; qr(d)$rank == ncol(d) }
  if (full_rank(covariates)) return(invisible(TRUE))
  nested <- Filter(function(cv) full_rank(setdiff(covariates, cv)) &&
                     cv %in% block_candidates(meta, setdiff(covariates, cv)), covariates)
  stop("Design matrix is not full rank — check for confounded covariates / empty groups.",
       if (length(nested)) sprintf(paste0(
         " %s is confounded with the groups but recurs across them: if it names the animal or ",
         "subject the samples came from, fit it as a random effect instead -- rename the column ",
         "(e.g. Mouse) and pass --block Mouse."), nested[1]) else "", call. = FALSE)
}

# The block first: a block that is also a covariate is the more specific error.
if (!is.null(block_col)) {
  block_validate_column(meta, block_col, covariates)
  block_check(trimws(as.character(meta[[block_col]])), build_design(meta, covariates)$design,
              block_col)
  # a requested fixed effect the design cannot carry (block nested in the groups) stops here
  if (identical(block_effect_req, "fixed"))
    block_choose_effect("fixed", trimws(as.character(meta[[block_col]])),
                        build_design(meta, covariates)$design, list(), block_col)
} else {
  .cand <- block_candidates(meta, covariates)
  if (length(.cand))
    message(sprintf(paste0("[run_de] note: metadata column(s) %s recur across groups (the same ",
                           "value in several conditions). If samples sharing a value come from one ",
                           "source (the same animal, patient or lysate), pass --block %s so they are ",
                           "not treated as independent."),
                    paste(.cand, collapse = ", "), .cand[1]))
}
check_design_rank(meta, covariates)

message(sprintf("[run_de] method=%s  q=%.3f  samples=%d  covariates=%s  block=%s  contaminants=%s",
                method, q_cutoff, nrow(meta),
                if (length(covariates)) paste(covariates, collapse = ",") else "none",
                if (is.null(block_col)) "none" else block_col,
                if (keep_contaminants) "kept (--keep-contaminants)" else "removed"))

# ---- build the protein-level object per pipeline ----------------------------
# Both branches produce:  E (proteins x samples, log2), run_names (cols of E),
# genes (annotation df), and a `descriptor` describing the method (provenance).

descriptor <- NULL
quantums_applied <- character(0)   # populated on the dpc path; kept defined for provenance
# Both branches set these: what the contaminant filter found (contaminants.R) and every
# filter the run applied, in order -- recorded in de_provenance.json and methods.txt.
cont_census <- NULL; cont_share <- NULL; cont_col <- NA_character_; cont_intensity <- NA_character_
filters_applied <- character(0)

# limpa/DPC is the DEFAULT path. It needs PRECURSOR-level input: readDIANN() keys on
# Precursor.Id + Precursor.Normalised. DIA-NN's native report.parquet has them; the
# adapters for Sage/FragPipe/Radiant collapse to one protein x run row, so limpa
# cannot run on their output. Check up front and say so in one line, rather than
# letting limpa fail deep inside with a column error the user cannot act on.
if (method == "dpc") {
  have_cols <- tryCatch({
    if (identical(format, "parquet")) names(arrow::open_dataset(input)$schema)
    else names(utils::read.delim(input, nrows = 1, check.names = FALSE))
  }, error = function(e) character(0))
  need_dpc <- c("Precursor.Id", "Precursor.Normalised")
  miss_dpc <- setdiff(need_dpc, have_cols)
  if (length(have_cols) > 0 && length(miss_dpc) > 0)
    stop(sprintf(paste0(
      "--method dpc (limpa) needs PRECURSOR-level input but %s has no %s.\n",
      "  Precursor-level reports that DO work with dpc:\n",
      "    DIA-NN   search_out/report.parquet          (native)\n",
      "    FragPipe dia-quant-output/report.tsv        (--format tsv; its DIA route\n",
      "             bundles DIA-NN, so this IS a DIA-NN report)\n",
      "    Radiant  radiant_to_delimp.py --out <x>.parquet\n",
      "  What does NOT work is the ADAPTED report.parquet from adapt_sage /\n",
      "  adapt_fragpipe / adapt_radiant -- those collapse to one row per protein x\n",
      "  run for the maxlfq path. Point --input at the precursor-level file above,\n",
      "  or use --method maxlfq."),
      basename(input), paste(miss_dpc, collapse = " / ")))

  if (!requireNamespace("limpa", quietly = TRUE))
    stop("limpa is required for --method dpc. BiocManager::install('limpa') (needs R 4.5+, Bioc 3.22+).")
  suppressMessages({ library(limpa); library(limma) })

  # readDIANN() prints "Q-value columns <x> not found." for a missing q-column and
  # then CARRIES ON with that filter simply not applied -- a message, not an error,
  # easy to lose in a long log. Resolving against the real header instead means the
  # run log states exactly which FDR filters were applied, and a report with no usable
  # q-column stops rather than producing unfiltered results that look filtered.
  # (FragPipe's own report.tsv does carry Lib.Q.Value / Lib.PG.Q.Value -- columns 39
  # and 40 -- so it needs no special handling; it works with limpa's defaults.)
  #
  # The columns in q_alt are applied WHENEVER THEY EXIST, not only as a fallback for a missing
  # run-level one. Per the DIA-NN README (master), they are not all the same kind of filter:
  #   Global.Q.Value / Global.PG.Q.Value  "global" -- experiment-wide
  #   PG.Q.Value                          "run-specific q-value for the protein group"
  # Q.Value / Lib.Q.Value / Lib.PG.Q.Value on their own do not control FDR across an
  # experiment: the union of run-level-passing IDs over N runs sits above the nominal cutoff
  # and grows with N. Lib.* is documented as "'global' if the library was created by DIA-NN",
  # but that cannot be relied on -- in some reports both Lib columns are identically 0 and
  # filter nothing -- which is why Global.* is applied explicitly rather than assumed.
  # Measured on identical input (373-run DIA-NN 2.6 report, mouse, 4 cohorts, contaminants
  # removed, Empirical.Quality >= 0.75 -- only the q-column set differs): the three extra
  # columns keep out 227 protein groups that the run-level three admit, 6,531 vs 6,758.
  # Those 227 carry a median of ONE precursor and appear in 9.4% of runs, against 12
  # precursors and 82.3% for the protein groups both sets retain.
  # Which column does the work is not what the name suggests: of the rows dropped,
  # PG.Q.Value accounts for 22,465 (20,447 of them uniquely), Global.PG.Q.Value 15,752 and
  # Global.Q.Value 7,873 -- i.e. the RUN-SPECIFIC column is the largest single contributor.
  # NOTE on cutoffs: DIA-NN recommends PG.Q.Value "at 0.01 to 0.05, typically 0.05 is
  # sufficient". This applies the single --q-cutoff (default 0.01) to it, a deliberate
  # tightening that matches build_maxlfq.R so the two --method values agree on FDR. On this
  # report the choice is nearly free: 0.01 vs 0.05 differs by 2 protein groups (6,564/6,566).
  q_want <- DIANN_FDR_REQUIRED
  q_alt  <- DIANN_FDR_OPTIONAL
  q_use  <- intersect(q_want, have_cols)
  q_extra <- intersect(setdiff(q_alt, q_use), have_cols)
  if (length(q_extra)) {
    q_use <- unique(c(q_use, q_extra))
    message(sprintf("[run_de] additional q-value columns present; filtering on %s",
                    paste(q_use, collapse = " / ")))
  }
  if (length(setdiff(q_want, have_cols)) > 0)
    message(sprintf(paste0("[run_de] this report has no %s; limpa would apply NO filter ",
                           "for those columns without saying so. Filtering on %s instead."),
                    paste(setdiff(q_want, have_cols), collapse = " / "),
                    paste(q_use, collapse = " / ")))
  if (length(q_use) == 0)
    stop("No usable q-value column found in ", basename(input),
         " -- cannot apply identification FDR. Columns present: ",
         paste(have_cols, collapse = ", "))

  # QuantUMS cutoffs on the dpc path.
  #
  # readDIANN() takes no QuantUMS argument, so --eq-cutoff / --pgq-cutoff used to be
  # parsed, silently dropped here, and then written into de_provenance.json anyway --
  # producing UNFILTERED limpa output whose provenance claimed a filter. Silent wrong
  # provenance is worse than a crash. Filter the report first, then hand limpa the
  # filtered file (this is what DE-LIMP's filter_quantums_parquet() does in the app).
  #
  # Note the two scores sit at different levels, and it matters here because limpa
  # RE-QUANTIFIES from precursors:
  #   Empirical.Quality  varies within a (Protein.Group, Run) cell -> precursor-level;
  #                      filtering it removes precursors and limpa re-rolls the protein
  #                      from the survivors. Coherent.
  #   PG.MaxLFQ.Quality  is constant within a cell -> protein x run level; filtering it
  #                      deletes whole cells, which limpa's detection model then imputes
  #                      back. Measured on a 274-run dataset: the fitted DPC slope
  #                      flattened 16% (0.713 -> 0.600), detected fraction halved
  #                      (73.5% -> 37.7%) and significant proteins fell 39%.
  # So pgQ is allowed but warned about; eQ is the sane knob for this path.
  dpc_input <- input
  quantums_applied <- character(0)
  # Does the report carry runs the design does not analyse?
  #
  # If so they MUST be removed before readDIANN(), not after. Subsetting columns
  # afterwards (dat[, .keep]) drops the excluded runs' data but keeps every
  # precursor row they justified -- rows that are then entirely NA in the
  # retained runs. Those rows still reach dpcCN()/dpcQuant(), so the excluded
  # runs go on deciding which proteins are eligible AND shifting the detection
  # model every retained protein is quantified against.
  #
  # Measured on a 12-run report analysed at 7 runs: 158 all-NA precursor rows
  # survived the subset; the protein set gained 3 proteins (exact superset, none
  # lost) -- and 100% of the 4,645 SHARED proteins came out with a different
  # logFC, flipping 54 significance calls. build_maxlfq.R has always applied
  # keep_runs to the arrow query before collect(); dpc now matches it.
  # arrow reads BOTH formats, so the restriction is not parquet-only. A tsv left
  # on the post-hoc dat[, .keep] path would keep exactly the defect this block
  # exists to remove -- and silently, since nothing errors.
  .open_report <- function(path, fmt) {
    if (identical(fmt, "parquet")) arrow::open_dataset(path)
    else arrow::read_tsv_arrow(path, as_data_frame = FALSE)
  }
  .rep_runs_pre <- tryCatch({
    if (requireNamespace("arrow", quietly = TRUE) &&
        requireNamespace("dplyr", quietly = TRUE))
      sort(unique(dplyr::collect(dplyr::select(.open_report(input, format), Run))$Run))
    else NULL
  }, error = function(e) NULL)
  .restrict_runs <- !is.null(.rep_runs_pre) && length(setdiff(.rep_runs_pre, meta$File.Name)) > 0

  if (.restrict_runs || (!is.na(eq_cutoff) && eq_cutoff > 0) ||
      (!is.na(pgq_cutoff) && pgq_cutoff > 0)) {
    if (!identical(format, "parquet") &&
        ((!is.na(eq_cutoff) && eq_cutoff > 0) || (!is.na(pgq_cutoff) && pgq_cutoff > 0)))
      stop("--eq-cutoff / --pgq-cutoff on --method dpc need parquet input.")
    if (!requireNamespace("arrow", quietly = TRUE) || !requireNamespace("dplyr", quietly = TRUE))
      stop("arrow and dplyr are required to pre-filter the report on --method dpc.")
    flt <- .open_report(input, format)
    if (.restrict_runs) {
      .keep_runs <- meta$File.Name
      flt <- dplyr::filter(flt, Run %in% !!.keep_runs)
      message(sprintf(paste0("[run_de] restricting to the %d of %d report runs in the metadata ",
                             "BEFORE quantification (excluded runs must not decide the protein set)"),
                      length(intersect(.rep_runs_pre, meta$File.Name)), length(.rep_runs_pre)))
      quantums_applied <- c(quantums_applied,
                            sprintf("restricted to %d analysed run(s) before rollup",
                                    length(intersect(.rep_runs_pre, meta$File.Name))))
    }
    if (!is.na(eq_cutoff) && eq_cutoff > 0) {
      if (!"Empirical.Quality" %in% have_cols)
        stop("--eq-cutoff given but this report has no Empirical.Quality column.")
      flt <- dplyr::filter(flt, Empirical.Quality >= !!eq_cutoff)
      quantums_applied <- c(quantums_applied, sprintf("Empirical.Quality >= %.2f", eq_cutoff))
    }
    if (!is.na(pgq_cutoff) && pgq_cutoff > 0) {
      if (!"PG.MaxLFQ.Quality" %in% have_cols)
        stop("--pgq-cutoff given but this report has no PG.MaxLFQ.Quality column.")
      warning("--pgq-cutoff on --method dpc deletes whole protein x run cells, which limpa's ",
              "detection model then imputes back. Prefer --eq-cutoff for this path.")
      flt <- dplyr::filter(flt, PG.MaxLFQ.Quality >= !!pgq_cutoff)
      quantums_applied <- c(quantums_applied, sprintf("PG.MaxLFQ.Quality >= %.2f", pgq_cutoff))
    }
    # Write back in the SAME format limpa is about to be told to read. Writing a
    # parquet and still passing format = "tsv" would fail inside readDIANN().
    if (identical(format, "parquet")) {
      dpc_input <- file.path(tempdir(), "quantums_filtered_for_dpc.parquet")
      arrow::write_parquet(dplyr::collect(flt), dpc_input)
    } else {
      dpc_input <- file.path(tempdir(), "runs_filtered_for_dpc.tsv")
      # base write.table, not arrow::write_csv_arrow -- the latter does not accept
      # a `delim` argument in arrow 24 ("not yet supported in Arrow"), so it
      # cannot write a TSV at all. collect() has already materialised the frame.
      utils::write.table(dplyr::collect(flt), dpc_input, sep = "\t",
                         quote = FALSE, row.names = FALSE, na = "")
    }
    message("[run_de] pre-filter for dpc: ", paste(quantums_applied, collapse = " | "))
  }

  # Per-column cutoffs. limpa recycles q.cutoffs against q.columns
  # (EListFromLongFormatFile: rep_len(q.cutoffs, length(q.columns)), then
  # Report[[q.columns[j]]] > q.cutoffs[j]), so a vector is honoured element-wise.
  # PG.Q.Value takes DIA-NN's own recommended 0.05 rather than the uniform
  # --q-cutoff -- see diann_q_columns.R.
  q_cuts <- vapply(q_use, diann_cutoff_for, numeric(1), q_cutoff = q_cutoff)
  # limpa's own annotation columns, plus the accession column(s) the contaminant filter
  # reads. Protein.Ids is dropped again after filtering, so nothing downstream changes.
  dpc_ann_default <- eval(formals(limpa::readDIANN)$annotation.columns)
  dpc_ann <- unique(c(dpc_ann_default, CONTAMINANT_ID_COLUMNS))
  dat <- limpa::readDIANN(dpc_input, format = format, q.cutoffs = unname(q_cuts),
                          q.columns = q_use, annotation.columns = dpc_ann)
  message(sprintf("[run_de] readDIANN: %d precursors x %d runs (FDR on %s)",
                  nrow(dat$E), ncol(dat$E),
                  paste(sprintf("%s<=%.3f", q_use, q_cuts), collapse = ", ")))

  # ---- FIX: reconcile against the metadata BEFORE quantifying -----------------
  # dpcQuant() on a 274-run report takes ~90 min. The run/metadata check used to sit
  # after it, so a metadata file covering a SUBSET of the report burned the full
  # quantification and then died with "Some report columns have no metadata row".
  # maxlfq already honours keep_runs; dpc now does too.
  .rep_runs <- colnames(dat$E)
  .keep     <- .rep_runs %in% meta$File.Name
  if (!any(.keep))
    stop("No report run matches File.Name in the metadata. First report run: ", .rep_runs[1])
  if (any(!.keep)) {
    message(sprintf("[run_de] restricting to the %d of %d report runs present in the metadata",
                    sum(.keep), length(.keep)))
    dat <- dat[, .keep]
  }
  .absent <- setdiff(meta$File.Name, .rep_runs)
  if (length(.absent))
    message(sprintf("[run_de] note: %d metadata row(s) have no run in the report, e.g. %s",
                    length(.absent), paste(utils::head(.absent, 2), collapse = ", ")))

  # ---- contaminants: out BEFORE dpcCN, or they shape the detection model and the DE ----
  # Rows of `dat` are precursors that passed the identification filters, so the counts
  # are exactly what limpa would have quantified. Rows left all-NA by the run subset
  # above carry no data in the analysed runs and are not counted.
  cont_col <- contaminant_id_column(names(dat$genes))
  if (!is.na(cont_col)) {
    .flag <- is_contaminant(dat$genes[[cont_col]])
    .seen <- rowSums(!is.na(dat$E)) > 0
    cont_census <- contaminant_census(group = dat$genes$Protein.Group[.seen],
                                      is_cont = .flag[.seen], feature = rownames(dat$E)[.seen],
                                      genes = dat$genes$Genes[.seen])
    # The share is on the measured signal (Precursor.Quantity), read for exactly the
    # precursor x run cells readDIANN kept; Precursor.Normalised -- what dat$E holds --
    # only when the report has no Precursor.Quantity.
    .qty <- tryCatch({
      .cols <- c("Run", "Precursor.Id", CONTAMINANT_SHARE_COLUMNS[1])
      .r <- if (identical(format, "parquet")) nanoparquet::read_parquet(dpc_input, col_select = .cols)
            else data.table::fread(dpc_input, select = .cols, data.table = FALSE, showProgress = FALSE)
      .i <- match(.r$Precursor.Id, rownames(dat$E)); .j <- match(.r$Run, colnames(dat$E))
      .ok <- !is.na(.i) & !is.na(.j)
      .ok[.ok] <- !is.na(dat$E[cbind(.i[.ok], .j[.ok])])
      .m <- matrix(NA_real_, nrow(dat$E), ncol(dat$E), dimnames = dimnames(dat$E))
      .m[cbind(.i[.ok], .j[.ok])] <- .r[[CONTAMINANT_SHARE_COLUMNS[1]]][.ok]
      .m
    }, error = function(e) {
      message("[run_de] contaminant share: ", CONTAMINANT_SHARE_COLUMNS[1], " not readable (",
              conditionMessage(e), "); using ", CONTAMINANT_SHARE_COLUMNS[2])
      NULL
    })
    cont_intensity <- if (is.null(.qty)) CONTAMINANT_SHARE_COLUMNS[2] else CONTAMINANT_SHARE_COLUMNS[1]
    cont_share <- contaminant_share(if (is.null(.qty)) 2^dat$E else .qty, .flag)
    if (!keep_contaminants && any(.flag)) dat <- dat[!.flag, ]
  }
  .extra <- setdiff(names(dat$genes), dpc_ann_default)
  if (length(.extra)) dat$genes <- dat$genes[, setdiff(names(dat$genes), .extra), drop = FALSE]

  dpcfit    <- limpa::dpcCN(dat)
  y_protein <- limpa::dpcQuant(dat, "Protein.Group", dpc = dpcfit)
  E         <- y_protein$E
  run_names <- colnames(E)
  genes     <- y_protein$genes
  descriptor <- list(
    pipeline_id   = "dpc",
    display_label = "DPC-Quant + limma (limpa)",
    rollup_method = "DPC-Quant (Detection Probability Curve quantification, dpcCN)",
    quantums_filter = if (length(quantums_applied)) paste(quantums_applied, collapse = " | ") else "none",
    de_engine     = "limpa::dpcDE (voomaLmFitWithImputation) -> contrasts.fit -> eBayes",
    missing_policy = "Missing precursors modelled via the detection probability curve; not imputed, not dropped.",
    citation      = "Li M, Cobbold SA, Smyth GK (2025) bioRxiv 10.1101/2025.04.28.651125; Li M, Smyth GK (2023) Bioinformatics 39(5):btad200"
  )
  filters_applied <- c(sprintf("%s <= %g", q_use, q_cuts), quantums_applied)

} else { # maxlfq
  if (!requireNamespace("arrow", quietly = TRUE) && identical(format, "parquet"))
    stop("arrow is required to read parquet for --method maxlfq. install.packages('arrow').")
  suppressMessages(library(limma))
  bm <- .sibling("build_maxlfq.R")             # defines build_maxlfq()
  if (is.null(bm)) stop("build_maxlfq.R not found next to run_de.R")
  source(bm)

  keep_runs <- meta$File.Name
  # Record the per-column cutoffs build_maxlfq() will apply, so the emitted
  # reproducibility script can state them instead of the scalar --q-cutoff.
  q_use  <- diann_fdr_columns(names(arrow::open_dataset(input, format = format)$schema))
  q_cuts <- vapply(q_use, diann_cutoff_for, numeric(1), q_cutoff = q_cutoff)
  ml <- build_maxlfq(input, format = format, q_cutoff = q_cutoff,
                     eq_cutoff = eq_cutoff, pgq_cutoff = pgq_cutoff,
                     keep_runs = keep_runs, drop_contaminants = !keep_contaminants)
  E         <- ml$E
  run_names <- colnames(E)
  genes     <- ml$genes
  descriptor <- ml$descriptor
  q_use     <- ml$q_columns
  cont_census    <- ml$contaminants$census
  cont_share     <- ml$contaminants$share
  cont_col       <- ml$contaminants$id_column
  cont_intensity <- ml$contaminants$intensity_column
  message(sprintf("[run_de] MaxLFQ matrix: %d proteins x %d runs (%.1f%% missing)",
                  nrow(E), ncol(E), 100 * mean(is.na(E))))

  # ---- coverage filter (ported from DE-LIMP server_data.R) -----------------
  # The MaxLFQ matrix keeps every protein seen in ANY run, so rows with 1-2
  # finite values reach limma. eBayes then moderates variance against rows whose
  # variance is barely estimable, which destabilises the whole fit -- not just
  # those rows. Require a protein to be quantified in at least `--coverage-min`
  # of samples before testing; the rest are an on/off observation, not a
  # differential-abundance result.
  n_samples     <- ncol(E)
  min_obs       <- max(2, ceiling(cov_min_frac * n_samples))
  n_obs_per_row <- rowSums(!is.na(E))
  keep_cov      <- n_obs_per_row >= min_obs
  n_dropped     <- sum(!keep_cov)
  message(sprintf("[run_de] coverage filter: keep proteins with >= %d/%d non-NA (%.0f%%). Kept %d, dropped %d.",
                  min_obs, n_samples, 100 * cov_min_frac, sum(keep_cov), n_dropped))
  if (sum(keep_cov) < 10)
    stop(sprintf(paste0("Coverage filter left only %d testable proteins (needed >= %d non-NA ",
                        "of %d samples). Loosen --coverage-min, or the QuantUMS cutoffs ",
                        "(--eq-cutoff / --pgq-cutoff) if those are doing the damage."),
                 sum(keep_cov), min_obs, n_samples))
  E <- E[keep_cov, , drop = FALSE]
  if (!is.null(genes) && nrow(genes) == length(keep_cov))
    genes <- genes[keep_cov, , drop = FALSE]
  # build_maxlfq() returns its filters beside the descriptor, not inside it -- reading
  # descriptor$filters_applied here used to drop every one of them from the record.
  filters_applied <- c(ml$filters_applied,
    sprintf("coverage >= %d/%d non-NA samples (%.0f%%); %d proteins dropped",
            min_obs, n_samples, 100 * cov_min_frac, n_dropped))
  descriptor$filters_applied <- filters_applied
}

# ---- the contaminant record (contaminants.R): one description, read by everything ----
cont_risk <- if (!is.null(cont_census) && cont_census$n_precursors > 0 && !keep_contaminants)
  contaminant_database_risk(fasta_meta, cont_census$n_groups_contaminant) else NULL
cont_rec <- contaminant_record(cont_census, cont_share, keep_contaminants, cont_col,
                               cont_intensity, risk = cont_risk,
                               fasta_meta = if (is.null(fasta_meta)) NULL
                                            else normalizePath(fasta_meta, mustWork = FALSE))
if (method == "dpc")   # build_maxlfq() records its own contaminant step in ml$filters_applied
  filters_applied <- c(filters_applied, switch(cont_rec$policy,
    removed = sprintf("contaminants removed: %d precursors mapping to a %s entry (%s); %d %s protein groups",
                      cont_rec$n_precursors, cont_rec$tag, cont_rec$id_column,
                      cont_rec$n_protein_groups, cont_rec$tag),
    kept = sprintf("contaminants kept (--keep-contaminants): %d %s protein groups",
                   cont_rec$n_protein_groups, cont_rec$tag),
    NULL))
for (.l in contaminant_methods_lines(cont_rec)) message("[run_de] ", sub("^ +", "", .l))
if (isTRUE(cont_rec$database_risk))
  warning("contaminant filter: ", cont_rec$database_note, call. = FALSE)
if (!is.null(cont_share)) {
  .qs <- cont_share
  .qs$Group <- meta$Group[match(.qs$Run, meta$File.Name)]
  .qs$Intensity.Column <- cont_intensity
  utils::write.csv(.qs[, c("Run", "Group", "Contaminant.Pct", "Contaminant.Intensity",
                           "Total.Intensity", "Intensity.Column")],
                   file.path(outdir, cont_rec$share_table), row.names = FALSE)
}
if (isTRUE(cont_rec$removed))
  utils::write.csv(cont_census$groups, file.path(outdir, cont_rec$removed_table), row.names = FALSE)

# ---- align metadata to the matrix columns -----------------------------------
meta <- meta[match(run_names, meta$File.Name), , drop = FALSE]
if (any(is.na(meta$Group)))
  stop("Some report columns have no metadata row. Report runs:\n  ",
       paste(run_names, collapse = "\n  "))

# ---- design matrix (~ 0 + groups [+ covariates]) ----------------------------
.dz <- build_design(meta, covariates)
design <- .dz$design; groups <- .dz$groups; formula_parts <- .dz$formula_parts

check_design_rank(meta, covariates)   # again, on the runs that survived quantification

# ---- blocking factor (--block): re-checked on the runs actually analysed ------
block <- NULL
if (!is.null(block_col)) {
  block <- trimws(as.character(meta[[block_col]]))
  block_check(block, design, block_col)
  meta[[block_col]] <- block   # the values fitted are the values the repro script and session carry
}

# ---- contrasts --------------------------------------------------------------
lvls <- levels(groups)
if (is.null(contrasts)) {
  ref  <- lvls[1]
  forms <- paste0(lvls[-1], "-", ref)        # all groups vs first level
} else {
  forms <- trimws(strsplit(contrasts, ",")[[1]])
}
message("[run_de] contrasts: ", paste(forms, collapse = ", "))
cmat <- limma::makeContrasts(contrasts = forms, levels = design)

# ---- fit --------------------------------------------------------------------
# The fit with samples as independent. With --block it is still needed for the contrasts
# --block-scope within reports from it (between-block ones); it reuses the quantification.
fit_independent <- function()
  if (method == "dpc") limpa::dpcDE(y_protein, design, plot = FALSE) else limma::lmFit(E, design)
block_rec <- block_none_record()
fit_ind <- NULL
# Fixed or random block (blocking.R: block_choose_effect). A fixed block is extra design
# columns, so the design and contrast matrix the fit uses change with it.
.eff <- if (!is.null(block))
  block_choose_effect(block_effect_req, block, design,
                      block_contrast_structure(block, groups, cmat), block_col) else NULL
if (is.null(block)) {
  fit_ind <- fit_independent()
} else if (identical(.eff$effect, "fixed")) {
  design <- .eff$design
  cmat <- limma::makeContrasts(contrasts = forms, levels = design)
  fit <- fit_independent()          # samples independent GIVEN the block columns
  block_rec <- block_record(block_col, block, method, NA_real_, numeric(0), nrow(E),
                            groups = groups, cmat = cmat, scope = block_scope,
                            effect = "fixed", effect_choice = .eff$choice)
} else if (method == "dpc") {
  # dpcDE(y, design, plot, ...) hands `block` to voomaLmFitWithImputation(), which
  # estimates the correlation with the vooma weights twice and fits lmFit(block =,
  # correlation =, weights =) -- see blocking.R. Its two estimates arrive as messages.
  .pass <- c(first = NA_real_, final = NA_real_)
  fit <- withCallingHandlers(
    limpa::dpcDE(y_protein, design, plot = FALSE, block = block),
    message = function(m) {
      .m <- regmatches(conditionMessage(m),
                       regexec("^(First|Final) intra-block correlation +(\\S+)", conditionMessage(m)))[[1]]
      if (length(.m) == 3) .pass[[tolower(.m[2])]] <<- as.numeric(.m[3])
    })
  # The per-protein estimates behind the consensus, for the record: the same call limpa's
  # final pass made (fit$EList holds its y and final weights), so it must agree with it.
  .dc <- suppressWarnings(limma::duplicateCorrelation(fit$EList, design, block = block,
                                                      weights = fit$EList$weights))
  .rho <- if (!is.null(fit$correlation)) fit$correlation else .pass[["final"]]
  if (!isTRUE(all.equal(.dc$consensus.correlation, .rho, tolerance = 1e-6)))
    warning(sprintf("--block: recomputed consensus correlation %.4f differs from limpa's %.4f",
                    .dc$consensus.correlation, .rho), call. = FALSE)
  block_rec <- block_record(block_col, block, method, .rho, .dc$atanh.correlations, nrow(E),
                            first_pass = .pass[["first"]], groups = groups, cmat = cmat,
                            scope = block_scope, effect_choice = .eff$choice)
} else {
  .dc <- limma::duplicateCorrelation(E, design, block = block)
  if (!is.finite(.dc$consensus.correlation))
    stop(sprintf(paste0("--block %s: no protein had enough quantified samples to estimate the ",
                        "within-block correlation. Drop --block, or loosen --coverage-min only if ",
                        "that is what emptied the matrix."), block_col))
  fit <- limma::lmFit(E, design, block = block, correlation = .dc$consensus.correlation)
  block_rec <- block_record(block_col, block, method, .dc$consensus.correlation,
                            .dc$atanh.correlations, nrow(E), groups = groups, cmat = cmat,
                            scope = block_scope, effect_choice = .eff$choice)
}
# Which fit reports each contrast: block.contrast_model, the one definition.
contrast_model <- if (is.null(block)) {
  setNames(rep("independent", length(forms)), forms)
} else unlist(block_rec$contrast_model)[forms]
if (!is.null(block)) {
  .fit_blocked <- limma::eBayes(limma::contrasts.fit(fit, cmat))
  if (any(contrast_model == "independent"))
    fit_ind <- fit_independent()
  descriptor$de_engine <- paste0(descriptor$de_engine, block_engine_suffix(block_rec))
  message(if (identical(block_rec$effect, "fixed"))
    sprintf("[run_de] --block %s: %d levels, FIXED effect (%s)", block_col, block_rec$n_blocks,
            block_rec$effect_choice)
  else sprintf("[run_de] --block %s (scope %s): %d levels; consensus within-%s correlation %.3f (%d of %d proteins estimated)",
               block_col, block_scope, block_rec$n_blocks, block_col,
               block_rec$consensus_correlation, block_rec$n_proteins_estimated,
               block_rec$n_proteins))
  for (.w in unlist(block_rec$warnings)) warning("--block: ", .w, call. = FALSE)
}
if (!is.null(fit_ind)) fit_ind <- limma::eBayes(limma::contrasts.fit(fit_ind, cmat))
# One fit object whose columns are each contrast's reporting fit (block_merge_fits), so the
# DE tables and the DE-LIMP session read the same numbers.
fit <- if (is.null(block)) fit_ind else block_merge_fits(.fit_blocked, fit_ind, contrast_model)

# The fitted design as text: a fixed block is a design term, a random one is not.
design_label <- paste0("~ 0 + ", paste(c(formula_parts,
                       if (identical(block_rec$effect, "fixed")) block_col), collapse = " + "))

# ---- write per-contrast results ---------------------------------------------
gene_cols <- intersect(c("Genes", "Protein.Names"), names(genes))
# limpa returns the protein identifier in rownames(genes), not as a column. Guarding
# the merge key on its presence silently dropped it, so the merges below failed with
# "'by' must specify a uniquely valid column" AFTER quantification had completed.
ann <- if (length(gene_cols)) {
  .a <- genes[, gene_cols, drop = FALSE]
  .a$Protein.Group <- if ("Protein.Group" %in% names(genes)) as.character(genes$Protein.Group) else rownames(genes)
  .a[, c("Protein.Group", gene_cols), drop = FALSE]
} else NULL

# Expression matrix (proteins x samples, log2) — feeds figures (PCA/heatmap) and
# the report; mirrors DE-LIMP's Expression_Matrix.csv export.
expr_df <- data.frame(Protein.Group = rownames(E), check.names = FALSE)
if (!is.null(ann)) expr_df <- merge(expr_df, ann, by = "Protein.Group", all.x = TRUE, sort = FALSE)
expr_df <- merge(expr_df, data.frame(Protein.Group = rownames(E), E, check.names = FALSE),
                 by = "Protein.Group", all.x = TRUE, sort = FALSE)
utils::write.csv(expr_df, file.path(outdir, "Expression_Matrix.csv"), row.names = FALSE)
message(sprintf("[run_de] Expression_Matrix.csv: %d proteins x %d samples", nrow(E), ncol(E)))

# ---- detection matrix: measured vs inferred/missing, per protein x sample ------
# The ONE definition of "was this value measured": QC_detected_vs_inferred.csv's totals
# below are column sums of this matrix, and the figures read it cell by cell (DE-LIMP
# colours each violin dot by it). dpc: precursors actually observed for that protein in
# that run (0 = value inferred by the detection-probability model). maxlfq: 1 where the
# matrix has a value, 0 where it is NA (missing). Rows and sample columns follow
# Expression_Matrix.csv exactly.
det_n <- if (method == "dpc" && exists("dat")) {
  .prot <- as.character(dat$genes[["Protein.Group"]])
  if (is.null(.prot) || !length(.prot)) .prot <- rownames(dat$E)
  .cnt <- rowsum((!is.na(dat$E)) * 1L, group = .prot, reorder = TRUE)
  .m <- matrix(NA_integer_, nrow(E), ncol(E), dimnames = dimnames(E))
  .ii <- match(rownames(E), rownames(.cnt))
  .m[!is.na(.ii), ] <- as.matrix(.cnt)[.ii[!is.na(.ii)], colnames(E), drop = FALSE]
  .m
} else {
  .m <- (!is.na(E)) * 1L
  dimnames(.m) <- dimnames(E)
  .m
}
detection_rec <- list(
  file = "Detection_Matrix.csv",
  zero_means = if (method == "dpc") "inferred" else "missing",
  values = if (method == "dpc")
    "precursors observed for the protein in that run; 0 = value inferred by the DPC model"
  else "1 = quantified in the MaxLFQ matrix, 0 = missing (NA)",
  na_means = "protein not in the precursor matrix: status not recorded")
utils::write.csv(data.frame(Protein.Group = expr_df$Protein.Group,
                            det_n[match(expr_df$Protein.Group, rownames(det_n)), , drop = FALSE],
                            check.names = FALSE),
                 file.path(outdir, detection_rec$file), row.names = FALSE)

all_sig <- list()
de_tables <- list()   # per DE table: its file and which fit it was reported from
for (cn in forms) {
  tt <- limma::topTable(fit, coef = cn, number = Inf, adjust.method = "BH")
  tt$Protein.Group <- rownames(tt)
  # topTable() already carries fit$genes, so merging ann in wholesale produced
  # Genes.x/Genes.y duplicates. Only bring across columns that aren't there yet.
  if (!is.null(ann)) {
    .miss <- setdiff(names(ann), c("Protein.Group", names(tt)))
    if (length(.miss)) tt <- merge(tt, ann[, c("Protein.Group", .miss), drop = FALSE],
                                   by = "Protein.Group", all.x = TRUE, sort = FALSE)
  }
  tt <- tt[order(tt$adj.P.Val), ]
  fn <- file.path(outdir, sprintf("DE_%s_%s.csv", method, make.names(cn)))
  utils::write.csv(tt, fn, row.names = FALSE)
  sig <- subset(tt, !is.na(adj.P.Val) & adj.P.Val < adjp_thr)
  all_sig[[cn]] <- nrow(sig)
  n_beyond <- sum(abs(sig$logFC) >= logfc_ref, na.rm = TRUE)   # descriptive, not a filter
  message(sprintf("[run_de] %-20s  %d proteins, %d significant (adj.P<%.2g); %d of those with |logFC|>=%.2g -> %s%s",
                  cn, nrow(tt), nrow(sig), adjp_thr, n_beyond, logfc_ref, basename(fn),
                  if (is.null(block)) "" else sprintf("  [%s fit]", contrast_model[[cn]])))
  de_tables[[cn]] <- list(file = basename(fn), model = contrast_model[[cn]], n_significant = nrow(sig))
}

# ---- the analysis as plain R ------------------------------------------------
# The bundle in reproducibility/ pins the whole environment, which is what you need
# to reproduce byte-for-byte -- but it is not something a reviewer can read. This is:
# flat, literal R with every value written out, runnable with nothing but R and the
# packages it loads. Same idea as DE-LIMP's reproducibility log. Generated from the
# objects that actually ran, so it cannot describe a different analysis than the one
# in the CSVs next to it.
repro_lines <- NULL
.rs <- .sibling("repro_script.R")
if (is.null(.rs)) {
  message("[run_de] repro_script.R not found next to run_de.R — reproducibility_log.R not written")
} else {
  tryCatch({
    source(.rs)
    repro_lines <- write_repro_script(
      path = file.path(outdir, "reproducibility_log.R"),
      method = method, input = input, format = format,
      # q_use / q_cuts exist only on the dpc branch; guard rather than rely on
      # the surrounding tryCatch, which would swallow the error and silently
      # skip writing reproducibility_log.R altogether.
      q_cutoff = q_cutoff,
      q_columns = if (exists("q_use")) q_use else NULL,
      q_cutoffs = if (exists("q_cuts")) unname(q_cuts) else NULL,
      eq_cutoff = eq_cutoff, pgq_cutoff = pgq_cutoff,
      cov_min_frac = cov_min_frac,
      meta = meta, covariates = covariates, formula_parts = formula_parts,
      forms = forms, adjp_thr = adjp_thr, logfc_ref = logfc_ref,
      ann_cols = gene_cols, descriptor = descriptor,
      contaminants = cont_rec, block = block_rec, contrast_model = contrast_model,
      dpc_annotation_columns = if (exists("dpc_ann")) dpc_ann else NULL)
    message(sprintf("[run_de] reproducibility_log.R: the analysis as %d lines of plain R (Rscript-runnable)",
                    length(repro_lines)))
  }, error = function(e) message("[run_de] reproducibility_log.R not written: ", e$message))
}

# ---- methods + reproducibility provenance -----------------------------------
methods_txt <- c(
  "Differential expression — methods",
  strrep("=", 40), "",
  sprintf("Pipeline      : %s", descriptor$display_label),
  sprintf("Quantification: %s", descriptor$rollup_method),
  sprintf("DE engine     : %s", descriptor$de_engine),
  sprintf("Missing values: %s", descriptor$missing_policy),
  sprintf("ID FDR cutoff : q <= %.3f", q_cutoff),
  sprintf("Filters       : %s", if (length(filters_applied)) filters_applied[1] else "none"),
  if (length(filters_applied) > 1) sprintf("                %s", filters_applied[-1]),
  contaminant_methods_lines(cont_rec),
  if (method == "maxlfq") sprintf("Normalization : quantile (limma::normalizeBetweenArrays)") else
                          sprintf("Normalization : DPC-CN (applied within dpcCN before dpcQuant)"),
  sprintf("Design        : %s", design_label),
  block_methods_lines(block_rec),
  sprintf("Contrasts     : %s", paste(forms, collapse = ", ")),
  sprintf("Significance  : adj.P.Val < %.3g (Benjamini-Hochberg), moderated t-test of", adjp_thr),
  "                H0: log2 fold change = 0. No fold-change filter is applied --",
  sprintf("                |log2FC| = %.3g is drawn on the volcano for reference only, so", logfc_ref),
  "                the reported error rate matches the hypothesis actually tested.",
  "",
  sprintf("Citation      : %s", descriptor$citation)
)
writeLines(methods_txt, file.path(outdir, "methods.txt"))

# ---- exact R provenance: sessionInfo + machine-readable record --------------
# Faithful capture of the stack that actually ran this DE (not re-derived later).
si <- file.path(outdir, "sessionInfo.txt")
con <- file(si, "w"); sink(con)
cat(R.version.string, "\n\n")
for (p in c("limpa", "limma", "arrow", "dplyr", "tidyr"))
  try(cat(sprintf("%-8s %s\n", p, as.character(packageVersion(p)))), silent = TRUE)
cat("\n"); print(utils::sessionInfo())
sink(); close(con)

# Minimal JSON writer — use jsonlite if present, else a dependency-free fallback.
jsonlite_or_manual <- function(x) {
  if (requireNamespace("jsonlite", quietly = TRUE))
    return(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"))
  esc <- function(s) gsub('"', '\\\\"', gsub('\\\\', '\\\\\\\\', as.character(s)))
  enc <- function(v) {
    if (is.null(v) || length(v) == 0) return("null")
    if (is.list(v)) {
      nm <- names(v)
      if (is.null(nm)) return(paste0("[", paste(vapply(v, enc, ""), collapse = ", "), "]"))
      return(paste0("{", paste(sprintf('"%s": %s', esc(nm), vapply(v, enc, "")), collapse = ", "), "}"))
    }
    if (length(v) > 1) return(paste0("[", paste(vapply(v, enc, ""), collapse = ", "), "]"))
    if (is.numeric(v) || is.logical(v)) return(tolower(as.character(v)))
    if (is.na(v)) return("null")
    sprintf('"%s"', esc(v))
  }
  enc(x)
}
pkg_ver <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
prov <- list(
  pipeline_id = descriptor$pipeline_id, display_label = descriptor$display_label,
  rollup_method = descriptor$rollup_method, de_engine = descriptor$de_engine,
  missing_policy = descriptor$missing_policy, citation = descriptor$citation,
  method = method, q_cutoff = q_cutoff,
  # WHICH columns were filtered and at WHAT cutoff. q_cutoff alone stopped being
  # the whole story once PG.Q.Value took DIA-NN's recommended 0.05: a record
  # saying "q_cutoff: 0.01" while six different cutoffs ran is the kind of
  # plausible-but-wrong provenance rule #2 exists to prevent.
  q_columns = if (exists("q_use")) as.list(q_use) else NULL,
  q_cutoffs = if (exists("q_cuts")) as.list(unname(q_cuts)) else NULL,
  # Ties this DE run to the EXACT search output it consumed. A path alone is not
  # traceability -- it does not survive the file being moved, regenerated, or
  # shipped in a bundle without its input.
  input_sha256 = tryCatch(
    if (file.exists(input) && requireNamespace("digest", quietly = TRUE))
      digest::digest(input, algo = "sha256", file = TRUE)
    else if (file.exists(input) && nzchar(Sys.which("shasum")))
      sub(" .*", "", system2("shasum", c("-a", "256", shQuote(input)), stdout = TRUE)[1])
    else NA_character_,
    error = function(e) NA_character_),
  input_bytes = tryCatch(as.numeric(file.info(input)$size), error = function(e) NA_real_),
  # Report the cutoffs as APPLIED. On the dpc path these were previously recorded
  # whether or not they had any effect, so a provenance file could assert a filter the
  # run never performed.
  eq_cutoff  = if (method == "dpc" && !any(grepl("Empirical", quantums_applied))) 0 else eq_cutoff,
  pgq_cutoff = if (method == "dpc" && !any(grepl("PG.MaxLFQ", quantums_applied))) 0 else pgq_cutoff,
  # Every filter the run applied, in order, and what the contaminant filter did. Both are
  # the record make_methods.py writes the Methods from -- it never assumes either.
  filters_applied = as.list(filters_applied),
  contaminants = cont_rec,
  logfc = logfc_ref, logfc_role = "reference_line_only", adjp = adjp_thr,
  significance_rule = "adj.P.Val < adjp (BH); no fold-change filter",
  design = design_label,
  # The random blocking factor (--block), if any: column, consensus correlation, how it
  # was estimated and fitted, and which contrasts are within vs between blocks.
  block = block_rec,
  # as.list: one contrast must still serialise as a list -- auto_unbox made it a bare
  # string, which readers then iterated character by character ("B, -, A").
  contrasts = as.list(forms), n_samples = nrow(meta), groups = as.list(table(groups)),
  significant_per_contrast = all_sig,
  # Each DE_*.csv and the fit it came from: "blocked", or "independent" (samples fitted as
  # independent -- every table of an unblocked run, between-block ones under --block-scope within).
  de_tables = de_tables,
  # Detection_Matrix.csv: what each 0 means depends on the pipeline -- say it here.
  detection_matrix = detection_rec,
  R_version = as.character(getRversion()),
  packages = list(limpa = pkg_ver("limpa"), limma = pkg_ver("limma"),
                  arrow = pkg_ver("arrow"), dplyr = pkg_ver("dplyr"), tidyr = pkg_ver("tidyr")),
  input = normalizePath(input, mustWork = FALSE), metadata = normalizePath(meta_path, mustWork = FALSE)
)
# The block's column name at the top level too, for readers that only need "was the design
# blocked, and on what" (the submission workflow). Present ONLY when blocked; `block` is the
# full record either way (applied = FALSE says independence was modelled).
if (isTRUE(block_rec$applied))
  prov <- append(prov, list(block_column = block_rec$column), after = which(names(prov) == "block"))
writeLines(jsonlite_or_manual(prov), file.path(outdir, "de_provenance.json"))

# ---- detected vs inferred QC ------------------------------------------------
# DPC models missing precursors rather than leaving holes, so the protein matrix
# is complete BY CONSTRUCTION. Counting non-missing cells therefore reports the
# same number for every sample and hides real depth differences. The meaningful
# split is detected (>=1 precursor actually observed in that run) vs inferred
# (value supplied by the detection-probability model) — the same view DE-LIMP's
# Data Completeness panel shows.
qc_di <- NULL
if (method == "dpc" && exists("dat") && exists("y_protein")) {
  tryCatch({
    det <- !is.na(det_n) & det_n > 0          # from Detection_Matrix.csv -- one definition

    qc_di <- data.frame(
      Sample   = colnames(E),
      Group    = if (exists("groups")) as.character(groups) else NA_character_,
      Detected = colSums(det),
      Inferred = nrow(E) - colSums(det),
      Total    = nrow(E),
      stringsAsFactors = FALSE
    )
    qc_di$PctDetected <- round(100 * qc_di$Detected / qc_di$Total, 1)
    qc_di$PctInferred <- round(100 * qc_di$Inferred / qc_di$Total, 1)
    qc_di <- qc_di[order(-qc_di$Detected), ]
    utils::write.csv(qc_di, file.path(outdir, "QC_detected_vs_inferred.csv"), row.names = FALSE)
    message(sprintf("[run_de] QC_detected_vs_inferred.csv: %.0f%%-%.0f%% detected across samples",
                    min(qc_di$PctDetected), max(qc_di$PctDetected)))
  }, error = function(e) message("[run_de] detected/inferred QC skipped: ", e$message))
}

# ---- DE-LIMP-loadable session ----------------------------------------------
# Everything the DE-LIMP Shiny app needs is already in memory here. Writing it in
# server_session.R's schema lets any result be dropped straight into the GUI at
# https://delimp.stan-proteomics.org/ for interactive exploration.
if (method == "dpc" && exists("dat") && exists("y_protein")) {
  tryCatch({
    session_data <- list(
      raw_data     = dat,
      metadata     = meta,
      fit          = fit,
      # --block-scope within: `fit` holds each contrast's reporting fit column by column
      # (fit$contrast_model); the independent fit whole, for anything needing a
      # consistent per-protein variance (s2.post) on its contrasts.
      fit_independent = if (!is.null(block) && any(contrast_model == "independent")) fit_ind else NULL,
      y_protein    = y_protein,
      dpc_fit      = if (exists("dpcfit")) dpcfit else NULL,
      design       = design,
      qc_stats     = list(detected_vs_inferred = qc_di),
      # The same plain-R script written to reproducibility_log.R. DE-LIMP's
      # Reproducibility tab renders `repro_log` verbatim, so loading this session
      # into the GUI shows the exact code that produced it — previously NULL, which
      # left that tab blank for every skill-produced session.
      repro_log    = repro_lines,
      contrast     = paste(forms, collapse = ","),
      # Key name is DE-LIMP's session schema (server_de.R draws its volcano lines from
      # it), so it stays -- but the value is our reference line, not a significance cut.
      logfc_cutoff = logfc_ref,
      q_cutoff     = adjp_thr,
      saved_at     = Sys.time(),
      app_version  = "DE-LIMP v2.5"
    )
    rds <- file.path(outdir, "DE-LIMP_session.rds")
    saveRDS(session_data, rds)
    message(sprintf("[run_de] DE-LIMP_session.rds (%.1f MB) — load at https://delimp.stan-proteomics.org/",
                    file.size(rds) / 1e6))
  }, error = function(e) message("[run_de] DE-LIMP session not written: ", e$message))
}

message("\n[run_de] done. Results + provenance (methods.txt, sessionInfo.txt, de_provenance.json) in ",
        normalizePath(outdir))
