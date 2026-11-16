#!/usr/bin/env Rscript
# =============================================================================
# make_screen_data.R -- the data behind the protein screen ("show me the data for <protein>").
#
# For a finished limpa (DPC-Quant) analysis, writes <session>/output/screen_data/:
#   manifest.json          what was read (sha256 of every source), what was written, what was skipped
#   index.tsv              Protein.Group, Genes, Protein.Names, part, n_precursors, then
#                          P.Value__<contrast>, adj.P.Val__<contrast>, logFC__<contrast> per contrast (empty = untested)
#   parts/part-NNNN.json.gz   {Protein.Group: record}, at most 250 proteins per part, y_protein order
#
# One record per protein group: the delivered values (y_protein), the contrasts (the DE tables,
# with a CI from t and df_total), the precursors behind it (raw_data) and the DIA-NN report rows
# with the checks each failed. Later parts of the screen add blocks to the record (build_record()
# runs RECORD_BLOCKS in order) and to the manifest; a block that cannot be built is recorded in
# manifest$skipped, never dropped silently.
#
# Shares: how much of each sample's protein value rests on each precursor. limpa's own
# quantification (dpcQuant.EList, limpa 1.4.0) is redone one protein at a time -- the hyperparameters
# once from the whole matrix, then a refit of the protein, which must reproduce the delivered value
# within 1e-6. A share is the derivative of a protein value with respect to one found precursor value
# (what nudging that value by a tiny step does), taken exactly at limpa's optimum from limpa's own
# Hessian (implicit-function theorem; design section 9b; screen_shares.R). The derivative mirrors
# limpa 1.4.0's objective; manifest$shares$limpa_version records the limpa that ran, and every run
# checks the derivatives of 3 proteins against nudges of fully converged fits
# (manifest$shares$nudge_check) -- if they differ by more than 1e-4, or the check cannot run, every
# protein's shares are withheld. A protein whose refit misses is withheld with the reason
# (record$shares$withheld, manifest$shares$withheld), never shipped.
#
# Detection chances (design section 5.7): the chance limpa's fit gives each precursor of being found
# in each sample, at the fitted level (pf_fit, every cell), and the detection curve at each found
# value (pf_obs). They come from the shares' own refit -- its BFGS stopping point, the one checked
# against the delivered value -- so wherever shares are withheld, detection is withheld with the same
# reason; and each protein's are checked against limpa's own objective value at that fit
# (screen_shares.R detection_at_fit; manifest$detection$withheld).
#
# Weights (design section 5.4): limma's precision weights of the delivered fit, as run_de.R kept them
# (set_test_inputs.rds weights$independent), each sample's share of its group's weight, and the
# protein's residual and limma's prior degrees of freedom (DE-LIMP_session.rds fit). Never refitted:
# when they are not there, weights are null and manifest$skipped says why.
#
# Usage:
#   Rscript make_screen_data.R --session <dir> [--out <dir>] [--cpus N] [--proteins PG1,PG2]
# Exit 0 = written; exit 2 = refused (the reason on stderr and in manifest.json: supported false, with
# the analysis's own description in `analysis`, so the page builder can name the route it cannot show).
# Which analyses part 1 covers is decided in one place, screen_support(), from what the analysis says
# it is: limpa DPC-Quant (de_provenance.json pipeline_id "dpc") over a DIA-NN report (one that declares
# no quantity of its own: build_maxlfq.R declared_quantity()), with a plain group design. Refused: every
# other route (MaxLFQ + limma, an engine's own protein quantities -- FragPipe/IonQuant, directLFQ,
# AlphaDIA, Radiant -- Sage's highest-peptide route, SILAC), a dpc DE over another engine's report, a
# blocked or covariate design, a contrast that is not one group against one other ((B+C)/2-A, an
# interaction), and an analysis limpa quantified row by row (dpcQuantByRow).
# The report is the one the DE's own record names (de_provenance.json input + input_sha256), never a search
# this script picks. Its run state goes to output/screen_data_status.json (session.py status shows it).
# The report is the one the DE read (de_provenance.json `input`, checked by its input_sha256), wherever
# the session keeps it under output/; report.stats.tsv is read from the same folder.
# --proteins writes a partial folder (manifest$partial true), into its own --out: never into
# <session>/output/screen_data, where it would replace the full screen data. A refusal for wrong
# arguments (that, or --proteins naming a protein group the analysis lacks) writes nothing: the
# folder keeps the manifest it had.
#
# Test hooks (tests only, never set them on a real run); when set, the manifest records them:
#   SCREEN_DATA_SIGMA_SCALE=<x> multiplies each protein's sigma before the refits, to make them miss
#     the delivered values as a CPU-family change would (shares$test_hook_sigma_scale).
#   SCREEN_DATA_DERIVATIVE_SCALE=<x> multiplies every derivative, to make the nudge check fail as a
#     change in limpa's objective would (shares$test_hook_derivative_scale).
#   SCREEN_DATA_DROP_RECORD=<k> empties the k-th protein's record after the workers return, as a worker
#     killed (out of memory) leaves it; the run then stops before any manifest is written.
# =============================================================================

suppressMessages({ library(jsonlite) })

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
session <- getarg("--session")
if (is.null(session)) stop("Required: --session <dir>  (the analysis folder holding output/)")
session <- normalizePath(session, mustWork = TRUE)
output  <- file.path(session, "output")
tables  <- file.path(output, "tables")
out_dir <- getarg("--out", file.path(output, "screen_data"))
cpus    <- suppressWarnings(as.integer(getarg("--cpus", Sys.getenv("SLURM_CPUS_PER_TASK", "1"))))
if (is.na(cpus) || cpus < 1L) cpus <- 1L
only_pgs <- getarg("--proteins"); if (!is.null(only_pgs)) only_pgs <- trimws(strsplit(only_pgs, ",")[[1]])
t_start <- Sys.time()

.script_dir <- local({
  f <- grep("^--file=", commandArgs(), value = TRUE)[1]
  if (is.na(f)) getwd() else dirname(normalizePath(sub("^--file=", "", f), mustWork = FALSE))
})
# the one contaminant rule, version reader, CPU family + sha256, shares, run names, and build_maxlfq.R
# for declared_quantity() (the one reading of what a report says its numbers are)
for (nm in c("contaminants.R", "skill_version.R", "compute_record.R", "screen_shares.R", "run_name.R",
             "build_maxlfq.R")) {
  f <- file.path(.script_dir, nm)
  if (!file.exists(f)) stop(nm, " not found next to make_screen_data.R")
  source(f, encoding = "UTF-8")
}
WRITTEN_BY <- paste("make_screen_data.R", skill_label(skill_version(.script_dir)))
SCHEMA <- "screen_data/1"
PART_SIZE <- 250L
NO_LATER <- "added in a later part (design section 12)"

iso_now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
to_json <- function(x, pretty = FALSE)
  toJSON(x, auto_unbox = TRUE, digits = NA, na = "null", null = "null", pretty = pretty,
         json_verbatim = TRUE)       # a pre-rendered `json` value (the rows) goes in as it is
write_manifest <- function(m) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  writeLines(as.character(to_json(m, pretty = TRUE)), file.path(out_dir, "manifest.json"), useBytes = TRUE)
}
# The run's status record, output/screen_data_status.json (beside the other step records, not inside
# screen_data/, which a rebuild empties): state "running" | "built" | "refused" | "failed" and a one-line
# reason. `session.py status` shows it as an optional note, never a finish item. It describes the
# session's own screen data only: a run with --out or --proteins writes none.
status_on <- is.null(getarg("--out")) && is.null(only_pgs)
write_status <- function(state, reason = NULL) {
  if (!status_on) return(invisible())
  one <- if (is.null(reason)) NULL else substr(gsub("[[:space:]]+", " ", paste(reason, collapse = " ")), 1L, 400L)
  tryCatch(writeLines(as.character(to_json(list(schema = "screen_data_status/1", state = state, reason = one,
                                                route = if (!is.null(analysis)) analysis$pipeline_id,
                                                written_by = WRITTEN_BY, at = iso_now()))),
                      file.path(output, "screen_data_status.json"), useBytes = TRUE),
           error = function(e) message("[make_screen_data] could not write screen_data_status.json: ",
                                       conditionMessage(e)))
}
# Refuse: reason on stderr and in manifest.json (supported false), exit 2. A refusal for wrong
# arguments (`args_wrong`) says nothing about the analysis, so it leaves manifest.json as it was: a
# typo must not turn screen data already written there into a refusal.
# What the page says the analysis did, copied from de_provenance.json (repo rule 1: never hard-coded):
# the pipeline and its label, the quantities and their normalisation, the q-value columns and cutoffs,
# every filter as applied, the roll-up, the DE engine and the significance rule. A key the analysis did
# not record is null. Set once de_provenance.json is read; a refusal after that carries it too.
ANALYSIS_KEYS <- c("pipeline_id", "display_label", "normalisation", "quantities", "q_columns", "q_cutoffs",
                   "filters_applied", "rollup_method", "de_engine", "significance_rule")
analysis <- NULL
refuse <- function(reason, args_wrong = FALSE) {
  if (!args_wrong)
    write_manifest(c(list(schema = SCHEMA, written_by = WRITTEN_BY, written_at = iso_now(), supported = FALSE,
                          unsupported_reason = reason, partial = !is.null(only_pgs)),
                     if (!is.null(analysis)) list(analysis = analysis)))
  if (!args_wrong) write_status("refused", reason)
  mf <- file.path(out_dir, "manifest.json")
  message("[make_screen_data] refused: ", reason,
          if (args_wrong) { if (file.exists(mf)) paste0(" (", mf, " left as it was)") else " (nothing written)" })
  quit(save = "no", status = 2L)
}
# The same folder, however it is spelled: one that exists is resolved whole (a trailing /., a link),
# one that does not yet by its parent.
same_dir <- function(a, b) {
  n <- function(p) if (dir.exists(p)) normalizePath(p) else file.path(normalizePath(dirname(p), mustWork = FALSE), basename(p))
  identical(n(a), n(b))
}
if (!is.null(only_pgs) && same_dir(out_dir, file.path(output, "screen_data")))
  refuse(paste("--proteins writes a partial folder; give it its own --out, not", file.path(output, "screen_data"),
               "(the session's full screen data, which it would replace)"), args_wrong = TRUE)
# ... and never into any folder that holds full screen data, however it is reached (a spelling the path check
# misses, another letter case on a case-insensitive disk, a copy): its manifest says supported, not partial.
if (!is.null(only_pgs) && file.exists(file.path(out_dir, "manifest.json"))) {
  old <- tryCatch(fromJSON(file.path(out_dir, "manifest.json"), simplifyVector = FALSE),
                  error = function(e) refuse(paste0(file.path(out_dir, "manifest.json"), " cannot be read (",
                    conditionMessage(e), "); remove it, or give --proteins its own --out"), args_wrong = TRUE))
  if (isTRUE(old$supported) && !isTRUE(old$partial))     # no `partial` key (older screen data) counts as full
    refuse(paste(out_dir, "holds full screen data (its manifest.json: supported, not partial); --proteins writes a",
                 "partial folder, so give it its own --out"), args_wrong = TRUE)
}

# Arguments are checked: from here a stop is about the analysis, so it is recorded.
write_status("running")
options(error = function() {
  write_status("failed", sub("^Error[^:]*: ", "", geterrmessage()))
  quit(save = "no", status = 1L)
})
# ---- what the analysis is ------------------------------------------------------
prov_f <- file.path(tables, "de_provenance.json")
sess_f <- file.path(tables, "DE-LIMP_session.rds")
if (!file.exists(prov_f)) refuse("tables/de_provenance.json is missing: this is not a finished skill analysis")
prov <- fromJSON(prov_f, simplifyVector = FALSE)
analysis <- setNames(lapply(ANALYSIS_KEYS, function(k) prov[[k]]), ANALYSIS_KEYS)
output_n <- normalizePath(output)
rel_output <- function(p) {         # a path under output/ as the manifest names it (relative to output/)
  pre <- paste0(output_n, "/")
  if (startsWith(p, pre)) substring(p, nchar(pre) + 1L) else p
}

# The report the DE read, from the DE's OWN record (de_provenance.json `input` path + `input_sha256`) and
# nothing else: the screen does not re-derive which search a session "chose" (session.py search_choice is that
# one definition) and does not look through the session for a report. A session copied or moved after its DE
# (to the delivery folder, rsynced between HIVE and a laptop) has a record that names the old location, so
# when the recorded file is not inside this session's output/, the SAME RELATIVE PATH under this session's
# output/ is tried (the part of the recorded path after an `/output/`). Whichever file is used must hash to
# input_sha256. A record without the input path or hash (an older analysis) is refused; so is a file that
# changed ("the report changed since the DE"), or one not found at either place.
de_report <- function() {
  want <- prov$input_sha256
  inp <- if (is.character(prov$input) && length(prov$input) == 1L && nzchar(prov$input)) prov$input else ""
  if (!is.character(want) || length(want) != 1L || !nzchar(want) || !nzchar(inp))
    refuse(paste("tables/de_provenance.json does not record the report this analysis read (its input path and",
                 "input_sha256): an older analysis; run the DE again, then the screen"))
  inside <- function(p) startsWith(normalizePath(p, mustWork = FALSE), paste0(output_n, "/"))
  cand <- if (file.exists(inp) && inside(inp)) inp else {
    # the recorded path's tail after each "/output/" (the last one first), under this session's output/
    m <- gregexpr("/output/", inp, fixed = TRUE)[[1]]
    tails <- if (m[1] > 0) rev(substring(inp, m + nchar("/output/"))) else character()
    rel <- file.path(output_n, tails)
    rel <- rel[file.exists(rel) & vapply(rel, inside, NA)]
    if (length(rel)) rel[1] else NA_character_
  }
  if (is.na(cand)) {
    if (file.exists(inp))
      refuse(sprintf(paste("the report this analysis read (%s) is outside this session's output/ folder, and the",
                           "same path under this session's output/ is not there; run the DE again, then the screen"), inp))
    refuse(sprintf(paste("the report this analysis read is not at %s, nor at the same path under this session's",
                         "output/ (the search folder was renamed or removed since the DE ran); run the DE again,",
                         "then the screen"), inp))
  }
  if (!identical(file_sha256(cand), want))
    refuse(sprintf(paste("%s is not the report this analysis read (input_sha256 %s): the report changed since the",
                         "DE ran; the screen reads only the report the DE read"), cand, substr(want, 1, 12)))
  normalizePath(cand)
}

# ---- which analyses part 1 covers: ONE decision, from what the analysis says it is ----------------
# (repo rule 1) The route is de_provenance.json pipeline_id, named by its own descriptor's display_label
# and rollup_method; whether the report is DIA-NN's is the report's own declaration (parquet metadata
# delimp.quantity.*, build_maxlfq.R declared_quantity(): an adapter -- Sage, FragPipe, AlphaDIA, Radiant
# -- declares what its numbers are; DIA-NN's report declares nothing). Every block of part 1 reads limpa's
# objects and DIA-NN's precursor rows, so anything else is refused here, with its own name. NULL = covered.
SCREEN_PIPELINES <- "dpc"
LATER_PART <- "the screen for other routes and designs is a later part (design section 12, part 5)"
# A contrast the page can describe is "one group minus one group" (numerator B, denominator A): a group name
# may itself hold "-", so every "-" is tried. NA/NA = anything else ((B+C)/2-A, an interaction, a sum).
split_contrast <- function(cn, lv) {
  pos <- gregexpr("-", cn, fixed = TRUE)[[1]]
  for (p in pos[pos > 0]) {
    l <- substr(cn, 1L, p - 1L); r <- substr(cn, p + 1L, nchar(cn))
    if (l %in% lv && r %in% lv && !identical(l, r)) return(c(l, r))
  }
  c(NA_character_, NA_character_)
}
screen_support <- function(prov, declared, groups = NULL) {
  id <- if (is.character(prov$pipeline_id)) prov$pipeline_id else "unknown"
  what <- paste0(if (is.character(prov$display_label)) prov$display_label else "an unnamed quantification",
                 " (pipeline_id ", id,
                 if (is.character(prov$rollup_method)) paste0("; ", prov$rollup_method) else "", ")")
  if (!id %in% SCREEN_PIPELINES)
    return(sprintf(paste("this analysis's quantification is %s; part 1 of the protein screen reads limpa",
                         "(DPC-Quant) analyses only, and %s"), what, LATER_PART))
  if (!is.null(declared))
    return(sprintf(paste("this analysis ran limpa (DPC-Quant) on a report written by %s (it declares its own",
                         "quantity: %s), not by DIA-NN; part 1 of the protein screen reads DIA-NN's precursor",
                         "rows, and %s"), declared$label %||% "another engine's adapter",
                   declared$value %||% "not stated", LATER_PART))
  if (isTRUE(prov$block$applied) || length(prov$covariates$terms))
    return(paste("the protein screen covers group comparisons; this design (blocked or with covariates) is",
                 "not one, and", LATER_PART))
  if (!is.null(groups)) {                  # the design's groups are known: every contrast must be one against one
    bad <- Filter(function(cn) anyNA(split_contrast(cn, unique(groups))), as.character(unlist(prov$contrasts)))
    if (length(bad))
      return(sprintf(paste("the protein screen covers a comparison of one group with one other; this analysis's",
                           "contrast %s is not one (a combination of groups, or an interaction), and %s"),
                     paste(shQuote(bad), collapse = ", "), LATER_PART))
  }
  NULL
}
if (!is.null(why <- screen_support(prov, NULL))) refuse(why)       # the route and the design first: no report read
if (!file.exists(sess_f)) refuse("tables/DE-LIMP_session.rds is missing")
report_f <- de_report()
declared <- tryCatch(declared_quantity(arrow::open_dataset(report_f, format = "parquet")),
                     error = function(e) refuse(paste0(rel_output(report_f), " cannot be read for what it declares (",
                                                       conditionMessage(e), ")")))
if (!is.null(why <- screen_support(prov, declared))) refuse(why)
report_rel <- rel_output(report_f)
stats_f <- file.path(dirname(report_f), "report.stats.tsv")
stats_rel <- rel_output(stats_f)

suppressMessages(library(limma))
sess <- readRDS(sess_f)
y    <- sess$y_protein
raw  <- sess$raw_data
# When every protein has one precursor, limpa's dpcQuant hands the analysis to dpcQuantByRow, which
# quantifies each precursor alone (y_protein's rows are then precursors): not the protein model here.
if (!anyDuplicated(as.character(raw$genes$Protein.Group)))
  refuse(paste("limpa quantified this analysis row by row (dpcQuantByRow: every protein has exactly one",
               "precursor); the protein screen does not cover that"))
meta <- sess$metadata
samples <- colnames(y$E)
# run_de.R matches the conditions' File.Name to the report's Run exactly (no extension is stripped),
# and aligns `metadata` to the matrix columns; so a sample's group is the metadata row of that name.
groups <- as.character(meta$Group[match(samples, meta$File.Name)])
if (anyNA(groups)) refuse("a sample of the analysis has no row in the session's metadata")

# ---- the contrasts ---------------------------------------------------------------
if (!is.null(why <- screen_support(prov, declared, groups))) refuse(why)   # the contrasts need the design's groups
contrast_names <- unlist(prov$contrasts)
set_f <- file.path(tables, "set_test_inputs.rds")
if (!file.exists(set_f)) refuse("tables/set_test_inputs.rds is missing")
set_in <- readRDS(set_f)
stat <- set_in$stat
sources <- c("tables/de_provenance.json", "tables/DE-LIMP_session.rds", "tables/set_test_inputs.rds")
de_tabs <- list(); contrasts <- list()
for (cn in contrast_names) {
  tinfo <- prov$de_tables[[cn]]
  f <- file.path(tables, tinfo$file)
  if (!file.exists(f)) refuse(sprintf("the DE table of %s (%s) is missing", cn, tinfo$file))
  sources <- c(sources, file.path("tables", tinfo$file))
  tb <- utils::read.csv(f, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
  de_tabs[[cn]] <- tb[match(rownames(y$E), tb$Protein.Group), , drop = FALSE]
  ng <- split_contrast(cn, unique(groups))
  contrasts[[length(contrasts) + 1L]] <- list(name = cn, numerator = ng[1], denominator = ng[2],
    table = tinfo$file, n_tested = tinfo$n_tested, n_sig = tinfo$n_significant, adjp = prov$adjp)
}

# ---- the DIA-NN report: only the columns a record uses --------------------------
q_cols <- unlist(prov$q_columns); q_cuts <- as.numeric(unlist(prov$q_cutoffs))
REPORT_COLS <- c("Run", "Precursor.Id", "Protein.Group", "Protein.Ids", "Proteotypic", "Precursor.Quantity",
                 "Normalisation.Factor", "Precursor.Normalised", "Quantity.Quality", "Empirical.Quality",
                 "Ms1.Area", "RT", q_cols, "PG.MaxLFQ", "PG.MaxLFQ.Quality")
have <- nanoparquet::read_parquet_schema(report_f)$name
rep <- as.data.frame(nanoparquet::read_parquet(report_f, col_select = intersect(unique(REPORT_COLS), have)))
sources <- c(sources, report_rel)
rep_idx <- split(seq_len(nrow(rep)), factor(rep$Protein.Group, levels = unique(rep$Protein.Group)))
raw_idx <- split(seq_len(nrow(raw$E)), factor(raw$genes$Protein.Group, levels = unique(raw$genes$Protein.Group)))
raw_ids <- rownames(raw$E)

# ---- per-run search statistics (optional) ----------------------------------------
skipped <- list()
note_skip <- function(block, reason) skipped[[length(skipped) + 1L]] <<- list(block = block, reason = reason)
# Run an optional block; on error record {block, reason} and carry on with `fallback`.
optional <- function(block, expr, fallback = NULL)
  tryCatch(expr, error = function(e) { note_skip(block, conditionMessage(e)); fallback })

run_stats <- optional("run_stats", {
  if (!file.exists(stats_f)) stop(stats_rel, " does not exist")
  st <- utils::read.delim(stats_f, check.names = FALSE, stringsAsFactors = FALSE)
  sources <- c(sources, stats_rel)
  # DIA-NN writes File.Name as the raw file's path (either slash, with its extension); the report's Run
  # is run_name() of it (run_name.R, the twin of check_report_runs.py's)
  key <- run_name(st$File.Name)
  lapply(setNames(samples, samples), function(s) {
    i <- match(s, key)
    if (is.na(i)) stop("no statistics for run ", s)
    as.list(st[i, setdiff(names(st), "File.Name"), drop = FALSE])
  })
})

t_read <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))   # session, report, DE tables, stats
# ---- limma's weights (design section 5.4) ------------------------------------------------
# The precision weights the delivered fit ran with -- the unblocked dpc fit's vooma weights, which
# run_de.R keeps in set_test_inputs.rds (weights$independent) -- and the degrees of freedom of the
# session's fit. Read, never refitted: when they are not there, every record's weights are null.
wts <- optional("weights", {
  w <- set_in$weights$independent
  if (is.null(w)) stop("tables/set_test_inputs.rds holds no weights of the delivered fit (weights$independent)")
  if (!is.matrix(w) || !setequal(rownames(w), rownames(y$E)) || !setequal(colnames(w), samples))
    stop(paste("the weights in tables/set_test_inputs.rds (weights$independent) are not one row per",
               "protein group and one column per sample of this analysis"))
  f <- sess$fit; fpg <- rownames(f$coefficients); dp <- as.numeric(f$df.prior)
  if (!setequal(fpg, rownames(y$E)) || length(f$df.residual) != length(fpg) || !length(dp) %in% c(1L, length(fpg)))
    stop(paste("the fit in tables/DE-LIMP_session.rds does not give residual and prior degrees of freedom",
               "for the protein groups of this analysis"))
  list(w = w[, samples, drop = FALSE], df_residual = setNames(as.numeric(f$df.residual), fpg),
       df_prior = if (length(dp) == 1L) dp else setNames(dp, fpg))
})
if (!is.null(wts) && any(is.infinite(wts$df_prior)))
  note_skip("weights.df_prior", paste("limma's prior degrees of freedom is infinite (its empirical Bayes step",
    "found no spread in the proteins' variances beyond chance); JSON has no infinity, so df_prior is null"))

# ---- the protein record ------------------------------------------------------------
r4 <- function(x) round(as.numeric(x), 4)
s6 <- function(x) signif(as.numeric(x), 6)
s3 <- function(x) signif(as.numeric(x), 3)
arr <- function(x) I(unname(x))          # a JSON array even at length 1

# limpa's count of each protein's precursors (its rows in raw_data): y_protein$genes$NPrec from limpa 1.4,
# NPeptides before (the same count, renamed). Recorded as n_precursors; null, with the reason in
# manifest$skipped, when the analysis has neither column.
nprec_col <- intersect(c("NPrec", "NPeptides"), names(y$genes))[1]
if (is.na(nprec_col))
  note_skip("n_precursors", paste("y_protein$genes has neither NPrec (limpa 1.4 and later) nor NPeptides (older",
                                  "limpa), so each record's n_precursors is null"))

# Context shared by the record blocks (read-only).
ctx <- list(samples = samples, groups = groups, y = y, raw = raw, rep = rep, rep_idx = rep_idx,
            raw_idx = raw_idx, raw_ids = raw_ids, de_tabs = de_tabs, stat = stat,
            contrasts = contrasts, q_cols = q_cols, q_cuts = q_cuts, wts = wts, nprec_col = nprec_col)

block_identity <- function(rec, pg, i, ctx) {
  g <- ctx$y$genes
  rec$pg <- pg
  rec$genes <- as.character(g$Genes[i]); rec$names <- as.character(g$Protein.Names[i])
  rec$n_precursors <- if (is.na(ctx$nprec_col)) NA_integer_ else as.integer(g[[ctx$nprec_col]][i])
  rec["sigma"] <- list(NULL)       # a later part fills it
  rec
}
block_de <- function(rec, pg, i, ctx) {
  rec$de <- setNames(lapply(ctx$contrasts, function(co) {
    cn <- co$name; tb <- ctx$de_tabs[[cn]]; st <- ctx$stat[[cn]]
    det <- grep("^Detected_", names(tb), value = TRUE)
    detected <- setNames(lapply(det, function(k) as.character(tb[[k]][i])), sub("^Detected_", "", det))
    lfc <- tb$logFC[i]; tt <- tb$t[i]
    df <- if (!is.null(st$df_total)) {
      d <- st$df_total      # per protein, in the order of the fit's t (named by protein group)
      if (length(d) == 1L) d else d[match(pg, names(st$t))]
    } else NA_real_
    tested <- is.finite(lfc) && is.finite(tt)
    se <- if (tested && tt != 0) lfc / tt else NA_real_
    q <- if (is.finite(df)) qt(0.975, df) else NA_real_
    list(tested = tested, logFC = if (tested) lfc else NA_real_,
         CI.L = if (tested) lfc - q * se else NA_real_, CI.R = if (tested) lfc + q * se else NA_real_,
         t = if (tested) tt else NA_real_, P.Value = if (tested) tb$P.Value[i] else NA_real_,
         adj.P.Val = if (tested) tb$adj.P.Val[i] else NA_real_, df_total = df,
         evidence = as.character(tb$Evidence[i]), detected = detected)
  }), vapply(ctx$contrasts, function(co) co$name, ""))
  rec
}
block_values <- function(rec, pg, i, ctx) {
  rec$values <- list(protein = arr(r4(ctx$y$E[i, ])), se = arr(r4(ctx$y$other$standard.error[i, ])),
                     n_obs = arr(as.integer(ctx$y$other$n.observations[i, ])))
  rec
}
block_precursors <- function(rec, pg, i, ctx) {
  rr <- ctx$raw_idx[[pg]]
  rec$precursors <- lapply(rr, function(k) list(id = ctx$raw_ids[k], log2 = arr(r4(ctx$raw$E[k, ]))))
  rec
}
block_rows <- function(rec, pg, i, ctx) {
  ri <- ctx$rep_idx[[pg]]; d <- ctx$rep[ri, , drop = FALSE]
  d <- d[order(match(d$Run, ctx$samples), d$Precursor.Id), , drop = FALSE]
  col <- function(nm, f = identity) if (nm %in% names(d)) f(d[[nm]]) else rep(NA, nrow(d))
  sample <- ifelse(d$Run %in% ctx$samples, d$Run, NA_character_)
  # used: a precursor of this protein in the analysed matrix with a value in this run
  rr <- ctx$raw_idx[[pg]]; ridx <- match(d$Precursor.Id, ctx$raw_ids[rr]); cidx <- match(sample, ctx$samples)
  used <- rep(FALSE, nrow(d)); ok <- !is.na(ridx) & !is.na(cidx)
  if (any(ok)) used[ok] <- !is.na(ctx$raw$E[cbind(rr[ridx[ok]], cidx[ok])])
  failed <- lapply(seq_len(nrow(d)), function(j) {
    f <- ctx$q_cols[vapply(seq_along(ctx$q_cols), function(k) {
      v <- if (ctx$q_cols[k] %in% names(d)) d[[ctx$q_cols[k]]][j] else NA_real_
      !is.na(v) && v > ctx$q_cuts[k] }, NA)]
    I(as.character(f))
  })
  cols <- list(Run = d$Run, sample = sample, Precursor.Id = d$Precursor.Id, Protein.Ids = col("Protein.Ids"),
               Proteotypic = col("Proteotypic"), Precursor.Quantity = col("Precursor.Quantity", s6),
               Normalisation.Factor = col("Normalisation.Factor", s6),
               Precursor.Normalised = col("Precursor.Normalised", s6),
               Quantity.Quality = col("Quantity.Quality", s6), Empirical.Quality = col("Empirical.Quality", s6),
               Ms1.Area = col("Ms1.Area", s6), RT = col("RT", s6))
  for (qc in ctx$q_cols) cols[[qc]] <- col(qc, s3)
  cols$PG.MaxLFQ <- col("PG.MaxLFQ", s6); cols$PG.MaxLFQ.Quality <- col("PG.MaxLFQ.Quality", s6)
  cols$failed <- failed
  cols$contaminant <- is_contaminant(col("Protein.Ids"))
  cols$used <- used
  # The rows go out as one data-frame call: jsonlite walks a data frame in C, but takes ~110 s on a
  # list of 700k per-row lists (HeLa, 2026-10-28). Same bytes: every cell keeps its own scalar form,
  # NA is null, `failed` is an array even when empty or of one name.
  rec$rows <- if (nrow(d)) structure(as.character(to_json(structure(cols, class = "data.frame",
                                                     row.names = .set_row_names(nrow(d))))), class = "json")
              else list()
  rec
}
# sigma and shares are computed before the records (section "shares" below); the record carries
# them. ctx$sigma is NULL when limpa's hyperparameters could not be reproduced.
block_shares <- function(rec, pg, i, ctx) {
  rec["sigma"] <- list(ctx$sigma[[pg]])
  rec$shares <- ctx$shares[[pg]]
  rec
}
# The detection chances, computed with the shares (same fit, same withholding).
block_detection <- function(rec, pg, i, ctx) {
  rec$detection <- ctx$detection[[pg]]
  rec
}
# limma's weights of this protein (wts, above): w per sample, each sample's share of its group's total
# weight, its residual df and limma's prior df. null when the weights were not there (manifest$skipped).
block_weights <- function(rec, pg, i, ctx) {
  if (is.null(ctx$wts)) { rec["weights"] <- list(NULL); return(rec) }
  w <- ctx$wts$w[pg, ]
  dp <- ctx$wts$df_prior; if (length(dp) > 1L) dp <- dp[[pg]]
  rec$weights <- list(w = arr(w), share_in_group = arr(w / stats::ave(w, ctx$groups, FUN = sum)),
                      df_residual = ctx$wts$df_residual[[pg]], df_prior = dp)
  rec
}
# Blocks run in this order.
RECORD_BLOCKS <- list(identity = block_identity, de = block_de, values = block_values,
                      precursors = block_precursors, rows = block_rows, shares = block_shares,
                      detection = block_detection, weights = block_weights)
build_record <- function(pg, i, ctx) {
  rec <- list()
  for (b in RECORD_BLOCKS) rec <- b(rec, pg, i, ctx)
  rec
}

# ---- which proteins ------------------------------------------------------------------
pgs_all <- rownames(y$E)
if (!is.null(only_pgs)) {
  miss <- setdiff(only_pgs, pgs_all)
  if (length(miss)) refuse(sprintf("--proteins names protein groups not in this analysis: %s",
                                   paste(utils::head(miss, 5), collapse = ", ")), args_wrong = TRUE)
  pgs_all <- pgs_all[pgs_all %in% only_pgs]
}
# Run f over x: forked workers on unix when --cpus > 1, else one at a time.
par_lapply <- function(x, f)
  if (cpus > 1L && .Platform$OS.type == "unix") parallel::mclapply(x, f, mc.cores = cpus) else lapply(x, f)
# The node this runs on: named in the manifest, and in the reason when limpa's numbers do not reproduce.
here_family <- tryCatch({
  fam <- slurm_cpu_family()$fam
  if (length(fam)) fam[1] else "unknown"
}, error = function(e) { note_skip("compute", conditionMessage(e)); "unknown" })
de_family <- if (length(prov$compute$cpu_family)) as.character(prov$compute$cpu_family)[1] else "not recorded"

# ---- shares --------------------------------------------------------------------------
# How much of each sample's protein value rests on each precursor: limpa's own quantification
# (dpcQuant.EList, limpa 1.4.0) redone one protein at a time, then differentiated at its optimum.
# Why not nudges: a refit stops where BFGS's tolerance lets it, short of the optimum, and a nudge
# reads that stopping point as a share -- plain refits at a step of 0.01 put one HeLa protein's total
# share at -5.7, and Newton-polished ones at 1e-3 were still off by up to 9e-4 (2026-10-28). The
# derivative below agrees with a tightly converged nudge reference to 1e-10 on all 3834 HeLa proteins.
SHARE_TOL <- 1e-6     # the refit must reproduce the delivered value to this
SETTLE_TOL <- 1e-5    # a Newton step from the second fit may move it at most this (HeLa: 1.7e-7)
CHECK_N <- 3L         # proteins in the run-time nudge check
CHECK_CELLS <- 24L    # found values nudged per checked protein, at most (evenly spaced when it has more)
CHECK_STEP <- 1e-3; CHECK_TOL <- 1e-4    # central-difference step; the largest difference allowed
DETECT_TOL <- 1e-8    # limpa's objective rebuilt from the detection chances vs its own value, relative
# A test hook (see the header): its value, NA when not set.
test_hook <- function(var, what) {
  v <- Sys.getenv(var, "")
  if (!nzchar(v)) return(NA_real_)
  x <- suppressWarnings(as.numeric(v))
  if (!is.finite(x) || x <= 0) stop(var, " (a test hook) must be a positive number, not ", v)
  message("[make_screen_data] TEST HOOK: ", var, "=", v, " ", what)
  x
}
sigma_scale <- test_hook("SCREEN_DATA_SIGMA_SCALE", "multiplies sigma before the refits")
deriv_scale <- test_hook("SCREEN_DATA_DERIVATIVE_SCALE", "multiplies every derivative")
drop_record <- test_hook("SCREEN_DATA_DROP_RECORD", "empties that protein's record, as a killed worker would")
t_shares <- Sys.time()
# limpa's hyperparameters as dpcQuant.EList computes them: rows ordered by protein, then prior.sd and
# prior.logFC floored at 1; they must match the ones y_protein carries. When they cannot be had,
# `hyper` is list(reason = ...) and every protein's shares are withheld with that reason.
hyper <- local({
  pid <- as.character(raw$genes$Protein.Group)     # some protein has 2+ precursors (refused above otherwise)
  if (length(y$dpc) != 2L) return(list(reason = "y_protein does not record the DPC limpa used"))
  o <- order(pid)
  h <- tryCatch(limpa::dpcQuantHyperparam(raw[o, ], protein.id = pid[o], dpc.slope = y$dpc[2]),
                error = function(e) e)
  if (inherits(h, "error"))
    return(list(reason = paste("limpa's hyperparameters could not be recomputed:", conditionMessage(h))))
  h$prior.sd <- max(h$prior.sd, 1); h$prior.logFC <- max(h$prior.logFC, 1)
  num <- function(x) if (length(x) == 1L) as.numeric(x) else NA_real_
  dh <- max(abs(c(h$prior.mean, h$prior.sd, h$prior.logFC) -
                c(num(y$prior.mean), num(y$prior.sd), num(y$prior.logFC))))
  if (!is.finite(dh) || dh > 1e-8)
    return(list(reason = sprintf(paste("limpa's hyperparameters recomputed here differ from the delivered",
      "ones by %s, so no refit here describes the analysis (compute.de %s, compute.here %s)"),
      format(signif(dh, 3)), de_family, here_family)))
  if (is.null(names(h$sigma))) names(h$sigma) <- unique(pid[o])
  # BFGS's start values, which peptides2Proteins imputes over the whole ordered matrix (its column
  # tilts use every protein); back in raw_data row order.
  h$yimp <- limpa::imputeByExpTilt(raw$E[o, , drop = FALSE], dpc.slope = y$dpc[2],
                                   prior.logfc = h$prior.logFC)[order(o), , drop = FALSE]
  h$dpc <- y$dpc
  h$scale <- if (is.na(sigma_scale)) 1 else sigma_scale
  h$dscale <- if (is.na(deriv_scale)) 1 else deriv_scale
  h
})
# One protein's derivatives (screen_shares.R): the refit from limpa's start values must reproduce the
# delivered value within SHARE_TOL, or the protein is withheld. Returns list(E, sg, f1, Wm, miss), f1
# the refit (shares_fit), Wm as shares_derivative's (column p + (k - 1) * np is d values / d E[p, k]),
# or list(withheld = <reason>).
deriv_of <- function(pg, h, ctx) {
  rr <- ctx$raw_idx[[pg]]
  E <- ctx$raw$E[rr, , drop = FALSE]
  sg <- h$sigma[[pg]] * h$scale
  f1 <- shares_fit(E, sg, h, shares_start(E, h$yimp[rr, , drop = FALSE]))
  miss <- max(abs(f1$par.bgfs[seq_len(ncol(E))] - ctx$y$E[pg, ]))
  if (!is.finite(miss) || miss > SHARE_TOL)
    return(list(withheld = sprintf("refit differs from the delivered value by %s", format(signif(miss, 3)))))
  d <- shares_derivative(E, sg, h, f1, SETTLE_TOL)
  if (!is.null(d$withheld)) return(d)
  list(E = E, sg = sg, f1 = f1, Wm = d$Wm * h$dscale, miss = miss)
}
# One protein's shares and detection chances, from the one refit (deriv_of): wherever its shares are
# withheld, so is its detection, with the same reason. Shares: own[p][j] = d value_j / d E[p, j] (null
# at a gap), total[j] = sum of value_j's derivatives, other[j] = total[j] - sum over p of own[p][j].
# Detection (screen_shares.R detection_at_fit, at the refit's BFGS stop): pf_fit[p][j], pf_obs[p][j]
# (null at a gap). 4 decimals. objective_diff: detection_at_fit's check (NA when not computed).
screen_of <- function(pg, h, ctx) {
  d <- deriv_of(pg, h, ctx)
  if (!is.null(d$withheld)) return(list(shares = d, detection = d, objective_diff = NA_real_))
  by_row <- function(m) lapply(seq_len(nrow(m)), function(p) arr(r4(m[p, ])))
  s <- shares_split(d$Wm, d$E)
  dt <- detection_at_fit(d$E, d$f1, d$sg, h, DETECT_TOL)
  list(shares = list(own = by_row(s$own), other = arr(r4(s$other)), total = arr(r4(s$total)),
                     refit_max_diff = signif(d$miss, 3)),
       detection = if (!is.null(dt$withheld)) dt else list(pf_fit = by_row(dt$pf_fit), pf_obs = by_row(dt$pf_obs)),
       objective_diff = if (is.null(dt$withheld)) dt$objective_diff else NA_real_)
}
# The run-time nudge check. The derivatives mirror limpa 1.4.0's objective, and nothing else ties them
# to the limpa that ran: a limpa that weighted the found values, or changed its design, would still
# reproduce the delivered values. So CHECK_N proteins' derivatives are compared with central
# differences of fully converged fits (shares_nudge). Chosen deterministically, first those with gaps
# and 3+ precursors (the detection curve and the prior both matter), with at most CHECK_CELLS found
# values (all of them nudged), then a sample with no found value, then more precursors, then name.
pick_check <- function(pgs, ctx) {
  f <- lapply(pgs, function(pg) !is.na(ctx$raw$E[ctx$raw_idx[[pg]], , drop = FALSE]))
  np <- vapply(f, nrow, 0L); n_found <- vapply(f, sum, 0L); gaps <- vapply(f, function(x) sum(!x), 0L)
  empty <- vapply(f, function(x) sum(colSums(x) == 0L), 0L)
  o <- order(-(gaps > 0), -(np >= 3L), -(n_found <= CHECK_CELLS), -(empty > 0), -np, pgs, method = "radix")
  utils::head(pgs[o], CHECK_N)
}
# One protein: the largest difference between its derivatives and the nudges, over the nudged
# columns, and over own/other/total when every found value was nudged.
nudge_check_of <- function(pg, h, ctx) {
  d <- deriv_of(pg, h, ctx)
  if (!is.null(d$withheld)) stop(d$withheld)
  found <- which(!is.na(as.vector(d$E)))
  cells <- if (length(found) > CHECK_CELLS) found[unique(round(seq(1, length(found), length.out = CHECK_CELLS)))] else found
  Wn <- shares_nudge(d$E, d$sg, h, cells, CHECK_STEP)
  diff <- max(abs(Wn - d$Wm[, cells, drop = FALSE]))
  if (length(cells) == length(found)) {
    Wf <- d$Wm; Wf[, cells] <- Wn
    a <- shares_split(d$Wm, d$E); b <- shares_split(Wf, d$E)
    diff <- max(diff, abs(a$own - b$own), abs(a$other - b$other), abs(a$total - b$total), na.rm = TRUE)
  }
  list(cells = length(cells), found = length(found), max_diff = diff)
}
nudge_check <- list(proteins = I(character(0)), max_diff = NULL, tolerance = CHECK_TOL, passed = NULL)
objective_max_diff <- NULL    # detection_at_fit's largest objective difference, over the proteins it passed
if (is.null(hyper$reason)) {
  quant <- par_lapply(pgs_all, function(pg) tryCatch(screen_of(pg, hyper, ctx), error = function(e) {
    w <- list(withheld = paste("the refit failed:", conditionMessage(e)))
    list(shares = w, detection = w)
  }))
  dead <- vapply(quant, function(s) is.null(s) || inherits(s, "try-error"), NA)
  if (any(dead)) stop("computing the shares of ", pgs_all[which(dead)[1]], " failed: its worker returned nothing")
  shares <- lapply(quant, `[[`, "shares"); detection <- lapply(quant, `[[`, "detection")
  names(shares) <- names(detection) <- pgs_all
  od <- unlist(lapply(quant, `[[`, "objective_diff"))
  objective_max_diff <- if (any(is.finite(od))) signif(max(od, na.rm = TRUE), 3)
  ctx$sigma <- setNames(as.list(unname(hyper$sigma[pgs_all])), pgs_all)   # NA (null) if limpa gave none
  computed <- pgs_all[vapply(shares, function(s) is.null(s$withheld), NA)]
  t_check <- Sys.time()
  if (!length(computed)) {
    nudge_check$reason <- "no protein's shares were computed, so there was nothing to check"
  } else {
    chk <- pick_check(computed, ctx)
    res <- par_lapply(chk, function(pg) tryCatch(nudge_check_of(pg, hyper, ctx),
                                                 error = function(e) list(error = conditionMessage(e))))
    got <- function(r, k) if (is.list(r) && is.numeric(r[[k]]) && length(r[[k]]) == 1L) r[[k]] else NA_real_
    err <- vapply(res, function(r) if (is.finite(got(r, "max_diff"))) NA_character_
                  else if (is.list(r) && is.character(r$error)) r$error else "no result", "")
    m <- if (all(is.na(err))) max(vapply(res, got, 0, "max_diff")) else NA_real_
    nudge_check <- list(proteins = I(chk), cells_checked = sum(vapply(res, got, 0, "cells"), na.rm = TRUE),
      cells_found = sum(vapply(res, got, 0, "found"), na.rm = TRUE), step = CHECK_STEP,
      max_diff = if (is.na(m)) NULL else signif(m, 3), tolerance = CHECK_TOL, passed = !is.na(m) && m <= CHECK_TOL,
      elapsed_s = round(as.numeric(difftime(Sys.time(), t_check, units = "secs")), 1))
    if (!nudge_check$passed) {
      why <- if (is.na(m))
        sprintf("the run-time nudge check of the shares could not be run (%s: %s)", chk[which(!is.na(err))[1]],
                err[!is.na(err)][1])
      else sprintf(paste("the run-time nudge check failed: limpa's derivatives differ from nudges of its converged",
                         "fits by %s (tolerance %s) on %s, so they do not describe limpa %s"),
                   format(signif(m, 3)), format(CHECK_TOL), paste(chk, collapse = ", "),
                   as.character(utils::packageVersion("limpa")))
      nudge_check$reason <- why
      message("[make_screen_data] shares and detection chances withheld for every protein: ", why)
      shares[computed] <- rep(list(list(withheld = why)), length(computed))
      detection[computed] <- shares[computed]
    }
  }
} else {
  message("[make_screen_data] shares and detection chances withheld for every protein: ", hyper$reason)
  shares <- rep(list(list(withheld = hyper$reason)), length(pgs_all))
  detection <- shares
  nudge_check$reason <- paste("no shares were computed:", hyper$reason)
}
names(shares) <- names(detection) <- pgs_all
ctx$shares <- shares; ctx$detection <- detection
withheld_of <- function(x) unname(Filter(Negate(is.null), lapply(pgs_all, function(pg)
  if (!is.null(x[[pg]]$withheld)) list(pg = pg, reason = x[[pg]]$withheld))))
shares_withheld <- withheld_of(shares); detection_withheld <- withheld_of(detection)
t_shares <- as.numeric(difftime(Sys.time(), t_shares, units = "secs"))
message(sprintf(paste("[make_screen_data] shares: %d of %d proteins withheld (manifest shares.withheld);",
                      "detection chances: %d withheld (manifest detection.withheld), %.1f s"),
                length(shares_withheld), length(pgs_all), length(detection_withheld), t_shares))

# ---- write the parts -----------------------------------------------------------------
dir.create(file.path(out_dir, "parts"), showWarnings = FALSE, recursive = TRUE)
# The old manifest and index go with the old parts: a run stopped from here on (a killed worker) then
# leaves a folder with no manifest.json -- "run make_screen_data.R" -- not one whose parts are missing.
unlink(file.path(out_dir, c("manifest.json", "index.tsv")))
unlink(list.files(file.path(out_dir, "parts"), full.names = TRUE))
chunks <- split(pgs_all, ceiling(seq_along(pgs_all) / PART_SIZE))
index <- list(); parts <- list()
t_records <- 0; t_write <- 0     # seconds: building + serialising the records, then writing the parts
for (ci in seq_along(chunks)) {
  pg_c <- chunks[[ci]]
  # each worker builds AND serialises its records (the JSON is most of the work), so --cpus speeds both
  one <- function(pg) as.character(to_json(build_record(pg, match(pg, rownames(y$E)), ctx)))
  t0 <- Sys.time()
  recs <- par_lapply(pg_c, one)
  if (!is.na(drop_record)) { k <- match(pgs_all[drop_record], pg_c); if (!is.na(k)) recs[k] <- list(NULL) }
  t_records <- t_records + as.numeric(difftime(Sys.time(), t0, units = "secs"))
  # A worker that threw gives a try-error; one that was killed (out of memory) gives NULL, which unlist()
  # below would drop, filing the next records under the wrong keys. Either stops the run.
  bad <- vapply(recs, function(r) is.null(r) || inherits(r, "try-error"), NA)
  if (any(bad)) {
    b <- which(bad)[1]
    stop("building the record of ", pg_c[b], " failed: ", if (is.null(recs[[b]]))
      "its worker returned nothing (killed, e.g. out of memory; rerun with fewer --cpus or more memory)"
      else as.character(recs[[b]]))
  }
  stopifnot(length(recs) == length(pg_c), length(unlist(recs)) == length(pg_c))
  t0 <- Sys.time()
  pname <- sprintf("part-%04d.json.gz", ci)
  con <- gzcon(file(file.path(out_dir, "parts", pname), "wb"))
  writeLines(paste0("{", paste0(vapply(pg_c, function(pg) as.character(to_json(pg)), ""), ":", unlist(recs),
                                collapse = ","), "}"), con, useBytes = TRUE)
  close(con)
  parts[[file.path("parts", pname)]] <- file_sha256(file.path(out_dir, "parts", pname))
  t_write <- t_write + as.numeric(difftime(Sys.time(), t0, units = "secs"))
  for (pg in pg_c) {
    k <- match(pg, rownames(y$E))
    index[[length(index) + 1L]] <- data.frame(Protein.Group = pg, Genes = y$genes$Genes[k],
      Protein.Names = y$genes$Protein.Names[k], part = pname, n_precursors = length(raw_idx[[pg]]),
      stringsAsFactors = FALSE)
  }
}
index <- do.call(rbind, index)
# One P.Value / adj.P.Val / logFC column per contrast so a list ("the top N changed proteins") can be ranked
# without opening the parts. NA (written empty) where block_de calls the protein untested.
for (co in contrasts) {
  tb <- de_tabs[[co$name]]; k <- match(index$Protein.Group, rownames(y$E))
  ok <- is.finite(tb$logFC[k]) & is.finite(tb$t[k])
  for (col in c("P.Value", "adj.P.Val", "logFC"))
    index[[paste0(col, "__", co$name)]] <- ifelse(ok, as.numeric(tb[[col]][k]), NA_real_)
}
utils::write.table(index, file.path(out_dir, "index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
parts[["index.tsv"]] <- file_sha256(file.path(out_dir, "index.tsv"))   # the index is read too: hashed like a part

# ---- manifest (last) --------------------------------------------------------------------
for (b in c("xic", "library_uniqueness", "check_routes", "fran", "report_text")) note_skip(b, NO_LATER)
srcs <- unique(sources)
manifest <- list(
  schema = SCHEMA, written_by = WRITTEN_BY, written_at = iso_now(), supported = TRUE,
  unsupported_reason = NULL, partial = !is.null(only_pgs),
  analysis = analysis,
  # where the files the page names are, relative to output/ (the report: the one the DE read)
  files = list(session_rds = "tables/DE-LIMP_session.rds", report = report_rel, stats = stats_rel),
  samples = arr(samples), groups = arr(groups), contrasts = contrasts,
  dpc = arr(as.numeric(sess$dpc_fit$dpc[c("beta0", "beta1")])),
  prior = list(prior.mean = y$prior.mean, prior.sd = y$prior.sd, prior.logFC = y$prior.logFC),
  n_precursors = nrow(raw$E), run_stats = run_stats,
  sources = setNames(lapply(srcs, function(r) file_sha256(file.path(output, r))), srcs),
  parts = parts, index = "index.tsv", skipped = skipped,
  shares = c(list(method = "derivative",
                  derivative = paste("exact derivative of limpa's protein values at its optimum, from limpa's",
                                     "own Hessian (implicit-function theorem): the limit of nudging a value"),
                  limpa_version = as.character(utils::packageVersion("limpa")),
                  tolerance = SHARE_TOL, settle_tolerance = SETTLE_TOL, nudge_check = nudge_check,
                  withheld = shares_withheld),
             if (!is.na(sigma_scale)) list(test_hook_sigma_scale = sigma_scale),
             if (!is.na(deriv_scale)) list(test_hook_derivative_scale = deriv_scale)),
  detection = list(
    method = paste("limpa's chance of finding each precursor in each sample (design section 5.7): pf_fit at",
                   "the level limpa's fit gives it, over limpa's normal Gauss quadrature with the protein's",
                   "sigma; pf_obs, the detection curve at a found value"),
    fit = paste("the shares' refit where BFGS stopped (par.bgfs), the point that reproduced the delivered",
                "value: dpcQuant runs no Newton polish"),
    dpc = arr(as.numeric(y$dpc)), quadrature_points = 16L, objective_tolerance = DETECT_TOL,
    objective_max_diff = objective_max_diff, withheld = detection_withheld),
  elapsed_s = list(total = round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1),
                   read = round(t_read, 1), shares = round(t_shares, 1), records = round(t_records, 1),
                   write = round(t_write, 1)),
  compute = list(de = prov$compute$cpu_family, here = here_family))
manifest["unsupported_reason"] <- list(NULL)
write_manifest(manifest)
write_status("built")
message(sprintf("[make_screen_data] %d proteins in %d part(s) -> %s", nrow(index), length(chunks), out_dir))
