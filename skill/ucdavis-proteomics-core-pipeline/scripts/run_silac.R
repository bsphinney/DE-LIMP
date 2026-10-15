#!/usr/bin/env Rscript
# =============================================================================
# run_silac.R -- SILAC (DIA-NN plexDIA) ratio analysis for the proteomics-pipeline skill.
#
# A SILAC run holds several samples, one per channel, so the comparison is a RATIO inside each
# run, not a difference between runs: run_de.R's group-vs-group model does not apply. This is the
# analysis a Core staff member had to write by hand for a triple-SILAC experiment (2026-10-14),
# made the skill's own: the statistics are silac_ratios.R (see there for the method), and this
# script reads the report and the channel design, writes the tables, figures, Methods and a
# de_provenance.json that describes itself (pipeline_id silac_ratio), so the report, the Methods,
# the audit, the delivery gate and the reproducibility bundle read it like any other DE record.
#
# Usage:
#   Rscript run_silac.R --input <search>/report.parquet --design input/silac_design.csv \
#       --outdir <S>/output/tables [--figures-dir <S>/output/figures] \
#       [--params <S>/input/wf/params.cfg] [--fasta <S>/input/search.fasta] \
#       [--comparisons "KO-WT,DKO-WT"] [--min-replicates 3] [--q-cutoff 0.01] \
#       [--pg-q-cutoff 0.05] [--channel-q-cutoff 0.01] [--adjp 0.05] [--logfc 1] \
#       [--scatter "KO-WT:DKO-KO"]
#
# The design (silac_design.csv; silac.py validates it -- the one definition of its rules):
#   File.Name,Replicate,Channel,Group      one row per run and channel that holds a sample
# Outputs (in --outdir): DE_silac_<comparison>.csv, SILAC_channel_specific_detections.csv,
#   SILAC_channel_summary.csv, SILAC_ratio_offsets.csv, SILAC_label_incorporation.csv (label-check
#   runs only), SILAC_sample_matrix.csv + SILAC_samples.csv (per replicate x channel, for
#   sample_quality.py), de_provenance.json, methods.txt, sessionInfo.txt, reproducibility_log.R;
#   figures and figures.json in --figures-dir.
# =============================================================================

for (p in c("limma", "data.table", "arrow", "jsonlite"))
  if (!requireNamespace(p, quietly = TRUE))
    stop(sprintf("run_silac.R needs the R package %s (the skill's environment has it: setup.sh)", p))
suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) == 0) return(default)
  if (i == length(args)) return(TRUE)
  val <- args[i + 1]
  if (startsWith(val, "--")) return(TRUE)
  val
}
given <- function(flag) flag %in% args
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

input       <- getarg("--input")
design_path <- getarg("--design")
outdir      <- getarg("--outdir", "silac_results")
figdir      <- getarg("--figures-dir",
                      if (basename(normalizePath(outdir, mustWork = FALSE)) == "tables")
                        file.path(dirname(normalizePath(outdir, mustWork = FALSE)), "figures")
                      else file.path(outdir, "figures"))
params_path <- getarg("--params")
fasta_path  <- getarg("--fasta")
comparisons <- getarg("--comparisons")
scatter_arg <- getarg("--scatter")
num_or <- function(flag, default) {
  v <- suppressWarnings(as.numeric(getarg(flag, as.character(default))))
  if (is.na(v)) stop(sprintf("%s must be a number", flag))
  v
}
min_reps  <- as.integer(num_or("--min-replicates", 3))
q_cut     <- num_or("--q-cutoff", 0.01)
pg_q_cut  <- num_or("--pg-q-cutoff", 0.05)
ch_q_cut  <- num_or("--channel-q-cutoff", 0.01)
adjp_thr  <- num_or("--adjp", 0.05)
logfc_ref <- num_or("--logfc", 1)
if (is.null(input) || is.null(design_path))
  stop("Required: --input <report.parquet> and --design <silac_design.csv>")
if (!file.exists(input)) stop("--input ", input, " does not exist")
if (!file.exists(design_path)) stop("--design ", design_path, " does not exist")
if (min_reps < 2) stop("--min-replicates must be at least 2: a one-sample test needs a variance")

.script_dir <- local({
  f <- grep("^--file=", commandArgs(), value = TRUE)[1]
  if (is.na(f)) getwd() else dirname(normalizePath(sub("^--file=", "", f), mustWork = FALSE))
})
.sibling <- function(nm) {
  for (p in c(file.path(.script_dir, nm), nm)) if (file.exists(p)) return(p)
  stop(nm, " not found next to run_silac.R")
}
source(.sibling("silac_ratios.R"))      # the statistics
source(.sibling("contaminants.R"))      # the contaminant rule (one definition) + python_json()
source(.sibling("skill_version.R"))
version_label <- paste("run_silac.R, ucdavis-proteomics-core-pipeline",
                       skill_label(skill_version(.script_dir)))

# The method's settings, each with where it came from (DE-LIMP rule 2)
setting <- function(flag, value, why)
  list(value = value, source = if (given(flag)) sprintf("given (%s)", flag) else
         sprintf("default of run_silac.R: %s", why))
settings <- list(
  min_replicates = setting("--min-replicates", min_reps,
                           "a ratio in at least 3 replicates, as validated on 4-replicate triple SILAC"),
  q_cutoff = setting("--q-cutoff", q_cut, "precursor Q.Value, DIA-NN's recommended 1%"),
  pg_q_cutoff = setting("--pg-q-cutoff", pg_q_cut,
                        "channel-specific PG.Q.Value; DIA-NN's README: 0.05 is typically sufficient"),
  channel_q_cutoff = setting("--channel-q-cutoff", ch_q_cut,
                             "Channel.Q.Value; DIA-NN's README recommends 0.01-0.5 for plexDIA"),
  adjp = setting("--adjp", adjp_thr, "Benjamini-Hochberg adjusted p below 0.05"),
  logfc = setting("--logfc", logfc_ref, "reference line on the volcano only, never a filter"))

# ---- the labels: the params sidecar (estimate_params.py --labels), else the cfg itself --------
labels <- NULL; labels_source <- "not given (--params): the channels are read from the report"
if (!is.null(params_path)) {
  side <- paste0(params_path, ".rationale.json")
  if (file.exists(side)) {
    s <- jsonlite::fromJSON(side, simplifyVector = FALSE)
    if (!is.null(s$labels)) { labels <- s$labels; labels_source <- basename(side) }
  }
  if (is.null(labels)) {
    r <- python_json("silac.py", "check-cfg", shQuote(params_path), "read the cfg's channels")
    if (!is.null(r$res$labels)) {
      labels <- jsonlite::fromJSON(jsonlite::toJSON(r$res$labels, auto_unbox = TRUE),
                                   simplifyVector = FALSE)
      labels_source <- paste(basename(params_path), "(its --channels; no --labels record)")
    } else labels_source <- paste("not recorded:", r$why %||% "the cfg has no --channels")
  }
}
cn <- labels$channel_norm
if (!is.null(cn$value) && identical(cn$value, "run"))
  stop("The search used --channel-run-norm: the channels are old and new protein of one culture ",
       "(pulsed SILAC / turnover). run_silac.R compares independent samples per channel; a ",
       "turnover analysis asks another question (references/silac.md). Nothing was written.")

# ---- the report ---------------------------------------------------------------------------------
need <- c("Run", "Channel", "Protein.Group", "Protein.Ids", "Precursor.Id", "Precursor.Quantity",
          "Q.Value", "PG.Q.Value", "Channel.Q.Value")
optional <- c("Protein.Names", "Genes")
have <- names(arrow::open_dataset(input)$schema)
miss <- setdiff(need, have)
if (length(miss))
  stop(input, " lacks ", paste(miss, collapse = ", "),
       if ("Channel" %in% miss) ": it is not a plexDIA (--channels) report" else "")
d_raw <- as.data.table(arrow::read_parquet(input, col_select = intersect(c(need, optional), have)))
for (o in optional) if (!o %in% names(d_raw)) d_raw[, (o) := NA_character_]
chan_order <- vapply(labels$channels %||% list(), function(c) c$name, "")
if (!length(chan_order)) {
  ch <- setdiff(unique(d_raw$Channel), "")
  chan_order <- ch[order(match(ch, c("L", "M", "H")), ch)]
}
light <- labels$light_channel %||% chan_order[1]

# ---- the design (silac.py: the one definition of its rules) -------------------------------------
runs_file <- tempfile(fileext = ".txt")
writeLines(unique(d_raw$Run), runs_file)
r <- python_json("silac.py", "design",
                 c("--validate", shQuote(design_path), "--runs-file", shQuote(runs_file),
                   "--channels", shQuote(paste(chan_order, collapse = ",")), "--json",
                   if (!is.null(comparisons)) c("--comparisons", shQuote(comparisons))),
                 "validate the SILAC design")
if (is.null(r$res)) stop("the design could not be checked: ", r$why)
if (!isTRUE(r$res$ok)) stop(r$res$error)
des <- r$res
asg <- as.data.table(des$assignments)
ratio_asg <- asg[role == "ratio"]
cmp <- as.data.table(des$comparisons)
groups <- des$groups
reps_all <- unique(ratio_asg$Replicate)

# ---- descriptions: the FASTA headers (UniProt and NCBI) --------------------------------------
fasta_desc <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  h <- sub("^>", "", grep("^>", readLines(path, warn = FALSE), value = TRUE))
  first <- sub("\\s.*$", "", h); rest <- trimws(sub("^\\S+\\s*", "", h))
  acc <- ifelse(grepl("^[a-z]{2}\\|[^|]+\\|", first), sub("^[a-z]{2}\\|([^|]+)\\|.*$", "\\1", first), first)
  gene <- ifelse(grepl("\\bGN=\\S+", rest), sub("^.*\\bGN=(\\S+).*$", "\\1", rest), NA_character_)
  desc <- sub("\\s+OS=.*$", "", rest)                     # UniProt: up to OS=
  desc <- sub("\\s*\\[[^]]*\\]\\s*$", "", desc)            # NCBI: the trailing [organism]
  desc <- sub("^MULTISPECIES:\\s*", "", desc)
  data.table(acc = c(acc, first), desc = c(desc, desc), gene = c(gene, gene))[!duplicated(acc)]
}
fd <- fasta_desc(fasta_path)
first_acc <- function(pg) sub(";.*$", "", pg)
annot <- unique(d_raw[, .(Protein.Group, Protein.Names, Genes)], by = "Protein.Group")
annot[, Description := if (is.null(fd)) NA_character_ else fd$desc[match(first_acc(Protein.Group), fd$acc)]]
annot[, Label := {
  g <- if (is.null(fd)) rep(NA_character_, .N) else fd$gene[match(first_acc(Protein.Group), fd$acc)]
  ifelse(!is.na(g) & nzchar(g), g,
         ifelse(!is.na(Description) & nzchar(Description), substr(sub(",.*$", "", Description), 1, 32),
                ifelse(!is.na(Genes) & nzchar(Genes), sub(";.*$", "", Genes), first_acc(Protein.Group))))
}]

# ---- the analysis -------------------------------------------------------------------------------
is_cont <- is_contaminant(d_raw$Protein.Ids)
prep <- silac_prepare(d_raw, ratio_asg, q_cut, pg_q_cut, ch_q_cut, is_cont)
dd <- prep$d
if (!nrow(dd)) stop("no precursor passed the filters in the runs the design holds")
w <- silac_wide(dd)
loaded <- lapply(setNames(groups, groups), function(g) unique(ratio_asg[Group == g]$Replicate))
results <- list(); offsets <- list(); tables_rec <- list()
for (i in seq_len(nrow(cmp))) {
  a <- cmp$numerator[i]; b <- cmp$denominator[i]; ct <- cmp$contrast[i]
  reps <- intersect(loaded[[a]], loaded[[b]])
  res <- silac_compare(w, a, b, reps, min_reps)
  tab <- merge(annot[, .(Protein.Group, Protein.Names, Genes, Description)], res$table,
               by = "Protein.Group", all.y = TRUE)
  det <- dd[ok == TRUE & Group %in% c(a, b), .(n = uniqueN(Replicate)), by = .(Protein.Group, Group)]
  for (g in c(a, b)) {
    k <- det[Group == g][match(tab$Protein.Group, Protein.Group)]$n
    tab[, (paste0("Detected_", g)) := sprintf("%d/%d", ifelse(is.na(k), 0L, k), length(loaded[[g]]))]
  }
  first <- c("Protein.Group", "Genes", "Protein.Names", "Description", "logFC", "t", "P.Value",
             "adj.P.Val", "n_replicates", paste0("Detected_", c(a, b)), "log2_ratio_uncentred_mean",
             "Evidence")
  setcolorder(tab, c(first, setdiff(names(tab), first)))
  setorderv(tab, c("P.Value", "Protein.Group"), na.last = TRUE)
  tab[, comparison := sprintf("%s vs %s (log2 %s/%s)", a, b, a, b)]
  results[[ct]] <- tab
  offsets[[ct]] <- res$offsets[, .(comparison = ct, Replicate, offset_log2 = offset,
                                   fold = 2^offset)]
  sig <- tab[!is.na(adj.P.Val) & adj.P.Val < adjp_thr]
  tables_rec[[ct]] <- list(
    file = sprintf("DE_silac_%s.csv", make.names(ct)), numerator = a, denominator = b,
    replicates = as.list(reps), n_replicates_available = length(reps),
    n_proteins = nrow(tab), n_tested = res$n_tested, n_significant = nrow(sig),
    n_up = sum(sig$logFC > 0), n_down = sum(sig$logFC < 0),
    df_prior = res$fit$df_prior, s2_prior = res$fit$s2_prior,
    detectable_log2 = res$fit$detectable, detectable_log2_alone = res$fit$detectable_alone,
    offset_log2 = setNames(as.list(res$offsets$offset), res$offsets$Replicate))
  message(sprintf("[run_silac] %s: %d proteins with a ratio; %d tested (ratio in >= %d of %d replicates); %d significant at adj.P < %g (%d up, %d down)",
                  ct, nrow(tab), res$n_tested, min_reps, length(reps), nrow(sig), adjp_thr,
                  sum(sig$logFC > 0), sum(sig$logFC < 0)))
}
detect <- silac_detection(dd, loaded, min_reps)
detect <- merge(annot[, .(Protein.Group, Genes, Protein.Names, Description)], detect, by = "Protein.Group")
specific <- detect[!is.na(pattern)]
setorderv(specific, c("pattern", "Protein.Group"))
chsum <- silac_channel_summary(dd)
chsum <- chsum[order(Replicate, match(Channel, chan_order))]
checks_runs <- asg[role == "label_check"]
incorp <- if (nrow(checks_runs)) {
  dchk <- d_raw[Run %in% checks_runs$Run][!is_contaminant(Protein.Ids)]
  silac_incorporation(dchk, setNames(as.list(checks_runs$Channel), checks_runs$Run), light,
                      q_cut, pg_q_cut, ch_q_cut)
} else NULL
lg <- labels$label_group %||% "SILAC"
doubled <- mean(grepl(sprintf("\\(%s\\)\\(%s\\)", lg, lg), d_raw$Precursor.Id, fixed = FALSE))

# ---- tables ------------------------------------------------------------------------------------
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
.prev <- file.path(outdir, "de_provenance.json")
if (file.exists(.prev)) {
  o <- tryCatch(jsonlite::fromJSON(.prev, simplifyVector = FALSE), error = function(e) list())
  if (!is.null(o$pipeline_id) && !identical(o$pipeline_id, "silac_ratio"))
    stop(outdir, " holds another analysis (", o$pipeline_id, "): write the SILAC analysis to ",
         "its own --outdir, or move that one first. Nothing was written.")
}
for (f in setdiff(list.files(outdir, pattern = "^DE_silac_.*\\.csv$"),
                  vapply(tables_rec, function(t) t$file, ""))) {
  unlink(file.path(outdir, f))
  message("[run_silac] removed ", f, ": a comparison of an earlier run, not of this one")
}
for (ct in names(results)) fwrite(results[[ct]], file.path(outdir, tables_rec[[ct]]$file))
fwrite(specific, file.path(outdir, "SILAC_channel_specific_detections.csv"))
fwrite(chsum, file.path(outdir, "SILAC_channel_summary.csv"))
fwrite(rbindlist(offsets), file.path(outdir, "SILAC_ratio_offsets.csv"))
# Per sample (replicate x channel): log2 of the summed channel-passing raw precursor quantity,
# median-centred within the sample -- for the sample-quality panels (sample_quality.py) and the
# audit's depth checks, which compare samples' composition. Not what the ratios test.
em <- dd[ok == TRUE, .(q = sum(Precursor.Quantity)), by = .(Protein.Group, Replicate, Group, Channel, Run)]
em[, Sample := paste(Replicate, Group, sep = "_")]
em[, v := log2(q) - stats::median(log2(q)), by = Sample]
samples_tab <- unique(em[, .(File.Name = Sample, Group, Replicate, Channel, Run)])
samples_tab <- samples_tab[order(Replicate, match(Channel, chan_order))]
EM <- dcast(em, Protein.Group ~ Sample, value.var = "v")
setcolorder(EM, c("Protein.Group", samples_tab$File.Name))
EM <- merge(annot[, .(Protein.Group, Genes, Protein.Names)], EM, by = "Protein.Group", all.y = TRUE)
fwrite(EM, file.path(outdir, "SILAC_sample_matrix.csv"))
fwrite(samples_tab, file.path(outdir, "SILAC_samples.csv"))
if (length(incorp)) fwrite(rbindlist(incorp), file.path(outdir, "SILAC_label_incorporation.csv"))
if (!length(incorp) && file.exists(file.path(outdir, "SILAC_label_incorporation.csv")))
  unlink(file.path(outdir, "SILAC_label_incorporation.csv"))

# ---- checks (audit_results.py quotes them: {check, status, message, detail}) -----------------
CHANNEL_SHARE_BAND <- c(0.5, 1.5)   # x the equal share (1/plex) before a channel is "unbalanced"
INCORPORATION_MIN <- 95             # % labelled; SILAC practice asks for > 95%
checks <- list()
chk <- function(check, status, message, detail = NULL)
  checks[[length(checks) + 1]] <<- list(check = check, status = status, message = message,
                                        detail = detail %||% list())
eq <- 100 / max(1, length(unique(ratio_asg$Channel)))
share_rng <- chsum[, .(lo = min(share_percent), hi = max(share_percent)), by = Group]
off_band <- share_rng[hi > CHANNEL_SHARE_BAND[2] * eq | lo < CHANNEL_SHARE_BAND[1] * eq]
if (nrow(off_band)) {
  chk("silac_channel_balance", "WARN", sprintf(paste(
    "The channels are not balanced: %s of the identified precursor signal per replicate, where an",
    "equal mix gives %.0f%%. The uncentred ratios (uncentred_* columns, SILAC_ratio_offsets.csv)",
    "carry this offset; the tested, median-centred ratios remove it, so they read 'relative to the",
    "typical protein'. The strongest decreases in a small channel can fall below detection there and",
    "then appear only as channel-specific detections, never as ratios. Confirm how the samples were",
    "mixed (equal volume is not equal protein)."),
    paste(sprintf("%s %.0f-%.0f%%", share_rng$Group, share_rng$lo, share_rng$hi), collapse = ", "), eq),
    list(equal_share_percent = eq, band = CHANNEL_SHARE_BAND,
         per_group = setNames(lapply(seq_len(nrow(share_rng)), function(i)
           list(min = share_rng$lo[i], max = share_rng$hi[i])), share_rng$Group)))
} else chk("silac_channel_balance", "PASS", sprintf(
  "Every channel holds %.0f-%.0f%% of the signal per replicate (an equal mix gives %.0f%%).",
  min(chsum$share_percent), max(chsum$share_percent), eq))
for (ct in names(tables_rec)) {
  n <- tables_rec[[ct]]$n_replicates_available
  if (n < min_reps)
    chk("silac_replicates", "WARN", sprintf(paste(
      "%s: only %d replicate(s) hold both channels, below the %d a ratio needs to be tested; no",
      "protein is tested in this comparison."), ct, n, min_reps))
  else if (n == min_reps)
    chk("silac_replicates", "WARN", sprintf(paste(
      "%s: %d replicates, the minimum: a protein is tested only when every replicate has its ratio,",
      "and the variance rests on %d values per protein (moderated by limma)."), ct, n, n))
  else chk("silac_replicates", "PASS", sprintf("%s: %d replicates; a protein is tested with a ratio in at least %d.",
                                               ct, n, min_reps))
}
if (length(incorp)) {
  low <- Filter(function(x) isTRUE(x$median_percent < INCORPORATION_MIN), incorp)
  txt <- paste(vapply(incorp, function(x) sprintf("%s %s %.1f%% (signal %.1f%%; %d precursors, %d with no light row)",
                                                  x$run, x$channel, x$median_percent, x$signal_percent,
                                                  x$precursors, x$precursors_without_light_row), ""),
               collapse = "; ")
  chk("silac_label_incorporation", if (length(low)) "WARN" else "PASS", paste0(
    "Label incorporation from the label-check runs (a lower bound): ", txt, ".",
    if (length(low)) sprintf(paste(" Below %g%%: unlabelled protein in the labelled cells adds light",
                                   "signal and pulls ratios towards the light channel."), INCORPORATION_MIN) else ""))
} else chk("silac_label_incorporation", "INFO", paste(
  "Label incorporation was not measured: every run is a mix that holds the light channel by design,",
  "and plexDIA gives each peptide one channel, so a partly labelled peptide is not searched. Measure",
  "it on a labelled-only injection of each labelled culture (a design row with that channel alone), or",
  "from a search with the labels as variable modifications (two-site method, references/silac.md)."))
if (doubled > 0)
  chk("silac_label_declaration", "WARN", sprintf(paste(
    "%.1f%% of the precursor ids carry the label twice (%s(%s)(%s)): the search declared it twice",
    "(--lib-fixed-mod on a library predicted from the FASTA), which DIA-NN reports as not supported.",
    "The channels were still assigned, but for a deliverable re-search with estimate_params.py",
    "--labels (references/silac.md)."), 100 * doubled, "K", lg, lg))
if (!is.null(cn) && identical(cn$given, FALSE))
  chk("silac_channel_norm", "INFO", paste(
    "The search's channel normalisation (", cn$flag %||% "?", ") is a DEFAULT, not user-confirmed.",
    "These ratios do not use DIA-NN's normalised quantities, so it does not change them; it changes",
    "the report's normalised columns (Precursor.Normalised, PG.MaxLFQ)."))
chk("silac_channel_specific", "INFO", sprintf(paste(
  "%d protein groups are present (>= %d replicates) in some channels and absent from others: no",
  "ratio, so not tested; listed in SILAC_channel_specific_detections.csv (%s). Check it for whole",
  "pathways before calling it a loading or detection-limit effect."), nrow(specific), min_reps,
  if (nrow(specific)) paste(sprintf("%s: %d", names(table(specific$pattern)), as.integer(table(specific$pattern))),
                            collapse = "; ") else "none"))

# ---- figures ------------------------------------------------------------------------------------
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
figs <- list(); failed <- list()
fig_ok <- function(file, type, caption) figs[[length(figs) + 1]] <<- list(file = file, type = type, caption = caption)
fig_fail <- function(file, type, reason) failed[[length(failed) + 1]] <<- list(file = file, type = type, reason = reason)
draw <- function(file, type, expr) tryCatch(expr, error = function(e) {
  fig_fail(file, type, conditionMessage(e))
  message("[run_silac] figure ", file, " not drawn: ", conditionMessage(e))
})
has_gg <- requireNamespace("ggplot2", quietly = TRUE)
has_repel <- requireNamespace("ggrepel", quietly = TRUE)
ink <- "#22302b"; muted <- "#5b6661"; up_c <- "#b0472a"; dn_c <- "#2f6db3"; ns_c <- "#b9c2be"
if (!has_gg) {
  for (ct in names(results)) fig_fail(sprintf("volcano_%s.png", make.names(ct)), "volcano",
                                      "ggplot2 is not installed")
} else {
  library(ggplot2)
  th <- theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", colour = ink), plot.subtitle = element_text(colour = muted),
    axis.text = element_text(colour = muted), legend.position = "top")
  all_lfc <- unlist(lapply(results, function(x) x$logFC))
  xl <- if (any(is.finite(all_lfc))) range(all_lfc, na.rm = TRUE) * 1.05 else c(-1, 1)
  for (ct in names(results)) {
    f <- sprintf("volcano_%s.png", make.names(ct)); a <- tables_rec[[ct]]
    draw(f, "volcano", {
      x <- merge(results[[ct]][!is.na(P.Value)], annot[, .(Protein.Group, Label)], by = "Protein.Group")
      if (!nrow(x)) stop("no protein was tested in this comparison")
      x[, cls := ifelse(adj.P.Val < adjp_thr & logFC > 0, "up", ifelse(adj.P.Val < adjp_thr, "down", "ns"))]
      lab <- if (a$n_significant) x[order(P.Value)][seq_len(min(12, .N))] else x[0]
      g <- ggplot(x, aes(logFC, -log10(P.Value))) +
        geom_point(aes(colour = cls), size = 1.3, alpha = 0.7) +
        geom_vline(xintercept = 0, colour = muted, linewidth = 0.3) +
        geom_vline(xintercept = c(-logfc_ref, logfc_ref), colour = muted, linewidth = 0.3, linetype = "dotted") +
        scale_colour_manual(values = c(up = up_c, down = dn_c, ns = ns_c), name = NULL,
                            labels = c(up = paste("higher in", a$numerator), down = paste("lower in", a$numerator),
                                       ns = "not significant")) +
        scale_x_continuous(limits = xl) +
        labs(x = sprintf("median-centred log2 %s/%s", a$numerator, a$denominator), y = "-log10 p-value",
             title = sprintf("%s vs %s", a$numerator, a$denominator),
             subtitle = sprintf("%s tested (ratio in >= %d of %d replicates); %s significant at adj. p < %g",
                                format(a$n_tested, big.mark = ","), min_reps, a$n_replicates_available,
                                format(a$n_significant, big.mark = ","), adjp_thr)) + th
      if (nrow(lab)) g <- g + (if (has_repel) ggrepel::geom_text_repel(data = lab, aes(label = Label), size = 2.7,
                                                                         colour = ink, max.overlaps = 20)
                               else geom_text(data = lab, aes(label = Label), size = 2.7, colour = ink,
                                              vjust = -0.8, check_overlap = TRUE))
      ggsave(file.path(figdir, f), g, width = 8, height = 6.2, dpi = 200, bg = "white")
      fig_ok(f, "volcano", sprintf(paste(
        "Volcano plot, %s vs %s: each point a protein group tested in at least %d of %d replicates;",
        "x the mean median-centred log2 ratio %s/%s, y the moderated one-sample t-test p-value.",
        "%d significant at Benjamini-Hochberg adj. p < %g (%d up, %d down); dotted lines mark |log2| = %g",
        "for reference only. All comparisons share the x-axis."), a$numerator, a$denominator, min_reps,
        a$n_replicates_available, a$numerator, a$denominator, a$n_significant, adjp_thr, a$n_up, a$n_down,
        logfc_ref))
    })
  }
  f <- "silac_channel_composition.png"
  draw(f, "silac_channel_composition", {
    x <- copy(chsum); x[, Group := factor(Group, levels = rev(groups))]
    g <- ggplot(x, aes(share_percent, Replicate, fill = Group)) + geom_col(width = 0.62, colour = "white") +
      geom_text(aes(label = sprintf("%.0f%%", share_percent)), position = position_stack(vjust = 0.5),
                colour = "white", size = 3.4, fontface = "bold") +
      geom_vline(xintercept = seq(eq, 100 - eq, by = eq), linetype = "dashed", colour = muted, linewidth = 0.4) +
      scale_fill_manual(values = setNames(c("#6b7a74", dn_c, up_c, "#8c6bb1")[seq_along(groups)], groups),
                        breaks = groups, name = NULL) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.03)), labels = function(v) paste0(v, "%")) +
      labs(x = "Share of the identified precursor signal", y = NULL, title = "Signal from each channel",
           subtitle = "Raw precursor quantities passing the channel q-value; dashed lines mark equal shares") + th
    ggsave(file.path(figdir, f), g, width = 8, height = 1.6 + 0.45 * length(unique(x$Replicate)), dpi = 200, bg = "white")
    fig_ok(f, "silac_channel_composition", sprintf(paste(
      "Each channel's share of the identified precursor signal (raw Precursor.Quantity of precursors",
      "passing channel q <= %g, contaminants removed), per replicate; dashed lines mark an equal mix."), ch_q_cut))
  })
  f <- "silac_channel_detections.png"
  draw(f, "silac_channel_detections", {
    x <- copy(chsum); x[, Group := factor(Group, levels = groups)]
    g <- ggplot(x, aes(Replicate, precursors, fill = Group)) + geom_col(position = position_dodge(0.8), width = 0.75) +
      scale_fill_manual(values = setNames(c("#6b7a74", dn_c, up_c, "#8c6bb1")[seq_along(groups)], groups), name = NULL) +
      scale_y_continuous(labels = function(v) format(v, big.mark = ",")) +
      labs(x = NULL, y = "Precursors", title = "Precursors identified in each channel",
           subtitle = sprintf("q <= %g, protein group q <= %g, channel q <= %g", q_cut, pg_q_cut, ch_q_cut)) + th
    ggsave(file.path(figdir, f), g, width = 8, height = 4.6, dpi = 200, bg = "white")
    fig_ok(f, "silac_channel_detections", sprintf(paste(
      "Precursors identified per replicate and channel (precursor q <= %g, protein group q <= %g,",
      "channel q <= %g, contaminants removed)."), q_cut, pg_q_cut, ch_q_cut))
  })
  sc <- if (!is.null(scatter_arg)) strsplit(scatter_arg, ":", fixed = TRUE)[[1]]
        else if (length(results) >= 3) names(results)[c(1, 3)]
        else if (length(results) == 2) names(results)[1:2] else NULL
  if (length(sc) == 2 && all(sc %in% names(results))) {
    f <- sprintf("silac_scatter_%s_vs_%s.png", make.names(sc[1]), make.names(sc[2]))
    draw(f, "silac_scatter", {
      x <- merge(results[[sc[1]]][!is.na(logFC), .(Protein.Group, x = logFC, px = adj.P.Val)],
                 results[[sc[2]]][!is.na(logFC), .(Protein.Group, y = logFC, py = adj.P.Val)], by = "Protein.Group")
      if (!nrow(x)) stop("no protein was tested in both comparisons")
      x[, cls := ifelse(py < adjp_thr, paste(sc[2], "significant"), ifelse(px < adjp_thr, paste(sc[1], "significant only"), "neither"))]
      g <- ggplot(x, aes(x, y)) + geom_hline(yintercept = 0, colour = muted, linewidth = 0.3) +
        geom_vline(xintercept = 0, colour = muted, linewidth = 0.3) +
        geom_point(aes(colour = cls), size = 1.3, alpha = 0.75) +
        scale_colour_manual(values = setNames(c(up_c, dn_c, ns_c), c(paste(sc[2], "significant"),
                                                                       paste(sc[1], "significant only"), "neither")), name = NULL) +
        labs(x = sprintf("centred log2 ratio, %s", gsub("-", "/", sc[1])), y = sprintf("centred log2 ratio, %s", gsub("-", "/", sc[2])),
             title = sprintf("%s against %s", gsub("-", "/", sc[2]), gsub("-", "/", sc[1])),
             subtitle = sprintf("%d proteins tested in both", nrow(x))) + th + guides(colour = guide_legend(nrow = 2))
      ggsave(file.path(figdir, f), g, width = 8, height = 6.6, dpi = 200, bg = "white")
      shared <- intersect(unlist(strsplit(sc[1], "-")), unlist(strsplit(sc[2], "-")))
      fig_ok(f, "silac_scatter", sprintf(paste(
        "Centred log2 ratios of the proteins tested in both comparisons: %s (x) against %s (y)%s."),
        gsub("-", "/", sc[1]), gsub("-", "/", sc[2]),
        if (length(shared)) sprintf(paste("; the two share the %s channel, so its measurement noise couples them",
                                          "and their agreement is not independent evidence"), shared[1]) else ""))
    })
  }
  f <- "qc_pvalue_panel.png"
  draw(f, "pvalue", {
    x <- rbindlist(lapply(names(results), function(ct) results[[ct]][!is.na(P.Value), .(comparison = ct, P.Value)]))
    if (!nrow(x)) stop("no p-values")
    g <- ggplot(x, aes(P.Value)) + geom_histogram(breaks = seq(0, 1, 0.05), fill = "#6b7a74", colour = "white") +
      facet_wrap(~comparison, scales = "free_y") + labs(x = "raw p-value", y = "proteins",
      title = "Raw p-value distributions", subtitle = "Flat with a peak near 0 when the test is calibrated") + th
    ggsave(file.path(figdir, f), g, width = 8, height = 3.4 + 0.6 * ceiling(length(results) / 3), dpi = 200, bg = "white")
    fig_ok(f, "pvalue", sprintf(paste("Raw p-value distributions for all %d comparison(s) (appendix; a",
                                      "calibration check, not a result)."), length(results)))
  })
}
writeLines(jsonlite::toJSON(list(figures = figs, failed = failed, written_by = "run_silac.R"),
                            auto_unbox = TRUE, pretty = TRUE, null = "null"),
           file.path(figdir, "figures.json"))

# ---- Methods, the brief, and the record --------------------------------------------------------
pkg <- function(p) tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
n_rep <- length(reps_all)
filters <- c(sprintf("Q.Value <= %g", q_cut), sprintf("PG.Q.Value <= %g (channel-specific)", pg_q_cut),
             sprintf("Channel.Q.Value <= %g in both channels of a ratio", ch_q_cut),
             "Precursor.Quantity > 0", "contaminant precursors removed (any accession tagged Cont_ / contam_)")
cmp_txt <- paste(vapply(names(tables_rec), function(ct) sprintf("%s/%s", tables_rec[[ct]]$numerator,
                                                                tables_rec[[ct]]$denominator), ""), collapse = ", ")
methods_paragraph <- paste0(
  sprintf("SILAC ratios were analysed with %s (limma %s; R %s). ", version_label, pkg("limma"),
          paste(R.version$major, R.version$minor, sep = ".")),
  sprintf(paste("From DIA-NN's main report, precursors were kept at Q.Value <= %g and channel-specific",
                "protein-group PG.Q.Value <= %g; precursors mapping to a contaminant entry were removed. A",
                "precursor ratio required both channels at Channel.Q.Value <= %g with a non-zero raw",
                "quantity (Precursor.Quantity; DIA-NN's channel-normalised quantities were not used). "),
          q_cut, pg_q_cut, ch_q_cut),
  sprintf(paste("For each comparison (%s), the log2 ratio of a protein group in each replicate was the median of",
                "its precursor log2 ratios, and each replicate's protein ratios were median-centred; the",
                "uncentred ratios, which carry the mixing offset, are reported alongside. "), cmp_txt),
  sprintf(paste("Protein groups with a ratio in at least %d of the %d replicates were tested against a log2",
                "ratio of 0 with a moderated one-sample t-test (limma lmFit with an intercept-only design and",
                "eBayes with robust = TRUE), and p-values were adjusted per comparison by the Benjamini-Hochberg",
                "method; protein groups with adj.P.Val < %g were called significant, with no fold-change filter. "),
          min_reps, n_rep, adjp_thr),
  sprintf(paste("Protein groups detected (Channel.Q.Value <= %g) in at least %d replicates of some channels and",
                "in no replicate of another have no ratio; they are listed separately as channel-specific",
                "detections. "), ch_q_cut, min_reps),
  if (length(incorp)) "Label incorporation was estimated from label-check runs as the labelled share of each precursor's signal (a lower bound). " else
    "Label incorporation was not measured from these data. ",
  "Citations: plexDIA, Derks et al., Nat Biotechnol 2023;41:50-59; limma, Ritchie et al., Nucleic Acids Res 2015;43:e47; ",
  "robust empirical Bayes, Phipson et al., Ann Appl Stat 2016;10:946-963; Benjamini and Hochberg, J R Stat Soc B 1995;57:289-300.")
off_txt <- paste(vapply(names(tables_rec), function(ct) {
  o <- unlist(tables_rec[[ct]]$offset_log2)
  sprintf("%s %+.2f to %+.2f", ct, min(o), max(o))
}, ""), collapse = "; ")
det_txt <- paste(vapply(names(tables_rec), function(ct) {
  v <- tables_rec[[ct]]$detectable_log2; v2 <- tables_rec[[ct]]$detectable_log2_alone
  if (is.finite(v)) sprintf("%s %.2f / %.2f", ct, v, v2) else sprintf("%s not estimated", ct)
}, ""), collapse = "; ")
brief_notes <- list(
  paste("logFC in the DE_silac tables is the MEDIAN-CENTRED log2 ratio (numerator/denominator): a protein's",
        "share relative to the typical protein of that replicate, so up and down counts come out close to",
        "even by construction. It is not 'amount per volume'."),
  sprintf(paste("The uncentred ratios (uncentred_<replicate>, log2_ratio_uncentred_mean) are amount per",
                "volume and carry each replicate's mixing offset; the median uncentred log2 ratio per replicate:",
                "%s (SILAC_ratio_offsets.csv)."), off_txt),
  paste("A ratio needs the precursor in both channels. A protein lost in the smaller channel has no ratio and",
        "appears only in SILAC_channel_specific_detections.csv, so the strongest decreases can be missing from",
        "the tables: read that list (whole pathways can hide there) before calling it a loading effect."),
  paste("Two comparisons that share a channel (B/A and C/A) are correlated by construction: their agreement is",
        "not independent evidence."),
  sprintf(paste("Not significant is not unchanged. From the median variance, a typical protein needs a centred",
                "|log2 ratio| of about the first value to reach p < 0.05 before the multiple-testing adjustment,",
                "and about the second to be significant after it when it is the only change: %s. A rough",
                "reading of the power, not a rule."), det_txt))
plain_language <- paste(
  "**SILAC ratios + limma**: each sample was grown with its own isotope label and the samples were mixed",
  "and measured together, so each protein is compared between labels inside the same run. The ratio is",
  "centred so that the typical protein reads 0, and a protein is called changed when its ratio differs",
  "from 0 consistently across the replicates (limma's moderated t-test, Benjamini-Hochberg adjusted).")
md_row <- function(...) paste0("| ", paste(..., sep = " | "), " |")
fmt_n <- function(x) format(x, big.mark = ",")
ch_md <- c(sprintf(paste("Each replicate is one run holding %s. Signal share = the channel's part of the",
                         "identified precursor signal (raw quantity, channel q <= %g); an equal mix gives %.0f%%."),
                   paste(sprintf("%s (%s)", unique(ratio_asg[, .(Channel, Group)])$Group,
                                 unique(ratio_asg[, .(Channel, Group)])$Channel), collapse = ", "), ch_q_cut, eq),
           "", md_row("Replicate", "Channel", "Group", "Precursors", "Protein groups", "Signal share"),
           "|---|---|---|---:|---:|---:|",
           vapply(seq_len(nrow(chsum)), function(i) md_row(chsum$Replicate[i], chsum$Channel[i], chsum$Group[i],
                                                           fmt_n(chsum$precursors[i]), fmt_n(chsum$protein_groups[i]),
                                                           sprintf("%.0f%%", chsum$share_percent[i])), ""),
           "", sprintf(paste("Median uncentred log2 ratio per replicate (the mixing offset the centred ratios",
                             "remove): %s. Label incorporation: %s"), off_txt,
                       if (length(incorp)) paste(vapply(incorp, function(x) sprintf("%s %.1f%%", x$run, x$median_percent), ""),
                                                 collapse = "; ")
                       else "not measured from these data (no labelled-only run)."))
sp_counts <- if (nrow(specific)) table(specific$pattern) else integer(0)
sp_md <- c(sprintf(paste("%s protein groups are present (channel q <= %g in at least %d replicates) in some",
                         "channels and absent from every replicate of another. They have no ratio, so they are not",
                         "in the tables above; the strongest losses in a small channel can only appear here. Full",
                         "list: `SILAC_channel_specific_detections.csv`."), fmt_n(nrow(specific)), ch_q_cut, min_reps),
           "", md_row("Pattern", "Protein groups"), "|---|---:|",
           vapply(names(sp_counts), function(k) md_row(k, fmt_n(as.integer(sp_counts[[k]]))), ""))
if (nrow(specific)) {
  top_sp <- merge(specific, annot[, .(Protein.Group, Label)], by = "Protein.Group")
  top_sp[, total := rowSums(.SD), .SDcols = paste0(groups, "_replicates")]
  setorderv(top_sp, c("pattern", "total", "Protein.Group"), c(1, -1, 1))
  top_sp <- top_sp[, head(.SD, 8), by = pattern]
  sp_md <- c(sp_md, "", sprintf("Up to 8 per pattern (replicates detected in %s):", paste(groups, collapse = " / ")), "",
             md_row("Protein", "Description", "Pattern", "Replicates"), "|---|---|---|---|",
             vapply(seq_len(nrow(top_sp)), function(i) md_row(
               top_sp$Label[i], gsub("\\|", "/", substr(ifelse(is.na(top_sp$Description[i]), "", top_sp$Description[i]), 1, 60)),
               top_sp$pattern[i],
               paste(vapply(groups, function(g) as.character(top_sp[[paste0(g, "_replicates")]][i]), ""), collapse = " / ")), ""))
}
report_sections <- list(list(title = "SILAC channels and mixing", markdown = paste(ch_md, collapse = "\n")),
                        list(title = "Channel-specific proteins", markdown = paste(sp_md, collapse = "\n")))
high_significance_note <- paste(
  "In a SILAC ratio test each protein is compared inside the same run, so run-to-run variation cancels and",
  "the test is powerful; between strains or treatments that differ broadly, a large share can be real. Check",
  "the channel assignment (a swapped label reverses a whole comparison) and the mixing (SILAC_channel_summary.csv);",
  "after median-centring, up and down counts close to even are expected.")
sig_counts <- lapply(tables_rec, function(t) list(up = t$n_up, down = t$n_down, total = t$n_significant))
input_sha256 <- tryCatch(if (requireNamespace("digest", quietly = TRUE))
  digest::digest(input, algo = "sha256", file = TRUE) else NA_character_, error = function(e) NA_character_)
replay <- c("--input", normalizePath(input), "--design", normalizePath(design_path),
            "--outdir", normalizePath(outdir), "--figures-dir", normalizePath(figdir),
            if (!is.null(params_path)) c("--params", normalizePath(params_path)),
            if (!is.null(fasta_path)) c("--fasta", normalizePath(fasta_path)),
            "--comparisons", paste(names(tables_rec), collapse = ","),
            "--min-replicates", min_reps, "--q-cutoff", q_cut, "--pg-q-cutoff", pg_q_cut,
            "--channel-q-cutoff", ch_q_cut, "--adjp", adjp_thr, "--logfc", logfc_ref)
prov <- list(
  pipeline_id = "silac_ratio", method = "silac",
  display_label = "SILAC ratios + limma (one-sample moderated t, robust eBayes)",
  methods_heading = "SILAC ratio analysis", methods_paragraph = methods_paragraph,
  plain_language = plain_language, brief_notes = brief_notes,
  rollup_method = "median of the precursor log2 ratios (raw Precursor.Quantity) per protein group and replicate",
  de_engine = "limma::lmFit (intercept only, per comparison) -> eBayes(robust = TRUE) -> BH",
  missing_policy = sprintf(paste("A ratio needs both channels at Channel.Q.Value <= %g; a protein group is tested",
                                 "with a ratio in at least %d replicates; one present in some channels and absent",
                                 "from others is listed in SILAC_channel_specific_detections.csv, not tested."),
                           ch_q_cut, min_reps),
  citation = paste("Derks J et al. (2023) Nat Biotechnol 41:50-59 (plexDIA); Ritchie ME et al. (2015) NAR 43:e47 (limma);",
                   "Phipson B et al. (2016) Ann Appl Stat 10:946-963 (robust eBayes)"),
  quantities = "raw",
  normalisation = paste("per-replicate median centring of the protein log2 ratios (uncentred ratios kept);",
                        "DIA-NN's channel-normalised quantities are not used (raw Precursor.Quantity)"),
  normalization_check = list(status = "not_applicable", reason = paste(
    "A SILAC ratio compares channels inside one run, so DIA-NN's between-run normalisation does not",
    "enter it; each replicate's ratios are median-centred instead, and the uncentred ratios are kept",
    "beside them. The channel-normalisation answer is the search's (references/silac.md)."),
    quantities_applied = "raw Precursor.Quantity, median-centred ratios"),
  set_tests = list(available = FALSE, suggested_for = list(), reason = paste(
    "run_sets.R reads run_de.R's group-vs-group fit; protein-set tests on SILAC ratios are not",
    "part of skill 2.13")),
  figures = list(written_by = "run_silac.R", figures_json = normalizePath(file.path(figdir, "figures.json"))),
  replay = list(script = "run_silac.R", args = as.list(as.character(replay))),
  contrasts = as.list(names(tables_rec)), adjp = adjp_thr, logfc = logfc_ref,
  logfc_role = "reference_line_only", significance_rule = "adj.P.Val < adjp (BH, per comparison); no fold-change filter",
  design = "an intercept-only (one-sample) model per comparison, one value per replicate (run)",
  n_samples = n_rep, n_replicates = n_rep,
  groups = setNames(lapply(loaded, length), names(loaded)),
  sample_matrix = list(file = "SILAC_sample_matrix.csv", samples = "SILAC_samples.csv", values = paste(
    "per sample (replicate x channel): log2 of the summed channel-passing raw precursor quantity,",
    "median-centred within the sample; for the sample-quality panels and depth checks, not what was tested")),
  report_sections = report_sections, high_significance_note = high_significance_note,
  de_columns = list(
    logFC = paste("the mean median-centred log2 ratio, numerator/denominator, over the replicates with a",
                  "ratio (limma's estimate): relative to the typical protein of each replicate; positive =",
                  "higher in the numerator (the first-named group); NA = not tested"),
    t = "moderated one-sample t-statistic (limma eBayes, robust)",
    P.Value = "raw p-value of the moderated one-sample t-test against a log2 ratio of 0; NA = not tested (Evidence says why)",
    adj.P.Val = "Benjamini-Hochberg adjusted p-value within the comparison -- THE significance column (NA = not tested)",
    Description = "the protein's description from the searched FASTA's header",
    n_replicates = "replicates (runs) in which the protein group had a ratio",
    `Detected_<group>` = paste("k/n: the replicates holding that group in which the protein group had at least",
                               "one precursor passing the channel q-value in that group's channel"),
    log2_ratio_uncentred_mean = paste("the mean UNCENTRED log2 ratio: amount per volume, carrying each",
                                      "replicate's mixing offset (SILAC_ratio_offsets.csv)"),
    Evidence = "tested, or why the protein group was not (a ratio in too few replicates)",
    `centred_<replicate>` = "that replicate's median-centred log2 ratio (median of the precursor log2 ratios, minus the replicate's median)",
    `uncentred_<replicate>` = "that replicate's log2 ratio before centring (the median of its precursor log2 ratios)",
    `n_precursors_<replicate>` = "precursors with both channels passing, behind that replicate's ratio",
    comparison = "the comparison, in words"),
  glance_tiles = list(list(n_rep, "replicates (runs)"), list(length(groups), "channels compared"),
                      list(nrow(specific), "channel-specific proteins")),
  q_columns = list("Q.Value", "PG.Q.Value", "Channel.Q.Value"), q_cutoffs = list(q_cut, pg_q_cut, ch_q_cut),
  filters_applied = as.list(filters), filter_counts = prep$counts,
  settings = settings,
  comparisons_source = des$comparisons_source,
  significant_per_contrast = sig_counts, de_tables = tables_rec,
  labels = labels, labels_source = labels_source,
  design_record = list(file = normalizePath(design_path), channels = as.list(chan_order),
                       label_swapped_groups = des$label_swapped_groups %||% list(),
                       label_check_runs = des$label_check_runs %||% list(),
                       unused_runs = des$unused_runs %||% list()),
  channel_summary = lapply(split(chsum, by = "Replicate"), function(x)
    setNames(lapply(seq_len(nrow(x)), function(i) list(channel = x$Channel[i], share_percent = x$share_percent[i],
                                                       precursors = x$precursors[i], protein_groups = x$protein_groups[i])),
             x$Group)),
  channel_specific = list(file = "SILAC_channel_specific_detections.csv", n = nrow(specific),
                          by_pattern = as.list(table(specific$pattern))),
  label_incorporation = if (length(incorp)) list(measured = TRUE, runs = unname(incorp),
                                                 basis = "labelled share H / (H + L) of the raw quantity; a lower bound")
                        else list(measured = FALSE, reason = paste(
                          "no label-check run (a run holding one labelled channel only): every run mixes the",
                          "light channel by design, and plexDIA gives each peptide one channel")),
  label_declaration = list(doubled_fraction = doubled,
                           note = if (doubled > 0) "the search declared the label twice (K(SILAC)(SILAC))" else "declared once"),
  checks = checks,
  input = normalizePath(input), input_sha256 = input_sha256, metadata = normalizePath(design_path),
  written = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  written_by = version_label, R_version = as.character(getRversion()),
  packages = list(limma = pkg("limma"), data.table = pkg("data.table"), arrow = pkg("arrow"),
                  ggplot2 = pkg("ggplot2"), jsonlite = pkg("jsonlite")),
  runtime = list(r_version = R.version.string, r_home = R.home(), lib_paths = as.list(.libPaths()),
                 conda_prefix = if (nzchar(Sys.getenv("CONDA_PREFIX"))) Sys.getenv("CONDA_PREFIX") else NULL,
                 session_info = "sessionInfo.txt", note = "recorded by run_silac.R in the process that ran it"))
writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA),
           file.path(outdir, "de_provenance.json"))
writeLines(c("SILAC ratio analysis (run_silac.R)", "", strwrap(methods_paragraph, 100), "",
             "Settings:", vapply(names(settings), function(k) sprintf("  %s = %s (%s)", k, settings[[k]]$value,
                                                                     settings[[k]]$source), "")),
           file.path(outdir, "methods.txt"))
con <- file(file.path(outdir, "sessionInfo.txt"), "w"); sink(con); print(utils::sessionInfo()); sink(); close(con)

# ---- reproducibility_log.R: the statistics as they ran, with every value written out ---------
rq <- function(x) paste0("'", gsub("'", "\\\\'", gsub("\\\\", "\\\\\\\\", as.character(x))), "'")
rvec <- function(x) paste0("c(", paste(vapply(x, rq, ""), collapse = ", "), ")")
core <- readLines(.sibling("silac_ratios.R"), warn = FALSE)
rl <- c("#!/usr/bin/env Rscript",
        sprintf("# The SILAC ratio analysis as it ran (%s, %s).", version_label, prov$written),
        "# Plain R: needs data.table, arrow and limma. Every value below is the one that ran.",
        "suppressMessages({library(data.table); library(limma)})",
        "", "# ---- the statistics (silac_ratios.R, verbatim) ----", core, "",
        "# ---- the inputs ----",
        sprintf("report <- %s", rq(normalizePath(input))),
        sprintf("design <- data.table(Run = %s,", rvec(ratio_asg$Run)),
        sprintf("                     Replicate = %s,", rvec(ratio_asg$Replicate)),
        sprintf("                     Channel = %s,", rvec(ratio_asg$Channel)),
        sprintf("                     Group = %s)", rvec(ratio_asg$Group)),
        sprintf("comparisons <- list(%s)", paste(vapply(names(tables_rec), function(ct)
          sprintf("%s = c(%s, %s)", rq(ct), rq(tables_rec[[ct]]$numerator), rq(tables_rec[[ct]]$denominator)), ""),
          collapse = ", ")),
        sprintf("q <- %s; pg_q <- %s; channel_q <- %s; min_reps <- %dL", q_cut, pg_q_cut, ch_q_cut, min_reps),
        sprintf("contaminant <- %s   # contaminants.R: any accession starting with a tag", rq(contaminant_regex())),
        "", "# ---- the run ----",
        sprintf("d <- as.data.table(arrow::read_parquet(report, col_select = %s))", rvec(need)),
        "prep <- silac_prepare(d, design, q, pg_q, channel_q, grepl(contaminant, d$Protein.Ids))",
        "w <- silac_wide(prep$d)",
        "loaded <- lapply(split(design$Replicate, design$Group), unique)",
        "for (ct in names(comparisons)) {",
        "  a <- comparisons[[ct]][1]; b <- comparisons[[ct]][2]",
        "  res <- silac_compare(w, a, b, intersect(loaded[[a]], loaded[[b]]), min_reps)",
        "  fwrite(res$table, sprintf('replay_DE_silac_%s.csv', make.names(ct)))",
        "  cat(ct, ': tested', res$n_tested, ' significant', sum(res$table$adj.P.Val < 0.05, na.rm = TRUE), '\\n')",
        "}")
writeLines(rl, file.path(outdir, "reproducibility_log.R"))

message(sprintf("[run_silac] done: %d comparison(s), %d channel-specific protein groups, %d figure(s)%s -> %s",
                length(results), nrow(specific), length(figs),
                if (length(failed)) sprintf(" (%d not drawn)", length(failed)) else "", normalizePath(outdir)))
