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
#       [--comparisons "KO-WT,DKO-WT"] [--min-replicates N] [--replicates biological] \
#       [--q-cutoff 0.01] [--pg-q-cutoff 0.05] [--channel-q-cutoff 0.01] [--run-level-only] \
#       [--adjp 0.05] [--logfc 1] [--scatter "KO-WT:DKO-WT"]
#
# --run-level-only drops the experiment-wide q-value columns (Global.*, Lib.*) and filters on the
#   run-specific Q.Value and channel-specific PG.Q.Value alone: the filter of the hand-made
#   analysis this method was validated against. Recorded; the default applies every column the
#   report has (diann_q_columns.R, the label-free path's set).
# --replicates biological: the user confirmed that each Replicate is a separate culture and mix.
#   Runs that share a Replicate are re-injections of one mix and are averaged, never counted twice.
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
run_level_only <- given("--run-level-only")
replicates_arg <- getarg("--replicates")
num_or <- function(flag, default) {
  v <- suppressWarnings(as.numeric(getarg(flag, as.character(default))))
  if (is.na(v)) stop(sprintf("%s must be a number", flag))
  v
}
min_reps_given <- if (given("--min-replicates")) as.integer(num_or("--min-replicates", 3)) else NULL
q_cut     <- num_or("--q-cutoff", 0.01)
pg_q_cut  <- num_or("--pg-q-cutoff", 0.05)
ch_q_cut  <- num_or("--channel-q-cutoff", 0.01)
adjp_thr  <- num_or("--adjp", 0.05)
logfc_ref <- num_or("--logfc", 1)
if (is.null(input) || is.null(design_path))
  stop("Required: --input <report.parquet> and --design <silac_design.csv>")
if (!file.exists(input)) stop("--input ", input, " does not exist")
if (!file.exists(design_path)) stop("--design ", design_path, " does not exist")
if (!is.null(min_reps_given) && min_reps_given < 2)
  stop("--min-replicates must be at least 2: a one-sample test needs a variance")
if (!is.null(replicates_arg) && !identical(replicates_arg, "biological"))
  stop("--replicates takes 'biological' (each Replicate is a separate culture and mix); runs that ",
       "are re-injections of one mix share a Replicate in the design")

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
source(.sibling("diann_q_columns.R"))   # the q-value columns and cutoffs (one definition)
source(.sibling("skill_version.R"))
source(.sibling("runtime_record.R"))     # the R this ran in, as run_de.R records it
version_label <- paste("run_silac.R, ucdavis-proteomics-core-pipeline",
                       skill_label(skill_version(.script_dir)))

# DE-LIMP rule 2: a value nobody chose says so -- in the record (sources start "DEFAULT") and in
# the client's text (TAG after the value)
DEFAULT_SRC <- "DEFAULT, not user-confirmed"
TAG <- " (DEFAULT \u2014 not user-confirmed)"
setting <- function(flag, value, why)
  list(value = value, given = given(flag),
       source = if (given(flag)) sprintf("given (%s)", flag) else sprintf("%s: %s", DEFAULT_SRC, why))
tag_if <- function(flag) if (given(flag)) "" else TAG
settings <- list(
  min_replicates = list(value = min_reps_given %||% "max(3, ceiling(n / 2)) per comparison",
                        given = !is.null(min_reps_given),
                        source = if (!is.null(min_reps_given)) "given (--min-replicates)" else paste0(
                          DEFAULT_SRC, ": max(3, ceiling(n / 2)) of the n replicates holding a ",
                          "comparison -- 3 of 4 as validated on 4-replicate triple SILAC; with more ",
                          "replicates a protein seen only where the small channel happened to be high ",
                          "would otherwise be tested on a biased few")),
  q_cutoff = setting("--q-cutoff", q_cut, "Q.Value and the experiment-wide columns at DIA-NN's recommended 1%"),
  pg_q_cutoff = setting("--pg-q-cutoff", pg_q_cut,
                        "channel-specific PG.Q.Value; DIA-NN's README: 0.05 is typically sufficient"),
  channel_q_cutoff = setting("--channel-q-cutoff", ch_q_cut,
                             "Channel.Q.Value; DIA-NN's README recommends 0.01-0.5 for plexDIA"),
  adjp = setting("--adjp", adjp_thr, "Benjamini-Hochberg adjusted p below 0.05"),
  logfc = setting("--logfc", logfc_ref, "reference line on the volcano only, never a filter"),
  replicates = list(value = if (!is.null(replicates_arg)) "biological" else "assumed biological",
                    given = !is.null(replicates_arg),
                    source = if (!is.null(replicates_arg))
                      "given (--replicates biological): each Replicate is a separate culture and mix"
                    else paste0(DEFAULT_SRC, ": each Replicate is taken to be a separate culture and ",
                                "mix; ask whether any are re-injections of one mix (references/silac.md)")),
  run_level_only = list(value = run_level_only, given = run_level_only,
                        source = if (run_level_only)
                          "given (--run-level-only): the experiment-wide Global.* / Lib.* q-values were not applied"
                        else paste0(DEFAULT_SRC, ": every q-value column the report has (diann_q_columns.R)")))

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
# The q-value filter: the label-free path's column set (diann_q_columns.R) restricted to what the
# report has -- the run-specific Q.Value, the channel-specific PG.Q.Value and, when present, the
# experiment-wide Global.Q.Value / Global.PG.Q.Value and Lib.* -- or, with --run-level-only, the
# run-specific two alone.
q_cols <- if (run_level_only) c("Q.Value", "PG.Q.Value") else diann_fdr_columns(have)
q_cols <- c("Q.Value", "PG.Q.Value", setdiff(q_cols, c("Q.Value", "PG.Q.Value")))
cutoffs <- vapply(q_cols, function(cc) if (cc == "PG.Q.Value") pg_q_cut
                  else diann_cutoff_for(cc, q_cut), numeric(1))
d_raw <- as.data.table(arrow::read_parquet(input, col_select = intersect(c(need, q_cols, optional), have)))
for (o in optional) if (!o %in% names(d_raw)) d_raw[, (o) := NA_character_]
chan_order <- vapply(labels$channels %||% list(), function(c) c$name, "")
if (!length(chan_order)) {
  ch <- setdiff(unique(d_raw$Channel), "")
  chan_order <- ch[order(match(ch, c("L", "M", "H")), ch)]
}
light <- labels$light_channel %||% chan_order[1]

# The decoy channel DIA-NN modelled channel confidence with, from the search's log when it is
# there (DIA-NN README: channel q-values are optimistic when the light channel is the carrier)
decoy_channel <- local({
  logs <- unique(c(sub("\\.parquet$", ".log.txt", input), file.path(dirname(input), "report.log.txt")))
  hit <- NULL
  for (lg in logs[file.exists(logs)]) {
    x <- grep("Registered \\S+ channel decoy with following masses:", readLines(lg, warn = FALSE), value = TRUE)
    if (length(x)) { hit <- list(masses = trimws(sub("^.*following masses:", "", x[1])), log = basename(lg)); break }
  }
  hit
})

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
runs_per_rep <- ratio_asg[, .(n = uniqueN(Run)), by = Replicate]
reinjected <- runs_per_rep[n > 1]$Replicate
# how a run is named where samples are listed: its replicate, or replicate + run for a re-injection
run_label <- unique(ratio_asg[, .(Run, Replicate)])
run_label[, label := ifelse(Replicate %in% reinjected, paste0(Replicate, "_", Run), Replicate)]

# ---- descriptions: the FASTA headers (UniProt and NCBI) --------------------------------------
fasta_desc <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  h <- sub("^>", "", grep("^>", readLines(path, warn = FALSE), value = TRUE))
  first <- sub("\\s.*$", "", h); rest <- trimws(sub("^\\S+\\s*", "", h))
  uniprot <- grepl("^(sp|tr)\\|[^|]+\\|", first) | grepl("\\bGN=\\S+", rest)
  acc <- ifelse(grepl("^[a-z]{2}\\|[^|]+\\|", first), sub("^[a-z]{2}\\|([^|]+)\\|.*$", "\\1", first), first)
  gene <- ifelse(grepl("\\bGN=\\S+", rest), sub("^.*\\bGN=(\\S+).*$", "\\1", rest),
                 ifelse(grepl("\\[gene=[^]]+\\]", rest), sub("^.*\\[gene=([^]]+)\\].*$", "\\1", rest),
                        NA_character_))
  desc <- sub("\\s+OS=.*$", "", rest)                     # UniProt: up to OS=
  desc <- ifelse(grepl("\\[protein=[^]]+\\]", desc),      # NCBI CDS headers: [protein=...]
                 sub("^.*\\[protein=([^]]+)\\].*$", "\\1", desc), desc)
  desc <- sub("\\s*\\[[^]]*\\]\\s*$", "", desc)            # NCBI: the trailing [organism]
  desc <- sub("^MULTISPECIES:\\s*", "", desc)
  data.table(acc = c(acc, first), desc = c(desc, desc), gene = c(gene, gene),
             uniprot = c(uniprot, uniprot))[!duplicated(acc)]
}
fd <- fasta_desc(fasta_path)
first_acc <- function(pg) sub(";.*$", "", pg)
annot <- unique(d_raw[, .(Protein.Group, Protein.Names, Genes)], by = "Protein.Group")
# DIA-NN reads gene and entry names from UniProt headers (GN=, sp|ACC|ENTRY_ORG). On any other
# header (NCBI: ">WP_000123.1 MULTISPECIES: hypothetical protein [Escherichia]") both columns hold
# its guess, the description's first word ("MULTISPECIES:", "hypothetical", "30S"), which a reader
# takes for a gene name. Those names are blanked -- a [gene=...] token is kept as the gene -- and
# Description names the protein. Without the FASTA the guess shows itself as Genes equal to
# Protein.Names, which a UniProt entry name (ALBU_HUMAN) never is for its gene (ALB).
names_record <- local({
  fa <- if (is.null(fd)) rep(NA_integer_, nrow(annot)) else match(first_acc(annot$Protein.Group), fd$acc)
  by_fasta <- !is.na(fa)
  if (any(by_fasta)) by_fasta[by_fasta] <- !fd$uniprot[fa[by_fasta]]
  guessed <- !by_fasta & !is.na(annot$Genes) & nzchar(annot$Genes) &
    !is.na(annot$Protein.Names) & annot$Genes == annot$Protein.Names
  token <- if (any(by_fasta)) fd$gene[fa[by_fasta]] else character(0)
  annot[by_fasta, `:=`(Genes = token, Protein.Names = NA_character_)]
  annot[guessed, `:=`(Genes = NA_character_, Protein.Names = NA_character_)]
  list(n_protein_groups = nrow(annot), n_not_uniprot = sum(by_fasta | guessed),
       n_gene_from_token = sum(!is.na(token)),
       source = if (!is.null(fd)) paste("the searched FASTA's headers (--fasta)", basename(fasta_path))
                else "no --fasta: Genes equal to Protein.Names marks DIA-NN's guess",
       rule = paste("Genes and Protein.Names are DIA-NN's, read from UniProt headers. For a protein whose",
                    "header is not UniProt's they hold DIA-NN's guess (the description's first word), so",
                    "they are left blank -- a [gene=...] token in the header is kept as the gene -- and",
                    "Description names the protein."))
})
annot[, Description := if (is.null(fd)) NA_character_ else fd$desc[match(first_acc(Protein.Group), fd$acc)]]
annot[, Label := {
  g <- if (is.null(fd)) rep(NA_character_, .N) else fd$gene[match(first_acc(Protein.Group), fd$acc)]
  ifelse(!is.na(g) & nzchar(g), g,
         ifelse(!is.na(Description) & nzchar(Description),
                {x <- sub(",.*$", "", Description); ifelse(nchar(x) > 32, paste0(substr(x, 1, 31), "\u2026"), x)},
                ifelse(!is.na(Genes) & nzchar(Genes), sub(";.*$", "", Genes), first_acc(Protein.Group))))
}]

# ---- the analysis -------------------------------------------------------------------------------
is_cont <- is_contaminant(d_raw$Protein.Ids)
prep <- silac_prepare(d_raw, ratio_asg, cutoffs, ch_q_cut, is_cont)
dd <- prep$d
if (!nrow(dd)) stop("no precursor passed the filters in the runs the design holds")
w <- silac_wide(dd)
loaded <- lapply(setNames(groups, groups), function(g) unique(ratio_asg[Group == g]$Replicate))
min_reps_for <- function(n) if (!is.null(min_reps_given)) min_reps_given else max(3L, as.integer(ceiling(n / 2)))
min_reps_all <- min_reps_for(length(reps_all))      # the detection patterns' threshold
results <- list(); offsets <- list(); tables_rec <- list()
for (i in seq_len(nrow(cmp))) {
  a <- cmp$numerator[i]; b <- cmp$denominator[i]; ct <- cmp$contrast[i]
  reps <- intersect(loaded[[a]], loaded[[b]])
  mr <- min_reps_for(length(reps))
  res <- silac_compare(w, a, b, reps, mr, adjp = adjp_thr)
  tab <- merge(annot[, .(Protein.Group, Protein.Names, Genes, Description)], res$table,
               by = "Protein.Group", all.y = TRUE)
  det <- dd[ok == TRUE & Group %in% c(a, b), .(n = uniqueN(Replicate)), by = .(Protein.Group, Group)]
  for (g in c(a, b)) {
    k <- det[Group == g][match(tab$Protein.Group, Protein.Group)]$n
    tab[, (paste0("Detected_", g)) := sprintf("%d/%d", ifelse(is.na(k), 0L, k), length(loaded[[g]]))]
  }
  first <- c("Protein.Group", "Genes", "Protein.Names", "Description", "logFC", "t", "P.Value",
             "adj.P.Val", "n_replicates", "min_precursors", paste0("Detected_", c(a, b)),
             "log2_ratio_uncentred_mean", "Evidence")
  setcolorder(tab, c(first, setdiff(names(tab), first)))
  setorderv(tab, c("P.Value", "Protein.Group"), na.last = TRUE)
  tab[, comparison := sprintf("%s vs %s (log2 %s/%s)", a, b, a, b)]
  results[[ct]] <- tab
  o <- merge(res$offsets, run_label[, .(Run, label)], by = "Run")
  offsets[[ct]] <- o[, .(comparison = ct, Replicate, Run, offset_log2 = offset, mode_log2 = mode,
                         fold = 2^offset)]
  sig <- tab[!is.na(adj.P.Val) & adj.P.Val < adjp_thr]
  tables_rec[[ct]] <- c(list(
    file = sprintf("DE_silac_%s.csv", make.names(ct)), numerator = a, denominator = b,
    replicates = as.list(reps), n_replicates_available = length(reps), min_replicates = mr,
    n_proteins = nrow(tab), n_tested = res$n_tested, n_significant = nrow(sig),
    n_up = sum(sig$logFC > 0), n_down = sum(sig$logFC < 0),
    n_significant_one_precursor = sum(sig$min_precursors == 1, na.rm = TRUE),
    offset_log2 = setNames(as.list(o$offset), o$label),
    mode_log2 = setNames(as.list(o$mode), o$label),
    zero_point = res$mode_check), res$fit)
  message(sprintf("[run_silac] %s: %d proteins with a ratio; %d tested (ratio in >= %d of %d replicates); %d significant at adj.P < %g (%d up, %d down)",
                  ct, nrow(tab), res$n_tested, mr, length(reps), nrow(sig), adjp_thr,
                  sum(sig$logFC > 0), sum(sig$logFC < 0)))
}
detect <- silac_detection(dd, loaded, min_reps_all)
inten <- silac_intensity(dd)
detect <- merge(annot[, .(Protein.Group, Genes, Protein.Names, Description)], detect, by = "Protein.Group")
specific <- merge(detect[!is.na(pattern)], inten, by = "Protein.Group", all.x = TRUE)
setorderv(specific, c("pattern", "Protein.Group"))
chsum <- silac_channel_summary(dd)
chsum <- merge(chsum, run_label[, .(Run, label)], by = "Run")
chsum <- chsum[order(Replicate, Run, match(Channel, chan_order))]
checks_runs <- asg[role == "label_check"]
incorp <- if (nrow(checks_runs)) {
  dchk <- d_raw[Run %in% checks_runs$Run][!is_contaminant(Protein.Ids)]
  silac_incorporation(dchk, setNames(as.list(checks_runs$Channel), checks_runs$Run), light,
                      cutoffs, ch_q_cut)
} else NULL
lg <- labels$label_group %||% "SILAC"
doubled <- mean(grepl(sprintf("\\(%s\\)\\(%s\\)", lg, lg), d_raw$Precursor.Id, fixed = FALSE))
# protein groups the analysis saw (a channel-passing precursor) that are in no DE table with a
# ratio and in no detection pattern: only in SILAC_sample_matrix.csv
in_tables <- unique(c(unlist(lapply(results, function(x) x$Protein.Group)), specific$Protein.Group))
unlisted <- setdiff(unique(dd[ok == TRUE]$Protein.Group), in_tables)

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
# Per sample (replicate x channel; a re-injection is its own column): log2 of the summed
# channel-passing raw precursor quantity, median-centred within the sample -- for the
# sample-quality panels (sample_quality.py) and the audit's depth checks. Not what the ratios test.
em <- dd[ok == TRUE, .(q = sum(Precursor.Quantity)), by = .(Protein.Group, Replicate, Group, Channel, Run)]
em <- merge(em, run_label[, .(Run, label)], by = "Run")
em[, Sample := paste(label, Group, sep = "_")]
em[, v := log2(q) - stats::median(log2(q)), by = Sample]
samples_tab <- unique(em[, .(File.Name = Sample, Group, Replicate, Channel, Run)])
samples_tab <- samples_tab[order(Replicate, Run, match(Channel, chan_order))]
EM <- dcast(em, Protein.Group ~ Sample, value.var = "v")
setcolorder(EM, c("Protein.Group", samples_tab$File.Name))
EM <- merge(annot[, .(Protein.Group, Genes, Protein.Names)], EM, by = "Protein.Group", all.y = TRUE)
fwrite(EM, file.path(outdir, "SILAC_sample_matrix.csv"))
fwrite(samples_tab, file.path(outdir, "SILAC_samples.csv"))
if (length(incorp)) fwrite(rbindlist(incorp), file.path(outdir, "SILAC_label_incorporation.csv"))
if (!length(incorp) && file.exists(file.path(outdir, "SILAC_label_incorporation.csv")))
  unlink(file.path(outdir, "SILAC_label_incorporation.csv"))

# ---- the comparisons, in words (the Methods, the report and the brief say the same) -----------
cmp_dir <- paste(vapply(names(tables_rec), function(ct) sprintf("%s/%s", tables_rec[[ct]]$numerator,
                                                                tables_rec[[ct]]$denominator), ""), collapse = ", ")
cmp_given <- !is.null(comparisons)
cmp_sentence <- sprintf(paste0("The comparisons were %s, each a log2 ratio numerator/denominator (a positive",
                               " value is higher in the first-named group)%s."),
                        cmp_dir, if (cmp_given) "" else paste0("; every pair of groups, the group in the heavier",
                                                                " channel over the lighter", TAG))

# ---- checks (audit_results.py quotes them: {check, status, message, detail}) -----------------
CHANNEL_SHARE_BAND <- c(0.5, 1.5)   # x the equal share (1/plex) before a channel is "unbalanced"
INCORPORATION_MIN <- 95             # % labelled; SILAC practice asks for > 95%
checks <- list()
# The report's Audit & caveats is a plain-language table (audit_results.client_section_md): each
# check carries its own label and, in `plain`, what its result means for the results in words a
# client reads -- the message stays the technical record (AUDIT.md). This pipeline writes them;
# the audit knows no pipeline by name.
SILAC_CHECK_LABELS <- c(
  silac_channel_balance = "SILAC channel mix", silac_replicates = "SILAC replicates",
  silac_zero_point = "Centring of the SILAC ratios", silac_single_precursor = "Hits resting on one peptide",
  silac_label_residue_bias = "Lysine and arginine labels agree",
  silac_label_incorporation = "Label incorporation",
  silac_label_declaration = "How the label was declared in the search",
  silac_channel_norm = "The search's channel normalisation",
  silac_replicate_independence = "Independent replicates",
  silac_detection_patterns = "Proteins seen in only some channels",
  silac_unlisted_proteins = "Proteins with no ratio")
chk <- function(check, status, message, detail = NULL, plain = NULL)
  checks[[length(checks) + 1]] <<- c(
    list(check = check, status = status, message = message, detail = detail %||% list()),
    if (!is.null(plain)) list(plain = plain),
    if (!is.na(SILAC_CHECK_LABELS[check])) list(label = unname(SILAC_CHECK_LABELS[check])))
eq <- 100 / max(1, length(unique(ratio_asg$Channel)))
share_rng <- chsum[, .(lo = min(share_percent), hi = max(share_percent)), by = Group]
off_band <- share_rng[hi > CHANNEL_SHARE_BAND[2] * eq | lo < CHANNEL_SHARE_BAND[1] * eq]
if (nrow(off_band)) {
  chk("silac_channel_balance", "WARN", sprintf(paste(
    "The channels are not balanced: %s of the identified precursor signal per run, where an",
    "equal mix gives %.0f%%. The uncentred ratios (uncentred_* columns, SILAC_ratio_offsets.csv)",
    "carry this offset; the tested, median-centred ratios remove it, so they read 'relative to the",
    "median detected protein'. The strongest decreases in a small channel can fall below detection there and",
    "then appear only in the detection patterns, never as ratios. Confirm how the samples were",
    "mixed (equal volume is not equal protein)."),
    paste(sprintf("%s %.0f-%.0f%%", share_rng$Group, share_rng$lo, share_rng$hi), collapse = ", "), eq),
    list(equal_share_percent = eq, band = CHANNEL_SHARE_BAND,
         per_group = setNames(lapply(seq_len(nrow(share_rng)), function(i)
           list(min = share_rng$lo[i], max = share_rng$hi[i])), share_rng$Group)),
    plain = paste("The labelled samples were not mixed in equal amounts. The tested ratios are",
                  "centred to remove that offset, so they read relative to the typical protein; the",
                  "strongest decreases in the smaller channel may show only as missing detections."))
} else chk("silac_channel_balance", "PASS", sprintf(
  "Every channel holds %.0f-%.0f%% of the signal per run (an equal mix gives %.0f%%).",
  min(chsum$share_percent), max(chsum$share_percent), eq),
  plain = "The labelled samples were mixed in similar amounts.")
for (ct in names(tables_rec)) {
  t <- tables_rec[[ct]]; n <- t$n_replicates_available; mr <- t$min_replicates
  if (n < mr)
    chk("silac_replicates", "WARN", sprintf(paste(
      "%s: only %d replicate(s) hold both channels, below the %d a ratio needs to be tested; no",
      "protein is tested in this comparison."), ct, n, mr),
      plain = "A comparison has too few replicates holding both channels, so no protein could be tested in it.")
  else if (n == mr)
    chk("silac_replicates", "WARN", sprintf(paste(
      "%s: %d replicates, the minimum: a protein is tested only when every replicate has its ratio,",
      "and the variance rests on %d values per protein (moderated by limma)."), ct, n, n),
      plain = paste("A comparison has the minimum number of replicates: a protein is tested only when",
                    "every replicate measured it, and the statistics rest on few values."))
  else chk("silac_replicates", "PASS", sprintf("%s: %d replicates; a protein is tested with a ratio in at least %d.",
                                               ct, n, mr),
                plain = "The comparison has more replicates than the minimum it needs.")
  z <- t$zero_point
  if (is.finite(z$max_abs_difference %||% NA_real_) && z$max_abs_difference > z$tolerance)
    chk("silac_zero_point", "WARN", sprintf(paste(
      "%s: in at least one run the median protein ratio (the zero of the centred ratios) is %.2f log2",
      "from the density peak (mode) of the ratios, above %.2f. A channel whose strongest losses fall",
      "below detection moves the median. Centred on the mode, %s proteins would be significant (%s up,",
      "%s down) against %d (%d up, %d down) as reported; read the up/down balance with that in mind."),
      ct, z$max_abs_difference, z$tolerance, format(z$n_significant_mode_centred %||% NA),
      format(z$n_up_mode_centred %||% NA), format(z$n_down_mode_centred %||% NA),
      t$n_significant, t$n_up, t$n_down), z,
      plain = paste("In at least one run the ratios' zero point may be shifted, so read the balance",
                    "of increases and decreases with care."))
  else chk("silac_zero_point", "PASS", sprintf(
    "%s: the median protein ratio (the zero) is within %.2f log2 of the density peak in every run.",
    ct, z$tolerance), z,
    plain = "The ratios' zero point agrees with where most proteins sit.")
  if (t$n_significant)
    chk("silac_single_precursor", "INFO", sprintf(paste(
      "%s: %d of the %d significant protein groups rest on a single precursor in at least one",
      "replicate (min_precursors = 1); limma's variance moderation treats them like the rest, so",
      "confirm the ones that matter."), ct, t$n_significant_one_precursor, t$n_significant),
      plain = if (t$n_significant_one_precursor > 0)
        paste("Some significant proteins rest on a single peptide in at least one replicate:",
              "confirm those that matter before relying on them.")
      else "No significant protein rests on a single peptide in any replicate.")
}
# A label that reads K- and R-ending precursors of one protein differently (silac_residue_bias):
# per CHANNEL pair, heavier over lighter, whatever groups they hold (a label swap must not cancel
# it), judged on the median and run by run; a deuterated label is named, since its retention
# shift is the known cause
RESIDUE_BIAS_MAX <- 0.15            # |log2| between K- and R-ending precursors of one protein
RESIDUE_BIAS_MIN_PAIRS <- 20        # protein-run pairs with both, before it is judged
deuterated <- Filter(function(ch) any(grepl("2H", c(ch$K$composition %||% "", ch$R$composition %||% ""),
                                           fixed = TRUE)), labels$channels %||% list())
wc <- silac_wide_channels(dd)
used_ch <- chan_order[chan_order %in% unique(ratio_asg$Channel)]
residue_bias <- list()
for (i in seq_along(used_ch)) for (j in seq_along(used_ch)) if (i < j) {
  lo <- used_ch[i]; hi <- used_ch[j]; key <- paste0(hi, "/", lo)
  kr <- silac_residue_bias(wc, hi, lo)
  holds <- function(ch) paste(unique(ratio_asg[Channel == ch]$Group), collapse = " or ")
  kr$groups <- sprintf("%s holds %s, %s holds %s", hi, holds(hi), lo, holds(lo))
  kr$basis <- paste("per run and protein group: median log2 ratio", key, "of the K-ending precursors",
                    "minus that of the R-ending ones; the median over protein-run pairs, and per run")
  residue_bias[[key]] <- kr
  if (kr$n_pairs < RESIDUE_BIAS_MIN_PAIRS) {
    chk("silac_label_residue_bias", "INFO", sprintf(paste(
      "%s: K- against R-ending precursors not compared: %d protein-run pairs hold both (fewer",
      "than %d)."), key, kr$n_pairs, RESIDUE_BIAS_MIN_PAIRS), kr,
      plain = "Not checked: too few proteins had both kinds of labelled peptide.")
    next
  }
  pr <- unlist(kr$per_run)
  deut <- Filter(function(ch) ch$name %in% c(hi, lo), deuterated)
  if (abs(kr$median) > RESIDUE_BIAS_MAX || stats::median(abs(pr)) > RESIDUE_BIAS_MAX)
    chk("silac_label_residue_bias", "WARN", sprintf(paste(
      "Channels %s (%s): lysine-ending precursors read %+.2f log2 against arginine-ending ones of",
      "the same protein and run (median of %d protein-run pairs, %d protein groups; per run %+.2f",
      "to %+.2f), above %.2f. A protein's ratio between these channels then depends on its mix of",
      "K- and R-ending peptides, by up to that much: read effects of that size with care. %s"),
      key, kr$groups, kr$median, kr$n_pairs, kr$n_proteins, min(pr), max(pr), RESIDUE_BIAS_MAX,
      if (length(deut)) sprintf(paste(
        "The %s channel's lysine label is deuterated (%s): deuterium shifts a peptide's retention",
        "time against its light partner, so the channels are not measured over the same elution",
        "-- the known cause. For a new design, a non-deuterated medium label (K6R6, 13C(6)",
        "lysine and arginine) avoids it (references/silac.md)."),
        paste(vapply(deut, function(ch) ch$name, ""), collapse = ", "),
        paste(unique(vapply(deut, function(ch) sprintf("%s, %s", ch$K$label %||% "?",
                                                         ch$K$composition %||% "?"), "")), collapse = "; "))
      else paste("Neither channel's label is recorded as deuterated: check the labels used, the label",
                 "incorporation and arginine-to-proline conversion (heavy arginine) before reading",
                 "these ratios.")), kr,
      plain = paste("Lysine- and arginine-ending peptides of the same protein give different ratios",
                    "between these channels, so effects of a similar size should be read with care."))
  else chk("silac_label_residue_bias", "PASS", sprintf(paste(
    "Channels %s: lysine- and arginine-ending precursors of the same protein agree (median difference",
    "%+.2f log2 over %d protein-run pairs, per run %+.2f to %+.2f; within %.2f)."), key, kr$median,
    kr$n_pairs, min(pr), max(pr), RESIDUE_BIAS_MAX), kr,
    plain = "Lysine- and arginine-ending peptides of the same protein agree.")
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
                                   "signal and pulls ratios towards the light channel."), INCORPORATION_MIN) else ""),
    plain = if (length(low)) paste("The label was not fully taken up in at least one labelled",
                                   "culture, which pulls ratios towards the unlabelled channel.")
            else sprintf("The label was taken up fully (above %g%%) in the labelled cultures.",
                         INCORPORATION_MIN))
} else chk("silac_label_incorporation", "INFO", paste(
  "Label incorporation was not measured: every run is a mix that holds the light channel by design,",
  "and plexDIA gives each peptide one channel, so a partly labelled peptide is not searched. Measure",
  "it on a labelled-only injection of each labelled culture (a design row with that channel alone), or",
  "from a search with the labels as variable modifications (two-site method, references/silac.md)."),
  plain = "Not measured: label uptake needs an injection of a labelled culture on its own.")
if (doubled > 0)
  chk("silac_label_declaration", "WARN", sprintf(paste(
    "%.1f%% of the precursor ids carry the label twice (%s(%s)(%s)): the search declared it twice",
    "(--lib-fixed-mod on a library predicted from the FASTA), which DIA-NN reports as not supported.",
    "The channels were still assigned, but for a deliverable re-search with estimate_params.py",
    "--labels (references/silac.md)."), 100 * doubled, "K", lg, lg),
    plain = paste("The search declared the label in a way the search engine does not support; the",
                  "channels were still assigned, but the search should be repeated with the label",
                  "declared once."))
if (!is.null(cn) && identical(cn$given, FALSE))
  chk("silac_channel_norm", "INFO", paste(
    "The search's channel normalisation (", cn$flag %||% "?", ") is a DEFAULT, not user-confirmed.",
    "These ratios do not use DIA-NN's normalised quantities, so it does not change them; it changes",
    "the report's normalised columns (Precursor.Normalised, PG.MaxLFQ)."),
    plain = "A search setting was left at its default; it does not change the ratios reported here.")
if (is.null(replicates_arg))
  chk("silac_replicate_independence", "INFO", paste(
    "Whether each replicate is a separate culture and mix was not confirmed: each Replicate in the",
    paste0("design is taken to be one", TAG, "."), "Re-injections of one mix must share a Replicate (they",
    "are then averaged); counted as separate replicates they would make the p-values too small."),
    plain = paste("It was not confirmed that each replicate is a separate culture and mix; if some are",
                  "re-injections of one mix, the statistics would be too optimistic."))
chk("silac_detection_patterns", "INFO", sprintf(paste(
  "%d protein groups are detected (channel q-value passed in >= %d replicates) in some channels and",
  "not detected in any replicate of another: no ratio, so not tested; listed in",
  "SILAC_channel_specific_detections.csv (%s). A detection pattern, not a statistical test: no error",
  "rate, and not detected is not absent, especially in the channel with less protein. Check the list",
  "for whole pathways before calling it a loading or detection-limit effect."), nrow(specific), min_reps_all,
  if (nrow(specific)) paste(sprintf("%s: %d", names(table(specific$pattern)), as.integer(table(specific$pattern))),
                            collapse = "; ") else "none"),
  plain = if (nrow(specific) > 0)
    paste("Some proteins are detected in some channels but not in others. They have no ratio,",
          "so they are listed separately as detection patterns, not tested.")
  else "No protein is detected in some channels while missing from every replicate of another.")
if (length(unlisted))
  chk("silac_unlisted_proteins", "INFO", sprintf(paste(
    "%d protein groups were detected but appear in neither a DE_silac table (no ratio) nor the detection",
    "patterns (detected in too few replicates of every channel): they are only in SILAC_sample_matrix.csv."),
    length(unlisted)),
    plain = "Some proteins were detected too rarely to get a ratio or a detection pattern; they are only in the sample matrix.")

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
q_short <- paste(sprintf("%s <= %g", names(cutoffs), cutoffs), collapse = ", ")
# A subtitle as lines that fit an 8-inch figure: the seven q-value filters ran off its edge
wrap_subtitle <- function(x, width = 90) paste(strwrap(x, width), collapse = "\n")
# the scatter: two comparisons over the same reference by default (M/L against H/L), and the
# correlation their shared channel builds in, with its sign
scatter_pair <- function() {
  if (!is.null(scatter_arg)) return(strsplit(scatter_arg, ":", fixed = TRUE)[[1]])
  nm <- names(tables_rec)
  if (length(nm) < 2) return(NULL)
  for (i in seq_along(nm)) for (j in seq_along(nm))
    if (i < j && tables_rec[[nm[i]]]$denominator == tables_rec[[nm[j]]]$denominator) return(nm[c(i, j)])
  nm[1:2]
}
coupling <- function(x, y) {
  tx <- tables_rec[[x]]; ty <- tables_rec[[y]]
  side <- function(t, g) if (t$numerator == g) "numerator" else if (t$denominator == g) "denominator" else NA
  shared <- intersect(c(tx$numerator, tx$denominator), c(ty$numerator, ty$denominator))
  if (!length(shared)) return("")
  same <- identical(side(tx, shared[1]), side(ty, shared[1]))
  sprintf(paste("; the two share %s as the %s, so its measurement noise alone gives a %s correlation",
                "(about %s0.5 at equal channel noise) even with no biology: their agreement is not",
                "independent evidence"), shared[1], if (same) paste(side(tx, shared[1]), "in both") else
                  "numerator of one and the denominator of the other", if (same) "positive" else "negative",
          if (same) "+" else "\u2212")
}
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
        labs(x = sprintf("log2 %s/%s, relative to the median detected protein", a$numerator, a$denominator),
             y = "-log10 p-value", title = sprintf("%s vs %s", a$numerator, a$denominator),
             subtitle = wrap_subtitle(sprintf("%s tested (ratio in >= %d of %d replicates); %s significant at adj. p < %g",
                                format(a$n_tested, big.mark = ","), a$min_replicates, a$n_replicates_available,
                                format(a$n_significant, big.mark = ","), adjp_thr))) + th
      if (nrow(lab)) g <- g + (if (has_repel) ggrepel::geom_text_repel(data = lab, aes(label = Label), size = 2.7,
                                                                         colour = ink, max.overlaps = 20)
                               else geom_text(data = lab, aes(label = Label), size = 2.7, colour = ink,
                                              vjust = -0.8, check_overlap = TRUE))
      ggsave(file.path(figdir, f), g, width = 8, height = 6.2, dpi = 200, bg = "white")
      fig_ok(f, "volcano", sprintf(paste(
        "Volcano plot, %s vs %s: each point a protein group tested in at least %d of %d replicates;",
        "x the mean log2 ratio %s/%s relative to the median detected protein (so the up/down split is",
        "close to even by construction), y the moderated one-sample t-test p-value.",
        "%d significant at Benjamini-Hochberg adj. p < %g (%d up, %d down); dotted lines mark |log2| = %g",
        "for reference only. All comparisons share the x-axis."), a$numerator, a$denominator, a$min_replicates,
        a$n_replicates_available, a$numerator, a$denominator, a$n_significant, adjp_thr, a$n_up, a$n_down,
        logfc_ref))
    })
  }
  f <- "silac_channel_composition.png"
  draw(f, "silac_channel_composition", {
    x <- copy(chsum); x[, Group := factor(Group, levels = rev(groups))]
    g <- ggplot(x, aes(share_percent, label, fill = Group)) + geom_col(width = 0.62, colour = "white") +
      geom_text(aes(label = sprintf("%.0f%%", share_percent)), position = position_stack(vjust = 0.5),
                colour = "white", size = 3.4, fontface = "bold") +
      geom_vline(xintercept = seq(eq, 100 - eq, by = eq), linetype = "dashed", colour = muted, linewidth = 0.4) +
      scale_fill_manual(values = setNames(c("#6b7a74", dn_c, up_c, "#8c6bb1")[seq_along(groups)], groups),
                        breaks = groups, name = NULL) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.03)), labels = function(v) paste0(v, "%")) +
      labs(x = "Share of the identified precursor signal", y = NULL, title = "Signal from each channel",
           subtitle = "Raw precursor quantities passing the channel q-value; dashed lines mark equal shares") + th
    ggsave(file.path(figdir, f), g, width = 8, height = 1.6 + 0.45 * length(unique(x$label)), dpi = 200, bg = "white")
    fig_ok(f, "silac_channel_composition", sprintf(paste(
      "Each channel's share of the identified precursor signal (raw Precursor.Quantity of precursors",
      "passing channel q <= %g, contaminants removed), per run; dashed lines mark an equal mix."), ch_q_cut))
  })
  f <- "silac_channel_detections.png"
  draw(f, "silac_channel_detections", {
    x <- copy(chsum); x[, Group := factor(Group, levels = groups)]
    g <- ggplot(x, aes(label, precursors, fill = Group)) + geom_col(position = position_dodge(0.8), width = 0.75) +
      scale_fill_manual(values = setNames(c("#6b7a74", dn_c, up_c, "#8c6bb1")[seq_along(groups)], groups), name = NULL) +
      scale_y_continuous(labels = function(v) format(v, big.mark = ",")) +
      labs(x = NULL, y = "Precursors", title = "Precursors identified in each channel",
           subtitle = wrap_subtitle(sprintf("%s, channel q <= %g", q_short, ch_q_cut))) + th
    ggsave(file.path(figdir, f), g, width = 8, height = 4.6, dpi = 200, bg = "white")
    fig_ok(f, "silac_channel_detections", sprintf(paste(
      "Precursors identified per run and channel (%s, channel q <= %g, contaminants removed)."),
      q_short, ch_q_cut))
  })
  sc <- scatter_pair()
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
      fig_ok(f, "silac_scatter", sprintf("Centred log2 ratios of the proteins tested in both comparisons: %s (x) against %s (y)%s.",
                                         gsub("-", "/", sc[1]), gsub("-", "/", sc[2]), coupling(sc[1], sc[2])))
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
q_text <- function(cc) {
  run_level <- c(Q.Value = "the run-specific Q.Value", PG.Q.Value = "the run- and channel-specific protein-group PG.Q.Value",
                 Global.Q.Value = "the experiment-wide Global.Q.Value", Global.PG.Q.Value = "the experiment-wide Global.PG.Q.Value",
                 Lib.Q.Value = "the library Lib.Q.Value", Lib.PG.Q.Value = "the library Lib.PG.Q.Value")
  flag <- if (cc == "PG.Q.Value") "--pg-q-cutoff" else "--q-cutoff"
  sprintf("%s <= %g%s", if (cc %in% names(run_level)) run_level[[cc]] else cc, cutoffs[[cc]], tag_if(flag))
}
filters <- c(sprintf("%s <= %g", names(cutoffs), cutoffs),
             sprintf("Channel.Q.Value <= %g in both channels of a ratio", ch_q_cut),
             "Precursor.Quantity > 0", "contaminant precursors removed (any accession tagged Cont_ / contam_)")
per_cmp_n <- paste(vapply(names(tables_rec), function(ct) sprintf("%s/%s: %d of %d", tables_rec[[ct]]$numerator,
                                                                 tables_rec[[ct]]$denominator, tables_rec[[ct]]$min_replicates,
                                                                 tables_rec[[ct]]$n_replicates_available), ""), collapse = "; ")
methods_paragraph <- paste0(
  sprintf("SILAC ratios were analysed with %s (limma %s; R %s). ", version_label, pkg("limma"),
          paste(R.version$major, R.version$minor, sep = ".")),
  sprintf("From DIA-NN's main report, precursors were kept at %s; precursors mapping to a contaminant entry were removed. ",
          paste(vapply(names(cutoffs), q_text, ""), collapse = ", ")),
  sprintf(paste0("A precursor ratio required both channels at Channel.Q.Value <= %g%s with a non-zero raw quantity ",
                 "(Precursor.Quantity; DIA-NN's channel-normalised quantities were not used). "),
          ch_q_cut, tag_if("--channel-q-cutoff")),
  cmp_sentence, " ",
  paste0("For each comparison, the log2 ratio of a protein group in each run was the median of its precursor ",
         "log2 ratios, median-centred within the run (on the median ratio of all its protein groups), so that ",
         "zero is the median detected protein and an even split of increases and decreases is built in; the ",
         "uncentred ratios, which carry the mixing offset, are reported alongside. "),
  if (length(reinjected)) sprintf("Runs that were re-injections of one mix were averaged into one replicate (%s). ",
                                  paste(sprintf("%s: %d runs", runs_per_rep[n > 1]$Replicate, runs_per_rep[n > 1]$n),
                                        collapse = "; ")) else "",
  if (!is.null(replicates_arg)) "Each replicate was a separate culture and mix. " else
    paste0("Each replicate was taken to be a separate culture and mix", TAG, ". "),
  sprintf(paste0("Protein groups with a ratio in at least the minimum number of replicates (%s%s) were tested ",
                 "against a log2 ratio of 0 with a moderated one-sample t-test (limma lmFit with an intercept-only ",
                 "design and eBayes with robust = TRUE), and p-values were adjusted per comparison by the ",
                 "Benjamini-Hochberg method; protein groups with adj.P.Val < %g%s were called significant, with no ",
                 "fold-change filter. "),
          per_cmp_n, if (is.null(min_reps_given)) paste0("; max(3, half the replicates)", TAG) else "",
          adjp_thr, tag_if("--adjp")),
  sprintf(paste0("Protein groups detected (Channel.Q.Value <= %g) in at least %d replicates of some channels and ",
                 "not detected in any replicate of another have no ratio; they are listed separately as a detection ",
                 "pattern, which is not a statistical test (no error rate; not detected is not absent). "),
          ch_q_cut, min_reps_all),
  if (!is.null(decoy_channel)) sprintf("Channel-specific q-values were estimated by DIA-NN against a decoy channel (%s, from the search log). ",
                                       decoy_channel$masses) else "",
  if (length(incorp)) "Label incorporation was estimated from label-check runs as the labelled share of each precursor's signal (a lower bound). " else
    "Label incorporation was not measured from these data. ",
  "Citations: plexDIA, Derks et al., Nat Biotechnol 2023;41:50-59 (doi:10.1038/s41587-022-01389-w); limma, Ritchie et al., ",
  "Nucleic Acids Res 2015;43:e47; robust empirical Bayes, Phipson et al., Ann Appl Stat 2016;10:946-963 ",
  "(doi:10.1214/16-AOAS920); Benjamini and Hochberg, J R Stat Soc B 1995;57:289-300.")
off_txt <- paste(vapply(names(tables_rec), function(ct) {
  o <- unlist(tables_rec[[ct]]$offset_log2)
  sprintf("%s %+.2f to %+.2f", ct, min(o), max(o))
}, ""), collapse = "; ")
# Precision, per comparison: a typical estimate's 95% interval, and the cutoff the correction
# applied -- said for what it is. With no call it is BH's first step, its strictest, which is a
# property of the cutoff and never set beside a comparison that had calls (silac_ratios.R).
pow_txt <- paste(vapply(names(tables_rec), function(ct) {
  t <- tables_rec[[ct]]
  if (!is.finite(t$critical_log2_raw_p05 %||% NA_real_)) return(sprintf("%s: not estimated", ct))
  sprintf("%s: a typical protein's estimate is within about \u00b1%.2f log2 (95%% interval; a change that size reaches raw p < 0.05); %s",
          ct, t$critical_log2_raw_p05,
          if (identical(t$p_cutoff_kind, "first_step"))
            sprintf(paste("because no protein passed, the correction applied its strictest cutoff (p < %.2g),",
                          "about %.2f log2 at typical variance for a protein changing alone -- a property of the",
                          "cutoff, not of the replicates or the variance: proteins changing together would be",
                          "called at a lower bar, so it is not comparable with a comparison that had calls"),
                    t$p_cutoff, t$critical_log2_significant)
          else sprintf("with %d protein groups called, the correction's cutoff was p < %.2g, about %.2f log2 at typical variance",
                       t$n_significant, t$p_cutoff, t$critical_log2_significant))
}, ""), collapse = "; ")
one_prec_txt <- paste(vapply(names(tables_rec), function(ct) sprintf("%s %d of %d", ct,
  tables_rec[[ct]]$n_significant_one_precursor, tables_rec[[ct]]$n_significant), ""), collapse = "; ")
brief_notes <- list(
  paste("logFC in the DE_silac tables is the MEDIAN-CENTRED log2 ratio (numerator/denominator): zero is the",
        "median detected protein of each run, so an even split of increases and decreases is built in. It is",
        "not 'amount per volume'."),
  cmp_sentence,
  sprintf(paste("The uncentred ratios (uncentred_<replicate>, log2_ratio_uncentred_mean) are amount per",
                "volume and carry each run's mixing offset; the median uncentred log2 ratio per run:",
                "%s (SILAC_ratio_offsets.csv)."), off_txt),
  paste("A ratio needs the precursor in both channels. A protein lost in the smaller channel has no ratio and",
        "appears only in SILAC_channel_specific_detections.csv, worded 'detected in X (k/n); not detected in Y (0/n)':",
        "a detection pattern, not a statistical test -- no error rate, and not detected is not absent, especially",
        "in the channel with less protein. Its intensity percentile in the channel where it is detected tells an",
        "abundant protein from one near the limit. Read the list for whole pathways before calling it a loading effect."),
  paste("Two comparisons that share a channel are correlated by construction (positive when the shared channel",
        "is the denominator of both, negative when it is the numerator of one and the denominator of the other):",
        "their agreement is not independent evidence."),
  sprintf(paste("Significant protein groups resting on a single precursor in at least one replicate: %s",
                "(min_precursors = 1; confirm the ones that matter)."), one_prec_txt),
  sprintf(paste("Not significant is not unchanged. From the median variance, a rough reading of precision:",
                "%s. These numbers rule nothing out: a protein not called may still have changed."), pow_txt),
  if (names_record$n_not_uniprot)
    sprintf(paste("%d of %d protein groups come from FASTA headers that are not UniProt's: their Genes and",
                  "Protein.Names are blank (DIA-NN's value there is the description's first word, such as",
                  "'MULTISPECIES:' or 'hypothetical') -- name them by Description."),
            names_record$n_not_uniprot, names_record$n_protein_groups))
brief_notes <- Filter(Negate(is.null), brief_notes)
plain_language <- paste(
  "**SILAC ratios + limma**: each sample was grown with its own isotope label and the samples were mixed",
  "and measured together, so each protein is compared between labels inside the same run. The ratio is",
  "centred so that the median detected protein reads 0, and a protein is called changed when its ratio",
  "differs from 0 consistently across the replicates (limma's moderated t-test, Benjamini-Hochberg adjusted).")
md_row <- function(...) paste0("| ", paste(..., sep = " | "), " |")
fmt_n <- function(x) format(x, big.mark = ",")
ch_md <- c(sprintf(paste("Each run holds %s. Signal share = the channel's part of the identified precursor",
                         "signal (raw quantity, channel q <= %g); an equal mix gives %.0f%%."),
                   paste(sprintf("%s (%s)", unique(ratio_asg[, .(Channel, Group)])$Group,
                                 unique(ratio_asg[, .(Channel, Group)])$Channel), collapse = ", "), ch_q_cut, eq),
           "", cmp_sentence,
           "", md_row("Replicate", "Channel", "Group", "Precursors", "Protein groups", "Signal share"),
           "|---|---|---|---:|---:|---:|",
           vapply(seq_len(nrow(chsum)), function(i) md_row(chsum$label[i], chsum$Channel[i], chsum$Group[i],
                                                           fmt_n(chsum$precursors[i]), fmt_n(chsum$protein_groups[i]),
                                                           sprintf("%.0f%%", chsum$share_percent[i])), ""),
           "", sprintf(paste("Median uncentred log2 ratio per run (the mixing offset the centred ratios",
                             "remove): %s. Label incorporation: %s"), off_txt,
                       if (length(incorp)) paste(vapply(incorp, function(x) sprintf("%s %.1f%%", x$run, x$median_percent), ""),
                                                 collapse = "; ")
                       else "not measured from these data (no labelled-only run)."),
           "", sprintf("Significant protein groups that rest on a single precursor in at least one replicate: %s.",
                       one_prec_txt))
sp_counts <- if (nrow(specific)) table(specific$pattern) else integer(0)
sp_md <- c(sprintf(paste("%s protein groups are detected (channel q <= %g in at least %d replicates) in some",
                         "channels and not detected in any replicate of another. They have no ratio, so they are",
                         "not in the tables above. This is a detection pattern, not a statistical test: it has no",
                         "error rate, and not detected is not absent, especially in the channel with less protein",
                         "-- the strongest losses in a small channel can only appear here. The percentile is the",
                         "protein's abundance among all protein groups detected in that channel: high means it is",
                         "abundant where present, low that it sits near the detection limit. Full list:",
                         "`SILAC_channel_specific_detections.csv`."), fmt_n(nrow(specific)), ch_q_cut, min_reps_all),
           "", md_row("Pattern", "Protein groups"), "|---|---:|",
           vapply(names(sp_counts), function(k) md_row(k, fmt_n(as.integer(sp_counts[[k]]))), ""))
if (nrow(specific)) {
  top_sp <- merge(specific, annot[, .(Protein.Group, Label)], by = "Protein.Group")
  top_sp[, total := rowSums(.SD), .SDcols = paste0(groups, "_replicates")]
  pct_where <- function(i) {
    present <- groups[vapply(groups, function(g) top_sp[[paste0(g, "_replicates")]][i] > 0, logical(1))]
    paste(vapply(present, function(g) {
      v <- top_sp[[paste0("percentile_", g)]][i]
      if (is.null(v) || is.na(v)) sprintf("%s: n/a", g) else sprintf("%s: %d", g, as.integer(v))
    }, ""), collapse = ", ")
  }
  top_sp[, where := vapply(seq_len(.N), pct_where, "")]
  setorderv(top_sp, c("pattern", "total", "Protein.Group"), c(1, -1, 1))
  top_sp <- top_sp[, head(.SD, 8), by = pattern]
  sp_md <- c(sp_md, "", "Up to 8 per pattern:", "",
             md_row("Protein", "Description", "Pattern", "Percentile where detected"), "|---|---|---|---|",
             vapply(seq_len(nrow(top_sp)), function(i) md_row(
               top_sp$Label[i], gsub("\\|", "/", substr(ifelse(is.na(top_sp$Description[i]), "", top_sp$Description[i]), 1, 60)),
               top_sp$pattern[i], top_sp$where[i]), ""))
}
report_sections <- list(list(title = "SILAC channels and mixing", markdown = paste(ch_md, collapse = "\n")),
                        list(title = "Channel-specific proteins", markdown = paste(sp_md, collapse = "\n")))
high_significance_note <- paste(
  "In a SILAC ratio test each protein is compared inside the same run, so run-to-run variation cancels and",
  "the test is powerful; between strains or treatments that differ broadly, a large share can be real. Check",
  "the channel assignment (a swapped label reverses a whole comparison) and the mixing (SILAC_channel_summary.csv);",
  "after median-centring, an even split of increases and decreases is built in.")
sig_counts <- lapply(tables_rec, function(t) list(up = t$n_up, down = t$n_down, total = t$n_significant))
# The contaminant step in run_de.R's record shape (contaminants.R contaminant_record): the brief,
# AGENTS.md and the Methods' readers say what was removed from it, not "not recorded".
n_cont_prec <- uniqueN(d_raw[is_cont & Run %in% ratio_asg$Run]$Precursor.Id)
contaminants_rec <- list(
  policy = if (n_cont_prec == 0) "none_present" else "removed", removed = n_cont_prec > 0,
  tag = CONTAMINANT_TAG, tags = as.list(CONTAMINANT_TAGS), id_column = "Protein.Ids",
  pattern = contaminant_regex(), unit = "precursors", unit_singular = "precursor",
  unit_level = "precursor",
  rule = sprintf("a precursor is a contaminant when any accession in Protein.Ids starts with %s -- applied by run_silac.R",
                 paste0("'", CONTAMINANT_TAGS, "'", collapse = " or ")),
  rule_origin = "applied by run_silac.R",
  counted_after = "the runs in the design, before the q-value filters",
  n_precursors = n_cont_prec, n_protein_groups = prep$counts$contaminant_protein_groups)
input_sha256 <- tryCatch(if (requireNamespace("digest", quietly = TRUE))
  digest::digest(input, algo = "sha256", file = TRUE) else NA_character_, error = function(e) NA_character_)
replay <- c("--input", normalizePath(input), "--design", normalizePath(design_path),
            "--outdir", normalizePath(outdir), "--figures-dir", normalizePath(figdir),
            if (!is.null(params_path)) c("--params", normalizePath(params_path)),
            if (!is.null(fasta_path)) c("--fasta", normalizePath(fasta_path)),
            "--comparisons", paste(names(tables_rec), collapse = ","),
            if (!is.null(min_reps_given)) c("--min-replicates", min_reps_given),
            if (!is.null(replicates_arg)) c("--replicates", replicates_arg),
            if (run_level_only) "--run-level-only",
            "--q-cutoff", q_cut, "--pg-q-cutoff", pg_q_cut,
            "--channel-q-cutoff", ch_q_cut, "--adjp", adjp_thr, "--logfc", logfc_ref)
prov <- list(
  pipeline_id = "silac_ratio", method = "silac",
  display_label = "SILAC ratios + limma (one-sample moderated t, robust eBayes)",
  methods_heading = "SILAC ratio analysis", methods_paragraph = methods_paragraph,
  plain_language = plain_language, brief_notes = brief_notes,
  rollup_method = paste("median of the precursor log2 ratios (raw Precursor.Quantity) per protein group and run,",
                        "averaged over a replicate's runs"),
  de_engine = "limma::lmFit (intercept only, per comparison) -> eBayes(robust = TRUE) -> BH",
  missing_policy = sprintf(paste("A ratio needs both channels at Channel.Q.Value <= %g; a protein group is tested",
                                 "with a ratio in at least the comparison's minimum of replicates (%s); one detected",
                                 "in some channels and not in another is listed as a detection pattern in",
                                 "SILAC_channel_specific_detections.csv, not tested."), ch_q_cut, per_cmp_n),
  citation = paste("Derks J et al. (2023) Nat Biotechnol 41:50-59 (plexDIA); Ritchie ME et al. (2015) NAR 43:e47 (limma);",
                   "Phipson B et al. (2016) Ann Appl Stat 10:946-963 (robust eBayes)"),
  quantities = "raw",
  normalisation = paste("per-run median centring of the protein log2 ratios: zero is the median detected protein",
                        "(uncentred ratios kept); DIA-NN's channel-normalised quantities are not used (raw Precursor.Quantity)"),
  normalization_check = list(status = "not_applicable", reason = paste(
    "A SILAC ratio compares channels inside one run, so DIA-NN's between-run normalisation does not",
    "enter it; each run's ratios are median-centred instead, and the uncentred ratios are kept",
    "beside them. The channel-normalisation answer is the search's (references/silac.md)."),
    quantities_applied = "raw Precursor.Quantity, median-centred ratios"),
  set_tests = list(available = FALSE, suggested_for = list(), reason = paste(
    "run_sets.R reads run_de.R's group-vs-group fit; protein-set tests on SILAC ratios are not",
    "part of skill 2.13")),
  figures = list(written_by = "run_silac.R", figures_json = normalizePath(file.path(figdir, "figures.json"))),
  replay = list(script = "run_silac.R", args = as.list(as.character(replay))),
  contrasts = as.list(names(tables_rec)), adjp = adjp_thr, logfc = logfc_ref,
  logfc_role = "reference_line_only", significance_rule = "adj.P.Val < adjp (BH, per comparison); no fold-change filter",
  design = "an intercept-only (one-sample) model per comparison, one value per replicate",
  n_samples = n_rep, n_replicates = n_rep,
  groups = setNames(lapply(loaded, length), names(loaded)),
  sample_matrix = list(file = "SILAC_sample_matrix.csv", samples = "SILAC_samples.csv", values = paste(
    "per sample (replicate x channel): log2 of the summed channel-passing raw precursor quantity,",
    "median-centred within the sample; for the sample-quality panels and depth checks, not what was tested"),
    proteins_only_here = length(unlisted)),
  report_sections = report_sections, high_significance_note = high_significance_note,
  contaminants = contaminants_rec,
  # how the brief (analysis_prompt.py) tells the writer to read a hit -- in place of run_de.R's
  # measured-versus-inferred tiers, CVs and next-step tiers, which do not apply to a ratio
  evidence_reading = paste(
    "Every tested protein group's logFC is a ratio measured inside each run, both channels",
    "passing, in at least the comparison's minimum of replicates (`n_replicates`;",
    "`Detected_<group>` gives k/n): there are no measured-versus-inferred tiers and no",
    "missing-value model. A row whose `Evidence` reads 'ratio in k of n replicates: not tested'",
    "has no p-value and is not a hit; count tested proteins as the rows with a P.Value. A protein",
    "group detected in some channels and not in another has no ratio: it is in",
    "`SILAC_channel_specific_detections.csv` as a detection pattern ('detected in X (k/n); not",
    "detected in Y (0/n)') -- report it that way, never with a fold change or as a test result."),
  next_steps = list(
    "**Follow up first** -- significant, with at least two precursors behind every replicate's ratio (`min_precursors` >= 2).",
    "**Check** -- significant but resting on a single precursor in at least one replicate (`min_precursors` = 1): confirm before quoting it.",
    "**Confirm presence** -- the channel-specific detections: a detection pattern, not a test; confirm with an orthogonal method before calling a protein absent."),
  spread = paste(
    "a protein's reproducibility is the spread of its replicate ratios (the `centred_<replicate>`",
    "columns), which limma's moderated variance summarises -- not a CV of intensities; do not",
    "compute CVs from SILAC_sample_matrix.csv"),
  # the README's line for output/tables (session_docs.py reads it instead of run_de.R's)
  tables_summary = paste("one `DE_silac_*.csv` per comparison (median-centred log2 ratios, with the",
                         "uncentred ones), `SILAC_channel_specific_detections.csv`,",
                         "`SILAC_channel_summary.csv`, `SILAC_ratio_offsets.csv`, `methods.txt`"),
  de_columns = list(
    logFC = paste("the mean median-centred log2 ratio, numerator/denominator, over the replicates with a",
                  "ratio (limma's estimate): zero is the median detected protein of each run; positive =",
                  "higher in the numerator (the first-named group); NA = not tested"),
    t = "moderated one-sample t-statistic (limma eBayes, robust)",
    P.Value = "raw p-value of the moderated one-sample t-test against a log2 ratio of 0; NA = not tested (Evidence says why)",
    adj.P.Val = "Benjamini-Hochberg adjusted p-value within the comparison -- THE significance column (NA = not tested)",
    Genes = paste("the gene name from a UniProt header (GN=; DIA-NN's) or a [gene=...] token; blank for",
                  "other headers, where DIA-NN's value is the description's first word"),
    Protein.Names = "the UniProt entry name (DIA-NN's); blank for a header that is not UniProt's",
    Description = "the protein's description from the searched FASTA's header",
    n_replicates = "replicates in which the protein group had a ratio",
    min_precursors = "the fewest precursors behind any one of its replicate ratios (1 = a single precursor)",
    `Detected_<group>` = paste("k/n: the replicates holding that group in which the protein group had at least",
                               "one precursor passing the channel q-value in that group's channel"),
    log2_ratio_uncentred_mean = paste("the mean UNCENTRED log2 ratio: amount per volume, carrying each",
                                      "run's mixing offset (SILAC_ratio_offsets.csv)"),
    Evidence = "tested, or why the protein group was not (a ratio in too few replicates)",
    `centred_<replicate>` = "that replicate's log2 ratio relative to the median detected protein (the mean over its runs)",
    `uncentred_<replicate>` = "that replicate's log2 ratio before centring (the median of its precursor log2 ratios)",
    `n_precursors_<replicate>` = "distinct precursors with both channels passing behind that replicate's ratio (over its runs)",
    comparison = "the comparison, in words"),
  glance_tiles = list(list(n_rep, if (length(reinjected)) sprintf("replicates (%d runs)", uniqueN(ratio_asg$Run))
                                  else "replicates (runs)"), list(length(groups), "channels compared"),
                      list(nrow(specific), "channel-specific proteins")),
  q_columns = as.list(c(names(cutoffs), "Channel.Q.Value")), q_cutoffs = as.list(c(unname(cutoffs), ch_q_cut)),
  filters_applied = as.list(filters), filter_counts = prep$counts,
  settings = settings,
  comparisons_source = des$comparisons_source, comparisons_sentence = cmp_sentence,
  significant_per_contrast = sig_counts, de_tables = tables_rec,
  labels = labels, labels_source = labels_source, protein_names = names_record,
  label_residue_bias = residue_bias,
  decoy_channel = decoy_channel %||% list(masses = NULL, note = "no DIA-NN log beside the report named it"),
  # searched runs the user decided to leave out (silac.py design --exclude), in run_de.R's shape
  runs_left_out = jsonlite::fromJSON(jsonlite::toJSON(des$runs_left_out, auto_unbox = TRUE, null = "null"),
                                     simplifyVector = FALSE),
  design_record = list(file = normalizePath(design_path), channels = as.list(chan_order),
                       label_swapped_groups = des$label_swapped_groups %||% list(),
                       reinjected_replicates = des$reinjected_replicates %||% list(),
                       label_check_runs = des$label_check_runs %||% list(),
                       unused_runs = des$unused_runs %||% list()),
  channel_summary = lapply(split(chsum, by = "Run"), function(x)
    setNames(lapply(seq_len(nrow(x)), function(i) list(channel = x$Channel[i], share_percent = x$share_percent[i],
                                                       precursors = x$precursors[i], protein_groups = x$protein_groups[i])),
             x$Group)),
  channel_specific = list(file = "SILAC_channel_specific_detections.csv", n = nrow(specific),
                          min_replicates = min_reps_all, by_pattern = as.list(table(specific$pattern)),
                          kind = "a detection pattern, not a statistical test: no error rate; not detected is not absent"),
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
  runtime = runtime_record("run_silac.R"))
writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA),
           file.path(outdir, "de_provenance.json"))
writeLines(c("SILAC ratio analysis (run_silac.R)", "", strwrap(methods_paragraph, 100), "",
             "Settings:", vapply(names(settings), function(k) sprintf("  %s = %s (%s)", k, format(settings[[k]]$value),
                                                                     settings[[k]]$source), "")),
           file.path(outdir, "methods.txt"))
con <- file(file.path(outdir, "sessionInfo.txt"), "w"); sink(con); print(utils::sessionInfo()); sink(); close(con)

# ---- reproducibility_log.R: the statistics as they ran, with every value written out ---------
rq <- function(x) paste0("'", gsub("'", "\\\\'", gsub("\\\\", "\\\\\\\\", as.character(x))), "'")
rvec <- function(x) paste0("c(", paste(vapply(x, rq, ""), collapse = ", "), ")")
rnum <- function(x) paste0("c(", paste(sprintf("%s = %s", vapply(names(x), rq, ""), format(unname(x), digits = 15)),
                                       collapse = ", "), ")")
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
        sprintf("min_reps <- c(%s)", paste(vapply(names(tables_rec), function(ct)
          sprintf("%s = %dL", rq(ct), as.integer(tables_rec[[ct]]$min_replicates)), ""), collapse = ", ")),
        sprintf("cutoffs <- %s   # the q-value columns, in order", rnum(cutoffs)),
        sprintf("channel_q <- %s; adjp <- %s", ch_q_cut, adjp_thr),
        sprintf("contaminant <- %s   # contaminants.R: any accession starting with a tag", rq(contaminant_regex())),
        "", "# ---- the run ----",
        sprintf("d <- as.data.table(arrow::read_parquet(report, col_select = %s))", rvec(intersect(c(need, q_cols), have))),
        "prep <- silac_prepare(d, design, cutoffs, channel_q, grepl(contaminant, d$Protein.Ids))",
        "w <- silac_wide(prep$d)",
        "loaded <- lapply(split(design$Replicate, design$Group), unique)",
        "for (ct in names(comparisons)) {",
        "  a <- comparisons[[ct]][1]; b <- comparisons[[ct]][2]",
        "  res <- silac_compare(w, a, b, intersect(loaded[[a]], loaded[[b]]), min_reps[[ct]], adjp = adjp)",
        "  fwrite(res$table, sprintf('replay_DE_silac_%s.csv', make.names(ct)))",
        "  cat(ct, ': tested', res$n_tested, ' significant', sum(res$table$adj.P.Val < adjp, na.rm = TRUE), '\\n')",
        "}")
writeLines(rl, file.path(outdir, "reproducibility_log.R"))

message(sprintf("[run_silac] done: %d comparison(s), %d channel-specific protein groups, %d figure(s)%s -> %s",
                length(results), nrow(specific), length(figs),
                if (length(failed)) sprintf(" (%d not drawn)", length(failed)) else "", normalizePath(outdir)))
