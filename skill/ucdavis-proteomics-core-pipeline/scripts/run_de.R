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
# metadata CSV columns: File.Name,Group[,Batch,Covariate1,Covariate2]
#   File.Name must match the Run / column names in the report.
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
outdir    <- getarg("--outdir", "de_results")
logfc_thr <- as.numeric(getarg("--logfc", "1.0"))
adjp_thr  <- as.numeric(getarg("--adjp", "0.05"))

if (is.null(input) || is.null(meta_path))
  stop("Required: --input <report> and --metadata <conditions.csv>")
if (!method %in% c("dpc", "maxlfq"))
  stop("--method must be 'dpc' or 'maxlfq'")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- read experimental design -----------------------------------------------
meta <- utils::read.csv(meta_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("File.Name", "Group") %in% names(meta)))
  stop("metadata CSV must have at least File.Name and Group columns")
covariates <- intersect(c("Batch", "Covariate1", "Covariate2"), names(meta))

message(sprintf("[run_de] method=%s  q=%.3f  samples=%d  covariates=%s",
                method, q_cutoff, nrow(meta),
                if (length(covariates)) paste(covariates, collapse = ",") else "none"))

# ---- build the protein-level object per pipeline ----------------------------
# Both branches produce:  E (proteins x samples, log2), run_names (cols of E),
# genes (annotation df), and a `descriptor` describing the method (provenance).

descriptor <- NULL

if (method == "dpc") {
  if (!requireNamespace("limpa", quietly = TRUE))
    stop("limpa is required for --method dpc. BiocManager::install('limpa') (needs R 4.5+, Bioc 3.22+).")
  suppressMessages({ library(limpa); library(limma) })

  dat <- limpa::readDIANN(input, format = format, q.cutoffs = q_cutoff)
  message(sprintf("[run_de] readDIANN: %d precursors x %d runs", nrow(dat$E), ncol(dat$E)))

  dpcfit    <- limpa::dpcCN(dat)
  y_protein <- limpa::dpcQuant(dat, "Protein.Group", dpc = dpcfit)
  E         <- y_protein$E
  run_names <- colnames(E)
  genes     <- y_protein$genes
  descriptor <- list(
    pipeline_id   = "dpc",
    display_label = "DPC-Quant + limma (limpa)",
    rollup_method = "DPC-Quant (Detection Probability Curve quantification, dpcCN)",
    de_engine     = "limpa::dpcDE (voomaLmFitWithImputation) -> contrasts.fit -> eBayes",
    missing_policy = "Missing precursors modelled via the detection probability curve; not imputed, not dropped.",
    citation      = "Li M, Cobbold SA, Smyth GK (2025) bioRxiv 10.1101/2025.04.28.651125; Li M, Smyth GK (2023) Bioinformatics 39(5):btad200"
  )

} else { # maxlfq
  if (!requireNamespace("arrow", quietly = TRUE) && identical(format, "parquet"))
    stop("arrow is required to read parquet for --method maxlfq. install.packages('arrow').")
  suppressMessages(library(limma))
  src <- file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "..")
  bm  <- file.path(src, "scripts", "build_maxlfq.R")
  if (file.exists(bm)) source(bm) else source("build_maxlfq.R")  # defines build_maxlfq()

  keep_runs <- meta$File.Name
  ml <- build_maxlfq(input, format = format, q_cutoff = q_cutoff,
                     eq_cutoff = eq_cutoff, pgq_cutoff = pgq_cutoff,
                     keep_runs = keep_runs)
  E         <- ml$E
  run_names <- colnames(E)
  genes     <- ml$genes
  descriptor <- ml$descriptor
  message(sprintf("[run_de] MaxLFQ matrix: %d proteins x %d runs (%.1f%% missing)",
                  nrow(E), ncol(E), 100 * mean(is.na(E))))
}

# ---- align metadata to the matrix columns -----------------------------------
meta <- meta[match(run_names, meta$File.Name), , drop = FALSE]
if (any(is.na(meta$Group)))
  stop("Some report columns have no metadata row. Report runs:\n  ",
       paste(run_names, collapse = "\n  "))
groups <- factor(meta$Group)

# ---- design matrix (~ 0 + groups [+ covariates]) ----------------------------
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

# rank check (DE-LIMP helpers.R guards against this before fitting)
if (qr(design)$rank < ncol(design))
  stop("Design matrix is not full rank — check for confounded covariates / empty groups.")

# ---- contrasts --------------------------------------------------------------
lvls <- levels(groups)
if (is.null(contrasts)) {
  ref  <- lvls[1]
  forms <- paste0(lvls[-1], "-", ref)        # all groups vs first level
} else {
  forms <- trimws(strsplit(contrasts, ",")[[1]])
}
message("[run_de] contrasts: ", paste(forms, collapse = ", "))

# ---- fit --------------------------------------------------------------------
if (method == "dpc") {
  fit <- limpa::dpcDE(y_protein, design, plot = FALSE)
} else {
  fit <- limma::lmFit(E, design)
}
fit <- limma::contrasts.fit(fit, limma::makeContrasts(contrasts = forms, levels = design))
fit <- limma::eBayes(fit)

# ---- write per-contrast results ---------------------------------------------
gene_cols <- intersect(c("Genes", "Protein.Names"), names(genes))
ann <- if (length(gene_cols)) genes[, c(if ("Protein.Group" %in% names(genes)) "Protein.Group", gene_cols), drop = FALSE] else NULL

# Expression matrix (proteins x samples, log2) — feeds figures (PCA/heatmap) and
# the report; mirrors DE-LIMP's Expression_Matrix.csv export.
expr_df <- data.frame(Protein.Group = rownames(E), check.names = FALSE)
if (!is.null(ann)) expr_df <- merge(expr_df, ann, by = "Protein.Group", all.x = TRUE, sort = FALSE)
expr_df <- merge(expr_df, data.frame(Protein.Group = rownames(E), E, check.names = FALSE),
                 by = "Protein.Group", all.x = TRUE, sort = FALSE)
utils::write.csv(expr_df, file.path(outdir, "Expression_Matrix.csv"), row.names = FALSE)
message(sprintf("[run_de] Expression_Matrix.csv: %d proteins x %d samples", nrow(E), ncol(E)))

all_sig <- list()
for (cn in forms) {
  tt <- limma::topTable(fit, coef = cn, number = Inf, adjust.method = "BH")
  tt$Protein.Group <- rownames(tt)
  if (!is.null(ann)) tt <- merge(tt, ann, by = "Protein.Group", all.x = TRUE, sort = FALSE)
  tt <- tt[order(tt$adj.P.Val), ]
  fn <- file.path(outdir, sprintf("DE_%s_%s.csv", method, make.names(cn)))
  utils::write.csv(tt, fn, row.names = FALSE)
  sig <- subset(tt, !is.na(adj.P.Val) & adj.P.Val < adjp_thr & abs(logFC) >= logfc_thr)
  all_sig[[cn]] <- nrow(sig)
  message(sprintf("[run_de] %-20s  %d proteins, %d significant (adj.P<%.2g, |logFC|>=%.2g) -> %s",
                  cn, nrow(tt), nrow(sig), adjp_thr, logfc_thr, basename(fn)))
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
  if (method == "maxlfq") sprintf("Normalization : quantile (limma::normalizeBetweenArrays)") else
                          sprintf("Normalization : DPC-CN (applied within dpcCN before dpcQuant)"),
  sprintf("Design        : ~ 0 + %s", paste(formula_parts, collapse = " + ")),
  sprintf("Contrasts     : %s", paste(forms, collapse = ", ")),
  sprintf("Thresholds    : adj.P.Val < %.3g and |logFC| >= %.3g", adjp_thr, logfc_thr),
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
  method = method, q_cutoff = q_cutoff, eq_cutoff = eq_cutoff, pgq_cutoff = pgq_cutoff,
  logfc = logfc_thr, adjp = adjp_thr,
  design = paste0("~ 0 + ", paste(formula_parts, collapse = " + ")),
  contrasts = forms, n_samples = nrow(meta), groups = as.list(table(groups)),
  significant_per_contrast = all_sig,
  R_version = as.character(getRversion()),
  packages = list(limpa = pkg_ver("limpa"), limma = pkg_ver("limma"),
                  arrow = pkg_ver("arrow"), dplyr = pkg_ver("dplyr"), tidyr = pkg_ver("tidyr")),
  input = normalizePath(input, mustWork = FALSE), metadata = normalizePath(meta_path, mustWork = FALSE)
)
writeLines(jsonlite_or_manual(prov), file.path(outdir, "de_provenance.json"))

message("\n[run_de] done. Results + provenance (methods.txt, sessionInfo.txt, de_provenance.json) in ",
        normalizePath(outdir))
