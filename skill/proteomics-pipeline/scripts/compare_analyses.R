#!/usr/bin/env Rscript
# =============================================================================
# compare_analyses.R  --  Compare two or more analyses of the SAME dataset
# (e.g. an original run and a re-analysis with different params/engine/version).
#
# This is a faithful port of DE-LIMP's Run Comparator core (R/server_comparator.R):
#   - normalize_protein_id()  : sp|P12345|GENE -> P12345, strip isoforms/groups
#   - classify_de()           : Up / Down / NS (NaN-safe), adj.P threshold
#   - the 3x3 concordance matrix on shared proteins
# Ported (not sourced) so the skill stays self-contained, the same way run_de.R
# and build_maxlfq.R are faithful ports of the DE-LIMP pipeline.
#
# For each contrast shared across analyses it reports, per analysis pair:
#   protein-universe overlap, the 3x3 Up/Down/NS concordance, direction
#   concordance on co-significant proteins, and logFC correlation on shared
#   proteins. Writes per-pair CSVs + a COMPARISON.md summary.
#
# Usage:
#   Rscript compare_analyses.R --out ./comparison --adjp 0.05 --logfc 0 \
#     --analysis "Original:/path/sessA/output/tables" \
#     --analysis "Reanalysis:/path/sessB/output/tables"
# Each path is a DE output dir holding DE_<method>_<contrast>.csv files.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
getall <- function(flag) { i <- which(args == flag); if (!length(i)) character(0) else args[i + 1] }
getone <- function(flag, default = NULL) { v <- getall(flag); if (!length(v)) default else v[1] }

out_dir  <- getone("--out", "comparison")
adjp_thr <- as.numeric(getone("--adjp", "0.05"))
logfc_thr<- as.numeric(getone("--logfc", "0"))
analyses_raw <- getall("--analysis")
if (length(analyses_raw) < 2)
  stop("Need >= 2 --analysis \"Label:/path/to/de_dir\" arguments to compare.")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- ported from DE-LIMP server_comparator.R --------------------------------
normalize_protein_id <- function(ids) {
  ids <- gsub(".*\\|([A-Z][0-9][A-Z0-9]{3}[0-9])(-[0-9]+)?\\|.*", "\\1", ids)
  ids <- gsub("-[0-9]+$", "", ids)
  ids <- gsub(";.*", "", ids)
  trimws(ids)
}
classify_de <- function(logfc, adjp, threshold = 0.05, lfc_min = 0) {
  ifelse(!is.finite(adjp) | adjp >= threshold, "NS",
         ifelse(is.finite(logfc) & logfc >  abs(lfc_min), "Up",
                ifelse(is.finite(logfc) & logfc < -abs(lfc_min), "Down", "NS")))
}

# ---- load analyses ----------------------------------------------------------
parse_spec <- function(s) {
  # "Label:/path" — split on the FIRST colon only (paths may contain none)
  i <- regexpr(":", s, fixed = TRUE)
  if (i < 0) stop("--analysis must be 'Label:/path/to/de_dir', got: ", s)
  list(label = substr(s, 1, i - 1), dir = substr(s, i + 1, nchar(s)))
}
specs <- lapply(analyses_raw, parse_spec)

# contrast key from a DE filename: DE_<method>_<contrast>.csv -> <contrast>
contrast_key <- function(fn) sub("\\.csv$", "", sub("^DE_[^_]+_", "", basename(fn)))

read_de_dir <- function(dir) {
  files <- list.files(dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
  if (!length(files)) stop("No DE_*.csv files in ", dir)
  out <- list()
  for (f in files) {
    df <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if (!"Protein.Group" %in% names(df)) {
      # tolerate a different id column name
      idc <- intersect(c("Protein.Group", "Protein.Ids", "Protein", "Protein.IDs"), names(df))
      if (!length(idc)) next
      names(df)[names(df) == idc[1]] <- "Protein.Group"
    }
    df$protein_id <- normalize_protein_id(df$Protein.Group)
    out[[contrast_key(f)]] <- df
  }
  out
}
loaded <- lapply(specs, function(s) list(label = s$label, de = read_de_dir(s$dir)))
labels <- vapply(loaded, function(x) x$label, "")

# ---- protein universe per analysis ------------------------------------------
universe_rows <- list()
for (a in loaded) {
  all_ids <- unique(unlist(lapply(a$de, function(df) df$protein_id)))
  n_sig <- length(unique(unlist(lapply(names(a$de), function(k) {
    df <- a$de[[k]]
    df$protein_id[classify_de(df$logFC, df$adj.P.Val, adjp_thr, logfc_thr) != "NS"]
  }))))
  universe_rows[[a$label]] <- data.frame(Analysis = a$label,
                                         n_proteins = length(all_ids),
                                         n_significant = n_sig)
}
universe_df <- do.call(rbind, universe_rows)
utils::write.csv(universe_df, file.path(out_dir, "protein_universe.csv"), row.names = FALSE)

# ---- pairwise concordance per shared contrast -------------------------------
pairs <- utils::combn(seq_along(loaded), 2, simplify = FALSE)
summary_rows <- list()
md <- c(sprintf("# Analysis comparison — %s", paste(labels, collapse = " vs ")),
        "",
        sprintf("Compared %d analyses of the same dataset. Significance: adj.P < %.3g, |logFC| > %.3g.",
                length(loaded), adjp_thr, logfc_thr),
        "",
        "## Protein universe", "",
        "| Analysis | Proteins | Significant |", "|---|---|---|")
for (i in seq_len(nrow(universe_df)))
  md <- c(md, sprintf("| %s | %d | %d |", universe_df$Analysis[i],
                      universe_df$n_proteins[i], universe_df$n_significant[i]))
md <- c(md, "", "## Pairwise concordance (shared contrasts)", "")

for (pr in pairs) {
  A <- loaded[[pr[1]]]; B <- loaded[[pr[2]]]
  shared_contrasts <- intersect(names(A$de), names(B$de))
  if (!length(shared_contrasts)) {
    md <- c(md, sprintf("- **%s vs %s**: no shared contrasts (%s | %s)",
                        A$label, B$label, paste(names(A$de), collapse=","),
                        paste(names(B$de), collapse=",")))
    next
  }
  for (ct in shared_contrasts) {
    da <- A$de[[ct]]; db <- B$de[[ct]]
    m <- merge(da[, c("protein_id", "logFC", "adj.P.Val")],
               db[, c("protein_id", "logFC", "adj.P.Val")],
               by = "protein_id", suffixes = c("_A", "_B"))
    if (!nrow(m)) next
    m$status_A <- classify_de(m$logFC_A, m$adj.P.Val_A, adjp_thr, logfc_thr)
    m$status_B <- classify_de(m$logFC_B, m$adj.P.Val_B, adjp_thr, logfc_thr)
    lev <- c("Up", "Down", "NS")
    conc <- table(factor(m$status_A, levels = lev), factor(m$status_B, levels = lev))

    # co-significant proteins (sig in both) and their direction agreement
    sig_both <- m[m$status_A != "NS" & m$status_B != "NS", ]
    dir_agree <- if (nrow(sig_both)) mean(sig_both$status_A == sig_both$status_B) else NA_real_
    # logFC correlation on all shared
    lfc_cor <- tryCatch(stats::cor(m$logFC_A, m$logFC_B, use = "complete.obs"),
                        error = function(e) NA_real_)
    sig_A <- sum(m$status_A != "NS"); sig_B <- sum(m$status_B != "NS")
    sig_shared <- sum(m$status_A != "NS" & m$status_B != "NS")

    tag <- sprintf("%s_vs_%s__%s", make.names(A$label), make.names(B$label), make.names(ct))
    cm <- as.data.frame.matrix(conc); cm <- cbind(`A\\B` = rownames(cm), cm)
    utils::write.csv(cm, file.path(out_dir, sprintf("concordance_%s.csv", tag)), row.names = FALSE)
    utils::write.csv(m[order(m$adj.P.Val_A), ],
                     file.path(out_dir, sprintf("merged_%s.csv", tag)), row.names = FALSE)

    summary_rows[[tag]] <- data.frame(
      contrast = ct, analysis_A = A$label, analysis_B = B$label,
      n_shared = nrow(m), sig_A = sig_A, sig_B = sig_B, sig_both = sig_shared,
      direction_concordance = round(dir_agree, 3), logFC_correlation = round(lfc_cor, 3))

    md <- c(md,
      sprintf("### %s — %s vs %s", ct, A$label, B$label),
      sprintf("- shared proteins: %d; significant: %s=%d, %s=%d, both=%d",
              nrow(m), A$label, sig_A, B$label, sig_B, sig_shared),
      sprintf("- direction concordance (co-significant): %s",
              ifelse(is.na(dir_agree), "n/a", sprintf("%.1f%%", 100*dir_agree))),
      sprintf("- logFC correlation: %s",
              ifelse(is.na(lfc_cor), "n/a", sprintf("%.3f", lfc_cor))),
      sprintf("- 3x3 table -> `concordance_%s.csv`; merged -> `merged_%s.csv`", tag, tag),
      "")
  }
}

if (length(summary_rows)) {
  summ <- do.call(rbind, summary_rows)
  utils::write.csv(summ, file.path(out_dir, "concordance_summary.csv"), row.names = FALSE)
}
md <- c(md, "---",
        "_Comparison generated by the proteomics-pipeline skill "
        , "(faithful port of DE-LIMP's Run Comparator concordance)._")
writeLines(md, file.path(out_dir, "COMPARISON.md"))

cat(sprintf("[compare] %d analyses, %d pair(s); outputs in %s\n",
            length(loaded), length(pairs), normalizePath(out_dir)))
