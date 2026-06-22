#!/usr/bin/env Rscript
# =============================================================================
# make_figures.R  --  Publication-quality proteomics figures for the report.
#
# Produces the standard figures a proteomics expert expects, from the DE output:
#   volcano_<contrast>.png   per-contrast volcano (top genes labelled)
#   pvalue_<contrast>.png    raw p-value distribution (calibration check)
#   pca.png                  sample PCA, coloured by group
#   heatmap_top.png          top differential proteins, z-scored, group-annotated
#   qc_protein_counts.png    proteins quantified per sample (loading/QC check)
#
# Inputs:  --de-dir <output/tables>  (DE_*.csv + Expression_Matrix.csv from run_de.R)
#          --conditions conditions.csv   --outdir output/figures
#          [--adjp 0.05] [--logfc 1] [--top 50]
# Writes the PNGs + figures.json (a list of {file, type, caption}) for the report.
# Each figure is wrapped in tryCatch so one failure never blocks the others.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
getone <- function(flag, d = NULL) { i <- which(args == flag); if (!length(i)) d else args[i + 1] }
de_dir   <- getone("--de-dir", "output/tables")
cond_path<- getone("--conditions")
outdir   <- getone("--outdir", "output/figures")
adjp_thr <- as.numeric(getone("--adjp", "0.05"))
logfc_thr<- as.numeric(getone("--logfc", "1"))
top_n    <- as.integer(getone("--top", "50"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

suppressWarnings(suppressMessages({
  has_ggplot <- requireNamespace("ggplot2", quietly = TRUE)
  has_repel  <- requireNamespace("ggrepel", quietly = TRUE)
  has_pheat  <- requireNamespace("pheatmap", quietly = TRUE)
}))
if (!has_ggplot) stop("ggplot2 is required for figures. Re-run setup.sh (installs r-ggplot2, r-ggrepel, r-pheatmap).")
library(ggplot2)

figs <- list()
add_fig <- function(file, type, caption) figs[[length(figs) + 1]] <<- list(file = basename(file), type = type, caption = caption)
THEME <- theme_bw(base_size = 13) + theme(panel.grid.minor = element_blank(),
                                          plot.title = element_text(face = "bold"))
gene_label <- function(df) {
  g <- if ("Genes" %in% names(df)) df$Genes else NA
  ifelse(is.na(g) | g == "", df$Protein.Group, g)
}

de_files <- list.files(de_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
contrast_of <- function(f) sub("\\.csv$", "", sub("^DE_[^_]+_", "", basename(f)))

# ---- volcano + p-value distribution, per contrast ---------------------------
for (f in de_files) {
  ct <- contrast_of(f)
  tryCatch({
    d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    d <- d[is.finite(d$logFC) & is.finite(d$adj.P.Val), ]
    d$sig <- ifelse(d$adj.P.Val < adjp_thr & abs(d$logFC) >= logfc_thr,
                    ifelse(d$logFC > 0, "Up", "Down"), "NS")
    d$lab <- gene_label(d)
    nlogp <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))
    d$nlogp <- nlogp
    top <- d[d$sig != "NS", ]; top <- top[order(top$adj.P.Val), ]; top <- head(top, 15)
    p <- ggplot(d, aes(logFC, nlogp, color = sig)) +
      geom_point(alpha = 0.6, size = 1.4) +
      scale_color_manual(values = c(Up = "#d6604d", Down = "#4393c3", NS = "grey75"), name = NULL) +
      geom_vline(xintercept = c(-logfc_thr, logfc_thr), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(adjp_thr), linetype = "dashed", color = "grey50") +
      labs(title = paste0("Volcano — ", ct),
           subtitle = sprintf("%d up, %d down (adj.P<%.2g, |logFC|>=%.2g)",
                              sum(d$sig == "Up"), sum(d$sig == "Down"), adjp_thr, logfc_thr),
           x = "log2 fold change", y = "-log10 adjusted p-value") + THEME
    if (has_repel && nrow(top)) p <- p + ggrepel::geom_text_repel(
      data = top, aes(label = lab), size = 3, max.overlaps = 20, show.legend = FALSE)
    fn <- file.path(outdir, sprintf("volcano_%s.png", make.names(ct)))
    ggsave(fn, p, width = 7, height = 6, dpi = 200)
    add_fig(fn, "volcano", sprintf("Volcano plot for %s: log2 fold change vs significance; coloured points pass adj.P<%.2g and |logFC|>=%.2g.", ct, adjp_thr, logfc_thr))

    if ("P.Value" %in% names(d)) {
      pp <- ggplot(d, aes(P.Value)) +
        geom_histogram(bins = 40, fill = "#4393c3", color = "white") +
        labs(title = paste0("p-value distribution — ", ct),
             x = "raw p-value", y = "proteins") + THEME
      fn2 <- file.path(outdir, sprintf("pvalue_%s.png", make.names(ct)))
      ggsave(fn2, pp, width = 6, height = 4.5, dpi = 200)
      add_fig(fn2, "pvalue", sprintf("Raw p-value distribution for %s. A peak near 0 over a flat background indicates real signal; a skew toward 1 or a spike mid-range suggests model/QC issues.", ct))
    }
  }, error = function(e) message("[figures] volcano/pvalue ", ct, " failed: ", e$message))
}

# ---- expression-matrix-based figures (PCA, heatmap, QC) ---------------------
em_path <- file.path(de_dir, "Expression_Matrix.csv")
meta <- if (!is.null(cond_path) && file.exists(cond_path))
  utils::read.csv(cond_path, stringsAsFactors = FALSE, check.names = FALSE) else NULL

if (file.exists(em_path)) {
  em <- utils::read.csv(em_path, stringsAsFactors = FALSE, check.names = FALSE)
  idcols <- intersect(c("Protein.Group", "Genes", "Protein.Names"), names(em))
  sample_cols <- setdiff(names(em), idcols)
  M <- as.matrix(em[, sample_cols, drop = FALSE]); rownames(M) <- em$Protein.Group
  storage.mode(M) <- "double"
  grp <- NULL
  if (!is.null(meta)) {
    g <- meta$Group[match(colnames(M), meta$File.Name)]
    if (all(!is.na(g))) grp <- factor(g)
  }

  # ---- QC: proteins quantified per sample ----
  tryCatch({
    cnt <- data.frame(Sample = colnames(M), n = colSums(!is.na(M)),
                      Group = if (!is.null(grp)) grp else "all")
    p <- ggplot(cnt, aes(reorder(Sample, n), n, fill = Group)) +
      geom_col() + coord_flip() +
      labs(title = "Proteins quantified per sample",
           x = NULL, y = "proteins (non-missing)") + THEME
    fn <- file.path(outdir, "qc_protein_counts.png")
    ggsave(fn, p, width = 7, height = max(3, 0.3 * ncol(M) + 1), dpi = 200)
    add_fig(fn, "qc", "Proteins quantified per sample — a loading/QC check. Large differences between samples (or systematic differences between groups) flag uneven input or sample-quality problems.")
  }, error = function(e) message("[figures] QC counts failed: ", e$message))

  # complete-ish matrix for PCA/heatmap: keep proteins seen in all samples; if too
  # few, mean-impute per protein (PCA/heatmap need no NAs).
  complete <- M[rowSums(is.na(M)) == 0, , drop = FALSE]
  Mi <- if (nrow(complete) >= 10) complete else {
    imp <- M; rm <- rowMeans(imp, na.rm = TRUE)
    imp[is.na(imp)] <- rm[row(imp)[is.na(imp)]]; imp[rowSums(is.na(imp)) == 0, , drop = FALSE]
  }

  # ---- PCA ----
  tryCatch({
    if (ncol(Mi) >= 3 && nrow(Mi) >= 5) {
      pc <- prcomp(t(Mi), scale. = TRUE)
      ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
      pdf <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], Sample = colnames(Mi),
                        Group = if (!is.null(grp)) grp else "all")
      p <- ggplot(pdf, aes(PC1, PC2, color = Group, label = Sample)) +
        geom_point(size = 3) +
        labs(title = "Sample PCA",
             x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) + THEME
      if (has_repel) p <- p + ggrepel::geom_text_repel(size = 3, show.legend = FALSE)
      fn <- file.path(outdir, "pca.png")
      ggsave(fn, p, width = 7, height = 5.5, dpi = 200)
      add_fig(fn, "pca", "Principal-component analysis of samples (top 2 PCs). Replicates of the same group should cluster; clear separation between groups indicates a strong global difference, while an outlier sample stands apart.")
    }
  }, error = function(e) message("[figures] PCA failed: ", e$message))

  # ---- heatmap of top differential proteins ----
  tryCatch({
    sig_ids <- character(0)
    for (f in de_files) {
      d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
      d <- d[is.finite(d$adj.P.Val) & is.finite(d$logFC), ]
      s <- d$Protein.Group[d$adj.P.Val < adjp_thr & abs(d$logFC) >= logfc_thr]
      sig_ids <- union(sig_ids, s)
    }
    pick <- intersect(sig_ids, rownames(Mi))
    if (length(pick) < 5) {  # fall back to most-variable proteins
      v <- apply(Mi, 1, var); pick <- names(sort(v, decreasing = TRUE))[seq_len(min(top_n, nrow(Mi)))]
    } else pick <- head(pick[order(match(pick, rownames(Mi)))], top_n)
    H <- Mi[pick, , drop = FALSE]
    lab <- gene_label(em[match(rownames(H), em$Protein.Group), , drop = FALSE])
    rownames(H) <- make.unique(ifelse(is.na(lab), rownames(H), lab))
    fn <- file.path(outdir, "heatmap_top.png")
    ann <- if (!is.null(grp)) data.frame(Group = grp, row.names = colnames(H)) else NA
    if (has_pheat) {
      grDevices::png(fn, width = 1700, height = 2200, res = 200)
      pheatmap::pheatmap(H, scale = "row", annotation_col = ann,
                         show_rownames = nrow(H) <= 60, show_colnames = TRUE,
                         main = sprintf("Top %d differential proteins (row z-score)", nrow(H)),
                         color = grDevices::colorRampPalette(c("#4393c3", "white", "#d6604d"))(100))
      grDevices::dev.off()
    } else {  # base-R fallback heatmap
      grDevices::png(fn, width = 1700, height = 2200, res = 200)
      stats::heatmap(t(scale(t(H))), col = grDevices::colorRampPalette(c("#4393c3","white","#d6604d"))(100),
                     margins = c(8, 8), main = sprintf("Top %d differential proteins", nrow(H)))
      grDevices::dev.off()
    }
    add_fig(fn, "heatmap", sprintf("Heatmap of the top %d differential proteins (row z-scored log2 abundance), samples annotated by group. Reveals which proteins drive the group separation and whether replicates behave consistently.", nrow(H)))
  }, error = function(e) message("[figures] heatmap failed: ", e$message))
} else {
  message("[figures] no Expression_Matrix.csv in ", de_dir, " — skipping PCA/heatmap/QC")
}

# ---- write figures.json -----------------------------------------------------
to_json <- function(figs) {
  if (requireNamespace("jsonlite", quietly = TRUE))
    return(jsonlite::toJSON(figs, auto_unbox = TRUE, pretty = TRUE))
  items <- vapply(figs, function(f) sprintf('  {"file": "%s", "type": "%s", "caption": "%s"}',
                  f$file, f$type, gsub('"', '\\\\"', f$caption)), "")
  paste0("[\n", paste(items, collapse = ",\n"), "\n]")
}
writeLines(to_json(figs), file.path(outdir, "figures.json"))
cat(sprintf("[figures] wrote %d figure(s) + figures.json to %s\n", length(figs), normalizePath(outdir)))
