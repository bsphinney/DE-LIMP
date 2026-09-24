#!/usr/bin/env Rscript
# =============================================================================
# make_figures.R  --  Publication-quality proteomics figures for the report.
#
# Produces the standard figures a proteomics expert expects, from the DE output:
#   volcano_<contrast>.png   per-contrast volcano (top genes labelled)
#   pvalue_<contrast>.png    raw p-value distribution (calibration check)
#   pca.png                  sample PCA, coloured by group
#   heatmap_top.png          top differential proteins, z-scored, group-annotated
#   qc_protein_counts.png    proteins quantified per sample (loading/QC check) -- only when the
#                            matrix has missing values; a DPC/limpa matrix is complete, so it is
#                            skipped and qc_detected_vs_inferred.png is the depth view instead
#   qc_detected_vs_inferred.png  detected vs DPC-inferred proteins per sample (DPC runs)
#   violin_top_<contrast>.png    the top proteins of each contrast as violins, one point per
#                                run, each point marked measured vs inferred/missing
#
# Inputs:  --de-dir <output/tables>  (DE_*.csv + Expression_Matrix.csv from run_de.R;
#                                     Detection_Matrix.csv + de_provenance.json if present)
#          --conditions conditions.csv   --outdir output/figures
#          [--adjp 0.05] [--logfc 1] [--top 50] [--violin-top 8]
# Writes the PNGs + figures.json (a list of {file, type, caption}) for the report.
# Each figure is wrapped in tryCatch so one failure never blocks the others.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
getone <- function(flag, d = NULL) { i <- which(args == flag); if (!length(i)) d else args[i + 1] }
de_dir   <- getone("--de-dir", "output/tables")
cond_path<- getone("--conditions")
outdir   <- getone("--outdir", "output/figures")
adjp_thr <- as.numeric(getone("--adjp", "0.05"))
# Reference line only -- drawn and labelled on the volcano, never used to call
# significance. See the note in run_de.R: significance is the BH adjusted p-value alone.
logfc_ref<- as.numeric(getone("--logfc", "1"))
top_n    <- as.integer(getone("--top", "50"))
violin_n <- as.integer(getone("--violin-top", "8"))   # proteins per violin_top_<contrast>.png
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

# Group colours: ONE definition shared by the PCA, the heatmap annotation and the violins,
# so a group is the same colour in every figure. pheatmap's default hue wheel spaced 10
# groups so closely that Old_JPH3 and Old_Kv21 came out as the same pink. These were picked
# with an OKLab colour-difference check over ALL pairs (groups sit side by side in any
# order): up to 6 groups stay apart for colour-blind readers (deltaE >= 10.7), all 10 stay
# apart for normal vision (deltaE >= 15), and none is close to the two detection-status
# colours below. Groups are always labelled too, so colour is never the only cue.
GROUP_PAL <- c("#2a78d6", "#fb6ca0", "#78aa36", "#642ea4", "#2bc4ea",
               "#83382e", "#a54d9f", "#c65e0b", "#a491fe", "#144d6e")
group_colours <- function(lvls) {
  n <- length(lvls)
  cols <- GROUP_PAL[seq_len(min(n, length(GROUP_PAL)))]
  if (n > length(GROUP_PAL)) {
    message("[figures] ", n, " groups: colours past the ", length(GROUP_PAL),
            "th cannot all be told apart; rely on the group labels")
    cols <- c(cols, grDevices::hcl.colors(n - length(GROUP_PAL), "Dark 3"))
  }
  stats::setNames(cols, lvls)
}

# Detection-status colours: ONE definition for qc_detected_vs_inferred.png and the violins,
# so "teal = measured in that run, amber = not measured" means the same thing report-wide.
STATUS_DETECTED <- "#12866f"
STATUS_ABSENT   <- "#e8a33d"
# Direction colours: the volcano's up/down, reused for the violins' direction arrow.
DIR_COLS <- c(Up = "#d6604d", Down = "#4393c3")

gene_label <- function(df) {
  # Must be one value PER ROW: a bare NA makes ifelse() return a single element, and the
  # caller's rownames(H) <- ... then dies with "dimnames [1] not equal to array extent".
  # run_de.R omits Genes whenever the report carries no annotation, so this is reachable.
  g <- if ("Genes" %in% names(df)) df$Genes else rep(NA_character_, nrow(df))
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
    d$sig <- ifelse(d$adj.P.Val < adjp_thr,
                    ifelse(d$logFC > 0, "Up", "Down"), "NS")
    d$lab <- gene_label(d)
    nlogp <- -log10(pmax(d$adj.P.Val, .Machine$double.xmin))
    d$nlogp <- nlogp
    top <- d[d$sig != "NS", ]; top <- top[order(top$adj.P.Val), ]; top <- head(top, 15)
    p <- ggplot(d, aes(logFC, nlogp, color = sig)) +
      geom_point(alpha = 0.6, size = 1.4) +
      scale_color_manual(values = c(DIR_COLS, NS = "grey75"), name = NULL) +
      geom_vline(xintercept = c(-logfc_ref, logfc_ref), linetype = "dashed", color = "grey50") +
      geom_hline(yintercept = -log10(adjp_thr), linetype = "dashed", color = "grey50") +
      # Label the reference lines in the MARGIN, as ticks on a secondary axis, never inside
      # the panel: labels at the top of the lines collided with the repelled gene labels,
      # and labels tucked inward at their foot overprinted each other ("2-fold2-fold")
      # whenever the lines were close. check.overlap drops a label rather than overprint.
      scale_x_continuous(sec.axis = dup_axis(name = NULL, breaks = c(-logfc_ref, logfc_ref),
                                             labels = sprintf("%.3g×", 2^c(-logfc_ref, logfc_ref)))) +
      scale_y_continuous(sec.axis = dup_axis(name = NULL, breaks = -log10(adjp_thr),
                                             labels = sprintf("adj.P %.2g", adjp_thr))) +
      guides(x.sec = guide_axis(check.overlap = TRUE)) +
      # Two short lines: one long line ran past the device edge and was clipped.
      labs(title = paste0("Volcano — ", ct),
           subtitle = sprintf("%d up, %d down at adj.P < %.2g (BH)\nDashed lines: %.3g-fold and adj.P %.2g (fold is a reference, not a cutoff)",
                              sum(d$sig == "Up"), sum(d$sig == "Down"), adjp_thr, 2^logfc_ref, adjp_thr),
           x = "log2 fold change", y = "-log10 adjusted p-value") + THEME +
      theme(plot.title.position = "plot",
            axis.text.x.top = element_text(color = "grey40", size = 9),
            axis.text.y.right = element_text(color = "grey40", size = 9),
            axis.ticks.x.top = element_line(color = "grey50"),
            axis.ticks.y.right = element_line(color = "grey50"))
    if (has_repel && nrow(top)) p <- p + ggrepel::geom_text_repel(
      data = top, aes(label = lab), size = 3, max.overlaps = 20, show.legend = FALSE)
    fn <- file.path(outdir, sprintf("volcano_%s.png", make.names(ct)))
    ggsave(fn, p, width = 7, height = 6, dpi = 200)
    add_fig(fn, "volcano", sprintf(paste("Volcano plot for %s: log2 fold change vs significance.",
      "Coloured points are significant at adj.P < %.2g (Benjamini-Hochberg); no fold-change",
      "filter is applied. Vertical dashed lines mark %.3g-fold (labelled on the top axis) for",
      "reference only, so a coloured point inside them is a confidently measured small change,",
      "not an error."),
      ct, adjp_thr, 2^logfc_ref))

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

# ---- helpers for the top-protein violins ------------------------------------
or_else <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x)) || identical(x, "")) y else x

# de_provenance.json is how run_de.R's pipeline describes itself (architectural rule 1:
# never hardcode which pipeline ran). jsonlite is installed by setup.sh, but only ggplot2
# is guaranteed, so fall back to reading the few fields used here.
read_provenance <- function(path) {
  if (!file.exists(path)) return(list())
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    p <- tryCatch(jsonlite::fromJSON(path, simplifyVector = TRUE), error = function(e) NULL)
    if (is.list(p)) return(p)
  }
  txt <- paste(readLines(path, warn = FALSE), collapse = " ")
  str <- function(key) {
    m <- regmatches(txt, regexec(sprintf('"%s"\\s*:\\s*"([^"]*)"', key), txt))[[1]]
    if (length(m) == 2) m[2] else NULL
  }
  arr <- regmatches(txt, regexec('"contrasts"\\s*:\\s*\\[([^]]*)\\]', txt))[[1]]
  list(pipeline_id = str("pipeline_id"), method = str("method"), rollup_method = str("rollup_method"),
       contrasts = if (length(arr) == 2) gsub('"', "", regmatches(arr[2], gregexpr('"[^"]*"', arr[2]))[[1]]),
       detection_matrix = list(zero_means = str("zero_means")))
}

# What a 0 in Detection_Matrix.csv means depends on the quantification, and the pipeline
# records which in de_provenance.json. This is the ONE place that turns that into words:
#   a value exists although no precursor was seen -> "Inferred" (DPC-Quant; drawn hollow)
#   no value exists at all                        -> "Missing"  (MaxLFQ; nothing to draw)
# An explicit detection_matrix.zero_means wins; otherwise pipeline_id decides. With no
# record the data decide (a zero cell that still holds a value can only have been
# inferred), and the wording says the pipeline was not recorded rather than guessing it.
detection_vocab <- function(prov, M, D) {
  zm  <- tolower(or_else(prov$detection_matrix$zero_means, ""))
  pid <- or_else(prov$pipeline_id, or_else(prov$method, ""))
  src <- sub("\\s*\\(.*$", "", or_else(prov$rollup_method, ""))   # "DPC-Quant (Detection ...)" -> "DPC-Quant"
  absent <- if (zm %in% c("inferred", "missing")) c(inferred = "Inferred", missing = "Missing")[[zm]] else
            if (identical(pid, "maxlfq")) "Missing" else if (identical(pid, "dpc")) "Inferred" else NA
  recorded <- !is.na(absent)
  if (!recorded) {
    r <- intersect(rownames(M), rownames(D)); cc <- intersect(colnames(M), colnames(D))
    z <- D[r, cc, drop = FALSE] == 0 & !is.na(M[r, cc, drop = FALSE])
    absent <- if (any(z, na.rm = TRUE)) "Inferred" else "Missing"
  }
  how <- if (absent == "Missing") "not quantified in that run, so there is no value to draw" else
         sprintf("no precursor observed in that run; value from %s",
                 if (nzchar(src)) src else "the quantification model")
  if (!recorded) how <- paste0(how, " (pipeline not recorded in de_provenance.json)")
  list(absent = absent, how = how, source = if (nzchar(src)) src else "the quantification model")
}

# The DE file name carries make.names(contrast) ("Old_JPH3-Old_IgG" -> "Old_JPH3.Old_IgG"),
# which cannot be split back reliably once group names contain dots or underscores. Take
# the formula run_de.R actually used from de_provenance.json; failing that, accept the one
# A-B pair of known groups whose make.names() is the file's contrast. Never guess past that.
contrast_groups <- function(ct, forms, lvls) {
  form <- forms[make.names(forms) == ct]
  if (length(form) != 1) {
    form <- character(0)
    for (a in lvls) for (b in lvls)
      if (a != b && make.names(paste0(a, "-", b)) == ct) form <- c(form, paste0(a, "-", b))
    if (length(form) != 1) return(NULL)
  }
  # Split at the first top-level minus: the left side is compared against the right.
  ch <- strsplit(form, "")[[1]]; depth <- 0; cut <- NA
  for (i in seq_along(ch)) {
    if (ch[i] == "(") depth <- depth + 1 else if (ch[i] == ")") depth <- depth - 1 else
    if (ch[i] == "-" && depth == 0 && i > 1) { cut <- i; break }
  }
  if (is.na(cut)) return(NULL)
  toks <- function(s) { t <- regmatches(s, gregexpr("[A-Za-z.][A-Za-z0-9._]*", s))[[1]]; unique(t[t %in% lvls]) }
  test <- toks(substr(form, 1, cut - 1)); ref <- toks(substr(form, cut + 1, nchar(form)))
  if (!length(test) || !length(ref)) return(NULL)
  simple <- length(test) == 1 && length(ref) == 1
  list(test = test, ref = ref, simple = simple, label = if (simple) paste(test, "vs", ref) else form)
}

# Wrap a long group name at its own separators (_ . - space) so axis labels never collide;
# a single unbroken run longer than the width is cut into pieces.
wrap_label <- function(x, width = 12) vapply(x, function(s) {
  if (nchar(s) <= width) return(s)
  parts <- regmatches(s, gregexpr("[^_. -]+[_. -]*", s))[[1]]
  parts <- unlist(lapply(parts, function(p) if (nchar(p) <= width) p else
    substring(p, seq(1, nchar(p), width), pmin(seq(width, nchar(p) + width - 1, width), nchar(p)))))
  lines <- character(0); cur <- ""
  for (p in parts) if (nzchar(cur) && nchar(cur) + nchar(p) > width) { lines <- c(lines, cur); cur <- p } else cur <- paste0(cur, p)
  paste(c(lines, cur), collapse = "\n")
}, character(1), USE.NAMES = FALSE)

fmt_p <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), formatC(p, format = "fg", digits = 2))
tint  <- function(col, f) grDevices::rgb(t(grDevices::col2rgb(col) / 255 * (1 - f) + f))

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
  # Only when the matrix HAS missing values (e.g. MaxLFQ). A DPC/limpa matrix is complete by
  # construction -- every protein gets a value in every run -- so every bar is identical and the
  # plot says nothing (a 30-run limpa report shipped 30 bars of 6,112 and had to explain them
  # away). Then the detected-vs-inferred plot below is the per-sample depth view instead.
  tryCatch({
    cnt <- data.frame(Sample = colnames(M), n = colSums(!is.na(M)),
                      Group = if (!is.null(grp)) grp else "all")
    if (length(unique(cnt$n)) <= 1) {
      message("[figures] proteins-per-sample plot skipped: the matrix is complete (",
              cnt$n[1], " proteins in every sample), so every bar would be identical",
              if (file.exists(file.path(de_dir, "QC_detected_vs_inferred.csv")))
                "; qc_detected_vs_inferred.png shows per-sample depth" else "")
    } else {
      p <- ggplot(cnt, aes(reorder(Sample, n), n, fill = Group)) +
        geom_col() + coord_flip() +
        labs(title = "Proteins quantified per sample",
             x = NULL, y = "proteins (non-missing)") + THEME
      fn <- file.path(outdir, "qc_protein_counts.png")
      ggsave(fn, p, width = 7, height = max(3, 0.3 * ncol(M) + 1), dpi = 200)
      add_fig(fn, "qc", "Proteins quantified per sample — a loading/QC check. Large differences between samples (or systematic differences between groups) flag uneven input or sample-quality problems.")
    }
  }, error = function(e) message("[figures] QC counts failed: ", e$message))

  # ---- detected vs inferred (the QC view that actually works after DPC) ----
  # The plot above counts non-missing cells, but a DPC matrix is complete by
  # construction, so every bar comes out identical and real depth differences
  # are invisible. run_de.R writes QC_detected_vs_inferred.csv for this reason;
  # plot it when present. Mirrors DE-LIMP's Data Completeness panel.
  tryCatch({
    qcf <- file.path(de_dir, "QC_detected_vs_inferred.csv")
    if (file.exists(qcf)) {
      q <- utils::read.csv(qcf, stringsAsFactors = FALSE, check.names = FALSE)
      q <- q[order(q$Detected), ]
      long <- rbind(
        data.frame(Sample = q$Sample, n = q$Detected, Kind = "Detected"),
        data.frame(Sample = q$Sample, n = q$Inferred, Kind = "Inferred"))
      long$Sample <- factor(long$Sample, levels = q$Sample)
      long$Kind   <- factor(long$Kind, levels = c("Detected", "Inferred"))
      p2 <- ggplot2::ggplot(long, ggplot2::aes(n, Sample, fill = Kind)) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::scale_fill_manual(values = c(Detected = STATUS_DETECTED, Inferred = STATUS_ABSENT)) +
        ggplot2::labs(title = "Detected vs inferred proteins per sample",
                      subtitle = "Inferred values come from the DPC detection model, not measurement",
                      x = "proteins", y = NULL, fill = NULL) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(legend.position = "top")
      fn2 <- file.path(outdir, "qc_detected_vs_inferred.png")
      ggsave(fn2, p2, width = 9, height = max(3, 0.34 * nrow(q) + 1.6), dpi = 200)
      add_fig(fn2, "qc", sprintf("Detected vs inferred proteins per sample (%.0f%%-%.0f%% detected). Detected means at least one precursor was actually observed in that run; inferred means the value came from the DPC detection-probability model. Samples with a large inferred fraction contribute weaker evidence, and fold-changes for proteins inferred in one whole group should be read as detection events rather than magnitudes.", min(q$PctDetected), max(q$PctDetected)))
    }
  }, error = function(e) message("[figures] detected/inferred QC failed: ", e$message))

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
      if (!is.null(grp)) p <- p + scale_color_manual(values = group_colours(levels(grp)))
      if (has_repel) p <- p + ggrepel::geom_text_repel(size = 3, show.legend = FALSE)
      fn <- file.path(outdir, "pca.png")
      ggsave(fn, p, width = 7, height = 5.5, dpi = 200)
      add_fig(fn, "pca", "Principal-component analysis of samples (top 2 PCs). Replicates of the same group should cluster; clear separation between groups indicates a strong global difference, while an outlier sample stands apart.")
    }
  }, error = function(e) message("[figures] PCA failed: ", e$message))

  # ---- heatmap of top differential proteins ----
  tryCatch({
    # Rank by strength of evidence (smallest adjusted p across contrasts). Selecting by
    # matrix row order instead would make "top differential proteins" an arbitrary slice
    # of the significant set -- which matters more now that no fold-change filter thins it.
    best_p <- numeric(0)
    for (f in de_files) {
      d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
      d <- d[is.finite(d$adj.P.Val) & is.finite(d$logFC), ]
      d <- d[d$adj.P.Val < adjp_thr, , drop = FALSE]
      if (!nrow(d)) next
      v <- stats::setNames(d$adj.P.Val, d$Protein.Group)
      shared <- intersect(names(v), names(best_p))
      if (length(shared)) best_p[shared] <- pmin(best_p[shared], v[shared])
      new <- setdiff(names(v), names(best_p))
      if (length(new)) best_p <- c(best_p, v[new])
    }
    sig_ids <- names(sort(best_p))              # most significant first
    pick <- intersect(sig_ids, rownames(Mi))    # intersect() preserves sig_ids' order
    if (length(pick) < 5) {  # fall back to most-variable proteins
      v <- apply(Mi, 1, var); pick <- names(sort(v, decreasing = TRUE))[seq_len(min(top_n, nrow(Mi)))]
    } else pick <- head(pick, top_n)
    H <- Mi[pick, , drop = FALSE]
    lab <- gene_label(em[match(rownames(H), em$Protein.Group), , drop = FALSE])
    lab <- ifelse(is.na(lab), rownames(H), lab)
    # Semicolon-joined protein groups (e.g. "H2ac12;H2ac13;H2ac15;...") overflow the
    # device and clip the title and sample columns. Keep the first symbol, cap length.
    lab <- vapply(strsplit(lab, ";"), function(x) x[1], character(1))
    lab <- ifelse(nchar(lab) > 16, paste0(substr(lab, 1, 15), "~"), lab)
    rownames(H) <- make.unique(lab)
    fn <- file.path(outdir, "heatmap_top.png")
    ann <- if (!is.null(grp)) data.frame(Group = grp, row.names = colnames(H)) else NA
    if (has_pheat) {
      # pheatmap manages its own device; passing filename= avoids the clipped output
      # produced by wrapping it in png()/dev.off().
      pheatmap::pheatmap(H, scale = "row", annotation_col = ann,
                         annotation_colors = if (!is.null(grp)) list(Group = group_colours(levels(grp))) else NA,
                         show_rownames = nrow(H) <= 60, show_colnames = TRUE,
                         fontsize_row = 8, fontsize_col = 9,
                         main = sprintf("Top %d differential proteins (row z-score)", nrow(H)),
                         color = grDevices::colorRampPalette(c("#4393c3", "white", "#d6604d"))(100),
                         width = 11, height = 11, filename = fn)
    } else {  # base-R fallback heatmap
      grDevices::png(fn, width = 1700, height = 2200, res = 200)
      stats::heatmap(t(scale(t(H))), col = grDevices::colorRampPalette(c("#4393c3","white","#d6604d"))(100),
                     margins = c(8, 8), main = sprintf("Top %d differential proteins", nrow(H)))
      grDevices::dev.off()
    }
    add_fig(fn, "heatmap", sprintf("Heatmap of the top %d differential proteins (row z-scored log2 abundance), samples annotated by group. Reveals which proteins drive the group separation and whether replicates behave consistently.", nrow(H)))
  }, error = function(e) message("[figures] heatmap failed: ", e$message))

  # ---- top-protein violins, one figure per contrast ----
  # The volcano says WHICH proteins changed; this says what the change rests on. Each
  # panel is one protein, the contrast's two groups as violins, one point per run, and
  # each point marked measured or not (DE-LIMP's expression-grid violin, server_viz.R).
  # That mark is the point of the figure: DPC-Quant gives every protein a value in every
  # run, so a large fold change can rest on values that were never measured, and a group
  # with no measured run at all turns the fold change into a detection event rather than a
  # magnitude -- invisible in the DE table. Other groups are left out on purpose: the
  # heatmap already shows the top proteins across every group, and adding them here would
  # bury the two groups the contrast is about.
  prov <- read_provenance(file.path(de_dir, "de_provenance.json"))
  det_f <- file.path(de_dir, "Detection_Matrix.csv")
  D <- if (file.exists(det_f)) tryCatch({
    dm <- utils::read.csv(det_f, stringsAsFactors = FALSE, check.names = FALSE)
    dmat <- as.matrix(dm[, setdiff(names(dm), c("Protein.Group", "Genes", "Protein.Names")), drop = FALSE])
    storage.mode(dmat) <- "double"; rownames(dmat) <- dm$Protein.Group
    dmat
  }, error = function(e) { message("[figures] Detection_Matrix.csv unreadable (", e$message,
                                   "); violins drawn without detection status"); NULL })
  vocab <- if (!is.null(D)) detection_vocab(prov, M, D) else NULL
  INK <- "#1f1f1d"; MUTED <- "#6b6a66"; EVENT_INK <- "#8a5300"

  for (f in de_files) {
    ct <- contrast_of(f)
    tryCatch({
      if (is.null(grp)) stop("no conditions.csv group for every sample")
      cg <- contrast_groups(ct, or_else(prov$contrasts, character(0)), levels(grp))
      if (is.null(cg)) stop("could not identify this contrast's groups from de_provenance.json or its name")
      gx <- c(cg$ref, cg$test)             # reference on the left, so "higher" reads as rising
      # Layout: ~1.3 in per group so labels never collide, as many columns as fit in ~11 in.
      pw <- max(2.6, 1.3 * length(gx))
      glab <- wrap_label(gx, width = max(10, floor(13 * pw / length(gx))))
      d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
      d <- d[is.finite(d$logFC) & is.finite(d$adj.P.Val) & d$Protein.Group %in% rownames(M), , drop = FALSE]
      if (!nrow(d)) stop("no tested protein is in Expression_Matrix.csv")
      d <- utils::head(d[order(d$adj.P.Val, -abs(d$logFC)), , drop = FALSE], violin_n)
      k <- nrow(d)
      ncv <- max(1, min(k, floor(11 / pw))); nr <- ceiling(k / ncv); fig_w <- pw * ncv + 0.7

      # panel title: first gene symbol, else the first accession; disambiguate repeats
      lab <- vapply(strsplit(gene_label(d), ";"), `[`, "", 1)
      lab <- ifelse(nchar(lab) > 18, paste0(substr(lab, 1, 17), "~"), lab)
      dup <- duplicated(lab) | duplicated(lab, fromLast = TRUE)
      lab[dup] <- paste0(lab[dup], " (", vapply(strsplit(d$Protein.Group[dup], ";"), `[`, "", 1), ")")
      lab <- make.unique(lab)

      runs <- stats::setNames(lapply(gx, function(g) intersect(colnames(M)[grp == g], colnames(M))), gx)
      long <- do.call(rbind, lapply(seq_len(k), function(i) do.call(rbind, lapply(gx, function(g) {
        r <- runs[[g]]; pid <- d$Protein.Group[i]
        nobs <- if (!is.null(D) && pid %in% rownames(D)) D[pid, match(r, colnames(D))] else rep(NA_real_, length(r))
        data.frame(panel = rep(lab[i], length(r)), group = rep(g, length(r)), value = unname(M[pid, r]),
                   nobs = unname(as.numeric(nobs)), stringsAsFactors = FALSE)
      }))))
      long$panel <- factor(long$panel, levels = lab)
      long$x <- match(long$group, gx)
      long$status <- if (is.null(D)) "Not recorded" else
        ifelse(is.na(long$nobs), "Not recorded", ifelse(long$nobs > 0, "Detected", vocab$absent))

      # per protein x group: runs, measured, values, mean -- drive the labels under each violin
      cnt <- do.call(rbind, lapply(split(long, list(long$panel, long$group), drop = TRUE), function(z) data.frame(
        panel = z$panel[1], group = z$group[1], x = z$x[1], n_runs = nrow(z),
        n_meas = sum(z$nobs > 0, na.rm = TRUE), n_rec = sum(!is.na(z$nobs)), n_val = sum(!is.na(z$value)),
        mean = if (any(!is.na(z$value))) mean(z$value, na.rm = TRUE) else NA_real_)))
      cnt$event <- !is.null(D) & cnt$n_runs > 0 & cnt$n_rec == cnt$n_runs & cnt$n_meas == 0
      # Label only the exceptions: "3/3 measured" under every violin is noise that hides
      # the one group that matters.
      cnt$lab <- if (is.null(D)) ifelse(cnt$n_val < cnt$n_runs, sprintf("%d/%d quantified", cnt$n_val, cnt$n_runs), "") else
        ifelse(cnt$n_meas < cnt$n_runs, paste0(sprintf("%d/%d measured", cnt$n_meas, cnt$n_runs),
                                               ifelse(cnt$event, paste0("\nall ", tolower(vocab$absent)), "")), "")

      pts <- long[!is.na(long$value), , drop = FALSE]
      n_grp <- max(table(long$group)) / k                 # runs per group
      # Spread each group's points evenly across the violin in a fixed shuffled order, not
      # random jitter: DPC can infer the SAME value for every run of a group, and random
      # jitter then stacks three hollow points into what looks like one.
      w <- min(0.24, 0.12 + 0.01 * n_grp)
      set.seed(1)
      pts$xj <- pts$x + stats::ave(pts$value, pts$panel, pts$group, FUN = function(v)
        if (length(v) < 2) 0 else sample(seq(-w, w, length.out = length(v))))
      psize <- if (n_grp <= 6) 2.7 else if (n_grp <= 15) 2.1 else 1.5

      # stats label per panel: the model's estimate, not the plotted means
      sig <- d$adj.P.Val < adjp_thr
      st <- data.frame(panel = factor(lab, levels = lab), x = (1 + length(gx)) / 2,
                       lab = sprintf("log2FC %+.2f   adj.P %s%s", d$logFC, fmt_p(d$adj.P.Val), ifelse(sig, "", "  (n.s.)")),
                       sig = sig)

      gcol <- group_colours(levels(grp))[gx]
      p <- ggplot()
      for (g in gx) {
        vg <- pts[pts$group == g, , drop = FALSE]
        if (nrow(vg)) p <- p + geom_violin(data = vg, aes(x = x, y = value, group = x), fill = tint(gcol[[g]], 0.74),
                                           colour = gcol[[g]], linewidth = 0.45, width = 0.74,
                                           trim = TRUE, scale = "width")
      }
      mb <- cnt[!is.na(cnt$mean), , drop = FALSE]
      p <- p + geom_segment(data = mb, aes(x = x - 0.22, xend = x + 0.22, y = mean, yend = mean),
                            colour = INK, linewidth = 0.9, lineend = "round")
      if (cg$simple) {                                  # arrow reference mean -> compared mean
        a <- merge(mb[mb$group == cg$ref, c("panel", "mean")], mb[mb$group == cg$test, c("panel", "mean")],
                   by = "panel", suffixes = c("_ref", "_test"))
        if (nrow(a)) {
          a$dir <- ifelse(a$mean_test >= a$mean_ref, "Up", "Down")
          for (dd in unique(a$dir))              # one layer per direction: a constant colour each
            p <- p + geom_segment(data = a[a$dir == dd, ], aes(x = 1.3, xend = 1.7, y = mean_ref, yend = mean_test),
                                  colour = DIR_COLS[[dd]], linewidth = 0.55, alpha = 0.9,
                                  arrow = grid::arrow(length = grid::unit(0.06, "in"), type = "closed"))
        }
      }
      st_keys <- intersect(c("Detected", "Inferred", "Missing", "Not recorded"), unique(pts$status))
      st_lab <- c(Detected = "Measured in that run (precursors observed)",
                  Inferred = sprintf("Inferred: not measured; value from %s", or_else(vocab$source, "the model")),
                  Missing = "Missing (not quantified)", `Not recorded` = "Detection status not recorded")
      p <- p + geom_point(data = pts, aes(x = xj, y = value, fill = status, colour = status, shape = status, stroke = status),
                          size = psize) +
        scale_fill_manual(NULL, values = c(Detected = STATUS_DETECTED, Inferred = "white", Missing = "white",
                                           `Not recorded` = "#5f5e5a"), breaks = st_keys, labels = st_lab[st_keys]) +
        scale_colour_manual(NULL, values = c(Detected = "white", Inferred = STATUS_ABSENT, Missing = STATUS_ABSENT,
                                             `Not recorded` = "white"), breaks = st_keys, labels = st_lab[st_keys]) +
        scale_shape_manual(NULL, values = c(Detected = 21, Inferred = 21, Missing = 21, `Not recorded` = 21),
                           breaks = st_keys, labels = st_lab[st_keys]) +
        scale_discrete_manual("stroke", name = NULL, values = c(Detected = 0.45, Inferred = 1.3, Missing = 1.3,
                              `Not recorded` = 0.45), breaks = st_keys, labels = st_lab[st_keys])
      # No detection record -> no status legend: one neutral point style, nothing invented.
      p <- p + if (is.null(D)) guides(fill = "none", colour = "none", shape = "none", stroke = "none") else
        guides(fill = guide_legend(override.aes = list(size = 3.6)))
      p <- p +
        geom_text(data = st[st$sig, ], aes(x = x, y = Inf, label = lab), vjust = 1.5, size = 3.15, colour = INK) +
        geom_text(data = st[!st$sig, ], aes(x = x, y = Inf, label = lab), vjust = 1.5, size = 3.15, colour = MUTED) +
        geom_text(data = cnt[!cnt$event & nzchar(cnt$lab), ], aes(x = x, y = -Inf, label = lab), vjust = -0.5, size = 2.75,
                  colour = MUTED, lineheight = 0.9) +
        geom_text(data = cnt[cnt$event, ], aes(x = x, y = -Inf, label = lab), vjust = -0.35, size = 2.75,
                  colour = EVENT_INK, fontface = "bold", lineheight = 0.9) +
        facet_wrap(~panel, ncol = ncv, scales = "free_y") +
        scale_x_continuous(breaks = seq_along(gx), labels = glab, limits = c(0.4, length(gx) + 0.6),
                           expand = c(0, 0)) +
        scale_y_continuous(expand = expansion(mult = c(if (any(cnt$event)) 0.34 else if (any(nzchar(cnt$lab))) 0.2 else 0.08, 0.2)))

      # ---- words: title, subtitle, caption (all from the data + provenance) ----
      nsig <- sum(sig); n_abs <- sum(pts$status %in% c("Inferred", "Missing")); n_missing <- sum(is.na(long$value))
      ev <- cnt[cnt$event, , drop = FALSE]; n_ev <- length(unique(ev$panel))
      rank_line <- paste0("Ranked by adjusted p-value (BH), then |log2FC|; ",
        if (nsig == k) sprintf("all %d significant at adj.P < %.2g.", k, adjp_thr) else
        if (nsig > 0) sprintf("%d of %d significant at adj.P < %.2g, the rest (n.s.) shown as context, not findings.", nsig, k, adjp_thr) else
        sprintf("none reaches adj.P < %.2g, so these are shown for inspection only, not as findings.", adjp_thr))
      mark_line <- paste0("Bar = group mean", if (cg$simple) sprintf("; arrow = %s mean to %s mean", cg$ref, cg$test), ".",
                          if (!is.null(D) && any(nzchar(cnt$lab))) " A label under a violin counts its measured runs." else "")
      status_line <- if (is.null(D))
        "Detection status was not recorded for this run (no Detection_Matrix.csv), so measured and inferred values cannot be told apart here." else
        if (vocab$absent == "Inferred")
          paste0(sprintf("Filled = measured in that run; hollow = inferred (%s). ", vocab$how),
                 if (n_abs) sprintf("%d of %d values shown are inferred.", n_abs, nrow(pts)) else
                 sprintf("All %d values shown were measured.", nrow(pts))) else
          paste0(sprintf("Filled = measured in that run. Missing values (%s) are not drawn", vocab$how),
                 if (n_missing) sprintf(": %d of %d runs here.", n_missing, nrow(long)) else
                 sprintf("; none among the %d runs here.", nrow(long)))
      event_line <- if (n_ev) {
        items <- sprintf("%s in %s", ev$panel, ev$group)
        if (length(items) > 4) items <- c(items[1:4], sprintf("%d more", length(items) - 4))
        sprintf("Never measured (amber): %s. There the fold change is a detection event, not a measured magnitude.",
                paste(items, collapse = ", "))
      }
      wrapw <- round(13 * fig_w)                         # ~characters per subtitle line at this width
      subtitle <- paste(unlist(lapply(c(rank_line, status_line, mark_line, event_line), strwrap, width = wrapw)), collapse = "\n")
      n_sub <- length(strsplit(subtitle, "\n")[[1]])
      message("[figures] violin ", ct, " subtitle: ", gsub("\n", " ", subtitle))
      message("[figures] violin ", ct, " status legend: ", if (is.null(D)) "none" else paste(st_keys, collapse = ", "))

      # spaces around a formula's operators give a long title somewhere to wrap
      title <- paste(strwrap(sprintf("Top %d proteins: %s", k, if (cg$simple) cg$label else
                                     gsub("\\s*([-+*/])\\s*", " \\1 ", cg$label)), width = round(7.5 * fig_w)),
                     collapse = "\n")
      n_tit <- length(strsplit(title, "\n")[[1]])
      p <- p + labs(title = title, subtitle = subtitle,
                    x = NULL, y = "log2 intensity") +
        theme_bw(base_size = 12) +
        theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_line(colour = "#ecebe7", linewidth = 0.35),
              panel.border = element_rect(colour = "#d6d5cf", fill = NA, linewidth = 0.5),
              strip.background = element_rect(fill = "#f3f2ee", colour = "#d6d5cf", linewidth = 0.5),
              strip.text = element_text(face = "bold", size = 12, colour = INK, margin = margin(4, 4, 4, 4)),
              axis.text.x = element_text(size = 9.5, colour = "#3a3936", lineheight = 0.9),
              axis.text.y = element_text(size = 9, colour = MUTED), axis.ticks.x = element_blank(),
              axis.title.y = element_text(size = 10.5, colour = "#3a3936"),
              plot.title = element_text(face = "bold", size = 15, colour = INK),
              plot.subtitle = element_text(size = 10.5, colour = "#52514e", lineheight = 1.15, margin = margin(b = 4)),
              plot.title.position = "plot",
              legend.position = "top", legend.justification = "left", legend.text = element_text(size = 10),
              legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(4, "pt"),
              plot.margin = margin(10, 14, 8, 10))

      n_xl <- max(lengths(strsplit(glab, "\n")))
      fn <- file.path(outdir, sprintf("violin_top_%s.png", make.names(ct)))
      suppressWarnings(ggsave(fn, p, width = fig_w, dpi = 200,
                              height = 2.75 * nr + 0.5 + 0.3 * n_tit + 0.2 * n_sub + 0.17 * (n_xl - 1) +
                                       if (is.null(D)) 0 else 0.35))
      cap_status <- if (is.null(D))
        "Detection status was not recorded for this run (no Detection_Matrix.csv), so points are not marked measured or inferred; do not read them as all measured." else
        if (vocab$absent == "Inferred")
          sprintf(paste("Filled teal points were measured in that run (at least one precursor observed); hollow amber points",
                        "are inferred values: %s. %s Where not every run of a group was measured,",
                        "a label under its violin counts the measured runs. A group with no measured run (amber 'all inferred') makes the fold change a detection",
                        "event -- the protein is seen in one group and not the other -- not a measured magnitude; confirm such",
                        "a protein before building on the size of its change."), vocab$how,
                  if (n_abs) sprintf("%d of %d values shown are inferred.", n_abs, nrow(pts)) else
                  sprintf("All %d values shown were measured.", nrow(pts))) else
          sprintf(paste("Filled teal points were measured in that run. Missing values (%s) are not drawn; %d of %d runs here",
                        "are missing; a label under a violin counts its measured runs. A group with no measured run",
                        "(amber 'all missing') makes any fold change a detection event, not a measured magnitude."),
                  vocab$how, n_missing, nrow(long))
      add_fig(fn, "violin", paste(sprintf(paste(
        "Top %d proteins for %s, ranked by adjusted p-value (BH) then |log2 fold change| -- %d significant at adj.P < %.2g.",
        "Each panel is one protein: log2 intensity, one point per run, the reference group (%s) on the left.",
        "Violins show each group's spread (with few replicates the outline is only a guide; the points are the data);",
        "the black bar is the group mean%s.",
        "The label gives the model's log2 fold change and adjusted p-value (n.s. = not significant)."),
        k, cg$label, nsig, adjp_thr, paste(cg$ref, collapse = ", "),
        if (cg$simple) sprintf(" and the arrow runs from the %s mean to the %s mean (red = higher, blue = lower)",
                               cg$ref, cg$test) else ""), cap_status))
    }, error = function(e) message("[figures] violin ", ct, " skipped: ", e$message))
  }
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
