#!/usr/bin/env Rscript
# =============================================================================
# make_figures.R  --  Publication-quality proteomics figures for the report.
#
# Produces the standard figures a proteomics expert expects, from the DE output:
#   volcano_<contrast>.png   per-contrast volcano (top genes labelled)
#   qc_pvalue_panel.png      raw p-value distributions, every contrast in one small-multiples
#                            panel (calibration check, for the report's appendix)
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
#          --adjp/--logfc are a FALLBACK: the cutoff the DE run recorded in de_provenance.json
#          wins, so the figures call significance exactly as the DE tables and the report do.
# Writes the PNGs + figures.json: {"figures": [{file, type, caption}],
#                                  "failed":  [{file, type, reason}]}
# Every PNG this script owns is deleted before drawing, and each figure is wrapped in
# tryCatch: one failure never blocks the others, and none is left stale or silently absent.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
getone <- function(flag, d = NULL) { i <- which(args == flag); if (!length(i)) d else args[i + 1] }
de_dir   <- getone("--de-dir", "output/tables")
cond_path<- getone("--conditions")
outdir   <- getone("--outdir", "output/figures")
cli_adjp <- getone("--adjp")    # fallbacks only -- see "significance cutoff" below
cli_logfc<- getone("--logfc")
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
# A figure this script meant to draw and did not -- an error, or too little data -- is listed in
# figures.json "failed" with the reason, so the report can say so instead of silently lacking it.
failed <- list()
add_failed <- function(file, type, reason) {
  failed[[length(failed) + 1]] <<- list(file = basename(file), type = type, reason = reason)
  message("[figures] ", basename(file), " not drawn: ", reason)
}
# Every PNG this script owns is removed before drawing. Otherwise a figure that now fails, or is
# no longer drawn (the per-sample count plot for a complete matrix, the old per-contrast p-value
# histograms), survives from an earlier run and is embedded as if current. Only these names.
OWNED_FIGS <- paste0("^(volcano_.+|pvalue_.+|violin_top_.+|pca|heatmap_top|qc_protein_counts|",
                     "qc_detected_vs_inferred|qc_pvalue_panel)\\.png$")
old_figs <- list.files(outdir, pattern = OWNED_FIGS, full.names = TRUE)
if (length(old_figs) && all(file.remove(old_figs)))
  message(sprintf("[figures] removed %d figure(s) left by an earlier run in %s", length(old_figs), outdir))
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

# Short sample names for figure labels. Raw run names ("08132026__60SPD_DIA-LRS-124_S3-F4_1_23659")
# took ~40% of the heatmap and buried the PCA in text. ONE definition for every figure that
# names samples; the full <-> short mapping goes to sample_labels.csv so nothing is ambiguous.
#   1. a Label (or Sample) column in conditions.csv, when the user supplied one;
#   2. otherwise derived: drop the tokens every run shares at the start and at the end, keep
#      the shortest leading stretch of what remains that tells the runs apart, and keep the
#      last shared token when that would leave a bare number ("LRS-124", not "124");
#   3. any label that comes out empty or not unique falls back to the full run name.
short_sample_names <- function(runs, meta = NULL) {
  runs <- as.character(runs)
  lab_col <- if (!is.null(meta)) intersect(c("Label", "Sample"), names(meta))[1] else NA
  if (!is.na(lab_col)) {
    lab <- trimws(as.character(meta[[lab_col]][match(runs, meta$File.Name)]))
    src <- rep(paste("conditions.csv", lab_col), length(runs))
  } else {
    src <- rep("derived from the run name", length(runs))
    tk <- lapply(runs, function(r) {
      m <- gregexpr("[^_. -]+", r)[[1]]
      if (m[1] < 0) return(list(s = 1L, e = nchar(r), t = r))
      s <- as.integer(m); e <- s + attr(m, "match.length") - 1L
      list(s = s, e = e, t = substring(r, s, e))
    })
    n_tok <- vapply(tk, function(x) length(x$t), 1L)
    shared <- function(pick) length(unique(vapply(tk, pick, ""))) == 1
    a <- 0L                                   # tokens shared at the start (keep >= 1 per run)
    while (a < min(n_tok) - 1L && shared(function(x) x$t[a + 1L])) a <- a + 1L
    b <- 0L                                   # tokens shared at the end
    while (a + b < min(n_tok) - 1L && shared(function(x) x$t[length(x$t) - b])) b <- b + 1L
    for (L in seq_len(max(n_tok - a - b))) {
      last <- pmin(a + L, n_tok - b)
      lab <- vapply(seq_along(runs), function(k) substring(runs[k], tk[[k]]$s[a + 1L], tk[[k]]$e[last[k]]), "")
      if (!anyDuplicated(lab)) break
    }
    if (a > 0L && all(grepl("^[0-9]", lab)))  # a bare number says little: keep its shared label
      lab <- vapply(seq_along(runs), function(k) substring(runs[k], tk[[k]]$s[a], tk[[k]]$e[last[k]]), "")
  }
  bad <- is.na(lab) | !nzchar(lab) | duplicated(lab) | duplicated(lab, fromLast = TRUE)
  lab[bad] <- runs[bad]
  src[bad] <- "full run name (label empty or not unique)"
  data.frame(File.Name = runs, Label = make.unique(lab), Source = src, stringsAsFactors = FALSE)
}

# Rounded hull around a group's samples: a small circle of points around every sample, then
# the convex hull of all of them. Works for any n -- a circle for 1 sample, a capsule for 2, a
# smooth rounded hull for 3+ -- and assumes no distribution, unlike a 95% normal ellipse, which
# with n = 3 is unstable and enormous. Returned closed (last vertex = first).
rounded_hull <- function(x, y, r, n_arc = 32) {
  th <- seq(0, 2 * pi, length.out = n_arc + 1)[-1]
  px <- as.vector(outer(x, r * cos(th), "+")); py <- as.vector(outer(y, r * sin(th), "+"))
  h <- grDevices::chull(px, py)
  data.frame(x = px[c(h, h[1])], y = py[c(h, h[1])])
}

# Do the group names cross two factors, e.g. Old_JPH3 / Young_IgG = Age x Bait? Only when the
# "_" split is unambiguous: every name has the same number of "_" tokens, exactly one split
# point gives two factors with >= 2 levels each crossing (near-)fully (>= 75% of cells), one of
# them has exactly 2 levels (drawn as point shape + outline style) and the other fits the
# palette (drawn as colour). Returns per-level colour/style factors, or NULL.
two_factor_split <- function(lvls) {
  tok <- strsplit(lvls, "_", fixed = TRUE); k <- unique(lengths(tok))
  if (length(lvls) < 4 || length(k) != 1 || k < 2) return(NULL)
  hits <- list()
  for (p in seq_len(k - 1)) {
    a <- vapply(tok, function(t) paste(t[seq_len(p)], collapse = "_"), "")
    b <- vapply(tok, function(t) paste(t[(p + 1):k], collapse = "_"), "")
    na <- length(unique(a)); nb <- length(unique(b))
    if (na < 2 || nb < 2 || length(lvls) / (na * nb) < 0.75 || !(na == 2 || nb == 2)) next
    split <- if (nb == 2 && na != 2) list(colour = a, style = b) else list(colour = b, style = a)
    if (length(unique(split$colour)) > length(GROUP_PAL)) next
    hits[[length(hits) + 1]] <- lapply(split, stats::setNames, lvls)
  }
  if (length(hits) == 1) hits[[1]] else NULL
}

de_files <- list.files(de_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)
contrast_of <- function(f) sub("\\.csv$", "", sub("^DE_[^_]+_", "", basename(f)))

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
  num <- function(key) {
    m <- regmatches(txt, regexec(sprintf('"%s"\\s*:\\s*([-0-9.eE+]+)', key), txt))[[1]]
    if (length(m) == 2) as.numeric(m[2]) else NULL
  }
  arr <- regmatches(txt, regexec('"contrasts"\\s*:\\s*\\[([^]]*)\\]', txt))[[1]]
  list(pipeline_id = str("pipeline_id"), method = str("method"), rollup_method = str("rollup_method"),
       adjp = num("adjp"), logfc = num("logfc"),
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

# ---- significance cutoff: ONE definition -------------------------------------
# The DE run records the cutoff it applied (de_provenance.json: adjp, logfc). The figures must
# call significance exactly as the DE tables, the audit and the report do, so that record wins;
# --adjp/--logfc only stand in when a run left no record, and every caption that depends on the
# cutoff then says where it came from (architectural rule 2: a default never passes silently
# as the analysis's own value).
prov <- read_provenance(file.path(de_dir, "de_provenance.json"))
num1 <- function(x) { v <- suppressWarnings(as.numeric(x)); if (length(v) == 1 && is.finite(v)) v else NA_real_ }
resolve_cut <- function(recorded, cli, default, flag, what) {
  if (!is.na(num1(recorded))) {
    if (!is.null(cli) && !isTRUE(all.equal(num1(cli), num1(recorded))))
      message(sprintf("[figures] %s %s ignored: the DE run used %s %s (de_provenance.json)",
                      flag, cli, what, format(num1(recorded))))
    return(list(value = num1(recorded), src = "record"))
  }
  if (!is.na(num1(cli))) return(list(value = num1(cli), src = "cli"))
  list(value = default, src = "default")
}
.adjp  <- resolve_cut(prov$adjp,  cli_adjp,  0.05, "--adjp",  "adj.P <")
.logfc <- resolve_cut(prov$logfc, cli_logfc, 1,    "--logfc", "a reference line at |log2FC| =")
adjp_thr  <- .adjp$value
# Reference line only -- drawn and labelled on the volcano, never used to call significance.
# See the note in run_de.R: significance is the BH adjusted p-value alone.
logfc_ref <- .logfc$value
cut_part <- function(label, src) if (src == "cli") sprintf("%s from the command line", label) else
  sprintf("%s is make_figures.R's DEFAULT, not user-confirmed", label)
.parts <- c(if (.adjp$src != "record") cut_part(sprintf("adj.P < %.2g", adjp_thr), .adjp$src),
            if (.logfc$src != "record") cut_part(sprintf("the %.3g-fold reference line", 2^logfc_ref), .logfc$src))
CUT_NOTE <- if (!length(.parts)) "" else
  sprintf(" (The DE run recorded no cutoff in de_provenance.json: %s.)", paste(.parts, collapse = "; "))
message(sprintf("[figures] significance cutoff: adj.P < %s (%s); fold reference %s (%s)",
                format(adjp_thr), .adjp$src, format(logfc_ref), .logfc$src))

# A contrast's display label from the formula the DE run recorded ("A-B" -> "A vs B").
contrast_display <- function(ct) {
  form <- or_else(prov$contrasts, character(0)); form <- form[make.names(form) == ct]
  if (length(form) == 1 && lengths(regmatches(form, gregexpr("-", form))) == 1 && !grepl("[()+*/]", form))
    sub("\\s*-\\s*", " vs ", form) else ct
}

# ---- volcano per contrast; p-values collected for one panel -------------------
pv <- list()
for (f in de_files) {
  ct <- contrast_of(f)
  fn <- file.path(outdir, sprintf("volcano_%s.png", make.names(ct)))
  tryCatch({
    d <- utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    if ("P.Value" %in% names(d) && any(is.finite(d$P.Value)))
      pv[[ct]] <- data.frame(contrast = contrast_display(ct), P.Value = d$P.Value[is.finite(d$P.Value)])
    d <- d[is.finite(d$logFC) & is.finite(d$adj.P.Val), ]
    if (!nrow(d)) stop("no protein has a finite logFC and adj.P.Val")
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
           subtitle = sprintf("%d up, %d down at adj.P < %.2g (BH)%s\nDashed lines: %.3g-fold and adj.P %.2g (fold is a reference, not a cutoff)",
                              sum(d$sig == "Up"), sum(d$sig == "Down"), adjp_thr,
                              if (nzchar(CUT_NOTE)) " -- cutoff not from the DE record" else "",
                              2^logfc_ref, adjp_thr),
           x = "log2 fold change", y = "-log10 adjusted p-value") + THEME +
      theme(plot.title.position = "plot",
            axis.text.x.top = element_text(color = "grey40", size = 9),
            axis.text.y.right = element_text(color = "grey40", size = 9),
            axis.ticks.x.top = element_line(color = "grey50"),
            axis.ticks.y.right = element_line(color = "grey50"))
    if (has_repel && nrow(top)) p <- p + ggrepel::geom_text_repel(
      data = top, aes(label = lab), size = 3, max.overlaps = 20, show.legend = FALSE)
    ggsave(fn, p, width = 7, height = 6, dpi = 200)
    add_fig(fn, "volcano", paste0(sprintf(paste("Volcano plot for %s: log2 fold change vs significance.",
      "Coloured points are significant at adj.P < %.2g (Benjamini-Hochberg); no fold-change",
      "filter is applied. Vertical dashed lines mark %.3g-fold (labelled on the top axis) for",
      "reference only, so a coloured point inside them is a confidently estimated small change,",
      "not an error."),
      ct, adjp_thr, 2^logfc_ref), CUT_NOTE))
  }, error = function(e) add_failed(fn, "volcano", conditionMessage(e)))
}

# ---- raw p-value distributions: ONE small-multiples panel ----------------------
# A calibration check, not a finding: one panel for the report's appendix instead of a
# separate figure per contrast that the narrative then interleaves with the results.
fn_pv <- file.path(outdir, "qc_pvalue_panel.png")
if (length(pv)) tryCatch({
  pvd <- do.call(rbind, pv)
  pvd$contrast <- factor(pvd$contrast, levels = unique(pvd$contrast))
  nc <- length(levels(pvd$contrast)); ncol_pv <- min(4, nc); nrow_pv <- ceiling(nc / ncol_pv)
  pp <- ggplot(pvd, aes(P.Value)) +
    geom_histogram(breaks = seq(0, 1, by = 0.025), fill = "#4393c3", colour = "white", linewidth = 0.15) +
    facet_wrap(~contrast, ncol = ncol_pv, scales = "free_y", labeller = label_wrap_gen(width = 28)) +
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"), expand = expansion(mult = 0.01)) +
    labs(title = "Raw p-value distributions",
         subtitle = "One panel per contrast. Well calibrated: flat, with a peak at 0 when there is real signal.",
         x = "raw p-value", y = "proteins") +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
          plot.title.position = "plot", strip.background = element_rect(fill = "#f3f2ee", colour = "#d6d5cf"),
          strip.text = element_text(size = 9), plot.subtitle = element_text(colour = "#52514e"))
  ggsave(fn_pv, pp, width = 2.6 * ncol_pv + 0.8, height = 2.1 * nrow_pv + 1.1, dpi = 200)
  add_fig(fn_pv, "pvalue", sprintf(paste(
    "Raw p-value distributions for all %d contrast%s (appendix; a calibration check, not a result).",
    "A flat background with a peak near 0 indicates real signal on a well-behaved model; a skew",
    "toward 1, a hump mid-range or a U shape suggests model or QC problems for that contrast."),
    nc, if (nc > 1) "s" else ""))
}, error = function(e) add_failed(fn_pv, "pvalue", conditionMessage(e)))

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
  # short sample labels for every figure that names samples (see short_sample_names)
  sl <- short_sample_names(colnames(M), meta)
  slab <- stats::setNames(sl$Label, sl$File.Name)
  short_of <- function(x) ifelse(is.na(slab[x]), x, slab[x])
  utils::write.csv(sl, file.path(outdir, "sample_labels.csv"), row.names = FALSE)
  LABELS_NOTE <- " Samples are labelled with short names; sample_labels.csv maps each one to its run file."
  # Which values were measured -- read once, used by the QC plot, PCA, heatmap and violins.
  # (prov, the pipeline's own record, was read above with the significance cutoff.)
  det_f <- file.path(de_dir, "Detection_Matrix.csv")
  D <- if (file.exists(det_f)) tryCatch({
    dm <- utils::read.csv(det_f, stringsAsFactors = FALSE, check.names = FALSE)
    dmat <- as.matrix(dm[, setdiff(names(dm), c("Protein.Group", "Genes", "Protein.Names")), drop = FALSE])
    storage.mode(dmat) <- "double"; rownames(dmat) <- dm$Protein.Group
    dmat
  }, error = function(e) { message("[figures] Detection_Matrix.csv unreadable (", e$message,
                                   "); figures drawn without detection status"); NULL })
  vocab <- if (!is.null(D)) detection_vocab(prov, M, D) else NULL

  # ---- QC: proteins quantified per sample ----
  # Only when the matrix HAS missing values (e.g. MaxLFQ). A DPC/limpa matrix is complete by
  # construction -- every protein gets a value in every run -- so every bar is identical and the
  # plot says nothing (a 30-run limpa report shipped 30 bars of 6,112 and had to explain them
  # away). Then the detected-vs-inferred plot below is the per-sample depth view instead.
  tryCatch({
    cnt <- data.frame(Sample = short_of(colnames(M)), n = colSums(!is.na(M)),
                      Group = if (!is.null(grp)) grp else "all")
    if (length(unique(cnt$n)) <= 1) {
      message("[figures] proteins-per-sample plot skipped: the matrix is complete (",
              cnt$n[1], " proteins in every sample), so every bar would be identical",
              if (file.exists(file.path(de_dir, "QC_detected_vs_inferred.csv")))
                "; qc_detected_vs_inferred.png shows per-sample depth" else "")
      # Not a failure -- nothing to show -- so not listed in "failed". A copy from an earlier
      # run was already removed by the sweep at the top (Silva08172026's re-rendered report
      # had shown the identical-bar plot again).
    } else {
      p <- ggplot(cnt, aes(reorder(Sample, n), n, fill = Group)) +
        geom_col() + coord_flip() +
        labs(title = "Proteins quantified per sample",
             x = NULL, y = "proteins (non-missing)") + THEME
      fn <- file.path(outdir, "qc_protein_counts.png")
      ggsave(fn, p, width = 7, height = max(3, 0.3 * ncol(M) + 1), dpi = 200)
      add_fig(fn, "qc", paste0("Proteins quantified per sample — a loading/QC check. Large differences between samples (or systematic differences between groups) flag uneven input or sample-quality problems.", LABELS_NOTE))
    }
  }, error = function(e) add_failed("qc_protein_counts.png", "qc", conditionMessage(e)))

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
      # "RUN · Group", so a group whose runs all sit at the bottom shows at a glance.
      q$Sample <- if ("Group" %in% names(q) && !all(is.na(q$Group)))
        paste(short_of(q$Sample), q$Group, sep = " \u00b7 ") else short_of(q$Sample)
      inf_by <- if (!is.null(vocab)) vocab$source else or_else(sub("\\s*\\(.*$", "", or_else(prov$rollup_method, "")),
                                                               "the detection-probability model")
      long <- rbind(
        data.frame(Sample = q$Sample, n = q$Detected, Kind = "Detected"),
        data.frame(Sample = q$Sample, n = q$Inferred, Kind = "Inferred"))
      long$Sample <- factor(long$Sample, levels = q$Sample)
      long$Kind   <- factor(long$Kind, levels = c("Detected", "Inferred"))
      p2 <- ggplot2::ggplot(long, ggplot2::aes(n, Sample, fill = Kind)) +
        ggplot2::geom_col(width = 0.72) +
        ggplot2::scale_fill_manual(values = c(Detected = STATUS_DETECTED, Inferred = STATUS_ABSENT)) +
        ggplot2::labs(title = "Detected vs inferred proteins per sample",
                      subtitle = sprintf("Inferred values come from %s, not measurement", inf_by),
                      x = "proteins", y = NULL, fill = NULL) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(legend.position = "top")
      fn2 <- file.path(outdir, "qc_detected_vs_inferred.png")
      ggsave(fn2, p2, width = 9, height = max(3, 0.34 * nrow(q) + 1.6), dpi = 200)
      add_fig(fn2, "qc", sprintf("Detected vs inferred proteins per sample (%.0f%%-%.0f%% detected), each bar labelled run \u00b7 group. Detected means at least one precursor was actually observed in that run; inferred means the value came from %s. Samples with a large inferred fraction contribute weaker evidence, and fold-changes for proteins inferred in one whole group should be read as detection events rather than magnitudes.%s", min(q$PctDetected), max(q$PctDetected), inf_by, LABELS_NOTE))
    }
  }, error = function(e) add_failed("qc_detected_vs_inferred.png", "qc", conditionMessage(e)))

  # complete-ish matrix for PCA/heatmap: keep proteins seen in all samples; if too
  # few, mean-impute per protein (PCA/heatmap need no NAs).
  complete <- M[rowSums(is.na(M)) == 0, , drop = FALSE]
  mi_imputed <- nrow(complete) < 10
  Mi <- if (!mi_imputed) complete else {
    imp <- M; rm <- rowMeans(imp, na.rm = TRUE)
    imp[is.na(imp)] <- rm[row(imp)[is.na(imp)]]; imp[rowSums(is.na(imp)) == 0, , drop = FALSE]
  }

  # ---- PCA ----
  # Every group is circled with a rounded hull (rounded_hull) and named on the plot itself, so
  # nobody bounces to a legend, and the axes share one scale, so distances between samples are
  # true. When the group names cross two factors (two_factor_split), one is colour and the other
  # point shape + outline style -- never filled vs hollow, which the violins reserve for
  # measured vs inferred. The subtitle states how the PCA was computed, from the code below.
  tryCatch({
    Mp <- Mi[apply(Mi, 1, stats::var) > 0, , drop = FALSE]    # prcomp cannot scale a constant protein
    n_const <- nrow(Mi) - nrow(Mp)
    if (ncol(Mp) < 3 || nrow(Mp) < 5) {
      add_failed("pca.png", "pca", sprintf(
        "a PCA needs at least 3 samples and 5 varying proteins; this matrix has %d sample%s and %d protein%s",
        ncol(Mp), if (ncol(Mp) == 1) "" else "s", nrow(Mp), if (nrow(Mp) == 1) "" else "s"))
    } else {
      pc <- prcomp(t(Mp), scale. = TRUE)
      ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
      g  <- if (!is.null(grp)) grp else factor(rep("all samples", ncol(Mp)))
      pdf <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], Sample = short_of(colnames(Mp)),
                        Group = as.character(g), stringsAsFactors = FALSE)
      tf <- two_factor_split(levels(g))
      pdf$Colour <- if (is.null(tf)) pdf$Group else unname(tf$colour[pdf$Group])
      pdf$Style  <- if (is.null(tf)) "all" else unname(tf$style[pdf$Group])
      col_lv <- unique(if (is.null(tf)) levels(g) else tf$colour[levels(g)])
      sty_lv <- unique(if (is.null(tf)) "all" else tf$style[levels(g)])
      pal  <- group_colours(col_lv)
      dark <- stats::setNames(grDevices::rgb(t(grDevices::col2rgb(pal) / 255 * 0.62)), paste(col_lv, "text"))
      message("[figures] PCA encoding: ", if (is.null(tf)) "one colour per group" else
              sprintf("two factors -- colour = %s; shape/outline = %s", paste(col_lv, collapse = "/"),
                      paste(sty_lv, collapse = "/")))

      span <- max(diff(range(pdf$PC1)), diff(range(pdf$PC2)))
      hull <- do.call(rbind, lapply(split(pdf, pdf$Group), function(z) cbind(
        rounded_hull(z$PC1, z$PC2, 0.035 * span), Group = z$Group[1], Colour = z$Colour[1], Style = z$Style[1])))
      cen <- do.call(rbind, lapply(split(pdf, pdf$Group), function(z) data.frame(
        Group = z$Group[1], Colour = z$Colour[1], PC1 = mean(z$PC1), PC2 = mean(z$PC2))))
      # samples far from their own group's centroid (> 2x the median such distance) get named
      m <- match(pdf$Group, cen$Group)
      dist <- sqrt((pdf$PC1 - cen$PC1[m])^2 + (pdf$PC2 - cen$PC2[m])^2)
      med <- stats::median(dist[dist > 0])
      out <- which(is.finite(med) & dist > 2 * med)
      out <- utils::head(out[order(-dist[out])], 5)

      # Crossed design: join each colour level's centroids across the two style levels (e.g. a
      # bait's Old and Young groups), so a consistent shift from the second factor shows at once.
      pair <- NULL
      if (!is.null(tf)) {
        cen$Style <- unname(tf$style[cen$Group])
        a1 <- cen[cen$Style == sty_lv[1], c("Colour", "PC1", "PC2")]
        a2 <- cen[cen$Style == sty_lv[2], c("Colour", "PC1", "PC2")]
        pair <- merge(a1, a2, by = "Colour", suffixes = c("_1", "_2"))
      }
      shapes <- stats::setNames(c(21, 24)[seq_along(sty_lv)], sty_lv)
      ltys   <- stats::setNames(c("solid", "22")[seq_along(sty_lv)], sty_lv)
      p <- ggplot() +
        geom_polygon(data = hull, aes(x, y, group = Group, fill = Colour, colour = Colour, linetype = Style),
                     alpha = 0.13, linewidth = 0.55, show.legend = FALSE) +
        (if (!is.null(pair) && nrow(pair)) geom_segment(data = pair, aes(x = PC1_1, y = PC2_1, xend = PC1_2, yend = PC2_2,
                                                                         colour = Colour), linewidth = 0.5, alpha = 0.75)) +
        geom_point(data = cen, aes(PC1, PC2, colour = Colour), shape = 3, size = 2, stroke = 0.6, alpha = 0.6) +
        geom_point(data = pdf, aes(PC1, PC2, fill = Colour, shape = Style), colour = "white", size = 3.1, stroke = 0.5) +
        scale_fill_manual(values = pal, guide = "none") +
        scale_colour_manual(values = c(pal, dark, outlier = "#52514e"), guide = "none") +
        scale_linetype_manual(values = ltys, guide = "none") +
        scale_shape_manual(values = shapes, name = NULL,
                           labels = if (length(sty_lv) == 2) paste0(sty_lv, c(" (circle, solid outline)",
                                                                              " (triangle, dashed outline)")) else sty_lv) +
        guides(shape = if (length(sty_lv) == 2) guide_legend(override.aes = list(fill = "#8f8d87", size = 3.4)) else "none") +
        coord_fixed()
      labs_df <- rbind(
        data.frame(x = cen$PC1, y = cen$PC2, label = cen$Group, key = paste(cen$Colour, "text"),
                   face = "bold", sz = 3.7),
        data.frame(x = pdf$PC1[out], y = pdf$PC2[out], label = pdf$Sample[out], key = rep("outlier", length(out)),
                   face = rep("plain", length(out)), sz = rep(2.9, length(out))))
      if (has_repel) {
        # the samples go in as empty labels: ggrepel keeps the names off the points
        keep <- setdiff(seq_len(nrow(pdf)), out)       # (x[-integer(0)] would drop every sample)
        obst <- data.frame(x = pdf$PC1[keep], y = pdf$PC2[keep], label = rep("", length(keep)),
                           key = "outlier", face = "plain", sz = 1)
        p <- p + ggrepel::geom_text_repel(data = rbind(labs_df, obst),
                   aes(x, y, label = label, colour = key, fontface = face, size = sz),
                   bg.color = "white", bg.r = 0.12, box.padding = 0.45, point.padding = 0.3,
                   min.segment.length = 0.3, segment.colour = "#9a9892", segment.size = 0.3,
                   max.overlaps = Inf, seed = 1, show.legend = FALSE)
      } else {
        p <- p + geom_text(data = labs_df, aes(x, y, label = label, colour = key, fontface = face, size = sz),
                           vjust = -1.1, show.legend = FALSE)
      }
      p <- p + scale_size_identity()

      pct_inf <- NA_real_
      if (!is.null(D) && !is.null(vocab) && identical(vocab$absent, "Inferred")) {
        z <- D[match(rownames(Mp), rownames(D)), match(colnames(Mp), colnames(D)), drop = FALSE]
        if (any(!is.na(z))) pct_inf <- 100 * mean(z == 0, na.rm = TRUE)
      }
      basis <- sprintf("%s proteins %s, log2 intensities centred and scaled to unit variance per protein%s.",
                       format(nrow(Mp), big.mark = ","),
                       if (!mi_imputed) "with a value in every sample" else
                         "with any value (missing values set to the protein's mean: fewer than 10 had a value in every sample)",
                       if (n_const) sprintf("; %d constant protein%s left out", n_const, if (n_const > 1) "s" else "") else "")
      inf_line <- if (is.finite(pct_inf) && pct_inf > 0)
        sprintf("%.0f%% of these values are inferred by %s (no precursor observed in that run), not measured.",
                pct_inf, vocab$source)
      enc_line <- paste0("Rounded hulls circle each group's samples; + = group centroid",
                         if (!is.null(pair) && nrow(pair)) sprintf("; a line joins the %s and %s centroids of each colour",
                           sty_lv[1], sty_lv[2]) else "",
                         if (length(out)) "; named samples lie > 2x the median distance from their group's centroid." else ".")
      subtitle <- paste(unlist(lapply(c(basis, inf_line, enc_line), strwrap, width = 118)), collapse = "\n")
      message("[figures] PCA subtitle: ", gsub("\n", " ", subtitle))
      p <- p + labs(title = "Sample PCA", subtitle = subtitle,
                    x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
        theme_bw(base_size = 12) +
        theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#ecebe7", linewidth = 0.35),
              panel.border = element_rect(colour = "#d6d5cf", fill = NA, linewidth = 0.5),
              plot.title = element_text(face = "bold", size = 15), plot.title.position = "plot",
              plot.subtitle = element_text(size = 10, colour = "#52514e", lineheight = 1.15),
              legend.position = "top", legend.justification = "left", legend.text = element_text(size = 10),
              legend.margin = margin(0, 0, 0, 0), axis.title = element_text(size = 11, colour = "#3a3936"))

      xr <- range(hull$x); yr <- range(hull$y)
      asp <- min(1.3, max(0.45, diff(yr) / diff(xr)))
      n_sub <- length(strsplit(subtitle, "\n")[[1]])
      fn <- file.path(outdir, "pca.png")
      ggsave(fn, p, width = 8.4, height = 6.9 * asp + 0.75 + 0.19 * n_sub + if (length(sty_lv) == 2) 0.3 else 0,
             dpi = 200)
      enc_cap <- if (!is.null(tf)) sprintf(paste(" The group names cross two factors, so colour = %s and point shape +",
        "hull outline = %s (%s: circle, solid outline; %s: triangle, dashed outline); a thin line joins the %s and %s",
        "centroids of each colour, so a consistent shift between them shows as parallel lines."),
        paste(col_lv, collapse = " / "), paste(sty_lv, collapse = " / "), sty_lv[1], sty_lv[2], sty_lv[1], sty_lv[2]) else ""
      add_fig(fn, "pca", paste0(
        "Principal-component analysis of samples (PC1 vs PC2, drawn to equal scale so distances compare). ", basis,
        if (!is.null(inf_line)) paste0(" ", inf_line) else "",
        " Each group is circled by a rounded hull around its samples -- the actual spread of its replicates, not a",
        " statistical confidence region -- and named on the plot; + marks the group centroid.", enc_cap,
        " Replicates of a group should sit together; well-separated hulls mean a strong global difference between",
        " groups, and a named sample lies more than twice the median distance from its group's centroid, worth checking.",
        if (length(ve) >= 3) sprintf(" PC3 explains a further %.1f%%.", ve[3]) else "", LABELS_NOTE))
    }
  }, error = function(e) add_failed("pca.png", "pca", conditionMessage(e)))

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
    hm_inf <- ""
    if (!is.null(D) && !is.null(vocab) && identical(vocab$absent, "Inferred")) {
      z <- D[match(rownames(H), rownames(D)), match(colnames(H), colnames(D)), drop = FALSE]
      if (any(!is.na(z))) hm_inf <- sprintf(paste(
        " %.0f%% of the cells shown are inferred by %s (no precursor observed in that run), not measured.",
        "A uniform block of colour across a group can therefore be the model's estimates for proteins",
        "never seen there, not agreement between replicates; the violins mark which values were measured."),
        100 * mean(z == 0, na.rm = TRUE), vocab$source)
    }
    lab <- gene_label(em[match(rownames(H), em$Protein.Group), , drop = FALSE])
    lab <- ifelse(is.na(lab), rownames(H), lab)
    # Semicolon-joined protein groups (e.g. "H2ac12;H2ac13;H2ac15;...") overflow the
    # device and clip the title and sample columns. Keep the first symbol, cap length.
    lab <- vapply(strsplit(lab, ";"), function(x) x[1], character(1))
    lab <- ifelse(nchar(lab) > 16, paste0(substr(lab, 1, 15), "~"), lab)
    rownames(H) <- make.unique(lab)
    fn <- file.path(outdir, "heatmap_top.png")
    colnames(H) <- short_of(colnames(H))
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
    add_fig(fn, "heatmap", sprintf("Heatmap of the top %d differential proteins (row z-scored log2 abundance), samples annotated by group. Reveals which proteins drive the group separation and whether replicates behave consistently.%s%s%s", nrow(H), hm_inf, LABELS_NOTE, CUT_NOTE))
  }, error = function(e) add_failed("heatmap_top.png", "heatmap", conditionMessage(e)))

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
      # BH gives many proteins the same adj.P.Val; break ties by the raw p-value (the DE table's
      # own topTable order), not by |logFC|, which would promote large but noisy changes.
      ord <- if ("P.Value" %in% names(d)) order(d$adj.P.Val, d$P.Value) else order(d$adj.P.Val)
      d <- utils::head(d[ord, , drop = FALSE], violin_n)
      message("[figures] violin ", ct, " proteins: ", paste(d$Protein.Group, collapse = ", "))
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
        data.frame(panel = rep(lab[i], length(r)), group = rep(g, length(r)), run = r, value = unname(M[pid, r]),
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
      # A density needs data: a violin drawn through 3 points is a shape the data do not have.
      # Groups with <= 4 runs in a panel are shown as their points and a mean bar only.
      MIN_VIOLIN <- 5
      n_violin <- 0
      for (g in gx) {
        vg <- pts[pts$group == g, , drop = FALSE]
        vg <- vg[stats::ave(vg$value, vg$panel, FUN = length) >= MIN_VIOLIN, , drop = FALSE]
        if (nrow(vg)) {
          n_violin <- n_violin + 1
          p <- p + geom_violin(data = vg, aes(x = x, y = value, group = x), fill = tint(gcol[[g]], 0.74),
                               colour = gcol[[g]], linewidth = 0.45, width = 0.74, trim = TRUE, scale = "width")
        }
      }
      # Within-block (paired) contrast: join each block's runs across the two groups, so the
      # reader sees the within-block changes the model actually tested, not two clouds.
      bcol <- or_else(prov$block$column, NA_character_)
      bstr <- unlist(prov$block$contrast_structure)
      paired <- cg$simple && !is.na(bcol) && !is.null(meta) && bcol %in% names(meta) &&
        identical(unname(bstr[make.names(names(bstr)) == ct]), "within")
      if (paired) {
        pts$block <- as.character(meta[[bcol]][match(pts$run, meta$File.Name)])
        one_each <- all(table(pts$panel, pts$block, pts$group) <= 1)
        bl <- if (one_each) pts[, c("panel", "block", "x", "xj", "value")] else {
          m <- stats::aggregate(value ~ panel + block + x, data = pts, FUN = mean); m$xj <- m$x; m }
        bl <- bl[order(bl$panel, bl$block, bl$x), ]
        p <- p + geom_line(data = bl, aes(x = xj, y = value, group = interaction(panel, block)),
                           colour = "#a8a69f", linewidth = 0.35, alpha = 0.9)
        message(sprintf("[figures] violin %s: paired lines join each %s (%s)", ct, bcol,
                        if (one_each) "one run per group" else "block means per group"))
      }
      mb <- cnt[!is.na(cnt$mean), , drop = FALSE]
      p <- p + geom_segment(data = mb, aes(x = x - 0.22, xend = x + 0.22, y = mean, yend = mean),
                            colour = INK, linewidth = 0.9, lineend = "round")
      # The arrow is the MODEL's log2 fold change, drawn up (or down) from the reference mean, so
      # it always agrees in sign with the label. The difference of plain means can disagree once
      # the model weights runs, models inferred values or removes a block effect.
      if (cg$simple) {
        a <- merge(mb[mb$group == cg$ref, c("panel", "mean")],
                   data.frame(panel = factor(lab, levels = lab), lfc = d$logFC), by = "panel")
        if (nrow(a)) {
          a$end <- a$mean + a$lfc
          message("[figures] violin ", ct, " arrows (model log2FC): ",
                  paste(sprintf("%s %+.2f", a$panel, a$end - a$mean), collapse = ", "))   # the span as drawn
          a$dir <- ifelse(a$lfc >= 0, "Up", "Down")
          for (dd in unique(a$dir))              # one layer per direction: a constant colour each
            p <- p + geom_segment(data = a[a$dir == dd, ], aes(x = 1.3, xend = 1.7, y = mean, yend = end),
                                  colour = DIR_COLS[[dd]], linewidth = 0.55, alpha = 0.9,
                                  arrow = grid::arrow(length = grid::unit(0.06, "in"), type = "closed"))
        }
      }
      message(sprintf("[figures] violin %s drawn as: %s", ct,
                      if (n_violin == length(gx)) "violins" else if (n_violin) "violins + points"
                      else sprintf("points + mean bars (<= %d runs per group)", MIN_VIOLIN - 1)))
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
      rank_line <- paste0("Ranked by adjusted p-value (BH), ties by raw p-value; ",
        if (nsig == k) sprintf("all %d significant at adj.P < %.2g.", k, adjp_thr) else
        if (nsig > 0) sprintf("%d of %d significant at adj.P < %.2g, the rest (n.s.) shown as context, not findings.", nsig, k, adjp_thr) else
        sprintf("none reaches adj.P < %.2g, so these are shown for inspection only, not as findings.", adjp_thr))
      mark_line <- paste0("Bar = group mean",
                          if (cg$simple) sprintf("; arrow = the model's log2FC, drawn from the %s mean", cg$ref) else "",
                          if (paired) sprintf("; grey lines join each %s's runs", bcol) else "",
                          if (n_violin < length(gx)) sprintf("; groups with <= %d runs are shown as points (no density)", MIN_VIOLIN - 1) else "",
                          ".",
                          if (!is.null(D) && any(nzchar(cnt$lab))) " A label under a group counts its measured runs." else "")
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
                        "a label under that group counts its measured runs. A group with no measured run (amber 'all inferred') makes the fold change a detection",
                        "event -- the protein is seen in one group and not the other -- not a measured magnitude; confirm such",
                        "a protein before building on the size of its change."), vocab$how,
                  if (n_abs) sprintf("%d of %d values shown are inferred.", n_abs, nrow(pts)) else
                  sprintf("All %d values shown were measured.", nrow(pts))) else
          sprintf(paste("Filled teal points were measured in that run. Missing values (%s) are not drawn; %d of %d runs here",
                        "are missing; a label under a group counts its measured runs. A group with no measured run",
                        "(amber 'all missing') makes any fold change a detection event, not a measured magnitude."),
                  vocab$how, n_missing, nrow(long))
      add_fig(fn, "violin", paste0(paste(sprintf(paste(
        "Top %d proteins for %s, ranked by adjusted p-value (BH), ties broken by raw p-value -- %d significant at adj.P < %.2g.",
        "Each panel is one protein: log2 intensity, one point per run, the reference group (%s) on the left.",
        "%s The black bar is the group mean%s.%s",
        "The label gives the model's log2 fold change and adjusted p-value (n.s. = not significant)."),
        k, cg$label, nsig, adjp_thr, paste(cg$ref, collapse = ", "),
        if (n_violin == length(gx)) "Violins show each group's spread; the points are the data." else if (n_violin)
          sprintf("A violin shows the spread of a group with %d or more runs; smaller groups are shown as their points only.", MIN_VIOLIN) else
          sprintf("With %d or fewer runs per group no density is drawn -- the points are the data.", MIN_VIOLIN - 1),
        if (cg$simple) sprintf(paste0("; the arrow is the model's log2 fold change, drawn from the %s mean (red = higher,",
                                      " blue = lower), so its tip is where the model puts %s, which need not be that",
                                      " group's plain mean%s"), cg$ref, cg$test,
                               {why <- c(if (identical(vocab$absent, "Inferred")) sprintf("values inferred by %s", vocab$source),
                                         if (paired) sprintf("the %s effect it removes", bcol))
                                if (length(why)) sprintf(" (here: %s)", paste(why, collapse = " and ")) else ""}) else "",
        if (paired) sprintf(paste(" This is a within-%s contrast: grey lines join each %s's runs, the paired changes the",
                                  "model tested."), bcol, bcol) else ""), cap_status), CUT_NOTE))
    }, error = function(e) add_failed(sprintf("violin_top_%s.png", make.names(ct)), "violin", conditionMessage(e)))
  }
} else {
  why <- sprintf("no Expression_Matrix.csv in %s", de_dir)
  add_failed("pca.png", "pca", why); add_failed("heatmap_top.png", "heatmap", why)
  for (f in de_files) add_failed(sprintf("violin_top_%s.png", make.names(contrast_of(f))), "violin", why)
}

# ---- write figures.json -----------------------------------------------------
to_json <- function(figs, failed) {
  if (requireNamespace("jsonlite", quietly = TRUE))
    return(jsonlite::toJSON(list(figures = figs, failed = failed), auto_unbox = TRUE, pretty = TRUE))
  esc <- function(x) gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", x))
  arr <- function(l, keys) if (!length(l)) "[]" else paste0("[\n", paste(vapply(l, function(f)
    paste0("    {", paste(sprintf('"%s": "%s"', keys, vapply(keys, function(k) esc(f[[k]]), "")), collapse = ", "), "}"),
    ""), collapse = ",\n"), "\n  ]")
  paste0('{\n  "figures": ', arr(figs, c("file", "type", "caption")),
         ',\n  "failed": ', arr(failed, c("file", "type", "reason")), "\n}")
}
writeLines(to_json(figs, failed), file.path(outdir, "figures.json"))
cat(sprintf("[figures] wrote %d figure(s) + figures.json to %s%s\n", length(figs), normalizePath(outdir),
            if (length(failed)) sprintf("; %d not drawn (listed under \"failed\")", length(failed)) else ""))
