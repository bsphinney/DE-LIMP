#!/usr/bin/env Rscript
# =============================================================================
# repro_script.R  --  Emit the analysis as a flat, runnable R script.
#
# This is the skill's answer to "just show me the R code". It is the same idea
# as DE-LIMP's reproducibility log (app.R `add_to_log()` -> R/server_session.R
# `download_repro_log`): plain top-to-bottom R with every value written out
# literally -- the actual report path, the actual FDR cutoff, the actual sample
# -> group assignment, the actual contrasts -- so a reviewer can read it, and a
# collaborator can run it, without installing the skill, conda, or any of the
# orchestration around it.
#
# It is generated FROM the objects that really ran (run_de.R passes them in),
# never re-derived from a config file, so it cannot drift from the result.
# DE-LIMP architectural rule #1: the pipeline self-describes.
#
# write_repro_script() returns the character vector of lines it wrote, so
# run_de.R can also drop them into the DE-LIMP session RDS as `repro_log` --
# the GUI then shows this same code in its Reproducibility tab.
# =============================================================================

# --- literal emitters: never let a value reach the script unescaped ----------
.rq <- function(x) {                          # quote a string for R source
  x <- gsub("\\\\", "\\\\\\\\", as.character(x))
  x <- gsub("'", "\\\\'", x)
  paste0("'", x, "'")
}
.rnum <- function(x) format(x, scientific = FALSE, trim = TRUE)
.rvec <- function(x) paste0("c(", paste(vapply(x, .rq, ""), collapse = ", "), ")")

# A named map, one entry per line, in DE-LIMP's reproducibility-log style.
.rmap <- function(varname, names_vec, values_vec) {
  c(sprintf("%s <- c(", varname),
    paste0(sprintf("  %s = %s", vapply(names_vec, .rq, ""), vapply(values_vec, .rq, "")),
           c(rep(",", max(0, length(names_vec) - 1)), "")),
    ")")
}

write_repro_script <- function(path,
                               method,           # "dpc" | "maxlfq"
                               input,            # report path as it was read
                               format,           # "parquet" | "tsv"
                               q_cutoff,
                               q_columns,        # q-value columns actually filtered on
                               eq_cutoff, pgq_cutoff,
                               cov_min_frac,     # maxlfq coverage filter
                               meta,             # design, already aligned to the matrix
                               covariates,       # covariate column names in `meta`
                               formula_parts,
                               forms,            # contrast strings
                               adjp_thr, logfc_ref,
                               ann_cols = character(0),   # Genes / Protein.Names, if present
                               descriptor = NULL,
                               timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                               # TRUE only when this is being generated for a run that
                               # predates the feature, from that run's recorded provenance
                               # rather than from live objects. The header must say so:
                               # the normal file's claim ("from the objects that actually
                               # ran") is precisely what makes it trustworthy, so a
                               # reconstruction must not borrow it. Architectural rule #2 --
                               # a value that reaches a reader is tagged or not asserted.
                               reconstructed = FALSE,
                               reconstructed_note = NULL) {

  is_dpc  <- identical(method, "dpc")
  rpt_abs <- normalizePath(input, mustWork = FALSE)
  eq_on   <- !is.na(eq_cutoff)  && eq_cutoff  > 0
  pgq_on  <- !is.na(pgq_cutoff) && pgq_cutoff > 0

  L <- c(
    "# =============================================================================",
    "# The analysis, as R code.",
    "#",
    sprintf("# Generated %s by the ucdavis-proteomics-core-pipeline skill.", timestamp),
    if (reconstructed) c(
    "#",
    "# RECONSTRUCTED AFTER THE FACT. This run predates the skill version that emits",
    "# this file, so it was not written by the run itself. Every value below was read",
    "# back from what the run recorded about itself (de_provenance.json, conditions.csv,",
    "# the expression matrix's column order) -- none of it was inferred. It is a faithful",
    "# description of the analysis, but it does not carry the generated-from-live-objects",
    "# guarantee that files from newer runs do, and package versions may have moved since.",
    if (!is.null(reconstructed_note)) paste0("# ", reconstructed_note) else NULL)
    else c(
    "# Written from the objects that actually produced the results -- not re-derived",
    "# from a config file, so it cannot describe a different analysis than the CSVs",
    "# beside it."),
    "#",
    "# This is the whole differential-expression analysis. It needs only R and the",
    "# packages loaded below; no conda environment, no skill install, no shell",
    "# scripts. Edit the two paths at the top and run:",
    "#",
    "#     Rscript reproducibility_log.R",
    "#",
    if (!is.null(descriptor)) sprintf("# Pipeline : %s", descriptor$display_label) else NULL,
    if (!is.null(descriptor)) sprintf("# Citation : %s", descriptor$citation) else NULL,
    "# =============================================================================",
    "",
    "# --- Paths (edit these) ------------------------------------------------------",
    sprintf("report <- %s", .rq(rpt_abs)),
    sprintf("outdir <- %s   # written here; the original results are left untouched",
            .rq("de_results_rerun")),
    "dir.create(outdir, showWarnings = FALSE, recursive = TRUE)",
    ""
  )

  # ---- libraries -------------------------------------------------------------
  libs <- if (is_dpc) {
    if (eq_on || pgq_on) "library(limpa); library(limma); library(arrow); library(dplyr)"
    else                 "library(limpa); library(limma)"
  } else "library(arrow); library(dplyr); library(tidyr); library(limma)"
  L <- c(L, "# --- Load Required Libraries -------------------------------------------------", libs, "")

  # ---- experimental design ----------------------------------------------------
  .hdr <- function(txt) {                      # "# --- txt ---...---" padded to 80 cols
    base <- paste0("# --- ", txt, " ")
    paste0(base, strrep("-", max(3, 80 - nchar(base))))
  }
  L <- c(L,
    .hdr(sprintf("Experimental design (%d samples%s)", nrow(meta),
                 if (length(covariates)) sprintf(", covariates: %s", paste(covariates, collapse = ", ")) else "")),
    .rmap("group_map", meta$File.Name, meta$Group))

  df_cols <- "Group = unname(group_map)"
  for (cv in covariates) {
    mapname <- paste0(tolower(cv), "_map")
    L <- c(L, .rmap(mapname, meta$File.Name, meta[[cv]]))
    df_cols <- paste0(df_cols, sprintf(", %s = unname(%s)", cv, mapname))
  }
  L <- c(L,
    sprintf("metadata <- data.frame(File.Name = names(group_map), %s, stringsAsFactors = FALSE)", df_cols),
    "")

  # ---- quantification ---------------------------------------------------------
  if (is_dpc) {
    src <- "report"
    if (eq_on || pgq_on) {
      filt <- character(0)
      if (eq_on)  filt <- c(filt, sprintf("Empirical.Quality >= %s", .rnum(eq_cutoff)))
      if (pgq_on) filt <- c(filt, sprintf("PG.MaxLFQ.Quality >= %s", .rnum(pgq_cutoff)))
      L <- c(L,
        "# --- QuantUMS pre-filter (applied to the report before limpa reads it) -------",
        "flt <- arrow::open_dataset(report)",
        sprintf("flt <- dplyr::filter(flt, %s)", paste(filt, collapse = ", ")),
        "report_filtered <- file.path(tempdir(), 'quantums_filtered.parquet')",
        "arrow::write_parquet(dplyr::collect(flt), report_filtered)",
        "")
      src <- "report_filtered"
    }
    L <- c(L,
      "# --- 1. Read the DIA-NN report, applying the identification FDR cutoff -------",
      sprintf("dat <- limpa::readDIANN(%s, format = %s, q.cutoffs = %s,",
              src, .rq(format), .rnum(q_cutoff)),
      sprintf("                        q.columns = %s)", .rvec(q_columns)),
      "",
      "# --- 2. Keep only the runs that appear in the design -------------------------",
      "dat <- dat[, colnames(dat$E) %in% metadata$File.Name]",
      "",
      "# --- 3. Normalise and roll precursors up to proteins (DPC-CN + DPC-Quant) ----",
      "dpcfit    <- limpa::dpcCN(dat)",
      "y_protein <- limpa::dpcQuant(dat, 'Protein.Group', dpc = dpcfit)",
      "E <- y_protein$E",
      "")
  } else {
    sel <- unique(c("Run", "Protein.Group", "PG.MaxLFQ", q_columns, ann_cols,
                    if (eq_on) "Empirical.Quality", if (pgq_on) "PG.MaxLFQ.Quality"))
    qfilt <- paste(sprintf("%s <= %s", q_columns, .rnum(q_cutoff)), collapse = ", ")
    if (eq_on)  qfilt <- paste0(qfilt, sprintf(", Empirical.Quality >= %s", .rnum(eq_cutoff)))
    if (pgq_on) qfilt <- paste0(qfilt, sprintf(", PG.MaxLFQ.Quality >= %s", .rnum(pgq_cutoff)))
    L <- c(L,
      "# --- 1. Read the DIA-NN report, applying FDR (+ QuantUMS) filters ------------",
      sprintf("rows <- arrow::open_dataset(report, format = %s) |>", .rq(format)),
      sprintf("  dplyr::select(dplyr::all_of(%s)) |>", .rvec(sel)),
      sprintf("  dplyr::filter(%s) |>", qfilt),
      "  dplyr::filter(Run %in% metadata$File.Name) |>",
      "  dplyr::collect()",
      "",
      "# --- 2. One PG.MaxLFQ per (protein, run); pivot wide; log2; quantile-normalise",
      "pg_run <- rows |>",
      "  dplyr::group_by(Protein.Group, Run) |>",
      "  dplyr::summarise(PG.MaxLFQ = max(PG.MaxLFQ, na.rm = TRUE), .groups = 'drop') |>",
      "  dplyr::mutate(PG.MaxLFQ = ifelse(is.finite(PG.MaxLFQ), PG.MaxLFQ, NA_real_))",
      "wide <- tidyr::pivot_wider(pg_run, id_cols = Protein.Group,",
      "                           names_from = Run, values_from = PG.MaxLFQ)",
      "E <- as.matrix(wide[, -1, drop = FALSE]); rownames(E) <- wide$Protein.Group",
      "E[E <= 0 | !is.finite(E)] <- NA_real_",
      "E <- log2(E)",
      "E <- E[, colSums(is.finite(E)) >= 2, drop = FALSE]   # drop runs with <2 proteins",
      "E <- limma::normalizeBetweenArrays(E, method = 'quantile')",
      "",
      sprintf("# --- 3. Coverage filter: a protein must be quantified in >= %.0f%% of samples",
              100 * cov_min_frac),
      sprintf("min_obs <- max(2, ceiling(%s * ncol(E)))", .rnum(cov_min_frac)),
      "E <- E[rowSums(!is.na(E)) >= min_obs, , drop = FALSE]",
      "")
  }

  # ---- design matrix ----------------------------------------------------------
  L <- c(L,
    "# --- 4. Align the design to the matrix columns, build the model matrix -------",
    "metadata <- metadata[match(colnames(E), metadata$File.Name), ]",
    "groups <- factor(metadata$Group)")
  for (cv in covariates)
    L <- c(L, sprintf("%s <- factor(metadata$%s)", cv, cv))
  L <- c(L,
    sprintf("design <- model.matrix(~ 0 + %s)", paste(formula_parts, collapse = " + ")),
    "colnames(design) <- sub('^groups', '', colnames(design))",
    "")

  # ---- fit --------------------------------------------------------------------
  L <- c(L,
    "# --- 5. Fit the model and test the contrasts ---------------------------------",
    if (is_dpc) "fit <- limpa::dpcDE(y_protein, design, plot = FALSE)"
    else        "fit <- limma::lmFit(E, design)",
    sprintf("contrast_matrix <- limma::makeContrasts(contrasts = %s, levels = design)", .rvec(forms)),
    "fit <- limma::contrasts.fit(fit, contrast_matrix)",
    "fit <- limma::eBayes(fit)",
    "")

  # ---- results ----------------------------------------------------------------
  # Significance is the BH adjusted p-value alone -- the hypothesis eBayes actually
  # tested. logfc_ref is reported as a description of the significant set, never as
  # a filter (see run_de.R for why).
  # Protein annotation (Genes / Protein.Names). It rides along in limpa's fit for the
  # dpc path, but lmFit(E, design) carries none, so on maxlfq it has to be merged back
  # in or the DE tables lose their gene symbols.
  has_ann <- length(ann_cols) > 0
  if (has_ann) {
    L <- c(L, "# --- Protein annotation (gene symbols / names) --------------------------------")
    if (is_dpc) {
      L <- c(L,
        sprintf("ann <- y_protein$genes[, %s, drop = FALSE]", .rvec(ann_cols)),
        "ann$Protein.Group <- if ('Protein.Group' %in% names(y_protein$genes))",
        "  as.character(y_protein$genes$Protein.Group) else rownames(y_protein$genes)")
    } else {
      L <- c(L,
        "ann <- rows |>",
        "  dplyr::group_by(Protein.Group) |>",
        sprintf("  dplyr::summarise(dplyr::across(dplyr::all_of(%s),", .rvec(ann_cols)),
        "    ~ names(sort(table(.x), decreasing = TRUE))[1]), .groups = 'drop') |>",
        "  as.data.frame()")
    }
    L <- c(L, sprintf("ann <- ann[, %s, drop = FALSE]", .rvec(c("Protein.Group", ann_cols))), "")
  }

  L <- c(L,
    "# --- 6. Write the expression matrix and one results table per contrast -------",
    "# Significance is the BH-adjusted p-value alone -- the hypothesis eBayes tested.",
    sprintf("# |log2FC| >= %s is only counted below for description; it filters nothing.", .rnum(logfc_ref)),
    "expr <- data.frame(Protein.Group = rownames(E), check.names = FALSE)",
    if (has_ann) "expr <- merge(expr, ann, by = 'Protein.Group', all.x = TRUE, sort = FALSE)" else NULL,
    "expr <- merge(expr, data.frame(Protein.Group = rownames(E), E, check.names = FALSE),",
    "              by = 'Protein.Group', all.x = TRUE, sort = FALSE)",
    "write.csv(expr, file.path(outdir, 'Expression_Matrix.csv'), row.names = FALSE)",
    "",
    sprintf("for (cn in %s) {", .rvec(forms)),
    "  tt <- limma::topTable(fit, coef = cn, number = Inf, adjust.method = 'BH')",
    "  tt$Protein.Group <- rownames(tt)",
    if (has_ann) c(
    "  miss <- setdiff(names(ann), c('Protein.Group', names(tt)))   # don't duplicate columns",
    "  if (length(miss)) tt <- merge(tt, ann[, c('Protein.Group', miss), drop = FALSE],",
    "                                by = 'Protein.Group', all.x = TRUE, sort = FALSE)") else NULL,
    "  tt <- tt[order(tt$adj.P.Val), ]",
    sprintf("  write.csv(tt, file.path(outdir, sprintf('DE_%s_%%s.csv', make.names(cn))), row.names = FALSE)",
            method),
    sprintf("  sig <- subset(tt, !is.na(adj.P.Val) & adj.P.Val < %s)", .rnum(adjp_thr)),
    sprintf("  cat(sprintf('%%-20s %%d proteins, %%d significant (adj.P < %s), %%d of those |log2FC| >= %s\\n',",
            .rnum(adjp_thr), .rnum(logfc_ref)),
    sprintf("              cn, nrow(tt), nrow(sig), sum(abs(sig$logFC) >= %s, na.rm = TRUE)))", .rnum(logfc_ref)),
    "}",
    "",
    "# --- 7. Record the exact package versions this run used ----------------------",
    "writeLines(capture.output(sessionInfo()), file.path(outdir, 'sessionInfo.txt'))",
    "sessionInfo()"
  )

  L <- L[!vapply(L, is.null, logical(1))]
  writeLines(L, path)
  invisible(L)
}
