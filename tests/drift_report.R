#!/usr/bin/env Rscript
# Claude Analysis Drift Report
#
# Reads all golden baseline .rds files and produces a summary table
# showing how Claude's analysis has changed over time, with trend
# assessment and gene stability analysis.
#
# Usage:
#   Rscript tests/drift_report.R
#   Rscript tests/drift_report.R --csv    # also save to tests/golden/drift_report.csv

args <- commandArgs(trailingOnly = TRUE)
save_csv <- "--csv" %in% args

golden_dir <- file.path(dirname(dirname(
  if (interactive()) rstudioapi::getActiveDocumentContext()$path else "."
)), "tests", "golden")

# If run from project root
if (!dir.exists(golden_dir)) golden_dir <- "tests/golden"
if (!dir.exists(golden_dir)) stop("Cannot find tests/golden/ directory. Run from project root.")

files <- sort(list.files(golden_dir, "^claude_findings_.*\\.rds$", full.names = TRUE))
if (length(files) == 0) stop("No golden baseline files found in ", golden_dir)

# Known ground truth
KEY_PROTEINS <- c("NAMPT", "ALDH3A1", "LGALS1", "C3", "ITGA2",
                   "PRSS23", "TIMP1", "GSTM3", "CGA")
STABLE_MARKERS <- c("MED14", "TUBA1B", "TBC1D17")

# Read all findings
all_findings <- lapply(files, readRDS)

# Build summary table
rows <- lapply(seq_along(all_findings), function(i) {
  f <- all_findings[[i]]
  fname <- basename(files[i])

  # Extract date from filename: claude_findings_YYYYMMDD_HHMMSS.rds
  date_str <- sub("^claude_findings_", "", sub("\\.rds$", "", fname))
  date_parsed <- tryCatch(
    as.POSIXct(date_str, format = "%Y%m%d_%H%M%S"),
    error = function(e) NA
  )
  date_display <- if (!is.na(date_parsed)) format(date_parsed, "%Y-%m-%d %H:%M") else date_str

  # Gene overlap with previous run
  if (i > 1) {
    prev <- all_findings[[i - 1]]
    overlap <- intersect(f$genes_mentioned, prev$genes_mentioned)
    overlap_pct <- round(length(overlap) / max(length(prev$genes_mentioned), 1) * 100)
    overlap_str <- paste0(overlap_pct, "%")
  } else {
    overlap_pct <- NA
    overlap_str <- "(baseline)"
  }

  # Gene overlap with baseline (first run)
  if (i > 1) {
    baseline_overlap <- intersect(f$genes_mentioned, all_findings[[1]]$genes_mentioned)
    baseline_pct <- round(length(baseline_overlap) / max(length(all_findings[[1]]$genes_mentioned), 1) * 100)
    baseline_str <- paste0(baseline_pct, "%")
  } else {
    baseline_pct <- NA
    baseline_str <- "(baseline)"
  }

  # Key proteins found
  n_key <- if (!is.null(f$n_key_found)) f$n_key_found else sum(f$key_proteins_found)
  missing_key <- KEY_PROTEINS[!f$key_proteins_found]

  # Stable markers
  n_stable <- if (!is.null(f$n_stable_found)) f$n_stable_found else sum(f$stable_markers_found)

  # Sections
  n_sections <- sum(f$sections_found)
  total_sections <- length(f$sections_found)

  # Model info (added in later versions — may be NULL in older baselines)
  model_id <- f$model %||% "unknown"
  # Shorten model ID for display: "claude-sonnet-4-20250514" -> "sonnet-4"
  model_short <- sub("^claude-", "", sub("-\\d{8}$", "", model_id))

  data.frame(
    Date = date_display,
    Model = model_short,
    Key_Proteins = sprintf("%d/%d", n_key, length(KEY_PROTEINS)),
    Missing = if (length(missing_key) > 0) paste(missing_key, collapse = ", ") else "",
    Stable = sprintf("%d/%d", n_stable, length(STABLE_MARKERS)),
    Genes = f$n_genes,
    Sections = sprintf("%d/%d", n_sections, total_sections),
    Words = f$word_count,
    Up = f$n_up_mentions,
    Down = f$n_down_mentions,
    Tokens_In = f$input_tokens %||% NA_integer_,
    Tokens_Out = f$output_tokens %||% NA_integer_,
    Resp_Sec = f$response_time_sec %||% NA_real_,
    Hedges = f$n_hedges %||% NA_integer_,
    Confident = f$n_confident %||% NA_integer_,
    Overlap_Prev = overlap_str,
    Overlap_Base = baseline_str,
    stringsAsFactors = FALSE
  )
})

report <- do.call(rbind, rows)

# ============================================================================
#   PRINT REPORT
# ============================================================================

cat("\n")
cat("=== Claude Analysis Drift Report ===\n")
cat(sprintf("Golden baselines: %d (from %s to %s)\n",
  nrow(report), report$Date[1], report$Date[nrow(report)]))
cat("\n")

# Print table
print(report, row.names = FALSE, right = FALSE)

# ============================================================================
#   ALERTS
# ============================================================================

cat("\n--- Alerts ---\n")
alerts <- 0

# Check for model changes
if (nrow(report) >= 2) {
  for (i in 2:nrow(report)) {
    if (report$Model[i] != report$Model[i - 1]) {
      cat(sprintf("  [*] %s: MODEL CHANGED: %s -> %s\n",
        report$Date[i], report$Model[i - 1], report$Model[i]))
      alerts <- alerts + 1
    }
  }
}

# Check for key protein drops
for (i in seq_len(nrow(report))) {
  if (nzchar(report$Missing[i])) {
    cat(sprintf("  [!] %s: Missing key proteins: %s\n", report$Date[i], report$Missing[i]))
    alerts <- alerts + 1
  }
}

# Check for overlap drops below 50%
if (nrow(report) >= 2) {
  for (i in 2:nrow(report)) {
    pct <- as.numeric(sub("%", "", report$Overlap_Prev[i]))
    if (!is.na(pct) && pct < 50) {
      cat(sprintf("  [!] %s: Gene overlap dropped to %s vs previous (threshold: 50%%)\n",
        report$Date[i], report$Overlap_Prev[i]))
      alerts <- alerts + 1
    }
  }

  # Check for baseline drift below 40%
  for (i in 2:nrow(report)) {
    pct <- as.numeric(sub("%", "", report$Overlap_Base[i]))
    if (!is.na(pct) && pct < 40) {
      cat(sprintf("  [!] %s: Gene overlap vs BASELINE dropped to %s (threshold: 40%%)\n",
        report$Date[i], report$Overlap_Base[i]))
      alerts <- alerts + 1
    }
  }

  # Check for word count anomalies (< 500 or > 2x baseline)
  baseline_words <- rows[[1]]$Words
  for (i in 2:nrow(report)) {
    w <- report$Words[i]
    if (w < 500) {
      cat(sprintf("  [!] %s: Response only %d words (minimum: 500)\n", report$Date[i], w))
      alerts <- alerts + 1
    } else if (w > baseline_words * 2.5) {
      cat(sprintf("  [!] %s: Response %d words (%.1fx baseline) — possible prompt regression\n",
        report$Date[i], w, w / baseline_words))
      alerts <- alerts + 1
    }
  }

  # Check for sudden gene count changes (>50% swing)
  for (i in 2:nrow(report)) {
    prev_genes <- report$Genes[i - 1]
    curr_genes <- report$Genes[i]
    if (prev_genes > 0 && abs(curr_genes - prev_genes) / prev_genes > 0.5) {
      cat(sprintf("  [!] %s: Gene count changed from %d to %d (>50%% swing)\n",
        report$Date[i], prev_genes, curr_genes))
      alerts <- alerts + 1
    }
  }
}

if (alerts == 0) cat("  None — all runs within normal range.\n")

# ============================================================================
#   GENE STABILITY ANALYSIS
# ============================================================================

if (length(all_findings) >= 2) {
  cat("\n--- Model History ---\n")
  model_runs <- table(report$Model)
  for (m in names(model_runs)) {
    cat(sprintf("  %s: %d runs\n", m, model_runs[m]))
  }
  n_changes <- sum(report$Model[-1] != report$Model[-nrow(report)])
  cat(sprintf("  Model changes detected: %d\n", n_changes))

  cat("\n--- Gene Stability ---\n")

  # Core genes: mentioned in ALL runs
  all_gene_sets <- lapply(all_findings, `[[`, "genes_mentioned")
  core_genes <- Reduce(intersect, all_gene_sets)
  cat(sprintf("  Core genes (in every run): %d", length(core_genes)))
  if (length(core_genes) > 0 && length(core_genes) <= 20) {
    cat(sprintf(" — %s", paste(sort(core_genes), collapse = ", ")))
  }
  cat("\n")

  # Frequent genes: mentioned in >= 75% of runs
  gene_freq <- table(unlist(all_gene_sets))
  threshold_75 <- ceiling(length(all_findings) * 0.75)
  frequent_genes <- sort(names(gene_freq[gene_freq >= threshold_75]))
  if (length(frequent_genes) > length(core_genes)) {
    near_core <- setdiff(frequent_genes, core_genes)
    cat(sprintf("  Frequent genes (>=75%% of runs): %d", length(frequent_genes)))
    if (length(near_core) > 0 && length(near_core) <= 15) {
      cat(sprintf(" — additional: %s", paste(near_core, collapse = ", ")))
    }
    cat("\n")
  }

  # Volatile genes: mentioned in exactly 1 run
  one_off_genes <- names(gene_freq[gene_freq == 1])
  cat(sprintf("  One-off genes (single run only): %d\n", length(one_off_genes)))

  # All unique genes ever mentioned
  all_unique <- unique(unlist(all_gene_sets))
  cat(sprintf("  Total unique genes across all runs: %d\n", length(all_unique)))

  # Key protein consistency
  key_found_per_run <- sapply(all_findings, function(f) {
    if (!is.null(f$n_key_found)) f$n_key_found else sum(f$key_proteins_found)
  })
  if (all(key_found_per_run == length(KEY_PROTEINS))) {
    cat(sprintf("  Key protein recall: PERFECT (%d/%d in every run)\n",
      length(KEY_PROTEINS), length(KEY_PROTEINS)))
  } else {
    cat(sprintf("  Key protein recall: %s (across %d runs)\n",
      paste(sprintf("%d/%d", key_found_per_run, length(KEY_PROTEINS)), collapse = ", "),
      length(all_findings)))
  }

  # Genes gained/lost between consecutive runs
  if (length(all_findings) >= 2) {
    cat("\n  Run-to-run gene changes:\n")
    for (i in 2:length(all_findings)) {
      prev_genes <- all_findings[[i - 1]]$genes_mentioned
      curr_genes <- all_findings[[i]]$genes_mentioned
      gained <- setdiff(curr_genes, prev_genes)
      lost <- setdiff(prev_genes, curr_genes)
      date_i <- report$Date[i]
      if (length(gained) > 0)
        cat(sprintf("    %s gained: %s\n", date_i, paste(gained, collapse = ", ")))
      if (length(lost) > 0)
        cat(sprintf("    %s lost:   %s\n", date_i, paste(lost, collapse = ", ")))
    }
  }

  # Per-key-protein detail (fold changes, direction, example quotes)
  has_details <- any(sapply(all_findings, function(f) !is.null(f$key_protein_details)))
  if (has_details) {
    cat("\n--- Key Protein Details ---\n")
    for (gene in KEY_PROTEINS) {
      cat(sprintf("\n  %s:\n", gene))
      for (i in seq_along(all_findings)) {
        f <- all_findings[[i]]
        if (is.null(f$key_protein_details) || is.null(f$key_protein_details[[gene]])) {
          cat(sprintf("    %s: (no detail recorded)\n", report$Date[i]))
          next
        }
        d <- f$key_protein_details[[gene]]
        fc_str <- if (length(d$fold_changes) > 0) paste(d$fold_changes, collapse = ", ") else "none cited"
        cat(sprintf("    %s: %s | FC: %s | mentions: %d\n",
          report$Date[i], d$direction, fc_str, d$n_mentions))
        if (!is.na(d$example_sentence) && nchar(d$example_sentence) < 200)
          cat(sprintf("      \"%s\"\n", d$example_sentence))
      }
    }
  }

  # Hedging/confident example quotes
  has_examples <- any(sapply(all_findings, function(f) length(f$hedge_examples %||% character()) > 0))
  if (has_examples) {
    cat("\n--- Language Examples (latest run) ---\n")
    latest <- all_findings[[length(all_findings)]]
    if (length(latest$hedge_examples %||% character()) > 0) {
      cat("  Hedging:\n")
      for (ex in head(latest$hedge_examples, 3))
        cat(sprintf("    \"%s\"\n", substr(ex, 1, 150)))
    }
    if (length(latest$confident_examples %||% character()) > 0) {
      cat("  Confident:\n")
      for (ex in head(latest$confident_examples, 3))
        cat(sprintf("    \"%s\"\n", substr(ex, 1, 150)))
    }
  }
}

# ============================================================================
#   TREND ASSESSMENT (requires 4+ data points)
# ============================================================================

if (nrow(report) >= 3) {
  cat("\n--- Trends ---\n")

  # Overlap trend (vs previous run)
  overlaps_prev <- as.numeric(sub("%", "", report$Overlap_Prev[-1]))
  overlaps_prev <- overlaps_prev[!is.na(overlaps_prev)]

  if (length(overlaps_prev) >= 2) {
    cat(sprintf("  Gene overlap (vs prev): mean %.0f%%, range %.0f%%–%.0f%%",
      mean(overlaps_prev), min(overlaps_prev), max(overlaps_prev)))

    if (length(overlaps_prev) >= 3) {
      # Simple linear trend: positive = improving, negative = degrading
      trend_coef <- coef(lm(overlaps_prev ~ seq_along(overlaps_prev)))[2]
      if (abs(trend_coef) < 2) {
        cat(" → STABLE")
      } else if (trend_coef > 0) {
        cat(sprintf(" → IMPROVING (+%.1f%%/run)", trend_coef))
      } else {
        cat(sprintf(" → DEGRADING (%.1f%%/run)", trend_coef))
      }
    }
    cat("\n")
  }

  # Overlap trend (vs baseline)
  overlaps_base <- as.numeric(sub("%", "", report$Overlap_Base[-1]))
  overlaps_base <- overlaps_base[!is.na(overlaps_base)]
  if (length(overlaps_base) >= 2) {
    cat(sprintf("  Gene overlap (vs baseline): mean %.0f%%, range %.0f%%–%.0f%%",
      mean(overlaps_base), min(overlaps_base), max(overlaps_base)))
    if (length(overlaps_base) >= 3) {
      trend_coef <- coef(lm(overlaps_base ~ seq_along(overlaps_base)))[2]
      if (abs(trend_coef) < 2) {
        cat(" → STABLE")
      } else if (trend_coef > 0) {
        cat(sprintf(" → CONVERGING (+%.1f%%/run)", trend_coef))
      } else {
        cat(sprintf(" → DRIFTING (%.1f%%/run)", trend_coef))
      }
    }
    cat("\n")
  }

  # Word count trend
  cat(sprintf("  Word count: mean %d, range %d–%d",
    round(mean(report$Words)), min(report$Words), max(report$Words)))
  if (nrow(report) >= 4) {
    wc_coef <- coef(lm(report$Words ~ seq_len(nrow(report))))[2]
    if (abs(wc_coef) < 20) {
      cat(" → STABLE")
    } else if (wc_coef > 0) {
      cat(sprintf(" → GROWING (+%.0f words/run)", wc_coef))
    } else {
      cat(sprintf(" → SHRINKING (%.0f words/run)", wc_coef))
    }
  }
  cat("\n")

  # Gene count trend
  cat(sprintf("  Genes mentioned: mean %d, range %d–%d",
    round(mean(report$Genes)), min(report$Genes), max(report$Genes)))
  if (nrow(report) >= 4) {
    gc_coef <- coef(lm(report$Genes ~ seq_len(nrow(report))))[2]
    if (abs(gc_coef) < 1) {
      cat(" → STABLE")
    } else if (gc_coef > 0) {
      cat(sprintf(" → EXPANDING (+%.1f genes/run)", gc_coef))
    } else {
      cat(sprintf(" → CONTRACTING (%.1f genes/run)", gc_coef))
    }
  }
  cat("\n")
}

# ============================================================================
#   OVERALL VERDICT (requires 4+ data points)
# ============================================================================

if (nrow(report) >= 4) {
  cat("\n--- Overall Verdict ---\n")

  overlaps <- as.numeric(sub("%", "", report$Overlap_Prev[-1]))
  overlaps <- overlaps[!is.na(overlaps)]
  key_recall <- sapply(all_findings, function(f) {
    if (!is.null(f$n_key_found)) f$n_key_found else sum(f$key_proteins_found)
  })

  issues <- character(0)
  if (any(overlaps < 50)) issues <- c(issues, "overlap below 50%")
  if (any(key_recall < length(KEY_PROTEINS))) issues <- c(issues, "missing key proteins")
  if (sd(report$Words) / mean(report$Words) > 0.3) issues <- c(issues, "high word count variance")

  mean_overlap <- mean(overlaps)
  if (length(issues) == 0 && mean_overlap >= 65) {
    cat("  ✓ HEALTHY — Claude's analysis is consistent and complete.\n")
    cat(sprintf("    Mean overlap: %.0f%%, Key recall: %d/%d in all runs\n",
      mean_overlap, length(KEY_PROTEINS), length(KEY_PROTEINS)))
  } else if (length(issues) == 0) {
    cat("  ~ ACCEPTABLE — No critical issues, but overlap could be higher.\n")
    cat(sprintf("    Mean overlap: %.0f%% (target: >65%%)\n", mean_overlap))
  } else {
    cat(sprintf("  ⚠ ATTENTION NEEDED — %s\n", paste(issues, collapse = "; ")))
  }
}

cat("\n")

# ============================================================================
#   SAVE CSV
# ============================================================================

if (save_csv) {
  csv_path <- file.path(golden_dir, "drift_report.csv")
  write.csv(report, csv_path, row.names = FALSE)
  cat(sprintf("Saved to: %s\n", csv_path))
}
