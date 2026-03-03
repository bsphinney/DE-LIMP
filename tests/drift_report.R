#!/usr/bin/env Rscript
# Claude Analysis Drift Report
#
# Reads all golden baseline .rds files and produces a summary table
# showing how Claude's analysis has changed over time.
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

  # Key proteins found
  n_key <- if (!is.null(f$n_key_found)) f$n_key_found else sum(f$key_proteins_found)
  missing_key <- KEY_PROTEINS[!f$key_proteins_found]

  # Stable markers
  n_stable <- if (!is.null(f$n_stable_found)) f$n_stable_found else sum(f$stable_markers_found)

  # Sections
  n_sections <- sum(f$sections_found)
  total_sections <- length(f$sections_found)

  data.frame(
    Date = date_display,
    Key_Proteins = sprintf("%d/%d", n_key, length(KEY_PROTEINS)),
    Missing = if (length(missing_key) > 0) paste(missing_key, collapse = ", ") else "",
    Stable = sprintf("%d/%d", n_stable, length(STABLE_MARKERS)),
    Genes = f$n_genes,
    Sections = sprintf("%d/%d", n_sections, total_sections),
    Words = f$word_count,
    Up = f$n_up_mentions,
    Down = f$n_down_mentions,
    Overlap = overlap_str,
    stringsAsFactors = FALSE
  )
})

report <- do.call(rbind, rows)

# Print report
cat("\n")
cat("=== Claude Analysis Drift Report ===\n")
cat(sprintf("Golden baselines: %d (from %s to %s)\n",
  nrow(report), report$Date[1], report$Date[nrow(report)]))
cat("\n")

# Print table
print(report, row.names = FALSE, right = FALSE)

# Print alerts
cat("\n--- Alerts ---\n")
alerts <- 0

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
    pct <- as.numeric(sub("%", "", report$Overlap[i]))
    if (!is.na(pct) && pct < 50) {
      cat(sprintf("  [!] %s: Gene overlap dropped to %s (threshold: 50%%)\n",
        report$Date[i], report$Overlap[i]))
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
}

if (alerts == 0) cat("  None — all runs within normal range.\n")

# Trend summary
if (nrow(report) >= 3) {
  cat("\n--- Trends ---\n")
  overlaps <- as.numeric(sub("%", "", report$Overlap[-1]))
  overlaps <- overlaps[!is.na(overlaps)]
  if (length(overlaps) >= 2) {
    cat(sprintf("  Gene overlap: mean %.0f%%, range %.0f%%–%.0f%%\n",
      mean(overlaps), min(overlaps), max(overlaps)))
  }
  cat(sprintf("  Word count: mean %d, range %d–%d\n",
    round(mean(report$Words)), min(report$Words), max(report$Words)))
  cat(sprintf("  Genes mentioned: mean %d, range %d–%d\n",
    round(mean(report$Genes)), min(report$Genes), max(report$Genes)))
}

cat("\n")

# Save CSV if requested
if (save_csv) {
  csv_path <- file.path(golden_dir, "drift_report.csv")
  write.csv(report, csv_path, row.names = FALSE)
  cat(sprintf("Saved to: %s\n", csv_path))
}
