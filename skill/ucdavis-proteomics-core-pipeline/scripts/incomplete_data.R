# incomplete_data.R -- "called on incomplete data": the per-call warning sign on EVERY DE route.
#
# A significant protein is the less trustworthy the fewer of the compared samples actually
# measured it. On the PXD063416 known-answer benchmark (human + yeast mix, 4 DIA windows,
# MaxLFQ + limma) calls made from proteins measured in every run were wrong 0.25-0.7% of the
# time, calls with missing values in a compared group 4-10%, and calls from a protein seen in
# only ONE run of a group 9-20%: flagged calls were about 10x more likely to be wrong. The signal
# needs no ground truth (filter-tuning prototype, Phinney lab, 2026-11-13,
# FILTER_TUNING_PROTOTYPE.md: "Decision for Brett" point 3).
#
# ONE DEFINITION (rule 3): what counts as MEASURED is not decided here. run_de.R's
# detection_counts() gives k (runs of a group in which the protein was measured) and n (runs of
# the group) from det_n -- the Detection_Matrix.csv definition (dpc: >= 1 precursor observed,
# 0 = inferred by the model; maxlfq: a finite value). This file only turns k and n into the flag.

INCOMPLETE_BENCHMARK <- list(
  data = "PXD063416 (human + yeast mix, 3-fold yeast change; four DIA windows)",
  wrong_pct_complete = "0.25-0.7",
  wrong_pct_incomplete = "4-10",
  wrong_pct_single_sample = "9-20",
  fold = "about 10x",
  citation = paste0("Phinney lab filter-tuning prototype, 2026-11-13 (FILTER_TUNING_PROTOTYPE.md), ",
                    "known-answer benchmark on PXD063416"))

INCOMPLETE_TIERS <- c(complete = "complete", incomplete = "incomplete", single = "single_sample")

# k: proteins x compared-groups matrix of runs with a measured value (NA = status unknown);
# n: runs per compared group (same order); grp: group names; sig: logical per protein, TRUE where
# the DE's own call rule calls it significant (NA/FALSE = not a call).
# -> data.frame(values_per_group, called_on_incomplete_data, incomplete_tier):
#   values_per_group          "Ctrl 3/3; Treat 1/3" (every compared group, the contrast's order)
#   called_on_incomplete_data TRUE: significant and some compared group has fewer than all its runs
#                             measured; FALSE: significant on complete data; NA: not significant
#   incomplete_tier           significant calls only: "complete" / "incomplete" / "single_sample"
#                             (a group with more than one run, measured in exactly one of them)
incomplete_flag_table <- function(k, n, grp, sig) {
  k <- as.matrix(k)
  stopifnot(ncol(k) == length(n), length(n) == length(grp), nrow(k) == length(sig))
  nn <- matrix(as.integer(n), nrow(k), length(n), byrow = TRUE)
  vpg <- vapply(seq_len(nrow(k)), function(i)
    if (anyNA(k[i, ])) NA_character_ else paste(sprintf("%s %d/%d", grp, as.integer(k[i, ]), nn[i, ]), collapse = "; "),
    "")
  unknown <- rowSums(is.na(k)) > 0
  incomplete <- rowSums(k < nn, na.rm = TRUE) > 0
  single <- rowSums(k == 1L & nn > 1L, na.rm = TRUE) > 0
  is_call <- !is.na(sig) & sig & !unknown
  tier <- ifelse(!is_call, NA_character_,
           ifelse(single, INCOMPLETE_TIERS[["single"]],
            ifelse(incomplete, INCOMPLETE_TIERS[["incomplete"]], INCOMPLETE_TIERS[["complete"]])))
  data.frame(values_per_group = vpg,
             called_on_incomplete_data = ifelse(is_call, incomplete, NA),
             incomplete_tier = unname(tier), stringsAsFactors = FALSE)
}

# The sentence every report reuses (de_provenance.json `incomplete_data.rule`).
incomplete_rule_text <- function(zero_means) {
  sprintf(paste0("for each contrast, values_per_group counts the runs of every group it compares in ",
                 "which the protein was measured (Detection_Matrix.csv > 0; 0 = %s). A significant ",
                 "protein (the DE's own call rule: adj.P.Val < the adjusted p-value threshold) is ",
                 "called_on_incomplete_data = TRUE when any compared group has fewer than all its ",
                 "runs measured, FALSE when it was measured in every run, NA when not significant; ",
                 "incomplete_tier splits out 'single_sample' (a group of several runs in which it ",
                 "was measured in one only)."), zero_means)
}

# The benchmark statement, from INCOMPLETE_BENCHMARK (never retyped by a consumer).
incomplete_benchmark_text <- function(b = INCOMPLETE_BENCHMARK) {
  sprintf(paste0("On the %s known-answer benchmark, calls from proteins measured in every run were ",
                 "wrong %s%% of the time, calls with missing values in a compared group %s%%, and ",
                 "calls from a protein seen in one run of a group %s%%: flagged calls were %s more ",
                 "likely to be wrong (%s)."),
          sub(" \\(.*$", "", b$data), b$wrong_pct_complete, b$wrong_pct_incomplete,
          b$wrong_pct_single_sample, b$fold, b$citation)
}
