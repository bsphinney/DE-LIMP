# =============================================================================
# diann_q_columns.R -- R mirror of scripts/diann_q_columns.py.
#
# The canonical definition and the full rationale live in the Python file; read
# that one. This exists because run_de.R and build_maxlfq.R are standalone
# Rscripts and treat jsonlite as OPTIONAL (run_de.R ships a manual JSON writer
# fallback), so making the FDR filter depend on parsing a shared JSON file would
# trade a duplication bug for a hard dependency on the identification filter --
# a much worse failure.
#
# The two files are asserted identical by tests/test_q_columns.py, which parses
# this file and compares it to the Python module. Drift is a CI failure, not a
# silent divergence, which is the property that actually matters here.
#
# Two DISTINCT concepts -- do not merge them:
#   DIANN_FDR_REQUIRED + DIANN_FDR_OPTIONAL
#       the FILTER SET: every one is applied, ANDed, at the q-cutoff.
#   DIANN_PROTEIN_Q_PREFERENCE
#       a PREFERENCE ORDER for reading ONE protein-level q-column; first
#       available wins and the rest are ignored.
#
# Column semantics (DIA-NN 2.6 README, "Main output reference"):
#   Q.Value            run-specific precursor q-value
#   Global.Q.Value     global (experiment-wide) precursor q-value
#   Lib.Q.Value        library-entry q-value; 'global' when DIA-NN built the lib
#   PG.Q.Value         RUN-SPECIFIC protein-group q-value, despite sitting beside
#                      the Global.* pair
#   Global.PG.Q.Value  global protein-group q-value
#   Lib.PG.Q.Value     library-entry protein-group q-value
# =============================================================================

DIANN_FDR_REQUIRED <- c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value")
DIANN_FDR_OPTIONAL <- c("PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value")
DIANN_PROTEIN_Q_PREFERENCE <- c("Global.PG.Q.Value", "Global.Q.Value",
                                "Lib.PG.Q.Value", "PG.Q.Value", "Q.Value")

#' Columns to FILTER on, restricted to those the report actually has.
#'
#' Never returns a column the report lacks: limpa::readDIANN() errors on an
#' unknown name, and arrow silently returns ZERO ROWS when filtering a column
#' that select() dropped -- the defect that made MaxLFQ produce nothing in
#' DE-LIMP v4.0.2-4.0.3.
diann_fdr_columns <- function(available = NULL) {
  all_cols <- c(DIANN_FDR_REQUIRED, DIANN_FDR_OPTIONAL)
  if (is.null(available)) return(all_cols)
  all_cols[all_cols %in% available]
}

#' The single best protein-level q-column present, or NA_character_.
diann_protein_q_column <- function(available) {
  hit <- DIANN_PROTEIN_Q_PREFERENCE[DIANN_PROTEIN_Q_PREFERENCE %in% available]
  if (length(hit)) hit[1] else NA_character_
}
