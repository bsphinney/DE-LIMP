# ==============================================================================
#  HELPER FUNCTIONS — General utilities
# ==============================================================================

# --- QC Stats Calculation ---
# Memory-optimized: reads only needed columns via Arrow col_select,
# then aggregates before collecting into R memory.
get_diann_stats_r <- function(file_path) {
  tryCatch({
    # Check which columns are available without reading data
    available_cols <- names(arrow::read_parquet(file_path, as_data_frame = FALSE))

    needed_cols <- c("Run", "Protein.Group", "Q.Value")
    has_pg_q <- "PG.Q.Value" %in% available_cols
    has_ms1  <- "Ms1.Apex.Area" %in% available_cols
    if (has_pg_q) needed_cols <- c(needed_cols, "PG.Q.Value")
    if (has_ms1)  needed_cols <- c(needed_cols, "Ms1.Apex.Area")

    # Read only the needed columns (saves 70-80% memory for large files)
    df <- arrow::read_parquet(file_path, col_select = dplyr::all_of(needed_cols))
    if ("Q.Value" %in% names(df)) df <- df %>% dplyr::filter(Q.Value <= 0.01)

    stats_df <- df %>%
      dplyr::group_by(Run) %>%
      dplyr::summarise(
        Precursors = dplyr::n(),
        Proteins = if(has_pg_q) {
          dplyr::n_distinct(Protein.Group[PG.Q.Value <= 0.01])
        } else {
          dplyr::n_distinct(Protein.Group)
        },
        MS1_Signal = if(has_ms1) sum(Ms1.Apex.Area, na.rm = TRUE) else NA_real_,
        .groups = 'drop'
      ) %>% dplyr::arrange(Run)

    # Free the large intermediate df immediately
    rm(df); gc(verbose = FALSE)

    return(stats_df)
  }, error = function(e) { data.frame(Run = "Error", Precursors = 0, Proteins = 0, MS1_Signal = 0) })
}

# --- Z-Score Utility ---
cal_z_score <- function(x) { (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE) }

# --- DPC-Quant Detection Class ---
# Classifies each protein based on n.observations across all samples.
# Returns a character vector: "Detected_All", "Detected_Partial", or "Inferred_All".
# n_obs_mat: matrix with proteins as rows, samples as columns (from y_protein$other$n.observations)
# protein_ids: character vector of protein IDs to classify (must match rownames of n_obs_mat)
compute_detection_class <- function(n_obs_mat, protein_ids) {
  if (is.null(n_obs_mat)) return(rep(NA_character_, length(protein_ids)))
  rn <- rownames(n_obs_mat)
  vapply(protein_ids, function(pid) {
    idx <- match(pid, rn)
    if (is.na(idx)) return(NA_character_)
    obs <- n_obs_mat[idx, ]
    if (all(obs > 0)) {
      "Detected_All"
    } else if (all(obs == 0)) {
      "Inferred_All"
    } else {
      "Detected_Partial"
    }
  }, character(1))
}

# --- Auto-detect Organism ---
detect_organism_db <- function(protein_ids) {
  ORGANISM_DB_MAP <- list(
    "_HUMAN" = "org.Hs.eg.db", "_MOUSE" = "org.Mm.eg.db", "_RAT"   = "org.Rn.eg.db",
    "_BOVIN" = "org.Bt.eg.db", "_CANLF" = "org.Cf.eg.db", "_CHICK" = "org.Gg.eg.db",
    "_DROME" = "org.Dm.eg.db", "_CAEEL" = "org.Ce.eg.db", "_DANRE" = "org.Dr.eg.db",
    "_YEAST" = "org.Sc.sgd.db", "_ARATH" = "org.At.tair.db", "_PIG"   = "org.Ss.eg.db"
  )
  for (suffix in names(ORGANISM_DB_MAP)) {
    if (any(grepl(suffix, protein_ids, ignore.case = TRUE))) {
      return(ORGANISM_DB_MAP[[suffix]])
    }
  }
  return("org.Hs.eg.db")
}
