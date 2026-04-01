# ==============================================================================
#  HELPERS — DDA Search (Sage pipeline)
#  Pure utility functions — no Shiny reactivity.
#  Called from: server_dda.R
# ==============================================================================

#' Generate a Sage search config JSON
#'
#' @param fasta_path Path to reference FASTA
#' @param raw_paths Character vector of .d directory paths
#' @param output_dir Where Sage writes results
#' @param preset One of "standard", "phospho", "tmt"
#' @param missed_cleavages Integer (default 2)
#' @param precursor_tol_ppm Precursor mass tolerance in ppm (default 20)
#' @param fragment_tol_da Fragment tolerance in Da (default 0.05)
#' @param min_peaks Min fragment peaks per spectrum (default 6)
#' @return Path to written sage.json file
generate_sage_config <- function(
  fasta_path,
  raw_paths,
  output_dir,
  preset           = "standard",
  missed_cleavages = 2,
  precursor_tol_ppm = 20,
  fragment_tol_da   = 0.05,
  min_peaks         = 6
) {
  # Preset-based modification definitions
  static_mods <- list("C" = 57.0215)  # carbamidomethyl always fixed

  variable_mods <- switch(preset,
    "standard" = list(
      list(
        sites    = list("M"),
        mass     = 15.9949,
        position = jsonlite::unbox(NA)
      )
    ),
    "phospho" = list(
      list(sites = list("M"), mass = 15.9949, position = jsonlite::unbox(NA)),
      list(sites = list("S", "T", "Y"), mass = 79.9663, position = jsonlite::unbox(NA))
    ),
    "tmt" = list(
      list(sites = list("K"), mass = 229.1629, position = jsonlite::unbox(NA)),
      list(sites = list("^"), mass = 229.1629, position = jsonlite::unbox(NA)),
      list(sites = list("M"), mass = 15.9949, position = jsonlite::unbox(NA))
    ),
    # default
    list(
      list(sites = list("M"), mass = 15.9949, position = jsonlite::unbox(NA))
    )
  )

  # TMT static mods: add TMT to n-term as static
  if (preset == "tmt") {
    static_mods[["^"]] <- 229.1629
  }

  config <- list(
    database = list(
      bucket_size = jsonlite::unbox(32768L),
      enzyme = list(
        missed_cleavages = jsonlite::unbox(as.integer(missed_cleavages)),
        min_len          = jsonlite::unbox(7L),
        max_len          = jsonlite::unbox(50L),
        cleave_at        = jsonlite::unbox("KR"),
        restrict         = jsonlite::unbox("P")
      ),
      fragment_min_mz  = jsonlite::unbox(150.0),
      fragment_max_mz  = jsonlite::unbox(2000.0),
      peptide_min_mass = jsonlite::unbox(500.0),
      peptide_max_mass = jsonlite::unbox(5000.0),
      ion_kinds        = c("b", "y"),
      min_ion_index    = jsonlite::unbox(2L),
      static_mods      = static_mods,
      variable_mods    = variable_mods,
      max_variable_mods = jsonlite::unbox(2L),
      decoy_tag        = jsonlite::unbox("rev_"),
      generate_decoys  = jsonlite::unbox(TRUE),
      fasta            = jsonlite::unbox(fasta_path)
    ),
    precursor_tol = list(ppm = c(-precursor_tol_ppm, precursor_tol_ppm)),
    fragment_tol  = list(ppm = c(-fragment_tol_da * 1000, fragment_tol_da * 1000)),
    report_psms   = jsonlite::unbox(1L),
    min_peaks     = jsonlite::unbox(as.integer(min_peaks)),
    max_peaks     = jsonlite::unbox(150L),
    max_fragment_charge = jsonlite::unbox(2L),
    chimera       = jsonlite::unbox(FALSE),
    predict_rt    = jsonlite::unbox(TRUE),
    mzml_paths    = raw_paths
  )

  # Add LFQ config for non-TMT

  if (preset != "tmt") {
    config$quant <- list(
      lfq = list(
        peak_scoring   = jsonlite::unbox("Hybrid"),
        integration    = jsonlite::unbox("Sum"),
        spectral_angle = jsonlite::unbox(0.7)
      )
    )
  } else {
    config$quant <- list(
      tmt = list(level = jsonlite::unbox(2L))
    )
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  config_path <- file.path(output_dir, "sage.json")
  jsonlite::write_json(config, config_path, auto_unbox = FALSE, pretty = TRUE, null = "null")
  message("[DDA] Sage config written to: ", config_path)
  config_path
}


#' Parse Sage search results into DE-LIMP-compatible format
#'
#' @param results_path Path to results.sage.tsv
#' @param lfq_path Path to lfq.tsv
#' @param fdr_threshold PSM-level FDR cutoff (default 0.01)
#' @param protein_fdr_threshold Protein-level FDR cutoff (default 0.01)
#' @return List with $psms, $lfq_wide (log2 matrix), $protein_meta
parse_sage_results <- function(
  results_path,
  lfq_path,
  fdr_threshold         = 0.01,
  protein_fdr_threshold = 0.01
) {
  if (!file.exists(results_path)) stop("Sage results file not found: ", results_path)
  if (!file.exists(lfq_path)) stop("Sage LFQ file not found: ", lfq_path)

  # PSM table
  psms <- data.table::fread(results_path)
  message("[DDA] Read ", nrow(psms), " total PSMs from Sage")

  # FDR filter
  psms_filtered <- psms[spectrum_q <= fdr_threshold & protein_q <= protein_fdr_threshold]
  message("[DDA] After FDR filter (spectrum_q <= ", fdr_threshold,
          ", protein_q <= ", protein_fdr_threshold, "): ", nrow(psms_filtered), " PSMs")

  # Per-protein peptide/spectra counts
  pep_counts <- psms_filtered[, .(
    NPeptides = data.table::uniqueN(peptide),
    NSpectra  = .N
  ), by = proteins]

  # LFQ matrix: long -> wide
  lfq <- data.table::fread(lfq_path)
  lfq_wide <- data.table::dcast(lfq, proteins ~ filename, value.var = "intensity")
  rownames_col <- lfq_wide$proteins
  lfq_mat <- as.matrix(lfq_wide[, -"proteins", with = FALSE])
  rownames(lfq_mat) <- rownames_col

  # Log2 transform (Sage outputs raw intensities)
  lfq_mat_log2 <- log2(lfq_mat)
  lfq_mat_log2[!is.finite(lfq_mat_log2)] <- NA

  # Strip path and .d extension from sample names
  colnames(lfq_mat_log2) <- gsub("\\.d$", "", basename(colnames(lfq_mat_log2)))

  # Protein meta (genes-equivalent in EList)
  protein_meta <- merge(
    data.frame(ProteinID = rownames(lfq_mat_log2), stringsAsFactors = FALSE),
    as.data.frame(pep_counts),
    by.x = "ProteinID", by.y = "proteins",
    all.x = TRUE
  )
  protein_meta$NPeptides[is.na(protein_meta$NPeptides)] <- 0L
  protein_meta$NSpectra[is.na(protein_meta$NSpectra)]   <- 0L

  message("[DDA] LFQ matrix: ", nrow(lfq_mat_log2), " proteins x ",
          ncol(lfq_mat_log2), " samples")

  list(
    psms         = psms_filtered,
    lfq_wide     = lfq_mat_log2,
    protein_meta = protein_meta
  )
}


#' Build a limma EList from Sage MaxLFQ output
#'
#' @param sage_parsed List from parse_sage_results()
#' @param metadata_df Data.frame with SampleID and Group columns
#' @return EList object compatible with limma DE pipeline
build_dda_elist <- function(sage_parsed, metadata_df) {
  mat  <- sage_parsed$lfq_wide
  meta <- sage_parsed$protein_meta

  # Match sample columns to metadata
  sample_ids <- metadata_df$SampleID
  available  <- intersect(sample_ids, colnames(mat))
  if (length(available) == 0) {
    stop("[DDA] No sample names match between LFQ matrix and metadata. ",
         "Matrix cols: ", paste(head(colnames(mat), 5), collapse = ", "),
         "; Metadata SampleIDs: ", paste(head(sample_ids, 5), collapse = ", "))
  }
  mat <- mat[, available, drop = FALSE]
  metadata_df <- metadata_df[metadata_df$SampleID %in% available, ]

  # Build protein-level gene annotation (Protein.Group for compatibility)
  genes_df <- data.frame(
    Protein.Group = rownames(mat),
    NPeptides     = meta$NPeptides[match(rownames(mat), meta$ProteinID)],
    NSpectra      = meta$NSpectra[match(rownames(mat), meta$ProteinID)],
    stringsAsFactors = FALSE,
    row.names     = rownames(mat)
  )

  elist <- new("EList",
    list(
      E       = mat,
      genes   = genes_df,
      targets = metadata_df
    )
  )

  message("[DDA] Built EList: ", nrow(elist$E), " proteins x ", ncol(elist$E), " samples")
  elist
}


#' Normalize DDA MaxLFQ log2 matrix
#'
#' @param mat log2 protein x sample matrix (with NAs)
#' @param method One of "cyclicloess" (default), "median", "quantile", "none"
#' @return Normalized matrix (same dimensions, NAs preserved)
normalize_dda_matrix <- function(mat, method = "cyclicloess") {
  switch(method,
    "cyclicloess" = {
      limma::normalizeCyclicLoess(mat, method = "fast")
    },
    "median" = {
      sample_medians <- apply(mat, 2, median, na.rm = TRUE)
      global_median  <- median(mat, na.rm = TRUE)
      sweep(mat, 2, sample_medians - global_median, "-")
    },
    "quantile" = {
      if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        warning("[DDA] preprocessCore not available, falling back to median normalization")
        return(normalize_dda_matrix(mat, method = "median"))
      }
      preprocessCore::normalize.quantiles(mat)
    },
    "none" = mat,
    stop("Unknown normalization method: ", method)
  )
}


#' Filter DDA protein matrix by valid value threshold per group
#'
#' @param mat log2 protein x sample matrix (with NAs, post-normalization)
#' @param group_vec Character/factor vector of group assignments (length = ncol(mat))
#' @param min_valid_fraction Minimum fraction of samples in EACH group that must
#'   have a valid (non-NA) value. Default 0.5.
#' @return Filtered matrix (rows passing threshold)
filter_dda_valid_values <- function(mat, group_vec, min_valid_fraction = 0.5) {
  groups <- unique(group_vec)

  passes <- apply(mat, 1, function(row) {
    all(vapply(groups, function(g) {
      group_vals <- row[group_vec == g]
      n_valid    <- sum(!is.na(group_vals))
      n_total    <- length(group_vals)
      (n_valid / n_total) >= min_valid_fraction
    }, logical(1)))
  })

  mat[passes, , drop = FALSE]
}


#' Impute missing values in DDA protein matrix
#'
#' @param mat log2 protein x sample matrix (post-filter, NAs remain)
#' @param method One of "perseus" (default), "minprob", "mindet", "none"
#' @param width Perseus width parameter (default 0.3)
#' @param shift Perseus downshift in SD units (default 1.8)
#' @param q MinProb/MinDet quantile (default 0.01)
#' @return Imputed matrix (no NAs unless method = "none")
impute_dda_matrix <- function(
  mat,
  method = "perseus",
  width  = 0.3,
  shift  = 1.8,
  q      = 0.01
) {
  set.seed(42)  # reproducibility for stochastic imputation

  switch(method,
    "perseus" = {
      apply(mat, 2, function(col) {
        missing <- is.na(col)
        if (!any(missing)) return(col)
        col_mean <- mean(col, na.rm = TRUE)
        col_sd   <- sd(col, na.rm = TRUE)
        if (is.na(col_sd) || col_sd == 0) col_sd <- 1
        col[missing] <- rnorm(
          n    = sum(missing),
          mean = col_mean - shift * col_sd,
          sd   = width * col_sd
        )
        col
      })
    },
    "minprob" = {
      apply(mat, 2, function(col) {
        missing <- is.na(col)
        if (!any(missing)) return(col)
        obs     <- col[!missing]
        center  <- quantile(obs, probs = q, na.rm = TRUE)
        col_sd  <- sd(obs, na.rm = TRUE)
        if (is.na(col_sd) || col_sd == 0) col_sd <- 1
        col[missing] <- rnorm(
          n    = sum(missing),
          mean = center,
          sd   = width * col_sd
        )
        col
      })
    },
    "mindet" = {
      apply(mat, 2, function(col) {
        missing <- is.na(col)
        if (!any(missing)) return(col)
        col[missing] <- quantile(col, probs = q, na.rm = TRUE)
        col
      })
    },
    "none" = mat,
    stop("Unknown imputation method: ", method)
  )
}


#' Run the full DDA pre-DE pipeline: normalize -> filter -> impute -> build EList
#'
#' @param lfq_wide log2 protein x sample matrix from parse_sage_results()
#' @param protein_meta protein metadata data.frame from parse_sage_results()
#' @param metadata_df data.frame with SampleID and Group columns
#' @param norm_method normalization method
#' @param min_valid_fraction valid value filter threshold
#' @param impute_method imputation method
#' @param perseus_width width parameter for Perseus imputation
#' @param perseus_shift shift parameter for Perseus imputation
#' @return List with $elist (EList), $n_prefilter, $n_postfilter
run_dda_pipeline <- function(
  lfq_wide,
  protein_meta,
  metadata_df,
  norm_method        = "cyclicloess",
  min_valid_fraction = 0.5,
  impute_method      = "perseus",
  perseus_width      = 0.3,
  perseus_shift      = 1.8
) {
  # Step 1: Normalize
  message("[DDA Pipeline] Normalizing with method: ", norm_method)
  mat_norm <- normalize_dda_matrix(lfq_wide, method = norm_method)

  # Step 2: Filter valid values
  group_vec <- metadata_df$Group[match(colnames(mat_norm), metadata_df$SampleID)]
  n_before <- nrow(mat_norm)
  mat_filtered <- filter_dda_valid_values(mat_norm, group_vec, min_valid_fraction)
  n_after <- nrow(mat_filtered)
  message(sprintf("[DDA Pipeline] Valid value filter: %d -> %d proteins (%.0f%% retained)",
    n_before, n_after, 100 * n_after / max(n_before, 1)))

  # Step 3: Impute
  message("[DDA Pipeline] Imputing with method: ", impute_method)
  mat_imputed <- impute_dda_matrix(mat_filtered,
    method = impute_method,
    width  = perseus_width,
    shift  = perseus_shift
  )

  # Step 4: Build EList
  # Update protein_meta to match filtered rows
  filtered_meta <- protein_meta[protein_meta$ProteinID %in% rownames(mat_imputed), ]
  elist <- build_dda_elist(
    list(lfq_wide = mat_imputed, protein_meta = filtered_meta),
    metadata_df
  )

  list(
    elist       = elist,
    n_prefilter = n_before,
    n_postfilter = n_after
  )
}


#' Generate sbatch script for Sage DDA search on Hive
#'
#' @param sage_bin Path to Sage binary on HPC
#' @param config_path Path to sage.json config on HPC
#' @param raw_dir Directory containing .d files on HPC
#' @param output_dir Output directory on HPC
#' @param experiment_name Name for the SLURM job
#' @param cpus Number of CPUs (default 32)
#' @param mem_gb Memory in GB (default 64)
#' @param time_limit Time limit string (default "02:00:00")
#' @param account SLURM account (default "genome-center-grp")
#' @param partition SLURM partition (default "high")
#' @return Character string: sbatch script content
generate_sage_sbatch <- function(
  sage_bin,
  config_path,
  raw_dir,
  output_dir,
  experiment_name = "sage_search",
  cpus            = 32,
  mem_gb          = 64,
  time_limit      = "02:00:00",
  account         = "genome-center-grp",
  partition       = "high"
) {
  # Sanitize experiment name for SLURM

  safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", experiment_name)

  script <- paste0(
'#!/bin/bash
#SBATCH --job-name=delimp_sage_', safe_name, '
#SBATCH --partition=', partition, '
#SBATCH --account=', account, '
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=', cpus, '
#SBATCH --mem=', mem_gb, 'G
#SBATCH --time=', time_limit, '
#SBATCH --output="', output_dir, '/logs/sage_%j.out"
#SBATCH --error="', output_dir, '/logs/sage_%j.err"

set -euo pipefail
echo "[DE-LIMP Sage] Job start: $(date)"
echo "[DE-LIMP Sage] Node: $(hostname)"
echo "[DE-LIMP Sage] CPUs: ', cpus, ', Memory: ', mem_gb, 'G"

SAGE_BIN="', sage_bin, '"
CONFIG="', config_path, '"
RAW_DIR="', raw_dir, '"
OUTPUT_DIR="', output_dir, '"

# Verify Sage binary
if [ ! -x "$SAGE_BIN" ]; then
  echo "[ERROR] Sage binary not found or not executable: $SAGE_BIN"
  exit 1
fi

# Collect .d directories
D_FILES=$(find "$RAW_DIR" -maxdepth 1 -name "*.d" -type d | sort)
N_FILES=$(echo "$D_FILES" | grep -c "." || true)
echo "[DE-LIMP Sage] Found $N_FILES .d files in $RAW_DIR"

if [ "$N_FILES" -eq 0 ]; then
  echo "[ERROR] No .d directories found in $RAW_DIR"
  exit 1
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run Sage
echo "[DE-LIMP Sage] Starting search..."
"$SAGE_BIN" \\
  --write-pin \\
  --output_directory "$OUTPUT_DIR" \\
  "$CONFIG" \\
  $D_FILES

echo "[DE-LIMP Sage] Search complete: $(date)"
echo "[DE-LIMP Sage] Output files:"
ls -lh "$OUTPUT_DIR"/*.tsv "$OUTPUT_DIR"/*.json 2>/dev/null || echo "(no output files found)"
echo "[DE-LIMP Sage] Done."
')
  script
}


#' Parse sage_report.json summary statistics
#'
#' @param json_path Path to sage_report.json or results.json
#' @return Named list: psms, unique_peptides, unique_proteins, sage_version
parse_sage_report <- function(json_path) {
  if (!file.exists(json_path)) {
    message("[DDA] sage report JSON not found: ", json_path)
    return(NULL)
  }

  tryCatch({
    report <- jsonlite::fromJSON(json_path)
    list(
      psms             = report$psms %||% NA_integer_,
      unique_peptides  = report$unique_peptides %||% NA_integer_,
      unique_proteins  = report$unique_proteins %||% NA_integer_,
      sage_version     = report$sage_version %||% report$version %||% "unknown",
      files            = report$files %||% NA_integer_
    )
  }, error = function(e) {
    message("[DDA] Failed to parse sage report: ", e$message)
    NULL
  })
}


#' Compute DDA QC metrics from PSM table
#'
#' @param psms Filtered PSM data.table from parse_sage_results()
#' @param lfq_wide log2 protein x sample matrix
#' @return Named list of QC metrics
compute_dda_qc_metrics <- function(psms, lfq_wide) {
  n_psms       <- nrow(psms)
  n_peptides   <- data.table::uniqueN(psms$peptide)
  n_proteins   <- nrow(lfq_wide)

  # Median peptides per protein
  pep_per_prot <- psms[, .(n = data.table::uniqueN(peptide)), by = proteins]
  med_pep_per_prot <- median(pep_per_prot$n, na.rm = TRUE)

  # Missed cleavages
  mc_pattern <- "[KR][^P]"  # tryptic missed cleavage sites
  n_mc <- sum(grepl(mc_pattern, psms$peptide))
  pct_missed_cleavage <- 100 * n_mc / max(n_psms, 1)

  # Precursor mass error (ppm) -- column may be called ppm_difference or similar
  ppm_col <- intersect(c("expmass_ppm", "ppm_difference", "delta_mass_ppm"), colnames(psms))
  if (length(ppm_col) > 0) {
    mass_error_median <- median(abs(psms[[ppm_col[1]]]), na.rm = TRUE)
  } else {
    mass_error_median <- NA_real_
  }

  list(
    n_psms              = n_psms,
    n_peptides          = n_peptides,
    n_proteins          = n_proteins,
    med_pep_per_prot    = med_pep_per_prot,
    pct_missed_cleavage = round(pct_missed_cleavage, 1),
    mass_error_ppm      = round(mass_error_median, 2)
  )
}
