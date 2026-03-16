# ============================================================================
# MOFA2 Test Data Setup for DE-LIMP
# ============================================================================
# This script downloads and prepares multi-omics test datasets for testing
# the MOFA2 integration in DE-LIMP.
#
# Three test datasets are provided:
#   1. TCGA Breast Cancer (mixOmics) - mRNA + miRNA + Proteomics (3 views)
#   2. Simulated Proteomics + Phosphoproteomics (2 views, realistic for DE-LIMP)
#   3. MOFA2 built-in test data (2 views)
#
# Output: CSV files ready for upload to DE-LIMP MOFA tab
# ============================================================================

# Install required packages if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages_needed <- c("mixOmics", "MOFA2")
for (pkg in packages_needed) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE)
  }
}

# Create output directory
output_dir <- "mofa_test_data"
dir.create(output_dir, showWarnings = FALSE)

cat("\n============================================================\n")
cat("MOFA2 Test Data Setup for DE-LIMP\n")
cat("============================================================\n\n")

# ============================================================================
# DATASET 1: TCGA Breast Cancer Multi-Omics (mixOmics)
# ============================================================================
# This is a real dataset with:
#   - 150 breast cancer samples
#   - 3 cancer subtypes: Basal, Her2, Luminal A
#   - 3 omics types: mRNA (520 genes), miRNA (184), Proteomics (142 proteins)
# ============================================================================

cat("=== Dataset 1: TCGA Breast Cancer Multi-Omics ===\n")

# Download from mixOmics GitHub
tcga_url <- "https://github.com/mixOmicsTeam/mixOmics/raw/master/data/breast.TCGA.rda"
tcga_file <- file.path(output_dir, "breast.TCGA.rda")

if (!file.exists(tcga_file)) {
  cat("Downloading TCGA breast cancer data...\n")
  download.file(tcga_url, tcga_file, mode = "wb")
}

load(tcga_file)

# Extract training data (has all 3 omics)
mrna_data <- breast.TCGA$data.train$mrna      # 150 samples × 520 genes
mirna_data <- breast.TCGA$data.train$mirna    # 150 samples × 184 miRNAs  
protein_data <- breast.TCGA$data.train$protein # 150 samples × 142 proteins
subtypes <- breast.TCGA$data.train$subtype    # Cancer subtypes

cat(sprintf("  mRNA: %d samples × %d features\n", nrow(mrna_data), ncol(mrna_data)))
cat(sprintf("  miRNA: %d samples × %d features\n", nrow(mirna_data), ncol(mirna_data)))
cat(sprintf("  Protein: %d samples × %d features\n", nrow(protein_data), ncol(protein_data)))
cat(sprintf("  Subtypes: %s\n", paste(table(subtypes), collapse=", ")))

# Transpose to features × samples format for MOFA/DE-LIMP
# (DE-LIMP expects features as rows, samples as columns)
mrna_matrix <- t(mrna_data)
mirna_matrix <- t(mirna_data)
protein_matrix <- t(protein_data)

# Save as CSV files (ready for DE-LIMP upload)
write.csv(mrna_matrix, file.path(output_dir, "tcga_breast_mRNA.csv"))
write.csv(mirna_matrix, file.path(output_dir, "tcga_breast_miRNA.csv"))
write.csv(protein_matrix, file.path(output_dir, "tcga_breast_protein.csv"))

# Save metadata
metadata <- data.frame(
  Sample = rownames(mrna_data),
  Subtype = as.character(subtypes),
  stringsAsFactors = FALSE
)
write.csv(metadata, file.path(output_dir, "tcga_breast_metadata.csv"), row.names = FALSE)

cat("  Saved: tcga_breast_mRNA.csv, tcga_breast_miRNA.csv, tcga_breast_protein.csv\n")
cat("  Saved: tcga_breast_metadata.csv\n\n")


# ============================================================================
# DATASET 2: Simulated Global Proteomics + Phosphoproteomics
# ============================================================================
# This simulates a realistic DE-LIMP use case:
#   - 40 samples (20 Control, 20 Treatment)
#   - Global proteomics: 2000 proteins
#   - Phosphoproteomics: 5000 phosphosites
#   - Some phospho changes are abundance-driven, others are true regulation
# ============================================================================

cat("=== Dataset 2: Simulated Proteomics + Phosphoproteomics ===\n")

set.seed(42)

n_samples <- 40
n_proteins <- 2000
n_phosphosites <- 5000

# Sample names and groups
sample_names <- paste0("Sample_", sprintf("%02d", 1:n_samples))
groups <- rep(c("Control", "Treatment"), each = n_samples/2)

# --- Global Proteomics ---
# Base expression + group effect + noise
base_expr <- rnorm(n_proteins, mean = 20, sd = 3)
group_effect <- c(rep(0, n_samples/2), rep(1, n_samples/2))  # Treatment effect

# 10% of proteins respond to treatment
responding_proteins <- sample(1:n_proteins, size = n_proteins * 0.10)
protein_matrix <- matrix(NA, nrow = n_proteins, ncol = n_samples)

for (i in 1:n_proteins) {
  effect_size <- ifelse(i %in% responding_proteins, rnorm(1, 0, 1.5), 0)
  protein_matrix[i, ] <- base_expr[i] + effect_size * group_effect + rnorm(n_samples, 0, 0.5)
}

rownames(protein_matrix) <- paste0("PROT_", sprintf("%04d", 1:n_proteins))
colnames(protein_matrix) <- sample_names

# --- Phosphoproteomics ---
# Three types of phosphosites:
#   1. Abundance-driven (follow their parent protein)
#   2. Truly regulated (independent of protein level)
#   3. Noise (no real signal)

phospho_matrix <- matrix(NA, nrow = n_phosphosites, ncol = n_samples)
site_info <- data.frame(
  SiteID = paste0("SITE_", sprintf("%05d", 1:n_phosphosites)),
  ParentProtein = paste0("PROT_", sprintf("%04d", sample(1:n_proteins, n_phosphosites, replace = TRUE))),
  Type = sample(c("abundance_driven", "regulated", "noise"), n_phosphosites, 
                replace = TRUE, prob = c(0.3, 0.2, 0.5)),
  stringsAsFactors = FALSE
)

for (i in 1:n_phosphosites) {
  parent_idx <- as.numeric(gsub("PROT_", "", site_info$ParentProtein[i]))
  parent_values <- protein_matrix[parent_idx, ]
  
  if (site_info$Type[i] == "abundance_driven") {
    # Phosphosite follows parent protein + small noise
    phospho_matrix[i, ] <- parent_values + rnorm(n_samples, 0, 0.3)
    
  } else if (site_info$Type[i] == "regulated") {
    # Independent regulation (could go opposite direction from protein)
    reg_effect <- rnorm(1, 0, 2)  # Stronger effect
    phospho_matrix[i, ] <- mean(parent_values) + reg_effect * group_effect + rnorm(n_samples, 0, 0.4)
    
  } else {
    # Noise - no real signal
    phospho_matrix[i, ] <- rnorm(n_samples, mean = 20, sd = 1)
  }
}

rownames(phospho_matrix) <- site_info$SiteID
colnames(phospho_matrix) <- sample_names

cat(sprintf("  Global proteomics: %d proteins × %d samples\n", nrow(protein_matrix), ncol(protein_matrix)))
cat(sprintf("  Phosphoproteomics: %d sites × %d samples\n", nrow(phospho_matrix), ncol(phospho_matrix)))
cat(sprintf("  Site types: %s\n", paste(names(table(site_info$Type)), table(site_info$Type), sep="=", collapse=", ")))

# Save as CSV
write.csv(protein_matrix, file.path(output_dir, "sim_global_proteomics.csv"))
write.csv(phospho_matrix, file.path(output_dir, "sim_phosphoproteomics.csv"))
write.csv(site_info, file.path(output_dir, "sim_phosphosite_info.csv"), row.names = FALSE)

# Save metadata
sim_metadata <- data.frame(
  Sample = sample_names,
  Group = groups,
  stringsAsFactors = FALSE
)
write.csv(sim_metadata, file.path(output_dir, "sim_metadata.csv"), row.names = FALSE)

cat("  Saved: sim_global_proteomics.csv, sim_phosphoproteomics.csv\n")
cat("  Saved: sim_phosphosite_info.csv, sim_metadata.csv\n\n")


# ============================================================================
# DATASET 3: MOFA2 Built-in Test Data
# ============================================================================
# Simple 2-view dataset included with MOFA2 package
# Good for quick testing of the integration
# ============================================================================

cat("=== Dataset 3: MOFA2 Built-in Test Data ===\n")

# Load MOFA2 test data
filepath <- system.file("extdata", "test_data.RData", package = "MOFA2")
if (file.exists(filepath)) {
  load(filepath)
  
  # 'dt' is a long-format data.table
  cat(sprintf("  Format: Long data.table with %d rows\n", nrow(dt)))
  cat(sprintf("  Views: %s\n", paste(unique(dt$view), collapse=", ")))
  cat(sprintf("  Samples: %d\n", length(unique(dt$sample))))
  
  # Convert to wide matrices for each view
  library(data.table)
  
  for (view_name in unique(dt$view)) {
    view_data <- dt[view == view_name]
    wide <- dcast(view_data, feature ~ sample, value.var = "value")
    mat <- as.matrix(wide[, -1, with = FALSE])
    rownames(mat) <- wide$feature
    
    filename <- sprintf("mofa2_builtin_%s.csv", gsub(" ", "_", view_name))
    write.csv(mat, file.path(output_dir, filename))
    cat(sprintf("  Saved: %s (%d features × %d samples)\n", filename, nrow(mat), ncol(mat)))
  }
} else {
  cat("  MOFA2 test data not found (package may need installation)\n")
}

cat("\n")

# ============================================================================
# DATASET 4: Create RDS files for testing DE-LIMP session import
# ============================================================================

cat("=== Dataset 4: RDS Files for Testing DE-LIMP Import ===\n")

# Simulate a DE-LIMP session RDS (what users would save from another analysis)
# This tests the RDS import functionality

# Create a mock DE-LIMP session with the simulated proteomics data
mock_delimp_session <- list(
  raw_data = protein_matrix,
  metadata = data.frame(
    Run = sample_names,
    Group = groups,
    stringsAsFactors = FALSE
  ),
  session_info = list(
    app_version = "2.1.1",
    timestamp = Sys.time(),
    description = "Mock DE-LIMP session for MOFA testing"
  )
)

# Add a mock limma fit object
library(limma)
design <- model.matrix(~ 0 + factor(groups))
colnames(design) <- c("Control", "Treatment")
fit <- lmFit(protein_matrix, design)
contrast_matrix <- makeContrasts(Treatment - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

mock_delimp_session$fit <- fit2

saveRDS(mock_delimp_session, file.path(output_dir, "mock_delimp_session.rds"))
cat("  Saved: mock_delimp_session.rds (full DE-LIMP session with fit object)\n")

# Save just a matrix as RDS
saveRDS(phospho_matrix, file.path(output_dir, "phospho_matrix_only.rds"))
cat("  Saved: phospho_matrix_only.rds (simple matrix)\n")

# Save as a named list (another common format)
saveRDS(
  list(
    matrix = phospho_matrix,
    metadata = sim_metadata,
    source = "Phosphoproteomics analysis"
  ),
  file.path(output_dir, "phospho_with_metadata.rds")
)
cat("  Saved: phospho_with_metadata.rds (named list with matrix)\n")

cat("\n")

# ============================================================================
# Summary
# ============================================================================

cat("============================================================\n")
cat("TEST DATA SETUP COMPLETE\n")
cat("============================================================\n\n")

cat("Files created in:", normalizePath(output_dir), "\n\n")

cat("DATASET 1: TCGA Breast Cancer (Real multi-omics)\n")
cat("  - tcga_breast_mRNA.csv (520 genes × 150 samples)\n")
cat("  - tcga_breast_miRNA.csv (184 miRNAs × 150 samples)\n")
cat("  - tcga_breast_protein.csv (142 proteins × 150 samples)\n")
cat("  - tcga_breast_metadata.csv (sample subtypes)\n")
cat("  USE CASE: True 3-view multi-omics integration\n\n")

cat("DATASET 2: Simulated Proteomics + Phospho (Realistic for DE-LIMP)\n")
cat("  - sim_global_proteomics.csv (2000 proteins × 40 samples)\n")
cat("  - sim_phosphoproteomics.csv (5000 sites × 40 samples)\n")
cat("  - sim_metadata.csv (Control vs Treatment)\n")
cat("  USE CASE: Deconvolute abundance vs true phospho-regulation\n\n")

cat("DATASET 3: MOFA2 Built-in (Quick testing)\n")
cat("  - mofa2_builtin_view_0.csv\n")
cat("  - mofa2_builtin_view_1.csv\n")
cat("  USE CASE: Verify MOFA training works\n\n")

cat("DATASET 4: RDS Files (Test import formats)\n")
cat("  - mock_delimp_session.rds (full session with fit object)\n")
cat("  - phospho_matrix_only.rds (simple matrix)\n")
cat("  - phospho_with_metadata.rds (named list)\n")
cat("  USE CASE: Test RDS import for second proteomics dataset\n\n")

cat("============================================================\n")
cat("TESTING WORKFLOW\n")
cat("============================================================\n\n")

cat("1. Start with Dataset 2 (simulated) for basic testing:\n")
cat("   - Load sim_global_proteomics.csv as View 1 (or use current DE-LIMP analysis)\n")
cat("   - Load sim_phosphoproteomics.csv as View 2\n")
cat("   - Train MOFA and look for factors that are:\n")
cat("     * Protein-heavy (abundance-driven changes)\n")
cat("     * Phospho-heavy (true kinase/signaling regulation)\n")
cat("     * Shared (coordinated response)\n\n")

cat("2. Test RDS import with Dataset 4:\n")
cat("   - Load mock_delimp_session.rds as View 2\n")
cat("   - Verify it extracts the matrix AND the fit object\n")
cat("   - Check Factor-DE correlation can use both fits\n\n")

cat("3. Test true multi-omics with Dataset 1 (TCGA):\n")
cat("   - Load all 3 views (mRNA, miRNA, protein)\n")
cat("   - Sample matching should find 150 common samples\n")
cat("   - Factors should separate cancer subtypes\n\n")

cat("============================================================\n")
