# Setup: source pure helper files for testing
# These have no Shiny dependency — safe to load standalone

# Find project root: go up from tests/testthat/ to project root
project_root <- normalizePath(file.path(getwd(), "..", ".."))

# Provide %||% operator (normally from rlang, loaded by Shiny)
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x

# Load helpers.R (cal_z_score, detect_organism_db)
source(file.path(project_root, "R", "helpers.R"))

# Load helpers_search.R (build_diann_flags, parse_sbatch_output, etc.)
source(file.path(project_root, "R", "helpers_search.R"))

# Load helpers_dda.R (build_dda_canonical_peptide, dda_blast_species,
# build_denovo_master, contaminant/skin-hair filters) — pure functions
source(file.path(project_root, "R", "helpers_dda.R"))

# Load helpers_phospho.R (parse_phospho_positions)
source(file.path(project_root, "R", "helpers_phospho.R"))

# Load proteogenomics helpers (no Shiny reactivity; pure functions)
source(file.path(project_root, "R", "helpers_proteogenomics.R"))
source(file.path(project_root, "R", "helpers_proteog_assembly.R"))
source(file.path(project_root, "R", "helpers_slims.R"))
source(file.path(project_root, "R", "helpers_rnaseq.R"))
source(file.path(project_root, "R", "helpers_proteog_qc.R"))
