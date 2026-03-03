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

# Load helpers_phospho.R (parse_phospho_positions)
source(file.path(project_root, "R", "helpers_phospho.R"))
