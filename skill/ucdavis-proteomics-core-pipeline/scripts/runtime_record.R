# =============================================================================
# runtime_record.R -- the R a DE ran in, as the process that ran it sees it: where R lives, its
# Rscript, the libraries it loaded packages from, the conda env and the container around it.
# de_provenance.json `runtime`; provenance.py builds the reproducibility bundle's environment
# from it (and sessionInfo.txt beside it). The ONE definition, sourced by run_de.R and
# run_silac.R: run_silac.R once wrote its own without `rscript`, and provenance.py told the user
# the DE had run in "a different environment" than setup.json's when it had not (2026-10).
# The container variables are the ones DE-LIMP's app.R checks (Apptainer/Singularity set them;
# Docker leaves /.dockerenv).
# =============================================================================

runtime_record <- function(writer) {
  env1 <- function(...) { for (v in c(...)) { x <- Sys.getenv(v, ""); if (nzchar(x)) return(x) }; NULL }
  rs <- file.path(R.home("bin"), "Rscript")
  list(r_version = R.version.string, r_home = R.home(),
       rscript = if (file.exists(rs)) rs else NULL,
       lib_paths = as.list(.libPaths()),
       conda_prefix = env1("CONDA_PREFIX"),
       container = env1("APPTAINER_CONTAINER", "SINGULARITY_CONTAINER"),
       container_name = env1("APPTAINER_NAME", "SINGULARITY_NAME"),
       docker = file.exists("/.dockerenv"),
       session_info = "sessionInfo.txt",
       note = sprintf("recorded by %s in the process that ran this DE", writer))
}
