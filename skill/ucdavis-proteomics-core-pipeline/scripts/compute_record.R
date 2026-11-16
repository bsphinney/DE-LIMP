# =============================================================================
# compute_record.R -- two small definitions shared by run_de.R and make_screen_data.R (repo
# rule 3: one definition each), sourced next to skill_version.R:
#   slurm_cpu_family()  the SLURM node's name, its AvailableFeatures and the CPU-family feature
#                       (what `sbatch --constraint=` selects); DPC-Quant results reproduce bit
#                       for bit only within one family (references/reproducibility.md)
#   file_sha256(path)   sha256 of a file; NA when it is missing or cannot be hashed
# =============================================================================

CPU_FAMILY_REGEX <- "^(zen[0-9]*|icelake|cascadelake|skylake|sapphirerapids|emeraldrapids|haswell|broadwell)$"

# list(node = "" off SLURM, feats = character(), fam = the family feature or character(0)).
slurm_cpu_family <- function() {
  sh <- function(cmd, args) tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
                                     error = function(e) character(0))
  node <- Sys.getenv("SLURMD_NODENAME", "")
  feats <- if (nzchar(node) && nzchar(Sys.which("scontrol"))) {
    out <- sh("scontrol", c("show", "node", node))
    f <- regmatches(out, regexpr("AvailableFeatures=[^ ]*", out))
    if (length(f)) strsplit(sub("^AvailableFeatures=", "", f[1]), ",")[[1]] else character(0)
  } else character(0)
  list(node = node, feats = feats, fam = grep(CPU_FAMILY_REGEX, feats, value = TRUE))
}

file_sha256 <- function(path) tryCatch(
  if (!file.exists(path)) NA_character_
  else if (requireNamespace("digest", quietly = TRUE))
    digest::digest(path, algo = "sha256", file = TRUE)
  else if (nzchar(Sys.which("sha256sum")))
    sub(" .*", "", system2("sha256sum", shQuote(path), stdout = TRUE)[1])
  else if (nzchar(Sys.which("shasum")))
    sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)[1])
  else NA_character_,
  error = function(e) NA_character_)
