#!/usr/bin/env Rscript
# Phase C dry-run: validate submit_proteogenomics_build() generates the right
# sbatch chain without actually submitting to SLURM.

project_root <- "/Users/brettphinney/Documents/claude"
if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path(project_root, "R", "helpers.R"))
source(file.path(project_root, "R", "helpers_search.R"))
source(file.path(project_root, "R", "helpers_proteogenomics.R"))
source(file.path(project_root, "R", "helpers_proteog_assembly.R"))
source(file.path(project_root, "R", "helpers_slims.R"))
source(file.path(project_root, "R", "helpers_rnaseq.R"))
source(file.path(project_root, "R", "helpers_proteog_qc.R"))
source(file.path(project_root, "R", "server_proteog_builder.R"))

# Mock sbatch: capture (script, dep) tuples; return synthetic incrementing jids.
captured <- list()
fake_id <- 0L
.sbatch_submit <- function(script_path, dep_jid = NULL) {
  fake_id <<- fake_id + 1L
  captured[[length(captured) + 1L]] <<- list(
    script = script_path,
    dep    = dep_jid %||% NA_character_
  )
  sprintf("%d", 90000L + fake_id)
}

# Mock load_reference_registry() — return a synthetic mm39 entry pointing at
# real Hive paths so script generation matches what would happen on Hive.
load_reference_registry <- function() {
  list(
    mm39_GRCm39 = list(
      organism      = "Mus musculus",
      build         = "GRCm39",
      annotation_source = "GENCODE",
      annotation_release = "vM38 basic",
      genome_fasta  = "/quobyte/proteomics-grp/de-limp/references/genomes/mm39_GRCm39_genome.fa",
      star_index    = "/quobyte/proteomics-grp/de-limp/references/star_index/mm39",
      gtf           = "/quobyte/proteomics-grp/de-limp/references/gtf/mm39.gtf",
      rrna_index    = "/quobyte/proteomics-grp/de-limp/references/rrna_index/mm39/rrna_mm39",
      completeness  = "complete"
    )
  )
}

# Create a fake rnaseq_dir with synthetic FASTQs (92bp, gzipped — to trigger
# significantly_relaxed tier matching the May 20 data).
fake_dir <- tempfile("fake_rnaseq_")
dir.create(fake_dir, recursive = TRUE)
for (s in c("SRR_FAKE1", "SRR_FAKE2")) {
  for (rd in c("R1", "R2")) {
    fq <- file.path(fake_dir, sprintf("%s_%s.fastq.gz", s, rd))
    con <- gzfile(fq, "wt")
    for (i in 1:10) {
      writeLines(c(
        sprintf("@read_%d", i),
        paste(rep("A", 92), collapse = ""),
        "+",
        paste(rep("I", 92), collapse = "")
      ), con)
    }
    close(con)
  }
}

# Project dir under a temp dir so we don't write to /quobyte from the laptop.
fake_rnaseq_root <- tempfile("fake_root_")
dir.create(fake_rnaseq_root, recursive = TRUE)

cat("=== DRY-RUN submit_proteogenomics_build ===\n")
result <- submit_proteogenomics_build(
  project_name   = "dryrun_test",
  rnaseq_dir     = fake_dir,
  reference_key  = "mm39_GRCm39",
  sample_names   = c("SRR_FAKE1", "SRR_FAKE2"),
  library_type   = "polyA",
  strand_flag    = "--rf",
  project_tag    = "DRYRUN",
  rnaseq_root    = fake_rnaseq_root
)

cat("\n=== Returned manifest ===\n")
cat("pipeline_id:", result$pipeline_id, "\n")
cat("project_dir:", result$project_dir, "\n")
cat("tier:", result$tier_params$tier, "(read length", result$build_metadata$detected_read_length, "bp)\n")
cat("methods_paragraph:\n")
cat("  ", result$methods_paragraph, "\n\n")

cat("=== Sbatch chain submitted (in order) ===\n")
for (i in seq_along(captured)) {
  cat(sprintf("  [%d] script=%-50s dep=%s\n",
              i, basename(captured[[i]]$script), captured[[i]]$dep))
}

cat("\n=== Files written to project_dir ===\n")
ls_out <- list.files(result$project_dir, recursive = TRUE)
for (f in ls_out) cat("  ", f, "\n")

cat("\n=== status.json contents ===\n")
status <- jsonlite::read_json(result$status_path)
cat("current_stage:", status$current_stage, "\n")
cat("stages:\n")
for (s in status$stages) {
  cat(sprintf("  %-12s status=%-8s job_id=%s\n", s$stage, s$status, s$job_id))
}

# Assertions
stopifnot(length(captured) == 10)                              # 10 sbatch submissions
stopifnot(is.na(captured[[1]]$dep))                            # first job has no dep
for (i in 2:10) {
  prev_jid <- sprintf("%d", 90000L + i - 1L)
  stopifnot(captured[[i]]$dep == prev_jid)                     # chained correctly
}
stopifnot(file.exists(result$status_path))
stopifnot(length(status$stages) == 11)                         # 11 stages declared (includes "assemble" placeholder)
# All 10 sbatch scripts on disk
sbatch_files <- list.files(file.path(result$project_dir, "sbatch"))
stopifnot(length(sbatch_files) == 10)
# Symlinks present for rnaseq inputs
fastq_links <- list.files(file.path(result$project_dir, "rnaseq"))
stopifnot(length(fastq_links) == 4)  # 2 samples × R1+R2

# Sample one of the actual sbatch script contents
star_script <- readLines(file.path(result$project_dir, "sbatch", "star.sbatch"))
stopifnot(any(grepl("--outFilterScoreMinOverLread 0.30", star_script)))  # 92bp tier
stopifnot(any(grepl("module load star/2.7.11a", star_script)))

# Tier choice was significantly_relaxed (92bp)
stopifnot(result$tier_params$tier == "significantly_relaxed")
stopifnot(result$tier_params$qc_gate_unique_pct == 25)

cat("\n=== DRY-RUN VALIDATION PASSED ===\n")
