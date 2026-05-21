#!/usr/bin/env Rscript
# Phase C step 2 validation: assemble_proteogenomics_fasta() against real
# inputs (UniProt mouse OPG + Phase C predicted ORFs + Mouse_Tissue_Contaminants).
#
# History: The original version of this script used
#   system2("grep", c("-c", "^>", f), stdout = TRUE)
# which bash parsed as `grep -c ^> /path/to/f` — interpreting `^>` as
# caret + redirect operator. The redirect truncated each input FASTA
# to a 2-byte "0\n" file before the assembly function was even called.
# See NOTES_spec_lessons #16 for the full anti-pattern.
#
# This version uses R native (grepl + readLines) per lesson #16: "prefer
# R native operations over system2() shell-outs."

# Adjust ROOT to where the helpers actually live on Hive.
ROOT  <- "/quobyte/proteomics-grp/de-limp/scratch/proteog_orchestrator"

if (!exists("%||%")) `%||%` <- function(x, y) if (is.null(x)) y else x
source(file.path(ROOT, "helpers.R"))
source(file.path(ROOT, "helpers_search.R"))
source(file.path(ROOT, "helpers_proteogenomics.R"))
source(file.path(ROOT, "helpers_proteog_assembly.R"))
source(file.path(ROOT, "helpers_slims.R"))
source(file.path(ROOT, "helpers_rnaseq.R"))
source(file.path(ROOT, "helpers_proteog_qc.R"))

UNIPROT  <- "/quobyte/proteomics-grp/de-limp/fasta/UP000000589_mus_musculus_opg_2026_05.fasta"
PREDORF  <- "/quobyte/proteomics-grp/de-limp/rnaseq/phase_c_proof_of_life/predicted_orfs.fasta"
MERGEDGTF<- "/quobyte/proteomics-grp/de-limp/rnaseq/phase_c_proof_of_life/stringtie_out/merged.gtf"
CONTAM   <- "/quobyte/proteomics-grp/de-limp/DE-LIMP/contaminants/Mouse_Tissue_Contaminants.fasta"
OUTDIR   <- "/quobyte/proteomics-grp/de-limp/databases/proteogenomics"

#' Count FASTA entries using R native — lesson #16 compliant
count_fasta_entries <- function(path) {
  if (!file.exists(path)) return(NA_integer_)
  sum(grepl("^>", readLines(path)))
}

cat("=== Input pre-flight ===\n")
inputs <- list(UniProt = UNIPROT, PredictedORFs = PREDORF, Contaminants = CONTAM)
for (lbl in names(inputs)) {
  f <- inputs[[lbl]]
  n <- count_fasta_entries(f)
  cat(sprintf("  %-15s %12s bytes  %s entries\n",
              lbl,
              format(file.info(f)$size, big.mark = ","),
              format(n, big.mark = ",")))
}

cat("\n=== Running assemble_proteogenomics_fasta() ===\n")
t0 <- Sys.time()
result <- assemble_proteogenomics_fasta(
  project_name         = "phase_c_assembly_test",
  uniprot_fasta        = UNIPROT,
  predicted_orfs_fasta = PREDORF,
  merged_gtf           = MERGEDGTF,
  contaminants_fasta   = CONTAM,
  output_dir           = OUTDIR,
  dedupe               = TRUE,
  build_metadata       = list(
    organism = "Mus musculus", build = "GRCm39",
    rnaseq_n_samples = 2L, read_length_tier = "significantly_relaxed",
    uniprot_release = "UP000000589 OPG 2026-05",
    contaminants_used = TRUE, contaminants_source = "HaoGroup"
  )
)
elapsed <- difftime(Sys.time(), t0, units = "secs")
cat(sprintf("\nElapsed: %.1f sec\n", as.numeric(elapsed)))

cat("\n=== Returned manifest ===\n")
cat("pipeline_id:    ", result$pipeline_id, "\n")
cat("path:           ", result$path, "\n")
cat("merged_gtf_path:", result$merged_gtf_path, "\n")
cat("composition:\n")
for (k in names(result$composition)) {
  cat(sprintf("  %-15s %s\n", k, format(result$composition[[k]], big.mark = ",")))
}
cat("\nmethods_paragraph:\n  ", result$methods_paragraph, "\n")

cat("\n=== Output FASTA on disk ===\n")
fi <- file.info(result$path)
cat(sprintf("  Size: %s bytes\n", format(fi$size, big.mark = ",")))
n_out <- count_fasta_entries(result$path)
cat(sprintf("  Entries (counted via R native): %s\n", format(n_out, big.mark = ",")))
stopifnot(n_out == result$composition$total)

cat("\n=== Merged GTF copied alongside ===\n")
cat("  Path:", result$merged_gtf_path, "\n")
cat("  Size:", format(file.info(result$merged_gtf_path)$size, big.mark = ","), "bytes\n")
stopifnot(file.exists(result$merged_gtf_path))

cat("\n=== Registry entry ===\n")
reg <- load_proteog_registry()
entry <- reg[["phase_c_assembly_test"]]
stopifnot(!is.null(entry))
cat("  registry path:        ", entry$path, "\n")
cat("  merged_gtf_path field:", entry$merged_gtf_path %||% "MISSING", "\n")
cat("  composition$total:    ", entry$composition$total, "\n")
cat("  created:              ", entry$created, "\n")
cat("  pipeline_version:     ", entry$pipeline_version, "\n")
stopifnot(!is.null(entry$merged_gtf_path))
stopifnot(entry$composition$total == result$composition$total)

cat("\n=== Sanity check: source class breakdown ===\n")
expected_predicted <- result$composition$REF + result$composition$NOVEL_ISOFORM +
                       result$composition$NOVEL_GENE
cat(sprintf("  Predicted-ORF total (REF+ISO+GENE):  %s\n",
            format(expected_predicted, big.mark = ",")))
cat(sprintf("  UNIPROT bucket (UniProt + Contam):    %s\n",
            format(result$composition$UNIPROT, big.mark = ",")))
cat(sprintf("  UNPARSED:                             %s (must be 0)\n",
            result$composition$UNPARSED))
stopifnot(result$composition$UNPARSED == 0)

cat("\n=========================================\n")
cat("ASSEMBLY VALIDATION PASSED\n")
cat("=========================================\n")
