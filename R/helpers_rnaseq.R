# helpers_rnaseq.R — Sbatch generators for the proteogenomics RNA-seq pipeline.
# No Shiny reactivity. Pure functions; each takes a project config and returns
# an sbatch script as a character string. Phase C composes and submits them.
#
# Pipeline order (each function below corresponds to one stage):
#   1) fastp       array, per sample
#   2) rrna_filter array, per sample (bowtie2 vs organism rRNA)
#   3) star        array, per sample (with adaptive thresholds)
#   3a) qc_gate    single job — halts dependency chain if unique% below tier
#   4) stringtie   array, per sample
#   5) merge       single (stringtie --merge)
#   5a) gffcompare single (class codes for NOVEL_ISOFORM detection)
#   6) gffread     single (extract transcript FASTA)
#   7) transdecoder single (LongOrfs + Predict)
#   8) rewrite     single (header rewriter Python)
#
# All sbatch templates use glue::glue() with .open="<<"/.close=">>" so the
# shell's ${VAR} expansion can pass through unmolested. Use <<R_VAR>> for
# R-side interpolation; bash variables stay literal.
#
# Module-load ordering is CRITICAL: load tool modules FIRST, then activate
# the conda env LAST. Otherwise module-loaded paths override conda's bin
# and you get the broken stringtie/2.2.1 module (bgzf assertion crash).
# Discovered during validation 2026-05-20.

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# =============================================================================
# Locked constants — single source of truth (CLAUDE.md rule #3)
# =============================================================================

PROTEOG_MODULES <- list(
  fastp        = "fastp/0.23.4",
  bowtie2      = "bowtie2/2.5.2",
  star         = "star/2.7.11a",
  samtools     = "samtools/1.19.2",
  gffread      = "gffread/0.12.7",
  transdecoder = "transdecoder/5.7.1",
  diamond      = "diamond/2.1.7",
  gffcompare   = "gffcompare/0.12.6"
)

PROTEOG_CONDA_ENV  <- "/quobyte/proteomics-grp/de-limp/envs/proteog_helpers"
PROTEOG_REWRITER   <- "/quobyte/proteomics-grp/de-limp/references/scripts/rewrite_transdecoder_headers.py"

# Default SLURM submission targets. The orchestrator overrides per build based
# on `sacctmgr show user $USER` output.
PROTEOG_DEFAULT_ACCOUNT   <- "genome-center-grp"
PROTEOG_DEFAULT_PARTITION <- "high"
PROTEOG_DEFAULT_QOS       <- "genome-center-grp-high-qos"

# =============================================================================
# Read-length adaptive STAR thresholds
# =============================================================================

#' Select STAR alignment parameters based on read length
#'
#' Locked per spec §5.3. Defaults assume DNA Tech Core's PE150 deliverable;
#' tiers below 130 bp progressively relax the score/match-length filters.
#' Reads <60 bp are refused outright (not suitable for proteogenomics).
#'
#' @param read_length integer or numeric — median post-fastp read length
#' @return list with $tier, $flags (character vector), $qc_gate_unique_pct,
#'   and optionally $warning or $error.
select_star_params <- function(read_length) {
  if (!is.numeric(read_length) || !is.finite(read_length) || read_length < 0) {
    return(list(tier = "refuse",
                error = sprintf("Invalid read length: %s", as.character(read_length))))
  }

  if (read_length >= 130) {
    list(
      tier = "default",
      flags = c(
        "--outFilterMismatchNoverLmax 0.04",
        "--outFilterScoreMinOverLread 0.66",
        "--outFilterMatchNminOverLread 0.66",
        "--outFilterMultimapNmax 20"
      ),
      qc_gate_unique_pct = 60
    )
  } else if (read_length >= 100) {
    list(
      tier = "mildly_relaxed",
      flags = c(
        "--outFilterMismatchNoverLmax 0.06",
        "--outFilterScoreMinOverLread 0.50",
        "--outFilterMatchNminOverLread 0.50",
        "--outFilterMultimapNmax 20"
      ),
      qc_gate_unique_pct = 45
    )
  } else if (read_length >= 60) {
    list(
      tier = "significantly_relaxed",
      flags = c(
        "--outFilterMismatchNoverLmax 0.10",
        "--outFilterScoreMinOverLread 0.30",
        "--outFilterMatchNminOverLread 0.30",
        "--outFilterMultimapNmax 20"
      ),
      qc_gate_unique_pct = 25,
      warning = "Read length below 100bp produces noisier transcript assembly. Consider commissioning a higher-quality RNA-seq run."
    )
  } else {
    list(
      tier = "refuse",
      error = sprintf(
        "Read length %.0fbp is too short (<60bp threshold) for proteogenomics novel-ORF discovery. Please consult the Proteomics Core about appropriate library design.",
        read_length
      )
    )
  }
}

#' Estimate median read length from a gzipped FASTQ
#'
#' Streams first N reads (default 100), parses the seq line of each record
#' (line 2 of every 4-line block), returns the median nchar. Used after
#' fastp to drive STAR tier selection.
#'
#' @param fastq_gz character — path to .fastq.gz
#' @param n_reads integer — how many reads to sample
#' @return numeric — median read length, or NA if file unreadable
detect_read_length <- function(fastq_gz, n_reads = 100L) {
  if (!file.exists(fastq_gz)) {
    warning("detect_read_length(): file not found: ", fastq_gz)
    return(NA_real_)
  }
  con <- gzfile(fastq_gz, open = "rt")
  on.exit(close(con), add = TRUE)
  lengths <- integer(n_reads)
  i <- 0L
  repeat {
    block <- readLines(con, n = 4L, warn = FALSE)
    if (length(block) < 2L) break
    i <- i + 1L
    lengths[i] <- nchar(block[2])
    if (i >= n_reads) break
  }
  if (i == 0L) return(NA_real_)
  median(lengths[seq_len(i)])
}

# =============================================================================
# Sbatch header — used by every generator
# =============================================================================

#' Construct the standard sbatch header
#'
#' @param job_name character
#' @param time character — e.g. "1:00:00"
#' @param mem character — e.g. "8G"
#' @param cpus integer
#' @param array character — e.g. "1-12" or NULL
#' @param slurm_account character
#' @param slurm_partition character
#' @param slurm_qos character — default depends on partition
#' @param out_dir character — directory for %j logs
#' @return character (multi-line)
slurm_header <- function(job_name,
                         time,
                         mem,
                         cpus,
                         array = NULL,
                         slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                         slurm_partition = PROTEOG_DEFAULT_PARTITION,
                         slurm_qos       = NULL,
                         out_dir) {
  qos <- slurm_qos %||% sprintf("%s-%s-qos", slurm_account, slurm_partition)
  array_line <- if (!is.null(array)) sprintf("#SBATCH --array=%s\n", array) else ""
  paste0(
    "#!/bin/bash -l\n",
    sprintf("#SBATCH --job-name=%s\n",  job_name),
    sprintf("#SBATCH --account=%s\n",   slurm_account),
    sprintf("#SBATCH --partition=%s\n", slurm_partition),
    sprintf("#SBATCH --qos=%s\n",       qos),
    array_line,
    sprintf("#SBATCH --time=%s\n",      time),
    sprintf("#SBATCH --mem=%s\n",       mem),
    sprintf("#SBATCH --cpus-per-task=%d\n", as.integer(cpus)),
    sprintf("#SBATCH -o %s/%s_%%j.out\n", out_dir, job_name),
    sprintf("#SBATCH -e %s/%s_%%j.err\n", out_dir, job_name)
  )
}

# =============================================================================
# Per-stage sbatch generators
# =============================================================================

#' fastp — adapter trim + quality filter, parallel array
generate_fastp_sbatch <- function(project_dir,
                                  sample_names,
                                  slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                  slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  n <- length(sample_names)
  if (n == 0) stop("generate_fastp_sbatch(): no samples")
  sample_list <- paste(shQuote(sample_names), collapse = " ")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_fastp",
      time = "1:00:00", mem = "8G", cpus = 8,
      array = sprintf("1-%d", n),
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    sprintf("module load %s\n\n", PROTEOG_MODULES$fastp),
    sprintf("PROJECT_DIR=%s\n", shQuote(project_dir)),
    sprintf("SAMPLES=(%s)\n",   sample_list),
    "SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}\n\n",
    "mkdir -p \"$PROJECT_DIR/fastp_out\"\n\n",
    "fastp \\\n",
    "  --in1  \"$PROJECT_DIR/rnaseq/${SAMPLE}_R1.fastq.gz\" \\\n",
    "  --in2  \"$PROJECT_DIR/rnaseq/${SAMPLE}_R2.fastq.gz\" \\\n",
    "  --out1 \"$PROJECT_DIR/fastp_out/${SAMPLE}_R1.fastq.gz\" \\\n",
    "  --out2 \"$PROJECT_DIR/fastp_out/${SAMPLE}_R2.fastq.gz\" \\\n",
    "  --json \"$PROJECT_DIR/fastp_out/${SAMPLE}.json\" \\\n",
    "  --html \"$PROJECT_DIR/fastp_out/${SAMPLE}.html\" \\\n",
    "  --thread 8 \\\n",
    "  --detect_adapter_for_pe \\\n",
    "  --length_required 36\n"
  )
}

#' rRNA pre-filter — bowtie2 array, keeps reads NOT mapped to rRNA
generate_rrna_sbatch <- function(project_dir,
                                 sample_names,
                                 rrna_index_prefix,
                                 slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                 slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  n <- length(sample_names)
  if (n == 0) stop("generate_rrna_sbatch(): no samples")
  if (!nzchar(rrna_index_prefix)) {
    stop("generate_rrna_sbatch(): rrna_index_prefix is empty")
  }
  sample_list <- paste(shQuote(sample_names), collapse = " ")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_rrna",
      time = "1:00:00", mem = "8G", cpus = 8,
      array = sprintf("1-%d", n),
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    sprintf("module load %s\n\n", PROTEOG_MODULES$bowtie2),
    sprintf("PROJECT_DIR=%s\n",    shQuote(project_dir)),
    sprintf("RRNA_INDEX=%s\n",      shQuote(rrna_index_prefix)),
    sprintf("SAMPLES=(%s)\n",       sample_list),
    "SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}\n\n",
    "mkdir -p \"$PROJECT_DIR/rrna_filt\"\n\n",
    "bowtie2 -x \"$RRNA_INDEX\" \\\n",
    "  -1 \"$PROJECT_DIR/fastp_out/${SAMPLE}_R1.fastq.gz\" \\\n",
    "  -2 \"$PROJECT_DIR/fastp_out/${SAMPLE}_R2.fastq.gz\" \\\n",
    "  --very-sensitive-local \\\n",
    "  --un-conc-gz \"$PROJECT_DIR/rrna_filt/${SAMPLE}_norrna_R%.fq.gz\" \\\n",
    "  -S /dev/null \\\n",
    "  -p 8 \\\n",
    "  --no-unal 2> \"$PROJECT_DIR/rrna_filt/${SAMPLE}_rrna_filter.log\"\n"
  )
}

#' STAR alignment array — uses tier-specific flags from select_star_params()
generate_star_sbatch <- function(project_dir,
                                 sample_names,
                                 star_index,
                                 tier_params,
                                 slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                 slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  n <- length(sample_names)
  if (n == 0) stop("generate_star_sbatch(): no samples")
  if (isTRUE(tier_params$tier == "refuse")) {
    stop("generate_star_sbatch(): refuse tier — pipeline must not start: ",
         tier_params$error %||% "unknown reason")
  }
  if (is.null(tier_params$flags) || length(tier_params$flags) == 0) {
    stop("generate_star_sbatch(): tier_params has no $flags")
  }
  sample_list <- paste(shQuote(sample_names), collapse = " ")
  adaptive_flags <- paste(tier_params$flags, collapse = " ")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_star",
      time = "4:00:00", mem = "48G", cpus = 16,
      array = sprintf("1-%d", n),
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    sprintf("module load %s\n\n", PROTEOG_MODULES$star),
    sprintf("PROJECT_DIR=%s\n", shQuote(project_dir)),
    sprintf("STAR_INDEX=%s\n",   shQuote(star_index)),
    sprintf("SAMPLES=(%s)\n",    sample_list),
    "SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}\n\n",
    "mkdir -p \"$PROJECT_DIR/star_out\"\n",
    "STAR_TMP=\"$PROJECT_DIR/star_out/_STARtmp_${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_ID}\"\n",
    "rm -rf \"$STAR_TMP\"\n\n",
    "STAR \\\n",
    "  --runMode alignReads \\\n",
    "  --runThreadN 16 \\\n",
    "  --genomeDir \"$STAR_INDEX\" \\\n",
    "  --readFilesIn \"$PROJECT_DIR/rrna_filt/${SAMPLE}_norrna_R1.fq.gz\" \\\n",
    "                \"$PROJECT_DIR/rrna_filt/${SAMPLE}_norrna_R2.fq.gz\" \\\n",
    "  --readFilesCommand zcat \\\n",
    "  --outSAMtype BAM SortedByCoordinate \\\n",
    "  --outSAMstrandField intronMotif \\\n",
    "  --outSAMattributes Standard XS NH \\\n",
    "  --outFileNamePrefix \"$PROJECT_DIR/star_out/${SAMPLE}_\" \\\n",
    "  --outTmpDir \"$STAR_TMP\" \\\n",
    sprintf("  %s\n\n", adaptive_flags),
    "rm -rf \"$STAR_TMP\"\n"
  )
}

#' QC gate — halts dependency chain if any sample's unique% below tier threshold.
#' Exits non-zero on failure so downstream `--dependency=afterok` jobs never run.
generate_qc_gate_sbatch <- function(project_dir,
                                    sample_names,
                                    qc_gate_unique_pct,
                                    slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                    slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  if (length(sample_names) == 0) stop("generate_qc_gate_sbatch(): no samples")
  sample_list <- paste(shQuote(sample_names), collapse = " ")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_qc_gate",
      time = "10:00", mem = "2G", cpus = 1,
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n\n",
    sprintf("PROJECT_DIR=%s\n", shQuote(project_dir)),
    sprintf("SAMPLES=(%s)\n",   sample_list),
    sprintf("THRESHOLD=%d\n",   as.integer(qc_gate_unique_pct)),
    "\nfail=0\n",
    "results=()\n",
    "for SAMPLE in \"${SAMPLES[@]}\"; do\n",
    "  log=\"$PROJECT_DIR/star_out/${SAMPLE}_Log.final.out\"\n",
    "  if [ ! -f \"$log\" ]; then\n",
    "    echo \"MISSING: $log\" >&2\n",
    "    fail=1\n",
    "    continue\n",
    "  fi\n",
    "  pct=$(grep 'Uniquely mapped reads %' \"$log\" | awk -F'|' '{print $2}' | tr -d ' %\\t')\n",
    "  pct_int=${pct%.*}\n",
    "  if [ -z \"$pct_int\" ]; then pct_int=0; fi\n",
    "  results+=(\"${SAMPLE}: ${pct}% (threshold ${THRESHOLD}%)\")\n",
    "  if [ \"$pct_int\" -lt \"$THRESHOLD\" ]; then fail=1; fi\n",
    "done\n\n",
    "printf 'QC gate results (uniquely mapped %% per sample):\\n' > \"$PROJECT_DIR/qc_gate_result.txt\"\n",
    "printf '  %s\\n' \"${results[@]}\" >> \"$PROJECT_DIR/qc_gate_result.txt\"\n",
    "cat \"$PROJECT_DIR/qc_gate_result.txt\"\n\n",
    "if [ \"$fail\" -ne 0 ]; then\n",
    "  echo >&2\n",
    "  echo 'QC GATE FAILED — alignment quality below threshold.' >&2\n",
    "  echo 'Likely causes (in order of frequency):' >&2\n",
    "  echo '  1. Wrong reference genome selected' >&2\n",
    "  echo '  2. Heavy non-rRNA contamination (mitochondrial, host cell line, bacterial)' >&2\n",
    "  echo '  3. Library type unsuited to this pipeline (Ribo-Seq, CLIP-Seq, 3-prime Tag-Seq)' >&2\n",
    "  echo '  4. Severe sample degradation' >&2\n",
    "  exit 1\n",
    "fi\n",
    "echo 'QC gate passed.'\n"
  )
}

#' Per-sample StringTie assembly array
#'
#' Uses the conda-installed stringtie 3.0.3 — the native stringtie/2.2.1 module
#' has a fatal bgzf assertion crash. Module loads BEFORE conda activation so
#' the conda bin/ wins PATH precedence.
generate_stringtie_sbatch <- function(project_dir,
                                      sample_names,
                                      ref_gtf,
                                      strand_flag = "",
                                      slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                      slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  n <- length(sample_names)
  if (n == 0) stop("generate_stringtie_sbatch(): no samples")
  if (!nzchar(ref_gtf)) stop("generate_stringtie_sbatch(): ref_gtf is empty")
  if (!strand_flag %in% c("", "--rf", "--fr")) {
    stop("generate_stringtie_sbatch(): strand_flag must be '', '--rf', or '--fr'; got: ", strand_flag)
  }
  sample_list <- paste(shQuote(sample_names), collapse = " ")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_stringtie",
      time = "2:00:00", mem = "16G", cpus = 8,
      array = sprintf("1-%d", n),
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n\n",
    "# Load conda LAST so its bin (with stringtie 3.0.3) wins PATH precedence.\n",
    "# (The native stringtie/2.2.1 module has a fatal bgzf assertion crash.)\n",
    "module load conda\n",
    sprintf("source activate %s\n\n", PROTEOG_CONDA_ENV),
    sprintf("PROJECT_DIR=%s\n", shQuote(project_dir)),
    sprintf("REF_GTF=%s\n",     shQuote(ref_gtf)),
    sprintf("SAMPLES=(%s)\n",   sample_list),
    "SAMPLE=${SAMPLES[$((SLURM_ARRAY_TASK_ID - 1))]}\n\n",
    "mkdir -p \"$PROJECT_DIR/stringtie_out\"\n\n",
    "stringtie \"$PROJECT_DIR/star_out/${SAMPLE}_Aligned.sortedByCoord.out.bam\" \\\n",
    "  -G \"$REF_GTF\" \\\n",
    "  -o \"$PROJECT_DIR/stringtie_out/${SAMPLE}.gtf\" \\\n",
    "  -p 8 \\\n",
    sprintf("  -l \"STRG_${SAMPLE}\" %s\n", strand_flag)
  )
}

#' StringTie --merge — combine per-sample GTFs into a unified transcript model
generate_merge_sbatch <- function(project_dir,
                                  ref_gtf,
                                  sample_names,
                                  slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                  slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  if (length(sample_names) == 0) stop("generate_merge_sbatch(): no samples")
  if (!nzchar(ref_gtf)) stop("generate_merge_sbatch(): ref_gtf is empty")
  sample_list <- paste(shQuote(sample_names), collapse = " ")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_merge",
      time = "1:00:00", mem = "32G", cpus = 8,
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    "module load conda\n",
    sprintf("source activate %s\n\n", PROTEOG_CONDA_ENV),
    sprintf("PROJECT_DIR=%s\n", shQuote(project_dir)),
    sprintf("REF_GTF=%s\n",     shQuote(ref_gtf)),
    sprintf("SAMPLES=(%s)\n",   sample_list),
    "\nGTF_LIST=\"$PROJECT_DIR/stringtie_out/gtf_list.txt\"\n",
    ": > \"$GTF_LIST\"\n",
    "for SAMPLE in \"${SAMPLES[@]}\"; do\n",
    "  echo \"$PROJECT_DIR/stringtie_out/${SAMPLE}.gtf\" >> \"$GTF_LIST\"\n",
    "done\n\n",
    "stringtie --merge \\\n",
    "  -G \"$REF_GTF\" \\\n",
    "  -o \"$PROJECT_DIR/stringtie_out/merged.gtf\" \\\n",
    "  \"$GTF_LIST\"\n"
  )
}

#' gffcompare — class codes (=/c/j/e/i/u/o/x/...) for NOVEL_ISOFORM detection
generate_gffcompare_sbatch <- function(project_dir,
                                       ref_gtf,
                                       slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                       slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  if (!nzchar(ref_gtf)) stop("generate_gffcompare_sbatch(): ref_gtf is empty")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_gffcompare",
      time = "30:00", mem = "8G", cpus = 1,
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    sprintf("module load %s\n\n", PROTEOG_MODULES$gffcompare),
    sprintf("PROJECT_DIR=%s\n", shQuote(project_dir)),
    sprintf("REF_GTF=%s\n",     shQuote(ref_gtf)),
    "\ncd \"$PROJECT_DIR/stringtie_out\"\n",
    "gffcompare \\\n",
    "  -r \"$REF_GTF\" \\\n",
    "  -o gffcmp \\\n",
    "  merged.gtf\n",
    "\necho 'gffcompare done; class code summary:'\n",
    "awk 'NR>1 {print $3}' gffcmp.merged.gtf.tmap | sort | uniq -c | sort -rn\n"
  )
}

#' gffread — extract transcript FASTA from genome + merged.gtf
generate_gffread_sbatch <- function(project_dir,
                                    genome_fasta,
                                    slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                    slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  if (!nzchar(genome_fasta)) stop("generate_gffread_sbatch(): genome_fasta is empty")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_gffread",
      time = "30:00", mem = "16G", cpus = 4,
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    sprintf("module load %s\n\n", PROTEOG_MODULES$gffread),
    sprintf("PROJECT_DIR=%s\n", shQuote(project_dir)),
    sprintf("GENOME=%s\n",       shQuote(genome_fasta)),
    "\ngffread -w \"$PROJECT_DIR/stringtie_out/merged_transcripts.fa\" \\\n",
    "  -g \"$GENOME\" \\\n",
    "  \"$PROJECT_DIR/stringtie_out/merged.gtf\"\n",
    "\necho 'Transcripts extracted:'\n",
    "grep -c '^>' \"$PROJECT_DIR/stringtie_out/merged_transcripts.fa\"\n"
  )
}

#' TransDecoder — LongOrfs + Predict (optionally DIAMOND-supported)
generate_transdecoder_sbatch <- function(project_dir,
                                         diamond_db = NULL,
                                         min_orf_len = 100L,
                                         slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                         slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  logs_dir <- file.path(project_dir, "logs")
  use_diamond <- !is.null(diamond_db) && nzchar(diamond_db)

  predict_block <- if (use_diamond) {
    paste0(
      "if [ -f \"$DIAMOND_DB\" ]; then\n",
      "  diamond blastp \\\n",
      "    --query \"merged_transcripts.fa.transdecoder_dir/longest_orfs.pep\" \\\n",
      "    --db \"$DIAMOND_DB\" \\\n",
      "    --max-target-seqs 1 --outfmt 6 --evalue 1e-5 --threads 16 \\\n",
      "    --out blastp.outfmt6\n",
      "  TransDecoder.Predict \\\n",
      "    -t \"$PROJECT_DIR/stringtie_out/merged_transcripts.fa\" \\\n",
      "    --retain_blastp_hits blastp.outfmt6 \\\n",
      "    --single_best_only\n",
      "else\n",
      "  TransDecoder.Predict \\\n",
      "    -t \"$PROJECT_DIR/stringtie_out/merged_transcripts.fa\" \\\n",
      "    --single_best_only\n",
      "fi\n"
    )
  } else {
    paste0(
      "TransDecoder.Predict \\\n",
      "  -t \"$PROJECT_DIR/stringtie_out/merged_transcripts.fa\" \\\n",
      "  --single_best_only\n"
    )
  }

  paste0(
    slurm_header(
      job_name = "proteog_transdecoder",
      time = "4:00:00", mem = "32G", cpus = 16,
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    sprintf("module load %s\n", PROTEOG_MODULES$transdecoder),
    if (use_diamond) sprintf("module load %s\n", PROTEOG_MODULES$diamond) else "",
    "\n",
    sprintf("PROJECT_DIR=%s\n", shQuote(project_dir)),
    if (use_diamond) sprintf("DIAMOND_DB=%s\n", shQuote(diamond_db)) else "",
    sprintf("MIN_ORF_LEN=%d\n", as.integer(min_orf_len)),
    "\nmkdir -p \"$PROJECT_DIR/transdecoder_out\"\n",
    "cd \"$PROJECT_DIR/transdecoder_out\"\n\n",
    "TransDecoder.LongOrfs \\\n",
    "  -t \"$PROJECT_DIR/stringtie_out/merged_transcripts.fa\" \\\n",
    "  -m \"$MIN_ORF_LEN\"\n\n",
    predict_block,
    "\necho 'TransDecoder done; ORF count:'\n",
    "grep -c '^>' \"merged_transcripts.fa.transdecoder.pep\"\n"
  )
}

#' Header rewriter — Python script in proteog_helpers conda env
#'
#' Single canonical signature (v1.0, locked Phase B.5):
#'   --transdecoder --merged-gtf --gffcompare-tmap --project-tag --output
#'
#' Produces `sp|<protein_id>|<sym>_<TAG>` headers with 7 metadata key=value
#' fields and source classes REF / NOVEL_GENE / NOVEL_ISOFORM / UNPARSED.
#' The rewriter itself exits non-zero on any UNPARSED entry; this sbatch
#' relies on that exit code (no shell-level recount needed).
generate_rewrite_sbatch <- function(project_dir,
                                    project_tag,
                                    rewriter_path = PROTEOG_REWRITER,
                                    conda_env_path = PROTEOG_CONDA_ENV,
                                    slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                    slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  if (!nzchar(project_tag)) stop("generate_rewrite_sbatch(): project_tag is empty")
  if (!grepl("^[A-Za-z0-9_-]+$", project_tag)) {
    stop("generate_rewrite_sbatch(): project_tag must match [A-Za-z0-9_-]+; got: ", project_tag)
  }
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_rewrite",
      time = "30:00", mem = "8G", cpus = 2,
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    "module load conda\n",
    sprintf("source activate %s\n\n", conda_env_path),
    sprintf("PROJECT_DIR=%s\n",  shQuote(project_dir)),
    sprintf("PROJECT_TAG=%s\n",  shQuote(project_tag)),
    sprintf("REWRITER=%s\n",     shQuote(rewriter_path)),
    "\npython3 \"$REWRITER\" \\\n",
    "  --transdecoder \"$PROJECT_DIR/transdecoder_out/merged_transcripts.fa.transdecoder.pep\" \\\n",
    "  --merged-gtf \"$PROJECT_DIR/stringtie_out/merged.gtf\" \\\n",
    "  --gffcompare-tmap \"$PROJECT_DIR/stringtie_out/gffcmp.merged.gtf.tmap\" \\\n",
    "  --project-tag \"$PROJECT_TAG\" \\\n",
    "  --output \"$PROJECT_DIR/predicted_orfs.fasta\"\n",
    "\necho 'Header rewrite done; ORF count:'\n",
    "grep -c '^>' \"$PROJECT_DIR/predicted_orfs.fasta\"\n"
  )
}

#' Final assembly: concatenate predicted_orfs.fasta + optional UniProt FASTA
#' and (if available) dedupe with seqkit. Output lands in the proteogenomics
#' databases dir as `<project>_proteogenomics_<YYYY_MM>.fasta` so DE-LIMP can
#' pick it up from the FASTA library scanner.
#'
#' UniProt FASTA is optional. If empty/missing, the output is just the
#' deduplicated predicted ORFs.
generate_assemble_sbatch <- function(project_dir,
                                     project_name,
                                     uniprot_fasta = "",
                                     databases_root = "/quobyte/proteomics-grp/de-limp/databases/proteogenomics",
                                     slurm_account   = PROTEOG_DEFAULT_ACCOUNT,
                                     slurm_partition = PROTEOG_DEFAULT_PARTITION) {
  if (!nzchar(project_name)) stop("generate_assemble_sbatch(): project_name is empty")
  logs_dir <- file.path(project_dir, "logs")

  paste0(
    slurm_header(
      job_name = "proteog_assemble",
      time = "30:00", mem = "4G", cpus = 2,
      slurm_account = slurm_account, slurm_partition = slurm_partition,
      out_dir = logs_dir
    ),
    "\nset -euo pipefail\n",
    sprintf("PROJECT_DIR=%s\n",   shQuote(project_dir)),
    sprintf("PROJECT_NAME=%s\n",  shQuote(project_name)),
    sprintf("UNIPROT=%s\n",       shQuote(uniprot_fasta %||% "")),
    sprintf("OUT_DIR=%s\n",       shQuote(databases_root)),
    "DATE_TAG=$(date +%Y_%m)\n",
    "OUT=\"$OUT_DIR/${PROJECT_NAME}_proteogenomics_${DATE_TAG}.fasta\"\n",
    "\nmkdir -p \"$OUT_DIR\"\n",
    "INPUTS=(\"$PROJECT_DIR/predicted_orfs.fasta\")\n",
    "if [ -n \"$UNIPROT\" ] && [ -f \"$UNIPROT\" ]; then INPUTS+=(\"$UNIPROT\"); fi\n",
    "echo \"Concatenating: ${INPUTS[@]} -> $OUT\"\n",
    "cat \"${INPUTS[@]}\" > \"$OUT\"\n",
    "\n# Optional dedupe — seqkit is not always available on Hive.\n",
    "if command -v seqkit >/dev/null 2>&1; then\n",
    "  seqkit rmdup -s -o \"${OUT}.dedup\" \"$OUT\" && mv \"${OUT}.dedup\" \"$OUT\"\n",
    "  echo \"Deduplicated with seqkit rmdup -s\"\n",
    "else\n",
    "  echo \"seqkit not on PATH — skipping dedupe (duplicates are harmless for DIA-NN).\"\n",
    "fi\n",
    "\necho \"Final FASTA: $OUT\"\n",
    "echo \"Sequence count: $(grep -c '^>' \"$OUT\")\"\n"
  )
}
