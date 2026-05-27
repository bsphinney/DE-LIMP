# server_proteog_builder.R — orchestrator for the proteogenomics RNA-seq pipeline.
#
# Three public functions:
#   submit_proteogenomics_build() — kick off the full SLURM dep chain
#   poll_proteog_build_status()  — read status.json + refresh per-stage via sacct
#   cancel_proteog_build()       — scancel the whole chain
#
# This is the FIRST file in the proteogenomics stack with side effects (writes
# sbatch scripts to disk, submits to SLURM, updates status.json). All pure
# script-generation lives in helpers_rnaseq.R. All pure log-parsing lives in
# helpers_proteog_qc.R. This file only orchestrates.
#
# Local vs remote execution:
#   When R runs on Hive (via apptainer de-limp.sif), sbatch/squeue/sacct are
#   on PATH (via SLURM proxy or login-shell). This file assumes local execution
#   via system2(). For Docker-from-laptop, a future iteration can layer an
#   ssh_exec wrapper. Phase C v1 is single-host.
#
# Self-describing build manifest (CLAUDE.md rule #1):
#   Every build object returned from submit_*() carries $pipeline_id,
#   $methods_paragraph, and $stages so downstream consumers (Claude export,
#   methods README, status panel) read these rather than hardcoding "what
#   the pipeline did."

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# Robust scalar-string coercion. JSON `{}` and `null` get parsed by jsonlite
# as list() and NULL respectively; nzchar() on those returns logical(0)/NA
# and crashes downstream `if` branches. Coerce all degenerate forms to "".
.empty_or_str <- function(v) {
  if (is.null(v)) return("")
  if (length(v) == 0) return("")
  if (length(v) == 1 && is.na(v[[1]])) return("")
  chr <- tryCatch(as.character(v)[1], error = function(e) "")
  if (length(chr) == 0 || is.na(chr)) "" else chr
}

PROTEOG_PIPELINE_ID    <- "proteogenomics_v1.1"
PROTEOG_RNASEQ_ROOT    <- "/quobyte/proteomics-grp/de-limp/rnaseq"
PROTEOG_DATABASES_ROOT <- "/quobyte/proteomics-grp/de-limp/databases/proteogenomics"

# Stage order is the contract — used by status.json schema, poller, and cancel.
PROTEOG_STAGE_ORDER <- c(
  "fastp", "rrna_filter", "star", "qc_gate", "stringtie",
  "merge", "gffcompare", "gffread", "transdecoder", "rewrite", "assemble"
)

# =============================================================================
# Input validation + project setup
# =============================================================================

.validate_build_inputs <- function(project_name, rnaseq_dir, sample_names,
                                   reference_key, library_type, strand_flag,
                                   ref_registry) {
  if (!nzchar(project_name) || !grepl("^[A-Za-z0-9_.-]+$", project_name)) {
    stop("submit_proteogenomics_build(): project_name must be non-empty and match [A-Za-z0-9_.-]+; got: ",
         project_name)
  }
  if (length(sample_names) == 0) {
    stop("submit_proteogenomics_build(): no sample_names provided")
  }
  if (any(!grepl("^[A-Za-z0-9_.-]+$", sample_names))) {
    bad <- sample_names[!grepl("^[A-Za-z0-9_.-]+$", sample_names)]
    stop("submit_proteogenomics_build(): invalid sample name(s): ",
         paste(bad, collapse = ", "))
  }
  if (!dir.exists(rnaseq_dir)) {
    stop("submit_proteogenomics_build(): rnaseq_dir does not exist: ", rnaseq_dir)
  }
  # Each sample must have R1 + R2 in rnaseq_dir
  for (s in sample_names) {
    for (rd in c("R1", "R2")) {
      f <- file.path(rnaseq_dir, sprintf("%s_%s.fastq.gz", s, rd))
      if (!file.exists(f)) {
        stop(sprintf(
          "submit_proteogenomics_build(): missing input FASTQ for sample %s: %s",
          s, f
        ))
      }
    }
  }
  if (!reference_key %in% names(ref_registry)) {
    stop("submit_proteogenomics_build(): reference_key not in registry: ",
         reference_key, " (known: ",
         paste(names(ref_registry), collapse = ", "), ")")
  }
  if (!library_type %in% c("polyA", "ribo_depleted", "stranded", "unstranded")) {
    stop("submit_proteogenomics_build(): library_type must be one of ",
         "polyA/ribo_depleted/stranded/unstranded; got: ", library_type)
  }
  if (!strand_flag %in% c("", "--rf", "--fr")) {
    stop("submit_proteogenomics_build(): strand_flag must be '', '--rf', or '--fr'; got: ",
         strand_flag)
  }
  invisible(TRUE)
}

#' Stage raw FASTQs into the project's expected layout
#'
#' The Phase B sbatch generators expect FASTQs at <project_dir>/rnaseq/<sample>_R{1,2}.fastq.gz.
#' If rnaseq_dir is different (e.g., user pointed at sra_data/), we symlink the
#' files into place. Idempotent.
.stage_rnaseq_inputs <- function(project_dir, rnaseq_dir, sample_names,
                                  ssh_config = NULL) {
  target <- file.path(project_dir, "rnaseq")
  .fs_mkdir(target, ssh_config = ssh_config)
  for (s in sample_names) {
    for (rd in c("R1", "R2")) {
      fname <- sprintf("%s_%s.fastq.gz", s, rd)
      src <- file.path(rnaseq_dir, fname)
      dst <- file.path(target, fname)
      .fs_symlink(src, dst, ssh_config = ssh_config)
    }
  }
  invisible(target)
}

# =============================================================================
# Filesystem helpers — local or via SSH depending on ssh_config
# =============================================================================
# These dispatch based on ssh_config: NULL → local (DE-LIMP-on-Hive case),
# non-NULL → remote via ssh_exec / scp_upload (DE-LIMP-on-Mac case).

#' Make a directory (with -p / recursive semantics) on local or remote host
.fs_mkdir <- function(path, ssh_config = NULL) {
  if (is.null(ssh_config)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  } else {
    res <- ssh_exec(ssh_config, sprintf("mkdir -p %s", shQuote(path)),
                    login_shell = FALSE, timeout = 10)
    if (!identical(res$status, 0L)) {
      stop("ssh mkdir failed for ", path, ": ",
           paste(res$stderr %||% character(), collapse = "; "))
    }
  }
  invisible(path)
}

#' Write a text file (e.g. sbatch script, status.json) — local or remote
.fs_write_text <- function(content, path, ssh_config = NULL, executable = FALSE) {
  if (is.null(ssh_config)) {
    writeLines(content, path)
    if (executable) Sys.chmod(path, "755")
  } else {
    tmp <- tempfile(pattern = "proteog_write_")
    on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
    writeLines(content, tmp)
    scp_res <- tryCatch(scp_upload(ssh_config, tmp, path, timeout = 60),
                        error = function(e) list(status = -1,
                                                  stderr = conditionMessage(e)))
    if (!identical(scp_res$status, 0L)) {
      stop("scp_upload failed for ", path, ": ",
           paste(scp_res$stderr %||% character(), collapse = "; "))
    }
    if (executable) {
      ssh_exec(ssh_config, sprintf("chmod 755 %s", shQuote(path)),
               login_shell = FALSE, timeout = 5)
    }
  }
  invisible(path)
}

#' Symlink a remote file into a remote project directory
.fs_symlink <- function(src, dst, ssh_config = NULL) {
  if (is.null(ssh_config)) {
    if (!file.exists(dst)) {
      file.symlink(normalizePath(src, mustWork = TRUE), dst)
    }
  } else {
    cmd <- sprintf("ln -sf %s %s", shQuote(src), shQuote(dst))
    res <- ssh_exec(ssh_config, cmd, login_shell = FALSE, timeout = 10)
    if (!identical(res$status, 0L)) {
      stop("ssh ln failed: ", paste(res$stderr %||% character(), collapse = "; "))
    }
  }
  invisible(dst)
}

#' Read a remote text file (e.g. status.json, log) — returns single string
.fs_read_text <- function(path, ssh_config = NULL) {
  if (is.null(ssh_config)) {
    if (!file.exists(path)) return(NULL)
    paste(readLines(path, warn = FALSE), collapse = "\n")
  } else {
    res <- ssh_exec(ssh_config, sprintf("cat %s", shQuote(path)),
                    login_shell = FALSE, timeout = 15)
    if (!identical(res$status, 0L) || length(res$stdout) == 0) return(NULL)
    paste(res$stdout, collapse = "\n")
  }
}

#' Detect median read length of the first n reads of a gzipped FASTQ
#'
#' Replaces the in-process detect_read_length() helper when running over SSH.
#' Both branches return the median nchar of the first n_reads seq lines, or NA.
.fs_detect_read_length <- function(fastq_gz, n_reads = 100L, ssh_config = NULL) {
  if (is.null(ssh_config)) {
    return(detect_read_length(fastq_gz, n_reads = n_reads))
  }
  # zcat | awk 'NR%4==2 {print length($1)}' | head -n N | sort -n | awk-median
  cmd <- sprintf(
    "zcat %s 2>/dev/null | awk 'NR%%4==2 {print length($1)}' | head -n %d",
    shQuote(fastq_gz), as.integer(n_reads)
  )
  res <- ssh_exec(ssh_config, cmd, login_shell = FALSE, timeout = 30)
  if (!identical(res$status, 0L) || length(res$stdout) == 0) return(NA_real_)
  lens <- suppressWarnings(as.integer(trimws(res$stdout)))
  lens <- lens[!is.na(lens) & lens > 0]
  if (length(lens) == 0) return(NA_real_)
  median(lens)
}

# =============================================================================
# Sbatch dispatch — local shell submission via system2() or SSH-relayed
# =============================================================================

#' Run sbatch on a script file; return parsed job_id or stop on failure.
#'
#' Dispatches based on ssh_config: NULL → direct local system2 call (Hive
#' apptainer case, where R inherits PATH from login shell); non-NULL →
#' remote sbatch via ssh_exec (Mac+SSH-to-Hive case).
.sbatch_submit <- function(script_path, dep_jid = NULL, ssh_config = NULL) {
  dep_arg <- if (!is.null(dep_jid) && nzchar(dep_jid)) {
    sprintf("--dependency=afterok:%s ", dep_jid)
  } else ""

  if (is.null(ssh_config)) {
    # Local execution
    args <- character()
    if (nzchar(dep_arg)) args <- c(args, trimws(dep_arg))
    args <- c(args, script_path)
    out <- tryCatch(
      suppressWarnings(system2("sbatch", args = args,
                                stdout = TRUE, stderr = TRUE)),
      error = function(e) stop("sbatch failed: ", conditionMessage(e))
    )
  } else {
    # Remote execution via SSH. Use login shell so SLURM tools are on PATH.
    cmd <- sprintf("sbatch %s%s", dep_arg, shQuote(script_path))
    res <- tryCatch(
      ssh_exec(ssh_config, cmd, login_shell = TRUE, timeout = 30),
      error = function(e) stop("ssh sbatch failed: ", conditionMessage(e))
    )
    if (!identical(res$status, 0L)) {
      stop("ssh sbatch returned non-zero. Output: ",
           paste(c(res$stdout, res$stderr %||% character()), collapse = "\n"))
    }
    out <- res$stdout
  }

  if (!is.character(out) || !any(grepl("Submitted batch job", out))) {
    stop("sbatch did not return a job id. Output: ",
         paste(out, collapse = "\n"))
  }
  jid <- parse_sbatch_output(out)
  if (is.null(jid) || !nzchar(jid)) {
    stop("could not parse job id from sbatch output: ",
         paste(out, collapse = "\n"))
  }
  jid
}

#' Query sacct for one job's current state. Returns one of:
#'   "pending" | "running" | "complete" | "failed" | "cancelled" | "unknown"
#'
#' Direct system2 call (same env story as .sbatch_submit above). Uses sacct
#' -X to exclude .extern/.batch substeps that confuse the standard
#' check_slurm_status() helper. Sacct works for both active and completed
#' jobs; squeue only returns active ones, which is why the SSH-centric
#' check_slurm_status() in helpers_search.R fails on this code path.
.sacct_state <- function(jid, ssh_config = NULL) {
  # Defensive: jid can arrive as NULL, NA, character(0), or "" depending on
  # how the status.json was last serialized by jsonlite.
  if (is.null(jid)) return("unknown")
  jid <- suppressWarnings(as.character(jid))
  if (length(jid) == 0) return("unknown")
  if (is.na(jid[1])) return("unknown")
  if (!nzchar(jid[1])) return("unknown")
  jid <- jid[1]

  out <- if (is.null(ssh_config)) {
    tryCatch(
      suppressWarnings(system2("sacct",
        args = c("-j", jid, "-X", "-n", "-o", "State"),
        stdout = TRUE, stderr = FALSE)),
      error = function(e) NULL
    )
  } else {
    res <- tryCatch(
      ssh_exec(ssh_config,
               sprintf("sacct -j %s -X -n -o State", shQuote(jid)),
               login_shell = TRUE, timeout = 15),
      error = function(e) NULL
    )
    if (is.null(res) || !identical(res$status, 0L)) NULL else res$stdout
  }

  if (is.null(out) || length(out) == 0 || !nzchar(trimws(out[1]))) {
    return("unknown")
  }
  state <- trimws(out[1])
  if (grepl("^COMPLETED",  state)) return("complete")
  if (grepl("^(RUNNING|COMPLETING)", state)) return("running")
  if (grepl("^PENDING",    state)) return("pending")
  if (grepl("^CANCELLED",  state)) return("cancelled")
  if (grepl("^(FAILED|TIMEOUT|OUT_OF_MEMORY|NODE_FAIL|BOOT_FAIL|PREEMPTED)",
            state)) return("failed")
  "unknown"
}

# =============================================================================
# status.json schema + I/O
# =============================================================================

.init_status_json <- function(project_dir, project_name, sample_names,
                              reference_key, tier_params, jids_by_stage,
                              build_metadata, ssh_config = NULL) {
  stages <- lapply(PROTEOG_STAGE_ORDER, function(s) {
    list(
      stage      = s,
      status     = "pending",
      job_id     = jids_by_stage[[s]] %||% NA_character_,
      started_at = NA_character_,
      finished_at = NA_character_
    )
  })

  status <- list(
    pipeline_id      = PROTEOG_PIPELINE_ID,
    project_name     = project_name,
    project_dir      = project_dir,
    sample_names     = as.list(sample_names),
    reference_key    = reference_key,
    read_length_tier = tier_params$tier,
    qc_gate_unique_pct = tier_params$qc_gate_unique_pct,
    submitted_at     = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    current_stage    = "fastp",
    stages           = stages,
    build_metadata   = build_metadata
  )

  status_path <- file.path(project_dir, "status.json")
  json_str <- jsonlite::toJSON(status, auto_unbox = TRUE, pretty = TRUE)
  .fs_write_text(as.character(json_str), status_path, ssh_config = ssh_config)
  invisible(status_path)
}

# =============================================================================
# Public: submit_proteogenomics_build
# =============================================================================

#' Submit the full proteogenomics build pipeline
#'
#' Generates sbatch scripts for every stage, submits them with
#' `--dependency=afterok` chaining, and writes a status.json that downstream
#' callers can poll.
#'
#' @param project_name      character — sanitized; used as project subdir name
#' @param rnaseq_dir        character — directory containing <sample>_R{1,2}.fastq.gz
#' @param reference_key     character — key in references/registry.json (e.g., "mm39_GRCm39")
#' @param sample_names      character vector
#' @param library_type      "polyA" | "ribo_depleted" | "stranded" | "unstranded"
#' @param strand_flag       "" | "--rf" | "--fr"
#' @param project_tag       character — passed to header rewriter; defaults to upper-cased project_name
#' @param uniprot_fasta     character or NULL — for Phase 4 assembly step; NULL skips merge
#' @param diamond_db        character or NULL — TransDecoder homology support
#' @param min_orf_len       integer — TransDecoder LongOrfs min length (default 100)
#' @param slurm_account     character (default "genome-center-grp")
#' @param slurm_partition   character (default "high")
#' @param ref_registry      list — if NULL, loaded from /quobyte/.../references/registry.json
#' @param rnaseq_root       character — base output dir
#' @return list with $project_dir, $status_path, $jids_by_stage, $tier_params,
#'         $pipeline_id, $methods_paragraph
submit_proteogenomics_build <- function(
  project_name,
  rnaseq_dir,
  reference_key,
  sample_names,
  library_type      = "polyA",
  strand_flag       = "",
  project_tag       = NULL,
  uniprot_fasta     = NULL,
  diamond_db        = NULL,
  min_orf_len       = 100L,
  slurm_account     = "genome-center-grp",
  slurm_partition   = "high",
  ref_registry      = NULL,
  rnaseq_root       = PROTEOG_RNASEQ_ROOT,
  ssh_config        = NULL
) {
  # ---- 1. Load reference registry + validate inputs --------------------------
  if (is.null(ref_registry)) {
    ref_registry <- load_reference_registry(ssh_config = ssh_config)
  }
  # NOTE: .validate_build_inputs() does file.exists() checks on the FASTQ
  # paths — these only work locally. Skip them when SSH-mode (the scan
  # handler already verified the directory + files exist on Hive).
  if (is.null(ssh_config)) {
    .validate_build_inputs(project_name, rnaseq_dir, sample_names,
                           reference_key, library_type, strand_flag,
                           ref_registry)
  }
  ref <- ref_registry[[reference_key]]

  # ---- 2. Set up project_dir -------------------------------------------------
  project_dir <- file.path(rnaseq_root, project_name)
  .fs_mkdir(project_dir, ssh_config = ssh_config)
  .fs_mkdir(file.path(project_dir, "logs"), ssh_config = ssh_config)
  .stage_rnaseq_inputs(project_dir, rnaseq_dir, sample_names,
                       ssh_config = ssh_config)

  # ---- 3. Detect read length on first R1 -------------------------------------
  first_r1 <- file.path(project_dir, "rnaseq",
                        sprintf("%s_R1.fastq.gz", sample_names[1]))
  read_len <- .fs_detect_read_length(first_r1, n_reads = 100L,
                                      ssh_config = ssh_config)
  if (is.na(read_len)) {
    stop("submit_proteogenomics_build(): could not detect read length from ", first_r1)
  }
  tier_params <- select_star_params(read_len)
  if (tier_params$tier == "refuse") {
    stop("submit_proteogenomics_build(): ", tier_params$error,
         " (detected median read length = ", read_len, " bp)")
  }

  if (is.null(project_tag)) {
    project_tag <- toupper(sanitize_project_name(project_name))
  }

  # ---- 4. Generate all sbatch scripts ----------------------------------------
  scripts <- list(
    fastp        = generate_fastp_sbatch(project_dir, sample_names,
                                          slurm_account, slurm_partition),
    rrna_filter  = generate_rrna_sbatch(project_dir, sample_names,
                                         ref$rrna_index,
                                         slurm_account, slurm_partition),
    star         = generate_star_sbatch(project_dir, sample_names,
                                         ref$star_index, tier_params,
                                         slurm_account, slurm_partition),
    qc_gate      = generate_qc_gate_sbatch(project_dir, sample_names,
                                            tier_params$qc_gate_unique_pct,
                                            slurm_account, slurm_partition),
    stringtie    = generate_stringtie_sbatch(project_dir, sample_names,
                                              ref$gtf, strand_flag,
                                              slurm_account, slurm_partition),
    merge        = generate_merge_sbatch(project_dir, ref$gtf, sample_names,
                                          slurm_account, slurm_partition),
    gffcompare   = generate_gffcompare_sbatch(project_dir, ref$gtf,
                                               slurm_account, slurm_partition),
    gffread      = generate_gffread_sbatch(project_dir, ref$genome_fasta,
                                            slurm_account, slurm_partition),
    transdecoder = generate_transdecoder_sbatch(project_dir, diamond_db,
                                                 min_orf_len,
                                                 slurm_account, slurm_partition),
    rewrite      = generate_rewrite_sbatch(project_dir, project_tag,
                                            slurm_account = slurm_account,
                                            slurm_partition = slurm_partition),
    assemble     = generate_assemble_sbatch(
                     project_dir, project_name,
                     uniprot_fasta = uniprot_fasta %||% "",
                     slurm_account = slurm_account,
                     slurm_partition = slurm_partition)
  )

  # Write each script to <project_dir>/sbatch/<stage>.sbatch
  sbatch_dir <- file.path(project_dir, "sbatch")
  .fs_mkdir(sbatch_dir, ssh_config = ssh_config)
  script_paths <- character()
  for (stage in names(scripts)) {
    p <- file.path(sbatch_dir, sprintf("%s.sbatch", stage))
    .fs_write_text(scripts[[stage]], p, ssh_config = ssh_config,
                   executable = TRUE)
    script_paths[[stage]] <- p
  }

  # ---- 5. Submit with afterok dependency chaining ---------------------------
  jids_by_stage <- list()
  prev <- NULL
  for (stage in names(scripts)) {
    jid <- .sbatch_submit(script_paths[[stage]], dep_jid = prev,
                          ssh_config = ssh_config)
    jids_by_stage[[stage]] <- jid
    prev <- jid
  }

  # ---- 6. Self-describing methods paragraph + status.json --------------------
  build_metadata <- list(
    rnaseq_dir       = rnaseq_dir,
    library_type     = library_type,
    strand_flag      = strand_flag,
    organism         = ref$organism %||% NA_character_,
    genome_build     = ref$build %||% NA_character_,
    annotation       = paste(c(ref$annotation_source, ref$annotation_release), collapse = " "),
    project_tag      = project_tag,
    detected_read_length = read_len
  )

  methods_paragraph <- sprintf(
    "Proteogenomics database built from %d samples (%s, %s; %s tier, read length %.0f bp). Pipeline: fastp → bowtie2 rRNA filter → STAR (%s) → stringtie → merge → gffcompare → gffread → TransDecoder → header rewrite (project tag %s). Reference: %s %s (%s %s).",
    length(sample_names), library_type, strand_flag,
    tier_params$tier, read_len, paste(tier_params$flags, collapse = " "),
    project_tag,
    build_metadata$organism, build_metadata$genome_build,
    ref$annotation_source %||% "", ref$annotation_release %||% ""
  )

  status_path <- .init_status_json(project_dir, project_name, sample_names,
                                    reference_key, tier_params,
                                    jids_by_stage, build_metadata,
                                    ssh_config = ssh_config)

  list(
    pipeline_id       = PROTEOG_PIPELINE_ID,
    project_dir       = project_dir,
    status_path       = status_path,
    jids_by_stage     = jids_by_stage,
    tier_params       = tier_params,
    build_metadata    = build_metadata,
    methods_paragraph = methods_paragraph
  )
}

# =============================================================================
# Public: submit_assemble_only — run JUST the assemble step for an existing
# build. Used by the per-row "Assemble" button in the Active Builds table
# for legacy builds that finished the SLURM chain before auto-assemble was
# wired up. Generates an assemble.sbatch, submits it, updates status.json.
# =============================================================================
submit_assemble_only <- function(project_dir,
                                 project_name,
                                 uniprot_fasta = "",
                                 slurm_account   = "genome-center-grp",
                                 slurm_partition = "high",
                                 ssh_config = NULL) {
  if (!nzchar(project_name)) stop("submit_assemble_only(): project_name required")

  sbatch_dir <- file.path(project_dir, "sbatch")
  .fs_mkdir(sbatch_dir, ssh_config = ssh_config)
  .fs_mkdir(file.path(project_dir, "logs"), ssh_config = ssh_config)

  script <- generate_assemble_sbatch(
    project_dir, project_name,
    uniprot_fasta = uniprot_fasta %||% "",
    slurm_account = slurm_account,
    slurm_partition = slurm_partition)
  script_path <- file.path(sbatch_dir, "assemble.sbatch")
  .fs_write_text(script, script_path, ssh_config = ssh_config,
                 executable = TRUE)

  jid <- .sbatch_submit(script_path, dep_jid = NULL, ssh_config = ssh_config)

  # Patch status.json: find the assemble stage, set job_id + status=pending
  status_path <- file.path(project_dir, "status.json")
  raw <- .fs_read_text(status_path, ssh_config = ssh_config)
  if (!is.null(raw) && nzchar(raw)) {
    status <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE),
                       error = function(e) NULL)
    if (!is.null(status) && is.list(status$stages)) {
      for (i in seq_along(status$stages)) {
        if (identical(status$stages[[i]]$stage, "assemble")) {
          status$stages[[i]]$job_id <- jid
          status$stages[[i]]$status <- "pending"
          status$stages[[i]]$started_at  <- NA_character_
          status$stages[[i]]$finished_at <- NA_character_
          break
        }
      }
      status$current_stage <- "assemble"
      # Record the uniprot input in build_metadata for traceability
      if (is.null(status$build_metadata)) status$build_metadata <- list()
      status$build_metadata$uniprot_fasta <- uniprot_fasta
      json_str <- jsonlite::toJSON(status, auto_unbox = TRUE, pretty = TRUE)
      .fs_write_text(as.character(json_str), status_path,
                     ssh_config = ssh_config)
    }
  }
  jid
}

# =============================================================================
# Public: poll_proteog_build_status
# =============================================================================

#' Refresh status.json by querying sacct for each stage's job_id
#'
#' @param project_dir character — from submit_*() result
#' @param ssh_config NULL (local) or ssh_config list (remote via ssh_exec)
#' @return updated status list with $current_stage and per-stage states
poll_proteog_build_status <- function(project_dir, ssh_config = NULL) {
  status_path <- file.path(project_dir, "status.json")
  raw <- .fs_read_text(status_path, ssh_config = ssh_config)
  if (is.null(raw)) {
    stop("poll_proteog_build_status(): status.json not found at ", status_path)
  }
  status <- tryCatch(jsonlite::fromJSON(raw, simplifyVector = FALSE),
                     error = function(e) NULL)
  if (is.null(status) || is.null(status$stages)) {
    stop("poll_proteog_build_status(): status.json missing $stages")
  }

  now <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  any_running    <- FALSE
  current_stage  <- "complete"
  any_failed     <- FALSE

  # See file-level .empty_or_str — re-declared here as a no-op for back-compat
  # (older code referenced this local). Kept so existing local references work.

  for (i in seq_along(status$stages)) {
    st <- status$stages[[i]]
    cur_status <- .empty_or_str(st$status)
    if (cur_status %in% c("complete", "failed", "cancelled")) next

    jid_str <- .empty_or_str(st$job_id)
    if (!nzchar(jid_str)) {
      # No SLURM job submitted for this stage (e.g., assemble is an
      # in-process consolidation step). Treat as still-pending so the
      # build doesn't get wrongly summarized as "complete".
      any_running <- TRUE
      if (current_stage == "complete") {
        current_stage <- .empty_or_str(st$stage)
        if (!nzchar(current_stage)) current_stage <- "pending"
      }
      next
    }

    new_state <- .sacct_state(jid_str, ssh_config = ssh_config)
    if (new_state == "running" && !nzchar(.empty_or_str(st$started_at))) {
      status$stages[[i]]$started_at <- now
    }
    if (new_state %in% c("complete", "failed", "cancelled") &&
        !nzchar(.empty_or_str(st$finished_at))) {
      status$stages[[i]]$finished_at <- now
    }
    status$stages[[i]]$status <- new_state
    if (new_state == "running" || new_state == "pending") {
      any_running <- TRUE
      if (current_stage == "complete") current_stage <- .empty_or_str(st$stage)
    }
    if (new_state == "failed" || new_state == "cancelled") {
      any_failed <- TRUE
    }
  }

  status$current_stage <- if (any_failed) "failed"
                          else if (any_running) current_stage
                          else "complete"
  status$last_polled_at <- now

  # ── Auto-register the assembled FASTA in the local FASTA library catalog ──
  # Fires once per build, on the transition from in-progress → complete.
  # Tracked via status$library_entry_id so we don't re-register on every poll.
  if (identical(status$current_stage, "complete") &&
      is.null(status$library_entry_id %||% NULL)) {
    entry_id <- tryCatch(
      .register_proteog_fasta_in_library(status, ssh_config = ssh_config),
      error = function(e) {
        message("[proteog] auto-register failed: ", conditionMessage(e))
        NULL
      })
    if (!is.null(entry_id) && is_scalar_char_safe(entry_id)) {
      status$library_entry_id <- entry_id
    }
  }

  json_str <- jsonlite::toJSON(status, auto_unbox = TRUE, pretty = TRUE)
  .fs_write_text(as.character(json_str), status_path, ssh_config = ssh_config)
  status
}

# Small scalar-char check used in poll + register helpers.
is_scalar_char_safe <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
}

#' Add the just-assembled proteogenomics FASTA to ~/.delimp_fasta_library/catalog.rds.
#'
#' Mirrors the schema used by the search-page FASTA library (see helpers_search.R
#' fasta_library_save) so proteog builds show up in the same picker modal.
#' Reads file size + sequence count from Hive via SSH so the entry is accurate.
.register_proteog_fasta_in_library <- function(status, ssh_config = NULL) {
  project_name <- status$project_name %||% basename(status$project_dir %||% "")
  if (!is_scalar_char_safe(project_name)) {
    stop("register: missing project_name")
  }
  date_tag <- format(Sys.Date(), "%Y_%m")
  fasta_name <- sprintf("%s_proteogenomics_%s.fasta", project_name, date_tag)
  fasta_path <- file.path(
    "/quobyte/proteomics-grp/de-limp/databases/proteogenomics", fasta_name)

  # File size + seq count from Hive (works for local too if path exists locally)
  size_bytes <- NA_integer_; seq_count <- NA_integer_
  if (is.null(ssh_config)) {
    if (file.exists(fasta_path)) {
      size_bytes <- as.integer(file.info(fasta_path)$size)
      seq_count  <- as.integer(length(grep("^>",
        readLines(fasta_path, warn = FALSE))))
    }
  } else {
    cmd <- sprintf(
      "stat -c%%s %s 2>/dev/null; echo SEP; grep -c '^>' %s 2>/dev/null",
      shQuote(fasta_path), shQuote(fasta_path))
    r <- tryCatch(ssh_exec(ssh_config, cmd, login_shell = FALSE, timeout = 30),
                  error = function(e) NULL)
    if (!is.null(r) && identical(r$status, 0L)) {
      txt <- paste(r$stdout %||% character(), collapse = "\n")
      parts <- strsplit(txt, "SEP", fixed = TRUE)[[1]]
      if (length(parts) >= 2) {
        size_bytes <- suppressWarnings(as.integer(trimws(parts[1])))
        seq_count  <- suppressWarnings(as.integer(trimws(parts[2])))
      }
    }
  }
  if (is.na(size_bytes) || size_bytes == 0L) {
    stop("register: FASTA not found or empty on Hive: ", fasta_path)
  }

  organism <- status$build_metadata$organism %||% NA_character_
  if (!is_scalar_char_safe(organism)) organism <- NA_character_

  entry_id <- sprintf("proteog_%s_%s",
                      project_name,
                      format(Sys.time(), "%Y%m%d_%H%M%S"))
  entry <- list(
    id                   = entry_id,
    name                 = sprintf("Proteogenomics: %s", project_name),
    organism             = organism,
    organism_common      = NA_character_,
    proteome_id          = NA_character_,
    content_type         = "proteogenomics",
    protein_count        = seq_count,
    file_size_bytes      = size_bytes,
    contaminant_library  = "None",
    contaminant_count    = 0L,
    custom_sequences     = NULL,
    custom_sequence_count = 0L,
    fasta_files          = fasta_name,
    fasta_dir            = "proteogenomics",
    remote_dir           = fasta_path,
    search_settings      = NULL,
    speclib_path         = NULL,
    speclib_search_mode  = NULL,
    created_at           = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),

    # Proteogenomics-specific metadata (visible in detail panel)
    proteog_pipeline_id        = PROTEOG_PIPELINE_ID,
    proteog_project_name       = project_name,
    proteog_project_dir        = status$project_dir %||% NA_character_,
    proteog_status_path        = file.path(status$project_dir %||% "",
                                           "status.json"),
    proteog_methods_paragraph  = status$build_metadata$methods_paragraph %||%
                                  NA_character_,
    proteog_sample_names       = status$sample_names %||% list(),
    proteog_reference_key      = status$reference_key %||% NA_character_,
    proteog_uniprot_fasta      = status$build_metadata$uniprot_fasta %||%
                                  NA_character_,
    proteog_read_length_tier   = status$read_length_tier %||% NA_character_
  )

  # Load catalog, de-dup by remote_dir (so re-poll doesn't re-add), save.
  catalog <- tryCatch(fasta_library_load(), error = function(e) list())
  if (!is.list(catalog)) catalog <- list()
  existing_idx <- which(vapply(catalog, function(c) {
    rd <- c$remote_dir %||% ""
    is_scalar_char_safe(rd) && identical(rd, fasta_path)
  }, logical(1)))
  if (length(existing_idx) > 0) {
    # Update in place rather than adding a duplicate
    catalog[[existing_idx[1]]] <- entry
  } else {
    catalog[[length(catalog) + 1L]] <- entry
  }
  fasta_library_save(catalog)
  message(sprintf("[proteog] registered FASTA in library: %s (%d seqs, %s bytes)",
                  fasta_name, seq_count, format(size_bytes, big.mark = ",")))
  entry_id
}

# =============================================================================
# Public: cancel_proteog_build
# =============================================================================

#' Scancel every non-terminal job in the build
cancel_proteog_build <- function(project_dir, ssh_config = NULL) {
  status <- poll_proteog_build_status(project_dir, ssh_config = ssh_config)
  to_cancel <- character()
  for (st in status$stages) {
    if (!is.null(st$job_id) && nzchar(st$job_id) &&
        !(st$status %in% c("complete", "failed", "cancelled"))) {
      to_cancel <- c(to_cancel, st$job_id)
    }
  }
  if (length(to_cancel) == 0) {
    return(invisible(list(cancelled = character(),
                          message = "no active jobs to cancel")))
  }
  out <- if (is.null(ssh_config)) {
    tryCatch(
      system2("scancel", args = to_cancel, stdout = TRUE, stderr = TRUE),
      error = function(e) stop("scancel failed: ", conditionMessage(e))
    )
  } else {
    cmd <- sprintf("scancel %s", paste(shQuote(to_cancel), collapse = " "))
    res <- tryCatch(
      ssh_exec(ssh_config, cmd, login_shell = TRUE, timeout = 30),
      error = function(e) stop("ssh scancel failed: ", conditionMessage(e))
    )
    res$stdout
  }
  invisible(list(
    cancelled = to_cancel,
    scancel_output = out
  ))
}


# =============================================================================
# Shiny module — UI renderer + server logic for the Build Database 🧬 tab
# =============================================================================
# Called from app.R as: server_proteog_builder(input, output, session, values)
#
# Hard-gated by hpc_available + !is_hf_space at the UI layer (see R/ui.R), so
# the server module assumes sbatch/squeue/sacct are reachable. If invoked on
# Docker-only, the UI panel never renders so the observers below are inert.

#' Render the Build Database tab body
#'
#' Five vertically-stacked accordion-style cards:
#'   1. Source — SLIMS URL OR comma-separated SRA accessions
#'   2. Sample scan / metadata verification
#'   3. Reference selection (from /quobyte/.../references/registry.json)
#'   4. Pipeline parameters (library type, strand, project name + tag)
#'   5. Submit
#'   plus an Active builds table polling status.json files
#'
#' All inputs are namespaced with `proteog_` to avoid colliding with the
#' existing Run Search inputs (analysis_name, fasta_source, etc.).
build_database_ui <- function() {
  div(
    style = "overflow-y: auto; max-height: calc(100vh - 150px); padding: 8px;",

    # ── Header
    div(
      style = "background: linear-gradient(135deg, #e8f5e9 0%, #c8e6c9 100%); padding: 12px 16px; border-radius: 8px; margin-bottom: 16px; display: flex; align-items: center; justify-content: space-between; gap: 12px;",
      div(
        h4(icon("dna"), " Proteogenomics — Build a sample-specific search database",
           style = "margin: 0; color: #1b5e20;"),
        p("Convert matched RNA-seq into a custom FASTA that contains your samples' ",
          "novel ORFs alongside the canonical reference proteome. The result appears ",
          "in the Run Search FASTA dropdown with a \U0001F9EC tag.",
          style = "margin: 4px 0 0 0; color: #2e7d32; font-size: 0.9em;")
      ),
      actionButton("proteog_explain_workflow_btn",
                   label = tagList(icon("circle-question"), " Explain this workflow"),
                   class = "btn-light",
                   style = "white-space: nowrap; border: 1px solid #1b5e20; color: #1b5e20;")
    ),

    # Helper: card header with a "?" info button on the right
    # local helper so we don't pollute the global namespace
    # (Shiny's `tagList` is fine here; bslib::card_header accepts arbitrary content)

    # ── Step 1: Source ──────────────────────────────────────────────────────
    bslib::card(
      bslib::card_header(div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        div(icon("upload"), " 1. RNA-seq Source"),
        actionButton("proteog_step1_info_btn", icon("question-circle"),
                     class = "btn-outline-info btn-sm",
                     title = "What does this step do?")
      )),
      bslib::card_body(
        radioButtons("proteog_source_mode", NULL,
                     choices = c("UC Davis DNA Tech Core (SLIMS URL)" = "slims",
                                 "Public archive (SRA/ENA accession)"  = "sra",
                                 "Folder of FASTQ files on the cluster" = "local"),
                     selected = "slims", inline = FALSE),
        conditionalPanel(
          "input.proteog_source_mode == 'slims'",
          textInput("proteog_slims_url", "SLIMS URL",
                    placeholder = "http://slimsdata.genomecenter.ucdavis.edu/Data/<id>/Unaligned/",
                    width = "100%"),
          helpText("For UC Davis users only. Paste the URL you received in your ",
                   "DNA Tech Core delivery email.")
        ),
        conditionalPanel(
          "input.proteog_source_mode == 'sra'",
          textInput("proteog_sra_accessions", "Accessions (comma-separated, max 24)",
                    placeholder = "SRR1303776, SRR1303777",
                    width = "100%"),
          checkboxInput("proteog_subsample",
                        "Stream-subsample to 5M read pairs per accession (fast test)",
                        value = FALSE),
          helpText("Use for re-analysis of published datasets or external SRA data. ",
                   "The pipeline streams FASTQs directly from ENA — no setup required.")
        ),
        conditionalPanel(
          "input.proteog_source_mode == 'local'",
          textInput("proteog_local_dir", "Directory path (on the cluster filesystem)",
                    placeholder = "/quobyte/proteomics-grp/myproject/rnaseq/",
                    width = "100%"),
          helpText("For data you already have on the cluster — from another sequencing ",
                   "core, a public download, or a collaborator. The folder must contain ",
                   "paired FASTQ files named ", tags$code("<sample>_R1.fastq.gz"), " / ",
                   tags$code("<sample>_R2.fastq.gz"), ". The Scan button below will list ",
                   "the samples it finds.")
        ),
        actionButton("proteog_scan_btn", "Scan / Verify",
                     icon = icon("magnifying-glass"),
                     class = "btn-outline-primary")
      )
    ),

    # ── Step 2: Scan / verification results ─────────────────────────────────
    bslib::card(
      bslib::card_header(div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        div(icon("clipboard-check"), " 2. Sample Verification"),
        actionButton("proteog_step2_info_btn", icon("question-circle"),
                     class = "btn-outline-info btn-sm",
                     title = "What is being verified?")
      )),
      bslib::card_body(uiOutput("proteog_scan_output"))
    ),

    # ── Step 3: Reference selection ─────────────────────────────────────────
    bslib::card(
      bslib::card_header(div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        div(icon("book-open"), " 3. Reference Genome"),
        actionButton("proteog_step3_info_btn", icon("question-circle"),
                     class = "btn-outline-info btn-sm",
                     title = "How do I pick the reference?")
      )),
      bslib::card_body(
        # Reference dropdown is rendered server-side via renderUI so its choices
        # are populated from load_reference_registry() with full reactive
        # context (the observe + updateSelectInput pattern raced renderUI on
        # initial session start, leaving "scanning…" stuck visible).
        uiOutput("proteog_reference_dropdown"),
        uiOutput("proteog_reference_info")
      )
    ),

    # ── Step 4: Pipeline parameters ─────────────────────────────────────────
    bslib::card(
      bslib::card_header(div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        div(icon("sliders"), " 4. Parameters"),
        actionButton("proteog_step4_info_btn", icon("question-circle"),
                     class = "btn-outline-info btn-sm",
                     title = "What do these parameters mean?")
      )),
      bslib::card_body(
        layout_columns(
          col_widths = c(6, 6),
          selectInput("proteog_library_type", "Library type",
                      choices = c("polyA mRNA-Seq"            = "polyA",
                                  "Total RNA + rRNA depletion" = "ribo_depleted",
                                  "Stranded RNA-Seq"           = "stranded",
                                  "Unstranded"                 = "unstranded"),
                      selected = "polyA"),
          selectInput("proteog_strand_flag", "Strand",
                      choices = c("Reverse stranded (TruSeq, --rf)" = "--rf",
                                  "Forward stranded (--fr)"          = "--fr",
                                  "Unstranded ()"                    = ""),
                      selected = "--rf")
        ),
        textInput("proteog_project_name", "Project name",
                  placeholder = "e.g. mouse_liver_pilot_2026_05",
                  width = "100%"),
        textInput("proteog_project_tag", "Project tag (uppercase, suffixes FASTA symbols)",
                  placeholder = "e.g. MOUSELIVER",
                  width = "100%"),
        numericInput("proteog_min_orf_len", "Minimum ORF length (aa)",
                     value = 100, min = 30, max = 300, step = 10, width = "50%"),
        tags$div(
          style = "border-top: 1px solid #e0e0e0; padding-top: 12px; margin-top: 12px;",
          tags$label("UniProt FASTA ", tags$em("(optional)"),
                     style = "font-weight: 600;"),
          tags$p(style = "color: #666; font-size: 0.85em; margin-bottom: 6px;",
                 "If provided, the final assemble step concatenates UniProt entries with the predicted ORFs ",
                 "so DIA-NN sees both. Leave on \"None\" to output predicted ORFs only."),
          selectInput("proteog_uniprot_source", label = NULL,
                      choices = c("None — predicted ORFs only" = "none",
                                  "Download from UniProt"      = "uniprot",
                                  "Download from NCBI"         = "ncbi",
                                  "Enter path on Hive"         = "path"),
                      selected = "none", width = "100%"),
          conditionalPanel("input.proteog_uniprot_source == 'uniprot'",
            actionButton("proteog_open_uniprot_modal", "Search UniProt",
                         class = "btn-info btn-sm w-100", icon = icon("search")),
            uiOutput("proteog_uniprot_selected_summary")
          ),
          conditionalPanel("input.proteog_uniprot_source == 'ncbi'",
            actionButton("proteog_open_ncbi_modal", "Search NCBI",
                         class = "btn-success btn-sm w-100", icon = icon("search")),
            uiOutput("proteog_ncbi_selected_summary")
          ),
          conditionalPanel("input.proteog_uniprot_source == 'path'",
            textInput("proteog_uniprot_fasta_path", label = NULL,
                      placeholder = "/quobyte/proteomics-grp/de-limp/databases/uniprot/UP000005640.fasta",
                      width = "100%")
          )
        )
      )
    ),

    # ── Step 5: Submit ──────────────────────────────────────────────────────
    bslib::card(
      bslib::card_header(div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        div(icon("rocket"), " 5. Submit"),
        actionButton("proteog_step5_info_btn", icon("question-circle"),
                     class = "btn-outline-info btn-sm",
                     title = "What happens when I click Submit?")
      )),
      bslib::card_body(
        uiOutput("proteog_submit_warnings"),
        actionButton("proteog_submit_btn", "Build Proteogenomics FASTA",
                     icon = icon("dna"),
                     class = "btn-success btn-lg w-100"),
        helpText("Estimated wall time: 3-6 hours for 12 samples × 30M PE150 reads. ",
                 "You can close the browser; the build continues on Hive.")
      )
    ),

    # ── Active builds table ─────────────────────────────────────────────────
    bslib::card(
      bslib::card_header(div(
        style = "display: flex; align-items: center; justify-content: space-between;",
        div(icon("list-check"), " Active & recent builds"),
        div(
          actionButton("proteog_restore_builds_btn",
                       label = tagList(icon("rotate"), "Restore from Hive"),
                       class = "btn-outline-secondary btn-sm",
                       title = "Scan Hive for in-progress builds and re-populate this list"),
          actionButton("proteog_builds_info_btn", icon("question-circle"),
                       class = "btn-outline-info btn-sm",
                       title = "What do the stage names mean?",
                       style = "margin-left: 6px;")
        )
      )),
      bslib::card_body(uiOutput("proteog_active_builds_table"))
    )
  )
}

#' Shiny server module for the Build Database tab
#'
#' Wires the UI inputs to:
#'   - load_reference_registry() / load_proteog_registry() for dropdowns
#'   - scan_slims_url() / verify_sra_accession() for Step 1 verify
#'   - submit_proteogenomics_build() for the submit button
#'   - reactivePoll on status.json files for the active-builds table
server_proteog_builder <- function(input, output, session, values) {

  ns <- session$ns %||% function(x) x  # tolerant of being called inside a moduleServer or not

  # ── Active SSH config for SSH-aware helpers ────────────────────────────────
  # When the user has connected to Hive via the Run Search tab,
  # `values$ssh_connected` is TRUE and the SSH inputs (input$ssh_host,
  # input$ssh_user, input$ssh_port, input$ssh_key_path, input$ssh_modules)
  # are populated. We re-construct the ssh_config list here rather than
  # depend on a cross-module reactive — matches the construction used in
  # `R/server_search.R:23` so the same connection state is reused.
  proteog_ssh_config <- function() {
    if (!isTRUE(values$ssh_connected)) return(NULL)
    if (is.null(input$ssh_host) || !nzchar(input$ssh_host %||% "")) return(NULL)
    list(
      host     = input$ssh_host,
      user     = input$ssh_user,
      port     = input$ssh_port %||% 22,
      key_path = input$ssh_key_path,
      modules  = input$ssh_modules %||% ""
    )
  }

  # ── Active-builds persistence ────────────────────────────────────────────
  # values$proteog_build_jobs lives in reactiveValues (in-memory only).
  # We mirror it to a local RDS so the Active Builds list survives Shiny
  # restarts and the user can resume tracking after closing R. Per CLAUDE.md
  # "never use mounted drives for app state" — local path only.
  proteog_builds_path <- function() {
    file.path(path.expand("~"), ".delimp_proteog_builds.rds")
  }
  proteog_save_builds <- function(jobs) {
    tryCatch(saveRDS(jobs, proteog_builds_path()),
             error = function(e) message("[proteog] save failed: ",
                                         conditionMessage(e)))
  }
  proteog_load_builds <- function() {
    p <- proteog_builds_path()
    if (!file.exists(p)) return(list())
    tryCatch(readRDS(p), error = function(e) {
      message("[proteog] load failed: ", conditionMessage(e)); list()
    })
  }

  # Restore on module init — only if not already populated (avoid overwriting
  # an in-progress session restore from server_session.R). Wrap in isolate()
  # because reactiveValues reads outside a reactive consumer throw.
  isolate({
    if (length(values$proteog_build_jobs %||% list()) == 0) {
      restored <- proteog_load_builds()
      if (length(restored) > 0) {
        values$proteog_build_jobs <- restored
        message(sprintf("[proteog] restored %d active build(s) from disk",
                        length(restored)))
      }
    }
  })

  # Persist on every change — but never clobber the file with an empty list,
  # which would happen on session startup before the user has submitted
  # anything in this session. The empty-on-startup case must NOT erase
  # builds from prior sessions.
  observe({
    jobs <- values$proteog_build_jobs %||% list()
    if (length(jobs) == 0 && file.exists(proteog_builds_path())) {
      existing <- tryCatch(readRDS(proteog_builds_path()),
                           error = function(e) NULL)
      if (length(existing) > 0) return()
    }
    proteog_save_builds(jobs)
  })

  # ── UniProt / NCBI FASTA download modals (proteog-specific) ─────────────
  # Reuses helper functions from helpers_search.R (search_uniprot_proteomes,
  # download_uniprot_fasta, download_ncbi_fasta) but routes results into
  # proteog state to avoid clobbering the main Search page's values$fasta_info.
  proteog_uniprot_state <- reactiveValues(
    results = NULL, hive_path = NULL, summary = NULL
  )
  proteog_ncbi_state <- reactiveValues(
    results = NULL, hive_path = NULL, summary = NULL
  )

  PROTEOG_UNIPROT_CACHE <- "/quobyte/proteomics-grp/de-limp/databases/uniprot"

  observeEvent(input$proteog_open_uniprot_modal, {
    showModal(modalDialog(
      title = tagList(icon("dna"), " UniProt — pick a proteome for the assemble step"),
      size = "l", easyClose = TRUE,
      div(style = "display: flex; gap: 8px; margin-bottom: 12px;",
        div(style = "flex: 1;",
          textInput("proteog_uniprot_query", NULL,
                    placeholder = "e.g., human, mouse, E. coli", width = "100%")),
        actionButton("proteog_search_uniprot_btn", "Search",
                     class = "btn-info", style = "margin-top: 0;")),
      DT::DTOutput("proteog_uniprot_results_table"),
      hr(),
      selectInput("proteog_uniprot_content_type", "Content:",
                  choices = c("One per gene (recommended)" = "one_per_gene",
                              "Swiss-Prot reviewed"        = "reviewed",
                              "Swiss-Prot + isoforms"      = "reviewed_isoforms",
                              "Full proteome"              = "full"),
                  selected = "one_per_gene", width = "100%"),
      footer = tagList(modalButton("Cancel"),
        actionButton("proteog_uniprot_download_btn", "Download + upload to Hive",
                     class = "btn-success", icon = icon("download")))))
  })

  observeEvent(input$proteog_search_uniprot_btn, {
    req(nzchar(input$proteog_uniprot_query %||% ""))
    withProgress(message = "Searching UniProt…", {
      proteog_uniprot_state$results <- tryCatch(
        search_uniprot_proteomes(input$proteog_uniprot_query),
        error = function(e) { showNotification(
          sprintf("UniProt search failed: %s", conditionMessage(e)),
          type = "error", duration = 8); data.frame() })
    })
    if (is.null(proteog_uniprot_state$results) ||
        nrow(proteog_uniprot_state$results) == 0) {
      showNotification("No proteomes found.", type = "warning", duration = 5)
    }
  })

  output$proteog_uniprot_results_table <- DT::renderDT({
    req(proteog_uniprot_state$results, nrow(proteog_uniprot_state$results) > 0)
    df <- proteog_uniprot_state$results[, c("upid", "organism", "common_name",
                                             "protein_count")]
    colnames(df) <- c("ID", "Organism", "Common Name", "Proteins")
    DT::datatable(df, selection = "single", rownames = FALSE,
      options = list(pageLength = 10, dom = "tip", scrollY = "300px"),
      class = "compact stripe")
  })

  observeEvent(input$proteog_uniprot_download_btn, {
    sel <- input$proteog_uniprot_results_table_rows_selected
    req(length(sel) > 0, proteog_uniprot_state$results)
    row <- proteog_uniprot_state$results[sel, ]
    sc <- proteog_ssh_config()
    if (is.null(sc)) {
      showNotification("Connect to Hive first.", type = "warning", duration = 5)
      return()
    }
    fname <- tryCatch(
      generate_fasta_filename(row$upid, row$organism,
                              input$proteog_uniprot_content_type),
      error = function(e) sprintf("%s_%s.fasta", row$upid,
                                  input$proteog_uniprot_content_type))
    hive_path <- file.path(PROTEOG_UNIPROT_CACHE, fname)
    # Reuse-if-cached: skip download when the file is already on Hive
    cached <- tryCatch(
      ssh_exec(sc, sprintf("test -s %s && echo OK", shQuote(hive_path)),
               login_shell = FALSE, timeout = 10),
      error = function(e) NULL)
    if (!is.null(cached) && identical(cached$status, 0L) &&
        any(grepl("OK", cached$stdout %||% character()))) {
      proteog_uniprot_state$hive_path <- hive_path
      proteog_uniprot_state$summary <- sprintf(
        "%s — cached on Hive (%s)", fname, row$organism)
      showNotification(sprintf("Already cached on Hive: %s", fname),
                       type = "default", duration = 6)
      removeModal()
      .maybe_auto_assemble_after_download(hive_path)
      return()
    }
    withProgress(message = "Downloading UniProt FASTA…", {
      tmp_local <- tempfile(pattern = "proteog_uniprot_", fileext = ".fasta")
      res <- tryCatch(download_uniprot_fasta(row$upid,
                                              input$proteog_uniprot_content_type,
                                              tmp_local),
                      error = function(e) {
                        showNotification(sprintf("Download failed: %s",
                                                  conditionMessage(e)),
                                          type = "error", duration = 10)
                        NULL
                      })
      if (is.null(res) || !file.exists(tmp_local)) return()
      setProgress(0.6, message = "Uploading to Hive…")
      mk <- ssh_exec(sc, sprintf("mkdir -p %s", shQuote(PROTEOG_UNIPROT_CACHE)),
                     login_shell = FALSE, timeout = 15)
      if (!identical(mk$status, 0L)) {
        showNotification("Failed to create Hive cache dir.",
                         type = "error", duration = 8); return()
      }
      up <- scp_upload(sc, tmp_local, hive_path, timeout = 300)
      if (!identical(up$status, 0L)) {
        showNotification("SCP upload failed.", type = "error", duration = 8)
        return()
      }
      proteog_uniprot_state$hive_path <- hive_path
      proteog_uniprot_state$summary <- sprintf("%s — %s (%d proteins)",
                                                fname, row$organism,
                                                as.integer(row$protein_count %||% 0L))
      file.remove(tmp_local)
    })
    showNotification(sprintf("Uploaded UniProt FASTA to %s", hive_path),
                     type = "default", duration = 8)
    removeModal()
    .maybe_auto_assemble_after_download(hive_path)
  })

  output$proteog_uniprot_selected_summary <- renderUI({
    if (is.null(proteog_uniprot_state$summary)) {
      return(helpText("No UniProt proteome selected yet."))
    }
    div(class = "alert alert-success py-2 px-3 mt-2",
        style = "font-size: 0.85em;",
        icon("check"), tags$strong(" Selected: "),
        proteog_uniprot_state$summary,
        tags$br(),
        tags$small(tags$code(proteog_uniprot_state$hive_path)))
  })

  # NCBI modal — same pattern, calls download_ncbi_fasta()
  observeEvent(input$proteog_open_ncbi_modal, {
    showModal(modalDialog(
      title = tagList(icon("dna"), " NCBI — pick a proteome for the assemble step"),
      size = "l", easyClose = TRUE,
      div(style = "display: flex; gap: 8px; margin-bottom: 12px;",
        div(style = "flex: 1;",
          textInput("proteog_ncbi_query", NULL,
                    placeholder = "e.g., Peromyscus, Bos taurus", width = "100%")),
        actionButton("proteog_search_ncbi_btn", "Search",
                     class = "btn-success", style = "margin-top: 0;")),
      DT::DTOutput("proteog_ncbi_results_table"),
      footer = tagList(modalButton("Cancel"),
        actionButton("proteog_ncbi_download_btn", "Download + upload to Hive",
                     class = "btn-success", icon = icon("download")))))
  })

  observeEvent(input$proteog_search_ncbi_btn, {
    req(nzchar(input$proteog_ncbi_query %||% ""))
    withProgress(message = "Searching NCBI…", {
      proteog_ncbi_state$results <- tryCatch(
        ncbi_search_assemblies(input$proteog_ncbi_query),
        error = function(e) { showNotification(
          sprintf("NCBI search failed: %s", conditionMessage(e)),
          type = "error", duration = 8); data.frame() })
    })
    if (is.null(proteog_ncbi_state$results) ||
        nrow(proteog_ncbi_state$results) == 0) {
      showNotification("No proteomes found.", type = "warning", duration = 5)
    }
  })

  output$proteog_ncbi_results_table <- DT::renderDT({
    req(proteog_ncbi_state$results, nrow(proteog_ncbi_state$results) > 0)
    df <- proteog_ncbi_state$results
    out <- data.frame(
      Accession = df$accession %||% "",
      Organism  = df$organism %||% "",
      Level     = df$assembly_level %||% "",
      Proteins  = format(as.integer(df$protein_count %||% 0L), big.mark = ","),
      Category  = df$refseq_category %||% "",
      stringsAsFactors = FALSE)
    DT::datatable(out, selection = "single", rownames = FALSE,
      options = list(pageLength = 10, dom = "tip", scrollY = "300px",
                     columnDefs = list(list(width = "120px", targets = 0))),
      class = "compact stripe")
  })

  observeEvent(input$proteog_ncbi_download_btn, {
    sel <- input$proteog_ncbi_results_table_rows_selected
    req(length(sel) > 0, proteog_ncbi_state$results)
    row <- proteog_ncbi_state$results[sel, ]
    sc <- proteog_ssh_config()
    if (is.null(sc)) {
      showNotification("Connect to Hive first.", type = "warning", duration = 5)
      return()
    }
    acc <- row$accession
    # ncbi_download_proteome takes a directory and generates the filename
    # plus a {basename}_gene_map.tsv (NCBI accessions don't carry gene symbols
    # in the FASTA header — DIA-NN's Genes column comes from this map TSV).
    tmp_local_dir <- tempfile(pattern = "proteog_ncbi_dl_")
    dir.create(tmp_local_dir, recursive = TRUE, showWarnings = FALSE)
    withProgress(message = sprintf("Downloading %s from NCBI…", row$organism), {
      local_fasta <- tryCatch(ncbi_download_proteome(acc, tmp_local_dir),
                              error = function(e) {
                                showNotification(sprintf("Download failed: %s",
                                                         conditionMessage(e)),
                                                 type = "error", duration = 10)
                                NULL
                              })
      if (is.null(local_fasta) || !file.exists(local_fasta)) return()
      # Find the gene_map.tsv generated alongside the FASTA (if any)
      gene_map <- sub("\\.fasta$", "_gene_map.tsv", local_fasta)
      fasta_name <- basename(local_fasta)
      gene_map_name <- basename(gene_map)
      hive_fasta <- file.path(PROTEOG_UNIPROT_CACHE, fasta_name)
      hive_gene_map <- file.path(PROTEOG_UNIPROT_CACHE, gene_map_name)

      setProgress(0.6, message = "Uploading FASTA + gene map to Hive…")
      ssh_exec(sc, sprintf("mkdir -p %s", shQuote(PROTEOG_UNIPROT_CACHE)),
               login_shell = FALSE, timeout = 15)
      up <- scp_upload(sc, local_fasta, hive_fasta, timeout = 300)
      if (!identical(up$status, 0L)) {
        showNotification("FASTA SCP upload failed.", type = "error", duration = 8)
        return()
      }
      if (file.exists(gene_map)) {
        scp_upload(sc, gene_map, hive_gene_map, timeout = 60)
      }
      proteog_ncbi_state$hive_path <- hive_fasta
      proteog_ncbi_state$summary <- sprintf("%s — %s (RefSeq accessions; gene map: %s)",
        fasta_name, row$organism,
        if (file.exists(gene_map)) "uploaded" else "missing")
      unlink(tmp_local_dir, recursive = TRUE)
    })
    showNotification(sprintf("Uploaded NCBI FASTA to %s", proteog_ncbi_state$hive_path),
                     type = "default", duration = 8)
    removeModal()
    .maybe_auto_assemble_after_download(proteog_ncbi_state$hive_path)
  })

  output$proteog_ncbi_selected_summary <- renderUI({
    if (is.null(proteog_ncbi_state$summary)) {
      return(helpText("No NCBI proteome selected yet."))
    }
    div(class = "alert alert-success py-2 px-3 mt-2",
        style = "font-size: 0.85em;",
        icon("check"), tags$strong(" Selected: "),
        proteog_ncbi_state$summary,
        tags$br(),
        tags$small(tags$code(proteog_ncbi_state$hive_path)))
  })

  # When a UniProt/NCBI download finishes WHILE the user is in the per-row
  # Assemble flow (proteog_assemble_target is set), auto-submit the assemble
  # job with the just-downloaded FASTA instead of asking the user to click
  # "Submit Assemble" with the path already filled in.
  .maybe_auto_assemble_after_download <- function(hive_path) {
    target <- proteog_assemble_target()
    if (is.null(target)) return(invisible())  # not in assemble flow → no-op
    if (!is_scalar_char_safe(hive_path)) return(invisible())
    sc <- proteog_ssh_config()
    if (is.null(sc)) return(invisible())
    tryCatch({
      jid <- submit_assemble_only(
        project_dir   = target$project_dir,
        project_name  = target$project_name,
        uniprot_fasta = hive_path,
        ssh_config    = sc)
      showNotification(
        sprintf("Auto-submitted assemble for %s with downloaded FASTA (SLURM %s)",
                target$project_name, jid),
        type = "default", duration = 10)
      proteog_assemble_target(NULL)
      removeModal()
    }, error = function(e) {
      showNotification(
        sprintf("Auto-assemble failed: %s — open Assemble modal manually to retry",
                conditionMessage(e)),
        type = "error", duration = 15)
    })
  }

  # Helper to resolve the chosen FASTA path based on the source dropdown.
  resolve_proteog_uniprot_path <- function() {
    src <- input$proteog_uniprot_source %||% "none"
    if (identical(src, "none"))     return("")
    if (identical(src, "path"))     return(input$proteog_uniprot_fasta_path %||% "")
    if (identical(src, "uniprot"))  return(proteog_uniprot_state$hive_path %||% "")
    if (identical(src, "ncbi"))     return(proteog_ncbi_state$hive_path %||% "")
    ""
  }

  # ── Restore-from-Hive: scan PROTEOG_RNASEQ_ROOT for status.json files,
  # read each, build a job entry, merge with current values$proteog_build_jobs.
  observeEvent(input$proteog_restore_builds_btn, {
    tryCatch({
    sc <- proteog_ssh_config()
    if (is.null(sc)) {
      showNotification("Connect to Hive first (Run Search → Test Connection).",
                       type = "warning", duration = 8)
      return()
    }
    # 1. find all status.json under PROTEOG_RNASEQ_ROOT (one per build)
    find_cmd <- sprintf(
      "find %s -maxdepth 2 -name status.json -type f 2>/dev/null",
      shQuote(PROTEOG_RNASEQ_ROOT))
    res <- tryCatch(
      ssh_exec(sc, find_cmd, login_shell = FALSE, timeout = 30),
      error = function(e) list(status = 1L, stdout = character(),
                               stderr = conditionMessage(e)))
    if (!identical(res$status, 0L)) {
      showNotification(sprintf("Scan failed: %s",
                               paste(res$stderr %||% "", collapse = "; ")),
                       type = "error", duration = 10)
      return()
    }
    raw_stdout <- paste(res$stdout %||% character(), collapse = "\n")
    paths <- trimws(unlist(strsplit(raw_stdout, "\n")))
    paths <- paths[!is.na(paths) & nzchar(paths)]
    if (length(paths) == 0) {
      showNotification("No builds found on Hive.",
                       type = "default", duration = 6)
      return()
    }
    # 2. read each, parse, derive jids_by_stage from stages list
    # Helper: TRUE iff x is a single non-NA, non-empty character.
    is_scalar_char <- function(x) {
      is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
    }
    restored <- list()
    message(sprintf("[proteog-restore] scanning %d path(s)", length(paths)))
    for (i in seq_along(paths)) {
      p <- paths[i]
      message(sprintf("[proteog-restore] [%d/%d] %s", i, length(paths), p))
      txt <- tryCatch(.fs_read_text(p, ssh_config = sc),
                      error = function(e) NULL)
      if (!is_scalar_char(txt)) { message("  skip: empty/NA txt"); next }
      parsed <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE),
                        error = function(e) {
                          message("  skip: fromJSON error: ", conditionMessage(e))
                          NULL
                        })
      if (is.null(parsed)) next
      pdir <- parsed$project_dir
      if (!is_scalar_char(pdir)) {
        message("  skip: project_dir not scalar char (class=",
                paste(class(pdir), collapse=","), ", len=", length(pdir), ")")
        next
      }
      jids_by_stage <- list()
      stages_in <- parsed$stages
      if (!is.list(stages_in)) stages_in <- list()
      for (s in stages_in) {
        if (!is.list(s)) next
        stage_name <- s$stage
        if (!is_scalar_char(stage_name)) next
        job_id_raw <- s$job_id
        if (is.null(job_id_raw)) next
        job_id_chr <- tryCatch(as.character(job_id_raw)[1],
                               error = function(e) NA_character_)
        if (!is_scalar_char(job_id_chr)) next
        jids_by_stage[[stage_name]] <- job_id_chr
      }
      proj_name <- parsed$project_name
      if (!is_scalar_char(proj_name)) proj_name <- basename(dirname(p))
      sub_at <- parsed$submitted_at
      if (!is_scalar_char(sub_at)) sub_at <- NA_character_
      restored[[length(restored) + 1L]] <- list(
        project_name      = proj_name,
        project_dir       = pdir,
        submitted_at      = sub_at,
        jids_by_stage     = jids_by_stage,
        methods_paragraph = parsed$build_metadata$methods_paragraph %||% NULL,
        download_pending  = FALSE
      )
      message("  added")
    }
    if (length(restored) == 0) {
      showNotification("Found status.json files but none parsed cleanly.",
                       type = "warning", duration = 8)
      return()
    }
    # 3. Merge into values$proteog_build_jobs — de-duplicate by project_dir
    current <- values$proteog_build_jobs %||% list()
    current_dirs <- vapply(current, function(j) {
      v <- j$project_dir
      if (is.null(v)) return("")
      v1 <- tryCatch(as.character(v)[1], error = function(e) "")
      if (length(v1) == 0 || is.na(v1) || !nzchar(v1)) "" else v1
    }, character(1))
    added <- 0L
    for (r in restored) {
      rd <- as.character(r$project_dir %||% "")
      if (is.na(rd) || !nzchar(rd)) next
      if (isTRUE(rd %in% current_dirs)) next
      current[[length(current) + 1L]] <- r
      added <- added + 1L
    }
    values$proteog_build_jobs <- current
    showNotification(sprintf("Restored %d build(s) from Hive (%d scanned).",
                             added, length(restored)),
                     type = "default", duration = 8)
    }, error = function(e) {
      showNotification(sprintf("Restore failed: %s",
                               conditionMessage(e)),
                       type = "error", duration = 15)
    })
  })

  # ── Skip-if-present helper: check whether the project_dir on Hive already
  # contains FASTQ files. If yes, skip download entirely and go straight to
  # the build submit.
  rnaseq_data_already_present <- function(project_dir, ssh_config) {
    if (is.null(ssh_config)) {
      if (!dir.exists(project_dir)) return(FALSE)
      files <- list.files(project_dir, pattern = "\\.(fastq|fq)\\.gz$",
                          recursive = TRUE, full.names = TRUE)
      return(length(files) > 0)
    }
    cmd <- sprintf(
      "test -d %s && find %s -maxdepth 3 \\( -name '*.fastq.gz' -o -name '*.fq.gz' \\) 2>/dev/null | head -1",
      shQuote(project_dir), shQuote(project_dir))
    res <- tryCatch(
      ssh_exec(ssh_config, cmd, login_shell = FALSE, timeout = 15),
      error = function(e) list(status = 1L, stdout = character()))
    nzchar(trimws(paste(res$stdout %||% character(), collapse = "")))
  }

  # ── render the Build Database body via uiOutput("build_database_content") ──
  output$build_database_content <- renderUI({ build_database_ui() })

  # ── render Reference Genome dropdown ────────────────────────────────────
  # Uses renderUI rather than the observe+updateSelectInput pattern so the
  # selectInput is built with the right choices at the moment the tab renders
  # (the prior observe approach raced the parent renderUI and left the
  # placeholder "scanning..." stuck visible).
  output$proteog_reference_dropdown <- renderUI({
    sc  <- proteog_ssh_config()
    message(sprintf(
      "[proteog] reference dropdown renderUI firing — ssh_connected=%s, sc_null=%s",
      isTRUE(values$ssh_connected), is.null(sc)
    ))
    reg <- tryCatch(load_reference_registry(ssh_config = sc),
                    error = function(e) {
                      message("[proteog] load_reference_registry error: ",
                              conditionMessage(e))
                      list()
                    })
    message(sprintf("[proteog] registry length = %d", length(reg)))

    if (length(reg) == 0) {
      msg <- if (is.null(sc)) {
        "Connect to HPC from the Run Search tab to see references"
      } else {
        "No references registered (registry empty or unreachable)"
      }
      empty_choices <- setNames("", msg)
      return(tagList(
        selectInput("proteog_reference_key", "Reference",
                    choices = empty_choices,
                    selected = "",
                    width = "100%"),
        tags$small(style = "color:#c62828;", msg)
      ))
    }

    labels <- vapply(names(reg), function(k) {
      e <- reg[[k]]
      sprintf("%s — %s %s", k,
              e$organism %||% "?",
              e$annotation_release %||% e$annotation_source %||% "")
    }, character(1), USE.NAMES = FALSE)
    choices <- setNames(names(reg), labels)
    selectInput("proteog_reference_key", "Reference",
                choices = choices, width = "100%")
  })

  output$proteog_reference_info <- renderUI({
    req(input$proteog_reference_key, nzchar(input$proteog_reference_key))
    reg <- load_reference_registry(ssh_config = proteog_ssh_config())
    entry <- reg[[input$proteog_reference_key]]
    if (is.null(entry)) return(NULL)
    tags$div(
      style = "background: #f5f5f5; padding: 8px; border-radius: 4px; font-size: 0.85em;",
      tags$div(strong("Genome FASTA:"), code(entry$genome_fasta %||% "(none)")),
      tags$div(strong("STAR index:"),   code(entry$star_index %||% "(none)")),
      tags$div(strong("GTF:"),          code(entry$gtf %||% "(none)")),
      tags$div(strong("rRNA index:"),   code(entry$rrna_index %||% "(none)")),
      tags$div(strong("Completeness:"), entry$completeness %||% "(unknown)")
    )
  })

  # ── Step 1: Scan / verify ────────────────────────────────────────────────
  # Stored in a reactiveVal so the UI in Step 2 can render it after the user
  # clicks Scan, and the submit handler can read it back to populate sample
  # names without re-fetching.
  scan_result <- reactiveVal(NULL)

  observeEvent(input$proteog_scan_btn, {
    if (input$proteog_source_mode == "slims") {
      url <- trimws(input$proteog_slims_url %||% "")
      if (!nzchar(url)) {
        scan_result(list(success = FALSE,
                         error = "Please enter a SLIMS URL."))
        return()
      }
      withProgress(message = "Scanning SLIMS URL…", value = 0.5, {
        res <- tryCatch(scan_slims_url(url),
                        error = function(e) list(success = FALSE,
                                                 error = conditionMessage(e)))
      })
      res$mode <- "slims"
      scan_result(res)

    } else if (input$proteog_source_mode == "local") {
      dir_path <- trimws(input$proteog_local_dir %||% "")
      if (!nzchar(dir_path)) {
        scan_result(list(success = FALSE,
                         error = "Please enter a directory path."))
        return()
      }
      sc <- proteog_ssh_config()

      # ── SSH path: list FASTQs via ssh_exec ───────────────────────────────
      if (!is.null(sc)) {
        # 1. Verify the directory exists on Hive
        test_cmd <- sprintf("test -d %s && echo OK", shQuote(dir_path))
        res_dir <- tryCatch(ssh_exec(sc, test_cmd, login_shell = FALSE, timeout = 10),
                            error = function(e) list(status = -1, stdout = character()))
        if (!identical(res_dir$status, 0L) ||
            !any(grepl("OK", res_dir$stdout, fixed = TRUE))) {
          scan_result(list(success = FALSE,
                           error = sprintf("Directory not found on Hive: %s", dir_path)))
          return()
        }
        # 2. List R1/R2 files; ls -1 outputs one filename per line
        ls_cmd <- sprintf("ls -1 %s 2>/dev/null | grep -E '_R[12]\\.fastq\\.gz$' || true",
                          shQuote(dir_path))
        res_ls <- tryCatch(ssh_exec(sc, ls_cmd, login_shell = FALSE, timeout = 15),
                           error = function(e) list(status = -1, stdout = character()))
        files <- if (identical(res_ls$status, 0L)) trimws(res_ls$stdout) else character()
        files <- files[nzchar(files)]
        r1_files <- files[grepl("_R1\\.fastq\\.gz$", files)]
        r2_files <- files[grepl("_R2\\.fastq\\.gz$", files)]
        if (length(r1_files) == 0) {
          scan_result(list(success = FALSE,
                           error = sprintf("No _R1.fastq.gz files in %s on Hive.",
                                           dir_path)))
          return()
        }
        sample_names <- sub("_R1\\.fastq\\.gz$", "", r1_files)
        r2_present <- sub("_R2\\.fastq\\.gz$", "", r2_files)
        missing_r2 <- setdiff(sample_names, r2_present)
        if (length(missing_r2) > 0) {
          scan_result(list(success = FALSE,
                           error = sprintf("%d sample(s) missing R2 on Hive: %s",
                                           length(missing_r2),
                                           paste(head(missing_r2, 5), collapse = ", "))))
          return()
        }
        scan_result(list(
          success = TRUE,
          mode = "local",
          local_dir = dir_path,
          n_samples = length(sample_names),
          sample_names = sample_names,
          is_paired = TRUE,
          has_md5 = FALSE  # not bothering to check via SSH
        ))
        return()
      }

      # ── Local-filesystem path (DE-LIMP-on-Hive case) ─────────────────────
      if (!dir.exists(dir_path)) {
        scan_result(list(success = FALSE,
                         error = sprintf("Directory not found: %s. Check the path is accessible from the cluster.", dir_path)))
        return()
      }
      r1_files <- list.files(dir_path, pattern = "_R1\\.fastq\\.gz$", full.names = FALSE)
      if (length(r1_files) == 0) {
        scan_result(list(success = FALSE,
                         error = sprintf(
                           "No _R1.fastq.gz files in %s.", dir_path)))
        return()
      }
      sample_names <- sub("_R1\\.fastq\\.gz$", "", r1_files)
      missing_r2 <- character()
      for (s in sample_names) {
        if (!file.exists(file.path(dir_path, sprintf("%s_R2.fastq.gz", s)))) {
          missing_r2 <- c(missing_r2, s)
        }
      }
      if (length(missing_r2) > 0) {
        scan_result(list(success = FALSE,
                         error = sprintf(
                           "%d sample(s) missing matching _R2.fastq.gz: %s.",
                           length(missing_r2),
                           paste(head(missing_r2, 5), collapse = ", "))))
        return()
      }
      scan_result(list(
        success = TRUE,
        mode = "local",
        local_dir = dir_path,
        n_samples = length(sample_names),
        sample_names = sample_names,
        is_paired = TRUE,
        has_md5 = file.exists(file.path(dir_path, "checksums.md5"))
      ))

    } else {  # SRA mode
      raw <- input$proteog_sra_accessions %||% ""
      accs <- trimws(strsplit(raw, "[,;[:space:]]+")[[1]])
      accs <- accs[nzchar(accs)]
      if (length(accs) == 0) {
        scan_result(list(success = FALSE,
                         error = "Please enter at least one SRA/ENA accession."))
        return()
      }
      if (length(accs) > 24) {
        scan_result(list(success = FALSE,
                         error = sprintf("Too many accessions (%d > 24).",
                                         length(accs))))
        return()
      }
      withProgress(message = "Verifying ENA metadata…", value = 0, {
        per_acc <- lapply(seq_along(accs), function(i) {
          incProgress(1 / length(accs), detail = accs[i])
          verify_sra_accession(accs[i])
        })
      })
      n_bad <- sum(!vapply(per_acc, function(r) isTRUE(r$success), logical(1)))
      scan_result(list(
        success = n_bad == 0,
        mode = "sra",
        accessions = accs,
        per_acc = per_acc,
        sample_names = accs,
        error = if (n_bad > 0) sprintf("%d accession(s) failed verification.", n_bad) else NULL
      ))
    }
  })

  output$proteog_scan_output <- renderUI({
    res <- scan_result()
    if (is.null(res)) {
      return(helpText("Click ", strong("Scan / Verify"),
                      " after filling in a source above."))
    }
    if (!isTRUE(res$success)) {
      return(tags$div(class = "alert alert-danger",
                      tags$strong("Error: "), res$error %||% "Unknown error."))
    }
    if (identical(res$mode, "slims")) {
      tagList(
        tags$div(class = "alert alert-success",
                 sprintf("SLIMS URL OK. %d samples (%s).",
                         res$n_samples,
                         if (isTRUE(res$is_paired)) "paired-end" else "single-end")),
        tags$details(
          tags$summary(sprintf("Sample list (%d)", length(res$sample_names))),
          tags$ul(lapply(res$sample_names, function(s) tags$li(code(s))))
        )
      )
    } else {
      # SRA mode — show per-accession metadata
      rows <- lapply(res$per_acc, function(r) {
        if (!isTRUE(r$success)) {
          return(tags$tr(tags$td(code(r$accession %||% "?")),
                         tags$td(colspan = 4,
                                 tags$span(style="color:#c62828;", r$error))))
        }
        tags$tr(
          tags$td(code(r$accession)),
          tags$td(r$scientific_name %||% "?"),
          tags$td(r$library_strategy %||% "?",
                  if (!isTRUE(r$suitable))
                    tags$span(style = "color:#c62828; font-weight:bold;",
                              " ⚠ unsuitable")),
          tags$td(r$instrument %||% "?"),
          tags$td(r$layout %||% "?")
        )
      })
      tagList(
        tags$table(class = "table table-sm",
                   tags$thead(tags$tr(
                     tags$th("Accession"), tags$th("Species"),
                     tags$th("Library"), tags$th("Instrument"), tags$th("Layout")
                   )),
                   tags$tbody(rows))
      )
    }
  })

  # ── Pre-submit warnings + button enable/disable ─────────────────────────────
  output$proteog_submit_warnings <- renderUI({
    warnings <- character()
    res <- scan_result()
    if (is.null(res) || !isTRUE(res$success)) {
      warnings <- c(warnings, "Run Scan / Verify before submitting.")
    }
    if (!nzchar(input$proteog_reference_key %||% "")) {
      warnings <- c(warnings, "Select a reference genome.")
    }
    if (!nzchar(input$proteog_project_name %||% "")) {
      warnings <- c(warnings, "Enter a project name.")
    }
    if (!nzchar(input$proteog_project_tag %||% "")) {
      warnings <- c(warnings, "Enter a project tag.")
    }
    # Species mismatch check (SRA mode only)
    if (!is.null(res) && identical(res$mode, "sra") &&
        nzchar(input$proteog_reference_key %||% "")) {
      reg <- load_reference_registry(ssh_config = proteog_ssh_config())
      ref_org <- reg[[input$proteog_reference_key]]$organism %||% ""
      accs_orgs <- vapply(res$per_acc,
        function(r) if (isTRUE(r$success)) r$scientific_name %||% "" else "",
        character(1))
      mismatched <- accs_orgs[nzchar(accs_orgs) & accs_orgs != ref_org]
      if (length(mismatched) > 0) {
        warnings <- c(warnings, sprintf(
          "Species mismatch: reference is %s; accessions report %s.",
          ref_org, paste(unique(mismatched), collapse = ", ")))
      }
    }
    if (length(warnings) == 0) return(NULL)
    tags$div(class = "alert alert-warning",
             tags$strong("Cannot submit yet:"),
             tags$ul(lapply(warnings, tags$li)))
  })

  # ── Submit handler ──────────────────────────────────────────────────────────
  observeEvent(input$proteog_submit_btn, {
    res <- scan_result()
    req(res, isTRUE(res$success))
    req(nzchar(input$proteog_reference_key %||% ""))
    req(nzchar(input$proteog_project_name %||% ""))
    req(nzchar(input$proteog_project_tag  %||% ""))

    # Resolve the rnaseq_dir + decide whether we need a download.
    #   local                  → use user-provided dir, no download, immediate submit
    #   sra/slims, files there → skip download, immediate submit (reuse existing data)
    #   sra/slims, fresh       → launch download, queue submit until status=complete
    sc <- proteog_ssh_config()
    pname <- sanitize_project_name(input$proteog_project_name)
    resolution <- tryCatch({
      if (identical(res$mode, "slims") || identical(res$mode, "sra")) {
        target_dir <- file.path(PROTEOG_RNASEQ_ROOT, pname)
        if (rnaseq_data_already_present(target_dir, sc)) {
          list(rnaseq_dir = target_dir, needs_download = FALSE,
               note = "FASTQ files already present on Hive — skipping download.")
        } else if (identical(res$mode, "slims")) {
          d <- launch_slims_download(res$url, pname, ssh_config = sc)
          list(rnaseq_dir = d$project_dir, needs_download = TRUE,
               note = "Download started — pipeline will auto-submit when files arrive.")
        } else {
          d <- launch_ena_download(
            res$accessions, pname,
            subsample_reads = if (isTRUE(input$proteog_subsample)) 5e6L else NULL,
            ssh_config = sc)
          list(rnaseq_dir = d$project_dir, needs_download = TRUE,
               note = "Download started — pipeline will auto-submit when files arrive.")
        }
      } else if (identical(res$mode, "local")) {
        list(rnaseq_dir = res$local_dir, needs_download = FALSE,
             note = "Submitting pipeline against on-cluster data.")
      } else {
        stop("Unknown source mode: ", res$mode)
      }
    }, error = function(e) {
      showNotification(sprintf("Source resolution failed: %s",
                               conditionMessage(e)),
                       type = "error", duration = 10)
      NULL
    })
    req(resolution)
    rnaseq_dir     <- resolution$rnaseq_dir
    needs_download <- resolution$needs_download

    showNotification(
      tags$div(
        tags$p(strong(resolution$note)),
        tags$p(tags$code(rnaseq_dir))),
      type = "message", duration = 10
    )

    # Build-args closure: capture everything submit_proteogenomics_build() needs
    # so we can either fire immediately (local mode) OR wait for the download
    # to finish (sra/slims modes, via the poll-and-submit observer below).
    build_args <- list(
      project_name    = pname,
      rnaseq_dir      = rnaseq_dir,
      reference_key   = input$proteog_reference_key,
      sample_names    = res$sample_names %||% character(0),
      library_type    = input$proteog_library_type,
      strand_flag     = input$proteog_strand_flag,
      project_tag     = input$proteog_project_tag,
      min_orf_len     = as.integer(input$proteog_min_orf_len %||% 100L),
      uniprot_fasta   = resolve_proteog_uniprot_path(),
      slurm_account   = "genome-center-grp",
      slurm_partition = "high",
      ssh_config      = proteog_ssh_config()
    )

    if (!needs_download) {
      # Data is already on Hive — submit immediately (covers local mode
      # AND the sra/slims-but-already-downloaded case)
      tryCatch({
        build <- do.call(submit_proteogenomics_build, build_args)
        jobs <- values$proteog_build_jobs %||% list()
        jobs[[length(jobs) + 1L]] <- list(
          project_name = build_args$project_name,
          project_dir  = build$project_dir,
          submitted_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
          jids_by_stage = build$jids_by_stage,
          methods_paragraph = build$methods_paragraph
        )
        values$proteog_build_jobs <- jobs
        showNotification(sprintf("Build submitted: %s", build_args$project_name),
                         type = "default", duration = 8)
      }, error = function(e) {
        showNotification(sprintf("Submit failed: %s", conditionMessage(e)),
                         type = "error", duration = 15)
      })
    } else {
      # sra/slims: download is running. Queue the build for the poll-and-submit
      # observer AND add a placeholder entry to proteog_build_jobs so the user
      # sees it in Active Builds immediately. When the download finishes, the
      # placeholder is upgraded with real SLURM job IDs (matched by project_dir).
      pending <- pending_build_submits()
      pending[[rnaseq_dir]] <- list(
        rnaseq_dir = rnaseq_dir,
        build_args = build_args,
        queued_at  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
      )
      pending_build_submits(pending)

      jobs <- values$proteog_build_jobs %||% list()
      jobs[[length(jobs) + 1L]] <- list(
        project_name      = pname,
        project_dir       = rnaseq_dir,
        submitted_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
        jids_by_stage     = list(),
        methods_paragraph = NULL,
        download_pending  = TRUE
      )
      values$proteog_build_jobs <- jobs

      showNotification(
        tags$div(
          tags$p(strong("Download running on Hive."),
                 " Pipeline will auto-submit when files arrive."),
          tags$p("You can close the browser; build resumes on Hive.")),
        type = "message", duration = 10)
    }
  })

  # ── Pending-build queue: waits for download to finish, then fires submit ──
  pending_build_submits <- reactiveVal(list())

  proteog_pending_poll <- reactivePoll(
    intervalMillis = 15000,
    session = session,
    checkFunc = function() {
      # Re-tick at least every 15 s while we have pending entries
      paste(names(pending_build_submits()), Sys.time(), collapse = "|")
    },
    valueFunc = function() {
      sc <- proteog_ssh_config()
      lapply(pending_build_submits(), function(p) {
        st <- tryCatch(poll_download_status(p$rnaseq_dir, ssh_config = sc),
                       error = function(e) list(state = "error",
                                                 message = conditionMessage(e)))
        list(rnaseq_dir = p$rnaseq_dir, state = st$state %||% "unknown",
             p = p)
      })
    }
  )

  observe({
    entries <- proteog_pending_poll()
    if (length(entries) == 0) return(invisible())
    pending <- pending_build_submits()
    changed <- FALSE
    for (e in entries) {
      st <- e$state
      pname <- e$p$build_args$project_name %||% basename(e$rnaseq_dir)
      if (identical(st, "complete")) {
        # Download finished — fire the build and upgrade the placeholder entry
        tryCatch({
          build <- do.call(submit_proteogenomics_build, e$p$build_args)
          jobs <- values$proteog_build_jobs %||% list()
          idx <- which(vapply(jobs, function(j)
                              identical(j$project_dir, e$rnaseq_dir),
                              logical(1)))
          replacement <- list(
            project_name      = pname,
            project_dir       = build$project_dir,
            submitted_at      = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
            jids_by_stage     = build$jids_by_stage,
            methods_paragraph = build$methods_paragraph,
            download_pending  = FALSE
          )
          if (length(idx) > 0) {
            jobs[[idx[1]]] <- replacement
          } else {
            jobs[[length(jobs) + 1L]] <- replacement
          }
          values$proteog_build_jobs <- jobs
          showNotification(sprintf("Download complete — build submitted: %s",
                                   pname),
                           type = "default", duration = 10)
        }, error = function(err) {
          showNotification(sprintf("Auto-submit failed for %s: %s",
                                   pname, conditionMessage(err)),
                           type = "error", duration = 15)
        })
        pending[[e$rnaseq_dir]] <- NULL; changed <- TRUE
      } else if (st %in% c("download_failed", "md5_failed", "error")) {
        # Mark the placeholder as failed so the user sees it; remove from queue.
        jobs <- values$proteog_build_jobs %||% list()
        idx <- which(vapply(jobs, function(j)
                            identical(j$project_dir, e$rnaseq_dir),
                            logical(1)))
        if (length(idx) > 0) {
          jobs[[idx[1]]]$download_pending <- FALSE
          jobs[[idx[1]]]$download_failed  <- st
          values$proteog_build_jobs <- jobs
        }
        showNotification(sprintf("Download failed for %s (%s) — build not submitted",
                                 pname, st),
                         type = "error", duration = 15)
        pending[[e$rnaseq_dir]] <- NULL; changed <- TRUE
      }
      # "running"/"missing"/"unknown" → keep waiting
    }
    if (changed) pending_build_submits(pending)
  })

  # ── Active builds table — polls status.json every 15 s ──────────────────────
  # Entries with download_pending=TRUE bypass the SLURM poll and render with a
  # "Downloading" stage instead.
  proteog_status_poll <- reactivePoll(
    intervalMillis = 15000,
    session = session,
    checkFunc = function() {
      # Tick on every interval regardless of file mtime — the SSH-mounted
      # case can't reliably check remote file mtime, and we want polled
      # refresh anyway.
      paste(length(values$proteog_build_jobs %||% list()), Sys.time())
    },
    valueFunc = function() {
      jobs <- values$proteog_build_jobs %||% list()
      lapply(jobs, function(j) {
        if (isTRUE(j$download_pending)) {
          return(list(
            project_name   = j$project_name,
            project_dir    = j$project_dir,
            submitted_at   = j$submitted_at,
            current_stage  = "downloading",
            stages         = list(),
            placeholder    = TRUE
          ))
        }
        if (!is.null(j$download_failed) && nzchar(as.character(j$download_failed))) {
          return(list(
            project_name   = j$project_name,
            project_dir    = j$project_dir,
            submitted_at   = j$submitted_at,
            current_stage  = sprintf("dl-%s", j$download_failed),
            stages         = list(),
            placeholder    = TRUE
          ))
        }
        tryCatch(poll_proteog_build_status(j$project_dir,
                                            ssh_config = proteog_ssh_config()),
                 error = function(e) {
                   message(sprintf("[proteog-poll] %s: %s",
                                   j$project_name %||% basename(j$project_dir %||% "?"),
                                   conditionMessage(e)))
                   list(project_name = j$project_name,
                        project_dir  = j$project_dir,
                        submitted_at = j$submitted_at,
                        current_stage = "?", stages = list(),
                        placeholder = TRUE)
                 })
      })
    }
  )

  output$proteog_active_builds_table <- renderUI({
    statuses <- proteog_status_poll()
    if (length(statuses) == 0) {
      return(helpText("No builds submitted in this session yet."))
    }
    rows <- lapply(seq_along(statuses), function(i) {
      st <- statuses[[i]]
      if (is.null(st)) return(NULL)
      current <- st$current_stage %||% "?"
      badge_color <- switch(current,
        "complete"    = "#27ae60",
        "failed"      = "#c0392b",
        "downloading" = "#3498db",
        "#f39c12")
      done <- if (length(st$stages))
        sum(vapply(st$stages, function(s) identical(s$status, "complete"), logical(1)))
      else 0L
      total <- length(st$stages)
      progress_txt <- if (total == 0L && identical(current, "downloading"))
        "downloading…"
      else sprintf("%d / %d", done, total)

      # Per-row Assemble button: shown when the assemble stage is the current
      # one AND it has no job_id yet (i.e., legacy builds that finished the
      # SLURM chain pre-auto-assemble).
      assemble_btn <- NULL
      asm_stage <- NULL
      if (length(st$stages) > 0) {
        # Find a stage named "assemble"
        for (s in st$stages) {
          if (identical(s$stage, "assemble")) { asm_stage <- s; break }
        }
      }
      needs_assemble <- !is.null(asm_stage) &&
        is_scalar_char_safe(.empty_or_str(asm_stage$status) %||% "") &&
        identical(.empty_or_str(asm_stage$status), "unknown") &&
        !is_scalar_char_safe(.empty_or_str(asm_stage$job_id))
      if (needs_assemble) {
        pdir_attr <- htmltools::htmlEscape(st$project_dir %||% "", attribute = TRUE)
        pname_attr <- htmltools::htmlEscape(st$project_name %||% "", attribute = TRUE)
        assemble_btn <- HTML(sprintf(
          '<button class="btn btn-warning btn-sm" onclick="Shiny.setInputValue(\'proteog_assemble_btn\', {dir: \'%s\', name: \'%s\', _ts: Date.now()}, {priority:\'event\'})"><i class="fa fa-cubes"></i> Assemble</button>',
          pdir_attr, pname_attr))
      }

      tags$tr(
        tags$td(st$project_name %||% "?"),
        tags$td(tags$span(style = sprintf("background:%s; color:white; padding:2px 8px; border-radius:4px;",
                                           badge_color),
                          current)),
        tags$td(progress_txt),
        tags$td(st$submitted_at %||% "?"),
        tags$td(code(basename(st$project_dir %||% ""))),
        tags$td(assemble_btn)
      )
    })
    tags$table(class = "table table-sm",
               tags$thead(tags$tr(
                 tags$th("Project"), tags$th("Stage"),
                 tags$th("Progress"), tags$th("Submitted"),
                 tags$th("Dir"), tags$th("Action")
               )),
               tags$tbody(rows))
  })

  # ── Per-row Assemble: opens modal with the same UniProt source dropdown
  # used in Step 4, lets user pick a UniProt path / search / NCBI, then
  # submits just the assemble.sbatch on Hive for that specific project.
  proteog_assemble_target <- reactiveVal(NULL)

  observeEvent(input$proteog_assemble_btn, {
    payload <- input$proteog_assemble_btn
    if (is.null(payload) || is.null(payload$dir)) return()
    proteog_assemble_target(list(project_dir = as.character(payload$dir),
                                  project_name = as.character(payload$name %||%
                                                              basename(payload$dir))))
    showModal(modalDialog(
      title = tagList(icon("cubes"), sprintf(" Assemble FASTA for %s",
                                              payload$name)),
      size = "l", easyClose = TRUE,
      tags$p("Generate the final proteogenomics FASTA for this build. ",
             "Pick whether to combine the predicted ORFs with a UniProt or NCBI proteome."),
      selectInput("proteog_assemble_uniprot_source", "UniProt FASTA",
                  choices = c("None — predicted ORFs only" = "none",
                              "Download from UniProt"      = "uniprot",
                              "Download from NCBI"         = "ncbi",
                              "Enter path on Hive"         = "path"),
                  selected = "none", width = "100%"),
      conditionalPanel("input.proteog_assemble_uniprot_source == 'uniprot'",
        actionButton("proteog_open_uniprot_modal", "Search UniProt",
                     class = "btn-info btn-sm", icon = icon("search")),
        uiOutput("proteog_uniprot_selected_summary")),
      conditionalPanel("input.proteog_assemble_uniprot_source == 'ncbi'",
        actionButton("proteog_open_ncbi_modal", "Search NCBI",
                     class = "btn-success btn-sm", icon = icon("search")),
        uiOutput("proteog_ncbi_selected_summary")),
      conditionalPanel("input.proteog_assemble_uniprot_source == 'path'",
        textInput("proteog_uniprot_fasta_path_modal", label = NULL,
                  placeholder = "/quobyte/proteomics-grp/de-limp/databases/uniprot/UP000005640.fasta",
                  width = "100%")),
      footer = tagList(modalButton("Cancel"),
        actionButton("proteog_submit_assemble_btn", "Submit Assemble",
                     class = "btn-warning", icon = icon("play")))
    ))
  })

  observeEvent(input$proteog_submit_assemble_btn, {
    target <- proteog_assemble_target()
    if (is.null(target)) return()
    sc <- proteog_ssh_config()
    if (is.null(sc)) {
      showNotification("Connect to Hive first.", type = "warning", duration = 5)
      return()
    }
    # Resolve UniProt FASTA path from the modal's dropdown
    src <- input$proteog_assemble_uniprot_source %||% "none"
    uniprot_path <- if (identical(src, "none")) ""
      else if (identical(src, "path"))
        input$proteog_uniprot_fasta_path_modal %||% ""
      else if (identical(src, "uniprot"))
        proteog_uniprot_state$hive_path %||% ""
      else if (identical(src, "ncbi"))
        proteog_ncbi_state$hive_path %||% ""
      else ""

    tryCatch({
      jid <- submit_assemble_only(
        project_dir   = target$project_dir,
        project_name  = target$project_name,
        uniprot_fasta = uniprot_path,
        ssh_config    = sc)
      showNotification(sprintf("Assemble job submitted: %s (SLURM %s)",
                               target$project_name, jid),
                       type = "default", duration = 8)
      removeModal()
    }, error = function(e) {
      showNotification(sprintf("Assemble submit failed: %s",
                               conditionMessage(e)),
                       type = "error", duration = 12)
    })
  })

  # ── Info modals ("?" buttons in each card header) ───────────────────────────
  # Pattern matches existing DE-LIMP info modals (CLAUDE.md):
  # actionButton(..._info_btn, icon("question-circle")) + observeEvent + showModal.

  # Workflow overview — aimed at a proteomics user who doesn't know
  # genomics / RNA-seq terminology. Long-form, sectioned.
  observeEvent(input$proteog_explain_workflow_btn, {
    showModal(modalDialog(
      title = tagList(icon("dna"), " Proteogenomics for proteomics people"),
      size = "xl", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.92em; line-height: 1.65; max-height: 70vh; overflow-y: auto; padding-right: 8px;",

        tags$h5("The 30-second pitch"),
        tags$p("Your DIA-NN search assumes the proteome lives in a FASTA file — usually UniProt's canonical set for the organism. That's great when your sample expresses ",
               "the same proteins everyone else does. But cancer cells, knock-outs, non-model organisms, tissue-specific splicing, and disease samples often produce ",
               tags$strong("proteins that don't exist in UniProt"), " — alternative splice forms, intergenic transcripts, novel ORFs. Your peptides for those proteins ",
               "land in DIA-NN as ", tags$em("unidentified"), " because the FASTA never had a matching sequence."),
        tags$p("This tool fixes that. You give it the RNA-seq from ", tags$em("the same samples"), " you're going to run mass-spec on. It figures out what your cells are ",
               "actually transcribing, predicts the proteins those transcripts encode, and merges them with the canonical UniProt FASTA. Now DIA-NN sees both the standard ",
               "proteome AND your sample-specific extras."),

        tags$hr(),
        tags$h5("What is RNA-seq, in proteomics terms?"),
        tags$p("RNA-seq is shotgun sequencing of mRNA. Think of it like bottom-up proteomics:"),
        tags$ul(
          tags$li(tags$strong("Sample: "), "Bulk RNA from your cells/tissue (analogous to your protein lysate)."),
          tags$li(tags$strong("Fragmentation: "), "RNA → cDNA → short fragments (~150–300 bp; analogous to tryptic digestion to peptides)."),
          tags$li(tags$strong("Detector: "), "Illumina sequencer reads 30M+ fragments per sample (analogous to a mass spectrometer producing MS2 spectra)."),
          tags$li(tags$strong("Output: "), "A pair of FASTQ files per sample — ", tags$code("sample_R1.fastq.gz"), " + ", tags$code("sample_R2.fastq.gz"),
                  " — each holding millions of 150-base reads with quality scores. (Analogous to the .raw/.d file from your instrument.)")
        ),
        tags$p("The expensive part is the sequencing run (~$50–200 per sample at a core facility). The compute is free once you have the FASTQs."),

        tags$hr(),
        tags$h5("What is the \"reference genome\"?"),
        tags$p("The genome is the full DNA sequence for the organism — chromosomes 1, 2, 3, ... laid out end-to-end. A canonical \"reference\" assembly (e.g. ",
               tags$code("GRCh38"), " for human, ", tags$code("GRCm39"), " for mouse) is the version everyone in the field agrees to use. We also need an ",
               tags$strong("annotation file (GTF)"), " — a separate text file that says \"chromosome 1 positions 1000–2500 are a gene called ACTB; positions 1000–1500 ",
               "are exon 1, 1700–2000 are exon 2, ...\" — basically the coordinate map of every known gene."),
        tags$p("DE-LIMP keeps pre-built references (genome FASTA + GTF + STAR alignment index + rRNA filter index) on Hive shared storage. You pick one from the ",
               "Reference dropdown — pick the species you're working in."),

        tags$hr(),
        tags$h5("What the pipeline does, stage by stage"),
        tags$p("It's a SLURM dependency chain — 11 stages, each waits for the previous one. ~3–6 hours wall time for a typical 12-sample run."),
        tags$ol(
          tags$li(tags$strong("fastp"), " — adapter trimming + quality filter on the raw FASTQ. Removes the technical sequencing artifacts. Detects read length to ",
                  "pick the right STAR alignment-stringency tier later (150 bp reads need stricter QC than 50 bp reads)."),
          tags$li(tags$strong("rRNA filter (bowtie2)"), " — most RNA in a cell is ribosomal RNA (rRNA), not the messenger RNA you care about. We map every read against ",
                  "an rRNA-only sequence library and keep ONLY the reads that don't match — those are your mRNA. (For proteomics analogy: it's like depleting albumin ",
                  "from serum before LC-MS.)"),
          tags$li(tags$strong("STAR"), " — spliced alignment of reads against the reference genome. mRNA reads usually span ", tags$em("exon junctions"),
                  " (gene segments stitched together with introns removed), so the aligner has to handle non-contiguous matches. This is the heavy step: ",
                  "~48 GB RAM, 16 cores, 10 min/sample. Output: BAM file with every read's genomic coordinates."),
          tags$li(tags$strong("QC gate"), " — checks STAR's uniquely-mapped rate against the read-length tier threshold (25% for short reads, 60% for long). ",
                  tags$strong("Halts the chain"), " if below — never produces a half-baked bad FASTA. If you fail QC, see the report and either re-check your library prep or relax the strand/library settings."),
          tags$li(tags$strong("stringtie"), " — per-sample transcript assembly. From the BAM, it figures out which reads belong to which transcript ",
                  "(\"this set of reads forms one mRNA, that set forms a different splice variant\"). Uses the reference GTF as a guide."),
          tags$li(tags$strong("merge"), " — combines all per-sample GTFs into a unified ", tags$em("non-redundant"), " transcript catalog for the study. ",
                  "Different samples may have caught slightly different transcripts; merge is how we agree on a single coordinate system for downstream steps."),
          tags$li(tags$strong("gffcompare"), " — classifies every transcript in the merged GTF vs the reference annotation. This is the key step that distinguishes:",
                  tags$ul(
                    tags$li(tags$code("REF"), " — already in the reference (canonical UniProt protein)"),
                    tags$li(tags$code("NOVEL_ISOFORM"), " — alternative splicing of a known gene (different exon combination, possibly different protein N/C-terminus)"),
                    tags$li(tags$code("NOVEL_GENE"), " — transcript in an intergenic region (could be a new protein-coding gene, or a long non-coding RNA — TransDecoder decides next)")
                  )),
          tags$li(tags$strong("gffread"), " — extracts the actual nucleotide sequence for every merged transcript from the genome FASTA + GTF coordinates. ",
                  "Now we have a FASTA of mRNA-like sequences, not protein sequences yet."),
          tags$li(tags$strong("TransDecoder"), " — predicts ORFs (open reading frames). For every transcript, finds the longest ATG-to-stop run and decides whether it ",
                  "looks like a real protein-coding region (based on length, Markov chain stats, and an optional DIAMOND BLAST hit against SwissProt). ",
                  tags$em("This is where transcripts become proteins."),
                  " Configurable minimum ORF length (default 100 aa)."),
          tags$li(tags$strong("rewrite headers"), " — the predicted ORFs come out with cryptic ", tags$code("TRINITY_DN123_c0_g1_i1.p1"),
                  "-style names. We rewrite them as ", tags$code("sp|<id>|<sym>_<TAG>"), " (mimicking UniProt's format) so DIA-NN's gene-name column works correctly. ",
                  "The ", tags$code("<TAG>"), " (your \"project tag\", e.g. ", tags$code("MOUSELIVER"),
                  ") tells you at a glance which proteins came from your custom FASTA vs UniProt."),
          tags$li(tags$strong("assemble"), " — final concatenation: predicted ORFs + (optional) UniProt FASTA → single combined FASTA, written to ",
                  tags$code("/quobyte/.../databases/proteogenomics/"), " on Hive. This is what DIA-NN searches against.")
        ),

        tags$hr(),
        tags$h5("What you get"),
        tags$ul(
          tags$li("A FASTA file named ", tags$code("<your_project>_proteogenomics_<YYYY_MM>.fasta"),
                  " on Hive, with REF / NOVEL_ISOFORM / NOVEL_GENE annotations baked into the FASTA headers."),
          tags$li("An entry in the Run Search → ", tags$strong("Proteogenomics DBs"), " dropdown on the main DE-LIMP page, ",
                  "with full provenance (samples, reference, methods paragraph)."),
          tags$li("Traceability: every protein in your DIA-NN results that ends with ", tags$code("_<TAG>"), " is one we discovered in YOUR data, not from UniProt.")
        ),

        tags$hr(),
        tags$h5("When this is worth it (and when it isn't)"),
        tags$ul(
          tags$li(tags$strong("Worth it: "), "cancer/disease samples with potential novel splicing; non-model organism (no good UniProt entry); ",
                  "knock-out/knock-in cell lines; tissue with known alternative splicing (brain, testis); discovering chimeric proteins; novel ORFs from long non-coding RNAs."),
          tags$li(tags$strong("Skip it: "), "you're running stock HEK293/HeLa with a well-curated UniProt proteome (you'll find ~5–20 novel proteins, mostly junk); ",
                  "you don't have matched RNA-seq from the same samples (mismatched RNA + protein from different conditions = false positives); your goal is just quantification, ",
                  "not discovery.")
        ),
        tags$h6("A quick reality check:"),
        tags$ul(
          tags$li("Cost: ~$50–200/sample sequencing + 3–6 h compute (free on Hive)."),
          tags$li("Sample volume: only 12 samples? Pool RNA from biological reps OR sequence one rep per condition. RNA-seq is per-condition, not per-sample-injection."),
          tags$li("False positives: NOVEL_GENE calls without DIAMOND BLAST support are statistically uncertain. Trust REF + NOVEL_ISOFORM most.")
        )
      )
    ))
  })

  observeEvent(input$proteog_step1_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("upload"), " Step 1 — RNA-seq Source"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Pick where your RNA-seq data lives"),
        tags$ul(
          tags$li(strong("DNA Tech Core (SLIMS URL):"),
                  " For UC Davis users whose sequencing was done at the DNA Technologies Core. ",
                  "Your delivery email contains a URL like ",
                  tags$code("http://slimsdata.genomecenter.ucdavis.edu/Data/<id>/Unaligned/"),
                  ". DE-LIMP will download all R1/R2 FASTQ files automatically + verify md5 checksums."),
          tags$li(strong("Public archive (SRA/ENA):"),
                  " For re-analyzing published datasets. Paste one or more accession IDs ",
                  "(e.g., ", tags$code("SRR1303776"), "). DE-LIMP queries ENA for metadata ",
                  "(species, library type) before download to catch mistakes like wrong-species ",
                  "accessions. Use the ", strong("subsample"), " checkbox for a fast test on the ",
                  "first 5 million read pairs."),
          tags$li(strong("Folder on the cluster:"),
                  " For data you already have on Hive — from another sequencing core, a ",
                  "previous download, or a collaborator. Paste the full path to a directory ",
                  "containing paired ", tags$code("<sample>_R1.fastq.gz"), " / ",
                  tags$code("<sample>_R2.fastq.gz"), " files.")
        ),
        tags$h6("Then click Scan / Verify"),
        p("The pipeline does basic checks before you submit anything heavy — wrong-species ",
          "verification, missing R2 detection, etc. Catching this here saves ~30 minutes of ",
          "wasted compute later."),
        tags$h6("Unsuitable libraries the pipeline will refuse"),
        p("Tag-Seq, miRNA-Seq, and other libraries that don't cover full transcripts can't ",
          "be used for proteogenomics novel-ORF discovery. If ENA reports one of these ",
          "for an accession, the Submit button stays disabled.")
      )
    ))
  })

  observeEvent(input$proteog_step2_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("clipboard-check"), " Step 2 — Sample Verification"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What's being verified"),
        tags$ul(
          tags$li(strong("Sample count + pairing:"),
                  " confirms every sample has both R1 and R2 (paired-end is required)."),
          tags$li(strong("Species (SRA only):"),
                  " ENA metadata is queried to confirm the organism. A common mistake is ",
                  "assuming an accession is one species when it's actually another — the ",
                  "validation that built DE-LIMP caught this exact bug with SRR1303776/77, ",
                  "claimed to be K562 human but actually mouse."),
          tags$li(strong("Library strategy (SRA only):"),
                  " RNA-Seq, polyA, total RNA — these are all OK. Tag-Seq, miRNA-Seq, ",
                  "Ribo-Seq, CLIP-Seq are flagged as unsuitable.")
        ),
        tags$h6("If verification fails"),
        p("The Submit button (Step 5) will stay disabled until everything passes. The ",
          "error message tells you what to fix. For species mismatches, change your reference ",
          "selection in Step 3 to match the actual organism.")
      )
    ))
  })

  observeEvent(input$proteog_step3_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("book-open"), " Step 3 — Reference Genome"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Pick the reference for your samples' organism"),
        p("The dropdown is populated from a curated registry of pre-staged references ",
          "(at ", tags$code("/quobyte/proteomics-grp/de-limp/references/registry.json"),
          "). Each reference includes the genome FASTA, the STAR index, the GTF annotation, ",
          "and the organism-specific rRNA filter sequences."),
        tags$h6("Currently available"),
        tags$ul(
          tags$li(strong("Mus musculus (mm39/GRCm39):"),
                  " GENCODE vM38 basic annotation. Full STAR index pre-built. Includes ",
                  "mouse-specific rRNA filter."),
          tags$li(strong("Homo sapiens (hg38/GRCh38.p14):"),
                  " RefSeq annotation. Full STAR index pre-built. Includes human-specific ",
                  "rRNA filter.")
        ),
        tags$h6("Adding new references"),
        p("Need a non-mouse, non-human reference (e.g., zebrafish, fly, plant)? Talk to ",
          "the Proteomics Core — staging a new reference is a one-time ~12 hour admin job ",
          "(download genome + GTF, build STAR index, build rRNA bowtie2 index, register).")
      )
    ))
  })

  observeEvent(input$proteog_step4_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("sliders"), " Step 4 — Parameters"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("Library type"),
        tags$ul(
          tags$li(strong("polyA mRNA-Seq:"), " standard mRNA-seq, polyA enriched. Most common."),
          tags$li(strong("Total RNA + rRNA depletion:"), " captures non-coding RNAs too. ",
                  "Preferred for proteogenomics because it captures more lincRNA and small ORF transcripts."),
          tags$li(strong("Stranded RNA-Seq:"), " directional library, tells the aligner which ",
                  "strand was transcribed.")
        ),
        tags$h6("Strand"),
        p("If you know your library prep kit, pick the right one:"),
        tags$ul(
          tags$li(strong("--rf (reverse stranded):"),
                  " TruSeq Stranded, Illumina Stranded mRNA. Most common in 2024+."),
          tags$li(strong("--fr (forward stranded):"),
                  " older prep kits, Lexogen QuantSeq 3' FWD."),
          tags$li(strong("Unstranded:"), " older non-stranded protocols (rare in 2024+).")
        ),
        p("If you're not sure, ", strong("--rf"), " is the right guess. The aligner is forgiving ",
          "of mis-specified strand; you'll lose ~10% of transcripts in the worst case."),
        tags$h6("Project name"),
        p("Short identifier for this build, e.g., ", tags$code("mouse_liver_pilot_2026_05"), ". ",
          "Used as the output directory name. Allowed: letters, numbers, dots, dashes, underscores."),
        tags$h6("Project tag"),
        p("An UPPERCASE tag suffixed to every predicted-protein symbol so they're traceable ",
          "back to which build produced them, e.g., ", tags$code("Gnai3_MOUSELIVER"),
          " instead of ", tags$code("Gnai3"),
          ". Keep it short — it shows up in every header."),
        tags$h6("Minimum ORF length"),
        p("TransDecoder won't predict ORFs shorter than this. Default 100 aa is the field ",
          "standard for novel-ORF discovery. Lowering to 50-70 aa picks up small ORFs (sORFs, ",
          "uORFs) at the cost of more false positives.")
      )
    ))
  })

  observeEvent(input$proteog_step5_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("rocket"), " Step 5 — Submit"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("What happens when you click Build"),
        tags$ol(
          tags$li(strong("Download (if SLIMS / SRA mode):"),
                  " a background process on the login node fetches your FASTQ files. ",
                  "Status appears under Active Builds below."),
          tags$li(strong("SLURM dependency chain (10 jobs):"),
                  " once data is on disk, the pipeline submits as 10 sbatch jobs with ",
                  tags$code("--dependency=afterok"), " chaining: fastp → bowtie2 rRNA filter → ",
                  "STAR → QC gate → stringtie → merge → gffcompare → gffread → ",
                  "TransDecoder → header rewrite."),
          tags$li(strong("Quality gates fire automatically:"),
                  " if the STAR uniquely-mapped rate is below threshold (25% for short reads, ",
                  "60% for ≥130bp reads), the chain halts and surfaces the cause to you. ",
                  "Doesn't silently produce bad data."),
          tags$li(strong("Final FASTA assembly:"),
                  " predicted ORFs + UniProt reference + contaminants → one FASTA, deduplicated."),
          tags$li(strong("Registry:"),
                  " the resulting FASTA appears in the Run Search FASTA dropdown with a ",
                  HTML("&#x1F9EC;"), " tag and composition breakdown.")
        ),
        tags$h6("Wall time"),
        p("Roughly 3-6 hours for 12 samples × 30M PE150 reads on the ", tags$code("high"),
          " partition. You can close the browser; the chain continues on Hive. ",
          "Check back via the Active Builds table below or the Output tab when complete."),
        tags$h6("What if it fails midway"),
        p("Every stage writes its own log under ", tags$code("logs/"),
          " in the project directory. The Active Builds table shows which stage failed. ",
          "Most failures are recoverable by re-running just the failed step — see the ",
          strong("Output"), " tab for resume options.")
      )
    ))
  })

  observeEvent(input$proteog_builds_info_btn, {
    showModal(modalDialog(
      title = tagList(icon("list-check"), " Pipeline Stages"),
      size = "l", easyClose = TRUE, footer = modalButton("Close"),
      div(style = "font-size: 0.9em; line-height: 1.7;",
        tags$h6("The 10 stages, in order"),
        tags$ol(
          tags$li(strong("fastp:"), " adapter trimming + quality filter. Detects read length to drive STAR threshold tier selection. ~5 min/sample, parallel array."),
          tags$li(strong("rrna_filter:"), " bowtie2 against organism-specific rRNA sequences. Reads NOT matching rRNA proceed to STAR. ~5 min/sample, parallel array."),
          tags$li(strong("star:"), " spliced alignment to the reference genome. ~10 min/sample on 16 cores. Highest memory use (~48 GB). Parallel array."),
          tags$li(strong("qc_gate:"), " checks STAR uniquely-mapped rate against tier threshold (25% / 45% / 60% based on read length). HALTS the chain if below — never produces a partial bad FASTA."),
          tags$li(strong("stringtie:"), " per-sample transcript assembly using STAR BAMs + reference GTF for guidance. ~1 min/sample, parallel array."),
          tags$li(strong("merge:"), " combines per-sample GTFs into a unified transcript model. ~1 min."),
          tags$li(strong("gffcompare:"), " classifies merged transcripts vs the reference annotation. This is what distinguishes REF (annotated) from NOVEL_ISOFORM (alternative splicing) from NOVEL_GENE (intergenic). ~30 sec."),
          tags$li(strong("gffread:"), " extracts transcript-level FASTA sequences from the genome FASTA + merged GTF. ~30 sec."),
          tags$li(strong("transdecoder:"), " predicts ORFs from transcripts. Longest single stage (~20 min) because it runs TransDecoder.LongOrfs + TransDecoder.Predict with start-codon refinement."),
          tags$li(strong("rewrite:"), " converts TransDecoder's idiosyncratic headers into the DE-LIMP ",
                  tags$code("sp|ID|SYM_TAG source=... ORF_type=..."),
                  " format. Fails non-zero if any UNPARSED entries are produced. ~1 min.")
        ),
        tags$h6("Stage states in the table"),
        tags$ul(
          tags$li(tags$span(style="background:#f39c12;color:white;padding:2px 8px;border-radius:4px;","pending"), " — submitted to SLURM, waiting for a node"),
          tags$li(tags$span(style="background:#f39c12;color:white;padding:2px 8px;border-radius:4px;","running"), " — actively executing"),
          tags$li(tags$span(style="background:#27ae60;color:white;padding:2px 8px;border-radius:4px;","complete"), " — succeeded; downstream stages will run"),
          tags$li(tags$span(style="background:#c0392b;color:white;padding:2px 8px;border-radius:4px;","failed"), " — non-zero exit; downstream stages will NOT run (afterok dependency holds them)")
        )
      )
    ))
  })

  invisible(NULL)
}
