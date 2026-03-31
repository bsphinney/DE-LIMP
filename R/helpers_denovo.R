# helpers_denovo.R — Pure helper functions for Cascadia De Novo Sequencing Integration
# No Shiny reactivity. All functions are testable standalone.
# Depends on: data.table, dplyr, stringr
# Remote operations use ssh_exec() from helpers_search.R

# =============================================================================
# SSL Parsing
# =============================================================================

# Parse Cascadia SSL output files
# ssl_paths: character vector of paths to .ssl files
# score_threshold: minimum confidence score to retain (default 0.8)
# Returns data.table with columns: file, scan, charge, sequence, score,
#   retention_time, source_file, seq_stripped, seq_norm
parse_cascadia_ssl <- function(ssl_paths, score_threshold = 0.8) {
  ssl_paths <- ssl_paths[file.exists(ssl_paths)]
  if (length(ssl_paths) == 0) {
    return(data.table::data.table(
      file = character(), scan = integer(), charge = integer(),
      sequence = character(), score = numeric(), retention_time = numeric(),
      source_file = character(), seq_stripped = character(), seq_norm = character()
    ))
  }

  ssl <- data.table::rbindlist(lapply(ssl_paths, function(path) {
    dt <- data.table::fread(path, sep = "\t", header = TRUE)
    # SSL headers use hyphens — normalize to underscores
    names(dt) <- gsub("-", "_", names(dt))
    dt$source_file <- basename(path)
    dt
  }), fill = TRUE)

  if (nrow(ssl) == 0) return(ssl)

  # Filter by score threshold
  if ("score" %in% names(ssl)) {
    ssl <- ssl[ssl$score >= score_threshold, ]
  }

  if (nrow(ssl) == 0) return(ssl)

  # Strip modifications for matching (e.g. SEQM[+15.99]ENCE -> SEQMENCE)
  ssl$seq_stripped <- gsub("\\[.*?\\]", "", ssl$sequence)

  # Normalize I/L ambiguity (de novo can't distinguish isoleucine/leucine)
  ssl$seq_norm <- gsub("I", "L", ssl$seq_stripped)

  # Convert RT to minutes if stored in seconds (Cascadia outputs seconds)
  if ("retention_time" %in% names(ssl)) {
    max_rt <- max(ssl$retention_time, na.rm = TRUE)
    if (is.finite(max_rt) && max_rt > 200) {
      ssl$retention_time <- ssl$retention_time / 60
    }
  }

  ssl
}


# =============================================================================
# Cross-referencing with DIA-NN
# =============================================================================

# Cross-reference Cascadia de novo results against DIA-NN peptides
# ssl_data: data.table from parse_cascadia_ssl()
# diann_report: data.table/data.frame with Stripped.Sequence and Protein.Group columns
# Returns list with $confirmed, $novel, $protein_summary
classify_denovo_peptides <- function(ssl_data, diann_report) {
  if (is.null(ssl_data) || nrow(ssl_data) == 0) {
    return(list(
      confirmed = data.table::data.table(),
      novel = data.table::data.table(),
      protein_summary = data.table::data.table(
        Protein.Group = character(),
        n_denovo_confirmed = integer(),
        denovo_max_score = numeric()
      )
    ))
  }

  # Build normalized DIA-NN peptide set
  diann_peptides <- unique(diann_report$Stripped.Sequence)
  diann_norm <- unique(gsub("I", "L", diann_peptides))

  # Classify each de novo peptide
  ssl_data$match_type <- ifelse(
    ssl_data$seq_norm %in% diann_norm,
    "confirmed",
    "novel"
  )

  # Build peptide-to-protein map from DIA-NN (I/L normalized)
  pep_to_protein <- dplyr::distinct(
    diann_report, Stripped.Sequence, Protein.Group
  )
  pep_to_protein$seq_norm <- gsub("I", "L", pep_to_protein$Stripped.Sequence)

  confirmed <- ssl_data[ssl_data$match_type == "confirmed", ]
  novel <- ssl_data[ssl_data$match_type == "novel", ]

  # Join confirmed peptides to protein groups
  if (nrow(confirmed) > 0) {
    confirmed <- merge(
      confirmed,
      pep_to_protein[, c("seq_norm", "Protein.Group")],
      by = "seq_norm", all.x = TRUE, allow.cartesian = TRUE
    )
  }

  # Per-protein summary of de novo confirmation
  if (nrow(confirmed) > 0 && "Protein.Group" %in% names(confirmed)) {
    protein_summary <- confirmed |>
      dplyr::filter(!is.na(Protein.Group)) |>
      dplyr::group_by(Protein.Group) |>
      dplyr::summarise(
        n_denovo_confirmed = dplyr::n_distinct(seq_norm),
        denovo_max_score = max(score, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    protein_summary <- data.table::data.table(
      Protein.Group = character(),
      n_denovo_confirmed = integer(),
      denovo_max_score = numeric()
    )
  }

  list(
    confirmed = confirmed,
    novel = novel,
    protein_summary = protein_summary
  )
}


# =============================================================================
# DIAMOND BLAST for Novel Peptides (runs on HPC via SSH)
# =============================================================================

# Map novel de novo peptides to proteins via DIAMOND BLAST on HPC
# novel_peptides: character vector of novel stripped sequences
# fasta_path: remote path to reference FASTA (same one used for DIA-NN search)
# ssh_config: list(host, user, port, key_path) for ssh_exec()
# sbatch_path: full path to diamond binary on HPC (or just "diamond" if in PATH)
# output_dir: remote directory for DIAMOND results
# diamond_db: remote path to pre-built DIAMOND database (auto-built if NULL)
# min_identity: minimum % identity to report a match (default 90)
# threads: number of CPU threads (default 4)
# Returns data.frame of BLAST hits
run_diamond_blast <- function(
  novel_peptides,
  fasta_path,
  ssh_config,
  sbatch_path = NULL,
  output_dir,
  diamond_db = NULL,
  min_identity = 90,
  threads = 4
) {
  if (length(novel_peptides) == 0) {
    return(data.frame(
      query_idx = character(), subject = character(), identity = numeric(),
      length = integer(), qlen = integer(), slen = integer(),
      evalue = numeric(), bitscore = numeric(),
      peptide_sequence = character(), protein = character(),
      stringsAsFactors = FALSE
    ))
  }

  diamond_bin <- sbatch_path %||% "diamond"
  query_fasta <- sprintf("%s/novel_denovo_queries.fasta", output_dir)
  blast_out <- sprintf("%s/novel_denovo_blast.tsv", output_dir)

  # Write query FASTA — build locally then upload via heredoc
  query_lines <- paste0(">", seq_along(novel_peptides), "\n", novel_peptides)
  query_content <- paste(query_lines, collapse = "\n")

  # Create output dir and write query FASTA on remote
  write_cmd <- sprintf(
    'mkdir -p "%s" && cat > "%s" << \'FASTA_EOF\'\n%s\nFASTA_EOF',
    output_dir, query_fasta, query_content
  )
  res <- ssh_exec(ssh_config, write_cmd, timeout = 30)
  if ((res$status %||% 0L) != 0L) {
    warning("Failed to write query FASTA to remote: ", paste(res$stdout, collapse = "\n"))
    return(data.frame())
  }

  # Build DIAMOND database if not provided
  if (is.null(diamond_db) || !nzchar(diamond_db)) {
    diamond_db <- sprintf("%s/ref_diamond.dmnd", output_dir)
    makedb_cmd <- sprintf(
      '%s makedb --in "%s" --db "%s" --threads %d --quiet',
      diamond_bin, fasta_path, diamond_db, threads
    )
    res <- ssh_exec(ssh_config, makedb_cmd, login_shell = TRUE, timeout = 300)
    if ((res$status %||% 0L) != 0L) {
      warning("DIAMOND makedb failed: ", paste(res$stdout, collapse = "\n"))
      return(data.frame())
    }
  }

  # Run DIAMOND blastp
  blastp_cmd <- sprintf(
    paste0(
      '%s blastp',
      ' --query "%s"',
      ' --db "%s"',
      ' --out "%s"',
      ' --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore',
      ' --id %d',
      ' --threads %d',
      ' --sensitive',
      ' --quiet'
    ),
    diamond_bin, query_fasta, diamond_db, blast_out,
    as.integer(min_identity), as.integer(threads)
  )
  res <- ssh_exec(ssh_config, blastp_cmd, login_shell = TRUE, timeout = 300)
  if ((res$status %||% 0L) != 0L) {
    warning("DIAMOND blastp failed: ", paste(res$stdout, collapse = "\n"))
    return(data.frame())
  }

  # Download results via SSH cat (small file)
  cat_cmd <- sprintf('cat "%s" 2>/dev/null', blast_out)
  res <- ssh_exec(ssh_config, cat_cmd, timeout = 30)
  if ((res$status %||% 0L) != 0L || length(res$stdout) == 0 ||
      all(!nzchar(trimws(res$stdout)))) {
    message("No DIAMOND hits found for novel de novo peptides")
    return(data.frame())
  }

  # Parse tab-delimited BLAST output from SSH stdout
  hits <- tryCatch({
    con <- textConnection(paste(res$stdout, collapse = "\n"))
    on.exit(close(con), add = TRUE)
    data.table::fread(text = paste(res$stdout, collapse = "\n"), header = FALSE)
  }, error = function(e) {
    warning("Failed to parse DIAMOND output: ", conditionMessage(e))
    return(data.frame())
  })

  if (nrow(hits) == 0) return(data.frame())

  names(hits) <- c("query_idx", "subject", "identity", "length",
                    "qlen", "slen", "evalue", "bitscore")

  # Map back to peptide sequences
  idx <- as.integer(hits$query_idx)
  hits$peptide_sequence <- ifelse(
    idx >= 1 & idx <= length(novel_peptides),
    novel_peptides[idx],
    NA_character_
  )

  # Extract clean protein accession from subject (handles sp|ACC|NAME format)
  hits$protein <- stringr::str_extract(hits$subject, "(?<=\\|)[^|]+(?=\\|)")
  no_match <- is.na(hits$protein)
  if (any(no_match)) {
    hits$protein[no_match] <- hits$subject[no_match]
  }

  as.data.frame(hits)
}


# =============================================================================
# Cascadia SLURM Job Script Generation
# =============================================================================

# Generate sbatch script for Cascadia GPU job
# Returns the script content as a character string (not written to file)
# Follows the same pattern as generate_sbatch_script() in helpers_search.R
generate_cascadia_sbatch <- function(
  analysis_name,
  raw_files,
  output_dir,
  cascadia_env = "cascadia",
  model_ckpt,
  gpu_partition = "gpu",
  gpu_account = "genome-center-grp",
  gpu_qos = NULL,
  min_score = 0.8,
  batch_size = 64,
  cpus = 8,
  mem_gb = 32,
  time_hours = 4,
  gpu_type = NULL,
  requeue = FALSE
) {
  # Sanitize analysis name for SLURM job name
  safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", analysis_name)

  # Detect file type
  file_exts <- unique(tolower(tools::file_ext(raw_files)))
  if (!all(file_exts %in% c("d", "mzml"))) {
    warning("Cascadia integration currently supports .d and .mzML files only. ",
            "Found: ", paste(file_exts, collapse = ", "))
  }

  denovo_dir <- sprintf("%s/denovo", output_dir)

  # GPU resource string

  gres <- if (!is.null(gpu_type) && nzchar(gpu_type)) {
    sprintf("gpu:%s:1", gpu_type)
  } else {
    "gpu:1"
  }

  # Build per-file cascadia commands
  file_cmds <- vapply(raw_files, function(f) {
    sample_name <- tools::file_path_sans_ext(basename(f))
    sprintf(
      paste0(
        '  echo "Processing: %s"\n',
        '  cascadia sequence \\\n',
        '    --input "%s" \\\n',
        '    --model "%s" \\\n',
        '    --output "%s/%s.ssl" \\\n',
        '    --batch-size %d \\\n',
        '    --min-score %s'
      ),
      basename(f), f, model_ckpt,
      denovo_dir, sample_name, batch_size, min_score
    )
  }, character(1))

  file_loop <- paste(file_cmds, collapse = "\n\n")

  # Assemble full sbatch script
  script <- paste0(
    '#!/bin/bash\n',
    sprintf('#SBATCH --job-name=cascadia_%s\n', safe_name),
    sprintf('#SBATCH --partition=%s\n', gpu_partition),
    sprintf('#SBATCH --account=%s\n', gpu_account),
    if (!is.null(gpu_qos) && nzchar(gpu_qos)) sprintf('#SBATCH --qos=%s\n', gpu_qos) else '',
    sprintf('#SBATCH --gres=%s\n', gres),
    sprintf('#SBATCH --cpus-per-task=%d\n', cpus),
    sprintf('#SBATCH --mem=%dG\n', mem_gb),
    sprintf('#SBATCH --time=%d:00:00\n', time_hours),
    sprintf('#SBATCH -o "%s/logs/cascadia_%%j.out"\n', output_dir),
    sprintf('#SBATCH -e "%s/logs/cascadia_%%j.err"\n', output_dir),
    if (isTRUE(requeue)) '#SBATCH --requeue\n' else '',
    '\n',
    '# Activate Cascadia conda environment\n',
    'source ~/.bashrc\n',
    sprintf('conda activate %s\n', cascadia_env),
    '\n',
    sprintf('echo "Cascadia de novo sequencing: %s"\n', analysis_name),
    'echo "Started: $(date)"\n',
    sprintf('echo "Output: %s"\n', denovo_dir),
    sprintf('echo "Model: %s"\n', model_ckpt),
    sprintf('echo "Min score: %s"\n', min_score),
    sprintf('echo "Files: %d"\n', length(raw_files)),
    '\n',
    '# Create output directory\n',
    sprintf('mkdir -p "%s"\n', denovo_dir),
    sprintf('mkdir -p "%s/logs"\n', output_dir),
    '\n',
    '# Verify GPU access\n',
    'echo "GPU info:"\n',
    'nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "WARNING: nvidia-smi failed"\n',
    '\n',
    '# Process each file\n',
    'N_SUCCESS=0\n',
    'N_FAIL=0\n',
    '\n',
    file_loop, '\n',
    '\n',
    '# Count SSL files produced\n',
    sprintf('N_SSL=$(ls "%s"/*.ssl 2>/dev/null | wc -l)\n', denovo_dir),
    'echo ""\n',
    sprintf('echo "Cascadia complete: ${N_SSL} / %d SSL files written"\n', length(raw_files)),
    'echo "Completed: $(date)"\n',
    '\n',
    sprintf('if [ "$N_SSL" -lt %d ]; then\n', length(raw_files)),
    '  echo "WARNING: Not all files produced SSL output"\n',
    '  exit 1\n',
    'fi\n',
    'exit 0\n'
  )

  script
}
