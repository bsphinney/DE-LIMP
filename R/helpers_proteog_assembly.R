# helpers_proteog_assembly.R — FASTA composition + registry I/O.
# No Shiny reactivity. Pure functions.
#
# Three responsibilities:
#  1) count_proteog_classes() — parse `source=` tags from FASTA headers and
#     count entries by class. Backward-compatible with v0.1 rewriter output
#     (which has no NOVEL_ISOFORM entries).
#  2) assemble_proteogenomics_fasta() — concatenate predicted ORFs + UniProt
#     reference (+ optional VARIANT + contaminants), dedup by sequence, and
#     register in the proteogenomics registry.
#  3) register_proteogenomics_fasta() / load_proteog_registry() — JSON I/O.

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# Single source of truth for the proteogenomics-database registry path.
.proteog_registry_path <- function() {
  Sys.getenv(
    "DELIMP_PROTEOG_REGISTRY",
    unset = "/quobyte/proteomics-grp/de-limp/databases/proteogenomics/registry.json"
  )
}

# Resolve seqkit binary path: prefer the proteog_helpers conda env (which is
# guaranteed to be present per Phase A), fall back to system PATH. This keeps
# assemble_proteogenomics_fasta() working both inside DE-LIMP's apptainer
# session and from a stock R invocation.
.find_seqkit <- function() {
  conda_env <- if (exists("PROTEOG_CONDA_ENV"))
    PROTEOG_CONDA_ENV
  else
    "/quobyte/proteomics-grp/de-limp/envs/proteog_helpers"
  candidate <- file.path(conda_env, "bin", "seqkit")
  if (file.exists(candidate)) return(candidate)
  on_path <- Sys.which("seqkit")
  if (nzchar(on_path)) return(unname(on_path))
  stop("assemble_proteogenomics_fasta(): seqkit not found. Tried ",
       candidate, " and system PATH. Install via `mamba install -p ",
       conda_env, " -c conda-forge -c bioconda seqkit`.")
}

# =============================================================================
# Input-integrity instrumentation (project convention — see NOTES_spec_lessons #14)
# =============================================================================
#
# Any helper that writes output FASTAs and takes input FASTAs MUST integrity-
# check its inputs. The check is cheap (one md5 of the first 1 KB per file,
# ~5ms each) and catches the entire class of "function corrupted my source
# data" bugs invisible at code-review time. CLAUDE.md Rule 4 applies — silent
# input mutation in an export-path helper is forbidden.

#' MD5 hash of the first n bytes of a file
#'
#' Uses the `digest` package when available; otherwise writes the head bytes
#' to a tempfile and uses base R's `tools::md5sum()`. The head-only hash
#' catches the corruption pattern we observed (whole-file truncation/
#' replacement) without paying the full-file scan cost on hundred-MB FASTAs.
.head_md5 <- function(path, n = 1024L) {
  if (!file.exists(path)) return(NA_character_)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  raw_bytes <- readBin(con, "raw", n = n)
  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(raw_bytes, algo = "md5", serialize = FALSE)
  } else {
    tmp <- tempfile(pattern = "head_md5_")
    on.exit(if (file.exists(tmp)) file.remove(tmp), add = TRUE)
    writeBin(raw_bytes, tmp)
    unname(tools::md5sum(tmp))
  }
}

#' Snapshot one input file's identity for later integrity verification
#'
#' Returns a list with $path, $size, $mtime, $head_hash — or NULL if the
#' input is NULL/missing (those are valid no-op skip cases for optional
#' inputs like variant_fasta).
.snapshot_input <- function(path) {
  if (is.null(path)) return(NULL)
  if (!file.exists(path)) return(NULL)
  list(
    path      = path,
    size      = file.size(path),
    mtime     = file.mtime(path),
    head_hash = .head_md5(path)
  )
}

#' Verify a snapshot against current file state
#'
#' @param snap snapshot from .snapshot_input(), or NULL (no-op)
#' @param strict if TRUE, stop() on mismatch; if FALSE, message() only
#' @return TRUE if unchanged, FALSE if changed (or stop()s when strict=TRUE)
.verify_input_unchanged <- function(snap, strict = TRUE) {
  if (is.null(snap)) return(TRUE)
  fail <- function(msg) {
    if (strict) stop(msg, call. = FALSE)
    message("[INPUT-INTEGRITY] ", msg)
    FALSE
  }
  if (!file.exists(snap$path)) {
    return(fail(sprintf("INPUT DELETED during execution: %s", snap$path)))
  }
  if (file.size(snap$path) != snap$size) {
    return(fail(sprintf(
      "INPUT SIZE CHANGED during execution: %s was %d bytes, now %d bytes",
      snap$path, snap$size, file.size(snap$path)
    )))
  }
  if (!identical(.head_md5(snap$path), snap$head_hash)) {
    return(fail(sprintf("INPUT CONTENT CHANGED during execution: %s", snap$path)))
  }
  TRUE
}

#' Log a disk-write operation to stderr
#'
#' Use this BEFORE every system2(..., stdout=path), file.copy(),
#' writeLines(), cat(file=path), jsonlite::write_json(...) call in the
#' assembly helpers. The first argument names the operation (for grep-ability
#' in the log); the second is the destination path string evaluated at the
#' moment of write.
.log_disk_write <- function(operation, path) {
  cat(sprintf("[DISK WRITE] %s -> %s\n", operation, path), file = stderr())
}

#' Validate that a file is a well-formed FASTA at function entry
#'
#' Companion to the integrity check: integrity catches MUTATION during the
#' function, this catches BAD-STATE-AT-ENTRY (e.g., a previously corrupted
#' input file). Together they cover both failure modes.
#'
#' Checks (cheap, ~5ms per file):
#'   - File exists
#'   - File size > 0
#'   - First non-empty line starts with `>` (FASTA header marker)
#'   - At least one `>` in the first 8 KB (basic sanity for multi-entry files)
#'
#' NULL/missing optional inputs are skipped (no-op return TRUE).
#'
#' @param path character path; NULL is OK (skip)
#' @param label character — argument name for the error message
#' @return TRUE on pass; stop() with clear message on fail
.validate_fasta_input <- function(path, label) {
  if (is.null(path) || !nzchar(path)) return(TRUE)
  if (!file.exists(path)) {
    stop(sprintf(
      "assemble_proteogenomics_fasta(): input %s does not exist: %s",
      label, path
    ), call. = FALSE)
  }
  sz <- file.size(path)
  if (is.na(sz) || sz <= 0) {
    stop(sprintf(
      "assemble_proteogenomics_fasta(): input %s is empty (size=%s): %s",
      label, as.character(sz), path
    ), call. = FALSE)
  }
  # Sniff first 8 KB for at least one FASTA header line.
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  raw_bytes <- readBin(con, "raw", n = 8192L)
  head_str <- rawToChar(raw_bytes[raw_bytes != as.raw(0)])
  if (!grepl("^>", head_str) && !grepl("\n>", head_str, fixed = TRUE)) {
    stop(sprintf(
      "assemble_proteogenomics_fasta(): input %s is not a valid FASTA (no '>' header found in first 8 KB): %s. File size = %d bytes. First 64 chars: %s",
      label, path, sz, substr(head_str, 1, 64)
    ), call. = FALSE)
  }
  TRUE
}

# =============================================================================
# Composition counting
# =============================================================================

#' Count proteins in a FASTA by `source=` class
#'
#' Reads only header lines (`^>`). For each header:
#'   - if it contains `source=...` → that class
#'   - else if accession matches `^(INDEL|SNV)_ENSP` → VARIANT
#'   - else → UNIPROT (canonical reference, no proteogenomic tag)
#'
#' Always returns a list with EVERY class field present (zero-filled), so
#' downstream consumers don't have to handle missing keys.
#'
#' @param fasta_path character — path to a FASTA file (plain or gzipped)
#' @return list with: total, UNIPROT, REF, NOVEL_GENE, NOVEL_ISOFORM, VARIANT, UNPARSED
count_proteog_classes <- function(fasta_path) {
  if (!file.exists(fasta_path)) {
    stop("count_proteog_classes(): file not found: ", fasta_path)
  }

  con <- if (grepl("\\.gz$", fasta_path, ignore.case = TRUE)) {
    gzfile(fasta_path, "rt")
  } else {
    file(fasta_path, "rt")
  }
  on.exit(close(con), add = TRUE)

  out <- list(
    total         = 0L,
    UNIPROT       = 0L,
    REF           = 0L,
    NOVEL_GENE    = 0L,
    NOVEL_ISOFORM = 0L,
    VARIANT       = 0L,
    UNPARSED      = 0L
  )

  # Stream line-by-line — FASTAs can be large (hundreds of MB)
  while (length(line <- readLines(con, n = 1L, warn = FALSE)) > 0) {
    if (!startsWith(line, ">")) next
    out$total <- out$total + 1L

    if (grepl("source=", line, fixed = TRUE)) {
      cls <- sub(".*source=([A-Z_]+).*", "\\1", line)
      if (cls %in% names(out)) {
        out[[cls]] <- out[[cls]] + 1L
      } else {
        out$UNPARSED <- out$UNPARSED + 1L
      }
    } else {
      # No source= tag — could be canonical UniProt or a VARIANT proteoform
      # (whose accession prefix marks it without needing a description tag).
      acc <- sub("^>([^[:space:]]+).*", "\\1", line)
      # Strip a leading sp|/tr| if present and grab the accession field.
      stripped <- sub("^(sp|tr)\\|([^|]+)\\|.*", "\\2", acc)
      if (grepl("^(INDEL|SNV)_ENSP", acc) || grepl("^(INDEL|SNV)_ENSP", stripped)) {
        out$VARIANT <- out$VARIANT + 1L
      } else {
        out$UNIPROT <- out$UNIPROT + 1L
      }
    }
  }
  out
}

# =============================================================================
# Registry I/O
# =============================================================================

#' Load the proteogenomics-database registry as a named list
#'
#' Empty/missing registry returns `list()`. Never throws — callers handle
#' an empty registry gracefully (no databases registered yet).
load_proteog_registry <- function() {
  path <- .proteog_registry_path()
  if (!file.exists(path)) return(list())
  raw <- tryCatch(jsonlite::read_json(path), error = function(e) NULL)
  if (is.null(raw)) {
    warning("load_proteog_registry(): could not parse ", path, " — returning empty registry")
    return(list())
  }
  if (length(raw) == 0) list() else raw
}

#' Write the proteogenomics-database registry from a named list
#'
#' @param registry named list (keys = project_name, values = entry list)
.save_proteog_registry <- function(registry) {
  path <- .proteog_registry_path()
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  .log_disk_write("jsonlite::write_json (registry)", path)
  jsonlite::write_json(registry, path, auto_unbox = TRUE, pretty = TRUE)
  invisible(path)
}

#' Register a built proteogenomics FASTA in the registry
#'
#' @param path             character — FASTA path
#' @param merged_gtf_path  character — preserved StringTie merged GTF alongside
#' @param project_name     character — user-supplied project identifier
#' @param composition      list — from count_proteog_classes()
#' @param build_metadata   list — optional extra fields (organism, build, n_samples, …)
#' @return invisible(path)
register_proteogenomics_fasta <- function(path,
                                          merged_gtf_path,
                                          project_name,
                                          composition,
                                          build_metadata = list()) {
  registry <- load_proteog_registry()

  entry <- c(
    list(
      path             = path,
      merged_gtf_path  = merged_gtf_path,
      project_name     = project_name,
      composition      = composition,
      created          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      created_by       = Sys.getenv("USER", unset = NA_character_),
      pipeline_version = "1.1"
    ),
    build_metadata
  )

  registry[[project_name]] <- entry
  .save_proteog_registry(registry)
  invisible(path)
}

# =============================================================================
# Assembly
# =============================================================================

#' Build a proteogenomics search database from component FASTAs
#'
#' Concatenates predicted ORFs + UniProt reference (+ optional variants +
#' optional cRAP contaminants), optionally dedups by sequence with `seqkit
#' rmdup -s`, copies the merged GTF alongside for coordinate lookups, and
#' registers the result.
#'
#' Self-describing: returned object carries `pipeline_id`, `methods_paragraph`,
#' `composition` — downstream consumers must read these rather than hardcoding
#' a description of "what we did" (CLAUDE.md architectural rule #1).
#'
#' @param project_name        character
#' @param uniprot_fasta       character path
#' @param predicted_orfs_fasta character path (from header rewriter)
#' @param merged_gtf          character path — StringTie merged.gtf
#' @param variant_fasta       character path or NULL
#' @param contaminants_fasta  character path or NULL
#' @param output_dir          character — where the final FASTA + GTF land
#' @param dedupe              logical — run `seqkit rmdup -s`
#' @param build_metadata      list — optional metadata for registry (organism, etc.)
#' @return list with $path, $merged_gtf_path, $composition, $pipeline_id, $methods_paragraph
assemble_proteogenomics_fasta <- function(project_name,
                                          uniprot_fasta,
                                          predicted_orfs_fasta,
                                          merged_gtf,
                                          variant_fasta = NULL,
                                          contaminants_fasta = NULL,
                                          output_dir,
                                          dedupe = TRUE,
                                          build_metadata = list()) {
  # Validate inputs (boundary check — CLAUDE.md rule on input validation).
  required_files <- c(uniprot_fasta, predicted_orfs_fasta, merged_gtf)
  missing <- required_files[!file.exists(required_files)]
  if (length(missing) > 0) {
    stop("assemble_proteogenomics_fasta(): missing required input(s): ",
         paste(missing, collapse = ", "))
  }
  if (!nzchar(project_name) || !grepl("^[A-Za-z0-9_.-]+$", project_name)) {
    stop("assemble_proteogenomics_fasta(): project_name must be non-empty and contain only [A-Za-z0-9_.-]; got: ", project_name)
  }

  # Structural validation: each FASTA input must be a well-formed FASTA at
  # entry. Catches the "previously-corrupted-input" class of bug (e.g., a
  # 2-byte "0\n" file from an upstream shell-meta accident) before any
  # downstream tool sees it. merged_gtf is excluded — it's a GTF, not FASTA.
  .validate_fasta_input(uniprot_fasta,        "uniprot_fasta")
  .validate_fasta_input(predicted_orfs_fasta, "predicted_orfs_fasta")
  .validate_fasta_input(variant_fasta,        "variant_fasta")
  .validate_fasta_input(contaminants_fasta,   "contaminants_fasta")

  # ── Input-integrity snapshot (project convention; NOTES_spec_lessons #14) ──
  # Any code path that mutates any of these snapshotted files mid-execution
  # is a bug. We verify both on the success path (final stop() if changed)
  # and on the failure path (on.exit message() with no stop, so the original
  # error still propagates).
  input_snapshots <- list(
    uniprot      = .snapshot_input(uniprot_fasta),
    predicted    = .snapshot_input(predicted_orfs_fasta),
    merged_gtf   = .snapshot_input(merged_gtf),
    variant      = .snapshot_input(variant_fasta),
    contaminants = .snapshot_input(contaminants_fasta)
  )
  on.exit({
    for (nm in names(input_snapshots)) {
      .verify_input_unchanged(input_snapshots[[nm]], strict = FALSE)
    }
  }, add = TRUE)

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  out_path <- file.path(
    output_dir,
    sprintf("%s_proteogenomics_%s.fasta",
            project_name,
            format(Sys.Date(), "%Y_%m"))
  )

  # Compose the component list in deterministic order: predicted ORFs first
  # so their headers come up early in the FASTA (useful for spot-checking),
  # then UniProt, then optional VARIANT, then optional contaminants.
  components <- c(predicted_orfs_fasta, uniprot_fasta)
  if (!is.null(variant_fasta) && file.exists(variant_fasta)) {
    components <- c(components, variant_fasta)
  }
  if (!is.null(contaminants_fasta) && file.exists(contaminants_fasta)) {
    components <- c(components, contaminants_fasta)
  }

  if (dedupe) {
    tmp_concat <- tempfile(pattern = "proteog_concat_", fileext = ".fasta")
    on.exit(if (file.exists(tmp_concat)) file.remove(tmp_concat), add = TRUE)

    .log_disk_write("system2('cat', stdout=tmp_concat)", tmp_concat)
    cat_status <- system2(
      "cat", args = shQuote(components),
      stdout = tmp_concat, stderr = FALSE
    )
    if (cat_status != 0) {
      stop("assemble_proteogenomics_fasta(): cat failed concatenating components")
    }

    seqkit_bin <- .find_seqkit()
    .log_disk_write("seqkit rmdup -o out_path", out_path)
    rmdup_status <- system2(
      seqkit_bin,
      args = c("rmdup", "-s", "-o", shQuote(out_path), shQuote(tmp_concat))
    )
    if (rmdup_status != 0) {
      stop("assemble_proteogenomics_fasta(): seqkit rmdup failed (exit ", rmdup_status, ")")
    }
  } else {
    .log_disk_write("system2('cat', stdout=out_path)", out_path)
    cat_status <- system2(
      "cat", args = shQuote(components),
      stdout = out_path, stderr = FALSE
    )
    if (cat_status != 0) {
      stop("assemble_proteogenomics_fasta(): cat failed concatenating components (no dedupe)")
    }
  }

  composition <- count_proteog_classes(out_path)

  # Preserve merged GTF alongside (coordinate-lookup reference) — registry
  # records the path so the Proteogenomics tab can resolve genomic coords on demand.
  gtf_dest <- sub("\\.fasta$", "_merged.gtf", out_path)
  .log_disk_write("file.copy(merged_gtf -> gtf_dest)", gtf_dest)
  file.copy(merged_gtf, gtf_dest, overwrite = TRUE)

  register_proteogenomics_fasta(
    path             = out_path,
    merged_gtf_path  = gtf_dest,
    project_name     = project_name,
    composition      = composition,
    build_metadata   = build_metadata
  )

  # Strict integrity check on the success path — stop() if any input changed.
  for (nm in names(input_snapshots)) {
    .verify_input_unchanged(input_snapshots[[nm]], strict = TRUE)
  }

  list(
    pipeline_id       = "proteogenomics_v1.1",
    path              = out_path,
    merged_gtf_path   = gtf_dest,
    composition       = composition,
    methods_paragraph = .proteog_methods_paragraph(composition, build_metadata)
  )
}

# Self-describing methods text used by downstream Claude export / methods readme.
# CLAUDE.md rule #1: pipeline objects must self-describe; downstream consumers
# read $methods_paragraph rather than hardcoding what the pipeline did.
.proteog_methods_paragraph <- function(composition, meta = list()) {
  fmt <- function(n) format(n %||% 0L, big.mark = ",")
  paste0(
    "The spectral search database was constructed by concatenating ",
    (meta$uniprot_release %||% "the UniProt reference proteome"),
    " with StringTie-predicted ORFs from sample-matched RNA-seq (",
    fmt(composition$REF), " reference-derived predicted, ",
    fmt(composition$NOVEL_GENE), " novel-gene candidates, ",
    fmt(composition$NOVEL_ISOFORM), " novel-isoform candidates",
    if ((composition$VARIANT %||% 0L) > 0) {
      paste0(", ", fmt(composition$VARIANT), " variant proteoforms")
    } else "",
    "). Total entries after deduplication: ", fmt(composition$total),
    " (", fmt(composition$UNIPROT), " canonical UniProt + ",
    fmt(composition$REF + composition$NOVEL_GENE + composition$NOVEL_ISOFORM),
    " proteogenomic).",
    if (!is.null(meta$rrna_pct_mean)) {
      sprintf(" Mean rRNA filter rate: %.1f%%.", meta$rrna_pct_mean)
    } else "",
    if (!is.null(meta$uniquely_mapped_pct_mean)) {
      sprintf(" Mean uniquely-mapped rate: %.1f%%.", meta$uniquely_mapped_pct_mean)
    } else "",
    if (!is.null(meta$read_length_tier)) {
      sprintf(" STAR threshold tier: %s.", meta$read_length_tier)
    } else "",
    if (isTRUE(meta$contaminants_used) || identical(meta$contaminants_source, "HaoGroup")) {
      paste0(
        " Protein contaminant sequences were obtained from the universal protein",
        " contaminant library (Frankenfield et al. 2022,",
        " doi:10.1021/acs.jproteome.2c00145;",
        " https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics)."
      )
    } else ""
  )
}

#' Load provenance metadata for a static asset (contaminant FASTA, etc.)
#'
#' Reads provenance.json from the same directory as the asset and returns
#' the entry for that asset, or NULL if no provenance is found.
#'
#' This is the "one filesystem lookup" implementation of NOTES_spec_lessons #15.
#' Any DE-LIMP code that needs the citation, license, or source URL for a
#' shipped static file should call this rather than hardcoding the metadata.
#'
#' @param asset_path character — absolute path to the asset file
#' @return named list with provenance fields, or NULL if not found
load_asset_provenance <- function(asset_path) {
  if (is.null(asset_path) || !nzchar(asset_path)) return(NULL)
  prov_path <- file.path(dirname(asset_path), "provenance.json")
  if (!file.exists(prov_path)) return(NULL)
  prov <- tryCatch(
    jsonlite::read_json(prov_path),
    error = function(e) NULL
  )
  if (is.null(prov)) return(NULL)
  entry_name <- basename(asset_path)
  prov[[entry_name]] %||% NULL
}
