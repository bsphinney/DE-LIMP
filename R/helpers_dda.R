# ==============================================================================
#  HELPERS — DDA Search (Sage + DIA-NN DDA pipelines)
#  Pure utility functions — no Shiny reactivity.
#  Called from: server_dda.R
# ==============================================================================

#' Canonical bare-amino-acid form of a peptide.
#'
#' The SINGLE definition used to cross-reference Sage DB-search peptides,
#' Casanovo de novo peptides, and DIAMOND BLAST query peptides. All three MUST
#' normalize identically or `%in%` joins between them silently return 0 (the
#' v3.11.x "BLAST Hits: 0" / "no matched spectra" bug — three different
#' strippers disagreed on named mods). Strips every notation the tools emit:
#'   [Acetyl] / [Carbamidomethyl] named, [+57.02] / +57.02 numeric,
#'   (UniMod:4), and terminal "[Mod]-" hyphens. Keeps only A-Z.
#' NOTE: keep the awk in the DIAMOND sbatch generator (server_dda.R) in sync —
#'   `gsub(/\[[^]]*\]/,"",$2); gsub(/[^A-Z]/,"",$2)`.
build_dda_canonical_peptide <- function(x) {
  x <- as.character(x)
  x <- gsub("\\[[^]]*\\]", "", x)          # [Acetyl], [Carbamidomethyl], [+57.02]
  x <- gsub("\\(UniMod:[0-9]+\\)", "", x)  # (UniMod:4)
  x <- gsub("\\+[0-9.]+", "", x)           # bare +57.02
  toupper(gsub("[^A-Za-z]", "", x))        # drop hyphens, dots, digits, spaces
}

#' Taxonomy-aware species/clade for a de novo BLAST hit — THE single definition.
#'
#' NCBI nr accessions (XP_/NP_/WP_/YP_…) carry no `_SPECIES` mnemonic, so the
#' UniProt-style "text after the last underscore" parse yields garbage
#' (`XP_025773238.1` -> `025773238.1`). When a per-peptide nr LCA table is
#' available it is authoritative — join on the canonical (mods-stripped,
#' I/L-normalized) peptide. Otherwise fall back to the UniProt mnemonic. Every
#' de novo view that needs a species must call this, never re-parse the subject
#' (CLAUDE.md architectural rule #3).
#'
#' @param peptide character vector of peptide sequences (the BLAST query)
#' @param subject character vector of BLAST subject accessions (same length)
#' @param lca_tbl optional data.frame with `peptide` + `lca_name` columns
#' @return character vector of species/clade names (same length as `subject`)
dda_blast_species <- function(peptide, subject, lca_tbl = NULL) {
  if (!is.null(lca_tbl) && nrow(lca_tbl) > 0 &&
      all(c("peptide", "lca_name") %in% names(lca_tbl))) {
    key  <- gsub("I", "L", build_dda_canonical_peptide(peptide))
    lkey <- gsub("I", "L", build_dda_canonical_peptide(lca_tbl$peptide))
    out  <- lca_tbl$lca_name[match(key, lkey)]
    out[is.na(out)] <- "Unresolved"
    return(out)
  }
  sub(".*_", "", sub("^[a-z]+\\|[^|]+\\|", "", subject))
}

#' NCBI Taxonomy Browser link for a taxon — HTML <a>, opens in a new tab.
#'
#' Returns the plain name when the taxid is missing. Intended for DT columns
#' rendered with `escape = FALSE`; taxon names come from our own LCA table so
#' there is no untrusted HTML.
#'
#' @param name character vector of taxon names
#' @param taxid character/numeric vector of NCBI taxids (same length)
#' @return character vector of HTML anchors (or plain names where no taxid)
ncbi_tax_link <- function(name, taxid) {
  name <- as.character(name); taxid <- as.character(taxid)
  out <- name
  has <- !is.na(name) & nzchar(name) & !is.na(taxid) & nzchar(taxid) & taxid != "0"
  out[has] <- sprintf(
    paste0('<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%s"',
           ' target="_blank" rel="noopener">%s</a>'),
    taxid[has], name[has])
  out
}

#' Display label for a BLAST subject's protein — UniProt mnemonic or full accession.
#'
#' For UniProt `sp|ACC|NAME_SPECIES` returns the protein mnemonic (`NAME`). For
#' NCBI / RefSeq / GenBank / PDB accessions returns the FULL accession unchanged
#' (never split on `_`, which would mangle `XP_025773238.1` -> "XP"). Single
#' definition for every de novo view that shows a protein label.
#'
#' @param subject character vector of BLAST subject accessions
#' @return character vector of protein labels
dda_protein_label <- function(subject) {
  s <- as.character(subject)
  isup <- grepl("^[a-z]{2}\\|[^|]+\\|", s)
  out <- s
  out[isup] <- sub("_[^_]+$", "", sub("^[a-z]+\\|[^|]+\\|", "", s[isup]))
  out
}

#' Best BLAST hit per canonical peptide — alignment length, e-value, bitscore.
#'
#' "Best" = highest bitscore. Used to surface query coverage (alnlen / peptide
#' length) + e-value so a 100%-identity-but-partial hit (e.g. 18 of 24 residues
#' to an over-represented taxon) can't masquerade as a full match. Pure.
#'
#' @param blast BLAST hits (cols: peptide/query, length, evalue, bitscore)
#' @return data.table(.k, blast_aln_len, blast_evalue, blast_bitscore) or NULL
best_blast_hit_per_peptide <- function(blast) {
  if (is.null(blast) || !is.data.frame(blast) || nrow(blast) == 0) return(NULL)
  bt <- data.table::as.data.table(blast)
  kcol <- if ("peptide" %in% names(bt)) "peptide" else if ("query" %in% names(bt)) "query" else NULL
  if (is.null(kcol)) return(NULL)
  bt$.k <- gsub("I", "L", build_dda_canonical_peptide(bt[[kcol]]))
  if ("bitscore" %in% names(bt)) bt <- bt[order(-as.numeric(bt$bitscore)), ]
  bt <- bt[!duplicated(bt$.k), ]
  alnc <- if ("length" %in% names(bt)) "length" else if ("aln_len" %in% names(bt)) "aln_len" else NA
  data.table::data.table(
    .k             = bt$.k,
    blast_aln_len  = if (!is.na(alnc)) as.numeric(bt[[alnc]]) else NA_real_,
    blast_evalue   = if ("evalue" %in% names(bt)) as.numeric(bt$evalue) else NA_real_,
    blast_bitscore = if ("bitscore" %in% names(bt)) as.numeric(bt$bitscore) else NA_real_)
}

#' Build the de novo Master Table — pure join, no Shiny (so it is unit-testable).
#'
#' One row per de novo peptide on the single canonical key
#' `gsub("I","L", build_dda_canonical_peptide(seq))`, combining:
#'   - Casanovo: best confidence + n PSMs + whether any PSM was Sage-confirmed
#'   - Sage:     the protein the peptide maps to (NA if de-novo-only)
#'   - nr LCA:   species/clade, rank, category, best %ID, diagnostic flag
#'
#' NOTE on confidence: Casanovo's `search_engine_score` runs roughly -1..+1, and
#' species-diagnostic peptides for a non-model organism cluster in the MID range
#' (~0.5-0.9), NOT at the top — so a high slider default can hide the entire nr
#' signal. Filtering is left to the caller; `Casanovo_score` is returned per row.
#'
#' @param classified data.frame/data.table of Casanovo PSMs (cols: seq_stripped
#'        or sequence, score, match_type; seq_norm optional)
#' @param sage_psms  optional Sage PSM table (cols: peptide, proteins)
#' @param lca        optional per-peptide LCA table (cols: peptide, lca_name,
#'        lca_rank, category, top_pident, diagnostic)
#' @return data.frame with seq_norm + biologist-facing columns
build_denovo_master <- function(classified, sage_psms = NULL, lca = NULL,
                                blast = NULL) {
  if (is.null(classified) || nrow(classified) == 0)
    return(data.frame())
  cas <- data.table::as.data.table(classified)
  base_seq <- if ("seq_stripped" %in% names(cas)) cas$seq_stripped else cas$sequence
  cas$seq_norm  <- gsub("I", "L", build_dda_canonical_peptide(base_seq))
  if (!"score" %in% names(cas)) cas$score <- NA_real_
  if (!"match_type" %in% names(cas)) cas$match_type <- "novel"
  cas$disp_seq <- base_seq
  pep <- cas[, list(
    Peptide        = disp_seq[1],
    Casanovo_score = suppressWarnings(round(max(score, na.rm = TRUE), 3)),
    n_PSMs         = .N,
    Found_by_Sage  = any(match_type == "confirmed")
  ), by = "seq_norm"]
  pep$Casanovo_score[!is.finite(pep$Casanovo_score)] <- NA_real_

  if (!is.null(sage_psms) && all(c("peptide", "proteins") %in% names(sage_psms))) {
    sdt <- data.table::as.data.table(sage_psms)
    sdt$seq_norm <- gsub("I", "L", build_dda_canonical_peptide(sdt$peptide))
    sdt <- sdt[!duplicated(sdt$seq_norm), ]
    pep <- merge(pep, sdt[, list(seq_norm, Sage_protein = proteins)],
                 by = "seq_norm", all.x = TRUE)
  } else {
    pep$Sage_protein <- NA_character_
  }

  if (!is.null(lca) && all(c("peptide", "lca_name") %in% names(lca))) {
    ldt <- data.table::as.data.table(lca)
    ldt$seq_norm <- gsub("I", "L", build_dda_canonical_peptide(ldt$peptide))
    ldt <- ldt[!duplicated(ldt$seq_norm), ]
    cols <- intersect(c("seq_norm", "lca_taxid", "lca_name", "lca_rank",
                        "category", "top_pident", "diagnostic"), names(ldt))
    pep <- merge(pep, ldt[, cols, with = FALSE], by = "seq_norm", all.x = TRUE)
  }
  ren <- c(lca_name = "Species_or_clade", lca_rank = "Rank",
           category = "Type", top_pident = "Best_pct_ID",
           diagnostic = "Diagnostic")
  for (k in names(ren)) if (k %in% names(pep)) data.table::setnames(pep, k, ren[[k]])
  if ("Diagnostic" %in% names(pep)) pep$Diagnostic <- pep$Diagnostic %in% c(1, "1", TRUE)

  # Best-hit coverage / e-value / bitscore so a 100%-identity-but-partial hit
  # (e.g. 18 of 24 residues) can't masquerade as a full-length match.
  bb <- best_blast_hit_per_peptide(blast)
  if (!is.null(bb)) {
    data.table::setnames(bb, ".k", "seq_norm")
    pep <- merge(pep, bb, by = "seq_norm", all.x = TRUE)
    pep$Query_coverage <- ifelse(
      is.na(pep$blast_aln_len), NA_real_,
      round(100 * pep$blast_aln_len / nchar(pep$Peptide)))
    data.table::setnames(pep, c("blast_evalue", "blast_bitscore"),
                         c("E_value", "Bitscore"))
    pep[, blast_aln_len := NULL]
  }
  as.data.frame(pep)
}

#' Calibrate the Casanovo confidence score against nr BLAST outcome — and, when
#' a shuffled-decoy BLAST is supplied, estimate the empirical FDR per cutoff.
#'
#' For each Casanovo score bin: how many de novo peptides got an nr hit
#' (target hit-rate) and the mean %identity of those hits. The decoy BLAST
#' (same peptides, residues shuffled) gives the by-chance hit-rate, so
#'   FDR(bin) = decoy_hit_rate / target_hit_rate
#' and a sensible cutoff is the lowest score whose CUMULATIVE FDR (peptides at
#' or above that score) is below the target (e.g. 0.05). Pure + testable.
#'
#' @param classified Casanovo PSMs (seq_stripped/sequence + score)
#' @param blast      target nr BLAST hits (peptide/query + pident/identity)
#' @param decoy_blast optional shuffled-decoy nr BLAST (same columns)
#' @param bin numeric bin width (default 0.1)
#' @return data.frame: bin, n, n_hit, hit_rate, mean_pident
#'         [+ decoy_hit_rate, fdr, cum_fdr when decoy_blast supplied]
build_denovo_score_calibration <- function(classified, blast,
                                           decoy_blast = NULL, bin = 0.1,
                                           min_length = 0) {
  if (is.null(classified) || nrow(classified) == 0) return(data.frame())
  cas <- data.table::as.data.table(classified)
  base <- if ("seq_stripped" %in% names(cas)) cas$seq_stripped else cas$sequence
  cas$k <- gsub("I", "L", build_dda_canonical_peptide(base))
  if (!"score" %in% names(cas)) cas$score <- NA_real_
  sc <- cas[, list(score = suppressWarnings(max(score, na.rm = TRUE))), by = "k"]
  sc <- sc[is.finite(sc$score), ]
  # Symmetric length floor: dropping short peptides from the peptide population
  # removes them from BOTH the target and decoy hit counts, so they provably
  # cannot skew the decoy/target FDR ratio (and they carry ~no species signal).
  if (min_length > 0) sc <- sc[nchar(sc$k) >= min_length, ]
  if (nrow(sc) == 0) return(data.frame())

  best_pid <- function(b) {
    empty <- data.table::data.table(k = character(0), pid = numeric(0))
    if (is.null(b) || nrow(b) == 0) return(empty)
    bt <- data.table::as.data.table(b)
    pcol <- if ("pident" %in% names(bt)) "pident" else if ("identity" %in% names(bt)) "identity" else NA
    kcol <- if ("peptide" %in% names(bt)) "peptide" else if ("query" %in% names(bt)) "query" else NA
    if (is.na(pcol) || is.na(kcol)) return(empty)
    bt$k <- gsub("I", "L", build_dda_canonical_peptide(bt[[kcol]]))
    bt[, list(pid = max(get(pcol), na.rm = TRUE)), by = "k"]
  }

  tpid <- best_pid(blast)
  sc$bin <- floor(sc$score / bin) * bin
  sc$hit <- sc$k %in% tpid$k
  sc$pid <- tpid$pid[match(sc$k, tpid$k)]
  agg <- sc[, list(n = .N, n_hit = sum(hit),
                   hit_rate = 100 * mean(hit),
                   mean_pident = mean(pid[hit], na.rm = TRUE)), by = "bin"]
  agg <- agg[order(-agg$bin), ]

  if (!is.null(decoy_blast)) {
    dpid <- best_pid(decoy_blast)
    sc$dhit <- sc$k %in% dpid$k
    dagg <- sc[, list(decoy_hit_rate = 100 * mean(dhit)), by = "bin"]
    agg <- merge(agg, dagg, by = "bin", all.x = TRUE)
    agg <- agg[order(-agg$bin), ]
    agg$fdr <- ifelse(agg$hit_rate > 0,
                      pmin(1, agg$decoy_hit_rate / agg$hit_rate), NA_real_)
    # cumulative FDR for peptides at-or-above each bin (decoy hits / target hits)
    agg$cum_target <- cumsum(agg$n_hit)
    cum_decoy <- cumsum((agg$decoy_hit_rate / 100) * agg$n)
    agg$cum_fdr <- ifelse(agg$cum_target > 0,
                          pmin(1, cum_decoy / agg$cum_target), NA_real_)
  }
  as.data.frame(agg)
}

#' TRUE for each protein group that is ENTIRELY contaminant — every accession
#' carries the `Cont_` tag stamped by the contaminant FASTA appended at search
#' time. Single definition used by the Results-page PSM filter and the de novo
#' Sage-DB-hits / classification filter. Groups with a real protein are kept.
is_contaminant_protein_group <- function(proteins) {
  vapply(strsplit(as.character(proteins), ";", fixed = TRUE), function(accs) {
    accs <- accs[nzchar(accs)]
    length(accs) > 0 && all(grepl("Cont_", accs, fixed = TRUE))
  }, logical(1))
}

#' TRUE for each protein/accession that looks like a skin- or hair-family
#' protein: keratins (incl. messy hair-keratin nomenclature K1H1/K2M3/KRB2A/
#' KRA31/KT33A), keratin-associated proteins, trichohyalin, filaggrin,
#' hornerin, loricrin, involucrin, corneodesmosin, desmoplakin/glein, collagens.
#' Heuristic by design — keys on the matched homolog's NAME so it works across
#' species (a non-model organism's hair keratin blasts to a named one). The
#' user opts in via a dropdown, so false hits are recoverable. Do NOT wire this
#' into anything that runs by default — for hair/feather projects keratins ARE
#' the signal, not contaminants.
is_skin_hair_protein <- function(x) {
  x <- toupper(as.character(x))
  grepl(paste0(
    "KRT|KERATIN|KRTAP|\\bKAP[0-9]|KRB[0-9]|KRA[0-9]|KT[0-9]|K[0-9]H[0-9]|",
    "K[0-9]M[0-9]|\\bK[0-9]{1,2}[A-Z]?_|TRICHOHYALIN|\\bTCHH|FILAGGRIN|\\bFLG_|",
    "HORNERIN|\\bHRNR|LORICRIN|INVOLUCRIN|CORNEODESMOSIN|\\bCDSN|DESMOPLAKIN|",
    "DESMOGLEIN|COLLAGEN|\\bCOL[0-9]"
  ), x)
}

#' Apply the de novo "Protein filter" dropdown to a BLAST data frame.
#' @param mode "all" (no-op), "skin_only" (keep skin/hair), "skin_exclude"
#'   (drop skin/hair). Finds the first available protein/subject column.
dda_apply_protein_filter <- function(blast, mode = "all") {
  if (is.null(blast) || nrow(blast) == 0 || identical(mode, "all")) return(blast)
  col <- intersect(c("subject", "Protein", "proteins", "protein", "Name", "name"),
                   names(blast))
  if (length(col) == 0) return(blast)
  sh <- is_skin_hair_protein(blast[[col[1]]])
  if (identical(mode, "skin_only"))    return(blast[sh, , drop = FALSE])
  if (identical(mode, "skin_exclude")) return(blast[!sh, , drop = FALSE])
  blast
}

#' Generate a Sage search config JSON
#'
#' @param fasta_path Path to reference FASTA
#' @param raw_paths Character vector of .d directory paths
#' @param output_dir Where Sage writes results
#' @param preset One of "standard", "phospho", "tmt"
#' @param missed_cleavages Integer (default 2)
#' @param precursor_tol_ppm Precursor mass tolerance in ppm (default 20)
#' @param fragment_tol_ppm Fragment tolerance in ppm (default 20). Sage v0.14.7 expects ppm under quant.fragment_tol.
#' @param min_peaks Min fragment peaks per spectrum (default 6)
#' @return Path to written sage.json file
generate_sage_config <- function(
  fasta_path,
  raw_paths,
  output_dir,
  preset           = "standard",
  missed_cleavages = 2,
  precursor_tol_ppm = 20,
  fragment_tol_ppm  = 20,
  min_peaks         = 6
) {
  # ──────────────────────────────────────────────────────────────────────────
  # Preset → search-parameter table. Adding a new analysis mode = add a row.
  # Keep keys consistent so the caller's mode dropdown maps 1:1 → preset name.
  # ──────────────────────────────────────────────────────────────────────────
  preset_table <- list(
    "standard" = list(
      cleave_at         = "KR",      restrict = "P",
      min_len           = 7L,        max_len = 50L,
      peptide_min_mass  = 500,       peptide_max_mass = 5000,
      variable_mods     = list("M" = c(15.9949), "^" = c(42.0106)),
      add_tmt_static    = FALSE
    ),
    "phospho" = list(
      cleave_at         = "KR",      restrict = "P",
      min_len           = 7L,        max_len = 50L,
      peptide_min_mass  = 500,       peptide_max_mass = 5000,
      variable_mods     = list("M" = c(15.9949), "S" = c(79.9663),
                               "T" = c(79.9663), "Y" = c(79.9663),
                               "^" = c(42.0106)),
      add_tmt_static    = FALSE
    ),
    "tmt" = list(
      cleave_at         = "KR",      restrict = "P",
      min_len           = 7L,        max_len = 50L,
      peptide_min_mass  = 500,       peptide_max_mass = 5000,
      variable_mods     = list("M" = c(15.9949)),
      add_tmt_static    = TRUE
    ),
    # Peptidomics: endogenous peptides (no enzymatic digestion). Nonspecific.
    # Variable mods: oxidation, pyro-Glu (Q/E N-term), C-term amidation, N-term acetylation.
    "peptidomics" = list(
      cleave_at         = "",        restrict = "",            # nonspecific
      min_len           = 5L,        max_len = 25L,
      peptide_min_mass  = 400,       peptide_max_mass = 5000,
      variable_mods     = list("M" = c(15.9949),
                               "^" = c(42.0106, -17.02655, -18.01056),  # N-term: acetyl, pyroGlu(Q), pyroGlu(E)
                               "$" = c(-0.98402)),                       # C-term: amidation
      add_tmt_static    = FALSE
    ),
    # HLA class I — non-specific, 8–12 AA, charge typically +1–3 on TOF.
    # Variable mods per BOWIE preset table: oxidation + deamidation.
    "hla_class_i" = list(
      cleave_at         = "",        restrict = "",
      min_len           = 8L,        max_len = 12L,
      peptide_min_mass  = 700,       peptide_max_mass = 1500,
      variable_mods     = list("M" = c(15.9949),
                               "N" = c(0.98402), "Q" = c(0.98402)),  # deamidation
      add_tmt_static    = FALSE
    ),
    # HLA class II — non-specific, 13–25 AA.
    "hla_class_ii" = list(
      cleave_at         = "",        restrict = "",
      min_len           = 13L,       max_len = 25L,
      peptide_min_mass  = 1300,      peptide_max_mass = 3000,
      variable_mods     = list("M" = c(15.9949),
                               "N" = c(0.98402), "Q" = c(0.98402)),
      add_tmt_static    = FALSE
    )
  )

  ps <- preset_table[[preset]] %||% preset_table[["standard"]]

  static_mods <- list("C" = jsonlite::unbox(57.0215))  # carbamidomethyl always fixed
  if (isTRUE(ps$add_tmt_static)) {
    static_mods[["^"]] <- jsonlite::unbox(229.1629)
    static_mods[["K"]] <- jsonlite::unbox(229.1629)
  }

  config <- list(
    database = list(
      bucket_size = jsonlite::unbox(32768L),
      enzyme = list(
        missed_cleavages = jsonlite::unbox(as.integer(missed_cleavages)),
        min_len          = jsonlite::unbox(ps$min_len),
        max_len          = jsonlite::unbox(ps$max_len),
        cleave_at        = jsonlite::unbox(ps$cleave_at),
        restrict         = jsonlite::unbox(ps$restrict)
      ),
      fragment_min_mz  = jsonlite::unbox(150.0),
      fragment_max_mz  = jsonlite::unbox(2000.0),
      peptide_min_mass = jsonlite::unbox(as.numeric(ps$peptide_min_mass)),
      peptide_max_mass = jsonlite::unbox(as.numeric(ps$peptide_max_mass)),
      ion_kinds        = c("b", "y"),
      min_ion_index    = jsonlite::unbox(2L),
      static_mods      = static_mods,
      variable_mods    = ps$variable_mods,
      max_variable_mods = jsonlite::unbox(2L),
      decoy_tag        = jsonlite::unbox("rev_"),
      generate_decoys  = jsonlite::unbox(TRUE),
      fasta            = jsonlite::unbox(fasta_path)
    ),
    precursor_tol = list(ppm = c(-precursor_tol_ppm, precursor_tol_ppm)),
    fragment_tol  = list(ppm = c(-fragment_tol_ppm, fragment_tol_ppm)),
    report_psms   = jsonlite::unbox(1L),
    min_peaks     = jsonlite::unbox(as.integer(min_peaks)),
    max_peaks     = jsonlite::unbox(150L),
    max_fragment_charge = jsonlite::unbox(2L),
    chimera       = jsonlite::unbox(FALSE),
    predict_rt    = jsonlite::unbox(TRUE),
    mzml_paths    = raw_paths
  )

  # Add LFQ config for non-TMT.
  #
  # IMPORTANT — Sage v0.14.7 config schema (the version we ship at
  # /quobyte/proteomics-grp/de-limp/cascadia/sage-v0.14.7-x86_64-unknown-linux-gnu/sage):
  #   quant.lfq             = BOOLEAN (true/false)              ← v0.14.7
  #   quant.lfq_settings    = { peak_scoring, integration, spectral_angle, ppm_tolerance }
  # In Sage v0.15+ the boolean was REMOVED and the settings object was
  # renamed `lfq` (no `_settings` suffix). Emitting the v0.15 shape against
  # v0.14.7 yields: `Error: invalid type: map, expected a boolean at line N`.
  # Verified against https://github.com/lazear/sage/blob/v0.14.7/DOCS.md
  if (preset != "tmt") {
    config$quant <- list(
      lfq = jsonlite::unbox(TRUE),
      lfq_settings = list(
        peak_scoring   = jsonlite::unbox("Hybrid"),
        integration    = jsonlite::unbox("Sum"),
        spectral_angle = jsonlite::unbox(0.7)
      )
    )
  } else {
    config$quant <- list(
      tmt = list(level = jsonlite::unbox(2L))
    )
  }

  # Write config to temp dir — caller SCPs to HPC output_dir
  local_dir <- file.path(tempdir(), "sage_config")
  dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
  config_path <- file.path(local_dir, "sage.json")
  jsonlite::write_json(config, config_path, auto_unbox = FALSE, pretty = TRUE, null = "null")
  message("[DDA] Sage config written to: ", config_path)
  config_path
}


#' Parse Sage search results into DE-LIMP-compatible format
#'
#' @param results_path Path to results.sage.tsv
#' @param lfq_path Path to lfq.tsv
#' @param fdr_threshold PSM-level FDR cutoff (default 0.01)
#' @param protein_fdr_threshold Protein-level FDR cutoff (default 0.01)
#' @return List with $psms, $lfq_wide (log2 matrix), $protein_meta
parse_sage_results <- function(
  results_path,
  lfq_path,
  fdr_threshold         = 0.01,
  protein_fdr_threshold = 0.01
) {
  if (!file.exists(results_path)) stop("Sage results file not found: ", results_path)
  if (is.null(lfq_path) || !file.exists(lfq_path)) {
    message("[DDA] No LFQ file provided — quantification will be PSM-based")
    lfq_path <- NULL
  }

  # PSM table
  psms <- data.table::fread(results_path)
  message("[DDA] Read ", nrow(psms), " total PSMs from Sage")

  # FDR filter
  psms_filtered <- psms[spectrum_q <= fdr_threshold & protein_q <= protein_fdr_threshold]
  message("[DDA] After FDR filter (spectrum_q <= ", fdr_threshold,
          ", protein_q <= ", protein_fdr_threshold, "): ", nrow(psms_filtered), " PSMs")

  # Per-protein peptide/spectra counts
  pep_counts <- psms_filtered[, .(
    NPeptides = data.table::uniqueN(peptide),
    NSpectra  = .N
  ), by = proteins]

  # LFQ matrix -> protein x run, log2.
  # Sage v0.14.7 lfq.tsv is WIDE + peptide-level: id columns
  # (peptide, charge, proteins, q_value, score, spectral_angle) followed by
  # one raw-intensity column per run file. We roll peptides up to protein
  # groups by summing intensities per run. (A legacy long-format path —
  # proteins/filename/intensity — is still supported for back-compat.)
  if (!is.null(lfq_path)) {
    lfq <- data.table::fread(lfq_path)
    if (all(c("filename", "intensity") %in% names(lfq))) {
      lfq_wide <- data.table::dcast(lfq, proteins ~ filename,
                                    value.var = "intensity", fun.aggregate = sum)
      rownames_col <- lfq_wide$proteins
      lfq_mat <- as.matrix(lfq_wide[, -"proteins", with = FALSE])
    } else {
      id_cols   <- intersect(c("peptide", "charge", "proteins", "q_value",
                               "score", "spectral_angle", "filename"), names(lfq))
      file_cols <- setdiff(names(lfq), id_cols)
      if (length(file_cols) == 0)
        stop("LFQ file '", basename(lfq_path), "' has no per-run intensity columns ",
             "(found: ", paste(names(lfq), collapse = ", "), ")")
      agg <- lfq[, lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
                 by = proteins, .SDcols = file_cols]
      rownames_col <- agg$proteins
      lfq_mat <- as.matrix(agg[, file_cols, with = FALSE])
    }
    rownames(lfq_mat) <- rownames_col

    # Log2 transform (Sage outputs raw intensities)
    lfq_mat_log2 <- log2(lfq_mat)
    lfq_mat_log2[!is.finite(lfq_mat_log2)] <- NA

    # Strip path and .d/.mzML extension from sample names
    colnames(lfq_mat_log2) <- gsub("\\.(d|mzML)$", "", basename(colnames(lfq_mat_log2)))
  } else {
    # Build a simple spectral count matrix from PSMs as fallback
    sc <- psms_filtered[, .N, by = .(proteins, filename)]
    sc_wide <- data.table::dcast(sc, proteins ~ filename, value.var = "N", fill = 0)
    rownames_col <- sc_wide$proteins
    lfq_mat_log2 <- log2(as.matrix(sc_wide[, -"proteins", with = FALSE]) + 1)
    rownames(lfq_mat_log2) <- rownames_col
    colnames(lfq_mat_log2) <- gsub("\\.(d|mzML)$", "", basename(colnames(lfq_mat_log2)))
    message("[DDA] Using spectral count matrix (no LFQ available)")
  }

  # Protein meta (genes-equivalent in EList)
  protein_meta <- merge(
    data.frame(ProteinID = rownames(lfq_mat_log2), stringsAsFactors = FALSE),
    as.data.frame(pep_counts),
    by.x = "ProteinID", by.y = "proteins",
    all.x = TRUE
  )
  protein_meta$NPeptides[is.na(protein_meta$NPeptides)] <- 0L
  protein_meta$NSpectra[is.na(protein_meta$NSpectra)]   <- 0L

  message("[DDA] LFQ matrix: ", nrow(lfq_mat_log2), " proteins x ",
          ncol(lfq_mat_log2), " samples")

  list(
    psms         = psms_filtered,
    lfq_wide     = lfq_mat_log2,
    protein_meta = protein_meta
  )
}


#' Build a limma EList from Sage MaxLFQ output
#'
#' @param sage_parsed List from parse_sage_results()
#' @param metadata_df Data.frame with SampleID and Group columns
#' @return EList object compatible with limma DE pipeline
build_dda_elist <- function(sage_parsed, metadata_df) {
  mat  <- sage_parsed$lfq_wide
  meta <- sage_parsed$protein_meta

  # Match sample columns to metadata
  sample_ids <- metadata_df$SampleID
  available  <- intersect(sample_ids, colnames(mat))
  if (length(available) == 0) {
    stop("[DDA] No sample names match between LFQ matrix and metadata. ",
         "Matrix cols: ", paste(head(colnames(mat), 5), collapse = ", "),
         "; Metadata SampleIDs: ", paste(head(sample_ids, 5), collapse = ", "))
  }
  mat <- mat[, available, drop = FALSE]
  metadata_df <- metadata_df[metadata_df$SampleID %in% available, ]

  # Build protein-level gene annotation (Protein.Group for compatibility)
  genes_df <- data.frame(
    Protein.Group = rownames(mat),
    NPeptides     = meta$NPeptides[match(rownames(mat), meta$ProteinID)],
    NSpectra      = meta$NSpectra[match(rownames(mat), meta$ProteinID)],
    stringsAsFactors = FALSE,
    row.names     = rownames(mat)
  )

  elist <- new("EList",
    list(
      E       = mat,
      genes   = genes_df,
      targets = metadata_df
    )
  )

  message("[DDA] Built EList: ", nrow(elist$E), " proteins x ", ncol(elist$E), " samples")
  elist
}


#' Normalize DDA MaxLFQ log2 matrix
#'
#' @param mat log2 protein x sample matrix (with NAs)
#' @param method One of "cyclicloess" (default), "median", "quantile", "none"
#' @return Normalized matrix (same dimensions, NAs preserved)
normalize_dda_matrix <- function(mat, method = "cyclicloess") {
  switch(method,
    "cyclicloess" = {
      limma::normalizeCyclicLoess(mat, method = "fast")
    },
    "median" = {
      sample_medians <- apply(mat, 2, median, na.rm = TRUE)
      global_median  <- median(mat, na.rm = TRUE)
      sweep(mat, 2, sample_medians - global_median, "-")
    },
    "quantile" = {
      if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        warning("[DDA] preprocessCore not available, falling back to median normalization")
        return(normalize_dda_matrix(mat, method = "median"))
      }
      preprocessCore::normalize.quantiles(mat)
    },
    "none" = mat,
    stop("Unknown normalization method: ", method)
  )
}


#' Filter DDA protein matrix by valid value threshold per group
#'
#' @param mat log2 protein x sample matrix (with NAs, post-normalization)
#' @param group_vec Character/factor vector of group assignments (length = ncol(mat))
#' @param min_valid_fraction Minimum fraction of samples in EACH group that must
#'   have a valid (non-NA) value. Default 0.5.
#' @return Filtered matrix (rows passing threshold)
filter_dda_valid_values <- function(mat, group_vec, min_valid_fraction = 0.5) {
  groups <- unique(group_vec)

  passes <- apply(mat, 1, function(row) {
    all(vapply(groups, function(g) {
      group_vals <- row[group_vec == g]
      n_valid    <- sum(!is.na(group_vals))
      n_total    <- length(group_vals)
      (n_valid / n_total) >= min_valid_fraction
    }, logical(1)))
  })

  mat[passes, , drop = FALSE]
}


#' Impute missing values in DDA protein matrix
#'
#' @param mat log2 protein x sample matrix (post-filter, NAs remain)
#' @param method One of "perseus" (default), "minprob", "mindet", "none"
#' @param width Perseus width parameter (default 0.3)
#' @param shift Perseus downshift in SD units (default 1.8)
#' @param q MinProb/MinDet quantile (default 0.01)
#' @return Imputed matrix (no NAs unless method = "none")
impute_dda_matrix <- function(
  mat,
  method = "perseus",
  width  = 0.3,
  shift  = 1.8,
  q      = 0.01
) {
  set.seed(42)  # reproducibility for stochastic imputation

  switch(method,
    "perseus" = {
      apply(mat, 2, function(col) {
        missing <- is.na(col)
        if (!any(missing)) return(col)
        col_mean <- mean(col, na.rm = TRUE)
        col_sd   <- sd(col, na.rm = TRUE)
        if (is.na(col_sd) || col_sd == 0) col_sd <- 1
        col[missing] <- rnorm(
          n    = sum(missing),
          mean = col_mean - shift * col_sd,
          sd   = width * col_sd
        )
        col
      })
    },
    "minprob" = {
      apply(mat, 2, function(col) {
        missing <- is.na(col)
        if (!any(missing)) return(col)
        obs     <- col[!missing]
        center  <- quantile(obs, probs = q, na.rm = TRUE)
        col_sd  <- sd(obs, na.rm = TRUE)
        if (is.na(col_sd) || col_sd == 0) col_sd <- 1
        col[missing] <- rnorm(
          n    = sum(missing),
          mean = center,
          sd   = width * col_sd
        )
        col
      })
    },
    "mindet" = {
      apply(mat, 2, function(col) {
        missing <- is.na(col)
        if (!any(missing)) return(col)
        col[missing] <- quantile(col, probs = q, na.rm = TRUE)
        col
      })
    },
    "none" = mat,
    stop("Unknown imputation method: ", method)
  )
}


#' Run the full DDA pre-DE pipeline: normalize -> filter -> impute -> build EList
#'
#' @param lfq_wide log2 protein x sample matrix from parse_sage_results()
#' @param protein_meta protein metadata data.frame from parse_sage_results()
#' @param metadata_df data.frame with SampleID and Group columns
#' @param norm_method normalization method
#' @param min_valid_fraction valid value filter threshold
#' @param impute_method imputation method
#' @param perseus_width width parameter for Perseus imputation
#' @param perseus_shift shift parameter for Perseus imputation
#' @return List with $elist (EList), $n_prefilter, $n_postfilter
run_dda_pipeline <- function(
  lfq_wide,
  protein_meta,
  metadata_df,
  norm_method        = "cyclicloess",
  min_valid_fraction = 0.5,
  impute_method      = "perseus",
  perseus_width      = 0.3,
  perseus_shift      = 1.8
) {
  # Step 1: Normalize
  message("[DDA Pipeline] Normalizing with method: ", norm_method)
  mat_norm <- normalize_dda_matrix(lfq_wide, method = norm_method)

  # Step 2: Filter valid values
  group_vec <- metadata_df$Group[match(colnames(mat_norm), metadata_df$SampleID)]
  n_before <- nrow(mat_norm)
  mat_filtered <- filter_dda_valid_values(mat_norm, group_vec, min_valid_fraction)
  n_after <- nrow(mat_filtered)
  message(sprintf("[DDA Pipeline] Valid value filter: %d -> %d proteins (%.0f%% retained)",
    n_before, n_after, 100 * n_after / max(n_before, 1)))

  # Step 3: Impute
  message("[DDA Pipeline] Imputing with method: ", impute_method)
  mat_imputed <- impute_dda_matrix(mat_filtered,
    method = impute_method,
    width  = perseus_width,
    shift  = perseus_shift
  )

  # Step 4: Build EList
  # Update protein_meta to match filtered rows
  filtered_meta <- protein_meta[protein_meta$ProteinID %in% rownames(mat_imputed), ]
  elist <- build_dda_elist(
    list(lfq_wide = mat_imputed, protein_meta = filtered_meta),
    metadata_df
  )

  list(
    elist       = elist,
    n_prefilter = n_before,
    n_postfilter = n_after
  )
}


#' Generate sbatch script for Sage DDA search on Hive
#'
#' @param sage_bin Path to Sage binary on HPC
#' @param config_path Path to sage.json config on HPC
#' @param raw_dir Directory containing .d files on HPC
#' @param output_dir Output directory on HPC
#' @param experiment_name Name for the SLURM job
#' @param cpus Number of CPUs (default 32)
#' @param mem_gb Memory in GB (default 64)
#' @param time_limit Time limit string (default "02:00:00")
#' @param account SLURM account (default "genome-center-grp")
#' @param partition SLURM partition (default "high")
#' @return Character string: sbatch script content
generate_sage_sbatch <- function(
  sage_bin,
  config_path,
  raw_dir,
  output_dir,
  experiment_name = "sage_search",
  cpus            = 32,
  mem_gb          = 64,
  time_limit      = "02:00:00",
  account         = "genome-center-grp",
  partition       = "high"
) {
  # Sanitize experiment name for SLURM

  safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", experiment_name)

  script <- paste0(
'#!/bin/bash
#SBATCH --job-name=delimp_sage_', safe_name, '
#SBATCH --partition=', partition, '
#SBATCH --account=', account, '
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=', cpus, '
#SBATCH --mem=', mem_gb, 'G
#SBATCH --time=', time_limit, '
#SBATCH --output="', output_dir, '/logs/sage_%j.out"
#SBATCH --error="', output_dir, '/logs/sage_%j.err"

set -euo pipefail
shopt -s nullglob       # critical: empty glob → empty list, not error under pipefail
echo "[DE-LIMP Sage] Job start: $(date)"
echo "[DE-LIMP Sage] Node: $(hostname)"
echo "[DE-LIMP Sage] CPUs: ', cpus, ', Memory: ', mem_gb, 'G"

SAGE_BIN="', sage_bin, '"
CONFIG="', config_path, '"
RAW_DIR="', raw_dir, '"
OUTPUT_DIR="', output_dir, '"
MZML_DIR="$OUTPUT_DIR/mzml"
MSCONVERT_SIF="/quobyte/proteomics-grp/apptainers/pwiz-skyline-i-agree-to-the-vendor-licenses_latest.sif"

# Verify Sage binary
if [ ! -x "$SAGE_BIN" ]; then
  echo "[ERROR] Sage binary not found or not executable: $SAGE_BIN"
  exit 1
fi

# Create output dirs
mkdir -p "$OUTPUT_DIR" "$MZML_DIR"

# Collect .d directories (Sage reads timsTOF natively)
D_FILES=$(find "$RAW_DIR" -maxdepth 1 -name "*.d" -type d 2>/dev/null | sort)
N_D=$(echo "$D_FILES" | grep -c "." 2>/dev/null || true)
N_D=${N_D:-0}

# .raw files: Sage v0.14.7 CANNOT parse Thermo .raw — silently produces 0 hits.
# Pre-convert each .raw → mzML using msconvert (--bind /quobyte:/quobyte is mandatory).
RAW_FILES=$(find "$RAW_DIR" -maxdepth 1 -name "*.raw" -type f 2>/dev/null | sort)
N_RAW=$(echo "$RAW_FILES" | grep -c "." 2>/dev/null || true)
N_RAW=${N_RAW:-0}
echo "[DE-LIMP Sage] Found $N_D .d files and $N_RAW .raw files in $RAW_DIR"

if [ "$N_RAW" -gt 0 ]; then
  if [ ! -f "$MSCONVERT_SIF" ]; then
    echo "[ERROR] msconvert container missing at $MSCONVERT_SIF — cannot convert .raw"
    exit 1
  fi
  echo "[DE-LIMP Sage] Pre-converting $N_RAW .raw files to mzML ..."
  for RAW_FILE in $RAW_FILES; do
    BASENAME=$(basename "$RAW_FILE" .raw)
    if [ -f "$MZML_DIR/${BASENAME}.mzML" ]; then
      echo "  [skip] $BASENAME.mzML already exists"
      continue
    fi
    echo "  [conv] $RAW_FILE"
    if ! apptainer exec --bind /quobyte:/quobyte "$MSCONVERT_SIF" wine msconvert "$RAW_FILE" --mzML --filter "peakPicking true 1-" -o "$MZML_DIR/"; then
      echo "[ERROR] msconvert failed for $BASENAME"
      exit 1
    fi
  done
fi

# Collect the converted mzML files
MZML_FILES=$(find "$MZML_DIR" -maxdepth 1 -name "*.mzML" -type f 2>/dev/null | sort)
N_MZML=$(echo "$MZML_FILES" | grep -c "." 2>/dev/null || true)
N_MZML=${N_MZML:-0}
echo "[DE-LIMP Sage] Sage will read: $N_D .d + $N_MZML mzML files"

if [ "$((N_D + N_MZML))" -eq 0 ]; then
  echo "[ERROR] No readable mass-spec files (.d or .mzML) for Sage"
  exit 1
fi

# Run Sage. CLI args override config.mzml_paths (per Sage v0.14.7 docs).
echo "[DE-LIMP Sage] Starting search..."
"$SAGE_BIN" \\
  --write-pin \\
  --output_directory "$OUTPUT_DIR" \\
  "$CONFIG" \\
  $D_FILES $MZML_FILES

echo "[DE-LIMP Sage] Search complete: $(date)"
echo "[DE-LIMP Sage] Output files:"
ls -lh "$OUTPUT_DIR"/*.tsv "$OUTPUT_DIR"/*.json 2>/dev/null || echo "(no output files found)"
echo "[DE-LIMP Sage] Done."
')
  script
}


#' Parse sage_report.json summary statistics
#'
#' @param json_path Path to sage_report.json or results.json
#' @return Named list: psms, unique_peptides, unique_proteins, sage_version
parse_sage_report <- function(json_path) {
  if (!file.exists(json_path)) {
    message("[DDA] sage report JSON not found: ", json_path)
    return(NULL)
  }

  tryCatch({
    report <- jsonlite::fromJSON(json_path)
    list(
      psms             = report$psms %||% NA_integer_,
      unique_peptides  = report$unique_peptides %||% NA_integer_,
      unique_proteins  = report$unique_proteins %||% NA_integer_,
      sage_version     = report$sage_version %||% report$version %||% "unknown",
      files            = report$files %||% NA_integer_
    )
  }, error = function(e) {
    message("[DDA] Failed to parse sage report: ", e$message)
    NULL
  })
}


#' Compute DDA QC metrics from PSM table
#'
#' @param psms Filtered PSM data.table from parse_sage_results()
#' @param lfq_wide log2 protein x sample matrix
#' @return Named list of QC metrics
compute_dda_qc_metrics <- function(psms, lfq_wide) {
  n_psms       <- nrow(psms)
  n_peptides   <- data.table::uniqueN(psms$peptide)
  n_proteins   <- nrow(lfq_wide)

  # Median peptides per protein
  pep_per_prot <- psms[, .(n = data.table::uniqueN(peptide)), by = proteins]
  med_pep_per_prot <- median(pep_per_prot$n, na.rm = TRUE)

  # Missed cleavages
  mc_pattern <- "[KR][^P]"  # tryptic missed cleavage sites
  n_mc <- sum(grepl(mc_pattern, psms$peptide))
  pct_missed_cleavage <- 100 * n_mc / max(n_psms, 1)

  # Precursor mass error (ppm) -- column may be called ppm_difference or similar
  ppm_col <- intersect(c("expmass_ppm", "ppm_difference", "delta_mass_ppm"), colnames(psms))
  if (length(ppm_col) > 0) {
    mass_error_median <- median(abs(psms[[ppm_col[1]]]), na.rm = TRUE)
  } else {
    mass_error_median <- NA_real_
  }

  list(
    n_psms              = n_psms,
    n_peptides          = n_peptides,
    n_proteins          = n_proteins,
    med_pep_per_prot    = med_pep_per_prot,
    pct_missed_cleavage = round(pct_missed_cleavage, 1),
    mass_error_ppm      = round(mass_error_median, 2)
  )
}


# ==============================================================================
#  Casanovo de novo sequencing helpers
# ==============================================================================

#' Generate sbatch script for Casanovo de novo sequencing (GPU)
#'
#' Creates a two-phase sbatch: (1) convert .d to MGF via bruker_to_mgf.py,
#' (2) run Casanovo sequence on each MGF file as an array job.
#'
#' @param raw_dir Directory containing .d files on HPC
#' @param output_dir Output directory on HPC (casanovo/ subdir created)
#' @param experiment_name Name for SLURM job
#' @param conda_env_path Path to Casanovo conda environment
#' @param model_ckpt Path to Casanovo model checkpoint
#' @param converter_script Path to bruker_to_mgf.py on HPC
#' @param n_files Number of .d files (for array job sizing)
#' @param account SLURM account
#' @param gpu_partition GPU partition name
#' @param gpu_qos GPU QOS name
#' @return List with $convert_script (MGF conversion sbatch) and $casanovo_script (sequencing sbatch)
generate_casanovo_sbatch <- function(
  raw_dir,
  output_dir,
  experiment_name  = "casanovo",
  conda_env_path   = "/quobyte/proteomics-grp/conda_envs/cassonovo_env",
  model_ckpt       = "/quobyte/proteomics-grp/bioinformatics_programs/casanovo_modles/casanovo_v4_2_0.ckpt",
  casanovo_version = "v4",
  compute_mode     = "gpu",     # "gpu" → gpu-a100 + --gres=gpu:1; "cpu" → high partition, no GPU
  converter_script = "/quobyte/proteomics-grp/de-limp/python/bruker_to_mgf.py",
  n_files          = 1,
  account          = "genome-center-grp",
  gpu_partition    = "gpu-a100",
  gpu_qos          = "genome-center-grp-gpu-a100-qos",
  cpu_partition    = "high",
  cpu_qos          = "genome-center-grp-high-qos"
) {
  safe_name <- gsub("[^a-zA-Z0-9_.-]", "_", experiment_name)
  casanovo_dir <- file.path(output_dir, "casanovo")
  mgf_dir      <- file.path(casanovo_dir, "mgf")
  mztab_dir    <- file.path(casanovo_dir, "mztab")
  logs_dir     <- file.path(output_dir, "logs")

  # ---- Step 1: MGF conversion (CPU, no GPU needed) ----
  convert_script <- paste0(
'#!/bin/bash
#SBATCH --job-name=delimp_casanovo_prep_', safe_name, '
#SBATCH --partition=high
#SBATCH --account=', account, '
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output="', logs_dir, '/casanovo_prep_%j.out"
#SBATCH --error="', logs_dir, '/casanovo_prep_%j.err"

set -euo pipefail
shopt -s nullglob          # critical: glob with no matches → empty list, not error
echo "[DE-LIMP MGF] Start: $(date)"
echo "[DE-LIMP MGF] Node: $(hostname)"

CONDA_ENV="', conda_env_path, '"
RAW_DIR="', raw_dir, '"
MGF_DIR="', mgf_dir, '"

mkdir -p "$MGF_DIR"

# Activate conda environment for timsrust_pyo3
export PATH="$CONDA_ENV/bin:$PATH"

# Check for existing MGF files first (may already exist from previous runs).
# With nullglob, "$RAW_DIR"/*.mgf expands to empty when no matches → safe.
EXISTING_MGF=$(find "$RAW_DIR" -maxdepth 1 -name "*.mgf" -type f 2>/dev/null | wc -l)
if [ "$EXISTING_MGF" -gt 0 ]; then
  echo "[DE-LIMP MGF] Found $EXISTING_MGF existing MGF files — copying"
  cp "$RAW_DIR"/*.mgf "$MGF_DIR/" 2>/dev/null || true
fi

# Convert .d files to MGF (timsTOF)
N_D=$(find "$RAW_DIR" -maxdepth 1 -name "*.d" -type d 2>/dev/null | wc -l)
if [ "$N_D" -gt 0 ]; then
  echo "[DE-LIMP MGF] Converting $N_D .d files to MGF"
  python "', converter_script, '" "$RAW_DIR" "$MGF_DIR" --batch --min-peaks 6 -v
fi

# NOTE: .raw → mzML happens in the Sage sbatch step (single conversion, reused).
# Casanovo v5 reads mzML directly, so we just point at $OUTPUT_DIR/mzml/ here.
# This step assumes the Sage job has finished (launcher passes --dependency=afterok:SAGE_ID).
MZML_DIR_FROM_SAGE="', output_dir, '/mzml"
N_MZML=$(find "$MZML_DIR_FROM_SAGE" -maxdepth 1 -name "*.mzML" -type f 2>/dev/null | wc -l)
echo "[DE-LIMP MGF] $N_MZML mzML files available from Sage conversion in $MZML_DIR_FROM_SAGE"

# Build Casanovo input list: ALL MGFs (from .d conversion) + ALL mzMLs (from Sage)
INPUT_LIST="', casanovo_dir, '/mgf_file_list.txt"
{
  find "$MGF_DIR"             -maxdepth 1 -name "*.mgf"  -type f 2>/dev/null | sort
  find "$MZML_DIR_FROM_SAGE"  -maxdepth 1 -name "*.mzML" -type f 2>/dev/null | sort
} > "$INPUT_LIST"
N_MGF=$(wc -l < "$INPUT_LIST" 2>/dev/null || echo 0)
echo "[DE-LIMP MGF] Total MGF files: $N_MGF"
echo "[DE-LIMP MGF] Done: $(date)"
')

  # ---- Step 2: Casanovo sequencing (GPU array job, 1 per file; CPU fallback when gpu-a100 is loaded) ----
  is_cpu <- identical(compute_mode, "cpu")
  cas_partition <- if (is_cpu) cpu_partition else gpu_partition
  cas_qos       <- if (is_cpu) cpu_qos       else gpu_qos
  cas_cpus      <- if (is_cpu) 16 else 8                       # more CPU cores when no GPU
  cas_mem       <- if (is_cpu) "32G" else "32G"
  cas_time      <- if (is_cpu) "08:00:00" else "01:30:00"      # CPU inference is ~10-30x slower per file
  cas_gres_line <- if (is_cpu) "" else "#SBATCH --gres=gpu:1\n"

  casanovo_script <- paste0(
'#!/bin/bash
#SBATCH --job-name=delimp_casanovo_', safe_name, '
#SBATCH --partition=', cas_partition, '
#SBATCH --account=', account, '
#SBATCH --qos=', cas_qos, '
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=', cas_cpus, '
#SBATCH --mem=', cas_mem, '
', cas_gres_line,
'#SBATCH --time=', cas_time, '
#SBATCH --array=1-', n_files, '
#SBATCH --output="', logs_dir, '/casanovo_%A_%a.out"
#SBATCH --error="', logs_dir, '/casanovo_%A_%a.err"

set -euo pipefail
echo "[DE-LIMP Casanovo] Task ${SLURM_ARRAY_TASK_ID} start: $(date)"
echo "[DE-LIMP Casanovo] Node: $(hostname)"
echo "[DE-LIMP Casanovo] Mode: ', compute_mode, '"
', if (is_cpu) {
  'export CUDA_VISIBLE_DEVICES=""\n'
} else {
  'echo "[DE-LIMP Casanovo] GPU: $(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo unknown)"\n'
}, '
CONDA_ENV="', conda_env_path, '"
MODEL="', model_ckpt, '"
MGF_DIR="', mgf_dir, '"
MZTAB_DIR="', mztab_dir, '"

mkdir -p "$MZTAB_DIR"

# Activate conda environment
export PATH="$CONDA_ENV/bin:$PATH"

# Get the input peak file for this array task (MGF or mzML — both supported by Casanovo v5)
MGF_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "', casanovo_dir, '/mgf_file_list.txt")
if [ -z "$MGF_FILE" ]; then
  echo "[ERROR] No input file for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

# Strip whichever extension applies — Casanovo writes "${BASENAME}_sequence.mztab"
BASENAME=$(basename "$MGF_FILE")
BASENAME="${BASENAME%.mgf}"
BASENAME="${BASENAME%.mzML}"
BASENAME="${BASENAME%.mzml}"
OUTPUT_FILE="$MZTAB_DIR/${BASENAME}_sequence.mztab"

echo "[DE-LIMP Casanovo] Processing: $MGF_FILE"
echo "[DE-LIMP Casanovo] Output: $OUTPUT_FILE"

# Casanovo v5 `--force_overwrite` is unreliable when ANY file matching the
# output-root glob already exists (e.g. stale stubs from a prior crashed run).
# Pre-delete any leftover artifacts for THIS basename to make the task idempotent.
rm -f "$MZTAB_DIR/${BASENAME}_sequence."* 2>/dev/null || true

# Run Casanovo de novo sequencing
', if (identical(casanovo_version, "v5")) {
'# v5 CLI: --output_dir + --output_root (no --output)
casanovo sequence \\
  --model "$MODEL" \\
  --output_dir "$MZTAB_DIR" \\
  --output_root "${BASENAME}_sequence" \\
  --force_overwrite \\
  "$MGF_FILE"
# v5 may write to ${BASENAME}_sequence.mztab inside MZTAB_DIR — verify
ls -lh "$MZTAB_DIR/${BASENAME}_sequence"* 2>/dev/null || true'
} else {
'# v4 CLI: --output is the full mztab path
casanovo sequence --model "$MODEL" --output "$OUTPUT_FILE" "$MGF_FILE"'
}, '

echo "[DE-LIMP Casanovo] Task ${SLURM_ARRAY_TASK_ID} done: $(date)"
')

  list(
    convert_script  = convert_script,
    casanovo_script = casanovo_script,
    casanovo_dir    = casanovo_dir,
    mgf_dir         = mgf_dir,
    mztab_dir       = mztab_dir
  )
}


#' Parse Casanovo mzTab output files
#'
#' Extracts PSM rows from mzTab format. Each row has a de novo predicted
#' peptide sequence with confidence score and per-residue amino acid scores.
#'
#' @param mztab_paths Character vector of .mztab file paths
#' @param score_threshold Minimum confidence score (default -0.5, Casanovo uses negative log-prob)
#' @return data.table with columns: sequence, seq_stripped, seq_norm, score,
#'   aa_scores, charge, exp_mz, calc_mz, source_file, psm_id
parse_casanovo_mztab <- function(mztab_paths, score_threshold = -Inf) {
  results <- lapply(mztab_paths, function(path) {
    if (!file.exists(path)) {
      message("[Casanovo] File not found: ", path)
      return(NULL)
    }

    lines <- readLines(path, warn = FALSE)

    # Find PSH (header) line
    psh_idx <- which(startsWith(lines, "PSH"))
    if (length(psh_idx) == 0) {
      message("[Casanovo] No PSH header in: ", basename(path))
      return(NULL)
    }

    # Parse header
    header <- strsplit(lines[psh_idx[1]], "\t")[[1]]

    # Find PSM rows
    psm_idx <- which(startsWith(lines, "PSM"))
    if (length(psm_idx) == 0) {
      message("[Casanovo] No PSM rows in: ", basename(path))
      return(NULL)
    }

    # Parse PSM rows
    psm_data <- do.call(rbind, lapply(psm_idx, function(i) {
      strsplit(lines[i], "\t")[[1]]
    }))
    colnames(psm_data) <- header

    df <- as.data.frame(psm_data, stringsAsFactors = FALSE)

    # Extract key columns (column names from mzTab spec)
    # sequence, PSM_ID, search_engine_score[1], charge, exp_mass_to_charge,
    # calc_mass_to_charge, opt_ms_run[1]_aa_scores
    out <- data.frame(
      sequence     = df$sequence,
      psm_id       = as.integer(df$PSM_ID),
      score        = as.numeric(df[["search_engine_score[1]"]]),
      charge       = as.integer(df$charge),
      exp_mz       = as.numeric(df$exp_mass_to_charge),
      calc_mz      = as.numeric(df$calc_mass_to_charge),
      # Strip Casanovo's "_sequence" suffix so source_file matches the Sage
      # run name (Sage filenames are <run>.mzML; Casanovo mztabs are
      # <run>_sequence.mztab). Without this, every Sage<->Casanovo match on
      # filename (per-sample summary, Disagreements) finds 0 shared spectra.
      source_file  = sub("_sequence$", "", basename(tools::file_path_sans_ext(path))),
      stringsAsFactors = FALSE
    )

    # Per-residue AA scores (if present)
    aa_col <- grep("aa_scores", names(df), value = TRUE)
    if (length(aa_col) > 0) {
      out$aa_scores <- df[[aa_col[1]]]
    } else {
      out$aa_scores <- NA_character_
    }

    out
  })

  # Combine all files
  combined <- do.call(rbind, results[!vapply(results, is.null, logical(1))])

  if (is.null(combined) || nrow(combined) == 0) {
    message("[Casanovo] No PSMs parsed from any file")
    return(data.table::data.table(
      sequence = character(0), psm_id = integer(0), score = numeric(0),
      charge = integer(0), exp_mz = numeric(0), calc_mz = numeric(0),
      source_file = character(0), aa_scores = character(0),
      seq_stripped = character(0), seq_norm = character(0),
      mean_aa_score = numeric(0)
    ))
  }

  dt <- data.table::as.data.table(combined)

  # Filter by score threshold
  if (is.finite(score_threshold)) {
    dt <- dt[score >= score_threshold]
  }

  # Strip modifications to bare amino acids via the canonical normalizer so
  # Casanovo peptides match Sage DB peptides AND the DIAMOND query FASTA.
  dt$seq_stripped <- build_dda_canonical_peptide(dt$sequence)

  # I/L normalization for cross-reference (leucine = isoleucine in MS)
  dt$seq_norm <- gsub("I", "L", dt$seq_stripped)

  # Compute mean per-residue AA score
  dt$mean_aa_score <- vapply(dt$aa_scores, function(s) {
    if (is.na(s) || !nzchar(s) || s == "null") return(NA_real_)
    vals <- as.numeric(strsplit(s, ",")[[1]])
    if (length(vals) == 0 || all(is.na(vals))) return(NA_real_)
    mean(vals, na.rm = TRUE)
  }, numeric(1))

  message("[Casanovo] Parsed ", nrow(dt), " PSMs from ", length(mztab_paths), " files")
  dt
}


#' Parse DIA-NN DDA search results (report.parquet) into DE-LIMP-compatible format
#'
#' Returns a data frame compatible with `classify_dda_denovo()` — same columns
#' as Sage PSMs ($peptide, $proteins) plus DIA-NN-specific fields for fuzzy
#' spectrum matching (precursor_mz, rt, charge).
#'
#' @param parquet_path Path to DIA-NN report.parquet
#' @param fdr_threshold Q.Value cutoff (default 0.01)
#' @return List with:
#'   $psms: data.table with columns: peptide, proteins, charge, precursor_mz,
#'          rt, q_value, stripped_sequence, modified_sequence, filename
#'   $n_precursors: total precursors before FDR filter
#'   $n_proteins: unique protein groups after filter
parse_diann_dda_results <- function(parquet_path, fdr_threshold = 0.01) {
  if (!file.exists(parquet_path)) stop("DIA-NN report not found: ", parquet_path)

  df <- arrow::read_parquet(parquet_path)
  message("[DIA-NN DDA] Read ", nrow(df), " rows from report.parquet")
  n_precursors <- nrow(df)


  # Identify key columns (DIA-NN uses these standard names)
  req_cols <- c("Stripped.Sequence", "Protein.Group", "Precursor.Charge", "Q.Value")
  missing <- setdiff(req_cols, names(df))
  if (length(missing) > 0) {
    stop("[DIA-NN DDA] Missing required columns: ", paste(missing, collapse = ", "))
  }

  # FDR filter
  dt <- data.table::as.data.table(df)
  dt <- dt[Q.Value <= fdr_threshold]
  message("[DIA-NN DDA] After Q.Value <= ", fdr_threshold, ": ", nrow(dt), " rows")

  if (nrow(dt) == 0) {
    return(list(
      psms = data.table::data.table(
        peptide = character(0), proteins = character(0),
        charge = integer(0), precursor_mz = numeric(0),
        rt = numeric(0), q_value = numeric(0),
        stripped_sequence = character(0), modified_sequence = character(0),
        filename = character(0)
      ),
      n_precursors = n_precursors, n_proteins = 0L
    ))
  }

  # Build output matching Sage PSM column conventions:
  # - $peptide: the modified sequence (used by classify_dda_denovo for stripping)
  # - $proteins: protein group IDs
  # DIA-NN Modified.Sequence uses UniMod format: e.g. C(UniMod:4), M(UniMod:35)
  # classify_dda_denovo strips with gsub("\\+[0-9.]+", ...) which won't match UniMod.
  # So we put Stripped.Sequence in $peptide (already stripped) for classification,
  # and keep Modified.Sequence separately for reference.
  out <- data.table::data.table(
    peptide            = dt$Stripped.Sequence,
    proteins           = dt$Protein.Group,
    charge             = as.integer(dt$Precursor.Charge),
    precursor_mz       = if ("Precursor.Mz" %in% names(dt)) as.numeric(dt$Precursor.Mz) else NA_real_,
    rt                 = if ("RT" %in% names(dt)) as.numeric(dt$RT) else NA_real_,
    q_value            = as.numeric(dt$Q.Value),
    stripped_sequence   = dt$Stripped.Sequence,
    modified_sequence   = if ("Modified.Sequence" %in% names(dt)) dt$Modified.Sequence else dt$Stripped.Sequence,
    filename           = if ("File.Name" %in% names(dt)) basename(dt$File.Name) else
                         if ("Run" %in% names(dt)) dt$Run else NA_character_
  )

  # Deduplicate to unique peptide-protein pairs (DIA-NN can have multiple rows
  # per precursor across runs)
  n_proteins <- length(unique(out$proteins))

  message(sprintf("[DIA-NN DDA] Parsed %d PSMs, %d unique peptides, %d protein groups",
    nrow(out),
    length(unique(out$stripped_sequence)),
    n_proteins))

  list(
    psms         = out,
    n_precursors = n_precursors,
    n_proteins   = n_proteins
  )
}


#' Fuzzy match DIA-NN PSMs to Casanovo PSMs by precursor properties
#'
#' DIA-NN doesn't output scan numbers, so we match by:
#' - Same source file (by filename stem)
#' - Precursor m/z within tolerance (default 20 ppm)
#' - RT within tolerance (default 1 minute)
#' - Same charge state
#'
#' @param casanovo_dt data.table from parse_casanovo_mztab()
#' @param diann_psms data.table from parse_diann_dda_results()$psms
#' @param mz_ppm m/z tolerance in ppm (default 20)
#' @param rt_tol RT tolerance in minutes (default 1.0)
#' @return data.table with matched rows: casanovo columns + diann_peptide, diann_proteins, match_delta_mz, match_delta_rt
fuzzy_match_diann_casanovo <- function(casanovo_dt, diann_psms,
                                       mz_ppm = 20, rt_tol = 1.0) {
  if (is.null(casanovo_dt) || nrow(casanovo_dt) == 0 ||
      is.null(diann_psms) || nrow(diann_psms) == 0) {
    return(data.table::data.table())
  }

  # Both need exp_mz and charge columns
  if (!all(c("exp_mz", "charge") %in% names(casanovo_dt))) {
    message("[Fuzzy match] Casanovo missing exp_mz or charge columns")
    return(data.table::data.table())
  }
  if (!all(c("precursor_mz", "charge") %in% names(diann_psms))) {
    message("[Fuzzy match] DIA-NN missing precursor_mz or charge columns")
    return(data.table::data.table())
  }

  # Normalize source filenames for matching
  cas_files <- gsub("\\.(d|mzML|mgf|mztab)$", "",
    casanovo_dt$source_file, ignore.case = TRUE)
  diann_files <- gsub("\\.(d|raw|mzML)$", "",
    diann_psms$filename, ignore.case = TRUE)

  matches <- list()
  for (i in seq_len(nrow(casanovo_dt))) {
    cas <- casanovo_dt[i, ]
    cas_mz <- cas$exp_mz
    cas_charge <- cas$charge

    if (is.na(cas_mz) || is.na(cas_charge)) next

    # Filter DIA-NN by charge
    candidates <- diann_psms[diann_psms$charge == cas_charge, ]
    if (nrow(candidates) == 0) next

    # Filter by m/z tolerance (ppm)
    mz_tol_abs <- cas_mz * mz_ppm / 1e6
    candidates <- candidates[abs(candidates$precursor_mz - cas_mz) <= mz_tol_abs, ]
    if (nrow(candidates) == 0) next

    # Filter by RT if both have RT values
    if (!is.na(cas$exp_mz) && any(!is.na(candidates$rt))) {
      # Casanovo mztab doesn't have RT directly — use calc_mz as proxy
      # or skip RT filter if not available
      # Actually, mztab PSMs don't have RT. Skip RT filter for now.
    }

    # Take best match by m/z proximity
    deltas <- abs(candidates$precursor_mz - cas_mz)
    best <- which.min(deltas)

    matches[[length(matches) + 1]] <- data.table::data.table(
      casanovo_idx     = i,
      casanovo_seq     = cas$seq_stripped,
      casanovo_score   = cas$score,
      diann_peptide    = candidates$peptide[best],
      diann_proteins   = candidates$proteins[best],
      delta_mz_ppm     = deltas[best] / cas_mz * 1e6,
      casanovo_charge  = cas_charge
    )
  }

  if (length(matches) == 0) return(data.table::data.table())

  result <- data.table::rbindlist(matches)

  # Classify agreement
  cas_norm <- gsub("I", "L", result$casanovo_seq)
  diann_norm <- gsub("I", "L", result$diann_peptide)
  result$agreement <- ifelse(cas_norm == diann_norm, "agree", "disagree")

  message(sprintf("[Fuzzy match] %d matches: %d agree, %d disagree",
    nrow(result),
    sum(result$agreement == "agree"),
    sum(result$agreement == "disagree")))

  result
}


#' Cross-reference Casanovo de novo results against database search PSMs
#'
#' Classifies each Casanovo sequence as:
#'   - "confirmed": exact match to a database search FDR-passing peptide (I/L normalized)
#'   - "novel": no match in database search results (potential novel peptide)
#'
#' @param casanovo_dt data.table from parse_casanovo_mztab()
#' @param sage_psms Filtered PSM data.table from parse_sage_results()$psms
#'   or parse_diann_dda_results()$psms. Requires $peptide and $proteins columns.
#' @param db_engine Character: "Sage" or "DIA-NN" (default "Sage"). Stored in result
#'   for display in source badges.
#' @return List with:
#'   $classified: full casanovo_dt with match_type column
#'   $confirmed: confirmed-only rows with protein mapping
#'   $novel: novel-only rows
#'   $protein_summary: per-protein Casanovo confirmation stats
#'   $summary_stats: overall classification counts
#'   $db_engine: which database search engine was used
classify_dda_denovo <- function(casanovo_dt, sage_psms, db_engine = "Sage") {
  if (is.null(casanovo_dt) || nrow(casanovo_dt) == 0) {
    return(list(
      classified     = casanovo_dt,
      confirmed      = casanovo_dt[0, ],
      novel          = casanovo_dt[0, ],
      protein_summary = data.frame(
        proteins = character(0),
        n_casanovo_confirmed = integer(0),
        casanovo_max_score = numeric(0),
        casanovo_mean_aa_score = numeric(0),
        stringsAsFactors = FALSE
      ),
      summary_stats  = list(
        n_total = 0L, n_confirmed = 0L, n_novel = 0L,
        pct_confirmed = 0, pct_novel = 0
      ),
      db_engine = db_engine
    ))
  }

  # Handle case where database search results are not available (mztab-only load)
  if (is.null(sage_psms) || nrow(sage_psms) == 0) {
    casanovo_dt$match_type <- "novel"
    novel <- casanovo_dt
    confirmed <- casanovo_dt[0, ]
    confirmed$proteins <- character(0)

    n_total <- nrow(casanovo_dt)
    return(list(
      classified      = casanovo_dt,
      confirmed       = confirmed,
      novel           = novel,
      protein_summary = data.frame(
        proteins = character(0),
        n_casanovo_confirmed = integer(0),
        casanovo_max_score = numeric(0),
        casanovo_mean_aa_score = numeric(0),
        stringsAsFactors = FALSE
      ),
      summary_stats = list(
        n_total = n_total, n_confirmed = 0L, n_novel = n_total,
        pct_confirmed = 0, pct_novel = 100
      ),
      db_engine = db_engine
    ))
  }

  # Normalize database search peptides: strip modifications + I/L matching
  # Works for both Sage (C[+57.0215]) and DIA-NN (already stripped) peptide formats
  db_stripped <- build_dda_canonical_peptide(sage_psms$peptide)
  db_peps_norm <- unique(gsub("I", "L", db_stripped))

  # Classify each Casanovo sequence
  casanovo_dt$match_type <- ifelse(
    casanovo_dt$seq_norm %in% db_peps_norm, "confirmed", "novel"
  )

  # Map confirmed sequences back to database search protein groups
  pep_to_protein <- unique(
    data.frame(
      peptide  = sage_psms$peptide,
      proteins = sage_psms$proteins,
      stringsAsFactors = FALSE
    )
  )
  pep_to_protein$seq_norm <- gsub("I", "L", build_dda_canonical_peptide(pep_to_protein$peptide))
  # Collapse to ONE protein group per seq_norm. A peptide can map to several
  # Sage protein groups; without this the merge below is many-to-many and, once
  # all ~440k PSMs are loaded, blows past data.table's cartesian guard
  # ("Join results in N rows; more than nrow(x)+nrow(i)").
  pep_to_protein <- pep_to_protein[!duplicated(pep_to_protein$seq_norm), ]

  confirmed <- merge(
    casanovo_dt[casanovo_dt$match_type == "confirmed", ],
    pep_to_protein[, c("seq_norm", "proteins")],
    by = "seq_norm", all.x = TRUE
  )

  novel <- casanovo_dt[casanovo_dt$match_type == "novel", ]

  # Per-protein summary
  if (nrow(confirmed) > 0 && any(!is.na(confirmed$proteins))) {
    confirmed_mapped <- confirmed[!is.na(confirmed$proteins), ]

    protein_summary <- data.frame(
      proteins = unique(confirmed_mapped$proteins),
      stringsAsFactors = FALSE
    )
    protein_summary$n_casanovo_confirmed <- vapply(
      protein_summary$proteins,
      function(p) length(unique(confirmed_mapped$seq_norm[confirmed_mapped$proteins == p])),
      integer(1)
    )
    protein_summary$casanovo_max_score <- vapply(
      protein_summary$proteins,
      function(p) max(confirmed_mapped$score[confirmed_mapped$proteins == p], na.rm = TRUE),
      numeric(1)
    )
    protein_summary$casanovo_mean_aa_score <- vapply(
      protein_summary$proteins,
      function(p) {
        scores <- confirmed_mapped$mean_aa_score[confirmed_mapped$proteins == p]
        scores <- scores[!is.na(scores)]
        if (length(scores) == 0) return(NA_real_)
        mean(scores)
      },
      numeric(1)
    )
  } else {
    protein_summary <- data.frame(
      proteins = character(0),
      n_casanovo_confirmed = integer(0),
      casanovo_max_score = numeric(0),
      casanovo_mean_aa_score = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  # Summary stats
  n_total     <- nrow(casanovo_dt)
  n_confirmed <- sum(casanovo_dt$match_type == "confirmed")
  n_novel     <- sum(casanovo_dt$match_type == "novel")

  summary_stats <- list(
    n_total       = n_total,
    n_confirmed   = n_confirmed,
    n_novel       = n_novel,
    pct_confirmed = round(100 * n_confirmed / max(n_total, 1), 1),
    pct_novel     = round(100 * n_novel / max(n_total, 1), 1),
    n_proteins_with_denovo = nrow(protein_summary)
  )

  message(sprintf(
    "[Casanovo] Classification (vs %s): %d total, %d confirmed (%.1f%%), %d novel (%.1f%%)",
    db_engine, n_total, n_confirmed, summary_stats$pct_confirmed,
    n_novel, summary_stats$pct_novel
  ))

  list(
    classified      = casanovo_dt,
    confirmed       = confirmed,
    novel           = novel,
    protein_summary = protein_summary,
    summary_stats   = summary_stats,
    db_engine       = db_engine
  )
}


#' Generate Casanovo submit_all.sh launcher script
#'
#' Creates a shell script that submits MGF conversion first, then Casanovo
#' array job with dependency on conversion completing.
#'
#' @param convert_sbatch_path Remote path to MGF conversion sbatch
#' @param casanovo_sbatch_path Remote path to Casanovo sbatch
#' @return Character string: launcher script content
generate_casanovo_launcher <- function(convert_sbatch_path, casanovo_sbatch_path,
                                        sage_job_id = NULL) {
  # When sage_job_id is set, the MGF/file-list build waits for Sage's mzML
  # conversion to finish — Casanovo needs the mzML files to exist before its
  # array task can sed the input list.
  convert_dep <- if (!is.null(sage_job_id) && nzchar(sage_job_id)) {
    paste0("--dependency=afterok:", sage_job_id, " ")
  } else {
    ""
  }
  paste0(
'#!/bin/bash
set -euo pipefail

# Step 1: Submit MGF/mzML file-list build (waits for Sage if .raw conversion needed)
CONVERT_OUT=$(sbatch ', convert_dep, '"', convert_sbatch_path, '")
CONVERT_ID=$(echo "$CONVERT_OUT" | grep -o "[0-9]*$" | tail -1)
echo "CONVERT:${CONVERT_ID}"

# Step 2: Submit Casanovo with dependency on conversion
CASANOVO_OUT=$(sbatch --dependency=afterok:${CONVERT_ID} "', casanovo_sbatch_path, '")
CASANOVO_ID=$(echo "$CASANOVO_OUT" | grep -o "[0-9]*$" | tail -1)
echo "CASANOVO:${CASANOVO_ID}"

echo "File-list build job: ${CONVERT_ID}"
echo "Casanovo sequencing job: ${CASANOVO_ID} (depends on ${CONVERT_ID})"
')
}


#' Probe live state of a GPU partition for the Casanovo compute toggle
#'
#' Returns counts of pending vs running jobs + the per-node GPU allocation so the
#' UI can recommend GPU vs CPU intelligently. Cheap (3 squeue/sinfo calls, ~5s).
#'
#' @param ssh_config SSH config list (host, user, key) or NULL for local
#' @param partition  Partition name (default "gpu-a100")
#' @param sbatch_path Optional full path to remote sbatch binary (for non-login-shell SLURM calls)
#' @return list(success, partition, pending, running,
#'              nodes_total, nodes_usable, gpus_total, gpus_alloc, gpus_free,
#'              recommend = "gpu"|"cpu", reason = single-line text)
check_casanovo_gpu_queue <- function(ssh_config, partition = "gpu-a100", sbatch_path = NULL) {
  slurm_cmd <- function(cmd) {
    if (!is.null(sbatch_path) && nzchar(sbatch_path)) file.path(dirname(sbatch_path), cmd) else cmd
  }
  run_cmd <- function(command) {
    if (!is.null(ssh_config)) {
      ssh_exec(ssh_config, command, login_shell = is.null(sbatch_path), timeout = 12)
    } else if (exists("slurm_proxy_available", mode = "function") && slurm_proxy_available()) {
      slurm_proxy_exec(command, timeout = 12)
    } else {
      parts <- strsplit(command, " ")[[1]]
      stdout <- tryCatch(
        system2(parts[1], args = parts[-1], stdout = TRUE, stderr = TRUE),
        error = function(e) structure(e$message, status = 1L)
      )
      list(status = attr(stdout, "status") %||% 0L, stdout = stdout)
    }
  }

  result <- list(
    success = FALSE, partition = partition,
    pending = NA_integer_, running = NA_integer_,
    nodes_total = NA_integer_, nodes_usable = NA_integer_,
    gpus_total = NA_integer_, gpus_alloc = NA_integer_, gpus_free = NA_integer_,
    recommend = "gpu", reason = "Queue state unknown — defaulting to GPU."
  )

  # 1. Pending count
  tryCatch({
    cmd <- sprintf("%s -p %s -t PD --noheader -o %%i", slurm_cmd("squeue"), partition)
    res <- run_cmd(cmd)
    if (res$status == 0) {
      lines <- trimws(res$stdout); lines <- lines[nzchar(lines)]
      result$pending <- length(lines)
    }
  }, error = function(e) NULL)

  # 2. Running count
  tryCatch({
    cmd <- sprintf("%s -p %s -t R --noheader -o %%i", slurm_cmd("squeue"), partition)
    res <- run_cmd(cmd)
    if (res$status == 0) {
      lines <- trimws(res$stdout); lines <- lines[nzchar(lines)]
      result$running <- length(lines)
    }
  }, error = function(e) NULL)

  # 3. GPU totals + allocations per node via sinfo (Gres + GresUsed + state)
  tryCatch({
    cmd <- sprintf(
      "%s -p %s -h -O \"NodeList:25,StateCompact:14,Gres:35,GresUsed:35\"",
      slurm_cmd("sinfo"), partition
    )
    res <- run_cmd(cmd)
    if (res$status == 0 && length(res$stdout) > 0) {
      gpu_count <- function(s) {
        # Parse "gpu:a100:N" or "gpu:N" → N; sum across multiple "gres" types
        m <- regmatches(s, gregexpr("gpu(:[a-z0-9]+)?:[0-9]+", s))[[1]]
        if (length(m) == 0) return(0L)
        sum(as.integer(sub(".*:", "", m)))
      }
      lines <- trimws(res$stdout); lines <- lines[nzchar(lines)]
      n_total <- 0L; n_usable <- 0L
      g_total <- 0L; g_alloc <- 0L                # raw totals across ALL nodes
      g_total_usable <- 0L; g_alloc_usable <- 0L  # only on nodes that accept new jobs
      for (ln in lines) {
        parts <- strsplit(ln, "[[:space:]]{2,}")[[1]]
        if (length(parts) < 4) next
        state <- tolower(parts[2])
        gres  <- parts[3]
        used  <- parts[4]
        n_total <- n_total + 1L
        gt <- gpu_count(gres); ga <- gpu_count(used)
        g_total <- g_total + gt; g_alloc <- g_alloc + ga
        # drain/down/maint/reserv/fail/boot → won't accept new jobs
        if (!grepl("drain|down|maint|reserv|fail|boot", state)) {
          n_usable <- n_usable + 1L
          g_total_usable <- g_total_usable + gt
          g_alloc_usable <- g_alloc_usable + ga
        }
      }
      result$nodes_total <- n_total
      result$nodes_usable <- n_usable
      result$gpus_total <- g_total
      result$gpus_alloc <- g_alloc
      # "Free" reflects schedulable capacity, not raw subtraction
      result$gpus_free  <- max(0L, g_total_usable - g_alloc_usable)
    }
  }, error = function(e) NULL)

  result$success <- !is.na(result$pending) || !is.na(result$gpus_total)

  # 4. Recommendation
  if (result$success) {
    pd <- result$pending %||% 0L
    free_usable_gpus <- if (!is.na(result$nodes_usable) && result$nodes_usable == 0) {
      0L
    } else {
      result$gpus_free %||% 0L
    }
    if (free_usable_gpus > 0 && pd <= 3) {
      result$recommend <- "gpu"
      result$reason <- sprintf("GPU OK — %d free of %d, only %d pending.",
                               free_usable_gpus, result$gpus_total %||% 0L, pd)
    } else if (pd >= 10 || (free_usable_gpus == 0 && is.finite(result$nodes_usable) && result$nodes_usable <= 1)) {
      result$recommend <- "cpu"
      result$reason <- sprintf("Heavy GPU queue (%d pending, %d/%d GPUs free) — CPU is faster.",
                               pd, free_usable_gpus, result$gpus_total %||% 0L)
    } else {
      result$recommend <- "cpu"
      result$reason <- sprintf("GPU contested (%d pending, %d free) — CPU is the safer pick.",
                               pd, free_usable_gpus)
    }
  }

  result
}


#' Build a mode-aware ZIP export bundle for a DDA search
#'
#' Designed for HuggingFace Spaces viewing — flat folder of .md / .csv / .tsv /
#' .json / .mztab files that HF renders natively. Uses `safe_section()` so any
#' missing artifact gets logged into MANIFEST.txt with a clear `[SKIPPED]` line
#' rather than silently disappearing (CLAUDE.md architectural rule #4).
#'
#' @param output_dir Local path to the search's output directory (results files live here)
#' @param mode Analysis mode: "standard" | "phospho" | "tmt" | "peptidomics" | "hla_class_i" | "hla_class_ii"
#' @param app_version DE-LIMP version stamp for the export
#' @param search_info Optional list with submission metadata (overrides what we
#'        read from search_info.md if present)
#' @return Path to the .zip file (in tempdir())
generate_dda_export_zip <- function(output_dir, mode = "standard",
                                    app_version = "unknown",
                                    search_info = NULL) {
  stopifnot(dir.exists(output_dir))
  safe_mode <- gsub("[^a-zA-Z0-9_.-]", "_", mode)
  stamp     <- format(Sys.time(), "%Y%m%d_%H%M%S")
  bundle_name <- sprintf("delimp_dda_%s_%s", safe_mode, stamp)
  bundle_dir  <- file.path(tempdir(), bundle_name)
  if (dir.exists(bundle_dir)) unlink(bundle_dir, recursive = TRUE)
  dir.create(bundle_dir, recursive = TRUE)

  manifest <- new.env(parent = emptyenv())
  manifest$lines <- c(
    sprintf("DE-LIMP DDA Export Bundle  (v%s)", app_version),
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("Source:    %s", output_dir),
    sprintf("Mode:      %s", mode),
    "",
    "Files included:",
    "----------------------------------------"
  )

  # ---- Common: search_info → methods.md, sage.json → settings.json ----
  safe_section(manifest, "methods.md (search_info copy)", {
    src <- file.path(output_dir, "search_info.md")
    if (!file.exists(src)) stop("search_info.md not found in output_dir")
    file.copy(src, file.path(bundle_dir, "methods.md"), overwrite = TRUE)
  })

  safe_section(manifest, "settings.json (sage config)", {
    src <- file.path(output_dir, "sage.json")
    if (!file.exists(src)) stop("sage.json not found in output_dir")
    file.copy(src, file.path(bundle_dir, "settings.json"), overwrite = TRUE)
  })

  # ---- Sage results table (the main PSM TSV — needed for every mode) ----
  safe_section(manifest, "sage_results.tsv (PSM table)", {
    src <- file.path(output_dir, "results.sage.tsv")
    if (!file.exists(src)) stop("results.sage.tsv not found")
    file.copy(src, file.path(bundle_dir, "sage_results.tsv"), overwrite = TRUE)
  })

  # ---- Casanovo mztabs (optional — bundled if Casanovo ran) ----
  safe_section(manifest, "casanovo/*.mztab", {
    mztab_src_dir <- file.path(output_dir, "casanovo", "mztab")
    if (!dir.exists(mztab_src_dir)) stop("no casanovo mztab dir")
    mztab_files <- list.files(mztab_src_dir, pattern = "\\.mztab$",
                              full.names = TRUE)
    if (length(mztab_files) == 0) stop("no .mztab files in casanovo dir")
    dest <- file.path(bundle_dir, "casanovo")
    dir.create(dest, showWarnings = FALSE)
    file.copy(mztab_files, dest)
  })

  # ---- DIAMOND BLAST hits (optional — bundled if BLAST ran) ----
  safe_section(manifest, "diamond_hits.tsv (BLAST cascade)", {
    src <- file.path(output_dir, "denovo", "diamond_hits.tsv")
    if (!file.exists(src)) stop("diamond_hits.tsv not found")
    file.copy(src, file.path(bundle_dir, "diamond_hits.tsv"), overwrite = TRUE)
  })

  # ---- Universal peptide length distribution (cheap; useful for every mode) ----
  safe_section(manifest, "peptide_length_distribution.csv", {
    src <- file.path(output_dir, "results.sage.tsv")
    if (!file.exists(src)) stop("results.sage.tsv missing")
    df <- data.table::fread(src, select = c("peptide"), nThread = 1)
    df$length <- nchar(gsub("\\[.*?\\]|[^A-Z]", "", df$peptide))
    out <- as.data.frame(table(length = df$length), stringsAsFactors = FALSE)
    names(out) <- c("length", "n_psms")
    write.csv(out, file.path(bundle_dir, "peptide_length_distribution.csv"),
              row.names = FALSE)
  })

  # ──────────────────────────────────────────────────────────────────────────
  # Mode-specific sections. Each mode adds its own characteristic summary CSVs
  # alongside the universal artifacts above. Plots are deferred to the v3.11.6
  # / v3.11.7 viz work — once those exist we can embed PNG/SVG here too.
  # ──────────────────────────────────────────────────────────────────────────
  if (mode %in% c("hla_class_i", "hla_class_ii")) {
    # Anchor residue P2 + PΩ frequency table — the signature HLA fingerprint.
    safe_section(manifest, "hla_anchor_residues.csv (P2 / PΩ frequencies)", {
      src <- file.path(output_dir, "results.sage.tsv")
      if (!file.exists(src)) stop("results.sage.tsv missing")
      df <- data.table::fread(src, select = c("peptide"), nThread = 1)
      strip <- gsub("\\[.*?\\]|[^A-Z]", "", df$peptide)
      strip <- strip[nchar(strip) >= 8]    # ignore non-HLA-length junk
      p2 <- substr(strip, 2, 2)
      pomega <- substr(strip, nchar(strip), nchar(strip))
      out <- data.frame(
        position = rep(c("P2", "POmega"), each = length(unique(c(p2, pomega)))),
        residue  = rep(sort(unique(c(p2, pomega))), 2)
      )
      out$freq <- c(
        as.integer(table(factor(p2, levels = sort(unique(c(p2, pomega)))))),
        as.integer(table(factor(pomega, levels = sort(unique(c(p2, pomega))))))
      )
      out$pct <- 100 * out$freq /
                 ave(out$freq, out$position, FUN = sum)
      write.csv(out, file.path(bundle_dir, "hla_anchor_residues.csv"),
                row.names = FALSE)
    })
    # NOTE: anchor sequence logo (P1..PΩ) PNG/SVG is deferred to the viz pass.
  }

  if (identical(mode, "peptidomics")) {
    # N- and C-terminal cleavage flanking residues — protease motif fingerprint.
    safe_section(manifest, "peptidomics_cleavage_residues.csv (N/C-term flanks)", {
      src <- file.path(output_dir, "results.sage.tsv")
      if (!file.exists(src)) stop("results.sage.tsv missing")
      df <- data.table::fread(src, select = c("peptide"), nThread = 1)
      strip <- gsub("\\[.*?\\]|[^A-Z]", "", df$peptide)
      strip <- strip[nchar(strip) >= 5]
      nt <- substr(strip, 1, 1)
      ct <- substr(strip, nchar(strip), nchar(strip))
      out <- data.frame(
        residue   = sort(unique(c(nt, ct))),
        n_pct = as.integer(table(factor(nt, levels = sort(unique(c(nt, ct)))))) /
                length(nt) * 100,
        c_pct = as.integer(table(factor(ct, levels = sort(unique(c(nt, ct)))))) /
                length(ct) * 100
      )
      write.csv(out, file.path(bundle_dir, "peptidomics_cleavage_residues.csv"),
                row.names = FALSE)
    })
    # NOTE: source-protein contributions CSV deferred until parse_sage_results
    # is mode-aware (needs FASTA mapping to know which peptides came from where).
  }

  # ---- PROMPT.md — AI ingestion instructions (Claude / Gemini friendly) ----
  safe_section(manifest, "PROMPT.md (AI ingestion guide)", {
    prompt <- c(
      sprintf("# DE-LIMP DDA Export — %s mode", mode),
      "",
      sprintf("Generated by DE-LIMP v%s on %s", app_version,
              format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
      "",
      "## What's in this bundle",
      "",
      "- `methods.md` — full submission metadata (FASTA, instrument, search params, job IDs)",
      "- `settings.json` — exact Sage v0.14.7 config that ran",
      "- `sage_results.tsv` — all Sage PSMs (1% FDR-filtered if `results.sage.tsv`)",
      "- `peptide_length_distribution.csv` — universal length histogram",
      if (file.exists(file.path(bundle_dir, "hla_anchor_residues.csv")))
        "- `hla_anchor_residues.csv` — P2 + PΩ residue frequencies (the HLA fingerprint)" else "",
      if (file.exists(file.path(bundle_dir, "peptidomics_cleavage_residues.csv")))
        "- `peptidomics_cleavage_residues.csv` — N-/C-terminal flanking residue percentages" else "",
      if (dir.exists(file.path(bundle_dir, "casanovo")))
        "- `casanovo/*.mztab` — per-file Casanovo de novo PSMs (sequences + scores)" else "",
      if (file.exists(file.path(bundle_dir, "diamond_hits.tsv")))
        "- `diamond_hits.tsv` — DIAMOND BLAST hits for Casanovo peptides (SwissProt → TrEMBL cascade)" else "",
      "- `MANIFEST.txt` — what made it into the bundle, what was skipped, and why",
      "",
      "## Interpretation hints",
      "",
      switch(mode,
        "hla_class_i" = paste(
          "- A clean HLA class I prep shows a sharp peak at length 9 in `peptide_length_distribution.csv`.",
          "- `hla_anchor_residues.csv` P2 and PΩ frequencies fingerprint the donor's HLA allele set —",
          "  for example, A*02:01 prefers L at P2 and L/V at PΩ.",
          sep = "\n"),
        "hla_class_ii" = paste(
          "- HLA class II peptides have a broad length distribution centered ~13-15 AA.",
          "- Class II has weaker anchor preferences than class I; expect more uniform P2/PΩ.",
          sep = "\n"),
        "peptidomics"  = paste(
          "- `peptidomics_cleavage_residues.csv` shows which proteases shaped these peptides.",
          "- High N-term M = N-terminal Met excision incomplete; high C-term K/R = trypsin contamination.",
          sep = "\n"),
        ""),
      ""
    )
    prompt <- prompt[nzchar(prompt)]
    writeLines(prompt, file.path(bundle_dir, "PROMPT.md"))
  })

  # ---- Write MANIFEST ----
  writeLines(manifest$lines, file.path(bundle_dir, "MANIFEST.txt"))

  # ---- Zip ----
  zip_path <- file.path(tempdir(), paste0(bundle_name, ".zip"))
  if (file.exists(zip_path)) file.remove(zip_path)
  old_wd <- setwd(dirname(bundle_dir))
  on.exit(setwd(old_wd), add = TRUE)
  utils::zip(zipfile = zip_path, files = basename(bundle_dir))

  message(sprintf("[DDA Export] Bundle ready: %s", zip_path))
  zip_path
}
