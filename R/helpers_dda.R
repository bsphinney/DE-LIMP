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

#' Hyperlink a protein accession to the right database — NCBI Protein for
#' RefSeq/GenBank (XP_, NP_, EGW08387.1, PNI49027.1 …), UniProt for UniProt.
#'
#' The discriminator: NCBI accessions carry an underscore (RefSeq) or a `.N`
#' version suffix (GenBank); UniProt accessions have neither. `sp|ACC|NAME` is
#' unwrapped to ACC. Label defaults to the FULL accession (never the mangled
#' prefix). For DT/HTML rendered with `escape = FALSE`.
#'
#' @param accession character vector of accessions
#' @param label optional display labels (default = accession)
#' @return character vector of HTML anchors
ncbi_protein_link <- function(accession, label = NULL) {
  acc <- as.character(accession)
  lab <- if (is.null(label)) acc else as.character(label)
  out <- lab
  ok  <- !is.na(acc) & nzchar(acc) & acc != "-"
  is_up_piped <- grepl("^[a-z]{2}\\|[^|]+\\|", acc)          # sp|ACC|NAME / tr|ACC|NAME
  upacc <- sub("^[a-z]{2}\\|([^|]+)\\|.*", "\\1", acc)        # -> ACC
  # NCBI only when NOT a piped UniProt id AND (RefSeq underscore or GenBank .version).
  # (the `_` test must come after the pipe check — UniProt NAMEs like KRT5_FELCA contain `_`)
  is_ncbi <- !is_up_piped & (grepl("_", acc) | grepl("\\.[0-9]+$", acc))
  url <- ifelse(is_ncbi,
                paste0("https://www.ncbi.nlm.nih.gov/protein/", acc),
                paste0("https://www.uniprot.org/uniprotkb/", ifelse(is_up_piped, upacc, acc)))
  out[ok] <- sprintf('<a href="%s" target="_blank" rel="noopener">%s</a>', url[ok], lab[ok])
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

#' Reconstruct the REAL gapped alignment from a DIAMOND/BLAST `btop` string + the
#' ungapped aligned query (`qseq`). btop tokens: an integer = that many matches;
#' a 2-char token = a non-match (query char, subject char), where `-` marks a gap.
#' Matched stretches are filled from the query (match ⇒ same residue). Returns the
#' aligned query + subject (equal length, `-` for gaps). No fabrication — purely
#' decodes what DIAMOND reported.
#'
#' @param btop traceback string (e.g. "20", "6QH13", "3-S4", "Q-")
#' @param qseq ungapped aligned query residues (DIAMOND qseq)
#' @return list(qaln, saln) equal-length aligned strings, or NULL if unparseable
parse_btop <- function(btop, qseq) {
  bt <- as.character(btop)
  if (length(bt) == 0 || is.na(bt) || !nzchar(bt)) return(NULL)
  qres <- strsplit(toupper(as.character(qseq)), "")[[1]]
  qi <- 0L; qa <- character(0); sa <- character(0)
  toks <- regmatches(bt, gregexpr("[0-9]+|[A-Za-z*-]{2}", bt))[[1]]
  for (t in toks) {
    if (grepl("^[0-9]+$", t)) {
      n <- as.integer(t)
      if (n > 0) { seg <- qres[qi + seq_len(n)]; qa <- c(qa, seg); sa <- c(sa, seg); qi <- qi + n }
    } else {
      qch <- substr(t, 1, 1); sch <- substr(t, 2, 2)
      if (qch == "-") {            # gap in query (insertion in subject)
        qa <- c(qa, "-"); sa <- c(sa, sch)
      } else {                      # mismatch, or gap in subject (sch == "-")
        qi <- qi + 1L; qa <- c(qa, if (qi <= length(qres)) qres[qi] else qch); sa <- c(sa, sch)
      }
    }
  }
  list(qaln = paste(qa, collapse = ""), saln = paste(sa, collapse = ""))
}

#' Confidence-weighted identity (CWI) for a de novo peptide vs a BLAST subject,
#' computed ONLY from the REAL DIAMOND alignment (qseq/sseq) — never from guessed
#' positions (CLAUDE.md / feedback: never fabricate viz data).
#'
#' Each aligned query residue is weighted by its Casanovo per-residue confidence.
#' A match at a confident residue counts fully; a mismatch where Casanovo was
#' unsure (a likely de novo sequencing error) barely penalizes. So a peptide
#' whose differences fall only at low-confidence residues scores high even at low
#' raw %identity — it formalizes "every residue the model was sure of agrees."
#'
#' @param qseq aligned query sequence (DIAMOND qseq; '-' = gap)
#' @param sseq aligned subject sequence (DIAMOND sseq; '-' = gap)
#' @param aa_scores numeric per-residue Casanovo confidence for the FULL peptide (0..1)
#' @param qstart 1-based start of the alignment within the full peptide (DIAMOND qstart)
#' @param hi_conf high-confidence threshold (default 0.95)
#' @return list(cwi, hca, n_hi, n_pos): CWI %, high-confidence-agreement %,
#'         # high-confidence aligned residues, # aligned query residues
compute_cwi <- function(qseq, sseq, aa_scores = NULL, qstart = 1, hi_conf = 0.95) {
  qc <- strsplit(toupper(as.character(qseq)), "")[[1]]
  sc <- strsplit(toupper(as.character(sseq)), "")[[1]]
  n <- min(length(qc), length(sc))
  if (n == 0) return(list(cwi = NA_real_, hca = NA_real_, n_hi = 0L, n_pos = 0L))
  qpos <- as.integer(qstart) - 1L
  num <- 0; den <- 0; hi_match <- 0L; hi_tot <- 0L; npos <- 0L
  for (i in seq_len(n)) {
    qa <- qc[i]; sa <- sc[i]
    if (qa == "-") next                          # insertion in subject; no query residue
    qpos <- qpos + 1L
    # 'X' = DIAMOND-masked (low-complexity) residue on either side — non-informative,
    # NOT a real substitution. Advance the query position (it maps to a real peptide
    # residue) but exclude it from scoring so masking neither rewards nor penalises.
    if (qa == "X" || sa == "X") next
    conf <- if (!is.null(aa_scores) && qpos <= length(aa_scores) && !is.na(aa_scores[qpos]))
              aa_scores[qpos] else NA_real_
    w <- if (is.na(conf)) 0.5 else max(0, conf)  # weight; clamp any negative to 0
    is_match <- sa != "-" && (qa == sa || (qa %in% c("I", "L") && sa %in% c("I", "L")))
    den <- den + w; npos <- npos + 1L
    if (is_match) num <- num + w
    if (!is.na(conf) && conf >= hi_conf) {
      hi_tot <- hi_tot + 1L; if (is_match) hi_match <- hi_match + 1L
    }
  }
  list(
    cwi  = if (den > 0)    round(100 * num / den, 1)         else NA_real_,
    hca  = if (hi_tot > 0) round(100 * hi_match / hi_tot, 1) else NA_real_,
    n_hi = hi_tot, n_pos = npos)
}

#' One-call CWI for a single BLAST hit row — the DRY entry point used by BOTH the
#' alignment viewer and the near-match table, so the metric has ONE definition
#' (CLAUDE.md architectural rule 3). Reconstructs the REAL gapped alignment from
#' btop (DIAMOND qseq/sseq are ungapped and unequal-length when the hit has gaps,
#' so they CANNOT be compared column-wise directly), then weights by per-residue
#' Casanovo confidence. If btop is absent the alignment cannot be reconstructed,
#' so CWI is reported as unavailable rather than fabricated.
#'
#' @param qseq,sseq  DIAMOND qseq/sseq (ungapped aligned subsequences)
#' @param btop       DIAMOND traceback string (needed to insert gaps)
#' @param qstart     1-based alignment start within the full peptide
#' @param aa_scores_str comma-separated per-residue Casanovo confidence (full peptide)
#' @param hi_conf    high-confidence threshold (default 0.95)
#' @return list(available, cwi, hca, n_hi, n_pos, qaln, saln, aa, qstart)
dda_alignment_cwi <- function(qseq, sseq, btop, qstart = 1,
                              aa_scores_str = NULL, hi_conf = 0.95) {
  na_out <- list(available = FALSE, cwi = NA_real_, hca = NA_real_,
                 n_hi = 0L, n_pos = 0L, qaln = NULL, saln = NULL,
                 aa = NULL, qstart = as.integer(qstart %||% 1L))
  if (is.null(qseq) || is.na(qseq) || !nzchar(as.character(qseq))) return(na_out)
  if (is.null(btop) || is.na(btop) || !nzchar(as.character(btop)))  return(na_out)
  aa <- NULL
  if (!is.null(aa_scores_str) && !is.na(aa_scores_str) &&
      nzchar(as.character(aa_scores_str)) && aa_scores_str != "null") {
    aa <- suppressWarnings(as.numeric(strsplit(as.character(aa_scores_str), ",")[[1]]))
  }
  aln <- parse_btop(btop, qseq)
  if (is.null(aln)) return(na_out)
  cw <- compute_cwi(aln$qaln, aln$saln, aa, qstart = qstart, hi_conf = hi_conf)
  list(available = TRUE, cwi = cw$cwi, hca = cw$hca, n_hi = cw$n_hi,
       n_pos = cw$n_pos, qaln = aln$qaln, saln = aln$saln,
       aa = aa, qstart = as.integer(qstart %||% 1L))
}

#' Three-tier biologist-facing call from the confidence-weighted metrics.
#' "Confident" = every residue Casanovo was sure of agrees with the reference
#' (so any disagreement is at de-novo-error residues); "Uncertain" = the model's
#' confident residues themselves disagree (a real divergence or a wrong hit).
#' Returns NA when there aren't enough high-confidence residues to judge.
#'
#' @param hca  high-confidence-agreement %  (from compute_cwi)
#' @param n_hi number of high-confidence aligned residues
#' @return one of "Confident", "Likely", "Uncertain", or NA_character_
cwi_call_label <- function(hca, n_hi) {
  if (is.na(hca) || is.na(n_hi) || n_hi < 2) return(NA_character_)
  if (hca >= 90 && n_hi >= 3) return("Confident")
  if (hca >= 75)              return("Likely")
  "Uncertain"
}

#' Per-peptide Casanovo aa_scores lookup (best PSM), keyed by the canonical
#' (I/L-normalized) peptide. ONE definition shared by the Master Table and the
#' BLAST alignment table so the CWI join uses the same key everywhere.
#'
#' @param psms Casanovo PSM table (needs seq_norm + aa_scores; score optional)
#' @return named character vector seq_norm -> aa_scores string, or NULL
build_casanovo_aa_lookup <- function(psms) {
  if (is.null(psms) || !is.data.frame(psms) ||
      !all(c("seq_norm", "aa_scores") %in% names(psms))) return(NULL)
  sc <- if ("score" %in% names(psms)) suppressWarnings(as.numeric(psms$score)) else rep(0, nrow(psms))
  o <- order(-ifelse(is.na(sc), -Inf, sc))
  p <- psms[o, , drop = FALSE]
  p <- p[!duplicated(p$seq_norm), , drop = FALSE]
  stats::setNames(as.character(p$aa_scores), p$seq_norm)
}

#' Best BLAST hit per canonical peptide — alignment length, e-value, bitscore.
#'
#' "Best" = highest bitscore. Used to surface query coverage (alnlen / peptide
#' length) + e-value so a 100%-identity-but-partial hit (e.g. 18 of 24 residues
#' to an over-represented taxon) can't masquerade as a full match. Pure.
#'
#' @param blast BLAST hits (cols: peptide/query, length, evalue, bitscore;
#'        optionally qseq/sseq/btop/qstart to carry the alignment for CWI)
#' @return data.table(.k, blast_aln_len, blast_evalue, blast_bitscore
#'         [+ blast_qseq/blast_sseq/blast_btop/blast_qstart when present]) or NULL
best_blast_hit_per_peptide <- function(blast) {
  if (is.null(blast) || !is.data.frame(blast) || nrow(blast) == 0) return(NULL)
  bt <- data.table::as.data.table(blast)
  kcol <- if ("peptide" %in% names(bt)) "peptide" else if ("query" %in% names(bt)) "query" else NULL
  if (is.null(kcol)) return(NULL)
  bt$.k <- gsub("I", "L", build_dda_canonical_peptide(bt[[kcol]]))
  if ("bitscore" %in% names(bt)) bt <- bt[order(-as.numeric(bt$bitscore)), ]
  bt <- bt[!duplicated(bt$.k), ]
  alnc <- if ("length" %in% names(bt)) "length" else if ("aln_len" %in% names(bt)) "aln_len" else NA
  out <- data.table::data.table(
    .k             = bt$.k,
    blast_aln_len  = if (!is.na(alnc)) as.numeric(bt[[alnc]]) else NA_real_,
    blast_evalue   = if ("evalue" %in% names(bt)) as.numeric(bt$evalue) else NA_real_,
    blast_bitscore = if ("bitscore" %in% names(bt)) as.numeric(bt$bitscore) else NA_real_)
  # carry the best hit's real alignment so CWI + protein coverage can be derived
  if ("subject" %in% names(bt)) out$blast_subject <- as.character(bt$subject)
  if ("qseq" %in% names(bt)) out$blast_qseq   <- as.character(bt$qseq)
  if ("sseq" %in% names(bt)) out$blast_sseq   <- as.character(bt$sseq)
  if ("btop" %in% names(bt)) out$blast_btop   <- as.character(bt$btop)
  if ("qstart" %in% names(bt)) out$blast_qstart <- suppressWarnings(as.integer(bt$qstart))
  if ("sstart" %in% names(bt)) out$blast_sstart <- suppressWarnings(as.integer(bt$sstart))
  if ("send" %in% names(bt))   out$blast_send   <- suppressWarnings(as.integer(bt$send))
  out
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
#' @param blast      optional BLAST hits; if they carry qseq/sseq/btop the best
#'        hit's real alignment drives the CWI columns
#' @param aa_lookup  optional named character vector (canonical seq_norm ->
#'        Casanovo aa_scores string) used to weight CWI/HCA by per-residue
#'        confidence. Without it, CWI still computes (unweighted) but HCA is NA.
#' @return data.frame with seq_norm + biologist-facing columns
build_denovo_master <- function(classified, sage_psms = NULL, lca = NULL,
                                blast = NULL, aa_lookup = NULL) {
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

    # Confidence-Weighted Identity from the best hit's REAL alignment (btop).
    # One definition (dda_alignment_cwi); HCA needs aa_scores via aa_lookup.
    # PERF: only a small fraction of peptides have a BLAST hit (btop); iterate
    # ONLY those rows over PRE-EXTRACTED vectors. Looping all N rows with per-row
    # data.table `$[i]` indexing was O(N) slow (42s on a 312k-peptide dataset).
    if (all(c("blast_qseq", "blast_sseq", "blast_btop") %in% names(pep))) {
      n <- nrow(pep)
      cwi <- rep(NA_real_, n); hca <- rep(NA_real_, n); call <- rep(NA_character_, n)
      btop_v  <- pep$blast_btop
      qseq_v  <- pep$blast_qseq
      sseq_v  <- pep$blast_sseq
      qst_v   <- if ("blast_qstart" %in% names(pep)) pep$blast_qstart else rep(1L, n)
      norm_v  <- pep$seq_norm
      hits <- which(!is.na(btop_v) & nzchar(btop_v))      # ~the BLAST-hit peptides only
      # Resolve aa_scores for the hit peptides in ONE vectorised lookup. Indexing a
      # 300k-element NAMED vector inside the loop is a linear name-scan per call
      # (O(hits x N) ~ billions of comparisons); doing it once is O(N).
      aa_hits <- if (!is.null(aa_lookup)) unname(aa_lookup[norm_v[hits]]) else rep(NA_character_, length(hits))
      for (j in seq_along(hits)) {
        i <- hits[j]
        rr <- dda_alignment_cwi(qseq_v[i], sseq_v[i], btop_v[i],
                                qstart = qst_v[i] %||% 1L, aa_scores_str = aa_hits[j])
        if (isTRUE(rr$available)) {
          cwi[i] <- rr$cwi; hca[i] <- rr$hca; call[i] <- cwi_call_label(rr$hca, rr$n_hi)
        }
      }
      pep$CWI <- cwi; pep$HCA <- hca; pep$Call <- call
    }
    # drop the raw alignment/subject carrier columns (not for the peptide display;
    # the protein rollup re-derives them straight from the BLAST table)
    for (cc in c("blast_qseq", "blast_sseq", "blast_btop", "blast_qstart",
                 "blast_subject", "blast_sstart", "blast_send"))
      if (cc %in% names(pep)) pep[, (cc) := NULL]
  }
  as.data.frame(pep)
}

#' FragPipe/Philosopher-style PARSIMONIOUS protein inference for the de novo
#' peptides, over the peptide->protein BLAST hit graph. Pure + testable.
#'
#' Mirrors the deterministic model FragPipe reports (parsimony + razor peptides
#' + indistinguishable-protein grouping):
#'   - Parsimony (Occam's razor): the minimal set of proteins explaining all
#'     peptides, found by greedy set-cover — repeatedly take the protein
#'     covering the most still-unexplained peptides.
#'   - Razor peptide: a shared peptide is credited to the single most-supported
#'     protein (the one that claims it during the greedy pass).
#'   - Unique peptides: map to only one protein. Total: every peptide that can
#'     map to the protein (here or elsewhere). Razor = unique + shared-assigned.
#'   - Indistinguishable proteins: proteins with the identical peptide set are
#'     collapsed under one lead and listed.
#' NOTE: this is the parsimony/razor MODEL, not ProteinProphet's probabilistic
#' EM — we have de novo + homology evidence, not DB-search PSM probabilities, so
#' there is no protein probability. Edges are BLAST hits with pident >= min_pident
#' (homology, not exact tryptic membership) — stated so it is not mistaken for a
#' database-search inference.
#'
#' @param blast_all the FULL BLAST hit graph (NOT best-hit-deduped): peptide/query
#'        + subject + pident (+ bitscore, sstart, send for coverage/ordering)
#' @param master    optional build_denovo_master() frame (seq_norm + per-peptide
#'        Casanovo_score / Species_or_clade / Diagnostic / CWI / Call / n_PSMs)
#' @param min_pident edge threshold for "peptide maps to protein" (default 50)
#' @return data.frame: Protein, Accession, Species, Unique_peptides,
#'         Razor_peptides, Total_peptides, N_indistinguishable, N_diagnostic,
#'         Best_CWI, Best_Casanovo, AA_covered, Found_by, Evidence, Subject
build_denovo_protein_groups <- function(blast_all, master = NULL, min_pident = 50) {
  if (is.null(blast_all) || !is.data.frame(blast_all) || nrow(blast_all) == 0)
    return(data.frame())
  bt <- data.table::as.data.table(blast_all)
  pcol <- if ("peptide" %in% names(bt)) "peptide" else if ("query" %in% names(bt)) "query" else NULL
  if (is.null(pcol) || !all(c("subject", "pident") %in% names(bt))) return(data.frame())
  bt <- bt[suppressWarnings(as.numeric(pident)) >= min_pident &
           !is.na(subject) & nzchar(as.character(subject))]
  if (nrow(bt) == 0) return(data.frame())
  bt$pep <- gsub("I", "L", build_dda_canonical_peptide(bt[[pcol]]))
  edges <- unique(bt[, list(pep, subject = as.character(subject))])

  peps  <- unique(edges$pep)
  prots <- unique(edges$subject)
  pep_i  <- stats::setNames(seq_along(peps),  peps)
  prot_i <- stats::setNames(seq_along(prots), prots)
  ep <- pep_i[edges$pep]; eq <- prot_i[edges$subject]
  prot2pep <- split(ep, eq)                       # protein idx -> peptide idxs
  pep2prot <- split(eq, ep)                        # peptide idx -> protein idxs
  total_per_prot <- lengths(prot2pep)              # Total peptides (all mappable)
  pep_degree     <- lengths(pep2prot)              # 1 == unique peptide

  # ---- greedy parsimony / razor assignment ----
  remaining <- total_per_prot                      # unexplained peptides per protein
  explained <- logical(length(peps))
  assigned  <- vector("list", length(prots))       # protein idx -> razor peptide idxs
  repeat {
    p <- which.max(remaining)
    if (length(p) == 0 || remaining[p] <= 0) break
    cand <- prot2pep[[p]]
    take <- cand[!explained[cand]]
    if (!length(take)) { remaining[p] <- 0L; next }
    assigned[[p]] <- take
    explained[take] <- TRUE
    for (pp in take) for (q in pep2prot[[pp]]) remaining[q] <- remaining[q] - 1L
  }
  lead <- which(lengths(assigned) > 0)
  if (!length(lead)) return(data.frame())

  # indistinguishable proteins: non-lead proteins whose full peptide set equals a
  # lead's (same evidence). PERF: hash each set to a signature string and count
  # via a table lookup — O(proteins), not the old O(proteins x leads) identical().
  sig_of <- function(p) paste0(sort(prot2pep[[p]]), collapse = ",")
  lead_sig <- vapply(lead, sig_of, character(1))
  nonlead <- setdiff(seq_along(prots), lead)
  n_indist <- integer(length(lead))
  if (length(nonlead)) {
    nl_tab <- table(vapply(nonlead, sig_of, character(1)))
    hit <- nl_tab[lead_sig]
    n_indist <- as.integer(ifelse(is.na(hit), 0L, hit))
  }

  # per-peptide attributes from the master (keyed by canonical seq_norm). PERF:
  # resolve peptide -> master-row ONCE (vectorised match); attr_for then does a
  # direct integer index instead of a match() over the full master per call.
  ma <- NULL; ma_row <- NULL
  if (!is.null(master) && is.data.frame(master) && "seq_norm" %in% names(master)) {
    ma <- data.table::as.data.table(master)
    ma_row <- match(peps, ma$seq_norm)
  }
  attr_for <- function(pep_idx, col, fun, default = NA) {
    if (is.null(ma) || !col %in% names(ma)) return(default)
    fun(ma[[col]][ma_row[pep_idx]])
  }
  mode_chr <- function(x) { x <- x[!is.na(x) & nzchar(as.character(x))]
    if (!length(x)) return(NA_character_); names(sort(table(x), decreasing = TRUE))[1] }
  # PERF: key the BLAST rows by subject once so per-protein coverage is a binary-
  # search lookup, not a full O(edges) scan per lead.
  have_span <- all(c("sstart", "send") %in% names(bt))
  bt_key <- NULL
  if (have_span) {
    bt_key <- data.table::data.table(
      pep = bt$pep, subject = as.character(bt$subject),
      sstart = suppressWarnings(as.integer(bt$sstart)),
      send   = suppressWarnings(as.integer(bt$send)))
    data.table::setkey(bt_key, subject)
  }
  span_union <- function(pep_idx, subj) {          # subject residues covered (union)
    if (!have_span) return(NA_integer_)
    e <- bt_key[list(subj), nomatch = 0L]
    if (!nrow(e)) return(NA_integer_)
    e <- e[e$pep %in% peps[pep_idx]]
    s <- e$sstart; en <- e$send
    ok <- !is.na(s) & !is.na(en); if (!any(ok)) return(NA_integer_)
    s <- s[ok]; en <- en[ok]; o <- order(s); s <- s[o]; en <- en[o]
    tot <- 0L; cs <- s[1]; ce <- en[1]
    for (i in seq_along(s)[-1]) {
      if (s[i] <= ce + 1L) ce <- max(ce, en[i])
      else { tot <- tot + (ce - cs + 1L); cs <- s[i]; ce <- en[i] } }
    as.integer(tot + (ce - cs + 1L))
  }

  rows <- lapply(seq_along(lead), function(li) {
    p <- lead[li]; subj <- prots[p]
    raz <- assigned[[p]]                            # razor peptides (unique + shared-assigned)
    allp <- prot2pep[[p]]
    data.frame(
      Subject          = subj,
      Protein          = dda_protein_label(subj),
      Accession        = sub("^[a-z]+\\|([^|]+)\\|.*", "\\1", subj),
      Species          = attr_for(raz, "Species_or_clade", mode_chr, NA_character_),
      Unique_peptides  = sum(pep_degree[allp] == 1L),
      Razor_peptides   = length(raz),
      Total_peptides   = length(allp),
      N_indistinguishable = n_indist[li],
      N_diagnostic     = attr_for(raz, "Diagnostic", function(v) sum(v %in% c(TRUE,1,"1"), na.rm=TRUE), NA_integer_),
      Best_CWI         = attr_for(raz, "CWI", function(v) { v<-suppressWarnings(max(as.numeric(v),na.rm=TRUE)); if(is.finite(v)) round(v,1) else NA_real_ }, NA_real_),
      Best_Casanovo    = attr_for(raz, "Casanovo_score", function(v) { v<-suppressWarnings(max(as.numeric(v),na.rm=TRUE)); if(is.finite(v)) round(v,3) else NA_real_ }, NA_real_),
      N_confident      = attr_for(raz, "Call", function(v) sum(v=="Confident", na.rm=TRUE), NA_integer_),
      AA_covered       = span_union(raz, subj),
      Found_by         = attr_for(raz, "Found_by_Sage", function(v) if (any(v %in% c(TRUE,1,"1"))) "Sage + de novo" else "de novo only", "de novo only"),
      stringsAsFactors = FALSE)
  })
  g <- do.call(rbind, rows)
  nconf <- ifelse(is.na(g$N_confident), 0L, g$N_confident)
  g$Evidence <- ifelse(g$Razor_peptides >= 3 | nconf >= 2, "Strong",
                ifelse(g$Razor_peptides >= 2 | nconf >= 1, "Moderate", "Weak"))
  g <- g[order(-g$Razor_peptides, -ifelse(is.na(g$Best_CWI), -1, g$Best_CWI)), , drop = FALSE]
  front <- c("Protein", "Accession", "Species", "Unique_peptides", "Razor_peptides",
             "Total_peptides", "N_indistinguishable", "N_diagnostic", "N_confident",
             "Best_CWI", "Best_Casanovo", "AA_covered", "Found_by", "Evidence", "Subject")
  g[, c(intersect(front, names(g)), setdiff(names(g), front)), drop = FALSE]
}

#' Per-residue protein COVERAGE track for one subject protein — lays every de
#' novo peptide that hits this protein onto the subject coordinate axis and
#' classifies each covered position. Pure + testable; the server turns it into
#' HTML. Uses ONLY the real btop alignment to this protein (no fabrication).
#'
#' state per subject position: "uncovered" | "match" | "variant" (a confident
#' de novo residue disagrees -> genuine difference) | "error" (low-confidence
#' disagreement -> likely de novo error) | "amber" (mid / gap / no score).
#'
#' @param subject   subject accession to map
#' @param blast     FULL BLAST hits (peptide/query + subject/qseq/sseq/btop/qstart/sstart/send)
#' @param aa_lookup optional seq_norm -> aa_scores string (for confidence)
#' @return data.frame(pos, ref, state, n_cov, n_match, n_mismatch) or NULL
build_protein_coverage_track <- function(subject, blast, aa_lookup = NULL,
                                         hi_conf = 0.95) {
  if (is.null(blast) || !is.data.frame(blast) || nrow(blast) == 0) return(NULL)
  need <- c("subject", "sstart", "send", "qseq", "sseq", "btop")
  if (!all(need %in% names(blast))) return(NULL)
  subj_id <- as.character(subject)
  bt <- data.table::as.data.table(blast)
  pcol <- if ("peptide" %in% names(bt)) "peptide" else if ("query" %in% names(bt)) "query" else NULL
  if (is.null(pcol)) return(NULL)
  bt$.kall <- gsub("I", "L", build_dda_canonical_peptide(bt[[pcol]]))
  # Peptide "uniqueness" = maps to only ONE protein in the full graph (FragPipe
  # unique-peptide sense). Computed across ALL subjects (>=50% identity edges)
  # before subsetting, so we can flag which residues are pinned to THIS protein
  # specifically vs covered only by peptides shared with homologs.
  deg <- if ("pident" %in% names(bt))
    bt[suppressWarnings(as.numeric(pident)) >= 50,
       list(.deg = data.table::uniqueN(subject)), by = ".kall"] else
    bt[, list(.deg = data.table::uniqueN(subject)), by = ".kall"]
  uniq_keys <- deg$.kall[deg$.deg == 1L]
  rows <- bt[as.character(bt$subject) == subj_id, ]   # base index — avoids data.table NSE name clash
  if (nrow(rows) == 0) return(NULL)
  rows$.k <- rows$.kall
  rows$.uniq <- rows$.k %in% uniq_keys
  if ("bitscore" %in% names(rows)) rows <- rows[order(-suppressWarnings(as.numeric(rows$bitscore))), ]
  rows <- rows[!duplicated(rows$.k), ]              # best hit per peptide TO this protein
  rows$sstart <- suppressWarnings(as.integer(rows$sstart))
  rows$send   <- suppressWarnings(as.integer(rows$send))
  rows <- rows[!is.na(rows$sstart) & !is.na(rows$send), ]
  if (nrow(rows) == 0) return(NULL)
  qst <- if ("qstart" %in% names(rows)) suppressWarnings(as.integer(rows$qstart)) else rep(1L, nrow(rows))
  # Show the WHOLE protein (residue 1 -> protein length), not just the covered
  # span, so the user can see where the matched peptides sit. The full length
  # comes from `slen` (subject length) when the BLAST output carries it; without
  # it we can only extend to the last covered residue (the C-terminal tail beyond
  # coverage is genuinely unknown — never fabricated).
  smin <- 1L
  slen_vec <- if ("slen" %in% names(rows)) suppressWarnings(as.integer(rows$slen)) else NA_integer_
  len_exact <- any(!is.na(slen_vec)) && max(slen_vec, na.rm = TRUE) >= max(rows$send)
  smax <- if (len_exact) max(slen_vec, na.rm = TRUE) else max(rows$send)
  L <- smax - smin + 1L
  if (!is.finite(L) || L <= 0 || L > 100000) return(NULL)
  ref <- rep(NA_character_, L); ncov <- integer(L); nuniq <- integer(L)
  obs_votes <- vector("list", L)                    # de novo residue calls per position
  conf_votes <- vector("list", L)                   # Casanovo confidence aligned to each call
  for (i in seq_len(nrow(rows))) {
    aln <- parse_btop(rows$btop[i], rows$qseq[i])
    if (is.null(aln)) next
    qa <- strsplit(aln$qaln, "")[[1]]; sa <- strsplit(aln$saln, "")[[1]]
    aa <- NULL
    if (!is.null(aa_lookup)) {
      s <- unname(aa_lookup[rows$.k[i]])
      if (!is.na(s) && nzchar(s)) aa <- suppressWarnings(as.numeric(strsplit(s, ",")[[1]]))
    }
    spos <- rows$sstart[i] - 1L
    qpos <- (qst[i] %||% 1L) - 1L
    for (cc in seq_along(sa)) {
      qch <- qa[cc]; sch <- sa[cc]
      if (qch != "-") qpos <- qpos + 1L
      if (sch == "-") next                          # gap in subject: no subject position
      spos <- spos + 1L
      idx <- spos - smin + 1L
      if (idx < 1L || idx > L) next
      ref[idx] <- sch; ncov[idx] <- ncov[idx] + 1L
      if (isTRUE(rows$.uniq[i])) nuniq[idx] <- nuniq[idx] + 1L
      conf <- if (qch != "-" && !is.null(aa) && qpos <= length(aa) && !is.na(aa[qpos])) aa[qpos] else NA_real_
      obs_votes[[idx]]  <- c(obs_votes[[idx]], qch)  # "-" = de novo lacks this residue
      conf_votes[[idx]] <- c(conf_votes[[idx]], conf)
    }
  }
  # CONSENSUS per position: the de novo residue most peptides call. A position is
  # a "difference" ONLY when that consensus differs from the reference — so the
  # coloured marker never sits over a residue that actually matches (I/L treated
  # as equivalent). Colour comes from the confidence of the consensus call.
  obs <- rep(NA_character_, L); state <- rep("uncovered", L)
  nmatch <- integer(L); nmis <- integer(L)
  for (idx in seq_len(L)) {
    if (ncov[idx] == 0L) next
    v <- obs_votes[[idx]]; cf <- conf_votes[[idx]]
    o <- names(sort(table(v), decreasing = TRUE))[1]; obs[idx] <- o
    rL <- gsub("I", "L", ref[idx]); vL <- gsub("I", "L", v)
    # 'X' = DIAMOND-masked low-complexity; non-informative (exclude from match/
    # mismatch counts; a masked consensus is shown as 'masked', not a difference).
    informative <- v != "-" & v != "X" & rL != "X"
    mv <- vL == rL & informative
    nmatch[idx] <- sum(mv); nmis[idx] <- sum(informative & !mv)
    if (rL == "X" || (!is.na(o) && o == "X")) { state[idx] <- "masked"; next }
    if (!is.na(o) && o != "-" && gsub("I", "L", o) == rL) { state[idx] <- "match"; next }
    if (is.na(o) || o == "-") { state[idx] <- "amber"; next }   # consensus deletion
    oc <- cf[v == o]; oc <- oc[!is.na(oc)]
    oconf <- if (length(oc)) max(oc) else NA_real_
    state[idx] <- if (!is.na(oconf) && oconf >= hi_conf) "variant" else
                  if (!is.na(oconf) && oconf < 0.70) "error" else "amber"
  }
  out <- data.frame(pos = smin:smax, ref = ref, obs = obs, state = state, n_cov = ncov,
             n_match = nmatch, n_mismatch = nmis, n_unique = nuniq,
             unique_support = nuniq > 0L, stringsAsFactors = FALSE)
  attr(out, "len_exact") <- isTRUE(len_exact)   # TRUE only when slen gave the true protein length
  # Per-peptide spans for the tiling view (PEAKS-style bars under the sequence):
  # one row per de novo peptide best-hitting this protein, with its subject span,
  # mean Casanovo confidence, and whether it is unique to this protein.
  pep_mc <- rep(NA_real_, nrow(rows))
  if (!is.null(aa_lookup)) {
    for (i in seq_len(nrow(rows))) {
      s <- unname(aa_lookup[rows$.k[i]])
      if (!is.na(s) && nzchar(s)) {
        v <- suppressWarnings(as.numeric(strsplit(s, ",")[[1]]))
        if (length(v)) pep_mc[i] <- mean(v, na.rm = TRUE)
      }
    }
  }
  attr(out, "peptides") <- data.frame(
    pep = rows$.k, sstart = rows$sstart, send = rows$send,
    uniq = rows$.uniq, mean_conf = pep_mc, stringsAsFactors = FALSE)
  out
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

#' Correct NovoBoard target-decoy SPECTRA competition + FDR, computed from the loaded
#' dataset. Unlike build_denovo_score_calibration() (which only works for shuffled-peptide
#' decoys), this bins BOTH the real and the decoy-spectra peptide populations by their OWN
#' Casanovo score, so it is valid for decoy spectra. See docs/DENOVO_FDR_VALIDATION.md.
#'
#' @param real  data.frame with numeric `score` (Casanovo peptide score) + logical `hit` (>=1 nr BLAST hit)
#' @param decoy data.frame with the same two columns, for the decoy-spectra population
#' @param bin   score bin width for the hit-rate-vs-score panel (default 0.1)
#' @param score_range plotted score range (default -1..1)
#' @return list(bins, fdr, enrichment, real_n, decoy_n, real_hit, decoy_hit,
#'   real_rate, decoy_rate, n_at_1pct) — or NULL if inputs are unusable
build_denovo_spectra_fdr <- function(real, decoy, bin = 0.1, score_range = c(-1, 1)) {
  ok <- function(d) is.data.frame(d) && all(c("score", "hit") %in% names(d)) && nrow(d) > 0
  if (!ok(real) || !ok(decoy)) return(NULL)
  real  <- real[is.finite(real$score), , drop = FALSE]
  decoy <- decoy[is.finite(decoy$score), , drop = FALSE]
  if (nrow(real) == 0 || nrow(decoy) == 0) return(NULL)
  rr <- mean(real$hit); dr <- mean(decoy$hit)

  brks <- seq(score_range[1], score_range[2], by = bin)
  mids <- (utils::head(brks, -1) + utils::tail(brks, -1)) / 2
  binhr <- function(d) {
    b <- cut(d$score, brks, include.lowest = TRUE)
    as.numeric(tapply(d$hit, b, function(x) 100 * mean(x)))[seq_along(mids)]
  }
  bins <- data.frame(mid = mids, real_hr = binhr(real), decoy_hr = binhr(decoy))

  # cumulative competition FDR on a fine score grid: for peptides with score >= t,
  # FDR = (# decoy hits >= t) / (# real hits >= t).
  ths <- seq(score_range[1], score_range[2], by = 0.02)
  n_real  <- vapply(ths, function(t) sum(real$hit  & real$score  >= t), numeric(1))
  n_decoy <- vapply(ths, function(t) sum(decoy$hit & decoy$score >= t), numeric(1))
  fdr <- data.frame(threshold = ths, n_real = n_real, n_decoy = n_decoy,
                    fdr = ifelse(n_real > 0, n_decoy / n_real, NA_real_))
  # q-value: min FDR over all equally-or-more-inclusive sets (lower threshold). ths is
  # ascending, so that is a prefix cummin — q is non-decreasing as the accepted set grows.
  fdr$q <- cummin(ifelse(is.na(fdr$fdr), 1, fdr$fdr))
  ok1 <- fdr[is.finite(fdr$q) & fdr$q <= 0.01, , drop = FALSE]
  list(bins = bins, fdr = fdr,
       enrichment = if (dr > 0) rr / dr else NA_real_,
       real_n = nrow(real), decoy_n = nrow(decoy),
       real_hit = sum(real$hit), decoy_hit = sum(decoy$hit),
       real_rate = 100 * rr, decoy_rate = 100 * dr,
       n_at_1pct = if (nrow(ok1)) max(ok1$n_real) else 0L)
}

#' De-novo-ONLY target-decoy FDR (the original NovoBoard FDR, BEFORE any BLAST): compete
#' real-spectra Casanovo PSMs against decoy-spectra PSMs by Casanovo score alone.
#' FDR(t) = (# decoy PSMs with score >= t) / (# real PSMs with score >= t). This is the pure
#' de-novo *sequencing* FDR — it answers "how many of my accepted de novo PSMs are likely
#' wrong", independent of homology. See docs/DENOVO_FDR_VALIDATION.md.
#'
#' @param real_scores  numeric Casanovo peptide scores for the real spectra
#' @param decoy_scores numeric Casanovo peptide scores for the decoy spectra
#' @param bin score bin width for the distribution panel; @param score_range plotted range
#' @return list(dist, fdr, n_at_1pct, real_n, decoy_n) or NULL
build_denovo_casanovo_fdr <- function(real_scores, decoy_scores,
                                      bin = 0.05, score_range = c(-1, 1)) {
  real_scores  <- real_scores[is.finite(real_scores)]
  decoy_scores <- decoy_scores[is.finite(decoy_scores)]
  if (length(real_scores) == 0 || length(decoy_scores) == 0) return(NULL)
  brks <- seq(score_range[1], score_range[2], by = bin)
  mids <- (utils::head(brks, -1) + utils::tail(brks, -1)) / 2
  frac <- function(s) as.numeric(table(cut(s, brks, include.lowest = TRUE)))[seq_along(mids)] / length(s) * 100
  dist <- data.frame(mid = mids, real_pct = frac(real_scores), decoy_pct = frac(decoy_scores))
  ths <- seq(score_range[1], score_range[2], by = 0.02)
  n_real  <- vapply(ths, function(t) sum(real_scores  >= t), numeric(1))
  n_decoy <- vapply(ths, function(t) sum(decoy_scores >= t), numeric(1))
  fdr <- data.frame(threshold = ths, n_real = n_real, n_decoy = n_decoy,
                    fdr = ifelse(n_real > 0, n_decoy / n_real, NA_real_))
  fdr$q <- cummin(ifelse(is.na(fdr$fdr), 1, fdr$fdr))
  ok1 <- fdr[is.finite(fdr$q) & fdr$q <= 0.01, , drop = FALSE]
  list(dist = dist, fdr = fdr, n_at_1pct = if (nrow(ok1)) max(ok1$n_real) else 0L,
       real_n = length(real_scores), decoy_n = length(decoy_scores))
}

#' Best (max) BLAST bitscore per de novo peptide from a DIAMOND hit table. Single
#' definition used by both FDR calibrations below. Returns data.frame(pep, b) or NULL.
denovo_best_bitscore_per_peptide <- function(b) {
  if (is.null(b) || !is.data.frame(b) || nrow(b) == 0) return(NULL)
  bt <- data.table::as.data.table(b)
  qcol <- intersect(c("qseqid", "query", "peptide", "Peptide", "seq_stripped"), names(bt))
  bcol <- intersect(c("bitscore", "bit_score"), names(bt))
  if (length(qcol) == 0 || length(bcol) == 0) return(NULL)
  bt$.q <- as.character(bt[[qcol[1]]])
  bt$.b <- suppressWarnings(as.numeric(bt[[bcol[1]]]))
  bt <- bt[is.finite(bt$.b) & bt$.b > 0 & !is.na(bt$.q) & nzchar(bt$.q), ]
  if (nrow(bt) == 0) return(NULL)
  agg <- bt[, list(b = max(.b)), by = ".q"]
  data.table::setnames(agg, ".q", "pep")
  as.data.frame(agg)
}

#' Confirmed-peptides-by-length table from a vector of accepted peptide ids. Only
#' computed when the ids look like sequences (so length = nchar is meaningful).
denovo_fdr_by_length <- function(peps) {
  if (length(peps) == 0) return(NULL)
  seqlike <- grepl("^[A-Za-z]{4,}$", peps)
  if (mean(seqlike) <= 0.8) return(NULL)
  brks <- c(0, 15, 20, 25, 30, Inf)
  labs <- c("≤15", "16-20", "21-25", "26-30", "31+")
  tb <- as.data.frame(table(length = cut(nchar(peps), brks, labs)))
  names(tb) <- c("length", "n")
  tb
}

#' RECOMMENDED de novo -> homology FDR: bitscore-calibrated decoy-SPECTRA null.
#' The ocelot validation (docs/denovo_decoy_method.html) showed the homology FDR
#' must be ranked by BLAST *bitscore* (or, for standard searches, E-value) and NOT
#' by Casanovo score — the Casanovo peptide score barely separates real from decoy
#' spectra, so it makes a poor threshold (only ~171 of 312,820 peptides clear 1%
#' FDR when gated on it, vs ~5,630 when gated on bitscore). Procedure: best bitscore
#' per peptide; FDR(t) = (decoy hits >= t, scaled by real_n/decoy_n) / (real hits >= t);
#' q-value monotonised over the larger (lower-threshold) sets.
#'
#' Only the BLAST hit tables are needed for the threshold counts (non-hit peptides
#' never enter n(>= t) for t > 0). real_n/decoy_n are used ONLY for the population
#' scale; when decoy_n is unknown we fall back to scale = 1, which over-estimates
#' the FDR when the decoy population is larger (i.e. it is conservative).
#'
#' @param real_blast  DIAMOND hits for REAL de novo peptides (query/peptide + bitscore)
#' @param decoy_blast DIAMOND hits for DECOY-SPECTRA peptides (same columns)
#' @param real_n  unique real de novo peptides searched (population denominator)
#' @param decoy_n unique decoy-spectra peptides searched; NULL -> scale = 1 (conservative)
#' @param fdr_level target FDR (default 0.01)
#' @return list(curve, n_at_fdr, bits_at_fdr, by_length, scale, real_hit, decoy_hit,
#'   global_fdr, fdr_level, mode) or NULL when inputs are unusable
build_denovo_bitscore_fdr <- function(real_blast, decoy_blast,
                                      real_n = NULL, decoy_n = NULL,
                                      fdr_level = 0.01) {
  rb <- denovo_best_bitscore_per_peptide(real_blast)
  db <- denovo_best_bitscore_per_peptide(decoy_blast)
  if (is.null(rb) || is.null(db)) return(NULL)
  scale <- if (!is.null(real_n) && !is.null(decoy_n) &&
               is.finite(real_n) && is.finite(decoy_n) && decoy_n > 0)
             real_n / decoy_n else 1
  ths <- sort(unique(rb$b))                          # ascending bitscore
  n_real  <- vapply(ths, function(t) sum(rb$b >= t), numeric(1))
  n_decoy <- vapply(ths, function(t) sum(db$b >= t), numeric(1))
  fdr <- (n_decoy * scale) / n_real
  q   <- cummin(ifelse(is.finite(fdr), fdr, 1))      # q over larger (lower-t) accepted sets
  curve <- data.frame(bitscore = ths, n_real = n_real, n_decoy = n_decoy, fdr = fdr, q = q)
  ok <- curve[is.finite(curve$q) & curve$q <= fdr_level, , drop = FALSE]
  n_at    <- if (nrow(ok)) max(ok$n_real)   else 0L
  bits_at <- if (nrow(ok)) min(ok$bitscore) else NA_real_
  by_length <- if (n_at > 0 && is.finite(bits_at)) denovo_fdr_by_length(rb$pep[rb$b >= bits_at]) else NULL
  list(curve = curve, n_at_fdr = n_at, bits_at_fdr = bits_at, by_length = by_length,
       scale = scale, real_hit = nrow(rb), decoy_hit = nrow(db),
       global_fdr = (nrow(db) * scale) / nrow(rb), fdr_level = fdr_level, mode = "spectra")
}

#' Decoy-DATABASE competition FDR (reversed-sequence null) — Mike Riffle's control.
#' Per real peptide, compete the best TARGET-DB bitscore against the best DECOY-DB
#' (reversed) bitscore: a peptide is a target-win if best-target >= best-decoy. Rank
#' by the winning bitscore; FDR(t) = (#decoy-wins with score >= t) / (#target-wins
#' with score >= t). This nulls *chance database matches* only — NOT de novo
#' sequencing error — so it is more permissive than the decoy-spectra null (it
#' over-accepts short peptides). Offered alongside the spectra null so the two can
#' be compared (and cross-referenced against the Sage IDs). Same return shape.
#'
#' @param target_blast best target-DB hit per real peptide (query/peptide + bitscore)
#' @param decoy_blast  best reversed-DB hit per real peptide (same columns)
#' @param fdr_level target FDR (default 0.01)
build_denovo_dbcompete_fdr <- function(target_blast, decoy_blast, fdr_level = 0.01) {
  tb <- denovo_best_bitscore_per_peptide(target_blast)
  if (is.null(tb)) return(NULL)
  db <- denovo_best_bitscore_per_peptide(decoy_blast)
  if (is.null(db)) db <- data.frame(pep = character(0), b = numeric(0))
  m <- merge(tb, db, by = "pep", all = TRUE, suffixes = c("_t", "_d"))
  m$b_t[is.na(m$b_t)] <- -Inf; m$b_d[is.na(m$b_d)] <- -Inf
  m <- m[is.finite(pmax(m$b_t, m$b_d)), , drop = FALSE]
  if (nrow(m) == 0) return(NULL)
  m$winT <- m$b_t >= m$b_d
  m$sc   <- pmax(m$b_t, m$b_d)
  ths <- sort(unique(m$sc))
  n_real  <- vapply(ths, function(t) sum(m$winT  & m$sc >= t), numeric(1))
  n_decoy <- vapply(ths, function(t) sum(!m$winT & m$sc >= t), numeric(1))
  fdr <- ifelse(n_real > 0, n_decoy / n_real, NA_real_)
  q   <- cummin(ifelse(is.finite(fdr), fdr, 1))
  curve <- data.frame(bitscore = ths, n_real = n_real, n_decoy = n_decoy, fdr = fdr, q = q)
  ok <- curve[is.finite(curve$q) & curve$q <= fdr_level, , drop = FALSE]
  n_at    <- if (nrow(ok)) max(ok$n_real)   else 0L
  bits_at <- if (nrow(ok)) min(ok$bitscore) else NA_real_
  by_length <- if (n_at > 0 && is.finite(bits_at)) denovo_fdr_by_length(m$pep[m$winT & m$sc >= bits_at]) else NULL
  list(curve = curve, n_at_fdr = n_at, bits_at_fdr = bits_at, by_length = by_length,
       scale = 1, real_hit = sum(m$winT), decoy_hit = sum(!m$winT),
       global_fdr = if (sum(m$winT) > 0) sum(!m$winT) / sum(m$winT) else NA_real_,
       fdr_level = fdr_level, mode = "database")
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

  # FDR filter — and drop DECOYS explicitly. Filtering on q alone lets the
  # occasional decoy (label -1) with q <= 0.01 slip into the analysed set; decoys
  # must never be counted as identifications.
  psms_filtered <- psms[spectrum_q <= fdr_threshold & protein_q <= protein_fdr_threshold]
  if ("label" %in% names(psms_filtered))
    psms_filtered <- psms_filtered[label == 1]
  else if ("is_decoy" %in% names(psms_filtered))
    psms_filtered <- psms_filtered[!as.logical(is_decoy)]
  message("[DDA] After FDR filter (spectrum_q <= ", fdr_threshold,
          ", protein_q <= ", protein_fdr_threshold, ", targets only): ", nrow(psms_filtered), " PSMs")

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


#' Generate the de novo decoy-spectra FDR sbatch (NovoBoard target-decoy null).
#'
#' Codifies the validated method (docs/DENOVO_FDR_VALIDATION.md): for every input
#' mzML, build a decoy spectrum (denovo_decoy_gen.py, FRAC>=0.8 peak replacement,
#' precursor preserved), sequence the decoys with the SAME Casanovo model, extract
#' bare peptides (mods stripped, >=7 aa, deduped), and DIAMOND blastp them against nr
#' with the identical 16-column format used for the real hits. Writes
#' `denovo/blast_results_decoy_spectra.tsv`, which the Target-Decoy FDR pane loads
#' (build_denovo_score_calibration is agnostic to which null produced the file).
#'
#' This is the SAME pipeline as the real arm with the spectra swapped for decoys, so
#' the FDR is a true 1:1 target-decoy competition. Submit after (or alongside) the
#' real Casanovo + nr-BLAST run.
#'
#' @param mzml_dir directory of input .mzML files (the real spectra)
#' @param output_dir search output dir (writes denovo/ + casanovo_decoy/ under it)
#' @param nr_db path to the DIAMOND nr database (no extension)
#' @param model_ckpt Casanovo model checkpoint
#' @param conda_env_path Casanovo v5 conda env
#' @param decoy_gen_script path to denovo_decoy_gen.py (shipped in scripts/)
#' @param frac peak-replacement fraction (default 0.8 — validated floor)
#' @param account,partition,qos SLURM placement for the GPU (Casanovo) step
#' @return character SLURM script (one job: decoy-gen + Casanovo + extract + blast)
generate_denovo_decoy_fdr_sbatch <- function(
    mzml_dir, output_dir, nr_db,
    model_ckpt     = "/quobyte/proteomics-grp/bioinformatics_programs/casanovo_modles/casanovo_v5_0_0.ckpt",
    conda_env_path = "/quobyte/proteomics-grp/conda_envs/casanovo5",
    decoy_gen_script = "/quobyte/proteomics-grp/de-limp/python/denovo_decoy_gen.py",
    frac           = 0.8,
    account        = "genome-center-grp",
    partition      = "gpu-a100",
    qos            = "genome-center-grp-gpu-a100-qos") {
  denovo_dir <- file.path(output_dir, "denovo")
  decoy_dir  <- file.path(output_dir, "casanovo_decoy")
  logs_dir   <- file.path(output_dir, "logs")
  fmt <- "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids qseq sseq btop"
  paste0(
'#!/bin/bash
#SBATCH --job-name=delimp_denovo_decoyfdr
#SBATCH --partition=', partition, '
#SBATCH --account=', account, '
#SBATCH --qos=', qos, '
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output="', logs_dir, '/decoyfdr_%j.out"
#SBATCH --error="', logs_dir, '/decoyfdr_%j.err"

set -euo pipefail
ENV="', conda_env_path, '"
MODEL="', model_ckpt, '"
MZML_DIR="', mzml_dir, '"
DECOY_DIR="', decoy_dir, '"
DENOVO_DIR="', denovo_dir, '"
NR="', nr_db, '"
GEN="', decoy_gen_script, '"
FRAC=', format(frac, nsmall = 1), '
mkdir -p "$DECOY_DIR" "$DENOVO_DIR" "', logs_dir, '"
export PATH="$ENV/bin:$PATH"

echo "[decoy-FDR] decoy-spectra + Casanovo start $(date) on $(hostname)"
for f in "$MZML_DIR"/*.mzML; do
  [ -e "$f" ] || continue
  b=$(basename "$f" .mzML)
  "$ENV/bin/python" "$GEN" "$f" "$DECOY_DIR/decoy_${b}.mgf" "$FRAC"
  casanovo sequence --model "$MODEL" --output_dir "$DECOY_DIR" \\
    --output_root "decoy_${b}_sequence" --force_overwrite "$DECOY_DIR/decoy_${b}.mgf"
done

echo "[decoy-FDR] extracting decoy peptides $(date)"
"$ENV/bin/python" - "$DECOY_DIR" "$DECOY_DIR/decoy_peptides.fasta" <<"PY"
import sys, re, glob, os
ddir, out = sys.argv[1], sys.argv[2]; seqs=set()
for mz in glob.glob(os.path.join(ddir, "decoy_*_sequence.mztab")):
    sc=None
    for line in open(mz):
        if line.startswith("PSH"): sc=line.rstrip("\\n").split("\\t").index("sequence"); continue
        if line.startswith("PSM") and sc is not None:
            s=line.rstrip("\\n").split("\\t")[sc]
            s=re.sub(r"\\[[^\\]]*\\]","",s); s=re.sub(r"\\([^)]*\\)","",s)
            s=re.sub(r"[+-]\\d+\\.?\\d*","",s); s=re.sub(r"[^A-Za-z]","",s).upper()
            if len(s)>=7: seqs.add(s)
open(out,"w").write("".join(f">{p}\\n{p}\\n" for p in sorted(seqs)))
print(f"decoy peptides >=7aa: {len(seqs)}")
PY

echo "[decoy-FDR] DIAMOND blastp vs nr $(date)"
module load diamond/2.1.7 2>/dev/null || true
diamond blastp -d "$NR" -q "$DECOY_DIR/decoy_peptides.fasta" \\
  -o "$DENOVO_DIR/blast_results_decoy_spectra.tsv" \\
  --outfmt ', fmt, ' \\
  --evalue 1 --max-target-seqs 25 --threads 16 --ignore-warnings
echo "[decoy-FDR] DONE $(date); decoy hits=$(wc -l < "$DENOVO_DIR/blast_results_decoy_spectra.tsv")"
')
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

    # Spectrum scan number from spectra_ref ("ms_run[1]:...scan=N"). This is the
    # real scan that aligns with Sage's scannr ("...scan=N") — psm_id is only a
    # row index and must NOT be used to match spectra (the Disagreements bug).
    sr <- if ("spectra_ref" %in% names(df)) df$spectra_ref else rep(NA_character_, nrow(df))
    scan_vec <- suppressWarnings(as.integer(sub(".*scan=([0-9]+).*", "\\1", sr)))

    # Extract key columns (column names from mzTab spec)
    # sequence, PSM_ID, search_engine_score[1], charge, exp_mass_to_charge,
    # calc_mass_to_charge, opt_ms_run[1]_aa_scores
    out <- data.frame(
      sequence     = df$sequence,
      psm_id       = as.integer(df$PSM_ID),
      scan         = scan_vec,
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
      sequence = character(0), psm_id = integer(0), scan = integer(0),
      score = numeric(0), charge = integer(0), exp_mz = numeric(0),
      calc_mz = numeric(0), source_file = character(0), aa_scores = character(0),
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
                                    search_info = NULL,
                                    denovo_cache = NULL) {
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
  safe_section(manifest, "sage_results.tsv (PSM table, 1% FDR, targets only)", {
    src <- file.path(output_dir, "results.sage.tsv")
    if (!file.exists(src)) stop("results.sage.tsv not found")
    # The RAW Sage output is UNFILTERED — it contains decoys (label -1) and the
    # full q-value range. Copying it verbatim ships decoy-inflated counts that a
    # downstream reader would mistake for identifications. Export the SAME table
    # the app analyses: targets only, 1% FDR (spectrum_q <= 0.01 & protein_q <= 0.01).
    ps <- data.table::fread(src)
    n_raw <- nrow(ps)
    keep <- rep(TRUE, nrow(ps))
    n_dec <- 0L
    if ("label" %in% names(ps))      { n_dec <- sum(ps$label == -1, na.rm = TRUE); keep <- keep & ps$label == 1 }
    else if ("is_decoy" %in% names(ps)) { n_dec <- sum(as.logical(ps$is_decoy), na.rm = TRUE); keep <- keep & !as.logical(ps$is_decoy) }
    if ("spectrum_q" %in% names(ps)) keep <- keep & ps$spectrum_q <= 0.01
    if ("protein_q"  %in% names(ps)) keep <- keep & ps$protein_q  <= 0.01
    ps <- ps[keep]
    if (nrow(ps) == 0) stop("no target PSMs pass 1% FDR")
    data.table::fwrite(ps, file.path(bundle_dir, "sage_results.tsv"), sep = "\t")
    message(sprintf("[DDA Export] sage_results.tsv: %d target PSMs at 1%% FDR (from %d raw rows, %d decoys removed)",
                    nrow(ps), n_raw, n_dec))
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
  safe_section(manifest, "diamond_hits.tsv (nr / UniProt BLAST)", {
    cand <- c(file.path(output_dir, "denovo", "blast_results.tsv"),
              file.path(output_dir, "denovo", "diamond_hits.tsv"))
    src <- cand[file.exists(cand)][1]
    if (is.na(src)) stop("no BLAST results (denovo/blast_results.tsv) found")
    file.copy(src, file.path(bundle_dir, "diamond_hits.tsv"), overwrite = TRUE)
  })

  # ---- nr LCA species attribution (per-peptide) ----
  safe_section(manifest, "peptide_lca.tsv (nr LCA species attribution)", {
    cand <- list.files(file.path(output_dir, "denovo"),
                       pattern = "_peptide_lca\\.tsv$", full.names = TRUE)
    if (length(cand) == 0) stop("no *_peptide_lca.tsv found")
    file.copy(tail(sort(cand), 1), file.path(bundle_dir, "peptide_lca.tsv"),
              overwrite = TRUE)   # prefer relaxed-evalue (*_e1) when both exist
  })

  # ---- Protein inference (parsimony + razor) from the BLAST graph + LCA ----
  safe_section(manifest, "protein_groups.tsv (parsimony / razor inference)", {
    bcand <- c(file.path(output_dir, "denovo", "blast_results.tsv"),
               file.path(output_dir, "denovo", "diamond_hits.tsv"))
    bsrc <- bcand[file.exists(bcand)][1]
    if (is.na(bsrc)) stop("no BLAST results for protein inference")
    bl <- data.table::fread(bsrc, header = FALSE)
    cn <- c("peptide","subject","pident","length","mismatch","gapopen",
            "qstart","qend","sstart","send","evalue","bitscore")
    if (ncol(bl) >= length(cn)) data.table::setnames(bl, seq_along(cn), cn)
    if (ncol(bl) >= 16) data.table::setnames(bl, 13:16, c("staxid","qseq","sseq","btop"))
    # Minimal master from the LCA file (species + diagnostic per peptide), so the
    # protein rows carry species attribution even in the headless export path.
    ma <- NULL
    lcand <- list.files(file.path(output_dir, "denovo"),
                        pattern = "_peptide_lca\\.tsv$", full.names = TRUE)
    if (length(lcand)) {
      lt <- tryCatch(data.table::fread(tail(sort(lcand), 1)), error = function(e) NULL)
      if (!is.null(lt) && "peptide" %in% names(lt)) {
        lt$seq_norm <- gsub("I", "L", build_dda_canonical_peptide(lt$peptide))
        keep <- intersect(c("seq_norm","lca_name","diagnostic"), names(lt))
        ma <- as.data.frame(lt[, keep, with = FALSE])
        if ("lca_name" %in% names(ma)) names(ma)[names(ma)=="lca_name"] <- "Species_or_clade"
        if ("diagnostic" %in% names(ma)) names(ma)[names(ma)=="diagnostic"] <- "Diagnostic"
      }
    }
    g <- build_denovo_protein_groups(as.data.frame(bl), master = ma, min_pident = 50)
    if (is.null(g) || nrow(g) == 0) stop("no protein groups inferred")
    data.table::fwrite(g, file.path(bundle_dir, "protein_groups.tsv"), sep = "\t")
  })

  # ---- Processed-data cache (skips the slow mztab re-parse on reload) ----
  # The expensive part of loading is parsing the Casanovo mztabs into ~hundreds
  # of thousands of PSMs + classifying them; BLAST/LCA are fast freads. Caching
  # the parsed PSMs + classification as an RDS lets a reload restore them
  # directly. Same idea as the activity-log session.rds.
  if (!is.null(denovo_cache)) {
    safe_section(manifest, "denovo/denovo_cache.rds (parsed-data cache for fast reload)", {
      dir.create(file.path(bundle_dir, "denovo"), showWarnings = FALSE, recursive = TRUE)
      payload <- list(
        schema         = "delimp-denovo-cache-v1",
        app_version    = app_version,
        casanovo_psms  = denovo_cache$casanovo_psms,
        classification = denovo_cache$classification)
      if (is.null(payload$casanovo_psms)) stop("no parsed PSMs to cache")
      saveRDS(payload, file.path(bundle_dir, "denovo", "denovo_cache.rds"),
              compress = "xz")
    })
  }

  # ---- Annotated-spectrum peaks (for the in-app MS/MS viewer) if pre-extracted ----
  if (all(file.exists(file.path(output_dir, "denovo",
          c("spectra_peaks.parquet", "spectra_meta.parquet"))))) {
    safe_section(manifest, "denovo/spectra_{peaks,meta}.parquet (annotated-spectrum viewer)", {
      dir.create(file.path(bundle_dir, "denovo"), showWarnings = FALSE, recursive = TRUE)
      file.copy(file.path(output_dir, "denovo", "spectra_peaks.parquet"),
                file.path(bundle_dir, "denovo", "spectra_peaks.parquet"), overwrite = TRUE)
      file.copy(file.path(output_dir, "denovo", "spectra_meta.parquet"),
                file.path(bundle_dir, "denovo", "spectra_meta.parquet"), overwrite = TRUE)
    })
  }

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
    # Species summary COMPUTED from this dataset's LCA table — never hardcoded,
    # so it is correct for the ocelot today and any other organism tomorrow.
    lca_lines <- character(0)
    lca_path <- file.path(bundle_dir, "peptide_lca.tsv")
    if (file.exists(lca_path)) {
      lt <- tryCatch(data.table::fread(lca_path), error = function(e) NULL)
      if (!is.null(lt) && nrow(lt) > 0 && "category" %in% names(lt)) {
        cat_tab <- sort(table(lt$category), decreasing = TRUE)
        cat_str <- paste(sprintf("%s (%d)", names(cat_tab), as.integer(cat_tab)),
                         collapse = ", ")
        is_diag <- if ("diagnostic" %in% names(lt)) lt$diagnostic %in% c(1, "1", TRUE) else TRUE
        host <- lt[lt$category == "host" & is_diag, ]
        top_taxa <- if ("lca_name" %in% names(host) && nrow(host) > 0) {
          tt <- sort(table(host$lca_name), decreasing = TRUE)
          k <- seq_len(min(8, length(tt)))
          paste(sprintf("%s (%d)", names(tt)[k], as.integer(tt)[k]), collapse = ", ")
        } else "none"
        lca_lines <- c(
          "",
          "## Species attribution (computed from this dataset's nr LCA)",
          "",
          sprintf("- %s de novo peptides placed by the lowest common ancestor of their nr BLAST hits.",
                  format(nrow(lt), big.mark = ",")),
          sprintf("- Category split: %s.", cat_str),
          sprintf("- Top diagnostic (species/genus) host taxa: %s.", top_taxa),
          "- Conserved peptides resolve only to family or higher and are NOT attributed to one species.",
          "- microbiome / viral hits (incl. nr over-represented taxa such as SARS-CoV-2 spike) are NOT host signal.",
          "- Weigh each hit by Casanovo score (-1..1; >=0 mass-consistent), query coverage, and e-value:",
          "  a 100% identity over partial coverage (e.g. 18 of 24 residues) is NOT a full-length match.",
          "",
          "## Confidence-Weighted Identity (CWI) — how to weigh a hit",
          "",
          "- Raw % identity penalizes a de novo peptide for its own sequencing errors. CWI weights each",
          "  aligned position by Casanovo's per-residue confidence: mismatches at low-confidence residues",
          "  (likely de-novo errors) are discounted; mismatches at high-confidence residues count fully.",
          "- HCA (High-Confidence Agreement) = % of >=0.95-confidence residues that match the reference.",
          "  HCA is the decisive number: a hit where every confident residue agrees is trustworthy even",
          "  at modest raw identity; a hit where confident residues disagree is a real divergence or wrong.",
          "- Call tier: Confident (HCA>=90 & >=3 hi-conf residues) | Likely (HCA>=75) | Uncertain (else).",
          "- This formalizes established de-novo+homology practice (SPIDER/PEAKS; ALPS) as a simple per-hit",
          "  score; it is computed from the real DIAMOND traceback, never from guessed positions.")
      }
    }
    # Protein-level summary, computed from this dataset's protein_groups.tsv.
    prot_lines <- character(0)
    pg_path <- file.path(bundle_dir, "protein_groups.tsv")
    if (file.exists(pg_path)) {
      pg <- tryCatch(data.table::fread(pg_path), error = function(e) NULL)
      if (!is.null(pg) && nrow(pg) > 0) {
        n_strong <- if ("Evidence" %in% names(pg)) sum(pg$Evidence == "Strong", na.rm = TRUE) else NA
        topn <- utils::head(pg, 8)
        lab <- if (all(c("Protein", "Razor_peptides") %in% names(topn)))
          paste(sprintf("%s (%d razor)", topn$Protein, as.integer(topn$Razor_peptides)),
                collapse = ", ") else "n/a"
        prot_lines <- c(
          "",
          "## Protein inference (parsimony + razor; FragPipe/IDPicker model)",
          "",
          sprintf("- %d protein groups inferred%s.", nrow(pg),
                  if (!is.na(n_strong)) sprintf(", %d with Strong evidence", n_strong) else ""),
          sprintf("- Top proteins by razor-peptide support: %s.", lab),
          "- Razor = peptides credited to a protein (unique + shared assigned by parsimony);",
          "  Total = all peptides mappable to it. Many razor peptides = strong support.",
          "- Edges are nr BLAST homology hits (>=50% identity), NOT exact tryptic membership, so",
          "  there is no protein probability. A protein backed by ONE peptide is tentative.",
          "- For species ID, read protein support alongside the per-peptide LCA: agreement between",
          "  many peptides on one protein AND a consistent LCA clade is the strongest call.")
      }
    }
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
      "- `sage_results.tsv` — Sage PSMs, **1% FDR-filtered, targets only** (decoys removed; spectrum_q & protein_q <= 0.01). NOT the raw Sage output — counts here are identifications, not decoy-inflated.",
      "- `peptide_length_distribution.csv` — universal length histogram",
      if (file.exists(file.path(bundle_dir, "hla_anchor_residues.csv")))
        "- `hla_anchor_residues.csv` — P2 + PΩ residue frequencies (the HLA fingerprint)" else "",
      if (file.exists(file.path(bundle_dir, "peptidomics_cleavage_residues.csv")))
        "- `peptidomics_cleavage_residues.csv` — N-/C-terminal flanking residue percentages" else "",
      if (dir.exists(file.path(bundle_dir, "casanovo")))
        "- `casanovo/*.mztab` — per-file Casanovo de novo PSMs (sequences + scores)" else "",
      if (file.exists(file.path(bundle_dir, "diamond_hits.tsv")))
        paste("- `diamond_hits.tsv` — DIAMOND BLAST hits for Casanovo peptides (nr or UniProt).",
              "When present, the `qseq`/`sseq`/`btop` columns hold the real aligned subsequences",
              "and traceback, from which Confidence-Weighted Identity is computed.") else "",
      if (file.exists(file.path(bundle_dir, "peptide_lca.tsv")))
        "- `peptide_lca.tsv` — per-peptide nr lowest-common-ancestor species/clade attribution" else "",
      if (file.exists(file.path(bundle_dir, "protein_groups.tsv")))
        paste("- `protein_groups.tsv` — parsimonious protein inference (FragPipe/IDPicker model):",
              "one row per reference protein with Unique/Razor/Total peptide counts,",
              "N_indistinguishable, N_diagnostic, AA_covered, and a Strong/Moderate/Weak Evidence tier.") else "",
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
      lca_lines,
      prot_lines,
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

# ==============================================================================
#  Annotated MS/MS spectrum support for the in-app de novo spectrum viewer.
#  Pure + testable: compute b/y fragment-ion m/z for a (modified) peptide, then
#  match observed peaks. Lets DE-LIMP render the spectrum behind a de novo call
#  alongside the per-residue Casanovo confidence — to confirm the BLAST alignment.
# ==============================================================================

# Monoisotopic residue masses (Da).
.DDA_AA_MASS <- c(
  G=57.02146, A=71.03711, S=87.03203, P=97.05276, V=99.06841, T=101.04768,
  C=103.00919, L=113.08406, I=113.08406, N=114.04293, D=115.02694, Q=128.05858,
  K=128.09496, E=129.04259, M=131.04049, H=137.05891, F=147.06841, R=156.10111,
  Y=163.06333, W=186.07931, U=150.95364)
.DDA_PROTON <- 1.0072764665
.DDA_WATER  <- 18.0105646863
# Common Casanovo/Sage modification masses by name (case-insensitive prefix).
.DDA_MOD_MASS <- c(carbamidomethyl=57.021464, oxidation=15.994915,
  deamidated=0.984016, acetyl=42.010565, phospho=79.966331,
  ammonialoss=-17.026549, "pyro-glu"=-17.026549, carbamyl=43.005814,
  "+57.021"=57.021464, "+15.995"=15.994915, "+0.984"=0.984016)

#' Parse a (modified) peptide into per-residue (aa, modmass). Handles Casanovo
#' style `C[Carbamidomethyl]`, `+mass` tags, and an optional N-terminal mod.
dda_parse_modseq <- function(raw_seq) {
  s <- as.character(raw_seq)
  aa <- character(0); mod <- numeric(0); nterm <- 0
  i <- 1L; n <- nchar(s)
  modval <- function(tok) {
    tok <- tolower(gsub("[][]", "", tok))    # strip [ and ] (]/[ first avoids escaping)
    if (grepl("^[+-]?[0-9.]+$", tok)) return(suppressWarnings(as.numeric(tok)))
    hit <- .DDA_MOD_MASS[which(startsWith(names(.DDA_MOD_MASS), substr(tok, 1, 6)))]
    if (length(hit)) as.numeric(hit[1]) else 0
  }
  while (i <= n) {
    ch <- substr(s, i, i)
    if (ch == "[") {
      j <- i + 1L; while (j <= n && substr(s, j, j) != "]") j <- j + 1L
      tok <- substr(s, i, j)
      if (length(aa) == 0) nterm <- nterm + modval(tok)        # N-terminal mod
      else mod[length(mod)] <- mod[length(mod)] + modval(tok)  # mod on prev residue
      i <- j + 1L
    } else if (grepl("[A-Za-z]", ch)) {
      aa <- c(aa, toupper(ch)); mod <- c(mod, 0); i <- i + 1L
    } else if (ch %in% c("+", "-")) {                          # bare +mass tag
      j <- i + 1L; while (j <= n && grepl("[0-9.]", substr(s, j, j))) j <- j + 1L
      tok <- substr(s, i, j - 1L)
      if (length(mod)) mod[length(mod)] <- mod[length(mod)] + suppressWarnings(as.numeric(tok))
      i <- j
    } else i <- i + 1L
  }
  list(aa = aa, mod = mod, nterm = nterm)
}

#' b/y fragment-ion m/z for a (modified) peptide. Pure.
#' @return data.frame(label, type, idx, charge, mz)
dda_ms2_fragments <- function(raw_seq, max_charge = 2L) {
  p <- dda_parse_modseq(raw_seq)
  res <- .DDA_AA_MASS[p$aa] + p$mod
  if (anyNA(res)) res[is.na(res)] <- 0      # unknown residue -> 0 (defensive)
  L <- length(res); if (L < 2) return(data.frame())
  bsum <- cumsum(res)[-L] + p$nterm         # b1..b(L-1) neutral (sum N-term residues)
  ysum <- cumsum(rev(res))[-L] + .DDA_WATER  # y1..y(L-1): y1 = C-term residue, etc.
  out <- list()
  for (z in seq_len(max_charge)) {
    out[[length(out)+1]] <- data.frame(
      label = sprintf("b%d%s", seq_along(bsum), if (z>1) paste0("^",z,"+") else ""),
      type = "b", idx = seq_along(bsum), charge = z,
      mz = (bsum + z * .DDA_PROTON) / z, stringsAsFactors = FALSE)
    out[[length(out)+1]] <- data.frame(
      label = sprintf("y%d%s", seq_along(ysum), if (z>1) paste0("^",z,"+") else ""),
      type = "y", idx = seq_along(ysum), charge = z,
      mz = (ysum + z * .DDA_PROTON) / z, stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

#' Match observed peaks to b/y fragments within a ppm tolerance. Pure.
#' @return data.frame(mz, intensity, label, type, idx, ppm) — label NA if unmatched
dda_annotate_peaks <- function(frag_df, mz, intensity, tol_ppm = 20) {
  mz <- as.numeric(mz); intensity <- as.numeric(intensity)
  lab <- rep(NA_character_, length(mz)); typ <- lab; ion_idx <- rep(NA_integer_, length(mz))
  ppm <- rep(NA_real_, length(mz))
  if (!is.null(frag_df) && nrow(frag_df)) {
    for (i in seq_along(mz)) {
      d <- abs(frag_df$mz - mz[i]) / mz[i] * 1e6
      k <- which.min(d)
      if (length(k) && d[k] <= tol_ppm) {
        lab[i] <- frag_df$label[k]; typ[i] <- frag_df$type[k]
        ion_idx[i] <- frag_df$idx[k]; ppm[i] <- round(d[k], 1)
      }
    }
  }
  data.frame(mz = mz, intensity = intensity, label = lab, type = typ,
             idx = ion_idx, ppm = ppm, stringsAsFactors = FALSE)
}

#' Trusted de novo substitutions from a real (gapped) alignment + per-residue
#' Casanovo confidence. A "trusted substitution" is a position where the de novo
#' residue DIFFERS from the reference AND the model was confident (>= hi_conf).
#' Calibration (ocelot, vs close nr references): residues at conf>=0.95 match a
#' close reference 99% of the time — so a confident mismatch is ~99% a GENUINE
#' sequence difference (variant), not a de novo error. `bracketed` flags the
#' strongest case: the flanking aligned residues are themselves confident matches
#' (the substitution is mass-pinned on both sides by the fragment-ion ladder).
#'
#' @param qaln,saln gapped aligned query/subject (from parse_btop); '-' = gap
#' @param aa numeric per-residue Casanovo confidence for the FULL peptide (0..1)
#' @param qstart 1-based alignment start within the peptide
#' @return data.frame(pos, ref, denovo, conf, bracketed) — one row per trusted sub
build_trusted_substitutions <- function(qaln, saln, aa = NULL, qstart = 1, hi_conf = 0.95) {
  empty <- data.frame(pos = integer(0), ref = character(0), denovo = character(0),
                      conf = numeric(0), bracketed = logical(0))
  if (is.null(qaln) || is.null(saln) || !nzchar(qaln)) return(empty)
  qc <- strsplit(qaln, "")[[1]]; sc <- strsplit(saln, "")[[1]]
  n <- min(length(qc), length(sc)); if (n == 0) return(empty)
  qpos <- as.integer(qstart) - 1L
  # per aligned column: query peptide position, confidence, is-match (NA at gaps/X)
  cols <- vector("list", n)
  for (i in seq_len(n)) {
    qa <- qc[i]; sa <- sc[i]
    if (qa == "-") { cols[[i]] <- list(qpos = NA, conf = NA, match = NA, q = qa, s = sa); next }
    qpos <- qpos + 1L
    conf <- if (!is.null(aa) && qpos <= length(aa) && !is.na(aa[qpos])) aa[qpos] else NA_real_
    il <- qa %in% c("I", "L") && sa %in% c("I", "L")
    m <- if (qa == "X" || sa == "X" || sa == "-") NA else (qa == sa || il)
    cols[[i]] <- list(qpos = qpos, conf = conf, match = m, q = qa, s = sa)
  }
  conf_match <- function(j) {                 # is column j a confident match?
    if (j < 1 || j > n) return(FALSE)
    cc <- cols[[j]]; isTRUE(cc$match) && !is.na(cc$conf) && cc$conf >= hi_conf
  }
  out <- empty
  for (i in seq_len(n)) {
    cc <- cols[[i]]
    if (isTRUE(is.na(cc$match)) || isTRUE(cc$match)) next   # only mismatches with real residues
    if (is.na(cc$conf) || cc$conf < hi_conf) next            # only confident ones
    out <- rbind(out, data.frame(pos = cc$qpos, ref = cc$s, denovo = cc$q,
                 conf = round(cc$conf, 3),
                 bracketed = conf_match(i - 1) && conf_match(i + 1),
                 stringsAsFactors = FALSE))
  }
  out
}
