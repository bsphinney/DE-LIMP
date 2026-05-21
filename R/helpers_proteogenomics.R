# helpers_proteogenomics.R — Pure helper functions for the Proteogenomics feature.
# No Shiny reactivity. All functions are testable standalone.
#
# Two responsibility groups:
#  1) classify_proteins() — turns a DIA-NN report's protein table into a
#     classification data.frame (source / orf_type / parent_gene).
#  2) build_proteog_*() — produce conditional text blocks for Claude export
#     prompts (Brief / Full / Manuscript). Each returns "" when
#     !isTRUE(values$is_proteogenomics) so they can be injected unconditionally.
#
# Source-tag contract (validated 2026-05-20):
#   sp|<protein_id>|<sym>_<TAG> source=<class> ORF_type=... strand=... len=... \
#     coords=... parent_gene=... transcript=...
# Five `source` classes used downstream:
#   REF | NOVEL_GENE | NOVEL_ISOFORM | UNPARSED | UNIPROT (no source= tag) | VARIANT (prefix-detected)
#
# A %||% B → A if !is.null(A) else B. Assumed defined in helpers.R.
# If not loaded, fall back to a local definition.
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# =============================================================================
# Classification
# =============================================================================

#' Pick the column in a DIA-NN genes table that contains FASTA descriptions
#'
#' DIA-NN's column naming has drifted across versions. We try the most likely
#' description-bearing columns in order, returning the first one present.
#'
#' @param tbl data.frame — the genes table (typically `values$raw_data$genes`)
#' @return character vector aligned with `tbl$Protein.Group`, or NULL if no
#'   description column is found.
.proteog_pick_description_col <- function(tbl) {
  candidates <- c(
    "Protein.Group.Description",
    "First.Protein.Description",
    "Protein.Names",
    "Description"
  )
  for (c in candidates) {
    if (c %in% names(tbl)) return(as.character(tbl[[c]]))
  }
  NULL
}

#' Classify proteins in a DIA-NN result by proteogenomic source
#'
#' Reads the `source=` tag from the FASTA description preserved in DIA-NN's
#' output table. Entries with no `source=` tag default to UNIPROT (canonical
#' reference). VARIANT entries are detected by the INDEL_ENSP*/SNV_ENSP*
#' accession prefix (optional Phase 3 output).
#'
#' Accepts EITHER a flat genes data.frame OR a limma/limpa EList-like object
#' that has a `$genes` slot — extracts `$genes` automatically.
#'
#' @param diann_report data.frame OR list with `$genes` data.frame
#' @return data.frame with columns: Protein.Group, source, orf_type, parent_gene
classify_proteins <- function(diann_report) {
  if (is.null(diann_report)) {
    return(.empty_classification())
  }

  tbl <- if (is.data.frame(diann_report)) {
    diann_report
  } else if (is.list(diann_report) && !is.null(diann_report$genes)) {
    diann_report$genes
  } else {
    warning("classify_proteins(): input is neither a data.frame nor a list with $genes; returning empty classification")
    return(.empty_classification())
  }

  if (!"Protein.Group" %in% names(tbl) || nrow(tbl) == 0) {
    return(.empty_classification())
  }

  ids <- as.character(tbl$Protein.Group)
  desc <- .proteog_pick_description_col(tbl)

  if (is.null(desc)) {
    # No description column → cannot read source= tags; everything defaults
    # to UNIPROT except prefix-detected VARIANTs. Surface a warning so we
    # don't silently produce misleading classifications (CLAUDE.md Rule 4).
    warning("classify_proteins(): no description column found in input (tried Protein.Group.Description, First.Protein.Description, Protein.Names, Description); defaulting all non-VARIANT entries to UNIPROT class")
    desc <- rep("", length(ids))
  }

  # Default everything to UNIPROT, then overwrite from explicit source= tags.
  src <- rep("UNIPROT", length(ids))
  has_src <- grepl("source=", desc, fixed = TRUE)
  src[has_src] <- sub(".*source=([A-Z_]+).*", "\\1", desc[has_src])

  # ORF_type and parent_gene only meaningful for proteogenomic-tagged entries.
  orf <- rep(NA_character_, length(ids))
  has_orf <- grepl("ORF_type=", desc, fixed = TRUE)
  orf[has_orf] <- sub(".*ORF_type=([^[:space:]]+).*", "\\1", desc[has_orf])

  pg <- rep(NA_character_, length(ids))
  has_pg <- grepl("parent_gene=", desc, fixed = TRUE)
  pg[has_pg] <- sub(".*parent_gene=([^[:space:]]+).*", "\\1", desc[has_pg])

  # VARIANT detection: accession prefix (Phase 3 optional output).
  is_variant <- grepl("^(INDEL|SNV)_ENSP", ids)
  src[is_variant] <- "VARIANT"

  # Deduplicate to one row per unique Protein.Group (the genes table can
  # have repeats if joined upstream)
  out <- data.frame(
    Protein.Group = ids,
    source        = src,
    orf_type      = orf,
    parent_gene   = pg,
    stringsAsFactors = FALSE
  )
  out <- out[!duplicated(out$Protein.Group), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.empty_classification <- function() {
  data.frame(
    Protein.Group = character(),
    source        = character(),
    orf_type      = character(),
    parent_gene   = character(),
    stringsAsFactors = FALSE
  )
}

#' Has a proteogenomic class been observed in the classification?
#'
#' Used by the data-load handler to set `values$is_proteogenomics`.
#'
#' @param classification data.frame from classify_proteins()
#' @return logical
is_proteogenomic_session <- function(classification) {
  if (is.null(classification) || nrow(classification) == 0) return(FALSE)
  any(classification$source %in% c("REF", "NOVEL_GENE", "NOVEL_ISOFORM", "VARIANT"))
}

# =============================================================================
# Claude export prompt helpers
# =============================================================================
# Each returns "" when !isTRUE(values$is_proteogenomics).
# Strictly additive — standard UniProt sessions get an empty string and the
# template renders unchanged.

#' Counts string for the experiment-overview block
#'
#' Format:
#'   Database type: PROTEOGENOMICS-EXPANDED
#'     (18,203 canonical UniProt + 51,289 reference + 5,821 novel genes + 119 novel isoforms)
#'   → See Proteogenomics_Glossary.txt for identifier decoding.
build_proteog_note <- function(values) {
  if (!isTRUE(values$is_proteogenomics)) return("")
  pc <- values$protein_classification
  if (is.null(pc) || nrow(pc) == 0) return("")

  fmt <- function(n) format(n, big.mark = ",")
  n_uniprot <- sum(pc$source == "UNIPROT")
  n_ref     <- sum(pc$source == "REF")
  n_novel_g <- sum(pc$source == "NOVEL_GENE")
  n_novel_i <- sum(pc$source == "NOVEL_ISOFORM")
  n_variant <- sum(pc$source == "VARIANT")

  parts <- c(
    sprintf("%s canonical UniProt", fmt(n_uniprot)),
    sprintf("%s reference (Ensembl)", fmt(n_ref))
  )
  if (n_novel_g > 0) parts <- c(parts, sprintf("%s novel genes", fmt(n_novel_g)))
  if (n_novel_i > 0) parts <- c(parts, sprintf("%s novel isoforms", fmt(n_novel_i)))
  if (n_variant > 0) parts <- c(parts, sprintf("%s variant proteoforms", fmt(n_variant)))

  paste0(
    "Database type: PROTEOGENOMICS-EXPANDED (",
    paste(parts, collapse = " + "),
    ")\n",
    "→ See Proteogenomics_Glossary.txt for identifier decoding.\n"
  )
}

#' File-list line for the glossary
build_proteog_file_note <- function(values) {
  if (!isTRUE(values$is_proteogenomics)) return("")
  "- Proteogenomics_Glossary.txt — How to decode REF/NOVEL_GENE/NOVEL_ISOFORM identifiers\n"
}

#' Template-specific prompt section
#'
#' @param values reactiveValues (uses $is_proteogenomics)
#' @param template_type "brief" | "full" | "manuscript"
build_proteog_section <- function(values, template_type) {
  if (!isTRUE(values$is_proteogenomics)) return("")
  switch(template_type,
    "brief"      = .proteog_section_brief(),
    "full"       = .proteog_section_full(),
    "manuscript" = .proteog_section_manuscript(),
    ""
  )
}

.proteog_section_brief <- function() {
  paste0(
    "## 8. Proteogenomic Discoveries\n",
    "This experiment used an expanded search database — see ",
    "Proteogenomics_Glossary.txt for full guidance.\n\n",
    "Generate a SEPARATE, clearly-labeled subsection for non-canonical hits. ",
    "Do NOT mix these into the canonical \"Key Findings\" section above.\n\n",
    "**Identify proteogenomic hits**: in DE_Results_Full.csv, look at the ",
    "`source` column (or parse from Protein.Group.Description). Non-canonical ",
    "classes are NOVEL_GENE, NOVEL_ISOFORM, and (if present) VARIANT.\n\n",
    "**For the Brief, produce only**:\n",
    "- A 2-3 sentence introduction explaining that this experiment combined ",
    "standard proteomics with sample-specific RNA-seq-derived ORFs.\n",
    "- A compact table of up to 10 significant non-canonical hits, ranked by ",
    "ORF_type confidence (`complete` first), then by |log2FC|:\n",
    "  | Identifier | Class | ORF_type | log2FC | FDR | Notes |\n",
    "- The Notes column should briefly describe each: e.g. ",
    "\"complete ORF, novel intergenic locus on chromosome from parent_gene MSTRG.10029\" ",
    "or \"novel isoform of Trim25, 5'-truncated\".\n",
    "- One closing sentence: \"These hits are candidate novel/variant proteins ",
    "and warrant orthogonal validation before biological interpretation.\"\n\n",
    "**Do NOT**:\n",
    "- speculate on biological function of NOVEL_GENE identifiers\n",
    "- claim a non-canonical hit is a known protein\n",
    "- include non-canonical hits in the volcano plot top-label list — label ",
    "only canonical (REF/UNIPROT) proteins on the volcano. Mention the ",
    "non-canonical count in the volcano caption as: \"An additional N non-canonical ",
    "(proteogenomic) hits are tabulated separately.\"\n"
  )
}

.proteog_section_full <- function() {
  paste0(
    "## 5b. Proteogenomic Discoveries\n\n",
    "This experiment used a proteogenomics-expanded database. See ",
    "Proteogenomics_Glossary.txt for identifier decoding rules; treat that ",
    "file as authoritative.\n\n",

    "### 5b-i. Database Composition Summary\n",
    "Report the breakdown of identified proteins by source class. Use the ",
    "counts from the experiment overview as a starting point, and note how ",
    "many of each class were significantly regulated in the contrasts.\n\n",

    "### 5b-ii. Significant Novel Genes (NOVEL_GENE)\n",
    "For each contrast, list all NOVEL_GENE hits passing significance. Format:\n",
    "| Identifier | ORF_type | log2FC | FDR | Peptides | Notes |\n\n",
    "Sort within the table by ORF_type confidence: ",
    "`complete` first, then `5prime_partial`, then `3prime_partial`, then ",
    "`internal`. Within each tier, sort by |log2FC|.\n\n",
    "For each hit, the Notes column should briefly describe the ORF ",
    "(e.g., \"novel intergenic locus, 285 aa, on minus strand of chromosome [parse from parent_gene if possible]\") ",
    "but do NOT invent a function.\n\n",

    "### 5b-iii. Significant Novel Isoforms (NOVEL_ISOFORM)\n",
    "For each contrast, list all NOVEL_ISOFORM hits passing significance. ",
    "Critically, for each one CHECK whether the canonical reference parent ",
    "protein is ALSO present in DE_Results_Full.csv. If both are present:\n",
    "- Generate a small comparison: canonical fold-change vs novel-isoform fold-change\n",
    "- Note whether they move concordantly (likely co-regulated) or differentially ",
    "  (novel isoform may have altered stability/regulation)\n",
    "- Differential behavior is the most biologically interesting outcome\n\n",
    "Table format:\n",
    "| Novel ID | Parent Gene Symbol | log2FC (novel) | log2FC (canonical, if detected) | FDR | Interpretation |\n\n",

    "### 5b-iv. Variant Proteoforms (if present)\n",
    "Only generate this subsection if any VARIANT-class entries are in the ",
    "results. For each:\n",
    "- Parse parent ENSP, translate to gene symbol\n",
    "- Report codon position and ref→alt change\n",
    "- If canonical parent protein is ALSO present, compare fold-changes\n",
    "- For tumor samples, mention candidate neoantigen status with appropriate ",
    "  caution (requires explicit clinical context to claim)\n\n",

    "### 5b-v. Interpretation Caveats\n",
    "Write 1 paragraph noting:\n",
    "- Expanded database increases multiple-testing burden\n",
    "- Single-peptide proteogenomic hits should be flagged as provisional\n",
    "- Orthogonal validation (targeted MS, de novo sequencing, matched RNA-seq) ",
    "is recommended for any biologically significant finding\n",
    "- For NOVEL_GENE entries specifically, the lack of prior annotation is ",
    "by definition — these are discoveries that warrant follow-up, not ",
    "established proteins\n\n",

    "### 5b-vi. Figure — Class-Stratified Volcano\n",
    "Generate a volcano plot where points are colored by `source` class:\n",
    "- canonical (REF/UNIPROT): grey, smaller points\n",
    "- NOVEL_GENE: orange, slightly larger points\n",
    "- NOVEL_ISOFORM: purple, slightly larger points\n",
    "- VARIANT (if present): teal, slightly larger points\n\n",
    "Use the same axes and significance lines as the main volcano. Label only ",
    "the top 5 significant hits PER CLASS (not just by overall |logFC|), so that ",
    "non-canonical discoveries are visible even when their fold-changes are ",
    "moderate. Use ORF_type as a label-decoration if space allows (e.g., ",
    "asterisk for `complete`, dash for partial).\n\n",
    "Caption: \"Figure X. Class-stratified volcano plot. Points colored by ",
    "protein source: canonical reference (grey, n=...), StringTie novel genes ",
    "(orange, n=...), novel isoforms of annotated genes (purple, n=...). ",
    "Top 5 per non-canonical class labeled, with ORF_type indicated.\"\n"
  )
}

.proteog_section_manuscript <- function() {
  paste0(
    "### Proteogenomic Identifications\n\n",
    "[Write a Results paragraph reporting non-canonical findings. Use the ",
    "following structure verbatim, filling in numbers from DE_Results_Full.csv:]\n\n",

    "\"To detect sample-specific proteins absent from standard reference ",
    "proteomes, mass spectra were searched against a proteogenomics-expanded ",
    "database. The database combined the canonical UniProt reference proteome ",
    "with [N] StringTie-predicted open reading frames derived from sample-",
    "matched RNA-seq, comprising [N] putative novel intergenic ORFs and [N] ",
    "novel splicing isoforms of annotated genes. At an FDR threshold of ",
    "{q_cutoff}, we identified [N] significantly regulated proteogenomic hits ",
    "(Supplementary Table SX): [N] novel-gene candidates and [N] novel-isoform ",
    "candidates. Of the novel-gene candidates, [N] were structurally complete ",
    "ORFs (containing both start and stop codons). [If VARIANT entries: ",
    "Additionally, [N] sequence-variant proteoforms of canonical proteins ",
    "were identified.] These proteogenomic identifications are provisional ",
    "and would benefit from orthogonal validation by targeted mass spectrometry ",
    "or matched transcript-level evidence.\"\n\n",

    "[Write a Methods paragraph:]\n\n",

    "\"For proteogenomics analysis, the spectral search database was constructed ",
    "by concatenating [reference proteome name and version] with predicted ORFs ",
    "from StringTie (Pertea et al., 2015) assembly of matched RNA-seq data. ",
    "RNA-seq reads were quality-trimmed with fastp (Chen et al., 2018), filtered ",
    "against organism-specific ribosomal RNA with bowtie2 (Langmead and Salzberg, ",
    "2012), aligned to the reference genome with STAR (Dobin et al., 2013), ",
    "assembled with StringTie into per-sample then merged transcript models, ",
    "and classified against the reference annotation using gffcompare ",
    "(Pertea and Pertea, 2020) to distinguish novel intergenic loci from novel ",
    "isoforms of annotated genes. Open reading frames were predicted with ",
    "TransDecoder, retaining the single best ORF per transcript. ",
    "[If the search database included a contaminant library, add this sentence: ",
    "'Protein contaminant sequences were obtained from the universal protein ",
    "contaminant library (Frankenfield et al. 2022, doi:10.1021/acs.jproteome.2c00145; ",
    "https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics).'] ",
    "Differential ",
    "expression statistics for proteogenomic identifications were computed using ",
    "the same limpa/limma framework as for canonical proteins; FDR thresholds ",
    "for novel-ORF and isoform hits should be interpreted in the context of the ",
    "expanded search space.\"\n\n",

    "**Critical formatting rules**:\n",
    "- Do NOT include MSTRG.* identifiers in main manuscript text body\n",
    "- Do NOT include MSTRG.* identifiers in top-protein lists in Results figures\n",
    "- Do NOT include MSTRG.* identifiers in main volcano plot labels\n",
    "- All non-canonical hits belong in SUPPLEMENTARY tables and figures\n",
    "- Reference them as \"candidate\" identifications\n",
    "- Do NOT propose biological functions for NOVEL_GENE candidates\n",
    "- For NOVEL_ISOFORM candidates, you MAY reference the parent gene's known ",
    "function but ONLY in the context of \"the canonical form is known to ... ; ",
    "this novel isoform may have altered ... pending functional characterization\"\n",
    "- Do NOT use any non-canonical identifier in a sentence claiming biological causation\n"
  )
}

#' Section-6 append: differentiate follow-up experiments by source class
build_biosynth_proteog_note <- function(values) {
  if (!isTRUE(values$is_proteogenomics)) return("")
  paste0(
    "\n- Discuss canonical and non-canonical findings SEPARATELY. ",
    "When suggesting follow-up experiments, distinguish:\n",
    "  (a) experiments targeting well-annotated proteins (functional validation)\n",
    "  (b) experiments targeting NOVEL_GENE candidates (identity validation — ",
    "      e.g., targeted MS for the specific peptide, matched RNA-seq, ",
    "      antibody generation if peptide is long enough)\n",
    "  (c) experiments targeting NOVEL_ISOFORM candidates (functional ",
    "      comparison to canonical isoform — alternative localization, ",
    "      stability, binding partner studies)\n",
    "- If both canonical and novel-isoform forms of the same gene appear in ",
    "  results, prioritize this finding in the discussion — it's the most ",
    "  interpretable proteogenomic signal."
  )
}

#' Inline data summaries — top 15 non-canonical sig hits per contrast.
#'
#' Sorted by source class → ORF_type confidence (complete > 5p_partial >
#' 3p_partial > internal) → |logFC|.
#'
#' @param values reactiveValues (uses $is_proteogenomics, $protein_classification, $fit)
#' @param input shiny input (uses $q_cutoff, $logfc_cutoff)
#' @return character — single string with one block per contrast, or "".
build_proteog_inline <- function(values, input) {
  if (!isTRUE(values$is_proteogenomics)) return("")
  if (is.null(values$fit) || is.null(values$protein_classification)) return("")

  pc <- values$protein_classification
  noncanon_ids <- pc$Protein.Group[pc$source %in% c("NOVEL_GENE", "NOVEL_ISOFORM", "VARIANT")]
  if (length(noncanon_ids) == 0) {
    return("\n--- PROTEOGENOMIC DISCOVERIES (non-canonical significant hits) ---\nNo non-canonical entries detected in the protein classification table.\n")
  }

  contrasts <- colnames(values$fit$contrasts)
  if (is.null(contrasts) || length(contrasts) == 0) return("")

  q_cutoff   <- input$q_cutoff   %||% 0.05
  fc_cutoff  <- input$logfc_cutoff %||% 0.6

  confidence_rank <- function(ot) {
    dplyr::case_when(
      ot == "complete"       ~ 1L,
      ot == "5prime_partial" ~ 2L,
      ot == "3prime_partial" ~ 3L,
      ot == "internal"       ~ 4L,
      TRUE                   ~ 5L
    )
  }

  blocks <- lapply(contrasts, function(ct) {
    tt <- tryCatch(
      limma::topTable(values$fit, coef = ct, number = Inf),
      error = function(e) NULL
    )
    if (is.null(tt) || nrow(tt) == 0) {
      return(sprintf("Contrast %s: topTable unavailable.", ct))
    }
    tt$Protein.Group <- rownames(tt)

    sig_nc <- tt[
      tt$Protein.Group %in% noncanon_ids &
        is.finite(tt$adj.P.Val) &
        is.finite(tt$logFC) &
        tt$adj.P.Val < q_cutoff &
        abs(tt$logFC) > fc_cutoff,
      , drop = FALSE
    ]
    if (nrow(sig_nc) == 0) {
      return(sprintf("Contrast %s: no significant non-canonical hits.", ct))
    }

    sig_nc <- merge(sig_nc, pc, by = "Protein.Group", all.x = TRUE, sort = FALSE)
    sig_nc$.rank <- confidence_rank(sig_nc$orf_type)
    sig_nc <- sig_nc[order(sig_nc$source, sig_nc$.rank, -abs(sig_nc$logFC)), , drop = FALSE]
    sig_nc <- utils::head(sig_nc, 15)

    out <- data.frame(
      Identifier = sig_nc$Protein.Group,
      Source     = sig_nc$source,
      ORF_type   = sig_nc$orf_type %||% NA_character_,
      ParentGene = sig_nc$parent_gene %||% NA_character_,
      logFC      = round(sig_nc$logFC, 2),
      FDR        = signif(sig_nc$adj.P.Val, 2),
      stringsAsFactors = FALSE
    )

    paste0(
      sprintf("Contrast %s — top non-canonical significant hits (sorted by class, then ORF_type confidence):\n", ct),
      paste(utils::capture.output(print(out, row.names = FALSE)), collapse = "\n")
    )
  })

  paste0(
    "\n--- PROTEOGENOMIC DISCOVERIES (non-canonical significant hits) ---\n",
    paste(unlist(blocks), collapse = "\n\n"),
    "\n\n(See Proteogenomics_Glossary.txt for identifier decoding.)"
  )
}
