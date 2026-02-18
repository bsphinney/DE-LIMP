# ==============================================================================
#  PHOSPHO HELPERS — Pure functions for phosphosite detection & extraction
#  No Shiny reactivity — called from server_phospho.R and server_data.R
# ==============================================================================

#' Detect phosphopeptides in a DIA-NN report parquet
#'
#' Scans the Modified.Sequence column for UniMod:21 (phospho STY).
#' Returns a list with detection status, counts, and enrichment flag.
detect_phospho <- function(report_path) {
  tryCatch({
    df <- arrow::read_parquet(report_path,
      col_select = "Modified.Sequence", as_data_frame = TRUE)
    phospho_hits <- grepl("UniMod:21", df$Modified.Sequence, fixed = TRUE)
    n_phospho <- sum(phospho_hits)
    n_total   <- nrow(df)
    pct       <- round(100 * n_phospho / n_total, 1)
    list(
      detected    = n_phospho > 100,
      n_phospho   = n_phospho,
      n_total     = n_total,
      pct_phospho = pct,
      is_enriched = pct > 30
    )
  }, error = function(e) list(detected = FALSE))
}

#' Parse phosphosite positions from a DIA-NN Modified.Sequence
#'
#' Walks through the modified sequence character-by-character, tracking
#' the position in the stripped sequence.  When "(UniMod:21)" is found,
#' records the preceding residue and its position.
#'
#' @return data.frame with columns: site_residue, site_peptide_pos
parse_phospho_positions <- function(mod_seq, stripped_seq) {
  positions <- integer(0)
  residues  <- character(0)
  stripped_pos <- 0L
  i <- 1L
  n <- nchar(mod_seq)

  while (i <= n) {
    char <- substr(mod_seq, i, i)
    if (char == "(") {
      # Find closing parenthesis
      end_paren <- regexpr("\\)", substr(mod_seq, i, n))
      if (end_paren == -1L) break
      end_paren <- end_paren + i - 1L
      mod_content <- substr(mod_seq, i + 1L, end_paren - 1L)
      if (grepl("UniMod:21", mod_content, fixed = TRUE)) {
        positions <- c(positions, stripped_pos)
        residues  <- c(residues, substr(stripped_seq, stripped_pos, stripped_pos))
      }
      i <- end_paren + 1L
    } else if (grepl("[A-Z]", char)) {
      stripped_pos <- stripped_pos + 1L
      i <- i + 1L
    } else {
      i <- i + 1L
    }
  }

  if (length(positions) == 0L) {
    return(data.frame(
      site_residue    = character(0),
      site_peptide_pos = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    site_residue     = residues,
    site_peptide_pos = positions,
    stringsAsFactors = FALSE
  )
}

#' Extract phosphosites from a DIA-NN report.parquet (Path B)
#'
#' Reimplements the core sitereport algorithm in R:
#'   1. Read parquet, filter by Q.Value and UniMod:21
#'   2. Optionally filter by PTM.Site.Confidence
#'   3. Parse phospho positions from Modified.Sequence
#'   4. Expand multiply-phosphorylated peptides to one row per site
#'   5. Aggregate per SiteID x Run (max intensity = "Top 1" method)
#'   6. Pivot to matrix (sites x samples), log2 transform
#'
#' @return list(matrix, info) — site_matrix and site_info data.frame
extract_phosphosites <- function(report_path, loc_threshold = 0.75, q_cutoff = 0.01) {

  # Determine which columns are available
  available_cols <- names(arrow::read_parquet(report_path, as_data_frame = FALSE))

  required_cols <- c("Run", "Protein.Group", "Genes",
                     "Modified.Sequence", "Stripped.Sequence",
                     "Precursor.Normalised", "Q.Value")

  # Optional columns
  optional_cols <- c("PTM.Site.Confidence", "Precursor.Quantity")
  cols_to_read  <- intersect(c(required_cols, optional_cols), available_cols)

  # Fall back to Precursor.Quantity if Precursor.Normalised missing
  if (!"Precursor.Normalised" %in% available_cols && "Precursor.Quantity" %in% available_cols) {
    cols_to_read <- union(cols_to_read, "Precursor.Quantity")
  }

  df <- arrow::read_parquet(report_path, col_select = cols_to_read, as_data_frame = TRUE)

  # Use Precursor.Normalised preferentially, fall back to Precursor.Quantity
  if ("Precursor.Normalised" %in% names(df)) {
    df$Intensity <- df$Precursor.Normalised
  } else if ("Precursor.Quantity" %in% names(df)) {
    df$Intensity <- df$Precursor.Quantity
  } else {
    stop("Neither Precursor.Normalised nor Precursor.Quantity found in report.")
  }

  # Filter: Q-value and phospho precursors
  df <- df[df$Q.Value <= q_cutoff & grepl("UniMod:21", df$Modified.Sequence, fixed = TRUE), ]

  if (nrow(df) == 0) {
    return(list(matrix = NULL, info = NULL,
                message = "No phosphoprecursors passed filters."))
  }

  # Filter by localization confidence (if column exists)
  has_loc_conf <- "PTM.Site.Confidence" %in% names(df)
  if (has_loc_conf) {
    df <- df[!is.na(df$PTM.Site.Confidence) & df$PTM.Site.Confidence >= loc_threshold, ]
    if (nrow(df) == 0) {
      return(list(matrix = NULL, info = NULL,
                  message = paste0("No phosphoprecursors with localization confidence >= ",
                                   loc_threshold, ".")))
    }
  }

  # Parse phospho positions — vectorised via mapply
  sites_list <- mapply(parse_phospho_positions,
                       df$Modified.Sequence, df$Stripped.Sequence,
                       SIMPLIFY = FALSE)

  # Expand to one row per site
  n_sites_per_row <- vapply(sites_list, nrow, integer(1))
  if (sum(n_sites_per_row) == 0) {
    return(list(matrix = NULL, info = NULL,
                message = "Could not parse any phosphosite positions."))
  }

  expanded <- do.call(rbind, sites_list[n_sites_per_row > 0])
  parent_rows <- rep(seq_along(n_sites_per_row), n_sites_per_row)

  df_sites <- data.frame(
    Run              = df$Run[parent_rows],
    Protein.Group    = df$Protein.Group[parent_rows],
    Genes            = df$Genes[parent_rows],
    Intensity        = df$Intensity[parent_rows],
    Residue          = expanded$site_residue,
    Position         = expanded$site_peptide_pos,
    stringsAsFactors = FALSE
  )

  if (has_loc_conf) {
    df_sites$Loc.Conf <- df$PTM.Site.Confidence[parent_rows]
  }

  # Build SiteID: ProteinGroup_Residue_Position (peptide-relative)
  df_sites$SiteID <- paste0(df_sites$Protein.Group, "_",
                            df_sites$Residue, df_sites$Position)

  # Aggregate: max intensity per SiteID x Run ("Top 1" method)
  agg <- stats::aggregate(
    Intensity ~ SiteID + Run + Protein.Group + Genes + Residue + Position,
    data = df_sites,
    FUN  = max, na.rm = TRUE
  )

  if (has_loc_conf) {
    agg_conf <- stats::aggregate(Loc.Conf ~ SiteID, data = df_sites, FUN = max, na.rm = TRUE)
    agg <- merge(agg, agg_conf, by = "SiteID", all.x = TRUE)
  }

  # Pivot to matrix: SiteID x Run
  site_wide <- tidyr::pivot_wider(
    agg[, c("SiteID", "Run", "Intensity")],
    names_from = "Run", values_from = "Intensity"
  )
  site_ids <- site_wide$SiteID
  site_matrix <- as.matrix(site_wide[, -1, drop = FALSE])
  rownames(site_matrix) <- site_ids

  # Log2 transform
  site_matrix <- log2(site_matrix)
  site_matrix[is.infinite(site_matrix)] <- NA

  # Site metadata (one row per unique SiteID)
  info_agg <- agg[!duplicated(agg$SiteID), ]
  site_info <- data.frame(
    SiteID        = info_agg$SiteID,
    Protein.Group = info_agg$Protein.Group,
    Genes         = info_agg$Genes,
    Residue       = info_agg$Residue,
    Position      = info_agg$Position,
    stringsAsFactors = FALSE
  )
  if (has_loc_conf) {
    site_info$Best.Loc.Conf <- info_agg$Loc.Conf
  }

  list(matrix = site_matrix, info = site_info)
}

# ==============================================================================
#  Phase 2: KSEA & Motif helpers
# ==============================================================================

#' Read FASTA file and return named list of protein sequences
#'
#' Parses standard FASTA format. Extracts UniProt accession from headers
#' like ">sp|P12345|NAME_HUMAN ..." or uses the first word otherwise.
#'
#' @return Named list: accession -> amino acid sequence string
read_fasta_sequences <- function(fasta_path) {
  lines <- readLines(fasta_path, warn = FALSE)
  sequences <- list()
  current_id <- NULL
  current_seq <- character(0)

  for (line in lines) {
    if (startsWith(line, ">")) {
      if (!is.null(current_id)) {
        sequences[[current_id]] <- paste(current_seq, collapse = "")
      }
      header <- sub("^>", "", line)
      # Try UniProt format: >sp|P12345|NAME_HUMAN ... or >tr|Q9Y6K1|...
      parts <- strsplit(header, "\\|")[[1]]
      if (length(parts) >= 2) {
        current_id <- parts[2]
      } else {
        current_id <- strsplit(trimws(header), "\\s+")[[1]][1]
      }
      current_seq <- character(0)
    } else {
      current_seq <- c(current_seq, trimws(line))
    }
  }
  if (!is.null(current_id)) {
    sequences[[current_id]] <- paste(current_seq, collapse = "")
  }

  sequences
}

#' Prepare phospho DE results for KSEAapp input format
#'
#' Transforms limma topTable + site_info into the data.frame required by
#' KSEAapp::KSEA.Scores(): columns Protein, Peptide, Residue.Both, p, FC
prepare_ksea_input <- function(phospho_fit, contrast, site_info) {
  de <- limma::topTable(phospho_fit, coef = contrast, number = Inf)
  de$SiteID <- rownames(de)
  de <- merge(de, site_info, by = "SiteID")

  ksea_input <- data.frame(
    Protein      = de$Genes,
    Peptide      = de$SiteID,
    Residue.Both = paste0(de$Residue, de$Position),
    p            = de$P.Value,
    FC           = 2^de$logFC,
    stringsAsFactors = FALSE
  )

  # Remove rows with missing gene names (KSEA needs gene symbols)
  ksea_input <- ksea_input[!is.na(ksea_input$Protein) & ksea_input$Protein != "", ]
  ksea_input
}

#' Extract flanking amino acid sequences around phosphosites
#'
#' For each significant phosphosite, extracts ±window residues from the
#' FASTA protein sequence. Pads with "_" near termini. Splits by direction.
#'
#' @return list(up = character(), down = character()) of fixed-width sequences
extract_flanking_sequences <- function(site_info, de_results, fasta_sequences, window = 7L) {
  sig <- de_results[de_results$adj.P.Val < 0.05, ]
  sig$SiteID <- rownames(sig)
  sig <- merge(sig, site_info, by = "SiteID")

  if (nrow(sig) == 0) return(list(up = character(0), down = character(0)))

  flanking <- vapply(seq_len(nrow(sig)), function(i) {
    pg <- sig$Protein.Group[i]
    accession <- strsplit(pg, "[;]")[[1]][1]
    pos <- sig$Position[i]
    seq_str <- fasta_sequences[[accession]]

    if (is.null(seq_str) || is.na(pos) || pos < 1 || pos > nchar(seq_str)) {
      return(NA_character_)
    }

    start_pos <- max(1L, pos - window)
    end_pos   <- min(nchar(seq_str), pos + window)
    flank <- substr(seq_str, start_pos, end_pos)

    left_pad  <- paste(rep("_", window - (pos - start_pos)), collapse = "")
    right_pad <- paste(rep("_", window - (end_pos - pos)), collapse = "")
    paste0(left_pad, flank, right_pad)
  }, character(1))

  valid <- !is.na(flanking)
  list(
    up   = flanking[valid & sig$logFC > 0],
    down = flanking[valid & sig$logFC < 0]
  )
}

# ==============================================================================
#  Phase 3: Protein correction & AI context
# ==============================================================================

#' Correct phosphosite logFC by subtracting protein-level logFC
#'
#' Isolates phosphorylation stoichiometry changes by removing the contribution
#' of total protein abundance changes. Requires both protein-level and
#' site-level DE to have been run on the same contrast.
#'
#' @return data.frame with original + corrected logFC columns
correct_phospho_for_protein <- function(phospho_fit, protein_fit, contrast, site_info) {
  protein_de <- limma::topTable(protein_fit, coef = contrast, number = Inf)
  protein_de$Protein.Group <- rownames(protein_de)

  site_de <- limma::topTable(phospho_fit, coef = contrast, number = Inf)
  site_de$SiteID <- rownames(site_de)
  site_de <- merge(site_de, site_info[, c("SiteID", "Protein.Group"), drop = FALSE],
                   by = "SiteID")

  site_de$protein_logFC <- protein_de$logFC[
    match(site_de$Protein.Group, protein_de$Protein.Group)
  ]
  site_de$corrected_logFC <- site_de$logFC - ifelse(
    is.na(site_de$protein_logFC), 0, site_de$protein_logFC
  )

  site_de
}

#' Build phospho-specific context string for AI chat
#'
#' Appends phosphosite DE results and KSEA kinase activities to the
#' system prompt when phospho analysis is active.
phospho_ai_context <- function(phospho_fit, contrast, ksea_results = NULL) {
  de <- limma::topTable(phospho_fit, coef = contrast, number = 20)
  de$SiteID <- rownames(de)
  all_de <- limma::topTable(phospho_fit, coef = contrast, number = Inf)

  context <- paste(
    "\n\n--- PHOSPHOPROTEOMICS RESULTS ---",
    sprintf("Phosphosite DE contrast: %s", contrast),
    sprintf("Total sites tested: %d", nrow(all_de)),
    sprintf("Significant sites (FDR < 0.05): %d", sum(all_de$adj.P.Val < 0.05)),
    "\nTop 20 regulated phosphosites:",
    paste(capture.output(print(de[, c("SiteID", "logFC", "adj.P.Val")])), collapse = "\n"),
    sep = "\n"
  )

  if (!is.null(ksea_results)) {
    sig_kinases <- ksea_results[ksea_results$FDR < 0.05, , drop = FALSE]
    if (nrow(sig_kinases) > 0) {
      context <- paste(context,
        "\n\nKSEA Kinase Activity (significant, FDR < 0.05):",
        paste(capture.output(
          print(sig_kinases[, c("Kinase.Gene", "z.score", "m", "FDR")])
        ), collapse = "\n"),
        sep = "\n"
      )
    }
  }

  context
}
