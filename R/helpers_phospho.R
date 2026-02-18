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
