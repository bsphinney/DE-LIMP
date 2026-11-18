# onoff.R -- "on/off" proteins for the MaxLFQ route: detected in one condition of a contrast and
# absent from the other. limma gives such a protein no logFC (all-NA on one side), and the
# coverage filter can drop it before limma ever sees it, so it is listed separately.
#
# ONE DEFINITION (DE-LIMP architectural rule 3): compute_onoff_proteins() below is a verbatim copy
# of DE-LIMP's compute_onoff_proteins() in R/helpers.R (added in commit b0683c04, DE-LIMP v3.9.x;
# the app calls it from R/server_data.R and R/server_de.R). Do not edit the function here
# without editing it there: tests/test_onoff_proteins.py runs both on the same fixtures and fails
# if they differ. Everything around it (contrast mapping, the CSV, the record) is the skill's.

# --- BEGIN verbatim copy of DE-LIMP R/helpers.R compute_onoff_proteins() (commit b0683c04) ---
# Surface proteins detected in >= n_min samples of one condition AND in 0
# samples of the other -- these get NA logFC from limma so they're invisible
# in the volcano. Returns a data.frame: Protein.Group, Gene, Contrast,
# Direction (one of "Group1_only" / "Group2_only"), detected_g1, detected_g2,
# total_g1, total_g2 (sample counts for each group).
compute_onoff_proteins <- function(E, group_factor, contrasts_list = NULL,
                                    n_min = 2, gene_lookup = NULL) {
  # contrasts_list can be a list of c(g1, g2) character pairs OR a 2-row matrix
  # where each column is a contrast. Normalize to list-of-pairs at function entry.
  if (is.matrix(contrasts_list)) {
    contrasts_list <- lapply(seq_len(ncol(contrasts_list)),
                              function(i) as.character(contrasts_list[, i]))
  }
  stopifnot(ncol(E) == length(group_factor))
  groups <- as.character(group_factor)
  unique_groups <- unique(groups[!is.na(groups) & nzchar(groups)])

  # If no contrasts given, generate all pairs
  if (is.null(contrasts_list) || length(contrasts_list) == 0) {
    if (length(unique_groups) < 2) return(NULL)
    pairs <- utils::combn(unique_groups, 2, simplify = FALSE)
    contrasts_list <- lapply(pairs, function(p) c(p[2], p[1]))
  }

  # Detection per cell
  detected <- !is.na(E)

  out <- list()
  for (con in contrasts_list) {
    g1 <- con[1]; g2 <- con[2]
    cols_g1 <- which(groups == g1)
    cols_g2 <- which(groups == g2)
    if (length(cols_g1) == 0 || length(cols_g2) == 0) next

    n1 <- rowSums(detected[, cols_g1, drop = FALSE])
    n2 <- rowSums(detected[, cols_g2, drop = FALSE])

    g1_only <- (n1 >= n_min) & (n2 == 0)
    g2_only <- (n2 >= n_min) & (n1 == 0)

    if (sum(g1_only) > 0) {
      df <- data.frame(
        Protein.Group = rownames(E)[g1_only],
        Contrast = paste0(g1, " - ", g2),
        Direction = paste0(g1, "_only"),
        detected_g1 = n1[g1_only],
        detected_g2 = n2[g1_only],
        total_g1 = length(cols_g1),
        total_g2 = length(cols_g2),
        stringsAsFactors = FALSE
      )
      out[[length(out) + 1]] <- df
    }
    if (sum(g2_only) > 0) {
      df <- data.frame(
        Protein.Group = rownames(E)[g2_only],
        Contrast = paste0(g1, " - ", g2),
        Direction = paste0(g2, "_only"),
        detected_g1 = n1[g2_only],
        detected_g2 = n2[g2_only],
        total_g1 = length(cols_g1),
        total_g2 = length(cols_g2),
        stringsAsFactors = FALSE
      )
      out[[length(out) + 1]] <- df
    }
  }
  if (length(out) == 0) return(NULL)
  result <- do.call(rbind, out)
  if (!is.null(gene_lookup)) {
    result$Gene <- gene_lookup[result$Protein.Group]
  }
  rownames(result) <- NULL
  result
}
# --- END verbatim copy ---

# The rule in words, recorded in de_provenance.json `onoff.rule` and read by the Methods, the
# analysis report and the AI brief (they never carry their own wording of it). `unit`: "sample"
# when the conditions name samples (Sample column), else "run"; `collapsed`: several runs of one
# sample were merged first (a sample counts as quantified if any of its runs has a value).
onoff_rule_text <- function(n_min, unit = "run", collapsed = FALSE) {
  sprintf(paste0("an on/off protein is quantified (has a value in the protein matrix, not missing) in at least %d %s%s of one ",
                 "group of a contrast and in none of the other group's %ss.%s It has no fold change ",
                 "and no p-value: this is a detection pattern, not a measured ratio. Counted on the ",
                 "protein matrix before the coverage filter, so proteins the filter drops are included."),
          n_min, unit, if (n_min == 1) "" else "s", unit,
          if (isTRUE(collapsed)) " A sample counts as quantified if any of its runs has a value (injections of one sample are collapsed)." else "")
}

# Runs -> samples, for counting on/off per SAMPLE. `sample_ids`: the conditions' Sample column, one
# per column of E (NULL / NA / blank = that run is its own sample). Returns list(E, group, collapsed,
# n_runs, n_samples): E has one column per sample, a value (1) where ANY of its runs has one, NA
# where none does; group is the sample's group (the groups of one sample's runs must agree --
# run_de.R refuses a Sample spanning groups before it gets here). Nothing to collapse (every sample
# is one run, or no Sample column): the input is returned unchanged with collapsed = FALSE.
# compute_onoff_proteins() above is NOT touched: it counts columns, and now gets sample columns.
onoff_collapse_samples <- function(E, group_factor, sample_ids = NULL) {
  ids <- if (is.null(sample_ids)) rep(NA_character_, ncol(E)) else trimws(as.character(sample_ids))
  blank <- is.na(ids) | !nzchar(ids)
  ids[blank] <- paste0("__run__", colnames(E)[blank])
  if (!anyDuplicated(ids))
    return(list(E = E, group = group_factor, collapsed = FALSE, n_runs = ncol(E), n_samples = ncol(E)))
  uniq <- unique(ids)
  det <- !is.na(E)
  Es <- vapply(uniq, function(s) as.numeric(rowSums(det[, ids == s, drop = FALSE]) > 0),
               numeric(nrow(E)))
  Es <- matrix(Es, nrow = nrow(E))
  Es[Es == 0] <- NA_real_
  dimnames(Es) <- list(rownames(E), sub("^__run__", "", uniq))
  g <- vapply(uniq, function(s) as.character(group_factor[ids == s][1]), "")
  list(E = Es, group = unname(g), collapsed = TRUE, n_runs = ncol(E), n_samples = length(uniq))
}

# Map each contrast column of limma's contrast matrix to a group pair, or say why it cannot be.
# cmat: makeContrasts() output (rows = design columns, one column per contrast); group_levels:
# the design's group names. A contrast is simple when it is exactly +1 on one group and -1 on
# another and 0 elsewhere (A-B, or -B+A); averages, interactions and anything with a covariate
# term are not (a detection pattern has no meaning for them) -> skipped, with the reason.
# Returns list(pairs = list(<name> = c(first, second)), skipped = list(<name> = reason)).
onoff_contrast_pairs <- function(cmat, group_levels) {
  pairs <- list(); skipped <- list()
  for (cn in colnames(cmat)) {
    w <- cmat[, cn]; names(w) <- rownames(cmat)
    nz <- w[!is.na(w) & w != 0]
    if (length(nz) == 2 && all(abs(nz) == 1) && sum(nz) == 0 && all(names(nz) %in% group_levels)) {
      pairs[[cn]] <- c(names(nz)[nz > 0], names(nz)[nz < 0])
    } else {
      skipped[[cn]] <- sprintf(
        "not a simple A-B contrast between two groups (%d term%s: %s); on/off needs one group against one group",
        length(nz), if (length(nz) == 1) "" else "s",
        paste(sprintf("%s%g*%s", ifelse(nz > 0, "+", ""), nz, names(nz)), collapse = " "))
    }
  }
  list(pairs = pairs, skipped = skipped)
}

# The table run_de.R writes (OnOff_Proteins.csv): compute_onoff_proteins() per simple contrast,
# in the shape the report reads. E must be the matrix BEFORE the coverage filter. `annot`: a
# data.frame with Protein.Group and Genes (or NULL). Always returns a data.frame (zero rows when
# nothing is on/off -- the file is still written, so "none found" is a result, not an absence).
onoff_table <- function(E, group_factor, pairs, n_min = 2, annot = NULL) {
  cols <- c("Protein.Group", "Genes", "Contrast", "Direction",
            "detected_g1", "detected_g2", "total_g1", "total_g2")
  res <- if (length(pairs))
    compute_onoff_proteins(E, group_factor, contrasts_list = unname(pairs), n_min = n_min) else NULL
  if (is.null(res)) {
    return(data.frame(Protein.Group = character(0), Genes = character(0), Contrast = character(0),
                      Direction = character(0), detected_g1 = integer(0), detected_g2 = integer(0),
                      total_g1 = integer(0), total_g2 = integer(0), stringsAsFactors = FALSE))
  }
  res$Genes <- NA_character_
  if (!is.null(annot) && all(c("Protein.Group", "Genes") %in% names(annot))) {
    lk <- stats::setNames(as.character(annot$Genes), as.character(annot$Protein.Group))
    res$Genes <- unname(lk[res$Protein.Group])
  }
  res[order(res$Contrast, res$Direction, -pmax(res$detected_g1, res$detected_g2)), cols]
}
