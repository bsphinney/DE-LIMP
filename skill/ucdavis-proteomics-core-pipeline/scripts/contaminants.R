# =============================================================================
# contaminants.R -- the contaminant filter run_de.R applies before quantification.
#
# Why (msalemi, 2026-09-24, Silva08172026): run_de.R re-quantifies from the DIA-NN
# report with limpa and had no contaminant filter. DIA-NN's --cont-quant-exclude
# Cont_ only shapes DIA-NN's OWN quantities, so all 121 Cont_ protein groups (bovine
# serum proteins from the antibody prep) entered the DE model and BH, and came out as
# significant hits (bovine HBB +10.5 log2 in the Kv2.1 IPs) -- while the Methods said
# contaminants "were excluded from quantification and normalisation".
#
# The rule is DIA-NN's own (README, --cont-quant-exclude): "peptides corresponding to
# protein sequence ids tagged with the specified tag will be excluded from
# normalisation as well as quantification of protein groups that do not include
# proteins with the tag". So a precursor is a contaminant when ANY accession it maps
# to (Protein.Ids) carries the tag -- not only when its inferred Protein.Group does.
# A peptide shared between a sample protein and a contaminant entry can carry the
# contaminant's signal, and DIA-NN already kept it out of that protein's quantity.
# Measured on that report (30 runs, mouse): 2,196 precursors removed, which takes out
# all 121 Cont_ protein groups, removes no sample protein group entirely, and trims
# conserved shared peptides from 47 sample protein groups (Eno1, Aldoa, Tubb3 ...).
#
# The tag mirrors fetch_fasta.py's CONT_TAG, sidecar_state() its sidecar_state(), and
# KEEP_TARGET_CONTAMINANTS_RULE / REBUILD_ADVICE its constants of those names;
# tests/test_run_de_contaminants.py asserts all of them, so a drift is a test failure
# rather than two filters that disagree.
# =============================================================================

CONTAMINANT_TAG <- "Cont_"
# Where a precursor's accessions are read from, most complete first. Protein.Ids lists
# every protein the precursor matches; Protein.Group only the inferred group. Adapted
# (protein-level) reports carry only Protein.Group.
CONTAMINANT_ID_COLUMNS <- c("Protein.Ids", "Protein.Group")

contaminant_id_column <- function(available) {
  hit <- CONTAMINANT_ID_COLUMNS[CONTAMINANT_ID_COLUMNS %in% available]
  if (length(hit)) hit[1] else NA_character_
}

# A ';'-separated accession list names a contaminant when one of its entries starts
# with the tag. Anchored, so an accession merely CONTAINING the tag does not count.
contaminant_regex <- function(tag = CONTAMINANT_TAG)
  paste0("(^|;)", gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", tag))

is_contaminant <- function(ids, tag = CONTAMINANT_TAG) {
  ids <- as.character(ids)
  !is.na(ids) & grepl(contaminant_regex(tag), ids)
}

# What the filter removes, counted per precursor (one entry per feature) so the numbers
# match what limpa quantifies. `feature` NULL = protein-level input: count rows.
contaminant_census <- function(group, is_cont, feature = NULL, genes = NULL) {
  if (is.null(feature)) feature <- seq_along(group)
  d <- data.frame(feature = feature, group = as.character(group), cont = is_cont,
                  genes = if (is.null(genes)) NA_character_ else as.character(genes),
                  stringsAsFactors = FALSE)
  d <- d[!duplicated(d$feature), , drop = FALSE]
  n      <- tapply(rep(1L, nrow(d)), d$group, sum)
  n_cont <- tapply(as.integer(d$cont), d$group, sum)
  hit    <- names(n_cont)[n_cont > 0]
  hit    <- hit[order(!is_contaminant(hit), -n_cont[hit])]
  groups <- data.frame(
    Protein.Group          = hit,
    Genes                  = d$genes[match(hit, d$group)],
    Contaminant.Group      = is_contaminant(hit),
    Precursors             = as.integer(n[hit]),
    Contaminant.Precursors = as.integer(n_cont[hit]),
    Removed.Entirely       = as.integer(n_cont[hit]) == as.integer(n[hit]),
    stringsAsFactors = FALSE)
  # A tagged group is always removed whole (its own accession carries the tag). A sample
  # group loses only the precursors it shares with a tagged entry -- all of them, rarely.
  list(n_precursors          = sum(d$cont),
       n_precursors_total    = nrow(d),
       n_groups_contaminant  = sum(groups$Contaminant.Group),
       n_groups_sample_whole = sum(!groups$Contaminant.Group & groups$Removed.Entirely),
       n_groups_sample_part  = sum(!groups$Contaminant.Group & !groups$Removed.Entirely),
       groups                = groups)
}

# The intensity the per-run share is computed on, most meaningful first. Precursor.Quantity
# is the measured signal; Precursor.Normalised has DIA-NN's RT-dependent normalisation on
# it, which does not preserve a run's signal fractions -- on Silva08172026 the contaminant
# share was median 12.7% of Precursor.Quantity but 28.3% of Precursor.Normalised.
CONTAMINANT_SHARE_COLUMNS <- c("Precursor.Quantity", "Precursor.Normalised")

# Per-run share of the precursor signal that maps to contaminant entries -- a sample
# QC number, computed whether or not the contaminants are then removed. `intensity` is
# either a features x runs matrix (run = NULL) or long vectors, one entry per feature x run.
contaminant_share <- function(intensity, is_cont, run = NULL) {
  if (is.null(run)) {
    m <- intensity
    m[!is.finite(m)] <- 0
    tot <- colSums(m)
    con <- colSums(m[is_cont, , drop = FALSE])
  } else {
    ok <- is.finite(intensity)
    s <- rowsum(cbind(intensity[ok], ifelse(is_cont[ok], intensity[ok], 0)),
                as.character(run[ok]))
    tot <- s[, 1]
    con <- s[, 2]
  }
  data.frame(Run = names(tot), Contaminant.Intensity = unname(con),
             Total.Intensity = unname(tot),
             Contaminant.Pct = round(100 * unname(con) / unname(tot), 2),
             stringsAsFactors = FALSE)
}

KEEP_TARGET_CONTAMINANTS_RULE <- "disabled (--keep-target-contaminants)"
REBUILD_ADVICE <- paste0("Rebuild the FASTA with this release's fetch_fasta.py (skill 2.8.0 or later) ",
                         "and re-search; the Core's shared human+contaminant FASTA was rebuilt with it ",
                         "on 2026-09-25 (MRS/UP000005640_9606_plus_universal_contam_2026-09.fasta).")

# Mirror of fetch_fasta.sidecar_state(): which rule built the search database.
#   "legacy"         no contaminant_target_rule: built before the overlap check
#   "identity_only"  the rule but no (or a 0) min_unique_peptides: built by the identity rule
#                    alone, so near-identical entries (bovine EEF1A1 / YWHAZ vs mouse) remain
#   "current"        both rules, --keep-target-contaminants, or no contaminants at all
sidecar_state <- function(meta) {
  truthy <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) &&
    !identical(x, FALSE) && !identical(x, 0L) && !identical(x, 0) && !identical(x, "")
  cs <- meta[["contaminant_set"]]
  holds <- truthy(meta[["n_contaminants_appended"]]) ||
    truthy(meta[["n_contaminants_already_present"]]) || (truthy(cs) && !identical(cs, "none"))
  if (!holds) return("current")
  if (!("contaminant_target_rule" %in% names(meta))) return("legacy")
  if (identical(meta[["contaminant_target_rule"]], KEEP_TARGET_CONTAMINANTS_RULE)) return("current")
  if (truthy(meta[["min_unique_peptides"]])) "current" else "identity_only"
}
sidecar_is_legacy <- function(meta) identical(sidecar_state(meta), "legacy")

# The measured near-identical set for this organism and contaminant set
# (near_identical_contaminants.json beside this file, shared with fetch_fasta.py), or NULL.
near_identical_measured <- function(meta) {
  f <- grep("^--file=", commandArgs(), value = TRUE)[1]
  dirs <- c(if (exists(".script_dir")) get(".script_dir"),
            if (!is.na(f)) dirname(normalizePath(sub("^--file=", "", f), mustWork = FALSE)),
            getwd())
  p <- file.path(dirs, "near_identical_contaminants.json")
  p <- p[file.exists(p)][1]
  if (is.na(p) || !requireNamespace("jsonlite", quietly = TRUE)) return(NULL)
  m <- tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
  tax <- suppressWarnings(as.integer(meta[["taxid"]]))
  if (is.null(m) || !identical(meta[["contaminant_set"]], m$contaminant_set) ||
      length(tax) != 1 || is.na(tax)) return(NULL)
  e <- m$by_taxid[[sprintf("%d", tax)]]
  if (is.null(e)) NULL else list(n = e$n, named = e$named, measured = m$measured)
}

# Does removing the Cont_ groups also remove real proteins of the searched organism?
# It does when the database carries contaminant entries IDENTICAL (or near-identical) to
# target proteins: DIA-NN then reports those proteins (ACTB, EEF1A1, YWHAZ, keratins ...)
# only as Cont_ groups. The sidecar records that -- a legacy or identity-only sidecar
# (sidecar_state), or a non-empty contaminants_identical_to_target_kept list.
# -> list(checked, risk, note):
# risk NA when the check could not run (note says why), FALSE when the sidecar is clean.
contaminant_database_risk <- function(sidecar, n_groups_removed) {
  if (is.null(sidecar) || !nzchar(sidecar))
    return(list(checked = FALSE, risk = NA, note = paste0(
      "not run: no FASTA sidecar (--fasta-meta) -- whether real proteins sat in the ",
      "database only as ", CONTAMINANT_TAG, " entries is unknown")))
  if (!requireNamespace("jsonlite", quietly = TRUE))
    return(list(checked = FALSE, risk = NA, note = paste0(
      "not run: jsonlite is not installed, so ", sidecar, " could not be read")))
  meta <- tryCatch(jsonlite::fromJSON(sidecar, simplifyVector = FALSE),
                   error = function(e) e)
  if (inherits(meta, "error"))
    return(list(checked = FALSE, risk = NA, note = sprintf(
      "not run: %s could not be read (%s)", sidecar, conditionMessage(meta))))
  org  <- if (is.null(meta[["organism"]]) || !nzchar(meta[["organism"]])) "target-organism"
          else meta[["organism"]]
  kept <- meta[["contaminants_identical_to_target_kept"]]
  genes <- unique(vapply(kept, function(r) {
    for (k in c("gene", "target_acc", "cont_acc"))
      if (is.character(r[[k]]) && length(r[[k]]) == 1 && nzchar(r[[k]])) return(r[[k]])
    "?"
  }, character(1)))
  state <- sidecar_state(meta)
  near <- if (identical(state, "identity_only")) near_identical_measured(meta) else NULL
  why <- c(
    if (identical(state, "legacy"))
      sprintf(paste0("the search database was built before fetch_fasta.py removed ",
                     "contaminant entries identical to %s proteins (no contaminant_target_rule ",
                     "in %s)"), org, basename(sidecar)),
    if (length(kept))
      sprintf("%d %s protein(s) are in the search database only as identical %s entries (%s)",
              length(kept), org, CONTAMINANT_TAG,
              paste(utils::head(genes, 12), collapse = ", ")),
    if (identical(state, "identity_only"))
      sprintf(paste0("the search database was built by fetch_fasta.py's identity rule alone ",
                     "(before skill 2.8.0, or --min-unique-peptides 0: %s has ",
                     "contaminant_target_rule but no min_unique_peptides), so contaminant ",
                     "entries NEAR-identical to %s proteins stayed in it%s"),
              basename(sidecar), org,
              if (is.null(near)) "" else sprintf(" -- with the %s set, %d of them: %s (measured %s)",
                                                 meta[["contaminant_set"]], near$n, near$named,
                                                 near$measured)))
  if (!length(why)) return(list(checked = TRUE, risk = FALSE, note = NULL))
  why <- paste(why, collapse = "; and ")
  list(checked = TRUE, risk = TRUE, note = sprintf(paste0(
    "%s, so some of the %d %s protein groups removed here are probably real %s proteins ",
    "(DIA-NN reports them only as %s groups), now missing from the DE. ",
    "audit_results.py --fasta-meta %s names them. %s Or re-run with --keep-contaminants to ",
    "test them together with the true contaminants."),
    why, n_groups_removed, CONTAMINANT_TAG, org, CONTAMINANT_TAG, sidecar, REBUILD_ADVICE))
}

# The record run_de.R writes into de_provenance.json ("contaminants") -- the ONE
# description of what the filter did. methods.txt below and make_methods.py read it;
# neither restates the policy on its own.
#   policy  "removed"      contaminant precursors dropped before quantification
#           "kept"         --keep-contaminants: quantified and tested with the sample
#           "none_present" no precursor maps to a tagged entry
#           "not_checked"  the report has no accession column to test
contaminant_record <- function(census, share, keep, id_column, intensity_column,
                               risk = NULL, fasta_meta = NULL,
                               share_table = "QC_contaminant_share.csv",
                               removed_table = "contaminants_removed.csv") {
  if (is.null(census))
    return(list(policy = "not_checked", removed = FALSE, tag = CONTAMINANT_TAG,
                note = "the report has neither Protein.Ids nor Protein.Group"))
  n <- census$n_precursors
  pct <- if (!is.null(share) && nrow(share)) share$Contaminant.Pct else numeric(0)
  list(
    policy = if (n == 0) "none_present" else if (keep) "kept" else "removed",
    removed = n > 0 && !keep,
    tag = CONTAMINANT_TAG,
    id_column = id_column,
    pattern = contaminant_regex(),
    rule = sprintf(paste0("a precursor is a contaminant when any accession in %s starts ",
                          "with '%s' (DIA-NN's --cont-quant-exclude rule)"),
                   id_column, CONTAMINANT_TAG),
    counted_after = "identification FDR filters, analysed runs only",
    n_precursors = n,
    n_precursors_total = census$n_precursors_total,
    n_protein_groups = census$n_groups_contaminant,
    n_sample_groups_all_shared = census$n_groups_sample_whole,
    n_sample_groups_sharing = census$n_groups_sample_part,
    intensity_column = intensity_column,
    share_median_pct = if (length(pct)) stats::median(pct) else NULL,
    share_max_pct = if (length(pct)) max(pct) else NULL,
    share_table = if (length(pct)) share_table else NULL,
    removed_table = if (n > 0 && !keep) removed_table else NULL,
    fasta_meta = fasta_meta,
    database_checked = if (is.null(risk)) NULL else risk$checked,
    database_risk = if (is.null(risk) || is.na(risk$risk)) NULL else risk$risk,
    database_note = if (is.null(risk)) NULL else risk$note)
}

# methods.txt lines, from the record.
contaminant_methods_lines <- function(rec) {
  pad <- "                "
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  share <- if (!is.null(rec$share_median_pct))
    sprintf("%sContaminant share of %s per run: median %.1f%%, max %.1f%% (%s).",
            pad, rec$intensity_column, rec$share_median_pct, rec$share_max_pct, rec$share_table)
  first <- switch(rec$policy,
    removed = c(
      sprintf("Contaminants  : REMOVED before quantification -- %s precursors mapping to a %s entry",
              fmt(rec$n_precursors), rec$tag),
      sprintf("%s(any accession in %s; DIA-NN's --cont-quant-exclude rule). %s %s protein",
              pad, rec$id_column, fmt(rec$n_protein_groups), rec$tag),
      sprintf("%sgroup(s) removed; they did not enter normalisation, the DE model or the BH", pad),
      sprintf("%scorrection. Sample protein groups that shared precursors with a %s entry:",
              pad, rec$tag),
      sprintf("%s%s lost some of them, %s lost all of them (%s).", pad,
              fmt(rec$n_sample_groups_sharing), fmt(rec$n_sample_groups_all_shared),
              rec$removed_table)),
    kept = c(
      sprintf("Contaminants  : KEPT (--keep-contaminants) -- %s precursors mapping to a %s entry",
              fmt(rec$n_precursors), rec$tag),
      sprintf("%s(%s %s protein group(s)) were quantified, normalised and tested together",
              pad, fmt(rec$n_protein_groups), rec$tag),
      sprintf("%swith the sample proteins.", pad)),
    none_present = sprintf("Contaminants  : none -- no precursor maps to a %s entry (%s).",
                           rec$tag, rec$id_column),
    sprintf("Contaminants  : NOT CHECKED -- %s", rec$note))
  c(first, share,
    if (isTRUE(rec$removed) && isTRUE(rec$database_risk))
      sprintf("%sCAUTION: %s", pad, rec$database_note)
    else if (isTRUE(rec$removed) && !is.null(rec$database_note))
      sprintf("%sDatabase check %s", pad, rec$database_note))
}
