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
# The tags mirror fetch_fasta.py's CONT_TAG and CONTAMINANT_TAGS, sidecar_state() its sidecar_state(), and
# KEEP_TARGET_CONTAMINANTS_RULE / REBUILD_ADVICE its constants of those names;
# tests/test_run_de_contaminants.py asserts all of them, so a drift is a test failure
# rather than two filters that disagree.
# =============================================================================

CONTAMINANT_TAG <- "Cont_"
# Every tag that marks a contaminant accession: the skill's own (CONTAMINANT_TAG -- what
# fetch_fasta.py writes and DIA-NN's --cont-quant-exclude is given) and FragPipe's `contam_`,
# which Philosopher puts in front of the contaminant entries of a database it builds with
# --contam, and run_search.py's FragPipe DDA adapter keeps on the accession (contam_P00167;
# protein_ids.py). Knowing only Cont_, the filter let those contaminants into the DE unfiltered.
# (FragPipe's DIA route drops the tag in library.tsv, so run_search.py refuses such a database for
# FragPipe.) Both are treated alike everywhere below, a keratin sample's exemption included.
CONTAMINANT_TAGS <- c(CONTAMINANT_TAG, "contam_")

# What the SEARCH ENGINE's own quantities did with the peptides a keratin sample's kept keratin
# entries share -- said by each pipeline's descriptor (`kept_contaminant_quant`, with
# `kept_keratin_under_quantified`), never assumed to be DIA-NN's: a FragPipe, AlphaDIA, Radiant
# or Sage DE was told "DIA-NN --cont-quant-exclude" left them out. This is DIA-NN's, for the
# pipelines that read a DIA-NN report (run_de.R dpc, build_maxlfq.R undeclared maxlfq), from `cq`
# = diann_cont_quant_exclude(): the flag as the run recorded it, or NOT RECORDED (sage-review N1).
# `requantified`: dpc re-derives protein quantities from the kept precursors (TRUE); maxlfq passes
# DIA-NN's PG.MaxLFQ through (FALSE). -> list(text, under_quantified): NA = not known.
diann_kept_quant <- function(cq, requantified) {
  if (!isTRUE(cq$recorded))
    return(list(text = sprintf(paste(
      "Whether DIA-NN's own contaminant handling (--cont-quant-exclude) left those peptides",
      "out of its %s: [not recorded -- confirm] (%s)."),
      if (requantified) "normalisation" else "quantities", cq$source),
      under_quantified = if (requantified) FALSE else NA))
  if (is.na(cq$tag))
    return(list(text = sprintf(paste(
      "DIA-NN was given no --cont-quant-exclude (%s), so those peptides also took part in its",
      "own %s."), cq$source,
      if (requantified) "normalisation" else "quantities (PG.MaxLFQ, used as reported)"),
      under_quantified = FALSE))
  list(text = if (requantified) sprintf(paste(
         "The search engine's own contaminant handling (DIA-NN --cont-quant-exclude %s, per %s)",
         "had still left those peptides out of its normalisation."), cq$tag, cq$source)
       else sprintf(paste(
         "But their quantities are the search engine's own (not re-derived from precursors",
         "here), so what it left out (DIA-NN --cont-quant-exclude %s, per %s: the peptides",
         "shared with those entries) stays out: keratin is still under-quantified."),
         cq$tag, cq$source),
       under_quantified = !requantified)
}
# Where a precursor's accessions are read from, most complete first. Protein.Ids lists
# every protein the precursor matches; Protein.Group only the inferred group. Adapted
# (protein-level) reports carry only Protein.Group.
CONTAMINANT_ID_COLUMNS <- c("Protein.Ids", "Protein.Group")

contaminant_id_column <- function(available) {
  hit <- CONTAMINANT_ID_COLUMNS[CONTAMINANT_ID_COLUMNS %in% available]
  if (length(hit)) hit[1] else NA_character_
}

# A ';'-separated accession list names a contaminant when one of its entries starts
# with a tag. Anchored, so an accession merely CONTAINING a tag does not count.
contaminant_tag_group <- function(tags = CONTAMINANT_TAGS)
  paste0("(", paste(gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", tags), collapse = "|"), ")")
contaminant_regex <- function(tags = CONTAMINANT_TAGS) paste0("(^|;)", contaminant_tag_group(tags))
# One accession (already split out of a list) carries a tag.
is_tagged_accession <- function(x, tags = CONTAMINANT_TAGS)
  grepl(paste0("^", contaminant_tag_group(tags)), x)

# `exempt`: a keratin sample's keratin-family contaminant accessions (keratin_exemption()), NULL
# for any other sample -- list(shared = every one, alone = the searched species' own). A list is
# NOT a contaminant when every tagged accession in it is in `shared` AND it also names an untagged
# (sample) accession or every tagged one is in `alone`: a peptide seen only in another species'
# keratin (mouse fur, sheep wool in a human hair sample) is not the sample's. Either tag counts.
# NULL = the rule exactly as before.
is_contaminant <- function(ids, tags = CONTAMINANT_TAGS, exempt = NULL) {
  ids <- as.character(ids)
  hit <- !is.na(ids) & grepl(contaminant_regex(tags), ids)
  if (!length(exempt$shared) || !any(hit)) return(hit)
  u <- unique(ids[hit])
  keep <- vapply(strsplit(u, ";", fixed = TRUE), function(x) {
    x <- trimws(x)
    tg <- is_tagged_accession(x, tags)
    all(x[tg] %in% exempt$shared) && (any(!tg) || all(x[tg] %in% exempt$alone))
  }, logical(1))
  hit[hit] <- !keep[match(ids[hit], u)]
  hit
}

# The flagged accession lists, one per precursor (`feature`; NULL = one per row).
.flagged_ids <- function(ids, flag, feature = NULL) {
  sel <- flag & !is.na(ids)
  d <- data.frame(f = if (is.null(feature)) seq_along(ids) else feature,
                  id = as.character(ids), stringsAsFactors = FALSE)[sel, , drop = FALSE]
  d$id[!duplicated(d$f)]
}
# Flagged precursors per tag -> named integer vector (a precursor naming both counts for both).
contaminant_tag_counts <- function(ids, flag, feature = NULL, tags = CONTAMINANT_TAGS) {
  x <- .flagged_ids(ids, flag, feature)
  vapply(tags, function(t) sum(grepl(contaminant_regex(t), x)), integer(1))
}
# The tag(s) the removed accessions carried, for the record and methods ("Cont_",
# "contam_", "Cont_/contam_"); the skill's own tag when none was seen.
contaminant_tag_label <- function(counts) {
  seen <- names(counts)[counts > 0]
  if (length(seen)) paste(seen, collapse = "/") else CONTAMINANT_TAG
}

# What ONE counted item of a report is -- the one definition, read from the report's declared
# quantity level (build_maxlfq.R declared_quantity; run_de.R's dpc path reads DIA-NN precursors):
#   precursor  DIA-NN's Precursor.Id                       "precursors"
#   peptide    a peptide-level adapter's Peptide.Id (Sage)  "peptides"
#   protein    a protein-level adapter's own groups          "protein groups" (the group is the item)
#   row        an adapted report from before 2.9, which declares nothing and names no feature:
#              its rows are counted as rows, and say so (re-adapt it: run_search.py --adapt-only)
# The census counts DISTINCT items, never item x run rows: a Sage DE printed 2,744 "precursors"
# for 686 peptides x 4 runs, a FragPipe DE 288 for 72 proteins x 4 runs (sage-review, b1b94f5).
CONTAMINANT_FEATURE_COLUMNS <- c("Precursor.Id", "Peptide.Id")
contaminant_feature_column <- function(available) {
  hit <- CONTAMINANT_FEATURE_COLUMNS[CONTAMINANT_FEATURE_COLUMNS %in% available]
  if (length(hit)) hit[1] else NA_character_
}
contaminant_unit <- function(level = NULL, has_feature = TRUE) {
  u <- if (identical(level, "peptide")) c("peptide", "peptides", "peptide", "Peptides")
       else if (identical(level, "protein")) c("protein", "protein groups", "protein group",
                                                "Protein.Groups")
       else if (has_feature) c("precursor", "precursors", "precursor", "Precursors")
       else c("row", "report rows", "report row", "Rows")
  list(level = u[1], plural = u[2], singular = u[3], column = u[4])
}
# "<n> <unit>" in the unit's own words -- the singular for exactly one ("1 protein group"; a DE
# printed "1 protein groups", sage-review) -- and "____ <plural>" for a count that is missing.
count_of <- function(n, unit = contaminant_unit()) {
  if (!length(n) || is.na(n[1])) return(paste("____", unit$plural))
  paste(format(n, big.mark = ",", scientific = FALSE, trim = TRUE),
        if (n[1] == 1) unit$singular else unit$plural)
}
# the verb agreeing with count_of()'s noun
count_verb <- function(n, plural = "were", singular = "was")
  if (length(n) && !is.na(n[1]) && n[1] == 1) singular else plural

# What the filter removes, counted per item of `unit` (one entry per feature) so the numbers
# match what is quantified. `feature` NULL = count rows (unit "row").
#
# Which SAMPLE proteins lost something: for DIA-NN precursors, per inferred Protein.Group (a
# sample group loses the precursors it shares with a tagged entry). For an adapted report the
# group is a peptide's protein set (Sage) or the engine's group, so a sample protein named in a
# removed group lost that peptide/group -- counted per sample accession: "part" = keeps other
# items, "whole" = has none left. A Sage DE used to print "0 lost some" for
# 101 mixed Cont_ groups (sage-review): the group-level count cannot see them there.
contaminant_census <- function(group, is_cont, feature = NULL, genes = NULL, exempt = NULL,
                               unit = contaminant_unit()) {
  if (is.null(feature)) feature <- seq_along(group)
  d <- data.frame(feature = feature, group = as.character(group), cont = is_cont,
                  genes = if (is.null(genes)) NA_character_ else as.character(genes),
                  stringsAsFactors = FALSE)
  d <- d[!duplicated(d$feature), , drop = FALSE]
  n      <- tapply(rep(1L, nrow(d)), d$group, sum)
  n_cont <- tapply(as.integer(d$cont), d$group, sum)
  hit    <- names(n_cont)[n_cont > 0]
  hit    <- hit[order(!is_contaminant(hit, exempt = exempt), -n_cont[hit])]
  cgroup <- is_contaminant(hit, exempt = exempt)
  mixed  <- cgroup & vapply(strsplit(hit, ";", fixed = TRUE),
                            function(x) any(!is_tagged_accession(trimws(x))), logical(1))
  groups <- data.frame(
    Protein.Group          = hit,
    Genes                  = d$genes[match(hit, d$group)],
    Contaminant.Group      = cgroup,
    n_all                  = as.integer(n[hit]),
    n_cont                 = as.integer(n_cont[hit]),
    Removed.Entirely       = as.integer(n_cont[hit]) == as.integer(n[hit]),
    stringsAsFactors = FALSE)
  names(groups)[4:5] <- c(unit$column, paste0("Contaminant.", unit$column))
  if (identical(unit$level, "precursor")) {
    # A tagged group is always removed whole (its own accession carries the tag). A sample
    # group loses only the precursors it shares with a tagged entry -- all of them, rarely.
    part  <- sum(!groups$Contaminant.Group & !groups$Removed.Entirely)
    whole <- sum(!groups$Contaminant.Group & groups$Removed.Entirely)
    scope <- "protein groups"
  } else {
    untag <- lapply(strsplit(d$group, ";", fixed = TRUE), function(x) {
      x <- trimws(x); x[nzchar(x) & !is_tagged_accession(x)] })
    rm_acc <- unique(unlist(untag[d$cont]))
    kept   <- unique(unlist(untag[!d$cont]))
    whole <- length(setdiff(rm_acc, kept))
    # peptide level: a sample protein that keeps other peptides lost some. Protein
    # level: the sharing happened inside the engine's own protein quantities, which this report
    # does not show -- NOT computed (NA, said so), never a 0 by default (rule 2)
    part  <- if (identical(unit$level, "protein")) NA_integer_ else length(intersect(rm_acc, kept))
    scope <- "proteins"
    if (identical(unit$level, "row")) part <- whole <- NA_integer_   # no feature: not computed
  }
  list(n_precursors          = sum(d$cont),
       n_precursors_total    = nrow(d),
       n_groups_contaminant  = sum(groups$Contaminant.Group),
       n_groups_mixed        = sum(mixed),
       n_groups_sample_whole = whole,
       n_groups_sample_part  = part,
       sample_scope          = scope,
       unit                  = unit,
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

# A file shipped beside this one (fetch_fasta.py, its JSON data), or NA.
contaminants_sibling <- function(name) {
  f <- grep("^--file=", commandArgs(), value = TRUE)[1]
  dirs <- c(if (exists(".script_dir")) get(".script_dir"),
            if (!is.na(f)) dirname(normalizePath(sub("^--file=", "", f), mustWork = FALSE)),
            getwd())
  p <- file.path(dirs, name)
  p[file.exists(p)][1]
}

# The measured near-identical set for this organism and contaminant set
# (near_identical_contaminants.json beside this file, shared with fetch_fasta.py), or NULL.
near_identical_measured <- function(meta) {
  p <- contaminants_sibling("near_identical_contaminants.json")
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
# With NO sidecar, `search_dir` (the report's folder) lets fetch_fasta.py check the FASTA the
# search itself names -- database_check_without_sidecar() below.
contaminant_database_risk <- function(sidecar, n_groups_removed, search_dir = NULL) {
  if (is.null(sidecar) || !nzchar(sidecar))
    return(database_check_without_sidecar(search_dir, n_groups_removed))
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
      sprintf(paste0("%d %s protein(s) are in the search database only as identical (or ",
                     "peptide-indistinguishable) %s entries (%s)"),
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
    "(the search reports them only as %s groups), now missing from the DE. ",
    "audit_results.py --fasta-meta %s names them. %s Or re-run with --keep-contaminants to ",
    "test them together with the true contaminants."),
    why, n_groups_removed, CONTAMINANT_TAG, org, CONTAMINANT_TAG, sidecar, REBUILD_ADVICE))
}

# The search's own folder -- where search_provenance.json and the engine log are: the report's
# folder, or its parent when only that holds search_provenance.json (FragPipe's DIA route writes
# its report to <workdir>/dia-quant-output/). The report's folder when neither does.
search_folder <- function(input) {
  d <- dirname(normalizePath(input, mustWork = FALSE))
  for (cand in c(d, dirname(d)))
    if (file.exists(file.path(cand, "search_provenance.json"))) return(cand)
  d
}

# <script> <sub> <args>'s JSON answer -> list(res, why): res NULL when there is none, and `why`
# then says why (python3, the script or jsonlite missing; no parseable output). The checks
# themselves are defined once, in the Python script; this only asks.
python_json <- function(script, sub, args, purpose) {
  none <- function(why) list(res = NULL, why = why)
  py <- Sys.which("python3")
  path <- contaminants_sibling(script)
  if (!nzchar(py)) return(none(paste("python3 is not on PATH to", purpose)))
  if (is.na(path)) return(none(sprintf("%s was not found beside contaminants.R", script)))
  if (!requireNamespace("jsonlite", quietly = TRUE))
    return(none(sprintf("jsonlite is not installed to read %s's answer", script)))
  out <- suppressWarnings(tryCatch(
    system2(py, c(shQuote(path), sub, args), stdout = TRUE, stderr = FALSE, timeout = 600),
    error = function(e) structure(character(0), status = conditionMessage(e))))
  res <- tryCatch(jsonlite::fromJSON(paste(out, collapse = "\n"), simplifyVector = TRUE),
                  error = function(e) NULL)
  if (is.null(res))
    return(none(sprintf("%s %s gave no answer (exit status %s)", script, sub,
                        if (is.null(attr(out, "status"))) "0"
                        else paste(attr(out, "status"), collapse = ""))))
  list(res = res, why = NULL)
}
fetch_fasta_json <- function(sub, args, purpose) python_json("fetch_fasta.py", sub, args, purpose)

# Did DIA-NN run with --cont-quant-exclude for this report? Asked of make_methods.py, THE reader
# of that flag (the command line in the DIA-NN log that wrote the report, then the parameters the
# search ran with). -> list(recorded, tag, source): tag NA when the flag was absent; recorded
# FALSE -- `source` then says why -- when nothing could be read. Never assumed (sage-review N1):
# FragPipe's DIA-NN step, or a DIA-NN run by hand, may not have been given it.
diann_cont_quant_exclude <- function(report_path) {
  a <- python_json("make_methods.py", "cont-quant-exclude",
                   shQuote(normalizePath(report_path, mustWork = FALSE)),
                   "read DIA-NN's --cont-quant-exclude")
  if (is.null(a$res)) return(list(recorded = FALSE, tag = NA_character_, source = a$why))
  r <- a$res
  if (!isTRUE(r$recorded))
    return(list(recorded = FALSE, tag = NA_character_,
                source = "no DIA-NN log beside the report and no DIA-NN parameters file was found"))
  list(recorded = TRUE, tag = if (length(r$value) && nzchar(r$value)) r$value else NA_character_,
       source = r$source)
}

# No sidecar: ask fetch_fasta.py about the FASTA the search itself names. Why (release review
# 2026-09-28): this used to say only "not run", so a search on the Core's superseded MRS human
# FASTA -- the Siegel entries staged 2026-09-25 -- lost ACTB, EEF1A1, KRT8 and ~150 other real
# proteins with no caveat. The check (find the FASTA in the search's provenance / log, re-check
# it, or recognise a known superseded database by md5 or name) is defined ONCE, in
# fetch_fasta.database_without_sidecar(); this only words its answer. Same list(checked, risk,
# note) as contaminant_database_risk().
database_check_without_sidecar <- function(search_dir, n_groups_removed) {
  not_run <- function(why) list(checked = FALSE, risk = NA, note = paste0(
    "not run: no FASTA sidecar (--fasta-meta), and ", why, " -- whether real proteins sat in ",
    "the database only as ", CONTAMINANT_TAG, " entries is unknown"))
  if (is.null(search_dir) || !nzchar(search_dir))
    return(not_run("the search's folder is not known"))
  q <- fetch_fasta_json("check-db", c("--search-dir", shQuote(search_dir)),
                        "check the search's FASTA")
  if (is.null(q$res)) return(not_run(q$why))
  res <- q$res
  if (!isTRUE(res$checked)) return(not_run(res$why))
  if (!isTRUE(res$risk))
    return(list(checked = TRUE, risk = FALSE, note = paste0("passed: ", res$why)))
  org <- if (is.character(res$organism) && nzchar(res$organism)) res$organism
         else "target-organism"
  genes <- as.character(unlist(res$genes))
  list(checked = TRUE, risk = TRUE, note = sprintf(paste0(
    "%s -- %s%s. So some of the %d %s protein groups removed here are probably real %s ",
    "proteins (the search reports them only as %s groups), now missing from the DE. %s Or re-run ",
    "with --keep-contaminants to test them together with the true contaminants."),
    res$why, paste(utils::head(genes, 15), collapse = ", "),
    if (length(genes) > 15) sprintf(", and %d more", length(genes) - 15) else "",
    n_groups_removed, CONTAMINANT_TAG, org, CONTAMINANT_TAG, REBUILD_ADVICE))
}

# Is this a keratin sample (hair, wool, feather, skin, nail ...: keratin is the ANALYTE), and
# which contaminant entries must the filter then leave in? Why (msalemi's hair benchmark SET28,
# skill 2.8.0): every precursor naming a Cont_ entry was removed, and ~46 keratin-family Cont_
# entries (KRT34, KRTAPs, mouse hair and sheep wool keratins) took the hair's own keratin
# peptides with them. Read, first that says so: run_de.R --keratin-sample; the FASTA sidecar's
# keratin_sample (fetch_fasta.py fetch --keratin-sample -- keratin_sample_recorded() there is the
# one reader); search_provenance.json's keratin_sample (run_search.py). With no --fasta-meta, the
# sidecar of the FASTA search_provenance.json names is read.
#   database "keratins_removed_at_build"  built with --keratin-sample: no keratin-family
#                                          contaminant entry is left, so nothing needs exempting
#            "keratins_in_database"       not built that way: the keratin-family entries
#                                          fetch_fasta.py keratin-db finds are exempted
#            "no_keratin_contaminants"    checked: the database holds none
#            "unknown"                    they could not be found: nothing is exempted -- the
#                                          filter removes them as before, and the note says so
# value FALSE: recorded as not a keratin sample. NA: recorded nowhere -- today's behaviour, and
# the note says so. `sample_source`: "user" when step 3's question was answered (a flag, or the
# sidecar's keratin_sample_source), "default" when nobody answered. Unless the user said "no", the
# database's keratin entries are looked up (keratin_accessions) so keratin_default_removed() can
# say how many keratin precursors the filter took on that unconfirmed default (rule 2).
# -> list(value, source, sample_source, database, exempt, exempt_alone, note, ...).
keratin_sample_status <- function(flag = FALSE, sidecar = NULL, search_dir = NULL) {
  read <- function(p) if (!is.null(p) && length(p) == 1 && nzchar(p) && file.exists(p) &&
                          requireNamespace("jsonlite", quietly = TRUE))
    tryCatch(jsonlite::fromJSON(p, simplifyVector = FALSE), error = function(e) NULL)
  lgl <- function(x) if (is.logical(x) && length(x) == 1 && !is.na(x)) x else NA
  answered <- function(x) if (is.character(x) && length(x) == 1 && nzchar(x)) x else "default"
  sp_path <- if (!is.null(search_dir) && nzchar(search_dir))
    file.path(search_dir, "search_provenance.json")
  sp <- read(sp_path)
  if (is.null(sidecar) && is.character(sp[["fasta"]]) && length(sp[["fasta"]]) == 1 &&
      file.exists(paste0(sp[["fasta"]], ".meta.json")))
    sidecar <- paste0(sp[["fasta"]], ".meta.json")
  meta <- read(sidecar)
  in_meta <- lgl(meta[["keratin_sample"]])
  in_sp   <- lgl(sp[["keratin_sample"]])
  value <- if (isTRUE(flag) || isTRUE(in_meta) || isTRUE(in_sp)) TRUE
           else if (identical(in_meta, FALSE) || identical(in_sp, FALSE)) FALSE else NA
  from_meta <- !isTRUE(flag) && !is.na(in_meta) && (isTRUE(in_meta) || !isTRUE(in_sp))
  from_sp   <- !isTRUE(flag) && !from_meta && !is.na(in_sp)
  rec <- list(value = value,
              source = if (isTRUE(flag)) "--keratin-sample" else if (from_meta) sidecar
                       else if (from_sp) sp_path else NULL,
              sample_source = if (isTRUE(flag)) "user"
                              else if (from_meta) answered(meta[["keratin_sample_source"]])
                              else if (from_sp) answered(sp[["keratin_sample_source"]])
                              else NA_character_,
              exempt = character(0), exempt_alone = character(0))
  if (is.na(value) || (!value && !identical(rec$sample_source, "user"))) {
    if (is.na(value))
      rec$note <- paste0(
        "whether these samples are keratin (hair, wool, feather, skin, nail ...) is not recorded ",
        "(no keratin_sample in the FASTA sidecar or search_provenance.json, and no ",
        "--keratin-sample), so keratin-family contaminant entries were treated as ",
        "contaminants, as for any sample. For a keratin sample, re-run with --keratin-sample.")
    listed <- meta[["keratin_contaminants_in_database"]]
    if (is.list(listed)) {
      rec$keratin_accessions <- as.character(unlist(listed))
    } else {
      q <- fetch_fasta_json("keratin-db",
                            c(if (!is.null(sidecar)) c("--fasta-meta", shQuote(sidecar)),
                              if (!is.null(search_dir)) c("--search-dir", shQuote(search_dir))),
                            "find the search database's keratin-family entries")
      if (!is.null(q$res) && isTRUE(q$res$checked))
        rec$keratin_accessions <- as.character(unlist(q$res$accessions))
    }
    return(rec)
  }
  if (!value) return(rec)
  # `tags_checked`: the contaminant tags the database that was checked holds. A tag the report
  # carries beyond them (contam_ entries from a database the skill did not build) was never
  # checked for keratins -- keratin_unchecked() below says so. A record without the field was
  # written by fetch_fasta.py / run_search.py before it existed; fetch_fasta.py tags what it
  # writes with CONTAMINANT_TAG, so that one is known.
  tags_of <- function(x) if (is.null(x)) CONTAMINANT_TAG else as.character(unlist(x))
  if (isTRUE(in_meta)) {
    n <- meta[["n_contaminants_dropped_keratin_sample"]]
    rec$database <- "keratins_removed_at_build"
    rec$n_removed_at_build <- if (is.numeric(n)) as.integer(n) else NULL
    rec$tags_checked <- tags_of(meta[["contaminant_tags_in_database"]])
    return(rec)
  }
  chk <- sp[["keratin_sample_check"]]
  if (isTRUE(in_sp) && identical(as.integer(chk[["keratin_contaminants_in_database"]]), 0L)) {
    rec$database <- "no_keratin_contaminants"
    rec$checked_in <- sp_path
    rec$tags_checked <- tags_of(chk[["contaminant_tags_in_database"]])
    return(rec)
  }
  q <- fetch_fasta_json("keratin-db",
                        c(if (!is.null(sidecar)) c("--fasta-meta", shQuote(sidecar)),
                          if (!is.null(search_dir)) c("--search-dir", shQuote(search_dir))),
                        "find the search database's keratin-family entries")
  why <- if (is.null(q$res)) q$why else if (!isTRUE(q$res$checked)) q$res$why
  if (!is.null(why)) {
    rec$database <- "unknown"
    rec$caution <- TRUE
    rec$note <- paste0(
      "a keratin sample, but which contaminant (", paste(CONTAMINANT_TAGS, collapse = "/"),
      ") entries of the search database are keratin-family could not be determined (", why,
      "), so what maps to them was REMOVED as contaminants, as for any sample: the ",
      "sample's keratins are under-counted. Rebuild the FASTA with fetch_fasta.py fetch ",
      "--keratin-sample and re-search.")
    return(rec)
  }
  rec$exempt <- as.character(unlist(q$res$accessions))
  # the searched species' own keratin entries: a precursor mapping only to one of these is the
  # sample's; one mapping only to another species' keratin (mouse fur, sheep wool in a human hair
  # sample) is kept only when it also names a sample (untagged) protein
  rec$exempt_alone <- as.character(unlist(q$res$own_species))
  rec$taxid <- q$res$taxid
  rec$database <- if (length(rec$exempt)) "keratins_in_database" else "no_keratin_contaminants"
  rec$checked_in <- q$res$source
  rec$tags_checked <- as.character(unlist(q$res$tags))      # read from the FASTA itself
  if (length(rec$exempt)) {
    rec$caution <- TRUE
    rec$note <- paste0(
      "the search database was not built with fetch_fasta.py --keratin-sample and holds ",
      length(rec$exempt), " keratin-family contaminant entries (", length(rec$exempt_alone),
      " of them the searched organism's own): what maps only to them was kept as sample ",
      "protein (another species' entry only when it also names a sample protein); what the ",
      "search engine's own quantities did with those peptides is stated with ",
      "the kept count. Rebuild the FASTA with fetch_fasta.py fetch --keratin-sample and re-search ",
      "for a database without them.")
  }
  rec
}

# The exemption is_contaminant() takes for a keratin record: NULL when nothing is exempt.
keratin_exemption <- function(k)
  if (length(k$exempt)) list(shared = k$exempt, alone = k$exempt_alone) else NULL

# A keratin sample whose report carries a contaminant tag the database check did not cover:
# contam_ entries come from a database built with FragPipe/Philosopher's --contam, not by the
# skill (run_search.py refuses such a database for FragPipe), so when the database checked holds
# none, which of them are keratins cannot be told. They are removed as contaminants (nothing to
# exempt them by), and the record says so -- never a silent loss.
keratin_unchecked <- function(k, ids, flag, feature = NULL, unit = contaminant_unit()) {
  # "unknown": nothing was checked, and its note already says every keratin entry was removed
  if (!isTRUE(k$value) || identical(k$database, "unknown") || !any(flag)) return(k)
  miss <- setdiff(CONTAMINANT_TAGS, k$tags_checked)
  if (!length(miss)) return(k)
  x <- .flagged_ids(ids, flag, feature)
  hit <- miss[vapply(miss, function(t) any(grepl(contaminant_regex(t), x)), logical(1))]
  if (!length(hit)) return(k)
  n <- sum(grepl(contaminant_regex(hit), x))
  k$unchecked_tags <- as.list(hit)
  k$n_precursors_unchecked <- n
  k$caution <- TRUE
  k$note <- paste(c(k$note, sprintf(paste0(
    "%s %s to %s-tagged contaminant entries the search database checked here does ",
    "not hold -- contaminants the skill did not add (a database built with FragPipe/Philosopher's ",
    "--contam), never checked for keratins: they were removed as contaminants, and the sample's ",
    "keratins may be under-counted."), count_of(n, unit), count_verb(n, "map", "maps"),
    paste(hit, collapse = "/"))), collapse = " Also: ")
  k
}

# Not a keratin sample by the user's answer: nothing to say. Otherwise -- not recorded, or a
# default nobody confirmed -- count the removed precursors that map to a keratin-family entry of
# the database (keratin_accessions): if any, the methods say they went on that default (rule 2).
keratin_default_removed <- function(k, ids, flag, feature = NULL) {
  if (isTRUE(k$value) || identical(k$sample_source, "user") || !length(k$keratin_accessions) ||
      !any(flag)) return(k)
  x <- .flagged_ids(ids, flag, feature)
  n <- sum(vapply(strsplit(x, ";", fixed = TRUE),
                  function(t) any(trimws(t) %in% k$keratin_accessions), logical(1)))
  if (n > 0) {
    k$default_removed <- TRUE
    k$n_keratin_precursors_removed <- n
  }
  k
}

# The record run_de.R writes into de_provenance.json ("contaminants") -- the ONE
# description of what the filter did. methods.txt below and make_methods.py read it;
# neither restates the policy on its own.
#   policy  "removed"      contaminant precursors dropped before quantification
#           "kept"         --keep-contaminants: quantified and tested with the sample
#           "none_present" no precursor maps to a tagged entry
#           "not_checked"  the report has no accession column to test
#   keratin_sample  keratin_sample_status(), plus n_precursors_kept (precursors mapping only to
#                   exempted keratin-family entries), requantified (the pipeline's descriptor:
#                   are protein quantities re-derived from precursors here?), and for a sample
#                   not confirmed as non-keratin, default_removed / n_keratin_precursors_removed:
#                   what a keratin sample kept, or why not
contaminant_record <- function(census, share, keep, id_column, intensity_column,
                               risk = NULL, fasta_meta = NULL,
                               share_table = "QC_contaminant_share.csv",
                               removed_table = "contaminants_removed.csv",
                               keratin = NULL, tag = CONTAMINANT_TAG) {
  ker <- if (!is.null(keratin)) {
    k <- keratin
    k$exempt_accessions <- as.list(k$exempt)     # always a JSON array
    k$exempt_alone_accessions <- as.list(k$exempt_alone)
    k$n_keratin_accessions <- if (!is.null(k$keratin_accessions)) length(k$keratin_accessions)
    k$exempt <- k$exempt_alone <- k$keratin_accessions <- NULL
    k
  }
  if (is.null(census))
    return(list(policy = "not_checked", removed = FALSE, tag = CONTAMINANT_TAG,
                note = "the report has neither Protein.Ids nor Protein.Group",
                keratin_sample = ker))
  n <- census$n_precursors
  u <- census$unit
  pct <- if (!is.null(share) && nrow(share)) share$Contaminant.Pct else numeric(0)
  list(
    policy = if (n == 0) "none_present" else if (keep) "kept" else "removed",
    removed = n > 0 && !keep,
    # the tag(s) the removed accessions carried (contaminant_tag_label); `tags`, every tag the
    # rule tests
    tag = tag,
    tags = as.list(CONTAMINANT_TAGS),
    id_column = id_column,
    pattern = contaminant_regex(),
    # what one counted item is (contaminant_unit): every count below is of distinct items
    unit = u$plural, unit_singular = u$singular, unit_level = u$level,
    # THE rule text: methods.txt and make_methods.py both quote it
    rule = sprintf(paste0("a %s is a contaminant when any accession in %s starts ",
                          "with %s -- the rule of DIA-NN's --cont-quant-exclude, applied here"),
                   u$singular, id_column, paste0("'", CONTAMINANT_TAGS, "'", collapse = " or ")),
    counted_after = "identification FDR filters, analysed runs only",
    n_precursors = n,
    n_precursors_total = census$n_precursors_total,
    n_protein_groups = census$n_groups_contaminant,
    n_protein_groups_mixed = census$n_groups_mixed,
    # NA (JSON null) = not computed for this unit -- never a 0 by default (rule 2)
    n_sample_groups_all_shared = census$n_groups_sample_whole,
    n_sample_groups_sharing = census$n_groups_sample_part,
    sample_scope = census$sample_scope,
    sample_loss = contaminant_sample_loss(census, tag),
    intensity_column = intensity_column,
    share_median_pct = if (length(pct)) stats::median(pct) else NULL,
    share_max_pct = if (length(pct)) max(pct) else NULL,
    share_table = if (length(pct)) share_table else NULL,
    removed_table = if (n > 0 && !keep) removed_table else NULL,
    fasta_meta = fasta_meta,
    database_checked = if (is.null(risk)) NULL else risk$checked,
    database_risk = if (is.null(risk) || is.na(risk$risk)) NULL else risk$risk,
    database_note = if (is.null(risk)) NULL else risk$note,
    keratin_sample = ker)
}

# THE sentence on which sample proteins lost something to the filter, per unit (census). methods.txt
# and make_methods.py both quote the record's copy.
contaminant_sample_loss <- function(census, tag) {
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  part <- census$n_groups_sample_part
  whole <- census$n_groups_sample_whole
  switch(census$unit$level,
    precursor = sprintf(paste0("Sample protein groups that shared precursors with a %s entry: ",
                               "%s lost some of them, %s lost all of them."),
                        tag, fmt(part), fmt(whole)),
    # "keep other peptides": in any kept group, their own or one shared with other proteins (on
    # gabrig's HeL50 Sage run 57 = 36 as their own group + 21 only inside a shared one)
    peptide = sprintf(paste0("Sample proteins named by a removed peptide (a peptide shared with a ",
                             "%s entry): %s keep other peptides, %s have no other peptide and are ",
                             "not quantified."), tag, fmt(part), fmt(whole)),
    protein = sprintf(paste0("%s sample protein(s) were removed because their protein group also ",
                             "named a %s entry; whether sample proteins lost peptides shared with a ",
                             "%s entry inside the search engine's own protein quantities is not ",
                             "computed (a protein-level report)."), fmt(whole), tag, tag),
    paste0("Which sample proteins lost report rows is not computed: this adapted report names no ",
           "precursor or peptide (made before skill 2.9; re-adapt it with run_search.py ",
           "--adapt-only)."))
}

# methods.txt lines, from the record.
contaminant_methods_lines <- function(rec) {
  pad <- "                "
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  # the record's unit (contaminant_unit); a record from before 2.9 counted DIA-NN precursors
  unit <- list(plural = if (is.null(rec$unit)) "precursors" else rec$unit,
               singular = if (is.null(rec$unit_singular)) "precursor" else rec$unit_singular)
  share <- if (!is.null(rec$share_median_pct))
    sprintf("%sContaminant share of %s per run: median %.1f%%, max %.1f%% (%s).",
            pad, rec$intensity_column, rec$share_median_pct, rec$share_max_pct, rec$share_table)
  mixed <- if (isTRUE(rec$n_protein_groups_mixed > 0))
    sprintf(" (%s of them also named a sample protein)", fmt(rec$n_protein_groups_mixed)) else ""
  first <- switch(rec$policy,
    removed = c(
      sprintf("Contaminants  : REMOVED before quantification -- %s mapping to a %s entry",
              count_of(rec$n_precursors, unit), rec$tag),
      sprintf("%s(%s). %s %s protein group(s) removed%s; they did not enter normalisation, the",
              pad, rec$rule, fmt(rec$n_protein_groups), rec$tag, mixed),
      sprintf("%sDE model or the BH correction.", pad),
      sprintf("%s%s (%s)", pad, rec$sample_loss, rec$removed_table)),
    kept = c(
      sprintf("Contaminants  : KEPT (--keep-contaminants) -- %s mapping to a %s entry",
              count_of(rec$n_precursors, unit), rec$tag),
      sprintf("%s(%s %s protein group(s)) were quantified, normalised and tested together",
              pad, fmt(rec$n_protein_groups), rec$tag),
      sprintf("%swith the sample proteins.", pad)),
    none_present = sprintf("Contaminants  : none -- no %s maps to a %s entry (%s).",
                           unit$singular, rec$tag, rec$id_column),
    sprintf("Contaminants  : NOT CHECKED -- %s", rec$note))
  c(first, share,
    if (isTRUE(rec$removed) && isTRUE(rec$database_risk))
      sprintf("%sCAUTION: %s", pad, rec$database_note)
    else if (isTRUE(rec$removed) && !is.null(rec$database_note))
      sprintf("%sDatabase check %s", pad, rec$database_note),
    keratin_methods_lines(rec$keratin_sample, pad, unit))
}

# methods.txt lines for a keratin sample (value TRUE); for a sample whose "not keratin" nobody
# confirmed, ONE tagged line when keratin precursors were removed on that default (rule 2).
# Nothing for a sample the user said is not keratin -- its lines above are exactly what they were
# before keratin samples were handled.
keratin_methods_lines <- function(k, pad = strrep(" ", 16), unit = contaminant_unit()) {
  if (is.null(k)) return(NULL)
  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  if (!isTRUE(k$value))
    return(if (isTRUE(k$default_removed)) sprintf(paste0(
      "Keratin sample: not confirmed -- %s mapping to keratin-family contaminant ",
      "entries %s removed as contaminants; whether samples are keratinous was not recorded ",
      "(DEFAULT \u2014 not user-confirmed)."), count_of(k$n_keratin_precursors_removed, unit),
      count_verb(k$n_keratin_precursors_removed)))
  head <- sprintf("Keratin sample: yes (%s) -- keratin is the analyte, not a contaminant.",
                  if (is.null(k$source)) "source not recorded" else k$source)
  kept <- count_of(k$n_precursors_kept, unit)
  were <- count_verb(k$n_precursors_kept)
  body <- switch(if (is.null(k$database)) "" else k$database,
    keratins_removed_at_build = sprintf(paste0(
      "%sThe FASTA was built with fetch_fasta.py --keratin-sample, which removed its %s ",
      "keratin-family contaminant entries, so keratin %s were quantified and tested as ",
      "sample proteins; every other contaminant entry (trypsin, BSA ...) was handled as above."),
      pad, if (is.null(k$n_removed_at_build)) "____" else fmt(k$n_removed_at_build),
      unit$plural),
    # Whether keeping them changes a quantity is the pipeline's to say (its descriptor's
    # requantified_from_precursors): re-quantified from precursors, the kept ones count; read from
    # the engine's own protein quantity, only the protein groups stay. What the ENGINE did with
    # those peptides is the descriptor's sentence too (kept_contaminant_quant), per engine.
    keratins_in_database = c(if (isTRUE(k$requantified)) c(
      sprintf(paste0("%s%s mapping only to keratin-family contaminant entries (%s ",
                     "such entries in the database) %s KEPT"), pad, kept,
              fmt(length(k$exempt_accessions)), were),
      sprintf(paste0("%sand quantified as sample proteins; every other contaminant entry ",
                     "(trypsin, BSA ...) was handled as above."), pad))
    else c(
      sprintf(paste0("%s%s mapping only to keratin-family contaminant entries (%s ",
                     "such entries in the database) %s KEPT,"), pad, kept,
              fmt(length(k$exempt_accessions)), were),
      sprintf(paste0("%sso their protein groups stay in the analysis; every other contaminant ",
                     "entry (trypsin, BSA ...) was handled as above."), pad)),
      # what the engine did with the kept ones -- moot when none was kept (sage-review N5)
      if (!is.null(k$kept_quant) && isTRUE(k$n_precursors_kept > 0))
        sprintf("%s%s", pad, k$kept_quant)),
    no_keratin_contaminants = sprintf(paste0(
      "%sThe search database holds no keratin-family contaminant entry, so keratin %s ",
      "were quantified and tested as sample proteins."), pad, unit$plural),
    NULL)
  c(head, body, if (isTRUE(k$caution)) sprintf("%sCAUTION: %s", pad, k$note))
}
