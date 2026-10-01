#!/usr/bin/env Rscript
# =============================================================================
# run_sets.R  --  protein-SET tests on run_de.R's model, for every contrast:
#   camera (competitive) and fry (self-contained) side by side; the same tests with log2 run
#   depth in the model (sets that do not survive are marked); the fraction of each set that
#   DPC presence calls make up, with a measured-only version of both tests; and, for
#   pulldowns (--ip-map), between-condition tests after each IP is offset by its bait's
#   interactome. set_tests.R has the methods and why; references/set-tests.md how to read them.
#
# Run after run_de.R, on its outdir (it reads set_test_inputs.rds and the DE tables):
#   Rscript run_sets.R --de-dir output/tables
#       [--sets go|reactome|go,reactome|none]   a-priori collections (default go)
#       [--gmt sets.gmt --gmt-origin "where the file came from and when" [--gmt-kind public|own]]
#       [--organism mouse]                       default: the search FASTA's taxid
#       [--ontologies BP,CC,MF] [--min-size 10] [--max-size 500]
#       [--contrasts A-B,C-D]                    default: every run_de contrast
#       [--camera-cor 0.01]                      default: estimated for each set
#       [--fdr 0.05]                             default: run_de's --adjp
#       [--ip-map ip_map.csv] [--ip-normaliser interactome|bait] [--ip-min-enrichment 2]
#       [--outdir <de-dir>] [--figdir <de-dir>/../figures]
#
# --ip-map: one row per group -- Group, Role (bait | control), Bait (the bait's name; empty
#   for controls), Condition (e.g. Old / Young), optional Bait_gene (the bait protein's gene
#   symbol, e.g. Ryr2: enables the bait-protein normaliser as a cross-check).
# --ip-min-enrichment: the interactome's minimum pooled log2 enrichment over the control (2 = 4-fold).
# --gmt-kind: public (a published collection, named with its version in --gmt-origin, e.g.
#   "MSigDB Hallmark v2024.1") or own (a list you drew up: valid only if drawn up before the
#   results, which the record checks against the first DE run). Default own.
#
# Writes to --outdir: Sets_<contrast>.csv, Sets_baitnorm_<contrast>.csv (pulldowns),
#   sets_provenance.json, sets_methods.txt; to --figdir: sets_<contrast>.png.
# =============================================================================

getarg <- function(flag, default = NULL) {
  a <- commandArgs(trailingOnly = TRUE); i <- match(flag, a)
  if (is.na(i)) return(default)
  if (i == length(a) || startsWith(a[i + 1], "--")) return(TRUE)
  a[i + 1]
}
.script_dir <- local({
  f <- grep("^--file=", commandArgs(), value = TRUE)[1]
  if (is.na(f)) getwd() else dirname(normalizePath(sub("^--file=", "", f), mustWork = FALSE))
})
.sibling <- function(nm) {
  for (p in c(file.path(.script_dir, nm), nm)) if (file.exists(p)) return(p)
  stop(nm, " not found next to run_sets.R")
}
source(.sibling("set_tests.R"))
source(.sibling("skill_version.R"), encoding = "UTF-8")
suppressMessages(library(limma))
check_limma_internals()
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("run_sets.R needs jsonlite")
say <- function(...) message("[run_sets] ", sprintf(...))

# ---- arguments; every value that reaches the record carries where it came from ----------
de_dir <- getarg("--de-dir")
if (is.null(de_dir) || isTRUE(de_dir)) stop("Required: --de-dir <run_de.R outdir>")
inp_path <- file.path(de_dir, "set_test_inputs.rds")
if (!file.exists(inp_path))
  stop("run_sets.R: ", inp_path, " not found. run_de.R writes it (the skill version that has run_sets.R); re-run run_de.R ",
       "on this analysis first -- the set tests must use the model run_de fitted.")
inp <- readRDS(inp_path)
dprov_path <- file.path(de_dir, "de_provenance.json")
dprov <- if (file.exists(dprov_path)) jsonlite::fromJSON(dprov_path, simplifyVector = FALSE) else NULL

DEFAULT_TAG <- "DEFAULT -- not user-confirmed"
setting <- function(value, given, source = if (given) "user" else DEFAULT_TAG)
  list(value = value, source = source)
arg_or <- function(flag, default) { v <- getarg(flag); list(value = if (is.null(v)) default else v, given = !is.null(v)) }

outdir <- getarg("--outdir", de_dir)
figdir <- getarg("--figdir", NULL)
if (is.null(figdir)) {
  fd <- file.path(dirname(normalizePath(de_dir)), "figures")
  figdir <- if (dir.exists(fd)) fd else outdir
}
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(figdir, showWarnings = FALSE, recursive = TRUE)
a_sets <- arg_or("--sets", "go")
sources <- setdiff(tolower(trimws(strsplit(a_sets$value, ",")[[1]])), "none")
if (length(setdiff(sources, c("go", "reactome"))))
  stop("--sets takes go, reactome, go,reactome or none")
gmt_path <- getarg("--gmt")
gmt_origin <- getarg("--gmt-origin")
if (!is.null(gmt_path) && (is.null(gmt_origin) || isTRUE(gmt_origin)))
  stop("--gmt needs --gmt-origin \"<where the sets came from, and when>\": only sets fixed before ",
       "the results can be tested, and the record has to say what this file is.")
a_gk <- arg_or("--gmt-kind", "own")
if (!a_gk$value %in% c("public", "own")) stop("--gmt-kind must be public or own")
if (!length(sources) && is.null(gmt_path)) stop("no sets: give --sets go/reactome and/or --gmt")
a_min <- arg_or("--min-size", SET_DEFAULTS$min_size); min_size <- as.integer(a_min$value)
a_max <- arg_or("--max-size", SET_DEFAULTS$max_size); max_size <- as.integer(a_max$value)
if (is.na(min_size) || is.na(max_size) || min_size < 2 || max_size < min_size)
  stop("--min-size must be >= 2 and --max-size >= --min-size")
a_ont <- arg_or("--ontologies", paste(SET_DEFAULTS$ontologies, collapse = ","))
ontologies <- toupper(trimws(strsplit(a_ont$value, ",")[[1]]))
if (length(setdiff(ontologies, c("BP", "CC", "MF")))) stop("--ontologies takes BP, CC, MF")
a_cor <- arg_or("--camera-cor", NA)
camera_cor <- if (a_cor$given) as.numeric(a_cor$value) else NA_real_
if (a_cor$given && (is.na(camera_cor) || abs(camera_cor) >= 1)) stop("--camera-cor must be in (-1, 1)")
a_fdr <- getarg("--fdr")
fdr <- if (is.null(a_fdr)) inp$adjp else as.numeric(a_fdr)
fdr_source <- if (is.null(a_fdr)) "run_de --adjp" else "user"
a_con <- getarg("--contrasts")
tested <- if (is.null(a_con)) inp$contrasts else trimws(strsplit(a_con, ",")[[1]])
if (length(setdiff(tested, inp$contrasts)))
  stop("--contrasts not fitted by run_de: ", paste(setdiff(tested, inp$contrasts), collapse = ", "),
       ". Add them to run_de.R --contrasts (the set tests use run_de's model).")
ip_map_path <- getarg("--ip-map")
a_ipn <- arg_or("--ip-normaliser", "interactome")
if (!a_ipn$value %in% c("interactome", "bait")) stop("--ip-normaliser must be interactome or bait")
a_ipe <- arg_or("--ip-min-enrichment", SET_DEFAULTS$ip_min_enrichment)
ip_min_enr <- as.numeric(a_ipe$value)
if (is.na(ip_min_enr) || ip_min_enr < 0) stop("--ip-min-enrichment must be a log2 value >= 0")
organism_arg <- getarg("--organism")

sha256_file <- function(p) {
  if (requireNamespace("digest", quietly = TRUE)) return(digest::digest(p, algo = "sha256", file = TRUE))
  if (nzchar(Sys.which("shasum"))) return(sub(" .*", "", system2("shasum", c("-a", "256", shQuote(p)), stdout = TRUE)[1]))
  if (nzchar(Sys.which("sha256sum"))) return(sub(" .*", "", system2("sha256sum", shQuote(p), stdout = TRUE)[1]))
  NA_character_
}

# ---- proteins, their genes, the tested universe ------------------------------------
E_all <- inp$E
prot <- rownames(E_all)
g <- inp$genes
gid <- if ("Protein.Group" %in% names(g)) as.character(g$Protein.Group) else rownames(g)
gene_field <- if ("Genes" %in% names(g)) as.character(g$Genes)[match(prot, gid)] else rep(NA_character_, length(prot))
gene_first <- vapply(split_gene_field(gene_field), function(s) if (length(s)) s[1] else NA_character_, "")
label <- ifelse(is.na(gene_first), prot, gene_first)
# fry and camera's residuals need a value in every run: complete rows (all of them for dpc,
# whose protein matrix is complete by construction) with a finite statistic in every contrast.
finite_t <- Reduce(`&`, lapply(inp$stat[tested], function(s) is.finite(s$t)))
uni <- which(rowSums(is.na(E_all)) == 0 & finite_t)
say("%d of %d proteins tested (%s)", length(uni), nrow(E_all),
    if (length(uni) < nrow(E_all)) "complete rows with a finite statistic" else "all")

# ---- set collections: a-priori only -------------------------------------------------
org <- NULL; route <- NULL; map_all <- NULL; collections <- list(); source_info <- list()
need_orgdb <- length(sources) > 0
gmt <- NULL
if (!is.null(gmt_path)) {
  gmt <- read_gmt(gmt_path)
  need_orgdb <- need_orgdb || identical(gmt$id_type, "Entrez Gene ID")
}
if (need_orgdb) {
  org <- resolve_organism(organism_arg, inp$fasta_meta)
  route <- annotation_route(org)
  map_all <- map_proteins_to_entrez(gene_field, route)
  say("%s: %d of %d protein groups with a gene symbol mapped (%.1f%%; %d via an unambiguous alias)",
      route$note, map_all$n_mapped, map_all$n_with_symbol, 100 * map_all$mapping_rate, map_all$n_via_alias)
  universe_entrez <- unique(unlist(map_all$ids))
  if ("go" %in% sources) {
    collections$go <- list(coll = go_sets(route$package, universe_entrez, ontologies), ids = map_all$ids)
    source_info$go <- c(list(name = paste0("Gene Ontology (", paste(ontologies, collapse = ", "),
                                           "), propagated to ancestor terms (GO2ALLEGS, as limma::goana)")),
                        orgdb_info(route$package), list(GO.db = pkg_version("GO.db")))
  }
  if ("reactome" %in% sources) {
    # Reactome's genes are the mapped Entrez IDs: human's when the symbols went to human
    rorg <- if (route$route == "native") org else as.list(SET_ORGANISMS[SET_ORGANISMS$code == "Hs", ])
    collections$reactome <- list(coll = reactome_sets(rorg, universe_entrez), ids = map_all$ids)
    source_info$reactome <- c(list(name = paste0("Reactome pathways (", rorg$latin, ")")), reactome_info())
  }
}
if (!is.null(gmt)) {
  ids <- if (identical(gmt$id_type, "Entrez Gene ID")) map_all$ids else
    lapply(split_gene_field(gene_field), toupper)
  collections$gmt <- list(coll = gmt, ids = ids)
  gmt_time <- file.info(gmt_path)$mtime
  # The FIRST time DE results existed for this input (run_de carries it forward on every re-run):
  # a re-run of run_de must not make a list drawn up after the results look older than them.
  first <- dprov$first_written
  de_time <- if (!is.null(first)) as.POSIXct(first, tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ") else
    file.info(if (file.exists(dprov_path)) dprov_path else inp_path)$mtime
  de_time_src <- if (!is.null(first)) paste("de_provenance.json first_written (the first DE run on this input in this DE",
                                             "folder; a run_de into a new folder starts it again)") else
    "the DE record's file time (this run_de does not record its first run; a re-run resets it)"
  public <- identical(a_gk$value, "public")
  source_info$gmt <- list(
    name = paste0(if (public) "public collection" else "user list", " GMT ", basename(gmt_path)),
    path = normalizePath(gmt_path), kind = setting(a_gk$value, a_gk$given),
    sha256 = sha256_file(gmt_path), origin = gmt_origin, identifiers = gmt$id_type,
    file_modified = format(gmt_time, tz = "UTC", usetz = TRUE),
    de_results_first = format(de_time, tz = "UTC", usetz = TRUE), de_results_first_source = de_time_src,
    predates_results = isTRUE(gmt_time < de_time),
    note = if (public) paste("a public collection, named with its version in its origin; its file date is",
                             "recorded, not a test of when the sets were chosen")
           else if (isTRUE(gmt_time < de_time)) "a user-made list whose file predates the first DE results"
           else paste("a user-made list modified AFTER the first DE results: if any set in it was",
                      "chosen with the results in view, its tests are not valid (selection)"))
  if (!public && !source_info$gmt$predates_results) warning("[run_sets] ", source_info$gmt$note, call. = FALSE)
}

# Every collection indexed on some protein rows, as one list of sets with their metadata.
index_all <- function(rows) {
  parts <- lapply(names(collections), function(nm) {
    x <- index_sets(collections[[nm]]$coll, collections[[nm]]$ids[rows], min_size, max_size)
    list(index = x$index, meta = data.frame(Set_ID = names(x$index), Set_Name = x$name,
                                            Source = x$source, stringsAsFactors = FALSE),
         counts = x$counts)
  })
  names(parts) <- names(collections)
  idx <- do.call(c, unname(lapply(parts, `[[`, "index")))
  meta <- do.call(rbind, lapply(parts, `[[`, "meta"))
  if (is.null(idx)) idx <- list()
  if (is.null(meta)) meta <- data.frame(Set_ID = character(0), Set_Name = character(0), Source = character(0))
  names(idx) <- make.unique(paste(meta$Source, meta$Set_ID, sep = ":"))
  list(index = idx, meta = meta, counts = lapply(parts, `[[`, "counts"))
}
sets_u <- index_all(uni)
for (nm in names(sets_u$counts)) {
  cc <- sets_u$counts[[nm]]
  say("%s: %d sets defined, %d with a tested protein, %d with %d-%d (tested); %d smaller, %d larger",
      nm, cc$defined, cc$with_members, cc$tested, min_size, max_size, cc$below_min, cc$above_max)
  source_info[[nm]]$counts <- cc
}
if (!length(sets_u$index)) stop("no set has ", min_size, "-", max_size, " tested proteins")

# ---- the model: run_de's, and a check that the refit route reproduces it ----------------
is_random <- function(mdl) identical(mdl, "blocked") && identical(inp$block_effect, "random")
block_for <- function(mdl, runs = seq_len(ncol(E_all))) if (is_random(mdl)) inp$block[runs]
models <- unique(unname(unlist(inp$contrast_model[tested])))
refit <- lapply(setNames(models, models), function(mdl)
  fit_like_run_de(inp$method, inp$y, inp$design, block_for(mdl)))
dt <- vapply(tested, function(cn) {
  cs <- contrast_stats(refit[[inp$contrast_model[[cn]]]]$fit, inp$cmat[, cn])
  max(abs(cs$t - inp$stat[[cn]]$t), na.rm = TRUE)
}, 0)
if (max(dt) > 1e-6)
  stop(sprintf(paste0("run_sets.R: refitting run_de's model gave different t-statistics (max |dt| = %.3g, ",
                      "%s): the depth analysis would not be the same model. Is this the environment ",
                      "run_de.R ran in (limpa/limma versions)?"), max(dt), names(which.max(dt))))
say("model check: refitting reproduces run_de's t for every contrast (max |dt| = %.2g)", max(dt))

primary_model <- function(cn, rows = uni) {
  mdl <- inp$contrast_model[[cn]]; w <- inp$weights[[mdl]]
  list(E = E_all[rows, , drop = FALSE], weights = if (!is.null(w)) w[rows, , drop = FALSE],
       design = inp$design, contrast = inp$cmat[, cn], block = block_for(mdl),
       correlation = if (is_random(mdl)) inp$correlation,
       z = z_from_t(inp$stat[[cn]]$t, inp$stat[[cn]]$df_total)[rows])
}

# ---- pulldowns (--ip-map): which contrast is which --------------------------------------------------------------
ip_rec <- NULL
if (!is.null(ip_map_path)) {
  ipm <- utils::read.csv(ip_map_path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("Group", "Role", "Bait", "Condition")
  if (length(setdiff(need, names(ipm)))) stop("--ip-map needs columns ", paste(need, collapse = ", "))
  ipm$Role <- tolower(trimws(ipm$Role)); ipm[is.na(ipm)] <- ""
  if (length(setdiff(ipm$Group, levels(inp$groups))))
    stop("--ip-map groups not in the analysis: ", paste(setdiff(ipm$Group, levels(inp$groups)), collapse = ", "))
  if (!all(ipm$Role %in% c("bait", "control"))) stop("--ip-map Role must be bait or control")
  role <- setNames(ipm$Role, ipm$Group); bait <- setNames(ipm$Bait, ipm$Group)
  cond <- setNames(ipm$Condition, ipm$Group)
  bgene <- if ("Bait_gene" %in% names(ipm)) setNames(ipm$Bait_gene, ipm$Group) else NULL
  sides <- function(cn) {
    v <- inp$cmat[, cn]; v <- v[v != 0]
    if (length(v) == 2 && all(sort(v) == c(-1, 1)) && all(names(v) %in% ipm$Group))
      list(pos = names(v)[v > 0], neg = names(v)[v < 0]) else NULL
  }
  classify <- function(cn) {
    s <- sides(cn); if (is.null(s)) return("other")
    r <- c(role[[s$pos]], role[[s$neg]])
    if (all(r == "control")) return("control_vs_control")
    if (r[1] == "bait" && r[2] == "control" && cond[[s$pos]] == cond[[s$neg]]) return("bait_vs_control")
    if (all(r == "bait") && bait[[s$pos]] == bait[[s$neg]] && cond[[s$pos]] != cond[[s$neg]])
      return("bait_between_conditions")
    "other"
  }
  roles <- vapply(inp$contrasts, classify, "")
  read_de <- function(cn) utils::read.csv(file.path(de_dir, inp$de_tables[[cn]]$file),
                                          stringsAsFactors = FALSE, check.names = FALSE)
  ip_rec <- list(map = normalizePath(ip_map_path), map_sha256 = sha256_file(ip_map_path),
                 contrast_roles = as.list(roles),
                 normaliser = setting(a_ipn$value, a_ipn$given),
                 interactome_rule = sprintf(paste0("proteins enriched over the control in the bait's ",
                   "bait-vs-control comparison pooled over conditions (the mean of its per-condition ",
                   "contrasts, from run_de's model): adj.P.Val < %.3g and log2 enrichment >= %g (%.0f-fold)"),
                   fdr, ip_min_enr, 2^ip_min_enr),
                 min_enrichment = setting(ip_min_enr, a_ipe$given),
                 direction_bias = ip_direction_bias(ip_min_enr),
                 weak_rule = sprintf(paste0("a set is 'mostly weak interactors' when over half its proteins are in ",
                   "the interactome's bottom %g%% of pooled enrichment over the control"), 100 * SET_DEFAULTS$ip_weak_quantile),
                 baits = list())
}
ip_role <- function(cn) if (is.null(ip_rec)) NA_character_ else ip_rec$contrast_roles[[cn]]

# depth: the same pipeline with log2 run depth in the design -- one slope per run type in a
# pulldown (bait / control IPs), one slope otherwise
run_role <- if (!is.null(ip_rec)) { r <- unname(role[as.character(inp$groups)]); r[is.na(r)] <- "other"; r }
dep <- run_depth(inp$det_n[, colnames(E_all), drop = FALSE], inp$detection$zero_means)
ad <- add_depth(inp$design, rep(0, ncol(inp$design)), dep$value, run_role)
dd <- ad$design
depth_ok <- qr(dd)$rank == ncol(dd) && nrow(dd) - ncol(dd) >= 1
depth_note <- if (!depth_ok) {
  "not run: with log2 run depth the design is not of full rank or has no residual df"
} else if (length(ad$columns) > 1) {
  sprintf(paste0("log2 run depth added to the design with one slope per run type (%s), each centred ",
                 "within its runs; everything else unchanged"), paste(sort(unique(run_role)), collapse = ", "))
} else "log2 run depth added to the design; everything else unchanged"
depth_fit <- if (depth_ok) lapply(setNames(models, models), function(mdl)
  fit_like_run_de(inp$method, inp$y, dd, block_for(mdl))) else NULL
depth_model <- function(cn, rows = uni) {
  mdl <- inp$contrast_model[[cn]]; f <- depth_fit[[mdl]]
  con <- c(inp$cmat[, cn], setNames(rep(0, length(ad$columns)), ad$columns))
  cs <- contrast_stats(f$fit, con)
  list(E = f$E[rows, , drop = FALSE], weights = if (!is.null(f$weights)) f$weights[rows, , drop = FALSE],
       design = dd, contrast = con, block = block_for(mdl), correlation = f$correlation,
       z = z_from_t(cs$t, cs$df_total)[rows], logfc = cs$logFC[rows])
}


# ---- one contrast ------------------------------------------------------------------
call_of <- function(cf, ff) set_call(cf, ff, fdr)
holds <- function(fdr0, dir0, fdr1, dir1) set_holds(fdr0, dir0, fdr1, dir1, fdr)
members_of <- function(index, rows) vapply(index, function(i)
  paste(sort(unique(label[rows[i]])), collapse = ";"), "")
set_effect <- function(index, logfc, rows) vapply(index, function(i) mean(logfc[rows[i]]), 0)

# A contrast between different subjects (block levels) that uses several samples per subject is
# reported from the blocked fit, with the DE record's CAUTION (blocking.R): anti-conservative for
# sets that share subject-level variation (review S5: 0.23 camera / 0.16 fry false positives).
between_block_caution <- function(cn) {
  st <- ((dprov$block %||% list())$contrast_structure %||% list())[[cn]]
  identical(st, "between") && is_random(inp$contrast_model[[cn]])
}

run_contrast <- function(cn) {
  s <- sets_u; idx <- s$index
  pm <- primary_model(cn)
  prim <- test_sets(pm, idx, camera_cor)
  tab <- cbind(s$meta, NGenes = prim$NGenes,
               Mean_logFC = set_effect(idx, inp$stat[[cn]]$logFC, uni),
               Call = call_of(prim$camera_FDR, prim$fry_FDR), prim[, -1], stringsAsFactors = FALSE)
  # fry's call on other reference bases: camera does not depend on the basis, fry does (S4)
  sig <- tab$fry_FDR < fdr
  tab$fry_Reference_Stable <- ""
  if (any(sig)) {
    rc <- fry_calls_on_references(pm, idx, tab$fry_Direction, fdr)
    tab$fry_Reference_Stable[sig] <- ifelse(rc$held[sig] == rc$of, "stable", sprintf("%d of %d", rc$held[sig], rc$of))
  }
  # Depth in the model separates a set's change from a run-depth difference between the groups.
  # In a bait-vs-control IP the depth difference IS the enrichment (the bait brings its
  # interactors, and their precursors, with it): adjusting for it would remove the effect
  # tested, so it is not run there.
  dcols <- c("depth_Mean_logFC", "depth_camera_Direction", "depth_camera_PValue", "depth_camera_FDR",
             "depth_fry_Direction", "depth_fry_PValue", "depth_fry_FDR")
  if (depth_ok && !identical(ip_role(cn), "bait_vs_control")) {
    dm <- depth_model(cn)
    d <- test_sets(dm, idx, camera_cor)
    names(d) <- paste0("depth_", names(d))
    d$depth_Mean_logFC <- set_effect(idx, dm$logfc, seq_along(uni))
    tab <- cbind(tab, d[, dcols])
    tab$camera_Depth <- holds(tab$camera_FDR, tab$camera_Direction, tab$depth_camera_FDR, tab$depth_camera_Direction)
    tab$fry_Depth <- holds(tab$fry_FDR, tab$fry_Direction, tab$depth_fry_FDR, tab$depth_fry_Direction)
  } else {
    for (m in dcols) tab[[m]] <- NA
    tab$camera_Depth <- tab$fry_Depth <- "not run"
  }
  # presence calls (DPC): run_de's Evidence for this contrast
  ev <- inp$evidence[[cn]][uni]
  tab$Presence_Call_Fraction <- round(presence_fraction(idx, ev), 3)
  tab$Mostly_Presence_Calls <- tab$Presence_Call_Fraction > SET_DEFAULTS$presence_majority
  meas <- which(!(ev %in% PRESENCE_CALL_LABEL))
  mcols <- c("measured_NGenes", "measured_camera_Direction", "measured_camera_PValue",
             "measured_camera_FDR", "measured_fry_Direction", "measured_fry_PValue", "measured_fry_FDR")
  for (m in mcols) tab[[m]] <- NA
  if (length(meas) == length(uni)) {
    tab[, mcols] <- tab[, c("NGenes", "camera_Direction", "camera_PValue", "camera_FDR",
                            "fry_Direction", "fry_PValue", "fry_FDR")]
  } else if (length(meas) >= min_size) {
    mi <- reindex(idx, meas)
    ok <- lengths(mi) >= min_size
    if (any(ok)) {
      mt <- test_sets(subset_model(pm, meas), mi[ok], camera_cor)
      tab[ok, mcols] <- mt[, c("NGenes", "camera_Direction", "camera_PValue", "camera_FDR",
                               "fry_Direction", "fry_PValue", "fry_FDR")]
    }
  }
  flags <- character(nrow(tab))
  add <- function(cond, txt) { cond[is.na(cond)] <- FALSE
    flags <<- ifelse(cond, ifelse(nzchar(flags), paste0(flags, "; ", txt), txt), flags) }
  with_depth <- function(test) sprintf("%s: not separable from run depth (with depth: FDR %s, effect %s)", test,
    signif(tab[[paste0("depth_", test, "_FDR")]], 2), round(tab$depth_Mean_logFC, 2))
  add(tab$camera_Depth == "lost", with_depth("camera"))
  add(tab$fry_Depth == "lost", with_depth("fry"))
  add(tab$fry_Reference_Stable != "" & tab$fry_Reference_Stable != "stable",
      sprintf("fry call depends on the residual basis (holds on %s references)", tab$fry_Reference_Stable))
  add(tab$Mostly_Presence_Calls, "mostly presence calls")
  meas_sig <- (!is.na(tab$measured_camera_FDR) & tab$measured_camera_FDR < fdr) |
              (!is.na(tab$measured_fry_FDR) & tab$measured_fry_FDR < fdr)
  add(nzchar(tab$Call) & is.na(tab$measured_NGenes), "too few measured proteins for the measured-only test")
  add(nzchar(tab$Call) & !is.na(tab$measured_NGenes) & !meas_sig, "not significant on measured proteins only")
  if (between_block_caution(cn))
    add(nzchar(tab$Call), sprintf("between-%s contrast from the blocked fit: anti-conservative for sets sharing %s-level variation",
                                  inp$block_column, inp$block_column))
  tab$Flags <- flags
  tab <- tab[order(pmin(tab$camera_FDR, tab$fry_FDR), tab$camera_PValue), ]
  rownames(tab) <- NULL
  attr(tab, "df_residual") <- attr(prim, "df_residual")
  attr(tab, "global_correlation") <- attr(prim, "global_correlation")
  tab
}

summarise_tab <- function(tab) list(
  n_sets = nrow(tab),
  camera = sum(tab$camera_FDR < fdr), fry = sum(tab$fry_FDR < fdr),
  fry_up = sum(tab$fry_FDR < fdr & tab$fry_Direction == "Up"),
  fry_down = sum(tab$fry_FDR < fdr & tab$fry_Direction == "Down"),
  fry_basis_dependent = sum(tab$fry_Reference_Stable != "" & tab$fry_Reference_Stable != "stable"),
  both = sum(tab$Call == "camera + fry"),
  both_depth_robust = sum(tab$Call == "camera + fry" & tab$camera_Depth == "holds" & tab$fry_Depth == "holds"),
  camera_depth_holds = sum(tab$camera_Depth == "holds"), fry_depth_holds = sum(tab$fry_Depth == "holds"),
  depth_lost = sum(grepl("not separable from run depth", tab$Flags)),
  mostly_presence_calls_significant = sum(nzchar(tab$Call) & tab$Mostly_Presence_Calls))

# ---- the figure: set effect against each test's FDR, camera and fry side by side ----------
# One panel per test, so a comparison where camera finds nothing still shows what fry found (a
# broad shift puts many sets high in the fry panel and none in the camera one). Colour: which
# test is significant; hollow: flagged. Labels: each panel's most significant sets.
plot_sets <- function(tab, title, path, effect_lab, n_lab = 6) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return("skipped: ggplot2 not installed")
  call <- ifelse(nzchar(tab$Call), tab$Call, "neither")
  long <- do.call(rbind, lapply(c("camera", "fry"), function(test) data.frame(
    test = test, effect = tab$Mean_logFC, FDR = tab[[paste0(test, "_FDR")]],
    Call = factor(call, levels = c("camera + fry", "camera only", "fry only", "neither")),
    Flagged = ifelse(nzchar(tab$Flags), "flagged (see table)", "no flag"),
    NGenes = tab$NGenes, name = tab$Set_Name, stringsAsFactors = FALSE)))
  long$y <- -log10(pmax(long$FDR, 1e-300))
  long$test <- factor(long$test, levels = c("camera", "fry"),
                      labels = c("camera: more changed than the other proteins?",
                                 "fry: changed at all?"))
  lab <- do.call(rbind, lapply(split(long, long$test), function(d)
    utils::head(d[d$FDR < fdr, ][order(d$FDR[d$FDR < fdr]), ], n_lab)))
  if (!is.null(lab) && nrow(lab))
    lab$name <- ifelse(nchar(lab$name) > 40, paste0(substr(lab$name, 1, 38), "..."), lab$name)
  p <- ggplot2::ggplot(long, ggplot2::aes(effect, y, colour = Call, shape = Flagged)) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = 2, colour = "grey50") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey85") +
    ggplot2::geom_point(ggplot2::aes(size = NGenes), alpha = 0.7) +
    ggplot2::facet_wrap(~test, nrow = 1) +
    ggplot2::scale_colour_manual(values = c("camera + fry" = "#b2182b", "camera only" = "#ef8a62",
                                            "fry only" = "#2166ac", "neither" = "grey70"), drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("no flag" = 16, "flagged (see table)" = 1)) +
    ggplot2::scale_size_area(max_size = 4) +
    ggplot2::labs(x = effect_lab, y = "-log10 FDR", title = title,
                  subtitle = sprintf("%d sets; dashed line FDR = %.2g", nrow(tab), fdr),
                  colour = NULL, shape = NULL, size = "proteins") +
    ggplot2::theme_bw(base_size = 11) + ggplot2::theme(legend.position = "bottom", legend.box = "vertical")
  if (!is.null(lab) && nrow(lab)) {
    if (requireNamespace("ggrepel", quietly = TRUE))
      p <- p + ggrepel::geom_text_repel(data = lab, ggplot2::aes(label = name), size = 2.6,
                                        show.legend = FALSE, max.overlaps = 30, min.segment.length = 0)
    else p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(label = name), size = 2.6, vjust = -0.8,
                                     show.legend = FALSE)
  }
  ggplot2::ggsave(path, p, width = 11, height = 6, dpi = 150)
  basename(path)
}

# ---- run the contrasts ---------------------------------------------------------------
# p-values and FDRs to 4 significant digits, effects and correlations to 4 decimals: the
# tables are read, and 15 digits of noise tripled their size.
tidy_numbers <- function(tab) {
  for (col in grep("PValue$|FDR$", names(tab), value = TRUE))
    tab[[col]] <- signif(suppressWarnings(as.numeric(tab[[col]])), 4)
  for (col in intersect(c("Mean_logFC", "camera_Correlation"), names(tab))) tab[[col]] <- round(tab[[col]], 4)
  tab
}
# Every set's tested proteins, once: the tested proteins are the same in every comparison
# (a pulldown's interactome tables carry their own Members column).
members_file <- "sets_members.csv"
utils::write.csv(data.frame(sets_u$meta, NGenes = lengths(sets_u$index),
                            Members = members_of(sets_u$index, uni), check.names = FALSE),
                 file.path(outdir, members_file), row.names = FALSE)
# How many independent units each side of a contrast has: subjects (block levels) for a
# between-block contrast, runs otherwise -- the "n" a reader needs next to "none found".
side_units <- function(cn) {
  v <- inp$cmat[, cn]; v <- v[names(v) %in% levels(inp$groups) & v != 0]
  runs <- lapply(list(pos = names(v)[v > 0], neg = names(v)[v < 0]),
                 function(gs) which(as.character(inp$groups) %in% gs))
  between <- !is.null(inp$block) &&
    identical(((dprov$block %||% list())$contrast_structure %||% list())[[cn]], "between")
  count <- function(r) if (between) length(unique(inp$block[r])) else length(r)
  per <- if (between) max(length(runs$pos) / max(1, count(runs$pos)), length(runs$neg) / max(1, count(runs$neg))) else 1
  list(pos = count(runs$pos), neg = count(runs$neg),
       words = if (!between) "runs" else if (per == 1) sprintf("samples, one per %s", inp$block_column)
               else sprintf("%s levels, %g samples each averaged", inp$block_column, per),
       pos_groups = names(v)[v > 0])
}

# Each comparison's reading: contrast_reading() (set_tests.R), written once the pulldown's recovery
# estimates are known (below).
contrast_rec <- list()
for (cn in tested) {
  tab <- run_contrast(cn)
  f <- file.path(outdir, sprintf("Sets_%s.csv", make.names(cn)))
  utils::write.csv(tidy_numbers(tab), f, row.names = FALSE)
  fig <- plot_sets(tab, sprintf("Protein-set tests: %s", cn),
                   file.path(figdir, sprintf("sets_%s.png", make.names(cn))),
                   "Set effect: mean log2 fold change of the set's proteins")
  sm <- summarise_tab(tab)
  # the contrast applied to the groups' mean log2 depth: how far apart the compared runs are in depth
  gm <- tapply(dep$value, as.character(inp$groups), mean)
  cw <- inp$cmat[intersect(rownames(inp$cmat), names(gm)), cn]
  su <- side_units(cn)
  rho0 <- attr(tab, "global_correlation")
  contrast_rec[[cn]] <- c(list(file = basename(f), figure = fig, model = inp$contrast_model[[cn]],
                               ip_role = ip_role(cn),
                               n_significant_proteins = inp$de_tables[[cn]]$n_significant,
                               n_pos = su$pos, n_neg = su$neg, n_words = su$words, pos_groups = as.list(su$pos_groups),
                               df_residual = attr(tab, "df_residual"),
                               camera_global_correlation = round(rho0, 3),
                               camera_vif_50 = round(max(1, 1 + 49 * rho0), 1),
                               camera_low_power = isTRUE(rho0 > SET_DEFAULTS$camera_low_power),
                               between_block_caution = between_block_caution(cn),
                               depth_difference_log2 = round(sum(cw * gm[names(cw)]), 3),
                               depth_tested = depth_ok && !identical(ip_role(cn), "bait_vs_control"),
                               # the proteins themselves: a "broad shift" must be one (review W1)
                               protein_share_up = round(mean(inp$stat[[cn]]$logFC[uni] > 0), 3),
                               protein_share_down = round(mean(inp$stat[[cn]]$logFC[uni] < 0), 3),
                               protein_median_logFC = round(stats::median(inp$stat[[cn]]$logFC[uni]), 4)), sm)
  say("%-24s %d sets, FDR < %.2g: camera %d, fry %d, both %d; %s -> %s", cn, sm$n_sets, fdr,
      sm$camera, sm$fry, sm$both,
      if (isTRUE(contrast_rec[[cn]]$depth_tested))
        sprintf("with run depth in the model camera %d and fry %d hold", sm$camera_depth_holds, sm$fry_depth_holds)
      else "run depth not tested", basename(f))
}

# ---- pulldowns (--ip-map): bait-normalised tests per bait ------------------------------
if (!is.null(ip_rec)) {
  for (b in unique(ipm$Bait[ipm$Role == "bait"])) {
    gb <- ipm$Group[ipm$Role == "bait" & ipm$Bait == b]
    rec <- list(groups = as.list(gb))
    vs_ctrl <- names(roles)[roles == "bait_vs_control"]
    vs_ctrl <- vs_ctrl[vapply(vs_ctrl, function(cn) sides(cn)$pos %in% gb, TRUE)]
    between <- names(roles)[roles == "bait_between_conditions"]
    between <- between[vapply(between, function(cn) sides(cn)$pos %in% gb, TRUE)]
    have_cond <- unique(cond[vapply(vs_ctrl, function(cn) sides(cn)$pos, "")])
    if (!length(between)) { rec$skipped <- "run_de has no between-condition contrast for this bait"; ip_rec$baits[[b]] <- rec; next }
    if (length(setdiff(unique(cond[gb]), have_cond))) {
      rec$skipped <- paste("no bait-vs-control contrast in run_de for condition(s)",
                           paste(setdiff(unique(cond[gb]), have_cond), collapse = ", "))
      ip_rec$baits[[b]] <- rec; next
    }
    # the interactome: enrichment over the control POOLED over conditions (set_tests.R says why)
    mdl_e <- unique(unname(unlist(inp$contrast_model[vs_ctrl])))
    if (length(mdl_e) != 1) {
      rec$skipped <- "the bait-vs-control contrasts come from different fits; the pooled enrichment is undefined"
      ip_rec$baits[[b]] <- rec; next
    }
    fe <- refit[[mdl_e]] %||% fit_like_run_de(inp$method, inp$y, inp$design, block_for(mdl_e))
    ce <- contrast_stats(fe$fit, pooled_enrichment_contrast(inp$cmat, vs_ctrl))
    adj_e <- stats::p.adjust(2 * stats::pt(-abs(ce$t), ce$df_total), "BH")
    enrich <- ce$logFC
    runs <- which(as.character(inp$groups) %in% gb)
    rows_b <- which(adj_e < fdr & enrich >= ip_min_enr)
    rows_b <- rows_b[rowSums(is.na(E_all[rows_b, runs, drop = FALSE])) == 0]
    rec$interactome_from <- as.list(vs_ctrl)
    rec$interactome_contrast <- paste0("mean of ", paste(vs_ctrl, collapse = " and "))
    rec$n_interactome <- length(rows_b)
    weak_row <- logical(length(enrich)); weak_row[rows_b] <- weak_interactors(enrich[rows_b])
    rec$n_pooled_enriched_any <- sum(adj_e < fdr & enrich > 0)
    # context: what "enriched at every condition" would have taken, one condition at a time
    rec$n_enriched_per_condition <- as.list(setNames(vapply(vs_ctrl, function(cn) {
      d <- read_de(cn); sum(d$logFC > 0 & d$adj.P.Val < fdr) }, 0L), vs_ctrl))
    say("bait %s: interactome %d proteins (pooled enrichment over the control >= %g log2 at FDR %.2g; %d at any enrichment)",
        b, length(rows_b), ip_min_enr, fdr, rec$n_pooled_enriched_any)
    design_b <- inp$design[runs, , drop = FALSE]
    design_b <- design_b[, colSums(abs(design_b)) > 0, drop = FALSE]
    if (qr(design_b)$rank < ncol(design_b) || nrow(design_b) - ncol(design_b) < 1) {
      rec$skipped <- "the bait's runs give no estimable design with residual df"; ip_rec$baits[[b]] <- rec; next
    }
    # the bait protein: its gene symbol among the interactome rows, if the map names it
    bg <- if (!is.null(bgene)) unique(bgene[gb][nzchar(bgene[gb])]) else character(0)
    bait_row <- if (length(bg) == 1) rows_b[vapply(split_gene_field(gene_field[rows_b]),
                                                   function(s) toupper(bg) %in% toupper(s), TRUE)] else integer(0)
    rec$bait_gene <- if (length(bg)) bg else NA_character_
    rec$bait_protein <- if (length(bait_row) == 1) prot[bait_row] else NA_character_
    hows <- a_ipn$value
    alt <- setdiff(c("interactome", "bait"), hows)
    if (alt == "interactome" || length(bait_row) == 1) hows <- c(hows, alt)
    else rec$bait_normaliser_note <- if (!length(bg)) "no Bait_gene in --ip-map: the bait-protein cross-check was not run"
      else sprintf("the bait gene %s is not (uniquely) in the interactome: the bait-protein cross-check was not run", bg)
    if (a_ipn$value == "bait" && length(bait_row) != 1)
      stop("--ip-normaliser bait: bait ", b, ": ", rec$bait_normaliser_note)
    res <- setNames(lapply(between, function(cn) list()), between); offs <- list()
    for (how in hows) {
      off <- ip_offsets(E_all, runs, rows_b, how, bait_row)
      uni_b <- if (how == "bait") setdiff(rows_b, bait_row) else rows_b
      offs[[how]] <- off
      if (how == hows[1] && ip_min_enr > 0) {   # how much a fall in one condition would cost (review W2)
        rec$losses <- ip_loss_sensitivity(enrich[uni_b], index_all(uni_b)$index, ip_min_enr)
        rec$losses_reading <- ip_losses_reading(rec$losses, ip_min_enr)
      }
      yb <- if (inp$method == "dpc") inp$y[uni_b, runs] else E_all[uni_b, runs, drop = FALSE]
      if (inp$method == "dpc") yb$E <- sweep(yb$E, 2, off) else yb <- sweep(yb, 2, off)
      fits <- list()
      for (cn in between) {
        mdl <- inp$contrast_model[[cn]]
        if (is.null(fits[[mdl]]))
          fits[[mdl]] <- tryCatch(fit_like_run_de(inp$method, yb, design_b, block_for(mdl, runs)),
                                  error = function(e) e)
        f <- fits[[mdl]]
        if (inherits(f, "error")) {
          res[[cn]][[how]] <- list(error = conditionMessage(f)); next
        }
        con <- inp$cmat[colnames(design_b), cn]
        cs <- contrast_stats(f$fit, con)
        s <- index_all(uni_b)
        if (!length(s$index)) { res[[cn]][[how]] <- list(error = sprintf("no set has %d-%d interactome proteins", min_size, max_size)); next }
        m <- list(E = f$E, weights = f$weights, design = design_b, contrast = con,
                  block = block_for(mdl, runs), correlation = f$correlation,
                  z = z_from_t(cs$t, cs$df_total))
        tt0 <- test_sets(m, s$index, camera_cor)
        tt <- cbind(s$meta, NGenes = tt0$NGenes, Mean_logFC = set_effect(s$index, cs$logFC, seq_along(uni_b)),
                    Call = call_of(tt0$camera_FDR, tt0$fry_FDR), tt0[, -1], stringsAsFactors = FALSE)
        tt$Weak_Fraction <- round(vapply(s$index, function(i) mean(weak_row[uni_b[i]]), 0), 3)
        tt$Presence_Call_Fraction <- round(presence_fraction(s$index, inp$evidence[[cn]][uni_b]), 3)
        tt$Members <- members_of(s$index, uni_b)
        res[[cn]][[how]] <- list(tab = tt, correlation = f$correlation, df = attr(tt0, "df_residual"),
                                 global = attr(tt0, "global_correlation"),
                                 n_rows = length(uni_b), n_runs = c(sum(inp$groups[runs] == sides(cn)$pos),
                                                                   sum(inp$groups[runs] == sides(cn)$neg)))
      }
    }
    # the recovery difference each reference implies: mean offset, first group minus second
    rec$offset_difference <- lapply(offs, function(o) lapply(setNames(between, between), function(cn) {
      s2 <- sides(cn)
      round(mean(o[inp$groups[runs] == s2$pos]) - mean(o[inp$groups[runs] == s2$neg]), 3)
    }))
    rec$contrasts <- list()
    for (cn in between) {
      r <- res[[cn]]
      prim <- r[[hows[1]]]
      if (!is.null(prim$error)) { rec$contrasts[[cn]] <- list(error = prim$error); next }
      tab <- prim$tab
      if (length(hows) > 1 && is.null(r[[hows[2]]]$error)) {
        other <- r[[hows[2]]]$tab
        k <- match(paste(tab$Source, tab$Set_ID), paste(other$Source, other$Set_ID))
        for (col in c("camera_Direction", "camera_FDR", "fry_Direction", "fry_FDR"))
          tab[[paste0(hows[2], "_", col)]] <- other[[col]][k]
        tab$Reference_Robust <- ifelse(!nzchar(tab$Call), "",
          ifelse((tab$camera_FDR >= fdr | (!is.na(tab[[paste0(hows[2], "_camera_FDR")]]) &
                   tab[[paste0(hows[2], "_camera_FDR")]] < fdr &
                   tab[[paste0(hows[2], "_camera_Direction")]] == tab$camera_Direction)) &
                 (tab$fry_FDR >= fdr | (!is.na(tab[[paste0(hows[2], "_fry_FDR")]]) &
                   tab[[paste0(hows[2], "_fry_FDR")]] < fdr &
                   tab[[paste0(hows[2], "_fry_Direction")]] == tab$fry_Direction)),
                 "holds with the other reference", "depends on the reference"))
      }
      rr <- if (is.null(tab$Reference_Robust)) rep("", nrow(tab)) else tab$Reference_Robust
      fl <- character(nrow(tab))
      addf <- function(cond, txt) { cond[is.na(cond)] <- FALSE
        fl <<- ifelse(cond, ifelse(nzchar(fl), paste0(fl, "; ", txt), txt), fl) }
      addf(rr == "depends on the reference", "depends on the reference")
      addf(tab$Weak_Fraction > 0.5, "mostly weak interactors (the interactome's least enriched quarter)")
      addf(tab$Presence_Call_Fraction > SET_DEFAULTS$presence_majority, "mostly presence calls")
      tab$Flags <- fl
      tab <- tab[order(pmin(tab$camera_FDR, tab$fry_FDR), tab$camera_PValue), ]
      f <- file.path(outdir, sprintf("Sets_baitnorm_%s.csv", make.names(cn)))
      utils::write.csv(tidy_numbers(tab), f, row.names = FALSE)
      fig <- plot_sets(tab, sprintf("Bait-normalised set tests: %s (%s reference)", cn, hows[1]),
                       file.path(figdir, sprintf("sets_baitnorm_%s.png", make.names(cn))),
                       "Set effect: mean log2 fold change relative to the reference")
      rc <- list(file = basename(f), figure = fig, model = inp$contrast_model[[cn]],
                 reference = hows[1], cross_check = if (length(hows) > 1) hows[2] else NA,
                 n_sets = nrow(tab), camera = sum(tab$camera_FDR < fdr),
                 fry = sum(tab$fry_FDR < fdr), both = sum(tab$Call == "camera + fry"),
                 depends_on_reference = if (!is.null(tab$Reference_Robust))
                   sum(tab$Reference_Robust == "depends on the reference") else NA,
                 mostly_weak_significant = sum(nzchar(tab$Call) & tab$Weak_Fraction > 0.5),
                 unflagged_significant = sum(nzchar(tab$Call) & !nzchar(tab$Flags)),
                 n_pos = prim$n_runs[1], n_neg = prim$n_runs[2], df_residual = prim$df,
                 n_prior_proteins = prim$n_rows, camera_global_correlation = round(prim$global, 3),
                 camera_low_power = isTRUE(prim$global > SET_DEFAULTS$camera_low_power))
      sg <- nzchar(tab$Call)
      if (length(hows) > 1 && !is.null(rec$offset_difference[[hows[2]]]))
        rc$recovery_difference_log2 <- round(abs(rec$offset_difference[[hows[1]]][[cn]] -
                                                 rec$offset_difference[[hows[2]]][[cn]]), 3)
      rc$reading <- bait_reading(utils::modifyList(rc, list(
        bait = b, pos_group = sides(cn)$pos, n_interactome = rec$n_interactome,
        cross_check = if (!is.null(tab$Reference_Robust)) hows[2] else NA,
        recovery_difference = rc$recovery_difference_log2 %||% NA_real_,
        losses = rec$losses_reading %||% NA_character_,
        sig = data.frame(direction = ifelse(tab$fry_FDR < fdr, tab$fry_Direction, tab$camera_Direction)[sg],
                         effect = tab$Mean_logFC[sg],
                         depends = if (is.null(tab$Reference_Robust)) rep(FALSE, sum(sg))
                                   else tab$Reference_Robust[sg] == "depends on the reference",
                         weak = tab$Weak_Fraction[sg] > 0.5, flagged = nzchar(tab$Flags[sg])))))
      rec$contrasts[[cn]] <- rc
      say("bait-normalised %-18s %d interactome sets: camera %d, fry %d, both %d -> %s", cn, nrow(tab),
          rec$contrasts[[cn]]$camera, rec$contrasts[[cn]]$fry, rec$contrasts[[cn]]$both, basename(f))
    }
    ip_rec$baits[[b]] <- rec
  }
}
for (cn in tested) {
  if (identical(ip_role(cn), "bait_between_conditions")) {
    ref <- a_ipn$value
    for (b in names(ip_rec$baits)) {
      d <- ((ip_rec$baits[[b]]$offset_difference %||% list())[[ref]] %||% list())[[cn]]
      if (!is.null(d)) contrast_rec[[cn]]$recovery <- list(bait = b, reference = ref, difference_log2 = d)
    }
  }
  contrast_rec[[cn]]$reading <- contrast_reading(contrast_rec[[cn]], fdr)
}

# ---- provenance + methods ------------------------------------------------------------------
blocked_any <- any(vapply(tested, function(cn) is_random(inp$contrast_model[[cn]]), TRUE))
suggested <- intersect(unlist((dprov$set_tests %||% list())$suggested_for), tested)
trigger_n <- (dprov$set_tests %||% list())$suggest_below %||% SET_DEFAULTS$trigger
low_power <- vapply(contrast_rec, function(r) isTRUE(r$camera_low_power), TRUE)
# rule 2: a default that reaches the Methods says so where it stands
tag <- function(given) if (given) "" else " (DEFAULT \u2014 not user-confirmed)"
fixed_any <- identical(inp$block_effect, "fixed")
src_words <- paste(vapply(source_info, function(x)
  if (!is.null(x$origin)) sprintf("%s (%s; sha256 %s; %s%s%s)", x$name, x$origin, substr(x$sha256, 1, 12),
                                  if (identical(x$kind$value, "own")) "a user-made list" else "a public collection",
                                  tag(x$kind$source != DEFAULT_TAG),
                                  if (!identical(x$kind$value, "own")) "" else if (isTRUE(x$predates_results))
                                    ", its file dated before the first DE results" else ", its file dated AFTER the first DE results")
  else if (!is.null(x$GOSOURCEDATE)) sprintf("%s, GO release %s (%s %s)", x$name, x$GOSOURCEDATE, x$package, x$version)
  else if (!is.null(x$SOURCEDATE)) sprintf("%s (%s %s)", x$name, x$package, x$version)
  else x$name, ""), collapse = "; ")
# the pulldown's losses sentence: each bait's own numbers (ip_loss_sensitivity), as a range over baits
losses_methods <- function() {
  L <- Filter(Negate(is.null), lapply(ip_rec$baits, function(r) r$losses))
  if (is.na(ip_rec$direction_bias) || !length(L)) return("")
  span <- function(k, j) { v <- 100 * vapply(L, function(l) l[[k]][j], 0)
    if (round(min(v)) == round(max(v))) sprintf("%.0f%%", min(v)) else sprintf("%.0f-%.0f%%", min(v), max(v)) }
  sprintf(paste0("Because the minimum applies to the enrichment pooled over conditions, a protein that falls in one ",
                 "condition can drop out before testing, so losses relative to the complex are harder to detect than ",
"gains: a %g log2 fall in one condition would have put about %s of %s below the minimum, and about %s ",
                 "of the tested sets would still have had enough members to be tested had all their proteins fallen ",
                 "that much (about %s for a %g log2 fall; estimated from the fitted enrichments, so slightly ",
                 "optimistic). "),
          L[[1]]$fall_log2[1], span("below_minimum", 1), if (length(L) > 1) "each bait's interactome" else "the interactome",
          span("sets_still_tested", 1), span("sets_still_tested", 2), L[[1]]$fall_log2[2])
}
methods_paragraph <- paste0(
  "Protein-set tests used the differential-expression model itself (run_sets.R, ",
  skill_label(skill_version(.script_dir)), "; limma ", pkg_version("limma"), "). For each contrast, ",
  "camera (competitive: are a set's proteins more changed than the others? Wu and Smyth 2012) ranked ",
  "that contrast's moderated t-statistics from the DE fit, converted to z-scores, with the ",
  "inter-protein correlation ", if (is.na(camera_cor)) "estimated for each set from the model's standardised residuals"
  else sprintf("preset to %g", camera_cor), tag(a_cor$given),
  ", and fry (self-contained: are they changed at all? Giner and Smyth 2016) used the same expression values",
  if (identical(inp$method, "dpc")) ", precision weights" else "", ", design and contrast",
  " (every protein's residual effects on one canonical basis, so results do not depend on the machine)",
  if (blocked_any) sprintf(paste0(" and, for contrasts fitted with %s as a random effect, the block and its ",
                                  "consensus correlation (%.3f); camera cannot take a block, so its statistics ",
                                  "came from the blocked fit and its correlation from the block-whitened residuals"),
                           inp$block_column, inp$correlation) else "",
  if (fixed_any) sprintf(" (%s as a fixed effect in the design)", inp$block_column) else "",
  ". Sets were ", src_words, ", with ", min_size, "-", max_size, " tested member proteins",
  tag(a_min$given && a_max$given),
  if (!is.null(map_all)) sprintf(" (%.1f%% of protein groups with a gene symbol mapped to %s Entrez Gene IDs%s)",
                                 100 * map_all$mapping_rate, if (route$route == "native") org$name else "human",
                                 if (route$route == "native") "" else ", by gene-symbol identity") else "",
  ". P-values were adjusted across the sets of each contrast (Benjamini-Hochberg), separately for each test; a ",
  "set's call joins the two lists, and no correction was made across the ", length(tested), " comparisons. ",
  if (length(suggested)) sprintf(paste0("Set tests were run because the differential-expression analysis found ",
                                        "fewer than %d significant proteins in %s. "), trigger_n, paste(suggested, collapse = ", ")) else "",
  "fry calls were re-tested on ", length(SET_DEFAULTS$ref_seeds_check), " other reference bases and those that did ",
  "not hold on every one are flagged. ",
  if (any(low_power)) sprintf(paste0("In %d of the %d comparisons random protein sets correlated across runs above ",
                                     "%g (whole runs moving together), which leaves camera little power; its counts ",
                                     "there are not interpreted. "), sum(low_power), length(low_power),
                              SET_DEFAULTS$camera_low_power) else "",
  if (depth_ok) paste0("As a sensitivity analysis both tests were repeated with ", sub(", centred$", "", dep$definition),
                       " added to the design",
                       if (length(ad$columns) > 1) " (one slope for each run type: bait and control IPs)" else "",
                       "; a set significant without it but not with it is marked as not separable from run depth",
                       if (!is.null(ip_rec)) " (not done for bait-vs-control contrasts, where the depth difference is the enrichment itself)" else "",
                       ". ")
  else "The run-depth sensitivity analysis could not be run (design not estimable with depth). ",
  if (identical(inp$method, "dpc")) paste0("Because DPC-Quant gives every protein a value in every run, the fraction ",
    "of each set's proteins never measured in one of the compared groups (presence calls) is reported, with both ",
    "tests repeated on measured proteins only. ") else "",
  if (!is.null(ip_rec)) sprintf(paste0("For the pulldown, each bait's IPs were offset by the median of the bait's ",
    "interactome -- proteins enriched over the control IP in the bait-vs-control comparison pooled over ",
    "conditions (adj.P.Val < %.3g, at least %.0f-fold%s), a selection independent of the between-condition ",
    "comparison -- and sets within that interactome were tested between conditions with the same pipeline ",
    "refitted on the bait's IPs%s. %sSets over half of whose proteins lay in the interactome's least-enriched ",
    "quarter were flagged as mostly weak interactors."), fdr, 2^ip_min_enr, tag(a_ipe$given),
    if (any(vapply(ip_rec$baits, function(r) !is.null(r$bait_protein) && !is.na(r$bait_protein), TRUE)))
      "; the bait protein's own value was used as a second reference, as a cross-check" else "",
    losses_methods()) else "")
writeLines(c("Protein-set tests -- methods", strrep("=", 40), "", strwrap(methods_paragraph, 95), "",
             "Citations: Wu D, Smyth GK (2012) Nucleic Acids Res 40(17):e133 (camera);",
             "           Giner G, Smyth GK (2016) F1000Research 5:2605 (fry);",
             "           Ritchie ME et al. (2015) Nucleic Acids Res 43(7):e47 (limma)."),
           file.path(outdir, "sets_methods.txt"))

prov <- list(
  tool = "run_sets.R", skill_version = skill_label(skill_version(.script_dir)),
  run_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  command = paste(c("Rscript run_sets.R", commandArgs(trailingOnly = TRUE)), collapse = " "),
  de_dir = normalizePath(de_dir),
  de_inputs = list(file = basename(inp_path), sha256 = sha256_file(inp_path), written_by = inp$written_by),
  method = inp$method,
  model = list(
    reuse = "each contrast's moderated t, expression, precision weights, design, contrasts and block from run_de.R (set_test_inputs.rds); no refit for the primary tests",
    refit_check = list(max_abs_t_difference = signif(max(dt), 3),
                       what = "fit_like_run_de() (used for the depth and bait-normalised refits) reproduces run_de's t for every tested contrast"),
    block = if (!is.null(inp$block)) list(column = inp$block_column, effect = inp$block_effect,
                                          correlation = inp$correlation) else NULL),
  tests = list(
    camera = list(what = "competitive: are the set's proteins more changed than the other tested proteins?",
                  statistic = "run_de's moderated t for the contrast, as z-scores (limma::zscoreT, Hill; df = the fit's df.total)",
                  inter_gene_correlation = if (is.na(camera_cor)) "estimated for each set from the model's standardised residual effects (camera's inter.gene.cor = NA mode; df = min(residual df, G - 2)); VIF >= 1"
                                           else sprintf("preset %g (camera's preset mode; df = G - 2)", camera_cor),
                  block = "camera cannot take a block (Smyth, support.bioconductor.org/p/78299): for a random-block contrast the statistics come from the blocked fit and the residual effects are whitened by the block correlation, as fry's are",
                  equivalence = "with no block and no weights this is limma::camera() exactly; with a preset correlation it is limma::cameraPR() on the z-scores (tests/test_set_tests.py)"),
    residual_basis = paste("every protein's residual effects on one canonical basis (set_tests.R lm_effects): the",
                           "polar factor of the design's residual projector applied to a fixed matrix",
                           sprintf("(seed %d), and with precision weights its polar projection into each", SET_DEFAULTS$ref_seed),
                           "protein's weighted residual space. limma's Householder basis is not unique: R's qr() of",
                           "the same design gave different reflection signs on different HIVE CPU types (PROT_0756:",
                           "p-values moved by up to 0.09 between otherwise identical runs). A change of reference",
                           "rotates every protein's basis by the same rotation, with weights too, and camera is",
                           "invariant to that: it does not depend on the reference. fry depends on it only through",
                           "its robust variance's largest-coordinate term (up to 0.03 in p without weights, 0.07",
                           "with). Apart from the reference: without weights camera is limma::camera() exactly and",
                           "fry is limma::fry() exactly when that term is taken from limma's basis; with weights",
                           "limma's per-protein bases are not aligned, and camera differs from weighted",
                           "limma::camera by up to 0.04 in p, fry from weighted limma::fry by up to 0.1 -- a",
                           "difference from limma, not a dependence on the reference. Each fry call is re-tested on",
                           sprintf("%d other references (seeds %s);", length(SET_DEFAULTS$ref_seeds_check),
                                   paste(SET_DEFAULTS$ref_seeds_check, collapse = ", ")),
                           "fry_Reference_Stable says whether it held on all of them."),
    fry = list(what = "self-contained: are the set's proteins changed at all?",
               call = "fry's computation (limma 3.68: .lmEffects for the contrast effect, standardize = posterior.sd, .fryEffects) on the fit's precision weights, design and contrast, and the block + correlation for a random-block contrast, on the canonical residual basis",
               dpc_note = if (identical(inp$method, "dpc")) "fry uses the DPC values of every run; limpa's reduced residual df for proteins wholly imputed in a group reaches camera's statistics but not fry: see the presence-call columns" else NULL),
    adjustment = "Benjamini-Hochberg across the sets of each contrast, separately for camera and fry",
    not_used = "geneSetTest (assumes independent proteins) and any list drawn up after seeing results"),
  settings = list(
    sets = setting(a_sets$value, a_sets$given),
    ontologies = if ("go" %in% sources) setting(as.list(ontologies), a_ont$given),
    min_size = setting(min_size, a_min$given), max_size = setting(max_size, a_max$given),
    camera_inter_gene_cor = setting(if (is.na(camera_cor)) "estimated per set" else camera_cor, a_cor$given),
    fdr = setting(fdr, !is.null(a_fdr), fdr_source),
    presence_majority = setting(SET_DEFAULTS$presence_majority, FALSE),
    camera_low_power_above = setting(SET_DEFAULTS$camera_low_power, FALSE),
    reference_seeds_checked = setting(as.list(SET_DEFAULTS$ref_seeds_check), FALSE),
    ip_min_enrichment = if (!is.null(ip_rec)) setting(ip_min_enr, a_ipe$given),
    organism = if (!is.null(org)) setting(org$name, !is.null(organism_arg), org$source) else NULL),
  annotation = if (!is.null(route)) list(route = route$route, note = route$note, package = orgdb_info(route$package),
    mapping = map_all[c("n_proteins", "n_with_symbol", "n_mapped", "mapping_rate", "n_symbols",
                        "n_symbols_mapped", "n_via_alias")],
    rule = "a protein group is in a set when any of its genes is") else NULL,
  sources = source_info,
  tested_proteins = list(n = length(uni), of = nrow(E_all),
                         rule = "complete rows with a finite statistic in every tested contrast"),
  depth = list(ran = depth_ok, definition = dep$definition, note = depth_note,
               correlation = if (depth_ok) lapply(depth_fit, `[[`, "correlation") else NULL,
               per_run = as.list(setNames(as.numeric(dep$counts), names(dep$counts))),
               columns = as.list(ad$columns),
               rule = "a set significant (FDR) in the primary test is 'lost' when the same test is not significant in the same direction with depth in the design; the depth_* columns give that test's p-value, FDR and the set's effect with depth"),
  presence_calls = list(label = PRESENCE_CALL_LABEL, source = "run_de's Evidence column for the contrast",
                        mostly = sprintf("more than %.0f%% of the set's tested proteins", 100 * SET_DEFAULTS$presence_majority),
                        measured_only = "both tests repeated with presence-call proteins removed from the sets and the background"),
  contrasts = contrast_rec,
  pulldown = ip_rec,
  # where the figures are, relative to this record's folder (the session moves as a whole)
  figure_dir = relative_path(figdir, outdir),
  columns = SET_COLUMNS,
  members = members_file,
  reading = SET_READING,
  methods_paragraph = methods_paragraph,
  packages = list(limma = pkg_version("limma"), limpa = pkg_version("limpa"),
                  AnnotationDbi = pkg_version("AnnotationDbi"), GO.db = pkg_version("GO.db"),
                  orgdb = if (!is.null(route)) pkg_version(route$package) else NULL,
                  reactome.db = pkg_version("reactome.db"), ggplot2 = pkg_version("ggplot2")),
  R_version = as.character(getRversion()))
writeLines(jsonlite::toJSON(prov, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null", digits = NA),
           file.path(outdir, "sets_provenance.json"))
say("done: sets_provenance.json, sets_methods.txt and %d Sets_*.csv in %s", length(contrast_rec) +
      sum(vapply(ip_rec$baits %||% list(), function(r) length(r$contrasts), 0L)), normalizePath(outdir))
