# =============================================================================
# set_tests.R  --  protein-SET tests on run_de.R's model. Sourced by run_sets.R (and the
# simulation sim_ip_normaliser.R); functions only, no command line.
#
# Two tests, always side by side, never one alone:
#   camera  COMPETITIVE: are the set's proteins more changed than the other proteins?
#           Wu D, Smyth GK (2012) Nucleic Acids Res 40(17):e133.
#   fry     SELF-CONTAINED: are the set's proteins changed at all? (roast with
#           nrot = Inf; Giner G, Smyth GK (2016) F1000Research 5:2605.)
# They answer different questions and do disagree: in the Dog CSF session (2026-09-21) a CNS
# set was camera-significant while fry gave p = 0.53 -- the set moved relative to the rest,
# the rest moved too. A reader needs both to know which is the case.
#
# What is NOT here, on purpose (the Dog CSF session's lessons):
#   geneSetTest  treats proteins as independent. Co-regulated proteins are not, so its
#                p-values for exactly the sets that matter are far too small.
#   hand lists   a list drawn up after looking at the results is selected on the outcome.
#                Testing it on the same data finds what selected it. Only sets fixed
#                BEFORE the results are tested: GO / Reactome releases, or a user GMT whose
#                file predates the DE results (run_sets.R records its sha256 and dates).
#   depth        a set can track run depth rather than biology (a complement-set p of 2e-4
#                became 0.81 once log run depth was a covariate). Every test is repeated
#                with depth in the model and a set that does not survive is marked.
#
# One model: the statistics camera ranks are run_de.R's own moderated t for the contrast --
# the numbers in the DE table -- and fry and the residual correlations use the same
# expression, precision weights, design, contrast, block and correlation. Nothing is refitted
# for the primary tests. The depth sensitivity analysis and the pulldown bait normalisation
# refit with fit_like_run_de(), which run_sets.R first checks reproduces run_de's t exactly.
# =============================================================================

SET_DEFAULTS <- list(
  min_size = 10L,          # members quantified in the tested proteins
  max_size = 500L,
  trigger = 10L,           # run_de suggests set tests below this many significant proteins
  presence_majority = 0.5, # a set is "mostly presence calls" above this fraction
  ontologies = c("BP", "CC", "MF"),
  ref_seed = 20260928L,    # the canonical residual basis's fixed matrix (lm_effects)
  ref_seeds_check = 1:4,   # other fixed matrices: a fry call must hold on every one of them
  ip_min_enrichment = 2,   # pulldown interactome: pooled log2 enrichment over the control (4-fold)
  ip_weak_quantile = 0.25, # "weak" interactor: the interactome's bottom quarter of pooled enrichment
  broad_protein_share = 0.6, # "most proteins moved together" needs this share of proteins that way
  camera_low_power = 0.05) # whole-run correlation of random sets above which camera has little power

PRESENCE_CALL_LABEL <- "presence call"   # run_de.R EVIDENCE_LABELS[["none_in_a_group"]]

`%||%` <- function(a, b) if (is.null(a)) b else a

# What each column of Sets_<contrast>.csv / Sets_baitnorm_<contrast>.csv holds. Written into
# sets_provenance.json; AGENTS.md, the analysis brief and the report quote it from there.
SET_COLUMNS <- list(
  Set_ID = "the set's identifier in its source (GO ID, Reactome stable ID, GMT set name)",
  Set_Name = "the set's name", Source = "GO_BP / GO_CC / GO_MF / Reactome / GMT:<file>",
  NGenes = "the set's proteins among those tested (protein groups; a group counts when any of its genes is in the set)",
  Mean_logFC = "the set's effect: mean log2 fold change of its proteins (Sets_baitnorm: relative to the reference)",
  Call = "'camera + fry', 'camera only', 'fry only', or empty: which test is significant (FDR)",
  camera_Direction = "Up / Down: the direction camera tested", camera_Correlation = "the inter-protein correlation camera used for the set",
  camera_PValue = "camera (competitive) p-value", camera_FDR = "Benjamini-Hochberg across this contrast's sets",
  fry_Direction = "Up / Down", fry_PValue = "fry (self-contained) p-value", fry_FDR = "Benjamini-Hochberg across this contrast's sets",
  "depth_*" = "the same tests with log2 run depth in the model (empty when not run); depth_Mean_logFC is the set's effect with depth in the model",
  camera_Depth = "for a camera-significant set: 'holds' with run depth in the model, or 'lost'",
  fry_Depth = "for a fry-significant set: 'holds' or 'lost'",
  fry_Reference_Stable = "for a fry-significant set: 'stable' when the call holds on every reference basis tried, else 'k of K'",
  Weak_Fraction = "Sets_baitnorm only: fraction of the set's proteins in the interactome's bottom quarter of pooled enrichment over the control",
  Presence_Call_Fraction = "fraction of the set's proteins never measured in one of the compared groups (the DE table's Evidence = presence call)",
  Mostly_Presence_Calls = "TRUE when that fraction is above one half",
  "measured_*" = "both tests with presence-call proteins removed from the sets and the background (NA: too few measured proteins)",
  "bait_* / interactome_*" = "Sets_baitnorm only: the same tests with the other reference (cross-check)",
  Reference_Robust = "Sets_baitnorm only: a significant set 'holds with the other reference' or 'depends on the reference'",
  Flags = "reasons to read the set with caution",
  Members = "Sets_baitnorm only: gene symbols of the set's interactome proteins (every other set's are in sets_members.csv)")

# How to read the results: one definition, quoted by the report, the analysis brief and AGENTS.md.
SET_READING <- c(
  paste("camera is competitive (are the set's proteins MORE changed than the other proteins?) and fry",
        "self-contained (are they changed at all?). A set's Call is the union of the two tests' BH lists;",
        "there is no correction across comparisons. fry alone: the set moved, but no more than proteins in",
        "general -- often part of a broad or global shift. camera alone: the set stands out from the rest,",
        "but the evidence that it moved at all is weaker."),
  paste("camera has little power when whole runs move together (a pulldown's recovery, a sample's",
        "loading): random protein sets then correlate across runs, and camera inflates each set's",
        "variance by that correlation. The record gives it per comparison (camera_global_correlation); where",
        "camera_low_power is true, 'camera found none' is not a result, and 'significant in both' is not a",
        "useful headline."),
  paste("A set test says a predefined group of proteins shifted together. It does not say a pathway is",
        "activated or inhibited, and not which proteins drive it: look at the set's proteins in the DE table."),
  paste("GO terms nest and overlap (a protein in a term is in all of its ancestors), so neighbouring",
        "significant terms are largely the same proteins, not independent confirmations."),
  paste("'not separable from run depth': significant without run depth in the model but not with it --",
        "the compared runs differ in depth and the set's change cannot be told apart from that (the depth_*",
        "columns give the p-value and effect with depth). It is not proof of an artefact: a depth difference",
        "can itself be biology. When most sets of a comparison move together and are not separable from",
        "depth, say exactly that; do not pick out the sets that happen to survive as a signature -- they are",
        "the largest movers of the same shift."),
  paste("'fry call depends on the residual basis': fry's robust variance takes each protein's largest",
        "residual coordinate, which depends on an arbitrary reference basis; the call did not hold on every",
        "reference tried. Treat it as borderline. camera does not depend on the reference."),
  paste("'mostly presence calls': most of the set's proteins were never measured in one of the compared",
        "groups, so its shift rests on the detection model -- read the measured_* columns."),
  paste("Pulldowns: a bait-vs-control set test is about what co-purifies with the bait, not abundance in the",
        "tissue; control-vs-control (e.g. IgG Old vs Young) is about the lysate background every IP carries."),
  paste("Bait-normalised tests (Sets_baitnorm_*) ask whether a set changed relative to the bait's complex",
        "(its interactome median; the bait protein as a cross-check). The interactome is chosen on the",
        "bait's enrichment over the control pooled over conditions, which is independent of the",
        "between-condition comparison -- but the minimum enrichment makes it one-sided: a protein that",
        "falls in one condition has a lower pooled enrichment and can drop out before the test, so losses",
        "from the complex are harder to detect than gains (each bait's reading says how much). Nothing",
        "found is never evidence that the complex lost nothing. 'depends on the reference': the two",
        "references disagree -- report the set as unresolved, with its direction, its effect and how far",
        "apart the references put the recovery. 'mostly weak interactors': most of the set's proteins are",
        "in the interactome's bottom quarter of enrichment over the control, where a change in the",
        "control background weighs most."),
  paste("'between-block contrast from the blocked fit': a comparison between different subjects (mice)",
        "that averages several samples per subject; its p-values are anti-conservative for sets sharing",
        "subject-level variation (the DE record's CAUTION)."),
  paste("Nothing significant means no coordinated shift was detectable with this design and these samples",
        "(the reading gives n per side and the residual df), not that nothing changed."))

# ---- organisms ------------------------------------------------------------------
# Bioconductor OrgDb code, name, NCBI taxid, Reactome species code.
SET_ORGANISMS <- data.frame(
  code     = c("Hs", "Mm", "Rn", "Cf", "Bt", "Ss", "Gg", "Dr", "Dm", "Ce"),
  name     = c("human", "mouse", "rat", "dog", "cow", "pig", "chicken", "zebrafish", "fly", "worm"),
  latin    = c("Homo sapiens", "Mus musculus", "Rattus norvegicus", "Canis familiaris",
               "Bos taurus", "Sus scrofa", "Gallus gallus", "Danio rerio",
               "Drosophila melanogaster", "Caenorhabditis elegans"),
  taxid    = c(9606, 10090, 10116, 9615, 9913, 9823, 9031, 7955, 7227, 6239),
  reactome = c("HSA", "MMU", "RNO", "CFA", "BTA", "SSC", "GGA", "DRE", "DME", "CEL"),
  stringsAsFactors = FALSE)

# The organism from --organism (name, code, latin name or taxid), else the taxid the search
# FASTA's .meta.json records. Neither: stop -- a guessed organism maps nothing, quietly.
resolve_organism <- function(arg = NULL, fasta_meta = NULL) {
  pick <- function(x) {
    x <- tolower(trimws(as.character(x)))
    i <- which(tolower(SET_ORGANISMS$name) == x | tolower(SET_ORGANISMS$code) == x |
               tolower(SET_ORGANISMS$latin) == x | as.character(SET_ORGANISMS$taxid) == x)
    if (length(i)) SET_ORGANISMS[i[1], ] else NULL
  }
  if (!is.null(arg)) {
    o <- pick(arg)
    if (!is.null(o)) return(c(as.list(o), source = "--organism"))
    # an organism without an OrgDb here: its gene symbols are matched to human's
    return(list(code = NA_character_, name = arg, latin = arg, taxid = NA_real_,
                reactome = NA_character_, source = "--organism"))
  }
  if (!is.null(fasta_meta) && file.exists(fasta_meta) && requireNamespace("jsonlite", quietly = TRUE)) {
    m <- jsonlite::fromJSON(fasta_meta, simplifyVector = FALSE)
    if (!is.null(m$taxid)) {
      o <- pick(m$taxid)
      if (!is.null(o)) return(c(as.list(o), source = paste0("taxid ", m$taxid, " in ", basename(fasta_meta))))
      return(list(code = NA_character_, name = m$organism %||% as.character(m$taxid),
                  latin = m$organism %||% NA_character_, taxid = m$taxid, reactome = NA_character_,
                  source = paste0("taxid ", m$taxid, " in ", basename(fasta_meta))))
    }
  }
  stop("run_sets.R: which organism? Give --organism (e.g. mouse, human, 10090), or run_de.R ",
       "--fasta-meta so the search FASTA's taxid is on record.")
}

# Which OrgDb maps the gene symbols. The organism's own when installed; otherwise human's,
# matching symbols by identity after upper-casing -- a proxy for 1:1 orthology that misses
# renamed genes and lineage-specific families, so the route and its mapping rate are recorded.
annotation_route <- function(org) {
  own <- if (!is.na(org$code)) paste0("org.", org$code, ".eg.db") else NA_character_
  if (!is.na(own) && requireNamespace(own, quietly = TRUE))
    return(list(package = own, route = "native", uppercase = FALSE,
                note = sprintf("%s gene symbols -> Entrez Gene IDs (%s)", org$name, own)))
  if (requireNamespace("org.Hs.eg.db", quietly = TRUE))
    return(list(package = "org.Hs.eg.db", route = "human symbol identity", uppercase = TRUE,
                note = sprintf(paste0("%s has no annotation package here%s: its gene symbols, ",
                                      "upper-cased, were matched to HUMAN symbols (org.Hs.eg.db). ",
                                      "Symbol identity stands in for 1:1 orthology; renamed genes ",
                                      "and lineage-expanded families do not map"),
                               org$name, if (is.na(own)) "" else paste0(" (", own, " not installed)"))))
  stop("run_sets.R: no annotation package. Install ", if (is.na(own)) "" else paste0(own, " or "),
       "org.Hs.eg.db (bash setup.sh installs org.Hs.eg.db, org.Mm.eg.db and GO.db).")
}

reactome_info <- function() {
  i <- reactome.db::reactome_dbInfo()
  list(package = "reactome.db", version = pkg_version("reactome.db"),
       SOURCEDATE = i$value[i$name == "SOURCEDATE"], SOURCEURL = i$value[i$name == "SOURCEURL"])
}

# a path relative to another folder (both exist), for records that move with the session
relative_path <- function(path, from) {
  p <- strsplit(normalizePath(path), "/", fixed = TRUE)[[1]]
  f <- strsplit(normalizePath(from), "/", fixed = TRUE)[[1]]
  n <- 0L
  while (n < min(length(p), length(f)) && p[n + 1] == f[n + 1]) n <- n + 1L
  paste(c(rep("..", length(f) - n), p[-seq_len(n)]), collapse = "/")
}

pkg_version <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)

# The OrgDb's own record of its sources (GO release date, Entrez date) -- the version of a
# set collection is its data date, not the package number alone.
orgdb_info <- function(pkg) {
  db <- getExportedValue(pkg, pkg)
  info <- AnnotationDbi::metadata(db)
  keep <- c("GOSOURCEDATE", "GOSOURCEURL", "EGSOURCEDATE", "GOEGSOURCEDATE", "ORGANISM")
  out <- setNames(as.list(info$value[match(keep, info$name)]), keep)
  c(list(package = pkg, version = pkg_version(pkg)), out[!vapply(out, is.na, TRUE)])
}

# ---- proteins -> gene identifiers --------------------------------------------------
# A protein group's Genes field may list several genes ("Jph3;Jph4"): the group belongs to a
# set when ANY of its genes does. Symbols first, then unambiguous aliases (one Entrez ID).
split_gene_field <- function(genes) {
  s <- strsplit(ifelse(is.na(genes), "", as.character(genes)), ";", fixed = TRUE)
  lapply(s, function(x) unique(trimws(x[nzchar(trimws(x))])))
}

map_proteins_to_entrez <- function(genes, route) {
  syms <- split_gene_field(genes)
  db <- getExportedValue(route$package, route$package)
  u <- unique(unlist(syms))
  key <- if (isTRUE(route$uppercase)) toupper(u) else u
  valid_sym <- AnnotationDbi::keys(db, "SYMBOL")
  eg <- setNames(rep(NA_character_, length(u)), u)
  ks <- key %in% valid_sym
  if (any(ks))
    eg[ks] <- suppressMessages(AnnotationDbi::mapIds(db, keys = key[ks], column = "ENTREZID",
                                                     keytype = "SYMBOL", multiVals = "first"))
  n_alias <- 0L
  miss <- is.na(eg) & key %in% AnnotationDbi::keys(db, "ALIAS")
  if (any(miss)) {
    al <- suppressMessages(AnnotationDbi::select(db, keys = unique(key[miss]), columns = "ENTREZID",
                                                 keytype = "ALIAS"))
    al <- unique(al[!is.na(al$ENTREZID), c("ALIAS", "ENTREZID")])
    one <- names(which(table(al$ALIAS) == 1))      # an alias naming two genes names neither
    hit <- miss & key %in% one
    eg[hit] <- al$ENTREZID[match(key[hit], al$ALIAS)]
    n_alias <- sum(hit)
  }
  per_protein <- lapply(syms, function(s) unique(stats::na.omit(unname(eg[s]))))
  has_sym <- lengths(syms) > 0
  list(ids = per_protein,
       n_proteins = length(syms), n_with_symbol = sum(has_sym),
       n_mapped = sum(lengths(per_protein) > 0),
       mapping_rate = if (sum(has_sym)) sum(lengths(per_protein) > 0) / sum(has_sym) else NA_real_,
       n_symbols = length(u), n_symbols_mapped = sum(!is.na(eg)), n_via_alias = n_alias)
}

# ---- set collections (a-priori only) ----------------------------------------------
# GO with ancestor propagation -- a gene annotated to a term belongs to every ancestor --
# from the OrgDb's GO2ALLEGS map, the map limma::goana uses. Term names from GO.db.
go_sets <- function(pkg, universe_entrez, ontologies = SET_DEFAULTS$ontologies) {
  if (!requireNamespace("GO.db", quietly = TRUE))
    stop("run_sets.R: GO.db is not installed (bash setup.sh installs it).")
  obj <- utils::getFromNamespace(paste0(sub("\\.db$", "", pkg), "GO2ALLEGS"), pkg)
  AnnotationDbi::Lkeys(obj) <- intersect(unique(universe_entrez), AnnotationDbi::Lkeys(obj))
  tab <- AnnotationDbi::toTable(obj)[, c("gene_id", "go_id", "Ontology")]
  tab <- unique(tab[tab$Ontology %in% ontologies, ])
  sets <- split(tab$gene_id, tab$go_id)
  ont <- tab$Ontology[match(names(sets), tab$go_id)]
  term <- suppressMessages(AnnotationDbi::Term(names(sets)))
  list(sets = sets, id = names(sets), name = unname(term[names(sets)]),
       source = paste0("GO_", ont))
}

# Reactome pathways for the organism (reactome.db, Entrez IDs), when installed.
reactome_sets <- function(org, universe_entrez) {
  if (!requireNamespace("reactome.db", quietly = TRUE))
    stop("run_sets.R --sets reactome: reactome.db is not installed.")
  if (is.na(org$reactome))
    stop("run_sets.R --sets reactome: no Reactome species code for ", org$name)
  p2g <- AnnotationDbi::as.list(reactome.db::reactomePATHID2EXTID)
  p2g <- p2g[startsWith(names(p2g), paste0("R-", org$reactome, "-"))]
  p2g <- lapply(p2g, function(g) intersect(as.character(g), universe_entrez))
  p2g <- p2g[lengths(p2g) > 0]
  nm <- AnnotationDbi::mget(names(p2g), reactome.db::reactomePATHID2NAME, ifnotfound = NA)
  nm <- sub(paste0("^", org$latin, ": "), "", vapply(nm, function(x) as.character(x[1]), ""))
  list(sets = p2g, id = names(p2g), name = unname(nm), source = rep("Reactome", length(p2g)))
}

# A user GMT (name <tab> description <tab> gene ...). Identifiers are Entrez IDs when every
# one is a number, else gene symbols, matched to the proteins' own symbols case-insensitively.
read_gmt <- function(path) {
  ln <- readLines(path, warn = FALSE, encoding = "UTF-8")
  ln <- ln[nzchar(trimws(ln))]
  f <- strsplit(ln, "\t", fixed = TRUE)
  bad <- lengths(f) < 3
  if (any(bad)) stop("GMT ", basename(path), ": line(s) ", paste(which(bad), collapse = ", "),
                     " have fewer than 3 tab-separated fields (name, description, genes)")
  ids <- lapply(f, function(x) unique(trimws(x[-(1:2)])))
  ids <- lapply(ids, function(x) x[nzchar(x)])
  entrez <- all(grepl("^[0-9]+$", unlist(ids)))
  nm <- vapply(f, `[`, "", 1)
  if (anyDuplicated(nm)) stop("GMT ", basename(path), ": duplicated set names: ",
                              paste(unique(nm[duplicated(nm)]), collapse = ", "))
  list(sets = setNames(if (entrez) ids else lapply(ids, function(x) unique(toupper(x))), nm), id = nm,
       name = nm, description = vapply(f, `[`, "", 2),     # MSigDB's description is a URL
       source = rep(paste0("GMT:", basename(path)), length(f)),
       id_type = if (entrez) "Entrez Gene ID" else "gene symbol (case-insensitive)")
}

# ---- sets -> row indices of the tested proteins ----------------------------------
# prot_ids: per protein row, the identifiers it carries (Entrez IDs, or upper-cased symbols
# for a symbol GMT). Size limits count the proteins a set has AMONG THE TESTED ONES, and a set
# must leave at least min_size tested proteins outside it: a competitive test of a set that is
# (nearly) the whole universe has nothing to compete with (camera's p is undefined at m = G) --
# it matters for a pulldown's interactome, a universe of a few dozen proteins.
index_sets <- function(collection, prot_ids, min_size, max_size) {
  max_size <- min(max_size, length(prot_ids) - min_size)
  long <- data.frame(row = rep(seq_along(prot_ids), lengths(prot_ids)),
                     id = unlist(prot_ids, use.names = FALSE), stringsAsFactors = FALSE)
  rows_by_id <- split(long$row, long$id)
  idx <- lapply(collection$sets, function(g)
    sort(unique(unlist(rows_by_id[intersect(g, names(rows_by_id))], use.names = FALSE))))
  n <- lengths(idx)
  keep <- n >= min_size & n <= max_size
  list(index = setNames(idx[keep], collection$id[keep]), name = collection$name[keep],
       source = collection$source[keep],
       counts = list(defined = length(idx), with_members = sum(n > 0), tested = sum(keep),
                     below_min = sum(n > 0 & n < min_size), above_max = sum(n > max_size),
                     max_size_applied = max_size))
}

# ---- the tests ---------------------------------------------------------------------
# camera, on the fit's own statistics. limma::camera() refits y ~ design itself and takes no
# block; with a random block that refit would not be run_de's model (and limpa's df
# correction for wholly imputed groups would be lost). This is camera.default's algorithm
# (limma 3.68) with its two inputs taken from run_de instead:
#   stat  the contrast's moderated t from run_de's fit, as camera converts it
#         (zscoreT, Hill's approximation, df = the fit's df.total);
#   U     the standardised residual effects of the same model: weights, design, contrast,
#         and -- for a random block -- the block-correlation whitening fry uses (lm_effects).
# inter_gene_cor = NA (default) estimates the correlation for each set from U: camera's
# "rigorous error rate control" mode, df = min(residual df, G - 2). A number is camera's
# preset mode (limma's default 0.01), df = G - 2. Either way VIF >= 1 (allow.neg.cor = FALSE).
# The Smyth-lab route for blocked designs (support.bioconductor.org/p/78299: "you can't use
# duplicateCorrelation() with camera()") is cameraPR on the blocked fit's t; for a preset
# correlation this function IS cameraPR on z-scores, and with no block and no weights it equals
# camera() exactly -- tests/test_set_tests.py checks both. (With weights, on lm_effects'
# canonical basis: see there.)
z_from_t <- function(t, df_total) limma::zscoreT(t, df = df_total, approx = TRUE, method = "hill")

# limma internals fry is built from, used here as fry uses them (checked against limma 3.68).
# A limma without them stops run_sets.R with this message rather than failing mid-run.
check_limma_internals <- function() {
  need <- c(".lmEffects", ".fryEffects", ".squeezeVar")
  miss <- need[!vapply(need, exists, TRUE, envir = asNamespace("limma"), inherits = FALSE)]
  if (length(miss))
    stop("run_sets.R: this limma (", utils::packageVersion("limma"), ") has no ",
         paste(miss, collapse = ", "), " -- set_tests.R follows limma 3.68's fry(); check it ",
         "against this limma's fry.default and camera.default before using it.")
  invisible(TRUE)
}

# The model's effects for one contrast: column 1 the contrast effect (limma:::.lmEffects, which
# fry uses), the rest the residual effects, per protein, on ONE canonical basis.
#
# camera's correlation estimate and fry's set statistic average a set's residual effects
# coordinate by coordinate, and fry's robust variance takes each protein's largest squared
# coordinate. limma takes the residual basis from a Householder QR -- of the design, or with
# precision weights of each protein's own weighted design -- and that basis is not unique: on
# PROT_0756 R's qr() of the SAME design matrix came out with different reflection signs on an
# AMD EPYC 9554 HIVE node than on EPYC 7532/7763 nodes (floating-point zeros in a 0/1 design
# landing either side of zero), so residual effects were rotated -- 1,184 of 6,231 proteins
# with weights -- and camera and fry p-values moved by up to 0.09 between otherwise identical
# runs. Here the basis is canonical:
#   reference  the polar factor of the residual projector (which depends only on the design's
#              column space) applied to a fixed generic matrix -- unique, so the same on every
#              machine and for every way of writing the design;
#   a protein  with weights: the polar factor of the reference projected into its own weighted
#              residual space -- the orthonormal basis closest to the reference, so coordinate
#              k is as nearly as possible the same sample contrast for every protein; with
#              constant weights, the reference itself.
# Residual sums of squares are unchanged. Changing the reference rotates every protein's basis by
# the SAME rotation, with weights too (the polar factor carries a common rotation through), and
# camera's statistics are invariant to a common rotation: camera does not depend on the reference
# (to 1e-15 with weights, review check1b). fry depends on it only through its robust variance's
# largest-coordinate term: up to 0.03 in p without weights, 0.07 with. Separately from the
# reference, both differ from limma's own functions where limma's basis differs: without weights
# camera is limma::camera() exactly and fry is limma::fry() exactly once that one term is taken
# from limma's basis (the tests); with weights limma's per-protein bases are not aligned, and
# camera differs from weighted limma::camera by up to 0.04 in p, fry from weighted limma::fry by
# up to 0.1 -- a difference from limma, not a dependence on the reference.
# "Canonical" is not "right": the reference is one fixed choice, so run_sets.R re-tests every fry
# call on other references (SET_DEFAULTS$ref_seeds_check) and flags calls that depend on it.
polar <- function(B) { sv <- svd(B); sv$u %*% t(sv$v) }
fixed_matrix <- function(n, k, seed = SET_DEFAULTS$ref_seed) {  # generic; R's RNG is portable
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
          assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
  matrix(stats::rnorm(n * k), n, k)
}
lm_effects <- function(E, design, contrast, weights = NULL, block = NULL, correlation = NULL,
                       ref_seed = SET_DEFAULTS$ref_seed) {
  eff <- limma:::.lmEffects(E, design = design, contrast = contrast, weights = weights,
                            block = block, correlation = correlation)
  n <- ncol(E); p <- ncol(design)
  whiten <- function(M) M
  if (!is.null(block)) {                       # as .lmEffects: weights first, then the block
    Z <- outer(block, unique(block), `==`)
    cm <- Z %*% (correlation * t(Z)); diag(cm) <- 1
    R <- chol(cm)
    whiten <- function(M) backsolve(R, M, transpose = TRUE)
  }
  X0 <- whiten(design)
  B0 <- qr.resid(qr(X0), fixed_matrix(n, n - p, ref_seed))
  if (min(svd(B0, nu = 0, nv = 0)$d) < 1e-6)
    stop("lm_effects: no canonical residual basis for this design (degenerate projection)")
  ref <- polar(B0)
  if (is.null(weights)) {
    U <- t(whiten(t(E))) %*% ref
  } else {
    sw <- sqrt(weights)
    U <- matrix(0, nrow(E), n - p)
    for (g in seq_len(nrow(E))) {
      q <- qr(whiten(design * sw[g, ]))
      U[g, ] <- crossprod(polar(qr.resid(q, ref)), qr.resid(q, whiten(E[g, ] * sw[g, ])))
    }
  }
  if (!isTRUE(all.equal(unname(rowSums(U^2)), unname(rowSums(eff[, -1, drop = FALSE]^2)), tolerance = 1e-8)))
    stop("lm_effects: the residual effects lost residual variance -- a protein's weighted design ",
         "is not of full rank")
  eff[, -1] <- U
  eff
}

# camera's residuals: each protein's residual effects, standardised (camera.default).
residual_U <- function(eff) {
  U <- eff[, -1, drop = FALSE]
  U / sqrt(pmax(rowMeans(U^2), 1e-8))
}

camera_fit <- function(stat, U, index, df_residual, inter_gene_cor = NA) {
  G <- length(stat)
  fixed <- !is.na(inter_gene_cor)
  df_camera <- if (fixed) G - 2L else min(df_residual, G - 2L)
  mean_stat <- mean(stat); var_stat <- stats::var(stat)
  res <- vapply(index, function(iset) {
    m <- length(iset); m2 <- G - m
    if (fixed) { cor <- inter_gene_cor; vif <- 1 + (m - 1) * cor }
    else if (m > 1) { vif <- m * mean(colMeans(U[iset, , drop = FALSE])^2); cor <- (vif - 1) / (m - 1) }
    else { vif <- 1; cor <- NA_real_ }
    vif <- max(1, vif)
    delta <- G / m2 * (mean(stat[iset]) - mean_stat)
    var_pooled <- ((G - 1) * var_stat - delta^2 * m * m2 / G) / (G - 2)
    tt <- delta / sqrt(var_pooled * (vif / m + 1 / m2))
    c(m, cor, stats::pt(tt, df_camera), stats::pt(tt, df_camera, lower.tail = FALSE))
  }, numeric(4))
  data.frame(NGenes = res[1, ], Correlation = res[2, ],
             Direction = ifelse(res[3, ] < res[4, ], "Down", "Up"),
             PValue = 2 * pmin(res[3, ], res[4, ]), row.names = names(index),
             stringsAsFactors = FALSE)
}

# fry on the model's effects: fry.default's standardize = "posterior.sd" (limma 3.68, no trend,
# not robust) and its set statistics (.fryEffects), given lm_effects()' effects -- so the
# weights, the block and its correlation enter as fry(..., weights, block, correlation) would
# take them, on the canonical basis above.
# u2_eff: the effects whose largest squared coordinate the robust variance uses (the only
# basis-dependent step without weights); lm_effects' own by default.
fry_effects <- function(eff, index, u2_eff = eff) {
  df_residual <- ncol(eff) - 1
  gq <- statmod::gauss.quad.prob(128, "uniform")
  Eu2max <- sum((df_residual + 1) * gq$nodes^df_residual * stats::qchisq(gq$nodes, df = 1) * gq$weights)
  s2_robust <- (rowSums(eff^2) - apply(u2_eff^2, 1, max)) / (df_residual + 1 - Eu2max)
  fit <- limma::fitFDist(rowMeans(eff[, -1, drop = FALSE]^2), df1 = df_residual)
  s2_robust <- limma:::.squeezeVar(s2_robust, df = 0.92 * df_residual, var.prior = fit$scale,
                                   df.prior = fit$df2)
  r <- limma:::.fryEffects(eff / sqrt(s2_robust), index = index, sort = "none")
  data.frame(NGenes = r$NGenes, Direction = as.character(r$Direction), PValue = r$PValue,
             row.names = names(index), stringsAsFactors = FALSE)
}
fry_fit <- function(E, index, design, contrast, weights = NULL, block = NULL, correlation = NULL)
  fry_effects(lm_effects(E, design, contrast, weights, block, correlation), index)

# Both tests on one model, BH across the sets of this contrast (one family per test).
test_sets <- function(model, index, inter_gene_cor = NA) {
  eff <- lm_effects(model$E, model$design, model$contrast, model$weights, model$block, model$correlation)
  U <- residual_U(eff)
  cam <- camera_fit(model$z, U, index, df_residual = ncol(U), inter_gene_cor = inter_gene_cor)
  fr <- fry_effects(eff, index)
  out <- data.frame(NGenes = cam$NGenes,
             camera_Direction = cam$Direction, camera_Correlation = cam$Correlation,
             camera_PValue = cam$PValue, camera_FDR = stats::p.adjust(cam$PValue, "BH"),
             fry_Direction = fr$Direction, fry_PValue = fr$PValue,
             fry_FDR = stats::p.adjust(fr$PValue, "BH"),
             row.names = names(index), stringsAsFactors = FALSE)
  attr(out, "df_residual") <- ncol(U)
  attr(out, "global_correlation") <- global_correlation(U)
  out
}

# fry on the same model with other reference bases (the call-stability check): for each set, on
# how many of them the call (FDR < fdr, same direction as `direction`) holds.
fry_calls_on_references <- function(model, index, direction, fdr, seeds = SET_DEFAULTS$ref_seeds_check) {
  held <- vapply(seeds, function(sd) {
    r <- fry_effects(lm_effects(model$E, model$design, model$contrast, model$weights, model$block,
                                model$correlation, ref_seed = sd), index)
    stats::p.adjust(r$PValue, "BH") < fdr & r$Direction == direction
  }, logical(length(index)))
  list(held = rowSums(matrix(held, nrow = length(index))), of = length(seeds))
}

# How strongly whole runs move together: camera's estimated correlation for RANDOM protein sets
# (sizes 20, 50, 100; 100 draws each, from a fixed seed). camera treats it as within-set
# correlation and inflates every set's variance by it -- the variance factor at 50 proteins is
# 1 + 49 x this -- so above SET_DEFAULTS$camera_low_power camera has little power.
global_correlation <- function(U, sizes = c(20L, 50L, 100L), draws = 100L) {
  G <- nrow(U)
  sizes <- sizes[sizes < G]
  if (!length(sizes)) return(NA_real_)
  idx <- local({
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
            assign(".Random.seed", old, envir = globalenv()))
    set.seed(SET_DEFAULTS$ref_seed, kind = "Mersenne-Twister", sample.kind = "Rejection")
    unlist(lapply(sizes, function(k) replicate(draws, sample(G, k), simplify = FALSE)), recursive = FALSE)
  })
  mean(vapply(idx, function(i) { m <- length(i); (m * mean(colMeans(U[i, , drop = FALSE])^2) - 1) / (m - 1) }, 0))
}

# The call for each set: which of the two tests is significant at `fdr`.
set_call <- function(camera_fdr, fry_fdr, fdr)
  ifelse(camera_fdr < fdr & fry_fdr < fdr, "camera + fry",
         ifelse(camera_fdr < fdr, "camera only", ifelse(fry_fdr < fdr, "fry only", "")))
# Does a significant result survive a sensitivity refit (depth)? "holds": the same test is
# significant in the same direction; "lost": it is not; "": not significant to begin with.
set_holds <- function(fdr0, dir0, fdr1, dir1, fdr)
  ifelse(!(fdr0 < fdr), "", ifelse(!is.na(fdr1) & fdr1 < fdr & dir1 == dir0, "holds", "lost"))

# The model for one contrast restricted to some protein rows (the measured-only variant, the
# pulldown interactome): same everything, fewer rows.
subset_model <- function(model, rows) {
  model$E <- model$E[rows, , drop = FALSE]
  if (!is.null(model$weights)) model$weights <- model$weights[rows, , drop = FALSE]
  model$z <- model$z[rows]
  model
}
reindex <- function(index, rows) lapply(index, function(i) { j <- match(i, rows); sort(j[!is.na(j)]) })

# ---- refitting like run_de (depth covariate, bait normalisation) --------------------
# The model run_de.R fits, as one call: dpc -> limpa::dpcDE (voomaLmFitWithImputation, block
# passed through); maxlfq -> lmFit, with duplicateCorrelation's consensus for a random block.
fit_like_run_de <- function(method, y, design, block = NULL) {
  if (method == "dpc") {
    fit <- suppressMessages(limpa::dpcDE(y, design, plot = FALSE, block = block))
    return(list(fit = fit, E = fit$EList$E, weights = fit$EList$weights,
                correlation = fit$correlation))
  }
  E <- if (is.list(y)) y$E else y
  if (is.null(block)) return(list(fit = limma::lmFit(E, design), E = E, weights = NULL, correlation = NULL))
  rho <- limma::duplicateCorrelation(E, design, block = block)$consensus.correlation
  list(fit = limma::lmFit(E, design, block = block, correlation = rho), E = E, weights = NULL,
       correlation = rho)
}
contrast_stats <- function(fit, contrast) {
  fc <- limma::eBayes(limma::contrasts.fit(fit, matrix(contrast, ncol = 1,
                                                       dimnames = list(names(contrast), "c"))))
  list(t = fc$t[, 1], df_total = fc$df.total, logFC = fc$coefficients[, 1])
}

# ---- depth ---------------------------------------------------------------------------
# Run depth, one definition: Detection_Matrix.csv's column sums -- precursors observed per
# run (dpc) or proteins quantified per run (maxlfq) -- on the log2 scale, centred.
run_depth <- function(det_n, zero_means) {
  n <- colSums(det_n, na.rm = TRUE)
  list(value = log2(pmax(n, 1)) - mean(log2(pmax(n, 1))), counts = n,
       definition = if (identical(zero_means, "inferred"))
         "log2 precursors observed per run (Detection_Matrix.csv column sums), centred"
       else "log2 proteins quantified per run (Detection_Matrix.csv column sums), centred")
}
# With `role` (a pulldown's run type per run: bait / control), one depth slope per run type,
# each centred within its runs: a control IP's depth and a bait IP's depth are not the same
# thing, and one shared slope let the bait runs' slope stand in for the controls' (review S3).
add_depth <- function(design, contrast, depth, role = NULL) {
  if (is.null(role) || length(unique(role)) < 2) {
    d <- cbind(design, run_depth = depth)
    return(list(design = d, contrast = c(contrast, run_depth = 0), columns = "run_depth"))
  }
  cols <- sapply(sort(unique(role)), function(r) ifelse(role == r, depth - mean(depth[role == r]), 0))
  colnames(cols) <- paste0("run_depth_", make.names(sort(unique(role))))
  list(design = cbind(design, cols),
       contrast = c(contrast, setNames(rep(0, ncol(cols)), colnames(cols))), columns = colnames(cols))
}

# ---- presence calls (DPC) -------------------------------------------------------------
# The fraction of a set's members that run_de's Evidence column calls a presence call for
# this contrast: never measured in one of the compared groups, so that protein's difference
# is the detection model's, not a measurement. One definition: the DE table's Evidence.
presence_fraction <- function(index, evidence)
  vapply(index, function(i) if (length(i)) mean(evidence[i] %in% PRESENCE_CALL_LABEL) else NA_real_, 0)

# ---- pulldowns: the interactome ----------------------------------------------------------
# The bait's interactome is chosen on its enrichment over the control POOLED over conditions --
# (mean of the bait's groups) - (mean of their controls) -- at FDR and at least
# SET_DEFAULTS$ip_min_enrichment log2. With equal runs per condition that statistic is
# uncorrelated with the between-condition comparison under the null (independent filtering,
# Bourgon et al. 2010 PNAS 107:9546). "Enriched at EVERY condition" is not: when recovery
# differs, it keeps a protein that fell in one condition only if it was still enriched there
# (review B1: 10.7 of 20 members of a set truly down 0.8 kept, against 19.4 pooled; PROT_0756
# JPH3 had 225 Old / 974 Young enriched). The minimum enrichment keeps out weak interactors,
# whose "change" relative to the complex is mostly the control background changing (review S2:
# null sets of the weakest interactors 0.29-0.68 false positives without it, 0.087 with it; the cost,
# losses left undetectable, is below).
pooled_enrichment_contrast <- function(cmat, vs_control)
  rowMeans(cmat[, vs_control, drop = FALSE])

# What the minimum costs (review M2, W2): it acts on the enrichment POOLED over conditions, so a
# protein that falls in one condition -- binds the bait less there, or leaves the complex -- has
# a lower pooled enrichment (by half the fall) and can drop below the minimum before any test.
# Losses are harder to detect than gains, by how much depends on how close the interactome sits
# to the minimum: in sim_set_tests.R --what ip (enrichment 0.3-3 log2) a falling set is almost
# never tested; on PROT_0756, whose interactomes have a median of 2.4-2.6 log2, 54-63% of the
# tested sets would survive a 0.8 log2 fall. So each bait's reading gives its own numbers.
ip_direction_bias <- function(min_enrichment) {   # the general statement (report, brief)
  if (!is.finite(min_enrichment) || min_enrichment <= 0) return(NA_character_)
  paste0("Losses from a complex are harder to detect than gains: the interactome is chosen on the ",
         "enrichment pooled over conditions, so a protein that falls in one condition can drop below the ",
         sprintf("%g-fold", signif(2^min_enrichment, 3)), " minimum before any test. Each bait's reading ",
         "says how much.")
}
# enrich: the pooled log2 enrichment of the interactome proteins the tested sets are indexed on;
# index: the tested sets. For each fall (log2, in one condition): the share of the interactome it
# would put below the minimum, and the share of the tested sets that would still have min_size
# members if ALL their proteins fell that much.
ip_loss_sensitivity <- function(enrich, index, min_enrichment, falls = c(0.8, 1.6),
                                min_size = SET_DEFAULTS$min_size) {
  list(fall_log2 = falls,
       below_minimum = round(vapply(falls, function(d) mean(enrich < min_enrichment + d / 2), 0), 3),
       sets_still_tested = round(vapply(falls, function(d)
         mean(vapply(index, function(i) sum(enrich[i] >= min_enrichment + d / 2) >= min_size, TRUE)), 0), 3),
       n_sets = length(index))
}
ip_losses_reading <- function(ls, min_enrichment) {   # one bait's sentence, from ip_loss_sensitivity
  if (!is.finite(min_enrichment) || min_enrichment <= 0 || !ls$n_sets) return(NA_character_)
  pc <- function(x) sprintf("%.0f%%", 100 * x)
  # "about": the estimates are taken as exact, so noise near the minimum is ignored (optimistic)
  sprintf(paste0("Losses are harder to detect than gains here: a %g log2 fall in one condition would put about %s ",
                 "of this interactome below the %g-fold minimum, and about %s of the sets tested here would keep ",
                 "enough members to be tested (about %s for a %g log2 fall)."),
          ls$fall_log2[1], pc(ls$below_minimum[1]), signif(2^min_enrichment, 3), pc(ls$sets_still_tested[1]),
          pc(ls$sets_still_tested[2]), ls$fall_log2[2])
}

# "Weak" interactors, relative to the interactome itself (review M3): its bottom quarter of pooled
# enrichment over the control, where a change in the control background weighs most. A random set
# of the interactome has a quarter of its members here and is flagged (over half) rarely; a
# margin above the minimum ("within 1 log2") flagged 73-88% of PROT_0756's tested sets and every
# simulated one. enrich: the interactome's pooled log2 enrichment; returns one TRUE/FALSE each.
weak_interactors <- function(enrich, q = SET_DEFAULTS$ip_weak_quantile)
  enrich <= stats::quantile(enrich, q, names = FALSE)

# ---- pulldowns: bait normalisation -------------------------------------------------------
# Each IP's values minus a per-run offset: the median of the bait's interactome in that run
# ("interactome"), or the bait protein's own value ("bait"). sim_set_tests.R --what ip compares
# them (references/set-tests.md).
ip_offsets <- function(E, runs, members, how = c("interactome", "bait"), bait_row = NULL) {
  how <- match.arg(how)
  if (how == "bait") return(E[bait_row, runs])
  apply(E[members, runs, drop = FALSE], 2, stats::median)
}

# ---- the readings -------------------------------------------------------------------------
# What the report, the brief and AGENTS.md say about each comparison, generated here so the
# wording has one definition and cannot be strengthened on the way.
signed <- function(x, d = 2) {       # "+0.41", "-0.12"; a value that rounds to zero gets no sign
  x <- round(x, d); ifelse(x == 0, sprintf(paste0("%.", d, "f"), 0), sprintf(paste0("%+.", d, "f"), x))
}

# One comparison. r: run_sets.R's per-contrast record. Every "none" carries its n and residual
# df; a camera count is not a result where camera has little power; "a broad shift" -- most
# proteins moved together -- is said only when the proteins did (at least broad_protein_share of
# them that way, and their median), not from the fry sets alone (review W1: Old_JPH3-Young_JPH3
# had 45 sets lower in Old while 53% of proteins were lower and the median was 0.000); a shift
# that cannot be separated from run depth is said to be so -- never a named signature (B2).
contrast_reading <- function(r, fdr) {
  n <- sprintf("%d vs %d %s, %d residual df", r$n_pos, r$n_neg, r$n_words, r$df_residual)
  pos <- paste(unlist(r$pos_groups), collapse = " + ")
  low <- isTRUE(r$camera_low_power)
  power <- sprintf(paste0("random protein sets correlate at %.2f across these runs, which inflates a 50-protein ",
                          "set's variance %.0f-fold"), r$camera_global_correlation, r$camera_vif_50)
  cam <- if (low) sprintf(" camera has little power here: %s, so camera's count is not a result.", power) else ""
  caution <- if (isTRUE(r$between_block_caution))
    " Between-subject comparison from the blocked fit: its p-values are anti-conservative (the DE record's CAUTION)." else ""
  recov <- if (identical(r$ip_role, "bait_between_conditions")) {
    rv <- r$recovery
    if (!is.null(rv) && is.numeric(rv$difference_log2) && is.finite(rv$difference_log2))
      sprintf(paste0(" Not corrected for bait recovery (the %s IPs brought down %.2f log2 %s of the %s complex, ",
                     "by its %s reference): see the bait-normalised tests."), pos, abs(rv$difference_log2),
              if (rv$difference_log2 < 0) "less" else "more", rv$bait, rv$reference)
    else " Not corrected for bait recovery: see the bait-normalised tests."
  } else ""
  if (r$camera == 0 && r$fry == 0)
    return(paste0(sprintf("No set change detectable in this design (%s); this is not evidence that nothing changed.", n),
                  recov, cam, caution))
  if (identical(r$ip_role, "bait_vs_control"))
    return(paste0(sprintf("What co-purifies with the bait: fry finds %d sets enriched in the IP over the control%s.",
                          r$fry_up, if (r$fry_down) sprintf(" and %d depleted", r$fry_down) else ""),
                  if (low) " Which of them stand out from the rest of the IP cannot be told here."
                  else sprintf(" camera: %d stand out from the rest of the IP.", r$camera), cam))
  up <- r$fry_up >= r$fry_down; dom <- if (up) r$fry_up else r$fry_down
  way <- if (up) "higher" else "lower"
  many_sets <- r$fry >= 20 && dom >= 0.9 * r$fry && r$camera == 0
  share <- if (up) r$protein_share_up else r$protein_share_down
  med <- r$protein_median_logFC
  proteins_moved <- is.numeric(share) && is.numeric(med) && isTRUE(share >= SET_DEFAULTS$broad_protein_share) &&
    isTRUE(if (up) med > 0 else med < 0)
  prot <- if (is.numeric(share) && is.numeric(med))
    sprintf("%.0f%% of the proteins %s, median log2 fold change %s", 100 * share, way, signed(med, 3)) else NULL
  as_are <- if (is.numeric(share) && is.numeric(med))
    sprintf("as are %.0f%% of the proteins (median log2 fold change %s)", 100 * share, signed(med, 3)) else NULL
  txt <- if (many_sets && proteins_moved && low)
    sprintf(paste0("A broad shift: %d of fry's %d significant sets are %s in %s, %s -- most proteins moved ",
                   "together; whether any category moved more than the rest cannot be told here (camera has ",
                   "little power: %s)."), dom, r$fry, way, pos, as_are, power)
  else if (many_sets && proteins_moved)
    sprintf(paste0("A broad shift: %d of fry's %d significant sets are %s in %s, %s, and none stands out from ",
                   "the rest (camera) -- most proteins moved together, not one category."),
            dom, r$fry, way, pos, as_are)
  else if (many_sets)
    sprintf(paste0("%d of fry's %d significant sets are %s in %s (FDR < %g; %s), but the proteins as a whole did not ",
                   "move that way (%s), so this is not a shift of the whole proteome."),
            dom, r$fry, way, pos, fdr, n, prot %||% "no protein-level summary recorded")
  else sprintf("Sets significant: camera %d, fry %d, both %d (FDR < %g; %s).", r$camera, r$fry, r$both, fdr, n)
  if (many_sets && proteins_moved) cam <- ""                    # said once, in the sentence above
  if (isTRUE(r$depth_tested) && r$fry > 0 && (r$fry - r$fry_depth_holds) >= 0.5 * r$fry)
    txt <- paste0(txt, sprintf(paste0(" It cannot be separated from run depth: the compared runs differ by %.2f log2 ",
                                      "in depth, and with depth in the model %d of the %d fry sets are no longer significant."),
                               abs(r$depth_difference_log2), r$fry - r$fry_depth_holds, r$fry))
  paste0(txt, recov, cam, caution)
}

# One bait-normalised comparison. x: the record's counts plus
#   pos_group, reference / cross_check ("interactome" | "bait" | NA), n_interactome, n_prior_proteins,
#   sig: data.frame(direction "Up"/"Down", effect = Mean_logFC, depends, weak, flagged) of the
#        significant sets, recovery_difference (|difference| of the two references' recovery
#        estimates, log2; NA without a cross-check), losses (ip_losses_reading or NA).
# Calls that depend on the reference give their direction, their effect range and how far apart
# the two references put the recovery -- the uncertainty the effects have to beat (review: JPH3's
# 6 sets were -0.12 to -0.09 log2 against a 0.19 log2 disagreement).
bait_reading <- function(x) {
  ref_words <- c(interactome = "the interactome median", bait = "the bait protein")
  other_words <- c(interactome = "interactome-median", bait = "bait-protein")
  tail_txt <- paste0(if (!is.na(x$losses %||% NA)) paste0(" ", x$losses) else "",
    if (x$n_prior_proteins < 100)
      sprintf(" The variance prior rests on only %d proteins, so these p-values are less stable.", x$n_prior_proteins)
    else "",
    if (isTRUE(x$camera_low_power))
      sprintf(" camera has little power here (random sets of the interactome correlate at %.2f), so its count is not a result.",
              x$camera_global_correlation) else "")
  if (x$camera == 0 && x$fry == 0)
    return(paste0(sprintf(paste0("No change relative to the %s complex detectable in this design (%d vs %d IPs, %d ",
                                 "residual df; %d interactome proteins chosen on the pooled enrichment over the ",
                                 "control); this is not evidence that the complex is unchanged."),
                          x$bait, x$n_pos, x$n_neg, x$df_residual, x$n_interactome), tail_txt))
  s <- x$sig
  rng <- function(e) if (length(e) == 1 || diff(range(e)) < 0.005) sprintf("%s log2", signed(e[1]))
                     else sprintf("%s to %s log2", signed(min(e)), signed(max(e)))
  part <- function(d, word, where) { e <- s$effect[s$direction == d]
    if (length(e)) sprintf("%d %s%s (%s)", length(e), word, where, rng(e)) else NULL }
  where <- paste(" in", x$pos_group)                         # named once: "6 lower in Old_JPH3 (...)"
  lo <- part("Down", "lower", where)
  dirs <- c(lo, part("Up", "higher", if (is.null(lo)) where else ""))
  txt <- sprintf("Relative to the %s complex (%d vs %d IPs, %d residual df): camera %d, fry %d, both %d sets -- %s, relative to %s.",
                 x$bait, x$n_pos, x$n_neg, x$df_residual, x$camera, x$fry, x$both,
                 paste(dirs, collapse = " and "), ref_words[[x$reference]])
  nd <- sum(s$depends); ns <- nrow(s)
  if (!is.na(x$cross_check %||% NA)) {
    other <- other_words[[x$cross_check]]
    rd <- x$recovery_difference
    rd_txt <- if (is.numeric(rd) && is.finite(rd)) sprintf("%.2f log2", rd) else "an unrecorded amount"
    txt <- paste0(txt,
      if (nd == ns && ns > 0) {
        within <- is.numeric(rd) && is.finite(rd) && max(abs(s$effect)) <= rd
        sprintf(paste0(" %s holds with the %s reference, whose recovery estimate differs by %s, so %s -- unresolved."),
                if (ns == 1) "It no longer" else "None", other, rd_txt,
                if (within) "the effect is within the reference uncertainty"
                else "the call rests on the choice of reference")
      } else if (nd > 0)
        sprintf(" %d of them do not hold with the %s reference (its recovery estimate differs by %s): read those as unresolved.",
                nd, other, rd_txt)
      else sprintf(" %s with the %s reference.", if (ns == 1) "It holds" else "All hold", other))
  }
  nw <- sum(s$weak)
  if (nw) txt <- paste0(txt, sprintf(" %d %s mostly weak interactors (the interactome's least-enriched quarter).",
                                     nw, if (nw == 1) "is" else "are"))
  if (ns > 0 && !any(!s$flagged) && nd < ns) txt <- paste0(txt, " Every one is flagged, so none stands as a finding.")
  paste0(txt, tail_txt)
}
