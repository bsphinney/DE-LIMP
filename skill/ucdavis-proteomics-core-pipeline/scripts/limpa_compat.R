# =============================================================================
# limpa_compat.R -- limpa::readDIANN() across limpa versions.
#
# Why (HIVE job 24154221, 2026-09-27): run_de.R passed readDIANN(annotation.columns =)
# unconditionally. That argument is limpa >= 1.4.0 (Bioconductor 3.23). bioconda's only
# build is limpa 1.2.5 (Bioconductor 3.22), which setup.sh installed -- so every dpc DE on a
# setup.sh env died with "unused argument (annotation.columns = dpc_ann)".
#
# Read from the sources (github.com/bioc/limpa):
#   RELEASE_3_22 (1.2.5)  readDIANN(..., q.columns, q.cutoffs, extra.columns =
#                         c("Protein.Group","Protein.Names","Genes","Proteotypic"))
#   RELEASE_3_23 (1.4.x)  the same argument renamed annotation.columns (same default);
#                         q.columns / q.cutoffs unchanged, cutoffs honoured element-wise.
# So an older limpa gets the SAME annotation through extra.columns. Only a readDIANN taking
# neither (none known) reads limpa's own annotation and joins the missing columns from the
# report by Precursor.Id -- they are precursor attributes (Protein.Ids lists every
# accession ONE precursor matches, which the contaminant rule reads), so a join by
# Protein.Group would give each precursor its group's accessions, not its own.
# setup.sh installs limpa >= 1.4.0; this keeps an older env working, and says which path ran.
# =============================================================================

LIMPA_ANNOTATION_DEFAULT <- c("Protein.Group", "Protein.Names", "Genes", "Proteotypic")

# The name this limpa's readDIANN() gives the annotation argument, or NA.
limpa_annotation_arg <- function(fn = limpa::readDIANN)
  intersect(c("annotation.columns", "extra.columns"), names(formals(fn)))[1]

limpa_annotation_default <- function(fn = limpa::readDIANN) {
  arg <- limpa_annotation_arg(fn)
  if (is.na(arg)) LIMPA_ANNOTATION_DEFAULT else eval(formals(fn)[[arg]])
}

limpa_version <- function()
  tryCatch(as.character(utils::packageVersion("limpa")), error = function(e) NA_character_)

# The name this limpa's readDIANN() gives the quantity column, or NA: limpa 1.4.x
# intensity.column, 1.2.x qty.column (RELEASE_3_22 R/readDIANN.R; both default to
# "Precursor.Normalised"). The non-normalised DE (run_de.R --quantities raw) reads DIA-NN's
# Precursor.Quantity through it; limpa itself applies no between-run normalisation.
LIMPA_INTENSITY_DEFAULT <- "Precursor.Normalised"
limpa_intensity_arg <- function(fn = limpa::readDIANN)
  intersect(c("intensity.column", "qty.column"), names(formals(fn)))[1]

# readDIANN() with the annotation columns `annotation`, whatever this limpa calls the
# argument. Returns list(dat, argument, path): `argument` is the name used (NA for the join)
# and `path` one line for de_provenance.json. `intensity`: the quantity column, passed under
# this limpa's name for it; a readDIANN() with no such argument can read only its default, so
# any other column is refused rather than silently replaced by the default.
read_diann_annotated <- function(file, format, q.cutoffs, q.columns, annotation,
                                 fn = limpa::readDIANN, intensity = LIMPA_INTENSITY_DEFAULT) {
  arg <- limpa_annotation_arg(fn)
  args <- list(file, format = format, q.cutoffs = q.cutoffs, q.columns = q.columns)
  if (!identical(intensity, LIMPA_INTENSITY_DEFAULT)) {
    iarg <- limpa_intensity_arg(fn)
    if (is.na(iarg))
      stop("this limpa's readDIANN() takes no quantity-column argument, so it cannot read ",
           intensity, " (it reads only ", LIMPA_INTENSITY_DEFAULT, "). Install limpa >= 1.2.")
    args[[iarg]] <- intensity
  }
  if (!is.na(arg)) {
    args[[arg]] <- annotation
    return(list(dat = do.call(fn, args), argument = arg,
                path = sprintf("readDIANN(%s = ...)", arg)))
  }
  dat <- do.call(fn, args)
  missing <- setdiff(annotation, names(dat$genes))
  if (length(missing)) {
    avail <- if (identical(format, "parquet")) names(arrow::open_dataset(file)$schema)
             else names(utils::read.delim(file, nrows = 1, check.names = FALSE))
    cols <- intersect(c("Precursor.Id", missing), avail)
    r <- if (identical(format, "parquet"))
      as.data.frame(nanoparquet::read_parquet(file, col_select = cols))
    else data.table::fread(file, select = cols, data.table = FALSE, showProgress = FALSE)
    r <- r[!duplicated(r$Precursor.Id), , drop = FALSE]
    i <- match(rownames(dat$E), r$Precursor.Id)
    for (cc in setdiff(cols, "Precursor.Id")) dat$genes[[cc]] <- r[[cc]][i]
  }
  list(dat = dat, argument = NA_character_,
       path = sprintf("readDIANN() without an annotation argument; %s joined from the report by Precursor.Id",
                      if (length(missing)) paste(missing, collapse = ", ") else "nothing"))
}
