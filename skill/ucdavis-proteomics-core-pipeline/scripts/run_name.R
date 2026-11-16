# run_name(): the name DIA-NN reports a run under -- the raw file's name, no folder, no extension.
# The R twin of check_report_runs.py run_name() (keep the extension list identical): both read a
# File.Name of report.stats.tsv or a raw-file path, with / or \ separators (DIA-NN on Windows).
RUN_NAME_EXTS <- c(".raw", ".d", ".mzml", ".wiff", ".wiff2", ".dia")
run_name <- function(path) {
  base <- basename(sub("/+$", "", gsub("\\", "/", as.character(path), fixed = TRUE)))
  ext <- tolower(sub("^.*?(\\.[^.]*)?$", "\\1", base, perl = TRUE))    # os.path.splitext: the last dot, not a leading one
  has <- nzchar(ext) & ext %in% RUN_NAME_EXTS & nchar(base) > nchar(ext)
  ifelse(has, substr(base, 1L, nchar(base) - nchar(ext)), base)
}
