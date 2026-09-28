# =============================================================================
# skill_version.R -- R mirror of skill_version.py, the ONE reader of this skill's version
# (.claude-plugin/plugin.json beside scripts/). run_de.R names the version with it in the
# DE-LIMP session ("App version: ..."); it used to read plugin.json with its own regex and its
# own fallback text. Kept equal to the Python original, fixture for fixture, by
# tests/test_skill_version.py. When plugin.json is not there -- a scripts/ folder copied up on
# its own -- the version is the tagged UNKNOWN, never a guessed one (CLAUDE.md rule 2).
# =============================================================================

# "(unknown -- plugin.json not found)" with skill_version.py's em dash, written as an escape so
# the file reads the same under any locale.
SKILL_VERSION_UNKNOWN <- "(unknown — plugin.json not found)"

# The version in <here>/../.claude-plugin/plugin.json, or SKILL_VERSION_UNKNOWN when it cannot
# be read: a missing file, no "version", a blank one, not JSON, not a JSON object.
skill_version <- function(here, use_jsonlite = requireNamespace("jsonlite", quietly = TRUE)) {
  f <- file.path(here, "..", ".claude-plugin", "plugin.json")
  v <- tryCatch({
    if (use_jsonlite) {
      m <- jsonlite::fromJSON(f, simplifyVector = FALSE)
      if (is.list(m) && !is.null(names(m))) m[["version"]] else NULL
    } else {
      # no jsonlite: the top-level "version": "<string>" of a JSON object
      txt <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
      if (!grepl("^[[:space:]]*\\{", txt)) NULL
      else {
        m <- regmatches(txt, regexpr('"version"[[:space:]]*:[[:space:]]*"[^"]*"', txt))
        if (length(m)) sub('^.*:[[:space:]]*"([^"]*)"$', "\\1", m) else NULL
      }
    }
  }, error = function(e) NULL, warning = function(w) NULL)
  if (is.character(v) && length(v) == 1L && nzchar(trimws(v))) trimws(v) else SKILL_VERSION_UNKNOWN
}

# How a version reads after the skill's name: "v2.8.0", or the UNKNOWN tag as it is.
skill_label <- function(version)
  if (identical(version, SKILL_VERSION_UNKNOWN)) version else paste0("v", version)
