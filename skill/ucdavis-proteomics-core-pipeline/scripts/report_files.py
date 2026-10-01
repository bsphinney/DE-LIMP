"""
report_files.py -- the names of report files that more than one script needs, defined once.
No imports: make_report.py (OUTPUT_FILES.md) reads it without pulling in make_podcast.py and
its dependencies, so a partial copy of scripts/ still catalogues.

  * SHARE_NAME  -- the report with the podcast's audio and transcript built in, one file to
    send (make_podcast.py share writes it; core_submission.py deliver ships it);
  * SHARE_STALE -- the name an out-of-date copy is renamed to, so it is never sent
    (scratch_files.py leaves every *.stale.html out).
"""
SHARE_NAME = "Analysis_Report_with_audio.html"
SHARE_STALE = "Analysis_Report_with_audio.stale.html"
