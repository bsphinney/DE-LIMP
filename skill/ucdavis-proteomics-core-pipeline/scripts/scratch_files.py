"""
scratch_files.py -- the ONE rule for files that are scratch: never zipped, catalogued or copied
to a record. session.py (the session zip) and make_report.py (OUTPUT_FILES.md) both use it.

  * every `.cache` directory -- make_podcast.py's per-chunk TTS audio (~48 KB per second of
    speech, ~60 MB for a 20-minute episode), kept only so a render can resume;
  * every `*.part` file -- a write still in progress, or one that died;
  * every `*.stale.pdf` -- a report PDF that no longer matches its HTML (html_to_pdf.print_report
    renames it rather than leave it looking current). It must never reach a collaborator;
    finalize deletes it once a current PDF exists (session.report_pdf_step);
  * every `*.stale.html` -- the same for the shareable report with the audio built in
    (make_podcast.ensure_share: Analysis_Report_with_audio.stale.html);
  * podcast.wav when podcast.m4a sits beside it (render --keep-wav): the .m4a is the episode.

They stay on disk where they are. Stdlib only.
"""
import os

SCRATCH_DIRS = frozenset({".cache"})
LABEL = ("scratch: .cache folders, *.part files, *.stale.pdf/.html, podcast.wav beside "
         "podcast.m4a (kept on disk)")


def is_scratch_dir(name):
    return name in SCRATCH_DIRS


def is_scratch_file(name, siblings=()):
    return (name.endswith((".part", ".stale.pdf", ".stale.html"))
            or (name == "podcast.wav" and "podcast.m4a" in siblings))


def prune(root, dirs, files):
    """One os.walk step: drops scratch folders from `dirs` in place (so the walk skips them) and
    returns (the files to keep, how many scratch files were left out -- everything under a
    dropped folder included)."""
    n = 0
    for d in [d for d in dirs if is_scratch_dir(d)]:
        dirs.remove(d)
        n += sum(len(fs) for _, _, fs in os.walk(os.path.join(root, d)))
    names = set(files)
    kept = [f for f in files if not is_scratch_file(f, names)]
    return kept, n + len(files) - len(kept)
