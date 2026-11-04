"""
scratch_files.py -- the ONE rule for files that are scratch: never zipped, catalogued or copied
to a record. session.py (the session zip) and make_report.py (OUTPUT_FILES.md) both use it.

  * every `.cache` directory -- make_podcast.py's per-chunk TTS audio (~48 KB per second of
    speech, ~60 MB for a 20-minute episode), kept only so a render can resume;
  * every `*.part` and `*.partial` file, and every `.partial` folder -- a write still in
    progress, or one that died (run_search.py converts a Sage search's raw files into
    <out>/mzml/.partial/ and moves each mzML out only once it is complete);
  * every `*.stale.pdf` -- a report PDF that no longer matches its HTML (html_to_pdf.print_report
    renames it rather than leave it looking current). It must never reach a collaborator;
    finalize deletes it once a current PDF exists (session.report_pdf_step);
  * every `*.stale.html` -- the same for the shareable report with the audio built in
    (make_podcast.ensure_share: Analysis_Report_with_audio.stale.html);
  * podcast.wav when podcast.m4a sits beside it (render --keep-wav): the .m4a is the episode;
  * every `skyline_work` directory -- the peak check's working folder (skyline_peaks.py): the
    Skyline document's caches, the per-file imports, the job scripts and logs. What a reader
    needs from it is in output/skyline/ (the pictures, the tables and the .sky.zip);
  * every `raw_stage` directory -- a FragPipe DDA search's staged inputs (run_search.py
    run_fragpipe: a copy of each Thermo .raw, a link to each .d / mzML) and what MSFragger writes
    beside them, a second copy of the raw data (~1 GB a run); and every `diatracer_stage`
    directory, the DIA route's links to the .d and diaTracer's pseudo-MS/MS mzML (L11 N1, N4);
  * MSFragger's `*_uncalibrated.mzML`, `*_uncalibrated.mzBIN` and `*.mzBIN` -- what it writes
    next to each input, wherever that was;
  * every `input_links` directory -- a DIA-NN search's links to its raw runs under names DIA-NN
    can read (diann_parallel.INPUT_LINKS_DIR, --rename-dashed). An interim delivery walks a search
    folder FOLLOWING links, and copied every linked .d a second time (C2 security L7).

They stay on disk where they are. Stdlib only.
"""
import os

# The folder run_search.py stages a FragPipe DDA search's inputs in, inside the search folder
FRAGPIPE_STAGE_DIR = "raw_stage"
DIATRACER_STAGE_DIR = "diatracer_stage"
INPUT_LINKS_DIR = "input_links"        # diann_parallel.INPUT_LINKS_DIR (a test ties the two)
SCRATCH_DIRS = frozenset({".cache", "skyline_work", ".partial", FRAGPIPE_STAGE_DIR,
                          DIATRACER_STAGE_DIR, INPUT_LINKS_DIR})
# Why a scratch folder stays out, where the generic reason would mislead (a delivery's SKIPPED)
DIR_REASONS = {
    FRAGPIPE_STAGE_DIR: "excluded: the search's staged copies of the raw files and MSFragger's "
                        "uncalibrated mzML/mzBIN beside them, a second copy of the raw data (the "
                        "raw files are in raw/)",
    DIATRACER_STAGE_DIR: "excluded: the search's links to the raw .d files and diaTracer's "
                         "pseudo-MS/MS mzML, a working copy of the raw data (the raw files are "
                         "in raw/)",
    INPUT_LINKS_DIR: "excluded: the search's links to the raw files under names DIA-NN can read, "
                     "not data of their own (the raw files are in raw/)"}
# MSFragger writes these beside each input (lower-cased for the test)
MSFRAGGER_SCRATCH_SUFFIXES = ("_uncalibrated.mzml", "_uncalibrated.mzbin", ".mzbin")
LABEL = ("scratch: .cache, skyline_work, .partial, raw_stage, diatracer_stage and input_links "
         "folders, "
         "*.part/*.partial files, *.stale.pdf/.html, MSFragger's *_uncalibrated.mzML/*.mzBIN, "
         "podcast.wav beside podcast.m4a (kept on disk)")


def is_scratch_dir(name):
    return name in SCRATCH_DIRS


def dir_reason(name, default):
    """Why scratch folder `name` stays out: its own reason (DIR_REASONS), else `default`."""
    return DIR_REASONS.get(name, default)


def is_scratch_file(name, siblings=()):
    return (name.endswith((".part", ".partial", ".stale.pdf", ".stale.html"))
            or name.lower().endswith(MSFRAGGER_SCRATCH_SUFFIXES)
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
