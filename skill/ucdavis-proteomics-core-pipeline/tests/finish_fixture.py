"""What a finished session holds before `session.py finalize`, for the tests that finalize one
(2.14.2: finalize refuses without the skill's own report and the reproducibility bundle; the
gates' own tests are test_finish_gates.py). Not a test module."""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

import report_files  # noqa: E402

# For a test where a missing finish item is the point: finalize goes ahead under this reason,
# which MANIFEST.txt records. A test of something else finishes the session first (finished()),
# so a regression in finalize's gate shows there too.
FIXTURE_OVERRIDE = "a test fixture: this test is about another part of finalize"


def skill_report(output_dir, text="<html><body><p>the skill's report</p></body></html>"):
    """output/Analysis_Report.html as make_analysis_html.py leaves it: the page and its record."""
    path = os.path.join(output_dir, "Analysis_Report.html")
    os.makedirs(output_dir, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    report_files.mark_skill_report(path, "make_analysis_html.py")
    return path


def bundle(output_dir):
    """output/reproducibility/run_manifest.json, as provenance.py (step 10) leaves it."""
    path = os.path.join(output_dir, "reproducibility", "run_manifest.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump({"files": []}, fh)
    return path


def finished(output_dir, report=True, methods=False):
    """What finalize's gate asks for: the skill's report (the one there, marked, or a new one;
    none with report=False) and the bundle (unless the test wrote its own). With `methods`, an
    output/methods.md too (step 9d), for a session finalize cannot write Methods for (no raw-file
    list): its stamp needs one."""
    if methods and not os.path.isfile(os.path.join(output_dir, "methods.md")):
        os.makedirs(output_dir, exist_ok=True)
        with open(os.path.join(output_dir, "methods.md"), "w", encoding="utf-8") as fh:
            fh.write("# Materials and Methods\n\nWritten at step 9d.\n")
    if report and not os.path.isfile(os.path.join(output_dir, "Analysis_Report.html")):
        skill_report(output_dir)
    elif report:
        report_files.mark_skill_report(os.path.join(output_dir, "Analysis_Report.html"),
                                       "make_analysis_html.py")
    if not os.path.isfile(os.path.join(output_dir, "reproducibility", "run_manifest.json")):
        bundle(output_dir)


def fake_pipeline_home(root):
    """<root>/pipeline_home: a skill home ($PROTEOMICS_PIPELINE_HOME) whose setup.json names a
    fake python and Rscript, so a provenance.py that finalize runs (it builds a missing bundle)
    lists no real conda env and starts no R. -> the folder, for the env."""
    home = os.path.join(root, "pipeline_home")
    bin_dir = os.path.join(home, "bin")
    os.makedirs(bin_dir, exist_ok=True)
    progs = {}
    for name, says in (("python3", "fake-package==1.0"), ("Rscript", "R version 4.5.1 (fake)")):
        path = os.path.join(bin_dir, name)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(f"#!/bin/sh\necho '{says}'\n")
        os.chmod(path, 0o755)
        progs[name] = path
    with open(os.path.join(home, "setup.json"), "w", encoding="utf-8") as fh:
        json.dump({"python": progs["python3"], "rscript": progs["Rscript"]}, fh)
    return home


def stamp(output_dir, **fields):
    """output/FINALIZED.json as session.py finalize leaves it (the fields its readers use), with
    the real sha256 of each stamped file that is there (report_files.STAMP_OUTPUTS), so a file
    changed after it makes it stale (report_files.stamp_drift), as it would a real one."""
    sd = report_files.session_of(output_dir)
    outputs = {rel: (report_files._sha256(os.path.join(sd, rel))
                     if os.path.isfile(os.path.join(sd, rel)) else "[SKIPPED] not in the session")
               for rel, _part, _required in report_files.STAMP_OUTPUTS}
    path = os.path.join(output_dir, report_files.FINALIZED_RECORD)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(dict({"schema": report_files.FINALIZED_SCHEMA,
                        "finalized_at": "2026-11-01T00:00:00+00:00", "skill_version": "2.14.2",
                        "outputs": outputs}, **fields), fh)
    return path
