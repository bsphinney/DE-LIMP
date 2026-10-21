#!/usr/bin/env python3
"""finalize writes output/OUTPUT_FILES.md itself when SKILL.md step 11 was not run.

Step 11 (make_report.py, the file catalog) was a command the agent had to run by hand before
finalize. When it was skipped, finalize only noted "[INFO] ... step 11 was not run" and the
delivered folder had no catalog. Now finalize builds it from the session's own folders -- step
11's inputs, plus the Methods and deposit package finalize has written by then -- before it adds
README, AGENTS.md and MANIFEST.txt, and records the outcome in MANIFEST.txt:
  * no catalog -> written, [OK] OUTPUT_FILES.md -- written at finalize (step 11 had not run)
  * a catalog step 11 wrote -> never rewritten, only the finalize block appended
  * make_report fails -> [SKIPPED] OUTPUT_FILES.md -- <reason>, and finalize still completes
stdlib only, no network, no browser (no Analysis_Report.html newer than its PDF, so
report_pdf_step never prints).
"""
import io
import json
import os
import re
import shutil
import sys
import tempfile
import time
import unittest
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

from job_env import job_env                                          # noqa: E402
from test_deposit_package import dia_session, finalize, manifest_lines, write  # noqa: E402
import make_report                                                   # noqa: E402
import session                                                       # noqa: E402

OK_NOTE = "written at finalize (step 11 had not run)"


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def catalog_lines(p):
    return [ln for ln in manifest_lines(p) if re.match(r"\[\w+\]\s+OUTPUT_FILES\.md\s+--", ln)]


def add_report_files(p):
    """What steps 8-9 leave in output/: the audit, a figure, the report and its PDF -- the PDF
    newer than the HTML, so finalize's report_pdf_step keeps it and never starts a browser."""
    out = p["output_dir"]
    write(os.path.join(out, "AUDIT.md"), "# Audit\n")
    write(os.path.join(out, "AI_Analysis_Report.md"), "# Report\n")
    write(os.path.join(p["figures_dir"], "volcano_Treated_vs_Control.png"), "png")
    html, pdf = (os.path.join(out, "Analysis_Report." + x) for x in ("html", "pdf"))
    write(html, "<html></html>")
    write(pdf, "%PDF-1.4")
    now = time.time()
    os.utime(html, (now - 100, now - 100))
    os.utime(pdf, (now, now))
    # raw archives prepare_upload.sbatch stages: never in the catalog, as never in the zip
    write(os.path.join(p["deposit_dir"], "upload_staging", "HeLa_ctrl_01.d.tar"), "tar")


class FinalizeWritesTheCatalog(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        self.p = dia_session(tmp)
        add_report_files(self.p)

    def test_no_catalog_finalize_writes_it_and_says_so(self):
        p = self.p
        self.assertFalse(os.path.exists(p["output_files_md"]))
        r = finalize(p["session_dir"], "--zip")
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        res = json.loads(r.stdout)
        self.assertTrue(os.path.isfile(p["output_files_md"]))
        self.assertFalse(os.path.exists(p["output_files_md"] + ".part"))
        md = read(p["output_files_md"])
        self.assertTrue(md.startswith("# Output files — what each one is"))
        # the session's own files, paths relative to the session folder (step 11's inputs)
        for rel in ("output/search/report.parquet", "output/tables/DE_Treated_vs_Control.csv",
                    "output/tables/de_provenance.json", "input/conditions.csv",
                    "input/search.fasta", "input/wf/params.cfg",
                    "input/wf/workflow.manifest.json",
                    "output/figures/volcano_Treated_vs_Control.png", "output/AUDIT.md",
                    "output/AI_Analysis_Report.md", "output/Analysis_Report.html",
                    "output/Analysis_Report.pdf",
                    # written by finalize before the catalog: the Methods, the deposit package
                    "output/methods.md", "output/DATA_SUBMISSION/sdrf.tsv",
                    "output/DATA_SUBMISSION/HOW_TO_SUBMIT.md"):
            self.assertIn(f"`{rel}`", md, rel)
        self.assertIn("output/search/quant_step2/*.quant", md)      # one line per kind
        self.assertNotIn("upload_staging", md)
        # README, AGENTS.md, the raw list: the finalize block, once each (MANIFEST.txt is
        # written after the block, so a first finalize cannot list it)
        self.assertEqual(md.count("written at finalize: BEGIN"), 1)
        block = md.partition("written at finalize: BEGIN")[2]
        for rel in ("README.html", "README.md", "AGENTS.md", "input/raw_files.txt"):
            self.assertEqual(md.count(f"`{rel}`"), 1, rel)
            self.assertIn(f"`{rel}`", block, rel)
        (line,) = catalog_lines(p)
        self.assertTrue(line.startswith("[OK]"), line)
        self.assertIn(f"-- {OK_NOTE}; ", line)
        (added,) = [ln for ln in manifest_lines(p)
                    if "OUTPUT_FILES.md: files written at finalize" in ln]
        self.assertTrue(added.startswith("[OK]"), added)
        # in the zip, and the catalog's [OK] line did not count as skipped
        base = os.path.basename(p["session_dir"])
        with zipfile.ZipFile(res["zip"]) as z:
            self.assertIn(f"{base}/output/OUTPUT_FILES.md", z.namelist())
            self.assertIn(OK_NOTE, z.read(f"{base}/MANIFEST.txt").decode())

    def test_a_catalog_step_11_wrote_is_kept_and_only_appended_to(self):
        p = self.p
        original = "# Output files — what each one is\n\nWritten by step 11, by hand.\n"
        write(p["output_files_md"], original)
        for _ in range(2):                                        # re-finalize: one block
            r = finalize(p["session_dir"])
            self.assertEqual(r.returncode, 0, r.stderr[-2000:])
            md = read(p["output_files_md"])
            self.assertTrue(md.startswith(original.rstrip("\n")), md[:200])
            self.assertEqual(md.count("written at finalize: BEGIN"), 1)
            self.assertIn("`README.html`", md)
            self.assertNotIn("`output/search/report.parquet`", md)   # not rebuilt
            self.assertEqual(catalog_lines(p), [])                   # nothing written: no line


class MakeReportFails(unittest.TestCase):
    def run_finalize(self, p):
        args = type("A", (), dict(dir=p["session_dir"], zip=True, no_deposit=True,
                                  no_notify=True, reanalysis_of=""))()
        out = io.StringIO()
        with mock.patch.dict(os.environ, job_env(os.path.dirname(p["session_dir"])),
                             clear=True), redirect_stdout(out), redirect_stderr(io.StringIO()), \
                mock.patch.object(make_report, "write_catalog",
                                  side_effect=OSError("No space left on device")):
            session.do_finalize(args)
        return json.loads(out.getvalue())

    def test_skipped_with_the_reason_and_finalize_completes(self):
        with tempfile.TemporaryDirectory() as d:
            p = dia_session(d)
            res = self.run_finalize(p)
            self.assertFalse(os.path.exists(p["output_files_md"]))
            (line,) = catalog_lines(p)
            self.assertTrue(line.startswith("[SKIPPED]"), line)
            self.assertIn("not written at finalize (step 11 had not run): OSError: No space "
                          "left on device", line)
            (added,) = [ln for ln in manifest_lines(p)
                        if "OUTPUT_FILES.md: files written at finalize" in ln]
            self.assertTrue(added.startswith("[INFO]"), added)
            self.assertGreaterEqual(res["skipped"], 1)
            # the rest of finalize still ran: README, AGENTS.md, the zip and its MANIFEST.txt
            for n in ("README.md", "README.html", "AGENTS.md"):
                self.assertTrue(os.path.isfile(os.path.join(p["session_dir"], n)), n)
            self.assertEqual(res["zip_docs"], "added")
            self.assertEqual(res["zip_manifest"], "added")

    def test_a_failed_write_leaves_no_half_catalog(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "OUTPUT_FILES.md")
            write(os.path.join(d, "tables", "DE_a_vs_b.csv"), "x")
            with mock.patch.object(make_report.os, "replace", side_effect=OSError("boom")):
                with self.assertRaises(OSError):
                    make_report.write_catalog(out, de_dir=os.path.join(d, "tables"), root=d)
            self.assertEqual(sorted(os.listdir(d)), ["tables"])
            got = make_report.write_catalog(out, de_dir=os.path.join(d, "tables"), root=d)
            self.assertEqual(got["n_files"], 1)
            self.assertIn("`tables/DE_a_vs_b.csv`", read(out))


if __name__ == "__main__":
    unittest.main()
