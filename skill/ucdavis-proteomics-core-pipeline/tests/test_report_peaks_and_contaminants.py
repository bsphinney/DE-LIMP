#!/usr/bin/env python3
"""
Two checks of skill 2.14.0 that both reach what the client reads, in one session.

The contaminant competition check (contaminant_competition.py) is a finding of the audit:
AUDIT.md for staff, and a row of the client's Audit & caveats table (from AUDIT.json). The
Skyline peak check (skyline_peaks.py) is a section of the report of its own, from
output/skyline/skyline_provenance.json. They were built on separate branches; here one session
has both, made by the real scripts (contaminant_competition.py, audit_results.py,
make_analysis_html.py --no-pdf): each says its piece once, in its own place, and nothing in the
page says the peak check confirmed anything. Stdlib only, as CI is.
"""
import base64
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import client_text  # noqa: E402
import contaminant_competition as cc  # noqa: E402
import skyline_peaks as sp  # noqa: E402
import test_contaminant_competition as tcc  # noqa: E402  fixtures only (no classes run here)
import test_skyline_peaks as tsp  # noqa: E402  fixtures only (no classes run here)

# a 1x1 PNG: the page reads the peak picture it names
PNG = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII=")


def section(md, title):
    m = re.search(r"^## " + re.escape(title) + r"\n(.*?)(?=^## |\Z)", md, re.S | re.M)
    return m.group(1) if m else None


class OneSessionWithBoth(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._tmp = tempfile.TemporaryDirectory()
        d = cls._tmp.name
        session = os.path.join(d, "session")
        out = os.path.join(session, "output")
        fasta = os.path.join(d, "search.fasta")
        tcc.write(fasta, tcc.FASTA)
        tcc.write(fasta + ".meta.json", json.dumps({"organism": "Homo sapiens", "fasta": fasta}))
        # the search: bovine actin took ACTB's peptides; its last job wrote the record
        search = os.path.join(d, "search_out")
        rep = tcc.write_tsv(os.path.join(search, "report.tsv"), tcc.CLEAN + tcc.COMPETED)
        tcc.write(os.path.join(search, "search_provenance.json"), json.dumps(
            {"engine": "diann", "fasta": fasta, "result": {"report": rep}}))
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "contaminant_competition.py"),
                            "--search-out", search], capture_output=True, text=True, timeout=120)
        assert r.returncode == 0, r.stderr
        # the DE's record, and the audit of it
        tables = os.path.join(out, "tables")
        tcc.write(os.path.join(tables, "de_provenance.json"),
                  json.dumps({"input": rep, "contaminants": tcc.DE_REMOVED}))
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                            "--out", os.path.join(out, "AUDIT.md"), "--de-dir", tables,
                            "--fasta-meta", fasta + ".meta.json", "--search-out", search],
                           capture_output=True, text=True, timeout=180, cwd=d)
        assert r.returncode == 0, r.stderr
        # the peak check's record and its picture
        sky = os.path.join(out, "skyline")
        os.makedirs(sky)
        with open(sp.record_path(sky), "w") as fh:
            json.dump(tsp.OK_REC, fh)
        with open(os.path.join(sky, "peaks_G1.png"), "wb") as fh:
            fh.write(PNG)
        tcc.write(os.path.join(out, "AI_Analysis_Report.md"),
                  "# A report\n\n## Key Findings\n\nThings changed.\n")
        cls.build = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"), "--session", session,
             "--out", os.path.join(out, "Analysis_Report.html"), "--no-pdf"],
            capture_output=True, text=True, timeout=180)
        cls.audit_md = tcc.read(os.path.join(out, "AUDIT.md"))
        cls.audit = json.loads(tcc.read(os.path.join(out, "AUDIT.json")))
        cls.page = cls.twin = ""
        if cls.build.returncode == 0:
            cls.page = tcc.read(os.path.join(out, "Analysis_Report.html"))
            cls.twin = tcc.read(os.path.join(out, "Analysis_Report.md"))

    @classmethod
    def tearDownClass(cls):
        cls._tmp.cleanup()

    def setUp(self):
        self.assertEqual(self.build.returncode, 0, self.build.stderr)

    def test_audit_md_carries_the_contaminant_finding_for_staff(self):
        self.assertIn("**contaminant_competition** — COMPETED OUT", self.audit_md)
        self.assertIn("--drop-contaminant Cont_P60712", self.audit_md)
        (f,) = [x for x in self.audit["findings"] if x["check"] == "contaminant_competition"]
        self.assertEqual(f["status"], "WARN")
        self.assertIn("ACTB is missing from the results", f["plain"])

    def test_the_client_table_has_the_contaminant_row_once(self):
        sec = section(self.twin, "Audit & caveats")
        self.assertIsNotNone(sec, self.twin)
        lines = [re.sub(r"^>\s?", "", ln) for ln in sec.splitlines()]   # the twin's callout
        rows = [ln for ln in lines if ln.startswith("| Sample proteins taken by the contaminant "
                                                    "list |")]
        self.assertEqual(len(rows), 1, sec)
        self.assertIn("ACTB is missing from the results", rows[0])
        self.assertNotIn("Skyline", sec)                      # the peak check is not an audit row
        for t in ("Cont_", "COMPETED", "--drop", "fetch_fasta"):
            self.assertNotIn(t, sec)

    def test_the_peak_check_has_its_own_section_once(self):
        self.assertEqual(self.twin.count(f"## {sp.SECTION_TITLE}\n"), 1, self.twin)
        self.assertEqual(self.twin.count("## Audit & caveats\n"), 1, self.twin)
        sec = section(self.twin, sp.SECTION_TITLE)
        self.assertIn("These pictures are a visual check of DIA-NN's peaks.", sec)
        self.assertIn("not independent confirmation", sec)
        self.assertIn("skyline/peaks_G1.png", sec)
        self.assertNotIn("contaminant list", sec)
        self.assertIn(sp.SECTION_TITLE, self.page)
        self.assertIn("Sample proteins taken by the contaminant list", self.page)

    def test_nothing_in_the_report_overclaims_the_peak_check(self):
        self.assertEqual(client_text.peak_check_overclaims(self.twin), [])
        self.assertEqual(client_text.peak_check_overclaims_loose(self.twin), [])
        self.assertNotIn("claims more for the Skyline peak check", self.build.stderr)

    def test_the_record_of_the_search_is_the_one_the_audit_read(self):
        rec = cc.load(os.path.join(self._tmp.name, "search_out"))
        self.assertEqual([g["gene"] for g in rec["competed_out"]], ["ACTB"])
        (f,) = [x for x in self.audit["findings"] if x["check"] == "contaminant_competition"]
        self.assertNotIn("checked_by", f.get("detail") or {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
