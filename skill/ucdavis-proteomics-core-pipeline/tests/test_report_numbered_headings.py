#!/usr/bin/env python3
"""A numbered report heading ("## 9. Data Quality Notes") replaces the tool's raw block (I12).

mstaff, 2026-10-22 (a Core submission): the report's own section was "## 9. Data Quality Notes", but
norm_title() kept the "9." and so did not see that the report already covered sample quality.
The page then also showed SAMPLE_QUALITY.md as it is: the z-score table, run file names and
developer notes reached the client. "10. Audit & caveats" failed the same way for AUDIT.md.

Guards (stdlib only):
  * "## 9. Data Quality Notes" removes the Sample-quality block: no z-score table, no run file
    name from SAMPLE_QUALITY.md;
  * "## 10. Audit & caveats" removes the raw audit block;
  * norm_title("9) Audit") == "audit";
  * an un-numbered heading is still matched (pin).
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import make_analysis_html as mah  # noqa: E402

RUN_FILE = "20251012_PROT9798_bovine_S07_rep2.d"
SAMPLE_QUALITY = ("# Sample quality\n\n## Panel abundance\n\n"
                  "| sample | group | z (panel abundance) |\n|---|---|---|\n"
                  f"| {RUN_FILE} | A | 3.4 |\n\nDEV NOTE: keratin panel v2 thresholds.\n")
AUDIT = "# Results audit\n\n- ⚠️ **B3** — RAW-TOOL-AUDIT-LINE not real biology\n"


class NumberedHeadings(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.session = self._td.name
        self.out = os.path.join(self.session, "output")
        os.makedirs(os.path.join(self.out, "figures"))
        self.write("SAMPLE_QUALITY.md", SAMPLE_QUALITY)
        self.write("AUDIT.md", AUDIT)

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, text):
        with open(os.path.join(self.out, name), "w", encoding="utf-8") as fh:
            fh.write(text)

    def absent(self, needle, page):
        self.assertFalse(needle in page, f"{needle!r} is on the page")

    def present(self, needle, page):
        self.assertTrue(needle in page, f"{needle!r} is not on the page")

    def build(self, report):
        self.write("AI_Analysis_Report.md", report)
        html_path = os.path.join(self.out, "Analysis_Report.html")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", self.session, "--out", html_path, "--no-pdf"],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        json.loads(r.stdout)
        with open(html_path, encoding="utf-8") as fh:
            return fh.read()

    def test_numbered_data_quality_heading_removes_the_sample_quality_block(self):
        page = self.build("# Study\n\n## 1. Overview\n\nText.\n\n"
                          "## 9. Data Quality Notes\n\nOWN-QUALITY-NOTES here.\n")
        self.present("OWN-QUALITY-NOTES", page)
        self.absent('id="sample-quality-notes"', page)
        self.absent("z (panel abundance)", page)
        self.absent(RUN_FILE, page)
        self.absent("DEV NOTE", page)

    def test_numbered_audit_heading_removes_the_raw_audit_block(self):
        page = self.build("# Study\n\n## 1. Overview\n\nText.\n\n"
                          "## 10. Audit & caveats\n\n- OWN-AUDIT line\n")
        self.present("OWN-AUDIT", page)
        self.absent('id="audit-caveats"', page)      # the tool's block, by its own title
        self.absent("<td>B3</td>", page)              # its check code
        self.absent("RAW-TOOL-AUDIT-LINE", page)

    def test_norm_title_strips_a_leading_section_number(self):
        self.assertEqual(mah.norm_title("9) Audit"), "audit")
        self.assertEqual(mah.norm_title("9. Data Quality Notes"), "data quality notes")
        self.assertEqual(mah.norm_title("10.2 Audit &amp; caveats"), "audit & caveats")
        self.assertEqual(mah.norm_title("IV. Audit"), "audit")
        self.assertEqual(mah.norm_title("iv) Audit"), "audit")

    def test_a_first_word_made_of_roman_letters_is_kept(self):
        """L5a review R9: a Roman numeral needs "." or ")" after it, so "CV analysis" and
        "LC MS" keep their first word."""
        for t in ("CV analysis", "LC MS", "CL test", "Mix design", "XL crosslinks"):
            self.assertEqual(mah.norm_title(t), t.lower(), t)
        self.assertEqual(mah.norm_title("9 Audit"), "audit")              # a bare number still

    def test_unnumbered_heading_still_matches(self):
        """(pin) an un-numbered heading was already matched; it still is."""
        self.assertEqual(mah.norm_title("Data Quality Notes"), "data quality notes")
        page = self.build("# Study\n\n## Overview\n\nText.\n\n"
                          "## Data Quality Notes\n\nOwn notes.\n\n## Audit & caveats\n\n- own\n")
        self.absent("z (panel abundance)", page)
        self.absent("<td>B3</td>", page)

    def test_numbered_heading_keeps_its_callout(self):
        """severity() reads norm_title too: a numbered audit section is still a callout."""
        self.assertEqual(mah.severity("10. Audit & caveats", "- fine"), "warning")
        self.assertEqual(mah.severity("9. Data Quality Notes", "- FAIL x"), "critical")


if __name__ == "__main__":
    unittest.main()
