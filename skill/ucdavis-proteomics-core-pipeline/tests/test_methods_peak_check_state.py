#!/usr/bin/env python3
"""The Methods describe the peak check (Skyline) only when the delivery can carry it (2.14.2,
L3 review).

make_methods.py wrote "Peak check (Skyline)" whenever the record said status "ok" -- also when
the check was STALE (the DE was run again after it). Lane L3 holds output/skyline/ back from a
delivery unless the check is finished and current, so a delivered methods.md could describe a
check the delivery does not hold. The Methods now read the state the delivery reads
(skyline_peaks.record_state, one definition): only "ok" gets the paragraph; a stale or unreadable
check gets one tagged note, a skipped one nothing.
stdlib only, no network.
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
sys.path.insert(0, HERE)

import skyline_peaks as sp                       # noqa: E402
from test_skyline_peaks import OK_REC            # noqa: E402

PARAGRAPH_WORDS = "not an independent identification"


class PeakCheckState(unittest.TestCase):

    def methods(self, **rec_extra):
        with tempfile.TemporaryDirectory() as td:
            tables, sky = (os.path.join(td, "output", "tables"),
                           os.path.join(td, "output", "skyline"))
            os.makedirs(tables)
            os.makedirs(sky)
            de = os.path.join(tables, "de_provenance.json")
            with open(de, "w") as fh:
                json.dump({"adjp": 0.05, "method": "dpc"}, fh)
            rec = dict(OK_REC, methods_paragraph=sp.methods_paragraph(dict(OK_REC)))
            rec.update({k: (sp.de_fingerprint(de) if v == "CURRENT" else v)
                        for k, v in rec_extra.items()})
            with open(os.path.join(sky, sp.RECORD), "w") as fh:
                json.dump(rec, fh)
            raw = os.path.join(td, "FL_sample01.raw")
            open(raw, "w").close()
            out = os.path.join(td, "methods.md")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                                "--raw", raw, "--de-dir", tables, "--out", out],
                               capture_output=True, text=True, timeout=300)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out, encoding="utf-8") as fh:
                return fh.read()

    def test_an_ok_but_stale_check_is_not_described(self):
        text = self.methods(de_record_sha256="0" * 64)        # its DE was re-run since
        self.assertNotIn(f"## {sp.SECTION_TITLE}", text)
        self.assertNotIn(PARAGRAPH_WORDS, text)
        self.assertIn("Peak check (Skyline) not described", text)
        self.assertIn("run again after this peak check", text)
        self.assertIn("— confirm]", text)

    def test_a_current_check_is_described(self):                        # (pin)
        text = self.methods(de_record_sha256="CURRENT")
        self.assertIn(f"## {sp.SECTION_TITLE}", text)
        self.assertIn(PARAGRAPH_WORDS, text)

    def test_a_skipped_check_is_not_mentioned(self):                    # (pin)
        text = self.methods(status="skipped", reason="turned off")
        self.assertNotIn("Peak check", text)

    def test_a_failed_check_is_not_described(self):                    # (pin)
        text = self.methods(status="failed")                  # record_state: "unreadable"
        self.assertNotIn(f"## {sp.SECTION_TITLE}", text)
        self.assertNotIn(PARAGRAPH_WORDS, text)


if __name__ == "__main__":
    unittest.main()
