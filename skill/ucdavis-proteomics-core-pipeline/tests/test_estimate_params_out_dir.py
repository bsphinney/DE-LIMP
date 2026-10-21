#!/usr/bin/env python3
"""estimate_params.py makes --out's folder when it does not exist yet (an eval run found it
crashing with a FileNotFoundError traceback on a new session's wf/)."""
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")


class OutFolder(unittest.TestCase):
    def test_a_missing_out_folder_is_made(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "input", "wf", "params.cfg")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                                "--engine", "diann", "--acquisition", "DIA", "--instrument", "timsTOF HT",
                                "--precursor-mz-range", "300", "1200", "--out", out],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertNotIn("Traceback", p.stderr)
            self.assertTrue(os.path.isfile(out))
            self.assertTrue(os.path.isfile(out + ".rationale.json"))


if __name__ == "__main__":
    unittest.main()
