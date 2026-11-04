#!/usr/bin/env python3
"""
estimate_params.py --ids-only writes a Sage config with LFQ off (a Core submission, 2026-10-29: the parameters side).

Sage 0.14.7 reads a timsTOF .d directly, but has no MS1 quantification on it: a search that needs
identifications only (no differential expression) can skip the mzML conversion, and its config
says so with `"quant": {"lfq": false}`. run_search.py reads that key (sage_lfq_enabled) to send
the .d straight to Sage and to skip the LFQ check; nothing else carries the choice.

A DIA-NN search always quantifies, so --ids-only is refused for it (exit 2), with nothing left at
--out.

Hermetic: estimate_params.py as a subprocess on a temp dir; no network.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")

WHY = ("identifications only (--ids-only: no differential expression); Sage 0.14.7 reads a "
       "timsTOF .d directly but has no MS1 quantification on it")


def run(d, engine, *extra):
    out = os.path.join(d, "params.cfg" if engine == "diann" else "sage.json")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                        "--engine", engine, "--acquisition", "DDA", "--instrument", "timsTOF HT",
                        "--precursor-mz-range", "100", "1700", "--out", out, *extra],
                       capture_output=True, text=True, timeout=60)
    return p, out


def load(out):
    with open(out) as fh:
        cfg = json.load(fh)
    with open(out + ".rationale.json") as fh:
        return cfg, json.load(fh)


class IdsOnly(unittest.TestCase):
    def test_ids_only_turns_lfq_off_and_says_why(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = run(d, "sage", "--ids-only")
            self.assertEqual(p.returncode, 0, p.stderr)
            cfg, side = load(out)
        self.assertIs(cfg["quant"]["lfq"], False)
        self.assertEqual(side["rationale"]["lfq"], {"value": False, "source": WHY})

    def test_ids_only_refused_for_diann(self):
        with tempfile.TemporaryDirectory() as d:
            stale = os.path.join(d, "params.cfg")
            for f in (stale, stale + ".rationale.json"):       # an earlier run's files
                with open(f, "w") as fh:
                    fh.write("--qvalue 0.01\n")
            p, out = run(d, "diann", "--ids-only")
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("--ids-only is for Sage; a DIA-NN search always quantifies", p.stderr)
            self.assertFalse(os.path.exists(out))
            self.assertFalse(os.path.exists(out + ".rationale.json"))

    def test_default_keeps_lfq_on(self):
        """(pin) Without --ids-only, Sage quantifies, as before."""
        with tempfile.TemporaryDirectory() as d:
            p, out = run(d, "sage")
            self.assertEqual(p.returncode, 0, p.stderr)
            cfg, side = load(out)
        self.assertIs(cfg["quant"]["lfq"], True)
        self.assertIs(side["rationale"]["lfq"]["value"], True)

    def test_a_quant_override_beside_ids_only_is_refused(self):
        """--overrides '{"quant": ...}' would replace the block --ids-only sets, and the rationale
        would describe a config that is not the one written."""
        with tempfile.TemporaryDirectory() as d:
            p, out = run(d, "sage", "--ids-only", "--overrides", '{"quant": {"lfq": true}}')
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("--ids-only", p.stderr)
            self.assertIn("quant", p.stderr)
            self.assertFalse(os.path.exists(out))


if __name__ == "__main__":
    unittest.main()
