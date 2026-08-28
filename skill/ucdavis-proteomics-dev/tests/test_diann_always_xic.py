#!/usr/bin/env python3
"""
Every DIA-NN search extracts chromatograms.

XICs are what let a person LOOK at an identification instead of trusting a q-value, and they
are the input to FRAN's XIC lane. `estimate_params.py` writes `--xic` into the cfg it
generates, but a cfg can also come from a workflow bundle or straight from the user — and if
that silently produced a search with no chromatograms, nobody would find out until they went
looking for a trace that was never extracted.
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import run_search  # noqa: E402


class GeneratedCfgTests(unittest.TestCase):
    def test_estimate_params_puts_xic_in_every_diann_cfg(self):
        """The generated cfg is the normal source, and it must carry XIC for every instrument —
        including the ones with no ion mobility, where DIA-NN simply ignores --mobilograms."""
        for instr in ("Orbitrap Astral", "timsTOF Pro", "Orbitrap Exploris 480"):
            with tempfile.TemporaryDirectory() as d:
                cfg = os.path.join(d, "diann.cfg")
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                                    "--engine", "diann", "--acquisition", "DIA",
                                    "--instrument", instr, "--out", cfg],
                                   capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, p.stderr)
                txt = open(cfg).read()
                self.assertIn("--xic", txt, instr)
                # --xic ALONE allocates mobilogram parquets and leaves them full of zeros,
                # silently and at plausible file size. The two flags travel together.
                self.assertIn("--mobilograms", txt, instr)


class EnforcementTests(unittest.TestCase):
    def test_a_cfg_without_xic_is_augmented_not_run_as_is(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "bundle.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n")
            out = os.path.join(d, "search_out")
            used = run_search.ensure_xic(cfg, out)
            self.assertNotEqual(used, cfg)
            txt = open(used).read()
            self.assertIn("--xic", txt)
            self.assertIn("--mobilograms", txt)
            self.assertIn("--mass-acc 15", txt)             # the original settings survive
            # The user's file is never edited -- the cfg that RAN is written beside the results.
            self.assertNotIn("--xic", open(cfg).read())
            self.assertEqual(os.path.dirname(used), out)

    def test_a_cfg_that_already_has_xic_is_left_exactly_alone(self):
        """Including its window: a cfg asking for --xic 30 must not be silently reset to 10."""
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--xic 30\n--mobilograms\n")
            self.assertEqual(run_search.ensure_xic(cfg, os.path.join(d, "o")), cfg)

    def test_mobilograms_is_added_when_only_xic_is_present(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n")
            txt = open(run_search.ensure_xic(cfg, os.path.join(d, "o"))).read()
            self.assertEqual(txt.count("--mobilograms"), 1)

    def test_an_unreadable_cfg_does_not_raise(self):
        """A missing params file is the caller's error to report; failing here would replace a
        clear "no such cfg" with a traceback out of an unrelated helper."""
        self.assertEqual(run_search.ensure_xic("/no/such/file.cfg", "/tmp/x"),
                         "/no/such/file.cfg")


class ParallelChainTests(unittest.TestCase):
    """The 5-step chain strips --xic from the shared flags and re-adds it to step 4 alone."""

    def test_xic_flag_is_recovered_from_the_cfg_with_its_window(self):
        import diann_parallel
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--xic 10\n--mobilograms\n")
            flag = diann_parallel.xic_flag(cfg)
            self.assertIn("--xic 10", flag)
            self.assertIn("--mobilograms", flag)

    def test_no_xic_in_the_cfg_means_no_flag(self):
        import diann_parallel
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n")
            self.assertEqual(diann_parallel.xic_flag(cfg), "")


if __name__ == "__main__":
    unittest.main(verbosity=2)
