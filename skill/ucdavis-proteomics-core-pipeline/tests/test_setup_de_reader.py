#!/usr/bin/env python3
"""setup.json's ready_for.de says whether the DE can read report.parquet with the limpa the env has.

HIVE, 2026-10-16: `setup.sh --check` (skill 2.13.1, in an srun job) wrote ready_for.de false for
the account's env -- limpa 1.2.5 (< 1.4.0) and no loadable nanoparquet, "so no DE can run" --
while every DE that week ran on that env and completed (de_provenance.json: limpa 1.2.5,
readDIANN(extra.columns = ...)). Both were half right. limpa 1.2.x reads report.parquet with
arrow (readDIANN: arrow::read_parquet, github.com/bioc/limpa RELEASE_3_22); only limpa >= 1.4.0
(RELEASE_3_23) reads it with nanoparquet. So that env runs the DE -- on an unpinned limpa, and
run_de.R's contaminant share, which reads Precursor.Quantity with nanoparquet, fell back to
Precursor.Normalised ("there is no package called 'nanoparquet'"; de_provenance.json
contaminants.intensity_column).

ready_for.de now follows the reader: arrow always (run_de.R filters the report with it on both
methods), plus nanoparquet when limpa >= 1.4.0. The nanoparquet note says what is lost without it
on limpa 1.2.x instead of "no DE can run". limpa.ok and the exit status still hold the env to the
pinned stack.
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import test_setup_limpa as limpa_t          # noqa: E402


class DeReadinessFollowsTheReader(limpa_t.SetupLimpaGate):
    def check(self, limpa, **env):
        self.prebuilt_env(limpa)
        return self.run_setup("--check", **env)

    def test_limpa_1_2_without_nanoparquet_reads_with_arrow_and_the_de_runs(self):
        r, s = self.check("1.2.5", **self.np_missing())
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertIs(s["ready_for"]["de"], True)
        self.assertEqual(s["de_reader"], {"limpa": "1.2.5", "reads_parquet_with": "arrow",
                                          "ok": True, "maxlfq": True})
        self.assertIs(s["limpa"]["ok"], False)                 # still not the pinned stack
        (note,) = [n for n in s["notes"] if n.startswith("nanoparquet")]
        self.assertNotIn("no DE can run", note)
        self.assertIn("reads report.parquet with arrow", note)
        self.assertIn("Precursor.Normalised", note)            # what the DE loses without it
        self.assertEqual(s["nanoparquet"]["fix"], note[note.index("Fix: ") + 5:].rstrip("."))
        self.assertIn("ERROR: nanoparquet", r.stderr)          # the pinned stack is not there

    def test_limpa_1_4_without_nanoparquet_cannot_read_the_report(self):
        r, s = self.check("1.4.1", **self.np_missing())
        self.assertIs(s["ready_for"]["de"], False)
        self.assertEqual(s["de_reader"], {"limpa": "1.4.1", "reads_parquet_with": "nanoparquet",
                                          "ok": False, "maxlfq": True})
        (note,) = [n for n in s["notes"] if n.startswith("nanoparquet")]
        self.assertIn("no DE can run", note)

    def test_without_arrow_no_de_runs_on_either_limpa(self):
        for limpa in ("1.2.5", "1.4.1"):
            with self.subTest(limpa=limpa):
                self.tearDown()
                self.setUp()
                _r, s = self.check(limpa, FAKE_ARROW_MISSING="1")
                self.assertIs(s["ready_for"]["de"], False)
                self.assertIs(s["de_reader"]["ok"], False)
                self.assertTrue([n for n in s["notes"] if n.startswith("arrow cannot be loaded")],
                                s["notes"])

    def test_the_pinned_stack_reads_with_nanoparquet(self):
        _r, s = self.check("1.4.1")
        self.assertIs(s["ready_for"]["de"], True)
        self.assertEqual(s["de_reader"]["reads_parquet_with"], "nanoparquet")
        self.assertFalse([n for n in s["notes"] if n.startswith(("nanoparquet", "arrow"))])

    def test_no_limpa_no_de_but_maxlfq_can_run(self):
        """run_de.R --method maxlfq needs arrow and limma, no limpa (review of this branch)."""
        _r, s = self.check("")
        self.assertIs(s["ready_for"]["de"], False)
        self.assertEqual(s["de_reader"]["reads_parquet_with"], "")
        self.assertIs(s["de_reader"]["maxlfq"], True)

    def test_a_devel_limpa_is_read_by_its_description_not_its_number(self):
        """A limpa between 1.2 and 1.4 (a devel build): its DESCRIPTION says which reader."""
        _r, s = self.check("1.3.9", FAKE_LIMPA_SUGGESTS="nanoparquet", **self.np_missing())
        self.assertEqual(s["de_reader"]["reads_parquet_with"], "nanoparquet")
        self.assertIs(s["ready_for"]["de"], False)
        (note,) = [n for n in s["notes"] if n.startswith("nanoparquet")]
        self.assertIn("no DE can run", note)


# the parent's own tests run in test_setup_limpa.py; not again here
for _name in [n for n in dir(DeReadinessFollowsTheReader) if n.startswith("test_")]:
    if _name not in DeReadinessFollowsTheReader.__dict__:
        setattr(DeReadinessFollowsTheReader, _name, None)
del _name

if __name__ == "__main__":
    unittest.main()
