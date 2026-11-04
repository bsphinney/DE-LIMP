#!/usr/bin/env python3
"""
A Sage search refused before anything is planned leaves no search folder behind (H2c N1).

The refusals in mzml_plan() (no converter for a .raw, no ProteoWizard container for a .d to
quantify, two inputs with one run name) say "Nothing was converted, written or submitted" -- but
run_sage() made the --out folder first, so an empty search_out/ was left behind, which the next
step can take for a search that started. mzml_plan creates nothing, so the folder is now made
only once it has passed. A folder that already existed (an earlier search) is left as it was.

Fakes on a temporary PATH (tests/test_sage_sbatch_conversion._Harness); nothing is searched.
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
# Module imports, not `from ... import`: a TestCase bound here would be collected twice.
import test_sage_sbatch_conversion as ssc  # noqa: E402
import test_sage_dotd_pwiz_conversion as pwz  # noqa: E402
# self.env() is job_env() (see _Harness.env); imported so the guard can see that.
from job_env import job_env  # noqa: E402,F401

NOTHING = "Nothing was converted, written or submitted"


class NoConverter(ssc._Harness):
    def test_the_refusal_leaves_no_search_folder(self):
        for n in ("ThermoRawFileParser", "msconvert"):
            os.remove(os.path.join(self.bin, n))
        p = self.run_search("--sbatch", os.path.join(self.d, "sage_job.sh"))
        self.assertNotEqual(p.returncode, 0)
        self.assertIn(NOTHING, p.stderr)
        self.assertFalse(os.path.exists(self.out), os.listdir(self.d))

    def test_an_earlier_search_folder_is_left_as_it_was(self):
        """(pin) the refusal removes nothing it did not make."""
        os.makedirs(self.out)
        keep = os.path.join(self.out, "report.parquet")
        with open(keep, "w") as fh:
            fh.write("an earlier search's report")
        for n in ("ThermoRawFileParser", "msconvert"):
            os.remove(os.path.join(self.bin, n))
        p = self.run_search("--sbatch", os.path.join(self.d, "sage_job.sh"))
        self.assertIn(NOTHING, p.stderr)
        self.assertEqual(os.listdir(self.out), ["report.parquet"])


class NoContainer(pwz._PwizHarness):
    def test_a_missing_container_leaves_no_search_folder(self):
        job, p = self.generate(sif=os.path.join(self.d, "nowhere", "pwiz.sif"), expect=None)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn(NOTHING, p.stderr)
        self.assertFalse(os.path.exists(self.out), os.listdir(self.d))


if __name__ == "__main__":
    unittest.main()
