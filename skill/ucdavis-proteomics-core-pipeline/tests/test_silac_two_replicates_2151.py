#!/usr/bin/env python3
"""
2.15.1 audit (2026-11-17) lane4 M24: SILAC with 2 replicates tested nothing, silently. The
default minimum, max(3, ceiling(n / 2)), is 3 -- above n = 2 -- so no protein reached the test and
the run "succeeded" with empty DE tables. Now run_silac.R refuses, naming the comparisons and the
choice (--min-replicates 2, said in the Methods as given; or drop them); with --min-replicates 2
it runs and tests. Synthetic triple-SILAC report (silac_fixture). Needs R.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import silac_fixture as fx  # noqa: E402
from test_silac_analysis import HAVE_R, NEEDS, run_silac  # noqa: E402


@unittest.skipUnless(HAVE_R, NEEDS)
class TwoReplicates(unittest.TestCase):
    def test_refused_then_run_with_the_users_minimum(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet"), n_runs=2)
            fx.write_design(os.path.join(d, "design.csv"), n_runs=2)
            p, _ = run_silac(d, "--comparisons", "KO-WT", ok=False)
            self.assertNotEqual(p.returncode, 0)
            err = " ".join(p.stderr.split())
            self.assertIn("too few replicates to test: KO-WT has 2 replicate(s)", err)
            self.assertIn("--min-replicates 2", err)
            p, out = run_silac(d, "--comparisons", "KO-WT", "--min-replicates", "2")
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            t = prov["de_tables"]["KO-WT"]
            self.assertEqual((t["min_replicates"], t["n_replicates_available"]), (2, 2))
            self.assertGreater(t["n_tested"], 0)
            self.assertEqual(prov["settings"]["min_replicates"]["source"],
                             "given (--min-replicates)")


if __name__ == "__main__":
    unittest.main()
