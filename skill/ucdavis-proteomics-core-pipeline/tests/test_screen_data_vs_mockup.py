#!/usr/bin/env python3
"""make_screen_data.R reproduces the numbers in the design mockup (trace_data13.json) on real data.

Runs only when both are set, and skips with the reason otherwise:
  SCREEN_MOCKUP_DIR   folder holding the mockup's trace_data13.json (read only)
  SCREEN_HELA_SESSION the finished HeLa limpa analysis the mockup was made from (read only)
Rscript with limpa is needed too. The script is run once per session (setUpClass), into a temp folder.
A mismatch here is a bug in the screen data, not in this test: do not widen a tolerance.

Tolerances: shares 1e-3, detection 1e-3, weights 1e-4, protein values 1e-6. Protein values are written
to 4 decimals (references/protein-screen.md), so they are compared with the mockup's value rounded the
same way; the unrounded difference is also bounded by that rounding (5e-5) plus the 1e-6.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import make_protein_screen as mps  # noqa: E402  (the builder's own reader: index.tsv -> part -> record)

MOCKUP_DIR = os.environ.get("SCREEN_MOCKUP_DIR")
SESSION = os.environ.get("SCREEN_HELA_SESSION")
PGS = ["P80723", "O75683", "Q6S8J3"]
TOL = {"shares": 1e-3, "detection": 1e-3, "weights": 1e-4, "protein": 1e-6}
ROUNDING = 5e-5          # half a unit of the 4th decimal

REASON = None
if not MOCKUP_DIR or not SESSION:
    REASON = "set SCREEN_MOCKUP_DIR and SCREEN_HELA_SESSION to run this (real HeLa data); both are unset or one is"
elif not os.path.isfile(os.path.join(MOCKUP_DIR, "trace_data13.json")):
    REASON = "no trace_data13.json in SCREEN_MOCKUP_DIR=%s" % MOCKUP_DIR
elif not os.path.isdir(os.path.join(SESSION, "output")):
    REASON = "SCREEN_HELA_SESSION=%s has no output/ folder" % SESSION
elif shutil.which("Rscript") is None:
    REASON = "Rscript is not on PATH"


def worst(a, b):
    """Largest |a - b| over equal-shaped nested lists; None in both places is skipped."""
    if isinstance(a, list):
        assert isinstance(b, list) and len(a) == len(b), "shape differs: %s vs %s" % (len(a), len(b))
        return max([worst(x, y) for x, y in zip(a, b)] or [0.0])
    if a is None and b is None:
        return 0.0
    assert a is not None and b is not None, "null in one only: ours %r, mockup %r" % (a, b)
    return abs(a - b)


@unittest.skipIf(REASON, REASON or "")
class ScreenDataVsMockup(unittest.TestCase):
    maxDiff = None

    @classmethod
    def setUpClass(cls):
        with open(os.path.join(MOCKUP_DIR, "trace_data13.json")) as fh:
            cls.mock = json.load(fh)
        cls.tmp = tempfile.mkdtemp(prefix="screen_vs_mockup_")
        cls.out = os.path.join(cls.tmp, "screen_data")      # never inside the read-only session
        cls.proc = subprocess.run(
            ["Rscript", os.path.join(SCRIPTS, "make_screen_data.R"), "--session", SESSION,
             "--proteins", ",".join(PGS), "--out", cls.out], capture_output=True, text=True)
        if cls.proc.returncode != 0:
            raise AssertionError("make_screen_data.R exited %s:\n%s" % (cls.proc.returncode, cls.proc.stderr))
        cls.manifest = mps.load_manifest(cls.out)
        cls.index = mps.parse_index(os.path.join(cls.out, "index.tsv"))
        cls.worst = {}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)
        for k in sorted(cls.worst):
            print("  largest difference %-30s %.3g" % (k, cls.worst[k]), file=sys.stderr)

    def record(self, pg):
        rows = [r for r in self.index if r["Protein.Group"] == pg]
        self.assertEqual(len(rows), 1, pg)
        return mps.read_record(self.out, rows[0], self.manifest)

    def reordered(self, rec, ids, ours):
        """Our per-precursor rows (`ours`, in rec.precursors order) in the mockup's order `ids`, by id."""
        mine = [p["id"] for p in rec["precursors"]]
        self.assertEqual(sorted(mine), sorted(ids))
        return [ours[mine.index(i)] for i in ids]

    def check(self, label, pg, got, want, tol):
        d = worst(got, want)
        key = "%s %s" % (label, pg)
        self.worst[key] = d
        self.assertLessEqual(d, tol, "%s: largest difference %.3g exceeds %g" % (key, d, tol))

    # ---- run level ----
    def test_samples_and_run_stats(self):
        self.assertEqual(self.manifest["samples"], self.mock["samples"])
        # the HeLa session has no search/report.stats.tsv
        self.assertIsNone(self.manifest["run_stats"])
        self.assertTrue([s for s in self.manifest["skipped"] if s["block"] == "run_stats"],
                        self.manifest["skipped"])
        self.assertEqual(self.manifest["shares"]["withheld"], [])      # a finding if not (cross-CPU refit)
        self.assertEqual(self.manifest["detection"]["withheld"], [])

    # ---- per protein ----
    def test_shares_match_contrib(self):
        for pg in PGS:
            rec, c = self.record(pg), self.mock["contrib"][pg]
            self.assertNotIn("withheld", rec["shares"], "%s shares withheld: %s" % (pg, rec["shares"]))
            own = self.reordered(rec, c["ids"], rec["shares"]["own"])
            # a gap is null here and 0 in the mockup
            for mine, theirs in zip(own, c["own"]):
                for a, b in zip(mine, theirs):
                    if a is None:
                        self.assertEqual(b, 0, pg)
            own0 = [[0 if a is None else a for a in row] for row in own]
            self.check("shares.own", pg, own0, c["own"], TOL["shares"])
            self.check("shares.other", pg, rec["shares"]["other"], c["other"], TOL["shares"])
            self.check("shares.total", pg, rec["shares"]["total"], c["total"], TOL["shares"])

    def test_detection_matches_detprob(self):
        for pg in PGS:
            rec, d = self.record(pg), self.mock["detprob"][pg]
            self.assertNotIn("withheld", rec["detection"], "%s detection withheld" % pg)
            self.check("detection.pf_fit", pg, self.reordered(rec, d["ids"], rec["detection"]["pf_fit"]),
                       d["pf_fit"], TOL["detection"])
            self.check("detection.pf_obs", pg, self.reordered(rec, d["ids"], rec["detection"]["pf_obs"]),
                       d["pf_obs"], TOL["detection"])

    def test_weights_match_routes(self):
        for pg in PGS:
            rec, r = self.record(pg), self.mock["routes"][pg]
            self.assertIsNotNone(rec["weights"], pg)
            self.check("weights.w", pg, rec["weights"]["w"], r["w"], TOL["weights"])
            self.check("weights.share_in_group", pg, rec["weights"]["share_in_group"], r["wshare"], TOL["weights"])
            self.assertEqual(rec["weights"]["df_residual"], r["df_res"], pg)

    def test_protein_values_match(self):
        for pg in PGS:
            rec, want = self.record(pg), self.mock["proteins"][pg]["prot"]
            self.check("values.protein (rounded to 4 dp)", pg, rec["values"]["protein"],
                       [None if x is None else round(x, 4) for x in want], TOL["protein"])
            self.check("values.protein (unrounded)", pg, rec["values"]["protein"], want, ROUNDING + TOL["protein"])


if __name__ == "__main__":
    unittest.main()
