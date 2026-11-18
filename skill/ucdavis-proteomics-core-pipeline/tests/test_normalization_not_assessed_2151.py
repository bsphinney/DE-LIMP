#!/usr/bin/env python3
"""
2.15.1 audit (2026-11-17), normalization_check.py: a check that could not run is not agreement.

L2-2  On a report the check cannot read (every Sage / FragPipe report: B1-B3 need DIA-NN's
      columns) all three data checks were "not assessed", nothing tripped, and the record said
      "the data check agreed" and decided the default by itself. Now: what could not be assessed
      is listed (`not_assessed`), data_check_agrees is false, and a decision is asked for.
L2-13 With no raw DE the recommendation read "(n anomalies vs 0 for raw)": a version with no DE
      is "not assessed", never 0 anomalies.

Pure Python: run()'s record-building pieces on hand-made records.
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402


def rec(default=et.NORMALISED, checks_ok=True, raw_volcano=True, etype="whole_proteome",
        anomalies_norm=0):
    na = {"status": "not_assessed", "why": "the report has no Precursor.Quantity column"}
    ok = {"status": "assessed"}
    vol = {"contrasts": {"B-A": {"anomalies": ["x"] * anomalies_norm}}}
    r = {"default": {"quantities": default, "why": "type default"},
         "experiment_type": {"type": etype},
         "checks": {"factor_confound": dict(ok) if checks_ok else dict(na),
                    "raw_shift": dict(ok) if checks_ok else dict(na),
                    "id_gap": dict(ok) if checks_ok else dict(na),
                    "instability": {"status": "not_assessed", "why": "no column"}},
         "volcano": {et.NORMALISED: vol,
                     et.RAW: {"contrasts": {"B-A": {"anomalies": []}}} if raw_volcano else None},
         "volcano_not_assessed": {} if raw_volcano else {et.RAW: "no --de-raw given"}}
    r["trips"] = nc.trips(r)
    r["tripped"] = bool(r["trips"])
    r["not_assessed"] = nc.not_assessed(r)
    r["data_check_agrees"] = not r["tripped"] and not r["not_assessed"]
    return r


class NotAssessedIsNotAgreement(unittest.TestCase):
    def test_unreadable_report_lists_b1_to_b3(self):
        r = rec(checks_ok=False)
        self.assertFalse(r["tripped"])
        self.assertEqual(len(r["not_assessed"]), 3)
        self.assertTrue(r["not_assessed"][0].startswith("B1 normalisation factors vs groups: "
                                                        "not assessed"))
        self.assertFalse(r["data_check_agrees"])
        pick, why = nc._recommend(r)
        self.assertEqual(pick, et.NORMALISED)
        self.assertNotIn("agreed", why)
        self.assertIn("not a data-checked one", why)

    def test_b4_instability_alone_is_not_required(self):
        r = rec()
        self.assertEqual(r["not_assessed"], [])
        self.assertTrue(r["data_check_agrees"])
        self.assertEqual(nc._recommend(r)[1],
                         "the default for this experiment type; the data check agreed")

    def test_the_defaults_own_volcano_is_required(self):
        r = rec(default=et.RAW, raw_volcano=False, etype="ip")
        self.assertEqual(r["not_assessed"], ["the raw volcano: not assessed (no --de-raw given)"])
        self.assertFalse(r["data_check_agrees"])

    def test_raw_default_needs_b1_only(self):
        r = rec(default=et.RAW, etype="ip")
        r["checks"]["raw_shift"] = {"status": "not_assessed", "why": "x"}
        self.assertEqual(nc.not_assessed(r), [])

    def test_summary_says_not_fully_assessed(self):
        r = rec(checks_ok=False)
        r.update(recommendation="x", decision=None, needs_decision=True, report="r",
                 conditions="c", groups={"A": 3, "B": 3}, design={}, bait=None,
                 thresholds={}, experiment_type={"type": "whole_proteome", "label": "Whole",
                                                  "source": "user"})
        md = nc.summary_md(r)
        self.assertIn("NOT FULLY ASSESSED", md)
        self.assertNotIn("agreed with the default", md)


class MissingVolcanoIsNotZero(unittest.TestCase):
    def test_counts(self):
        r = rec(raw_volcano=False, anomalies_norm=2)
        n = nc._anomaly_counts(r)
        self.assertEqual(n, {et.NORMALISED: 2, et.RAW: None})

    def test_recommendation_says_not_assessed(self):
        r = rec(raw_volcano=False, anomalies_norm=2)
        self.assertTrue(r["tripped"])
        pick, why = nc._recommend(r)
        self.assertEqual(pick, et.NORMALISED)
        self.assertIn("raw: not assessed", why)
        self.assertNotIn("vs 0", why)

    def test_no_type_and_one_volcano_cannot_decide(self):
        r = rec(raw_volcano=False, etype=None, anomalies_norm=5)
        pick, why = nc._recommend(r)
        self.assertIn("the volcanoes cannot decide", why)


if __name__ == "__main__":
    unittest.main()
