#!/usr/bin/env python3
"""
2.15.1 audit (2026-11-17) L2-11: the audit's confounding check fired only when EVERY group had
its own single batch. Groups A and B in batch 1 and C in batch 2 was reported "'Batch' is not
confounded with Group" -- false for the A-C and B-C comparisons. Now each pair of groups that
shares no batch value is named (WARN). Pure Python.
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import audit_results as ar  # noqa: E402


def rows(spec):
    return [{"File.Name": f"r{i}", "Group": g, "Batch": b} for i, (g, b) in enumerate(spec)]


def confounding(spec):
    f = []
    ar.audit_conditions(f, rows(spec), 0.05)
    return [x for x in f if x["check"] == "confounding"]


class PartialConfounding(unittest.TestCase):
    def test_partial_is_warned_per_pair(self):
        out = confounding([("A", "1"), ("A", "1"), ("B", "1"), ("B", "1"), ("C", "2"), ("C", "2")])
        self.assertEqual([x["status"] for x in out], ["WARN"])
        m = out[0]["message"]
        self.assertIn("A (1) and C (2) share no 'Batch' value", m)
        self.assertIn("B (1) and C (2) share no 'Batch' value", m)
        self.assertNotIn("A (1) and B (1)", m)

    def test_perfect_still_fails(self):
        out = confounding([("A", "1"), ("A", "1"), ("B", "2"), ("B", "2")])
        self.assertEqual([x["status"] for x in out], ["FAIL"])

    def test_crossed_passes(self):
        out = confounding([("A", "1"), ("A", "2"), ("B", "1"), ("B", "2"), ("C", "1"), ("C", "2")])
        self.assertEqual([x["status"] for x in out], ["PASS"])


if __name__ == "__main__":
    unittest.main()
