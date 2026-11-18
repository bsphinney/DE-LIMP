#!/usr/bin/env python3
"""
2.15.1 audit (2026-11-17), Sage cfg: what the cfg says it searched is what Sage searches.

L2-5  "chimera": true was written for every DDA search with report_psms unset (1). Sage's chimeric
      search needs report_psms > 1 to do anything (Sage docs, configuration/additional; v0.14.7
      scoring.rs score_chimera_fast), so it was a no-op the rationale called "on". Now off, and
      the rationale says why.
L2-8  max_fragment_charge was 2 (Sage's LOW-res recommendation) for high-res MS2 too. Sage's docs
      recommend 1 for high-res deisotoped MS2, 2 for low-res: now by MS2 analyser.
L2-9  the 7-30 length cap stays (a facility choice: the broad-digest memory was measured with it),
      and the enzyme rationale says it is the facility's, not an engine default.

Hermetic: build_sage() in-process.
"""
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import estimate_params as ep  # noqa: E402


def build(acq="DDA", analyzer=None, overrides=None):
    text, rationale = ep.build_sage(acq, "orbitrap_generic", None, overrides or {},
                                    ms2_analyzer=analyzer)
    return json.loads(text), rationale


class ChimeraOff(unittest.TestCase):
    def test_dda_cfg_has_no_chimera_without_report_psms(self):
        cfg, rat = build("DDA")
        self.assertIs(cfg["chimera"], False)
        self.assertNotIn("report_psms", cfg)
        self.assertIs(rat["chimera"]["value"], False)
        self.assertIn("report_psms > 1", rat["chimera"]["source"])

    def test_no_cfg_claims_chimera_with_one_psm(self):
        """Whatever the acquisition, never chimera true with report_psms 1 (Sage's default)."""
        for acq in ("DDA", "DIA"):
            cfg, _ = build(acq)
            if cfg.get("chimera"):
                self.assertGreater(cfg.get("report_psms", 1), 1, acq)

    def test_user_can_still_turn_it_on(self):
        cfg, rat = build("DDA", overrides={"chimera": True, "report_psms": 2})
        self.assertIs(cfg["chimera"], True)
        self.assertEqual(cfg["report_psms"], 2)
        self.assertNotIn("report_psms > 1, which this skill does not use",
                         rat["chimera"]["source"])


class FragmentCharge(unittest.TestCase):
    def test_high_res_ms2_is_one(self):
        cfg, rat = build("DDA", analyzer=None)
        self.assertTrue(cfg["deisotope"])
        self.assertEqual(cfg["max_fragment_charge"], 1)
        self.assertEqual(rat["max_fragment_charge"]["value"], 1)
        self.assertIn("high-res", rat["max_fragment_charge"]["source"])

    def test_ion_trap_ms2_is_two(self):
        for an in sorted(ep.ION_TRAP_ANALYZERS):
            cfg, rat = build("DDA", analyzer=an)
            self.assertFalse(cfg["deisotope"], an)
            self.assertEqual(cfg["max_fragment_charge"], 2, an)
            self.assertIn("low-res", rat["max_fragment_charge"]["source"])

    def test_override_is_named(self):
        cfg, rat = build("DDA", overrides={"max_fragment_charge": 3})
        self.assertEqual(cfg["max_fragment_charge"], 3)
        self.assertEqual(rat["max_fragment_charge"]["value"], 3)


class LengthCapIsTheFacilitys(unittest.TestCase):
    def test_enzyme_rationale_says_so(self):
        cfg, rat = build("DDA")
        self.assertEqual((cfg["database"]["enzyme"]["min_len"], cfg["database"]["enzyme"]["max_len"]),
                         (7, 30))
        self.assertIn("facility's DDA cap", rat["enzyme"]["source"])
        self.assertIn("5-50", rat["enzyme"]["source"])

    def test_not_said_when_lengths_overridden(self):
        cfg, rat = build("DDA", overrides={"database": {"enzyme": {"max_len": 50}}})
        self.assertEqual(cfg["database"]["enzyme"]["max_len"], 50)
        self.assertNotIn("facility's DDA cap", json.dumps(rat))


if __name__ == "__main__":
    unittest.main()
