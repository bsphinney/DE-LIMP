#!/usr/bin/env python3
"""Two mass-accuracy refusals that could not change the pin.

1. A ONE-FILE search (a QC run, a method test) always failed the probe: "only 1 run logged
   everything asked; mass accuracy is pinned from the median of at least 2 representative runs"
   (Exploris 480 60k/15k, msalemi 2026-10-05: MS2 measured 12 ppm, job FAILED after the library,
   nothing searched). A cohort of one run has no second run to measure, so its own measurement is
   pinned -- inside the band, floored at the SOP like any other. Cohorts of 2+ keep the rule.
2. Step 1b refused an MS2 spread over 50% on 18 streptavidin pull-down runs (20/24/10 ppm,
   2026-10-08) although the median, 20, IS the SOP floor, so 20 was the pin either way. A spread
   whose median is at or below the floor now WARNS -- when every run is inside the band -- and the
   refusal names enrichment cohorts with sparse runs as a known benign cause.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import probe_window  # noqa: E402
from test_step1b_mass_accuracy_probe import (FAKE_DIANN, GB, _cohort, _exe, _fasta_lib, _run,  # noqa: E402
                                             _stem, ms2_as, no_ms2, real_log)


def probe(d, raws):
    """probe_window.py measuring MS2 with MS1 documented at 7, as the single-shot job runs it."""
    fasta, lib = _fasta_lib(d)
    diann = _exe(os.path.join(d, "diann"), FAKE_DIANN)
    cfg = os.path.join(d, "resolved.cfg")
    with open(cfg, "w") as fh:
        fh.write("--qvalue 0.01")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "probe_window.py"),
                        "--diann", diann, "--raw", *raws, "--fasta", fasta, "--lib", lib,
                        "--threads", "8", "--timeout", "60", "--write-cfg", cfg,
                        "--workdir", os.path.join(d, "w"), "--measure", "mass-acc",
                        "--ms1-ppm", "7"],
                       capture_output=True, text=True, timeout=240,
                       env=dict(os.environ, DOTNET_ROOT="/opt/fake-dotnet"))
    with open(cfg) as fh:
        return p, json.loads(p.stdout), fh.read()


def with_ms2(raws, values):
    for r, v in zip(raws, values):
        with open(r + ".log", "w") as fh:
            fh.write(ms2_as(real_log(_stem(r)), v))
    return raws


class OneRunCohort(unittest.TestCase):
    def test_the_rule_is_one_run_only_for_a_cohort_of_one(self):
        self.assertEqual(probe_window.mass_acc_min_runs(1), 1)
        for n in (2, 3, 40):
            self.assertEqual(probe_window.mass_acc_min_runs(n), probe_window.MASS_ACC_MIN_RUNS)

    def test_a_one_file_search_pins_its_own_measurement_floored_at_the_sop(self):
        with tempfile.TemporaryDirectory() as d:
            raws = [_run(d, "Ex01162023_12_TT34", 13 * GB // 10)]      # MS2 14 ppm, MS1 4.1
            p, out, cfg = probe(d, raws)
            self.assertEqual(p.returncode, 0, p.stderr)
            ma = out["mass_acc"]
            self.assertTrue(ma["single_run_cohort"])
            self.assertEqual((ma["measured_ms2_ppm"], ma["pinned_ms2_ppm"]), (14.0, 20.0))
            self.assertEqual(ma["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
            self.assertIn("the cohort is ONE run", ma["rule"])
            self.assertIn("--mass-acc 20\n", cfg)
            self.assertIn("the cohort is ONE run", p.stderr)

    def test_a_one_file_search_still_has_to_be_plausible(self):
        with tempfile.TemporaryDirectory() as d:
            raws = with_ms2([_run(d, "Ex01162023_12_TT34", 13 * GB // 10)], [60])
            p, out, cfg = probe(d, raws)
            self.assertEqual(p.returncode, probe_window.EXIT_REFUSED, p.stderr)
            self.assertIsNone(out["mass_acc"])
            self.assertNotIn("--mass-acc", cfg)

    def test_a_cohort_of_two_with_one_answer_is_still_refused(self):
        with tempfile.TemporaryDirectory() as d:
            raws = _cohort(d)[:2]
            with open(raws[1] + ".log", "w") as fh:
                fh.write(no_ms2(real_log(_stem(raws[1]))))
            p, out, cfg = probe(d, raws)
            self.assertNotEqual(p.returncode, 0)
            self.assertTrue(any("at least 2" in r for r in out["mass_acc_refused"]),
                            out["mass_acc_refused"])
            self.assertNotIn("--mass-acc", cfg)


class SpreadAtTheFloor(unittest.TestCase):
    def test_a_spread_whose_median_is_the_floor_warns_and_pins_the_floor(self):
        """The pull-down cohort's own numbers: 20/24/10 ppm, 70% of the median, median 20."""
        pin = probe_window.pin_mass_acc([{"ms2_ppm": v} for v in (20.0, 24.0, 10.0)],
                                        documented={"ms1_ppm": 7})
        self.assertEqual(pin["rejected"], [])
        self.assertEqual(pin["pin_as"], "--mass-acc 20 --mass-acc-ms1 7")
        self.assertEqual(len(pin["warnings"]), 1)
        self.assertIn("disagree on MS2", pin["warnings"][0])
        self.assertIn("NOT refused", pin["warnings"][0])
        self.assertIn("pull-down", pin["warnings"][0])
        # below the floor too: 4/4/14 is pinned at 20, warned
        low = probe_window.pin_mass_acc([{"ms2_ppm": v} for v in (4.0, 4.0, 14.0)],
                                        documented={"ms1_ppm": 7})
        self.assertEqual((low["rejected"], low["pin_as"]), ([], "--mass-acc 20 --mass-acc-ms1 7"))
        self.assertTrue(low["floored"]["ms2_ppm"])

    def test_a_spread_still_refuses_above_the_floor_or_with_an_implausible_run(self):
        above = probe_window.pin_mass_acc([{"ms2_ppm": v} for v in (24.0, 24.0, 10.0)],
                                          documented={"ms1_ppm": 7})
        self.assertIsNone(above["pin_as"])
        self.assertTrue(any("disagree on MS2" in r for r in above["rejected"]))
        self.assertEqual(above["warnings"], [])
        # median 14 is under the floor, but one run's 45 ppm is no mass accuracy at all
        wild = probe_window.pin_mass_acc([{"ms2_ppm": v} for v in (14.0, 14.0, 45.0)],
                                         documented={"ms1_ppm": 7})
        self.assertIsNone(wild["pin_as"])
        self.assertTrue(any("disagree on MS2" in r for r in wild["rejected"]))

    def test_end_to_end_the_warning_is_in_the_job_log_and_the_floor_in_the_cfg(self):
        with tempfile.TemporaryDirectory() as d:
            p, out, cfg = probe(d, with_ms2(_cohort(d), [20, 24, 10]))
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("--mass-acc 20\n", cfg)
            self.assertIn("[probe_window] WARNING: the probed runs disagree on MS2", p.stderr)
            self.assertEqual(len(out["mass_acc"]["warnings"]), 1)

    def test_the_refusal_names_enrichment_cohorts_as_a_benign_cause(self):
        with tempfile.TemporaryDirectory() as d:
            p, out, cfg = probe(d, with_ms2(_cohort(d), [24, 24, 10]))
            self.assertEqual(p.returncode, probe_window.EXIT_REFUSED, p.stderr)
            self.assertIn("Refusing to pin", p.stderr)
            self.assertIn("pull-down", p.stderr)
            self.assertIn("ask before anything is pinned by hand", p.stderr)
            self.assertNotIn("--mass-acc", cfg)


if __name__ == "__main__":
    unittest.main()
