#!/usr/bin/env python3
"""A normalisation check made before 2.14.2, and a final MaxLFQ DE made with 2.14.2 (C1 review, F1).

A session whose two check DEs ran under 2.14.1 reaches its final DE under 2.14.2. The check record
names no between-run step, so it counts as quantile -- 2.14.1 quantile-normalised every
normalised MaxLFQ DE (normalization_check.LEGACY_MAXLFQ_BETWEEN_RUN). The final DE asks for the
default, auto. The gate compared the REQUESTED steps (quantile vs auto) and refused, and its first
way forward was "Run this DE with --between-run quantile": on a DIA-NN report that is 2.14.1's
quantile step on top of DIA-NN's normalisation again (I20: 61-77% of the significant calls wrong
on the PXD063416 benchmark), and the Methods then said the check agreed. On a Sage/FragPipe report auto IS
quantile, so the refusal stopped a step identical to the one the check's DE had.

Now, for a check from before 2.14.2, the gate compares the step this DE APPLIES with the old
rule's quantile. For auto that is known only once run_de.R has read the report (build_maxlfq.R
maxlfq_between_run, the one rule), so the gate first says the comparison is pending and run_de.R
asks it again with the applied step (--between-run-applied):

  * DIA-NN report, auto -> none: refused; the first way forward is to run the check again with
    this skill; --between-run quantile is said to repeat 2.14.1's quantile step on top of
    DIA-NN's normalisation, for reproducing a 2.14.1 analysis only. The harm is said to be the
    quantile step itself, not normalising twice (C2 statistician L1: alone, on un-normalised
    quantities, it gave 46-72% wrong calls);
  * Sage report, auto -> quantile: the same step as the check's DE -- it runs, "design matches";
  * an explicit other step (none, median): refused at once, the same way forward first.

A check made by 2.14.2 records its normalised DE's requested step, and is held to it as before.

Unit tests are stdlib; the end-to-end class needs pyarrow and R with limma/arrow/dplyr/tidyr/
jsonlite (skips otherwise). Hermetic.
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
sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402
from test_maxlfq_between_run import HAVE_ARROW, NEEDS, SAGE_META, write_report  # noqa: E402
from test_run_de_contaminants import r_has  # noqa: E402

CHECK = os.path.join(SCRIPTS, "normalization_check.py")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")

MAXLFQ = {"method": "maxlfq", "formula": "~ 0 + groups", "numeric_covariates": [],
          "block_column": None, "block_effect": None, "block_scope": None,
          "contrasts": ["B-A"], "flags": ["--contrasts", "B-A"]}
OLD = {"check_design": dict(MAXLFQ, status="recorded")}      # 2.14.1: no between_run
RERUN = "run the normalisation check again with this skill"


def old_record(path, report=None, conditions=None):
    """A decided check record as 2.14.1 wrote it (no between_run in check_design) -> its path."""
    rec = dict(OLD, schema=nc.SCHEMA, schema_version=nc.SCHEMA_VERSION, tripped=False, trips=[],
               de_method="maxlfq",
               decision={"quantities": et.NORMALISED, "source": "default", "by": "the skill",
                         "reason": "the default for this experiment type"},
               default={"quantities": et.NORMALISED, "why": "x"},
               experiment_type={"type": "whole_proteome", "label": "whole proteome",
                                "source": "user"})
    if report:
        rec.update(report=report, report_sha256=nc.sha256(report))
    if conditions:
        rec.update(conditions=conditions, conditions_sha256=nc.sha256(conditions))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(rec, fh)
    return path


class TheResolvedStepIsCompared(unittest.TestCase):
    def check(self, final, applied=None):
        return nc.design_check(OLD, dict(MAXLFQ, between_run=final), et.NORMALISED,
                               applied=applied)

    def test_auto_is_compared_once_the_report_is_read(self):
        dc = self.check("auto")
        self.assertIsNone(dc["refusal"])
        self.assertEqual(dc["between_run_pending"], "quantile")
        # said, so a record that never got the second answer does not read as "matches"
        self.assertNotEqual(dc["match"], "matches")
        self.assertIn("between-run step", dc["match"])

    def test_dia_nn_auto_none_is_refused_and_running_the_check_again_comes_first(self):
        r = self.check("auto", applied="none")["refusal"]
        self.assertIsNotNone(r)
        self.assertIn(RERUN, r)
        self.assertIn("--between-run quantile", r)
        self.assertLess(r.index(RERUN), r.index("--between-run quantile"))
        self.assertNotIn("Run this DE with --between-run quantile", r)
        self.assertIn("quantile step on top of DIA-NN's normalisation", r)
        self.assertNotIn("double normalisation", r)          # C2 statistician L1: not the mechanism
        self.assertIn("not normalising twice", r)
        self.assertIn("46-72%", r)
        self.assertIn("only to reproduce a 2.14.1 analysis", r.lower())
        self.assertIn("before 2.14.2", r)

    def test_sage_auto_quantile_is_the_checks_own_step(self):
        dc = self.check("auto", applied="quantile")
        self.assertEqual((dc["match"], dc["refusal"]), ("matches", None))
        self.assertNotIn("between_run_pending", dc)

    def test_an_explicit_other_step_is_refused_at_once_the_same_way(self):
        for step in ("none", "median"):
            with self.subTest(step=step):
                r = self.check(step)["refusal"]
                self.assertIsNotNone(r)
                self.assertIn(f"--between-run {step}", r)
                self.assertLess(r.index(RERUN), r.index("--between-run quantile"))
                self.assertNotIn("Run this DE with --between-run quantile", r)
        # ... and quantile asked for explicitly is the check's own step
        self.assertEqual(self.check("quantile")["match"], "matches")

    def test_contrasts_are_still_warned_while_the_step_is_pending(self):
        dc = nc.design_check(OLD, dict(MAXLFQ, between_run="auto", contrasts=["A-B"]),
                             et.NORMALISED)
        self.assertEqual(dc["between_run_pending"], "quantile")
        self.assertTrue(any("the check did not look at A-B" in w for w in dc["warnings"]))

    def test_a_2_14_2_check_is_held_to_its_requested_step_as_before(self):
        new = {"check_design": dict(MAXLFQ, status="recorded", between_run="auto")}
        dc = nc.design_check(new, dict(MAXLFQ, between_run="auto"), et.NORMALISED)
        self.assertEqual((dc["match"], dc["refusal"]), ("matches", None))
        self.assertNotIn("between_run_pending", dc)
        dc = nc.design_check(new, dict(MAXLFQ, between_run="quantile"), et.NORMALISED,
                             applied="quantile")
        self.assertIsNotNone(dc["refusal"])
        self.assertNotIn(RERUN, dc["refusal"])

    def test_a_raw_final_de_has_no_step_to_compare(self):
        dc = nc.design_check(OLD, dict(MAXLFQ, between_run="auto"), et.RAW)
        self.assertEqual((dc["match"], dc["refusal"]), ("matches", None))
        self.assertNotIn("between_run_pending", dc)


class TheGateAsksTwice(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.path = old_record(os.path.join(self.d, "nc", nc.RECORD))
        self.design = os.path.join(self.d, "design.json")
        with open(self.design, "w") as fh:
            json.dump(dict(MAXLFQ, between_run="auto"), fh)

    def gate(self, *extra):
        p = subprocess.run([sys.executable, CHECK, "gate", "--check", self.path, "--quantities",
                            et.NORMALISED, "--design-json", self.design, *extra],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        return json.loads(p.stdout)

    def test_first_pending_then_decided_by_the_applied_step(self):
        first = self.gate()
        self.assertTrue(first["ok"], first["message"])
        self.assertEqual(first["between_run_pending"], "quantile")
        dia_nn = self.gate("--between-run-applied", "none")
        self.assertFalse(dia_nn["ok"])
        self.assertIn(RERUN, dia_nn["message"])
        self.assertIsNone(dia_nn["between_run_pending"])
        sage = self.gate("--between-run-applied", "quantile")
        self.assertTrue(sage["ok"], sage["message"])
        self.assertEqual(sage["design_match"], "matches")
        self.assertIsNone(sage["between_run_pending"])

    def test_gate_function(self):
        ok, msg, _ = nc.gate(self.path, et.NORMALISED, design=dict(MAXLFQ, between_run="auto"),
                             applied="none")
        self.assertFalse(ok)
        self.assertIn(RERUN, msg)
        ok, msg, _ = nc.gate(self.path, et.NORMALISED, design=dict(MAXLFQ, between_run="auto"),
                             applied="quantile")
        self.assertTrue(ok, msg)
        self.assertIn("design matches", msg)


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class RunDeWithA2141Check(unittest.TestCase):
    """run_de.R --method maxlfq under a check record 2.14.1 wrote."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def final(self, name, **report):
        rep, cond, _ = write_report(os.path.join(self.d, name), **report)
        path = old_record(os.path.join(self.d, name, "nc", nc.RECORD), rep, cond)
        out = os.path.join(self.d, name, "final")
        p = subprocess.run(["Rscript", RUN_DE, "--input", rep, "--metadata", cond, "--method",
                            "maxlfq", "--quantities", et.NORMALISED, "--contrasts", "B-A",
                            "--normalization-check", path, "--outdir", out],
                           capture_output=True, text=True, timeout=900, cwd=self.d)
        return p, out

    def test_a_dia_nn_report_is_refused_and_told_to_run_the_check_again_first(self):
        p, out = self.final("dn")
        self.assertNotEqual(p.returncode, 0)
        flat = " ".join(p.stderr.split())
        self.assertIn(RERUN, flat)
        self.assertLess(flat.index(RERUN), flat.index("--between-run quantile"))
        self.assertIn("quantile step on top of DIA-NN's normalisation", flat)
        self.assertNotIn("double normalisation", flat)
        self.assertFalse(os.path.exists(os.path.join(out, "de_provenance.json")))

    def test_a_sage_report_runs_its_quantile_step_and_matches(self):
        p, out = self.final("sage", declared=SAGE_META)
        self.assertEqual(p.returncode, 0, p.stderr[-3000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["between_run"]["method"], "quantile")
        self.assertEqual(prov["normalization_check"]["design_match"], "matches")
        self.assertNotIn("between_run_pending", prov["normalization_check"])


if __name__ == "__main__":
    unittest.main()
