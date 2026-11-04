#!/usr/bin/env python3
"""The normalisation check is held to the model its DEs were fitted with (2.12 review, M1).

Before: step8_commands built the two check DEs with no --numeric-covariates / --block /
--contrasts, run() never looked at the design or the inputs its DEs used, and gate() compared
only the report and conditions sha256 and the quantities -- so a final DE with another model (a
covariate, a block, another method) went through on a decision made about a different one. Now:

  * run_de.R records the design it was asked for before fitting (de_provenance.json
    `design_request`: method, formula, numeric covariates, block column/effect/scope, contrasts,
    the flags given) and the conditions' sha256 (`metadata_sha256`);
  * normalization_check.py run refuses two check DEs fitted with different designs, or that read
    another report or conditions file than the check, and records `check_design`;
  * the final DE passes its design to gate (--design-json): another method, formula, numeric
    covariate or block is refused; other contrasts are a warning; a check from before 2.12 (no
    check_design) is NOT CHECKED, said, and the Methods line carries [confirm];
  * the printed step-8 commands carry the asking DE's design flags.

Unit tests are stdlib; the end-to-end class needs pyarrow and R with limpa (skips otherwise).
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
from test_normalization_check import (HAVE_ARROW, NEEDS, make_report, no_staff_list,  # noqa: E402
                                     session)
from test_run_de_contaminants import r_has  # noqa: E402

CHECK = os.path.join(SCRIPTS, "normalization_check.py")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")

DESIGN = {"method": "dpc", "formula": "~ 0 + groups + Covariate1", "numeric_covariates": ["Covariate1"],
          "block_column": "Mouse", "block_effect": "auto", "block_scope": "within",
          "contrasts": ["B-A"], "flags": ["--numeric-covariates", "Covariate1", "--block", "Mouse"]}


def de(design=DESIGN, input_sha="r" * 64, metadata_sha="c" * 64):
    prov = {"input": "/x/report.parquet", "metadata": "/x/conditions.csv",
            "input_sha256": input_sha, "metadata_sha256": metadata_sha}
    if design is not None:
        prov["design_request"] = design
    return {"prov": prov, "tables": {}}


class TheCheckRecordsOneDesign(unittest.TestCase):
    def test_two_check_des_with_one_design_are_recorded(self):
        cd = nc.check_design_of(de(), de(), "r" * 64, None, "c" * 64)
        self.assertEqual(cd["status"], "recorded")
        for k in nc.DESIGN_KEYS + ("contrasts", "flags"):
            self.assertEqual(cd[k], DESIGN[k], k)
        self.assertEqual(cd["files"]["raw"]["metadata_sha256"], "c" * 64)

    def test_different_designs_are_refused(self):
        for change in ({"numeric_covariates": []}, {"block_column": None}, {"method": "maxlfq"},
                       {"formula": "~ 0 + groups"}, {"contrasts": ["A-B"]}):
            with self.subTest(change=change):
                with self.assertRaises(ValueError) as e:
                    nc.check_design_of(de(), de(dict(DESIGN, **change)), "r" * 64, None, "c" * 64)
                self.assertIn("the two check DEs were fitted with different designs", str(e.exception))
                self.assertIn(next(iter(change)), str(e.exception))

    def test_another_report_or_conditions_is_refused(self):
        with self.assertRaises(ValueError) as e:
            nc.check_design_of(de(), de(input_sha="x" * 64), "r" * 64, None, "c" * 64)
        self.assertIn("the raw check DE (/x/report.parquet) read another report than the check's "
                      "--report", str(e.exception))
        with self.assertRaises(ValueError) as e:
            nc.check_design_of(de(metadata_sha="x" * 64), de(), "r" * 64, None, "c" * 64)
        self.assertIn("read another conditions than the check's --conditions", str(e.exception))
        # MaxLFQ's raw DE reads its own --no-norm report: held to --report-raw, not --report
        cd = nc.check_design_of(de(), de(input_sha="n" * 64), "r" * 64, "n" * 64, "c" * 64)
        self.assertEqual(cd["status"], "recorded")

    def test_check_des_from_before_2_12_record_no_design(self):
        cd = nc.check_design_of(de(None), de(), "r" * 64, None, "c" * 64)
        self.assertEqual(cd["status"], "not_recorded")
        self.assertIn("normalised check DE(s) record no design", cd["note"])
        old = de(None)
        del old["prov"]["metadata_sha256"]                  # an older run_de.R: unverified, said
        cd = nc.check_design_of(old, de(None), "r" * 64, None, "c" * 64)
        self.assertEqual(cd["status"], "not_recorded")


class TheGateHoldsTheFinalDeToIt(unittest.TestCase):
    REC = {"check_design": dict(DESIGN, status="recorded")}

    def test_the_same_design_matches(self):
        self.assertEqual(nc.design_check(self.REC, DESIGN),
                         {"match": "matches", "warnings": [], "refusal": None})

    def test_another_model_is_refused(self):
        for change in ({"numeric_covariates": []}, {"block_column": None, "block_effect": None,
                                                    "block_scope": None},
                       {"method": "maxlfq"}, {"formula": "~ 0 + groups + Batch + Covariate1"},
                       {"block_effect": "random"}):
            with self.subTest(change=change):
                dc = nc.design_check(self.REC, dict(DESIGN, **change))
                self.assertEqual(dc["match"], "differs")
                self.assertTrue(dc["refusal"].startswith("REFUSED: this DE's design differs"))
                self.assertIn("--numeric-covariates Covariate1 --block Mouse", dc["refusal"])

    def test_other_contrasts_are_a_warning(self):
        dc = nc.design_check(self.REC, dict(DESIGN, contrasts=["B-A", "C-A"]))
        self.assertEqual((dc["match"], dc["refusal"]), ("matches; contrasts differ", None))
        self.assertIn("the check did not look at C-A", dc["warnings"][0])

    def test_a_check_from_before_2_12_is_not_checked_and_says_so(self):
        for rec in ({}, {"check_design": {"status": "not_recorded", "note": "the raw check DE(s) "
                                                                            "record no design"}}):
            dc = nc.design_check(rec, DESIGN)
            self.assertEqual((dc["match"], dc["refusal"]), ("not recorded", None))
            self.assertIn("NOT RECORDED", dc["warnings"][0])
            self.assertIn("NOT CHECKED", dc["warnings"][0])
            self.assertIn("run step 8's check again", dc["warnings"][0])

    def test_gate_refuses_before_reading_the_decision(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, nc.RECORD)
            with open(path, "w") as fh:
                json.dump(dict(self.REC, schema=nc.SCHEMA, decision={"quantities": et.RAW, "by": "x",
                                                                     "reason": "r"}, tripped=False,
                               trips=[]), fh)
            ok, msg, _ = nc.gate(path, et.RAW, design=dict(DESIGN, method="maxlfq"))
            self.assertFalse(ok)
            self.assertIn("method", msg)
            ok, msg, _ = nc.gate(path, et.RAW, design=DESIGN)
            self.assertTrue(ok, msg)
            self.assertIn("design matches", msg)


class TheCheckCoversTheMaxlfqRoute(unittest.TestCase):
    """2.14.2 (I20; Brett's Q1): the check also covers the MaxLFQ route. On the PXD063416
    benchmark its record said "method dpc" while the MaxLFQ DE that went wrong was never checked:
    the check's DEs, the final DE and the summary must all be about one route."""
    REPORT, COND, S = "/s/output/search/report.parquet", "/s/input/conditions.csv", "/s"

    def test_step8_commands_maxlfq_runs_all_three_des_with_maxlfq(self):
        """(pin) step8_commands already carried --method into every DE it prints."""
        cmds = nc.step8_commands(self.REPORT, self.COND, self.S, method="maxlfq")
        for c in cmds["check"][:2] + [cmds["final"]]:
            self.assertIn("--method maxlfq", c)
        self.assertNotIn("--method dpc", cmds["next"])

    def test_maxlfq_final_de_is_refused_against_a_dpc_check(self):
        """(pin) design_check already refused another --method than the check's DEs ran."""
        rec = {"check_design": dict(DESIGN, status="recorded", method="dpc")}
        dc = nc.design_check(rec, dict(DESIGN, method="maxlfq"))
        self.assertIsNotNone(dc["refusal"])
        self.assertIn("'dpc'", dc["refusal"])
        self.assertIn("'maxlfq'", dc["refusal"])

    def test_check_summary_names_the_route(self):
        """The summary's header names the route the check's DEs ran, even when their design was
        not recorded (DEs from before 2.12)."""
        rec = {"default": {"quantities": et.NORMALISED}, "de_method": "maxlfq"}
        head = nc.summary_md(rec).split("## Checks")[0]
        self.assertIn("method maxlfq", head)
        self.assertIn("--method maxlfq", head)
        self.assertIn("DE route checked:** NOT RECORDED",
                      nc.summary_md({"default": {"quantities": et.NORMALISED}}))

    @unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
    def test_the_record_names_the_des_method(self):
        with tempfile.TemporaryDirectory() as d:
            report, cond, _ = make_report(os.path.join(d, "w"), "whole")
            s = session(os.path.join(d, "w"), "whole_proteome")
            dirs = {}
            for q in (et.NORMALISED, et.RAW):
                dirs[q] = os.path.join(d, "nc", q)
                os.makedirs(dirs[q])
                with open(os.path.join(dirs[q], "de_provenance.json"), "w") as fh:
                    json.dump({"method": "maxlfq", "contrasts": ["B-A"], "de_tables": {},
                               "input_sha256": nc.sha256(report),
                               "metadata_sha256": nc.sha256(cond),
                               "design_request": dict(DESIGN, method="maxlfq")}, fh)
            rec, out = nc.run(report, cond, dirs[et.NORMALISED], dirs[et.RAW], session=s,
                              out=os.path.join(d, "nc"))
            self.assertEqual(rec["de_method"], "maxlfq")
            with open(os.path.join(out, nc.SUMMARY)) as fh:
                self.assertIn("method maxlfq", fh.read().split("## Checks")[0])


class TheGateSeesTheBetweenRunStep(unittest.TestCase):
    """2.14.2 review: run_de.R --between-run was outside the gate, so a final MaxLFQ DE with
    another between-run step than the check's normalised DE passed as "design matches" -- and its
    Methods said the data check (volcano shape) agreed about a volcano drawn without that step.
    The step changes the volcano (I20), so it is held like the method: the check records its
    normalised DE's requested step (check_design `between_run`); a check from before 2.14.2
    counts as quantile (normalization_check.LEGACY_MAXLFQ_BETWEEN_RUN); a final DE reading
    normalised quantities with another step is refused. A raw final DE has no step to compare."""
    MAXLFQ = dict(DESIGN, method="maxlfq")

    def check(self, **cd):
        return {"check_design": dict(self.MAXLFQ, status="recorded", **cd)}

    def test_check_design_records_the_normalised_des_step(self):
        norm, raw = de(dict(self.MAXLFQ, between_run="quantile")), de(dict(self.MAXLFQ,
                                                                            between_run="auto"))
        norm["prov"].update(method="maxlfq", quantities=et.NORMALISED,
                            between_run={"method": "quantile", "requested": "quantile", "why": "x"})
        raw["prov"].update(method="maxlfq", quantities=et.RAW,
                           between_run={"method": "none", "requested": "auto", "why": "x"})
        cd = nc.check_design_of(norm, raw, "r" * 64, None, "c" * 64)   # not refused: raw has none
        self.assertEqual(cd["between_run"], "quantile")
        old = de(self.MAXLFQ)                                 # a 2.14.1 maxlfq DE: no between_run
        old["prov"].update(method="maxlfq", quantities=et.NORMALISED)
        self.assertEqual(nc.check_design_of(old, de(self.MAXLFQ), "r" * 64, None,
                                            "c" * 64)["between_run"], "quantile")

    def test_scenario_a_auto_check_then_a_quantile_final_is_refused(self):
        dc = nc.design_check(self.check(between_run="auto"),
                             dict(self.MAXLFQ, between_run="quantile"), et.NORMALISED)
        self.assertEqual(dc["match"], "differs")
        self.assertIn("--between-run quantile", dc["refusal"])
        self.assertIn("auto", dc["refusal"])

    def test_scenario_b_a_2_14_1_check_then_an_auto_final_is_refused(self):
        # what auto gives is known once the report is read: the step that runs is compared (C1
        # review F1, tests/test_check_before_2142.py); on a DIA-NN report it is none
        dc = nc.design_check(self.check(), dict(self.MAXLFQ, between_run="auto"), et.NORMALISED)
        self.assertEqual(dc["between_run_pending"], "quantile")
        dc = nc.design_check(self.check(), dict(self.MAXLFQ, between_run="auto"), et.NORMALISED,
                             applied="none")
        self.assertIsNotNone(dc["refusal"])
        self.assertIn("before 2.14.2", dc["refusal"])
        self.assertIn("--between-run quantile", dc["refusal"])
        # ... and with the step that check's DEs had, it matches
        self.assertEqual(nc.design_check(self.check(), dict(self.MAXLFQ, between_run="quantile"),
                                         et.NORMALISED)["match"], "matches")

    def test_the_same_step_matches_and_raw_or_dpc_compare_none(self):
        self.assertEqual(nc.design_check(self.check(between_run="auto"),
                                         dict(self.MAXLFQ, between_run="auto"),
                                         et.NORMALISED)["refusal"], None)
        # a raw final DE has no between-run step: the normalised DE's is not its to match
        self.assertIsNone(nc.design_check(self.check(between_run="quantile"),
                                          dict(self.MAXLFQ, between_run="auto"), et.RAW)["refusal"])
        self.assertEqual(nc.design_check({"check_design": dict(DESIGN, status="recorded")},
                                         DESIGN, et.NORMALISED)["match"], "matches")

    def test_gate_refuses_scenario_a(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, nc.RECORD)
            with open(path, "w") as fh:
                json.dump(dict(self.check(between_run="auto"), schema=nc.SCHEMA, tripped=False,
                               trips=[], decision={"quantities": et.NORMALISED, "by": "x",
                                                   "reason": "r"}), fh)
            ok, msg, _ = nc.gate(path, et.NORMALISED,
                                 design=dict(self.MAXLFQ, between_run="quantile"))
            self.assertFalse(ok)
            self.assertIn("--between-run", msg)
            ok, msg, _ = nc.gate(path, et.NORMALISED, design=dict(self.MAXLFQ, between_run="auto"))
            self.assertTrue(ok, msg)

    def test_step8_commands_carry_the_opt_in_to_the_normalised_and_final_des_only(self):
        cmds = nc.step8_commands("/s/output/search/report.parquet", "/s/input/conditions.csv", "/s",
                                 method="maxlfq", design=dict(self.MAXLFQ, between_run="median"))
        norm, raw, _ = cmds["check"]
        self.assertIn("--quantities normalised", norm)
        self.assertIn("--between-run median", norm)
        self.assertNotIn("--between-run", raw)                # the raw arm refuses a step
        self.assertIn("--between-run median", cmds["final"])
        auto = nc.step8_commands("/s/output/search/report.parquet", "/s/input/conditions.csv", "/s",
                                 method="maxlfq", design=dict(self.MAXLFQ, between_run="auto"))
        self.assertNotIn("--between-run", auto["next"])
        with tempfile.TemporaryDirectory() as d:              # ... and required() prints them
            s = session(d, "ip")
            r = nc.required(os.path.join(s, "output", "search", "report.parquet"),
                            os.path.join(s, "input", "conditions.csv"), method="maxlfq",
                            design=dict(self.MAXLFQ, between_run="quantile"))
        self.assertEqual(r["message"].count("--between-run quantile"), 2)

    def test_the_legacy_rule_is_one_definition(self):
        self.assertEqual(nc.LEGACY_MAXLFQ_BETWEEN_RUN, "quantile")
        self.assertEqual(nc.between_run_of({"method": "maxlfq", "quantities": et.NORMALISED}),
                         "quantile")
        self.assertEqual(nc.between_run_of({"method": "maxlfq", "quantities": et.RAW}), "none")
        self.assertIsNone(nc.between_run_of({"method": "dpc"}))
        self.assertEqual(nc.between_run_of({"method": "maxlfq", "between_run": {
            "method": "none", "requested": "auto"}}, "method"), "none")


class TheCommandsCarryTheDesign(unittest.TestCase):
    def test_step8_commands_carry_the_flags_to_all_three_des(self):
        cmds = nc.step8_commands("/s/output/search/report.parquet", "/s/input/conditions.csv", "/s",
                                 design=dict(DESIGN, flags=DESIGN["flags"] + ["--contrasts", "B-A,C-A"]))
        flags = "--numeric-covariates Covariate1 --block Mouse --contrasts B-A,C-A"
        for c in cmds["check"][:2] + [cmds["final"]]:
            self.assertIn(flags, c)
        self.assertNotIn("--numeric-covariates", cmds["check"][2])     # the check itself takes none

    def test_required_prints_them(self):
        with tempfile.TemporaryDirectory() as d:
            s = session(d, "ip")
            r = nc.required(os.path.join(s, "output", "search", "report.parquet"),
                            os.path.join(s, "input", "conditions.csv"), design=DESIGN)
        self.assertEqual(r["message"].count("--numeric-covariates Covariate1 --block Mouse"), 3)

    def test_the_cli_reads_the_design_file(self):
        with tempfile.TemporaryDirectory() as d:
            s = session(d, "ip")
            dj = os.path.join(d, "design.json")
            with open(dj, "w") as fh:
                json.dump(DESIGN, fh)
            p = subprocess.run([sys.executable, CHECK, "required", "--input",
                                os.path.join(s, "output", "search", "report.parquet"), "--metadata",
                                os.path.join(s, "input", "conditions.csv"), "--design-json", dj],
                               capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("--numeric-covariates Covariate1 --block Mouse", json.loads(p.stdout)["message"])


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class EndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        cls.report, cls.cond, cls.contrast = make_report(os.path.join(cls.d, "w"), "whole")
        cls.s = session(os.path.join(cls.d, "w"), "whole_proteome")
        cls.nc = os.path.join(cls.d, "nc")
        cls.check_des = {q: cls.de(q, os.path.join(cls.nc, q), "--check-input")
                         for q in (et.NORMALISED, et.RAW)}
        cls.check_run = subprocess.run([sys.executable, CHECK, "run", "--report", cls.report, "--conditions",
                                  cls.cond, "--de-normalised", os.path.join(cls.nc, et.NORMALISED),
                                  "--de-raw", os.path.join(cls.nc, et.RAW), "--session", cls.s],
                                 capture_output=True, text=True, timeout=300)
        cls.path = os.path.join(cls.nc, nc.RECORD)

    @classmethod
    def de(cls, quantities, out, *extra, contrast=None, cond=None):
        return subprocess.run(["Rscript", RUN_DE, "--input", cls.report, "--metadata", cond or cls.cond,
                               "--contrasts", contrast or cls.contrast, "--quantities", quantities,
                               "--outdir", out, *extra], capture_output=True, text=True,
                              timeout=900, cwd=cls.d)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def setUp(self):
        no_staff_list(self, self.d)
        for q, p in self.check_des.items():
            self.assertEqual(p.returncode, 0, (q, p.stderr[-2000:]))
        self.assertEqual(self.check_run.returncode, 0, self.check_run.stdout + self.check_run.stderr)

    def prov(self, out):
        with open(os.path.join(out, "de_provenance.json")) as fh:
            return json.load(fh)

    def test_run_de_records_its_design_and_conditions(self):
        p = self.prov(os.path.join(self.nc, et.RAW))
        self.assertEqual(p["design_request"]["formula"], "~ 0 + groups")
        self.assertEqual(p["design_request"]["contrasts"], [self.contrast])
        self.assertEqual(p["design_request"]["flags"], ["--contrasts", self.contrast])
        self.assertIsNone(p["design_request"]["block_column"])
        self.assertEqual(p["metadata_sha256"], nc.sha256(self.cond))
        with open(self.path) as fh:
            cd = json.load(fh)["check_design"]
        self.assertEqual((cd["status"], cd["method"], cd["contrasts"]), ("recorded", "dpc", [self.contrast]))
        with open(os.path.join(self.nc, nc.SUMMARY)) as fh:
            self.assertIn("**Model of the two DEs compared:** ~ 0 + groups (method dpc)", fh.read())

    def test_the_final_de_with_the_checks_design_matches(self):
        out = os.path.join(self.d, "final_same")
        p = self.de(et.NORMALISED, out, "--normalization-check", self.path)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        self.assertEqual(self.prov(out)["normalization_check"]["design_match"], "matches")

    def test_another_method_is_refused_and_other_contrasts_warned(self):
        p = self.de(et.NORMALISED, os.path.join(self.d, "final_ml"), "--method", "maxlfq",
                    "--normalization-check", self.path)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("REFUSED: this DE's design differs from the one the normalisation check's two "
                      "DEs were fitted with (method: 'dpc' vs 'maxlfq')", p.stderr)
        a, b = self.contrast.split("-")
        out = os.path.join(self.d, "final_rev")
        p = self.de(et.NORMALISED, out, "--normalization-check", self.path, contrast=f"{b}-{a}")
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        self.assertIn(f"the check did not look at {b}-{a}", p.stderr)
        n = self.prov(out)["normalization_check"]
        self.assertEqual(n["design_match"], "matches; contrasts differ")
        with open(os.path.join(out, "methods.txt")) as fh:
            self.assertIn("Design vs the check's DEs: matches; contrasts differ", fh.read())

    def test_a_check_from_before_2_12_runs_but_is_not_checked(self):
        with open(self.path) as fh:
            rec = json.load(fh)
        del rec["check_design"]
        old = os.path.join(self.d, "old", nc.RECORD)
        os.makedirs(os.path.dirname(old))
        with open(old, "w") as fh:
            json.dump(rec, fh)
        out = os.path.join(self.d, "final_old")
        p = self.de(et.NORMALISED, out, "--normalization-check", old)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        self.assertIn("NOT CHECKED", p.stderr)
        self.assertEqual(self.prov(out)["normalization_check"]["design_match"], "not recorded")
        with open(os.path.join(out, "methods.txt")) as fh:
            self.assertIn("Design vs the check's DEs: not recorded [confirm]", fh.read())

    def test_check_des_with_different_designs_or_conditions_are_refused(self):
        a, b = self.contrast.split("-")
        other = os.path.join(self.d, "nc_other")
        p = self.de(et.RAW, os.path.join(other, et.RAW), "--check-input", contrast=f"{b}-{a}")
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        q = subprocess.run([sys.executable, CHECK, "run", "--report", self.report, "--conditions",
                            self.cond, "--de-normalised", os.path.join(self.nc, et.NORMALISED),
                            "--de-raw", os.path.join(other, et.RAW), "--out", other],
                           capture_output=True, text=True, timeout=300)
        self.assertNotEqual(q.returncode, 0)
        self.assertIn("the two check DEs were fitted with different designs", q.stderr)
        self.assertFalse(os.path.exists(os.path.join(other, nc.RECORD)))
        # the same design on a conditions file with other content (its rows in another order)
        cond2 = os.path.join(self.d, "conditions_copy.csv")
        with open(self.cond, "rb") as fh:
            head, *rows = fh.read().splitlines(keepends=True)
        with open(cond2, "wb") as out:
            out.write(head + b"".join(reversed(rows)))
        self.assertNotEqual(nc.sha256(cond2), nc.sha256(self.cond))
        third = os.path.join(self.d, "nc_third")
        p = self.de(et.RAW, os.path.join(third, et.RAW), "--check-input", cond=cond2)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        q = subprocess.run([sys.executable, CHECK, "run", "--report", self.report, "--conditions",
                            self.cond, "--de-normalised", os.path.join(self.nc, et.NORMALISED),
                            "--de-raw", os.path.join(third, et.RAW), "--out", third],
                           capture_output=True, text=True, timeout=300)
        self.assertNotEqual(q.returncode, 0)
        self.assertIn("read another conditions than the check's --conditions", q.stderr)


if __name__ == "__main__":
    unittest.main()
