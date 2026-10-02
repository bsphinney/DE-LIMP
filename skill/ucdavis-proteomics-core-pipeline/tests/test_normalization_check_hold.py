#!/usr/bin/env python3
"""
The normalisation check cannot be skipped, is tied to its data, and never decides for a user
who said nothing (2.10 review of 40c5b9f: HIGH 1, MED 5; tests/test_review210_results.py holds
the reviewer's own cases).

  * a session that records an experiment type refuses a DE with neither --normalization-check
    nor --check-input, printing step 8's exact commands (normalization_check.required);
  * with NO type recorded the check asks (exit 3, no decision), and the audit never says PASS;
  * a DE without a decided check is an AUDIT.md WARN, and core_submission.py deliver holds the
    delivery (normalization_check.delivery_gate, the second gate beside the QC one);
  * the gate refuses another report or another conditions file than the check's (sha256);
  * reproduce.sh says plainly that a maxlfq raw DE needs a --no-norm re-quantification.

Synthetic reports and DE folders on a temp dir; the R parts skip without R. Hermetic.
"""
import csv
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
AUDIT = os.path.join(SCRIPTS, "audit_results.py")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")


def fake_de(d, contrast, rows):
    """A run_de.R output folder: de_provenance.json naming one DE table, and the table."""
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "DE_x.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein.Group", "Genes", "logFC", "P.Value", "adj.P.Val"])
        w.writerows(rows)
    with open(os.path.join(d, "de_provenance.json"), "w") as fh:
        json.dump({"contrasts": [contrast], "adjp": 0.05,
                   "de_tables": {contrast: {"file": "DE_x.csv"}}}, fh)
    return d


BALANCED = [(f"P{i}", f"G{i}", 0.02 * ((i % 7) - 3), 0.2 + (i % 50) / 62.5, 0.9)
            for i in range(200)]


@unittest.skipUnless(HAVE_ARROW, "pyarrow is needed to write the report")
class NoTypeAsks(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        no_staff_list(self, self.d)
        self.report, self.cond, self.contrast = make_report(os.path.join(self.d, "w"), "whole")
        self.dn = fake_de(os.path.join(self.d, "nc", "normalised"), self.contrast, BALANCED)
        self.dr = fake_de(os.path.join(self.d, "nc", "raw"), self.contrast, BALANCED)

    def run_check(self, s):
        return subprocess.run([sys.executable, CHECK, "run", "--report", self.report,
                               "--conditions", self.cond, "--de-normalised", self.dn,
                               "--de-raw", self.dr, "--session", s],
                              capture_output=True, text=True, timeout=300)

    def test_a_session_with_no_type_is_asked_never_auto_decided(self):
        p = self.run_check(session(os.path.join(self.d, "untyped")))
        self.assertEqual(p.returncode, nc.TRIPPED_EXIT, p.stdout + p.stderr)
        self.assertIn("NO EXPERIMENT TYPE", p.stderr)
        j = json.loads(p.stdout)
        self.assertIsNone(j["decision"])
        self.assertFalse(j["tripped"], "the data are fine -- it is the missing type that asks")
        path = j["check"]
        with open(path) as fh:
            rec = json.load(fh)
        self.assertIn("no experiment type is recorded", rec["needs_decision_why"])
        self.assertIn("no experiment type is recorded", rec["recommendation"])
        with open(os.path.join(os.path.dirname(path), nc.SUMMARY)) as fh:
            self.assertIn("**Ask the user:** no experiment type is recorded", fh.read())
        self.assertFalse(nc.gate(path, et.NORMALISED)[0])
        # once someone chooses, the audit quotes it as a WARN -- never a PASS
        rec = nc.decide(path, et.NORMALISED, "canalyst", "a whole-lysate comparison")
        tables = os.path.join(self.d, "tables")
        os.makedirs(tables)
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump({"normalization_check": dict(rec, status="decided",
                                                   quantities_applied=et.NORMALISED)}, fh)
        a = subprocess.run([sys.executable, AUDIT, "--out", "AUDIT.md", "--de-dir", tables],
                           capture_output=True, text=True, timeout=120, cwd=self.d)
        self.assertEqual(a.returncode, 0, a.stderr)
        with open(os.path.join(self.d, "AUDIT.json")) as fh:
            f = [x for x in json.load(fh)["findings"] if x["check"] == "normalization_check"]
        self.assertEqual([x["status"] for x in f], ["WARN"])
        self.assertIn("No experiment type was recorded", f[0]["message"])

    def test_the_gate_refuses_other_conditions_than_the_checks(self):
        s = session(os.path.join(self.d, "typed"), "whole_proteome")
        p = self.run_check(s)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        path = json.loads(p.stdout)["check"]
        with open(path) as fh:
            rec = json.load(fh)
        self.assertEqual(rec["report_sha256"], nc.sha256(self.report))
        ok, msg, _ = nc.gate(path, et.NORMALISED, self.report, self.cond)
        self.assertTrue(ok, msg)
        edited = os.path.join(self.d, "conditions_v2.csv")
        with open(self.cond) as fh:
            text = fh.read()
        with open(edited, "w") as fh:
            fh.write(text.replace("B_3,B", "B_3,A"))
        ok, msg, _ = nc.gate(path, et.NORMALISED, self.report, edited)
        self.assertFalse(ok)
        self.assertIn("never authorises another", msg)


class StepEightIsRequired(unittest.TestCase):
    def test_a_typed_session_requires_the_check_and_names_its_commands(self):
        with tempfile.TemporaryDirectory() as d:
            s = session(d, "ip")
            meta = os.path.join(s, "input", "conditions.csv")
            report = os.path.join(s, "output", "search", "report.parquet")
            r = nc.required(report, meta)
            self.assertTrue(r["required"])
            self.assertEqual((r["type"], r["default"]), ("ip", et.RAW))
            msg = r["message"]
            self.assertIn("REFUSED", msg)
            self.assertEqual(msg.count("--check-input"), 3)        # two commands + the note
            self.assertIn("normalization_check.py", msg)
            self.assertIn("--normalization-check", msg)
            self.assertIn(os.path.join(s, "output", "norm_check", "raw"), msg)
            # found from the report as well as from the metadata; nothing found elsewhere
            self.assertEqual(et.session_of(report), s)
            self.assertFalse(nc.required(os.path.join(d, "x.parquet"),
                                         os.path.join(d, "c.csv"))["required"])

    def test_the_chains_recovery_step_is_the_check(self):
        nxt = nc.step8_commands("/s/output/search/report.parquet", "/s/input/conditions.csv",
                                "/s")["next"]
        first, rest = nxt.split(" && ", 1)
        self.assertIn("--check-input", first)
        self.assertIn("normalization_check.py", rest)
        self.assertIn("exit 3: stop", nxt)


class DeliveryIsHeld(unittest.TestCase):
    def gate(self, prov):
        with tempfile.TemporaryDirectory() as d:
            t = os.path.join(d, "output", "tables")
            os.makedirs(t)
            if prov is not None:
                with open(os.path.join(t, "de_provenance.json"), "w") as fh:
                    json.dump(prov, fh)
            return nc.delivery_gate(d)

    def test_only_a_decided_check_proceeds(self):
        decided = {"normalization_check": {"status": "decided",
                                           "decision": {"quantities": "raw", "by": "x"}}}
        self.assertTrue(self.gate(decided)["proceed"])
        for prov, words in (({"normalization_check": {"status": "not_run"}}, "without"),
                            ({"normalization_check": {"status": "check_input"}}, "check's two inputs"),
                            ({"method": "dpc"}, "predates"),
                            (None, "could not be read")):
            g = self.gate(prov)
            self.assertFalse(g["proceed"], prov)
            self.assertIn(words, g["reason"])
            self.assertIn("step 8", g["hint"])

    def test_core_submission_deliver_holds_an_unchecked_de(self):
        from test_core_submission import DeliverBase

        class Held(DeliverBase):
            def runTest(self):
                pass
        t = Held()
        t.setUp()
        try:
            with open(os.path.join(t.output, "tables", "de_provenance.json"), "w") as fh:
                json.dump({"normalization_check": {"status": "not_run",
                                                   "note": "no check"}}, fh)
            rc, out, p = t.deliver()
            self.assertNotEqual(rc, 0, p.stdout)
            self.assertIn("held: the DE ran without the normalisation check", out["error"])
            self.assertFalse(out["normalization_gate"]["proceed"])
            self.assertFalse(os.path.exists(t.delivery), "nothing may be delivered")
        finally:
            t.tearDown()


class ReproduceSaysNoNorm(unittest.TestCase):
    def test_a_maxlfq_raw_de_needs_a_no_norm_requantification(self):
        with tempfile.TemporaryDirectory() as d:
            de = os.path.join(d, "tables")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump({"method": "maxlfq", "quantities": "raw"}, fh)
            out = os.path.join(d, "repro")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir",
                                out, "--de-dir", de, "--de-method", "maxlfq"],
                               capture_output=True, text=True, timeout=300)
            self.assertEqual(p.returncode, 0, p.stderr[-1500:])
            with open(os.path.join(out, "reproduce.sh")) as fh:
                sh = fh.read()
        self.assertIn("the DE read a --no-norm report's PG.MaxLFQ", sh)
        self.assertIn("diann_parallel.py\n#    --no-norm", sh)
        self.assertIn("--quantities raw", sh)


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class RunDeRefusesAnUncheckedDe(unittest.TestCase):
    def test_refused_with_the_commands_and_check_input_allowed(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        report, cond, contrast = make_report(os.path.join(d, "ip"), "ip")
        s = session(os.path.join(d, "ip"), "ip", bait="BAITX")
        meta = os.path.join(s, "input", "conditions.csv")
        shutil.copy(cond, meta)

        def de(*extra):
            return subprocess.run(["Rscript", RUN_DE, "--input", report, "--metadata", meta,
                                   "--contrasts", contrast, *extra], capture_output=True,
                                  text=True, timeout=900, cwd=d)
        r = de("--outdir", os.path.join(s, "output", "tables"))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("REFUSED: this session records the experiment type", r.stderr)
        self.assertIn("Run step 8's check first", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(s, "output", "tables", "de_provenance.json")))
        out = os.path.join(s, "output", "norm_check", "raw")
        r = de("--quantities", "raw", "--check-input", "--outdir", out)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            self.assertEqual(json.load(fh)["normalization_check"]["status"], "check_input")
        with open(os.path.join(out, "methods.txt")) as fh:
            self.assertIn("CHECK INPUT, NOT A FINAL DE", fh.read())


if __name__ == "__main__":
    unittest.main()
