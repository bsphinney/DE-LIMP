#!/usr/bin/env python3
"""
A Sage LFQ window that does not fit the runs' MS1 mass error is a GATE with a way through, not
only a warning (SKILL_OPEN_DEFECTS, 2.9 backlog; staff report 2026-09-29: +7.1/+7.7 ppm runs
against Sage's 5 ppm window kept 0 MS1 peaks, and the quantities went on into the DE).

  * the corrected window is computed from the measured offsets and written into a copy of the
    config the search ran with (<out>/sage_config.lfq_ppm<N>.json);
  * the exact run_search.py command that repeats the search with it, into <out>_lfq_ppm<N>, is
    recorded in sage_lfq_check.json and search_provenance.json, and printed by the job;
  * run_search.py --adapt-only and an inline search refuse to build report.parquet, and run_de.R
    refuses one built before the gate, until the re-run is used or the user accepts the
    quantities (--accept-lfq-window "<who, why>", recorded, and named in AUDIT.md);
  * running that command, and its job, gives a search the check passes.

Sage, the converters and SLURM are stand-ins (tests/test_sage_sbatch_conversion.py's harness;
test_sage_lfq_check.py's Sage-shaped fixtures). Hermetic.
"""
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import sage_lfq_check as slc  # noqa: E402
import test_sage_lfq_check as fx  # noqa: E402
from job_env import job_env  # noqa: E402,F401  (the jobs below run under the harness's job_env)
from test_sage_sbatch_conversion import _Harness  # noqa: E402
from test_run_de_contaminants import r_has  # noqa: E402

AUDIT = os.path.join(SCRIPTS, "audit_results.py")
OFF = {"HeLa_1.mzML": 7.1, "HeLa_2.mzML": 7.7}
REASON = "core analyst: identifications only, these quantities are not reported"


@unittest.skipUnless(fx.HAVE_ARROW, "pyarrow is needed to write the Sage fixtures")
class TheCheckComputesTheFix(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.out = os.path.join(self.d, "search")
        self.cfg = os.path.join(self.d, "sage_config.json")
        with open(self.cfg, "w") as fh:
            json.dump({"precursor_tol": {"ppm": [-10.0, 10.0]},
                       "quant": {"lfq": True, "lfq_settings": {"peak_scoring": "Hybrid"}}}, fh)

    def test_the_corrected_config_is_the_run_config_with_the_measured_window(self):
        fx.write_sage_outputs(self.out, OFF, peaks=0)
        rec = slc.check(self.out, self.cfg)
        self.assertEqual(rec["suggested_ppm_tolerance"], 10)                 # ceil(7.7 + 2)
        self.assertEqual(rec["corrected_params"],
                         os.path.join(self.out, "sage_config.lfq_ppm10.json"))
        with open(rec["corrected_params"]) as fh:
            cfg = json.load(fh)
        self.assertEqual(cfg["quant"]["lfq_settings"],
                         {"peak_scoring": "Hybrid", "ppm_tolerance": 10.0})
        self.assertEqual(cfg["precursor_tol"], {"ppm": [-10.0, 10.0]}, "nothing else changes")
        self.assertTrue(slc.gated(rec))

    def test_without_a_config_it_says_why(self):
        fx.write_sage_outputs(self.out, OFF, peaks=0)
        rec = slc.check(self.out)
        self.assertIsNone(rec["corrected_params"])
        self.assertIn("not known", rec["corrected_params_why"])
        self.assertIn("not known", slc.refusal(self.out, rec))

    def test_a_low_peak_count_with_a_fitting_window_is_not_gated(self):
        """Its cause is unknown, so there is no corrected run to send the user to."""
        fx.write_sage_outputs(self.out, {"a.mzML": 0.2}, peaks=0)
        rec = slc.check(self.out, self.cfg)
        self.assertEqual(rec["status"], "warn")
        self.assertFalse(slc.gated(rec))
        self.assertIsNone(slc.refusal(self.out, rec))
        self.assertNotIn("corrected_params", rec)

    def test_the_rerun_command_is_built_from_the_record(self):
        fx.write_sage_outputs(self.out, OFF, peaks=0)
        raws = ["/data/raw/HeLa_1.raw", "/data/raw/HeLa 2.raw"]
        with open(os.path.join(self.out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "sage", "tools": "/t/tools.json", "bundle": "/w/wf.json",
                       "fasta": "/db/search.fasta", "files": raws, "threads": 16,
                       "keratin_sample": True, "submitted_sbatch": "/s/sage_job.sh",
                       "queue": {"partition": "low", "account": "publicgrp", "qos": None}}, fh)
        rec = slc.check(self.out, self.cfg)
        msg = slc.refusal(self.out, rec)
        plan = rec["rerun"]
        self.assertEqual(plan["out"], self.out + "_lfq_ppm10")
        argv = shlex.split(plan["command"])
        self.assertEqual(argv[1], os.path.join(SCRIPTS, "run_search.py"))

        def opt(flag):
            return argv[argv.index(flag) + 1]
        self.assertEqual((opt("--params"), opt("--out"), opt("--tools"), opt("--bundle"),
                          opt("--fasta"), opt("--threads"), opt("--engine"), opt("--sbatch"),
                          opt("--partition"), opt("--account")),
                         (rec["corrected_params"], self.out + "_lfq_ppm10", "/t/tools.json",
                          "/w/wf.json", "/db/search.fasta", "16", "sage", "/s/sage_job_lfq_ppm10.sh",
                          "low", "publicgrp"))
        self.assertNotIn("--qos", argv)
        self.assertIn("--keratin-sample", argv)
        self.assertEqual(argv[argv.index("--files") + 1:], raws)
        self.assertIn(plan["command"], msg)
        self.assertIn("+/-10 ppm (the worst run's median 7.7 ppm + 2 ppm, rounded up)", msg)
        self.assertIn("--accept-lfq-window", msg)
        # recorded: the record and search_provenance.json carry the plan and the refusal
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)["sage_lfq_check"]
        self.assertEqual(prov["rerun"]["command"], plan["command"])
        self.assertEqual(prov["refusal"], msg)

    def test_a_record_without_what_the_command_needs_says_so(self):
        fx.write_sage_outputs(self.out, OFF, peaks=0)
        with open(os.path.join(self.out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "sage", "fasta": "/db/search.fasta"}, fh)
        rec = slc.check(self.out, self.cfg)
        slc.refusal(self.out, rec)
        self.assertIsNone(rec["rerun"]["command"])
        self.assertIn("not recorded: tools, bundle, files", rec["rerun"]["why"])
        self.assertIn(f"--params {rec['corrected_params']}", rec["rerun"]["why"])

    def test_acceptance_is_recorded_kept_on_a_recheck_and_named_in_the_audit(self):
        fx.write_sage_outputs(self.out, OFF, peaks=0)
        rec = slc.check(self.out, self.cfg)
        slc.accept(self.out, rec, REASON, by="analyst1")
        self.assertFalse(slc.gated(rec))
        again = slc.carry_acceptance(slc.load(self.out), slc.check(self.out, self.cfg))
        self.assertEqual(again["accepted"]["reason"], REASON)
        # what was accepted changed (another window): the acceptance does not carry
        fx.write_sage_outputs(self.out, OFF, peaks=0, lfq_settings={"ppm_tolerance": 7.0})
        changed = slc.carry_acceptance(slc.load(self.out), slc.check(self.out, self.cfg))
        self.assertNotIn("accepted", changed)
        status, msg, detail = slc.audit_finding(again)
        self.assertEqual(status, "WARN")
        self.assertIn(f"The user accepted these quantities as they are (analyst1, ", msg)
        self.assertIn(REASON, msg)
        with self.assertRaises(ValueError):
            slc.accept(self.out, rec, "  ")

    def test_the_gate_cli_answers_run_de(self):
        fx.write_sage_outputs(self.out, OFF, peaks=0)
        slc.record(self.out, slc.check(self.out, self.cfg))
        gate = [sys.executable, os.path.join(SCRIPTS, "sage_lfq_check.py"), "gate", "--out",
                self.out]
        j = json.loads(subprocess.run(gate, capture_output=True, text=True, timeout=60).stdout)
        self.assertEqual((j["checked"], j["refused"]), (True, True))
        self.assertIn("REFUSED", j["message"])
        rec = slc.load(self.out)
        slc.accept(self.out, rec, REASON)
        j = json.loads(subprocess.run(gate, capture_output=True, text=True, timeout=60).stdout)
        self.assertEqual((j["checked"], j["refused"], j["message"]), (True, False, None))


@unittest.skipUnless(fx.HAVE_ARROW, "pyarrow is needed to write the Sage fixtures")
class TheGateEndToEnd(_Harness):
    def adapt(self, out, *extra):
        return subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--tools", self.tools,
             "--bundle", self.bundle, "--params", self.cfg, "--fasta", self.fasta, "--out", out,
             "--files", *self.raws, "--engine", "sage", "--adapt-only", *extra],
            capture_output=True, text=True, timeout=180, env=self.env(), cwd=self.d)

    def test_job_refusal_rerun_and_acceptance(self):
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        r = subprocess.run(["bash", job], capture_output=True, text=True, timeout=180,
                           env=self.env(FAKE_SAGE_OFFSET="7.4"), cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)        # Sage did its part
        self.assertIn("[sage_lfq_check] REFUSED", r.stderr)            # ... and the job says why
        rec = slc.load(self.out)
        command = rec["rerun"]["command"]
        self.assertIsNotNone(command, rec["rerun"])
        self.assertIn(command, r.stderr)

        # the adapt step refuses: no report.parquet from these quantities
        a = self.adapt(self.out)
        self.assertNotEqual(a.returncode, 0)
        self.assertIn("REFUSED: report.parquet is not built from this Sage search", a.stderr)
        self.assertIn(command, a.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "report.parquet")))

        # the recorded command, as printed, repeats the search with the corrected window
        argv = shlex.split(command)
        g = subprocess.run([sys.executable] + argv[1:], capture_output=True, text=True,
                           timeout=180, env=self.env(), cwd=self.d)
        self.assertEqual(g.returncode, 0, g.stdout + g.stderr)
        new_out = argv[argv.index("--out") + 1]
        new_job = argv[argv.index("--sbatch") + 1]
        r2 = subprocess.run(["bash", new_job], capture_output=True, text=True, timeout=180,
                            env=self.env(FAKE_SAGE_OFFSET="7.4", FAKE_SAGE_PEAKS="380"),
                            cwd=self.d)
        self.assertEqual(r2.returncode, 0, r2.stdout + r2.stderr)
        rerun = slc.load(new_out)
        self.assertEqual((rerun["status"], rerun["ppm_tolerance"]), ("ok", 10.0), rerun)
        a2 = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--tools", self.tools,
             "--bundle", self.bundle, "--params", argv[argv.index("--params") + 1],
             "--fasta", self.fasta, "--out", new_out, "--files", *self.raws, "--engine", "sage",
             "--adapt-only"], capture_output=True, text=True, timeout=180, env=self.env(),
            cwd=self.d)
        self.assertEqual(a2.returncode, 0, a2.stderr)
        self.assertTrue(os.path.isfile(os.path.join(new_out, "report.parquet")))

        # or the user accepts the first search's quantities as they are: recorded, and adapted
        a3 = self.adapt(self.out, "--accept-lfq-window", REASON)
        self.assertEqual(a3.returncode, 0, a3.stderr)
        self.assertTrue(os.path.isfile(os.path.join(self.out, "report.parquet")))
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            self.assertEqual(json.load(fh)["sage_lfq_check"]["accepted"]["reason"], REASON)
        au = subprocess.run([sys.executable, AUDIT, "--out", "AUDIT.md", "--search-out",
                             self.out], capture_output=True, text=True, timeout=120, cwd=self.d)
        self.assertEqual(au.returncode, 0, au.stderr)
        with open(os.path.join(self.d, "AUDIT.md")) as fh:
            md = fh.read()
        self.assertIn("The user accepted these quantities as they are", md)
        self.assertIn(REASON, md)

    def test_an_inline_search_records_itself_and_refuses(self):
        p = self.run_search(env=self.env(FAKE_SAGE_OFFSET="7.4"))
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("REFUSED: report.parquet is not built from this Sage search", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "report.parquet")))
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["sage_lfq_check"]["status"], "warn")
        self.assertIn("--params", prov["sage_lfq_check"]["rerun"]["command"])
        self.assertIn(prov["sage_lfq_check"]["rerun"]["command"], p.stderr)

    @unittest.skipUnless(r_has("jsonlite"), "needs R with jsonlite")
    def test_run_de_refuses_a_report_built_before_the_gate(self):
        """A report.parquet adapted by 2.9 sits beside a check that now refuses: the DE stops."""
        fx.write_sage_outputs(self.out, {f"S{i}.mzML": 7.4 for i in range(1, 5)}, peaks=0)
        slc.record(self.out, slc.check(self.out, self.cfg))
        with open(os.path.join(self.out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "sage"}, fh)
        report = os.path.join(self.out, "report.parquet")
        open(report, "w").close()                         # the gate stops before reading it
        cond = os.path.join(self.d, "conditions.csv")
        with open(cond, "w") as fh:
            fh.write("File.Name,Group\nS1,A\nS2,A\nS3,B\nS4,B\n")
        # R and this python3 on PATH: run_de.R asks the gate of sage_lfq_check.py through it
        env = self.env()
        env["PATH"] = os.pathsep.join([os.path.dirname(sys.executable),
                                       os.path.dirname(shutil.which("Rscript")), "/usr/bin", "/bin"])
        r = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"), "--input", report,
                            "--metadata", cond, "--method", "maxlfq",
                            "--outdir", os.path.join(self.d, "de")],
                           capture_output=True, text=True, timeout=300, cwd=self.d, env=env)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("REFUSED: report.parquet is not built from this Sage search", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.d, "de", "de_provenance.json")))


if __name__ == "__main__":
    unittest.main()
