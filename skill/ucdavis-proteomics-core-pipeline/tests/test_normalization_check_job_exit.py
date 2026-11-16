#!/usr/bin/env python3
"""
A normalisation check that stops to ask is not a failed Slurm job (2.15 acceptance, lane X4:
three runs saw `sacct` say FAILED for a check that had done its work, indistinguishable from a
crash for a watcher).

  * `normalization_check.py run` still exits 3 for a direct caller (TRIPPED_EXIT: stop and ask);
  * `run --job-exit0` -- what a job script uses -- exits 0 and leaves NEEDS_DECISION.json beside
    the check record, naming the exit code it would have had; a check that needs nothing leaves no
    marker, and a stale one is removed;
  * `decide` takes the marker away: the question is answered;
  * session.py status says a decision is owed while the marker is there (a note, not a missing
    finish item).

Hermetic: synthetic reports and DE folders on a temp dir; pyarrow needed for the report.
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
from test_normalization_check import HAVE_ARROW, make_report, no_staff_list, session  # noqa: E402
from test_normalization_check_hold import BALANCED, CHECK, fake_de  # noqa: E402

MARKER = "NEEDS_DECISION.json"


@unittest.skipUnless(HAVE_ARROW, "pyarrow is needed to write the report")
class JobExit(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        no_staff_list(self, self.d)
        self.report, self.cond, self.contrast = make_report(os.path.join(self.d, "w"), "whole")
        self.dn = fake_de(os.path.join(self.d, "nc", "normalised"), self.contrast, BALANCED)
        self.dr = fake_de(os.path.join(self.d, "nc", "raw"), self.contrast, BALANCED)

    def run_check(self, s, *extra):
        return subprocess.run([sys.executable, CHECK, "run", "--report", self.report,
                               "--conditions", self.cond, "--de-normalised", self.dn,
                               "--de-raw", self.dr, "--session", s, *extra],
                              capture_output=True, text=True, timeout=300)

    def untyped(self):
        return session(os.path.join(self.d, "untyped"))        # no experiment type: asks

    def test_direct_callers_still_get_exit_3(self):
        p = self.run_check(self.untyped())
        self.assertEqual(p.returncode, nc.TRIPPED_EXIT, p.stdout + p.stderr)

    def test_a_job_exits_0_and_leaves_the_marker(self):
        s = self.untyped()
        p = self.run_check(s, "--job-exit0")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("NO EXPERIMENT TYPE", p.stderr)          # still said, in the job's log
        marker = os.path.join(s, "output", "norm_check", MARKER)
        with open(marker) as fh:
            m = json.load(fh)
        self.assertTrue(m["needs_decision"])
        self.assertEqual(m["exit_code"], nc.TRIPPED_EXIT)
        self.assertIn("no experiment type", m["why"])
        self.assertTrue(m["summary"].endswith(nc.SUMMARY))
        self.assertTrue(os.path.isfile(m["check"]))
        # the check json still says a decision is owed, and the gate still refuses
        self.assertFalse(nc.gate(m["check"], et.NORMALISED)[0])

    def test_the_marker_is_written_for_a_direct_caller_too(self):
        s = self.untyped()
        self.assertEqual(self.run_check(s).returncode, nc.TRIPPED_EXIT)
        self.assertTrue(os.path.isfile(os.path.join(s, "output", "norm_check", MARKER)))

    def test_a_check_that_needs_nothing_leaves_no_marker_and_clears_a_stale_one(self):
        s = session(os.path.join(self.d, "typed"), "whole_proteome")
        out = os.path.join(s, "output", "norm_check")
        os.makedirs(out)
        with open(os.path.join(out, MARKER), "w") as fh:
            fh.write("{}")                                      # left by an earlier run
        p = self.run_check(s, "--job-exit0")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertFalse(os.path.exists(os.path.join(out, MARKER)))

    def test_deciding_takes_the_marker_away(self):
        s = self.untyped()
        p = self.run_check(s, "--job-exit0")
        marker = os.path.join(s, "output", "norm_check", MARKER)      # where status looks
        self.assertTrue(os.path.isfile(marker))
        nc.decide(json.loads(p.stdout)["check"], et.NORMALISED, "canalyst",
                  "a whole-lysate comparison")
        self.assertFalse(os.path.exists(marker))


class StatusSaysADecisionIsOwed(unittest.TestCase):
    def test_status_notes_the_marker(self):
        import session as sm
        with tempfile.TemporaryDirectory() as d:
            sd = os.path.join(d, "S")
            nd = os.path.join(sd, "output", "norm_check")
            os.makedirs(nd)
            os.makedirs(os.path.join(sd, "input"))
            with open(os.path.join(nd, MARKER), "w") as fh:
                json.dump({"needs_decision": True, "exit_code": 3, "why": "it tripped",
                           "summary": os.path.join(nd, nc.SUMMARY)}, fh)
            notes = sm.status_notes(sd)
        self.assertTrue(any("normalisation decision is owed" in n and "NORMALIZATION_CHECK.md" in n
                            for n in notes), notes)


class StatusSaysASageJobFailed(unittest.TestCase):
    def test_status_notes_sage_failure_json(self):
        import session as sm
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "output", "search"))
            with open(os.path.join(d, "output", "search", "sage_failure.json"), "w") as fh:
                json.dump({"failure": "sage_stack_overflow", "exit_code": 134,
                           "hint": "re-run this job once"}, fh)
            notes = sm.status_notes(d)
        self.assertTrue(any("sage_stack_overflow" in n and "re-run this job once" in n
                            for n in notes), notes)


if __name__ == "__main__":
    unittest.main()
