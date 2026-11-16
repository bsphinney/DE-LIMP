#!/usr/bin/env python3
"""
session.py status and finalize tell the same truth about logs/commands.log, and print a bundle
command that works (2.15 acceptance, lane X4).

  * logs/commands.log missing: `status` carries the same warning finalize prints (the one
    definition, session.commands_log_warning) as a note -- WARN, never a missing finish item, so
    the exit code and "Complete" are as before;
  * logs/commands.log written: neither says it;
  * `session.py log-command` is how an agent writes it (nothing in the skill's own scripts does
    when they run over plain ssh): one timestamped line per command, appended;
  * the reproducibility command `status` and finalize's gate print carries `--session <S>`: run
    without it, provenance.py built a bundle that claimed no DE ran (pxd4-215).
"""
import json
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import report_files as rf  # noqa: E402
import session as sm  # noqa: E402
import test_finish_gates as tfg  # noqa: E402

PY = sys.executable
SESSION_PY = os.path.join(SCRIPTS, "session.py")


class Base(tfg.FinalizeBase):
    def status(self, *extra):
        r = subprocess.run([PY, SESSION_PY, "status", "--session", self.sd, *extra],
                           capture_output=True, text=True, env=self.env)
        return r

    def log_path(self):
        return self.p["commands_log"]


class StatusAndFinalizeAgree(Base):
    def test_a_missing_commands_log_is_a_status_note_and_not_a_missing_item(self):
        tfg.skill_report(self.out)
        self.assertEqual(self.finalize().returncode, 0)
        self.assertFalse(os.path.exists(self.log_path()))
        r = self.status()
        want = sm.commands_log_warning(self.p)
        self.assertTrue(want and "commands.log" in want and "missing" in want, want)
        self.assertIn("note     " + want, r.stdout)
        self.assertEqual(r.returncode, sm.STATUS_COMPLETE, r.stdout)     # a WARN, not a gap
        self.assertIn("Complete: every item is done.", r.stdout)
        j = json.loads(self.status("--json").stdout)
        self.assertIn(want, j["notes"])
        self.assertTrue(j["complete"])
        # finalize said the same thing, in the stamp
        stamp = json.load(open(self.stamp_path()))
        self.assertIn(want, stamp["warnings"])

    def test_a_written_commands_log_is_not_a_note(self):
        tfg.skill_report(self.out)
        r = subprocess.run([PY, SESSION_PY, "log-command", "--session", self.sd,
                            "python3 scripts/run_search.py --engine diann"],
                           capture_output=True, text=True, env=self.env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.finalize().returncode, 0)
        r = self.status()
        self.assertNotIn("commands.log", r.stdout)

    def test_log_command_appends_one_line_per_command(self):
        for c in ("first command", "second command"):
            r = subprocess.run([PY, SESSION_PY, "log-command", "--session", self.sd, c],
                               capture_output=True, text=True, env=self.env)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        lines = open(self.log_path()).read().splitlines()
        self.assertEqual(len(lines), 2, lines)
        self.assertTrue(lines[0].endswith("first command") and lines[1].endswith("second command"))
        self.assertRegex(lines[0], r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ ")
        self.assertEqual(sm.commands_log_warning(self.p) is None or
                         "older" in sm.commands_log_warning(self.p), True)


class TheBundleCommandNamesTheSession(Base):
    def test_status_prints_provenance_with_session(self):
        tfg.skill_report(self.out)
        cmd = rf._bundle_command(self.sd)
        self.assertIsNotNone(cmd)
        self.assertIn("--session " + self.sd, cmd, cmd)
        self.assertTrue(cmd.startswith("python3 scripts/provenance.py --outdir "), cmd)

    def test_the_finish_gate_prints_it_too(self):
        gap = rf.finish_gaps(self.out, ("reproducibility",))[0]
        self.assertIn("--session " + self.sd, gap, gap)


if __name__ == "__main__":
    unittest.main()
