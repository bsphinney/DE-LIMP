#!/usr/bin/env python3
"""finish_stop_hook.py -- the plugin's Stop hook (hooks/hooks.json): blocks an agent's stop ONCE
while a skill analysis it ran has finish items missing, and is silent in every other case.

Guards (stdlib; synthetic sessions in a temp folder, the hook run as Claude Code runs it: a
subprocess reading JSON on stdin):
  * no recorded analysis for this Claude session: silent allow (exit 0, no output);
  * stop_hook_active true: allow (never loops);
  * results + missing items: a block whose reason names each missing item;
  * the same missing set is not reminded twice; a changed set is;
  * a finished session, a held session, a session with no results yet, a HIVE-only record: allow;
  * malformed / empty / non-object stdin: allow, never a failed turn;
  * hooks.json registers the Stop hook through ${CLAUDE_PLUGIN_ROOT}.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import finish_fixture  # noqa: E402
from job_env import job_env  # noqa: E402

HOOK = os.path.join(SCRIPTS, "finish_stop_hook.py")
SID = "aaaaaaaa-1111-4222-8333-bbbbbbbbbbbb"


class StopHook(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = self._td.name
        self.cfg = os.path.join(self.tmp, "skill_config")
        os.makedirs(self.cfg)
        self.sdir = os.path.join(self.tmp, "analysis_session")
        os.makedirs(os.path.join(self.sdir, "output", "tables"))
        os.makedirs(os.path.join(self.sdir, "logs"))
        self.addCleanup(self._td.cleanup)

    def record(self, sdir=None, hive=None, sid=SID):
        key = f"hive:{hive}" if hive else f"local:{sdir}"
        m = {}
        if os.path.isfile(os.path.join(self.cfg, "transcript_sessions.json")):
            with open(os.path.join(self.cfg, "transcript_sessions.json")) as fh:
                m = json.load(fh)
        m[sid] = {"sessions": {key: {"session_dir": sdir, "hive": hive, "recorded": "t"}},
                  "segments": [{"session": key, "start": {"line": 0, "after_uuid": None}}]}
        with open(os.path.join(self.cfg, "transcript_sessions.json"), "w") as fh:
            json.dump(m, fh)

    def results(self):
        with open(os.path.join(self.sdir, "output", "tables", "de.csv"), "w") as fh:
            fh.write("protein,logFC\n")

    def run_hook(self, payload=None, stdin=None):
        data = stdin if stdin is not None else json.dumps(
            {"session_id": SID, "hook_event_name": "Stop", "stop_hook_active": False,
             "cwd": self.tmp, **(payload or {})})
        r = subprocess.run([sys.executable, HOOK], input=data, capture_output=True, text=True,
                           env=job_env(self.tmp), cwd=self.tmp, timeout=60)
        return r

    def assertAllowed(self, r):
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(r.stdout, "", r.stdout)

    def assertBlocked(self, r):
        self.assertEqual(r.returncode, 0, r.stderr)
        out = json.loads(r.stdout)
        # additionalContext continues the agent WITHOUT the user's "Stop hook error occurred"
        # (live-verified); decision:block would show that error
        self.assertNotIn("decision", out)
        hso = out["hookSpecificOutput"]
        self.assertEqual(hso["hookEventName"], "Stop")
        return hso["additionalContext"]

    def test_no_recorded_session_is_silent(self):
        self.results()
        self.assertAllowed(self.run_hook())

    def test_stop_hook_active_allows(self):
        self.results()
        self.record(self.sdir)
        self.assertAllowed(self.run_hook({"stop_hook_active": True}))

    def test_missing_items_block_and_are_named(self):
        self.results()
        self.record(self.sdir)
        why = self.assertBlocked(self.run_hook())
        for label in ("Methods", "reproducibility", "AGENTS.md", "finalize stamp"):
            self.assertIn(label.lower(), why.lower())
        self.assertIn(self.sdir, why)
        self.assertIn("session.py hold", why)

    def test_reminded_once_per_set(self):
        self.results()
        self.record(self.sdir)
        self.assertBlocked(self.run_hook())
        self.assertAllowed(self.run_hook())               # same set: not again
        finish_fixture.finished(os.path.join(self.sdir, "output"))   # the set changes
        self.assertBlocked(self.run_hook())

    def test_finished_session_allows(self):
        # a session with every item done (finish_status's own test is test_finish_gates.py)
        self.results()
        self.record(self.sdir)
        import session
        import finish_stop_hook
        from unittest import mock
        done = lambda d: [{"item": "report", "label": "Report", "done": True,   # noqa: E731
                           "state": "skill", "detail": "ok", "command": None}]
        with mock.patch.dict(os.environ, {"SKILL_CONFIG_DIR": self.cfg}), \
                mock.patch.object(session, "finish_status", done):
            self.assertIsNone(finish_stop_hook.reminder({"session_id": SID}))
        with mock.patch.dict(os.environ, {"SKILL_CONFIG_DIR": self.cfg}):
            self.assertIsNotNone(finish_stop_hook.reminder({"session_id": SID}))

    def test_held_session_allows(self):
        self.results()
        self.record(self.sdir)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "hold",
                            "--session", self.sdir, "--why", "the user said stop"],
                           capture_output=True, text=True, env=job_env(self.tmp))
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertAllowed(self.run_hook())
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "hold", "--session",
                        self.sdir, "--release"], capture_output=True, env=job_env(self.tmp))
        self.assertBlocked(self.run_hook())

    def test_no_results_yet_allows(self):
        self.record(self.sdir)                            # output/tables is empty
        self.assertAllowed(self.run_hook())

    def _only(self, relpath):
        """a session whose ONLY output is this file (no output/tables content)"""
        full = os.path.join(self.sdir, "output", relpath)
        os.makedirs(os.path.dirname(full), exist_ok=True)
        with open(full, "w") as fh:
            fh.write("x")
        self.record(self.sdir)

    def test_fetched_search_report_counts_as_results(self):
        # hive_remote: `hive_exec.sh --get '<out>/report.parquet' ./<session>/output/`
        for rel in ("report.parquet", os.path.join("out", "report.parquet"),
                    os.path.join("a", "b", "report.tsv")):
            with self.subTest(rel=rel):
                self.setUp()
                self._only(rel)
                self.assertBlocked(self.run_hook())

    def test_search_branches_count_as_results(self):
        # a search record in a non-default folder (output/search_dda/), and a DE table
        self._only(os.path.join("search_dda", "search_provenance.json"))
        self.assertBlocked(self.run_hook())

    def test_unrelated_file_is_not_results(self):
        self._only(os.path.join("notes", "scratch.txt"))
        self.assertAllowed(self.run_hook())

    def test_finalized_stamp_counts_as_results(self):
        self._only("FINALIZED.json")
        self.assertBlocked(self.run_hook())

    def test_corrupt_reminder_state_is_ignored(self):
        self.results()
        self.record(self.sdir)
        with open(os.path.join(self.cfg, "finish_reminders.json"), "w") as fh:
            fh.write("{broken")
        self.assertBlocked(self.run_hook())

    def test_reminder_state_write_is_atomic(self):
        self.results()
        self.record(self.sdir)
        self.assertBlocked(self.run_hook())
        left = [n for n in os.listdir(self.cfg) if n.endswith(".tmp")]
        self.assertEqual(left, [])
        with open(os.path.join(self.cfg, "finish_reminders.json")) as fh:
            self.assertIn(SID, json.load(fh))

    def test_system_exit_in_the_path_allows(self):
        # a sys.exit(2) anywhere would be a block + `|| python` double run: caught, exit 0
        self.results()
        self.record(self.sdir)
        code = ("import sys; sys.path.insert(0, %r); import session; import finish_stop_hook;"
                "session.finish_status = lambda d: sys.exit(2);"
                "sys.exit(finish_stop_hook.main())" % SCRIPTS)
        r = subprocess.run([sys.executable, "-c", code], input=json.dumps(
            {"session_id": SID, "stop_hook_active": False}), capture_output=True, text=True,
            env=job_env(self.tmp), cwd=self.tmp, timeout=60)
        self.assertAllowed(r)

    def test_keyboard_interrupt_is_not_swallowed(self):
        self.results()
        self.record(self.sdir)
        code = ("import sys; sys.path.insert(0, %r); import session; import finish_stop_hook;"
                "session.finish_status = lambda d: (_ for _ in ()).throw(KeyboardInterrupt());"
                "sys.exit(finish_stop_hook.main())" % SCRIPTS)
        r = subprocess.run([sys.executable, "-c", code], input=json.dumps(
            {"session_id": SID, "stop_hook_active": False}), capture_output=True, text=True,
            env=job_env(self.tmp), cwd=self.tmp, timeout=60)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("KeyboardInterrupt", r.stderr)

    @unittest.skipIf(os.name == "nt", "SIGALRM is POSIX only")
    def test_time_budget_ends_silently(self):
        self.results()
        self.record(self.sdir)
        code = ("import sys, time; sys.path.insert(0, %r); import session; import finish_stop_hook;"
                "finish_stop_hook.TIME_BUDGET_S = 1;"
                "session.finish_status = lambda d: time.sleep(30);"
                "sys.exit(finish_stop_hook.main())" % SCRIPTS)
        r = subprocess.run([sys.executable, "-c", code], input=json.dumps(
            {"session_id": SID, "stop_hook_active": False}), capture_output=True, text=True,
            env=job_env(self.tmp), cwd=self.tmp, timeout=20)
        self.assertAllowed(r)

    def test_reason_tells_how_to_hold_for_a_board_pause(self):
        self.results()
        self.record(self.sdir)
        why = self.assertBlocked(self.run_hook())
        self.assertIn("board", why)
        self.assertIn("session.py hold", why)

    def test_hive_record_allows(self):
        self.results()
        self.record(self.sdir, hive="/quobyte/x/session")   # a HIVE entry is checked on HIVE
        self.assertAllowed(self.run_hook())

    def test_malformed_stdin_allows(self):
        self.results()
        self.record(self.sdir)
        for bad in ("{not json", "", "[1, 2]", "null", '"text"'):
            self.assertAllowed(self.run_hook(stdin=bad))

    def test_internal_error_allows(self):
        self.results()
        self.record(self.sdir)
        with open(os.path.join(self.cfg, "transcript_sessions.json"), "w") as fh:
            fh.write('{"%s": {"sessions": 7, "segments": "x"}}' % SID)
        self.assertAllowed(self.run_hook())

    def _reason_for(self, items):
        """reminder() in-process, for a session whose finish items are `items`."""
        from unittest import mock
        sys.path.insert(0, SCRIPTS)
        import finish_stop_hook as fsh
        import save_transcript
        import session
        with mock.patch.object(save_transcript, "lookup",
                               return_value=[{"session_dir": self.sdir}]), \
                mock.patch.object(session, "finish_hold", return_value=None), \
                mock.patch.object(session, "has_results", return_value=True), \
                mock.patch.object(session, "finish_status", return_value=items), \
                mock.patch.object(fsh, "_read_state", return_value={}), \
                mock.patch.object(fsh, "_write_state"):
            return fsh.reminder({"session_id": SID})

    def test_the_delivery_is_asked_about_never_ordered(self):
        """C2 review ux M5: the delivery (core_submission.py deliver --apply: Bioshare links, the
        client email plan) is outward-facing -- the reminder says to ask, never "Do them now"."""
        deliver = {"item": "delivery", "label": "The delivery of the submission", "done": False,
                   "state": "missing", "detail": "not delivered yet",
                   "command": "python3 scripts/core_submission.py deliver ... --apply (step 12b)"}
        methods = {"item": "methods", "label": "Methods", "done": False, "state": "missing",
                   "detail": "not written", "command": "make_methods.py ..."}
        why = self._reason_for([deliver])
        self.assertIn("ASK the person before doing this", why)
        self.assertIn("Do NOT run it unasked", why)
        self.assertNotIn("Do them now", why)
        self.assertNotIn("Do the others now", why)
        self.assertNotIn("-> python3 scripts/core_submission.py deliver", why)
        both = self._reason_for([methods, deliver])
        self.assertIn("Do the others now", both)
        self.assertIn("-> make_methods.py", both)
        self.assertIn("ASK the person", both)
        self.assertIn("Do them now", self._reason_for([methods]))

    def test_hooks_json_registers_stop(self):
        with open(os.path.join(SKILL, "hooks", "hooks.json")) as fh:
            hooks = json.load(fh)["hooks"]
        cmds = [h["command"] for g in hooks["Stop"] for h in g["hooks"]]
        self.assertTrue(cmds and all("${CLAUDE_PLUGIN_ROOT}/scripts/finish_stop_hook.py" in c
                                     for c in cmds))
        self.assertIn("PreCompact", hooks)                 # the transcript hooks stay


if __name__ == "__main__":
    unittest.main()
