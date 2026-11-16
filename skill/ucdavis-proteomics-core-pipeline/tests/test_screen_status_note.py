#!/usr/bin/env python3
"""`session.py status` shows the protein screen data's state (step 8g, make_screen_data.R's
output/screen_data_status.json) as an OPTIONAL note: never a finish item, so the Stop hook and
finalize never ask for it. Stdlib only; no R."""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import session as session_mod  # noqa: E402
from test_core_client_survey import NOT_STAFF, cli  # noqa: E402
from test_submission_report import make_session  # noqa: E402

PROV = {"pipeline_id": "dpc", "contrasts": ["B-A"]}


def status_file(s, **rec):
    with open(os.path.join(s, "output", "screen_data_status.json"), "w") as fh:
        json.dump(dict({"schema": "screen_data_status/1", "route": "dpc"}, **rec), fh)


class ScreenStatusNote(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.s = make_session(self._td.name, prov=PROV)
        self.p = session_mod.paths_for(self.s)
        # the tests here are about the screen note, not the commands log (golden rule 5): keep one
        os.makedirs(os.path.dirname(self.p["commands_log"]), exist_ok=True)
        with open(self.p["commands_log"], "w") as fh:
            fh.write("# a test's command\n")

    def note(self):
        return session_mod.screen_data_note(self.p)

    def test_not_built_by_default(self):
        self.assertEqual(self.note(), "protein screen data: not built (optional)")

    def test_no_de_no_note(self):
        d = os.path.join(self._td.name, "nodE")
        os.makedirs(d)
        s = make_session(d)                                  # no output/tables/de_provenance.json
        self.assertIsNone(session_mod.screen_data_note(session_mod.paths_for(s)))

    def test_built(self):
        status_file(self.s, state="built", reason=None)
        self.assertEqual(self.note(), "protein screen data: built (output/screen_data; optional)")

    def test_refused_names_route_and_reason(self):
        status_file(self.s, state="refused", route="maxlfq", reason="this analysis's quantification is MaxLFQ")
        self.assertEqual(self.note(), "protein screen data: not built, refused (maxlfq): this analysis's "
                                      "quantification is MaxLFQ (optional)")

    def test_failed_says_the_line(self):
        status_file(self.s, state="failed", reason="building the record of P1 failed:\n it was killed")
        n = self.note()
        self.assertTrue(n.startswith("protein screen data: failed: building the record of P1 failed: it was killed"), n)
        self.assertNotIn("\n", n)

    def test_started_and_not_finished(self):
        status_file(self.s, state="running")
        self.assertIn("started and not finished", self.note())

    def test_a_status_line_is_never_a_finish_item(self):
        before = session_mod.finish_status(self.s)
        for rec in ({"state": "failed", "reason": "x"}, {"state": "refused", "reason": "y"},
                    {"state": "running"}):
            status_file(self.s, **rec)
            self.assertEqual(session_mod.finish_status(self.s), before)

    def test_status_shows_it_and_completeness_ignores_it(self):
        def run(*extra):
            return cli(session_mod, ["status", "--session", self.s, *extra], NOT_STAFF)
        rc0, out0, _ = run()
        self.assertIn("  note     protein screen data: not built (optional)", out0)
        status_file(self.s, state="failed", reason="boom")
        rc1, out1, _ = run()
        self.assertEqual(rc1, rc0)                       # the exit code (complete or not) is unchanged
        self.assertIn("  note     protein screen data: failed: boom (optional; step 8g)", out1)
        rc, out, _ = run("--json")
        j = json.loads(out)
        self.assertEqual(j["notes"][-1], "protein screen data: failed: boom (optional; step 8g)")
        self.assertEqual(rc, rc0)
        self.assertNotIn("screen", json.dumps(j["items"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
