#!/usr/bin/env python3
"""
2.15.1 review, finding 4 (and minor 11): the Core's instrument acknowledgment was printed only
with an attached CoreOmics record, a core-client mark or --acquired-at-core, so Core runs lost it:
a Core staff session started with --internal, and one started with --submission PROT_... whose
record was not attached yet. Now submission_report.core_acquisition() is THE test, and it reads
session.json's analysis_for (init --submission / --internal) too. A standalone Methods (no
session) stays a tagged blank -- SKILL.md 9d tells the agent to ask the user and pass
--acquired-at-core -- and a session with no Core mark, a missing or a damaged session.json never
gets the sentence. The Core's LC defaults (column, emitter, column temperature, mobile phases)
follow the same record (minor 11).

Synthetic sessions and a synthetic timsTOF .d only. stdlib.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
SKILL = os.path.dirname(HERE)
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_methods as mm                                        # noqa: E402
import submission_report as sr                                   # noqa: E402
from test_instrument_methods import make_timstof_d, section     # noqa: E402

HHMI = "Howard Hughes Medical Institute"


def session(tmp, sj):
    os.makedirs(os.path.join(tmp, "input"), exist_ok=True)
    os.makedirs(os.path.join(tmp, "output"), exist_ok=True)
    if sj is not None:
        with open(os.path.join(tmp, "session.json"), "w") as fh:
            fh.write(sj if isinstance(sj, str) else json.dumps(sj))
    return tmp


class TheOneTest(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.tmp = self._t.name

    def tearDown(self):
        self._t.cleanup()

    def why(self, sj):
        return sr.core_acquisition(session(self.tmp, sj))

    def test_internal_core_work_is_core(self):
        self.assertIn("internal work", self.why({"analysis_for": {"internal": "method test"}}))

    def test_an_announced_submission_without_its_record_is_core(self):
        self.assertIn("PROT_9001", self.why({"analysis_for": {"submissions": ["PROT_9001"]}}))

    def test_an_attached_coreomics_block_is_core(self):
        self.assertIn("PROT_9001", self.why({"coreomics": {"internal_id": "PROT_9001"}}))

    def test_nothing_on_record_is_not_core(self):
        self.assertIsNone(self.why({"schema": 1}))
        self.assertIsNone(self.why({"analysis_for": {}}))
        self.assertIsNone(self.why({"analysis_for": {"internal": "   "}}))

    def test_missing_or_damaged_session_json_is_not_core(self):
        self.assertIsNone(self.why(None))
        self.assertIsNone(self.why("{not json"))
        self.assertIsNone(self.why(["a list"]))
        self.assertIsNone(sr.core_acquisition(None))

    def test_make_methods_reads_the_same_test(self):
        session(self.tmp, {"analysis_for": {"internal": "method test"}})
        a = SimpleNamespace(acquired_at_core=False, session=self.tmp, submission=None)
        self.assertTrue(mm.core_acquisition_recorded(a))
        a = SimpleNamespace(acquired_at_core=True, session=None, submission=None)
        self.assertIn("--acquired-at-core", mm.core_acquisition_recorded(a))
        a = SimpleNamespace(acquired_at_core=False, session=None, submission=None)
        self.assertIsNone(mm.core_acquisition_recorded(a))


class EndToEnd(unittest.TestCase):
    """make_methods.py on a synthetic timsTOF .d, with each kind of session."""

    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.tmp = self._t.name
        self.d = make_timstof_d(os.path.join(self.tmp, "data", "r1.d"))

    def tearDown(self):
        self._t.cleanup()

    def run_mm(self, sj="__none__", extra=()):
        out = os.path.join(self.tmp, "m.md")
        cmd = [sys.executable, os.path.join(SCRIPTS, "make_methods.py"), "--raw", self.d,
               "--out", out, *extra]
        if sj != "__none__":
            sd = session(os.path.join(self.tmp, "s"), sj)
            cmd += ["--session", sd]
        r = subprocess.run(cmd, capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out, encoding="utf-8") as fh:
            md = fh.read()
        with open(os.path.splitext(out)[0] + "_params.json", encoding="utf-8") as fh:
            return md, json.load(fh)

    def assertCore(self, md, params):
        self.assertIn(HHMI, section(md, "Acknowledgments"))
        self.assertIn(f"{mm.LC_COLUMN_DEFAULT} {mm.DEF}", section(md, "Liquid chromatography"))
        self.assertIn(f"{mm.EMITTER_DEFAULT} {mm.DEF}", section(md, "Mass spectrometry"))
        self.assertTrue(params["core_acquisition"])

    def assertNotCore(self, md, params):
        ack = section(md, "Acknowledgments")
        self.assertNotIn(HHMI, ack)
        self.assertIn("add your own facility's", ack)
        self.assertNotIn(mm.LC_COLUMN_DEFAULT, md)
        self.assertNotIn(mm.EMITTER_DEFAULT, md)
        self.assertNotIn(mm.COLUMN_TEMP_DEFAULT, md)
        self.assertNotIn(mm.DEF, section(md, "Liquid chromatography"))
        self.assertIsNone(params["core_acquisition"])

    def test_internal_session_gets_the_acknowledgment(self):
        self.assertCore(*self.run_mm({"analysis_for": {"internal": "method test"}}))

    def test_announced_submission_gets_the_acknowledgment(self):
        self.assertCore(*self.run_mm({"analysis_for": {"submissions": ["PROT_9001"]}}))

    def test_standalone_with_the_users_yes_gets_it(self):
        self.assertCore(*self.run_mm(extra=["--acquired-at-core"]))

    def test_standalone_without_an_answer_is_a_tagged_blank(self):
        self.assertNotCore(*self.run_mm())

    def test_a_session_with_no_core_mark_is_a_tagged_blank(self):
        self.assertNotCore(*self.run_mm({"schema": 1}))

    def test_a_missing_session_json_is_a_tagged_blank(self):
        self.assertNotCore(*self.run_mm(None))


class SkillTellsTheAgentToAsk(unittest.TestCase):
    def test_9d_asks_whether_the_core_acquired_the_data(self):
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            s = fh.read()
        step = " ".join(s[s.index("### 9d."):s.index("### 10.")].split())
        self.assertIn("ask the user", step)
        self.assertIn("UC Davis Proteomics Core?", step)
        self.assertIn("--acquired-at-core", step)
        fm = s.split("\n---", 1)[0]
        self.assertIn("write the Methods for its analysis or from raw files", " ".join(fm.split()))


if __name__ == "__main__":
    unittest.main()
