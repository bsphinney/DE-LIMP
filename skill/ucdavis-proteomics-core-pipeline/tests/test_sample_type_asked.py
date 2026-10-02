#!/usr/bin/env python3
"""
The sample type is ASKED, never inferred (a staff report, 2.10).

An abbreviation in a cohort's file names was written up as a tissue through the draft report and
Methods, and the data did not support it. SKILL.md asked what the samples are only to decide
keratin.
Now step 3 asks it for every analysis, sample_type.py records the user's words in the session,
and the analysis brief (analysis_prompt.py) quotes them -- or says NOT STATED and forbids naming
a tissue, cell type or biofluid, from file names or otherwise.

Synthetic sessions on a temp dir; no network.
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

import sample_type as st  # noqa: E402


def run(script, *args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args],
                          capture_output=True, text=True, timeout=120)


def session(tmp):
    """A session as session.py init leaves one: input/ and session.json (no submission)."""
    s = os.path.join(tmp, "2026-10-01_cohort")
    os.makedirs(os.path.join(s, "input"))
    with open(os.path.join(s, "session.json"), "w") as fh:
        json.dump({"name": "cohort"}, fh)
    return s


class TheRecord(unittest.TestCase):
    def test_the_users_words_are_recorded_as_given(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            p = run("sample_type.py", "set", "--session", s, "--stated",
                    "  mouse   hind-limb muscle ")
            self.assertEqual(p.returncode, 0, p.stderr)
            rec = st.load(s)
            self.assertEqual(rec["stated"], "mouse hind-limb muscle")
            self.assertEqual(rec["source"], "user")
            self.assertIn("never infer it from file or folder names", rec["rule"])
            show = json.loads(run("sample_type.py", "show", "--session", s).stdout)
            self.assertEqual(show["sample_type"], "mouse hind-limb muscle (as the user stated it)")

    def test_a_user_who_does_not_know_is_recorded_as_such(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            self.assertEqual(run("sample_type.py", "set", "--session", s,
                                 "--not-stated").returncode, 0)
            self.assertIsNone(st.load(s)["stated"])
            self.assertIn("NOT STATED", st.describe(st.load(s)))

    def test_nothing_is_written_without_an_answer(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            p = run("sample_type.py", "set", "--session", s, "--stated", "   ")
            self.assertNotEqual(p.returncode, 0)
            self.assertFalse(os.path.exists(st.path_for(s)))
            p = run("sample_type.py", "set", "--session", tmp, "--stated", "plasma")
            self.assertNotEqual(p.returncode, 0, "not a session directory")
            self.assertIn("not a session directory", p.stderr)
        self.assertEqual(st.describe(None),
                         "NOT STATED by the user (not recorded in this session -- ask the user)")


class TheBrief(unittest.TestCase):
    def brief(self, tmp, *extra):
        de = os.path.join(tmp, "tables")
        os.makedirs(de, exist_ok=True)
        out = os.path.join(tmp, "PROMPT.md")
        p = run("analysis_prompt.py", "--de-dir", de, "--out", out, *extra)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out, encoding="utf-8") as fh:
            return fh.read(), p.stderr

    def test_a_stated_sample_type_is_quoted(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            st.record(s, "mouse hind-limb muscle")
            brief, _ = self.brief(tmp, "--session", s)
        self.assertIn("## What the samples are — as the user stated it", brief)
        self.assertIn("**Sample type: mouse hind-limb muscle (as the user stated it)**", brief)
        self.assertIn("adding nothing they do not state", brief)
        self.assertIn("as a question for the user", brief)

    def test_without_a_record_no_tissue_may_be_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            brief, _ = self.brief(tmp, "--session", s)
        self.assertIn("**Sample type: NOT STATED by the user (not recorded in this session -- "
                      "ask the user)**", brief)
        self.assertIn("Do not name a tissue, cell type or biofluid anywhere", brief)
        self.assertIn("never infer it from file or folder names", brief)
        with tempfile.TemporaryDirectory() as tmp:
            brief, _ = self.brief(tmp)                    # no session at all: the same rule
        self.assertIn("Sample type: NOT STATED", brief)

    def test_the_submission_session_is_read_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            st.record(s, "HeLa cells")
            brief, _ = self.brief(tmp, "--submission", s)
        self.assertIn("**Sample type: HeLa cells (as the user stated it)**", brief)

    def test_an_unreadable_record_is_said_never_skipped(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            with open(st.path_for(s), "w") as fh:
                fh.write("{not json")
            brief, err = self.brief(tmp, "--session", s)
        self.assertIn("Sample type: NOT STATED", brief)
        self.assertIn("its record could not be read", brief)
        self.assertIn("could not be read", err)


class SkillMdAsksIt(unittest.TestCase):
    def test_step_3_asks_it_for_every_analysis_and_step_9_passes_the_session(self):
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            md = fh.read()
        self.assertIn("**Sample type — ask, for every analysis, and never infer it.**", md)
        self.assertIn("scripts/sample_type.py set --session <session> --stated", md)
        self.assertIn("--workflow-manifest ./wf/workflow.manifest.json --session <session>", md)
        self.assertIn("**sample type** (tissue,\n   cell line, biofluid, IP …) is likewise "
                      "**asked, never inferred**", md)


if __name__ == "__main__":
    unittest.main()
