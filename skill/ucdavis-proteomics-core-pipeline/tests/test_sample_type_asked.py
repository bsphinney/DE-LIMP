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

    def test_core_staff_not_asked_is_recorded_as_not_asked(self):
        """A staff report (2026-10-08): "stop asking what tissue! irrelevant!". Core staff are asked
        only when a step depends on it, and their record says it was NOT ASKED -- never that they
        were asked and did not say."""
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            p = run("sample_type.py", "set", "--session", s, "--not-asked")
            self.assertEqual(p.returncode, 0, p.stderr)
            rec = st.load(s)
            self.assertIsNone(rec["stated"])
            self.assertEqual(rec["source"], st.NOT_ASKED)
            self.assertEqual(st.describe(rec), f"NOT STATED by the user ({st.NOT_ASKED})")
            self.assertNotIn("asked; the user did not say", st.describe(rec))
            p = run("sample_type.py", "set", "--session", s, "--not-asked", "--not-stated")
            self.assertNotEqual(p.returncode, 0, "one of the three, never two")
            # a later answer (a flagged panel, 8c) replaces it
            self.assertEqual(run("sample_type.py", "set", "--session", s, "--stated",
                                 "mouse brain IP").returncode, 0)
            self.assertEqual(st.load(s)["source"], "user")
        with self.assertRaises(ValueError):
            st.record(tmp, None, user_did_not_say=True, not_asked=True)

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

    def test_core_staff_not_asked_adds_no_data_quality_note(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp)
            st.record(s, not_asked=True)
            brief, _ = self.brief(tmp, "--session", s)
        self.assertIn("Core staff were not asked", brief)
        self.assertIn("as the attached submission record does, when there is one, or else by "
                      "their groups only", brief)
        self.assertIn("never infer it from file or folder names", brief)
        self.assertNotIn("Say in **Data Quality Notes** that the sample type was not stated", brief)

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
    def test_step_3_asks_outside_users_and_core_staff_only_when_it_matters(self):
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            md = fh.read()
        self.assertIn("**Sample type — never infer it; ask outside users, and Core staff only when "
                      "a choice\n  depends on it.**", md)
        self.assertIn("**Outside users: ask**", md)
        # any Core staff analysis, internal work too (staff review of 2.12): not asked, or what
        # they already said recorded as stated
        self.assertIn("**Core staff, any analysis (a CoreOmics submission or internal work): do "
                      "not ask it.**\n    Record `--not-asked` (step 3b) — or `--stated` with what "
                      "they already said (\"HeLa digest\n    standard\")", md)
        self.assertIn("a submission record, when there is one, describes the samples", md)
        self.assertIn("though Core staff (any analysis: a\n   submission or internal work) are "
                      "asked it only when a step depends on it", md)
        self.assertIn("scripts/sample_type.py set --session <session> --stated", md)
        self.assertIn("scripts/sample_type.py set --session <session> --not-asked", md)
        self.assertIn("--workflow-manifest ./wf/workflow.manifest.json --session <session>", md)
        self.assertIn("**sample type** (tissue,\n   cell line, biofluid, IP …) is likewise "
                      "**asked, never inferred**", md)
        self.assertIn("ask what the samples are only when a panel flags", md)


if __name__ == "__main__":
    unittest.main()
