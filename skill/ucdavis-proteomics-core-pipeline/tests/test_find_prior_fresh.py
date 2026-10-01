#!/usr/bin/env python3
"""
A new session is a FRESH analysis by default (SKILL.md step 1b).

`session.py find-prior` used to answer a match with "re-analysis of <session>", and step 1b
said a match "is a re-analysis" -- so a Core staff member starting a fresh run was steered into
nesting under, and reading from, an old session. A match is now a fact to mention in one line;
nothing is read from the earlier session unless the user asks to re-analyse it.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SESSION = os.path.join(SKILL, "scripts", "session.py")
sys.path.insert(0, HERE)
from job_env import job_env     # noqa: E402


class FindPriorIsFreshByDefault(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.raw = os.path.join(self.d, "raw")
        self.files = [os.path.join(self.raw, f"S{i}.d") for i in (1, 2, 3)]
        for f in self.files:
            os.makedirs(f)
        self.env = job_env(self.d, HOME=self.d)

    def tearDown(self):
        self._tmp.cleanup()

    def session(self, *args):
        r = subprocess.run([sys.executable, SESSION, *args], capture_output=True, text=True,
                           env=self.env, cwd=self.d, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        return json.loads(r.stdout)

    def earlier_session(self):
        """An earlier, finished analysis of the same files, beside them (init's default)."""
        prior = os.path.join(self.raw, "2026-09-01_Earlier")
        os.makedirs(os.path.join(prior, "input"))
        with open(os.path.join(prior, "input", "raw_files.txt"), "w") as fh:
            fh.write("# Raw MS files used in this analysis\n" + "\n".join(self.files) + "\n")
        with open(os.path.join(prior, "input", "conditions.csv"), "w") as fh:
            fh.write("File.Name,Group\nS1,A\nS2,B\nS3,B\n")
        return prior

    def test_a_match_is_mentioned_not_reused(self):
        prior = self.earlier_session()
        j = self.session("find-prior", "--raw", os.path.join(self.raw, "*.d"))
        self.assertEqual([m["session"] for m in j["matches"]], [prior])
        self.assertTrue(j["matches"][0]["same_dataset"])
        self.assertEqual(j["default"], "fresh")
        self.assertFalse(j["suggestion"].startswith("re-analysis"), j["suggestion"])
        self.assertIn("fresh analysis", j["suggestion"])
        self.assertIn(prior, j["suggestion"])
        self.assertIn("only if the user asks", j["suggestion"])

    def test_no_match_is_fresh_too(self):
        j = self.session("find-prior", "--raw", os.path.join(self.raw, "*.d"))
        self.assertEqual((j["matches"], j["default"]), ([], "fresh"))
        self.assertIn("fresh analysis", j["suggestion"])

    def test_init_without_reanalysis_of_takes_nothing_from_the_earlier_session(self):
        prior = self.earlier_session()
        j = self.session("init", "--name", "Fresh run", "--raw", *self.files)
        new = j["created"]
        self.assertEqual((j["placement"], j["reanalysis_of"]), ("with-raw-data", None))
        self.assertNotIn(os.path.join(prior, "reanalysis"), new)
        self.assertFalse(os.path.exists(os.path.join(new, ".reanalysis_of")))
        self.assertFalse(os.path.exists(os.path.join(new, "input", "conditions.csv")))

    def test_skill_md_0c_a_new_request_is_not_a_resume(self):
        """Review: "only unfinished work is a resume" made an old session with a COMPLETED search
        and no DE into a resume of a NEW request -- skip the search, continue from its
        next_commands. The condition is the user coming back, and it comes before the bullets."""
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            md = fh.read()
        step = md[md.index("### 0c."):md.index("### 1. ")]
        rule = step.index("resume only when it is unfinished work the user is coming back to")
        self.assertLess(rule, step.index("**COMPLETED"))
        self.assertIn("A new request (new settings, or \"analyse\nthese files\") is fresh", step)
        self.assertIn("never resubmit one that is RUNNING without asking", step)
        self.assertNotIn("Only unfinished work is a resume", step)

    def test_skill_md_step_1b_says_fresh_by_default(self):
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            md = fh.read()
        step = md[md.index("### 1b."):md.index("### 1c.")]
        self.assertIn("fresh analysis by default", step)
        self.assertIn("read nothing from that session", step)
        self.assertIn("only when the user asks", step)
        self.assertNotIn("this is a **re-analysis**", step)
        with open(os.path.join(SKILL, "scripts", "session.py"), encoding="utf-8") as fh:
            self.assertNotIn('("re-analysis of " + hits[0]', fh.read())


if __name__ == "__main__":
    unittest.main(verbosity=2)
