#!/usr/bin/env python3
"""What a Core staff analysis is for (skill 2.12, SKILL.md step 0a).

Brett (2026-10-15): "automatically create an entry on the project board for every analysis my
staff does and the skill should ask for a submission number when they start the analysis."
Staff answer with one or more CoreOmics submission numbers, or "none" for internal work. The
answer goes on the Core Project Board before the long steps, and `session.py init
--submission ... | --internal "<title>"` keeps it in session.json's "analysis_for", so later
steps and a resume (`checkpoint.py status`) know it without asking again. An analysis by
anyone else records nothing, and the record never makes a run the Core's (its feedback
survey stays with the attached CoreOmics record).
stdlib only.
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
from job_env import job_env     # noqa: E402
import session                  # noqa: E402
import submission_report as sr  # noqa: E402


class AnalysisFor(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.raw = os.path.join(self.d, "raw")
        self.files = [os.path.join(self.raw, f"S{i}.d") for i in (1, 2)]
        for f in self.files:
            os.makedirs(f)
        self.env = job_env(self.d, HOME=self.d)

    def tearDown(self):
        self._tmp.cleanup()

    def run_py(self, script, *args):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args],
                              capture_output=True, text=True, env=self.env, cwd=self.d, timeout=120)

    def init(self, *extra, name="PROT_0807"):
        r = self.run_py("session.py", "init", "--name", name, "--date", "2026-10-15",
                        "--raw", *self.files, *extra)
        self.assertEqual(r.returncode, 0, r.stderr)
        return json.loads(r.stdout)

    def session_json(self, s):
        with open(os.path.join(s, "session.json"), encoding="utf-8") as fh:
            return json.load(fh)

    def test_submission_numbers_are_kept_in_order_the_first_the_main_one(self):
        j = self.init("--submission", "807", "--submission", "PROT_0822", "--submission", "prot-807")
        want = {"submissions": ["PROT_0807", "PROT_0822"]}
        self.assertEqual(j["analysis_for"], want)
        self.assertEqual(self.session_json(j["created"])["analysis_for"], want)
        self.assertEqual(session.analysis_for(j["created"]), want)

    def test_internal_work_keeps_its_board_title(self):
        j = self.init("--internal", "  Astral  HeLa gradient test ", name="Astral HeLa gradient test")
        self.assertEqual(self.session_json(j["created"])["analysis_for"],
                         {"internal": "Astral HeLa gradient test"})

    def test_without_either_nothing_is_recorded(self):
        """An outside user's analysis (never asked) is exactly as before: no session.json."""
        j = self.init()
        self.assertIsNone(j["analysis_for"])
        self.assertFalse(os.path.exists(os.path.join(j["created"], "session.json")))
        self.assertIsNone(session.analysis_for(j["created"]))

    def test_a_wrong_answer_is_refused_before_anything_is_made(self):
        for extra, said in ((["--submission", "807", "--internal", "Method test"], "not both"),
                            (["--submission", "99922f5337f8"], "PROT number"),
                            (["--submission", "abc"], "not a submission number"),
                            (["--internal", "   "], "title")):
            with self.subTest(extra=extra):
                r = self.run_py("session.py", "init", "--name", "Bad", "--raw", *self.files, *extra)
                self.assertNotEqual(r.returncode, 0)
                self.assertIn(said, r.stderr)
        self.assertEqual(sorted(os.listdir(self.raw)), ["S1.d", "S2.d"], "no session was made")

    def test_attach_and_init_keep_each_others_records(self):
        j = self.init("--submission", "807")
        s = j["created"]
        sr.attach(s, sr.sanitize({"internal_id": "PROT_0807"}))
        sj = self.session_json(s)
        self.assertEqual(sj["analysis_for"], {"submissions": ["PROT_0807"]})
        self.assertEqual(sj["coreomics"]["internal_id"], "PROT_0807")
        self.init("--submission", "807", "--submission", "822")       # init again: same folder
        sj = self.session_json(s)
        self.assertEqual(sj["coreomics"]["internal_id"], "PROT_0807")
        self.assertEqual(sj["analysis_for"]["submissions"], ["PROT_0807", "PROT_0822"])

    def test_a_damaged_session_json_is_not_overwritten(self):
        s = self.init()["created"]
        with open(os.path.join(s, "session.json"), "w") as fh:
            fh.write("{not json")
        r = self.run_py("session.py", "init", "--name", "PROT_0807", "--date", "2026-10-15",
                        "--raw", *self.files, "--submission", "807")
        self.assertNotEqual(r.returncode, 0)
        with open(os.path.join(s, "session.json")) as fh:
            self.assertEqual(fh.read(), "{not json")

    def test_the_record_never_makes_a_run_the_cores(self):
        """The feedback survey keys on the attached CoreOmics record alone (core_run), so internal
        work, or a number staff gave with no record attached, never shows it."""
        for extra in (["--internal", "Method test"], ["--submission", "807"]):
            with self.subTest(extra=extra):
                s = self.init(*extra, name="Run " + extra[0])["created"]
                self.assertIsNone(sr.core_run(s))
                self.assertIsNone(sr.load(s))

    def test_a_resume_reads_it_from_checkpoint_status(self):
        for extra, want in ((["--internal", "Method test"], {"internal": "Method test"}), ([], None)):
            with self.subTest(extra=extra):
                s = self.init(*extra, name="Resume " + str(bool(extra)))["created"]
                with open(os.path.join(s, ".recovery.json"), "w") as fh:
                    json.dump({"stages": [{"stage": "search", "jobs": [], "state": "submitted"}]}, fh)
                r = self.run_py("checkpoint.py", "status", "--session", s)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual(json.loads(r.stdout).get("analysis_for"), want)
                self.assertFalse([w for w in json.loads(r.stdout).get("warnings", [])
                                  if "NOT KNOWN" in w])

    def test_a_damaged_session_json_is_said_on_resume_not_dropped(self):
        """2.12 review: a damaged session.json dropped analysis_for from `status` without a word,
        so a resume lost which submission and board thread it was."""
        s = self.init("--submission", "807")["created"]
        with open(os.path.join(s, ".recovery.json"), "w") as fh:
            json.dump({"stages": [{"stage": "search", "jobs": [], "state": "submitted"}]}, fh)
        sj = os.path.join(s, "session.json")
        for text, said in (("{not json", "could not be read (JSONDecodeError"),
                           ('["a list"]', "is not a JSON object"),
                           ('{"analysis_for": "PROT_0807"}', "that is not an object")):
            with self.subTest(text=text):
                with open(sj, "w", encoding="utf-8") as fh:
                    fh.write(text)
                self.assertIsNone(session.analysis_for(s))
                self.assertIn(said, session.read_analysis_for(s)[1])
                r = self.run_py("checkpoint.py", "status", "--session", s)
                self.assertEqual(r.returncode, 0, r.stderr)
                out = json.loads(r.stdout)
                self.assertNotIn("analysis_for", out)
                (w,) = [w for w in out.get("warnings", []) if "NOT KNOWN" in w]
                self.assertIn(sj, w)
                self.assertIn(said, w)
                self.assertIn("--submission <PROT number>", w)
                with open(sj, encoding="utf-8") as fh:
                    self.assertEqual(fh.read(), text, "status rewrote session.json")

    def test_no_session_json_is_not_a_warning(self):
        s = self.init()["created"]
        self.assertEqual(session.read_analysis_for(s), (None, None))


def _read(*parts):
    with open(os.path.join(SKILL, *parts), encoding="utf-8") as fh:
        return fh.read()


def _flat(text):
    """One line: no wrapping, no bold, no shell line continuations."""
    return " ".join(text.replace("\\\n", " ").replace("**", "").split())


class TheDocsSayIt(unittest.TestCase):
    """SKILL.md asks the question of Core staff only and passes the answer on; the references
    say what each answer does on the board and in the session."""

    @classmethod
    def setUpClass(cls):
        cls.skill = _read("SKILL.md")
        cls.board = _read("references", "board.md")
        cls.core = _read("references", "core-submissions.md")

    def test_step_0a_asks_core_staff_only(self):
        step = _flat(self.skill[self.skill.index("### 0a."):self.skill.index("### 0. ")])
        self.assertIn("Core staff only, as soon as they say yes to 2: \"Which CoreOmics submission "
                      "is this analysis for?\"", step)
        self.assertIn("\"none\" for internal work", step)
        self.assertIn("on a resume (step 0c)", step)
        self.assertIn("Anyone else: never ask it, and never mention the board", step)
        self.assertIn("--submission PROT_0807", step)
        self.assertIn("--internal \"<its title on the board>\"", step)

    def test_a_resume_finds_its_thread_without_asking(self):
        step = _flat(self.skill[self.skill.index("### 0c."):self.skill.index("### 1. ")])
        self.assertIn("a resume asks no submission question", step)
        self.assertIn("analysis_for", step)
        self.assertIn("bash scripts/board.sh threads", step)

    def test_internal_work_goes_on_the_board_in_step_1(self):
        step = _flat(self.skill[self.skill.index("### 1. "):self.skill.index("### 1a.")])
        self.assertIn("None (internal work) → no `identify` and no submission. Put it on the board "
                      "now, before step 2", step)
        self.assertIn("bash scripts/board.sh start --internal --title", step)
        self.assertIn("--hours 336 --with bsphinney@ucdavis.edu", step)

    def test_the_core_session_records_its_submission(self):
        step = self.skill[self.skill.index("### 1c."):self.skill.index("### 2. ")]
        self.assertIn("--raw $(cat ~/core/PROT_0807/files.txt) --submission PROT_0807'", step)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init", "--help"],
                           capture_output=True, text=True, timeout=60)
        self.assertIn("--submission", r.stdout)
        self.assertIn("--internal", r.stdout)

    def test_board_md_starts_every_kind_of_analysis(self):
        ref = _flat(self.board)
        self.assertIn("### One submission (step 1c, right after `fetch`)", self.board)
        self.assertIn("### Several submissions in one analysis", self.board)
        self.assertIn("### Internal work (staff said \"none\")", self.board)
        self.assertIn("bash scripts/board.sh start --internal --title \"<the analysis's short name>\" "
                      "--goal", ref)
        self.assertIn("What the first `start` writes is kept.", ref)
        self.assertIn("Every record in section 3 goes to the main thread only", ref)
        self.assertIn("Internal work has no stage and no delivery, and no QC row", ref)
        self.assertIn("--submission PROT_0807", ref)
        self.assertIn("--internal \"<the title it has on the board>\"", ref)

    def test_several_submissions_attach_the_main_record_only(self):
        sec = _flat(self.core[self.core.index("## Several submissions in one analysis"):
                              self.core.index("## Troubleshooting by exit code")])
        self.assertIn("Attach the main submission's record only.", sec)
        self.assertIn("--submission PROT_0807 --submission PROT_0811", sec)
        self.assertIn("deliver takes the combined session into the main submission's share only",
                      sec.replace("`", ""))
        self.assertIn("never work around that", sec)


if __name__ == "__main__":
    unittest.main(verbosity=2)
