#!/usr/bin/env python3
"""finalize --zip leaves logs/ out, except the redacted logs/commands.log (I11).

gstaff, 2026-10-21 (a Core analysis): the session zip packed most of logs/ -- the expert reviews, the
search job logs and an agent-named qc_bracket output carrying the staff-only verdict line with
no marker -- although SKILL.md says logs/ stays out. Now logs/ is Core-internal as a whole: the
zip holds only logs/commands.log, with staff names replaced by the role (staff.redact). The
existing logs/conversation and logs/decisions.md labels keep their own counts; the new logs/
label counts the rest. Stdlib only.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

LOGS_LABEL = ("logs/ (Core-internal: reviews, job logs, decisions -- kept on disk; "
              "commands.log is in, redacted)")


def finish_fixture(sdir):
    """What a finished session holds before finalize: the skill's own report, the bundle and the
    Methods, so finalize's finish gate and its stamp pass (test_finish_gates.py). Nothing here is
    under logs/."""
    import finish_fixture as ff
    out = os.path.join(sdir, "output")
    ff.finished(out)
    with open(os.path.join(out, "methods.md"), "w", encoding="utf-8") as fh:
        fh.write("# Methods\n")


class ZipLeavesLogsOut(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.tmp = self._td.name
        self.env = dict(os.environ, RECORD_RUN="off",
                        SKILL_CONFIG_DIR=os.path.join(self.tmp, "cfg"),
                        CLAUDE_CONFIG_DIR=os.path.join(self.tmp, "cc"))
        self.env.pop("CLAUDE_CODE_SESSION_ID", None)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init",
                            "--name", "Zip logs", "--base", os.path.join(self.tmp, "base")],
                           capture_output=True, text=True, env=self.env)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.sdir = json.loads(r.stdout)["paths"]["session_dir"]
        finish_fixture(self.sdir)

    def tearDown(self):
        self._td.cleanup()

    def write(self, rel, text):
        path = os.path.join(self.sdir, rel)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)

    def test_zip_holds_only_the_redacted_commands_log_under_logs(self):
        self.write("logs/review_expert1.md", "# Expert review\n\nThe statistician said ...\n")
        self.write("logs/search_25168141.log", "DIA-NN 2.3.0 ... job log\n")
        # an agent-named copy of qc_bracket's output: the verdict line, but no staff-only marker
        self.write("logs/qc_out.json", json.dumps({"verdict": "concern", "note": "Lumos drift"}))
        self.write("logs/commands.log",
                   "python3 normalization_check.py decide --by alice --reason \"low IgG\"\n")
        # the decisions log keeps its own label and count
        self.write("logs/decisions.md", "# Decisions log\n")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "finalize",
                            "--dir", self.sdir, "--zip", "--no-deposit", "--no-notify"],
                           capture_output=True, text=True, env=self.env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        res = json.loads(r.stdout)
        base = os.path.basename(self.sdir)
        names = zipfile.ZipFile(res["zip"]).namelist()
        under_logs = [n for n in names if n.startswith(base + "/logs/")]
        self.assertEqual(under_logs, [base + "/logs/commands.log"], names)
        text = zipfile.ZipFile(res["zip"]).read(base + "/logs/commands.log").decode()
        self.assertIn("<Core staff>", text)
        self.assertNotIn("alice", text)
        excluded = res["zip_excluded"]
        # the three above, and logs/skill_docs.json (what the skill last wrote to README.md,
        # README.html and AGENTS.md: report_files.SKILL_DOCS_RECORD, L5b review M1)
        self.assertEqual(excluded[LOGS_LABEL], 4, excluded)
        self.assertTrue(os.path.isfile(os.path.join(self.sdir, "logs", "skill_docs.json")))
        self.assertEqual(excluded["logs/decisions.md (the decisions log: Core-internal, kept on "
                                  "disk)"], 1)
        # kept on disk where they are
        for rel in ("logs/review_expert1.md", "logs/search_25168141.log", "logs/qc_out.json"):
            self.assertTrue(os.path.isfile(os.path.join(self.sdir, rel)), rel)
        with open(os.path.join(self.sdir, "logs", "commands.log")) as fh:
            self.assertIn("--by alice", fh.read(), "the session's own log keeps the command")

    def test_skill_md_says_logs_stays_out_except_commands_log(self):
        with open(os.path.join(os.path.dirname(HERE), "SKILL.md"), encoding="utf-8") as fh:
            skill = " ".join(fh.read().split())
        self.assertIn("`logs/` is Core-internal: it is never delivered, and the session zip leaves "
                      "it out except `logs/commands.log`, whose staff names are replaced by the "
                      "role.", skill)


if __name__ == "__main__":
    unittest.main()
