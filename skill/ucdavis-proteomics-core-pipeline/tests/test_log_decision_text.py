#!/usr/bin/env python3
"""log_decision.py takes its text as DATA (2.12 security review).

The decisions log's text is the user's own words and the evidence: quotes, apostrophes, an error
text with backticks or $(...). In `--what "..."` / `--why "..."` bash runs those -- on this
computer, and again on HIVE through hive_exec.sh's `bash -l -c` (the review reproduced a backtick
running a command through report_issue.sh's --what). Now --what-file / --why-file F|- read it
from a file or standard input (both -: the what, a line `---`, the why), and SKILL.md's examples
send it in a heredoc quoted <<'EOF'. These run the documented commands themselves.
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
LOG_DECISION = os.path.join(SCRIPTS, "log_decision.py")


def _read(*parts):
    with open(os.path.join(SKILL, *parts), encoding="utf-8") as fh:
        return fh.read()


class TextAsData(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._tmp.name)
        self.session = os.path.join(self.d, "session")
        os.makedirs(self.session)
        self.ran = [os.path.join(self.d, n) for n in ("RAN_BACKTICK", "RAN_DOLLAR")]

    def tearDown(self):
        self._tmp.cleanup()

    def run_py(self, *args, stdin=None):
        env = {k: v for k, v in os.environ.items() if k != "CLAUDE_CODE_SESSION_ID"}
        return subprocess.run([sys.executable, LOG_DECISION, "--session", self.session, *args],
                              input=stdin, capture_output=True, text=True, timeout=60, env=env)

    def log(self):
        with open(os.path.join(self.session, "logs", "decisions.md"), encoding="utf-8") as fh:
            return fh.read()

    def hostile(self):
        return (f"Contrasts: `touch {self.ran[0]}` Old - Young",
                f"the user's words: \"yes, $(touch {self.ran[1]}) that's right\" \\ $HOME")

    def assert_never_run(self):
        for f in self.ran:
            self.assertFalse(os.path.exists(f), "the text was run as a command")

    def test_both_from_standard_input(self):
        what, why = self.hostile()
        r = self.run_py("--step", "3. design", "--what-file", "-", "--why-file", "-",
                        stdin=f"{what}\n---\n{why}\nsecond line\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIs(json.loads(r.stdout)["logged"], True)
        self.assert_never_run()
        log = self.log()
        self.assertIn(f"-- {what}\n- **Why:** {why} second line\n- **Step:** 3. design", log)

    def test_each_from_a_file_or_one_from_stdin(self):
        what, why = self.hostile()
        fw, fy = os.path.join(self.d, "what.txt"), os.path.join(self.d, "why.txt")
        with open(fw, "w", encoding="utf-8") as fh:
            fh.write("﻿" + what + "\r\n")                  # a Notepad BOM, a CRLF
        with open(fy, "w", encoding="utf-8") as fh:
            fh.write(why)
        self.assertEqual(self.run_py("--what-file", fw, "--why-file", fy).returncode, 0)
        r = self.run_py("--what", "Groups: ctrl vs trt", "--why-file", "-", stdin=why)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assert_never_run()
        log = self.log()
        self.assertIn(f"-- {what}\n- **Why:** {why}\n", log)
        self.assertIn(f"-- Groups: ctrl vs trt\n- **Why:** {why}\n", log)
        self.assertNotIn("﻿", log)

    def test_the_old_arguments_still_work(self):
        r = self.run_py("--what", "Groups: Old vs Young", "--why", "the user confirmed")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("-- Groups: Old vs Young\n- **Why:** the user confirmed", self.log())

    def test_wrong_uses_say_why_and_log_nothing(self):
        f = os.path.join(self.d, "w.txt")
        with open(f, "w") as fh:
            fh.write("w")
        for args, stdin, said in (
                (("--what-file", "-", "--why-file", "-"), "no separator here\n", "a line '---'"),
                (("--what", "w", "--what-file", f, "--why", "y"), None, "not allowed with"),
                (("--what-file", os.path.join(self.d, "nope"), "--why", "y"), None, "could not read"),
                (("--what-file", "-", "--why", "y"), "  \n", "the what is empty"),
                (("--what-file", "-", "--why-file", "-"), "w\n---\n\n", "the why is empty"),
                (("--why", "y"), None, "one of the arguments --what --what-file is required")):
            with self.subTest(args=args):
                r = self.run_py(*args, stdin=stdin)
                self.assertEqual(r.returncode, 2, r.stdout)
                self.assertIn(said, r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.session, "logs", "decisions.md")))

    def test_a_secret_in_a_file_is_redacted_as_before(self):
        tok = "ghp_" + "a1B2c3D4e5" * 4
        r = self.run_py("--what", "Overrides", "--why-file", "-", stdin=f"the user pasted {tok}\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn(tok, self.log())

    # ---- the documented commands, run as written --------------------------------------------
    def block(self, doc, marker):
        """The fenced code block of `doc` that holds `marker`."""
        i = doc.index(marker)
        return doc[doc.rindex("```", 0, i) + 3:doc.index("```", i)].strip("\n")

    def test_skill_md_commands_run_as_written(self):
        """The decisions-log block: the local command, and the one through hive_exec.sh (a
        stand-in that runs the command as HIVE's login shell would, with the skill at
        ~/proteomics-pipeline/scripts and standard input passed through, as ssh does)."""
        home = os.path.join(self.d, "home")
        os.makedirs(os.path.join(home, "proteomics-pipeline"))
        os.symlink(SCRIPTS, os.path.join(home, "proteomics-pipeline", "scripts"))
        stub = os.path.join(self.d, "hive_exec_stub.sh")
        with open(stub, "w") as fh:
            fh.write('#!/bin/bash\nexec bash -c "$1"\n')
        skill = _read("SKILL.md")
        marker = "# no usable local Python: the same on HIVE (the text goes over ssh"
        block = self.block(skill, marker)
        local, hive = block.split(marker)
        self.assertIn("<<'EOF'", local)
        self.assertIn("<<'EOF'", hive)
        for cmd in (local, hive.split("\n", 1)[1]):
            with self.subTest(cmd=cmd.split("\n")[0]):
                cmd = (cmd.replace("<S>", self.session)
                       .replace("python3 scripts/", f"{sys.executable} {SCRIPTS}/")
                       .replace("python3 ~/", f"{sys.executable} ~/")
                       .replace("bash scripts/hive_exec.sh", f"bash {stub}"))
                r = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, timeout=60,
                                   env=dict(os.environ, HOME=home), cwd=self.d)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertIs(json.loads(r.stdout.strip().splitlines()[-1])["logged"], True)
        log = self.log()
        self.assertEqual(log.count("-- <what was decided>\n- **Why:** <the user's words, or the "
                                   "evidence>\n- **Step:** 3. design"), 2)

    def test_every_documented_example_reads_standard_input(self):
        examples = [(where, line) for where, doc in (("SKILL.md", _read("SKILL.md")),
                                                     ("core-submissions.md",
                                                      _read("references", "core-submissions.md")))
                    for line in doc.splitlines() if "log_decision.py --session" in line]
        self.assertGreaterEqual(len(examples), 4, examples)
        for where, line in examples:
            with self.subTest(where=where, line=line):
                self.assertIn("--what-file - --why-file -", line)
                self.assertTrue(line.rstrip().endswith("<<'EOF'") or line.rstrip().endswith("-'"),
                                line)
                self.assertNotRegex(line, r'--(what|why) ')
        r = subprocess.run([sys.executable, LOG_DECISION, "--help"], capture_output=True, text=True)
        self.assertIn("--what-file F|-", r.stdout)
        self.assertIn("--why-file F|-", r.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
