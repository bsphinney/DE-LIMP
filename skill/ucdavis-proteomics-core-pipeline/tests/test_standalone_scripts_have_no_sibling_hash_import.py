"""Scripts that must run without the rest of scripts/ keep their own helpers (cleanup audit E).

record_run.py, notify_slack.py, notes.py and slack_collab.py are piped to a remote python on stdin
(hive_exec.sh `python3 -`); report_files.py is loaded from a partial copy by make_report.py;
core_submission.py must reach its "sync the whole scripts/ directory" refusal from a partial copy.
Importing safe_write (the shared sha256_file / load_json) at module top would turn that into a
traceback, so none of them may.
"""
import os
import re
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
STANDALONE = ("record_run.py", "notify_slack.py", "notes.py", "slack_collab.py",
              "report_files.py", "core_submission.py")


class NoTopLevelSafeWriteImport(unittest.TestCase):
    def test_none_imports_safe_write_at_module_level(self):
        for name in STANDALONE:
            with open(os.path.join(SCRIPTS, name), encoding="utf-8") as fh:
                text = fh.read()
            with self.subTest(name):
                self.assertIsNone(re.search(r"^(from safe_write import|import safe_write)", text, re.M))


if __name__ == "__main__":
    unittest.main()
