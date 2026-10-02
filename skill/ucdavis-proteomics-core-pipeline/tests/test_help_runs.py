#!/usr/bin/env python3
"""`--help` exits 0, for every subcommand of the scripts the 2.10 FASTA / submission / repeats /
re-injection work changed.

2.10 safety review: `fetch_fasta.py --help` and `fetch --help` crashed (exit 1, "TypeError: %o
format") -- an f-string put "50% of" into argparse help, which %-formats it. A user (or agent)
asking how to use the script got a traceback. Subcommands are read from each script's own
--help ({a,b,...}), so a new one is covered without editing this list.
"""
import os
import re
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")

TOUCHED = ["fetch_fasta.py", "core_submission.py", "collect_conditions.py", "ht_manifest.py",
           "run_search.py", "diann_parallel.py", "radiant_parallel.py", "check_report_runs.py",
           "make_methods.py", "make_analysis_html.py", "provenance.py", "session.py",
           "submission_report.py", "experiment_type.py"]
_SUBS = re.compile(r"^\s*\{([A-Za-z0-9_,-]+)\}", re.M)
# fetch_fasta.py reads a bare flag as its old flat CLI (= `fetch`), so its top-level --help is
# fetch's and lists no subcommands: named here.
KNOWN_SUBS = {"fetch_fasta.py": ["resolve", "check-db", "keratin-db", "fetch"]}


def run_help(script, *sub):
    env = {k: v for k, v in os.environ.items() if k not in ("CLAUDE_CODE_SESSION_ID",
                                                             "ANTHROPIC_API_KEY")}
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *sub, "--help"],
                          capture_output=True, text=True, timeout=60, env=env)


class HelpExitsZero(unittest.TestCase):
    def test_every_subcommand_of_every_touched_script(self):
        for script in TOUCHED:
            with self.subTest(script=script):
                p = run_help(script)
                self.assertEqual(p.returncode, 0, p.stderr[-800:])
                self.assertIn("usage:", p.stdout)
                subs = [x for g in _SUBS.findall(p.stdout) for x in g.split(",")]
                for sub in dict.fromkeys(subs + KNOWN_SUBS.get(script, [])):
                    with self.subTest(script=script, sub=sub):
                        q = run_help(script, sub)
                        self.assertEqual(q.returncode, 0, q.stderr[-800:])
                        self.assertIn("usage:", q.stdout)

    def test_the_reported_crash(self):
        for sub in ((), ("fetch",)):
            p = run_help("fetch_fasta.py", *sub)
            self.assertEqual(p.returncode, 0, p.stderr[-800:])
            self.assertNotIn("Traceback", p.stderr)


if __name__ == "__main__":
    unittest.main()
