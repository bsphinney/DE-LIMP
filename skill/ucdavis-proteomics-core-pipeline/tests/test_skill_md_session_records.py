#!/usr/bin/env python3
"""SKILL.md puts the session records in an order that lets them work (2.10 safety review).

  * step 3b attaches the CoreOmics submission BEFORE `experiment_type.py set`: set with
    --source submission reads the attached record (its ASK compares the submission's
    normalisation with the type's default), and set first meant the ASK never fired;
  * the Core flow's copy-back list brings the records the local report reads:
    sample_type.json, experiment_type.json and locate_decisions.json;
  * after step 6 the bait is asked explicitly (ask_bait) and recorded with set-bait.
"""
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.join(os.path.dirname(HERE), "SKILL.md")


def text():
    with open(SKILL, encoding="utf-8") as fh:
        return fh.read()


class SessionRecords(unittest.TestCase):
    def test_step_3b_attaches_before_the_experiment_type(self):
        s = text()
        b = s.index("### 3b. Create the analysis session")
        attach = s.index("python3 scripts/submission_report.py attach --session <session>", b)
        etype = s.index("python3 scripts/experiment_type.py set --session <session>", b)
        self.assertLess(attach, etype)

    def test_the_copy_back_list_brings_the_session_records(self):
        line = next(ln for ln in text().splitlines()
                    if ln.strip().startswith("for f in submission.json samples.tsv"))
        for f in ("sample_type.json", "experiment_type.json", "locate_decisions.json"):
            self.assertIn(f, line)

    def test_the_bait_is_asked_after_step_6_and_recorded(self):
        s = text()
        step6, step7 = s.index("### 6. Build the FASTA"), s.index("### 7. Run the search")
        part = s[step6:step7]
        self.assertIn("ask_bait", part)
        self.assertTrue(re.search(r"experiment_type\.py set-bait --session <session> --bait", part))


if __name__ == "__main__":
    unittest.main()
