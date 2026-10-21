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
        """The pull list lives in scripts/core_sync.sh since the slim-down (SKILL.md step 1c.7
        and core-submissions.md call it instead of a hand-copied loop)."""
        self.assertIn("core_sync.sh --get --session", " ".join(text().split()))
        with open(os.path.join(os.path.dirname(SKILL), "scripts", "core_sync.sh"),
                  encoding="utf-8") as fh:
            sync = fh.read()
        records = re.search(r'GET_RECORDS="([^"]*)"', sync).group(1)
        for f in ("sample_type.json", "experiment_type.json", "locate_decisions.json"):
            self.assertIn("input/" + f, records.split())

    def test_the_bait_is_asked_after_step_6_and_recorded(self):
        s = text()
        step6, step7 = s.index("### 6. Build the FASTA"), s.index("### 7. Run the search")
        part = s[step6:step7]
        self.assertIn("ask_bait", part)
        self.assertTrue(re.search(r"experiment_type\.py set-bait --session <session> --bait", part))


if __name__ == "__main__":
    unittest.main()


class CallIdIsExplained(unittest.TestCase):
    def test_step_3b_says_when_and_how_to_give_a_call_id(self):
        """Review of 2.13.1, round 4: the token is needed whenever a subagent ran in the last
        6 hours (finalize after 8d's reviewers too), typed literally, fresh, in the foreground."""
        s = text()
        router = " ".join(s[s.index("### 3b."):s.index("### 4.")].split())
        self.assertIn("When a subagent ran in the last 6 hours, `--call-id <token>` on `init`, "
                      "`finalize`, `checkpoint.py status` and `save_transcript.py`", router)
        for want in ("typed literally", "fresh each call", "in the foreground"):
            self.assertIn(want, router)
        # and the finalize line carries it (review of the slim-down, H3)
        twelve = " ".join(s[s.index("### 12. "):s.index("### 12b.")].split())
        self.assertIn("finalize --dir <session> --zip [--call-id <token>]", twelve)
        with open(os.path.join(os.path.dirname(SKILL), "references", "outputs.md"),
                  encoding="utf-8") as fh:
            ref = fh.read()
        b = ref[ref.index("## 3b. Creating the analysis session"):
                ref.index("**Core data with a CoreOmics submission")]
        for want in ("`--call-id` whenever any subagent or teammate of this conversation ran in "
                     "the last 6 hours", "never `$(uuidgen)`", "fresh for every call",
                     "a background call is never counted"):
            self.assertIn(want, " ".join(b.split()))


class ReviewsAreKept(unittest.TestCase):
    def test_step_8d_reviewers_write_their_findings_to_logs(self):
        """2026-10: three step-8d reviewers, run from a teammate, all ended without handing
        their report back (twice); told to also write it to a file, they delivered."""
        s = text()
        d8 = " ".join(s[s.index("### 8d."):s.index("### 8e.")].split())     # layout-free
        self.assertIn("<S>/logs/review_<proteomics|biology|statistics>.md", d8)
        self.assertIn("never delivered", d8)
        self.assertIn("hive_remote: locally, then `hive_exec.sh --put`", d8)
