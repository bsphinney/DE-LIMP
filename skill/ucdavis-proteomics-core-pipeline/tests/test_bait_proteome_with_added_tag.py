#!/usr/bin/env python3
"""A TurboID bait from the proteome survives "none of the added sequences is the bait" (2.14.2, I15).

A proximity study whose bait is a protein of the searched proteome (recorded at step 3 with
`experiment_type.py set --bait VPS35`) and whose added sequence is only the TurboID enzyme had no
stable state: `propose` asked whether an added sequence is the bait; the right answer, `set-bait
--bait none`, ERASED the recorded bait; and a later `set --bait VPS35` dropped that answer, so
`propose` asked again. Now "none" with a recorded proteome bait keeps the bait and records that
the added sequences are not it, nothing is asked again, and `set --bait` with the same bait keeps
the answer. Stdlib only; hermetic.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402
from test_bait_question_many_added import bait_session  # noqa: E402

ETYPE = os.path.join(SCRIPTS, "experiment_type.py")
# the one added sequence: the TurboID enzyme, not the bait
TURBO_ONLY = [("1000001", "IN|1000001|TurboID", "IN|1000001|TurboID biotin ligase GN=birA")]
PROTEOME_SOURCE = "set --bait (a proteome protein); the added sequences are not the bait"


def propose(s):
    p = subprocess.run([sys.executable, ETYPE, "propose", "--session", s], capture_output=True,
                       text=True, timeout=60)
    assert p.returncode == 0, p.stderr
    return json.loads(p.stdout)


class ProteomeBaitWithAnAddedTag(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.s = bait_session(self.d, TURBO_ONLY)
        et.record(self.s, "proximity", "user", bait="VPS35", controls="Ctrl")

    def test_the_question_says_none_keeps_the_recorded_bait(self):
        q = propose(self.s)["ask_bait"]
        self.assertTrue(q.startswith("The bait was recorded at step 3 as VPS35, which names none "
                                     "of the sequences the user added."), q)
        self.assertIn("--bait none keeps VPS35 as the bait", q)

    def test_none_for_added_tag_keeps_proteome_bait(self):
        et.record_bait(self.s, "none")
        rec = et.load(self.s)
        self.assertEqual(rec["bait"], "VPS35")
        self.assertEqual(rec["bait_confirmed"]["answer"], et.NO_BAIT)
        self.assertIs(rec["bait_confirmed"]["added_sequences_are_not_the_bait"], True)
        self.assertTrue(rec["bait_confirmed"]["by"] and rec["bait_confirmed"]["at"])
        # the data check uses the proteome bait
        self.assertEqual(et.bait_for_check(self.s, rec=rec)[:2],
                         ("VPS35", "the experiment-type record"))

    def test_no_question_asked_again(self):
        et.record_bait(self.s, "none")
        out = propose(self.s)
        self.assertIsNone(out["ask_bait"])
        self.assertEqual(out["bait_recorded"], {"bait": "VPS35", "source": PROTEOME_SOURCE})

    def test_set_bait_keeps_confirmation(self):
        et.record_bait(self.s, "none")
        et.record(self.s, "proximity", "user", bait="VPS35", controls="Ctrl")   # re-recorded
        rec = et.load(self.s)
        self.assertEqual(rec["bait"], "VPS35")
        self.assertIs(rec["bait_confirmed"]["added_sequences_are_not_the_bait"], True)
        self.assertIsNone(propose(self.s)["ask_bait"])
        # ... and through the CLI, as step 3 runs it
        p = subprocess.run([sys.executable, ETYPE, "set", "--session", self.s, "--type",
                            "proximity", "--source", "user", "--bait", "VPS35"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIsNone(propose(self.s)["ask_bait"])

    def test_none_first_then_a_proteome_bait_is_that_bait(self):
        """Review of 2.14.2: no bait at step 3, set-bait --bait none, THEN set --bait VPS35. The
        record said VPS35 while propose answered bait "none" and asked nothing -- the agent was
        told there is no bait. The kept answer now says the added sequences are not the bait, and
        propose returns the recorded proteome bait."""
        s = bait_session(self.d, TURBO_ONLY)
        et.record(s, "proximity", "user", controls="Ctrl")
        et.record_bait(s, "none")
        et.record(s, "proximity", "user", bait="VPS35", controls="Ctrl")
        rec = et.load(s)
        self.assertEqual(rec["bait"], "VPS35")
        self.assertIs(rec["bait_confirmed"]["added_sequences_are_not_the_bait"], True)
        out = propose(s)
        self.assertIsNone(out["ask_bait"])
        self.assertEqual(out["bait_recorded"], {"bait": "VPS35", "source": PROTEOME_SOURCE})
        # an older record (the answer without the flag) is read the same way
        with open(et.path_for(s)) as fh:
            raw = json.load(fh)
        del raw["bait_confirmed"]["added_sequences_are_not_the_bait"]
        self.assertEqual(et.bait_question(raw, et.bait_candidates(s))[1],
                         {"bait": "VPS35", "source": PROTEOME_SOURCE})

    def test_set_bait_naming_an_added_sequence_drops_a_contradicted_none(self):
        """(pin) 'none of the added sequences is the bait' does not survive a --bait that names
        one (set --bait always dropped the answer; it still does here)."""
        et.record_bait(self.s, "none")
        et.record(self.s, "proximity", "user", bait="1000001")
        rec = et.load(self.s)
        self.assertIsNone(rec.get("bait_confirmed"))
        self.assertEqual(propose(self.s)["bait_recorded"]["accession"], "1000001")

    def test_set_bait_keeps_an_accession_answer_for_the_same_bait(self):
        s = bait_session(self.d, TURBO_ONLY)
        et.record(s, "proximity", "user")
        et.record_bait(s, "1000001")
        et.record(s, "proximity", "user", bait="1000001")
        self.assertEqual(et.load(s)["bait_confirmed"]["answer"], "1000001")
        et.record(s, "proximity", "user", bait="VPS35")     # another bait: the answer was not about it
        self.assertIsNone(et.load(s).get("bait_confirmed"))

    def test_none_with_no_recorded_bait_still_erases(self):
        """(pin) no bait recorded: "none" leaves no bait, as before."""
        s = bait_session(self.d, TURBO_ONLY)
        et.record(s, "proximity", "user")
        et.record_bait(s, "none")
        rec = et.load(s)
        self.assertIsNone(rec["bait"])
        self.assertEqual(rec["bait_confirmed"]["answer"], et.NO_BAIT)
        self.assertNotIn("added_sequences_are_not_the_bait", rec["bait_confirmed"])
        self.assertEqual(propose(s)["bait_recorded"], {"bait": et.NO_BAIT,
                                                       "source": "set-bait (asked after step 6)"})

    def test_none_for_a_bait_that_is_an_added_sequence_still_erases(self):
        """(pin) the recorded bait IS an added sequence and the user says none is: no bait."""
        s = bait_session(self.d, TURBO_ONLY)
        et.record(s, "proximity", "user", bait="1000001")
        et.record_bait(s, "none")
        self.assertIsNone(et.load(s)["bait"])


if __name__ == "__main__":
    unittest.main()
