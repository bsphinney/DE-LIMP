#!/usr/bin/env python3
"""reproduce.sh replays the contrasts the DE ran, from the DE's own record (C1 review, F2).

reproduce.sh's step 6 took its --contrasts only from `provenance.py --contrasts`, which SKILL.md
step 10 shows as a placeholder. 2.14.2 made the default contrasts "group minus control" when
controls are recorded (T2: Bait-IgG), while run_de.R's own default is every group against the
first level alphabetically (IgG-Bait). So a blank --contrasts replayed IgG-Bait, and a stale one
replayed whatever it said: every fold change's sign flipped, silently.

Now the contrasts come from the DE record (de_provenance.json `contrasts`, what ran), as the
between-run step does (L1-F3: replay what ran). A --contrasts that differs from the record is not
replayed; MANIFEST.txt has a NOTE naming both. run_manifest.json and REPRODUCE.md say the same
contrasts as reproduce.sh. A DE record with no contrasts (an older run_de.R) falls back to
--contrasts; with neither, the replay's run_de.R default is said in a NOTE.

Stdlib only: reproduce.sh's steps 5 and 6 run with stub run_search.py / Rscript.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from test_reproduce_sh_runs import arg, provenance, run_steps_5_and_6  # noqa: E402


class ReplayWhatRan(unittest.TestCase):
    def bundle(self, recorded, *given):
        """-> (the replayed --contrasts, MANIFEST.txt, run_manifest.json, REPRODUCE.md)."""
        with tempfile.TemporaryDirectory() as d:
            de = os.path.join(d, "tables")
            os.makedirs(de)
            rec = {"method": "dpc", "quantities": "normalised"}
            if recorded is not None:
                rec["contrasts"] = recorded
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump(rec, fh)
            out = os.path.join(d, "repro")
            sh = provenance(out, de, "--de-method", "dpc", *given)
            calls, _ = run_steps_5_and_6(d, out, sh)
            with open(os.path.join(out, "MANIFEST.txt")) as fh:
                manifest = fh.read()
            with open(os.path.join(out, "run_manifest.json")) as fh:
                run_manifest = json.load(fh)
            with open(os.path.join(out, "REPRODUCE.md")) as fh:
                md = fh.read()
        de_call = calls[1]["argv"]
        self.assertTrue(de_call[0].endswith("/scripts/run_de.R"), de_call)
        return arg(de_call, "--contrasts"), manifest, run_manifest, md

    def test_no_contrasts_given_replays_the_record(self):
        flag, manifest, rm, md = self.bundle(["Bait-IgG"])
        self.assertEqual(flag, "Bait-IgG")
        self.assertEqual(rm["de"]["contrasts"], "Bait-IgG")
        self.assertIn("contrasts=`Bait-IgG`", md)
        self.assertFalse([ln for ln in manifest.splitlines()
                          if ln.startswith("[NOTE]") and "contrasts" in ln], manifest)

    def test_matching_contrasts_are_replayed_without_a_note(self):
        flag, manifest, rm, _ = self.bundle(["B-A", "C-A"], "--contrasts", "B-A, C-A")
        self.assertEqual(flag, "B-A,C-A")
        self.assertEqual(rm["de"]["contrasts"], "B-A,C-A")
        self.assertFalse([ln for ln in manifest.splitlines()
                          if ln.startswith("[NOTE]") and "contrasts" in ln], manifest)

    def test_stale_contrasts_are_not_replayed_and_the_note_names_both(self):
        flag, manifest, rm, md = self.bundle(["Bait-IgG"], "--contrasts", "IgG-Bait")
        self.assertEqual(flag, "Bait-IgG", "the record is replayed, not the stale value")
        self.assertEqual(rm["de"]["contrasts"], "Bait-IgG")
        self.assertIn("contrasts=`Bait-IgG`", md)
        notes = [ln for ln in manifest.splitlines()
                 if ln.startswith("[NOTE]") and "contrasts" in ln]
        self.assertEqual(len(notes), 1, manifest)
        self.assertIn("IgG-Bait", notes[0])
        self.assertIn("Bait-IgG", notes[0])

    def test_a_one_contrast_record_written_as_a_bare_string_is_read(self):
        self.assertEqual(self.bundle("Bait-IgG")[0], "Bait-IgG")

    def test_a_record_without_contrasts_falls_back_to_the_flag(self):
        flag, manifest, rm, _ = self.bundle(None, "--contrasts", "B-A")
        self.assertEqual(flag, "B-A")
        self.assertEqual(rm["de"]["contrasts"], "B-A")
        self.assertTrue([ln for ln in manifest.splitlines()
                         if ln.startswith("[NOTE]") and "--contrasts" in ln], manifest)

    def test_neither_says_the_replay_runs_run_de_defaults(self):
        flag, manifest, _, _ = self.bundle(None)
        self.assertIsNone(flag)
        notes = [ln for ln in manifest.splitlines()
                 if ln.startswith("[NOTE]") and "contrasts" in ln]
        self.assertEqual(len(notes), 1, manifest)
        self.assertIn("default", notes[0])


if __name__ == "__main__":
    unittest.main()
