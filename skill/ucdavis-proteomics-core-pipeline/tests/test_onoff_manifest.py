#!/usr/bin/env python3
"""finalize's MANIFEST.txt carries one line for the MaxLFQ route's on/off list (CLAUDE.md rule 4:
no silent catch in an export path): [OK] with counts when OnOff_Proteins.csv is there, [SKIPPED]
with the reason when run_de.R could not write it / the file is gone / no contrast was a simple A-B,
a [NOTE] per skipped contrast, and nothing for a DE with no `onoff` record (dpc). Python only."""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import make_deposit  # noqa: E402
import session  # noqa: E402


class OnOffManifest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.p = {"session_dir": self.tmp.name, "de_dir": os.path.join(self.tmp.name, "output", "tables")}
        os.makedirs(self.p["de_dir"])

    def record(self, onoff, file=True, key="onoff"):
        prov = {"method": "maxlfq"}
        if onoff is not None or key:
            prov[key] = onoff
        with open(os.path.join(self.p["de_dir"], "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        if file:
            with open(os.path.join(self.p["de_dir"], "OnOff_Proteins.csv"), "w") as fh:
                fh.write("Protein.Group\n")

    def lines(self):
        man = make_deposit.Manifest()
        session.onoff_step(self.p, man)
        return man.lines

    def test_ok_with_counts_and_a_note_per_skipped_contrast(self):
        self.record({"file": "OnOff_Proteins.csv", "status": "written", "n_min": 2,
                     "contrasts": {"B-A": {"n_onoff": 7}}, "skipped": {"C-(A+B)/2": "not a simple A-B"}})
        ls = self.lines()
        self.assertTrue(ls[0].startswith("[OK]") and "B-A 7" in ls[0] and "min. 2 runs" in ls[0], ls)
        self.assertTrue(ls[1].startswith("[NOTE]") and "C-(A+B)/2" in ls[1], ls)

    def test_failure_is_skipped_with_the_reason(self):
        self.record({"file": "OnOff_Proteins.csv", "status": "failed", "reason": "disk full"}, file=False)
        ls = self.lines()
        self.assertEqual(len(ls), 1)
        self.assertTrue(ls[0].startswith("[SKIPPED]") and "disk full" in ls[0], ls)

    def test_missing_file_is_skipped_not_ok(self):
        self.record({"file": "OnOff_Proteins.csv", "status": "written", "n_min": 2,
                     "contrasts": {"B-A": {"n_onoff": 1}}}, file=False)
        ls = self.lines()
        self.assertTrue(ls[0].startswith("[SKIPPED]") and "is not in" in ls[0], ls)

    def test_no_simple_contrast_is_skipped_with_that_reason(self):
        self.record({"file": "OnOff_Proteins.csv", "status": "no_simple_contrasts", "n_min": 2,
                     "contrasts": [], "skipped": {"X": "why"}})
        ls = self.lines()
        self.assertTrue(ls[0].startswith("[SKIPPED]") and "simple A-B" in ls[0], ls)
        self.assertTrue(ls[1].startswith("[NOTE]"), ls)

    def test_a_record_without_onoff_adds_nothing(self):
        self.record(None, file=False)
        self.assertEqual(self.lines(), [])
        os.remove(os.path.join(self.p["de_dir"], "de_provenance.json"))
        self.assertEqual(self.lines(), [])


if __name__ == "__main__":
    unittest.main()
