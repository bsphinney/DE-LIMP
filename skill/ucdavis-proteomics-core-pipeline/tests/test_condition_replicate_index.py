#!/usr/bin/env python3
"""A condition name keeps the digits that are part of it (Core review of the client IP-MS guide,
2026-10-20).

Before: submission_report's REPLICATE_SUFFIX took any trailing number after a letter as a
replicate index, so the sheet check read ATG5 as ATG replicate 5 and HSP90 as HSP replicate 90.
A sheet with ATG5 and ATG7 groups became one "ATG" condition, and the report said the analysed
groups split it ("conditions differ"), which was false. collect_conditions.split_replicate was a
second copy of the same rule. Now there is one definition (collect_conditions.split_replicate):
only a MARKED index is a replicate (a separator then digits, or rep / replicate / R / #), and a
word that merely ends in "rep" is no marker (Strep_2, SamplePrep_2). Numbers glued to a name
(WT1..3, Ctl1..3) are read as replicates only where collect_conditions.replicate_reading's checks
pass -- collapse_replicates, which also asks the user, and the sheet check when the marked reading
leaves every sample its own condition -- so ATG5 / ATG7 / CD45 / HSP90 keep their names. A name
several samples carry is a condition as written, whatever its ending, unless it is a slip in a
numbered series (WT_1, WT_2, WT_2, WT_3).

Stdlib only.
"""
import csv
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import collect_conditions as cc  # noqa: E402
import submission_report as sr  # noqa: E402
from test_submission_report import fixture, notes_by_id  # noqa: E402

KEEP = ("ATG5", "HSP90", "CD45", "IgG1", "ATG5_IP", "CXCR4", "Kv2.1", "Prep2", "Strep2",
        "Biorep2")
STRIP = {"ATG5_2": "ATG5", "ATG5-rep2": "ATG5", "ATG5_R2": "ATG5", "WT_2": "WT", "KO-3": "KO",
         "Ctrl rep 2": "Ctrl", "WT#2": "WT", "WTrep2": "WT", "KORep2": "KO", "WT_replicate_3": "WT",
         "WT_BioRep2": "WT", "WT - Tech Rep 2": "WT", "WT \u2013 2": "WT", "WT \u20132": "WT",
         "Control 3": "Control", "Old - A - Rep 2": "Old - A", "Old - A - Bio Rep 2": "Old - A",
         # a word that ends in "rep" keeps it
         "Strep_2": "Strep", "SamplePrep_2": "SamplePrep", "NoPrep-1": "NoPrep",
         "FastPrep 3": "FastPrep"}


class OneMarkedIndex(unittest.TestCase):
    def test_digits_that_are_part_of_the_name_are_kept(self):
        for name in KEEP:
            self.assertIsNone(sr.replicate_base(name), name)
            self.assertIsNone(cc.split_replicate(name), name)

    def test_a_marked_index_is_stripped(self):
        for name, base in STRIP.items():
            self.assertEqual(sr.replicate_base(name), base, name)

    def test_the_glued_reading_is_only_for_the_asked_collapse(self):
        self.assertIsNone(cc.split_replicate("Control1"))
        self.assertEqual(cc.split_replicate("Control1", glued=True), ("Control", 1))
        self.assertEqual(cc.split_replicate("CtrlA_mix_3", glued=True), ("CtrlA_mix", 3))

    def test_the_known_condition_list_decides(self):
        self.assertEqual(sr.replicate_base("Mutant-1"), "Mutant")
        self.assertIsNone(sr.replicate_base("Mutant-1", {"mutant-1", "mutant-2"}))
        self.assertEqual(sr.replicate_base("Mutant-1_2", {"mutant-1"}), "Mutant-1")
        self.assertEqual(sr.known_conditions(["Mutant-1", "Mutant-1", "WT_1", "WT_2"]), {"mutant-1"})
        # a slip in a numbered series is not a condition of its own
        self.assertEqual(sr.known_conditions(["WT_1", "WT_2", "WT_2", "WT_3"]), set())


def record(names):
    r = fixture()
    r["submission_data"]["samples"] = [
        {"unique_id": f"SMP{100 + i}", "sample_name": f"sample {i}", "condition_name": c}
        for i, c in enumerate(names, 1)]
    return sr.sanitize(r)


class TheSheetsGroups(unittest.TestCase):
    def test_gene_named_groups_stay_apart(self):
        rec = record(["ATG5"] * 3 + ["ATG7"] * 3 + ["HSP90"] * 3)
        self.assertEqual(set(sr.sheet_groups(rec).values()), {"ATG5", "ATG7", "HSP90"})

    def test_replicate_suffixes_after_a_gene_name_collapse(self):
        rec = record([f"{g}_{n}" for g in ("ATG5", "ATG7") for n in (1, 2, 3)])
        self.assertEqual(set(sr.sheet_groups(rec).values()), {"ATG5", "ATG7"})

    def test_names_several_samples_carry_are_the_conditions(self):
        rec = record(["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["WT"] * 3)
        self.assertEqual(set(sr.sheet_groups(rec).values()), {"Mutant-1", "Mutant-2", "WT"})

    def test_a_sheet_of_gene_names_with_one_sample_each_has_no_groups(self):
        rec = record(["ATG5", "ATG7", "CD45", "HSP90"])
        self.assertIn("sheet_conditions_unique", notes_by_id(rec))

    def test_numbers_glued_to_the_names_are_replicates_when_the_checks_pass(self):
        for names, groups in (([f"{g}{n}" for g in ("WT", "KO") for n in (1, 2, 3)], {"WT", "KO"}),
                              ([f"{g}{n}" for g in ("WT", "KO") for n in (1, 2, 3)] + ["Pool"],
                               {"WT", "KO", "Pool"})):
            rec = record(names)
            self.assertEqual(set(sr.sheet_groups(rec).values()), groups, names)
            self.assertNotIn("sheet_conditions_unique", notes_by_id(rec), names)
        # a numbered series is not a set of replicates: Day1..Day4 stay four conditions
        rec = record(["Day1", "Day2", "Day3", "Day4"])
        self.assertIn("sheet_conditions_unique", notes_by_id(rec))

    def test_a_slip_in_a_numbered_series_joins_its_condition(self):
        rec = record(["WT_1", "WT_2", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"])
        self.assertEqual(set(sr.sheet_groups(rec).values()), {"WT", "KO"})


class NoFalseConditionsDiffer(unittest.TestCase):
    """End to end: the analysed groups match the sheet, so no "conditions differ" note."""

    def session(self, tmp, names):
        s = os.path.join(tmp, "2026-10-20_conditions")
        os.makedirs(os.path.join(s, "input"))
        os.makedirs(os.path.join(s, "output"))
        rec = record(names)
        runs = {s_["unique_id"]: f"10202026_DIA_{s_['unique_id']}_1_{i}"
                for i, s_ in enumerate(rec["samples"], 1)}
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw\n" + "".join(f"/data/OCT26/{r}.d\n" for r in runs.values()))
        with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["File.Name", "Group"])
            for smp in rec["samples"]:
                w.writerow([runs[smp["unique_id"]], smp["condition_name"].replace("-", "_")])
        with open(os.path.join(s, "input", "locate_decisions.json"), "w") as fh:
            json.dump({"assignments": [{"sample": u, "file": f"{r}.d"} for u, r in runs.items()]},
                      fh)
        sr.attach(s, rec)
        return s

    def test_gene_named_groups(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3)
            self.assertNotIn("conditions_differ", notes_by_id(sr.load(s), s))

    def test_numbered_conditions_the_sheet_repeats(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, ["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["WT"] * 3)
            self.assertNotIn("conditions_differ", notes_by_id(sr.load(s), s))

    def test_a_real_difference_is_still_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3)
            p = os.path.join(s, "input", "conditions.csv")
            with open(p, newline="") as fh:
                rows = list(csv.reader(fh))
            for r in rows[1:]:
                if r[1] in ("ATG5", "ATG7"):
                    r[1] = "ATG"
            with open(p, "w", newline="") as fh:
                csv.writer(fh).writerows(rows)
            self.assertIn("pools the sheet's ATG5, ATG7",
                          notes_by_id(sr.load(s), s)["conditions_differ"])


if __name__ == "__main__":
    unittest.main()
