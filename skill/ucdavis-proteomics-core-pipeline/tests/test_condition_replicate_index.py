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
         "Old - A - Bio 2": "Old - A", "WT - Biological 2": "WT", "WT - Tech 2": "WT",
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

    def test_an_r_beside_another_letter_is_a_level(self):
        """Brain_L / Brain_R, Cis_R / Cis_S: left and right, resistant and sensitive."""
        for a, b in (("Brain_L", "Brain_R"), ("Cis_R", "Cis_S")):
            for sep in ("", "_", "-", " "):        # Cis_R1, Cis_R_1, Cis_R-1, "Cis_R 1"
                labels = [f"{g}{sep}{n}" for g in (a, b) for n in (1, 2, 3)]
                rec = cc.collapse_replicates([{"sample": f"S{i}", "group": g, "extras": {}}
                                              for i, g in enumerate(labels)])[1]
                self.assertEqual(sorted(rec["conditions"]), [a, b], labels)
                self.assertEqual(set(sr.sheet_groups(record(labels)).values()), {a, b}, labels)
        self.assertEqual(cc.split_replicate("ATG5_R2", among=["ATG5_R1", "ATG7_R1"]), ("ATG5", 2))

    def test_the_known_condition_list_decides(self):
        self.assertEqual(sr.replicate_base("Mutant-1"), "Mutant")
        self.assertIsNone(sr.replicate_base("Mutant-1", {"mutant-1", "mutant-2"}))
        self.assertEqual(sr.replicate_base("Mutant-1_2", {"mutant-1"}), "Mutant-1")
        self.assertEqual(sr.known_conditions(["Mutant-1", "Mutant-1", "WT_1", "WT_2"]),
                         {"mutant-1"})
        # a slip in a numbered series is not a condition of its own; a third level is no slip
        self.assertEqual(sr.known_conditions(["WT_1", "WT_2", "WT_2", "WT_3"]), set())
        self.assertEqual(sr.known_conditions(["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["Mutant-3"]),
                         {"mutant-1", "mutant-2"})


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
        # mixed formats, and case the sheet check ignores anyway
        for names in ([f"WT{n}" for n in (1, 2, 3)] + [f"KO_{n}" for n in (1, 2, 3)],
                      ["wt1", "wt2", "WT3", "KO1", "KO2", "KO3"]):
            rec = record(names)
            self.assertEqual({sr._key(g) for g in sr.sheet_groups(rec).values()}, {"wt", "ko"})
            self.assertNotIn("sheet_conditions_unique", notes_by_id(rec), names)
        # as 2.13.2 read them: a gap is a missing replicate; gene names beside a series keep theirs
        for names, groups in (
                (["WT1", "WT2", "WT3", "KO1", "KO2", "KO4"], {"wt", "ko"}),
                (["HSP70", "HSP90", "WT1", "WT2"], {"hsp70", "hsp90", "wt"})):
            self.assertEqual({sr._key(g) for g in sr.sheet_groups(record(names)).values()}, groups)
        # numbered conditions the sheet repeats are not merged, with a third level beside them
        groups = set(sr.sheet_groups(record(["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["Mutant-3"]))
                     .values())
        self.assertTrue({"Mutant-1", "Mutant-2"} <= groups, groups)
        # a numbered series is not a set of replicates: Day1..Day4 stay four conditions
        rec = record(["Day1", "Day2", "Day3", "Day4"])
        self.assertIn("sheet_conditions_unique", notes_by_id(rec))

    def test_a_slip_in_a_numbered_series_joins_its_condition(self):
        for names in (["WT_1", "WT_2", "WT_2", "WT_3", "KO_1", "KO_2", "KO_3"],
                      ["WT1", "WT2", "WT2", "WT3", "KO1", "KO2", "KO3"],
                      ["WT_1", "WT_2", "WT_2", "WT_3", "KO1", "KO2", "KO3"]):
            self.assertEqual(set(sr.sheet_groups(record(names)).values()), {"WT", "KO"}, names)
        # a gene name two samples share, beside two other genes, keeps it (start: ATG 5, 7, 8)
        self.assertEqual(set(sr.sheet_groups(record(["ATG5", "ATG5", "ATG7", "ATG8"])).values()),
                         {"ATG5", "ATG7", "ATG8"})

    def test_a_repeated_name_leaves_the_rest_of_the_sheet_as_2_13_2_read_it(self):
        for names, groups in (
                (["WT1", "WT2", "WT3", "KO", "KO", "KO"], {"WT", "KO"}),
                (["Ctrl"] * 3 + ["Treated1", "Treated2", "Treated3"], {"Ctrl", "Treated"}),
                (["HSP90"] * 3 + ["WT1", "WT2", "WT3"], {"HSP90", "WT"}),
                (["WT"] * 3 + ["WT4"] + ["KO"] * 3, {"WT", "KO"}),
                # the series of a numbered condition the sheet repeats stays apart
                (["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["Mutant-3", "Mutant-4"],
                 {"Mutant-1", "Mutant-2", "Mutant-3", "Mutant-4"}),
                # deliberate (2.13.2 merged them): repeated numbered names are conditions as written
                (["WT1"] * 2 + ["WT2"] * 2, {"WT1", "WT2"})):
            self.assertEqual(set(sr.sheet_groups(record(names)).values()), groups, names)

    def test_numbering_that_does_not_start_at_one_is_not_a_series(self):
        """Kept on purpose (2.13.2 read WT/KO): numbers from 2 read like gene names."""
        self.assertIn("sheet_conditions_unique",
                      notes_by_id(record([f"{g}{n}" for g in ("WT", "KO") for n in (2, 3, 4)])))

    def test_results_that_changed_for_the_better(self):
        """2.13.2 said "no groups" for these sheets, and proposed "WT_R", "WTrep" as conditions."""
        for names, groups in (([f"{g}_{n}" for g in ("T0", "T24") for n in (1, 2)], {"T0", "T24"}),
                              ([f"{g} \u2013 {n}" for g in ("WT", "KO") for n in (1, 2)],
                               {"WT", "KO"})):
            self.assertEqual(set(sr.sheet_groups(record(names)).values()), groups, names)
        # never a number after a dot, as 2.13.2: Kv2.1 / Kv2.2 are two channels; WT.1 stays WT.1
        for names, groups in ((["Kv2.1", "Kv2.2", "WT_1", "WT_2"], {"Kv2.1", "Kv2.2", "WT"}),
                              (["WT.1", "WT.2", "KO.1", "KO.2"], {"WT.1", "WT.2", "KO.1", "KO.2"})):
            self.assertEqual(set(sr.sheet_groups(record(names)).values()), groups, names)
        for fmt in ("{g}_R{n}", "{g}rep{n}", "{g}#{n}", "{g} - Bio {n}"):
            labels = [fmt.format(g=g, n=n) for g in ("WT", "KO") for n in (1, 2)]
            rec = cc.collapse_replicates([{"sample": f"S{i}", "group": g, "extras": {}}
                                          for i, g in enumerate(labels)])[1]
            self.assertEqual((rec["collapsed"], sorted(rec["conditions"])), (True, ["KO", "WT"]),
                             fmt)


class NoFalseConditionsDiffer(unittest.TestCase):
    """End to end: the analysed groups match the sheet, so no "conditions differ" note."""

    def session(self, tmp, names, group=lambda c: c.replace("-", "_")):
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
                w.writerow([runs[smp["unique_id"]], group(smp["condition_name"])])
        with open(os.path.join(s, "input", "locate_decisions.json"), "w") as fh:
            json.dump({"assignments": [{"sample": u, "file": f"{r}.d"} for u, r in runs.items()]},
                      fh)
        sr.attach(s, rec)
        return s

    def test_gene_named_groups(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3)
            self.assertNotIn("conditions_differ", notes_by_id(sr.load(s), s))

    def test_mixed_formats_and_left_right_levels(self):
        for names in ([f"WT{n}" for n in (1, 2, 3)] + [f"KO_{n}" for n in (1, 2, 3)],
                      [f"Brain_{s}{n}" for s in "LR" for n in (1, 2, 3)]):
            with tempfile.TemporaryDirectory() as tmp:
                s = self.session(tmp, names, group=lambda c: cc.split_replicate(
                    c, glued=True, among=names)[0])
                self.assertNotIn("conditions_differ", notes_by_id(sr.load(s), s), names)

    def test_a_repeated_name_beside_numbered_ones(self):
        for names, analysed in (
                (["WT1", "WT2", "WT3", "KO", "KO", "KO"], None),
                (["WT1", "WT2", "WT2", "WT3", "KO1", "KO2", "KO3"], None),
                (["WT_1", "WT_2", "WT_2", "WT_3", "KO1", "KO2", "KO3"], None),
                (["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["Mutant-3", "Mutant-4"],
                 lambda c: c.replace("-", "_"))):
            with tempfile.TemporaryDirectory() as tmp:
                s = self.session(tmp, names, group=analysed or (
                    lambda c: (cc.split_replicate(c, glued=True) or (c,))[0].rstrip("_")))
                self.assertNotIn("conditions_differ", notes_by_id(sr.load(s), s), names)

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
