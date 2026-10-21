#!/usr/bin/env python3
"""A sheet's own condition names are never merged by taking a "replicate index" off them (Core
review of the client IP-MS guide, 2026-10-20).

Before: submission_report's REPLICATE_SUFFIX read ATG5 as ATG replicate 5 and HSP90 as HSP
replicate 90, so a sheet with ATG5 x3, ATG7 x3 and WT x3 became "ATG" and "WT", and the report
said the analysed groups split the sheet's "ATG" ("conditions differ"), which was false. Now
(sheet_groups, _own_conditions_apart): two names that several samples carry each are never
merged into one base when the numbers taken off are glued to letters and look like gene numbers,
not sample or replicate numbers (ATG 5, 7) -- every name with that base keeps its name.
Everything else reads exactly as 2.13.2 read it (test_replicate_parity holds that on a corpus,
against 2.13.2's answers and the truth); the sheets in AsTwoThirteenTwo are the ones the reviews
of broader rules found broken. Known limits, as 2.13.2: one-sample names (ATG5 and ATG7 alone
are "ATG"), numbers after a separator (Mutant-1 x3 / Mutant-2 x3 are "Mutant"), and a sheet that
numbers its other samples (ATG5 / ATG7 beside WT_1..3).

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

import submission_report as sr  # noqa: E402
from test_submission_report import fixture, notes_by_id  # noqa: E402


def record(names):
    r = fixture()
    r["submission_data"]["samples"] = [
        {"unique_id": f"SMP{100 + i}", "sample_name": f"sample {i}", "condition_name": c}
        for i, c in enumerate(names, 1)]
    return sr.sanitize(r)


def groups(names):
    return set(sr.sheet_groups(record(names)).values())


# (the sheet's condition names, the groups the analysis used -- written by hand, one per sample,
# 2.13.2's sheet groups)
SAME_AS_2_13_2 = (
    ([f"Cis_{x}_{n}" for x in "RS" for n in (1, 2, 3)], ["Cis_R"] * 3 + ["Cis_S"] * 3,
     {"Cis_R", "Cis_S"}),
    ([f"Brain_{x}-{n}" for x in "RL" for n in (1, 2, 3)], ["Brain_R"] * 3 + ["Brain_L"] * 3,
     {"Brain_R", "Brain_L"}),
    ([f"Brain_{x} {n}" for x in "RL" for n in (1, 2, 3)], ["Brain_R"] * 3 + ["Brain_L"] * 3,
     {"Brain_R", "Brain_L"}),
    (["WT1", "WT2", "WT3", "KO", "KO", "KO"], ["WT"] * 3 + ["KO"] * 3, {"WT", "KO"}),
    (["Ctrl"] * 3 + ["Treated1", "Treated2", "Treated3"], ["Ctrl"] * 3 + ["Treated"] * 3,
     {"Ctrl", "Treated"}),
    (["HSP90"] * 3 + ["WT1", "WT2", "WT3"], ["HSP90"] * 3 + ["WT"] * 3, {"HSP", "WT"}),
    (["WT1", "WT2", "WT2", "WT3", "KO1", "KO2", "KO3"], ["WT"] * 4 + ["KO"] * 3, {"WT", "KO"}),
    (["WT_1", "WT_2", "WT_2", "WT_3", "KO1", "KO2", "KO3"], ["WT"] * 4 + ["KO"] * 3, {"WT", "KO"}),
    (["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["Mutant-3", "Mutant-4"],
     ["Mutant_1"] * 3 + ["Mutant_2"] * 3 + ["Mutant_3", "Mutant_4"],
     {"Mutant-1", "Mutant-2", "Mutant-3", "Mutant-4"}),
    ([f"{g}{n}" for g in ("WT", "KO") for n in (2, 3, 4)], ["WT"] * 3 + ["KO"] * 3, {"WT", "KO"}),
)


class TheReportedCase(unittest.TestCase):
    def test_two_conditions_the_sheet_repeats_are_not_merged(self):
        for names, want in (
                (["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3, {"ATG5", "ATG7", "WT"}),
                (["HSP70"] * 3 + ["HSP90"] * 3 + ["WT"] * 3, {"HSP70", "HSP90", "WT"}),
                (["ATG5"] * 3 + ["ATG7"] * 3, {"ATG5", "ATG7"}),
                # every name with the merged base keeps its own
                (["ATG5"] * 3 + ["ATG7"] * 3 + ["ATG8", "WT", "WT"],
                 {"ATG5", "ATG7", "ATG8", "WT"})):
            self.assertEqual(groups(names), want, names)

    def test_numbers_after_a_separator_are_never_kept(self):
        """Sample numbers, not gene numbers: as 2.13.2 read them, listed twice or not."""
        self.assertEqual(groups([x for x in ("Parental_12", "Parental_51") for _ in (0, 1)]
                                + ["Dox(+)"] * 4), {"Parental", "Dox(+)"})

    def test_every_replicate_listed_twice_reads_as_2_13_2(self):
        """IP-MS sheets often list each replicate twice (IP and input rows): a run is replicates."""
        for names in (["WT1"] * 2 + ["WT2"] * 2 + ["KO"] * 2,
                      ["WT1"] * 2 + ["WT2"] * 2 + ["KO1"] * 2 + ["KO2"] * 2,
                      ["WT_1"] * 2 + ["WT_2"] * 2 + ["KO_1"] * 2 + ["KO_2"] * 2):
            self.assertEqual(groups(names), {"WT", "KO"}, names)

    def test_known_limits_read_as_before(self):
        """As 2.13.2 read them: no name repeated (nothing says ATG5 and ATG7 are conditions); a
        number after a separator (Mutant-1, -2); a sheet that numbers its other samples, so its
        numbers are sample numbers too (WT_1..3 beside ATG5 / ATG7)."""
        self.assertEqual(groups(["ATG5", "ATG7", "WT_1", "WT_2"]), {"ATG", "WT"})
        self.assertEqual(groups(["Mutant-1"] * 3 + ["Mutant-2"] * 3 + ["WT"] * 3),
                         {"Mutant", "WT"})
        self.assertEqual(groups(["ATG5"] * 3 + ["ATG7"] * 3 + ["WT_1", "WT_2", "WT_3"]),
                         {"ATG", "WT"})


class AsTwoThirteenTwo(unittest.TestCase):
    def test_the_sheets_reviews_found_read_as_2_13_2(self):
        for names, _analysed, want in SAME_AS_2_13_2:
            self.assertEqual(groups(names), want, names)


class NoFalseConditionsDiffer(unittest.TestCase):
    """End to end: the analysed groups (written by hand) match the sheet, so no "conditions
    differ" note."""

    def session(self, tmp, names, group=lambda c: c.replace("-", "_")):
        """`group`: the analysed group of each sample -- a list in the sheet's order, or a
        function of its condition name."""
        s = os.path.join(tmp, "2026-10-20_conditions")
        os.makedirs(os.path.join(s, "input"))
        os.makedirs(os.path.join(s, "output"))
        rec = record(names)
        runs = {s_["unique_id"]: f"10202026_DIA_{s_['unique_id']}_1_{i}"
                for i, s_ in enumerate(rec["samples"], 1)}
        analysed = group if isinstance(group, list) else [group(c) for c in names]
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw\n" + "".join(f"/data/OCT26/{r}.d\n" for r in runs.values()))
        with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["File.Name", "Group"])
            for smp, g in zip(rec["samples"], analysed):
                w.writerow([runs[smp["unique_id"]], g])
        with open(os.path.join(s, "input", "locate_decisions.json"), "w") as fh:
            json.dump({"assignments": [{"sample": u, "file": f"{r}.d"} for u, r in runs.items()]},
                      fh)
        sr.attach(s, rec)
        return s

    def note(self, names, analysed):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, names, analysed)
            return notes_by_id(sr.load(s), s).get("conditions_differ")

    def test_the_reported_sheet(self):
        names = ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3
        self.assertIsNone(self.note(names, list(names)))

    def test_the_sheets_reviews_found(self):
        for names, analysed, _want in SAME_AS_2_13_2:
            self.assertIsNone(self.note(names, analysed), names)

    def test_every_replicate_listed_twice_gets_no_note(self):
        names = ["WT_1"] * 2 + ["WT_2"] * 2 + ["KO_1"] * 2 + ["KO_2"] * 2
        self.assertIsNone(self.note(names, ["WT"] * 4 + ["KO"] * 4))

    def test_a_real_difference_is_still_reported(self):
        names = ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3
        self.assertIn("pools the sheet's ATG5, ATG7", self.note(names, ["ATG"] * 6 + ["WT"] * 3))


if __name__ == "__main__":
    unittest.main()
