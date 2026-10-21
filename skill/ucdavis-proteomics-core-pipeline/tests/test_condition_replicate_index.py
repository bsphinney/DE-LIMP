#!/usr/bin/env python3
"""Names with a number glued to them: a possible difference for staff, never a finding in the
report (Core review of the client IP-MS guide, 2026-10-20).

Before: submission_report's sheet check read ATG5 as ATG replicate 5 and HSP90 as HSP replicate
90, so a sheet with ATG5 x3, ATG7 x3 and WT x3 became "ATG" and "WT", and the report told the
client the analysed groups split the sheet's "ATG" ("conditions differ"), which was false. Three
reviews showed that no rule reading names alone can tell ATG5 / ATG7 (genes) from Parental12 /
Parental51 (sample IDs), so the sheet's groups stay exactly 2.13.2's (test_replicate_parity).
What changed is the note: a part of "Conditions analysed vs the sample sheet" that rests on a
group made by taking a number glued to the names off them is withheld from the report and said
to staff instead, as a possibility (submission_report.sheet_cautions), in the brief and the
decisions log. Every other note is 2.13.2's.

Stdlib only.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import submission_report as sr  # noqa: E402
from test_submission_report import fixture, notes_by_id  # noqa: E402


def twice(names):
    return [x for x in names for _ in (0, 1)]


def record(names):
    r = fixture()
    r["submission_data"]["samples"] = [
        {"unique_id": f"SMP{100 + i}", "sample_name": f"sample {i}", "condition_name": c}
        for i, c in enumerate(names, 1)]
    return sr.sanitize(r)


def session(tmp, names, analysed):
    """A session whose conditions.csv holds `analysed` (by hand), each run assigned to its
    sample, with the sheet `names` attached."""
    s = os.path.join(tmp, "2026-10-20_conditions")
    os.makedirs(os.path.join(s, "input"))
    os.makedirs(os.path.join(s, "output"))
    rec = record(names)
    runs = {x["unique_id"]: f"10202026_DIA_{x['unique_id']}_1_{i}"
            for i, x in enumerate(rec["samples"], 1)}
    with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
        fh.write("# raw\n" + "".join(f"/data/OCT26/{r}.d\n" for r in runs.values()))
    with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        for x, g in zip(rec["samples"], analysed):
            w.writerow([runs[x["unique_id"]], g])
    with open(os.path.join(s, "input", "locate_decisions.json"), "w") as fh:
        json.dump({"assignments": [{"sample": u, "file": f"{r}.d"} for u, r in runs.items()]}, fh)
    sr.attach(s, rec)
    return s


# (the sheet, the analysis's groups by hand: the names are separate conditions)
GLUED = (
    (["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3, ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3, "ATG"),
    (twice(["WT3", "WT4", "WT6"]) + ["KO"] * 6, twice(["WT3", "WT4", "WT6"]) + ["KO"] * 6, "WT"),
    (twice(["Parental12", "Parental51"]) + ["Dox(+)"] * 4,
     twice(["Parental12", "Parental51"]) + ["Dox(+)"] * 4, "Parental"),
)


class GluedNumbersAreForStaff(unittest.TestCase):
    def test_the_groups_are_2_13_2s(self):
        for names, _analysed, base in GLUED:
            self.assertIn(base, set(sr.sheet_groups(record(names)).values()), names)

    def test_the_client_note_is_withheld_and_staff_are_told(self):
        for names, analysed, base in GLUED:
            with tempfile.TemporaryDirectory() as tmp:
                s = session(tmp, names, analysed)
                rec = sr.load(s)
                self.assertNotIn("conditions_differ", notes_by_id(rec, s), names)
                cautions = sr.sheet_cautions(rec, s)
                self.assertEqual([c["id"] for c in cautions], ["conditions_glued_names"], names)
                text = cautions[0]["text"]
                self.assertIn(f"the sheet's “{base}” is split across groups", text)
                self.assertIn(f"were read as replicates of “{base}”", text)
                self.assertIn("staff only, a possibility, never in the report", text)
                # not in the report's Submission section either
                self.assertNotIn("split across", sr.report_section(session=s)["md"])

    def test_a_note_not_resting_on_glued_numbers_is_2_13_2s(self):
        """Marked replicate numbers (WT_1): the client note stays, and no caution."""
        names = [f"{g}_{n}" for g in ("WT", "KO") for n in (1, 2, 3)]
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp, names, ["WT"] * 6)
            rec = sr.load(s)
            self.assertIn("group “WT” pools the sheet's KO, WT",
                          notes_by_id(rec, s)["conditions_differ"])
            self.assertEqual(sr.sheet_cautions(rec, s), [])

    def test_glued_and_other_parts_are_split_between_client_and_staff(self):
        names = ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 3 + ["KO"] * 3
        analysed = ["ATG5"] * 3 + ["ATG7"] * 3 + ["WT"] * 6        # WT pools the sheet's KO
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp, names, analysed)
            rec = sr.load(s)
            note = notes_by_id(rec, s)["conditions_differ"]
            self.assertIn("group “WT” pools the sheet's KO, WT", note)
            self.assertNotIn("ATG", note)
            self.assertIn("“ATG” is split across", sr.sheet_cautions(rec, s)[0]["text"])

    def test_the_caution_is_logged_once_and_in_the_brief(self):
        names, analysed, _base = GLUED[0]
        with tempfile.TemporaryDirectory() as tmp:
            s = session(tmp, names, analysed)
            for _ in range(2):                           # `notes` run twice: logged once
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "submission_report.py"),
                                    "notes", "--session", s], capture_output=True, text=True,
                                   timeout=120)
                self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertNotIn("conditions_differ", [n["id"] for n in out["notes"]])
            self.assertEqual([c["id"] for c in out["staff_cautions"]], ["conditions_glued_names"])
            with open(os.path.join(s, "logs", "decisions.md"), encoding="utf-8") as fh:
                self.assertEqual(fh.read().count("“ATG” is split across"), 1)
            brief = []
            import analysis_prompt
            analysis_prompt.submission_brief(brief.append, s)
            brief = "\n".join(brief)
            self.assertIn("For the Core staff member only -- never in the report", brief)
            self.assertIn("“ATG” is split across", brief)


if __name__ == "__main__":
    unittest.main()
