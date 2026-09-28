#!/usr/bin/env python3
"""collect_conditions.py keeps a subject / animal column under its own name.

PROT_0756 (Dickson lab): five IPs from each of six mice. --map filed any extra sample-sheet
column under Covariate1/2 -- a fixed effect in run_de.R -- so a Mouse column became
Covariate1; nested in the Old/Young groups that made the design rank-deficient, and the
pairing could not be passed to run_de.R --block. Now a column named for the unit samples
come from (Mouse, Animal, Subject, Patient, Donor ...) is written as that column, and the
--map report names it and says whether --block applies.

The header alone is not trusted (review C1): a 'Subject' of M/F, a 'Patient' of Yes/No or a
'Donor' that is the Group relabelled would silently drop that factor from the model as a
block. Their VALUES must look like subjects -- >= 3 of them, not the groups relabelled, and
recurring across groups where they span them -- or the column stays a covariate and --map
asks (subject_ambiguous) until --subject-column confirms it.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "collect_conditions.py")
sys.path.insert(0, os.path.dirname(SCRIPT))

import collect_conditions as cc   # noqa: E402

RUNS = [f"08132026_DIA-LRS-{n}_S3" for n in range(96, 106)]


def run_map(tmp, sheet_rows, header, runs=RUNS, extra=(), **kw):
    sheet = os.path.join(tmp, "sheet.csv")
    with open(sheet, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(sheet_rows)
    out = os.path.join(tmp, "conditions.csv")
    args = ["python3", SCRIPT, "--map", out, "--runs", ",".join(runs), *extra]
    args += ["--mapping-json", kw["mapping_json"]] if "mapping_json" in kw else ["--from-file", sheet]
    p = subprocess.run(args, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    with open(out, newline="") as fh:
        return json.loads(p.stdout), list(csv.DictReader(fh))


def ip_sheet(subject_header="Mouse", unique=False):
    # LRS96-105: two baits x five mice, one IP per bait per mouse
    rows = []
    for i, n in enumerate(range(96, 106)):
        bait = "JPH3" if i < 5 else "IgG"
        mouse = f"Mouse {i + 1}" if unique else f"Mouse {i % 5 + 1}"
        rows.append([f"LRS-{n}_", bait, "plate1", "F" if i % 2 else "M", mouse])
    return rows, ["sample", "condition", "batch", "Sex", subject_header]


class SubjectColumnKeptUnderItsOwnName(unittest.TestCase):
    def test_mouse_column_is_kept_not_a_covariate(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet()
            rep, out = run_map(tmp, rows, header)
        self.assertEqual(list(out[0]), ["File.Name", "Group", "Batch", "Covariate1", "Mouse"])
        self.assertEqual([r["Mouse"] for r in out], [f"Mouse {i % 5 + 1}" for i in range(10)])
        self.assertEqual({r["Covariate1"] for r in out}, {"M", "F"})   # Sex, not the mouse
        self.assertNotIn("Covariate2", out[0])
        self.assertEqual(rep["block_column"], "Mouse")
        self.assertTrue(rep["block_suggested"])
        self.assertEqual(rep["subjects"], {f"Mouse {k}": 2 for k in range(1, 6)})
        self.assertIn("--block Mouse", rep["guidance"])
        self.assertFalse(rep["needs_confirmation"])

    def test_header_variants(self):
        for header, col in (("Animal ID", "Animal_ID"), ("patient", "patient"), ("Donor#", "Donor"),
                            ("subject_no", "subject_no"), ("Mice", "Mice")):
            with self.subTest(header=header), tempfile.TemporaryDirectory() as tmp:
                rows, hdr = ip_sheet(header)
                rep, out = run_map(tmp, rows, hdr)
                self.assertEqual(rep["block_column"], col)
                self.assertIn(col, out[0])
                self.assertNotIn(header if header != col else "Covariate2", out[0])

    def test_not_a_subject_column_stays_a_covariate(self):
        for header in ("Sex", "Tissue", "Treatment_time", "Plate"):
            with self.subTest(header=header):
                self.assertFalse(cc.subject_header(header))
        for header in ("Mouse", "mouse_id", "Animal Number", "Patient #", "Pair", "Block"):
            with self.subTest(header=header):
                self.assertTrue(cc.subject_header(header))

    def test_one_run_per_subject_is_an_id_not_a_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet(unique=True)
            rep, out = run_map(tmp, rows, header)
        self.assertIn("Mouse", out[0])            # still carried, under its name
        self.assertFalse(rep["block_suggested"])  # but every mouse has one run: no --block
        self.assertNotIn("--block", rep["guidance"])
        self.assertEqual(rep["ambiguities"]["single_run_subjects"], [])

    def test_some_single_run_subjects_need_confirmation(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet()
            rows[-1][-1] = "Mouse 9"   # Mouse 5 and Mouse 9 now hold one run each
            rep, _ = run_map(tmp, rows, header)
        self.assertFalse(rep["block_suggested"])
        self.assertEqual(rep["ambiguities"]["single_run_subjects"], ["Mouse 5", "Mouse 9"])
        self.assertTrue(rep["needs_confirmation"])

    def test_agent_built_mapping_can_carry_subjects(self):
        groups = {"JPH3": [f"LRS-{n}_" for n in range(96, 101)],
                  "IgG": [f"LRS-{n}_" for n in range(101, 106)]}
        subjects = {f"LRS-{n}_": f"M{(n - 96) % 5 + 1}" for n in range(96, 106)}
        with tempfile.TemporaryDirectory() as tmp:
            rep, out = run_map(tmp, [], ["x"], mapping_json=json.dumps(
                {"groups": groups, "subjects": subjects, "subject_column": "Mouse"}))
        self.assertEqual(list(out[0]), ["File.Name", "Group", "Mouse"])
        self.assertEqual(rep["block_column"], "Mouse")
        self.assertTrue(rep["block_suggested"])

    # -- review C1: a subject header over values that are not subjects ----------------
    def assert_ambiguous(self, rep, out, header, reason, slot):
        self.assertIsNone(rep["block_column"])
        self.assertFalse(rep["block_suggested"])
        self.assertTrue(rep["needs_confirmation"])
        amb = rep["ambiguities"]["subject_ambiguous"]
        self.assertEqual(amb["column"], header)
        self.assertTrue(any(reason in r for r in amb["reasons"]), amb["reasons"])
        self.assertEqual(amb["written_as"], slot)
        self.assertIn(f"--subject-column '{header}'", amb["to_confirm"])
        self.assertNotIn(header, out[0])                      # not written as a block column
        self.assertIn(slot, out[0])                           # kept, as a covariate
        self.assertNotIn("--block", rep["guidance"])

    def test_subject_that_is_sex_is_not_a_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet("Subject")
            for i, r in enumerate(rows):
                r[-1] = "M" if i % 2 == 0 else "F"
            rep, out = run_map(tmp, rows, header)
        self.assert_ambiguous(rep, out, "Subject", "only 2 distinct value(s) (F, M)", "Covariate2")
        self.assertEqual({r["Covariate2"] for r in out}, {"M", "F"})

    def test_patient_yes_no_is_not_a_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet("Patient")
            for i, r in enumerate(rows):
                r[-1] = "Yes" if i < 3 or i > 7 else "No"
            rep, out = run_map(tmp, rows, header)
        self.assert_ambiguous(rep, out, "Patient", "only 2 distinct value(s)", "Covariate2")

    def test_donor_that_is_the_group_is_not_a_block(self):
        runs = [f"R{i:02d}" for i in range(1, 10)]
        rows = [[f"R{i:02d}", g, f"D_{g}"] for i, g in zip(range(1, 10), ["A", "B", "C"] * 3)]
        with tempfile.TemporaryDirectory() as tmp:
            rep, out = run_map(tmp, rows, ["sample", "group", "Donor"], runs=runs)
        self.assert_ambiguous(rep, out, "Donor", "the groups relabelled", "Covariate1")

    def test_ambiguous_subject_with_no_free_slot_says_it_was_not_written(self):
        # two other extra columns fill Covariate1/2: the ambiguous subject cannot be kept as a
        # covariate -- the report must say so, not drop it silently
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet("Subject")
            for i, r in enumerate(rows):
                r[-1] = "M" if i % 2 == 0 else "F"
                r.insert(4, f"T{i % 3}")
            header = header[:4] + ["Tissue"] + header[4:]
            rep, out = run_map(tmp, rows, header)
        amb = rep["ambiguities"]["subject_ambiguous"]
        self.assertIsNone(amb["written_as"])
        self.assertFalse(amb["written"])
        self.assertIn("it was NOT written, because both covariate slots are taken "
                      "(Covariate1 = Sex, Covariate2 = Tissue)", amb["to_confirm"])
        self.assertIn("--subject-column 'Subject' to keep it as the block", amb["to_confirm"])
        self.assertNotIn("Subject", out[0])
        self.assertEqual(rep["covariate_columns"], {"Batch": "batch", "Covariate1": "Sex",
                                                    "Covariate2": "Tissue"})
        self.assertTrue(rep["needs_confirmation"])

    def test_more_than_two_extra_columns_are_reported_not_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet()
            rows = [r[:4] + [f"T{i % 3}", f"L{i % 2}"] for i, r in enumerate(rows)]
            rep, out = run_map(tmp, rows, header[:4] + ["Tissue", "Lane"])
        self.assertEqual(list(out[0]), ["File.Name", "Group", "Batch", "Covariate1", "Covariate2"])
        nw = rep["ambiguities"]["columns_not_written"]
        self.assertEqual(nw["columns"], ["Lane"])
        self.assertIn("NOT written", nw["why"])
        self.assertTrue(rep["needs_confirmation"])

    def test_confirmed_subject_column_is_kept(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet("Subject")
            for i, r in enumerate(rows):
                r[-1] = "M" if i % 2 == 0 else "F"
            rep, out = run_map(tmp, rows, header, extra=["--subject-column", "Subject"])
        self.assertEqual(rep["block_column"], "Subject")
        self.assertTrue(rep["subject_confirmed"])
        self.assertIn("Subject", out[0])
        self.assertNotIn("subject_ambiguous", rep["ambiguities"])

    def test_subject_column_none_turns_detection_off(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet()
            rep, out = run_map(tmp, rows, header, extra=["--subject-column", "none"])
        self.assertIsNone(rep["block_column"])
        self.assertEqual(list(out[0]), ["File.Name", "Group", "Batch", "Covariate1", "Covariate2"])

    def test_assessment_directly(self):
        grp = {f"r{i}": g for i, g in enumerate("AABBCC")}
        self.assertEqual(cc.subject_assessment({f"r{i}": s for i, s in enumerate("123123")}, grp), [])
        # one subject per group and one group per subject: the groups relabelled
        relab = {"r0": "1", "r1": "1", "r2": "2", "r3": "2", "r4": "3", "r5": "3"}
        self.assertTrue(any("relabelled" in x for x in cc.subject_assessment(relab, grp)))
        # nested replicates (two subjects per group, two runs each) are a valid block
        grp2 = {f"r{i}": g for i, g in enumerate("AAAABBBB")}
        reps = {f"r{i}": s for i, s in enumerate("11223344")}
        self.assertEqual(cc.subject_assessment(reps, grp2), [])
        # most values confined to one group, one crossing: not a subject design
        few = {f"r{i}": s for i, s in enumerate("11234156")}
        self.assertTrue(any("recur across groups" in x for x in cc.subject_assessment(few, grp2)))

    def test_no_subject_column_no_change(self):
        with tempfile.TemporaryDirectory() as tmp:
            rows, header = ip_sheet()
            rep, out = run_map(tmp, [r[:4] for r in rows], header[:4])
        self.assertEqual(list(out[0]), ["File.Name", "Group", "Batch", "Covariate1"])
        self.assertIsNone(rep["block_column"])
        self.assertFalse(rep["block_suggested"])
        self.assertNotIn("single_run_subjects", rep["ambiguities"])


if __name__ == "__main__":
    unittest.main()
