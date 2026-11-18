#!/usr/bin/env python3
"""
2.15.1 review, finding 3: on Sage data the on/off list counted identifications while limma tested
the values Sage's LFQ measured by matching between runs, so a protein could sit in
OnOff_Proteins.csv ("no fold change") AND be significant in the DE table for the same contrast; a
protein with matched values in both groups was set to "not measured in either group" although
Expression_Matrix had its values.

Now one definition of "present in a run" (a value in the matrix limma tests) for Detection_Matrix,
the Evidence, the not-measured rule, the incomplete flag and the on/off list; matched values are
named beside it (DE: matched_between_runs_only_in; on/off: matched_between_runs_only).

2.15.1 wording correction: Sage's LFQ (lfq.rs, v0.14.7) MEASURES MS1 signal in every run at the
aligned retention time -- nothing is "transferred". A run where it found no signal has intensity 0:
such a row is never present and never matched (ZERO, IDZERO below).

The report follows the reviewer's synthetic case (review/tmp/e3/gen.py), three groups of three.
Needs pyarrow and R with limma/arrow/dplyr/tidyr/jsonlite.
"""
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)

from test_run_de_contaminants import r_has  # noqa: E402

sys.path.insert(0, SCRIPTS)
import sage_lfq_check as slc  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

NEEDS = ("limma", "arrow", "dplyr", "tidyr", "jsonlite")
RUNS = [(f"{g.lower()}{i}", g) for g in "ABC" for i in range(1, 4)]
COLS = ["Run", "Protein.Group", "Protein.Names", "Peptide.Id", "PG.MaxLFQ", "Q.Value",
        "Lib.Q.Value", "Lib.PG.Q.Value", "Identified.In.Run"]
META = {"delimp.quantity.level": "peptide", "delimp.quantity.label": "Sage",
        "delimp.quantity.value": "Sage synthetic", "delimp.quantity.identification": "synthetic",
        "delimp.quantity.q_columns": "placeholder", "delimp.quantity.citation": "synthetic",
        "delimp.quantity.engine_normalised": "false",
        "delimp.quantity.matched_between_runs": slc.MATCHED_BETWEEN_RUNS}


def write_report(out):
    rnd = random.Random(3)
    load = {r: rnd.uniform(-0.6, 0.6) for r, _ in RUNS}
    rows = {c: [] for c in COLS}

    def add(r, acc, pep, lv, ident=1):
        for c, x in zip(COLS, (r, acc, "N" + acc, pep, None if lv is None else 2 ** lv,
                               0.0, 0.0, 0.0, ident)):
            rows[c].append(x)
    for i in range(300):
        acc, base = f"P{i:05d}", rnd.uniform(18, 24)
        for k in range(3):
            off = rnd.gauss(0, 1)
            for r, g in RUNS:
                add(r, acc, f"PEP{i}K{k}", base + off + load[r] + rnd.gauss(0, 0.15))
    for r, g in RUNS:
        n = load[r] + rnd.gauss(0, 0.1)
        if g == "A":            # identified in A; B values all matched between runs, 3 log2 lower
            add(r, "XFERLOW", "PEPY", 20 + n, 1)
        if g == "B":
            add(r, "XFERLOW", "PEPY", 17 + n, 0)
        if g in "AB":           # matched everywhere in A and B, identified only in C
            add(r, "NOID", "PEPZ", 20 + n + (2 if g == "B" else 0), 0)
        if g == "C":
            add(r, "NOID", "PEPZ", 20 + n, 1)
        if g == "C":            # values only in C, every one matched between runs
            add(r, "TONLY", "PEPT", 21 + n, 0)
        # ZERO: identified in A; in B Sage found no MS1 signal (0, one null): absent in B, never
        # matched, so on/off A_only in B - A, no fold change
        if g == "A":
            add(r, "ZERO", "PEPZE", 21 + n, 1)
        if g == "B":
            add(r, "ZERO", "PEPZE", None if r == "b3" else float("-inf"), 0)
        # IDZERO: in B its MS2-identified peptide has no signal (0) and only a matched peptide has
        # signal, so B's values rest on matching alone; identified with signal in A
        if g in "AB":
            add(r, "IDZERO", "PEPIA", 21 + n if g == "A" else float("-inf"), 1)
            add(r, "IDZERO", "PEPIB", 20 + n, 1 if g == "A" else 0)
    pq.write_table(pa.table(rows).replace_schema_metadata(META), os.path.join(out, "report.parquet"))
    with open(os.path.join(out, "conditions.csv"), "w") as fh:
        fh.write("File.Name,Group\n" + "".join(f"{r},{g}\n" for r, g in RUNS))


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class OnOffAndDeAgree(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        write_report(cls.d)
        cls.out = os.path.join(cls.d, "de")
        cls.p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"),
                                "--input", "report.parquet", "--metadata", "conditions.csv",
                                "--method", "maxlfq", "--outdir", cls.out],
                               capture_output=True, text=True, cwd=cls.d, timeout=600)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def setUp(self):
        self.assertEqual(self.p.returncode, 0, self.p.stderr[-3000:])

    def de(self, name):
        with open(os.path.join(self.out, f"DE_maxlfq_{name}.csv"), newline="") as fh:
            return {r["Protein.Group"]: r for r in csv.DictReader(fh)}

    def onoff(self):
        with open(os.path.join(self.out, "OnOff_Proteins.csv"), newline="") as fh:
            return list(csv.DictReader(fh))

    @staticmethod
    def has_p(row):
        return row is not None and row["P.Value"] not in ("", "NA")

    def test_no_protein_is_on_off_and_tested_in_the_same_contrast(self):
        on = self.onoff()
        self.assertTrue(on)
        for name in ("B.A", "C.A"):
            tab = self.de(name)
            label = name.replace(".", " - ")
            for r in on:
                if r["Contrast"] == label:
                    self.assertFalse(self.has_p(tab.get(r["Protein.Group"])),
                                     f"{r['Protein.Group']} on/off and tested in {label}")

    def test_matched_values_are_tested_and_named(self):
        tab = self.de("B.A")
        self.assertTrue(self.has_p(tab["XFERLOW"]))
        self.assertEqual(tab["XFERLOW"]["matched_between_runs_only_in"], "B")
        self.assertNotIn("XFERLOW", {r["Protein.Group"] for r in self.onoff()
                                     if r["Contrast"] == "B - A"})
        # values in both groups (all transferred): tested, never "not measured in either group"
        self.assertTrue(self.has_p(tab["NOID"]))
        self.assertEqual(tab["NOID"]["Evidence"], "measured in both")
        self.assertEqual(set(tab["NOID"]["matched_between_runs_only_in"].split(";")), {"A", "B"})
        self.assertEqual(tab["P00001"]["matched_between_runs_only_in"], "")

    def test_an_on_off_protein_present_only_by_matching_is_flagged(self):
        rows = {(r["Protein.Group"], r["Contrast"]): r for r in self.onoff()}
        self.assertEqual(rows[("TONLY", "C - A")]["Direction"], "C_only")
        self.assertEqual(rows[("TONLY", "C - A")]["matched_between_runs_only"], "TRUE")
        self.assertEqual(rows[("XFERLOW", "C - A")]["matched_between_runs_only"], "FALSE")
        with open(os.path.join(self.out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertIn("counts as present", prov["onoff"]["presence"])
        self.assertIn(slc.MATCHED_BETWEEN_RUNS, prov["onoff"]["presence"])
        self.assertNotIn("transferred", prov["onoff"]["presence"] + prov["onoff"]["rule"])

    def test_no_ms1_signal_is_absent_never_matched(self):
        """A row with intensity 0 / null (Sage found no MS1 signal) is a missing value: absent in
        Detection_Matrix, not counted as matched, and the on/off list treats the group as
        absent."""
        with open(os.path.join(self.out, "Detection_Matrix.csv"), newline="") as fh:
            det = {r[next(iter(r))]: r for r in csv.DictReader(fh)}
        self.assertEqual({det["ZERO"][r] for r, g in RUNS if g == "B"}, {"0"})
        self.assertEqual({det["ZERO"][r] for r, g in RUNS if g == "A"}, {"1"})
        on = {(r["Protein.Group"], r["Contrast"]): r for r in self.onoff()}
        self.assertEqual(on[("ZERO", "B - A")]["Direction"], "A_only")
        self.assertEqual(on[("ZERO", "B - A")]["matched_between_runs_only"], "FALSE")
        tab = self.de("B.A")
        self.assertFalse(self.has_p(tab.get("ZERO")))
        if tab.get("ZERO") is not None:
            self.assertNotIn("B", tab["ZERO"]["matched_between_runs_only_in"].split(";"))
        # an identified peptide with no signal does not make the cell "identified": IDZERO's B
        # values rest on its matched peptide alone
        self.assertEqual(tab["IDZERO"]["matched_between_runs_only_in"], "B")
        # matched_cells counts cells with a value and no identified peptide with signal: B's
        # XFERLOW, NOID, IDZERO; A's NOID; C's TONLY -- 3 runs each. ZERO's B runs are not cells.
        with open(os.path.join(self.out, "de_provenance.json")) as fh:
            dm = json.load(fh)["detection_matrix"]
        self.assertEqual(dm["matched_cells"], 5 * 3)
        self.assertEqual(dm["matched_definition"], slc.MATCHED_BETWEEN_RUNS)
        self.assertIn(slc.MATCHED_BETWEEN_RUNS, dm["values"])

    def test_the_methods_say_on_off_has_no_fold_change_and_mean_it(self):
        with open(os.path.join(self.out, "methods.txt")) as fh:
            txt = fh.read()
        self.assertNotIn("not identified", txt)


if __name__ == "__main__":
    unittest.main()
