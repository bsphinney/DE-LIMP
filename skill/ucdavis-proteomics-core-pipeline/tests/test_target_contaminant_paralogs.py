#!/usr/bin/env python3
"""A peptide shared across a paralog family does not make the family a contaminant.

PROT_0756 v2 (2026-09-28): AUDIT.md's target_contaminants listed Ywhab/Ywhae/Ywhag/Ywhah/Ywhaq,
Sfn, Tuba8, Tubal3 and Eef1a2 -- abundant endogenous brain proteins -- as possible contamination.
The shared_peptides records for 1433Z_BOVIN / TBA1D_BOVIN / EF1A1_BOVIN carry target_accs = every
target sharing >= 1 peptide (the paralogs), and target_contaminants() flagged all of them. Only a
record's matched pair -- the target it is identical to, contained in, or (shared_peptides) shares
the most peptides with -- is the protein it cannot be told apart from.
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

import fetch_fasta as ff   # noqa: E402

# A 14-3-3 family built from tryptic blocks (the skill's digest: K*,R*, 1 missed cleavage,
# 7-30 aa). All three share one conserved peptide; the bovine stand-in is zeta with one residue
# changed, so the search cannot tell it from zeta -- and shares that one peptide with beta/eps.
C1 = "AEDLLAEYK"
ZETA = "G" + "WDKNELVQK" + C1 + "GSVTEQGAELSNEERNLLSVAYKNVVGARRSSWRVVSSLEQKTEGAEKKQQMAR"
BETA = "G" + "TMDKSELVQK" + C1 + "LAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVPYK"
EPS = "G" + "DDREDLVYQAK" + C1 + "LAEQAERYDEMVESMKKVAGMDVELTVEER"
TARGET = (f">sp|P63101|1433Z_MOUSE 14-3-3 zeta OS=Mus musculus OX=10090 GN=Ywhaz PE=1 SV=1\n{ZETA}\n"
          f">sp|Q9CQV8|1433B_MOUSE 14-3-3 beta OS=Mus musculus OX=10090 GN=Ywhab PE=1 SV=3\n{BETA}\n"
          f">sp|P62259|1433E_MOUSE 14-3-3 epsilon OS=Mus musculus OX=10090 GN=Ywhae PE=1 SV=1\n{EPS}\n")
CONT = (f">sp|Cont_P63103|1433Z_BOVIN 14-3-3 zeta OS=Bos taurus OX=9913 GN=YWHAZ PE=1 SV=1\n"
        f"{ZETA[:-3]}W{ZETA[-2:]}\n")
SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]


class MatchedPairOnly(unittest.TestCase):
    def test_shared_peptides_record_flags_its_matched_target_not_the_family(self):
        _kept, dropped, _enz = ff.drop_target_contaminants(CONT, TARGET)
        (rec,) = dropped
        self.assertEqual((rec["reason"], rec["target_acc"], rec["gene"]),
                         ("shared_peptides", "P63101", "Ywhaz"))
        # The family is still ON RECORD (information) ...
        self.assertEqual(rec["target_accs"], ["P63101", "Q9CQV8", "P62259"])
        # ... but only the matched pair is flagged.
        self.assertEqual(ff.matched_target_accs(rec), ["P63101"])
        tc = ff.target_contaminants({"contaminants_dropped_as_target": dropped})
        self.assertEqual((tc["accessions"], tc["genes"]), ({"P63101"}, {"YWHAZ"}))

    def test_identical_record_flags_every_identical_target(self):
        """Two targets with one sequence are both the contaminant: both stay flagged."""
        twin = TARGET + f">sp|P99999|1433Z2_MOUSE zeta copy OS=Mus musculus OX=10090 GN=Ywhaz2 PE=1\n{ZETA}\n"
        cont = f">sp|Cont_P63103|1433Z_BOVIN copy OS=Bos taurus GN=YWHAZ PE=1\n{ZETA}\n"
        _kept, (rec,), _enz = ff.drop_target_contaminants(cont, twin)
        self.assertEqual(rec["reason"], "identical")
        self.assertEqual(sorted(ff.matched_target_accs(rec)), ["P63101", "P99999"])


class AuditorsDoNotFlagParalogs(unittest.TestCase):
    """End to end: a real fetch builds the sidecar; the auditors read it."""

    @classmethod
    def setUpClass(cls):
        cls._td = tempfile.TemporaryDirectory()
        root = cls.root = cls._td.name
        tgt, cp = os.path.join(root, "mouse.fasta"), os.path.join(root, "cont.fasta")
        with open(tgt, "w") as fh:
            fh.write(TARGET)
        with open(cp, "w") as fh:
            fh.write(CONT)
        out = os.path.join(root, "search.fasta")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fetch_fasta.py"), "fetch",
                            "--path", tgt, "--contaminants", "universal", "--contaminants-path", cp,
                            "--out", out], capture_output=True, text=True, cwd=root)
        assert p.returncode == 0, p.stderr
        cls.meta = out + ".meta.json"
        de = cls.de = os.path.join(root, "de_results")
        os.makedirs(de)
        # run_de.R's column order: Protein.Group, Genes, Protein.Names, then samples. Zeta is
        # up in S6 so sample_quality's panel has something to see.
        with open(os.path.join(de, "Expression_Matrix.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Protein.Group", "Genes", "Protein.Names"] + SAMPLES)
            w.writerow(["P63101", "Ywhaz", "1433Z_MOUSE"] + [20, 20, 20, 20, 20, 24])
            w.writerow(["Q9CQV8", "Ywhab", "1433B_MOUSE"] + [19] * 6)
            w.writerow(["P62259", "Ywhae", "1433E_MOUSE"] + [21] * 6)
            w.writerow(["P04406", "Gapdh", "G3P_MOUSE"] + [22] * 6)

    @classmethod
    def tearDownClass(cls):
        cls._td.cleanup()

    def test_audit_flags_only_the_matched_protein(self):
        out = os.path.join(self.root, "AUDIT.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"), "--out", out,
                            "--de-dir", self.de, "--fasta-meta", self.meta],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(self.root, "AUDIT.json")) as fh:
            f = next(x for x in json.load(fh)["findings"] if x["check"] == "target_contaminants")
        self.assertEqual(f["detail"]["quantified"], ["Ywhaz"])
        for paralog in ("Ywhab", "Ywhae"):
            self.assertNotIn(paralog, f["message"])

    def test_sample_quality_panel_holds_only_the_matched_protein(self):
        out = os.path.join(self.root, "SAMPLE_QUALITY.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                            "--matrix", os.path.join(self.de, "Expression_Matrix.csv"),
                            "--fasta-meta", self.meta, "--out", out],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out.replace(".md", ".json")) as fh:
            panel = json.load(fh)["panels"]["CONTAMINANT_IDENTICAL"]
        self.assertEqual(panel["n_panel_proteins"], 1)
        self.assertNotIn("YWHAB", panel["matched_genes"])
        self.assertNotIn("Q9CQV8", panel["matched_genes"])


if __name__ == "__main__":
    unittest.main()
