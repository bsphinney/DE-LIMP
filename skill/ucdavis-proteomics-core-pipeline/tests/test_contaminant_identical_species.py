#!/usr/bin/env python3
"""Contaminant-identical proteins carry their own species; gene symbols match only within one
(2.14.2 triage I21).

PXD063416 (human + yeast, 2026-10-27): the Universal contaminant set's yeast ADH1 and TPM2 are
identical to the yeast proteins of the searched database, so fetch_fasta.py dropped them and the
auditors flagged the yeast proteins -- calling them "Homo sapiens" proteins (the database's primary
organism) and matching by gene symbol, so human TPM2 was flagged too.

What these tests pin:
  * every dropped record carries target_species (the entry name's species code) and
    target_organism (OS=), both None when the header has neither
  * audit_results.species_code() is THE reader of a species code; sample_quality uses the same
  * the rule (fetch_fasta.identical_target_matcher, both auditors): accession first; gene symbol
    only within the record's own species; no species (NCBI, an old sidecar) -> accession only,
    and an old sidecar's finding says so
  * each hit is named with its own organism, not the database's primary one
stdlib only, no network.
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

import audit_results as ar      # noqa: E402
import fetch_fasta as ff        # noqa: E402
import sample_quality as sq     # noqa: E402

YEAST_OS = "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"
# synthetic sequences from tryptic blocks (the skill's digest), no two sharing a peptide
ADH1_Y = "M" + "SIPETQK" + "GVIFYESHGK" + "LEYKDIPVPK" + "ANELLINVK" + "YSGVCHTDLHAWHGDWPLPVK"
TPM2_Y = "M" + "EAIFDTLK" + "SVQEELTSAEGR" + "VSSLEQQLTETSK" + "NIEDLEFATDR" + "LQEELTTVSNELESK"
TPM2_H = "M" + "DAIKQMNPLK" + "LDKENAIDRAEQAEADK" + "WWGPCQSLEDHQGK" + "YTVHNIQDGMEK" + "FPYSAGVR"
G3P_H = "M" + "GKVGVNGFGR" + "IGRLVTRAAFNSGK" + "VDIVAINDPFIDLNYMVYMFQYDSTHGK" + "FHGTVK"
TARGET = (f">sp|P07951|TPM2_HUMAN Tropomyosin beta chain OS=Homo sapiens OX=9606 GN=TPM2 PE=1 SV=1\n"
          f"{TPM2_H}\n"
          f">sp|P04406|G3P_HUMAN Glyceraldehyde-3-phosphate dehydrogenase OS=Homo sapiens OX=9606 "
          f"GN=GAPDH PE=1 SV=3\n{G3P_H}\n"
          f">sp|P00330|ADH1_YEAST Alcohol dehydrogenase 1 OS={YEAST_OS} OX=559292 GN=ADH1 PE=1 "
          f"SV=5\n{ADH1_Y}\n"
          f">sp|P40414|TPM2_YEAST Tropomyosin-2 OS={YEAST_OS} OX=559292 GN=TPM2 PE=1 SV=1\n"
          f"{TPM2_Y}\n")
CONT = (f">sp|Cont_P00330|ADH1_YEAST Alcohol dehydrogenase 1 OS={YEAST_OS} GN=ADH1\n{ADH1_Y}\n"
        f">sp|Cont_P40414|TPM2_YEAST Tropomyosin-2 OS={YEAST_OS} GN=TPM2\n{TPM2_Y}\n")
SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]


def dropped_records():
    _kept, dropped, _enz = ff.drop_target_contaminants(CONT, TARGET)
    return dropped


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


class SidecarRecords(unittest.TestCase):

    def test_each_dropped_record_names_its_own_species(self):
        by = {r["target_acc"]: r for r in dropped_records()}
        self.assertEqual(set(by), {"P00330", "P40414"})
        for r in by.values():
            self.assertEqual(r["target_species"], "YEAST")
            self.assertEqual(r["target_organism"], YEAST_OS)

    def test_an_ncbi_header_has_neither(self):
        seq = "M" + "AEDLLAEYK" + "GSVTEQGAELSNEER"
        _kept, (r,), _enz = ff.drop_target_contaminants(
            f">Cont_XP_1.1 keratin\n{seq}\n", f">XP_000001.1 keratin [Homo sapiens]\n{seq}\n")
        self.assertIsNone(r["target_species"])
        self.assertIsNone(r["target_organism"])


class SpeciesCode(unittest.TestCase):

    def test_reads_the_entry_names_suffix(self):
        self.assertEqual(ar.species_code("ADH1_YEAST"), "YEAST")
        self.assertEqual(ar.species_code("A0A024R161_HUMAN"), "HUMAN")
        self.assertEqual(ar.species_code("K2C1_HUMAN;K2C1B_HUMAN"), "HUMAN")

    def test_none_without_one_or_when_they_disagree(self):
        for text in (None, "", "P00330", "Cont_P00330", "XP_000001.1", "TPM2_HUMAN;TPM2_YEAST"):
            self.assertIsNone(ar.species_code(text), text)

    def test_sample_quality_reads_it_with_the_same_function(self):
        self.assertIs(sq.species_code, ar.species_code)
        self.assertIs(sq.identical_target_matcher, ar.identical_target_matcher)


class Auditors(unittest.TestCase):

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.de = os.path.join(self.root, "de_results")
        os.makedirs(self.de)

    def tearDown(self):
        self._td.cleanup()

    def meta(self, old=False, no_entry=False):
        dropped = dropped_records()
        if old:                                        # a sidecar written before 2.14.2
            for r in dropped:
                r.pop("target_species")
                r.pop("target_organism")
                if no_entry:                           # ... and before 2.5.1: no entry name
                    r.pop("target_entry")
        p = os.path.join(self.root, "search.fasta.meta.json")
        with open(p, "w") as fh:
            json.dump({"organism": "Homo sapiens", "taxid": 9606,
                       "contaminants_dropped_as_target": dropped,
                       "contaminants_identical_to_target_kept": []}, fh)
        return p

    def write_tables(self, matrix_rows, de_rows=()):
        write_csv(os.path.join(self.de, "Expression_Matrix.csv"),
                  ["Protein.Group", "Genes", "Protein.Names"] + SAMPLES, matrix_rows)
        write_csv(os.path.join(self.de, "DE_limpa_B.vs.A.csv"),
                  ["Protein.Group", "Genes", "logFC", "adj.P.Val"], list(de_rows))

    def audit(self, matrix_rows, de_rows=(), old=False, no_entry=False):
        self.write_tables(matrix_rows, de_rows)
        out = os.path.join(self.root, "AUDIT.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                            "--out", out, "--de-dir", self.de, "--fasta-meta",
                            self.meta(old, no_entry)],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(self.root, "AUDIT.json")) as fh:
            return next(f for f in json.load(fh)["findings"] if f["check"] == "target_contaminants")

    def test_yeast_hit_is_labelled_yeast_not_the_primary_organism(self):
        f = self.audit([["P00330", "ADH1", "ADH1_YEAST"] + [20] * 6])
        self.assertEqual(f["status"], "WARN")
        self.assertIn(YEAST_OS, f["message"])
        self.assertNotIn("Homo sapiens", f["message"])
        self.assertEqual(f["detail"]["quantified"], ["ADH1"])

    def test_human_symbol_twin_is_not_flagged(self):
        f = self.audit([["P07951", "TPM2", "TPM2_HUMAN"] + [20] * 6,
                        ["P04406", "GAPDH", "G3P_HUMAN"] + [18] * 6],
                       de_rows=[["P07951", "TPM2", "2.0", "0.001"]])
        self.assertEqual(f["status"], "PASS")
        self.assertNotIn("TPM2", f["message"])

    def test_a_same_species_row_is_matched_by_symbol(self):
        # another accession of the yeast protein (an isoform) is still the same yeast TPM2
        f = self.audit([["P40414-2", "TPM2", "TPM2_YEAST"] + [20] * 6],
                       de_rows=[["P40414-2", "TPM2", "2.0", "0.001"]])
        self.assertEqual(f["detail"]["quantified"], ["TPM2"])
        self.assertEqual(f["detail"]["significant"], {"DE_limpa_B.vs.A.csv": ["TPM2"]})

    def test_mixed_species_hits_each_carry_their_own(self):
        rec = ff.target_contaminants({"organism": "Homo sapiens", "contaminants_dropped_as_target": [
            {"cont_acc": "Cont_P04264", "target_acc": "P04264", "target_entry": "K2C1_HUMAN",
             "gene": "KRT1", "reason": "identical", "target_accs": ["P04264"],
             "target_species": "HUMAN", "target_organism": "Homo sapiens"}] + dropped_records()})
        p = os.path.join(self.root, "search.fasta.meta.json")
        self.write_tables([["P00330", "ADH1", "ADH1_YEAST"] + [20] * 6,
                           ["P04264", "KRT1", "K2C1_HUMAN"] + [20] * 6])
        with open(p, "w") as fh:
            json.dump({"organism": "Homo sapiens",
                       "contaminants_dropped_as_target": rec["dropped"]}, fh)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                            "--out", os.path.join(self.root, "AUDIT.md"), "--de-dir", self.de,
                            "--fasta-meta", p], capture_output=True, text=True, cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(self.root, "AUDIT.json")) as fh:
            f = next(x for x in json.load(fh)["findings"] if x["check"] == "target_contaminants")
        self.assertIn(f"{YEAST_OS}: ADH1", f["message"])
        self.assertIn("Homo sapiens: KRT1", f["message"])

    def test_old_sidecar_takes_the_species_from_its_entry_name(self):
        # every sidecar since 2.5.1 records target_entry (ADH1_YEAST): PXD063416's is one
        f = self.audit([["P00330", "ADH1", "ADH1_YEAST"] + [20] * 6,
                        ["P40414-2", "TPM2", "TPM2_YEAST"] + [20] * 6,
                        ["P07951", "TPM2", "TPM2_HUMAN"] + [20] * 6], old=True)
        self.assertEqual(f["detail"]["quantified"], ["ADH1", "TPM2"])   # yeast TPM2 by symbol
        self.assertEqual(f["detail"]["organism_of"], {"ADH1": "YEAST", "TPM2": "YEAST"})
        self.assertIn("2 quantified YEAST protein group(s)", f["message"])
        self.assertNotIn("Homo sapiens", f["message"])
        self.assertNotIn("accession only", f["message"])

    def test_sidecar_without_species_or_entry_name_matches_by_accession_only(self):
        f = self.audit([["P00330", "ADH1", "ADH1_YEAST"] + [20] * 6,
                        ["P40414-2", "TPM2", "TPM2_YEAST"] + [20] * 6,
                        ["P07951", "TPM2", "TPM2_HUMAN"] + [20] * 6], old=True, no_entry=True)
        self.assertEqual(f["detail"]["quantified"], ["ADH1"])      # the accession match only
        self.assertIn("accession only", f["message"])

    def test_sample_quality_panel_uses_the_same_rule(self):
        # human TPM2 is elevated in S6: it must not make the contaminant-identical panel move
        self.write_tables([["P00330", "ADH1", "ADH1_YEAST"] + [20] * 6,
                           ["P07951", "TPM2", "TPM2_HUMAN"] + [10, 10, 10, 10, 10, 16],
                           ["P04406", "GAPDH", "G3P_HUMAN"] + [18] * 6])
        out = os.path.join(self.root, "SAMPLE_QUALITY.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                            "--matrix", os.path.join(self.de, "Expression_Matrix.csv"),
                            "--fasta-meta", self.meta(), "--out", out],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out.replace(".md", ".json")) as fh:
            panel = json.load(fh)["panels"]["CONTAMINANT_IDENTICAL"]
        self.assertEqual(panel["n_panel_proteins"], 1)
        self.assertNotIn("TPM2", panel["matched_genes"])
        self.assertIn("ADH1", panel["matched_genes"])
        self.assertEqual(panel["elevated_samples"], [])
        with open(out) as fh:
            md = fh.read()
        self.assertIn(f"{YEAST_OS} proteins whose sequence is also a common-contaminant", md)
        self.assertNotIn("Homo sapiens proteins whose sequence", md)
        self.assertNotIn("matched_by", panel)

    def sample_quality_old(self, no_entry):
        self.write_tables([["P00330", "ADH1", "ADH1_YEAST"] + [20] * 6,
                           ["P40414-2", "TPM2", "TPM2_YEAST"] + [20] * 6,
                           ["P07951", "TPM2", "TPM2_HUMAN"] + [20] * 6])
        out = os.path.join(self.root, "SAMPLE_QUALITY.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                            "--matrix", os.path.join(self.de, "Expression_Matrix.csv"),
                            "--fasta-meta", self.meta(old=True, no_entry=no_entry), "--out", out],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out.replace(".md", ".json")) as fh:
            panel = json.load(fh)["panels"]["CONTAMINANT_IDENTICAL"]
        with open(out) as fh:
            return panel, fh.read()

    def test_sample_quality_takes_an_old_sidecars_species_from_its_entry_name(self):
        panel, md = self.sample_quality_old(no_entry=False)
        self.assertEqual(panel["n_panel_proteins"], 2)        # yeast ADH1 and TPM2, not human
        self.assertNotIn("matched_by", panel)
        self.assertIn("_YEAST proteins whose sequence", md)
        self.assertNotIn("Homo sapiens proteins whose sequence", md)

    def test_sample_quality_says_a_sidecar_without_names_matches_by_accession_only(self):
        panel, md = self.sample_quality_old(no_entry=True)
        self.assertEqual(panel["n_panel_proteins"], 1)                 # ADH1 by accession
        self.assertIn("accession only", panel["matched_by"])
        self.assertIn("Matched by accession only", md)


if __name__ == "__main__":
    unittest.main()
