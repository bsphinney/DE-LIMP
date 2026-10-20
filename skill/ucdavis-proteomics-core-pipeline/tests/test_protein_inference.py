#!/usr/bin/env python3
"""
DIA-NN's protein inference level is recorded and stated, never left to be guessed
(diann_parallel.protein_inference: search_provenance.json `protein_inference`, the Methods).

A staff report (PROT_0000, skill 2.11.2): a transcriptome FASTA was rewritten with UniProt-style
GN= tags so that ORFs of one gene would group; no step of the 5-step chain passed --pg-level,
nothing recorded which level ran, and the report overstated what the tags did to the protein
groups. Measured on HIVE with DIA-NN 2.7.0 (2026-10-16, a 6,752-entry
subset of that FASTA and one of its runs):
  * no --pg-level = --pg-level 2 (genes): identical Protein.Group on all 12,190 precursors, and
    the log counts "genes"; 0 (isoform IDs) and 1 (protein names) regroup 242 of them;
  * at every level a gene whose ORFs have peptides of their own keeps one protein group per ORF
    (35 genes, the tagged one included): the gene level is the gene-group matrix (Genes);
  * step 5 re-run at 0, 1 and 2 on the same .quant files and empirical library: identical groups.
    A level other than the default is set from step 1 (a new search), or it does nothing.

Stdlib only.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import diann_parallel as dp  # noqa: E402
import make_methods  # noqa: E402

BASE = "--qvalue 0.01\n--cut K*,R*\n--missed-cleavages 1\n--met-excision\n--unimod4\n"


def cfg(d, extra=""):
    p = os.path.join(d, "params.cfg")
    with open(p, "w") as fh:
        fh.write(BASE + extra)
    return p


def groups(path):
    return dp.cfg_groups(dp.cfg_tokens(path))


class Record(unittest.TestCase):
    def test_no_flag_is_dianns_default_genes_as_measured(self):
        with tempfile.TemporaryDirectory() as d:
            r = dp.protein_inference(groups(cfg(d)))
        self.assertTrue(r["inference"])
        self.assertEqual((r["pg_level"], r["basis"]), (2, "genes"))
        self.assertIsNone(r["passed"])
        self.assertIn("default", r["source"])
        self.assertIn("2.7.0", r["default_measured"])

    def test_a_level_in_the_cfg_is_what_ran(self):
        for n, basis in ((0, "isoform IDs"), (1, "protein names"), (2, "genes")):
            with self.subTest(level=n), tempfile.TemporaryDirectory() as d:
                r = dp.protein_inference(groups(cfg(d, f"--pg-level {n}\n")))
                self.assertEqual((r["pg_level"], r["basis"], r["passed"]),
                                 (n, basis, f"--pg-level {n}"))

    def test_inference_switched_off(self):
        with tempfile.TemporaryDirectory() as d:
            r = dp.protein_inference(groups(cfg(d, "--no-prot-inf\n")))
        self.assertFalse(r["inference"])
        self.assertEqual(r["passed"], "--no-prot-inf")

    def test_an_unknown_level_is_said_not_guessed(self):
        with tempfile.TemporaryDirectory() as d:
            r = dp.protein_inference(groups(cfg(d, "--pg-level 7\n")))
        self.assertEqual(r["pg_level"], 7)
        self.assertIsNone(r["basis"])


class Methods(unittest.TestCase):
    def paragraph(self, d, extra="", version="2.7.0"):
        prov = os.path.join(d, "search_provenance.json")
        with open(prov, "w") as fh:
            json.dump({"engine": "diann", "version": version, "params_file": cfg(d, extra),
                       "search_mode": "parallel_5step"}, fh)
        return make_methods.search_paragraph(make_methods.search_record(search_prov=prov))

    def test_the_default_is_named_as_the_default(self):
        with tempfile.TemporaryDirectory() as d:
            s = self.paragraph(d)
        self.assertIn("implicit protein grouping by genes", s)
        self.assertIn("DIA-NN's default; no --pg-level was set", s)

    def test_another_version_does_not_borrow_the_measurement(self):
        with tempfile.TemporaryDirectory() as d:
            s = self.paragraph(d, version="2.6.1")
        self.assertIn("DIA-NN's default implicit protein grouping (no --pg-level was set; "
                      "measured as genes for DIA-NN 2.7.0)", s)

    def test_a_chosen_level_is_stated(self):
        with tempfile.TemporaryDirectory() as d:
            s = self.paragraph(d, "--pg-level 0\n")
        self.assertIn("implicit protein grouping by isoform IDs (--pg-level 0)", s)
        with tempfile.TemporaryDirectory() as d:
            s = self.paragraph(d, "--no-prot-inf\n")
        self.assertIn("Protein inference was switched off (--no-prot-inf)", s)


class Provenance(unittest.TestCase):
    def test_search_provenance_records_it(self):
        import test_single_shot_mass_accuracy as ss
        h = ss.SingleShotMassAccTests("test_sbatch_single_shot_measures_before_the_search_and_pins_it")
        with tempfile.TemporaryDirectory() as d:
            p, out = h._run_search(d, *h._setup(d), "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            with open(os.path.join(out, "search_provenance.json")) as fh:
                rec = json.load(fh)["protein_inference"]
        self.assertEqual((rec["pg_level"], rec["basis"], rec["passed"]), (2, "genes", None))


class Documented(unittest.TestCase):
    def test_the_engine_reference_says_what_the_level_does_and_when_to_set_it(self):
        with open(os.path.join(SKILL, "references", "search-engines.md"), encoding="utf-8") as fh:
            t = " ".join(fh.read().split())
        for s in ("--pg-level", "genes", "from step 1", "gg_matrix", "--overrides"):
            self.assertIn(s, t)

    def test_the_step5_warning_is_explained(self):
        with open(os.path.join(SKILL, "references", "watcher.md"), encoding="utf-8") as fh:
            t = " ".join(fh.read().split())
        self.assertIn("protein inference is enabled while empirical library is used", t)
        self.assertIn("Protein.Ids", t)


if __name__ == "__main__":
    unittest.main(verbosity=2)
