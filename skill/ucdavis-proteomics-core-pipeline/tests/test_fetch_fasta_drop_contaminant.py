#!/usr/bin/env python3
"""fetch_fasta.py --drop-contaminant leaves a named contaminant entry out of the database.

audit_results.py's contaminant_evidence check (2.12) names a contaminant group that rests only on
the sample's own peptides -- horse myoglobin Cont_P68082 in human heart, whose peptides there were
human myoglobin's (one differing only by I/L). The build keeps that entry by its rule (9 peptides
of its own), so the fix was a hand-made copy of the contaminant set passed with
--contaminants-path, which loses the set's name and is not replayed. Now:

  * --drop-contaminant <accession> (repeatable; Cont_P68082 or P68082) removes the entry from the
    contaminant set and from a supplied database's own contaminant entries, recorded in the
    sidecar (drop_contaminant_requested, contaminants_dropped_by_name with each one's source);
  * an accession that matches no contaminant entry stops the build;
  * the Methods' database paragraph names them (counted above ADDED_SEQUENCES_NAMED_MAX), and
    reproduce.sh replays the flag. Offline: --path and --contaminants-path.
"""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fetch_fasta as ff  # noqa: E402
import make_methods as mm  # noqa: E402

# UniProt REST, 2026-10-15
MYG_HUMAN = ("MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE"
             "DLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKH"
             "PGDFGADAQGAMNKALELFRKDMASNYKELGFQG")
MYG_HORSE = ("MGLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASE"
             "DLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKH"
             "PGDFGADAQGAMTKALELFRNDIAAKYKELGFQG")
PROTEOME = (f">sp|P02144|MYG_HUMAN Myoglobin OS=Homo sapiens OX=9606 GN=MB PE=1 SV=2\n{MYG_HUMAN}\n"
            ">sp|P04406|G3P_HUMAN GAPDH OS=Homo sapiens OX=9606 GN=GAPDH\n"
            "MGKVKVGVNGFGRIGRLVTRAAFNSGKVDIVAINDPFIDLNYMVYMFQYDSTHGKFHGTVK\n")
CONT = (f">sp|Cont_P68082|MYG_HORSE Myoglobin OS=Equus caballus OX=9796 GN=MB PE=1 SV=2\n"
        f"{MYG_HORSE}\n"
        ">sp|Cont_P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB\n"
        "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVK\n"
        ">sp|Cont_P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823\nFPTDDDDKIVGGYTCAANSIPYQVSLNSG\n")


class DropContaminant(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.proteome = self.write("proteome.fasta", PROTEOME)
        self.cont = self.write("cont.fasta", CONT)
        self.out = os.path.join(self.root, "out", "search.fasta")

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, text):
        p = os.path.join(self.root, name)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def fetch(self, *args, path=None, contaminants=True):
        cont = ["--contaminants", "universal", "--contaminants-path", self.cont] if contaminants \
            else ["--contaminants", "none"]
        argv = ["fetch_fasta.py", "fetch", "--path", path or self.proteome, "--organism",
                "Homo sapiens", "--taxid", "9606", *cont, *args, "--out", self.out]
        err = io.StringIO()
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(err):
            try:
                rc = ff.main()
            except SystemExit as e:
                rc = e.code
        if rc not in (0, None):
            return rc, None, None, err.getvalue()
        with open(self.out + ".meta.json") as fh, open(self.out) as fa:
            return rc, json.load(fh), fa.read(), err.getvalue()

    def test_without_the_flag_horse_myoglobin_stays(self):
        rc, m, fasta, err = self.fetch()
        self.assertIn(rc, (0, None), err)
        self.assertIn("Cont_P68082", fasta)
        self.assertEqual((m["drop_contaminant_requested"], m["contaminants_dropped_by_name"]), ([], []))
        self.assertEqual(mm.named_contaminants_sentence(m), "")

    def test_a_named_entry_leaves_the_set_and_is_recorded(self):
        for given in ("Cont_P68082", "P68082", "cont_p68082"):
            with self.subTest(given=given):
                rc, m, fasta, err = self.fetch("--drop-contaminant", given)
                self.assertIn(rc, (0, None), err)
                self.assertNotIn("Cont_P68082", fasta)
                self.assertIn("Cont_P02769", fasta)
                self.assertIn(">sp|P02144|MYG_HUMAN", fasta)
                self.assertEqual(m["drop_contaminant_requested"], [given])
                self.assertEqual(m["contaminants_dropped_by_name"],
                                 [{"cont_acc": "Cont_P68082", "cont_entry": "MYG_HORSE",
                                   "cont_gene": "MB", "source": "contaminant_set"}])
                self.assertEqual((m["n_contaminants_in_set"], m["n_contaminants_appended"]), (3, 2))
                self.assertIn("--drop-contaminant: left out Cont_P68082 (MYG_HORSE, contaminant set)",
                              err)

    def test_several_and_from_a_supplied_database(self):
        db = self.write("db.fasta", PROTEOME + CONT)
        rc, m, fasta, err = self.fetch("--drop-contaminant", "Cont_P68082", "--drop-contaminant",
                                       "P02769", path=db, contaminants=False)
        self.assertIn(rc, (0, None), err)
        self.assertNotIn("Cont_P68082", fasta)
        self.assertNotIn("Cont_P02769", fasta)
        self.assertIn("Cont_P00761", fasta)
        self.assertEqual([(r["cont_acc"], r["source"]) for r in m["contaminants_dropped_by_name"]],
                         [("Cont_P68082", "supplied_database"), ("Cont_P02769", "supplied_database")])
        self.assertEqual(m["n_contaminants_already_present"], 1)

    def test_an_accession_that_matches_nothing_stops_the_build(self):
        rc, m, _f, err = self.fetch("--drop-contaminant", "Cont_P68082", "--drop-contaminant",
                                    "Cont_Q99999")
        self.assertIsNone(m)
        self.assertIn("--drop-contaminant Cont_Q99999: no contaminant entry in this build has "
                      "that accession", str(rc))
        self.assertFalse(os.path.exists(self.out))

    def test_no_contaminant_set_is_said_as_such(self):
        rc, m, _f, _e = self.fetch("--drop-contaminant", "Cont_P68082", contaminants=False)
        self.assertIsNone(m)
        self.assertIn("no contaminant set was fetched for this build (--contaminants none, and the "
                      "supplied database's own contaminant entries do not include it)", str(rc))
        self.assertNotIn("has that accession", str(rc))

    def test_a_target_accession_is_not_a_contaminant(self):
        rc, m, _f, _e = self.fetch("--drop-contaminant", "P02144")
        self.assertIn("no contaminant entry in this build has that accession", str(rc))

    def test_the_methods_name_them(self):
        _, m, _, _ = self.fetch("--drop-contaminant", "Cont_P68082")
        self.assertEqual(mm.named_contaminants_sentence(m),
                         " 1 contaminant entry (Cont_P68082 (MYG_HORSE)) was left out of the "
                         "database by name.")
        many = dict(m, contaminants_dropped_by_name=[
            {"cont_acc": f"Cont_Q{i:05d}", "cont_entry": f"C{i}_BOVIN"} for i in range(11)])
        self.assertEqual(mm.named_contaminants_sentence(many),
                         " 11 contaminant entries were left out of the database by name.")
        md = os.path.join(self.root, "methods.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"), "--raw",
                            os.path.join(self.root, "FL_sample01.raw"), "--fasta-meta",
                            self.out + ".meta.json", "--out", md],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(md, encoding="utf-8") as fh:
            self.assertIn("1 contaminant entry (Cont_P68082 (MYG_HORSE)) was left out of the "
                          "database by name.", fh.read())

    def test_reproduce_sh_replays_it(self):
        _, m, _, _ = self.fetch("--drop-contaminant", "Cont_P68082")
        wm = self.write("workflow.manifest.json", "{}")
        out = os.path.join(self.root, "repro")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir", out,
                            "--workflow-manifest", wm, "--engine", "diann", "--acquisition", "DIA",
                            "--setup-json", os.path.join(self.root, "none.json"),
                            "--fasta-info-file", self.out + ".meta.json"],
                           capture_output=True, text=True, cwd=self.root, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out, "reproduce.sh")) as fh:
            sh = fh.read()
        line = next(ln for ln in sh.splitlines() if "--contaminants " in ln and "--out" in ln)
        self.assertIn("--drop-contaminant Cont_P68082", line)
        self.assertIn("Contaminant entry Cont_P68082 was left out by name (--drop-contaminant).", sh)


if __name__ == "__main__":
    unittest.main()
