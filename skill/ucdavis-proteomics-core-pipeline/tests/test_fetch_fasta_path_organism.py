#!/usr/bin/env python3
"""A --path database records the organism the user named -- and asks when nobody did.

A staff build (2026-10-01, skill 2.9.0): `fetch --path <a supplied FASTA> --ncbi-organism '...'
--ncbi-taxid ...` exited 0 with no warning, but the sidecar said organism "", taxid 0,
organism_source "none". FRAN files a search by that organism, and Methods and the report
read it too, so the search landed under no species and the Methods named none -- while the
flags looked accepted.

Guards, all offline:
  * --path honours --organism/--taxid and the --ncbi-* spellings, recorded as the user's;
  * --path with no organism stops and asks; `--organism none` is the answer for a database
    with no single organism;
  * an organism flag the chosen source would ignore (--proteome) is refused, as are two
    spellings that disagree;
  * Methods names the organism of a supplied database.
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

DB = (">XP_000001.1 synthetic protein one\nMPEPTIDEKAAAGGGLLLR\n"
      ">XP_000002.1 synthetic protein two\nMSTVWYKCDEFGHIKLMNPR\n")


class PathOrganism(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.db = os.path.join(self.root, "assembly_one_per_gene.fasta")
        with open(self.db, "w") as fh:
            fh.write(DB)
        self.out = os.path.join(self.root, "out", "search.fasta")

    def tearDown(self):
        self._td.cleanup()

    def fetch(self, *args):
        argv = ["fetch_fasta.py", "fetch", *args, "--contaminants", "none", "--out", self.out]
        err = io.StringIO()
        # No organism lookup may reach the network here: a call means the flag was ignored.
        with mock.patch.object(sys, "argv", argv), \
                mock.patch.object(ff, "_get_json", side_effect=AssertionError("network")), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            try:
                rc = ff.main()
            except SystemExit as e:
                rc = e.code
        meta = None
        if os.path.exists(self.out + ".meta.json"):
            with open(self.out + ".meta.json") as fh:
                meta = json.load(fh)
        return rc, meta

    def test_the_reported_spelling_is_honoured(self):
        rc, m = self.fetch("--path", self.db, "--ncbi-organism", "Synthetica exempli",
                           "--ncbi-taxid", "424242")
        self.assertEqual(rc, 0)
        self.assertEqual((m["organism"], m["taxid"]), ("Synthetica exempli", 424242))
        self.assertEqual(m["organism_source"], "user (--ncbi-organism/--ncbi-taxid)")
        self.assertEqual(m["proteome_type"], "user-supplied FASTA (--path)")
        self.assertEqual(m["source"], f"override:{self.db}")

    def test_organism_and_taxid_are_recorded_as_the_users(self):
        rc, m = self.fetch("--path", self.db, "--organism", "Synthetica exempli", "--taxid", "424242")
        self.assertEqual(rc, 0)
        self.assertEqual((m["organism"], m["taxid"]), ("Synthetica exempli", 424242))
        self.assertEqual(m["organism_source"], "user (--organism/--taxid)")

    def test_no_organism_stops_and_asks(self):
        rc, m = self.fetch("--path", self.db)
        self.assertNotEqual(rc, 0)
        self.assertIn("--organism", str(rc))
        self.assertIn("--organism none", str(rc))
        self.assertIsNone(m)
        self.assertFalse(os.path.exists(self.out))

    def test_none_is_the_answer_for_no_single_organism(self):
        rc, m = self.fetch("--path", self.db, "--organism", "none")
        self.assertEqual(rc, 0)
        self.assertEqual((m["organism"], m["taxid"]), ("", 0))
        self.assertEqual(m["organism_source"], "user: no single organism (--organism none)")
        rc, _ = self.fetch("--path", self.db, "--organism", "none", "--taxid", "9606")
        self.assertIn("contradict", str(rc))

    def test_a_name_and_taxid_of_different_organisms_are_refused(self):
        rc, m = self.fetch("--path", self.db, "--organism", "Homo sapiens", "--taxid", "10090")
        self.assertIn("contradict each other: taxid 10090 is Mus musculus", str(rc))
        self.assertIsNone(m)
        # a strain or a shorter name of the same organism is not a contradiction
        rc, m = self.fetch("--path", self.db, "--organism", "Escherichia coli", "--taxid", "83333")
        self.assertEqual(rc, 0)
        self.assertEqual((m["organism"], m["taxid"]), ("Escherichia coli", 83333))

    def test_a_curated_taxid_alone_fills_the_name_and_says_so(self):
        rc, m = self.fetch("--path", self.db, "--taxid", "10090")
        self.assertEqual(rc, 0)
        self.assertEqual((m["organism"], m["taxid"]), ("Mus musculus", 10090))
        self.assertEqual(m["organism_source"], "user (--taxid) + curated_table name for the taxid")
        rc, _ = self.fetch("--path", self.db, "--taxid", "424242")
        self.assertIn("give --organism", str(rc))

    def test_a_flag_the_source_ignores_is_refused(self):
        for flags in (["--organism", "Homo sapiens"], ["--taxid", "9606"],
                      ["--ncbi-organism", "Homo sapiens"]):
            rc, _ = self.fetch("--proteome", "UP000005640", *flags)
            self.assertIn("would be ignored", str(rc), flags)

    def test_two_spellings_that_disagree_are_refused(self):
        rc, _ = self.fetch("--path", self.db, "--organism", "Synthetica exempli",
                           "--ncbi-organism", "Synthetica altera")
        self.assertIn("disagree", str(rc))
        rc, m = self.fetch("--path", self.db, "--organism", "Synthetica exempli",
                           "--ncbi-organism", "synthetica EXEMPLI")
        self.assertEqual(rc, 0)

    def test_ncbi_accession_takes_organism_as_well(self):
        with mock.patch.object(ff, "ncbi_download_proteome",
                               return_value=(DB, os.path.join(self.root, "p.faa"), "https://x")):
            rc, m = self.fetch("--ncbi-accession", "GCF_000000001.1", "--organism",
                               "Synthetica exempli", "--taxid", "424242")
        self.assertEqual(rc, 0)
        self.assertEqual((m["organism"], m["taxid"]), ("Synthetica exempli", 424242))
        self.assertTrue(m["proteome_type"].startswith("NCBI RefSeq assembly proteins"))

    def methods(self, meta):
        mp = os.path.join(self.root, "meta.json")
        with open(mp, "w") as fh:
            json.dump(meta, fh)
        md = os.path.join(self.root, "methods.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                            "--raw", os.path.join(self.root, "FL_run.raw"), "--fasta-meta", mp,
                            "--out", md], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(md) as fh:
            return fh.read()

    def test_methods_names_the_organism_of_a_supplied_database(self):
        _, m = self.fetch("--path", self.db, "--organism", "Synthetica exempli", "--taxid", "424242")
        text = self.methods(m)
        self.assertIn("searched against a supplied Synthetica exempli sequence database "
                      "(search.fasta; taxid 424242; 2 sequences)", text)
        _, m = self.fetch("--path", self.db, "--organism", "none")
        self.assertIn("a supplied sequence database with no single source organism "
                      "(search.fasta; 2 sequences)", self.methods(m))


if __name__ == "__main__":
    unittest.main()
