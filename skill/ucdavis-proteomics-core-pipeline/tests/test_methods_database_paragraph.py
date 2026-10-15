#!/usr/bin/env python3
"""The Sequence-database paragraph of methods.md: added sequences counted, not listed, when there
are many; the database's build warnings and every server path kept out of the client's Methods.

Staff report (2026-10-10): a database of an NCBI proteome plus a second organism's 4,506-entry
UniProt proteome added with fetch_fasta.py --add-fasta. make_methods.py listed all 4,506
accessions in the Sequence-database paragraph (a 145 KB line; methods.md 154 KB) and appended a
block quote "Database build warnings (resolve before publication)" carrying an /nfs server path.
It was nearly sent to the client.

Now: added sequences are named one by one only up to ADDED_SEQUENCES_NAMED_MAX (the bait / tag
case), else counted with the file name and SHA-256 -- likewise the contaminant entries removed or
kept beside them; build warnings go to stderr and the printed JSON (database_build_warnings), never
to methods.md; and no absolute path from any record reaches methods.md. stdlib only, no network.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import make_methods as mm  # noqa: E402

SERVER_PATHS = ("/quobyte", "/nfs", "/home/")
ADDED_FILE = "/nfs/lssc0/flinders/Data/lab/service/PROT_0000/second_organism.fasta"
SHA = "c0ffee" * 10 + "abcd"


def entries(n, prefix="A0A"):
    return [{"accession": f"{prefix}{i:07d}", "name": f"tr|{prefix}{i:07d}|P{i}_METSP",
             "header": f"tr|{prefix}{i:07d}|P{i}_METSP protein", "length": 250,
             "sequence_sha256": f"{i:064x}"} for i in range(n)]


def cont(i):
    return {"cont_acc": f"Cont_Q{i:05d}", "cont_entry": f"C{i}_BOVIN", "target_acc": f"A0A{i:07d}",
            "shared_peptides": [f"PEPTIDE{i}K", f"SHARED{i}R"]}


def meta(n_added=4506, n_gone=0, n_kept=0, warnings=()):
    return {"fasta": "/quobyte/proteomics-grp/users/someone/PROT_0000/search.fasta",
            "organism": "Methylococcus sp.", "proteome": "UP000000001", "uniprot_release": "2026_04",
            "proteome_type": "Reference proteome", "content_used": "one_per_gene",
            "n_proteome": 2900, "n_contaminants_appended": 381, "contaminant_set": "universal",
            "n_contaminants_dropped_as_target": 0, "min_unique_peptides": 2,
            "added_sequences": [{"file": ADDED_FILE, "sha256": SHA, "n_entries": n_added,
                                 "entries": entries(n_added)}],
            "contaminants_dropped_for_added_sequences": [cont(i) for i in range(n_gone)],
            "contaminants_sharing_peptides_with_added_sequences": [cont(100 + i) for i in range(n_kept)],
            "warnings": list(warnings)}


class AddedSequencesSentence(unittest.TestCase):
    def test_a_whole_added_proteome_is_counted_with_its_file_and_hash(self):
        s = mm.added_sequences_sentence(meta())
        self.assertIn(f"4,506 user-supplied sequences were added to the database as target entries "
                      f"(second_organism.fasta, SHA-256 {SHA}).", s)
        self.assertNotIn("A0A0000001", s)
        self.assertLess(len(s), 400)

    def test_a_bait_and_a_tag_are_named(self):
        m = meta(2)
        s = mm.added_sequences_sentence(m)
        self.assertIn("2 user-supplied sequences (tr|A0A0000000|P0_METSP, tr|A0A0000001|P1_METSP) "
                      "were added", s)

    def test_the_limit_is_ten(self):
        self.assertEqual(mm.ADDED_SEQUENCES_NAMED_MAX, 10)
        self.assertIn("tr|A0A0000009|P9_METSP", mm.added_sequences_sentence(meta(10)))
        s = mm.added_sequences_sentence(meta(11))
        self.assertNotIn("tr|A0A0000000|", s)
        self.assertIn("11 user-supplied sequences were added", s)

    def test_a_sidecar_without_its_entries_is_counted_from_n_entries(self):
        m = meta(4506)
        m["added_sequences"] = [{k: v for k, v in m["added_sequences"][0].items() if k != "entries"}]
        self.assertIn("4,506 user-supplied sequences were added", mm.added_sequences_sentence(m))

    def test_many_removed_or_kept_contaminants_are_counted(self):
        s = mm.added_sequences_sentence(meta(4506, n_gone=12, n_kept=11))
        self.assertIn("12 contaminant entries sharing most of their peptides with them were removed "
                      "from the library", s)
        self.assertIn("11 contaminant entries sharing fewer of their peptides with them were kept, "
                      "and the peptides shared are ambiguous between the contaminant and the added "
                      "sequence (22 peptides).", s)
        self.assertNotIn("Cont_Q00000", s)
        self.assertNotIn("PEPTIDE100K", s)

    def test_a_few_are_still_named(self):
        s = mm.added_sequences_sentence(meta(1, n_gone=1, n_kept=1))
        self.assertIn("(Cont_Q00000 (C0_BOVIN)) was removed", s)
        self.assertIn("Cont_Q00100 (C100_BOVIN) and A0A0000100: PEPTIDE100K, SHARED100R.", s)


class NoWarningsOrServerPathsInMethods(unittest.TestCase):
    WARNINGS = (f"added sequence A0A0000003 ({ADDED_FILE}) is identical to database entry Q00003: "
                f"the two share every peptide, so the search cannot tell them apart",
                "a pre-staged copy was used: /quobyte/proteomics-grp/MRS/UP000000001_0000.fasta "
                "(release unknown)",
                "the user's copy at /home/someone/core/PROT_0000/old.fasta was not used")

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.meta_path = os.path.join(self.root, "search.fasta.meta.json")
        with open(self.meta_path, "w") as fh:
            json.dump(meta(4506, n_gone=2, warnings=self.WARNINGS), fh)
        self.out = os.path.join(self.root, "methods.md")
        self.r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                                 "--raw", os.path.join(self.root, "FL_sample01.raw"),
                                 "--fasta-meta", self.meta_path, "--out", self.out],
                                capture_output=True, text=True, timeout=300)
        self.assertEqual(self.r.returncode, 0, self.r.stderr)
        with open(self.out, encoding="utf-8") as fh:
            self.text = fh.read()

    def tearDown(self):
        self._td.cleanup()

    def test_methods_md_has_no_server_path_and_no_build_warnings(self):
        for p in SERVER_PATHS:
            self.assertNotIn(p, self.text)
        self.assertNotIn(self.root, self.text)
        self.assertNotIn("Database build warnings", self.text)
        self.assertNotIn("identical to database entry", self.text)
        self.assertIn("(second_organism.fasta, SHA-256", self.text)
        self.assertLess(os.path.getsize(self.out), 20 * 1024)

    def test_the_warnings_go_to_the_agent(self):
        out = json.loads(self.r.stdout)
        self.assertEqual(out["database_build_warnings"], list(self.WARNINGS))
        self.assertTrue(out["next"].startswith("Tell the user the database build warnings"))
        for w in self.WARNINGS:
            self.assertIn(f"database build warning (not in methods.md): {w}", self.r.stderr)

    def test_no_warnings_no_mention(self):
        with open(self.meta_path, "w") as fh:
            json.dump(meta(2), fh)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                            "--raw", os.path.join(self.root, "FL_sample01.raw"),
                            "--fasta-meta", self.meta_path, "--out", self.out],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        out = json.loads(r.stdout)
        self.assertEqual(out["database_build_warnings"], [])
        self.assertTrue(out["next"].startswith("Verify the draft"))

    def test_the_drop_notes_described_in_the_paragraph_are_not_warnings(self):
        m = meta(2, warnings=("dropped 3 contaminant entries ...",))
        m["contaminants_dropped_note"] = "dropped 3 contaminant entries ..."
        self.assertEqual(mm.database_build_warnings(m), [])


class ClientText(unittest.TestCase):
    def test_paths_are_cut_to_the_file_name(self):
        for raw, want in (("from /quobyte/proteomics-grp/MRS/x.fasta.", "from x.fasta."),
                          ("(/nfs/lssc0/flinders/Data/a.d)", "(a.d)"),
                          ("/home/someone/core/PROT_0000/search.fasta", "search.fasta"),
                          ("~/core/PROT_0000/in.fasta", "in.fasta"),
                          ("C:\\Users\\someone\\bait.fasta", "bait.fasta")):
            self.assertEqual(mm.client_text(raw), (want, 1), raw)

    def test_units_ratios_and_urls_are_not_paths(self):
        s = ("m/z 300-1200, 1/K0 0.6-1.6 V·s/cm², Trypsin/P, 300 nL/min, 0.1% (v/v), and/or, "
             "https://proteomics.ucdavis.edu/instrument-grant-acknowledgments, "
             "doi:10.1021/acs.jproteome.2c00145, ~20 min/sample, a / b")
        self.assertEqual(mm.client_text(s), (s, 0))


if __name__ == "__main__":
    unittest.main()
