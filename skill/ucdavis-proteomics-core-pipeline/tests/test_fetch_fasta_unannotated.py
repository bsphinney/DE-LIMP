#!/usr/bin/env python3
"""A supplied FASTA whose headers name no protein is said so at fetch (2.14.2 triage I07, the
warning only).

The coral study's session (2026-10-22) searched a gene-model FASTA whose headers are bare IDs
(">adig_s0001.g1.t1"): no GN=, no description. Every table, the volcano and the report then named
proteins by those IDs, and nothing said so until the report was read. fetch --path now warns and
records `unannotated: true` in the sidecar. (Annotating them -- DIAMOND against Swiss-Prot -- is
2.15.)

All offline.
"""
import contextlib
import io
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fetch_fasta as ff  # noqa: E402

CORAL = (">adig_s0001.g1.t1\nMPEPTIDEKAAAGGGLLLR\n"
         ">adig_s0002.g7.t2\nMSTVWYKCDEFGHIKLMNPR\n"
         ">adig_s0013.g21.t1\nMKLVNEATEFAKTCVADESAENCDK\n")
NCBI = (">XP_000001.1 synthetic protein one\nMPEPTIDEKAAAGGGLLLR\n"
        ">XP_000002.1 synthetic protein two\nMSTVWYKCDEFGHIKLMNPR\n")
UNIPROT_BARE = (">sp|P12345|ABC1_HUMAN\nMPEPTIDEKAAAGGGLLLR\n"
                ">tr|A0A024R161|A0A024R161_HUMAN\nMSTVWYKCDEFGHIKLMNPR\n")


class Unannotated(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.out = os.path.join(self.root, "out", "search.fasta")

    def tearDown(self):
        self._td.cleanup()

    def fetch(self, text):
        db = os.path.join(self.root, "db.fasta")
        with open(db, "w") as fh:
            fh.write(text)
        argv = ["fetch_fasta.py", "fetch", "--path", db, "--organism", "none",
                "--contaminants", "none", "--out", self.out]
        err = io.StringIO()
        with mock.patch.object(sys, "argv", argv), \
                mock.patch.object(ff, "_get_json", side_effect=AssertionError("network")), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            try:
                rc = ff.main()
            except SystemExit as e:
                rc = e.code
        self.assertIn(rc, (0, None), err.getvalue())
        with open(self.out + ".meta.json") as fh:
            return json.load(fh), err.getvalue()

    def test_bare_gene_model_ids_warn_and_are_recorded(self):
        meta, err = self.fetch(CORAL)
        self.assertIs(meta["unannotated"], True)
        self.assertIn("no gene names (GN=) and no protein descriptions", err)
        self.assertTrue(any("no gene names (GN=)" in w for w in meta["warnings"]), meta["warnings"])
        self.assertIn(">adig_s0001.g1.t1", err)

    def test_described_headers_are_annotated(self):                     # (pin)
        meta, err = self.fetch(NCBI)
        self.assertIs(meta["unannotated"], False)
        self.assertNotIn("no gene names (GN=)", err)

    def test_a_uniprot_entry_name_counts_as_a_name(self):
        meta, _ = self.fetch(UNIPROT_BARE)
        self.assertIs(meta["unannotated"], False)

    def test_the_rule_on_its_own(self):
        self.assertTrue(ff.fasta_unannotated(CORAL)[0])
        self.assertFalse(ff.fasta_unannotated(NCBI)[0])
        self.assertFalse(ff.fasta_unannotated("")[0])                   # nothing to judge

    # 2.14.2 review: one named entry -- an added bait -- silenced the warning, and key=value
    # tokens (Trinity's len=, path=) passed for a description
    def test_one_named_bait_does_not_silence_it(self):
        flag, n_unnamed, n, _ = ff.fasta_unannotated(
            CORAL + ">EGFP enhanced green fluorescent protein GN=EGFP\nMVSKGEELFTGVVPILVELDGDVNGHK\n")
        self.assertTrue(flag)
        self.assertEqual((n_unnamed, n), (3, 4))

    def test_key_value_tokens_are_not_a_description(self):
        trinity = (">TRINITY_DN1_c0_g1_i1 len=300 path=[0:0-299] gene=TRINITY_DN1_c0_g1\nMPEPTIDEK\n"
                   ">TRINITY_DN2_c0_g1_i1 len=210 path=[0:0-209]\nMSTVWYKCDEFGHIKLMNPR\n")
        self.assertTrue(ff.fasta_unannotated(trinity)[0])

    def test_a_mostly_named_database_is_annotated(self):
        self.assertFalse(ff.fasta_unannotated(NCBI * 3 + ">adig_s0001.g1.t1\nMK\n")[0])

    def test_the_sidecar_says_how_many(self):
        meta, err = self.fetch(CORAL + ">EGFP enhanced green fluorescent protein\nMVSKGEELFTGVVPILVELDGDVNGHK\n")
        self.assertIs(meta["unannotated"], True)
        self.assertIn("3 of 4 entries", err)


if __name__ == "__main__":
    unittest.main()
