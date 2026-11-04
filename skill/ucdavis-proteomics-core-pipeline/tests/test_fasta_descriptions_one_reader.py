#!/usr/bin/env python3
"""One reader of a FASTA header's accession and description (C1 review, F5a; rule 3).

normalization_check.fasta_descriptions (the carboxylase panel's names, 2.14.2 I26) parsed FASTA
headers itself, beside fetch_fasta._header_ids, and the two read a `>Cont_X|ENTRY` header
differently: the second field only for three or more |-fields, so normalization_check's
accession was "Cont_X|ENTRY" where fetch_fasta's (what DIA-NN reports) is "Cont_X". Now
fetch_fasta.fasta_descriptions reads them -- the accession by _header_ids, the description by
_protein_name -- and normalization_check's is that one.

Stdlib only; hermetic (a synthetic FASTA).
"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fetch_fasta as ff  # noqa: E402
import normalization_check as nc  # noqa: E402

HEADERS = {
    ">Cont_X|ENTRY Serum albumin OS=Bos taurus OX=9913 GN=ALB": "Serum albumin",
    ">sp|P12345|NAME_HUMAN Pyruvate carboxylase, mitochondrial OS=Homo sapiens OX=9606 GN=PC "
    "PE=1 SV=2": "Pyruvate carboxylase, mitochondrial",
    ">XP_012345.1 keratin, type II cytoskeletal": "keratin, type II cytoskeletal",
}


class OneReader(unittest.TestCase):
    def setUp(self):
        fd, self.fasta = tempfile.mkstemp(suffix=".fasta")
        with os.fdopen(fd, "w") as fh:
            for h in HEADERS:
                fh.write(f"{h}\nMSTNPEPTIDEK\n")
        self.addCleanup(os.unlink, self.fasta)
        self.want = {ff._header_ids(h)[0]: d for h, d in HEADERS.items()}

    def test_the_accessions_are_the_ones_fetch_fasta_reads(self):
        self.assertEqual(sorted(self.want), ["Cont_X", "P12345", "XP_012345.1"])
        self.assertEqual(ff.fasta_descriptions(self.fasta, set(self.want)), self.want)

    def test_normalization_check_reads_them_the_same_way(self):
        self.assertEqual(nc.fasta_descriptions(self.fasta, set(self.want)), self.want)

    def test_only_the_asked_accessions_and_no_fasta(self):
        self.assertEqual(ff.fasta_descriptions(self.fasta, {"P12345"}),
                         {"P12345": "Pyruvate carboxylase, mitochondrial"})
        self.assertEqual(nc.fasta_descriptions(None, {"P12345"}), {})


if __name__ == "__main__":
    unittest.main()
