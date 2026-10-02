#!/usr/bin/env python3
"""Independent 2.10 review (wrong-results risks): fetch_fasta.py --add-fasta and --path organism.

Each test here FAILS on release/skill-2.10.0 (40c5b9f) and describes the behaviour a fix must give.

  * A supplied database made by FragPipe tags its contaminants contam_ (FRAGPIPE_CONT_TAG), and
    the --add-fasta rule looks only at Cont_ entries: wild-type GFP stays beside an EGFP bait,
    neither removed nor flagged, and nothing in the sidecar or Methods says so.
  * Text before the first '>' of an --add-fasta file is written verbatim after the proteome, so it
    becomes part of the LAST proteome entry's sequence in the searched database.
  * The --add-fasta rule and the ordinary peptide rule (MIN_UNIQUE_PEPTIDES) are judged apart: a
    contaminant whose peptides are split between a proteome protein and the bait keeps ONE peptide
    of its own in the whole database and still stays, taking peptides from both.
  * --path --organism / --taxid that contradict the curated organism table are recorded as given
    (FRAN files the search under a taxid that is not the organism Methods names).

Offline; no network.
"""
import argparse
import contextlib
import io
import json
import os
import random
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))
sys.path.insert(0, HERE)

import fetch_fasta as ff  # noqa: E402
from test_fetch_fasta_add_fasta import EGFP, WT_GFP  # noqa: E402

ACTB = (">sp|P60709|ACTB_HUMAN Actin OS=Homo sapiens OX=9606 GN=ACTB\n"
        "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPR\n")


class _Fetch(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.out = os.path.join(self.root, "out", "search.fasta")

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, text):
        p = os.path.join(self.root, name)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def fetch(self, proteome_text, added_text):
        argv = ["fetch_fasta.py", "fetch", "--path", self.write("p.fasta", proteome_text),
                "--organism", "Homo sapiens", "--taxid", "9606", "--contaminants", "none",
                "--add-fasta", self.write("a.fasta", added_text), "--out", self.out]
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            try:
                rc = ff.main()
            except SystemExit as e:
                rc = e.code
        if rc not in (0, None):
            return rc, None, None
        with open(self.out + ".meta.json") as fh, open(self.out) as fa:
            return rc, json.load(fh), fa.read()


class FragPipeTaggedSuppliedDatabase(_Fetch):
    def test_a_contam_tagged_wild_type_gfp_is_removed_or_flagged_beside_an_egfp_bait(self):
        # The same database with the Core's tag works (Cont_P42212 removed); FragPipe's tag is
        # ignored by the --add-fasta rule -- every other contaminant rule here reads both tags.
        db = ACTB + f">contam_sp|P42212|GFP_AEQVI Green fluorescent protein\n{WT_GFP}\n"
        rc, meta, fasta = self.fetch(db, f">sp|EGFP01|EGFP_SYN EGFP bait\n{EGFP}\n")
        self.assertIn(rc, (0, None))
        seen = ([r["cont_acc"] for r in meta["contaminants_dropped_for_added_sequences"]]
                + [r["cont_acc"] for r in
                   meta["contaminants_sharing_peptides_with_added_sequences"]])
        self.assertTrue(any("P42212" in acc for acc in seen),
                        "contam_sp|P42212 (wild-type GFP, 75% of its peptides in EGFP) was "
                        "neither removed nor flagged; the sidecar and Methods say nothing")


class TextBeforeTheFirstHeader(_Fetch):
    def test_a_preamble_line_never_joins_the_last_proteome_entry(self):
        rc, meta, fasta = self.fetch(
            ACTB, "EGFP construct from plasmid pX330\n>EGFP_construct\nMVSKGEELFTGVVPILVELDGDVNGHK\n")
        if rc not in (0, None):
            return                                    # refusing the file is a correct fix too
        recs = ff._fasta_records(fasta)
        actb = next(lines for h, lines in recs if "P60709" in h)
        self.assertEqual(ff._record_seq(actb), "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPR",
                         "the --add-fasta file's first line was appended to ACTB's sequence in "
                         "the searched database")


class RulesJudgedApart(unittest.TestCase):
    def test_a_contaminant_with_one_own_peptide_in_the_whole_database_leaves(self):
        random.seed(7)
        aa = "ACDEFGHLMNQSTVWY"                       # no K/R/P/I: one tryptic peptide each
        pep = ["".join(random.choice(aa) for _ in range(16)) + "K" for _ in range(7)]
        proteome = f">sp|Q00001|HUM_HUMAN human protein GN=HUM\n{''.join(pep[0:3])}\n"
        added = f">sp|BAIT01|BAIT_SYN bait\n{''.join(pep[3:6])}\n"
        contam = f">sp|Cont_Q99999|CON_BOVIN bovine protein\n{''.join(pep)}\n"
        used = ff.parse_enzymes(ff.DEFAULT_ENZYMES)
        # the order cmd_fetch applies them: the ordinary rule vs EVERY target -- the proteome and
        # the added sequences (the fix: judged together) -- then --add-fasta
        kept, _d, _e = ff.drop_target_contaminants(
            contam, proteome + ff.added_as_targets(ff._fasta_records(added)), used,
            ff.MIN_UNIQUE_PEPTIDES)
        kept, _d2, _flagged, _e2 = ff.drop_contaminants_near_added(
            kept, ff._fasta_records(added), used, "contaminant_set", only_tagged=False)
        own = ff.peptide_overlap(ff._fasta_records(contam),
                                 ff._fasta_records(proteome + added))[0]["n_unique_peptides"]
        self.assertEqual(own, 1)
        self.assertNotIn("Cont_Q99999", kept,
                         f"Cont_Q99999 has {own} peptide of its own among proteome + added "
                         f"sequences (MIN_UNIQUE_PEPTIDES = {ff.MIN_UNIQUE_PEPTIDES}) but stays, "
                         f"taking 3 peptides from Q00001 and 3 from the bait")


class PathOrganismContradiction(unittest.TestCase):
    def test_a_name_and_taxid_the_curated_table_contradicts_are_refused(self):
        a = argparse.Namespace(organism="Homo sapiens", taxid=10090, ncbi_organism="",
                               ncbi_taxid=0, path="/x.fasta", ncbi_accession="")
        with contextlib.redirect_stderr(io.StringIO()):
            try:
                got = ff.user_organism(a)
            except SystemExit:
                return
        self.fail(f"--organism 'Homo sapiens' --taxid 10090 (Mus musculus in ORGANISM_TAXIDS) "
                  f"was recorded as {got}")


if __name__ == "__main__":
    unittest.main()
