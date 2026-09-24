#!/usr/bin/env python3
"""A pre-staged HIVE proteome must still describe itself in <fasta>.meta.json.

gabrig, 2026-09-23 (human HeLa, HIVE Core member): `fetch --proteome UP000005640 --hive`
reused /quobyte/proteomics-grp/MRS/UP000005640_9606.fasta and wrote organism "", taxid 0,
uniprot_release "" -- although the proteome ID was known. make_methods.py (methods text) and
fran_deposit.py (the FRAN corpus row) both read organism/taxid from that sidecar.

Guards, all offline (the UniProt lookup is mocked):
  * organism/taxid come from UniProt when reachable, else from the filename's taxid + the
    curated table, and the sidecar says which;
  * the live UniProt release is NEVER recorded as the staged copy's release;
  * the one-per-gene inference passes 20,663 vs geneCount 20,652 and fails 147,520;
  * a staged file whose name carries no taxid still works.
"""
import contextlib
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
import urllib.error
from datetime import datetime, timezone
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fetch_fasta as ff  # noqa: E402

UP = "UP000005640"
LIVE_RELEASE = "2026_03"          # what UniProt's header says TODAY -- not the staged copy's
STAGED_MTIME = datetime(2025, 4, 25, 23, 25, 36, tzinfo=timezone.utc)   # the real file's date


def fasta(n, isoforms=0):
    recs = [f">sp|P{i:05d}|G{i}_HUMAN Protein {i} OS=Homo sapiens OX=9606 GN=G{i}\nPEPTIDEK\n"
            for i in range(n)]
    recs += [f">sp|P{i:05d}-2|G{i}_HUMAN Isoform 2 of Protein {i}\nPEPTIDER\n"
             for i in range(isoforms)]
    return "".join(recs)


def uniprot(gene_count, taxid=9606, organism="Homo sapiens"):
    """Stand-in for _get_json on /proteomes/UP000005640 -- shape verified against the live
    endpoint 2026-09-23 (proteinCount 147520, geneCount 20652, x-uniprot-release 2026_03)."""
    data = {"id": UP, "taxonomy": {"scientificName": organism, "taxonId": taxid},
            "proteomeType": "Reference proteome", "proteinCount": 147520,
            "geneCount": gene_count, "superkingdom": "eukaryota"}
    headers = {"x-uniprot-release": LIVE_RELEASE, "x-uniprot-release-date": "02-September-2026"}
    return mock.patch.object(ff, "_get_json", return_value=(data, headers))


def offline():
    return mock.patch.object(ff, "_get_json",
                             side_effect=urllib.error.URLError("[Errno 101] Network is unreachable"))


class StagedMeta(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.mrs = os.path.join(self.root, "MRS")
        os.makedirs(self.mrs)
        self.cont = os.path.join(self.root, "cont.fasta")
        # Sequences no synthetic target contains: a contaminant that IS a target protein is
        # dropped (test_fetch_fasta_contaminant_overlap.py), which is not what these test.
        with open(self.cont, "w") as fh:
            fh.write(">Cont_P02768|ALBU_HUMAN\nMKWVTFISLLFLFSSAYS\n"
                     ">Cont_P00761|TRYP_PIG\nFPTDDDDKIVGGYTCAANSIPYQVSLNSG\n")

    def tearDown(self):
        self._td.cleanup()

    def stage(self, name, text):
        p = os.path.join(self.mrs, name)
        with open(p, "w") as fh:
            fh.write(text)
        os.utime(p, (STAGED_MTIME.timestamp(), STAGED_MTIME.timestamp()))
        return p

    def fetch(self, *extra, proteome=UP):
        out = os.path.join(self.root, "out", "search.fasta")
        argv = ["fetch_fasta.py", "fetch", "--proteome", proteome, "--hive",
                "--contaminants", "universal", "--contaminants-path", self.cont,
                "--out", out, *extra]
        err = io.StringIO()
        with mock.patch.object(ff, "HIVE_MRS", self.mrs), mock.patch.object(sys, "argv", argv), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            rc = ff.main()
        self.assertEqual(rc, 0, err.getvalue())
        with open(out + ".meta.json") as fh:
            return json.load(fh), out + ".meta.json"

    # -- organism / taxid ------------------------------------------------------------------
    def test_staged_file_fills_organism_and_taxid_from_uniprot(self):
        self.stage("UP000005640_9606.fasta", fasta(50))
        with uniprot(gene_count=50):
            m, _ = self.fetch()
        self.assertEqual(m["organism"], "Homo sapiens")
        self.assertEqual(m["taxid"], 9606)
        self.assertEqual(m["proteome_type"], "Reference proteome")
        self.assertEqual(m["organism_source"], "uniprot_api:UP000005640")
        # Still the literal enum: provenance.py writes it into reproduce.sh as --content.
        self.assertEqual(m["content_used"], "as_staged")
        self.assertEqual(m["warnings"], [])

    def test_offline_fills_organism_from_filename_taxid_and_curated_table(self):
        self.stage("UP000005640_9606.fasta", fasta(50))
        with offline():
            m, _ = self.fetch()
        self.assertEqual(m["organism"], "Homo sapiens")
        self.assertEqual(m["taxid"], 9606)
        self.assertEqual(m["organism_source"], "filename_taxid+curated_table")
        self.assertIn("URLError", m["staged_file"]["uniprot_lookup_error"])
        # "Reference proteome" is UniProt's claim; offline we do not make it for them.
        self.assertEqual(m["proteome_type"], "")
        # No geneCount offline -> composition unchecked, and said so rather than assumed.
        self.assertEqual(m["content_check"]["verdict"], "unchecked")
        self.assertIsNone(m["content_inferred"])
        self.assertTrue(any("could not confirm" in w for w in m["warnings"]), m["warnings"])

    def test_filename_taxid_disagreeing_with_uniprot_is_flagged(self):
        self.stage("UP000005640_10090.fasta", fasta(50))
        with uniprot(gene_count=50):
            m, _ = self.fetch()
        self.assertEqual(m["taxid"], 9606)      # UniProt's answer wins...
        self.assertTrue(any("taxid 10090" in w for w in m["warnings"]), m["warnings"])  # ...loudly

    # -- release ---------------------------------------------------------------------------
    def test_live_release_is_never_recorded_as_the_staged_copys(self):
        staged = self.stage("UP000005640_9606.fasta", fasta(50))
        with uniprot(gene_count=50):
            m, meta_path = self.fetch()
        self.assertEqual(m["uniprot_release"], "")
        self.assertEqual(m["uniprot_release_date"], "")
        with open(meta_path) as fh:
            self.assertNotIn(LIVE_RELEASE, fh.read())
        self.assertTrue(m["staged_release_unknown"])
        st = m["staged_file"]
        self.assertTrue(st["release_unknown"])
        self.assertEqual(st["path"], staged)
        self.assertEqual(st["mtime_utc"], "2025-04-25T23:25:36Z")
        self.assertIn("2025-04-25", st["release_note"])
        # The STAGED file's hash: the output also holds the appended contaminants.
        with open(staged, "rb") as fh:
            self.assertEqual(st["sha256"], hashlib.sha256(fh.read()).hexdigest())
        self.assertNotEqual(st["sha256"], m["sha256"])

    def test_non_staged_builds_carry_no_staged_record(self):
        # --path override: not staged, so no stand-in fields claiming otherwise.
        p = os.path.join(self.root, "mine.fasta")
        with open(p, "w") as fh:
            fh.write(fasta(5))
        out = os.path.join(self.root, "o2", "search.fasta")
        argv = ["fetch_fasta.py", "fetch", "--path", p, "--contaminants", "none", "--out", out]
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(ff.main(), 0)
        with open(out + ".meta.json") as fh:
            m = json.load(fh)
        self.assertIsNone(m["staged_file"])
        self.assertFalse(m["staged_release_unknown"])
        self.assertIsNone(m["content_inferred"])
        self.assertEqual(m["organism_source"], "none")

    # -- one-per-gene inference ------------------------------------------------------------
    def test_gene_count_inference_on_the_real_human_numbers(self):
        # The Core's staged copy vs UniProt's live geneCount (both measured 2026-09-23).
        ok = ff.infer_staged_content(20663, 20652)
        self.assertEqual(ok["verdict"], "consistent_with_one_per_gene")
        self.assertIn("inferred", ok["note"])
        # The FULL human proteome (proteinCount) must never pass as one-per-gene.
        self.assertEqual(ff.infer_staged_content(147520, 20652)["verdict"],
                         "larger_than_one_per_gene")
        # Arabidopsis full/gene = 1.43x, the closest full set among the curated organisms.
        self.assertEqual(ff.infer_staged_content(39272, 27496)["verdict"],
                         "larger_than_one_per_gene")
        self.assertEqual(ff.infer_staged_content(20663, 20652, n_isoform=12)["verdict"],
                         "contains_isoforms")
        self.assertEqual(ff.infer_staged_content(20663, 0)["verdict"], "unchecked")

    def test_one_per_gene_inferred_is_labelled_and_separate_from_content_used(self):
        self.stage("UP000005640_9606.fasta", fasta(50))
        with uniprot(gene_count=50):
            m, _ = self.fetch()
        self.assertEqual(m["content_inferred"], "one_per_gene")
        self.assertEqual(m["content_check"]["n_entries"], 50)
        self.assertEqual(m["content_check"]["uniprot_gene_count"], 50)
        self.assertTrue(m["content_check"]["note"].startswith("one_per_gene (inferred:"))

    def test_full_set_staged_while_one_per_gene_requested_warns(self):
        self.stage("UP000005640_9606.fasta", fasta(60))
        with uniprot(gene_count=10):          # 6x the gene count: a full set
            m, _ = self.fetch("--content", "one_per_gene")
        self.assertIsNone(m["content_inferred"])
        self.assertEqual(m["content_check"]["verdict"], "larger_than_one_per_gene")
        self.assertTrue(any("does not look like the one_per_gene" in w for w in m["warnings"]),
                        m["warnings"])

    def test_undersized_staged_copy_warns(self):
        """A truncated / partial / reviewed-only copy was searched with warnings: [] (review
        2026-09-23): the smaller verdict raised nothing."""
        self.stage("UP000005640_9606.fasta", fasta(30))
        with uniprot(gene_count=50):          # 0.6x the gene count
            m, _ = self.fetch("--content", "one_per_gene")
        self.assertEqual(m["content_check"]["verdict"], "smaller_than_one_per_gene")
        self.assertIsNone(m["content_inferred"])
        self.assertTrue(any("FEWER entries" in w for w in m["warnings"]), m["warnings"])

    def test_isoform_accessions_fail_the_inference_even_at_the_right_count(self):
        self.stage("UP000005640_9606.fasta", fasta(45, isoforms=5))
        with uniprot(gene_count=50):
            m, _ = self.fetch()
        self.assertEqual(m["content_check"]["verdict"], "contains_isoforms")
        self.assertIsNone(m["content_inferred"])
        self.assertTrue(m["warnings"])

    def test_other_content_requested_says_it_was_not_applied(self):
        self.stage("UP000005640_9606.fasta", fasta(50))
        with uniprot(gene_count=50):
            m, _ = self.fetch("--content", "full")
        self.assertTrue(any("--content full was requested" in w for w in m["warnings"]),
                        m["warnings"])

    # -- filename without a taxid ----------------------------------------------------------
    def test_staged_name_without_taxid_uses_uniprot(self):
        self.stage("UP000005640.fasta", fasta(50))
        with uniprot(gene_count=50):
            m, _ = self.fetch()
        self.assertEqual((m["organism"], m["taxid"]), ("Homo sapiens", 9606))
        self.assertEqual(m["organism_source"], "uniprot_api:UP000005640")

    def test_staged_name_without_taxid_offline_uses_curated_accession(self):
        self.stage("UP000005640.fasta", fasta(50))
        with offline():
            m, _ = self.fetch()
        self.assertEqual((m["organism"], m["taxid"]), ("Homo sapiens", 9606))
        self.assertEqual(m["organism_source"], "curated_table_by_accession")

    def test_unknown_staged_proteome_offline_records_nothing_rather_than_guessing(self):
        self.stage("UP000999999.fasta", fasta(50))
        with offline():
            m, _ = self.fetch(proteome="UP000999999")
        self.assertEqual((m["organism"], m["taxid"]), ("", 0))
        self.assertEqual(m["organism_source"], "none")
        self.assertEqual(m["n_proteome"], 50)            # the build itself still succeeded

    # -- downstream: the methods sentence --------------------------------------------------
    def test_methods_text_names_the_copy_not_a_release(self):
        self.stage("UP000005640_9606.fasta", fasta(50))
        with uniprot(gene_count=50):
            _, meta_path = self.fetch()
        md = os.path.join(self.root, "methods.md")
        # Thermo is identified by filename alone, so a bare name is a valid --raw here.
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                            "--raw", os.path.join(self.root, "FL_hela.raw"),
                            "--fasta-meta", meta_path, "--out", md],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(md) as fh:
            text = fh.read()
        self.assertIn("pre-staged copy of the UniProt Homo sapiens reference proteome "
                      "(UP000005640, taxid 9606; copy dated 2025-04-25", text)
        self.assertIn("inferred, not verified: UniProt lists 50 genes", text)
        self.assertNotIn(LIVE_RELEASE, text)


if __name__ == "__main__":
    unittest.main()
