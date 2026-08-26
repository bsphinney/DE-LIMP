#!/usr/bin/env python3
"""The search DATABASE must reach FRAN's manifest.

FRAN's delimp_searches has carried fasta_path / fasta_md5 / fasta_n_proteins since its schema
was written, and as of 2026-08-25 they were populated for 157 / 0 / 0 of 2,014 searches --
corpus_ingest.py never wrote them. FRAN can now parse them out of an engine log as a fallback,
but a search that came through this skill knows its database exactly, so it should hand it over
rather than making the corpus guess.

The failure this guards is silent: a manifest with no database still stages, still ingests, and
leaves the corpus unable to say whether two searches used comparable proteomes. Entries-per-gene
is what separates a real depth difference from database redundancy.
"""
import hashlib
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fran_deposit as fd  # noqa: E402
from fetch_fasta import _count_entries  # noqa: E402  -- the ONE definition; fran_deposit imports it

# 3 entries; the middle '>' is placed so a naive per-chunk count would still see it, while
# test_entry_count_survives_chunk_boundary below forces the boundary case explicitly.
FASTA = ">sp|P1|A\nPEPTIDEK\n>sp|P2|B\nPEPTIDER\n>sp|P3|C\nPEPTIDEK\n"


class FastaInManifest(unittest.TestCase):
    def _meta(self, root, **extra):
        fa = os.path.join(root, "db.fasta")
        with open(fa, "w") as fh:
            fh.write(FASTA)
        meta = {"fasta": fa, "organism": "Homo sapiens", "taxid": 9606, **extra}
        mp = os.path.join(root, "db.fasta.meta.json")
        with open(mp, "w") as fh:
            json.dump(meta, fh)
        return fa, mp

    def test_reads_path_md5_and_count_from_meta(self):
        with tempfile.TemporaryDirectory() as root:
            fa, mp = self._meta(root, md5="deadbeef", n_entries=3)
            path, md5, n = fd.fasta_from_meta(root, mp)
            self.assertEqual(path, fa)
            self.assertEqual(md5, "deadbeef")   # trusted, not recomputed
            self.assertEqual(n, 3)

    def test_recomputes_when_meta_predates_the_fields(self):
        """Older meta.json files have `fasta` but no md5/n_entries. Those must be filled from
        the file rather than left NULL, or every search staged before this change stays blind."""
        with tempfile.TemporaryDirectory() as root:
            fa, mp = self._meta(root)
            path, md5, n = fd.fasta_from_meta(root, mp)
            self.assertEqual(path, fa)
            self.assertEqual(md5, hashlib.md5(FASTA.encode()).hexdigest())  # noqa: S324
            self.assertEqual(n, 3)

    def test_the_count_already_in_older_sidecars_is_used_rather_than_rescanned(self):
        """`fetch_fasta.py` has ALWAYS written `n_sequences` (proteome + appended contaminants),
        and it is the same number as the header count -- verified against a real sidecar,
        n_sequences 34,306 against 34,306 headers. Reading it keeps `check` a stat rather than a
        parse: no sidecar predating the md5/n_entries fields has to be re-scanned, which on a
        login node costs 0.6-1.3 s per hundred MB of FASTA."""
        with tempfile.TemporaryDirectory() as root:
            fa, mp = self._meta(root)
            with open(mp) as fh:
                m = json.load(fh)
            m["n_sequences"] = 4321               # deliberately NOT the real header count
            with open(mp, "w") as fh:
                json.dump(m, fh)
            _, _, n = fd.fasta_from_meta(root, mp)
            self.assertEqual(n, 4321, "the count already in the sidecar was ignored and the "
                                      "FASTA re-scanned instead")

    def test_the_helpers_are_imported_not_copied(self):
        """One definition. The chunk-boundary counter is exactly the logic you do not want two
        copies of -- a bug fixed in one would silently survive in the other."""
        src = open(os.path.join(SCRIPTS, "fran_deposit.py")).read()
        self.assertNotIn("def _count_entries(", src)
        self.assertIn("from fetch_fasta import", src)

    def test_absent_not_guessed_when_no_meta(self):
        with tempfile.TemporaryDirectory() as root:
            self.assertEqual(fd.fasta_from_meta(root, None), (None, None, None))

    def test_unreadable_fasta_leaves_md5_null_rather_than_inventing(self):
        with tempfile.TemporaryDirectory() as root:
            _, mp = self._meta(root)
            os.remove(os.path.join(root, "db.fasta"))
            path, md5, n = fd.fasta_from_meta(root, mp)
            self.assertTrue(path)          # the path is still what the search used
            self.assertIsNone(md5)         # but nothing is fabricated about its contents
            self.assertIsNone(n)

    def test_entry_count_survives_chunk_boundary(self):
        """A per-chunk count(b"\\n>") drops any header landing exactly on a 1 MiB boundary --
        invisible on small files, wrong on real proteomes."""
        with tempfile.TemporaryDirectory() as root:
            fa = os.path.join(root, "big.fasta")
            chunk = 1 << 20
            with open(fa, "w") as fh:
                head = ">a\n"
                fh.write(head + "P" * (chunk - len(head) - 1) + "\n")  # next '>' starts at 1 MiB
                fh.write(">b\nPEPTIDEK\n>c\nPEPTIDER\n")
            self.assertEqual(_count_entries(fa), 3)


if __name__ == "__main__":
    unittest.main(verbosity=2)
