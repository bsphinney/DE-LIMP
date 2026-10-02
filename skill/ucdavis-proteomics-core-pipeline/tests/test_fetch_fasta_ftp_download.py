#!/usr/bin/env python3
"""A one-per-gene download that does not arrive whole must never become the full proteome.

A staff search (2026-10-01, human, skill 2.9.0): the FTP transfer of UP000005640_9606.fasta.gz
was cut short ("Compressed file ended before the end-of-stream marker was reached"), and
fetch_fasta fell back to the REST full proteome -- 147,520 entries for the 20,652 asked for --
with exit 0. A re-run a minute later was whole.

Guards, all offline (the HTTP layer is faked):
  * a body shorter than its Content-Length, a cut gzip stream and a 5xx are retried, and a
    whole file on a later attempt is used, with the attempts recorded in the sidecar;
  * still broken after the retries -> exit != 0, and the REST full proteome is never asked for;
  * only a 404 falls back to the REST full set (loudly, as before);
  * a whole file whose entry count is far from UniProt's geneCount stops too.
"""
import contextlib
import gzip
import http.client
import io
import json
import os
import sys
import tempfile
import unittest
import urllib.error
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import fetch_fasta as ff  # noqa: E402

UP, TAXID = "UP000005640", 9606
FTP_URL = f"{ff.FTP_REF}/Eukaryota/{UP}/{UP}_{TAXID}.fasta.gz"


def fasta(n):
    return "".join(f">sp|P{i:05d}|G{i}_HUMAN Protein {i} OS=Homo sapiens OX=9606 GN=G{i}\n"
                   f"MPEPTIDEKAAAGGG{i}\n" for i in range(n))


def gz(text):
    return gzip.compress(text.encode())


def uniprot(gene_count):
    data = {"id": UP, "taxonomy": {"scientificName": "Homo sapiens", "taxonId": TAXID},
            "proteomeType": "Reference proteome", "proteinCount": gene_count * 7,
            "geneCount": gene_count, "superkingdom": "eukaryota"}
    return mock.patch.object(ff, "_get_json", return_value=(data, {"x-uniprot-release": "2026_03"}))


class FakeResponse(io.BytesIO):
    """What urlopen returns, as far as _download reads it: a body and case-insensitive headers."""

    def __init__(self, body, content_length=None):
        super().__init__(body)
        self.headers = http.client.HTTPMessage()
        if content_length is not None:
            self.headers["Content-Length"] = str(content_length)


class FakeNet:
    """Answers each URL from a script: a list of (body, Content-Length) tuples or HTTP codes,
    one per request, in order. Records every URL asked for."""

    def __init__(self, script):
        self.script = {u: list(v) for u, v in script.items()}
        self.asked = []

    def __call__(self, url, timeout=300):
        self.asked.append(url)
        key = next((u for u in self.script if url.startswith(u)), None)
        if key is None or not self.script[key]:
            raise urllib.error.URLError(f"unexpected request {url}")
        step = self.script[key].pop(0)
        if isinstance(step, int):
            raise urllib.error.HTTPError(url, step, "status", http.client.HTTPMessage(), None)
        body, clen = step
        return FakeResponse(body, clen)


class FtpDownload(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.out = os.path.join(self._td.name, "search.fasta")
        self.slept = []

    def tearDown(self):
        self._td.cleanup()

    def fetch(self, net, gene_count=40, *extra):
        argv = ["fetch_fasta.py", "fetch", "--proteome", UP, "--contaminants", "none",
                "--out", self.out, *extra]
        stdout, stderr = io.StringIO(), io.StringIO()
        with uniprot(gene_count), mock.patch.object(ff, "_open", net), \
                mock.patch.object(ff, "_sleep", self.slept.append), \
                mock.patch.object(sys, "argv", argv), \
                contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            try:
                rc = ff.main()
            except SystemExit as e:
                rc = e.code
        return rc, stdout.getvalue(), stderr.getvalue()

    def sidecar(self):
        with open(self.out + ".meta.json") as fh:
            return json.load(fh)

    def test_a_short_body_is_retried_and_the_whole_file_used(self):
        whole = gz(fasta(40))
        net = FakeNet({FTP_URL: [(whole[: len(whole) // 2], len(whole)), (whole, len(whole))]})
        rc, _out, err = self.fetch(net)
        self.assertIn(rc, (0, None), err)
        m = self.sidecar()
        self.assertEqual(m["content_used"], "one_per_gene")
        self.assertEqual(m["source"], f"uniprot_ftp:{UP}")
        self.assertEqual(m["n_proteome"], 40)
        chk = m["download_check"]
        self.assertEqual(chk["attempts"], 2)
        self.assertEqual(chk["bytes"], len(whole))
        self.assertEqual(chk["content_length"], len(whole))
        self.assertEqual(chk["n_entries"], 40)
        self.assertEqual(chk["uniprot_gene_count"], 40)
        self.assertIn("IncompleteDownload", chk["failed_attempts"][0])
        self.assertEqual(self.slept, [ff.FTP_RETRY_DELAYS[0]])
        self.assertIn("retrying", err)

    def test_a_cut_gzip_stream_with_no_content_length_never_becomes_the_full_proteome(self):
        """The reported failure: the gzip ends early. Every attempt cut -> exit != 0, and the
        REST full proteome is not even asked for."""
        whole = gz(fasta(40))
        cut = (whole[: len(whole) - 30], None)
        net = FakeNet({FTP_URL: [cut] * (1 + len(ff.FTP_RETRY_DELAYS)),
                       ff.UNIPROT_REST + "/uniprotkb/stream": [(fasta(280).encode(), None)]})
        rc, _out, err = self.fetch(net)
        self.assertNotIn(rc, (0, None))
        self.assertIn("NOT falling back to the REST full proteome", str(rc))
        self.assertIn("EOFError", str(rc))
        self.assertIn("--content full", str(rc))
        self.assertFalse(any("/uniprotkb/stream" in u for u in net.asked), net.asked)
        self.assertEqual(len([u for u in net.asked if u == FTP_URL]), 1 + len(ff.FTP_RETRY_DELAYS))
        self.assertEqual(self.slept, list(ff.FTP_RETRY_DELAYS))
        self.assertFalse(os.path.exists(self.out), "nothing may be written on a failed build")
        self.assertFalse(os.path.exists(os.path.join(self._td.name, f"{UP}.fasta.gz")),
                         "no partial download may be left beside the output")

    def test_a_server_error_is_retried(self):
        whole = gz(fasta(40))
        net = FakeNet({FTP_URL: [503, (whole, len(whole))]})
        rc, _out, err = self.fetch(net)
        self.assertIn(rc, (0, None), err)
        self.assertEqual(self.sidecar()["download_check"]["failed_attempts"], ["HTTP 503"])

    def test_only_a_404_falls_back_to_the_rest_full_set_loudly(self):
        net = FakeNet({FTP_URL: [404],
                       ff.UNIPROT_REST + "/uniprotkb/stream": [(fasta(280).encode(), None)]})
        rc, _out, err = self.fetch(net)
        self.assertIn(rc, (0, None), err)
        m = self.sidecar()
        self.assertEqual(m["content_used"], "full")
        self.assertEqual(m["source"], f"uniprot_rest:{UP}")
        self.assertEqual(m["download_check"]["fallback"], "uniprot_rest full")
        self.assertIn("404", m["download_check"]["one_per_gene_ftp"])
        self.assertTrue(any("Falling back to the REST full proteome" in w for w in m["warnings"]))
        self.assertEqual(self.slept, [], "a 404 is an answer, not a transient failure")

    def test_a_whole_file_with_the_wrong_entry_count_stops(self):
        whole = gz(fasta(30))
        net = FakeNet({FTP_URL: [(whole, len(whole))]})
        rc, _out, _err = self.fetch(net, 40)
        self.assertNotIn(rc, (0, None))
        self.assertIn("geneCount", str(rc))
        self.assertEqual(len(net.asked), 1, "a count that does not fit is not retried")

    def test_an_unmappable_superkingdom_stops_instead_of_guessing(self):
        meta = {"superkingdom": "", "taxid": TAXID, "gene_count": 40}
        with self.assertRaises(RuntimeError) as cm:
            ff.download_ftp_one_per_gene(UP, meta, self._td.name)
        self.assertIn("could not be looked for", str(cm.exception))


class ContentLength(unittest.TestCase):
    def test_every_download_checks_the_declared_length(self):
        with tempfile.TemporaryDirectory() as td:
            dest = os.path.join(td, "x")
            with mock.patch.object(ff, "_open", lambda u, t=0: FakeResponse(b"abc", 5)):
                with self.assertRaises(ff.IncompleteDownload):
                    ff._download("https://example.invalid/x", dest)
            with mock.patch.object(ff, "_open", lambda u, t=0: FakeResponse(b"abcde", 5)):
                self.assertEqual(ff._download("https://example.invalid/x", dest), (5, 5))
            with mock.patch.object(ff, "_open", lambda u, t=0: FakeResponse(b"abc")):
                self.assertEqual(ff._download("https://example.invalid/x", dest), (3, None))


if __name__ == "__main__":
    unittest.main()
