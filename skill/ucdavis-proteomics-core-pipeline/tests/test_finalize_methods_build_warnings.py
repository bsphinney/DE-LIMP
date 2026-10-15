#!/usr/bin/env python3
"""When finalize builds the Methods, the database's build warnings reach the agent, not the zip.

Since 2.12 make_methods.py keeps fetch_fasta.py's build warnings out of methods.md (a client's
Methods once carried them with an /nfs path) and lists them under database_build_warnings in the
JSON it prints. session.py finalize builds the Methods through make_deposit.ensure_methods when
the session has none, and that JSON was discarded -- so on that route the warnings went nowhere.
Now make_deposit prints them on stderr and hands them to finalize, which adds them to the JSON it
prints; MANIFEST.txt, which travels in the zip, gets only their count.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from test_deposit_package import dia_session, finalize, manifest_lines, write  # noqa: E402

WARNINGS = ["one-per-gene FTP file not found for UP000005640: used the REST full proteome from "
            "/quobyte/proteomics-grp/MRS/UP000005640_9606.fasta instead",
            "added sequence A0A0000003 (/nfs/lssc0/flinders/Data/lab/service/PROT_0000/x.fasta) "
            "is identical to database entry Q00003"]


class FinalizeSurfacesBuildWarnings(unittest.TestCase):
    def run_finalize(self, warnings):
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        p = dia_session(tmp)
        with open(p["fasta_meta"]) as fh:
            meta = json.load(fh)
        write(p["fasta_meta"], json.dumps(dict(meta, warnings=warnings)))
        r = finalize(p["session_dir"], "--zip")
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        return p, r, json.loads(r.stdout)

    def test_the_agent_gets_them_and_the_zip_does_not(self):
        p, r, res = self.run_finalize(WARNINGS)
        self.assertEqual(res["database_build_warnings"], WARNINGS)
        for w in WARNINGS:
            self.assertIn(f"[make_deposit] database build warning (not in methods.md; tell the "
                          f"user): {w}", r.stderr)
        (line,) = [ln for ln in manifest_lines(p)
                   if "Publication methods (output/methods.md)" in ln]
        self.assertIn("2 database build warnings left out of the Methods, for you to tell the user",
                      line)
        with open(p["manifest_txt"]) as fh:
            manifest = fh.read()
        with open(res["methods"], encoding="utf-8") as fh:
            methods = fh.read()
        for text in (manifest, methods):
            self.assertNotIn("/quobyte", text)
            self.assertNotIn("/nfs", text)
            self.assertNotIn("identical to database entry", text)

    def test_none_no_key(self):
        p, r, res = self.run_finalize([])
        self.assertNotIn("database_build_warnings", res)
        self.assertNotIn("database build warning", r.stderr)


if __name__ == "__main__":
    unittest.main()
