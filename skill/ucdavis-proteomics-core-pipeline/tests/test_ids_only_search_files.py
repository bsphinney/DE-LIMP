#!/usr/bin/env python3
"""An identifications-only search's results travel with its record (2.14.2 C2, E2E F4 and
check 9).

`session.py which-search --paths` -- the list the laptop pull (SKILL.md "Pull what the report
needs") fetches -- named the search record, `<folder>/report.parquet` and the workflow manifest.
A Sage search for identifications only writes no report.parquet: its results are
identifications.tsv and sage_identifications.json (how the rows were kept). So the E2E run's pull
brought back the record and the manifest, and a laptop-side report writer had nothing to read.

session.search_result_files() is now THE list of a search's results: the report it quantified
into, or the identifications and their record; with engine_tables, the engine's own tables the
record's adapter entries name (what the bundle fingerprints, test_reproduce_replays_what_ran.py).
What these tests pin:
  * which-search --paths for a Sage identifications-only search: record, identifications.tsv,
    sage_identifications.json, manifest -- and the JSON's `pull` says the same;
  * a FragPipe identifications-only search: identifications.tsv and fragpipe_identifications.json;
  * a quantifying search is unchanged (pin): record, report.parquet, manifest;
  * engine_tables finds Sage's results.sage.parquet through the record, also in a copied session
    whose record names the original path.
Stdlib only, no network.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import two_search_fixture as fx                 # noqa: E402
import session                                  # noqa: E402

write = fx.write


def ids_only_session(root, engine="sage", record_root=None):
    """A session whose one search (output/search/) was for identifications only. `record_root`:
    the folder the record names its files under (another, for a copied session). -> (p, folder)."""
    p, raws = fx._layout(root, "ids_only", True)
    folder = p["search_out"]
    named = record_root or folder
    write(os.path.join(folder, "identifications.tsv"),
          "run\tpeptide\tproteins\tcharge\tspectrum_q\tpeptide_q\tprotein_q\n")
    rec = {"engine": engine, "version": "0.14.7" if engine == "sage" else "24.0",
           "files": raws, "quantified": False,
           "result": {"engine": engine, "report": None, "quantified": False}}
    if engine == "sage":
        write(os.path.join(folder, "results.sage.parquet"), "PAR1")
        ids = {"source": os.path.join(named, "results.sage.parquet"),
               "identifications": os.path.join(named, "identifications.tsv"),
               "row_unit": "one row per (run, modified peptide, charge)"}
        write(os.path.join(folder, "sage_identifications.json"), json.dumps(ids))
        rec["sage_identifications"] = ids
        rec["sage_input"] = {"dotd": "direct"}
    else:
        write(os.path.join(folder, "fragpipe_identifications.json"), json.dumps({"rows": 0}))
        rec["quantification"] = {"status": "not_applicable"}
    write(os.path.join(folder, "search_provenance.json"), json.dumps(rec))
    return p, folder


class Tmp(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def which(self, p, *extra):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "which-search",
                            "--dir", p["session_dir"], *extra],
                           capture_output=True, text=True, timeout=120)
        return r.returncode, r.stdout, r.stderr


class PullTheIdentifications(Tmp):
    def test_a_sage_ids_only_search_pulls_its_identifications_and_their_record(self):
        p, _ = ids_only_session(self.tmp)
        rc, out, err = self.which(p, "--paths")
        self.assertEqual(rc, 0, err)
        want = ["output/search/search_provenance.json", "output/search/identifications.tsv",
                "output/search/sage_identifications.json", "input/wf/workflow.manifest.json"]
        self.assertEqual(out.split(), want)
        rc, out, err = self.which(p)
        self.assertEqual(rc, 0, err)
        self.assertEqual(json.loads(out)["pull"], want)

    def test_a_fragpipe_ids_only_search_pulls_its_identifications_and_their_record(self):
        p, _ = ids_only_session(self.tmp, engine="fragpipe")
        rc, out, err = self.which(p, "--paths")
        self.assertEqual(rc, 0, err)
        self.assertEqual(out.split(), ["output/search/search_provenance.json",
                                       "output/search/identifications.tsv",
                                       "output/search/fragpipe_identifications.json",
                                       "input/wf/workflow.manifest.json"])

    def test_a_quantifying_search_pulls_what_it_did(self):
        """(pin) record, the report the DE read, the manifest -- as before."""
        p, _ = fx.one_search(self.tmp)
        rc, out, err = self.which(p, "--paths")
        self.assertEqual(rc, 0, err)
        self.assertEqual(out.split(), ["output/search/search_provenance.json",
                                       "output/search/report.parquet",
                                       "input/wf/workflow.manifest.json"])


class EngineTables(Tmp):
    def test_sages_own_table_is_found_through_the_record(self):
        p, folder = ids_only_session(self.tmp)
        got = session.search_result_files(p, engine_tables=True)
        self.assertEqual([os.path.basename(f) for f in got],
                         ["identifications.tsv", "sage_identifications.json",
                          "results.sage.parquet"])
        self.assertEqual(session.search_result_files(p),
                         [os.path.join(folder, "identifications.tsv"),
                          os.path.join(folder, "sage_identifications.json")],
                         "the pull leaves the engine's own table out")

    def test_a_copied_session_finds_it_in_its_own_folder(self):
        p, folder = ids_only_session(self.tmp, record_root="/elsewhere/output/search")
        got = session.search_result_files(p, engine_tables=True)
        self.assertIn(os.path.join(folder, "results.sage.parquet"), got)
        self.assertFalse([f for f in got if f.startswith("/elsewhere")])


if __name__ == "__main__":
    unittest.main()
