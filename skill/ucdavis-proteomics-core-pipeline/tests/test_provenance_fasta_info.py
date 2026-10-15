#!/usr/bin/env python3
"""provenance.py reads the FASTA sidecar from a FILE, checks its shape, and keeps the bundle small.

Staff report (2026-10-09, an NCBI proteome plus a second organism's UniProt proteome added with
--add-fasta): search.fasta.meta.json was 2 MB, 1.69 MB of it the 8,685 per-entry records of
added_sequences[0].entries. `--fasta-info "$(cat search.fasta.meta.json)"` failed with "Argument
list too long" (exit 126). A hand-trimmed copy with added_sequences as a dict then crashed
provenance.py (`'str' object has no attribute 'get'`), after the environment capture had started.

Now: --fasta-info-file <path> and --fasta-info @<path> read the sidecar from the file (the JSON
text still works); a sidecar that is not JSON, not an object, or whose added_sequences is not a
list of objects stops with exit 2 and a message naming the problem, before anything is written;
and the bundle records each added file by path, sha256 and entry count, never the entries.
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

import provenance  # noqa: E402

N_ADDED = 9000


def sidecar(added_file):
    """A sidecar shaped as fetch_fasta.py writes it, with a big --add-fasta file."""
    entries = [{"accession": f"A0A{i:07d}", "name": f"tr|A0A{i:07d}|X{i}_METOB",
                "header": f"tr|A0A{i:07d}|X{i}_METOB Uncharacterized protein OS=Methylocystis "
                          f"sp. OX=0000 GN=g{i} PE=4 SV=1", "length": 300 + i % 200,
                "sequence_sha256": f"{i:064x}"} for i in range(N_ADDED)]
    return {"fasta": "/data/search.fasta", "proteome": None, "content_used": "unknown",
            "organism": "Methylosinus sp.", "taxid": 0, "contaminant_set": "universal",
            "n_contaminants_appended": 381, "digestion_enzymes_used": ["trypsin", "lysc"],
            "n_added_sequences": N_ADDED,
            "added_sequences": [{"file": added_file, "sha256": "ab" * 32,
                                 "n_entries": N_ADDED, "entries": entries}]}


class FastaInfoFromAFile(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.added = os.path.join(self.root, "second_organism.fasta")
        self.meta = os.path.join(self.root, "search.fasta.meta.json")
        with open(self.meta, "w") as fh:
            json.dump(sidecar(self.added), fh)
        self.wm = os.path.join(self.root, "workflow.manifest.json")
        with open(self.wm, "w") as fh:
            fh.write("{}")
        self.out = os.path.join(self.root, "repro")

    def tearDown(self):
        self._td.cleanup()

    def run_prov(self, *fasta_args):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                               "--outdir", self.out, "--workflow-manifest", self.wm,
                               "--setup-json", os.path.join(self.root, "no_setup.json"),
                               "--engine", "diann", "--acquisition", "DIA",
                               "--instrument", "timsTOF HT", *fasta_args],
                              capture_output=True, text=True, cwd=self.root, timeout=300)

    def bundle(self):
        with open(os.path.join(self.out, "run_manifest.json")) as fh:
            man = json.load(fh)
        with open(os.path.join(self.out, "reproduce.sh")) as fh:
            return man, fh.read()

    def test_the_sidecar_is_big_enough_to_have_failed_on_the_command_line(self):
        self.assertGreater(os.path.getsize(self.meta), 2 * 1024 ** 2)

    def test_file_form_records_the_added_file_without_its_entries(self):
        for args in (["--fasta-info-file", self.meta], ["--fasta-info", "@" + self.meta]):
            with self.subTest(args=args[0]):
                r = self.run_prov(*args)
                self.assertEqual(r.returncode, 0, r.stderr[-2000:])
                man, sh = self.bundle()
                self.assertEqual(man["inputs"]["fasta_info"]["added_sequences"],
                                 [{"file": self.added, "sha256": "ab" * 32, "n_entries": N_ADDED}])
                self.assertEqual(man["inputs"]["fasta_info"]["n_added_sequences"], N_ADDED)
                self.assertLess(os.path.getsize(os.path.join(self.out, "run_manifest.json")),
                                100 * 1024)
                self.assertIn(f"--add-fasta {self.added}", sh)
                self.assertIn(f"({N_ADDED} entries, sha256 {'ab' * 6}...)", sh)

    def test_the_json_text_still_works(self):
        small = dict(sidecar(self.added))
        small["added_sequences"] = [dict(small["added_sequences"][0], entries=[])]
        r = self.run_prov("--fasta-info", json.dumps(small))
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        man, sh = self.bundle()
        self.assertEqual(man["inputs"]["fasta_info"]["added_sequences"][0]["n_entries"], N_ADDED)
        self.assertIn(f"--add-fasta {self.added}", sh)

    def refused(self, *args):
        r = self.run_prov(*args)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertFalse(os.path.exists(self.out), "nothing may be written before the check")
        return r.stderr

    def test_added_sequences_as_a_dict_is_a_clear_error(self):
        m = sidecar(self.added)
        m["added_sequences"] = {"file": self.added, "sha256": "ab" * 32, "n_entries": N_ADDED}
        with open(self.meta, "w") as fh:
            json.dump(m, fh)
        err = self.refused("--fasta-info-file", self.meta)
        self.assertIn("added_sequences must be a list of objects", err)
        self.assertNotIn("Traceback", err)

    def test_a_list_of_strings_is_refused_too(self):
        m = sidecar(self.added)
        m["added_sequences"] = [self.added]
        with open(self.meta, "w") as fh:
            json.dump(m, fh)
        self.assertIn("added_sequences must be a list of objects",
                      self.refused("--fasta-info-file", self.meta))

    def test_not_json_not_an_object_or_missing(self):
        with open(self.meta, "w") as fh:
            fh.write("{not json")
        self.assertIn("not valid JSON", self.refused("--fasta-info-file", self.meta))
        with open(self.meta, "w") as fh:
            json.dump([1, 2], fh)
        self.assertIn("not fetch_fasta.py's sidecar", self.refused("--fasta-info", "@" + self.meta))
        self.assertIn("cannot read", self.refused("--fasta-info-file",
                                                  os.path.join(self.root, "nope.json")))

    def test_both_forms_at_once_is_refused(self):
        self.assertIn("give the FASTA sidecar once",
                      self.refused("--fasta-info-file", self.meta, "--fasta-info", "@" + self.meta))


class Helpers(unittest.TestCase):
    def test_trim_keeps_everything_but_the_entries(self):
        m = sidecar("/x/added.fasta")
        t = provenance.fasta_info_for_bundle(m)
        self.assertEqual(t["added_sequences"], [{"file": "/x/added.fasta", "sha256": "ab" * 32,
                                                 "n_entries": N_ADDED}])
        self.assertEqual({k: v for k, v in t.items() if k != "added_sequences"},
                         {k: v for k, v in m.items() if k != "added_sequences"})
        self.assertIn("entries", m["added_sequences"][0])          # the sidecar is not changed
        self.assertEqual(provenance.fasta_info_for_bundle({"proteome": "UP1"}), {"proteome": "UP1"})


if __name__ == "__main__":
    unittest.main()
