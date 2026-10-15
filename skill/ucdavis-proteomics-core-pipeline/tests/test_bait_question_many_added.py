#!/usr/bin/env python3
"""The bait question for a whole added proteome counts the sequences and asks for the bait by name.

experiment_type.py propose lists every sequence added with fetch_fasta.py --add-fasta as a bait
candidate and asks "is one of them the bait?" -- right for an EGFP bait or a tag. A staff search
(2026-10-10) added a second organism's UniProt proteome, 4,506 entries: the question would have
listed 4,506 accessions, and set-bait's refusal would have repeated them. Now, above
fetch_fasta.FEW_ADDED_SEQUENCES (10, the same limit the Methods use), propose gives the count and
the file, lists no candidate, and asks for the bait by name; set-bait then takes the bait's
accession, entry name or gene, recording the accession it names and what the user said. At 10 or
fewer nothing changes. Stdlib only.
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

import experiment_type as et  # noqa: E402
import fetch_fasta as ff  # noqa: E402
import make_methods as mm  # noqa: E402

ETYPE = os.path.join(SCRIPTS, "experiment_type.py")


def session(root, n, genes=None):
    """A session whose FASTA record has n added sequences (fetch_fasta.py --add-fasta's shape)."""
    s = os.path.join(root, f"s{n}")
    os.makedirs(os.path.join(s, "input"))
    genes = genes or {}
    entries = [{"accession": f"A0A{i:07d}", "name": f"tr|A0A{i:07d}|P{i}_METSP",
                "header": f"tr|A0A{i:07d}|P{i}_METSP protein {i} OS=Methylocystis sp. OX=0000"
                          + (f" GN={genes[i]}" if i in genes else "") + " PE=4 SV=1",
                "length": 300, "sequence_sha256": f"{i:064x}"} for i in range(n)]
    with open(os.path.join(s, et.FASTA_META), "w", encoding="utf-8") as fh:
        json.dump({"added_sequences": [{"file": "/data/PROT_0000/second_organism.fasta",
                                        "sha256": "1" * 64, "n_entries": n, "entries": entries}]}, fh)
    return s


def run(*args):
    return subprocess.run([sys.executable, ETYPE, *args], capture_output=True, text=True,
                          timeout=60)


class BaitQuestion(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_one_limit_for_the_methods_and_the_bait_question(self):
        self.assertEqual(ff.FEW_ADDED_SEQUENCES, 10)
        self.assertEqual(mm.ADDED_SEQUENCES_NAMED_MAX, ff.FEW_ADDED_SEQUENCES)

    def test_a_whole_added_proteome_is_counted_not_listed(self):
        s = session(self.d, 4506)
        p = run("propose", "--session", s)
        self.assertEqual(p.returncode, 0, p.stderr)
        out = json.loads(p.stdout)
        self.assertEqual(out["bait_candidates"], [])
        self.assertEqual(out["n_added_sequences"], 4506)
        self.assertEqual(out["ask_bait"],
                         "The search database has 4,506 sequences the user added "
                         "(second_organism.fasta): a whole proteome, too many to list. If this is "
                         "an enrichment, ask the user which protein is the bait -- by gene, entry "
                         "name or accession -- and record it with set-bait --bait <the bait> (or "
                         "--bait none).")
        self.assertLess(len(p.stdout), 4000)

    def test_ten_are_still_listed(self):
        out = json.loads(run("propose", "--session", session(self.d, 10)).stdout)
        self.assertEqual(len(out["bait_candidates"]), 10)
        self.assertIn("A0A0000009", out["ask_bait"])
        self.assertIn("set-bait --bait <accession> (or --bait none)", out["ask_bait"])
        out = json.loads(run("propose", "--session", session(self.d, 11)).stdout)
        self.assertEqual(out["bait_candidates"], [])
        self.assertIn("11 sequences the user added", out["ask_bait"])

    def test_set_bait_takes_a_name_in_a_whole_proteome(self):
        s = session(self.d, 4506, genes={17: "pmoA", 18: "mmoX", 19: "mmoX"})
        et.record(s, "ip", "user")
        for given in ("pmoA", "PMOA", "P17_METSP", "A0A0000017"):
            with self.subTest(given=given):
                p = run("set-bait", "--session", s, "--bait", given)
                self.assertEqual(p.returncode, 0, p.stderr)
                r = et.load(s)
                self.assertEqual((r["bait"], r["bait_confirmed"]["answer"]), ("A0A0000017", "A0A0000017"))
                self.assertEqual(r["bait_confirmed"].get("given_as"),
                                 None if given == "A0A0000017" else given)
        p = run("set-bait", "--session", s, "--bait", "mmoX")          # two entries carry it
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("names 2 of the added sequences (A0A0000018 tr|A0A0000018|P18_METSP, "
                      "A0A0000019 tr|A0A0000019|P19_METSP); give the bait's accession", p.stderr)
        p = run("set-bait", "--session", s, "--bait", "EGFP")
        self.assertIn("no accession, entry name or gene among the 4,506 added sequences "
                      "(second_organism.fasta) is EGFP", p.stderr)
        self.assertLess(len(p.stderr), 1000)
        self.assertEqual(run("set-bait", "--session", s, "--bait", "none").returncode, 0)
        self.assertIsNone(et.load(s)["bait"])

    def test_at_ten_or_fewer_the_accession_is_required_as_before(self):
        s = session(self.d, 2, genes={1: "pmoA"})
        et.record(s, "ip", "user")
        p = run("set-bait", "--session", s, "--bait", "pmoA")
        self.assertIn("not one of the added sequences (A0A0000000, A0A0000001)", p.stderr)
        self.assertEqual(run("set-bait", "--session", s, "--bait", "A0A0000001").returncode, 0)

    def test_the_data_check_says_how_many(self):
        bait, src, cands = et.bait_for_check(session(self.d, 4506))
        self.assertIsNone(bait)
        self.assertEqual(len(cands), 4506)
        self.assertIn("the user added 4,506 sequences (--add-fasta, second_organism.fasta)", src)


if __name__ == "__main__":
    unittest.main()
