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


def bait_session(root, entries):
    """A session whose FASTA record has the added entries [(accession, name, header)]."""
    s = os.path.join(root, "bait_" + str(len(os.listdir(root))))
    os.makedirs(os.path.join(s, "input"))
    with open(os.path.join(s, et.FASTA_META), "w", encoding="utf-8") as fh:
        json.dump({"added_sequences": [{"file": "/data/bait.fasta", "sha256": "1" * 64,
                                        "n_entries": len(entries), "entries": [
            {"accession": a, "name": n, "header": h, "length": 240, "sequence_sha256": "0" * 64}
            for a, n, h in entries]}]}, fh)
    return s


TURBO = [("1000001", "IN|1000001|TurboID-Ryr2", "IN|1000001|TurboID-Ryr2 fusion GN=Ryr2"),
         ("1000002", "IN|1000002|TAG", "IN|1000002|TAG synthetic tag")]


class BaitAskedOnce(unittest.TestCase):
    """Staff review (2.12): a TurboID study with its fusion added by --add-fasta was asked for the
    bait at step 3 (`set --bait`) and again by `propose` after step 6. Now propose asks nothing
    when the recorded bait names one of the added sequences or set-bait was answered."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def propose(self, s):
        p = run("propose", "--session", s)
        self.assertEqual(p.returncode, 0, p.stderr)
        return json.loads(p.stdout)

    def test_a_step_3_bait_that_names_an_added_sequence_is_not_asked_again(self):
        for told in ("Ryr2", "TurboID-Ryr2", "1000001"):
            with self.subTest(told=told):
                s = bait_session(self.d, TURBO)
                et.record(s, "proximity", "user", bait=told)
                out = self.propose(s)
                self.assertIsNone(out["ask_bait"])
                self.assertEqual(out["bait_recorded"]["accession"], "1000001")
                self.assertEqual(out["bait_recorded"]["bait"], told)

    def test_a_set_bait_answer_is_not_asked_again(self):
        s = bait_session(self.d, TURBO)
        et.record(s, "proximity", "user")
        self.assertIsNotNone(self.propose(s)["ask_bait"])          # not asked yet: asked
        et.record_bait(s, "1000001")
        out = self.propose(s)
        self.assertIsNone(out["ask_bait"])
        self.assertEqual(out["bait_recorded"], {"bait": "1000001",
                                                "source": "set-bait (asked after step 6)"})
        et.record_bait(s, "none")
        self.assertIsNone(self.propose(s)["ask_bait"])

    def test_a_step_3_bait_that_names_no_added_sequence_is_confirmed_not_re_asked_blind(self):
        s = bait_session(self.d, [TURBO[1]])
        et.record(s, "ip", "user", bait="Ryr2")
        out = self.propose(s)
        self.assertTrue(out["ask_bait"].startswith(
            "The bait was recorded at step 3 as Ryr2, which names none of the sequences the user "
            "added. The search database has sequence(s) the user added: 1000002."), out["ask_bait"])
        self.assertIsNone(out["bait_recorded"])

    def test_no_type_recorded_asks_as_before(self):
        out = self.propose(bait_session(self.d, TURBO))
        self.assertIn("is one of them the bait?", out["ask_bait"])

    def test_an_unreadable_record_is_said_and_the_bait_asked(self):
        s = bait_session(self.d, TURBO)
        with open(et.path_for(s), "w") as fh:
            fh.write("{}")
        p = run("propose", "--session", s)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("is not a experiment_type/1 record: the bait is asked", p.stderr)
        self.assertIsNotNone(json.loads(p.stdout)["ask_bait"])


if __name__ == "__main__":
    unittest.main()
