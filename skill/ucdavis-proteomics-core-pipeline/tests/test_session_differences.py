#!/usr/bin/env python3
"""DIFFERENCES.md -- what a re-analysis changed from the session it re-analyses.

PROT_0756's v2 re-ran the DE with --block and nothing else, and DIFFERENCES.md read "identical
settings; difference is data/environment only": the block, the contaminant policy, the database's
sidecar state, the design, the skill and limpa versions were never compared, and a dead
"Validated workflow commit" row was. A setting neither session records is said to be not
compared, never shown as "None" or counted as identical.
stdlib only.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))
sys.path.insert(0, HERE)

import session                          # noqa: E402
import test_deposit_package as tdp      # noqa: E402

UNBLOCKED = {"applied": False, "note": "no --block: samples modelled as independent"}
BLOCKED = {"applied": True, "column": "Mouse", "effect": "random", "scope": "within",
           "consensus_correlation": 0.41}


def write(path, obj):
    tdp.write(path, obj if isinstance(obj, str) else json.dumps(obj))


def pair(d, old_de, new_de):
    """Two sessions over the same raw files, differing only in de_provenance.json."""
    a = tdp.dia_session(os.path.join(d, "a"))
    b = tdp.dia_session(os.path.join(d, "b"))
    shutil.copyfile(a["raw_list"], b["raw_list"])
    for p, de in ((a, old_de), (b, new_de)):
        write(os.path.join(p["de_dir"], "de_provenance.json"), de)
    return a, b


def differences(a, b):
    with open(session.write_differences(a["session_dir"], b["session_dir"]),
              encoding="utf-8") as fh:
        return fh.read()


def rows(md):
    return [ln for ln in md.splitlines() if ln.startswith("| ") and not ln.startswith("| Aspect")]


class Differences(unittest.TestCase):
    def test_a_block_only_change_is_a_difference(self):
        with tempfile.TemporaryDirectory() as d:
            a, b = pair(d, dict(tdp.DE_PROV, block=UNBLOCKED), dict(tdp.DE_PROV, block=BLOCKED))
            md = differences(a, b)
            self.assertIn("| Blocking (column, effect, scope) | none (samples modelled as "
                          "independent) | Mouse (random effect, scope within) |", md)
            self.assertEqual(len(rows(md)), 1, rows(md))
            self.assertNotIn("identical", md)
            self.assertNotIn("none of the settings", md)
            self.assertIn("Same raw files (3, input/raw_files.txt).", md)

    def test_every_setting_the_records_carry_is_compared(self):
        with tempfile.TemporaryDirectory() as d:
            old = dict(tdp.DE_PROV, design="~ 0 + groups", packages={"limpa": "1.4.0"},
                       contaminants={"policy": "kept", "tag": "Cont_"})
            new = dict(tdp.DE_PROV, design="~ 0 + groups + Batch", packages={"limpa": "1.5.0"},
                       contaminants={"policy": "removed", "tag": "Cont_"})
            a, b = pair(d, old, new)
            for p, ver, reg in ((a, "2.7.0", "abc123"), (b, "2.8.0", "def456")):
                write(os.path.join(p["repro_dir"], "run_manifest.json"),
                      {"skill": {"version": ver}, "registry": {"commit": reg}})
            meta = json.load(open(b["fasta_meta"], encoding="utf-8"))
            meta.update(contaminant_target_rule="identity", min_unique_peptides=2)
            write(b["fasta_meta"], meta)
            md = differences(a, b)
            for want in ("| Design (with covariates) | ~ 0 + groups | ~ 0 + groups + Batch |",
                         "| Contaminant policy (DE) | kept (Cont_ entries) | removed (Cont_ "
                         "entries) |",
                         "| FASTA sidecar state (fetch_fasta.sidecar_state) | legacy | current |",
                         "| Skill version | 2.7.0 | 2.8.0 |",
                         "| limpa version | 1.4.0 | 1.5.0 |"):
                self.assertIn(want, md)
            self.assertNotIn("workflow commit", md)
            self.assertNotIn("abc123", md)

    def test_what_neither_session_records_is_said_not_compared(self):
        with tempfile.TemporaryDirectory() as d:
            a, b = pair(d, tdp.DE_PROV, tdp.DE_PROV)
            md = differences(a, b)
            self.assertIn("none of the settings compared here differ", md)
            self.assertIn("Not recorded in either session, so not compared:", md)
            self.assertIn("Blocking (column, effect, scope)", md)
            self.assertIn("Skill version", md)
            self.assertNotIn("None", md)

    def test_one_side_not_recorded_is_a_difference_and_says_so(self):
        with tempfile.TemporaryDirectory() as d:
            a, b = pair(d, tdp.DE_PROV, dict(tdp.DE_PROV, block=BLOCKED))
            md = differences(a, b)
            self.assertIn("| Blocking (column, effect, scope) | not recorded | Mouse (random "
                          "effect, scope within) |", md)

    def test_different_raw_files_are_not_called_the_same(self):
        with tempfile.TemporaryDirectory() as d:
            a, b = pair(d, tdp.DE_PROV, tdp.DE_PROV)
            with open(b["raw_list"], "a", encoding="utf-8") as fh:
                fh.write("/data/extra_run.d\n")
            md = differences(a, b)
            self.assertNotIn("Same raw files", md)
            self.assertIn("| Raw files (input/raw_files.txt) | 3 | 4 (3 in common) |", md)

    def test_the_search_parameters_are_found_in_input_wf(self):
        """The params diff looked only at input/params.* -- the workflow step stages them in
        input/wf/, so it never ran."""
        with tempfile.TemporaryDirectory() as d:
            a, b = pair(d, tdp.DE_PROV, tdp.DE_PROV)
            cfg = os.path.join(b["workflow_dir"], "params.cfg")
            with open(cfg, encoding="utf-8") as fh:
                text = fh.read()
            write(cfg, text + "--mass-acc 12\n")
            md = differences(a, b)
            self.assertIn("## Search-parameter diff", md)
            self.assertIn("+--mass-acc 12", md)


if __name__ == "__main__":
    unittest.main(verbosity=2)
