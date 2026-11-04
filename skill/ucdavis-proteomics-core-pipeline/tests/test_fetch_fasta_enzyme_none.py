#!/usr/bin/env python3
"""
`fetch_fasta.py fetch --enzyme none`: a sample with no digestion enzyme (2.14.2, L11 N3).

A non-specific search of peptides that were never digested (a hydrolysate, a peptidome) was built
with the default --enzyme trypsin,lysc, because "no enzyme" could not be said: parse_enzymes
refused an empty or "none" value. So the bovine trypsin contaminant entry stayed Cont_ -- kept as
a digestion reagent -- although it is identical to a target protein of the bovine database, and
the search had to be cancelled and rebuilt with --drop-contaminant by hand.

With --enzyme none no enzyme was added, so no protease is a digestion contaminant: every protease
entry identical to a target is dropped as a target, like any other protein. The sidecar records
`enzyme: none` (and digestion_enzymes_used []), and every replay of the build says --enzyme none.
"""
import argparse
import contextlib
import io
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import fetch_fasta as ff  # noqa: E402
from test_fetch_fasta_contaminant_overlap import CONT, PIG_TARGET, _organism, headers  # noqa: E402


def build(target, cont, *args):
    with tempfile.TemporaryDirectory() as root:
        tgt, cp = os.path.join(root, "target.fasta"), os.path.join(root, "cont.fasta")
        with open(tgt, "w") as fh:
            fh.write(target)
        with open(cp, "w") as fh:
            fh.write(cont)
        out = os.path.join(root, "out", "search.fasta")
        argv = ["fetch_fasta.py", "fetch", "--path", tgt, "--contaminants", "universal",
                "--contaminants-path", cp, *args]
        with mock.patch.object(sys, "argv", argv + _organism(argv) + ["--out", out]), \
                contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            assert ff.main() == 0
        with open(out + ".meta.json") as fh, open(out) as fa:
            return json.load(fh), fa.read()


class EnzymeNone(unittest.TestCase):
    def test_none_is_no_enzyme(self):
        self.assertEqual(ff.parse_enzymes("none"), ())
        self.assertEqual(ff.parse_enzymes(" None "), ())
        with self.assertRaises(argparse.ArgumentTypeError):
            ff.parse_enzymes("none,trypsin")                # none is alone or not at all
        with self.assertRaises(argparse.ArgumentTypeError):
            ff.parse_enzymes("")                            # still not a way to say it

    def test_the_trypsin_entry_identical_to_a_target_is_dropped(self):
        m, fasta = build(PIG_TARGET, CONT, "--enzyme", "none")
        self.assertEqual(m["enzyme"], "none")
        self.assertEqual(m["digestion_enzymes_used"], [])
        self.assertEqual(m["contaminants_kept_despite_target_match"], [])
        drop = {r["cont_acc"]: r for r in m["contaminants_dropped_as_target"]}
        self.assertEqual(drop["Cont_P00761"]["target_acc"], "P00761")
        self.assertEqual(drop["Cont_P00761"]["enzyme_not_used"], "trypsin")
        self.assertNotIn(">sp|Cont_P00761|TRYP_PIG", headers(fasta))
        self.assertIn(">sp|P00761|TRYP_PIG", headers(fasta))         # the target stays
        self.assertFalse([w for w in m["warnings"] if w.startswith("kept 1 digestion-enzyme")])

    def test_the_default_is_unchanged(self):
        """(pin) no --enzyme: trypsin,lysc, and porcine trypsin stays a contaminant"""
        m, fasta = build(PIG_TARGET, CONT)
        self.assertEqual(m["digestion_enzymes_used"], ["lysc", "trypsin"])
        self.assertEqual(m["enzyme"], "lysc,trypsin")
        self.assertEqual([r["cont_acc"] for r in m["contaminants_kept_despite_target_match"]],
                         ["Cont_P00761"])
        self.assertIn(">sp|Cont_P00761|TRYP_PIG", headers(fasta))

    def test_the_recorded_none_is_read_back_as_none(self):
        """A sidecar's [] is the build's 'none', not a sidecar from before --enzyme (default)."""
        self.assertEqual(ff.enzymes_recorded({"digestion_enzymes_used": []}), ())
        self.assertEqual(ff.enzyme_arg({"digestion_enzymes_used": []}), "none")
        self.assertEqual(ff.enzymes_recorded({}), ff.parse_enzymes(ff.DEFAULT_ENZYMES))
        self.assertEqual(ff.enzyme_arg({}), ff.DEFAULT_ENZYMES)     # a sidecar from before
        self.assertEqual(ff.enzyme_arg({"digestion_enzymes_used": ["gluc"]}), "gluc")

    def test_the_replay_says_enzyme_none(self):
        """reproduce.sh rebuilds the database with --enzyme none, not the default"""
        import subprocess
        from test_provenance_replay import build_bundle
        with tempfile.TemporaryDirectory() as tmp:
            out, wf = build_bundle(tmp)
            out2 = os.path.join(tmp, "repro2")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                                "--outdir", out2, "--workflow-manifest",
                                os.path.join(wf, "workflow.manifest.json"),
                                "--engine", "diann", "--acquisition", "DIA",
                                "--instrument", "timsTOF HT", "--organism-taxid", "10090",
                                "--fasta-info", json.dumps(
                                    {"proteome": "UP000000589", "content_used": "one_per_gene",
                                     "contaminant_set": "universal",
                                     "digestion_enzymes_used": [], "enzyme": "none"})],
                               capture_output=True, text=True, timeout=300)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(out2, "reproduce.sh")) as fh:
                sh = fh.read()
        self.assertIn("--enzyme none", sh)
        self.assertNotIn("--enzyme trypsin,lysc", sh)

    def test_the_keratin_refusal_s_rebuild_command_keeps_none(self):
        """run_search's keratin refusal prints the fetch command to rebuild the database: it
        keeps the build's --enzyme none"""
        import run_search
        with tempfile.TemporaryDirectory() as d:
            fasta = os.path.join(d, "search.fasta")
            with open(fasta, "w") as fh:
                fh.write(">sp|P1|X_BOVIN X\nPEPTIDEK\n"
                         ">sp|Cont_P35908|K22E_HUMAN Keratin, type II cytoskeletal 2 epidermal "
                         "GN=KRT2\nSGGGFSSGSAGIINYQR\n")
            with open(fasta + ".meta.json", "w") as fh:
                json.dump({"proteome": "UP000009136", "contaminant_set": "universal",
                           "digestion_enzymes_used": [], "enzyme": "none",
                           "keratin_sample": True, "keratin_sample_source": "user"}, fh)
            with self.assertRaises(SystemExit) as cm:
                run_search.keratin_sample_check(fasta, True, engine="sage")
        self.assertIn("--enzyme none", str(cm.exception.code))
        self.assertNotIn("--enzyme trypsin,lysc", str(cm.exception.code))


class TheDigestIsSaid(unittest.TestCase):
    """Follow-up 3 (H7 bug 5): a non-specific search was compared with the tryptic contaminant
    digest without a warning (cleave_at "" dropped out of the comparison), and an --enzyme none
    build called its tryptic stand-in "the search's digest"."""

    def test_a_non_specific_sage_search_is_warned(self):
        import run_search
        with tempfile.TemporaryDirectory() as d:
            fasta = os.path.join(d, "search.fasta")
            with open(fasta, "w") as fh:
                fh.write(">sp|P1|X\nPEPTIDEK\n")
            with open(fasta + ".meta.json", "w") as fh:
                json.dump({"contaminant_digest": dict(ff.DIANN_DIGEST),
                           "min_unique_peptides": 2}, fh)
            cfg = os.path.join(d, "sage.json")
            with open(cfg, "w") as fh:
                json.dump({"database": {"enzyme": {"cleave_at": "", "missed_cleavages": 2,
                                                   "restrict": "P"}}}, fh)
            with contextlib.redirect_stderr(io.StringIO()):
                msg = run_search.warn_contaminant_digest("sage", cfg, fasta)
        self.assertTrue(msg.startswith("[run_search] WARNING"), msg)
        self.assertIn("non-specific", msg)

    def test_an_enzyme_none_build_says_its_rule_used_a_tryptic_stand_in(self):
        m, _ = build(PIG_TARGET, CONT, "--enzyme", "none")
        self.assertIn("tryptic stand-in", m["contaminant_target_rule"])
        self.assertNotIn("in the search's digest", m["contaminant_target_rule"])
        self.assertIn("tryptic stand-in", m["contaminant_digest_note"])
        self.assertTrue([w for w in m["warnings"] if "tryptic stand-in" in w])

    def test_a_build_with_an_enzyme_is_unchanged(self):
        """(pin)"""
        m, _ = build(PIG_TARGET, CONT)
        self.assertIn("in the search's digest", m["contaminant_target_rule"])
        self.assertNotIn("contaminant_digest_note", m)


if __name__ == "__main__":
    unittest.main()
