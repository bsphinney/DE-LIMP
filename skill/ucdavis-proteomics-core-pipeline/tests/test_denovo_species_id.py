#!/usr/bin/env python3
"""
The de novo -> homology species-ID workflow.

Guards the settings that fail SILENTLY rather than loudly. Every number in the
assertions below came from measurement on the ocelot cohort + HeLa entrapment;
the failure modes are all "runs fine, answer is wrong".
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

from denovo_homology_fdr import (  # noqa: E402
    rate_normalised_accept, FEATURE_SETS, DEFAULT_FEATURES, read_hits,
)


def run(script, *args):
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args],
                       capture_output=True, text=True)
    return p


class TestSearchSettings(unittest.TestCase):
    def _emit(self, *extra):
        p = run("denovo_search_cmd.py", "--query", "q.fa", "--db", "nr",
                "--out", "o.tsv", *extra)
        self.assertEqual(p.returncode, 0, p.stderr)
        return json.loads(p.stdout)

    def test_short_peptide_flags_are_present_by_default(self):
        # Stock DIAMOND finds almost nothing on 7-30 aa queries.
        cmd = self._emit()["command"]
        for flag in ("--ultra-sensitive", "--short-query-ungapped-bitscore",
                     "--masking 0", "--gapped-filter-evalue 0", "--min-score 1"):
            self.assertIn(flag, cmd, f"missing short-peptide flag: {flag}")

    def test_max_target_seqs_defaults_to_something_lca_can_use(self):
        d = self._emit()
        self.assertGreaterEqual(d["max_target_seqs"], 2)
        self.assertTrue(d["lca_possible"])

    def test_max_target_seqs_1_is_flagged_not_silently_accepted(self):
        # The production ocelot search used 1 and could not run the LCA step its
        # own method recommends. That must not happen quietly again.
        d = self._emit("--max-target-seqs", "1")
        self.assertFalse(d["lca_possible"])
        self.assertTrue(any("LCA" in w for w in d["warnings"]))

    def test_standard_mode_warns_about_lost_short_peptides(self):
        d = self._emit("--mode", "standard")
        self.assertTrue(any("short" in w.lower() or "9,801" in w for w in d["warnings"]))

    def test_taxid_and_outfmt_carry_staxids(self):
        # LCA needs staxids in the output or it has nothing to work with.
        self.assertIn("staxids", self._emit()["command"])


class TestFdrMath(unittest.TestCase):
    def test_rate_normalisation_uses_the_peptide_universes(self):
        # Equal hit counts but unequal universes must NOT give 100% FDR.
        scores = [10.0, 9.0, 8.0, 7.0]
        tgt = [True, False, True, False]
        n, thr = rate_normalised_accept(scores, tgt, n_target=1000, n_decoy=4000)
        self.assertGreater(n, 0)

    def test_ties_move_together(self):
        # Bitscore is coarse and heavily tied; accepting a partial tie group
        # inflated the ocelot count 10,404 -> 11,558 (+11%).
        scores = [5.0, 5.0, 5.0, 5.0]
        tgt = [True, True, False, False]
        n, thr = rate_normalised_accept(scores, tgt, n_target=100, n_decoy=100, fdr=1.0)
        self.assertIn(n, (0, 2), "a tie group must be accepted whole or not at all")

    def test_default_model_is_the_five_feature_set(self):
        self.assertEqual(FEATURE_SETS["default"], DEFAULT_FEATURES)
        self.assertEqual(len(DEFAULT_FEATURES), 5)
        self.assertIn("mean_conf", DEFAULT_FEATURES,
                      "confidence is load-bearing: the homology-only model does "
                      "not separate at all on a standard search")
        self.assertIn("log_evalue", DEFAULT_FEATURES,
                      "calibrate on E-value, never raw bitscore")

    def test_bigger_feature_sets_exist_but_are_not_the_default(self):
        # They accept more peptides, less precisely, and recover FEWER correct.
        self.assertGreater(len(FEATURE_SETS["12feat"]), len(FEATURE_SETS["default"]))
        self.assertNotEqual(FEATURE_SETS["default"], FEATURE_SETS["12feat"])


class TestHitParsing(unittest.TestCase):
    def _write(self, tmp, rows):
        p = os.path.join(tmp, "h.tsv")
        with open(p, "w") as fh:
            for r in rows:
                fh.write("\t".join(str(x) for x in r) + "\n")
        return p

    def test_reads_the_7_column_layout(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self._write(tmp, [["PEPTIDEK", "sp1", 90.0, 8, "1e-5", 40.0, "9606"]])
            h = read_hits(p)
            self.assertEqual(len(h), 1)
            self.assertAlmostEqual(h[0]["bitscore"], 40.0)
            self.assertEqual(h[0]["staxids"], "9606")

    def test_reads_the_16_column_extended_layout(self):
        # Both layouts occur in real output; picking the wrong one silently reads
        # mismatch counts as E-values.
        with tempfile.TemporaryDirectory() as tmp:
            p = self._write(tmp, [["PEPTIDEK", "sp1", 90.0, 8, 1, 0, 1, 8, 10, 17,
                                   "1e-5", 40.0, "9606", "PEP", "PEP", "8"]])
            h = read_hits(p)
            self.assertEqual(len(h), 1)
            self.assertAlmostEqual(h[0]["bitscore"], 40.0)
            self.assertEqual(h[0]["staxids"], "9606")

    def test_keeps_only_the_best_hit_per_peptide(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self._write(tmp, [
                ["PEPTIDEK", "a", 90.0, 8, "1e-2", 30.0, "1"],
                ["PEPTIDEK", "b", 95.0, 8, "1e-9", 55.0, "2"]])
            h = read_hits(p)
            self.assertEqual(len(h), 1)
            self.assertAlmostEqual(h[0]["bitscore"], 55.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
