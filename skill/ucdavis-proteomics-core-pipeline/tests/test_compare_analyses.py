#!/usr/bin/env python3
"""compare_analyses.R counts TESTED proteins: the rows with a P.Value.

A SILAC re-analysis (2026-10) compared two DE_silac tables whose untested rows (Evidence "ratio in
k of n replicates: not tested", no P.Value) the Comparator counted as part of the protein
universe and as "not significant" in the 3x3 concordance. run_de.R's "not measured in either
group" rows (2.12) are the same. SKILL.md: count tested proteins as the rows with a P.Value.
Needs Rscript only (base R).
"""
import csv
import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
HEAD = ["Protein.Group", "logFC", "P.Value", "adj.P.Val", "Evidence"]


def table(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEAD)
        w.writerows(rows)


@unittest.skipUnless(shutil.which("Rscript"), "needs Rscript")
class TestedProteinsOnly(unittest.TestCase):
    def test_untested_rows_are_not_in_the_universe_or_the_concordance(self):
        tested = [["P01", 2.0, 1e-6, 1e-4, "tested"], ["P02", -1.5, 1e-5, 1e-3, "tested"],
                  ["P03", 0.1, 0.5, 0.8, "tested"]]
        untested = [["P04", "", "", "", "ratio in 2 of 4 replicates: not tested (minimum 3)"],
                    ["P05", "", "", "", "ratio in 1 of 4 replicates: not tested (minimum 3)"]]
        with tempfile.TemporaryDirectory() as d:
            a, b, out = (os.path.join(d, x) for x in ("a", "b", "cmp"))
            table(os.path.join(a, "DE_silac_KO.WT.csv"), tested + untested)
            # P04 is tested in B only: shared as a listed row, but not tested in both
            table(os.path.join(b, "DE_silac_KO.WT.csv"),
                  tested + [["P04", 0.2, 0.4, 0.7, "tested"], untested[1]])
            r = subprocess.run(["Rscript", os.path.join(SCRIPTS, "compare_analyses.R"),
                                "--out", out, "--analysis", f"Original:{a}",
                                "--analysis", f"Redo:{b}"], capture_output=True, text=True,
                               timeout=300)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(out, "protein_universe.csv")) as fh:
                uni = {x["Analysis"]: x for x in csv.DictReader(fh)}
            self.assertEqual((uni["Original"]["n_proteins"], uni["Original"]["n_tested"],
                              uni["Original"]["n_significant"]), ("5", "3", "2"))
            self.assertEqual(uni["Redo"]["n_tested"], "4")
            with open(os.path.join(out, "concordance_summary.csv")) as fh:
                summ = next(csv.DictReader(fh))
            self.assertEqual((summ["n_shared"], summ["sig_both"]), ("3", "2"))
            with open(os.path.join(out, "COMPARISON.md")) as fh:
                md = fh.read()
            self.assertIn("| Original | 5 | 3 | 2 |", md)
            self.assertIn("proteins tested in both: 3", md)


if __name__ == "__main__":
    unittest.main()
