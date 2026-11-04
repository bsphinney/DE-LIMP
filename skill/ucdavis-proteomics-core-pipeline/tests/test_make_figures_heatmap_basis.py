"""make_figures.R: the heatmap says when it shows the most variable proteins, not differential ones
(2.14.2 triage I10).

A Core analysis (a Core member, 2026-10-21): no protein passed, so the heatmap fell back to the 50 most variable
proteins -- and was still titled and captioned "top 50 differential proteins". The basis is now
recorded where the proteins are picked, and the title and the figures.json caption come from it.

Needs Rscript with ggplot2; skipped otherwise (CI has no R).
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "make_figures.R")

from test_make_figures_release import have_r_ggplot, write_inputs   # noqa: E402


def de_rows(n_pass, n_prot=20):
    """(protein, logFC, P.Value, adj.P.Val): the first n_pass pass adj.P < 0.05, the rest not."""
    return [(f"P{p:02d}", 1.0, 1e-6 * (p + 1), 1e-3 * (p + 1) if p < n_pass else 0.5)
            for p in range(n_prot)]


@unittest.skipUnless(have_r_ggplot(), "Rscript with ggplot2 not available")
class HeatmapBasis(unittest.TestCase):
    def caption(self, n_pass, incomplete=()):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        tables = write_inputs(d, provenance={"adjp": 0.05, "logfc": 1}, de_rows=de_rows(n_pass))
        if incomplete:                      # these proteins miss a value in one sample
            path = os.path.join(tables, "Expression_Matrix.csv")
            with open(path, newline="") as fh:
                rows = list(csv.reader(fh))
            for r in rows[1:]:
                if r[0] in incomplete:
                    r[2] = "NA"
            with open(path, "w", newline="") as fh:
                csv.writer(fh).writerows(rows)
        out = os.path.join(d, "figs")
        r = subprocess.run(["Rscript", SCRIPT, "--de-dir", tables, "--conditions",
                            os.path.join(d, "conditions.csv"), "--outdir", out],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out, "figures.json")) as fh:
            fj = json.load(fh)
        caps = {f["file"]: f["caption"] for f in fj["figures"]}
        self.assertIn("heatmap_top.png", caps, fj.get("failed"))
        return caps["heatmap_top.png"]

    def test_fallback_caption_says_most_variable(self):
        cap = self.caption(3)
        self.assertIn("most variable", cap)
        self.assertIn("only 3 passed adj.P < 0.05", cap)
        self.assertIn("chosen by variance, not by the DE test", cap)
        self.assertNotIn("differential proteins (row", cap)
        self.assertNotIn("top 20 differential", cap)

    def test_passing_proteins_without_a_complete_row_are_counted_as_passing(self):
        # 2.14.2 review: the heatmap draws complete rows only; 8 passed but 4 have a gap, so the
        # fallback was taken -- and the caption said "fewer than 5 passed"
        cap = self.caption(8, incomplete={"P04", "P05", "P06", "P07"})
        self.assertIn("most variable", cap)
        self.assertIn("8 passed adj.P < 0.05", cap)
        self.assertIn("fewer than 5 of them have a value in every sample", cap)
        self.assertNotIn("fewer than 5 passed", cap)
        self.assertNotIn("not differential proteins", cap)

    def test_five_or_more_keep_differential_wording(self):              # (pin)
        cap = self.caption(6)
        self.assertIn("top 6 differential", cap)
        self.assertNotIn("most variable", cap)


if __name__ == "__main__":
    unittest.main()
