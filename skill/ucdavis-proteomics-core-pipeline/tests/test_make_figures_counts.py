"""make_figures.R: the "Proteins quantified per sample" plot only when it can say something.

A DPC/limpa expression matrix is complete by construction -- every protein gets a value in every
run -- so the old plot drew one identical bar per sample (Silva08172026, 2026-09-24: 30 bars of
6,112) and the report had to explain it away. It is skipped for a complete matrix, where
qc_detected_vs_inferred.png is the per-sample depth view, and kept for a matrix with missing
values (MaxLFQ), where the counts differ and mean something.

Needs Rscript with ggplot2; skipped otherwise (CI has no R).
"""
import csv
import json
import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "make_figures.R")


def have_r_ggplot():
    if not shutil.which("Rscript"):
        return False
    r = subprocess.run(["Rscript", "-e", "quit(status = !requireNamespace('ggplot2', quietly=TRUE))"],
                       capture_output=True)
    return r.returncode == 0


def write_inputs(d, missing):
    samples = [f"S{i}" for i in range(1, 7)]
    groups = ["A", "A", "A", "B", "B", "B"]
    with open(os.path.join(d, "conditions.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        w.writerows(zip(samples, groups))
    tables = os.path.join(d, "tables")
    os.makedirs(tables)
    with open(os.path.join(tables, "Expression_Matrix.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein.Group", "Genes"] + samples)
        for p in range(40):
            row = [20 + ((p * 7 + s * 3) % 11) / 3 for s in range(6)]
            if missing:
                row = [("NA" if p < 3 * s else v) for s, v in enumerate(row)]  # S6 misses 15
            w.writerow([f"P{p:05d}", f"G{p}"] + row)
    return tables


@unittest.skipUnless(have_r_ggplot(), "Rscript with ggplot2 not available")
class ProteinCountsPlot(unittest.TestCase):
    def run_figures(self, missing, stale=False):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        tables = write_inputs(d, missing)
        out = os.path.join(d, "figs")
        if stale:                                  # a copy left behind by an earlier run
            os.makedirs(out)
            with open(os.path.join(out, "qc_protein_counts.png"), "wb") as fh:
                fh.write(b"stale")
        r = subprocess.run(["Rscript", SCRIPT, "--de-dir", tables,
                            "--conditions", os.path.join(d, "conditions.csv"), "--outdir", out],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out, "figures.json")) as fh:
            listed = [f["file"] for f in json.load(fh)["figures"]]
        return out, listed, r.stderr + r.stdout

    def test_a_complete_matrix_gets_no_identical_bar_plot(self):
        out, listed, log = self.run_figures(missing=False)
        self.assertFalse(os.path.exists(os.path.join(out, "qc_protein_counts.png")))
        self.assertNotIn("qc_protein_counts.png", listed)
        self.assertIn("proteins-per-sample plot skipped", log)
        self.assertIn("pca.png", listed)                 # the other figures are unaffected

    def test_a_stale_copy_is_removed_when_the_plot_is_skipped(self):
        out, listed, log = self.run_figures(missing=False, stale=True)
        self.assertFalse(os.path.exists(os.path.join(out, "qc_protein_counts.png")))
        self.assertNotIn("qc_protein_counts.png", listed)
        self.assertIn("removed 1 figure(s) left by an earlier run", log)

    def test_a_stale_copy_is_replaced_when_the_plot_is_drawn(self):
        out, listed, log = self.run_figures(missing=True, stale=True)
        path = os.path.join(out, "qc_protein_counts.png")
        with open(path, "rb") as fh:
            self.assertTrue(fh.read(8).startswith(b"\x89PNG"))   # regenerated, not the stale bytes
        self.assertIn("qc_protein_counts.png", listed)
        self.assertNotIn("not drawn", log)

    def test_a_matrix_with_missing_values_keeps_the_plot(self):
        out, listed, _ = self.run_figures(missing=True)
        self.assertTrue(os.path.isfile(os.path.join(out, "qc_protein_counts.png")))
        self.assertIn("qc_protein_counts.png", listed)


if __name__ == "__main__":
    unittest.main()
