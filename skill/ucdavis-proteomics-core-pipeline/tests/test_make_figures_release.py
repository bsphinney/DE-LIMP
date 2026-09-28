"""make_figures.R release-review fixes (2.8.0).

- Nothing stale, nothing silently absent: every PNG the script owns is deleted before drawing
  (only those names), and a figure it meant to draw and did not is listed in figures.json
  "failed" with the reason -- including a PCA skipped for fewer than 3 samples.
- One significance cutoff: de_provenance.json's adjp/logfc win; --adjp/--logfc are only a
  fallback, and the captions then say where the cutoff came from.
- Violins: ties ranked by raw p-value (topTable order), the arrow is the model's log2FC (never
  the plain-mean difference, which can disagree in sign), groups with <= 4 runs are points and
  a mean bar (no density), and a within-block contrast joins each block's runs.
- The p-value histograms are one panel, qc_pvalue_panel.png; the heatmap caption gives the
  inferred share; the volcano says "estimated", not "measured".

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


def write_inputs(d, n=3, provenance=None, block=False, de_rows=None, detection=True, n_prot=20, de_name="B.A"):
    """Two groups A/B (n runs each). de_rows: list of (protein, logFC, P.Value, adj.P.Val)."""
    runs = [(f"run_{g}{i}", g, f"S{i}") for g in ("A", "B") for i in range(1, n + 1)]
    with open(os.path.join(d, "conditions.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"] + (["Subject"] if block else []))
        w.writerows([(r, g, s) if block else (r, g) for r, g, s in runs])
    tables = os.path.join(d, "tables")
    os.makedirs(tables)
    with open(os.path.join(tables, "Expression_Matrix.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein.Group", "Genes"] + [r for r, _, _ in runs])
        for p in range(n_prot):
            # group B's plain mean is LOWER than A's for every protein; the DE rows below can say
            # otherwise, as a model with weights or a block effect may
            w.writerow([f"P{p:02d}", f"G{p}"] + [20 + (p % 5) / 2 + (0 if g == "A" else -0.3) + i / 10
                                                  for i, (_, g, _) in enumerate(runs)])
    if detection:
        with open(os.path.join(tables, "Detection_Matrix.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Protein.Group"] + [r for r, _, _ in runs])
            for p in range(n_prot):
                w.writerow([f"P{p:02d}"] + [0 if (p == 0 and g == "A") else 2 for _, g, _ in runs])
    rows = de_rows or [(f"P{p:02d}", 1.0, 1e-6 * (p + 1), 1e-4 * (p + 1)) for p in range(n_prot)]
    with open(os.path.join(tables, f"DE_dpc_{de_name}.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein.Group", "Genes", "logFC", "P.Value", "adj.P.Val"])
        for pid, lfc, p, adj in rows:
            w.writerow([pid, "G" + pid[1:].lstrip("0") if pid[1:].lstrip("0") else "G0", lfc, p, adj])
    if provenance is not None:
        prov = {"pipeline_id": "dpc", "rollup_method": "DPC-Quant (Detection Probability Curve)",
                "contrasts": ["B-A"]}
        prov.update(provenance)
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
    return tables


@unittest.skipUnless(have_r_ggplot(), "Rscript with ggplot2 not available")
class ReleaseFixes(unittest.TestCase):
    def run_figures(self, extra=(), pre_existing=(), **kw):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        tables = write_inputs(d, **kw)
        out = os.path.join(d, "figs")
        os.makedirs(out)
        for name in pre_existing:
            with open(os.path.join(out, name), "wb") as fh:
                fh.write(b"stale")
        r = subprocess.run(["Rscript", SCRIPT, "--de-dir", tables, "--conditions",
                            os.path.join(d, "conditions.csv"), "--outdir", out, *extra],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out, "figures.json")) as fh:
            fj = json.load(fh)
        caps = {f["file"]: f["caption"] for f in fj["figures"]}
        return out, caps, fj["failed"], r.stderr + r.stdout

    def log_line(self, log, prefix):
        hits = [ln[len(prefix):] for ln in log.splitlines() if ln.startswith(prefix)]
        self.assertEqual(len(hits), 1, f"no single '{prefix}' line in:\n{log}")
        return hits[0]

    # ---- 1. stale figures + failed list ------------------------------------
    def test_owned_figures_are_swept_and_foreign_files_kept(self):
        out, caps, failed, log = self.run_figures(
            provenance={"adjp": 0.05, "logfc": 1},
            pre_existing=("pvalue_B.A.png", "violin_top_old.png", "volcano_gone.png", "notes.png"))
        left = set(os.listdir(out))
        for gone in ("pvalue_B.A.png", "violin_top_old.png", "volcano_gone.png"):
            self.assertNotIn(gone, left)
        self.assertIn("notes.png", left)                      # not ours: never touched
        self.assertIn("removed 3 figure(s) left by an earlier run", log)
        self.assertEqual(failed, [])

    def test_pca_skipped_for_two_samples_is_recorded_with_a_reason(self):
        out, caps, failed, _ = self.run_figures(n=1, provenance={"adjp": 0.05, "logfc": 1})
        self.assertNotIn("pca.png", caps)
        pca = [f for f in failed if f["file"] == "pca.png"]
        self.assertEqual(len(pca), 1)
        self.assertIn("at least 3 samples", pca[0]["reason"])
        self.assertIn("2 samples", pca[0]["reason"])

    def test_a_violin_that_cannot_be_drawn_is_recorded(self):
        # a contrast between groups that are not in conditions.csv (nor in the record)
        out, caps, failed, _ = self.run_figures(provenance={"adjp": 0.05, "logfc": 1, "contrasts": ["X-Y"]},
                                                de_name="X.Y")
        self.assertIn("volcano_X.Y.png", caps)                # the volcano needs no groups
        v = [f for f in failed if f["file"] == "violin_top_X.Y.png"]
        self.assertEqual(len(v), 1)
        self.assertIn("could not identify", v[0]["reason"])

    # ---- 2. one cutoff ------------------------------------------------------
    def test_the_de_record_cutoff_wins_over_the_flag(self):
        _, caps, _, log = self.run_figures(extra=("--adjp", "0.2"), provenance={"adjp": 0.01, "logfc": 1})
        self.assertIn("--adjp 0.2 ignored: the DE run used adj.P < 0.01", log)
        self.assertIn("adj.P < 0.01", caps["volcano_B.A.png"])
        self.assertNotIn("DEFAULT", caps["volcano_B.A.png"])

    def test_a_fallback_cutoff_is_tagged(self):
        _, caps, _, _ = self.run_figures()                     # no record, no flag
        self.assertIn("adj.P < 0.05 is make_figures.R's DEFAULT, not user-confirmed", caps["volcano_B.A.png"])
        self.assertIn("not user-confirmed", caps["violin_top_B.A.png"])
        _, caps, _, _ = self.run_figures(extra=("--adjp", "0.1"))     # flag given, no record
        cap = caps["volcano_B.A.png"]
        self.assertIn("adj.P < 0.1 from the command line", cap)
        self.assertIn("2-fold reference line is make_figures.R's DEFAULT", cap)   # --logfc was not given

    # ---- 3. violins ---------------------------------------------------------
    def test_ties_are_ranked_by_raw_p_not_fold_change(self):
        # all tie on adj.P; |logFC| grows while the raw p-value gets WORSE
        rows = [(f"P{p:02d}", 0.5 + p, 1e-6 * (p + 1), 1e-3) for p in range(12)]
        _, _, _, log = self.run_figures(provenance={"adjp": 0.05, "logfc": 1}, de_rows=rows, n_prot=12)
        chosen = self.log_line(log, "[figures] violin B.A proteins: ").split(", ")
        self.assertEqual(chosen, [f"P{p:02d}" for p in range(8)])

    def test_the_arrow_is_the_models_log2fc_not_the_plain_means(self):
        # B's plain mean is below A's for every protein; the model says +1.5
        rows = [(f"P{p:02d}", 1.5, 1e-6 * (p + 1), 1e-4 * (p + 1)) for p in range(20)]
        _, caps, _, log = self.run_figures(provenance={"adjp": 0.05, "logfc": 1}, de_rows=rows)
        arrows = self.log_line(log, "[figures] violin B.A arrows (model log2FC): ")
        self.assertTrue(all(a.endswith("+1.50") for a in arrows.split(", ")), arrows)
        self.assertIn("the arrow is the model's log2 fold change", caps["violin_top_B.A.png"])

    def test_small_groups_are_points_and_larger_ones_violins(self):
        _, caps, _, log = self.run_figures(n=3, provenance={"adjp": 0.05, "logfc": 1})
        self.assertIn("points + mean bars", self.log_line(log, "[figures] violin B.A drawn as: "))
        self.assertIn("no density is drawn", caps["violin_top_B.A.png"])
        _, caps, _, log = self.run_figures(n=5, provenance={"adjp": 0.05, "logfc": 1})
        self.assertEqual(self.log_line(log, "[figures] violin B.A drawn as: "), "violins")

    def test_a_within_block_contrast_joins_each_blocks_runs(self):
        block = {"applied": True, "column": "Subject", "contrast_structure": {"B-A": "within"}}
        _, caps, _, log = self.run_figures(block=True, provenance={"adjp": 0.05, "logfc": 1, "block": block})
        self.assertIn("paired lines join each Subject (one run per group)", log)
        self.assertIn("within-Subject contrast", caps["violin_top_B.A.png"])
        # a between-block contrast gets no lines
        block["contrast_structure"] = {"B-A": "between"}
        _, _, _, log = self.run_figures(block=True, provenance={"adjp": 0.05, "logfc": 1, "block": block})
        self.assertNotIn("paired lines", log)

    # ---- 4-7. captions + the p-value panel ----------------------------------
    def test_one_pvalue_panel_and_honest_captions(self):
        out, caps, _, _ = self.run_figures(provenance={"adjp": 0.05, "logfc": 1})
        self.assertIn("qc_pvalue_panel.png", caps)
        self.assertFalse([f for f in os.listdir(out) if f.startswith("pvalue_")])
        self.assertIn("confidently estimated small change", caps["volcano_B.A.png"])
        self.assertNotIn("confidently measured", caps["volcano_B.A.png"])
        self.assertRegex(caps["heatmap_top.png"], r"\d+% of the cells shown are inferred")
        self.assertIn("not agreement between replicates", caps["heatmap_top.png"])


if __name__ == "__main__":
    unittest.main()
