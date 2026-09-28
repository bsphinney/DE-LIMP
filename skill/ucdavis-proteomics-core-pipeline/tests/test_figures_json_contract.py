#!/usr/bin/env python3
"""Contract: make_figures.R's figures.json -> analysis_prompt.py's brief (and the HTML page).

make_figures.R writes figures.json; analysis_prompt.py turns it into the figure list the report
writer embeds. When the two drift -- 2.8.0 changed the file from a bare list to
{"figures": [...], "failed": [...]} -- the brief silently lists ZERO figures and every other
test stays green. These tests pin the hand-off:
  * the 2.8.0 object: the brief lists every figure and names every failed one with its reason;
  * the bare list sessions before 2.8.0 wrote: still read, every figure listed;
  * a figures.json that is unreadable or of an unknown shape is said (brief + stderr), never an
    empty figure list;
  * the brief and the HTML page read figures.json through ONE function;
  * end to end (needs Rscript + ggplot2): whatever make_figures.R in this tree writes, the brief
    lists all of it. CONTRACT_MAKE_FIGURES_R points it at another copy of the script.
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
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
MAKE_FIGURES = os.environ.get("CONTRACT_MAKE_FIGURES_R") or os.path.join(SCRIPTS, "make_figures.R")
sys.path.insert(0, SCRIPTS)

# The 2.8.0 shape, as make_figures.R writes it (fix/280-figures 3f2a0f6), reasons verbatim.
FIG_280 = {
    "figures": [
        {"file": "volcano_B.A.png", "type": "volcano", "caption": "Volcano plot for B vs A."},
        {"file": "violin_top_B.A.png", "type": "violin", "caption": "Top 12 proteins for B vs A."},
        {"file": "heatmap_top.png", "type": "heatmap", "caption": "Top proteins, z-scored."},
        {"file": "qc_detected_vs_inferred.png", "type": "qc", "caption": "Detected vs inferred."},
        {"file": "qc_pvalue_panel.png", "type": "pvalue",
         "caption": "Raw p-value distributions for all 1 contrast (appendix; a calibration "
                    "check, not a result)."}],
    "failed": [
        {"file": "pca.png", "type": "pca",
         "reason": "a PCA needs at least 3 samples and 5 varying proteins; this matrix has 2 "
                   "samples and 20 proteins"},
        {"file": "violin_top_X.Y.png", "type": "violin",
         "reason": "could not identify this contrast's groups from de_provenance.json or its "
                   "name"}]}
# Before 2.8.0: a bare list of figures, one p-value histogram per contrast.
FIG_LEGACY = [
    {"file": "volcano_B.A.png", "type": "volcano", "caption": "Volcano plot for B vs A."},
    {"file": "pvalue_B.A.png", "type": "pvalue", "caption": "Raw p-values for B vs A."},
    {"file": "pca.png", "type": "pca", "caption": "PCA of all samples."}]


def de_dir(root):
    """The smallest DE folder the brief reads: a record and one DE table."""
    d = os.path.join(root, "tables")
    os.makedirs(d)
    with open(os.path.join(d, "de_provenance.json"), "w") as fh:
        json.dump({"pipeline_id": "maxlfq", "method": "maxlfq", "adjp": 0.05,
                   "contrasts": ["B-A"]}, fh)
    with open(os.path.join(d, "DE_maxlfq_B.A.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein.Group", "Genes", "logFC", "P.Value", "adj.P.Val"])
        w.writerow(["P1", "G1", 2.0, 1e-6, 1e-4])
    return d


def run_brief(root, figures_json=None, raw=None):
    """-> (brief text, the JSON analysis_prompt prints, stderr)."""
    figs = os.path.join(root, "figures")
    os.makedirs(figs, exist_ok=True)
    if figures_json is not None or raw is not None:
        with open(os.path.join(figs, "figures.json"), "w") as fh:
            fh.write(raw if raw is not None else json.dumps(figures_json))
    tables = os.path.join(root, "tables")
    if not os.path.isdir(tables):
        de_dir(root)
    out = os.path.join(root, "brief.md")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                        "--out", out, "--de-dir", tables, "--figures-dir", figs],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(out, encoding="utf-8") as fh:
        return fh.read(), json.loads(r.stdout), r.stderr


def shape(fj):
    """(figures, failed) of a figures.json in either shape make_figures.R has written; any other
    shape is a contract break, not an empty list."""
    if isinstance(fj, list):
        return fj, []
    if isinstance(fj, dict) and isinstance(fj.get("figures"), list):
        return fj["figures"], fj.get("failed") or []
    raise AssertionError(f"figures.json shape drifted: {type(fj).__name__} "
                         f"{sorted(fj) if isinstance(fj, dict) else ''}")


class Carries:
    def assert_brief_carries(self, brief, out, fj):
        figures, failed = shape(fj)
        self.assertTrue(figures, "a contract test with no figures proves nothing")
        expected = [f["file"] for f in figures]
        # the one figure left out by design: per-sample counts of a complete-by-construction
        # matrix -- and then the brief says so
        if "Do NOT embed, reference or describe a proteins-quantified-per-sample plot" in brief:
            expected = [x for x in expected if not x.startswith("qc_protein_counts")]
        for x in expected:
            self.assertIn(f"`figures/{x}`", brief, f"figure dropped from the brief: {x}")
        for f in failed:
            self.assertIn(f"`figures/{f['file']}`", brief, f"failed figure not named: {f}")
            if f.get("reason"):
                self.assertIn(f["reason"], brief)
        self.assertEqual(out["figures"], expected)
        self.assertEqual(out["figures_failed"], [f["file"] for f in failed])
        self.assertIsNone(out["figures_error"])


class Contract(Carries, unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name

    def tearDown(self):
        self._td.cleanup()

    def test_280_object_lists_every_figure_and_names_each_failed_one(self):
        brief, out, _ = run_brief(self.root, FIG_280)
        self.assert_brief_carries(brief, out, FIG_280)
        self.assertIn("**Not drawn in this run** (make_figures.R's `failed` list)", brief)
        # four to embed; the p-value panel is the HTML report's appendix
        self.assertEqual(out["n_figures"], 4)
        self.assertEqual(out["figures_appendix"], ["qc_pvalue_panel.png"])
        self.assertIn("embed all 4 figures", out["next"])

    def test_legacy_bare_list_still_lists_every_figure(self):
        brief, out, err = run_brief(self.root, FIG_LEGACY)
        self.assert_brief_carries(brief, out, FIG_LEGACY)
        self.assertEqual(out["n_figures"], 3)
        self.assertNotIn("Not drawn in this run", brief)
        self.assertNotIn("WARNING", err)

    def test_a_drifted_or_unreadable_figures_json_is_said_not_silently_empty(self):
        for i, (fj, raw) in enumerate((({"plots": FIG_LEGACY}, None), (None, "{not json"),
                                       ({"figures": "volcano_B.A.png"}, None))):
            with self.subTest(case=i):
                root = tempfile.mkdtemp(dir=self.root)
                brief, out, err = run_brief(root, fj, raw)
                self.assertIn("WARNING", err)
                self.assertIn("**No figure list:** `figures/figures.json`", brief)
                self.assertTrue(out["figures_error"])
                self.assertEqual(out["figures"], [])

    def test_one_reader_for_the_brief_and_the_page(self):
        import analysis_prompt
        import make_analysis_html
        self.assertIs(analysis_prompt.read_figures_json, make_analysis_html.read_figures_json)
        for fj in (FIG_280, FIG_LEGACY):
            d = tempfile.mkdtemp(dir=self.root)
            with open(os.path.join(d, "figures.json"), "w") as fh:
                json.dump(fj, fh)
            got = make_analysis_html.read_figures_json(d)
            figures, failed = shape(fj)
            self.assertEqual([f["file"] for f in got["figures"]], [f["file"] for f in figures])
            self.assertEqual([(f["file"], f["reason"]) for f in got["failed"]],
                             [(f["file"], f["reason"]) for f in failed])

    def test_the_page_names_each_failed_figure_too(self):
        # no report: the page is the figure gallery, built from the same figures.json
        tables = de_dir(self.root)
        figs = os.path.join(self.root, "figures")
        os.makedirs(figs)
        with open(os.path.join(figs, "figures.json"), "w") as fh:
            json.dump(FIG_280, fh)
        for f in FIG_280["figures"]:
            with open(os.path.join(figs, f["file"]), "wb") as fh:
                fh.write(b"\x89PNG\r\n\x1a\n" + f["file"].encode())
        html_out = os.path.join(self.root, "Analysis_Report.html")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--figures", figs, "--tables", tables, "--out", html_out, "--no-pdf"],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        res = json.loads(r.stdout)
        self.assertEqual(res["figures_embedded"], len(FIG_280["figures"]))
        with open(html_out[:-5] + ".md", encoding="utf-8") as fh:
            md = fh.read()
        self.assertIn("2 figures could not be drawn in this run", md)
        for f in FIG_280["failed"]:
            self.assertIn(f"`{f['file']}`: {f['reason']}", md)
        self.assertIn("## Appendix: p-value calibration", md)


def have_r_ggplot():
    if not shutil.which("Rscript"):
        return False
    r = subprocess.run(["Rscript", "-e", "quit(status = !requireNamespace('ggplot2', quietly=TRUE))"],
                       capture_output=True)
    return r.returncode == 0


@unittest.skipUnless(have_r_ggplot(), "Rscript with ggplot2 not available")
class EndToEnd(Carries, unittest.TestCase):
    """The real writer: make_figures.R on a two-sample study (so a 2.8.0 make_figures cannot
    draw the PCA and says so), then analysis_prompt.py on what it wrote."""

    def test_the_brief_carries_whatever_make_figures_writes(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [("run_A1", "A"), ("run_B1", "B")]
            with open(os.path.join(d, "conditions.csv"), "w", newline="") as fh:
                csv.writer(fh).writerows([("File.Name", "Group")] + runs)
            tables = os.path.join(d, "tables")
            os.makedirs(tables)
            with open(os.path.join(tables, "Expression_Matrix.csv"), "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["Protein.Group", "Genes"] + [r for r, _ in runs])
                for p in range(20):
                    w.writerow([f"P{p:02d}", f"G{p}", 20 + p / 7, 20.5 + p / 9])
            with open(os.path.join(tables, "Detection_Matrix.csv"), "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["Protein.Group"] + [r for r, _ in runs])
                for p in range(20):
                    w.writerow([f"P{p:02d}", 2, 0 if p == 0 else 2])
            with open(os.path.join(tables, "DE_dpc_B.A.csv"), "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["Protein.Group", "Genes", "logFC", "P.Value", "adj.P.Val"])
                for p in range(20):
                    w.writerow([f"P{p:02d}", f"G{p}", 1.0, 1e-6 * (p + 1), 1e-4 * (p + 1)])
            with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
                json.dump({"pipeline_id": "dpc", "method": "dpc", "adjp": 0.05, "logfc": 1,
                           "contrasts": ["B-A"], "rollup_method": "DPC-Quant"}, fh)
            figs = os.path.join(d, "figures")
            os.makedirs(figs)
            r = subprocess.run(["Rscript", MAKE_FIGURES, "--de-dir", tables, "--conditions",
                                os.path.join(d, "conditions.csv"), "--outdir", figs],
                               capture_output=True, text=True, timeout=300)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(figs, "figures.json")) as fh:
                fj = json.load(fh)
            brief, out, _ = run_brief(d)
            self.assert_brief_carries(brief, out, fj)


if __name__ == "__main__":
    unittest.main(verbosity=2)
