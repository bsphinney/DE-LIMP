#!/usr/bin/env python3
"""make_analysis_html.py embeds the figures the report is about -- not every file in figures/.

msalemi, 2026-09-24 (Silva08172026): after redrawing pca.png she kept the old image as
figures/pca_original_labels.png; Analysis_Report.html embedded 29 figures although
AI_Analysis_Report.md references 28, so the stale PCA went out as an extra figure.

Guards (stdlib only): the report's image references decide; figures.json decides when there
is no report (or it references none); an image in neither is NOT embedded and one warning
line names it; a reference to a missing image is warned about; with no list at all every
image is embedded and that is said, not done silently.
"""
import base64
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import make_analysis_html as mah  # noqa: E402

FIGS = ["pca.png", "pca_original_labels.png", "volcano_B.A.png", "qc_protein_counts.png"]
REPORT = ("# Analysis\n\n![PCA of the samples](figures/pca.png)\n\n"
          "Text.\n\n![Volcano — B vs A](figures/volcano_B.A.png \"volcano\")\n"
          "<img src=\"figures/qc_protein_counts.png\" alt=\"counts\">\n")


def png_bytes(name):
    return b"\x89PNG\r\n\x1a\n" + name.encode() * 4       # distinct per file


def b64(name):
    return base64.b64encode(png_bytes(name)).decode("ascii")


class FigureSelection(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.figs = os.path.join(self.root, "figures")
        os.makedirs(self.figs)
        for f in FIGS:
            with open(os.path.join(self.figs, f), "wb") as fh:
                fh.write(png_bytes(f))

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, text):
        path = os.path.join(self.root, name)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return path

    def build(self, *args):
        out = os.path.join(self.root, "report.html")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--figures", self.figs, "--out", out, *args],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out, encoding="utf-8") as fh:
            return json.loads(r.stdout), fh.read(), r.stderr

    def warnings(self, stderr):
        return [ln for ln in stderr.splitlines() if "WARNING" in ln]

    def test_only_figures_the_report_references_are_embedded(self):
        res, page, err = self.build("--report", self.write("AI_Analysis_Report.md", REPORT))
        self.assertEqual(res["figures_embedded"], 3)
        self.assertEqual(res["figures_not_embedded"], ["pca_original_labels.png"])
        self.assertEqual(res["figure_list"], "AI_Analysis_Report.md")
        self.assertNotIn(b64("pca_original_labels.png"), page)
        for f in ("pca.png", "volcano_B.A.png", "qc_protein_counts.png"):
            self.assertIn(b64(f), page, f)
        w = self.warnings(err)
        self.assertEqual(len(w), 1, err)                     # ONE line, naming the image
        self.assertIn("pca_original_labels.png", w[0])
        self.assertIn("NOT embedded", w[0])

    def test_figures_json_decides_without_a_report(self):
        with open(os.path.join(self.figs, "figures.json"), "w") as fh:
            json.dump([{"file": "pca.png", "caption": "PCA"},
                       {"file": "volcano_B.A.png", "caption": "Volcano"}], fh)
        res, page, err = self.build()
        self.assertEqual(res["figure_list"], "figures.json")
        self.assertEqual(sorted(res["figures_not_embedded"]),
                         ["pca_original_labels.png", "qc_protein_counts.png"])
        self.assertNotIn(b64("pca_original_labels.png"), page)
        self.assertEqual(len(self.warnings(err)), 1)

    def test_a_report_without_images_falls_back_to_figures_json(self):
        with open(os.path.join(self.figs, "figures.json"), "w") as fh:
            json.dump([{"file": "pca.png"}], fh)
        res, _, _ = self.build("--report", self.write("r.md", "# Analysis\n\nNo figures.\n"))
        self.assertEqual(res["figure_list"], "figures.json")
        self.assertEqual(res["figures_embedded"], 1)

    def test_a_missing_referenced_image_is_warned_about(self):
        res, _, err = self.build("--report", self.write(
            "r.md", REPORT + "\n![Heatmap](figures/heatmap_top.png)\n"))
        self.assertEqual(res["figures_missing"], ["heatmap_top.png"])
        self.assertTrue(any("heatmap_top.png" in w and "not in" in w for w in self.warnings(err)))

    def test_with_no_list_everything_is_embedded_and_said(self):
        res, _, err = self.build()
        self.assertEqual(res["figures_embedded"], len(FIGS))
        self.assertEqual(res["figure_list"], "all (no list)")
        self.assertTrue(any("embedding all" in w for w in self.warnings(err)))

    def test_reference_parser(self):
        md = ('![a](figures/pca.png) ![b](<figures/my%20fig.png> "t") <IMG SRC="figures/v.png"> '
              '![again](./figures/pca.png) ![remote](https://x.org/y.png) ![c](figs/z.svg?v=2)')
        self.assertEqual(mah.report_figures(md), ["pca.png", "my fig.png", "v.png", "z.svg"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
