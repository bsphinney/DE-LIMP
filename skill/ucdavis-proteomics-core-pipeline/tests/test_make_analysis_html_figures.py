#!/usr/bin/env python3
"""make_analysis_html.py puts each figure where the report discusses it, and nothing else.

msalemi, 2026-09-24 (Silva08172026): Analysis_Report.html printed all 28 figure references
of AI_Analysis_Report.md as literal Markdown (`![Proteins quantified per sample](figures/
qc_protein_counts.png)`), dumped every image in figures/ -- a stale pca_original_labels.png
included -- into galleries at the top, and showed "Overview" and "Audit & caveats" twice.

Guards (stdlib only):
  * an image reference becomes an embedded figure AT ITS POSITION; `<img` count equals the
    referenced images; no literal `![` survives;
  * an unreferenced image is left out with one warning line; a missing one is a visible
    "figure missing" note; a path outside the session or a URL is never embedded;
  * galleries only without a report, and then figures.json's figures only;
  * no h2 title appears twice.
"""
import base64
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import make_analysis_html as mah  # noqa: E402

FIGS = ["pca.png", "pca_original_labels.png", "volcano_B.A.png", "qc_protein_counts.png"]
REPORT = ("# Study title\n\n## Overview\n\nBEFORE-PCA text.\n\n![PCA of the samples](figures/pca.png)\n\n"
          "AFTER-PCA text.\n\n## Key findings\n\nSee ![Volcano — B vs A](figures/volcano_B.A.png \"v\") "
          "inline.\n\n<img src=\"figures/qc_protein_counts.png\" alt=\"counts\">\n\n"
          "## Audit & caveats\n\n- the report's own audit\n")


def png_bytes(name):
    return b"\x89PNG\r\n\x1a\n" + name.encode() * 4       # distinct per file


def b64(name):
    return base64.b64encode(png_bytes(name)).decode("ascii")


class FigureSelection(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.session = self._td.name
        self.out_dir = os.path.join(self.session, "output")
        self.figs = os.path.join(self.out_dir, "figures")
        os.makedirs(self.figs)
        for f in FIGS:
            with open(os.path.join(self.figs, f), "wb") as fh:
                fh.write(png_bytes(f))

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, text):
        path = os.path.join(self.out_dir, name)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return path

    def build(self, *args, session=True):
        out = os.path.join(self.session, "report.html")
        base = ["--session", self.session] if session else ["--figures", self.figs]
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            *base, "--out", out, *args], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out, encoding="utf-8") as fh:
            return json.loads(r.stdout), fh.read(), r.stderr

    @staticmethod
    def warnings(stderr):
        return [ln for ln in stderr.splitlines() if "WARNING" in ln]

    def test_figures_are_embedded_where_the_text_references_them(self):
        self.write("AI_Analysis_Report.md", REPORT)
        res, page, err = self.build()
        self.assertEqual(res["figures_embedded"], 3)
        self.assertEqual(page.count("<img"), 3)                 # one per referenced image
        self.assertNotIn("![", page)                            # never literal Markdown
        self.assertNotIn("&lt;img", page)
        # at its position: between the paragraphs around the reference
        self.assertLess(page.index("BEFORE-PCA"), page.index(b64("pca.png")))
        self.assertLess(page.index(b64("pca.png")), page.index("AFTER-PCA"))
        self.assertIn('Figure 1.</span><span class="figt">PCA of the samples', page)
        self.assertIn('Figure 2.</span><span class="figt">Volcano — B vs A', page)

    def test_unreferenced_image_is_left_out_with_one_warning(self):
        self.write("AI_Analysis_Report.md", REPORT)
        res, page, err = self.build()
        self.assertEqual(res["figures_not_embedded"], ["pca_original_labels.png"])
        self.assertNotIn(b64("pca_original_labels.png"), page)
        w = self.warnings(err)
        self.assertEqual(len(w), 1, err)
        self.assertIn("pca_original_labels.png", w[0])
        self.assertIn("NOT embedded", w[0])

    def test_a_missing_image_is_a_visible_note(self):
        os.remove(os.path.join(self.figs, "qc_protein_counts.png"))   # make_figures no longer draws it
        self.write("AI_Analysis_Report.md", REPORT)
        res, page, err = self.build()
        self.assertEqual(res["figures_missing"], ["qc_protein_counts.png"])
        self.assertIn("figure missing: qc_protein_counts.png", page)
        self.assertNotIn("![", page)
        self.assertTrue(any("missing" in w and "qc_protein_counts.png" in w
                            for w in self.warnings(err)))

    def test_paths_outside_the_session_and_urls_are_never_embedded(self):
        outside = os.path.join(os.path.dirname(self.session), "secret.png")
        self.write("AI_Analysis_Report.md",
                   "# T\n\n![a](../../secret.png)\n\n![b](https://example.org/x.png)\n\n"
                   f"![c]({outside})\n\n![d](figures/pca.png)\n")
        res, page, _ = self.build()
        self.assertEqual(res["figures_embedded"], 1)
        self.assertEqual(len(res["figures_rejected"]) + len(res["figures_missing"]), 3)
        self.assertNotIn("https://example.org", page.split("<body")[1].split("figure not embedded")[0])

    def test_no_h2_title_appears_twice(self):
        self.write("AI_Analysis_Report.md", REPORT)
        self.write("AUDIT.md", "# Results audit\n\n- ⚠️ **x** — the tool's audit\n")
        self.write("SAMPLE_QUALITY.md", "# Sample quality\n\nNo panel confounded.\n")
        _, page, _ = self.build()
        titles = [re.sub(r"<[^>]+>", "", t) for t in re.findall(r"<h2[^>]*>(.*?)</h2>", page)]
        self.assertEqual(len(titles), len(set(titles)), titles)
        self.assertNotIn("the tool's audit", page)        # the report's own audit wins
        ids = re.findall(r'\bid="([^"]+)"', page)
        self.assertEqual(len(ids), len(set(ids)))

    def test_galleries_only_without_a_report_and_only_figures_json(self):
        with open(os.path.join(self.figs, "figures.json"), "w") as fh:
            json.dump([{"file": "pca.png", "caption": "PCA"},
                       {"file": "volcano_B.A.png", "caption": "Volcano"}], fh)
        res, page, err = self.build()                    # no AI_Analysis_Report.md
        self.assertEqual(res["figures_embedded"], 2)
        self.assertEqual(sorted(res["figures_not_embedded"]),
                         ["pca_original_labels.png", "qc_protein_counts.png"])
        self.assertEqual(len(self.warnings(err)), 1)

    def test_a_report_without_images_embeds_none(self):
        with open(os.path.join(self.figs, "figures.json"), "w") as fh:
            json.dump([{"file": "pca.png"}], fh)
        self.write("AI_Analysis_Report.md", "# Analysis\n\nNo figures.\n")
        res, page, _ = self.build()
        self.assertEqual(res["figures_embedded"], 0)
        self.assertEqual(page.count("<img"), 0)

    def test_results_at_a_glance_counts_adj_p_only(self):
        tables = os.path.join(self.out_dir, "tables")
        os.makedirs(tables)
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump({"adjp": 0.05, "logfc": 1, "logfc_role": "reference_line_only"}, fh)
        with open(os.path.join(tables, "DE_dpc_B.A.csv"), "w") as fh:
            fh.write("Protein.Group,logFC,adj.P.Val\nP1,0.3,0.01\nP2,-2.0,0.01\n"
                     "P3,3.0,0.2\nP4,0.1,0.049\n")
        self.write("AI_Analysis_Report.md", "# T\n\nText.\n")
        _, page, _ = self.build()
        tile = re.search(r'<div class="v">(\d+)</div><div class="l">B vs A</div>'
                         r'<div class="s">(.*?)</div>', page)
        self.assertIsNotNone(tile, "no stat tile for the contrast")
        self.assertEqual(tile.group(1), "3")                        # 0.3-log2FC counts, 3.0 does not
        self.assertEqual(re.findall(r"(\d+) (?:up|down)", tile.group(2)), ["2", "1"])
        self.assertIn("no fold-change filter", page)

    def test_reference_parser(self):
        md = ('![a](figures/pca.png) ![b](<figures/my%20fig.png> "t") <IMG SRC="figures/v.png"> '
              '![again](./figures/pca.png) ![remote](https://x.org/y.png) ![c](figs/z.svg?v=2)')
        self.assertEqual(mah.report_figures(md), ["pca.png", "my fig.png", "v.png", "z.svg"])


class Restyle(unittest.TestCase):
    """report_style.py: one self-contained look -- no fetched fonts/CSS/JS, working contents
    links, a real dark mode, print styles, and callouts for what a PI must not miss."""

    @classmethod
    def setUpClass(cls):
        cls._td = tempfile.TemporaryDirectory()
        s = cls._td.name
        figs = os.path.join(s, "output", "figures")
        tables = os.path.join(s, "output", "tables")
        os.makedirs(figs)
        os.makedirs(tables)
        with open(os.path.join(figs, "pca.png"), "wb") as fh:
            fh.write(png_bytes("pca.png"))
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump({"adjp": 0.05, "n_samples": 6, "groups": {"A": 3, "B": 3},
                       "contrasts": ["B-A"], "display_label": "DPC-Quant + limma (limpa)"}, fh)
        with open(os.path.join(tables, "DE_dpc_B.A.csv"), "w") as fh:
            fh.write("Protein.Group,logFC,adj.P.Val\nP1,1,0.01\n")
        with open(os.path.join(tables, "QC_detected_vs_inferred.csv"), "w") as fh:
            fh.write("Sample,Group,Detected,Inferred,Total,PctDetected,PctInferred\n"
                     "S1,A,40,60,100,40,60\nS2,B,80,20,100,80,20\n")
        with open(os.path.join(s, "output", "AI_Analysis_Report.md"), "w") as fh:
            fh.write("# Study\n\n*Standfirst.*\n\n## Overview\n\n![PCA](figures/pca.png)\n\n"
                     "| Protein | logFC |\n|---|---|\n| P1 | 1.5 |\n| P2 | -0.3 |\n\n"
                     "## Data Quality Notes\n\n1. one\n\n- detail\n\n2. two\n\n"
                     "## Expert Review Notes\n\n- **Critical:** something\n")
        # the header's Submission fact is the attached record's label (submission_report.py)
        import submission_report
        os.makedirs(os.path.join(s, "input"))
        submission_report.attach(s, {"internal_id": "PROT_0001"})
        out = os.path.join(s, "r.html")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", s, "--out", out],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        with open(out, encoding="utf-8") as fh:
            cls.page = fh.read()

    @classmethod
    def tearDownClass(cls):
        cls._td.cleanup()

    def test_self_contained(self):
        p = self.page
        self.assertNotRegex(p, r"<link\b")
        self.assertNotRegex(p, r"<script[^>]+\bsrc=")
        self.assertNotRegex(p, r"@import|@font-face")
        self.assertNotRegex(p, r"url\(\s*['\"]?(https?:)?//")
        self.assertNotRegex(p, r"(?:src|href)=[\"']?(https?:)?//")

    def test_toc_links_resolve(self):
        ids = set(re.findall(r'\bid="([^"]+)"', self.page))
        hrefs = re.findall(r'href="#([^"]+)"', self.page)
        self.assertGreaterEqual(len(hrefs), 4)
        self.assertEqual([h for h in hrefs if h not in ids], [])

    def test_dark_mode_tokens_and_print_styles(self):
        p = self.page
        self.assertRegex(p, r"@media \(prefers-color-scheme:dark\)\{:root:not\(\[data-theme=\"light\"\]\)\{\s*--bg:")
        self.assertIn(':root[data-theme="dark"]{', p)
        self.assertRegex(p, r"body\{[^}]*background:var\(--bg\)")
        self.assertRegex(p, r"figure\.fig \.imgbox\{[^}]*background:#fff")   # plots stay on white
        self.assertRegex(p, r"@media print\{[^@]*\.toc[^}]*display:none")

    def test_header_tiles_callouts_and_figure_cards(self):
        p = self.page
        self.assertIn("<b>Submission</b>PROT_0001", p)
        self.assertIn('class="tiles"', p)
        self.assertIn("Some values are inferred, not measured", p)          # 60% inferred -> warning
        self.assertRegex(p, r'callout warning[^>]*>.*?Some values are inferred', )
        self.assertRegex(p, r'id="data-quality-notes">Data Quality Notes</h2><div class="callout warning')
        self.assertRegex(p, r'id="expert-review-notes">Expert Review Notes</h2><div class="callout critical')
        self.assertIn('class="imgbox"', p)
        self.assertIn('id="lb"', p)                                          # the lightbox
        self.assertIn('<ol start="2">', p)                                   # numbering survives
        self.assertIn('<td class="num">-0.3</td>', p)                       # numbers right-aligned
        self.assertIn("<p class=\"subtitle\"><em>Standfirst.</em></p>", p)


if __name__ == "__main__":
    unittest.main(verbosity=2)
