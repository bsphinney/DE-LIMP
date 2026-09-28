#!/usr/bin/env python3
"""Analysis_Report.md: the HTML report's plain-text twin, for NotebookLM and other AI notebooks.

Brett, 2026-09-25: "can you make an md version of the output too so I can feed it to google
notebook llm". make_analysis_html.py assembles the page once and renders it twice, so the
Markdown must carry the same sections as the HTML, no HTML tags, every figure's caption as
text plus a data summary taken from the tables, and numbers that match the DE tables.
Stdlib only.
"""
import csv
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

GENES = ["Jph3", "Ryr2", "Stim2", "Cacna1a", "Fv1", "Kcnb1", "Actb", "Gapdh", "Tubb5", "Eno1",
         "Aldoa", "Pkm", "Hspa8", "Ywhaz", "Cfl1", "Pfn1", "Vim", "Des", "Myh9", "Tln1",
         "Vcl", "Flna", "Iqgap1", "Ahnak", "Spta1"]
SAMPLES = ["a1", "a2", "a3", "b1", "b2", "b3"]


def png(name):
    return b"\x89PNG\r\n\x1a\n" + name.encode() * 3


class MarkdownTwin(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._td = tempfile.TemporaryDirectory()
        s = cls.session = cls._td.name
        out = os.path.join(s, "output")
        figs, tables = os.path.join(out, "figures"), os.path.join(out, "tables")
        for d in (figs, tables, os.path.join(s, "input")):
            os.makedirs(d, exist_ok=True)
        for f in ("volcano_B.A.png", "pvalue_B.A.png", "pca.png"):
            with open(os.path.join(figs, f), "wb") as fh:
                fh.write(png(f))
        with open(os.path.join(figs, "figures.json"), "w") as fh:
            json.dump([{"file": "volcano_B.A.png", "caption": "Volcano caption from figures.json."},
                       {"file": "pvalue_B.A.png", "caption": "P-value caption."},
                       {"file": "pca.png", "caption": "PCA caption."}], fh)
        # DE table: adj.P rises with the row; alternate signs; P.Value half of adj.P
        cls.de = []
        for i, g in enumerate(GENES):
            p = 1e-6 * (3 ** i) if i < 12 else 0.2 + 0.03 * i
            lfc = (5.0 - 0.3 * i) * (1 if i % 3 else -1)
            cls.de.append({"Protein.Group": f"P{i:05d}", "Genes": g, "logFC": lfc,
                           "P.Value": min(1.0, p / 2), "adj.P.Val": min(1.0, p)})
        with open(os.path.join(tables, "DE_dpc_B.A.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, ["Protein.Group", "Genes", "logFC", "P.Value", "adj.P.Val"])
            w.writeheader()
            w.writerows(cls.de)
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump({"adjp": 0.05, "n_samples": 6, "groups": {"A": 3, "B": 3},
                       "contrasts": ["B-A"], "display_label": "DPC-Quant + limma (limpa)",
                       "pipeline_id": "dpc",
                       "detection_matrix": {"file": "Detection_Matrix.csv",
                                            "zero_means": "inferred"}}, fh)
        with open(os.path.join(tables, "QC_detected_vs_inferred.csv"), "w") as fh:
            fh.write("Sample,Group,Detected,Inferred,Total,PctDetected,PctInferred\n"
                     + "".join(f"{x},{'A' if x[0] == 'a' else 'B'},{20 + k},{5 - k},25,"
                               f"{(20 + k) * 4},{(5 - k) * 4}\n" for k, x in enumerate(SAMPLES)))
        with open(os.path.join(tables, "Detection_Matrix.csv"), "w") as fh:
            fh.write("Protein.Group," + ",".join(SAMPLES) + "\n")
            for i in range(len(GENES)):       # P00000: seen in all B, no A -> "detected 3/3 B, 0/3 A"
                fh.write(f"P{i:05d}," + ",".join(("0" if (i == 0 and x[0] == "a") else "2")
                                                 for x in SAMPLES) + "\n")
        with open(os.path.join(s, "input", "conditions.csv"), "w") as fh:
            fh.write("File.Name,Group\n" + "".join(f"{x},{'A' if x[0] == 'a' else 'B'}\n"
                                                   for x in SAMPLES))
        with open(os.path.join(out, "AI_Analysis_Report.md"), "w", encoding="utf-8") as fh:
            fh.write("# Twin study\n\n*One-line standfirst.*\n\n## Overview\n\nText <b>bold</b>.\n\n"
                     "![PCA of the samples](figures/pca.png)\n\n## Key findings\n\n"
                     "![Volcano — B vs A](figures/volcano_B.A.png)\n"
                     "![p-values — B vs A](figures/pvalue_B.A.png)\n- a list item\n\n"
                     "![Gone](figures/heatmap_top.png)\n\n"
                     "## Audit & caveats\n\n- ⚠️ one warning\n\n## Data Quality Notes\n\n"
                     "1. **Something.**\n   - detail\n")
        # the Submission section and header fact come from the attached record
        import submission_report
        submission_report.attach(s, {"internal_id": "PROT_0001", "organism": "Mouse"})
        cls.html_path = os.path.join(out, "Analysis_Report.html")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", s, "--out", cls.html_path, "--no-pdf"],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        cls.res = json.loads(r.stdout)
        with open(cls.html_path, encoding="utf-8") as fh:
            cls.html = fh.read()
        with open(os.path.join(out, "Analysis_Report.md"), encoding="utf-8") as fh:
            cls.md = fh.read()

    @classmethod
    def tearDownClass(cls):
        cls._td.cleanup()

    def test_written_beside_the_html(self):
        self.assertEqual(self.res["markdown_twin"],
                         os.path.join(self.session, "output", "Analysis_Report.md"))
        self.assertLess(len(self.md.encode("utf-8")), 1_000_000)

    def test_same_sections_as_the_html(self):
        html_h2 = [re.sub(r"<[^>]+>", "", t).replace("&amp;", "&")
                   for t in re.findall(r'<h2 id="[^"]+">(.*?)</h2>', self.html)]
        md_h2 = re.findall(r"^## (.*)$", self.md, re.M)
        self.assertEqual(md_h2, html_h2)
        self.assertIn("Results at a glance", md_h2)
        self.assertIn("Top proteins per contrast", md_h2)
        self.assertTrue(self.md.startswith("# Twin study\n"))
        self.assertIn("**Submission:** PROT_0001", self.md)
        self.assertEqual(md_h2[0], "Submission")          # the record, first, in both renderings
        self.assertIn("| Submission | PROT_0001 |", self.md)

    def test_no_html_tags(self):
        self.assertIsNone(re.search(r"</?[A-Za-z][A-Za-z0-9]*\b[^>]*>", self.md), self.md[:500])

    def test_figures_carry_their_captions_and_data(self):
        self.assertEqual(len(re.findall(r"^!\[Figure \d+\. ", self.md, re.M)), 3)
        self.assertEqual(self.html.count("<img"), 3)
        for cap in ("Volcano caption from figures.json.", "P-value caption.", "PCA caption."):
            self.assertIn(cap, self.md)
        self.assertIn("![Figure 2. Volcano — B vs A](figures/volcano_B.A.png)", self.md)
        self.assertIn("> **Note:** figure missing: heatmap_top.png", self.md)
        self.assertIn("figure missing: heatmap_top.png", self.html)

    def test_numbers_match_the_de_table(self):
        sig = [r for r in self.de if r["adj.P.Val"] < 0.05]
        up, dn = sum(r["logFC"] > 0 for r in sig), sum(r["logFC"] < 0 for r in sig)
        self.assertIn(f"| B vs A | {up + dn} | {up} | {dn} | {len(self.de)} |", self.md)
        self.assertIn(f"{up + dn} of {len(self.de)} proteins significant at adj.P < 0.05 "
                      f"({up} up, {dn} down", self.md)
        # volcano summary: the top 5 by adj.P, in order, with their values
        top5 = sorted(self.de, key=lambda r: r["adj.P.Val"])[:5]
        m = re.search(r"Top 5 by adj\.P: (.*?)\.\n", self.md)
        got = re.findall(r"(\w+) \(log2FC ([+-][\d.]+), adj\.P ([\d.e+-]+)\)", m.group(1))
        self.assertEqual([g for g, _, _ in got], [r["Genes"] for r in top5])
        for (g, lfc, p), r in zip(got, top5):
            self.assertAlmostEqual(float(lfc), r["logFC"], places=2)
            self.assertAlmostEqual(float(p) / r["adj.P.Val"], 1, delta=0.006)
        # p-value summary
        below = sum(r["P.Value"] < 0.05 for r in self.de)
        self.assertIn(f"{len(self.de)} p-values, {below} ", self.md)

    def test_top_protein_table(self):
        block = self.md.split("## Top proteins per contrast", 1)[1]
        rows = re.findall(r"^\| (P\d+) \| (\w+) \| ([+-][\d.]+) \| ([\d.e+-]+) \| (.*?) \|$",
                          block, re.M)
        self.assertEqual(len(rows), 20)
        want = sorted(self.de, key=lambda r: r["adj.P.Val"])[:20]
        self.assertEqual([r[0] for r in rows], [w["Protein.Group"] for w in want])
        for (prot, gene, lfc, p, det), w in zip(rows, want):
            self.assertEqual(gene, w["Genes"])
            self.assertAlmostEqual(float(lfc), w["logFC"], places=2)
            self.assertAlmostEqual(float(p) / w["adj.P.Val"], 1, delta=0.006)
        self.assertEqual(rows[0][4], "detected 3/3 B, 0/3 A")       # from Detection_Matrix.csv
        self.assertEqual(rows[1][4], "detected 3/3 B, 3/3 A")

    def test_callouts_are_labelled_blockquotes(self):
        self.assertRegex(self.md, r"## Data Quality Notes\n\n> \*\*Warning\.\*\*\n>\n> 1\. \*\*Something")
        self.assertIn("> **Note:** Some values are inferred, not measured.", self.md)  # <50%
        self.assertIn("Significant = adjusted p < 0.05", self.md)

    def test_standfirst_and_text_survive_without_tags(self):
        self.assertIn("*One-line standfirst.*", self.md)
        self.assertIn("Text bold.", self.md)


if __name__ == "__main__":
    unittest.main(verbosity=2)
