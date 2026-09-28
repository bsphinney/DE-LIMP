#!/usr/bin/env python3
"""2.8.0 release review: the analysis brief, the report page and the PDF step.

Blockers guarded here:
  * the brief no longer makes never-measured detection events the "most stable" or "top"
    proteins: CV from measured values only (linear scale), top by adj.P, Detection_Matrix
    attached, the multi-comparison claim softened;
  * PropObs is defined as limpa computes it (study-wide) and never used to tier hits; tiers
    come from per-group detection (the DE tables' Detected_<group>, else Detection_Matrix);
  * run_de.R's database_risk note is a FIXED callout on the page (HTML, .md, PDF), naming
    the proteins -- it never depends on the report writer.
Should-fixes: pull-down specificity + never-in-control names + Next steps; top tables with an
n.s. divider and Ig/contaminant flags; rule-2 tags when nothing recorded the rule; figures
not in this run's figures.json (stale / failed) are notes; the stale PDF is retired; the
inferred note only for dpc, in the record's own words; one contrast form + direction line;
the p-value definition; a scroll hint for wide tables.
"""
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import html_to_pdf as h2p  # noqa: E402

RUNS = {"Bait": ["b1", "b2", "b3"], "IgG": ["c1", "c2", "c3"]}
POLICY = "Missing precursors modelled via the detection probability curve; not imputed, not dropped."
# (protein, gene, logFC, adj.P, detection per run b1..b3 c1..c3)
ROWS = [("P1", "Jph3", 9.1, 1e-9, [3, 3, 3, 0, 0, 0]),      # never measured in IgG
        ("P2", "Ryr2", 5.0, 2e-8, [4, 4, 4, 2, 2, 1]),      # measured in both
        ("P3", "Ighg2c", 6.0, 3e-8, [5, 5, 5, 0, 0, 0]),    # the antibody's own chain
        ("Cont_P02070", "HBB", 4.0, 4e-7, [2, 2, 2, 0, 0, 0]),
        ("P5", "Stim2", 3.0, 0.004, [1, 0, 0, 0, 0, 0]),    # mostly unmeasured
        ("P6", "Actb", 0.2, 0.3, [9, 9, 9, 9, 9, 9]),       # n.s.
        ("P7", "Gapdh", -0.1, 0.8, [9, 9, 9, 9, 9, 9])]


def session(root, prov=None, detmat=True, det_cols=False, figures_json=None, report=True,
            removed=None):
    out = os.path.join(root, "output")
    tables, figs = os.path.join(out, "tables"), os.path.join(out, "figures")
    for d in (tables, figs, os.path.join(root, "input")):
        os.makedirs(d, exist_ok=True)
    prov = prov if prov is not None else {
        "adjp": 0.05, "q_cutoff": 0.01, "logfc": 1.0, "pipeline_id": "dpc", "method": "dpc",
        "missing_policy": POLICY, "n_samples": 6, "groups": {"Bait": 3, "IgG": 3},
        "contrasts": ["Bait-IgG"], "detection_matrix": {"zero_means": "inferred"}}
    with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
        json.dump(prov, fh)
    head = ["Protein.Group", "Genes", "logFC", "P.Value", "adj.P.Val", "PropObs"]
    if det_cols:
        head += ["Detected_Bait", "Detected_IgG"]
    with open(os.path.join(tables, "DE_dpc_Bait.IgG.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        for p, g, lfc, adj, d in ROWS:
            row = [p, g, lfc, adj / 2, adj, 0.5]
            if det_cols:
                row += [f"{sum(x > 0 for x in d[:3])}/3", f"{sum(x > 0 for x in d[3:])}/3"]
            w.writerow(row)
    if detmat:
        with open(os.path.join(tables, "Detection_Matrix.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Protein.Group"] + RUNS["Bait"] + RUNS["IgG"])
            for p, _, _, _, d in ROWS:
                w.writerow([p] + d)
    with open(os.path.join(tables, "QC_detected_vs_inferred.csv"), "w") as fh:
        fh.write("Sample,Group,Detected,Inferred,Total,PctDetected,PctInferred\n"
                 "b1,Bait,6,1,7,86,14\nc1,IgG,3,4,7,43,57\n")
    if removed:
        with open(os.path.join(tables, "contaminants_removed.csv"), "w") as fh:
            fh.write("Protein.Group,Genes,Contaminant.Group,Precursors,Contaminant.Precursors,"
                     "Removed.Entirely\n" + "".join(f"{p},{g},TRUE,3,3,TRUE\n" for p, g in removed))
    with open(os.path.join(root, "input", "conditions.csv"), "w") as fh:
        fh.write("File.Name,Group\n" + "".join(f"{r},{g}\n" for g, rs in RUNS.items() for r in rs))
    for f in ("volcano_Bait.IgG.png", "old_plot.png", "pvalue_Bait.IgG.png"):
        with open(os.path.join(figs, f), "wb") as fh:
            fh.write(b"\x89PNG\r\n\x1a\n" + f.encode())
    if figures_json is not None:
        with open(os.path.join(figs, "figures.json"), "w") as fh:
            json.dump(figures_json, fh)
    if report:
        with open(os.path.join(out, "AI_Analysis_Report.md"), "w") as fh:
            fh.write("# Pull-down\n\n## Overview\n\nText.\n\n![Volcano](figures/volcano_Bait.IgG.png)"
                     "\n\n![Old](figures/old_plot.png)\n\n![P](figures/pvalue_Bait.IgG.png)\n")
    return out, tables


def brief(root, **kw):
    out, tables = session(root, **kw)
    b = os.path.join(root, "brief.md")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"), "--out", b,
                        "--de-dir", tables, "--figures-dir", os.path.join(out, "figures"),
                        "--conditions", os.path.join(root, "input", "conditions.csv")],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(b) as fh:
        return fh.read()


def page(root, **kw):
    out, _ = session(root, **kw)
    html_out = os.path.join(out, "Analysis_Report.html")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                        "--session", root, "--out", html_out, "--no-pdf"],
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    with open(html_out, encoding="utf-8") as fh, open(html_out[:-5] + ".md", encoding="utf-8") as m:
        return json.loads(r.stdout), fh.read(), m.read(), r.stderr


class Brief(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name

    def tearDown(self):
        self._td.cleanup()

    def test_measured_only_cv_top_by_adjp_and_detection_matrix_attached(self):
        b = brief(self.root)
        self.assertIn("`tables/Detection_Matrix.csv`", b)
        self.assertIn("**0 = inferred**", b)
        self.assertIn("on the **linear** scale (2^value), from **measured values only**", b)
        self.assertIn("**ranked by adj.P.Val (ties by P.Value) — never by |logFC|**", b)
        self.assertNotIn("by fold change (use gene names", b)
        self.assertNotIn("most stable significant proteins (lowest CV)", b)
        self.assertNotIn("are the highest-confidence candidates", b)

    def test_propobs_is_study_wide_and_never_a_tier(self):
        b = brief(self.root)
        self.assertIn("rowMeans(n.observations) / NPeptides over **all runs of the study**", b)
        self.assertIn("**Never tier hits by PropObs", b)
        self.assertIn("count, from Detection_Matrix.csv + conditions.csv", b)
        b2 = brief(tempfile.mkdtemp(dir=self.root), det_cols=True)
        self.assertIn("the DE tables' `Detected_<group>` columns (k/n samples measured)", b2)
        self.assertIn("`Detected_Bait`, `Detected_IgG`", b2)

    def test_pulldown_specificity_never_in_control_and_next_steps(self):
        b = brief(self.root)
        self.assertIn("## Pull-down design — enrichment over a control", b)
        self.assertIn("### Specificity: bait-specific, shared, background", b)
        self.assertNotIn("### Cross-Comparison Biomarkers", b)
        self.assertNotIn("### High-Confidence", b)
        self.assertIn('Do NOT write a blanket sentence such as "the IgG value is a model-inferred '
                      'floor"', b)
        # Jph3 and Stim2 are significant, up and 0/3 in IgG; the Ig chain and the Cont_ entry
        # are background; Ryr2 was measured in IgG
        self.assertIn("Bait vs IgG: 2 significant enriched proteins never measured in IgG "
                      "(top by adj.P: Jph3, Stim2).", b)
        self.assertIn("### Next steps", b)
        self.assertIn("1. **Follow up first**", b)
        self.assertIn("these are the interactor candidates", b)
        self.assertIn("`README.html`", b)

    def test_nothing_recorded_is_said_not_assumed(self):
        b = brief(self.root, prov={"contrasts": ["Bait-IgG"], "groups": {"Bait": 3, "IgG": 3}})
        self.assertIn("Significance rule: **not recorded**", b)
        self.assertIn("0.05 (DEFAULT — not user-confirmed)", b)
        self.assertNotIn("q ≤ 0.01", b)
        self.assertNotIn("|log2FC| = 1", b)

    def test_pvalue_definition_contrast_names_and_direction(self):
        b = brief(self.root)
        self.assertIn("if the protein did not truly change, the probability of a difference at "
                      "least this large from measurement noise alone", b)
        self.assertIn("BH within each contrast", b)
        self.assertIn("`Bait-IgG` = **Bait vs IgG**", b)
        self.assertIn("Positive log2FC = higher in the first-named group of the contrast.", b)
        self.assertIn("**do not add a second per-contrast table**", b)


class Page(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name

    def tearDown(self):
        self._td.cleanup()

    def test_database_risk_is_a_fixed_callout_naming_the_proteins(self):
        prov = {"adjp": 0.05, "pipeline_id": "dpc", "missing_policy": POLICY,
                "contrasts": ["Bait-IgG"],
                "contaminants": {"policy": "removed", "removed": True,
                                 "removed_table": "contaminants_removed.csv",
                                 "database_risk": True,
                                 "database_note": "the search database was built before the "
                                                  "overlap check."}}
        res, html, md, _ = page(self.root, prov=prov,
                                removed=[("Cont_P60712", "ACTB"), ("Cont_P68103", "EEF1A1")])
        for doc in (html, md):
            self.assertIn("Real proteins may have been removed as contaminants", doc)
            self.assertIn("The search database was built before the overlap check.", doc)
            self.assertIn("ACTB, EEF1A1", doc)
        self.assertRegex(html, r'callout warning[^>]*>.*?Real proteins may have been removed')

    def test_inferred_note_only_for_dpc_in_the_records_words(self):
        _, html, md, _ = page(self.root)
        self.assertIn(POLICY.rstrip("."), html)
        self.assertIn("Some values are inferred, not measured", md)
        _, html2, _, _ = page(tempfile.mkdtemp(dir=self.root),
                              prov={"adjp": 0.05, "pipeline_id": "maxlfq", "method": "maxlfq",
                                    "contrasts": ["Bait-IgG"]})
        self.assertNotIn("Some values are inferred", html2)
        self.assertNotIn("limpa", html2)

    def test_glance_without_a_record_is_tagged(self):
        _, html, md, _ = page(self.root, prov={"contrasts": ["Bait-IgG"]})
        self.assertIn("(DEFAULT — not user-confirmed)", md)
        self.assertNotIn("the only rule the DE applied", md)
        self.assertIn("default — not recorded", html)

    def test_top_table_divider_and_flags(self):
        _, html, md, _ = page(self.root)
        self.assertIn('<tr class="divider"><td colspan="6">not significant below this line '
                      '(adj.P ≥ 0.05)</td></tr>', html)
        # the divider sits between the last significant row (Stim2) and the first n.s. (Actb)
        self.assertLess(html.index("Stim2"), html.index('class="divider"'))
        self.assertLess(html.index('class="divider"'), html.index("Actb"))
        self.assertIn("<b>Ig chain</b>", html)
        self.assertIn("<b>contaminant</b>", html)
        self.assertIn("| *not significant below this line (adj.P ≥ 0.05)* |", md)
        self.assertIn("**Ig chain**", md)
        self.assertIn("detected 3/3 Bait, 0/3 IgG", md)
        self.assertIn("Positive log2FC = higher in the first-named group", md)

    def test_only_figures_in_this_runs_figures_json(self):
        fj = {"figures": [{"file": "volcano_Bait.IgG.png", "caption": "v"}],
              "failed": [{"file": "pvalue_Bait.IgG.png", "reason": "ggplot error"}]}
        res, html, md, err = page(self.root, figures_json=fj)
        self.assertEqual(res["figures_embedded"], 1)
        self.assertEqual(res["figures_stale"], ["old_plot.png"])
        self.assertIn("figure not shown: old_plot.png is not in this run", html)
        self.assertIn("figure could not be drawn in this run: pvalue_Bait.IgG.png (ggplot error)",
                      html)
        self.assertIn("figure could not be drawn in this run", md)
        self.assertEqual(html.count("<img"), 1)
        self.assertIn("not in this run's figures.json", err)

    def test_wide_tables_get_a_scroll_hint(self):
        _, html, _, _ = page(self.root)
        self.assertIn("swipe the table sideways", html)
        self.assertIn(".scrollhint", html)


class StalePdf(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        d = self._td.name
        self.html, self.pdf = os.path.join(d, "Analysis_Report.html"), os.path.join(d, "Analysis_Report.pdf")

    def tearDown(self):
        self._td.cleanup()

    def _write(self, path, when):
        with open(path, "w") as fh:
            fh.write("x")
        os.utime(path, (when, when))

    def test_failed_print_retires_an_older_pdf(self):
        now = time.time()
        self._write(self.pdf, now - 100)
        self._write(self.html, now)
        with mock.patch.object(h2p, "convert", return_value=(False, "no browser")):
            status, note = h2p.print_report(self.html, self.pdf)
        self.assertEqual(status, "SKIPPED")
        self.assertFalse(os.path.exists(self.pdf))
        self.assertTrue(os.path.exists(self.pdf[:-4] + ".stale.pdf"))     # kept, never deleted
        self.assertIn("renamed Analysis_Report.stale.pdf", note)

    def test_failed_print_without_a_pdf_is_info(self):
        self._write(self.html, time.time())
        with mock.patch.object(h2p, "convert", return_value=(False, "no browser")):
            self.assertEqual(h2p.print_report(self.html, self.pdf)[0], "INFO")

    def test_success_is_ok(self):
        self._write(self.html, time.time())
        with mock.patch.object(h2p, "convert", return_value=(True, "3 pages")):
            self.assertEqual(h2p.print_report(self.html, self.pdf), ("OK", "3 pages"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
