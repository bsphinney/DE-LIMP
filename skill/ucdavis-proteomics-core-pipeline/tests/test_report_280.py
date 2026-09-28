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
Then (figures.json 2.8.0, fix/280-figures; DE detection columns, fix/280-stats): the tiers are
run_de.R's Evidence categories, described in the record's own words; the top-protein plot is
not always a violin; the p-value panel is the page's appendix, not narrative; every figure
make_figures.R could not draw is a fixed callout.
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
PROV = {"adjp": 0.05, "q_cutoff": 0.01, "logfc": 1.0, "pipeline_id": "dpc", "method": "dpc",
        "missing_policy": POLICY, "n_samples": 6, "groups": {"Bait": 3, "IgG": 3},
        "contrasts": ["Bait-IgG"], "detection_matrix": {"zero_means": "inferred"}}
# run_de.R's de_provenance.json detection_matrix.de_columns (fix/280-stats 4d5bb69), verbatim
DE_COLUMNS_REC = {
    "Detected_<group>": "k/n: the runs of that group (n) in which the protein was measured (k; "
                        "Detection_Matrix.csv > 0, 0 = inferred)",
    "Evidence": "'measured in both' = measured in every run of every group in the contrast; "
                "'presence call' = never measured in at least one group, so the difference there "
                "is the model's, not a measurement; 'partly inferred' = otherwise"}
# (protein, gene, logFC, adj.P, detection per run b1..b3 c1..c3)
ROWS = [("P1", "Jph3", 9.1, 1e-9, [3, 3, 3, 0, 0, 0]),      # never measured in IgG
        ("P2", "Ryr2", 5.0, 2e-8, [4, 4, 4, 2, 2, 1]),      # measured in both
        ("P3", "Ighg2c", 6.0, 3e-8, [5, 5, 5, 0, 0, 0]),    # the antibody's own chain
        ("Cont_P02070", "HBB", 4.0, 4e-7, [2, 2, 2, 0, 0, 0]),
        ("P5", "Stim2", 3.0, 0.004, [1, 0, 0, 0, 0, 0]),    # mostly unmeasured
        ("P6", "Actb", 0.2, 0.3, [9, 9, 9, 9, 9, 9]),       # n.s.
        ("P7", "Gapdh", -0.1, 0.8, [9, 9, 9, 9, 9, 9])]


def evidence(d):
    """run_de.R's Evidence for one row's detection (b1..b3, c1..c3)."""
    k = [sum(x > 0 for x in d[:3]), sum(x > 0 for x in d[3:])]
    return ("presence call" if 0 in k else "measured in both" if k == [3, 3]
            else "partly inferred")


def session(root, prov=None, detmat=True, det_cols=False, figures_json=None, report=True,
            removed=None, evidence_col=False, extra_figs=(), omit=(), report_extra="",
            rows=None, prov_bytes=None):
    """rows: (protein, gene, logFC, adj.P, detection b1..b3 c1..c3[, raw P]) -- ROWS by default,
    raw P = adj.P / 2 unless given. prov_bytes: de_provenance.json exactly as these bytes."""
    rows = rows or ROWS
    out = os.path.join(root, "output")
    tables, figs = os.path.join(out, "tables"), os.path.join(out, "figures")
    for d in (tables, figs, os.path.join(root, "input")):
        os.makedirs(d, exist_ok=True)
    prov = prov if prov is not None else dict(PROV)
    with open(os.path.join(tables, "de_provenance.json"), "wb") as fh:
        fh.write(prov_bytes if prov_bytes is not None else json.dumps(prov).encode())
    head = ["Protein.Group", "Genes", "logFC", "P.Value", "adj.P.Val", "PropObs"]
    if det_cols:
        head += ["Detected_Bait", "Detected_IgG"] + (["Evidence"] if evidence_col else [])
    with open(os.path.join(tables, "DE_dpc_Bait.IgG.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        for p, g, lfc, adj, d, *raw in rows:
            row = [p, g, lfc, raw[0] if raw else adj / 2, adj, 0.5]
            if det_cols:
                row += [f"{sum(x > 0 for x in d[:3])}/3", f"{sum(x > 0 for x in d[3:])}/3"]
                row += [evidence(d)] if evidence_col else []
            w.writerow(row)
    if detmat:
        with open(os.path.join(tables, "Detection_Matrix.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Protein.Group"] + RUNS["Bait"] + RUNS["IgG"])
            for p, _, _, _, d, *_ in rows:
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
    for f in ("volcano_Bait.IgG.png", "old_plot.png", "pvalue_Bait.IgG.png") + tuple(extra_figs):
        if f in omit:
            continue
        with open(os.path.join(figs, f), "wb") as fh:
            fh.write(b"\x89PNG\r\n\x1a\n" + f.encode())
    if figures_json is not None:
        with open(os.path.join(figs, "figures.json"), "w") as fh:
            json.dump(figures_json, fh)
    if report:
        with open(os.path.join(out, "AI_Analysis_Report.md"), "w") as fh:
            fh.write("# Pull-down\n\n## Overview\n\nText.\n\n![Volcano](figures/volcano_Bait.IgG.png)"
                     "\n\n![Old](figures/old_plot.png)\n\n![P](figures/pvalue_Bait.IgG.png)\n"
                     + report_extra)
    return out, tables


def brief_run(root, env=None, **kw):
    """-> (brief text, the JSON analysis_prompt prints, stderr)."""
    out, tables = session(root, **kw)
    b = os.path.join(root, "brief.md")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"), "--out", b,
                        "--de-dir", tables, "--figures-dir", os.path.join(out, "figures"),
                        "--conditions", os.path.join(root, "input", "conditions.csv")],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    with open(b, encoding="utf-8") as fh:
        return fh.read(), json.loads(r.stdout), r.stderr


def brief(root, **kw):
    return brief_run(root, **kw)[0]


def page(root, env=None, **kw):
    out, _ = session(root, **kw)
    html_out = os.path.join(out, "Analysis_Report.html")
    r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                        "--session", root, "--out", html_out, "--no-pdf"],
                       capture_output=True, text=True, env=env)
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
        self.assertIn("from Detection_Matrix.csv + conditions.csv (count each group's runs", b)
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

    def test_tiers_are_run_de_evidence_quoted_from_the_record(self):
        prov = dict(PROV, detection_matrix={"zero_means": "inferred",
                                            "de_columns": DE_COLUMNS_REC})
        b = brief(self.root, prov=prov, det_cols=True, evidence_col=True)
        self.assertIn("`Detected_Bait`, `Detected_IgG` (k/n samples measured), `Evidence`.", b)
        self.assertIn(f"  - `Evidence`: {DE_COLUMNS_REC['Evidence']} (as recorded in "
                      "`de_provenance.json`).", b)
        self.assertIn(f"  - `Detected_<group>`: {DE_COLUMNS_REC['Detected_<group>']}", b)
        self.assertIn("from the DE tables' `Evidence` column:", b)
        # the record defines the categories: named here, not restated
        self.assertIn("  - **measured in both**: the fold change is a measured ratio;", b)
        self.assertIn("  - **presence call**: a detection event", b)
        self.assertIn("  - **partly inferred**: a measured ratio with gaps", b)
        self.assertNotIn("(measured in every run of every group in the contrast)", b)
        self.assertIn("1. **Follow up first** — significant and *measured in both*", b)
        self.assertIn("significant with a *presence call*", b)
        self.assertNotIn("≥ half the samples", b)
        # an older run without the column: the same three categories, defined here
        b0 = brief(tempfile.mkdtemp(dir=self.root))
        self.assertIn("sorted into the three categories run_de.R writes as `Evidence`", b0)
        self.assertIn("**measured in both** (measured in every run of every group in the "
                      "contrast)", b0)

    def test_top_protein_plots_and_the_pvalue_panel_goes_to_the_appendix(self):
        fj = {"figures": [{"file": "volcano_Bait.IgG.png", "type": "volcano", "caption": "V."},
                          {"file": "violin_top_Bait.IgG.png", "type": "violin",
                           "caption": "Top 12 proteins; points + mean bar."},
                          {"file": "qc_pvalue_panel.png", "type": "pvalue",
                           "caption": "Raw p-value distributions for all 1 contrast (appendix; a "
                                      "calibration check, not a result)."}],
              "failed": []}
        b = brief(self.root, figures_json=fj)
        self.assertIn("Call it the **top-protein plot**", b)
        self.assertNotIn("hollow points were not", b)
        self.assertIn("**Appendix figure — do not embed:** `figures/qc_pvalue_panel.png`", b)
        self.assertNotIn("- `figures/qc_pvalue_panel.png`", b)      # not in the embed list
        self.assertIn("*Appendix: p-value calibration*", b)
        self.assertNotIn("volcano + p-value figures", b)
        self.assertNotIn("Use the p-value distribution", b)
        # an older session's per-contrast histograms: together, in an appendix of the report
        b2 = brief(tempfile.mkdtemp(dir=self.root),
                   figures_json=[{"file": "pvalue_Bait.IgG.png", "type": "pvalue", "caption": "P."}])
        self.assertIn("- `figures/pvalue_Bait.IgG.png` (pvalue) — P.", b2)
        self.assertIn("embed them together in a final `## Appendix: p-value calibration` "
                      "section, never in the per-comparison sections", b2)


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

    def test_failed_figures_are_a_fixed_callout_and_the_panel_is_the_appendix(self):
        why = ("a PCA needs at least 3 samples and 5 varying proteins; this matrix has 2 samples "
               "and 20 proteins")
        fj = {"figures": [{"file": "volcano_Bait.IgG.png", "type": "volcano", "caption": "v"},
                          {"file": "qc_pvalue_panel.png", "type": "pvalue", "caption": "Panel."}],
              "failed": [{"file": "pca.png", "type": "pca", "reason": why}]}
        res, html, md, _ = page(self.root, figures_json=fj, extra_figs=("qc_pvalue_panel.png",),
                                omit=("pvalue_Bait.IgG.png",))
        self.assertEqual(res["figures_embedded"], 2)                 # volcano + the appendix
        for doc in (html, md):
            self.assertIn("1 figure could not be drawn in this run", doc)
        self.assertIn(f"<code>pca.png</code>: {why}", html)
        self.assertIn(f"`pca.png`: {why}", md)
        # the appendix: the page's last section, after the top-protein tables, with its summary
        self.assertLess(html.rindex("Top proteins per contrast"),
                        html.rindex("Appendix: p-value calibration"))
        self.assertIn("## Appendix: p-value calibration", md)
        self.assertIn("Share of raw p-values below 0.05 per contrast (about 5% when nothing "
                      "changed): Bait vs IgG 71%.", md)
        # a report written for per-contrast histograms says where they went
        self.assertIn("figure missing: pvalue_Bait.IgG.png (this run draws every contrast's "
                      "p-values in one panel: see Appendix: p-value calibration)", md)
        # a report that embeds the panel itself gets no second copy
        res2, html2, md2, _ = page(tempfile.mkdtemp(dir=self.root), figures_json=fj,
                                   extra_figs=("qc_pvalue_panel.png",),
                                   report_extra="\n## QC\n\n![Panel](figures/qc_pvalue_panel.png)\n")
        self.assertNotIn("## Appendix: p-value calibration", md2)
        self.assertEqual(md2.count("qc_pvalue_panel.png"), 1)

    def test_top_table_note_carries_run_de_evidence(self):
        _, _, md, _ = page(self.root, det_cols=True, evidence_col=True)
        self.assertIn("detected 3/3 Bait, 0/3 IgG — presence call", md)
        self.assertIn("detected 3/3 Bait, 3/3 IgG — measured in both", md)

    def test_wide_tables_get_a_scroll_hint(self):
        _, html, _, _ = page(self.root)
        self.assertIn("swipe the table sideways", html)
        self.assertIn(".scrollhint", html)


# Round 2 (verifier on 9dc5aa2)
TIES = [("PI", "Inf1", 9.0, 1e-6, [3, 3, 3, 0, 0, 0], 1e-8),     # presence call, big |logFC|
        ("PM", "Meas1", 1.0, 1e-6, [3, 3, 3, 3, 3, 3], 1e-9),    # same adj.P, smaller raw P
        ("PX", "Next1", 2.0, 5e-4, [3, 3, 3, 3, 3, 3], 2e-4)]
# A latin-1 locale: open() without an encoding mis-decodes (or, under cp1252, rejects) UTF-8.
LATIN1 = dict(os.environ, LC_ALL="en_US.ISO8859-1", LANG="en_US.ISO8859-1", PYTHONUTF8="0")
NOTE = "the search database — built before the overlap check — held 5 µg of bovine entries."


class TopTies(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name

    def tearDown(self):
        self._td.cleanup()

    def test_adjp_ties_break_by_raw_p_not_fold_change(self):
        import make_analysis_html as mah
        _, tables = session(self.root, rows=TIES)
        t = mah.Tables(tables, dict(PROV), 0.05, "de_provenance.json")
        self.assertEqual([r["Genes"] for r in t.top("Bait.IgG", 3)], ["Meas1", "Inf1", "Next1"])
        _, _, md, _ = page(tempfile.mkdtemp(dir=self.root), rows=TIES)
        top = md[md.index("## Top proteins per contrast"):]
        self.assertLess(top.index("Meas1"), top.index("Inf1"))

    def test_a_tie_without_a_raw_p_keeps_its_table_order(self):
        import make_analysis_html as mah
        _, tables = session(self.root, rows=[("PA", "A1", 5.0, 1e-3, [1] * 6, ""),
                                             ("PB", "B1", 9.0, 1e-3, [1] * 6, "")])
        t = mah.Tables(tables, dict(PROV), 0.05, "de_provenance.json")
        self.assertEqual([r["Genes"] for r in t.top("Bait.IgG", 2)], ["A1", "B1"])


class Encoding(unittest.TestCase):
    """Records are UTF-8 whatever the locale, and one that cannot be read is said -- on stderr,
    on the page and in the brief -- never a silent {} that drops what it held."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name

    def tearDown(self):
        self._td.cleanup()

    def prov(self, note=NOTE):
        return dict(PROV, contaminants={"policy": "removed", "removed": True,
                                        "database_risk": True, "database_note": note,
                                        "removed_table": "contaminants_removed.csv"})

    def test_utf8_records_read_right_under_a_latin1_locale(self):
        raw = json.dumps(self.prov(), ensure_ascii=False).encode("utf-8")
        res, html, md, err = page(self.root, env=LATIN1, prov_bytes=raw,
                                  removed=[("Cont_P60712", "ACTB")])
        for doc in (html, md):
            self.assertIn("Real proteins may have been removed as contaminants", doc)
            self.assertIn("The search database — built before the overlap check — held 5 µg", doc)
            self.assertNotIn("\u00e2\u0080\u0094", doc)                  # "—" read as latin-1
        self.assertEqual(res["records_unreadable"], {})
        self.assertNotIn("WARNING", err)
        b, out, _ = brief_run(tempfile.mkdtemp(dir=self.root), env=LATIN1, prov_bytes=raw)
        self.assertIn("built before the overlap check — held 5 µg", b)
        self.assertIn("≥2 contrasts", b)                              # written as UTF-8
        self.assertEqual(out["records_unreadable"], {})

    def test_a_record_that_is_not_utf8_is_loud_and_still_used(self):
        raw = json.dumps(self.prov(note="5 \u00b5g"), ensure_ascii=False).encode("latin-1")
        res, html, md, err = page(self.root, prov_bytes=raw, removed=[("Cont_P60712", "ACTB")])
        path = os.path.join(self.root, "output", "tables", "de_provenance.json")
        self.assertIn(f"WARNING: {path} is not UTF-8", err)
        self.assertEqual(list(res["records_unreadable"]), [path])
        for doc in (html, md):
            self.assertIn("1 record of this run could not be read", doc)
            self.assertIn("Real proteins may have been removed as contaminants", doc)  # kept

    def test_an_unparseable_record_is_said_not_dropped(self):
        raw = json.dumps(self.prov()).encode()[:-20]                  # truncated
        res, html, md, err = page(self.root, prov_bytes=raw)
        self.assertIn("de_provenance.json could not be read (JSONDecodeError", err)
        for doc in (html, md):
            self.assertIn("1 record of this run could not be read", doc)
            self.assertIn("missing from the page, not from the run", doc)
        self.assertIn("(DEFAULT — not user-confirmed)", md)            # the cutoff: tagged
        b, out, err2 = brief_run(tempfile.mkdtemp(dir=self.root), prov_bytes=raw)
        self.assertIn("## Records that could not be read — say so, do not fill in", b)
        self.assertIn("de_provenance.json` could not be read (JSONDecodeError", b)
        self.assertEqual(len(out["records_unreadable"]), 1)
        self.assertIn("WARNING", err2)


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
