#!/usr/bin/env python3
"""
A SILAC session through the rest of the pipeline: the standard report, the Methods, the audit,
the analysis brief, the delivery gate, the reproducibility bundle, the set tests and the figures.

A triple-SILAC client report had to be made outside the skill (2026-10-14), so it lacked the
Submission section, the QC panels, the audit and the Core's survey. Every consumer below reads
the SILAC analysis through its self-describing record (run_silac.R's de_provenance.json:
pipeline_id silac_ratio, methods_paragraph, normalization_check not_applicable, checks,
report_sections, glance_tiles, brief_notes, replay, figures.written_by) -- none of them tests for
SILAC by name (DE-LIMP rule 1). The session is synthetic (tests/silac_fixture.py), its submission
a placeholder.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import silac_fixture as fx  # noqa: E402
from test_silac_analysis import HAVE_R, NEEDS, RSCRIPT  # noqa: E402

PROT = "PROT_0807"


def run(args, ok=True, **kw):
    p = subprocess.run(args, capture_output=True, text=True, timeout=600, **kw)
    if ok:
        assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    return p


@unittest.skipUnless(HAVE_R, NEEDS)
class SilacSession(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        s = cls.s = cls.tmp.name
        for d in ("input/wf", "output/search", "output/tables", "output/figures", "logs"):
            os.makedirs(os.path.join(s, d))
        cls.report = fx.make_report(os.path.join(s, "output", "search", "report.parquet"))
        cls.design = fx.write_design(os.path.join(s, "input", "silac_design.csv"))
        cls.cfg = os.path.join(s, "input", "wf", "params.cfg")
        run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"), "--engine", "diann",
             "--acquisition", "DIA", "--instrument", "timsTOF HT", "--precursor-mz-range", "299.5",
             "1200.5", "--labels", "silac3", "--label-preset", "K4R6,K8R10", "--channel-norm",
             "spec", "--out", cls.cfg])
        with open(os.path.join(s, "output", "search", "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "version": "2.7.0", "params_file": cls.cfg,
                       "search_mode": "two_job_libfree",
                       "result": {"report": cls.report}}, fh)
        cls.tables = os.path.join(s, "output", "tables")
        run([RSCRIPT, os.path.join(SCRIPTS, "run_silac.R"), "--input", cls.report,
             "--design", cls.design, "--outdir", cls.tables, "--params", cls.cfg])
        import submission_report
        submission_report.attach(s, {"internal_id": PROT})
        # the audit (AUDIT.md becomes the report's "Audit & caveats")
        cls.audit = run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                         "--out", os.path.join(s, "output", "AUDIT.md"), "--de-dir", cls.tables],
                        cwd=s)
        cls.audit_json = json.load(open(os.path.join(s, "output", "AUDIT.json")))
        with open(os.path.join(s, "output", "AI_Analysis_Report.md"), "w") as fh:
            fh.write("# Triple SILAC example\n\n*Three groups, four replicates.*\n\n"
                     "## Key Findings\n\n![KO vs WT](figures/volcano_KO.WT.png)\n\n"
                     "## Data Quality Notes\n\n![Channels](figures/silac_channel_composition.png)\n\n"
                     "WT dominates the signal.\n")
        out = os.path.join(s, "output", "Analysis_Report.html")
        cls.html_run = run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", s, "--out", out, "--no-pdf"])
        cls.page = open(out, encoding="utf-8").read()
        cls.md = open(os.path.join(s, "output", "Analysis_Report.md"), encoding="utf-8").read()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    # ---- the report ---------------------------------------------------------------------------
    def test_report_is_the_standard_one(self):
        p = self.page
        self.assertIn(">Submission<", p)
        self.assertIn(f"prot={PROT}", p, "the Core's survey")
        self.assertIn("Audit &amp; caveats", p)
        self.assertIn("silac_channel_balance", p)
        for title in ("SILAC channels and mixing", "Channel-specific proteins",
                      "Top proteins per contrast", "Results at a glance"):
            self.assertIn(title, p)
        self.assertIn("replicates (runs)", p)
        self.assertIn("channels compared", p)
        self.assertIn("SILAC ratios + limma", p)
        self.assertIn("KO vs WT", p)
        self.assertIn("detected 4/4 KO, 4/4 WT", p)
        n_img = len(re.findall(r'<img[^>]+src="data:image/png;base64,', p))
        self.assertGreaterEqual(n_img, 3, "two referenced figures + the p-value appendix")
        self.assertIn("Appendix: p-value calibration", p)
        self.assertIn("SILAC channels and mixing", self.md)

    def test_glance_counts_are_the_tables(self):
        # 20 of 52 tested significant in KO-WT: 10 up, 10 down (P40 fails Global.PG.Q.Value)
        text = re.sub(r"<[^>]+>", " ", self.page)
        self.assertRegex(text, r"20\s+KO vs WT\s+&#9650;&thinsp;10 up &nbsp;&#9660;&thinsp;10 down")
        self.assertRegex(text, r"0\s+DKO vs KO\s+&#9650;&thinsp;0 up")
        self.assertRegex(text, r"52\s+proteins tested")

    # ---- Methods ------------------------------------------------------------------------------
    def test_methods(self):
        out = os.path.join(self.s, "output", "methods.md")
        run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"), "--raw",
             os.path.join(self.s, "nowhere", "*.d"), "--instrument", "timsTOF HT",
             "--acquisition", "DIA", "--de-dir", self.tables, "--params", self.cfg,
             "--search-prov", os.path.join(self.s, "output", "search", "search_provenance.json"),
             "--out", out])
        m = open(out, encoding="utf-8").read()
        self.assertIn("## SILAC ratio analysis", m)
        self.assertNotIn("## Differential expression", m)
        for want in ("plexDIA", "Lys +4.025107 Da (2H(4), UniMod:481)",
                     "Arg +10.008269 Da (13C(6)15N(4), UniMod:267)",
                     "--channel-spec-norm", "median-centred", "robust = TRUE",
                     "| SILAC channels (plexDIA) | L K0R0; M K4R6; H K8R10 |",
                     "| Channel normalisation | --channel-spec-norm |"):
            self.assertIn(want, m)
        self.assertNotIn("Fixed modifications: Carbamidomethyl (C); SILAC", m)
        self.assertNotIn("| Fixed modification | SILAC", m)

    def test_methods_tag_a_defaulted_normalisation(self):
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "params.cfg")
            run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"), "--engine", "diann",
                 "--acquisition", "DIA", "--instrument", "timsTOF HT", "--precursor-mz-range",
                 "299.5", "1200.5", "--labels", "silac2", "--label-preset", "K8R10", "--out", cfg])
            sp = os.path.join(d, "search_provenance.json")
            json.dump({"engine": "diann", "version": "2.7.0", "params_file": cfg}, open(sp, "w"))
            out = os.path.join(d, "methods.md")
            run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"), "--raw",
                 os.path.join(d, "x.d"), "--instrument", "timsTOF HT", "--acquisition", "DIA",
                 "--params", cfg, "--search-prov", sp, "--out", out])
            m = open(out, encoding="utf-8").read()
            self.assertIn("--channel-spec-norm (independent samples in each channel", m)
            self.assertIn("DEFAULT, not user-confirmed", m)

    # ---- the audit, the brief, the gate --------------------------------------------------------
    def test_audit_quotes_the_pipeline(self):
        by = {}
        for f in self.audit_json["findings"]:
            by.setdefault(f["check"], []).append(f)
        self.assertEqual(by["normalization_check"][0]["status"], "PASS")
        self.assertIn("Not applicable", by["normalization_check"][0]["message"])
        self.assertEqual(by["silac_channel_balance"][0]["status"], "WARN")
        self.assertIn("silac_label_incorporation", by)
        if "normalisation_stability" in by:
            self.assertEqual(by["normalisation_stability"][0]["status"], "INFO")

    def test_brief(self):
        out = os.path.join(self.s, "brief.md")
        run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"), "--out", out,
             "--de-dir", self.tables, "--figures-dir", os.path.join(self.s, "output", "figures")])
        b = open(out, encoding="utf-8").read()
        self.assertIn("not applicable to this pipeline", b)
        self.assertNotIn("NOT RUN", b)
        self.assertIn("MEDIAN-CENTRED", b)
        self.assertIn("SILAC_channel_specific_detections.csv", b)
        self.assertIn("volcano_KO.WT.png", b)

    def test_delivery_gate_passes_on_the_records_own_reason(self):
        import normalization_check
        g = normalization_check.delivery_gate(self.s)
        self.assertTrue(g["proceed"], g)
        self.assertEqual(g["status"], "not_applicable")

    def test_a_not_applicable_without_a_reason_does_not_pass(self):
        import normalization_check
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "output", "tables"))
            json.dump({"normalization_check": {"status": "not_applicable"}},
                      open(os.path.join(d, "output", "tables", "de_provenance.json"), "w"))
            self.assertFalse(normalization_check.delivery_gate(d)["proceed"])

    def test_sdrf_has_a_row_per_run_and_channel_with_the_pride_label(self):
        import make_deposit as md
        with open(os.path.join(self.s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw\n" + "".join(f"/data/{r}.d\n" for r in fx.RUNS))
        try:
            cols, rows, sources, _ = md.build_sdrf(md.gather(self.s))
        finally:
            os.unlink(os.path.join(self.s, "input", "raw_files.txt"))
        c = {name: i for i, name in reversed(list(enumerate(cols)))}
        self.assertEqual(len(rows), 12, "4 runs x 3 channels")
        first = [r for r in rows if r[c["assay name"]] == "run_A1"]
        self.assertEqual([(r[c["source name"]], r[c["comment[label]"]], r[-1]) for r in first],
                         [("Rep1_WT", "SILAC light", "WT"), ("Rep1_KO", "SILAC medium", "KO"),
                          ("Rep1_DKO", "SILAC heavy", "DKO")])
        self.assertEqual({r[c["comment[data file]"]] for r in first}, {"run_A1.d"})
        keys = {(r[c["source name"]], r[c["assay name"]], r[c["comment[label]"]]) for r in rows}
        self.assertEqual(len(keys), len(rows), "source + assay + label must be unique (spec)")
        self.assertEqual(sorted(r[c["characteristics[biological replicate]"]] for r in rows
                                if r[-1] == "KO"), ["1", "2", "3", "4"])
        label_src = next(why for col, _, why in sources if col == "comment[label]")
        self.assertIn("PRIDE:0000327", label_src)

    def test_agents_md_reads_the_pipelines_own_columns(self):
        run([sys.executable, os.path.join(SCRIPTS, "session.py"), "docs", "--dir", self.s])
        a = open(os.path.join(self.s, "AGENTS.md"), encoding="utf-8").read()
        self.assertIn("- `logFC` — the mean median-centred log2 ratio", a)
        self.assertIn("- `uncentred_Rep2` — that replicate's log2 ratio before centring", a)
        self.assertNotIn("no description recorded", a.split("## Table columns")[1].split("##")[0])
        r = open(os.path.join(self.s, "README.md"), encoding="utf-8").read()
        self.assertIn("one `DE_silac_*.csv` per comparison", r)
        self.assertNotIn("Expression_Matrix.csv", r)

    # ---- reproducibility, set tests, figures ---------------------------------------------------
    def test_reproduce_sh_replays_the_silac_analysis(self):
        out = os.path.join(self.s, "repro")
        run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir", out,
             "--params", self.cfg, "--de-dir", self.tables, "--engine", "diann",
             "--acquisition", "DIA", "--instrument", "timsTOF HT", "--report", self.report])
        sh = open(os.path.join(out, "reproduce.sh")).read()
        self.assertIn('Rscript "$SKILL/scripts/run_silac.R" --input ./search_out/report.parquet', sh)
        self.assertIn("--design inputs/silac_design.csv", sh)
        self.assertIn("--params inputs/params.cfg", sh)
        self.assertIn("--comparisons KO-WT,DKO-WT,DKO-KO", sh)
        self.assertNotIn("run_de.R", sh)
        self.assertIn("--params ./inputs/params.cfg", sh, "the search replays the recorded labels")
        self.assertTrue(os.path.exists(os.path.join(out, "inputs", "silac_design.csv")))
        n = subprocess.run(["bash", "-n", os.path.join(out, "reproduce.sh")],
                           capture_output=True, text=True)
        self.assertEqual(n.returncode, 0, n.stderr)

    def test_set_tests_say_why_not(self):
        p = run([RSCRIPT, os.path.join(SCRIPTS, "run_sets.R"), "--de-dir", self.tables], ok=False)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("no protein-set tests for this analysis", p.stderr)

    def test_make_figures_leaves_the_silac_figures_alone(self):
        figs = os.path.join(self.s, "output", "figures")
        before = sorted(os.listdir(figs))
        fj = open(os.path.join(figs, "figures.json")).read()
        p = run([RSCRIPT, os.path.join(SCRIPTS, "make_figures.R"), "--de-dir", self.tables,
                 "--outdir", figs])
        self.assertIn("run_silac.R draws this analysis's figures", p.stderr + p.stdout)
        self.assertEqual(sorted(os.listdir(figs)), before)
        self.assertEqual(open(os.path.join(figs, "figures.json")).read(), fj)

    def test_sample_quality_reads_the_sample_matrix(self):
        out = os.path.join(self.s, "SAMPLE_QUALITY.md")
        run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"), "--matrix",
             os.path.join(self.tables, "SILAC_sample_matrix.csv"), "--conditions",
             os.path.join(self.tables, "SILAC_samples.csv"), "--out", out])
        self.assertTrue(os.path.getsize(out) > 0)


if __name__ == "__main__":
    unittest.main()
