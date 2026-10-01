#!/usr/bin/env python3
"""OUTPUT_FILES.md (make_report.py) describes every file the skill writes.

2.8.0 release review: the detection and contaminant tables, the DE-LIMP session, the audit and
sample-quality records, the CoreOmics submission, the deposit package and DIFFERENCES.md all came
out as "unrecognized output"; run_manifest.json and reproduce.sh were still described as the
removed workflow registry; and the Figures category was missing from CATEGORY_ORDER, so every
figure was dropped from the catalog without a word.
stdlib only.
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

import make_report      # noqa: E402

SKILL_FILES = [
    # run_de.R / contaminants.R / make_figures.R
    "Detection_Matrix.csv", "QC_detected_vs_inferred.csv", "QC_contaminant_share.csv",
    "contaminants_removed.csv", "DE-LIMP_session.rds", "sample_labels.csv",
    "qc_pvalue_panel.png",
    # audit_results.py / sample_quality.py
    "AUDIT.md", "AUDIT.json", "SAMPLE_QUALITY.md", "SAMPLE_QUALITY.json",
    # sage_lfq_check.py (via run_search.py)
    "sage_lfq_check.json", "sage_adapt.json",
    # submission_report.py
    "submission.json", "samples.tsv", "session.json",
    # make_deposit.py
    "HOW_TO_SUBMIT.md", "HOW_TO_SUBMIT.html", "sdrf.tsv", "protocols.txt",
    "files_to_upload.tsv", "prepare_upload.sbatch", "methods_complete_draft.md",
    # session.py
    "DIFFERENCES.md", "README.html", "README.md", "AGENTS.md", "raw_files.txt",
    # make_podcast.py
    "podcast.m4a", "podcast_script.md", "transcript.html", "podcast.json", "check.txt",
    "verify.txt", "verify_transcript.txt", "Analysis_Report_with_audio.html",
    # save_transcript.py / log_decision.py
    "conversation.md", "11111111-2222-4333-8444-555555555555.jsonl", "decisions.md",
]


class Catalog(unittest.TestCase):
    def test_the_shareable_reports_name_is_defined_once(self):
        # review of f18d95e: the CATALOG builds it from SHARE_NAME, which lives in the
        # import-free report_files.py that make_podcast and core_submission read too
        import core_submission
        import make_podcast
        import report_files
        self.assertEqual(make_report.describe(report_files.SHARE_NAME)[0], "Analysis report")
        self.assertIs(make_podcast.SHARE_NAME, report_files.SHARE_NAME)
        self.assertIs(core_submission.SHARE_FILE, report_files.SHARE_NAME)
        with open(make_report.__file__, encoding="utf-8") as fh:
            self.assertNotIn(report_files.SHARE_NAME.split(".")[0], fh.read())

    def test_a_partial_scripts_copy_still_catalogues(self):
        # verification of cba30b0: make_report imported make_podcast (and notify_slack) at load
        with tempfile.TemporaryDirectory() as d:
            lone = os.path.join(d, "scripts")
            os.makedirs(lone)
            for f in ("make_report.py", "report_files.py", "scratch_files.py"):
                shutil.copy(os.path.join(SCRIPTS, f), lone)
            out = os.path.join(d, "OUTPUT_FILES.md")
            with open(os.path.join(d, "Analysis_Report_with_audio.html"), "w") as fh:
                fh.write("x")
            r = subprocess.run([sys.executable, os.path.join(lone, "make_report.py"), "--out", out,
                                "--extra", os.path.join(d, "Analysis_Report_with_audio.html"),
                                "--root", d], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out, encoding="utf-8") as fh:
                self.assertIn("The report with the podcast's audio and transcript built in",
                              fh.read())

    def test_every_file_the_skill_writes_is_described(self):
        for name in SKILL_FILES:
            with self.subTest(name):
                cat, desc = make_report.describe(name)
                self.assertNotIn("unrecognized", desc)
                self.assertNotEqual(cat, "Other")

    def test_the_specific_figure_wins_over_the_generic_png(self):
        self.assertIn("p-value", make_report.describe("qc_pvalue_panel.png")[1])
        self.assertIn("figure", make_report.describe("volcano_A-B.png")[1])

    def test_no_description_names_the_removed_registry(self):
        for pat, _, desc in make_report.CATALOG:
            with self.subTest(pat):
                self.assertNotIn("registry commit", desc)
                self.assertNotIn("pinned workflow", desc)
                self.assertNotIn("Cont_", desc, "the tag is fetch_fasta.CONT_TAG's to spell")

    def test_the_predicted_library_is_described_as_left_out_of_the_zip(self):
        self.assertIn("Left out of the session zip",
                      make_report.describe("step1.predicted.speclib")[1])
        self.assertNotIn("predicted", make_report.describe("custom.speclib")[1])

    def test_every_category_is_printed_figures_included(self):
        self.assertEqual({c for _, c, _ in make_report.CATALOG} - set(make_report.CATEGORY_ORDER),
                         set())
        with tempfile.TemporaryDirectory() as d:
            figs = os.path.join(d, "figures")
            os.makedirs(figs)
            for n in ("volcano.png", "figures.json"):
                open(os.path.join(figs, n), "w").close()
            out = os.path.join(d, "OUTPUT_FILES.md")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--out",
                                out, "--extra", figs, "--root", d], capture_output=True,
                               text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out, encoding="utf-8") as fh:
                text = fh.read()
            self.assertIn("## Figures", text)
            self.assertIn("volcano.png", text)
            self.assertIn("figures.json", text)


class RealSessionLayout(unittest.TestCase):
    """PROT_0756 v2 (a 5-step DIA-NN chain, 30 runs): OUTPUT_FILES.md left out output/podcast/
    and marked 298 of the chain's files "unrecognized"."""

    CHAIN = ["quant_step2/r1.quant", "quant_step2/r2.quant", "quant_step4/r1.quant",
             "quant_step2_orig/r1.quant", "xic/t0.parquet", "xic/t0.log.txt",
             "xic/t0.stats.tsv", "xic/t0_xic/r1.xic.parquet",
             "xic/t0_xic/r1.ms1_mobilogram.parquet", "empirical.parquet",
             "step3_assembly.parquet", "step3_assembly.log.txt", "step3_assembly.stats.tsv",
             "step1.log.txt", "file_list.txt", "parallel_input_files.txt", "jobs.txt",
             "fran_deposit.json", "window.json", "window.txt", "step1_libpred.sbatch",
             "step5_report.sbatch", "submit.sh", "report.parquet", "report.pg_matrix.tsv",
             "report.pr_matrix.tsv", "report.unique_genes_matrix.tsv", "report.manifest.txt",
             "report.protein_description.tsv"]
    PODCAST = ["podcast.m4a", "transcript.html", "podcast_script.md", "check.txt", "verify.txt",
               "verify_transcript.txt", "podcast.json"]

    def build(self, d):
        out = os.path.join(d, "output")
        for rel in self.CHAIN:
            path = os.path.join(out, "search", rel)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as fh:
                fh.write("x" * 10)
        for n in self.PODCAST + ["podcast.wav", ".cache/chunk_001.wav"]:
            path = os.path.join(out, "podcast", n)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            open(path, "w").close()
        md = os.path.join(out, "OUTPUT_FILES.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--out", md,
                            "--search-out", os.path.join(out, "search"), "--root", d],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(md, encoding="utf-8") as fh:
            return fh.read(), json.loads(r.stdout)

    def test_nothing_the_chain_writes_is_unrecognized(self):
        with tempfile.TemporaryDirectory() as d:
            text, res = self.build(d)
            self.assertEqual(res["unrecognized"], 0, text)
            self.assertNotIn("unrecognized", text)

    def test_the_per_run_files_are_one_line_per_kind(self):
        with tempfile.TemporaryDirectory() as d:
            text, res = self.build(d)
            internals = text.split("## Search internals")[1].split("\n## ")[0]
            for line, n in (("output/search/xic/", 5), ("output/search/quant_step2/*.quant", 2),
                            ("output/search/quant_step4/*.quant", 1),
                            ("output/search/quant_step2_orig/*.quant", 1)):
                self.assertRegex(internals, rf"\| `{re.escape(line)}` \| {10 * n} B in {n} "
                                            rf"files? \|")
            self.assertNotIn("t0_xic/r1.xic.parquet", text)
            self.assertNotIn("r1.quant`", text)
            self.assertIn("The empirical spectral library", text)
            self.assertIn("| `output/search/step3_assembly.parquet` |", internals)
            # the header counts every file, the collapsed ones included
            self.assertIn(f"This run produced {len(self.CHAIN) + len(self.PODCAST)} file(s)", text)
            self.assertEqual(res["n_files"], len(self.CHAIN) + len(self.PODCAST))

    def test_the_podcast_beside_output_files_is_listed(self):
        with tempfile.TemporaryDirectory() as d:
            text, _ = self.build(d)
            for n in self.PODCAST:
                self.assertIn(f"| `output/podcast/{n}` |", text)
            self.assertNotIn("podcast.wav", text)            # scratch: beside podcast.m4a
            self.assertNotIn(".cache", text)


if __name__ == "__main__":
    unittest.main(verbosity=2)
