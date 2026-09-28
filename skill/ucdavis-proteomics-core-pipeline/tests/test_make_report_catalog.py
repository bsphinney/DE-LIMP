#!/usr/bin/env python3
"""OUTPUT_FILES.md (make_report.py) describes every file the skill writes.

2.8.0 release review: the detection and contaminant tables, the DE-LIMP session, the audit and
sample-quality records, the CoreOmics submission, the deposit package and DIFFERENCES.md all came
out as "unrecognized output"; run_manifest.json and reproduce.sh were still described as the
removed workflow registry; and the Figures category was missing from CATEGORY_ORDER, so every
figure was dropped from the catalog without a word.
stdlib only.
"""
import os
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
    # submission_report.py
    "submission.json", "samples.tsv", "session.json",
    # make_deposit.py
    "HOW_TO_SUBMIT.md", "HOW_TO_SUBMIT.html", "sdrf.tsv", "protocols.txt",
    "files_to_upload.tsv", "prepare_upload.sbatch", "methods_complete_draft.md",
    # session.py
    "DIFFERENCES.md", "README.html", "README.md", "AGENTS.md", "raw_files.txt",
    # make_podcast.py
    "podcast.m4a", "podcast_script.md", "transcript.html", "podcast.json", "check.txt",
    "verify.txt", "verify_transcript.txt",
]


class Catalog(unittest.TestCase):
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
