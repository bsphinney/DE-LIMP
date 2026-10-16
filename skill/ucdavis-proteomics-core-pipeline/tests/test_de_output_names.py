#!/usr/bin/env python3
"""A delivery without differential expression (core_submission.py deliver --mode no-de) refuses
the DE's output, and report_files.is_de_output is the one definition of it: the names run_de.R
(with contaminants.R's tables), run_sets.R and make_figures.R write. These tests read the R
scripts, so an output added there and not to report_files fails here instead of slipping into a
delivery that says it has no DE. stdlib only.
"""
import os
import re
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import report_files as rf  # noqa: E402

NAME = r"[A-Za-z][\w.-]*\.(?:csv|png|rds|json|txt|R)"


def source(script):
    with open(os.path.join(SCRIPTS, script), encoding="utf-8") as fh:
        return fh.read()


class DeOutputNames(unittest.TestCase):
    def test_every_file_the_de_scripts_write_is_de_output(self):
        for script, d in (("run_de.R", "tables"), ("contaminants.R", "tables"),
                          ("run_sets.R", "tables"), ("make_figures.R", "figures"),
                          ("run_silac.R", "tables")):
            src = source(script)
            # file.path(outdir, "x.csv"), and names given through a variable
            names = re.findall(r'file\.path\((?:outdir|figdir),\s*"(' + NAME + r')"\)', src)
            names += re.findall(r'(?:\bfile|set_inputs_file|members_file|share_table|removed_table)'
                                r'\s*(?:=|<-)\s*"(' + NAME + r')"', src)
            # one per contrast: sprintf("volcano_%s.png", ...)
            per_contrast = re.findall(r'sprintf\("([A-Za-z][\w-]*_)%s[^"]*\.(csv|png)"', src)
            self.assertTrue(names or per_contrast, f"{script}: no output names found -- has it moved?")
            for n in names:
                self.assertTrue(rf.is_de_output(f"{d}/{n}"), f"{script} writes {n}")
            for prefix, ext in per_contrast:
                self.assertTrue(rf.is_de_output(f"{d}/{prefix}A.vs.B.{ext}"),
                                f"{script} writes {prefix}<contrast>.{ext}")
        # run_silac.R's figures: named by a variable (f <- "x.png") and drawn into figdir
        for n in re.findall(r'\bf <- "(' + NAME + r')"', source("run_silac.R")):
            self.assertTrue(rf.is_de_output(f"figures/{n}"), f"run_silac.R draws {n}")
        # make_figures.R removes, as its own, the per-contrast p-value plots of sessions before 2.8.0
        self.assertIn("pvalue_.+", source("make_figures.R"))
        self.assertTrue(rf.is_de_output("figures/pvalue_A.vs.B.png"))

    def test_the_de_record_is_de_output(self):
        self.assertTrue(rf.is_de_output("tables/" + rf.DE_RECORD))
        self.assertTrue(rf.is_de_output("tables\\DE_dpc_A.vs.B.csv"))

    def test_a_qc_report_session_is_not_de_output(self):
        for rel in ("Analysis_Report.html", "Analysis_Report.md", "methods.md", "AUDIT.md",
                    "SAMPLE_QUALITY.md", "OUTPUT_FILES.md", "tables/qc_pooled_standard.csv",
                    "figures/qc_standard_trend.png", "reproducibility/REPRODUCE.md",
                    "search/report.parquet", "search/search_provenance.json", "podcast/podcast.m4a",
                    # a DE name outside tables/ and figures/ is something else
                    "methods.txt", "reproducibility/methods.txt"):
            self.assertFalse(rf.is_de_output(rel), rel)


if __name__ == "__main__":
    unittest.main()
