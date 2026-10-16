#!/usr/bin/env python3
"""
A client deliverable never narrates the drafts the client never saw (client_text.py).

A staff report (PROT_0000, skill 2.11.2): after the step-8d expert review, the report's Expert
Review Notes read as a changelog ("checked the first draft", "two protein names were wrong ...
corrected"), the audit said an earlier normalisation reason was "corrected", and the podcast
said "the first draft of this report fell for it". The client never saw a draft, and
`make_podcast.py check` passed all of it. The phrase list lives in client_text.py, and both
the podcast check (FAIL) and make_analysis_html.py (a warning) read it.

Stdlib only, no network.
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import client_text  # noqa: E402
from test_make_podcast import SEG1, SEG2, Workspace, read, run, script_text, write  # noqa: E402

# Draft history: each must be caught.
HISTORY = (
    "The first draft of this report fell for it.",
    "We checked the first draft against the tables.",
    "Two protein names were wrong and have been corrected.",
    "Two protein names were wrong; both are now corrected.",
    "The categories were corrected.",
    "An earlier normalisation reason has since been corrected.",
    "In an earlier version, the volcano was mislabelled.",
    "A previous version of this report listed 12 proteins.",
    "Honest correction: the bait is not enriched.",
    "We corrected the protein names after the review.",
    "We got that wrong at first.",
    "The draft said the opposite.",
    "This version of the report fixes the legend.",
    "The report now says the controls were IgG.",
    "We fixed the mistake in the legend.",
    "Our earlier draft overstated it.",
    "The initial draft called it significant.",
    "That trap is easy to fall into, and we fell into it.",
    "The error was corrected before delivery.",
)
# Ordinary science: never caught.
SCIENCE = (
    "P-values were corrected for multiple testing (Benjamini-Hochberg).",
    "The p-values were corrected using the Benjamini-Hochberg procedure.",
    "Batch-corrected intensities are in the matrix.",
    "BH-corrected p-values below 0.05 were called significant.",
    "Corrected p-values are in the adj.P.Val column.",
    "We corrected for the batch effect by adding Batch to the design.",
    "The mass accuracy was fixed at 15 ppm.",
    "Carbamidomethylation was a fixed modification, and Batch a fixed effect.",
    "The signal fell for most proteins in the knockout.",
    "Proteins were annotated from a draft genome assembly.",
    "The first draft genome of this species was used for the database.",
    "These are previously reported interactors of the bait.",
    "As previously described, peptides were desalted.",
    "An earlier version of UniProt listed this protein under another name.",
    "The previous version of DIA-NN used a different scoring model.",
    "The first pass identified 62,000 precursors per run.",
    "The label on the tube was wrong, so the sheet lists ME03 twice.",
    "Expression now increases with temperature.",
    "Cells were drafted into the experiment in two batches.",
    "The mass error was corrected by recalibration.",
    "Mass errors were corrected using the lock mass.",
    "Use this table for your first draft of the manuscript.",
    "The draft manuscript is attached.",
    "Values were corrected (Benjamini-Hochberg).",
)


class Phrases(unittest.TestCase):
    def test_draft_history_is_caught(self):
        for text in HISTORY:
            with self.subTest(text=text):
                self.assertTrue(client_text.draft_history(text), text)

    def test_ordinary_science_is_not(self):
        for text in SCIENCE:
            with self.subTest(text=text):
                self.assertEqual(client_text.draft_history(text), [], text)

    def test_each_hit_is_a_readable_snippet(self):
        hits = client_text.draft_history("All good. The first draft of this report fell for it. "
                                         "Then more.")
        self.assertTrue(all(isinstance(h, str) and h for h in hits))
        self.assertTrue(any("first draft" in h for h in hits))


class Podcast(Workspace):
    def test_draft_history_in_the_transcript_fails_the_check(self):
        seg = list(SEG2) + [("LEO", "And the first draft of this report fell for it, which we "
                                    "corrected.")]
        write(self.script, script_text(segs=(SEG1, seg)))
        rc, _out, _err = run("check", self.script, "--source", self.report)
        txt = read(os.path.join(self.pod, "check.txt"))
        self.assertEqual(rc, 1, txt)
        self.assertIn("draft history", txt)
        self.assertIn("first draft", txt)

    def test_a_clean_script_still_passes(self):
        write(self.script, script_text())
        rc, _out, _err = run("check", self.script, "--source", self.report)
        self.assertEqual(rc, 0, read(os.path.join(self.pod, "check.txt")))

    def test_the_brief_says_so(self):
        with open(os.path.join(os.path.dirname(HERE), "references", "podcast.md"),
                  encoding="utf-8") as fh:
            self.assertIn("draft history", fh.read().lower())


class ReportLint(unittest.TestCase):
    def build(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            rep = os.path.join(tmp, "AI_Analysis_Report.md")
            with open(rep, "w", encoding="utf-8") as fh:
                fh.write(text)
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                                "--report", rep, "--out", os.path.join(tmp, "Analysis_Report.html"),
                                "--no-pdf"], capture_output=True, text=True, timeout=120)
            return r

    def test_draft_history_in_the_report_is_a_warning_not_a_failure(self):
        r = self.build("# Report\n\n## Expert Review Notes\n\nWe checked the first draft; two "
                       "protein names were wrong and have been corrected.\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("draft history", r.stderr)
        self.assertIn("first draft", r.stderr)

    def test_a_report_without_it_says_nothing(self):
        r = self.build("# Report\n\n## Expert Review Notes\n\nThe review checked the protein "
                       "names against UniProt (Key Findings) and the normalisation (Audit & "
                       "caveats). P-values were corrected for multiple testing.\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn("draft history", r.stderr)


class Documented(unittest.TestCase):
    def read(self, *rel):
        with open(os.path.join(os.path.dirname(HERE), *rel), encoding="utf-8") as fh:
            return " ".join(fh.read().split())

    def test_expert_review_notes_state_what_was_checked_never_draft_history(self):
        for rel in (("SKILL.md",), ("references", "anomaly-checks.md")):
            with self.subTest(rel=rel):
                t = self.read(*rel)
                self.assertIn("present tense", t)
                self.assertIn("logs/decisions.md", t)
                self.assertIn("draft history", t.lower())


if __name__ == "__main__":
    unittest.main(verbosity=2)
