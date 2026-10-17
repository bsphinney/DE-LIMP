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

# Draft history that only narrates a draft or a fix: the STRICT tier (podcast check FAILS).
STRICT = (
    "The first draft of this report fell for it.",
    "We checked the first draft against the tables.",
    "Two protein names were wrong and have been corrected.",
    "An earlier normalisation reason has since been corrected.",
    "In an earlier version, the volcano was mislabelled.",
    "A previous version of this report listed 12 proteins.",
    "Honest correction: the bait is not enriched.",
    "We got that wrong at first.",
    "This version of the report fixes the legend.",
    "The report now says the controls were IgG.",
    "Our earlier draft overstated it.",
    "The initial draft called it significant.",
    # review of 2.13.2 (M5): the false negatives
    "The expert review caught two mislabelled proteins, which we then fixed.",
    "Originally the report said Ryr2 was down.",
    "The report originally said Ryr2 was down.",
    "Two protein names were mislabelled; they are fixed now.",
    "We've updated the wording after the review.",
    "The review flagged that the normalisation reason was incorrect, and it has been changed.",
    "Our first pass of this analysis missed it.",
)
# Draft history only in context: the LOOSE tier (podcast check WARNS).
LOOSE = (
    "Two protein names were wrong; both are now corrected.",
    "We corrected the protein names after the review.",
    "The draft said the opposite.",
    "We fixed the mistake in the legend.",
    "That trap is easy to fall into, and we fell into it.",
    "The error was corrected before delivery.",
    "We had that backwards before.",
    "An earlier pass named the wrong protein.",
    "This was initially reported as significant.",
)
# Left to the meaning check's (f): no phrase tells these from science.
MEANING_ONLY = ("The categories were corrected.",)
# Ordinary science: neither tier.
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
    "The label on the tube was wrong, so the sheet lists S03 twice.",
    "Expression now increases with temperature.",
    "Cells were drafted into the experiment in two batches.",
    "The mass error was corrected by recalibration.",
    "Mass errors were corrected using the lock mass.",
    "Use this table for your first draft of the manuscript.",
    "The draft manuscript is attached.",
    "Values were corrected (Benjamini-Hochberg).",
    # review of 2.13.2 (M4): the false positives
    "Its abundance fell for this group of mice.",
    "Most proteins fell into this category.",
    "These are corrected p-values, so a value below 0.05 means something.",
    "The p-values are corrected.",
    "Retention times were corrected, then the runs were aligned.",
    "Isotope impurities were corrected.",
    "We corrected these for multiple testing.",
    "We corrected the systematic mass errors by recalibrating each run.",
    "In earlier versions, DIA-NN used a different scoring.",
    "This draft proteome map is from 2014.",
    "Batch effects were corrected, and the groups then separated.",
    "The intensity values are corrected, not raw.",
    # and more
    "Retention time drift was corrected before alignment.",
    "We updated the UniProt database to the 2026_03 release.",
    "The search was repeated with an updated FASTA.",
    "The first pass of DIA-NN identified fewer precursors than the second.",
    "An earlier study reported this interaction.",
    "The protein was initially described in yeast.",
    "Values were fixed at 15 ppm for both MS levels.",
    "The false discovery rate was 1%, and the error rate of the test is controlled.",
    "Ratios were corrected for label incorporation.",
    "We changed the normalisation to raw quantities because the loading differed.",
)


class Phrases(unittest.TestCase):
    def test_the_strict_tier_catches_what_only_narrates_a_draft(self):
        for text in STRICT:
            with self.subTest(text=text):
                self.assertTrue(client_text.draft_history(text), text)

    def test_the_loose_tier_catches_the_rest_by_context(self):
        for text in LOOSE:
            with self.subTest(text=text):
                self.assertEqual(client_text.draft_history(text), [], text)
                self.assertTrue(client_text.draft_history_loose(text), text)

    def test_ordinary_science_trips_neither_tier(self):
        for text in SCIENCE + MEANING_ONLY:
            with self.subTest(text=text):
                self.assertEqual(client_text.draft_history(text), [], text)
                self.assertEqual(client_text.draft_history_loose(text), [], text)

    def test_each_hit_is_a_readable_snippet(self):
        hits = client_text.draft_history("All good. The first draft of this report fell for it. "
                                         "Then more.")
        self.assertTrue(all(isinstance(h, str) and h for h in hits))
        self.assertTrue(any("first draft" in h for h in hits))

    def test_the_meaning_check_asks_for_draft_history_too(self):
        import make_podcast
        self.assertIn("f", [k for k, _ in make_podcast.MEANING_CHECKS])
        self.assertIn("draft", dict(make_podcast.MEANING_CHECKS)["f"])


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

    def test_a_loose_phrase_is_a_warning_not_a_failure(self):
        seg = list(SEG2) + [("LEO", "We corrected the protein names after the review.")]
        write(self.script, script_text(segs=(SEG1, seg)))
        rc, _out, _err = run("check", self.script, "--source", self.report)
        txt = read(os.path.join(self.pod, "check.txt"))
        self.assertEqual(rc, 0, txt)
        self.assertIn("possible draft history", txt)

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
