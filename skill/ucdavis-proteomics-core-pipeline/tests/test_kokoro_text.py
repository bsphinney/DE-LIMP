#!/usr/bin/env python3
"""kokoro_text.py: each text rule for the local podcast voice does its job, and ordinary
sentences are left alone.

The rules come from Core episodes rendered with Kokoro on HIVE (2026-10), each misreading found
by two Whisper transcriptions that agreed: spelled-out letters read as the word "a", ISO dates,
versions, a PROT number's leading zero, signs and number ranges. The examples here are made up.
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

from kokoro_text import apply, fix_table  # noqa: E402


def fix(text):
    return apply(text)[0]


class Letters(unittest.TestCase):
    def test_spelled_symbols_become_capitals(self):
        self.assertEqual(fix("The channel subunit k c n b two is at plus 4.10."),
                         "The channel subunit K C N B two is at plus 4.10.")
        self.assertEqual(fix("delta j p h against the wild type"), "delta J P H against the wild type")
        self.assertEqual(fix("like c s g a, f i m a or w z a."), "like C S G A, F I M A or W Z A.")
        self.assertEqual(fix("stress proteins, h s p a and Slp"), "stress proteins, H S P A and Slp")
        self.assertEqual(fix("we were looking for a r y r effect."),
                         "we were looking for a R Y R effect.")

    def test_article_before_a_capital_is_a_word(self):
        for s in ("So a b-cell is not a T cell.", "a B cell and a T cell", "is it a G or a C?"):
            self.assertEqual(fix(s), s)

    def test_already_capital_is_untouched_and_unreported(self):
        text, changes = apply("H S P A and S L P")
        self.assertEqual(text, "H S P A and S L P")
        self.assertEqual(changes, [])

    def test_table_form(self):
        self.assertEqual(fix_table([("CsgA", "c s g a"), ("OmpT", "omp t"), ("SILAC", "sigh lack")]),
                         [("CsgA", "C S G A"), ("OmpT", "omp t"), ("SILAC", "sigh lack")])

    def test_symbol_starting_with_a(self):
        # in finished text a leading "a" reads as the article; the table form spells it out
        self.assertEqual(fix("the kinase a k t"), "the kinase a K T")
        self.assertEqual(fix_table([("AKT1", "a k t one")]), [("AKT1", "A K T one")])


class Dates(unittest.TestCase):
    def test_iso_date_in_words(self):
        self.assertEqual(fix("The same data were analysed on 2026-10-14 with a hand-made analysis."),
                         "The same data were analysed on October 14, 2026 with a hand-made analysis.")
        self.assertEqual(fix("(2026-01-05)"), "(January 5, 2026)")

    def test_not_a_date(self):
        for s in ("run 2026-13-40 failed", "folder 2026-10-15_PROT_0807 is", "a 1-2-3 count"):
            self.assertEqual(fix(s), s.replace("PROT_0807", "P R O T zero eight zero seven"))


class Versions(unittest.TestCase):
    def test_versions_spelled(self):
        self.assertEqual(fix("built-in SILAC support, version 2.13.0."),
                         "built-in SILAC support, version 2 point 13 point 0.")
        self.assertEqual(fix("DIA-NN 2.7.0 was used"), "DIA-NN 2 point 7 point 0 was used")
        self.assertEqual(fix("skill version 2.13 and v2.9"), "skill version 2 point 13 and version 2 point 9")
        self.assertEqual(fix("v1.2.3"), "1 point 2 point 3")

    def test_decimals_stay(self):
        for s in ("Ryr2 is at plus 5.97, which is about 60-fold.", "C1qa at plus 2.30.",
                  "p below 0.05", "a 1.5-fold change", "median minus 0.86, and 4.10."):
            self.assertEqual(fix(s), s)


class ProtNumbers(unittest.TestCase):
    def test_prot_keeps_its_zero(self):
        self.assertEqual(fix("an AI-generated discussion of the PROT_0807 report."),
                         "an AI-generated discussion of the P R O T zero eight zero seven report.")
        self.assertEqual(fix("PROT 0807 and prot-0012"),
                         "P R O T zero eight zero seven and P R O T zero zero one two")
        self.assertEqual(fix("PROT_0807's results"), "P R O T zero eight zero seven's results")

    def test_words_starting_with_prot_stay(self):
        for s in ("protein groups", "PROTEOMICS core", "a protease", "Protocol 2"):
            self.assertEqual(fix(s), s)


class PlusAndRanges(unittest.TestCase):
    def test_plus_sign(self):
        self.assertEqual(fix("proteasome +0.61, translation factors +0.32"),
                         "proteasome plus 0.61, translation factors plus 0.32")
        self.assertEqual(fix("it's +5.89. But"), "it's plus 5.89. But")
        self.assertEqual(fix("(+7.14)"), "(plus 7.14)")

    def test_plus_inside_expressions_stays(self):
        for s in ("1+1 is 2", "log2(x+1)", "C++ code", "Na+ and K+", "a +b"):
            self.assertEqual(fix(s), s)

    def test_ranges(self):
        self.assertEqual(fix("detected 51,204–58,977 precursors"), "detected 51,204 to 58,977 precursors")
        self.assertEqual(fix("91–97% of the 8,315"), "91 to 97% of the 8,315")
        self.assertEqual(fix("between 18 — 24 degrees"), "between 18 to 24 degrees")

    def test_hyphen_ranges(self):
        self.assertEqual(fix("each mutant carried 13-17% of the signal"),
                         "each mutant carried 13 to 17% of the signal")
        self.assertEqual(fix("between 1,200-1,500 proteins"), "between 1,200 to 1,500 proteins")
        self.assertEqual(fix("a 0.5-1.5 window"), "a 0.5 to 1.5 window")
        self.assertEqual(fix("between 2020-2021"), "between 2020 to 2021")     # a year range

    def test_hyphens_that_are_not_ranges_stay(self):
        for s in ("X1234-C1", "a–b", "pre–post", "Lys8-labelled", "COVID-19 samples", "1-2-3 order",
                  "in 2026-10 we ran", "the T-cell", "anti-5"):
            self.assertEqual(fix(s), s)

    def test_minus_sign(self):
        self.assertEqual(fix("Jph3 is at −1.82, Kcnb2 at -3.42"),
                         "Jph3 is at minus 1.82, Kcnb2 at minus 3.42")
        self.assertEqual(fix("median (−0.86)"), "median (minus 0.86)")

    def test_minus_inside_words_and_expressions_stays(self):
        for s in ("x-1", "COVID-19", "anti-5", "log2(x)-1"):
            self.assertEqual(fix(s), s)


class OrdinarySentences(unittest.TestCase):
    SENTENCES = [
        "Imagine you could listen to three bacterial strains at once, on one tape.",
        "This is a test, and I think a sample is a sample.",
        "I'm Maya, the biologist of the pair.",
        "Plan A or plan B? Vitamin A and B12.",
        "1,412 protein groups had a ratio in at least three of four replicates.",
        "615 differ significantly: 300 relatively enriched and 315 relatively depleted.",
        "In this log2 fold change, or logFC, plus 1 means a doubling and minus 1 a halving.",
        "The three supernatants were combined 1:1:1 by volume.",
        "Lys4 and Arg6, a few neutrons heavier; Lys8 and Arg10.",
        "An x-ray, a T-cell and an e-mail.",
        "So a b-cell is not a T cell.",
        "It ran on 2026 hardware in 36 isolation windows.",
        "A 1% false-discovery rate, the FDR, means 1% are expected wrong.",
    ]

    def test_unchanged(self):
        for s in self.SENTENCES:
            text, changes = apply(s)
            self.assertEqual(text, s, f"changed: {changes}")
            self.assertEqual(changes, [])


class Reporting(unittest.TestCase):
    def test_changes_listed(self):
        _, changes = apply("PROT_0807 on 2026-10-14, version 2.13.0, 91–97%, +0.61, c s g a")
        self.assertEqual([c[0] for c in changes],
                         ["prot", "date", "version", "range", "plus", "letters"])

    def test_signs_and_hyphen_ranges_listed(self):
        _, changes = apply("68-72% and −0.86")
        self.assertEqual([c[0] for c in changes], ["range", "minus"])


if __name__ == "__main__":
    unittest.main()
