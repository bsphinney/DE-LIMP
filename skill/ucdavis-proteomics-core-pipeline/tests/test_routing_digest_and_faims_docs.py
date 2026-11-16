#!/usr/bin/env python3
"""SKILL.md carries the digest answer through the steps, and routes FragPipe DDA to its own docs
(2.15 C2 fix round, docs H1/H2/MED/LOW).

The flags named here are checked against the scripts' argparse, so a rename breaks this test
instead of leaving the router describing a flag that no longer exists.
"""
import os
import re
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)


def read(*parts):
    with open(os.path.join(SKILL, *parts), encoding="utf-8") as fh:
        return fh.read()


def step(md, number):
    """The text of one '### <number>.' step of SKILL.md."""
    m = re.search(r"^### %s\.? .*?(?=^### |^## )" % re.escape(number), md, re.M | re.S)
    assert m, number
    return m.group(0)


def helptext(script):
    out = subprocess.run([sys.executable, "-I", os.path.join(SKILL, "scripts", script), "--help"],
                         capture_output=True, text=True, timeout=120)
    return out.stdout + out.stderr


class DigestAnswerIsAskedAndCarried(unittest.TestCase):
    def setUp(self):
        self.md = read("SKILL.md")

    def test_step3_asks_the_enzyme(self):
        s3 = step(self.md, "3")
        self.assertIn("Digestion enzyme", s3)
        self.assertIn("never infer", s3)

    def test_step3_asks_it_in_study_questions(self):
        sq = read("references", "study-questions.md")
        self.assertRegex(sq, r"(?m)^## Digestion enzyme and specificity")
        self.assertIn("--semi-specific", sq)

    def test_steps_4_6_6b_6c_carry_the_flags(self):
        self.assertIn("--semi-specific", step(self.md, "4"))
        self.assertIn("--enzyme", step(self.md, "6"))
        self.assertIn("--enzyme", step(self.md, "6b"))
        self.assertIn("--semi-specific", step(self.md, "6c"))
        self.assertIn("--semi-specific", read("references", "parameters.md").split("## 6b.")[1])

    def test_the_flags_exist_where_the_router_sends_them(self):
        self.assertIn("--semi-specific", helptext("resolve_defaults.py"))
        self.assertIn("--enzyme", read("scripts", "fetch_fasta.py"))
        for script in ("estimate_params.py", "make_presets.py"):
            h = helptext(script)
            self.assertIn("--enzyme", h, script)
            self.assertIn("--semi-specific", h, script)
        self.assertNotIn("--enzyme", helptext("run_search.py").replace("--enzyme-", ""))


class FragPipeDdaIsRouted(unittest.TestCase):
    def test_step7_cites_the_dda_sections_and_keeps_the_dia_one(self):
        s7 = step(read("SKILL.md"), "7")
        self.assertIn('"FragPipe (the timsTOF DDA-PASEF default"', s7)
        self.assertIn('parameters.md "FragPipe DDA: digest', s7)
        self.assertIn('fragpipe-diatracer.md "7."', s7)
        self.assertIn("--faims-split", s7)
        self.assertIn("--semi-specific", s7)

    def test_a_faims_split_has_a_written_finish(self):
        se = read("references", "search-engines.md")
        self.assertRegex(se, r"(?m)^## Finishing a FAIMS-split session")
        part = se.split("## Finishing a FAIMS-split session")[1].split("\n## ")[0]
        for word in ("choose-search", "faims_split.json", "--mode interim"):
            self.assertIn(word, part)
        self.assertIn("Finishing a FAIMS-split session", step(read("SKILL.md"), "12"))


class SmallStaleStatements(unittest.TestCase):
    def test_no_pointer_to_a_script_that_does_not_exist(self):
        self.assertFalse(os.path.exists(os.path.join(SKILL, "scripts", "sage_protein_groups.py")))
        self.assertNotIn("sage_protein_groups.py", read("references", "search-engines.md"))

    def test_session_comment_does_not_cite_step_12b_for_choosing(self):
        src = read("scripts", "session.py")
        self.assertNotIn("the refusals and SKILL.md step 12b say it", src)

    def test_transcript_flag_names_its_script(self):
        md = read("SKILL.md")
        self.assertNotIn("always use `--transcript-from-now`", md)
        self.assertIn("save_transcript.py --from-now", md)
        self.assertIn("--from-now", read("scripts", "save_transcript.py"))
        self.assertIn("--transcript-from-now", read("scripts", "session.py"))

    def test_step2_reference_says_faims_cvs(self):
        se = read("references", "search-engines.md")
        part = se.split("## 2. Detect acquisition + instrument")[1].split("\n## ")[0]
        self.assertIn("faims_cvs", part)


if __name__ == "__main__":
    unittest.main()
