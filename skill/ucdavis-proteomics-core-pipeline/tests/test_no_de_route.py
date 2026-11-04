#!/usr/bin/env python3
"""
The no-DE route is signposted where the agent decides it, and again where it would stop (T14:
I03; Rulings 37 and 38).

Three sessions with no differential expression left the skill's path after the search: one
wrote its own report and delivered by hand; one (an independent re-search, no traditional DE)
ended with a hand-made reproducibility zip and no AGENTS.md, README, bundle or finalize; one
filed everything numbered after step 8 as "after the DE / later" and delivered raw + search
results with its own copy script. Nothing in steps 3, 4 or 8 said that a session with no DE
still runs the report, the Methods, finalize and the delivery.

Now: a short block at the end of step 3 (the design questions), one line in step 4, the
branch at the TOP of step 8, and the route itself in references/no-de-route.md. The analysis
brief of a session with no DE is not titled "differential expression".

Hermetic: SKILL.md and the reference are read; analysis_prompt.py runs on a temp dir.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
REF = os.path.join(SKILL, "references", "no-de-route.md")


def skill_md():
    with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
        return fh.read()


def between(text, start, end):
    a = text.index(start)
    return text[a:text.index(end, a)]


def flat(text):
    return " ".join(text.split())


class TheSignposts(unittest.TestCase):
    def test_skill_md_names_the_no_de_route_before_step_4(self):
        s3 = flat(between(skill_md(), "\n### 3. ", "\n### 4. "))
        self.assertIn("deliver --mode no-de", s3)
        self.assertIn("--ids-only", s3)
        self.assertIn("references/no-de-route.md", s3)
        self.assertIn("NOT done", s3)
        self.assertIn("finalize", s3)

    def test_step_3_block_comes_before_3b(self):
        text = skill_md()
        block = text.index("**No differential expression?**")
        self.assertLess(text.index("\n### 3. "), block)
        self.assertLess(block, text.index("\n### 3b. "))

    def test_step_4_confirms_it_with_the_defaults(self):
        s4 = flat(between(skill_md(), "\n### 4. ", "\n### 5. "))
        self.assertIn("If step 3 chose the no-DE route, confirm it with the defaults", s4)

    def test_step_8_opens_with_the_no_de_branch(self):
        """Ruling 38: the branch is at the TOP of step 8, before anything a DE needs."""
        s8 = between(skill_md(), "\n### 8. Differential expression", "\n### 8a. ")
        top = flat("\n".join(s8.splitlines()[1:5]))
        self.assertIn("No DE wanted?", top)
        self.assertIn("NOT done", top)
        self.assertIn("deliver --mode no-de", top)
        self.assertIn("references/no-de-route.md", top)
        self.assertLess(s8.index("No DE wanted?"), s8.index("SILAC"))

    def test_the_route_matches_the_finish_steps(self):
        """The reference: the steps in order, finalize making the documents for every analysis
        (golden rule 9 and step 12 as they are), the interim delivery, Ruling 37's two lines."""
        with open(REF, encoding="utf-8") as fh:
            ref = flat(fh.read())
        for want in ("estimate_params.py --ids-only", "make_analysis_html.py", "9d",
                     "session.py finalize", "deliver --mode no-de", "deliver --mode interim",
                     "qc_bracket.py", "--code-dir", "IN ADDITION to the skill's bundle",
                     "session.py status", "AGENTS.md", "re-analysis", "benchmark"):
            self.assertIn(want, ref, want)
        self.assertNotIn("provenance.py without --de-dir", ref,
                         "since L5 finalize builds the bundle from the session's records")
        order = [ref.index(x) for x in ("(step 7)", "(step 9)", "(step 9d)", "(step 12)",
                                        "(step 12b)")]
        self.assertEqual(order, sorted(order))


class TheReviewsGaps(unittest.TestCase):
    """Review F6, F7, F11: finalize only adds to an existing OUTPUT_FILES.md, so a route that
    skips step 11 leaves none (and 11a is the step that answers a benchmark of two searches);
    the "in addition to the bundle" sentence must be where the agent reads, in SKILL.md; the
    step-3 chain named no podcast (9e) while step 8's did."""

    def test_both_blocks_and_the_reference_name_step_11_and_11a(self):
        text = skill_md()
        s3 = flat(text[text.index("**No differential expression?**"):text.index("\n### 3b. ")])
        s8 = flat("\n".join(between(text, "\n### 8. Differential expression",
                                     "\n### 8a. ").splitlines()[1:6]))
        with open(REF, encoding="utf-8") as fh:
            ref = flat(fh.read())
        for where, t in (("step 3", s3), ("step 8", s8), ("reference", ref)):
            self.assertIn("OUTPUT_FILES.md", t, where)
            self.assertIn("11a", t, where)
        order = [ref.index(x) for x in ("(step 9)", "(step 11)", "(step 12)")]
        self.assertEqual(order, sorted(order))

    def test_step_3_says_in_addition_to_the_bundle_and_names_9e(self):
        text = skill_md()
        s3 = flat(text[text.index("**No differential expression?**"):text.index("\n### 3b. ")])
        self.assertIn("in addition to the bundle, never instead", s3)
        self.assertIn("9e (Core)", s3)
        self.assertLess(s3.index("9d"), s3.index("9e (Core)"))
        self.assertLess(s3.index("9e (Core)"), s3.index("12 `finalize`"))


class TheBrief(unittest.TestCase):
    def brief(self, tables):
        with tempfile.TemporaryDirectory() as d:
            de = os.path.join(d, "tables")
            os.makedirs(de)
            for name, text in tables.items():
                with open(os.path.join(de, name), "w") as fh:
                    fh.write(text)
            out = os.path.join(d, "ANALYSIS_PROMPT.md")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                                "--de-dir", de, "--out", out], capture_output=True, text=True,
                               timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out, encoding="utf-8") as fh:
                return fh.read()

    def test_brief_without_de_is_not_titled_differential_expression(self):
        b = self.brief({})
        first = next(ln for ln in b.splitlines() if ln.startswith("# "))
        self.assertIn("identification and QC", first)
        self.assertNotIn("differential expression results", first.replace(
            "(no differential expression)", ""))
        self.assertIn("(no differential expression)", first)
        self.assertNotIn("analysis of the differential expression results", flat(b))

    def test_a_de_brief_keeps_its_title(self):
        """(pin)"""
        b = self.brief({"DE_B_vs_A.csv": "Protein.Group,logFC,adj.P.Val\nP1,1.0,0.01\n",
                        "de_provenance.json": json.dumps({"contrasts": ["B-A"]})})
        first = next(ln for ln in b.splitlines() if ln.startswith("# "))
        self.assertEqual(first, "# Analysis brief — proteomics differential expression")
        self.assertTrue(re.search(r"analysis of the differential expression\s+results", b))


if __name__ == "__main__":
    unittest.main()
