#!/usr/bin/env python3
"""
A podcast is checked for MEANING before it is rendered (make_podcast.py meaning; check's
robustness and superlative lines).

A staff episode (PROT_0000, skill 2.11.2, 2026-10-16) passed `check` (0 problems, 0 warnings),
yet a separate read of every line against the report found a hit from one comparison told as
part of another comparison's story and called "nearly absent" (detected in 2 of 3 cultures),
proteins put right after "It survives everything" whose robustness table says robust = no or
"not testable", "not a main driver" spoken as "aren't driving", and a tool presented as one the
client already uses. check
matches tokens and numbers only. Now a reader who did not write the script records a verdict for
the script's exact sha256 (meaning_check.md), and render, link, share and job refuse without its
PASS; what can be caught mechanically is caught by check.

Stdlib only; the TTS is a fake.
"""
import json
import os
import sys
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
sys.path.insert(0, HERE)

from test_make_podcast import (REPORT, SEG1, SEG2, FakeTTS, Workspace, mp,  # noqa: E402
                               read, run, script_text, write)

ROBUST_CSV = ('"Protein.Group","Genes","logFC","adj.P.Val","robust"\n'
              '"P1","Ryr2",2.1,0.001,"no"\n'
              '"P2","Jph3",-1.3,0.002,"yes"\n'
              '"P3","C1qa",0.9,0.01,"not testable in S1"\n')
ROBUST_CSV_2 = ('"Protein.Group","Genes","robust"\n'
                '"P3","C1qa","yes"\n')


class Base(Workspace):
    def setUp(self):
        super().setUp()
        mp.BACKENDS["fake"] = FakeTTS
        FakeTTS.calls = []
        self.addCleanup(mp.BACKENDS.pop, "fake", None)
        enc = mock.patch.object(mp, "encode_aac", return_value=(None, "test"))
        enc.start()
        self.addCleanup(enc.stop)
        self.findings = os.path.join(self.pod, "findings.txt")
        write(self.findings, "- line 31: the hit is A vs B, not C vs B -> fixed (resolved)\n")

    def checked(self, text=None):
        write(self.script, text if text is not None else script_text())
        rc, out, err = run("check", self.script, "--source", self.report)
        return rc, read(os.path.join(self.pod, "check.txt"))

    def meaning(self, *extra):
        return run("meaning", self.script, *extra)

    def record(self, verdict="PASS", reviewer="a reviewer agent that did not write the script"):
        return self.meaning("--verdict", verdict, "--reviewer", reviewer, "--findings",
                            self.findings)

    def render(self, *extra):
        return run("render", self.script, "--tts", "fake", *extra)


class Gate(Base):
    def test_render_needs_a_pass_for_this_script(self):
        self.assertEqual(self.checked()[0], 0)
        rc, out, err = self.render()
        self.assertEqual(rc, 2)
        self.assertIn("no meaning_check.md beside podcast_script.md", err)
        self.assertIn("make_podcast.py meaning", err)
        self.assertEqual(FakeTTS.calls, [])
        # a FAIL verdict is recorded (exit 1) and still refused
        self.assertEqual(self.record("FAIL")[0], 1)
        rc, out, err = self.render()
        self.assertEqual(rc, 2)
        self.assertIn("meaning_check.md does not say PASS", err)
        # PASS for this script: render goes on and records it
        self.assertEqual(self.record()[0], 0)
        mc = read(os.path.join(self.pod, "meaning_check.md"))
        self.assertIn(f"script_sha256: {mp.parse_script(self.script).sha256}", mc)
        self.assertIn("verdict: PASS", mc)
        for k, _ in mp.MEANING_CHECKS:
            self.assertIn(f"- ({k}) ", mc)
        self.assertIn("A vs B, not C vs B", mc)                                   # the findings
        rc, out, err = self.render()
        self.assertEqual(rc, 0, err)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertEqual(man["meaning_check"]["status"], "PASS")
        self.assertFalse(man["meaning_check"]["overridden"])
        # an edit after the meaning check: the check AND the meaning check are stale
        write(self.script, script_text() + "\n")
        self.assertEqual(self.checked(read(self.script))[0], 0)
        rc, out, err = self.render()
        self.assertEqual(rc, 2)
        self.assertIn("the script changed after meaning_check.md was recorded", err)

    def test_an_override_is_recorded_with_its_reason(self):
        self.assertEqual(self.checked()[0], 0)
        rc, out, err = self.render("--no-meaning-check", "  the director read it  ")
        self.assertEqual(rc, 0, err)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertTrue(man["meaning_check"]["overridden"])
        self.assertEqual(man["meaning_check"]["status"], "missing")
        self.assertEqual(man["meaning_check"]["override_reason"], "the director read it")

    def test_the_writer_cannot_be_the_reviewer(self):
        head = ["# Signal to Noise — x", "", "Show: Signal to Noise", "Title: x",
                "Hosts: Maya (cell biologist), Leo (statistician)", "Voices: Maya=Kore, Leo=Charon",
                "Author: Claude (the session's writer)"]
        self.assertEqual(self.checked(script_text(header=head))[0], 0)
        rc, out, err = self.record(reviewer="claude (the session's writer)")
        self.assertEqual(rc, 2)
        self.assertIn("is the script's author", err)
        self.assertFalse(os.path.exists(os.path.join(self.pod, "meaning_check.md")))

    def test_findings_and_a_reviewer_are_required(self):
        self.assertEqual(self.checked()[0], 0)
        self.assertEqual(self.meaning("--verdict", "PASS", "--findings", self.findings)[0], 2)
        write(self.findings, "  \n")
        self.assertEqual(self.record()[0], 2)

    def test_the_brief_lists_what_to_check_and_every_sentence(self):
        self.assertEqual(self.checked()[0], 0)
        rc, out, err = self.meaning("--brief", "--source", self.report)
        self.assertEqual(rc, 0, err)
        brief = read(os.path.join(self.pod, "meaning_brief.md"))
        self.assertEqual(brief.strip(), out.strip())                     # printed and written
        self.assertIn("You did NOT write this script", brief)
        for k, _ in mp.MEANING_CHECKS:
            self.assertIn(f"- ({k}) ", brief)
        self.assertIn("Leo: Kv2.1 has 215 proteins at both ages; 77 of them are Kv2.1-only.",
                      brief)
        self.assertIn("../AI_Analysis_Report.md", brief)


class LinkShareJob(Base):
    def rendered(self):
        self.assertEqual(self.checked()[0], 0)
        self.assertEqual(self.record()[0], 0)
        self.assertEqual(self.render()[0], 0)

    def test_link_refuses_once_the_meaning_check_no_longer_holds(self):
        self.rendered()
        self.assertEqual(run("link", self.out)[0], 0)
        os.remove(os.path.join(self.pod, "meaning_check.md"))
        rc, out, err = run("link", self.out)
        self.assertEqual(rc, 2)
        self.assertIn("no meaning_check.md", err)

    def test_an_episode_rendered_before_the_meaning_check_only_warns(self):
        self.rendered()
        man_path = os.path.join(self.pod, "podcast.json")
        man = json.loads(read(man_path))
        del man["meaning_check"]
        write(man_path, json.dumps(man))
        os.remove(os.path.join(self.pod, "meaning_check.md"))
        mp._WARNED.clear()
        rc, out, err = run("link", self.out)
        self.assertEqual(rc, 0, err)
        self.assertIn("rendered before the meaning check existed", err)

    def test_job_refuses_without_it(self):
        self.assertEqual(self.checked()[0], 0)
        rc, out, err = run("job", self.script, "--cloud-ok", "PI, today")
        self.assertEqual(rc, 2)
        self.assertIn("no meaning_check.md", err)


class Mechanical(Base):
    def tables(self, *extra):
        write(os.path.join(self.out, "tables", "robustness", "Robustness_A.B.csv"), ROBUST_CSV)
        for i, text in enumerate(extra):
            write(os.path.join(self.out, "tables", "robustness", f"Robustness_C{i}.D.csv"), text)

    def script_saying(self, line):
        seg = list(SEG2) + [("LEO", line)]
        return script_text(segs=(SEG1, seg))

    def test_a_robust_claim_next_to_a_protein_the_tables_call_not_robust_fails(self):
        self.tables()
        rc, txt = self.checked(self.script_saying("It survives everything: Ryr2 and Jph3 both."))
        self.assertEqual(rc, 1, txt)
        self.assertIn("Ryr2 is named next to 'survives', but the robustness tables say otherwise "
                      "(tables/robustness/Robustness_A.B.csv: robust = no)", txt)
        self.assertNotIn("Jph3 is named next to", txt)
        self.assertIn("robustness tables: tables/robustness/Robustness_A.B.csv", txt)
        # and meaning --verdict PASS is refused while it stands
        rc, out, err = self.record()
        self.assertEqual(rc, 2)
        self.assertIn("refusing PASS", err)

    def test_a_negated_or_distant_claim_is_not_one(self):
        self.tables()
        for line in ("Ryr2 is not robust, and that is said in the report.",
                     "Ryr2 does not survive the second analysis."):
            with self.subTest(line=line):
                rc, txt = self.checked(self.script_saying(line))
                self.assertNotIn("Ryr2 is named next to", txt)

    def test_robust_in_one_comparison_only_is_a_warning(self):
        self.tables(ROBUST_CSV_2)
        rc, txt = self.checked(self.script_saying("C1qa is robust, says the report."))
        self.assertIn("C1qa is named next to 'robust'", txt)
        self.assertIn("robust in tables/robustness/Robustness_C0.D.csv; not in "
                      "tables/robustness/Robustness_A.B.csv: robust = not testable in S1", txt)
        self.assertNotIn("FAIL\n- turn 9", txt)

    def test_no_tables_means_nothing_to_check(self):
        rc, txt = self.checked(self.script_saying("It survives everything: Ryr2 and Jph3."))
        self.assertNotIn("Ryr2 is named next to", txt)
        self.assertIn("robustness tables: none", txt)

    def test_a_superlative_the_report_does_not_make_is_a_warning(self):
        write(self.report, REPORT + "\nRyr2 is the strongest hit of all.\n")
        rc, txt = self.checked(self.script_saying(
            "Ryr2 is the strongest hit, and Jph3 is the only partner that matters."))
        self.assertEqual(rc, 0, txt)                                   # a warning, not a fail
        self.assertIn("'the only partner' -- a superlative the sources do not say", txt)
        self.assertNotIn("'the strongest hit'", txt)


class Documented(unittest.TestCase):
    def read(self, *rel):
        with open(os.path.join(SKILL, *rel), encoding="utf-8") as fh:
            return " ".join(fh.read().split())

    def test_skill_md_and_the_brief_describe_the_reviewer_pass(self):
        for rel in (("SKILL.md",), ("references", "podcast.md")):
            with self.subTest(rel=rel):
                t = self.read(*rel)
                self.assertIn("make_podcast.py meaning", t)
                self.assertIn("meaning_check.md", t)
                self.assertIn("did not write", t)


if __name__ == "__main__":
    unittest.main(verbosity=2)
