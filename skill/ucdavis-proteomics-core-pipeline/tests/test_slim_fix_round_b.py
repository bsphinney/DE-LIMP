#!/usr/bin/env python3
"""The 2.15 slim port's router and reference fixes (review REVIEW_215_port.md, Ruling 93, fixer B).

The slim-down's router shortened the routes and moved text into references; an independent
review found the places where that lost a rule an agent acts on:

* H2: the one-search no-DE route listed "0a, 0, 2, 3 (organism), 3b, 4-7b, 9-12", so a "run Sage
  on these hair files" request skipped step 3's other questions (sample type, keratin, the
  identifications-only choice), step 1's hive_path.sh, 7c and 8e. It is "every step, with no DE"
  now, and no-de-route.md's route names the steps before the search.
* H3: the podcast's "a free-tier key is for published or synthetic data only" lived only in a
  section bannered "Superseded"; step 9e said only "suggest a paid-tier key".
* M1: the resume route skipped 0a (no core_member, no mode). M2: "(Core: 8e)" was dropped from
  step 3 and the step-8 header, and ht-submissions.md never said a plate with no DE runs 8e.
  M3: step 7's "commercial -> AlphaDIA" pointed nowhere near AlphaDIA's timsTOF limit.
* LOW 4 (Ruling 93): the rechecks before steps 7 and 12 are `skill_version.sh --check-update`
  alone, as access.md and notes.md say; `--session-start` there could run `next` (--put-skill)
  unasked in the middle of an analysis.
* The smaller ones: the description's triggers, golden rules 8 and 11, 0b in the route, the
  CPU-family pointer at 11b, install.md "When something is missing", 3b's --remember, the
  --slurm watch form's --hive, stale "SKILL.md step 12b" pointers, and the rest below.

Hermetic: SKILL.md and the references are read; nothing runs.
"""
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
REFS = os.path.join(SKILL, "references")


def read(*parts):
    with open(os.path.join(SKILL, *parts), encoding="utf-8") as fh:
        return fh.read()


def flat(text):
    return " ".join(text.split())


def between(text, start, end):
    a = text.index(start)
    return text[a:text.index(end, a + len(start))]


def router():
    return read("SKILL.md")


def step(head, nxt):
    return flat(between(router(), head, nxt))


def section(rel, head):
    text = read("references", rel)
    a = text.index(head)
    b = text.find("\n## ", a + 3)
    return flat(text[a:b if b >= 0 else len(text)])


def route():
    return step("## Route the request", "## Flow")


def description():
    fm = router().split("\n---\n", 1)[0]
    body = fm.split("description: >\n", 1)[1]
    return " ".join(ln.strip() for ln in body.splitlines())


class TheRoutes(unittest.TestCase):
    def test_one_search_with_no_de_runs_every_step(self):
        """H2."""
        r = route()
        self.assertIn('One search, no DE ("run DIA-NN on these"): every step, with no DE '
                      '(no-de-route.md)', r)
        self.assertNotIn("3 (organism)", r)
        self.assertNotIn("4–7b, 9–12", r)


class Step7(unittest.TestCase):
    def test_the_commercial_clause_points_at_alphadias_limit(self):
        """M3: AlphaDIA does not work on timsTOF whole-proteome directDIA."""
        s7 = step("\n### 7. Run the search", "\n### 7b. ")
        self.assertIn('(commercial: search-engines.md "7. Licensing")', s7)
        self.assertNotIn("commercial → AlphaDIA", s7)
        lic = section("search-engines.md", "## 7. Licensing, AlphaDIA, chromatograms")
        self.assertIn("directDIA", lic)
        self.assertIn("timsTOF", lic)


class TheMidAnalysisRechecks(unittest.TestCase):
    """LOW 4 (Ruling 93): `--check-update` alone before steps 7 and 12, in the router and in both
    references, never an unasked --put-skill in the middle of an analysis."""

    def test_the_router_says_check_update_alone_at_all_three_places(self):
        s0 = step("\n### 0. ", "\n### 0b. ")
        again = s0[s0.index("**Check again**"):]
        self.assertIn("`skill_version.sh --check-update` alone", again)
        self.assertNotIn("skill_version.sh --session-start", again)
        head = s0[:s0.index("bash scripts/skill_version.sh --session-start")]
        self.assertIn("before steps 7 and 12, `--check-update` alone", head)
        s7 = step("\n### 7. Run the search", "\n### 7b. ")
        before = s7[s7.index("Before submitting:"):s7.index("notifications.md")]
        self.assertIn("`--check-update` alone (step 0)", before)
        self.assertNotIn("--session-start", before)
        s12 = step("\n### 12. ", "\n### 12b. ")
        self.assertIn("`--check-update` alone (step 0)", s12)
        self.assertNotIn("--session-start", s12)

    def test_the_references_agree(self):
        access = flat(read("references", "access.md"))
        self.assertIn('`--check-update` alone is what step 0\'s "Check again" runs before steps 7 '
                      'and 12', access)
        notes = section("notes.md", "## Step 0: the notes check")
        self.assertIn("**The version too**, at the same two points: "
                      "`bash scripts/skill_version.sh --check-update`", notes)


class ThePodcastsFreeTierRule(unittest.TestCase):
    """H3."""

    def test_9e_states_it(self):
        s9e = section("podcast.md", "## 9e. The audio discussion (step 9e)")
        self.assertIn("a free-tier key is for published or synthetic data only; for unpublished "
                      "data a paid-tier key or `--tts say`, and say what the tier means when you "
                      "ask", s9e)
        self.assertIn('"Privacy and consent (read before any cloud voice)" above', s9e)

    def test_the_privacy_section_is_not_bannered_and_keeps_its_rules(self):
        priv = section("podcast.md", "## Privacy and consent")
        self.assertNotIn("Superseded", priv)
        self.assertIn("`--core-run` and `--cloud-ok` cannot be given together", priv)
        self.assertIn("A free-tier key is for published or synthetic data only", priv)
        self.assertIn("`core_submission.key_mode_problem`", priv)       # the key-file rules


class SupersededBlocksMoved(unittest.TestCase):
    """M4's named instances (the banner rule itself: test_skill_structure.BannersAreTrue)."""

    def test_methods_9d_carries_the_which_search_block(self):
        s9d = section("methods.md", "## 9d. The Methods section")
        for want in ("**Which search, and what its record says (2.14.2):**",
                     "`session.search_choice`", "optimised automatically per run",
                     "`results.json`", "`fragger.params`", "An identifications-only search",
                     "between-run normalisation", "`--search-prov <file>` names one by hand"):
            self.assertIn(want, s9d, want)
        self.assertEqual(flat(read("references", "methods.md")).count("Which search, and what its "
                                                                     "record says"), 1)


class StalePointers(unittest.TestCase):
    """LOW 11: "SKILL.md step 12b" text lives in core-submissions.md "12b." (staff only)."""

    def test_no_reference_points_at_router_text_that_moved(self):
        for f in sorted(os.listdir(REFS)):
            if f.endswith(".md"):
                t = flat(read("references", f))
                with self.subTest(ref=f):
                    self.assertNotIn("SKILL.md step 12b", t)
                    self.assertNotIn("SKILL.md 12b", t)
        self.assertNotIn("→ detail: `references/methods.md`",
                         section("methods.md", "## 9d. The Methods section"))


class SmallerRouterFixes(unittest.TestCase):
    def test_the_description_keeps_its_triggers(self):
        """LOW 7: a hand-written job would escape golden rule 6 if the skill does not trigger."""
        d = description()
        self.assertLessEqual(len(d), 1024)
        for want in ("submit a search to the cluster / SLURM / HPC", "extract XICs",
                     "diann-linux / sage command line"):
            self.assertIn(want, d)


    def test_3b_registers_the_conversation_on_hive(self):
        """LOW 14."""
        self.assertIn('hive_remote: `save_transcript.py --hive "$S" --remember` right after '
                      '`init`', step("\n### 3b. ", "\n### 4. "))

    def test_the_slurm_watch_form_has_hive(self):
        """LOW 16."""
        self.assertIn("`--slurm <id> --log <log> --out <out> --poll <n> --hive`",
                      step("\n### 7b. ", "\n### 7c. "))

    def test_the_router_stays_within_its_lines(self):
        self.assertLessEqual(router().count("\n"), 499)


class SmallerReferenceFixes(unittest.TestCase):
    def test_an_override_names_the_role(self):
        """LOW 3: the S4 rule, "Core staff" only for a login on the staff list."""
        s4 = section("parameters.md", "## 4. Resolving the defaults")
        self.assertNotIn('"set by Core staff"', s4)
        self.assertIn('"set by <Core staff|the user>" (the role, `staff.role_for`)', s4)

    def test_6b_points_at_the_one_rule(self):
        """LOW 2: the identifications-only criteria are said once, in no-de-route.md step 1."""
        s6b = section("parameters.md", "## 6b. Estimating the search parameters")
        self.assertIn("**unless amounts are wanted** (the one rule: "
                      "`references/no-de-route.md`, step 1)", s6b)
        self.assertNotIn("abundance between runs, an intensity table, a DE later", s6b)


    def test_step_10_says_what_session_fills_and_the_one_replay_rule(self):
        """Ruling 96 (fixer A's H1): --session takes finalize's flags, and the bundle replays the
        params file and the FASTA as the search ran them."""
        s10 = section("reproducibility.md", "## 10. The bundle")
        first = s10[s10.index("Inside an analysis session"):s10.index("An explicit argument wins")]
        for want in ("the value finalize gives it (`make_deposit.bundle_flags`, one definition",
                     "an identifications-only search's `identifications.tsv`",
                     "the params file and the FASTA as the search ran them (the search record's "
                     "first, with a note when the session's own file differs)",
                     "the peak check's record"):
            self.assertIn(want, first, want)
        self.assertNotIn("`input/acquisition.json`", first)
        self.assertNotIn("the files the search ran on", first)

    def test_audit_uses_the_sessions_acquisition_json(self):
        """LOW 16."""
        s8c = section("audit.md", "## 8c.")
        self.assertIn("--acquisition-json <session>/input/acquisition.json", s8c)
        self.assertNotIn("/tmp/acq.json", s8c)


if __name__ == "__main__":
    unittest.main(verbosity=2)
