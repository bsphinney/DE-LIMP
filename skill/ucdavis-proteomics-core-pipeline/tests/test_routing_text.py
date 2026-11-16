#!/usr/bin/env python3
"""The text says the routing the code does (2.15, F5; Rulings 55, 56, 67, 99).

engine_choice.ROUTES is the one rule: DIA -> DIA-NN; timsTOF DDA-PASEF -> FragPipe, with Sage as
the automatic fallback where FragPipe cannot run (its reason staff-only); Orbitrap and other DDA
-> Sage, FragPipe on request. 2.14.2's text said "DDA -> Sage" flatly, and the plan's first brief
said "DDA -> FragPipe where licensed"; both are wrong now. So:
  * no "DDA -> Sage" (any arrow) anywhere in the router, the references or the scripts without
    timsTOF or FragPipe on the same line to say which DDA it is;
  * no client-type text on the route (Ruling 56: no commercial/for-profit gate, no question) in
    the router's route paragraph, search-engines.md "4c." or the client sentence;
  * the router stays under 500 lines (499 at 2.15's start: the rewording is net zero);
  * one pin per new rule sentence: the router's DDA sentence, "4c." pointing at the one table,
    the no-DE rule for both engines, the gluc line (Ruling 105.3), the Stellar rule (105.2), FRAN's
    dda_search, the quantification record's new name, log_decision --from-json in step 4.
"""
import glob
import os
import re
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)

DDA_TO_SAGE = re.compile(r"\bDDA\s*(?:→|->|=>|-->)\s*Sage\b", re.I)
WHICH_DDA = re.compile(r"timsTOF|FragPipe|unclassified|other DDA|Orbitrap", re.I)
CLIENT_TYPE = re.compile(r"\b(commercial|for-profit|non-profit|nonprofit|academic|industry|"
                         r"company|companies|client type)\b", re.I)


def read(*parts):
    with open(os.path.join(SKILL, *parts), encoding="utf-8") as fh:
        return fh.read()


def flat(text):
    return " ".join(text.split())


def section(name, heading):
    text = read("references", name)
    i = text.index(heading)
    j = text.find("\n## ", i + len(heading))
    return text[i:j if j > 0 else None]


def text_files():
    yield "SKILL.md"
    for pat in ("references/*.md", "scripts/*.py", "scripts/*.sh", "scripts/*.R"):
        for p in sorted(glob.glob(os.path.join(SKILL, pat))):
            yield os.path.relpath(p, SKILL)


class NoBareDdaToSage(unittest.TestCase):
    def test_every_dda_to_sage_says_which_dda(self):
        bad = []
        for rel in text_files():
            for n, line in enumerate(read(rel).splitlines(), 1):
                if DDA_TO_SAGE.search(line) and not WHICH_DDA.search(line):
                    bad.append(f"{rel}:{n}: {line.strip()}")
        self.assertEqual(bad, [], "a bare 'DDA -> Sage' (no timsTOF/FragPipe/Orbitrap beside it)")


class TheRouter(unittest.TestCase):
    def setUp(self):
        self.text = read("SKILL.md")

    def test_under_500_lines(self):
        self.assertLessEqual(self.text.count("\n"), 499)

    def route_paragraph(self):
        i = self.text.index("#### 4a. / 4b. / 4c. Engine routes and the bake-off")
        return flat(self.text[i:self.text.index("### 5.", i)])

    def test_the_dda_rule_sentence(self):
        p = self.route_paragraph()
        # Ruling 120: a non-specific Orbitrap DDA search routes to FragPipe too
        self.assertIn("**DDA: timsTOF DDA-PASEF and non-specific Orbitrap DDA (`--enzyme none`/"
                      "pepsin at step 4) → FragPipe (Sage where it cannot run, or on request); "
                      "other DDA → Sage unless FragPipe is asked for.**", p)
        self.assertIn('search-engines.md "4a."/"4c."', p)
        self.assertIsNone(CLIENT_TYPE.search(p))

    def test_step_4_flags_selection_sentence_and_its_log(self):
        i = self.text.index("### 4. Resolve the defaults")
        s4 = flat(self.text[i:self.text.index("#### 4a.", i)])
        self.assertIn("[--tools <tools.json>] [--ids-only] [--enzyme <name|none>] [--semi-specific]", s4)
        self.assertIn("say `selection.sentence` as written when it is not null", s4)
        self.assertIn("`log_decision.py --from-json ./wf/workflow.manifest.json`", s4)
        self.assertIn('never say "read from the data" about numbers the user gave', s4)

    def test_step_5_pins_the_manifests_engine(self):
        self.assertIn("`PIN_ENGINE=<the manifest's engine> PIN_VERSION=<its version>", self.text)


class TheReferences(unittest.TestCase):
    def test_4c_points_at_the_one_table_without_client_type_text(self):
        s = flat(section("search-engines.md", "## 4c. DDA: the routes (step 4c)"))
        self.assertIn("The rule is `engine_choice.ROUTES` (one table;", s)
        self.assertIn("**timsTOF DDA-PASEF → FragPipe**", s)
        # Ruling 120: the non-specific Orbitrap row, and Sage for the rest
        self.assertIn("**Orbitrap DDA with no enzyme rule → FragPipe** (Ruling 120)", s)
        self.assertIn("`engine_choice.NONSPECIFIC_ROUTES`", s)
        self.assertIn("**Other Orbitrap and other DDA → Sage**", s)
        self.assertIn("go only to the staff record", s)
        self.assertIsNone(CLIENT_TYPE.search(s), CLIENT_TYPE.search(s))
        routing = flat(section("search-engines.md", "## Routing"))
        self.assertIn("**DIA → DIA-NN; timsTOF DDA-PASEF, and Orbitrap DDA with no enzyme rule "
                      "(`--enzyme none` or pepsin), → FragPipe, with Sage as the automatic "
                      "fallback where FragPipe cannot run; other Orbitrap and other DDA → Sage, "
                      "FragPipe when it is asked for.**", routing)
        self.assertNotIn("opt-in", routing)

    def test_the_licence_table_is_unchanged(self):
        # Ruling 56: no new licence or commercial-use text on the FragPipe row
        s = read("references", "search-engines.md")
        self.assertIn("| **FragPipe** (MSFragger/IonQuant) | license-gated | depends | requires the "
                      "user's own license |", s)

    def test_the_client_sentence_has_no_client_type(self):
        import engine_choice
        sel, _ = engine_choice.choose("DDA", "timstof", env={"engines": {"fragpipe": {
            "available": False, "note": "no macOS build"}}})
        self.assertEqual(sel["how"], "fallback")
        self.assertIsNone(CLIENT_TYPE.search(sel["sentence"]))

    def test_no_de_rule_for_both_engines(self):
        s = flat(read("references", "no-de-route.md"))
        self.assertIn("One rule, both DDA engines: a no-DE DDA search is identifications only "
                      "unless the user", s)
        self.assertIn("A DIA search (DIA-NN, Radiant, or FragPipe's DIA route) always quantifies: "
                      "`--ids-only` is for DDA searches only (Sage or FragPipe), and step 4 "
                      "refuses it otherwise -- leave it out; the search's own report is the "
                      "result.", s)

    def test_gluc_and_stellar_lines(self):
        s = flat(read("references", "parameters.md"))
        self.assertIn("`--enzyme gluc` searches FragPipe's `gluc` rule (cut after D and E, not "
                      "before P: the phosphate-buffer specificity); for a Glu-C digest in ammonium "
                      "bicarbonate, which cuts after E only, give the rule with `--overrides` "
                      "(Sage: `{\"database\": {\"enzyme\": {\"cleave_at\": \"E\"}}}`; FragPipe: "
                      "`{\"msfragger.search_enzyme_cut_1\": \"E\"}`).", s)
        self.assertIn("is searched for **identifications only** (`--ids-only`)", s)
        self.assertIn("targeted quantification (PRM, with Skyline: skyline.md)", s)
        self.assertIn("the non-specific workflows keep FragPipe's own 3", s)

    def test_fragpipe_dia_quantifies_with_its_bundled_diann_said_once(self):
        # Ruling 107 M4: the skill never points FragPipe at a DIA-NN 2.x (run_search passes the
        # bundled one; each search reads its own settings folder). One statement, pointed at.
        rule = "FragPipe's DIA quantification is its bundled DIA-NN"
        se = read("references", "search-engines.md")
        self.assertEqual(se.count(f"**{rule} (2.15; the one statement of this rule).**"), 1)
        for name in ("fragpipe-diatracer.md", "cross-tool-comparison.md"):
            text = flat(read("references", name))
            self.assertIn(f'`references/search-engines.md`, "{rule}"', text, name)
            self.assertNotIn("fragpipe_bootstrap.py --diann <2.x", text, name)
            self.assertNotIn("passes no `--config-diann`", text, name)
        self.assertNotIn("Do not choose a DIA-NN 2.x for FragPipe yourself", se)

    def test_fran_says_dda_search(self):
        s = read("references", "fran.md")
        self.assertIn("| `dda_search` | a DDA search by an engine that also searches DIA "
                      "(FragPipe DDA, DIA-NN `--dda`)", s)

    def test_the_quantification_record_is_named_by_its_new_name(self):
        for rel in ("SKILL.md",) + tuple(os.path.relpath(p, SKILL) for p in
                                         glob.glob(os.path.join(SKILL, "references", "*.md"))):
            for n, line in enumerate(read(rel).splitlines(), 1):
                if "sage_lfq_check.json" in line:
                    self.assertRegex(line, r"2\.14\.x|still read", f"{rel}:{n}")
        self.assertIn("quantification_check.json", read("references", "outputs.md"))


if __name__ == "__main__":
    unittest.main()
