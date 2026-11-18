#!/usr/bin/env python3
"""The instructions agree with themselves and with the scripts (2.15.1 design audit, lane 4).

Each check pins one finding of the 2026-11-17 design audit that the text now gets right:
  * C1  a file that stalls twice is the user's decision: no reference tells the agent to drop
        it or to edit file_list.txt / the chain's run list by hand;
  * C3  step 3 asks about re-injections for everyone, conditions.md documents the
        `technical_replicates` key, and that key really writes the Sample column run_de.R
        blocks on (run here, end to end);
  * C4  a combined CoreOmics session names another submission's PI only for the same lab;
  * C5  no reference calls a group-readable STAN credential a fix;
  * M1/M2 the licence question sits in step 3 (before the confirmation and the download), and
        no reference offers AlphaDIA, which steps 4 and 6b refuse (checked against the script);
  * M7  step 8 takes the DE method from the manifest, and the text's engine -> method pairs are
        engine_choice.de_method's;
  * M8  "chymotrypsin changes the engine" is gone, and the route still agrees (engine_choice);
  * M13 hive_local, which check_access.sh emits, is defined in SKILL.md;
  * M22 the Radiant preset file the references name is the one resolve_defaults.py writes;
  * M3  the Stop-hook sentence says it does not see a session on HIVE;
  * M11/M29 the router says what the skill does not do, and the description names exclusions
        and stays within the frontmatter limit.
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
sys.path.insert(0, SCRIPTS)


def read(*parts):
    with open(os.path.join(SKILL, *parts), encoding="utf-8") as fh:
        return fh.read()


def flat(text):
    return " ".join(text.split())


def step(text, head, nxt):
    i = text.index(head)
    return flat(text[i:text.index(nxt, i)])


SKILL_MD = read("SKILL.md")


class StalledFileIsTheUsersDecision(unittest.TestCase):          # C1, M5
    DROP = re.compile(r"\bdrop (?:it|that (?:one )?file)\b|remove it from step 3", re.I)

    def test_no_reference_drops_a_file_on_its_own(self):
        for ref in ("watcher.md", "ht-submissions.md"):
            text = read("references", ref)
            with self.subTest(ref=ref):
                self.assertIsNone(self.DROP.search(text), self.DROP.search(text))
                self.assertNotIn("is *success*", text)

    def test_the_playbook_asks_and_never_edits_the_run_list(self):
        w = flat(read("references", "watcher.md"))
        self.assertIn("**If it stalls again, stop and ask the user**", w)
        self.assertIn("never edit `file_list.txt`, a step's `--f` list or the run count", w)
        self.assertNotIn("When resubmitting after a fix, edit the sbatch", w)

    def test_7b_lists_the_exceptions(self):
        s7b = step(SKILL_MD, "### 7b.", "### 7c.")
        for want in ("`cancelled` (find out who first)", "`held`", "a failed run check",
                     "a file that stalled twice (ask the user"):
            self.assertIn(want, s7b)
        self.assertIn("record_run.py search-done --out <dir> --status failed", s7b)


class TechnicalReplicates(unittest.TestCase):                    # C3
    def test_step_3_asks_everyone(self):
        s3 = step(SKILL_MD, "### 3. ", "### 3b.")
        self.assertIn('"are any runs repeat injections of one sample?"', s3)
        self.assertIn("`technical_replicates`", s3)

    def test_conditions_md_documents_the_key(self):
        c = flat(read("references", "conditions.md"))
        self.assertIn('**"Are any runs repeat injections of the same sample?"**', c)
        self.assertIn('`{"technical_replicates": {"<run>": "<sample>", ...}}`', c)

    def test_the_documented_key_writes_the_sample_column(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "conditions.csv")
            intent = {"groups": {"control": ["c1", "c2", "c2b"], "treated": ["t1", "t2", "t3"]},
                      "technical_replicates": {"c2": "c2", "c2b": "c2"}}
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "collect_conditions.py"),
                                "--map", out, "--runs", "c1,c2,c2b,t1,t2,t3",
                                "--mapping-json", json.dumps(intent)],
                               capture_output=True, text=True, cwd=tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            res = json.loads(r.stdout)
            self.assertEqual(res["technical_replicates"]["column"], "Sample")
            self.assertEqual(res["technical_replicates"]["samples"], {"c2": ["c2", "c2b"]})
            with open(out) as fh:
                rows = fh.read().splitlines()
            self.assertEqual(rows[0], "File.Name,Group,Sample")
            self.assertIn("c2b,control,c2", rows)


class CombinedSubmissions(unittest.TestCase):                    # C4
    def test_another_pi_is_named_only_for_the_same_lab(self):
        c = flat(read("references", "core-submissions.md"))
        i = c.index("## Several submissions in one analysis")
        sec = c[i:c.index("## Core client work", i)]
        self.assertIn("**Same lab only, by default.**", sec)
        self.assertIn("its number, its PI and the samples it brought when it is the same lab; "
                      "only the number of samples it brought when it is not", sec)
        self.assertIn("another lab never gets it", sec)


class StanCredential(unittest.TestCase):                         # C5, M18
    def test_no_group_readable_token_is_called_a_fix(self):
        h = flat(read("references", "ht-submissions.md"))
        self.assertNotIn("durable fix for the CLI path is a group-readable token", h)
        self.assertIn("**Never copy it, and never make it group-readable**", h)

    def test_the_share_token_never_passes_through_the_agent(self):
        h = read("references", "ht-submissions.md")
        self.assertNotIn("printf '%s\\n' '<tok>'", h)
        self.assertIn("read -rs T && printf '%s\\n' \"$T\" > ~/.stan/share_0807; unset T", h)


class Licensing(unittest.TestCase):                              # M1, M2, m19
    def test_the_question_is_asked_at_step_3_not_step_7(self):
        self.assertIn('**"Academic/non-profit, or commercial?"**', step(SKILL_MD, "### 3. ", "### 3b."))
        self.assertNotIn("commercial?", step(SKILL_MD, "### 7. ", "### 7b."))

    def test_alphadia_is_never_offered(self):
        s = flat(read("references", "search-engines.md"))
        for bad in ("route to AlphaDIA", "offer AlphaDIA", "commercial → AlphaDIA",
                    "(`--engine alphadia`)"):
            self.assertNotIn(bad, s)
        self.assertIn("AlphaDIA is not a skill route", s)

    def test_step_4_really_refuses_alphadia(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                                "--acquisition", "DIA", "--instrument", "timsTOF HT",
                                "--engine", "alphadia", "--dest", tmp],
                               capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("unknown engine 'alphadia'", r.stderr + r.stdout)


class DeMethodFromTheManifest(unittest.TestCase):                # M7
    def test_step_8_does_not_hardcode_dpc(self):
        s8 = step(SKILL_MD, "### 8. ", "### 8a.")
        self.assertNotIn("--method dpc", s8)
        self.assertIn("--method <the manifest's de.method>", s8)

    def test_the_pairs_named_are_engine_choices(self):
        import engine_choice
        self.assertEqual(engine_choice.de_method("diann", "DIA"), "dpc")
        self.assertEqual(engine_choice.de_method("sage", "DDA"), "maxlfq")
        self.assertEqual(engine_choice.de_method("fragpipe", "DDA"), "maxlfq")
        # 2.15.1 follow-up: Radiant's adapted report is a protein table too
        self.assertEqual(engine_choice.de_method("radiant", "DIA"), "maxlfq")
        self.assertIn("`dpc` for DIA-NN, `maxlfq` for Sage, FragPipe DDA and Radiant",
                      step(SKILL_MD, "### 8. ", "### 8a."))


class Enzyme(unittest.TestCase):                                 # M8, M9
    def test_chymotrypsin_does_not_change_the_engine(self):
        self.assertNotIn("chymotrypsin change the engine", SKILL_MD)
        import engine_choice
        sel, _ = engine_choice.choose("DDA", "orbitrap", enzyme="chymotrypsin")
        self.assertEqual(sel["engine"], "sage")

    def test_the_fasta_enzyme_default_is_not_a_reason_to_skip_asking(self):
        e = flat(read("references", "environment.md"))
        self.assertNotIn("leave it unless the user says otherwise", e)
        self.assertIn("**Always pass step 3's answer**", e)


class ModesAndHooks(unittest.TestCase):                          # M13, M3, M12
    def test_hive_local_is_defined(self):
        self.assertIn('MODE="hive_local"', read("scripts", "check_access.sh"))
        self.assertIn("`hive_local` (Claude Code already runs on HIVE", SKILL_MD)

    def test_stop_hook_limit_is_stated(self):
        self.assertIn("**for a session on this computer only**: it cannot see a session on HIVE",
                      SKILL_MD)

    def test_one_order_of_the_first_checks(self):
        self.assertIn("**The order is fixed: 0a's questions and `check_access.sh` (the mode) → "
                      "this version check → setup → the notes → 0b → 0c.**", flat(SKILL_MD))


class PresetNames(unittest.TestCase):                            # M22
    def test_the_radiant_preset_is_the_one_step_4_writes(self):
        import resolve_defaults
        name = "radiant" + resolve_defaults.PRESET_ENGINES["radiant"]
        s = read("references", "search-engines.md")
        self.assertNotIn("default.radiantConfig", s)
        self.assertIn("./wf/" + name, s)
        self.assertIn("./wf/" + name, SKILL_MD)


class Scope(unittest.TestCase):                                  # M11, M29
    def test_router_says_what_it_does_not_do(self):
        sec = flat(SKILL_MD[SKILL_MD.index("## What this skill does not do"):])
        for want in ("TMT / iTRAQ", "PTM-site work", "SCIEX `.wiff`", "Spectronaut, MaxQuant"):
            self.assertIn(want, sec)

    def test_only_two_variable_mods_are_searched(self):
        # the scope section says only Ox(M) and deamidation are searched as variable mods
        import estimate_params
        self.assertEqual(sorted(estimate_params.VAR_MODS), ["deam", "ox"])

    def test_description_names_exclusions_and_fits(self):
        fm = SKILL_MD.split("---")[1]
        desc = flat(fm.split("description: >", 1)[1])
        self.assertLessEqual(len(desc), 1024)
        self.assertIn("Not for editing the DE-LIMP Shiny app", desc)
        self.assertNotIn("Spectronaut or FragPipe", desc)      # no Spectronaut reader in the skill


if __name__ == "__main__":
    unittest.main()
