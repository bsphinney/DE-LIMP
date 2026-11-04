#!/usr/bin/env python3
"""A README.md or AGENTS.md written by hand is never changed (2.14.2, Amendment 2 to Task 13,
item 3).

The beer re-analysis avoided `session.py docs` because it would OVERWRITE the README.md it had
written by hand, so it hand-wrote AGENTS.md too and never finalized. The rule (report_files):
a README.md, README.html or AGENTS.md without the skill's marker (SKILL_DOC_MARKER, or the words
a pre-2.14.2 skill document carries) is kept byte for byte; the skill's own copy goes beside it
as README.skill.md / README.skill.html / AGENTS.skill.md, each saying so and linking the other,
and MANIFEST.txt records it. finalize still stamps (the skill's copies are what it requires), the
zip carries both, `session.py status` reads the skill's copies, `session.py docs` follows the
same rule, and so does `make_podcast.py link`.
Stdlib only, no network.
"""
import json
import os
import subprocess
import sys
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import report_files as rf                      # noqa: E402
import test_finish_gates as tfg                # noqa: E402  (FinalizeBase: a no-DE DIA session)
import test_make_podcast as tmp_               # noqa: E402  (the podcast workspace)

read, write = tfg.read, tfg.write
PY = sys.executable
HAND_README = "# Beer species composition\n\nWritten by hand: what we did and why.\n"
HAND_AGENTS = "# AGENTS\n\nThe species table is species.csv. Do not re-run the search.\n"


def raw_bytes(path):
    with open(path, "rb") as fh:
        return fh.read()


class HandWrittenDocs(tfg.FinalizeBase):
    def setUp(self):
        super().setUp()
        self.finish()
        self.readme = write(self.p["readme"], HAND_README)            # over init's starter
        self.agents = write(os.path.join(self.sd, "AGENTS.md"), HAND_AGENTS)

    def at(self, name):
        return os.path.join(self.sd, name)

    def docs(self):
        return subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "docs", "--dir", self.sd],
                              capture_output=True, text=True, env=self.env)

    def status(self):
        return subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "status", "--session",
                               self.sd, "--json"], capture_output=True, text=True, env=self.env)

    def test_finalize_keeps_them_writes_the_skills_beside_them_and_stamps(self):
        before = raw_bytes(self.readme), raw_bytes(self.agents)
        r = self.finalize("--zip")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual((raw_bytes(self.readme), raw_bytes(self.agents)), before, "byte-identical")
        skill_readme, skill_agents = read(self.at("README.skill.md")), read(self.at("AGENTS.skill.md"))
        for path in ("README.skill.md", "AGENTS.skill.md", "README.html"):
            self.assertTrue(rf.is_skill_doc(self.at(path)), path)
        self.assertFalse(rf.is_skill_doc(self.readme))
        # each says what it is beside, and links the other
        self.assertIn("`README.md` in this folder was written by hand", skill_readme)
        self.assertIn("- `AGENTS.skill.md` — for an AI assistant", skill_readme)  # no link: T33
        self.assertIn("`AGENTS.md` in this folder was written by hand", skill_agents)
        self.assertIn("README.skill.md", skill_agents)
        self.assertIn("**Open `README.html`**", skill_readme)           # no hand-made README.html
        man = self.manifest()
        self.assertRegex(man, r"\[INFO\] +README\.md \(written by hand\) +-- kept as written: it "
                              r"has no skill marker; the skill's README is README\.skill\.md")
        self.assertRegex(man, r"\[INFO\] +AGENTS\.md \(written by hand\) +-- kept as written: it "
                              r"has no skill marker; the skill's guide is AGENTS\.skill\.md")
        st = json.loads(read(self.stamp_path()))
        for rel in ("README.skill.md", "AGENTS.skill.md", "README.html", "README.md", "AGENTS.md"):
            self.assertEqual(st["outputs"][rel], rf._sha256(self.at(rel)), rel)
        self.assertEqual(rf.stamp_drift(self.sd, st), [])
        names = {n.split("/", 1)[1] for n in zipfile.ZipFile(self.sd + ".zip").namelist()}
        for rel in ("README.md", "AGENTS.md", "README.skill.md", "AGENTS.skill.md", "README.html"):
            self.assertIn(rel, names)
        r = self.status()
        self.assertEqual(r.returncode, 0, r.stdout)
        # the hand-written README is part of what was stamped: edited later, the stamp is stale
        write(self.readme, HAND_README + "\nOne more line.\n")
        items = {i["item"]: i for i in json.loads(self.status().stdout)["items"]}
        self.assertEqual(items["finalized"]["state"], "stale")
        self.assertIn("README.md", items["finalized"]["detail"])

    def test_docs_follows_the_same_rule(self):
        before = raw_bytes(self.readme)
        r = self.docs()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(raw_bytes(self.readme), before)
        self.assertEqual(read(self.agents), HAND_AGENTS)
        out = json.loads(r.stdout)
        self.assertEqual(out["readme"], self.at("README.skill.md"))
        self.assertEqual(out["agents"], self.at("AGENTS.skill.md"))
        self.assertTrue(any(x["part"] == "README.md (written by hand)" for x in out["parts"]))

    def test_status_wants_the_skills_copies_when_the_folder_has_hand_written_ones(self):
        items = {i["item"]: i for i in json.loads(self.status().stdout)["items"]}
        for item in ("agents", "readme"):
            self.assertFalse(items[item]["done"], items[item])
            self.assertIn("written by hand", items[item]["detail"])

    def test_a_hand_made_readme_html_is_kept_too(self):
        html = write(self.at("README.html"), "<html>by hand</html>")
        self.assertEqual(self.finalize().returncode, 0)
        self.assertEqual(read(html), "<html>by hand</html>")
        self.assertTrue(rf.is_skill_doc(self.at("README.skill.html")))
        self.assertIn("**Open `README.skill.html`**", read(self.at("README.skill.md")))


class TheSkillsOwnDocs(tfg.FinalizeBase):
    """(pin, and the marker) the skill's own documents -- init's starter README, a README or
    AGENTS.md from an earlier finalize, from before 2.14.2 too -- are written over as before."""

    def setUp(self):
        super().setUp()
        self.finish()

    def test_init_marks_its_starter_readme_and_finalize_replaces_it(self):
        self.assertTrue(rf.is_skill_doc(self.p["readme"]), read(self.p["readme"]))
        self.assertEqual(self.finalize().returncode, 0)
        self.assertFalse(os.path.exists(os.path.join(self.sd, "README.skill.md")))
        self.assertIn("## Where this lives on HIVE", read(self.p["readme"]))
        # the marker is the last line: README.md's first stays `# <name>` (fran_deposit reads it)
        self.assertTrue(read(self.p["readme"]).startswith("# 2026-09-24_demo\n"))
        self.assertTrue(read(self.p["readme"]).rstrip().endswith(rf.SKILL_DOC_MARKER))
        self.assertIn(rf.SKILL_DOC_MARKER, read(os.path.join(self.sd, "README.html")))
        self.assertTrue(read(os.path.join(self.sd, "AGENTS.md")).rstrip().endswith(
            rf.SKILL_DOC_MARKER))
        self.assertEqual(self.finalize().returncode, 0)                   # again: written over
        self.assertFalse(os.path.exists(os.path.join(self.sd, "README.skill.md")))

    def test_documents_from_before_the_marker_are_the_skills(self):
        # a session from before 2.14.2 has no record of what the skill wrote (M1)
        os.remove(os.path.join(self.sd, rf.SKILL_DOCS_RECORD))
        old = {"README.md": "# 2026-09-24_demo\n\n- Date: 2026-09-24\n- Status: in progress\n\n"
                            "Layout: `input/` (conditions, FASTA, params), `output/` (search, "
                            "tables, figures, reproducibility, report), `scripts/`, `logs/`.\n",
               "AGENTS.md": "# AGENTS.md — guide to this folder for an AI agent\n\nGenerated by "
                            "`session.py` from this session's own records. Every fact here.\n",
               "README.html": "<p><em>Written 2026-09-30 by <code>session.py</code> from this "
                              "session's records.</em></p>"}
        for name, text in old.items():
            write(os.path.join(self.sd, name), text)
            self.assertTrue(rf.is_skill_doc(os.path.join(self.sd, name)), name)
        write(os.path.join(self.sd, "README.md"),
              "# x\n\n*Written 2026-09-30 by `session.py` from this session's records.*\n")
        self.assertTrue(rf.is_skill_doc(os.path.join(self.sd, "README.md")))
        self.assertEqual(self.finalize().returncode, 0)
        for name in old:
            self.assertIn(rf.SKILL_DOC_MARKER, read(os.path.join(self.sd, name)), name)
        for name in ("README.skill.md", "AGENTS.skill.md", "README.skill.html"):
            self.assertFalse(os.path.exists(os.path.join(self.sd, name)), name)


class EditedSkillDocs(tfg.FinalizeBase):
    """L5b review M1: a document the skill wrote and someone then EDITED is hand-written: init's
    starter README.md with notes added (the beer and HT sessions, held for 2.14.2), a 2.14.1
    starter (no marker, no record) edited, or a finalized README.md edited. The skill records
    the sha256 of each document it writes (report_files.record_skill_doc, logs/skill_docs.json);
    one that differs is kept, and the skill's copy goes beside it. A starter with no record is the
    skill's only while it holds nothing but the starter's lines."""

    NOTES = "\n## My own notes (hand-written)\nSpecies composition was done in ~/beer_local.\n"

    def setUp(self):
        super().setUp()
        self.finish()

    def test_a_hand_edited_starter_readme_survives_finalize(self):
        self.assertIn(rf.STARTER_README_SIGN, read(self.p["readme"]))
        with open(self.p["readme"], "a", encoding="utf-8") as fh:
            fh.write(self.NOTES)
        before = raw_bytes(self.p["readme"])
        self.assertFalse(rf.is_skill_doc(self.p["readme"]))
        r = self.finalize("--zip")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(raw_bytes(self.p["readme"]), before, "kept byte for byte")
        self.assertTrue(rf.is_skill_doc(os.path.join(self.sd, "README.skill.md")))
        self.assertRegex(self.manifest(), r"\[INFO\] +README\.md \(written by hand\)")

    def test_a_pre_2142_starter_edited_by_hand_survives_finalize(self):
        os.remove(os.path.join(self.sd, rf.SKILL_DOCS_RECORD))     # 2.14.1 recorded nothing
        txt = read(self.p["readme"]).replace(rf.SKILL_DOC_MARKER, "")
        write(self.p["readme"], txt + "\n## Hand-written results summary\nSee my notebook.\n")
        self.assertEqual(self.finalize().returncode, 0)
        self.assertIn("Hand-written results summary", read(self.p["readme"]))
        self.assertTrue(os.path.isfile(os.path.join(self.sd, "README.skill.md")))

    def test_an_untouched_pre_2142_starter_is_still_the_skills(self):
        """(pin) the starter's own lines only, with or without the marker and the re-analysis
        line: replaced by finalize, as before."""
        os.remove(os.path.join(self.sd, rf.SKILL_DOCS_RECORD))
        for extra in ("", "- **Re-analysis of:** `/x/prior` — see `DIFFERENCES.md` (written at "
                          "finalize) for exactly what changed.\n"):
            write(self.p["readme"], "# demo\n\n- Date: 2026-09-24\n- Status: in progress\n"
                                    + extra + "\n" + rf.STARTER_README_LAYOUT + "\n")
            self.assertTrue(rf.is_skill_doc(self.p["readme"]), extra)
        write(self.p["readme"], read(self.p["readme"]) + "\n" + rf.SKILL_DOC_MARKER + "\n")
        self.assertTrue(rf.is_skill_doc(self.p["readme"]))
        write(self.p["readme"], read(self.p["readme"]) + "# A second title\n")
        self.assertFalse(rf.is_skill_doc(self.p["readme"]))

    def test_a_finalized_readme_edited_afterwards_is_kept_by_the_next_finalize(self):
        self.assertEqual(self.finalize().returncode, 0)
        for name in ("README.md", "AGENTS.md"):
            self.assertTrue(rf.is_skill_doc(os.path.join(self.sd, name)), name)
        agents = os.path.join(self.sd, "AGENTS.md")
        write(agents, read(agents) + "\n- Also: the species table is in ~/beer_local.\n")
        before = raw_bytes(agents)
        self.assertEqual(self.finalize().returncode, 0)
        self.assertEqual(raw_bytes(agents), before)
        self.assertTrue(os.path.isfile(os.path.join(self.sd, "AGENTS.skill.md")))
        self.assertFalse(os.path.exists(os.path.join(self.sd, "README.skill.md")),
                         "README.md was not edited: still the skill's")

    def test_a_skill_copy_left_from_an_earlier_finalize_is_not_zipped(self):
        """N4: README.skill.md from a finalize when README.md was hand-written; the hand README.md
        is then removed, so README.md is the skill's again -- the old copy stays out of the zip."""
        write(self.p["readme"], HAND_README)
        self.assertEqual(self.finalize().returncode, 0)
        self.assertTrue(os.path.isfile(os.path.join(self.sd, "README.skill.md")))
        os.remove(self.p["readme"])
        r = self.finalize("--zip")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        names = {n.split("/", 1)[1] for n in zipfile.ZipFile(self.sd + ".zip").namelist()}
        self.assertIn("README.md", names)
        self.assertNotIn("README.skill.md", names)
        self.assertIn("README.skill.md", json.dumps(json.loads(r.stdout)["zip_excluded"]))


class PodcastLink(tmp_.Workspace):
    """make_podcast.py link puts its line into the skill's README and AGENTS.md: beside a
    hand-written one, into the .skill.md copy, never the hand-written file."""

    def test_link_leaves_a_hand_written_readme_and_agents_alone(self):
        hand = write(os.path.join(self.d, "README.md"), HAND_README)
        hand_agents = write(os.path.join(self.d, "AGENTS.md"), HAND_AGENTS)
        skill = write(os.path.join(self.d, "README.skill.md"),
                      rf.SKILL_DOC_MARKER + "\n# x\n\n## Start here\n\n- [Analysis report]"
                      "(output/Analysis_Report.html) — the report\n")
        skill_agents = write(os.path.join(self.d, "AGENTS.skill.md"),
                             rf.SKILL_DOC_MARKER + "\n# AGENTS.md\n\n## Do not\n\n- x\n")
        tfg.skill_report(self.out)
        tmp_.checked_podcast(self.out, self.report)
        rc, _, err = tmp_.run("link", self.out)
        self.assertEqual(rc, 0, err)
        self.assertEqual((read(hand), read(hand_agents)), (HAND_README, HAND_AGENTS))
        self.assertIn("podcast", read(skill).lower())
        self.assertIn("podcast", read(skill_agents).lower())

    def test_links_edit_of_the_skills_readme_keeps_it_the_skills(self):
        """M1: link's own edit is the skill's; it records it, so the next finalize still writes
        README.md and AGENTS.md (never mistakes the link line for a hand edit)."""
        readme = write(os.path.join(self.d, "README.md"),
                       "# x\n\n## Start here\n\n- [Analysis report](output/Analysis_Report.html) — "
                       "the report\n\n" + rf.SKILL_DOC_MARKER + "\n")
        agents = write(os.path.join(self.d, "AGENTS.md"),
                       "# AGENTS.md\n\n## Do not\n\n- x\n\n" + rf.SKILL_DOC_MARKER + "\n")
        for path in (readme, agents):
            rf.record_skill_doc(path, "session_docs.py")
        tfg.skill_report(self.out)
        tmp_.checked_podcast(self.out, self.report)
        rc, _, err = tmp_.run("link", self.out)
        self.assertEqual(rc, 0, err)
        self.assertIn("podcast", read(readme).lower())
        for path in (readme, agents):
            self.assertTrue(rf.is_skill_doc(path), path)


if __name__ == "__main__":
    unittest.main()
