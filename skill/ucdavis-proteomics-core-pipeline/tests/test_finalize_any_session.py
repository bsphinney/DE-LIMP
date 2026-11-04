#!/usr/bin/env python3
"""`session.py finalize` for ANY session (2.14.2, Amendment 2 to Task 13; Brett, 2026-10-31).

The beer re-analysis (species composition, no DE) and three earlier HT submissions (search QC)
never ran finalize: the agent read it as the DE path's last step, because finalize "reads
de_provenance.json and writes the DE README". It hand-wrote its README.md and AGENTS.md instead,
avoided `session.py docs` because docs would overwrite them, and kept its own analysis code in a
second folder outside the session. What these tests pin:
  * a session with no DE -- a search-only one, one whose search is in output/search_dda/, one
    with no search at all -- still gets AGENTS.md, the README, MANIFEST.txt and the Methods, with
    a [SKIPPED] line for each part only a DE makes (never a silent gap, DE-LIMP rule 4);
  * finalize --report/--methods/--code-dir/--extra-dir fold work done outside the session in,
    under output/custom/, with checksums in the bundle's manifest and a line in AGENTS.md and the
    README; a report the skill did not write goes in only under --override-finish;
  * a README.md or AGENTS.md written by hand (no skill marker) is never changed: the skill's
    copy is README.skill.md / AGENTS.skill.md, and the session still finalizes and stamps.
Stdlib only, no network.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import report_files as rf                      # noqa: E402
import session                                 # noqa: E402
import test_finish_gates as tfg                # noqa: E402  (FinalizeBase: a no-DE DIA session)

read, write = tfg.read, tfg.write


def skipped(manifest, part):
    """The reason of MANIFEST.txt's [SKIPPED] line for `part`, or None."""
    for ln in manifest.splitlines():
        if ln.startswith("[SKIPPED]") and ln[10:].startswith(part) and " -- " in ln:
            return ln.split(" -- ", 1)[1]
    return None


class AnySessionType(tfg.FinalizeBase):
    """Amendment 2 item 1: finalize works for ANY session type."""

    def setUp(self):
        super().setUp()
        self.finish()

    def done(self, r):
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        for rel in ("AGENTS.md", "README.md", "README.html", "MANIFEST.txt", "output/methods.md"):
            self.assertTrue(os.path.isfile(os.path.join(self.sd, rel)), rel)
        st = json.loads(read(self.stamp_path()))
        self.assertEqual(st["outputs"]["AGENTS.md"], rf._sha256(os.path.join(self.sd, "AGENTS.md")))
        return st

    def test_a_search_only_session_has_a_skipped_line_for_each_de_only_part(self):
        self.done(self.finalize("--zip"))
        man = self.manifest()
        self.assertTrue(session.DE_ONLY_PARTS)
        for part in session.DE_ONLY_PARTS + (session.PEAK_CHECK_PART,):
            why = skipped(man, part)
            self.assertIsNotNone(why, f"no [SKIPPED] line for {part}:\n{man}")
            self.assertIn("no differential expression in this session", why)
            self.assertIn("no output/tables/de_provenance.json", why)
        self.assertIsNone(skipped(man, session.SEARCH_RECORD_PART), "its search is found")
        # client words, not a pointer at a missing record (C2 biologist L5)
        self.assertIn("Differential expression: none in this analysis", read(self.p["readme"]))
        # L5b review L2: nothing in the README or AGENTS.md describes a DE this session has not
        readme, agents = read(self.p["readme"]), read(os.path.join(self.sd, "AGENTS.md"))
        for doc, name in ((readme, "README.md"), (agents, "AGENTS.md")):
            for words in ("DE results (DE_*.csv)", "reproducibility_log.R (the analysis as R)",
                          "the database, the DE and", "(LC-MS + search + DE + acknowledgment)",
                          "the contrasts in `output/tables/de_provenance.json`"):
                self.assertNotIn(words, doc, f"{name}: {words}")
        self.assertIn("LC-MS acquisition, the search, the database and the UC Davis "
                      "instrument-grant acknowledgment", readme)
        self.assertIn("Publication Methods (LC-MS + search + acknowledgment)", agents)
        self.assertIn("this session has no differential expression", agents)

    def test_a_session_with_a_de_has_no_de_only_lines(self):
        """(pin) a DE session's MANIFEST.txt says nothing about missing DE parts."""
        write(os.path.join(self.p["de_dir"], "de_provenance.json"),
              json.dumps(tfg.tdp.DE_PROV))
        write(os.path.join(self.p["de_dir"], "DE_Treated_vs_Control.csv"), "Protein,logFC\nP1,1\n")
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        man = self.manifest()
        for part in session.DE_ONLY_PARTS:
            self.assertIsNone(skipped(man, part), part)
        # ... and its README / AGENTS.md describe the DE (L5b review L2, the other side)
        readme, agents = read(self.p["readme"]), read(os.path.join(self.sd, "AGENTS.md"))
        self.assertIn("the database, the DE and the UC Davis", readme)
        self.assertIn("(LC-MS + search + DE + acknowledgment)", agents)
        self.assertIn("the contrasts in `output/tables/de_provenance.json`", agents)

    def test_a_search_outside_output_search_is_the_one_described(self):
        # A Core member's submission: the searches were in output/search_dda/ (and output/search_sage/)
        dda = os.path.join(self.out, "search_dda")
        shutil.move(self.p["search_out"], dda)
        st = self.done(self.finalize("--zip"))
        self.assertTrue(st["outputs"]["output/methods.md"])
        readme, agents = read(self.p["readme"]), read(os.path.join(self.sd, "AGENTS.md"))
        for doc in (readme, agents):
            self.assertIn("Search engine: DIA-NN 2.7.0", doc)
            self.assertNotIn("the search's own record is not here", doc)
            self.assertNotIn("search_provenance.json names no output folder", doc)
        self.assertIn("| Search engine, the version that ran, exact command | "
                      "`output/search_dda/search_provenance.json` |", agents)
        self.assertIn("| `output/search_dda` | the search output", readme)
        methods = read(os.path.join(self.out, "methods.md"))
        self.assertIn("DIA-NN 2.7.0", methods)
        self.assertNotIn("[pinned version — confirm it is what ran]", methods)
        self.assertIsNone(skipped(self.manifest(), session.SEARCH_RECORD_PART))
        names = zipfile.ZipFile(self.sd + ".zip").namelist()
        self.assertTrue(any(n.endswith("output/search_dda/report.parquet") for n in names), names)
        self.assertFalse([n for n in names if n.endswith(".quant")])

    def test_a_raw_list_comes_from_a_search_outside_output_search(self):
        # a hive_remote session is initialised without --raw: the list comes from the search
        os.remove(self.p["raw_list"])
        dda = os.path.join(self.out, "search_dda")
        shutil.move(self.p["search_out"], dda)
        prov = json.loads(read(os.path.join(dda, "search_provenance.json")))
        prov["files"] = ["/data/runs/HeLa_ctrl_01.d", "/data/runs/HeLa_ctrl_02.d"]
        write(os.path.join(dda, "search_provenance.json"), json.dumps(prov))
        self.done(self.finalize())
        self.assertIn("/data/runs/HeLa_ctrl_01.d", read(self.p["raw_list"]))
        self.assertIn("written from output/search_dda/search_provenance.json", self.manifest())

    def test_a_session_with_no_search_still_gets_its_documents_and_says_so(self):
        # no search at all here (a re-analysis of published results, a search run elsewhere):
        # finalize writes what it can; part a's stamp rules decide the stamp -- AGENTS.md, the
        # README and the Methods were written, so it stamps
        shutil.rmtree(self.p["search_out"])
        self.done(self.finalize())
        why = skipped(self.manifest(), session.SEARCH_RECORD_PART)
        self.assertIsNotNone(why, self.manifest())
        self.assertIn("no search record in this session", why)
        self.assertIn("Search engine: DIA-NN 2.7.0 (workflow manifest pin -- the search's own "
                      "record is not here)", read(self.p["readme"]))
        self.assertIn("[pinned version — confirm it is what ran]",
                      read(os.path.join(self.out, "methods.md")))

    def test_two_searches_and_no_tie_break_say_which_one_the_documents_describe(self):
        # output/search (the session's own) and a second one, and no DE input to choose by:
        # one answer everywhere (T44 review F3, Ruling 78) -- the README, AGENTS.md and the
        # Methods describe neither, never output/search's by default, and MANIFEST.txt says so
        other = os.path.join(self.out, "search_sage")
        shutil.copytree(self.p["search_out"], other)
        self.done(self.finalize())
        why = skipped(self.manifest(), session.SEARCH_RECORD_PART)
        self.assertIsNotNone(why, self.manifest())
        self.assertIn("2 search records", why)
        self.assertIn("output/search_sage/search_provenance.json", why)
        self.assertIn("describe none of them", why)


if __name__ == "__main__":
    unittest.main()
