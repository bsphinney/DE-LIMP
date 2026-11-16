#!/usr/bin/env python3
"""
A session whose search was split by FAIMS group (run_search --faims-split; 2.15 lane X10).

The split writes one search folder per group (output/<out>/faims_none, .../faims_faims), each with
its own search_provenance.json carrying `faims`. That was "N search records, none chosen", the
brief said "No search output was found", and the only documented way out renamed the other group
to output/superseded_<name>/ -- calling a valid result superseded. What these pin:
  * the session says it holds a FAIMS-split search, names each group and its runs, and never says
    "No search output" or "superseded" about it;
  * `session.py choose-search` names one group WITHOUT renaming any folder; the other group stays
    where it is, still a record, and every reader then gives the chosen one;
  * a folder that is not a search of the session is refused, with the ones that are;
  * the finalize MANIFEST keeps the per-group [NOTE] lines (faims_step) and the how-to-choose
    sentence leads with choose-search.
Hermetic: temp sessions.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import faims  # noqa: E402
import session  # noqa: E402
from test_no_de_brief import Base, write  # noqa: E402


def split_session(s, out="search"):
    """Two group searches under output/<out>/, as run_search --faims-split writes them."""
    others = {"no FAIMS": 58, "FAIMS": 23}
    for key, name, n in (("none", "no FAIMS", 58), ("faims", "FAIMS", 23)):
        rec = faims.split_record(name, [] if key == "none" else [-55.0], n,
                                 {k: v for k, v in others.items() if k != name})
        write(os.path.join(s, "output", out, f"faims_{key}", "search_provenance.json"),
              {"engine": "fragpipe", "version": "24.0", "faims": rec})
    write(os.path.join(s, "output", out, "faims_split.json"), {"comparable_across_groups": False})


class TheChoice(Base):
    def p(self):
        return session.paths_for(self.s)

    def test_a_split_is_said_to_be_a_split_not_none_found_or_superseded(self):
        split_session(self.s)
        c = session.search_choice(self.p())
        self.assertEqual(c["how"], "faims_split")
        self.assertTrue(session.search_unchosen(c))
        self.assertIsNone(c["path"])
        for word in ("faims_none", "faims_faims", "no FAIMS", "58 run", "23 run",
                     "not comparable across the groups"):
            self.assertIn(word, c["why"])
        for bad in ("No search output", "superseded"):
            self.assertNotIn(bad, c["why"] + (c["short"] or ""))
        self.assertTrue(session.SEARCH_HOW_TO_CHOOSE.startswith("to choose one, run session.py "
                                                                "choose-search"))

    def test_choose_search_names_a_group_without_renaming_anything(self):
        split_session(self.s)
        out = os.path.join(self.s, "output")
        before = sorted(os.listdir(os.path.join(out, "search")))
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "choose-search",
                            "--dir", self.s, "--folder", "search/faims_faims"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(sorted(os.listdir(os.path.join(out, "search"))), before)
        c = session.search_choice(self.p())
        self.assertEqual(c["how"], "chosen")
        self.assertEqual(os.path.basename(os.path.dirname(c["path"])), "faims_faims")
        self.assertFalse(session.search_unchosen(c))
        self.assertEqual(session.search_folder(self.p(), c), os.path.dirname(c["path"]))
        self.assertEqual(len(c["candidates"]), 2)         # the other group is still a record
        self.assertEqual(session.search_engine(self.p())["engine"], "fragpipe")

    def test_a_folder_that_is_no_search_is_refused_with_the_ones_that_are(self):
        split_session(self.s)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "choose-search",
                            "--dir", self.s, "--folder", "nowhere"],
                           capture_output=True, text=True, timeout=60)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("search/faims_none", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.s, "output", "search_choice.json")))

    def choose(self, folder, flag="--dir"):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"),
                               "choose-search", flag, self.s, "--folder", folder],
                              capture_output=True, text=True, timeout=60)

    def de_read(self, folder):
        rep = os.path.join(self.s, "output", folder, "report.parquet")
        write(rep, "x")
        write(os.path.join(self.p()["de_dir"], "de_provenance.json"), {"input": rep})

    def test_a_choice_the_des_input_overrides_is_refused_not_printed_chosen(self):
        """C2 review ux M1: the DE read faims_none; choosing faims_faims printed "chosen" while
        every reader still described faims_none."""
        split_session(self.s)
        self.de_read("search/faims_none")
        p = self.choose("search/faims_faims")
        self.assertNotEqual(p.returncode, 0, p.stdout)
        self.assertNotIn("chosen:", p.stdout)
        self.assertIn("NOT chosen", p.stderr)
        self.assertIn("output/search/faims_none", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.s, "output", "search_choice.json")))
        c = session.search_choice(self.p())
        self.assertEqual((c["how"], os.path.basename(os.path.dirname(c["path"]))),
                         ("de_input", "faims_none"))

    def test_choosing_the_search_the_de_read_is_fine_and_says_so(self):
        split_session(self.s)
        self.de_read("search/faims_none")
        p = self.choose("search/faims_none", flag="--session")        # --session is an alias
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("chosen: output/search/faims_none", p.stdout)
        self.assertIn("the search the DE read", p.stdout)

    def test_a_de_on_a_report_no_record_describes_refuses_any_choice(self):
        split_session(self.s)
        rep = os.path.join(self.s, "output", "loose", "report.parquet")
        write(rep, "x")
        write(os.path.join(self.p()["de_dir"], "de_provenance.json"), {"input": rep})
        p = self.choose("search/faims_faims")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("no search record describes it", p.stderr)

    def test_a_split_straight_under_output_is_found_too(self):
        split_session(self.s, out=".")
        self.assertEqual(session.search_choice(self.p())["how"], "faims_split")

    def test_two_unrelated_searches_are_still_none_chosen(self):
        for f, eng in (("search_a", "sage"), ("search_b", "diann")):
            write(os.path.join(self.s, "output", f, "search_provenance.json"),
                  {"engine": eng, "version": "1"})
        self.assertEqual(session.search_choice(self.p())["how"], "ambiguous")

    def test_the_manifest_notes_name_each_group(self):
        split_session(self.s)

        class Man:
            notes = []

            def note(self, name, text):
                self.notes.append((name, text))
        m = Man()
        session.faims_step(self.p(), m)
        self.assertEqual(len(m.notes), 2)
        self.assertTrue(all("not comparable across the groups" in t for _, t in m.notes))


class TheBrief(Base):
    def test_the_brief_of_a_split_session_does_not_say_no_search_output(self):
        split_session(self.s)
        b = self.brief("--acquisition", "DDA")
        self.assertNotIn("No search output was found", b)
        self.assertNotIn("superseded", b)
        self.assertIn("split by FAIMS group", b)
        self.assertIn("choose-search", b)


if __name__ == "__main__":
    unittest.main()
