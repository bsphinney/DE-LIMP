#!/usr/bin/env python3
"""T44 follow-up (2.14.2, Ruling 80): the working copy of a Core session, the manifest's pin in the
report, and how to choose between searches.

N4: a Core submission runs on HIVE and the report is written in a working copy pulled to the
laptop (SKILL.md "Pull what the report needs", references/core-submissions.md). The pull took
output/search/report.parquet only -- no search record, no workflow manifest -- so the copy's
report dropped the engine line and the scan-window / mass-accuracy caution, silently, and for a
session searched into output/search_sage/ it took nothing. What these tests pin:
  * `session.py which-search` prints the session's search (session.search_choice): the folder,
    how it was chosen or why none was, and the files the copy needs (--paths: the record, the
    report the DE read, the workflow manifest); none chosen -> nothing to pull, why, exit 3;
  * the documented pull takes that list, not output/search/report.parquet;
  * a report whose session holds no search record says "not recorded", never omits the line.
Stdlib only, no network.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import two_search_fixture as fx                 # noqa: E402

SKILL = os.path.dirname(HERE)


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


class Tmp(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def which(self, p, *extra):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "which-search",
                            "--dir", p["session_dir"], *extra],
                           capture_output=True, text=True, timeout=120)
        return r.returncode, r.stdout, r.stderr

    def report_facts(self, p):
        import make_analysis_html as mah
        try:
            de = json.loads(read(os.path.join(p["de_dir"], "de_provenance.json")))
        except OSError:
            de = {}
        a = argparse.Namespace(session=p["session_dir"], tables=p["de_dir"],
                               out=os.path.join(p["output_dir"], "Analysis_Report.html"))
        header = dict(mah.session_facts(a, de)).get("Search")
        made = " ".join(b[1] for b in mah.how_made(a, de, None, []) if b[0] == "md")
        return header, made


# ------------------------------------------------------------------------------------- N4 --
class WhichSearch(Tmp):
    """`session.py which-search`: the one answer, printed on HIVE for the working copy's pull."""

    def test_the_real_case_names_the_search_the_de_read_and_what_to_pull(self):
        p, d = fx.two_searches(self.tmp)
        rc, out, err = self.which(p)
        self.assertEqual(rc, 0, err)
        js = json.loads(out)
        self.assertEqual((js["how"], js["folder"], js["record"]),
                         ("de_input", "output/search_sage",
                          "output/search_sage/search_provenance.json"))
        rc, out, err = self.which(p, "--paths")
        self.assertEqual(rc, 0, err)
        self.assertEqual(out.split(), ["output/search_sage/search_provenance.json",
                                       "output/search_sage/report.parquet",
                                       "input/wf/workflow.manifest.json"])

    def test_fragpipe_dia_pulls_the_report_the_de_read(self):
        p, _d = fx.layout(self.tmp, {
            "search_fragpipe": (fx.FRAGPIPE, fx.FRAGPIPE_DIA)},
            de_read=("search_fragpipe", fx.FRAGPIPE_DIA))
        rc, out, err = self.which(p, "--paths")
        self.assertEqual(rc, 0, err)
        self.assertIn("output/search_fragpipe/dia-quant-output/report.tsv", out.split())

    def test_none_chosen_pulls_nothing_and_says_why(self):
        p, _d = fx.layout(self.tmp, {
            "search": (fx.DIANN, "report.parquet"),
            "search_sage": (fx.SAGE, "report.parquet")})
        rc, out, err = self.which(p, "--paths")
        self.assertEqual((rc, out), (3, ""))
        self.assertIn("2 search records", err)
        rc, out, _err = self.which(p)
        self.assertEqual(rc, 3)
        js = json.loads(out)
        self.assertEqual((js["how"], js["folder"]), ("ambiguous", None))

    def test_pin_output_search(self):
        p, _d = fx.one_search(self.tmp)
        rc, out, err = self.which(p, "--paths")
        self.assertEqual(rc, 0, err)
        self.assertEqual(out.split(), ["output/search/search_provenance.json",
                                       "output/search/report.parquet",
                                       "input/wf/workflow.manifest.json"])


class ThePullTakesTheList(Tmp):
    """The documented pull asks which-search, and no longer takes output/search/report.parquet."""

    def test_the_pull_blocks(self):
        for rel in ("SKILL.md", os.path.join("references", "core-submissions.md")):
            text = read(os.path.join(SKILL, rel))
            self.assertNotIn('--get "$S/output/search/report.parquet"', text, rel)
            # $S quoted for the remote shell: a lab folder can hold a space (verify V1)
            self.assertIn('session.py which-search --dir $(printf %q "$S") --paths', text, rel)


class ACopyWithNoRecordSaysSo(Tmp):
    """A working copy that holds the report and no search record: the report says "not
    recorded", never omits the engine line."""

    def test_header_and_how_made_say_not_recorded(self):
        p, d = fx.one_search(self.tmp)
        os.remove(p["search_prov"])                    # the old pull: report.parquet only
        os.remove(p["workflow_manifest"])
        header, made = self.report_facts(p)
        self.assertEqual(header, "not recorded (no search record)")
        self.assertIn("Which search these numbers come from was not recorded (no search "
                      "record).", made)



# ------------------------------------------------------------------------------------- N2 --
class ThePinIsTagged(Tmp):
    """N2 (rule 2): with no search record, the report's engine is the workflow manifest's pin --
    a request, not a record of what ran -- and says so, as the README, Slack and the Methods do."""

    def test_no_record_the_header_and_how_made_tag_the_pin(self):
        p, _d = fx.one_search(self.tmp)
        os.remove(p["search_prov"])
        header, made = self.report_facts(p)
        self.assertEqual(header, "DIA-NN 2.7.0 (workflow manifest pin)")
        self.assertIn("DIA-NN 2.7.0 (the version planned for this workflow; the search's own "
                      "record was not available) searched the raw instrument files", made)

    def test_pin_a_record_is_not_tagged(self):
        p, _d = fx.one_search(self.tmp)
        header, made = self.report_facts(p)
        self.assertEqual(header, "DIA-NN 2.7.0")
        self.assertIn("DIA-NN 2.7.0 searched the raw instrument files", made)



# ------------------------------------------------------------------------------------- N1 --
class HowToChoose(Tmp):
    """N1: several searches and no DE -- the refusal and the docs say how to choose one."""

    def layout_E(self):
        return fx.layout(self.tmp, {
            "search": (fx.DIANN, "report.parquet"),
            "search_sage": (fx.SAGE, "report.parquet")})

    def test_deliver_says_how_to_choose(self):
        import core_submission
        p, _d = self.layout_E()
        why = core_submission.delivery_search(p["output_dir"])[1]
        self.assertIn("superseded", why)
        self.assertIn("--mode interim", why)
        self.assertIn("output/search_sage/search_provenance.json", why)

    def test_which_search_and_the_methods_note_say_how_to_choose(self):
        import make_deposit as md
        p, _d = self.layout_E()
        _rc, _out, err = self.which(p, "--paths")
        self.assertIn("superseded", err)
        f = md.gather(p["session_dir"])
        man = md.Manifest()
        md.search_ambiguity(f, man)
        note = " ".join(man.for_agent.get("methods_to_fix") or [])
        self.assertIn("superseded", note)

    def test_renaming_the_other_folder_chooses_the_one_left(self):
        import session
        p, d = self.layout_E()
        os.rename(d["search"], os.path.join(p["output_dir"], "superseded_search"))
        c = session.search_choice(p)
        self.assertEqual((c["how"], os.path.basename(os.path.dirname(c["path"]))),
                         ("only", "search_sage"))

    def test_the_docs_say_how_to_choose(self):
        for rel in ("SKILL.md", os.path.join("references", "core-submissions.md")):
            text = read(os.path.join(SKILL, rel))
            self.assertIn("superseded_<name>", text, rel)


if __name__ == "__main__":
    unittest.main()
