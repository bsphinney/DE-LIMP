#!/usr/bin/env python3
"""The engine choice is logged and replayed in its record's own words (2.15, F5).

* `log_decision.py --from-json` writes step 4's decision from the record engine_choice.choose()
  made (a workflow manifest's engine.selection, resolve_defaults.py's printed `selection`, or a
  search_provenance.json's engine_selection): the client sentence quoted as written, never
  re-worded; the staff record beside it is not read; typed --what/--why beside it are refused,
  and a document with no selection is refused.
* REPRODUCE.md's engine line is the engine's name, version and how it was chosen -- not the
  manifest's record printed raw (2.15's engine record became a dict, and the line printed
  "{'name': 'fragpipe', ...}"); a FragPipe search's REPRODUCE.md carries the MSBooster re-run
  sentence (Ruling 106; fragpipe_params.msbooster_rerun_sentence, the one sentence).
"""
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

import provenance  # noqa: E402

FALLBACK = ("FragPipe, the pipeline's default search engine for timsTOF DDA-PASEF data, was not "
            "available where this search ran, so the data were searched with Sage.")
SEL_FALLBACK = {"engine": "sage", "default": "fragpipe", "how": "fallback",
                "data_type": "timsTOF DDA-PASEF data", "sentence": FALLBACK, "ids_only": False,
                "staff_record": "workflow.manifest.json.staff.json",
                "decided_by": "engine_choice.py"}
SEL_DEFAULT = {"engine": "fragpipe", "default": "fragpipe", "how": "default",
               "data_type": "timsTOF DDA-PASEF data", "sentence": None, "ids_only": True,
               "staff_record": None, "decided_by": "engine_choice.py"}


class LogDecisionFromJson(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.session = os.path.join(self.d, "S")
        os.makedirs(self.session)
        # a staff record beside the manifest, which must never be read into the log
        with open(os.path.join(self.d, "workflow.manifest.json.staff.json"), "w") as fh:
            json.dump({"why": "SECRET-STAFF-DETAIL darwin/arm64 /home/someone"}, fh)

    def run_log(self, doc, *extra):
        src = os.path.join(self.d, "workflow.manifest.json")
        with open(src, "w") as fh:
            json.dump(doc, fh)
        env = dict(os.environ)
        env.pop("CLAUDE_CODE_SESSION_ID", None)
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "log_decision.py"),
                               "--session", self.session, "--from-json", src, *extra],
                              capture_output=True, text=True, env=env)

    def logged(self):
        with open(os.path.join(self.session, "logs", "decisions.md"), encoding="utf-8") as fh:
            return fh.read()

    def test_a_fallback_is_logged_with_its_sentence_as_written(self):
        p = self.run_log({"engine": {"name": "sage", "version": "0.14.7",
                                     "selection": SEL_FALLBACK}})
        self.assertEqual(p.returncode, 0, p.stderr)
        log = self.logged()
        self.assertIn("-- Search engine: Sage (fallback)", log)
        self.assertIn(f"- **Why:** {FALLBACK} (decided by engine_choice.py)", log)
        self.assertIn("- **Step:** 4. defaults (engine)", log)
        self.assertNotIn("SECRET-STAFF-DETAIL", log)

    def test_the_default_and_ids_only_from_each_record_shape(self):
        for doc in ({"selection": SEL_DEFAULT},                # resolve_defaults.py's summary
                    {"engine": "fragpipe", "engine_selection": SEL_DEFAULT}):   # search record
            shutil.rmtree(os.path.join(self.session, "logs"), True)
            p = self.run_log(doc)
            self.assertEqual(p.returncode, 0, p.stderr)
            log = self.logged()
            self.assertIn("-- Search engine: FragPipe (default); identifications only", log)
            self.assertIn("the default for timsTOF DDA-PASEF data is FragPipe; how: default", log)

    def test_typed_text_beside_it_and_a_record_without_a_choice_are_refused(self):
        p = self.run_log({"engine": {"name": "sage", "selection": SEL_FALLBACK}},
                         "--what", "something else")
        self.assertEqual(p.returncode, 2)
        self.assertIn("leave out --what/--why", p.stderr)
        p = self.run_log({"engine": {"name": "diann", "version": "2.7.0"}})
        self.assertEqual(p.returncode, 2)
        self.assertIn("holds no engine selection", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.session, "logs", "decisions.md")))

    def test_typed_decisions_still_need_both_parts(self):
        env = dict(os.environ)
        env.pop("CLAUDE_CODE_SESSION_ID", None)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "log_decision.py"),
                            "--session", self.session, "--what", "x"],
                           capture_output=True, text=True, env=env)
        self.assertEqual(p.returncode, 2)
        self.assertIn("one of the arguments --why --why-file is required", p.stderr)


class ReproduceMdEngine(unittest.TestCase):
    def test_the_engine_line_is_words_not_a_dict(self):
        line = provenance._engine_line({"engine": {"name": "sage", "version": "0.14.7",
                                                   "selection": SEL_FALLBACK}})
        self.assertEqual(line, f"sage 0.14.7 (chosen: fallback). {FALLBACK}")
        self.assertEqual(provenance._engine_line({"engine": {"name": "fragpipe", "version": "24.0",
                                                             "selection": SEL_DEFAULT}}),
                         "fragpipe 24.0 (chosen: default)")
        self.assertEqual(provenance._engine_line({"engine": "diann"}), "diann")
        self.assertNotIn("{", provenance._engine_line({"engine": {"name": "x"}}))

    def test_a_fragpipe_rerun_note_from_the_record(self):
        bw = {"runs": [{"run": "S1-a", "bandwidth": 0.1, "source": "chosen by MSBooster"}]}
        note = provenance.rerun_note({"engine": "fragpipe", "msbooster_rt_bandwidth": bw})
        self.assertIn("A re-run of this FragPipe search can differ slightly", note)
        self.assertIn("this search ran with 0.1 for one run", note)
        self.assertNotIn("S1-a", note)
        self.assertEqual(provenance.rerun_note({"engine": "sage"}), "")
        self.assertEqual(provenance.rerun_note({"engine": "fragpipe"}), "")


class AgentsAndCatalog(unittest.TestCase):
    """E2E reading (2.15 F5): AGENTS.md named input/wf/fragpipe.workflow (the default preset, the
    first *.workflow there) as "Search parameters as run" when another workflow had run, and the
    catalog called FragPipe's saved fragpipe.workflow "parameters actually used"."""

    def test_agents_md_names_the_workflow_fragpipe_ran(self):
        import session_docs
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        wfdir, search = os.path.join(d, "input", "wf"), os.path.join(d, "output", "search")
        os.makedirs(wfdir)
        os.makedirs(search)
        for n in (os.path.join(wfdir, "fragpipe.workflow"), os.path.join(wfdir, "custom.workflow"),
                  os.path.join(search, "fragpipe_input.workflow")):
            open(n, "w").close()
        p = {"workflow_dir": wfdir, "input_dir": os.path.join(d, "input"),
             "output_dir": os.path.join(d, "output"), "search_out": search}
        sp = {"fragpipe_workflow": {"given": "/elsewhere/custom.workflow",
                                    "run": "/elsewhere/search/fragpipe_input.workflow"}}
        self.assertEqual(session_docs._params_as_run(p, sp, search),
                         os.path.join(search, "fragpipe_input.workflow"))

    def test_the_catalog_tells_fragpipes_workflows_apart(self):
        import make_report
        self.assertIn("The workflow FragPipe ran",
                      make_report.describe("fragpipe_input.workflow")[1])
        self.assertIn("never replayed", make_report.describe("fragpipe.workflow")[1])
        self.assertIn("make_presets.py's record",
                      make_report.describe("x.workflow.provenance.json")[1])
        self.assertIn("actually used", make_report.describe("params.cfg")[1])


class ReportHowMade(unittest.TestCase):
    """E2E reading (2.15 F5): a no-DE FragPipe DDA report said only "FragPipe 24.0 searched the raw
    instrument files"; its protein amounts are directLFQ's, and the record says so."""

    def session_with(self, quantity):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        os.makedirs(os.path.join(d, "output", "search"))
        with open(os.path.join(d, "output", "search", "search_provenance.json"), "w") as fh:
            json.dump({"engine": "fragpipe", "version": "24.0", "fragpipe_quantity": quantity}, fh)
        return d

    def test_the_quantity_is_named_from_the_record(self):
        import make_analysis_html as mah
        self.assertEqual(mah.search_quantity_clause(self.session_with(
            {"quantity": "directLFQ", "status": "ok"})),
            "directLFQ estimated the protein amounts from FragPipe's ion intensities")
        self.assertEqual(mah.search_quantity_clause(self.session_with(
            {"quantity": "directLFQ", "status": "ok", "n_experiments": 5,
             "why_not_maxlfq": "directLFQ was requested instead of IonQuant's MaxLFQ (the "
                               "workflow sets ionquant.maxlfq=0)."})),
            "directLFQ estimated the protein amounts from FragPipe's ion intensities "
            "(requested instead of IonQuant's MaxLFQ)")
        self.assertIn("MaxLFQ", mah.search_quantity_clause(self.session_with(
            {"quantity": "IonQuant MaxLFQ Intensity", "status": "ok"})))
        # C2 bio M1: one experiment's plain Intensity carries the between-run caveat
        self.assertIn("not normalised between runs; compare runs with care",
                      mah.search_quantity_clause(self.session_with(
                          {"quantity": "IonQuant Intensity", "status": "ok"})))
        self.assertIsNone(mah.search_quantity_clause(self.session_with(
            {"quantity": "directLFQ", "status": "refused"})))
        self.assertIsNone(mah.search_quantity_clause(None))


if __name__ == "__main__":
    unittest.main()
