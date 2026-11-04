#!/usr/bin/env python3
"""finalize's Methods find the search wherever it is, and never claim "optimised automatically"
without a record of the run (2.14.2 triage I04).

A Core member's lateElu re-analysis (2026-10-20) searched into output/search_dda/. finalize looked only in
output/search/ and output/search/*/, found no search record, read the unresolved params file the
search was given -- which does not carry the mass accuracy the search pinned at 10/20 ppm -- and
wrote "mass accuracy was not fixed in the parameters; DIA-NN optimised it automatically per run".

What these tests pin:
  * make_deposit.find_search_provenance() looks in output/*/ and output/*/*/, by real path, and
    skips any folder whose name says "superseded"; one record is THE record, several are named
    and none is chosen
  * gather() uses it; with several, finalize's Methods carry a tagged blank search paragraph and
    MANIFEST.txt a [SKIPPED] line naming every record
  * search_record() writes the "optimised automatically" sentence only with a run record
  * (pin) a session whose only search is output/search/ gives the Methods it gave in 2.14.1
stdlib only, no network.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_deposit as md                                        # noqa: E402
import make_methods as mm                                        # noqa: E402
from test_deposit_package import dia_session, write              # noqa: E402

CFG_NO_MASS_ACC = "--qvalue 0.01\n--cut K*,R*\n--missed-cleavages 1\n--unimod4\n"
CFG_PINNED = CFG_NO_MASS_ACC + "--mass-acc 20\n--mass-acc-ms1 10\n"


def out_with(root, *rels):
    """An output/ folder holding a search_provenance.json at each relative path."""
    out = os.path.join(root, "output")
    for rel in rels:
        write(os.path.join(out, rel), json.dumps({"engine": "diann", "version": "2.7.0"}))
    os.makedirs(out, exist_ok=True)
    return out


class FindSearchProvenance(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_search_in_output_search_dda_is_found(self):                # not output/search
        path, cands = md.find_search_provenance(
            out_with(self.tmp, "search_dda/search_provenance.json"))
        self.assertTrue(path.endswith(os.path.join("search_dda", "search_provenance.json")))
        self.assertEqual(len(cands), 1)

    def test_one_level_deeper_is_found(self):
        path, _ = md.find_search_provenance(
            out_with(self.tmp, "search/rerun/search_provenance.json"))
        self.assertTrue(path.endswith(os.path.join("rerun", "search_provenance.json")))

    def test_two_searches_choose_neither_and_say_so(self):
        path, cands = md.find_search_provenance(out_with(
            self.tmp, "search_a/search_provenance.json", "search_b/search_provenance.json"))
        self.assertIsNone(path)
        self.assertEqual(len(cands), 2)

    def test_a_superseded_folder_is_not_a_candidate(self):
        path, cands = md.find_search_provenance(out_with(
            self.tmp, "search/search_provenance.json",
            "search_superseded_2026-10-19/search_provenance.json",
            "search/superseded/search_provenance.json"))
        self.assertTrue(path.endswith(os.path.join("output", "search", "search_provenance.json")))
        self.assertEqual(len(cands), 1)

    def test_the_same_file_twice_is_one_candidate(self):
        out = out_with(self.tmp, "search/search_provenance.json")
        os.symlink(os.path.join(out, "search"), os.path.join(out, "search_link"))
        path, cands = md.find_search_provenance(out)
        self.assertIsNotNone(path)
        self.assertEqual(len(cands), 1)

    def test_none(self):
        self.assertEqual(md.find_search_provenance(out_with(self.tmp)), (None, []))

    # 2.14.2 review (ruling 35): the finder lives in session.py, beside paths_for, so the
    # readers that import session (session_docs, make_analysis_html, notify_slack) need not pull
    # in make_deposit -> make_methods -> fetch_fasta; make_deposit re-exports it
    def test_it_lives_in_session_and_make_deposit_re_exports_it(self):
        import session
        self.assertIs(md.find_search_provenance, session.find_search_provenance)

    def test_the_record_beside_the_de_input_breaks_a_tie(self):
        import session
        out = out_with(self.tmp, "search_a/search_provenance.json",
                       "search_b/search_provenance.json")
        de_input = os.path.join(out, "search_b", "report.parquet")
        path, cands = session.find_search_provenance(out, de_input=de_input)
        self.assertTrue(path.endswith(os.path.join("search_b", "search_provenance.json")))
        self.assertEqual(len(cands), 2)
        # a DE input that is beside neither breaks no tie
        self.assertIsNone(session.find_search_provenance(
            out, de_input=os.path.join(self.tmp, "elsewhere", "report.parquet"))[0])


class ToleranceClaim(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.cfg = os.path.join(self.tmp, "params.cfg")
        write(self.cfg, CFG_NO_MASS_ACC)
        self.wf = os.path.join(self.tmp, "workflow.manifest.json")
        write(self.wf, json.dumps({"acquisition": "DIA",
                                   "engine": {"name": "diann", "version": "2.7.0"}}))

    def test_no_run_record_never_claims_automatic_mass_accuracy(self):
        rec = mm.search_record(self.cfg, None, self.wf)
        self.assertNotIn("optimised it automatically", rec.get("tol_note") or "")
        self.assertIn("not recorded", rec["tol_note"])
        self.assertIn(mm.NR_TAG, mm.search_paragraph(rec))

    def prov(self, **extra):
        path = os.path.join(self.tmp, "search", "search_provenance.json")
        write(path, json.dumps(dict({"engine": "diann", "version": "2.7.0",
                                     "params_file": self.cfg}, **extra)))
        return path

    def test_recorded_unpinned_search_still_says_automatic(self):        # (pin)
        rec = mm.search_record(None, self.prov(mass_acc={"mode": "auto", "fixed": False}),
                               self.wf)
        self.assertIn("optimised it automatically", rec["tol_note"])

    def test_the_resolved_cfg_that_ran_without_mass_acc_is_automatic(self):   # (pin)
        resolved = os.path.join(self.tmp, "run", "params.resolved.cfg")
        write(resolved, CFG_NO_MASS_ACC)
        rec = mm.search_record(None, self.prov(resolved_params_file=resolved), self.wf)
        self.assertIn("optimised it automatically", rec["tol_note"])

    # 2.14.2 review: the 5-step chain and the single-shot probe record result.mass_acc
    # {mode: measured, fixed: True, ms1/ms2: None} and a resolved cfg on HIVE; read away from
    # HIVE, the given params (no --mass-acc) wrote the automatic claim for a pinned search
    MEASURED = {"mode": "measured", "fixed": True, "measured": True, "ms1": None, "ms2": None}

    def test_a_measured_record_whose_value_is_not_readable_is_not_automatic(self):
        prov = self.prov(resolved_params_file="/quobyte/elsewhere/params.resolved.cfg",
                         result={"mass_acc": self.MEASURED})
        rec = mm.search_record(None, prov, self.wf)
        self.assertNotIn("optimised it automatically", rec["tol_note"])
        self.assertIn("measured before the search and pinned", rec["tol_note"])
        self.assertIn(mm.NR_TAG, rec["tol_note"])

    def test_the_top_level_record_is_read_too(self):
        rec = mm.search_record(None, self.prov(mass_acc=self.MEASURED), self.wf)
        self.assertIn("measured before the search and pinned", rec["tol_note"])

    def test_a_record_that_does_not_say_is_not_automatic(self):
        rec = mm.search_record(None, self.prov(), self.wf)       # the GIVEN params, no mass_acc
        self.assertNotIn("optimised it automatically", rec["tol_note"])
        self.assertIn("not recorded", rec["tol_note"])

    def test_the_resolved_cfg_beside_the_record_is_read(self):
        prov = self.prov(resolved_params_file="/quobyte/elsewhere/params.resolved.cfg",
                         result={"mass_acc": self.MEASURED})
        write(os.path.join(os.path.dirname(prov), "params.resolved.cfg"), CFG_PINNED)
        rec = mm.search_record(None, prov, self.wf)
        self.assertEqual((rec["ms1_tol"]["value"], rec["ms2_tol"]["value"]), (10.0, 20.0))
        self.assertIn("beside search_provenance.json", rec["params_source"])


class FinalizeMethods(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.p = dia_session(self.tmp)

    def methods(self):
        man = md.Manifest()
        path = md.ensure_methods(self.p["session_dir"], man)
        with open(path, encoding="utf-8") as fh:
            return md.md_sections(fh.read()).get("Database search", ""), man.lines

    def test_only_output_search_gives_the_same_methods(self):           # (pin)
        self.assertEqual(md.gather(self.p["session_dir"])["search_prov_path"],
                         self.p["search_prov"])
        para, lines = self.methods()
        self.assertIn("DIA-NN 2.7.0", para)
        self.assertNotIn("pinned version", para)
        self.assertIn("precursor (MS1) 15 ppm and fragment (MS2) 15 ppm", para)
        self.assertFalse([ln for ln in lines if "Search paragraph" in ln])

    def test_a_search_in_output_search_dda_is_the_one_described(self):
        # the session's own params file (what the search was given) does not fix the mass
        # accuracy; the search's resolved parameters pinned it at 10/20 ppm
        write(os.path.join(self.p["workflow_dir"], "params.cfg"), CFG_NO_MASS_ACC)
        dda = os.path.join(self.p["output_dir"], "search_dda")
        resolved = os.path.join(dda, "params.resolved.cfg")
        write(resolved, CFG_PINNED)
        write(os.path.join(dda, "search_provenance.json"), json.dumps({
            "engine": "diann", "version": "2.7.0", "resolved_params_file": resolved}))
        os.remove(self.p["search_prov"])
        para, _ = self.methods()
        self.assertIn("precursor (MS1) 10 ppm and fragment (MS2) 20 ppm", para)
        self.assertNotIn("optimised it automatically", para)
        self.assertNotIn("pinned version", para)

    def test_two_searches_are_named_and_neither_described(self):
        other = os.path.join(self.p["output_dir"], "search_dda", "search_provenance.json")
        write(other, json.dumps({"engine": "sage", "version": "0.14.7"}))
        para, lines = self.methods()
        self.assertIn(mm.NR_TAG, para)
        self.assertNotIn("DIA-NN 2.7.0", para)
        self.assertNotIn("Sage", para)
        self.assertNotIn("optimised it automatically", para)
        (line,) = [ln for ln in lines if "Search paragraph" in ln]
        self.assertTrue(line.startswith("[SKIPPED]"), line)
        self.assertIn("2 search records", line)
        self.assertIn(os.path.join("output", "search", "search_provenance.json"), line)
        self.assertIn(os.path.join("output", "search_dda", "search_provenance.json"), line)
        self.assertIn("make_methods.py --search-prov", line)

    def test_the_upload_list_walks_the_search_described(self):
        # 2.14.2 review: gather() described output/search_dda/, but the PRIDE file list walked
        # only output/search/
        dda = os.path.join(self.p["output_dir"], "search_dda")
        os.rename(self.p["search_out"], dda)
        rows = md.plan_uploads(md.gather(self.p["session_dir"]))
        by = {r["upload_name"]: r for r in rows}
        self.assertIn("report.parquet", by)
        self.assertEqual(os.path.dirname(by["report.parquet"]["source_path"]), dda)

    def test_gather_with_two_searches_reads_neither(self):
        write(os.path.join(self.p["output_dir"], "search_dda", "search_provenance.json"),
              json.dumps({"engine": "sage", "version": "0.14.7"}))
        f = md.gather(self.p["session_dir"])
        self.assertIsNone(f["search_prov_path"])
        self.assertEqual(len(f["search_prov_candidates"]), 2)
        self.assertIsNone(f["srec"].get("ms1_tol"))       # the SDRF is not filled from a guess


if __name__ == "__main__":
    unittest.main()
