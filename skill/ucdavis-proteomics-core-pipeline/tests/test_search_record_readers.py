#!/usr/bin/env python3
"""Every reader of a session's search record finds it wherever the search is (2.14.2, T44).

A 2.14.1 session searched two timsTOF DDA runs into output/search_dda/ (DIA-NN --dda, whose
second pass failed) and output/search_sage/ (Sage), never into output/search/, and its DE read
the Sage report. Readers that looked only in output/search/ missed both searches, and the
session left the skill's path. session.find_search_provenance is the one definition of where a
session's search records are, tie-broken by the report the DE read; session.search_record is
its wrapper for a session. What these tests pin, for the real case (two_search_fixture):
  * session.py: the resolved parameters (params_file), a re-analysis's link to the prior's search
    record, and the raw list taken from the search record are the searches the DE read;
  * the report's engine line and its search record (make_analysis_html);
  * the Slack posts: a search's session, and the analysis post's engine (notify_slack);
  * the analysis delivery's search files and the README line about them (core_submission);
  * the run log's search folder (record_run.session_search_out);
  * finalize's Methods describe the search the DE read (make_deposit.gather);
  * the peak check's default search folder (skyline_peaks.resolve_inputs);
and (pin) a session whose one search is output/search/ reads exactly as before.
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

import session                                  # noqa: E402
import session_docs                             # noqa: E402
import two_search_fixture as fx                 # noqa: E402
from job_env import job_env                     # noqa: E402


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


class Tmp(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)


# ------------------------------------------------------------------------------ session.py --
class SessionReaders(Tmp):
    """session.params_file (resolved), a re-analysis's prior_records, session_docs.raw_record."""

    def init_reanalysis(self, prior):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init",
                            "--name", "again", "--date", "2026-11-03", "--reanalysis-of", prior,
                            "--beside"], capture_output=True, text=True, env=job_env(self.tmp),
                           timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        new = json.loads(r.stdout)["created"]
        return json.loads(read(os.path.join(new, session.REANALYSIS_JSON)))

    def test_the_resolved_parameters_are_the_search_the_de_read(self):
        p, d = fx.two_searches(self.tmp)
        self.assertEqual(session.params_file(p), os.path.join(d["sage"], "params.resolved.cfg"))

    def test_a_reanalysis_links_the_priors_search_record_wherever_it_is(self):
        p, d = fx.two_searches(self.tmp)
        link = self.init_reanalysis(p["session_dir"])
        self.assertEqual(link["prior_records"].get("search_prov"),
                         os.path.join(d["sage"], "search_provenance.json"))

    def test_the_raw_list_comes_from_the_search_the_de_read(self):
        p, d = fx.two_searches(self.tmp, raw_list=False)
        raws, src = session_docs.raw_record(p)
        self.assertEqual(src, "output/search_sage/search_provenance.json")
        self.assertEqual([os.path.basename(r) for r in raws], ["runA_01.d", "runB_01.d"])

    def test_pin_output_search_reads_as_before(self):
        p, d = fx.one_search(self.tmp, raw_list=False)
        self.assertEqual(session.params_file(p),
                         os.path.join(p["search_out"], "params.resolved.cfg"))
        self.assertEqual(session_docs.raw_record(p)[1], "output/search/search_provenance.json")
        link = self.init_reanalysis(p["session_dir"])
        self.assertEqual(link["prior_records"]["search_prov"], p["search_prov"])


# ---------------------------------------------------------------------------- notify_slack --
class SlackPosts(Tmp):
    """A search's session (search-done) and the analysis post's engine (analysis-done)."""

    def setUp(self):
        super().setUp()
        import notify_slack
        self.ns = notify_slack

    def test_a_search_in_output_search_sage_names_its_session(self):
        p, d = fx.two_searches(self.tmp)
        for folder in (d["sage"], d["dda"]):
            f = self.ns.search_facts(folder, exit_code=0)
            self.assertEqual(f["session"], os.path.basename(p["session_dir"]), folder)
            self.assertIsNone(f["folder"])

    def test_the_analysis_post_names_the_engine_of_the_search_the_de_read(self):
        p, _d = fx.two_searches(self.tmp)
        f = self.ns.analysis_facts(p["session_dir"])
        self.assertEqual((f["engine"], f["engine_version"]), ("Sage", "0.14.7"))

    def test_a_folder_outside_a_session_names_none(self):
        out = fx.write(os.path.join(self.tmp, "job", "output", "search_x",
                                    "search_provenance.json"), json.dumps(fx.SAGE))
        self.assertIsNone(self.ns.search_facts(os.path.dirname(out), exit_code=0)["session"])

    def test_pin_output_search_reads_as_before(self):
        p, d = fx.one_search(self.tmp)
        f = self.ns.search_facts(d["search"], exit_code=0)
        self.assertEqual((f["session"], f["folder"]), (os.path.basename(p["session_dir"]), None))
        a = self.ns.analysis_facts(p["session_dir"])
        self.assertEqual((a["engine"], a["engine_version"]), ("DIA-NN", "2.7.0"))
        # <x>/output/search names <x> as it always has, a session or not
        bare = os.path.join(self.tmp, "bare", "output", "search")
        os.makedirs(bare)
        self.assertEqual(self.ns.search_facts(bare, exit_code=0)["session"], "bare")


# ---------------------------------------------------------------------- make_analysis_html --
class ReportSearch(Tmp):
    """The header's Search fact, "How this report was made", and the report's search record."""

    def setUp(self):
        super().setUp()
        import make_analysis_html
        self.mah = make_analysis_html

    def page_args(self, p):
        return argparse.Namespace(session=p["session_dir"], tables=p["de_dir"],
                                  out=os.path.join(p["output_dir"], "Analysis_Report.html"))

    def de(self, p):
        return json.loads(read(os.path.join(p["de_dir"], "de_provenance.json")))

    def facts_and_how_made(self, p):
        a, prov = self.page_args(p), self.de(p)
        facts = dict(self.mah.session_facts(a, prov))
        made = " ".join(b[1] for b in self.mah.how_made(a, prov, None, []) if b[0] == "md")
        return facts, made

    def test_the_engine_line_names_the_engine_of_the_search_the_de_read(self):
        p, _d = fx.two_searches(self.tmp)
        facts, made = self.facts_and_how_made(p)
        self.assertEqual(facts["Search"], "Sage 0.14.7")
        self.assertIn("Sage 0.14.7 searched the raw instrument files", made)
        self.assertNotIn("DIA-NN", made)

    def test_the_search_record_is_found_outside_output_search(self):
        # an older DE record names no input: the session's one search record, in search_dda/
        p, d = fx.two_searches(self.tmp, de=False)
        shutil.rmtree(d["sage"])
        rec = self.mah.search_record({}, p["session_dir"])
        self.assertEqual((rec.get("engine"), rec.get("search_mode")), ("diann", "dda"))

    def test_pin_output_search_reads_as_before(self):
        p, _d = fx.one_search(self.tmp)
        facts, made = self.facts_and_how_made(p)
        self.assertEqual(facts["Search"], "DIA-NN 2.7.0")
        self.assertIn("DIA-NN 2.7.0 searched the raw instrument files", made)
        self.assertEqual(self.mah.search_record({}, p["session_dir"]).get("engine"), "diann")


# ------------------------------------------------------------------------------ record_run --
class RunLogSearch(Tmp):
    """The search folder the run log records and copies from (session_search_out)."""

    def setUp(self):
        super().setUp()
        import record_run
        self.rr = record_run

    def test_the_run_log_copies_from_the_search_the_de_read(self):
        p, d = fx.two_searches(self.tmp)
        self.assertEqual(self.rr.session_search_out(p["session_dir"], None), d["sage"])

    def test_a_lone_search_outside_output_search_is_the_one(self):
        p, d = fx.two_searches(self.tmp, de=False)
        shutil.rmtree(d["sage"])
        self.assertEqual(self.rr.session_search_out(p["session_dir"], None), d["dda"])

    def test_pin_output_search_reads_as_before(self):
        p, d = fx.one_search(self.tmp)
        self.assertEqual(self.rr.session_search_out(p["session_dir"], None), d["search"])
        self.assertEqual(self.rr.session_search_out(p["session_dir"], "/x/out"), "/x/out")
        # output/search/ with a report and no record (an older search) is still the one
        os.remove(p["search_prov"])
        self.assertEqual(self.rr.session_search_out(p["session_dir"], None), d["search"])


# --------------------------------------------------------------------------- skyline_peaks --
class PeakCheckSearch(Tmp):
    """The peak check's search folder with --session and no --search (resolve_inputs)."""

    def setUp(self):
        super().setUp()
        import skyline_peaks
        self.sp = skyline_peaks

    def args(self, sd, search=None):
        return argparse.Namespace(session=sd, search=search, out=None, de_dir=None,
                                  manifest=None, conditions=None, name=None)

    def test_the_peak_check_reads_the_search_the_de_read(self):
        p, d = fx.two_searches(self.tmp)
        a = self.args(p["session_dir"])
        inp = self.sp.resolve_inputs(a)
        self.assertEqual(inp["search"], d["sage"])
        # and its reason not to run is the true one (DIA-NN only), not "no search folder"
        with self.assertRaises(self.sp.Skip) as cm:
            self.sp.make_plan(a, inp)
        # fixer A (C2 bio H2) names the engine as recorded ("Sage"), not lower-cased
        self.assertIn("the search was run with sage", str(cm.exception).lower())

    def test_pin_output_search_reads_as_before(self):
        p, _d = fx.one_search(self.tmp)
        self.assertEqual(self.sp.resolve_inputs(self.args(p["session_dir"]))["search"],
                         p["search_out"])
        given = os.path.join(self.tmp, "elsewhere")
        self.assertEqual(self.sp.resolve_inputs(self.args(p["session_dir"], given))["search"],
                         given)


# ---------------------------------------------------------------------------- fran_deposit --
class FranSession(Tmp):
    """fran_deposit.session_for already reads any <session>/output/<folder> (its docstring said
    output/search): kept so."""

    def test_a_search_in_output_search_sage_has_its_session(self):
        import fran_deposit
        p, d = fx.two_searches(self.tmp)
        for folder in (d["sage"], d["dda"]):
            self.assertEqual(fran_deposit.session_for(folder), p["session_dir"])


# ---------------------------------------------------------------------------- make_deposit --
class MethodsSearch(Tmp):
    """finalize's Methods (make_deposit.gather) describe the search the DE read (Ruling RT-F7):
    the same choice as its bundle command (_chosen_search) and the README."""

    def setUp(self):
        super().setUp()
        import make_deposit
        self.md = make_deposit

    def described(self, p):
        f = self.md.gather(p["session_dir"])
        man = self.md.Manifest()
        self.md.search_ambiguity(f, man)
        cmd = self.md.methods_command(f, os.path.join(self.tmp, "methods.md"))
        provs = [cmd[i + 1] for i, x in enumerate(cmd) if x == "--search-prov"]
        self.by_session = "--session" in cmd
        return f, [ln for ln in man.lines if "Search paragraph" in ln], provs

    def test_finalize_describes_the_search_the_de_read(self):
        p, d = fx.two_searches(self.tmp)
        f, lines, provs = self.described(p)
        sage = os.path.join(d["sage"], "search_provenance.json")
        self.assertEqual(f["search_prov_path"], sage)
        self.assertEqual(f["srec"].get("engine"), "sage")
        self.assertFalse(f["srec"].get("ambiguous"))
        self.assertEqual((lines, provs), ([], [sage]))

    def test_with_no_tie_break_neither_is_described_as_before(self):
        p, _d = fx.two_searches(self.tmp, de=False)
        f, lines, provs = self.described(p)
        self.assertIsNone(f["search_prov_path"])
        self.assertEqual(f["srec"].get("ambiguous"), 2)
        # make_methods.py makes the same choice from --session (Ruling 78 F3): no record named
        self.assertEqual((len(lines), len(provs), self.by_session), (1, 0, True))

    def test_pin_output_search_reads_as_before(self):
        p, _d = fx.one_search(self.tmp)
        f, lines, provs = self.described(p)
        self.assertEqual((f["search_prov_path"], lines, provs),
                         (p["search_prov"], [], [p["search_prov"]]))


# ------------------------------------------------------------------------- core_submission --
import test_core_submission as tcs              # noqa: E402  (the deliver fixture, the fake API)


class DeliverSearch(tcs.DeliverBase):
    """An analysis delivery carries the search the DE read, under its own folder name (the
    delivery drops `output/`: AGENTS_DELIVERY_NOTE), and its README says where."""

    def setUp(self):
        super().setUp()
        import finish_fixture
        shutil.rmtree(os.path.join(self.output, "search"))
        self.sage = os.path.join(self.output, "search_sage")
        dda = os.path.join(self.output, "search_dda")
        for folder, rec, report in ((self.sage, fx.SAGE, "report.parquet"),
                                    (dda, fx.DIANN, "report-first-pass.parquet")):
            fx.write(os.path.join(folder, report), "parquet")
            fx.write(os.path.join(folder, "search_provenance.json"), json.dumps(rec))
        self.de_record(os.path.join(self.sage, "report.parquet"))
        self.restamp = lambda: finish_fixture.stamp(self.output)
        self.restamp()

    def de_record(self, report):
        fx.write(os.path.join(self.output, "tables", "de_provenance.json"), json.dumps(
            dict(tcs.NORM_DECIDED, **({"input": report} if report else {}))))

    def items(self):
        rc, out, p = self.deliver()
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        return {i["name"]: i for i in out["items"]}

    def test_deliver_plans_the_search_the_de_read(self):
        items = self.items()
        self.assertEqual(items["search_sage/report.parquet"]["status"], "OK")
        self.assertEqual(items["search_sage/search_provenance.json"]["status"], "OK")
        self.assertFalse([n for n in items if n.startswith(("search/", "search_dda/"))], items)

    def test_the_readme_names_the_delivered_search_folder(self):
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertTrue(os.path.isfile(os.path.join(self.delivery, "search_sage", "report.parquet")))
        readme = read(os.path.join(self.delivery, "README.md"))
        self.assertIn("the search output is in `search_sage/`", readme)
        self.assertIn("| `search_sage/` | the search engine's own output", readme)
        self.assertNotIn("`search/`", readme)

    def test_two_searches_and_no_tie_break_are_said_not_guessed(self):
        self.de_record(None)                          # a DE record that names no input
        self.restamp()
        items = self.items()
        self.assertFalse([n for n in items if n.startswith(("search/", "search_sage/",
                                                            "search_dda/")) and n != "search/"])
        why = items["search/"]
        self.assertEqual(why["status"], "SKIPPED")
        self.assertIn("output/search_dda/search_provenance.json", why["reason"])
        self.assertIn("output/search_sage/search_provenance.json", why["reason"])


class DeliverSearchPin(tcs.DeliverBase):
    """(pin) output/search/ -- with no record, as the deliver fixture has it, or with one -- is
    delivered as search/, as before."""

    def test_pin_output_search_is_delivered_as_before(self):
        items = {i["name"]: i for i in self.deliver()[1]["items"]}
        self.assertEqual(items["search/report.parquet"]["status"], "OK")
        self.assertEqual(items["search/search_provenance.json"]["status"], "SKIPPED")
        fx.write(os.path.join(self.output, "search", "search_provenance.json"),
                 json.dumps(fx.DIANN))
        items = {i["name"]: i for i in self.deliver()[1]["items"]}
        self.assertEqual(items["search/search_provenance.json"]["status"], "OK")
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertIn("the search output is in `search/`", read(os.path.join(self.delivery,
                                                                             "README.md")))


if __name__ == "__main__":
    unittest.main()
