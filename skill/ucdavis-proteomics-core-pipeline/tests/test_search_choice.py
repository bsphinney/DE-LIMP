#!/usr/bin/env python3
"""One answer to "which search", everywhere (2.14.2, T44 fix round, Ruling 78).

T44 sent every reader of a session's search record through session.search_record, and the review
found layouts where they still disagreed, or agreed on a search the DE never read:
  J/K  the DE read FragPipe's DIA report (<workdir>/dia-quant-output/report.tsv): the tie-break
       looked only beside the report, so output/search/ (DIA-NN) was described (F1);
  H    the DE read a recordless output/search/ beside a recorded search: the lone record won (F2);
  E    output/search/ and a second search, no DE input: some readers took output/search/, the
       Methods and the bundle said none was chosen (F3, F10);
  I    a Sage record with no version (every conda-env Sage): Slack, the README and the Methods
       gave it the manifest's DIA-NN pin, 2.7.0 (F4);
  E    with none chosen, the report, Slack and the bundle named the manifest's engine (F5);
and delivery.json did not say which search folder went (F6), the peak check gave a false reason
(F7), SKILL.md's commands named output/search/ and make_methods took a missing --search-prov as
no record (F8), and the run log said nothing when it could not choose (F9).
What these tests pin: session.search_choice is the one answer -- search_folder_of(the DE's
input) ties, a DE input that is here and no record describes chooses none, and several with no
tie choose none -- and every reader gives it; one engine/version rule (session.search_engine);
and (pin) output/search/ alone, and the real case, read as before.
Stdlib only, no network.
"""
import argparse
import io
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import session                                  # noqa: E402
import two_search_fixture as fx                 # noqa: E402

SKILL = os.path.dirname(HERE)
# the fixture's raw folders are empty: make_methods.py reads the instrument and acquisition from
# the flags finalize passes from the session's records (make_deposit.methods_command)
FINALIZE_FLAGS = ("--instrument", "timsTOF HT", "--acquisition", "DDA")


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def tagged(rec, folder):
    return dict(rec, tag=folder)


class Tmp(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def de(self, p):
        try:
            return json.loads(read(os.path.join(p["de_dir"], "de_provenance.json")))
        except OSError:
            return {}

    def picks(self, p):
        """{reader: the folder (under output/) or engine label it took, None for none}."""
        import audit_results
        import core_submission
        import make_analysis_html as mah
        import make_deposit as md
        import notify_slack
        import record_run
        import session_docs
        import skyline_peaks
        out, sd, de = p["output_dir"], p["session_dir"], self.de(p)

        def folder(path):
            if not path:
                return None
            rel = os.path.relpath(os.path.dirname(os.path.abspath(path)), out)
            return None if rel.startswith("..") else rel.replace(os.sep, "/")

        f = md.gather(sd)
        cmd = md.bundle_command(f)[0]
        flag = lambda name: cmd[cmd.index(name) + 1] if name in cmd else None  # noqa: E731
        rows = md.plan_uploads(f)
        up = sorted({os.path.relpath(r["source_path"], out).split(os.sep)[0] for r in rows
                     if r["source_path"].startswith(out + os.sep)} - {"tables", "DATA_SUBMISSION"})
        inp = skyline_peaks.resolve_inputs(argparse.Namespace(
            session=sd, search=None, out=None, de_dir=None, manifest=None, conditions=None,
            name=None))
        peak = folder(os.path.join(inp["search"], "x")) if inp["search"] else None
        docs = session_docs.gather(sd)
        run_log = record_run.session_search_out(sd, None)
        return {
            "search_record": folder(session.search_record(p)[0]),
            "report_record": mah.search_record(de, sd).get("tag"),
            "report_engine": mah.search_engine(sd, de)[0],
            "slack_engine": notify_slack.analysis_facts(sd)["engine"],
            "docs_record": docs["sp"].get("tag"),
            "docs_engine": docs["engine"],
            "deliver": core_submission.delivery_search(out)[0],
            "methods": folder(f["search_prov_path"]),
            "methods_engine": f["srec"].get("engine"),
            "bundle_search": folder(md._chosen_search(f)[0]),
            "bundle_engine": flag("--engine"),
            "uploads": up[0] if len(up) == 1 else (up or None),
            "run_log": folder(os.path.join(run_log, "x")) if run_log else None,
            "peak_check": peak,
            "params": folder(session.params_file(p)),
            "audit": (audit_results.load_search_prov(de)[1] or {}).get("tag"),
        }


# ------------------------------------------------------------------------------------- F1 --
class FragPipeDiaTie(Tmp):
    """F1: the DE read <workdir>/dia-quant-output/report.tsv; its search is the workdir's."""

    def check_every_reader_takes_fragpipe(self, p):
        got = self.picks(p)
        want = {"search_record": "search_fragpipe", "report_record": "search_fragpipe",
                "report_engine": "FragPipe", "slack_engine": "FragPipe",
                "docs_record": "search_fragpipe", "docs_engine": "FragPipe",
                "deliver": "search_fragpipe", "methods": "search_fragpipe",
                "methods_engine": "fragpipe", "bundle_search": "search_fragpipe",
                "bundle_engine": "fragpipe", "uploads": "search_fragpipe",
                "run_log": "search_fragpipe", "params": "search_fragpipe",
                "audit": "search_fragpipe", "peak_check": "search_fragpipe"}
        self.assertEqual(got, want)

    def test_J_beside_output_search(self):
        p, _d = fx.layout(self.tmp, {
            "search": (tagged(fx.DIANN, "search"), "report.parquet"),
            "search_fragpipe": (tagged(fx.FRAGPIPE, "search_fragpipe"), fx.FRAGPIPE_DIA)},
            de_read=("search_fragpipe", fx.FRAGPIPE_DIA))
        self.check_every_reader_takes_fragpipe(p)

    def test_K_beside_another_search(self):
        p, _d = fx.layout(self.tmp, {
            "search_dda": (tagged(fx.DIANN, "search_dda"), "report-first-pass.parquet"),
            "search_fragpipe": (tagged(fx.FRAGPIPE, "search_fragpipe"), fx.FRAGPIPE_DIA)},
            de_read=("search_fragpipe", fx.FRAGPIPE_DIA))
        self.check_every_reader_takes_fragpipe(p)

    def test_search_folder_of(self):
        _p, d = fx.layout(self.tmp, {
            "search_fragpipe": (tagged(fx.FRAGPIPE, "search_fragpipe"), fx.FRAGPIPE_DIA),
            "search_x": (None, "report.parquet")})
        self.assertEqual(session.search_folder_of(os.path.join(d["search_fragpipe"],
                                                               fx.FRAGPIPE_DIA)),
                         d["search_fragpipe"])
        self.assertEqual(session.search_folder_of(os.path.join(d["search_x"], "report.parquet")),
                         d["search_x"])



NONE = {k: None for k in ("search_record", "report_record", "report_engine", "slack_engine",
                          "docs_record", "docs_engine", "deliver", "methods", "methods_engine",
                          "bundle_search", "bundle_engine", "uploads", "run_log", "peak_check",
                          "params", "audit")}


# ------------------------------------------------------------------------------------- F2 --
class UnrecordedInput(Tmp):
    """F2: the DE read a recordless output/search/; the recorded search beside it is not the
    DE's, and no reader names it."""

    def test_H_no_reader_names_the_lone_recorded_search(self):
        p, _d = fx.layout(self.tmp, {
            "search": (None, "report.parquet"),
            "search_sage": (tagged(fx.SAGE, "search_sage"), "report.parquet")},
            de_read=("search", "report.parquet"))
        got = self.picks(p)
        self.assertEqual(dict(got, methods_engine=None), NONE)
        c = session.search_choice(p)
        self.assertEqual(c["how"], "unrecorded_input")
        self.assertIn("output/search/report.parquet", c["why"])
        self.assertIn("output/search_sage/search_provenance.json", c["why"])

    def test_G_an_input_that_is_not_here_keeps_the_lone_record(self):
        # a laptop copy of a HIVE session: the DE's input names the HIVE path
        p, _d = fx.layout(self.tmp, {
            "search_sage": (tagged(fx.SAGE, "search_sage"), "report.parquet")})
        fx._de(p, "/hive/elsewhere/output/search_sage/report.parquet")
        self.assertEqual(session.search_choice(p)["how"], "only")
        self.assertEqual(self.picks(p)["deliver"], "search_sage")


# ------------------------------------------------------------------------------- F3 / F10 --
class OneAnswer(Tmp):
    """F3/F10: output/search/ and a second search, no DE input: no reader takes output/search/
    by default, and the client-facing ones say so."""

    def layout_E(self, de_read=None):
        return fx.layout(self.tmp, {
            "search": (tagged(fx.DIANN, "search"), "report.parquet"),
            "search_sage": (tagged(fx.SAGE, "search_sage"), "report.parquet")}, de_read=de_read)

    def test_E_no_reader_takes_output_search_by_default(self):
        for de_read in (None, "no-input"):
            shutil.rmtree(self.tmp)
            os.makedirs(self.tmp)
            p, _d = self.layout_E(de_read)
            self.assertEqual(self.picks(p), NONE, de_read)
            self.assertEqual(session.search_choice(p)["how"], "ambiguous")

    def test_E_deliver_skips_the_search_naming_both_records(self):
        import core_submission
        p, _d = self.layout_E()
        why = core_submission.delivery_search(p["output_dir"])[1]
        self.assertIn("output/search/search_provenance.json", why)
        self.assertIn("output/search_sage/search_provenance.json", why)
        items = {i["name"]: i for i in core_submission.plan_delivery(p["output_dir"])}
        self.assertEqual(items["search/"]["status"], "SKIPPED")

    def test_E_the_deposit_list_says_it(self):
        import make_deposit as md
        p, _d = self.layout_E()
        man = md.Manifest()
        md.build(p["session_dir"], man)
        (line,) = [ln for ln in man.lines if md.SEARCH_UPLOAD_PART in ln]
        self.assertTrue(line.startswith("[SKIPPED]"), line)
        self.assertIn("2 search records", line)

    def test_E_the_finalize_line_says_the_documents_describe_none(self):
        p, _d = self.layout_E()
        man = __import__("make_deposit").Manifest()
        session.search_record_step(p, man)
        (line,) = man.lines
        self.assertIn("2 search records", line)
        self.assertIn("describe none of them", line)

    def test_F10_the_bundle_reads_gathers_choice(self):
        import make_deposit as md
        p, d = fx.two_searches(self.tmp)
        f = md.gather(p["session_dir"])
        f["de_prov"] = {}                     # a reader that re-derived the choice would lose it
        self.assertEqual(md._chosen_search(f)[0],
                         os.path.join(d["sage"], "search_provenance.json"))


# --------------------------------------------------------------------------------- F4 / F5 --
class EngineVersion(Tmp):
    """F4: one engine/version rule -- a record's null version (conda-env Sage) is never the
    manifest's pin; F5: records and none chosen name no engine, and say so."""

    def run_methods(self, p):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                            "--raw", *session.read_raw_list(p["session_dir"]),
                            *FINALIZE_FLAGS, "--session", p["session_dir"],
                            "--workflow-manifest", p["workflow_manifest"],
                            "--out", os.path.join(self.tmp, "m.md")],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return read(os.path.join(self.tmp, "m.md"))

    def texts(self, p):
        import make_analysis_html as mah
        import notify_slack
        import session_docs
        de = self.de(p)
        a = argparse.Namespace(session=p["session_dir"], tables=p["de_dir"],
                               out=os.path.join(p["output_dir"], "Analysis_Report.html"))
        header = dict(mah.session_facts(a, de))["Search"]
        made = " ".join(b[1] for b in mah.how_made(a, de, None, []) if b[0] == "md")
        post = notify_slack.render(notify_slack.analysis_facts(p["session_dir"]))
        readme = session_docs.readme_md(session_docs.gather(p["session_dir"]))
        return header, made, json.dumps(post), readme

    def test_I_a_sage_record_without_a_version_is_never_given_the_pin(self):
        p, d = fx.two_searches(self.tmp)
        rec = json.loads(read(os.path.join(d["sage"], "search_provenance.json")))
        fx.write(os.path.join(d["sage"], "search_provenance.json"), json.dumps(
            dict(rec, version=None)))
        header, made, post, readme = self.texts(p)
        self.assertEqual(header, "Sage")
        for t in (made, post, readme):
            self.assertIn("Sage", t)
            self.assertNotIn("2.7.0", t)
        methods = self.run_methods(p)
        self.assertIn("Sage", methods)
        self.assertNotIn("2.7.0", methods)
        self.assertNotIn("pinned version", methods)

    def test_E_none_chosen_names_no_engine_and_says_so(self):
        p, _d = fx.layout(self.tmp, {
            "search": (tagged(fx.DIANN, "search"), "report.parquet"),
            "search_sage": (tagged(fx.SAGE, "search_sage"), "report.parquet")},
            de_read="no-input")
        header, made, post, readme = self.texts(p)
        self.assertEqual(header, "not recorded (2 search records, none chosen)")
        self.assertIn("Which search these numbers come from was not recorded (2 search "
                      "records, none chosen).", made)
        self.assertIn("2 search records, none chosen", post)
        self.assertIn("Search engine: not recorded (2 search records, none chosen)", readme)
        for t in (header, made, post, readme):
            self.assertNotIn("DIA-NN", t)

    def test_pin_one_search_and_no_record_read_as_before(self):
        p, d = fx.one_search(self.tmp)
        header, made, post, readme = self.texts(p)
        self.assertEqual(header, "DIA-NN 2.7.0")
        self.assertIn("Search engine: DIA-NN 2.7.0", readme)
        shutil.rmtree(d["search"])                # no record at all: the manifest's pin, tagged
        header, made, post, readme = self.texts(p)
        self.assertEqual(header, "DIA-NN 2.7.0 (workflow manifest pin)")    # Ruling 80 N2
        self.assertIn("Search engine: DIA-NN 2.7.0 (workflow manifest pin -- the search's own "
                      "record is not here)", readme)
        self.assertIn("workflow manifest pin", post)


# ------------------------------------------------------------------------------------- F6 --
import test_core_submission as tcs              # noqa: E402  (the deliver fixture, the fake API)


class DeliveryRecordsTheSearchFolder(tcs.DeliverBase):
    """F6: delivery.json says which search folder went, and the email draft names it."""

    def setUp(self):
        super().setUp()
        import finish_fixture
        shutil.rmtree(os.path.join(self.output, "search"))
        sage = os.path.join(self.output, "search_sage")
        fx.write(os.path.join(sage, "report.parquet"), "parquet")
        fx.write(os.path.join(sage, "search_provenance.json"), json.dumps(fx.SAGE))
        fx.write(os.path.join(self.output, "tables", "de_provenance.json"), json.dumps(
            dict(tcs.NORM_DECIDED, input=os.path.join(sage, "report.parquet"))))
        finish_fixture.stamp(self.output)

    def test_the_email_names_the_search_output(self):
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        dj = os.path.join(self.session, "delivery.json")
        self.assertEqual(tcs.load(dj)["search_folder"], {"folder": "search_sage", "why": None})
        draft = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        rc, out, p = tcs.run(["email-draft", "--summary", self.summary, "--delivery", dj,
                              "--share-url", "https://bioshare.example/bsx/", "--out", draft],
                             self.env)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertIn("  - the search engine output (in search_sage/)", read(draft))


# ------------------------------------------------------------------------- F3 (Methods) --
class MethodsTakeTheSession(Tmp):
    """make_methods.py --session gives the session's answer (finalize passes it when none is
    chosen): the search the DE read, or tagged blanks and why."""

    def methods(self, p):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                            "--raw", *session.read_raw_list(p["session_dir"]), *FINALIZE_FLAGS,
                            "--session", p["session_dir"], "--out", os.path.join(self.tmp, "m.md")],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return read(os.path.join(self.tmp, "m.md"))

    def test_the_search_the_de_read(self):
        p, _d = fx.two_searches(self.tmp)
        self.assertIn("Sage 0.14.7", self.methods(p))

    def test_H_the_de_read_a_search_with_no_record(self):
        p, _d = fx.layout(self.tmp, {
            "search": (None, "report.parquet"),
            "search_sage": (tagged(fx.SAGE, "search_sage"), "report.parquet")},
            de_read=("search", "report.parquet"))
        text = self.methods(p)
        self.assertIn("the search the differential expression read has no search record", text)
        self.assertNotIn("Sage", text)


# ------------------------------------------------------------------------------------- F7 --
class PeakCheckReason(Tmp):
    """F7: with none chosen, the peak check skips with the true reason."""

    def test_E_the_skip_names_both_records(self):
        import skyline_peaks
        p, _d = fx.layout(self.tmp, {
            "search": (tagged(fx.DIANN, "search"), "report.parquet"),
            "search_sage": (tagged(fx.SAGE, "search_sage"), "report.parquet")})
        a = argparse.Namespace(session=p["session_dir"], search=None, out=None, de_dir=None,
                               manifest=None, conditions=None, name=None)
        with self.assertRaises(skyline_peaks.Skip) as cm:
            skyline_peaks.make_plan(a, skyline_peaks.resolve_inputs(a))
        why = str(cm.exception)
        self.assertIn("2 search records", why)
        self.assertIn("output/search_sage/search_provenance.json", why)
        self.assertIn("--search", why)
        self.assertNotIn("there is no search folder", why)


# ------------------------------------------------------------------------------------- F8 --
class CommandsNameTheRule(Tmp):
    """F8: the skill's commands take the session's search, and a missing --search-prov file is
    refused."""

    def test_the_commands_do_not_hard_code_output_search(self):
        for rel in ("SKILL.md", os.path.join("references", "core-submissions.md")):
            text = read(os.path.join(SKILL, rel))
            self.assertNotIn("--search-out output/search", text, rel)
            self.assertNotIn("--search-prov <session>/output/search/", text, rel)

    def test_make_methods_refuses_a_missing_search_prov(self):
        p, _d = fx.two_searches(self.tmp)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                            "--raw", *session.read_raw_list(p["session_dir"]), *FINALIZE_FLAGS,
                            "--search-prov",
                            os.path.join(p["output_dir"], "search", "search_provenance.json"),
                            "--out", os.path.join(self.tmp, "m.md")],
                           capture_output=True, text=True, timeout=300)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("no such file", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "m.md")))

    def test_make_report_session_catalogs_the_search_the_de_read(self):
        p, _d = fx.two_searches(self.tmp)
        out = os.path.join(p["output_dir"], "OUTPUT_FILES.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_report.py"), "--out", out,
                            "--session", ".", "--de-dir", "output/tables"],
                           capture_output=True, text=True, timeout=120, cwd=p["session_dir"])
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        text = read(out)
        self.assertIn("output/search_sage/report.parquet", text)
        self.assertNotIn("output/search_dda/", text)


# ------------------------------------------------------------------------------------- F9 --
import test_record_run as trr                   # noqa: E402  (the run-log harness)


class RunLogSaysItCannotChoose(trr.Base):
    """F9: a session whose records leave the choice open: the run record says so."""

    def test_the_analysis_record_says_no_search_was_chosen(self):
        sess, _ = trr.make_session(self.d, quant_in_zip=False)
        for folder, rec in (("search_a", fx.DIANN), ("search_b", fx.SAGE)):
            fx.write(os.path.join(sess, "output", folder, "search_provenance.json"),
                     json.dumps(rec))
        self.run_it("analysis-done", "--session", sess)
        rec = self.read(self.only_folder(), "run_record.json")
        warns = " ".join(rec["analysis"]["warnings"])
        self.assertIn("2 search records", warns)
        self.assertIn("pass --out", warns)


if __name__ == "__main__":
    unittest.main()
