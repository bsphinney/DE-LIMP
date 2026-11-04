#!/usr/bin/env python3
"""The finish gates (2.14.2, I35 + the 2026-10-31 amendments).

A Core submission (mstaff, 2026-10-29): the agent made the report with pandoc and COPIED it over
output/Analysis_Report.html so make_podcast.py would take it (it did); Methods, the
reproducibility bundle, finalize and deliver never ran, and the collaborator got a hand-made
folder. The beer re-analysis and three earlier HT submissions ended the same way: every one an
analysis with no traditional DE, and none with AGENTS.md, a README, Methods or the bundle.

What these tests pin:
  * report_files: the skill's record of its own report (report_provenance.json); a pre-2.14.2
    report (footer, no record) is accepted with a warning; anything else is refused;
  * make_podcast.py check / link / share / job refuse a report the skill did not write, or no
    Methods; link re-marks only a report that was the skill's before its edit; share and link
    WARN when the session is not finalized;
  * core_submission.py deliver (analysis, no-de) needs the report, Methods, the bundle and the
    finalize stamp; raw-only is never gated; --override-finish is recorded; bioshare send and
    email-draft refuse a delivery made from a session that was not finalized;
  * notify_slack.py analysis-done needs the stamp; board.sh only warns;
  * session.py finalize refuses without the report or the bundle, and a Claude Code analysis
    with no saved conversation (Review Focus 5); it writes output/FINALIZED.json last, and
    refuses to stamp when AGENTS.md cannot be written;
  * session.py status: the checklist, exit 3 until complete, 0 after.
Stdlib only, no network: the CoreOmics API is test_core_submission's local fake.
"""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import finish_fixture                          # noqa: E402  (a finished session's pieces)
import report_files as rf                      # noqa: E402
import test_core_submission as tcs             # noqa: E402  (the deliver fixture, the fake API)
import test_deposit_package as tdp             # noqa: E402  (a finished DIA-NN session)
import test_make_podcast as tmp_               # noqa: E402  (the podcast workspace)
import test_record_run as tdp_rr              # noqa: E402  (record_run's run-log fixture)
from job_env import job_env                    # noqa: E402

PY = sys.executable
SID = "0b6c2c1e-1111-4222-8333-444455556666"   # a Claude Code session id (shape only)


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def skill_report(out, text="<html><body><p>the skill's report</p></body></html>"):
    """output/Analysis_Report.html as make_analysis_html.py leaves it: the page and its record."""
    return finish_fixture.skill_report(out, text)


bundle = finish_fixture.bundle            # output/reproducibility/run_manifest.json (step 10)
stamp = finish_fixture.stamp              # output/FINALIZED.json (step 12)


# --------------------------------------------------------------------------- the record --
class SkillReportRecord(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.out = os.path.join(self._td.name, "output")

    def tearDown(self):
        self._td.cleanup()

    def test_the_skills_own_report_a_legacy_one_and_a_hand_made_one(self):
        html = skill_report(self.out)
        self.assertEqual(rf.skill_report_state(html)[0], "skill")
        write(html, "<html>pandoc</html>")                 # a Core submission: a copy over it
        state, why = rf.skill_report_state(html)
        self.assertEqual(state, "not_skill")
        self.assertIn("changed after make_analysis_html.py wrote it", why)
        os.remove(os.path.join(self.out, rf.SKILL_REPORT_RECORD))
        state, why = rf.skill_report_state(html)
        self.assertEqual(state, "not_skill")
        self.assertIn("not written by make_analysis_html.py", why)
        write(html, f"<html><footer>{rf.LEGACY_FOOTER}. Click a figure.</footer></html>")
        self.assertEqual(rf.skill_report_state(html)[0], "legacy_skill")   # Review Focus 1
        os.remove(html)
        self.assertEqual(rf.skill_report_state(html)[0], "missing")

    def test_make_analysis_html_marks_the_report_it_writes(self):
        session = self._td.name
        write(os.path.join(self.out, "AI_Analysis_Report.md"), "# Study\n\n## Overview\n\nText.\n")
        html = os.path.join(self.out, "Analysis_Report.html")
        r = subprocess.run([PY, os.path.join(SCRIPTS, "make_analysis_html.py"), "--session",
                            session, "--out", html, "--no-pdf"], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(json.loads(r.stdout)["skill_record"],
                         os.path.join(self.out, rf.SKILL_REPORT_RECORD))
        self.assertEqual(rf.skill_report_state(html)[0], "skill")
        rec = json.loads(read(os.path.join(self.out, rf.SKILL_REPORT_RECORD)))
        self.assertEqual((rec["schema"], rec["report"], rec["writer"]),
                         (rf.SKILL_REPORT_SCHEMA, "Analysis_Report.html", "make_analysis_html.py"))
        self.assertIn(rf.LEGACY_FOOTER, read(html), "the footer the legacy test reads is still there")

    def test_the_needs_are_one_definition(self):
        self.assertEqual(rf.FINISH_NEEDS["podcast"], ("report", "methods"))
        self.assertEqual(rf.FINISH_NEEDS["finalize"], ("report",))   # it makes the bundle (H1)
        self.assertEqual(rf.FINISH_NEEDS["deliver"],
                         ("report", "methods", "reproducibility", "finalized"))
        gaps = rf.finish_gaps(self.out, rf.FINISH_NEEDS["deliver"])
        self.assertEqual(len(gaps), 4)
        for g, cmd in zip(gaps, ("make_analysis_html.py --session", "make_methods.py",
                                 "provenance.py --outdir", "session.py finalize --dir")):
            self.assertIn(cmd, g)
        self.assertTrue(rf.refusal("x", gaps).endswith(
            'or re-run with --override-finish "<reason>"'))


# --------------------------------------------------------------------------- the podcast --
class PodcastGate(tmp_.Workspace):
    def setUp(self):
        super().setUp()
        skill_report(self.out)
        write(os.path.join(self.out, "methods.md"), "# Methods\n")

    def test_hand_made_report_is_refused_by_the_podcast(self):           # a Core submission
        write(os.path.join(self.out, "Analysis_Report.html"), "<html>pandoc</html>")
        os.remove(os.path.join(self.out, rf.SKILL_REPORT_RECORD))
        write(self.script, tmp_.script_text())
        rc, _, err = tmp_.run("check", self.script, "--source", self.report)
        self.assertEqual(rc, 2, err)
        self.assertIn("not written by make_analysis_html.py", err)
        self.assertIn("make_analysis_html.py --session", err)
        self.assertIn('or re-run with --override-finish "<reason>"', err)
        self.assertFalse(os.path.exists(os.path.join(self.pod, "check.txt")))

    def test_missing_methods_is_refused(self):
        os.remove(os.path.join(self.out, "methods.md"))
        write(self.script, tmp_.script_text())
        rc, _, err = tmp_.run("check", self.script, "--source", self.report)
        self.assertEqual(rc, 2, err)
        self.assertIn("output/methods.md does not exist", err)

    def test_legacy_report_is_accepted_with_a_warning(self):              # Review Focus 1
        os.remove(os.path.join(self.out, rf.SKILL_REPORT_RECORD))
        write(os.path.join(self.out, "Analysis_Report.html"),
              f"<html><body>{tmp_.REPORT_HTML}<footer>{rf.LEGACY_FOOTER}.</footer></body></html>")
        write(self.script, tmp_.script_text())
        rc, _, err = tmp_.run("check", self.script, "--source", self.report)
        self.assertEqual(rc, 0, err)
        self.assertIn("WARNING: Analysis_Report.html was made before skill 2.14.2", err)

    def test_check_override_is_recorded_for_render(self):
        write(os.path.join(self.out, "Analysis_Report.html"), "<html>pandoc</html>")
        write(self.script, tmp_.script_text())
        rc, _, err = tmp_.run("check", self.script, "--source", self.report,
                              "--override-finish", "staff checked the hand-made report")
        self.assertEqual(rc, 0, err)
        rec = tmp_.mp.check_finish_override(tmp_.mp.parse_script(self.script))
        self.assertEqual(rec["reason"], "staff checked the hand-made report")
        self.assertEqual(rec["tool"], "make_podcast.py check")

    def test_skill_report_passes_and_link_keeps_it_marked(self):
        tmp_.checked_podcast(self.out, self.report)
        html = os.path.join(self.out, "Analysis_Report.html")
        before = read(html)
        rc, _, err = tmp_.run("link", self.out)
        self.assertEqual(rc, 0, err)
        self.assertNotEqual(read(html), before, "link put its Listen card in")
        self.assertEqual(rf.skill_report_state(html)[0], "skill")
        rec = json.loads(read(os.path.join(self.out, rf.SKILL_REPORT_RECORD)))
        self.assertEqual(rec["writer"], "make_analysis_html.py, then make_podcast.py link")
        self.assertEqual(tmp_.run("link", self.out)[0], 0)              # again: still marked
        self.assertEqual(rf.skill_report_state(html)[0], "skill")

    def test_link_never_blesses_a_hand_copy(self):
        tmp_.checked_podcast(self.out, self.report)
        html = write(os.path.join(self.out, "Analysis_Report.html"), tmp_.REPORT_HTML)  # by hand
        rc, _, err = tmp_.run("link", self.out)
        self.assertEqual(rc, 2, err)
        self.assertIn("changed after make_analysis_html.py wrote it", err)
        rc, _, err = tmp_.run("link", self.out, "--override-finish", "the Core's own page")
        self.assertEqual(rc, 0, err)
        self.assertEqual(rf.skill_report_state(html)[0], "not_skill")
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertEqual(man["finish_override"]["reason"], "the Core's own page")

    def test_link_and_share_warn_when_the_session_is_not_finalized(self):   # Amendment 2.4
        tmp_.checked_podcast(self.out, self.report)
        rc, _, err = tmp_.run("link", self.out)
        self.assertEqual(rc, 0, err)
        self.assertIn("WARNING: this session is not finalized: run session.py finalize before "
                      "sending (AGENTS.md, README, Methods and the bundle come from it)", err)
        _, _, err = tmp_.run("share", self.out, "--no-browser")
        self.assertIn("this session is not finalized", err)
        stamp(self.out)
        rc, _, err = tmp_.run("link", self.out)
        self.assertEqual(rc, 0, err)
        self.assertNotIn("not finalized", err)

    def test_share_refuses_a_hand_made_report(self):
        tmp_.checked_podcast(self.out, self.report)
        write(os.path.join(self.pod, "podcast.m4a"), "not decoded: refused before")
        write(os.path.join(self.out, "Analysis_Report.html"), "<html>pandoc</html>")
        rc, _, err = tmp_.run("share", self.out, "--no-browser")
        self.assertEqual(rc, 2, err)
        self.assertIn("[make_podcast.py share] refusing", err)

    def test_job_refuses_a_hand_made_report(self):
        tmp_.checked_podcast(self.out, self.report)
        write(os.path.join(self.out, "Analysis_Report.html"), "<html>pandoc</html>")
        rc, _, err = tmp_.run("job", self.script, "--tts", "local")
        self.assertEqual(rc, 2, err)
        self.assertIn("[make_podcast.py job] refusing", err)


    def test_a_corrupt_override_line_is_kept_and_said(self):                 # F9
        write(os.path.join(self.out, "Analysis_Report.html"), "<html>pandoc</html>")
        write(self.script, tmp_.script_text())
        rc, _, err = tmp_.run("check", self.script, "--source", self.report,
                              "--override-finish", "staff checked the hand-made report")
        self.assertEqual(rc, 0, err)
        check = os.path.join(self.pod, "check.txt")
        text = read(check)
        start = text.index(tmp_.mp.FINISH_OVERRIDE_LINE)
        end = text.index("\n", start)
        write(check, text[:start] + tmp_.mp.FINISH_OVERRIDE_LINE + '{"reason": "staff ch'
              + text[end:])                                         # cut short: not JSON
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rec = tmp_.mp.check_finish_override(tmp_.mp.parse_script(self.script))
        self.assertIsNotNone(rec, "kept, not dropped")
        self.assertIn('{"reason": "staff ch', rec["reason"])
        self.assertIn("cannot be read", rec["unreadable"])


# --------------------------------------------------------------------------- deliver --
class DeliverGates(tcs.DeliverBase):
    """A finished session -- the skill's report, Methods, the bundle, the stamp -- that each test
    takes one piece from (Review Focus 4)."""

    def setUp(self):
        super().setUp()
        skill_report(self.output)
        write(os.path.join(self.output, "methods.md"), "# Methods\n")
        bundle(self.output)
        stamp(self.output)

    def delivery_json(self):
        return tcs.load(os.path.join(self.session, "delivery.json"))

    def test_a_finished_session_delivers_and_records_its_stamp(self):
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        d = self.delivery_json()
        self.assertEqual(d["finalized"]["finalized_at"], "2026-11-01T00:00:00+00:00")
        self.assertIsNone(d["finish_override"])
        self.assertEqual(rf.delivery_finish_state(d)[0], True)

    def test_deliver_analysis_refuses_without_the_bundle(self):
        os.remove(os.path.join(self.output, "reproducibility", "run_manifest.json"))
        for extra in ((), ("--apply",)):                                  # dry run included
            rc, out, p = self.deliver(*extra)
            self.assertEqual(rc, 2, p.stdout + p.stderr)
            self.assertIn("the reproducibility bundle", out["error"])
            self.assertIn("provenance.py --outdir", p.stderr)
        self.assertFalse(os.path.exists(self.delivery), "nothing copied")

    def test_deliver_refuses_a_hand_made_report_and_no_stamp(self):
        write(os.path.join(self.output, "Analysis_Report.html"), "<html>pandoc</html>")
        os.remove(os.path.join(self.output, rf.FINALIZED_RECORD))
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 2, p.stdout + p.stderr)
        self.assertEqual(len(out["finish_gaps"]), 2, out["finish_gaps"])
        self.assertIn("session.py finalize --dir", p.stderr)
        self.assertFalse(os.path.exists(self.delivery))

    def test_deliver_no_de_passes_with_report_methods_bundle_and_no_tables(self):
        import shutil
        shutil.rmtree(os.path.join(self.output, "tables"))
        os.remove(os.path.join(self.output, "figures", "pca.png"))
        rc, out, p = self.deliver("--mode", "no-de", "--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertTrue(out["verified"])
        self.assertEqual(self.delivery_json()["mode"], "no-de")

    def test_deliver_raw_only_is_not_gated(self):
        for rel in ("Analysis_Report.html", "methods.md", rf.FINALIZED_RECORD,
                    os.path.join("reproducibility", "run_manifest.json")):
            os.remove(os.path.join(self.output, rel))
        rc, out, p = self.deliver("--mode", "raw-only", "--include-raw", "yes", "--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertNotIn("finish_gaps", out)
        self.assertEqual(rf.delivery_finish_state(self.delivery_json())[0], True)

    def test_override_is_recorded(self):
        os.remove(os.path.join(self.output, rf.FINALIZED_RECORD))
        rc, out, p = self.deliver("--override-finish", "staff hand-checked the report", "--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        d = self.delivery_json()
        self.assertEqual(d["finish_override"]["reason"], "staff hand-checked the report")
        self.assertEqual(d["finish_override"]["tool"], "core_submission.py deliver")
        self.assertIsNone(d["finalized"])
        manifest = tcs.read(os.path.join(self.delivery, "MANIFEST.txt"))
        self.assertIn("[OVERRIDE] the skill's finish steps -- staff hand-checked the report; "
                      "delivered without: the finalize stamp (output/FINALIZED.json)", manifest)
        self.assertEqual(rf.delivery_finish_state(d)[0], True)


# ---------------------------------------------------------- the outward steps after deliver --
class OutwardSteps(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.summary, self.s = tcs.write_summary(self.tmp)

    def tearDown(self):
        self._tmp.cleanup()

    def delivery(self, **over):
        d = dict({"verified": True, "internal_id": "PROT_0807", "share_dir": tcs.SERVER_SHARE,
                  "mode": "analysis", "contents": ["Analysis_Report.html"]}, **over)
        return tcs.write(self.tmp, "delivery.json", json.dumps(d))

    def send(self, fake, *args):
        share_url = tcs.SHARES + "41/share/"
        fake.routes[("GET", tcs.SHARES)] = lambda q, b: (200, [{
            "id": 41, "submission": tcs.HEX, "bioshare_id": "bsx", "name": "x", "notes": "",
            "sub_folder": None, "link_to_path": self.s["share_dir"] + "/",
            "url": "https://bioshare.example/bsx/"}])
        fake.routes[("POST", share_url)] = lambda q, b: (200, {"ok": True})
        return tcs.run(["bioshare", "send", "--summary", self.summary, "--apply", *args],
                       tcs.child_env(self.tmp, fake.base))

    def test_bioshare_send_refuses_a_delivery_from_a_session_not_finalized(self):
        with tcs.FakeCoreOmics() as fake:
            rc, out, p = self.send(fake, "--delivery", self.delivery(finalized=None,
                                                                     finish_override=None))
            self.assertEqual(rc, 2, p.stdout + p.stderr)
            self.assertIn("not finalized", out["error"])
            self.assertEqual(fake.posts(), [])
            rc, out, p = self.send(fake, "--delivery", self.delivery(finalized=None,
                                                                     finish_override=None),
                                   "--override-finish", "the PI asked for it today")
            self.assertEqual(rc, 0, p.stdout + p.stderr)
            self.assertEqual(out["finish_override"]["reason"], "the PI asked for it today")
            rc, _, p = self.send(fake, "--delivery", self.delivery(
                finalized={"finalized_at": "2026-11-01T00:00:00+00:00"}, finish_override=None))
            self.assertEqual(rc, 0, p.stdout + p.stderr)
            rc, out, p = self.send(fake, "--delivery", self.delivery())        # before 2.14.2
            self.assertEqual(rc, 0, p.stdout + p.stderr)
            self.assertIn("from before skill 2.14.2", out["finish_warning"])
            self.assertEqual(len(fake.posts()), 3)

    def draft(self, delivery, *extra):
        out = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        if os.path.exists(out):
            os.remove(out)
        rc, js, p = tcs.run(["email-draft", "--summary", self.summary, "--delivery", delivery,
                             "--share-url", "https://bioshare.example/bsx/", "--out", out, *extra],
                            tcs.child_env(self.tmp))
        return rc, js, p, (read(out) if os.path.exists(out) else None)

    def test_email_draft_promises_only_a_finalized_delivery(self):
        rc, js, p, text = self.draft(self.delivery(finalized=None, finish_override=None))
        self.assertEqual(rc, 2, p.stdout + p.stderr)
        self.assertIsNone(text, "no draft for a delivery from a session that was not finalized")
        self.assertIn("not finalized", js["error"])
        rc, js, p, text = self.draft(self.delivery(finalized=None, finish_override=None),
                                     "--override-finish", "an interim look, agreed by phone")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertIn("NOT FINALIZED: drafted under --override-finish (an interim look, agreed "
                      "by phone)", text)
        rc, js, p, text = self.draft(self.delivery(
            finalized={"finalized_at": "2026-11-01T00:00:00+00:00"}))
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertNotIn("NOT FINALIZED", text)
        rc, js, p, text = self.draft(self.delivery(mode="raw-only", finalized=None,
                                                   contents=["MANIFEST.txt"]))
        self.assertEqual(rc, 0, p.stdout + p.stderr)                      # raw-only: not gated


    def test_a_draft_without_delivery_json_says_it_was_not_checked(self):     # F6
        out = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        rc, js, p = tcs.run(["email-draft", "--summary", self.summary, "--share-url",
                             "https://bioshare.example/bsx/", "--out", out],
                            tcs.child_env(self.tmp))
        self.assertIn("NOT CHECKED: no delivery.json, finalize state unknown", read(out))


# ---------------------------------------------------------------- finalize, status, Slack --
class FinalizeBase(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.p = tdp.dia_session(self.root, de=False)        # no DE: the sessions that went wrong
        self.sd = self.p["session_dir"]
        self.out = self.p["output_dir"]
        self.env = job_env(self.root)
        self.env.pop("CLAUDE_CODE_SESSION_ID", None)
        self.env["SKILL_CONFIG_DIR"] = os.path.join(self.root, "cfg")
        self.env["CLAUDE_CONFIG_DIR"] = os.path.join(self.root, "cc")
        # finalize builds a missing bundle (provenance.py): against a fake setup.json
        self.env["PROTEOMICS_PIPELINE_HOME"] = finish_fixture.fake_pipeline_home(self.root)

    def tearDown(self):
        self._td.cleanup()

    def finish(self):
        skill_report(self.out)
        bundle(self.out)

    def finalize(self, *extra, env=None):
        return subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "finalize", "--dir",
                               self.sd, "--no-deposit", *extra], capture_output=True, text=True,
                              env=env or self.env)

    def manifest(self):
        return read(self.p["manifest_txt"])

    def stamp_path(self):
        return os.path.join(self.out, rf.FINALIZED_RECORD)


class FinalizeGate(FinalizeBase):
    def test_finalize_refuses_without_the_skills_report_and_would_build_the_bundle(self):
        write(os.path.join(self.out, "Analysis_Report.html"), "<html>pandoc</html>")
        r = self.finalize("--zip")
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("not written by make_analysis_html.py", r.stderr)
        # the bundle is no gap: finalize builds it from the session's records (L5b review H1)
        refusal = r.stderr.split("refusing:", 1)[1]
        self.assertNotIn("the reproducibility bundle", refusal)
        self.assertIn("finalize builds it", r.stderr)
        self.assertFalse(os.path.exists(self.sd + ".zip"))
        self.assertFalse(os.path.exists(self.p["manifest_txt"]), "nothing written")
        self.assertFalse(os.path.exists(os.path.join(self.out, "reproducibility",
                                                     "run_manifest.json")),
                         "nothing built before a refusal")

    def test_finalize_writes_the_stamp(self):
        self.finish()
        r = self.finalize("--zip")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        res = json.loads(r.stdout)
        self.assertEqual(res["finalized"], self.stamp_path())
        st = json.loads(read(self.stamp_path()))
        self.assertEqual(st["schema"], rf.FINALIZED_SCHEMA)
        for rel in ("AGENTS.md", "README.md", "README.html", "output/methods.md", "MANIFEST.txt",
                    "output/reproducibility/run_manifest.json", "output/Analysis_Report.html"):
            self.assertEqual(st["outputs"][rel], rf._sha256(os.path.join(self.sd, rel)), rel)
        self.assertTrue(st["outputs"]["output/methods.docx"].startswith(("[SKIPPED]", ))
                        or len(st["outputs"]["output/methods.docx"]) == 64)
        self.assertEqual(st["zip"], self.sd + ".zip")
        self.assertIsNone(st["finish_override"])
        self.assertRegex(self.manifest(), r"\[OK\] +Finalize stamp \(output/FINALIZED\.json\)")
        self.assertEqual(rf.stamp_drift(self.sd, st), [])

    def test_finalize_refuses_to_stamp_when_agents_md_cannot_be_written(self):
        self.finish()
        stamp(self.out)                                  # an earlier finalize's stamp
        os.makedirs(os.path.join(self.sd, "AGENTS.md"))  # a folder in its place: not writable
        r = self.finalize()
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertFalse(os.path.exists(self.stamp_path()), "the earlier stamp is gone too")
        res = json.loads(r.stdout)
        self.assertIsNone(res["finalized"])
        self.assertTrue(any(x.startswith("AGENTS.md: ") for x in res["not_finalized"]), res)
        self.assertRegex(self.manifest(), r"\[SKIPPED\] Finalize stamp \(output/FINALIZED\.json\)"
                                          r" +-- not written: AGENTS\.md was not written")
        self.assertIn("NOT FINALIZED", r.stderr)
        self.assertRegex(self.manifest(), r"\[INFO\] +Slack notification|\[SKIPPED\] +Slack")

    def test_override_is_recorded_in_the_manifest_and_the_stamp(self):
        r = self.finalize("--override-finish", "a QC run: no report was asked for")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        man = self.manifest()
        self.assertRegex(man, r"\[OVERRIDE\] Finish check \(session\.py finalize\) +-- a QC run: "
                              r"no report was asked for -- went ahead without: the skill's report")
        st = json.loads(read(self.stamp_path()))
        self.assertEqual(st["finish_override"]["reason"], "a QC run: no report was asked for")
        # the report only: finalize built the bundle (L5b review H1)
        self.assertEqual(len(st["finish_override"]["missing"]), 1, st["finish_override"])
        self.assertTrue(st["outputs"]["output/Analysis_Report.html"].startswith("[SKIPPED]"))
        self.assertEqual(len(st["outputs"]["output/reproducibility/run_manifest.json"]), 64)

    def test_an_empty_reason_is_no_override(self):
        r = self.finalize("--override-finish", "  ")
        self.assertEqual(r.returncode, 2)
        self.assertIn("--override-finish needs a reason", r.stderr)


class FinalizeTranscript(FinalizeBase):
    """Review Focus 5 and I35: a Claude Code analysis owes its conversation."""

    def setUp(self):
        super().setUp()
        self.finish()

    def decisions(self, conversation=None):
        text = "# Decisions log\n\n## 2026-10-30 10:00 UTC -- Groups: ctrl vs trt\n- **Why:** asked\n"
        if conversation:
            text += f"- **Conversation:** `{conversation}`\n"
        write(os.path.join(self.p["logs_dir"], "decisions.md"), text)

    def test_finalize_refuses_a_claude_code_session_without_a_transcript(self):
        env = dict(self.env, CLAUDE_CODE_SESSION_ID=SID)        # nothing under CLAUDE_CONFIG_DIR
        r = self.finalize("--zip", env=env)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("this is a Claude Code analysis (CLAUDE_CODE_SESSION_ID is set here)",
                      r.stderr)
        self.assertIn("save_transcript.py --hive", r.stderr)
        self.assertIn('--no-transcript "<why>"', r.stderr)
        self.assertFalse(os.path.exists(self.sd + ".zip"))

    def test_finalize_on_hive_with_a_saved_transcript_passes(self):
        # hive_remote: no CLAUDE_CODE_SESSION_ID on HIVE; the laptop saved the conversation, and
        # the decisions log names it
        self.decisions(SID)
        write(os.path.join(self.sd, "logs", "conversation", SID + ".jsonl"), '{"type": "user"}\n')
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertRegex(self.manifest(), r"\[OK\] +Analysis conversation \(logs/conversation\)")

    def test_finalize_with_only_the_decisions_log_passes_with_info(self):
        self.decisions()                                 # no record of a Claude Code session
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertRegex(self.manifest(), r"\[INFO\] +Analysis conversation \(logs/conversation\)"
                                          r" +-- .*the decisions log \(logs/decisions\.md\) is the "
                                          r"record")

    def test_a_record_naming_a_claude_code_session_with_nothing_saved_is_refused(self):
        self.decisions(SID)                              # the record names one; nothing saved
        r = self.finalize()
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("logs/decisions.md names a Claude Code conversation", r.stderr)

    def test_the_transcript_hooks_record_from_init_counts(self):
        # session.py init, run inside Claude Code, recorded this folder in the hook's map; this
        # finalize runs where CLAUDE_CODE_SESSION_ID is not set
        key = "local:" + os.path.abspath(self.sd)
        write(os.path.join(self.env["SKILL_CONFIG_DIR"], "transcript_sessions.json"), json.dumps(
            {SID: {"sessions": {key: {"session_dir": os.path.abspath(self.sd), "hive": None,
                                      "recorded": "2026-11-01T00:00:00+00:00"}},
                   "segments": [{"session": key, "start": {"line": 0, "after_uuid": None}}],
                   "updated": "2026-11-01T00:00:00+00:00"}}))
        r = self.finalize()
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("the transcript hook's record of this folder", r.stderr)

    def test_no_transcript_reason_goes_to_the_manifest(self):
        env = dict(self.env, CLAUDE_CODE_SESSION_ID=SID)
        r = self.finalize("--no-transcript", "the laptop that ran it was reimaged", env=env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertRegex(self.manifest(), r"\[SKIPPED\] Analysis conversation \(logs/conversation\)"
                                          r" +-- --no-transcript: the laptop that ran it was "
                                          r"reimaged")
        self.assertEqual(json.loads(read(self.stamp_path()))["no_transcript"],
                         "the laptop that ran it was reimaged")


class Status(FinalizeBase):
    def status(self, *extra):
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "status", "--session",
                            self.sd, *extra], capture_output=True, text=True, env=self.env)
        return r

    def test_status_lists_what_a_no_de_session_that_never_finalized_still_owes(self):
        skill_report(self.out)
        r = self.status()
        self.assertEqual(r.returncode, 3, r.stdout + r.stderr)
        for line in ("  done     The skill's report (output/Analysis_Report.html)",
                     "  missing  The Methods (output/methods.md)",
                     "  missing  The reproducibility bundle",
                     "  missing  AGENTS.md (a guide for an AI agent)",
                     "  missing  README.html + README.md (for collaborators)",
                     "  missing  The finalize stamp (output/FINALIZED.json)",
                     "  done     The analysis conversation (logs/conversation) -- not a Claude "
                     "Code analysis",
                     "  done     The delivery (core_submission.py deliver) -- not a CoreOmics "
                     "submission"):
            self.assertIn(line, r.stdout)
        self.assertIn("provenance.py --outdir", r.stdout)
        self.assertIn("session.py finalize --dir", r.stdout)
        self.assertIn("NOT complete: 5 item(s) missing", r.stdout)
        bundle(self.out)
        self.assertEqual(self.finalize().returncode, 0)
        r = self.status()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("Complete: every item is done.", r.stdout)
        write(os.path.join(self.sd, "AGENTS.md"), "# hand-edited after finalize\n")
        r = self.status("--json")
        self.assertEqual(r.returncode, 3)
        items = {i["item"]: i for i in json.loads(r.stdout)["items"]}
        self.assertEqual(items["finalized"]["state"], "stale")
        self.assertIn("AGENTS.md", items["finalized"]["detail"])

    def test_a_submission_owes_its_delivery_and_an_interim_one_keeps_it_open(self):
        import submission_report
        submission_report.attach(self.sd, {"internal_id": "PROT_9812", "organism": "Mouse"})
        self.finish()
        self.assertEqual(self.finalize().returncode, 0)
        r = self.status("--json")
        self.assertEqual(r.returncode, 3, r.stdout)
        items = {i["item"]: i for i in json.loads(r.stdout)["items"]}
        self.assertEqual(items["delivery"]["state"], "missing")
        self.assertIn("core_submission.py deliver", items["delivery"]["command"])
        write(os.path.join(self.sd, "delivery.json"), json.dumps(
            {"mode": "interim", "verified": True, "owed": ["the report", "Methods"]}))
        r = self.status()
        self.assertEqual(r.returncode, 3)
        self.assertIn("interim delivered; owed: the report, Methods", r.stdout)
        write(os.path.join(self.sd, "delivery.json"), json.dumps(
            {"mode": "no-de", "verified": True, "delivered_at": "2026-11-01T00:00:00+00:00"}))
        r = self.status()
        self.assertEqual(r.returncode, 0, r.stdout)

    def test_a_missing_session_is_exit_2(self):
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "status", "--session",
                            os.path.join(self.root, "nope")], capture_output=True, text=True)
        self.assertEqual(r.returncode, 2)


class SlackAndBoard(FinalizeBase):
    def test_analysis_done_announces_only_a_finalized_session(self):
        env = dict(self.env, SKILL_SLACK="0")
        cmd = [PY, os.path.join(SCRIPTS, "notify_slack.py"), "analysis-done", "--session",
               self.sd, "--dry-run"]
        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("not sent", r.stderr)
        self.assertIn("finalize stamp", r.stderr)
        r = subprocess.run(cmd + ["--override-finish", "the PI wants a heads-up"],
                           capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("NOT FINALIZED", r.stdout)
        self.assertIn("the PI wants a heads-up", r.stdout)
        stamp(self.out)
        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertNotIn("NOT FINALIZED", r.stdout)

    def test_the_board_only_warns(self):
        env = {k: v for k, v in os.environ.items() if k not in ("BOARD_URL", "BOARD_KEY")}
        env["HOME"] = env["USERPROFILE"] = os.path.join(self.root, "home")
        os.makedirs(env["HOME"])
        cmd = ["bash", os.path.join(SCRIPTS, "board.sh"), "where", "7", "--session-folder",
               self.sd, "--report", os.path.join(self.out, "Analysis_Report.html")]
        r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=60)
        self.assertIn("[board] WARNING: this session is not finalized", r.stderr)
        self.assertIn("No key", r.stderr, "and it still runs the client (here: no key)")
        stamp(self.out)
        r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=60)
        self.assertNotIn("not finalized", r.stderr)
        r = subprocess.run(cmd[:5] + [self.sd], capture_output=True, text=True, env=env,
                           timeout=60)
        self.assertNotIn("not finalized", r.stderr, "the 'session made' row: no warning")


# ---------------------------------------------------- the L5a review's fix round (2.14.2) --
class StaleStampDeliver(tcs.DeliverBase):
    """F1 (a): a stamp counts only while what it stamped is unchanged (report_files.stamp_drift):
    finalize, then make_analysis_html again, then deliver used to pass."""

    def setUp(self):
        super().setUp()
        skill_report(self.output)
        write(os.path.join(self.output, "methods.md"), "# Methods\n")
        bundle(self.output)
        stamp(self.output)

    def test_a_report_made_again_after_finalize_holds_the_delivery(self):
        skill_report(self.output, "<html><body><p>made again after finalize</p></body></html>")
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 2, p.stdout + p.stderr)
        self.assertEqual(len(out["finish_gaps"]), 1, out["finish_gaps"])
        self.assertIn("changed since: output/Analysis_Report.html -- run session.py finalize "
                      "again", out["finish_gaps"][0])
        self.assertFalse(os.path.exists(self.delivery), "nothing copied")

    def test_a_stamp_with_no_outputs_is_refused(self):
        """L5b review N2: {schema: finalized/1, outputs: {}} -- the old DeliverBase fixture, or
        one written by hand -- never went stale, so it vouched for anything."""
        for outputs in ({}, None, {"AGENTS.md": "[SKIPPED] x"}):
            write(os.path.join(self.output, rf.FINALIZED_RECORD),
                  json.dumps({"schema": rf.FINALIZED_SCHEMA, "finalized_at": "2026-11-01",
                              "outputs": outputs}))
            rec, why = rf.finalized_record(self.output)
            self.assertIsNone(rec, outputs)
            self.assertIn("names none of", why) if not outputs else self.assertIn("README.md", why)
            rc, out, p = self.deliver()
            self.assertEqual(rc, 2, p.stdout + p.stderr)
            self.assertTrue(any("not a stamp session.py finalize wrote" in g
                                for g in out["finish_gaps"]), out["finish_gaps"])

    def test_no_stamp_says_it_may_be_an_older_finalize(self):
        os.remove(os.path.join(self.output, rf.FINALIZED_RECORD))
        rc, out, p = self.deliver()
        self.assertEqual(rc, 2, p.stdout + p.stderr)
        self.assertIn("no 2.14.2 finalize stamp, output/FINALIZED.json (finalized by an older "
                      "skill version, or not at all): run session.py finalize", p.stderr)


class DeliverFixtureStamp(tcs.DeliverBase):
    def test_the_fixture_stamps_with_real_checksums(self):
        """N2: DeliverBase's stamp is finish_fixture.stamp's, so it goes stale like a real one."""
        st = json.loads(read(os.path.join(self.output, rf.FINALIZED_RECORD)))
        self.assertEqual(st["outputs"]["output/methods.md"],
                         rf._sha256(os.path.join(self.output, "methods.md")))
        for rel in rf.STAMP_REQUIRED:
            self.assertIn(rel, st["outputs"])


class StaleStampFinalize(FinalizeBase):
    """F1 (b): a refused finalize used to exit before it removed the earlier stamp, so the Slack
    post, the board and share went on trusting it."""

    def announce(self):
        env = dict(self.env, SKILL_SLACK="0")
        return subprocess.run([PY, os.path.join(SCRIPTS, "notify_slack.py"), "analysis-done",
                               "--session", self.sd, "--dry-run"], capture_output=True,
                              text=True, env=env)

    def test_a_refused_finalize_removes_the_earlier_stamp(self):
        self.finish()
        self.assertEqual(self.finalize().returncode, 0)
        self.assertEqual(self.announce().returncode, 0)
        # a report copied over the skill's: finalize refuses (a missing bundle it would build)
        write(os.path.join(self.out, "Analysis_Report.html"), "<html>pandoc</html>")
        r = self.finalize()
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertFalse(os.path.exists(self.stamp_path()), "the earlier stamp is gone")
        r = self.announce()
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("no 2.14.2 finalize stamp", r.stderr)

    def test_an_argument_error_keeps_the_earlier_stamp(self):
        """L5b review L4: a typo in --code-dir (or an empty reason) is not a finalize that ran:
        "nothing written" must be true, the earlier stamp included."""
        self.finish()
        self.assertEqual(self.finalize().returncode, 0)
        before = read(self.stamp_path())
        for extra, words in ((("--code-dir", "/no/such/dir"), "--code-dir /no/such/dir: not found"),
                             (("--override-finish", "  "), "--override-finish needs a reason"),
                             (("--no-transcript", ""), "--no-transcript needs a reason")):
            r = self.finalize(*extra)
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            self.assertIn(words, r.stderr)
            self.assertEqual(read(self.stamp_path()), before, f"{extra}: the stamp was removed")
        self.assertEqual(self.announce().returncode, 0)

    def test_a_report_made_again_after_finalize_is_not_announced(self):
        self.finish()
        self.assertEqual(self.finalize().returncode, 0)
        skill_report(self.out, "<html><body><p>again</p></body></html>")
        r = self.announce()
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("changed since: output/Analysis_Report.html", r.stderr)


class StatusWaivers(FinalizeBase):
    """F2: `finalize --no-transcript "<why>"` and `--override-finish "<why>"` exit 0, so status
    must not then stay open on what they waived: it reads the stamp and says so."""

    def status(self, env=None):
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "status", "--session",
                            self.sd, "--json"], capture_output=True, text=True,
                           env=env or self.env)
        return r.returncode, {i["item"]: i for i in json.loads(r.stdout)["items"]}

    def test_no_transcript_is_waived_in_status(self):
        self.finish()
        env = dict(self.env, CLAUDE_CODE_SESSION_ID=SID)
        r = self.finalize("--no-transcript", "the laptop that ran it was reimaged", env=env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rc, items = self.status(env)
        self.assertEqual(rc, 0, json.dumps(items, indent=1))
        self.assertEqual(items["conversation"]["state"], "waived")
        self.assertEqual(items["conversation"]["detail"],
                         "waived (the laptop that ran it was reimaged)")

    def test_overridden_items_are_done_with_the_reason(self):
        r = self.finalize("--override-finish", "a QC run: no report was asked for")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rc, items = self.status()
        self.assertEqual(rc, 0, json.dumps(items, indent=1))
        self.assertEqual(items["report"]["state"], "overridden", items["report"])
        self.assertEqual(items["report"]["detail"], "overridden (a QC run: no report was asked for)")
        # made, not waived: finalize builds the Methods and (L5b review H1) the bundle
        self.assertEqual(items["reproducibility"]["state"], "present", items["reproducibility"])
        self.assertEqual(items["methods"]["state"], "present")
        self.assertIn("under --override-finish", items["finalized"]["detail"])

    def text_status(self, env=None):
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "status", "--session",
                            self.sd], capture_output=True, text=True, env=env or self.env)
        return r.returncode, r.stdout

    def test_status_never_says_every_item_is_done_while_one_is_overridden(self):
        """L5b review H1: an override read as "Complete: every item is done." """
        r = self.finalize("--override-finish", "a QC run: no report was asked for")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rc, out = self.text_status()
        self.assertEqual(rc, 0, out)
        self.assertNotIn("every item is done", out)
        last = out.strip().splitlines()[-1]
        self.assertRegex(last, r"^Complete, with \d+ overridden/waived: ")
        self.assertIn("The skill's report (output/Analysis_Report.html) (overridden: a QC run: no "
                      "report was asked for)", last)
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "status", "--session",
                            self.sd, "--json"], capture_output=True, text=True, env=self.env)
        js = json.loads(r.stdout)
        self.assertTrue(js["complete"])
        self.assertIn("report", [i["item"] for i in js["overridden_or_waived"]])

    def test_a_waived_conversation_is_named_in_the_summary(self):
        self.finish()
        env = dict(self.env, CLAUDE_CODE_SESSION_ID=SID)
        r = self.finalize("--no-transcript", "the laptop that ran it was reimaged", env=env)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rc, out = self.text_status(env)
        self.assertEqual(rc, 0, out)
        self.assertEqual(out.strip().splitlines()[-1],
                         "Complete, with 1 overridden/waived: The analysis conversation "
                         "(logs/conversation) (waived: the laptop that ran it was reimaged)")

    def test_a_complete_session_with_nothing_waived_says_so(self):
        self.finish()
        self.assertEqual(self.finalize().returncode, 0)
        rc, out = self.text_status()
        self.assertEqual((rc, out.strip().splitlines()[-1]), (0, "Complete: every item is done."))


class OverrideAnnounced(FinalizeBase):
    """F3: a finalize that went ahead under --override-finish was announced as complete: finalize
    posted without the reason, the stamp's own override was not read by analysis-done, and the
    [OVERRIDE] lines counted nowhere."""

    def test_finalize_under_override_posts_it_with_the_reason(self):
        import contextlib
        import io
        from unittest import mock
        import notify_slack
        import session
        argv = ["session.py", "finalize", "--dir", self.sd, "--no-deposit", "--override-finish",
                "a QC run: no report was asked for"]
        with mock.patch.object(sys, "argv", argv), mock.patch.dict(os.environ, self.env, clear=True), \
                mock.patch.object(notify_slack, "analysis_done",
                                  return_value=(False, "off (test)")) as posted, \
                contextlib.redirect_stdout(io.StringIO()):
            session.main()
        self.assertEqual(posted.call_args.kwargs.get("finish_override"),
                         "a QC run: no report was asked for")
        self.finish()
        with mock.patch.object(sys, "argv", argv[:5]), mock.patch.dict(os.environ, self.env, clear=True), \
                mock.patch.object(notify_slack, "analysis_done",
                                  return_value=(False, "off (test)")) as posted, \
                contextlib.redirect_stdout(io.StringIO()):
            session.main()
        self.assertIsNone(posted.call_args.kwargs.get("finish_override"), "nothing overridden")

    def test_analysis_done_reads_the_stamps_override(self):
        stamp(self.out, finish_override={"tool": "session.py finalize", "missing": ["x"],
                                         "reason": "the PI asked for it today", "at": "t"})
        r = subprocess.run([PY, os.path.join(SCRIPTS, "notify_slack.py"), "analysis-done",
                            "--session", self.sd, "--dry-run"], capture_output=True, text=True,
                           env=dict(self.env, SKILL_SLACK="0"))
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("NOT FINALIZED", r.stdout)
        self.assertIn("the PI asked for it today", r.stdout)

    def test_override_lines_are_counted_in_the_post(self):
        import notify_slack
        write(self.p["manifest_txt"], "Session export manifest\n=======\n\n"
              "[OVERRIDE] Finish check (session.py finalize)              -- why -- went ahead\n"
              "[SKIPPED] Deposit package                                 -- --no-deposit\n")
        f = notify_slack.analysis_facts(self.sd)
        self.assertEqual((f["export_skipped"], f["export_overridden"]), (1, 1))


class OverrideInTheRunLog(tdp_rr.Base):
    """F3: record_run.py's run log lists [OVERRIDE] lines beside the [SKIPPED] ones."""

    def test_the_run_log_lists_override_lines(self):
        out = tdp_rr.make_search(self.d)
        sess, _ = tdp_rr.make_session(self.d, quant_in_zip=True)
        with open(os.path.join(sess, "MANIFEST.txt"), "a", encoding="utf-8") as fh:
            fh.write("[OVERRIDE] Finish check (session.py finalize)  -- a QC run -- went ahead "
                     "without: the skill's report\n")
        self.run_it("analysis-done", "--session", sess, "--out", out)
        log = self.read(self.only_folder())
        self.assertIn("1 overridden (--override-finish)", log)
        self.assertIn("[OVERRIDE] Finish check (session.py finalize)", log)


class EmptyFinishFiles(unittest.TestCase):
    """F4: a 0-byte output/methods.md, or an empty or broken run_manifest.json, passed the
    podcast and deliver gates. The empty-file rule is make_deposit.ensure_methods' own
    (report_files.has_text); the bundle's record must be a JSON object with something in it."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.out = os.path.join(self._td.name, "output")
        skill_report(self.out)

    def tearDown(self):
        self._td.cleanup()

    def gap(self, item):
        return [g for g in rf.finish_gaps(self.out, (item,))]

    def test_an_empty_methods_md_is_no_methods(self):
        import make_deposit
        for text in ("", "  \n\n"):
            write(os.path.join(self.out, "methods.md"), text)
            self.assertEqual(len(self.gap("methods")), 1)
            self.assertIn("output/methods.md is empty", self.gap("methods")[0])
            self.assertFalse(make_deposit._has_text(os.path.join(self.out, "methods.md")))
        write(os.path.join(self.out, "methods.md"), "# Methods\n")
        self.assertEqual(self.gap("methods"), [])
        self.assertIs(make_deposit._has_text, rf.has_text, "one definition")

    def test_an_empty_or_broken_bundle_record_is_no_bundle(self):
        path = os.path.join(self.out, "reproducibility", "run_manifest.json")
        for text, words in (("", "cannot be read"), ("{}", "is empty"), ("[1]", "not a JSON "
                                                                               "object"),
                            ('{"files": [', "cannot be read")):
            write(path, text)
            gaps = self.gap("reproducibility")
            self.assertEqual(len(gaps), 1, text)
            self.assertIn(words, gaps[0], text)
        write(path, '{"files": []}')
        self.assertEqual(self.gap("reproducibility"), [])


class EmptyMethodsPodcast(tmp_.Workspace):
    def setUp(self):
        super().setUp()
        skill_report(self.out)

    def test_a_0_byte_methods_md_is_refused_by_the_podcast(self):
        write(os.path.join(self.out, "methods.md"), "")
        write(self.script, tmp_.script_text())
        rc, _, err = tmp_.run("check", self.script, "--source", self.report)
        self.assertEqual(rc, 2, err)
        self.assertIn("output/methods.md is empty", err)


class Nits(FinalizeBase):
    """F7, F8, F10 of the L5a review."""

    def test_the_gate_speaks_once_even_after_the_compute_node_re_exec(self):
        # F7: on a login node finalize runs itself again on a compute node, with
        # SKILL_FINISH_GATE_DONE set: the gate's warnings were printed there a second time
        skill_report(self.out)
        write(os.path.join(self.out, "Analysis_Report.html"),
              f"<html><footer>{rf.LEGACY_FOOTER}.</footer></html>")      # a pre-2.14.2 report
        os.remove(os.path.join(self.out, rf.SKILL_REPORT_RECORD))
        bundle(self.out)
        words = "WARNING: Analysis_Report.html was made before skill 2.14.2"
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(r.stderr.count(words), 1, r.stderr)
        r = self.finalize(env=dict(self.env, SKILL_FINISH_GATE_DONE=self.sd))   # the re-exec
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(r.stderr.count(words), 0, r.stderr)

    def test_an_unreadable_transcript_map_is_said(self):
        # F8: a map that cannot be read silently meant "no Claude Code record"
        import contextlib
        import io
        from unittest import mock
        import session
        write(os.path.join(self.env["SKILL_CONFIG_DIR"], "transcript_sessions.json"), "{broken")
        err = io.StringIO()
        with mock.patch.dict(os.environ, self.env, clear=True), contextlib.redirect_stderr(err):
            self.assertEqual(session.claude_code_record(self.sd), (None, None))
        self.assertIn("transcript_sessions.json", err.getvalue())
        self.assertIn("cannot be read", err.getvalue())

    def test_status_requires_what_the_stamp_requires(self):
        # F10: README.html is optional to the stamp (STAMP_OUTPUTS), so not required by status
        self.finish()
        os.makedirs(os.path.join(self.sd, "README.html"))      # cannot be written
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "status", "--session",
                            self.sd, "--json"], capture_output=True, text=True, env=self.env)
        self.assertEqual(r.returncode, 0, r.stdout)
        item = {i["item"]: i for i in json.loads(r.stdout)["items"]}["readme"]
        self.assertTrue(item["done"], item)
        self.assertIn("README.html not written (optional)", item["detail"])


class ExactCommands(FinalizeBase):
    """F11: a refusal names the real command where the session gives it -- no "..." -- and says
    where the rest comes from where it does not."""

    def test_the_methods_command_is_the_one_finalize_runs(self):
        import make_deposit
        gap = rf.finish_gaps(self.out, ("methods",))[0]
        cmd = gap.split(" -- make it with: ", 1)[1]
        self.assertNotIn("...", cmd)
        want = make_deposit.methods_command(make_deposit.gather(self.sd),
                                            os.path.join(self.out, "methods.md"))
        for arg in want[2:]:
            self.assertIn(arg, cmd)
        self.assertTrue(cmd.startswith("python3 scripts/make_methods.py --raw "), cmd)

    def test_a_long_raw_list_is_read_from_the_session(self):
        raws = [f"/data/runs/S{i:02d}.d" for i in range(1, 9)]
        write(self.p["raw_list"], "# raw\n" + "\n".join(raws) + "\n")
        gap = rf.finish_gaps(self.out, ("methods",))[0]
        self.assertIn("--raw $(grep -v '^#' " + self.p["raw_list"] + ")", gap)
        self.assertNotIn("S08.d", gap)

    def test_no_raw_list_says_so(self):
        os.remove(self.p["raw_list"])
        gap = rf.finish_gaps(self.out, ("methods",))[0]
        self.assertIn("(the session has no raw-file list to name them)", gap)
        gap = rf.finish_gaps(self.out, ("reproducibility",))[0]
        # every flag the session records, a placeholder for the raw files (L5b review H1)
        self.assertIn("--raw <the raw files>", gap)
        self.assertNotIn("SKILL.md step 10 lists", gap)
        self.assertNotIn("...", gap)


class ReferenceText(unittest.TestCase):
    """L5b review N1 and the [NOTE] line (M3): what the references say matches the code."""

    def test_outputs_md_says_docs_makes_the_stamp_stale_and_finalize_makes_the_bundle(self):
        with open(os.path.join(SKILL, "references", "outputs.md"), encoding="utf-8") as fh:
            text = " ".join(fh.read().split())
        self.assertIn("`session.py docs` after finalize makes the finalize stamp stale "
                      "(README.md, README.html and AGENTS.md are stamped): `deliver`, the Slack "
                      "post and `email-draft` then refuse until `session.py finalize` runs again.",
                      text)
        self.assertIn("builds the reproducibility bundle when it is missing or out of date", text)

    def test_step_12_relays_note_lines(self):
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            text = " ".join(fh.read().split())
        self.assertIn("**Read `MANIFEST.txt` and relay every `[SKIPPED]` and `[NOTE]` line with "
                      "its reason**", text)


class SkillMdText(unittest.TestCase):
    def setUp(self):
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            self.text = " ".join(fh.read().split())

    def test_golden_rule_9(self):
        self.assertIn(
            "9. **Stay on the standard path; leave it only when it is clearly better.** When you "
            "must (a step fails, the data need something no script does), record why with "
            "`log_decision.py` at that moment, and still run the finish steps: the report (step "
            "9, `make_analysis_html.py`), Methods (9d), the reproducibility bundle (10), finalize "
            "(12) and, for a Core submission, `core_submission.py deliver` (12b). Each refuses "
            "without what it needs: `make_podcast.py` the skill's own report and Methods; "
            "`finalize` the report (it makes the Methods and the bundle itself from the session's "
            "records); `deliver` the report, Methods, the bundle and "
            "the finalize stamp, which `bioshare send`, `email-draft` and the Slack post "
            "(`notify_slack.py analysis-done`) need too. A stamp counts only while nothing it "
            "stamped has changed. `--override-finish \"<reason>\"` goes ahead and records the "
            "reason. Before you tell the user an analysis is finished, run "
            "`python3 scripts/session.py status --session <session>` and show its result: exit 0 "
            "is complete, exit 3 lists each missing item with the command that makes it. Finish "
            "every missing item first.", self.text)

    def test_golden_rule_9_names_what_each_tool_needs(self):
        """(pin) F11: the rule says what report_files.FINISH_NEEDS says, tool by tool."""
        self.assertEqual(rf.FINISH_NEEDS["podcast"], ("report", "methods"))
        # finalize makes the bundle itself (L5b review H1): it needs only the report
        self.assertEqual(rf.FINISH_NEEDS["finalize"], ("report",))
        self.assertEqual(rf.FINISH_NEEDS["deliver"],
                         ("report", "methods", "reproducibility", "finalized"))
        self.assertEqual(rf.FINISH_NEEDS["announce"], ("finalized",))

    def test_step_12_says_what_finalize_refuses(self):
        self.assertIn(
            "finalize refuses when the report is missing, or when a Claude Code "
            "analysis has no saved conversation: save it (`save_transcript.py`; from a laptop for "
            "a HIVE session: `save_transcript.py --hive <S>`), or pass `--no-transcript "
            "\"<why>\"`. When the bundle (step 10) is missing or out of date, finalize builds it, "
            "with step 10's flags read from the session's records; it refuses, naming the exact "
            "`provenance.py` command, only when a record it needs is absent (no raw-file list).",
            self.text)


if __name__ == "__main__":
    unittest.main()
