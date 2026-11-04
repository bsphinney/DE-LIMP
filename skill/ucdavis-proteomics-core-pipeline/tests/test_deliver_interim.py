#!/usr/bin/env python3
"""`core_submission.py deliver --mode interim` (2.14.2, Amendment 3 to Task 13; three
earlier HT submissions, 2026-10-31).

One of those sessions followed the skill through the search, then left it at search QC: it
delivered the raw files with its own corrected-name script (deliver used the wrong original
names, I43) and the search folders with its own `deliver_results_*.sh` copy -- which bypassed the
qc_bracket gate, so step 8e was never prompted -- and "everything after step 8" became a backlog
line that survived compaction as a to-do. An interim delivery is now a mode, so an agent never
hand-copies search results. What these tests pin:
  * interim links the raw files under their staged (corrected) names and copies every search
    folder (session.find_search_provenance: output/search_dda/ ...), without chromatograms,
    engine intermediates or job scripts; it needs no report, Methods, bundle or finalize stamp,
    and refuses a session with no search record;
  * its README and MANIFEST say what is still owed (the report, the Methods, the bundle), and
    delivery.json records it as `owed`: `session.py status` stays open after it, and a later
    no-de (or analysis) delivery closes it;
  * the qc_bracket gate holds it like the other modes;
  * bioshare send and email-draft accept it, and the draft says what is to follow.
Stdlib only, no network: the CoreOmics API is test_core_submission's local fake.
"""
import json
import os
import shutil
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import core_submission as cs                   # noqa: E402
import finish_fixture                          # noqa: E402
import report_files as rf                      # noqa: E402
import test_core_submission as tcs             # noqa: E402  (the deliver fixture, the fake API)

OWED = list(cs.DELIVERY_MODES["interim"]["owed"])


class InterimBase(tcs.DeliverBase):
    """A staged submission whose session has searched (output/search_dda/) and stopped: no
    report, no Methods, no bundle, no finalize."""

    def setUp(self):
        super().setUp()
        for rel in ("Analysis_Report.html", "methods.md", rf.FINALIZED_RECORD,
                    os.path.join("reproducibility", "run_manifest.json")):
            os.remove(os.path.join(self.output, rel))
        dda = os.path.join(self.output, "search_dda")
        for rel, text in (("report.parquet", "parquet"), ("report.pg_matrix.tsv", "pg"),
                          ("report.log.txt", "log"), ("params.resolved.cfg", "--qvalue 0.01"),
                          ("quant/run1.quant", "quant"), ("lib.predicted.speclib", "lib"),
                          ("report_xic/x.parquet", "xic"), ("slurm-4242.out", "node log"),
                          ("run_search.sbatch", "#!/bin/bash"), ("submit.sh", "sbatch x")):
            tcs.write(dda, rel, text)
        tcs.write(dda, "search_provenance.json", json.dumps({"engine": "diann", "version": "2.7.0"}))
        self.dda = dda
        self.env.pop("CLAUDE_CODE_SESSION_ID", None)
        self.env["SKILL_CONFIG_DIR"] = os.path.join(self.tmp, "cfg")

    def delivery_json(self):
        return tcs.load(os.path.join(self.session, "delivery.json"))

    def interim(self, *extra, label="t1"):
        args = ["deliver", "--summary", self.summary, "--session", self.session, "--mode",
                "interim", "--label", label]
        return tcs.run(args + list(extra), self.env)

    def status(self):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "status",
                            "--session", self.session, "--json"], capture_output=True, text=True,
                           env=self.env)
        return r.returncode, {i["item"]: i for i in json.loads(r.stdout)["items"]}


class Interim(InterimBase):
    def test_interim_links_the_raw_and_copies_the_search_folders_and_writes_owed(self):
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertTrue(out["verified"], out)
        self.assertEqual(out["owed"], OWED)
        self.assertNotIn("finish_gaps", out)
        for rel in ("search_dda/report.parquet", "search_dda/report.pg_matrix.tsv",
                    "search_dda/report.log.txt", "search_dda/search_provenance.json",
                    "search_dda/params.resolved.cfg"):
            self.assertTrue(os.path.isfile(os.path.join(self.delivery, rel)), rel)
        for rel in ("search_dda/quant", "search_dda/lib.predicted.speclib", "search_dda/report_xic",
                    "search_dda/slurm-4242.out", "search_dda/run_search.sbatch",
                    "search_dda/submit.sh", "Analysis_Report.html", "tables", "figures",
                    "search/report.parquet"):
            self.assertFalse(os.path.lexists(os.path.join(self.delivery, rel)), rel)
        share_raw = os.path.join(self.share, "raw")
        self.assertEqual(sorted(os.listdir(share_raw)),
                         sorted(os.path.basename(f) for f in self.raw_files))
        d = self.delivery_json()
        self.assertEqual((d["mode"], d["owed"], d["differential_expression"]),
                         ("interim", OWED, False))
        self.assertEqual(rf.delivery_finish_state(d)[0], True)
        man = tcs.read(os.path.join(self.delivery, "MANIFEST.txt"))
        self.assertIn("[OK] search_dda/report.parquet", man)
        # the scheduler's logs and the Core's records: staff.is_core_internal (C2 security M1)
        self.assertIn("[SKIPPED] search_dda/slurm-4242.out -- Core-internal: a scheduler job log",
                      man)
        self.assertIn("[SKIPPED] search_dda/run_search.sbatch -- Core-internal: a job script", man)
        for item in OWED:
            self.assertIn(f"[SKIPPED] {item} -- not in this interim delivery: it follows in a "
                          "later one (--mode no-de or analysis)", man)
        readme = tcs.read(os.path.join(self.delivery, "README.md"))
        self.assertIn("**Still to come, in a later delivery:** " + "; ".join(OWED), readme)
        self.assertIn("`search_dda/`", readme)
        self.assertNotIn("Open `Analysis_Report.html` first", readme)

    def test_the_raw_links_keep_their_staged_names(self):
        # T8 (I43): a link staff made to correct a mis-named run keeps the corrected name
        wrong = tcs.raw(self.tmp, "tTOF_HT", "mar25", "wrong_sample.d")
        fixed = os.path.join(self.tmp, "corrected", "S07.d")
        os.makedirs(os.path.dirname(fixed))
        os.symlink(wrong, fixed)
        flist = tcs.write(self.tmp, "files_corrected.txt", "\n".join(self.raw_files + [fixed]))
        rc, _, p = tcs.run(["stage", "--summary", self.summary, "--files", flist, "--apply"],
                           self.env)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        share_raw = os.path.join(self.share, "raw")
        self.assertTrue(os.path.islink(os.path.join(share_raw, "S07.d")))
        self.assertFalse(os.path.lexists(os.path.join(share_raw, "wrong_sample.d")))

    def test_raw_can_be_left_out(self):
        rc, out, p = self.interim("--include-raw", "no", "--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.share, "raw")))
        self.assertIn("[SKIPPED] ../raw/ -- not requested (--include-raw no)",
                      tcs.read(os.path.join(self.delivery, "MANIFEST.txt")))
        # L5b review L3: the README promises only what was delivered
        title = tcs.read(os.path.join(self.delivery, "README.md")).splitlines()[0]
        self.assertNotIn("raw data", title.lower())
        self.assertIn("search results", title)

    def test_no_stage_record_hints_include_raw_no(self):
        """L5b review N3: an interim delivery carries the raw files by default, so with no stage
        record it said "raw data was requested" to someone who never asked: say why, and how to
        deliver the search results alone."""
        os.remove(os.path.join(self.service_project, cs.MARKER))
        os.remove(os.path.join(self.work, cs.MARKER))
        rc, out, p = self.interim("--force")
        self.assertEqual(rc, 2, p.stdout + p.stderr)
        self.assertIn("an interim delivery links the raw files unless --include-raw no",
                      out["error"])
        self.assertIn("--include-raw no", out["hint"])
        rc, out, p = self.interim("--force", "--include-raw", "no")
        self.assertEqual(rc, 0, p.stdout + p.stderr)

    def test_every_search_folder_goes(self):
        sage = os.path.join(self.output, "search_sage")
        tcs.write(sage, "results.sage.tsv", "psm")
        tcs.write(sage, "search_provenance.json", json.dumps({"engine": "sage"}))
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertTrue(os.path.isfile(os.path.join(self.delivery, "search_sage",
                                                    "results.sage.tsv")))
        self.assertTrue(os.path.isfile(os.path.join(self.delivery, "search_dda",
                                                    "report.parquet")))

    def test_a_session_with_no_search_record_is_refused(self):
        shutil.rmtree(self.dda)
        rc, out, p = self.interim()
        self.assertEqual(rc, 2, p.stdout + p.stderr)
        self.assertIn("no search record", out["error"])
        self.assertFalse(os.path.exists(self.delivery))

    def test_the_qc_bracket_gate_holds_it(self):
        self.set_qc("concern")
        for extra in ((), ("--apply",)):
            rc, out, p = self.interim(*extra)
            self.assertEqual(rc, 2, p.stdout + p.stderr)
            self.assertIn("held for staff review", out["error"])
        self.assertFalse(os.path.exists(self.delivery))
        self.ack()
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)

    def test_status_stays_open_after_interim_and_closes_after_a_no_de_delivery(self):
        import submission_report
        os.makedirs(os.path.join(self.session, "input"), exist_ok=True)
        submission_report.attach(self.session, {"internal_id": "PROT_0807", "organism": "Human"})
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        rc, items = self.status()
        self.assertEqual(rc, 3)
        self.assertEqual(items["delivery"]["state"], "owed")
        self.assertIn("interim delivered; owed: " + ", ".join(OWED), items["delivery"]["detail"])
        # the rest of the finish steps, then the no-de delivery beside the interim one
        shutil.rmtree(os.path.join(self.output, "tables"))
        os.remove(os.path.join(self.output, "figures", "pca.png"))
        finish_fixture.skill_report(self.output)
        tcs.write(self.output, "methods.md", "# Methods\n")
        finish_fixture.bundle(self.output)
        finish_fixture.stamp(self.output)
        for name in ("AGENTS.md", "README.md", "README.html"):
            tcs.write(self.session, name, rf.SKILL_DOC_MARKER + "\nfrom finalize\n")
        rc, out, p = tcs.run(["deliver", "--summary", self.summary, "--session", self.session,
                              "--mode", "no-de", "--label", "t2", "--include-raw", "no",
                              "--apply"], self.env)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(self.delivery_json()["owed"], [])
        rc, items = self.status()
        self.assertEqual(items["delivery"]["state"], "delivered", items["delivery"])
        self.assertEqual(rc, 0, json.dumps(items, indent=1))


class InterimLeavesConversionsOut(InterimBase):
    """L5b review L1: an interim delivery walks whole search folders, so a Sage search's mzML
    conversion of the raw files (a second copy of the raw data, often tens of GB) and the
    half-written files beside it were listed as OK."""

    def test_a_sage_searchs_mzml_folder_and_partial_files_stay_behind(self):
        sage = os.path.join(self.output, "search_sage")
        tcs.write(sage, "results.sage.tsv", "psm")
        tcs.write(sage, "report.parquet", "parquet")
        tcs.write(sage, "search_provenance.json", json.dumps({"engine": "sage"}))
        tcs.write(sage, "mzml/S01.mzML", "<mzML>" + "x" * 1000)
        tcs.write(sage, "mzml/.partial/S02.mzML", "<mzML")
        tcs.write(sage, "S03.mzML.partial", "x")
        tcs.write(sage, "lfq.parquet.part", "x")
        tcs.write(self.dda, "window_probe/probe.tsv", "x")
        tcs.write(self.dda, "radiant_lib/radiant_library.tsv", "x")
        rc, out, p = self.interim()                  # the dry run: the plan
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        ok = [i["name"] for i in out["items"] if i.get("status") == "OK"]
        self.assertIn(os.path.join("search_sage", "results.sage.tsv"), ok)
        bad = [n for n in ok if n.lower().endswith((".mzml", ".partial", ".part"))
               or "window_probe" in n or "radiant_lib" in n]
        self.assertEqual(bad, [], ok)
        skipped = {i["name"]: i["reason"] for i in out["items"] if i.get("status") == "SKIPPED"}
        self.assertIn("mzML conversion", skipped[os.path.join("search_sage", "mzml") + "/"])
        self.assertIn("half-written", skipped[os.path.join("search_sage", "S03.mzML.partial")])


class OwedCarriesForward(InterimBase):
    """L5b review M2: any later delivery cleared `owed` -- a raw-only re-send after an interim one
    (a late plate, the T8 name fix) made status say "delivered" although the report, Methods and
    bundle never went. Only a no-de or analysis delivery delivers them (references/
    core-submissions.md, SKILL.md 12b): any other carries what is owed forward."""

    def finish(self):
        """The session finalized after the interim delivery: report, Methods, bundle, stamp."""
        shutil.rmtree(os.path.join(self.output, "tables"))
        os.remove(os.path.join(self.output, "figures", "pca.png"))
        finish_fixture.skill_report(self.output)
        tcs.write(self.output, "methods.md", "# Methods\n")
        finish_fixture.bundle(self.output)
        for name in ("AGENTS.md", "README.md", "README.html"):
            tcs.write(self.session, name, rf.SKILL_DOC_MARKER + "\nfrom finalize\n")
        finish_fixture.stamp(self.output)

    def test_a_raw_only_delivery_after_an_interim_one_leaves_it_owed(self):
        import submission_report
        os.makedirs(os.path.join(self.session, "input"), exist_ok=True)
        submission_report.attach(self.session, {"internal_id": "PROT_0807", "organism": "Human"})
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.finish()
        rc, out, p = tcs.run(["deliver", "--summary", self.summary, "--session", self.session,
                              "--mode", "raw-only", "--label", "t3", "--apply"], self.env)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(out["owed"], OWED, "carried forward from the interim delivery")
        self.assertEqual(self.delivery_json()["owed"], OWED)
        rc, items = self.status()
        self.assertEqual(rc, 3, json.dumps(items, indent=1))
        self.assertEqual(items["delivery"]["state"], "owed")
        self.assertIn("raw-only delivered; owed: ", items["delivery"]["detail"])
        # the report's delivery closes it
        rc, out, p = tcs.run(["deliver", "--summary", self.summary, "--session", self.session,
                              "--mode", "no-de", "--label", "t4", "--include-raw", "no",
                              "--apply"], self.env)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(self.delivery_json()["owed"], [])
        rc, items = self.status()
        self.assertEqual(items["delivery"]["state"], "delivered", items["delivery"])

    def test_a_second_interim_delivery_owes_the_same_once(self):
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        rc, out, p = self.interim("--apply", label="t2")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(out["owed"], OWED)


class InterimOutward(unittest.TestCase):
    """bioshare send and email-draft run where the key is, from delivery.json alone."""

    def setUp(self):
        import tempfile
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.summary, self.s = tcs.write_summary(self.tmp)
        self.delivery = tcs.write(self.tmp, "delivery.json", json.dumps({
            "verified": True, "internal_id": "PROT_0807", "share_dir": tcs.SERVER_SHARE,
            "mode": "interim", "contents": ["MANIFEST.txt", "README.md", "search_dda"],
            "finalized": None, "finish_override": None, "owed": OWED,
            "raw_links": {"created": 2}}))

    def tearDown(self):
        self._tmp.cleanup()

    def test_bioshare_send_accepts_an_interim_delivery(self):
        with tcs.FakeCoreOmics() as fake:
            fake.routes[("GET", tcs.SHARES)] = lambda q, b: (200, [{
                "id": 41, "submission": tcs.HEX, "bioshare_id": "bsx", "name": "x", "notes": "",
                "sub_folder": None, "link_to_path": self.s["share_dir"] + "/",
                "url": "https://bioshare.example/bsx/"}])
            fake.routes[("POST", tcs.SHARES + "41/share/")] = lambda q, b: (200, {"ok": True})
            rc, out, p = tcs.run(["bioshare", "send", "--summary", self.summary, "--delivery",
                                  self.delivery, "--apply"], tcs.child_env(self.tmp, fake.base))
            self.assertEqual(rc, 0, p.stdout + p.stderr)
            self.assertEqual(len(fake.posts()), 1)

    def test_the_email_draft_says_what_is_to_follow(self):
        out_md = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        rc, js, p = tcs.run(["email-draft", "--summary", self.summary, "--delivery",
                             self.delivery, "--share-url", "https://bioshare.example/bsx/",
                             "--out", out_md], tcs.child_env(self.tmp))
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        text = tcs.read(out_md)
        self.assertNotIn("NOT FINALIZED", text)
        self.assertIn("Subject: Raw data and search results for PROT_0807", text)
        self.assertIn("Still to come, in a later delivery:", text)
        for item in OWED:
            self.assertIn(item, text)
        self.assertIn("Your raw instrument files are in the raw folder", text)
        self.assertNotIn("Please start with Analysis_Report.html", text)
        self.assertNotIn("the wording, including the instrument grant, is at the end of "
                         "methods.md", text)

    def draft(self, **fields):
        rec = json.loads(tcs.read(self.delivery))
        rec.update(fields)
        tcs.write(self.tmp, "delivery.json", json.dumps(rec))
        out_md = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        rc, js, p = tcs.run(["email-draft", "--summary", self.summary, "--delivery",
                             self.delivery, "--share-url", "https://bioshare.example/bsx/",
                             "--out", out_md], tcs.child_env(self.tmp))
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        return tcs.read(out_md)

    def test_an_interim_draft_without_raw_files_does_not_promise_them(self):
        """L5b review L3: --include-raw no -- the subject said "Raw data and search results"."""
        text = self.draft(raw_links={})
        self.assertIn("Subject: Search results for PROT_0807", text)
        self.assertNotIn("raw data", text.split("Subject:", 1)[1].splitlines()[0].lower())
        self.assertNotIn("Your raw instrument files", text)

    def test_a_raw_only_draft_after_an_interim_one_says_what_is_still_owed(self):
        """M2's carried `owed`: a raw-only re-send must not say no search was run."""
        text = self.draft(mode="raw-only", raw_links={"created": 2})
        self.assertIn("Subject: Raw data for PROT_0807 is ready", text)
        self.assertNotIn("no search or analysis was run", text)
        self.assertIn("Still to come, in a later delivery:", text)
        for item in OWED:
            self.assertIn(item, text)


if __name__ == "__main__":
    unittest.main()
