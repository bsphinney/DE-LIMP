#!/usr/bin/env python3
"""C2 errors_ux M1, M2, L2 (2.14.2): what "finalized" means after deliver, said the same way by
every step.

M1. report_files.delivery_finish_state returned ok on delivery.json `finalized` before it looked
at `finish_override`, and deliver filled `finalized` from finalized_record, which does not check
drift. So: finalize, re-make the report (stale stamp), `deliver --override-finish "client
waiting"` -- delivery.json held the stale stamp AND the override, and the email draft said
nothing, bioshare send warned nothing, and after a later finalize `status` said "Complete:
every item is done". A stamp written under `finalize --override-finish` was as silent. Slack
alone said NOT FINALIZED. Now the override ranks first (delivery_finish_override: deliver's,
else the stamp's), deliver records `finalized` only for a CURRENT stamp and carries the stamp's
own override, and email-draft, bioshare send and status all say NOT FINALIZED / overridden.

M2. A finalize that refused to stamp still ran record_run analysis-done, which logged the
analysis as finalized and complete. Now it is not logged (a [SKIPPED] line says why).

L2. `status` tested `owed` before `verified`: an interim delivery that failed verification read
"interim delivered". Now an unverified delivery is "unverified" first.
Stdlib only, no network.
"""
import json
import os
import shutil
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import report_files as rf                      # noqa: E402
import test_core_submission as tcs             # noqa: E402  (the deliver fixture, the fake API)
import test_deliver_interim as tdi             # noqa: E402  (a searched session, status())
import test_finish_gates as tfg                # noqa: E402  (FinalizeBase, the fixtures)

STALE = ("the finalize stamp (output/FINALIZED.json): finalized at 2026-11-01, but changed "
         "since: output/methods.md -- run session.py finalize again -- make it with: ...")


def override(reason, missing=(STALE,), tool="core_submission.py deliver"):
    return {"tool": tool, "missing": list(missing), "reason": reason, "at": "2026-11-02"}


class TheOverrideRanksFirst(unittest.TestCase):
    def test_a_stamp_recorded_beside_an_override_is_not_finalized(self):
        d = {"mode": "analysis", "finalized": {"finalized_at": "2026-11-01"},
             "finish_override": override("client waiting")}
        ok, why = rf.delivery_finish_state(d)
        self.assertIs(ok, True)
        self.assertTrue(why.startswith("NOT FINALIZED: delivered under --override-finish "
                                       "(client waiting)"), why)
        self.assertIn("the finalize stamp (output/FINALIZED.json)", why)
        self.assertEqual(rf.delivery_finish_override(d)["reason"], "client waiting")

    def test_a_stamp_written_under_an_override_is_not_finalized(self):
        d = {"mode": "no-de", "finish_override": None,
             "finalized": {"finalized_at": "2026-11-01",
                           "finish_override": override("Methods by hand",
                                                       ("output/methods.md: not written",),
                                                       "session.py finalize")}}
        ok, why = rf.delivery_finish_state(d)
        self.assertIs(ok, True)
        self.assertTrue(why.startswith("NOT FINALIZED"), why)
        self.assertIn("Methods by hand", why)

    def test_a_finalized_delivery_and_the_ungated_modes_are_unchanged(self):
        d = {"mode": "analysis", "finalized": {"finalized_at": "2026-11-01"},
             "finish_override": None}
        self.assertEqual(rf.delivery_finish_state(d),
                         (True, "the session was finalized at 2026-11-01"))
        self.assertIsNone(rf.delivery_finish_override(d))
        for mode in rf.UNGATED_DELIVERY_MODES:
            self.assertIsNone(rf.delivery_finish_override(
                {"mode": mode, "finish_override": override("x")}))


class DeliverRecordsOnlyACurrentStamp(tcs.DeliverBase):
    def setUp(self):
        super().setUp()
        tfg.skill_report(self.output)
        tfg.write(os.path.join(self.output, "methods.md"), "# Methods\n")
        tfg.bundle(self.output)
        tfg.stamp(self.output)

    def delivery_json(self):
        return tcs.load(os.path.join(self.session, "delivery.json"))

    def test_a_stale_stamp_delivered_under_override_is_not_recorded_as_finalized(self):
        tfg.write(os.path.join(self.output, "methods.md"), "# Methods\n\nedited after finalize\n")
        rc, out, p = self.deliver("--override-finish", "client waiting", "--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        d = self.delivery_json()
        self.assertIsNone(d["finalized"], "a stale stamp is not finalized")
        self.assertEqual(d["finish_override"]["reason"], "client waiting")
        ok, why = rf.delivery_finish_state(d)
        self.assertTrue(why.startswith("NOT FINALIZED"), why)

    def test_a_stamp_written_under_override_is_carried(self):
        tfg.stamp(self.output, finish_override=override("Methods by hand",
                                                        ("output/methods.md: by hand",),
                                                        "session.py finalize"))
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        d = self.delivery_json()
        self.assertEqual(d["finalized"]["finish_override"]["reason"], "Methods by hand")
        self.assertTrue(rf.delivery_finish_state(d)[1].startswith("NOT FINALIZED"))


class OutwardStepsSayIt(unittest.TestCase):
    """email-draft and bioshare send on a delivery made under an override (deliver's or the
    stamp's): NOT FINALIZED, never silent (as the Slack post). The fixture is OutwardSteps' (its
    methods borrowed, so its own tests do not run twice)."""
    setUp = tfg.OutwardSteps.setUp
    tearDown = tfg.OutwardSteps.tearDown
    delivery = tfg.OutwardSteps.delivery
    send = tfg.OutwardSteps.send
    draft = tfg.OutwardSteps.draft

    def test_the_draft_says_delivered_under_override(self):
        for d in (self.delivery(finalized={"finalized_at": "2026-11-01"},
                                finish_override=override("client waiting")),
                  self.delivery(finish_override=None, finalized={
                      "finalized_at": "2026-11-01",
                      "finish_override": override("client waiting",
                                                  ("output/methods.md: by hand",),
                                                  "session.py finalize")})):
            rc, js, p, text = self.draft(d)
            self.assertEqual(rc, 0, p.stdout + p.stderr)
            self.assertIn("NOT FINALIZED: DELIVERED UNDER --override-finish (client waiting) "
                          "without: ", text)

    def test_bioshare_send_warns_and_says_the_override(self):
        with tcs.FakeCoreOmics() as fake:
            rc, out, p = self.send(fake, "--delivery", self.delivery(
                finalized={"finalized_at": "2026-11-01"},
                finish_override=override("client waiting")))
            self.assertEqual(rc, 0, p.stdout + p.stderr)
            self.assertTrue(out["finish_warning"].startswith("NOT FINALIZED"), out)
            self.assertEqual(out["finish_override"]["reason"], "client waiting")
            self.assertIn("NOT FINALIZED", p.stderr)


class StatusReadsTheDelivery(tdi.InterimBase):
    def setUp(self):
        super().setUp()
        import submission_report
        os.makedirs(os.path.join(self.session, "input"), exist_ok=True)
        submission_report.attach(self.session, {"internal_id": "PROT_9999", "organism": "Human"})

    def put(self, **d):
        base = {"mode": "analysis", "verified": True, "owed": [],
                "delivered_at": "2026-11-02T00:00:00+00:00"}
        tcs.write(self.session, "delivery.json", json.dumps(dict(base, **d)))

    def test_a_delivery_under_override_is_overridden_never_simply_delivered(self):
        self.put(finalized={"finalized_at": "2026-11-01"},
                 finish_override=override("client waiting"))
        _, items = self.status()
        self.assertEqual(items["delivery"]["state"], "overridden", items["delivery"])
        self.assertTrue(items["delivery"]["done"])
        self.assertIn("NOT FINALIZED", items["delivery"]["detail"])
        self.assertIn("client waiting", items["delivery"]["detail"])

    def test_a_delivery_of_an_earlier_finalize_is_stale(self):
        tfg.stamp(self.output)                       # the session finalized again since
        self.put(finalized={"finalized_at": "2026-10-01", "sha256": "0" * 64},
                 finish_override=None)
        rc, items = self.status()
        self.assertEqual(rc, 3)
        self.assertEqual(items["delivery"]["state"], "stale", items["delivery"])
        self.assertIn("deliver again", items["delivery"]["detail"])

    def test_the_current_finalize_delivered_is_delivered(self):
        path = tfg.stamp(self.output)
        self.put(finalized={"finalized_at": "2026-11-01", "sha256": rf._sha256(path)},
                 finish_override=None)
        _, items = self.status()
        self.assertEqual(items["delivery"]["state"], "delivered", items["delivery"])

    def test_an_unverified_interim_is_unverified_not_delivered(self):
        """L2: verified is read before owed"""
        self.put(mode="interim", verified=False, owed=["the report"])
        rc, items = self.status()
        self.assertEqual(rc, 3)
        self.assertEqual(items["delivery"]["state"], "unverified", items["delivery"])
        self.assertIn("do not share it", items["delivery"]["detail"])


class ARefusedFinalizeIsNotLogged(tfg.FinalizeBase):
    """M2: no raw-file list, so finalize cannot write the Methods and refuses to stamp (NOT
    FINALIZED, exit 2): the Core run log must not record a finished analysis."""

    def test_no_analysis_done_for_a_finalize_that_did_not_stamp(self):
        os.remove(self.p["raw_list"])
        shutil.rmtree(self.p["search_out"])
        mine = tfg.write(os.path.join(self.root, "my_methods.md"), "# Methods\nBy hand.\n")
        self.finish()
        r = self.finalize("--methods", mine)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("NOT FINALIZED", r.stderr)
        self.assertFalse(os.path.exists(self.stamp_path()))
        man = self.manifest()
        self.assertRegex(man, r"\[SKIPPED\] Core run log +-- not logged as a finished analysis: "
                              r"the session is not finalized")


class StatusTakesDirAndHelpSaysWhatFinalizeDoes(tfg.FinalizeBase):
    """C2 docs M2 (code): `status` was the one subcommand without --dir (finalize, docs and
    which-search take it), so `session.py status --dir $S` was an argparse error; and --help
    said finalize refuses without the bundle, which it builds."""

    def test_status_takes_dir(self):
        import subprocess
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "status",
                            "--dir", self.sd, "--json"], capture_output=True, text=True,
                           env=self.env, timeout=120)
        self.assertIn(r.returncode, (0, 3), r.stderr)
        self.assertEqual(json.loads(r.stdout)["session_dir"], os.path.abspath(self.sd))

    def test_help_says_finalize_builds_the_bundle(self):
        import subprocess
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "--help"],
                           capture_output=True, text=True, timeout=60)
        flat = " ".join(r.stdout.replace("#", " ").split())
        self.assertIn("it BUILDS a missing or stale reproducibility bundle", flat)
        self.assertNotIn("own report or the reproducibility bundle", flat)


if __name__ == "__main__":
    unittest.main()
