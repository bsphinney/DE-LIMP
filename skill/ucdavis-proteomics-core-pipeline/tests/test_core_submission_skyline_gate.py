"""
The peak check's folder (output/skyline/, skyline_peaks.py) goes to a collaborator only when the
check FINISHED (I31, I29; 2.14.2).

A skipped check writes its record into output/skyline/ with the reason. Before 2.14.2, deliver
shipped any skyline/ folder that existed: the delivery held a skyline/ folder with only the skip
record in it, marked [OK]; the README described pictures and a Skyline document that were not
there; and the client email promised "pictures of the peaks ... drawn with Skyline". Now
`core_submission.skyline_deliverable()` (skyline_peaks.record_state "ok", nothing else) decides
whether the folder goes, a folder that does not go is a [SKIPPED] line with the reason, and the
email names Skyline only when deliver recorded a delivered Skyline document (*_top_hits.sky.zip)
in delivery.json.

Stdlib only, no network: the staged submission and session come from test_core_submission's
DeliverBase.
"""
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))
sys.path.insert(0, HERE)

import core_submission as cs  # noqa: E402
import skyline_peaks as sp  # noqa: E402
from test_core_submission import DeliverBase, read, load, run, write  # noqa: E402

SKIP_REASON = "the search has no fixed carbamidomethyl cysteine, so Skyline cannot match its peptides"
EMAIL_SKYLINE = "drawn with Skyline"


def ok_record(**extra):
    rec = {"schema": sp.SCHEMA, "status": "ok", "name": "T",
           "proteins": [{"protein_group": "P1", "gene": "G1"}],
           "runs": {"imported": ["a", "b"], "failed": []},
           "results": {"compared": 4, "flagged": 0},
           "files": {"share_zip": "T_top_hits.sky.zip"}}
    rec.update(extra)
    return rec


def put_record(output_dir, rec):
    write(output_dir, os.path.join("skyline", sp.RECORD), json.dumps(rec))


class TestSkylineDeliverable(unittest.TestCase):
    """The one rule: True only for a finished check; otherwise the reason, in skyline_peaks' words."""

    def test_states(self):
        with tempfile.TemporaryDirectory() as td:
            os.makedirs(os.path.join(td, "skyline"))
            ok, why = cs.skyline_deliverable(td)
            self.assertFalse(ok)
            self.assertIn(sp.RECORD, why)
            put_record(td, sp.skipped_record(SKIP_REASON, "plan"))
            ok, why = cs.skyline_deliverable(td)
            self.assertFalse(ok)
            self.assertIn(SKIP_REASON, why)
            put_record(td, dict(ok_record(), status="planned"))
            ok, why = cs.skyline_deliverable(td)
            self.assertFalse(ok)
            self.assertIn(sp.STATE_WORDS["planned"], why)
            put_record(td, ok_record())
            self.assertEqual(cs.skyline_deliverable(td), (True, None))

    def test_an_unknown_state_is_said_not_a_crash(self):
        """Review of d3737cfd: STATE_WORDS[state] raised KeyError for a state it has no words for."""
        with tempfile.TemporaryDirectory() as td, \
                mock.patch.object(sp, "record_state", return_value=("archived", {})):
            ok, why = cs.skyline_deliverable(td)
        self.assertFalse(ok)
        self.assertIn("its record's state is 'archived'", why)


class TestPlanDelivery(unittest.TestCase):
    def test_skipped_peak_check_is_not_delivered(self):
        with tempfile.TemporaryDirectory() as td:
            put_record(td, sp.skipped_record(SKIP_REASON, "plan"))
            items = cs.plan_delivery(td)
        under = [i for i in items if i["name"].startswith("skyline/")]
        self.assertEqual([i for i in under if i["status"] == "OK"], [], under)
        folder = [i for i in under if i["name"] == "skyline/"]
        self.assertEqual(len(folder), 1, under)
        self.assertEqual(folder[0]["status"], "SKIPPED")
        self.assertIn(SKIP_REASON, folder[0]["reason"])

    def test_no_skyline_folder_says_nothing(self):
        """(pin) A session without the step has no skyline/ line at all, as before."""
        with tempfile.TemporaryDirectory() as td:
            items = cs.plan_delivery(td)
        self.assertEqual([i for i in items if i["name"].startswith("skyline")], [])


class TestDeliverSkyline(DeliverBase):
    def draft(self, delivery_json):
        out_md = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        rc, js, p = run(["email-draft", "--summary", self.summary, "--delivery", delivery_json,
                         "--share-url", "https://bioshare.example/bsx/", "--out", out_md], self.env)
        self.assertIn(rc, (0, 2), p.stdout + p.stderr)
        return read(out_md), js

    def deliver_ok(self):
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        return load(os.path.join(self.session, "delivery.json"))

    def test_readme_has_no_skyline_row_when_skipped(self):
        put_record(self.output, sp.skipped_record(SKIP_REASON, "plan"))
        dj = self.deliver_ok()
        readme = read(os.path.join(self.delivery, "README.md"))
        self.assertNotIn("skyline/", readme)
        self.assertNotIn("Skyline", readme)
        self.assertFalse(os.path.exists(os.path.join(self.delivery, "skyline")))
        manifest = read(os.path.join(self.delivery, "MANIFEST.txt"))
        line = [x for x in manifest.splitlines() if x.startswith("[SKIPPED] skyline/ -- ")]
        self.assertEqual(line, [f"[SKIPPED] skyline/ -- the peak check was skipped: {SKIP_REASON}"],
                         manifest)                   # said once, not "not delivered: ... skipped"
        self.assertNotIn("skyline", dj["contents"])
        self.assertEqual(dj["skyline"]["delivered"], False)
        self.assertIn(SKIP_REASON, dj["skyline"]["reason"])

    def test_email_does_not_name_skyline_for_a_skipped_check(self):
        put_record(self.output, sp.skipped_record(SKIP_REASON, "plan"))
        self.deliver_ok()
        text, _ = self.draft(os.path.join(self.session, "delivery.json"))
        self.assertNotIn(EMAIL_SKYLINE, text)

    def test_email_names_skyline_only_with_a_sky_zip(self):
        put_record(self.output, ok_record())
        write(self.output, "skyline/peaks_G1.png", "png")
        dj = self.deliver_ok()
        self.assertIn("skyline", dj["contents"])
        self.assertIsNone(dj["skyline"]["document"])
        text, _ = self.draft(os.path.join(self.session, "delivery.json"))
        self.assertNotIn(EMAIL_SKYLINE, text)
        # the README does not promise a document that is not there either
        readme = read(os.path.join(self.delivery, "README.md"))
        self.assertIn("`skyline/`", readme)
        self.assertNotIn("Skyline document", readme)

    def test_ok_peak_check_is_delivered(self):
        """(pin) A finished check goes, with its document, and the README and email name it."""
        put_record(self.output, ok_record())
        write(self.output, "skyline/peaks_G1.png", "png")
        write(self.output, "skyline/T_top_hits.sky.zip", "zip")
        items = {i["name"]: i["status"] for i in cs.plan_delivery(self.output)}
        self.assertEqual(items.get("skyline/peaks_G1.png"), "OK")
        self.assertEqual(items.get("skyline/T_top_hits.sky.zip"), "OK")
        self.assertNotIn("skyline/", items)
        self.deliver_ok()
        self.assertTrue(os.path.isfile(os.path.join(self.delivery, "skyline", "T_top_hits.sky.zip")))
        readme = read(os.path.join(self.delivery, "README.md"))
        self.assertIn("`skyline/T_top_hits.sky.zip`", readme)
        self.assertIn(sp.open_document_text(), readme)
        self.assertIn("drawn with Skyline, its tables and the Skyline document", readme)
        text, _ = self.draft(os.path.join(self.session, "delivery.json"))
        self.assertIn(EMAIL_SKYLINE, text)

    def test_delivery_json_records_the_delivered_document(self):
        """What the email (run on a laptop, from delivery.json alone) needs to know."""
        put_record(self.output, ok_record())
        write(self.output, "skyline/T_top_hits.sky.zip", "zip")
        dj = self.deliver_ok()
        self.assertEqual(dj["skyline"], {"delivered": True, "reason": None,
                                         "document": "skyline/T_top_hits.sky.zip"})
        # a session with no peak check: no `skyline` record at all
        os.remove(os.path.join(self.output, "skyline", "T_top_hits.sky.zip"))
        os.remove(os.path.join(self.output, "skyline", sp.RECORD))
        os.rmdir(os.path.join(self.output, "skyline"))
        rc, _, p = self.deliver("--apply", "--label", "t10")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertIsNone(load(os.path.join(self.session, "delivery.json"))["skyline"])

    def test_a_delivery_record_from_before_the_gate_names_no_skyline_and_says_why(self):
        """A delivery.json written before 2.14.2 has `skyline` in contents and no record of a
        document: the draft cannot know one was delivered, so it does not promise one, and says so."""
        dj = write(self.tmp, "old_delivery.json", json.dumps(
            {"mode": "analysis", "contents": ["Analysis_Report.html", "skyline"]}))
        text, js = self.draft(dj)
        self.assertNotIn(EMAIL_SKYLINE, text)
        self.assertTrue(any("skyline" in w for w in js["warnings"]), js["warnings"])


if __name__ == "__main__":
    unittest.main()
