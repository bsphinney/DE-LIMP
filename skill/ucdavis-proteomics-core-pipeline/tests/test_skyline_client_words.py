#!/usr/bin/env python3
"""
The peak check's words in a client report (2.14.2 C2 review: biologist H2, M5, L3).

H2: a skipped check on a Sage or FragPipe DE told the client "the quantities in this report are
    DIA-NN's", beside a reason naming "sage". Now the skipped section states the reason (its first
    letter capital, the engine by its label) and only that the peaks were not drawn.
L3: the reason reaches the page and the README row without server paths (make_methods.client_text,
    the Methods' rule).
M5: a check still running when the report was made told the client to "make the report again once
    they finish" -- words for staff. The page now says only that the check was not part of the
    report (REPORT_WORDS). The report's record (report_provenance.json) keeps the state it was
    made with, and a delivery never ships a finished skyline/ beside a report that says the check
    was not part of it (skyline_peaks.report_disagrees, through core_submission.
    skyline_deliverable): the report and the folder never contradict each other.
Hermetic: temp dirs.
"""
import io
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import core_submission as cs  # noqa: E402
import report_files as rf  # noqa: E402
import skyline_peaks as sp  # noqa: E402
from test_skyline_peaks import OK_REC  # noqa: E402


class Skipped(unittest.TestCase):
    def test_a_sage_search_names_its_engine_and_no_dianns_quantities(self):
        with tempfile.TemporaryDirectory() as td:
            search = os.path.join(td, "output", "search_sage")
            os.makedirs(search)
            with open(os.path.join(search, "search_provenance.json"), "w") as fh:
                json.dump({"engine": "sage", "version": "0.14.7"}, fh)
            out = os.path.join(td, "output", "skyline")
            with mock.patch("sys.stdout", new=io.StringIO()):
                self.assertEqual(sp.main(["plan", "--search", search, "--out", out]), 3)
            rec = sp.load_record(out)
            self.assertIn("the search was run with Sage;", rec["reason"])
            md = sp.report_section(rec)["md"]
            self.assertTrue(md.startswith("**Not run.** The search was run with Sage;"), md)
            self.assertNotIn("quantities in this report are DIA-NN's", md)
            self.assertNotIn("quantities", md)

    def test_the_reason_reaches_the_page_without_a_server_path(self):
        rec = sp.skipped_record("there is no search folder at /quobyte/somegroup/sessions/x/"
                                "output/search", "plan")
        md = sp.report_section(rec)["md"]
        self.assertNotIn("/quobyte/", md)
        self.assertIn("**Not run.** There is no search folder at search.", md)


class Unfinished(unittest.TestCase):
    def test_the_page_never_gives_staff_a_task(self):
        for state, words in sp.REPORT_WORDS.items():
            with self.subTest(state=state):
                self.assertTrue(words.startswith("**Not part of this report.**"), words)
                for task in ("make the report again", "run again", "needs to be",
                             "skyline_provenance.json", "jobs"):
                    self.assertNotIn(task, words)


class TheReportAndTheFolderAgree(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.out = os.path.join(self._t.name, "output")
        sky = os.path.join(self.out, "skyline")
        os.makedirs(sky)
        with open(sp.record_path(sky), "w") as fh:
            json.dump(dict(OK_REC, de_record=None, de_record_sha256=None), fh)
        self.assertEqual(sp.record_state(self.out)[0], "ok")
        self.html = os.path.join(self.out, "Analysis_Report.html")
        with open(self.html, "w") as fh:
            fh.write("<html>report</html>")

    def tearDown(self):
        self._t.cleanup()

    def test_a_report_made_while_the_jobs_ran_keeps_skyline_out_of_the_delivery(self):
        rf.mark_skill_report(self.html, "make_analysis_html.py", "x",
                             facts={"peak_check_state": "submitted"})
        why = sp.report_disagrees(self.out)
        self.assertIn("make the report again", why)
        ok, reason = cs.skyline_deliverable(self.out)
        self.assertFalse(ok)
        self.assertEqual(reason, why)

    def test_the_state_survives_a_later_writer(self):
        """make_podcast.py link records its own edit without the facts: they are kept"""
        rf.mark_skill_report(self.html, "make_analysis_html.py", "x",
                             facts={"peak_check_state": "planned"})
        rec = rf.mark_skill_report(self.html, "make_podcast.py link", "x")
        self.assertEqual(rec["peak_check_state"], "planned")
        self.assertFalse(cs.skyline_deliverable(self.out)[0])

    def test_a_report_made_with_the_finished_check_delivers_it(self):
        rf.mark_skill_report(self.html, "make_analysis_html.py", "x",
                             facts={"peak_check_state": "ok"})
        self.assertIsNone(sp.report_disagrees(self.out))
        self.assertEqual(cs.skyline_deliverable(self.out), (True, None))

    def test_a_report_record_without_the_state_claims_nothing(self):
        rf.mark_skill_report(self.html, "make_analysis_html.py", "x")
        self.assertEqual(cs.skyline_deliverable(self.out), (True, None))

    def test_report_data_says_its_state(self):
        self.assertEqual(sp.report_data(self.out)["state"], "ok")


if __name__ == "__main__":
    unittest.main()
