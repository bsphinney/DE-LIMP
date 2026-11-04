#!/usr/bin/env python3
"""C2 proteomics M2 (2.14.2): the Sage LFQ check said "ok" while the DE saw a small part of what
was identified.

Two timsTOF DDA-PASEF runs through the default quantifying Sage route (the .d converted in the
job): Sage identified 753 peptides and logged 289 target MS1 peaks (38%: the check passed), but
the adapter kept 65 precursors in 34 protein groups -- 8.6% of the identified peptides, below
the check's own 10% "LFQ collapsed" line -- and the DE ran on them with nothing said.

Now, after the adapter writes sage_adapt.json, sage_lfq_check.add_coverage puts
quantified_coverage (sage_adapt's target_precursors_kept against target_peptides_1pct, the one
place both counts are made) into the check's record and search_provenance.json, and WARNS with
the numbers when it is under LOW_PEAK_FRACTION. A warning, never a gate. AUDIT.md reads it.
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

import sage_lfq_check as slc                   # noqa: E402

try:
    import pyarrow  # noqa: F401
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False


class TheNumbers(unittest.TestCase):
    def test_the_reviewed_case_is_low(self):
        cov = slc.quantified_coverage({"target_precursors_kept": 65, "target_peptides_1pct": 753})
        self.assertEqual((cov["target_precursors_quantified"],
                          cov["target_peptides_identified_1pct"]), (65, 753))
        self.assertAlmostEqual(cov["fraction"], 0.0863, places=4)
        self.assertTrue(cov["low"])
        self.assertEqual(cov["low_below"], slc.LOW_PEAK_FRACTION)       # the one line
        note = slc.coverage_note(cov)
        for words in ("Only 65 target precursors were quantified", "753 target peptides",
                      "8.6%", "not a refusal"):
            self.assertIn(words, note)

    def test_a_healthy_search_is_not_low_and_missing_counts_say_nothing(self):
        self.assertFalse(slc.quantified_coverage({"target_precursors_kept": 400,
                                                  "target_peptides_1pct": 753})["low"])
        for adapt in (None, {}, {"target_precursors_kept": 5}, {"target_peptides_1pct": 0,
                                                                "target_precursors_kept": 5}):
            self.assertIsNone(slc.quantified_coverage(adapt), adapt)


class TheRecordWarns(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "sage"}, fh)

    def ok_record(self):
        slc.record(self.d, {"check": "sage_lfq_window", "status": "ok",
                            "message": "Sage LFQ window fits the data.",
                            "runs_outside_lfq_window": [], "target_peptides_1pct": 753})

    def test_low_coverage_turns_ok_into_a_warning_with_the_numbers(self):
        self.ok_record()
        rec = slc.add_coverage(self.d, {"target_precursors_kept": 65,
                                        "target_peptides_1pct": 753})
        self.assertEqual(rec["status"], "warn")
        self.assertIn("Only 65 target precursors were quantified", rec["message"])
        self.assertFalse(slc.gated(rec), "a warning, never a gate")
        on_disk = slc.load(self.d)
        self.assertEqual(on_disk["quantified_coverage"]["target_precursors_quantified"], 65)
        with open(os.path.join(self.d, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["sage_lfq_check"]["status"], "warn")
        status, msg, detail, plain = slc.audit_finding(on_disk)
        self.assertEqual(status, "WARN")
        self.assertIn("753 target peptides", msg)
        # the client's words: the settings DID fit; few peptides were quantified
        self.assertIn("The quantification settings fit these runs", plain)
        self.assertIn("65 quantified peptide signals for 753 identified peptides", plain)
        self.assertEqual(detail["quantified_coverage"]["target_precursors_quantified"], 65)
        again = slc.add_coverage(self.d, {"target_precursors_kept": 65,
                                          "target_peptides_1pct": 753})
        self.assertEqual(again["message"].count("Only 65"), 1, "said once")

    def test_a_healthy_coverage_is_recorded_and_stays_ok(self):
        self.ok_record()
        rec = slc.add_coverage(self.d, {"target_precursors_kept": 600,
                                        "target_peptides_1pct": 753})
        self.assertEqual(rec["status"], "ok")
        self.assertFalse(rec["quantified_coverage"]["low"])
        self.assertNotIn("coverage_warning", rec)

    def test_an_unchecked_record_is_left_alone(self):
        slc.record(self.d, {"status": "unchecked", "message": "x"})
        self.assertIsNone(slc.add_coverage(self.d, {"target_precursors_kept": 1,
                                                    "target_peptides_1pct": 753}))
        self.assertNotIn("quantified_coverage", slc.load(self.d))


@unittest.skipUnless(HAVE_ARROW, "pyarrow is needed for the Sage fixtures")
class TheAdapterSaysIt(unittest.TestCase):
    """run_search.adapt_sage on a search that identified 200 peptides per run and whose
    lfq.parquet keeps 2 target precursors: the record warns after the adapter."""

    def test_adapt_sage_records_and_warns(self):
        import run_search
        import test_sage_adapt as tsa
        import test_sage_lfq_check as fx
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        fx.write_sage_outputs(d, {"a.mzML": 0.5, "b.mzML": 0.5}, peaks=150)
        slc.record(d, slc.check(d))
        self.assertEqual(slc.load(d)["status"], "ok", slc.load(d).get("message"))
        tsa.write_lfq(os.path.join(d, "lfq.parquet"))           # 2 target precursors pass
        run_search.adapt_sage(d)
        rec = slc.load(d)
        self.assertEqual(rec["status"], "warn", rec.get("message"))
        cov = rec["quantified_coverage"]
        self.assertEqual(cov["target_precursors_quantified"], 2)
        self.assertEqual(cov["protein_groups_quantified"], 2)     # P1, P2 (sage_adapt's count)
        with open(os.path.join(d, "sage_adapt.json")) as fh:
            self.assertEqual(json.load(fh)["protein_groups_kept"], 2)
        self.assertTrue(cov["low"])
        self.assertIn("Only 2 target precursors in 2 protein groups were quantified",
                      rec["message"])
        self.assertFalse(slc.gated(rec))


class TheWordsForEachReader(unittest.TestCase):
    """Verify round (proteomics): the client's words said the amounts "rest on few peptides",
    when the bigger effect is that most identified proteins get no amount at all; and AUDIT.md,
    which goes into the client's bundle, ended with a staff task. Now the client's words give
    the protein groups with amounts (sage_adapt's protein_groups_kept), and the task is said to
    staff only (stderr), never in the record's message, AUDIT.md or the client table."""

    LOW = {"target_precursors_kept": 65, "target_peptides_1pct": 753, "protein_groups_kept": 34}

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        slc.record(self.d, {"check": "sage_lfq_window", "status": "ok",
                            "message": "Sage LFQ window fits the data.",
                            "runs_outside_lfq_window": []})

    def test_the_client_words_say_most_proteins_have_no_amount(self):
        import contextlib
        import io
        with contextlib.redirect_stderr(io.StringIO()):
            slc.add_coverage(self.d, self.LOW)
        _s, _m, _d, plain = slc.audit_finding(slc.load(self.d))
        self.assertIn("the comparison has amounts for 34 protein groups", plain)
        self.assertIn("from 65 quantified peptide signals for the 753 peptides identified", plain)
        self.assertIn("Many identified proteins therefore have no amount in the comparison",
                      plain)
        self.assertIn("the identifications are not affected", plain)

    def test_the_staff_task_is_said_to_staff_only(self):
        import contextlib
        import io
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            rec = slc.add_coverage(self.d, self.LOW)
        self.assertIn("STAFF: " + slc.COVERAGE_STAFF_TASK, err.getvalue())
        status, msg, detail, plain = slc.audit_finding(slc.load(self.d))
        for text in (rec["message"], rec["coverage_warning"], msg, plain, json.dumps(detail)):
            self.assertNotIn("another engine", text)
            self.assertNotIn("compare", text.lower())
        self.assertIn("Only 65 target precursors in 34 protein groups were quantified", msg)


if __name__ == "__main__":
    unittest.main()
