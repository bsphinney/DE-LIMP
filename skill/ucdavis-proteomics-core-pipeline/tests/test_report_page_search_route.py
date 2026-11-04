#!/usr/bin/env python3
"""
The report page itself says when its search was a first pass or for identifications only
(2.14.2 C2 review: biologist M1).

On the no-DE route the audit (step 8c) does not run, so the first-pass caveat -- which reached
the page only through AUDIT.md -- never did, and an identifications-only report never said on the
page that there were no quantities. Both are now fixed callouts in "Results at a glance", like
search_measurement_note, never left to the writer:
  - first_pass_note: dda_first_pass.caution decides, in AUDIT.md's plain-language words;
  - identifications_only_note: sage_lfq_check.search_quantified (the predicate the DE gate
    shares) decides, with the identifications per run from the search's own record.
Hermetic: temp sessions, in-process.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import audit_results  # noqa: E402
import make_analysis_html as mah  # noqa: E402
import sage_lfq_check as slc  # noqa: E402


def write(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh)


class Base(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.s = self._t.name
        self.tables = os.path.join(self.s, "output", "tables")
        os.makedirs(self.tables)

    def tearDown(self):
        self._t.cleanup()

    def notes(self):
        t = mah.Tables(self.tables, {}, 0.05, "default")
        return mah.glance_data({}, t, self.tables, self.s)["notes"]


class IdentificationsOnly(Base):
    def test_the_page_says_it_with_the_runs(self):
        d = os.path.join(self.s, "output", "search_sage")
        write(os.path.join(d, "search_provenance.json"), {"engine": "sage", "quantified": False})
        write(os.path.join(d, slc.RECORD), slc.not_applicable())
        write(os.path.join(d, "sage_identifications.json"),
              {"per_run": {"runA": {"precursors": 812}, "runB": {"precursors": 1790}}})
        note = mah.identifications_only_note({}, self.s)
        self.assertEqual(note["title"], "Identifications only: no quantities")
        self.assertIn("no differential expression", note["text"])
        self.assertIn("runA: 812 precursors; runB: 1,790 precursors", note["text"])
        self.assertIn(note, self.notes())

    def test_a_fragpipe_record_counts_peptides(self):
        d = os.path.join(self.s, "output", "search_fragpipe")
        write(os.path.join(d, "search_provenance.json"), {
            "engine": "fragpipe", "quantification": slc.not_applicable(engine="fragpipe")})
        write(os.path.join(d, "fragpipe_identifications.json"), {"per_run": {"r1": 40}})
        self.assertIn("r1: 40 peptides", mah.identifications_only_note({}, self.s)["text"])

    def test_a_quantifying_search_has_no_such_note(self):
        write(os.path.join(self.s, "output", "search", "search_provenance.json"),
              {"engine": "diann", "version": "2.7.0"})
        self.assertIsNone(mah.identifications_only_note({}, self.s))
        self.assertIsNone(mah.first_pass_note({}, self.s))


class FirstPass(Base):
    def test_the_page_says_it_in_the_audits_plain_words(self):
        write(os.path.join(self.s, "output", "search_dda", "search_provenance.json"),
              {"engine": "diann", "version": "2.7.0", "report_used": "first_pass",
               "dda_first_pass": {"report_used": "first_pass"}})
        note = mah.first_pass_note({}, self.s)
        self.assertEqual(note["kind"], "warning")
        self.assertEqual(note["text"], audit_results.PLAIN["dda_first_pass"][1]["WARN"])
        self.assertIn(note, self.notes())


if __name__ == "__main__":
    unittest.main()
