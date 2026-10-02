#!/usr/bin/env python3
"""
Independent 2.10.0 review (wrong-results risks), qc_bracket.py. Each test below FAILS on
release/skill-2.10.0 @ 40c5b9f and states what should hold instead.

Hermetic: STAN is a fake (qc_bracket.get_json replaced; urllib's urlopen made to fail).
Everything is synthetic. Stdlib only.
"""
import contextlib
import io
import json
import os
import sys
import unittest
import urllib.parse
from datetime import timedelta
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import qc_bracket as qb  # noqa: E402
from test_qc_bracket import (Base, FakeStan, FIRST, LAST, NOW, TIMS, EXPLORIS,  # noqa: E402
                             baseline, entry, files_json, run)


class ReviewBase(Base):
    def session_with(self, rows, entries):
        """A session dir with input/raw_files.txt + acquisition.json, checked in --session mode."""
        s = os.path.join(self.tmp, f"session_{len(os.listdir(self.tmp))}")
        for d in ("input", "logs"):
            os.makedirs(os.path.join(s, d), exist_ok=True)
        self.write_inputs(s, entries)
        with mock.patch.object(qb, "get_json", FakeStan(rows)), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            rc = qb.main(["--session", s, "--stan-url", "http://stan.invalid"], now=NOW)
        return s, rc

    @staticmethod
    def write_inputs(s, entries):
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw files\n" + "".join(e["file"] + "\n" for e in entries))
        files_json(os.path.join(s, "input", "acquisition.json"), entries)


class UnplacedFilesHoldTheGate(ReviewBase):
    def test_files_that_could_not_be_placed_do_not_let_good_through(self):
        """3 timsTOF files are placed and their QC is normal; 4 Exploris .raw could not be read
        (ThermoRawFileParser failed), so the Exploris QC was never looked at. The record says
        'good' and the gate lets the client report out with no staff look -- a partial check
        is a check that did not run for those samples (the gate holds when NONE could be
        placed, but not when most could not)."""
        unread = [dict(entry(f"E{i}.raw", None, None), acquired_at_note="not read: TRFP failed")
                  for i in range(4)]
        s, rc = self.session_with(baseline() + [run(FIRST - timedelta(hours=2))],
                                  self.project() + unread)
        self.assertEqual(rc, qb.EXIT_OK)
        with open(qb.record_paths(s)["json"]) as fh:
            rec = json.load(fh)
        self.assertEqual(len(rec["not_checked"]), 4)
        g = qb.delivery_gate(s)
        self.assertFalse(g["proceed"], "4 of 7 files were never checked, yet the gate passes: "
                                       f"verdict={g['verdict']!r}")


class NoComparisonIsNotNormal(ReviewBase):
    def test_a_qc_run_with_no_baseline_to_compare_is_not_graded_normal(self):
        """Fewer than MIN_BASELINE same-mode QC runs in the 30 days before the project -- a
        sparse Exploris schedule, or the month after a logged downtime -- leaves nothing to
        compare with, and grade_run() then grades the run 'ok' whatever its IDs. A QC run during
        the project at 10% of the instrument's usual IDs reads GOOD ("within this instrument's
        normal range") and the delivery goes out unreviewed."""
        kw = dict(instrument=EXPLORIS, mode="DIA", spd=38)
        cases = {
            # two QC runs in the 30 days before (40,000 each), then 4,000 during the project
            "sparse": [run(FIRST - timedelta(days=d), ids=40000, **kw) for d in (9, 20)],
            # a 5-week downtime: the last normal QC runs are 36-50 days before the project
            "after_downtime": [run(FIRST - timedelta(days=d), ids=40000, **kw)
                               for d in (36, 40, 44, 50)],
        }
        for name, rows in cases.items():
            with self.subTest(name):
                rows = rows + [run(FIRST + timedelta(hours=3), ids=4000, **kw)]
                _rc, rec, _ = self.check(rows, self.project(EXPLORIS))
                during = rec["instruments"][0]["during"][0]
                self.assertIn("no comparison", during["baseline"]["runs"])
                self.assertNotEqual(rec["verdict"], "good",
                                    f"a QC run at 10% of usual IDs, with no baseline to compare "
                                    f"against, reads {rec['verdict']!r}")


class IsolationNeedsNearbyNeighbours(ReviewBase):
    def test_a_flagged_run_is_not_isolated_by_neighbours_weeks_away(self):
        """The only QC run during the project identified 60% of the usual (a check). Its
        neighbours in STAN are 20 days before and 20 days after the project, both normal, so it
        is called 'isolated' (a bad injection, 'the instrument working on both sides of it'),
        counts as ok, and the verdict is GOOD -- the delivery goes out unreviewed, although
        nothing measured the instrument as normal within weeks of the samples."""
        rows = [run(FIRST - timedelta(days=d)) for d in (20, 23, 26, 29)] + [
            run(FIRST + timedelta(hours=3), ids=24000),           # 60% of 40,000: check
            run(LAST + timedelta(days=20))]
        _rc, rec, _ = self.check(rows, self.project())
        bad = rec["instruments"][0]["during"][0]
        self.assertEqual(bad["grade"], "check")
        self.assertNotEqual(rec["verdict"], "good",
                            "a check on the only QC run near the project is discounted as "
                            "'isolated' by neighbours 20 days away")


class MethodsSentenceIsTrue(ReviewBase):
    def test_methods_do_not_claim_qc_on_an_instrument_with_none(self):
        """The client's Methods say QC standards were 'acquired regularly on each instrument
        used' whenever ANY instrument has QC on record. Here half the samples ran on an
        instrument STAN has no QC for (and in the second case on an instrument that could not
        even be read): the sentence is false for those samples."""
        cases = {
            "no_qc_on_record": self.project() + [entry("A1.raw", "Orbitrap Astral", FIRST),
                                                 entry("A2.raw", "Orbitrap Astral", LAST)],
            "not_checked": self.project() + [dict(entry("E1.raw", None, None),
                                                  acquired_at_note="not read: TRFP failed")],
        }
        for name, ents in cases.items():
            with self.subTest(name):
                _rc, rec, _ = self.check(baseline() + [run(FIRST - timedelta(hours=2))], ents)
                s = rec["methods_sentence"] or ""
                self.assertNotIn("each instrument used", s,
                                 f"{name}: Methods claim QC on every instrument used")


class TheRecordMustCoverTheFilesDelivered(ReviewBase):
    def test_a_record_for_another_file_list_does_not_pass_the_gate(self):
        """qc_bracket ran at step 2 on 3 files (good). Later the session's raw list grows (a
        late batch, re-injections kept with --reinjections all, a re-search with more files)
        and one new file ran just after a failed QC injection. The record does not store which
        files it checked (only a count), and delivery_gate() never compares it with the
        session's raw list: the old 'good' releases the new data."""
        rows = baseline() + [run(FIRST - timedelta(hours=2)), run(LAST + timedelta(hours=1)),
                             run(LAST + timedelta(days=2), ids=0)]
        s, rc = self.session_with(rows, self.project())
        self.assertEqual(rc, qb.EXIT_OK)
        self.assertTrue(qb.delivery_gate(s)["proceed"])
        self.write_inputs(s, self.project() + [entry("S04.d", TIMS, LAST + timedelta(days=2,
                                                                                    hours=1))])
        g = qb.delivery_gate(s)
        self.assertFalse(g["proceed"], "the QC record covers 3 files; the session now has 4, "
                                       "and the gate still says the report may go out")


class PagingStopsAtAnEmptyFilteredPage(ReviewBase):
    def test_an_empty_qc_page_is_not_the_end_of_stans_history(self):
        """STAN's /api/runs?qc_only=true (stan/db.py get_runs, db_pg.get_runs_pg) fetches
        limit*3 raw rows at `offset`, keeps the QC-named ones and returns the first `limit`. A
        page is therefore EMPTY whenever 3*limit consecutive rows of the instrument are non-QC
        (legacy baseline rows), not only at the end of the table. stan_runs() stops at the
        first empty page, so every QC run older than such a block -- the 30-day baseline, the
        nearest QC before -- is lost. (Mechanism checked in STAN's source; whether STAN's live
        table holds such a block was NOT checked: no network in this review.)"""
        LIMIT = 10
        qc_new = [run(FIRST + timedelta(days=40) + timedelta(hours=k)) for k in range(5)]
        legacy = [dict(run(FIRST + timedelta(days=39) - timedelta(minutes=k)),
                       run_name=f"Sample_{k:04d}.d") for k in range(3 * LIMIT + 5)]
        qc_old = baseline()                                     # the 30 days before FIRST
        table = sorted(qc_new + legacy + qc_old, key=lambda r: r["run_date"], reverse=True)

        def stan_like(url, timeout):
            q = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
            o, lim = int(q["offset"][0]), int(q["limit"][0])
            raw = table[o:o + 3 * lim]
            return [r for r in raw if "QC" in r["run_name"]][:lim]

        with mock.patch.object(qb, "PAGE", LIMIT), mock.patch.object(qb, "get_json", stan_like):
            got, _pages = qb.stan_runs("http://stan.invalid", TIMS, FIRST - timedelta(days=30),
                                       5, lambda m: None)
        ids = {r["id"] for r in got}
        self.assertTrue({r["id"] for r in qc_old} <= ids,
                        f"{len([r for r in qc_old if r['id'] not in ids])} of {len(qc_old)} "
                        f"baseline QC runs were never fetched: paging stopped at an empty page")


if __name__ == "__main__":
    unittest.main()
