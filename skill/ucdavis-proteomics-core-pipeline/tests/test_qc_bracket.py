#!/usr/bin/env python3
"""
qc_bracket.py: the QC runs around a project, from STAN, and one honest verdict.

What it must never do: call a project "good" when STAN could not be asked; grade on STAN's
`gate_result` (a constant "pass"); read a Thermo creation date as UTC; hide a file it could not
place in time; leave the report silent when the check did not run.

Hermetic: STAN is a fake (qc_bracket.get_json replaced, and urllib's urlopen made to fail, so
nothing here can reach the network). Every QC run, file name and instrument time below is
synthetic. Stdlib only (unittest).
"""
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
import urllib.error
import urllib.parse
from datetime import datetime, timedelta, timezone
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import acq_time  # noqa: E402
import qc_bracket as qb  # noqa: E402
import staff  # noqa: E402

TIMS, LUMOS, EXPLORIS = "timsTOF HT", "Orbitrap Fusion Lumos", "Orbitrap Exploris 480"
# The project: Tuesday 2026-09-15 09:00 to Wednesday 2026-09-16 17:00, Pacific (PDT, UTC-7)
FIRST = datetime(2026, 9, 15, 16, 0, tzinfo=timezone.utc)
LAST = datetime(2026, 9, 17, 0, 0, tzinfo=timezone.utc)
NOW = datetime(2026, 10, 17, tzinfo=timezone.utc)      # every check runs "today" = this


def iso(t):
    return t.isoformat()


def run(when, ids=40000, ips=70, instrument=TIMS, mode="diaPASEF", spd=60, ms1=1.2, ms2=2.5,
        fwhm_min=0.06, name=None, n=[0]):
    """A synthetic STAN /api/runs row, with every field this script reads -- and gate_result
    "pass", always, as STAN records it."""
    n[0] += 1
    return {"id": f"run-{n[0]:05d}", "instrument": instrument,
            "run_name": name or f"QC_HeLa_example_{n[0]:03d}.d", "run_date": iso(when),
            "mode": mode, "spd": spd, "gradient_length_min": 21,
            "n_precursors": ids if "dia" in mode.lower() else None,
            "n_psms": ids if "dda" in mode.lower() else None,
            "ips_score": ips, "median_mass_acc_ms1_ppm": ms1, "median_mass_acc_ms2_ppm": ms2,
            "fwhm_rt_min": fwhm_min, "gate_result": "pass", "failed_gates": "[]"}


def event(etype, date, end=None, instrument=TIMS, first_run=None, n=[0]):
    """A synthetic STAN maintenance event. notes/operator carry canaries: they name people in
    real logs, and nothing may copy them out."""
    n[0] += 1
    return {"id": f"ev-{n[0]:04d}", "instrument": instrument, "event_type": etype,
            "event_date": date, "end_date": end, "notes": "CANARY-NOTES Dr Example said so",
            "operator": "CANARY-OPERATOR", "created_by": "canary@example.edu",
            "column_vendor": None, "column_model": None, "column_serial": None,
            "first_run": first_run, "part_spec": None}


def baseline(instrument=TIMS, days=30, every=2, **kw):
    """Normal QC runs every `every` days over the 30 days before FIRST."""
    return [run(FIRST - timedelta(days=d, hours=3), instrument=instrument, **kw)
            for d in range(every, days, every)]


def files_json(path, entries):
    with open(path, "w") as fh:
        json.dump({"files": entries}, fh)
    return path


def entry(name, instrument, when, note=None):
    return {"file": f"/data/PROT_0000/{name}", "instrument": instrument,
            "acquired_at": acq_time.iso_utc(when) if when else None,
            "acquired_at_source": "synthetic" if when else None,
            "acquired_at_note": note}


class FakeStan:
    """/api/runs as STAN serves it (one instrument, newest first, limit/offset paging) and
    /api/instruments/{instrument}/events (the maintenance log, newest first)."""

    def __init__(self, rows, fail=None, answer=None, events=None, events_fail=None):
        self.rows, self.fail, self.answer, self.calls = rows, fail, answer, []
        self.events, self.events_fail, self.event_calls = events or {}, events_fail, []

    def __call__(self, url, timeout):
        u = urllib.parse.urlparse(url)
        q = urllib.parse.parse_qs(u.query)
        if self.fail:
            raise self.fail
        if u.path.endswith("/events"):
            inst = urllib.parse.unquote(u.path.split("/")[-2])
            self.event_calls.append(inst)
            if self.events_fail:
                raise self.events_fail
            return sorted(self.events.get(inst, []), key=lambda e: e["event_date"], reverse=True)
        self.calls.append(q)
        if self.answer is not None:
            return self.answer
        assert q["qc_only"] == ["true"]
        rows = sorted((r for r in self.rows if r["instrument"] == q["instrument"][0]),
                      key=lambda r: r["run_date"], reverse=True)
        o, lim = int(q["offset"][0]), int(q["limit"][0])
        return rows[o:o + lim]


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="qcb_")
        # no network, whatever happens
        p = mock.patch("urllib.request.urlopen",
                       side_effect=AssertionError("a test reached for the network"))
        p.start()
        self.addCleanup(p.stop)
        # never the Core's real staff list (it exists on HIVE): a missing one unless a test says
        e = mock.patch.dict(os.environ, {staff.STAFF_FILE_ENV: os.path.join(self.tmp, "no_staff.txt")})
        e.start()
        self.addCleanup(e.stop)

    def check(self, rows, entries, fake=None, extra=()):
        fake = fake or FakeStan(rows)
        acq = files_json(os.path.join(self.tmp, "acq.json"), entries)
        out = os.path.join(self.tmp, "qc_bracket.json")
        err = io.StringIO()
        with mock.patch.object(qb, "get_json", fake), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(err):
            rc = qb.main(["--files", *[e["file"] for e in entries], "--acquisition-json", acq,
                          "--out", out, "--stan-url", "http://stan.invalid", *extra], now=NOW)
        with open(out) as fh:
            rec = json.load(fh)
        self.stderr = err.getvalue()
        return rc, rec, fake

    def project(self, instrument=TIMS):
        return [entry("S01.d", instrument, FIRST), entry("S02.d", instrument, FIRST + timedelta(hours=10)),
                entry("S03.d", instrument, LAST)]


class Bracketing(Base):
    def test_qc_before_during_and_after_all_normal_is_good(self):
        rows = baseline() + [run(FIRST - timedelta(hours=2)), run(FIRST + timedelta(hours=20)),
                             run(LAST + timedelta(hours=3))]
        rc, rec, _ = self.check(rows, self.project())
        self.assertEqual(rc, qb.EXIT_OK)
        self.assertEqual((rec["status"], rec["verdict"]), ("ok", "good"))
        o = rec["instruments"][0]
        self.assertEqual(o["before"]["days_from_project"], round(2 / 24, 1))
        self.assertEqual(len(o["during"]), 1)
        self.assertEqual(o["after"]["days_from_project"], round(3 / 24, 1))
        self.assertFalse(o["desert"])
        self.assertIn("GOOD", rec["summary"])
        self.assertIn("1 QC run during the project, all within normal range.", rec["summary"])
        self.assertEqual(rec["methods_sentence"], qb.METHODS_SENTENCE)

    def test_a_qc_desert_is_said_plainly_and_is_a_check(self):
        """Real deserts exist (a Fusion Lumos went 27 days without QC): the nearest QC on each
        side, with how far away, and no 'good'."""
        rows = baseline(LUMOS, mode="DIA", spd=32) + [
            run(FIRST - timedelta(days=12), instrument=LUMOS, mode="DIA", spd=32),
            run(LAST + timedelta(days=15), instrument=LUMOS, mode="DIA", spd=32)]
        # the baseline's newest run is 2 days before FIRST: drop those within 12 days
        rows = [r for r in rows if acq_time.parse_iso(r["run_date"]) <= FIRST - timedelta(days=12)
                or acq_time.parse_iso(r["run_date"]) > LAST]
        rc, rec, _ = self.check(rows, self.project(LUMOS))
        o = rec["instruments"][0]
        self.assertEqual(rc, qb.EXIT_OK)
        self.assertTrue(o["desert"])
        self.assertEqual(o["verdict"], "check")
        self.assertEqual(o["during"], [])
        self.assertEqual(round(o["before"]["days_from_project"]), 12)
        self.assertEqual(round(o["after"]["days_from_project"]), 15)
        self.assertIn("No QC run during the project.", rec["summary"])
        self.assertIn("No QC run within 7 days of the project", rec["summary"])
        self.assertIn("12 days before the first sample", rec["summary"])

    def test_no_qc_after_yet(self):
        rows = baseline() + [run(FIRST - timedelta(hours=2))]
        rc, rec, _ = self.check(rows, self.project())
        o = rec["instruments"][0]
        self.assertIsNone(o["after"])
        self.assertTrue(o["after_note"].startswith("none yet"))
        self.assertIn("Nearest QC after: none yet", rec["summary"])
        self.assertEqual(rec["verdict"], "good")        # the run before was near and normal

    def test_a_zero_id_qc_run_is_a_concern_and_names_the_samples_around_it(self):
        rows = baseline() + [run(FIRST - timedelta(hours=2)),
                             run(FIRST + timedelta(hours=5), ids=0, ips=0)]
        rc, rec, _ = self.check(rows, self.project())
        o = rec["instruments"][0]
        bad = o["during"][0]
        self.assertEqual(bad["grade"], "concern")
        self.assertEqual(bad["flags"][0]["code"], "zero_ids")
        self.assertFalse(bad["isolated"])               # no normal QC after it
        self.assertEqual(rec["verdict"], "concern")
        # the samples a staff member should look at: after the last normal QC, no normal one since
        self.assertEqual(bad["samples_between_normal_qc"]["n_files"], 3)
        self.assertEqual(bad["samples_between_normal_qc"]["to"], "no normal QC run after it yet")

    def test_an_isolated_failed_injection_counts_one_level_lower(self):
        rows = baseline() + [run(FIRST - timedelta(hours=2)),
                             run(FIRST + timedelta(hours=5), ids=0, ips=0),
                             run(FIRST + timedelta(hours=8)), run(LAST + timedelta(hours=2))]
        rc, rec, _ = self.check(rows, self.project())
        o = rec["instruments"][0]
        bad = [x for x in o["during"] if x["grade"] != "ok"][0]
        self.assertTrue(bad["isolated"])
        self.assertEqual(bad["counts_as"], "check")
        self.assertEqual(rec["verdict"], "check")
        self.assertEqual(bad["samples_between_normal_qc"]["n_files"], 1)   # only S01 ran between
        self.assertIn("Isolated", rec["summary"])

    def test_ids_below_80_and_50_percent_of_the_30_day_median(self):
        rows = baseline() + [run(FIRST + timedelta(hours=1), ids=30000),    # 75%: check
                             run(FIRST + timedelta(hours=2), ids=19000)]    # 47.5%: concern
        rc, rec, _ = self.check(rows, self.project())
        d = rec["instruments"][0]["during"]
        self.assertEqual([x["grade"] for x in d], ["check", "concern"])
        self.assertEqual([x["ids_vs_median"] for x in d], [0.75, 0.475])

    def test_mass_error_and_peak_width_against_the_baseline(self):
        rows = baseline() + [run(FIRST + timedelta(hours=1), ms1=2.0),     # +0.8 ppm: under floor
                             run(FIRST + timedelta(hours=2), ms1=2.5),     # 2.1x, +1.3: check
                             run(FIRST + timedelta(hours=3), fwhm_min=0.1),  # 1.67x: check
                             run(FIRST + timedelta(hours=4))]
        rc, rec, _ = self.check(rows, self.project())
        d = rec["instruments"][0]["during"]
        self.assertEqual([x["grade"] for x in d], ["ok", "check", "check", "ok"])
        self.assertEqual(d[1]["flags"][-1]["code"], "ms1_ppm_high")
        self.assertEqual(d[2]["flags"][-1]["code"], "fwhm_wide")

    def test_baseline_is_same_throughput_unless_too_few(self):
        rows = baseline(spd=60) + baseline(spd=100, ids=30000, every=5) + [
            run(FIRST + timedelta(hours=1), spd=100, ids=29000),
            run(FIRST + timedelta(hours=2), spd=30, ids=40000)]
        rc, rec, _ = self.check(rows, self.project())
        d = rec["instruments"][0]["during"]
        self.assertIn("at 100 samples/day", d[0]["baseline"]["runs"])
        self.assertEqual(d[0]["grade"], "ok")         # 97% of the 100-spd median, not 72% of all
        self.assertIn("all throughputs", d[1]["baseline"]["runs"])


class IpsIsShownNotGraded(Base):
    def test_bands_match_stans_dashboard(self):
        self.assertEqual([qb.ips_band(x) for x in (100, 80, 79, 60, 59, 0, None)],
                         ["green", "green", "amber", "amber", "red", "red", None])

    def test_red_amber_and_green_ips_with_usual_ids_are_all_good(self):
        """IPS compares with STAN's reference cohort, not the instrument's recent self: an
        instrument whose QC usually scores red is not a concern for scoring red again."""
        rows = baseline(ips=25) + [run(FIRST + timedelta(hours=1), ips=25),
                                   run(FIRST + timedelta(hours=2), ips=65),
                                   run(FIRST + timedelta(hours=3), ips=85)]
        rc, rec, _ = self.check(rows, self.project())
        d = rec["instruments"][0]["during"]
        self.assertEqual([x["ips_band"] for x in d], ["red", "amber", "green"])
        self.assertEqual([x["grade"] for x in d], ["ok", "ok", "ok"])
        self.assertEqual(rec["verdict"], "good")
        self.assertIn("'something is wrong' level", d[0]["flags"][0]["text"])  # still named
        self.assertIn("median IPS of 25 (red)", rec["summary"])

    def test_gate_result_is_never_read(self):
        import inspect
        self.assertNotIn('"gate_result"', inspect.getsource(qb.grade_run))
        self.assertNotIn("'gate_result'", inspect.getsource(qb))


class InstrumentsAndFiles(Base):
    def test_unknown_instrument_has_no_qc_on_record(self):
        rc, rec, _ = self.check([], [entry("A.raw", "Orbitrap Astral", FIRST)])
        self.assertEqual(rc, qb.EXIT_OK)
        self.assertEqual(rec["verdict"], "no_qc_on_record")
        self.assertIn("STAN has no QC run on record", rec["summary"])

    def test_mixed_instruments_get_one_bracket_each_and_the_worst_verdict(self):
        rows = (baseline() + [run(FIRST - timedelta(hours=2))]
                + baseline(EXPLORIS, mode="DIA", spd=38)
                + [run(FIRST + timedelta(hours=1), instrument=EXPLORIS, mode="DIA", spd=38, ids=0)])
        ents = self.project() + [entry("E01.raw", EXPLORIS, FIRST + timedelta(hours=2))]
        rc, rec, fake = self.check(rows, ents)
        by = {o["instrument"]: o for o in rec["instruments"]}
        self.assertEqual(set(by), {TIMS, EXPLORIS})
        self.assertEqual(by[TIMS]["verdict"], "good")
        self.assertEqual(by[EXPLORIS]["verdict"], "concern")
        self.assertEqual(rec["verdict"], "concern")
        self.assertEqual({q["instrument"][0] for q in fake.calls}, {TIMS, EXPLORIS})

    def test_no_qc_on_record_ranks_above_check_and_below_concern(self):
        self.assertEqual(max(["good", "check", "no_qc_on_record"], key=qb.RANK.get),
                         "no_qc_on_record")
        self.assertEqual(max(["no_qc_on_record", "concern"], key=qb.RANK.get), "concern")

    def test_a_file_with_no_readable_time_is_listed_not_guessed(self):
        bad = acq_time.from_thermo_value("30.09.2026 03:06:18")      # a day-first culture
        self.assertIsNone(bad["acquired_at"])
        rows = baseline() + [run(FIRST - timedelta(hours=2))]
        ents = self.project() + [dict(entry("X.raw", EXPLORIS, None),
                                      acquired_at_note=bad["acquired_at_note"])]
        rc, rec, _ = self.check(rows, ents)
        self.assertEqual(rc, qb.EXIT_OK)
        self.assertEqual(rec["files"]["n_checked"], 3)
        self.assertEqual(rec["not_checked"][0]["file"], "/data/PROT_0000/X.raw")
        self.assertIn("no format this reads", rec["not_checked"][0]["why"])
        self.assertIn("1 file could not be placed in time", rec["summary"])

    def test_a_time_with_no_zone_or_no_sense_is_not_checked(self):
        rows = baseline() + [run(FIRST - timedelta(hours=2))]
        ents = self.project() + [dict(entry("N.raw", EXPLORIS, None), acquired_at="2026-09-15T10:00:00"),
                                 dict(entry("G.raw", EXPLORIS, None), acquired_at="soon")]
        rc, rec, _ = self.check(rows, ents)
        self.assertEqual(rc, qb.EXIT_OK)
        self.assertEqual(sorted(os.path.basename(f["file"]) for f in rec["not_checked"]),
                         ["G.raw", "N.raw"])
        self.assertTrue(all("with a time zone" in f["why"] for f in rec["not_checked"]))

    def test_an_unreadable_step2_json_is_said_and_the_files_are_read(self):
        bad = os.path.join(self.tmp, "broken.json")
        with open(bad, "w") as fh:
            fh.write("{")
        err = io.StringIO()
        with mock.patch.object(qb, "get_json", FakeStan([])), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(err):
            rc = qb.main(["--files", "/data/PROT_0000/none.d", "--acquisition-json", bad,
                          "--out", os.path.join(self.tmp, "o.json")], now=NOW)
        self.assertEqual(rc, qb.EXIT_NOTHING)
        self.assertIn("WARNING", err.getvalue())
        with open(os.path.join(self.tmp, "o.json")) as fh:
            rec = json.load(fh)
        self.assertIn("could not be read", rec["files"]["acquisition_json_problem"])
        self.assertEqual(rec["not_checked"][0]["why"], "not found on this computer")

    def test_nothing_placeable_exits_4_with_a_record(self):
        rc, rec, fake = self.check([], [entry("X.raw", EXPLORIS, None, note="not read: synthetic")])
        self.assertEqual(rc, qb.EXIT_NOTHING)
        self.assertEqual(rec["status"], "nothing_checked")
        self.assertIsNone(rec["verdict"])
        self.assertEqual(fake.calls, [])              # STAN is not asked about nothing

    def test_blank_runs_and_impossible_dates_are_left_out(self):
        rows = baseline() + [run(FIRST + timedelta(hours=1), ids=0, name="blankDia_before_HeLa.d"),
                             run(datetime(1980, 1, 2, tzinfo=timezone.utc)),
                             run(FIRST - timedelta(hours=2))]
        rc, rec, _ = self.check(rows, self.project())
        o = rec["instruments"][0]
        self.assertEqual(o["during"], [])
        self.assertEqual(sorted(x["why"] for x in o["ignored_rows"]),
                         ["named as a blank: not a QC standard",
                          "run_date is not a plausible acquisition time"])
        self.assertEqual(rec["verdict"], "good")


class MaintenanceLog(Base):
    """Long QC gaps were often the instrument being down, not QC missed. STAN's maintenance log
    explains a gap when it has an event there; when it has none, the report says so plainly and
    guesses no cause."""

    def lumos_desert(self, events=None, events_fail=None):
        """A Fusion Lumos project 2026-09-15..16 with no QC for 27 days around it."""
        kw = dict(instrument=LUMOS, mode="DIA", spd=32)
        rows = ([run(FIRST - timedelta(days=d, hours=3), **kw) for d in range(12, 40, 2)]
                + [run(LAST + timedelta(days=14), **kw)])
        fake = FakeStan(rows, events={LUMOS: events or []}, events_fail=events_fail)
        return self.check(rows, self.project(LUMOS), fake=fake)

    def test_a_gap_explained_by_a_logged_downtime(self):
        rc, rec, fake = self.lumos_desert([event("downtime", "2026-09-04", "2026-09-29",
                                                 instrument=LUMOS)])
        o = rec["instruments"][0]
        self.assertEqual(fake.event_calls, [LUMOS])
        self.assertEqual(len(o["gaps"]), 1)
        g = o["gaps"][0]
        self.assertEqual(g["days"], 27)
        self.assertEqual(g["events"], ["instrument downtime 2026-09-04 to 2026-09-29"])
        self.assertIn("STAN's maintenance log has, in this period: instrument downtime "
                      "2026-09-04 to 2026-09-29.", rec["summary"])
        # samples acquired inside a logged downtime: a contradiction worth a look
        m = o["maintenance"]
        self.assertEqual(m["files_during_downtime"], 3)
        self.assertIn("acquired while STAN's log has the instrument down", rec["summary"])
        self.assertEqual(rec["verdict"], "check")

    def test_a_gap_with_no_event_is_said_plainly_with_no_cause(self):
        rc, rec, _ = self.lumos_desert([event("column_change", "2026-01-10", instrument=LUMOS)])
        g = rec["instruments"][0]["gaps"][0]
        self.assertEqual(g["events"], [])
        self.assertTrue(g["text"].startswith("No QC for 27 days (2026-09-03 to 2026-09-30)"))
        self.assertTrue(g["text"].endswith("STAN has no maintenance record for this period."))
        self.assertIn("nothing for this instrument within 30 days", rec["summary"])

    def test_an_instrument_with_an_empty_log_says_so(self):
        rc, rec, _ = self.lumos_desert([])
        self.assertIn("(none is logged for this instrument at all)",
                      rec["instruments"][0]["gaps"][0]["text"])
        self.assertIn("STAN has no maintenance events logged for this instrument.", rec["summary"])

    def test_an_unreadable_log_leaves_the_verdict_and_says_the_gap_is_unexplained(self):
        rc, rec, _ = self.lumos_desert(events_fail=urllib.error.HTTPError("u", 404, "nf", {}, None))
        self.assertEqual(rc, qb.EXIT_OK)
        self.assertEqual(rec["verdict"], "check")       # the desert, from the QC runs alone
        self.assertEqual(rec["instruments"][0]["maintenance"]["log"], "unreadable")
        self.assertIn("could not be read (http://stan.invalid answered HTTP 404",
                      rec["instruments"][0]["gaps"][0]["text"])
        self.assertIn("WARNING", self.stderr)

    def test_samples_straddling_a_column_change_and_the_post_repair_qc(self):
        """Column changed on 2026-09-15 (a date only, as STAN's Trends form stores one): the
        project's first file ran that morning, the others after. The first QC on or after that
        day is the check for them and speaks for the project; files after the change and
        before any QC are a check."""
        rows = baseline() + [run(FIRST - timedelta(days=3)),              # before the change
                             run(LAST + timedelta(hours=3))]               # the first QC after
        ents = [entry("S01.d", TIMS, FIRST),                               # 09:00 on the day
                entry("S02.d", TIMS, FIRST + timedelta(hours=26)),         # the next day
                entry("S03.d", TIMS, LAST)]
        fake = FakeStan(rows, events={TIMS: [event("column_change", "2026-09-15T12:00:00Z")]})
        rc, rec, _ = self.check(rows, ents, fake=fake)
        o = rec["instruments"][0]
        ev = o["maintenance"]["events"][0]
        self.assertEqual(ev["when"], "2026-09-15 (date only)")
        self.assertTrue(ev["date_only"])
        self.assertEqual(ev["files_same_day"], ["S01.d"])
        self.assertEqual(ev["files_after"], 2)
        self.assertEqual(ev["files_unchecked_after"], ["S02.d", "S03.d"])
        self.assertEqual(ev["first_qc_after"]["stan_id"], o["after"]["stan_id"])
        self.assertEqual(o["after"]["after_maintenance"], "column change")
        self.assertIn("whether they ran before or after it is not known", rec["summary"])
        self.assertIn("acquired after the column change and before the first QC run after it",
                      rec["summary"])
        self.assertEqual(rec["verdict"], "check")

    def test_first_run_anchors_the_change_and_a_failed_post_repair_qc_is_never_discounted(self):
        after_change = run(FIRST + timedelta(hours=2), ids=0, ips=0,
                           name="QC_HeLa_example_newcolumn.d")
        rows = baseline() + [run(FIRST - timedelta(hours=1)), after_change,
                             run(FIRST + timedelta(hours=4)), run(LAST + timedelta(hours=2))]
        fake = FakeStan(rows, events={TIMS: [event("column_change", "2026-09-15",
                                                   first_run="QC_HeLa_example_newcolumn")]})
        ents = [entry("S01.d", TIMS, FIRST),                       # on the old column
                entry("S02.d", TIMS, FIRST + timedelta(hours=3)),  # after it, before the next QC
                entry("S03.d", TIMS, LAST)]
        rc, rec, _ = self.check(rows, ents, fake=fake)
        o = rec["instruments"][0]
        ev = o["maintenance"]["events"][0]
        self.assertEqual(ev["anchored_by"], "a QC run in STAN")
        self.assertFalse(ev["date_only"])
        self.assertEqual(ev["files_same_day"], [])
        self.assertEqual(ev["files_after"], 2)               # S01 ran before the new column
        self.assertEqual(ev["files_checked_by_first_qc"], 1)  # S02: no QC between it and the check
        bad = [x for x in o["during"] if x["stan_id"] == after_change["id"]][0]
        self.assertEqual(bad["after_maintenance"], "column change")
        self.assertTrue(bad["isolated"])                    # normal QC on both sides ...
        self.assertEqual(bad["counts_as"], "concern")       # ... but it is the post-repair check
        self.assertEqual(rec["verdict"], "concern")

    def test_a_far_post_repair_qc_still_speaks_and_is_listed(self):
        rows = [r for r in baseline()                       # nothing in the 9 days before
                if acq_time.parse_iso(r["run_date"]) < FIRST - timedelta(days=10)]
        rows.append(run(FIRST - timedelta(days=9), ids=0, ips=0))   # the first after cleaning
        fake = FakeStan(rows, events={TIMS: [event("source_clean", "2026-09-05T18:30:00-07:00")]})
        rc, rec, _ = self.check(rows, self.project(), fake=fake)
        o = rec["instruments"][0]
        self.assertTrue(o["desert"])                        # nothing within 7 days ...
        self.assertEqual(rec["verdict"], "concern")         # ... the post-cleaning QC failed
        self.assertIn("First QC run after the ion source cleaning", rec["summary"])
        ev = o["maintenance"]["events"][0]
        self.assertEqual(ev["files_checked_by_first_qc"], 3)

    def test_a_change_long_before_with_qc_since_does_not_speak_for_the_project(self):
        """A column changed 20 days before, with normal QC every other day since: its failed
        first QC run is old news, not the check for these samples."""
        rows = baseline() + [run(FIRST - timedelta(days=20) + timedelta(hours=5), ids=0, ips=0),
                             run(FIRST - timedelta(hours=2))]
        fake = FakeStan(rows, events={TIMS: [event("column_change", "2026-08-26")]})
        rc, rec, _ = self.check(rows, self.project(), fake=fake)
        ev = rec["instruments"][0]["maintenance"]["events"][0]
        self.assertEqual(ev["files_after"], 3)
        self.assertEqual(ev["files_checked_by_first_qc"], 0)
        self.assertEqual(rec["verdict"], "good")
        self.assertIn("column change 2026-08-26 (date only)", rec["summary"])
        self.assertNotIn("First QC run", rec["summary"])

    def test_notes_and_people_are_never_copied_out(self):
        rows = baseline() + [run(FIRST - timedelta(hours=2))]
        fake = FakeStan(rows, events={TIMS: [event("calibration", "2026-09-14")]})
        rc, rec, _ = self.check(rows, self.project(), fake=fake)
        text = json.dumps(rec)
        self.assertIn("mass calibration 2026-09-14", text)
        for canary in ("CANARY", "Dr Example", "canary@example.edu"):
            self.assertNotIn(canary, text)

    def test_event_times(self):
        day, d_only = qb._event_time("2026-07-30T12:00:00Z")
        self.assertTrue(d_only)
        self.assertEqual(acq_time.iso_utc(day), "2026-07-30T07:00:00Z")      # PDT midnight
        end, _ = qb._event_time("2026-07-30", end=True)
        self.assertEqual(acq_time.iso_utc(end), "2026-07-31T07:00:00Z")
        inst, d_only = qb._event_time("2026-09-16T19:09:26Z")
        self.assertFalse(d_only)
        local, _ = qb._event_time("2026-01-10T09:00:00")                    # no offset: Pacific
        self.assertEqual(acq_time.iso_utc(local), "2026-01-10T17:00:00Z")
        self.assertEqual(qb._event_time("someday"), (None, None))


class StanUnreachable(Base):
    def test_unreachable_is_exit_3_and_never_a_verdict(self):
        fake = FakeStan([], fail=urllib.error.URLError("Name or service not known"))
        rc, rec, _ = self.check([], self.project(), fake=fake)
        self.assertEqual(rc, qb.EXIT_STAN)
        self.assertEqual(rec["status"], "stan_unreachable")
        self.assertIsNone(rec["verdict"])
        self.assertIn("could NOT be checked", rec["summary"])
        self.assertIn("not a pass", rec["summary"])
        self.assertNotIn("GOOD", rec["summary"])

    def test_an_answer_that_is_not_a_run_list_is_unreachable_too(self):
        for answer in ({"detail": "maintenance"}, ["not a run"]):
            rc, rec, _ = self.check([], self.project(), fake=FakeStan([], answer=answer))
            self.assertEqual(rc, qb.EXIT_STAN, answer)
            self.assertIsNone(rec["verdict"])

    def test_http_error_is_unreachable(self):
        fake = FakeStan([], fail=urllib.error.HTTPError("u", 502, "Bad Gateway", {}, None))
        rc, rec, _ = self.check([], self.project(), fake=fake)
        self.assertEqual(rc, qb.EXIT_STAN)
        self.assertIn("HTTP 502", rec["stan"]["error"])


class Paging(Base):
    def test_pages_step_by_page_size_dedupe_and_stop_past_the_baseline(self):
        """With qc_only, STAN filters each page after fetching it, so a page can repeat rows of
        the next one: offsets step by the page size, never by the rows returned, and every
        run is kept once. Paging stops at the first page older than the baseline."""
        rows = [run(FIRST + timedelta(days=60) - timedelta(hours=6 * k)) for k in range(900)]
        plain = FakeStan(rows)

        def overlapping(url, timeout):
            page = plain(url, timeout)
            o = int(plain.calls[-1]["offset"][0])
            nxt = sorted(rows, key=lambda r: r["run_date"], reverse=True)[o + 100:o + 103]
            return page + nxt if page else page
        with mock.patch.object(qb, "PAGE", 100), mock.patch.object(qb, "get_json", overlapping):
            got, pages = qb.stan_runs("http://stan.invalid", TIMS, FIRST - timedelta(days=30), 5,
                                      lambda m: None)
        # 900 rows, 4 a day, the newest 60 days after FIRST: FIRST - 30 days is row 360, so
        # the 4th page of 100 (rows 300-399) is the first to reach past it
        self.assertEqual([int(q["offset"][0]) for q in plain.calls], [0, 100, 200, 300])
        self.assertEqual(pages, 4)
        self.assertEqual(len({r["id"] for r in got}), len(got))
        self.assertLess(min(acq_time.parse_iso(r["run_date"]) for r in got),
                        FIRST - timedelta(days=30))


class TimeZones(Base):
    def test_bruker_offset_and_thermo_wall_clock(self):
        self.assertEqual(acq_time.from_bruker_value("2026-09-30T13:04:03.063-07:00")["acquired_at"],
                         "2026-09-30T20:04:03Z")
        # TRFP on HIVE (invariant culture) and on a US-English Windows: the same instant
        for text in ("09/30/2026 03:06:18", "9/30/2026 3:06:18 AM"):
            r = acq_time.from_thermo_value(text)
            self.assertEqual(r["acquired_at"], "2026-09-30T10:06:18Z", text)     # PDT, UTC-7
            self.assertIn("America/Los_Angeles", r["acquired_at_note"])
        self.assertEqual(acq_time.from_thermo_value("01/15/2026 03:06:18")["acquired_at"],
                         "2026-01-15T11:06:18Z")                                 # PST, UTC-8
        r = acq_time.from_bruker_value("2026-09-30T13:04:03")                    # no offset
        self.assertEqual(r["acquired_at"], "2026-09-30T20:04:03Z")
        self.assertIn("assumption", r["acquired_at_note"])

    def test_iso_forms_python_39_cannot_parse_alone(self):
        for text, want in (("2026-09-29T15:32:57.1234567-07:00", "2026-09-29T22:32:57Z"),
                           ("2026-09-30T12:06:20Z", "2026-09-30T12:06:20Z"),
                           ("2026-09-29T23:15:50.077161+00:00", "2026-09-29T23:15:50Z"),
                           ("2026-09-29 23:15:50+0000", "2026-09-29T23:15:50Z")):
            self.assertEqual(acq_time.iso_utc(acq_time.parse_iso(text)), want, text)
        self.assertIsNone(acq_time.parse_iso("yesterday"))

    def test_the_built_in_pacific_rule_matches_zoneinfo(self):
        """Windows Pythons often have no time-zone database; the stand-in must agree with it
        wherever it exists (the repeated November hour aside: both read its first occurrence)."""
        try:
            from zoneinfo import ZoneInfo
            z = ZoneInfo("America/Los_Angeles")
        except Exception:
            self.skipTest("no zoneinfo database here to compare with")
        p = acq_time._USPacific()
        t = datetime(2024, 1, 1, tzinfo=timezone.utc)
        while t.year < 2028:
            self.assertEqual(t.astimezone(p).replace(tzinfo=None),
                             t.astimezone(z).replace(tzinfo=None), t)
            t += timedelta(hours=7)
        for wall in (datetime(2026, 3, 8, 1, 59), datetime(2026, 3, 8, 3, 0),
                     datetime(2026, 11, 1, 0, 59), datetime(2026, 11, 1, 2, 0)):
            self.assertEqual(wall.replace(tzinfo=p).utcoffset(), wall.replace(tzinfo=z).utcoffset())

    def test_without_zoneinfo_core_tz_still_reads_and_says_how(self):
        with mock.patch.object(acq_time, "_zone", side_effect=lambda name: (
                (acq_time._USPacific(), "built-in US Pacific rule (no time-zone database on this "
                 "computer)") if name == acq_time.CORE_TZ else (None, "no database"))):
            r = acq_time.from_thermo_value("09/30/2026 03:06:18")
            self.assertEqual(r["acquired_at"], "2026-09-30T10:06:18Z")
            self.assertIn("built-in US Pacific rule", r["acquired_at_note"])
            self.assertIsNone(acq_time.from_thermo_value("09/30/2026 03:06:18",
                                                         "Europe/Berlin")["acquired_at"])

    def test_a_late_evening_qc_is_before_an_early_morning_thermo_sample(self):
        """23:30 PDT on the 14th is 06:30 UTC on the 15th: a sample created at 00:30 on the
        15th (PDT, from the .raw) comes after it. Read as UTC, the QC would land 'after'."""
        sample = acq_time.parse_iso(acq_time.from_thermo_value("09/15/2026 00:30:00")["acquired_at"])
        qc = run(datetime(2026, 9, 15, 6, 30, tzinfo=timezone.utc), instrument=LUMOS,
                 mode="DIA", spd=32)
        rows = baseline(LUMOS, mode="DIA", spd=32) + [qc]
        rc, rec, _ = self.check(rows, [entry("L1.raw", LUMOS, sample)])
        o = rec["instruments"][0]
        self.assertEqual(o["before"]["stan_id"], qc["id"])
        self.assertEqual(o["before"]["local_time"], "2026-09-14 23:30 PDT")
        self.assertEqual(o["window"]["first_local"], "2026-09-15 00:30 PDT")

    def test_an_acquisition_time_after_the_file_was_written_is_flagged(self):
        f = os.path.join(self.tmp, "x.raw")
        open(f, "w").close()
        os.utime(f, (datetime(2026, 9, 1, tzinfo=timezone.utc).timestamp(),) * 2)
        rec = acq_time.from_thermo_value("09/05/2026 10:00:00")      # 4 days after the mtime
        self.assertIn("after the file was last written", acq_time.check_against_mtime(rec, f))
        self.assertIsNone(acq_time.check_against_mtime(
            acq_time.from_thermo_value("08/31/2026 10:00:00"), f))


class SessionAndReport(Base):
    def session(self, core=False):
        s = os.path.join(self.tmp, "2026-09-17_Example")
        for d in ("input", "output/tables", "output/figures"):
            os.makedirs(os.path.join(s, d), exist_ok=True)
        ents = self.project()
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw files\n" + "".join(e["file"] + "\n" for e in ents))
        files_json(os.path.join(s, "input", "acquisition.json"), ents)
        with open(os.path.join(s, "output", "AI_Analysis_Report.md"), "w") as fh:
            fh.write("# Example study\n\nIntro.\n\n## Overview\n\nText.\n")
        if core:
            with open(os.path.join(s, "session.json"), "w") as fh:
                json.dump({"coreomics": {"internal_id": "PROT_9999"}}, fh)   # 0000 is not a valid number
        return s

    def render(self, s, extra=()):
        out = os.path.join(s, "output", "Analysis_Report.html")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                            "--session", s, "--out", out, "--no-pdf", *extra],
                           capture_output=True, text=True, env=dict(
                               os.environ, PATH=os.environ.get("PATH", "")))
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out, encoding="utf-8") as fh, \
                open(os.path.splitext(out)[0] + ".md", encoding="utf-8") as fm:
            return fh.read(), fm.read()

    def check_session(self, rows, fake=None):
        s = self.session(core=True)
        with mock.patch.object(qb, "get_json", fake or FakeStan(rows)), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            rc = qb.main(["--session", s, "--stan-url", "http://stan.invalid"], now=NOW)
        return s, rc

    def test_session_mode_writes_the_staff_record_into_logs(self):
        s, rc = self.check_session(baseline() + [run(FIRST - timedelta(hours=2))])
        self.assertEqual(rc, qb.EXIT_OK)
        p = qb.record_paths(s)
        with open(p["json"]) as fh:
            rec = json.load(fh)
        self.assertEqual((rec["schema"], rec["schema_version"], rec["staff_only"]),
                         ("qc_bracket/1", 1, True))
        self.assertEqual(rec["files"]["from"], "input/raw_files.txt")
        self.assertTrue(rec["files"]["acquisition_json"].endswith("acquisition.json"))
        self.assertEqual(rec["verdict"], "good")
        with open(p["md"]) as fh:
            self.assertIn("(staff only -- never delivered)", fh.read())
        for f in (p["json"], p["md"]):
            self.assertEqual(os.stat(f).st_mode & 0o777, 0o640, f)
        self.assertFalse(os.path.exists(os.path.join(s, "output", "qc_bracket.json")))

    # -- Brett, 2026-10-01: the verdict is STAFF-ONLY -------------------------------------
    VERDICT_WORDS = ("Verdict", "CONCERN", "CHECK.", "NO QC ON RECORD", "GOOD.", "QC runs around",
                     "qc_bracket", "IPS", "Not checked")

    def test_the_client_report_carries_no_qc_section_or_verdict_word(self):
        """Whatever the verdict, and whether or not the check could run, the client's report and
        its .md twin say nothing about it."""
        cases = {"concern": baseline() + [run(FIRST - timedelta(hours=2)),
                                          run(FIRST + timedelta(hours=1), ids=0)],
                 "good": baseline() + [run(FIRST - timedelta(hours=2))]}
        for verdict, rows in cases.items():
            s, _rc = self.check_session(rows)
            with open(qb.record_paths(s)["json"]) as fh:
                self.assertEqual(json.load(fh)["verdict"], verdict)
            html_doc, md = self.render(s)
            for doc in (html_doc, md):
                for w in self.VERDICT_WORDS:
                    self.assertNotIn(w, doc, (verdict, w))
            os.rename(s, s + "_" + verdict)
        s, rc = self.check_session([], fake=FakeStan([], fail=urllib.error.URLError("timed out")))
        self.assertEqual(rc, qb.EXIT_STAN)
        for doc in self.render(s):
            for w in self.VERDICT_WORDS:
                self.assertNotIn(w, doc, ("unreachable", w))

    def test_the_methods_sentence_is_identical_across_verdicts(self):
        import make_methods
        rows = {"good": baseline() + [run(FIRST - timedelta(hours=2))],
                "check": baseline() + [run(FIRST + timedelta(hours=1), ids=30000)],
                "concern": baseline() + [run(FIRST + timedelta(hours=1), ids=0)]}
        lines = {}
        for verdict, r in rows.items():
            rc, rec, _ = self.check(r, self.project())
            self.assertEqual(rec["verdict"], verdict)
            lines[verdict] = make_methods.qc_lines(rec, "qc_bracket.json")
        self.assertEqual(lines["good"], lines["check"])
        self.assertEqual(lines["good"], lines["concern"])
        self.assertEqual(lines["good"], ["## Instrument performance", "", qb.METHODS_SENTENCE, ""])
        for w in ("normal", "within", "verdict", "concern", "check"):
            self.assertNotIn(w, qb.METHODS_SENTENCE.lower())
        # no QC on record, a check that could not run, an unreadable record: the SAME sentence --
        # one that came and went with the outcome would tell the client (2.10 safety review)
        _rc, rec, _ = self.check([], [entry("A.raw", "Orbitrap Astral", FIRST)])
        self.assertEqual(make_methods.qc_lines(rec, "x"), lines["good"])
        _rc, rec, _ = self.check([], self.project(),
                                 fake=FakeStan([], fail=urllib.error.URLError("down")))
        self.assertEqual(rec["status"], "stan_unreachable")
        self.assertEqual(make_methods.qc_lines(rec, "x"), lines["good"])
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(make_methods.qc_lines(None, "x.json"), lines["good"])

    def test_provenance_points_to_the_record_without_its_verdict(self):
        _rc, rec, _ = self.check(baseline() + [run(FIRST + timedelta(hours=1), ids=0)],
                                 self.project())
        path = os.path.join(self.tmp, "qc_bracket.json")
        out = os.path.join(self.tmp, "repro")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir", out,
                            "--qc-bracket", path, "--setup-json", os.path.join(self.tmp, "none.json")],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out, "run_manifest.json")) as fh:
            man = json.load(fh)
        self.assertEqual(man["instrument_qc"]["sha256"], qb._sha256(path))
        self.assertNotIn("verdict", man["instrument_qc"])
        self.assertFalse(os.path.exists(os.path.join(out, "inputs", "qc_bracket.json")))
        for root, _d, files in os.walk(out):           # the bundle is delivered: no verdict word
            for f in files:
                with open(os.path.join(root, f), encoding="utf-8", errors="replace") as fh:
                    text = fh.read()
                for w in ("CONCERN", "concern", "Verdict"):
                    self.assertNotIn(w, text, f)


class DeliveryGate(Base):
    """A check or concern holds the client delivery until a staff member records who, when and
    a one-line note, for THIS record; good and no-QC-on-record go through."""

    def session_with(self, rows=None, entries=None, fake=None):
        s = os.path.join(self.tmp, f"session_{len(os.listdir(self.tmp))}")
        os.makedirs(os.path.join(s, "logs"))
        os.makedirs(os.path.join(s, "input"))
        entries = entries or self.project()
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:   # the list it covers
            fh.write("".join(e["file"] + "\n" for e in entries))
        _rc, _rec, _ = self.check(rows if rows is not None else [], entries, fake=fake)
        os.replace(os.path.join(self.tmp, "qc_bracket.json"), qb.record_paths(s)["json"])
        os.replace(os.path.join(self.tmp, "qc_bracket.md"), qb.record_paths(s)["md"])
        return s

    def test_check_and_concern_block_until_acknowledged(self):
        for ids, verdict in ((30000, "check"), (0, "concern")):
            s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=ids)])
            g = qb.delivery_gate(s)
            self.assertEqual((g["verdict"], g["proceed"], g["needs_ack"]), (verdict, False, True))
            self.assertIn(qb.LABEL[verdict], g["reason"])
            g = qb.acknowledge(s, "staffexample", "reviewed; re-injection not needed",
                               now=NOW)
            self.assertTrue(g["proceed"])
            self.assertEqual((g["ack"]["by"], g["ack"]["at"], g["ack"]["verdict"]),
                             ("staffexample", "2026-10-17T00:00:00Z", verdict))
            self.assertEqual(os.stat(qb.record_paths(s)["ack"]).st_mode & 0o777, 0o640)

    def test_good_and_no_qc_on_record_proceed(self):
        s = self.session_with(baseline() + [run(FIRST - timedelta(hours=2))],
                              self.project() + [entry("A.raw", "Orbitrap Astral", FIRST)])
        g = qb.delivery_gate(s)
        self.assertEqual(g["verdict"], "no_qc_on_record")
        self.assertTrue(g["proceed"])
        self.assertFalse(g["needs_ack"])
        self.assertEqual(g["no_qc_on_record"], ["Orbitrap Astral"])

    def test_a_new_record_needs_a_new_acknowledgement(self):
        s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=0)])
        qb.acknowledge(s, "staffexample", "looked", now=NOW)
        self.assertTrue(qb.delivery_gate(s)["proceed"])
        with open(qb.record_paths(s)["json"]) as fh:
            rec = json.load(fh)
        rec["checked_at"] = "2026-10-18T00:00:00Z"               # re-run: another record
        with open(qb.record_paths(s)["json"], "w") as fh:
            json.dump(rec, fh)
        self.assertFalse(qb.delivery_gate(s)["proceed"])

    def test_no_record_and_unchecked_need_an_acknowledgement(self):
        s = os.path.join(self.tmp, "bare")
        os.makedirs(os.path.join(s, "logs"))
        os.makedirs(os.path.join(s, "input"))
        g = qb.delivery_gate(s)
        self.assertFalse(g["proceed"])
        self.assertIn("were not checked", g["reason"])
        s = self.session_with(fake=FakeStan([], fail=urllib.error.URLError("down")))
        g = qb.delivery_gate(s)
        self.assertFalse(g["proceed"])
        self.assertIn("stan_unreachable", g["reason"])

    def test_a_stale_record_cannot_be_acknowledged_only_re_run(self):
        s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=0)])
        with open(os.path.join(s, "input", "raw_files.txt"), "a") as fh:
            fh.write("/data/PROT_0000/S99.d\n")                    # a file added since
        g = qb.delivery_gate(s)
        self.assertTrue(g["needs_rerun"])
        self.assertIn("re-run step 8e", g["reason"])
        with self.assertRaises(ValueError):
            qb.acknowledge(s, "staffexample", "looked")
        self.assertFalse(os.path.exists(qb.record_paths(s)["ack"]))

    def test_every_staff_only_output_carries_the_marker(self):
        s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=0)])
        qb.acknowledge(s, "staffexample", "looked", now=NOW)
        qb.write_gate_snapshot(s, qb.delivery_gate(s))
        p = qb.record_paths(s)
        for f in (p["json"], p["md"], p["ack"], p["gate"]):
            self.assertTrue(qb.is_staff_only_file(f), f)
            self.assertEqual(os.stat(f).st_mode & 0o777, 0o640, f)
        with open(p["json"]) as fh:
            self.assertEqual(next(iter(json.load(fh))), "staff_only_marker")
        bare = os.path.join(self.tmp, "logs_copy", "qc_bracket.json")   # 2.10.0: no marker
        os.makedirs(os.path.dirname(bare))
        with open(bare, "w") as fh:
            json.dump({"verdict": "concern"}, fh)
        self.assertTrue(qb.is_staff_only_file(bare))
        self.assertFalse(qb.is_staff_only_file(os.path.join(SCRIPTS, "qc_bracket.py")))
        other = os.path.join(self.tmp, "plain.json")
        with open(other, "w") as fh:
            json.dump({"verdict": "good"}, fh)
        self.assertFalse(qb.is_staff_only_file(other))
        self.assertEqual(qb.public_gate(qb.delivery_gate(s)),
                         {"proceed": True, "record_sha256": qb._sha256(p["json"]),
                          "ack_at": "2026-10-17T00:00:00Z"})

    def test_an_agent_or_a_full_name_cannot_acknowledge(self):
        """2.10 safety review: `ack --by "Claude"` passed."""
        s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=0)])
        for by in ("Claude", "claude", "claude-code", "assistant", "qc-agent", "my-bot", "ai",
                   "Staff Example", "Brett Phinney"):
            with self.assertRaises(ValueError, msg=by):
                qb.acknowledge(s, by, "looked")
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            g = qb.acknowledge(s, "talbot", "looked")         # a login, not an agent's name
        self.assertTrue(g["proceed"])
        self.assertTrue(g["ack"]["staff_list"].startswith("NOT VERIFIED"))
        self.assertIn("staff list is not set up", err.getvalue())

    def test_a_trusted_staff_list_decides_who_may_acknowledge(self):
        """staff.py is the one check (staff.require): its list, trusted under the test switch
        with this account as the Core admin, as test_staff_decisions does."""
        import getpass
        s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=0)])
        d = os.path.join(self.tmp, "config")
        os.makedirs(d, mode=0o755)
        os.chmod(d, 0o755)
        path = os.path.join(d, "core_staff.txt")
        with open(path, "w") as fh:
            fh.write("# Core staff\nstaffexample\n")
        os.chmod(path, 0o644)
        me = getpass.getuser()
        knobs = {staff.STAFF_FILE_ENV: path, "SKILL_SLACK_TEST_LOOPBACK": "1",
                 "SKILL_CORE_ADMINS": me}
        with mock.patch.dict(os.environ, knobs):
            with self.assertRaises(ValueError) as cm:
                qb.acknowledge(s, "talbot", "looked")
            self.assertIn("not on the Core staff list", str(cm.exception))
            g = qb.acknowledge(s, "staffexample", "looked")
            self.assertTrue(g["ack"]["staff_list"].startswith("on the Core staff list"))
            os.chmod(path, 0o664)                       # group-writable: not trusted
            self.assertIsNone(staff.load_staff()[0])
        with mock.patch.dict(os.environ, dict(knobs, SKILL_CORE_ADMINS="someoneelse")):
            os.chmod(path, 0o644)
            self.assertIn("not to a Core admin", staff.load_staff()[1])

    def test_no_record_needs_a_reason_and_its_ack_never_releases_a_later_record(self):
        s = os.path.join(self.tmp, "norec")
        for d in ("logs", "input"):
            os.makedirs(os.path.join(s, d))
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("".join(e["file"] + "\n" for e in self.project()))
        with self.assertRaises(ValueError) as cm:
            qb.acknowledge(s, "staffexample", "looked")
        self.assertIn("--no-record-reason", str(cm.exception))
        with contextlib.redirect_stderr(io.StringIO()):
            g = qb.acknowledge(s, "staffexample", "looked",
                               no_record_reason="STAN down all week; Brett approved")
        self.assertTrue(g["proceed"])
        self.assertEqual(g["ack"]["no_record_reason"], "STAN down all week; Brett approved")
        _rc, rec, _ = self.check(baseline() + [run(FIRST + timedelta(hours=1), ids=0)],
                                 self.project())
        os.replace(os.path.join(self.tmp, "qc_bracket.json"), qb.record_paths(s)["json"])
        self.assertFalse(qb.delivery_gate(s)["proceed"])      # the record written later is held
        with self.assertRaises(ValueError):                    # and has a record: no reason now
            qb.acknowledge(s, "staffexample", "looked", no_record_reason="x")

    def test_an_acknowledgement_needs_a_name_and_one_line(self):
        s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=0)])
        for by, note in (("", "ok"), ("staffexample", ""), ("staffexample", "two\nlines"),
                         ("staffexample", "x" * (qb.NOTE_MAX + 1))):
            with self.assertRaises(ValueError):
                qb.acknowledge(s, by, note)
        self.assertFalse(os.path.exists(qb.record_paths(s)["ack"]))
        with open(qb.record_paths(s)["ack"], "w") as fh:          # a damaged file never counts
            fh.write("{")
        self.assertFalse(qb.delivery_gate(s)["proceed"])

    def test_gate_and_ack_cli_exit_codes(self):
        s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=0)])
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(qb.main(["gate", "--session", s]), qb.EXIT_USAGE)
            self.assertEqual(qb.main(["ack", "--session", s, "--by", "staffexample",
                                      "--note", "reviewed"]), qb.EXIT_OK)
            self.assertEqual(qb.main(["gate", "--session", s]), qb.EXIT_OK)

    def test_the_run_registry_carries_the_staff_section_and_the_gate(self):
        import record_run
        s = self.session_with(baseline() + [run(FIRST + timedelta(hours=1), ids=0)])
        q = record_run.qc_facts(s)
        self.assertEqual((q["verdict"], q["schema_version"]), ("concern", 1))
        self.assertFalse(q["gate"]["proceed"])
        text = "\n".join(record_run.render_qc(q))
        self.assertIn("## Instrument QC around this project (staff only -- never delivered)", text)
        self.assertIn("**Verdict:** CONCERN", text)
        self.assertIn("HELD until a staff member acknowledges", text)
        qb.acknowledge(s, "staffexample", "reviewed the affected files", now=NOW)
        text = "\n".join(record_run.render_qc(record_run.qc_facts(s)))
        self.assertIn('acknowledged by staffexample at 2026-10-17T00:00:00Z: "reviewed the '
                      'affected files"', text)
        bare = os.path.join(self.tmp, "bare_reg")
        os.makedirs(bare)
        self.assertIn("was not run", "\n".join(record_run.render_qc(record_run.qc_facts(bare))))

    def test_the_session_zip_leaves_every_qc_record_and_delivery_json_out(self):
        """2.10 safety review: delivery.json carried the gate (verdict, who, note) into the zip,
        and a record written with --out under another name escaped the exact-name rule."""
        import zipfile
        env = dict(os.environ, RECORD_RUN="off", SKILL_CONFIG_DIR=os.path.join(self.tmp, "cfg"),
                   CLAUDE_CONFIG_DIR=os.path.join(self.tmp, "cc"))
        env.pop("CLAUDE_CODE_SESSION_ID", None)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init", "--name",
                            "QC zip", "--base", os.path.join(self.tmp, "base")],
                           capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        sdir = json.loads(r.stdout)["paths"]["session_dir"]
        _rc, rec, _ = self.check(baseline() + [run(FIRST + timedelta(hours=1), ids=0)],
                                 self.project(), extra=())
        p = qb.record_paths(sdir)
        with open(p["json"], "w") as fh:                       # a concern record ...
            json.dump(rec, fh)
        with open(p["md"], "w") as fh:
            fh.write(qb.staff_markdown(rec) + "\n")
        os.makedirs(os.path.join(sdir, "output", "tables"), exist_ok=True)
        with open(os.path.join(sdir, "output", "tables", "qc.json"), "w") as fh:
            json.dump(rec, fh)                                 # ... and one under another name
        with open(os.path.join(sdir, "logs", "qc_bracket_ack.json"), "w") as fh:
            json.dump({"acks": [{"by": "staffexample", "verdict": "concern"}]}, fh)   # 2.10.0: no marker
        gate = {"verdict": "concern", "ack": {"by": "staffexample", "note": "CANARYNOTE"}}
        qb.write_gate_snapshot(sdir, gate)
        with open(os.path.join(sdir, "delivery.json"), "w") as fh:     # an old full delivery.json
            json.dump({"applied": True, "qc_gate": gate}, fh)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "finalize",
                            "--dir", sdir, "--zip", "--no-deposit", "--no-notify"],
                           capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        z = zipfile.ZipFile(sdir + ".zip")
        leaks = [(n, c) for n in z.namelist() if "/scripts/" not in n
                 for c in ('"verdict": "concern"', "CANARYNOTE", "staffexample", "CONCERN")
                 if c in z.read(n).decode("utf-8", "replace")]
        self.assertEqual(leaks, [])
        excluded = json.loads(r.stdout)["zip_excluded"]
        label = [k for k in excluded if k.startswith("staff-only QC records")][0]
        self.assertEqual(excluded[label], 6)    # json, md, ack, qc_gate.json, tables/qc.json, delivery.json
        self.assertTrue([n for n in z.namelist() if n.endswith("/scripts/qc_bracket.py")])

if __name__ == "__main__":
    unittest.main()
