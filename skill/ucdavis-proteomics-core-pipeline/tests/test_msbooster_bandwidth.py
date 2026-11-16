#!/usr/bin/env python3
"""MSBooster's RT LOESS bandwidth, read from FragPipe's log and recorded as MSBooster's choice
(2.15, Ruling 106).

MSBooster (FragPipe 24.0 runs v1.4.14) fits each run's retention times with a LOESS whose
bandwidth it CHOOSES per run by a grid search. The same input searched twice by the same skill,
workflow and install gave 0.10 in one search and 0.11 in another for one run of the two-run
timsTOF DDA reference (T12), and that alone moved the result: 6,797 PSMs / 187 protein groups
against 6,762 / 186 (everything before MSBooster byte-identical). So the value that ran is
recorded, from FragPipe's own log, as chosen by MSBooster -- never as a setting of the skill or
the workflow -- and "not stated in the log" where the log does not give it.

What these tests pin:
  * the reader (fragpipe_params.msbooster_rt_bandwidths) takes each run's RT bandwidth -- the
    line after "RT regression", never an IM regression's -- for 0.10 and 0.11, a run whose RT
    regression did not run (not stated), a log with no MSBooster, no log; the newest log wins;
    only file names, never a path;
  * FragPipe's logs are listed once (fragpipe_params.logs), newest first, for every reader;
  * the record reaches search_provenance.json, fragpipe_quantity.json and
    fragpipe_identifications.json.
The fixtures are invented runs in the line shapes of a real FragPipe 24.0 log. Stdlib only, plus
pyarrow for the quantified adapter.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
FIX = os.path.join(HERE, "fixtures", "fragpipe24_msbooster_log")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import fragpipe_params as fpp  # noqa: E402

try:
    import pyarrow  # noqa: F401
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

KEY = "msbooster_rt_bandwidth"


def folder_with(d, *logs):
    """`d` holding FragPipe logs: (fixture name, the name FragPipe gives it), in that order."""
    os.makedirs(d, exist_ok=True)
    for fixture, name in logs:
        shutil.copy(os.path.join(FIX, fixture), os.path.join(d, name))
    return d


def by_run(rec):
    return {r["run"]: r for r in rec["runs"]}


class TheReader(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def read(self, fixture):
        folder_with(self.d, (fixture, "log_2026-11-05_00-27-23.txt"))
        return fpp.msbooster_rt_bandwidths(self.d)

    def test_each_runs_rt_bandwidth_as_msbooster_chose_it(self):
        for fixture, want in (("log_bw010.txt", 0.1), ("log_bw011.txt", 0.11)):
            shutil.rmtree(self.d, True)
            rec = self.read(fixture)
            runs = by_run(rec)
            self.assertEqual(sorted(runs), ["S1-a", "S2-b"], fixture)
            self.assertEqual(runs["S1-a"]["bandwidth"], want, fixture)
            # the RT line, never the IM regression's (1.0 / 0.14 in the fixture)
            self.assertEqual(runs["S2-b"]["bandwidth"], 0.17, fixture)
            self.assertEqual(runs["S1-a"]["rt_regression_psms"], 126)
            self.assertEqual(runs["S1-a"]["experiment"], "S1_a")
            self.assertEqual(runs["S1-a"]["grid"], "0.01,0.05,0.1,0.2")
            for r in rec["runs"]:
                self.assertEqual(r["source"], fpp.MSBOOSTER_CHOSEN)
            self.assertEqual(rec["status"], "stated")
            self.assertEqual(rec["log"], "log_2026-11-05_00-27-23.txt")

    def test_it_is_said_to_be_msboosters_choice_never_a_setting(self):
        rec = self.read("log_bw011.txt")
        self.assertEqual(fpp.MSBOOSTER_CHOSEN, "chosen by MSBooster")
        self.assertTrue(rec["source"].startswith("chosen by MSBooster"), rec["source"])
        self.assertIn("not a setting", rec["source"])

    def test_a_run_whose_rt_regression_did_not_run_is_not_stated(self):
        runs = by_run(self.read("log_one_unstated.txt"))
        self.assertEqual(runs["S1-a"]["bandwidth"], fpp.NOT_STATED)
        self.assertEqual(runs["S1-a"]["source"], fpp.NOT_STATED)
        self.assertEqual(fpp.NOT_STATED, "not stated in the log")
        self.assertEqual(runs["S2-b"]["bandwidth"], 0.17)      # the IM line was not taken

    def test_a_log_without_msbooster_states_none(self):
        rec = self.read("log_no_msbooster.txt")
        self.assertEqual(rec["runs"], [])
        self.assertEqual(rec["status"], fpp.NOT_STATED)
        self.assertEqual(rec["log"], "log_2026-11-05_00-27-23.txt")

    def test_no_log_at_all(self):
        rec = fpp.msbooster_rt_bandwidths(self.d)
        self.assertEqual((rec["runs"], rec["status"], rec["log"]), ([], fpp.NOT_STATED, None))
        self.assertEqual(fpp.msbooster_rt_bandwidths(None)["status"], fpp.NOT_STATED)

    def test_the_newest_log_is_the_search_that_left_these_files(self):
        folder_with(self.d, ("log_bw011.txt", "log_2026-11-04_21-33-20.txt"),
                    ("log_bw010.txt", "log_2026-11-05_00-27-23.txt"))
        rec = fpp.msbooster_rt_bandwidths(self.d)
        self.assertEqual(by_run(rec)["S1-a"]["bandwidth"], 0.1)
        self.assertEqual(rec["log"], "log_2026-11-05_00-27-23.txt")

    def test_no_path_reaches_the_record(self):
        rec = self.read("log_bw010.txt")
        self.assertNotIn("/data/", json.dumps(rec))
        self.assertNotIn(self.d, json.dumps(rec))


class OneListOfLogs(unittest.TestCase):
    def test_every_reader_walks_the_same_list_newest_first(self):
        with tempfile.TemporaryDirectory() as d:
            folder_with(d, ("log_no_msbooster.txt", "log_2026-11-04_21-33-20.txt"),
                        ("log_no_msbooster.txt", "log_2026-11-05_00-27-23.txt"))
            self.assertEqual([os.path.basename(p) for p in fpp.logs(d)],
                             ["log_2026-11-05_00-27-23.txt", "log_2026-11-04_21-33-20.txt"])
            self.assertEqual(fpp.logs(os.path.join(d, "none")), [])
            import fragpipe_quant
            seen = []
            real = fpp.logs
            with mock.patch.object(fpp, "logs", side_effect=lambda f: seen.append(f) or real(f)):
                fpp.find_filter_runs(d)
                fragpipe_quant.ionquant_maxlfq_line(d)
                fpp.msbooster_rt_bandwidths(d)
            self.assertEqual(seen, [d, d, d])


class TheRecords(unittest.TestCase):
    """search_provenance.json, fragpipe_identifications.json and fragpipe_quantity.json carry it"""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_identifications_only_adapt_only(self):
        from test_fragpipe_ids_only import SHIPPED, IdentificationsOnly
        h = IdentificationsOnly("setUp")
        h.setUp()
        self.addCleanup(shutil.rmtree, h.d, True)
        self.assertEqual(h.generate().returncode, 0)
        self.assertEqual(h.run_job().returncode, 0)
        with open(os.path.join(h.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(SHIPPED)
        folder_with(h.out, ("log_bw011.txt", "log_2026-11-05_00-27-23.txt"))
        p = h.generate("--adapt-only")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        with open(os.path.join(h.out, "search_provenance.json")) as fh:
            prov = json.load(fh)[KEY]
        with open(os.path.join(h.out, "fragpipe_identifications.json")) as fh:
            ids = json.load(fh)[KEY]
        for rec in (prov, ids):
            self.assertEqual(by_run(rec)["S1-a"]["bandwidth"], 0.11)
            self.assertEqual(by_run(rec)["S1-a"]["source"], "chosen by MSBooster")

    def test_no_log_is_recorded_as_not_stated(self):
        from test_fragpipe_ids_only import SHIPPED, IdentificationsOnly
        h = IdentificationsOnly("setUp")
        h.setUp()
        self.addCleanup(shutil.rmtree, h.d, True)
        self.assertEqual(h.generate().returncode, 0)
        self.assertEqual(h.run_job().returncode, 0)
        with open(os.path.join(h.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(SHIPPED)
        self.assertEqual(h.generate("--adapt-only").returncode, 0)
        with open(os.path.join(h.out, "search_provenance.json")) as fh:
            self.assertEqual(json.load(fh)[KEY]["status"], "not stated in the log")

    @unittest.skipUnless(HAVE_ARROW, "pyarrow")
    def test_the_quantity_record(self):
        import fragpipe_quant as fq
        import run_search
        from test_fragpipe_dda_quantity import write_output
        write_output(self.d, ["S-1", "S-2", "S-3"])                 # MaxLFQ: no venv needed
        folder_with(self.d, ("log_bw010.txt", "log_2026-11-05_00-27-23.txt"))
        with mock.patch.dict(os.environ, {fq.VENV_ENV: os.path.join(self.d, "no_venv")}):
            run_search.adapt_fragpipe_dda(self.d)
        with open(os.path.join(self.d, fq.RECORD)) as fh:
            rec = json.load(fh)[KEY]
        with open(os.path.join(self.d, "search_provenance.json")) as fh:
            prov = json.load(fh)[fq.PROVENANCE_KEY][KEY]
        for r in (rec, prov):
            self.assertEqual(by_run(r)["S1-a"]["bandwidth"], 0.1)
            self.assertEqual(r["status"], "stated")


if __name__ == "__main__":
    unittest.main()
