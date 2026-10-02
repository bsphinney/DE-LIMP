#!/usr/bin/env python3
"""Repeats in a search's file list are okay, but they are flagged (Brett, 2026-10-01) -- by every
route, with one rule and one wording (check_report_runs).

A staff HT plate (2026-10-01): STAN's manifest listed one run 4 times and another twice -- 100
lines for 96 files -- with every gate PASS. Fed to run_search.py, DIA-NN refused them as "inputs
share a run name ... Rename them", which does not apply to one file listed twice; Sage refused
them as sharing an mzML name; the other engines took them, so those runs would have counted 2-4
times in the experiment and the empirical library.

  * the same path listed more than once: searched ONCE, flagged on stderr and recorded with its
    count in search_provenance.json `repeated_paths`, which the analysis report shows;
  * different files sharing a run name (one name in two folders): flagged, and the search stops
    before anything is written, naming them -- an engine would merge them into one run;
  * a mix of both: both flags.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import check_report_runs as crr  # noqa: E402
import make_analysis_html  # noqa: E402
from test_single_shot_sbatch import LIBFREE_CFG, _Workspace  # noqa: E402
from test_sage_sbatch_conversion import _Harness  # noqa: E402


def twin_in(d, path, folder="rerun"):
    """The same file name in another folder: a rerun, or one name on two plates."""
    os.makedirs(os.path.join(d, folder), exist_ok=True)
    twin = os.path.join(d, folder, os.path.basename(path))
    open(twin, "w").close()
    return twin


class TheRule(unittest.TestCase):
    def test_one_file_once_in_first_seen_order_judged_on_the_resolved_path(self):
        with tempfile.TemporaryDirectory() as d:
            a, b = os.path.join(d, "a.d"), os.path.join(d, "b.d")
            os.makedirs(a)
            os.makedirs(b)
            link = os.path.join(d, "link_to_b.d")
            os.symlink(b, link)
            kept, repeats = crr.distinct_inputs([b, a, b + "/", link, os.path.join(d, ".", "a.d"), b])
            self.assertEqual(kept, [b, a])
            self.assertEqual(repeats, [
                {"path": b, "times": 4, "also_listed_as": sorted([b + "/", link])},
                {"path": a, "times": 2, "also_listed_as": [os.path.join(d, ".", "a.d")]}])

    def test_one_name_in_two_folders_both_kept(self):
        files = ["/plate1/s1.d", "/plate2/s1.d", "/plate1/s2.d"]
        self.assertEqual(crr.distinct_inputs(files)[0], files)
        self.assertEqual(crr.repeated_names(files),
                         [{"run_name": "s1", "paths": ["/plate1/s1.d", "/plate2/s1.d"]}])

    def test_a_mix_is_worded_both_ways(self):
        files = ["/plate1/s1.d", "/plate2/s1.d", "/plate1/s1.d", "/plate1/s2.d"]
        kept, paths = crr.distinct_inputs(files)
        names = crr.repeated_names(kept)
        self.assertEqual([(r["path"], r["times"]) for r in paths], [("/plate1/s1.d", 2)])
        self.assertEqual([r["run_name"] for r in names], ["s1"])
        text = crr.repeats_text(paths, names)
        self.assertEqual(text, [
            "1 input file(s) were listed more than once; each is searched once: /plate1/s1.d "
            "(listed 2 times)",
            "1 run name(s) are shared by different files, all kept: s1: /plate1/s1.d, /plate2/s1.d"])
        note = crr.repeats_note(paths, "x", names)
        self.assertEqual(note.count("[x] FLAG: "), 2)

    def test_the_report_check_counts_a_repeated_path_once(self):
        files = ["/x/s1.d", "/x/s2.d", "/x/s1.d"]
        ok, msg = crr.verify("report.parquet", files, parquet_runs=lambda _r: {"s1", "s2"})
        self.assertTrue(ok, msg)
        self.assertIn("all 2 runs", msg)


class RunSearchDiann(unittest.TestCase):
    def prov(self, w):
        with open(os.path.join(w.out, "search_provenance.json")) as fh:
            return json.load(fh)

    def test_one_path_twice_is_searched_once_and_flagged(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, n_runs=2)
            link = os.path.join(d, "alias_1.mzML")
            os.symlink(w.runs[1], link)
            w.runs = [w.runs[0], w.runs[1], w.runs[0], link]
            p = w.generate()
            self.assertNotIn("share a run name", p.stderr)
            self.assertIn(f"[run_search] FLAG: 2 input file(s) were listed more than once; each is "
                          f"searched once: {w.runs[0]} (listed 2 times)", p.stderr)
            self.assertIn(f"(listed 2 times, also as {link})", p.stderr)
            prov = self.prov(w)
            self.assertEqual(prov["n_files"], 2)
            self.assertEqual(prov["files"], [w.runs[0], w.runs[1]])
            self.assertEqual([(r["path"], r["times"]) for r in prov["repeated_paths"]],
                             [(w.runs[0], 2), (w.runs[1], 2)])
            with open(w.search_job) as fh:
                self.assertEqual(fh.read().count(w.runs[0]), 1)

    def test_one_name_in_two_folders_is_flagged_and_stops_naming_both(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, n_runs=2)
            twin = twin_in(d, w.runs[0])
            w.runs = [w.runs[0], w.runs[1], twin]
            p = w.generate(check=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn(f"[run_search] FLAG: 1 run name(s) are shared by different files, all "
                          f"kept: sample_0: {w.runs[0]}, {twin}", p.stderr)
            self.assertIn(f"STOPPED before searching: different input files share a run name -- "
                          f"sample_0: {w.runs[0]}, {twin}. They are kept in the file list and "
                          f"flagged (repeated_names)", p.stderr)
            self.assertFalse(os.path.exists(w.out), "nothing may be written")
            self.assertFalse(os.path.exists(w.lib_job))

    def test_a_mix_flags_both_and_stops(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, n_runs=2)
            twin = twin_in(d, w.runs[0])
            w.runs = [w.runs[0], w.runs[1], w.runs[1], twin]
            p = w.generate(check=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertEqual(p.stderr.count("[run_search] FLAG:"), 2)
            self.assertIn(f"{w.runs[1]} (listed 2 times)", p.stderr)
            self.assertIn(f"sample_0: {w.runs[0]}, {twin}", p.stderr)

    def test_diann_parallel_run_directly_follows_the_same_rule(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
            out = os.path.join(d, "chain_out")
            argv = [sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"), "--diann", w.diann,
                    "--fasta", w.fasta, "--out", out, "--cfg", w.cfg, "--raw"]
            p = subprocess.run(argv + w.runs + [w.runs[2]], capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn(f"[diann_parallel] FLAG: 1 input file(s) were listed more than once", p.stderr)
            with open(os.path.join(out, "file_list.txt")) as fh:
                self.assertEqual([ln.strip() for ln in fh if ln.strip()], w.runs)
            twin = os.path.join(d, "plate2", "sample_0.d")
            os.makedirs(twin)
            p = subprocess.run(argv + w.runs + [twin], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn(f"[diann_parallel] STOPPED before searching", p.stderr)
            self.assertIn(twin, p.stderr)


class RunSearchFragPipe(unittest.TestCase):
    def test_fragpipe_holds_to_the_same_rule(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, n_runs=2)
            twin = twin_in(d, w.runs[0])
            tools = os.path.join(d, "tools_fp.json")
            with open(tools, "w") as fh:
                json.dump({"fragpipe": os.path.join(d, "fake-fragpipe"),
                           "versions": {"fragpipe": "24.0"}}, fh)
            p = subprocess.run(
                [sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--tools", tools,
                 "--bundle", w.bundle, "--engine", "fragpipe", "--params", w.cfg,
                 "--fasta", w.fasta, "--out", w.out, "--files", w.runs[0], twin,
                 "--threads", "4", "--sbatch", w.job], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("STOPPED before searching", p.stderr)
            self.assertFalse(os.path.exists(w.out))


class RunSearchSage(_Harness):
    def test_one_path_twice_is_converted_and_searched_once(self):
        once = os.path.join(self.d, "sage_once.sh")
        p = self.run_search("--sbatch", once, files=self.raws)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        twice = os.path.join(self.d, "sage_twice.sh")
        p = self.run_search("--sbatch", twice, files=[self.raws[0], self.raws[1], self.raws[0]])
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertNotIn("share an mzML name", p.stderr)
        self.assertIn(f"FLAG: 1 input file(s) were listed more than once; each is searched once: "
                      f"{self.raws[0]} (listed 2 times)", p.stderr)
        with open(once) as a, open(twice) as b:
            self.assertEqual(b.read().count(self.raws[0]), a.read().count(self.raws[0]))
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["n_files"], 2)
        self.assertEqual(prov["repeated_paths"][0], {"path": self.raws[0], "times": 2,
                                                      "also_listed_as": []})

    def test_one_name_in_two_folders_stops_before_converting(self):
        twin = twin_in(self.d, self.raws[0])
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job, files=[self.raws[0], twin])
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("STOPPED before searching", p.stderr)
        self.assertFalse(os.path.exists(job))
        self.assertEqual(self.calls(), [], "nothing may be converted")


class TheReportSaysSo(unittest.TestCase):
    def test_the_analysis_report_shows_the_flag(self):
        with tempfile.TemporaryDirectory() as session:
            search = os.path.join(session, "output", "search")
            os.makedirs(search)
            rec = {"engine": "diann", "repeated_paths": [{"path": "/x/s1.d", "times": 4,
                                                          "also_listed_as": []}]}
            with open(os.path.join(search, "search_provenance.json"), "w") as fh:
                json.dump(rec, fh)
            note = make_analysis_html.search_inputs_note({}, session)
            self.assertEqual(note["kind"], "info")
            self.assertEqual(note["text"], "1 input file(s) were listed more than once; each is "
                                           "searched once: /x/s1.d (listed 4 times).")
            with open(os.path.join(search, "search_provenance.json"), "w") as fh:
                json.dump({"engine": "diann", "repeated_paths": []}, fh)
            self.assertIsNone(make_analysis_html.search_inputs_note({}, session))


if __name__ == "__main__":
    unittest.main()
