#!/usr/bin/env python3
"""
run_search reuses step 2's analysis.tdf check instead of re-reading every .d on the login node
(I05).

detect_acquisition.py (step 2) already ran bruker_tdf.tdf_integrity() on every .d and the session
keeps its output as input/acquisition.json. run_search.py then read every analysis.tdf again --
on HIVE, from the login node, for a few hundred .d. Now each record carries a fingerprint
(bruker_tdf.fingerprint: the sizes and mtimes of analysis.tdf and analysis.tdf_bin, from os.stat),
and a .d whose fingerprint is unchanged is not read again; its recorded status still decides the
gate. A record from before 2.14.2 has no fingerprint and counts as not recorded (Review Focus 1).
More than TDF_LOGIN_MAX .d still to read on a login node: the whole command moves to a compute
node (on_compute.ensure_compute_node), the skill's established pattern.

Synthetic .d only (the one deliberate writer, synthetic_tdf_write_uri); on_compute is patched.
"""
import json
import os
import sqlite3
import struct
import sys
import tempfile
import time
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

from synthetic_tdf import synthetic_tdf_write_uri  # noqa: E402  (fixture writes only)

import bruker_tdf  # noqa: E402
import on_compute  # noqa: E402
import run_search  # noqa: E402

BLOCK, N = 256, 400


def make_d(tmp, name):
    """A finished, intact .d: a frame index that reaches the end of its tdf_bin."""
    d = os.path.join(tmp, name)
    os.makedirs(d)
    with open(os.path.join(d, "analysis.tdf_bin"), "wb") as fh:
        for _ in range(N):
            fh.write(struct.pack("<II", BLOCK, 5) + b"\0" * (BLOCK - 8))
    tdf = os.path.join(d, "analysis.tdf")
    con = sqlite3.connect(synthetic_tdf_write_uri(tdf), uri=True)
    con.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, TimsId INTEGER)")
    con.executemany("INSERT INTO Frames VALUES (?,?)", [(i + 1, i * BLOCK) for i in range(N)])
    con.commit()
    con.close()
    for side in ("-wal", "-shm"):
        if os.path.exists(tdf + side):
            os.remove(tdf + side)
    return d


def session(tmp, n, record=True, fingerprint=True):
    """<tmp>/input/acquisition.json for n .d (as step 2 saves it) -> (files, out dir)."""
    files = [make_d(tmp, f"run{i}.d") for i in range(n)]
    if record:
        recs = []
        for f in files:
            integ = bruker_tdf.tdf_integrity(f)
            if not fingerprint:
                integ.pop("fingerprint", None)
            recs.append({"file": f, "vendor": "Bruker", "tdf_integrity": integ})
        os.makedirs(os.path.join(tmp, "input"))
        with open(os.path.join(tmp, "input", "acquisition.json"), "w") as fh:
            json.dump({"files": recs}, fh)
    return files, os.path.join(tmp, "output", "search")


def refuse_reads(*a, **k):
    raise AssertionError("tdf_integrity() was called: a recorded .d was read again")


class Fingerprint(unittest.TestCase):
    def test_from_stat_and_in_every_integrity_record(self):
        with tempfile.TemporaryDirectory() as t:
            d = make_d(t, "a.d")
            fp = bruker_tdf.fingerprint(d)
            # sizes, mtimes, and what no copy tool can set: ctime and inode (L8R-F6)
            self.assertEqual(set(fp), {"tdf_bytes", "tdf_mtime_ns", "tdf_ctime_ns", "tdf_ino",
                                       "tdf_bin_bytes", "tdf_bin_mtime_ns", "tdf_bin_ctime_ns",
                                       "tdf_bin_ino"})
            self.assertEqual(bruker_tdf.tdf_integrity(d)["rules"], bruker_tdf.INTEGRITY_RULES)
            self.assertEqual(fp["tdf_bin_bytes"], N * BLOCK)
            self.assertEqual(bruker_tdf.tdf_integrity(d)["fingerprint"], fp)
            self.assertEqual(bruker_tdf.tdf_integrity(d)["status"], "ok")


class Reuse(unittest.TestCase):
    def test_find_acquisition_json_walks_up_from_out(self):
        with tempfile.TemporaryDirectory() as t:
            _, out = session(t, 1)
            want = os.path.join(t, "input", "acquisition.json")
            self.assertEqual(run_search.find_acquisition_json(out), want)
            self.assertEqual(run_search.find_acquisition_json(out, "/given.json"), "/given.json")
            self.assertIsNone(run_search.find_acquisition_json(os.path.join(t, "a", "b", "c",
                                                                            "d", "e")))

    def test_matching_fingerprint_is_not_read_again(self):
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 3)
            rec = run_search.recorded_integrity(files, run_search.find_acquisition_json(out))
            self.assertEqual(sorted(rec), sorted(files))
            with mock.patch.object(run_search, "tdf_integrity", refuse_reads):
                self.assertEqual(run_search.tdf_problems(files, rec), [])

    def test_a_recorded_problem_still_stops_the_search(self):
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 2)
            acq = run_search.find_acquisition_json(out)
            with open(acq) as fh:
                doc = json.load(fh)
            doc["files"][0]["tdf_integrity"].update(status="truncated", statuses=["truncated"],
                                                    problems=["analysis.tdf is truncated"])
            with open(acq, "w") as fh:
                json.dump(doc, fh)
            rec = run_search.recorded_integrity(files, acq)
            with mock.patch.object(run_search, "tdf_integrity", refuse_reads):
                bad = run_search.tdf_problems(files, rec)
            self.assertEqual([p for p, _ in bad], [files[0]])

    def test_changed_mtime_is_read_again(self):
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 2)
            later = time.time() + 120
            os.utime(os.path.join(files[1], "analysis.tdf"), (later, later))
            rec = run_search.recorded_integrity(files, run_search.find_acquisition_json(out))
            self.assertEqual(list(rec), [files[0]])
            calls = []
            real = run_search.tdf_integrity
            with mock.patch.object(run_search, "tdf_integrity",
                                   lambda p: calls.append(p) or real(p)):
                run_search.tdf_problems(files, rec)
            self.assertEqual(calls, [files[1]])

    def test_a_copy_with_preserved_times_is_read_again(self):
        """L8R-F6: a .d replaced at the same path by `rsync -a` / `cp -p` keeps sizes and mtimes;
        its inode and ctime change, and step 2's verdict on the old bytes is not reused."""
        import shutil
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 2)
            for name in ("analysis.tdf", "analysis.tdf_bin"):
                src = os.path.join(files[1], name)
                tmp = src + ".copy"
                shutil.copy2(src, tmp)                 # contents, mode and times preserved
                os.replace(tmp, src)
                self.assertEqual(os.stat(src).st_mtime_ns, os.stat(src).st_mtime_ns)
            rec = run_search.recorded_integrity(files, run_search.find_acquisition_json(out))
            self.assertEqual(list(rec), [files[0]])

    def test_a_record_made_under_other_rules_is_read_again(self):
        """L8R-F6: a later skill with stricter tdf_integrity() rules never trusts an older
        verdict: a record with no `rules`, or other rules, is re-read."""
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 2)
            acq = run_search.find_acquisition_json(out)
            with open(acq) as fh:
                doc = json.load(fh)
            doc["files"][0]["tdf_integrity"].pop("rules", None)
            doc["files"][1]["tdf_integrity"]["rules"] = bruker_tdf.INTEGRITY_RULES - 1
            with open(acq, "w") as fh:
                json.dump(doc, fh)
            self.assertEqual(run_search.recorded_integrity(files, acq), {})

    def test_a_side_file_that_appeared_since_is_read_again(self):
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 1)
            with open(os.path.join(files[0], "analysis.tdf-wal"), "wb") as fh:
                fh.write(b"\0" * 4096)
            self.assertEqual(run_search.recorded_integrity(
                files, run_search.find_acquisition_json(out)), {})

    def test_old_acquisition_json_without_fingerprint_is_read_again(self):
        """Review Focus 1: a record from before 2.14.2 is a re-read, not a crash."""
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 2, fingerprint=False)
            self.assertEqual(run_search.recorded_integrity(
                files, run_search.find_acquisition_json(out)), {})
            calls = []
            real = run_search.tdf_integrity
            with mock.patch.object(run_search, "tdf_integrity",
                                   lambda p: calls.append(p) or real(p)):
                self.assertEqual(run_search.tdf_problems(files, {}), [])
            self.assertEqual(calls, files)

    def test_an_unreadable_acquisition_json_is_a_re_read(self):
        with tempfile.TemporaryDirectory() as t:
            files, _ = session(t, 1, record=False)
            bad = os.path.join(t, "acq.json")
            with open(bad, "w") as fh:
                fh.write("{not json")
            self.assertEqual(run_search.recorded_integrity(files, bad), {})
            self.assertEqual(run_search.recorded_integrity(files, None), {})


class LoginNode(unittest.TestCase):
    def _gate(self, files, out, login=True, queue=None):
        moved = []

        def ensure(why, **kw):
            moved.append((why, kw))
            raise SystemExit(0)                # the real one never returns
        with mock.patch.object(on_compute, "on_login_node", lambda env=None: login), \
                mock.patch.object(on_compute, "ensure_compute_node", ensure), \
                mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(run_search.TDF_REEXEC_ENV, None)
            try:
                run_search.check_tdf(files, out, None, allow=False, queue=queue)
            except SystemExit:
                pass
            marked = os.environ.get(run_search.TDF_REEXEC_ENV)
        return moved, marked

    def test_six_unchecked_on_a_login_node_moves_to_a_compute_node(self):
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 6, record=False)
            with mock.patch.object(run_search, "tdf_integrity", refuse_reads):
                moved, marked = self._gate(files, out)
            self.assertEqual(len(moved), 1)
            why, kw = moved[0]
            self.assertEqual(why, "run_search.py reads the analysis.tdf of 6 .d")
            self.assertEqual(kw["minutes"], 10)
            self.assertEqual(kw["mem"], "2G")
            self.assertEqual(os.path.realpath(kw["script"]),
                             os.path.realpath(run_search.__file__))
            self.assertEqual(marked, "1", "the search under that srun must know why it is there")

    def test_the_move_takes_run_search_s_own_queue(self):
        """L8R-F7: --account / --partition / --qos given to run_search go to the srun; the check
        only reads, so a preemptible queue may take it."""
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 6, record=False)
            with mock.patch.object(run_search, "tdf_integrity", refuse_reads):
                moved, _ = self._gate(files, out, queue={"partition": "high",
                                                         "account": "labbgrp", "qos": None,
                                                         "exclude": None, "mem": "96G"})
            (why, kw), = moved
            self.assertEqual((kw["account"], kw["partition"], kw["qos"]),
                             ("labbgrp", "high", None))
            self.assertIs(kw["preemptible_ok"], True)
            self.assertEqual(kw["mem"], "2G", "the search's --mem is not the check's")

    def test_six_recorded_are_not_moved(self):
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 6)
            with mock.patch.object(run_search, "tdf_integrity", refuse_reads):
                moved, _ = self._gate(files, out)
            self.assertEqual(moved, [])

    def test_five_or_fewer_are_read_here(self):
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 5, record=False)
            calls = []
            real = run_search.tdf_integrity
            with mock.patch.object(run_search, "tdf_integrity",
                                   lambda p: calls.append(p) or real(p)):
                moved, _ = self._gate(files, out)
            self.assertEqual(moved, [])
            self.assertEqual(sorted(calls), sorted(files))

    def test_off_a_login_node_they_are_read_here(self):
        with tempfile.TemporaryDirectory() as t:
            files, out = session(t, 6, record=False)
            moved, _ = self._gate(files, out, login=False)
            self.assertEqual(moved, [])

    def test_minutes_grow_with_the_count(self):
        self.assertEqual(run_search.tdf_check_minutes(6), 10)
        self.assertEqual(run_search.tdf_check_minutes(300), 25)

    def test_no_search_runs_inline_in_the_srun_it_moved_to(self):
        """That srun is 2G and minutes long: generating jobs there is fine, searching is not."""
        import subprocess
        import test_single_shot_mass_accuracy as ss
        with tempfile.TemporaryDirectory() as d:
            case = ss.SingleShotMassAccTests()
            raws, fasta, cfg, tools, bundle = case._setup(d)
            env = dict(case._env(d), **{run_search.TDF_REEXEC_ENV: "1", "SLURM_JOB_ID": "77"})
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", tools, "--bundle", bundle, "--params", cfg,
                                "--fasta", fasta, "--out", os.path.join(d, "out"),
                                "--files", *raws, "--engine", "diann", "--threads", "8"],
                               capture_output=True, text=True, env=env, timeout=240)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("REFUSING to run the search inline: this command was moved", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(d, "argv.txt")))


if __name__ == "__main__":
    unittest.main()
