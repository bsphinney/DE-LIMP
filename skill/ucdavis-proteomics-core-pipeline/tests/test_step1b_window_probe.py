#!/usr/bin/env python3
"""
Step 1b -- the scan-window probe -- has to measure runs that are representative of the cohort,
and must never measure a run whose data DIA-NN cannot actually read.

Upstream (PR #70, c67e658/484cb9f) already gives step 1b the .NET 8 export, a per-probe
watchdog, a process-group stop, and a fall-through "try the first file, then the next, up to 3".
What was still missing, each measured on HIVE with DIA-NN 2.7.0 (FRAN re-search pilot,
2026-09-16):

  1. WHICH runs. DIA-NN's README says auto-optimised values "depend on which run is first in the
     list" and recommends "several representative runs". The first file of a listing is whatever
     sorts first. On a timsTOF cohort one-file probes pinned 10, 11 or 14 depending on the run;
     the median of three representative runs was stable. So the probe is handed the whole
     cohort, ranks it, keeps blanks/washes/failed acquisitions out, probes the median and
     quartile runs, and pins the MEDIAN radius. A run that logs no radius is replaced by the
     next representative run (the fall-through upstream introduced), not by the next file.

  2. A TRUNCATED Bruker index. 342 .d on HIVE have an analysis.tdf whose frame index stops early
     while analysis.tdf_bin is complete: a read-write sqlite open checkpointed a stale
     mid-acquisition -wal into the finished file. Step 1b's size rule picked one in the pilot
     (tdf_bin full size, index 0.7%): DIA-NN read 121 cycles, 80 precursors. So a .d is never
     probed when its tdf header is in WAL mode, a non-empty -wal/-journal sits beside it, or its
     last frame block ends short of the end of tdf_bin -- and the tdf is only ever opened with
     `?mode=ro&immutable=1` (mode=ro alone still reads the stale WAL and leaves -wal/-shm files).
     Bruker runs are ranked by acquisition time (max Frames.Time), not bytes.

  3. Evidence and time. A failure still writes window.json with every probe; the job log shows
     each probe live (a timsTOF step 1b took 1083 s, past watch_run.sh's 15-minute stall rule);
     one hung run is cut at --timeout so the next gets its turn, and --budget keeps all probes
     inside the job's wall clock so the evidence is written before SLURM kills the job.

The fake DIA-NN below mimics the real one where it matters: it refuses .raw without DOTNET_ROOT
(printing the real error and exiting 0), prints the real "Scan window radius set to N" line, then
keeps running until the probe stops it.
"""
import ast
import glob
import hashlib
import json
import os
import random
import re
import shutil
import signal
import sqlite3
import stat
import struct
import subprocess
import sys
import tempfile
import threading
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import probe_window  # noqa: E402
import diann_parallel as dp  # noqa: E402
# Building a synthetic .d is a deliberate WRITE to an analysis.tdf, and it is named so
# that tests/test_tdf_readonly_open.py can tell it from a reader that must never write.
from synthetic_tdf import synthetic_tdf_write_uri  # noqa: E402

PROBE = os.path.join(SCRIPTS, "probe_window.py")
GB = 1024 ** 3
MB = 1024 ** 2

# Radius is encoded in the run name (`..._w8.raw`) so a test controls what each run "measures".
# `nowin` never logs a radius (a run DIA-NN cannot calibrate); `hang` goes silent.
FAKE_DIANN = r"""#!/bin/bash
f=""; th=""
while [ $# -gt 0 ]; do
  case "$1" in --f) f="$2"; shift;; --threads) th="$2"; shift;; esac; shift
done
echo "fake DIA-NN: file=$f threads=$th"
case "$f" in
  *.raw|*.RAW)
    if [ -z "${DOTNET_ROOT:-}" ]; then
      echo "ERROR: cannot read .raw files, please download and install .NET Runtime 8: 8.0.17 or later https://dotnet.microsoft.com/en-us/download/dotnet/8.0 : 1"
      exit 0
    fi;;
esac
case "$f" in *nowin*) echo "Finished"; exit 0;; esac
case "$f" in *hang*) exec sleep 60;; esac
r=$(basename "$f" | sed -n 's/.*_w\([0-9][0-9]*\).*/\1/p')
echo "[0:04] Scan window radius set to ${r:-7}"
exec sleep 20
"""

FAKE_DOTNET_ENV = {"DOTNET_ROOT": "/opt/fake-dotnet"}


def _exe(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def _read(path):
    with open(path) as fh:
        return fh.read()


def _sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def _raw(d, name, size):
    """A sparse file of the given logical size: getsize() sees `size`, the disk sees ~0."""
    p = os.path.join(d, name)
    with open(p, "wb") as fh:
        fh.truncate(size)
    return p


def _bruker(d, name, minutes=60.0, gb=2.4, nframes=140, indexed_frames=None):
    """A Bruker timsTOF .d as far as a reader of the index cares.

    analysis.tdf_bin holds `nframes` frame blocks; each block starts with its uint32 size, as in
    a real file. analysis.tdf is an SQLite database whose Frames table indexes the first
    `indexed_frames` of them (all by default) with TimsId = the block's offset and Time in
    seconds. indexed_frames=1 of 140 is the pilot's truncated file: tdf_bin full, index 0.7%."""
    p = os.path.join(d, name)
    os.makedirs(p)
    total = int(gb * GB)
    block = total // nframes
    offsets = [i * block for i in range(nframes)]
    sizes = [block] * (nframes - 1) + [total - block * (nframes - 1)]
    with open(os.path.join(p, "analysis.tdf_bin"), "wb") as fh:
        for off, size in zip(offsets, sizes):
            fh.seek(off)
            fh.write(struct.pack("<I", size))
        fh.truncate(total)
    k = nframes if indexed_frames is None else indexed_frames
    con = sqlite3.connect(synthetic_tdf_write_uri(os.path.join(p, "analysis.tdf")),
                          uri=True)
    con.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL NOT NULL, "
                "TimsId INTEGER NOT NULL, NumScans INTEGER NOT NULL)")
    con.executemany("INSERT INTO Frames VALUES (?, ?, ?, ?)",
                    [(i + 1, (i + 1) * minutes * 60.0 / nframes, offsets[i], 918)
                     for i in range(k)])
    con.commit()
    con.close()
    return p


def _set_wal_header(dotd):
    """Bytes 18-19 of the SQLite header = 2,2: the file is in WAL mode, as every one of the 342
    truncated tdfs on HIVE is (intact ones are 1,1)."""
    with open(os.path.join(dotd, "analysis.tdf"), "r+b") as fh:
        fh.seek(18)
        fh.write(bytes([2, 2]))


def _stale_wal_copy(d, name):
    """Reproduce the damage mechanism's precondition for real: a writer in WAL mode with frames
    still in its -wal, copied mid-acquisition (as the instrument->HIVE copier did). A READ-WRITE
    open of the copy checkpoints those pages into it; a correct reader must leave it alone."""
    src = os.path.join(d, "_writer_" + name)
    os.makedirs(src)
    tdf = os.path.join(src, "analysis.tdf")
    con = sqlite3.connect(synthetic_tdf_write_uri(tdf), uri=True, isolation_level=None)
    con.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL NOT NULL, "
                "TimsId INTEGER NOT NULL, NumScans INTEGER NOT NULL)")
    con.execute("INSERT INTO Frames VALUES (1, 0.1, 0, 918)")
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA wal_autocheckpoint=0")
    con.executemany("INSERT INTO Frames VALUES (?, ?, ?, 918)",
                    [(i, i * 0.1, (i - 1) * 1000) for i in range(2, 200)])
    dest = os.path.join(d, name)
    os.makedirs(dest)
    shutil.copyfile(tdf, os.path.join(dest, "analysis.tdf"))
    shutil.copyfile(tdf + "-wal", os.path.join(dest, "analysis.tdf-wal"))
    con.close()
    with open(os.path.join(dest, "analysis.tdf_bin"), "wb") as fh:
        fh.truncate(199 * 1000)
    return dest


def _frames_without_usable_time(dotd, keep_column=False):
    """Rebuild Frames so no acquisition time can be read out of it: the Time column is gone
    (another vendor's schema, or an older one), or it is NULL for every frame. The frame INDEX --
    Id, TimsId, NumScans -- is untouched, so every spectrum is exactly where it was. PR #73's
    integrity checker calls such a run `ok`; ranking by time is what it loses, not its data."""
    cols = ("Id INTEGER PRIMARY KEY, Time REAL, TimsId INTEGER NOT NULL, NumScans INTEGER NOT NULL"
            if keep_column else
            "Id INTEGER PRIMARY KEY, TimsId INTEGER NOT NULL, NumScans INTEGER NOT NULL")
    pick = "Id, NULL, TimsId, NumScans" if keep_column else "Id, TimsId, NumScans"
    con = sqlite3.connect(synthetic_tdf_write_uri(os.path.join(dotd, "analysis.tdf")),
                          uri=True)
    con.execute("CREATE TABLE F2 (%s)" % cols)
    con.execute("INSERT INTO F2 SELECT %s FROM Frames" % pick)
    con.execute("DROP TABLE Frames")
    con.execute("ALTER TABLE F2 RENAME TO Frames")
    con.commit()
    con.close()
    return dotd


def _fasta_lib(d):
    fasta = os.path.join(d, "db.fasta")
    with open(fasta, "w") as fh:
        fh.write(">sp|P1|X\nPEPTIDER\n")
    lib = os.path.join(d, "lib.speclib")
    with open(lib, "w") as fh:
        fh.write("lib")
    return fasta, lib


def _fake_dotnet_root(d):
    """A directory ensure_dotnet8.sh accepts as a .NET 8 >= 8.0.17 install (it also needs
    AspNetCore 8), so the REAL dotnet_prefix() code path runs in a test without downloading
    anything."""
    root = os.path.join(d, "dotnet8")
    os.makedirs(root)
    _exe(os.path.join(root, "dotnet"),
         '#!/bin/bash\necho "Microsoft.AspNetCore.App 8.0.28 [%s/shared/Microsoft.AspNetCore.App]"\n'
         'echo "Microsoft.NETCore.App 8.0.28 [%s/shared/Microsoft.NETCore.App]"\n' % (root, root))
    return root


def _names(entries):
    return [os.path.basename(x["file"].rstrip("/")) for x in entries]


# --------------------------------------------------------------------------------------------
class BrukerIndexIntegrityTests(unittest.TestCase):
    """A .d is only as long as its analysis.tdf frame index says. A truncated index reads as a
    perfectly valid, short run -- DIA-NN reports no error, it just searches 0.7% of the data."""

    def test_an_intact_run_reports_its_acquisition_time_and_a_complete_index(self):
        with tempfile.TemporaryDirectory() as d:
            p = _bruker(d, "run.d", minutes=60, gb=2.4)
            st = probe_window.bruker_tdf_status(p)
            self.assertEqual(st["problems"], [])
            self.assertAlmostEqual(st["time_s"], 3600.0, places=3)
            self.assertEqual(st["indexed_bytes"], st["tdf_bin_bytes"])
            self.assertAlmostEqual(st["index_coverage"], 1.0)

    def test_the_pilots_truncated_index_is_a_problem_even_with_a_full_tdf_bin(self):
        with tempfile.TemporaryDirectory() as d:
            p = _bruker(d, "trunc.d", minutes=21, gb=2.4, nframes=140, indexed_frames=1)
            st = probe_window.bruker_tdf_status(p)
            self.assertTrue(st["problems"], st)
            self.assertLess(st["index_coverage"], 0.01)
            self.assertRegex(" ".join(st["problems"]), r"0\.7%")

    def test_an_index_that_stops_short_of_the_end_by_more_than_a_tenth_of_a_percent(self):
        with tempfile.TemporaryDirectory() as d:
            p = _bruker(d, "almost.d", nframes=2000, indexed_frames=1996)   # 99.8%
            self.assertTrue(probe_window.bruker_tdf_status(p)["problems"])
            q = _bruker(d, "whole.d", nframes=2000)
            self.assertEqual(probe_window.bruker_tdf_status(q)["problems"], [])

    def test_a_wal_mode_header_is_a_problem_and_the_tdf_is_not_opened(self):
        with tempfile.TemporaryDirectory() as d:
            p = _bruker(d, "wal.d")
            _set_wal_header(p)
            st = probe_window.bruker_tdf_status(p)
            self.assertRegex(" ".join(st["problems"]), r"WAL")

    def test_a_stale_wal_copy_is_refused_and_left_exactly_as_it_was(self):
        with tempfile.TemporaryDirectory() as d:
            p = _stale_wal_copy(d, "stale.d")
            tdf = os.path.join(p, "analysis.tdf")
            with open(tdf, "rb") as fh:
                self.assertEqual(fh.read(20)[18:20], bytes([2, 2]), "fixture: header not WAL")
            before = (_sha(tdf), os.path.getsize(tdf + "-wal"), sorted(os.listdir(p)))
            st = probe_window.bruker_tdf_status(p)
            self.assertTrue(st["problems"])
            probe_window.measure_run(p)
            self.assertEqual((_sha(tdf), os.path.getsize(tdf + "-wal"), sorted(os.listdir(p))),
                             before, "reading the .d changed it (a checkpoint, or new -shm)")

            # The fixture reproduces the hazard: a plain READ-WRITE open rewrites the copy.
            victim = os.path.join(d, "victim.d")
            shutil.copytree(p, victim)
            victim_tdf = os.path.join(victim, "analysis.tdf")
            con = sqlite3.connect(synthetic_tdf_write_uri(victim_tdf), uri=True)
            con.execute("SELECT COUNT(*) FROM Frames").fetchone()
            con.close()
            self.assertNotEqual(_sha(os.path.join(victim, "analysis.tdf")), before[0])

    def test_a_non_empty_wal_or_journal_beside_an_intact_tdf_is_a_problem(self):
        """Intact tdfs with a stale -wal still exist on HIVE (e.g. apr26/blankDia_S1-H11): the
        next read-write open by anything truncates them."""
        for side in ("analysis.tdf-wal", "analysis.tdf-journal"):
            with self.subTest(side), tempfile.TemporaryDirectory() as d:
                p = _bruker(d, "run.d")
                with open(os.path.join(p, side), "wb") as fh:
                    fh.write(b"\x37\x7f\x06\x82" + b"\0" * 28)
                self.assertRegex(" ".join(probe_window.bruker_tdf_status(p)["problems"]),
                                 re.escape(side))
                open(os.path.join(p, side), "wb").close()          # empty is harmless
                self.assertEqual(probe_window.bruker_tdf_status(p)["problems"], [])

    def test_the_tdf_is_opened_immutable_and_nothing_appears_beside_it(self):
        with tempfile.TemporaryDirectory() as d:
            p = _bruker(d, "has space #1 ?x.d")
            tdf = os.path.join(p, "analysis.tdf")
            before = (_sha(tdf), sorted(os.listdir(p)))
            st = probe_window.bruker_tdf_status(p)
            self.assertEqual(st["problems"], [])
            self.assertEqual((_sha(tdf), sorted(os.listdir(p))), before)
            self.assertIn("mode=ro&immutable=1", probe_window.tdf_uri(tdf))

    def test_every_sqlite_open_in_the_probe_goes_through_the_immutable_uri(self):
        src = _read(PROBE)
        calls = re.findall(r"sqlite3\.connect\(([^\n]*)", src)
        self.assertTrue(calls, "no sqlite3.connect found -- the tdf reader moved?")
        for args in calls:
            self.assertIn("tdf_uri(", args)
            self.assertIn("uri=True", args)

    def test_a_tdf_bin_without_its_tdf_and_the_reverse_are_problems(self):
        with tempfile.TemporaryDirectory() as d:
            p = _bruker(d, "no_tdf.d")
            os.remove(os.path.join(p, "analysis.tdf"))
            self.assertTrue(probe_window.bruker_tdf_status(p)["problems"])
            q = _bruker(d, "no_bin.d")
            os.remove(os.path.join(q, "analysis.tdf_bin"))
            self.assertTrue(probe_window.bruker_tdf_status(q)["problems"])

    def test_a_frames_table_without_a_time_column_is_a_schema_not_a_damaged_run(self):
        """Asking for MAX(Time), MAX(TimsId) and COUNT(*) in one statement reported a Frames
        table with no Time column as "frame index cannot be read" -- damage -- and the run was
        never probed. The index is intact; only the RANKING key is missing, so the cohort falls
        back to size. (PR #73's checker calls this run ok; this is the disagreement.)"""
        with tempfile.TemporaryDirectory() as d:
            p = _frames_without_usable_time(_bruker(d, "notime.d", minutes=60, gb=2.4))
            st = probe_window.bruker_tdf_status(p)
            self.assertEqual(st["problems"], [], st)
            self.assertAlmostEqual(st["index_coverage"], 1.0)
            self.assertEqual(st["n_frames"], 140)
            self.assertIsNone(st["time_s"])
            self.assertRegex(st["time_problem"], r"Time")
            self.assertRegex(st["time_problem"], r"ranked by acquisition time")

    def test_a_frames_time_that_is_null_for_every_frame_is_not_damage_either(self):
        with tempfile.TemporaryDirectory() as d:
            p = _frames_without_usable_time(_bruker(d, "nulltime.d", gb=0.5), keep_column=True)
            st = probe_window.bruker_tdf_status(p)
            self.assertEqual(st["problems"], [], st)
            self.assertIsNone(st["time_s"])
            self.assertRegex(st["time_problem"], r"NULL")

    def test_a_frames_table_that_is_really_unreadable_is_still_damage(self):
        """The split must not turn a missing index into a ranking footnote."""
        with tempfile.TemporaryDirectory() as d:
            p = _bruker(d, "noframes.d", gb=0.5)
            noframes_tdf = os.path.join(p, "analysis.tdf")
            con = sqlite3.connect(synthetic_tdf_write_uri(noframes_tdf), uri=True)
            con.execute("DROP TABLE Frames")
            con.commit()
            con.close()
            self.assertRegex(" ".join(probe_window.bruker_tdf_status(p)["problems"]),
                             r"frame index cannot be read")

    def test_a_directory_that_is_not_a_bruker_run_has_no_tdf_status(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "other.d"))
            self.assertIsNone(probe_window.bruker_tdf_status(os.path.join(d, "other.d")))


# --------------------------------------------------------------------------------------------
class RepresentativeSelectionTests(unittest.TestCase):
    """Which runs step 1b measures."""

    def test_raw_size_is_the_file_size(self):
        with tempfile.TemporaryDirectory() as d:
            m = probe_window.measure_run(_raw(d, "a.raw", 3 * GB))
            self.assertEqual(m["size_bytes"], 3 * GB)
            self.assertIsNone(m["time_s"])

    def test_a_missing_input_is_unreadable(self):
        self.assertFalse(probe_window.measure_run("/no/such/run.raw")["readable"])

    def test_three_or_more_runs_probe_the_median_first_then_the_quartiles(self):
        with tempfile.TemporaryDirectory() as d:
            sizes = {"r1": 1.00, "r2": 1.10, "r3": 1.20, "r4": 1.30, "r5": 1.40}
            paths = [_raw(d, f"{k}.raw", int(v * GB)) for k, v in sizes.items()]
            sel = probe_window.select_representative(paths)
            self.assertEqual([(os.path.basename(c["file"]), c["role"]) for c in sel["chosen"]],
                             [("r3.raw", "median"), ("r2.raw", "lower_quartile"),
                              ("r4.raw", "upper_quartile")])
            self.assertEqual(sel["rank_by"], "size")

    def test_reserves_are_the_remaining_runs_nearest_the_median_first(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"r{i}.raw", (10 + i) * GB // 10) for i in range(9)]
            sel = probe_window.select_representative(paths)
            chosen = set(_names(sel["chosen"]))
            reserves = _names(sel["reserves"])
            self.assertEqual(set(reserves) | chosen, {f"r{i}.raw" for i in range(9)})
            self.assertFalse(set(reserves) & chosen)
            # median r4; quartiles r2, r6; nearest remaining: r3/r5 (larger first), then r1/r7
            self.assertEqual(reserves[:4], ["r5.raw", "r3.raw", "r7.raw", "r1.raw"])

    def test_selection_does_not_depend_on_input_order(self):
        """raws[0] was the whole bug: the answer must be a property of the cohort."""
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"r{i}.raw", (10 + i) * GB // 10) for i in range(8)]
            want = probe_window.select_representative(paths)
            rng = random.Random(7)
            for _ in range(5):
                shuffled = paths[:]
                rng.shuffle(shuffled)
                got = probe_window.select_representative(shuffled)
                self.assertEqual(_names(got["chosen"]), _names(want["chosen"]))
                self.assertEqual(_names(got["reserves"]), _names(want["reserves"]))

    def test_a_blank_first_file_is_never_probed_or_held_in_reserve(self):
        with tempfile.TemporaryDirectory() as d:
            blank = _raw(d, "00_blank.raw", 40 * MB)
            runs = [_raw(d, f"s{i}.raw", (12 + i) * GB // 10) for i in range(5)]
            sel = probe_window.select_representative([blank] + runs)
            self.assertNotIn("00_blank.raw", _names(sel["chosen"]) + _names(sel["reserves"]))
            self.assertIn("00_blank.raw", _names(sel["excluded_small"]))

    def test_blanks_injected_between_every_sample_are_never_probed(self):
        """With a blank after every sample, blanks are the MAJORITY: the median of all runs is a
        blank's size and a floor measured from it excludes nothing. The reference is the median
        of the LARGER half."""
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for i in range(6):
                paths.append(_raw(d, f"blank{i}.raw", 50 * MB))
                paths.append(_raw(d, f"s{i}.raw", 3 * GB + i))
            paths.append(_raw(d, "blank6.raw", 50 * MB))
            sel = probe_window.select_representative(paths)
            picked = _names(sel["chosen"]) + _names(sel["reserves"])
            self.assertEqual(len(sel["chosen"]), 3)
            self.assertTrue(all(c.startswith("s") for c in picked), picked)
            self.assertEqual(len(sel["excluded_small"]), 7)

    def test_a_wash_beside_a_few_dia_runs_is_not_probed(self):
        """Real sizes, HIVE raw_data/Lumos1/noi25: three 90-min DIA runs and a 60-min wash."""
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, "FL011125_DWang-Dia_90m_440.raw", 1753040624),
                     _raw(d, "FL011125_DWang-Dia_90m_441.raw", 1954554235),
                     _raw(d, "FL011125_DWang-Dia_90m_442.raw", 1852751309),
                     _raw(d, "FL24noi25_wa60m-cleanAPInewAnlypreC.raw", 550651027)]
            sel = probe_window.select_representative(paths)
            self.assertNotIn("FL24noi25_wa60m-cleanAPInewAnlypreC.raw",
                             _names(sel["chosen"]) + _names(sel["reserves"]))
            self.assertEqual(len(sel["chosen"]), 3)

    def test_fewer_than_three_runs_probe_the_median_and_keep_the_other_in_reserve(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, "a.raw", GB), _raw(d, "b.raw", GB + GB // 10)]
            sel = probe_window.select_representative(paths)
            self.assertEqual([(os.path.basename(c["file"]), c["role"]) for c in sel["chosen"]],
                             [("b.raw", "median")])
            self.assertEqual(_names(sel["reserves"]), ["a.raw"])

    def test_max_probes_one_keeps_the_single_median_run(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"r{i}.raw", (10 + i) * GB // 10) for i in range(5)]
            sel = probe_window.select_representative(paths, max_probes=1)
            self.assertEqual(_names(sel["chosen"]), ["r2.raw"])

    def test_unreadable_inputs_are_reported_and_never_chosen(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"r{i}.raw", GB) for i in range(3)] + ["/no/such/run.raw"]
            sel = probe_window.select_representative(paths)
            self.assertEqual(sel["unreadable"], ["/no/such/run.raw"])
            self.assertNotIn("/no/such/run.raw", [c["file"] for c in sel["chosen"]])

    def test_no_probeable_input_is_an_error_not_a_guess(self):
        with self.assertRaises(ValueError):
            probe_window.select_representative(["/no/such/a.raw", "/no/such/b.raw"])

    # ---- Bruker ----------------------------------------------------------------------------
    def test_timstof_runs_are_ranked_by_acquisition_time_not_bytes(self):
        """Bytes are what misled the pilot. By bytes the median here is the 80-min run; by the
        index's max(Time) it is the 60-min run."""
        with tempfile.TemporaryDirectory() as d:
            spec = [(40, 3.0), (50, 3.4), (60, 3.1), (70, 3.3), (80, 3.2)]
            paths = [_bruker(d, f"t{mins}.d", minutes=mins, gb=gb) for mins, gb in spec]
            sel = probe_window.select_representative(paths)
            self.assertEqual(sel["rank_by"], "time")
            self.assertEqual(sel["chosen"][0]["role"], "median")
            self.assertEqual(os.path.basename(sel["chosen"][0]["file"]), "t60.d")
            self.assertAlmostEqual(sel["chosen"][0]["time_s"], 3600.0, places=3)

    def test_the_pilots_truncated_run_is_never_probed_though_its_tdf_bin_is_largest(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_bruker(d, f"ok{i}.d", minutes=21, gb=2.3 + i / 100) for i in range(4)]
            trunc = _bruker(d, "8aug25_Koganti_trunc.d", minutes=21, gb=2.5, nframes=140,
                            indexed_frames=1)
            sel = probe_window.select_representative([trunc] + paths)
            self.assertNotIn("8aug25_Koganti_trunc.d",
                             _names(sel["chosen"]) + _names(sel["reserves"]))
            damaged = {os.path.basename(x["file"]): x for x in sel["excluded_damaged"]}
            self.assertIn("8aug25_Koganti_trunc.d", damaged)
            self.assertTrue(damaged["8aug25_Koganti_trunc.d"]["problems"])
            self.assertEqual(len(sel["chosen"]), 3)

    def test_wal_mode_and_stale_wal_runs_are_never_probed(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_bruker(d, f"ok{i}.d", gb=2.0 + i / 10) for i in range(3)]
            wal = _bruker(d, "walmode.d", gb=2.2)
            _set_wal_header(wal)
            stale = _bruker(d, "stalewal.d", gb=2.25)
            with open(os.path.join(stale, "analysis.tdf-wal"), "wb") as fh:
                fh.write(b"x" * 4096)
            sel = probe_window.select_representative(paths + [wal, stale])
            self.assertEqual(sorted(_names(sel["excluded_damaged"])), ["stalewal.d", "walmode.d"])
            self.assertEqual(sorted(_names(sel["chosen"])), ["ok0.d", "ok1.d", "ok2.d"])

    def test_a_short_bruker_acquisition_is_not_probed(self):
        """A 10-min wash on the same bytes as 60-min runs: excluded by acquisition time."""
        with tempfile.TemporaryDirectory() as d:
            paths = [_bruker(d, f"s{i}.d", minutes=60 + i, gb=2.0) for i in range(4)]
            wash = _bruker(d, "wash.d", minutes=10, gb=2.0)
            sel = probe_window.select_representative(paths + [wash])
            small = {os.path.basename(x["file"]): x for x in sel["excluded_small"]}
            self.assertIn("wash.d", small)
            self.assertRegex(small["wash.d"]["why"], r"time")

    def test_a_full_length_bruker_blank_is_not_probed(self):
        """Same gradient, little signal: excluded by the indexed bytes."""
        with tempfile.TemporaryDirectory() as d:
            paths = [_bruker(d, f"s{i}.d", minutes=60, gb=2.0 + i / 10) for i in range(4)]
            blank = _bruker(d, "blank.d", minutes=60, gb=0.2)
            sel = probe_window.select_representative(paths + [blank])
            small = {os.path.basename(x["file"]): x for x in sel["excluded_small"]}
            self.assertIn("blank.d", small)
            self.assertRegex(small["blank.d"]["why"], r"size")

    def test_an_empty_failed_acquisition_does_not_switch_a_timstof_cohort_to_size(self):
        """HIVE, garg cohort (srun job 23537237): the folder's failed acquisition is a .d with
        neither analysis.tdf nor tdf_bin. It has no acquisition time, and deciding the ranking
        before the size floor dropped it made all eight real runs rank by bytes."""
        with tempfile.TemporaryDirectory() as d:
            failed = os.path.join(d, "260822_He200ng_failed.d")
            os.makedirs(failed)
            spec = [(40, 3.0), (50, 3.4), (60, 3.1), (70, 3.3), (80, 3.2)]
            paths = [_bruker(d, f"t{mins}.d", minutes=mins, gb=gb) for mins, gb in spec]
            sel = probe_window.select_representative([failed] + paths)
            self.assertEqual(sel["rank_by"], "time")
            self.assertIn("260822_He200ng_failed.d", _names(sel["excluded_small"]))
            self.assertEqual(os.path.basename(sel["chosen"][0]["file"]), "t60.d")

    def test_a_bruker_run_without_a_usable_time_falls_back_to_size_not_to_damaged(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_bruker(d, f"t{i}.d", minutes=60, gb=2.0 + i / 10) for i in range(4)]
            _frames_without_usable_time(paths[1])
            sel = probe_window.select_representative(paths)
            self.assertEqual(sel["excluded_damaged"], [])
            self.assertEqual(sel["rank_by"], "size")
            self.assertEqual(len(sel["chosen"]), 3)
            self.assertEqual(_names(sel["no_acquisition_time"]), ["t1.d"])

    # ---- one cohort, one kind of run -------------------------------------------------------
    def test_a_thermo_and_bruker_cohort_is_warned_about_not_pruned_to_one_instrument(self):
        """size_bytes is not one quantity: a .raw's is the file's bytes, a .d's the INDEXED bytes
        of analysis.tdf_bin. One floor over both excluded all three 0.9 GB Orbitrap .raw here as
        "likely a blank, wash or failed injection", flipped the ranking to time because only .d
        were left, and pinned a two-instrument cohort from one instrument."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"Ex_{i}.raw", 9 * GB // 10) for i in range(3)]
            dotds = [_bruker(d, f"tt{i}.d", minutes=60, gb=8.0 + i / 10) for i in range(3)]
            sel = probe_window.select_representative(raws + dotds)
            self.assertEqual(sel["excluded_small"], [])
            self.assertEqual(sel["n_eligible"], 6)
            self.assertEqual(sel["rank_by"], "size")
            self.assertRegex(sel["mixed_cohort"], r"mixes Bruker \.d and Thermo \.raw")
            self.assertRegex(sel["mixed_cohort"], r"not valid for both")
            self.assertEqual(sorted(sel["size_reference_by_kind"]), ["Bruker .d", "Thermo .raw"])
            self.assertIsNone(sel["reference_size_bytes"], "one floor for two kinds of run")

    def test_one_kind_of_run_still_reports_the_single_floor_it_measured(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"r{i}.raw", (10 + i) * GB // 10) for i in range(5)]
            sel = probe_window.select_representative(paths)
            self.assertIsNone(sel["mixed_cohort"])
            self.assertEqual(list(sel["size_reference_by_kind"]), ["Thermo .raw"])
            self.assertEqual(sel["min_size_bytes"], 0.5 * sel["reference_size_bytes"])

    def test_a_failed_acquisition_is_measured_against_the_real_d_beside_it(self):
        """A failed Bruker acquisition is a .d with neither analysis.tdf nor tdf_bin in it. It
        has to share the size floor with the real .d -- give it a partition of its own and it
        survives, has no acquisition time, and drags the cohort back to size ranking (the HIVE
        garg case, in the other direction)."""
        with tempfile.TemporaryDirectory() as d:
            failed = os.path.join(d, "260822_He200ng_failed.d")
            os.makedirs(failed)
            paths = [_bruker(d, f"t{m}.d", minutes=m, gb=3.0) for m in (40, 50, 60, 70, 80)]
            sel = probe_window.select_representative([failed] + paths)
            self.assertIsNone(sel["mixed_cohort"])
            self.assertEqual(sel["rank_by"], "time")
            self.assertIn("260822_He200ng_failed.d", _names(sel["excluded_small"]))

    # ---- a cohort that is mostly blanks ----------------------------------------------------
    def test_a_cohort_that_is_mostly_washes_says_the_median_run_looks_like_one(self):
        """The floor is measured from the cohort, and here the washes ARE the cohort: nothing is
        excluded and the median run picked is wash07 (verified). Nothing in this script can tell
        a wash from a sample, so it says what the shape of the cohort is instead of pretending."""
        with tempfile.TemporaryDirectory() as d:
            washes = [_raw(d, f"wash{i:02d}.raw", 150 * MB + i) for i in range(12)]
            samples = [_raw(d, f"sample{i}.raw", 2 * GB + i) for i in range(2)]
            sel = probe_window.select_representative(washes + samples)
            self.assertEqual(sel["excluded_small"], [])
            self.assertEqual(os.path.basename(sel["chosen"][0]["file"]), "wash07.raw")
            self.assertRegex(sel["median_much_smaller"], r"median run wash07\.raw")
            self.assertRegex(sel["median_much_smaller"], r"under half the largest run left")

    def test_nine_tiny_files_beside_one_real_run_say_so_too(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"t{i}.raw", 1) for i in range(9)] + [_raw(d, "big.raw", 10 * GB)]
            sel = probe_window.select_representative(paths)
            self.assertTrue(all(n.startswith("t") for n in _names(sel["chosen"])))
            self.assertRegex(sel["median_much_smaller"], r"big\.raw")

    def test_a_cohort_of_real_runs_is_not_warned_about_its_median(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"s{i}.raw", (20 + i) * GB // 10) for i in range(6)]
            self.assertIsNone(probe_window.select_representative(paths)["median_much_smaller"])

    # ---- the last resort: a cohort that is entirely in WAL mode -----------------------------
    def test_a_cohort_whose_tdfs_are_all_in_wal_mode_is_probed_with_a_warning(self):
        """Reproduced 6 of 6 in review: every .d in WAL mode failed step 1b outright, and steps
        2-5 -- which SEARCH those very runs -- sat in DependencyNeverSatisfied. The header is a
        sign of the damage, not the damage: these indexes are complete and nothing stale sits
        beside them."""
        with tempfile.TemporaryDirectory() as d:
            paths = [_bruker(d, f"w{i}.d", minutes=60, gb=2.0 + i / 10) for i in range(6)]
            for p in paths:
                _set_wal_header(p)
            sel = probe_window.select_representative(paths)
            self.assertEqual(len(sel["chosen"]), 3)
            self.assertEqual(sel["excluded_damaged"], [])
            self.assertEqual(_names(sel["probed_despite_wal_header"]),
                             [f"w{i}.d" for i in range(6)])
            self.assertRegex(sel["probed_despite_wal_header"][0]["problems"][0], r"WAL")
            self.assertEqual(sel["rank_by"], "time")

    def test_a_wal_mode_run_is_still_refused_while_any_intact_run_can_be_probed(self):
        with tempfile.TemporaryDirectory() as d:
            ok = [_bruker(d, f"ok{i}.d", minutes=60, gb=2.0 + i / 10) for i in range(3)]
            wal = _bruker(d, "wal.d", minutes=60, gb=2.05)
            _set_wal_header(wal)
            sel = probe_window.select_representative(ok + [wal])
            self.assertEqual(_names(sel["excluded_damaged"]), ["wal.d"])
            self.assertEqual(sel["probed_despite_wal_header"], [])
            self.assertNotIn("wal.d", _names(sel["chosen"]) + _names(sel["reserves"]))

    def test_a_wal_mode_run_with_a_truncated_index_is_never_the_last_resort(self):
        """The last resort is the HEADER alone. A short index is the damage the header warns
        about, and it stays refused even when it is all there is."""
        with tempfile.TemporaryDirectory() as d:
            paths = [_bruker(d, f"t{i}.d", minutes=21, gb=2.4, nframes=140, indexed_frames=1)
                     for i in range(3)]
            for p in paths:
                _set_wal_header(p)
            with self.assertRaises(probe_window.NoProbeableRun):
                probe_window.select_representative(paths)

    def test_a_stale_wal_beside_the_tdf_is_never_the_last_resort_either(self):
        with tempfile.TemporaryDirectory() as d:
            paths = [_bruker(d, f"s{i}.d", minutes=60, gb=2.0) for i in range(3)]
            for p in paths:
                _set_wal_header(p)
                with open(os.path.join(p, "analysis.tdf-wal"), "wb") as fh:
                    fh.write(b"x" * 4096)
            with self.assertRaises(probe_window.NoProbeableRun) as caught:
                probe_window.select_representative(paths)
            self.assertRegex(str(caught.exception), r"--raw does not override these checks")

    # ---- the floors, and what "nearest the median" means ------------------------------------
    def test_the_floors_can_never_empty_the_cohort(self):
        """The largest run of every kind clears both floors on its own, so the "nothing left"
        branch was dead code. It is an assert now; only an impossible floor reaches it."""
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"r{i}.raw", (10 + i) * GB // 10) for i in range(5)]
            sel = probe_window.select_representative(paths, min_fraction=1.0)
            self.assertTrue(sel["chosen"])
            with self.assertRaises(AssertionError):
                probe_window.select_representative(paths, min_fraction=2.0)

    def test_reserves_are_nearest_the_median_by_rank_position_and_it_says_so(self):
        """|i - mid| in the ranking, not |value - the median's value|: the two differ wherever
        the ranked sizes are spaced unevenly, and the docs have to say which one it is."""
        with tempfile.TemporaryDirectory() as d:
            paths = [_raw(d, f"r{i}.raw", v * GB // 10)
                     for i, v in enumerate((10, 11, 12, 13, 19, 20, 21))]
            sel = probe_window.select_representative(paths)
            self.assertEqual(_names(sel["chosen"]), ["r3.raw", "r1.raw", "r5.raw"])
            # r4 (1.9 GB) is six times further from r3's 1.3 GB than r2 (1.2 GB) is, and is
            # still the first reserve: one rank away, not one value away.
            self.assertEqual(_names(sel["reserves"]), ["r4.raw", "r2.raw", "r6.raw", "r0.raw"])
            self.assertRegex(probe_window.select_representative.__doc__, r"RANK POSITION")
            self.assertRegex(probe_window.SELECTION_RULE,
                             r"nearest the median in rank position")

    def test_a_directory_run_is_sized_by_every_file_under_it_not_its_top_level(self):
        """An Agilent .d keeps its data in AcqData/, a Waters .raw directory in _FUNC*.DAT.
        Counted at the top level only they all tied at ~0 bytes and the ranking fell back to the
        path tie-break -- alphabetical, which is the first-file order this probe exists to
        remove."""
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for name, mb in (("c_small.d", 20), ("a_big.d", 24), ("b_mid.d", 22)):
                acq = os.path.join(d, name, "AcqData")
                os.makedirs(acq)
                _raw(acq, "MSScan.bin", mb * MB)
                paths.append(os.path.join(d, name))
            self.assertEqual(probe_window.measure_run(paths[1])["size_bytes"], 24 * MB)
            sel = probe_window.select_representative(paths)
            self.assertEqual(_names(sel["chosen"]), ["b_mid.d", "c_small.d", "a_big.d"])

    def test_directories_without_a_tdf_are_ranked_by_size(self):
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for i in range(3):
                p = os.path.join(d, f"x{i}.d")
                os.makedirs(p)
                _raw(p, "data.bin", (10 + i) * MB)
                paths.append(p)
            sel = probe_window.select_representative(paths)
            self.assertEqual(sel["rank_by"], "size")
            self.assertEqual(len(sel["chosen"]), 3)


# --------------------------------------------------------------------------------------------
def _run_probe(d, raws, *extra, env_extra=None, fake=FAKE_DIANN, timeout=240):
    fasta, lib = _fasta_lib(d)
    diann = _exe(os.path.join(d, "diann"), fake)
    cfg = os.path.join(d, "resolved.cfg")
    if not os.path.exists(cfg):
        with open(cfg, "w") as fh:
            fh.write("--mass-acc 20\n--mass-acc-ms1 7\n")
    env = {k: v for k, v in os.environ.items() if k != "DOTNET_ROOT"}
    env.update(env_extra or {})
    t0 = time.time()
    p = subprocess.run([sys.executable, PROBE, "--diann", diann, "--raw", *raws,
                        "--fasta", fasta, "--lib", lib, "--threads", "12",
                        "--write-cfg", cfg, "--workdir", os.path.join(d, "probe_work"), *extra],
                       capture_output=True, text=True, env=env, timeout=timeout)
    return p, cfg, time.time() - t0


class ProbeManyRunsTests(unittest.TestCase):
    """probe_window.py end to end, against the fake DIA-NN."""

    def test_the_median_radius_is_pinned_and_every_run_is_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, "a_w9.raw", 11 * GB // 10), _raw(d, "b_w7.raw", 12 * GB // 10),
                    _raw(d, "c_w8.raw", 13 * GB // 10)]
            p, cfg, _ = _run_probe(d, raws, env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["window_radius"], 8)
            self.assertEqual(out["pin_as"], "--window 8")
            self.assertFalse(out["incomplete"])
            per_run = {os.path.basename(x["file"]): x["radius"] for x in out["probes"]}
            self.assertEqual(per_run, {"a_w9.raw": 9, "b_w7.raw": 7, "c_w8.raw": 8})
            for x in out["probes"]:
                for key in ("role", "size_bytes", "seconds", "timed_out", "log"):
                    self.assertIn(key, x)
            self.assertIn("rule", out["selection"])
            txt = _read(cfg)
            self.assertEqual(txt.count("--window"), 1, txt)
            self.assertIn("--window 8", txt)

    def test_a_run_that_logs_no_radius_is_replaced_by_the_next_representative_run(self):
        """Upstream's fall-through, kept: one bad run must not fail the cohort -- but its
        replacement is the next run nearest the median, not the next file in the listing, and
        the pin is still the median of three measured runs."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, "r1_w7.raw", 10 * GB // 10), _raw(d, "r2_w8.raw", 11 * GB // 10),
                    _raw(d, "r3_nowin.raw", 12 * GB // 10), _raw(d, "r4_w9.raw", 13 * GB // 10),
                    _raw(d, "r5_w8.raw", 14 * GB // 10)]
            p, cfg, _ = _run_probe(d, raws, env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(_names(out["probes"]),
                             ["r3_nowin.raw", "r2_w8.raw", "r4_w9.raw", "r5_w8.raw"])
            self.assertEqual(out["failed"], ["r3_nowin.raw"])
            self.assertEqual(out["probes"][3]["role"], "reserve")
            self.assertEqual(os.path.basename(out["probes"][3]["replaces"]), "r3_nowin.raw")
            self.assertEqual(out["radii"], [8, 9, 8])
            self.assertEqual(out["window_radius"], 8)
            self.assertFalse(out["incomplete"])
            self.assertIn("--window 8", _read(cfg))

    def test_fewer_radii_than_planned_are_pinned_with_a_warning_not_hidden(self):
        """Never worse than upstream, which pinned the first file that answered: one measured
        radius is still pinned -- but the JSON and the job log say it is incomplete."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, "a_w7.raw", 11 * GB // 10), _raw(d, "b_nowin.raw", 12 * GB // 10),
                    _raw(d, "c_nowin.raw", 13 * GB // 10)]
            p, cfg, _ = _run_probe(d, raws, env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["window_radius"], 7)
            self.assertTrue(out["incomplete"])
            self.assertEqual(sorted(out["failed"]), ["b_nowin.raw", "c_nowin.raw"])
            self.assertIn("WARNING", p.stderr)
            self.assertIn("--window 7", _read(cfg))

    def test_after_max_failures_with_no_radius_it_fails_and_still_writes_its_evidence(self):
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"r{i}_nowin.raw", (10 + i) * GB // 10) for i in range(6)]
            p, cfg, _ = _run_probe(d, raws, "--max-failures", "3", env_extra=FAKE_DOTNET_ENV)
            self.assertNotEqual(p.returncode, 0)
            out = json.loads(p.stdout)
            self.assertIsNone(out["window_radius"])
            self.assertEqual(len(out["probes"]), 3)
            self.assertEqual(out["stopped_because"], "max_failures")
            self.assertEqual(len(out["failed"]), 3)
            self.assertIn("excluded_small", out["selection"])
            for x in out["probes"]:
                self.assertTrue(os.path.isfile(x["log"]), x["log"])
            self.assertIn("Could not read the scan-window radius", p.stderr)
            self.assertNotIn("--window", _read(cfg))

    def test_a_missing_dotnet_stops_at_once_and_says_how_to_fix_it(self):
        """No other run can succeed without .NET, so falling through would only spend time."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"r{i}_w6.raw", (10 + i) * GB // 10) for i in range(5)]
            p, _, _ = _run_probe(d, raws)
            self.assertNotEqual(p.returncode, 0)
            out = json.loads(p.stdout)
            self.assertEqual(len(out["probes"]), 1)
            self.assertEqual(out["stopped_because"], "environment")
            self.assertIn("DOTNET_ROOT", p.stderr)
            m = re.search(r"bash (\S*ensure_dotnet8\.sh)", p.stderr)
            self.assertIsNotNone(m, p.stderr)
            helper = m.group(1).strip("'\"")
            self.assertTrue(os.path.isabs(helper), helper)
            self.assertTrue(os.path.isfile(helper), helper)

    def test_a_hung_run_is_cut_at_its_timeout_and_replaced(self):
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, "r1_w7.raw", 10 * GB // 10), _raw(d, "r2_w7.raw", 11 * GB // 10),
                    _raw(d, "r3_hang.raw", 12 * GB // 10), _raw(d, "r4_w7.raw", 13 * GB // 10)]
            p, _, elapsed = _run_probe(d, raws, "--timeout", "2", env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["probes"][0]["timed_out"], True)
            self.assertEqual(os.path.basename(out["probes"][0]["file"]), "r3_hang.raw")
            self.assertEqual(out["window_radius"], 7)
            self.assertEqual(len(out["radii"]), 3)
            self.assertLess(elapsed, 45, "the hung run was not cut at --timeout")

    def test_the_budget_bounds_all_probes_together(self):
        """--timeout is per probe; --budget is for everything, so the probes cannot outlive the
        step's own wall clock (SLURM would kill the job before window.json is written). Each
        probe here takes ~1.5 s against a 2 s budget: per-probe, all three would finish."""
        slow = FAKE_DIANN.replace('echo "[0:04] Scan window radius',
                                  'sleep 1.5; echo "[0:04] Scan window radius')
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"r{i}_w7.raw", (10 + i) * GB // 10) for i in range(3)]
            p, _, elapsed = _run_probe(d, raws, "--timeout", "60", "--budget", "2",
                                       env_extra=FAKE_DOTNET_ENV, fake=slow)
            out = json.loads(p.stdout)
            self.assertEqual(out["stopped_because"], "budget")
            self.assertLess(len(out["radii"]), 3)
            self.assertLess(elapsed, 30)
            if out["radii"]:
                self.assertEqual(p.returncode, 0, p.stderr)
                self.assertTrue(out["incomplete"])
            self.assertIn("[probe_window] budget", p.stderr)

    def test_a_budget_spent_before_any_radius_fails(self):
        slow = FAKE_DIANN.replace('echo "[0:04] Scan window radius',
                                  'sleep 5; echo "[0:04] Scan window radius')
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"r{i}_w7.raw", (10 + i) * GB // 10) for i in range(3)]
            p, _, _ = _run_probe(d, raws, "--budget", "1", env_extra=FAKE_DOTNET_ENV, fake=slow)
            self.assertNotEqual(p.returncode, 0)
            self.assertIsNone(json.loads(p.stdout)["window_radius"])

    def test_each_probe_gets_its_own_temp_and_out_and_cfg_flags_stay_last(self):
        """Measured on HIVE: a step-1b DIA-NN that never logged a radius wrote report.parquet
        into the chain's output directory -- the path step 5's output check reads -- and with no
        --temp a completed probe writes its .quant next to the raw. The cfg flags passed after
        `--` must still be the tail of DIA-NN's argv, exactly as in steps 2-5."""
        with tempfile.TemporaryDirectory() as d:
            argv_log = os.path.join(d, "argv.txt")
            fake = FAKE_DIANN.replace("#!/bin/bash\n",
                                      f'#!/bin/bash\nprintf "%s\\n" "$*" >> {argv_log}\n', 1)
            raws = [_raw(d, f"r{i}_w7.raw", (10 + i) * GB // 10) for i in range(3)]
            p, _, _ = _run_probe(d, raws, "--", "--qvalue", "0.01", "--var-mod",
                                 "UniMod:35,15.994915,M", env_extra=FAKE_DOTNET_ENV, fake=fake)
            self.assertEqual(p.returncode, 0, p.stderr)
            calls = [ln.split() for ln in _read(argv_log).splitlines()]
            self.assertEqual(len(calls), 3)
            work = os.path.join(d, "probe_work")
            temps = set()
            for argv in calls:
                temp = argv[argv.index("--temp") + 1]
                outp = argv[argv.index("--out") + 1]
                self.assertTrue(temp.startswith(work + os.sep), temp)
                self.assertTrue(outp.startswith(temp + os.sep), outp)
                temps.add(temp)
                k = argv.index("--threads") + 2
                self.assertEqual(argv[k:], ["--qvalue", "0.01", "--var-mod",
                                            "UniMod:35,15.994915,M"])
            self.assertEqual(len(temps), 3, "probes shared a --temp")

    def test_raw_list_is_accepted(self):
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"r{i}_w7.raw", (10 + i) * GB // 10) for i in range(4)]
            lst = os.path.join(d, "file list.txt")
            with open(lst, "w") as fh:
                fh.write("\n".join(raws) + "\n")
            fasta, lib = _fasta_lib(d)
            p = subprocess.run([sys.executable, PROBE, "--diann",
                                _exe(os.path.join(d, "diann"), FAKE_DIANN), "--raw-list", lst,
                                "--fasta", fasta, "--lib", lib, "--threads", "8",
                                "--workdir", os.path.join(d, "w")],
                               capture_output=True, text=True, timeout=240,
                               env=dict(os.environ, **FAKE_DOTNET_ENV))
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(len(json.loads(p.stdout)["probes"]), 3)

    def test_damaged_runs_are_named_in_the_job_log_because_they_are_still_searched(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_bruker(d, f"ok{i}_w7.d", gb=0.02 + i / 1000) for i in range(3)]
            trunc = _bruker(d, "trunc_w3.d", gb=0.03, nframes=140, indexed_frames=1)
            p, _, _ = _run_probe(d, runs + [trunc])
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["window_radius"], 7)
            self.assertNotIn("trunc_w3.d", _names(out["probes"]))
            self.assertEqual(_names(out["selection"]["excluded_damaged"]), ["trunc_w3.d"])
            self.assertRegex(p.stderr, r"WARNING.*trunc_w3\.d")

    def test_one_measured_radius_never_claims_the_radii_agree(self):
        """With two runs there is only ever ONE probe, and `radii_agree: true` on it asserted a
        corroboration that was never measured: radii [12], radii_agree true, n_probes 1."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, "a_w12.raw", 11 * GB // 10), _raw(d, "b_w12.raw", 12 * GB // 10)]
            p, _, _ = _run_probe(d, raws, env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["radii"], [12])
            self.assertEqual(len(out["probes"]), 1)
            self.assertIsNone(out["radii_agree"])
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"s{i}_w7.raw", (11 + i) * GB // 10) for i in range(3)]
            p, _, _ = _run_probe(d, raws, env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["radii"], [7, 7, 7])
            self.assertTrue(out["radii_agree"], "three runs that agreed are not recorded as such")

    def test_radii_that_disagree_by_more_than_two_are_called_out_in_the_job_log(self):
        """radii [7, 13] pinned 7 with only `radii_agree: false` buried in window.json."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, "a_w7.raw", 11 * GB // 10), _raw(d, "b_w13.raw", 12 * GB // 10),
                    _raw(d, "c_w7.raw", 13 * GB // 10)]
            p, _, _ = _run_probe(d, raws, env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(sorted(out["radii"]), [7, 7, 13])
            self.assertEqual(out["window_radius"], 7)
            self.assertFalse(out["radii_agree"])
            self.assertRegex(p.stderr, r"WARNING: the measured radii disagree by more than 2")

    def test_a_cohort_entirely_in_wal_mode_is_measured_rather_than_failing_the_chain(self):
        """Verified: --raw <one WAL-mode .d> gave rc 1, stopped_because no_probeable_run, so a
        whole cohort in WAL mode failed step 1b and left steps 2-5 in DependencyNeverSatisfied --
        over runs those steps would have searched anyway."""
        with tempfile.TemporaryDirectory() as d:
            runs = [_bruker(d, f"w{i}_w7.d", minutes=60, gb=0.02 + i / 1000) for i in range(6)]
            for r in runs:
                _set_wal_header(r)
            p, cfg, _ = _run_probe(d, runs)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["window_radius"], 7)
            self.assertEqual(len(out["probes"]), 3)
            self.assertEqual(_names(out["selection"]["probed_despite_wal_header"]),
                             [f"w{i}_w7.d" for i in range(6)])
            self.assertRegex(p.stderr, r"WARNING: probing w\d_w7\.d anyway")
            self.assertIn("--window 7", _read(cfg))

    def test_the_failure_does_not_promise_that_raw_overrides_the_checks(self):
        """`--raw` goes through the same filter, so "re-run with --raw naming runs to use" sent
        the reader in a circle -- --raw <a refused run> is refused again."""
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"r{i}_nowin.raw", (10 + i) * GB // 10) for i in range(4)]
            p, _, _ = _run_probe(d, raws, "--max-failures", "3", env_extra=FAKE_DOTNET_ENV)
            self.assertNotEqual(p.returncode, 0)
            self.assertNotIn("re-run with --raw naming runs to use", p.stderr)
            self.assertIn("--raw is not a way past it", p.stderr)

    def test_the_job_log_warns_that_a_mixed_cohort_has_no_single_radius(self):
        with tempfile.TemporaryDirectory() as d:
            raws = [_raw(d, f"Ex_{i}_w7.raw", (20 + i) * MB) for i in range(3)]
            dotds = [_bruker(d, f"tt{i}_w9.d", minutes=60, gb=0.1) for i in range(3)]
            p, _, _ = _run_probe(d, raws + dotds, env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertRegex(p.stderr, r"WARNING: cohort mixes Bruker \.d and Thermo \.raw")
            self.assertNotRegex(p.stderr, r"not probing Ex_\d_w7\.raw")
            self.assertEqual(json.loads(p.stdout)["selection"]["excluded_small"], [])

    def test_the_job_log_says_when_the_median_run_looks_like_a_wash(self):
        """12 washes and 2 samples probes wash07, wash03 and wash10 and pins a wash's radius.
        Nothing here can tell them apart -- but it must not do it quietly."""
        with tempfile.TemporaryDirectory() as d:
            washes = [_raw(d, f"wash{i:02d}_w3.raw", 15 * MB + i) for i in range(12)]
            samples = [_raw(d, f"sample{i}_w7.raw", 200 * MB + i) for i in range(2)]
            p, _, _ = _run_probe(d, washes + samples, env_extra=FAKE_DOTNET_ENV)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(_names(out["probes"]),
                             ["wash07_w3.raw", "wash03_w3.raw", "wash10_w3.raw"])
            self.assertEqual(out["window_radius"], 3)
            self.assertRegex(p.stderr, r"WARNING: the median run wash07_w3\.raw")

    def test_runs_on_python_3_8(self):
        """Step 1b runs `python3` on the COMPUTE node, not the interpreter that generated the
        chain; Ubuntu 20.04's is 3.8. No 3.9+ syntax, and no argparse.BooleanOptionalAction."""
        ast.parse(_read(PROBE), feature_version=(3, 8))
        code = ("import argparse, runpy, sys\n"
                "if hasattr(argparse, 'BooleanOptionalAction'): del argparse.BooleanOptionalAction\n"
                "script = sys.argv[1]; sys.argv = [script, '--help']; "
                "runpy.run_path(script, run_name='__main__')")
        p = subprocess.run([sys.executable, "-c", code, PROBE],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("--budget", p.stdout)


# --------------------------------------------------------------------------------------------
class ProbeLifecycleTests(unittest.TestCase):
    """What the probe does to the processes it starts, and what the job log shows."""

    # A wrapper that does NOT exec: the engine is this shell's child, as under apptainer.
    WRAPPER = """#!/bin/bash
bash -c 'echo started > "$1/started"; sleep 3; echo alive > "$1/alive"; echo "[0:04] Scan window radius set to 9"; sleep 20' _ "{d}"
"""

    def _argv(self, d, wrapper_body, timeout):
        fasta, lib = _fasta_lib(d)
        raw = _raw(d, "only_w9.raw", GB)
        diann = _exe(os.path.join(d, "diann_wrapper"), wrapper_body.format(d=d))
        cfg = os.path.join(d, "resolved.cfg")
        open(cfg, "w").close()
        argv = [sys.executable, PROBE, "--diann", diann, "--raw", raw, "--fasta", fasta,
                "--lib", lib, "--timeout", str(timeout), "--write-cfg", cfg,
                "--workdir", os.path.join(d, "w")]
        return argv, cfg, dict(os.environ, **FAKE_DOTNET_ENV)

    def _wait_for(self, path, seconds):
        end = time.time() + seconds
        while time.time() < end and not os.path.exists(path):
            time.sleep(0.05)
        return os.path.exists(path)

    def test_a_timeout_stops_the_engine_behind_a_wrapper_that_does_not_exec_it(self):
        with tempfile.TemporaryDirectory() as d:
            argv, cfg, env = self._argv(d, self.WRAPPER, timeout=1)
            t0 = time.time()
            p = subprocess.run(argv, capture_output=True, text=True, env=env, timeout=120)
            elapsed = time.time() - t0
            time.sleep(max(0.0, 4.0 - elapsed))
            self.assertFalse(os.path.exists(os.path.join(d, "alive")),
                             "the engine outlived the probe's timeout")
            self.assertLess(elapsed, 3.0, "the probe waited for the engine, not its timeout")
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIsNone(json.loads(p.stdout)["window_radius"])
            self.assertNotIn("--window", _read(cfg))

    def test_a_signal_to_the_probe_stops_the_engine_too(self):
        """DIA-NN runs in its own process group so a stop reaches all of it -- which also takes
        it out of the job's group. SIGTERM and Ctrl-C to the probe must still reach it."""
        with tempfile.TemporaryDirectory() as d:
            argv, cfg, env = self._argv(d, self.WRAPPER, timeout=60)
            p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, env=env)
            t_started = time.time()
            try:
                self.assertTrue(self._wait_for(os.path.join(d, "started"), 30))
                t_started = time.time()
                p.send_signal(signal.SIGTERM)
                p.communicate(timeout=30)
            finally:
                if p.poll() is None:
                    p.kill()
            time.sleep(max(0.0, 4.0 - (time.time() - t_started)))
            self.assertFalse(os.path.exists(os.path.join(d, "alive")),
                             "the engine kept running after the probe was terminated")
            self.assertNotEqual(p.returncode, 0)
            self.assertNotIn("--window", _read(cfg))

    def test_a_signal_during_selection_still_writes_the_evidence(self):
        """The handlers went on AFTER select_representative(), which stats every run in the
        cohort -- thousands of files over the cluster's storage. A SLURM kill in there ended the
        probe on the default disposition and window.json, the chain's only evidence, was never
        written."""
        with tempfile.TemporaryDirectory() as d:
            marker = os.path.join(d, "selecting")
            fasta, lib = _fasta_lib(d)
            raw = _raw(d, "only_w9.raw", GB)
            diann = _exe(os.path.join(d, "diann"), FAKE_DIANN)
            code = ("import sys, time\n"
                    "sys.path.insert(0, %r)\n"
                    "import probe_window\n"
                    "real = probe_window.measure_run\n"
                    "def slow(path, **kw):\n"
                    "    open(%r, 'w').close()\n"
                    "    time.sleep(5)\n"
                    "    return real(path, **kw)\n"
                    "probe_window.measure_run = slow\n"
                    "sys.argv = sys.argv[1:]\n"
                    "probe_window.main()\n") % (SCRIPTS, marker)
            argv = [sys.executable, "-c", code, "probe_window.py", "--diann", diann,
                    "--raw", raw, "--fasta", fasta, "--lib", lib,
                    "--workdir", os.path.join(d, "w")]
            p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                                 env=dict(os.environ, **FAKE_DOTNET_ENV))
            try:
                self.assertTrue(self._wait_for(marker, 30), "selection never started")
                p.send_signal(signal.SIGTERM)
                out, err = p.communicate(timeout=60)
            finally:
                if p.poll() is None:
                    p.kill()
            self.assertEqual(p.returncode, 128 + signal.SIGTERM, out + err)
            w = json.loads(out)
            self.assertIsNone(w["window_radius"])
            self.assertEqual(w["stopped_because"], "signal")

    def test_each_probe_is_announced_and_diann_output_reaches_the_job_log_as_it_runs(self):
        """A timsTOF step 1b took 1083 s with its job log empty until the end -- past
        watch_run.sh's 15-minute stall rule, whose playbook is to scancel and drop the file."""
        slow = """#!/bin/bash
echo "[0:00] Loading run $2"
sleep 20
"""
        with tempfile.TemporaryDirectory() as d:
            argv, _, env = self._argv(d, slow, timeout=60)
            p = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                 text=True, env=env)
            seen = []
            reader = threading.Thread(target=lambda: seen.extend(p.stderr), daemon=True)
            reader.start()
            try:
                end = time.time() + 15
                while time.time() < end:
                    txt = "".join(seen)
                    if "probe 1/1" in txt and "Loading run" in txt:
                        break
                    time.sleep(0.1)
                self.assertIsNone(p.poll(), "the probe ended before the test could look")
            finally:
                p.terminate()
                try:
                    p.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    p.kill()
            txt = "".join(seen)
            self.assertRegex(txt, r"probe 1/1: only_w9\.raw")
            self.assertIn("Loading run", txt, "DIA-NN's output was not in the job log while it ran")


# --------------------------------------------------------------------------------------------
class Step1bChainTests(unittest.TestCase):
    """The generated step1b_window.sbatch, executed with bash exactly as a compute node would
    run it -- a login-shell environment with no .NET on it."""

    def _chain(self, d, raws, dotnet_root=None):
        fasta = os.path.join(d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDER\n")
        cfg = os.path.join(d, "p.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--mass-acc 20\n--mass-acc-ms1 7\n")
        diann = _exe(os.path.join(d, "diann"), FAKE_DIANN)
        out = os.path.join(d, "out")
        env = {k: v for k, v in os.environ.items() if k != "DOTNET_ROOT"}
        if dotnet_root:
            env["PROTEOMICS_DOTNET_DIR"] = dotnet_root
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", diann, "--raw", *raws, "--fasta", fasta,
                            "--out", out, "--cfg", cfg, "--threads-per-file", "12"],
                           capture_output=True, text=True, env=env, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(out, "step1.predicted.speclib"), "w") as fh:
            fh.write("lib")                           # stands in for step 1's library
        return out, env, json.loads(p.stdout)

    def _run_step1b(self, out, env):
        env = {k: v for k, v in env.items() if k not in ("DOTNET_ROOT", "PROTEOMICS_DOTNET_DIR")}
        return subprocess.run(["bash", os.path.join(out, "step1b_window.sbatch")],
                              cwd=out, capture_output=True, text=True, env=env, timeout=240)

    def test_step1b_on_thermo_raw_carries_the_dotnet_prefix_and_measures(self):
        with tempfile.TemporaryDirectory() as d:
            root = _fake_dotnet_root(d)
            raws = [_raw(d, f"Ex_{i}_w8.raw", (12 + i) * GB // 10) for i in range(4)]
            out, env, _ = self._chain(d, raws, dotnet_root=root)
            body = _read(os.path.join(out, "step1b_window.sbatch"))
            self.assertIn(f"export DOTNET_ROOT={root}", body)
            self.assertLess(body.index("export DOTNET_ROOT"), body.index("probe_window.py"))
            p = self._run_step1b(out, env)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(_read(os.path.join(out, "window.txt")).strip(), "8")

    def test_step1b_probes_representative_runs_not_the_first_file(self):
        with tempfile.TemporaryDirectory() as d:
            root = _fake_dotnet_root(d)
            blank = _raw(d, "000_blank_w3.raw", 20 * MB)
            runs = [_raw(d, "Ex_a_w7.raw", 12 * GB // 10), _raw(d, "Ex_b_w8.raw", 13 * GB // 10),
                    _raw(d, "Ex_c_w8.raw", 14 * GB // 10), _raw(d, "Ex_d_w9.raw", 15 * GB // 10)]
            out, env, _ = self._chain(d, [blank] + runs, dotnet_root=root)
            p = self._run_step1b(out, env)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            w = json.load(open(os.path.join(out, "window.json")))
            probed = _names(w["probes"])
            self.assertNotIn("000_blank_w3.raw", probed)
            self.assertEqual(len(probed), 3)
            self.assertEqual(_read(os.path.join(out, "window.txt")).strip(),
                             str(w["window_radius"]))
            self.assertNotEqual(w["window_radius"], 3, "the blank's radius was pinned")
            resolved = dp.cfg_tokens(os.path.join(out, "params.resolved.cfg"))
            self.assertEqual(resolved.count("--window"), 1)
            self.assertEqual(resolved[resolved.index("--window") + 1], str(w["window_radius"]))
            # every run's radius is in the job log, not only in window.json
            for x in w["probes"]:
                self.assertRegex(p.stdout + p.stderr,
                                 re.escape(os.path.basename(x["file"])) + r".*radius %d" % x["radius"])

    def test_step1b_never_probes_a_truncated_timstof_run(self):
        """The pilot: the truncated run had the largest tdf_bin, so a bytes rule chose it."""
        with tempfile.TemporaryDirectory() as d:
            runs = [_bruker(d, f"A{i}_w7.d", minutes=21, gb=0.02 + i / 1000) for i in range(4)]
            trunc = _bruker(d, "Koganti_trunc_w3.d", minutes=21, gb=0.03, nframes=140,
                            indexed_frames=1)
            out, env, _ = self._chain(d, [trunc] + runs)
            p = self._run_step1b(out, env)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            w = json.load(open(os.path.join(out, "window.json")))
            self.assertNotIn("Koganti_trunc_w3.d", _names(w["probes"]))
            self.assertEqual(_names(w["selection"]["excluded_damaged"]), ["Koganti_trunc_w3.d"])
            self.assertEqual(_read(os.path.join(out, "window.txt")).strip(), "7")
            self.assertRegex(p.stderr, r"WARNING.*Koganti_trunc_w3\.d")

    def test_step1b_timeout_is_per_probe_and_the_budget_fits_the_job(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_raw(d, f"run{i}_w7.mzML", (30 + i) * MB) for i in range(3)]
            out, _, _ = self._chain(d, runs)
            body = _read(os.path.join(out, "step1b_window.sbatch"))
            hours = int(re.search(r"#SBATCH --time=(\d+):00:00", body).group(1))
            self.assertEqual(hours, dp.PROBE_WALL_HOURS)
            self.assertIn(f"--timeout {dp.PROBE_TIMEOUT_S}", body)
            budget = int(re.search(r"--budget (\d+)", body).group(1))
            self.assertLess(budget, hours * 3600, "probes could outlive the job's wall clock")
            self.assertGreaterEqual(budget, dp.PROBE_CANDIDATES * dp.PROBE_TIMEOUT_S,
                                    "three probes at their full timeout no longer fit")
            self.assertIn("--raw-list", body)

    def test_step1b_probes_in_its_own_workdir_under_out(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_raw(d, f"run{i}_w7.mzML", (30 + i) * MB) for i in range(3)]
            out, env, _ = self._chain(d, runs)
            p = self._run_step1b(out, env)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            logs = glob.glob(os.path.join(out, "window_probe", "probe*", "probe.log"))
            self.assertEqual(len(logs), 3, logs)
            for x in json.load(open(os.path.join(out, "window.json")))["probes"]:
                self.assertTrue(x["log"].startswith(os.path.join(out, "window_probe") + os.sep))

    def test_a_failed_step1b_stops_cleanly_keeps_its_evidence_and_leaves_no_radius(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_raw(d, f"run{i}_nowin.mzML", (30 + i) * MB) for i in range(5)]
            out, env, _ = self._chain(d, runs)
            with open(os.path.join(out, "window.txt"), "w") as fh:
                fh.write("7\n")                    # from an earlier run of this chain
            p = self._run_step1b(out, env)
            log = p.stdout + p.stderr
            self.assertNotEqual(p.returncode, 0, log)
            self.assertNotIn("Traceback", log)
            self.assertIn("FAILED", log)
            self.assertIn("DependencyNeverSatisfied", log)
            self.assertFalse(os.path.exists(os.path.join(out, "window.txt")))
            self.assertFalse(os.path.exists(os.path.join(out, "params.resolved.cfg")))
            w = json.load(open(os.path.join(out, "window.json")))
            self.assertIsNone(w["window_radius"])
            self.assertEqual(len(w["probes"]), dp.PROBE_MAX_FAILURES)

    def test_timstof_only_step1b_has_no_dotnet_export(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_bruker(d, f"run{i}_w7.d", gb=0.02 + i / 1000) for i in range(3)]
            out, env, _ = self._chain(d, runs)
            body = _read(os.path.join(out, "step1b_window.sbatch"))
            self.assertNotIn("DOTNET_ROOT", body)
            p = self._run_step1b(out, env)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(_read(os.path.join(out, "window.txt")).strip(), "7")

    def test_provenance_points_at_the_evidence_not_at_the_first_files(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_raw(d, f"run{i}_w7.mzML", (30 + i) * MB) for i in range(6)]
            out, _, info = self._chain(d, runs)
            sw = info["scan_window"]
            self.assertEqual(sw["evidence_file"], os.path.join(out, "window.json"))
            self.assertNotIn("probe_candidates", sw)
            self.assertRegex(sw["probe_rule"], r"median")


if __name__ == "__main__":
    unittest.main(verbosity=2)
