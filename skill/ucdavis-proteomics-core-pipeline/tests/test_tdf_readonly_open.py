#!/usr/bin/env python3
"""
A Bruker analysis.tdf is opened so nothing can write to it or replay a stale -wal,
and a truncated one is reported before anything searches it.

Written for the Hive incident (2026-09-16): 342 of 39,374 Bruker .d had an
analysis.tdf whose frame index covered only part of a complete analysis.tdf_bin.
Every one was still in WAL mode (SQLite header bytes 18-19 = 2,2; intact files are
1,1). The instrument->storage copy can leave a stale, mid-acquisition
analysis.tdf-wal beside the finished file, and a READ-WRITE sqlite open replays it:
it checkpoints the mid-run pages into the finished file and truncates it. One run's
index ended up covering 0.7% of its 2.4 GB tdf_bin, and a search read 121 cycles of
it without an error anywhere.

`mode=ro` alone does not truncate, but it still READS the stale -wal (so the reader
sees the mid-run database) and drops -shm files beside the tdf. `immutable=1` reads
the file as it is on disk and never writes. The fixtures below reproduce the
mechanism with sqlite itself, and each hazard test first proves the fixture really
has the hazard, so the immutable tests cannot pass vacuously.

Stdlib only (unittest + sqlite3), like the rest of this suite.
"""
import hashlib
import json
import os
import re
import shutil
import sqlite3
import struct
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import bruker_tdf                           # noqa: E402
import detect_acquisition as da             # noqa: E402
import make_methods                         # noqa: E402

BLOCK = 64                                  # bytes per synthetic frame block in tdf_bin
N_FRAMES = 1000
N_STALE = 10                                # frames the stale mid-acquisition -wal knows about


def _write_bin(d, n_blocks):
    """analysis.tdf_bin: one block per frame, each starting with its uint32 block size
    (then a uint32 scan count), as in a real TDF."""
    with open(os.path.join(d, "analysis.tdf_bin"), "wb") as fh:
        for _ in range(n_blocks):
            fh.write(struct.pack("<II", BLOCK, 5) + b"\0" * (BLOCK - 8))


def _schema(con):
    con.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
    con.executemany("INSERT INTO GlobalMetadata VALUES (?,?)",
                    [("InstrumentName", "timsTOF HT"), ("AcquisitionSoftware", "timsControl")])
    con.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, MsMsType INTEGER, "
                "TimsId INTEGER, NumScans INTEGER, AccumulationTime REAL, RampTime REAL)")
    con.execute("CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER, WindowGroup INTEGER)")
    con.execute("CREATE TABLE DiaFrameMsMsWindowGroups (Id INTEGER)")
    con.execute("CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, "
                "ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, "
                "CollisionEnergy REAL)")


def _frames(con, lo, hi):
    con.executemany("INSERT INTO Frames VALUES (?,?,?,?,?,?,?)",
                    [(i + 1, i * 0.1, 9 if i % 2 else 0, i * BLOCK, 5, 100.0, 100.0)
                     for i in range(lo, hi)])


def _windows(con, lo=299.5, hi=1200.5, width=25.0, n=6):
    step = (hi - lo - width) / (n - 1)
    con.executemany("INSERT INTO DiaFrameMsMsWindows VALUES (?,?,?,?,?,?)",
                    [(1, 0, 100, lo + width / 2 + i * step, width, 30.0) for i in range(n)])


def make_intact_d(tmp, name="intact.d", n_indexed=N_FRAMES, n_bin=N_FRAMES, wal=False):
    """A finished run in rollback-journal mode (header 1,1) unless `wal`. `n_indexed` <
    `n_bin` gives a truncated index; `n_bin` < `n_indexed` an incomplete tdf_bin."""
    d = os.path.join(tmp, name)
    os.makedirs(d)
    _write_bin(d, n_bin)
    con = sqlite3.connect(os.path.join(d, "analysis.tdf"))
    if wal:
        con.execute("PRAGMA journal_mode=WAL")
    _schema(con)
    _frames(con, 0, n_indexed)
    _windows(con)
    con.commit()
    con.close()
    return d


def make_stale_wal_d(tmp, name="stale.d"):
    """A finished run with a stale mid-acquisition -wal beside it -- the Hive state.

    The -wal is captured after N_STALE frames; the acquisition then finishes, is
    checkpointed into the main file and closed (which deletes the live -wal but leaves
    the header in WAL mode); the stale -wal is then put back. The DIA windows are
    written after the snapshot, so a reader that replays the -wal cannot see them."""
    d = os.path.join(tmp, name)
    os.makedirs(d)
    _write_bin(d, N_FRAMES)
    tdf = os.path.join(d, "analysis.tdf")
    con = sqlite3.connect(tdf)
    con.execute("PRAGMA journal_mode=WAL")
    con.execute("PRAGMA wal_autocheckpoint=0")
    _schema(con)
    _frames(con, 0, N_STALE)
    con.commit()
    with open(tdf + "-wal", "rb") as fh:
        stale = fh.read()
    _frames(con, N_STALE, N_FRAMES)
    _windows(con)
    con.commit()
    con.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    con.close()
    with open(tdf + "-wal", "wb") as fh:
        fh.write(stale)
    return d


def truncate_like_hive(d):
    """What happened on Hive: one plain read-write open of the tdf, then close."""
    con = sqlite3.connect(os.path.join(d, "analysis.tdf"))
    con.execute("SELECT COUNT(*) FROM Frames").fetchone()
    con.close()
    return d


def snapshot(d):
    """Names and content hashes of everything in the .d, to prove a reader wrote nothing."""
    out = {}
    for f in sorted(os.listdir(d)):
        with open(os.path.join(d, f), "rb") as fh:
            out[f] = hashlib.sha256(fh.read()).hexdigest()
    return out


def frame_count(tdf, query):
    con = sqlite3.connect(f"file:{tdf}?{query}", uri=True)
    try:
        return con.execute("SELECT COUNT(*) FROM Frames").fetchone()[0]
    finally:
        con.close()


class TestFixturesReallyHaveTheHazard(unittest.TestCase):
    """Without these, every 'immutable reads the finished file' test below could pass
    against a fixture that never had a stale -wal in it."""

    def test_mode_ro_alone_replays_the_stale_wal_and_writes_a_shm(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            tdf = os.path.join(d, "analysis.tdf")
            self.assertEqual(frame_count(tdf, "mode=ro"), N_STALE)
            self.assertIn("analysis.tdf-shm", os.listdir(d))

    def test_a_read_write_open_truncates_the_finished_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            tdf = os.path.join(d, "analysis.tdf")
            size_before = os.path.getsize(tdf)
            truncate_like_hive(d)
            self.assertLess(os.path.getsize(tdf), size_before)
            self.assertEqual(frame_count(tdf, "mode=ro&immutable=1"), N_STALE)
            self.assertNotIn("analysis.tdf-wal", os.listdir(d))
            with open(tdf, "rb") as fh:
                hdr = fh.read(20)
            self.assertEqual((hdr[18], hdr[19]), (2, 2), "the header stays in WAL mode")


class TestConnectTdf(unittest.TestCase):
    def test_reads_the_finished_file_not_the_stale_wal(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            con = bruker_tdf.connect_tdf(os.path.join(d, "analysis.tdf"))
            try:
                self.assertEqual(con.execute("SELECT COUNT(*) FROM Frames").fetchone()[0],
                                 N_FRAMES)
            finally:
                con.close()

    def test_writes_nothing_beside_the_tdf(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            before = snapshot(d)
            con = bruker_tdf.connect_tdf(os.path.join(d, "analysis.tdf"))
            con.execute("SELECT * FROM Frames").fetchall()
            con.close()
            self.assertEqual(snapshot(d), before)

    def test_refuses_writes(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            con = bruker_tdf.connect_tdf(os.path.join(d, "analysis.tdf"))
            try:
                with self.assertRaises(sqlite3.OperationalError):
                    con.execute("CREATE TABLE scribble (x)")
            finally:
                con.close()

    def test_uri_is_read_only_and_immutable(self):
        uri = bruker_tdf.tdf_uri("/data/run.d/analysis.tdf")
        self.assertTrue(uri.startswith("file:"))
        self.assertIn("mode=ro", uri)
        self.assertIn("immutable=1", uri)

    def test_uri_survives_characters_that_end_a_uri_path(self):
        # `?` and `#` end the path of a URI and `%` starts an escape; spliced in raw, the
        # open would silently target a different file name.
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp, name="odd #1? 50%.d")
            con = bruker_tdf.connect_tdf(os.path.join(d, "analysis.tdf"))
            try:
                self.assertEqual(con.execute("SELECT COUNT(*) FROM Frames").fetchone()[0],
                                 N_FRAMES)
            finally:
                con.close()

    def test_missing_file_is_an_error_not_a_new_empty_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(sqlite3.OperationalError):
                bruker_tdf.connect_tdf(os.path.join(tmp, "nope.d", "analysis.tdf"))
            self.assertEqual(os.listdir(tmp), [])


class TestTdfIntegrity(unittest.TestCase):
    def test_intact_rollback_file_is_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_intact_d(tmp))
            self.assertEqual(r["status"], "ok", r)
            self.assertEqual(r["problems"], [])
            self.assertFalse(r["sqlite_header_wal"])
            self.assertAlmostEqual(r["index_coverage"], 1.0)
            self.assertEqual(r["index_end_bytes"], N_FRAMES * BLOCK)
            self.assertEqual(r["tdf_bin_bytes"], N_FRAMES * BLOCK)

    def test_truncated_by_a_read_write_open_is_truncated(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(truncate_like_hive(make_stale_wal_d(tmp)))
            self.assertEqual(r["status"], "truncated", r)
            self.assertAlmostEqual(r["index_coverage"], N_STALE / N_FRAMES)
            self.assertTrue(r["sqlite_header_wal"])
            self.assertTrue(any("1.0%" in p for p in r["problems"]), r["problems"])

    def test_truncated_index_is_caught_without_any_wal_signal(self):
        # the coverage check stands on its own: header 1,1, no -wal, index 1% of tdf_bin
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_intact_d(tmp, n_indexed=N_STALE))
            self.assertEqual(r["status"], "truncated", r)
            self.assertFalse(r["sqlite_header_wal"])
            self.assertEqual(r["wal_bytes"], 0)

    def test_stale_wal_beside_a_complete_index_is_at_risk(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            before = snapshot(d)
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "at_risk", r)
            self.assertAlmostEqual(r["index_coverage"], 1.0, "read from the finished file")
            self.assertGreater(r["wal_bytes"], 0)
            self.assertTrue(any("read-write" in p for p in r["problems"]), r["problems"])
            self.assertEqual(snapshot(d), before, "the check itself must write nothing")

    def test_wal_header_alone_is_at_risk(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp, wal=True)
            self.assertNotIn("analysis.tdf-wal", os.listdir(d))
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "at_risk", r)
            self.assertTrue(r["sqlite_header_wal"])

    def test_non_empty_wal_alone_is_at_risk(self):
        # header 1,1, so only the -wal size can raise it
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            with open(os.path.join(d, "analysis.tdf-wal"), "wb") as fh:
                fh.write(b"\x37\x7f\x06\x82" + b"\0" * 28)
            r = bruker_tdf.tdf_integrity(d)
            self.assertFalse(r["sqlite_header_wal"])
            self.assertEqual(r["status"], "at_risk", r)
            self.assertEqual(r["wal_bytes"], 32)
            self.assertTrue(any("analysis.tdf-wal" in p for p in r["problems"]), r["problems"])

    def test_empty_wal_beside_a_rollback_file_is_not_a_problem(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            open(os.path.join(d, "analysis.tdf-wal"), "wb").close()
            self.assertEqual(bruker_tdf.tdf_integrity(d)["status"], "ok")

    def test_hot_rollback_journal_is_at_risk(self):
        # a read-write open would roll a non-empty -journal back into the file, too
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            with open(os.path.join(d, "analysis.tdf-journal"), "wb") as fh:
                fh.write(b"\xd9\xd5\x05\xf9\x20\xa1\x63\xd7" + b"\0" * 504)
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "at_risk", r)
            self.assertGreater(r["journal_bytes"], 0)

    def test_index_past_the_end_of_tdf_bin_is_bin_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_intact_d(tmp, n_bin=N_FRAMES // 2))
            self.assertEqual(r["status"], "bin_incomplete", r)

    def test_missing_tdf_bin_is_bin_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            os.remove(os.path.join(d, "analysis.tdf_bin"))
            self.assertEqual(bruker_tdf.tdf_integrity(d)["status"], "bin_incomplete")

    def test_unreadable_index_is_unverified_not_ok(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = os.path.join(tmp, "noframes.d")
            os.makedirs(d)
            _write_bin(d, 4)
            con = sqlite3.connect(os.path.join(d, "analysis.tdf"))
            con.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
            con.commit()
            con.close()
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "unverified", r)
            self.assertTrue(r["problems"])

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0, "root reads mode-000 files")
    def test_unreadable_tdf_bin_is_unverified_not_a_crash(self):
        # one unreadable file must not take down detection for the whole cohort
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            b = os.path.join(d, "analysis.tdf_bin")
            os.chmod(b, 0)
            try:
                r = bruker_tdf.tdf_integrity(d)
            finally:
                os.chmod(b, 0o644)
            self.assertEqual(r["status"], "unverified", r)
            self.assertTrue(any("analysis.tdf_bin" in p for p in r["problems"]), r["problems"])

    def test_not_an_sqlite_file_is_unverified(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = os.path.join(tmp, "junk.d")
            os.makedirs(d)
            _write_bin(d, 4)
            with open(os.path.join(d, "analysis.tdf"), "wb") as fh:
                fh.write(b"not a database" * 10)
            self.assertEqual(bruker_tdf.tdf_integrity(d)["status"], "unverified")

    def test_truncation_outranks_the_wal_signals(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp, n_indexed=N_STALE, wal=True)
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "truncated")
            # the finding that decides what to do is the one read first
            self.assertTrue(r["problems"][0].startswith("analysis.tdf is truncated"), r["problems"])
            self.assertEqual(len(r["problems"]), 2)


class TestDetectAcquisition(unittest.TestCase):
    """Warnings are per file (`files[i].warnings`) and echoed to stderr, and any warning sets
    needs_confirmation -- the same shape the Thermo .raw reader uses for its read failures,
    so a caller has one place to look."""

    def _run(self, *paths):
        res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                              *paths], capture_output=True, text=True)
        self.assertEqual(res.returncode, 0, res.stderr)
        return json.loads(res.stdout), res.stderr

    def _cli(self, *paths):
        return self._run(*paths)[0]

    def test_reads_the_finished_file_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            before = snapshot(d)
            out = self._cli(d)
            self.assertEqual(snapshot(d), before, "detect_acquisition wrote beside the tdf")
            # the DIA windows exist only in the finished file, not in the stale -wal
            self.assertEqual(out["precursor_mz_range"], [299.5, 1200.5])
            self.assertEqual(out["instrument"], "timsTOF HT")

    def test_intact_run_needs_no_confirmation(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, err = self._run(make_intact_d(tmp))
            f = out["files"][0]
            self.assertEqual(f["tdf_integrity"]["status"], "ok")
            self.assertEqual(f["warnings"], [])
            self.assertEqual(out["tdf_integrity_problem_files"], [])
            self.assertFalse(out["needs_confirmation"], out)
            self.assertNotIn("WARNING", err)

    def test_truncated_run_warns_and_needs_confirmation(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = make_intact_d(tmp)
            bad = truncate_like_hive(make_stale_wal_d(tmp, name="bad.d"))
            out, err = self._run(good, bad)
            self.assertTrue(out["needs_confirmation"])
            self.assertEqual(out["tdf_integrity_problem_files"], [bad])
            by_file = {f["file"]: f for f in out["files"]}
            self.assertEqual(by_file[good]["warnings"], [])
            self.assertEqual(len(by_file[bad]["warnings"]), 1, by_file[bad]["warnings"])
            self.assertIn("truncated", by_file[bad]["warnings"][0])
            self.assertIn("Do NOT search", by_file[bad]["warnings"][0])
            # a person reading the terminal sees it too, with the file named
            self.assertIn(f"WARNING: {bad}: ", err)
            # acquisition detection itself is unchanged -- the warning is what gates
            self.assertEqual(out["overall"], "DIA")
            self.assertEqual(out["low_confidence_files"], [])

    def test_at_risk_run_warns_and_needs_confirmation(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            out = self._cli(d)
            self.assertTrue(out["needs_confirmation"])
            f = out["files"][0]
            self.assertEqual(f["tdf_integrity"]["status"], "at_risk")
            self.assertIn("read-write", f["warnings"][0])

    def test_non_bruker_files_carry_no_tdf_check(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "run.mzML")
            open(p, "w").close()
            r = da.classify(p)
            self.assertIsNone(r["tdf_integrity"])
            self.assertIsInstance(r["warnings"], list)


class TestMakeMethods(unittest.TestCase):
    def test_bruker_meta_reads_the_finished_file_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            before = snapshot(d)
            m = make_methods.bruker_meta(d)
            self.assertEqual(snapshot(d), before, "make_methods wrote beside the tdf")
            self.assertNotIn("error", m)
            self.assertEqual(m["n_windows"], 6)
            self.assertEqual(m["mode"], "dia-PASEF")


class TestNoOtherTdfOpens(unittest.TestCase):
    def test_every_tdf_open_in_the_scripts_is_immutable(self):
        # `mode=ro` without `immutable=1` got into two scripts because each spelled its own
        # sqlite3.connect. Any script that reads an analysis.tdf must open it through
        # bruker_tdf.connect_tdf, or at least through a URI helper / literal that carries
        # immutable=1 -- never a bare path and never mode=ro alone.
        offenders = []
        for name in sorted(os.listdir(SCRIPTS)):
            if not name.endswith(".py") or name == "bruker_tdf.py":
                continue
            with open(os.path.join(SCRIPTS, name), encoding="utf-8") as fh:
                src = fh.read()
            if "analysis.tdf" not in src:
                continue
            for m in re.finditer(r"sqlite3\.connect\s*\((.*)", src):
                arg = m.group(1)
                if "immutable=1" not in arg and not re.search(r"\b\w*tdf_uri\s*\(", arg):
                    offenders.append(f"{name}: sqlite3.connect({arg.strip()}")
        self.assertEqual(offenders, [],
                         "open analysis.tdf with bruker_tdf.connect_tdf() (mode=ro&immutable=1)")

    def test_the_skill_scripts_here_use_the_shared_helper(self):
        for name in ("detect_acquisition.py", "make_methods.py"):
            with open(os.path.join(SCRIPTS, name), encoding="utf-8") as fh:
                src = fh.read()
            self.assertIn("from bruker_tdf import", src, name)
            self.assertNotRegex(src, r"sqlite3\.connect\s*\(", name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
