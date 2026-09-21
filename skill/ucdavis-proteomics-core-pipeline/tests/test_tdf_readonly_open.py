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


def _repo_root():
    """The top of the checkout, or this skill package when it is installed on its own.

    marketplace.json is what makes a directory the checkout: it is the file that decides
    which skills get installed, so everything below it is shipped code and everything below
    it gets scanned. Installed from the marketplace there is no repo above the package, and
    the package alone is then the whole of what shipped."""
    d = HERE
    for _ in range(8):
        d = os.path.dirname(d)
        if os.path.isfile(os.path.join(d, ".claude-plugin", "marketplace.json")):
            return d
    return os.path.dirname(HERE)


ROOT = _repo_root()

# Directories with no shipped source of ours in them.
SKIP_DIRS = {".git", ".github", ".pytest_cache", "__pycache__", ".mypy_cache", ".ruff_cache",
             "node_modules", ".venv", "venv", "env", "renv", "build", "dist", "site-packages",
             ".Rproj.user", ".quarto", "_snaps"}

# (parent directory, file name) of the files allowed to spell an unsafe tdf open, because
# writing one IS their job: bruker_tdf.py and helpers_instrument.R define the safe opens and
# document the unsafe ones beside them, and the fixture writers here have to CREATE a tdf
# (and, in this file, prove a mode=ro open really replays a stale -wal, or the immutable
# tests below would pass against a fixture that never had the hazard). Matching on the
# parent directory rather than a full path keeps the entry valid in the dev fork and in a
# standalone install, where the path above the package differs.
ALLOW_UNSAFE_PY = {("scripts", "bruker_tdf.py"),
                   ("tests", "test_tdf_readonly_open.py"),
                   ("tests", "test_precursor_mz_range.py")}
ALLOW_UNSAFE_R = {("R", "helpers_instrument.R"),
                  ("testthat", "test-helpers_instrument_tdf.R")}


def _sources(root, suffixes):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for name in sorted(filenames):
            if not name.endswith(suffixes):
                continue
            path = os.path.join(dirpath, name)
            try:
                with open(path, encoding="utf-8") as fh:
                    text = fh.read()
            except (OSError, UnicodeDecodeError):
                continue
            if "analysis.tdf" in text:       # only files that read a Bruker tdf
                yield path, text


def scanned_files(root):
    """Every file the two scans below look at -- what proves the walk is not empty."""
    return [p for p, _ in _sources(root, (".py", ".R", ".r"))]


def _offenders(root, suffixes, allow, call, safe):
    out = []
    for path, text in _sources(root, suffixes):
        if not path.endswith(suffixes):
            continue
        parent = os.path.basename(os.path.dirname(path))
        if (parent, os.path.basename(path)) in allow:
            continue
        rel = os.path.relpath(path, root)
        for m in re.finditer(call, text):
            arg = m.group(1)
            if not re.search(safe, arg):
                line = text.count("\n", 0, m.start()) + 1
                out.append(f"{rel}:{line}: {m.group(0).strip()}")
    return sorted(out)


def python_tdf_offenders(root):
    """Every sqlite3.connect of an analysis.tdf under `root` that is not immutable.

    `mode=ro` without `immutable=1` got into five scripts across two skills because each
    spelled its own sqlite3.connect. Any file that reads an analysis.tdf must open it
    through bruker_tdf.connect_tdf, or at least through a URI helper / literal that carries
    immutable=1 -- never a bare path and never mode=ro alone.

    The one act that legitimately writes is BUILDING a synthetic .d for a fixture (or
    damaging a throwaway copy of one on purpose, to prove the fixture carries the hazard).
    It says so by name: synthetic_tdf.synthetic_tdf_write_uri(). The exemption is that
    NAME, not the file it appears in -- it is granted per CALL, so a bare connect in the
    very same module is still an offender, and `grep` lists every deliberate writer there
    is. Widening ALLOW_UNSAFE_PY instead is how this guard stops being one: a per-directory
    allowlist is exactly what failed the first time round (TestNoOtherTdfOpens)."""
    return _offenders(root, (".py",), ALLOW_UNSAFE_PY,
                      r"sqlite3\.connect\s*\((.*)",
                      r"immutable=1|\b\w*tdf_uri\s*\(|\bsynthetic_tdf_write_uri\s*\(")


def r_tdf_offenders(root):
    """The same for RSQLite: every dbConnect of an analysis.tdf that is not immutable."""
    return _offenders(root, (".R", ".r"), ALLOW_UNSAFE_R,
                      r"dbConnect\s*\((.*)",
                      r"immutable=1|\bsqlite_immutable_uri\s*\(")

sys.path.insert(0, HERE)

import bruker_tdf                           # noqa: E402
from synthetic_tdf import synthetic_tdf_write_uri   # noqa: E402  (the deliberate writer)
import detect_acquisition as da             # noqa: E402
import make_methods                         # noqa: E402
import run_search                           # noqa: E402

# Bytes per synthetic frame block in tdf_bin. N_FRAMES * BLOCK has to clear
# bruker_tdf.SLACK_BYTES (64 KiB) by a wide margin, or every "this is truncated"
# fixture below would sit inside the absolute slack the check allows a tiny tdf_bin
# and read as intact. A real tdf_bin is hundreds of MB to GB, so 250 KiB here is
# still orders of magnitude on the small side.
BLOCK = 256
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


def make_stale_wal_d(tmp, name="stale.d", rollback_header=False):
    """A finished run with a stale mid-acquisition -wal beside it -- the Hive state.

    The -wal is captured after N_STALE frames; the acquisition then finishes, is
    checkpointed into the main file and closed (which deletes the live -wal but leaves
    the header in WAL mode); the stale -wal is then put back. The DIA windows are
    written after the snapshot, so a reader that replays the -wal cannot see them.

    `rollback_header` switches the finished file out of WAL mode before it is closed
    (header 1,1) -- the state of the two intact HIVE blank runs with a ~4 MB stale -wal."""
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
    if rollback_header:
        con.execute("PRAGMA journal_mode=DELETE")
    con.close()
    with open(tdf + "-wal", "wb") as fh:
        fh.write(stale)
    return d


def header_bytes(tdf):
    with open(tdf, "rb") as fh:
        hdr = fh.read(20)
    return hdr[18], hdr[19]


def move_side_files_out(d, backup):
    """The remedy a stale_side_file warning gives: copy the side files to a backup outside
    the .d, then remove them from the .d."""
    os.makedirs(backup, exist_ok=True)
    for side in ("-wal", "-shm", "-journal"):
        p = os.path.join(d, "analysis.tdf" + side)
        if os.path.exists(p):
            shutil.copy2(p, backup)
            os.remove(p)


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
            self.assertEqual(header_bytes(tdf), (2, 2), "the header stays in WAL mode")

    def test_a_rollback_header_does_not_stop_a_stale_wal_being_replayed(self):
        # SQLite opens WAL mode whenever a -wal exists, whatever header bytes 18-19 say. So a
        # finished 1,1 file with a stale -wal beside it (the HIVE blank runs) is read as the
        # mid-acquisition database by a read-only open, and truncated by a read-write one.
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp, rollback_header=True)
            tdf = os.path.join(d, "analysis.tdf")
            self.assertEqual(header_bytes(tdf), (1, 1))
            self.assertEqual(frame_count(tdf, "mode=ro&immutable=1"), N_FRAMES)
            self.assertEqual(frame_count(tdf, "mode=ro"), N_STALE)
            size_before = os.path.getsize(tdf)
            truncate_like_hive(d)
            self.assertLess(os.path.getsize(tdf), size_before)
            self.assertEqual(header_bytes(tdf), (2, 2), "the replay brings the WAL header back")
            self.assertEqual(frame_count(tdf, "mode=ro&immutable=1"), N_STALE)
            self.assertEqual(bruker_tdf.tdf_integrity(d)["status"], "truncated")

    def test_with_the_side_files_moved_out_every_open_reads_the_whole_run(self):
        # what the stale_side_file remedy relies on: the finished file indexes the whole run
        # on its own, and once nothing stale sits beside it no kind of open can change that
        for rollback_header in (True, False):
            with self.subTest(rollback_header=rollback_header), \
                    tempfile.TemporaryDirectory() as tmp:
                d = make_stale_wal_d(tmp, rollback_header=rollback_header)
                tdf = os.path.join(d, "analysis.tdf")
                frame_count(tdf, "mode=ro")         # leaves an -shm beside it, as on HIVE
                move_side_files_out(d, os.path.join(tmp, "backup"))
                with open(tdf, "rb") as fh:
                    finished = fh.read()
                self.assertEqual(frame_count(tdf, "mode=ro"), N_FRAMES)
                truncate_like_hive(d)
                with open(tdf, "rb") as fh:
                    self.assertEqual(fh.read(), finished, "a read-write open changed the file")
                self.assertEqual(frame_count(tdf, "mode=ro&immutable=1"), N_FRAMES)
                self.assertIn(bruker_tdf.tdf_integrity(d)["status"], ("ok", "at_risk"))

    def test_a_wal_header_alone_is_read_whole_by_every_open(self):
        # header 2,2 and nothing beside it: nothing to replay, so it is safe to search
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp, wal=True)
            tdf = os.path.join(d, "analysis.tdf")
            with open(tdf, "rb") as fh:
                finished = fh.read()
            self.assertEqual(frame_count(tdf, "mode=ro"), N_FRAMES)
            for side in ("-wal", "-shm"):             # mode=ro leaves them; both are harmless
                if os.path.exists(tdf + side):
                    os.remove(tdf + side)
            truncate_like_hive(d)
            with open(tdf, "rb") as fh:
                self.assertEqual(fh.read(), finished, "a read-write open changed the file")
            self.assertEqual(sorted(os.listdir(d)), ["analysis.tdf", "analysis.tdf_bin"])


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
            # the shortfall in bytes, so a few KB of trailing slack on a tiny run is recognisable
            short = N_FRAMES * BLOCK - N_STALE * BLOCK
            self.assertTrue(any(f"{short:,} bytes short" in p for p in r["problems"]),
                            r["problems"])

    def test_truncated_index_is_caught_without_any_wal_signal(self):
        # the coverage check stands on its own: header 1,1, no -wal, index 1% of tdf_bin
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_intact_d(tmp, n_indexed=N_STALE))
            self.assertEqual(r["status"], "truncated", r)
            self.assertFalse(r["sqlite_header_wal"])
            self.assertEqual(r["wal_bytes"], 0)

    def test_stale_wal_beside_a_complete_wal_mode_index_is_stale_side_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp)
            before = snapshot(d)
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "stale_side_file", r)
            self.assertAlmostEqual(r["index_coverage"], 1.0, "read from the finished file")
            self.assertGreater(r["wal_bytes"], 0)
            self.assertEqual(snapshot(d), before, "the check itself must write nothing")

    def test_stale_wal_beside_a_rollback_header_is_stale_side_file(self):
        # the two intact HIVE blank runs: header 1,1, complete index, ~4 MB stale -wal. The
        # 1,1 header protects nothing -- a read-only open that is not immutable replays it.
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp, rollback_header=True)
            before = snapshot(d)
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "stale_side_file", r)
            self.assertFalse(r["sqlite_header_wal"])
            self.assertAlmostEqual(r["index_coverage"], 1.0)
            self.assertEqual(len(r["problems"]), 1, r["problems"])
            self.assertIn("read-only", r["problems"][0])
            self.assertIn("not immutable", r["problems"][0])
            self.assertEqual(snapshot(d), before, "the check itself must write nothing")

    def test_wal_header_alone_is_at_risk(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp, wal=True)
            self.assertNotIn("analysis.tdf-wal", os.listdir(d))
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "at_risk", r)
            self.assertTrue(r["sqlite_header_wal"])

    def test_non_empty_wal_alone_is_stale_side_file(self):
        # header 1,1, so only the -wal size can raise it
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            with open(os.path.join(d, "analysis.tdf-wal"), "wb") as fh:
                fh.write(b"\x37\x7f\x06\x82" + b"\0" * 28)
            r = bruker_tdf.tdf_integrity(d)
            self.assertFalse(r["sqlite_header_wal"])
            self.assertEqual(r["status"], "stale_side_file", r)
            self.assertEqual(r["wal_bytes"], 32)
            self.assertTrue(any("analysis.tdf-wal" in p for p in r["problems"]), r["problems"])

    def test_empty_wal_beside_a_rollback_file_is_not_a_problem(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            open(os.path.join(d, "analysis.tdf-wal"), "wb").close()
            self.assertEqual(bruker_tdf.tdf_integrity(d)["status"], "ok")

    def test_hot_rollback_journal_is_stale_side_file(self):
        # a read-write open would roll a non-empty -journal back into the file
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            with open(os.path.join(d, "analysis.tdf-journal"), "wb") as fh:
                fh.write(b"\xd9\xd5\x05\xf9\x20\xa1\x63\xd7" + b"\0" * 504)
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "stale_side_file", r)
            self.assertGreater(r["journal_bytes"], 0)
            self.assertTrue(any("analysis.tdf-journal" in p for p in r["problems"]),
                            r["problems"])

    def test_side_file_outranks_the_wal_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_stale_wal_d(tmp))
            self.assertEqual(r["status"], "stale_side_file")
            self.assertEqual(len(r["problems"]), 2, r["problems"])
            self.assertIn("analysis.tdf-wal", r["problems"][0])
            self.assertIn("WAL mode", r["problems"][1])

    def test_index_past_the_end_of_tdf_bin_is_bin_incomplete(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_intact_d(tmp, n_bin=N_FRAMES // 2))
            self.assertEqual(r["status"], "bin_incomplete", r)

    def test_last_block_running_past_the_end_of_tdf_bin_is_bin_incomplete(self):
        # the last indexed block STARTS inside tdf_bin but ends past it: a copy cut off
        # inside its final frame
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp)
            b = os.path.join(d, "analysis.tdf_bin")
            with open(b, "r+b") as fh:
                fh.truncate(N_FRAMES * BLOCK - BLOCK // 2)
            r = bruker_tdf.tdf_integrity(d)
            self.assertEqual(r["status"], "bin_incomplete", r)
            self.assertEqual(r["index_end_bytes"], N_FRAMES * BLOCK)
            self.assertTrue(any("ends at byte" in p for p in r["problems"]), r["problems"])

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


class TestIntegrityWarning(unittest.TestCase):
    """The remedy is what the user acts on, so its wording is pinned per state."""

    def test_every_status_but_ok_has_a_remedy(self):
        for status in bruker_tdf.STATUSES:
            if status != "ok":
                self.assertIn(status, bruker_tdf._REMEDY)

    def test_stale_side_file_says_do_not_search_and_how_to_clear_it(self):
        # a read-only engine that is not immutable would search the mid-acquisition index; a
        # read-write one would truncate the file. Backing up analysis.tdf prevents neither.
        for rollback_header in (True, False):
            with self.subTest(rollback_header=rollback_header), \
                    tempfile.TemporaryDirectory() as tmp:
                w = bruker_tdf.integrity_warning(bruker_tdf.tdf_integrity(
                    make_stale_wal_d(tmp, rollback_header=rollback_header)))
                self.assertIn("Do NOT search", w)
                self.assertNotIn("can be searched", w)
                self.assertIn("backup outside the .d", w)
                self.assertIn("-wal", w)
                self.assertIn("-shm", w)
                self.assertIn("re-check", w)

    def test_wal_header_alone_can_be_searched(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = bruker_tdf.integrity_warning(bruker_tdf.tdf_integrity(
                make_intact_d(tmp, wal=True)))
            self.assertIn("can be searched", w)
            self.assertNotIn("Do NOT", w)
            self.assertIn("immutable", w)


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
            self.assertFalse(out["needs_confirmation"], out)
            self.assertNotIn("WARNING", err)

    def test_truncated_run_warns_and_needs_confirmation(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = make_intact_d(tmp)
            bad = truncate_like_hive(make_stale_wal_d(tmp, name="bad.d"))
            out, err = self._run(good, bad)
            self.assertTrue(out["needs_confirmation"])
            by_file = {f["file"]: f for f in out["files"]}
            self.assertEqual(by_file[bad]["tdf_integrity"]["status"], "truncated")
            self.assertEqual(by_file[good]["warnings"], [])
            self.assertEqual(len(by_file[bad]["warnings"]), 1, by_file[bad]["warnings"])
            self.assertIn("truncated", by_file[bad]["warnings"][0])
            self.assertIn("Do NOT search", by_file[bad]["warnings"][0])
            # a person reading the terminal sees it too, with the file named
            self.assertIn(f"WARNING: {bad}: ", err)
            # acquisition detection itself is unchanged -- the warning is what gates
            self.assertEqual(out["overall"], "DIA")
            self.assertEqual(out["low_confidence_files"], [])

    def test_stale_wal_run_says_do_not_search_and_needs_confirmation(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp, rollback_header=True)
            out, err = self._run(d)
            self.assertTrue(out["needs_confirmation"])
            f = out["files"][0]
            self.assertEqual(f["tdf_integrity"]["status"], "stale_side_file")
            self.assertIn("Do NOT search", f["warnings"][0])
            self.assertNotIn("can be searched", f["warnings"][0])
            self.assertIn(f"WARNING: {d}: analysis.tdf stale_side_file", err)

    def test_wal_header_only_run_warns_but_is_searchable(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = self._cli(make_intact_d(tmp, wal=True))
            self.assertTrue(out["needs_confirmation"])
            f = out["files"][0]
            self.assertEqual(f["tdf_integrity"]["status"], "at_risk")
            self.assertIn("can be searched", f["warnings"][0])

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
    """The whole checkout, not this one scripts/ directory.

    The first version of this test walked `os.listdir(SCRIPTS)` only. Three raw
    `sqlite3.connect(f"file:{tdf}?mode=ro", uri=True)` calls sat one directory away the
    whole time, in skill/ucdavis-proteomics-dev/ -- a second plugin in the same
    marketplace.json, so it installs for everyone who adds the marketplace, points at the
    same real .d, and its SKILL.md saying "not for production" is an instruction to a
    model, not packaging. A per-directory allowlist is the thing that failed; the scan is
    the repo.
    """

    def test_no_python_in_the_repo_opens_a_tdf_without_immutable(self):
        self.assertEqual(python_tdf_offenders(ROOT), [],
                         "open analysis.tdf with bruker_tdf.connect_tdf() (mode=ro&immutable=1)")

    def test_no_r_in_the_repo_opens_a_tdf_without_immutable(self):
        # the R side has the same hazard and its own helper (R/helpers_instrument.R:
        # tdf_dbconnect_ro / sqlite_immutable_uri)
        self.assertEqual(r_tdf_offenders(ROOT), [],
                         "open analysis.tdf with tdf_dbconnect_ro() (mode=ro&immutable=1)")

    def test_the_scan_reaches_past_this_skill(self):
        """A scan that silently covered nothing would pass the two tests above for ever."""
        scanned = scanned_files(ROOT)
        self.assertIn(os.path.join(SCRIPTS, "detect_acquisition.py"), scanned)
        if ROOT != os.path.dirname(HERE):        # a full checkout, not a lone installed skill
            for rel in ("skill/ucdavis-proteomics-dev/scripts/detect_acquisition.py",
                        "skill/ucdavis-proteomics-dev/scripts/make_methods.py",
                        "cascadia/bruker_mobility_filter.py",
                        "R/helpers_instrument.R"):
                self.assertIn(os.path.join(ROOT, *rel.split("/")), scanned, rel)

    def test_the_scan_catches_a_planted_offender(self):
        """...and a scan that found files but waved every open through would too."""
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "scripts"))
            os.makedirs(os.path.join(tmp, "R"))
            plant = {
                "scripts/mode_ro_only.py":
                    'tdf = "analysis.tdf"\ncon = sqlite3.connect(f"file:{tdf}?mode=ro", uri=True)\n',
                "scripts/bare_path.py":
                    'tdf = "analysis.tdf"\ncon = sqlite3.connect(tdf)\n',
                "scripts/fine.py":
                    'tdf = "analysis.tdf"\ncon = connect_tdf(tdf)\n',
                "scripts/also_fine.py":
                    'tdf = "analysis.tdf"\ncon = sqlite3.connect(tdf_uri(tdf), uri=True)\n',
                "scripts/not_a_tdf.py":
                    'con = sqlite3.connect("results.sqlite")\n',
                # The exemption is a NAME and it is granted per CALL: the third line here
                # is a fixture being built and passes, the fourth is a plain read-write
                # open in the very same file and must still be caught.
                "scripts/fixture_writer.py":
                    'from synthetic_tdf import synthetic_tdf_write_uri\n'
                    'tdf = "analysis.tdf"\n'
                    'con = sqlite3.connect(synthetic_tdf_write_uri(tdf), uri=True)\n'
                    'sneaked = sqlite3.connect(tdf)\n',
                "R/mode_ro_only.R":
                    'tdf <- "analysis.tdf"\ncon <- DBI::dbConnect(RSQLite::SQLite(), tdf)\n',
                "R/fine.R":
                    'tdf <- "analysis.tdf"\ncon <- tdf_dbconnect_ro(tdf)\n',
            }
            for rel, text in plant.items():
                with open(os.path.join(tmp, *rel.split("/")), "w", encoding="utf-8") as fh:
                    fh.write(text)
            py = python_tdf_offenders(tmp)
            self.assertEqual(sorted(o.split(":")[0] for o in py),
                             ["scripts/bare_path.py", "scripts/fixture_writer.py",
                              "scripts/mode_ro_only.py"], py)
            self.assertEqual([o for o in py if o.startswith("scripts/fixture_writer.py")],
                             ["scripts/fixture_writer.py:4: sqlite3.connect(tdf)"], py)
            r = r_tdf_offenders(tmp)
            self.assertEqual([o.split(":")[0] for o in r], ["R/mode_ro_only.R"], r)

    def test_the_named_fixture_writer_really_opens_a_writable_database(self):
        """The exemption above is a NAME, so the name has to do what it says.

        If synthetic_tdf_write_uri() ever stopped granting write access the fixtures would
        break loudly -- but if it quietly started handing out a READING uri, the guard
        would go on waving every call through while the calls did something else. Pin both
        halves: it writes, and it is not immutable, so nothing can mistake it for the
        reader."""
        uri = synthetic_tdf_write_uri("/some/run.d/analysis.tdf")
        self.assertIn("mode=rwc", uri)
        self.assertNotIn("immutable", uri)
        with tempfile.TemporaryDirectory() as tmp:
            tdf = os.path.join(tmp, "analysis.tdf")
            con = sqlite3.connect(synthetic_tdf_write_uri(tdf), uri=True)
            con.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY)")
            con.execute("INSERT INTO Frames VALUES (1)")
            con.commit()
            con.close()
            self.assertGreater(os.path.getsize(tdf), 0, "the fixture writer wrote nothing")

    def test_the_skill_scripts_here_use_the_shared_helper(self):
        for name in ("detect_acquisition.py", "make_methods.py"):
            with open(os.path.join(SCRIPTS, name), encoding="utf-8") as fh:
                src = fh.read()
            self.assertIn("from bruker_tdf import", src, name)
            self.assertNotRegex(src, r"sqlite3\.connect\s*\(", name)

    def test_the_dev_fork_uses_the_same_helper_file(self):
        """The fork ships its own scripts/, so it needs its own copy -- kept byte-identical.

        A fork that edits this one module on its own is how the safety fix gets lost again:
        the fork is what a user who typed "dev skill" points at real .d with."""
        fork = os.path.join(ROOT, "skill", "ucdavis-proteomics-dev", "scripts")
        if not os.path.isdir(fork):
            self.skipTest("no dev fork in this tree (skill installed on its own)")
        with open(os.path.join(SCRIPTS, "bruker_tdf.py"), "rb") as fh:
            stable = fh.read()
        with open(os.path.join(fork, "bruker_tdf.py"), "rb") as fh:
            forked = fh.read()
        self.assertEqual(forked, stable,
                         "copy scripts/bruker_tdf.py over skill/ucdavis-proteomics-dev/"
                         "scripts/bruker_tdf.py -- they must not drift")
        for name in ("detect_acquisition.py", "make_methods.py"):
            with open(os.path.join(fork, name), encoding="utf-8") as fh:
                src = fh.read()
            self.assertIn("from bruker_tdf import", src, name)
            self.assertNotRegex(src, r"sqlite3\.connect\s*\(", name)


class TestTrailingSlackOnASmallRun(unittest.TestCase):
    """A relative coverage margin alone is wrong at the small end.

    Four intact HIVE runs' frame indexes ended 0.8-3.5 KB short of their tdf_bin. 0.1% of a
    3.5 MB tdf_bin is 3.5 KB, so every blank and QC injection under about that size read as
    `truncated` -- a warning on a healthy file, in the one place the user is being asked to
    believe a warning."""

    SLACK_BLOCKS = 14                        # 14 * 256 B = 3,584 B, the largest slack measured

    def test_a_few_kb_of_slack_on_a_small_tdf_bin_is_not_truncation(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_intact_d(tmp, n_indexed=N_FRAMES, n_bin=N_FRAMES + self.SLACK_BLOCKS)
            r = bruker_tdf.tdf_integrity(d)
            self.assertLess(r["index_coverage"], bruker_tdf.COVERAGE_MIN,
                            "fixture must be below the RELATIVE margin, or this proves nothing")
            self.assertEqual(r["status"], "ok", r)

    def test_the_same_shortfall_on_a_big_tdf_bin_is_still_truncation(self):
        # the absolute floor must not swallow real damage: the relative margin decides once
        # tdf_bin is big enough for 0.1% of it to exceed SLACK_BYTES
        size = 400 * 1024 * 1024
        self.assertGreater((1 - bruker_tdf.COVERAGE_MIN) * size, bruker_tdf.SLACK_BYTES)
        self.assertEqual(bruker_tdf.allowed_slack(size), (1 - bruker_tdf.COVERAGE_MIN) * size)

    def test_a_real_short_index_on_a_small_run_is_still_caught(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_intact_d(tmp, n_indexed=N_STALE))
            self.assertEqual(r["status"], "truncated", r)


class TestEveryStatusGetsItsRemedy(unittest.TestCase):
    """The worst status decides the headline; it must not be the only thing acted on.

    A .d with a hot -journal AND an unreadable index is two separate problems with two
    separate fixes, and one of them names a file the user can move."""

    def _journal_and_unreadable_index(self, tmp):
        d = make_intact_d(tmp)
        with open(os.path.join(d, "analysis.tdf"), "wb") as fh:
            fh.write(b"not a database" * 10)
        with open(os.path.join(d, "analysis.tdf-journal"), "wb") as fh:
            fh.write(b"\xd9\xd5\x05\xf9\x20\xa1\x63\xd7" + b"\0" * 504)
        return d

    def test_a_movable_side_file_outranks_could_not_read_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(self._journal_and_unreadable_index(tmp))
            self.assertEqual(r["status"], "stale_side_file", r)
            self.assertEqual(r["statuses"], ["stale_side_file", "unverified"], r)

    def test_both_remedies_reach_the_user(self):
        with tempfile.TemporaryDirectory() as tmp:
            w = bruker_tdf.integrity_warning(
                bruker_tdf.tdf_integrity(self._journal_and_unreadable_index(tmp)))
            self.assertIn("backup outside the .d", w)          # the stale_side_file remedy
            self.assertIn(bruker_tdf._REMEDY["unverified"], w)  # and the unverified one

    def test_searchable_as_it_is_is_never_said_alongside_a_worse_finding(self):
        # at_risk's remedy is a verdict on the whole file; a stale -wal makes it false
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_stale_wal_d(tmp))
            self.assertIn("at_risk", r["statuses"])
            self.assertIn("can be searched", bruker_tdf._REMEDY["at_risk"])
            self.assertNotIn("can be searched", bruker_tdf.integrity_warning(r))

    def test_an_ok_file_gets_an_empty_warning_not_a_crash(self):
        # guarded at the one caller today, but this is a shared module: a caller that maps
        # it over every file must not blow up on the healthy ones
        with tempfile.TemporaryDirectory() as tmp:
            r = bruker_tdf.tdf_integrity(make_intact_d(tmp))
            self.assertEqual(r["status"], "ok")
            self.assertEqual(bruker_tdf.integrity_warning(r), "")


class TestRemediesNameTheStateTheyApplyTo(unittest.TestCase):
    def test_truncated_names_a_run_that_is_still_being_acquired(self):
        """A .d being written right now looks exactly like a truncated one. Sending that user
        to hunt for another copy of a file that is simply not finished is the wrong errand."""
        w = bruker_tdf._REMEDY["truncated"]
        self.assertIn("ACQUIRED", w)
        self.assertIn("wait", w)

    def test_stale_side_file_does_not_tell_a_live_acquisition_to_move_its_own_files(self):
        """The remedy moves files OUT of a raw .d. Beside a running acquisition those files
        are in use, not stale, and moving them is how you lose the run."""
        w = bruker_tdf._REMEDY["stale_side_file"]
        self.assertIn("only once acquisition has finished", w)
        self.assertIn("nothing holds the .d open", w)

    def test_at_risk_states_both_caveats_next_to_searchable_as_it_is(self):
        """True for the truncation mechanism, and only for that: a WAL-mode header still
        means a read-write open writes into the raw .d, and that WAL needs locks a network
        mount does not have."""
        w = bruker_tdf._REMEDY["at_risk"]
        self.assertIn("can be searched", w)
        self.assertIn("analysis.tdf-shm", w)
        for mount in ("NFS", "SMB", "Quobyte"):
            self.assertIn(mount, w)


class TestRunSearchChecksTheTdfItself(unittest.TestCase):
    """The integrity check lived only in detect_acquisition.py -- step 2 of the skill.

    Every "just search these" entry point in SKILL.md's own description ("re-run this
    search", "re-search with different parameters", "the paths are already known") arrives
    at run_search.py with the files in hand and no step 2 behind it. A truncated .d was
    searched there with nothing said."""

    def _args(self, tmp, *files, extra=()):
        for name, text in (("tools.json", "{}"), ("bundle.json", "{}"),
                           ("p.cfg", "--qvalue 0.01\n"), ("db.fasta", ">sp|P1|X\nPEPTIDER\n")):
            with open(os.path.join(tmp, name), "w") as fh:
                fh.write(text)
        return [sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                "--tools", os.path.join(tmp, "tools.json"),
                "--bundle", os.path.join(tmp, "bundle.json"),
                "--params", os.path.join(tmp, "p.cfg"),
                "--fasta", os.path.join(tmp, "db.fasta"),
                "--out", os.path.join(tmp, "out"),
                *extra, "--files", *files]

    def _run(self, tmp, *files, extra=()):
        return subprocess.run(self._args(tmp, *files, extra=extra),
                              capture_output=True, text=True)

    def test_tdf_problems_finds_the_damaged_d_and_leaves_the_rest_alone(self):
        with tempfile.TemporaryDirectory() as tmp:
            good = make_intact_d(tmp)
            bad = truncate_like_hive(make_stale_wal_d(tmp, name="bad.d"))
            mzml = os.path.join(tmp, "run.mzML")
            open(mzml, "w").close()
            found = run_search.tdf_problems([good, bad, mzml])
            self.assertEqual([p for p, _ in found], [bad])
            self.assertEqual(found[0][1]["status"], "truncated")

    def test_an_intact_cohort_has_no_problems(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(run_search.tdf_problems([make_intact_d(tmp)]), [])

    def test_a_trailing_slash_does_not_hide_a_d(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = truncate_like_hive(make_stale_wal_d(tmp, name="bad.d"))
            self.assertEqual([p for p, _ in run_search.tdf_problems([bad + "/"])], [bad])

    def test_the_search_is_refused_before_anything_is_provisioned(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = truncate_like_hive(make_stale_wal_d(tmp, name="bad.d"))
            before = snapshot(bad)
            res = self._run(tmp, bad)
            self.assertNotEqual(res.returncode, 0, res.stdout)
            self.assertIn("REFUSING to search", res.stderr)
            self.assertIn(bad, res.stderr)
            self.assertIn("analysis.tdf truncated", res.stderr)
            self.assertIn("Do NOT search", res.stderr)
            self.assertIn("--allow-damaged-tdf", res.stderr)
            # it stopped before it got as far as looking for an engine
            self.assertNotIn("has no command for engine", res.stderr)
            self.assertEqual(snapshot(bad), before, "run_search wrote beside the tdf")
            self.assertFalse(os.path.exists(os.path.join(tmp, "out")))

    def test_a_stale_side_file_is_refused_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_stale_wal_d(tmp, rollback_header=True)
            res = self._run(tmp, d)
            self.assertNotEqual(res.returncode, 0, res.stdout)
            self.assertIn("analysis.tdf stale_side_file", res.stderr)

    def test_an_intact_run_is_not_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            res = self._run(tmp, make_intact_d(tmp))
            self.assertNotIn("REFUSING to search", res.stderr)
            # it got past the gate and on to the next thing run_search does
            self.assertIn("has no command for engine", res.stderr)

    def test_the_override_says_what_it_is_letting_through(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = truncate_like_hive(make_stale_wal_d(tmp, name="bad.d"))
            res = self._run(tmp, bad, extra=("--allow-damaged-tdf",))
            self.assertNotIn("REFUSING to search", res.stderr)
            self.assertIn("--allow-damaged-tdf", res.stdout)
            self.assertIn("analysis.tdf truncated", res.stdout)
            self.assertIn("has no command for engine", res.stderr)

    def test_a_cohort_of_non_bruker_files_is_never_gated(self):
        with tempfile.TemporaryDirectory() as tmp:
            mzml = os.path.join(tmp, "run.mzML")
            open(mzml, "w").close()
            res = self._run(tmp, mzml)
            self.assertNotIn("REFUSING to search", res.stderr)



if __name__ == "__main__":
    unittest.main(verbosity=2)
