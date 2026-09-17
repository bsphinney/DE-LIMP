"""
bruker_tdf.py -- the ONE way the skill opens a Bruker analysis.tdf, and a cheap check
that the tdf still indexes the whole run.

Why a helper and not a one-line `sqlite3.connect`: on UC Davis HIVE (measured
2026-09-16) 342 of 39,374 Bruker .d had an analysis.tdf whose frame index covered only
part of a complete analysis.tdf_bin -- one run's index ended at 0.7% of its tdf_bin, and
a DIA-NN search of it read 121 cycles without an error anywhere. Every truncated file was
still in WAL mode (SQLite header bytes 18-19 = 2,2; intact files are 1,1). The copy from
the instrument can leave a stale, mid-acquisition analysis.tdf-wal beside the finished
file, and a READ-WRITE sqlite open replays it: SQLite checkpoints the mid-run pages into
the finished file, cuts it down to the size the -wal recorded, and deletes the -wal. The
file keeps its WAL header and otherwise looks normal.

    sqlite3.connect(tdf)                       read-write: truncates the file
    sqlite3.connect("file:<tdf>?mode=ro")      no truncation, but READS the stale -wal (the
                                               reader sees the mid-run database) and drops
                                               an analysis.tdf-shm beside the tdf
    connect_tdf(tdf)  (mode=ro&immutable=1)    reads the file as it is on disk, takes no
                                               locks, writes nothing

immutable=1 is right for a tdf because nothing writes one after acquisition, and taking no
locks is what a network mount wants anyway.

tdf_integrity() is the check to run before anything searches a .d. It needs no Bruker
library: the SQLite header, the size of any -wal/-journal beside the tdf, and where the
last frame's block ends in tdf_bin (Frames.TimsId is that block's byte offset, and a block
starts with its own uint32 size).

Stdlib only; imported by detect_acquisition.py and make_methods.py.
"""
import os
import pathlib
import sqlite3
import struct

# An index that ends before this fraction of tdf_bin is reported as truncated. The last
# block normally ends exactly at the end of the file; the margin only keeps trailing
# padding from ever reading as damage. The truncated files on HIVE were far below it.
COVERAGE_MIN = 0.999

SQLITE_MAGIC = b"SQLite format 3\x00"

# worst first; tdf_integrity() reports the worst status it found
STATUSES = ("truncated", "bin_incomplete", "unverified", "at_risk", "ok")


def tdf_uri(path):
    """The sqlite URI that opens `path` read-only and immutable.

    pathlib percent-encodes the path: a `?`, `#` or `%` in a folder name spliced into
    `file:` raw would end the path (or start an escape), and the open would silently target
    a different name."""
    return pathlib.Path(os.path.abspath(path)).as_uri() + "?mode=ro&immutable=1"


def connect_tdf(path):
    """Open a Bruker analysis.tdf so nothing is written to it or beside it and no stale -wal
    is replayed. Raises sqlite3.OperationalError for a missing file (mode=ro never creates
    one) and for any attempt to write."""
    return sqlite3.connect(tdf_uri(path), uri=True)


def _size(path):
    try:
        return os.path.getsize(path)
    except OSError:
        return 0


def tdf_integrity(d_path):
    """Does this .d's analysis.tdf still index its whole analysis.tdf_bin, and is it exposed
    to a read-write open? Reads only; writes nothing. Returns

      status             ok | at_risk | truncated | bin_incomplete | unverified
      problems           one sentence per finding ([] when ok)
      sqlite_header_wal  header bytes 18-19 say WAL mode (None when the header is unreadable)
      wal_bytes          size of analysis.tdf-wal (0 when absent)
      journal_bytes      size of analysis.tdf-journal (0 when absent)
      tdf_bin_bytes      size of analysis.tdf_bin (None when absent)
      index_end_bytes    the last frame's TimsId + its block size (None when not read)
      index_coverage     index_end_bytes / tdf_bin_bytes (None when not read)

    `status` is the worst finding: truncated (the index ends before COVERAGE_MIN of tdf_bin,
    so a search reads only part of the run) > bin_incomplete (tdf_bin missing, or shorter
    than the index) > unverified (the index could not be read) > at_risk (the index is
    complete, but the header is in WAL mode or a non-empty -wal/-journal sits beside it, so
    one read-write open can truncate it) > ok."""
    tdf = os.path.join(d_path, "analysis.tdf")
    tdf_bin = os.path.join(d_path, "analysis.tdf_bin")
    r = {"status": "ok", "problems": [], "sqlite_header_wal": None,
         "wal_bytes": _size(tdf + "-wal"), "journal_bytes": _size(tdf + "-journal"),
         "tdf_bin_bytes": os.path.getsize(tdf_bin) if os.path.isfile(tdf_bin) else None,
         "index_end_bytes": None, "index_coverage": None}
    findings = []                           # (status, sentence)

    def problem(status, text):
        findings.append((status, text))

    # 1. the header, read as bytes -- never through sqlite
    try:
        with open(tdf, "rb") as fh:
            hdr = fh.read(100)
    except OSError as e:
        hdr = None
        problem("unverified", f"analysis.tdf cannot be read ({e})")
    if hdr and not hdr.startswith(SQLITE_MAGIC):
        problem("unverified", "analysis.tdf is not an SQLite database")
    elif hdr and len(hdr) >= 20:
        r["sqlite_header_wal"] = 2 in (hdr[18], hdr[19])
        if r["sqlite_header_wal"]:
            problem("at_risk",
                    f"analysis.tdf is still in WAL mode (SQLite header bytes 18-19 = "
                    f"{hdr[18]},{hdr[19]}; a finished tdf is 1,1) -- every truncated tdf found "
                    f"on HIVE was, and a read-write sqlite open is what truncated them")

    # 2. what a read-write open would replay into it
    if r["wal_bytes"]:
        problem("at_risk",
                f"a non-empty analysis.tdf-wal ({r['wal_bytes']:,} bytes) sits beside it, which "
                f"a read-write sqlite open would checkpoint into the finished file")
    if r["journal_bytes"]:
        problem("at_risk",
                f"a non-empty analysis.tdf-journal ({r['journal_bytes']:,} bytes) sits beside "
                f"it, which a read-write sqlite open would roll back into the finished file")

    # 3. does the index reach the end of tdf_bin?
    if not any(s == "unverified" for s, _ in findings):
        last = None
        try:
            con = connect_tdf(tdf)
            try:
                last = con.execute(
                    "SELECT TimsId FROM Frames ORDER BY TimsId DESC LIMIT 1").fetchone()
            finally:
                con.close()
            if last is None:
                problem("unverified", "Frames has no rows, so the index cannot be checked")
        except sqlite3.Error as e:
            problem("unverified", f"the frame index could not be read ({e})")
        size = r["tdf_bin_bytes"]
        offset = last[0] if last else None
        block = None
        if last is None:
            pass                            # already reported as unverified
        elif size is None:
            problem("bin_incomplete", "analysis.tdf_bin is missing, so the run has no spectra")
        elif not isinstance(offset, int) or offset < 0:
            problem("unverified", f"the last frame's TimsId is not a byte offset ({offset!r})")
        elif offset + 4 > size:
            problem("bin_incomplete",
                    f"the index points at byte {offset:,} but analysis.tdf_bin is only "
                    f"{size:,} bytes, so the binary is incomplete")
        else:
            try:
                with open(tdf_bin, "rb") as fh:
                    fh.seek(offset)
                    block = struct.unpack("<I", fh.read(4))[0]
            except (OSError, struct.error) as e:
                problem("unverified", f"analysis.tdf_bin cannot be read at byte {offset:,} ({e})")
        if block is not None:
            end = offset + block
            r["index_end_bytes"] = end
            r["index_coverage"] = end / size
            if end > size:
                problem("bin_incomplete",
                        f"the last frame's block ends at byte {end:,} but analysis.tdf_bin is "
                        f"only {size:,} bytes, so the binary is incomplete")
            elif end < COVERAGE_MIN * size:
                problem("truncated",
                        f"analysis.tdf is truncated: its frame index ends at byte {end:,} of the "
                        f"{size:,}-byte analysis.tdf_bin ({100 * end / size:.1f}%), so a search "
                        f"reads only that part of the run, silently")

    # worst first, so the finding that decides what to do is the one read first
    findings.sort(key=lambda f: STATUSES.index(f[0]))
    r["problems"] = [text for _, text in findings]
    r["status"] = findings[0][0] if findings else "ok"
    return r


_REMEDY = {
    "truncated": "Do NOT search this run as it is: find an intact analysis.tdf (another copy "
                 "of the run with the same analysis.tdf_bin, or the instrument PC), then re-check",
    "bin_incomplete": "Do NOT search this run: re-copy the .d from its source, then re-check",
    "unverified": "Confirm this .d is readable before searching it",
    "at_risk": "Its index is complete, so it can be searched, but never open this tdf "
               "read-write (e.g. the sqlite3 CLI on the bare path), and back up "
               "analysis.tdf (it is small) first -- whether each search engine opens it "
               "read-only is unverified",
}


def integrity_warning(r):
    """One line for a file's `warnings` list: what was found and what to do."""
    return f"analysis.tdf {r['status']} -- {'; '.join(r['problems'])}. {_REMEDY[r['status']]}."
