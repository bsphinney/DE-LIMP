#!/usr/bin/env python3
"""The one deliberate WRITE to a Bruker analysis.tdf in this repo: building a fixture.

test_tdf_readonly_open.py walks the whole checkout and fails on every `sqlite3.connect`
of an analysis.tdf that is not immutable. That is the Hive incident (2026-09-16): a stale
mid-acquisition `-wal` beside a finished run is REPLAYED by a read-write open, which
checkpoints the mid-run pages into the file and truncates it -- 342 runs, one of them
indexing 0.7% of its own 2.4 GB of spectra. Nothing that READS a run may write to it.

A test that MAKES a synthetic .d has to write: there is no file yet, or it is a throwaway
copy being damaged on purpose to prove the fixture really carries the hazard. That is a
different act from reading a run, and it is named here so that it says so at the call site:

    con = sqlite3.connect(synthetic_tdf_write_uri(tdf), uri=True)

The guard exempts THAT NAME, not the files that use it. A per-file allowlist is the thing
that already failed once (see TestNoOtherTdfOpens) and it gains an entry every time someone
adds a fixture, until it is no longer a guard. A name cannot be written by accident, it is
granted per CALL rather than per file -- a bare connect on the path, in this very module,
would still be an offender -- and `grep -rn synthetic_tdf_write_uri` lists every deliberate
writer in the repo. Never point it at a real run.

The guard reads one line at a time, so the name has to sit on the same line as the
connect it makes safe; a call split across lines is reported, which is the right way
round for a guard to be wrong.
"""
import os
import pathlib


def synthetic_tdf_write_uri(path):
    """A sqlite URI that opens `path` READ-WRITE, creating it if need be -- fixtures only.

    Percent-encoded exactly as bruker_tdf.tdf_uri() does it, so a `?`, `#` or `%` in a
    temporary directory name cannot end the path and send the open at another file."""
    return pathlib.Path(os.path.abspath(path)).as_uri() + "?mode=rwc"
