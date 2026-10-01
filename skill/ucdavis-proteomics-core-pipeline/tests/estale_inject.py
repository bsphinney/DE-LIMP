"""Inject a stale NFS file handle (ESTALE) into the probe's own process, as HIVE's Flinders NFS did.

fran-5b, 2026-09-29/30: 62 of 261 step-1b jobs died in probe_window.run_probe with
`OSError: [Errno 116] Stale file handle` at `pending += src.read()` -- the tail of the probe's
DIA-NN log. The injection is LOCATION-INDEPENDENT: any read of a file named `probe.log`, opened
for binary reading, wherever it is and whoever reads it -- so the same test crashes the probe
of skill <= 979dfd5 (which tails <workdir>/probe.log) and must not crash the fixed one (which
tails a node-local copy). It rides in on a sitecustomize.py on PYTHONPATH, so it reaches the
probe inside a generated job script too.

    env.update(estale_env(tmpdir, "transient"))

  transient  the 2nd read of each probe.log, in every process, fails ONCE (a handle that went
             stale mid-tail and would have worked when reopened); give DIA-NN's first lines a
             moment (the fakes' FAKE_FIRST) so there IS a 2nd read
  once       every read of a probe.log fails -- in the FIRST process that reads one only
             (attempt 1 of a job; the retry is clean)
  always     every read of a probe.log fails, in every process
"""
import os

SITECUSTOMIZE = r'''
import builtins, errno, os
_real_open = builtins.open
_reads = {}
_cursed = []


class _Stale:
    """The real handle, whose read() goes stale as FAKE_ESTALE says."""

    def __init__(self, fh, path):
        self._fh, self._path = fh, path

    def read(self, *a):
        n = _reads[self._path] = _reads.get(self._path, 0) + 1
        how, at = os.environ.get("FAKE_ESTALE"), int(os.environ.get("FAKE_ESTALE_AT", "2"))
        if how == "transient" and n == at or how in ("once", "always") and n >= at and _cursed[0]:
            raise OSError(errno.ESTALE, os.strerror(errno.ESTALE))
        return self._fh.read(*a)

    def __getattr__(self, name):
        return getattr(self._fh, name)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self._fh.close()

    def __iter__(self):
        return iter(self._fh)


def _open(file, mode="r", *a, **k):
    fh = _real_open(file, mode, *a, **k)
    how = os.environ.get("FAKE_ESTALE")
    if how and "r" in mode and "b" in mode and "+" not in mode \
            and os.path.basename(str(file)) == "probe.log":
        if not _cursed:
            marker = os.environ.get("FAKE_ESTALE_MARKER", "")
            _cursed.append(how == "always" or not (marker and os.path.exists(marker)))
            if marker and _cursed[0]:
                _real_open(marker, "w").close()
        return _Stale(fh, os.path.realpath(str(file)))
    return fh


builtins.open = _open
'''


def estale_env(tmpdir, mode, at=None):
    """Environment additions that inject ESTALE (`mode`: transient | once | always) into every
    Python process started with them. Writes <tmpdir>/estale_site/sitecustomize.py."""
    site = os.path.join(tmpdir, "estale_site")
    os.makedirs(site, exist_ok=True)
    with open(os.path.join(site, "sitecustomize.py"), "w") as fh:
        fh.write(SITECUSTOMIZE)
    at = at if at is not None else (2 if mode == "transient" else 1)
    return {"PYTHONPATH": site, "FAKE_ESTALE": mode, "FAKE_ESTALE_AT": str(at),
            "FAKE_ESTALE_MARKER": os.path.join(tmpdir, "estale_marker")}
