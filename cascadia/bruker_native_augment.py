"""Cascadia's own augmentation, fed from a Bruker diaPASEF .d.

Cascadia was trained on spectra built by ``cascadia.augment`` from Orbitrap mzML. Each augmented
spectrum there carries the centre MS2 scan's fragments (top 150, sqrt-scaled) PLUS the nearest
MS1 survey scan's peaks around the precursor (top 150, 4th-root-scaled), and every peak's
retention-time offset from the centre scan; the model encodes both the RT column and the
MS-level column (model.py rt_encoder / level_encoder).

``cascadia.bruker_augment`` (the .d route patched into cascadia.py) writes co-added fragments
only, every RT offset 0.0, no MS1 peaks, and the CLI then hard-codes a 25 m/z window and 1.0 s
cycle. This module instead hands ``cascadia.augment``'s UNMODIFIED functions a stream of
mzML-shaped records built from the .d, so windowing, neighbour selection, normalisation, MS1
filtering, charge enumeration and the returned window size / cycle time are Cascadia's own:

* MS2: one record per timsrust DIA spectrum (one diaPASEF window, collapsed over its mobility
  range), isolation window from ``isolation_mz`` / ``isolation_width``.
* MS1: one record per MS1 frame, summed over all mobility scans by TOF index, smoothed over +/-1
  TOF bin, centroided at local maxima (intensity-weighted TOF), converted to m/z with the same
  boundary formula timsrust uses for its MS2 m/z (MzAcqRangeLower/Upper, DigitizerNumSamples).
* RT: timsrust/Bruker ``Frames.Time`` is SECONDS; mzML 'scan start time' is minutes, so records
  carry rt/60 and ``augment`` multiplies back by 60.

m/z calibration: the boundary formula matches the instrument's stored calibration (Bruker
libtimsdata on analysis.tdf's MzCalibration) to within +/-3 ppm. What it cannot see is drift
since that calibration: usually a few ppm, but -25 ppm on the mouse benchmark run. ``calib``
applies a per-run linear ppm correction fitted by prepare_training_native.py (RUN.calib.json):
true = obs * (1 - (a + b*obs)*1e-6), separately for MS2 and MS1.

Peaks are pre-truncated to the ``top_n`` most intense (augment takes the same top_n again), which
keeps memory bounded on long runs. The raw .d is only read: analysis.tdf is opened
``mode=ro&immutable=1``.
"""
from __future__ import annotations

import contextlib
import logging
import os
import sqlite3
import time

import numpy as np

logger = logging.getLogger("bruker_native_augment")


def _tdf_meta(d_folder):
    tdf = os.path.join(str(d_folder), "analysis.tdf")
    conn = sqlite3.connect(f"file:{tdf}?mode=ro&immutable=1", uri=True)
    try:
        meta = dict(conn.execute("SELECT Key, Value FROM GlobalMetadata").fetchall())
    finally:
        conn.close()
    return float(meta["MzAcqRangeLower"]), float(meta["MzAcqRangeUpper"]), int(meta["DigitizerNumSamples"])


def tof_to_mz(tof, mz_lower, mz_upper, dns):
    """timsrust Tof2MzConverter::from_boundaries: sqrt(mz) is linear in the TOF index."""
    a = np.sqrt(mz_lower)
    b = (np.sqrt(mz_upper) - a) / dns
    return (a + b * np.asarray(tof, dtype=np.float64)) ** 2


def centroid_ms1_frame(tof_indices, intensities, conv, top_n=None):
    """Sum one MS1 frame over mobility by TOF index, smooth +/-1 bin, keep local maxima.

    Returns (mz, intensity) sorted by m/z; the top_n most intense if top_n is given."""
    tof = np.asarray(tof_indices, dtype=np.int64)
    if tof.size == 0:
        return np.empty(0), np.empty(0)
    w = np.asarray(intensities, dtype=np.float64)
    lo = int(tof.min()) - 2
    c = np.bincount(tof - lo, weights=w)
    c = np.concatenate([c, [0.0, 0.0]])
    left, right = np.roll(c, 1), np.roll(c, -1)
    left[0] = right[-1] = 0.0
    s = left + c + right                                   # +/-1 bin smoothing
    sl, sr = np.roll(s, 1), np.roll(s, -1)
    sl[0] = sr[-1] = 0.0
    idx = np.nonzero((s > 0) & (s >= sl) & (s > sr) & (c > 0))[0]
    if idx.size == 0:
        return np.empty(0), np.empty(0)
    # intensity-weighted TOF over the 3 bins, for a better centroid m/z
    num = (idx - 1) * left[idx] + idx * c[idx] + (idx + 1) * right[idx]
    tof_c = num / s[idx] + lo
    inten = s[idx]
    if top_n is not None and inten.size > top_n:
        keep = np.argsort(inten)[-top_n:]
        tof_c, inten = tof_c[keep], inten[keep]
    mz = conv(tof_c)
    o = np.argsort(mz)
    return mz[o], inten[o]


def _top(mz, inten, top_n):
    if top_n is not None and inten.size > top_n:
        keep = np.argsort(inten)[-top_n:]
        mz, inten = mz[keep], inten[keep]
    o = np.argsort(mz)
    return mz[o], inten[o]


def _correct(mz, coef):
    if not coef:
        return mz
    return mz * (1.0 - (float(coef["a"]) + float(coef["b"]) * mz) * 1e-6)


@contextlib.contextmanager
def _without_appledouble(d_folder):
    """Yield a .d path timsrust can open. A Mac copy leaves AppleDouble files such as
    ``._analysis.tdf`` beside the real one; timsrust picks the first ``*.tdf`` it finds, opens the
    4 KB metadata stub and fails with "file is not a database" (35 of 378 FRAN runs, 2026-11-16).
    When such files exist, read through a temporary folder of symlinks that leaves them out."""
    import shutil
    import tempfile
    names = os.listdir(d_folder)
    if not any(n.startswith("._") for n in names):
        yield d_folder
        return
    tmp = tempfile.mkdtemp(prefix="dclean_")
    clean = os.path.join(tmp, os.path.basename(d_folder.rstrip("/")))
    os.mkdir(clean)
    for n in names:
        if not n.startswith("._"):
            os.symlink(os.path.join(os.path.abspath(d_folder), n), os.path.join(clean, n))
    try:
        yield clean
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def read_bruker_records(d_folder, top_n=150, calib=None):
    """mzML-shaped records (the subset of keys cascadia.augment reads), in RT order.
    calib: {"ms2": {"a","b"}, "ms1": {"a","b"}} ppm correction, or None."""
    c2 = (calib or {}).get("ms2")
    c1 = (calib or {}).get("ms1")
    d_folder = str(d_folder)
    with _without_appledouble(d_folder) as d_read:
        return _read_bruker_records(d_read, top_n, c1, c2)


def _read_bruker_records(d_folder, top_n, c1, c2):
    import timsrust_pyo3 as timsrust

    mz_lower, mz_upper, dns = _tdf_meta(d_folder)
    conv = lambda t: tof_to_mz(t, mz_lower, mz_upper, dns)  # noqa: E731
    t0 = time.time()
    records = []
    n_empty2 = 0
    for sp in timsrust.read_all_spectra(d_folder):
        mz = np.asarray(sp.mz_values, dtype=np.float64)
        it = np.asarray(sp.intensities, dtype=np.float64)
        if mz.size == 0:
            n_empty2 += 1
            continue
        mz, it = _top(_correct(mz, c2), it, top_n)
        half = float(sp.isolation_width) / 2.0
        records.append((float(sp.precursor.rt), 2, {
            "ms level": 2,
            "scanList": {"scan": [{"scan start time": float(sp.precursor.rt) / 60.0}]},
            "precursorList": {"precursor": [{"isolationWindow": {
                "isolation window target m/z": float(sp.isolation_mz),
                "isolation window lower offset": half,
                "isolation window upper offset": half}}]},
            "m/z array": mz, "intensity array": it}))
    n2 = len(records)
    n_empty1 = 0
    for fr in timsrust.FrameReader(d_folder).read_ms1_frames():
        mz, it = centroid_ms1_frame(fr.tof_indices, fr.intensities, conv, top_n)
        mz = _correct(mz, c1)
        if mz.size == 0:
            n_empty1 += 1
            continue
        records.append((float(fr.rt), 1, {
            "ms level": 1,
            "scanList": {"scan": [{"scan start time": float(fr.rt) / 60.0}]},
            "m/z array": mz, "intensity array": it}))
    records.sort(key=lambda r: (r[0], r[1]))
    logger.info("%s: %d MS2 window spectra (%d empty), %d MS1 frames (%d empty) in %.0fs",
                os.path.basename(d_folder.rstrip("/")), n2, n_empty2, len(records) - n2, n_empty1,
                time.time() - t0)
    if n2 == 0 or len(records) == n2:
        raise RuntimeError(f"{d_folder}: need both MS2 window spectra and MS1 frames (got {n2} MS2, "
                           f"{len(records) - n2} MS1); is this a diaPASEF acquisition?")
    return [r[2] for r in records]


class _RecordReader:
    """Stands in for pyteomics.mzml inside cascadia.augment: read(path, ...) yields the records."""

    def __init__(self, records):
        self._records = records

    def read(self, *_a, **_k):
        return contextlib.nullcontext(iter(self._records))


def augment_spectra_bruker_native(d_folder, temp_path, top_n=150, scan_width=1, max_pep_length=30,
                                  max_charge=3, calib=None):
    """Drop-in for cascadia.augment.augment_spectra on a .d: returns (asf_path, window_size, cycle_time)."""
    from cascadia import augment as aug

    records = read_bruker_records(d_folder, top_n=top_n, calib=calib)
    out = os.path.join(temp_path, "temp.asf")
    if os.path.exists(out):
        os.remove(out)  # augment appends
    real = aug.mzml
    aug.mzml = _RecordReader(records)
    try:
        return aug.augment_spectra(str(d_folder), temp_path, top_n=top_n, scan_width=scan_width,
                                   max_pep_length=max_pep_length, max_charge=max_charge)
    finally:
        aug.mzml = real
