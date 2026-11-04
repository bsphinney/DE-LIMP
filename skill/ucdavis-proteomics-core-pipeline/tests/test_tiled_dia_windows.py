#!/usr/bin/env python3
"""
A tiled, repeated set of MS2 windows with no data-dependent scan is DIA "high", whatever its
width or window count (S2, T41).

detect_acquisition.py judged a run by window WIDTH first: >= 3 m/z with at most 80 distinct
centres was DIA, <= 2 m/z DDA. Narrow-window DIA (Orbitrap Astral, e.g. 2 m/z x 300 windows over
380-980) therefore read as DDA by width, rescued to DIA "medium" by the Thermo data-dependent
flag; 3-4 m/z methods with more than 80 windows (Astral 4 m/z over 380-980 is 150) read as
"unknown"; an mzML, which carries no flag the reader uses, read 2 m/z DIA as DDA "high".

The shape a DIA method has and a DDA run does not: ONE fixed set of windows that tiles the
precursor range, each window acquired again every cycle at exactly the same edges
(detect_acquisition.tiled_window_set). A DDA run isolates each precursor at its own measured m/z,
so its windows do not repeat; a PRM list repeats but does not tile. Now:
  - a tiled, repeated set is DIA, past the 80-window cap;
  - with the Thermo flag read and no MS2 scan data-dependent, it is "high";
  - without the flag (mzML), a set of narrow (< 3 m/z) or 3-4 m/z windows is "medium" (the
    user confirms), and a >= 4 m/z set "high" as before.

Synthetic window tables; no parser involved. Hermetic.
"""
import os
import random
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import detect_acquisition as da  # noqa: E402


def tiled(width, lo, n, cycles, dependent=False, overlap=0.0):
    """A Thermo window table (thermo_isolation_windows' shape): `n` windows of `width` m/z from
    `lo`, abutting (or overlapping by `overlap`), every one acquired each of `cycles` cycles."""
    step = width - overlap
    centres = [lo + width / 2 + step * i for i in range(n)]
    w = {"widths": [], "centres": [], "los": [], "his": [],
         "n_ms2": 0, "n_filter": 0, "n_dependent": 0}
    for _ in range(cycles):
        for c in centres:
            w["widths"].append(width); w["centres"].append(c)
            w["los"].append(c - width / 2); w["his"].append(c + width / 2)
            w["n_ms2"] += 1; w["n_filter"] += 1; w["n_dependent"] += int(dependent)
    return w


def dda(n=1200, width=1.6, seed=3, reselect=0.25):
    """A DDA table: each MS2 isolates a precursor at its measured m/z (4 decimals), a quarter of
    them re-selected later at that m/z +/- 3 ppm (dynamic exclusion over), all data-dependent."""
    rnd = random.Random(seed)
    picked = []
    for _ in range(n):
        if picked and rnd.random() < reselect:
            mz = rnd.choice(picked) * (1 + rnd.uniform(-3e-6, 3e-6))
        else:
            mz = rnd.uniform(400.0, 1100.0)
        picked.append(mz)
    w = {"widths": [], "centres": [], "los": [], "his": [],
         "n_ms2": n, "n_filter": n, "n_dependent": n}
    for mz in picked:
        c = round(mz, 4)
        w["widths"].append(width); w["centres"].append(c)
        w["los"].append(c - width / 2); w["his"].append(c + width / 2)
    return w


def shared(w):
    """The mzML reader's verdict: the width/centre rule alone, no data-dependent flag."""
    return da.classify_isolation_windows(w["widths"], w["centres"], w["los"], w["his"])


class TiledDIA(unittest.TestCase):
    def test_astral_2mz_x_300_windows_is_dia_high(self):
        kind, conf, why, rng = da.classify_thermo_windows(tiled(2.0, 380.0, 300, 3))
        self.assertEqual((kind, conf), ("DIA", "high"), why)
        self.assertEqual(rng, (380.0, 980.0))
        self.assertIn("300 windows", why)

    def test_astral_4mz_x_150_windows_is_dia_high(self):
        """Past the old 80-window cap: it read as "unknown" by width, DIA "medium" by flag."""
        kind, conf, why, rng = da.classify_thermo_windows(tiled(4.0, 380.0, 150, 4))
        self.assertEqual((kind, conf), ("DIA", "high"), why)
        self.assertEqual(rng, (380.0, 980.0))

    def test_3mz_x_200_overlapping_windows_is_dia_high(self):
        kind, conf, why, rng = da.classify_thermo_windows(
            tiled(3.0, 400.0, 200, 3, overlap=0.5))
        self.assertEqual((kind, conf), ("DIA", "high"), why)
        self.assertEqual(rng[0], 400.0)

    def test_an_mzml_of_2mz_dia_is_no_longer_called_dda(self):
        """No flag to read: DIA, for the user to confirm -- it read as DDA "high", and a DDA
        verdict routes the run to Sage without asking."""
        kind, conf, why, rng = shared(tiled(2.0, 380.0, 300, 3))
        self.assertEqual((kind, conf), ("DIA", "medium"), why)
        self.assertEqual(rng, (380.0, 980.0))

    def test_an_mzml_of_4mz_x_150_windows_is_dia(self):
        kind, conf, why, rng = shared(tiled(4.0, 380.0, 150, 4))
        self.assertEqual((kind, conf), ("DIA", "high"), why)

    def test_the_shape_itself(self):
        self.assertEqual(da.tiled_window_set(*_edges(tiled(2.0, 380.0, 300, 3))), 300)
        self.assertIsNone(da.tiled_window_set(*_edges(tiled(2.0, 380.0, 300, 1))),
                          "one cycle: nothing repeated")
        self.assertIsNone(da.tiled_window_set(*_edges(dda())))


class Pins(unittest.TestCase):
    def test_a_4mz_x_75_window_table_stays_dia_high(self):
        """(pin) inside the old cap: high by width, with or without the flag."""
        w = tiled(4.0, 400.0, 75, 4)
        self.assertEqual(da.classify_thermo_windows(w)[:2], ("DIA", "high"))
        self.assertEqual(shared(w)[:2], ("DIA", "high"))

    def test_a_dda_table_stays_dda(self):
        """(pin) re-selected precursors repeat only within a few ppm, never at the same edges,
        so the windows are not one repeated set."""
        w = dda()
        self.assertEqual(da.classify_thermo_windows(w)[:2], ("DDA", "high"))
        self.assertEqual(shared(w)[:2], ("DDA", "high"))
        self.assertIsNone(da.classify_thermo_windows(w)[3])

    def test_a_dda_table_with_two_activations_per_precursor_stays_dda(self):
        """A DDA method that fragments each precursor twice (HCD then EThcD, or CID and HCD)
        acquires two MS2 at the same edges, one after the other: its windows recur, and dense
        precursors abut -- but not once per cycle. Without the flag (mzML) it stays DDA."""
        w = dda(n=1500, reselect=0.0)       # 3,000 MS2: the mzML reader's cap; dense enough to abut
        for k in ("widths", "centres", "los", "his"):
            w[k] = [x for x in w[k] for _ in (0, 1)]
        w["n_ms2"] = w["n_filter"] = w["n_dependent"] = len(w["widths"])
        self.assertIsNone(da.tiled_window_set(w["los"], w["his"]))
        self.assertEqual(shared(w)[0], "DDA")
        self.assertEqual(da.classify_thermo_windows(w)[:2], ("DDA", "high"))

    def test_repeated_targets_that_do_not_tile_are_still_asked_about(self):
        """(pin) a PRM-like list -- 2 m/z windows on 40 scattered targets, every cycle, nothing
        data-dependent -- is not DIA "high": "narrow-window DIA or PRM?", medium."""
        rnd = random.Random(5)
        targets = sorted(rnd.uniform(400.0, 1000.0) for _ in range(40))
        w = {"widths": [], "centres": [], "los": [], "his": [],
             "n_ms2": 0, "n_filter": 0, "n_dependent": 0}
        for _ in range(5):
            for c in targets:
                w["widths"].append(2.0); w["centres"].append(c)
                w["los"].append(c - 1.0); w["his"].append(c + 1.0)
                w["n_ms2"] += 1; w["n_filter"] += 1
        kind, conf, why, _ = da.classify_thermo_windows(w)
        self.assertEqual((kind, conf), ("DIA", "medium"), why)
        self.assertIn("PRM", why)

    def test_a_tiled_set_that_is_all_data_dependent_is_not_called_dia(self):
        """Contradicting evidence is asked about: the instrument's own flag outranks the window
        shape, and the verdict is not "high" (by width alone it was DDA "high")."""
        kind, conf, why, _ = da.classify_thermo_windows(tiled(2.0, 380.0, 300, 3,
                                                              dependent=True))
        self.assertEqual(kind, "DDA", why)
        self.assertNotEqual(conf, "high", why)


class TheReasonDoesNotContradictItself(unittest.TestCase):
    """Review F9: the shared rule said "no data-dependent flag was read to confirm it", and the
    Thermo reader then added "0/900 MS2 scans flagged data-dependent" to the same reason. The
    clause belongs to the reader that reads no flag (mzML) only."""

    NO_FLAG = "no data-dependent flag was read"

    def test_the_thermo_reason_states_the_flag_once(self):
        for dep in (False, True):
            why = da.classify_thermo_windows(tiled(2.0, 380.0, 300, 3, dependent=dep))[2]
            self.assertNotIn(self.NO_FLAG, why, why)
            self.assertIn("MS2 scans flagged data-dependent", why)

    def test_the_mzml_reason_says_no_flag_was_read(self):
        why = shared(tiled(2.0, 380.0, 300, 3))[2]
        self.assertIn(self.NO_FLAG, why)


def _edges(w):
    return w["los"], w["his"]


if __name__ == "__main__":
    unittest.main()
