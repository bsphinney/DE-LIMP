#!/usr/bin/env python3
"""
The precursor m/z search range must follow the acquisition.

Written with the fix for SKILL_OPEN_DEFECTS #1: estimate_params.py emitted a
hardcoded 380-980 for every dataset, tagged as a universal default. On a UC Davis
two-platform pilot that matched neither instrument --

    timsTOF HT dia-PASEF   acquired 299.5-1200.5  -> lost 300-380 and 980-1200
    Orbitrap Lumos DIA     acquired  350 -1200    -> lost 350-380 and 980-1200

-- on a run whose stated purpose was measuring proteome depth. It does not error
and the loss is invisible in the output; a human caught it by reading the emitted
parameters against the instrument method. Nothing in CI would have.

Stdlib only (unittest + sqlite3), so CI needs no new dependencies.
Run: python3 -m unittest discover -s skill/ucdavis-proteomics-core-pipeline/tests
"""
import json
import os
import sqlite3
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import detect_acquisition as da          # noqa: E402
from estimate_params import build_diann  # noqa: E402


def make_bruker_d(tmp, lo=299.5, hi=1200.5, width=25.0):
    """A minimal dia-PASEF .d: the tables the detector keys on, nothing else.

    Schema mirrors a real analysis.tdf, verified 2026-08-17: the isolation
    windows live in DiaFrameMsMsWindows (IsolationMz, IsolationWidth), while
    DiaFrameMsMsWindowGroups carries only an Id. Reading the range from the
    wrong one of those two returns nothing.
    """
    d = os.path.join(tmp, "synthetic.d")
    os.makedirs(d, exist_ok=True)
    con = sqlite3.connect(os.path.join(d, "analysis.tdf"))
    cur = con.cursor()
    cur.execute("CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER, WindowGroup INTEGER)")
    cur.execute("CREATE TABLE DiaFrameMsMsWindowGroups (Id INTEGER)")
    cur.execute("CREATE TABLE DiaFrameMsMsWindows ("
                "WindowGroup INTEGER, ScanNumBegin INTEGER, ScanNumEnd INTEGER, "
                "IsolationMz REAL, IsolationWidth REAL, CollisionEnergy REAL)")
    centre_lo, centre_hi = lo + width / 2, hi - width / 2
    n = 6
    step = (centre_hi - centre_lo) / (n - 1)
    cur.executemany(
        "INSERT INTO DiaFrameMsMsWindows VALUES (?,?,?,?,?,?)",
        [(1, 0, 100, centre_lo + i * step, width, 30.0) for i in range(n)])
    con.commit()
    con.close()
    return d


MZML = """<?xml version="1.0"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
 <run>
  <spectrumList count="{n}">
  {spectra}
  </spectrumList>
 </run>
</mzML>
"""
SPEC = """
   <spectrum index="{i}">
    <cvParam accession="MS:1000511" name="ms level" value="2"/>
    <precursorList>
     <precursor><isolationWindow>
      <cvParam accession="MS:1000827" name="isolation window target m/z" value="{tgt}"/>
      <cvParam accession="MS:1000828" name="isolation window lower offset" value="{lo}"/>
      <cvParam accession="MS:1000829" name="isolation window upper offset" value="{hi}"/>
     </isolationWindow></precursor>
    </precursorList>
   </spectrum>
"""


def make_mzml(tmp, lo=350.0, hi=1200.0, width=20.0):
    centres = [lo + width / 2 + i * (hi - lo - width) / 9 for i in range(10)]
    spectra = "".join(SPEC.format(i=i, tgt=c, lo=width / 2, hi=width / 2)
                      for i, c in enumerate(centres))
    p = os.path.join(tmp, "synthetic.mzML")
    with open(p, "w") as fh:
        fh.write(MZML.format(n=len(centres), spectra=spectra))
    return p


class TestDetectionReturnsRange(unittest.TestCase):
    def test_bruker_d_reports_the_acquired_range(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_bruker_d(tmp, lo=299.5, hi=1200.5)
            kind, conf, why, rng = da.detect_bruker_d(d)
            self.assertEqual(kind, "DIA")
            self.assertIsNotNone(rng, "the range must be read, not discarded")
            self.assertAlmostEqual(rng[0], 299.5, places=3)
            self.assertAlmostEqual(rng[1], 1200.5, places=3)

    def test_mzml_reports_the_acquired_range(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = make_mzml(tmp, lo=350.0, hi=1200.0)
            kind, conf, why, rng = da.detect_mzml(p)
            self.assertEqual(kind, "DIA")
            self.assertIsNotNone(rng)
            self.assertAlmostEqual(rng[0], 350.0, places=1)
            self.assertAlmostEqual(rng[1], 1200.0, places=1)

    def test_classify_surfaces_the_range(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_bruker_d(tmp)
            self.assertEqual(da.classify(d)["precursor_mz_range"], [299.5, 1200.5])

    def test_every_detector_returns_four_elements(self):
        # The detectors are unpacked positionally in classify(); a 3-tuple from
        # any branch is a TypeError at runtime that no import-time check catches.
        with tempfile.TemporaryDirectory() as tmp:
            empty = os.path.join(tmp, "empty.d")
            os.makedirs(empty)
            for res in (da.detect_bruker_d(empty),
                        da.detect_mzml(os.path.join(tmp, "missing.mzML")),
                        da.detect_thermo_raw(os.path.join(tmp, "missing.raw"))):
                self.assertEqual(len(res), 4, f"detector returned {len(res)} elements")

    def test_range_is_none_when_unreadable_rather_than_guessed(self):
        with tempfile.TemporaryDirectory() as tmp:
            empty = os.path.join(tmp, "empty.d")
            os.makedirs(empty)
            self.assertIsNone(da.detect_bruker_d(empty)[3])


class TestEstimateParamsUsesRange(unittest.TestCase):
    def _rationale(self, mz_range):
        _text, r = build_diann("DIA", "timstof", 15, 15, "Bruker timsTOF",
                               "test", [], {}, None, mz_range)
        return r

    def test_measured_range_is_used_and_rounded_outward(self):
        r = self._rationale([299.5, 1200.5])
        # Outward, never inward: rounding 299.5 up to 300 would clip the edge window.
        self.assertEqual(r["--min-pr-mz"]["value"], 299)
        self.assertEqual(r["--max-pr-mz"]["value"], 1201)
        self.assertIn("measured", r["--min-pr-mz"]["source"])

    def test_measured_range_is_not_tagged_as_a_universal_default(self):
        # The original defect was not only the number -- it was labelling a guess
        # 'universal trypsin/LFQ default', which reads as deliberate.
        r = self._rationale([299.5, 1200.5])
        self.assertNotIn("universal", r["--min-pr-mz"]["source"].lower())

    def test_fallback_is_tagged_as_a_fallback(self):
        r = self._rationale(None)
        self.assertEqual(r["--min-pr-mz"]["value"], 380)
        self.assertEqual(r["--max-pr-mz"]["value"], 980)
        for k in ("--min-pr-mz", "--max-pr-mz"):
            self.assertIn("FALLBACK", r[k]["source"])
            self.assertNotIn("universal", r[k]["source"].lower())

    def test_orbitrap_range_survives_too(self):
        r = self._rationale([350.0, 1200.0])
        self.assertEqual(r["--min-pr-mz"]["value"], 350)
        self.assertEqual(r["--max-pr-mz"]["value"], 1200)


class TestEndToEndCLI(unittest.TestCase):
    def test_cli_emits_the_measured_range(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "diann.cfg")
            res = subprocess.run(
                [sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                 "--engine", "diann", "--acquisition", "DIA",
                 "--instrument", "timsTOF HT",
                 "--precursor-mz-range", "299.5", "1200.5", "--out", out],
                capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, res.stderr)
            with open(out) as fh:
                cfg = fh.read()
            self.assertIn("--min-pr-mz 299", cfg)
            self.assertIn("--max-pr-mz 1201", cfg)
            self.assertNotIn("--min-pr-mz 380", cfg)

    def test_detect_cli_reports_range_in_json(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_bruker_d(tmp)
            res = subprocess.run(
                [sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"), d],
                capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, res.stderr)
            payload = json.loads(res.stdout)
            self.assertEqual(payload["precursor_mz_range"], [299.5, 1200.5])
            self.assertFalse(payload["precursor_mz_range_mixed"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
