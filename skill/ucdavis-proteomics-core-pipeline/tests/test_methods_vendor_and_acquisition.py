#!/usr/bin/env python3
"""Methods name the instrument's own vendor, the recorded DIA/DDA, and no Core LC defaults for
data from elsewhere (2.14.2 triage S7, I09).

The coral study's session (2026-10-22) wrote Methods for another lab's mzML from a Bruker impact II
that called it "(Thermo Fisher Scientific)" and printed the Core's PepSep column and CaptiveSpray
emitter. A Thermo cohort got the Bruker PepSep column and "not recorded in the .d", and the step
9d run said "[DDA/DIA — confirm]" although step 2's detection was in the workflow manifest.

What these tests pin:
  * instrument_vendor() is THE vendor rule: words of the instrument name, Bruker or Thermo, else
    None -- and a None vendor is printed "[vendor — confirm]", never guessed
  * a Bruker that is not a timsTOF gets the plain MS sentence with "(Bruker Daltonics)", no TIMS
  * the Core's LC defaults (PepSep column, CaptiveSpray emitter) are for a Core timsTOF .d only;
    a Thermo .raw has the column blank and tagged; outside data gets tagged blanks, no defaults
  * the acquisition mode is --acquisition, else the workflow manifest's, else the session's
    input/acquisition.json, else the tagged "[DDA/DIA — confirm]"
stdlib only, no network.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_methods as mm                                        # noqa: E402
from test_instrument_methods import make_timstof_d, section     # noqa: E402

PY = sys.executable


def touch(path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        fh.write(b"\0" * 64)
    return path


def methods(tmp, raw, instrument=None, acquisition=None, workflow_manifest=None, out=None):
    """make_methods.py end to end; returns (methods.md text, methods_params.json)."""
    out = out or os.path.join(tmp, "methods.md")
    cmd = [PY, os.path.join(SCRIPTS, "make_methods.py"), "--raw", *raw, "--out", out]
    if instrument:
        cmd += ["--instrument", instrument]
    if acquisition:
        cmd += ["--acquisition", acquisition]
    if workflow_manifest is not None:
        wf = os.path.join(tmp, "workflow.manifest.json")
        with open(wf, "w", encoding="utf-8") as fh:
            json.dump(workflow_manifest, fh)
        cmd += ["--workflow-manifest", wf]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(r.stderr)
    with open(out, encoding="utf-8") as fh:
        text = fh.read()
    with open(os.path.splitext(out)[0] + "_params.json", encoding="utf-8") as fh:
        return text, json.load(fh)


class VendorRule(unittest.TestCase):

    def test_bruker_names(self):
        for name in ("timsTOF HT", "timsTOF Pro 2", "impact II", "maXis 4G", "compact",
                     "timsTOF-flex"):
            self.assertEqual(mm.instrument_vendor(name), "Bruker", name)

    def test_thermo_names(self):
        for name in ("Orbitrap Exploris 480", "Q Exactive HF-X", "Orbitrap Fusion Lumos",
                     "Orbitrap Eclipse", "Orbitrap Astral", "Orbitrap Ascend", "LTQ Orbitrap Velos",
                     "Exploris480"):
            self.assertEqual(mm.instrument_vendor(name), "Thermo", name)

    def test_unknown_and_empty_are_none(self):
        for name in (None, "", "Mystery 9000", "ZenoTOF 7600", "compactness tester"):
            self.assertIsNone(mm.instrument_vendor(name), name)


class Vendor(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.mzml = touch(os.path.join(self.tmp, "data", "coral_oocyte_1.mzML"))
        self.raw = touch(os.path.join(self.tmp, "data", "sample1.raw"))

    def tearDown(self):
        self._tmp.cleanup()

    def test_mzml_from_a_bruker_impact_is_not_called_thermo(self):
        md, params = methods(self.tmp, [self.mzml], instrument="impact II")
        ms = section(md, "Mass spectrometry")
        self.assertIn("on an impact II mass spectrometer (Bruker Daltonics)", ms)
        self.assertNotIn("Thermo Fisher Scientific", md)
        self.assertNotIn("TIMS", ms)                          # an impact II has no TIMS
        self.assertEqual(params["representative"]["vendor"], "Bruker")

    def test_unknown_vendor_is_tagged_not_guessed(self):
        md, _ = methods(self.tmp, [self.mzml], instrument="Mystery 9000")
        self.assertIn("[vendor — confirm]", md)
        self.assertNotIn("Thermo Fisher", md)
        self.assertNotIn("Bruker Daltonics", md)

    def test_outside_data_gets_no_core_lc_defaults(self):
        md, params = methods(self.tmp, [self.mzml], instrument="impact II")
        self.assertNotIn("PepSep", md)
        self.assertNotIn("CaptiveSpray emitter", md)
        lc = section(md, "Liquid chromatography")
        self.assertNotIn(mm.DEF, lc)                          # tagged blanks, no Core defaults
        self.assertIn(mm.NR_TAG, lc)
        self.assertNotEqual(params["column"]["tag"], mm.DEF)

    def test_thermo_gets_no_bruker_column_and_says_raw_file(self):
        md, params = methods(self.tmp, [self.raw], instrument="Orbitrap Exploris 480")
        self.assertNotIn("PepSep", md)
        self.assertNotIn("not recorded in the .d", md)
        self.assertIn("not recorded in the raw file", md)
        self.assertIn("(Thermo Fisher Scientific)", section(md, "Mass spectrometry"))
        self.assertEqual(params["column"]["tag"], mm.NR_TAG)

    def test_a_bruker_d_read_from_the_record_is_not_a_timstof_unless_named(self):
        # an impact .d holds analysis.baf, not analysis.tdf: it is not read here, so the
        # instrument and vendor come from the record -- and it gets no TIMS sentence
        d = os.path.join(self.tmp, "data", "run1.d")
        os.makedirs(d)
        md, params = methods(self.tmp, [d], instrument="impact II")
        ms = section(md, "Mass spectrometry")
        self.assertIn("(Bruker Daltonics)", ms)
        self.assertNotIn("TIMS", ms)
        self.assertNotIn("PepSep", md)

    def test_core_timstof_d_keeps_the_facility_defaults(self):          # (pin)
        d = make_timstof_d(os.path.join(self.tmp, "data", "r1.d"))
        md, params = methods(self.tmp, [d])
        self.assertIn(f"{mm.LC_COLUMN_DEFAULT} {mm.DEF}", section(md, "Liquid chromatography"))
        self.assertIn(f"{mm.EMITTER_DEFAULT} {mm.DEF}", section(md, "Mass spectrometry"))
        self.assertIn("not recorded in the .d", md)

    def test_core_timstof_from_the_record_keeps_the_facility_defaults(self):   # (pin)
        d = os.path.join(self.tmp, "data", "r1.d")
        os.makedirs(d)
        md, _ = methods(self.tmp, [d], instrument="timsTOF HT")
        self.assertIn(f"{mm.LC_COLUMN_DEFAULT} {mm.DEF}", section(md, "Liquid chromatography"))
        self.assertIn("TIMS", section(md, "Mass spectrometry"))


class TimsTofOutsideData(unittest.TestCase):
    """2.14.2 review: a timsTOF mzML is described with the timsTOF paragraph, which read its mode
    from the raw file only; and outside data printed the Core's acknowledgment untagged."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.mzml = touch(os.path.join(self.tmp, "data", "hela_1.mzML"))

    def tearDown(self):
        self._tmp.cleanup()

    def test_a_timstof_mzml_states_the_recorded_mode(self):
        md, params = methods(self.tmp, [self.mzml], instrument="timsTOF Pro 2",
                             workflow_manifest={"acquisition": "DIA"})
        ms = section(md, "Mass spectrometry")
        self.assertIn("operated in DIA mode", ms)
        self.assertNotIn(f"operated in ____ {mm.NR_TAG} mode", ms)

    def test_outside_data_tags_the_core_acknowledgment(self):
        md, _ = methods(self.tmp, [self.mzml], instrument="timsTOF Pro 2")
        ack = section(md, "Acknowledgments")
        self.assertIn("Howard Hughes Medical Institute", ack)          # kept, but tagged
        self.assertIn(mm.ACK_OUTSIDE_TAG, ack)

    def test_an_outside_thermo_mzml_is_tagged_too(self):
        md, _ = methods(self.tmp, [self.mzml], instrument="Orbitrap Exploris 480")
        ack = section(md, "Acknowledgments")
        self.assertIn("S10OD026918-01A1", ack)
        self.assertIn(mm.ACK_OUTSIDE_TAG, ack)

    def test_a_core_raw_is_not_tagged(self):                             # (pin)
        raw = touch(os.path.join(self.tmp, "data", "s1.raw"))
        md, _ = methods(self.tmp, [raw], instrument="Orbitrap Exploris 480")
        ack = section(md, "Acknowledgments")
        self.assertIn("S10OD026918-01A1", ack)
        self.assertNotIn(mm.ACK_OUTSIDE_TAG, ack)

    def test_a_core_timstof_d_is_not_tagged(self):                       # (pin)
        md, _ = methods(self.tmp, [make_timstof_d(os.path.join(self.tmp, "data", "r1.d"))])
        self.assertNotIn(mm.ACK_OUTSIDE_TAG, section(md, "Acknowledgments"))


class Acquisition(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.raw = touch(os.path.join(self.tmp, "data", "sample1.raw"))

    def tearDown(self):
        self._tmp.cleanup()

    def test_recorded_acquisition_is_used_without_the_flag(self):
        md, params = methods(self.tmp, [self.raw], workflow_manifest={"acquisition": "DIA"})
        self.assertIn("DIA mode", md)
        self.assertNotIn("[DDA/DIA — confirm]", md)
        self.assertEqual(params["acquisition"], "DIA")
        self.assertIn("workflow manifest", params["acquisition_source"])

    def test_the_flag_outranks_the_manifest(self):
        md, params = methods(self.tmp, [self.raw], acquisition="DDA",
                             workflow_manifest={"acquisition": "DIA"})
        self.assertIn("DDA mode", md)
        self.assertEqual(params["acquisition"], "DDA")
        self.assertIn("--acquisition", params["acquisition_source"])

    def test_the_sessions_acquisition_json_is_read(self):
        sess = os.path.join(self.tmp, "session")
        os.makedirs(os.path.join(sess, "output"))
        with open(touch(os.path.join(sess, "input", "acquisition.json")), "w",
                  encoding="utf-8") as fh:
            json.dump({"overall": "DDA", "instrument": "Orbitrap Exploris 480"}, fh)
        md, params = methods(self.tmp, [self.raw],
                             out=os.path.join(sess, "output", "methods.md"))
        self.assertIn("DDA mode", md)
        self.assertNotIn("[DDA/DIA — confirm]", md)
        self.assertIn("acquisition.json", params["acquisition_source"])

    def test_a_mixed_or_missing_record_stays_tagged(self):
        md, params = methods(self.tmp, [self.raw], workflow_manifest={"acquisition": "mixed"})
        self.assertIn("[DDA/DIA — confirm]", md)
        self.assertIsNone(params["acquisition"])

    def test_acquisition_from_records_order(self):
        class A:
            acquisition = None
            workflow_manifest = None
            out = os.path.join(self.tmp, "nowhere", "methods.md")
            raw = [self.raw]
            de_dir = params = search_prov = None
        self.assertEqual(mm.acquisition_from_records(A()), (None, None))
        A.acquisition = "dia"
        self.assertEqual(mm.acquisition_from_records(A())[0], "DIA")


if __name__ == "__main__":
    unittest.main()
