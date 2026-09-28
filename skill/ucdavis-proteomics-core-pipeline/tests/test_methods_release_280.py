#!/usr/bin/env python3
"""2.8.0 release review of make_methods.py: the acknowledgment, unreadable runs, tags, encoding.

What these tests pin:
  * BLOCKER: the acknowledgment follows the instrument NAME, across every registry entry, before
    any filename rule. A real timsTOF HT run renamed FLAG_IP_1.d got the Fusion Lumos S10 grant,
    and Exp3_HeLa.d / ExcisedBand1.d the Exploris one. The Thermo prefix (FL*, Ex*) is only a
    fallback for .raw files whose instrument nothing names -- never for a .d, and never over a
    named instrument (a foreign FLAG_pulldown.raw from an Astral).
  * a .d that cannot be read is named on stderr, in the header and in the note under Mass
    spectrometry, and its blanks are tagged "[raw file not readable here — confirm]"
  * a value no record holds is tagged "[not recorded — confirm]", never a facility default
  * methods.md is written UTF-8 whatever the locale ("1/K₀" fails under cp1252/latin-1) and
    atomically: a failed write leaves the old file, not a 0-byte one
  * nits: one "CaptiveSpray" spelling, "Thermo" not named twice, a ddaPASEF precursor-selection
    placeholder, STAN's 50 °C column temperature as a tagged facility default
stdlib only, no network.
"""
import os
import re
import sqlite3
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_methods as mm                            # noqa: E402
from synthetic_tdf import synthetic_tdf_write_uri   # noqa: E402  (the deliberate writer)
from test_instrument_methods import make_timstof_d, run_methods, section   # noqa: E402

PY = sys.executable
LUMOS, EXPLORIS, TIMS = ("Thermo Orbitrap Fusion Lumos", "Thermo Orbitrap Exploris 480",
                         "Bruker timsTOF")


def touch(path):
    with open(path, "wb") as fh:
        fh.write(b"\0" * 64)
    return path


def run(tmp, raws, *extra, env=None):
    out = os.path.join(tmp, "methods.md")
    return subprocess.run([PY, os.path.join(SCRIPTS, "make_methods.py"), "--raw", *raws,
                           "--out", out, *extra], capture_output=True, text=True, env=env)


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


class Acknowledgment(unittest.TestCase):

    def test_the_instrument_name_wins_over_a_filename(self):
        for name in ("FLAG_IP_1.d", "Exp3_HeLa.d", "ExcisedBand1.d"):
            self.assertEqual(mm.pick_ack("timsTOF HT", [f"/data/{name}"])[0], TIMS, name)

    def test_a_d_never_takes_the_thermo_prefix(self):
        for name in ("FLAG_IP_1.d", "Exp3_HeLa.d", "FL_run.d"):
            label, text = mm.pick_ack(None, [f"/data/{name}"])
            self.assertIsNone(label, name)
            self.assertIn("not in the UC Davis acknowledgment registry", text)

    def test_a_foreign_fl_raw_with_a_named_instrument_is_not_a_lumos(self):
        label, text = mm.pick_ack("Orbitrap Astral", ["/data/FLAG_pulldown.raw"])
        self.assertIsNone(label)
        self.assertNotIn("S10OD021801", text)

    def test_the_prefix_names_a_facility_raw_whose_instrument_nothing_names(self):
        self.assertEqual(mm.pick_ack(None, ["/d/FL20240101_a.raw", "/d/FL20240101_b.raw"])[0],
                         LUMOS)
        self.assertEqual(mm.pick_ack(None, ["/d/Ex_A1.raw"])[0], EXPLORIS)
        self.assertIsNone(mm.pick_ack(None, ["/d/FL_a.raw", "/d/Ex_b.raw"])[0])  # ambiguous
        self.assertIsNone(mm.pick_ack(None, ["/d/FL_a.raw", "/d/FL_b.d"])[0])     # not all .raw

    def test_a_renamed_timstof_run_gets_the_timstof_acknowledgment(self):
        with tempfile.TemporaryDirectory() as tmp:
            text, params = run_methods(tmp, [make_timstof_d(os.path.join(tmp, "FLAG_IP_1.d"))])
            ack = section(text, "Acknowledgments")
            self.assertIn("Howard Hughes Medical Institute", ack)
            self.assertNotIn("S10OD021801", ack)
            self.assertEqual(params["acknowledgment_for"], TIMS)

    def test_a_foreign_raw_keeps_its_recorded_instrument(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw = touch(os.path.join(tmp, "FLAG_pulldown.raw"))
            r = run(tmp, [raw], "--instrument", "Orbitrap Astral")
            self.assertEqual(r.returncode, 0, r.stderr)
            text = read(os.path.join(tmp, "methods.md"))
            self.assertIn("on an Orbitrap Astral mass spectrometer", text)
            self.assertNotIn("Lumos", text)
            self.assertNotIn("S10OD021801", section(text, "Acknowledgments"))

    def test_a_facility_raw_is_named_by_prefix_as_a_guess(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = run(tmp, [touch(os.path.join(tmp, "FL20260101_a.raw"))])
            self.assertEqual(r.returncode, 0, r.stderr)
            text = read(os.path.join(tmp, "methods.md"))
            ms = section(text, "Mass spectrometry")
            self.assertIn("on an Orbitrap Fusion Lumos mass spectrometer (Thermo Fisher "
                          "Scientific)", ms)                   # "Thermo" not named twice
            self.assertIn("S10OD021801", section(text, "Acknowledgments"))
            self.assertIn("a guess — confirm", text)


class UnreadableRun(unittest.TestCase):

    def bad_d(self, tmp, name="broken.d"):
        d = os.path.join(tmp, name)
        os.makedirs(d)
        with open(os.path.join(d, "analysis.tdf"), "wb") as fh:
            fh.write(b"this is not a sqlite database " * 200)
        return d

    def test_the_only_run_unreadable_is_said_and_tagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = run(tmp, [self.bad_d(tmp)], "--instrument", "timsTOF HT")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("broken.d could not be read", r.stderr)
            text = read(os.path.join(tmp, "methods.md"))
            self.assertIn("1 could not be read", text.splitlines()[2])
            ms = section(text, "Mass spectrometry")
            self.assertIn(mm.UNREADABLE_TAG, ms)
            self.assertIn("broken.d could not be read", ms)          # the resolve-before note
            self.assertNotIn(f"____ {mm.DEF}", text)

    def test_a_readable_run_represents_the_series_and_the_bad_one_is_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [self.bad_d(tmp, "a_broken.d"), make_timstof_d(os.path.join(tmp, "b.d"))]
            r = run(tmp, raws)
            self.assertEqual(r.returncode, 0, r.stderr)
            text = read(os.path.join(tmp, "methods.md"))
            self.assertIn("(2 file(s); 1 could not be read", text)
            ms = section(text, "Mass spectrometry")
            self.assertIn("m/z 100–1700", ms)                         # from the good run
            self.assertIn("a_broken.d could not be read", ms)
            self.assertNotIn("all other values were extracted", text)


class Tags(unittest.TestCase):

    def test_a_value_no_record_holds_is_not_recorded_not_a_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_timstof_d(os.path.join(tmp, "x.d"))
            con = sqlite3.connect(synthetic_tdf_write_uri(os.path.join(d, "analysis.tdf")), uri=True)
            con.execute("DELETE FROM GlobalMetadata WHERE Key LIKE 'OneOverK0AcqRange%'")
            con.commit()
            con.close()
            text, _ = run_methods(tmp, [d])
            ms = section(text, "Mass spectrometry")
            self.assertIn(f"1/K₀ ____ {mm.NR_TAG}–____ {mm.NR_TAG} V·s/cm²", ms)
            self.assertNotIn(f"____ {mm.DEF}", text)

    def test_captivespray_one_spelling_and_the_column_temperature_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            text, _ = run_methods(tmp, [make_timstof_d(os.path.join(tmp, "x.d"))])
            prose = section(text, "Liquid chromatography") + section(text, "Mass spectrometry")
            self.assertNotIn("Captive Spray", prose)
            self.assertIn("via a CaptiveSpray ion source", prose)
            self.assertIn(f"column temperature of 50 °C {mm.DEF}", prose)

    def test_ddapasef_says_precursor_selection_is_not_extracted(self):
        ms = mm.ms_paragraph({"mode": "ddaPASEF", "instrument": "timsTOF HT"},
                             lambda x, unit="", default=None: x or f"____ {mm.NR_TAG}")
        self.assertIn(f"Precursor selection (PASEF ramps per cycle, target intensity, charge and "
                      f"mobility filters, dynamic exclusion) was ____ {mm.DDA_TAG}.", ms)
        dia = mm.ms_paragraph({"mode": "dia-PASEF"}, lambda x, unit="", default=None: x or "")
        self.assertNotIn("Precursor selection", dia)


class Encoding(unittest.TestCase):

    def test_every_text_write_names_utf8(self):
        with open(os.path.join(SCRIPTS, "make_methods.py"), encoding="utf-8") as fh:
            src = fh.read()
        for m in re.finditer(r"open\(([^()]*(\([^()]*\))?[^()]*)\)", src):
            args = m.group(1)
            if re.search(r"""["'][wa]b?["']""", args) and '"wb"' not in args:
                self.assertIn('encoding="utf-8"', args, m.group(0))

    def test_a_latin1_locale_still_writes_the_file_whole(self):
        env = dict(os.environ, LC_ALL="en_US.ISO8859-1", LANG="en_US.ISO8859-1", PYTHONUTF8="0")
        probe = subprocess.run([PY, "-c", "import locale;print(locale.getpreferredencoding(False))"],
                               capture_output=True, text=True, env=env).stdout.strip().lower()
        if "utf" in probe:
            self.skipTest(f"this system gives Python {probe} under a latin-1 locale")
        with tempfile.TemporaryDirectory() as tmp:
            r = run(tmp, [make_timstof_d(os.path.join(tmp, "x.d"))], env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            text = read(os.path.join(tmp, "methods.md"))
            self.assertIn("1/K₀ 0.70–1.30 V·s/cm²", text)
            self.assertIn("## Acknowledgments", text)
            self.assertFalse(os.path.exists(os.path.join(tmp, "methods.md.part")))

    def test_a_failed_write_keeps_the_old_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "methods.md")
            mm._write_text(p, "old ✓\n")
            with self.assertRaises(UnicodeEncodeError):
                mm._write_text(p, "new \ud800\n")                    # not encodable
            self.assertEqual(read(p), "old ✓\n")
            self.assertEqual(os.listdir(tmp), ["methods.md"])


if __name__ == "__main__":
    unittest.main()
