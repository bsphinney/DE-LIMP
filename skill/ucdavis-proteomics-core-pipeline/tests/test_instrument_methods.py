#!/usr/bin/env python3
"""The LC-MS paragraphs of the Methods must come from the run, in the published order.

PROT_0756 (Dickson lab, timsTOF HT + Evosep One, 2026-08-13) got a Methods section that left
"[LC system / gradient — confirm]" although every .d names the Evosep One and its "60 samples
per day" method; reported m/z 99.993933-1700; said the collision energy was "ramped from ≈26.3
to ≈49.4 eV" when the method ramps 20 eV at 1/K0 0.6 to 65 eV at 1/K0 1.6 (26.3-49.4 are the
per-window values); and never gave the window overlap, the 1/K0 coverage, the cycle time or
the ion-source settings, all of which are recorded. bruker_method.py now reads them.

What these tests pin:
  * the LC system, vendor, method, run time and Evotip tray come from HyStarMetadata.xml /
    SampleInfo.xml; the ion source is named by the .d's own code table, not from memory
  * the dia-PASEF scheme (windows, ramps, width, spacing, overlap, m/z and 1/K0 coverage),
    cycle time and duty cycle are computed from what the .d records
  * the collision-energy ramp is stated only when the per-window energies match it; when
    they do not, the per-window range is stated and tagged
  * the column: --lc-column > HyStar ColumnInfo > a column log matched to the acquisition
    dates > the facility default, tagged; a change inside the series is never papered over
  * what no file records (column temperature, emitter, mobile phases, Evotip loading) is
    always tagged (DE-LIMP rule 2), and runs acquired differently are flagged
  * reading a .d writes nothing into it, and an old .d without side files still works
stdlib only, no network.
"""
import csv
import json
import os
import sqlite3
import struct
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import bruker_method as bm                           # noqa: E402
import make_deposit as md                            # noqa: E402  md_sections(): one reader
import make_methods as mm                            # noqa: E402
from synthetic_tdf import synthetic_tdf_write_uri   # noqa: E402  (the deliberate writer)

PY = sys.executable
NS = "https://www.bruker.com/compass/metadata"
SOURCE_CODES = ("0:No Source;1:ESI;2:APCI;3:Nano ESI Offline;4:Nano ESI Online;5:APPI;"
                "6:Multi Mode;9:Nano Flow ESI;10:ionBooster;11:Captive Spray;12:GC-APCI;"
                "13:VIP-HESI;14:VIP-HESI-APCI;15:DART")


def scheme_windows():
    """The PROT_0756 dia-PASEF scheme: 11 ramps, 36 windows of 26 Th at 25 Th spacing over
    m/z 299.5-1200.5 -- (group, centre)."""
    out = []
    for g in range(1, 12):
        for k in range(4):
            c = 1187.5 - 25 * (g - 1) - 275 * k
            if c >= 312.5:
                out.append((g, c))
    return out


def window_im(mz):
    """A 1/K0 band per window, inside 0.70-1.30, rising with m/z (a synthetic stand-in)."""
    mid = 0.70 + (mz - 300.0) / 900.0 * 0.60
    return max(0.70, mid - 0.06), min(1.30, mid + 0.06)


def blob(*xs):
    return struct.pack("<%dd" % len(xs), *xs)


def make_timstof_d(path, method="DIA_11x3-k07t13Ra85.m", ce_hi=65.0, ramp_hi=65.0,
                   acquired="2026-08-14T01:30:04.726-07:00", column_info="",
                   source_codes=SOURCE_CODES, side_files=True):
    """A synthetic timsTOF HT + Evosep One dia-PASEF .d carrying the metadata PROT_0756's
    runs carry. `ce_hi` sets the window energies' ramp end, `ramp_hi` the method's."""
    os.makedirs(path)
    tdf = os.path.join(path, "analysis.tdf")
    con = sqlite3.connect(synthetic_tdf_write_uri(tdf), uri=True)
    con.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
    con.executemany("INSERT INTO GlobalMetadata VALUES (?,?)", [
        ("InstrumentName", "timsTOF HT"), ("InstrumentSerialNumber", "1895883.10878"),
        ("AcquisitionSoftware", "timsTOF"), ("AcquisitionSoftwareVersion", "6.0.6"),
        ("AcquisitionDateTime", acquired), ("MethodName", method),
        ("MzAcqRangeLower", "99.993933"), ("MzAcqRangeUpper", "1700.000000"),
        ("OneOverK0AcqRangeLower", "0.700000"), ("OneOverK0AcqRangeUpper", "1.300000")])
    con.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, Polarity TEXT, "
                "MsMsType INTEGER, AccumulationTime REAL, RampTime REAL)")
    frames, t = [], 0.0
    for i in range(1, 12 * 5 + 1):                   # five cycles: 1 MS1 + 11 dia-PASEF
        frames.append((i, round(t, 4), "+", 0 if (i - 1) % 12 == 0 else 9, 85.05, 85.05))
        t += 0.0916
    con.executemany("INSERT INTO Frames VALUES (?,?,?,?,?,?)", frames)
    con.execute("CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, "
                "ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, "
                "CollisionEnergy REAL)")
    rows = []
    for g, c in scheme_windows():
        lo, hi = window_im(c)
        rows.append((g, 0, 0, c, 26.0, 20.0 + (ce_hi - 20.0) * ((lo + hi) / 2 - 0.6)))
    con.executemany("INSERT INTO DiaFrameMsMsWindows VALUES (?,?,?,?,?,?)", rows)
    con.execute("CREATE TABLE PropertyDefinitions (Id INTEGER PRIMARY KEY, PermanentName TEXT, "
                "Type INTEGER, DisplayGroupName TEXT, DisplayName TEXT, DisplayValueText TEXT, "
                "DisplayFormat TEXT, DisplayDimension TEXT, Description TEXT)")
    props = [("Source_Type", "11", source_codes, ""),
             ("Source_CapillarySetValue", "1700", "", "V"),
             ("Source_DryGasSetValue", "3.0", "", "l/min"),
             ("Source_DryHeaterSetValue", "200.0", "", "°C"),
             ("Mode_IonPolarity", "0", "0:Positive;1:Negative", ""),
             ("Energy_Ramping_Collision_Energy_Active", "1", "", ""),
             ("Energy_Ramping_Advanced_Settings_Active", "0", "", ""),
             ("Energy_Ramping_Mobility_StartEnd", blob(0.6, 1.6), "", "V·s/cm²"),
             ("Energy_Ramping_Collision_Energy_StartEnd", blob(20.0, ramp_hi), "", "eV"),
             # the inactive advanced list is a different ramp: it must not be read
             ("Energy_Ramping_Advanced_ListMobilityValues", blob(1.6, 0.6), "", "V·s/cm²"),
             ("Energy_Ramping_Advanced_ListCollisionEnergyValues", blob(59.0, 20.0), "", "eV")]
    con.executemany("INSERT INTO PropertyDefinitions VALUES (?,?,0,'','',?,'',?,'')",
                    [(i + 1, n, enum, dim) for i, (n, _, enum, dim) in enumerate(props)])
    con.execute("CREATE TABLE Properties (Frame INTEGER, Property INTEGER, Value)")
    con.executemany("INSERT INTO Properties VALUES (2,?,?)",
                    [(i + 1, val) for i, (_, val, _, _) in enumerate(props)])
    con.commit()
    con.close()
    if not side_files:
        return path
    mdir = os.path.join(path, "23648.m")
    os.makedirs(os.path.join(mdir, "backup-2026-08-14.m"))
    con = sqlite3.connect(synthetic_tdf_write_uri(os.path.join(mdir, "diaSettings.diasqlite")),
                          uri=True)
    con.execute("CREATE TABLE DiaWindowsSpecification (Id INTEGER, Type INTEGER, CycleId INTEGER, "
                "OneOverK0Start REAL, OneOverK0End REAL, IsolationMz REAL, IsolationWidth REAL, "
                "CollisionEnergy REAL)")
    con.execute("INSERT INTO DiaWindowsSpecification VALUES (1,0,0,NULL,NULL,NULL,NULL,NULL)")
    con.executemany("INSERT INTO DiaWindowsSpecification VALUES (?,1,?,?,?,?,26.0,NULL)",
                    [(i + 2, g, *window_im(c), c) for i, (g, c) in enumerate(scheme_windows())])
    con.commit()
    con.close()
    with open(os.path.join(mdir, "microTOFQImpacTemAcquisition.method"), "w",
              encoding="latin-1") as fh:
        fh.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n<root>\n    <fileinfo '
                 'type="ImpacTEM Pro Method" appname="Bruker timsControl" appversion="6.0.0.4" '
                 'createdate="2026-08-14T01:51:05.293-07:00"/>\n</root>\n')
    with open(os.path.join(mdir, "hystar.method"), "w", encoding="utf-8") as fh:
        fh.write('<?xml version="1.0" encoding="utf-8"?>\n<root>\n  <LCMethodData Version="2">'
                 '<TotalRunTime>21</TotalRunTime></LCMethodData>\n'
                 f'  <ColumnInfo>{column_info}\n</ColumnInfo>\n</root>\n')
    with open(os.path.join(path, "HyStarMetadata.xml"), "w", encoding="utf-8") as fh:
        fh.write(f'''<?xml version="1.0"?>
<Metadata xmlns="{NS}">
  <Plugin ID="ms" Name="Bruker OTOF MS"><PluginMetadata>
    <ParameterGroup ID="Connection" Name="Connection">
      <Parameter ID="Product Name" Name="Product Name"><Value>timsTOF HT</Value></Parameter>
      <Parameter ID="MS Control" Name="MS Control"><Value>timsControl</Value></Parameter>
    </ParameterGroup></PluginMetadata></Plugin>
  <Plugin ID="icf" Name="ICF System"><PluginMetadata>
    <Module ID="SAMPLER0" Name="Evosep One (Sampler0)" Type="ALS">
      <ParameterGroup ID="Configuration" Name="Configuration">
        <Parameter ID="ConfigurationObject_DeviceName" Name="Device name"><Value>Evosep One</Value></Parameter>
        <Parameter ID="ConfigurationObject_SerialNumber" Name="Serial number"><Value>S00230</Value></Parameter>
      </ParameterGroup>
      <ParameterGroup ID="Method" Name="Method">
        <Parameter ID="MethodDataObject_RunTimeMinutes" Name="Run time" Unit="min"><Value>21</Value></Parameter>
        <Parameter ID="MethodDataObject_MethodName" Name="Method name"><Value>60 samples per day</Value></Parameter>
      </ParameterGroup>
    </Module>
    <ParameterGroup ID="ICF.Configuration" Name="Configuration">
      <Table ID="ICF.Configuration.Modules.Table" Name="Modules">
        <Col ColIndex="0" Name="Name" /><Col ColIndex="1" Name="Module type" />
        <Col ColIndex="5" Name="Vendor" />
        <Row RowIndex="0"><Cell ColIndex="0"><Value>Evosep One</Value></Cell>
          <Cell ColIndex="1"><Value>EVOSEP_ONE</Value></Cell>
          <Cell ColIndex="5"><Value>Evosep Biosystems</Value></Cell></Row>
      </Table>
    </ParameterGroup></PluginMetadata></Plugin>
</Metadata>
''')
    sample_info = ('<?xml version="1.0" encoding="utf-16"?><SampleTable><SampleTableHeader '
                   'Version="3" HyStarVersion="6.3.1.8" /><Sample Position="S3-F1" '
                   'Method="D:\\Methods\\EvoSepLCmeth\\60spd.m?HyStar_LC" /><SampleTableProperties>'
                   '<Property PropertyNumber="1" Name="TrayType" Value="96Evotip" />'
                   '<Property PropertyNumber="2" Name="HyStar_LC_Method_Name" '
                   'Value="60 samples per day" /></SampleTableProperties></SampleTable>')
    with open(os.path.join(path, "SampleInfo.xml"), "wb") as fh:
        fh.write(sample_info.encode("utf-16"))                  # BOM + UTF-16, as HyStar writes
    return path


def snapshot(d):
    out = {}
    for root, _, files in os.walk(d):
        for f in files:
            p = os.path.join(root, f)
            st = os.stat(p)
            out[os.path.relpath(p, d)] = (st.st_size, st.st_mtime_ns)
    return out


def run_methods(tmp, raws, *extra):
    out = os.path.join(tmp, "methods.md")
    r = subprocess.run([PY, os.path.join(SCRIPTS, "make_methods.py"), "--raw", *raws,
                        "--out", out, *extra], capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(r.stderr)
    with open(out, encoding="utf-8") as fh:
        text = fh.read()
    with open(os.path.join(tmp, "methods_params.json"), encoding="utf-8") as fh:
        params = json.load(fh)
    return text, params


def section(text, name):
    return md.md_sections(text).get(name, "")


class ReadRun(unittest.TestCase):
    """bruker_method.read_run: what the .d records, and where each value came from."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.d = make_timstof_d(os.path.join(cls.tmp.name, "run.d"))
        cls.r = bm.read_run(cls.d)
        cls.v, cls.src = cls.r["values"], cls.r["sources"]

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_lc_system_method_and_tray(self):
        self.assertEqual(self.v["lc_system"], "Evosep One")
        self.assertEqual(self.v["lc_vendor"], "Evosep Biosystems")
        self.assertEqual(self.v["lc_method"], "60 samples per day")
        self.assertEqual(self.v["lc_run_min"], 21.0)
        self.assertEqual(self.v["tray_type"], "96Evotip")
        self.assertEqual(self.v["hystar_version"], "6.3.1.8")
        self.assertIn("HyStarMetadata.xml", self.src["lc_system"])

    def test_software(self):
        self.assertEqual(self.v["ms_control"], "timsControl")
        self.assertEqual(self.v["acquisition_software_version"], "6.0.6")
        self.assertEqual(self.v["control_software"], "Bruker timsControl")
        self.assertEqual(self.v["control_software_version"], "6.0.0.4")

    def test_ion_source_and_polarity(self):
        self.assertEqual(self.v["source_type"], "Captive Spray")
        self.assertEqual(self.v["capillary_v"], {"value": 1700.0, "unit": "V"})
        self.assertEqual(self.v["dry_gas"], {"value": 3.0, "unit": "l/min"})
        self.assertEqual(self.v["dry_temp"]["value"], 200.0)
        self.assertEqual(self.v["polarity"], "positive")

    def test_source_is_named_by_the_files_own_code_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_timstof_d(os.path.join(tmp, "x.d"), source_codes="11:Whatever The File Says")
            self.assertEqual(bm.read_run(d)["values"]["source_type"], "Whatever The File Says")
            d2 = make_timstof_d(os.path.join(tmp, "y.d"), source_codes="1:ESI")
            self.assertEqual(bm.read_run(d2)["values"]["source_type"], "source type code 11")

    def test_timing(self):
        self.assertEqual(self.v["mode"], "dia-PASEF")
        self.assertEqual(self.v["frames_per_cycle"], 12)
        self.assertAlmostEqual(self.v["cycle_s"], 1.10, places=2)
        self.assertEqual((self.v["ramp_ms"], self.v["accumulation_ms"]), (85.05, 85.05))

    def test_window_scheme(self):
        s = bm.window_scheme(self.v)
        self.assertEqual(s["n_windows"], 36)
        self.assertEqual(s["n_ramps"], 11)
        self.assertEqual(s["per_ramp"], (3, 4))
        self.assertEqual(s["width"], 26.0)
        self.assertEqual((s["spacing"], s["overlap"]), (25.0, 1.0))
        self.assertEqual((s["mz_lo"], s["mz_hi"]), (299.5, 1200.5))
        self.assertEqual((round(s["im_lo"], 2), round(s["im_hi"], 2)), (0.70, 1.30))

    def test_ce_ramp_is_the_active_one_and_matches_the_windows(self):
        self.assertEqual(self.v["ce_ramp"]["points"], [(0.6, 20.0), (1.6, 65.0)])
        self.assertFalse(self.v["ce_ramp"]["advanced"])
        self.assertEqual(bm.check_ce_ramp(self.v)["status"], "ok")

    def test_ce_ramp_that_does_not_match_the_windows_is_a_mismatch(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_timstof_d(os.path.join(tmp, "x.d"), ce_hi=59.0)
            chk = bm.check_ce_ramp(bm.read_run(d)["values"])
            self.assertEqual(chk["status"], "mismatch")
            self.assertGreater(chk["max_diff_ev"], bm.CE_CHECK_EV)

    def test_empty_column_info_is_not_a_column(self):
        self.assertIsNone(self.v["column_info"])
        self.assertIn("empty", self.src["column_info"])

    def test_reading_writes_nothing_into_the_d(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_timstof_d(os.path.join(tmp, "x.d"))
            before = snapshot(d)
            bm.read_run(d)
            mm.bruker_meta(d)
            self.assertEqual(snapshot(d), before)

    def test_a_d_without_side_files_loses_only_those_values(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = make_timstof_d(os.path.join(tmp, "x.d"), side_files=False)
            v = bm.read_run(d)["values"]
            self.assertEqual(v["instrument"], "timsTOF HT")
            self.assertEqual(bm.window_scheme(v)["n_windows"], 36)
            for key in ("lc_system", "window_im", "hystar_version", "column_info"):
                self.assertNotIn(key, v)
            self.assertEqual(bm.check_ce_ramp(v)["status"], "unchecked")


class Paragraphs(unittest.TestCase):
    """make_methods.py end to end on synthetic runs."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        raws = [make_timstof_d(os.path.join(cls.tmp.name, f"r{i}.d"),
                               acquired=f"2026-08-13T1{8 + i}:00:00.000-07:00")
                for i in range(2)]
        cls.text, cls.params = run_methods(cls.tmp.name, raws)
        cls.lc = section(cls.text, "Liquid chromatography")
        cls.ms = section(cls.text, "Mass spectrometry")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_lc_names_the_recorded_system_and_method(self):
        self.assertIn("Evosep One LC system (Evosep Biosystems)", self.lc)
        self.assertIn("the 60 samples per day (60 SPD) method (run time 21 min)", self.lc)
        self.assertIn(f"loaded onto Evotips {mm.EVOTIP_TAG}", self.lc)
        self.assertNotIn("[LC system / gradient — confirm]", self.lc)

    def test_lc_paragraph_in_published_order(self):
        # literature survey 2026-09-24: LC coupled to the timsTOF via its source; column and
        # its temperature; the LC method; mobile phases
        order = ["Evosep One LC system", "coupled online to a timsTOF HT mass spectrometer "
                 "(Bruker Daltonics) via a Captive Spray ion source", mm.LC_COLUMN_DEFAULT,
                 "column temperature", "(60 SPD) method", "Mobile phase A"]
        pos = [self.lc.find(x) for x in order]
        self.assertNotIn(-1, pos, dict(zip(order, pos)))
        self.assertEqual(pos, sorted(pos))

    def test_ms_paragraph_in_published_order(self):
        order = ["dia-PASEF mode", "m/z 100–1700", "1/K₀ 0.70–1.30 V·s/cm²",
                 "85 ms each (100% duty cycle)", "one MS1 frame and 11 dia-PASEF frames",
                 "36 isolation windows of 26 Th", "25 Th spacing, 1 Th overlap",
                 "3–4 per TIMS ramp", "m/z 299.5–1200.5", "The cycle time was 1.10 s.",
                 "from 20 eV at 1/K₀ 0.60 V·s/cm² to 65 eV at 1/K₀ 1.60 V·s/cm²",
                 "capillary voltage of 1700 V", "dry gas flow of 3 L/min",
                 "dry temperature of 200 °C", "timsControl", "HyStar 6.3.1.8"]
        pos = [self.ms.find(x) for x in order]
        self.assertNotIn(-1, pos, dict(zip(order, pos)))
        self.assertEqual(pos, sorted(pos))
        self.assertNotIn("99.993933", self.ms)                  # the table keeps the raw value
        self.assertNotIn("≈26.3", self.ms)                      # the old per-window "ramp"

    def test_one_unit_form_and_no_cycle_time_called_duty_cycle(self):
        text = self.lc + self.ms
        self.assertNotIn("1/K₀ =", text)
        for m in __import__("re").finditer(r"1/K₀ [0-9.]+(–[0-9.]+)?", text):
            self.assertTrue(text[m.end():].startswith(" V·s/cm²"), text[m.start():m.end() + 12])
        self.assertNotRegex(text, r"duty cycle (of|was) [0-9.]+ s")

    def test_what_no_file_records_is_tagged(self):
        self.assertIn(f"{mm.LC_COLUMN_DEFAULT} {mm.DEF}", self.lc)
        self.assertIn(f"column temperature of ____ °C {mm.NR_TAG}", self.lc)
        self.assertIn(f"acetonitrile {mm.DEF}", self.lc)
        self.assertIn(f"{mm.EMITTER_DEFAULT} {mm.DEF}", self.ms)

    def test_table_and_params_carry_sources(self):
        tab = section(self.text, "Acquisition parameters (extracted from the raw data)")
        self.assertIn("| LC system | Evosep One (Evosep Biosystems) S/N S00230 |", tab)
        self.assertIn("matches all 36 window energies", tab)
        self.assertIn("| Series consistency | all runs share the acquisition method |", tab)
        self.assertEqual(self.params["column"]["tag"], mm.DEF)
        self.assertEqual(self.params["series"]["differences"], {})
        self.assertTrue(self.params["series"]["acquired_first"].startswith("2026-08-13T18"))
        self.assertNotIn("windows", self.params["representative"])

    def test_mismatched_ramp_is_not_stated(self):
        with tempfile.TemporaryDirectory() as tmp:
            text, _ = run_methods(tmp, [make_timstof_d(os.path.join(tmp, "x.d"), ce_hi=59.0)])
            ms = section(text, "Mass spectrometry")
            self.assertNotIn("ramped linearly", ms)
            self.assertIn("across the isolation windows [they do not match the method's "
                          "recorded ramp — confirm]", ms)

    def test_runs_acquired_differently_are_flagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [make_timstof_d(os.path.join(tmp, "a.d")),
                    make_timstof_d(os.path.join(tmp, "b.d"), method="DIA_other.m")]
            text, params = run_methods(tmp, raws)
            self.assertIn("One paragraph cannot describe every run", text)
            self.assertIn("ms_method", params["series"]["differences"])

    def test_old_d_without_side_files_keeps_the_tagged_lc_sentence(self):
        with tempfile.TemporaryDirectory() as tmp:
            text, _ = run_methods(tmp, [make_timstof_d(os.path.join(tmp, "x.d"),
                                                       side_files=False)])
            lc = section(text, "Liquid chromatography")
            self.assertIn("[LC system / gradient — confirm]", lc)
            self.assertIn("36 isolation windows", section(text, "Mass spectrometry"))


class RenderFromDict(unittest.TestCase):
    """The paragraph writers on a plain .d-metadata dict (what bruker_meta returns)."""

    REP = {"instrument": "timsTOF HT", "mode": "dia-PASEF", "polarity": "positive",
           "mz_low": 99.993933, "mz_high": 1700.0, "im_low": 0.7, "im_high": 1.3,
           "ramp_ms": 85.05, "accumulation_ms": 85.05, "frames_per_cycle": 12, "cycle_s": 1.1,
           "scheme": {"n_windows": 36, "n_ramps": 11, "per_ramp": (3, 4), "width": 26.0,
                      "spacing": 25.0, "overlap": 1.0, "mz_lo": 299.5, "mz_hi": 1200.5,
                      "im_lo": 0.7, "im_hi": 1.3},
           "ce_ramp": {"points": [(0.6, 20.0), (1.6, 65.0)], "advanced": False},
           "ce_check": {"status": "ok"}, "lc_system": "Evosep One",
           "lc_vendor": "Evosep Biosystems", "lc_method": "60 samples per day",
           "lc_run_min": 21.0, "tray_type": "96Evotip", "source_type": "Captive Spray"}
    COL = {"text": "PepSep MAX C18 (Bruker PepSep)", "tag": None}

    @staticmethod
    def v(x, unit="", default=None):
        return f"____ {mm.DEF}" if x is None else f"{x}{unit}"

    def test_expected_sentences(self):
        lc = mm.lc_paragraph(self.REP, self.COL, True, True)
        self.assertIn("analysed on an Evosep One LC system (Evosep Biosystems) coupled online to "
                      "a timsTOF HT mass spectrometer (Bruker Daltonics) via a Captive Spray ion "
                      "source. Peptides were separated on a PepSep MAX C18 (Bruker PepSep), at a "
                      "column temperature of ____ °C [not recorded — confirm], with the 60 "
                      "samples per day (60 SPD) method (run time 21 min).", lc)
        ms = mm.ms_paragraph(self.REP, self.v, coupled=True)
        self.assertTrue(ms.startswith("The mass spectrometer was operated in positive-ion "
                                      "dia-PASEF mode. MS1 and MS2 spectra were recorded over "
                                      "m/z 100–1700. The trapped ion mobility (TIMS) ramp "
                                      "spanned 1/K₀ 0.70–1.30 V·s/cm², with ramp and "
                                      "accumulation times of 85 ms each (100% duty cycle)."), ms)

    def test_an_unknown_value_is_marked_never_filled(self):
        rep = dict(self.REP, im_low=None, im_high=None, polarity=None, cycle_s=None,
                   capillary_v=None, instrument=None)
        ms = mm.ms_paragraph(rep, self.v, coupled=False)
        self.assertIn(f"1/K₀ ____ {mm.DEF}–____ {mm.DEF} V·s/cm²", ms)
        self.assertIn(f"on a ____ {mm.DEF} mass spectrometer", ms)
        self.assertNotIn("positive-ion", ms)                    # polarity not recorded: omitted
        self.assertNotIn("cycle time", ms)
        self.assertNotIn("capillary", ms)
        lc = mm.lc_paragraph(dict(self.REP, lc_method=None, instrument=None), self.COL, True, True)
        self.assertIn(f"with the LC method ____ {mm.NR_TAG}.", lc)
        self.assertIn(f"coupled online to a {mm.NOT_RECORDED} mass spectrometer", lc)


def write_log(path, rows, fmt="csv"):
    keys = ["instrument", "event_type", "event_date", "column_vendor", "column_model",
            "column_serial", "first_run"]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        if fmt == "json":
            json.dump(rows, fh)
        else:
            wr = csv.DictWriter(fh, fieldnames=keys)
            wr.writeheader()
            for r in rows:
                wr.writerow({k: r.get(k, "") for k in keys})
    return path


HT = "timsTOF HT"
CHANGES = [
    {"instrument": HT, "event_type": "column_change", "event_date": "2026-07-31T18:00:00",
     "column_vendor": "Bruker PepSep", "column_model": "PepSep Max C18 10cm x 150um, 1.5um",
     "first_run": "23229"},
    {"instrument": HT, "event_type": "source_clean", "event_date": "2026-08-20T09:00:00"},
    {"instrument": HT, "event_type": "column_change", "event_date": "2026-09-02T05:00:00",
     "column_vendor": "Bruker PepSep", "column_model": "PepSep Max C18 10cm x 150um, 1.5um",
     "column_serial": "0000507372"},
]
FIRST, LAST = "2026-08-13T18:23:03.532-07:00", "2026-08-14T05:51:10.241-07:00"


class Column(unittest.TestCase):
    """column_record(): which column the Methods names, and why."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def rec(self, rows=None, fmt="csv", **kw):
        log = write_log(os.path.join(self.dir, "log." + fmt), rows, fmt) if rows else None
        return mm.column_record(kw.pop("user", None), kw.pop("runs", ()), log,
                                kw.pop("instrument", None), kw.pop("first", FIRST),
                                kw.pop("last", LAST))

    def test_default_is_tagged(self):
        c = self.rec()
        self.assertEqual((c["text"], c["tag"]), (mm.LC_COLUMN_DEFAULT, mm.DEF))

    def test_user_column_wins(self):
        c = self.rec(CHANGES, user="Aurora 25 cm", runs=[{"column_info": "X"}])
        self.assertEqual((c["text"], c["tag"], c["source"]),
                         ("Aurora 25 cm", None, "--lc-column (user-given)"))

    def test_hystar_column_info_beats_the_log(self):
        c = self.rec(CHANGES, runs=[{"column_info": "PepSep 15 cm"}] * 2)
        self.assertEqual(c["text"], "PepSep 15 cm")
        self.assertIn("ColumnInfo", c["source"])

    def test_different_column_infos_are_not_chosen_between(self):
        c = self.rec(runs=[{"column_info": "A"}, {"column_info": "B"}])
        self.assertEqual(c["tag"], mm.DEF)
        self.assertIn("2 different columns", c["warnings"][0])

    def test_log_gives_the_column_in_place_at_the_first_run(self):
        for fmt in ("csv", "json"):
            c = self.rec(CHANGES, fmt=fmt)
            self.assertEqual(c["text"], "PepSep Max C18 10cm x 150um, 1.5um (Bruker PepSep)", fmt)
            self.assertIsNone(c["tag"])
            self.assertIn("column_change of 2026-07-31 18:00:00", c["source"])
            self.assertIn("first run on it 23229", c["source"])
            self.assertEqual(c["warnings"], [])

    def test_serial_and_vendor_are_carried(self):
        c = self.rec(CHANGES, first="2026-09-03T10:00:00-07:00", last="2026-09-03T12:00:00-07:00")
        self.assertTrue(c["text"].endswith("(Bruker PepSep; serial/LOT 0000507372)"), c["text"])
        rows = [dict(CHANGES[0], column_model="Aurora Ultimate 25cm", column_vendor="IonOpticks")]
        self.assertEqual(self.rec(rows)["text"], "Aurora Ultimate 25cm (IonOpticks)")
        rows = [dict(CHANGES[0], column_model="Bruker PepSep MAX 10 cm")]   # vendor not repeated
        self.assertEqual(self.rec(rows)["text"], "Bruker PepSep MAX 10 cm")

    def test_a_change_inside_the_series_falls_back_tagged(self):
        rows = CHANGES + [dict(CHANGES[0], event_date="2026-08-14T00:00:00")]
        c = self.rec(rows)
        self.assertEqual(c["tag"], mm.DEF)
        self.assertIn("during this series", " ".join(c["warnings"]))

    def test_an_offset_aware_event_is_compared_on_the_acquisition_clock(self):
        # 08:00 UTC is 01:00 at the instrument's -07:00: inside the 18:23-05:51 series
        rows = CHANGES + [dict(CHANGES[0], event_date="2026-08-14T08:00:00Z")]
        self.assertIn("during this series", " ".join(self.rec(rows)["warnings"]))

    def test_no_change_before_the_first_run(self):
        c = self.rec(CHANGES[2:])
        self.assertEqual(c["tag"], mm.DEF)
        self.assertIn("no column change before the first run", c["warnings"][0])

    def test_a_change_that_names_no_column(self):
        c = self.rec([dict(CHANGES[0], column_model="")])
        self.assertEqual(c["tag"], mm.DEF)
        self.assertIn("does not name the column", c["warnings"][0])

    def test_several_instruments_need_a_choice(self):
        rows = CHANGES + [dict(CHANGES[0], instrument="Orbitrap Exploris 480",
                               column_model="Aurora")]
        c = self.rec(rows)
        self.assertEqual(c["tag"], mm.DEF)
        self.assertIn("--column-log-instrument", c["warnings"][0])
        self.assertEqual(self.rec(rows, instrument=HT)["text"],
                         "PepSep Max C18 10cm x 150um, 1.5um (Bruker PepSep)")

    def test_log_flows_into_the_methods(self):
        raws = [make_timstof_d(os.path.join(self.dir, "r.d"), acquired=FIRST)]
        log = write_log(os.path.join(self.dir, "stan.csv"), CHANGES)
        text, params = run_methods(self.dir, raws, "--column-log", log)
        lc = section(text, "Liquid chromatography")
        self.assertIn("on a PepSep Max C18 10cm x 150um, 1.5um (Bruker PepSep), at a column "
                      "temperature", lc)
        self.assertNotIn(mm.LC_COLUMN_DEFAULT, lc)
        self.assertIn("column log stan.csv", params["column"]["source"])


if __name__ == "__main__":
    unittest.main()
