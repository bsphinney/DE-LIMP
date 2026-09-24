#!/usr/bin/env python3
"""
The Orbitrap MS1/MS2 resolution, read from the .raw scan trailer.

estimate_params.py pins DIA-NN's documented Orbitrap tolerances from the resolution; without it
an Orbitrap is `orbitrap_generic`, DIA-NN calibrates per run and the 5-step parallel chain
declines (gabrig 2026-09-23, 15 Fusion Lumos .raw). ThermoRawFileParser never outputs the
resolution, but every scan's trailer carries it ('Orbitrap Resolution:' on a Fusion Lumos,
'FT Resolution:' on an Exploris 480; HIVE jobs 23989081, 23989170: Lumos 60000/15000, Exploris
480 DIA 120000/15000). scripts/thermo_resolution.py reads it through pythonnet and Thermo's
RawFileReader DLLs, in a subprocess of detect_acquisition.py.

No .NET and no DLLs here: fixtures/fake_rawfilereader/ stands in for pythonnet, `clr` and the
ThermoFisher namespaces, and the REAL thermo_resolution.py runs on it under `python -S` (so a
pythonnet installed on the machine cannot take over). Each fake .raw's scans are described by
the `<raw>.scans.json` beside it. The .raw itself is read by the ThermoRawFileParser stand-in of
test_thermo_raw_detection.py, so acquisition, instrument and range come from real captured
output. POSIX only.
"""
import json
import os
import subprocess
import sys
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
FAKES = os.path.join(HERE, "fixtures", "fake_rawfilereader")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import detect_acquisition as da           # noqa: E402
import thermo_resolution as tr            # noqa: E402
import estimate_params                    # noqa: E402
# Module imports, not `from ... import`: a TestCase class bound in this namespace would be
# collected and run a second time here.
import test_thermo_raw_detection as trfp  # noqa: E402
import test_trfp_dotnet_discovery as dn   # noqa: E402

MS1_LUMOS = "FTMS + p NSI Full ms [350.0000-1500.0000]"
MS2_LUMOS = "FTMS + p NSI Full ms2 368.0000@hcd35.00 [250.0000-746.0000]"
MS2_ITMS = "ITMS + c NSI r d Full ms2 572.3181@cid35.00 [150.0000-2000.0000]"


def cycle(ms1, ms2, key="Orbitrap Resolution:", ms2_filter=MS2_LUMOS, n_ms2=4):
    """One DIA-like cycle: an MS1 then n_ms2 MS2 scans, each trailer with its resolution."""
    scans = [{"filter": MS1_LUMOS, "trailer": {"Ion Injection Time (ms):": "50",
                                               key: str(ms1)}}]
    for _ in range(n_ms2):
        trailer = {"Ion Injection Time (ms):": "22"}
        if not ms2_filter.startswith("ITMS"):
            trailer[key] = str(ms2)
        scans.append({"filter": ms2_filter, "trailer": trailer})
    return scans


class _ResolutionCase(trfp._FakeParserCase):
    """The TRFP stand-in on PATH with the two RawFileReader DLLs (empty files) beside it, a
    .NET root with Microsoft.NETCore.App 8 (a directory tree), and the reader run by `python -S`
    with the fakes on PYTHONPATH."""

    def setUp(self):
        super().setUp()
        for name in da.RAWFILEREADER_DLLS:
            open(os.path.join(self.bin, name), "w").close()
        os.environ["PATH"] = os.pathsep.join(
            d for d in os.environ["PATH"].split(os.pathsep)
            if not os.path.exists(os.path.join(d, "dotnet")))
        for k in ("DOTNET_ROOT", "DOTNET_CORE_SDK_ROOT", "DOTNET_ROOT_X64", "DOTNET_ROOT_ARM64",
                  "FAKE_PYTHONNET_LOAD_FAIL", "FAKE_RFR_SLEEP_ON", "FAKE_RFR_CRASH_ON"):
            os.environ.pop(k, None)
        os.environ["PROTEOMICS_DOTNET_SYSTEM_ROOTS"] = ""
        self.root = dn.fake_dotnet_root(os.path.join(self.tmp, "dotnet8"), (dn.NETCORE,))
        os.environ["PROTEOMICS_DOTNET_DIR"] = self.root
        os.environ["THERMO_RESOLUTION_PYTHON"] = dn._exe(
            os.path.join(self.tmp, "py_S"), f'#!/bin/sh\nexec "{sys.executable}" -S "$@"\n')
        os.environ["PYTHONPATH"] = FAKES
        self.pn_log = os.path.join(self.tmp, "pythonnet.log")
        self.rfr_log = os.path.join(self.tmp, "rfr.log")
        os.environ["FAKE_PYTHONNET_LOG"] = self.pn_log
        os.environ["FAKE_RFR_LOG"] = self.rfr_log

    def lumos(self, name, ms1=60000, ms2=15000, spec=None, **kw):
        """A Fusion Lumos DIA .raw (TRFP fixture) whose scan trailers say ms1/ms2."""
        raw = self.raw(trfp.LUMOS, subdir=name)
        with open(raw + ".scans.json", "w") as fh:
            json.dump(spec if spec is not None else
                      {"n": 5000, "cycle": cycle(ms1, ms2, **kw)}, fh)
        return raw

    def detect(self, *paths):
        res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                              *paths], capture_output=True, text=True, env=os.environ.copy())
        self.assertEqual(res.returncode, 0, res.stderr)
        return json.loads(res.stdout), res.stderr

    def lines(self, path):
        return dn._lines(path)


# --------------------------------------------------------------------------------------------
class ResolutionIsReadFromTheTrailer(_ResolutionCase):
    def test_a_cohort_that_agrees_gets_one_value_at_the_top(self):
        raws = [self.lumos(f"r{i}") for i in range(3)]
        out, _ = self.detect(*raws)
        self.assertEqual((out["ms1_resolution"], out["ms2_resolution"]), (60000, 15000))
        self.assertEqual(out["resolution_mixed"], [])
        self.assertIsNone(out["orbitrap_resolution_unknown"])
        self.assertFalse(out["needs_confirmation"])
        self.assertIsNone(out["ms2_ion_trap"], "pure Orbitrap MS2 is not ion trap")
        for f in out["files"]:
            self.assertEqual((f["ms1_resolution"], f["ms2_resolution"]), (60000, 15000))
            self.assertEqual(f["ms2_analyzer"], "FTMS")
            self.assertIn("Orbitrap Resolution:", f["resolution_note"])
            self.assertIn("RawFileReader 8.0.6.0", f["resolution_note"])
            self.assertEqual(f["warnings"], [])
        # one reader process for the whole cohort, told exactly which .NET to use
        loads = [json.loads(ln) for ln in self.lines(self.pn_log)]
        self.assertEqual(len(loads), 1, loads)
        self.assertEqual(loads[0]["runtime"], "coreclr")
        self.assertEqual(loads[0]["params"]["dotnet_root"], self.root)
        self.assertEqual(loads[0]["DOTNET_ROOT"], self.root)
        self.assertTrue(loads[0]["params"]["runtime_config"].endswith(".runtimeconfig.json"))
        self.assertFalse(os.path.exists(loads[0]["params"]["runtime_config"]),
                         "the temp runtimeconfig was left behind")
        self.assertEqual(sorted(self.lines(self.rfr_log)), sorted(raws))

    def test_the_reader_shows_which_scans_the_values_came_from(self):
        """The evidence that the values are the method's own MS1 and MS2: the full filter
        strings of the scans read, and every scan type in the walk (an extra SIM or survey scan
        type would show up there)."""
        raw = self.lumos("r0")
        p = subprocess.run([os.environ["THERMO_RESOLUTION_PYTHON"], da.RES_READER, "--dll-dir",
                            self.bin, raw], capture_output=True, text=True,
                           env=dict(os.environ, DOTNET_ROOT=self.root))
        rec = json.loads(p.stdout)
        self.assertEqual((rec["ms1_filter"], rec["ms2_filter"]), (MS1_LUMOS, MS2_LUMOS))
        self.assertEqual(rec["scan_types"], {"FTMS + p NSI Full ms": 80,
                                             "FTMS + p NSI Full ms2 @hcd35.00": 320})
        self.assertEqual(rec["scans_walked"], "2300-2699")          # mid-run, 400 scans

    def test_the_exploris_key_is_read_too(self):
        raw = self.lumos("ex", ms1=120000, ms2=15000, key="FT Resolution:")
        out, _ = self.detect(raw)
        self.assertEqual((out["ms1_resolution"], out["ms2_resolution"]), (120000, 15000))
        self.assertIn("FT Resolution:", out["files"][0]["resolution_note"])

    def test_the_values_feed_estimate_params(self):
        """What step 6b does with them: a Lumos at 60k/15k is no longer orbitrap_generic."""
        out, _ = self.detect(self.lumos("r0"))
        cls, ms1, ms2, _label, _src = estimate_params.classify_instrument(
            out["instrument"], out["ms1_resolution"], out["ms2_resolution"])
        self.assertNotEqual(cls, "orbitrap_generic")
        self.assertEqual(ms1, 10)                       # 60k MS1: DIA-NN's table
        self.assertEqual(estimate_params.classify_instrument(out["instrument"])[0],
                         "orbitrap_generic", "without the values it would have been")


class AMixedCohortGetsNoValue(_ResolutionCase):
    def test_two_methods_are_listed_not_merged(self):
        a = [self.lumos(f"hi{i}", ms1=120000, ms2=30000) for i in range(2)]
        b = self.lumos("lo", ms1=60000, ms2=15000)
        out, _ = self.detect(*a, b)
        self.assertEqual((out["ms1_resolution"], out["ms2_resolution"]), (None, None))
        groups = {(g["ms1_resolution"], g["ms2_resolution"]): sorted(g["files"])
                  for g in out["resolution_mixed"]}
        self.assertEqual(groups, {(120000, 30000): sorted(a), (60000, 15000): [b]})
        self.assertTrue(out["needs_confirmation"], "a mixed cohort went out unconfirmed")
        self.assertIsNone(out["orbitrap_resolution_unknown"], "every file WAS read")

    def test_a_level_that_agrees_keeps_its_value(self):
        a = self.lumos("a", ms1=60000, ms2=15000)
        b = self.lumos("b", ms1=60000, ms2=30000)
        out, _ = self.detect(a, b)
        self.assertEqual((out["ms1_resolution"], out["ms2_resolution"]), (60000, None))
        self.assertEqual(len(out["resolution_mixed"]), 2)

    def test_two_ms2_resolutions_in_one_file_give_no_value(self):
        spec = {"n": 5000, "cycle": cycle(60000, 15000) + [
            {"filter": MS2_LUMOS, "trailer": {"Orbitrap Resolution:": "30000"}}]}
        out, _ = self.detect(self.lumos("two", spec=spec))
        f = out["files"][0]
        self.assertEqual((f["ms1_resolution"], f["ms2_resolution"]), (60000, None))
        self.assertIn("several resolutions", f["resolution_note"])
        self.assertEqual(out["orbitrap_resolution_unknown"]["files"], [f["file"]])


class WhatCannotBeReadIsAskedFor(_ResolutionCase):
    def test_only_the_unread_file_is_asked_about(self):
        good = [self.lumos(f"g{i}") for i in range(2)]
        bad = self.lumos("bad", spec={"error": "The file is corrupt"})
        out, _ = self.detect(*good, bad)
        unknown = out["orbitrap_resolution_unknown"]
        self.assertEqual(unknown["files"], [bad])
        self.assertEqual(unknown["read_in_other_files"], [[60000, 15000]])
        self.assertTrue(any("The file is corrupt" in r for r in unknown["reasons"]),
                        unknown["reasons"])
        self.assertIn("ASK the user", unknown["ask"])
        # not every Orbitrap file has a value, so none goes to the top ...
        self.assertEqual((out["ms1_resolution"], out["ms2_resolution"]), (None, None))
        # ... but nothing disagrees, and the reads themselves were clean
        self.assertEqual(out["resolution_mixed"], [])
        self.assertFalse(out["needs_confirmation"])

    def test_an_ion_trap_ms2_is_not_asked_for(self):
        """No MS2 Orbitrap resolution exists to ask for: only MS1 applies, and the MS2
        tolerance cannot come from DIA-NN's Orbitrap table."""
        raw = self.lumos("itms", ms1=120000, ms2=None, ms2_filter=MS2_ITMS)
        out, _ = self.detect(raw)
        f = out["files"][0]
        self.assertEqual((f["ms1_resolution"], f["ms2_resolution"], f["ms2_analyzer"]),
                         (120000, None, "ITMS"))
        self.assertIn("ion trap (ITMS", f["resolution_note"])
        self.assertIsNone(out["orbitrap_resolution_unknown"], "asked for an MS2 that is not there")
        itms = out["ms2_ion_trap"]
        self.assertEqual((itms["files"], itms["ms1_resolution"]), ([raw], 120000))
        for words in ("only the MS1 Orbitrap resolution applies", "do not pass --ms2-resolution",
                      "Orbitrap resolution table", "--ms2-analyzer ITMS"):
            self.assertIn(words, itms["note"])
        self.assertEqual((itms["mixed_files"], itms["mixed_note"]), ([], None))
        self.assertEqual((out["ms1_resolution"], out["ms2_resolution"]), (120000, None))
        self.assertFalse(out["needs_confirmation"], "pure ion trap is a known method")

    def test_orbitrap_and_ion_trap_ms2_in_one_run_is_mixed_not_ftms(self):
        """A Tribrid decision-tree method: HCD-OT for high charge, CID-IT for 2+. Reported as
        FTMS it would have Sage search the ion-trap spectra at +/-10 ppm and lose them."""
        spec = {"n": 5000, "cycle": cycle(120000, 30000, n_ms2=3) + [
            {"filter": MS2_ITMS, "trailer": {"Ion Injection Time (ms):": "35"}}] * 2}
        raw = self.lumos("tree", spec=spec)
        out, _ = self.detect(raw)
        f = out["files"][0]
        self.assertEqual(f["ms2_analyzer"], "mixed")
        # 400 walked scans of a 6-scan cycle (1 MS1, 3 OT MS2, 2 IT MS2)
        self.assertIn("MS2 is MIXED: 201 Orbitrap (FTMS) and 133 ion-trap (ITMS)",
                      f["resolution_note"])
        self.assertIn("x 'ITMS + c NSI r d Full ms2 @cid35.00'", f["resolution_note"])
        itms = out["ms2_ion_trap"]
        self.assertEqual((itms["files"], itms["mixed_files"]), ([raw], [raw]))
        self.assertEqual(itms["ms2_analyzer"], {raw: "mixed"})
        for words in ("CONFIRM the method", "--ms2-analyzer mixed", "no --ms2-resolution",
                      "0.4 Da"):
            self.assertIn(words, itms["mixed_note"])
        self.assertTrue(out["needs_confirmation"], "a mixed-analyzer run went out unconfirmed")
        self.assertIsNone(out["ms2_resolution"], "an Orbitrap MS2 value for a mixed cohort")
        self.assertEqual(out["ms1_resolution"], 120000)
        self.assertIsNone(out["orbitrap_resolution_unknown"], "nothing to ask: all was read")

    def test_an_ion_trap_ms2_with_an_unread_ms1_asks_for_ms1_only(self):
        spec = {"n": 500, "cycle": [{"filter": MS1_LUMOS, "trailer": {}},
                                    {"filter": MS2_ITMS, "trailer": {}}]}
        raw = self.lumos("itms_nores", spec=spec)
        out, _ = self.detect(raw)
        unknown = out["orbitrap_resolution_unknown"]
        self.assertEqual(unknown["levels"], {raw: ["MS1"]})

    def test_ion_trap_and_orbitrap_ms2_runs_are_a_mixed_cohort(self):
        a = self.lumos("ot", ms1=120000, ms2=30000)
        b = self.lumos("it", ms1=120000, ms2=None, ms2_filter=MS2_ITMS)
        out, _ = self.detect(a, b)
        self.assertEqual(out["ms1_resolution"], 120000)
        self.assertIsNone(out["ms2_resolution"])
        analyzers = {g["ms2_analyzer"]: g["files"] for g in out["resolution_mixed"]}
        self.assertEqual(analyzers, {"FTMS": [a], "ITMS": [b]})
        self.assertTrue(out["needs_confirmation"])

    def test_no_pythonnet_leaves_detection_alone(self):
        """The REAL reader under `python -S` with no fakes: pythonnet cannot be imported."""
        os.environ.pop("PYTHONPATH")
        raws = [self.lumos(f"r{i}") for i in range(2)]
        out, err = self.detect(*raws)
        self.assertEqual(out["overall"], "DIA")
        self.assertEqual(out["precursor_mz_range"], [350.05, 1200.95])
        self.assertFalse(out["needs_confirmation"])
        for f in out["files"]:
            self.assertEqual((f["acquisition"], f["confidence"]), ("DIA", "high"))
            self.assertEqual(f["warnings"], [])
            self.assertIsNone(f["ms1_resolution"])
            self.assertIn("pythonnet is not installed", f["resolution_note"])
            self.assertIn("setup.sh", f["resolution_note"])
        self.assertEqual(sorted(out["orbitrap_resolution_unknown"]["files"]), sorted(raws))
        self.assertNotIn("Traceback", err)

    def test_no_interpreter_with_pythonnet_says_the_same(self):
        os.environ["THERMO_RESOLUTION_PYTHON"] = ""
        out, _ = self.detect(self.lumos("r0"))
        self.assertIn("pythonnet is not installed", out["files"][0]["resolution_note"])
        self.assertEqual(self.lines(self.pn_log), [], "nothing should have been started")

    def test_no_dlls_anywhere_is_a_note_naming_the_override(self):
        for name in da.RAWFILEREADER_DLLS:
            os.remove(os.path.join(self.bin, name))
        out, _ = self.detect(self.lumos("r0"))
        note = out["files"][0]["resolution_note"]
        self.assertIn("THERMO_RAWFILEREADER_DIR", note)
        self.assertEqual(out["files"][0]["acquisition"], "DIA")

    def test_the_dll_folder_can_be_named(self):
        elsewhere = os.path.join(self.tmp, "rfr")
        os.makedirs(elsewhere)
        for name in da.RAWFILEREADER_DLLS:
            os.rename(os.path.join(self.bin, name), os.path.join(elsewhere, name))
        os.environ["THERMO_RAWFILEREADER_DIR"] = elsewhere
        out, _ = self.detect(self.lumos("r0"))
        self.assertEqual(out["ms1_resolution"], 60000)
        self.assertIn(elsewhere, out["files"][0]["resolution_note"])

    def test_no_dotnet_for_pythonnet_names_the_fix(self):
        os.environ["PROTEOMICS_DOTNET_DIR"] = os.path.join(self.tmp, "none")
        out, _ = self.detect(self.lumos("r0"))
        note = out["files"][0]["resolution_note"]
        self.assertIn("Microsoft.NETCore.App 8+", note)
        self.assertIn("ensure_dotnet8.sh", note)
        self.assertIn("self-contained", note)

    def test_a_runtime_that_will_not_load_is_a_note(self):
        os.environ["FAKE_PYTHONNET_LOAD_FAIL"] = "1"
        out, _ = self.detect(self.lumos("r0"))
        self.assertIn("pythonnet could not load .NET", out["files"][0]["resolution_note"])
        self.assertEqual(out["files"][0]["acquisition"], "DIA")


class TheReaderCannotHangOrCrashDetection(_ResolutionCase):
    def test_a_reader_that_does_not_answer_is_cut_off_and_keeps_what_it_said(self):
        first, stuck = self.lumos("a"), self.lumos("b")
        os.environ["FAKE_RFR_SLEEP_ON"] = os.path.basename(stuck)   # both are lumos_dia.raw
        os.environ["FAKE_RFR_SLEEP"] = "20"
        launch = da.trfp_launch(da.find_trfp())
        with mock.patch.object(da, "RES_TIMEOUT_BASE_S", 3), \
                mock.patch.object(da, "RES_TIMEOUT_EACH_S", 0):
            got = da.read_resolutions([first, stuck], launch)
        # same basename, so the first file sleeps too: nothing answered in time
        for p in (first, stuck):
            self.assertIsNone(got[p]["ms1_resolution"])
            self.assertIn("no answer within 3 s", got[p]["note"])

    def test_answers_before_a_hang_are_kept(self):
        first = self.lumos("a")
        stuck = self.raw(trfp.EXPLORIS, subdir="b")
        with open(stuck + ".scans.json", "w") as fh:
            json.dump({"n": 50, "cycle": cycle(120000, 15000, key="FT Resolution:")}, fh)
        os.environ["FAKE_RFR_SLEEP_ON"] = os.path.basename(stuck)
        os.environ["FAKE_RFR_SLEEP"] = "20"
        launch = da.trfp_launch(da.find_trfp())
        with mock.patch.object(da, "RES_TIMEOUT_BASE_S", 4), \
                mock.patch.object(da, "RES_TIMEOUT_EACH_S", 0):
            got = da.read_resolutions([first, stuck], launch)
        self.assertEqual(got[first]["ms1_resolution"], 60000)
        self.assertIsNone(got[stuck]["ms1_resolution"])
        self.assertIn("no answer within 4 s", got[stuck]["note"])

    def test_a_native_crash_costs_only_the_files_it_had_not_answered(self):
        first = self.lumos("a")
        crash = self.raw(trfp.EXPLORIS, subdir="b")
        os.environ["FAKE_RFR_CRASH_ON"] = os.path.basename(crash)
        out, _ = self.detect(first, crash)
        by = {f["file"]: f for f in out["files"]}
        self.assertEqual(by[first]["ms1_resolution"], 60000)
        self.assertIsNone(by[crash]["ms1_resolution"])
        self.assertIn("exited 134", by[crash]["resolution_note"])
        self.assertEqual(by[crash]["acquisition"], "DIA", "detection went down with the reader")


class AHangCostsOneBatch(_ResolutionCase):
    """RES_BATCH files per reader process, each with its own timeout: a hang inside .NET costs
    its batch, and two hung batches in a row (a stalled mount) stop the rest."""

    def _read(self, raws, hang, batch=3):
        os.environ["FAKE_RFR_SLEEP_ON"] = ",".join(hang)
        os.environ["FAKE_RFR_SLEEP"] = "30"
        launch = da.trfp_launch(da.find_trfp())
        with mock.patch.object(da, "RES_BATCH", batch), \
                mock.patch.object(da, "RES_TIMEOUT_BASE_S", 3), \
                mock.patch.object(da, "RES_TIMEOUT_EACH_S", 0):
            return da.read_resolutions(raws, launch)

    def test_a_hang_on_one_file_costs_only_its_batch(self):
        raws = [self.lumos(f"f{i}") for i in range(7)]          # batches: 0-2, 3-5, 6
        got = self._read(raws, ["f1"])
        self.assertEqual(got[raws[0]]["ms1_resolution"], 60000, "answered before the hang")
        for p in raws[1:3]:
            self.assertIsNone(got[p]["ms1_resolution"])
            self.assertIn("no answer within 3 s", got[p]["note"])
        for p in raws[3:]:
            self.assertEqual(got[p]["ms1_resolution"], 60000, p)

    def test_two_hung_batches_in_a_row_stop_the_rest(self):
        raws = [self.lumos(f"f{i}") for i in range(9)]          # batches: 0-2, 3-5, 6-8
        got = self._read(raws, ["f0", "f3"])
        for p in raws[:6]:
            self.assertIsNone(got[p]["ms1_resolution"], p)
        for p in raws[6:]:
            self.assertIn("not attempted: the reader hung on 2 batches in a row", got[p]["note"])
        self.assertNotIn(raws[6], self.lines(self.rfr_log), "a batch after the stop was read")

    def test_the_budgets_bound_a_200_file_cohort(self):
        """Worst case: two hung batches, then the stop -- minutes, not the ~5 h of one timeout
        for everything."""
        worst = da.RES_HUNG_BATCHES_STOP * (da.RES_TIMEOUT_BASE_S
                                            + da.RES_TIMEOUT_EACH_S * da.RES_BATCH)
        self.assertLessEqual(worst, 15 * 60)
        self.assertLess(tr.WALK_SECONDS, da.RES_TIMEOUT_EACH_S)


class DetectionNeverFailsOverTheResolution(_ResolutionCase):
    def locked_root(self, frameworks):
        """A .NET root whose host/fxr cannot be listed (restored by each test before tearDown)."""
        bad = dn.fake_dotnet_root(os.path.join(self.tmp, "locked"), frameworks)
        os.chmod(os.path.join(bad, "host", "fxr"), 0)
        return bad

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0, "root reads anything")
    def test_an_unreadable_dotnet_root_is_a_reason_not_a_traceback(self):
        """dotnet_root_lacks listed <root>/host/fxr unguarded: a PermissionError there aborted
        detection with a traceback and no JSON."""
        bad = self.locked_root((dn.NETCORE,))
        os.environ["DOTNET_ROOT"] = bad
        os.environ["PROTEOMICS_DOTNET_DIR"] = os.path.join(self.tmp, "none")
        raw = self.lumos("r0")
        try:
            out, err = self.detect(raw)
        finally:
            os.chmod(os.path.join(bad, "host", "fxr"), 0o755)
        self.assertNotIn("Traceback", err)
        f = out["files"][0]
        self.assertEqual((f["acquisition"], f["instrument"]), ("DIA", "Orbitrap Fusion Lumos"))
        self.assertIn("unreadable", f["resolution_note"])
        self.assertIn(bad, f["resolution_note"])
        self.assertEqual(out["orbitrap_resolution_unknown"]["files"], [f["file"]])

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0, "root reads anything")
    def test_an_unreadable_root_does_not_stop_the_parser_either(self):
        """The same unguarded listdir sat on the ThermoRawFileParser path (trfp_launch)."""
        bad = self.locked_root((dn.NETCORE, dn.ASPNET))
        try:
            lacks = da.dotnet_root_lacks(bad, [(dn.NETCORE, (8, 0, 0), False)])
        finally:
            os.chmod(os.path.join(bad, "host", "fxr"), 0o755)
        self.assertTrue(lacks[0].startswith("unreadable (PermissionError"), lacks)

    def test_anything_the_reader_raises_becomes_the_reason(self):
        import contextlib
        import io
        raw = self.lumos("r0")
        buf = io.StringIO()
        with mock.patch.object(da, "read_resolutions", side_effect=RuntimeError("boom")), \
                contextlib.redirect_stdout(buf), contextlib.redirect_stderr(io.StringIO()):
            da.main([raw])
        out = json.loads(buf.getvalue())
        self.assertEqual(out["files"][0]["acquisition"], "DIA")
        self.assertIn("resolution reader failed (RuntimeError: boom)",
                      out["files"][0]["resolution_note"])
        self.assertEqual(out["orbitrap_resolution_unknown"]["files"], [raw])
        self.assertIn("RuntimeError: boom", out["orbitrap_resolution_unknown"]["reasons"][0])


class OnlyOrbitrapsOfUnknownResolutionAreRead(unittest.TestCase):
    def test_an_astral_is_not_read_and_says_why(self):
        recs = [{"file": "a.raw", "vendor": "Thermo", "instrument": "Orbitrap Astral",
                 "ms1_resolution": None, "ms2_resolution": None, "resolution_note": None},
                {"file": "b.d", "vendor": "Bruker", "instrument": "timsTOF HT",
                 "ms1_resolution": None, "ms2_resolution": None, "resolution_note": None}]
        with mock.patch.object(da, "read_resolutions") as rr:
            da.add_resolutions(recs, {"cmd": ["x"], "dotnet_root": None})
        rr.assert_not_called()
        self.assertIn("orbitrap_astral", recs[0]["resolution_note"])
        self.assertIn("MS1 4 / MS2 10 ppm", recs[0]["resolution_note"])
        self.assertIsNone(recs[1]["resolution_note"])
        top = da.resolution_summary(recs)
        self.assertEqual((top["ms1_resolution"], top["resolution_mixed"],
                          top["orbitrap_resolution_unknown"]), (None, [], None))

    def test_the_reader_is_given_only_the_orbitraps(self):
        recs = [{"file": f"{n}.raw", "vendor": "Thermo", "instrument": inst,
                 "ms1_resolution": None, "ms2_resolution": None, "resolution_note": None}
                for n, inst in (("l", "Orbitrap Fusion Lumos"), ("a", "Orbitrap Astral"),
                                ("x", None))]
        seen = {}

        def fake(paths, launch):
            seen["paths"] = paths
            return {p: {"file": p, "ms1_resolution": 60000, "ms2_resolution": 15000,
                        "ms1_key": "Orbitrap Resolution:", "ms2_key": "Orbitrap Resolution:",
                        "ms1_scan": 1, "ms2_scan": 2, "reader": "r", "note": None}
                    for p in paths}
        with mock.patch.object(da, "read_resolutions", side_effect=fake):
            da.add_resolutions(recs, {"cmd": ["x"], "dotnet_root": None})
        self.assertEqual(seen["paths"], ["l.raw"])
        self.assertIn("instrument is unknown", recs[2]["resolution_note"])


class ScanFilterGrammar(unittest.TestCase):
    def test_levels_and_analyzers(self):
        self.assertEqual(tr.scan_level(MS1_LUMOS), ("FTMS", 1, True))
        self.assertEqual(tr.scan_level(MS2_LUMOS), ("FTMS", 2, True))
        self.assertEqual(tr.scan_level(MS2_ITMS), ("ITMS", 2, True))
        self.assertEqual(tr.scan_level("FTMS + p NSI Full msx ms2 500.00@hcd30.00 [200-2000]"),
                         ("FTMS", 2, True))
        self.assertEqual(tr.scan_level("FTMS + p NSI SIM ms [500.0000-520.0000]"),
                         ("FTMS", 1, False))
        self.assertEqual(tr.scan_level(""), ("", None, False))

    def test_trailer_keys(self):
        self.assertEqual(tr.trailer_resolution(["Orbitrap Resolution:"], ["60000"]),
                         ("Orbitrap Resolution:", 60000))
        self.assertEqual(tr.trailer_resolution([" FT Resolution: "], ["120000.0 "]),
                         ("FT Resolution:", 120000))
        self.assertEqual(tr.trailer_resolution(["Mass Resolution:", "FT Resolution:"],
                                               ["0.5", "0"]), (None, None))


class ReaderReadinessForSetup(_ResolutionCase):
    def _check(self):
        res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                              "--check-reader"], capture_output=True, text=True,
                             env=os.environ.copy())
        return json.loads(res.stdout)

    def test_ready_when_it_loads(self):
        rr = self._check()["resolution_reader"]
        self.assertIs(rr["ready"], True, rr)
        self.assertEqual((rr["dll_dir"], rr["dotnet_root"]), (os.path.realpath(self.bin),
                                                              self.root))
        self.assertIn("RawFileReader 8.0.6.0", rr["reader"])
        self.assertEqual(self.lines(self.rfr_log), [], "readiness must not open a .raw")

    def test_not_ready_without_pythonnet_and_the_parser_is_still_ready(self):
        os.environ.pop("PYTHONPATH")
        st = self._check()
        self.assertIs(st["ready"], True, "the parser's readiness does not depend on it")
        self.assertIs(st["resolution_reader"]["ready"], False)
        self.assertIn("pythonnet is not installed", st["resolution_reader"]["note"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
