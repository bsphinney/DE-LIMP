#!/usr/bin/env python3
"""
A Thermo .raw must yield its acquisition, its instrument and its ACQUIRED precursor m/z bounds.

Before this fix none of the three came back for any .raw. detect_acquisition.py called
ThermoRawFileParser with a command line the parser does not accept -- `metadata -i <raw>`
(there is no `metadata` subcommand: exit 255, "Unexpected extra arguments") and
`query -i <raw>` (query requires -n: exit 255, "specify a valid scan range") -- ignored the
exit codes, and grepped the empty stdout for "dia". Every .raw therefore came back
acquisition unknown / confidence low / instrument null / precursor_mz_range null, so
estimate_params.py searched its 380-980 FALLBACK. Measured on HIVE (srun job 23510571)
against both FRAN pilot Orbitraps, whose methods acquire much wider than that:

    Orbitrap Exploris 480   25 x 35.0 m/z windows   350.0 -1201.0
    Orbitrap Fusion Lumos   19 x 45.7 m/z windows   350.05-1200.95

The fixtures in fixtures/trfp/ are that run's real ThermoRawFileParser 2.0.0.0 output
(`-i=<raw> -m=0 -f=4 -o=<dir>` and `query -i=<raw> -n=<scans> -b=<file>`), cut to one
acquisition cycle -- DIA: an MS1, every window, the next MS1 (scans 1-27 / 1-21); DDA: 30
mid-run scans -- with the peak arrays emptied, the file/method paths replaced (they name
clients), and the instrument serial/slot, creation date and vial/row/sample number replaced
by placeholders. Every attribute the detector reads is verbatim. They carry a trap worth
keeping: the metadata's "MS min MZ"/"MS max MZ" (367.5/1183.5, 372.9/1178.1) are isolation
window CENTRES. Searching those would clip half a window off each end of the range.

fixtures/trfp/fake_trfp.py replays them behind the real parser's command-line grammar and
rejects the two calls that shipped, so these tests fail on the old code for the right
reason. POSIX only: the stand-in is put on PATH as an executable shell shim.
"""
import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
FIXTURES = os.path.join(HERE, "fixtures", "trfp")
FAKE = os.path.join(FIXTURES, "fake_trfp.py")
sys.path.insert(0, SCRIPTS)

import detect_acquisition as da  # noqa: E402

EXPLORIS = "exploris480_dia"   # Orbitrap Exploris 480, DIA, 350.0-1201.0
LUMOS = "lumos_dia"            # Orbitrap Fusion Lumos, DIA, 350.05-1200.95
DDA = "exploris480_dda"        # Orbitrap Exploris 480, DDA (top-N, 1.6 m/z isolation)


def _meta_value(stem, accession):
    with open(os.path.join(FIXTURES, stem + "-metadata.json")) as fh:
        meta = json.load(fh)
    for section in meta.values():
        for term in section:
            if term.get("accession") == accession:
                return term.get("value")
    return None


@unittest.skipUnless(os.name == "posix", "the ThermoRawFileParser stand-in is a POSIX shim")
class _FakeParserCase(unittest.TestCase):
    """A temp dir holding placeholder .raw files and, on PATH, the replaying stand-in."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.bin = os.path.join(self.tmp, "bin")
        os.makedirs(self.bin)
        shim = os.path.join(self.bin, "ThermoRawFileParser")
        with open(shim, "w") as fh:
            fh.write(f'#!/bin/sh\nexec "{sys.executable}" "{FAKE}" "$@"\n')
        os.chmod(shim, os.stat(shim).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
        self.log = os.path.join(self.tmp, "trfp_calls.jsonl")
        self._env = os.environ.copy()
        os.environ["PATH"] = self.bin + os.pathsep + self._env.get("PATH", "")
        os.environ["FAKE_TRFP_LOG"] = self.log
        for k in ("THERMORAWFILEPARSER", "FAKE_TRFP_FAIL", "FAKE_TRFP_GARBAGE",
                  "FAKE_TRFP_NO_FILTER", "FAKE_TRFP_FILTER_ACCESSION", "FAKE_TRFP_SLEEP"):
            os.environ.pop(k, None)

    def tearDown(self):
        os.environ.clear()
        os.environ.update(self._env)
        self._tmp.cleanup()

    def raw(self, stem, subdir=""):
        d = os.path.join(self.tmp, subdir)
        os.makedirs(d, exist_ok=True)
        p = os.path.join(d, stem + ".raw")
        with open(p, "wb") as fh:
            fh.write(b"\x01\xa1F\x00i\x00n\x00n\x00i\x00g\x00a\x00n")  # Finnigan magic, nothing more
        return p

    def hide_parser(self):
        """PATH without the stand-in, or any real ThermoRawFileParser a developer has installed."""
        keep = []
        for d in os.environ.get("PATH", "").split(os.pathsep):
            if d == self.bin or any(os.path.exists(os.path.join(d, n))
                                    for n in ("ThermoRawFileParser", "thermorawfileparser")):
                continue
            keep.append(d)
        os.environ["PATH"] = os.pathsep.join(keep)

    def calls(self):
        if not os.path.exists(self.log):
            return []
        with open(self.log) as fh:
            return [json.loads(line) for line in fh if line.strip()]


class ThermoRawIsReadThroughTheRealCommandLine(_FakeParserCase):
    def test_exploris_dia_yields_acquisition_instrument_and_bounds(self):
        r = da.classify(self.raw(EXPLORIS))
        self.assertEqual(r["acquisition"], "DIA", r["reason"])
        self.assertEqual(r["confidence"], "high", r["reason"])
        self.assertEqual(r["instrument"], "Orbitrap Exploris 480")
        self.assertIsNotNone(r["precursor_mz_range"], r["reason"])
        self.assertAlmostEqual(r["precursor_mz_range"][0], 350.0, places=3)
        self.assertAlmostEqual(r["precursor_mz_range"][1], 1201.0, places=3)

    def test_lumos_bounds_are_window_edges_not_window_centres(self):
        stem = LUMOS
        r = da.classify(self.raw(stem))
        self.assertEqual(r["acquisition"], "DIA", r["reason"])
        self.assertEqual(r["instrument"], "Orbitrap Fusion Lumos")
        lo, hi = r["precursor_mz_range"]
        # exact: the parser's float32 noise (350.0499938964844) must not reach the JSON
        self.assertEqual([lo, hi], [350.05, 1200.95])
        # The metadata JSON offers a tempting shortcut that is wrong by half a window.
        centre_lo = float(_meta_value(stem, "PRIDE:0000476"))
        centre_hi = float(_meta_value(stem, "PRIDE:0000477"))
        self.assertLess(lo, centre_lo - 20, "used the lowest window CENTRE as the bound")
        self.assertGreater(hi, centre_hi + 20, "used the highest window CENTRE as the bound")

    def test_dda_is_detected_with_confidence_and_no_range(self):
        r = da.classify(self.raw(DDA))
        self.assertEqual(r["acquisition"], "DDA", r["reason"])
        self.assertEqual(r["confidence"], "high", r["reason"])
        self.assertEqual(r["instrument"], "Orbitrap Exploris 480")
        # DDA isolation windows are precursor picks, not an acquired range.
        self.assertIsNone(r["precursor_mz_range"])

    def test_dependent_scan_flag_is_used_as_evidence(self):
        """The filter string's `d` token is the instrument saying 'data-dependent'. It is
        stronger evidence than any width heuristic and the reason should say it was read."""
        dia = da.classify(self.raw(EXPLORIS))
        dda = da.classify(self.raw(DDA))
        self.assertIn("data-dependent", dia["reason"])
        self.assertIn("data-dependent", dda["reason"])

    def test_the_calls_use_grammar_the_parser_accepts(self):
        da.classify(self.raw(EXPLORIS))
        calls = self.calls()
        self.assertTrue(calls, "ThermoRawFileParser was never invoked")
        for argv in calls:
            self.assertNotEqual(argv[:1], ["metadata"], "there is no `metadata` subcommand")
        meta = [a for a in calls if a[:1] != ["query"] and any(x.startswith("-m=") for x in a)]
        query = [a for a in calls if a[:1] == ["query"]]
        self.assertEqual(len(meta), 1, calls)
        self.assertEqual(len(query), 1, calls)
        self.assertIn("-m=0", meta[0])
        self.assertTrue(any(x.startswith("-n=") for x in query[0]), "query needs -n=<scans>")
        # -option=value is the only form the README promises for every option
        for argv in meta + query:
            self.assertTrue(any(x.startswith("-i=") for x in argv), argv)

    def test_scans_are_sampled_mid_run_and_cover_several_cycles(self):
        da.classify(self.raw(EXPLORIS))
        query = [a for a in self.calls() if a[:1] == ["query"]][0]
        spec = next(x for x in query if x.startswith("-n="))[3:]
        a, b = (int(v) for v in spec.split("-"))
        first, last = (int(v) for v in _meta_value(EXPLORIS, "PRIDE:0000479").split(":"))
        self.assertTrue(first <= a < b <= last, spec)
        self.assertLessEqual(a, (first + last) // 2)
        self.assertGreaterEqual(b, (first + last) // 2)
        n_ms1 = int(_meta_value(EXPLORIS, "PRIDE:0000481"))
        n_ms2 = int(_meta_value(EXPLORIS, "PRIDE:0000482"))
        cycle = n_ms2 / n_ms1 + 1
        self.assertGreaterEqual(b - a + 1, 3 * cycle, "too few scans to see every window")

    def test_a_waters_raw_folder_is_not_handed_to_the_thermo_reader(self):
        folder = os.path.join(self.tmp, "waters_run.raw")
        os.makedirs(folder)
        r = da.classify(folder)
        self.assertEqual(r["vendor"], "Waters")
        self.assertEqual(self.calls(), [], "ran a Thermo parser on a Waters folder")

    def test_detect_instrument_alone_reads_the_model(self):
        self.assertEqual(da.detect_instrument(self.raw(LUMOS)), "Orbitrap Fusion Lumos")

    def test_filter_string_under_the_pre_1_4_5_accession_is_still_read(self):
        """v1.3.0-v1.4.4 DO write the filter string, as "MS:10000512" (one zero too many;
        Query/ProxiSpectrumReader.cs at every one of those tags, corrected in v1.4.5), and
        bioconda still serves all of them. Reading only MS:1000512 threw the data-dependent
        flag away on those builds -- and with it the check that keeps narrow-window DIA
        from being called DDA."""
        os.environ["FAKE_TRFP_FILTER_ACCESSION"] = "MS:10000512"
        dia = da.classify(self.raw(EXPLORIS))
        self.assertEqual((dia["acquisition"], dia["confidence"]), ("DIA", "high"), dia["reason"])
        self.assertRegex(dia["reason"], r"\b0/\d+ MS2 scans flagged data-dependent")
        self.assertNotIn("not available", dia["reason"])
        dda = da.classify(self.raw(DDA))
        self.assertEqual((dda["acquisition"], dda["confidence"]), ("DDA", "high"), dda["reason"])
        self.assertRegex(dda["reason"], r"\b(\d+)/\1 MS2 scans flagged data-dependent")

    def test_a_parser_that_reports_no_filter_string_still_classifies_from_windows(self):
        """No release omits it, but if one does the width rule alone must still carry DIA
        with its bounds, and the reason must say the flag was missing -- without blaming a
        parser version for it."""
        os.environ["FAKE_TRFP_NO_FILTER"] = "1"
        r = da.classify(self.raw(EXPLORIS))
        self.assertEqual((r["acquisition"], r["confidence"]), ("DIA", "high"), r["reason"])
        self.assertEqual([round(x, 3) for x in r["precursor_mz_range"]], [350.0, 1201.0])
        self.assertIn("not available", r["reason"])
        self.assertNotIn("1.4.5", r["reason"])
        d = da.classify(self.raw(DDA))
        self.assertEqual(d["acquisition"], "DDA", d["reason"])

    def test_reader_field_records_the_parser_version_and_command(self):
        """Which parser build read the file decides what could be read (see the filter
        string accession above), so it belongs in the per-file record."""
        r = da.classify(self.raw(EXPLORIS))
        self.assertIsNotNone(r["reader"])
        self.assertTrue(r["reader"].startswith("ThermoRawFileParser 2.0.0.0 "), r["reader"])
        self.assertIn(os.path.join(self.bin, "ThermoRawFileParser"), r["reader"])
        self.hide_parser()
        self.assertIsNone(da.classify(self.raw(EXPLORIS))["reader"])

    def test_parser_named_by_environment_variable_is_used(self):
        """The framework-dependent release is a DLL (`dotnet ThermoRawFileParser.dll`) and the
        1.4.x builds need mono: neither is a single executable on PATH."""
        self.hide_parser()
        os.environ["THERMORAWFILEPARSER"] = f'"{sys.executable}" "{FAKE}"'
        r = da.classify(self.raw(LUMOS))
        self.assertEqual(r["acquisition"], "DIA", r["reason"])
        self.assertEqual(r["instrument"], "Orbitrap Fusion Lumos")


class DataDependentFlagOutranksWindowShape(unittest.TestCase):
    """Synthetic windows, because no pilot run has these shapes -- the rule still needs pinning."""

    @staticmethod
    def windows(width, centres, cycles, dependent):
        w = {"widths": [], "centres": [], "los": [], "his": [],
             "n_ms2": 0, "n_filter": 0, "n_dependent": 0}
        for _ in range(cycles):
            for c in centres:
                w["widths"].append(width); w["centres"].append(c)
                w["los"].append(c - width / 2); w["his"].append(c + width / 2)
                w["n_ms2"] += 1; w["n_filter"] += 1; w["n_dependent"] += int(dependent)
        return w

    def test_narrow_windows_that_are_not_data_dependent_are_not_called_dda(self):
        # Astral-style 2 m/z DIA: the width rule alone would call this DDA.
        centres = [381.0 + 2.0 * i for i in range(300)]
        kind, conf, why, rng = da.classify_thermo_windows(self.windows(2.0, centres, 3, False))
        self.assertEqual((kind, conf), ("DIA", "medium"), why)
        self.assertEqual(rng, (380.0, 980.0))

    def test_wide_windows_that_are_all_data_dependent_are_not_called_dia(self):
        kind, conf, why, rng = da.classify_thermo_windows(
            self.windows(4.0, [400.0 + 50 * i for i in range(10)], 5, True))
        self.assertEqual((kind, conf), ("DDA", "medium"), why)
        self.assertIsNone(rng)

    def test_mixed_dependent_and_independent_scans_are_never_high(self):
        """A hybrid method (DIA windows plus data-dependent scans) must reach the user."""
        centres = [367.5 + 34.0 * i for i in range(25)]
        w = self.windows(35.0, centres, 4, False)
        w["n_dependent"] = w["n_filter"] // 2
        kind, conf, why, rng = da.classify_thermo_windows(w)
        self.assertEqual((kind, conf), ("DIA", "medium"), why)
        self.assertIn("mixed", why)
        self.assertIsNotNone(rng)

    @staticmethod
    def proxi(width, centres, cycles, filter_accession):
        """TRFP `query` JSON as the parser writes it: an MS1, then one MS2 per window."""
        spectra = []
        for _ in range(cycles):
            spectra.append({"attributes": [
                {"accession": "MS:1000511", "value": "1"},
                {"accession": filter_accession,
                 "value": "FTMS + p NSI Full ms [350.0000-1500.0000]"}]})
            for c in centres:
                spectra.append({"attributes": [
                    {"accession": "MS:1000511", "value": "2"},
                    {"accession": "MS:1000827", "value": str(c)},
                    {"accession": "MS:1000828", "value": str(width / 2)},
                    {"accession": "MS:1000829", "value": str(width / 2)},
                    {"accession": filter_accession,
                     "value": f"FTMS + p NSI Full ms2 {c:.4f}@hcd30.00 [150.0000-2000.0000]"}]})
        return spectra

    def test_narrow_window_dia_is_not_called_dda_on_a_pre_1_4_5_parser(self):
        """End to end from query JSON: 300 x 2 m/z windows, nothing data-dependent, the
        filter string under v1.3.0-v1.4.4's "MS:10000512". Reading only the correct
        accession returned DDA/high with no range and no confirmation -- DIA data routed to
        Sage without asking."""
        centres = [381.0 + 2.0 * i for i in range(300)]
        for acc in ("MS:1000512", "MS:10000512"):
            w = da.thermo_isolation_windows(self.proxi(2.0, centres, 3, acc))
            self.assertEqual((w["n_filter"], w["n_dependent"]), (900, 0), acc)
            kind, conf, why, rng = da.classify_thermo_windows(w)
            self.assertEqual((kind, conf), ("DIA", "medium"), f"{acc}: {why}")
            self.assertEqual(rng, (380.0, 980.0), acc)


class ThermoRawFailuresAreLoud(_FakeParserCase):
    def test_parser_not_found_says_so_and_points_at_the_public_source(self):
        self.hide_parser()
        r = da.classify(self.raw(EXPLORIS))
        self.assertEqual(r["acquisition"], "unknown")
        self.assertEqual(r["confidence"], "low")
        self.assertIsNone(r["precursor_mz_range"])
        self.assertIn("not found", r["reason"])
        self.assertIn("github.com/compomics/ThermoRawFileParser", r["reason"])
        # The downstream consequence must be spelled out, not left for the user to discover.
        self.assertIn("FALLBACK", r["reason"])

    def test_query_failure_reports_exit_code_and_parser_message(self):
        os.environ["FAKE_TRFP_FAIL"] = "query"
        r = da.classify(self.raw(EXPLORIS))
        self.assertEqual(r["acquisition"], "unknown")
        self.assertEqual(r["confidence"], "low")
        self.assertIsNone(r["precursor_mz_range"])
        self.assertIn("exit 1", r["reason"])
        self.assertIn("file is corrupt", r["reason"])
        self.assertIn("FALLBACK", r["reason"])

    def test_unreadable_query_output_is_a_failure_not_a_guess(self):
        os.environ["FAKE_TRFP_GARBAGE"] = "query"
        r = da.classify(self.raw(EXPLORIS))
        self.assertEqual(r["acquisition"], "unknown")
        self.assertIsNone(r["precursor_mz_range"])
        self.assertIn("FALLBACK", r["reason"])

    def test_a_parser_that_does_not_answer_is_a_failure_not_a_hang(self):
        """A raw on a stalled mount must come back unknown, loudly -- not block forever,
        and not look like a result."""
        self.addCleanup(setattr, da, "TRFP_TIMEOUT_S", da.TRFP_TIMEOUT_S)
        da.TRFP_TIMEOUT_S = 1
        os.environ["FAKE_TRFP_SLEEP"] = "4"
        r = da.classify(self.raw(EXPLORIS))
        self.assertEqual((r["acquisition"], r["confidence"]), ("unknown", "low"))
        self.assertIsNone(r["precursor_mz_range"])
        self.assertIn("no answer after 1 s", r["reason"])
        self.assertIn("FALLBACK", r["reason"])
        self.assertTrue(r["warnings"])

    def test_dia_without_readable_window_edges_is_a_warning(self):
        """DIA with no range would silently search the 380-980 FALLBACK. It must warn, and
        a warning is what sets needs_confirmation."""
        with mock.patch.object(da, "classify_thermo_windows",
                               return_value=("DIA", "high", "stub", None)):
            r = da.classify(self.raw(EXPLORIS))
        self.assertEqual(r["acquisition"], "DIA")
        self.assertIsNone(r["precursor_mz_range"])
        self.assertTrue(any("no isolation window edges" in w and "FALLBACK" in w
                            for w in r["warnings"]), r["warnings"])

    def test_metadata_failure_still_classifies_but_flags_the_instrument(self):
        os.environ["FAKE_TRFP_FAIL"] = "metadata"
        r = da.classify(self.raw(EXPLORIS))
        self.assertEqual(r["acquisition"], "DIA", r["reason"])
        self.assertAlmostEqual(r["precursor_mz_range"][0], 350.0, places=3)
        self.assertIsNone(r["instrument"])
        self.assertTrue(any("instrument" in w and "exit 1" in w for w in r["warnings"]),
                        r["warnings"])


class ThermoRawCLI(_FakeParserCase):
    def _run(self, *raws):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                               *raws], capture_output=True, text=True, env=os.environ.copy())

    def test_cli_json_carries_the_measured_range_for_estimate_params(self):
        res = self._run(self.raw(EXPLORIS))
        self.assertEqual(res.returncode, 0, res.stderr)
        payload = json.loads(res.stdout)
        self.assertEqual(payload["overall"], "DIA")
        self.assertEqual(payload["instrument"], "Orbitrap Exploris 480")
        self.assertEqual(payload["precursor_mz_range"], [350.0, 1201.0])
        self.assertEqual(payload["precursor_mz_range_files_without"], [])
        self.assertFalse(payload["needs_confirmation"])

    def test_cli_unions_two_orbitraps(self):
        res = self._run(self.raw(EXPLORIS), self.raw(LUMOS))
        self.assertEqual(res.returncode, 0, res.stderr)
        payload = json.loads(res.stdout)
        self.assertEqual(payload["overall"], "DIA")
        lo, hi = payload["precursor_mz_range"]
        self.assertAlmostEqual(lo, 350.0, places=3)
        self.assertAlmostEqual(hi, 1201.0, places=3)
        self.assertTrue(payload["precursor_mz_range_mixed"])
        self.assertTrue(payload["needs_confirmation"], "two instruments must be confirmed")

    def test_cli_warns_on_stderr_when_the_parser_is_missing(self):
        self.hide_parser()
        res = self._run(self.raw(EXPLORIS))
        self.assertEqual(res.returncode, 0, res.stderr)
        self.assertIn("WARNING", res.stderr)
        self.assertIn("ThermoRawFileParser", res.stderr)
        payload = json.loads(res.stdout)          # stdout stays pure JSON
        self.assertTrue(payload["needs_confirmation"])
        self.assertIsNone(payload["precursor_mz_range"])

    def test_cli_metadata_failure_requires_confirmation(self):
        os.environ["FAKE_TRFP_FAIL"] = "metadata"
        res = self._run(self.raw(EXPLORIS))
        payload = json.loads(res.stdout)
        self.assertTrue(payload["needs_confirmation"],
                        "an unread instrument decides mass accuracy; the user must see it")
        self.assertIn("WARNING", res.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
