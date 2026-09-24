#!/usr/bin/env python3
"""
Where an Orbitrap's resolution came from, and what happens when nobody knows it.

From a real Fusion Lumos session (gabrig, skill 2.5.0):

  #11 resolve_defaults.py took --ms1-res/--ms2-res and estimate_params.py
      --ms1-resolution/--ms2-resolution, so the same numbers needed two spellings. And values
      the USER typed in were labelled "MS1 60,000 / MS2 15,000 resolution read from the data"
      in workflow.manifest.json and <cfg>.rationale.json -- a value presented as something it
      is not (DE-LIMP rule #2). Both spellings now work in both scripts, and the label says
      where the numbers came from: the .raw scan trailer (--resolution-source detected, what
      detect_acquisition.py reads), the user (the default), the mzML, or a saved configuration.
  #12 A Fusion Lumos with NO resolution exited 0 with "resolution unknown" in the manifest and
      nothing asked the user. It still exits 0, but sets needs_confirmation + ask_user and says
      so on stderr -- for the engines the resolution changes (DIA-NN, Radiant), not for routes
      that ignore it (FragPipe keeps vendor tolerances; Sage goes by instrument class).

Stdlib only; runs the scripts as the orchestrator does, in temp dirs, with no network.
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

import estimate_params as ep  # noqa: E402

LUMOS = "Orbitrap Fusion Lumos"
USER = ep.RESOLUTION_SOURCES["user"]
DETECTED = ep.RESOLUTION_SOURCES["detected"]
FALSE_CLAIM = "read from the data"


class Scripts(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def resolve(self, *args, name="wf", instrument=LUMOS, acquisition="DIA"):
        """resolve_defaults.py -> (stdout JSON, manifest, stderr)"""
        dest = os.path.join(self.d, name)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--acquisition", acquisition, "--instrument", instrument,
                            "--dest", dest, *args], capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(dest, "workflow.manifest.json")) as fh:
            return json.loads(p.stdout), json.load(fh), p.stderr

    def estimate(self, *args, name="params.cfg", instrument=LUMOS, engine="diann",
                 acquisition="DIA"):
        """estimate_params.py -> (stdout JSON, <cfg>.rationale.json, cfg text, stderr)"""
        out = os.path.join(self.d, name)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", engine, "--acquisition", acquisition,
                            "--instrument", instrument, "--out", out, *args],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out + ".rationale.json") as fh:
            side = json.load(fh)
        with open(out) as fh:
            return json.loads(p.stdout), side, fh.read(), p.stderr


class BothSpellingsTests(Scripts):
    SHORT = ("--ms1-res", "120000", "--ms2-res", "30000")
    LONG = ("--ms1-resolution", "120000", "--ms2-resolution", "30000")

    def test_resolve_defaults_takes_both_and_they_mean_the_same(self):
        _, short, _ = self.resolve(*self.SHORT, name="short")
        _, long_, _ = self.resolve(*self.LONG, name="long")
        self.assertEqual(short["search"], long_["search"])
        self.assertEqual(short["instrument_label"], long_["instrument_label"])
        self.assertEqual(short["resolution"], long_["resolution"])
        self.assertEqual((short["search"]["ms1_ppm"], short["search"]["ms2_ppm"]), (7, 15))

    def test_estimate_params_takes_both_and_they_mean_the_same(self):
        s_out, s_side, s_cfg, _ = self.estimate(*self.SHORT, name="short.cfg")
        l_out, l_side, l_cfg, _ = self.estimate(*self.LONG, name="long.cfg")
        self.assertEqual(s_cfg, l_cfg)
        self.assertIn("--mass-acc 15", s_cfg.splitlines())
        for k in ("class_label", "resolution", "rationale", "mass_accuracy_plan"):
            self.assertEqual(s_side[k], l_side[k], k)

    def test_the_short_spelling_is_a_real_alias_not_an_abbreviation(self):
        """argparse happens to accept --ms1-res as a PREFIX of --ms1-resolution; one more
        --ms1-res... option would make that ambiguous and break every existing command."""
        import argparse
        ap = argparse.ArgumentParser(allow_abbrev=False)
        ep.add_resolution_args(ap)
        a = ap.parse_args(["--ms1-res", "120000", "--ms2-res", "30000"])
        b = ap.parse_args(["--ms1-resolution", "120000", "--ms2-resolution", "30000"])
        self.assertEqual(vars(a), vars(b))
        self.assertEqual((a.ms1_resolution, a.ms2_resolution), (120000.0, 30000.0))

    def test_make_presets_takes_both_too(self):
        """resolve_defaults.py hands the Radiant/FragPipe config to make_presets.py by the short
        spelling; the long one must not be the one that breaks there instead."""
        for flags in (self.SHORT, self.LONG):
            out = os.path.join(self.d, f"r{len(flags[0])}.radiantConfig")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_presets.py"),
                                "--engine", "radiant", "--instrument", LUMOS, *flags,
                                "--out", out], capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("120,000", json.loads(p.stdout)["instrument_label"])


class SourceLabelTests(Scripts):
    RES = ("--ms1-res", "60000", "--ms2-res", "15000")

    def test_numbers_on_the_command_line_are_the_users_by_default(self):
        """gabrig's exact case: typed-in 60k/15k came out 'read from the data'."""
        out, man, _ = self.resolve(*self.RES)
        self.assertEqual(man["resolution"],
                         {"ms1": 60000, "ms2": 15000, "ms2_analyzer": None, "source": "user",
                          "source_label": USER})
        self.assertIn(USER, man["instrument_label"])
        self.assertIn(USER, man["name"])
        self.assertEqual(out["resolution"]["source"], "user")
        self.assertNotIn(FALSE_CLAIM, json.dumps(man))

        e_out, side, _, _ = self.estimate(*self.RES)
        self.assertIn(USER, side["class_label"])
        self.assertEqual(side["resolution"]["source"], "user")
        # the per-flag rationale carries the label too -- that is what the methods text quotes
        self.assertIn(USER, side["rationale"]["--mass-acc-ms1"]["source"])
        self.assertNotIn(FALSE_CLAIM, json.dumps(side))
        self.assertEqual(e_out, side)

    def test_detected_says_the_scan_trailer(self):
        _, man, _ = self.resolve(*self.RES, "--resolution-source", "detected")
        self.assertEqual(man["resolution"]["source"], "detected")
        self.assertIn(DETECTED, man["instrument_label"])
        self.assertNotIn(USER, json.dumps(man))

        _, side, _, _ = self.estimate(*self.RES, "--resolution-source", "detected")
        self.assertIn(DETECTED, side["class_label"])
        self.assertIn(DETECTED, side["rationale"]["--mass-acc-ms1"]["source"])

    def test_a_saved_configuration_is_not_called_data(self):
        _, man, _ = self.resolve(*self.RES, "--resolution-source", "cfg")
        self.assertIn("not read from the data", man["instrument_label"])

    def test_radiant_preset_provenance_carries_the_same_source(self):
        """make_presets.py writes the label into every narrowed extraction width's reason."""
        _, man, _ = self.resolve("--engine", "radiant", "--ms1-res", "120000",
                                 "--ms2-res", "30000", "--resolution-source", "detected")
        prov = man["search"]["preset_provenance"]
        self.assertIn(DETECTED, prov["instrument_label"])
        self.assertIn(DETECTED, json.dumps(prov["keys_changed"]))
        self.assertNotIn(FALSE_CLAIM, json.dumps(prov))

    def test_from_mzml_is_labelled_as_the_mzml(self):
        mzml = os.path.join(self.d, "run.mzML")
        with open(mzml, "w") as fh:
            for res in (120000, 30000, 30000):
                fh.write(f'<cvParam cvRef="MS" accession="MS:1000800" '
                         f'name="mass resolving power" value="{res}"/>\n')
        _, side, _, err = self.estimate("--from-mzml", mzml)
        self.assertEqual(side["resolution"]["source"], "mzml")
        self.assertIn(ep.RESOLUTION_SOURCES["mzml"], side["class_label"])
        self.assertIn("read resolution from run.mzML", err)

    def test_no_resolution_records_no_source(self):
        _, man, _ = self.resolve(instrument="timsTOF HT")
        self.assertIsNone(man["resolution"])

    def test_an_unstated_source_is_never_guessed_in_the_function(self):
        """A direct caller that passes resolutions without a source gets 'not recorded', not a
        claim either way."""
        label = ep.classify_instrument(LUMOS, 120000, 30000)[3]
        self.assertIn(ep.RESOLUTION_SOURCE_UNRECORDED, label)
        for phrase in ep.RESOLUTION_SOURCES.values():
            self.assertNotIn(phrase, label)

    def test_the_phrases_live_in_one_file(self):
        """DE-LIMP rule 3: one definition. A second copy of a label is how 'read from the data'
        survived in one place after being fixed in another."""
        def text(f):
            with open(os.path.join(SCRIPTS, f)) as fh:
                return fh.read()
        for phrase in (DETECTED, USER):
            holders = [f for f in sorted(os.listdir(SCRIPTS))
                       if f.endswith(".py") and phrase in text(f)]
            self.assertEqual(holders, ["estimate_params.py"], phrase)


class NeedsConfirmationTests(Scripts):
    def test_an_orbitrap_with_no_resolution_asks(self):
        """Exit 0 as before (resolve() asserts it), but flagged at the top of the manifest and
        stdout, and printed to stderr."""
        out, man, err = self.resolve()
        self.assertIs(man["needs_confirmation"], True)
        self.assertIs(out["needs_confirmation"], True)
        for words in ("Ask the user", "MS1 and MS2", "detect_acquisition.py", "--ms1-resolution"):
            self.assertIn(words, man["ask_user"])
        self.assertEqual(out["ask_user"], man["ask_user"])
        self.assertIn("NEEDS CONFIRMATION", err)

        e_out, side, _, e_err = self.estimate()
        self.assertIs(side["needs_confirmation"], True)
        self.assertEqual(side["ask_user"], man["ask_user"])
        self.assertIn("NEEDS CONFIRMATION", e_err)

    def test_only_the_missing_level_is_asked_for(self):
        _, man, _ = self.resolve("--ms1-res", "120000")
        self.assertIn("MS2 resolution", man["ask_user"])
        self.assertNotIn("MS1 and", man["ask_user"])

    def test_radiant_asks_too(self):
        _, man, _ = self.resolve("--engine", "radiant")
        self.assertIs(man["needs_confirmation"], True)

    def test_astral_and_timstof_do_not(self):
        for instrument in ("Orbitrap Astral", "timsTOF HT"):
            with self.subTest(instrument=instrument):
                out, man, err = self.resolve(instrument=instrument, name=instrument[:5])
                self.assertIs(man["needs_confirmation"], False)
                self.assertIsNone(man["ask_user"])
                self.assertIs(out["needs_confirmation"], False)
                self.assertNotIn("NEEDS CONFIRMATION", err)
                e_out, side, _, e_err = self.estimate(instrument=instrument,
                                                      name=instrument[:5] + ".cfg")
                self.assertIs(side["needs_confirmation"], False)
                self.assertNotIn("NEEDS CONFIRMATION", e_err)

    def test_a_known_resolution_does_not(self):
        _, man, err = self.resolve("--ms1-res", "120000", "--ms2-res", "15000")
        self.assertIs(man["needs_confirmation"], False)
        self.assertNotIn("NEEDS CONFIRMATION", err)

    def test_routes_that_ignore_the_resolution_do_not_ask(self):
        """FragPipe keeps its vendor tolerances and Sage goes by class: asking would change
        nothing about the search."""
        _, man, err = self.resolve("--engine", "fragpipe")
        self.assertIs(man["needs_confirmation"], False)
        _, man, _ = self.resolve(acquisition="DDA", name="dda")
        self.assertEqual(man["engine"]["name"], "sage")
        self.assertIs(man["needs_confirmation"], False)
        _, side, _, err = self.estimate(engine="sage", acquisition="DDA", name="sage.json")
        self.assertIs(side["needs_confirmation"], False)
        self.assertNotIn("NEEDS CONFIRMATION", err)


class SageOrbitrapTests(Scripts):
    def test_giving_sage_the_resolution_does_not_make_the_instrument_unknown(self):
        """orbitrap_measured used to fall through sage_ppm() to 20/20 'instrument not
        identified' -- for the instrument that got 10/10 when its resolution was unknown."""
        _, without, _, _ = self.estimate(engine="sage", acquisition="DDA", name="a.json")
        _, with_res, _, _ = self.estimate("--ms1-res", "120000", "--ms2-res", "30000",
                                          engine="sage", acquisition="DDA", name="b.json")
        for side in (without, with_res):
            self.assertEqual(side["rationale"]["precursor_tol_ppm"]["value"], 10)
            self.assertNotIn("not identified", side["rationale"]["precursor_tol_ppm"]["source"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
