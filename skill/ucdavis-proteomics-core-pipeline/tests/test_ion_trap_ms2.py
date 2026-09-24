#!/usr/bin/env python3
"""
Orbitrap MS1 / ion-trap MS2 (a Fusion Lumos OT/IT method): MS2 is read in the ion trap, so it has
no Orbitrap resolution at all. detect_acquisition.py reports it per file as ms2_analyzer "ITMS"
and lists the files under ms2_ion_trap; the orchestrator passes --ms2-analyzer ITMS.

Before this:
  * Sage got a 10 ppm FRAGMENT window for every Orbitrap class. An ion-trap spectrum's fragment
    error is Da-scale, so almost nothing matched. Now: Sage's own documented low-res MS/MS
    setting, fragment_tol {"da": [-0.4, 0.4]} (sage-docs.vercel.app/docs/configuration/tolerance),
    with deisotope off and a larger bucket_size, as the same docs recommend.
  * DIA-NN, given only the MS1 resolution, fell to orbitrap_generic and was labelled as if the
    resolution were unknown, and asked for an MS2 resolution that does not exist. Now: class
    orbitrap_iontrap, both mass-accuracy flags omitted (either one fixes both levels), labelled
    as an ion-trap MS2, and never asked for an MS2 resolution.

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
ITMS = ("--ms1-resolution", "120000", "--ms2-analyzer", "ITMS", "--resolution-source", "detected")


class Scripts(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def run_script(self, name, *args):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, name), *args],
                              capture_output=True, text=True, timeout=120)

    def estimate(self, *args, engine="sage", acquisition="DDA", name="p"):
        """-> (cfg text, rationale sidecar, stderr)"""
        out = os.path.join(self.d, name)
        p = self.run_script("estimate_params.py", "--engine", engine, "--acquisition",
                            acquisition, "--instrument", LUMOS, "--out", out, *args)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out) as fh, open(out + ".rationale.json") as side:
            return fh.read(), json.load(side), p.stderr

    def resolve(self, *args, engine="diann", acquisition="DIA", name="wf"):
        """-> (manifest, stderr)"""
        dest = os.path.join(self.d, name)
        p = self.run_script("resolve_defaults.py", "--acquisition", acquisition, "--engine",
                            engine, "--instrument", LUMOS, "--dest", dest, *args)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(dest, "workflow.manifest.json")) as fh:
            return json.load(fh), p.stderr


class SageIonTrapTests(Scripts):
    def test_an_ion_trap_ms2_gets_sages_documented_da_window(self):
        text, side, _ = self.estimate(*ITMS)
        cfg = json.loads(text)
        self.assertEqual(cfg["fragment_tol"], {"da": [-0.4, 0.4]})
        self.assertIs(cfg["deisotope"], False)
        self.assertEqual(cfg["database"]["bucket_size"], 32768)
        # the precursor is still the Orbitrap MS1's
        self.assertEqual(cfg["precursor_tol"], {"ppm": [-10.0, 10.0]})
        r = side["rationale"]
        self.assertEqual(r["fragment_tol_da"]["value"], 0.4)
        self.assertIn("Sage docs", r["fragment_tol_da"]["source"])
        self.assertIn("ion trap", r["fragment_tol_da"]["source"])
        self.assertNotIn("fragment_tol_ppm", r, "a ppm fragment tag beside a Da window")
        self.assertEqual(side["resolution"]["ms2_analyzer"], "ITMS")
        self.assertEqual(side["resolution"]["source"], "detected")

    def test_the_value_is_cited_where_it_is_defined(self):
        """The lead's rule: from Sage's own docs, not memory -- so the source stays next to it."""
        self.assertEqual(ep.SAGE_ITMS_FRAGMENT_DA, 0.4)
        with open(os.path.join(SCRIPTS, "estimate_params.py")) as fh:
            src = fh.read()
        self.assertIn("https://sage-docs.vercel.app/docs/configuration/tolerance", src)
        self.assertIn('{ "da": [-0.4, 0.4] }', src)
        self.assertIn("pub enum Tolerance { Ppm(f32, f32), Da(f32, f32) }", src)

    def test_ftms_is_unchanged(self):
        """--ms2-analyzer FTMS must write exactly what the same run wrote without the flag."""
        res = ("--ms1-resolution", "120000", "--ms2-resolution", "30000")
        with_ftms, s1, _ = self.estimate(*res, "--ms2-analyzer", "FTMS", name="ftms")
        without, s2, _ = self.estimate(*res, name="none")
        self.assertEqual(with_ftms, without)
        self.assertEqual(json.loads(without)["fragment_tol"], {"ppm": [-10.0, 10.0]})
        self.assertIs(json.loads(without)["deisotope"], True)
        self.assertEqual(s1["rationale"], s2["rationale"])


class DiannIonTrapTests(Scripts):
    def test_nothing_is_pinned_and_the_label_says_why(self):
        text, side, err = self.estimate(*ITMS, engine="diann", acquisition="DIA", name="d.cfg")
        self.assertEqual(side["instrument_class"], "orbitrap_iontrap")
        lines = text.splitlines()
        self.assertFalse([ln for ln in lines if ln.startswith("--mass-acc")],
                         "a lone flag would fix both levels, the other at 20 ppm")
        self.assertEqual(side["mass_accuracy_plan"], "auto")
        plan = side["rationale"]["mass_accuracy_plan"]["source"]
        self.assertIn("ion trap", plan)
        self.assertIn("left to DIA-NN calibration", plan)
        blob = json.dumps(side)
        self.assertNotIn("resolution unknown", blob)
        self.assertNotIn("instrument not identified", blob)
        self.assertIn("read in the ion trap", side["class_label"])
        # there is nothing to ask for
        self.assertIs(side["needs_confirmation"], False)
        self.assertIsNone(side["ask_user"])
        self.assertNotIn("NEEDS CONFIRMATION", err)

    def test_no_ms2_resolution_is_asked_for_even_without_ms1(self):
        _, side, _ = self.estimate("--ms2-analyzer", "ITMS", engine="diann", acquisition="DIA",
                                   name="n.cfg")
        self.assertIs(side["needs_confirmation"], False)

    def test_radiant_asks_for_ms1_only(self):
        man, err = self.resolve("--ms2-analyzer", "ITMS", engine="radiant")
        self.assertIs(man["needs_confirmation"], True)
        ask = man["ask_user"]
        self.assertIn("MS1 resolution", ask)
        self.assertIn("--ms1-resolution", ask)
        self.assertIn("--ms2-analyzer ITMS", ask)
        self.assertNotIn("--ms2-resolution", ask)
        self.assertNotIn("MS2 resolution", ask)
        tol = man["search"]["preset_provenance"]["tolerances"]
        self.assertIn("ion trap", tol)

    def test_both_at_once_is_refused(self):
        for script, extra in (("estimate_params.py", ("--engine", "diann", "--out",
                                                      os.path.join(self.d, "x.cfg"))),
                              ("resolve_defaults.py", ("--dest", os.path.join(self.d, "x")))):
            with self.subTest(script=script):
                p = self.run_script(script, "--acquisition", "DIA", "--instrument", LUMOS,
                                    "--ms2-resolution", "30000", "--ms2-analyzer", "ITMS", *extra)
                self.assertNotEqual(p.returncode, 0)
                self.assertIn("contradict", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.d, "x.cfg")))


class ResolveDefaultsTests(Scripts):
    def test_both_flag_spellings_carry_the_analyzer(self):
        a, _ = self.resolve("--ms1-res", "120000", "--ms2-analyzer", "ITMS",
                            "--resolution-source", "detected", name="short")
        b, _ = self.resolve("--ms1-resolution", "120000", "--ms2-analyzer", "itms",
                            "--resolution-source", "detected", name="long")
        self.assertEqual(a["resolution"], b["resolution"])
        self.assertEqual(a["resolution"],
                         {"ms1": 120000, "ms2": None, "ms2_analyzer": "ITMS", "source": "detected",
                          "source_label": ep.RESOLUTION_SOURCES["detected"]})
        self.assertEqual(a["search"], b["search"])
        self.assertEqual(a["instrument_class"], "orbitrap_iontrap")
        self.assertEqual((a["search"]["ms1_ppm"], a["search"]["ms2_ppm"]), (7, None))
        self.assertIn("left to DIA-NN calibration", a["search"]["ppm_source"])
        self.assertIs(a["needs_confirmation"], False)

    def test_sage_route_names_the_da_window(self):
        man, err = self.resolve(*ITMS, engine="sage", acquisition="DDA", name="sage")
        self.assertIn("0.4 Da", man["search"]["ppm_source"])
        self.assertIs(man["needs_confirmation"], False)
        self.assertNotIn("NEEDS CONFIRMATION", err)


class MixedAnalyzerTests(Scripts):
    """A file with MS2 from both analyzers: one fragment window must fit the ion-trap spectra."""

    def test_mixed_gets_the_ion_trap_window_and_says_mixed(self):
        for spelling in ("mixed", "MIXED"):
            with self.subTest(spelling=spelling):
                text, side, _ = self.estimate("--ms1-resolution", "120000", "--ms2-analyzer",
                                              spelling, name=f"m_{spelling}")
                self.assertEqual(json.loads(text)["fragment_tol"], {"da": [-0.4, 0.4]})
                self.assertEqual(side["resolution"]["ms2_analyzer"], "mixed")
                self.assertIn("mixed", side["class_label"])
                self.assertIn("mixed", side["rationale"]["fragment_tol_da"]["source"])

    def test_mixed_leaves_dia_nn_to_calibrate(self):
        text, side, _ = self.estimate("--ms1-resolution", "120000", "--ms2-analyzer", "mixed",
                                      engine="diann", acquisition="DIA", name="m.cfg")
        self.assertFalse([ln for ln in text.splitlines() if ln.startswith("--mass-acc")])
        self.assertEqual(side["mass_accuracy_plan"], "auto")
        self.assertIs(side["needs_confirmation"], False)


class MakePresetsGuardTests(Scripts):
    def test_make_presets_refuses_an_ion_trap_ms2_with_a_resolution(self):
        out = os.path.join(self.d, "r.radiantConfig")
        p = self.run_script("make_presets.py", "--engine", "radiant", "--instrument", LUMOS,
                            "--ms1-res", "120000", "--ms2-res", "30000", "--ms2-analyzer", "ITMS",
                            "--out", out)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("contradict", p.stderr)
        self.assertFalse(os.path.exists(out))


class RunSearchCrossCheckTests(Scripts):
    """SKILL.md passes --ms2-analyzer to resolve_defaults.py (step 4, the manifest) AND to
    estimate_params.py (step 6b, the Sage cfg); only the second changes the search. run_search.py
    checks the two agree before it converts, writes or submits anything."""

    def setUp(self):
        super().setUp()
        self.cfg = {}
        for kind, extra in (("itms", ("--ms2-analyzer", "ITMS")),
                            ("ftms", ("--ms2-resolution", "30000"))):
            self.estimate("--ms1-resolution", "120000", *extra, name=f"{kind}.json")
            self.cfg[kind] = os.path.join(self.d, f"{kind}.json")
        self.tools = os.path.join(self.d, "tools.json")
        with open(self.tools, "w") as fh:
            json.dump({"sage": "/opt/sage/sage", "versions": {"sage": "v0.14.7"}}, fh)
        self.fasta = os.path.join(self.d, "db.fasta")
        with open(self.fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEK\n")
        self.mzml = os.path.join(self.d, "a.mzML")
        open(self.mzml, "w").close()

    def search(self, manifest_analyzer, cfg_kind):
        bundle = os.path.join(self.d, "workflow.manifest.json")
        res = (None if manifest_analyzer == "absent" else
               {"ms1": 120000, "ms2": None, "ms2_analyzer": manifest_analyzer,
                "source": "detected"})
        with open(bundle, "w") as fh:
            json.dump({"acquisition": "DDA", "engine": {"name": "sage", "version": "0.14.7"},
                       "resolution": res}, fh)
        job = os.path.join(self.d, f"job_{manifest_analyzer}_{cfg_kind}.sh")
        out = os.path.join(self.d, f"out_{manifest_analyzer}_{cfg_kind}")
        # --sbatch only WRITES a script; PATH without sbatch keeps the login-node guard away
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                            "--tools", self.tools, "--bundle", bundle,
                            "--params", self.cfg[cfg_kind], "--fasta", self.fasta, "--out", out,
                            "--files", self.mzml, "--engine", "sage", "--sbatch", job],
                           capture_output=True, text=True, timeout=120,
                           env=dict(os.environ, PATH="/usr/bin:/bin"))
        return p, os.path.exists(job), os.path.exists(out)

    def test_an_ion_trap_manifest_with_a_ppm_cfg_is_refused_before_anything_is_written(self):
        for analyzer in ("ITMS", "mixed"):
            with self.subTest(analyzer=analyzer):
                p, job, out = self.search(analyzer, "ftms")
                self.assertNotEqual(p.returncode, 0)
                self.assertIn("REFUSED", p.stderr)
                self.assertIn("re-run estimate_params.py", p.stderr)
                self.assertIn(f"--ms2-analyzer {analyzer}", p.stderr)
                self.assertFalse(job or out, "something was generated before the refusal")

    def test_an_ftms_manifest_with_a_da_cfg_is_refused_too(self):
        p, job, out = self.search("FTMS", "itms")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("REFUSED", p.stderr)
        self.assertIn("--ms2-analyzer ITMS", p.stderr)
        self.assertFalse(job or out)

    def test_matching_pairs_run(self):
        for analyzer, kind in (("ITMS", "itms"), ("mixed", "itms"), ("FTMS", "ftms"),
                               ("absent", "ftms")):
            with self.subTest(analyzer=analyzer, cfg=kind):
                p, job, _ = self.search(analyzer, kind)
                self.assertEqual(p.returncode, 0, p.stderr)
                self.assertTrue(job)
                self.assertNotIn("REFUSED", p.stderr)
                self.assertNotIn("WARNING: " + self.cfg[kind], p.stderr)

    def test_a_da_cfg_with_no_recorded_analyzer_runs_with_a_warning(self):
        """Nothing to check it against: say so, do not block."""
        p, job, _ = self.search("absent", "itms")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertTrue(job)
        self.assertIn("WARNING", p.stderr)
        self.assertIn("records no ms2_analyzer", p.stderr)


class Ms1OnlyLabelTests(unittest.TestCase):
    def test_a_known_ms1_is_not_called_unknown(self):
        """The other half of the bug: MS1-only (no analyzer) fell to 'resolution unknown -- pass
        --ms1-resolution/--ms2-resolution' although MS1 had been given."""
        cls, ms1, ms2, label, src = ep.classify_instrument(LUMOS, 120000, None, "user")
        self.assertEqual((cls, ms1, ms2), ("orbitrap_generic", None, None))
        self.assertIn("MS1 120,000 resolution supplied by the user", label)
        self.assertIn("MS2 resolution unknown", label)
        self.assertNotIn("--ms1-resolution", label)

    def test_an_unknown_resolution_orbitrap_is_not_called_unidentified(self):
        _, rationale = ep.build_diann("DIA", "orbitrap_generic", None, None, "Orbitrap", "x",
                                      [], {})
        self.assertNotIn("instrument not identified",
                         rationale["mass_accuracy_plan"]["source"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
