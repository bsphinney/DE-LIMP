#!/usr/bin/env python3
"""
The reproducibility bundle describes the environment the DE RAN in (staff report, 2026-09-28).

A DE run in a separate R 4.6 env (limpa 1.4.0, which the default env's limpa 1.2.5 could not
replace) would have been bundled with setup.json's env: provenance.py took env_prefix / rscript
from setup.json, so conda-explicit.txt and r-sessionInfo.txt would have named the wrong R and
limpa. run_de.R now records `runtime` in de_provenance.json (R home, Rscript, library paths, conda
env, container), and provenance.py builds environment/ from it, copies the DE's own
sessionInfo.txt, says when setup.json names another env, and says when a DE record has no
runtime at all instead of passing setup.json's off as the DE's.

The two "envs" are fake prefixes with fake Rscripts that print what they are. Hermetic.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")

DE_SESSIONINFO = ("R version 4.6.0 (2026-04-24) \n\nR home   {home}\nlimpa    1.4.0\n"
                  "limma    3.68.0\n")


def fake_env(root, name, says):
    """<root>/<name>: a conda-looking prefix whose bin/Rscript prints `says`."""
    prefix = os.path.join(root, name)
    os.makedirs(os.path.join(prefix, "conda-meta"))
    os.makedirs(os.path.join(prefix, "lib", "R", "library"))
    os.makedirs(os.path.join(prefix, "bin"))
    rs = os.path.join(prefix, "bin", "Rscript")
    with open(rs, "w") as fh:
        fh.write(f"#!/bin/sh\necho '{says}'\n")
    os.chmod(rs, 0o755)
    return prefix, rs


class BundleDescribesTheDe(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.setup_prefix, self.setup_rs = fake_env(self.tmp, "default_env",
                                                    "R version 4.5.3 limpa 1.2.5")
        self.de_prefix, self.de_rs = fake_env(self.tmp, "r46_env", "R version 4.6.0 limpa 1.4.0")
        self.setup_json = os.path.join(self.tmp, "setup.json")
        with open(self.setup_json, "w") as fh:
            json.dump({"rscript": self.setup_rs, "env_prefix": self.setup_prefix}, fh)

    def de_dir(self, runtime=True, session_info=True):
        de = os.path.join(self.tmp, "tables")
        os.makedirs(de, exist_ok=True)
        rec = {"method": "dpc", "R_version": "4.6.0",
               "packages": {"limpa": "1.4.0", "limma": "3.68.0", "arrow": None}}
        if runtime:
            rec["runtime"] = {"r_version": "R version 4.6.0 (2026-04-24)",
                              "r_home": os.path.join(self.de_prefix, "lib", "R"),
                              "rscript": self.de_rs,
                              "lib_paths": [os.path.join(self.de_prefix, "lib", "R", "library")],
                              "conda_prefix": None, "container": None, "container_name": None,
                              "docker": False, "session_info": "sessionInfo.txt"}
        with open(os.path.join(de, "de_provenance.json"), "w") as fh:
            json.dump(rec, fh)
        if session_info:
            with open(os.path.join(de, "sessionInfo.txt"), "w") as fh:
                fh.write(DE_SESSIONINFO.format(home=os.path.join(self.de_prefix, "lib", "R")))
        return de

    def bundle(self, de):
        out = os.path.join(self.tmp, "repro")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                            "--outdir", out, "--de-dir", de, "--setup-json", self.setup_json,
                            "--timestamp", "2026-10-01T00:00:00Z"],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)

        def read(*parts):
            with open(os.path.join(out, *parts), encoding="utf-8") as fh:
                return fh.read()
        return (json.loads(read("run_manifest.json")), read("environment", "r-sessionInfo.txt"),
                read("MANIFEST.txt"), read("REPRODUCE.md"),
                json.loads(read("environment", "versions.txt")))

    def test_the_bundle_is_the_des_environment_not_setup_jsons(self):
        man, si, manifest, md, versions = self.bundle(self.de_dir())
        env = man["environment"]
        self.assertEqual(env["rscript"], self.de_rs)
        self.assertEqual(env["env_prefix"], self.de_prefix, "derived from R home: <prefix>/lib/R")
        self.assertTrue(env["source"].startswith("de_provenance.json runtime"))
        self.assertEqual(env["de_runtime"]["r_version"], "R version 4.6.0 (2026-04-24)")
        # the DE's own sessionInfo.txt, not setup.json's R
        self.assertIn("limpa    1.4.0", si)
        self.assertNotIn("1.2.5", si)
        self.assertEqual(versions["de_r"]["packages"]["limpa"], "1.4.0")
        # and that setup.json names another env is said, not hidden
        self.assertEqual(env["setup_json_differs"]["rscript"],
                         {"setup_json": self.setup_rs, "de": self.de_rs})
        self.assertIn("[NOTE]    setup.json names a different environment", manifest)
        self.assertIn("## The environment the DE ran in", md)
        self.assertIn("limpa 1.4.0", md)
        self.assertIn(f"conda env: `{self.de_prefix}`", md)

    def test_without_its_sessioninfo_the_des_rscript_is_asked(self):
        _, si, manifest, _, _ = self.bundle(self.de_dir(session_info=False))
        self.assertIn("limpa 1.4.0", si)
        self.assertIn(f"from the DE's Rscript ({self.de_rs})", manifest)

    def test_an_older_de_record_says_the_environment_is_setup_jsons(self):
        man, si, manifest, md, _ = self.bundle(self.de_dir(runtime=False))
        env = man["environment"]
        self.assertEqual(env["rscript"], self.setup_rs)
        self.assertIn("NOT necessarily the environment the DE ran in", env["source"])
        self.assertIsNone(env["de_runtime"])
        self.assertIn("[SKIPPED] DE runtime -- de_provenance.json records none", manifest)
        self.assertIn("NOT RECORDED -- this DE record has no `runtime`", md)
        # the DE's own sessionInfo.txt is still the DE's: it is used
        self.assertIn("limpa    1.4.0", si)

    def test_a_container_is_named_for_the_re_run(self):
        de = self.de_dir()
        with open(os.path.join(de, "de_provenance.json")) as fh:
            rec = json.load(fh)
        rec["runtime"]["container"] = "/containers/delimp-r.sif"
        with open(os.path.join(de, "de_provenance.json"), "w") as fh:
            json.dump(rec, fh)
        _, _, _, md, versions = self.bundle(de)
        self.assertIn("inside the container `/containers/delimp-r.sif`", md)
        self.assertEqual(versions["de_r"]["container"], "/containers/delimp-r.sif")
        with open(os.path.join(self.tmp, "repro", "reproduce.sh")) as fh:
            self.assertIn("The DE ran inside the container /containers/delimp-r.sif", fh.read())



class AMissingValueIsNotADifference(unittest.TestCase):
    """2026-10: a SILAC DE's record had no `rscript`, and MANIFEST.txt said "[NOTE] setup.json
    names a different environment than the one the DE ran in (rscript: setup.json <env>/bin/
    Rscript, DE None)" -- the DE had run in that same env, and SKILL.md relays the NOTE."""
    setUp, de_dir, bundle = (BundleDescribesTheDe.setUp, BundleDescribesTheDe.de_dir,
                             BundleDescribesTheDe.bundle)

    def same_env_record(self, rscript):
        de = self.de_dir()
        with open(os.path.join(de, "de_provenance.json")) as fh:
            rec = json.load(fh)
        rec["runtime"]["r_home"] = os.path.join(self.setup_prefix, "lib", "R")
        if rscript is None:
            rec["runtime"].pop("rscript")
        else:
            rec["runtime"]["rscript"] = rscript
        with open(os.path.join(de, "de_provenance.json"), "w") as fh:
            json.dump(rec, fh)
        return de

    def test_an_unrecorded_rscript_is_not_compared(self):
        man, _, manifest, _, _ = self.bundle(self.same_env_record(None))
        self.assertIsNone(man["environment"]["setup_json_differs"])
        self.assertNotIn("different environment", manifest)
        self.assertIn("the DE's record has no rscript, so setup.json's is not compared", manifest)

    def test_one_envs_two_rscripts_are_the_same_environment(self):
        # R.home("bin") is <env>/lib/R/bin; setup.json names <env>/bin/Rscript
        lib_rs = os.path.join(self.setup_prefix, "lib", "R", "bin", "Rscript")
        os.makedirs(os.path.dirname(lib_rs))
        shutil.copy(self.setup_rs, lib_rs)
        man, _, manifest, _, _ = self.bundle(self.same_env_record(lib_rs))
        self.assertIsNone(man["environment"]["setup_json_differs"])
        self.assertNotIn("different environment", manifest)


if __name__ == "__main__":
    unittest.main()
