#!/usr/bin/env python3
"""setup.sh builds an env whose limpa is >= 1.4.0, solves with its own channels, and says
plainly when limpa is too old.

HIVE brettsp, 2026-09-27: (1) sbatch 24154212 -- ~/.condarc with `defaults` and
channel_priority strict left libmamba only defaults' r-statmod (R <= 4.3), and bioconductor-limma
3.66 was unsolvable: setup.sh now passes --override-channels. (2) sbatch 24154221 -- the env had
bioconda's only limpa, 1.2.5, and run_de.R died in readDIANN(annotation.columns =), limpa >= 1.4.0.
setup.sh now installs limpa from the Bioconductor 3.23 source repository (pure R), checks
packageVersion("limpa") >= 1.4.0 last, reports it in setup.json, and exits 1 when it is short.

A stub micromamba and a stub Rscript stand in (no conda, no network): the Rscript stub keeps
"the installed limpa version" in a file, answers the version check from it, and "installs"
1.4.0 unless told the download fails.
"""
import json
import os
import shutil
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import test_setup_diann_readiness as setup_t          # noqa: E402
import test_thermo_raw_detection as trfp             # noqa: E402
from test_trfp_dotnet_discovery import (FAKE_ENSURE, FAKE_MICROMAMBA, _exe, _lines,   # noqa: E402
                                        fake_dotnet_root)

FAKE_RSCRIPT = r'''#!/bin/bash
echo "$*" >> "$FAKE_R_LOG"
v="$(cat "$FAKE_LIMPA_STATE" 2>/dev/null)"
case "$*" in
  *install.packages*)
    [ -n "${FAKE_LIMPA_INSTALL_FAIL:-}" ] && { echo "cannot open URL bioconductor.org" >&2; exit 1; }
    echo "1.4.0" > "$FAKE_LIMPA_STATE"; exit 0 ;;
  *"cat(tryCatch"*) printf '%s' "$v"; exit 0 ;;
  *requireNamespace*)
    case "$v" in 1.[4-9]*|1.[1-9][0-9]*|[2-9].*) exit 0 ;; *) exit 1 ;; esac ;;
esac
exit 0
'''


class SetupLimpaGate(setup_t.SetupCheckHarness):
    def setUp(self):
        super().setUp()
        for tool in ("chmod", "ln"):
            os.symlink(shutil.which(tool), os.path.join(self.sys, tool))
        _exe(os.path.join(self.sys, "micromamba"), FAKE_MICROMAMBA)
        self.ensure = _exe(os.path.join(self.d, "ensure_dotnet8_stub.sh"), FAKE_ENSURE)
        self.dotnet_root = fake_dotnet_root(os.path.join(self.d, ".proteomics-pipeline", "dotnet8"))
        self.mm_log = os.path.join(self.d, "mm.log")
        self.r_log = os.path.join(self.d, "r.log")
        self.state = os.path.join(self.d, "limpa_version")
        self.prefix = os.path.join(self.home, "micromamba", "envs", "proteomics-pipeline")

    def prebuilt_env(self, limpa="1.2.5"):
        """An env that exists already (create_env will not run): python + the Rscript stub."""
        os.makedirs(os.path.join(self.prefix, "bin"))
        os.makedirs(os.path.join(self.prefix, "conda-meta"))
        os.symlink(sys.executable, os.path.join(self.prefix, "bin", "python"))
        _exe(os.path.join(self.prefix, "bin", "Rscript"), FAKE_RSCRIPT)
        for meta in ("pythonnet-3.1.0-pyhd8ed1ab_0.json", "pandas-2.3.2-py311h0000000_0.json"):
            with open(os.path.join(self.prefix, "conda-meta", meta), "w") as fh:
                fh.write("{}")
        with open(self.state, "w") as fh:
            fh.write(limpa)

    def run_setup(self, *args, **extra):
        env = {"PATH": self.sys, "HOME": self.d, "PP_HOME": self.home, "QUOBYTE_DIR": "",
               "THERMORAWFILEPARSER_SHARED": "", "PROTEOMICS_DOTNET_SYSTEM_ROOTS": "",
               "FAKE_MM_LOG": self.mm_log, "FAKE_PY": sys.executable, "FAKE_TRFP": trfp.FAKE,
               "PROTEOMICS_ENSURE_DOTNET8": self.ensure, "FAKE_ENSURE_LOG": os.path.join(self.d, "e.log"),
               "FAKE_ENSURE_ROOT": self.dotnet_root, "FAKE_R_LOG": self.r_log,
               "FAKE_LIMPA_STATE": self.state}
        env.update(extra)
        r = subprocess.run([self.bash, setup_t.SETUP, *args], capture_output=True, text=True,
                           env=env, timeout=300)
        with open(os.path.join(self.home, "setup.json")) as fh:
            return r, json.load(fh)

    def test_an_old_limpa_is_upgraded_from_bioconductor_3_23(self):
        self.prebuilt_env("1.2.5")
        r, s = self.run_setup()
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        installs = [c for c in _lines(self.r_log) if "install.packages" in c]
        self.assertEqual(len(installs), 1, _lines(self.r_log))
        self.assertIn("https://bioconductor.org/packages/3.23/bioc", installs[0])
        self.assertIn("type = 'source'", installs[0])
        self.assertEqual(s["limpa"], {"version": "1.4.0", "required": ">= 1.4.0", "ok": True,
                                      "source": "Bioconductor 3.23"})
        self.assertNotIn("ERROR", r.stderr)

    def test_a_current_limpa_is_left_alone(self):
        self.prebuilt_env("1.4.1")
        r, s = self.run_setup()
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertFalse(any("install.packages" in c for c in _lines(self.r_log)))
        self.assertIs(s["limpa"]["ok"], True)
        self.assertEqual(s["limpa"]["version"], "1.4.1")

    def test_a_failed_upgrade_is_an_error_with_the_fix(self):
        self.prebuilt_env("1.2.5")
        r, s = self.run_setup(FAKE_LIMPA_INSTALL_FAIL="1")
        self.assertEqual(r.returncode, 1, r.stderr[-2000:])
        self.assertIn("[setup] ERROR: limpa 1.2.5 in", r.stderr)
        self.assertIn("needs limpa >= 1.4.0 (Bioconductor 3.23)", r.stderr)
        self.assertIn("install.packages('limpa', repos = 'https://bioconductor.org/packages/3.23/bioc'",
                      r.stderr)
        self.assertIs(s["limpa"]["ok"], False)
        self.assertEqual(s["limpa"]["version"], "1.2.5")
        self.assertTrue(any("limpa 1.2.5 is older than 1.4.0" in n for n in s["notes"]), s["notes"])

    def test_check_only_reports_and_installs_nothing(self):
        self.prebuilt_env("1.2.5")
        r, s = self.run_setup("--check")
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertFalse(any("install.packages" in c for c in _lines(self.r_log)))
        self.assertIs(s["limpa"]["ok"], False)
        self.assertIn("ERROR: limpa 1.2.5", r.stderr)

    def test_the_solve_names_its_own_channels(self):
        # a fresh env: the stub micromamba creates it (its Rscript answers every check "yes")
        r, s = self.run_setup()
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        calls = _lines(self.mm_log)
        creates = [c for c in calls if c.startswith("create")]
        installs = [c for c in calls if c.startswith("install")]
        self.assertEqual(len(creates), 1, calls)
        for c in creates + installs:
            self.assertIn("--override-channels -c conda-forge -c bioconda", c)
        for pkg in ("bioconductor-limma", "r-statmod", "r-data.table", "r-nanoparquet"):
            self.assertIn(f" {pkg} ", f" {creates[0]} ")


if __name__ == "__main__":
    unittest.main()
