#!/usr/bin/env python3
"""setup.sh builds an env whose limpa is >= 1.4.0, solves with its own channels, and says
plainly when limpa is too old.

HIVE bstaff, 2026-09-27: (1) sbatch 24154212 -- ~/.condarc with `defaults` and
channel_priority strict left libmamba only defaults' r-statmod (R <= 4.3), and bioconductor-limma
3.66 was unsolvable: setup.sh now passes --override-channels. (2) sbatch 24154221 -- the env had
bioconda's only limpa, 1.2.5, and run_de.R died in readDIANN(annotation.columns =), limpa >= 1.4.0.
setup.sh now installs limpa from the Bioconductor 3.23 source repository (pure R), checks
packageVersion("limpa") >= 1.4.0 last, reports it in setup.json, and exits 1 when it is short.

A stub micromamba and a stub Rscript stand in (no conda, no network): the Rscript stub keeps
"the installed limpa version" in a file, answers the version check from it, and "installs"
1.4.0 unless told the download fails.

nanoparquet (2.8.0 report): limpa reads DIA-NN's report.parquet with it (EListFromLongFormatFile:
"nanoparquet package required but is not installed (or can't be loaded)"), but it is only in
limpa's Suggests, so installing limpa never brings it, and r-nanoparquet was in create_env from
2.8.0 only -- which runs for a MISSING env. An older env never got it, `setup.sh --check` did not
look, and CRAN's nanoparquet compiles C++, which a login node without a compiler cannot. setup.sh
now checks it loads, installs conda-forge's r-nanoparquet into an existing env that lacks it, and
says so with that fix. The Rscript stub answers the nanoparquet check from FAKE_NP_STATE.
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
  *getRversion*) printf '%s' "${FAKE_R_VERSION:-}"; exit 0 ;;
  *install.packages*nanoparquet*)
    [ -n "${FAKE_NP_CRAN_FAIL:-}" ] && { echo "compilation failed for package 'nanoparquet'" >&2; exit 1; }
    echo present > "$FAKE_NP_STATE"; exit 0 ;;
  *install.packages*)
    [ -n "${FAKE_LIMPA_INSTALL_FAIL:-}" ] && { echo "cannot open URL bioconductor.org" >&2; exit 1; }
    echo "1.4.0" > "$FAKE_LIMPA_STATE"; exit 0 ;;
  *"cat(tryCatch"*) printf '%s' "$v"; exit 0 ;;
  *packageDescription*)     # what this limpa's DESCRIPTION Suggests for reading parquet
    [ -n "${FAKE_LIMPA_SUGGESTS:-}" ] && { printf '%s' "$FAKE_LIMPA_SUGGESTS"; exit 0; }
    case "$v" in 1.[4-9]*|1.[1-9][0-9]*|[2-9].*) printf nanoparquet ;; 1.*) printf arrow ;; esac
    exit 0 ;;
  *nanoparquet*)
    [ -z "${FAKE_NP_STATE:-}" ] && exit 0
    [ "$(cat "$FAKE_NP_STATE" 2>/dev/null)" = present ] && exit 0 || exit 1 ;;
  *"'arrow'"*) [ -n "${FAKE_ARROW_MISSING:-}" ] && exit 1; exit 0 ;;
  *requireNamespace*)
    case "$v" in 1.[4-9]*|1.[1-9][0-9]*|[2-9].*) exit 0 ;; *) exit 1 ;; esac ;;
esac
exit 0
'''


# conda's R names its compilers in its Makeconf (`R CMD config CC`): conda-forge's
# x86_64-conda-linux-gnu-cc etc. Stand-ins: the names are answered here; whether they RUN depends
# on the env's bin (compilers()).
FAKE_R = r'''#!/bin/bash
[ "$1 $2" = "CMD config" ] || exit 1
case "$3" in CC) echo "fakecc -std=gnu23" ;; CXX) echo "fakecxx -std=gnu++20" ;; MAKE) echo make ;; *) exit 1 ;; esac
'''


# the stub micromamba, plus: installing r-nanoparquet makes the Rscript stub's check pass
FAKE_MICROMAMBA_NP = FAKE_MICROMAMBA.replace(
    "             pandas) ",
    '             r-nanoparquet) echo present > "$FAKE_NP_STATE" ;;\n             pandas) ')


class SetupLimpaGate(setup_t.SetupCheckHarness):
    def setUp(self):
        super().setUp()
        for tool in ("chmod", "ln"):
            os.symlink(shutil.which(tool), os.path.join(self.sys, tool))
        _exe(os.path.join(self.sys, "micromamba"), FAKE_MICROMAMBA_NP)
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
        for meta in ("pythonnet-3.1.0-pyhd8ed1ab_0.json", "pandas-2.3.2-py311h0000000_0.json",
                     "google-genai-2.31.0-pyh0000000_0.json"):     # setup.sh step 2d's
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
               "FAKE_LIMPA_STATE": self.state, "FAKE_NP_STATE": ""}
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

    # ---- nanoparquet -------------------------------------------------------------------------
    def np_missing(self):
        np_state = os.path.join(self.d, "nanoparquet_state")
        with open(np_state, "w") as fh:
            fh.write("missing")
        return {"FAKE_NP_STATE": np_state}

    def test_check_says_nanoparquet_is_missing_with_the_conda_forge_fix(self):
        self.prebuilt_env("1.4.1")
        r, s = self.run_setup("--check", **self.np_missing())
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertIs(s["nanoparquet"]["ok"], False)
        self.assertIs(s["ready_for"]["de"], False)
        note = [n for n in s["notes"] if "nanoparquet" in n]
        self.assertEqual(len(note), 1, s["notes"])
        self.assertIn("--freeze-installed --override-channels -c conda-forge -c bioconda "
                      "r-nanoparquet", note[0])
        self.assertIn(self.prefix, note[0])
        self.assertNotIn("install.packages", note[0][note[0].index("Fix: "):])   # never CRAN
        self.assertEqual(s["nanoparquet"]["fix"], note[0][note[0].index("Fix: ") + 5:].rstrip("."))
        self.assertFalse([c for c in _lines(self.mm_log) if c.startswith("install")])
        self.assertIn("ERROR: nanoparquet", r.stderr)

    def test_setup_installs_it_from_conda_forge_into_an_existing_env(self):
        self.prebuilt_env("1.4.1")
        r, s = self.run_setup(**self.np_missing())
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        installs = [c for c in _lines(self.mm_log) if c.startswith("install")]
        (np,) = [c for c in installs if " r-nanoparquet" in f" {c} "]
        self.assertIn("--freeze-installed --override-channels -c conda-forge", np)
        self.assertTrue(np.endswith(" -c bioconda r-nanoparquet"), np)    # its own call
        self.assertFalse(any("nanoparquet" in c and "install.packages" in c
                             for c in _lines(self.r_log)))
        self.assertIs(s["nanoparquet"]["ok"], True)
        self.assertIs(s["ready_for"]["de"], True)
        self.assertFalse([n for n in s["notes"] if "nanoparquet" in n], s["notes"])

    def test_a_loadable_nanoparquet_is_left_alone(self):
        self.prebuilt_env("1.4.1")
        r, s = self.run_setup()
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertFalse(any("r-nanoparquet" in c for c in _lines(self.mm_log)))
        self.assertIs(s["nanoparquet"]["ok"], True)

    def test_when_it_cannot_be_installed_setup_fails_with_the_fix(self):
        self.prebuilt_env("1.4.1")
        r, s = self.run_setup(FAKE_MM_INSTALL_FAIL="1", **self.np_missing())
        self.assertEqual(r.returncode, 1, r.stderr[-2000:])
        self.assertIn("ERROR: nanoparquet", r.stderr)
        self.assertIn("--freeze-installed --override-channels -c conda-forge -c bioconda r-nanoparquet",
                      r.stderr)
        self.assertIs(s["nanoparquet"]["ok"], False)
        self.assertIs(s["ready_for"]["de"], False)


    # ---- nanoparquet, 2.12 review: never a solve that can downgrade R; R 4.6 has no conda build --
    def r46_env(self, compilers=True):
        """An env with R 4.6.1 (FAKE_R_VERSION) and nanoparquet missing; its R names compilers that
        run here only when `compilers`."""
        self.prebuilt_env("1.4.1")
        _exe(os.path.join(self.prefix, "bin", "R"), FAKE_R)
        if compilers:
            for c in ("fakecc", "fakecxx", "make"):
                _exe(os.path.join(self.prefix, "bin", c), "#!/bin/sh\nexit 0\n")
        return dict(self.np_missing(), FAKE_R_VERSION="4.6.1")

    def test_the_frozen_install_is_the_fix_it_prints(self):
        """NANOPARQUET_FIX is the command setup.sh runs: --check prints it, a run runs it."""
        self.prebuilt_env("1.4.1")
        env = dict(self.np_missing(), FAKE_MM_INSTALL_FAIL="1")
        _r, s = self.run_setup("--check", **env)
        fix = s["nanoparquet"]["fix"]
        self.run_setup(**env)
        (np,) = [c for c in _lines(self.mm_log) if " r-nanoparquet" in f" {c} "]
        self.assertTrue(fix.endswith(" " + np), (fix, np))
        self.assertNotIn("thermorawfileparser", np)

    def test_r46_is_never_given_to_conda_and_compiles_from_cran(self):
        r, s = self.run_setup(**self.r46_env())
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertFalse([c for c in _lines(self.mm_log) if "r-nanoparquet" in c],
                         "an R 4.6 env was handed to conda for r-nanoparquet (R 4.6.1 -> 4.5.3)")
        (cran,) = [c for c in _lines(self.r_log) if "install.packages('nanoparquet'" in c]
        self.assertIn("repos = 'https://cloud.r-project.org', type = 'source'", cran)
        self.assertIs(s["nanoparquet"]["ok"], True)
        self.assertIs(s["ready_for"]["de"], True)
        self.assertEqual(s["nanoparquet"]["fix"], "")

    def test_r46_without_compilers_says_install_them_then_cran(self):
        r, s = self.run_setup(**self.r46_env(compilers=False))
        self.assertEqual(r.returncode, 1, r.stderr[-2000:])
        self.assertFalse([c for c in _lines(self.mm_log) if "r-nanoparquet" in c])
        self.assertFalse([c for c in _lines(self.r_log) if "install.packages('nanoparquet'" in c],
                         "CRAN tried without a compiler (x86_64-conda-linux-gnu-c++: not found)")
        fix = s["nanoparquet"]["fix"]
        self.assertIn("--freeze-installed --override-channels -c conda-forge -c bioconda c-compiler "
                      "cxx-compiler make && ", fix)
        self.assertIn(f"{self.prefix}/bin/Rscript -e \"install.packages('nanoparquet', repos = "
                      f"'https://cloud.r-project.org', type = 'source')\"", fix)
        self.assertIn(f'PATH="{self.prefix}/bin:$PATH"', fix)
        note = [n for n in s["notes"] if "nanoparquet" in n]
        self.assertIn("R 4.6.1: conda-forge has no r-nanoparquet for it", note[0])
        self.assertIn(fix, r.stderr)

    def test_r46_check_prints_the_cran_line_when_compilers_run(self):
        r, s = self.run_setup("--check", **self.r46_env())
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertTrue(s["nanoparquet"]["fix"].startswith(f'PATH="{self.prefix}/bin:$PATH" '),
                        s["nanoparquet"]["fix"])
        self.assertFalse([c for c in _lines(self.r_log) if "install.packages" in c])

    def test_r45_falls_back_to_cran_when_conda_cannot_add_it(self):
        env = dict(self.r46_env(), FAKE_R_VERSION="4.5.3", FAKE_MM_INSTALL_FAIL="1")
        r, s = self.run_setup(**env)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertTrue([c for c in _lines(self.mm_log) if "r-nanoparquet" in c])
        self.assertTrue([c for c in _lines(self.r_log) if "install.packages('nanoparquet'" in c])
        self.assertIs(s["nanoparquet"]["ok"], True)

    # ---- ready_for: a search can run, the DE cannot -- and says why -------------------------
    def test_no_env_with_a_search_ready_says_why_de_is_not(self):
        for tool in ("Rscript", "sage"):            # outside the env: enough for a search
            _exe(os.path.join(self.sys, tool), "#!/bin/sh\nexit 0\n")
        r, s = self.run_setup("--check")
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertEqual((s["ready_for"]["dda"], s["ready_for"]["de"]), (True, False))
        note = [n for n in s["notes"] if n.startswith("ready_for.de is false")]
        self.assertEqual(len(note), 1, s["notes"])
        self.assertIn(self.prefix, note[0])
        self.assertIn("setup.sh (without --check)", note[0])

    def test_a_ready_de_has_no_such_note(self):
        self.prebuilt_env("1.4.1")
        _r, s = self.run_setup("--check")
        self.assertIs(s["ready_for"]["de"], True)
        self.assertFalse([n for n in s["notes"] if n.startswith("ready_for.de is false")])


if __name__ == "__main__":
    unittest.main()
