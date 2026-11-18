#!/usr/bin/env python3
"""
2.15.1 review (blocker): setup.sh's conda pins must be solvable together, and limpa must come
from the route that can actually give limpa >= 1.4.0.

The 2.15.1 audit pinned bioconductor-limma=3.68.5 while create_env still asked for
bioconductor-limpa -- bioconda's only build, 1.2.5, which requires bioconductor-limma
>=3.66.0,<3.67.0. Every fresh install failed at the solve (conda --dry-run, osx-arm64 and
linux-64: LibMambaUnsatisfiableError). limpa >= 1.4.0 comes only from step 2a (Bioconductor
3.23 source, pure R). These static checks keep the conda list and that route consistent:

  * create_env's list names no bioconductor-limpa (bioconda has nothing >= 1.4.0, and its 1.2.5
    cannot be solved with limma 3.68.5);
  * limpa's own Depends/Imports (limma, data.table, statmod) are all in the conda list, so the
    source install needs nothing from CRAN that would compile;
  * step 2a installs limpa from Bioconductor and checks packageVersion >= 1.4.0, and the pin is
    >= that floor.

A real solve runs only when asked (PROTEOMICS_TEST_CONDA_SOLVE=1 and conda on PATH): it needs the
network and a minute per platform.
"""
import json
import os
import re
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
SETUP = os.path.join(SCRIPTS, "setup.sh")

# bioconda's bioconductor-limpa builds and the limma range each requires (anaconda.org file
# metadata, 2026-11-17). If a build >= 1.4.0 with a compatible limma ever appears, the conda
# route may come back -- add it here and change the test, not the other way round.
BIOCONDA_LIMPA_LIMMA_RANGE = {"1.2.5": ((3, 66, 0), (3, 67, 0))}
# limpa 1.4.0's DESCRIPTION (Bioconductor 3.23): Depends limma; Imports methods, stats,
# data.table, statmod (methods/stats ship with R)
LIMPA_NEEDS_FROM_CONDA = ("bioconductor-limma", "r-data.table", "r-statmod")


def _src():
    with open(SETUP) as fh:
        return fh.read()


def _pins(src):
    return {m.group(1): m.group(2) for m in
            re.finditer(r'^(PIN_[A-Z_]+)="\$\{PIN_[A-Z_]+:-([^}]+)\}"', src, re.M)}


def create_env_packages(src=None):
    """create_env's `pkgs` array as conda specs, PIN_* expanded (comments dropped)."""
    src = src or _src()
    body = re.search(r"create_env\(\) \{(.*?)\n\}", src, re.S).group(1)
    # the array's lines, comments dropped, up to the line that closes it
    lines = []
    for ln in body[body.index("local pkgs=(") + len("local pkgs=("):].splitlines():
        code = ln.split("#", 1)[0].rstrip()
        if code.endswith(")"):
            lines.append(code[:-1])
            break
        lines.append(code)
    arr = "\n".join(lines)
    pins = _pins(src)
    out = []
    for tok in re.findall(r'"[^"]*"|\S+', arr):
        tok = tok.strip('"')
        for k, v in pins.items():
            tok = tok.replace("$" + k, v)
        out.append(tok)
    linux = re.search(r'if \[ "\$OS" = "linux" \]; then pkgs\+=\(([^)]*)\)', body)
    return out, (linux.group(1).split() if linux else [])


def _ver(s):
    return tuple(int(x) for x in re.findall(r"\d+", s)[:3])


class SetupPinsAreSolvable(unittest.TestCase):
    def setUp(self):
        self.src = _src()
        self.pkgs, self.linux_pkgs = create_env_packages(self.src)
        self.names = {re.split(r"[=<>]", p, maxsplit=1)[0] for p in self.pkgs}

    def test_the_parser_reads_the_list(self):
        self.assertIn("bioconductor-limma=3.68.5", self.pkgs)
        self.assertIn("sage-proteomics=0.14.7", self.pkgs)
        self.assertIn("proteowizard", self.linux_pkgs)

    def test_no_conda_limpa_beside_the_limma_pin(self):
        self.assertNotIn("bioconductor-limpa", self.names,
                         "bioconda's limpa (1.2.5) cannot be solved with the limma pin; limpa "
                         "comes from step 2a")
        limma = _ver(_pins(self.src)["PIN_LIMMA"])
        for lv, (lo, hi) in BIOCONDA_LIMPA_LIMMA_RANGE.items():
            # the reason it is out: every bioconda limpa needs a limma the pin is outside of
            self.assertFalse(lo <= limma < hi, f"limpa {lv} would solve with limma {limma}: "
                             "re-check whether the conda route is wanted again")

    def test_limpa_dependencies_come_from_conda(self):
        for dep in LIMPA_NEEDS_FROM_CONDA:
            self.assertIn(dep, self.names, f"limpa needs {dep}; the source install must not "
                          "have to fetch (and compile) it from CRAN")

    def test_step_2a_is_the_limpa_route_and_meets_the_floor(self):
        self.assertRegex(self.src, r'LIMPA_MIN="1\.4\.0"')
        self.assertIn('LIMPA_PIN="$PIN_LIMPA"', self.src)
        self.assertGreaterEqual(_ver(_pins(self.src)["PIN_LIMPA"]), (1, 4, 0))
        self.assertIn("install.packages('limpa', repos = c(BioCsoft = '$LIMPA_REPO'", self.src)
        self.assertIn("packageVersion('limpa') >= '$LIMPA_MIN'", self.src)
        self.assertRegex(self.src, r'LIMPA_BIOC="\$\{PROTEOMICS_LIMPA_BIOC:-3\.23\}"')

    def test_the_fallback_hint_carries_the_pin_and_no_limpa(self):
        hint = re.search(r'Environment solve failed\. Try: ([^"]*)"', self.src).group(1)
        self.assertNotIn("bioconductor-limpa", hint)
        self.assertIn("bioconductor-limma=$PIN_LIMMA", hint)

    @unittest.skipUnless(os.environ.get("PROTEOMICS_TEST_CONDA_SOLVE") == "1"
                         and shutil.which("conda"),
                         "set PROTEOMICS_TEST_CONDA_SOLVE=1 (needs conda + network)")
    def test_real_dry_run_solve(self):
        for sub in ("linux-64", "osx-arm64"):
            pk = self.pkgs + (self.linux_pkgs if sub.startswith("linux") else [])
            env = dict(os.environ, CONDA_SUBDIR=sub, CONDA_OVERRIDE_GLIBC="2.28",
                       CONDA_OVERRIDE_OSX="13.0")
            prefix = os.path.join(tempfile.mkdtemp(), "never")
            r = subprocess.run(["conda", "create", "--dry-run", "--json", "-p", prefix,
                                "--override-channels", "-c", "conda-forge", "-c", "bioconda"]
                               + pk, capture_output=True, text=True, env=env, timeout=900)
            d = json.loads(r.stdout or "{}")
            self.assertTrue("actions" in d, f"{sub}: {d.get('error', r.stderr)[:500]}")
            got = {p["name"]: p["version"] for p in d["actions"]["LINK"]}
            self.assertEqual(got.get("bioconductor-limma"), "3.68.5", sub)


if __name__ == "__main__":
    unittest.main()
