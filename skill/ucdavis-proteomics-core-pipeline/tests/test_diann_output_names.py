#!/usr/bin/env python3
"""
A search folder with a '.' in its name (a re-injection's `.01dup2`, a `v1.2`).

DIA-NN names everything it writes beside an output -- `<name>.predicted.speclib`, `<name>.log.txt`,
`<name>.stats.tsv` -- by cutting the --out / --out-lib path at its LAST '.', over the whole path
(remove_extension(), diann.cpp at tag 1.8). On 2026-10-02 (DIA-NN 2.7.0) a one-file search into
`.../species_check/S07_A01.01dup2` passed `--out-lib <out>/diann_lib`: DIA-NN wrote
`species_check/S07_A01.predicted.speclib` and `S07_A01.log.txt` one folder UP, the
library job's guard failed, and the search job never ran. Re-run on HIVE (DIA-NN 2.7.0, job
24560319): the same with `--out-lib <out>/diann_lib.speclib` writes both inside the folder.

The fake DIA-NN here derives names exactly that way (bash `${path%.*}` is remove_extension), so a
generated job that hands DIA-NN an output file name without an extension misplaces its output
here as it did on HIVE.
"""
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
from job_env import job_env  # noqa: E402  (env for running job scripts)

import diann_parallel as dp  # noqa: E402
import test_single_shot_sbatch as sss  # noqa: E402  (the single-shot workspace)

# Prediction: --out-lib and no --f / --lib (the cfg, read by DIA-NN itself, says --predictor).
# Conversion (make_radiant_library.py): --lib + --out-lib, no --f -> a .parquet library.
# The log: from --out when runs are searched, else from --out-lib (diann.cpp 1.8, "Log saved").
FAKE = r"""#!/bin/bash
[ -n "${FAKE_ARGV:-}" ] && printf '%s\n' "$*" >> "$FAKE_ARGV"
out_lib=""; out=""; lib=""; nf=0
while [ $# -gt 0 ]; do
  case "$1" in
    --out-lib) out_lib="$2"; shift ;;
    --out) out="$2"; shift ;;
    --lib) lib="$2"; shift ;;
    --f) nf=$((nf + 1)); shift ;;
  esac
  shift
done
if [ -n "$out_lib" ] && [ "$nf" = 0 ] && [ -z "$lib" ]; then
  head -c 2000 /dev/zero | tr '\0' 'p' > "${out_lib%.*}.predicted.speclib"
elif [ -n "$out_lib" ] && [ "$nf" = 0 ]; then
  head -c 2000 /dev/zero | tr '\0' 'q' > "${out_lib%.*}.parquet"
fi
if [ -n "$out" ] && [ "$nf" -gt 0 ]; then echo log > "${out%.*}.log.txt"
elif [ -n "$out_lib" ]; then echo log > "${out_lib%.*}.log.txt"; fi
exit 0
"""


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)
    return path


def _strays(parent, folder):
    """What DIA-NN put one level up: a name in `parent` that is `folder` cut at a '.'."""
    stem = folder.split(".")[0]
    return sorted(n for n in os.listdir(parent) if n.startswith(stem) and n != folder)


def _diann_outputs(script, diann):
    """(flag, value) for every --out / --out-lib on a DIA-NN command line in a job script."""
    text = open(script).read().replace("\\\n", " ")
    found = []
    for line in text.splitlines():
        if diann not in line:
            continue
        try:
            words = shlex.split(line.strip(), comments=False)
        except ValueError:
            continue
        for i, w in enumerate(words[:-1]):
            if w in ("--out", "--out-lib"):
                found.append((w, words[i + 1]))
    return found


class Helpers(unittest.TestCase):
    def test_the_predicted_library_is_named_as_diann_names_it(self):
        self.assertEqual(dp.predicted_speclib("/s/P2_X.01dup2/diann_lib.speclib"),
                         "/s/P2_X.01dup2/diann_lib.predicted.speclib")
        self.assertEqual(dp.diann_stem("/s/P2_X.01dup2/report.parquet"),
                         "/s/P2_X.01dup2/report")

    def test_an_output_file_name_without_an_extension_is_refused(self):
        for bad in ("/s/P2_X.01dup2/diann_lib", "/s/v1.2/lib", "diann_lib", "/s/.lib"):
            with self.subTest(bad):
                with self.assertRaises(ValueError) as cm:
                    dp.diann_output_path(bad, "--out-lib")
                self.assertIn("--out-lib", str(cm.exception))
        self.assertEqual(dp.diann_output_path("/s/P2_X.01dup2/t${SLURM_ARRAY_TASK_ID}.parquet"),
                         "/s/P2_X.01dup2/t${SLURM_ARRAY_TASK_ID}.parquet")


class SingleShotInADotFolder(unittest.TestCase):
    def test_the_library_job_writes_its_library_inside_the_search_folder(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d, n_runs=1)
            _exe(w.diann, FAKE)
            parent = os.path.join(d, "species_check")
            w.out = os.path.join(parent, "P2_X.01dup2")
            w.generate()
            p = w.run_job(w.lib_job)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertTrue(os.path.isfile(os.path.join(w.out, "diann_lib.predicted.speclib")),
                            os.listdir(w.out))
            self.assertTrue(os.path.isfile(os.path.join(w.out, "diann_lib.log.txt")))
            self.assertEqual(_strays(parent, "P2_X.01dup2"), [])
            # ...and the search job reads the library the library job wrote
            srch = open(w.search_job).read()
            self.assertIn(shlex.quote(os.path.join(w.out, "diann_lib.predicted.speclib")), srch)


class EveryGeneratedDiannOutputHasAnExtension(unittest.TestCase):
    """Every --out / --out-lib a generated job hands DIA-NN, on every route, in a dot folder."""

    def _check(self, scripts, diann):
        seen = []
        for s in scripts:
            for flag, value in _diann_outputs(s, diann):
                seen.append((os.path.basename(s), flag, value))
                name = value.rsplit("/", 1)[-1]
                self.assertIn(".", name[1:], f"{os.path.basename(s)}: {flag} {value}")
        self.assertTrue(seen, scripts)
        return seen

    def test_single_shot_two_job_and_one_step(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d, n_runs=1)
            w.out = os.path.join(d, "species_check", "P2_X.01dup2")
            w.generate()
            seen = self._check([w.lib_job, w.search_job], w.diann)
            self.assertIn("--out-lib", [f for _, f, _ in seen])
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d, n_runs=1)
            w.out = os.path.join(d, "species_check", "P2_X.01dup2")
            w.generate("--one-step")
            self._check([w.job], w.diann)

    def test_the_five_step_chain(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d, n_runs=3)
            out = os.path.join(d, "proj.v1.2", "out.01dup2")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", w.diann, "--raw", *w.runs, "--fasta", w.fasta,
                                "--out", out, "--cfg", w.cfg, "--partition", "low",
                                "--account", "publicgrp"],
                               capture_output=True, text=True, env=job_env(d), timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            scripts = [os.path.join(out, n) for n in sorted(os.listdir(out))
                       if n.endswith(".sbatch")]
            seen = self._check(scripts, w.diann)
            self.assertTrue({"--out", "--out-lib"} <= {f for _, f, _ in seen}, seen)

    def test_the_radiant_chain_predicts_its_library_inside_the_folder(self):
        with tempfile.TemporaryDirectory() as d:
            diann = _exe(os.path.join(d, "fake-diann"), FAKE)
            mzml = [os.path.join(d, f"r{i}.mzML") for i in range(2)]
            for p in mzml:
                open(p, "w").close()
            listing = os.path.join(d, "files.txt")
            with open(listing, "w") as fh:
                fh.write("\n".join(mzml) + "\n")
            fasta, config = os.path.join(d, "db.fasta"), os.path.join(d, "x.radiantConfig")
            for p in (fasta, config):
                open(p, "w").close()
            out = os.path.join(d, "P2_X.01dup2")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "radiant_parallel.py"),
                                "--runtime", "apptainer", "--image", "/x/radiant.sif",
                                "--raw-list", listing, "--fasta", fasta, "--config", config,
                                "--out", out, "--diann", diann, "--partition", "low",
                                "--account", "publicgrp"],
                               capture_output=True, text=True, env=job_env(d), timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            step1 = os.path.join(out, "radiant_parallel", "step1_libpred.sbatch")
            self._check([step1], diann)
            p = subprocess.run(["bash", step1], capture_output=True, text=True,
                               env=job_env(d), cwd=d, timeout=120)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(json.loads(g.stdout)["library"],
                             os.path.join(out, "radiant_parallel",
                                          "diann_predicted_lib.predicted.speclib"))
            self.assertEqual(_strays(d, "P2_X.01dup2"), [])


class RadiantLibraryInADotFolder(unittest.TestCase):
    def test_both_diann_steps_land_in_the_out_dir(self):
        """Prediction, then the speclib -> parquet conversion: both outputs are found where
        make_radiant_library.py looks (it reaches step 3/3), none one level up."""
        with tempfile.TemporaryDirectory() as d:
            diann = _exe(os.path.join(d, "fake-diann"), FAKE)
            fasta = os.path.join(d, "db.fasta")
            open(fasta, "w").close()
            parent = os.path.join(d, "radiant")
            libdir = os.path.join(parent, "lib.v1.2")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_radiant_library.py"),
                                "--diann", diann, "--fasta", fasta, "--out-dir", libdir],
                               capture_output=True, text=True, timeout=60)
            # step 3 reads the (fake) parquet and fails on it; steps 1 and 2 are what is tested
            self.assertIn("[3/3]", p.stdout, p.stdout + p.stderr)
            self.assertEqual(_strays(parent, "lib.v1.2"), [])
            names = os.listdir(libdir)
            self.assertIn("diann_predicted_lib.predicted.speclib", names)
            self.assertIn("radiant_library_src.parquet", names)


if __name__ == "__main__":
    unittest.main()
