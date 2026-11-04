#!/usr/bin/env python3
"""
The parallel Radiant route gives Radiant a TSV library, not DIA-NN 2.x's .predicted.speclib (S8).

radiant_parallel.py's step 1 ran DIA-NN's predictor and handed step 2 the
`diann_predicted_lib.predicted.speclib` it wrote. Radiant rejects every DIA-NN 2.x speclib
(format v-10/-11; its reader supports v >= -3: "ERROR: version is not supported11"), so every
array task of a multi-file Radiant search on SLURM would fail. The single-job route already
built a TSV with make_radiant_library.py (predict -> parquet -> renamed TSV).

Now step 1 runs the same job lines as run_search.radiant_library_job() -- one argv builder,
run_search.radiant_library_lines(), never a copy (rule 3) -- into the same <out>/radiant_lib/, and
steps 2 and 3 read the .tsv it writes.

Hermetic: radiant_parallel.py writes scripts only; nothing is submitted.
"""
import json
import os
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
from job_env import job_env  # noqa: E402

import run_search  # noqa: E402

DIANN = "/opt/diann-2.7.0/diann-linux"


class ParallelRadiantLibrary(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        d = self.d = os.path.realpath(self._tmp.name)
        runs = [os.path.join(d, "raw", f"r{i}.mzML") for i in range(3)]
        os.makedirs(os.path.dirname(runs[0]))
        for r in runs:
            open(r, "w").close()
        self.listing = os.path.join(d, "files.txt")
        with open(self.listing, "w") as fh:
            fh.write("\n".join(runs) + "\n")
        self.fasta = os.path.join(d, "db", "search.fasta")
        os.makedirs(os.path.dirname(self.fasta))
        with open(self.fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEK\n")
        self.config = os.path.join(d, "default.radiantConfig")
        open(self.config, "w").close()
        self.out = os.path.join(d, "search out")

    def generate(self, *extra):
        g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "radiant_parallel.py"),
                            "--runtime", "apptainer", "--image", "/x/radiant.sif",
                            "--raw-list", self.listing, "--fasta", self.fasta,
                            "--config", self.config, "--out", self.out, "--diann", DIANN,
                            "--cys-alkylation", "nem", "--partition", "high",
                            "--account", "genome-center-grp", *extra],
                           capture_output=True, text=True, env=job_env(self.d), timeout=120)
        self.assertEqual(g.returncode, 0, g.stderr)
        D = os.path.join(self.out, "radiant_parallel")

        def read(p):
            with open(p) as fh:
                return fh.read()
        return (json.loads(g.stdout), read(os.path.join(D, "step1_libpred.sbatch")),
                read(os.path.join(D, "step2_search.sbatch")),
                read(os.path.join(self.out, "radiant_results", "fulcrum_rescore.toml")))

    def tsv(self):
        return os.path.join(self.out, "radiant_lib", "radiant_library.tsv")

    def test_step1_builds_a_tsv_library_with_make_radiant_library(self):
        info, step1, _, _ = self.generate()
        self.assertEqual(info["library"], self.tsv())
        self.assertIn("make_radiant_library.py", step1)
        self.assertIn("--out-dir " + shlex.quote(os.path.join(self.out, "radiant_lib")), step1)
        self.assertIn("--cys-alkylation nem", step1)
        self.assertNotIn("--out-lib", step1, "DIA-NN is run by make_radiant_library.py")
        self.assertNotIn(".predicted.speclib", step1)
        # the guard checks for the TSV it hands step 2
        self.assertIn(f"[ -s {shlex.quote(self.tsv())} ]", step1)

    def test_step1_is_the_single_job_routes_lines(self):
        """Rule 3: one argv builder. Every line run_search.radiant_library_job() writes for
        this library (threads as step 1 runs it) is in step 1, in order."""
        _, step1, _, _ = self.generate()
        tsv, lines = run_search.radiant_library_lines(DIANN, self.fasta,
                                                      os.path.join(self.out, "radiant_lib"),
                                                      16, "nem")
        self.assertEqual(tsv, self.tsv())
        body = step1.splitlines()
        at = [body.index(ln) for ln in lines]
        self.assertEqual(at, sorted(at), "in the builder's order")
        # ...and the single-job route's lines are that builder's, with its own cysteine choice
        # (a config with no make_presets.py record: cam)
        self.assertEqual(run_search.radiant_library_job({"diann": DIANN}, self.fasta, self.out,
                                                        16, self.config),
                         run_search.radiant_library_lines(
                             DIANN, self.fasta, os.path.join(self.out, "radiant_lib"), 16,
                             "cam"))

    def test_step2_reads_the_tsv(self):
        _, _, step2, toml = self.generate()
        self.assertIn("/radiant_library.tsv ", step2)
        self.assertNotIn("speclib", step2)
        lib = [ln for ln in toml.splitlines() if ln.startswith("library = ")]
        self.assertEqual(len(lib), 1, toml)
        self.assertTrue(json.loads(lib[0].split("= ", 1)[1]).endswith("/radiant_library.tsv"))
        self.assertNotIn("speclib", toml)

    def test_a_given_library_is_used_as_given(self):
        """(pin) --library skips step 1, and steps 2 and 3 read that file."""
        given = os.path.join(self.d, "lib", "my_library.tsv")
        os.makedirs(os.path.dirname(given))
        open(given, "w").close()
        g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "radiant_parallel.py"),
                            "--runtime", "apptainer", "--image", "/x/radiant.sif",
                            "--raw-list", self.listing, "--fasta", self.fasta,
                            "--config", self.config, "--out", self.out, "--library", given,
                            "--partition", "high", "--account", "genome-center-grp"],
                           capture_output=True, text=True, env=job_env(self.d), timeout=120)
        self.assertEqual(g.returncode, 0, g.stderr)
        info = json.loads(g.stdout)
        self.assertEqual(info["library"], given)
        self.assertNotIn("step1_libpred", info["scripts"])


if __name__ == "__main__":
    unittest.main()
