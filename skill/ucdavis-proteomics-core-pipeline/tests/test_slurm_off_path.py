#!/usr/bin/env python3
"""
sbatch off PATH (I13, I44): run_search still finds SLURM, and still refuses an inline search on a
login node.

A non-login shell on HIVE (an agent's `ssh hive 'cmd'`, or a shell after `module purge`) has no
sbatch on PATH. run_search.slurm_available() was `shutil.which("sbatch")`, so there:
  * a 6-file DIA-NN cohort was routed to ONE single-shot search ("no SLURM here"), not the 5-step
    chain -- a different order of magnitude in wall time;
  * the hard stop on an inline search ("sbatch present + no SLURM_JOB_ID") did not fire, and the
    search could start on the login node.
SLURM lives at a fixed path on all three HPC@UCD clusters (on_compute.SLURM_BIN_DEFAULT), and
on_compute.on_login_node() is THE definition of a login node; both now decide here too.

Fakes only: a fake sbatch in a temporary HPCCF_SLURM_BIN, a temporary HPCCF cluster file.
"""
import contextlib
import io
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import run_search  # noqa: E402
# Module import, not `from ... import`: a TestCase bound here would be collected twice.
import test_single_shot_mass_accuracy as ss  # noqa: E402  (single-shot harness, no SLURM)

NO_SBATCH_PATH = os.pathsep.join([os.path.dirname(sys.executable), "/usr/bin", "/bin"])
PINNED = "--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n"


class Args:
    """The subset of the argparse namespace parallel_decision reads."""
    parallel_threshold = 5
    no_parallel = False


def fake_slurm_bin(d):
    b = os.path.join(d, "slurm", "bin")
    os.makedirs(b)
    p = os.path.join(b, "sbatch")
    with open(p, "w") as fh:
        fh.write("#!/bin/sh\necho 1\n")
    os.chmod(p, 0o755)
    return b


class SbatchAtTheFixedPath(unittest.TestCase):
    def setUp(self):
        run_search._SBATCH_OFF_PATH_WARNED = False

    def test_fixed_path_sbatch_routes_six_files_to_the_chain(self):
        with tempfile.TemporaryDirectory() as d:
            b = fake_slurm_bin(d)
            cfg = os.path.join(d, "params.cfg")
            with open(cfg, "w") as fh:
                fh.write(PINNED)
            err = io.StringIO()
            with mock.patch.dict(os.environ, {"PATH": NO_SBATCH_PATH, "HPCCF_SLURM_BIN": b}), \
                    contextlib.redirect_stderr(err):
                self.assertIsNone(run_search.shutil.which("sbatch"))
                self.assertEqual(run_search.sbatch_path(), os.path.join(b, "sbatch"))
                use, why = run_search.parallel_decision(
                    "diann", [f"f{i}.d" for i in range(6)], cfg, Args())
                run_search.sbatch_path()
            self.assertTrue(use, why)
            msg = err.getvalue()
            self.assertEqual(msg.count("WARNING"), 1, msg)          # one warning, not one per call
            self.assertIn("sbatch is not on PATH in this shell", msg)
            self.assertIn(f"PATH={b}:$PATH bash", msg)

    def test_no_sbatch_anywhere_is_no_slurm(self):
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.dict(os.environ, {"PATH": NO_SBATCH_PATH,
                                              "HPCCF_SLURM_BIN": os.path.join(d, "none")}):
                self.assertIsNone(run_search.sbatch_path())
                self.assertFalse(run_search.slurm_available())
                cfg = os.path.join(d, "params.cfg")
                with open(cfg, "w") as fh:
                    fh.write(PINNED)
                use, why = run_search.parallel_decision(
                    "diann", [f"f{i}.d" for i in range(6)], cfg, Args())
            self.assertFalse(use)
            self.assertIn("no SLURM here", why)
            self.assertIn(os.path.join(d, "none"), why, "it says where it looked")


class InlineRefusal(unittest.TestCase):
    def _run(self, d, cluster_conf):
        case = ss.SingleShotMassAccTests()
        raws, fasta, cfg, tools, bundle = case._setup(d)
        env = dict(case._env(d), UCD_HPC_CLUSTER_CONF=cluster_conf,
                   HPCCF_SLURM_BIN=os.path.join(d, "no_slurm_here"))
        self.assertNotIn("SLURM_JOB_ID", env)
        out = os.path.join(d, "out")
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                               "--tools", tools, "--bundle", bundle, "--params", cfg,
                               "--fasta", fasta, "--out", out, "--files", *raws,
                               "--engine", "diann", "--threads", "8"],
                              capture_output=True, text=True, env=env, timeout=240)

    def test_inline_search_on_a_login_node_is_refused_without_sbatch_on_path(self):
        with tempfile.TemporaryDirectory() as d:
            conf = os.path.join(d, "cluster_name.conf")
            with open(conf, "w") as fh:
                fh.write("hive\n")
            p = self._run(d, conf)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("REFUSING to run the search inline", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(d, "argv.txt")),
                             "DIA-NN must not have been started")

    def test_laptop_still_runs_inline(self):
        """(pin) No SLURM and no HPCCF cluster file: a laptop searches inline, as before."""
        with tempfile.TemporaryDirectory() as d:
            p = self._run(d, os.path.join(d, "no_cluster_file"))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertNotIn("REFUSING", p.stderr)


if __name__ == "__main__":
    unittest.main()
