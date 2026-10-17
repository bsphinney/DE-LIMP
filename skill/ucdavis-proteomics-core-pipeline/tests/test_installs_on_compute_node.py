#!/usr/bin/env python3
"""setup.sh and ensure_dotnet8.sh install on a compute node, never on a cluster login node.

HPC@UCD's ucd-hpc skill (mandatory on Hive since 2026-10-16) names "large conda solves" as work
that may not run on a login node; downloads and package builds go "inside a job". SKILL.md ran
`hive_exec.sh 'bash ~/proteomics-pipeline/scripts/setup.sh'` -- a micromamba solve, R source
installs and a .NET download on HIVE's login node -- and its printed fixes said to run them
"where there is internet (a login node)". HIVE's compute nodes reach conda-forge, Bioconductor,
CRAN and Microsoft's .NET downloads (srun job 24827776, 2026-10-16), so that reason is gone.

On a login node (sbatch on PATH, SLURM_JOB_ID unset) both scripts now run themselves again
through on_compute.py: one blocking srun with the queue spelled out. ensure_dotnet8.sh only when
it has something to download: reusing an install is a check, and stays where it is.

The stubs: sbatch (a "login node"), sacctmgr (the queue), and an srun that records its arguments
and runs the command inside a pretend job (test_on_compute.FAKE_SRUN).
"""
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import test_setup_limpa as limpa_t                   # noqa: E402
import test_trfp_dotnet_discovery as dn_t            # noqa: E402
from test_on_compute import FAKE_SACCTMGR, FAKE_SRUN, LAB, _exe   # noqa: E402


def slurm_stubs(bindir, with_python=True):
    _exe(os.path.join(bindir, "sbatch"), "#!/bin/sh\nexit 0\n")
    _exe(os.path.join(bindir, "srun"), FAKE_SRUN)
    _exe(os.path.join(bindir, "sacctmgr"), FAKE_SACCTMGR)
    if with_python:
        os.symlink(sys.executable, os.path.join(bindir, "python3"))


def srun_calls(log):
    if not os.path.exists(log):
        return []
    with open(log) as fh:
        return [json.loads(ln) for ln in fh if ln.strip()]


class SetupOnAComputeNode(limpa_t.SetupLimpaGate):
    def login_node(self, with_python=True):
        slurm_stubs(self.sys, with_python)
        for tool in ("tr",):                  # FAKE_SACCTMGR's
            if not os.path.exists(os.path.join(self.sys, tool)):
                os.symlink(limpa_t.shutil.which(tool), os.path.join(self.sys, tool))
        self.srun_log = os.path.join(self.d, "srun.log")
        return {"FAKE_SRUN_LOG": self.srun_log, "FAKE_ASSOC": LAB, "USER": "staff1"}

    def test_an_install_runs_under_srun_with_the_queue_spelled_out(self):
        self.prebuilt_env("1.4.1")
        r, s = self.run_setup(**self.login_node())
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        (call,) = srun_calls(self.srun_log)
        for flag in ("--account=genome-center-grp", "--partition=high", "--time=120", "--mem=8G",
                     "--cpus-per-task=2", "--job-name=pp_setup"):
            self.assertIn(flag, call)
        self.assertEqual(call[-1], limpa_t.setup_t.SETUP)   # setup.sh itself, again
        self.assertIs(s["check_only"], False)                # and it ran, in the "job"
        self.assertIn("compute node", r.stderr)

    def test_the_check_runs_under_srun_too_with_less(self):
        """--check loads R and limpa, nanoparquet, AnnotationDbi, GO.db, org.*.eg.db and starts
        the .NET parser: more than a few seconds of CPU."""
        self.prebuilt_env("1.4.1")
        r, s = self.run_setup("--check", **self.login_node())
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        (call,) = srun_calls(self.srun_log)
        self.assertIn("--time=20", call)
        self.assertIn("--mem=4G", call)
        self.assertEqual(call[-2:], [limpa_t.setup_t.SETUP, "--check"])
        self.assertIs(s["check_only"], True)

    def test_inside_a_job_it_runs_where_it_is(self):
        self.prebuilt_env("1.4.1")
        r, _s = self.run_setup(SLURM_JOB_ID="31", **self.login_node())
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertEqual(srun_calls(self.srun_log), [])

    def test_without_a_way_to_a_compute_node_it_refuses_and_installs_nothing(self):
        self.prebuilt_env("1.2.5")
        env = dict(self.login_node(with_python=False))
        r = limpa_t.subprocess.run([self.bash, limpa_t.setup_t.SETUP], capture_output=True,
                                   text=True, timeout=120,
                                   env={"PATH": self.sys, "HOME": self.d, "PP_HOME": self.home,
                                        "FAKE_R_LOG": self.r_log, **env})
        self.assertEqual(r.returncode, 2, r.stderr)
        self.assertIn("REFUSING", r.stderr)
        self.assertIn("srun", r.stderr)
        self.assertEqual(limpa_t._lines(self.r_log), [])
        self.assertFalse(os.path.exists(os.path.join(self.home, "setup.json")))

    def test_a_compute_node_run_that_did_not_finish_is_exit_3_not_a_stale_setup_json(self):
        """Review of this branch: srun refused (an invalid account) or killed (a time limit) exited 1,
        which already means "the R stack is not OK", and left the old setup.json, whose ready_for
        the agent would then gate on."""
        self.prebuilt_env("1.4.1")
        env = self.login_node()
        r, s = self.run_setup("--check", **env)                  # a first, good run
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertTrue(s["run_id"])                             # the guard's id, written there
        r = limpa_t.subprocess.run([self.bash, limpa_t.setup_t.SETUP, "--check"], capture_output=True,
                                   text=True, timeout=120,
                                   env={"PATH": self.sys, "HOME": self.d, "PP_HOME": self.home,
                                        "FAKE_SRUN_FAIL": "1", **env})
        self.assertEqual(r.returncode, 3, r.stderr)
        self.assertIn("did not finish", r.stderr)
        self.assertIn("not rewritten", r.stderr)

    def test_no_fix_says_to_install_on_a_login_node(self):
        self.prebuilt_env("1.2.5")
        np_state = os.path.join(self.d, "np")
        with open(np_state, "w") as fh:
            fh.write("missing")
        r, s = self.run_setup("--check", FAKE_NP_STATE=np_state)
        text = json.dumps(s) + r.stderr
        self.assertNotIn("(a login node)", text)
        self.assertNotIn("where there is internet", text)
        (note,) = [n for n in s["notes"] if n.startswith("limpa ")]
        self.assertIn("inside a job", note)


class EnsureDotnet8OnAComputeNode(dn_t.EnsureDotnet8):
    TOOLS = dn_t.EnsureDotnet8.TOOLS + ("tr",)

    def setUp(self):
        super().setUp()
        slurm_stubs(self.bin)
        self.srun_log = os.path.join(self.d, "srun.log")
        self.env.update({"FAKE_SRUN_LOG": self.srun_log, "FAKE_ASSOC": LAB, "USER": "staff1"})

    def test_a_download_runs_under_srun_and_the_root_is_still_the_last_line(self):
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], self.dest)
        self.assertEqual(log, ["curl", "install dotnet", "curl", "install aspnetcore"])
        (call,) = srun_calls(self.srun_log)
        for flag in ("--time=30", "--mem=4G", "--job-name=pp_dotnet8", "--account=genome-center-grp"):
            self.assertIn(flag, call)

    def test_reusing_an_install_is_a_check_and_stays_here(self):
        dn_t.fake_dotnet_root(self.dest)
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], self.dest)
        self.assertEqual(srun_calls(self.srun_log), [])
        self.assertEqual(log, [])

    def test_without_python_on_a_cluster_nothing_is_downloaded_and_a_netcore_root_still_serves(self):
        """Review of this branch: with no python3 + on_compute.py the guard exited 1 before step 4,
        so a NETCore >= 8.0.17 root without AspNetCore no longer reached DIA-NN (the header's
        promise: never lose a .raw search that worked before AspNetCore was asked for)."""
        os.remove(os.path.join(self.bin, "python3"))
        dn_t.fake_dotnet_root(self.dest, (dn_t.NETCORE,))
        p, log = self.run_it()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], self.dest)
        self.assertEqual(log, [])                               # no download on the login node
        self.assertIn("nothing is downloaded here", p.stderr)
        self.assertEqual(srun_calls(self.srun_log), [])

    def test_inside_a_job_it_downloads_where_it_is(self):
        p, log = self.run_it(SLURM_JOB_ID="12")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(srun_calls(self.srun_log), [])
        self.assertIn("install dotnet", log)


# the parent classes' own tests run in their own modules; not again here
for _cls in (SetupOnAComputeNode, EnsureDotnet8OnAComputeNode):
    for _name in [n for n in dir(_cls) if n.startswith("test_")]:
        if _name not in _cls.__dict__:
            setattr(_cls, _name, None)
del _cls, _name                      # a TestCase bound at module level is collected again

if __name__ == "__main__":
    unittest.main()
