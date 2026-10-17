#!/usr/bin/env python3
"""on_compute.py runs a command on a SLURM compute node, never on the cluster's login node.

HPC@UCD's agent skill (ucd-hpc, mandatory on Hive, Farm and Franklin since 2026-10-16): "No
analysis on the login node ... anything past a few seconds of CPU or much I/O -- runs in srun
(blocking) or sbatch, and you wait", and every request names --account, --partition, --time and
memory. In hive_remote the skill ran nearly every script on HIVE's login node through
hive_exec.sh: core_submission.py deliver's dry run walked the Flinders service tree there for
about 4 minutes (a re-analysis, 2026-10-16), and setup.sh solved conda environments there.

on_compute.py is the one way to send a command to a compute node: one blocking srun, its queue
from run_search.slurm_queue (the skill's one queue rule) and spelled out. Inside a job, or with
no SLURM (a laptop), the command runs where it is. `hive_exec.sh --job '<cmd>'` runs it on HIVE
from a laptop; ensure_compute_node() lets a heavy script move itself.

No test needs SLURM: sbatch, srun and sacctmgr are stubs on PATH. The srun stub records its
arguments and then runs the command it was given, inside a pretend job (SLURM_JOB_ID set).
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
ON_COMPUTE = os.path.join(SCRIPTS, "on_compute.py")
HIVE_EXEC = os.path.join(SCRIPTS, "hive_exec.sh")
sys.path.insert(0, SCRIPTS)

# srun: one line of JSON per call (its argv), then the command after the options, in a "job"
FAKE_SRUN = r'''#!/usr/bin/env python3
import json, os, sys
with open(os.environ["FAKE_SRUN_LOG"], "a") as fh:
    fh.write(json.dumps(sys.argv[1:]) + "\n")
if os.environ.get("FAKE_SRUN_ENV_LOG"):
    with open(os.environ["FAKE_SRUN_ENV_LOG"], "a") as fh:
        fh.write(os.environ.get("PROTEOMICS_ON_COMPUTE", "") + "\n")
args = sys.argv[1:]
while args and args[0].startswith("--"):
    args.pop(0)
if os.environ.get("FAKE_SRUN_NORUN"):          # record the call only
    sys.exit(0)
if os.environ.get("FAKE_SRUN_FAIL"):
    sys.stderr.write("srun: error: Unable to allocate resources: Invalid account\n")
    sys.exit(1)
os.environ["SLURM_JOB_ID"] = "4242"
os.execvp(args[0], args)
'''
FAKE_SACCTMGR = '#!/bin/sh\nprintf "%s\\n" "$FAKE_ASSOC" | tr ";" "\\n"\n'
LAB = "genome-center-grp|high|genome-center-grp-high-qos"
PUBLIC = "publicgrp|low|publicgrp-low-qos"


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)
    return path


class Harness(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        self.srun_log = os.path.join(self.d, "srun.log")
        _exe(os.path.join(self.bin, "srun"), FAKE_SRUN)
        _exe(os.path.join(self.bin, "sacctmgr"), FAKE_SACCTMGR)
        self.out = os.path.join(self.d, "ran.txt")

    def tearDown(self):
        self._tmp.cleanup()

    def login_node(self):
        _exe(os.path.join(self.bin, "sbatch"), "#!/bin/sh\nexit 0\n")

    def env(self, **extra):
        e = dict(os.environ)
        for k in ("SLURM_JOB_ID", "FAKE_SRUN_FAIL", "PROTEOMICS_ON_COMPUTE"):
            e.pop(k, None)
        # hermetic on a real cluster node too: no HPCCF cluster file, no CVMFS SLURM
        e.update({"PATH": self.bin + os.pathsep + e["PATH"], "FAKE_SRUN_LOG": self.srun_log,
                  "FAKE_ASSOC": LAB, "USER": "staff1",
                  "UCD_HPC_CLUSTER_CONF": os.path.join(self.d, "no-cluster-file"),
                  "HPCCF_SLURM_BIN": os.path.join(self.d, "no-slurm-bin")})
        e.update(extra)
        return e

    def srun_calls(self):
        if not os.path.exists(self.srun_log):
            return []
        with open(self.srun_log) as fh:
            return [json.loads(ln) for ln in fh if ln.strip()]

    def probe(self):
        """A command that records whether it ran inside a job."""
        return [sys.executable, "-c",
                "import os,sys; open(sys.argv[1],'w').write(os.environ.get('SLURM_JOB_ID','none'))",
                self.out]

    def ran_in(self):
        with open(self.out) as fh:
            return fh.read()


class OnComputeCommandLine(Harness):
    def run_oc(self, *args, cmd=None, **env):
        return subprocess.run([sys.executable, ON_COMPUTE, *args, "--", *(cmd or self.probe())],
                              capture_output=True, text=True, env=self.env(**env), timeout=60)

    def test_on_a_login_node_the_command_runs_under_srun_with_the_queue_spelled_out(self):
        self.login_node()
        r = self.run_oc("--time", "20", "--mem", "6G", "--cpus", "2", "--name", "pp_walk")
        self.assertEqual(r.returncode, 0, r.stderr)
        (call,) = self.srun_calls()
        opts = [a for a in call if a.startswith("--")]
        self.assertEqual(opts, ["--account=genome-center-grp", "--partition=high",
                                "--qos=genome-center-grp-high-qos", "--time=20", "--mem=6G",
                                "--ntasks=1", "--cpus-per-task=2", "--job-name=pp_walk"])
        self.assertEqual(call[len(opts):], self.probe())
        self.assertEqual(self.ran_in(), "4242")             # inside the (pretend) job
        self.assertIn("srun --account=genome-center-grp", r.stderr)   # said out loud

    def test_a_public_account_gets_low_and_its_qos(self):
        self.login_node()
        r = self.run_oc(FAKE_ASSOC=PUBLIC)
        self.assertEqual(r.returncode, 0, r.stderr)
        (call,) = self.srun_calls()
        self.assertIn("--account=publicgrp", call)
        self.assertIn("--partition=low", call)
        self.assertIn("--qos=publicgrp-low-qos", call)

    def test_defaults_are_explicit_time_memory_and_one_task(self):
        self.login_node()
        self.run_oc()
        (call,) = self.srun_calls()
        for flag in ("--time=30", "--mem=4G", "--ntasks=1", "--cpus-per-task=1"):
            self.assertIn(flag, call)

    def test_inside_a_job_it_runs_in_place(self):
        self.login_node()
        r = self.run_oc(SLURM_JOB_ID="77")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.srun_calls(), [])
        self.assertEqual(self.ran_in(), "77")

    def test_with_no_slurm_it_runs_in_place(self):
        r = self.run_oc()                                   # no sbatch on PATH: a laptop
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.srun_calls(), [])
        self.assertEqual(self.ran_in(), "none")

    def test_the_commands_exit_status_and_output_come_through(self):
        self.login_node()
        r = self.run_oc(cmd=[sys.executable, "-c", "import sys; print('hello'); sys.exit(5)"])
        self.assertEqual(r.returncode, 5)
        self.assertEqual(r.stdout.strip(), "hello")

    def test_print_shows_the_line_and_runs_nothing(self):
        self.login_node()
        r = self.run_oc("--print")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.srun_calls(), [])
        self.assertFalse(os.path.exists(self.out))
        line = shlex.split(r.stdout.strip())
        self.assertEqual(line[0], "srun")
        self.assertIn("--account=genome-center-grp", line)
        self.assertEqual(line[-len(self.probe()):], self.probe())

    def test_a_memory_size_needs_a_unit(self):
        """A bare number is megabytes to SLURM (ucd-hpc): refused, never guessed."""
        self.login_node()
        r = self.run_oc("--mem", "8")
        self.assertEqual(r.returncode, 2)
        self.assertIn("unit", r.stderr)
        self.assertEqual(self.srun_calls(), [])

    def test_no_command_is_a_usage_error(self):
        r = subprocess.run([sys.executable, ON_COMPUTE, "--time", "5"], capture_output=True,
                           text=True, env=self.env(), timeout=60)
        self.assertEqual(r.returncode, 2)

    def test_a_named_queue_is_used(self):
        self.login_node()
        self.run_oc("--partition", "low", "--account", "publicgrp", FAKE_ASSOC=f"{LAB};{PUBLIC}")
        (call,) = self.srun_calls()
        self.assertIn("--partition=low", call)
        self.assertIn("--account=publicgrp", call)

    def test_no_srun_beside_sbatch_is_said_not_a_traceback(self):
        self.login_node()
        os.remove(os.path.join(self.bin, "srun"))
        r = self.run_oc()
        self.assertEqual(r.returncode, 127)
        self.assertIn("cannot run", r.stderr)
        self.assertNotIn("Traceback", r.stderr)
        self.assertFalse(os.path.exists(self.out))

    def test_an_allocation_that_fails_is_a_failure(self):
        self.login_node()
        r = self.run_oc(FAKE_SRUN_FAIL="1")
        self.assertEqual(r.returncode, 1)
        self.assertIn("Invalid account", r.stderr)
        self.assertFalse(os.path.exists(self.out))


class SbatchHeader(Harness):
    """--header: the top of a job script the skill does not generate (a DE, set tests). The DE
    jobs of 2026-10 were written by hand each time; the header gives them the skill's queue rule
    and ucd-hpc's canonical lines."""

    def header(self, *args, **env):
        r = subprocess.run([sys.executable, ON_COMPUTE, "--header", os.path.join(self.d, "logs"),
                            "--name", "de_final", *args], capture_output=True, text=True,
                           env=self.env(**env), timeout=60)
        return r

    def test_every_request_is_explicit_and_the_log_is_per_job(self):
        r = self.header("--time", "120", "--mem", "48G", "--cpus", "8")
        self.assertEqual(r.returncode, 0, r.stderr)
        lines = r.stdout.splitlines()
        self.assertEqual(lines[0], "#!/bin/bash -l")
        for want in ("#SBATCH --account=genome-center-grp", "#SBATCH --partition=high",
                     "#SBATCH --qos=genome-center-grp-high-qos", "#SBATCH --time=120",
                     "#SBATCH --ntasks=1", "#SBATCH --cpus-per-task=8", "#SBATCH --mem=48G",
                     f'#SBATCH --output="{self.d}/logs/de_final_%j.log"', "set -euo pipefail"):
            self.assertIn(want, lines)
        self.assertNotIn("#SBATCH --requeue", lines)
        self.assertFalse([ln for ln in lines if "--exclusive" in ln])
        self.assertEqual(self.srun_calls(), [])

    def test_a_preemptible_queue_is_requeued(self):
        r = self.header(FAKE_ASSOC=PUBLIC)
        self.assertIn("#SBATCH --requeue", r.stdout.splitlines())
        self.assertIn("#SBATCH --partition=low", r.stdout.splitlines())

    def test_a_size_without_a_unit_is_refused(self):
        self.assertEqual(self.header("--mem", "48").returncode, 2)


class QueueNeverPreempts(Harness):
    """Review of this branch: run_search.slurm_queue sends work to publicgrp/low when the account
    is at its `high` cap -- fine for a requeued sbatch job, not for an srun, which a preemption
    cancels half way through an install or a copy. An srun now takes a queue that does not
    preempt; --preemptible lets the job rule choose."""

    def at_capacity(self):
        _exe(os.path.join(self.bin, "squeue"), "#!/bin/sh\necho 64\n")          # 64 CPUs running
        _exe(os.path.join(self.bin, "sinfo"), "#!/bin/sh\necho 0/4000/0/4000\n")  # low is idle
        return {"FAKE_ASSOC": f"{LAB};{PUBLIC}"}

    def run_oc(self, *args, **env):
        return subprocess.run([sys.executable, ON_COMPUTE, *args, "--", *self.probe()],
                              capture_output=True, text=True, env=self.env(**env), timeout=60)

    def test_the_job_rule_would_pick_low_here(self):
        self.login_node()
        env = self.env(**self.at_capacity())
        r = subprocess.run([sys.executable, "-c", "import run_search; "
                            "print(run_search.slurm_queue(peak_cpus=1, preemptible_ok=False))"],
                           cwd=SCRIPTS, capture_output=True, text=True, env=env, timeout=60)
        self.assertIn("'low'", r.stdout, r.stderr)

    def test_an_srun_stays_on_high_at_capacity(self):
        self.login_node()
        r = self.run_oc(**self.at_capacity())
        self.assertEqual(r.returncode, 0, r.stderr)
        (call,) = self.srun_calls()
        self.assertIn("--partition=high", call)
        self.assertIn("--account=genome-center-grp", call)
        self.assertNotIn("--requeue", r.stderr)

    def test_preemptible_lets_the_job_rule_choose(self):
        self.login_node()
        self.run_oc("--preemptible", **self.at_capacity())
        (call,) = self.srun_calls()
        self.assertIn("--partition=low", call)

    def test_a_free_tier_account_takes_its_high(self):
        self.login_node()
        self.run_oc(FAKE_ASSOC=f"publicgrp|high|publicgrp-high-qos;{PUBLIC}")
        (call,) = self.srun_calls()
        self.assertIn("--partition=high", call)
        self.assertIn("--qos=publicgrp-high-qos", call)

    def test_but_not_past_its_8_cpu_cap(self):
        self.login_node()
        self.run_oc("--cpus", "16", FAKE_ASSOC=f"publicgrp|high|publicgrp-high-qos;{PUBLIC}")
        (call,) = self.srun_calls()
        self.assertIn("--partition=low", call)

    def test_only_low_is_said(self):
        self.login_node()
        r = self.run_oc(FAKE_ASSOC=PUBLIC)
        self.assertIn("preemptible", r.stderr)


class LoginNodeWithoutSbatchOnPath(Harness):
    """Review of this branch: the login-node test needed sbatch on PATH. A non-login shell, or one
    after `module purge`, has none -- ucd-hpc: SLURM is at a fixed CVMFS path -- and every heavy
    step then ran on the login node. HPCCF's cluster file, or SLURM at that path, counts too."""

    def env_without_sbatch(self, **extra):
        e = self.env(**extra)
        e["PATH"] = os.pathsep.join(p for p in e["PATH"].split(os.pathsep) if p != self.bin)
        return e

    def run_oc(self, env):
        return subprocess.run([sys.executable, ON_COMPUTE, "--", *self.probe()],
                              capture_output=True, text=True, env=env, timeout=60)

    def test_slurm_at_the_fixed_path_is_a_login_node_and_its_srun_is_used(self):
        fixed = os.path.join(self.d, "cvmfs-slurm-bin")
        os.makedirs(fixed)
        _exe(os.path.join(fixed, "sbatch"), "#!/bin/sh\nexit 0\n")
        os.rename(os.path.join(self.bin, "srun"), os.path.join(fixed, "srun"))
        os.rename(os.path.join(self.bin, "sacctmgr"), os.path.join(fixed, "sacctmgr"))
        env = self.env_without_sbatch(HPCCF_SLURM_BIN=fixed)
        r = self.run_oc(env)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.ran_in(), "4242")
        (call,) = self.srun_calls()
        self.assertIn("--time=30", call)

    def test_hpccfs_cluster_file_alone_is_a_login_node(self):
        conf = os.path.join(self.d, "cluster_name.conf")
        with open(conf, "w") as fh:
            fh.write('export CLUSTER_NAME="hive"\n')
        env = self.env_without_sbatch(UCD_HPC_CLUSTER_CONF=conf)
        r = subprocess.run([sys.executable, ON_COMPUTE, "--is-login-node"], env=env, timeout=60)
        self.assertEqual(r.returncode, 0)
        env["SLURM_JOB_ID"] = "5"
        r = subprocess.run([sys.executable, ON_COMPUTE, "--is-login-node"], env=env, timeout=60)
        self.assertEqual(r.returncode, 1)

    def test_a_laptop_is_not(self):
        env = self.env_without_sbatch(UCD_HPC_CLUSTER_CONF=os.path.join(self.d, "nope"),
                                      HPCCF_SLURM_BIN=os.path.join(self.d, "nope"))
        r = subprocess.run([sys.executable, ON_COMPUTE, "--is-login-node"], env=env, timeout=60)
        self.assertEqual(r.returncode, 1)


# a script that moves itself: ensure_compute_node() first, then its work
MOVER = r'''import os, sys
sys.path.insert(0, {scripts!r})
import on_compute
on_compute.ensure_compute_node("this walks a share", minutes=15, mem="3G", script=__file__)
open(sys.argv[1], "w").write(os.environ.get("SLURM_JOB_ID", "none") + " " + " ".join(sys.argv[2:]))
'''


class EnsureComputeNode(Harness):
    def run_mover(self, *args, **env):
        script = os.path.join(self.d, "mover.py")
        with open(script, "w") as fh:
            fh.write(MOVER.format(scripts=SCRIPTS))
        return subprocess.run([sys.executable, script, self.out, *args], capture_output=True,
                              text=True, env=self.env(**env), timeout=60)

    def test_on_a_login_node_the_script_runs_again_under_srun_with_its_arguments(self):
        self.login_node()
        r = self.run_mover("--summary", "a b.json")
        self.assertEqual(r.returncode, 0, r.stderr)
        (call,) = self.srun_calls()
        self.assertIn("--time=15", call)
        self.assertIn("--mem=3G", call)
        self.assertEqual(self.ran_in(), "4242 --summary a b.json")
        self.assertIn("this walks a share", r.stderr)

    def test_whatever_runs_under_it_is_marked_and_never_sent_again(self):
        self.login_node()
        log = os.path.join(self.d, "env.log")
        self.run_mover("x", FAKE_SRUN_ENV_LOG=log)
        with open(log) as fh:
            self.assertEqual(fh.read().strip(), "1")
        os.remove(self.out)
        self.run_mover("y", PROTEOMICS_ON_COMPUTE="1")      # marked: runs here, no second srun
        self.assertEqual(len(self.srun_calls()), 1)
        self.assertEqual(self.ran_in(), "none y")

    def test_a_main_called_in_process_is_never_sent(self):
        """A test (or another script) that imports a script and calls its main() must not be
        re-run under srun: sys.argv[0] is then pytest, not the script."""
        self.login_node()
        r = subprocess.run([sys.executable, "-c",
                            "import sys; sys.path.insert(0, sys.argv[1]); import on_compute; "
                            "on_compute.ensure_compute_node('x', script=sys.argv[1] + '/some_script.py'); "
                            "print('carried on')", SCRIPTS],
                           capture_output=True, text=True, env=self.env(), timeout=60)
        self.assertEqual(r.stdout.strip(), "carried on", r.stderr)
        self.assertEqual(self.srun_calls(), [])

    def test_inside_a_job_or_without_slurm_it_carries_on_here(self):
        self.login_node()
        self.run_mover("x", SLURM_JOB_ID="9")
        self.assertEqual(self.ran_in(), "9 x")
        os.remove(os.path.join(self.bin, "sbatch"))
        self.run_mover("y")
        self.assertEqual(self.ran_in(), "none y")
        self.assertEqual(self.srun_calls(), [])


class LoginNodeIsOneDefinition(unittest.TestCase):
    def test_detect_acquisition_asks_on_compute(self):
        import detect_acquisition
        import on_compute
        from unittest import mock
        with mock.patch.object(on_compute, "on_login_node", return_value=True):
            self.assertTrue(detect_acquisition.on_cluster_login_node())
        with mock.patch.object(on_compute, "on_login_node", return_value=False):
            self.assertFalse(detect_acquisition.on_cluster_login_node())


class HiveExecJob(Harness):
    """hive_exec.sh --job '<cmd>': the command goes to HIVE wrapped in on_compute.py. The ssh
    stub runs what it is sent here, in a pretend HIVE home whose ~/proteomics-pipeline/scripts
    is this checkout."""

    def setUp(self):
        super().setUp()
        self.home = os.path.join(self.d, "home")
        os.makedirs(os.path.join(self.home, "proteomics-pipeline"))
        os.symlink(SCRIPTS, os.path.join(self.home, "proteomics-pipeline", "scripts"))
        # the remote side is a LOGIN shell (bash -l), whose system profile may rebuild PATH
        # (macOS path_helper): the pretend HIVE home puts the stubs and this python first again
        os.symlink(sys.executable, os.path.join(self.bin, "python3"))
        with open(os.path.join(self.home, ".bash_profile"), "w") as fh:
            fh.write(f'export PATH={shlex.quote(self.bin)}:"$PATH"\n')
        self.key = _exe(os.path.join(self.d, "id_test"), "not a key\n")
        self.ssh_log = os.path.join(self.d, "ssh.log")
        _exe(os.path.join(self.bin, "ssh"),
             '#!/bin/bash\nprintf "%s\\n" "${!#}" >> "$FAKE_SSH_LOG"\n'
             'HOME="$FAKE_HIVE_HOME" exec bash -c "${!#}"\n')

    def hive_exec(self, *args, **env):
        e = self.env(HIVE_USER="staff1", HIVE_KEY=self.key, HIVE_SSH_MUX="0",
                     FAKE_SSH_LOG=self.ssh_log, FAKE_HIVE_HOME=self.home,
                     HIVE_ENV_FILE=os.path.join(self.d, "none.env"), **env)
        return subprocess.run(["bash", HIVE_EXEC, *args], capture_output=True, text=True,
                              env=e, timeout=60)

    def test_job_runs_the_command_on_a_compute_node(self):
        self.login_node()
        cmd = (f"{shlex.quote(sys.executable)} -c \"import os,sys; "
               f"open(sys.argv[1],'w').write(os.environ.get('SLURM_JOB_ID','none'))\" "
               f"{shlex.quote(self.out)}")
        r = self.hive_exec("--job", "--time", "45", "--mem", "8G", cmd)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.ran_in(), "4242")
        (call,) = self.srun_calls()
        self.assertIn("--time=45", call)
        self.assertIn("--mem=8G", call)
        self.assertIn("--account=genome-center-grp", call)
        with open(self.ssh_log) as fh:
            self.assertIn("on_compute.py", fh.read())

    def test_job_keeps_a_quoted_argument_whole(self):
        self.login_node()
        cmd = (f"{shlex.quote(sys.executable)} -c \"import sys; "
               f"open(sys.argv[1],'w').write(sys.argv[2])\" {shlex.quote(self.out)} 'two words'")
        r = self.hive_exec("--job", cmd)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.ran_in(), "two words")

    def test_job_refuses_a_bad_size_before_connecting(self):
        r = self.hive_exec("--job", "--mem", "lots", "true")
        self.assertEqual(r.returncode, 2)
        self.assertFalse(os.path.exists(self.ssh_log))

    def test_job_needs_a_command(self):
        r = self.hive_exec("--job", "--time", "10")
        self.assertEqual(r.returncode, 2)


if __name__ == "__main__":
    unittest.main()
