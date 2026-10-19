#!/usr/bin/env python3
"""
watch_run.sh --hive with no usable python3 on the laptop (a staff report, skill 2.13.1, Git
Bash on Windows): `python3` there was the Microsoft Store stub, which prints "Python was not
found ..." and exits 9009. `--all <out> --hive` then returned {failed: true, done: true,
err_class: watcher_query_failed} for a healthy search -- a false failure an agent could act on
by resubmitting a running search. Without MSYS_NO_PATHCONV=1, Git Bash also rewrote the /nfs
path and the watcher reported no_jobs_file.

  * no usable local python3 + --hive  -> the whole poll runs ON HIVE (hive_exec.sh), whose
                                         python3 is real, and its JSON is printed as is;
  * HIVE cannot be reached that way   -> watcher_query_failed, failed false, done false;
  * no usable python3 and no --hive   -> watcher_no_python, failed false, done false;
  * a usable local python3            -> nothing changes: no poll is handed to HIVE.
With --hive, MSYS_NO_PATHCONV=1 is set for every HIVE path.

A fake hive_exec.sh stands in for HIVE (it runs the command here, with the real python3) and
a fake sacct/squeue for SLURM; nothing contacts a cluster.
"""
import json
import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
WATCH = os.path.join(SCRIPTS, "watch_run.sh")

STORE_STUB = ('#!/bin/bash\n'
              'echo "Python was not found; run without arguments to install from the Microsoft '
              'Store, or disable this shortcut from Settings > Apps > Advanced app settings > App '
              'execution aliases." >&2\nexit 49\n')
# "HIVE": logs what it was asked to run and the path-conversion setting it was called with,
# then runs it with the real PATH and the skill's scripts where ~/proteomics-pipeline is.
FAKE_HIVE_EXEC = r'''#!/bin/bash
printf 'MSYS_NO_PATHCONV=%s\n%s\n' "${MSYS_NO_PATHCONV:-}" "$1" >> "$HIVE_LOG"
[ -n "${HIVE_DOWN:-}" ] && { echo "ssh: connect to host hive.hpc.ucdavis.edu port 22: Connection timed out" >&2; exit 255; }
cmd="${1//\~\/proteomics-pipeline\/scripts/$REAL_SCRIPTS}"
PATH="$REAL_PATH" bash -c "$cmd"
'''


def _exe(path, body):
    with open(path, "w") as fh:
        fh.write(body)
    os.chmod(path, 0o755)


class NoLocalPython(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        # the laptop's copy of the skill: watch_run.sh beside the fake hive_exec.sh
        self.local = os.path.join(self.tmp, "laptop_scripts")
        os.makedirs(self.local)
        shutil.copy(WATCH, self.local)
        _exe(os.path.join(self.local, "hive_exec.sh"), FAKE_HIVE_EXEC)
        # SLURM on "HIVE"
        self.slurm = os.path.join(self.tmp, "slurm_bin")
        os.makedirs(self.slurm)
        _exe(os.path.join(self.slurm, "sacct"), '#!/bin/bash\ncase "$*" in *"-j 201 "*) '
             'echo "  COMPLETED" ;; *"-j 202 "*|*"-j 999 "*) echo "  RUNNING" ;; esac\n')
        _exe(os.path.join(self.slurm, "squeue"), "#!/bin/bash\n")
        # the laptop's python3: the Store stub
        self.stub = os.path.join(self.tmp, "stub_bin")
        os.makedirs(self.stub)
        _exe(os.path.join(self.stub, "python3"), STORE_STUB)
        self.log = os.path.join(self.tmp, "hive_exec.log")
        self.out = os.path.join(self.tmp, "search_out")
        os.makedirs(self.out)
        with open(os.path.join(self.out, "jobs.txt"), "w") as fh:
            fh.write("201\n202\n")
        with open(os.path.join(self.out, "diann_search_202.log"), "w") as fh:
            fh.write("[0:17] File #1/2\n")

    def watch(self, *args, python="stub", **env_extra):
        real = self.slurm + os.pathsep + os.environ.get("PATH", "")
        path = (self.stub + os.pathsep + real) if python == "stub" else real
        if python == "none":
            path = self.slurm + os.pathsep + "/usr/bin:/bin"
            if shutil.which("python3", path=path):
                path = self.slurm + os.pathsep + self.stub + os.pathsep + "/usr/bin:/bin"
        env = dict(os.environ, PATH=path, REAL_PATH=real, REAL_SCRIPTS=SCRIPTS,
                   HIVE_LOG=self.log, **env_extra)
        env.pop("MSYS_NO_PATHCONV", None)
        r = subprocess.run(["bash", os.path.join(self.local, "watch_run.sh"), *args],
                           capture_output=True, text=True, env=env, timeout=120)
        return r, json.loads(r.stdout)

    def hive_calls(self):
        if not os.path.exists(self.log):
            return ""
        with open(self.log) as fh:
            return fh.read()

    def test_a_chain_is_watched_on_hive_and_never_reads_as_failed(self):
        r, out = self.watch("--all", self.out, "--hive", "--poll", "3")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertFalse(out["failed"], out)
        self.assertFalse(out["done"], out)
        self.assertEqual(out["chain"], "201=COMPLETED 202=RUNNING")
        self.assertIn("DIA-NN is on file 1 of 2", out["progress"]["summary"])
        calls = self.hive_calls()
        self.assertIn("proteomics-pipeline/scripts/watch_run.sh --all", calls)
        self.assertNotIn("--hive", calls, "the poll on HIVE is a local one there")
        self.assertIn("MSYS_NO_PATHCONV=1", calls)

    def test_a_single_job_is_watched_on_hive(self):
        r, out = self.watch("--slurm", "999", "--hive")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out["state"], "RUNNING")
        self.assertFalse(out["failed"], out)
        self.assertIn("watch_run.sh --slurm 999", self.hive_calls())

    def test_hive_out_of_reach_is_a_query_failure_not_a_job_failure(self):
        r, out = self.watch("--all", self.out, "--hive", HIVE_DOWN="1")
        self.assertEqual(r.returncode, 3)
        self.assertEqual(out["error_class"], "watcher_query_failed")
        self.assertFalse(out["failed"], out)
        self.assertFalse(out["done"], out)
        self.assertIn("NOT a job failure", out["fix"])
        self.assertIn("Connection timed out", out["detail"])

    def test_hive_out_of_reach_with_a_local_python_is_not_failed_either(self):
        """The chain's own reachability check said failed: true, done: true -- the same false
        failure, from the same cause the fix calls "NOT a job failure"."""
        r, out = self.watch("--all", self.out, "--hive", python="real", HIVE_DOWN="1")
        self.assertEqual(r.returncode, 3)
        self.assertEqual(out["err_class"], "watcher_query_failed")
        self.assertFalse(out["failed"], out)
        self.assertFalse(out["done"], out)

    def test_no_python_and_no_hive_says_so_and_is_not_failed(self):
        r, out = self.watch("--slurm", "999", python="none")
        self.assertEqual(r.returncode, 3)
        self.assertEqual(out["error_class"], "watcher_no_python")
        self.assertFalse(out["failed"], out)
        self.assertFalse(out["done"], out)

    def test_a_usable_local_python_hands_nothing_to_hive_whole(self):
        r, out = self.watch("--all", self.out, "--hive", "--poll", "3", python="real")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(out["chain"], "201=COMPLETED 202=RUNNING")
        self.assertNotIn("watch_run.sh", self.hive_calls())
        self.assertIn("MSYS_NO_PATHCONV=1", self.hive_calls())


if __name__ == "__main__":
    unittest.main()
