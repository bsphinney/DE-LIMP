#!/usr/bin/env python3
"""
`watch_run.sh` says "step 1b/5 measuring" while step 1b runs, not "step 2/5" (S6).

The watcher read a 5-step chain's stage from its output folder, and called it step 2 as soon as
step 1's predicted library existed -- while step 1b (step1b_window.sbatch) was still measuring the
scan window on representative runs and no file had been searched. Step 1b writes
params.resolved.cfg as its last act (diann_parallel.py), so until that exists, in a chain that has
a step 1b, the stage is step 1b. A chain with no step 1b (a pinned window, a DDA cfg) goes from
step 1 to step 2 as before.

A fake sacct/squeue stands in for SLURM; nothing contacts a cluster.
"""
import json
import os
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
WATCH = os.path.join(os.path.dirname(HERE), "scripts", "watch_run.sh")


def chain_dir(o, step1b=True, resolved=False):
    with open(os.path.join(o, "file_list.txt"), "w") as fh:
        fh.write("".join(f"run{i}.d\n" for i in range(66)))
    with open(os.path.join(o, "step1.predicted.speclib"), "w") as fh:
        fh.write("lib")
    if step1b:
        with open(os.path.join(o, "step1b_window.sbatch"), "w") as fh:
            fh.write("#!/bin/bash\n")
    if resolved:
        with open(os.path.join(o, "params.resolved.cfg"), "w") as fh:
            fh.write("--window 7\n")


def watch(o, mode):
    with tempfile.TemporaryDirectory() as d:
        if mode == "chain":
            ids = ["301", "302", "303", "304", "305", "306"]
            with open(os.path.join(o, "jobs.txt"), "w") as fh:
                fh.write("\n".join(ids) + "\n")
            sacct = ('case "$*" in *"-j 301 "*) echo "  COMPLETED" ;; *"-j 302 "*) '
                     'echo "  RUNNING" ;; *) echo "  PENDING" ;; esac')
            args = ["bash", WATCH, "--all", o]
        else:
            sacct = 'echo "  RUNNING"'
            args = ["bash", WATCH, "--slurm", "302", "--out", o]
        for name, body in (("sacct", sacct), ("squeue", "")):
            p = os.path.join(d, name)
            with open(p, "w") as fh:
                fh.write(f"#!/bin/bash\n{body}\n")
            os.chmod(p, 0o755)
        env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
        r = subprocess.run(args, capture_output=True, text=True, env=env, timeout=60)
        assert r.returncode == 0, r.stderr
        return json.loads(r.stdout)


class Step1b(unittest.TestCase):
    def test_a_chain_dir_with_the_speclib_and_no_resolved_cfg_is_step_1b(self):
        for mode in ("slurm", "chain"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as o:
                chain_dir(o)
                out = watch(o, mode)
                self.assertEqual(out["progress"]["stage"], "step1b")
                self.assertIn("step 1b/5 measuring", out["progress"]["summary"])
                self.assertNotIn("step 2/5", out["progress"]["summary"])
                self.assertIn("Measuring", out.get("doing", ""))

    def test_once_the_resolved_cfg_exists_it_is_step_2(self):
        for mode in ("slurm", "chain"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as o:
                chain_dir(o, resolved=True)
                out = watch(o, mode)
                self.assertEqual(out["progress"]["stage"], "step2")
                self.assertIn("step 2/5", out["progress"]["summary"])

    def test_a_chain_with_no_step_1b_goes_to_step_2(self):
        """(pin) a pinned window or a DDA cfg: no step1b_window.sbatch, no resolved cfg."""
        with tempfile.TemporaryDirectory() as o:
            chain_dir(o, step1b=False)
            out = watch(o, "slurm")
            self.assertEqual(out["progress"]["stage"], "step2")


if __name__ == "__main__":
    unittest.main()
