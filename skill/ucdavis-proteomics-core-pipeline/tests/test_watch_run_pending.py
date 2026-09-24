#!/usr/bin/env python3
"""
watch_run.sh's `progress.summary` and `doing` are what the orchestrator hands the user
verbatim, and `error_class` is what it acts on. What a PENDING job means depends on why:

  * nothing has run yet            -> "queued, not started" (a held search was reported as
                                      "search running" / "Searching your files" -- HIVE e2e
                                      test 2026-09-23);
  * waiting on an earlier step     -> the chain's step-5 job is PENDING (Dependency) for most
                                      of a run while steps 2-4 finish files: the file count
                                      must survive (review 2026-09-23: it was overwritten
                                      with "not started");
  * part of an array already done  -> "partly done", not "not started";
  * held (JobHeldUser/Admin)       -> error_class "held": it never starts by itself, and it is
                                      NOT failed (a failed job gets resubmitted = a duplicate).

A fake sacct/squeue stands in for SLURM; nothing contacts a cluster.
"""
import json
import os
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
WATCH = os.path.join(os.path.dirname(HERE), "scripts", "watch_run.sh")


class PendingIsNotRunning(unittest.TestCase):
    def watch(self, sacct_lines, reason="", outdir=None):
        with tempfile.TemporaryDirectory() as d:
            sacct = "\n".join(f'echo "  {ln}"' for ln in sacct_lines)
            squeue = f'case "$*" in *%R*) echo "{reason}" ;; esac'
            for name, body in (("sacct", sacct), ("squeue", squeue)):
                p = os.path.join(d, name)
                with open(p, "w") as fh:
                    fh.write(f"#!/bin/bash\n{body}\n")
                os.chmod(p, 0o755)
            env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
            args = ["bash", WATCH, "--slurm", "999"] + (["--out", outdir] if outdir else [])
            r = subprocess.run(args, capture_output=True, text=True, env=env, timeout=60)
            self.assertEqual(r.returncode, 0, r.stderr)
            return json.loads(r.stdout)

    def test_nothing_run_yet_says_queued(self):
        out = self.watch(["PENDING"], reason="Priority")
        self.assertEqual(out["state"], "PENDING")
        self.assertIn("queued, not started", out["progress"]["summary"])
        self.assertIn("not started", out["doing"])
        self.assertFalse(out["failed"])

    def test_running_is_unchanged(self):
        out = self.watch(["RUNNING"])
        self.assertEqual(out["progress"]["summary"], "search running")
        self.assertNotIn("queue", out.get("doing", ""))

    def test_chain_step5_pending_keeps_the_file_count(self):
        with tempfile.TemporaryDirectory() as o:
            with open(os.path.join(o, "file_list.txt"), "w") as fh:
                fh.write("".join(f"run{i}.raw\n" for i in range(66)))
            os.makedirs(os.path.join(o, "quant_step2"))
            for i in range(34):
                open(os.path.join(o, "quant_step2", f"run{i}.quant"), "w").close()
            with open(os.path.join(o, "step1.predicted.speclib"), "w") as fh:
                fh.write("lib")
            out = self.watch(["PENDING"], reason="Dependency", outdir=o)
        self.assertIn("34/66", out["progress"]["summary"])
        self.assertNotIn("not started", out["progress"]["summary"])
        self.assertNotIn("not started", out.get("doing", ""))

    def test_partly_done_array_is_not_not_started(self):
        out = self.watch(["3 COMPLETED", "5 PENDING"], reason="QOSMaxCpuPerUserLimit")
        self.assertEqual(out["state"], "PENDING")
        self.assertIn("partly done", out["progress"]["summary"])
        self.assertNotIn("not started", out.get("doing", ""))

    def test_waiting_on_an_earlier_step_says_so(self):
        out = self.watch(["PENDING"], reason="Dependency")
        self.assertIn("waiting for earlier steps", out["progress"]["summary"])

    def test_held_job_is_actionable_but_not_failed(self):
        out = self.watch(["PENDING"], reason="JobHeldUser")
        self.assertEqual(out["error_class"], "held")
        self.assertIn("scontrol release", out["fix"])
        self.assertFalse(out["failed"], "a held job must not be resubmitted as a failure")
        self.assertIn("HELD", out["progress"]["summary"])


if __name__ == "__main__":
    unittest.main()
