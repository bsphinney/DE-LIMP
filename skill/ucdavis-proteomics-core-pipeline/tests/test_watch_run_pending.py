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



class ChainModeSaysWhatItIsDoing(unittest.TestCase):
    """2026-10: `watch_run.sh --all <out>` on a two-job search printed only
    {"chain": "101=COMPLETED 102=RUNNING"}: a 66-minute search had nothing to relay
    but RUNNING, though SKILL.md says to relay progress.summary, doing and note on every poll."""

    def watch(self, out, states):
        with tempfile.TemporaryDirectory() as d:
            body = 'case "$*" in ' + " ".join(f'*"-j {j} "*) echo "  {s}" ;;'
                                              for j, s in states.items()) + " esac"
            for name, text in (("sacct", body), ("squeue", "")):
                p = os.path.join(d, name)
                with open(p, "w") as fh:
                    fh.write(f"#!/bin/bash\n{text}\n")
                os.chmod(p, 0o755)
            with open(os.path.join(out, "jobs.txt"), "w") as fh:
                fh.write("".join(f"{j}\n" for j in states))
            env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
            r = subprocess.run(["bash", WATCH, "--all", out, "--poll", "3"], capture_output=True,
                               text=True, env=env, timeout=60)
            self.assertEqual(r.returncode, 0, r.stderr)
            return json.loads(r.stdout)

    def test_a_two_job_search_reports_the_file_diann_is_on(self):
        with tempfile.TemporaryDirectory() as o:
            with open(os.path.join(o, "file_list.txt"), "w") as fh:
                fh.write("".join(f"run{i}.d\n" for i in range(4)))
            with open(os.path.join(o, "diann_search_202.log"), "w") as fh:
                fh.write("[0:17] File #1/4\n[0:40] File #2/4\n[1:02] File #3/4\n"
                         "[1:05] Loading run /x/run2.d\n")
            out = self.watch(o, {"201": "COMPLETED", "202": "RUNNING"})
        self.assertFalse(out["done"])
        self.assertEqual(out["chain"], "201=COMPLETED 202=RUNNING")
        self.assertEqual((out["progress"]["file"], out["progress"]["files"]), (3, 4))
        self.assertIn("DIA-NN is on file 3 of 4", out["progress"]["summary"])
        self.assertIn("diann_search_202.log", out["progress"]["summary"])
        self.assertTrue(out["doing"])
        self.assertTrue(out["note"])

    def test_an_array_jobs_task_log_is_read(self):
        """Review of 2.13.1 (note 7): an array task logs to <job>_<id>_<task>.log."""
        with tempfile.TemporaryDirectory() as o:
            with open(os.path.join(o, "search_202_3.log"), "w") as fh:
                fh.write("[0:17] File #2/6\n")
            out = self.watch(o, {"201": "COMPLETED", "202": "RUNNING"})
        self.assertEqual((out["progress"]["file"], out["progress"]["files"]), (2, 6))

    def test_a_five_step_chain_counts_finished_files(self):
        with tempfile.TemporaryDirectory() as o:
            with open(os.path.join(o, "file_list.txt"), "w") as fh:
                fh.write("".join(f"run{i}.d\n" for i in range(4)))
            open(os.path.join(o, "step1.predicted.speclib"), "w").write("x")
            os.makedirs(os.path.join(o, "quant_step2"))
            for i in range(2):
                open(os.path.join(o, "quant_step2", f"run{i}.quant"), "w").write("x")
            out = self.watch(o, {"301": "COMPLETED", "302": "RUNNING", "303": "PENDING",
                                 "304": "PENDING", "305": "PENDING"})
        self.assertEqual(out["progress"]["summary"], "step 2/5: 2/4 files done (50%)")
        self.assertEqual(out["progress"]["stage"], "step2")
        self.assertTrue(out["doing"].startswith("First pass"))

    def test_finished_and_failed_chains_are_as_before(self):
        with tempfile.TemporaryDirectory() as o:
            out = self.watch(o, {"401": "COMPLETED", "402": "COMPLETED"})
        self.assertTrue(out["done"])
        self.assertNotIn("progress", out)


if __name__ == "__main__":
    unittest.main()
