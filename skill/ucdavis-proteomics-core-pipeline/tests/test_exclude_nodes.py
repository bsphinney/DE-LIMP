#!/usr/bin/env python3
"""
Keeping a search off nodes known to be bad: `run_search.py --exclude <node[,node...]>`.

2026-10-05: two projects had to keep their jobs off HIVE nodes with storage faults (one job hung
67 min at 'Loading run'; two were stuck 'Kill task failed'). run_search.py had no --exclude, the
generated job scripts took none, and `export SBATCH_EXCLUDE=...` did nothing (ExcNodeList=(null)):
sbatch reads no such variable (man sbatch, INPUT ENVIRONMENT VARIABLES). The way round was
`scontrol update JobId=<id> ExcNodeList=...` on 64 pending jobs by hand.

These pin:
  * the value is checked as SLURM node names / ranges before anything is written -- it lands on a
    #SBATCH line of every job;
  * `#SBATCH --exclude=` in EVERY generated script: the single-shot jobs, the 5-step chain (also
    when run_search.py routes there), the Radiant chain, the one-job routes (Sage, ...);
  * it is recorded in search_provenance.json, and a Sage re-run command carries it;
  * node_fault.py retry passes the union of the script's nodes and the failed ones (a command-line
    --exclude replaces the script's -- tests/test_node_fault.py).
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
from job_env import job_env  # noqa: E402

import node_fault  # noqa: E402
import run_search  # noqa: E402
import test_single_shot_sbatch as sss  # noqa: E402  (the single-shot workspace)


def _header(path):
    """The #SBATCH lines SLURM reads: those before the first command."""
    lines = []
    with open(path) as fh:
        for ln in fh:
            s = ln.strip()
            if s and not s.startswith("#"):
                break
            if s.startswith("#SBATCH"):
                lines.append(s)
    return lines


def _excludes(path):
    return [ln for ln in _header(path) if re.match(r"#SBATCH\s+(--exclude|-x)\b", ln)]


class ParseTests(unittest.TestCase):
    def test_names_and_ranges_become_one_sorted_list_without_repeats(self):
        self.assertEqual(node_fault.exclude_nodes("hive-dc-7-4-54,hive-dc-7-4-[50-51],hive-dc-7-4-54"),
                         ["hive-dc-7-4-50", "hive-dc-7-4-51", "hive-dc-7-4-54"])
        self.assertEqual(node_fault.exclude_nodes("cpu[01-02,09]"), ["cpu01", "cpu02", "cpu09"])
        self.assertEqual(node_fault.exclude_nodes(" n1 , n2 "), ["n1", "n2"])

    def test_anything_that_is_not_a_node_list_is_refused(self):
        for bad in ("", ",", "n1,,n2", "n1 n2", "n1;rm -rf x", "$(touch x)", "`x`", "n'1",
                    "/etc/bad_nodes", "-x", "--exclude=n1", "n[3-1]", "n[1-", "n[a-b]", "n]1[",
                    "n1\nn2", "n[1-100000]"):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    node_fault.exclude_nodes(bad)

    def test_the_directive(self):
        self.assertEqual(node_fault.exclude_directive(["n1", "n2"]), "#SBATCH --exclude=n1,n2")
        self.assertIsNone(node_fault.exclude_directive([]))


class SingleShot(unittest.TestCase):
    def test_both_jobs_carry_it_and_the_record_says_so(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d)
            w.generate("--exclude", "hive-dc-7-4-[50-51],hive-as-11-3-62")
            want = "#SBATCH --exclude=hive-as-11-3-62,hive-dc-7-4-50,hive-dc-7-4-51"
            for job in (w.lib_job, w.search_job):
                self.assertEqual(_excludes(job), [want], os.path.basename(job))
            with open(os.path.join(w.out, "search_provenance.json")) as fh:
                prov = json.load(fh)
            self.assertEqual(prov["queue"]["exclude"],
                             ["hive-as-11-3-62", "hive-dc-7-4-50", "hive-dc-7-4-51"])

    def test_without_it_no_job_excludes_anything(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d)
            w.generate()
            for job in (w.lib_job, w.search_job):
                self.assertEqual(_excludes(job), [], os.path.basename(job))
            with open(os.path.join(w.out, "search_provenance.json")) as fh:
                self.assertIsNone(json.load(fh)["queue"].get("exclude"))

    def test_a_bad_value_stops_before_anything_is_written(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d)
            p = w.generate("--exclude", "n1;touch pwned", check=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("--exclude", p.stderr)
            self.assertFalse(os.path.exists(w.lib_job))
            self.assertFalse(os.path.exists(os.path.join(w.out, "submit.sh")))


class Chains(unittest.TestCase):
    def test_every_step_of_the_five_step_chain(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d, n_runs=3)
            out = os.path.join(d, "out")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", w.diann, "--raw", *w.runs, "--fasta", w.fasta,
                                "--out", out, "--cfg", w.cfg, "--partition", "low",
                                "--account", "publicgrp", "--exclude", "nB,nA"],
                               capture_output=True, text=True, env=job_env(d), timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            scripts = sorted(n for n in os.listdir(out) if n.endswith(".sbatch"))
            self.assertGreaterEqual(len(scripts), 6, scripts)
            for s in scripts:
                self.assertEqual(_excludes(os.path.join(out, s)), ["#SBATCH --exclude=nA,nB"], s)
            self.assertEqual(json.loads(g.stdout)["exclude"], ["nA", "nB"])

    def test_run_search_hands_it_to_the_chain(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d, n_runs=6)
            bindir = os.path.join(d, "bin")
            os.makedirs(bindir)
            sss._write(os.path.join(bindir, "sbatch"), "#!/bin/sh\necho 1\n", 0o755)
            env = job_env(d, PATH=bindir + os.pathsep + "/usr/bin:/bin")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", w.tools, "--bundle", w.bundle, "--params", w.cfg,
                                "--fasta", w.fasta, "--out", w.out, "--files", *w.runs,
                                "--partition", "low", "--account", "publicgrp",
                                "--exclude", "hive-dc-7-4-50"],
                               capture_output=True, text=True, env=env, timeout=180)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            with open(os.path.join(w.out, "search_provenance.json")) as fh:
                prov = json.load(fh)
            self.assertEqual(prov["search_mode"], "parallel_5step")
            self.assertEqual(prov["queue"]["exclude"], ["hive-dc-7-4-50"])
            scripts = [n for n in os.listdir(w.out) if n.endswith(".sbatch")]
            self.assertTrue(scripts)
            for s in scripts:
                self.assertEqual(_excludes(os.path.join(w.out, s)),
                                 ["#SBATCH --exclude=hive-dc-7-4-50"], s)

    def test_every_step_of_the_radiant_chain(self):
        with tempfile.TemporaryDirectory() as d:
            mzml = [os.path.join(d, f"r{i}.mzML") for i in range(2)]
            for f in mzml:
                open(f, "w").close()
            listing = sss._write(os.path.join(d, "files.txt"), "\n".join(mzml) + "\n")
            fasta = sss._write(os.path.join(d, "db.fasta"), "")
            config = sss._write(os.path.join(d, "x.radiantConfig"), "")
            out = os.path.join(d, "out")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "radiant_parallel.py"),
                                "--runtime", "apptainer", "--image", "/x/radiant.sif",
                                "--raw-list", listing, "--fasta", fasta, "--config", config,
                                "--out", out, "--diann", "/opt/diann-linux", "--partition", "low",
                                "--account", "publicgrp", "--exclude", "n[1-2]"],
                               capture_output=True, text=True, env=job_env(d), timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            D = os.path.join(out, "radiant_parallel")
            scripts = sorted(n for n in os.listdir(D) if n.endswith(".sbatch"))
            self.assertEqual(len(scripts), 3, scripts)
            for s in scripts:
                self.assertEqual(_excludes(os.path.join(D, s)), ["#SBATCH --exclude=n1,n2"], s)

    def test_a_chain_generator_refuses_a_bad_value_too(self):
        with tempfile.TemporaryDirectory() as d:
            w = sss._Workspace(d, n_runs=3)
            out = os.path.join(d, "out")
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", w.diann, "--raw", *w.runs, "--fasta", w.fasta,
                                "--out", out, "--cfg", w.cfg, "--partition", "low",
                                "--account", "publicgrp", "--exclude", "/tmp/nodes.txt"],
                               capture_output=True, text=True, env=job_env(d), timeout=120)
            self.assertNotEqual(g.returncode, 0)
            self.assertIn("--exclude", g.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "submit.sh")))


class OneJobRoutes(unittest.TestCase):
    def test_emit_sbatch_writes_it_into_the_header(self):
        """Sage, AlphaDIA, FragPipe and single-job Radiant all write their job with
        emit_sbatch(..., **queue); the queue carries the nodes."""
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            run_search.emit_sbatch(job, "echo search", d, 4, job="sage_search",
                                   partition="low", account="publicgrp",
                                   qos="publicgrp-low-qos", exclude=["nA", "nB"], notify=None)
            self.assertEqual(_excludes(job), ["#SBATCH --exclude=nA,nB"])
            run_search.emit_sbatch(job, "echo search", d, 4, job="sage_search",
                                   partition="low", account="publicgrp", notify=None)
            self.assertEqual(_excludes(job), [])

    def test_the_job_end_hook_keeps_it_in_the_header(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            run_search.emit_sbatch(job, "echo search", d, 4, job="sage_search",
                                   partition="low", account="publicgrp", exclude=["nA"])
            self.assertEqual(_excludes(job), ["#SBATCH --exclude=nA"])


class SageRerunCarriesIt(unittest.TestCase):
    def test_the_corrected_window_rerun_keeps_the_nodes_out(self):
        import sage_lfq_check
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "sage_out")
            os.makedirs(out)
            with open(os.path.join(out, "search_provenance.json"), "w") as fh:
                json.dump({"tools": "/t/tools.json", "bundle": "/t/b.json", "fasta": "/t/db.fasta",
                           "files": ["/t/a.mzML"], "threads": 8,
                           "queue": {"partition": "low", "account": "publicgrp", "qos": None,
                                     "exclude": ["nA", "nB"]}}, fh)
            plan = sage_lfq_check.rerun_plan(out, "/t/corrected.json", 20)
            self.assertIn("--exclude nA,nB", plan["command"])


if __name__ == "__main__":
    unittest.main()
