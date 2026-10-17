#!/usr/bin/env python3
"""The skill's heavy non-search steps move themselves off a cluster login node.

HPC@UCD's ucd-hpc skill (mandatory on Hive since 2026-10-16): on a login node, nothing "past a few
seconds of CPU or much I/O"; "every directory that find, du, ls -R ... enters is a round trip to
a metadata server shared by the whole cluster". In hive_remote, SKILL.md ran these through
hive_exec.sh, i.e. on HIVE's login node:

  core_submission.py locate   lists every instrument's month folders on Flinders (list_raw)
  core_submission.py stage    lists the service tree's PI and institution folders, links raw files
  core_submission.py deliver  globs the whole service tree for the stage record
                              (find_service_record: ~4 min on a re-analysis) and copies
  fetch_fasta.py fetch/check-db  downloads, digests and reads whole FASTAs
  provenance.py               sha256 of every input raw file (up to 5 GB each)
  session.py finalize         Methods, Word, and a zip of the whole session
  checkpoint.py find          a recursive glob under --base (SKILL.md: "e.g. ~")

On a login node (sbatch on PATH, SLURM_JOB_ID unset) each now runs itself again under one
blocking srun (on_compute.ensure_compute_node). The light subcommands stay where they are.
The srun stub here only records the call (FAKE_SRUN_NORUN): what matters is that the work was
sent, with its own time and memory, and that a light command was not.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)

from test_on_compute import FAKE_SACCTMGR, FAKE_SRUN, LAB, _exe   # noqa: E402

HEAVY = [  # (script, args, --time, --mem)
    ("core_submission.py", ["locate", "--summary", "s.json", "--out", "o"], 30, "4G"),
    ("core_submission.py", ["stage", "--summary", "s.json", "--files", "f.txt"], 60, "4G"),
    ("core_submission.py", ["deliver", "--summary", "s.json", "--session", "S"], 60, "4G"),
    ("fetch_fasta.py", ["fetch", "--proteome", "UP000005640", "--out", "x.fasta"], 60, "8G"),
    ("fetch_fasta.py", ["check-db", "--search-dir", "S"], 30, "8G"),
    ("provenance.py", ["--outdir", "repro"], 60, "4G"),
    ("session.py", ["finalize", "--dir", "S", "--zip"], 60, "8G"),
    ("checkpoint.py", ["find", "--base", "."], 15, "2G"),
]
LIGHT = [
    ("core_submission.py", ["check"]),
    ("session.py", ["find-prior", "--raw", "nothing.raw"]),
    ("checkpoint.py", ["status", "--session", "nowhere"]),
    ("fetch_fasta.py", ["resolve", "--organism", "9606", "--offline"]),
]


class HeavyStepsMove(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        _exe(os.path.join(self.bin, "sbatch"), "#!/bin/sh\nexit 0\n")
        _exe(os.path.join(self.bin, "srun"), FAKE_SRUN)
        _exe(os.path.join(self.bin, "sacctmgr"), FAKE_SACCTMGR)
        self.log = os.path.join(self.d, "srun.log")

    def tearDown(self):
        self._tmp.cleanup()

    def run_script(self, script, args):
        env = dict(os.environ)
        env.pop("SLURM_JOB_ID", None)
        env.update({"PATH": self.bin + os.pathsep + env["PATH"], "FAKE_SRUN_LOG": self.log,
                    "FAKE_SRUN_NORUN": "1", "FAKE_ASSOC": LAB, "USER": "staff1",
                    "HOME": self.d, "CORE_FLINDERS_ROOT": self.d, "CORE_WORK_ROOT": self.d})
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args], cwd=self.d,
                              capture_output=True, text=True, env=env, timeout=120)

    def calls(self):
        if not os.path.exists(self.log):
            return []
        with open(self.log) as fh:
            return [json.loads(ln) for ln in fh if ln.strip()]

    def test_each_heavy_step_runs_itself_on_a_compute_node(self):
        for script, args, minutes, mem in HEAVY:
            with self.subTest(script=script, cmd=args[0]):
                if os.path.exists(self.log):
                    os.remove(self.log)
                r = self.run_script(script, args)
                self.assertEqual(r.returncode, 0, r.stderr[-1500:])
                (call,) = self.calls()
                self.assertIn(f"--time={minutes}", call)
                self.assertIn(f"--mem={mem}", call)
                self.assertIn("--account=genome-center-grp", call)
                i = call.index(os.path.join(SCRIPTS, script))
                self.assertEqual(call[i + 1:], args)          # the same command, unchanged
                self.assertIn("compute node", r.stderr)

    def test_light_steps_stay_where_they_are(self):
        for script, args in LIGHT:
            with self.subTest(script=script, cmd=args[0]):
                self.run_script(script, args)
                self.assertEqual(self.calls(), [])


class CheckpointFindIsBounded(unittest.TestCase):
    """checkpoint.py find walks --base on a compute node now, and still never everything: it
    stops at --max-depth (default 6) and never enters a search's per-run or raw-data folders."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.base = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def plant(self, *parts):
        d = os.path.join(self.base, *parts)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, ".recovery.json"), "w") as fh:
            json.dump({"updated": "2026-10-16", "stages": [{"stage": "search"}]}, fh)
        return d

    def find(self, *extra):
        env = {k: v for k, v in os.environ.items() if k != "SLURM_JOB_ID"}
        env["SLURM_JOB_ID"] = "1"                         # as if in the job already
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "checkpoint.py"), "find",
                            "--base", self.base, *extra], capture_output=True, text=True,
                           env=env, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        return sorted(h["session"] for h in json.loads(r.stdout)["sessions"])

    def test_sessions_within_the_depth_are_found(self):
        a = self.plant("PROT_0807", "sessions", "s1")
        b = self.plant("lab", "PROT_9001", "sessions", "s2")
        self.assertEqual(self.find(), sorted([a, b]))

    def test_deeper_than_max_depth_is_not_walked(self):
        self.plant(*[f"d{i}" for i in range(8)])
        self.assertEqual(self.find(), [])
        self.assertEqual(len(self.find("--max-depth", "9")), 1)

    def test_a_searchs_per_run_folders_are_not_entered(self):
        self.plant("S", "output", "search", "quant_step2", "x")
        self.plant("S", "output", "search", "quant", "x")
        self.plant("S", "output", "search", "report_xic", "x")
        self.plant("S", "raw", "run1.d", "y")
        self.assertEqual(self.find(), [])

    def test_a_session_beside_the_raw_data_is_found(self):
        """session.py init's default puts the session beside the raw files (<raw dir>/<date>_<name>),
        and raw folders are often called raw/: never skipped (review of this branch: the first
        version did, so a resume could have resubmitted a running search)."""
        a = self.plant("proj", "raw", "2026-10-16_hela")
        b = self.plant("proj", "Quant_study", "sessions", "s1")
        c = self.plant("proj", "quantification", "s2")
        self.assertEqual(self.find(), sorted([a, b, c]))

    def test_a_skipped_folder_that_is_itself_a_session_is_found(self):
        a = self.plant("proj", "temp")
        self.assertEqual(self.find(), [a])

    def test_a_linked_folder_is_followed_once_and_a_loop_ends(self):
        elsewhere = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        self.addCleanup(lambda: __import__("shutil").rmtree(elsewhere, ignore_errors=True))
        os.makedirs(os.path.join(elsewhere, "sessions", "s3"))
        with open(os.path.join(elsewhere, "sessions", "s3", ".recovery.json"), "w") as fh:
            json.dump({"updated": "x", "stages": []}, fh)
        os.makedirs(os.path.join(self.base, "proj"))
        os.symlink(elsewhere, os.path.join(self.base, "proj", "linked"))
        os.symlink(self.base, os.path.join(self.base, "proj", "loop"))     # back to the top
        found = self.find()
        self.assertEqual(found, [os.path.join(self.base, "proj", "linked", "sessions", "s3")])


if __name__ == "__main__":
    unittest.main()
