#!/usr/bin/env python3
"""
CPUs per array task in the DIA-NN 5-step chain (run_search.array_task_cpus).

genome-center-grp-high-qos caps a user at 64 CPUs (MaxTRESPU cpu=64,mem=1T; `sacctmgr show qos`
on HIVE, 2026-10-01). A chain generated with --threads 32 asked 32 CPUs per array task, so 2 of
29 files ran at a time and the rest sat PENDING with QOSMaxCpuPerUserLimit (2026-09-25).
Under a per-user cap, fewer CPUs per task never lowers throughput while files are waiting, so the
array tasks are now sized to the cap; a queue without one (publicgrp/low) keeps what was asked,
and --threads-per-file pins it. Everything SLURM is stubbed; PATH holds only the stub directory.
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
import run_search  # noqa: E402

GC_HIGH = {"qos": "genome-center-grp-high-qos", "cpu": 64, "mem_gb": 1024, "per_job_cpu": None,
           "per_job_mem_gb": None, "source": "test"}
PUB_LOW = {"qos": "publicgrp-low-qos", "cpu": None, "mem_gb": None, "per_job_cpu": None,
           "per_job_mem_gb": None, "source": "test"}
PUB_HIGH = {"qos": "publicgrp-high-qos", "cpu": None, "mem_gb": None, "per_job_cpu": 8,
            "per_job_mem_gb": 128, "source": "test"}

# What sacctmgr says on HIVE (2026-10-01), for the commands run_search.py runs.
SACCTMGR = """#!/bin/sh
case "$*" in
  *"show assoc"*)
    echo 'genome-center-grp|high|genome-center-grp-high-qos'
    echo 'publicgrp|low|publicgrp-low-qos' ;;
  *"show qos"*genome-center-grp-high-qos*)
    echo 'genome-center-grp-high-qos|cpu=64,gres/gpu=0,mem=1T|' ;;
  *"show qos"*publicgrp-low-qos*)
    echo 'publicgrp-low-qos||' ;;
esac
"""


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)


def _read(path):
    with open(path) as fh:
        return fh.read()


def size(n, ceiling, limits, pinned=None, mem=64, max_sim=20):
    return run_search.array_task_cpus(n, ceiling, limits=limits, mem_per_task_gb=mem,
                                      max_simultaneous=max_sim, pinned=pinned)


class SizingRuleTests(unittest.TestCase):
    def test_a_29_file_chain_runs_eight_files_at_once_not_two(self):
        r = size(29, 32, GC_HIGH)
        self.assertEqual((r["cpus"], r["concurrent"]), (8, 8))
        self.assertEqual(r["concurrent_at_requested"], 2)
        self.assertIn("32 per file would run 2 at once", r["summary"])

    def test_few_files_get_more_cpus_each_all_running_at_once(self):
        r = size(6, 32, GC_HIGH)
        self.assertEqual((r["cpus"], r["concurrent"]), (10, 6))

    def test_never_more_than_was_asked(self):
        self.assertEqual(size(6, 8, GC_HIGH)["cpus"], 8)
        self.assertEqual(size(29, 4, GC_HIGH)["cpus"], 4)      # below the floor: as asked

    def test_no_per_user_cap_keeps_what_was_asked(self):
        r = size(29, 32, PUB_LOW)
        self.assertEqual((r["cpus"], r["concurrent"]), (32, 20))   # --max-simultaneous
        self.assertIn("no per-user CPU cap", r["reason"])

    def test_a_pin_is_honoured_and_its_cost_said(self):
        r = size(29, 8, GC_HIGH, pinned=32)
        self.assertEqual((r["cpus"], r["concurrent"]), (32, 2))
        self.assertTrue(r["pinned"])
        self.assertIn("pinned", r["reason"])

    def test_a_pin_above_the_cap_is_lowered_or_it_would_never_start(self):
        r = size(29, 8, GC_HIGH, pinned=128)
        self.assertEqual(r["cpus"], 64)
        self.assertIn("lowered to 64", r["reason"])

    def test_a_per_job_cap_is_respected(self):
        self.assertEqual(size(29, 16, PUB_HIGH)["cpus"], 8)

    def test_the_memory_cap_limits_concurrency_too(self):
        r = size(29, 32, GC_HIGH, mem=256)
        self.assertEqual(r["concurrent"], 4)                     # 1024 GB / 256 GB

    def test_single_jobs_are_fitted_to_the_queue(self):
        self.assertEqual(run_search.fit_to_queue(64, GC_HIGH), 64)
        self.assertEqual(run_search.fit_to_queue(96, GC_HIGH), 64)
        self.assertEqual(run_search.fit_to_queue(64, PUB_HIGH), 8)
        self.assertEqual(run_search.fit_to_queue(64, PUB_LOW), 64)

    def test_tres_strings(self):
        self.assertEqual(run_search._tres("cpu=64,gres/gpu=0,mem=1T", "cpu"), 64)
        self.assertEqual(run_search._tres("cpu=64,gres/gpu=0,mem=1T", "mem"), 1024)
        self.assertEqual(run_search._tres("cpu=8,gres/gpu=1,mem=128G", "mem"), 128)
        self.assertEqual(run_search._tres("mem=512000M", "mem"), 500)
        self.assertEqual(run_search._tres("mem=2048", "mem"), 2)          # MB when no unit
        self.assertIsNone(run_search._tres("", "cpu"))
        self.assertIsNone(run_search._tres("gres/gpu=0", "cpu"))


class StubbedSlurm(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        _exe(os.path.join(self.bin, "sinfo"), "#!/bin/sh\necho '0/5000/0/5000'\n")
        _exe(os.path.join(self.bin, "squeue"), "#!/bin/sh\necho 0\n")
        self.env = {"PATH": self.bin, "HOME": self.d, "USER": "someone"}

    def tearDown(self):
        self._tmp.cleanup()

    def limits(self, *queue, no_sacctmgr=False):
        # no_sacctmgr: also hide the install locations _sacctmgr_path() falls back to -- on HIVE
        # /usr/bin/sacctmgr exists, and a PATH of stubs alone does not hide it
        code = ("import json, sys; sys.path.insert(0, sys.argv[1]); import run_search; "
                + ("run_search._sacctmgr_path = lambda: None; " if no_sacctmgr else "")
                + "print(json.dumps(run_search.user_limits(*sys.argv[2:])))")
        r = subprocess.run([sys.executable, "-c", code, SCRIPTS, *queue], env=self.env,
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        return json.loads(r.stdout)


class UserLimitsTests(StubbedSlurm):
    def test_read_from_the_qos_of_your_association(self):
        _exe(os.path.join(self.bin, "sacctmgr"), SACCTMGR)
        lim = self.limits("high", "genome-center-grp")
        self.assertEqual((lim["qos"], lim["cpu"], lim["mem_gb"]),
                         ("genome-center-grp-high-qos", 64, 1024))
        self.assertIn("sacctmgr", lim["source"])
        lim = self.limits("low", "publicgrp", "publicgrp-low-qos")
        self.assertEqual((lim["cpu"], lim["mem_gb"]), (None, None))

    def test_without_sacctmgr_only_the_known_cap_is_assumed(self):
        lim = self.limits("high", "genome-center-grp", no_sacctmgr=True)
        self.assertEqual(lim["cpu"], run_search.HIVE_USER_CPU_CAP)
        self.assertIn("HIVE_USER_CPU_CAP", lim["source"])
        lim = self.limits("low", "publicgrp", "publicgrp-low-qos", no_sacctmgr=True)
        self.assertIsNone(lim["cpu"])
        self.assertIn("unknown", lim["source"])


class ChainHeaderTests(StubbedSlurm):
    def setUp(self):
        super().setUp()
        _exe(os.path.join(self.bin, "sacctmgr"), SACCTMGR)

    def chain(self, *extra, n=29, window=True):
        raws = []
        for i in range(n):
            raws.append(os.path.join(self.d, f"f{i:02d}.d"))
            os.makedirs(raws[-1], exist_ok=True)
        cfg = os.path.join(self.d, "diann.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n"
                     + ("--window 7\n" if window else ""))
        fasta = os.path.join(self.d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEK\n")
        out = os.path.join(self.d, "out")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", "/opt/diann-2.7.0/diann-linux", "--raw", *raws,
                            "--fasta", fasta, "--out", out, "--cfg", cfg,
                            *([] if not window else ["--no-probe-window"]), *extra],
                           cwd=self.d, capture_output=True, text=True, env=self.env, timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return out, json.loads(r.stdout), r.stderr

    @staticmethod
    def cpus(path):
        return int(re.search(r"^#SBATCH --cpus-per-task=(\d+)$", _read(path), re.M).group(1))

    def test_array_steps_are_sized_to_the_cap_and_say_so(self):
        out, info, err = self.chain("--threads-max", "32")
        for step in ("step2_firstpass", "step4_finalpass"):
            path = os.path.join(out, f"{step}.sbatch")
            self.assertEqual(self.cpus(path), 8, step)
            self.assertIn("--threads 8 ", _read(path), step)
        self.assertEqual(self.cpus(os.path.join(out, "step3_assembly.sbatch")), 64)
        for step in ("step2_firstpass", "step4_finalpass"):
            text = _read(os.path.join(out, f"{step}.sbatch"))
            self.assertIn("#SBATCH --time=8:00:00", text)    # 4 h at 16 CPUs, scaled to 8
            self.assertIn("#SBATCH --mem=64G" if step == "step2_firstpass" else
                          "#SBATCH --mem=48G", text)          # MEM_PER_FILE_GB / FINAL_PASS
        rec = info["cpu_sizing"]
        self.assertEqual((rec["cpus"], rec["concurrent"], rec["per_user_cpu_cap"]), (8, 8, 64))
        self.assertIn("up to 8 of 29 files at once", err)

    def test_step1b_keeps_the_cpus_asked_for(self):
        out, info, _ = self.chain("--threads-max", "32", window=False)
        s1b = os.path.join(out, "step1b_window.sbatch")
        self.assertEqual(self.cpus(s1b), 32)
        self.assertIn("--threads 32", _read(s1b))
        self.assertEqual(info["cpu_sizing"]["step1b_cpus"], 32)

    def test_a_pin_is_written_as_given(self):
        out, info, _ = self.chain("--threads-per-file", "32")
        self.assertEqual(self.cpus(os.path.join(out, "step2_firstpass.sbatch")), 32)
        self.assertTrue(info["cpu_sizing"]["pinned"])

    def test_the_wall_clock_scales_with_fewer_cpus_and_a_pin_is_kept(self):
        import diann_parallel as dp
        self.assertEqual([dp.array_task_hours(c) for c in (4, 8, 10, 16, 32)], [16, 8, 7, 4, 4])
        # real runs at 16 CPUs (sacct): p99 160 min, p99.5 248 min -- the base clears p99 with
        # room; the 32 of 6,182 tasks over 240 min take --time-per-file
        self.assertGreaterEqual(dp.TIME_PER_FILE_HOURS * 60, 240)
        out, info, _ = self.chain("--threads-max", "32", "--time-per-file", "3")
        self.assertIn("#SBATCH --time=3:00:00", _read(os.path.join(out, "step2_firstpass.sbatch")))
        self.assertEqual(info["cpu_sizing"]["time_per_file_hours"]["rule"],
                         "--time-per-file, as given")

    def test_the_memory_request_can_be_raised(self):
        out, info, _ = self.chain("--threads-max", "32", "--mem-per-file", "96")
        for step in ("step2_firstpass", "step4_finalpass"):
            self.assertIn("#SBATCH --mem=96G", _read(os.path.join(out, f"{step}.sbatch")))
        self.assertEqual(info["cpu_sizing"]["mem_per_task_gb"], 96)
        self.assertEqual(info["cpu_sizing"]["mem_gb"]["rule"], "--mem-per-file, as given")

    def test_publicgrp_low_keeps_what_was_asked(self):
        out, info, _ = self.chain("--threads-max", "32", "--partition", "low",
                                  "--account", "publicgrp")
        self.assertEqual(self.cpus(os.path.join(out, "step2_firstpass.sbatch")), 32)
        self.assertIsNone(info["cpu_sizing"]["per_user_cpu_cap"])


class MeasuredDefaultsTests(unittest.TestCase):
    """The floor and the per-file wall clock come from the HIVE measurement (2026-10-01) that
    references/diann_parallel.md tabulates; the doc and the constants must say the same."""

    def test_the_doc_states_the_constants(self):
        import diann_parallel
        doc = _read(os.path.join(os.path.dirname(HERE), "references", "diann_parallel.md"))
        self.assertIn(f"Hence the floor of {run_search.ARRAY_TASK_MIN_CPUS}\n"
                      "(`run_search.ARRAY_TASK_MIN_CPUS`)", doc)
        self.assertIn(f"(`diann_parallel.TIME_PER_FILE_HOURS` = {diann_parallel.TIME_PER_FILE_HOURS}"
                      f" h at `TIME_REFERENCE_CPUS` = {diann_parallel.TIME_REFERENCE_CPUS} CPUs", doc)
        self.assertIn("jobs 24276579–82", doc)
        self.assertIn(f"(default {diann_parallel.MEM_PER_FILE_GB} GB, "
                      "`diann_parallel.MEM_PER_FILE_GB`;", doc)
        self.assertIn(f"({diann_parallel.MEM_FINAL_PASS_GB} GB, "
                      "`diann_parallel.MEM_FINAL_PASS_GB`;", doc)


class RunSearchForwardingTests(unittest.TestCase):
    def test_threads_is_the_ceiling_and_threads_per_file_the_pin(self):
        seen = {}

        class Done(Exception):
            pass

        def fake_run(argv, **kw):
            seen["argv"] = argv
            raise Done()
        a = type("A", (), dict(partition=None, account=None, qos=None, max_simultaneous=None,
                               threads_per_file=None, no_notify=True, no_fran=True,
                               fran_name=None, qc=False, not_qc=False))()
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n")
            orig = run_search.subprocess.run
            run_search.subprocess.run = fake_run
            try:
                for pin in (None, 12):
                    a.threads_per_file = pin
                    a.mem_per_file = 48 if pin else None
                    with self.assertRaises(Done):
                        run_search.run_diann_parallel("diann", cfg, ["/x/a.d", "/x/b.d"],
                                                      "/x/db.fasta", os.path.join(d, "o"), 32, a)
                    argv = seen["argv"]
                    self.assertEqual(argv[argv.index("--threads-max") + 1], "32")
                    if pin:
                        self.assertEqual(argv[argv.index("--mem-per-file") + 1], "48")
                    else:                                  # diann_parallel's default stands
                        self.assertNotIn("--mem-per-file", argv)
                    if pin:
                        self.assertEqual(argv[argv.index("--threads-per-file") + 1], "12")
                    else:
                        self.assertNotIn("--threads-per-file", argv)
            finally:
                run_search.subprocess.run = orig


if __name__ == "__main__":
    unittest.main()
