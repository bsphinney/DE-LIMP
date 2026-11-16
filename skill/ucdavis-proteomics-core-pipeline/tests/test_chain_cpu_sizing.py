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

    def test_single_jobs_leave_headroom_under_a_per_user_cap(self):
        self.assertEqual(run_search.single_job_cpus(64, GC_HIGH), 48)    # 64 - 64 // 4
        self.assertEqual(run_search.single_job_cpus(96, GC_HIGH), 48)
        self.assertEqual(run_search.single_job_cpus(16, GC_HIGH), 16)    # below it: as asked
        self.assertEqual(run_search.single_job_cpus(32, GC_HIGH), 32)
        self.assertEqual(run_search.single_job_cpus(64, PUB_HIGH), 8)    # per-job cap still wins
        self.assertEqual(run_search.single_job_cpus(64, PUB_LOW), 64)    # no per-user cap

    def test_tres_strings(self):
        self.assertEqual(run_search._tres("cpu=64,gres/gpu=0,mem=1T", "cpu"), 64)
        self.assertEqual(run_search._tres("cpu=64,gres/gpu=0,mem=1T", "mem"), 1024)
        self.assertEqual(run_search._tres("cpu=8,gres/gpu=1,mem=128G", "mem"), 128)
        self.assertEqual(run_search._tres("mem=512000M", "mem"), 500)
        self.assertEqual(run_search._tres("mem=2048", "mem"), 2)          # MB when no unit
        self.assertIsNone(run_search._tres("", "cpu"))
        self.assertIsNone(run_search._tres("gres/gpu=0", "cpu"))


class SizeCpusTests(unittest.TestCase):
    """diann_parallel.size_cpus: the chain's single jobs (steps 1, 1b, 3, 5) on a capped queue."""

    def run_size(self, limits, **kw):
        import argparse
        import contextlib
        import io
        import diann_parallel
        a = argparse.Namespace(threads_per_file=None, threads_max=16, max_simultaneous=None,
                               mem_per_file=None, time_per_file=None, libpred_cpus=16,
                               assembly_cpus=64)
        for k, v in kw.items():
            setattr(a, k, v)
        orig = run_search.user_limits
        run_search.user_limits = lambda *q: dict(limits)
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err):
                rec, _probe = diann_parallel.size_cpus(a, 29, ("high", "x", None),
                                                       ("high", "x", None))
        finally:
            run_search.user_limits = orig
        return a, rec, err.getvalue()

    def test_default_assembly_asks_for_48_of_a_64_cap(self):
        a, rec, err = self.run_size(GC_HIGH)
        self.assertEqual(a.assembly_cpus, 48)
        self.assertEqual(a.libpred_cpus, 16)
        low = rec["single_jobs_lowered"]["--assembly-cpus"]
        self.assertEqual((low["asked"], low["used"]), (64, 48))
        self.assertIn("1/4 of the per-user cap", low["reason"])
        self.assertIn("--assembly-cpus lowered 64 -> 48", err)

    def test_a_per_job_cap_of_8_gives_8(self):
        a, rec, _ = self.run_size(PUB_HIGH)
        self.assertEqual(a.assembly_cpus, 8)
        self.assertNotIn("1/4", rec["single_jobs_lowered"]["--assembly-cpus"]["reason"])

    def test_no_per_user_cap_keeps_64(self):
        a, rec, _ = self.run_size(PUB_LOW)
        self.assertEqual(a.assembly_cpus, 64)
        self.assertNotIn("single_jobs_lowered", rec)

    def test_an_explicit_32_is_kept(self):
        a, rec, _ = self.run_size(GC_HIGH, assembly_cpus=32)
        self.assertEqual(a.assembly_cpus, 32)
        self.assertNotIn("single_jobs_lowered", rec)

    def test_an_explicit_64_is_lowered_and_said(self):
        a, rec, err = self.run_size(GC_HIGH, assembly_cpus=64)
        self.assertEqual(a.assembly_cpus, 48)
        self.assertIn("lowered 64 -> 48", err)

    def test_step_1b_probe_leaves_headroom_too(self):
        _, rec, _ = self.run_size(GC_HIGH, threads_max=64)
        self.assertEqual(rec["step1b_cpus"], 48)


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
        # X2: step 3 is sized from the files (29 -> 8 CPUs); step 5 keeps 64 asked, lowered to 48
        # (a quarter of the cap stays free)
        self.assertEqual(self.cpus(os.path.join(out, "step3_assembly.sbatch")), 8)
        self.assertEqual(self.cpus(os.path.join(out, "step5_report.sbatch")), 48)
        for step in ("step2_firstpass", "step4_finalpass"):
            text = _read(os.path.join(out, f"{step}.sbatch"))
            self.assertIn("#SBATCH --time=8:00:00", text)    # 4 h at 16 CPUs, scaled to 8
            self.assertIn("#SBATCH --mem=64G" if step == "step2_firstpass" else
                          "#SBATCH --mem=48G", text)          # MEM_PER_FILE_GB / FINAL_PASS
        rec = info["cpu_sizing"]
        self.assertEqual((rec["cpus"], rec["concurrent"], rec["per_user_cpu_cap"]), (8, 8, 64))
        self.assertIn("up to 8 of 29 files at once", err)

    def test_step_3_is_sized_from_the_files_and_step_5_keeps_its_own_size(self):
        out, info, _ = self.chain("--threads-max", "32")       # 29 files
        s3 = _read(os.path.join(out, "step3_assembly.sbatch"))
        s5 = _read(os.path.join(out, "step5_report.sbatch"))
        self.assertIn("#SBATCH --mem=32G", s3)
        self.assertIn("--threads 8 ", s3)
        self.assertNotIn("#SBATCH --time=12:00:00", s3)
        self.assertIn("#SBATCH --mem=128G", s5)
        self.assertIn("#SBATCH --time=12:00:00", s5)
        self.assertEqual(info["cpu_sizing"]["step3_sizing"]["mem"], 32)

    def test_step_3_takes_the_users_assembly_flags(self):
        out, _, _ = self.chain("--threads-max", "32", "--assembly-cpus", "24",
                               "--assembly-mem", "96", "--assembly-time", "20")
        s3 = _read(os.path.join(out, "step3_assembly.sbatch"))
        self.assertIn("#SBATCH --mem=96G", s3)
        self.assertIn("#SBATCH --time=20:00:00", s3)
        self.assertEqual(self.cpus(os.path.join(out, "step3_assembly.sbatch")), 24)

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
        self.assertIn("(HIVE, 2026-10-01, four test jobs)", doc)
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


class AssemblySizing(unittest.TestCase):
    """X2 5/5: step 3 (empirical library assembly) asked for 48 CPU / 128G / 12 h and used 2:24
    and 5 GB on 77 files, after waiting hours on the per-user CPU cap."""

    def test_the_measured_case_gets_a_small_job(self):
        import diann_parallel as dp
        cpus, mem, hours, rule = dp.assembly_sizing(77)
        self.assertEqual((cpus, mem), (16, 32))
        self.assertGreaterEqual(hours, 3)             # measured 2:24 at 48 CPUs; at 16 it is longer
        self.assertLessEqual(hours, 24)
        self.assertIn("77 files", rule)
        self.assertIn("2:24", rule)

    def test_it_scales_with_the_files_and_stays_inside_the_old_request(self):
        import diann_parallel as dp
        small, big, huge = (dp.assembly_sizing(n) for n in (5, 300, 5000))
        self.assertEqual(small[0], 8)
        self.assertLessEqual(huge[0], 16)
        self.assertLessEqual(huge[1], 128)
        self.assertLessEqual(huge[2], 48)
        self.assertGreater(big[1], small[1])
        self.assertGreater(big[2], small[2])
        self.assertGreaterEqual(small[2], 2)

    def run_size(self, n, **kw):
        import argparse
        import contextlib
        import io
        import diann_parallel
        a = argparse.Namespace(threads_per_file=None, threads_max=16, max_simultaneous=None,
                               mem_per_file=None, time_per_file=None, libpred_cpus=16,
                               assembly_cpus=64, assembly_mem=128, assembly_time=12,
                               assembly_given={"cpus": False, "mem": False, "time": False})
        for k, v in kw.items():
            setattr(a, k, v)
        orig = run_search.user_limits
        run_search.user_limits = lambda *q: dict(GC_HIGH)
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                rec, _ = diann_parallel.size_cpus(a, n, ("high", "x", None), ("high", "x", None))
        finally:
            run_search.user_limits = orig
        return a, rec

    def test_size_cpus_sets_step_3_and_records_the_rule(self):
        a, rec = self.run_size(77)
        self.assertEqual((a.s3_cpus, a.s3_mem), (16, 32))
        self.assertEqual(a.assembly_cpus, 48)         # step 5: the old ask, under the cap rule
        self.assertEqual(a.assembly_mem, 128)
        r = rec["step3_sizing"]
        self.assertEqual((r["cpus"], r["mem"]), (16, 32))
        self.assertIn("default for 77 files", r["rule"])

    def test_the_users_assembly_flags_win_for_step_3(self):
        a, rec = self.run_size(77, assembly_cpus=32, assembly_mem=96, assembly_time=20,
                               assembly_given={"cpus": True, "mem": True, "time": True})
        self.assertEqual((a.s3_cpus, a.s3_mem, a.s3_time), (32, 96, 20))
        self.assertIn("as given", rec["step3_sizing"]["rule"])


if __name__ == "__main__":
    unittest.main()
