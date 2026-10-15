#!/usr/bin/env python3
"""Every search states the most CPU-hours it can use, from the job scripts it actually wrote.

Brett (2026-10-15): on the Core board a request under 10,000 CPU-hours is approved automatically;
10,000+ needs a person. cpu_hours.py (the one rule) sums CPUs per task x time limit x array tasks
over every job script a search submits; run_search.py, diann_parallel.py and radiant_parallel.py
print it as `cpu_hours_reserved`, and each submit.sh carries it in a header comment. The expected
figures below are the arithmetic of the sizes each test asks for.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import cpu_hours  # noqa: E402
import run_search  # noqa: E402
from test_single_shot_sbatch import HIGH, LIBFREE_CFG, _Workspace, _read, _write  # noqa: E402


def last_json(stdout):
    return json.loads(stdout[stdout.rindex("\n{") + 1:] if "\n{" in stdout else stdout)


class TheRule(unittest.TestCase):
    def test_slurm_time_formats(self):
        for value, hours in (("12:00:00", 12), ("90", 1.5), ("30:00", 0.5), ("1-06", 30),
                             ("2-00:30", 48.5), ("1-00:00:00", 24), ("4:00:00", 4)):
            self.assertAlmostEqual(cpu_hours.slurm_hours(value), hours, msg=value)

    def test_array_sizes(self):
        for spec, n in (("0-9%5", 10), ("0-2", 3), ("1,3,5-7", 5), ("0-10:2", 6), ("4", 1)):
            self.assertEqual(cpu_hours.array_tasks(spec), n, spec)

    def test_a_job_is_cpus_times_hours_times_tasks(self):
        with tempfile.TemporaryDirectory() as d:
            a = _write(os.path.join(d, "a.sbatch"), "#!/bin/bash\n#SBATCH --cpus-per-task=8\n"
                       "#SBATCH --time=4:00:00\n#SBATCH --array=0-9%3\necho hi\n")
            b = _write(os.path.join(d, "b.sbatch"), "#!/bin/bash\n#SBATCH -c 32\n#SBATCH -t 12:00:00\n")
            rec = cpu_hours.reservation([a, b])
            self.assertEqual(rec["cpu_hours_reserved"], 8 * 4 * 10 + 32 * 12)
            self.assertIn("the most this search can use", rec["what"])
            c = _write(os.path.join(d, "c.sbatch"), "#!/bin/bash\n#SBATCH --cpus-per-task=4\n")
            rec = cpu_hours.reservation([a, c])
            self.assertIsNone(rec["cpu_hours_reserved"], "no time limit: nothing bounds it")
            self.assertEqual(rec["no_time_limit"], [c])


class Routes(unittest.TestCase):
    def test_single_shot_one_job(self):
        cfg = [ln for ln in LIBFREE_CFG if ln not in ("--fasta-search", "--predictor")]
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=cfg)
            p = w.generate(*HIGH)
            res = last_json(p.stdout)
            self.assertEqual(res["cpu_hours_reserved"], 4 * run_search.SEARCH_WALL_HOURS)
            self.assertIn(f"[run_search] cpu_hours_reserved: {res['cpu_hours_reserved']} -- the "
                          "most this search can use", p.stdout)
            with open(os.path.join(w.out, "search_provenance.json")) as fh:
                self.assertEqual(json.load(fh)["cpu_hours_reserved"], res["cpu_hours_reserved"])

    def test_single_shot_two_jobs_and_its_submit_sh(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            p = w.generate(*HIGH)
            res = last_json(p.stdout)
            jobs = res["cpu_hours_reservation"]["jobs"]
            self.assertEqual([os.path.basename(j["script"]) for j in jobs],
                             ["diann_job_1_lib.sh", "diann_job_2_search.sh"])
            self.assertEqual(res["cpu_hours_reserved"], sum(4 * j["hours"] for j in jobs))
            self.assertTrue(_read(os.path.join(w.out, "submit.sh")).splitlines()[1].startswith(
                f"# cpu_hours_reserved: {res['cpu_hours_reserved']} -- the most this search"))

    def test_the_five_step_chain(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
            out = os.path.join(d, "chain")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", w.diann, "--fasta", w.fasta, "--out", out, "--cfg",
                                w.cfg, "--partition", "high", "--account", "genome-center-grp",
                                "--threads-per-file", "8", "--time-per-file", "4",
                                "--assembly-cpus", "32", "--assembly-time", "12",
                                "--libpred-cpus", "16", "--libpred-time", "4", "--raw", *w.runs],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            res = json.loads(p.stdout)
            # step 1 16x4, steps 2 and 4 8x4x3 runs, steps 3 and 5 32x12
            self.assertEqual(res["cpu_hours_reserved"], 16 * 4 + 2 * (8 * 4 * 3) + 2 * (32 * 12))
            self.assertEqual(_read(os.path.join(out, "submit.sh")).splitlines()[1],
                             cpu_hours.submit_comment(res["cpu_hours_reservation"]))

    def test_the_radiant_chain(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_write(os.path.join(d, f"r{i}.mzML"), "") for i in range(4)]
            listing = _write(os.path.join(d, "list.txt"), "\n".join(runs) + "\n")
            fasta = _write(os.path.join(d, "db.fasta"), ">sp|P1|X\nPEPTIDEK\n")
            config = _write(os.path.join(d, "x.radiantConfig"), "")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "radiant_parallel.py"),
                                "--runtime", "apptainer", "--image", "/x/r.sif", "--raw-list",
                                listing, "--fasta", fasta, "--config", config, "--out",
                                os.path.join(d, "out"), "--diann", "/opt/diann",
                                "--libpred-cpus", "16", "--libpred-time", "4",
                                "--threads-per-file", "8", "--time-per-file", "2",
                                "--fulcrum-cpus", "32", "--fulcrum-time", "8",
                                "--partition", "high", "--account", "genome-center-grp"],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(json.loads(p.stdout)["cpu_hours_reserved"],
                             16 * 4 + 8 * 2 * 4 + 32 * 8)

    def test_fragpipe(self):
        with tempfile.TemporaryDirectory() as d:
            tools = _write(os.path.join(d, "tools.json"), json.dumps(
                {"fragpipe": os.path.join(d, "fake-fragpipe"), "versions": {"fragpipe": "24.0"}}))
            bundle = _write(os.path.join(d, "manifest.json"), json.dumps(
                {"engine": {"name": "fragpipe"}, "acquisition": "DDA"}))
            wf = _write(os.path.join(d, "LFQ.workflow"), "diann.run-dia-nn=false\n")
            fasta = _write(os.path.join(d, "db.fasta"), ">sp|P1|X\nPEPTIDEK\n")
            runs = [_write(os.path.join(d, f"r{i}.mzML"), "") for i in range(2)]
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", tools, "--bundle", bundle, "--engine", "fragpipe",
                                "--params", wf, "--fasta", fasta, "--out", os.path.join(d, "o"),
                                "--files", *runs, "--threads", "6", "--sbatch",
                                os.path.join(d, "fp.sh"), *HIGH],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(last_json(p.stdout)["cpu_hours_reserved"],
                             6 * run_search.SEARCH_WALL_HOURS)

    def test_inline_reserves_nothing(self):
        self.assertEqual(run_search.search_cpu_hours({"ran": True}), {"ran": True})


if __name__ == "__main__":
    unittest.main()
