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


class Gpus(unittest.TestCase):
    """Brett (2026-10-15): CPU OR GPU requests under 10,000 hours are approved automatically."""

    def test_every_way_slurm_asks_for_gpus(self):
        head = "#!/bin/bash\n#SBATCH --cpus-per-task=4\n#SBATCH --time=10:00:00\n"
        with tempfile.TemporaryDirectory() as d:
            for i, (lines, gpus) in enumerate((
                    ("#SBATCH --gres=gpu:a100:2\n", 2), ("#SBATCH --gres=gpu:1\n", 1),
                    ("#SBATCH --gres=gpu\n", 1), ("#SBATCH --gpus=3\n", 3),
                    ("#SBATCH -G a100:2\n", 2),
                    ("#SBATCH --ntasks=2\n#SBATCH --gpus-per-task=1\n", 2),
                    ("#SBATCH --nodes=2\n#SBATCH --gpus-per-node=4\n", 8),
                    ("#SBATCH --gres=tmpfs:10G\n", 0), ("", 0))):
                path = _write(os.path.join(d, f"j{i}.sh"), head + lines)
                job = cpu_hours.job_reservation(path)
                self.assertEqual(job["gpus"], gpus, lines)
                self.assertEqual(job["gpu_hours"], gpus * 10.0, lines)
            arr = _write(os.path.join(d, "arr.sh"), head + "#SBATCH --gres=gpu:1\n#SBATCH --array=0-4\n")
            rec = cpu_hours.reservation([arr, os.path.join(d, "j8.sh")])
            self.assertEqual(rec["gpu_hours_reserved"], 1 * 10 * 5)
            self.assertEqual(rec["cpu_hours_reserved"], 4 * 10 * 5 + 4 * 10)


class Routes(unittest.TestCase):
    def test_single_shot_one_job(self):
        cfg = [ln for ln in LIBFREE_CFG if ln not in ("--fasta-search", "--predictor")]
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=cfg)
            p = w.generate(*HIGH)
            res = last_json(p.stdout)
            self.assertEqual(res["cpu_hours_reserved"], 4 * run_search.SEARCH_WALL_HOURS)
            self.assertIn(f"[run_search] cpu_hours_reserved: {res['cpu_hours_reserved']}  "
                          "gpu_hours_reserved: 0.0 -- the most this search can use", p.stdout)
            self.assertEqual(res["gpu_hours_reserved"], 0.0)
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
            self.assertEqual(_read(os.path.join(w.out, "submit.sh")).splitlines()[1],
                             cpu_hours.submit_comment(res["cpu_hours_reservation"]))
            self.assertTrue(cpu_hours.submit_comment(res["cpu_hours_reservation"]).startswith(
                f"# cpu_hours_reserved: {res['cpu_hours_reserved']}  gpu_hours_reserved: 0.0 -- "))

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


class EveryGenerator(unittest.TestCase):
    """The board requires cpu_hours >= 0.1 on every job (2026-10-15), so every job script any
    generator writes carries an #SBATCH --time and gives a non-null figure. Searches are covered
    by Routes; these are the rest, and the generators the Routes do not reach. Their log paths are
    double-quoted (an unquoted space makes sbatch refuse the script; HIVE, 2026-10-15)."""

    def assert_bounded(self, path):
        rec = cpu_hours.reservation([path])
        self.assertIsNotNone(rec["cpu_hours_reserved"], path)
        self.assertGreaterEqual(rec["cpu_hours_reserved"], 0.1, path)
        self.assertEqual(rec["gpu_hours_reserved"], 0.0, path)
        return _read(path)

    def test_diatracer_array(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "stage", "r1.d"))
            out = os.path.join(d, "out")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diatracer_parallel.py"),
                                "--stage", os.path.join(d, "stage"), "--out", out, "--jar",
                                "/x/diatracer.jar", "--partition", "high", "--account",
                                "genome-center-grp"], capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assert_bounded(os.path.join(out, "diatracer_array.sbatch"))

    def test_fran_backfill(self):
        import fran_deposit as fd
        with tempfile.TemporaryDirectory() as d:
            class A:
                sbatch_dir, time_budget, max_depth = os.path.join(d, "with space"), 1200, 9
                roots, list, no_homes, apply = None, None, False, False
                partition, account, qos = "high", "genome-center-grp", None
                sbatch_minutes = 30
            res = fd.write_backfill_sbatch(A())
            text = self.assert_bounded(res["sbatch"])
            self.assertIn(f'#SBATCH --output="{A.sbatch_dir}/fran_backfill_%j.log"', text)

    def test_deposit_prep(self):
        import make_deposit as md
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "deposit with space")
            os.makedirs(out)
            md.write_prep_script({"session": "s"}, out, [
                {"kind": "file", "source_path": "/x/a.raw", "upload_name": "a.raw"}])
            text = self.assert_bounded(os.path.join(out, "prepare_upload.sbatch"))
            self.assertIn(f'#SBATCH --output="{out}/prepare_upload_%j.log"', text)

    def test_core_delivery(self):
        import argparse
        import core_submission as cs
        with tempfile.TemporaryDirectory() as d:
            where = os.path.join(d, "PROT_0000 service")
            os.makedirs(where)
            a = argparse.Namespace(summary=os.path.join(d, "s.json"), mode="analysis",
                                   include_raw="no", label=None, service_project=None, force=False)
            path = cs.write_deliver_job(a, None, where, "PROT_0000", 10 ** 12, [])
            text = self.assert_bounded(path)
            self.assertIn(f'#SBATCH --output="{where}/deliver_%j.log"', text)

    def test_podcast_job(self):
        import make_podcast as mp
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "podcast_job.sh")
            mp.write_job(path, sys.executable, os.path.join(d, "script.json"), d, [], [],
                         ("high", "genome-center-grp", None, False),
                         os.path.join(d, "logs", "podcast_job_%j.log"), "t")
            self.assert_bounded(path)

    def test_sage_and_alphadia_single_jobs(self):
        with tempfile.TemporaryDirectory() as d:
            fasta = _write(os.path.join(d, "db.fasta"), ">sp|P1|X\nPEPTIDEK\n")
            runs = [_write(os.path.join(d, f"r{i}.mzML"), "") for i in range(2)]
            for engine, params, acq in (("sage", "sage_config.json", "DDA"),
                                        ("alphadia", "alphadia.yaml", "DIA")):
                tools = _write(os.path.join(d, f"tools_{engine}.json"), json.dumps(
                    {engine: os.path.join(d, f"fake-{engine}"), "versions": {engine: "1.0.0"}}))
                bundle = _write(os.path.join(d, f"m_{engine}.json"), json.dumps(
                    {"engine": {"name": engine}, "acquisition": acq}))
                cfg = _write(os.path.join(d, params), "{}\n")
                job = os.path.join(d, f"{engine}.sh")
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                    "--tools", tools, "--bundle", bundle, "--engine", engine,
                                    "--params", cfg, "--fasta", fasta, "--out",
                                    os.path.join(d, f"o_{engine}"), "--files", *runs,
                                    "--threads", "4", "--sbatch", job, *HIGH],
                                   capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, engine + p.stdout + p.stderr)
                self.assertEqual(last_json(p.stdout)["cpu_hours_reserved"],
                                 4 * run_search.SEARCH_WALL_HOURS, engine)
                self.assert_bounded(job)


if __name__ == "__main__":
    unittest.main()
