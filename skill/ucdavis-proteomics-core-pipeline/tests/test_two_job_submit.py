#!/usr/bin/env python3
"""
A library-free DIA-NN search on SLURM is TWO jobs (predict the library, then search), chained
by <out>/submit.sh -- the route every <=5-file search from a FASTA takes. The live HIVE test of
2026-09-23 (3 Fusion Lumos runs, DIA-NN 2.7.0) found three things wrong with it:

  * submit.sh named the job scripts by the RELATIVE path --sbatch was given as, while it lives
    in <out> and run_search.py prints `bash <out>/submit.sh`: run from anywhere but the --sbatch
    folder, sbatch failed "Unable to open file job_1_lib.sh".
  * it wrote no jobs.txt, so `watch_run.sh --all <out>` answered failed/no_jobs_file for a
    healthy running search -- and step 7b says to resubmit on failed.
  * search_provenance.json recorded "submitted_sbatch": "job.sh", a file that was never
    written, and a job.sh left over from an earlier search stayed where `sbatch job.sh` would
    resubmit it.

run_search.py runs as the orchestrator runs it (relative --sbatch, from the --sbatch folder),
against the fake DIA-NN of test_single_shot_mass_accuracy.py; submit.sh then runs from a
different directory against a fake sbatch that refuses a path it cannot open.
"""
import glob
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)

import test_single_shot_mass_accuracy as ss  # noqa: E402  (fake DIA-NN + cohort fixtures)
from job_env import job_env  # noqa: E402  (env for running generated scripts)

FAKE_SBATCH = r"""#!/bin/bash
for last; do :; done
[ -f "$last" ] || { echo "sbatch: error: Unable to open file $last" >&2; exit 1; }
n=$(( $(cat "$SB_COUNTER" 2>/dev/null || echo 900) + 1 )); echo "$n" > "$SB_COUNTER"
echo "$n"
"""


class TwoJobSubmit(unittest.TestCase):
    # Borrow the harness, not the tests (those run in their own module).
    _setup = ss.SingleShotMassAccTests._setup
    _env = ss.SingleShotMassAccTests._env

    def _generate(self, d, stale=False):
        raws, fasta, cfg, tools, bundle = self._setup(d)
        if stale:
            with open(os.path.join(d, "job.sh"), "w") as fh:
                fh.write("#!/bin/bash\necho an EARLIER search\n")
        out = os.path.join(d, "out")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                            "--tools", tools, "--bundle", bundle, "--params", cfg,
                            "--fasta", fasta, "--out", out, "--files", *raws,
                            "--engine", "diann", "--threads", "8", "--sbatch", "job.sh"],
                           capture_output=True, text=True, env=self._env(d), cwd=d,
                           timeout=240)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        return out, p

    def test_submit_sh_works_from_any_directory_and_writes_jobs_txt(self):
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._generate(d)
            submit = os.path.join(out, "submit.sh")
            text = open(submit).read()
            self.assertIn(os.path.join(os.path.realpath(d), "job_1_lib.sh"),
                          text.replace(d, os.path.realpath(d)))
            fake = os.path.join(d, "fakebin")
            os.makedirs(fake)
            ss.fx._exe(os.path.join(fake, "sbatch"), FAKE_SBATCH)
            env = job_env(d, PATH=fake + os.pathsep + os.environ.get("PATH", ""),
                          SB_COUNTER=os.path.join(d, "counter"))
            elsewhere = tempfile.mkdtemp(dir=d)
            r = subprocess.run(["bash", submit], capture_output=True, text=True, env=env,
                               cwd=elsewhere, timeout=60)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            with open(os.path.join(out, "jobs.txt")) as fh:
                self.assertEqual(fh.read().split(), ["901", "902"])
            self.assertIn("watch_run.sh --all", r.stdout)

    def test_provenance_names_what_was_written_and_a_stale_job_sh_is_moved(self):
        with tempfile.TemporaryDirectory() as d:
            out, p = self._generate(d, stale=True)
            self.assertFalse(os.path.exists(os.path.join(d, "job.sh")),
                             "a job.sh from an earlier search must not stay where "
                             "`sbatch job.sh` would resubmit it")
            moved = glob.glob(os.path.join(d, "job.sh.stale-*"))
            self.assertEqual(len(moved), 1, os.listdir(d))
            self.assertIn("EARLIER", open(moved[0]).read())          # renamed, not deleted
            with open(os.path.join(out, "search_provenance.json")) as fh:
                prov = json.load(fh)
            self.assertEqual(os.path.realpath(prov["submitted_sbatch"]),
                             os.path.realpath(os.path.join(out, "submit.sh")))
            self.assertIn("bash", p.stdout + p.stderr)

    def test_search_names_its_own_temp_unless_the_cfg_does(self):
        """Without --temp, DIA-NN writes each run's .quant beside the raw file -- on HIVE, the
        instrument archive every user searches from. One --temp, never two."""
        with tempfile.TemporaryDirectory() as d:
            out, _ = self._generate(d)
            srch = open(os.path.join(d, "job_2_search.sh")).read()
            self.assertEqual(srch.count("--temp "), 1, srch)
            self.assertIn("--temp " + os.path.join(out, "quant"), srch)
            self.assertTrue(os.path.isdir(os.path.join(out, "quant")))
        with tempfile.TemporaryDirectory() as d:
            raws, fasta, cfg, tools, bundle = self._setup(d)
            mine = os.path.join(d, "mytemp")
            with open(cfg, "a") as fh:
                fh.write(f"\n--temp {mine}\n")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", tools, "--bundle", bundle, "--params", cfg,
                                "--fasta", fasta, "--out", os.path.join(d, "out"),
                                "--files", *raws, "--engine", "diann", "--threads", "8",
                                "--sbatch", "job.sh"],
                               capture_output=True, text=True, env=self._env(d), cwd=d,
                               timeout=240)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            srch = open(os.path.join(d, "job_2_search.sh")).read()
            self.assertEqual(srch.count("--temp "), 1, srch)
            self.assertIn(mine, srch)

    def test_sbatch_naming_the_generated_submit_sh_is_not_moved_aside(self):
        """--sbatch <out>/submit.sh: the run writes that very file, and set_aside must not then
        rename the fresh file to .stale-* (review 2026-09-23)."""
        with tempfile.TemporaryDirectory() as d:
            raws, fasta, cfg, tools, bundle = self._setup(d)
            out = os.path.join(d, "out")
            submit = os.path.join(out, "submit.sh")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                                "--tools", tools, "--bundle", bundle, "--params", cfg,
                                "--fasta", fasta, "--out", out, "--files", *raws,
                                "--engine", "diann", "--threads", "8", "--sbatch", submit],
                               capture_output=True, text=True, env=self._env(d), cwd=d,
                               timeout=240)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertTrue(os.path.isfile(submit), os.listdir(out))
            self.assertEqual(glob.glob(submit + ".stale-*"), [])


if __name__ == "__main__":
    unittest.main()
