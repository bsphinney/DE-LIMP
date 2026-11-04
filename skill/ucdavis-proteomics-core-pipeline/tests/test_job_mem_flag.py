#!/usr/bin/env python3
"""
`run_search.py --mem GB` sets the memory of every single-job search (I46).

Every job emit_sbatch() writes asked for --mem=64G, with no way to change it from the command
line: a Sage search of a large DDA cohort, or a DIA-NN search of a large library, had to be
hand-edited before submission. --mem now reaches the header of every single-job route (Sage,
DIA-NN single shot -- both of its jobs -- FragPipe, AlphaDIA, Radiant single job) and is recorded
as search_provenance.json `job_mem`; without it the job asks for 64G, as before, recorded as
"64G (default)" (DE-LIMP rule 2). The 5-step chain sizes its own jobs (--mem-per-file,
--assembly-mem), so there --mem changes nothing, and run_search says so.

Fakes only: the engines are /bin/true, sbatch is a stub, nothing is submitted.
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

PINNED = ("--qvalue 0.01\n--xic 10\n--mobilograms\n--fasta-search\n--gen-spec-lib\n--predictor\n"
          "--reanalyse\n--mass-acc 15\n--mass-acc-ms1 15\n")


def _write(path, text, mode=None):
    with open(path, "w") as fh:
        fh.write(text)
    if mode:
        os.chmod(path, mode)
    return path


def mem_of(path):
    m = re.findall(r"^#SBATCH --mem=(\S+)$", open(path).read(), re.M)
    assert len(m) == 1, m
    return m[0]


def search(d, engine, n_files, *extra, slurm=False):
    """run_search.py for `engine` on n_files inputs -> (CompletedProcess, out dir)."""
    files = []
    for i in range(n_files):
        p = os.path.join(d, f"run{i}.d")
        os.makedirs(p)
        files.append(p)
    fasta = _write(os.path.join(d, "db.fasta"), ">sp|P1|X\nPEPTIDER\n")
    params = (_write(os.path.join(d, "diann.cfg"), PINNED) if engine == "diann" else
              _write(os.path.join(d, "sage.json"), json.dumps({"database": {"fasta": "x"}})))
    if engine == "sage":
        files = [_write(os.path.join(d, f"run{i}.mzML"), "") for i in range(n_files)]
    tools = _write(os.path.join(d, "tools.json"), json.dumps({engine: "/bin/true"}))
    bundle = _write(os.path.join(d, "bundle.json"), json.dumps(
        {"engine": {"name": engine}, "acquisition": "DIA" if engine == "diann" else "DDA"}))
    env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
    path = [os.path.dirname(sys.executable), "/usr/bin", "/bin"]
    if slurm:
        b = os.path.join(d, "bin")
        os.makedirs(b)
        _write(os.path.join(b, "sbatch"), "#!/bin/sh\necho 1\n", 0o755)
        path.insert(0, b)
    env.update(PATH=os.pathsep.join(path), HPCCF_SLURM_BIN=os.path.join(d, "no_slurm"),
               UCD_HPC_CLUSTER_CONF=os.path.join(d, "no_cluster_file"))
    out = os.path.join(d, "out")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                        "--tools", tools, "--bundle", bundle, "--params", params,
                        "--fasta", fasta, "--out", out, "--files", *files,
                        "--engine", engine, "--threads", "8", *extra],
                       capture_output=True, text=True, env=env, cwd=d, timeout=240)
    return p, out


def provenance(out):
    with open(os.path.join(out, "search_provenance.json")) as fh:
        return json.load(fh)


class JobMem(unittest.TestCase):
    def test_sage_job_header_takes_mem(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            p, out = search(d, "sage", 2, "--sbatch", job, "--mem", "128")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(mem_of(job), "128G")
            self.assertEqual(provenance(out)["job_mem"], "128G")

    def test_diann_single_shot_takes_mem(self):
        """Both jobs of the library-free search: the library job and the search job."""
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "job.sh")
            p, out = search(d, "diann", 2, "--sbatch", job, "--mem", "96")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            for s in ("job_1_lib.sh", "job_2_search.sh"):
                self.assertEqual(mem_of(os.path.join(d, s)), "96G", s)
            self.assertEqual(provenance(out)["job_mem"], "96G")
            self.assertEqual(provenance(out)["queue"]["mem"], "96G")

    def test_default_is_recorded_as_default(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            p, out = search(d, "sage", 2, "--sbatch", job)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(mem_of(job), run_search.JOB_MEM_DEFAULT)
            self.assertEqual(provenance(out)["job_mem"], f"{run_search.JOB_MEM_DEFAULT} (default)")
            self.assertNotIn("mem", provenance(out)["queue"])

    def test_chain_notes_that_mem_is_not_used(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "diann", 6, "--mem", "77", slurm=True)    # no chain default
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("parallel routing: YES", p.stdout)
            self.assertIn("--mem sets the memory of a single-job search", p.stderr)
            self.assertIn("--mem-per-file", p.stderr)
            self.assertIn("--assembly-mem", p.stderr)
            self.assertIsNone(provenance(out)["job_mem"])
            for f in os.listdir(out):
                if f.endswith(".sbatch"):
                    self.assertNotIn("--mem=77G", open(os.path.join(out, f)).read(), f)

    def test_the_radiant_chain_note_names_radiant_s_flags(self):
        """L8R-F9: Radiant's per-file array is sized by radiant_parallel.py, not by the DIA-NN
        chain's --mem-per-file / --assembly-mem that the note used to name."""
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "radiant", 2, "--mem", "77", "--sbatch",
                            os.path.join(d, "radiant_job.sh"), slurm=True)
            self.assertIn("--mem sets the memory of a single-job search", p.stderr)
            self.assertIn("radiant_parallel.py", p.stderr)
            self.assertIn("--fulcrum-mem", p.stderr)
            self.assertNotIn("--assembly-mem", p.stderr)

    def test_an_inline_search_says_mem_changes_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "sage", 1, "--mem", "128")
            self.assertIn("--mem sets the memory of a SLURM job", p.stderr)

    def test_mem_must_be_a_positive_whole_number(self):
        with tempfile.TemporaryDirectory() as d:
            p, _ = search(d, "sage", 1, "--sbatch", os.path.join(d, "j.sh"), "--mem", "0")
            self.assertEqual(p.returncode, 2)
            self.assertIn("--mem", p.stderr)


if __name__ == "__main__":
    unittest.main()
