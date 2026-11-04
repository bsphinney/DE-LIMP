#!/usr/bin/env python3
"""
`run_search.py --time HOURS`, and the job size of a non-specific Sage or FragPipe search
(2.14.2, L11 N1/N2).

A non-specific Sage search had its generated job.sh hand-edited (64G -> 192G, 12 h -> 48 h)
before submission: there was --mem (L8) but no --time, and nothing sized the job from the
digest. Now:
  * --time HOURS reaches the header of every single-job route (as --mem does), is recorded as
    search_provenance.json job_time, and the Sage LFQ gate's re-run command keeps it;
  * a Sage or FragPipe search whose digest is non-specific or semi-specific gets a larger
    DEFAULT, from measured use, recorded as a default with what it was measured on (DE-LIMP
    rule 2) and never as the user's choice -- so a re-run recomputes it rather than replaying it;
  * a tryptic search keeps 64G / 12 h (pin).
FragPipe's Java heap follows the job's memory (--ram), and its threads the job's CPUs.

Fakes only: the engines are /bin/true, nothing is submitted.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import run_search  # noqa: E402

FP_WORKFLOW = os.path.join(HERE, "fixtures", "fragpipe24_nonspecific_peptidome",
                           "Nonspecific-peptidome.workflow")
PINNED = ("--qvalue 0.01\n--xic 10\n--mobilograms\n--fasta-search\n--gen-spec-lib\n--predictor\n"
          "--reanalyse\n--mass-acc 15\n--mass-acc-ms1 15\n")
TRYPTIC = {"database": {"enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 2,
                                   "min_len": 7, "max_len": 30}, "fasta": "x"},
           "quant": {"lfq": True}}


def _write(path, text, mode=None):
    with open(path, "w") as fh:
        fh.write(text)
    if mode:
        os.chmod(path, mode)
    return path


def header(path, key):
    m = re.findall(rf"^#SBATCH --{key}=(\S+)$", open(path).read(), re.M)
    assert len(m) == 1, (key, m)
    return m[0]


def search(d, engine, *extra, cfg=None, n_files=2, slurm=False):
    """run_search.py on n_files inputs -> (CompletedProcess, out dir)."""
    fasta = _write(os.path.join(d, "db.fasta"), ">sp|P1|X\nPEPTIDER\n")
    if engine == "diann":
        params = _write(os.path.join(d, "diann.cfg"), PINNED)
        files = []
        for i in range(n_files):
            os.makedirs(os.path.join(d, f"run{i}.d"))
            files.append(os.path.join(d, f"run{i}.d"))
    elif engine == "sage":
        params = _write(os.path.join(d, "sage.json"), json.dumps(cfg or TRYPTIC))
        files = [_write(os.path.join(d, f"run{i}.mzML"), "") for i in range(n_files)]
    else:
        params = os.path.join(d, "search.workflow")
        if cfg is None:
            shutil.copy(FP_WORKFLOW, params)
        else:
            _write(params, cfg)
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
    env.pop("FRAGPIPE_TOOLS_FOLDER", None)
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


NONSPEC = {"database": {"enzyme": {"cleave_at": "", "min_len": 7, "max_len": 30}, "fasta": "x"},
           "quant": {"lfq": True}}


class TimeFlag(unittest.TestCase):
    def test_sage_job_header_takes_time(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            p, out = search(d, "sage", "--sbatch", job, "--time", "36")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(header(job, "time"), "36:00:00")
            prov = provenance(out)
            self.assertEqual(prov["job_time"], "36h")
            self.assertEqual(prov["job_time_source"], "user (--time)")
            self.assertEqual(prov["queue"]["time"], 36)

    def test_diann_single_shot_takes_time_for_both_jobs(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "job.sh")
            p, out = search(d, "diann", "--sbatch", job, "--time", "30")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            for s in ("job_1_lib.sh", "job_2_search.sh"):
                self.assertEqual(header(os.path.join(d, s), "time"), "30:00:00", s)

    def test_a_tryptic_search_keeps_64g_and_12h(self):
        """(pin)"""
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            p, out = search(d, "sage", "--sbatch", job)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(header(job, "mem"), "64G")
            self.assertEqual(header(job, "time"), "12:00:00")
            prov = provenance(out)
            self.assertEqual(prov["job_mem"], "64G (default)")
            self.assertEqual(prov["job_time"], "12h (default)")
            self.assertNotIn("mem", prov["queue"])
            self.assertNotIn("time", prov["queue"])

    def test_the_chain_says_time_is_not_used(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "diann", "--time", "30", n_files=6, slurm=True)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("--time sets the wall time of a single-job search", p.stderr)
            self.assertIn("--time-per-file", p.stderr)
            self.assertIsNone(provenance(out)["job_time"])

    def test_an_inline_search_says_time_changes_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "sage", "--time", "30", n_files=1)
            self.assertIn("--time sets the wall time of a SLURM job", p.stderr)

    def test_the_lfq_gate_s_rerun_keeps_time(self):
        import sage_lfq_check
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "search_provenance.json"), "w") as fh:
                json.dump({"tools": "/t.json", "bundle": "/b.json", "fasta": "/f.fasta",
                           "files": ["/r1.mzML"], "queue": {"mem": "128G", "time": 36},
                           "job_mem": "128G", "job_time": "36h"}, fh)
            cmd = sage_lfq_check.rerun_plan(d, "/c.json", 7.0)["command"]
            self.assertIn("--time 36", cmd)
            self.assertIn("--mem 128", cmd)
            # a digest-sized default is not the user's: a re-run recomputes it
            with open(os.path.join(d, "search_provenance.json"), "w") as fh:
                json.dump({"tools": "/t.json", "bundle": "/b.json", "fasta": "/f.fasta",
                           "files": ["/r1.mzML"], "queue": {},
                           "job_mem": "192G (default for a non-specific digest)",
                           "job_time": "24h (default for a non-specific digest)"}, fh)
            cmd = sage_lfq_check.rerun_plan(d, "/c.json", 7.0)["command"]
            self.assertNotIn("--time", cmd)
            self.assertNotIn("--mem", cmd)


class DigestSizedDefault(unittest.TestCase):
    def test_a_non_specific_sage_search(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            p, out = search(d, "sage", "--sbatch", job, cfg=NONSPEC)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(header(job, "mem"), f"{run_search.BROAD_DIGEST_SAGE_MEM_GB}G")
            self.assertEqual(header(job, "time"), f"{run_search.BROAD_DIGEST_HOURS}:00:00")
            prov = provenance(out)
            self.assertEqual(prov["job_mem_source"],
                             "default for a non-specific digest (measured on a 38-run "
                             "hydrolysate search)")
            self.assertEqual(prov["job_time_source"], prov["job_mem_source"])
            self.assertIn("(default for a non-specific digest)", prov["job_mem"])
            self.assertNotIn("mem", prov["queue"])               # never the user's choice
            self.assertNotIn("time", prov["queue"])
            self.assertIn("non-specific digest", p.stdout + p.stderr)    # said at generation

    def test_a_semi_specific_sage_search(self):
        cfg = json.loads(json.dumps(TRYPTIC))
        cfg["database"]["enzyme"]["semi_enzymatic"] = True
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            p, out = search(d, "sage", "--sbatch", job, cfg=cfg)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(header(job, "mem"), f"{run_search.BROAD_DIGEST_SAGE_MEM_GB}G")
            self.assertIn("semi-specific", provenance(out)["job_mem_source"])

    def test_the_user_s_values_win(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "sage_job.sh")
            p, out = search(d, "sage", "--sbatch", job, "--mem", "300", cfg=NONSPEC)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(header(job, "mem"), "300G")
            self.assertEqual(header(job, "time"), f"{run_search.BROAD_DIGEST_HOURS}:00:00")
            prov = provenance(out)
            self.assertEqual(prov["job_mem_source"], "user (--mem)")
            self.assertIn("non-specific", prov["job_time_source"])

    def test_a_non_specific_fragpipe_search(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "fp_job.sh")
            p, out = search(d, "fragpipe", "--sbatch", job)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            mem = run_search.BROAD_DIGEST_FRAGPIPE_MEM_GB
            self.assertEqual(header(job, "mem"), f"{mem}G")
            self.assertEqual(header(job, "time"), f"{run_search.BROAD_DIGEST_HOURS}:00:00")
            text = open(job).read()
            self.assertIn(f"--ram {run_search.fragpipe_ram_gb(mem)}", text)
            self.assertIn('--threads "${SLURM_CPUS_PER_TASK:-8}"', text)
            self.assertIn("non-specific digest", provenance(out)["job_mem_source"])

    def test_a_tryptic_fragpipe_search_keeps_the_default(self):
        """(pin) and its Java heap still follows the job's memory"""
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "fp_job.sh")
            wf = open(FP_WORKFLOW, newline="").read()
            wf = (wf.replace("msfragger.num_enzyme_termini=0", "msfragger.num_enzyme_termini=2")
                    .replace("msfragger.search_enzyme_name_1=nonspecific",
                             "msfragger.search_enzyme_name_1=stricttrypsin")
                    .replace("msfragger.search_enzyme_cut_1=@", "msfragger.search_enzyme_cut_1=KR"))
            p, out = search(d, "fragpipe", "--sbatch", job, cfg=wf)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(header(job, "mem"), "64G")
            self.assertEqual(header(job, "time"), "12:00:00")
            self.assertIn(f"--ram {run_search.fragpipe_ram_gb(64)}", open(job).read())

    def test_the_ram_rule(self):
        self.assertEqual(run_search.fragpipe_ram_gb(200), 190)     # the measured job's own pair
        self.assertEqual(run_search.fragpipe_ram_gb(64), 54)
        self.assertGreater(run_search.fragpipe_ram_gb(8), 0)

    def test_the_measured_numbers_are_recorded(self):
        m = run_search.BROAD_DIGEST_MEASURED
        self.assertIn("126.2 GiB", m)
        self.assertIn("6,300", m)
        # review F14: 126.2 GiB was the 0:45 job; the 1:02 one peaked at 121.9 GiB
        self.assertNotIn("126.2 GiB in 1:02", m)
        self.assertIn("0:45", m)


class ReviewSizing(unittest.TestCase):
    """Review F11, F15 (of 9abc8390)."""
    PUB_HIGH = {"qos": "publicgrp-high-qos", "cpu": None, "mem_gb": None, "per_job_cpu": 8,
                "per_job_mem_gb": 128, "source": "test"}

    def test_a_job_over_the_queue_s_per_job_cap_is_refused_with_the_options(self):
        with self.assertRaises(SystemExit) as cm:
            run_search.refuse_job_over_queue_cap("192G", self.PUB_HIGH, "default for a "
                                                 "non-specific digest (measured ...)")
        msg = str(cm.exception.code)
        for words in ("REFUSED", "192G", "128G", "publicgrp-high-qos", "--mem", "--partition",
                      "never start"):
            self.assertIn(words, msg)

    def test_under_the_cap_or_no_known_cap_passes(self):
        self.assertIsNone(run_search.refuse_job_over_queue_cap("64G", self.PUB_HIGH, "default"))
        self.assertIsNone(run_search.refuse_job_over_queue_cap(
            "192G", dict(self.PUB_HIGH, per_job_mem_gb=None), "default"))

    def test_the_digest_reading_is_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "sage", "--sbatch", os.path.join(d, "j.sh"), cfg=NONSPEC)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            prov = provenance(out)
            self.assertEqual(prov["job_digest"], "non-specific")
            self.assertEqual(prov["job_digest_source"], "sage.json")

    def test_an_unreadable_digest_is_said(self):
        z = run_search.single_job_sizing("fragpipe", "/nonexistent/x.workflow")
        self.assertIsNone(z["digest"])
        self.assertIn("could not be read", z["digest_source"])
        self.assertIn("could not be read", run_search.sizing_note(z))


if __name__ == "__main__":
    unittest.main()
