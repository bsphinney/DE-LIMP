#!/usr/bin/env python3
"""DIA-NN cannot be handed a path with '--' in it, and a run it cannot read must fail the search.

A staff search (2026-10-14, DIA-NN 2.7.0): a run named `<date>__60SPD_DIA--asp3_<well>.d` logged
"unrecognised option [--asp3_<well>.d]" and "skipping .../<date>__60SPD_DIA - invalid raw MS data
format", and the search went on with the other runs. DIA-NN's reader joins every argument with a
space and cuts the result at every '--' (arguments(), diann.cpp at tag 1.8), so no quoting helps.

  * '--' in a folder only: the run is searched through a link under <out>/input_links/ with the
    same file name, so its run name -- and conditions.csv -- do not change;
  * '--' in the file name: refused, naming the run name it would get, unless
    --rename-dashed-runs; then linked under that name;
  * --out / --fasta / --params with '--': refused before anything is written;
  * the 5-step chain's step 5 checks the report's Run column against the input list, as the
    single-shot search job already did, so a run DIA-NN left out fails the search.

The fake DIA-NN below reads its command line the way diann.cpp does, so the links are tested
against the parsing that dropped the run, not against a stand-in that takes argv as given.
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

import diann_parallel as dp  # noqa: E402
from test_single_shot_sbatch import HIGH, LIBFREE_CFG, _Workspace, _read, _write  # noqa: E402

# DIA-NN's arguments() at tag 1.8: argv joined with spaces, cut at every "--"; "--f X" adds a file;
# a file whose extension is not a raw format is skipped. The report holds every file it kept.
FAKE_DIANN_PARSER = r'''
import os, sys
args = " ".join(sys.argv[1:]) + " "
files, out = [], None
start = args.find("--")
while start != -1:
    start += 2
    nxt = args.find("--", start)
    opt = args[start:len(args) if nxt == -1 else nxt]
    if opt.startswith("f "):
        files.append(opt[2:].strip())
    elif opt.startswith("out "):
        out = opt[4:].strip()
    elif opt.split(" ")[0] not in ("fasta", "lib", "threads", "temp", "reanalyse", "matrices",
                                   "qvalue", "xic", "mobilograms", "rt-profiling", "cut",
                                   "var-mod", "mass-acc", "mass-acc-ms1", "window"):
        print("WARNING: unrecognised option [--%s]" % opt.strip())
    start = nxt
kept = []
for f in files:
    if os.path.splitext(f)[1].lower() in (".mzml", ".d", ".raw"):
        kept.append(f)
    else:
        print("WARNING: skipping %s - invalid raw MS data format" % f)
if out:
    runs = [os.path.splitext(os.path.basename(f.rstrip("/")))[0] for f in kept]
    import pyarrow as pa, pyarrow.parquet as pq
    pq.write_table(pa.table({"Run": runs, "Protein.Group": ["P1"] * len(runs)}), out)
'''


def parser_fake(path):
    _write(path, f"#!{sys.executable}\n" + FAKE_DIANN_PARSER, 0o755)


def search_f_args(job):
    """The --f paths on the search job's DIA-NN line."""
    words = [w for line in _read(job).splitlines() if "--f " in line for w in line.split()]
    return [w.strip("'\"") for i, w in enumerate(words) if i and words[i - 1] == "--f"]


def dashed_cohort(w, d, folder="plate--1", names=("sample_0.mzML", "sample_1.mzML")):
    os.makedirs(os.path.join(d, folder), exist_ok=True)
    w.runs = [_write(os.path.join(d, folder, n), "") for n in names]
    return w.runs


class ThePlan(unittest.TestCase):
    def test_a_dashed_folder_links_under_the_same_name(self):
        plan = dp.dashed_input_plan(["/raw/p--1/s1.d", "/raw/ok/s2.d"], "/out")
        self.assertIsNone(plan["refused"])
        self.assertEqual(plan["links"], [{"original": "/raw/p--1/s1.d",
                                          "link": "/out/input_links/s1.d",
                                          "run_name": "s1", "renamed": False}])

    def test_a_dashed_file_name_is_refused_unless_renaming_is_asked_for(self):
        files = ["/raw/x_DIA--asp3_S3.d", "/raw/y.d"]
        plan = dp.dashed_input_plan(files, "/out")
        self.assertIn("1 input file name(s) contain '--'", plan["refused"])
        self.assertIn("would be searched as x_DIA-asp3_S3", plan["refused"])
        self.assertIn(dp.RENAME_DASHED_FLAG, plan["refused"])
        plan = dp.dashed_input_plan(files, "/out", rename=True)
        self.assertIsNone(plan["refused"])
        self.assertEqual(plan["links"][0]["link"], "/out/input_links/x_DIA-asp3_S3.d")
        self.assertTrue(plan["links"][0]["renamed"])

    def test_a_new_name_another_input_has_is_refused(self):
        plan = dp.dashed_input_plan(["/raw/a--b.d", "/raw/a-b.d"], "/out", rename=True)
        self.assertIn("a-b", plan["refused"])
        self.assertIn("merge", plan["refused"])
        self.assertIn("a-b", dp.dashed_input_plan(["/r1/a--b.d", "/r2/a---b.d"], "/out",
                                                  rename=True)["refused"])

    def test_one_dash_is_not_an_option(self):
        """Only '--' starts an option; every path handed to DIA-NN is absolute."""
        self.assertEqual(dp.dashed_input_plan(["/raw/-old/s-1.d", "/raw/a-/-b.d"], "/out"),
                         {"links": [], "refused": None})

    def test_links_are_made_and_a_stale_one_replaced(self):
        with tempfile.TemporaryDirectory() as d:
            src = _write(os.path.join(d, "s1.mzML"), "")
            plan = {"links": [{"original": src, "link": os.path.join(d, "o", "input_links",
                                                                     "s1.mzML")}]}
            os.makedirs(os.path.join(d, "o", "input_links"))
            os.symlink("/elsewhere", plan["links"][0]["link"])
            self.assertEqual(dp.make_input_links(plan), {src: plan["links"][0]["link"]})
            self.assertEqual(os.readlink(plan["links"][0]["link"]), src)


class SingleShot(unittest.TestCase):
    def test_a_dashed_folder_is_searched_through_links_and_every_run_is_reported(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            originals = dashed_cohort(w, d)
            p = w.generate(*HIGH)
            self.assertIn("[run_search] FLAG: 2 input path(s) contain '--'", p.stderr)
            self.assertIn("same file name: the run keeps its name", p.stderr)
            links = [os.path.join(w.out, "input_links", os.path.basename(f)) for f in originals]
            self.assertEqual(search_f_args(w.search_job), links)
            for link, orig in zip(links, originals):
                self.assertEqual(os.readlink(link), orig)
            listing = os.path.join(w.out, "diann_job_input_files.txt")
            self.assertEqual(_read(listing).split(), links)
            with open(os.path.join(w.out, "search_provenance.json")) as fh:
                prov = json.load(fh)
            self.assertEqual(prov["files"], originals, "the record keeps the inputs as given")
            self.assertEqual([(r["original"], r["link"], r["renamed"])
                              for r in prov["input_links"]],
                             [(o, l, False) for o, l in zip(originals, links)])
            # the node check looks at the raw data's own folder, not the links'
            self.assertIn(os.path.join(d, "plate--1"), _read(w.search_job))
            # DIA-NN's own parsing, end to end: both runs reach the report and the check passes
            parser_fake(w.diann)
            r = w.run_job(w.search_job)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertNotIn("unrecognised option", r.stdout)
            self.assertIn("all 2 runs", r.stdout)

    def test_the_bug_itself_a_dashed_path_drops_the_run_and_the_job_fails(self):
        """Handed straight to the parsing fake (what 2.11 did), the run is skipped -- and the
        search job's report check fails it rather than finishing on the other run."""
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.generate(*HIGH)
            dashed = os.path.join(d, "x_DIA--asp3_S3.mzML")
            job = _read(w.search_job).replace(w.runs[1], dashed)
            listing = os.path.join(w.out, "diann_job_input_files.txt")
            _write(listing, f"{w.runs[0]}\n{dashed}\n")
            _write(w.search_job, job)
            parser_fake(w.diann)
            r = w.run_job(w.search_job)
            self.assertIn("unrecognised option [--asp3_S3.mzML]", r.stdout)
            self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertIn("1 of 2 runs", r.stderr)
            self.assertIn("x_DIA--asp3_S3", r.stderr)

    def test_a_dashed_file_name_is_refused_before_anything_is_written(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.runs = [w.runs[0], _write(os.path.join(d, "x_DIA--asp3_S3.mzML"), "")]
            p = w.generate(*HIGH, check=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("REFUSING: 1 input file name(s) contain '--'", p.stderr)
            self.assertIn("would be searched as x_DIA-asp3_S3", p.stderr)
            self.assertIn("--rename-dashed-runs", p.stderr)
            self.assertFalse(os.path.exists(w.out), "nothing may be written")
            self.assertFalse(os.path.exists(w.lib_job))

    def test_rename_dashed_runs_searches_it_under_the_new_name(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            dashed = _write(os.path.join(d, "x_DIA--asp3_S3.mzML"), "")
            w.runs = [w.runs[0], dashed]
            p = w.generate(*HIGH, "--rename-dashed-runs")
            link = os.path.join(w.out, "input_links", "x_DIA-asp3_S3.mzML")
            self.assertIn("(the run is named x_DIA-asp3_S3)", p.stderr)
            self.assertEqual(search_f_args(w.search_job), [w.runs[0], link])
            parser_fake(w.diann)
            r = w.run_job(w.search_job)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertIn("all 2 runs", r.stdout)

    def test_an_out_with_a_double_dash_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.out = os.path.join(d, "search--1")
            p = w.generate(*HIGH, check=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("REFUSING --out", p.stderr)
            self.assertIn("contains '--'", p.stderr)
            self.assertFalse(os.path.exists(w.out))


class Chain(unittest.TestCase):
    def _generate(self, d, runs, *more):
        w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
        out = os.path.join(d, "chain_out")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", w.diann, "--fasta", w.fasta, "--out", out, "--cfg", w.cfg,
                            "--partition", "high", "--account", "genome-center-grp",
                            *more, "--raw", *(runs or w.runs)], capture_output=True, text=True)
        return w, out, p

    def test_diann_parallel_run_directly_links_a_dashed_folder(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "p--1"))
            runs = [_write(os.path.join(d, "p--1", f"s{i}.d"), "") for i in range(3)]
            w, out, p = self._generate(d, runs)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("[diann_parallel] FLAG: 3 input path(s) contain '--'", p.stderr)
            links = [os.path.join(out, "input_links", f"s{i}.d") for i in range(3)]
            self.assertEqual(_read(os.path.join(out, "file_list.txt")).split(), links)
            self.assertEqual([r["link"] for r in json.loads(p.stdout)["input_links"]], links)
            self.assertNotIn("p--1/s", _read(os.path.join(out, "step5_report.sbatch"))
                             .split("--f ", 1)[1])

    def test_diann_parallel_refuses_a_dashed_file_name_and_a_dashed_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_write(os.path.join(d, n), "") for n in ("a.d", "b--c.d")]
            w, out, p = self._generate(d, runs)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("would be searched as b-c", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "step5_report.sbatch")))
            w, out, p = self._generate(d, runs, "--rename-dashed-runs")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn(os.path.join(out, "input_links", "b-c.d"),
                          _read(os.path.join(out, "file_list.txt")))
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=2, ext=".d")
            fasta = _write(os.path.join(d, "db--x.fasta"), ">sp|P1|X\nPEPTIDER\n")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", w.diann, "--fasta", fasta, "--out",
                                os.path.join(d, "o"), "--cfg", w.cfg, "--raw", *w.runs],
                               capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("REFUSING --fasta", p.stderr)

    def test_step5_fails_when_a_run_is_missing_from_the_report(self):
        """The .quant count passes (every run was searched) and DIA-NN still leaves one out of
        the report: step 5 must fail, naming it, not print OK."""
        with tempfile.TemporaryDirectory() as d:
            w, out, p = self._generate(d, None)
            self.assertEqual(p.returncode, 0, p.stderr)
            s5 = os.path.join(out, "step5_report.sbatch")
            self.assertIn("check_report_runs.py", _read(s5))
            self.assertIn(os.path.join(out, "file_list.txt"), _read(s5))
            os.makedirs(os.path.join(out, "quant_step4"), exist_ok=True)
            for i in range(3):
                _write(os.path.join(out, "quant_step4", f"sample_{i}.quant"), "q")
            r = w.run_job(s5, FAKE_REPORT_RUNS="sample_0,sample_2")
            self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertIn("2 of 3 runs", r.stderr)
            self.assertIn("sample_1", r.stderr)
            self.assertNotIn("OK: report built from all", r.stdout)
            r = w.run_job(s5, FAKE_REPORT_RUNS="sample_0,sample_1,sample_2")
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertIn("OK: report holds all 3 runs", r.stdout)


if __name__ == "__main__":
    unittest.main()
