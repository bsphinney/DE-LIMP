#!/usr/bin/env python3
"""A run whose FILE NAME holds '--' is decided at step 2, not discovered at step 7.

DIA-NN reads '--' as the start of an option, so such a run is searched only under another name
(run_search.py --rename-dashed-runs: each run of dashes made one). conditions.csv is built at
step 3 from the raw names, so learning that at step 7 meant rebuilding the design afterwards.

  * check_report_runs.renamed_dashed_name / dashed_file_names: THE rule, which diann_parallel's
    links use too;
  * detect_acquisition.py (step 2) lists each such run with the name it would get, as a warning
    (needs_confirmation) and as `dashed_run_names`;
  * collect_conditions.py --rename-dashed-runs (and core_submission.py conditions) write those
    names into conditions.csv; --validate passes them against the report and names a row that
    still holds the raw name;
  * FragPipe's DIA-NN step reads its runs from a --cfg file list that DIA-NN also cuts at '--'
    (HIVE srun 24807938), and has no renaming: run_search.py refuses such a path there;
  * a run missing from the report says what to do, and not to resubmit.
"""
import csv
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

import check_report_runs as crr  # noqa: E402
import diann_parallel as dp  # noqa: E402
from test_core_submission import child_env, run as cs_run, write_summary  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_PYARROW = True
except ImportError:
    HAVE_PYARROW = False

DASHED = "x_DIA--asp3_S3-H1_1_23896"
RENAMED = "x_DIA-asp3_S3-H1_1_23896"


def write(path, text=""):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


class TheRule(unittest.TestCase):
    def test_one_rule_for_every_step(self):
        self.assertEqual(crr.renamed_dashed_name(DASHED + ".d"), RENAMED + ".d")
        self.assertEqual(crr.renamed_dashed_name("a---b"), "a-b")
        recs = crr.dashed_file_names([f"/raw/{DASHED}.d", "/raw/p--1/s1.d", "/raw/ok.d"])
        self.assertEqual(recs, [{"file": f"/raw/{DASHED}.d", "run_name": DASHED,
                                 "renamed_run_name": RENAMED}])
        # the chain's links are named by the same rule
        plan = dp.dashed_input_plan([f"/raw/{DASHED}.d"], "/out", rename=True)
        self.assertEqual(plan["links"][0]["run_name"], RENAMED)
        self.assertIs(dp.renamed_dashed_name, crr.renamed_dashed_name)
        self.assertEqual(dp.RENAME_DASHED_FLAG, crr.RENAME_DASHED_FLAG)


class DetectAcquisition(unittest.TestCase):
    def test_step2_names_the_run_and_its_new_name(self):
        with tempfile.TemporaryDirectory() as d:
            dashed = write(os.path.join(d, DASHED + ".mzML"))
            folder = write(os.path.join(d, "plate--1", "s1.mzML"))
            plain = write(os.path.join(d, "s2.mzML"))
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                                dashed, folder, plain], capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual(out["dashed_run_names"], [{"file": dashed, "run_name": DASHED,
                                                        "renamed_run_name": RENAMED}])
            self.assertTrue(out["needs_confirmation"])
            warned = {f["file"]: f.get("warnings") or [] for f in out["files"]}
            self.assertTrue(any(f"searches it as '{RENAMED}'" in w and "--rename-dashed-runs" in w
                                for w in warned[dashed]), warned[dashed])
            self.assertFalse(any("'--'" in w for w in warned[folder] + warned[plain]))
            self.assertIn(f"WARNING: {dashed}: its file name holds '--'", p.stderr)


class Conditions(unittest.TestCase):
    def cc(self, *args):
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "collect_conditions.py"), *args],
                           capture_output=True, text=True)
        return p, (json.loads(p.stdout) if p.stdout.strip().startswith("{") else None)

    def test_list_runs_and_map_write_the_searched_names(self):
        runs = f"{DASHED},ctrl_2,trt_1,trt_2"
        p, out = self.cc("--list-runs", "--runs", runs, "--rename-dashed-runs")
        self.assertEqual(out["runs"], [RENAMED, "ctrl_2", "trt_1", "trt_2"])
        p, out = self.cc("--list-runs", "--runs", runs)
        self.assertEqual(out["runs"][0], DASHED, "without the flag the raw name stands")
        with tempfile.TemporaryDirectory() as d:
            csv_out = os.path.join(d, "conditions.csv")
            # the user's sheet names the RAW file: matching ignores punctuation, so it still finds it
            intent = {"mapping": {DASHED: "ctrl", "ctrl_2": "ctrl", "trt_1": "trt", "trt_2": "trt"}}
            p, out = self.cc("--map", csv_out, "--runs", runs, "--rename-dashed-runs",
                             "--mapping-json", json.dumps(intent))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            with open(csv_out, newline="") as fh:
                rows = {r["File.Name"]: r["Group"] for r in csv.DictReader(fh)}
            self.assertEqual(rows[RENAMED], "ctrl")
            self.assertNotIn(DASHED, rows)

    @unittest.skipUnless(HAVE_PYARROW, "the report is read with pyarrow")
    def test_validate_matches_the_renamed_names_and_names_a_raw_one(self):
        with tempfile.TemporaryDirectory() as d:
            report = os.path.join(d, "report.parquet")
            pq.write_table(pa.table({"Run": [RENAMED, "ctrl_2", "trt_1", "trt_2"]}), report)
            good = write(os.path.join(d, "good.csv"), "File.Name,Group\n"
                         f"{RENAMED},ctrl\nctrl_2,ctrl\ntrt_1,trt\ntrt_2,trt\n")
            p, out = self.cc("--validate", good, "--against", report)
            self.assertEqual(p.returncode, 0, p.stdout)
            self.assertTrue(out["valid"])
            raw = write(os.path.join(d, "raw.csv"), "File.Name,Group\n"
                        f"{DASHED},ctrl\nctrl_2,ctrl\ntrt_1,trt\ntrt_2,trt\n")
            p, out = self.cc("--validate", raw, "--against", report)
            self.assertNotEqual(p.returncode, 0)
            hint = [x for x in out["problems"] if "--rename-dashed-runs" in x]
            self.assertEqual(len(hint), 1, out["problems"])
            self.assertIn(f"{DASHED} -> {RENAMED}", hint[0])

    def test_core_submission_conditions_takes_the_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = child_env(tmp)
            summary, _ = write_summary(tmp)
            tsv = os.path.join(tmp, "sample_files.tsv")
            names = {"KG1": DASHED, "KG2": "03122025__60SPD_DIA-KG2_S3-A1_1_24001",
                     "KG13": "03122025__60SPD_DIA-KG13_S3-A2_1_24002",
                     "KG14": "03122025__60SPD_DIA-KG14_S3-A3_1_24003"}
            conds = {"KG1": "ctrl", "KG2": "ctrl", "KG13": "treat", "KG14": "treat"}
            with open(tsv, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t", lineterminator="\n")
                w.writerow(["unique_id", "sample_name", "condition_name", "file", "acquired",
                            "date_source", "status", "alternates", "note"])
                for uid, name in names.items():
                    f = os.path.join(tmp, "rawdir", name + ".d")
                    os.makedirs(f)
                    w.writerow([uid, f"sample {uid}", conds[uid], f, "2025-03-12", "filename",
                                "matched", "", ""])
            out = os.path.join(tmp, "conditions.csv")
            rc, js, p = cs_run(["conditions", "--summary", summary, "--sample-files", tsv,
                                "--out", out, "--rename-dashed-runs"], env)
            self.assertEqual(rc, 0, p.stdout + p.stderr)
            with open(out, newline="") as fh:
                rows = {r["File.Name"]: r["Group"] for r in csv.DictReader(fh)}
            self.assertEqual(rows.get(RENAMED), "ctrl", rows)
            self.assertNotIn(DASHED, rows)


class FragPipeRoute(unittest.TestCase):
    """FragPipe 24.0's CmdDiann writes `--f <path>` lines to filelist_diann.txt and hands DIA-NN
    `--cfg <filelist>--`; DIA-NN cuts that text at every '--' (HIVE srun 24807938)."""

    def generate(self, d, files, runs_diann=True, acquisition="DIA"):
        tools = write(os.path.join(d, "tools.json"), json.dumps(
            {"fragpipe": os.path.join(d, "fake-fragpipe"), "versions": {"fragpipe": "24.0"}}))
        bundle = write(os.path.join(d, "manifest.json"), json.dumps(
            {"engine": {"name": "fragpipe"}, "acquisition": acquisition}))
        wf = write(os.path.join(d, "DIA_SpecLib_Quant_diaPASEF.workflow"),
                   f"workflow.misc.save-sdrf=false\ndiann.run-dia-nn={str(runs_diann).lower()}\n")
        fasta = write(os.path.join(d, "search.fasta"), ">sp|P1|X\nPEPTIDER\n")
        out = os.path.join(d, "search_out")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                            "--tools", tools, "--bundle", bundle, "--engine", "fragpipe",
                            "--params", wf, "--fasta", fasta, "--out", out, "--files", *files,
                            "--threads", "4", "--sbatch", os.path.join(d, "fp_job.sh")],
                           capture_output=True, text=True)
        return p, out

    def dotd(self, d, *parts):
        path = os.path.join(d, *parts)
        os.makedirs(path)
        return path

    def test_a_dashed_file_name_is_refused_when_the_workflow_runs_diann(self):
        with tempfile.TemporaryDirectory() as d:
            files = [self.dotd(d, "raw", DASHED + ".d"), self.dotd(d, "raw", "ok.d")]
            p, out = self.generate(d, files)
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("REFUSING: this FragPipe workflow runs DIA-NN", p.stderr)
            self.assertIn(f"would be searched as {RENAMED}", p.stderr)
            self.assertIn("--engine diann", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "fragpipe.fp-manifest")))

    def test_a_dashed_folder_is_fine_on_the_diatracer_route_but_not_off_it(self):
        with tempfile.TemporaryDirectory() as d:
            # dia-PASEF .d: DIA-NN gets the staged link in <out>/diatracer_stage/, no '--' in it
            files = [self.dotd(d, "p--1", "s1.d"), self.dotd(d, "p--1", "s2.d")]
            p, out = self.generate(d, files)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            with open(os.path.join(out, "fragpipe.fp-manifest")) as fh:
                self.assertNotIn("--", fh.read())
        with tempfile.TemporaryDirectory() as d:
            # mzML goes to FragPipe as given, folder and all
            files = [write(os.path.join(d, "p--1", f"s{i}.mzML")) for i in range(2)]
            p, out = self.generate(d, files)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("(a folder holds '--')", p.stderr)

    def test_no_diann_in_the_workflow_no_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            files = [write(os.path.join(d, "raw", DASHED + ".mzML")),
                     write(os.path.join(d, "raw", "ok.mzML"))]
            p, out = self.generate(d, files, runs_diann=False, acquisition="DDA")
            self.assertNotIn("REFUSING: this FragPipe workflow", p.stderr)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)


class MissingRunMessage(unittest.TestCase):
    def test_it_says_what_to_do_and_not_to_resubmit(self):
        ok, msg = crr.verify("/s/report.parquet", ["/r/s1.d", "/r/blank.d"],
                             parquet_runs=lambda _r: {"s1"})
        self.assertFalse(ok)
        self.assertIn("Missing: blank.", msg)
        self.assertIn("do NOT resubmit this job", msg)
        self.assertIn("The report was written without the run(s) above (/s/report.parquet)", msg)
        self.assertIn("leave it out of conditions.csv, go on with this report", msg)
        self.assertIn("ask Core staff", msg)


if __name__ == "__main__":
    unittest.main()
