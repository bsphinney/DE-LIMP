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
        self.assertIn("--accept-missing <run> --reason", msg)
        self.assertIn("leave the run out of conditions.csv and go on with this report", msg)
        self.assertIn("ask Core staff", msg)


@unittest.skipUnless(HAVE_PYARROW, "the report is read with pyarrow")
class FragPipeDiaReportRuns(unittest.TestCase):
    """Review of 2.12 (P1): a run FragPipe's DIA-NN step left out fails the adapt, as on the
    DIA-NN routes (check_report_runs.verify)."""

    def out_with(self, d, runs_in_report, manifest=True):
        import run_search  # noqa: F401
        out = os.path.join(d, "fp_out")
        dq = os.path.join(out, "group", "dia-quant-output")
        os.makedirs(dq)
        pq.write_table(pa.table({"Run": runs_in_report, "Protein.Group": ["P1"] * len(runs_in_report)}),
                       os.path.join(dq, "report.parquet"))
        if manifest:
            stage = os.path.join(out, "diatracer_stage")
            write(os.path.join(out, "fragpipe.fp-manifest"),
                  f"{stage}/s1_diatracer.mzML\t\t\tDDA\n"
                  + "".join(f"{stage}/s{i}.d\t\t\t{'DIA-Quant' if i == 1 else 'DIA'}\n"
                            for i in (1, 2, 3)))
        return out

    def test_a_missing_run_fails_and_a_complete_report_passes(self):
        import run_search
        with tempfile.TemporaryDirectory() as d:
            out = self.out_with(d, ["s1", "s3"])
            with self.assertRaises(SystemExit) as cm:
                run_search.adapt_fragpipe_dia(out)
            self.assertIn("FragPipe DIA: FAILED: the report holds 2 of 3 runs", str(cm.exception))
            self.assertIn("Missing: s2.", str(cm.exception))
        with tempfile.TemporaryDirectory() as d:
            out = self.out_with(d, ["s1", "s2", "s3"])
            self.assertTrue(run_search.adapt_fragpipe_dia(out).endswith("report.parquet"))
        with tempfile.TemporaryDirectory() as d:       # no manifest of ours: said, not checked
            out = self.out_with(d, ["s1"], manifest=False)
            self.assertTrue(run_search.adapt_fragpipe_dia(out).endswith("report.parquet"))


@unittest.skipUnless(HAVE_PYARROW, "the report is read with pyarrow")
class AcceptMissing(unittest.TestCase):
    """Review of 2.12 (P5): dropping a blank from file_list.txt and resubmitting step 5 does NOT
    work -- the step's --f list and run count are fixed at generation (tested below) -- so the
    report already written, without the blank, is accepted on the record instead."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.d = self._td.name
        self.report = os.path.join(self.d, "report.parquet")
        pq.write_table(pa.table({"Run": ["s0", "s2"]}), self.report)
        self.listing = write(os.path.join(self.d, "file_list.txt"),
                             "".join(f"/raw/s{i}.d\n" for i in range(3)))

    def tearDown(self):
        self._td.cleanup()

    def crr(self, *more):
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "check_report_runs.py"),
                               "--report", self.report, "--files-list", self.listing, *more],
                              capture_output=True, text=True)

    def test_the_failure_gives_the_command_and_accepting_records_it(self):
        p = self.crr()
        self.assertEqual(p.returncode, 1)
        self.assertIn("do NOT resubmit this job, and do not edit the file list", p.stderr)
        self.assertIn(f"--report {self.report} --files-list {self.listing} --accept-missing <run>",
                      p.stderr)
        p = self.crr("--accept-missing", "s1", "--reason", "wash injection, staff agreed")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("left out by decision: s1", p.stdout)
        with open(crr.accepted_path(self.report)) as fh:
            rec = json.load(fh)
        self.assertEqual((rec["runs"], rec["reason"]), (["s1"], "wash injection, staff agreed"))

    def test_only_a_missing_input_can_be_accepted_and_only_with_a_reason(self):
        self.assertEqual(self.crr("--accept-missing", "s1").returncode, 2)
        p = self.crr("--accept-missing", "s0", "--reason", "x")
        self.assertEqual(p.returncode, 1)
        self.assertIn("in the report: s0", p.stderr)
        p = self.crr("--accept-missing", "s9", "--reason", "x")
        self.assertIn("not an input: s9", p.stderr)
        self.assertFalse(os.path.exists(crr.accepted_path(self.report)))

    def test_editing_file_list_and_resubmitting_step5_fails(self):
        from job_env import job_env
        from test_single_shot_sbatch import LIBFREE_CFG, _Workspace
        w = _Workspace(self.d, cfg_lines=LIBFREE_CFG + ["--window 7"], n_runs=3, ext=".d")
        out = os.path.join(self.d, "chain")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", w.diann, "--fasta", w.fasta, "--out", out, "--cfg", w.cfg,
                            "--partition", "high", "--account", "genome-center-grp",
                            "--raw", *w.runs], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        os.makedirs(os.path.join(out, "quant_step4"))
        for i in range(3):
            write(os.path.join(out, "quant_step4", f"sample_{i}.quant"), "q")
        fl = os.path.join(out, "file_list.txt")
        with open(fl) as fh:
            kept = [x for x in fh.read().split() if "sample_1" not in x]
        write(fl, "\n".join(kept) + "\n")
        env = job_env(self.d, FAKE_PY=sys.executable, FAKE_REPORT_RUNS="sample_0,sample_2")
        r = subprocess.run(["bash", os.path.join(out, "step5_report.sbatch")], cwd=out, env=env,
                           capture_output=True, text=True)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("3 inputs map to only 2 distinct run names", r.stderr)
        self.assertIn("check_report_runs.py --accept-missing", r.stderr)


class ReplayAndChecks(unittest.TestCase):
    def test_reproduce_sh_replays_rename_dashed_runs(self):
        """Review of 2.12 (P6)."""
        with tempfile.TemporaryDirectory() as d:
            search = os.path.join(d, "search")
            write(os.path.join(search, "report.parquet"))
            for renamed, want in ((True, True), (False, False)):
                write(os.path.join(search, "search_provenance.json"), json.dumps(
                    {"engine": "diann", "input_links": [{"original": f"/r/{DASHED}.d",
                                                         "link": f"/o/input_links/{RENAMED}.d",
                                                         "run_name": RENAMED,
                                                         "renamed": renamed}]}))
                out = os.path.join(d, f"repro_{renamed}")
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                                    "--outdir", out, "--engine", "diann", "--acquisition", "DIA",
                                    "--report", os.path.join(search, "report.parquet")],
                                   capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, p.stderr)
                with open(os.path.join(out, "reproduce.sh")) as fh:
                    line = next(ln for ln in fh if "--files" in ln)
                self.assertEqual("--rename-dashed-runs" in line, want, line)

    def test_a_cfg_value_with_dashes_is_refused(self):
        """Review of 2.12 (P7): a --lib or --temp path in the cfg is read by DIA-NN like argv."""
        from test_single_shot_sbatch import HIGH, LIBFREE_CFG, _Workspace
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d, cfg_lines=LIBFREE_CFG + [f"--temp {d}/tmp--1"])
            p = w.generate(*HIGH, check=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn(f"--temp '{d}/tmp--1'", p.stderr)
            self.assertIn("reads a cfg like a command line", p.stderr)
            self.assertFalse(os.path.exists(w.out))

    def test_probe_window_by_hand_refuses_dashed_paths(self):
        """Review of 2.12 (P7): probe_window.py run by hand (with --seed-lib, as documented)."""
        with tempfile.TemporaryDirectory() as d:
            raw = write(os.path.join(d, "plate--1", "s1.mzML"))
            fasta = write(os.path.join(d, "db.fasta"), ">sp|P1|X\nPEPTIDEK\n")
            lib = write(os.path.join(d, "lib.speclib"), "lib")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "probe_window.py"),
                                "--diann", "/nonexistent/diann", "--raw", raw, "--fasta", fasta,
                                "--lib", lib], capture_output=True, text=True)
            import probe_window
            self.assertEqual(p.returncode, probe_window.EXIT_CONFIG, p.stderr)
            self.assertIn(f"run {raw}", p.stderr)
            self.assertIn("Link each under a name and folder without '--'", p.stderr)

    def test_the_node_check_looks_where_a_link_points(self):
        """Error audit of 2.12 (E1): with '--' inputs the job got links, and its node check tested
        <out>/input_links -- not the raw data's mount, which is what a node can lose."""
        import node_fault
        with tempfile.TemporaryDirectory() as d:
            real = write(os.path.join(d, "raw--data", "s1.d"))
            link = os.path.join(d, "out", "input_links", "s1.d")
            os.makedirs(os.path.dirname(link))
            os.symlink(real, link)
            checks = node_fault.chain_checks("/x/diann", os.path.join(d, "out"), real, [link])
            raw_dirs = [p for f, p, what in checks if what == "the raw data"]
            self.assertEqual(raw_dirs, [os.path.dirname(os.path.realpath(real))])


if __name__ == "__main__":
    unittest.main()
