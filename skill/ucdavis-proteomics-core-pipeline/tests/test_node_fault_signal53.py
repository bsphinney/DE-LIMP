#!/usr/bin/env python3
"""
Black-hole nodes (I08): a job that died as it started -- sacct FAILED, ExitCode 0:53 (signal 53),
Elapsed a few seconds, and no log at all -- is a NODE fault that `node_fault.py retry --exclude`
acts on.

A node that cannot start jobs kills every job SLURM sends it before the job script runs, so
there is no log line to read and no exit 75 from the node check. classify() knew only the log's
I/O errors, exit 75 and NODE_FAIL, so these read as an ordinary failure: the watcher said
"unknown_failure", and retry refused them ("failures a node retry cannot fix"). A job that dies
by signal 53 AFTER writing a log is not one of these: it ran, so the log decides.

sacct's Elapsed is now read beside the other four fields (node_fault.sacct_rows, watch_run.sh
NF_FORMAT) and named in the evidence. SLURM is stubbed; nothing contacts a cluster.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
NODE_FAULT = os.path.join(SCRIPTS, "node_fault.py")
WATCH = os.path.join(SCRIPTS, "watch_run.sh")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import node_fault  # noqa: E402
# Module import, not `from ... import`: a TestCase bound here would be collected twice.
import test_node_fault as tnf  # noqa: E402  (the retry harness's fake SLURM)


class Classify(unittest.TestCase):
    def test_signal_53_without_a_log_is_a_node_fault(self):
        for log in ("", "\n  \n", "[probe] retrying the probe once\n"):
            with self.subTest(log=log):
                rec = node_fault.classify(log, "FAILED", ["0:53"], ["hive-dc-7-4-60"],
                                          elapsed=["00:00:02"])
                self.assertIsNotNone(rec)
                self.assertEqual(rec["signature"], "slurm_signal53")
                self.assertEqual(rec["meaning"], "the job died as it started (signal 53) and "
                                                 "wrote no log: a node that cannot start jobs")
                self.assertEqual(rec["evidence"], "sacct ExitCode 0:53, Elapsed 00:00:02, no log")
                self.assertEqual(rec["nodes"], ["hive-dc-7-4-60"])
                self.assertIn("hive-dc-7-4-60", rec["say"])

    def test_without_elapsed_it_says_so(self):
        rec = node_fault.classify("", "FAILED", ["0:53"], ["n1"])
        self.assertEqual(rec["evidence"], "sacct ExitCode 0:53, Elapsed not recorded, no log")

    def test_signal_53_with_a_log_is_not(self):
        self.assertIsNone(node_fault.classify("DIA-NN 2.7.0\nLoading run x.d\n", "FAILED",
                                              ["0:53"], ["n1"], elapsed=["00:41:10"]))

    def test_other_signals_and_states_are_not(self):
        self.assertIsNone(node_fault.classify("", "FAILED", ["0:9"], ["n1"]))
        self.assertIsNone(node_fault.classify("", "FAILED", ["53:0"], ["n1"]))
        self.assertIsNone(node_fault.classify("", "CANCELLED", ["0:53"], ["n1"]))
        self.assertIsNone(node_fault.classify("", "OUT_OF_MEMORY", ["0:53"], ["n1"]))

    def test_the_cli_takes_elapsed(self):
        r = subprocess.run([sys.executable, NODE_FAULT, "classify", "--state", "FAILED",
                            "--exit-codes", "0:53", "--nodes", "n7", "--elapsed", "00:00:01"],
                           input="", capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        out = json.loads(r.stdout)
        self.assertTrue(out["node_fault"])
        self.assertIn("Elapsed 00:00:01", out["evidence"])


class SacctRows(unittest.TestCase):
    def _rows(self, printed):
        with tempfile.TemporaryDirectory() as d:
            argv_log = os.path.join(d, "argv")
            tnf._exe(os.path.join(d, "sacct"), "#!/bin/bash\n"
                     f'echo "$*" > {argv_log}\n' + "".join(f'echo "{ln}"\n' for ln in printed))
            with mock.patch.dict(os.environ, {"PATH": d + os.pathsep + os.environ["PATH"]}):
                rows = node_fault.sacct_rows("103")
            return rows, open(argv_log).read()

    def test_sacct_rows_accepts_elapsed(self):
        rows, argv = self._rows(["103|FAILED|0:53|hive-dc-7-4-60|00:00:02"])
        self.assertIn("JobID,State,ExitCode,NodeList,Elapsed", argv)
        self.assertEqual(rows, [{"id": "103", "task": None, "state": "FAILED", "exit": "0:53",
                                 "node": "hive-dc-7-4-60", "elapsed": "00:00:02"}])

    def test_four_fields_still_parse(self):
        rows, _ = self._rows(["103_2|FAILED|75:0|n1"])
        self.assertEqual(rows, [{"id": "103", "task": "2", "state": "FAILED", "exit": "75:0",
                                 "node": "n1", "elapsed": None}])


class Retry(unittest.TestCase):
    """The retry harness of test_node_fault.py (a 5-step chain, jobs 101..105; borrowed, not its
    tests), with step 3 killed by a black-hole node: five sacct fields, no log."""
    setUp = tnf.RetryTests.setUp
    tearDown = tnf.RetryTests.tearDown
    set_state = tnf.RetryTests.set_state
    retry = tnf.RetryTests.retry

    def test_retry_accepts_a_signal53_task(self):
        self.set_state({
            "101": [["101", "COMPLETED", "0:0", "n0", "00:10:00"]],
            "102": [[f"102_{i}", "COMPLETED", "0:0", "n0", "01:00:00"] for i in range(6)],
            "103": [["103", "FAILED", "0:53", "hive-dc-7-4-60", "00:00:02"]],
            "104": [["104", "PENDING", "0:0", "None assigned", "00:00:00"]],
            "105": [["105", "PENDING", "0:0", "None assigned", "00:00:00"]]})
        r, out = self.retry("--dry-run")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(out["dry_run"], out)
        self.assertNotIn("cannot fix", json.dumps(out))
        self.assertEqual(out["node_fault"]["signature"], "slurm_signal53")
        self.assertIn("--exclude=hive-dc-7-4-60", out["plan"][0]["argv"])

    def test_retry_refuses_a_signal53_task_whose_log_cannot_be_read(self):
        """L8R-F5: a log there but unreadable is not "no log": not a node fault, so refused."""
        self.set_state({
            "101": [["101", "COMPLETED", "0:0", "n0", "00:10:00"]],
            "102": [[f"102_{i}", "COMPLETED", "0:0", "n0", "01:00:00"] for i in range(6)],
            "103": [["103", "FAILED", "0:53", "hive-dc-7-4-60", "00:00:02"]],
            "104": [["104", "PENDING", "0:0", "None assigned", "00:00:00"]],
            "105": [["105", "PENDING", "0:0", "None assigned", "00:00:00"]]})
        log = os.path.join(self.out, "s3_assembly_103.log")
        with open(log, "w") as fh:
            fh.write("Step 3/5 empirical library assembly\n")
        os.chmod(log, 0)
        try:
            if os.access(log, os.R_OK):
                self.skipTest("running as root: a mode-0 file is still readable")
            r, out = self.retry("--dry-run")
        finally:
            os.chmod(log, 0o644)
        self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("cannot fix", out["say"])


class LogTail(unittest.TestCase):
    def test_a_folder_that_cannot_be_read_is_none_not_empty(self):
        """L8R-F5: retry judges each task on its log; "" means no log was written, None that it
        could not be read -- which is never a black hole."""
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(node_fault._log_tail(d, "103"), "")
            self.assertIsNone(node_fault._log_tail(os.path.join(d, "missing"), "103"))
            log = os.path.join(d, "s3_assembly_103.log")
            with open(log, "w") as fh:
                fh.write("Step 3/5\n")
            os.chmod(log, 0)
            try:
                if os.access(log, os.R_OK):
                    self.skipTest("running as root: a mode-0 file is still readable")
                self.assertIsNone(node_fault._log_tail(d, "103"))
            finally:
                os.chmod(log, 0o644)


class Watcher(unittest.TestCase):
    def watch(self, rows, chain=False, out="", log_in_out=None, log=None):
        """out "": an empty search folder that exists (a job that wrote no log); a path: that
        path as given; None: no --out. log: a --log path as given."""
        with tempfile.TemporaryDirectory() as d:
            st = "\n".join(f'echo "{r}"' for r in rows)
            plain = "\n".join(f'echo "  {r.split("|")[1]}"' for r in rows)
            tnf._exe(os.path.join(d, "sacct"),
                     '#!/bin/bash\ncase "$*" in *ExitCode,NodeList,Elapsed*) '
                     f"{st} ;; *ExitCode*) echo FORMAT-WITHOUT-ELAPSED ;; *) {plain} ;; esac\n")
            tnf._exe(os.path.join(d, "squeue"), "#!/bin/bash\n")
            env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
            if chain:
                with open(os.path.join(d, "jobs.txt"), "w") as fh:
                    fh.write("103\n")
                args = ["bash", WATCH, "--all", d]
            else:
                if out == "" or log_in_out is not None:
                    out = os.path.join(d, "search out")
                    os.makedirs(out)
                if log_in_out is not None:
                    with open(os.path.join(out, "diann_search_103.log"), "w") as fh:
                        fh.write(log_in_out)
                args = (["bash", WATCH, "--slurm", "103"] + (["--out", out] if out else [])
                        + (["--log", log] if log else []))
            r = subprocess.run(args, capture_output=True, text=True, env=env, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            return json.loads(r.stdout)

    def test_the_watcher_names_a_black_hole_node(self):
        for chain in (False, True):
            with self.subTest(chain=chain):
                out = self.watch(["103|FAILED|0:53|hive-dc-7-4-60|00:00:02"], chain=chain)
                self.assertEqual(out["error_class"], "node_fault")
                self.assertEqual(out["node_fault"]["signature"], "slurm_signal53")
                self.assertIn("Elapsed 00:00:02", out["node_fault"]["evidence"])
                self.assertIn("node_fault.py retry", out["fix"])


    def test_a_job_whose_log_is_in_out_is_judged_on_it(self):
        """No --log: the job's own log in --out is read; it ran, so it is not a black hole."""
        out = self.watch(["103|FAILED|0:53|n1|00:12:00"], log_in_out="DIA-NN 2.7.0\nFile #1/2\n")
        self.assertNotEqual(out["error_class"], "node_fault")
        out = self.watch(["103|FAILED|0:53|n1|00:00:01"], log_in_out="")
        self.assertEqual(out["error_class"], "node_fault")

    def test_an_out_folder_that_cannot_be_read_is_not_a_black_hole(self):
        """L8R-F5: a wrong --out, or one that ssh failed to read under --hive, is a log NOT READ,
        never a log not written: unknown, and the watcher says the log could not be read."""
        out = self.watch(["103|FAILED|0:53|hive-dc-7-4-60|00:00:02"], out="/x/search")
        self.assertEqual(out["error_class"], "unknown_failure")
        self.assertNotIn("node_fault", out)
        self.assertIn("could not be read", out["fix"])
        self.assertIn("/x/search", out["fix"])

    def test_a_log_that_cannot_be_read_is_not_a_black_hole(self):
        """--log naming a file that is not there, and no --out to look in: not read."""
        out = self.watch(["103|FAILED|0:53|n1|00:00:01"], out=None, log="/x/search/job_103.log")
        self.assertEqual(out["error_class"], "unknown_failure")
        self.assertIn("could not be read", out["fix"])

    def test_a_log_file_that_cannot_be_read_is_not_a_black_hole(self):
        """The job's log is there but cannot be read: not read, in either mode."""
        for chain in (False, True):
            with self.subTest(chain=chain), tempfile.TemporaryDirectory() as d:
                rows = ["103|FAILED|0:53|hive-dc-7-4-60|00:00:02"]
                st = "\n".join(f'echo "{r}"' for r in rows)
                tnf._exe(os.path.join(d, "sacct"),
                         '#!/bin/bash\ncase "$*" in *ExitCode,NodeList,Elapsed*) '
                         f"{st} ;; *) echo '  FAILED' ;; esac\n")
                tnf._exe(os.path.join(d, "squeue"), "#!/bin/bash\n")
                folder = os.path.join(d, "search out")
                os.makedirs(folder)
                with open(os.path.join(folder, "jobs.txt"), "w") as fh:
                    fh.write("103\n")
                log = os.path.join(folder, "s3_assembly_103.log")
                with open(log, "w") as fh:
                    fh.write("Step 3/5 empirical library assembly\n")
                os.chmod(log, 0)
                try:
                    if os.access(log, os.R_OK):
                        self.skipTest("running as root: a mode-0 file is still readable")
                    env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
                    args = (["bash", WATCH, "--all", folder] if chain else
                            ["bash", WATCH, "--slurm", "103", "--out", folder])
                    r = subprocess.run(args, capture_output=True, text=True, env=env, timeout=120)
                finally:
                    os.chmod(log, 0o644)
                o = json.loads(r.stdout)
                self.assertNotEqual(o.get("error_class"), "node_fault", o)
                self.assertNotIn("node_fault", o)
                self.assertIn("could not be read", json.dumps(o))

    def test_with_no_log_and_no_out_it_is_not_called_a_black_hole(self):
        """Nothing says whether the job wrote a log, so signal 53 alone decides nothing."""
        out = self.watch(["103|FAILED|0:53|n1|00:00:01"], out=None)
        self.assertNotEqual(out["error_class"], "node_fault")
        self.assertEqual(out["error_class"], "unknown_failure")


if __name__ == "__main__":
    unittest.main()
