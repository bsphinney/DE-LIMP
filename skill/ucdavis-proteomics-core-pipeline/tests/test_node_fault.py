#!/usr/bin/env python3
"""
A compute node that cannot reach the storage is a NODE fault (node_fault.py), not a search
failure.

2026-09-28: the chain's step 3 died in 2 s on a node whose /quobyte mount had dropped
(`diann-linux: Transport endpoint is not connected`), the job said `FAILED: DIA-NN exited 0 but did
not write the empirical spectral library`, and steps 4-5 sat DependencyNeverSatisfied. These pin:
  * the node check every chain job runs first -- bounded, exit 75, NODE_FAULT named, DIA-NN never
    blamed;
  * classify(): storage I/O errors, exit 75 and NODE_FAIL are node faults; OOM / TIMEOUT are not,
    and neither is an ordinary failure;
  * retry: the failed step and every step waiting on it resubmitted with --exclude, at most
    MAX_NODE_RETRIES times per step, recorded; refused when it is not one bad node;
  * watch_run.sh reports error_class node_fault with the retry command.
SLURM is stubbed (sacct/sbatch/scancel/squeue on PATH); nothing contacts a cluster.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
NODE_FAULT = os.path.join(SCRIPTS, "node_fault.py")
WATCH = os.path.join(SCRIPTS, "watch_run.sh")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
import node_fault  # noqa: E402
from job_env import job_env  # noqa: E402

ENOTCONN_LOG = ("Step 3/5 empirical library assembly\n"
                "/var/spool/slurm/job24176598/slurm_script: line 16: /quobyte/proteomics-grp/"
                "dia-nn/build_270/diann-2.7.0/diann-linux: Transport endpoint is not connected\n"
                "FAILED: DIA-NN exited 0 but did not write the empirical spectral library\n")


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)


def _read(path):
    with open(path) as fh:
        return fh.read()


class ClassifyTests(unittest.TestCase):
    def test_a_dropped_mount_is_a_node_fault_not_diann(self):
        r = node_fault.classify(ENOTCONN_LOG, "FAILED", ["1:0"], ["hive-dc-7-4-50"])
        self.assertEqual(r["signature"], "enotconn")
        self.assertIn("Transport endpoint is not connected", r["evidence"])
        self.assertIn("hive-dc-7-4-50", r["say"])
        self.assertIn("not a DIA-NN or data problem", r["say"])

    def test_the_node_check_line_is_the_evidence(self):
        log = ("NODE_FAULT: node=hive-dc-7-4-50 check=the DIA-NN engine path=/q/diann-linux "
               "detail=ls: cannot access '/q/diann-linux': Transport endpoint is not connected\n")
        r = node_fault.classify(log, "FAILED", ["75:0"], ["hive-dc-7-4-50"])
        self.assertEqual(r["signature"], "preflight")
        self.assertIn("check=the DIA-NN engine", r["evidence"])

    def test_exit_75_needs_no_log(self):
        r = node_fault.classify("", "FAILED", ["75:0"], ["n1"])
        self.assertEqual(r["signature"], "preflight")

    def test_slurm_node_fail_needs_no_log(self):
        self.assertEqual(node_fault.classify("", "NODE_FAIL", ["0:0"], ["n1"])["signature"],
                         "slurm_node_fail")

    def test_an_estale_storm_and_eio(self):
        self.assertEqual(node_fault.classify("x: Stale file handle\n" * 4, "FAILED")["count"], 4)
        self.assertEqual(node_fault.classify("read: Input/output error", "FAILED")["signature"],
                         "eio")

    def test_oom_timeout_and_cancel_keep_their_own_class(self):
        for st in ("OUT_OF_MEMORY", "TIMEOUT", "CANCELLED by 123"):
            self.assertIsNone(node_fault.classify(ENOTCONN_LOG, st, ["0:125"]), st)

    def test_an_ordinary_failure_is_not_a_node_fault(self):
        log = "ERROR: no spectra\nFAILED: DIA-NN exited 0 but did not write the report\n"
        self.assertIsNone(node_fault.classify(log, "FAILED", ["1:0"], ["n1"]))

    def test_what_the_probe_absorbed_is_not_the_jobs_failure(self):
        """Step 1b's probe reread its log on ESTALE and went on; the job then failed because a
        run logged no radius. That is not a node fault (review of 2.10)."""
        log = ("[probe_window] the probe log could not be read (OSError: [Errno 116] Stale file "
               "handle); this run is abandoned\n"
               "[probe 2] [0:05] Loading run x.d\n"
               "FAILED: step 1b measured no scan-window radius -- the reason and each probe's "
               "DIA-NN log tail are above\n")
        self.assertIsNone(node_fault.classify(log, "FAILED", ["6:0"], ["n1"]))
        self.assertEqual(node_fault.classify("x: Stale file handle\n", "FAILED", ["1:0"],
                                             ["n1"])["signature"], "estale")

    def test_the_cli(self):
        r = subprocess.run([sys.executable, NODE_FAULT, "classify", "--state", "FAILED",
                            "--exit-codes", "1:0", "--nodes", "n1"], input=ENOTCONN_LOG,
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue(json.loads(r.stdout)["node_fault"])

    def test_the_job_end_post_names_it(self):
        import notify_slack
        st, why = notify_slack.classify(node_fault.NODE_FAULT_EXIT)
        self.assertEqual(st, "failed")
        self.assertIn("node problem", why)


class NodeCheckTests(unittest.TestCase):
    def run_check(self, checks, timeout_s=5, path_prefix=None):
        env = dict(os.environ)
        if path_prefix:
            env["PATH"] = path_prefix + os.pathsep + env["PATH"]
        script = "\n".join(node_fault.preflight_lines(checks, timeout_s=timeout_s)
                           + ['echo "work ran"'])
        return subprocess.run(["bash", "-c", script], capture_output=True, text=True, env=env,
                              timeout=120)

    def test_a_healthy_node_goes_on(self):
        with tempfile.TemporaryDirectory() as d:
            r = self.run_check([("-w", d, "the search folder")])
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("node check:", r.stdout)
        self.assertIn("work ran", r.stdout)

    def test_a_path_this_node_cannot_see_is_exit_75(self):
        with tempfile.TemporaryDirectory() as d:
            gone = os.path.join(d, "diann-linux")
            r = self.run_check([("-w", d, "the search folder"),
                                ("-x", gone, "the DIA-NN engine")])
        self.assertEqual(r.returncode, node_fault.NODE_FAULT_EXIT)
        self.assertNotIn("work ran", r.stdout)
        self.assertRegex(r.stderr, r"NODE_FAULT: node=\S+ check=the DIA-NN engine path=" + re.escape(gone))
        self.assertIn("not with DIA-NN or the data", r.stderr)

    def test_a_hung_mount_is_bounded(self):
        with tempfile.TemporaryDirectory() as d:
            _exe(os.path.join(d, "ls"), "#!/bin/sh\nsleep 30\n")      # a stat that never answers
            r = self.run_check([("-r", d, "the raw data")], timeout_s=2, path_prefix=d)
        self.assertEqual(r.returncode, node_fault.NODE_FAULT_EXIT)
        self.assertIn("no answer in 2 s (a hung mount)", r.stderr)

    def test_chain_checks_take_only_paths_that_exist(self):
        with tempfile.TemporaryDirectory() as d:
            binary = os.path.join(d, "diann-linux")
            _exe(binary, "#!/bin/sh\n")
            raw = os.path.join(d, "raw", "a.d")
            os.makedirs(raw)
            fasta = os.path.join(d, "db.fasta")
            _exe(fasta, ">x\nK\n")
            checks = node_fault.chain_checks(f"{binary} --verbose 1 /not/there", d, fasta, [raw])
        self.assertEqual([c[2] for c in checks], ["the DIA-NN engine", "the search folder",
                                                  "the FASTA", "the raw data"])
        self.assertEqual(checks[0][0], "-x")


class ChainJobTests(unittest.TestCase):
    """The generated chain runs the node check before DIA-NN, inside the job-end hook."""

    def test_step3_on_a_node_without_the_binary_fails_as_a_node_fault(self):
        with tempfile.TemporaryDirectory() as d:
            binary = os.path.join(d, "bin", "diann-linux")
            os.makedirs(os.path.dirname(binary))
            _exe(binary, "#!/bin/sh\nexit 0\n")
            raws = []
            for i in range(6):
                raws.append(os.path.join(d, "raw", f"r{i}.d"))
                os.makedirs(raws[-1])
            cfg = os.path.join(d, "p.cfg")
            with open(cfg, "w") as fh:
                fh.write("--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n")
            fasta = os.path.join(d, "db.fasta")
            with open(fasta, "w") as fh:
                fh.write(">sp|P1|X\nPEPTIDEK\n")
            out = os.path.join(d, "out")
            env = job_env(d, PATH=os.environ["PATH"])
            g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", binary, "--raw", *raws, "--fasta", fasta,
                                "--out", out, "--cfg", cfg, "--no-probe-window",
                                "--partition", "low", "--account", "publicgrp"],
                               capture_output=True, text=True, env=env, timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            s3 = os.path.join(out, "step3_assembly.sbatch")
            text = _read(s3)
            self.assertLess(text.index("_nf_check"), text.index("Step 3/5"))
            self.assertIn(f"_nf_check -x {binary}", text)
            os.remove(binary)                               # this node cannot see it
            r = subprocess.run(["bash", s3], capture_output=True, text=True, env=env,
                               cwd=out, timeout=120)
        self.assertEqual(r.returncode, node_fault.NODE_FAULT_EXIT, r.stdout + r.stderr)
        self.assertIn("NODE_FAULT: node=", r.stderr)
        self.assertNotIn("DIA-NN exited 0", r.stdout + r.stderr)

    def test_single_shot_jobs_check_their_node_first_too(self):
        import run_search
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "job.sh")
            run_search.emit_sbatch(job, "diann-linux --threads 4", d, 4, job="diann_search",
                                   partition="low", account="publicgrp", qos="publicgrp-low-qos",
                                   preflight=node_fault.preflight_lines([("-w", d, "x")]))
            text = _read(job)
        self.assertLess(text.index("_nf_check -w"), text.index("diann-linux --threads 4"))
        self.assertLess(text.index("#SBATCH --qos"), text.index("_nf_check -w"))


class Step5WritesNothingOnAnIncompleteFinalPass(unittest.TestCase):
    """Step 5 counted the final-pass .quant files only AFTER DIA-NN had written report.parquet:
    a chain resubmitted without one step-4 task (review of 2.10) left a report from N-1 runs at
    the path run_de.R reads, and then failed. It counts first now."""

    def chain(self, d):
        diann = os.path.join(d, "fake-diann")
        _exe(diann, "#!/bin/sh\nwhile [ $# -gt 0 ]; do [ \"$1\" = --out ] && echo report > \"$2\"; "
                    "shift; done\n")
        raws = []
        for i in range(3):
            raws.append(os.path.join(d, "raw", f"r{i}.d"))
            os.makedirs(raws[-1])
        cfg = os.path.join(d, "p.cfg")
        with open(cfg, "w") as fh:
            fh.write("--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n")
        fasta = os.path.join(d, "db.fasta")
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEK\n")
        out = os.path.join(d, "out")
        env = job_env(d, PATH=os.environ["PATH"])
        g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", diann, "--raw", *raws, "--fasta", fasta, "--out", out,
                            "--cfg", cfg, "--no-probe-window", "--partition", "low",
                            "--account", "publicgrp"],
                           capture_output=True, text=True, env=env, timeout=120)
        self.assertEqual(g.returncode, 0, g.stderr)
        os.makedirs(os.path.join(out, "quant_step4"), exist_ok=True)
        with open(os.path.join(out, "empirical.parquet"), "w") as fh:
            fh.write("lib")
        return out, env

    def run_step5(self, out, env, quants):
        for q in quants:
            with open(os.path.join(out, "quant_step4", f"{q}.quant"), "w") as fh:
                fh.write("q")
        return subprocess.run(["bash", os.path.join(out, "step5_report.sbatch")],
                              capture_output=True, text=True, env=env, cwd=out, timeout=120)

    def test_a_missing_final_pass_quant_means_no_report(self):
        with tempfile.TemporaryDirectory() as d:
            out, env = self.chain(d)
            p = self.run_step5(out, env, ["r0", "r1"])
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("MISSING: r2.quant", p.stderr)
            self.assertIn("no report was written", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "report.parquet")))

    def test_all_runs_present_builds_the_report(self):
        with tempfile.TemporaryDirectory() as d:
            out, env = self.chain(d)
            p = self.run_step5(out, env, ["r0", "r1", "r2"])
            self.assertIn("OK: report built from all 3 runs", p.stdout)
            self.assertTrue(os.path.exists(os.path.join(out, "report.parquet")))


# ------------------------------------------------------------------------------------- retry
SACCT = r"""#!/usr/bin/env python3
import json, os, sys
st = json.load(open(os.environ["NF_STATE"]))
job = sys.argv[sys.argv.index("-j") + 1]
for row in st.get(job, []):
    print("|".join(row))
"""
SBATCH = r"""#!/usr/bin/env python3
import json, os, sys
log = os.environ["NF_LOG"]
n = sum(1 for _ in open(log)) if os.path.exists(log) else 0
with open(log, "a") as fh:
    fh.write(json.dumps({"cmd": "sbatch", "argv": sys.argv[1:], "cwd": os.getcwd()}) + "\n")
print(201 + n)
"""
SCANCEL = r"""#!/usr/bin/env python3
import json, os, sys
with open(os.environ["NF_LOG"] + ".scancel", "a") as fh:
    fh.write(sys.argv[1] + "\n")
"""


class RetryTests(unittest.TestCase):
    """A 5-step chain (no step 1b): jobs 101..105."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        for name, body in (("sacct", SACCT), ("sbatch", SBATCH), ("scancel", SCANCEL)):
            _exe(os.path.join(self.bin, name), body)
        self.sess = os.path.join(self.d, "2026-10-01_Test")
        self.out = os.path.join(self.sess, "output", "search")
        os.makedirs(self.out)
        with open(os.path.join(self.out, "submit.sh"), "w") as fh:
            fh.write("#!/bin/bash\nset -euo pipefail\n"
                     f'cd "{self.out}"\n'
                     "jid1=$(sbatch --parsable step1_libpred.sbatch)\n"
                     "jid2=$(sbatch --parsable --dependency=afterok:$jid1 step2_firstpass.sbatch)\n"
                     "jid3=$(sbatch --parsable --dependency=afterok:$jid2 step3_assembly.sbatch)\n"
                     "jid4=$(sbatch --parsable --dependency=afterok:$jid3 step4_finalpass.sbatch)\n"
                     "jid5=$(sbatch --parsable --dependency=afterok:$jid4 step5_report.sbatch)\n"
                     f'printf "%s\\n" $jid1 $jid2 $jid3 $jid4 $jid5 > "{self.out}/jobs.txt"\n')
        with open(os.path.join(self.out, "jobs.txt"), "w") as fh:
            fh.write("101\n102\n103\n104\n105\n")
        with open(os.path.join(self.out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann"}, fh)
        self.state = os.path.join(self.d, "state.json")
        self.log = os.path.join(self.d, "calls.jsonl")
        self.env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ["PATH"],
                        NF_STATE=self.state, NF_LOG=self.log)

    def tearDown(self):
        self._tmp.cleanup()

    def set_state(self, st):
        with open(self.state, "w") as fh:
            json.dump(st, fh)

    def step3_failed_on(self, node, code="75:0"):
        self.set_state({
            "101": [["101", "COMPLETED", "0:0", "n0"]],
            "102": [[f"102_{i}", "COMPLETED", "0:0", "n0"] for i in range(6)],
            "103": [["103", "FAILED", code, node]],
            "104": [["104", "PENDING", "0:0", "None assigned"]],
            "105": [["105", "PENDING", "0:0", "None assigned"]]})

    def retry(self, *extra):
        r = subprocess.run([sys.executable, NODE_FAULT, "retry", "--out", self.out, *extra],
                           capture_output=True, text=True, env=self.env, timeout=120)
        return r, (json.loads(r.stdout) if r.stdout.strip().startswith("{") else None)

    def calls(self):
        if not os.path.exists(self.log):
            return []
        with open(self.log) as fh:
            return [json.loads(ln) for ln in fh]

    def test_the_failed_step_and_its_dependants_go_elsewhere(self):
        self.step3_failed_on("hive-dc-7-4-50")
        r, out = self.retry()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        argvs = [c["argv"] for c in self.calls()]
        self.assertEqual(argvs[0], ["--parsable", "--exclude=hive-dc-7-4-50",
                                    "step3_assembly.sbatch"])
        self.assertEqual(argvs[1], ["--parsable", "--dependency=afterok:201",
                                    "--exclude=hive-dc-7-4-50", "step4_finalpass.sbatch"])
        self.assertEqual(argvs[2][1], "--dependency=afterok:202")
        self.assertTrue(all(c["cwd"] == os.path.realpath(self.out) or c["cwd"] == self.out
                            for c in self.calls()))
        self.assertEqual(_read(self.log + ".scancel").split(), ["104", "105"])
        self.assertEqual(_read(os.path.join(self.out, "jobs.txt")).split(),
                         ["101", "102", "201", "202", "203"])
        self.assertTrue(out["complete"])
        self.assertIn("Node problem, not a DIA-NN or data problem", out["say"])
        rec = json.loads(_read(os.path.join(self.out, "node_faults.json")))["retries"]
        self.assertEqual((rec[0]["attempt"], rec[0]["nodes"]), (1, ["hive-dc-7-4-50"]))
        prov = json.loads(_read(os.path.join(self.out, "search_provenance.json")))
        self.assertEqual(prov["node_faults"][0]["failed_job"], "103")

    def test_bounded_per_step_and_every_bad_node_stays_excluded(self):
        self.step3_failed_on("nA")
        self.assertEqual(self.retry()[0].returncode, 0)
        # the retried step 3 (201) fails on another node
        with open(os.path.join(self.out, "jobs.txt"), "w") as fh:
            fh.write("101\n102\n203\n204\n205\n")
        st = json.loads(_read(self.state))
        st["203"] = [["203", "FAILED", "75:0", "nB"]]
        st["204"] = st["205"] = [["x", "PENDING", "0:0", "None assigned"]]
        self.set_state(st)
        r, out = self.retry()
        self.assertEqual(r.returncode, 0, r.stdout)
        self.assertIn("--exclude=nA,nB", self.calls()[-3]["argv"])
        with open(os.path.join(self.out, "jobs.txt"), "w") as fh:
            fh.write("101\n102\n303\n304\n305\n")
        st["303"] = [["303", "FAILED", "75:0", "nC"]]
        st["304"] = st["305"] = [["x", "PENDING", "0:0", "None assigned"]]
        self.set_state(st)
        n = len(self.calls())
        r, out = self.retry()
        self.assertEqual(r.returncode, 3)
        self.assertIn("not one bad node", out["say"])
        self.assertEqual(len(self.calls()), n)                    # nothing submitted

    def test_an_array_retries_only_its_failed_tasks(self):
        self.set_state({
            "101": [["101", "COMPLETED", "0:0", "n0"]],
            "102": [[f"102_{i}", "COMPLETED", "0:0", "n0"] for i in (0, 1, 2, 4, 5)]
                   + [["102_3", "FAILED", "75:0", "hive-dc-7-4-54"]],
            "103": [["103", "PENDING", "0:0", "None assigned"]],
            "104": [["104", "PENDING", "0:0", "None assigned"]],
            "105": [["105", "PENDING", "0:0", "None assigned"]]})
        r, out = self.retry("--job", "102_3")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.calls()[0]["argv"], ["--parsable", "--exclude=hive-dc-7-4-54",
                                                   "--array=3", "step2_firstpass.sbatch"])
        self.assertEqual(_read(self.log + ".scancel").split(), ["103", "104", "105"])

    def array4(self, task4_state, task4_exit):
        """Step 2 an array of 6: task 3 a node fault on nX, task 4 as given on n2."""
        self.set_state({
            "101": [["101", "COMPLETED", "0:0", "n0"]],
            "102": [[f"102_{i}", "COMPLETED", "0:0", "n0"] for i in (0, 1, 2, 5)]
                   + [["102_3", "FAILED", "75:0", "nX"], ["102_4", task4_state, task4_exit, "n2"]],
            "103": [["103", "PENDING", "0:0", "None assigned"]],
            "104": [["104", "PENDING", "0:0", "None assigned"]],
            "105": [["105", "PENDING", "0:0", "None assigned"]]})

    def test_an_oom_task_beside_a_node_fault_is_said_never_left_out(self):
        """Review of 2.10: the OOM task was neither retried nor mentioned, step 5 was resubmitted
        afterok on the node-fault task alone, and the retry said complete."""
        self.array4("OUT_OF_MEMORY", "0:125")
        r, out = self.retry()
        self.assertEqual(r.returncode, 3)
        self.assertEqual(self.calls(), [])                       # nothing resubmitted
        self.assertFalse(out.get("complete"))
        oom = out["not_node_faults"]
        self.assertEqual((oom[0]["task"], oom[0]["state"]), ("4", "OUT_OF_MEMORY"))
        self.assertIn("raise the memory", oom[0]["fix"])
        self.assertEqual([t["task"] for t in out["node_faults"]], ["3"])

    def test_each_task_is_judged_on_its_own_log(self):
        """The first task's NODE_FAULT once stood for all: task 4's DIA-NN error got its node
        excluded and blamed."""
        self.array4("FAILED", "1:0")
        with open(os.path.join(self.out, "s2_firstpass_102_4.log"), "w") as fh:
            fh.write("ERROR: cannot load file: run4.raw\n")
        with open(os.path.join(self.out, "s2_firstpass_102_3.log"), "w") as fh:
            fh.write("NODE_FAULT: node=nX check=the raw data path=/q/x detail=timeout\n")
        r, out = self.retry()
        self.assertEqual(r.returncode, 3)
        self.assertEqual(self.calls(), [])
        self.assertEqual(out["not_node_faults"][0]["node"], "n2")
        self.assertEqual(out["node_faults"][0]["signature"], "preflight")

    def test_two_node_fault_tasks_are_both_retried_and_complete(self):
        self.array4("FAILED", "75:0")
        r, out = self.retry()
        self.assertEqual(r.returncode, 0, r.stdout)
        self.assertEqual(self.calls()[0]["argv"], ["--parsable", "--exclude=n2,nX",
                                                   "--array=3,4", "step2_firstpass.sbatch"])
        self.assertTrue(out["complete"])
        rec = json.loads(_read(os.path.join(self.out, "node_faults.json")))["retries"][0]
        self.assertEqual(sorted(rec["evidence_per_task"]), ["3", "4"])

    def test_a_retry_holds_a_lock_and_leaves_none(self):
        self.step3_failed_on("n1")
        lock = os.path.join(self.out, node_fault.LOCK)
        os.mkdir(lock)                                     # another retry is running
        r, out = self.retry_env(NODE_FAULT_LOCK_WAIT_S="1")
        self.assertEqual(r.returncode, 3)
        self.assertIn("another node_fault.py retry", out["say"])
        self.assertEqual(self.calls(), [])
        old = time.time() - node_fault.LOCK_STALE_S - 60  # ...or a dead one's, long ago
        os.utime(lock, (old, old))
        r, out = self.retry()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertFalse(os.path.exists(lock))

    def retry_env(self, **extra):
        env = dict(self.env, **extra)
        r = subprocess.run([sys.executable, NODE_FAULT, "retry", "--out", self.out],
                           capture_output=True, text=True, env=env, timeout=120)
        return r, (json.loads(r.stdout) if r.stdout.strip().startswith("{") else None)

    def test_an_array_still_running_is_waited_for(self):
        self.set_state({"101": [["101", "COMPLETED", "0:0", "n0"]],
                        "102": [["102_3", "FAILED", "75:0", "nX"], ["102_4", "RUNNING", "0:0", "n1"]],
                        "103": [], "104": [], "105": []})
        r, out = self.retry()
        self.assertEqual(r.returncode, 3)
        self.assertIn("still running", out["say"])
        self.assertEqual(self.calls(), [])

    def test_an_ordinary_failure_is_refused_unless_forced(self):
        self.step3_failed_on("n1", code="1:0")
        with open(os.path.join(self.out, "s3_assembly_103.log"), "w") as fh:
            fh.write("FAILED: DIA-NN exited 0 but did not write the empirical spectral library\n")
        r, out = self.retry()
        self.assertEqual(r.returncode, 3)
        self.assertIn("a node retry cannot fix", out["say"])
        self.assertEqual(out["not_node_faults"][0]["state"], "FAILED")
        self.assertEqual(self.calls(), [])
        r, out = self.retry("--force")
        self.assertEqual(r.returncode, 0, r.stdout)

    def test_the_log_evidence_is_read(self):
        self.step3_failed_on("hive-dc-7-4-50", code="1:0")
        with open(os.path.join(self.out, "s3_assembly_103.log"), "w") as fh:
            fh.write(ENOTCONN_LOG)
        r, out = self.retry()
        self.assertEqual(r.returncode, 0, r.stdout)
        self.assertEqual(out["node_fault"]["signature"], "enotconn")

    def test_failures_on_many_nodes_are_the_storage(self):
        self.set_state({"101": [["101", "COMPLETED", "0:0", "n0"]],
                        "102": [[f"102_{i}", "FAILED", "75:0", f"n{i}"] for i in range(5)],
                        "103": [], "104": [], "105": []})
        r, out = self.retry()
        self.assertEqual(r.returncode, 3)
        self.assertIn("not one bad node", out["say"])

    def test_a_dependant_already_running_is_refused(self):
        self.step3_failed_on("n1")
        st = json.loads(_read(self.state))
        st["104"] = [["104", "RUNNING", "0:0", "n2"]]
        self.set_state(st)
        r, out = self.retry()
        self.assertEqual(r.returncode, 3)
        self.assertEqual(self.calls(), [])

    def test_dry_run_submits_nothing(self):
        self.step3_failed_on("n1")
        r, out = self.retry("--dry-run")
        self.assertEqual(r.returncode, 0)
        self.assertEqual(len(out["plan"]), 3)
        self.assertEqual(self.calls(), [])
        self.assertFalse(os.path.exists(os.path.join(self.out, "node_faults.json")))

    def test_the_recovery_checkpoint_follows_the_new_ids(self):
        import checkpoint
        checkpoint._save(self.sess, {"session": self.sess, "stages": [{
            "stage": "search", "jobs": ["101", "102", "103", "104", "105"], "desc": "chain",
            "state": "submitted", "next": None, "watch_job": "105",
            "watch_log": f"{self.out}/s5_report_105.log"}]})
        self.step3_failed_on("n1")
        r, out = self.retry()
        self.assertEqual(r.returncode, 0, r.stdout)
        st = json.loads(_read(os.path.join(self.sess, ".recovery.json")))["stages"][0]
        self.assertEqual(st["jobs"], ["101", "102", "201", "202", "203"])
        self.assertEqual(st["watch_job"], "203")
        self.assertTrue(st["watch_log"].endswith("s5_report_203.log"))

    def test_the_two_job_search_is_understood_too(self):
        lib, srch = os.path.join(self.d, "job_1_lib.sh"), os.path.join(self.d, "job_2_search.sh")
        with open(os.path.join(self.out, "submit.sh"), "w") as fh:
            fh.write("#!/bin/bash -l\nset -euo pipefail\n"
                     f"j1=$(sbatch --parsable {lib})\n"
                     f"j2=$(sbatch --parsable --dependency=afterok:$j1 {srch})\n"
                     f'printf "%s\\n" "$j1" "$j2" > {self.out}/jobs.txt\n')
        with open(os.path.join(self.out, "jobs.txt"), "w") as fh:
            fh.write("101\n102\n")
        self.set_state({"101": [["101", "FAILED", "75:0", "n9"]],
                        "102": [["102", "PENDING", "0:0", "None assigned"]]})
        r, out = self.retry()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.calls()[0]["argv"][-1], lib)
        self.assertEqual(self.calls()[1]["argv"][-1], srch)


class WatchRunTests(unittest.TestCase):
    def watch(self, rows, tail="", chain=False):
        with tempfile.TemporaryDirectory() as d:
            st = "\n".join(f'echo "{r}"' for r in rows)
            plain = "\n".join(f'echo "  {r.split("|")[1]}"' for r in rows)
            _exe(os.path.join(d, "sacct"),
                 f'#!/bin/bash\ncase "$*" in *ExitCode*) {st} ;; *) {plain} ;; esac\n')
            _exe(os.path.join(d, "squeue"), "#!/bin/bash\n")
            env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
            if chain:
                with open(os.path.join(d, "jobs.txt"), "w") as fh:
                    fh.write("103\n")
                if tail:
                    with open(os.path.join(d, "s3_assembly_103.log"), "w") as fh:
                        fh.write(tail)
                args = ["bash", WATCH, "--all", d]
            else:
                log = os.path.join(d, "job.log")
                with open(log, "w") as fh:
                    fh.write(tail)
                args = ["bash", WATCH, "--slurm", "103", "--log", log, "--out", "/x/search"]
            r = subprocess.run(args, capture_output=True, text=True, env=env, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            return json.loads(r.stdout)

    def test_slurm_mode_names_the_node_and_the_retry(self):
        out = self.watch(["103|FAILED|75:0|hive-dc-7-4-50"])
        self.assertEqual(out["error_class"], "node_fault")
        self.assertIn("hive-dc-7-4-50", out["fix"])
        self.assertIn("node_fault.py retry --out /x/search --job 103", out["fix"])
        self.assertEqual(out["node_fault"]["nodes"], ["hive-dc-7-4-50"])

    def test_a_dropped_mount_in_the_log(self):
        out = self.watch(["103|FAILED|1:0|hive-dc-7-4-50"], tail=ENOTCONN_LOG)
        self.assertEqual(out["error_class"], "node_fault")
        self.assertNotIn("DIA-NN exited 0", out["fix"])

    def test_oom_stays_oom(self):
        out = self.watch(["103|OUT_OF_MEMORY|0:125|n1"], tail="Stale file handle\n")
        self.assertEqual(out["error_class"], "out_of_memory")

    def test_chain_mode_reads_the_failed_tasks_log_not_the_first(self):
        """Task 0 of the array succeeded; task 3 lost its mount. The first task's log (all a
        `ls | head -1` saw) says nothing."""
        with tempfile.TemporaryDirectory() as d:
            folder = os.path.join(d, "search out")                 # a space, too
            os.makedirs(folder)
            rows = ["104_0|COMPLETED|0:0|n0", "104_3|FAILED|1:0|hive-dc-7-4-50"]
            _exe(os.path.join(d, "sacct"), "#!/bin/bash\ncase \"$*\" in *ExitCode*) "
                 + " ".join(f'echo "{r}";' for r in rows) + " ;; *) echo '  COMPLETED'; "
                 "echo '  FAILED' ;; esac\n")
            _exe(os.path.join(d, "squeue"), "#!/bin/bash\n")
            with open(os.path.join(folder, "jobs.txt"), "w") as fh:
                fh.write("104\n")
            with open(os.path.join(folder, "s4_finalpass_104_0.log"), "w") as fh:
                fh.write("Step 4/5 final pass, task 0 of 6\nall fine\n")
            with open(os.path.join(folder, "s4_finalpass_104_3.log"), "w") as fh:
                fh.write(ENOTCONN_LOG)
            env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
            r = subprocess.run(["bash", WATCH, "--all", folder], capture_output=True, text=True,
                               env=env, timeout=120)
            out = json.loads(r.stdout)
        self.assertEqual(out["error_class"], "node_fault")
        self.assertEqual(out["node_fault"]["signature"], "enotconn")
        self.assertIn("search\\ out", out["fix"])                  # printf %q: one word

    def test_chain_mode(self):
        out = self.watch(["103|FAILED|75:0|hive-dc-7-4-50"], chain=True)
        self.assertEqual(out["error_class"], "node_fault")
        self.assertIn("node_fault.py retry", out["fix"])
        self.assertEqual(out["first_failed_job"], "103")


if __name__ == "__main__":
    unittest.main()
