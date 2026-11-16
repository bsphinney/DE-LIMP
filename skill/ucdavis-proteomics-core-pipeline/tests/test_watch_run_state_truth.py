#!/usr/bin/env python3
"""
watch_run.sh tells the truth about the Slurm state of a job (2.15 acceptance, lane X4).

  * a job sacct reports FAILED / OUT_OF_MEMORY / TIMEOUT / CANCELLED is reported as that, with the
    exit code and the log's tail, and its one-line progress summary never says "search running";
  * `--out <dir>` alone (no --slurm, no --all) finds the job ids the search folder holds
    (jobs.txt, <step>_<jobid>.log) and reads their state from sacct -- it used to read none, and a
    FAILED job passed for a running one (an agent's watcher never woke);
  * a Sage "fatal runtime error: stack overflow" (exit 134), or the job's own sage_failure.json,
    is the named class sage_stack_overflow with a re-run hint, never "out of memory";
  * a normalisation check that stopped to ask a question (its job exits 0 and leaves
    NEEDS_DECISION.json) is a COMPLETED job that says a decision is owed -- not a failed one.

A fake sacct/squeue stands in for SLURM; nothing contacts a cluster.
"""
import json
import os
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
WATCH = os.path.join(os.path.dirname(HERE), "scripts", "watch_run.sh")

FAKE_SACCT = r'''#!/bin/bash
# -P: JobID|State|ExitCode|NodeList|Elapsed rows (the node-fault query and the exit-code query);
# otherwise the bare State column, as `sacct -X --noheader -o State` prints it.
STATE="${FAKE_STATE:-FAILED}"; CODE="${FAKE_EXIT:-1:0}"
case "$*" in
  *-P*) echo "${FAKE_JOB:-9900001}|$STATE|$CODE|hive-node-1|01:02:03" ;;
  *)    echo "  $STATE" ;;
esac
'''
FAKE_SQUEUE = '#!/bin/bash\ncase "$*" in *%R*) echo "${FAKE_REASON:-}" ;; *) echo "$FAKE_STATE" ;; esac\n'


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.bin = os.path.join(self.tmp.name, "bin")
        os.makedirs(self.bin)
        for name, body in (("sacct", FAKE_SACCT), ("squeue", FAKE_SQUEUE)):
            p = os.path.join(self.bin, name)
            with open(p, "w") as fh:
                fh.write(body)
            os.chmod(p, 0o755)
        self.out = os.path.join(self.tmp.name, "search_out")
        os.makedirs(self.out)

    def write(self, rel, text):
        p = os.path.join(self.out, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def watch(self, *args, state="FAILED", code="1:0", job="9900001", rc=0):
        env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ.get("PATH", ""),
                   FAKE_STATE=state, FAKE_EXIT=code, FAKE_JOB=job)
        r = subprocess.run(["bash", WATCH, *args], capture_output=True, text=True, env=env,
                           timeout=90)
        self.assertEqual(r.returncode, rc, r.stdout + r.stderr)
        return json.loads(r.stdout)


class FailedIsFailed(Base):
    def test_failed_job_is_reported_failed_with_code_and_log_tail(self):
        self.write("search_9900001.log", "starting\nSomething broke in the engine\n")
        o = self.watch("--slurm", "9900001", "--out", self.out)
        self.assertTrue(o["failed"] and o["done"])
        self.assertEqual(o["state"], "FAILED")
        self.assertEqual(o["exit_code"], "1:0")
        self.assertIn("Something broke", o["log_tail"])        # read from --out, no --log needed
        self.assertNotIn("running", o["progress"]["summary"].lower())
        self.assertIn("FAILED", o["progress"]["summary"])
        self.assertIn("1:0", o["progress"]["summary"])

    def test_each_bad_state_is_named(self):
        for state in ("OUT_OF_MEMORY", "TIMEOUT", "CANCELLED", "NODE_FAIL"):
            with self.subTest(state=state):
                o = self.watch("--slurm", "9900001", "--out", self.out, state=state, code="0:9")
                self.assertTrue(o["failed"] and o["done"], o)
                self.assertEqual(o["state"], state)
                self.assertIn(state, o["progress"]["summary"])
                self.assertNotIn("running", o["progress"]["summary"].lower())

    def test_a_running_job_still_says_running(self):
        o = self.watch("--slurm", "9900001", "--out", self.out, state="RUNNING", code="0:0")
        self.assertFalse(o["failed"] or o["done"])
        self.assertEqual(o["progress"]["summary"], "search running")

    def test_completed_is_not_called_running(self):
        o = self.watch("--slurm", "9900001", "--out", self.out, state="COMPLETED", code="0:0")
        self.assertTrue(o["done"] and not o["failed"])
        self.assertNotIn("running", o["progress"]["summary"].lower())


class OutAloneFindsTheJob(Base):
    def test_out_with_a_job_log_reads_that_jobs_state(self):
        self.write("search_9900001.log", "boom\n")
        o = self.watch("--out", self.out, state="FAILED", code="2:0")
        self.assertTrue(o["failed"] and o["done"], o)
        self.assertEqual(o["state"], "FAILED")
        self.assertEqual(o["exit_code"], "2:0")

    def test_out_with_a_log_flag_too(self):
        log = self.write("search_9900001.log", "boom\n")
        o = self.watch("--out", self.out, "--log", log, state="FAILED")
        self.assertTrue(o["failed"], o)

    def test_out_with_jobs_txt_reads_the_whole_chain(self):
        self.write("jobs.txt", "9900001\n9900002\n")
        o = self.watch("--out", self.out, state="FAILED", code="1:0")
        self.assertTrue(o["failed"] and o["done"], o)
        self.assertEqual(o["mode"], "chain")
        self.assertEqual(o["first_failed_job"], "9900002")           # the last failed id it scans
        self.assertEqual(o.get("exit_code"), "1:0")

    def test_out_with_no_job_ids_stays_a_local_poll(self):
        o = self.watch("--out", self.out, state="FAILED")
        self.assertEqual(o["mode"], "local")
        self.assertFalse(o["failed"])


class SageStackOverflow(Base):
    TEXT = "thread 'main' has overflowed its stack\nfatal runtime error: stack overflow\n"

    def test_log_text_names_it(self):
        self.write("sage_9900001.log", self.TEXT)
        o = self.watch("--slurm", "9900001", "--out", self.out, state="FAILED", code="134:0")
        self.assertEqual(o["error_class"], "sage_stack_overflow", o)
        self.assertIn("RUST_MIN_STACK", o["fix"])
        self.assertIn("re-run", o["fix"].lower())

    def test_the_jobs_own_record_names_it(self):
        self.write("sage_failure.json", json.dumps({
            "failure": "sage_stack_overflow", "exit_code": 134,
            "hint": "re-run this job once; RUST_MIN_STACK is raised in jobs from 2.15"}))
        o = self.watch("--slurm", "9900001", "--out", self.out, state="FAILED", code="134:0")
        self.assertEqual(o["error_class"], "sage_stack_overflow", o)
        self.assertIn("re-run this job once", o["fix"])

    def test_a_chain_names_it_too(self):
        self.write("jobs.txt", "9900001\n")
        self.write("sage_9900001.log", self.TEXT)
        o = self.watch("--all", self.out, state="FAILED", code="134:0")
        self.assertEqual(o["error_class"], "sage_stack_overflow", o)

    def test_a_plain_failure_is_not_a_stack_overflow(self):
        self.write("sage_9900001.log", "error: no such file\n")
        o = self.watch("--slurm", "9900001", "--out", self.out, state="FAILED", code="1:0")
        self.assertNotEqual(o["error_class"], "sage_stack_overflow")


class NeedsDecision(Base):
    def marker(self, d, job_id="9900001"):
        # as normalization_check.write records it: the job that wrote it (SLURM_JOB_ID)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "NEEDS_DECISION.json"), "w") as fh:
            json.dump({"needs_decision": True, "exit_code": 3, "why": "the check tripped",
                       "summary": os.path.join(d, "NORMALIZATION_CHECK.md"),
                       "job_id": job_id}, fh)

    def test_a_completed_check_job_says_a_decision_is_owed(self):
        s = os.path.join(self.tmp.name, "S")
        self.marker(os.path.join(s, "output", "norm_check"))
        os.makedirs(os.path.join(s, "logs"))
        log = os.path.join(s, "logs", "de_check_9900001.log")
        with open(log, "w") as fh:
            fh.write("[normalization_check] TRIPPED\n")
        o = self.watch("--slurm", "9900001", "--log", log, "--session", s, state="COMPLETED",
                       code="0:0")
        self.assertTrue(o["done"])
        self.assertFalse(o["failed"])
        self.assertTrue(o["needs_decision"], o)
        self.assertEqual(o["state"], "COMPLETED")
        self.assertIn("decision", o["progress"]["summary"].lower())
        self.assertIn("NORMALIZATION_CHECK.md", o["fix"])

    def test_the_session_is_found_from_the_log_folder(self):
        s = os.path.join(self.tmp.name, "S2")
        self.marker(os.path.join(s, "output", "norm_check"))
        os.makedirs(os.path.join(s, "logs"))
        log = os.path.join(s, "logs", "de_check_9900001.log")
        open(log, "w").close()
        o = self.watch("--slurm", "9900001", "--log", log, state="COMPLETED", code="0:0")
        self.assertTrue(o.get("needs_decision"), o)

    def test_another_jobs_marker_is_a_warning_not_this_jobs_summary(self):
        """C2 review ux M8: a later, unrelated COMPLETED job of the session (a search step, a DE)
        was reported as "a normalisation DECISION is owed ... do not resubmit"."""
        s = os.path.join(self.tmp.name, "S3")
        self.marker(os.path.join(s, "output", "norm_check"), job_id="9900777")
        o = self.watch("--slurm", "9900001", "--out", self.out, "--session", s,
                       state="COMPLETED", code="0:0")
        self.assertTrue(o["done"] and not o["failed"])
        self.assertFalse(o.get("needs_decision", False), o)
        self.assertNotIn("decision", o["progress"]["summary"].lower())
        self.assertTrue(any("job 9900777" in w and "decision is owed" in w
                            for w in o.get("warnings", [])), o.get("warnings"))

    def test_a_marker_with_no_job_id_is_a_warning(self):
        s = os.path.join(self.tmp.name, "S4")
        self.marker(os.path.join(s, "output", "norm_check"), job_id=None)
        o = self.watch("--slurm", "9900001", "--out", self.out, "--session", s,
                       state="COMPLETED", code="0:0")
        self.assertFalse(o.get("needs_decision", False), o)
        self.assertTrue(any("no Slurm job id" in w for w in o.get("warnings", [])), o)

    def test_an_unreadable_marker_is_said(self):
        if os.geteuid() == 0:
            self.skipTest("root reads anything")
        s = os.path.join(self.tmp.name, "S5")
        d = os.path.join(s, "output", "norm_check")
        self.marker(d)
        os.chmod(os.path.join(d, "NEEDS_DECISION.json"), 0)
        self.addCleanup(os.chmod, os.path.join(d, "NEEDS_DECISION.json"), 0o644)
        o = self.watch("--slurm", "9900001", "--out", self.out, "--session", s,
                       state="COMPLETED", code="0:0")
        self.assertTrue(any("could not be read" in w for w in o.get("warnings", [])), o)

    def test_no_marker_no_question(self):
        o = self.watch("--slurm", "9900001", "--out", self.out, state="COMPLETED", code="0:0")
        self.assertFalse(o.get("needs_decision", False))


if __name__ == "__main__":
    unittest.main()
