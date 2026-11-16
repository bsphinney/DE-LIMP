#!/usr/bin/env python3
"""
watch_run.sh reads every Slurm state from ONE list (slurm_states.py), in chain mode and single-job
mode alike (2.15, C2 review ux H1 and M4).

  * PREEMPTED (not requeued), DEADLINE, BOOT_FAIL, REVOKED and SPECIAL_EXIT are failures: done and
    failed, in a chain (which read them as "done, not failed") and for one job (which read them as
    running forever);
  * a REQUEUED job (preempted, then put back by Slurm) is not done and not failed;
  * a chain whose job ids neither sacct nor squeue answers for is watcher_query_failed, with the
    ids -- never a quiet "running";
  * record_run.py and checkpoint.py read the same list;
  * a CANCELLED job is "cancelled", never a class read from its log; "empty_results" needs a failed
    job and the engines' own words; a chain of Skyline jobs (labelled ids, logs/ folder) reads its
    ids and its failed job's log (pm-215x check d).

A fake sacct/squeue answers per job id; nothing contacts a cluster.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
WATCH = os.path.join(SCRIPTS, "watch_run.sh")
sys.path.insert(0, SCRIPTS)

import slurm_states  # noqa: E402

# FAKE_STATES="<id>=<STATE>,<id>=<STATE>": sacct answers that state for that id, nothing for others.
# FAKE_SQ the same for squeue. -P rows for the exit-code / node-fault queries.
FAKE_SACCT = r'''#!/bin/bash
id=""; prev=""
for a in "$@"; do [ "$prev" = "-j" ] && id="$a"; prev="$a"; done
st="$(printf '%s' "$FAKE_STATES" | tr ',' '\n' | sed -n "s/^$id=//p" | head -n1)"
[ -z "$st" ] && exit 0
case "$*" in
  *-P*) echo "$id|$st|0:0|hive-node-1|00:10:00" ;;
  *)    echo "  $st" ;;
esac
'''
FAKE_SQUEUE = r'''#!/bin/bash
id=""; prev=""
for a in "$@"; do [ "$prev" = "-j" ] && id="$a"; prev="$a"; done
case "$*" in *%R*) exit 0 ;; esac
printf '%s' "${FAKE_SQ:-}" | tr ',' '\n' | sed -n "s/^$id=//p" | head -n1
'''


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

    def watch(self, *args, states="", sq=""):
        env = dict(os.environ, PATH=self.bin + os.pathsep + os.environ.get("PATH", ""),
                   FAKE_STATES=states, FAKE_SQ=sq)
        r = subprocess.run(["bash", WATCH, *args], capture_output=True, text=True, env=env,
                           timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return json.loads(r.stdout)

    def chain(self, ids):
        with open(os.path.join(self.out, "jobs.txt"), "w") as fh:
            fh.write("\n".join(ids) + "\n")


ENDED = ("PREEMPTED", "DEADLINE", "BOOT_FAIL", "REVOKED", "SPECIAL_EXIT")


class ChainFailures(Base):
    def test_every_ended_state_fails_the_chain(self):
        self.chain(["9900001", "9900002"])
        for st in ENDED + ("OUT_OF_ME+",):
            with self.subTest(state=st):
                o = self.watch("--all", self.out, states=f"9900001=COMPLETED,9900002={st}")
                self.assertTrue(o["failed"] and o["done"], o)
                self.assertEqual(o["first_failed_job"], "9900002")

    def test_a_preempted_step_says_resubmit_on_high(self):
        self.chain(["9900001"])
        o = self.watch("--all", self.out, states="9900001=PREEMPTED")
        self.assertEqual(o["error_class"], "preempted", o)
        self.assertIn("high queue", o["fix"])

    def test_requeued_is_not_done_and_not_failed(self):
        self.chain(["9900001", "9900002"])
        o = self.watch("--all", self.out, states="9900001=COMPLETED,9900002=REQUEUED")
        self.assertFalse(o["failed"] or o["done"], o)
        self.assertEqual(o["progress"]["requeued_jobs"], ["9900002"])
        self.assertIn("REQUEUED", o["progress"]["summary"])

    def test_completed_chain_is_still_done(self):
        self.chain(["9900001", "9900002"])
        o = self.watch("--all", self.out, states="9900001=COMPLETED,9900002=COMPLETED")
        self.assertTrue(o["done"] and not o["failed"], o)

    def test_a_state_the_list_does_not_know_is_not_done(self):
        self.chain(["9900001"])
        o = self.watch("--all", self.out, states="9900001=SOMETHING_NEW")
        self.assertFalse(o["done"], o)


class SingleFailures(Base):
    def test_every_ended_state_fails_the_job(self):
        for st in ENDED:
            with self.subTest(state=st):
                o = self.watch("--slurm", "9900001", "--out", self.out, states=f"9900001={st}")
                self.assertTrue(o["failed"] and o["done"], o)
                self.assertNotIn("running", o["progress"]["summary"].lower())

    def test_preempted_class(self):
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=PREEMPTED")
        self.assertEqual(o["error_class"], "preempted", o)
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=DEADLINE")
        self.assertEqual(o["error_class"], "slurm_ended", o)
        # BOOT_FAIL is a node that failed to start: node_fault.py's class, with its retry
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=BOOT_FAIL")
        self.assertEqual(o["error_class"], "node_fault", o)

    def test_a_failed_job_never_says_searching(self):
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=PREEMPTED")
        self.assertTrue(o["doing"].startswith("Stopped"), o["doing"])

    def test_requeued_job_waits(self):
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=REQUEUED")
        self.assertFalse(o["failed"] or o["done"], o)
        self.assertEqual(o["state"], "REQUEUED")
        self.assertIn("run again", o["progress"]["summary"])


class ChainUnanswered(Base):
    def test_no_answer_for_any_job_is_query_failed(self):
        self.chain(["9900001", "9900002"])
        o = self.watch("--all", self.out, states="")
        self.assertEqual(o["error_class"], "watcher_query_failed", o)
        self.assertEqual(o["state"], "query_failed")
        self.assertEqual(o["unknown_jobs"], ["9900001", "9900002"])
        self.assertFalse(o["failed"] or o["done"])
        self.assertIn("NOT a job failure", o["fix"])

    def test_one_unanswered_job_is_named(self):
        self.chain(["9900001", "9900002"])
        o = self.watch("--all", self.out, states="9900001=COMPLETED")
        self.assertEqual(o["error_class"], "watcher_query_failed", o)
        self.assertEqual(o["unknown_jobs"], ["9900002"])
        self.assertFalse(o["done"])

    def test_squeue_answers_when_sacct_has_no_row_yet(self):
        self.chain(["9900001"])
        o = self.watch("--all", self.out, states="", sq="9900001=PENDING")
        self.assertNotIn("error_class", o)
        self.assertFalse(o["done"])


# A FragPipe log as the real cancelled acceptance job left it (anonymised): its count lines hold
# "...=1340 proteins=..." and the Slurm CANCELLED line (pm-215x check d).
FRAGPIPE_CANCELLED_LOG = """\
time="16:15:59" level=info msg="Final report numbers after FDR filtering, and post-processing" ions=1478 peptides=1340 proteins=443 psms=2155
time="16:16:12" level=info msg="Final report numbers after FDR filtering, and post-processing" ions=2341 peptides=2080 proteins=597 psms=3484
an empty line in the workflow is skipped
[2026-11-11T16:16:15.498] error: *** JOB 9900001 ON hive-node-1 CANCELLED AT 2026-11-11T16:16:15 DUE to SIGNAL Terminated ***
"""


class ErrorClassFromTheLog(Base):
    def log(self, name, text, sub=""):
        d = os.path.join(self.out, sub) if sub else self.out
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, name), "w") as fh:
            fh.write(text)

    def test_a_cancelled_fragpipe_job_is_cancelled_not_empty_results(self):
        self.log("fragpipe_search_9900001.log", FRAGPIPE_CANCELLED_LOG)
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=CANCELLED")
        self.assertTrue(o["failed"] and o["done"], o)
        self.assertIn("proteins=443", o["log_tail"])        # the log was read
        self.assertEqual(o["error_class"], "cancelled", o)
        self.assertIn("scancel", o["fix"])
        self.assertNotIn("FASTA", o["fix"])

    def test_a_failed_job_with_counts_and_the_word_empty_is_not_empty_results(self):
        self.log("fragpipe_search_9900001.log", FRAGPIPE_CANCELLED_LOG.replace("CANCELLED", "ended"))
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=FAILED")
        self.assertNotEqual(o["error_class"], "empty_results", o)

    def test_a_failed_job_that_found_0_proteins_is_empty_results(self):
        self.log("search_9900001.log", "Protein inference: 0 proteins after FDR\n")
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=FAILED")
        self.assertEqual(o["error_class"], "empty_results", o)

    def test_a_completed_job_is_never_empty_results(self):
        self.log("search_9900001.log", "No precursors for this charge state, skipped\n")
        o = self.watch("--slurm", "9900001", "--out", self.out, states="9900001=COMPLETED")
        self.assertTrue(o["done"] and not o["failed"], o)
        self.assertNotEqual(o["error_class"], "empty_results", o)


class SkylineChain(Base):
    """skyline_peaks.py's submit.sh: jobs.txt lines are "<label><TAB><id>", logs in logs/."""

    def setUp(self):
        super().setUp()
        with open(os.path.join(self.out, "jobs.txt"), "w") as fh:
            fh.write("build\t9900001\nimport\t9900002\njoin\t9900003\n")
        os.makedirs(os.path.join(self.out, "logs"))
        with open(os.path.join(self.out, "logs", "sky_build_9900001.log"), "w") as fh:
            fh.write("SkylineCmd: the isolation scheme file could not be read\n")

    def test_labels_are_not_read_as_job_ids(self):
        o = self.watch("--all", self.out,
                       states="9900001=FAILED,9900002=CANCELLED,9900003=COMPLETED")
        self.assertNotIn("UNKNOWN", o["chain"], o)
        self.assertNotIn("build", o["chain"], o)
        self.assertEqual(o["chain"], "9900001=FAILED 9900002=CANCELLED 9900003=COMPLETED")
        self.assertEqual(o["first_failed_job"], "9900001")

    def test_the_failed_skyline_job_log_is_read_from_logs(self):
        o = self.watch("--all", self.out,
                       states="9900001=FAILED,9900002=CANCELLED,9900003=COMPLETED")
        self.assertIn("isolation scheme file could not be read", o.get("log_tail", ""), o)
        self.assertNotIn("log_unread", o)

    def test_a_label_only_chain_is_not_watched_as_running(self):
        o = self.watch("--all", self.out, states="9900001=COMPLETED,9900002=COMPLETED,"
                                                 "9900003=COMPLETED")
        self.assertTrue(o["done"] and not o["failed"], o)
        self.assertNotIn("unknown_jobs", o)


class OneList(unittest.TestCase):
    def test_python_readers_share_it(self):
        import checkpoint
        import record_run
        # record_run.py keeps a copy (it is piped to a remote python, no sibling import there)
        self.assertEqual(record_run.BAD_STATES, slurm_states.BAD_STATES)
        self.assertEqual(record_run.LIVE_STATES, slurm_states.LIVE_STATES)
        self.assertIs(checkpoint.TERMINAL_BAD, slurm_states.BAD_STATES)

    def test_groups(self):
        for st in ENDED + ("FAILED", "CANCELLED", "OUT_OF_ME+", "CANCELLED by 1"):
            self.assertEqual(slurm_states.classify(st), "bad", st)
        for st in ("REQUEUED", "REQUEUE_H", "SUSPENDED", "PENDING", "RUNNING"):
            self.assertEqual(slurm_states.classify(st), "live", st)
        self.assertEqual(slurm_states.classify("COMPLETED"), "ok")
        self.assertEqual(slurm_states.classify(""), "unknown")


if __name__ == "__main__":
    unittest.main()
