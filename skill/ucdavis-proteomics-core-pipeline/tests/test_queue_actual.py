#!/usr/bin/env python3
"""
search_provenance.json records where the jobs actually ran, beside the queue that was requested
(2.15 acceptance, lane X4: a pending job moved to publicgrp/low was recorded as genome-center-grp/
high). A fake sacct stands in for Slurm; nothing contacts a cluster. Hermetic.
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

import queue_actual as qa  # noqa: E402

TOOL = os.path.join(SCRIPTS, "queue_actual.py")
REQUESTED = {"partition": "high", "account": "genome-center-grp", "qos": None, "exclude": None}


def fake_sacct(d, rows):
    p = os.path.join(d, "sacct")
    with open(p, "w") as fh:
        fh.write("#!/bin/bash\n" + "".join(f"echo '{r}'\n" for r in rows))
    os.chmod(p, 0o755)
    return p


class Base(unittest.TestCase):
    def setUp(self):
        self.td = tempfile.TemporaryDirectory()
        self.addCleanup(self.td.cleanup)
        self.d = self.td.name
        self.prov = os.path.join(self.d, "search_provenance.json")
        with open(self.prov, "w") as fh:
            json.dump({"engine": "diann", "queue": REQUESTED}, fh)
        with open(os.path.join(self.d, "jobs.txt"), "w") as fh:
            fh.write("9900001\n9900002\n")

    def run_tool(self, sacct_rows, *extra, jobs_file=True):
        env = dict(os.environ, QUEUE_ACTUAL_SACCT=fake_sacct(self.d, sacct_rows))
        env.pop("SLURM_JOB_ID", None)
        args = [sys.executable, TOOL, "--provenance", self.prov]
        if jobs_file:
            args += ["--jobs-file", os.path.join(self.d, "jobs.txt")]
        r = subprocess.run(args + list(extra), capture_output=True, text=True, env=env, timeout=60)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(self.prov) as fh:
            return json.load(fh)


class Records(Base):
    def test_a_moved_job_is_recorded_where_it_ran_and_the_request_is_kept(self):
        prov = self.run_tool(["9900001|publicgrp|publicgrp|publicgrp-low-qos|COMPLETED",
                              "9900002|publicgrp|publicgrp|publicgrp-low-qos|COMPLETED"])
        qa_rec = prov["queue_actual"]
        self.assertEqual(prov["queue"], REQUESTED)                  # the request is untouched
        self.assertEqual(qa_rec["partition"], ["publicgrp"])
        self.assertEqual(qa_rec["account"], ["publicgrp"])
        self.assertEqual(qa_rec["qos"], ["publicgrp-low-qos"])
        self.assertTrue(qa_rec["moved"])
        self.assertEqual(qa_rec["requested"], {"partition": "high", "account": "genome-center-grp"})
        self.assertEqual([j["id"] for j in qa_rec["jobs"]], ["9900001", "9900002"])

    def test_jobs_that_ran_where_requested_are_not_moved(self):
        prov = self.run_tool(["9900001|high|genome-center-grp|genome-center-grp-high-qos|COMPLETED"])
        self.assertFalse(prov["queue_actual"]["moved"])

    def test_an_array_is_one_job(self):
        prov = self.run_tool(["9900001_1|high|genome-center-grp|q|COMPLETED",
                              "9900001_2|high|genome-center-grp|q|COMPLETED"])
        self.assertEqual(len(prov["queue_actual"]["jobs"]), 1)

    def test_no_request_to_compare_is_null_not_false(self):
        with open(self.prov, "w") as fh:
            json.dump({"engine": "sage"}, fh)
        prov = self.run_tool(["9900001|high|genome-center-grp|q|COMPLETED"])
        self.assertIsNone(prov["queue_actual"]["moved"])

    def test_no_sacct_is_said_never_left_out(self):
        env = dict(os.environ, QUEUE_ACTUAL_SACCT="/nonexistent/sacct", PATH="/nonexistent")
        r = subprocess.run([sys.executable, TOOL, "--provenance", self.prov, "--jobs-file",
                            os.path.join(self.d, "jobs.txt")], capture_output=True, text=True,
                           env=env, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(self.prov) as fh:
            rec = json.load(fh)["queue_actual"]
        self.assertIn("unavailable", rec)

    def test_a_job_id_alone_is_enough(self):
        os.remove(os.path.join(self.d, "jobs.txt"))
        prov = self.run_tool(["9900009|low|publicgrp|q|RUNNING"], "--job", "9900009",
                             jobs_file=False)
        self.assertEqual(prov["queue_actual"]["jobs"][0]["id"], "9900009")


class WiredIntoTheJobEnd(unittest.TestCase):
    def test_the_job_end_hook_records_the_queue(self):
        src = open(os.path.join(SCRIPTS, "notify_slack.py")).read()
        body = src[src.index("def search_done("):src.index("def announce_gate(")]
        self.assertIn("queue_actual", body)
        self.assertLess(body.index("queue_actual"), body.index('record_run("search-done"'))


if __name__ == "__main__":
    unittest.main()
