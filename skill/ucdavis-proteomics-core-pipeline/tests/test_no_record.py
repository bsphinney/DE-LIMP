#!/usr/bin/env python3
"""
A test run must be able to stay out of the Core's shared run records (lane X7; pep22 F16, ochro
F5, pxd4, bat4, PROT_9812 F11). Before this the only opt-out was the environment variable
RECORD_RUN=off, which a hive_remote job never saw and which agents had to know about. Now:

  * record_run.no_record_reason() is the ONE definition of "should this run be recorded": the
    --no-record flag, RECORD_RUN=off / SKILL_RUNS_DIR=off, or a /skill_test/ folder in the
    session or search path;
  * record_run.py search-done / analysis-done honour --no-record and a /skill_test/ path before
    anything is read or written (the registry folder stays untouched);
  * the job-end hook (notify_slack.wrap_job_script / job_end_plan / job_end_flags) bakes the
    choice into the job script as --no-record, and the plan records why;
  * notify_slack.record_run() -- what finalize and the job hook call -- never starts record_run.py
    for such a run.

Nothing here can reach the real registry: SKILL_RUNS_DIR is a temp dir, and the scripts are
imported or run with RECORD_RUN-family variables removed.
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
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
import notify_slack  # noqa: E402
import record_run  # noqa: E402

CLEAN = {k: v for k, v in os.environ.items()
         if not k.startswith(("RECORD_RUN", "SKILL_RUNS", "SLURM_", "HIVE_"))}


def mkdirs(*p):
    q = os.path.join(*p)
    os.makedirs(q, exist_ok=True)
    return q


class OneDefinition(unittest.TestCase):
    def setUp(self):
        self._p = mock.patch.dict(os.environ, CLEAN, clear=True)
        self._p.start()
        self.addCleanup(self._p.stop)

    def test_an_ordinary_run_is_recorded(self):
        self.assertIsNone(record_run.no_record_reason("/quobyte/x/sessions/2026-01-01_a"))
        self.assertIsNone(record_run.no_record_reason(None, None))

    def test_the_flag_wins_and_says_so(self):
        self.assertIn("--no-record", record_run.no_record_reason("/a/b", flag=True))

    def test_the_environment_switch_still_works(self):
        for k, v in (("RECORD_RUN", "off"), ("SKILL_RUNS_DIR", "off")):
            with mock.patch.dict(os.environ, {k: v}):
                self.assertIn("RECORD_RUN=off", record_run.no_record_reason("/a/b"))

    def test_a_skill_test_folder_anywhere_in_the_path(self):
        for path in ("/q/skill_test/s1", "/q/skill_test/s1/output/search", "/q/skill_test",
                     "/q/skill_test/"):
            why = record_run.no_record_reason(path)
            self.assertIn("/skill_test/", why, path)
            self.assertIn("never recorded", why)

    def test_a_look_alike_is_still_recorded(self):
        for path in ("/q/skill_tests/s1", "/q/my_skill_test/s1", "/q/skill_test_2/s1"):
            self.assertIsNone(record_run.no_record_reason(path), path)

    def test_any_of_several_paths_counts(self):
        self.assertIsNotNone(record_run.no_record_reason("/ok/a", "/x/skill_test/b"))

    def test_a_link_into_skill_test_counts(self):
        with tempfile.TemporaryDirectory() as d:
            real = mkdirs(d, "skill_test", "s")
            link = os.path.join(d, "alias")
            os.symlink(real, link)
            self.assertIsNotNone(record_run.no_record_reason(link))


class RecordRunCli(unittest.TestCase):
    def run_rr(self, d, *args, **env):
        e = dict(CLEAN, SKILL_RUNS_DIR=os.path.join(d, "runs"), RECORD_RUN_SACCT="none",
                 TMPDIR=d, **env)
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "record_run.py"), *args],
                              capture_output=True, text=True, env=e, timeout=60)

    def check_not_recorded(self, d, r, needle):
        self.assertEqual(r.returncode, 0, r.stderr)
        res = json.loads(r.stdout)
        self.assertEqual((res["recorded"], res["reason"]), (False, "disabled"), res)
        self.assertIn(needle, res["detail"])
        self.assertIn("not recorded", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(d, "runs")), "the registry was touched")

    def test_no_record_flag_search_done(self):
        with tempfile.TemporaryDirectory() as d:
            out = mkdirs(d, "s", "output", "search")
            self.check_not_recorded(d, self.run_rr(d, "search-done", "--out", out, "--no-record"),
                                    "--no-record")

    def test_no_record_flag_analysis_done(self):
        with tempfile.TemporaryDirectory() as d:
            ses = mkdirs(d, "s")
            self.check_not_recorded(
                d, self.run_rr(d, "analysis-done", "--session", ses, "--no-record"), "--no-record")

    def test_skill_test_path_search_done_and_analysis_done(self):
        with tempfile.TemporaryDirectory() as d:
            ses = mkdirs(d, "skill_test", "s")
            out = mkdirs(ses, "output", "search")
            self.check_not_recorded(d, self.run_rr(d, "search-done", "--out", out), "/skill_test/")
            self.check_not_recorded(d, self.run_rr(d, "analysis-done", "--session", ses),
                                    "/skill_test/")

    def test_a_normal_path_is_not_skipped_by_the_new_check(self):
        with tempfile.TemporaryDirectory() as d:
            out = mkdirs(d, "s", "output", "search")
            r = self.run_rr(d, "search-done", "--out", out)
            res = json.loads(r.stdout)
            self.assertNotIn("--no-record", res.get("detail") or "")
            self.assertNotIn("skill_test", res.get("detail") or "")


class JobEndHook(unittest.TestCase):
    def setUp(self):
        self._p = mock.patch.dict(os.environ, CLEAN, clear=True)
        self._p.start()
        self.addCleanup(self._p.stop)

    def wrap(self, out="/q/s/output/search", **kw):
        return notify_slack.wrap_job_script("#!/bin/bash\n#SBATCH --x\nwork\n", out, final=True,
                                            time_limit_h=1, stage="s", scripts_dir="/sc", **kw)

    def test_default_hook_records(self):
        self.assertNotIn("--no-record", self.wrap())

    def test_flag_is_baked_into_the_hook(self):
        self.assertIn(" --no-record", self.wrap(record=False))

    def test_a_skill_test_out_is_baked_in_without_the_flag(self):
        self.assertIn(" --no-record", self.wrap(out="/q/skill_test/s/output/search"))

    def test_the_environment_is_read_where_the_job_runs_not_baked_in(self):
        with mock.patch.dict(os.environ, {"RECORD_RUN": "off"}):
            self.assertNotIn(" --no-record", self.wrap())

    def test_plan_says_why(self):
        plan = notify_slack.job_end_plan(record=False)
        self.assertEqual(plan["run_log"], "off")
        self.assertIn("--no-record", plan["run_log_off_because"])
        plan = notify_slack.job_end_plan(out="/q/skill_test/s")
        self.assertEqual(plan["run_log"], "off")
        self.assertIn("/skill_test/", plan["run_log_off_because"])
        plan = notify_slack.job_end_plan(out="/q/s")
        self.assertEqual(plan["run_log"], "on")
        self.assertNotIn("run_log_off_because", plan)

    def test_flags_and_a_recorded_plan_round_trip(self):
        self.assertIn("--no-record", notify_slack.job_end_flags(record=False))
        self.assertNotIn("--no-record", notify_slack.job_end_flags())
        plan = notify_slack.job_end_plan(record=False)
        self.assertIn("--no-record", notify_slack.job_end_argv(plan))
        self.assertNotIn("--no-record", notify_slack.job_end_argv(notify_slack.job_end_plan()))

    def test_search_done_never_starts_record_run_for_such_a_run(self):
        calls = []
        with mock.patch.object(notify_slack, "_run_helper",
                               side_effect=lambda *a, **k: calls.append(a) or (0, "{}", "")):
            r = notify_slack.record_run("search-done", out="/q/s", session=None, status="ok",
                                        exit_code=0, no_record=True)
            r2 = notify_slack.record_run("analysis-done", session="/q/skill_test/s")
        self.assertEqual(calls, [])
        for x in (r, r2):
            self.assertTrue(x["off"] and x["logged"] is False, x)
        self.assertIn("--no-record", r["detail"])
        self.assertIn("/skill_test/", r2["detail"])
        self.assertEqual(notify_slack.run_log_manifest(r)[0], "INFO")



class TestRunsStayOutOfFranAndSlack(unittest.TestCase):
    """C2 review ux H2: the ONE no-record reason (--no-record, a /skill_test/ folder) keeps a test
    run out of FRAN staging and the Core Slack channel as well as the run log -- at generation
    (job_end_plan, baked into the job script), where the job ends (search_done, read again) and
    in fran_deposit.py check itself. RECORD_RUN=off is the registry's switch only."""

    def setUp(self):
        self._p = mock.patch.dict(os.environ, CLEAN, clear=True)
        self._p.start()
        self.addCleanup(self._p.stop)

    def test_the_plan_turns_fran_and_slack_off_and_says_why(self):
        for kw in ({"record": False}, {"out": "/q/skill_test/s/output/search"}):
            plan = notify_slack.job_end_plan(fran_guarded=True, **kw)
            self.assertEqual((plan["run_log"], plan["slack"], plan["fran"]), ("off",) * 3, kw)
            self.assertEqual(plan["fran_off_because"], plan["run_log_off_because"])
            self.assertEqual(plan["slack_off_because"], plan["run_log_off_because"])
        plan = notify_slack.job_end_plan(fran_guarded=True, out="/q/s/output/search")
        self.assertEqual((plan["run_log"], plan["slack"], plan["fran"]), ("on", "on", "stage"))
        self.assertNotIn("fran_off_because", plan)

    def test_the_environment_switch_is_the_registry_only(self):
        with mock.patch.dict(os.environ, {"RECORD_RUN": "off"}):
            plan = notify_slack.job_end_plan(fran_guarded=True, out="/q/s/output/search")
        self.assertEqual((plan["slack"], plan["fran"]), ("on", "stage"))

    def test_the_job_script_carries_all_three(self):
        js = notify_slack.wrap_job_script("#!/bin/bash\n#SBATCH --x\nwork\n",
                                          "/q/skill_test/s/output/search", final=True,
                                          time_limit_h=1, stage="s", scripts_dir="/sc",
                                          fran_guarded=True)
        for flag in (" --no-record", " --no-slack", " --no-fran"):
            self.assertIn(flag, js)
        self.assertNotIn("--fran-stage", js)

    def _search_done(self, out, **kw):
        seen = {"fran": [], "posts": []}
        facts = {"status": "ok", "location": out, "kind": "search"}
        with mock.patch.object(notify_slack, "search_facts", return_value=dict(facts)), \
                mock.patch.object(notify_slack, "_session_of", return_value=None), \
                mock.patch.object(notify_slack, "_fran_step",
                                  side_effect=lambda o, s, mode, n, q: seen["fran"].append(mode)
                                  or {"staged": False, "reason": "opted_out"}), \
                mock.patch.object(notify_slack, "record_run", return_value={"logged": False,
                                                                           "off": True}), \
                mock.patch.object(notify_slack, "deliver",
                                  side_effect=lambda f, **k: seen["posts"].append(f)
                                  or (True, "sent")):
            sent, msg = notify_slack.search_done(out, exit_code=0, final=True, fran="stage",
                                                 slack=True, **kw)
        return sent, msg, seen

    def test_the_job_end_hook_reads_it_again(self):
        """An older job script (generated without the plan) still keeps a test run out."""
        for out, kw in (("/q/skill_test/s/output/search", {}),
                        ("/q/s/output/search", {"record": False})):
            sent, msg, seen = self._search_done(out, **kw)
            self.assertEqual(seen["fran"], ["off"], out)
            self.assertEqual(seen["posts"], [], out)
            self.assertFalse(sent)
            self.assertIn("a test run", msg)
            self.assertIn("FRAN: skipped: a test run", msg)
            self.assertNotIn("/q/", msg)                      # words, not the folder

    def test_an_ordinary_search_still_stages_and_posts(self):
        sent, msg, seen = self._search_done("/q/s/output/search")
        self.assertEqual(seen["fran"], ["stage"])
        self.assertEqual(len(seen["posts"]), 1)
        self.assertTrue(sent)

    def test_fran_deposit_refuses_a_skill_test_search_and_records_it(self):
        import argparse
        import fran_deposit
        self.assertIn("test_run", fran_deposit.DECISION_STATUS)
        with tempfile.TemporaryDirectory() as d:
            out = mkdirs(d, "skill_test", "s", "output", "search")
            c = fran_deposit.check(argparse.Namespace(out=out, skip=False, name=None, qc=False,
                                                      not_qc=False))
            self.assertEqual(c["reason"], "test_run", c)
            self.assertFalse(c["eligible"])
            ordinary = mkdirs(d, "s", "output", "search")
            c = fran_deposit.check(argparse.Namespace(out=ordinary, skip=False, name=None,
                                                      qc=False, not_qc=False))
            self.assertNotEqual(c["reason"], "test_run", c)


if __name__ == "__main__":
    unittest.main()
