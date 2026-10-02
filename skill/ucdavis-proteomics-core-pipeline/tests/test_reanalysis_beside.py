#!/usr/bin/env python3
"""
A re-analysis BESIDE the session it re-does, and a finalize that can be run again.

A Core re-analysis (2026-09-25): the prior session was a delivered, read-only-by-policy folder;
`session.py init --reanalysis-of` could only nest the new run INSIDE it, so the agent patched
session._raw_dir from a wrapper and wrote .reanalysis_of by hand. `--beside` puts it next to the
prior as <date>_<name>_v<N>, writes nothing into the prior, and links the two (reanalysis.json,
the README / AGENTS.md line, DIFFERENCES.md, the run log).

2026-09-28: two finalize runs on one session (before and after adding the podcast) left two
identical "analysis complete" blocks in the Core run log. A re-finalize of an unchanged analysis
now logs nothing new and posts nothing (record_run analysis_digest; tests/test_record_run.py
pins the log side).
"""
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)
import session  # noqa: E402
import session_docs  # noqa: E402
import test_deposit_package as tdp  # noqa: E402
from job_env import job_env  # noqa: E402


def _read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def snapshot(root):
    """Every path under `root` with its size and mtime: proof nothing was written there."""
    out = {}
    for d, dirs, files in os.walk(root):
        for n in dirs + files:
            p = os.path.join(d, n)
            st = os.lstat(p)
            out[os.path.relpath(p, root)] = (st.st_size, st.st_mtime_ns)
    return out


class Beside(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.p = tdp.dia_session(self.d)                    # <d>/sessions/2026-09-24_demo
        self.prior = self.p["session_dir"]

    def tearDown(self):
        self._tmp.cleanup()

    def init(self, *extra, name="demo", check=True):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init",
                            "--name", name, "--date", "2026-10-01", *extra],
                           capture_output=True, text=True, env=job_env(self.d), timeout=120)
        if check:
            self.assertEqual(r.returncode, 0, r.stderr)
            return json.loads(r.stdout)
        return r

    def test_v2_beside_the_prior_and_nothing_written_into_it(self):
        before = snapshot(self.prior)
        res = self.init("--reanalysis-of", self.prior, "--beside")
        new = res["created"]
        self.assertEqual(os.path.dirname(new), os.path.dirname(self.prior))
        self.assertEqual(os.path.basename(new), "2026-10-01_demo_v2")
        self.assertEqual((res["placement"], res["version"]), ("beside-prior", 2))
        self.assertEqual(snapshot(self.prior), before)
        self.assertFalse(os.path.exists(os.path.join(self.prior, "reanalysis")))
        link = json.loads(_read(os.path.join(new, session.REANALYSIS_JSON)))
        self.assertEqual((link["reanalysis_of"], link["version"], link["prior_version"]),
                         (self.prior, 2, 1))
        self.assertEqual(link["prior_records"]["search_prov"], self.p["search_prov"])
        self.assertEqual(_read(os.path.join(new, ".reanalysis_of")).strip(), self.prior)
        self.assertIn("this is version 2", _read(os.path.join(new, "README.md")))

    def test_a_second_re_analysis_never_reuses_a_version(self):
        a = self.init("--reanalysis-of", self.prior, "--beside")["created"]
        b = self.init("--reanalysis-of", self.prior, "--beside")["created"]
        self.assertEqual([os.path.basename(x) for x in (a, b)],
                         ["2026-10-01_demo_v2", "2026-10-01_demo_v3"])
        c = self.init("--reanalysis-of", b, "--beside", name="demo v9")   # a name saying v9
        self.assertEqual(os.path.basename(c["created"]), "2026-10-01_demo_v4")
        self.assertEqual(c["version"], 4)

    def test_beside_needs_the_prior(self):
        r = self.init("--beside", check=False)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("--reanalysis-of", r.stderr)

    def test_nested_re_analysis_still_works_and_is_versioned(self):
        res = self.init("--reanalysis-of", self.prior)
        self.assertEqual(res["placement"], "reanalysis")
        self.assertTrue(res["created"].startswith(os.path.join(self.prior, "reanalysis")))
        self.assertEqual(res["version"], 2)

    def test_the_link_is_read_back_and_legacy_markers_still_count(self):
        new = self.init("--reanalysis-of", self.prior, "--beside")["created"]
        self.assertEqual(session.reanalysis_of(new)["session"], self.prior)
        os.remove(os.path.join(new, session.REANALYSIS_JSON))          # a 2.9 session
        rec = session.reanalysis_of(new)
        self.assertEqual((rec["session"], rec["version"]), (self.prior, None))
        self.assertIsNone(session.reanalysis_of(self.prior))

    def test_finalize_writes_differences_in_the_new_session_only(self):
        new = self.init("--reanalysis-of", self.prior, "--beside")["created"]
        before = snapshot(self.prior)
        args = type("A", (), dict(dir=new, zip=False, no_deposit=True, no_notify=True,
                                  reanalysis_of=""))()
        out = io.StringIO()
        with mock.patch.dict(os.environ, job_env(self.d), clear=True), \
                redirect_stdout(out), redirect_stderr(io.StringIO()):
            session.do_finalize(args)
        res = json.loads(out.getvalue())
        self.assertEqual(res["reanalysis_of"], self.prior)
        self.assertEqual(res["differences"], os.path.join(new, "DIFFERENCES.md"))
        self.assertIn(f"Re-analysis of `{self.prior}`", _read(res["differences"]))
        self.assertEqual(snapshot(self.prior), before)
        for doc in ("README.md", "AGENTS.md"):
            self.assertIn(f"Re-analysis of `{self.prior}` (this is version 2)",
                          _read(os.path.join(new, doc)), doc)

    def test_the_docs_line_comes_from_the_one_reader(self):
        new = self.init("--reanalysis-of", self.prior, "--beside")["created"]
        f = session_docs.gather(new)
        self.assertEqual(f["reanalysis_of"]["version"], 2)
        self.assertNotIn("Re-analysis", "\n".join(session_docs.summary_lines(
            session_docs.gather(self.prior))))


class RefinalizeDoesNotRepost(unittest.TestCase):
    def hooks(self, run_log):
        import notify_slack
        args = type("A", (), dict(no_notify=False))()
        with mock.patch.object(notify_slack, "record_run", return_value=run_log), \
                mock.patch.object(notify_slack, "analysis_done",
                                  return_value=(True, "sent")) as post, \
                redirect_stderr(io.StringIO()):
            res = session._finish_hooks(args, "/x/session", None)
        return res, post

    def test_an_unchanged_analysis_is_not_posted_again(self):
        res, post = self.hooks({"logged": True, "path": "/r", "analysis_logged": "unchanged"})
        post.assert_not_called()
        self.assertFalse(res["slack"]["sent"])
        self.assertEqual(res["slack"]["level"], "INFO")
        self.assertIn("not posted again", res["slack"]["detail"])

    def test_a_new_or_changed_analysis_is_posted(self):
        for state in ("new", "changed", None):
            res, post = self.hooks({"logged": True, "path": "/r", "analysis_logged": state})
            post.assert_called_once()
            self.assertTrue(res["slack"]["sent"], state)


if __name__ == "__main__":
    unittest.main()
