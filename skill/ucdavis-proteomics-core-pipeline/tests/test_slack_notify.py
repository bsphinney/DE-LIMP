#!/usr/bin/env python3
"""
notify_slack.py is the job-end hook of every search job (and the Slack half of finalize): when
a search job ends it logs the run (record_run.py), hands a finished Core search to FRAN
(fran_deposit.py stage, success only) and posts to the Core's Slack channel, in that order.
Two rules it must never break, ported from STAN's stan/notify.py:

  1. A notifier never takes down its caller: a missing webhook, an HTTP 500, a hang or a DNS
     failure returns False and raises nothing, and the job trap never changes a job's exit
     status or what the report guard decides.
  2. The webhook URL is a bearer credential: it never appears in stdout, stderr, a returned
     message, an exception or a generated job script.

NOTHING here talks to Slack. Every webhook is a loopback mock (http.server in this process);
the one real-looking hooks.slack.com URL is only ever handed to a patched getaddrinfo that
fails, so no packet leaves the machine. The Core's group webhook file is redirected to a path
that does not exist, and relays go through a fake hive_exec.sh. Nothing reaches the real run
log or FRAN's drop directory either: job and finalize tests run a COPY of the scripts with a
fake record_run.py, FRAN_DEPOSIT=off, and FRAN_DROP_DIR inside the test's temp dir.
"""
import base64
import glob
import http.server
import io
import json
import os
import shutil
import signal
import socket
import stat
import subprocess
import sys
import tempfile
import threading
import time
import unittest
import zipfile
from contextlib import redirect_stderr, redirect_stdout
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import notify_slack as ns  # noqa: E402
import run_search  # noqa: E402
from job_env import job_env  # noqa: E402  (the shared isolation; base_env builds on it)

PY = sys.executable
# The CLI tests run a COPY of scripts/ (setUpModule) whose record_run.py is a fake, so no test
# can write to the Core's real run log once record_run.py ships.
NOTIFIER = None
_MODULE_TMP = None
# Shaped like a real webhook, with a token that is obviously fake. Used only where the network
# is patched out (DNS test) or as a value that must never be printed.
FAKE_SLACK = "https://hooks.slack.com/services/T0FAKE000/B0FAKE000/NotARealTokenXYZ123"
SECRET_BIT = "NotARealTokenXYZ123"


class Mock:
    """A loopback webhook. behaviour: ok | 500 | slow."""

    def __init__(self, behaviour="ok", delay=0):
        self.bodies, self.times = [], []
        mock_self = self

        class H(http.server.BaseHTTPRequestHandler):
            def do_POST(self):
                n = int(self.headers.get("Content-Length", 0))
                mock_self.bodies.append(json.loads(self.rfile.read(n)))
                mock_self.times.append(time.time())
                if behaviour == "slow":
                    time.sleep(delay)
                code = 500 if behaviour == "500" else 200
                self.send_response(code)
                self.end_headers()
                # Slack answers with a word; a hostile or buggy server might echo the URL back.
                self.wfile.write(b"ok" if code == 200 else
                                 ("no_service " + self.path).encode())

            def log_message(self, *a):
                pass

        self.srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
        self.srv.daemon_threads = True
        self.url = f"http://127.0.0.1:{self.srv.server_port}/services/T1/B1/{SECRET_BIT}"
        threading.Thread(target=self.srv.serve_forever, daemon=True).start()

    def close(self):
        self.srv.shutdown()
        self.srv.server_close()


def base_env(d, **extra):
    """job_env() -- the shared isolation: no SLURM job, FRAN off with its drop directory in the
    temp dir, no HIVE login -- with Slack and the run log switched back ON, because these tests
    exercise them: the webhook is a loopback mock and record_run.py a fake (scripts_copy), and
    the group webhook file does not exist."""
    env = job_env(d)
    env.pop("SKILL_SLACK")
    env.pop("RECORD_RUN")
    env.update(HOME=d, SKILL_SLACK_TEST_LOOPBACK="1",
               SKILL_SLACK_GROUP_FILE=os.path.join(d, "no_group_webhook"),
               SKILL_CORE_GROUP_DIR=os.path.join(d, "no_core_group"),
               SKILL_ISSUES_DIR=os.path.join(d, "issues"),
               SKILL_ISSUES_LOCAL_DIR=os.path.join(d, "local_issues"),
               HIVE_ENV_FILE=os.path.join(d, "no_hive.env"))
    env.update(extra)
    return env


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def read(path):
    with open(path) as fh:
        return fh.read()


def search_session(root, name="2026-09-24_HeLa_liver_DIA"):
    """A finished single-shot DIA-NN search inside a session layout."""
    sess = os.path.join(root, name)
    out = os.path.join(sess, "output", "search")
    wf = os.path.join(sess, "input", "wf", "workflow.manifest.json")
    write(wf, json.dumps({"acquisition": "DIA", "instruments": ["Orbitrap Astral"],
                          "engine": {"name": "diann", "version": "2.7.0"}}))
    write(os.path.join(out, "search_provenance.json"), json.dumps({
        "engine": "diann", "version": "2.7.0", "n_files": 4, "bundle": wf,
        "search_mode": "single_shot",
        "files": ["/raw/PatientA_tumour_01.raw", "/raw/PatientB_normal_02.raw"],
        "result": {"report": os.path.join(out, "report.parquet"), "mode": "two_job_libfree"}}))
    write(os.path.join(out, "report.stats.tsv"),
          "File.Name\tPrecursors.Identified\tProteins.Identified\n"
          "/raw/PatientA_tumour_01.raw\t45210\t6120\n/raw/PatientB_normal_02.raw\t44001\t6050\n"
          "/raw/PatientC_tumour_03.raw\t46900\t6200\n/raw/blank.raw\t0\t0\n")
    write(os.path.join(out, "report.parquet"), "PAR1" + "x" * 200)
    write(os.path.join(sess, "input", "search.fasta.meta.json"),
          json.dumps({"organism": "Homo sapiens", "taxid": 9606}))
    return sess, out


FAKE_RECORD_RUN = r"""import json, os, sys, time
with open(os.environ.get("FAKE_ORDER_LOG", os.devnull), "a") as fh:
    fh.write(json.dumps({"who": "record_run", "argv": sys.argv[1:], "t": time.time()}) + "\n")
if os.environ.get("FAKE_RECORD_MODE") == "crash":
    sys.exit("record_run exploded")
if os.environ.get("FAKE_RECORD_REASON"):
    print(json.dumps({"recorded": False, "reason": os.environ["FAKE_RECORD_REASON"],
                      "detail": "NotADirectoryError: sessions is a file"}))
    sys.exit(0)
print(json.dumps({"recorded": True, "path": "/runs/skill_runs.tsv"}))
"""

FAKE_FRAN = r"""import json, os, sys, time
# Guarded like the real fran_deposit.py: the notifier IMPORTS it to find stage_argv().
def main():
    with open(os.environ.get("FAKE_ORDER_LOG", os.devnull), "a") as fh:
        fh.write(json.dumps({"who": "fran_deposit", "argv": sys.argv[1:], "t": time.time()}) + "\n")
    if os.environ.get("FAKE_FRAN_MODE") == "crash":
        raise RuntimeError("drop dir exploded")
    out = {"staged": True, "reason": "ok", "entry": "/drop/x__1234",
           "fran_health": {"verdict": "healthy", "summary": "ok"}}
    if "--skip" in sys.argv:
        out = {"staged": False, "reason": "opted_out", "detail": "--skip"}
    elif "--qc" in sys.argv:
        out = {"staged": False, "reason": "qc_run", "detail": "--qc"}
    if os.environ.get("FAKE_FRAN_MODE") == "stuck":
        out["fran_health"] = {"verdict": "stuck", "summary": "nothing ingested for 7 days"}
        out["health_warning"] = "FRAN ingest stuck: nothing ingested since 2026-09-17"
        print("[fran_deposit] WARNING: " + out["health_warning"], file=sys.stderr)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
"""


def scripts_copy(d, record_run=FAKE_RECORD_RUN, fran=None):
    """A copy of the skill's scripts/ (so notify_slack.py finds its helpers beside itself), with
    record_run.py and fran_deposit.py set to a fake's text, removed (False), or kept (None)."""
    sd = os.path.join(d, "skill_scripts")
    if os.path.isdir(sd):
        shutil.rmtree(sd)
    os.makedirs(sd)
    for f in os.listdir(SCRIPTS):
        if os.path.isfile(os.path.join(SCRIPTS, f)):
            shutil.copy2(os.path.join(SCRIPTS, f), sd)
    for name, text in (("record_run.py", record_run), ("fran_deposit.py", fran)):
        path = os.path.join(sd, name)
        if text is False:
            if os.path.exists(path):
                os.remove(path)
        elif text is not None:
            write(path, text)
    return sd


def setUpModule():
    global NOTIFIER, _MODULE_TMP
    _MODULE_TMP = tempfile.mkdtemp()
    NOTIFIER = os.path.join(scripts_copy(_MODULE_TMP), "notify_slack.py")


def tearDownModule():
    shutil.rmtree(_MODULE_TMP, ignore_errors=True)


def order_log(path):
    if not os.path.exists(path):
        return []
    return [json.loads(ln) for ln in read(path).splitlines() if ln.strip()]


def env_patch(d, **extra):
    return mock.patch.dict(os.environ, base_env(d, **extra), clear=True)


class Payloads(unittest.TestCase):
    def test_search_done_payload_names_who_what_where_and_numbers(self):
        with tempfile.TemporaryDirectory() as d:
            sess, out = search_session(d)
            with env_patch(d, SLURM_JOB_ID="4242", SLURM_JOB_USER="alice",
                           SLURM_JOB_NAME="diann_search"):
                f = ns.search_facts(out, exit_code=0, stage="search (job 2 of 2)")
                p = ns.render(f)
            blob = json.dumps(p)
            self.assertIn("Search finished", p["text"])
            # the median is over EVERY run searched -- the blank included, and flagged
            for want in ("2026-09-24_HeLa_liver_DIA", "alice", "DIA-NN 2.7.0", "Orbitrap Astral",
                         "DIA", "44,605 precursors", "6,085 proteins", out, "job 4242",
                         "1 run(s) identified nothing", "library + search jobs"):
                self.assertIn(want, blob, want)
            self.assertIn("*Runs*\\n4", blob)
            # no raw sample annotations: not one input file name reaches Slack
            for leak in ("PatientA", "PatientB", "tumour", "blank.raw"):
                self.assertNotIn(leak, blob)

    def test_failed_intermediate_step_says_the_chain_waits(self):
        with tempfile.TemporaryDirectory() as d:
            _, out = search_session(d)
            with env_patch(d, SLURM_JOB_ID="7"):
                p = ns.render(ns.search_facts(out, exit_code=1, final=False,
                                              stage="step 3/5 empirical library assembly"))
            self.assertIn("FAILED at step 3/5", p["text"])
            self.assertIn("exit 1", p["text"])
            self.assertIn("afterok", json.dumps(p))

    def test_analysis_done_payload_counts_per_contrast_and_points_at_the_package(self):
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            prov = dict(tdp.DE_PROV, n_samples=3,
                        significant_per_contrast={"Treated-Control": 312},
                        significance_rule="adj.P.Val < adjp (BH); no fold-change filter")
            write(os.path.join(p["de_dir"], "de_provenance.json"), json.dumps(prov))
            write(os.path.join(p["deposit_dir"], "HOW_TO_SUBMIT.md"), "# how\n")
            write(p["manifest_txt"], "M\n=\n[OK]      a\n[SKIPPED] b -- why\n")
            zp = p["session_dir"] + ".zip"
            write(zp, "zip")
            with env_patch(d, HIVE_USER="bob"):
                body = json.dumps(ns.render(ns.analysis_facts(p["session_dir"], zp)))
            for want in ("Analysis complete", "Treated-Control", "*312*", "no fold-change filter",
                         "timsTOF HT", "DIA-NN 2.7.0", "DPC-Quant + limma (limpa)",
                         "HOW_TO_SUBMIT.md", zp, "1 export part(s) skipped", "bob"):
                self.assertIn(want, body, want)
            # conditions.csv's sample-to-group annotations never go out
            for leak in ("HeLa_ctrl_01", "HeLa_trt_01"):
                self.assertNotIn(leak, body)

    def test_issues_recorded_by_report_issue_are_counted(self):
        with tempfile.TemporaryDirectory() as d:
            _, out = search_session(d)
            today = time.strftime("%Y-%m-%d")
            write(os.path.join(d, "issues", f"{today}_alice_sess.md"),
                  "# Skill issues -- alice\n\n<!-- end of header -->\n\n## one [bug, high]\n"
                  "- x\n\n## two [docs, low]\n- y\n")
            with env_patch(d, SLURM_JOB_ID="1", SLURM_JOB_USER="alice"):
                f = ns.add_issues(ns.search_facts(out, exit_code=0))
            self.assertEqual(f["issues"]["entries"], 2)
            self.assertIn("2 skill issue(s)", json.dumps(ns.render(f)))


class NeverBreaksTheCaller(unittest.TestCase):
    payload = {"text": "t"}

    def test_no_webhook_is_a_silent_no_op_with_one_info_line(self):
        with tempfile.TemporaryDirectory() as d:
            _, out = search_session(d)
            r = subprocess.run([PY, NOTIFIER, "search-done", "--out", out, "--exit-code", "0"],
                               capture_output=True, text=True, env=base_env(d), timeout=60)
            self.assertEqual(r.returncode, 0)
            self.assertEqual(r.stdout, "")
            # one line for Slack (the other job-end lines are the run log and FRAN's result),
            # in public words: a collaborator's job log must not describe the Core's internals
            lines = [ln for ln in r.stderr.splitlines() if "not configured" in ln]
            self.assertEqual(len(lines), 1, r.stderr)
            self.assertIn("off: Core notification not configured for this user", lines[0])
            for internal in ("quobyte", "proteomics-grp", "no_group_webhook", ".config"):
                self.assertNotIn(internal, r.stderr)

    def test_http_500_returns_false_and_raises_nothing(self):
        m = Mock("500")
        try:
            ok, detail = ns.post(self.payload, m.url)
        finally:
            m.close()
        self.assertFalse(ok)
        self.assertIn("HTTP 500", detail)
        self.assertNotIn(SECRET_BIT, detail)
        self.assertNotIn(m.url, detail)

    def test_a_hung_slack_times_out_quickly(self):
        m = Mock("slow", delay=6)
        try:
            t = time.time()
            ok, detail = ns.post(self.payload, m.url, timeout=1)
            took = time.time() - t
        finally:
            m.close()
        self.assertFalse(ok)
        self.assertLess(took, 1 + ns._DEADLINE_PAD + 1)
        self.assertNotIn(SECRET_BIT, detail)

    def test_dns_failure_returns_false(self):
        def no_dns(*a, **k):
            raise socket.gaierror(-2, "Name or service not known")
        with mock.patch("socket.getaddrinfo", side_effect=no_dns), \
                mock.patch.dict(os.environ, {"https_proxy": "", "HTTPS_PROXY": "",
                                             "no_proxy": "*"}):
            ok, detail = ns.post(self.payload, FAKE_SLACK, timeout=2)
        self.assertFalse(ok)
        self.assertIn("URLError", detail)
        self.assertNotIn(SECRET_BIT, detail)

    def test_a_dns_lookup_that_hangs_is_cut_off_by_the_deadline(self):
        def hang(*a, **k):
            time.sleep(30)
        with mock.patch("socket.getaddrinfo", side_effect=hang):
            t = time.time()
            ok, detail = ns.post(self.payload, FAKE_SLACK, timeout=1)
        self.assertFalse(ok)
        self.assertLess(time.time() - t, 1 + ns._DEADLINE_PAD + 1)
        self.assertIn("no answer", detail)

    def test_an_exception_inside_fact_gathering_is_contained(self):
        with tempfile.TemporaryDirectory() as d, env_patch(d, SLURM_JOB_ID="1"):
            with mock.patch.object(ns, "search_facts", side_effect=RuntimeError("boom " + FAKE_SLACK)):
                ok, detail = ns.search_done(d, exit_code=1)
        self.assertFalse(ok)
        self.assertNotIn(SECRET_BIT, detail)


class AlertApi(unittest.TestCase):
    """send_alert() is what other skill scripts call (fran_deposit.py health --alert): a bool,
    nothing on stdout or stderr -- their stdout is a JSON contract."""

    def _call(self, d, **env):
        out, err = io.StringIO(), io.StringIO()
        with env_patch(d, **env), redirect_stdout(out), redirect_stderr(err):
            res = (ns.send_alert("FRAN ingest stalled: 3 staged searches older than 24 h",
                                 title="FRAN health"),
                   ns.send_alert("again", with_status=True))
        return res, out.getvalue(), err.getvalue()

    def test_sent_off_and_failed(self):
        with tempfile.TemporaryDirectory() as d:
            m = Mock("ok")
            try:
                (sent, (sent2, status)), out, err = self._call(d, SKILL_SLACK_WEBHOOK=m.url)
            finally:
                m.close()
            self.assertEqual((sent, sent2, status, out, err), (True, True, "sent", "", ""))
            self.assertEqual(m.bodies[0]["text"], ":warning: FRAN health: FRAN ingest stalled: "
                                                  "3 staged searches older than 24 h")
            (sent, (sent2, status)), out, err = self._call(d)             # no webhook
            self.assertEqual((sent, sent2, out, err), (False, False, "", ""))
            self.assertEqual(status, "off: Core notification not configured for this user")
            (sent, _), out, err = self._call(d, SKILL_SLACK="0")
            self.assertEqual((sent, out, err), (False, "", ""))
            m = Mock("500")
            try:
                (sent, (sent2, status)), out, err = self._call(d, SKILL_SLACK_WEBHOOK=m.url)
            finally:
                m.close()
            self.assertEqual((sent, out, err), (False, "", ""))
            self.assertIn("HTTP 500", status)
            self.assertNotIn(SECRET_BIT, status)


class TheUrlNeverLeaks(unittest.TestCase):
    def _cli(self, d, url, *args):
        r = subprocess.run([PY, NOTIFIER, *args], capture_output=True, text=True, timeout=90,
                           env=base_env(d, SKILL_SLACK_WEBHOOK=url, SLURM_JOB_ID="9"))
        return r

    def test_cli_output_never_carries_the_url(self):
        with tempfile.TemporaryDirectory() as d:
            _, out = search_session(d)
            for behaviour in ("ok", "500"):
                m = Mock(behaviour)
                try:
                    for args in (["search-done", "--out", out, "--exit-code", "1"],
                                 ["--test"], ["--test", "--dry-run"],
                                 ["search-done", "--out", out, "--dry-run"]):
                        r = self._cli(d, m.url, *args)
                        for text in (r.stdout, r.stderr):
                            self.assertNotIn(SECRET_BIT, text, (behaviour, args))
                            self.assertNotIn(m.url, text)
                finally:
                    m.close()
            # a dead port: connection refused
            r = self._cli(d, "http://127.0.0.1:9/services/T/B/" + SECRET_BIT, "--test")
            self.assertEqual(r.returncode, 1)
            self.assertTrue(r.stdout.startswith("not sent"), r.stdout)
            self.assertNotIn(SECRET_BIT, r.stdout + r.stderr)

    def test_test_message_prints_only_sent(self):
        with tempfile.TemporaryDirectory() as d:
            m = Mock("ok")
            try:
                r = self._cli(d, m.url, "--test")
            finally:
                m.close()
            self.assertEqual((r.returncode, r.stdout.strip()), (0, "sent"))
            self.assertEqual(len(m.bodies), 1)
            self.assertIn("notifications work", m.bodies[0]["text"])

    def test_dry_run_prints_the_payload_and_sends_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            _, out = search_session(d)
            m = Mock("ok")
            try:
                r = self._cli(d, m.url, "search-done", "--out", out, "--dry-run")
            finally:
                m.close()
            self.assertEqual(m.bodies, [])
            self.assertIn("Search finished", json.loads(r.stdout)["text"])
            self.assertIn("would send via $SKILL_SLACK_WEBHOOK", r.stderr)

    def test_scrub(self):
        self.assertEqual(ns._scrub(f"err {FAKE_SLACK} x"), "err <webhook> x")
        self.assertEqual(ns._scrub("a secret-url b", "secret-url"), "a <webhook> b")

    def test_generated_job_scripts_carry_no_url(self):
        """A webhook in the generating environment must not end up in any job it writes."""
        with tempfile.TemporaryDirectory() as d, \
                env_patch(d, SKILL_SLACK_WEBHOOK=FAKE_SLACK):
            job = os.path.join(d, "job.sh")
            run_search.emit_sbatch(job, "echo search", d, 4, job="diann_search",
                                   partition="high", account="acct", qos="q")
            text = read(job)
            self.assertIn("notify_slack.py", text)
            self.assertNotIn("hooks.slack.com", text)
            self.assertNotIn(SECRET_BIT, text)


class Resolution(unittest.TestCase):
    def test_order_env_then_user_file_then_group_file(self):
        with tempfile.TemporaryDirectory() as d:
            user = os.path.join(d, ".config", "ucdavis-proteomics", "slack_webhook")
            group = os.path.join(d, "group_webhook")
            write(user, "http://127.0.0.1:1/user\n")
            write(group, "http://127.0.0.1:1/group\n")
            with env_patch(d, SKILL_SLACK_WEBHOOK="http://127.0.0.1:1/env"), \
                    mock.patch.object(ns, "GROUP_FILE", group):
                self.assertEqual(ns.resolve_webhook(), ("http://127.0.0.1:1/env",
                                                        "$SKILL_SLACK_WEBHOOK"))
            with env_patch(d), mock.patch.object(ns, "GROUP_FILE", group):
                self.assertEqual(ns.resolve_webhook()[0], "http://127.0.0.1:1/user")
                os.remove(user)
                self.assertEqual(ns.resolve_webhook(), ("http://127.0.0.1:1/group", group))

    def test_a_non_slack_url_switches_notifications_off(self):
        with tempfile.TemporaryDirectory() as d:
            with env_patch(d, SKILL_SLACK_TEST_LOOPBACK="0",
                           SKILL_SLACK_WEBHOOK="https://evil.example.com/collect"):
                url, why = ns.resolve_webhook()
            self.assertIsNone(url)
            self.assertIn("not a https://hooks.slack.com/ URL", why)
            self.assertNotIn("evil", why)
            # a bad user file does not fall through to the group file
            write(os.path.join(d, ".config", "ucdavis-proteomics", "slack_webhook"),
                  "https://evil.example.com/x")
            group = os.path.join(d, "g")
            write(group, FAKE_SLACK)
            with env_patch(d, SKILL_SLACK_TEST_LOOPBACK="0"), \
                    mock.patch.object(ns, "GROUP_FILE", group):
                self.assertIsNone(ns.resolve_webhook()[0])

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0, "root reads anything")
    def test_an_unreadable_group_file_says_why(self):
        with tempfile.TemporaryDirectory() as d:
            gdir = os.path.join(d, "cfg")
            os.makedirs(gdir)
            group = os.path.join(gdir, "skill_slack_webhook")
            write(group, FAKE_SLACK)
            os.chmod(gdir, 0)                      # like a 2750 dir to someone outside the group
            try:
                with env_patch(d), mock.patch.object(ns, "GROUP_FILE", group):
                    url, why = ns.resolve_webhook()
            finally:
                os.chmod(gdir, stat.S_IRWXU)
            self.assertIsNone(url)
            self.assertIn("not readable", why)
            self.assertNotIn(SECRET_BIT, why)

    def test_skill_slack_0_disables_everything(self):
        with tempfile.TemporaryDirectory() as d:
            m = Mock("ok")
            try:
                with env_patch(d, SKILL_SLACK="0", SKILL_SLACK_WEBHOOK=m.url):
                    ok, detail = ns.send_test()
                    script = ns.wrap_job_script("#!/bin/bash\necho hi\n", d, final=True,
                                                time_limit_h=1, stage="search")
            finally:
                m.close()
            self.assertFalse(ok)
            self.assertEqual(detail, "off: SKILL_SLACK=0")
            self.assertEqual(m.bodies, [])
            # the job-end hook stays (run log + FRAN); only the post is switched off
            self.assertIn("notify_slack.py", script)
            self.assertIn("--no-slack", script)
            self.assertIn("echo hi", script)


class JobTrap(unittest.TestCase):
    """The generated wrapper, run with bash exactly as SLURM runs a batch script, from a copy of
    the scripts with a fake record_run.py (and, where named, a fake fran_deposit.py)."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.sess, self.out = search_session(self.d)
        self.m = Mock("ok")
        self.log = os.path.join(self.d, "order.log")

    def tearDown(self):
        self.m.close()
        shutil.rmtree(self.d, ignore_errors=True)

    def _job(self, body, final=True, minutes=60, record_run=FAKE_RECORD_RUN, fran=None,
             slack=True, fran_on=True, fran_guarded=True, fran_name=None, qc=None, **env):
        """fran_guarded=True by default: the DIA-NN routes, the only ones that stage from a job."""
        sd = scripts_copy(self.d, record_run=record_run, fran=fran)
        with mock.patch.dict(os.environ, {"FRAN_DEPOSIT": "", "SKILL_SLACK": ""}):
            script = ns.wrap_job_script("#!/bin/bash -l\n#SBATCH --job-name=t\n" + body,
                                        self.out, final=final, time_limit_h=minutes / 60,
                                        stage="search", slack=slack, fran=fran_on,
                                        fran_guarded=fran_guarded, fran_name=fran_name, qc=qc,
                                        scripts_dir=sd)
        path = os.path.join(self.d, "job.sh")
        write(path, script)
        chk = subprocess.run(["bash", "-n", path], capture_output=True, text=True)
        self.assertEqual(chk.returncode, 0, chk.stderr)
        e = base_env(self.d, SKILL_SLACK_WEBHOOK=self.m.url, SLURM_JOB_ID="123",
                     FAKE_ORDER_LOG=self.log)
        e.update(env)
        return path, e

    def _run(self, body, **kw):
        path, e = self._job(body, **kw)
        return subprocess.run(["bash", path], capture_output=True, text=True, env=e, timeout=120)

    def _calls(self, who):
        return [c for c in order_log(self.log) if c["who"] == who]

    def test_exit_codes_are_kept_and_reported(self):
        for body, rc, word in (("set -euo pipefail\necho ok\n", 0, "Search finished"),
                               ("set -euo pipefail\nexit 7\n", 7, "exit 7"),
                               ("set -euo pipefail\nfalse\necho NOT REACHED\n", 1, "exit 1")):
            before = len(self.m.bodies)
            r = self._run(body)
            self.assertEqual(r.returncode, rc, r.stderr)
            self.assertNotIn("NOT REACHED", r.stdout)
            self.assertEqual(len(self.m.bodies), before + 1, r.stderr)
            self.assertIn(word, self.m.bodies[-1]["text"])
            self.assertNotIn(SECRET_BIT, r.stdout + r.stderr)

    def test_success_stages_then_logs_then_posts(self):
        """On success FRAN goes first, so the run record carries the real FRAN outcome (the stage
        receipt) instead of "not staged (yet)"; then the run log; Slack last."""
        r = self._run("true\n", fran=FAKE_FRAN)
        self.assertEqual(r.returncode, 0, r.stderr)
        calls = order_log(self.log)
        self.assertEqual([c["who"] for c in calls], ["fran_deposit", "record_run"], calls)
        meta = os.path.join(self.sess, "input", "search.fasta.meta.json")
        self.assertEqual(calls[0]["argv"], ["stage", "--out", self.out, "--fasta-meta", meta])
        self.assertEqual(calls[1]["argv"], ["search-done", "--out", self.out, "--status",
                                            "completed", "--exit-code", "0",
                                            "--session", self.sess])
        self.assertGreaterEqual(calls[1]["t"], calls[0]["t"])
        self.assertEqual(len(self.m.bodies), 1)
        self.assertGreater(self.m.times[0], calls[1]["t"])          # Slack last
        body = json.dumps(self.m.bodies[0])
        self.assertIn("*Run log:* yes", body)
        self.assertIn("*Staged for FRAN:* yes", body)
        # stage's result, one line, in the job log
        self.assertIn('fran_deposit stage: {"staged":true,"reason":"ok"', r.stderr)
        self.assertIn("run log: yes; FRAN: yes; sent", r.stderr)

    def test_failure_is_logged_but_never_staged(self):
        r = self._run("exit 7\n", fran=FAKE_FRAN)
        self.assertEqual(r.returncode, 7)
        self.assertEqual(self._calls("fran_deposit"), [])
        rr = self._calls("record_run")
        self.assertEqual(len(rr), 1)
        self.assertEqual(rr[0]["argv"][3:7], ["--status", "failed", "--exit-code", "7"])
        self.assertNotIn("Staged for FRAN", json.dumps(self.m.bodies[-1]))

    def test_staged_but_fran_ingest_unhealthy_is_said(self):
        r = self._run("true\n", fran=FAKE_FRAN, FAKE_FRAN_MODE="stuck")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("*Staged for FRAN:* yes, but FRAN's ingest is unhealthy: FRAN ingest "
                      "stuck: nothing ingested since 2026-09-17", json.dumps(self.m.bodies[-1]))
        self.assertIn('"health_warning":"FRAN ingest stuck', r.stderr)   # the job log line

    def test_a_crashing_stage_does_not_change_the_exit_code(self):
        r = self._run("true\n", fran=FAKE_FRAN, FAKE_FRAN_MODE="crash")
        self.assertEqual(r.returncode, 0, r.stderr)
        body = json.dumps(self.m.bodies[-1])
        self.assertIn("*Staged for FRAN:* error: fran_deposit.py stage exit 1: "
                      "RuntimeError: drop dir exploded", body)

    def test_the_real_stage_hands_a_finished_search_to_fran(self):
        drop = os.path.join(self.d, "fran_drop")
        r = self._run("true\n", fran=None, FRAN_DEPOSIT="", FRAN_DROP_DIR=drop)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn('fran_deposit stage: {"staged":true', r.stderr)
        self.assertIn("*Staged for FRAN:* yes", json.dumps(self.m.bodies[-1]))
        entries = os.listdir(drop)
        self.assertEqual(len(entries), 1, entries)
        man = json.loads(read(os.path.join(drop, entries[0], "fran_manifest.json")))
        self.assertEqual((man["engine"], man["organism"]), ("diann", "Homo sapiens"))
        self.assertTrue(os.path.islink(os.path.join(drop, entries[0], "report.parquet")))

    def test_fran_deposit_off_is_reported_as_skipped(self):
        r = self._run("true\n", fran=None)                      # base_env: FRAN_DEPOSIT=off
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("*Staged for FRAN:* skipped: off for this search",
                      json.dumps(self.m.bodies[-1]))
        self.assertFalse(os.path.exists(os.path.join(self.d, "fran_drop")))

    def test_record_run_absent_is_skipped_cleanly(self):
        r = self._run("true\n", record_run=False, fran=FAKE_FRAN)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([c["who"] for c in order_log(self.log)], ["fran_deposit"])
        body = json.dumps(self.m.bodies[-1])
        self.assertNotIn("Run log", body)
        self.assertIn("*Staged for FRAN:* yes", body)

    def test_a_crashing_record_run_does_not_change_the_exit_code(self):
        for body, rc in (("true\n", 0), ("exit 6\n", 6)):
            r = self._run(body, fran=FAKE_FRAN, FAKE_RECORD_MODE="crash")
            self.assertEqual(r.returncode, rc, r.stderr)
            self.assertIn("*Run log:* error -- record_run.py exit 1: record_run exploded",
                          json.dumps(self.m.bodies[-1]))

    def test_no_slack_still_logs_and_stages(self):
        r = self._run("true\n", fran=FAKE_FRAN, slack=False)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self.m.bodies, [])
        self.assertEqual([c["who"] for c in order_log(self.log)], ["fran_deposit", "record_run"])
        self.assertIn("Slack: off for this job", r.stderr)

    def test_the_report_guard_still_decides(self):
        """must_exist() failing = the job fails with 1, and the post says failed."""
        import diann_parallel as dp
        guard = dp.must_exist(os.path.join(self.out, "never_written.parquet"), "the report")
        r = self._run("set -euo pipefail\necho searched\n" + guard + "\necho NOT REACHED\n",
                      fran=FAKE_FRAN)
        self.assertEqual(r.returncode, 1)
        self.assertIn("FAILED: DIA-NN exited 0 but did not write the report", r.stderr)
        self.assertIn("FAILED", self.m.bodies[-1]["text"])
        self.assertEqual(self._calls("fran_deposit"), [])

    def test_an_unreachable_slack_does_not_change_the_exit_code(self):
        r = self._run("exit 3\n", SKILL_SLACK_WEBHOOK="http://127.0.0.1:9/dead")
        self.assertEqual(r.returncode, 3)
        self.assertIn("not sent", r.stderr)

    def test_a_non_final_job_acts_only_on_failure(self):
        r = self._run("true\n", final=False, fran=FAKE_FRAN)
        self.assertEqual(r.returncode, 0)
        self.assertEqual((self.m.bodies, order_log(self.log)), ([], []))
        r = self._run("exit 4\n", final=False, fran=FAKE_FRAN)
        self.assertEqual(r.returncode, 4)
        self.assertEqual(len(self.m.bodies), 1)
        self.assertEqual([c["who"] for c in order_log(self.log)], ["record_run"])

    def test_an_array_reports_its_first_failure_only(self):
        for task in range(3):
            r = self._run("exit 2\n", final=False, SLURM_ARRAY_JOB_ID="555",
                          SLURM_ARRAY_TASK_ID=str(task), SLURM_JOB_ID=str(556 + task))
            self.assertEqual(r.returncode, 2)
        self.assertEqual(len(self.m.bodies), 1)
        self.assertEqual(len(self._calls("record_run")), 1)
        # mkdir, not an O_EXCL file: the primitive measured atomic across HIVE nodes on Quobyte
        self.assertTrue(os.path.isdir(os.path.join(self.out, ".slack_failed_555")))

    def test_outside_a_slurm_job_nothing_is_done(self):
        path, e = self._job("exit 5\n", fran=FAKE_FRAN)
        e.pop("SLURM_JOB_ID")
        r = subprocess.run(["bash", path], capture_output=True, text=True, env=e, timeout=60)
        self.assertEqual(r.returncode, 5)
        self.assertEqual((self.m.bodies, order_log(self.log)), ([], []))

    def _term(self, minutes):
        """SIGTERM to the batch shell ONLY, as SLURM does at the time limit (HIVE, measured)."""
        path, e = self._job("set -euo pipefail\necho started\nsleep 60\necho NOT REACHED\n",
                            minutes=minutes, fran=FAKE_FRAN)
        errf = os.path.join(self.d, "job.err")
        with open(errf, "w") as eh:
            p = subprocess.Popen(["bash", path], stdout=subprocess.DEVNULL, stderr=eh,
                                 env=e, start_new_session=True)
            time.sleep(1.5)
            t = time.time()
            p.send_signal(signal.SIGTERM)
            rc = p.wait(timeout=45)
            took = time.time() - t
        try:
            os.killpg(p.pid, signal.SIGKILL)       # the orphaned `sleep`, as SLURM's KillWait would
        except OSError:
            pass
        return rc, took, read(errf)

    def test_time_limit_term_is_handled_at_once_logged_and_posted(self):
        rc, took, err = self._term(minutes=1)      # 1-minute limit: any TERM is the limit
        self.assertEqual(rc, -signal.SIGTERM, err)  # still dies by the signal
        self.assertLess(took, 20)                   # not after the 60 s foreground sleep
        self.assertEqual(len(self.m.bodies), 1, err)
        self.assertIn("time limit", self.m.bodies[0]["text"])
        rr = self._calls("record_run")
        self.assertEqual(rr[0]["argv"][3:7], ["--status", "failed", "--exit-code", "143"])
        self.assertEqual(self._calls("fran_deposit"), [])

    def test_term_before_the_limit_does_nothing(self):
        rc, took, err = self._term(minutes=600)    # cancelled, or preempted and requeued
        self.assertEqual(rc, -signal.SIGTERM)
        self.assertLess(took, 20)
        self.assertEqual((self.m.bodies, order_log(self.log)), ([], []))
        self.assertIn("nothing done: stopped by SIGTERM", err)


class Generators(unittest.TestCase):
    """Every route's LAST job reports; a job that others wait on reports only a failure."""
    import test_single_shot_mass_accuracy as _ss
    # Borrow the harness, not the tests (as test_two_job_submit.py does).
    _setup = _ss.SingleShotMassAccTests._setup
    _env = _ss.SingleShotMassAccTests._env

    def test_emit_sbatch_final_and_fail_only(self):
        with tempfile.TemporaryDirectory() as d, env_patch(d):
            for notify, fail_only in (("final", False), ("fail", True)):
                job = os.path.join(d, f"{notify}.sh")
                run_search.emit_sbatch(job, "echo search", d, 4, job="x", partition="high",
                                       account="a", qos="q", notify=notify, stage="s")
                text = read(job)
                self.assertIn("search-done --from-job", text)
                self.assertEqual("--fail-only" in text, fail_only)
                self.assertTrue(text.startswith("#!/bin/bash -l\n#SBATCH --job-name=x"))
                self.assertEqual(subprocess.run(["bash", "-n", job]).returncode, 0)

    def test_only_guarded_routes_may_stage_for_fran(self):
        """emit_sbatch's default is unguarded (Sage, AlphaDIA, FragPipe, single Radiant): those
        jobs never stage; the Radiant chain's step 3 counts files, not a complete report."""
        with tempfile.TemporaryDirectory() as d, env_patch(d, FRAN_DEPOSIT=""):
            for guarded in (False, True):
                job = os.path.join(d, f"g{guarded}.sh")
                run_search.emit_sbatch(job, "echo search", d, 4, job="x", partition="high",
                                       account="a", qos="q", fran_guarded=guarded)
                self.assertEqual("--fran-stage" in read(job), guarded)
            files = []
            for i in range(2):
                files.append(os.path.join(d, f"run{i}.mzML"))
                write(files[-1], "mzml")
            for name in ("lib.tsv", "db.fasta", "cfg.radiantConfig"):
                write(os.path.join(d, name), "x")
            write(os.path.join(d, "list.txt"), "\n".join(files) + "\n")
            out = os.path.join(d, "radiant_out")
            g = subprocess.run([PY, os.path.join(SCRIPTS, "radiant_parallel.py"), "--runtime",
                                "apptainer", "--image", os.path.join(d, "r.sif"), "--raw-list",
                                os.path.join(d, "list.txt"), "--fasta", os.path.join(d, "db.fasta"),
                                "--config", os.path.join(d, "cfg.radiantConfig"), "--out", out,
                                "--library", os.path.join(d, "lib.tsv"), "--partition", "high",
                                "--account", "a", "--qos", "q"],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(g.returncode, 0, g.stderr)
            steps = glob.glob(os.path.join(out, "radiant_parallel", "*.sbatch"))
            self.assertTrue(steps)
            for st in steps:
                text = read(st)
                self.assertIn("search-done --from-job", text, st)
                self.assertNotIn("--fran-stage", text, st)

    def test_run_search_two_job_route_and_no_notify(self):
        for extra, no_slack, no_fran in (([], False, False), (["--no-notify"], True, False),
                                         (["--no-fran"], False, True),
                                         (["--no-notify", "--no-fran"], True, True)):
            with tempfile.TemporaryDirectory() as d:
                raws, fasta, cfg, tools, bundle = self._setup(d)
                out = os.path.join(d, "out")
                env = self._env(d)
                env.pop("SKILL_SLACK", None)            # job_env's opt-outs would be baked in
                env.pop("FRAN_DEPOSIT", None)
                p = subprocess.run([PY, os.path.join(SCRIPTS, "run_search.py"),
                                    "--tools", tools, "--bundle", bundle, "--params", cfg,
                                    "--fasta", fasta, "--out", out, "--files", *raws,
                                    "--engine", "diann", "--threads", "8", "--sbatch", "job.sh",
                                    *extra], capture_output=True, text=True, env=env, cwd=d,
                                   timeout=240)
                self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
                lib, srch = (read(os.path.join(d, f"job_{s}.sh"))
                             for s in ("1_lib", "2_search"))
                for text in (lib, srch):            # the hook is always there (run log, FRAN)
                    self.assertIn("notify_slack.py search-done --from-job", text)
                    self.assertEqual("--no-slack" in text, no_slack)
                    self.assertEqual("--no-fran" in text, no_fran)
                self.assertIn("--fail-only", lib)
                self.assertNotIn("--fail-only", srch)
                # the search job ends with report_guard: a guarded route, so it may stage
                self.assertEqual("--fran-stage" in srch, not no_fran)
                prov = json.loads(read(os.path.join(out, "search_provenance.json")))
                self.assertEqual(prov["bundle"], os.path.abspath(bundle))
                self.assertEqual(prov["job_end_hook"],
                                 {"run_log": "on", "slack": "off" if no_slack else "on",
                                  "fran": "off" if no_fran else "stage",
                                  "fran_name": None, "qc": None})

    def test_run_search_records_the_name_and_the_qc_decision(self):
        for extra, qc in ((["--fran-name", "Lumos QC 2026-09-24", "--qc"], True),
                          (["--fran-name", "Lumos QC 2026-09-24", "--not-qc"], False)):
            with tempfile.TemporaryDirectory() as d:
                raws, fasta, cfg, tools, bundle = self._setup(d)
                out = os.path.join(d, "out")
                env = self._env(d)
                env.pop("SKILL_SLACK", None)
                env.pop("FRAN_DEPOSIT", None)
                p = subprocess.run([PY, os.path.join(SCRIPTS, "run_search.py"),
                                    "--tools", tools, "--bundle", bundle, "--params", cfg,
                                    "--fasta", fasta, "--out", out, "--files", *raws,
                                    "--engine", "diann", "--threads", "8", "--sbatch", "job.sh",
                                    *extra], capture_output=True, text=True, env=env, cwd=d,
                                   timeout=240)
                self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
                srch = read(os.path.join(d, "job_2_search.sh"))
                self.assertIn("--fran-name 'Lumos QC 2026-09-24'", srch)
                self.assertIn("--qc" if qc else "--not-qc", srch)
                self.assertEqual("--fran-stage" in srch, not qc)     # a QC run never stages
                prov = json.loads(read(os.path.join(out, "search_provenance.json")))
                self.assertEqual(prov["job_end_hook"],
                                 {"run_log": "on", "slack": "on",
                                  "fran": "off" if qc else "stage",
                                  "fran_name": "Lumos QC 2026-09-24", "qc": qc})

    def test_diann_chain_step5_final_the_rest_fail_only(self):
        import test_parallel_routing_window as prw
        with tempfile.TemporaryDirectory() as d, env_patch(d, FRAN_DEPOSIT=""):
            raws, fasta = prw._inputs(d)
            cfg = prw._cfg(d, prw.PINNED_MA + "--window 7\n")
            p, out = prw._generate(d, cfg, raws, fasta, "--partition", "high",
                                   "--account", "a", "--qos", "q")
            self.assertEqual(p.returncode, 0, p.stderr)
            steps = sorted(f for f in os.listdir(out) if f.endswith(".sbatch"))
            self.assertIn("step5_report.sbatch", steps)
            for s in steps:
                text = read(os.path.join(out, s))
                self.assertIn("search-done --from-job", text, s)
                self.assertEqual("--fail-only" in text, s != "step5_report.sbatch", s)
                # step 5's quant count is the chain's completeness guard: only it stages
                self.assertEqual("--fran-stage" in text, s == "step5_report.sbatch", s)
                self.assertEqual(subprocess.run(["bash", "-n", os.path.join(out, s)]).returncode, 0)
            # step 4 keeps its own private-library clean-up, inside the work
            s4 = read(os.path.join(out, "step4_finalpass.sbatch"))
            self.assertLess(s4.index("trap _job_end_term TERM"),
                            s4.index("trap 'rm -rf \"$LIBPRIV\"' EXIT"))
            self.assertIn("submit.sh", os.listdir(out))
            self.assertNotIn("notify_slack", read(os.path.join(out, "submit.sh")))
            d2 = os.path.join(d, "second")
            os.makedirs(d2)
            p2, out2 = prw._generate(d2, cfg, raws, fasta, "--partition", "high",
                                     "--account", "a", "--qos", "q", "--no-notify")
            self.assertEqual(p2.returncode, 0, p2.stderr)
            for s in (f for f in os.listdir(out2) if f.endswith(".sbatch")):
                self.assertIn("--no-slack", read(os.path.join(out2, s)))
            d3 = os.path.join(d, "third")
            os.makedirs(d3)
            p3, out3 = prw._generate(d3, cfg, raws, fasta, "--partition", "high",
                                     "--account", "a", "--qos", "q", "--no-fran")
            self.assertEqual(p3.returncode, 0, p3.stderr)
            s5 = read(os.path.join(out3, "step5_report.sbatch"))
            self.assertIn("--no-fran", s5)
            self.assertNotIn("--fran-stage", s5)
            d4 = os.path.join(d, "fourth")
            os.makedirs(d4)
            p4, out4 = prw._generate(d4, cfg, raws, fasta, "--partition", "high",
                                     "--account", "a", "--qos", "q",
                                     "--fran-name", "Poplar drought", "--not-qc")
            self.assertEqual(p4.returncode, 0, p4.stderr)
            s5 = read(os.path.join(out4, "step5_report.sbatch"))
            for want in ("--fran-stage", "--not-qc", "--fran-name 'Poplar drought'"):
                self.assertIn(want, s5)


class Finalize(unittest.TestCase):
    """finalize, run from a copy of the scripts with a fake record_run.py (see scripts_copy)."""

    def _finalize(self, sd, env, *extra, record_run=FAKE_RECORD_RUN):
        root = os.path.dirname(sd)
        scripts = scripts_copy(os.path.join(root, "copy"), record_run=record_run)
        env = dict(env, FAKE_ORDER_LOG=os.path.join(root, "order.log"))
        return subprocess.run([PY, os.path.join(scripts, "session.py"), "finalize", "--dir", sd,
                               *extra], capture_output=True, text=True, env=env, timeout=240)

    def test_the_run_is_logged_before_the_post(self):
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            m = Mock("ok")
            try:
                r = self._finalize(p["session_dir"], base_env(d, SKILL_SLACK_WEBHOOK=m.url),
                                   "--zip")
            finally:
                m.close()
            self.assertEqual(r.returncode, 0, r.stderr)
            res = json.loads(r.stdout)
            self.assertIs(res["run_log"]["logged"], True)
            calls = order_log(os.path.join(os.path.dirname(p["session_dir"]), "order.log"))
            self.assertEqual([c["argv"] for c in calls],
                             [["analysis-done", "--timeout", "300", "--session",
                               p["session_dir"]]])
            self.assertGreater(m.times[0], calls[0]["t"])            # logged, then posted
            self.assertIn("*Run log:* yes", json.dumps(m.bodies[0]))
            # both outcomes, run log first, are the last lines of MANIFEST.txt -- and in the zip
            lines = read(p["manifest_txt"]).splitlines()
            self.assertTrue(lines[-2].startswith("[OK]      Core run log"), lines[-2:])
            self.assertTrue(lines[-2].endswith("-- logged"), lines[-2])
            self.assertNotIn("/runs/", lines[-2])          # never a path: the zip is shared
            self.assertTrue(lines[-1].startswith("[OK]      Slack notification"), lines[-1])
            self.assertEqual(res["zip_manifest"], "added")
            base = os.path.basename(p["session_dir"])
            with zipfile.ZipFile(res["zip"]) as z:
                zipped = z.read(f"{base}/MANIFEST.txt").decode()
            self.assertIn("Core run log", zipped)

    def test_record_run_absent_or_crashing_never_fails_finalize(self):
        import test_deposit_package as tdp
        for record_run, env_extra, want in ((False, {}, None),
                                            (FAKE_RECORD_RUN, {"FAKE_RECORD_MODE": "crash"},
                                             False)):
            with tempfile.TemporaryDirectory() as d:
                p = tdp.dia_session(d)
                r = self._finalize(p["session_dir"], base_env(d, **env_extra), "--zip",
                                   record_run=record_run)
                self.assertEqual(r.returncode, 0, r.stderr)
                res = json.loads(r.stdout)
                self.assertTrue(os.path.isfile(res["zip"]))
                line = [ln for ln in read(p["manifest_txt"]).splitlines()
                        if "Core run log" in ln]
                self.assertEqual(len(line), 1, line)
                if want is None:                   # not installed: a notice, not a failure
                    self.assertIsNone(res["run_log"])
                    self.assertTrue(line[0].startswith("[INFO]"), line)
                    self.assertTrue(line[0].endswith("-- not configured for this user"), line)
                else:                              # attempted and failed: [SKIPPED]
                    self.assertTrue(line[0].startswith("[SKIPPED]"), line)
                    self.assertIs(res["run_log"]["logged"], want)
                    self.assertIn("record_run exploded", res["run_log"]["detail"])
                    self.assertIn("record_run.py exit 1: record_run exploded", line[0])

    def test_a_real_record_run_failure_is_skipped_with_its_error(self):
        """record_run.py always exits 0; a failed record is {"recorded": false, "reason":
        "error", "detail": "<exception>"}. That is [SKIPPED] with the text, never [INFO]."""
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            runs = os.path.join(d, "runs")
            os.makedirs(runs)
            write(os.path.join(runs, "sessions"), "a FILE where the sessions/ directory goes")
            r = self._finalize(p["session_dir"], base_env(d, SKILL_RUNS_DIR=runs), "--zip",
                               record_run=None)                       # the REAL record_run.py
            self.assertEqual(r.returncode, 0, r.stderr)
            res = json.loads(r.stdout)
            self.assertIs(res["run_log"]["error"], True, res["run_log"])
            line = [ln for ln in read(p["manifest_txt"]).splitlines() if "Core run log" in ln]
            self.assertEqual(len(line), 1, line)
            self.assertTrue(line[0].startswith("[SKIPPED]"), line)
            self.assertIn("-- error: ", line[0])
            self.assertRegex(line[0], r"(NotADirectoryError|FileExistsError|OSError|Errno)", line)

    def test_no_notify_still_logs_the_run(self):
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            r = self._finalize(p["session_dir"], base_env(d), "--no-notify")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIs(json.loads(r.stdout)["run_log"]["logged"], True)

    def _manifest_slack(self, path):
        with open(path) as fh:
            return [ln for ln in fh.read().splitlines() if "notification" in ln]

    def test_sent_is_recorded_in_manifest_and_in_the_zip(self):
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            m = Mock("ok")
            try:
                r = self._finalize(p["session_dir"], base_env(d, SKILL_SLACK_WEBHOOK=m.url),
                                   "--zip")
            finally:
                m.close()
            self.assertEqual(r.returncode, 0, r.stderr)
            res = json.loads(r.stdout)
            self.assertEqual(res["slack"], {"sent": True, "level": "OK", "detail": "sent"})
            self.assertEqual(len(m.bodies), 1)
            self.assertIn("Analysis complete", m.bodies[0]["text"])
            self.assertIn(res["zip"], json.dumps(m.bodies[0]))
            line = self._manifest_slack(p["manifest_txt"])
            self.assertEqual(len(line), 1)
            self.assertTrue(line[0].startswith("[OK]"), line)
            base = os.path.basename(p["session_dir"])
            with zipfile.ZipFile(res["zip"]) as z:
                self.assertEqual(z.namelist().count(f"{base}/MANIFEST.txt"), 1)
                self.assertIn("Slack notification", z.read(f"{base}/MANIFEST.txt").decode())
            self.assertNotIn(SECRET_BIT, r.stdout + r.stderr)

    def test_off_and_no_notify_are_info_not_skipped(self):
        """Not configured / opted out is a notice ([INFO], not relayed as a missing export part),
        in words that name no path and no group -- the zip goes to collaborators."""
        import test_deposit_package as tdp
        for env_extra, flag, part, note in (
                ({}, [], "Core notification", "not configured for this user"),
                ({"SKILL_SLACK": "0"}, [], "Slack notification (Core channel)",
                 "off (SKILL_SLACK=0)"),
                ({}, ["--no-notify"], "Slack notification (Core channel)", "off (--no-notify)")):
            with tempfile.TemporaryDirectory() as d:
                p = tdp.dia_session(d)
                r = self._finalize(p["session_dir"], base_env(d, **env_extra), *flag)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertFalse(json.loads(r.stdout)["slack"]["sent"])
                line = self._manifest_slack(p["manifest_txt"])
                self.assertEqual(len(line), 1, line)
                self.assertTrue(line[0].startswith("[INFO]"), line)
                self.assertEqual(line[0][10:].split(" -- ")[0].strip(), part)
                self.assertEqual(line[0].split(" -- ", 1)[1], note)
                self.assertEqual(json.loads(r.stdout)["slack"]["detail"], note)
                text = read(p["manifest_txt"])
                self.assertFalse(any(ln.startswith("[SKIPPED]") for ln in text.splitlines()[-2:]))
                for internal in ("quobyte", "proteomics-grp", "no_group_webhook", ".config"):
                    self.assertNotIn(internal, text)

    def test_a_500_is_recorded_and_finalize_still_succeeds(self):
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            m = Mock("500")
            try:
                r = self._finalize(p["session_dir"], base_env(d, SKILL_SLACK_WEBHOOK=m.url))
            finally:
                m.close()
            self.assertEqual(r.returncode, 0, r.stderr)
            line = self._manifest_slack(p["manifest_txt"])[0]
            self.assertTrue(line.startswith("[SKIPPED]"), line)     # attempted and failed
            self.assertIn("HTTP 500", line)
            self.assertNotIn(SECRET_BIT, line + r.stdout + r.stderr)


FAKE_HIVE_EXEC = r"""#!/bin/bash
# Stands in for hive_exec.sh: runs the command "on HIVE" -- here, locally, in an environment
# where the Core's webhook (the mock) is configured and the laptop's is not.
echo "$1" > "$RELAY_CMD_LOG"
export SKILL_SLACK_WEBHOOK="$HIVE_SIDE_WEBHOOK"
exec bash -c "$1"
"""


class RelayThroughHive(unittest.TestCase):
    def test_a_laptop_without_a_webhook_posts_through_hive(self):
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            hx = os.path.join(d, "fake_hive_exec.sh")
            write(hx, FAKE_HIVE_EXEC)
            os.chmod(hx, 0o755)
            m = Mock("ok")
            try:
                env = base_env(d, HIVE_USER="carol", HIVE_EXEC=hx, HIVE_SIDE_WEBHOOK=m.url,
                               RELAY_CMD_LOG=os.path.join(d, "relay_cmd"))
                r = subprocess.run([PY, NOTIFIER, "analysis-done", "--session",
                                    p["session_dir"]], capture_output=True, text=True,
                                   env=env, timeout=120)
            finally:
                m.close()
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("sent through HIVE (carol)", r.stderr)
            self.assertEqual(len(m.bodies), 1)
            self.assertIn("Analysis complete", m.bodies[0]["text"])
            self.assertIn("carol", m.bodies[0]["text"])
            cmd = read(os.path.join(d, "relay_cmd"))
            self.assertTrue(cmd.startswith("python3 - relay --facts-b64 "), cmd)
            self.assertNotIn(SECRET_BIT, cmd + r.stdout + r.stderr)
            facts = json.loads(base64.b64decode(cmd.split()[-1]))
            self.assertEqual(facts["kind"], "analysis")

    def test_no_relay_when_opted_out_or_already_on_hive(self):
        with tempfile.TemporaryDirectory() as d:
            hx = os.path.join(d, "fake_hive_exec.sh")
            write(hx, "#!/bin/bash\ntouch \"$RELAY_CMD_LOG\"\n")
            flag = os.path.join(d, "called")
            for extra in ({"SKILL_SLACK": "0"},
                          {"SKILL_CORE_GROUP_DIR": d}):   # the group folder exists: on HIVE
                env = base_env(d, HIVE_USER="carol", HIVE_EXEC=hx, RELAY_CMD_LOG=flag, **extra)
                r = subprocess.run([PY, NOTIFIER, "--test"], capture_output=True, text=True,
                                   env=env, timeout=60)
                self.assertEqual(r.returncode, 1)
                self.assertFalse(os.path.exists(flag), extra)


FAKE_SQUEUE = r"""#!/bin/bash
# squeue -h -j <id> -o %e -> this job's end time; FAKE_SQUEUE_FAIL=1 makes it fail like a
# missing or unreachable squeue.
[ -n "${FAKE_SQUEUE_FAIL:-}" ] && exit 1
echo "$FAKE_SQUEUE_END"
"""


def _race_worker(out, rounds, barrier, q):
    """One array task: at each round, claim the failure report for array `race<r>`."""
    wins = []
    for r in range(rounds):
        os.environ["SLURM_ARRAY_JOB_ID"] = f"race{r}"
        barrier.wait()
        wins.append(ns._array_first_failure(out))
    q.put(wins)


class RunLogOutcomes(unittest.TestCase):
    """record_run.py's JSON, read by its reason codes (it always exits 0)."""

    def _rr(self, payload, kind="search-done"):
        seen = {}

        def fake(argv, timeout):
            seen.update(argv=argv, timeout=timeout)
            return 0, json.dumps(payload), ""
        with mock.patch.object(ns, "_helper", return_value="/s/record_run.py"), \
                mock.patch.object(ns, "_run_helper", side_effect=fake), \
                mock.patch.dict(os.environ, {"RECORD_RUN": ""}):
            res = ns.record_run(kind, out="/o", session="/s", status="failed", exit_code=1)
        return res, seen

    def test_failure_codes_are_errors_with_their_text(self):
        for code in ("error", "timeout", "ssh_failed", "bad_input", "nothing_to_record",
                     "out_not_found", "session_not_found"):
            res, _ = self._rr({"recorded": False, "reason": code, "detail": "the text"})
            self.assertIs(res.get("error"), True, code)
            self.assertEqual(res["detail"], f"{code}: the text")
            self.assertEqual(ns.run_log_manifest(res),
                             ("SKIPPED", "Core run log", f"{code}: the text"))

    def test_not_configured_off_and_unknown_are_notices(self):
        for code in ("not_core_member", "not_on_hive"):
            res, _ = self._rr({"recorded": False, "reason": code,
                               "detail": "/quobyte/proteomics-grp/skill_runs is not writable"})
            self.assertFalse(res.get("error"))
            self.assertEqual(ns.run_log_manifest(res),
                             ("INFO", "Core run log", "not configured for this user"))
        res, _ = self._rr({"recorded": False, "reason": "disabled", "detail": "RECORD_RUN=off"})
        self.assertEqual(ns.run_log_manifest(res), ("INFO", "Core run log", "off (RECORD_RUN=off)"))
        res, _ = self._rr({"recorded": False, "reason": "dry_run", "detail": "nothing written"})
        self.assertEqual(ns.run_log_manifest(res), ("INFO", "Core run log", "not recorded (dry_run)"))
        res, _ = self._rr({"recorded": True, "path": "/quobyte/x"})
        self.assertEqual(ns.run_log_manifest(res), ("OK", "Core run log", "logged"))

    def test_budgets(self):
        _, seen = self._rr({"recorded": True}, kind="search-done")
        self.assertEqual(seen["timeout"], 60)                       # the job hook stays 60 s
        self.assertNotIn("--timeout", seen["argv"])
        _, seen = self._rr({"recorded": True}, kind="analysis-done")
        self.assertEqual(seen["argv"][2:5], ["analysis-done", "--timeout", "300"])
        self.assertGreaterEqual(seen["timeout"], 310)               # room to report its own timeout


class ArrayMarkerRace(unittest.TestCase):
    def test_two_racing_tasks_exactly_one_wins_each_round(self):
        """mkdir is the first-wins primitive (atomic across HIVE nodes on Quobyte, where flock
        is not): two processes released together by a barrier, 100 rounds."""
        import multiprocessing as mp
        ctx = mp.get_context("spawn")
        rounds = 100
        with tempfile.TemporaryDirectory() as d:
            barrier, q = ctx.Barrier(2), ctx.Queue()
            procs = [ctx.Process(target=_race_worker, args=(d, rounds, barrier, q))
                     for _ in range(2)]
            for pr in procs:
                pr.start()
            results = [q.get(timeout=120) for _ in procs]
            for pr in procs:
                pr.join(timeout=30)
            for r in range(rounds):
                self.assertEqual(results[0][r] + results[1][r], 1, f"round {r}")
                self.assertTrue(os.path.isdir(os.path.join(d, f".slack_failed_race{r}")))


class JobEndChoices(unittest.TestCase):
    """--no-fran, routes without a completeness guard, and the time-limit guard on staging."""
    setUp, tearDown = JobTrap.setUp, JobTrap.tearDown
    _job, _run, _calls = JobTrap._job, JobTrap._run, JobTrap._calls

    def _squeue(self, end=None, fail=False):
        fake = os.path.join(self.d, "fakebin")
        os.makedirs(fake, exist_ok=True)
        write(os.path.join(fake, "squeue"), FAKE_SQUEUE)
        os.chmod(os.path.join(fake, "squeue"), 0o755)
        env = {"PATH": fake + os.pathsep + os.environ.get("PATH", "")}
        if fail:
            env["FAKE_SQUEUE_FAIL"] = "1"
        if end is not None:
            env["FAKE_SQUEUE_END"] = time.strftime("%Y-%m-%dT%H:%M:%S",
                                                   time.localtime(time.time() + end))
        return env

    def test_no_fran_at_generation_stages_nothing_in_a_clean_env(self):
        path, e = self._job("true\n", fran=FAKE_FRAN, fran_on=False, FRAN_DEPOSIT="")
        self.assertIn("--no-fran", read(path))
        self.assertNotIn("--fran-stage", read(path))
        r = subprocess.run(["bash", path], capture_output=True, text=True, env=e, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        # stage is still CALLED -- with --skip, which stages nothing but records `opted_out`, so a
        # later `fran_deposit.py backfill` leaves the search alone
        calls = self._calls("fran_deposit")
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["argv"][-1], "--skip")
        self.assertEqual(len(self._calls("record_run")), 1)      # the run log still happens
        self.assertIn("*Staged for FRAN:* skipped: off for this search",
                      json.dumps(self.m.bodies[-1]))

    def test_fran_deposit_off_when_generating_is_baked_into_the_job(self):
        """hive_remote: the generating shell's FRAN_DEPOSIT=off never reaches the job's own
        environment -- so it must be in the script."""
        sd = scripts_copy(self.d, fran=FAKE_FRAN)
        with mock.patch.dict(os.environ, {"FRAN_DEPOSIT": "off", "SKILL_SLACK": ""}):
            script = ns.wrap_job_script("#!/bin/bash -l\ntrue\n", self.out, final=True,
                                        time_limit_h=1, stage="search", fran_guarded=True,
                                        scripts_dir=sd)
        self.assertIn("--no-fran", script)
        path = os.path.join(self.d, "job.sh")
        write(path, script)
        e = base_env(self.d, SKILL_SLACK_WEBHOOK=self.m.url, SLURM_JOB_ID="1",
                     FAKE_ORDER_LOG=self.log, FRAN_DEPOSIT="")      # clean at run time
        r = subprocess.run(["bash", path], capture_output=True, text=True, env=e, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual([c["argv"][-1] for c in self._calls("fran_deposit")], ["--skip"])

    def test_a_route_without_a_completeness_guard_leaves_fran_to_the_agent(self):
        r = self._run("true\n", fran=FAKE_FRAN, fran_guarded=False)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self._calls("fran_deposit"), [])
        self.assertIn("fran: left_to_agent (no completeness guard on this route)", r.stderr)
        self.assertIn("*Staged for FRAN:* skipped: left to the agent (no completeness guard "
                      "on this route)", json.dumps(self.m.bodies[-1]))

    def test_near_the_time_limit_staging_is_left_to_the_agent(self):
        r = self._run("true\n", fran=FAKE_FRAN, **self._squeue(end=60))
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertEqual(self._calls("fran_deposit"), [])
        self.assertIn('"reason":"near_time_limit"', r.stderr)
        self.assertIn("left to the agent (too close to the job's time limit)",
                      json.dumps(self.m.bodies[-1]))

    def test_with_time_to_spare_or_no_squeue_it_stages(self):
        for extra in (self._squeue(end=3600), self._squeue(fail=True)):
            os.path.exists(self.log) and os.remove(self.log)
            r = self._run("true\n", fran=FAKE_FRAN, **extra)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(len(self._calls("fran_deposit")), 1, (extra, r.stderr))

    def test_a_qc_run_is_never_staged_from_the_job(self):
        """--qc at generation: the race it closes is a QC session named without "QC" being
        staged by the hook before the agent's own `stage --name "... QC"`."""
        path, e = self._job("true\n", fran=FAKE_FRAN, qc=True, fran_name="HeLa standard",
                            FRAN_DEPOSIT="")
        text = read(path)
        self.assertIn("--no-fran", text)
        self.assertIn("--qc", text)
        self.assertNotIn("--fran-stage", text)
        r = subprocess.run(["bash", path], capture_output=True, text=True, env=e, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        # stage --qc stages nothing and records the decision (a later stage/backfill honours it)
        calls = self._calls("fran_deposit")
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0]["argv"][-3:], ["--name", "HeLa standard", "--qc"])
        self.assertNotIn("--skip", calls[0]["argv"])
        self.assertIn("*Staged for FRAN:* skipped: instrument QC / standard run",
                      json.dumps(self.m.bodies[-1]))

    def test_the_stage_argv_matches_fran_deposit_stage_argv(self):
        """The shape of fran_deposit.stage_argv() (the FRAN side's builder, on its branch):
        stage --out O [--fasta-meta M] [--name N] [--qc | --not-qc]; a blank name dropped."""
        py = sys.executable
        self.assertEqual(ns._stage_argv("/s/fran_deposit.py", "/o", name=" PROT_0793 mouse — x ",
                                        qc=True, fasta_meta="/m"),
                         [py, "/s/fran_deposit.py", "stage", "--out", "/o", "--fasta-meta", "/m",
                          "--name", "PROT_0793 mouse — x", "--qc"])
        self.assertEqual(ns._stage_argv("/f", "/o", name="  ", qc=False),
                         [py, "/f", "stage", "--out", "/o", "--not-qc"])
        self.assertEqual(ns._stage_argv("/f", "/o"), [py, "/f", "stage", "--out", "/o"])

    def test_the_hook_builder_and_fran_deposit_stage_argv_cannot_drift(self):
        """Two builders of one argv (architectural rule 3): the hook keeps its own so it never
        imports fran_deposit.py inside a job, and fran_deposit.stage_argv() is the FRAN side's.
        Both ship in this skill, so pin them to each other on every input shape."""
        import fran_deposit as fd
        path = os.path.abspath(fd.__file__)
        for kw in ({}, {"name": " PROT_0793 mouse — x "}, {"name": "  "}, {"qc": True},
                   {"qc": False}, {"fasta_meta": "/m"},
                   {"name": "Lumos QC", "qc": True, "fasta_meta": "/m"},
                   {"name": "HeLa study", "qc": False}):
            self.assertEqual(ns._stage_argv(path, "/o", **kw),
                             fd.stage_argv("/o", python=sys.executable, **kw), kw)

    def test_record_only_calls_have_the_same_time_bounds_as_stage(self):
        seen = []

        def fake_stage(out, session=None, timeout=None, **kw):
            seen.append((timeout, kw))
            return {"staged": False, "reason": "opted_out"}
        for qc, mode, skip in ((True, "off", False), (None, "off", True)):
            seen.clear()
            with mock.patch.object(ns, "_seconds_left", return_value=60), \
                    mock.patch.object(ns, "fran_stage", side_effect=fake_stage):
                r = ns._fran_step("/x", None, mode, "N", qc)
            self.assertEqual(seen, [])                           # under NEAR_LIMIT_S: not called
            self.assertIn("not recorded", r["detail"])
            with mock.patch.object(ns, "_seconds_left", return_value=200), \
                    mock.patch.object(ns, "fran_stage", side_effect=fake_stage):
                ns._fran_step("/x", None, mode, "N", qc)
            self.assertEqual(seen, [(170, {"name": "N", "qc": qc, "skip": skip})])

    def test_the_name_and_not_qc_reach_stage(self):
        for qc, tail in ((False, ["--name", "Mouse liver KO vs WT", "--not-qc"]),
                         (None, ["--name", "Mouse liver KO vs WT"])):
            os.path.exists(self.log) and os.remove(self.log)
            r = self._run("true\n", fran=FAKE_FRAN, qc=qc, fran_name="Mouse liver KO vs WT")
            self.assertEqual(r.returncode, 0, r.stderr)
            argv = self._calls("fran_deposit")[0]["argv"]
            self.assertEqual(argv[-len(tail):], tail, argv)
            self.assertNotIn("--qc", argv)

    def test_stage_timeout_is_bounded_by_the_time_left(self):
        seen = {}

        def fake_stage(out, session=None, timeout=None, **_):
            seen["timeout"] = timeout
            return {"staged": True, "reason": "ok"}
        for left, want in ((200, 170), (None, ns.FRAN_STAGE_TIMEOUT_S),
                           (5000, ns.FRAN_STAGE_TIMEOUT_S)):
            with mock.patch.object(ns, "_seconds_left", return_value=left), \
                    mock.patch.object(ns, "fran_stage", side_effect=fake_stage):
                ns._fran_step("/x", None, "stage")
            self.assertEqual(seen["timeout"], want, left)


class ScriptAndPayloadSafety(unittest.TestCase):
    def test_a_comment_between_sbatch_lines_stays_in_the_header(self):
        script = ("#!/bin/bash -l\n#SBATCH --job-name=t\n# the queue, chosen by slurm_queue()\n"
                  "#SBATCH --time=1:00:00\n\nset -euo pipefail\necho hi\n")
        with tempfile.TemporaryDirectory() as d:
            out = ns.wrap_job_script(script, d, final=True, time_limit_h=1, stage="search")
        lines = out.splitlines()
        hook = lines.index("_JOB_T0=$(date +%s)")
        self.assertLess(lines.index("#SBATCH --time=1:00:00"), hook)
        self.assertLess(lines.index("#SBATCH --job-name=t"), hook)
        self.assertGreater(lines.index("set -euo pipefail"), lines.index("("))

    def test_the_notification_text_is_escaped_too(self):
        with tempfile.TemporaryDirectory() as d:
            _, out = search_session(d)
            with env_patch(d, SLURM_JOB_ID="1"):
                f = ns.search_facts(out, exit_code=0)
            f["session"] = "<!channel> <https://evil.example|click here> & co"
            p = ns.payload(f)
        self.assertNotIn("<!channel>", json.dumps(p))
        self.assertNotIn("<https://evil", json.dumps(p))
        self.assertIn("&lt;!channel&gt;", p["text"])
        self.assertIn("&amp; co", p["text"])

    SECRETS = [
        ("-----BEGIN OPENSSH PRIVATE KEY-----\nb3BlbnNzaC1rZXk\n-----END OPENSSH PRIVATE KEY-----",
         "b3BlbnNzaC1rZXk"),
        ("ghp_" + "A1b2C3d4E5f6G7h8I9j0K1l2", "ghp_A1b2"),
        ("github_pat_11ABCDEFG0_abcdefXYZ", "github_pat_11"),
        ("hf_" + "abcdefghijklmnopqrstuvwxyz", "hf_abcdefgh"),
        ("Authorization: Bearer eyJhbGciOiJIUzI1NiJ9.x.y", "eyJhbGci"),
        ("password=hunter2", "hunter2"),
        ("postgresql://fran:s3cretpw@pgfarm.ucdavis.edu:5432/fran", "s3cretpw"),
        ("postgres://u:pw1234@h/db", "pw1234"),
        (FAKE_SLACK, SECRET_BIT),
    ]

    def test_send_alert_redacts_every_secret_pattern(self):
        with tempfile.TemporaryDirectory() as d:
            m = Mock("ok")
            try:
                with env_patch(d, SKILL_SLACK_WEBHOOK=m.url):
                    for secret, bit in self.SECRETS:
                        sent = ns.send_alert(f"ingest failed: {secret} (retrying)",
                                             title=f"FRAN {secret}")
                        self.assertTrue(sent)
                        posted = json.dumps(m.bodies[-1])
                        self.assertNotIn(bit, posted, secret)
                        self.assertIn("[redacted]", posted)
                        self.assertIn("(retrying)", posted)       # the alert still goes out
            finally:
                m.close()

    def test_job_end_error_tails_are_redacted_in_the_post_and_the_log(self):
        t = JobTrap("test_exit_codes_are_kept_and_reported")
        t.setUp()
        try:
            fran = FAKE_FRAN.replace(
                'raise RuntimeError("drop dir exploded")',
                'raise RuntimeError("could not reach postgresql://fran:s3cretpw@db/x with '
                'ghp_A1b2C3d4E5f6G7h8I9j0K1l2")')
            r = t._run("true\n", fran=fran, FAKE_FRAN_MODE="crash")
            self.assertEqual(r.returncode, 0, r.stderr)
            for text in (json.dumps(t.m.bodies[-1]), r.stderr):
                self.assertNotIn("s3cretpw", text)
                self.assertNotIn("ghp_A1b2", text)
                self.assertIn("[redacted]", text)
        finally:
            t.tearDown()

    def test_a_stalled_webhook_lookup_cannot_hold_send_alert(self):
        def stalled():
            time.sleep(20)                          # an open() on a hung quobyte mount
        with tempfile.TemporaryDirectory() as d, env_patch(d), \
                mock.patch.object(ns, "DELIVER_DEADLINE_S", 1), \
                mock.patch.object(ns, "resolve_webhook", side_effect=stalled):
            t = time.time()
            sent, status = ns.send_alert("x", with_status=True)
            took = time.time() - t
        self.assertFalse(sent)
        self.assertLess(took, 5)
        self.assertIn("no answer within", status)

    def test_an_alert_never_counts_skill_issues(self):
        with tempfile.TemporaryDirectory() as d:
            m = Mock("ok")
            try:
                with env_patch(d, SKILL_SLACK_WEBHOOK=m.url), \
                        mock.patch.object(ns, "add_issues") as ai:
                    self.assertTrue(ns.send_alert("x"))
            finally:
                m.close()
        ai.assert_not_called()

    def test_every_status_and_manifest_line_is_one_line(self):
        import session
        self.assertNotIn("\n", ns.slack_manifest(False, "not sent: a\nb\r\nc")[2])
        with mock.patch.object(ns, "deliver", return_value=(False, "not sent: x\ny\rz")), \
                tempfile.TemporaryDirectory() as d, env_patch(d, SLURM_JOB_ID="1"):
            _, out = search_session(d)
            _, status = ns.search_done(out, exit_code=1, from_job=True)
        self.assertNotIn("\n", status)
        self.assertNotIn("\r", status)
        with tempfile.TemporaryDirectory() as d:
            man = os.path.join(d, "MANIFEST.txt")
            session._append_manifest(man, "SKIPPED", "Slack notification", "bad\r\nsecond\nthird")
            self.assertEqual(len(read(man).splitlines()), 1)
            self.assertIn("bad second third", read(man))


class ZipManifestFallback(unittest.TestCase):
    """If MANIFEST.txt cannot go into the zip after the hooks, the pre-hook one does -- a zip with
    no manifest at all would hide every [SKIPPED] part (CLAUDE.md rule 4)."""

    def _finalize(self, d, fail_writes):
        import session
        import test_deposit_package as tdp
        p = tdp.dia_session(d)
        real = zipfile.ZipFile.writestr
        calls = {"n": 0}

        def flaky(zself, arcname, data, *a, **k):
            if str(arcname).endswith("MANIFEST.txt"):
                calls["n"] += 1
                if calls["n"] <= fail_writes:
                    raise OSError(28, "No space left on device")
            return real(zself, arcname, data, *a, **k)
        args = type("A", (), dict(dir=p["session_dir"], zip=True, no_deposit=True,
                                  no_notify=False, reanalysis_of=""))()
        out = io.StringIO()
        with mock.patch.dict(os.environ, job_env(d), clear=True), \
                mock.patch.object(zipfile.ZipFile, "writestr", flaky), \
                redirect_stdout(out), redirect_stderr(io.StringIO()):
            session.do_finalize(args)
        res = json.loads(out.getvalue())
        with zipfile.ZipFile(res["zip"]) as z:
            names = z.namelist()
            arc = f"{os.path.basename(p['session_dir'])}/MANIFEST.txt"
            zipped = z.read(arc).decode() if arc in names else None
        return res, zipped

    def test_a_failed_append_falls_back_to_the_pre_hook_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            res, zipped = self._finalize(d, fail_writes=1)
        self.assertTrue(res["zip_manifest"].startswith("added without the run-log/Slack lines"),
                        res["zip_manifest"])
        self.assertIsNotNone(zipped)
        self.assertIn("[SKIPPED] Deposit package", zipped)          # the export parts are there
        self.assertNotIn("Core run log", zipped)

    def test_both_attempts_failing_is_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            res, zipped = self._finalize(d, fail_writes=2)
        self.assertTrue(res["zip_manifest"].startswith("MISSING:"), res["zip_manifest"])
        self.assertIsNone(zipped)


if __name__ == "__main__":
    unittest.main()
