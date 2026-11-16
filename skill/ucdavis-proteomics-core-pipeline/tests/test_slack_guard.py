#!/usr/bin/env python3
"""No test reaches the Core's real Slack (review A10, Ruling 109.1).

On HIVE, test_slack_notify.AlertApi.test_sent_off_and_failed's "no webhook" call posted
":warning: FRAN health: FRAN ingest stalled: ..." to the Core channel on every suite run:
notify_slack reads its group webhook path (GROUP_FILE) from the environment at IMPORT, the test
patched the environment later, and the Core's real group file was read.

These tests rebuild that situation in a fresh python whose environment names a group webhook file
from the START (as HIVE's import-time default does) -- a fake one, holding a loopback URL to a
server the test listens on -- and check that:
  * test_slack_notify's own env_patch() moves the import-time paths: its no-webhook call is off and
    the server hears nothing (fails on the tree before the fix: the server got the alert);
  * tests/slack_guard.py (installed by job_env.py for every test process) stops an in-process
    send that does NOT go through env_patch, and records it;
  * the guard refuses a webhook file outside a temp dir, a real Slack URL, and a relay through a
    hive_exec outside a temp dir.
Stdlib only; nothing leaves the machine (loopback only).
"""
import http.server
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import unittest
import urllib.request
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

try:
    import slack_guard  # noqa: E402
except ImportError:         # the tree before the fix: the tests that need it fail on their own
    slack_guard = None

PY = sys.executable


class Ear:
    """A loopback HTTP server that counts what is posted to it."""

    def __init__(self):
        bodies = self.bodies = []

        class H(http.server.BaseHTTPRequestHandler):
            def do_POST(self):
                n = int(self.headers.get("Content-Length") or 0)
                bodies.append(self.rfile.read(n).decode("utf-8", "replace"))
                self.send_response(200)
                self.end_headers()
                self.wfile.write(b"ok")

            def log_message(self, *a):
                pass

        self.srv = http.server.HTTPServer(("127.0.0.1", 0), H)
        self.url = f"http://127.0.0.1:{self.srv.server_address[1]}/hook"
        threading.Thread(target=self.srv.serve_forever, daemon=True).start()

    def close(self):
        self.srv.shutdown()
        self.srv.server_close()


# Run in a fresh python: SKILL_SLACK_GROUP_FILE is set from the start, as HIVE's real default is.
# argv[1] = "env_patch" (test_slack_notify's helper, the A10 call) or "bare" (an in-process call
# with no helper at all: only slack_guard, through job_env, stands between it and the file).
CHILD = r"""
import contextlib, io, json, os, sys, tempfile
sys.path[:0] = [sys.argv[2], sys.argv[3]]
mode = sys.argv[1]
try:
    import slack_guard
except ImportError:          # the tree before the fix
    slack_guard = None
if mode == "env_patch":
    import test_slack_notify as t
    t.setUpModule()
    if slack_guard:
        # env_patch() ALONE must move the import-time paths: the process-wide guard is lifted,
        # so GROUP_FILE is back to the import-time value, the probe's group file
        slack_guard.uninstall()
    try:
        with tempfile.TemporaryDirectory() as d, t.env_patch(d), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            sent, status = t.ns.send_alert("FRAN ingest stalled (slack_guard proof)",
                                           title="FRAN health", with_status=True)
    finally:
        t._MODULE_TMP and __import__("shutil").rmtree(t._MODULE_TMP, ignore_errors=True)
    print(json.dumps({"sent": sent, "status": status}))
else:
    import job_env, notify_slack as ns
    with tempfile.TemporaryDirectory() as d:
        os.environ.pop("SKILL_SLACK", None)
        os.environ.update(HOME=d, SKILL_SLACK_TEST_LOOPBACK="1",
                          HIVE_ENV_FILE=os.path.join(d, "none"))
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            sent, status = ns.send_alert("bare call (slack_guard proof)", with_status=True)
    print(json.dumps({"sent": sent, "status": status, "group_file": ns.GROUP_FILE,
                      "violations": slack_guard.violations() if slack_guard else None}))
"""


class ImportTimeGroupFile(unittest.TestCase):
    """A group webhook file named in the environment before notify_slack is imported."""

    def setUp(self):
        # NOT under a temp dir where it can be helped: the guard's file check must see it as the
        # real group file. Beside the skill tree when that is not itself in a temp dir (the Mac),
        # else (HIVE's suite unpacks the tree under /tmp) in a temp dir -- the env_patch and
        # GROUP_FILE checks below do not depend on where it is.
        base = os.path.dirname(os.path.dirname(HERE))
        tmp = os.path.realpath(tempfile.gettempdir())
        root = (tempfile.gettempdir() if os.path.realpath(base).startswith(tmp + os.sep)
                or os.path.realpath(base).startswith(("/tmp/", "/private/tmp/")) else base)
        self.d = tempfile.mkdtemp(prefix=".slack_guard_probe_", dir=root)
        self.addCleanup(shutil.rmtree, self.d, True)
        self.ear = Ear()
        self.addCleanup(self.ear.close)
        self.group = os.path.join(self.d, "skill_slack_webhook")
        with open(self.group, "w") as fh:
            fh.write(self.ear.url + "\n")

    def child(self, mode):
        env = {k: v for k, v in os.environ.items()
               if not k.startswith(("SKILL_SLACK", "HIVE_", "SLURM_"))}
        env.update(SKILL_SLACK_GROUP_FILE=self.group, SKILL_SLACK_TEST_LOOPBACK="1",
                   SKILL_CORE_GROUP_DIR=self.d, PYTHONDONTWRITEBYTECODE="1")
        r = subprocess.run([PY, "-c", CHILD, mode, SCRIPTS, HERE], capture_output=True,
                           text=True, env=env, timeout=120)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return json.loads(r.stdout.strip().splitlines()[-1])

    def test_env_patch_moves_the_import_time_paths(self):
        """the A10 call: test_slack_notify's no-webhook send_alert, as AlertApi makes it"""
        out = self.child("env_patch")
        self.assertEqual(self.ear.bodies, [], "the test alert reached the group webhook")
        self.assertEqual((out["sent"], out["status"]),
                         (False, "off: Core notification not configured for this user"))

    def test_the_guard_stops_a_bare_in_process_call(self):
        """no env_patch at all: job_env's slack_guard has moved GROUP_FILE into its temp folder"""
        out = self.child("bare")
        self.assertEqual(self.ear.bodies, [], "a bare in-process call reached the group file")
        self.assertIs(out["sent"], False)
        self.assertNotEqual(out["group_file"], self.group)


class TheGuardRefuses(unittest.TestCase):
    """slack_guard's three stops, each recorded (check() raises on them)."""

    def setUp(self):
        self.assertIsNotNone(slack_guard, "tests/slack_guard.py is missing")
        slack_guard.install()
        slack_guard.clear()
        self.addCleanup(slack_guard.clear)
        import notify_slack as ns
        self.ns = ns

    def test_a_webhook_file_outside_a_temp_dir(self):
        with tempfile.TemporaryDirectory() as d:
            group = os.path.join(d, "skill_slack_webhook")
            with open(group, "w") as fh:
                fh.write("https://hooks.slack.com/services/T0FAKE000/B0FAKE000/x\n")
            env = {"HOME": d, "SKILL_SLACK_TEST_LOOPBACK": "0"}
            # the file as if it were not in a temp dir (the Core's real one is on /quobyte)
            with mock.patch.dict(os.environ, env), \
                    mock.patch.object(self.ns, "GROUP_FILE", group), \
                    mock.patch.object(slack_guard, "_temp_roots", lambda: {"/nowhere-temp"}):
                os.environ.pop("SKILL_SLACK_WEBHOOK", None)
                url, why = self.ns.resolve_webhook()
        self.assertIsNone(url)
        self.assertIn("slack_guard: refused a webhook read from", why)
        with self.assertRaises(AssertionError):
            slack_guard.check()

    def test_a_real_slack_url(self):
        with self.assertRaises(OSError):
            urllib.request.urlopen("https://hooks.slack.com/services/TREAL/BREAL/realtoken",
                                   timeout=1)
        self.assertTrue(slack_guard.violations())

    def test_a_relay_through_the_real_hive(self):
        ok, why = self.ns._relay({"kind": "test"}, "/opt/real/hive_exec.sh")
        self.assertFalse(ok)
        self.assertIn("slack_guard: refused a relay", why)
        self.assertTrue(slack_guard.violations())

    def test_job_env_installs_it(self):
        import job_env  # noqa: F401
        self.assertTrue(slack_guard._STATE["installed"])
        self.assertTrue(slack_guard.inside_temp(self.ns.GROUP_FILE), self.ns.GROUP_FILE)


if __name__ == "__main__":
    unittest.main()
