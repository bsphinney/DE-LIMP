#!/usr/bin/env python3
"""
One list of secret-shaped strings, and the places a secret leaked past it (release review, 2.8.0).

  * notify_slack._SECRET_PATTERNS is THE list. record_run.py imports it; report_issue.sh must run
    with no Python (Git Bash on Windows), so it carries an ERE mirror. Before this there were three
    lists that disagreed: report_issue.sh accepted Slack webhooks, xox tokens and Google API keys
    verbatim, and record_run.py had no share-token pattern at all.
  * ht_manifest.py printed the full URL -- `token=<share token>` included -- on any HTTP error but
    401/403, and the documented `--share-token <tok>` put the token in commands.log, which the run
    registry copied into a folder the whole Core group can read.

Every secret in this file is FAKE, built by concatenation so no scanner mistakes it for a real one.
Nothing here leaves the machine: the HTTP server is on 127.0.0.1.
"""
import http.server
import json
import os
import subprocess
import sys
import tempfile
import threading
import unittest
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import notify_slack as ns  # noqa: E402
import record_run as rr  # noqa: E402

FAKE_T = "Fk7" + "Ab3dEf6hIj9kLm2nOp5q"          # a fake STAN share token
FAKE_JWT = "eyJhbGciOiJIUzI1NiJ9" + "." + "eyJzdWIiOiJmYWtlIn0" + "." + "fakesignaturefake"
EXAMPLES = {
    "private key": "-----BEGIN " + "OPENSSH PRIVATE KEY-----\nAAAAfake\n-----END OPENSSH PRIVATE KEY-----",
    "github": "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1l2",
    "github oauth": "gho_" + "A1b2C3d4E5f6G7h8I9j0K1l2",           # and ghu_, ghs_, ghr_
    "github server": "ghs_" + "A1b2C3d4E5f6G7h8I9j0K1l2",
    "github pat": "github_pat_" + "11AAAAAAA0fakefakefake",
    "huggingface": "hf_" + "abcdefghijklmnopqrstuvwxyz",
    "authorization": "Authorization: Bearer " + "fakefakefake123",
    "password": "password=" + "hunter2",
    "postgres dsn": "postgresql://brett:" + "s3cretpw@pgfarm.example.edu/db",
    "slack webhook": "curl said: https://hooks.slack.com/services/" + "T0FAKE000/B0FAKE000/abcdefghijklmnopqrstuvwx",
    "google key": "key AIza" + "SyFAKEFAKEFAKEFAKEFAKEFAKEFAKEFAKE12",
    "google key AQ": "key AQ." + "Ab8RN6FAKEFAKEFAKEFAKEFAKE",
    "key=": "https://api.example.org/v1?key=" + "abcdefghijklmnopqrstuvwx",
    "slack token": "token xoxb-" + "1234567890-abcdefghij",
    "share token in a URL": "https://ucd.stan-proteomics.org/api/ht/manifest?q=0793&token=" + FAKE_T,
    "share token argv": "ht_manifest.py fetch 0793 --share-token " + FAKE_T + " --out ~/ht0793",
    "share token argv, quoted": "--share-token '" + FAKE_T + "'",
    "share token argv, =": "--share-token=" + FAKE_T,
    "share token env": "STAN_HT_SHARE_TOKEN=" + FAKE_T + " python3 ht_manifest.py fetch 0793",
    "cookie argv": "--cookie 'AppServiceAuthSession=" + FAKE_T + "'",
    # release verification, 2.8.0: report_issue.sh recorded all of these
    "anthropic key": "export KEY=sk-ant-" + "api03-FAKEfakeFAKEfake0123456789abcdef",
    "openai project key": "sk-proj-" + "FAKEfakeFAKEfake0123456789_-ab",
    "openai key": "sk-" + "FAKEfakeFAKEfake0123456789abcdefFAKEfake01234567",
    "aws access key id": "aws_access_key_id = AKIA" + "IOSFODNN7EXAMPLE",
    "bare bearer": "curl -H 'Bearer " + "FAKEfakeFAKEfake0123456789'",
    "jwt": "the PG Farm token is " + FAKE_JWT,
    "share token argv, two spaces": "--share-token  " + FAKE_T,
    "share token argv, tab": "--share-token\t" + FAKE_T,
    "cookie argv, tab": "--cookie\tAppServiceAuthSession=" + FAKE_T,
}
# Look-alikes that must pass: file NAMES (the documented forms), placeholders, prose.
BENIGN = [
    "--share-token-file ~/.stan/share_0793",
    "--share-token-file=/home/u/.stan/share_0793",
    "STAN_PG_TOKEN=/quobyte/proteomics-grp/etc/pgfarm_token",
    "--token /quobyte/proteomics-grp/etc/pgfarm_token",
    "--share-token <tok>",
    "--cookie-file ~/.stan/cookie",
    "the query takes token=~/.stan/share_0793",
    "STAN_HT_SHARE_TOKEN=$HOME/.stan/share_0793 python3 ht_manifest.py",
    "a token= field with no value",
    "key=short",
    "set the webhook up at hooks.slack.com",
    "password",
    "AQ.short",
    "detect_acquisition.py returned 'not found' for all 15 .raw",
    "sk-learn",
    "risk-2026-09-28_session_name_with_many_parts_xyz",
    "sk-2026-09-28_session_name_foo_bar_baz_qux",
    "AKIA is the AWS access-key prefix",
    "Bearer tokens are sent in a header",
    'eyJ is base64 for {"',
    "--share-token\t<tok>",
]


class OneListTests(unittest.TestCase):
    def test_every_canonical_pattern_has_an_example_here(self):
        """A pattern added to notify_slack without a fake example here fails -- and with one,
        the mirror test below fails until report_issue.sh catches it too."""
        for pat in ns._SECRET_PATTERNS:
            self.assertTrue(any(pat.search(t) for t in EXAMPLES.values()),
                            f"no example in EXAMPLES exercises {pat.pattern!r}")

    def test_the_canonical_list_catches_every_example_and_spares_the_look_alikes(self):
        for name, text in EXAMPLES.items():
            with self.subTest(name=name):
                self.assertTrue(ns.contains_secret(text))
                self.assertIn(ns.REDACTED, ns.redact(text))
        self.assertNotIn(FAKE_T, ns.redact(EXAMPLES["share token in a URL"]))
        for text in BENIGN:
            with self.subTest(benign=text):
                self.assertFalse(ns.contains_secret(text))
                self.assertEqual(ns.redact(text), text)

    def test_record_run_uses_that_list_and_keeps_none_of_its_own(self):
        self.assertFalse(hasattr(rr, "SECRET_TEXT_RE"))
        self.assertIs(rr.contains_secret, ns.contains_secret)


class ReportIssueMirrorTests(unittest.TestCase):
    """report_issue.sh's ERE mirror must refuse exactly what the canonical list catches."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.shared = os.path.join(self.d, "shared")
        os.makedirs(self.shared)

    def tearDown(self):
        self._tmp.cleanup()

    def report(self, text, *more):
        env = {k: v for k, v in os.environ.items() if not k.startswith("HIVE_")}
        env.update(HOME=self.d, HIVE_ENV_FILE=os.path.join(self.d, "none"), TMPDIR=self.d,
                   SKILL_ISSUES_DIR=self.shared,
                   SKILL_ISSUES_LOCAL_DIR=os.path.join(self.d, "local"))
        return subprocess.run(["bash", os.path.join(SCRIPTS, "report_issue.sh"), "--title", "t",
                               "--what", text, *more], capture_output=True, text=True, env=env,
                              timeout=60)

    def test_every_example_is_refused_and_nothing_is_written(self):
        for name, text in EXAMPLES.items():
            with self.subTest(name=name):
                p = self.report(text)
                self.assertEqual(p.returncode, 2, p.stdout + p.stderr)
                self.assertNotIn(FAKE_T, p.stdout + p.stderr)
        self.assertEqual(os.listdir(self.shared), [])

    def test_mode_is_checked_too(self):
        """--mode is written into the file, and was the one free-text field left unchecked."""
        p = self.report("t", "--mode", "hive_remote " + EXAMPLES["openai key"])
        self.assertEqual(p.returncode, 2, p.stdout + p.stderr)
        self.assertEqual(os.listdir(self.shared), [])

    def test_the_look_alikes_are_recorded(self):
        for text in BENIGN:
            with self.subTest(benign=text):
                p = self.report(text)
                self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(len(os.listdir(self.shared)), len(BENIGN))


class RegistryTests(unittest.TestCase):
    """record_run.py refuses (and lists) a file whose contents match: a commands.log with the
    token on a command line never reaches the group-readable registry."""

    def check(self, text):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "commands.log")
            with open(p, "w") as fh:
                fh.write("bash scripts/hive_exec.sh 'python3 ht_manifest.py fetch 0793 "
                         f"--http https://ucd.stan-proteomics.org {text} --out ~/ht0793'\n")
            return rr.secret_reason(p)

    def test_a_logged_share_token_is_refused(self):
        self.assertIn("token", self.check(f"--share-token {FAKE_T}"))

    def test_the_documented_file_form_is_copied(self):
        self.assertIsNone(self.check("--share-token-file ~/.stan/share_0793"))

    def test_without_the_list_nothing_is_copied_unchecked(self):
        saved = rr.contains_secret
        rr.contains_secret = None
        try:
            self.assertIn("cannot be checked", self.check("--share-token-file ~/.stan/x"))
        finally:
            rr.contains_secret = saved


class _Echo(http.server.BaseHTTPRequestHandler):
    """A STAN stand-in that fails, and echoes the request line in its error body -- what a
    server's error page may well do."""
    status = 500
    seen = []

    def do_GET(self):  # noqa: N802
        type(self).seen.append({"path": self.path, "cookie": self.headers.get("Cookie")})
        if type(self).status == 200:
            # a reply that echoes what it was sent into every field fetch() prints or keeps
            tok = (urllib.parse.parse_qs(urllib.parse.urlparse(self.path).query).get("token")
                   or [""])[0]
            body = json.dumps({"submission": "0793", "include": "samples", "files": [],
                               "plates": [f"P1-{tok}"], "counts": {"sample": 1, f"x{tok}": 2},
                               "missing_paths": [f"/raw/{tok}.d"], "note": self.path}).encode()
        else:
            body = f"internal error handling GET {self.path}".encode()
        self.send_response(type(self).status)
        self.send_header("Content-Type", "text/plain")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *a):
        pass


class HtManifestTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        _Echo.seen = []
        self.srv = http.server.HTTPServer(("127.0.0.1", 0), _Echo)
        threading.Thread(target=self.srv.serve_forever, daemon=True).start()
        self.base = f"http://127.0.0.1:{self.srv.server_port}"
        self.tokfile = os.path.join(self.d, "share_0793")
        with open(self.tokfile, "w") as fh:
            fh.write(FAKE_T + "\n")
        os.chmod(self.tokfile, 0o600)

    def tearDown(self):
        self.srv.shutdown()
        self.srv.server_close()
        self._tmp.cleanup()

    def fetch(self, *args, status=500):
        _Echo.status = status
        env = {k: v for k, v in os.environ.items() if not k.startswith("STAN_")}
        return subprocess.run([sys.executable, os.path.join(SCRIPTS, "ht_manifest.py"), "fetch",
                               "0793", "--http", self.base, "--out", self.d, *args],
                              capture_output=True, text=True, env=env, timeout=60)

    def test_an_http_error_never_prints_the_token(self):
        """The reviewer's repro: any code but 401/403 printed the whole URL, token included --
        and here the server echoes it back in the body as well."""
        p = self.fetch("--share-token-file", self.tokfile)
        self.assertEqual(p.returncode, 3, p.stderr)
        self.assertIn("HTTP 500", p.stderr)
        self.assertNotIn(FAKE_T, p.stdout + p.stderr)
        self.assertIn("[redacted]", p.stderr)
        # the token from the file WAS sent: only the printing is masked
        sent = urllib.parse.parse_qs(urllib.parse.urlparse(_Echo.seen[0]["path"]).query)
        self.assertEqual(sent["token"], [FAKE_T])

    def test_a_refusal_points_at_the_file_forms(self):
        p = self.fetch("--share-token-file", self.tokfile, status=403)
        self.assertEqual(p.returncode, 3)
        self.assertIn("--share-token-file", p.stderr)
        self.assertNotIn(FAKE_T, p.stdout + p.stderr)

    def test_the_argv_form_still_works_but_says_why_not(self):
        p = self.fetch("--share-token", FAKE_T)
        self.assertIn("use --share-token-file", p.stderr)
        self.assertNotIn(FAKE_T, p.stdout + p.stderr)

    def test_a_token_file_others_can_read_is_flagged(self):
        os.chmod(self.tokfile, 0o644)
        self.assertIn("chmod 600", self.fetch("--share-token-file", self.tokfile).stderr)

    def test_the_cookie_file_is_sent_and_never_printed(self):
        cookie = os.path.join(self.d, "cookie")
        with open(cookie, "w") as fh:
            fh.write("AppServiceAuthSession=" + FAKE_T)
        os.chmod(cookie, 0o600)
        p = self.fetch("--cookie-file", cookie)
        self.assertEqual(_Echo.seen[0]["cookie"], "AppServiceAuthSession=" + FAKE_T)
        self.assertNotIn(FAKE_T, p.stdout + p.stderr)

    def test_a_200_reply_that_echoes_the_token_is_masked_in_print_and_in_the_json(self):
        """fetch() printed plates, counts and the missing-path examples straight from the reply,
        and wrote ht_manifest.json verbatim (release verification, 2.8.0)."""
        p = self.fetch("--share-token-file", self.tokfile, status=200)
        self.assertEqual(p.returncode, 2, p.stdout + p.stderr)       # the missing path gate
        self.assertIn("plates", p.stdout)
        self.assertNotIn(FAKE_T, p.stdout + p.stderr)
        self.assertIn("[redacted]", p.stdout)
        with open(os.path.join(self.d, "ht_manifest.json")) as fh:
            text = fh.read()
        self.assertNotIn(FAKE_T, text)
        self.assertIn("[redacted]", text)
        self.assertIn("gates", json.loads(text))

    def test_the_cli_path_masks_the_pg_credential(self):
        """The STAN CLI route reads the PG Farm credential; a STAN that echoes it must not get it
        printed or kept. Not JWT-shaped on purpose: the value itself is masked, not a pattern."""
        fake_pg = "pgFAKE" + "Zz9yXw8vUt7sRq6pOn5m"
        pg = os.path.join(self.d, "pg_token")
        with open(pg, "w") as fh:
            fh.write(fake_pg + "\n")
        stan = os.path.join(self.d, "fake_stan")
        with open(stan, "w") as fh:
            fh.write("#!/bin/sh\nprintf '{\"files\": [], \"plates\": [\"%s\"], "
                     "\"missing_paths\": [\"/raw/%s.d\"]}' \"$PGPASSWORD\" \"$PGPASSWORD\"\n")
        os.chmod(stan, 0o755)
        env = {k: v for k, v in os.environ.items() if not k.startswith(("STAN_", "PG"))}
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "ht_manifest.py"), "--stan", stan,
                            "--token", pg, "fetch", "0793", "--out", self.d],
                           capture_output=True, text=True, env=env, timeout=60)
        self.assertEqual(p.returncode, 2, p.stdout + p.stderr)
        self.assertNotIn(fake_pg, p.stdout + p.stderr)
        with open(os.path.join(self.d, "ht_manifest.json")) as fh:
            self.assertNotIn(fake_pg, fh.read())

    def test_nothing_it_writes_holds_the_token(self):
        self.fetch("--share-token-file", self.tokfile)
        for name in os.listdir(self.d):
            path = os.path.join(self.d, name)
            if os.path.isfile(path) and path != self.tokfile:
                with open(path, errors="replace") as fh:
                    self.assertNotIn(FAKE_T, fh.read(), name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
