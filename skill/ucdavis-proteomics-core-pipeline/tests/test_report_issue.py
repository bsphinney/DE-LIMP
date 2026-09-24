#!/usr/bin/env python3
"""
report_issue.sh is how a session's skill defects reach the people who fix them. The ten
defects of 2026-09-23 (a Windows Core member in hive_remote mode) were only fixed because the
agent happened to leave a notes file in the user's folder; this script makes that routine.

What must hold, because a failure here is silent -- nobody is waiting for the report:
  * three routes -- write to the Core's shared folder directly (on HIVE), through
    hive_exec.sh over SSH (a laptop in hive_remote mode), or into a local folder with a note
    to send it on -- and a failed shared write falls back to local instead of losing it;
  * one file per user/day/session, the header written once, entries appended in order;
  * an entry without optional fields is still written (a `[ -n "$FIX" ] && printf` as the
    last line of the block once made every such entry look like a failed write);
  * text that looks like a key, token or password is refused, and nothing is written;
  * the session name cannot steer the file out of the folder.

Nothing here contacts HIVE: HIVE_EXEC is replaced by a stub that runs the remote command
against a temp directory standing in for /quobyte.
"""
import os
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "report_issue.sh")
FAKE_ROOT = "/nonexistent_hive_root_for_tests"     # never exists locally -> not "direct"


class ReportIssue(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.shared = os.path.join(self.d, "shared")
        self.local = os.path.join(self.d, "local")
        self.remote = os.path.join(self.d, "remote")          # stands in for FAKE_ROOT on "HIVE"
        self.key = os.path.join(self.d, "id_test")
        open(self.key, "w").close()
        # hive_exec.sh stand-in: run the command locally with FAKE_ROOT mapped to self.remote,
        # stdin passed through exactly as ssh would.
        self.fake_exec = os.path.join(self.d, "fake_hive_exec.sh")
        with open(self.fake_exec, "w") as fh:
            fh.write('#!/usr/bin/env bash\ncmd="${1//%s/%s}"\nexec bash -c "$cmd"\n'
                     % (FAKE_ROOT, self.remote))
        os.chmod(self.fake_exec, 0o755)

    def tearDown(self):
        self._tmp.cleanup()

    def run_it(self, *args, route="direct", extra_env=None):
        env = {k: v for k, v in os.environ.items() if not k.startswith("HIVE_")}
        env.update(HOME=self.d, HIVE_ENV_FILE=os.path.join(self.d, "no-hive.env"),
                   SKILL_ISSUES_LOCAL_DIR=self.local, TMPDIR=self.d)
        if route == "direct":
            env["SKILL_ISSUES_DIR"] = self.shared
        else:
            env["SKILL_ISSUES_DIR"] = FAKE_ROOT + "/skill_issues"
        if route == "ssh":
            env.update(HIVE_USER="gabrig", HIVE_KEY=self.key, HIVE_EXEC=self.fake_exec)
        env.update(extra_env or {})
        return subprocess.run(["bash", SCRIPT, *args], capture_output=True, text=True, env=env)

    def files(self, d):
        return sorted(os.listdir(d)) if os.path.isdir(d) else []

    def read_only_file(self, d):
        names = [n for n in self.files(d) if n.endswith(".md")]
        self.assertEqual(len(names), 1, names)
        with open(os.path.join(d, names[0])) as fh:
            return names[0], fh.read()

    # ------------------------------------------------------------------ routes
    def test_direct_header_once_entries_in_order(self):
        os.makedirs(self.shared)
        for title in ("first", "second"):
            p = self.run_it("--title", title, "--what", "w", "--session", "chkLUppm")
            self.assertEqual(p.returncode, 0, p.stderr)
        name, text = self.read_only_file(self.shared)
        self.assertTrue(name.endswith("_chkLUppm.md"), name)
        self.assertEqual(text.count("# Skill issues"), 1, text)
        self.assertLess(text.index("## first"), text.index("## second"))
        self.assertIn("**Skill version:**", text)

    def test_entry_without_optional_fields_is_written(self):
        os.makedirs(self.shared)
        p = self.run_it("--title", "bare", "--what", "only the required fields")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertNotIn("cannot write", p.stderr)
        _, text = self.read_only_file(self.shared)
        self.assertIn("## bare", text)
        self.assertNotIn("Proposed fix", text)

    def test_ssh_route_writes_on_hive(self):
        os.makedirs(os.path.join(self.remote, "skill_issues"))
        p = self.run_it("--title", "over ssh", "--what", "w", "--fix", "f", route="ssh")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("recorded: gabrig@hive:", p.stdout)
        name, text = self.read_only_file(os.path.join(self.remote, "skill_issues"))
        self.assertTrue(name.startswith(tuple("0123456789")) and "_gabrig" in name, name)
        self.assertIn("hive_remote", text)
        self.assertEqual(self.files(self.local), [])

    def test_ssh_failure_falls_back_to_local(self):
        # no skill_issues folder on "HIVE" (a non-Core account): kept locally, not lost
        p = self.run_it("--title", "kept", "--what", "w", route="ssh")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("recorded locally", p.stdout)
        _, text = self.read_only_file(self.local)
        self.assertIn("## kept", text)

    def test_no_hive_at_all_goes_local_with_instructions(self):
        p = self.run_it("--title", "laptop", "--what", "w", route="local")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("send that file", p.stdout)
        self.read_only_file(self.local)

    def test_where_names_the_route(self):
        os.makedirs(self.shared)
        self.assertTrue(self.run_it("--where").stdout.startswith("direct:"))
        self.assertTrue(self.run_it("--where", route="ssh").stdout.startswith("ssh:"))
        self.assertTrue(self.run_it("--where", route="local").stdout.startswith("local:"))

    # ------------------------------------------------------------------ guards
    def test_secrets_are_refused_and_nothing_written(self):
        os.makedirs(self.shared)
        for bad in ("-----BEGIN OPENSSH PRIVATE KEY-----", "ghp_" + "a" * 30,
                    "Authorization: Token 0123456789abcdef", "password=hunter2"):
            p = self.run_it("--title", "t", "--what", bad)
            self.assertEqual(p.returncode, 2, bad)
            self.assertIn("key, token or password", p.stderr)
        self.assertEqual(self.files(self.shared), [])

    def test_session_name_cannot_escape_the_folder(self):
        os.makedirs(self.shared)
        p = self.run_it("--title", "t", "--what", "w", "--session", "../../etc/x y")
        self.assertEqual(p.returncode, 0, p.stderr)
        name, _ = self.read_only_file(self.shared)
        self.assertNotIn("/", name)
        self.assertNotIn(" ", name)

    def test_required_fields_and_enums(self):
        self.assertEqual(self.run_it("--title", "t").returncode, 2)
        self.assertEqual(self.run_it("--title", "t", "--what", "w", "--kind", "oops").returncode, 2)
        self.assertEqual(self.run_it("--title", "t", "--what", "w",
                                     "--severity", "urgent").returncode, 2)


if __name__ == "__main__":
    unittest.main()
