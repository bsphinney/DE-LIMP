#!/usr/bin/env python3
"""
A script piped to `python3 -` has no folder of its own: Python sets __file__ to "<stdin>" (3.9,
3.10 on HIVE, 3.13 -- checked) and puts the CURRENT directory first on sys.path. notify_slack.py
(`python3 - relay`, `python3 - --test`) and record_run.py (`python3 - list`) go to HIVE that
way, where the current directory is the home folder. Treating "__file__ exists" as "running
from a file" made HERE the home folder: sibling files read from there, and a stray
~/core_submission.py -- or ~/statistics.py, shadowing the standard library -- imported.

Now HERE comes from __file__ only when it is a real file; otherwise there is no HERE and the
current directory comes off sys.path before any other import. Each test pipes the script to
`python3 -` from a temporary folder of decoys -- modules that leave a marker and raise when
imported, and a hive_exec.sh that leaves a marker when run -- and checks that none is touched
and the command still does its job.
"""
import base64
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
from job_env import job_env     # noqa: E402

DECOY_MODULE = ('import os\n'
                'open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "IMPORTED_{name}"),'
                ' "w").close()\n'
                'raise RuntimeError("decoy {name} was imported from the current directory")\n')
DECOY_HIVE_EXEC = '#!/bin/bash\ntouch "$(dirname "$0")/RAN_hive_exec"\necho SKILLDECOY\n'
# siblings both scripts import, and a standard-library name both import at the top
MODULES = ("core_submission", "notify_slack", "record_run", "skill_version", "fetch_fasta",
           "make_methods", "session", "check_report_runs", "fran_deposit", "statistics")
# Piped to `python3 -` itself: runs each script's source as __main__ with __file__ "<stdin>"
# (as a piped run has it; `--help` stops main), then asks the functions that find things
# beside the script what they found. The commands' own output does not show all of them: the
# --test payload has no version, and `list` names hive_exec.sh only on the SSH route.
PROBE = r"""
import contextlib, io, json, os, sys
out = {}
path0 = list(sys.path)
for name in ("notify_slack.py", "record_run.py"):
    with open(os.path.join(os.environ["PROBE_SCRIPTS"], name), encoding="utf-8") as fh:
        src = fh.read()
    sys.path[:] = path0                   # each starts with the current directory first
    g = {"__name__": "__main__", "__file__": "<stdin>"}
    sys.argv = ["-", "--help"]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            exec(compile(src, "<stdin>", "exec"), g)
        except SystemExit:
            pass
    if name == "notify_slack.py":
        out.update(ns_HERE=g["HERE"], relay_route=g["_relay_route"](),
                   skill_version=g["_skill_version"](), helper=g["_helper"]("record_run.py"),
                   unknown=g["SKILL_VERSION_UNKNOWN"])
    else:
        out.update(rr_HERE=g["HERE"], hive_exec_path=g["hive_exec_path"](),
                   rr_skill_version=g["skill_version"]())
    out[name + " cwd on sys.path"] = any(p in ("", ".", os.getcwd()) for p in sys.path)
print(json.dumps(out))
"""


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


class Harness(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._tmp.name)
        self.home = os.path.join(self.d, "home")            # the HIVE home: the current dir
        os.makedirs(self.home)
        for name in MODULES:
            with open(os.path.join(self.home, f"{name}.py"), "w") as fh:
                fh.write(DECOY_MODULE.format(name=name))
        with open(os.path.join(self.home, "hive_exec.sh"), "w") as fh:
            fh.write(DECOY_HIVE_EXEC)
        os.chmod(os.path.join(self.home, "hive_exec.sh"), 0o755)
        # a plugin.json where a HERE of "the current directory" would look for the version
        os.makedirs(os.path.join(self.d, ".claude-plugin"))
        with open(os.path.join(self.d, ".claude-plugin", "plugin.json"), "w") as fh:
            json.dump({"version": "6.6.6-decoy"}, fh)
        self.key = os.path.join(self.d, "id_test")
        open(self.key, "w").close()

    def tearDown(self):
        self._tmp.cleanup()

    def piped(self, script, *args, **env):
        """`python3 - <args>` with `script` on stdin, run from the decoy home."""
        with open(os.path.join(SCRIPTS, script), "rb") as fh:
            src = fh.read()
        e = job_env(self.d, HOME=self.home, **{k: v for k, v in env.items() if v is not None})
        for k in [k for k, v in env.items() if v is None]:
            e.pop(k, None)
        return subprocess.run([sys.executable, "-", *args], input=src, cwd=self.home, env=e,
                              capture_output=True, timeout=120)

    def assert_untouched(self, r):
        touched = sorted(n for n in os.listdir(self.home)
                         if n.startswith(("IMPORTED_", "RAN_")))
        self.assertEqual(touched, [], r.stderr.decode("utf-8", "replace"))
        self.assertNotIn(b"decoy", r.stdout + r.stderr)


class RecordRunList(Harness):
    def registry(self):
        reg = os.path.join(self.d, "skill_runs")
        rec = os.path.join(reg, "sessions", "2026-09-29_Test_Run")
        os.makedirs(rec)
        with open(os.path.join(rec, "run_record.json"), "w") as fh:
            json.dump({"schema_version": 2, "name": "2026-09-29_Test_Run", "status": "completed",
                       "user": "tester", "engine": "diann", "submitted": "2026-09-29T10:00"}, fh)
        return reg

    def test_the_hive_side_of_list_imports_nothing_from_the_current_directory(self):
        """What the laptop runs on HIVE: `python3 - list --json --remote-hop`."""
        reg = self.registry()
        r = self.piped("record_run.py", "list", "--json", "--remote-hop",
                       RECORD_RUN="on", SKILL_RUNS_DIR=reg)
        self.assertEqual(r.returncode, 0, r.stderr.decode())
        self.assert_untouched(r)
        rows = json.loads(r.stdout.decode()[r.stdout.decode().index("["):])
        self.assertEqual([row.get("name") or row.get("session") for row in rows],
                         ["2026-09-29_Test_Run"])

    def test_no_hive_exec_is_found_beside_nothing(self):
        """No HERE: no hive_exec.sh beside it, so a HIVE login does not turn into an SSH route
        through ~/hive_exec.sh, and open("<stdin>") is never tried."""
        r = self.piped("record_run.py", "list", RECORD_RUN="on",
                       SKILL_RUNS_DIR=os.path.join(self.d, "not_here", "skill_runs"),
                       HIVE_USER="tester", HIVE_KEY=self.key)
        self.assertEqual(r.returncode, 0, r.stderr.decode())
        self.assert_untouched(r)
        out = r.stdout + r.stderr
        self.assertNotIn(b"via hive_exec.sh", out)            # the SSH route's where_text
        self.assertNotIn(b"could not list over SSH", out)
        self.assertIn(b"cannot list -- not on HIVE", out)


class NotifySlack(Harness):
    def test_a_relayed_post_imports_nothing_from_the_current_directory(self):
        """The HIVE side of a relay: `python3 - relay --facts-b64 ...`, to a loopback webhook."""
        sys.path.insert(0, HERE)
        from test_slack_notify import Mock
        m = Mock("ok")
        try:
            facts = {"kind": "test", "status": "test", "who": "tester", "host": "laptop"}
            b64 = base64.b64encode(json.dumps(facts).encode()).decode()
            env = {k: v for k, v in dict(
                SKILL_SLACK_WEBHOOK=m.url, SKILL_SLACK_TEST_LOOPBACK="1",
                SKILL_SLACK_GROUP_FILE=os.path.join(self.d, "no_group_webhook"),
                SKILL_CORE_GROUP_DIR=os.path.join(self.d, "no_core_group")).items()}
            e = job_env(self.d, HOME=self.home, **env)
            e.pop("SKILL_SLACK")
            with open(os.path.join(SCRIPTS, "notify_slack.py"), "rb") as fh:
                src = fh.read()
            r = subprocess.run([sys.executable, "-", "relay", "--facts-b64", b64], input=src,
                               cwd=self.home, env=e, capture_output=True, timeout=120)
        finally:
            m.close()
        self.assertEqual(r.returncode, 0, r.stderr.decode())
        self.assert_untouched(r)
        self.assertEqual(r.stdout.decode().strip().splitlines()[-1], "sent")
        self.assertEqual(len(m.bodies), 1)

    def test_no_hive_exec_and_no_version_are_found_beside_nothing(self):
        """`python3 - --test --dry-run` with a HIVE login: HERE used to be the current directory,
        so ~/hive_exec.sh became the relay and ~/../.claude-plugin/plugin.json the version."""
        # SKILL_SLACK on (job_env turns it off): the relay decision must be reached
        r = self.piped("notify_slack.py", "--test", "--dry-run", SKILL_SLACK=None,
                       HIVE_USER="tester", HIVE_KEY=self.key,
                       SKILL_SLACK_GROUP_FILE=os.path.join(self.d, "none"),
                       SKILL_CORE_GROUP_DIR=os.path.join(self.d, "no_core_group"))
        self.assertEqual(r.returncode, 0, r.stderr.decode())
        self.assert_untouched(r)
        out = (r.stdout + r.stderr).decode()
        self.assertIn("dry run, nothing sent", out)
        self.assertNotIn("would relay through HIVE", out)
        self.assertNotIn("6.6.6-decoy", out)


class NothingFoundBesideNothing(Harness):
    def test_each_finder_finds_nothing_when_piped(self):
        """Mutants that put a current-directory fallback back into _relay_route,
        _skill_version, _helper or hive_exec_path passed every command-level test above."""
        r = subprocess.run([sys.executable, "-"], input=PROBE.encode(), cwd=self.home,
                           capture_output=True, timeout=120,
                           env=job_env(self.d, HOME=self.home, PROBE_SCRIPTS=SCRIPTS,
                                       HIVE_USER="tester", HIVE_KEY=self.key,
                                       SKILL_CORE_GROUP_DIR=os.path.join(self.d, "no_core")))
        self.assertEqual(r.returncode, 0, r.stderr.decode())
        self.assert_untouched(r)
        got = json.loads(r.stdout.decode().strip().splitlines()[-1])
        self.assertEqual((got["ns_HERE"], got["rr_HERE"]), (None, None))
        self.assertIsNone(got["relay_route"])
        self.assertEqual(got["skill_version"], got["unknown"])
        self.assertIsNone(got["helper"])
        self.assertEqual(got["hive_exec_path"], "")
        self.assertIsNone(got["rr_skill_version"])
        self.assertFalse(got["notify_slack.py cwd on sys.path"])
        self.assertFalse(got["record_run.py cwd on sys.path"])


class FromAFile(unittest.TestCase):
    def test_run_from_a_file_here_is_its_folder(self):
        for script in ("record_run.py", "notify_slack.py"):
            with self.subTest(script=script):
                code = (f"import runpy, sys; sys.argv=['x']; "
                        f"g = runpy.run_path({os.path.join(SCRIPTS, script)!r}, "
                        f"run_name='not_main'); print(g['HERE'])")
                with tempfile.TemporaryDirectory() as tmp:
                    r = subprocess.run([sys.executable, "-c", code], cwd=tmp,
                                       capture_output=True, text=True, timeout=60,
                                       env=job_env(tmp))
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual(os.path.realpath(r.stdout.strip()), os.path.realpath(SCRIPTS))

    def test_the_two_copies_of_the_block_are_identical(self):
        """Both files go to HIVE alone, so neither can import the other's: one text, twice."""
        def block(name):
            src = read(os.path.join(SCRIPTS, name))
            return src[src.index("# ---- where this file is (begin"):
                       src.index("# ---- where this file is (end)")]
        self.assertEqual(block("record_run.py"), block("notify_slack.py"))
        for name in ("record_run.py", "notify_slack.py"):
            src = read(os.path.join(SCRIPTS, name))
            self.assertLess(src.index("# ---- where this file is (end)"),
                            src.index("\nimport argparse"), f"{name}: the block runs first")
            self.assertNotIn('"__file__" in globals()', src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
