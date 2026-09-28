"""
The HT gates are the entire safety value of step 1a, so they get a test.

Each hard gate exists because it otherwise produces a search that SUCCEEDS while covering
the wrong files -- the failure mode that costs a plate-sized SLURM run and is invisible in
the results. A regression here does not throw; it quietly searches a subset. So these
assert on the exit code the orchestrator branches on, not on log text.

Stdlib only (unittest), matching the rest of the suite -- CI installs nothing.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

SCRIPTS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
HT = os.path.join(SCRIPTS, "ht_manifest.py")

# A stand-in for the `stan` binary: prints whatever manifest JSON we hand it, so the gate
# logic can be exercised without STAN, its database, or a credential.
FAKE_STAN = """#!/usr/bin/env python3
import json, sys, os
if "ht-manifest" not in sys.argv:
    sys.stderr.write("Usage: stan\\nNo such command '%s'.\\n" % sys.argv[1]); sys.exit(2)
sys.stdout.write(os.environ["FAKE_MANIFEST"])
"""


def _fake_stan(tmp):
    p = os.path.join(tmp, "stan")
    with open(p, "w") as fh:
        fh.write(FAKE_STAN)
    os.chmod(p, 0o755)
    return p


def _run(tmp, manifest, extra=None):
    env = dict(os.environ)
    env["FAKE_MANIFEST"] = json.dumps(manifest)
    env["PGPASSWORD"] = "not-a-real-token"          # skips the credential path entirely
    argv = [sys.executable, HT, "--stan", _fake_stan(tmp), "fetch", "0793", "--out", tmp]
    return subprocess.run(argv + (extra or []), capture_output=True, text=True, env=env)


def _ok(tmp, n=24):
    """A manifest that should pass every gate. Paths must really exist -- the paths_exist
    gate checks the filesystem, so the fixture creates them."""
    files = []
    for i in range(n):
        f = os.path.join(tmp, f"run{i}.d")
        open(f, "w").close()
        files.append(f)
    return {"submission": "0793", "include": "samples", "files": files, "n_files": n,
            "plates": ["S5", "S6"], "counts": {"sample": n, "standard": 8, "blank": 7},
            "n_needs_rerun": 0, "missing_paths": []}


class TestHardGates(unittest.TestCase):
    def test_clean_manifest_passes_and_writes_both_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = _run(tmp, _ok(tmp))
            self.assertEqual(r.returncode, 0, r.stderr)
            files_txt = os.path.join(tmp, "files.txt")
            self.assertTrue(os.path.exists(files_txt))
            self.assertTrue(os.path.exists(os.path.join(tmp, "ht_manifest.json")))
            # one absolute path per line, no blank last line surprises
            lines = [l for l in open(files_txt).read().split("\n") if l]
            self.assertEqual(len(lines), 24)
            self.assertTrue(all(os.path.isabs(l) for l in lines))

    def test_missing_paths_is_a_hard_fail(self):
        """Runs STAN knows about but has no path for are EXCLUDED from `files`, so
        proceeding searches a subset and reports success."""
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            m["missing_paths"] = ["/nfs/.../lost_1.d", "/nfs/.../lost_2.d"]
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 2, "missing_paths must block the search")

    def test_empty_file_list_is_a_hard_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            m["files"], m["n_files"], m["counts"] = [], 0, {"sample": 0}
            self.assertEqual(_run(tmp, m).returncode, 2)

    def test_nonexistent_path_is_a_hard_fail(self):
        """Otherwise a 120-file array dies partway through, hours in."""
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            m["files"].append(os.path.join(tmp, "never_written.d"))
            self.assertEqual(_run(tmp, m).returncode, 2)

    def test_stale_stan_exits_3_not_2(self):
        """A STAN too old to have ht-manifest must be distinguishable from a bad
        submission -- 3 means 'fix STAN', 2 means 'the data is wrong'."""
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, PGPASSWORD="x", FAKE_MANIFEST="{}")
            fake = _fake_stan(tmp)
            r = subprocess.run([sys.executable, HT, "--stan", fake, "fetch", "0793",
                                "--out", tmp, "--include", "standards"],
                               capture_output=True, text=True,
                               env=dict(env, FAKE_MANIFEST=""))
            # the fake only answers ht-manifest; with no manifest it emits invalid JSON
            self.assertEqual(r.returncode, 3)


class TestWarnGates(unittest.TestCase):
    """Warnings must NOT block -- a genuinely small or wide submission is legal. They must
    still be visible, because the orchestrator is told to surface every non-PASS gate."""

    def _gates(self, tmp, m):
        _run(tmp, m)
        return {g["gate"]: g for g in
                json.load(open(os.path.join(tmp, "ht_manifest.json")))["gates"]}

    def test_too_many_plates_warns_but_proceeds(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            m["plates"] = ["S1", "S2", "S3", "S4"]
            self.assertEqual(_run(tmp, m).returncode, 0)
            self.assertEqual(self._gates(tmp, m)["plates"]["status"], "WARN")

    def test_implausibly_few_samples_warns_but_proceeds(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp, n=3)
            m["counts"] = {"sample": 3}
            self.assertEqual(_run(tmp, m).returncode, 0)
            self.assertEqual(self._gates(tmp, m)["counts"]["status"], "WARN")

    def test_needs_rerun_is_reported_because_those_samples_are_included(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            m["n_needs_rerun"] = 6
            self.assertEqual(_run(tmp, m).returncode, 0)
            g = self._gates(tmp, m)
            self.assertIn("needs_rerun", g)
            self.assertEqual(g["needs_rerun"]["n"], 6)


if __name__ == "__main__":
    unittest.main()
