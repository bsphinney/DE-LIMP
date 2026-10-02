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


class CredentialIsNeverShared(unittest.TestCase):
    """The owner's .pgfarm_token is the service account's long-lived SECRET (STAN CLAUDE.md:
    512 bytes, mode 0600), not a 7-day token. Skill 2.9 called it safe to copy group-readable to
    /quobyte/proteomics-grp/etc/pgfarm_token and searched there first."""

    def setUp(self):
        sys.path.insert(0, SCRIPTS)
        import ht_manifest
        self.ht = ht_manifest

    def test_no_shared_copy_is_searched(self):
        self.assertFalse(hasattr(self.ht, "SHARED_TOKEN"))
        self.assertFalse(any("/etc/" in c for c in self.ht.TOKEN_CANDIDATES))
        with open(HT, encoding="utf-8") as fh:
            src = fh.read()
        # the old wording may be quoted as history, never stated
        self.assertFalse("7-day token" in src.replace('"7-day token"', ""),
                         "ht_manifest.py calls the secret a 7-day token")
        self.assertFalse("group-readable copy at" in src.split("TOKEN_CANDIDATES =")[1],
                         "ht_manifest.py still suggests a group-readable copy")

    def test_no_credential_says_never_copy_and_points_at_http(self):
        from unittest import mock
        with tempfile.TemporaryDirectory() as tmp:
            env = {k: v for k, v in os.environ.items() if k not in ("PGPASSWORD", "STAN_PG_TOKEN")}
            with mock.patch.object(self.ht, "TOKEN_CANDIDATES", (os.path.join(tmp, "none"),)), \
                    mock.patch.dict(os.environ, env, clear=True):
                with self.assertRaises(SystemExit) as cm:
                    self.ht._env(None)
        msg = str(cm.exception.code)
        self.assertIn("NEVER copy it or make it group-readable", msg)
        self.assertIn("--http https://ucd.stan-proteomics.org --share-token-file", msg)
        for wrong in ("7-day", "publish a group-readable", "etc/pgfarm_token"):
            self.assertNotIn(wrong, msg)

    def test_the_reference_no_longer_suggests_sharing_it(self):
        ref = os.path.join(os.path.dirname(SCRIPTS), "references", "ht-submissions.md")
        with open(ref, encoding="utf-8") as fh:
            doc = fh.read()
        self.assertNotIn("7-day token", doc.replace('"7-day token"', ""))
        self.assertNotIn("publish here", doc)
        self.assertIn("Never copy it, and never make it group-readable", doc)


def _entries(files, well=lambda i: f"A{i + 1}"):
    """STAN's per-run entries, shaped like the live payload (2026-10-01): one per line of
    `files`, a repeated file repeating its entry exactly."""
    first = {}
    out = []
    for f in files:
        i = first.setdefault(f, len(first))
        out.append({"run_name": os.path.basename(f)[:-2], "raw_path": f, "class": "sample",
                    "plate": "S5", "well": well(i), "injection": 24000 + i,
                    "needs_rerun": False, "verdict": "pass"})
    return out


class TestRepeatedRuns(unittest.TestCase):
    """Brett (2026-10-01): repeats in a run list are okay, but they are flagged. STAN listed one
    run 4 times and another twice -- 100 lines for 96 files -- with every gate PASS (a staff
    plate). The same file goes into files.txt once and is flagged with its count; different
    files sharing a run name are all kept and flagged. Neither stops ht_manifest."""

    def manifest(self, tmp):
        with open(os.path.join(tmp, "ht_manifest.json")) as fh:
            return json.load(fh)

    def listed(self, tmp):
        with open(os.path.join(tmp, "files.txt")) as fh:
            return [l for l in fh.read().split("\n") if l]

    def rerun_in_another_folder(self, tmp, f):
        os.makedirs(os.path.join(tmp, "rerun"))
        twin = os.path.join(tmp, "rerun", os.path.basename(f))
        open(twin, "w").close()
        return twin

    def test_one_path_listed_twice_is_listed_once_and_flagged(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            f = m["files"]
            m["files"] = f + [f[3], f[3], f[3], f[7]]
            m["entries"], m["n_files"] = _entries(m["files"]), len(m["files"])
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(self.listed(tmp), f)
            out = self.manifest(tmp)
            self.assertEqual((out["n_files_listed"], out["n_files"]), (28, 24))
            self.assertEqual([(d["path"], d["times"]) for d in out["repeated_paths"]],
                             [(f[3], 4), (f[7], 2)])
            self.assertEqual(out["repeated_names"], [])
            g = {x["gate"]: x for x in out["gates"]}
            self.assertEqual((g["repeated_paths"]["status"], g["repeated_paths"]["n"]), ("WARN", 4))
            self.assertEqual(g["repeated_names"]["status"], "PASS")
            self.assertEqual(g["n_files"]["n"], 24)
            self.assertIn(f"[ht_manifest] FLAG: 2 input file(s) were listed more than once; each "
                          f"is searched once: {f[3]} (listed 4 times)", r.stderr)
            self.assertIn("files      : 24 (28 listed by STAN; repeats removed)", r.stdout)

    def test_one_name_in_two_folders_keeps_both_and_flags_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            twin = self.rerun_in_another_folder(tmp, m["files"][2])
            m["files"].append(twin)
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn(twin, self.listed(tmp))
            self.assertIn(m["files"][2], self.listed(tmp))
            out = self.manifest(tmp)
            self.assertEqual(out["repeated_names"], [{"run_name": "run2",
                                                      "paths": [m["files"][2], twin]}])
            g = {x["gate"]: x for x in out["gates"]}
            self.assertEqual(g["repeated_names"]["status"], "WARN")
            self.assertIn("run_search.py stops before searching them as given",
                          g["repeated_names"]["detail"])
            self.assertIn(f"FLAG: 1 run name(s) are shared by different files, all kept: run2: "
                          f"{m['files'][2]}, {twin}", r.stderr)

    def test_a_mix_of_both_is_flagged_both_ways(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            f = list(m["files"])
            twin = self.rerun_in_another_folder(tmp, f[2])
            m["files"] = f + [f[5], twin, twin]
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(self.listed(tmp), f + [twin])
            out = self.manifest(tmp)
            self.assertEqual([(d["path"], d["times"]) for d in out["repeated_paths"]],
                             [(f[5], 2), (twin, 2)])
            self.assertEqual([d["run_name"] for d in out["repeated_names"]], ["run2"])
            self.assertEqual(r.stderr.count("[ht_manifest] FLAG:"), 2)

    def test_a_symlink_to_a_listed_run_is_the_same_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            alias = os.path.join(tmp, "alias.d")
            os.symlink(m["files"][0], alias)
            m["files"].append(alias)
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertNotIn(alias, self.listed(tmp))
            self.assertEqual(self.manifest(tmp)["repeated_paths"][0]["also_listed_as"], [alias])

    def test_entries_that_disagree_about_one_file_are_a_hard_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            f = m["files"]
            m["files"] = f + [f[5]]
            m["entries"] = _entries(m["files"])
            m["entries"][-1] = dict(m["entries"][-1], well="H12")
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 2)
            g = {x["gate"]: x for x in self.manifest(tmp)["gates"]}
            self.assertEqual(g["repeated_paths"]["status"], "FAIL")
            self.assertIn("not a repeat but a contradiction about which sample the file is",
                          g["repeated_paths"]["detail"])
            self.assertEqual(g["repeated_paths"]["examples"][0], f"{f[5]}: well=A6 vs well=H12")
            self.assertEqual(g["repeated_paths"]["conflicts"][0]["entries"],
                             [{"well": "A6"}, {"well": "H12"}])

    def test_no_repeat_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            m["entries"] = _entries(m["files"])
            self.assertEqual(_run(tmp, m).returncode, 0)
            g = {x["gate"]: x for x in self.manifest(tmp)["gates"]}
            self.assertEqual((g["repeated_paths"]["status"], g["repeated_names"]["status"]),
                             ("PASS", "PASS"))

if __name__ == "__main__":
    unittest.main()
