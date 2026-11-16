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
import shutil
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


def _run(tmp, manifest, extra=None, ht=HT):
    env = dict(os.environ)
    env["FAKE_MANIFEST"] = json.dumps(manifest)
    env["PGPASSWORD"] = "not-a-real-token"          # skips the credential path entirely
    argv = [sys.executable, ht, "--stan", _fake_stan(tmp), "fetch", "9001", "--out", tmp]
    return subprocess.run(argv + (extra or []), capture_output=True, text=True, env=env)


def _well(i):
    return "ABCDEFGH"[i // 12 % 8] + str(i % 12 + 1)


def _touch(tmp, name):
    f = os.path.join(tmp, name)
    open(f, "w").close()
    return f


def _plate_run(tmp, i, prefix="20260827_9001_100spd", plate="S5", inj=24000):
    """One of the submission's own plate runs, named as the Core names them:
    <date>_[PROT_]<n>_<speed>_<sample name>_S<plate>-<well>_1_<acquisition #>.d"""
    return _touch(tmp, f"{prefix}_XY{i + 1}_{plate}-{_well(i)}_1_{inj + i}.d")


def _ok(tmp, n=24):
    """A manifest that should pass every gate. Paths must really exist -- the paths_exist
    gate checks the filesystem, so the fixture creates them."""
    files = [_plate_run(tmp, i) for i in range(n)]
    return {"submission": "9001", "include": "samples", "files": files, "n_files": n,
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
            m["files"].append(os.path.join(tmp, "20260827_9001_100spd_XY99_S6-A1_1_24999.d"))
            self.assertEqual(_run(tmp, m).returncode, 2)
            g = {x["gate"]: x for x in
                 json.load(open(os.path.join(tmp, "ht_manifest.json")))["gates"]}
            self.assertEqual(g["paths_exist"]["status"], "FAIL")

    def test_stale_stan_exits_3_not_2(self):
        """A STAN too old to have ht-manifest must be distinguishable from a bad
        submission -- 3 means 'fix STAN', 2 means 'the data is wrong'."""
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(os.environ, PGPASSWORD="x", FAKE_MANIFEST="{}")
            fake = _fake_stan(tmp)
            r = subprocess.run([sys.executable, HT, "--stan", fake, "fetch", "9001",
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
            name = os.path.basename(twin)[:-2]
            self.assertEqual(out["repeated_names"], [{"run_name": name,
                                                      "paths": [m["files"][2], twin]}])
            g = {x["gate"]: x for x in out["gates"]}
            self.assertEqual(g["repeated_names"]["status"], "WARN")
            self.assertIn("run_search.py stops before searching them as given",
                          g["repeated_names"]["detail"])
            self.assertIn(f"FLAG: 1 run name(s) are shared by different files, all kept: {name}: "
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
            self.assertEqual([d["run_name"] for d in out["repeated_names"]],
                             [os.path.basename(twin)[:-2]])
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


def _summary(tmp, names, internal_id="PROT_9001"):
    """core_submission.py fetch's submission_summary.json, as much of it as the check reads: the
    CoreOmics record's own sample list (an HT plate run carries the sample_name)."""
    p = os.path.join(tmp, "submission_summary.json")
    with open(p, "w") as fh:
        json.dump({"internal_id": internal_id, "id": "0123456789ab",
                   "samples": [{"unique_id": f"U{i + 1:03d}", "sample_name": n,
                                "condition_name": "b" if i % 2 else "a"}
                               for i, n in enumerate(names)]}, fh)
    return p


class TheFilesAreThisSubmissions(unittest.TestCase):
    """2026-10-05 (2.10.0): a 96-sample HT plate came back as 136 files with every gate PASS: 95
    of its 96 runs (one real run missing) and 41 runs of OTHER projects that had used the same
    tray on earlier days -- named another way, other dates, other methods, 5 of them DDA. STAN
    assigns runs to a submission by tray and injection window, whoever's they are. Searched as
    given, other clients' data enters the search, the report and FRAN, and a sample is silently
    left out. Caught by hand. These assert on the exit code and on what files.txt holds."""

    def gates(self, tmp):
        with open(os.path.join(tmp, "ht_manifest.json")) as fh:
            return {g["gate"]: g for g in json.load(fh)["gates"]}

    def listed(self, tmp):
        with open(os.path.join(tmp, "files.txt")) as fh:
            return [l for l in fh.read().split("\n") if l]

    def shared_tray(self, tmp, missing=68):
        """The plate: 96 samples on tray S3, one of them never recorded by STAN; earlier on the
        same tray, 36 DIA and 5 DDA runs of another project, named another way."""
        dia = [_touch(tmp, f"09282026__60SPD_DIA-QZ{k + 1}_S3-{_well(k)}_1_{24700 + k}.d")
               for k in range(36)]
        dda = [_touch(tmp, f"09292026__30SPD_DDAQZ{k + 1}_S3-{_well(40 + k)}_1_{24800 + k}.d")
               for k in range(5)]
        own = [_plate_run(tmp, i, prefix="20261005_PROT_9001_100spd", plate="S3", inj=25000)
               for i in range(96) if i != missing]
        files = dia + dda + own
        entries = [{"run_name": os.path.basename(f)[:-2], "raw_path": f, "class": "sample",
                    "plate": "S3", "well": None, "injection": int(f[:-2].rsplit("_", 1)[1]),
                    "needs_rerun": False, "verdict": "pass"} for f in files]
        m = {"submission": "9001", "include": "samples", "files": files, "entries": entries,
             "n_files": len(files), "plates": ["S3"], "counts": {"sample": len(files)},
             "n_needs_rerun": 0, "missing_paths": []}
        return m, own, dia + dda

    def test_a_tray_shared_with_earlier_projects_stops_and_names_every_foreign_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            m, own, foreign = self.shared_tray(tmp)
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 2, r.stdout)
            g = self.gates(tmp)
            f = g["named_for_submission"]
            self.assertEqual((f["status"], f["n"]), ("FAIL", 41))
            self.assertEqual(sorted(x["file"] for x in f["files"]), sorted(foreign))
            self.assertEqual({x["class"] for x in f["files"]}, {"sample"})
            self.assertIn("STAN assigns runs to a submission by tray", f["detail"])
            self.assertIn("41 of them ran before", f["detail"])
            for x in foreign:
                self.assertIn(x, r.stdout)                      # every one, not five examples
            self.assertEqual(self.listed(tmp), own)             # never another project's run
            self.assertEqual(g["n_files"]["n"], 95)
            meth = g["one_method"]
            self.assertEqual(meth["status"], "FAIL")
            self.assertEqual({x["method"]: (x["n"], x["n_named_for_submission"])
                              for x in meth["groups"]},
                             {"100spd": (95, 95), "60spd DIA": (36, 0), "30spd": (5, 0)})
            self.assertIn("files      : 95 (136 listed by STAN; 41 not named for submission 9001)",
                          r.stdout)
            self.assertIn("HARD GATE FAILED", r.stderr)

    def test_the_records_sample_list_names_the_missing_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            m, _own, _foreign = self.shared_tray(tmp)
            s = _summary(tmp, [f"XY{i + 1}" for i in range(96)])
            r = _run(tmp, m, ["--summary", s])
            self.assertEqual(r.returncode, 2)
            g = self.gates(tmp)["sample_list"]
            self.assertEqual(g["status"], "FAIL")
            self.assertEqual([(x["unique_id"], x["sample_name"]) for x in g["missing"]],
                             [("U069", "XY69")])
            self.assertEqual(g["extra"], [])
            self.assertIn("U069 (XY69)", r.stdout)

    def test_a_plate_that_is_all_its_own_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            r = _run(tmp, m, ["--summary", _summary(tmp, [f"XY{i + 1}" for i in range(24)])])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = self.gates(tmp)
            for k in ("named_for_submission", "one_method", "sample_list"):
                self.assertEqual(g[k]["status"], "PASS", g[k])
            self.assertEqual(self.listed(tmp), m["files"])

    def test_a_run_of_no_sample_on_the_record_is_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            r = _run(tmp, m, ["--summary", _summary(tmp, [f"XY{i + 1}" for i in range(23)])])
            self.assertEqual(r.returncode, 2)
            g = self.gates(tmp)["sample_list"]
            self.assertEqual((g["status"], g["missing"]), ("FAIL", []))
            self.assertEqual([x["file"] for x in g["extra"]], [m["files"][23]])

    def test_the_record_of_another_submission_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            s = _summary(tmp, [f"XY{i + 1}" for i in range(24)], internal_id="PROT_9002")
            r = _run(tmp, m, ["--summary", s])
            self.assertEqual(r.returncode, 2)
            g = self.gates(tmp)["sample_list"]
            self.assertEqual(g["status"], "FAIL")
            self.assertIn("PROT_9002", g["detail"])

    def test_without_the_record_the_sample_list_is_not_checked_and_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(_run(tmp, _ok(tmp)).returncode, 0)
            g = self.gates(tmp)["sample_list"]
            self.assertEqual(g["status"], "WARN")
            self.assertIn("--summary", g["detail"])

    def test_a_rerun_list_is_not_held_to_the_whole_sample_list(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = dict(_ok(tmp, n=3), include="rerun")
            s = _summary(tmp, [f"XY{i + 1}" for i in range(24)])
            self.assertEqual(_run(tmp, m, ["--include", "rerun", "--summary", s]).returncode, 0)
            self.assertEqual(self.gates(tmp)["sample_list"]["status"], "INFO")

    def test_what_counts_as_named_for_the_submission(self):
        """core_submission.ht_pattern, the one definition (2.9.1): the number right after the
        leading date, bare or PROT_, any case, and ending there."""
        with tempfile.TemporaryDirectory() as tmp:
            ours = [_touch(tmp, n) for n in (
                "20260827_9001_100spd_XY1_S5-A1_1_24000.d",
                "20260827_PROT_9001_100spd_XY2_S5-A2_1_24001.d",
                "20260827_prot9001_100spd_XY3_S5-A3_1_24002.d",
                "20260901_9001_rerun_XY4_S1-A1_1_24200.d")]
            not_ours = [_touch(tmp, n) for n in (
                "20260827_90010_100spd_XY5_S5-A5_1_24003.d",       # another submission
                "20260828_100spd_XY6_S5-A6_1_29001.d")]            # 9001 only in the counter
            m = dict(_ok(tmp, n=0), files=ours + not_ours, counts={"sample": 6})
            self.assertEqual(_run(tmp, m).returncode, 2)
            self.assertEqual(self.listed(tmp), ours)
            self.assertEqual(sorted(x["file"] for x in self.gates(tmp)["named_for_submission"]["files"]),
                             sorted(not_ours))

    def test_runs_that_declare_different_methods_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp, n=12)
            m["files"] += [_plate_run(tmp, i, prefix="20260828_9001_60SPD", inj=24100)
                           for i in range(12, 16)]
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 2)
            g = self.gates(tmp)
            self.assertEqual(g["named_for_submission"]["status"], "PASS")
            self.assertEqual(g["one_method"]["status"], "FAIL")
            self.assertEqual({x["method"]: x["n"] for x in g["one_method"]["groups"]},
                             {"100spd": 12, "60spd": 4})

    def test_dia_beside_dda_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            files = [_touch(tmp, "20260827_9001_60SPD_DIA-XY1_S5-A1_1_24000.d"),
                     _touch(tmp, "20260827_9001_60-spd-dda-XY2_S5-A2_1_24001.d")]
            m = dict(_ok(tmp, n=0), files=files, counts={"sample": 2})
            self.assertEqual(_run(tmp, m).returncode, 2)
            g = self.gates(tmp)["one_method"]
            self.assertEqual(g["status"], "FAIL")
            self.assertEqual({x["method"] for x in g["groups"]}, {"60spd DIA", "60spd DDA"})

    def test_a_run_whose_name_declares_no_method_is_flagged_not_stopped(self):
        """A re-injection is named <date>_<n>_rerun_<sample>_...: no speed in the name."""
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp)
            rerun = _touch(tmp, "20260901_9001_rerun_XY1_S1-A1_1_24200.d")
            m["files"].append(rerun)
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = self.gates(tmp)["one_method"]
            self.assertEqual(g["status"], "WARN")
            self.assertEqual(g["undeclared"], [rerun])

    def test_no_file_named_for_the_submission_points_at_a_curated_list(self):
        """A continuation tray named without the number: STAN joined it by the acquisition
        counter, and the name cannot say whose it is. Staff choose the files instead."""
        with tempfile.TemporaryDirectory() as tmp:
            files = [_touch(tmp, f"20260828_100spd_XY{i + 1}_S5-{_well(i)}_1_{24126 + i}.d")
                     for i in range(4)]
            m = dict(_ok(tmp, n=0), files=files, counts={"sample": 4})
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 2)
            g = self.gates(tmp)["n_files"]
            self.assertEqual((g["status"], g["n"]), ("FAIL", 0))
            self.assertIn("core_submission.py locate --files-from", g["detail"])
            self.assertEqual(self.listed(tmp), [])

    def test_without_the_name_rule_nothing_passes(self):
        """A partial copy of scripts/ (no core_submission.py, no fran_deposit.py): the names cannot
        be checked, so the check FAILS -- it never passes for want of its rule."""
        with tempfile.TemporaryDirectory() as tmp:
            part = os.path.join(tmp, "partial")
            os.makedirs(part)
            for f in ("ht_manifest.py", "check_report_runs.py"):
                shutil.copy(os.path.join(SCRIPTS, f), part)
            r = _run(tmp, _ok(tmp), ht=os.path.join(part, "ht_manifest.py"))
            self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
            g = self.gates(tmp)
            self.assertEqual(g["named_for_submission"]["status"], "FAIL")
            self.assertIn("core_submission.py could not be imported",
                          g["named_for_submission"]["detail"])
            self.assertEqual(g["one_method"]["status"], "FAIL")
            self.assertEqual(self.listed(tmp), [])


class Review212(unittest.TestCase):
    """Skill 2.12 review (error handling): three ways a sound plate was stopped, or stopped with
    the wrong advice."""

    def gates(self, tmp):
        with open(os.path.join(tmp, "ht_manifest.json")) as fh:
            return {g["gate"]: g for g in json.load(fh)["gates"]}

    # ---- runs that name several samples are staff's choice, not runs of no sample -----------
    def two_runs_named_xy1(self, tmp):
        """Two samples on the record share the sample_name XY1 (a submitter's duplicate), and
        the plate has two runs named XY1."""
        m = _ok(tmp, n=4)
        new = m["files"][1].replace("_XY2_", "_XY1_")
        os.rename(m["files"][1], new)
        m["files"][1] = new
        return m, _summary(tmp, ["XY1", "XY1", "XY3", "XY4"])

    def test_runs_naming_two_samples_ask_staff_to_choose(self):
        with tempfile.TemporaryDirectory() as tmp:
            m, s = self.two_runs_named_xy1(tmp)
            r = _run(tmp, m, ["--summary", s])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = self.gates(tmp)["sample_list"]
            self.assertEqual((g["status"], g["missing"], g["extra"]), ("INFO", [], []))
            self.assertEqual(g["several_runs"], ["U001", "U002"])
            self.assertEqual([(x["file"], x["samples"]) for x in g["runs_of_several_samples"]],
                             [(f, ["U001", "U002"]) for f in m["files"][:2]])
            self.assertIn("2 run(s) name several samples", g["detail"])
            self.assertIn("core_submission.py locate --files-from", g["detail"])
            self.assertIn("--choose <unique_id>=<file>", g["detail"])
            self.assertNotIn("belong to no sample", g["detail"])
            self.assertIn(f"several samples (U001, U002): {m['files'][0]}", r.stdout)

    def test_a_run_of_no_sample_still_stops_beside_them(self):
        with tempfile.TemporaryDirectory() as tmp:
            m, _s = self.two_runs_named_xy1(tmp)
            stray = _touch(tmp, "20260827_9001_100spd_QQ9_S5-A9_1_24008.d")
            m["files"].append(stray)
            r = _run(tmp, m, ["--summary", _summary(tmp, ["XY1", "XY1", "XY3", "XY4"])])
            self.assertEqual(r.returncode, 2)
            g = self.gates(tmp)["sample_list"]
            self.assertEqual(g["status"], "FAIL")
            self.assertEqual([x["file"] for x in g["extra"]], [stray])
            self.assertEqual(len(g["runs_of_several_samples"]), 2)
            self.assertIn("1 run(s) named for it belong to no sample", g["detail"])
            self.assertIn("Separately, 2 run(s) name several samples", g["detail"])

    # ---- the method is read right after the plate token, never in the sample name ------------
    def test_a_timepoint_in_a_reinjected_sample_name_is_not_a_method(self):
        """An aging study: `..._rerun_Liver_24m_3_...` declared method "24m" and the plate was
        stopped as "2 acquisition methods (100spd, 24m)"."""
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp, n=12)
            rr = _touch(tmp, "20260901_9001_rerun_Liver_24m_3_S5-B2_1_24100.d")
            m["files"].append(rr)
            r = _run(tmp, m)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            g = self.gates(tmp)["one_method"]
            self.assertEqual((g["status"], g["undeclared"]), ("WARN", [rr]))
            self.assertEqual({x["method"]: x["n"] for x in g["groups"]}, {"100spd": 12})

    def test_declared_method_reads_only_the_token_after_the_plate(self):
        sys.path.insert(0, SCRIPTS)
        import core_submission
        import ht_manifest
        plate = core_submission.ht_pattern("PROT_9001")
        cases = {
            "20260901_9001_rerun_Liver_24m_3_S5-B2_1_24100.d": None,
            "20260901_9793_rerun_XY4_S5-A4_1_24100.d": None,               # another plate: whole name
            "20261005_PROT_9001_100spd_XY1_S3-A1_1_25000.d": "100spd",
            "20260827_9001_100spd_Liver_24m_3_S5-A1_1_24000.d": "100spd",
            "20260827_9001_60SPD_DIA-XY1_S5-A1_1_24000.d": "60spd DIA",
            "20260827_9001_60-spd-dda-XY2_S5-A2_1_24001.d": "60spd DDA",
            "20260827_9001_Liver_90m_S5-A3_1_24002.d": None,               # a method in the sample
            # not named for the submission: read anywhere (never searched; it only describes it)
            "09282026__60SPD_DIA-QZ1_S3-A1_1_24700.d": "60spd DIA",
            "20260827_9002_Liver_90m_S5-A3_1_24002.d": "90m",
        }
        for name, want in cases.items():
            with self.subTest(name=name):
                self.assertEqual(ht_manifest.declared_method(name, plate), want)

    # ---- a missing --summary is a file to put on HIVE, not a question for the operator --------
    def test_a_missing_summary_says_put_it_on_hive_not_ask_the_operator(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = _run(tmp, _ok(tmp, n=12), ["--summary", os.path.join(tmp, "nope", "s.json")])
            self.assertEqual(r.returncode, 2)
            g = self.gates(tmp)["sample_list"]
            self.assertEqual((g["status"], g["summary_unreadable"]), ("FAIL", True))
            self.assertIn("Put submission_summary.json on HIVE and pass --summary", g["detail"])
            self.assertIn("~/core/PROT_9001/hive/submission_summary.json", g["detail"])
            self.assertNotIn("operator", r.stderr)
            self.assertIn("put submission_summary.json on HIVE and pass --summary", r.stderr)
            self.assertNotIn("Once the operator has said", r.stderr)

    def test_with_another_failed_gate_both_are_said(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = _ok(tmp, n=12)
            m["files"].append(_touch(tmp, "09282026__60SPD_DIA-QZ1_S3-A1_1_24700.d"))
            r = _run(tmp, m, ["--summary", os.path.join(tmp, "nope.json")])
            self.assertEqual(r.returncode, 2)
            self.assertIn("Resolve the above with the operator first.", r.stderr)
            self.assertIn("put submission_summary.json on HIVE and pass --summary", r.stderr)

    def test_a_summary_path_given_with_a_tilde_is_read(self):
        from unittest import mock
        with tempfile.TemporaryDirectory() as tmp:
            s = _summary(tmp, [f"XY{i + 1}" for i in range(12)])
            with mock.patch.dict(os.environ, {"HOME": tmp}):
                r = _run(tmp, _ok(tmp, n=12), ["--summary", "~/" + os.path.basename(s)])
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertEqual(self.gates(tmp)["sample_list"]["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
