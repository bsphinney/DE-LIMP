#!/usr/bin/env python3
"""
detect_acquisition.py has ONE place for "the user must hear this before a search starts":
each file's `warnings` list. Every entry is printed to stderr with the file named, and any
entry sets `needs_confirmation`.

Two readers fill it, and they were written on separate branches:
  - a Thermo .raw the parser could not read, or read only in part (no instrument, no
    isolation window edges, parser missing) -- test_thermo_raw_detection.py;
  - a Bruker .d whose analysis.tdf is not `ok` (truncated, stale side file, WAL header,
    ...) -- bruker_tdf.tdf_integrity(), test_tdf_readonly_open.py.

Merged naively the two overwrite each other: the Bruker side assigned
`warnings = [integrity warning] or []` for every file, which empties the Thermo reader's
list for a .raw -- the parser failure then prints nothing and no longer gates. And the
Bruker side's top-level `tdf_integrity_problem_files` was a second, vendor-specific list
to check, which says "nothing wrong" for a cohort whose only problem is a .raw.

These tests read a mixed cohort through the CLI, so both readers run in one call.
POSIX only, like the Thermo tests: the parser stand-in is a shell shim on PATH.
"""
import json
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import detect_acquisition as da          # noqa: E402
# Module imports, not `from ... import`: a TestCase class bound in this namespace would be
# collected and run a second time here.
import test_thermo_raw_detection as trfp  # noqa: E402
import test_tdf_readonly_open as tdf      # noqa: E402

FILE_KEYS = {"file", "vendor", "acquisition", "confidence", "reason", "instrument",
             "precursor_mz_range", "tdf_integrity", "warnings", "reader"}


class _MixedCohort(trfp._FakeParserCase):
    def _run(self, *paths):
        res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                              *paths], capture_output=True, text=True, env=os.environ.copy())
        self.assertEqual(res.returncode, 0, res.stderr)
        return json.loads(res.stdout), res.stderr


class EveryProblemIsInItsFilesWarnings(_MixedCohort):
    def test_a_raw_failure_and_a_truncated_tdf_each_warn_on_their_own_file(self):
        os.environ["FAKE_TRFP_FAIL"] = "metadata"      # the .raw classifies, instrument unread
        raw = self.raw(trfp.EXPLORIS)
        good = tdf.make_intact_d(self.tmp)
        bad = tdf.truncate_like_hive(tdf.make_stale_wal_d(self.tmp, name="bad.d"))
        out, err = self._run(raw, good, bad)
        by_file = {f["file"]: f for f in out["files"]}

        self.assertTrue(any("instrument unknown" in w for w in by_file[raw]["warnings"]),
                        by_file[raw]["warnings"])
        self.assertFalse(any("analysis.tdf" in w for w in by_file[raw]["warnings"]))
        self.assertEqual(len(by_file[bad]["warnings"]), 1, by_file[bad]["warnings"])
        self.assertIn("analysis.tdf truncated", by_file[bad]["warnings"][0])
        self.assertEqual(by_file[good]["warnings"], [])

        # a person watching the terminal sees both, each with its file named
        self.assertIn(f"WARNING: {raw}: instrument unknown", err)
        self.assertIn(f"WARNING: {bad}: analysis.tdf truncated", err)
        self.assertNotIn(f"WARNING: {good}:", err)

        # nothing else would have asked: one DIA acquisition, one instrument, all high
        self.assertEqual(out["overall"], "DIA")
        self.assertEqual(out["instruments_seen"], ["timsTOF HT"])
        self.assertEqual(out["low_confidence_files"], [])
        self.assertTrue(out["needs_confirmation"])

    def test_a_raw_warning_is_kept_beside_healthy_bruker_runs(self):
        """The overwrite bug: a .d-only integrity assignment emptied the .raw's warnings.
        (needs_confirmation is set here by the unknown acquisition too; the test above is
        the one where only a warning can set it.)"""
        self.hide_parser()
        raw = self.raw(trfp.EXPLORIS)
        good = tdf.make_intact_d(self.tmp)
        out, err = self._run(good, raw)
        by_file = {f["file"]: f for f in out["files"]}
        self.assertTrue(any("not found" in w for w in by_file[raw]["warnings"]),
                        by_file[raw]["warnings"])
        self.assertIn(f"WARNING: {raw}: ThermoRawFileParser not found", err)
        self.assertEqual(by_file[good]["tdf_integrity"]["status"], "ok")
        self.assertEqual(by_file[good]["warnings"], [])

    def test_a_tdf_warning_is_kept_beside_a_raw_that_read_cleanly(self):
        raw = self.raw(trfp.EXPLORIS)
        d = tdf.make_stale_wal_d(self.tmp, rollback_header=True)
        out, err = self._run(raw, d)
        by_file = {f["file"]: f for f in out["files"]}
        self.assertEqual(by_file[raw]["warnings"], [])
        self.assertEqual(by_file[d]["tdf_integrity"]["status"], "stale_side_file")
        self.assertIn(f"WARNING: {d}: analysis.tdf stale_side_file", err)
        self.assertNotIn(f"WARNING: {raw}:", err)


class OneShapeForEveryFile(_MixedCohort):
    def test_every_file_record_has_the_same_fields_whatever_its_vendor(self):
        raw = self.raw(trfp.EXPLORIS)
        d = tdf.make_intact_d(self.tmp)
        mzml = os.path.join(self.tmp, "run.mzML")
        open(mzml, "w").close()
        out, _ = self._run(raw, d, mzml)
        by_file = {f["file"]: f for f in out["files"]}
        for f in out["files"]:
            self.assertEqual(set(f), FILE_KEYS, f["file"])
            self.assertIsInstance(f["warnings"], list, f["file"])
        # each reader's own detail stays on the files it applies to
        self.assertIsNone(by_file[raw]["tdf_integrity"])
        self.assertIsNotNone(by_file[raw]["reader"])
        self.assertEqual(by_file[d]["tdf_integrity"]["status"], "ok")
        self.assertIsNone(by_file[d]["reader"])
        self.assertIsNone(by_file[mzml]["tdf_integrity"])
        self.assertIsNone(by_file[mzml]["reader"])

    def test_there_is_no_second_top_level_list_of_problems(self):
        """`files[].warnings` is the list to read. A top-level list for one vendor's problems
        reads as "nothing wrong" whenever the problem is the other vendor's."""
        self.hide_parser()
        raw = self.raw(trfp.EXPLORIS)
        bad = tdf.truncate_like_hive(tdf.make_stale_wal_d(self.tmp, name="bad.d"))
        out, _ = self._run(raw, bad)
        self.assertNotIn("warnings", out)
        self.assertNotIn("tdf_integrity_problem_files", out)
        self.assertEqual(sorted(f["file"] for f in out["files"] if f["warnings"]),
                         sorted([raw, bad]))


class BrukerChecksOnALoginNode(_MixedCohort):
    """The login-node refusal counts Thermo .raw only. The tdf check is a header read, one
    indexed query and a 4-byte read, so it must neither be refused nor skipped there."""

    def setUp(self):
        super().setUp()
        sbin = os.path.join(self.tmp, "sbin")
        os.makedirs(sbin)
        with open(os.path.join(sbin, "sbatch"), "w") as fh:
            fh.write("#!/bin/sh\nexit 0\n")
        os.chmod(os.path.join(sbin, "sbatch"), 0o755)
        os.environ["PATH"] = sbin + os.pathsep + os.environ["PATH"]
        os.environ.pop("SLURM_JOB_ID", None)

    def test_a_cohort_of_bruker_runs_is_checked_and_warns_on_a_login_node(self):
        self.assertTrue(da.on_cluster_login_node())
        runs = [tdf.truncate_like_hive(tdf.make_stale_wal_d(self.tmp, name=f"run{i}.d"))
                for i in range(da.LOGIN_NODE_MAX_RAW + 1)]
        out, err = self._run(*runs, self.raw(trfp.EXPLORIS))
        self.assertEqual(len(out["files"]), len(runs) + 1)
        for d in runs:
            self.assertIn(f"WARNING: {d}: analysis.tdf truncated", err)
        self.assertTrue(out["needs_confirmation"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
