#!/usr/bin/env python3
"""2.15.1 HIVE real-data checks (2026-11-17, audit_20261117/hive_checks/HIVE_CHECKS_2151.md):
the findings fixed after them.

F1  --between-run auto on a DIA-NN --no-norm report gave quantile (46-72% wrong calls at
    4-10 m/z on PXD063416); it now gives median_ratio, like every quantity no engine normalised
    (build_maxlfq.R maxlfq_between_run; tests/test_maxlfq_between_run.py runs it on reports).
    Here: the rule's whole table, and the normalisation check's gate, which now compares the
    step auto APPLIED -- a check made before 2.15.1 drew its volcano with the old rule's step.
F2  the check summary's "Data check:" line said TRIPPED and hid NOT FULLY ASSESSED.
F4  Sage's parsimony groups with no peptide of their own are listed in a CSV beside the report.
F5  ensure_dotnet8.sh's default folder follows the skill's per-user folder, not always $HOME.
F6  Microsoft's dotnet-install.sh fails under a TMPDIR with a version-like folder (2.15.1).
F7  setup.json: ready_for.dia_raw says whether DIA-NN can read .raw, from what the file records.
F8  an unpinned conda Sage records its package version, not "env" (its binary says 0.14.6).
F9  setup.sh's conda calls get an absolute package cache and run from inside the prefix.
"""
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import experiment_type as et               # noqa: E402
import normalization_check as nc           # noqa: E402
from job_env import no_cluster_env         # noqa: E402

ENSURE = os.path.join(SCRIPTS, "ensure_dotnet8.sh")
SETUP = os.path.join(SCRIPTS, "setup.sh")
# dotnet-install.sh's own regex (extract_dotnet_package, folders_with_version_regex)
DOTNET_INSTALL_VERSION_RE = re.compile(r"^.*/[0-9]+\.[0-9]+[^/]+/")
BASH = shutil.which("bash")


def _exe(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)
    return path


# ------------------------------------------------------------------------------------- F1 rule
RULE_R = r'''
suppressMessages({source("build_maxlfq.R")})
d <- function(...) list(...)
cases <- list(
  dia_normalised   = maxlfq_between_run("auto", "normalised", NULL, d(no_norm = FALSE, source = "x")),
  dia_no_norm      = maxlfq_between_run("auto", "normalised", NULL, d(no_norm = TRUE, source = "x")),
  dia_unknown      = maxlfq_between_run("auto", "normalised", NULL, d(no_norm = NA, source = "x")),
  sage             = maxlfq_between_run("auto", "normalised", d(level = "peptide", engine_normalised = "false")),
  fragpipe_normed  = maxlfq_between_run("auto", "normalised", d(level = "protein", engine_normalised = "true")),
  adapted_unknown  = maxlfq_between_run("auto", "normalised", d(level = "protein", engine_normalised = "unknown")),
  no_norm_quantile = maxlfq_between_run("quantile", "normalised", NULL, d(no_norm = TRUE, source = "x")),
  raw              = maxlfq_between_run("auto", "raw", NULL))
cat(jsonlite::toJSON(lapply(cases, function(b) list(method = b$method, why = b$why,
  has_en = !is.null(b$engine_normalised))), auto_unbox = TRUE))
'''


@unittest.skipUnless(shutil.which("Rscript"), "needs R")
class TheAutoRuleAsksOneQuestion(unittest.TestCase):
    """Did the engine already normalise these quantities across runs? yes -> none; no ->
    median_ratio, whichever engine (DIA-NN --no-norm, Sage); not known -> DIA-NN none
    [confirm] (it normalises unless --no-norm), an adapted report median_ratio [confirm]."""

    @classmethod
    def setUpClass(cls):
        p = subprocess.run(["Rscript", "-e", RULE_R], cwd=SCRIPTS, capture_output=True,
                           text=True, timeout=120)
        if p.returncode != 0 and "jsonlite" in p.stderr:
            raise unittest.SkipTest("needs jsonlite")
        assert p.returncode == 0, p.stderr
        cls.r = json.loads(p.stdout)

    def test_the_table(self):
        got = {k: v["method"] for k, v in self.r.items()}
        self.assertEqual(got, {"dia_normalised": "none", "dia_no_norm": "median_ratio",
                               "dia_unknown": "none", "sage": "median_ratio",
                               "fragpipe_normed": "none", "adapted_unknown": "median_ratio",
                               "no_norm_quantile": "quantile", "raw": "none"})

    def test_the_why_says_what_ran_and_what_ran_before(self):
        self.assertIn("--no-norm", self.r["dia_no_norm"]["why"])
        self.assertIn("quantile up to 2.15.0", self.r["dia_no_norm"]["why"])
        self.assertIn("[confirm]", self.r["dia_unknown"]["why"])
        self.assertIn("[confirm]", self.r["adapted_unknown"]["why"])
        self.assertEqual(self.r["no_norm_quantile"]["why"], "--between-run quantile (chosen for this DE)")

    def test_a_dia_nn_answer_carries_no_engine_declaration(self):
        for k in ("dia_normalised", "dia_no_norm", "dia_unknown"):
            self.assertFalse(self.r[k]["has_en"], k)
        for k in ("sage", "fragpipe_normed", "adapted_unknown"):
            self.assertTrue(self.r[k]["has_en"], k)


# ------------------------------------------------------------------------------------- F1 gate
DESIGN = {"method": "maxlfq", "formula": "~ 0 + groups", "numeric_covariates": [],
          "block_column": None, "block_effect": None, "block_scope": None,
          "contrasts": ["B-A"], "flags": []}


def de(br):
    prov = {"input": "/x/report.parquet", "metadata": "/x/conditions.csv",
            "input_sha256": "r" * 64, "metadata_sha256": "c" * 64, "design_request": DESIGN,
            "method": "maxlfq", "quantities": et.NORMALISED, "between_run": br}
    return {"prov": prov, "tables": {}}


class TheGateComparesTheStepAutoApplied(unittest.TestCase):
    def check(self, applied=None):
        cd = dict(DESIGN, status="recorded", between_run="auto")
        if applied:
            cd["between_run_applied"] = applied
        return {"check_design": cd}

    def auto(self):
        return dict(DESIGN, between_run="auto")

    def test_the_check_records_the_step_its_normalised_de_applied(self):
        raw = de({"method": "none", "requested": "auto", "why": "x"})
        raw["prov"]["quantities"] = et.RAW
        cd = nc.check_design_of(de({"method": "median_ratio", "requested": "auto", "why": "x"}),
                                raw, "r" * 64, None, "c" * 64)
        self.assertEqual((cd["between_run"], cd["between_run_applied"]), ("auto", "median_ratio"))

    def test_a_new_check_is_pending_then_held_to_its_step(self):
        dc = nc.design_check(self.check("median_ratio"), self.auto(), et.NORMALISED)
        self.assertIsNone(dc["refusal"])
        self.assertEqual(dc["between_run_pending"], "median_ratio")
        self.assertIn(nc.AUTO_STEP_PENDING, dc["match"])
        ok = nc.design_check(self.check("median_ratio"), self.auto(), et.NORMALISED,
                             applied="median_ratio")
        self.assertEqual((ok["refusal"], ok["match"]), (None, "matches"))
        bad = nc.design_check(self.check("median_ratio"), self.auto(), et.NORMALISED,
                              applied="none")
        self.assertIn("gave median_ratio", bad["refusal"])
        self.assertIn(nc.RERUN_CHECK, bad["refusal"])
        self.assertNotIn("before 2.15.1", bad["refusal"])

    def test_a_check_from_before_2_15_1_on_a_no_norm_dia_nn_report_is_refused(self):
        """It drew its volcano with quantile (the old auto); this DE's auto gives median_ratio."""
        dc = nc.design_check(self.check(), self.auto(), et.NORMALISED)
        self.assertIsNone(dc["refusal"])
        self.assertTrue(dc["between_run_pending"])
        dc = nc.design_check(self.check(), self.auto(), et.NORMALISED, applied="median_ratio")
        self.assertIn("gave quantile", dc["refusal"])
        self.assertIn("before 2.15.1", dc["refusal"])
        self.assertIn("gives median_ratio", dc["refusal"])

    def test_a_check_from_before_2_15_1_on_dia_nns_normalised_report_still_matches(self):
        dc = nc.design_check(self.check(), self.auto(), et.NORMALISED, applied="none")
        self.assertEqual((dc["refusal"], dc["match"]), (None, "matches"))

    def test_a_check_from_before_2_15_1_on_an_adapted_report_is_refused(self):
        """2.14.2-2.15.0 quantile-normalised every adapted report; 2.15.1 gives none to an
        engine that normalised (FragPipe IonQuant) and median_ratio to Sage."""
        for applied in ("none", "median_ratio"):
            dc = nc.design_check(self.check(), self.auto(), et.NORMALISED, applied=applied,
                                 adapted=True)
            self.assertIn("gave quantile", dc["refusal"] or "", applied)

    def test_the_old_rule_is_one_definition(self):
        self.assertEqual(nc.auto_step_before_2151("none", False), "none")
        self.assertEqual(nc.auto_step_before_2151("median_ratio", False), "quantile")
        self.assertEqual(nc.auto_step_before_2151("none", True), "quantile")

    def test_raw_and_explicit_steps_are_unchanged(self):
        self.assertIsNone(nc.design_check(self.check(), self.auto(), et.RAW,
                                          applied="none")["refusal"])
        dc = nc.design_check(self.check("quantile"), dict(DESIGN, between_run="quantile"),
                             et.NORMALISED)
        self.assertIn("differs from the one", dc["refusal"])

    def test_the_gate_cli_takes_the_adapted_flag(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, nc.RECORD)
            with open(path, "w") as fh:
                json.dump(dict(self.check(), schema=nc.SCHEMA, tripped=False, trips=[],
                               decision={"quantities": et.NORMALISED, "by": "x", "reason": "r"}),
                          fh)
            dj = os.path.join(d, "design.json")
            with open(dj, "w") as fh:
                json.dump(self.auto(), fh)
            base = [sys.executable, os.path.join(SCRIPTS, "normalization_check.py"), "gate",
                    "--check", path, "--quantities", et.NORMALISED, "--design-json", dj,
                    "--between-run-applied", "none"]
            ok = json.loads(subprocess.run(base, capture_output=True, text=True).stdout)
            self.assertTrue(ok["ok"], ok["message"])
            ad = json.loads(subprocess.run(base + ["--adapted-report"], capture_output=True,
                                           text=True).stdout)
            self.assertFalse(ad["ok"])
            self.assertIn("gave quantile", ad["message"])


# ------------------------------------------------------------------------------------- F2
class TrippedNeverHidesNotAssessed(unittest.TestCase):
    def rec(self):
        from test_normalization_not_assessed_2151 import rec
        r = rec(checks_ok=False, anomalies_norm=2)
        r.update(recommendation="x", decision=None, needs_decision=True, report="r",
                 conditions="c", groups={"A": 3, "B": 3}, design={}, bait=None, thresholds={},
                 experiment_type={"type": "whole_proteome", "label": "Whole", "source": "user"})
        self.assertTrue(r["tripped"] and r["not_assessed"])
        return r

    def test_the_summary_shows_both(self):
        md = nc.summary_md(self.rec())
        line = next(ln for ln in md.splitlines() if ln.startswith("- **Data check:**"))
        self.assertIn("TRIPPED", line)
        self.assertIn("NOT FULLY ASSESSED", line)
        self.assertIn("  - not assessed: B1", md)
        self.assertIn("  - tripped: ", md)

    def test_the_gate_message_says_both(self):
        r = self.rec()
        r.update(schema=nc.SCHEMA, decision={"quantities": et.NORMALISED, "by": "x",
                                             "reason": "r"})
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, nc.RECORD)
            with open(path, "w") as fh:
                json.dump(r, fh)
            ok, msg, _ = nc.gate(path, et.NORMALISED)
        self.assertTrue(ok, msg)
        self.assertIn("tripped: ", msg)
        self.assertIn(" and not fully assessed: ", msg)


# ------------------------------------------------------------------------------------- F4
class SageGroupsWithNoPeptideAreListed(unittest.TestCase):
    def test_the_csv(self):
        import run_search
        # three proteins: P2's only peptides are razor peptides of the larger P1 group
        _assigned, grouping = run_search.sage_protein_groups({
            "PEPA": "P1", "PEPB": "P1;P2", "PEPC": "P1;P2;P3", "PEPD": "P3", "PEPE": "P2;P3",
            "PEPF": "P1"})
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(run_search.write_sage_not_quantified(d, grouping),
                             run_search.SAGE_NOT_QUANTIFIED_CSV)
            with open(os.path.join(d, run_search.SAGE_NOT_QUANTIFIED_CSV)) as fh:
                rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), grouping["groups_not_quantified"])
        self.assertEqual([r["protein_group"] for r in rows],
                         [g["group"] for g in grouping["groups_not_quantified_list"]])
        for r in rows:
            self.assertIn("razor", r["why"])

    def test_none_is_a_header(self):
        import run_search
        with tempfile.TemporaryDirectory() as d:
            run_search.write_sage_not_quantified(d, {"groups_not_quantified_list": []})
            with open(os.path.join(d, run_search.SAGE_NOT_QUANTIFIED_CSV)) as fh:
                self.assertEqual(fh.read().strip(), ",".join(run_search.SAGE_NOT_QUANTIFIED_COLUMNS))

    def test_the_report_describes_it(self):
        import make_report
        with open(make_report.__file__, encoding="utf-8") as fh:
            src = fh.read()
        self.assertIn("sage_groups_not_quantified", src)


# ------------------------------------------------------------------------------------- F5, F6
FAKE_DOTNET = r'''#!/bin/bash
[ "$1" = --list-runtimes ] || exit 0
echo "Microsoft.NETCore.App 8.0.28 [x]"
echo "Microsoft.AspNetCore.App 8.0.28 [x]"
'''


@unittest.skipUnless(BASH, "needs bash")
class DotnetFolderAndInstallerTmp(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._t.name)
        self.addCleanup(self._t.cleanup)

    def run_ensure(self, *args, **env):
        e = no_cluster_env(self.d, base={"PATH": "/usr/bin:/bin", "HOME": self.d,
                                         "PROTEOMICS_ON_COMPUTE": "1"})
        e.update(env)
        return subprocess.run([BASH, ENSURE, *args], capture_output=True, text=True, env=e,
                              timeout=60)

    def test_the_default_folder_is_the_pipeline_homes(self):              # F5
        root = os.path.join(self.d, "pp", "dotnet8")
        _exe(os.path.join(root, "dotnet"), FAKE_DOTNET)
        p = self.run_ensure(PROTEOMICS_PIPELINE_HOME=os.path.join(self.d, "pp"))
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(p.stdout.splitlines()[-1], root)

    def test_a_root_at_the_old_default_is_reused(self):                    # F5
        old = os.path.join(self.d, ".proteomics-pipeline", "dotnet8")
        _exe(os.path.join(old, "dotnet"), FAKE_DOTNET)
        p = self.run_ensure(PROTEOMICS_PIPELINE_HOME=os.path.join(self.d, "pp"))
        self.assertEqual(p.stdout.splitlines()[-1], old, p.stderr)
        # ... but PROTEOMICS_DOTNET_DIR, when set, is the only place looked at
        p = self.run_ensure("--print-install-tmp",
                            PROTEOMICS_DOTNET_DIR=os.path.join(self.d, "mine", "dotnet8"))
        self.assertEqual(p.returncode, 0, p.stderr)

    def test_detect_acquisition_reads_the_same_default(self):              # F5
        import detect_acquisition as da
        old = os.environ.get(da.PIPELINE_HOME_ENV)
        try:
            os.environ[da.PIPELINE_HOME_ENV] = "/x/pp"
            self.assertEqual(da.dotnet_dir_default(), "/x/pp/dotnet8")
            os.environ.pop(da.PIPELINE_HOME_ENV)
            self.assertEqual(da.dotnet_dir_default(), os.path.expanduser(da.DOTNET_DIR_DEFAULT))
        finally:
            if old is not None:
                os.environ[da.PIPELINE_HOME_ENV] = old

    def test_a_version_like_tmpdir_is_not_given_to_the_installer(self):     # F6
        bad = os.path.join(self.d, "skill_test", "2.15.1", "jobtmp")
        os.makedirs(bad)
        # the bug, reproduced on dotnet-install.sh's own regex: a file with no version folder
        # of its own (the `dotnet` host) matches on TMPDIR's "2.15.1/"
        self.assertTrue(DOTNET_INSTALL_VERSION_RE.match(f"{bad}/dotnet.ABCDEFGHI/dotnet"))
        dest = os.path.join(self.d, "pp", "dotnet8")
        os.makedirs(os.path.dirname(dest))
        p = self.run_ensure("--print-install-tmp", TMPDIR=bad, PROTEOMICS_DOTNET_DIR=dest)
        self.assertEqual(p.returncode, 0, p.stderr)
        chosen = p.stdout.strip()
        self.assertNotIn("2.15.1", chosen)
        self.assertTrue(os.path.isdir(chosen))
        self.assertIsNone(DOTNET_INSTALL_VERSION_RE.match(f"{chosen}/ppdotnet.abc123/dotnet.ABCDEFGHI/dotnet"))
        # a plain TMPDIR is used as it is
        plain = os.path.join(self.d, "plain")
        os.makedirs(plain)
        self.assertEqual(self.run_ensure("--print-install-tmp", TMPDIR=plain,
                                         PROTEOMICS_DOTNET_DIR=dest).stdout.strip(), plain)

    def test_the_installer_runs_with_that_tmpdir(self):                     # F6
        with open(ENSURE, encoding="utf-8") as fh:
            src = fh.read()
        self.assertIn('parent="$(installer_tmp_parent)"', src)
        self.assertIn('TMPDIR="$work" bash "$tmp" --channel 8.0', src)


# ------------------------------------------------------------------------------------- F7, F9
@unittest.skipUnless(BASH, "needs bash")
class SetupJsonAgreesWithItself(unittest.TestCase):
    """setup.sh (not --check) through test_trfp_dotnet_discovery's stand-ins, on a pretend
    Linux (a stub `uname`)."""

    def setUp(self):
        import test_trfp_dotnet_discovery as td
        self.td = td

        class H(td.SetupInstallsTheParserSeparately):
            def runTest(self):
                pass
        self.h = H()
        self.h.setUp()
        self.addCleanup(self.h.tearDown)
        os.remove(os.path.join(self.h.sys, "uname"))
        _exe(os.path.join(self.h.sys, "uname"),
             '#!/bin/sh\ncase "$1" in -m) echo x86_64 ;; *) echo Linux ;; esac\n')

    def test_dia_raw_is_true_with_a_root_and_false_without(self):          # F7
        s, _ = self.h.setup()
        self.assertIs(s["ready_for"]["dia"], True)
        self.assertIs(s["ready_for"]["dia_raw"], True)
        self.assertNotIn("dia_raw", s["next_steps"])
        s, _ = self.h.setup(FAKE_ENSURE_FAIL="1")
        self.assertIs(s["ready_for"]["dia"], True, ".d and mzML need no .NET")
        self.assertIs(s["ready_for"]["dia_raw"], False)
        self.assertEqual(s["dotnet8"]["root"], "")
        self.assertIn("dia_raw", s["next_steps"])
        note = next(n for n in s["notes"] if "could not provide .NET 8" in n)
        self.assertIn("DIA-NN cannot read .raw", note)
        # the resolution reader is said to fail only when its own record says so
        rr = s["thermo_raw_reader"].get("resolution_reader") or {}
        self.assertEqual("Orbitrap resolution cannot be read" in note, rr.get("ready") is not True)

    def test_conda_runs_with_an_absolute_cache_from_the_prefix(self):      # F9
        log = os.path.join(self.h.d, "cwd.log")
        wrap = _exe(os.path.join(self.h.d, "wrapbin", "micromamba"),
                    f'#!/bin/bash\necho "$PWD|${{CONDA_PKGS_DIRS:-}}" >> "{log}"\n'
                    f'exec "{os.path.join(self.h.sys, "micromamba")}" "$@"\n')
        self.h.setup(PATH=os.path.dirname(wrap) + os.pathsep + self.h.sys,
                     CONDA_PKGS_DIRS="relative/pkgs")
        with open(log) as fh:
            calls = [ln.strip().split("|") for ln in fh if ln.strip()]
        self.assertTrue(calls)
        root = os.path.realpath(os.path.join(self.h.home, "micromamba"))
        for cwd, pkgs in calls:
            self.assertEqual(os.path.realpath(cwd), root)
            self.assertEqual(pkgs, os.path.join(self.h.home, "micromamba", "pkgs"))

    def test_absolute_pkgs_dirs_the_user_set_are_kept(self):               # F9
        log = os.path.join(self.h.d, "cwd.log")
        wrap = _exe(os.path.join(self.h.d, "wrapbin", "micromamba"),
                    f'#!/bin/bash\necho "${{CONDA_PKGS_DIRS:-}}" >> "{log}"\n'
                    f'exec "{os.path.join(self.h.sys, "micromamba")}" "$@"\n')
        mine = os.path.join(self.h.d, "mypkgs")
        self.h.setup(PATH=os.path.dirname(wrap) + os.pathsep + self.h.sys, CONDA_PKGS_DIRS=mine)
        with open(log) as fh:
            self.assertEqual({ln.strip() for ln in fh if ln.strip()}, {mine})


# ------------------------------------------------------------------------------------- F8
@unittest.skipUnless(BASH, "needs bash")
class AnUnpinnedCondaSageRecordsItsPackage(unittest.TestCase):
    def test_package_version_not_env_and_not_its_self_report(self):
        import test_acquire_tools_versions as at
        h = at.SageTests("test_a_sage_on_PATH_is_not_a_release")
        h.setUp()
        self.addCleanup(h.tearDown)
        # the harness's bin folder as a conda env's bin: <prefix>/bin/sage + <prefix>/conda-meta
        _exe(os.path.join(h.bin, "sage"), "#!/bin/sh\necho 'sage 0.14.6'\n")
        meta = os.path.join(os.path.dirname(h.bin), "conda-meta")
        os.makedirs(meta)
        with open(os.path.join(meta, "sage-proteomics-0.14.7-h0000000_2.json"), "w") as fh:
            fh.write("{}")
        t = h.acquire("linux")
        self.assertEqual(t["versions"]["sage"], "0.14.7")
        self.assertTrue(any("0.14.7 (conda package sage-proteomics)" in n for n in t["notes"]),
                        t["notes"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
