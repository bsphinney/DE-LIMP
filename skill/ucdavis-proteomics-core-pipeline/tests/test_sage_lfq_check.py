#!/usr/bin/env python3
"""
Sage's LFQ integrates MS1 only within +/- quant.lfq_settings.ppm_tolerance (default 5 ppm) of the
theoretical mass. gabrig 2026-09-29 (HeLa, Fusion Lumos, Sage 0.14.7): runs whose MS1 sat +7.1 /
+7.7 ppm off searched normally (29,846 PSMs) but LFQ logged "discovered 0 target MS1 peaks at 5%
FDR", and nothing warned. sage_lfq_check.py now reads each run's median precursor mass error and
Sage's MS1-peak count after every Sage search and WARNS; the record reaches search_provenance.json,
watch_run.sh --out, checkpoint.py status and audit_results.py --search-out (-> AUDIT.md -> the
report's "Audit & caveats").

The fixtures are Sage 0.14.7's own output shapes (crates/sage-cloudpath/src/parquet.rs: the
results.sage.parquet and lfq.parquet columns; precursor_ppm is ABSOLUTE, scoring.rs .abs()), with
per-run mass offsets set by hand: +7 ppm must warn, +/-1 ppm must not.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import sage_lfq_check as slc  # noqa: E402
from job_env import job_env  # noqa: E402

WATCH = os.path.join(SCRIPTS, "watch_run.sh")
AUDIT = os.path.join(SCRIPTS, "audit_results.py")
CHECKPOINT = os.path.join(SCRIPTS, "checkpoint.py")
RUN_SEARCH = os.path.join(SCRIPTS, "run_search.py")

try:
    import pyarrow  # noqa: F401
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False


def write_sage_outputs(out, offsets, peaks=None, lfq_q=0.01, lfq_settings=None, n=200,
                       spread=0.8, log=True):
    """Sage 0.14.7-shaped output in `out`: results.sage.parquet (n confident target PSMs per run
    at `offsets[run]` ppm, +/- `spread`, plus decoys at a wild error that must be ignored),
    lfq.parquet (targets and their +11.06 Da decoys, all at q_value `lfq_q`), results.json (the
    parameters as run; lfq_settings only if given, as Sage writes what the config set) and, if
    `log`, sage.log with Sage's "discovered N target MS1 peaks at 5% FDR"."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    os.makedirs(out, exist_ok=True)
    cols = {k: [] for k in ("filename", "is_decoy", "rank", "peptide", "proteins", "peptide_q",
                            "precursor_ppm", "expmass", "calcmass", "isotope_error")}

    def psm(run, i, ppm, decoy=False, q=0.001):
        calc = 1000.0 + 7.3 * i
        cols["filename"].append(run)
        cols["is_decoy"].append(decoy)
        cols["rank"].append(1)
        cols["peptide"].append(f"PEPTIDE{i}K")
        cols["proteins"].append(f"sp|P{i:05d}|PROT{i}_HUMAN")
        cols["peptide_q"].append(q)
        cols["precursor_ppm"].append(abs(ppm))                     # Sage writes it absolute
        cols["expmass"].append(calc * (1 + ppm * 1e-6))
        cols["calcmass"].append(calc)
        cols["isotope_error"].append(0.0)

    for run, off in offsets.items():
        for i in range(n):
            psm(run, i, off + spread * ((i % 5) - 2) / 2.0)
        for i in range(n // 4):
            psm(run, i, 40.0, decoy=True)                           # never part of the median
            psm(run, i, 30.0, q=0.5)                                # not confident: ignored too
    f32 = pa.float32()
    pq.write_table(pa.table({
        "filename": cols["filename"], "is_decoy": cols["is_decoy"],
        "rank": pa.array(cols["rank"], pa.int32()), "peptide": cols["peptide"],
        "proteins": cols["proteins"], "peptide_q": pa.array(cols["peptide_q"], f32),
        "precursor_ppm": pa.array(cols["precursor_ppm"], f32),
        "expmass": pa.array(cols["expmass"], pa.float64()),
        "calcmass": pa.array(cols["calcmass"], pa.float64()),
        "isotope_error": pa.array(cols["isotope_error"], f32)}),
        os.path.join(out, "results.sage.parquet"))
    lf = {k: [] for k in ("peptide", "stripped_peptide", "charge", "proteins", "is_decoy",
                          "q_value", "filename", "intensity")}
    for i in range(n):
        for decoy in (False, True):
            for run in offsets:
                lf["peptide"].append(f"PEPTIDE{i}K")
                lf["stripped_peptide"].append(f"PEPTIDE{i}K")
                lf["charge"].append(None)
                lf["proteins"].append(f"sp|P{i:05d}|PROT{i}_HUMAN")
                lf["is_decoy"].append(decoy)
                lf["q_value"].append(lfq_q)
                lf["filename"].append(run)
                lf["intensity"].append(1.85e5)
    pq.write_table(pa.table({**{k: v for k, v in lf.items() if k != "charge"},
                             "charge": pa.array(lf["charge"], pa.int32()),
                             "q_value": pa.array(lf["q_value"], f32),
                             "intensity": pa.array(lf["intensity"], f32)}),
                   os.path.join(out, "lfq.parquet"))
    quant = {"tmt": None, "lfq": True}
    if lfq_settings is not None:
        quant["lfq_settings"] = lfq_settings
    with open(os.path.join(out, "results.json"), "w") as fh:
        json.dump({"quant": quant, "mzml_paths": list(offsets)}, fh)
    if log:
        with open(os.path.join(out, "sage.log"), "w") as fh:
            fh.write("[2026-09-29T11:59:00Z INFO  sage] tracing MS1 features\n")
            if peaks is not None:
                fh.write(f"[2026-09-29T11:59:01Z INFO  sage] discovered {peaks} target MS1 "
                         f"peaks at 5% FDR\n")
            fh.write(f"[2026-09-29T11:59:01Z INFO  sage] discovered {len(offsets) * n} target "
                     f"peptide-spectrum matches at 1% FDR\n")


@unittest.skipUnless(HAVE_ARROW, "pyarrow is needed to write the Sage fixtures")
class LfqWindowCheck(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, self.d, True)
        self.out = os.path.join(self.d, "search")

    def test_plus_seven_ppm_warns_and_names_the_fix(self):
        """gabrig's cohort: +7.1/+7.7 ppm against Sage's default 5 ppm window, 0 peaks kept."""
        write_sage_outputs(self.out, {"HeLa_1.mzML": 7.1, "HeLa_2.mzML": 7.7}, peaks=0,
                           lfq_q=0.128)
        rec = slc.check(self.out)
        self.assertEqual(rec["status"], "warn", rec)
        self.assertEqual(rec["runs_outside_lfq_window"], ["HeLa_1.mzML", "HeLa_2.mzML"])
        self.assertEqual(rec["target_ms1_peaks_5pct_fdr"], 0)
        self.assertTrue(rec["low_ms1_peaks"])
        self.assertEqual(rec["ppm_tolerance"], 5.0)
        # the default is TAGGED as a default (architectural rule 2), not presented as set
        self.assertIn("DEFAULT", rec["ppm_tolerance_source"])
        self.assertEqual(rec["suggested_ppm_tolerance"], 10)          # ceil(7.7 + 2)
        self.assertAlmostEqual(rec["runs"]["HeLa_2.mzML"]["median_offset_ppm"], 7.7, delta=0.1)
        msg = rec["message"]
        self.assertIn('"lfq_settings": {"ppm_tolerance": 10}', msg)
        self.assertIn("0 target MS1 peaks at 5% FDR", msg)
        self.assertIn("identifications are not affected, only the MS1 quantities", msg)
        self.assertIn("Nothing was re-run", msg)

    def test_plus_minus_one_ppm_does_not_warn(self):
        write_sage_outputs(self.out, {"a.mzML": 1.0, "b.mzML": -1.0}, peaks=380)
        rec = slc.check(self.out)
        self.assertEqual(rec["status"], "ok", rec)
        self.assertEqual(rec["runs_outside_lfq_window"], [])
        self.assertFalse(rec["low_ms1_peaks"])
        # a negative offset is kept signed for direction; the rule uses Sage's absolute error
        self.assertLess(rec["runs"]["b.mzML"]["median_offset_ppm"], 0)
        self.assertGreater(rec["runs"]["b.mzML"]["median_abs_ppm"], 0)

    def test_one_run_off_is_named_alone(self):
        write_sage_outputs(self.out, {"good.mzML": 0.5, "off.mzML": 3.6}, peaks=390)
        rec = slc.check(self.out)
        self.assertEqual(rec["status"], "warn")
        self.assertEqual(rec["runs_outside_lfq_window"], ["off.mzML"])   # 3.6 + 2 > 5
        self.assertIn("1 of 2 run(s)", rec["message"])

    def test_a_wider_window_that_sage_ran_with_clears_it(self):
        """results.json is what Sage actually ran: a 10 ppm window fits +7 ppm."""
        write_sage_outputs(self.out, {"a.mzML": 7.0}, peaks=150,
                           lfq_settings={"ppm_tolerance": 10.0})
        rec = slc.check(self.out)
        self.assertEqual(rec["status"], "ok", rec)
        self.assertEqual(rec["ppm_tolerance"], 10.0)
        self.assertIn("results.json", rec["ppm_tolerance_source"])
        self.assertNotIn("DEFAULT", rec["ppm_tolerance_source"])

    def test_zero_peaks_warns_even_when_the_mass_error_fits(self):
        write_sage_outputs(self.out, {"a.mzML": 0.2}, peaks=0)
        rec = slc.check(self.out)
        self.assertEqual(rec["status"], "warn")
        self.assertEqual(rec["runs_outside_lfq_window"], [])
        self.assertIn("cannot name the cause", rec["message"])
        self.assertNotIn("ppm_tolerance\": ", rec["message"])       # no window advice here

    def test_peak_count_falls_back_to_lfq_parquet_without_a_log(self):
        write_sage_outputs(self.out, {"a.mzML": 7.0}, lfq_q=0.128, log=False)
        rec = slc.check(self.out)
        self.assertEqual(rec["target_ms1_peaks_5pct_fdr"], 0)
        self.assertIn("lfq.parquet", rec["target_ms1_peaks_source"])
        self.assertIn("lower bound", rec["target_ms1_peaks_source"])   # the stored q is shared
        write_sage_outputs(self.out, {"a.mzML": 0.1}, lfq_q=0.01, log=False)
        rec = slc.check(self.out)
        # targets only: the +11.06 Da decoy rows are not peaks
        self.assertEqual(rec["target_ms1_peaks_5pct_fdr"], 200)
        self.assertEqual(rec["status"], "ok")

    def test_no_results_is_unchecked_and_says_so(self):
        os.makedirs(self.out)
        rec = slc.check(self.out)
        self.assertEqual(rec["status"], "unchecked")
        self.assertIn("results.sage.parquet", rec["message"])
        self.assertIn("MS1 peaks at 5% FDR", rec["message"])

    def test_too_few_psms_to_judge_is_said_not_passed(self):
        write_sage_outputs(self.out, {"tiny.mzML": 9.0}, peaks=30, n=20)
        rec = slc.check(self.out)
        self.assertEqual(rec["status"], "unchecked", rec)          # never an "ok" on nothing
        self.assertIn("50 confident target PSMs", rec["message"])
        write_sage_outputs(self.out, {"good.mzML": 0.5}, peaks=180)
        rec = slc.check(self.out)
        self.assertEqual(rec["status"], "ok")
        self.assertNotIn("not judged", rec["message"])

    def test_lfq_off_is_not_applicable(self):
        write_sage_outputs(self.out, {"a.mzML": 9.0}, peaks=0)
        with open(os.path.join(self.out, "results.json"), "w") as fh:
            json.dump({"quant": {"lfq": False}}, fh)
        self.assertEqual(slc.check(self.out)["status"], "not_applicable")

    def test_record_goes_into_search_provenance(self):
        write_sage_outputs(self.out, {"a.mzML": 7.1}, peaks=0)
        with open(os.path.join(self.out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "sage", "version": "0.14.7"}, fh)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sage_lfq_check.py"),
                            "--out", self.out], capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)             # a warning is not a failure
        self.assertIn("[sage_lfq_check] WARNING:", r.stderr)
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["engine"], "sage")                # merged, not replaced
        self.assertEqual(prov["sage_lfq_check"]["status"], "warn")
        with open(os.path.join(self.out, "sage_lfq_check.json")) as fh:
            self.assertEqual(json.load(fh)["status"], "warn")


@unittest.skipUnless(HAVE_ARROW, "pyarrow is needed to write the Sage fixtures")
class WarningReachesEveryReader(unittest.TestCase):
    """The same record, read by the watcher, the recovery status and the audit (-> report)."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(__import__("shutil").rmtree, self.d, True)
        self.out = os.path.join(self.d, "search")

    def _record(self, offsets, peaks):
        write_sage_outputs(self.out, offsets, peaks=peaks)
        rec = slc.check(self.out)
        slc.record(self.out, rec)
        return rec

    def test_watcher_reports_the_warning_and_the_peak_count(self):
        self._record({"a.mzML": 7.3}, 0)
        r = subprocess.run(["bash", WATCH, "--log", os.path.join(self.out, "sage.log"),
                            "--out", self.out], capture_output=True, text=True, timeout=120,
                           env=job_env(self.d))
        res = json.loads(r.stdout)
        self.assertEqual(res["sage_lfq"]["status"], "warn")
        self.assertEqual(res["sage_lfq"]["target_ms1_peaks_5pct_fdr"], 0)
        self.assertTrue(any("label-free quantification is unreliable" in w
                            for w in res["warnings"]), res)
        self.assertFalse(res["failed"])                         # a warning, not a failure
        self.assertEqual(res["error_class"], "")

    def test_watcher_without_out_reads_the_log(self):
        write_sage_outputs(self.out, {"a.mzML": 0.3}, peaks=4242)
        r = subprocess.run(["bash", WATCH, "--log", os.path.join(self.out, "sage.log")],
                           capture_output=True, text=True, timeout=120, env=job_env(self.d))
        res = json.loads(r.stdout)
        self.assertEqual(res["sage_lfq"]["target_ms1_peaks_5pct_fdr"], 4242)
        self.assertNotIn("warnings", res)

    def test_ok_record_adds_no_watcher_warning(self):
        self._record({"a.mzML": 0.4}, 180)
        r = subprocess.run(["bash", WATCH, "--out", self.out], capture_output=True, text=True,
                           timeout=120, env=job_env(self.d))
        res = json.loads(r.stdout)
        self.assertEqual(res["sage_lfq"]["status"], "ok")
        self.assertNotIn("warnings", res)

    def test_audit_and_so_the_report_carry_it(self):
        self._record({"a.mzML": 7.3, "b.mzML": 7.9}, 0)
        audit = os.path.join(self.d, "AUDIT.md")
        r = subprocess.run([sys.executable, AUDIT, "--out", audit, "--search-out", self.out],
                           capture_output=True, text=True, timeout=120, cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        text = open(audit).read()
        self.assertIn("**sage_lfq**", text)
        self.assertIn("⚠️ **sage_lfq**", text)
        self.assertIn("label-free quantification is unreliable", text)
        with open(os.path.join(self.d, "AUDIT.json")) as fh:
            self.assertEqual(json.load(fh)["overall"], "WARN")

    def test_audit_passes_a_fitting_window(self):
        self._record({"a.mzML": 1.0}, 190)
        audit = os.path.join(self.d, "AUDIT.md")
        subprocess.run([sys.executable, AUDIT, "--out", audit, "--search-out", self.out],
                       capture_output=True, text=True, timeout=120, cwd=self.d, check=True)
        self.assertIn("✅ **sage_lfq**", open(audit).read())

    def test_recovery_status_says_it(self):
        self._record({"a.mzML": 7.3}, 0)
        session = os.path.join(self.d, "session")
        cfg = os.path.join(self.d, "claude_config")
        os.makedirs(cfg)
        env = job_env(self.d, CLAUDE_CONFIG_DIR=cfg)
        subprocess.run([sys.executable, CHECKPOINT, "record", "--session", session, "--stage",
                        "search", "--report", os.path.join(self.out, "report.parquet")],
                       capture_output=True, text=True, timeout=60, env=env, check=True)
        r = subprocess.run([sys.executable, CHECKPOINT, "status", "--session", session],
                           capture_output=True, text=True, timeout=60, env=env)
        res = json.loads(r.stdout)
        self.assertTrue(any("Sage LFQ:" in w for w in res.get("warnings", [])), res)

    def test_adapt_only_rechecks_prints_and_records(self):
        """After a Sage job, `run_search.py --adapt-only` is the step the orchestrator runs: the
        warning is on its screen, and in search_provenance.json, before any DE."""
        self._record({"a.mzML": 7.3}, 0)
        os.remove(os.path.join(self.out, "sage_lfq_check.json"))
        with open(os.path.join(self.out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "sage"}, fh)
        tools, bundle, cfg, fasta = (os.path.join(self.d, n) for n in
                                     ("tools.json", "bundle.json", "sage.json", "db.fasta"))
        for p, v in ((tools, {"sage": "/opt/sage/sage"}),
                     (bundle, {"acquisition": "DDA", "engine": {"name": "sage"}}),
                     (cfg, {"quant": {"lfq": True}})):
            with open(p, "w") as fh:
                json.dump(v, fh)
        open(fasta, "w").close()
        argv = [sys.executable, RUN_SEARCH, "--tools", tools, "--bundle", bundle,
                "--params", cfg, "--fasta", fasta, "--out", self.out,
                "--files", "a.mzML", "--engine", "sage", "--adapt-only"]
        env = job_env(self.d, PATH="/usr/bin:/bin")
        r = subprocess.run(argv, capture_output=True, text=True, timeout=120, env=env)
        # since 2.10 a window that does not fit the mass error is a gate, not only a warning:
        # tests/test_sage_lfq_gate.py has the rest
        self.assertNotEqual(r.returncode, 0, r.stderr)
        self.assertIn("[sage_lfq_check] WARNING:", r.stderr)
        self.assertIn("REFUSED: report.parquet is not built from this Sage search", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "report.parquet")))
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            self.assertEqual(json.load(fh)["sage_lfq_check"]["status"], "warn")
        r = subprocess.run(argv + ["--accept-lfq-window", "core analyst: IDs only, quantities "
                                   "not reported"], capture_output=True, text=True, timeout=120,
                           env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        res = json.loads(r.stdout[r.stdout.index("{"):])        # after adapt_sage's own line
        self.assertEqual(res["sage_lfq_check"]["status"], "warn")
        self.assertIn("IDs only", res["sage_lfq_check"]["accepted"]["reason"])
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            self.assertEqual(json.load(fh)["sage_lfq_check"]["status"], "warn")


if __name__ == "__main__":
    unittest.main()
