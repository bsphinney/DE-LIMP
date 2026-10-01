#!/usr/bin/env python3
"""
DIA-NN on DDA data, through every route -- the four defects msalemi hit on SET28 (PROT_0000, 28
Exploris 480 DDA .raw, 60k/15k, skill 2.6.0, 2026-09-25):

  (a) no --dda in the 5-step chain. run_search.py appended it to the SINGLE-SHOT command only;
      estimate_params.py never wrote it and diann_parallel.py strips nothing -- so a DDA cohort of
      more than 5 files was searched as DIA, with no error. Caught by a manual grep.
  (b) the precursor range fell back to 380-980. detect_acquisition.py reported no range for DDA
      and SKILL.md said to omit it; her survey scans were 350-1500 and her FragPipe PSMs spanned
      360-1315 m/z, so about a third would have been cut from the predicted library.
  (c) step 1b's probe cannot measure under --dda: DIA-NN logs no scan-window radius there, and
      each probe ran to its 3600 s timeout -- ~3 h, then a failed step 1b.
  (11) submit.sh lost jobs.txt and RECOVERY.md when its stdout was piped into `head -1`: under
      `set -euo pipefail` the second echo died of SIGPIPE, after every sbatch had succeeded.

Each fix is tested where it lives: the cfg (estimate_params.py), the chain's gate and scripts
(diann_parallel.py), both routes' refusal (run_search.py), the probe (probe_window.py), the range
reader (detect_acquisition.py, three formats) and submit.sh run with a closed stdout.
"""
import json
import os
import sqlite3
import struct
import subprocess
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

from job_env import job_env  # noqa: E402  (env for running generated scripts)
from synthetic_tdf import synthetic_tdf_write_uri  # noqa: E402  (fixture writes only)

import detect_acquisition as da  # noqa: E402
import diann_parallel as dp  # noqa: E402
import estimate_params as ep  # noqa: E402
import run_search  # noqa: E402
# Module imports, not `from ... import`: a TestCase class bound in this namespace would be
# collected and run a second time here.
import test_single_shot_mass_accuracy as ss  # noqa: E402  (fake DIA-NN + no-SLURM env)
import test_step1b_mass_accuracy_probe as fx  # noqa: E402  (fake-executable helper)
import test_thermo_raw_detection as trfp  # noqa: E402  (fake ThermoRawFileParser harness)

ESTIMATE = os.path.join(SCRIPTS, "estimate_params.py")
GENERATE = os.path.join(SCRIPTS, "diann_parallel.py")
PROBE = os.path.join(SCRIPTS, "probe_window.py")

# SET28: Orbitrap Exploris 480, MS1 60,000 / MS2 15,000, survey scans 350-1500 m/z
SET28 = ("--instrument", "Orbitrap Exploris 480", "--ms1-resolution", "60000",
        "--ms2-resolution", "15000", "--resolution-source", "detected")
SURVEY = ("--precursor-mz-range", "350", "1500")

# sbatch that prints a new job id per call (and refuses a script it cannot open)
FAKE_SBATCH = r"""#!/bin/bash
for last; do :; done
[ -f "$last" ] || { echo "sbatch: error: Unable to open file $last" >&2; exit 1; }
n=$(( $(cat "$SB_COUNTER" 2>/dev/null || echo 700) + 1 )); echo "$n" > "$SB_COUNTER"
echo "$n"
"""


def estimate(d, acquisition, *args, name="params.cfg", engine="diann"):
    """-> (CompletedProcess, cfg path)."""
    out = os.path.join(d, name)
    p = subprocess.run([sys.executable, ESTIMATE, "--engine", engine, "--acquisition",
                        acquisition, *args, "--out", out], capture_output=True, text=True,
                       timeout=60)
    return p, out


def read(path):
    with open(path) as fh:
        return fh.read()


def cfg_flags(path):
    return [f for f, _ in dp.cfg_groups(dp.cfg_tokens(path))]


def side(path):
    with open(path + ".rationale.json") as fh:
        return json.load(fh)


def cohort(d, n=6, ext=".mzML"):
    """Input files the generator accepts without reading them (mzML: no .NET prefix)."""
    paths = []
    for i in range(n):
        p = os.path.join(d, f"run{i:02d}{ext}")
        with open(p, "w") as fh:
            fh.write("<mzML/>")
        paths.append(p)
    return paths


def fasta(d):
    p = os.path.join(d, "db.fasta")
    with open(p, "w") as fh:
        fh.write(">sp|P1|X\nPEPTIDER\n")
    return p


def generate(d, cfg, out, *more):
    return subprocess.run([sys.executable, GENERATE, "--diann", "/bin/true", "--raw",
                           *cohort(d), "--fasta", fasta(d), "--out", out, "--cfg", cfg, *more],
                          capture_output=True, text=True, timeout=120)


# ------------------------------------------------------------------------------------------
# (a)(b)(c) the cfg estimate_params.py writes
# ------------------------------------------------------------------------------------------
class DdaCfg(unittest.TestCase):
    def test_a_dda_cfg_carries_dda_and_a_dia_cfg_does_not(self):
        with tempfile.TemporaryDirectory() as d:
            p, dda = estimate(d, "DDA", *SET28, *SURVEY, name="dda.cfg")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("--dda", cfg_flags(dda))
            self.assertIn("must not be used with DIA data", side(dda)["rationale"]["--dda"]["source"])
            p, dia = estimate(d, "DIA", *SET28, *SURVEY, name="dia.cfg")
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertNotIn("--dda", cfg_flags(dia))
            self.assertNotIn("--dda", side(dia)["rationale"])

    def test_the_dda_range_is_the_survey_scan_and_says_so(self):
        with tempfile.TemporaryDirectory() as d:
            p, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            self.assertEqual(p.returncode, 0, p.stderr)
            text = read(cfg)
            self.assertIn("--min-pr-mz 350\n", text)
            self.assertIn("--max-pr-mz 1500\n", text)
            self.assertNotIn("--min-pr-mz 380", text)
            src = side(cfg)["rationale"]["--min-pr-mz"]["source"]
            self.assertIn("MS1 survey scan range", src)
            self.assertNotIn("isolation windows", src)

    def test_a_diann_dda_cfg_without_a_range_is_refused_and_nothing_is_left_behind(self):
        """The 380-980 fallback is gone for DDA. A cfg (and sidecar) from an earlier run of the
        same command must not survive the refusal: the next step reads them by path."""
        with tempfile.TemporaryDirectory() as d:
            p, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertTrue(os.path.exists(cfg + ".rationale.json"))
            p, cfg = estimate(d, "DDA", *SET28)                  # same --out, no range
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("--precursor-mz-range", p.stderr)
            self.assertIn("survey scan range", p.stderr)
            self.assertIn("No cfg was written", p.stderr)
            self.assertFalse(os.path.exists(cfg))
            self.assertFalse(os.path.exists(cfg + ".rationale.json"))

    def test_dia_without_a_range_still_falls_back_and_sage_dda_needs_none(self):
        with tempfile.TemporaryDirectory() as d:
            p, cfg = estimate(d, "DIA", *SET28)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("FALLBACK", side(cfg)["rationale"]["--min-pr-mz"]["source"])
            p, _ = estimate(d, "DDA", *SET28, engine="sage", name="sage.json")
            self.assertEqual(p.returncode, 0, p.stderr)

    def test_dda_pins_mass_accuracy_the_table_level_and_the_sop_for_the_other(self):
        """(c): no measure_with_diann for DDA. 60k MS1 -> DIA-NN's documented 10 ppm; the 15k
        MS2 (outside the table, which DIA would measure) -> the SOP 20 ppm, tagged DEFAULT."""
        with tempfile.TemporaryDirectory() as d:
            p, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            self.assertEqual(p.returncode, 0, p.stderr)
            text = read(cfg)
            self.assertIn(f"--mass-acc {ep.SOP_MASS_ACC['ms2_ppm']:g}\n", text)
            self.assertIn("--mass-acc-ms1 10\n", text)
            s = side(cfg)
            self.assertEqual(s["mass_accuracy_plan"], ep.PLAN_PINNED)
            self.assertEqual(s["mass_accuracy_default"], {"--mass-acc": 20})
            self.assertIn("DEFAULT, not user-confirmed", s["rationale"]["--mass-acc"]["source"])
            self.assertNotIn("DEFAULT", s["rationale"]["--mass-acc-ms1"]["source"])
            self.assertIn("DEFAULT, not user-confirmed",
                          s["rationale"]["mass_accuracy_plan"]["source"])
            # ...while the same instrument on DIA is still measured
            p, dia = estimate(d, "DIA", *SET28, *SURVEY, name="dia.cfg")
            self.assertEqual(side(dia)["mass_accuracy_plan"], ep.MEASURE_WITH_DIANN)
            self.assertEqual(side(dia)["mass_accuracy_default"], {})

    def test_an_overridden_level_is_not_called_a_default(self):
        with tempfile.TemporaryDirectory() as d:
            p, cfg = estimate(d, "DDA", *SET28, *SURVEY, "--overrides", '{"--mass-acc": 25}')
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("--mass-acc 25\n", read(cfg))
            self.assertEqual(side(cfg)["mass_accuracy_default"], {})

    def test_the_dda_window_rationale_says_it_is_not_measured(self):
        with tempfile.TemporaryDirectory() as d:
            p, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            self.assertEqual(side(cfg)["rationale"]["--window"]["source"], ep.DDA_WINDOW_NOTE)
            self.assertNotIn("--window", cfg_flags(cfg))


# ------------------------------------------------------------------------------------------
# (a)(c) the 5-step chain
# ------------------------------------------------------------------------------------------
class DdaChain(unittest.TestCase):
    def test_a_dda_cfg_is_parallel_safe_without_a_probe(self):
        with tempfile.TemporaryDirectory() as d:
            _, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            v = dp.parallel_safe(cfg)
            self.assertTrue(v["ok"], v["reason"])
            self.assertFalse(v["probe"])
            self.assertEqual(v["measure"], [])
            self.assertEqual(v["code"], "dda_window_unset")

    def test_every_step_carries_dda_and_there_is_no_step_1b(self):
        with tempfile.TemporaryDirectory() as d:
            _, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            out = os.path.join(d, "out")
            p = generate(d, cfg, out)
            self.assertEqual(p.returncode, 0, p.stderr)
            info = json.loads(p.stdout)
            self.assertFalse(os.path.exists(os.path.join(out, "step1b_window.sbatch")))
            self.assertEqual(info["step1b_measures"], [])
            for step in ("step1_libpred", "step2_firstpass", "step3_assembly",
                         "step4_finalpass", "step5_report"):
                body = read(os.path.join(out, step + ".sbatch"))
                runs = [ln for ln in body.splitlines() if "--threads" in ln]
                self.assertTrue(runs, step)
                self.assertTrue(all(" --dda" in ln for ln in runs), f"{step}: {runs}")
            self.assertNotIn("probe_window.py", read(os.path.join(out, "submit.sh")))
            # the provenance says what ran: the window unmeasured, MS2 a DEFAULT
            self.assertEqual(info["scan_window"]["source"], ep.DDA_WINDOW_NOTE)
            self.assertEqual(info["mass_acc"]["default"], {"--mass-acc": 20})
            self.assertIn("DEFAULT", info["mass_acc"]["default_note"])

    def test_a_pre_2_9_measure_plan_beside_dda_is_refused_not_probed(self):
        """A cfg planned `measure_with_diann` beside --dda -- what estimate_params.py wrote for
        DDA before 2.9 (acquisition DDA, the plan, --dda added by hand as the docs then said):
        step 1b would probe under --dda and burn its timeouts. Refused, loudly."""
        with tempfile.TemporaryDirectory() as d:
            _, cfg = estimate(d, "DIA", *SET28, *SURVEY)
            with open(cfg, "a") as fh:
                fh.write("--dda\n")
            sidecar = side(cfg)
            sidecar["acquisition"] = "DDA"
            with open(cfg + ".rationale.json", "w") as fh:
                json.dump(sidecar, fh)
            v = dp.parallel_safe(cfg)
            self.assertFalse(v["ok"])
            self.assertEqual(v["code"], "mass_acc_dda")
            self.assertIn("--acquisition DDA", v["remedy"])
            out = os.path.join(d, "out")
            p = generate(d, cfg, out)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("mass accuracy is to be measured", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out, "step1b_window.sbatch")))

    def test_the_generator_itself_refuses_a_cfg_that_disagrees_with_its_sidecar(self):
        """F3: the documented direct call of diann_parallel.py, and a --seed-lib phase 2, never
        pass through run_search.py's check. The generator reads the acquisition the cfg was
        estimated for from <cfg>.rationale.json and refuses a mismatch itself."""
        with tempfile.TemporaryDirectory() as d:
            _, dda = estimate(d, "DDA", *SET28, *SURVEY, name="dda.cfg")
            with open(dda) as fh:                   # --dda removed by hand
                text = fh.read().replace("--dda\n", "")
            with open(dda, "w") as fh:
                fh.write(text)
            _, dia = estimate(d, "DIA", "--instrument", "Orbitrap Exploris 480",
                              "--ms1-resolution", "120000", "--ms2-resolution", "30000",
                              *SURVEY, name="dia.cfg")
            with open(dia, "a") as fh:              # --dda added by hand to a DIA cfg
                fh.write("--dda\n")
            for cfg, words in ((dda, "acquisition is DDA, but"), (dia, "acquisition is DIA, but")):
                out = os.path.join(d, "out_" + os.path.basename(cfg))
                p = generate(d, cfg, out)
                self.assertNotEqual(p.returncode, 0, cfg)
                self.assertIn("REFUSING", p.stderr)
                self.assertIn(words, p.stderr)
                self.assertFalse(os.path.exists(os.path.join(out, "submit.sh")), cfg)
            # --seed-lib is covered the same way: the check runs before anything else
            seed = os.path.join(d, "seed.parquet")
            with open(seed, "w") as fh:
                fh.write("lib")
            p = generate(d, dda, os.path.join(d, "out_seed"), "--seed-lib", seed)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("acquisition is DDA, but", p.stderr)

    def test_window_zero_in_a_dda_cfg_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            _, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            with open(cfg, "a") as fh:
                fh.write("--window 0\n")
            v = dp.parallel_safe(cfg)
            self.assertEqual((v["ok"], v["code"]), (False, "window_invalid"))
            self.assertIn("DDA", v["remedy"])

    def test_dda_mismatch_both_ways(self):
        with tempfile.TemporaryDirectory() as d:
            _, dda = estimate(d, "DDA", *SET28, *SURVEY, name="dda.cfg")
            _, dia = estimate(d, "DIA", *SET28, *SURVEY, name="dia.cfg")
            self.assertIsNone(dp.dda_mismatch(dda, "DDA"))
            self.assertIsNone(dp.dda_mismatch(dia, "dia"))
            self.assertIn("no --dda", dp.dda_mismatch(dia, "DDA"))
            self.assertIn("must not be used with DIA data", dp.dda_mismatch(dda, "DIA"))
            for unknown in ("", None, "unknown", "mixed"):
                self.assertIsNone(dp.dda_mismatch(dia, unknown))


# ------------------------------------------------------------------------------------------
# (a) run_search.py refuses a mismatch on BOTH routes, and the single-shot search takes --dda
#     from the cfg (once)
# ------------------------------------------------------------------------------------------
class RunSearchRefusesMismatch(unittest.TestCase):
    # Borrow the single-shot harness (fake DIA-NN, no SLURM on PATH), not its tests.
    _setup = ss.SingleShotMassAccTests._setup
    _env = ss.SingleShotMassAccTests._env
    _run_search = ss.SingleShotMassAccTests._run_search
    def test_the_chain_route_refuses_a_dda_bundle_whose_cfg_lacks_dda(self):
        with tempfile.TemporaryDirectory() as d:
            _, dia = estimate(d, "DIA", *SET28, *SURVEY)
            out = os.path.join(d, "out")
            with self.assertRaises(SystemExit) as cm:
                run_search.run_diann_parallel("/bin/true", dia, cohort(d), fasta(d), out, 8,
                                              None, acquisition="DDA")
            self.assertIn("no --dda", str(cm.exception))
            self.assertFalse(os.path.exists(out), "nothing may be generated")

    def test_the_single_shot_route_refuses_a_dia_bundle_whose_cfg_has_dda(self):
        with tempfile.TemporaryDirectory() as d:
            _, dda = estimate(d, "DDA", *SET28, *SURVEY)
            out = os.path.join(d, "out")
            with self.assertRaises(SystemExit) as cm:
                run_search.run_diann("/bin/true", dda, cohort(d, 2), fasta(d), out, 8,
                                     os.path.join(d, "job.sh"), acquisition="DIA")
            self.assertIn("must not be used with DIA data", str(cm.exception))
            self.assertFalse(os.path.exists(out), "nothing may be generated")

    def test_a_single_shot_dda_cfg_with_a_leftover_measure_plan_is_refused_at_generation(self):
        """F4: it used to be carried into the job, where probe_window.py refused --dda at run
        time -- after the library job had run. The chain's own rule (mass_acc_dda_refusal)
        refuses it now, before anything is written."""
        with tempfile.TemporaryDirectory() as d:
            raws, fasta_, _, tools, bundle = self._setup(d)
            p, cfg = estimate(d, "DIA", *SET28, *SURVEY, name="plan.cfg")    # measure_with_diann
            self.assertEqual(side(cfg)["mass_accuracy_plan"], ep.MEASURE_WITH_DIANN)
            with open(cfg, "a") as fh:
                fh.write("--dda\n")
            with open(bundle, "w") as fh:
                json.dump({"acquisition": "DDA", "engine": {"name": "diann"}}, fh)
            self.assertEqual(dp.mass_acc_dda_refusal(cfg)["code"], "mass_acc_dda")
            self.assertEqual(dp.parallel_safe(cfg)["code"], "mass_acc_dda")
            p, out = self._run_search(d, raws, fasta_, cfg, tools, bundle,
                                      "--sbatch", os.path.join(d, "job.sh"))
            self.assertNotEqual(p.returncode, 0, p.stdout)
            self.assertIn("REFUSING the search", p.stderr)
            self.assertIn(dp.MASS_ACC_DDA_REASON, p.stderr)
            for f in ("job_1_lib.sh", "job_2_search.sh"):
                self.assertFalse(os.path.exists(os.path.join(d, f)), f)

    def test_the_single_shot_search_takes_dda_from_the_cfg_once(self):
        """It used to append --dda to the command; with --dda now in the cfg that would put it
        on the command line twice."""
        with tempfile.TemporaryDirectory() as d:
            raws, fasta_, _, tools, bundle = self._setup(d)
            p, cfg = estimate(d, "DDA", *SET28, *SURVEY, name="dda.cfg")
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(bundle, "w") as fh:
                json.dump({"acquisition": "DDA", "engine": {"name": "diann"}}, fh)
            p, out = self._run_search(d, raws, fasta_, cfg, tools, bundle,
                                      "--sbatch", os.path.join(d, "job.sh"))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            search = [ln for ln in read(os.path.join(d, "job_2_search.sh")).splitlines()
                      if "--reanalyse" in ln]
            self.assertEqual(len(search), 1, search)
            self.assertEqual(search[0].split().count("--dda"), 1, search[0])
            with open(os.path.join(out, "search_provenance.json")) as fh:
                prov = json.load(fh)
            self.assertIs(prov["result"]["dda"], True)
            self.assertEqual(prov["result"]["mass_acc"]["default"], {"--mass-acc": 20})
            self.assertEqual(prov["scan_window"]["source"], ep.DDA_WINDOW_NOTE)


# ------------------------------------------------------------------------------------------
# (c) rule 2 in published text: the Methods say which mass accuracy is the SOP default
# ------------------------------------------------------------------------------------------
import make_methods as mm  # noqa: E402

SOP_NOTE = ep.SOP_DEFAULT_PUBLISHED


class MethodsNameTheSopDefault(unittest.TestCase):
    """The 20 ppm a DDA cfg pins for a 15k MS2 is the Core SOP, not a value anyone chose for
    these data or measured on them. The rationale, the sidecar and search_provenance.json said
    so; the Methods paragraph -- the user-facing text -- printed a bare "fragment (MS2) 20 ppm".
    The tag is read from the run's provenance, not from a branch on DDA."""
    _setup = ss.SingleShotMassAccTests._setup
    _env = ss.SingleShotMassAccTests._env
    _run_search = ss.SingleShotMassAccTests._run_search

    def _search(self, d, acquisition, *res):
        raws, fasta_, _, tools, bundle = self._setup(d)
        p, cfg = estimate(d, acquisition, "--instrument", "Orbitrap Exploris 480", *res,
                          "--resolution-source", "detected", *SURVEY, name=f"{acquisition}.cfg")
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(bundle, "w") as fh:
            json.dump({"acquisition": acquisition, "engine": {"name": "diann"}}, fh)
        p, out = self._run_search(d, raws, fasta_, cfg, tools, bundle,
                                  "--sbatch", os.path.join(d, "job.sh"))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        return raws, cfg, os.path.join(out, "search_provenance.json")

    def test_the_methods_sentence_and_table_tag_the_sop_level_only(self):
        with tempfile.TemporaryDirectory() as d:
            raws, cfg, prov = self._search(d, "DDA", "--ms1-resolution", "60000",
                                           "--ms2-resolution", "15000")
            md = os.path.join(d, "methods.md")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                                "--raw", raws[0], "--params", cfg, "--search-prov", prov,
                                "--instrument", "Orbitrap Exploris 480", "--acquisition", "DDA",
                                "--out", md], capture_output=True, text=True, timeout=120,
                               env=dict(os.environ, THERMORAWFILEPARSER_SHARED="",
                                        THERMO_RESOLUTION_PYTHON=""))
            self.assertEqual(r.returncode, 0, r.stderr)
            text = read(md)
            self.assertIn(f"precursor (MS1) 10 ppm and fragment (MS2) 20 ppm ({SOP_NOTE})", text)
            self.assertNotIn(f"10 ppm ({SOP_NOTE})", text, "the table-derived MS1 is not a default")
            self.assertIn(f"| Fragment (MS2) tolerance | 20 ppm ({SOP_NOTE}) |", text)
            self.assertIn("| Precursor (MS1) tolerance | 10 ppm |", text)
            self.assertIn("DEFAULT per search_provenance.json result.mass_acc.default", text)
            # ...and the paragraph says DIA-NN searched the spectra as DDA
            self.assertIn("in its DDA mode (--dda", text)
            self.assertIn("| Spectra searched as | DDA (--dda) |", text)

    def test_a_dia_cfg_pinned_from_the_table_is_not_tagged(self):
        with tempfile.TemporaryDirectory() as d:
            _, cfg, prov = self._search(d, "DIA", "--ms1-resolution", "120000",
                                        "--ms2-resolution", "30000")
            rec = mm.search_record(cfg, prov)
            self.assertEqual((rec["ms1_tol"]["value"], rec["ms2_tol"]["value"]), (7, 15))
            self.assertNotIn("default", rec["ms1_tol"])
            self.assertNotIn("default", rec["ms2_tol"])
            self.assertNotIn(SOP_NOTE, mm.search_paragraph(rec))
            self.assertNotIn("DDA mode", mm.search_paragraph(rec))

    def test_the_run_records_parameter_table_carries_the_tag_too(self):
        import record_run
        with tempfile.TemporaryDirectory() as d:
            _, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            table = "\n".join(record_run.render_parameters(
                {"parameters": {"record": mm.search_record(cfg), "why": {}}}))
        self.assertIn(f"| MS2 mass accuracy | 20 ppm ({SOP_NOTE}) |", table)
        self.assertNotIn(f"10 ppm ({SOP_NOTE})", table)

    def test_with_no_run_record_the_cfgs_sidecar_is_read_and_a_stale_one_is_not(self):
        with tempfile.TemporaryDirectory() as d:
            _, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            rec = mm.search_record(cfg)
            self.assertEqual(rec["ms2_tol"].get("default"), SOP_NOTE)
            self.assertNotIn("default", rec["ms1_tol"])
            # the cfg edited to another MS2 value: the sidecar no longer describes it
            with open(cfg) as fh:
                text = fh.read().replace("--mass-acc 20\n", "--mass-acc 25\n")
            with open(cfg, "w") as fh:
                fh.write(text)
            rec = mm.search_record(cfg)
            self.assertEqual(rec["ms2_tol"]["value"], 25)
            self.assertNotIn("default", rec["ms2_tol"])


# ------------------------------------------------------------------------------------------
# (c) probe_window.py stops at once under --dda
# ------------------------------------------------------------------------------------------
class ProbeRefusesDda(unittest.TestCase):
    def test_no_dia_nn_is_started_and_it_returns_at_once(self):
        with tempfile.TemporaryDirectory() as d:
            marker = os.path.join(d, "diann_started")
            diann = fx._exe(os.path.join(d, "diann"),
                            f"#!/bin/bash\ntouch {marker}\nsleep 600\n")
            lib = os.path.join(d, "lib.speclib")
            with open(lib, "w") as fh:
                fh.write("lib")
            t0 = time.monotonic()
            p = subprocess.run([sys.executable, PROBE, "--diann", diann, "--raw", *cohort(d, 3),
                                "--fasta", fasta(d), "--lib", lib, "--timeout", "3600",
                                "--workdir", os.path.join(d, "w"),
                                "--", "--dda", "--qvalue", "0.01"],
                               capture_output=True, text=True, timeout=60)
            self.assertLess(time.monotonic() - t0, 30)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("--dda", p.stderr)
            self.assertIn("this probe cannot measure a DDA search", p.stderr)
            self.assertNotIn("nothing can be measured", p.stderr)
            self.assertIn("no scan-window radius is logged in DDA mode", p.stderr)
            self.assertIn("did not finish within 3600 s on the SET28 runs", p.stderr)
            self.assertFalse(os.path.exists(marker), "a DIA-NN probe was started")


# ------------------------------------------------------------------------------------------
# (b) the survey scan range, read by detect_acquisition.py from each format
# ------------------------------------------------------------------------------------------
class ThermoFilterRange(unittest.TestCase):
    def test_full_ms_filter_strings(self):
        f = da.ms1_filter_range
        self.assertEqual(f("FTMS + p NSI Full ms [350.0000-1500.0000]"), (350.0, 1500.0))
        self.assertEqual(f("FTMS + p NSI Full ms [350.0000-800.0000, 800.0000-1500.0000]"),
                         (350.0, 1500.0))
        self.assertIsNone(f("FTMS + c NSI d Full ms2 572.3181@hcd30.00 [118.9161-1189.1609]"))
        self.assertIsNone(f("FTMS + p NSI SIM ms [400.0000-420.0000]"))
        self.assertIsNone(f(None))



class ThermoDdaSurveyRange(trfp._FakeParserCase):
    def test_the_exploris_dda_survey_filter_is_the_range(self):
        r = da.classify(self.raw(trfp.DDA))
        self.assertEqual((r["acquisition"], r["confidence"]), ("DDA", "high"), r["reason"])
        self.assertEqual(r["precursor_mz_range"], [350.0, 1500.0])
        self.assertEqual(r["precursor_mz_range_source"], "ms1_survey_scan")

    def test_the_cli_hands_it_to_estimate_params(self):
        res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                              self.raw(trfp.DDA)], capture_output=True, text=True,
                             env=os.environ.copy(), timeout=120)
        self.assertEqual(res.returncode, 0, res.stderr)
        payload = json.loads(res.stdout)
        self.assertEqual(payload["overall"], "DDA")
        self.assertEqual(payload["precursor_mz_range"], [350.0, 1500.0])
        self.assertEqual(payload["precursor_mz_range_source"], "ms1_survey_scan")
        self.assertEqual(payload["precursor_mz_range_files_without"], [])

    def test_a_dda_slice_with_no_readable_survey_filter_says_so_in_a_note(self):
        """No filter string at all (FAKE_TRFP_NO_FILTER): still DDA, no range -- and a NOTE in
        the reason, not a DDA_NO_RANGE warning: Sage, the DDA default, does not use the range,
        and estimate_params.py refuses a DIA-NN DDA cfg without one."""
        os.environ["FAKE_TRFP_NO_FILTER"] = "1"
        r = da.classify(self.raw(trfp.DDA))
        self.assertEqual(r["acquisition"], "DDA", r["reason"])
        self.assertIsNone(r["precursor_mz_range"])
        self.assertIsNone(r["precursor_mz_range_source"])
        self.assertIn("NOTE: " + da.DDA_NO_RANGE, r["reason"])
        self.assertNotIn(da.DDA_NO_RANGE, r["warnings"])


MZML = """<?xml version="1.0"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml">
 <run><spectrumList count="{n}">{spectra}</spectrumList></run>
</mzML>
"""
MS1 = """
  <spectrum index="{i}">
   <cvParam accession="MS:1000511" name="ms level" value="1"/>
   <scanList><scan><scanWindowList count="1"><scanWindow>{window}</scanWindow></scanWindowList>
   </scan></scanList>
  </spectrum>"""
WINDOW = ('<cvParam accession="MS:1000501" name="scan window lower limit" value="{lo}"/>'
          '<cvParam accession="MS:1000500" name="scan window upper limit" value="{hi}"/>')
MS2 = """
  <spectrum index="{i}">
   <cvParam accession="MS:1000511" name="ms level" value="2"/>
   <scanList><scan><scanWindowList count="1"><scanWindow>
    <cvParam accession="MS:1000501" name="scan window lower limit" value="120"/>
    <cvParam accession="MS:1000500" name="scan window upper limit" value="1900"/>
   </scanWindow></scanWindowList></scan></scanList>
   <precursorList><precursor><isolationWindow>
    <cvParam accession="MS:1000827" name="isolation window target m/z" value="{tgt}"/>
    <cvParam accession="MS:1000828" name="isolation window lower offset" value="0.8"/>
    <cvParam accession="MS:1000829" name="isolation window upper offset" value="0.8"/>
   </isolationWindow></precursor></precursorList>
  </spectrum>"""


def dda_mzml(d, window=True):
    """Top-10 DDA: an MS1 over 350-1500 (msconvert's scan window terms, as in a real Thermo DDA
    mzML), then ten 1.6 m/z MS2 at distinct precursors, three cycles."""
    spectra, i = [], 0
    for cycle in range(3):
        spectra.append(MS1.format(i=i, window=WINDOW.format(lo=350, hi=1500) if window else ""))
        i += 1
        for k in range(10):
            spectra.append(MS2.format(i=i, tgt=400.0 + 37.1 * k + 3.3 * cycle))
            i += 1
    p = os.path.join(d, "dda.mzML")
    with open(p, "w") as fh:
        fh.write(MZML.format(n=len(spectra), spectra="".join(spectra)))
    return p


class MzmlSurveyRange(unittest.TestCase):
    def test_the_ms1_scan_window_is_the_range_and_ms2_windows_are_not(self):
        with tempfile.TemporaryDirectory() as d:
            kind, conf, why, rng = da.detect_mzml(dda_mzml(d))
            self.assertEqual(kind, "DDA", why)
            # 350-1500, not the MS2 fragment scan window 120-1900 or the precursor picks
            self.assertEqual(rng, (350.0, 1500.0))
            self.assertIn("survey scan range", why)

    def test_no_ms1_scan_window_is_no_range_and_a_note(self):
        with tempfile.TemporaryDirectory() as d:
            r = da.classify(dda_mzml(d, window=False))
            self.assertEqual(r["acquisition"], "DDA", r["reason"])
            self.assertIsNone(r["precursor_mz_range"])
            self.assertIn("NOTE: " + da.DDA_NO_RANGE, r["reason"])
            self.assertEqual(r["warnings"], [])

    def test_a_missing_dda_range_is_a_note_that_does_not_ask_for_confirmation(self):
        """F5: a Sage search -- the DDA default -- never uses the range, so its absence must not
        stop the run with needs_confirmation. DIA-NN is gated by estimate_params.py instead."""
        with tempfile.TemporaryDirectory() as d:
            p = dda_mzml(d, window=False)
            res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                                  p], capture_output=True, text=True, timeout=60)
            self.assertEqual(res.returncode, 0, res.stderr)
            payload = json.loads(res.stdout)
            self.assertFalse(payload["needs_confirmation"], payload)
            self.assertEqual(payload["dda_range_note"]["files"], [p])
            self.assertEqual(payload["dda_range_note"]["note"], da.DDA_NO_RANGE)
            with_range = subprocess.run([sys.executable, os.path.join(SCRIPTS,
                                         "detect_acquisition.py"), dda_mzml(d)],
                                        capture_output=True, text=True, timeout=60)
            self.assertIsNone(json.loads(with_range.stdout)["dda_range_note"])


BLOCK, N_FRAMES = 256, 400
MZ_ACQ = (99.993933, 1700.0)      # both SET1-28 runs below record exactly this


def synthetic_d(tmp, name, kind, acq_range=MZ_ACQ, dia_info_rows=None, pasef_rows=None,
                both_frame_types=False):
    """A finished timsTOF run shaped like the real ones read on HIVE, 2026-09-29 (PI_Example/
    SET1-28, timsTOF HT, immutable open):

        kind="dda"  09292026__30SPD_DDANS27: MsMsType {0, 8}; PasefFrameMsMsInfo rows; the three
                    dia-PASEF tables PRESENT with 0 rows
        kind="dia"  09112026__60SPD_DIA-NS10: MsMsType {0, 9}; DiaFrameMsMsInfo one row per
                    dia-PASEF frame, DiaFrameMsMsWindows and -WindowGroups filled; no Pasef tables

    `dia_info_rows` / `pasef_rows` override whether those frame tables get rows (for runs whose
    evidence disagrees); `both_frame_types` makes every tenth MS2 frame the other type. The
    tdf_bin matches the index, so tdf_integrity() calls it ok."""
    dia = kind == "dia"
    dia_info_rows = dia if dia_info_rows is None else dia_info_rows
    pasef_rows = (not dia) if pasef_rows is None else pasef_rows
    ms2, other = (9, 8) if dia else (8, 9)
    d = os.path.join(tmp, name)
    os.makedirs(d)
    with open(os.path.join(d, "analysis.tdf_bin"), "wb") as fh:
        for _ in range(N_FRAMES):
            fh.write(struct.pack("<II", BLOCK, 5) + b"\0" * (BLOCK - 8))
    tdf = os.path.join(d, "analysis.tdf")
    con = sqlite3.connect(synthetic_tdf_write_uri(tdf), uri=True)
    meta = [("InstrumentName", "timsTOF HT")]
    if acq_range:
        meta += [("MzAcqRangeLower", str(acq_range[0])), ("MzAcqRangeUpper", str(acq_range[1]))]
    con.execute("CREATE TABLE GlobalMetadata (Key TEXT PRIMARY KEY, Value TEXT)")
    con.executemany("INSERT INTO GlobalMetadata VALUES (?,?)", meta)
    con.execute("CREATE TABLE Frames (Id INTEGER PRIMARY KEY, Time REAL, MsMsType INTEGER, "
                "TimsId INTEGER, NumScans INTEGER, AccumulationTime REAL, RampTime REAL)")
    frames = [(i + 1, i * 0.1, 0 if i % 5 == 0 else
               other if both_frame_types and i % 10 == 1 else ms2,
               i * BLOCK, 5, 100.0, 100.0) for i in range(N_FRAMES)]
    con.executemany("INSERT INTO Frames VALUES (?,?,?,?,?,?,?)", frames)
    ms2_ids = [f[0] for f in frames if f[2] in (8, 9)]
    # The dia-PASEF tables: in every dia-PASEF run, and EMPTY in the real ddaPASEF run
    con.execute("CREATE TABLE DiaFrameMsMsInfo (Frame INTEGER, WindowGroup INTEGER)")
    con.execute("CREATE TABLE DiaFrameMsMsWindowGroups (Id INTEGER)")
    con.execute("CREATE TABLE DiaFrameMsMsWindows (WindowGroup INTEGER, ScanNumBegin INTEGER, "
                "ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, "
                "CollisionEnergy REAL)")
    if dia:
        con.executemany("INSERT INTO DiaFrameMsMsWindowGroups VALUES (?)", [(g,) for g in (1, 2)])
        con.executemany("INSERT INTO DiaFrameMsMsWindows VALUES (?,?,?,?,?,?)",
                        [(1 + k % 2, 0, 100, 412.5 + 25.0 * k, 25.0, 30.0) for k in range(32)])
    if dia_info_rows:
        con.executemany("INSERT INTO DiaFrameMsMsInfo VALUES (?,?)",
                        [(f, 1 + f % 2) for f in ms2_ids])
    if not dia or pasef_rows:
        con.execute("CREATE TABLE PasefFrameMsMsInfo (Frame INTEGER, ScanNumBegin INTEGER, "
                    "ScanNumEnd INTEGER, IsolationMz REAL, IsolationWidth REAL, "
                    "CollisionEnergy REAL, Precursor INTEGER)")
    if pasef_rows:
        con.executemany("INSERT INTO PasefFrameMsMsInfo VALUES (?,?,?,?,?,?,?)",
                        [(f, 0, 100, 500.0 + k, 2.0, 30.0, k) for k, f in enumerate(ms2_ids)])
    con.commit()
    con.close()
    for side_file in ("-wal", "-shm"):
        if os.path.exists(tdf + side_file):
            os.remove(tdf + side_file)
    return d


class BrukerAcquisitionAndRange(unittest.TestCase):
    def classify(self, *args, **kw):
        with tempfile.TemporaryDirectory() as d:
            r = da.classify(synthetic_d(d, *args, **kw))
        self.assertEqual(r["tdf_integrity"]["status"], "ok", r["tdf_integrity"])
        return r

    def test_the_real_ddapasef_shape_is_dda_with_the_ms1_acquisition_range(self):
        """msalemi 2026-09-29: the three dia-PASEF tables present but EMPTY, beside PASEF rows
        and MsMsType-8 frames, came back "DIA, high confidence". Table existence is not
        acquisition."""
        r = self.classify("ddapasef.d", "dda")
        self.assertEqual((r["acquisition"], r["confidence"]), ("DDA", "high"), r["reason"])
        self.assertEqual(r["precursor_mz_range"], list(MZ_ACQ))
        self.assertEqual(r["precursor_mz_range_source"], "ms1_survey_scan")
        self.assertIn("MzAcqRange", r["reason"])
        self.assertEqual(r["warnings"], [])

    def test_the_real_dia_pasef_shape_is_still_dia_with_its_window_range(self):
        r = self.classify("diapasef.d", "dia")
        self.assertEqual((r["acquisition"], r["confidence"]), ("DIA", "high"), r["reason"])
        # the isolation windows (32 x 25 m/z from 400 to 1200), not the MS1 acquisition range
        self.assertEqual(r["precursor_mz_range"], [400.0, 1200.0])
        self.assertEqual(r["precursor_mz_range_source"], "isolation_windows")

    def test_frames_of_both_types_are_not_guessed(self):
        for kind in ("dda", "dia"):
            r = self.classify(f"{kind}.d", kind, both_frame_types=True)
            self.assertEqual((r["acquisition"], r["confidence"]), ("unknown", "low"), r["reason"])
            self.assertIn("disagrees", r["reason"])
            self.assertIsNone(r["precursor_mz_range"])

    def test_typed_frames_beside_the_other_kinds_frame_table_are_not_guessed(self):
        dda_frames_dia_rows = self.classify("a.d", "dda", dia_info_rows=True)
        dia_frames_pasef_rows = self.classify("b.d", "dia", pasef_rows=True)
        for r in (dda_frames_dia_rows, dia_frames_pasef_rows):
            self.assertEqual((r["acquisition"], r["confidence"]), ("unknown", "low"), r["reason"])
            self.assertIn("confirm the acquisition method", r["reason"])

    def test_disagreement_asks_the_user(self):
        with tempfile.TemporaryDirectory() as d:
            p = synthetic_d(d, "x.d", "dda", dia_info_rows=True)
            res = subprocess.run([sys.executable, os.path.join(SCRIPTS, "detect_acquisition.py"),
                                  p], capture_output=True, text=True, timeout=60)
        self.assertEqual(res.returncode, 0, res.stderr)
        payload = json.loads(res.stdout)
        self.assertTrue(payload["needs_confirmation"])
        self.assertEqual(payload["overall"], "unknown")

    def test_no_acquisition_range_is_no_range_and_a_note(self):
        r = self.classify("dda.d", "dda", acq_range=None)
        self.assertEqual(r["acquisition"], "DDA", r["reason"])
        self.assertIsNone(r["precursor_mz_range"])
        self.assertIn("NOTE: " + da.DDA_NO_RANGE, r["reason"])
        self.assertEqual(r["warnings"], [])


# ------------------------------------------------------------------------------------------
# (11) submit.sh with a reader that stops early
# ------------------------------------------------------------------------------------------
class SubmitShSurvivesAClosedStdout(unittest.TestCase):
    def test_jobs_txt_and_the_checkpoint_are_written_when_stdout_is_closed(self):
        """`bash submit.sh | head -1`: the reader is gone after one line. Here the pipe's read
        end is closed before the script starts, so its FIRST echo already hits a broken pipe --
        the worst case. Every job id must still reach jobs.txt, and RECOVERY.md must exist."""
        with tempfile.TemporaryDirectory() as d:
            _, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            sess = os.path.join(d, "sess")
            out = os.path.join(sess, "output", "search")
            p = generate(d, cfg, out)
            self.assertEqual(p.returncode, 0, p.stderr)
            fake = os.path.join(d, "fakebin")
            os.makedirs(fake)
            fx._exe(os.path.join(fake, "sbatch"), FAKE_SBATCH)
            env = job_env(d, PATH=fake + os.pathsep + os.environ.get("PATH", ""),
                          SB_COUNTER=os.path.join(d, "counter"))
            r_end, w_end = os.pipe()
            os.close(r_end)
            try:
                res = subprocess.run(["bash", os.path.join(out, "submit.sh")], stdout=w_end,
                                     stderr=subprocess.PIPE, text=True, env=env, timeout=60)
            finally:
                os.close(w_end)
            self.assertEqual(res.returncode, 0, res.stderr)
            with open(os.path.join(out, "jobs.txt")) as fh:
                self.assertEqual(fh.read().split(), ["701", "702", "703", "704", "705"])
            self.assertTrue(os.path.exists(os.path.join(sess, "RECOVERY.md")))

    def test_the_record_comes_before_any_message(self):
        with tempfile.TemporaryDirectory() as d:
            _, cfg = estimate(d, "DDA", *SET28, *SURVEY)
            out = os.path.join(d, "sess", "output", "search")
            self.assertEqual(generate(d, cfg, out).returncode, 0)
            lines = read(os.path.join(out, "submit.sh")).splitlines()
            last_sbatch = max(i for i, ln in enumerate(lines) if "sbatch --parsable" in ln)
            first_message = min(i for i, ln in enumerate(lines)
                                if ln.startswith(("echo ", "say ")) and i > last_sbatch)
            jobs_txt = next(i for i, ln in enumerate(lines) if "jobs.txt\"" in ln
                            and ln.startswith("printf"))
            record = next(i for i, ln in enumerate(lines) if "checkpoint.py" in ln)
            self.assertLess(jobs_txt, first_message)
            self.assertLess(record, first_message)
            self.assertFalse([ln for ln in lines if ln.startswith("echo ")],
                             "a bare echo can still die of SIGPIPE under set -e")


if __name__ == "__main__":
    unittest.main()
