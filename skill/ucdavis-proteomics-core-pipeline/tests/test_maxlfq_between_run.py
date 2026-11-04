#!/usr/bin/env python3
"""The MaxLFQ route's between-run step (2.14.2, I20; Brett's Q1).

Up to 2.14.1, `run_de.R --method maxlfq --quantities normalised` quantile-normalised every
report's PG.MaxLFQ -- on top of DIA-NN's own cross-run normalisation for a normal DIA-NN report.
On the PXD063416 yeast-spike benchmark (4 m/z) that extra step made 61% of the significant calls
wrong; the default limpa route got 2%. Now ONE rule, build_maxlfq.R maxlfq_between_run(), decides
the step (run_de.R --between-run auto|none|quantile|median, default auto):

  * DIA-NN report, PG.MaxLFQ normalised by DIA-NN          -> none (no second step)
  * DIA-NN report searched with --no-norm                   -> quantile (the only normalisation)
  * DIA-NN report, whether DIA-NN normalised NOT RECORDED   -> none, said "[confirm]"
  * a report that declares its own quantity (Sage, FragPipe) -> quantile, as in 2.14.1
  * --quantities raw                                        -> none; any other step is refused
  * --between-run quantile | median | none                  -> as given, said "chosen for this DE"

and the descriptor's `normalisation` text, its new `between_run` field and the reproducibility
log are all built from that answer. --between-run is refused with --method dpc.

Synthetic DIA-NN-shaped reports with no missing values (so a quantile step gives identical column
distributions). Needs pyarrow and R with limma/arrow/dplyr/tidyr/jsonlite; skips otherwise.
Hermetic.
"""
import csv
import json
import math
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

from test_run_de_contaminants import r_has  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

NEEDS = ("limma", "arrow", "dplyr", "tidyr", "jsonlite")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")
RUNS = [("A_1", "A"), ("A_2", "A"), ("A_3", "A"), ("B_1", "B"), ("B_2", "B"), ("B_3", "B")]
Q_COLS = ("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.Q.Value", "Global.Q.Value",
          "Global.PG.Q.Value")
SAGE_META = {"delimp.quantity.level": "peptide", "delimp.quantity.label": "Sage",
             "delimp.quantity.value": "Sage label-free peptide intensity",
             "delimp.quantity.identification": "Sage's own q-values, 1% FDR",
             "delimp.quantity.q_columns": "placeholder",
             "delimp.quantity.citation": "Sage (Lazear 2023)"}


def write_report(out, normalised=True, quant_columns=True, declared=None, n_prot=40, seed=3):
    """A complete (no missing values) DIA-NN-shaped report.parquet + conditions.csv in `out`.
    Each run carries a loading offset; `normalised`: DIA-NN's normalisation removed it from
    Precursor.Normalised and PG.MaxLFQ (PG.MaxLFQ = sum of the normalised precursors); otherwise
    the report is a --no-norm one (Normalised == Quantity). `quant_columns` False drops
    Precursor.Quantity / Precursor.Normalised, so nothing in the data says which it was.
    `declared`: parquet schema metadata (an adapted report declaring its quantity).
    -> (report, conditions, {(protein, run): PG.MaxLFQ})."""
    rnd = random.Random(seed)
    os.makedirs(out, exist_ok=True)
    load = {r: rnd.uniform(-0.8, 0.8) for r, _ in RUNS}
    rows = {k: [] for k in ("Run", "Precursor.Id", "Protein.Group", "Protein.Ids",
                            "Protein.Names", "Genes", "Proteotypic", "Precursor.Quantity",
                            "Precursor.Normalised", "PG.MaxLFQ") + Q_COLS}
    pg = {}
    for i in range(n_prot):
        acc, base = f"P{i:05d}", rnd.uniform(16, 22)
        for r, grp in RUNS:
            lv = base + (1.0 if grp == "B" and i < 8 else 0.0) + rnd.gauss(0, 0.2)
            precs = []
            for k in range(2):
                q = 2 ** (lv + 0.5 * k + load[r] + rnd.gauss(0, 0.1))
                n = q / (2 ** load[r]) if normalised else q
                precs.append((f"PEP{i}K{k}2", q, n))
            pg[(acc, r)] = sum(n for _, _, n in precs)
            for pid, q, n in precs:
                for c, x in (("Run", r), ("Precursor.Id", pid), ("Protein.Group", acc),
                             ("Protein.Ids", acc), ("Protein.Names", f"G{i}_HUMAN"),
                             ("Genes", f"G{i}"), ("Proteotypic", 1),
                             ("Precursor.Quantity", q), ("Precursor.Normalised", n),
                             ("PG.MaxLFQ", pg[(acc, r)])):
                    rows[c].append(x)
                for c in Q_COLS:
                    rows[c].append(0.001)
    if not quant_columns:
        del rows["Precursor.Quantity"], rows["Precursor.Normalised"]
    t = pa.table(rows)
    if declared:
        t = t.replace_schema_metadata(declared)
    report = os.path.join(out, "report.parquet")
    pq.write_table(t, report)
    cond = os.path.join(out, "conditions.csv")
    with open(cond, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        w.writerows(RUNS)
    return report, cond, pg


R_CALL = r'''
setwd("{scripts}")
source("build_maxlfq.R")
m <- build_maxlfq("{report}", quantities = "{quantities}", between_run = "{between_run}")
write.csv(m$E, file.path("{out}", "E.csv"))
jsonlite::write_json(m$descriptor, file.path("{out}", "descriptor.json"), auto_unbox = TRUE,
                     null = "null", na = "null", digits = NA)
'''


def run_maxlfq(report, out, between_run="auto", quantities="normalised"):
    """build_maxlfq() on `report` -> ({protein: {run: log2 value}}, descriptor)."""
    os.makedirs(out, exist_ok=True)
    p = subprocess.run(["Rscript", "-e", R_CALL.format(scripts=SCRIPTS, report=report, out=out,
                                                       quantities=quantities,
                                                       between_run=between_run)],
                       capture_output=True, text=True, timeout=300)
    if p.returncode != 0:
        raise AssertionError(p.stderr[-2000:])
    with open(os.path.join(out, "E.csv"), newline="") as fh:
        r = csv.reader(fh)
        runs = next(r)[1:]
        E = {row[0]: {run: float(v) for run, v in zip(runs, row[1:])} for row in r}
    with open(os.path.join(out, "descriptor.json")) as fh:
        return E, json.load(fh)


def columns(E):
    runs = sorted(next(iter(E.values())))
    return {r: [E[p][r] for p in E] for r in runs}


def median(xs):
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class TheRule(unittest.TestCase):
    """build_maxlfq(between_run = ...) and the descriptor it writes (maxlfq_between_run)."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_dia_nn_normalised_report_gets_no_second_step(self):
        report, _, pg = write_report(os.path.join(self.d, "n"))
        E, desc = run_maxlfq(report, os.path.join(self.d, "o"))      # default --between-run auto
        self.assertEqual(desc["between_run"]["method"], "none")
        self.assertEqual(desc["between_run"]["requested"], "auto")
        self.assertIn("already carries DIA-NN's cross-run normalisation", desc["between_run"]["why"])
        for (acc, run), v in pg.items():                             # untouched log2 PG.MaxLFQ
            self.assertAlmostEqual(E[acc][run], math.log2(v), places=9)
        self.assertIn("no second step", desc["normalisation"])
        self.assertNotIn("quantile", desc["normalisation"])

    def test_explicit_quantile_is_still_available(self):              # the opt-in kept (Q1)
        report, _, _ = write_report(os.path.join(self.d, "n"))
        E, desc = run_maxlfq(report, os.path.join(self.d, "o"), between_run="quantile")
        cols = list(columns(E).values())
        for c in cols[1:]:                                            # identical distributions
            for a, b in zip(sorted(c), sorted(cols[0])):
                self.assertAlmostEqual(a, b, places=9)
        self.assertEqual(desc["between_run"]["method"], "quantile")
        self.assertEqual(desc["between_run"]["why"], "--between-run quantile (chosen for this DE)")
        self.assertIn("DIA-NN's cross-run normalisation (PG.MaxLFQ), then quantile normalisation "
                      "(limma::normalizeBetweenArrays)", desc["normalisation"])
        self.assertIn("chosen for this DE", desc["normalisation"])

    def test_a_no_norm_report_keeps_its_quantile_step(self):
        report, _, _ = write_report(os.path.join(self.d, "nn"), normalised=False)
        _, desc = run_maxlfq(report, os.path.join(self.d, "o"))
        self.assertEqual(desc["between_run"]["method"], "quantile")
        self.assertEqual(desc["between_run"]["why"],
                         "the only between-run normalisation: the search ran with --no-norm")
        self.assertTrue(desc["normalisation"].startswith(
            "quantile normalisation (limma::normalizeBetweenArrays) only: the DIA-NN search ran "
            "with --no-norm"))

    def test_not_recorded_gets_no_step_and_says_confirm(self):
        report, _, _ = write_report(os.path.join(self.d, "u"), quant_columns=False)
        _, desc = run_maxlfq(report, os.path.join(self.d, "o"))
        self.assertEqual(desc["between_run"]["method"], "none")
        self.assertIn("NOT RECORDED [confirm]", desc["between_run"]["why"])
        self.assertIn("no second step", desc["normalisation"])
        self.assertIn("NOT RECORDED [confirm]", desc["normalisation"])

    def test_declared_sage_quantities_keep_quantile(self):            # unchanged until measured
        report, _, _ = write_report(os.path.join(self.d, "s"), declared=SAGE_META)
        _, desc = run_maxlfq(report, os.path.join(self.d, "o"))
        self.assertEqual(desc["pipeline_id"], "peptide_max")
        self.assertEqual(desc["between_run"]["method"], "quantile")
        self.assertEqual(desc["normalisation"],
                         "quantile normalisation (limma::normalizeBetweenArrays) of the engine's "
                         "quantities")

    def test_median_centres_columns(self):
        report, _, pg = write_report(os.path.join(self.d, "nn"), normalised=False)
        E, desc = run_maxlfq(report, os.path.join(self.d, "o"), between_run="median")
        before = [median([math.log2(v) for (a, r), v in pg.items() if r == run]) for run, _ in RUNS]
        after = [median(c) for c in columns(E).values()]
        for m in after:                                              # every run's median ...
            self.assertAlmostEqual(m, median(before), places=9)      # ... the median of the medians
        self.assertEqual(desc["between_run"]["method"], "median")
        self.assertIn("median", desc["normalisation"])
        self.assertIn("--between-run median, chosen for this DE", desc["normalisation"])

    def test_explicit_none_on_a_no_norm_report_is_no_normalisation_at_all(self):
        report, _, pg = write_report(os.path.join(self.d, "nn"), normalised=False)
        E, desc = run_maxlfq(report, os.path.join(self.d, "o"), between_run="none")
        for (acc, run), v in pg.items():
            self.assertAlmostEqual(E[acc][run], math.log2(v), places=9)
        self.assertEqual(desc["between_run"]["method"], "none")
        self.assertTrue(desc["normalisation"].startswith("none:"), desc["normalisation"])
        self.assertIn("--no-norm", desc["normalisation"])

    def test_raw_with_a_step_is_refused_by_the_rule(self):
        report, _, _ = write_report(os.path.join(self.d, "nn"), normalised=False)
        with self.assertRaises(AssertionError) as e:
            run_maxlfq(report, os.path.join(self.d, "o"), between_run="quantile", quantities="raw")
        self.assertIn("--between-run quantile", str(e.exception))
        self.assertIn("--quantities raw", str(e.exception))
        _, desc = run_maxlfq(report, os.path.join(self.d, "o2"), between_run="none",
                             quantities="raw")
        self.assertEqual(desc["between_run"]["method"], "none")


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class RunDe(unittest.TestCase):
    """run_de.R --between-run: refused where it has no meaning, recorded and reproduced where it
    has."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def de(self, report, cond, *extra):
        out = os.path.join(self.d, "de")
        p = subprocess.run(["Rscript", RUN_DE, "--input", report, "--metadata", cond,
                            "--contrasts", "B-A", "--outdir", out, *extra],
                           capture_output=True, text=True, timeout=900, cwd=self.d)
        return p, out

    def test_raw_refuses_a_between_run_step(self):
        report, cond, _ = write_report(os.path.join(self.d, "nn"), normalised=False)
        p, _ = self.de(report, cond, "--method", "maxlfq", "--quantities", "raw",
                       "--between-run", "quantile")
        self.assertNotEqual(p.returncode, 0)
        flat = " ".join(p.stderr.split())
        self.assertIn("--between-run quantile", flat)
        self.assertIn("--quantities raw", flat)

    def test_dpc_refuses_between_run(self):
        report, cond, _ = write_report(os.path.join(self.d, "n"))
        p, _ = self.de(report, cond, "--method", "dpc", "--between-run", "none")
        self.assertNotEqual(p.returncode, 0)
        flat = " ".join(p.stderr.split())
        self.assertIn("--between-run", flat)
        self.assertIn("--method maxlfq", flat)

    def test_an_unknown_value_is_refused(self):
        report, cond, _ = write_report(os.path.join(self.d, "n"))
        p, _ = self.de(report, cond, "--method", "maxlfq", "--between-run", "loess")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("auto, none, quantile or median", " ".join(p.stderr.split()))

    def test_the_default_records_no_second_step_everywhere(self):
        report, cond, _ = write_report(os.path.join(self.d, "n"))
        p, out = self.de(report, cond, "--method", "maxlfq")
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["between_run"]["method"], "none")
        # the requested step is part of the design the normalisation gate holds a final DE to
        self.assertEqual(prov["design_request"]["between_run"], "auto")
        self.assertIn("no second step", prov["normalisation"])
        with open(os.path.join(out, "methods.txt")) as fh:
            self.assertIn("no second step", " ".join(fh.read().split()))
        with open(os.path.join(out, "reproducibility_log.R")) as fh:
            src = fh.read()
        self.assertNotIn("normalizeBetweenArrays", src)
        self.assertIn("no between-run step", src)

    def test_median_is_recorded_and_reproduced(self):
        report, cond, _ = write_report(os.path.join(self.d, "nn"), normalised=False)
        p, out = self.de(report, cond, "--method", "maxlfq", "--between-run", "median")
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual((prov["between_run"]["method"], prov["between_run"]["requested"]),
                         ("median", "median"))
        self.assertEqual(prov["design_request"]["between_run"], "median")
        with open(os.path.join(out, "reproducibility_log.R")) as fh:
            src = fh.read()
        self.assertNotIn("normalizeBetweenArrays", src)
        self.assertIn("stats::median", src)


@unittest.skipUnless(r_has("jsonlite"), "needs Rscript")
class ReproScriptNamesTheStep(unittest.TestCase):
    """repro_script.R write_repro_script() emits the step the run's descriptor recorded. Given
    neither the descriptor's between_run nor between_run =, a maxlfq script is refused: it once
    fell back to quantile silently -- a default reaching the reproducibility log (rule 2), and a
    second copy of provenance.py's legacy rule (normalization_check.between_run_of)."""
    R = r'''
      source("{scripts}/diann_q_columns.R"); source("{scripts}/repro_script.R")
      meta <- data.frame(File.Name = c("r1","r2","r3","r4"), Group = c("A","A","B","B"))
      invisible(write_repro_script(path = "{out}", method = "maxlfq", input = "report.parquet",
        format = "parquet", q_cutoff = 0.01, q_columns = "Q.Value", q_cutoffs = 0.01,
        eq_cutoff = 0, pgq_cutoff = 0, cov_min_frac = 0.5, meta = meta,
        covariates = character(0), formula_parts = "groups", forms = "B-A",
        adjp_thr = 0.05, logfc_ref = 1{extra}))
      cat(readLines("{out}"), sep = "\n")'''

    def emit(self, extra=""):
        with tempfile.TemporaryDirectory() as d:
            return subprocess.run(["Rscript", "-e", self.R.format(
                scripts=SCRIPTS, out=os.path.join(d, "r.R"), extra=extra)],
                capture_output=True, text=True, timeout=120)

    def test_no_step_named_is_refused(self):
        p = self.emit()
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("between-run step", " ".join(p.stderr.split()))

    def test_the_descriptors_step_is_emitted(self):
        p = self.emit(', descriptor = list(between_run = list(method = "median", '
                      'requested = "median", why = "x"))')
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertIn("E <- sweep(E, 2, cm - stats::median(cm, na.rm = TRUE))", p.stdout)
        p = self.emit(', between_run = "quantile"')
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertIn("normalizeBetweenArrays", p.stdout)


class ReproduceSh(unittest.TestCase):
    """reproduce.sh (provenance.py) replays the between-run step the DE APPLIED, always, for a
    normalised maxlfq DE: a later skill's auto may give another step (2.14.2 changed it once), so a
    replay relying on it could silently be another analysis. A DE that read a --no-norm report's
    PG.MaxLFQ is noted: the replayed search is DIA-NN-normalised, so its step would sit on top of
    DIA-NN's. Stdlib only."""

    def replay(self, between_run, quantities="normalised", **more):
        """-> (the replayed --between-run, the replay's stderr, reproduce.sh)."""
        from test_reproduce_sh_runs import arg, provenance, run_steps_5_and_6
        with tempfile.TemporaryDirectory() as d:
            de = os.path.join(d, "tables")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump(dict({"method": "maxlfq", "quantities": quantities,
                                "between_run": between_run}, **more), fh)
            out = os.path.join(d, "repro")
            sh = provenance(out, de, "--de-method", "maxlfq")
            calls, err = run_steps_5_and_6(d, out, sh)
        de_call = calls[1]["argv"]
        self.assertTrue(de_call[0].endswith("/scripts/run_de.R"))
        return arg(de_call, "--between-run"), err, sh

    def test_an_explicit_step_is_replayed(self):
        for m in ("median", "quantile", "none"):
            with self.subTest(m=m):
                self.assertEqual(self.replay({"method": m, "requested": m,
                                              "why": f"--between-run {m} (chosen for this DE)"})[0],
                                 m)

    def test_the_step_that_ran_is_spelled_out_for_auto_too(self):
        for applied in ("none", "quantile"):
            with self.subTest(applied=applied):
                self.assertEqual(self.replay({"method": applied, "requested": "auto",
                                              "why": "x"})[0], applied)

    def test_a_raw_de_replays_no_step(self):
        self.assertIsNone(self.replay({"method": "none", "requested": "auto", "why": "x"},
                                      quantities="raw")[0])

    def test_a_de_from_before_2_14_2_replays_its_quantile_step(self):
        """It records no between_run, and every normalised maxlfq DE was quantile-normalised then:
        the replay asks for that step, which the current run_de.R would not apply by itself."""
        self.assertEqual(self.replay(None)[0], "quantile")

    def test_a_no_norm_reports_normalised_de_is_noted(self):
        flag, err, sh = self.replay({"method": "quantile", "requested": "auto", "why": "x"},
                                    dia_nn_normalisation={"no_norm": True, "source": "the name"})
        self.assertEqual(flag, "quantile")
        self.assertIn("re-quantify with --no-norm", err)
        self.assertIn("--no-norm report's PG.MaxLFQ (--method maxlfq --quantities normalised", sh)
        flag, err, sh = self.replay({"method": "none", "requested": "none", "why": "x"},
                                    dia_nn_normalisation={"no_norm": True, "source": "the name"})
        self.assertEqual(flag, "none")
        self.assertIn("re-quantify with --no-norm", err)
        self.assertIn("no between-run normalisation at all", sh)
        _, err, _ = self.replay({"method": "none", "requested": "auto", "why": "x"},
                                dia_nn_normalisation={"no_norm": False, "source": "x"})
        self.assertNotIn("re-quantify", err)


if __name__ == "__main__":
    unittest.main()
