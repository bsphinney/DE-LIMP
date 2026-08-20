#!/usr/bin/env python3
"""
CONTRACT: every analysis path the skill can run must emit a reproducibility
script that reproduces it.

This is a guard for FUTURE work, not a fix for a past bug. Two reproducibility
defects shipped in a row, both with the same shape:

  v2.0.1  dpc     emitted script gave 4,624 proteins / 1,846 significant
                  against the run's 4,648 / 1,859
  v2.0.2  maxlfq  emitted script gave 4,444 / 1,281
                  against the run's 4,526 / 1,326

Neither errored. Both produced correct CSVs beside a script that disagreed with
them. The second happened *because* the first was fixed on one branch only --
`repro_script.R` has a separate emitter per --method, and the maxlfq branch was
missed. That is the failure this file exists to make impossible:

  1. Every --method run_de.R accepts must be listed in COVERED. Adding a third
     method fails this test until someone extends the coverage, rather than
     shipping an unverified emitter.
  2. Every cutoff a run APPLIES must appear in the script it emits, for each
     method, whatever shape that method's emitter uses (readDIANN q.cutoffs
     vector for dpc, a dplyr::filter chain for maxlfq).

Structural checks need Rscript but NOT limpa/arrow, so they run in skill-checks
CI. The end-to-end check -- run the pipeline for real, execute the emitted
script, diff the numbers -- needs the full stack and skips when it is absent.
Run it locally with:

    DELIMP_TEST_REPORT=/path/to/report.parquet \\
      python3 tests/test_reproducibility_contract.py TestEndToEndReproduction
"""
import os
import re
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPTS = os.path.join(ROOT, "scripts")

# Every --method that must have verified reproducibility. Extend this ONLY
# together with real coverage below; the first test fails otherwise.
COVERED = {"dpc", "maxlfq"}

ALL_SIX = ["Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
           "PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value"]
MIXED = [0.01, 0.01, 0.01, 0.05, 0.01, 0.01]


def have(*progs):
    return all(shutil.which(p) for p in progs)


def r_has(*pkgs):
    if not have("Rscript"):
        return False
    expr = ";".join(f'q(status=if(requireNamespace("{p}",quietly=TRUE)) 0 else 1)'
                    for p in pkgs[:1])
    for p in pkgs:
        r = subprocess.run(["Rscript", "-e",
                            f'quit(status = if (requireNamespace("{p}", quietly=TRUE)) 0 else 1)'],
                           capture_output=True)
        if r.returncode != 0:
            return False
    return True


def declared_methods():
    """The methods run_de.R accepts, read from its own validation line.

    Parsed rather than hardcoded so a new method shows up here automatically —
    that is the whole point: the test must notice work nobody told it about.
    """
    with open(os.path.join(SCRIPTS, "run_de.R")) as fh:
        src = fh.read()
    m = re.search(r'method\s*%in%\s*c\(([^)]*)\)', src)
    assert m, "could not find run_de.R's --method validation; update this parser"
    return set(re.findall(r'"([^"]+)"', m.group(1)))


def emit(method, q_columns=ALL_SIX, q_cutoffs=MIXED, q_cutoff=0.01,
         eq_cutoff=0, pgq_cutoff=0):
    """Generate a reproducibility script and return its source."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "reproducibility_log.R")
        cols = ", ".join(f'"{c}"' for c in q_columns)
        cuts = "NULL" if q_cutoffs is None else "c(" + ", ".join(map(str, q_cutoffs)) + ")"
        r = f'''
          source("{SCRIPTS}/diann_q_columns.R"); source("{SCRIPTS}/repro_script.R")
          meta <- data.frame(File.Name = c("r1","r2","r3","r4"),
                             Group = c("A","A","B","B"), stringsAsFactors = FALSE)
          invisible(write_repro_script(
            path = "{out}", method = "{method}", input = "report.parquet",
            format = "parquet", q_cutoff = {q_cutoff},
            q_columns = c({cols}), q_cutoffs = {cuts},
            eq_cutoff = {eq_cutoff}, pgq_cutoff = {pgq_cutoff}, cov_min_frac = 0.5,
            meta = meta, covariates = character(0),
            formula_parts = "groups", forms = "A-B",
            adjp_thr = 0.05, logfc_ref = 1))
          cat(readLines("{out}"), sep = "\\n")
        '''
        p = subprocess.run(["Rscript", "-e", r], capture_output=True, text=True)
        if p.returncode != 0:
            raise AssertionError(f"emitter failed for method={method}:\n{p.stderr[-1500:]}")
        return p.stdout


class TestEveryMethodIsCovered(unittest.TestCase):
    def test_no_method_ships_without_reproducibility_coverage(self):
        declared = declared_methods()
        uncovered = declared - COVERED
        self.assertEqual(
            uncovered, set(),
            f"run_de.R accepts --method {sorted(uncovered)} with no reproducibility "
            f"coverage. Add it to COVERED here and verify its emitted script "
            f"reproduces a real run before shipping — the maxlfq emitter was "
            f"missed exactly this way and shipped wrong for a release.")

    def test_covered_list_has_not_gone_stale(self):
        # The inverse: a method removed from run_de.R should not linger here
        # pretending to be tested.
        self.assertEqual(COVERED - declared_methods(), set())


class TestEmittedScriptDeclaresEveryAppliedCutoff(unittest.TestCase):
    """The invariant both shipped bugs violated, checked for every method."""

    def setUp(self):
        if not have("Rscript"):
            self.skipTest("Rscript not available")

    def test_every_method_emits_the_per_column_cutoffs_it_applied(self):
        for method in sorted(COVERED):
            with self.subTest(method=method):
                src = emit(method)
                for col, cut in zip(ALL_SIX, MIXED):
                    # Two emitter shapes to satisfy:
                    #   dpc     c('PG.Q.Value' = 0.05, ...)   quoted, "="
                    #   maxlfq  dplyr::filter(PG.Q.Value <= 0.05, ...)  bare, "<="
                    # Anchored on the preceding separator because
                    # "Lib.PG.Q.Value" ENDS WITH "PG.Q.Value" -- an unanchored
                    # check silently passes on the wrong column.
                    pat = (r"[,(]\s*'?" + re.escape(col) + r"'?\s*(<=|=)\s*" +
                           re.escape(str(cut)))
                    self.assertRegex(
                        src, pat,
                        f"method={method}: emitted script does not apply "
                        f"{col} at {cut}; it will not reproduce the run")

    def test_no_method_emits_a_bare_scalar_for_mixed_cutoffs(self):
        # The exact shape of both bugs: one number standing in for six.
        for method in sorted(COVERED):
            with self.subTest(method=method):
                src = emit(method)
                self.assertNotRegex(
                    src, r"q\.cutoffs\s*=\s*0\.01\s*,",
                    f"method={method}: emitted a scalar q.cutoffs alongside "
                    f"multiple columns; limpa recycles it and the script will "
                    f"silently filter on the wrong thresholds")

    def test_every_method_emits_quantums_filters_when_applied(self):
        # Generalises past q-values: any filter the run applies must appear.
        for method in sorted(COVERED):
            with self.subTest(method=method):
                src = emit(method, eq_cutoff=0.75, pgq_cutoff=0.5)
                self.assertIn("Empirical.Quality", src)
                self.assertIn("0.75", src)

    def test_every_method_emits_a_runnable_script(self):
        for method in sorted(COVERED):
            with self.subTest(method=method):
                src = emit(method)
                with tempfile.TemporaryDirectory() as tmp:
                    f = os.path.join(tmp, "s.R")
                    with open(f, "w") as fh:
                        fh.write(src)
                    # parse() only — executing needs the data and the full stack.
                    p = subprocess.run(
                        ["Rscript", "-e", f'invisible(parse("{f}"))'],
                        capture_output=True, text=True)
                    self.assertEqual(p.returncode, 0,
                                     f"method={method}: emitted script is not "
                                     f"valid R:\n{p.stderr[-800:]}")


class TestEndToEndReproduction(unittest.TestCase):
    """The real check: run it, re-run the emitted script, diff the numbers.

    Needs limpa + arrow + a DIA-NN report, so it skips in skill-checks CI. It is
    the only test that would have caught either shipped bug on its own — the
    structural checks above encode what these runs taught us.
    """

    def setUp(self):
        if not r_has("limpa", "arrow", "limma", "dplyr", "tidyr"):
            self.skipTest("needs R with limpa/arrow/limma/dplyr/tidyr")
        self.report = os.environ.get("DELIMP_TEST_REPORT")
        if not self.report or not os.path.exists(self.report):
            self.skipTest("set DELIMP_TEST_REPORT to a DIA-NN report.parquet")

    def _run(self, method, tmp):
        import csv
        cond = os.path.join(tmp, "conditions.csv")
        runs = subprocess.run(
            ["Rscript", "-e",
             f'suppressMessages(library(arrow));'
             f'r <- sort(unique(as.data.frame(arrow::open_dataset("{self.report}",'
             f' format="parquet") |> dplyr::select(Run) |> dplyr::collect())$Run));'
             f'cat(r, sep="\\n")'],
            capture_output=True, text=True).stdout.split()
        with open(cond, "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["File.Name", "Group"])
            for i, r in enumerate(runs):
                w.writerow([r, "A" if i % 2 == 0 else "B"])
        out = os.path.join(tmp, f"de_{method}")
        p = subprocess.run(
            ["Rscript", os.path.join(SCRIPTS, "run_de.R"),
             "--input", self.report, "--metadata", cond, "--method", method,
             "--q-cutoff", "0.01", "--contrasts", "A-B", "--outdir", out],
            capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        return out

    @staticmethod
    def _summarise(csv_path):
        import csv
        with open(csv_path) as fh:
            rows = list(csv.DictReader(fh))
        def num(v):
            try: return float(v)
            except (TypeError, ValueError): return None
        return (len(rows),
                sum(1 for r in rows if (num(r.get("adj.P.Val")) or 1) < 0.05),
                {list(r.values())[0]: r.get("logFC") for r in rows})

    def test_emitted_script_reproduces_the_run(self):
        import glob
        for method in sorted(COVERED):
            with self.subTest(method=method):
                with tempfile.TemporaryDirectory() as tmp:
                    out = self._run(method, tmp)
                    rerun = os.path.join(tmp, f"rerun_{method}")
                    os.makedirs(rerun)
                    shutil.copy(os.path.join(out, "reproducibility_log.R"), rerun)
                    p = subprocess.run(["Rscript", "reproducibility_log.R"],
                                       cwd=rerun, capture_output=True, text=True)
                    self.assertEqual(p.returncode, 0, p.stderr[-1500:])
                    a = self._summarise(glob.glob(os.path.join(out, "DE_*.csv"))[0])
                    b = self._summarise(glob.glob(
                        os.path.join(rerun, "de_results_rerun", "DE_*.csv"))[0])
                    # Guard against a hollow pass: 0 == 0 would satisfy every
                    # assertion below while proving nothing ran.
                    self.assertGreater(a[0], 100,
                                       f"{method}: only {a[0]} proteins — the run "
                                       f"itself failed, so the comparison is vacuous")
                    self.assertEqual(a[0], b[0], f"{method}: protein count differs")
                    self.assertEqual(a[1], b[1], f"{method}: significant count differs")
                    diff = [k for k in a[2] if k in b[2] and a[2][k] != b[2][k]]
                    self.assertEqual(diff[:5], [], f"{method}: logFC differs for "
                                                   f"{len(diff)} protein(s)")


if __name__ == "__main__":
    unittest.main(verbosity=2)
