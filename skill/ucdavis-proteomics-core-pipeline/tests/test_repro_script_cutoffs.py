#!/usr/bin/env python3
"""
reproducibility_log.R must describe the analysis that actually ran.

Caught by auditing the skill's own output. Once PG.Q.Value took DIA-NN's
recommended 0.05 while the other five kept --q-cutoff, the emitter still wrote

    limpa::readDIANN(report, q.cutoffs = 0.01, q.columns = c(...six...))

a SCALAR. limpa recycles q.cutoffs across q.columns, so the emitted script ran
clean and produced different results from the run it claims to reproduce:

    original run   4,648 proteins, 1,859 significant
    emitted script 4,624 proteins, 1,846 significant

Nothing errored. The CSVs beside it were right. Only the script disagreed --
DE-LIMP architectural rule #1 ("exports must describe what actually ran").

Needs Rscript but NOT limpa: the emitter only builds strings.
"""
import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")


def have_rscript():
    return shutil.which("Rscript") is not None


def emit(q_columns, q_cutoffs, q_cutoff=0.01):
    """Run write_repro_script() and return the generated R source."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "reproducibility_log.R")
        cols = ", ".join(f'"{c}"' for c in q_columns)
        cuts = "NULL" if q_cutoffs is None else "c(" + ", ".join(str(x) for x in q_cutoffs) + ")"
        r = f'''
          source("{SCRIPTS}/diann_q_columns.R"); source("{SCRIPTS}/repro_script.R")
          meta <- data.frame(File.Name = c("r1","r2","r3","r4"),
                             Group = c("A","A","B","B"), stringsAsFactors = FALSE)
          invisible(write_repro_script(
            path = "{out}", method = "dpc", input = "report.parquet",
            format = "parquet", q_cutoff = {q_cutoff},
            q_columns = c({cols}), q_cutoffs = {cuts},
            eq_cutoff = 0, pgq_cutoff = 0, cov_min_frac = 0,
            meta = meta, covariates = character(0), formula_parts = character(0),
            forms = "A-B", adjp_thr = 0.05, logfc_ref = 1))
          cat(readLines("{out}"), sep = "\\n")
        '''
        p = subprocess.run(["Rscript", "-e", r], capture_output=True, text=True)
        if p.returncode != 0:
            raise AssertionError(p.stderr[-2000:])
        return p.stdout


ALL_SIX = ["Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
           "PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value"]


class TestEmittedCutoffs(unittest.TestCase):
    def setUp(self):
        if not have_rscript():
            self.skipTest("Rscript not available")

    def test_differing_cutoffs_are_emitted_as_a_named_vector(self):
        src = emit(ALL_SIX, [0.01, 0.01, 0.01, 0.05, 0.01, 0.01])
        self.assertIn("'PG.Q.Value' = 0.05", src)
        # The exact shape of the old bug: a bare scalar alongside six columns.
        self.assertNotIn("q.cutoffs = 0.01,", src)

    def test_every_column_appears_in_the_emitted_cutoff_vector(self):
        src = emit(ALL_SIX, [0.01, 0.01, 0.01, 0.05, 0.01, 0.01])
        for c in ALL_SIX:
            self.assertIn(f"'{c}' = ", src)

    def test_uniform_cutoffs_still_emit_a_scalar(self):
        # Don't make the common case ugly: one value for all six stays scalar.
        src = emit(ALL_SIX, [0.01] * 6)
        self.assertIn("q.cutoffs = 0.01", src)
        self.assertNotIn("'PG.Q.Value' = ", src)

    def test_cutoffs_default_to_the_shared_map_when_not_passed(self):
        # A caller that omits q_cutoffs must still get the real cutoffs, not the
        # scalar -- otherwise the old bug returns through the default path.
        src = emit(ALL_SIX, None)
        self.assertIn("'PG.Q.Value' = 0.05", src)

    def test_emitted_script_reads_the_recorded_input(self):
        src = emit(ALL_SIX, [0.01, 0.01, 0.01, 0.05, 0.01, 0.01])
        self.assertIn("report.parquet", src)
        self.assertIn("limpa::readDIANN", src)


if __name__ == "__main__":
    unittest.main(verbosity=2)
