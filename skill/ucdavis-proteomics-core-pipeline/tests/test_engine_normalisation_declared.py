#!/usr/bin/env python3
"""
The Methods name the engine's own cross-run normalisation before the DE's step (2.14.2 C2 review:
statistician M2).

A replicated FragPipe DDA LFQ study: IonQuant normalised across runs (ionquant.normalization=1 in
FragPipe's shipped workflows), then the DE's `auto` between-run rule applies quantile on top. The
Methods said only "quantile normalisation (limma::normalizeBetweenArrays) of the engine's
quantities", leaving out the first normalisation. AlphaDIA with directLFQ is the same. The
default is unchanged (Brett's decision, deferred); only the record and the words come forward:
  - run_search.py's FragPipe and AlphaDIA declarations carry `engine_normalisation`, read from
    the engine's config (ionquant.normalization; AlphaDIA's search_output.normalize_directlfq),
    or NOT RECORDED [confirm] -- never assumed;
  - build_maxlfq.R maxlfq_normalisation_text says it first: "<engine's>, then <step>".
Hermetic: temp folders; the R part needs Rscript.
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import run_search  # noqa: E402
from test_run_de_contaminants import r_has  # noqa: E402
from test_adapters_describe_themselves import FROZEN_2_0_1, FROZEN_2_1_2  # noqa: E402


class Declared(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = self._t.name

    def tearDown(self):
        self._t.cleanup()

    def workflow(self, extra):
        with open(os.path.join(self.d, "fragpipe.workflow"), "w") as fh:
            fh.write("# FragPipe version 23.1\n# IonQuant version 1.11.20\nionquant.mbr=1\n"
                     "phi-report.filter=--sequential --prot 0.01 --picked\n" + extra)

    def test_ionquant_as_its_workflow_says(self):
        for line, words in (("ionquant.normalization=1\n", "IonQuant's cross-run normalisation "
                                                           "(ionquant.normalization=1"),
                            ("ionquant.normalization=0\n", "no cross-run normalisation by "
                                                           "IonQuant"),
                            ("", "whether IonQuant normalised across runs is NOT RECORDED")):
            with self.subTest(line=line):
                self.workflow(line)
                dec = run_search.fragpipe_quantity_declaration(self.d, True)
                self.assertIn(words, dec["engine_normalisation"])
        self.assertIn("[confirm]", dec["engine_normalisation"])

    def test_alphadia_as_its_frozen_config_says(self):
        for frozen, words in ((FROZEN_2_1_2, "directLFQ's between-run normalisation "
                                             "(search_output.normalize_directlfq: true"),
                              (FROZEN_2_0_1, "NOT RECORDED")):
            with self.subTest(words=words):
                with open(os.path.join(self.d, "frozen_config.yaml"), "w") as fh:
                    fh.write(frozen)
                dec = run_search.alphadia_quantity_declaration(self.d, "pg.matrix.parquet")
                self.assertIn(words, dec["engine_normalisation"])


@unittest.skipUnless(r_has("jsonlite"), "needs Rscript")
class TheMethodsWords(unittest.TestCase):
    # build_maxlfq.R finds diann_q_columns.R in --file's folder, else the working directory
    # (its refusal otherwise is deliberate): run from scripts/, as test_sage_de_describes_itself
    # does -- HIVE's runner calls this file from tests/
    R = r'''setwd("{scripts}"); source("build_maxlfq.R")
      cat(maxlfq_normalisation_text(list(method = "{m}", requested = "auto"), "normalised",
          list(level = "protein"{eng})))'''

    def text(self, method, eng=None):
        p = subprocess.run(["Rscript", "-e", self.R.format(
            scripts=SCRIPTS, m=method,
            eng=f', engine_normalisation = "{eng}"' if eng else "")],
            capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        return p.stdout

    def test_the_engines_step_comes_first(self):
        eng = "IonQuant cross-run normalisation (ionquant.normalization=1)"
        self.assertEqual(self.text("quantile", eng),
                         f"{eng}, then quantile normalisation (limma::normalizeBetweenArrays) of "
                         "the engine's quantities")
        self.assertEqual(self.text("none", eng), f"{eng}; no between-run step here")

    def test_a_declaration_without_it_reads_as_before(self):
        """(pin)"""
        self.assertEqual(self.text("quantile"), "quantile normalisation "
                                                "(limma::normalizeBetweenArrays) of the engine's "
                                                "quantities")


if __name__ == "__main__":
    unittest.main()
