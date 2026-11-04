#!/usr/bin/env python3
"""
A FragPipe DE's Methods state the Philosopher filter as it ran, and over which runs its rates
hold (2.14.2 C2 review: proteomics M1).

For a quantified FragPipe DDA DE (LFQ-MBR, then run_de.R --method maxlfq), the DE paragraph quoted
fragpipe_quantity_declaration's identification_fdr, built from the workflow's phi-report.filter
line -- but FragPipe adds --razor when it runs the filter, and with LFQ on it can change --prot --
while the search paragraph's own FDR sentence (the filter as run, from FragPipe's log, with its
scope) was dropped because the adapted report's placeholder q-columns set de_prov.q_columns. Now
fragpipe_params.fdr_statement is the one statement: the declaration and the Methods both read it,
and the search paragraph keeps it when the DE's q-columns are placeholders.
Hermetic: a temp FragPipe output folder, in-process.
"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fragpipe_params as fpp  # noqa: E402
import make_methods  # noqa: E402
import run_search  # noqa: E402

WORKFLOW_FILTER = "--sequential --prot 0.01 --picked"
LOG = """FragPipe log
PhilosopherFilter [Work dir: {d}/exp_runA]
/opt/philosopher/philosopher filter --sequential --prot 0.01 --picked --razor --pepxml {d}/exp_runA --protxml {d}/combined.prot.xml --tag rev_
PhilosopherFilter [Work dir: {d}/exp_runB]
/opt/philosopher/philosopher filter --sequential --prot 0.01 --picked --razor --pepxml {d}/exp_runB --protxml {d}/combined.prot.xml --tag rev_
"""


class FilterAsRun(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = self._t.name
        with open(os.path.join(self.d, fpp.SAVED_WORKFLOW), "w") as fh:
            fh.write("# FragPipe version 23.1\n# IonQuant version 1.11.20\n"
                     f"{fpp.PHILOSOPHER_FILTER_KEY}={WORKFLOW_FILTER}\n"
                     "phi-report.run-report=true\nquantitation.run-label-free-quant=true\n"
                     "ionquant.mbr=1\n")
        with open(os.path.join(self.d, "fragpipe.fp-manifest"), "w") as fh:
            fh.write(f"{self.d}/runA.d\trunA\t\tDDA\n{self.d}/runB.d\trunB\t\tDDA\n")

    def tearDown(self):
        self._t.cleanup()

    def log(self):
        with open(os.path.join(self.d, "log_2026-01-01_00-00-00.txt"), "w") as fh:
            fh.write(LOG.format(d=self.d))

    def test_the_declaration_quotes_the_filter_as_run_with_its_scope(self):
        self.log()
        ident = run_search.fragpipe_quantity_declaration(self.d, True)["identification"]
        self.assertIn("--razor", ident)
        self.assertIn("FragPipe's log, the filter as run", ident)
        self.assertNotIn("(fragpipe.workflow)", ident)
        self.assertIn("within each run", ident)
        self.assertIn("over all runs together", ident)          # the protein FDR's scope
        self.assertTrue(ident.endswith(run_search.PLACEHOLDER_Q), ident)

    def test_without_a_log_the_workflow_line_is_tagged(self):
        ident = run_search.fragpipe_quantity_declaration(self.d, True)["identification"]
        self.assertIn(f"Philosopher filter {WORKFLOW_FILTER} (fragpipe.workflow)", ident)
        self.assertIn("This is the workflow's filter line", ident)
        self.assertIn("confirm", ident)

    def test_one_statement_for_the_de_and_the_search_paragraph(self):
        self.log()
        st = fpp.fdr_statement(self.d, WORKFLOW_FILTER)
        self.assertTrue(st["from_log"])
        self.assertIn(st["sentence"],
                      run_search.fragpipe_quantity_declaration(self.d, True)["identification"])
        with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
            fh.write('{"engine": "fragpipe", "version": "23.1"}')
        rec = make_methods.search_record(
            params=os.path.join(self.d, fpp.SAVED_WORKFLOW),
            search_prov=os.path.join(self.d, "search_provenance.json"))
        self.assertEqual(rec["fdr_sentence"], st["sentence"])
        # the DE's q-columns are 0.0 placeholders: the search paragraph keeps the filter as run
        para = make_methods.search_paragraph(rec, {"q_columns": ["Q.Value"], "q_cutoffs": [0.01],
                                                   "q_columns_role": "placeholder"})
        self.assertIn(st["sentence"], para)
        # ...and a DE whose q-columns really filtered states its own, as before
        para = make_methods.search_paragraph(rec, {"q_columns": ["Q.Value"], "q_cutoffs": [0.01],
                                                   "q_columns_role": "real"})
        self.assertNotIn(st["sentence"], para)

    def test_a_fragpipe_dia_de_reads_as_before(self):
        """(C2 verify, proteomics nit) FragPipe's DIA route: the DE reads DIA-NN's report, whose
        q-values are real and whose record names no q_columns_role. The Philosopher filter only
        built the library there, so the search paragraph does not print it beside DIA-NN's
        q-values, where it would read as the final filter."""
        self.log()
        with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
            fh.write('{"engine": "fragpipe", "version": "23.1"}')
        rec = make_methods.search_record(
            params=os.path.join(self.d, fpp.SAVED_WORKFLOW),
            search_prov=os.path.join(self.d, "search_provenance.json"))
        st = fpp.fdr_statement(self.d, WORKFLOW_FILTER)
        dia_de = {"q_columns": ["Q.Value", "Lib.Q.Value"], "q_cutoffs": [0.01, 0.01],
                  "pipeline_id": "dpc"}                   # no q_columns_role: DIA-NN's own
        para = make_methods.search_paragraph(rec, dia_de)
        self.assertNotIn(st["sentence"], para)
        self.assertNotIn("Philosopher", para)


if __name__ == "__main__":
    unittest.main()
