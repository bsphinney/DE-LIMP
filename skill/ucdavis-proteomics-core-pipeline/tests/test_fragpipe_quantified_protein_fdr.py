#!/usr/bin/env python3
"""A QUANTIFIED FragPipe search with no protein-level FDR (C2 proteomics H1, skill 2.15).

FragPipe's Nonspecific-peptidome and -HLA workflows ship Philosopher's filter
`--sequential --prot 1`: no protein-level FDR. make_presets switches label-free quantification on
in them when amounts are wanted, so adapt_fragpipe_dda writes every protein group of
combined_protein.tsv into report.parquet. Before this fix the quantified adapter recorded no
`protein_fdr` (only the identifications-only adapter did), and the analysis brief described the
report as "already FDR-filtered by FragPipe" -- the protein level was not filtered at all.

Now:
  * adapt_fragpipe_dda records `protein_fdr` / `protein_filter` (the one reading,
    run_search.fragpipe_filter_as_run -> fragpipe_params.protein_fdr_record) in
    fragpipe_quantity.json and search_provenance.json `fragpipe_quantity`;
  * analysis_prompt.report_parquet_what states the protein level from that record, never
    "already FDR-filtered";
  * the DE brief and the no-DE brief ask for no protein-group count as a result;
  * the report page carries a fixed warning callout (make_analysis_html.protein_fdr_note);
  * with a protein FDR in the filter (LFQ-MBR's --prot 0.01) none of this fires.
One predicate: fragpipe_params.no_protein_fdr. Hermetic; stdlib + pyarrow.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import fragpipe_params as fpp                   # noqa: E402
import fragpipe_quant as fq                     # noqa: E402
from test_fragpipe_dda_quantity import write_output   # noqa: E402

try:
    import pyarrow  # noqa: F401
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

NO_PROT = "--sequential --prot 1"
WITH_PROT = "--sequential --prot 0.01 --picked"
RUNS = ["r-1", "r-2", "r-3"]


def set_filter(d, filt):
    """fragpipe.workflow with this Philosopher filter line (no FragPipe log: the workflow's line
    is the filter read, as fragpipe_params.filter_as_recorded says)."""
    path = os.path.join(d, fpp.SAVED_WORKFLOW)
    with open(path) as fh:
        lines = [ln for ln in fh if not ln.startswith(fpp.PHILOSOPHER_FILTER_KEY + "=")]
    with open(path, "w") as fh:
        fh.write("".join(lines) + f"{fpp.PHILOSOPHER_FILTER_KEY}={filt}\n")


class Predicate(unittest.TestCase):
    def test_the_one_predicate(self):
        self.assertTrue(fpp.no_protein_fdr(fpp.protein_fdr_record(NO_PROT)))
        self.assertTrue(fpp.no_protein_fdr(fpp.protein_fdr_record("")))      # not recorded
        self.assertFalse(fpp.no_protein_fdr(fpp.protein_fdr_record(WITH_PROT)))
        self.assertFalse(fpp.no_protein_fdr(None))                           # older record
        self.assertFalse(fpp.no_protein_fdr({}))

    def test_recorded_protein_fdr_takes_the_first_record_that_has_one(self):
        p = fpp.protein_fdr_record(NO_PROT)
        self.assertIs(fpp.recorded_protein_fdr(None, {"x": 1}, {"protein_fdr": p}), p)
        self.assertIsNone(fpp.recorded_protein_fdr(None, {}))


@unittest.skipUnless(HAVE_ARROW, "pyarrow")
class _Search(unittest.TestCase):
    FILTER = NO_PROT

    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.addCleanup(self._t.cleanup)
        env = mock.patch.dict(os.environ, {fq.VENV_ENV: os.path.join(self._t.name, "no_venv")})
        env.start()
        self.addCleanup(env.stop)
        self.s = os.path.join(self._t.name, "session")
        self.d = os.path.join(self.s, "output", "search_fragpipe")
        os.makedirs(self.d)
        write_output(self.d, RUNS)                    # 3 experiments: IonQuant MaxLFQ
        set_filter(self.d, self.FILTER)
        import run_search
        self.report = run_search.adapt_fragpipe_dda(self.d)
        with open(os.path.join(self.d, fq.RECORD)) as fh:
            self.qrec = json.load(fh)
        with open(os.path.join(self.d, "search_provenance.json")) as fh:
            self.sp = json.load(fh)


class NoProteinFdr(_Search):
    def test_the_quantity_record_says_no_protein_fdr(self):
        self.assertIs(self.qrec["protein_fdr"]["applied"], False)
        self.assertIn("not FDR-controlled", self.qrec["protein_filter"])
        self.assertIs(self.sp["fragpipe_quantity"]["protein_fdr"]["applied"], False)

    def test_the_file_list_never_says_already_fdr_filtered(self):
        import analysis_prompt
        what = analysis_prompt.report_parquet_what(self.sp)
        self.assertNotIn("already FDR-filtered", what)
        self.assertIn("NO protein-level FDR", what)
        self.assertIn("--prot 1", what)

    def test_the_page_warns(self):
        import make_analysis_html as mah
        note = mah.protein_fdr_note({"input": self.report})
        self.assertIsNotNone(note)
        self.assertEqual(note["kind"], "warning")
        self.assertIn("not a controlled result", note["text"])

    def _brief(self, de):
        tables = os.path.join(self.s, "output", "tables")
        os.makedirs(tables, exist_ok=True)
        cond = os.path.join(self.s, "input", "conditions.csv")
        os.makedirs(os.path.dirname(cond), exist_ok=True)
        with open(cond, "w") as fh:
            fh.write("File.Name,Group\nr-1,a\nr-2,b\nr-3,b\n")
        if de:
            with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
                json.dump({"contrasts": ["b-a"], "adjp": 0.05, "input": self.report,
                           "groups": {"a": ["r-1"], "b": ["r-2", "r-3"]}}, fh)
        out = os.path.join(self._t.name, "brief.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                            "--de-dir", tables, "--out", out, "--conditions", cond,
                            "--session", self.s, "--report", self.report,
                            "--acquisition", "DDA"],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out, encoding="utf-8") as fh:
            return " ".join(fh.read().split())

    def test_the_de_brief_asks_for_no_protein_count(self):
        b = self._brief(de=True)
        self.assertIn("No protein-level FDR (say this in Data Quality Notes)", b)
        self.assertIn("Never give the number of protein groups", b)

    def test_the_no_de_brief_asks_for_no_protein_count(self):
        b = self._brief(de=False)
        self.assertNotIn("already FDR-filtered", b)
        self.assertIn("never a protein-group count", b)
        self.assertNotIn("distinct peptides and protein groups", b)


class WithProteinFdr(_Search):
    FILTER = WITH_PROT

    def test_a_protein_fdr_is_recorded_and_nothing_warns(self):
        import analysis_prompt
        import make_analysis_html as mah
        self.assertIs(self.qrec["protein_fdr"]["applied"], True)
        what = analysis_prompt.report_parquet_what(self.sp)
        self.assertIn("protein level included", what)
        self.assertNotIn("NO protein-level FDR", what)
        self.assertIsNone(mah.protein_fdr_note({"input": self.report}))


if __name__ == "__main__":
    unittest.main()
