#!/usr/bin/env python3
"""
A normalised MaxLFQ DE says what its PG.MaxLFQ carries (2.10 review; CLAUDE.md rule 1).

`run_de.R --method maxlfq --quantities normalised` on a report searched with --no-norm described
itself as "DIA-NN's cross-run normalisation (PG.MaxLFQ), then quantile normalisation" -- a wrong
Methods sentence: nothing but the quantile step normalised it. Now the search's provenance decides
(normalization_check.searched_no_norm: the chain's no_norm_report.parquet name, else
search_provenance.json -- run_search.py's 5-step chain never has --no-norm, a single-shot search
has it when its cfg does), and only when that records nothing are the data read; when neither
can tell, the description claims no DIA-NN normalisation and says NOT RECORDED.

Synthetic reports on a temp dir (tests/test_normalization_check.make_report); the R parts skip
without R. Hermetic.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import diann_parallel as dp  # noqa: E402
import normalization_check as nc  # noqa: E402
from test_normalization_check import HAVE_ARROW, NEEDS, make_report  # noqa: E402
from test_run_de_contaminants import r_has  # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")
QUANTILE_ONLY = ("quantile normalisation (limma::normalizeBetweenArrays) only: the DIA-NN search "
                 "ran with --no-norm, so PG.MaxLFQ carries no DIA-NN cross-run normalisation")
DIANN_THEN_QUANTILE = ("DIA-NN's cross-run normalisation (PG.MaxLFQ), then quantile normalisation "
                       "(limma::normalizeBetweenArrays)")


def provenance(d, **rec):
    with open(os.path.join(d, "search_provenance.json"), "w") as fh:
        json.dump(dict({"engine": "diann"}, **rec), fh)


class FromProvenance(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.report = os.path.join(self.d, "report.parquet")

    def test_the_chains_no_norm_report_name_says_it(self):
        r = nc.searched_no_norm(os.path.join(self.d, dp.NO_NORM_REPORT))
        self.assertIs(r["no_norm"], True)
        self.assertIn("diann_parallel.py --no-norm", r["source"])

    def test_run_searchs_chain_never_has_no_norm(self):
        provenance(self.d, search_mode="parallel_5step")
        self.assertIs(nc.searched_no_norm(self.report)["no_norm"], False)

    def test_a_single_shot_search_has_it_when_its_cfg_does(self):
        cfg = os.path.join(self.d, "params.cfg")
        for text, want in (("--qvalue 0.01\n--no-norm\n", True), ("--qvalue 0.01\n", False),
                           ("--qvalue 0.01\n# --no-norm\n", False)):
            with open(cfg, "w") as fh:
                fh.write(text)
            provenance(self.d, search_mode="single_shot", params_file=cfg)
            r = nc.searched_no_norm(self.report)
            self.assertIs(r["no_norm"], want, text)
            self.assertIn(cfg, r["source"])

    def test_not_recorded_is_none_never_a_guess(self):
        self.assertIsNone(nc.searched_no_norm(self.report)["no_norm"])       # no provenance
        for rec in ({"engine": "sage", "search_mode": "single_shot"},
                    {"search_mode": "single_shot"},                          # no cfg named
                    {"search_mode": "single_shot", "params_file": os.path.join(self.d, "gone.cfg")},
                    {"search_mode": "something_else"}):
            with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
                json.dump(dict({"engine": "diann"}, **rec), fh)
            r = nc.searched_no_norm(self.report)
            self.assertIsNone(r["no_norm"], rec)
            self.assertTrue(r["source"].startswith("not recorded"), r)

    def test_the_search_folder_run_de_found_is_the_one_read(self):
        sub = os.path.join(self.d, "dia-quant-output")
        os.makedirs(sub)
        provenance(self.d, search_mode="parallel_5step")
        self.assertIsNone(nc.searched_no_norm(os.path.join(sub, "report.parquet"))["no_norm"])
        self.assertIs(nc.searched_no_norm(os.path.join(sub, "report.parquet"),
                                          search_dir=self.d)["no_norm"], False)


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class TheDESaysWhatRan(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def de(self, report, cond, contrast, name):
        out = os.path.join(self.d, name)
        p = subprocess.run(["Rscript", RUN_DE, "--input", report, "--metadata", cond,
                            "--contrasts", contrast, "--method", "maxlfq", "--quantities",
                            "normalised", "--outdir", out], capture_output=True, text=True,
                           timeout=900, cwd=self.d)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        with open(os.path.join(out, "methods.txt")) as fh:
            return prov, " ".join(fh.read().split())

    def test_a_no_norm_report_is_never_called_dia_nn_normalised(self):
        report, cond, contrast = make_report(os.path.join(self.d, "nn"), "whole", no_norm=True)
        # no provenance: the data say it
        prov, methods = self.de(report, cond, contrast, "data")
        self.assertEqual(prov["normalisation"], QUANTILE_ONLY)
        self.assertIn(QUANTILE_ONLY, methods)
        self.assertNotIn("DIA-NN's cross-run normalisation (PG.MaxLFQ)", methods)
        self.assertTrue(prov["dia_nn_normalisation"]["no_norm"])
        self.assertTrue(prov["dia_nn_normalisation"]["source"].startswith("read from the data"))
        # the chain's own name: the provenance says it, the data are not read
        named = os.path.join(os.path.dirname(report), dp.NO_NORM_REPORT)
        os.rename(report, named)
        prov, _ = self.de(named, cond, contrast, "named")
        self.assertEqual(prov["normalisation"], QUANTILE_ONLY)
        self.assertIn("the report's name", prov["dia_nn_normalisation"]["source"])

    def test_a_normal_report_says_dia_nn_then_quantile_from_its_provenance(self):
        report, cond, contrast = make_report(os.path.join(self.d, "n"), "whole")
        provenance(os.path.dirname(report), search_mode="parallel_5step")
        prov, methods = self.de(report, cond, contrast, "chain")
        self.assertEqual(prov["normalisation"], DIANN_THEN_QUANTILE)
        self.assertFalse(prov["dia_nn_normalisation"]["no_norm"])
        self.assertIn("search_provenance.json", prov["dia_nn_normalisation"]["source"])
        self.assertIn(DIANN_THEN_QUANTILE, methods)

    def test_neither_provenance_nor_data_claims_no_dia_nn_normalisation(self):
        import pyarrow.parquet as pq
        report, cond, contrast = make_report(os.path.join(self.d, "u"), "whole")
        t = pq.read_table(report)
        pq.write_table(t.drop(["Precursor.Quantity", "Precursor.Normalised"]), report)
        prov, methods = self.de(report, cond, contrast, "unknown")
        self.assertIn("whether DIA-NN also normalised it is NOT RECORDED [confirm]",
                      prov["normalisation"])
        self.assertIsNone(prov["dia_nn_normalisation"]["no_norm"])
        self.assertNotIn("DIA-NN's cross-run normalisation (PG.MaxLFQ)", methods)


if __name__ == "__main__":
    unittest.main()
