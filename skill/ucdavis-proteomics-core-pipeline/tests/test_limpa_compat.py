#!/usr/bin/env python3
"""run_de.R reads the DIA-NN report whatever this limpa calls readDIANN()'s annotation argument.

HIVE sbatch 24154221 (2026-09-27): the 2.8.0 run_de.R passed readDIANN(annotation.columns =),
an argument of limpa >= 1.4.0 (Bioconductor 3.23). setup.sh's env had bioconda's only build,
limpa 1.2.5, whose readDIANN() calls it extra.columns -- so every dpc DE died with "unused
argument". limpa_compat.R now passes the annotation under the name this limpa has, and only a
readDIANN() with neither (none known) joins the columns from the report by Precursor.Id.

Guards (skip without R + limpa/arrow/nanoparquet):
  * the three paths -- annotation.columns (the installed limpa), extra.columns (limpa 1.2.x's
    signature, as a wrapper) and neither (the join) -- give the same precursor matrix and the
    same per-precursor accessions, and each says which path ran;
  * run_de.R records the limpa version and the path in de_provenance.json;
  * the emitted reproducibility script picks the argument name at run time, so it runs on
    either limpa;
  * OPT-IN, with a real limpa 1.2.x: LIMPA_OLD_RLIBS=<R_LIBS holding limpa 1.2.x (+ limma)>
    runs run_de.R on it end to end and compares the tables with the installed limpa's.
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)

from test_run_de_contaminants import SYNTH_R, r_has, rscript   # noqa: E402

COMPAT_R = os.path.join(SCRIPTS, "limpa_compat.R")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")
NEEDS = ("limpa", "limma", "arrow", "dplyr", "tidyr", "jsonlite", "nanoparquet")


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@unittest.skipUnless(r_has(*NEEDS), "needs R with " + "/".join(NEEDS))
class ReadDiannAcrossLimpaVersions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        cls.report = os.path.join(cls.tmp, "report.parquet")
        cls.res = json.loads(rscript(f'''
          source("{COMPAT_R}")
          q <- c("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value")
          ann <- c(LIMPA_ANNOTATION_DEFAULT, "Protein.Ids")
          # limpa 1.2.x's signature, delegating to the installed readDIANN
          old <- function(file = "report.parquet", path = NULL, format = NULL, sep = "\\t",
                          log = TRUE, q.columns = q, q.cutoffs = 0.01,
                          extra.columns = LIMPA_ANNOTATION_DEFAULT)
            limpa::readDIANN(file, format = format, q.columns = q.columns,
                             q.cutoffs = q.cutoffs, annotation.columns = extra.columns)
          # a readDIANN() with neither argument: limpa's default annotation only
          neither <- function(file, format = NULL, q.columns = q, q.cutoffs = 0.01)
            limpa::readDIANN(file, format = format, q.columns = q.columns, q.cutoffs = q.cutoffs)
          runs <- list(native = read_diann_annotated("{cls.report}", "parquet", 0.01, q, ann),
                       extra = read_diann_annotated("{cls.report}", "parquet", 0.01, q, ann, fn = old),
                       join = read_diann_annotated("{cls.report}", "parquet", 0.01, q, ann, fn = neither))
          key <- function(r) paste(rownames(r$dat$E), r$dat$genes$Protein.Ids)
          E <- lapply(runs, function(r) r$dat$E)
          cat(jsonlite::toJSON(list(
            argument = lapply(runs, function(r) r$argument), path = lapply(runs, function(r) r$path),
            args_old = limpa_annotation_arg(old), args_neither = limpa_annotation_arg(neither),
            default_old = limpa_annotation_default(old),
            same_E = c(identical(E$native, E$extra), isTRUE(all.equal(E$native, E$join))),
            same_ids = c(identical(key(runs$native), key(runs$extra)),
                         identical(key(runs$native), key(runs$join))),
            n_ids = sum(!is.na(runs$join$dat$genes$Protein.Ids)), n = nrow(E$join)),
            auto_unbox = TRUE, na = "null"))'''))
        cls.out = os.path.join(cls.tmp, "de")
        cls.proc = subprocess.run(["Rscript", RUN_DE, "--input", cls.report,
                                  "--metadata", os.path.join(cls.tmp, "conditions.csv"),
                                  "--method", "dpc", "--outdir", cls.out],
                                 capture_output=True, text=True, cwd=cls.tmp)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_the_argument_name_is_read_from_the_signature(self):
        self.assertEqual(self.res["args_old"], "extra.columns")
        self.assertIsNone(self.res["args_neither"])
        self.assertEqual(self.res["default_old"],
                         ["Protein.Group", "Protein.Names", "Genes", "Proteotypic"])

    def test_each_path_says_which_ran(self):
        self.assertEqual(self.res["argument"], {"native": "annotation.columns",
                                                "extra": "extra.columns", "join": None})
        self.assertEqual(self.res["path"]["extra"], "readDIANN(extra.columns = ...)")
        self.assertIn("Protein.Ids joined from the report by Precursor.Id", self.res["path"]["join"])

    def test_all_three_read_the_same_precursors_and_accessions(self):
        self.assertEqual(self.res["same_E"], [True, True])
        self.assertEqual(self.res["same_ids"], [True, True])
        self.assertEqual(self.res["n_ids"], self.res["n"])   # every precursor got its own

    def test_run_de_records_the_limpa_and_the_path(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr[-2000:])
        with open(os.path.join(self.out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        lr = prov["limpa_read"]
        self.assertEqual(lr["annotation_argument"], "annotation.columns")
        self.assertEqual(lr["limpa_version"], prov["packages"]["limpa"])
        self.assertIn("readDIANN (limpa ", self.proc.stderr)

    def test_the_repro_script_picks_the_argument_at_run_time(self):
        with open(os.path.join(self.out, "reproducibility_log.R")) as fh:
            src = fh.read()
        self.assertIn("ann_arg <- intersect(c('annotation.columns', 'extra.columns'), "
                      "names(formals(limpa::readDIANN)))[1]", src)
        self.assertIn("setNames(list(ann_cols), ann_arg)", src)
        self.assertNotIn("annotation.columns = c(", src)

    @unittest.skipUnless(os.environ.get("LIMPA_OLD_RLIBS"),
                         "set LIMPA_OLD_RLIBS to an R_LIBS path holding limpa 1.2.x (+ its limma)")
    def test_a_real_limpa_1_2_runs_end_to_end(self):
        out = os.path.join(self.tmp, "de_old")
        env = dict(os.environ, R_LIBS=os.environ["LIMPA_OLD_RLIBS"])
        p = subprocess.run(["Rscript", RUN_DE, "--input", self.report,
                            "--metadata", os.path.join(self.tmp, "conditions.csv"),
                            "--method", "dpc", "--outdir", out],
                           capture_output=True, text=True, cwd=self.tmp, env=env)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            lr = json.load(fh)["limpa_read"]
        self.assertTrue(lr["limpa_version"].startswith("1.2"), lr)
        self.assertEqual(lr["annotation_argument"], "extra.columns")
        a = {r["Protein.Group"]: r for r in read_csv(os.path.join(self.out, "DE_dpc_B.A.csv"))}
        b = {r["Protein.Group"]: r for r in read_csv(os.path.join(out, "DE_dpc_B.A.csv"))}
        self.assertEqual(set(a), set(b))
        for k in a:
            for col in ("logFC", "P.Value"):
                self.assertAlmostEqual(float(a[k][col]), float(b[k][col]), places=8)


if __name__ == "__main__":
    unittest.main()
