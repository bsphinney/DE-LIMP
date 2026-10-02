#!/usr/bin/env python3
"""run_de.R blocks technical replicates on the Sample column (2.10 review, MED 7).

core_submission.py locate --reinjections all keeps every injection of a sample. Before this fix
the conditions offered no column naming each run's biological sample, so two injections of one
sample were fitted as two independent samples -- n doubled, p-values too small. Now
`conditions` writes the Sample column (collect_conditions.TECH_REPLICATE_COLUMN) and run_de.R,
seeing a Sample shared by several runs, fits it as a random blocking factor
(duplicateCorrelation) and reports every contrast from that fit; a sample injected once is a
block of one. The Methods say which samples this applied to.

End to end on a synthetic DIA-NN report; skipped without R + limpa/limma/statmod/arrow/...
NOTE: CI runs Python only and skips every R-invoking test -- run them locally.
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
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import collect_conditions as cc  # noqa: E402
import make_methods as mm  # noqa: E402
from test_run_de_block import r_has, rscript  # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")

# 6 biological samples, 3 per group; C1, C3, T1 and T3 injected twice (10 runs). Each sample has
# its own level (sd 0.5) and each injection a little noise (sd 0.1); P00001-P00020 go up 0.6.
SYNTH_R = r'''
set.seed(5)
samples <- c("C1", "C1", "C2", "C3", "C3", "T1", "T1", "T2", "T3", "T3")
grp <- ifelse(substr(samples, 1, 1) == "C", "ctrl", "treat")
runs <- sprintf("run%02d", seq_along(samples))
out <- list(); pgq <- list()
for (i in 1:120) {
  pg <- sprintf("P%05d", i); base <- runif(1, 12, 20)
  s_eff <- setNames(rnorm(6, 0, 0.5), c("C1", "C2", "C3", "T1", "T2", "T3"))
  eff <- if (i <= 20) 0.6 else 0
  run_lv <- base + s_eff[samples] + (grp == "treat") * eff + rnorm(length(runs), 0, 0.1)
  pgq[[i]] <- data.frame(Protein.Group = pg, Run = runs, PG.MaxLFQ = 2^run_lv)
  for (j in 1:3) {
    lv <- run_lv + rnorm(1, 0, 1)
    int <- 2^(lv + rnorm(length(runs), 0, 0.1))
    out[[length(out) + 1]] <- data.frame(Run = runs, Precursor.Id = sprintf("PEP%dK%d2", i, j),
      Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0("G", i, "_X"),
      Genes = paste0("G", i), Proteotypic = 1L, Precursor.Normalised = int,
      Precursor.Quantity = int, Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001,
      PG.Q.Value = 0.001, Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001,
      stringsAsFactors = FALSE)
  }
}
d <- merge(do.call(rbind, out), do.call(rbind, pgq), by = c("Protein.Group", "Run"))
arrow::write_parquet(d, file.path(OUT, "report.parquet"))
meta <- data.frame(File.Name = runs, Group = grp, Sample = samples, Other = samples)
write.csv(meta, file.path(OUT, "conditions.csv"), row.names = FALSE)
m2 <- meta; m2$Sample[10] <- "C1"          # C1 now has a run in the treat group
write.csv(m2, file.path(OUT, "conditions_span.csv"), row.names = FALSE)
'''


def collect_header_matches():
    """run_de.R's column name is collect_conditions.TECH_REPLICATE_COLUMN: one name, two files."""
    with open(RUN_DE) as fh:
        return f'TECH_REPLICATE_COLUMN <- "{cc.TECH_REPLICATE_COLUMN}"' in fh.read()


class OneColumnName(unittest.TestCase):
    def test_run_de_and_collect_conditions_name_the_same_column(self):
        self.assertTrue(collect_header_matches())


@unittest.skipUnless(r_has("limpa", "limma", "statmod", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limpa/limma/statmod/arrow/dplyr/tidyr/jsonlite")
class TechnicalReplicatesAreBlocked(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        cls.runs = {}
        for key, meta, args in (("dpc", "conditions.csv", ["--method", "dpc"]),
                                ("ml", "conditions.csv", ["--method", "maxlfq"]),
                                ("other_block", "conditions.csv", ["--method", "maxlfq",
                                                                   "--block", "Other"]),
                                ("scope", "conditions.csv", ["--method", "maxlfq", "--block",
                                                             "Sample", "--block-scope", "within"]),
                                ("span", "conditions_span.csv", ["--method", "maxlfq"])):
            out = os.path.join(cls.tmp, "out_" + key)
            p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(cls.tmp, "report.parquet"),
                                "--metadata", os.path.join(cls.tmp, meta), "--outdir", out, *args],
                               capture_output=True, text=True, cwd=cls.tmp)
            cls.runs[key] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def prov(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            return json.load(fh)

    def test_the_sample_is_a_random_block_for_every_contrast(self):
        for key in ("dpc", "ml"):
            b = self.prov(key)["block"]
            self.assertTrue(b["applied"], key)
            self.assertEqual((b["column"], b["effect"]), ("Sample", "random"), key)
            self.assertEqual(set(b["contrast_model"].values()), {"blocked"}, key)
            tr = b["technical_replicates"]
            self.assertEqual((tr["n_samples"], tr["n_runs"]), (4, 8), key)
            self.assertEqual(tr["samples"], {"C1": 2, "C3": 2, "T1": 2, "T3": 2})

    def test_methods_say_so(self):
        _p, out = self.runs["dpc"]
        with open(os.path.join(out, "methods.txt")) as fh:
            text = fh.read()
        self.assertIn("Replicates    : 4 sample(s) were injected more than once (8 runs; Sample "
                      "column): technical replicates, blocked on Sample below, never counted as "
                      "independent samples", text)
        sent = mm.de_block_sentence(self.prov("dpc"))
        self.assertTrue(sent.startswith("4 sample(s) were injected more than once (8 runs); these "
                                        "technical replicates were not counted as independent "
                                        "samples. Samples sharing a Sample were modelled as "
                                        "correlated"), sent)

    def test_another_block_or_scope_or_a_sample_in_two_groups_stops(self):
        for key, why in (("other_block", "cannot both be blocked"),
                         ("scope", "--block-scope within cannot apply"),
                         ("span", "has runs in more than one Group")):
            p, _out = self.runs[key]
            self.assertNotEqual(p.returncode, 0, key)
            self.assertIn(why, p.stderr, key)

    def test_the_conditions_csv_carries_the_column(self):
        with open(os.path.join(self.tmp, "conditions.csv"), newline="") as fh:
            self.assertIn("Sample", next(csv.reader(fh)))


if __name__ == "__main__":
    unittest.main()
