#!/usr/bin/env python3
"""run_de.R --method maxlfq lists the proteins detected in one group and absent from the other
(OnOff_Proteins.csv), 2.15.x. Brett, 2026-11-13: "make sure we find these." (An islet-organoid
client cares most: ATF6 seen with isotype, not with mAb43, under cytokines.)

Before: such a protein either fell below --coverage-min and was dropped without a trace, or passed
it, was all-NA on one side and came back with NA logFC, never listed. Now, on the MaxLFQ route only:
  * counted on the matrix BEFORE the coverage filter (so a coverage-dropped protein still appears);
  * simple A-B contrasts only; an average / interaction contrast is skipped WITH its reason in
    de_provenance.json `onoff.skipped`, never silently;
  * --onoff-min-n (default 2) = runs needed in the group that has it;
  * de_provenance.json `onoff` (rule, n_min, per-contrast counts, skipped + why) feeds methods.txt,
    make_methods, the analysis page note and the AI brief -- none carries its own wording;
  * session.py finalize records it in MANIFEST.txt ([OK] / [SKIPPED] with the reason).
dpc (limpa infers missing values) is untouched: no file, `onoff` null, the flag refused.

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

import analysis_prompt as ap  # noqa: E402,F401
import make_analysis_html as mah  # noqa: E402
import make_methods as mm  # noqa: E402
from test_run_de_contaminants import r_has, rscript  # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")

# 3 groups x 3 runs (9 runs; --coverage-min 0.5 keeps a protein seen in >= 5 of 9).
#   P00001-P00050  every run (P00001-P00005 up 1.5 log2 in B)
#   P00051-P00055  in B and C, not A: 6/9 runs -> PASSES coverage; B - A and C - A on/off
#   P00061-P00063  in B only: 3/9 runs -> kept by the per-group coverage rule (2.15.1; dropped
#                  before); B - A on/off
#   P00071         in 1 run of B only: below n_min 2; on/off at --onoff-min-n 1 (and dropped by coverage)
#   P00081-P00082  in A only: 3/9 -> kept (per group); A is the SECOND group of B - A, so "A_only" there
# Genes: G<n> except P00061, whose Genes is empty (NA gene must not drop the protein).
SYNTH_R = r'''
set.seed(11)
runs <- sprintf("run%02d", 1:9); grp <- rep(c("A", "B", "C"), each = 3)
rows <- list()
add <- function(pg, prec, base, only = NULL, only_run = NULL, effB = 0, gene = NULL) {
  for (r in seq_along(runs)) {
    if (!is.null(only) && !(grp[r] %in% only)) next
    if (!is.null(only_run) && r != only_run) next
    lv <- base + rnorm(1, 0, 0.2) + (grp[r] == "B") * effB
    int <- 2^(lv + rnorm(1, 0, 0.3))
    rows[[length(rows) + 1]] <<- data.frame(Run = runs[r], Precursor.Id = prec,
      Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0(pg, "_X"),
      Genes = if (is.null(gene)) sub("P0*", "G", pg) else gene,
      Proteotypic = 1L, Precursor.Normalised = int, Precursor.Quantity = int,
      Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001,
      Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
  }
}
spec <- function(i) {
  if (i >= 51 && i <= 55) list(only = c("B", "C")) else
  if (i >= 61 && i <= 63) list(only = "B") else
  if (i == 71) list(only_run = 4) else
  if (i >= 81 && i <= 82) list(only = "A") else list()
}
for (i in c(1:55, 61:63, 71, 81:82)) {
  pg <- sprintf("P%05d", i); s <- spec(i)
  base <- runif(1, 17, 20)
  for (k in 1:3)
    add(pg, sprintf("PEP%dK%d2", i, k), base + rnorm(1, 0, 0.5), only = s[["only"]], only_run = s[["only_run"]],
        effB = if (i <= 5) 1.5 else 0, gene = if (i == 61) "" else NULL)
}
d <- do.call(rbind, rows)
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum)
names(pgm)[3] <- "PG.MaxLFQ"
d <- merge(d, pgm, by = c("Protein.Group", "Run"))
arrow::write_parquet(d, file.path(OUT, "report.parquet"))
write.csv(data.frame(File.Name = runs, Group = grp), file.path(OUT, "conditions.csv"), row.names = FALSE)
'''

PASSES_COVERAGE = {f"P{i:05d}" for i in range(51, 56)}
COVERAGE_DROPPED_B = {f"P{i:05d}" for i in (61, 62, 63)}
A_ONLY = {"P00081", "P00082"}


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@unittest.skipUnless(r_has("limma", "arrow", "dplyr", "tidyr", "jsonlite", "limpa", "statmod"),
                     "needs R with limma/limpa/statmod/arrow/dplyr/tidyr/jsonlite")
class OnOffMaxlfq(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        cls.runs = {}
        for key, method, extra in (
                ("default", "maxlfq", []),                                        # B-A, C-A
                ("n1", "maxlfq", ["--onoff-min-n", "1"]),
                ("complex", "maxlfq", ["--contrasts", "B-A,C-(A+B)/2"]),
                ("cov0", "maxlfq", ["--coverage-min", "0.1"]),
                ("dpc", "dpc", [])):
            out = os.path.join(cls.tmp, "out_" + key)
            p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(cls.tmp, "report.parquet"),
                                "--metadata", os.path.join(cls.tmp, "conditions.csv"),
                                "--method", method, "--outdir", out] + extra,
                               capture_output=True, text=True, cwd=cls.tmp)
            cls.runs[key] = (p, out)
        cls.refused = subprocess.run(
            ["Rscript", RUN_DE, "--input", os.path.join(cls.tmp, "report.parquet"),
             "--metadata", os.path.join(cls.tmp, "conditions.csv"), "--method", "dpc",
             "--onoff-min-n", "2", "--outdir", os.path.join(cls.tmp, "out_refused")],
            capture_output=True, text=True, cwd=cls.tmp)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-3000:])
        return out

    def onoff(self, key):
        return read_csv(os.path.join(self.out(key), "OnOff_Proteins.csv"))

    def prov(self, key):
        with open(os.path.join(self.out(key), "de_provenance.json")) as fh:
            return json.load(fh)

    def test_file_written_with_the_columns_and_the_default_n_min_2(self):
        rows = self.onoff("default")
        self.assertEqual(list(rows[0]), ["Protein.Group", "Genes", "Contrast", "Direction",
                                         "detected_g1", "detected_g2", "total_g1", "total_g2"])
        got = {(r["Protein.Group"], r["Contrast"], r["Direction"]) for r in rows}
        # B - A and C - A; B is the first group of "B - A"
        want = {(pg, "B - A", "B_only") for pg in PASSES_COVERAGE | COVERAGE_DROPPED_B}
        want |= {(pg, "B - A", "A_only") for pg in A_ONLY}
        want |= {(pg, "C - A", "C_only") for pg in PASSES_COVERAGE}
        want |= {(pg, "C - A", "A_only") for pg in A_ONLY}
        self.assertEqual(got, want)
        r = next(x for x in rows if x["Protein.Group"] == "P00061" and x["Contrast"] == "B - A")
        self.assertEqual((r["detected_g1"], r["detected_g2"], r["total_g1"], r["total_g2"]),
                         ("3", "0", "3", "3"))
        # one run only (P00071) is below n_min 2
        self.assertNotIn("P00071", {x["Protein.Group"] for x in rows})
        # the other direction's protein is never listed under the wrong group
        self.assertTrue(all(x["Protein.Group"] not in A_ONLY or x["Direction"] == "A_only" for x in rows))

    def test_a_group_specific_protein_now_reaches_limma(self):
        """2.15.1 (audit L2-1): the coverage rule counts any one group too, so a protein in every
        run of B and none of A (3 of 9 runs) is no longer deleted: it is in the DE table (no
        fold change for B-A, Evidence 'presence call') AND listed as on/off. Only a protein in
        fewer than 2 runs of every group is still dropped (P00071)."""
        out = self.out("default")
        expr = {r["Protein.Group"] for r in read_csv(os.path.join(out, "Expression_Matrix.csv"))}
        self.assertTrue((COVERAGE_DROPPED_B | A_ONLY | PASSES_COVERAGE) <= expr)
        self.assertNotIn("P00071", expr)
        listed = {r["Protein.Group"] for r in self.onoff("default")}
        self.assertTrue(COVERAGE_DROPPED_B <= listed)
        t = self.prov("default")["onoff"]["contrasts"]["B-A"]
        self.assertEqual(t["n_dropped_by_coverage_filter"], 0)
        de = {r["Protein.Group"]: r for r in read_csv(os.path.join(out, "DE_maxlfq_B.A.csv"))}
        for pg in COVERAGE_DROPPED_B | A_ONLY:
            self.assertIn(pg, de)
            self.assertIn(de[pg]["logFC"], ("", "NA"))
            self.assertEqual(de[pg]["Evidence"], "presence call")
        cov = self.prov("default")["coverage"]
        self.assertEqual(cov["min_all"], 5)
        self.assertEqual(cov["min_group"], {"A": 2, "B": 2, "C": 2})
        self.assertEqual(cov["n_kept_by_group_only"], len(COVERAGE_DROPPED_B | A_ONLY))
        self.assertEqual(cov["n_dropped"], 1)                         # P00071

    def test_na_gene_keeps_the_protein(self):
        rows = {(r["Protein.Group"], r["Contrast"]): r for r in self.onoff("default")}
        self.assertIn(("P00061", "B - A"), rows)
        self.assertIn(rows[("P00061", "B - A")]["Genes"], ("", "NA"))
        self.assertEqual(rows[("P00062", "B - A")]["Genes"], "G62")

    def test_onoff_min_n_1_adds_the_single_run_protein(self):
        rows = self.onoff("n1")
        self.assertIn(("P00071", "B - A", "B_only"),
                      {(r["Protein.Group"], r["Contrast"], r["Direction"]) for r in rows})
        rec = self.prov("n1")["onoff"]
        self.assertEqual(rec["n_min"], 1)
        self.assertIn("at least 1 run ", rec["rule"])

    def test_provenance_record(self):
        rec = self.prov("default")["onoff"]
        self.assertEqual(rec["status"], "written")
        self.assertEqual((rec["file"], rec["n_min"]), ("OnOff_Proteins.csv", 2))
        for w in ("not a measured ratio", "before the coverage filter", "at least 2 runs"):
            self.assertIn(w, rec["rule"])
        self.assertEqual(set(rec["contrasts"]), {"B-A", "C-A"})
        ba = rec["contrasts"]["B-A"]
        self.assertEqual((ba["group_1"], ba["group_2"]), ("B", "A"))
        self.assertEqual((ba["n_only_group_1"], ba["n_only_group_2"], ba["n_onoff"]),
                         (len(PASSES_COVERAGE | COVERAGE_DROPPED_B), len(A_ONLY),
                          len(PASSES_COVERAGE | COVERAGE_DROPPED_B | A_ONLY)))
        self.assertEqual(rec["n_onoff_total"], len(self.onoff("default")))
        self.assertFalse(rec["skipped"])      # R writes an empty list as []
        # methods.txt says it, from the same record
        with open(os.path.join(self.out("default"), "methods.txt")) as fh:
            txt = fh.read()
        self.assertIn("OnOff_Proteins.csv", txt)
        self.assertIn("not as fold changes", txt)

    def test_complex_contrast_is_skipped_with_the_reason_not_silently(self):
        rec = self.prov("complex")["onoff"]
        self.assertEqual(set(rec["contrasts"]), {"B-A"})
        self.assertEqual(set(rec["skipped"]), {"C-(A+B)/2"})
        self.assertIn("not a simple A-B contrast", rec["skipped"]["C-(A+B)/2"])
        self.assertIn("3 terms", rec["skipped"]["C-(A+B)/2"])
        # nothing under the skipped contrast in the table
        self.assertTrue(all(r["Contrast"] == "B - A" for r in self.onoff("complex")))
        with open(os.path.join(self.out("complex"), "methods.txt")) as fh:
            self.assertIn("Skipped (not a simple A-B contrast): C-(A+B)/2", fh.read())

    def test_loosening_coverage_does_not_change_the_list(self):
        # the list is counted before the filter, so --coverage-min does not move it
        key = lambda rows: sorted((r["Protein.Group"], r["Contrast"], r["Direction"]) for r in rows)
        self.assertEqual(key(self.onoff("cov0")), key(self.onoff("default")))

    def test_dpc_missing_values_are_called_estimates(self):          # 2.15.1, lane3 #15
        mp = self.prov("dpc")["missing_policy"]
        self.assertIn("model-based estimate", mp)
        self.assertNotIn("not imputed", mp)

    def test_onoff_sentence_only_when_the_table_was_written(self):   # 2.15.1, audit L2-1
        self.assertIn("OnOff_Proteins.csv", self.prov("default")["missing_policy"])
        self.assertNotIn("OnOff_Proteins.csv", self.prov("dpc")["missing_policy"])

    def test_dpc_is_untouched(self):
        out = self.out("dpc")
        self.assertFalse(os.path.exists(os.path.join(out, "OnOff_Proteins.csv")))
        self.assertIsNone(self.prov("dpc")["onoff"])
        self.assertNotEqual(self.refused.returncode, 0)
        self.assertIn("--onoff-min-n applies to --method maxlfq only", self.refused.stderr)

    def test_text_comes_from_the_record_only(self):
        """Methods, analysis note: built from de_provenance.json; absent on a record without `onoff`
        (dpc), and a failed table is said, not hidden."""
        prov = self.prov("default")
        s = mm.de_onoff_sentence(prov)
        self.assertIn("detection patterns, not fold changes", s)
        self.assertIn("B-A: " + str(prov["onoff"]["contrasts"]["B-A"]["n_onoff"]), s)
        self.assertIsNone(mm.de_onoff_sentence(self.prov("dpc")))
        note = mah.onoff_note(prov)
        self.assertEqual(note["kind"], "info")
        self.assertIn("OnOff_Proteins.csv", note["text"])
        self.assertIsNone(mah.onoff_note(self.prov("dpc")))
        bad = {"onoff": {"file": "OnOff_Proteins.csv", "status": "failed", "reason": "disk full"}}
        self.assertIn("could not be written", mm.de_onoff_sentence(bad))
        self.assertEqual(mah.onoff_note(bad)["kind"], "warning")
        skipped_only = {"onoff": {"status": "no_simple_contrasts", "file": "OnOff_Proteins.csv",
                                  "contrasts": {}, "skipped": {"X": "why"}}}
        self.assertIn("Not assessed", mah.onoff_note(skipped_only)["text"])

    def test_the_ai_brief_says_it(self):
        """analysis_prompt's block reads the record: counts per contrast, never a fold change."""
        import inspect
        src = inspect.getsource(ap)
        self.assertIn("On/off proteins", src)
        self.assertIn("Never give them a fold", src)


if __name__ == "__main__":
    unittest.main()
