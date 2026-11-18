#!/usr/bin/env python3
"""Three DE additions (2.15.x), each end to end through run_de.R:

1. `values_per_group` / `called_on_incomplete_data` / `incomplete_tier` on EVERY route's DE tables
   (dpc and maxlfq). "Measured" is the existing Detection_Matrix definition (detection_counts() in
   run_de.R reads det_n; incomplete_data.R only turns k and n into the flag). Counts in
   de_provenance.json (`incomplete_data`, `de_tables`), with the PXD063416 benchmark statement.
2. QuantUMS_tradeoff.csv (maxlfq only): the DE re-fitted at five cutoffs from the report read once.
   The row at the cutoffs actually used must equal the DE's own numbers (that is the test that the
   re-fit is the run's fit), the used row says where its cutoffs came from, and the skip paths
   (dpc, no QuantUMS columns, --no-quantums-tradeoff) are recorded, never silent.
3. On/off counted per SAMPLE where the conditions carry a Sample column.

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

from test_run_de_contaminants import r_has, rscript  # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")

# 3 groups (A, B, C) x 3 runs. Every protein has 3 precursors.
#   P00001-P00050  every run, 5 of them up 1.5 log2 in B
#   P00091         B runs 4 and 5 only (ONE sample when run4/run5 share a Sample)
#   P00092         B runs 4 and 6 only (TWO samples: on/off per sample as well as per run)
#   P00101-P00110  every run, +3 log2 in B; but PG.MaxLFQ.Quality 0.6 in run 1 (group A): at pgQ >= 0.75
#                  that one cell is censored, so A has 2/3 measured -> a call on incomplete data
#   P00111-P00115  every run, +3 log2 in B; Empirical.Quality 0.6 on every precursor of run 2 (A):
#                  censored only by an eQ cutoff
#   P00121-P00125  every run, PG.MaxLFQ.Quality 0.4 everywhere: gone at any pgQ >= 0.5 (proteins lost)
#   P00201-P00205  A: run 1 only (value 12); B and C every run at 18 -> if called, SINGLE-sample
SYNTH_R = r'''
set.seed(21)
runs <- sprintf("run%02d", 1:9); grp <- rep(c("A", "B", "C"), each = 3)
rows <- list()
add <- function(pg, prec, base, only_run = NULL, effB = 0, pgq = 0.95, eq = 0.95, fixed = NULL) {
  for (r in seq_along(runs)) {
    if (!is.null(only_run) && !(r %in% only_run)) next
    lv <- if (!is.null(fixed)) fixed[grp[r]] else base + (grp[r] == "B") * effB
    lv <- lv + rnorm(1, 0, 0.2)
    int <- 2^(lv + rnorm(1, 0, 0.1))
    rows[[length(rows) + 1]] <<- data.frame(Run = runs[r], Precursor.Id = prec,
      Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0(pg, "_X"), Genes = sub("P0*", "G", pg),
      Proteotypic = 1L, Precursor.Normalised = int, Precursor.Quantity = int,
      Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001,
      Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001,
      Empirical.Quality = eq, PG.MaxLFQ.Quality = pgq, stringsAsFactors = FALSE)
  }
}
for (i in c(1:50, 91, 92, 101:115, 121:125, 201:205)) {
  pg <- sprintf("P%05d", i); base <- runif(1, 16, 20)
  for (k in 1:3) {
    prec <- sprintf("PEP%dK%d2", i, k); b <- base + rnorm(1, 0, 0.3)
    if (i <= 50)       add(pg, prec, b, effB = if (i <= 5) 1.5 else 0)
    else if (i == 91)  add(pg, prec, b, only_run = 4:5)
    else if (i == 92)  add(pg, prec, b, only_run = c(4, 6))
    else if (i <= 110) add(pg, prec, b, effB = 3, pgq = 0)            # patched below for run 1
    else if (i <= 115) add(pg, prec, b, effB = 3)
    else if (i <= 125) add(pg, prec, b)
    else               add(pg, prec, 18, fixed = c(A = 12, B = 18, C = 18))
  }
}
d <- do.call(rbind, rows)
# P00201-205: A measured in run 1 only
d <- d[!(grepl("^P002", d$Protein.Group) & d$Run %in% c("run02", "run03")), ]
# quality: PG-level 0.6 in run 1 for P00101-110 (0.95 elsewhere); precursor-level 0.6 in run 2 for P00111-115
is110 <- d$Protein.Group %in% sprintf("P%05d", 101:110)
d$PG.MaxLFQ.Quality <- ifelse(is110 & d$Run == "run01", 0.6,
                              ifelse(d$Protein.Group %in% sprintf("P%05d", 121:125), 0.4, 0.95))
d$Empirical.Quality <- ifelse(d$Protein.Group %in% sprintf("P%05d", 111:115) & d$Run == "run02", 0.6, 0.95)
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum)
names(pgm)[3] <- "PG.MaxLFQ"
d <- merge(d, pgm, by = c("Protein.Group", "Run"))
arrow::write_parquet(d, file.path(OUT, "report.parquet"))
arrow::write_parquet(d[, setdiff(names(d), c("Empirical.Quality", "PG.MaxLFQ.Quality"))],
                     file.path(OUT, "report_noq.parquet"))
write.csv(data.frame(File.Name = runs, Group = grp), file.path(OUT, "conditions.csv"), row.names = FALSE)
smp <- paste0(tolower(grp), c(1, 1, 2))
write.csv(data.frame(File.Name = runs, Group = grp, Sample = smp), file.path(OUT, "conditions_tr.csv"),
          row.names = FALSE)
'''

LOW_PGQ = {f"P{i:05d}" for i in range(101, 111)}
SINGLE = {f"P{i:05d}" for i in range(201, 206)}
DEFAULT_TAG = "DEFAULT — not user-confirmed"


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def parse_vpg(s):
    out = {}
    for part in s.split("; "):
        g, kn = part.rsplit(" ", 1)
        k, n = kn.split("/")
        out[g] = (int(k), int(n))
    return out


@unittest.skipUnless(r_has("limma", "arrow", "dplyr", "tidyr", "jsonlite", "limpa", "statmod"),
                     "needs R with limma/limpa/statmod/arrow/dplyr/tidyr/jsonlite")
class IncompleteQuantumsOnOff(unittest.TestCase):
    RUNS = (
        ("ml", "report.parquet", "conditions.csv", "maxlfq", []),
        ("ml_cmplx", "report.parquet", "conditions.csv", "maxlfq", ["--contrasts", "B-A,C-(A+B)/2"]),
        ("ml_user", "report.parquet", "conditions.csv", "maxlfq",
         ["--pgq-cutoff", "0.75", "--eq-cutoff", "0.75"]),
        ("ml_pgq6", "report.parquet", "conditions.csv", "maxlfq", ["--pgq-cutoff", "0.6"]),
        ("ml_tr", "report.parquet", "conditions_tr.csv", "maxlfq", []),
        ("ml_noq", "report_noq.parquet", "conditions.csv", "maxlfq", []),
        ("ml_off", "report.parquet", "conditions.csv", "maxlfq", ["--no-quantums-tradeoff"]),
        ("dpc", "report.parquet", "conditions.csv", "dpc", []),
    )

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        cls.runs = {}
        for key, rep, meta, method, extra in cls.RUNS:
            out = os.path.join(cls.tmp, "out_" + key)
            p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(cls.tmp, rep),
                                "--metadata", os.path.join(cls.tmp, meta), "--method", method,
                                "--outdir", out] + extra, capture_output=True, text=True, cwd=cls.tmp)
            cls.runs[key] = (p, out, method)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, key):
        p, out, _ = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-3000:])
        return out

    def prov(self, key):
        with open(os.path.join(self.out(key), "de_provenance.json"), encoding="utf-8") as fh:
            return json.load(fh)

    def table(self, key, contrast):
        _, _, method = self.runs[key]
        return read_csv(os.path.join(self.out(key), f"DE_{method}_{contrast}.csv"))

    # ---- 1. the incomplete-data flag, every route ------------------------------------------
    def check_flags(self, key, contrast):
        rows = self.table(key, contrast)
        n_sig = 0
        for r in rows:
            vpg = parse_vpg(r["values_per_group"])
            sig = r["adj.P.Val"] not in ("NA", "") and float(r["adj.P.Val"]) < 0.05
            # values_per_group says exactly what Detected_<group> says (one definition of measured)
            for g, (k, n) in vpg.items():
                self.assertEqual(r[f"Detected_{g}"], f"{k}/{n}", (key, r["Protein.Group"]))
            if not sig:
                self.assertEqual((r["called_on_incomplete_data"], r["incomplete_tier"]), ("NA", "NA"),
                                 (key, r["Protein.Group"]))
                continue
            n_sig += 1
            incomplete = any(k < n for k, n in vpg.values())
            single = any(k == 1 and n > 1 for k, n in vpg.values())
            self.assertEqual(r["called_on_incomplete_data"], "TRUE" if incomplete else "FALSE", r)
            self.assertEqual(r["incomplete_tier"],
                             "single_sample" if single else "incomplete" if incomplete else "complete", r)
        return rows, n_sig

    def test_flags_on_maxlfq(self):
        rows, n_sig = self.check_flags("ml", "B.A")
        self.assertGreater(n_sig, 10)
        by = {r["Protein.Group"]: r for r in rows}
        # complete data, significant: FALSE
        self.assertEqual(by["P00001"]["called_on_incomplete_data"], "FALSE")
        self.assertEqual(by["P00001"]["values_per_group"], "B 3/3; A 3/3")
        # a protein seen in one run of A only, called: the single-sample tier
        hit = [by[p] for p in SINGLE if by[p]["called_on_incomplete_data"] == "TRUE"]
        self.assertTrue(hit, "fixture: at least one single-run protein should be called")
        self.assertTrue(all(r["incomplete_tier"] == "single_sample" for r in hit))
        self.assertTrue(all(r["values_per_group"].endswith("A 1/3") for r in hit))

    def test_flags_on_dpc(self):
        _, n_sig = self.check_flags("dpc", "B.A")
        self.assertGreater(n_sig, 0)
        prov = self.prov("dpc")
        self.assertIn("columns", prov["incomplete_data"])
        self.assertIn("0 = inferred", prov["incomplete_data"]["rule"])     # dpc's own meaning of "not measured"

    def test_a_complex_contrast_lists_every_group_it_uses(self):
        rows, _ = self.check_flags("ml_cmplx", "C..A.B..2")
        self.assertEqual(set(parse_vpg(rows[0]["values_per_group"])), {"A", "B", "C"})

    def test_flags_on_the_blocked_route_count_runs(self):
        rows, _ = self.check_flags("ml_tr", "B.A")
        self.assertEqual(parse_vpg(rows[0]["values_per_group"])["A"][1], 3)

    def test_provenance_counts_and_benchmark_statement(self):
        prov = self.prov("ml")
        inc = prov["incomplete_data"]
        self.assertIn("about 10x", inc["benchmark_text"])
        self.assertIn("PXD063416", inc["benchmark_text"])
        self.assertIn("FILTER_TUNING_PROTOTYPE", inc["citation"])
        self.assertIn("called_on_incomplete_data", inc["rule"])
        for cn, t in prov["de_tables"].items():
            rows = self.table("ml", cn.replace("-", "."))
            self.assertEqual(t["n_called_incomplete"],
                             sum(r["called_on_incomplete_data"] == "TRUE" for r in rows))
            self.assertEqual(t["n_called_single_sample"],
                             sum(r["incomplete_tier"] == "single_sample" for r in rows))
        with open(os.path.join(self.out("ml"), "methods.txt"), encoding="utf-8") as fh:
            self.assertIn("calls on incomplete data", fh.read().lower())

    # ---- 2. the QuantUMS trade-off ---------------------------------------------------------
    def tradeoff(self, key):
        return read_csv(os.path.join(self.out(key), "QuantUMS_tradeoff.csv"))

    def test_table_rows_and_the_default_tag(self):
        rows = self.tradeoff("ml")
        self.assertEqual({r["Cutoff"] for r in rows},
                         {"none", "pgQ 0.5", "pgQ 0.75", "pgQ 0.75 + eQ 0.75", "pgQ 0.85"})
        self.assertEqual(len(rows), 5 * 2)                      # x B-A, C-A
        used = [r for r in rows if r["used_in_this_run"] == "TRUE"]
        self.assertEqual({r["Cutoff"] for r in used}, {"none"})
        self.assertTrue(all(r["cutoff_source"] == DEFAULT_TAG for r in used))
        rec = [r for r in rows if "benchmark-recommended (PXD063416)" in r["note"]]
        self.assertEqual({r["Cutoff"] for r in rec}, {"pgQ 0.75 + eQ 0.75"})
        # "none" loses nothing relative to itself, a cutoff only ever loses depth
        for r in rows:
            if r["Cutoff"] == "none":
                self.assertEqual(r["proteins_lost_vs_none"], "0")
        for cn in ("B-A", "C-A"):
            t = {r["Cutoff"]: r for r in rows if r["Contrast"] == cn}
            self.assertLessEqual(int(t["pgQ 0.85"]["proteins_tested"]), int(t["pgQ 0.5"]["proteins_tested"]))
            self.assertGreaterEqual(int(t["pgQ 0.5"]["proteins_lost_vs_none"]), 5)

    def test_the_cutoff_censors_cells_so_calls_become_incomplete(self):
        t = {r["Cutoff"]: r for r in self.tradeoff("ml") if r["Contrast"] == "B-A"}
        # P00101-110: complete at none, one A cell censored at pgQ 0.75 -> 2/3 measured
        self.assertGreaterEqual(int(t["pgQ 0.75"]["called_on_incomplete_data"]),
                                int(t["none"]["called_on_incomplete_data"]) + len(LOW_PGQ) - 2)
        # the eQ cutoff does what the pgQ-only rows do not (P00111-115 lose a cell only to eQ)
        self.assertGreater(int(t["pgQ 0.75 + eQ 0.75"]["called_on_incomplete_data"]),
                           int(t["pgQ 0.75"]["called_on_incomplete_data"]))

    def used_row_matches_the_de(self, key, label):
        prov = self.prov(key)
        rows = [r for r in self.tradeoff(key) if r["used_in_this_run"] == "TRUE"]
        self.assertEqual({r["Cutoff"] for r in rows}, {label}, key)
        for r in rows:
            t = prov["de_tables"][r["Contrast"]]
            self.assertEqual(int(r["proteins_tested"]), t["n_tested"], (key, r))
            self.assertEqual(int(r["significant"]), t["n_significant"], (key, r))
            self.assertEqual(int(r["called_on_incomplete_data"]), t["n_called_incomplete"], (key, r))
            self.assertEqual(int(r["called_single_sample"]), t["n_called_single_sample"], (key, r))

    def test_the_row_at_the_cutoffs_used_is_the_des_own_fit(self):
        self.used_row_matches_the_de("ml", "none")
        self.used_row_matches_the_de("ml_user", "pgQ 0.75 + eQ 0.75")
        self.used_row_matches_the_de("ml_tr", "none")             # the blocked (duplicateCorrelation) fit
        self.used_row_matches_the_de("ml_cmplx", "none")

    def test_user_chosen_cutoffs_are_not_tagged_default(self):
        used = [r for r in self.tradeoff("ml_user") if r["used_in_this_run"] == "TRUE"]
        self.assertTrue(all(r["cutoff_source"] == "pgQ: user-chosen; eQ: user-chosen" for r in used))
        self.assertEqual(self.prov("ml_user")["quantums_tradeoff"]["cutoff_source"],
                         "pgQ: user-chosen; eQ: user-chosen")

    def test_a_cutoff_outside_the_five_is_added_as_run_and_partly_default(self):
        rows = self.tradeoff("ml_pgq6")
        used = [r for r in rows if r["used_in_this_run"] == "TRUE"]
        self.assertEqual({r["Cutoff"] for r in used}, {"as run (pgQ 0.6, eQ 0)"})
        self.assertTrue(all(r["cutoff_source"] == f"pgQ: user-chosen; eQ: {DEFAULT_TAG}" for r in used))
        self.assertEqual(len({r["Cutoff"] for r in rows}), 6)
        self.used_row_matches_the_de("ml_pgq6", "as run (pgQ 0.6, eQ 0)")

    def test_provenance_record_and_defaults_unchanged(self):
        q = self.prov("ml")["quantums_tradeoff"]
        self.assertEqual(q["status"], "written")
        self.assertEqual(q["used"], {"pgQ": 0, "eQ": 0})                 # the skill's default is still off
        self.assertEqual(q["cutoff_source"], DEFAULT_TAG)
        self.assertEqual(q["recommended"]["tag"], "benchmark-recommended (PXD063416)")
        self.assertIn("FILTER_TUNING_PROTOTYPE", q["recommended"]["citation"])
        self.assertIn("unchanged", q["defaults_unchanged"])
        self.assertEqual(len(q["table"]), 10)
        self.assertEqual(self.prov("ml")["pgq_cutoff"], 0)

    def test_skip_paths_are_recorded_not_silent(self):
        for key, want in (("dpc", "imputes"), ("ml_noq", "no QuantUMS quality columns"),
                          ("ml_off", "--no-quantums-tradeoff")):
            with self.subTest(run=key):
                q = self.prov(key)["quantums_tradeoff"]
                self.assertEqual(q["status"], "skipped")
                self.assertIn(want, q["reason"])
                self.assertFalse(os.path.exists(os.path.join(self.out(key), "QuantUMS_tradeoff.csv")))

    # ---- 3. on/off per sample --------------------------------------------------------------
    def onoff(self, key):
        return read_csv(os.path.join(self.out(key), "OnOff_Proteins.csv"))

    def test_onoff_per_sample_collapses_injections(self):
        rec = self.prov("ml_tr")["onoff"]
        self.assertTrue(rec["collapsed_runs_to_samples"])
        self.assertEqual((rec["n_runs"], rec["n_samples"], rec["n_min_unit"]), (9, 6, "sample"))
        self.assertIn("conditions' Sample column: 9 runs collapsed to 6 samples", rec["unit_source"])
        self.assertIn("at least 2 samples", rec["rule"])
        self.assertIn("any of its runs", rec["rule"])
        # P00091 is in B runs 4+5 = ONE sample: per run it is B_only, per sample it is not
        per_run = {r["Protein.Group"] for r in self.onoff("ml")}
        per_smp = {r["Protein.Group"] for r in self.onoff("ml_tr")}
        self.assertIn("P00091", per_run)
        self.assertNotIn("P00091", per_smp)
        self.assertIn("P00092", per_smp)
        # counts are in samples: B has 2 samples, A 2, C 2
        rows = self.onoff("ml_tr")
        self.assertTrue(rows)
        self.assertEqual({x["total_g1"] for x in rows} | {x["total_g2"] for x in rows}, {"2"})

    def test_onoff_without_a_sample_column_stays_per_run_and_says_so(self):
        rec = self.prov("ml")["onoff"]
        self.assertFalse(rec["collapsed_runs_to_samples"])
        self.assertEqual(rec["n_min_unit"], "run")
        self.assertIn("no Sample column", rec["unit_source"])
        self.assertIn("at least 2 runs", rec["rule"])


if __name__ == "__main__":
    unittest.main()
