#!/usr/bin/env python3
"""run_de.R fits a continuous covariate as one numeric column when asked (2.12).

A staff report (2026-10-08): build_design() wrapped every Batch/Covariate1/Covariate2 in factor(),
so a continuous covariate -- a centred log2 TurboID intensity with 11 distinct values -- became an
11-level factor and "Design matrix is not full rank". Now `--numeric-covariates Covariate1[,...]`
fits it as one numeric column: one slope shared by all groups, so each contrast between groups
compares them at the same covariate value. Never decided from the values (a batch numbered 1/2/3
stays a factor; Batch always is). A value that is not a finite number is refused, never filled
in. Group explaining >= 80% of its variance is warned about. de_provenance.json `covariates`,
methods.txt, the Methods paragraph (make_methods.de_covariate_sentence), the analysis brief, the
report's callout and reproducibility_log.R all say how it was fitted.

End to end on a synthetic DIA-NN report (each protein its own slope on the covariate, 12 proteins
up 1 log2 in B); skipped without R + limpa/limma/arrow/dplyr/tidyr/jsonlite.
NOTE: CI runs Python only and skips every R-invoking test -- run them locally.
"""
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_analysis_html as mah  # noqa: E402
import make_methods as mm  # noqa: E402
from test_run_de_block import r_has, rscript  # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")
NEEDS = ("limpa", "limma", "statmod", "arrow", "dplyr", "tidyr", "jsonlite")

# 12 runs, A x6 and B x6; Covariate1 a continuous measurement, 12 distinct values, barely related
# to Group. Every protein has its own slope on it (mean ~0, so normalisation leaves it in place);
# P00001-P00012 are up 1 log2 in B. Low-intensity precursors go missing more often (limpa's
# detection model needs some missing values).
SYNTH_R = r'''
set.seed(11)
runs <- sprintf("run%02d", 1:12); grp <- rep(c("A", "B"), each = 6)
cov <- c(-1.1, -0.6, -0.2, 0.3, 0.7, 1.2, -0.9, -0.4, 0.1, 0.5, 0.9, 1.4)
out <- list()
for (i in 1:120) {
  pg <- sprintf("P%05d", i); base <- runif(1, 12, 20); slope <- runif(1, -1.5, 1.5)
  run_lv <- base + slope * cov + (grp == "B") * (if (i <= 12) 1 else 0) + rnorm(12, 0, 0.15)
  for (j in 1:3) {
    lv <- run_lv + rnorm(1, 0, 1) + rnorm(12, 0, 0.1)
    seen <- runif(12) >= plogis(-(lv - 13) * 2)
    if (!any(seen)) next
    out[[length(out) + 1]] <- data.frame(Run = runs[seen], Precursor.Id = sprintf("PEP%dK%d2", i, j),
      Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0("G", i, "_X"),
      Genes = paste0("G", i), Proteotypic = 1L, Precursor.Normalised = 2^lv[seen],
      Precursor.Quantity = 2^lv[seen], Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001,
      PG.Q.Value = 0.001, Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001,
      stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, out)
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum)
names(pgm)[3] <- "PG.MaxLFQ"
arrow::write_parquet(merge(d, pgm, by = c("Protein.Group", "Run")), file.path(OUT, "report.parquet"))
w <- function(m, f) write.csv(m, file.path(OUT, f), row.names = FALSE)
meta <- data.frame(File.Name = runs, Group = grp, Covariate1 = cov)
w(meta, "conditions.csv")
w(meta[, c("File.Name", "Group")], "conditions_nocov.csv")
w(transform(meta, Batch = rep(c("b1", "b2", "b3"), 4))[, c("File.Name", "Group", "Batch", "Covariate1")],
  "conditions_batch.csv")
# Covariate1 nearly the grouping: A around 0, B around 1 (R^2 with Group ~0.95)
w(transform(meta, Covariate1 = c(0.05, 0.1, -0.08, 0.2, 0.0, 0.12, 0.95, 1.1, 0.9, 1.2, 1.0, 0.85)),
  "conditions_conf.csv")
# paired: run01 and run07 from patient P1 ...; Covariate1 = the patient's age (constant per patient)
w(transform(meta, Patient = paste0("P", c(1:6, 1:6)), Covariate1 = c(30, 41, 52, 35, 60, 47)[c(1:6, 1:6)]),
  "conditions_paired.csv")
# technical replicates: run01/run02 are one sample, run07/run08 another
w(transform(meta, Sample = c("S1", "S1", sprintf("S%d", 3:6), "S7", "S7", sprintf("S%d", 9:12))),
  "conditions_tech.csv")
'''


def run_de(tmp, meta, *args):
    out = os.path.join(tmp, "out_" + os.path.splitext(meta)[0] + "_" + "_".join(
        re.sub(r"[^A-Za-z0-9]+", "", a) for a in args))
    p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(tmp, "report.parquet"),
                        "--metadata", os.path.join(tmp, meta), "--outdir", out, *args],
                       capture_output=True, text=True, cwd=tmp)
    return p, out


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@unittest.skipUnless(r_has(*NEEDS), "needs R with " + "/".join(NEEDS))
class NumericCovariateEndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        num = ("--numeric-covariates", "Covariate1")
        cls.runs = {
            "factor": run_de(cls.tmp, "conditions.csv", "--method", "maxlfq"),
            "nocov": run_de(cls.tmp, "conditions_nocov.csv", "--method", "maxlfq"),
            "ml": run_de(cls.tmp, "conditions.csv", "--method", "maxlfq", *num),
            "dpc": run_de(cls.tmp, "conditions.csv", "--method", "dpc", *num),
            "batch": run_de(cls.tmp, "conditions_batch.csv", "--method", "maxlfq", *num),
            "conf": run_de(cls.tmp, "conditions_conf.csv", "--method", "maxlfq", *num),
            "tech_ml": run_de(cls.tmp, "conditions_tech.csv", "--method", "maxlfq", *num),
            "tech_dpc": run_de(cls.tmp, "conditions_tech.csv", "--method", "dpc", *num),
            "paired": run_de(cls.tmp, "conditions_paired.csv", "--method", "maxlfq", *num,
                             "--block", "Patient"),
        }

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def ok(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-2500:])
        return out

    def prov(self, key):
        with open(os.path.join(self.ok(key), "de_provenance.json")) as fh:
            return json.load(fh)

    def test_without_the_flag_it_stays_a_factor_and_the_error_names_the_flag(self):
        p, _ = self.runs["factor"]
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("not full rank", p.stderr)
        self.assertIn("Covariate1 holds only numbers (12 distinct values) and is fitted as a factor",
                      p.stderr)
        self.assertIn("--numeric-covariates Covariate1", p.stderr)

    def test_it_is_one_numeric_column_and_the_record_says_so(self):
        for key, method in (("ml", "maxlfq"), ("dpc", "dpc")):
            prov = self.prov(key)
            self.assertEqual(prov["design"], "~ 0 + groups + Covariate1", key)
            (t,) = prov["covariates"]["terms"]
            self.assertEqual((t["column"], t["type"], t["source"], t["n_distinct"]),
                             ("Covariate1", "numeric", "--numeric-covariates", 12), key)
            self.assertEqual((t["min"], t["max"]), (-1.1, 1.4), key)
            self.assertLess(t["group_r2"], 0.1, key)
            self.assertFalse(t["confounded"], key)
            self.assertEqual(list(t["se_inflation"]), ["B-A"], key)
            self.assertEqual(prov["covariates"]["numeric_requested"], ["Covariate1"])
            self.assertEqual(prov["covariates"]["warnings"], [])
            with open(os.path.join(self.ok(key), "methods.txt")) as fh:
                methods = fh.read()
            self.assertIn("Covariates    : Covariate1 numeric (--numeric-covariates; -1.1 to 1.4): "
                          "one slope shared by all groups", methods, key)

    def test_the_contrast_is_between_groups_with_the_covariate_held_fixed(self):
        """The design is A, B, Covariate1; contrast B-A puts 0 on the covariate's column, so its
        logFC is the B-A difference adjusted for Covariate1 -- what lm(y ~ 0 + group + cov) gives
        for every protein (maxlfq: unweighted lmFit, so exactly)."""
        for key in ("ml", "dpc"):
            out = self.ok(key)
            r = rscript(f'''
              x <- readRDS("{out}/set_test_inputs.rds")
              cat(colnames(x$design), "|", x$cmat[, "B-A"], "|", x$design[, "Covariate1"])''')
            cols, weights, values = (s.split() for s in r.split("|"))
            self.assertEqual(cols, ["A", "B", "Covariate1"], key)
            self.assertEqual([float(v) for v in weights], [-1, 1, 0], key)
            self.assertEqual([float(v) for v in values],
                             [-1.1, -0.6, -0.2, 0.3, 0.7, 1.2, -0.9, -0.4, 0.1, 0.5, 0.9, 1.4], key)
        out = self.ok("ml")
        diff = rscript(f'''
          em <- read.csv("{out}/Expression_Matrix.csv", check.names = FALSE)
          meta <- read.csv("{self.tmp}/conditions.csv")
          E <- as.matrix(em[, meta$File.Name]); rownames(E) <- em$Protein.Group
          g <- factor(meta$Group); cv <- meta$Covariate1
          ols <- apply(E, 1, function(y) {{ b <- coef(lm(y ~ 0 + g + cv)); b[["gB"]] - b[["gA"]] }})
          tt <- read.csv("{out}/DE_maxlfq_B.A.csv")
          cat(max(abs(tt$logFC - ols[tt$Protein.Group])))''')
        self.assertLess(float(diff), 1e-8)

    def test_the_covariate_takes_its_variance_out_of_the_comparison(self):
        """Each protein has its own slope on Covariate1: fitted as a covariate, most of the twelve
        proteins that really change are found, and more than with no covariate at all; Covariate1
        barely relates to Group, so the comparison's standard error is hardly widened."""
        def found(key):
            rows = read_csv(os.path.join(self.ok(key), "DE_maxlfq_B.A.csv"))
            sig = {r["Protein.Group"] for r in rows if r["adj.P.Val"] not in ("", "NA")
                   and float(r["adj.P.Val"]) < 0.05}
            return len(sig & {f"P{i:05d}" for i in range(1, 13)})
        self.assertGreaterEqual(found("ml"), 10)
        self.assertGreater(found("ml"), found("nocov"))
        self.assertLess(self.prov("ml")["covariates"]["terms"][0]["se_inflation"]["B-A"], 1.05)
        self.assertEqual(self.prov("nocov")["covariates"]["terms"], [])

    def test_the_reproducibility_script_fits_the_same_numeric_covariate(self):
        for key, method in (("ml", "maxlfq"), ("dpc", "dpc")):
            out = self.ok(key)
            with open(os.path.join(out, "reproducibility_log.R")) as fh:
                src = fh.read()
            self.assertIn("covariates: Covariate1 (numeric)", src)
            self.assertIn("  'run01' = -1.1,", src)
            self.assertIn("Covariate1 <- as.numeric(metadata$Covariate1)", src)
            self.assertNotIn("Covariate1 <- factor(", src)
            rerun = os.path.join(self.tmp, "rerun_" + key)
            os.makedirs(rerun)
            p = subprocess.run(["Rscript", os.path.join(out, "reproducibility_log.R")],
                               capture_output=True, text=True, cwd=rerun)
            self.assertEqual(p.returncode, 0, p.stderr[-2000:])
            a = {r["Protein.Group"]: r for r in read_csv(os.path.join(out, f"DE_{method}_B.A.csv"))}
            b = {r["Protein.Group"]: r for r in read_csv(
                os.path.join(rerun, "de_results_rerun", f"DE_{method}_B.A.csv"))}
            self.assertEqual(set(a), set(b), key)
            for pg in a:
                for col in ("logFC", "P.Value"):
                    self.assertAlmostEqual(float(a[pg][col]), float(b[pg][col]), places=10,
                                           msg=(key, pg, col))

    def test_a_factor_beside_it_stays_a_factor(self):
        prov = self.prov("batch")
        self.assertEqual(prov["design"], "~ 0 + groups + Batch + Covariate1")
        terms = {t["column"]: t for t in prov["covariates"]["terms"]}
        self.assertEqual((terms["Batch"]["type"], terms["Batch"]["n_levels"]), ("factor", 3))
        self.assertEqual(terms["Covariate1"]["type"], "numeric")
        sent = mm.de_covariate_sentence(prov)
        self.assertIn("Covariate1 was fitted as a continuous covariate (values -1.1 to 1.4", sent)
        self.assertIn("Batch was a categorical factor (3 levels).", sent)

    def test_a_covariate_that_nearly_matches_the_groups_is_warned_about(self):
        p, out = self.runs["conf"]
        prov = self.prov("conf")
        t = prov["covariates"]["terms"][0]
        self.assertTrue(t["confounded"])
        self.assertGreater(t["group_r2"], 0.9)
        self.assertGreater(t["se_inflation"]["B-A"], 3)
        self.assertEqual(t["confounded_contrasts"], ["B-A"])
        pc = t["per_contrast"]["B-A"]
        self.assertEqual((pc["groups"], pc["overlap"], pc["confounded"]), (["B", "A"], False, True))
        self.assertIn("Covariate1 is strongly confounded with the groups in B-A (", p.stderr)
        self.assertIn("the design alone widens its standard error", p.stderr)
        self.assertEqual(len(prov["covariates"]["warnings"]), 1)
        self.assertIn("removes part of the effect being tested", prov["covariates"]["warnings"][0])
        with open(os.path.join(out, "methods.txt")) as fh:
            self.assertIn("CAUTION: Covariate1 is strongly confounded with the groups in B-A", fh.read())
        self.assertIn("In B-A, Covariate1 nearly matched the groups compared (", mm.de_covariate_sentence(prov))

    def test_the_audit_surfaces_it(self):
        """audit_results.py's numeric_covariate check quotes run_de.R's caution (WARN), and passes
        the covariate that does not match the groups."""
        def audit(key):
            out = self.ok(key)
            md = os.path.join(self.tmp, f"AUDIT_{key}.md")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                                "--out", md, "--de-dir", out], capture_output=True, text=True,
                               cwd=self.tmp, timeout=120)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(md[:-3] + ".json") as fh:
                return [f for f in json.load(fh)["findings"] if f["check"] == "numeric_covariate"]
        (f,) = audit("conf")
        self.assertEqual(f["status"], "WARN")
        self.assertTrue(f["message"].startswith(self.prov("conf")["covariates"]["warnings"][0]))
        (f,) = audit("ml")
        self.assertEqual(f["status"], "PASS")
        self.assertEqual(audit("nocov"), [])

    def test_technical_replicates_are_blocked_with_it_on_both_methods(self):
        for key in ("tech_ml", "tech_dpc"):
            prov = self.prov(key)
            self.assertTrue(prov["block"]["applied"], key)
            self.assertEqual((prov["block"]["column"], prov["block"]["effect"]),
                             ("Sample", "random"), key)
            self.assertEqual(prov["design"], "~ 0 + groups + Covariate1", key)
            self.assertEqual(prov["covariates"]["terms"][0]["type"], "numeric", key)
            self.assertEqual(set(prov["block"]["contrast_model"].values()), {"blocked"}, key)

    def test_a_fixed_block_absorbs_a_covariate_constant_within_it(self):
        """Each patient's age is constant within the patient: the fixed Patient effect absorbs it
        (blocking.R block_choose_effect), as it would a factor, and the record says so."""
        prov = self.prov("paired")
        self.assertEqual((prov["block"]["effect"], prov["design"]), ("fixed", "~ 0 + groups + Patient"))
        self.assertEqual(prov["block"]["absorbed_covariates"], ["Covariate1"])
        self.assertEqual(prov["covariates"]["terms"], [])
        self.assertEqual(prov["covariates"]["absorbed_by_block"], ["Covariate1"])
        self.assertIsNone(mm.de_covariate_sentence(prov))


@unittest.skipUnless(r_has("limma", "jsonlite"), "needs R with limma/jsonlite")
class NumericCovariateRefusals(unittest.TestCase):
    """Refused before quantification: a value that is not a number is never filled in."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        meta = [("run01", "A", "-1.1"), ("run02", "A", "0.3"), ("run03", "A", "0.7"),
                ("run04", "B", "0.2"), ("run05", "B", "0.9"), ("run06", "B", "1.4")]

        def write(name, rows, header=("File.Name", "Group", "Covariate1")):
            with open(os.path.join(cls.tmp, name), "w", newline="") as fh:
                csv.writer(fh).writerows([header, *rows])
        write("ok.csv", meta)
        write("blank.csv", [m if m[0] != "run03" else (m[0], m[1], "") for m in meta])
        write("text.csv", [m if m[0] != "run05" else (m[0], m[1], "high") for m in meta])
        write("inf.csv", [m if m[0] != "run05" else (m[0], m[1], "Inf") for m in meta])
        write("const.csv", [(m[0], m[1], "0.5") for m in meta])
        write("batch.csv", [(*m, "1") for m in meta], ("File.Name", "Group", "Covariate1", "Batch"))
        open(os.path.join(cls.tmp, "report.parquet"), "w").close()   # never read: refused first

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def refused(self, meta, *args):
        p, _ = run_de(self.tmp, meta, "--method", "maxlfq", *args)
        self.assertNotEqual(p.returncode, 0, p.stderr)
        return p.stderr

    def test_each_bad_value_or_column_is_refused_with_why(self):
        num = ("--numeric-covariates", "Covariate1")
        self.assertIn("1 sample(s) have no number in Covariate1: run03 (blank)",
                      self.refused("blank.csv", *num))
        self.assertIn("run05 ('high')", self.refused("text.csv", *num))
        self.assertIn("run05 ('Inf')", self.refused("inf.csv", *num))
        self.assertIn("never a made-up value", self.refused("text.csv", *num))
        self.assertIn("every sample has the same value", self.refused("const.csv", *num))
        self.assertIn("Batch is a categorical label and is always fitted as a factor",
                      self.refused("batch.csv", "--numeric-covariates", "Batch"))
        self.assertIn("has no Covariate2 column",
                      self.refused("ok.csv", "--numeric-covariates", "Covariate1,Covariate2"))
        self.assertIn("only Covariate1 or Covariate2 can be numeric covariates",
                      self.refused("ok.csv", "--numeric-covariates", "Age"))
        self.assertIn("--numeric-covariates needs column names",
                      self.refused("ok.csv", "--numeric-covariates", "--keep-contaminants"))


@unittest.skipUnless(shutil.which("Rscript"), "needs Rscript")
class ReproEmitter(unittest.TestCase):
    """Structural: the emitter writes a numeric covariate as numbers, exactly, and the script parses."""

    def test_numbers_round_trip_and_the_script_parses(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "r.R")
            src = rscript(f'''
              source("{SCRIPTS}/diann_q_columns.R"); source("{SCRIPTS}/repro_script.R")
              meta <- data.frame(File.Name = c("r1","r2","r3","r4"), Group = c("A","A","B","B"),
                                 Batch = c("1","2","1","2"), Covariate1 = c(0.1 + 0.2, 1/3, -2, 5))
              invisible(write_repro_script(path = "{out}", method = "maxlfq", input = "report.parquet",
                format = "parquet", q_cutoff = 0.01, q_columns = "Q.Value", q_cutoffs = 0.01,
                eq_cutoff = 0, pgq_cutoff = 0, cov_min_frac = 0.5, meta = meta,
                covariates = c("Batch", "Covariate1"), numeric_covariates = "Covariate1",
                formula_parts = c("groups", "Batch", "Covariate1"), forms = "B-A",
                adjp_thr = 0.05, logfc_ref = 1))
              e <- new.env(); invisible(parse("{out}"))
              L <- readLines("{out}"); i <- grep("^covariate1_map", L)
              eval(parse(text = paste(L[i:(i + 5)], collapse = "\\n")), e)
              cat(identical(unname(e$covariate1_map), meta$Covariate1), "\\n")
              cat(L, sep = "\\n")''')
        self.assertTrue(src.startswith("TRUE"), src[:200])
        self.assertIn("covariates: Batch, Covariate1 (numeric)", src)
        self.assertIn("  'r1' = 0.30000000000000004,", src)
        self.assertIn("  'r1' = '1',", src)                    # Batch: still a quoted label
        self.assertIn("Batch <- factor(metadata$Batch)", src)
        self.assertIn("Covariate1 <- as.numeric(metadata$Covariate1)", src)


class ReproduceShReplaysIt(unittest.TestCase):
    def test_reproduce_sh_passes_the_numeric_covariates(self):
        with tempfile.TemporaryDirectory() as d:
            de = os.path.join(d, "tables")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump({"method": "maxlfq", "quantities": "normalised",
                           "covariates": {"terms": [], "numeric_requested": ["Covariate1"]}}, fh)
            out = os.path.join(d, "repro")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir",
                                out, "--de-dir", de, "--de-method", "maxlfq",
                                "--timestamp", "2026-10-15T00:00:00Z"],
                               capture_output=True, text=True, timeout=300)
            self.assertEqual(p.returncode, 0, p.stderr[-1500:])
            with open(os.path.join(out, "reproduce.sh")) as fh:
                self.assertIn("--quantities normalised \\\n  --numeric-covariates Covariate1\n",
                              fh.read())


class TheWordsReadTheRecord(unittest.TestCase):
    """make_methods.de_covariate_sentence is the one description; the brief and the page quote it."""
    NUM = {"column": "Covariate1", "type": "numeric", "source": "--numeric-covariates",
           "n_distinct": 11, "min": -1.25, "max": 2.5, "group_r2": 0.21, "confounded": False,
           "se_inflation": {"B-A": 1.12}}

    def prov(self, *terms, warnings=()):
        return {"design": "~ 0 + groups + " + " + ".join(t["column"] for t in terms),
                "contrasts": ["B-A"],
                "covariates": {"terms": list(terms), "numeric_requested": ["Covariate1"],
                               "warnings": list(warnings)}}

    def test_nothing_is_said_without_a_numeric_covariate(self):
        self.assertIsNone(mm.de_covariate_sentence({}))                     # before 2.12
        self.assertIsNone(mm.de_covariate_sentence(self.prov(
            {"column": "Batch", "type": "factor", "n_levels": 3})))
        self.assertIsNone(mah.covariate_note({"covariates": {"terms": []}}))

    def test_a_numeric_covariate_is_described_in_the_methods(self):
        prov = self.prov(self.NUM)
        sent = mm.de_covariate_sentence(prov)
        self.assertEqual(sent, "Covariate1 was fitted as a continuous covariate (values -1.25 to "
                               "2.5, as given in the conditions file), with one slope shared by "
                               "all groups, so each group comparison is adjusted for it and "
                               "compares the groups at the same covariate value.")
        self.assertIn(sent, mm.de_paragraph(prov))
        note = mah.covariate_note(prov)
        self.assertEqual((note["kind"], note["text"]), ("info", sent))

    def test_a_confounded_one_is_a_warning_everywhere(self):
        t = dict(self.NUM, group_r2=0.93, confounded=True, se_inflation={"B-A": 3.8})
        prov = self.prov(t, warnings=["Covariate1 is strongly confounded with Group: ..."])
        sent = mm.de_covariate_sentence(prov)
        # a record with no per-contrast detail: the overall share, the design-only factor said so
        self.assertTrue(sent.endswith("Group explained 93% of the variance of Covariate1, so the "
                                      "adjusted group comparisons rest largely on its fitted slope "
                                      "(the design alone widens their standard errors up to 3.8x)."),
                        sent)
        self.assertEqual(mah.covariate_note(prov)["kind"], "warning")
        with tempfile.TemporaryDirectory() as tmp:
            de = os.path.join(tmp, "tables")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump(prov, fh)
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                                "--de-dir", de, "--out", os.path.join(tmp, "P.md")],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(os.path.join(tmp, "P.md"), encoding="utf-8") as fh:
                brief = fh.read()
        self.assertIn(f"- Continuous covariate: {sent} Describe every group difference as "
                      f"adjusted for it.", brief)
        self.assertIn("- **Covariate caveat (say this in Data Quality Notes):** Covariate1 is "
                      "strongly confounded with Group", brief)


# Four groups x 4 runs; Covariate1 separates B from A completely but spreads C and D over the whole
# range: overall R^2 0.595 -- below the caution -- yet 0.988 within A and B (the 2.12 review's
# example, rev-stats/r2.R). Every protein has its own slope; P00001-P00010 are up 1 log2 in B.
SYNTH4_R = r'''
set.seed(3)
grp <- rep(c("A", "B", "C", "D"), each = 4); runs <- sprintf("run%02d", seq_along(grp))
cov <- c(0.0, 0.1, 0.2, 0.3,  2.0, 2.1, 2.2, 2.3,  0, 0.8, 1.5, 2.3,  0.1, 0.9, 1.6, 2.2)
out <- list()
for (i in 1:100) {
  pg <- sprintf("P%05d", i); base <- runif(1, 14, 20); slope <- runif(1, -1, 1)
  run_lv <- base + slope * cov + (grp == "B") * (if (i <= 10) 1 else 0) + rnorm(16, 0, 0.15)
  for (j in 1:3) {
    lv <- run_lv + rnorm(1, 0, 0.8) + rnorm(16, 0, 0.1)
    out[[length(out) + 1]] <- data.frame(Run = runs, Precursor.Id = sprintf("PEP%dK%d2", i, j),
      Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0("G", i, "_X"),
      Genes = paste0("G", i), Proteotypic = 1L, Precursor.Normalised = 2^lv,
      Precursor.Quantity = 2^lv, Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001,
      PG.Q.Value = 0.001, Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001,
      stringsAsFactors = FALSE)
  }
}
d <- do.call(rbind, out)
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum)
names(pgm)[3] <- "PG.MaxLFQ"
arrow::write_parquet(merge(d, pgm, by = c("Protein.Group", "Run")), file.path(OUT, "report.parquet"))
write.csv(data.frame(File.Name = runs, Group = grp, Covariate1 = cov), file.path(OUT, "conditions.csv"),
          row.names = FALSE)
'''


@unittest.skipUnless(r_has("limma", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limma/arrow/dplyr/tidyr/jsonlite")
class ConfoundingIsPerContrast(unittest.TestCase):
    """The 2.12 review (M2): the caution used the overall R^2 only, so with 3+ groups a pair the
    covariate separates went unflagged. Now each contrast is judged on its own groups' samples
    (R^2 within them, overlapping ranges, the design's standard-error inflation), and the
    caution, methods.txt, the Methods, the page callout, the brief and the audit name it."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH4_R)
        cls.p, cls.out = run_de(cls.tmp, "conditions.csv", "--method", "maxlfq", "--numeric-covariates",
                                "Covariate1", "--contrasts", "B-A,C-A,D-A,C-D")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def prov(self):
        self.assertEqual(self.p.returncode, 0, self.p.stderr[-2500:])
        with open(os.path.join(self.out, "de_provenance.json")) as fh:
            return json.load(fh)

    def test_the_pair_it_separates_is_flagged_though_the_overall_share_is_modest(self):
        t = self.prov()["covariates"]["terms"][0]
        self.assertAlmostEqual(t["group_r2"], 0.595, places=2)          # under 0.8 overall
        self.assertTrue(t["confounded"])
        self.assertEqual(t["confounded_contrasts"], ["B-A"])
        pc = t["per_contrast"]
        self.assertAlmostEqual(pc["B-A"]["r2"], 0.988, places=2)
        self.assertFalse(pc["B-A"]["overlap"])
        self.assertAlmostEqual(pc["B-A"]["se_inflation"], 1.57, places=2)
        self.assertEqual(len(pc["B-A"]["why"]), 2)                     # R^2 and no overlap; 1.57x < 2.24x
        for cn in ("C-A", "D-A", "C-D"):
            self.assertFalse(pc[cn]["confounded"], cn)
            self.assertTrue(pc[cn]["overlap"], cn)
        self.assertEqual(pc["C-D"]["groups"], ["C", "D"])

    def test_every_reader_names_the_contrast(self):
        prov = self.prov()
        (w,) = prov["covariates"]["warnings"]
        self.assertTrue(w.startswith("Covariate1 is strongly confounded with the groups in B-A (the "
                                     "compared groups explain 99% of its variance among their samples"),
                        w)
        self.assertIn("their values do not overlap (B 2 to 2.3; A 0 to 0.3)", w)
        self.assertIn("That comparison is adjusted", w)
        self.assertNotIn("C-A", w)
        self.assertIn(w, self.p.stderr.replace("\n", " "))
        with open(os.path.join(self.out, "methods.txt")) as fh:
            self.assertIn("CAUTION: " + w, " ".join(fh.read().split()))
        sent = mm.de_covariate_sentence(prov)
        self.assertIn("In B-A, Covariate1 nearly matched the groups compared (the compared groups "
                      "explain 99% of its variance", sent)
        self.assertNotIn("C-A", sent)
        note = mah.covariate_note(prov)
        self.assertEqual((note["kind"], note["text"]), ("warning", sent))
        md = os.path.join(self.tmp, "AUDIT.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"), "--out", md,
                            "--de-dir", self.out], capture_output=True, text=True, cwd=self.tmp,
                           timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(md[:-3] + ".json") as fh:
            (f,) = [x for x in json.load(fh)["findings"] if x["check"] == "numeric_covariate"]
        self.assertEqual(f["status"], "WARN")
        self.assertTrue(f["message"].startswith(w), f["message"])
        self.assertEqual(f["detail"]["confounded_contrasts"], ["B-A"])


if __name__ == "__main__":
    unittest.main()
