#!/usr/bin/env python3
"""run_de.R --block: a random blocking factor (the mouse several IPs came from).

PROT_0756 (Dickson lab, Silva08172026): 30 IPs from 6 mouse brains, one JPH3/JPH4/Kv2.1/
RyR/IgG IP per mouse, mice 1-3 Old and 4-6 Young. Every bait-vs-IgG contrast compares IPs
from the same mice, but run_de.R fitted all 30 as independent samples.

Guards (end to end on synthetic DIA-NN reports; skip without R + limpa/limma/statmod/...):
  * a paired design where subjects differ a lot: blocking finds the true changes the
    unblocked fit misses, on dpc and maxlfq;
  * de_provenance.json / methods.txt / make_methods.py record the block, its consensus
    correlation and how it was fitted; an unblocked run says it modelled independence;
  * the nested design (subject within age) fits, and each contrast is labelled within- or
    between-subject;
  * the block cannot also be a fixed covariate, the Group itself, encoded in the design,
    or hold a single sample -- each stops BEFORE quantification with a clear message; a
    subject column fitted as a fixed covariate inside the groups points at --block;
  * reproducibility_log.R refits the blocked model and reproduces the DE tables;
  * the record's warnings: correlation <= 0, few proteins, a moving estimate, 2 blocks.
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

import make_methods as mm   # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")
BLOCK_R = os.path.join(SCRIPTS, "blocking.R")


def r_has(*pkgs):
    if not shutil.which("Rscript"):
        return False
    expr = "quit(status = if (all(vapply(c(%s), requireNamespace, logical(1), quietly = TRUE))) 0 else 1)" \
        % ", ".join(f'"{p}"' for p in pkgs)
    return subprocess.run(["Rscript", "-e", expr], capture_output=True).returncode == 0


def rscript(expr):
    p = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError(p.stderr[-1500:])
    return p.stdout


# 6 subjects x 2 conditions. Each protein has a subject effect (sd 0.5) larger than the
# run-to-run noise (sd 0.25): paired differences are precise, unpaired ones are not.
# P00001-P00020 go up 0.5 log2 in the second condition, P00021-P00040 down 0.5 (balanced,
# so normalisation does not shift the unchanged proteins). `nested` relabels the same
# layout as subject-within-age: S1-S3 Old, S4-S6 Young, each with one IgG and one Bait IP.
SYNTH_R = r'''
make <- function(OUT, nested) {
  set.seed(11)
  nsub <- 6; runs <- sprintf("run%02d", 1:12)
  subj <- rep(sprintf("S%d", 1:nsub), times = 2); cond <- rep(c("A", "B"), each = nsub)
  grp <- if (nested) paste0(ifelse(subj %in% c("S1", "S2", "S3"), "Old", "Young"), "_",
                            ifelse(cond == "A", "IgG", "Bait")) else cond
  out <- list(); pgq <- list()
  for (i in 1:150) {
    pg <- sprintf("P%05d", i); base <- runif(1, 12, 20)
    s_eff <- setNames(rnorm(nsub, 0, 0.5), sprintf("S%d", 1:nsub))
    eff <- if (i <= 20) 0.5 else if (i <= 40) -0.5 else 0
    run_lv <- base + s_eff[subj] + (cond == "B") * eff + rnorm(12, 0, 0.25)
    # PG.MaxLFQ: the protein level itself (MaxLFQ is built to survive missing precursors)
    pgq[[i]] <- data.frame(Protein.Group = pg, Run = runs, PG.MaxLFQ = 2^run_lv)
    for (j in 1:4) {
      lv <- run_lv + rnorm(1, 0, 1)
      int <- 2^(lv + rnorm(12, 0, 0.1))
      keep <- runif(12) >= plogis(-(lv - 12) * 2)
      if (!any(keep)) next
      out[[length(out) + 1]] <- data.frame(Run = runs[keep], Precursor.Id = sprintf("PEP%dK%d2", i, j),
        Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0("G", i, "_X"),
        Genes = paste0("G", i), Proteotypic = 1L, Precursor.Normalised = int[keep],
        Precursor.Quantity = int[keep], Q.Value = 0.001, Lib.Q.Value = 0.001,
        Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001, Global.Q.Value = 0.001,
        Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
    }
  }
  d <- do.call(rbind, out)
  d <- merge(d, do.call(rbind, pgq), by = c("Protein.Group", "Run"))
  dir.create(OUT, showWarnings = FALSE)
  arrow::write_parquet(d, file.path(OUT, "report.parquet"))
  meta <- data.frame(File.Name = runs, Group = grp, Subject = subj)
  write.csv(meta, file.path(OUT, "conditions.csv"), row.names = FALSE)
  # the same column as a fixed covariate (run_de.R fits Covariate1 as one)
  write.csv(transform(meta, Covariate1 = Subject), file.path(OUT, "conditions_cov.csv"), row.names = FALSE)
  # a block that is just the grouping under another name
  write.csv(transform(meta, Cage = Group), file.path(OUT, "conditions_cage.csv"), row.names = FALSE)
  # one subject left with a single sample
  m1 <- meta; m1$Subject[12] <- "S7"
  write.csv(m1, file.path(OUT, "conditions_single.csv"), row.names = FALSE)
}
make(file.path(OUT, "paired"), FALSE)
make(file.path(OUT, "nested"), TRUE)
'''

TRUE_HITS = {f"P{i:05d}" for i in range(1, 41)}
# within-subject (bait vs IgG), between with one sample per subject (Old vs Young for one
# IP type), and between POOLING two samples per subject (Old vs Young over both IP types --
# unblocked, that one would be pseudo-replicated)
POOLED = "(Old_Bait+Old_IgG)/2-(Young_Bait+Young_IgG)/2"
NESTED_CONTRASTS = ",".join(["Old_Bait-Old_IgG", "Young_Bait-Young_IgG", "Old_Bait-Young_Bait",
                             "Old_IgG-Young_IgG", POOLED])
NESTED_MODEL = {"Old_Bait-Old_IgG": "blocked", "Young_Bait-Young_IgG": "blocked",
                "Old_Bait-Young_Bait": "independent", "Old_IgG-Young_IgG": "independent",
                POOLED: "blocked"}


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def significant(path, adjp=0.05):
    return {r["Protein.Group"] for r in read_csv(path)
            if r["adj.P.Val"] not in ("", "NA") and float(r["adj.P.Val"]) < adjp}


@unittest.skipUnless(r_has("limpa", "limma", "statmod", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limpa/limma/statmod/arrow/dplyr/tidyr/jsonlite")
class RunDeBlock(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        cls.runs = {}
        plans = {
            "dpc":         ("paired", "conditions.csv", ["--method", "dpc"]),
            "dpc_block":   ("paired", "conditions.csv", ["--method", "dpc", "--block", "Subject"]),
            "ml":          ("paired", "conditions.csv", ["--method", "maxlfq"]),
            "ml_block":    ("paired", "conditions.csv", ["--method", "maxlfq", "--block", "Subject"]),
            # --block-scope within (the default), all, and no block at all -- same contrasts
            "nested":      ("nested", "conditions.csv", ["--method", "dpc", "--block", "Subject",
                                                         "--contrasts", NESTED_CONTRASTS]),
            "nested_all":  ("nested", "conditions.csv", ["--method", "dpc", "--block", "Subject",
                                                         "--block-scope", "all",
                                                         "--contrasts", NESTED_CONTRASTS]),
            "nested_plain": ("nested", "conditions.csv", ["--method", "dpc",
                                                          "--contrasts", NESTED_CONTRASTS]),
            "nested_ml":   ("nested", "conditions.csv", ["--method", "maxlfq", "--block", "Subject",
                                                         "--contrasts", NESTED_CONTRASTS]),
            "nested_ml_all": ("nested", "conditions.csv", ["--method", "maxlfq", "--block", "Subject",
                                                           "--block-scope", "all",
                                                           "--contrasts", NESTED_CONTRASTS]),
            "nested_ml_plain": ("nested", "conditions.csv", ["--method", "maxlfq",
                                                             "--contrasts", NESTED_CONTRASTS]),
            "scope_alone": ("paired", "conditions.csv", ["--method", "dpc", "--block-scope", "all"]),
            "scope_bad":   ("paired", "conditions.csv", ["--method", "dpc", "--block", "Subject",
                                                         "--block-scope", "between"]),
            "conflict":    ("nested", "conditions_cov.csv", ["--method", "dpc", "--block", "Covariate1"]),
            # the subject as a FIXED covariate, nested in the groups: rank-deficient
            "fixed":       ("nested", "conditions_cov.csv", ["--method", "dpc"]),
            "group":       ("paired", "conditions.csv", ["--method", "dpc", "--block", "Group"]),
            "encoded":     ("paired", "conditions_cage.csv", ["--method", "dpc", "--block", "Cage"]),
            "single":      ("paired", "conditions_single.csv", ["--method", "dpc", "--block", "Subject"]),
        }
        for key, (ds, meta, args) in plans.items():
            d = os.path.join(cls.tmp, ds)
            out = os.path.join(cls.tmp, "out_" + key)
            p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(d, "report.parquet"),
                                "--metadata", os.path.join(d, meta), "--outdir", out, *args],
                               capture_output=True, text=True, cwd=d)
            cls.runs[key] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        return out

    def prov(self, key):
        with open(os.path.join(self.out(key), "de_provenance.json")) as fh:
            return json.load(fh)

    def methods(self, key):
        with open(os.path.join(self.out(key), "methods.txt")) as fh:
            return fh.read()

    def failed(self, key):
        p, _ = self.runs[key]
        self.assertNotEqual(p.returncode, 0, f"{key} should have stopped")
        # the design checks run before quantification: limpa never got as far as dpcQuant
        self.assertNotIn("Quantifying proteins", p.stdout + p.stderr)
        return p.stderr

    # -- the point of it ----------------------------------------------------------------
    def test_blocking_increases_power(self):
        for plain, blocked, de in (("dpc", "dpc_block", "DE_dpc_B.A.csv"),
                                   ("ml", "ml_block", "DE_maxlfq_B.A.csv")):
            with self.subTest(method=plain):
                s0 = significant(os.path.join(self.out(plain), de))
                s1 = significant(os.path.join(self.out(blocked), de))
                self.assertGreater(len(s1 & TRUE_HITS), len(s0 & TRUE_HITS) + 10,
                                   f"unblocked {len(s0)}, blocked {len(s1)}")
                # the gain is real changes, not a looser test
                self.assertLessEqual(len(s1 - TRUE_HITS), max(2, len(s1) // 10))

    def test_logfc_is_the_same_estimate(self):
        # complete, balanced pairs: GLS and OLS give the same group means (maxlfq is
        # unweighted), so blocking changes the standard errors, not the fold changes
        complete = {r["Protein.Group"] for r in
                    read_csv(os.path.join(self.out("ml"), "Expression_Matrix.csv"))
                    if all(v not in ("", "NA") for k, v in r.items() if k.startswith("run"))}
        a = {r["Protein.Group"]: float(r["logFC"]) for r in
             read_csv(os.path.join(self.out("ml"), "DE_maxlfq_B.A.csv"))}
        b = {r["Protein.Group"]: float(r["logFC"]) for r in
             read_csv(os.path.join(self.out("ml_block"), "DE_maxlfq_B.A.csv"))}
        shared = [k for k in a if k in b and k in complete]
        self.assertGreater(len(shared), 100)
        self.assertLess(max(abs(a[k] - b[k]) for k in shared), 1e-6)

    # -- the record ---------------------------------------------------------------------
    def test_provenance_records_the_block(self):
        for key, fit in (("dpc_block", "limpa::dpcDE(block =)"),
                         ("ml_block", "limma::duplicateCorrelation(E, design, block =)")):
            with self.subTest(run=key):
                p = self.prov(key)
                b = p["block"]
                self.assertTrue(b["applied"])
                self.assertEqual(b["column"], "Subject")
                self.assertEqual(p["block_column"], "Subject")   # top-level, for simple readers
                self.assertEqual(b["n_blocks"], 6)
                self.assertEqual(b["block_sizes"], {f"S{i}": 2 for i in range(1, 7)})
                self.assertGreater(b["consensus_correlation"], 0.5)
                self.assertLess(b["consensus_correlation"], 1)
                self.assertEqual(b["n_proteins"], 150)
                self.assertGreater(b["n_proteins_estimated"], 100)
                self.assertTrue(b["fit"].startswith(fit), b["fit"])
                self.assertIn("duplicateCorrelation", b["estimator"])
                self.assertEqual(b["contrast_structure"], {"B-A": "within"})
                self.assertEqual(b["scope"], "within")
                self.assertEqual(b["contrast_model"], {"B-A": "blocked"})
                self.assertEqual(p["de_tables"]["B-A"]["model"], "blocked")
                self.assertEqual(b["warnings"], [])
                # the fixed part of the model is unchanged: Subject is not a design term
                self.assertEqual(p["design"], "~ 0 + groups")
                self.assertIn("block = Subject", p["de_engine"])
                self.assertEqual(p["contrasts"], ["B-A"])   # a list even when there is one
        # dpc: limpa estimates twice; both are kept
        self.assertIsInstance(self.prov("dpc_block")["block"]["first_pass_correlation"], float)

    def test_unblocked_run_says_independent(self):
        b = self.prov("dpc")["block"]
        self.assertFalse(b["applied"])
        self.assertNotIn("block_column", self.prov("dpc"))   # present only when blocked
        self.assertIn("Blocking      : none", self.methods("dpc"))
        self.assertNotIn("block", self.prov("dpc")["de_engine"])
        # ... and points at the column that looks like a blocking unit
        self.assertIn("pass --block Subject", self.runs["dpc"][0].stderr)

    def test_methods_txt_states_the_block(self):
        rho = self.prov("dpc_block")["block"]["consensus_correlation"]
        txt = self.methods("dpc_block")
        self.assertIn("Blocking      : Subject as a random effect (6 levels, 2 samples each)", txt)
        self.assertIn(f"consensus within-Subject correlation {rho:.3f}", txt)
        self.assertIn("Within-Subject contrasts, blocked fit: B-A", txt)
        self.assertIn("Scope: all contrasts from the blocked fit (scope within)", txt)
        self.assertNotIn("CAUTION", txt)

    def test_make_methods_sentence_reads_the_record(self):
        p = self.prov("dpc_block")
        para = mm.de_paragraph(p)
        self.assertIn("Subject (6 levels) was fitted as a random blocking factor", para)
        self.assertIn(f"{p['block']['consensus_correlation']:.3f}", para)
        self.assertIn("limpa::dpcDE(block =)", para)
        self.assertIn("with contrasts B-A.", para)
        self.assertIn("All contrasts were reported from this fit.", para)
        self.assertNotIn("blocking factor", mm.de_paragraph(self.prov("dpc")))
        # a record from before --block existed: no sentence, nothing assumed
        old = dict(p); old.pop("block")
        self.assertNotIn("blocking factor", mm.de_paragraph(old))

    # -- nested: subject within age ------------------------------------------------------
    def test_nested_design_fits_and_labels_contrasts(self):
        p = self.prov("nested")
        self.assertEqual(p["design"], "~ 0 + groups")
        self.assertEqual(p["block"]["contrast_structure"],
                         {"Old_Bait-Old_IgG": "within", "Young_Bait-Young_IgG": "within",
                          "Old_Bait-Young_Bait": "between", "Old_IgG-Young_IgG": "between",
                          POOLED: "between"})
        self.assertGreater(p["block"]["consensus_correlation"], 0.5)
        txt = self.methods("nested")
        self.assertIn("Within-Subject contrasts, blocked fit: Old_Bait-Old_IgG, Young_Bait-Young_IgG", txt)
        self.assertIn("Between-Subject contrasts, independent fit (judged on the Subject levels, "
                      "not the samples): Old_Bait-Young_Bait, Old_IgG-Young_IgG", txt)
        self.assertIn(f"Between-Subject contrasts, blocked fit (judged on the Subject levels, "
                      f"not the samples): {POOLED}", txt)
        self.assertIn("Scope: between-Subject contrasts using at most one sample per Subject from "
                      "the fit with samples independent, all others from the blocked fit "
                      "(scope within)", txt)
        for t in p["de_tables"].values():
            self.assertTrue(os.path.exists(os.path.join(self.out("nested"), t["file"])))

    # -- --block-scope: one run, each contrast from its fit ----------------------------------
    def test_scope_within_picks_the_fit_per_contrast(self):
        for key in ("nested", "nested_ml"):
            with self.subTest(run=key):
                p = self.prov(key)
                b = p["block"]
                self.assertEqual(b["scope"], "within")
                self.assertTrue(b["applied"])                 # it still reports within contrasts
                self.assertEqual(p["block_column"], "Subject")
                self.assertEqual(b["contrast_model"], NESTED_MODEL)
                self.assertEqual({cn: t["model"] for cn, t in p["de_tables"].items()}, NESTED_MODEL)
                self.assertIn("between-Subject contrasts using at most one sample per Subject "
                              "from the fit with samples independent", p["de_engine"])
                self.assertIn("at most one sample per Subject", b["contrast_model_rule"])
        all_ = self.prov("nested_all")["block"]
        self.assertEqual(all_["scope"], "all")
        self.assertEqual(set(all_["contrast_model"].values()), {"blocked"})
        self.assertEqual({t["model"] for t in self.prov("nested_plain")["de_tables"].values()},
                         {"independent"})

    def test_scope_within_tables_equal_the_separate_runs(self):
        cols = ("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
        for key, blocked, plain in (("nested", "nested_all", "nested_plain"),
                                    ("nested_ml", "nested_ml_all", "nested_ml_plain")):
            p = self.prov(key)
            for cn, t in p["de_tables"].items():
                with self.subTest(run=key, contrast=cn, model=t["model"]):
                    ref = blocked if t["model"] == "blocked" else plain
                    a = {r["Protein.Group"]: r for r in read_csv(os.path.join(self.out(key), t["file"]))}
                    b = {r["Protein.Group"]: r for r in read_csv(os.path.join(self.out(ref), t["file"]))}
                    self.assertEqual(set(a), set(b))
                    for k in a:
                        for c in cols:
                            self.assertAlmostEqual(float(a[k][c]), float(b[k][c]), places=8,
                                                   msg=f"{k} {c}")

    def test_make_methods_names_the_fit_per_contrast(self):
        para = mm.de_paragraph(self.prov("nested"))
        self.assertIn(f"This fit reported Old_Bait-Old_IgG, Young_Bait-Young_IgG, {POOLED};", para)
        self.assertIn("Old_Bait-Young_Bait, Old_IgG-Young_IgG -- contrasts between different "
                      "Subject levels using at most one sample per Subject -- were reported from "
                      "the same data fitted with samples as independent", para)
        # a blocked record from before --block-scope: which fit reported what is not assumed
        old = self.prov("dpc_block"); old["block"].pop("contrast_model")
        self.assertIn("not recorded", mm.de_paragraph(old))

    def test_block_scope_arguments(self):
        self.assertIn("--block-scope needs --block", self.failed("scope_alone"))
        self.assertIn("--block-scope must be 'within' (default) or 'all'", self.failed("scope_bad"))

    # -- what must stop ------------------------------------------------------------------
    def test_block_that_is_also_a_covariate_stops(self):
        err = self.failed("conflict")
        self.assertIn("--block Covariate1", err)
        self.assertIn("also a fixed-effect covariate", err)

    def test_nested_subject_as_fixed_covariate_points_at_block(self):
        err = self.failed("fixed")
        self.assertIn("Design matrix is not full rank", err)
        self.assertIn("Covariate1 is confounded with the groups but recurs across them", err)
        self.assertIn("--block", err)

    def test_block_that_is_the_grouping_stops(self):
        self.assertIn("fixed grouping itself", self.failed("group"))
        err = self.failed("encoded")
        self.assertIn("--block Cage is already a fixed effect in the design", err)

    def test_single_sample_block_stops(self):
        err = self.failed("single")
        self.assertIn("block(s) hold a single analysed sample", err)
        self.assertIn("S6", err)
        self.assertIn("S7", err)

    # -- reproducible ----------------------------------------------------------------------
    def test_reproducibility_log_refits_the_blocked_model(self):
        # paired: one blocked fit; nested: the blocked fit AND the independent one (scope within)
        for key in ("dpc_block", "ml_block", "nested", "nested_ml"):
            with self.subTest(run=key):
                out = self.out(key)
                with open(os.path.join(out, "reproducibility_log.R")) as fh:
                    src = fh.read()
                self.assertIn("block <- metadata[['Subject']]", src)
                self.assertIn("block = block", src)
                if key.startswith("nested"):
                    self.assertIn("fit_independent <- limma::eBayes(", src)
                    self.assertIn("'Old_Bait-Young_Bait' = 'independent'", src)
                p = subprocess.run(["Rscript", "reproducibility_log.R"], capture_output=True,
                                   text=True, cwd=out)
                self.assertEqual(p.returncode, 0, p.stderr[-2000:])
                for t in self.prov(key)["de_tables"].values():
                    a = {r["Protein.Group"]: r for r in read_csv(os.path.join(out, t["file"]))}
                    b = {r["Protein.Group"]: r for r in
                         read_csv(os.path.join(out, "de_results_rerun", t["file"]))}
                    self.assertEqual(set(a), set(b))
                    for k in a:
                        for col in ("logFC", "P.Value"):
                            self.assertAlmostEqual(float(a[k][col]), float(b[k][col]), places=8)


@unittest.skipUnless(r_has("limma", "jsonlite"), "needs Rscript + limma/jsonlite")
class BlockRecordWarnings(unittest.TestCase):
    """blocking.R's record, directly: what it warns about, and the contrast labels."""

    @classmethod
    def setUpClass(cls):
        cls.res = json.loads(rscript(f'''
          source("{BLOCK_R}")
          b6 <- rep(sprintf("M%d", 1:6), times = 2); b2 <- rep(c("M1", "M2"), times = 3)
          many <- atanh(rep(0.4, 500))
          r <- list(
            ok       = block_record("Mouse", b6, "maxlfq", 0.4, many, 500),
            negative = block_record("Mouse", b6, "maxlfq", -0.05, many, 500),
            few      = block_record("Mouse", b6, "maxlfq", 0.4, atanh(rep(0.4, 10)), 500),
            moving   = block_record("Mouse", b6, "dpc", 0.45, many, 500, first_pass = 0.2),
            two      = block_record("Mouse", b2, "maxlfq", 0.4, many, 500))
          warn <- lapply(r, function(x) x$warnings)   # a list stays a JSON array
          g <- factor(c("Old_IgG", "Old_Bait", "Young_IgG", "Young_Bait", "Old_IgG", "Old_Bait"))
          blk <- c("M1", "M1", "M2", "M2", "M3", "M4")
          design <- model.matrix(~ 0 + g); colnames(design) <- levels(g)
          cm <- limma::makeContrasts(contrasts = c("Old_Bait-Old_IgG", "Old_Bait-Young_Bait"),
                                     levels = design)
          st <- block_contrast_structure(blk, g, cm)
          # the reporting-fit rule: 3 mice per age, one Bait + one IgG IP each
          g2 <- factor(rep(c("Old_IgG", "Old_Bait", "Young_IgG", "Young_Bait"), each = 3))
          b2 <- c("M1", "M2", "M3", "M1", "M2", "M3", "M4", "M5", "M6", "M4", "M5", "M6")
          d2 <- model.matrix(~ 0 + g2); colnames(d2) <- levels(g2)
          cm2 <- limma::makeContrasts(contrasts = c("Old_Bait-Old_IgG", "Old_Bait-Young_Bait",
            "(Old_Bait+Old_IgG)/2-(Young_Bait+Young_IgG)/2"), levels = d2)
          b3 <- b2; b3[2] <- "M1"   # M1 now holds two Old_IgG samples: a technical replicate
          models <- list(within = block_contrast_model(b2, g2, cm2, "within"),
                         all = block_contrast_model(b2, g2, cm2, "all"),
                         reps = block_contrast_model(b3, g2, cm2, "within"))
          cat(jsonlite::toJSON(list(warn = warn, structure = st, n_ok = length(r$ok$warnings),
                                    models = models), auto_unbox = TRUE))'''))

    def test_clean_estimate_has_no_warning(self):
        self.assertEqual(self.res["n_ok"], 0)

    def test_correlation_at_or_below_zero_warns(self):
        self.assertTrue(any("(<= 0)" in w for w in self.res["warn"]["negative"]))

    def test_unstable_estimates_warn(self):
        self.assertTrue(any("unstable" in w and "only 10 protein" in w
                            for w in self.res["warn"]["few"]))
        self.assertTrue(any("unstable" in w and "two estimation passes" in w
                            for w in self.res["warn"]["moving"]))
        self.assertTrue(any("only 2 Mouse levels" in w for w in self.res["warn"]["two"]))

    def test_reporting_fit_rule(self):
        m = self.res["models"]
        self.assertEqual(m["within"], {"Old_Bait-Old_IgG": "blocked", "Old_Bait-Young_Bait": "independent",
                                       "(Old_Bait+Old_IgG)/2-(Young_Bait+Young_IgG)/2": "blocked"})
        self.assertEqual(set(m["all"].values()), {"blocked"})
        # a subject with two samples in one group pseudo-replicates the independent fit too
        self.assertEqual(set(m["reps"].values()), {"blocked"})

    def test_contrast_structure(self):
        # M1 has both an Old IgG and an Old Bait IP, M3/M4 only one each -> partial
        self.assertEqual(self.res["structure"],
                         {"Old_Bait-Old_IgG": "partial", "Old_Bait-Young_Bait": "between"})


if __name__ == "__main__":
    unittest.main()
