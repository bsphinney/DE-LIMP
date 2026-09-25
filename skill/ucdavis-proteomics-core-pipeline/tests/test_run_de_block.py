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
NESTED_CONTRASTS = "Old_Bait-Old_IgG,Young_Bait-Young_IgG,Old_Bait-Young_Bait,Old_IgG-Young_IgG"


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
            "nested":      ("nested", "conditions.csv", ["--method", "dpc", "--block", "Subject",
                                                         "--contrasts", NESTED_CONTRASTS]),
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
        self.assertIn("Within-Subject contrasts: B-A", txt)
        self.assertNotIn("CAUTION", txt)

    def test_make_methods_sentence_reads_the_record(self):
        p = self.prov("dpc_block")
        para = mm.de_paragraph(p)
        self.assertIn("Subject (6 levels) was fitted as a random blocking factor", para)
        self.assertIn(f"{p['block']['consensus_correlation']:.3f}", para)
        self.assertIn("limpa::dpcDE(block =)", para)
        self.assertIn("with contrasts B-A.", para)
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
                          "Old_Bait-Young_Bait": "between", "Old_IgG-Young_IgG": "between"})
        self.assertGreater(p["block"]["consensus_correlation"], 0.5)
        txt = self.methods("nested")
        self.assertIn("Within-Subject contrasts: Old_Bait-Old_IgG, Young_Bait-Young_IgG", txt)
        self.assertIn("Between-Subject contrasts (judged on the Subject levels, not the samples): "
                      "Old_Bait-Young_Bait, Old_IgG-Young_IgG", txt)
        for cn in ("Old_Bait.Old_IgG", "Old_Bait.Young_Bait"):
            self.assertTrue(os.path.exists(os.path.join(self.out("nested"), f"DE_dpc_{cn}.csv")))

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
        for key, de in (("dpc_block", "DE_dpc_B.A.csv"), ("ml_block", "DE_maxlfq_B.A.csv")):
            with self.subTest(run=key):
                out = self.out(key)
                with open(os.path.join(out, "reproducibility_log.R")) as fh:
                    src = fh.read()
                self.assertIn("block <- metadata[['Subject']]", src)
                self.assertIn("block = block", src)
                p = subprocess.run(["Rscript", "reproducibility_log.R"], capture_output=True,
                                   text=True, cwd=out)
                self.assertEqual(p.returncode, 0, p.stderr[-2000:])
                a = {r["Protein.Group"]: r for r in read_csv(os.path.join(out, de))}
                b = {r["Protein.Group"]: r for r in
                     read_csv(os.path.join(out, "de_results_rerun", de))}
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
          cat(jsonlite::toJSON(list(warn = warn, structure = st, n_ok = length(r$ok$warnings)),
                               auto_unbox = TRUE))'''))

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

    def test_contrast_structure(self):
        # M1 has both an Old IgG and an Old Bait IP, M3/M4 only one each -> partial
        self.assertEqual(self.res["structure"],
                         {"Old_Bait-Old_IgG": "partial", "Old_Bait-Young_Bait": "between"})


if __name__ == "__main__":
    unittest.main()
