#!/usr/bin/env python3
"""On/off proteins (2.15.x, MaxLFQ route): onoff.R's compute_onoff_proteins() is DE-LIMP's, and
the two must not drift (DE-LIMP architectural rule 3: a concept has ONE definition).

  * parity: the same fixtures through DE-LIMP's R/helpers.R function and the skill's copy in
    onoff.R give identical results, and the two functions deparse to the same code. Skipped when
    this skill is installed without the DE-LIMP repo around it (R/helpers.R absent).
  * the rule itself: one group missing entirely, n_min, orientation (g1 - g2), NA genes, zero
    results, NA group labels, a protein in 1 run only.
  * the contrast mapper: A-B and -B+A are simple; averages and interactions are skipped WITH a
    reason; the table is always a data.frame (zero rows when nothing is on/off).

Needs only base R (no limma/limpa). NOTE: CI runs Python only and skips R tests -- run locally.
"""
import json
import os
import shutil
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
ONOFF_R = os.path.join(SCRIPTS, "onoff.R")
HELPERS_R = os.path.normpath(os.path.join(HERE, "..", "..", "..", "R", "helpers.R"))

HAVE_R = bool(shutil.which("Rscript"))


def rscript(code):
    p = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError(p.stderr[-2500:])
    return p.stdout


# 8 proteins x 7 runs: groups A (3 runs), B (3 runs), C (1 run), plus one run with NA group.
# rows (A: 3 runs, B: 3 runs, C: 1 run, NA-group run):
#   p_Aonly   A 3/3, B 0/3            -> A_only
#   p_Bonly   A 0/3, B 2/3            -> B_only (n_min 2); 1/3 would need n_min 1
#   p_B1      A 0/3, B 1/3            -> only at n_min 1
#   p_both    A 3/3, B 3/3            -> never
#   p_part    A 1/3, B 2/3            -> never (present on both sides)
#   p_none    nowhere                 -> never
#   p_Conly   C 1/1 only              -> with C vs A: C_only at n_min 1 only
#   p_NAgrp   only in the run whose group is NA -> never
FIXTURE_R = r'''
cols <- c("a1","a2","a3","b1","b2","b3","c1","x1")
E <- matrix(NA_real_, 8, 8, dimnames = list(
  c("p_Aonly","p_Bonly","p_B1","p_both","p_part","p_none","p_Conly","p_NAgrp"), cols))
E["p_Aonly", c("a1","a2","a3")] <- 20
E["p_Bonly", c("b1","b2")] <- 21
E["p_B1", "b3"] <- 22
E["p_both", c("a1","a2","a3","b1","b2","b3")] <- 19
E["p_part", c("a1","b2","b3")] <- 18
E["p_Conly", "c1"] <- 17
E["p_NAgrp", "x1"] <- 16
grp <- c("A","A","A","B","B","B","C", NA)
genes <- c(p_Aonly = "GA", p_Bonly = NA, p_B1 = "GB1", p_both = "GBOTH", p_part = "GP",
           p_none = "GN", p_Conly = "GC", p_NAgrp = "GX")
'''

# the scenarios both functions are run through: (label, call args as R text)
SCENARIOS = {
    "pairs_AB_n2": 'list(c("A","B")), n_min = 2',
    "pairs_BA_n2": 'list(c("B","A")), n_min = 2',
    "pairs_AB_n1": 'list(c("A","B")), n_min = 1',
    "pairs_AB_n3": 'list(c("A","B")), n_min = 3',
    "matrix_two_contrasts": 'cbind(c("A","B"), c("B","A")), n_min = 2',
    "all_pairs_default": 'NULL, n_min = 1',
    "one_group_missing_entirely": 'list(c("A","Z")), n_min = 1',
    "with_gene_lookup": 'list(c("A","B")), n_min = 1, gene_lookup = genes',
    "three_pairs": 'list(c("A","B"), c("C","A"), c("B","C")), n_min = 1',
    "zero_results": 'list(c("A","B")), n_min = 99',
}


@unittest.skipUnless(HAVE_R, "needs Rscript")
class OnOffRule(unittest.TestCase):
    """The rule, through the skill's copy."""

    @classmethod
    def setUpClass(cls):
        cls.res = {}
        body = FIXTURE_R + f'source("{ONOFF_R}")\n'
        body += "out <- list()\n"
        for k, args in SCENARIOS.items():
            body += f'out[["{k}"]] <- compute_onoff_proteins(E, grp, {args.replace("NULL, ", "contrasts_list = NULL, ", 1) if args.startswith("NULL") else "contrasts_list = " + args})\n'
        body += 'cat(jsonlite::toJSON(out, auto_unbox = TRUE, null = "null", na = "null", dataframe = "rows"))'
        cls.has_jsonlite = subprocess.run(["Rscript", "-e", "quit(status=!requireNamespace('jsonlite',quietly=TRUE))"]).returncode == 0
        if cls.has_jsonlite:
            cls.res = json.loads(rscript(body))

    def need(self):
        if not self.has_jsonlite:
            self.skipTest("needs jsonlite")

    def test_orientation_is_first_minus_second_and_direction_names_the_group(self):
        self.need()
        ab = self.res["pairs_AB_n2"]
        self.assertEqual({(r["Protein.Group"], r["Contrast"], r["Direction"]) for r in ab},
                         {("p_Aonly", "A - B", "A_only"), ("p_Bonly", "A - B", "B_only")})
        ba = self.res["pairs_BA_n2"]
        self.assertEqual({(r["Protein.Group"], r["Contrast"], r["Direction"]) for r in ba},
                         {("p_Aonly", "B - A", "A_only"), ("p_Bonly", "B - A", "B_only")})
        # detected_g1 is the FIRST group of the contrast
        r = {x["Protein.Group"]: x for x in ba}["p_Aonly"]
        self.assertEqual((r["detected_g1"], r["detected_g2"], r["total_g1"], r["total_g2"]), (0, 3, 3, 3))

    def test_n_min_is_the_runs_needed_in_the_group_that_has_it(self):
        self.need()
        n1 = {r["Protein.Group"] for r in self.res["pairs_AB_n1"]}
        self.assertEqual(n1, {"p_Aonly", "p_Bonly", "p_B1"})
        n2 = {r["Protein.Group"] for r in self.res["pairs_AB_n2"]}
        self.assertEqual(n2, {"p_Aonly", "p_Bonly"})
        n3 = {r["Protein.Group"] for r in self.res["pairs_AB_n3"]}
        self.assertEqual(n3, {"p_Aonly"})                    # B has only 2 runs for p_Bonly

    def test_protein_on_both_sides_or_nowhere_is_never_listed(self):
        self.need()
        # A vs B scenarios only (C has one run, so p_both is "A only" against C, correctly)
        every = {r["Protein.Group"] for k in SCENARIOS if k.startswith("pairs_")
                 and isinstance(self.res.get(k), list) for r in self.res[k]}
        for pg in ("p_both", "p_part", "p_none", "p_NAgrp"):
            self.assertNotIn(pg, every)

    def test_a_group_missing_entirely_gives_nothing_not_everything(self):
        self.need()
        # no run is labelled Z: the contrast is skipped (not "A has it, Z lacks it" for all of A)
        self.assertIsNone(self.res.get("one_group_missing_entirely"))   # NULL: key absent in JSON

    def test_zero_results_is_null(self):
        self.need()
        self.assertIsNone(self.res.get("zero_results"))

    def test_default_contrasts_are_all_pairs_second_minus_first(self):
        self.need()
        con = {r["Contrast"] for r in self.res["all_pairs_default"]}
        self.assertEqual(con, {"B - A", "C - A", "C - B"})

    def test_na_gene_stays_na_and_is_not_dropped(self):
        self.need()
        rows = {r["Protein.Group"]: r for r in self.res["with_gene_lookup"]}
        self.assertEqual(rows["p_Aonly"]["Gene"], "GA")
        self.assertIsNone(rows["p_Bonly"]["Gene"])           # NA gene: the protein is still there
        self.assertEqual(rows["p_B1"]["Gene"], "GB1")

    def test_a_run_with_no_group_label_counts_for_no_group(self):
        self.need()
        # p_NAgrp is quantified only in the unlabelled run
        self.assertNotIn("p_NAgrp", {r["Protein.Group"] for r in self.res["three_pairs"]})


@unittest.skipUnless(HAVE_R, "needs Rscript")
class ContrastMapperAndTable(unittest.TestCase):
    def run_r(self, code):
        out = rscript(f'source("{ONOFF_R}")\n' + code)
        return json.loads(out)

    def setUp(self):
        if subprocess.run(["Rscript", "-e", "quit(status=!requireNamespace('jsonlite',quietly=TRUE))"]).returncode:
            self.skipTest("needs jsonlite")

    def test_simple_contrasts_map_to_pairs_complex_ones_are_skipped_with_a_reason(self):
        r = self.run_r(r'''
          lv <- c("A","B","C","D")
          cm <- cbind("B-A" = c(-1,1,0,0), "-B+A" = c(1,-1,0,0),
                      "C-(A+B)/2" = c(-.5,-.5,1,0), "(B-A)-(D-C)" = c(-1,1,1,-1),
                      "2B-2A" = c(-2,2,0,0))
          rownames(cm) <- lv
          x <- onoff_contrast_pairs(cm, lv)
          cat(jsonlite::toJSON(x, auto_unbox = TRUE))''')
        self.assertEqual(r["pairs"], {"B-A": ["B", "A"], "-B+A": ["A", "B"]})
        self.assertEqual(set(r["skipped"]), {"C-(A+B)/2", "(B-A)-(D-C)", "2B-2A"})
        for why in r["skipped"].values():
            self.assertIn("not a simple A-B contrast", why)
        self.assertIn("3 terms", r["skipped"]["C-(A+B)/2"])

    def test_a_covariate_term_makes_a_contrast_not_simple(self):
        r = self.run_r(r'''
          cm <- cbind("B-A" = c(-1,1,0)); rownames(cm) <- c("A","B","age")
          x <- onoff_contrast_pairs(cm, c("A","B"))
          cm2 <- cbind("B-age" = c(0,1,-1)); rownames(cm2) <- c("A","B","age")
          y <- onoff_contrast_pairs(cm2, c("A","B"))
          cat(jsonlite::toJSON(list(x = x, y = y), auto_unbox = TRUE))''')
        self.assertEqual(r["x"]["pairs"], {"B-A": ["B", "A"]})
        self.assertEqual(list(r["y"]["pairs"]), [])
        self.assertIn("B-age", r["y"]["skipped"])

    def test_table_columns_genes_and_zero_rows(self):
        r = self.run_r(FIXTURE_R + r'''
          annot <- data.frame(Protein.Group = names(genes), Genes = unname(genes))
          t1 <- onoff_table(E, grp, list("A-B" = c("A","B")), n_min = 2, annot = annot)
          t0 <- onoff_table(E, grp, list("A-B" = c("A","B")), n_min = 99, annot = annot)
          tn <- onoff_table(E, grp, list(), n_min = 2, annot = annot)
          t_ng <- onoff_table(E, grp, list("A-B" = c("A","B")), n_min = 2, annot = NULL)
          cat(jsonlite::toJSON(list(t1 = t1, t0_n = nrow(t0), t0_cols = names(t0), tn_n = nrow(tn),
                                    t_ng = t_ng), auto_unbox = TRUE, na = "null", dataframe = "rows"))''')
        cols = ["Protein.Group", "Genes", "Contrast", "Direction",
                "detected_g1", "detected_g2", "total_g1", "total_g2"]
        self.assertEqual(list(r["t1"][0]), cols)
        rows = {x["Protein.Group"]: x for x in r["t1"]}
        self.assertEqual(rows["p_Aonly"]["Genes"], "GA")
        self.assertIsNone(rows["p_Bonly"]["Genes"])
        self.assertEqual((r["t0_n"], r["t0_cols"]), (0, cols))     # zero rows, same columns
        self.assertEqual(r["tn_n"], 0)                             # no simple contrast at all
        self.assertTrue(all(x["Genes"] is None for x in r["t_ng"]))  # no annotation: NA, not dropped
        self.assertEqual(len(r["t_ng"]), 2)

    def test_rule_text_names_n_min_and_says_not_a_fold_change(self):
        r = self.run_r('cat(jsonlite::toJSON(list(a = onoff_rule_text(2), b = onoff_rule_text(1)), auto_unbox = TRUE))')
        self.assertIn("at least 2 runs", r["a"])
        self.assertIn("at least 1 run ", r["b"])
        for t in r.values():
            self.assertIn("not a measured ratio", t)
            self.assertIn("before the coverage filter", t)


@unittest.skipUnless(HAVE_R and os.path.isfile(HELPERS_R),
                     "needs Rscript and DE-LIMP's R/helpers.R (this skill is not inside the repo)")
class ParityWithDeLimp(unittest.TestCase):
    """The skill's copy and DE-LIMP's function agree: same code, same answers."""

    LOAD = f'''
      src <- parse("{HELPERS_R}", keep.source = FALSE)
      delimp <- new.env()
      for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-")) &&
                         identical(as.character(e[[2]]), "compute_onoff_proteins")) eval(e, delimp)
      stopifnot(is.function(delimp$compute_onoff_proteins))
      skill <- new.env(); sys.source("{ONOFF_R}", skill)
    '''

    def test_same_code(self):
        out = rscript(self.LOAD + '''
          a <- deparse(delimp$compute_onoff_proteins); b <- deparse(skill$compute_onoff_proteins)
          cat(identical(a, b))''')
        self.assertEqual(out.strip(), "TRUE",
                         "onoff.R's compute_onoff_proteins() no longer matches R/helpers.R: copy "
                         "the DE-LIMP function over verbatim (or fix both) -- one definition")

    def test_same_answers_on_every_scenario(self):
        calls = []
        for k, args in SCENARIOS.items():
            a = ("contrasts_list = NULL, " + args[len("NULL, "):]) if args.startswith("NULL") \
                else "contrasts_list = " + args
            calls.append(f'"{k}" = function(f) f(E, grp, {a})')
        out = rscript(self.LOAD + FIXTURE_R + f'''
          sc <- list({", ".join(calls)})
          bad <- character(0)
          for (k in names(sc)) {{
            a <- sc[[k]](delimp$compute_onoff_proteins); b <- sc[[k]](skill$compute_onoff_proteins)
            if (!identical(a, b)) bad <- c(bad, k)
          }}
          cat(length(sc), "|", paste(bad, collapse = ","))''')
        n, bad = [x.strip() for x in out.split("|")]
        self.assertEqual(int(n), len(SCENARIOS))
        self.assertEqual(bad, "", f"DE-LIMP and the skill's copy differ on: {bad}")


if __name__ == "__main__":
    unittest.main()
