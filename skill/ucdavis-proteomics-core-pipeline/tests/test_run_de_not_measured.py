#!/usr/bin/env python3
"""A protein measured in neither group of a contrast is not tested in that contrast (2.12.0).

Staff report, 2026-10-03 (skill 2.9.1): in a several-group DE, proteins measured in 0 of n runs
of BOTH groups of a contrast -- measured only in the other groups -- were labelled `presence
call` (a label that means present in one group, absent in the other), tested, and counted in that
contrast's Benjamini-Hochberg adjustment. On dpc both sides are the detection model's inferred
value, so logFC = 0 and p = 1 exactly: 186 of the 249 proteins with p > 0.975 in one contrast
were 0/4 vs 0/4, every p-value histogram had a spike at 1, and the BH correction counted
hypotheses that carried no information.

Guards, for dpc and maxlfq, with and without technical replicates (the Sample column,
duplicateCorrelation):
  * such a protein's Evidence is `not measured in either group`; its logFC, t, P.Value,
    adj.P.Val and B are NA (the row stays); `presence call` still means measured in one group;
  * adj.P.Val is BH over the tested proteins only, and the significant counts are those;
  * de_provenance.json records the rule and, per contrast, how many were not tested; methods.txt
    and make_methods.py say so; reproducibility_log.R reproduces the same NA rows;
  * the consumers (analysis page counts, normalisation check volcano, audit, set tests' measured-
    only version) read the NA rows as not tested.

End to end on a synthetic DIA-NN report (needs R with limpa/limma/statmod/arrow/...; skips
otherwise). NOTE: CI runs Python only and skips every R-invoking test -- run them locally.
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

import analysis_prompt as ap  # noqa: E402
import audit_results as ar  # noqa: E402
import make_analysis_html as mah  # noqa: E402
import make_methods as mm  # noqa: E402
import normalization_check as nc  # noqa: E402
from test_run_de_contaminants import r_has, rscript  # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")
LABEL = "not measured in either group"
NA_COLS = ("logFC", "t", "P.Value", "adj.P.Val", "B")

# 4 groups x 3 runs. P00001-P00050 in every run (P00001-P00005 up 1.5 log2 in B); P00051-P00070
# only in C and D (6 of 12 runs: they pass maxlfq's 50% coverage filter), so in B vs A they are
# 0/3 vs 0/3 -- and presence calls in C vs A; P00071-P00075 only in B (dpc: presence calls in
# B vs A; maxlfq's coverage filter drops them). The Sample column makes two injections of one
# sample in every group: technical replicates.
SYNTH_R = r'''
set.seed(7)
runs <- sprintf("run%02d", 1:12); grp <- rep(c("A", "B", "C", "D"), each = 3)
smp <- paste0(tolower(grp), c(1, 1, 2))
rows <- list()
add <- function(pg, prec, base, only = NULL, effB = 0) {
  for (r in seq_along(runs)) {
    if (!is.null(only) && !(grp[r] %in% only)) next
    lv <- base + rnorm(1, 0, 0.2) + (grp[r] == "B") * effB
    int <- 2^(lv + rnorm(1, 0, 0.3))
    if (runif(1) < plogis(-(lv - 13) * 2)) next
    rows[[length(rows) + 1]] <<- data.frame(Run = runs[r], Precursor.Id = prec,
      Protein.Group = pg, Protein.Ids = pg, Protein.Names = paste0(pg, "_X"), Genes = sub("P0*", "G", pg),
      Proteotypic = 1L, Precursor.Normalised = int, Precursor.Quantity = int,
      Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001,
      Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
  }
}
for (i in 1:75) {
  pg <- sprintf("P%05d", i)
  base <- if (i > 50) runif(1, 17, 20) else runif(1, 14, 20)
  only <- if (i > 70) "B" else if (i > 50) c("C", "D") else NULL
  for (k in 1:4) add(pg, sprintf("PEP%dK%d2", i, k), base + rnorm(1, 0, 0.5), only = only,
                     effB = if (i <= 5) 1.5 else 0)
}
d <- do.call(rbind, rows)
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum)
names(pgm)[3] <- "PG.MaxLFQ"
d <- merge(d, pgm, by = c("Protein.Group", "Run"))
arrow::write_parquet(d, file.path(OUT, "report.parquet"))
write.csv(data.frame(File.Name = runs, Group = grp), file.path(OUT, "conditions.csv"), row.names = FALSE)
write.csv(data.frame(File.Name = runs, Group = grp, Sample = smp), file.path(OUT, "conditions_tr.csv"),
          row.names = FALSE)
'''
CD_ONLY = {f"P{i:05d}" for i in range(51, 71)}
B_ONLY = {f"P{i:05d}" for i in range(71, 76)}      # dpc only: maxlfq's coverage filter drops them


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if v != v else v


def bh(ps):
    """Benjamini-Hochberg, as p.adjust(method = "BH")."""
    n = len(ps)
    order = sorted(range(n), key=lambda i: ps[i], reverse=True)
    out, run = [0.0] * n, 1.0
    for rank, i in enumerate(order):
        run = min(run, ps[i] * n / (n - rank))
        out[i] = run
    return out


def zeros(row):
    """True when every Detected_<group> of the row is 0/n."""
    det = [v for k, v in row.items() if k.startswith("Detected_")]
    return bool(det) and all(v.split("/")[0] == "0" for v in det)


@unittest.skipUnless(r_has("limpa", "limma", "statmod", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limpa/limma/statmod/arrow/dplyr/tidyr/jsonlite")
class NotMeasuredInEitherGroup(unittest.TestCase):
    RUNS = (("dpc", "conditions.csv", "dpc"), ("ml", "conditions.csv", "maxlfq"),
            ("dpc_tr", "conditions_tr.csv", "dpc"), ("ml_tr", "conditions_tr.csv", "maxlfq"))

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        cls.runs = {}
        for key, meta, method in cls.RUNS:
            out = os.path.join(cls.tmp, "out_" + key)
            p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(cls.tmp, "report.parquet"),
                                "--metadata", os.path.join(cls.tmp, meta), "--method", method,
                                "--outdir", out], capture_output=True, text=True, cwd=cls.tmp)
            cls.runs[key] = (p, out, method)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, key):
        p, out, _ = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-3000:])
        return out

    def table(self, key, contrast):
        _, _, method = self.runs[key]
        return read_csv(os.path.join(self.out(key), f"DE_{method}_{contrast}.csv"))

    def prov(self, key):
        with open(os.path.join(self.out(key), "de_provenance.json")) as fh:
            return json.load(fh)

    def test_measured_in_neither_group_is_labelled_and_not_tested(self):
        for key, _, _ in self.RUNS:
            with self.subTest(run=key):
                rows = {r["Protein.Group"]: r for r in self.table(key, "B.A")}
                self.assertTrue(CD_ONLY <= set(rows), sorted(CD_ONLY - set(rows)))
                for pg in CD_ONLY:
                    r = rows[pg]
                    self.assertEqual((r["Detected_B"], r["Detected_A"]), ("0/3", "0/3"), pg)
                    self.assertEqual(r["Evidence"], LABEL, pg)
                    for c in NA_COLS:
                        self.assertEqual(r[c], "NA", f"{pg} {c}")
                # the label is exactly "every compared group 0/n"
                for pg, r in rows.items():
                    self.assertEqual(r["Evidence"] == LABEL, zeros(r), pg)

    def test_the_delimp_session_leaves_them_untested(self):
        """2.12 review (M4): DE-LIMP draws its volcano and tables with topTable(session$fit), which
        gave an untested protein its p-value back and counted it in BH. In the session's fit its t,
        p-value and B are NA for that contrast, so topTable reproduces DE_*.csv."""
        for key in ("dpc", "dpc_tr"):
            out = self.out(key)
            res = json.loads(rscript(f'''
              suppressMessages(library(limma)); s <- readRDS("{out}/DE-LIMP_session.rds")
              r <- lapply(setNames(colnames(s$fit$p.value), colnames(s$fit$p.value)), function(cn) {{
                tt <- topTable(s$fit, coef = cn, number = Inf)
                de <- read.csv(file.path("{out}", sprintf("DE_dpc_%s.csv", make.names(cn))))
                de <- de[match(rownames(tt), de$Protein.Group), ]
                list(na_p = sum(is.na(tt$P.Value)), na_same = identical(is.na(tt$adj.P.Val), is.na(de$adj.P.Val)),
                     dp = max(abs(tt$P.Value - de$P.Value), na.rm = TRUE),
                     dadj = max(abs(tt$adj.P.Val - de$adj.P.Val), na.rm = TRUE),
                     na_t = sum(is.na(s$fit$t[, cn])), na_lods = sum(is.na(s$fit$lods[, cn])),
                     n_listed = length(s$not_measured$proteins[[cn]]))
              }})
              cat(jsonlite::toJSON(list(per = r, na_F = sum(is.na(s$fit$F)), label = s$not_measured$label),
                                   auto_unbox = TRUE, digits = NA))'''))
            self.assertEqual(res["label"], LABEL)
            n_nm = {c: t["n_not_measured"] for c, t in self.prov(key)["de_tables"].items()}
            self.assertEqual(n_nm["B-A"], 20, key)
            for cn, r in res["per"].items():
                self.assertEqual((r["na_p"], r["na_t"], r["na_lods"], r["n_listed"]),
                                 (n_nm[cn],) * 4, (key, cn))
                self.assertTrue(r["na_same"], (key, cn))
                self.assertLess(max(r["dp"], r["dadj"]), 1e-12, (key, cn))
            # F (all contrasts together): NA for a protein any contrast does not test
            self.assertEqual(res["na_F"], len(CD_ONLY | B_ONLY), key)

    def test_presence_call_still_means_measured_in_one_group(self):
        for key, _, _ in self.RUNS:
            with self.subTest(run=key):
                rows = {r["Protein.Group"]: r for r in self.table(key, "C.A")}
                for pg in CD_ONLY:
                    self.assertEqual(rows[pg]["Evidence"], "presence call", pg)
                    self.assertEqual(rows[pg]["Detected_A"], "0/3", pg)
                    self.assertNotEqual(rows[pg]["Detected_C"].split("/")[0], "0", pg)
                    if key.startswith("dpc"):           # dpc tests them; maxlfq cannot (A all NA)
                        self.assertIsNotNone(num(rows[pg]["P.Value"]), pg)
                # the B-only proteins are measured in neither C nor A
                want = B_ONLY if key.startswith("dpc") else set()
                self.assertEqual({pg for pg, r in rows.items() if r["Evidence"] == LABEL}, want)
                if key.startswith("dpc"):
                    ba = {r["Protein.Group"]: r for r in self.table(key, "B.A")}
                    for pg in B_ONLY:
                        self.assertEqual(ba[pg]["Evidence"], "presence call", pg)
                        self.assertIsNotNone(num(ba[pg]["P.Value"]), pg)

    def test_bh_is_over_the_tested_proteins_only(self):
        for key, _, _ in self.RUNS:
            for con in ("B.A", "C.A", "D.A"):
                with self.subTest(run=key, contrast=con):
                    rows = self.table(key, con)
                    tested = [r for r in rows if num(r["P.Value"]) is not None]
                    want = bh([num(r["P.Value"]) for r in tested])
                    for r, w in zip(tested, want):
                        self.assertAlmostEqual(num(r["adj.P.Val"]), w, places=10, msg=r["Protein.Group"])
                    for r in rows:
                        if num(r["P.Value"]) is None:
                            self.assertIsNone(num(r["adj.P.Val"]), r["Protein.Group"])

    def test_untested_proteins_would_have_inflated_bh(self):
        # before 2.12.0 the 20 sat in B vs A's BH at p = 1 (dpc): every adjusted p-value of a
        # tested protein is now no larger than it was with them in
        rows = self.table("dpc", "B.A")
        ps = [num(r["P.Value"]) for r in rows if num(r["P.Value"]) is not None]
        before = bh(ps + [1.0] * len(CD_ONLY))[:len(ps)]
        after = [num(r["adj.P.Val"]) for r in rows if num(r["P.Value"]) is not None]
        self.assertTrue(all(a <= b + 1e-12 for a, b in zip(after, before)))
        self.assertLess(sum(after), sum(before))

    def test_the_record_says_what_was_not_tested(self):
        for key, _, method in self.RUNS:
            with self.subTest(run=key):
                prov = self.prov(key)
                rec = prov["not_measured"]
                self.assertEqual(rec["label"], LABEL)
                self.assertIn("Benjamini-Hochberg", rec["rule"])
                self.assertEqual(rec["columns_na"], list(NA_COLS))
                self.assertIn(f"'{LABEL}'", prov["detection_matrix"]["de_columns"]["Evidence"])
                for con, t in prov["de_tables"].items():
                    rows = read_csv(os.path.join(self.out(key), t["file"]))
                    self.assertEqual(t["n_proteins"], len(rows), con)
                    self.assertEqual(t["n_tested"],
                                     sum(num(r["P.Value"]) is not None for r in rows), con)
                    self.assertEqual(t["n_not_measured"],
                                     sum(r["Evidence"] == LABEL for r in rows), con)
                    sig = sum(num(r["adj.P.Val"]) is not None and num(r["adj.P.Val"]) < 0.05
                              for r in rows)
                    self.assertEqual(t["n_significant"], sig, con)
                    self.assertEqual(prov["significant_per_contrast"][con], sig, con)
                self.assertEqual(prov["de_tables"]["B-A"]["n_not_measured"], len(CD_ONLY))
                self.assertEqual(prov["de_tables"]["C-A"]["n_not_measured"],
                                 len(B_ONLY) if method == "dpc" else 0)

    def test_methods_txt_and_make_methods_say_so(self):
        for key, _, _ in self.RUNS:
            with self.subTest(run=key):
                with open(os.path.join(self.out(key), "methods.txt")) as fh:
                    txt = fh.read()
                self.assertIn("Not tested", txt)
                self.assertIn(f"B-A {len(CD_ONLY)}", txt)
                para = mm.de_paragraph(self.prov(key))
                self.assertIn("was not tested in that contrast", para)
                self.assertIn(f"B-A: {len(CD_ONLY)} of ", para)
                self.assertIn(f"C-A: {len(B_ONLY) if key.startswith('dpc') else 0} of ", para)

    def test_the_run_log_counts_tested_proteins(self):
        p, _, _ = self.runs["dpc"]
        self.assertRegex(p.stderr, r"B-A +75 proteins, 55 tested \(20 not measured in either "
                                   r"group\), \d+ significant")

    def test_reproducibility_log_reproduces_the_untested_rows(self):
        for key in ("dpc", "ml", "dpc_tr", "ml_tr"):
            with self.subTest(run=key):
                out = self.out(key)
                p = subprocess.run(["Rscript", "reproducibility_log.R"], capture_output=True,
                                   text=True, cwd=out)
                self.assertEqual(p.returncode, 0, p.stderr[-2000:])
                for t in self.prov(key)["de_tables"].values():
                    a = {r["Protein.Group"]: r for r in read_csv(os.path.join(out, t["file"]))}
                    b = {r["Protein.Group"]: r for r in
                         read_csv(os.path.join(out, "de_results_rerun", t["file"]))}
                    self.assertEqual(set(a), set(b))
                    for k in a:
                        for col in ("logFC", "P.Value", "adj.P.Val"):
                            x, y = num(a[k][col]), num(b[k][col])
                            if x is None:
                                self.assertIsNone(y, f"{t['file']} {k} {col}")
                            else:
                                self.assertAlmostEqual(x, y, places=8, msg=f"{t['file']} {k} {col}")

    def test_figures_leave_untested_proteins_out(self):
        if not r_has("ggplot2"):
            self.skipTest("needs ggplot2")
        out = self.out("dpc")
        figs = os.path.join(self.tmp, "figs_dpc")
        p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "make_figures.R"), "--de-dir", out,
                            "--outdir", figs], capture_output=True, text=True, cwd=self.tmp)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(figs, "figures.json")) as fh:
            fj = json.load(fh)
        cap = json.dumps(fj)
        self.assertIn(f"{len(CD_ONLY)} proteins measured in neither group of B vs A", cap)


# ---- the consumers, without R ----------------------------------------------------------------
NOT_MEASURED_REC = {"label": LABEL, "columns_na": list(NA_COLS),
                    "rule": "a protein measured in no run of the groups a contrast compares is not "
                            "tested in it: its logFC, t, P.Value, adj.P.Val and B are NA and it is "
                            "left out of that contrast's Benjamini-Hochberg adjustment"}
PROV = {"display_label": "DPC-Quant + limma (limpa)", "pipeline_id": "dpc", "method": "dpc",
        "adjp": 0.05, "logfc": 1.0, "logfc_role": "reference_line_only",
        "de_engine": "limpa::dpcDE", "contrasts": ["B-A", "C-A"],
        "not_measured": NOT_MEASURED_REC,
        "de_tables": {"B-A": {"file": "DE_dpc_B.A.csv", "n_proteins": 75, "n_tested": 55,
                              "n_not_measured": 20, "n_significant": 3},
                      "C-A": {"file": "DE_dpc_C.A.csv", "n_proteins": 75, "n_tested": 75,
                              "n_not_measured": 0, "n_significant": 9}}}


class MethodsSentence(unittest.TestCase):
    def test_states_the_rule_and_the_counts(self):
        s = mm.de_paragraph(PROV)
        self.assertIn("was not tested in that contrast", s)
        self.assertIn("left out of that contrast's Benjamini–Hochberg adjustment", s)
        self.assertIn("B-A: 20 of 75; C-A: 0 of 75", s)

    def test_an_older_record_says_nothing_it_did_not_do(self):
        old = {k: v for k, v in PROV.items() if k != "not_measured"}
        self.assertNotIn("not tested", mm.de_paragraph(old))


def write_de(root, rows, name="DE_dpc_B.A.csv"):
    os.makedirs(root, exist_ok=True)
    head = ["Protein.Group", "Genes", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B",
            "Detected_B", "Detected_A", "Evidence"]
    with open(os.path.join(root, name), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(head)
        for r in rows:
            w.writerow([r.get(h, "") for h in head])
    return root


def synthetic_rows(n_tested=60, n_untested=20, n_sig=4):
    rows = []
    for i in range(n_tested):
        p = 1e-6 * (i + 1) if i < n_sig else (i + 1) / (n_tested + 1)
        rows.append({"Protein.Group": f"P{i:05d}", "Genes": f"G{i}", "logFC": 1.5 if i < n_sig else 0.1,
                     "AveExpr": 15, "t": 3, "P.Value": p, "adj.P.Val": min(1.0, p * n_tested / (i + 1)),
                     "B": 0, "Detected_B": "3/3", "Detected_A": "3/3", "Evidence": "measured in both"})
    for i in range(n_untested):
        rows.append({"Protein.Group": f"Q{i:05d}", "Genes": f"H{i}", "logFC": "NA", "AveExpr": 18,
                     "t": "NA", "P.Value": "NA", "adj.P.Val": "NA", "B": "NA",
                     "Detected_B": "0/3", "Detected_A": "0/3", "Evidence": LABEL})
    return rows


class ConsumersReadUntestedRowsAsUntested(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_analysis_page_counts_tested_proteins(self):
        d = write_de(self.tmp, synthetic_rows())
        prov = dict(PROV, contrasts=["B-A"])
        t = mah.Tables(d, prov, 0.05, "de_provenance.json")
        self.assertEqual(t.counts("B.A"), (60, 4, 0))
        g = mah.glance_data(prov, t, d)
        self.assertEqual(g["contrasts"][0]["tested"], 60)
        self.assertTrue(any(LABEL in json.dumps(n) and "B vs A 20" in json.dumps(n)
                            for n in g["notes"]), g["notes"])

    def test_normalisation_check_volcano_has_no_spike_at_one(self):
        rows = [{k: str(v) for k, v in r.items()} for r in synthetic_rows(n_tested=60, n_untested=40)]
        v = nc.volcano(rows, 0.05, False, None)
        self.assertEqual(v["n"], 60)
        self.assertFalse(v["pvalue"]["shape"].startswith("odd"), v["pvalue"])

    def test_audit_counts_tested_proteins(self):
        write_de(self.tmp, synthetic_rows())
        findings = []
        ar.audit_de(findings, self.tmp, 0.05, 1.0)
        self.assertEqual(len(findings), 1, findings)
        self.assertIn("4/60 proteins significant", json.dumps(findings))

    def test_the_brief_explains_the_label(self):
        self.assertIn(LABEL, ap.EVIDENCE)


class BriefAndPage(unittest.TestCase):
    """The analysis brief and the report page, end to end on test_report_280's session."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name

    def tearDown(self):
        self._td.cleanup()

    def prov(self):
        from test_report_280 import PROV as P280
        return dict(P280, not_measured=NOT_MEASURED_REC,
                    de_tables={"Bait-IgG": {"file": "DE_dpc_Bait.IgG.csv", "n_proteins": 10,
                                            "n_tested": 7, "n_not_measured": 3,
                                            "n_significant": 2}})

    def test_the_brief_says_untested_rows_are_never_hits(self):
        from test_report_280 import brief
        b = brief(self.root, prov=self.prov(), det_cols=True, evidence_col=True)
        self.assertIn(f"**{LABEL}**: measured only in other groups of the study, so NOT TESTED", b)
        self.assertIn("count tested proteins per comparison as the rows with a P.Value", b)

    def test_the_page_says_how_many_were_not_tested(self):
        from test_report_280 import page
        _, html, md, _ = page(self.root, prov=self.prov(), det_cols=True, evidence_col=True)
        for text in (html, md):
            self.assertIn("Proteins measured only in other groups are not tested in a comparison", text)
            self.assertIn("Bait vs IgG 3", text)


class OneVocabulary(unittest.TestCase):
    """run_de.R's EVIDENCE_LABELS is the definition; set_tests.R and analysis_prompt.py mirror it."""

    def labels(self):
        with open(RUN_DE) as fh:
            src = fh.read()
        m = re.search(r"EVIDENCE_LABELS <- c\((.*?)\)\n", src, re.S)
        self.assertIsNotNone(m)
        return dict(re.findall(r'(\w+) = "([^"]+)"', m.group(1)))

    def test_set_tests_mirror_run_de(self):
        lab = self.labels()
        with open(os.path.join(SCRIPTS, "set_tests.R")) as fh:
            st = fh.read()
        self.assertIn(f'PRESENCE_CALL_LABEL <- "{lab["none_in_a_group"]}"', st)
        self.assertIn(f'NOT_MEASURED_LABEL <- "{lab["none_measured"]}"', st)

    def test_analysis_prompt_mirrors_run_de(self):
        lab = self.labels()
        self.assertEqual(ap.EVIDENCE[:2], (lab["all"], lab["none_in_a_group"]))
        self.assertIn(lab["none_measured"], ap.EVIDENCE)


# The set tests' primary universe (2.12 review, L1): the GMT's sets mix proteins tested in B vs A
# with C/D-only ones (not measured in either group of B-A).
SETS_GMT = "\n".join("\t".join([name, "synthetic"] + genes) for name, genes in (
    ("b_up", [f"G{i}" for i in (1, 2, 3, 4, 5, 11, 12, 13, 14, 15)]),
    ("cd_mix", [f"G{i}" for i in range(51, 61)] + [f"G{i}" for i in range(6, 11)]),
    ("cd_only", [f"G{i}" for i in range(61, 71)]),
    ("bg_a", [f"G{i}" for i in range(16, 31)]),
    ("bg_b", [f"G{i}" for i in range(31, 46)]))) + "\n"


@unittest.skipUnless(r_has("limpa", "limma", "statmod", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limpa/limma/statmod/arrow/dplyr/tidyr/jsonlite")
class SetTestsAndLostRuns(unittest.TestCase):
    """L1: run_sets.R's primary tests leave out a contrast's untested proteins, as its DE table
    does (they had t ~ 0 on dpc and joined camera's background). L2: a design run that has no
    quantities (no precursor passed the identification filters) is recorded as not quantified,
    with that reason -- not as a left-out run whose reason is missing."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R + r'''
d <- arrow::read_parquet(file.path(OUT, "report.parquet"))
d$Q.Value[d$Run == "run12"] <- 0.5        # nothing of run12 passes the identification filter
arrow::write_parquet(d, file.path(OUT, "report_lost.parquet"))
''')
        with open(os.path.join(cls.tmp, "sets.gmt"), "w") as fh:
            fh.write(SETS_GMT)
        cls.de = os.path.join(cls.tmp, "de")
        run = lambda *a: subprocess.run(["Rscript", *a], capture_output=True, text=True, cwd=cls.tmp)
        cls.de_p = run(RUN_DE, "--input", "report.parquet", "--metadata", "conditions.csv",
                       "--method", "dpc", "--outdir", cls.de)
        cls.sets_p = run(os.path.join(SCRIPTS, "run_sets.R"), "--de-dir", cls.de, "--sets", "none",
                         "--gmt", "sets.gmt", "--gmt-origin", "synthetic sets written before run_de",
                         "--min-size", "5")
        with open(os.path.join(cls.tmp, "cd_only.gmt"), "w") as fh:
            fh.write("cd_only\tsynthetic\t" + "\t".join(f"G{i}" for i in range(61, 71)) + "\n")
        cls.cd_p = run(os.path.join(SCRIPTS, "run_sets.R"), "--de-dir", cls.de, "--sets", "none",
                       "--gmt", "cd_only.gmt", "--gmt-origin", "synthetic sets written before run_de",
                       "--min-size", "5", "--contrasts", "B-A", "--outdir", os.path.join(cls.tmp, "cd"))
        cls.lost = {m: (run(RUN_DE, "--input", "report_lost.parquet", "--metadata", "conditions.csv",
                            "--method", m, "--outdir", os.path.join(cls.tmp, "lost_" + m)),
                        os.path.join(cls.tmp, "lost_" + m)) for m in ("dpc", "maxlfq")}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def sets(self, contrast):
        self.assertEqual(self.de_p.returncode, 0, self.de_p.stderr[-2000:])
        self.assertEqual(self.sets_p.returncode, 0, self.sets_p.stderr[-2000:])
        return {r["Set_ID"]: r for r in read_csv(os.path.join(self.de, f"Sets_{contrast}.csv"))}

    def test_set_tests_leave_out_what_the_de_table_does_not_test(self):
        ba, ca = self.sets("B.A"), self.sets("C.A")
        self.assertEqual(ba["cd_mix"]["NGenes"], "5")               # its 10 C/D-only proteins out
        self.assertNotIn("cd_only", ba)                             # nothing of it tested in B-A
        self.assertEqual((ca["cd_mix"]["NGenes"], ca["cd_only"]["NGenes"]), ("15", "10"))
        with open(os.path.join(self.de, "sets_provenance.json")) as fh:
            rec = json.load(fh)
        c = rec["contrasts"]
        self.assertEqual((c["B-A"]["n_not_measured_left_out"], c["B-A"]["n_sets_not_tested"]), (20, 1))
        self.assertEqual(c["C-A"]["n_not_measured_left_out"], 5)    # B-only: untested in C vs A
        self.assertEqual(c["B-A"]["n_proteins_tested"] + 20, c["C-A"]["n_proteins_tested"] + 5)
        self.assertIn("proteins measured in neither compared group -- not tested in its "
                      "differential-expression table -- were left out of both tests",
                      " ".join(rec["methods_paragraph"].split()))

    def test_a_comparison_that_tests_none_of_a_sets_proteins_stops_with_the_fix(self):
        self.assertNotEqual(self.cd_p.returncode, 0)
        self.assertIn("no set keeps 5 proteins tested in B-A", self.cd_p.stderr)
        self.assertIn("leave it out with --contrasts, or lower --min-size", self.cd_p.stderr)

    def test_a_design_run_with_no_quantities_is_not_a_left_out_run(self):
        for m, (p, out) in self.lost.items():
            with self.subTest(method=m):
                self.assertEqual(p.returncode, 0, p.stderr[-2000:])
                with open(os.path.join(out, "de_provenance.json")) as fh:
                    r = json.load(fh)["runs_left_out"]
                self.assertEqual(r["runs"], [])                 # nothing left out by the design
                (nq,) = r["not_quantified"]
                self.assertEqual(nq["run"], "run12")
                self.assertIn("in the report, but", nq["reason"])
                self.assertIn("no precursor of it passed" if m == "dpc" else "fewer than 2 of its proteins",
                              nq["reason"])
                self.assertEqual((r["n_report_runs"], r["n_analysed"]), (12, 11))
                with open(os.path.join(out, "methods.txt")) as fh:
                    methods = fh.read()
                self.assertIn("Not quantified: runs in the design with no quantities, not analysed:", methods)
                self.assertNotIn("NOT RECORDED [confirm]", methods)
                self.assertIn("1 run(s) of the design have no quantities", p.stderr)
                prov = {"runs_left_out": r}
                sent = mm.de_runs_left_out_sentence(prov)
                self.assertTrue(sent.startswith("1 run was in the design but had no quantities, so it "
                                                "was not analysed: run12 (not quantified: in the report"),
                                sent)
                self.assertEqual(mah.runs_left_out_note(prov)["kind"], "warning")


@unittest.skipUnless(r_has("limma", "jsonlite"), "needs R with limma and jsonlite")
class SetTestsMeasuredOnly(unittest.TestCase):
    def test_not_measured_counts_as_unmeasured(self):
        out = rscript(f'source("{os.path.join(SCRIPTS, "set_tests.R")}")\n' + r'''
ev <- c("measured in both", "presence call", "not measured in either group", "partly inferred", NA)
idx <- list(a = 1:3, b = c(1, 4, 5))
cat(jsonlite::toJSON(list(f = unname(presence_fraction(idx, ev)), m = measured_rows(ev)), digits = NA))
''')
        r = json.loads(out)
        self.assertAlmostEqual(r["f"][0], 2 / 3, places=6)
        self.assertEqual(r["f"][1], 0)
        self.assertEqual(r["m"], [1, 4, 5])


if __name__ == "__main__":
    unittest.main()
