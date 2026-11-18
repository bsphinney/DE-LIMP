#!/usr/bin/env python3
"""Unit tests (no limma/limpa) for the 2.15.x DE additions:

  * incomplete_data.R incomplete_flag_table(): tiers, NA rules, complex contrasts' group lists;
  * quantums_tradeoff.R: cutoff set, DEFAULT source text (rule 2), the as-run row;
  * onoff.R onoff_collapse_samples(): runs -> samples, with compute_onoff_proteins() untouched
    (the parity test in test_onoff_proteins.py keeps passing);
  * the Python consumers (methods sentences, analysis notes, MANIFEST lines) read the records only.
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

import make_analysis_html as mah  # noqa: E402
import make_deposit  # noqa: E402
import make_methods as mm  # noqa: E402
import session  # noqa: E402

HAVE_R = bool(shutil.which("Rscript"))


def rjson(code):
    p = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError(p.stderr[-2000:])
    return json.loads(p.stdout)


def has_jsonlite():
    return HAVE_R and subprocess.run(
        ["Rscript", "-e", "quit(status=!requireNamespace('jsonlite',quietly=TRUE))"]).returncode == 0


@unittest.skipUnless(has_jsonlite(), "needs Rscript with jsonlite")
class FlagTable(unittest.TestCase):
    def flags(self, k, n, grp, sig):
        k_r = "rbind(" + ",".join("c(" + ",".join("NA" if v is None else str(v) for v in row) + ")"
                                 for row in k) + ")"
        sig_r = "c(" + ",".join({True: "TRUE", False: "FALSE", None: "NA"}[s] for s in sig) + ")"
        grp_r = "c(" + ",".join(f'"{g}"' for g in grp) + ")"
        return rjson(f'source("{SCRIPTS}/incomplete_data.R"); '
                     f'x <- incomplete_flag_table({k_r}, c({",".join(map(str, n))}), {grp_r}, {sig_r}); '
                     f'cat(jsonlite::toJSON(x, na = "null", dataframe = "rows"))')

    def test_tiers_and_na_rules(self):
        r = self.flags(k=[[3, 3], [3, 2], [3, 1], [0, 3], [2, 1], [3, 3], [3, 2]], n=[3, 3],
                       grp=["Ctrl", "Treat"], sig=[True, True, True, True, True, False, None])
        self.assertEqual([x["values_per_group"] for x in r][:3],
                         ["Ctrl 3/3; Treat 3/3", "Ctrl 3/3; Treat 2/3", "Ctrl 3/3; Treat 1/3"])
        self.assertEqual([x["called_on_incomplete_data"] for x in r],
                         [False, True, True, True, True, None, None])
        self.assertEqual([x["incomplete_tier"] for x in r],
                         ["complete", "incomplete", "single_sample", "incomplete", "single_sample", None, None])

    def test_a_group_of_one_run_measured_once_is_complete_not_single(self):
        r = self.flags(k=[[1, 3]], n=[1, 3], grp=["A", "B"], sig=[True])
        self.assertEqual((r[0]["called_on_incomplete_data"], r[0]["incomplete_tier"]), (False, "complete"))

    def test_unknown_detection_gives_no_flag(self):
        r = self.flags(k=[[3, None]], n=[3, 3], grp=["A", "B"], sig=[True])
        self.assertIsNone(r[0].get("values_per_group"))
        self.assertIsNone(r[0].get("called_on_incomplete_data"))

    def test_complex_contrast_lists_all_its_groups(self):
        r = self.flags(k=[[3, 2, 3]], n=[3, 3, 3], grp=["C", "A", "B"], sig=[True])
        self.assertEqual(r[0]["values_per_group"], "C 3/3; A 2/3; B 3/3")
        self.assertTrue(r[0]["called_on_incomplete_data"])


@unittest.skipUnless(has_jsonlite(), "needs Rscript with jsonlite")
class QuantumsHelpers(unittest.TestCase):
    def test_cutoffs_source_and_as_run_row(self):
        r = rjson(f'source("{SCRIPTS}/quantums_tradeoff.R", encoding = "UTF-8"); '
                  'a <- quantums_cutoffs_for_run(0, 0); b <- quantums_cutoffs_for_run(0.6, 0.25); '
                  'cat(jsonlite::toJSON(list(a = vapply(a, function(c) c$label, ""), '
                  'b = vapply(b, function(c) c$label, ""), '
                  's0 = quantums_source_text(FALSE, FALSE), s1 = quantums_source_text(TRUE, FALSE), '
                  's2 = quantums_source_text(TRUE, TRUE), rec = QUANTUMS_RECOMMENDED$tag), auto_unbox = TRUE))')
        self.assertEqual(r["a"], ["none", "pgQ 0.5", "pgQ 0.75", "pgQ 0.75 + eQ 0.75", "pgQ 0.85"])
        self.assertEqual(r["b"][-1], "as run (pgQ 0.6, eQ 0.25)")
        self.assertEqual(len(r["b"]), 6)
        self.assertEqual(r["s0"], "DEFAULT — not user-confirmed")
        self.assertEqual(r["s1"], "pgQ: user-chosen; eQ: DEFAULT — not user-confirmed")
        self.assertEqual(r["s2"], "pgQ: user-chosen; eQ: user-chosen")
        self.assertEqual(r["rec"], "benchmark-recommended (PXD063416)")


@unittest.skipUnless(has_jsonlite(), "needs Rscript with jsonlite")
class CollapseSamples(unittest.TestCase):
    CODE = f'''source("{SCRIPTS}/onoff.R")
      E <- matrix(NA_real_, 3, 6, dimnames = list(c("p1","p2","p3"), paste0("r", 1:6)))
      E["p1", c("r1","r2")] <- 1          # both injections of sample s1 (group A)
      E["p2", "r3"] <- 1                  # one injection of sample s2 (group A)
      E["p3", c("r4","r5","r6")] <- 1     # group B, three samples
      grp <- factor(c("A","A","A","B","B","B"))
      smp <- c("s1","s1","s2","s3","s4","s5")
      cs <- onoff_collapse_samples(E, grp, smp)
      per_run <- compute_onoff_proteins(E, grp, list(c("A","B")), n_min = 2)
      per_smp <- compute_onoff_proteins(cs$E, cs$group, list(c("A","B")), n_min = 2)
      none <- onoff_collapse_samples(E, grp, NULL)
      uniq <- onoff_collapse_samples(E, grp, paste0("x", 1:6))
      cat(jsonlite::toJSON(list(coll = cs$collapsed, nr = cs$n_runs, ns = cs$n_samples,
            cols = colnames(cs$E), grp = cs$group, run = per_run$Protein.Group, smp = per_smp$Protein.Group,
            smp_tot = per_smp$total_g1, none = none$collapsed, uniq = uniq$collapsed,
            p1 = unname(cs$E["p1", ]), blank = onoff_collapse_samples(E, grp, c("s1","s1","","","",""))$n_samples),
            auto_unbox = FALSE, na = "null"))'''

    def test_runs_collapse_to_samples_any_run_counts(self):
        r = rjson(self.CODE)
        self.assertEqual((r["coll"], r["nr"], r["ns"]), ([True], [6], [5]))
        self.assertEqual(r["cols"], ["s1", "s2", "s3", "s4", "s5"])
        self.assertEqual(r["grp"], ["A", "A", "B", "B", "B"])
        # p1: two injections of ONE sample -> per run it is A_only (2 runs), per sample it is not (1 sample)
        self.assertIn("p1", r["run"])
        self.assertNotIn("p1", r["smp"])
        self.assertIn("p3", r["smp"])
        self.assertEqual(r["p1"][0], 1)
        self.assertEqual(r["p1"][1:], [None] * 4)

    def test_no_sample_info_or_nothing_to_collapse_leaves_runs_alone(self):
        r = rjson(self.CODE)
        self.assertEqual((r["none"], r["uniq"]), ([False], [False]))
        self.assertEqual(r["blank"], [5])        # blank Sample cells: each such run is its own sample


class TextFromRecords(unittest.TestCase):
    INC = {"columns": ["values_per_group", "called_on_incomplete_data", "incomplete_tier"],
           "rule": "RULE-WORDS", "benchmark_text": "BENCH-WORDS about 10x more likely to be wrong"}
    PROV = {"incomplete_data": INC,
            "de_tables": {"B-A": {"n_called_incomplete": 3, "n_significant": 20, "n_proteins": 100}},
            "quantums_tradeoff": {"status": "written", "file": "QuantUMS_tradeoff.csv",
                                  "used": {"pgQ": 0, "eQ": 0}, "cutoff_source": "DEFAULT — not user-confirmed",
                                  "recommended": {"label": "pgQ 0.75 + eQ 0.75", "tag": "benchmark-recommended (PXD063416)",
                                                  "citation": "CITE"},
                                  "defaults_unchanged": "defaults stay"}}

    def test_methods_sentences_use_the_records_words(self):
        s = mm.de_incomplete_sentence(self.PROV)
        self.assertIn("RULE-WORDS", s)
        self.assertIn("B-A: 3 of 20", s)
        self.assertIn("about 10x", s)
        q = mm.de_quantums_sentence(self.PROV)
        self.assertIn("DEFAULT — not user-confirmed", q)
        self.assertIn("benchmark-recommended (PXD063416)", q)
        self.assertIn("CITE", q)
        self.assertIsNone(mm.de_incomplete_sentence({}))
        self.assertIsNone(mm.de_quantums_sentence({}))
        skipped = {"quantums_tradeoff": {"status": "skipped", "file": "QuantUMS_tradeoff.csv", "reason": "WHY"}}
        self.assertIn("was not made (skipped): WHY", mm.de_quantums_sentence(skipped))

    def test_analysis_notes(self):
        n = mah.incomplete_note(self.PROV)
        self.assertEqual(n["kind"], "warning")
        self.assertIn("3 of 20", n["text"])
        self.assertIn("BENCH-WORDS", n["text"])
        none = {"incomplete_data": self.INC, "de_tables": {"B-A": {"n_called_incomplete": 0, "n_significant": 5}}}
        self.assertIsNone(mah.incomplete_note(none))          # nothing flagged: no callout
        self.assertIsNone(mah.incomplete_note({}))
        q = mah.quantums_note(self.PROV)
        self.assertIn("DEFAULT", q["text"])
        self.assertIn("your choice", q["text"])
        self.assertEqual(mah.quantums_note({"quantums_tradeoff": {"status": "failed", "reason": "boom"}})["kind"], "warning")
        self.assertIsNone(mah.quantums_note({"quantums_tradeoff": {"status": "skipped", "reason": "x"}}))


class QuantumsManifest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.p = {"session_dir": self.tmp.name, "de_dir": os.path.join(self.tmp.name, "output", "tables")}
        os.makedirs(self.p["de_dir"])

    def go(self, prov, file=False):
        with open(os.path.join(self.p["de_dir"], "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        if file:
            with open(os.path.join(self.p["de_dir"], "QuantUMS_tradeoff.csv"), "w") as fh:
                fh.write("x\n")
        man = make_deposit.Manifest()
        session.quantums_step(self.p, man)
        return man.lines

    def test_ok_skipped_failed_and_missing_file(self):
        inc = {"columns": []}
        tabs = {"B-A": {"n_called_incomplete": 2, "n_significant": 9}}
        ok = self.go({"incomplete_data": inc, "de_tables": tabs, "quantums_tradeoff": {
            "status": "written", "file": "QuantUMS_tradeoff.csv", "used": {"pgQ": 0, "eQ": 0},
            "cutoff_source": "DEFAULT"}}, file=True)
        self.assertTrue(ok[0].startswith("[OK]") and "B-A 2/9" in ok[0], ok)
        self.assertTrue(ok[1].startswith("[OK]") and "DEFAULT" in ok[1], ok)
        sk = self.go({"incomplete_data": inc, "de_tables": tabs, "quantums_tradeoff": {
            "status": "skipped", "reason": "limpa imputes"}})
        self.assertTrue(sk[1].startswith("[SKIPPED]") and "limpa imputes" in sk[1], sk)
        fail = self.go({"quantums_tradeoff": {"status": "failed", "reason": "boom"}})
        self.assertTrue(fail[0].startswith("[SKIPPED]") and "failed: boom" in fail[0], fail)
        os.remove(os.path.join(self.p["de_dir"], "QuantUMS_tradeoff.csv"))
        gone = self.go({"quantums_tradeoff": {"status": "written", "file": "QuantUMS_tradeoff.csv"}})
        self.assertTrue(gone[0].startswith("[SKIPPED]") and "is not in" in gone[0], gone)

    def test_old_record_without_the_flag_is_a_skipped_line_not_silence(self):
        ls = self.go({"de_tables": {"B-A": {}}})
        self.assertTrue(ls[0].startswith("[SKIPPED]") and "incomplete_data" in ls[0], ls)
        self.assertEqual(self.go({}), [])


if __name__ == "__main__":
    unittest.main()
