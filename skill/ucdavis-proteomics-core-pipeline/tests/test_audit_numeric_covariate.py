#!/usr/bin/env python3
"""AUDIT.md says when a covariate fitted as a number nearly matches the groups (2.12).

run_de.R --numeric-covariates records each covariate in de_provenance.json `covariates`: for a
numeric one, how much of its variance Group explains (`group_r2`), whether that reaches the
caution (`confounded`, from `confound_r2_threshold`), how much it widens each contrast's standard
error (`se_inflation`), and its caution in `warnings`. The Methods, the brief and the report
already said so; the audit -- the list of what to surface to the user -- did not. Now the
numeric_covariate check quotes run_de.R's caution as a WARN (composed from the term only when a
record has no caution text), passes a numeric covariate that is not confounded, and says nothing
when there is none (a factor covariate is the `confounding` check's). Stdlib only.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")

CAUTION = ("Covariate1 is strongly confounded with Group: Group explains 93% of its variance "
           "(R^2 0.93, warned from 0.8). The group comparisons are adjusted for a covariate that "
           "nearly coincides with the groups, so they rest largely on its fitted straight-line "
           "slope, and their standard errors are up to 3.8x those without it (B-A 3.80x). If "
           "Covariate1 differs between the groups because of the grouping itself (a consequence "
           "of the treatment), adjusting for it removes part of the effect being tested: ask "
           "whether it belongs in the model.")
NUM = {"column": "Covariate1", "type": "numeric", "source": "--numeric-covariates",
       "n_distinct": 12, "min": -1.1, "max": 1.4, "group_r2": 0.21, "confounded": False,
       "se_inflation": {"B-A": 1.12, "C-A": 1.05}}


def audit(covariates):
    with tempfile.TemporaryDirectory() as d:
        de = os.path.join(d, "tables")
        os.makedirs(de)
        prov = {"method": "maxlfq", "adjp": 0.05, "logfc": 1.0}
        if covariates is not None:
            prov["covariates"] = covariates
        with open(os.path.join(de, "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        out = os.path.join(d, "AUDIT.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"), "--out", out,
                            "--de-dir", de], capture_output=True, text=True, cwd=d, timeout=120)
        assert r.returncode == 0, r.stderr
        with open(os.path.join(d, "AUDIT.json")) as fh:
            res = json.load(fh)
        with open(out, encoding="utf-8") as fh:
            md = fh.read()
    return [f for f in res["findings"] if f["check"] == "numeric_covariate"], res, md


class NumericCovariateAudit(unittest.TestCase):
    def test_a_confounded_one_is_a_warn_quoting_run_de(self):
        t = dict(NUM, group_r2=0.93, confounded=True, se_inflation={"B-A": 3.8})
        (f,), res, md = audit({"terms": [t], "numeric_requested": ["Covariate1"],
                               "confound_r2_threshold": 0.8, "warnings": [CAUTION]})
        self.assertEqual(f["status"], "WARN")
        self.assertEqual(f["message"], CAUTION + " Run the DE without Covariate1 as well, and "
                         "compare, before trusting the adjusted comparisons.")
        self.assertEqual(f["detail"]["group_r2"], 0.93)
        self.assertEqual(res["overall"], "WARN")
        self.assertIn("removes part of the effect being tested", md)

    def test_without_run_de_s_words_the_caution_is_composed_from_the_term(self):
        t = dict(NUM, group_r2=0.88, confounded=True, se_inflation={"B-A": 2.5, "C-A": 1.9})
        (f,), _, _ = audit({"terms": [t], "warnings": []})
        self.assertEqual(f["status"], "WARN")
        self.assertIn("Group explains 88% of its variance", f["message"])
        self.assertIn("standard errors are up to 2.5x those without it (B-A 2.50x, C-A 1.90x)",
                      f["message"])
        self.assertIn("removes part of the effect being tested", f["message"])

    def test_one_that_is_not_confounded_passes(self):
        (f,), res, _ = audit({"terms": [NUM], "confound_r2_threshold": 0.8, "warnings": []})
        self.assertEqual(f["status"], "PASS")
        self.assertEqual(f["message"], "Covariate1 was fitted as a number; Group explains 21% of "
                                       "its variance (cautioned from 80%), and it widens the "
                                       "comparisons' standard errors by at most 1.12x.")
        self.assertEqual(res["overall"], "PASS")

    def test_each_numeric_covariate_is_judged_on_its_own(self):
        conf = dict(NUM, column="Covariate2", group_r2=0.95, confounded=True,
                    se_inflation="not computed: singular")
        found, _, _ = audit({"terms": [NUM, conf], "warnings": []})
        self.assertEqual([(f["detail"]["column"], f["status"]) for f in found],
                         [("Covariate1", "PASS"), ("Covariate2", "WARN")])
        self.assertNotIn("standard errors are up to", found[1]["message"])

    def test_nothing_without_a_numeric_covariate(self):
        for cov in (None, {"terms": []}, {"terms": [{"column": "Batch", "type": "factor",
                                                     "n_levels": 3}], "warnings": []},
                    {"terms": [], "absorbed_by_block": ["Covariate1"]}):
            with self.subTest(cov=cov):
                found, _, _ = audit(cov)
                self.assertEqual(found, [])


if __name__ == "__main__":
    unittest.main()
