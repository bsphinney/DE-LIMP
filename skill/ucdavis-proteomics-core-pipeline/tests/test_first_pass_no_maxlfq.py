#!/usr/bin/env python3
"""
A DIA-NN --dda first-pass report has no engine protein quantities, and the skill says so
(Ruling 53, T23b).

DIA-NN 2.7.0's --dda first-pass report has PG.MaxLFQ = 0 in every row (two timsTOF DDA-PASEF
runs, every row of both; Precursor.Normalised > 0 everywhere). The default DE (--method dpc,
precursor-level) is unaffected, but `run_de.R --method maxlfq` set every value <= 0 to NA, dropped
every run and stopped with "MaxLFQ: fewer than 2 usable runs survived the filters." -- which sends
the user to the filters, not to the report.

Now ONE definition, dda_first_pass.pg_maxlfq_state() (is any PG.MaxLFQ > 0?), is:
  - recorded by the fallback: dda_first_pass.json and search_provenance.json dda_first_pass
    carry pg_maxlfq "absent" / "present";
  - asked by build_maxlfq.R (dda_first_pass.py pg-maxlfq), which refuses --method maxlfq on an
    "absent" report with a message naming the cause and --method dpc;
  - worded once in the caveat (dda_first_pass.caution): when the record says "absent", the
    Methods and AUDIT.md say protein quantities come from the precursors (limpa DPC-Quant).

Synthetic reports on a temp dir; the R parts skip without R. Hermetic.
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

import dda_first_pass  # noqa: E402
import audit_results as ar  # noqa: E402
import make_methods as mm  # noqa: E402
from test_normalization_check import HAVE_ARROW, NEEDS, make_report  # noqa: E402
from test_run_de_contaminants import r_has  # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")
CAUSE = "no engine protein quantities (PG.MaxLFQ is 0 or missing in every row)"


def write_parquet(path, cols):
    import pyarrow as pa
    import pyarrow.parquet as pq
    pq.write_table(pa.table(cols), path)
    return path


def zero_maxlfq(report):
    """`report` with PG.MaxLFQ set to 0 in every row: a DIA-NN --dda first pass."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    t = pq.read_table(report)
    i = t.schema.get_field_index("PG.MaxLFQ")
    t = t.set_column(i, "PG.MaxLFQ", pa.array([0.0] * t.num_rows))
    pq.write_table(t, report)
    return report


def fallback_prov(pg_maxlfq):
    rec = {"report_used": "first_pass", "exit_code": 132, "caveat": dda_first_pass.CAVEAT}
    if pg_maxlfq is not None:
        rec["pg_maxlfq"] = pg_maxlfq
    return {"engine": "diann", "version": "2.7.0", "dda": True, "report_used": "first_pass",
            "dda_first_pass": rec}


@unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
class TheOneDefinition(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def p(self, name):
        return os.path.join(self.d, name)

    def test_all_zero_is_absent(self):
        r = write_parquet(self.p("a.parquet"), {"Run": ["a", "b"], "PG.MaxLFQ": [0.0, 0.0]})
        self.assertEqual(dda_first_pass.pg_maxlfq_state(r), "absent")

    def test_zero_or_missing_in_every_row_is_absent(self):
        r = write_parquet(self.p("n.parquet"), {"Run": ["a", "b", "c"],
                                                "PG.MaxLFQ": [0.0, None, float("nan")]})
        self.assertEqual(dda_first_pass.pg_maxlfq_state(r), "absent")
        r = write_parquet(self.p("none.parquet"), {"Run": ["a", "b"]})
        self.assertEqual(dda_first_pass.pg_maxlfq_state(r), "absent", "no column at all")

    def test_one_positive_value_is_present(self):
        r = write_parquet(self.p("p.parquet"), {"Run": ["a", "b"], "PG.MaxLFQ": [0.0, 1250.0]})
        self.assertEqual(dda_first_pass.pg_maxlfq_state(r), "present")

    def test_a_tsv_report_is_read_too(self):
        path = self.p("report.tsv")
        for vals, want in ((["0", "0"], "absent"), (["0", "3.5e4"], "present")):
            with open(path, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["Run", "PG.MaxLFQ"])
                w.writerows([["a", vals[0]], ["b", vals[1]]])
            self.assertEqual(dda_first_pass.pg_maxlfq_state(path), want)

    def test_the_command_answers_in_json_for_run_de(self):
        r = write_parquet(self.p("a.parquet"), {"Run": ["a"], "PG.MaxLFQ": [0.0]})
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "dda_first_pass.py"),
                            "pg-maxlfq", "--report", r], capture_output=True, text=True,
                           timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        res = json.loads(p.stdout)
        self.assertEqual(res["pg_maxlfq"], "absent")
        self.assertIn(CAUSE, res["refusal"])
        self.assertIn("--method dpc", res["refusal"])


@unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
class TheRecord(unittest.TestCase):
    def run_record(self, maxlfq):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        report = os.path.join(d, "report.parquet")
        fp = write_parquet(dda_first_pass.first_pass_path(report),
                           {"Run": ["a", "b"], "PG.MaxLFQ": maxlfq})
        prov = os.path.join(d, "search_provenance.json")
        with open(prov, "w") as fh:
            json.dump({"engine": "diann"}, fh)
        rec = dda_first_pass.record(report, fp, prov, 132)
        with open(prov) as fh:
            sp = json.load(fh)
        with open(dda_first_pass.record_path(report)) as fh:
            side = json.load(fh)
        return rec, sp, side

    def test_an_all_zero_first_pass_is_recorded_absent(self):
        rec, sp, side = self.run_record([0.0, 0.0])
        self.assertEqual(rec["pg_maxlfq"], "absent")
        self.assertEqual(sp["dda_first_pass"]["pg_maxlfq"], "absent")
        self.assertEqual(side["pg_maxlfq"], "absent")

    def test_a_quantified_first_pass_is_recorded_present(self):
        rec, sp, side = self.run_record([0.0, 812.5])
        self.assertEqual(sp["dda_first_pass"]["pg_maxlfq"], "present")
        self.assertEqual(side["pg_maxlfq"], "present")


class TheCaveat(unittest.TestCase):
    def test_absent_says_protein_quantities_come_from_the_precursors(self):
        c = dda_first_pass.caution(fallback_prov("absent"))
        self.assertTrue(c.startswith(dda_first_pass.CAVEAT), c)
        self.assertIn(CAUSE.split(" (")[0], c)
        self.assertIn("precursors (limpa DPC-Quant)", c)
        for claim in ("fewer", "exit", "132"):
            self.assertNotIn(claim, c)

    def test_present_or_not_recorded_keeps_the_caveat_alone(self):
        """(pin) a record that says nothing about PG.MaxLFQ -- one from before T23b -- or says it
        is present gets the caveat as before."""
        for state in (None, "present"):
            self.assertEqual(dda_first_pass.caution(fallback_prov(state)), dda_first_pass.CAVEAT)

    def test_methods_and_audit_say_it(self):
        """The search paragraph states the fact only; DPC-Quant is the DE paragraph's, and only
        for a DE that ran dpc (C2 biologist B2 = statistician H2: a no-DE Methods said limpa had
        estimated protein quantities that nothing computed). AUDIT.md runs only after a DE."""
        with tempfile.TemporaryDirectory() as d:
            sp = os.path.join(d, "search_provenance.json")
            with open(sp, "w") as fh:
                json.dump(fallback_prov("absent"), fh)
            srec = mm.search_record(None, sp, None)
            para = mm.search_paragraph(srec)
            self.assertIn(dda_first_pass.NO_PG_MAXLFQ_FACT, para)
            self.assertNotIn("DPC-Quant", para)
            self.assertNotIn("limpa", para)
            de = mm.de_paragraph({"pipeline_id": "dpc", "display_label": "DPC-Quant + limma",
                                  "adjp": 0.05}, srec)
            self.assertIn(dda_first_pass.DE_DPC_SENTENCE, de)
            # another pipeline, or no search record: no such sentence
            self.assertNotIn(dda_first_pass.DE_DPC_SENTENCE,
                             mm.de_paragraph({"pipeline_id": "peptide_max", "adjp": 0.05}, srec))
            self.assertNotIn(dda_first_pass.DE_DPC_SENTENCE,
                             mm.de_paragraph({"pipeline_id": "dpc", "adjp": 0.05}))
            findings = []
            ar.audit_dda_first_pass(findings, {"input": os.path.join(d, "report.parquet")})
            self.assertIn("precursors (limpa DPC-Quant)", findings[0]["message"])
            self.assertEqual(findings[0]["detail"]["pg_maxlfq"], "absent")


class NoDEMethods(unittest.TestCase):
    def test_a_no_de_methods_md_never_names_dpc_quant(self):
        """(C2 statistician H2) two DDA runs, the --dda second pass crashed, no DE: methods.md
        says the report has no engine protein quantities, and never that limpa estimated any."""
        with tempfile.TemporaryDirectory() as d:
            sp = os.path.join(d, "search_provenance.json")
            with open(sp, "w") as fh:
                json.dump(fallback_prov("absent"), fh)
            raw = os.path.join(d, "run1.raw")
            open(raw, "w").close()
            md = os.path.join(d, "methods.md")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                                "--raw", raw, "--search-prov", sp, "--acquisition", "DDA",
                                "--out", md], capture_output=True, text=True, timeout=300,
                               env=dict(os.environ, THERMORAWFILEPARSER_SHARED="",
                                        THERMO_RESOLUTION_PYTHON=""))
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(md, encoding="utf-8") as fh:
                text = " ".join(fh.read().split())
            self.assertIn(dda_first_pass.NO_PG_MAXLFQ_FACT, text)
            self.assertNotIn("DPC-Quant", text)


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class RunDE(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def de(self, report, cond, contrast, method, name):
        out = os.path.join(self.d, name)
        return subprocess.run(["Rscript", RUN_DE, "--input", report, "--metadata", cond,
                               "--contrasts", contrast, "--method", method, "--outdir", out],
                              capture_output=True, text=True, timeout=900, cwd=self.d), out

    def test_maxlfq_is_refused_on_a_report_with_no_pg_maxlfq(self):
        report, cond, contrast = make_report(os.path.join(self.d, "fp"), "whole")
        zero_maxlfq(report)
        p, out = self.de(report, cond, contrast, "maxlfq", "maxlfq")
        self.assertEqual(p.returncode, 1, p.stdout[-1500:] + p.stderr[-1500:])
        err = " ".join(p.stderr.split())
        self.assertIn(CAUSE, err)
        self.assertIn("DIA-NN --dda first-pass report", err)
        self.assertIn("--method dpc", err)
        self.assertNotIn("fewer than 2 usable runs", err)
        self.assertFalse(os.path.exists(os.path.join(out, "de_provenance.json")))

    def test_dpc_runs_on_the_same_report(self):
        """(pin) the route the refusal names works on that report."""
        report, cond, contrast = make_report(os.path.join(self.d, "fp"), "whole")
        zero_maxlfq(report)
        p, out = self.de(report, cond, contrast, "dpc", "dpc")
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        self.assertTrue(os.path.isfile(os.path.join(out, "de_provenance.json")))

    def test_maxlfq_runs_on_a_normal_report(self):
        """(pin) a report with PG.MaxLFQ is not refused."""
        report, cond, contrast = make_report(os.path.join(self.d, "ok"), "whole")
        p, out = self.de(report, cond, contrast, "maxlfq", "maxlfq")
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        self.assertNotIn(CAUSE, p.stderr)


if __name__ == "__main__":
    unittest.main()
