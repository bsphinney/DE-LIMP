#!/usr/bin/env python3
"""
DIA-NN's per-run Normalisation.Instability was never read (a staff report): a timsTOF tissue
cohort had 0.67-1.00 in every run of report.stats.tsv (an E. coli cohort: 0.04-0.07), and that normalisation passed silently into the DE, which uses DIA-NN's
normalised quantities. check_report_runs.normalisation_instability() now reads it (the one
reader and threshold), and audit_results.py warns on it, naming the runs, with the alternatives
DIA-NN's README documents (--global-norm, --no-norm, the non-normalised Precursor.Quantity).

Stats files are fixtures with DIA-NN 2.7.0's column layout (read off a real one on HIVE) and
generic run names. Hermetic: no network, no real data.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import check_report_runs as crr  # noqa: E402

AUDIT = os.path.join(SCRIPTS, "audit_results.py")
# DIA-NN 2.7.0's report.stats.tsv header, in its order
HEADER = ["File.Name", "Precursors.Identified", "Proteins.Identified", "Total.Quantity",
          "MS1.Signal", "MS2.Signal", "FWHM.Scans", "FWHM.RT", "Median.Mass.Acc.MS1",
          "Median.Mass.Acc.MS1.Corrected", "Median.Mass.Acc.MS2", "Median.Mass.Acc.MS2.Corrected",
          "Normalisation.Instability", "Median.RT.Prediction.Acc", "Average.Peptide.Length",
          "Average.Peptide.Charge", "Average.Missed.Tryptic.Cleavages"]


def write_stats(path, rows, header=HEADER):
    """rows: (run, precursors, instability)."""
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for run, prec, ni in rows:
            vals = {"File.Name": f"/data/raw/{run}.d", "Precursors.Identified": prec,
                    "Proteins.Identified": prec // 10, "Normalisation.Instability": ni}
            fh.write("\t".join(str(vals.get(h, 1)) for h in header) + "\n")


UNSTABLE = [("s01", 9000, 0.95), ("s02", 8800, 0.67), ("s03", 9100, 1.0), ("s04", 8700, 0.71)]
STABLE = [("s01", 9000, 0.04), ("s02", 8800, 0.07), ("s03", 9100, 0.05)]


class TheReader(unittest.TestCase):
    def test_values_per_run_and_the_runs_above_the_threshold(self):
        with tempfile.TemporaryDirectory() as d:
            write_stats(os.path.join(d, "report.stats.tsv"),
                        UNSTABLE + [("s05", 9000, 0.05)])
            ni = crr.normalisation_instability(os.path.join(d, "report.parquet"))
        self.assertTrue(ni["column"])
        self.assertEqual(ni["per_run"]["s01"], 0.95)
        self.assertEqual(set(ni["flagged"]), {"s01", "s02", "s03", "s04"})
        self.assertEqual(ni["median"], 0.71)
        self.assertEqual(ni["max"], 1.0)
        self.assertEqual(ni["threshold"], crr.NORM_INSTABILITY_WARN)
        self.assertIn("documents no cutoff", ni["basis"])

    def test_only_the_analysed_runs_and_never_a_run_with_no_identifications(self):
        with tempfile.TemporaryDirectory() as d:
            write_stats(os.path.join(d, "report.stats.tsv"),
                        STABLE + [("blank1", 0, 0.0), ("wash1", 12, 1.0)])
            ni = crr.normalisation_instability(os.path.join(d, "report.parquet"),
                                               runs=["s01", "s02.d", "blank1"])
        self.assertEqual(set(ni["per_run"]), {"s01", "s02"},
                         "wash1 was not analysed; blank1 identified nothing")
        self.assertEqual(ni["flagged"], {})

    def test_no_stats_file_and_no_column(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(crr.normalisation_instability(os.path.join(d, "report.parquet")))
            write_stats(os.path.join(d, "report.stats.tsv"), STABLE,
                        header=[h for h in HEADER if h != "Normalisation.Instability"])
            ni = crr.normalisation_instability(os.path.join(d, "report.parquet"))
        self.assertFalse(ni["column"])
        self.assertEqual(ni["per_run"], {})

    def test_a_diann_1x_tsv_report_names_its_stats_file(self):
        self.assertEqual(crr.stats_path("/x/report.tsv"), "/x/report.stats.tsv")
        self.assertEqual(crr.stats_path("/x/report.parquet"), "/x/report.stats.tsv")


class TheAudit(unittest.TestCase):
    def audit(self, d, rows=UNSTABLE, engine="diann", stats=True, conditions=None):
        out = os.path.join(d, "search")
        os.makedirs(out)
        if stats:
            write_stats(os.path.join(out, "report.stats.tsv"), rows)
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": engine, "search_mode": "single_shot"}, fh)
        tables = os.path.join(d, "tables")
        os.makedirs(tables)
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump({"input": os.path.join(out, "report.parquet"), "adjp": 0.05,
                       "logfc": 1}, fh)
        argv = [sys.executable, AUDIT, "--out", "AUDIT.md", "--de-dir", tables]
        if conditions:
            cpath = os.path.join(d, "conditions.csv")
            with open(cpath, "w") as fh:
                fh.write("File.Name,Group\n" + "".join(f"{r},{g}\n" for r, g in conditions))
            argv += ["--conditions", cpath]
        res = subprocess.run(argv, capture_output=True, text=True, timeout=120, cwd=d)
        self.assertEqual(res.returncode, 0, res.stderr)
        with open(os.path.join(d, "AUDIT.json")) as fh:
            found = [f for f in json.load(fh)["findings"]
                     if f["check"] == "normalisation_stability"]
        with open(os.path.join(d, "AUDIT.md")) as fh:
            return found, fh.read()

    def test_unstable_runs_are_a_warning_naming_each_run_and_value(self):
        with tempfile.TemporaryDirectory() as d:
            found, md = self.audit(d)
        self.assertEqual([f["status"] for f in found], ["WARN"])
        msg = found[0]["message"]
        self.assertIn("above 0.3 in 4 of 4 runs", msg)
        self.assertIn("s03 1.00", msg)
        self.assertIn("s02 0.67", msg)
        self.assertIn("median of all 4: 0.83", msg)
        self.assertIn("normalisation_stability", md)
        # the value is recorded, every run of it
        self.assertEqual(found[0]["detail"]["per_run"]["s04"], 0.71)
        self.assertEqual(found[0]["detail"]["threshold"], crr.NORM_INSTABILITY_WARN)

    def test_the_warning_offers_only_what_diann_documents(self):
        with tempfile.TemporaryDirectory() as d:
            found, _ = self.audit(d)
        msg = found[0]["message"]
        for documented in ("--global-norm", "--no-norm", "Precursor.Quantity"):
            self.assertIn(documented, msg)
        # DIA-NN's flags are named as DIA-NN's, and no invented one appears
        flags = {w.strip("(),.;") for w in msg.split() if w.startswith("--")}
        self.assertEqual(flags, {"--global-norm", "--no-norm"})
        self.assertIn("documents no cutoff", msg)

    def test_the_conditions_decide_which_runs_count(self):
        with tempfile.TemporaryDirectory() as d:
            found, _ = self.audit(d, rows=STABLE + [("unused", 9000, 0.9)],
                                  conditions=[("s01", "A"), ("s02", "A"), ("s03", "B")])
        self.assertEqual([f["status"] for f in found], ["PASS"])
        self.assertIn("all 3 runs", found[0]["message"])

    def test_a_diann_search_without_its_stats_file_says_not_assessed(self):
        with tempfile.TemporaryDirectory() as d:
            found, _ = self.audit(d, stats=False)
        self.assertEqual([f["status"] for f in found], ["INFO"])
        self.assertIn("report.stats.tsv", found[0]["message"])

    def test_another_engines_report_is_not_asked(self):
        with tempfile.TemporaryDirectory() as d:
            found, _ = self.audit(d, engine="sage", stats=False)
        self.assertEqual(found, [])


if __name__ == "__main__":
    unittest.main()
