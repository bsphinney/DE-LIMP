#!/usr/bin/env python3
"""run_de.R fixes from the 2.8.0 release review (R/stat core).

  * every DE table carries per-group detection -- Detected_<group> = "k/n" measured runs, from
    Detection_Matrix.csv -- and an Evidence label (measured in both / presence call / partly
    inferred|missing), described in de_provenance.json;
  * no replicates (no residual degrees of freedom) stops BEFORE quantification, with a message;
  * a non-dpc run removes the dpc-only files an earlier dpc run left in the same outdir;
  * the DE-LIMP session names the skill version that wrote it, not a hard-coded "DE-LIMP v2.5";
  * a column that splits a group into units of several runs (technical replicates) is pointed
    at --block.

End to end on the synthetic report of test_run_de_contaminants (A x3, B x3); skips without R.
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
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, HERE)

from test_run_de_contaminants import SYNTH_R, r_has, rscript   # noqa: E402

RUN_DE = os.path.join(SCRIPTS, "run_de.R")
NEEDS = ("limpa", "limma", "arrow", "dplyr", "tidyr", "jsonlite")


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


@unittest.skipUnless(r_has(*NEEDS), "needs R with " + "/".join(NEEDS))
class RunDeReleaseFixes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        cond = read_csv(os.path.join(cls.tmp, "conditions.csv"))
        # one run per group: no residual degrees of freedom
        with open(os.path.join(cls.tmp, "norep.csv"), "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["File.Name", "Group"])
            w.writerow(["run01", "A"]); w.writerow(["run04", "B"])
        # technical replicates: two runs of M1 and one of M2 in A, the same in B
        with open(os.path.join(cls.tmp, "techrep.csv"), "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["File.Name", "Group", "Mouse"])
            for r, m in zip(cond, ["M1", "M1", "M2", "M3", "M3", "M4"]):
                w.writerow([r["File.Name"], r["Group"], m])
        cls.out = os.path.join(cls.tmp, "de")
        cls.runs = {}
        for key, meta, method in (("dpc", "conditions.csv", "dpc"),
                                  ("maxlfq_same_dir", "conditions.csv", "maxlfq"),
                                  ("norep", "norep.csv", "dpc"),
                                  ("techrep", "techrep.csv", "maxlfq")):
            out = cls.out if key in ("dpc", "maxlfq_same_dir") else os.path.join(cls.tmp, key)
            if key == "dpc":
                cls.dpc_files = None
            p = subprocess.run(["Rscript", RUN_DE, "--input", os.path.join(cls.tmp, "report.parquet"),
                                "--metadata", os.path.join(cls.tmp, meta), "--method", method,
                                "--outdir", out], capture_output=True, text=True, cwd=cls.tmp)
            if key == "dpc" and p.returncode == 0:
                cls.dpc_tables = {f: read_csv(os.path.join(out, f)) for f in
                                  ("DE_dpc_B.A.csv", "Detection_Matrix.csv")}
                with open(os.path.join(out, "de_provenance.json")) as fh:
                    cls.dpc_prov = json.load(fh)
                cls.dpc_files = sorted(os.listdir(out))
                cls.app_version = rscript(
                    f'cat(readRDS("{os.path.join(out, "DE-LIMP_session.rds")}")$app_version)')
            cls.runs[key] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def ok(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        return p, out

    # -- per-group detection ---------------------------------------------------------------
    def test_de_tables_carry_per_group_detection(self):
        self.ok("dpc")
        rows = self.dpc_tables["DE_dpc_B.A.csv"]
        self.assertEqual(list(rows[0])[-3:], ["Detected_B", "Detected_A", "Evidence"])
        det = {r["Protein.Group"]: r for r in self.dpc_tables["Detection_Matrix.csv"]}
        a_runs, b_runs = ["run01", "run02", "run03"], ["run04", "run05", "run06"]
        seen = set()
        for r in rows:
            d = det[r["Protein.Group"]]
            ka = sum(float(d[c]) > 0 for c in a_runs)
            kb = sum(float(d[c]) > 0 for c in b_runs)
            self.assertEqual(r["Detected_A"], f"{ka}/3")
            self.assertEqual(r["Detected_B"], f"{kb}/3")
            want = ("presence call" if min(ka, kb) == 0 else
                    "measured in both" if ka == kb == 3 else "partly inferred")
            self.assertEqual(r["Evidence"], want, r["Protein.Group"])
            seen.add(want)
        self.assertEqual(seen, {"presence call", "measured in both", "partly inferred"})

    def test_the_columns_are_described_in_the_record(self):
        self.ok("dpc")
        dc = self.dpc_prov["detection_matrix"]["de_columns"]
        self.assertIn("k/n", dc["Detected_<group>"])
        self.assertIn("'presence call' = never measured in at least one group", dc["Evidence"])
        self.assertIn("'partly inferred'", dc["Evidence"])

    def test_maxlfq_says_missing_not_inferred(self):
        _p, out = self.ok("maxlfq_same_dir")
        rows = read_csv(os.path.join(out, "DE_maxlfq_B.A.csv"))
        labels = {r["Evidence"] for r in rows}
        self.assertTrue(labels <= {"measured in both", "presence call", "partly missing", ""}, labels)
        with open(os.path.join(out, "de_provenance.json")) as fh:
            self.assertIn("'partly missing'",
                          json.load(fh)["detection_matrix"]["de_columns"]["Evidence"])

    # -- stale dpc files -------------------------------------------------------------------
    def test_a_maxlfq_run_removes_the_dpc_only_files(self):
        self.ok("dpc")
        self.assertIn("QC_detected_vs_inferred.csv", self.dpc_files)
        self.assertIn("DE-LIMP_session.rds", self.dpc_files)
        p, out = self.ok("maxlfq_same_dir")
        self.assertFalse(os.path.exists(os.path.join(out, "QC_detected_vs_inferred.csv")))
        self.assertFalse(os.path.exists(os.path.join(out, "DE-LIMP_session.rds")))
        self.assertIn("removed QC_detected_vs_inferred.csv: left by an earlier --method dpc run",
                      p.stderr)

    # -- no replicates ---------------------------------------------------------------------
    def test_no_residual_df_stops_before_quantification(self):
        p, _out = self.runs["norep"]
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("No residual degrees of freedom in this design: 2 samples, 2 model "
                      "coefficients", p.stderr)
        self.assertIn("one sample in: A, B", p.stderr)
        self.assertNotIn("Quantifying proteins", p.stdout + p.stderr)

    # -- session version -------------------------------------------------------------------
    def test_the_session_names_the_skill_version(self):
        self.ok("dpc")
        with open(os.path.join(SKILL, ".claude-plugin", "plugin.json")) as fh:
            v = json.load(fh)["version"]
        self.assertEqual(self.app_version, f"ucdavis-proteomics-core-pipeline run_de.R v{v}")

    # -- technical replicates hint ---------------------------------------------------------
    def test_technical_replicates_point_at_block(self):
        p, _out = self.ok("techrep")
        self.assertIn("metadata column(s) Mouse repeat within a group", p.stderr)
        self.assertIn("pass --block Mouse", p.stderr)


if __name__ == "__main__":
    unittest.main()
