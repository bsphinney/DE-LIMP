#!/usr/bin/env python3
"""Cysteine that was never alkylated must not be searched as carbamidomethyl.

Staff reports (2026-10-09): estimate_params.py always wrote --unimod4 (fixed carbamidomethyl C),
and --overrides could only add or replace a flag, so on-bead / IP / pull-down digests with no
reduction and alkylation were searched with every cysteine peptide unsearchable -- 0.03-0.5% of
identified peptides held a C, against 7.8% on an S-Trap prep with DTT/IAA. The fix used was a
hand `sed` of params.cfg plus a hand-edited rationale sidecar.

  * estimate_params.py --cys-alkylation cam|none: none leaves --unimod4 out (Sage: no static C),
    and the sidecar records which, and on whose word -- cam unasked is tagged DEFAULT;
  * the Methods name carbamidomethylation only when it was searched, and say why it was not;
  * sample_quality.py --report warns when carbamidomethyl was fixed and under 1% of the identified
    precursors hold a C.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import estimate_params as ep  # noqa: E402
import make_methods as mm  # noqa: E402
import sample_quality as sq  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_PYARROW = True
except ImportError:
    HAVE_PYARROW = False


def estimate(d, engine="diann", *more):
    out = os.path.join(d, "params.cfg" if engine == "diann" else "sage_config.json")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                        "--engine", engine, "--acquisition", "DIA", "--instrument", "timsTOF HT",
                        "--precursor-mz-range", "300", "1200", "--out", out, *more],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError(p.stderr)
    with open(out) as fh, open(out + ".rationale.json") as side:
        return out, fh.read(), json.load(side)


class EstimateParams(unittest.TestCase):
    def test_the_default_is_fixed_carbamidomethyl_tagged_as_a_default(self):
        with tempfile.TemporaryDirectory() as d:
            _, cfg, side = estimate(d)
            self.assertIn("--unimod4\n", cfg)
            self.assertEqual(side["cys_alkylation"]["value"], "cam")
            self.assertFalse(side["cys_alkylation"]["given"])
            self.assertTrue(side["cys_alkylation"]["source"].startswith("DEFAULT, not user-confirmed"))
            self.assertIn("DEFAULT", side["rationale"]["--unimod4"]["source"])

    def test_none_leaves_cysteine_unmodified_and_says_who_said_so(self):
        with tempfile.TemporaryDirectory() as d:
            path, cfg, side = estimate(d, "diann", "--cys-alkylation", "none")
            self.assertNotIn("--unimod4", cfg)
            self.assertEqual(side["cys_alkylation"]["value"], "none")
            self.assertTrue(side["cys_alkylation"]["given"])
            self.assertIn("not reduced and alkylated", side["cys_alkylation"]["source"])
            self.assertIn("flag omitted", side["rationale"]["--unimod4"]["value"])
            self.assertFalse(ep.cam_searched(path))
            _, cfg, side = estimate(d, "diann", "--cys-alkylation", "cam")
            self.assertIn("--unimod4\n", cfg)
            self.assertTrue(side["cys_alkylation"]["given"])

    def test_sage_none_has_no_static_cysteine(self):
        with tempfile.TemporaryDirectory() as d:
            path, text, side = estimate(d, "sage", "--cys-alkylation", "none")
            self.assertEqual(json.loads(text)["database"]["static_mods"], {})
            self.assertIn("unmodified", side["rationale"]["static_mods"]["source"])
            self.assertFalse(ep.cam_searched(path))
            path, text, _ = estimate(d, "sage")
            self.assertEqual(json.loads(text)["database"]["static_mods"], {"C": 57.0215})
            self.assertTrue(ep.cam_searched(path))

    def test_cam_searched_reads_every_way_a_cfg_can_fix_it(self):
        with tempfile.TemporaryDirectory() as d:
            def cfg(text):
                p = os.path.join(d, "x.cfg")
                with open(p, "w") as fh:
                    fh.write(text)
                return p
            self.assertTrue(ep.cam_searched(cfg("--qvalue 0.01\n--unimod4\n")))
            self.assertTrue(ep.cam_searched(cfg("--fixed-mod UniMod:4,57.021464,C\n")))
            self.assertTrue(ep.cam_searched(cfg("--fixed-mod Carbamidomethyl,57.021464,C\n")))
            self.assertFalse(ep.cam_searched(cfg("--qvalue 0.01\n# --unimod4\n")))
            self.assertIsNone(ep.cam_searched(os.path.join(d, "missing.cfg")))


class Methods(unittest.TestCase):
    def paragraph(self, d, *more):
        path, _, _ = estimate(d, "diann", *more)
        prov = os.path.join(d, "search_provenance.json")
        with open(prov, "w") as fh:
            json.dump({"engine": "diann", "version": "2.7.0", "params_file": path,
                       "result": {}}, fh)
        return mm.search_paragraph(mm.search_record(path, prov))

    def test_carbamidomethylation_is_named_only_when_it_was_searched(self):
        with tempfile.TemporaryDirectory() as d:
            text = self.paragraph(d)
            self.assertIn("Carbamidomethyl", text)
            self.assertNotIn("searched unmodified", text)
        with tempfile.TemporaryDirectory() as d:
            text = self.paragraph(d, "--cys-alkylation", "none")
            self.assertNotIn("Carbamidomethyl (", text)
            self.assertIn("No fixed modifications were set.", text)
            self.assertIn("Cysteine was searched unmodified (no carbamidomethylation), as the "
                          "samples were not reduced and alkylated.", text)


def write_report(path, seqs):
    """A DIA-NN-shaped report: one precursor per sequence, all at 1% FDR."""
    pq.write_table(pa.table({"Run": ["r1"] * len(seqs), "Protein.Group": ["P1"] * len(seqs),
                             "Precursor.Id": [f"{s}2" for s in seqs],
                             "Stripped.Sequence": seqs, "Q.Value": [0.001] * len(seqs)}), path)


@unittest.skipUnless(HAVE_PYARROW, "the cysteine check reads report.parquet with pyarrow")
class SampleQuality(unittest.TestCase):
    # 1 of 200 precursors with a C: 0.5%, the top of what the on-bead searches gave
    LOW = ["PEPTCDEK"] + [f"PEPTIDE{'K' * (i % 5 + 1)}{i}" for i in range(199)]
    NORMAL = [f"PEPCTIDE{i}" for i in range(16)] + [f"PEPTIDE{i}" for i in range(184)]

    def search(self, d, seqs, cys="cam"):
        params, _, _ = estimate(d, "diann", "--cys-alkylation", cys)
        report = os.path.join(d, "report.parquet")
        write_report(report, seqs)
        with open(os.path.join(d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "params_file": params}, fh)
        return report

    def test_fixed_cam_with_almost_no_cysteine_warns_and_says_the_fix(self):
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, self.LOW))
            self.assertEqual(rec["status"], "WARN")
            self.assertTrue(rec["cam_searched"])
            self.assertEqual((rec["n_cys"], rec["n"]), (1, 200))
            self.assertIn("NOT reduced and alkylated", rec["message"])
            self.assertIn("--cys-alkylation none", rec["message"])

    def test_a_normal_share_or_an_unmodified_search_does_not_warn(self):
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, self.NORMAL))
            self.assertEqual(rec["status"], "INFO")
            self.assertAlmostEqual(rec["fraction"], 0.08)
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, self.LOW, cys="none"))
            self.assertEqual(rec["status"], "INFO")
            self.assertFalse(rec["cam_searched"])
            self.assertIn("left cysteine unmodified", rec["message"])

    def test_no_params_found_is_said_not_guessed(self):
        with tempfile.TemporaryDirectory() as d:
            report = os.path.join(d, "report.parquet")
            write_report(report, self.LOW)
            rec = sq.cys_check(report)
            self.assertEqual(rec["status"], "INFO")
            self.assertIsNone(rec["cam_searched"])
            self.assertIn("pass --params", rec["message"])

    def test_the_warning_reaches_sample_quality_md(self):
        with tempfile.TemporaryDirectory() as d:
            report = self.search(d, self.LOW)
            matrix = os.path.join(d, "Expression_Matrix.csv")
            with open(matrix, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["Protein.Group", "Genes", "Protein.Names", "S1", "S2", "S3"])
                w.writerow(["P1", "ACTB", "ACTB_HUMAN", 16, 16, 16])
            out = os.path.join(d, "SAMPLE_QUALITY.md")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                                "--matrix", matrix, "--report", report, "--out", out],
                               capture_output=True, text=True, cwd=d)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                md = fh.read()
            self.assertIn("## Cysteine", md)
            self.assertIn("- CYSTEINE: only 0.50% of the identified precursors", md)
            self.assertEqual(json.loads(p.stdout)["cysteine"]["status"], "WARN")


if __name__ == "__main__":
    unittest.main()
