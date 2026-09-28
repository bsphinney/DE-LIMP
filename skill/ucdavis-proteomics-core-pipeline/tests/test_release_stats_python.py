#!/usr/bin/env python3
"""2.8.0 release review, the Python side of the R/stat core (pure Python: runs in CI).

  * sample_quality.py maps samples to groups by EXACT stem first; a substring only when
    exactly one row matches ('S1' is in 'S10'..'S19': S10 used to get S1's group);
  * its confound permutation test respects the block (the animal): blocks nested in the
    groups -> the test runs on block means; blocks spanning groups -> labels are relabelled
    within each block. It reads the block from --block or run_de.R's de_provenance.json;
  * audit_results.py counts replicates as distinct blocks (mice), not runs, when the DE was
    blocked or conditions.csv has a subject column whose values look like subjects;
  * audit_results.py takes adjp / logfc from de_provenance.json (one cutoff definition) and
    says where they came from;
  * collect_conditions.py: a "Block" column is a blocking unit, not Batch (as Batch,
    `--block Batch` is refused).
"""
import csv
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import sample_quality as sq   # noqa: E402
from test_collect_conditions_subject import run_map   # noqa: E402

TECHREP = [("Ctrl_m1_inj1", "Ctrl", "Ctrl_M1"), ("Ctrl_m1_inj2", "Ctrl", "Ctrl_M1"),
           ("Ctrl_m1_inj3", "Ctrl", "Ctrl_M1"), ("Ctrl_m2_inj1", "Ctrl", "Ctrl_M2"),
           ("Ctrl_m2_inj2", "Ctrl", "Ctrl_M2"), ("Ctrl_m2_inj3", "Ctrl", "Ctrl_M2"),
           ("Trt_m1_inj1", "Trt", "Trt_M1"), ("Trt_m1_inj2", "Trt", "Trt_M1"),
           ("Trt_m1_inj3", "Trt", "Trt_M1"), ("Trt_m2_inj1", "Trt", "Trt_M2"),
           ("Trt_m2_inj2", "Trt", "Trt_M2"), ("Trt_m2_inj3", "Trt", "Trt_M2")]


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(header); w.writerows(rows)


class SampleQualityMapping(unittest.TestCase):
    def test_exact_match_beats_a_substring(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "conditions.csv")
            write_csv(p, ["File.Name", "Group"],
                      [(f"exp_S{i}", "Ctrl" if i <= 6 else "Trt") for i in range(1, 13)])
            g = sq.load_conditions(p, [f"exp_S{i}" for i in range(1, 13)])
        self.assertEqual(g["exp_S10"], "Trt")          # used to be S1's Ctrl
        self.assertEqual(g["exp_S1"], "Ctrl")
        self.assertEqual(sum(v == "Trt" for v in g.values()), 6)

    def test_a_substring_is_used_only_when_unique(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "conditions.csv")
            write_csv(p, ["File.Name", "Group"], [("S1", "Ctrl"), ("S10", "Trt")])
            g = sq.load_conditions(p, ["2026_S10_run", "2026_S1_run"])
        # '2026_s10_run' contains both 's1' and 's10': ambiguous -> not assigned by guess
        self.assertNotIn("2026_S10_run", g)
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "conditions.csv")
            write_csv(p, ["File.Name", "Group"], [("Ctrl_A", "Ctrl"), ("Trt_B", "Trt")])
            g = sq.load_conditions(p, ["2026_Ctrl_A.raw", "2026_Trt_B.raw"])
        self.assertEqual(g, {"2026_Ctrl_A.raw": "Ctrl", "2026_Trt_B.raw": "Trt"})

    def test_the_block_column_maps_like_the_group(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "conditions.csv")
            write_csv(p, ["File.Name", "Group", "Mouse"], TECHREP)
            b = sq.load_conditions(p, [r[0] for r in TECHREP], column="Mouse")
        self.assertEqual(b["Trt_m2_inj3"], "Trt_M2")

    def test_the_block_is_read_from_de_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            m = os.path.join(tmp, "Expression_Matrix.csv")
            open(m, "w").close()
            self.assertIsNone(sq.block_column_for(m))
            with open(os.path.join(tmp, "de_provenance.json"), "w") as fh:
                json.dump({"block_column": "Mouse"}, fh)
            self.assertEqual(sq.block_column_for(m), "Mouse")
            self.assertEqual(sq.block_column_for(m, "Subject"), "Subject")   # --block wins


class ConfoundTestRespectsTheBlock(unittest.TestCase):
    @staticmethod
    def techrep(mice_per_group, seed):
        # a big MOUSE effect, three injections per mouse, no group effect
        rng = random.Random(seed)
        samples, gmap, bmap, vals = [], {}, {}, {}
        for g in ("Ctrl", "Trt"):
            for m in range(mice_per_group):
                eff = rng.gauss(0, 2)
                for i in range(3):
                    s = f"{g}_m{m}_{i}"
                    samples.append(s); gmap[s] = g; bmap[s] = f"{g}_M{m}"
                    vals[s] = eff + rng.gauss(0, 0.1)
        return sq.zscore(vals, samples), gmap, samples, bmap

    def test_technical_replicates_are_tested_on_mouse_means(self):
        z, gmap, samples, bmap = self.techrep(6, 5)
        c, detail, p = sq.confound_check(z, gmap, samples, 1.5, bmap=bmap, block_name="Mouse",
                                         n_perm=999)
        self.assertIn("on Mouse means (one value per Mouse)", detail)
        self.assertIn("(n=6 Mouse)", detail)
        self.assertIsNotNone(p)
        self.assertFalse(c, detail)

    def test_three_mice_per_group_is_too_few_for_a_test(self):
        z, gmap, samples, bmap = self.techrep(3, 5)
        _c, detail, p = sq.confound_check(z, gmap, samples, 1.5, bmap=bmap, block_name="Mouse")
        self.assertIsNone(p)                      # 20 arrangements of 6 means: no 1% test
        self.assertIn("too few Mouse means for a test", detail)
        self.assertIn("(n=3 Mouse)", detail)

    def test_blocks_spanning_groups_are_relabelled_within_block(self):
        # 8 mice, each giving one IP to each of 3 groups; mice differ a lot, groups do not
        rng = random.Random(9)
        samples, gmap, bmap, vals = [], {}, {}, {}
        for m in range(8):
            eff = rng.gauss(0, 3)
            for g in ("A", "B", "C"):
                s = f"m{m}_{g}"
                samples.append(s); gmap[s] = g; bmap[s] = f"M{m}"
                vals[s] = eff + rng.gauss(0, 0.3)
        z = sq.zscore(vals, samples)
        c, detail, p = sq.confound_check(z, gmap, samples, 1.5, bmap=bmap, block_name="Mouse",
                                         n_perm=999)
        self.assertIn("relabelling within each Mouse", detail)
        self.assertIsNotNone(p)
        self.assertFalse(c, detail)
        # a real group effect inside every mouse is still found
        for s in samples:
            if gmap[s] == "C":
                vals[s] += 1.5
        z = sq.zscore(vals, samples)
        c, detail, p = sq.confound_check(z, gmap, samples, 1.5, bmap=bmap, block_name="Mouse",
                                         n_perm=999)
        self.assertTrue(c, detail)

    def test_without_a_block_nothing_changes(self):
        rng = random.Random(3)
        samples = [f"S{i}" for i in range(30)]
        vals = {s: rng.gauss(0, 1) + (2.5 if (i // 3) % 2 else 0) for i, s in enumerate(samples)}
        z = sq.zscore(vals, samples)
        gmap = {s: f"G{i // 3}" for i, s in enumerate(samples)}
        a = sq.confound_check(z, gmap, samples, 1.5)
        b = sq.confound_check(z, gmap, samples, 1.5, bmap={})
        self.assertEqual(a, b)


class AuditReplicationAndCutoffs(unittest.TestCase):
    def audit(self, tmp, cond_rows, header, prov=None, extra=()):
        cond = os.path.join(tmp, "conditions.csv")
        write_csv(cond, header, cond_rows)
        de = os.path.join(tmp, "tables"); os.makedirs(de, exist_ok=True)
        if prov is not None:
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump(prov, fh)
        out = os.path.join(tmp, "AUDIT.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"), "--out", out,
                            "--conditions", cond, "--de-dir", de, *extra],
                           capture_output=True, text=True, cwd=tmp)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(tmp, "AUDIT.json")) as fh:
            return json.load(fh)

    def replication(self, res):
        return next(f for f in res["findings"] if f["check"] == "replication")

    def test_injections_of_one_mouse_are_one_replicate(self):
        with tempfile.TemporaryDirectory() as tmp:
            r = self.replication(self.audit(tmp, TECHREP, ["File.Name", "Group", "Mouse"]))
        self.assertEqual(r["status"], "WARN")             # 2 mice per group, not 6 runs
        self.assertIn("only 2 distinct Mouse (n = Mouse, not runs)", r["message"])
        self.assertEqual(r["detail"]["unit"], "Mouse")
        self.assertEqual(r["detail"]["runs_per_group"], {"Ctrl": 6, "Trt": 6})

    def test_the_block_run_de_recorded_is_the_unit(self):
        rows = [(f"r{i}", "A" if i < 6 else "B", f"U{i % 3}") for i in range(12)]
        with tempfile.TemporaryDirectory() as tmp:
            r = self.replication(self.audit(tmp, rows, ["File.Name", "Group", "Unit"],
                                            prov={"block_column": "Unit"}))
        self.assertEqual(r["detail"]["unit"], "Unit")
        self.assertEqual(r["detail"]["group_sizes"], {"A": 3, "B": 3})

    def test_a_subject_header_over_non_subject_values_is_not_the_unit(self):
        rows = [(f"r{i}", "A" if i < 6 else "B", "M" if i % 2 else "F") for i in range(12)]
        with tempfile.TemporaryDirectory() as tmp:
            r = self.replication(self.audit(tmp, rows, ["File.Name", "Group", "Subject"]))
        self.assertEqual(r["detail"]["unit"], "run")
        self.assertEqual(r["status"], "PASS")

    def test_cutoffs_come_from_de_provenance(self):
        rows = [(f"r{i}", "A" if i < 3 else "B") for i in range(6)]
        with tempfile.TemporaryDirectory() as tmp:
            res = self.audit(tmp, rows, ["File.Name", "Group"], prov={"adjp": 0.1, "logfc": 0.58})
            self.assertEqual(res["cutoffs"], {"adjp": 0.1, "logfc": 0.58, "source": {
                "adjp": "de_provenance.json", "logfc": "de_provenance.json"}})
        with tempfile.TemporaryDirectory() as tmp:
            res = self.audit(tmp, rows, ["File.Name", "Group"], prov={"adjp": 0.1},
                             extra=["--adjp", "0.01"])
            self.assertEqual(res["cutoffs"]["adjp"], 0.01)
            self.assertEqual(res["cutoffs"]["source"]["adjp"], "command line")
            self.assertEqual(res["cutoffs"]["source"]["logfc"],
                             "DEFAULT -- not recorded in de_provenance.json")


class BlockHeaderIsNotBatch(unittest.TestCase):
    def test_a_block_column_is_kept_for_block_not_filed_as_batch(self):
        runs = [f"R{i:02d}" for i in range(1, 13)]
        rows = [[f"R{i:02d}", "Ctrl" if i <= 6 else "Trt", f"B{(i - 1) % 6 + 1}"] for i in range(1, 13)]
        with tempfile.TemporaryDirectory() as tmp:
            rep, out = run_map(tmp, rows, ["sample", "group", "Block"], runs=runs)
        self.assertEqual(rep["block_column"], "Block")
        self.assertTrue(rep["block_suggested"])
        self.assertEqual(list(out[0]), ["File.Name", "Group", "Block"])   # not Batch


if __name__ == "__main__":
    unittest.main()
