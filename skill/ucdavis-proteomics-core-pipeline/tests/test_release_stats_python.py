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
        # each mouse in one group: the between-block factor IS the group, tested on the 12
        # mouse means by exact permutation; no within-mouse question to ask
        self.assertIn("[between block (Mouse) levels: Ctrl vs Trt]", detail)
        self.assertIn("one value per block (Mouse) mean (exact over all 924 arrangements)", detail)
        self.assertIn("Ctrl: mean z", detail)
        self.assertIn("(n = 6)", detail)
        self.assertNotIn("within one Mouse", detail)
        self.assertIsNotNone(p)
        self.assertFalse(c, detail)

    def test_three_mice_per_group_is_too_few_for_a_test(self):
        z, gmap, samples, bmap = self.techrep(3, 5)
        _c, detail, p = sq.confound_check(z, gmap, samples, 1.5, bmap=bmap, block_name="Mouse")
        self.assertIsNone(p)                      # 20 arrangements of 6 means: no 1% test
        self.assertIn("best possible p 0.1 > 0.01 (exact over all 20 arrangements), so no test", detail)
        self.assertIn("(n = 3)", detail)

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
        self.assertIn("labels relabelled within each block", detail)
        self.assertIn("[within one block (Mouse)]", detail)
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


class WithinAndBetweenBlockQuestions(unittest.TestCase):
    """PROT_0756 v2 (michelle-rerun, job 24175676): 6 mice x 5 IPs, age between mice. The
    within-mouse relabelling test flagged hemolysis at p = 0.001 -- a BAIT effect (Kv2.1 IPs) --
    and the message read like an age confound it can never detect. Now the two questions are
    separate tests with separate words."""

    BAITS = ("JPH3", "JPH4", "Kv21", "RyR", "IgG")

    def prot0756(self, age_blood=3.0, bait_blood=2.0, mice_per_age=3, seed=7):
        rng = random.Random(seed)
        samples, gmap, bmap, vals = [], {}, {}, {}
        for m in range(1, 2 * mice_per_age + 1):
            age = "Old" if m <= mice_per_age else "Young"
            mouse = rng.gauss(0, 0.3) + (age_blood if age == "Old" else 0.0)   # per-mouse blood
            for b in self.BAITS:
                s = f"M{m}_{b}"
                samples.append(s); gmap[s] = f"{age}_{b}"; bmap[s] = f"Mouse{m}"
                vals[s] = mouse + (bait_blood if b == "Kv21" else 0.0) + rng.gauss(0, 0.2)
        return samples, gmap, bmap, vals

    def tests(self, **kw):
        samples, gmap, bmap, vals = self.prot0756(**kw)
        z = sq.zscore(vals, samples)
        return {t["question"]: t for t in sq.block_confound_checks(z, gmap, bmap, samples, 1.5, "Mouse")}

    def test_the_bait_effect_is_the_within_question(self):
        t = self.tests(age_blood=0.0)
        self.assertTrue(t["within"]["flag"])
        self.assertLess(t["within"]["p"], sq.CONFOUND_P)
        self.assertTrue(t["within"]["detail"].startswith("[within one block (Mouse)]"))
        self.assertFalse(t["between"]["flag"])       # no age effect, and 3 vs 3 cannot test

    def test_age_is_the_between_question_on_mouse_means(self):
        t = self.tests(bait_blood=0.0)
        b = t["between"]
        self.assertEqual(b["factor"], "Old vs Young")
        self.assertIsNone(b["p"])                     # 3 vs 3 mice: best p 2/20 = 0.1
        self.assertAlmostEqual(b["min_p"], 0.1)
        self.assertFalse(b["flag"])
        self.assertTrue(b["caution"])
        self.assertIn("complete separation (all Old above all Young)", b["detail"])
        self.assertFalse(t["within"]["flag"])        # no bait effect

    def test_enough_mice_make_age_testable(self):
        t = self.tests(bait_blood=0.0, mice_per_age=8)
        b = t["between"]
        self.assertTrue(b["flag"])
        self.assertLess(b["p"], sq.CONFOUND_P)
        self.assertIn("one value per block (Mouse) mean (exact over all 12,870 arrangements)",
                      b["detail"])

    def test_a_paired_design_too_small_to_test_gets_a_within_caution(self):
        # review round 2: 6 mice x (Bait, IgG) = 64 arrangements, best p 2/64 = 0.031 > 0.01,
        # so no test -- and a clean +2 SD bait effect used to be SILENT. Now: a caution.
        rng = random.Random(4)
        samples, gmap, bmap, vals = [], {}, {}, {}
        for m in range(6):
            eff = rng.gauss(0, 0.3)
            for g in ("IgG", "Bait"):
                s = f"M{m}_{g}"
                samples.append(s); gmap[s] = g; bmap[s] = f"M{m}"
                vals[s] = eff + (2.0 if g == "Bait" else 0.0) + rng.gauss(0, 0.1)
        z = sq.zscore(vals, samples)
        (t,) = sq.confound_tests(z, gmap, samples, 1.5, bmap=bmap, block_name="Mouse")
        self.assertGreaterEqual(t["gap"], 1.5)       # the same |z| gap rule as everywhere else
        self.assertEqual(t["question"], "within")
        self.assertFalse(t["flag"])
        self.assertTrue(t["caution"])
        self.assertAlmostEqual(t["min_p"], 2 / 64)
        self.assertIn("every block orders the groups the same way", t["detail"])
        self.assertEqual((t["hi"], t["lo"]), ("Bait", "IgG"))

    def test_a_paired_design_has_no_between_question(self):
        rng = random.Random(1)
        samples, gmap, bmap, vals = [], {}, {}, {}
        for p in range(6):
            eff = rng.gauss(0, 2)
            for g in ("Ctrl", "Trt"):
                s = f"P{p}_{g}"
                samples.append(s); gmap[s] = g; bmap[s] = f"P{p}"; vals[s] = eff + rng.gauss(0, 0.2)
        z = sq.zscore(vals, samples)
        qs = [t["question"] for t in sq.block_confound_checks(z, gmap, bmap, samples, 1.5, "Patient")]
        self.assertEqual(qs, ["within"])

    def test_the_report_words_each_question_apart(self):
        samples, gmap, bmap, vals = self.prot0756()        # age AND bait blood
        genes = ["HBB", "HBA1", "CA1", "CA2", "CAT", "AHSP"]
        rng = random.Random(3)
        with tempfile.TemporaryDirectory() as tmp:
            m = os.path.join(tmp, "Expression_Matrix.csv")
            rows = [[f"P{i}", g] + [f"{20 + vals[s] + rng.gauss(0, 0.1):.3f}" for s in samples]
                    for i, g in enumerate(genes)]
            rows += [[f"Q{i}", f"G{i}"] + [f"{18 + rng.gauss(0, 0.3):.3f}" for _ in samples]
                     for i in range(60)]
            write_csv(m, ["Protein.Group", "Genes"] + samples, rows)
            c = os.path.join(tmp, "conditions.csv")
            write_csv(c, ["File.Name", "Group", "Mouse"], [(s, gmap[s], bmap[s]) for s in samples])
            with open(os.path.join(tmp, "de_provenance.json"), "w") as fh:
                json.dump({"block_column": "Mouse"}, fh)          # the block comes from here
            out = os.path.join(tmp, "SAMPLE_QUALITY.md")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                                "--matrix", m, "--conditions", c, "--out", out, "--taxid", "9606"],
                               capture_output=True, text=True, cwd=tmp)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                md = fh.read()
        self.assertIn("**HEMOLYSIS differs among the samples of one block (Mouse)**", md)
        self.assertIn("This answers the WITHIN-Mouse question only", md)
        self.assertIn("HEMOLYSIS: complete separation between block (Mouse) levels -- all Old above all "
                      "Young", md)
        self.assertIn("this is a caution, not a finding", md)
        self.assertNotIn("CONFOUNDED WITH GROUP", md)
        self.assertIn("_differs among the samples of one block? (within-block contrasts)_", md)
        self.assertIn("_differs between block levels, Old vs Young? (between-block contrasts)_", md)
        # review round 2, grammar: no "Mouses", no "too few sample, relabelled ..."
        self.assertNotIn("Mouses", md)
        self.assertNotIn("too few sample", md)


class OneRuleForWeakEvidence(unittest.TestCase):
    """Review round 2: 3 vs 3 complete separation happens in ~4% of null panels (1 in 10
    arrangements, halved by the direction): it must never be a flag, blocked or not."""

    def test_null_three_vs_three_never_flags(self):
        rng = random.Random(11)
        samples = [f"S{i}" for i in range(6)]
        gmap = {s: ("A" if i < 3 else "B") for i, s in enumerate(samples)}
        flags = cautions = 0
        for _ in range(300):
            z = sq.zscore({s: rng.gauss(0, 1) for s in samples}, samples)
            (t,) = sq.confound_tests(z, gmap, samples, 1.5)
            flags += t["flag"]; cautions += t["caution"]
        self.assertEqual(flags, 0)
        self.assertGreater(cautions, 0)            # it does happen -- as a caution

    def test_the_report_says_caution_not_confounded(self):
        rows_z = [0.1, 0.2, 0.0, 2.1, 2.3, 2.2]
        samples = [f"S{i}" for i in range(6)]
        genes = ["HBB", "HBA1", "CA1", "CA2", "CAT", "AHSP"]
        rng = random.Random(2)
        with tempfile.TemporaryDirectory() as tmp:
            m = os.path.join(tmp, "Expression_Matrix.csv")
            rows = [[f"P{i}", g] + [f"{20 + 1.5 * v + rng.gauss(0, 0.05):.3f}" for v in rows_z]
                    for i, g in enumerate(genes)]
            rows += [[f"Q{i}", f"G{i}"] + [f"{18 + rng.gauss(0, 0.3):.3f}" for _ in samples]
                     for i in range(60)]
            write_csv(m, ["Protein.Group", "Genes"] + samples, rows)
            c = os.path.join(tmp, "conditions.csv")
            write_csv(c, ["File.Name", "Group"], [(s, "A" if i < 3 else "B") for i, s in enumerate(samples)])
            out = os.path.join(tmp, "SAMPLE_QUALITY.md")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                                "--matrix", m, "--conditions", c, "--out", out, "--taxid", "9606"],
                               capture_output=True, text=True, cwd=tmp)
            self.assertEqual(p.returncode, 0, p.stderr)
            md = open(out).read()
            summary = json.loads(p.stdout)
        self.assertIn("HEMOLYSIS: complete separation between groups -- all B above all A", md)
        self.assertIn("so this is a caution, not a finding", md)
        self.assertNotIn("CONFOUNDED WITH GROUP", md)
        self.assertEqual(summary["group_confounded"], [])


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
