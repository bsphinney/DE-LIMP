#!/usr/bin/env python3
"""sample_quality.py's marker panels see the searched organism's own gene names.

msalemi, 2026-09-24 (Silva08172026, mouse brain IPs, UP000000589): SAMPLE_QUALITY.md said
HEMOLYSIS / SKELETAL_MUSCLE / EPIDERMIS "0 panel proteins detected - not assessable" while
Hba (P01942) and Hbb-bs (A8DUK4) were quantified at ~15 log2. The panels are human symbols;
mouse Hbb-bs / Hbb-b1 / Hba-a1 never matched, bovine serum haemoglobin (Cont_P02070) did,
and "not assessable" hid that the check had not run at all.

Guards (stdlib only):
  * the ortholog table is ONE definition keyed on the panels' own human symbols;
  * the mouse haemoglobins of that real matrix (rows copied verbatim) match HEMOLYSIS;
  * a Cont_-tagged marker is reported as "contaminant, not sample" and never scored --
    also when run_de.R already removed it (read from its record beside the matrix);
  * a matrix no panel matches at all says "check could not run", while a single empty
    panel on a matrix the panels do match is a real absence.
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

import sample_quality as sq  # noqa: E402

SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]
HEADER = ["Protein.Group", "Genes", "Protein.Names"] + SAMPLES   # run_de.R's column order
FLAT = [16.0] * 6
# Verbatim identifiers from Silva08172026's Expression_Matrix.csv (values are synthetic).
MOUSE_ROWS = [["A8DUK4", "Hbb-bs", "A8DUK4_MOUSE"] + [20, 15, 15, 15, 15, 15],
              ["P01942", "Hba", "HBA_MOUSE"] + [20, 15, 15, 15, 15, 15],
              ["F6UYE3", "Hbq1b", "F6UYE3_MOUSE"] + FLAT,          # theta: not an adult marker
              ["P24270", "Cat", "CATA_MOUSE"] + FLAT,
              ["P58252", "Eef2", "EF2_MOUSE"] + FLAT]


def write_matrix(path, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        w.writerows(rows)


class Base(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.matrix = os.path.join(self.root, "Expression_Matrix.csv")

    def tearDown(self):
        self._td.cleanup()

    def run_sq(self, rows, *extra):
        write_matrix(self.matrix, rows)
        out = os.path.join(self.root, "SAMPLE_QUALITY.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                            "--matrix", self.matrix, "--out", out, *extra],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out.replace(".md", ".json")) as fh, open(out) as md:
            return json.load(fh), md.read()


class OrthologTable(unittest.TestCase):
    def test_keys_are_the_panels_own_human_symbols(self):
        human = {g for v in sq.PANELS.values() for g in v}
        for taxid, table in sq.PANEL_ORTHOLOGS.items():
            self.assertIn(taxid, sq.PANEL_ORGANISMS)
            self.assertLessEqual(set(table), human, f"taxid {taxid} aliases a non-panel gene")
            for g, names in table.items():
                # Only names that differ beyond case belong here; the rest match already.
                self.assertTrue(all(n.upper() != g for n in names), (taxid, g, names))

    def test_organism_decides_the_names(self):
        self.assertIn("HBB-BS", sq.panel_genes("HEMOLYSIS", 10090))
        self.assertIn("HBB-BS", sq.panel_genes("HEMOLYSIS", None))     # organism unknown
        self.assertNotIn("HBB-BS", sq.panel_genes("HEMOLYSIS", 9606))  # human: symbols only
        self.assertIn("HBBL1", sq.panel_genes("HEMOLYSIS", 10116))
        self.assertNotIn("HBB-BH1", sq.panel_genes("HEMOLYSIS", 10090))  # embryonic, not lysis

    def test_loricrin_is_matched_by_its_current_name(self):
        self.assertIn("LORICRIN", sq.panel_genes("EPIDERMIS", 9606))


class MousePanels(Base):
    def test_mouse_haemoglobins_match_hemolysis(self):
        res, md = self.run_sq(MOUSE_ROWS, "--taxid", "10090")
        h = res["panels"]["HEMOLYSIS"]
        self.assertEqual(h["matched_genes"], ["CAT", "HBA", "HBB-BS"])
        self.assertEqual(h["n_panel_proteins"], 3)
        self.assertEqual(h["elevated_samples"], ["S1"])
        self.assertEqual(h["status"], "assessed")
        self.assertIn("markers: CAT, HBA, HBB-BS", md)

    def test_sidecar_taxid_is_used(self):
        meta = os.path.join(self.root, "search.fasta.meta.json")
        with open(meta, "w") as fh:
            json.dump({"organism": "Mus musculus", "taxid": 10090, "contaminant_set": "universal",
                       "n_contaminants_appended": 381, "contaminant_target_rule": "x"}, fh)
        res, _ = self.run_sq(MOUSE_ROWS, "--fasta-meta", meta)
        self.assertEqual(res["taxid"], 10090)
        self.assertIn("HBB-BS", res["panels"]["HEMOLYSIS"]["matched_genes"])


class ContaminantNotSample(Base):
    def test_cont_tagged_haemoglobin_is_listed_not_scored(self):
        rows = [["Cont_P02070", "HBB", "HBB_BOVIN"] + [26, 15, 15, 15, 15, 15],   # antibody prep
                ["P01942", "Hba", "HBA_MOUSE"] + FLAT, ["P24270", "Cat", "CATA_MOUSE"] + FLAT]
        res, md = self.run_sq(rows, "--taxid", "10090")
        h = res["panels"]["HEMOLYSIS"]
        self.assertEqual(h["contaminant_matches"], ["HBB (Cont_P02070)"])
        self.assertEqual(h["matched_genes"], ["CAT", "HBA"])     # the bovine HBB is not scored
        self.assertEqual(h["elevated_samples"], [])              # so S1 is not "haemolysed"
        self.assertTrue(any("contaminant, not sample" in f and "Cont_P02070" in f
                            for f in res["flags"]), res["flags"])
        self.assertIn("contaminant, not sample", md)

    def test_contaminants_run_de_removed_are_still_reported(self):
        with open(os.path.join(self.root, "de_provenance.json"), "w") as fh:
            json.dump({"contaminants": {"policy": "removed", "removed": True,
                                        "removed_table": "contaminants_removed.csv"}}, fh)
        with open(os.path.join(self.root, "contaminants_removed.csv"), "w") as fh:
            fh.write("Protein.Group,Genes,Contaminant.Group,Precursors,Contaminant.Precursors,"
                     "Removed.Entirely\nCont_P02070,HBB,TRUE,12,12,TRUE\n"
                     "P17182,Eno1,FALSE,30,1,FALSE\n")
        res, _ = self.run_sq(MOUSE_ROWS, "--taxid", "10090")
        self.assertEqual(res["panels"]["HEMOLYSIS"]["contaminant_matches"],
                         ["HBB (Cont_P02070) -- removed before DE"])


    def test_legacy_database_gets_a_caveat(self):
        # Silva08172026's sidecar: contaminants appended, no contaminant_target_rule, FASTA gone.
        meta = os.path.join(self.root, "search.fasta.meta.json")
        with open(meta, "w") as fh:
            json.dump({"organism": "Mus musculus", "taxid": 10090, "contaminant_set": "universal",
                       "n_contaminants_appended": 381,
                       "fasta": os.path.join(self.root, "gone.fasta")}, fh)
        rows = MOUSE_ROWS + [["Cont_P02070", "HBB", "HBB_BOVIN"] + FLAT]
        res, _ = self.run_sq(rows, "--fasta-meta", meta)
        f = next(x for x in res["flags"] if x.startswith("HEMOLYSIS: 1 `Cont_`"))
        self.assertIn("built by an older contaminant rule", f)

    def test_keratin_sample_keratins_are_not_called_not_sample(self):
        rows = [["Cont_P04264", "KRT1", "K2C1_HUMAN"] + FLAT, ["P15924", "DSP", "DESP_HUMAN"] + FLAT,
                ["P04040", "CAT", "CATA_HUMAN"] + FLAT]
        res, _ = self.run_sq(rows, "--taxid", "9606", "--keratin-sample")
        self.assertEqual(res["panels"]["EPIDERMIS"]["contaminant_matches"], [])
        self.assertFalse(any(f.startswith("EPIDERMIS") for f in res["flags"]), res["flags"])


class DidTheCheckRun(Base):
    def test_no_panel_matching_anything_says_the_check_could_not_run(self):
        # Accessions only, no gene names: the shape that read "not assessable" before.
        rows = [[f"Q{i:05d}", "", ""] + FLAT for i in range(30)]
        res, md = self.run_sq(rows, "--taxid", "10090")
        for name in sq.PANELS:
            self.assertEqual(res["panels"][name]["status"], "not_run", name)
        self.assertIn("check could not run (no panel gene matched this organism's symbols)", md)
        self.assertNotIn("not assessable", md)
        self.assertTrue(any("could not run" in f for f in res["flags"]))

    def test_one_empty_panel_on_a_matching_matrix_is_an_absence(self):
        rows = [["P04040", "CAT", "CATA_HUMAN"] + FLAT, ["P32119", "PRDX2", "PRDX2_HUMAN"] + FLAT,
                ["P15924", "DSP", "DESP_HUMAN"] + FLAT]
        res, md = self.run_sq(rows, "--taxid", "9606")
        m = res["panels"]["SKELETAL_MUSCLE"]
        self.assertEqual(m["status"], "none_detected")
        self.assertIn("the check ran", m["status_note"])
        self.assertFalse(any("could not run" in f for f in res["flags"]))

    def test_unverified_organism_is_not_called_an_absence(self):
        rows = [["P04040", "CAT", "CATA_X"] + FLAT]
        res, _ = self.run_sq(rows, "--taxid", "9031")      # chicken: names not verified
        self.assertEqual(res["panels"]["SKELETAL_MUSCLE"]["status"], "none_detected_unverified")
        self.assertFalse(res["panel_organisms_verified"])


class ConfoundTest(unittest.TestCase):
    """Silva08172026 (10 groups x 3) flagged all three panels "CONFOUNDED WITH GROUP": the old
    rule compared the highest and lowest of 10 noisy group means, whose spread grows with the
    number of groups. The permutation F-test asks whether the labels explain the score."""

    @staticmethod
    def design(values, per_group=3):
        samples = [f"S{i}" for i in range(len(values))]
        z = sq.zscore(dict(zip(samples, values)), samples)
        return z, {s: f"G{i // per_group}" for i, s in enumerate(samples)}, samples

    def test_ten_noise_groups_are_not_confounded(self):
        import random
        rng = random.Random(7)
        z, gmap, samples = self.design([rng.gauss(0, 1) for _ in range(30)])
        means = {}
        for s in samples:
            means.setdefault(gmap[s], []).append(z[s])
        gap = max(sum(v) / 3 for v in means.values()) - min(sum(v) / 3 for v in means.values())
        self.assertGreaterEqual(gap, 1.5)            # the old rule's trigger
        confounded, detail, p = sq.confound_check(z, gmap, samples, 1.5)
        self.assertFalse(confounded, detail)
        self.assertGreater(p, sq.CONFOUND_P)
        self.assertIn("permutation F-test across 10 groups", detail)

    def test_a_real_group_effect_is_flagged(self):
        import random
        rng = random.Random(3)
        vals = [rng.gauss(0, 1) + (2.5 if (i // 3) % 2 else 0) for i in range(30)]
        confounded, detail, p = sq.confound_check(*self.design(vals), 1.5)
        self.assertTrue(confounded, detail)
        self.assertLess(p, sq.CONFOUND_P)

    def test_three_vs_three_falls_back_to_separation(self):
        # 20 labellings: no test can reach 1%, so complete separation is reported instead.
        confounded, detail, p = sq.confound_check(*self.design([0.1, 0.2, 0.0, 2.1, 2.3, 2.2]), 1.5)
        self.assertIsNone(p)
        self.assertTrue(confounded)
        self.assertIn("too few samples for a test", detail)
        self.assertAlmostEqual(sq._min_p([3, 3]), 0.1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
