#!/usr/bin/env python3
"""The biotinylated-carboxylase reference panel also finds bacterial biotin proteins (2.14.2, I26).

The panel (experiment_type.CARBOXYLASES) is mammalian gene symbols: PC, PCCA, MCCC1, ACACA,
ACACB. A V. cholerae TurboID study (a Core submission) has them under locus tags (Genes = VC_0296), so
NORMALIZATION_CHECK.md said "none found" for both options while the DE held the biotin carboxyl
carrier protein, the biotin carboxylase and the oxaloacetate decarboxylase alpha subunit,
quantified in every sample. Now the panel also matches the search FASTA's protein descriptions
against experiment_type.CARBOXYLASE_NAMES (specific biotin-enzyme names, never a bare
"carboxylase", which would also take decarboxylases and RuBisCO); when nothing is found the summary
says what the panel looks for. Stdlib only; hermetic (synthetic matrix and FASTA).
"""
import csv
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402

SAMPLES = ["T1", "T2", "T3", "C1", "C2", "C3"]
# (accession, gene, description, log2 values): three biotin proteins under locus tags, and three
# that only LOOK like carboxylases by a bare substring
BACTERIA = [
    ("Q9KV19", "VC_0296", "Biotin carboxyl carrier protein of acetyl-CoA carboxylase",
     [20.1, 20.3, 20.2, 21.5, 21.6, 21.4]),
    ("Q9KV20", "VC_0295", "Biotin carboxylase", [19.0, 19.2, 19.1, 20.5, 20.6, 20.4]),
    ("Q9KUX1", "VC_0550", "Oxaloacetate decarboxylase alpha subunit",
     [18.0, 18.1, 18.2, 19.0, 19.1, 19.2]),
    ("Q9KQ01", "VC_1111", "Glutamate decarboxylase", [17.0, 17.5, 17.2, 17.1, 17.3, 17.4]),
    ("Q9KQ02", "VC_2222", "Ribulose bisphosphate carboxylase large chain",
     [16.0, 16.4, 16.2, 16.1, 16.3, 16.5]),
    ("Q9KQ03", "VC_3333", "Phosphoenolpyruvate carboxylase", [15.0, 15.2, 15.4, 15.1, 15.3, 15.5]),
]


def write_de(de_dir, proteins):
    os.makedirs(de_dir, exist_ok=True)
    with open(os.path.join(de_dir, "Expression_Matrix.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein.Group", "Genes", "Protein.Names"] + SAMPLES)
        for acc, gene, _, vals in proteins:
            w.writerow([acc, gene, f"{acc}_VIBCH"] + vals)


def write_fasta(path, proteins):
    with open(path, "w") as fh:
        for acc, gene, desc, _ in proteins:
            fh.write(f">tr|{acc}|{acc}_VIBCH {desc} OS=Vibrio cholerae OX=243277 "
                     f"OrderedLocusNames={gene} PE=3 SV=1\nMSTNPEPTIDEK\n")


class ThePanelByName(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.de = os.path.join(self.d, "de")
        write_de(self.de, BACTERIA)
        self.fasta = os.path.join(self.d, "search.fasta")
        write_fasta(self.fasta, BACTERIA)

    def test_a_bacterial_matrix_is_found_by_name(self):
        desc = nc.fasta_descriptions(self.fasta, {p[0] for p in BACTERIA})
        cx = nc.carboxylase_stability(self.de, desc)
        self.assertEqual(cx["n_found"], 3)
        self.assertEqual(sorted(cx["sd_log2"]), ["VC_0295", "VC_0296", "VC_0550"])
        self.assertEqual(cx["matched_by"]["VC_0296"], "biotin carboxyl carrier protein")

    def test_a_bare_carboxylase_is_not_the_panel(self):
        for acc, _, d, _ in BACTERIA[3:]:
            with self.subTest(d=d):
                self.assertIsNone(et.carboxylase_of("", d))

    def test_the_gene_symbols_still_match(self):
        """The mammalian symbols match as before (the carboxylase_stability half is a pin)."""
        self.assertEqual(et.carboxylase_of("PC;PCX", None), "PC")
        write_de(self.de, [("P11498", "PC", "", [20, 20.1, 20.2, 20.3, 20.4, 20.5])])
        self.assertEqual(nc.carboxylase_stability(self.de)["n_found"], 1)

    def test_the_session_fasta_is_read_by_run(self):
        """run() reads the descriptions from the session's FASTA record (input/search.fasta.meta.json
        `fasta`), so a proximity check finds the bacterial panel."""
        s = os.path.join(self.d, "session")
        os.makedirs(os.path.join(s, "input"))
        with open(os.path.join(s, et.FASTA_META), "w") as fh:
            json.dump({"fasta": self.fasta}, fh)
        self.assertEqual(nc.session_fasta(s), self.fasta)
        self.assertIsNone(nc.session_fasta(os.path.join(self.d, "nowhere")))
        # the whole check, on a proximity session: the record carries the panel found by name
        et.record(s, "proximity", "user", controls="C")
        cond = os.path.join(s, "input", "conditions.csv")
        with open(cond, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["File.Name", "Group"])
            w.writerows((x, x[0]) for x in SAMPLES)
        report = os.path.join(self.d, "report.parquet")
        with open(report, "wb") as fh:
            fh.write(b"not a parquet: B1-B3 are not assessed, the panel still is")
        dirs = {}
        for q in (et.NORMALISED, et.RAW):
            dirs[q] = os.path.join(self.d, "nc", q)
            write_de(dirs[q], BACTERIA)
            with open(os.path.join(dirs[q], "de_provenance.json"), "w") as fh:
                json.dump({"method": "dpc", "contrasts": ["T-C"], "de_tables": {}}, fh)
        rec, _ = nc.run(report, cond, dirs[et.NORMALISED], dirs[et.RAW], session=s,
                        out=os.path.join(self.d, "nc"))
        cx = rec["carboxylases"]
        self.assertEqual(cx["fasta"], self.fasta)
        for q in (et.NORMALISED, et.RAW):
            self.assertEqual(cx[q]["n_found"], 3, cx[q])

    def test_none_found_says_what_the_panel_looks_for(self):
        md = nc.summary_md({"default": {"quantities": et.NORMALISED},
                            "carboxylases": {"normalised": {"n_found": 0}, "raw": {"n_found": 0},
                                             "genes": list(et.CARBOXYLASES)}})
        line = [ln for ln in md.splitlines() if ln.startswith("- normalised:")][0]
        self.assertIn("none found", line)
        self.assertIn("mammalian gene symbols", line)
        self.assertIn("biotin carboxyl carrier protein", line)


# A human matrix: the five panel proteins and the two beta chains of the same enzymes, with their
# UniProt descriptions. The beta chains carry no biotin (UniProt KW-0092 is on the alpha chains).
HUMAN = [
    ("P11498", "PC", "Pyruvate carboxylase, mitochondrial", [20.0, 20.2, 20.1, 20.3, 20.4, 20.2]),
    ("P05165", "PCCA", "Propionyl-CoA carboxylase alpha chain, mitochondrial",
     [19.0, 19.3, 19.1, 19.2, 19.4, 19.0]),
    ("P05166", "PCCB", "Propionyl-CoA carboxylase beta chain, mitochondrial",
     [18.0, 19.5, 18.2, 19.9, 18.4, 19.1]),
    ("Q96RQ3", "MCCC1", "Methylcrotonoyl-CoA carboxylase subunit alpha, mitochondrial",
     [17.0, 17.1, 17.3, 17.2, 17.0, 17.4]),
    ("Q9HCC0", "MCCC2", "Methylcrotonoyl-CoA carboxylase beta chain, mitochondrial",
     [16.0, 17.8, 16.1, 17.5, 16.3, 17.9]),
    ("Q13085", "ACACA", "Acetyl-CoA carboxylase 1", [21.0, 21.1, 21.2, 21.0, 21.3, 21.1]),
    ("O00763", "ACACB", "Acetyl-CoA carboxylase 2", [15.0, 15.2, 15.1, 15.3, 15.0, 15.2]),
]


class NoBetaChains(unittest.TestCase):
    """Review of 2.14.2: the name fallback put PCCB and MCCC2 -- the non-biotinylated beta chains
    -- into a human session's panel, against its KW-0092 definition, and moved the median SD from
    2.14.1's. Descriptions naming a beta chain/subunit or a carboxyltransferase are not the panel."""

    def test_pccb_and_mccc2_are_excluded(self):
        for gene, desc in (("PCCB", HUMAN[2][2]), ("MCCC2", HUMAN[4][2]),
                           ("pccB", "propionyl-CoA carboxylase carboxyltransferase subunit"),
                           ("accD", "acetyl-CoA carboxylase, carboxyl transferase subunit beta"),
                           ("mccB", "methylcrotonoyl-CoA carboxylase subunit beta")):
            with self.subTest(gene=gene):
                self.assertIsNone(et.carboxylase_of(gene, desc))

    def test_the_mammalian_panel_is_unchanged_from_2_14_1(self):
        with tempfile.TemporaryDirectory() as d:
            write_de(d, HUMAN)
            by_symbol = nc.carboxylase_stability(d)          # 2.14.1: the symbols only
            with_names = nc.carboxylase_stability(d, {p[0]: p[2] for p in HUMAN})
        self.assertEqual(sorted(with_names["sd_log2"]), sorted(et.CARBOXYLASES))
        for k in ("sd_log2", "median_sd_log2", "n_found"):
            self.assertEqual(with_names[k], by_symbol[k], k)

    def test_ncbi_pgaps_alpha_subunit_name_is_matched(self):
        self.assertEqual(et.carboxylase_of("VC_0550", "oxaloacetate decarboxylase subunit alpha"),
                         "oxaloacetate decarboxylase subunit alpha")


if __name__ == "__main__":
    unittest.main()
