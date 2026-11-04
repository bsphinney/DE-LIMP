#!/usr/bin/env python3
"""make_methods.UNIMOD's variable modifications come from estimate_params.VAR_MODS (C1 review,
F5b; rule 3).

make_methods.UNIMOD wrote the Oxidation and Deamidated masses out a second time beside
estimate_params.VAR_MODS, the table the search is built from -- as it once did the cysteine
alkylations, which it reads from estimate_params.CYS_REAGENTS. VAR_MODS now carries each
modification's Unimod record (unimod, mod, mass, sites) and builds its DIA-NN and Sage forms
from it; make_methods builds those UNIMOD entries from it. The Methods text is unchanged (pin).

Hermetic: estimate_params.py as a subprocess on a temp dir; make_methods.py in-process.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import estimate_params as ep  # noqa: E402
import make_methods as mm  # noqa: E402

# make_methods.UNIMOD as it was before (2.14.2 at ba6dc244)
BEFORE = {1: ("Acetyl", 42.010565), 4: ("Carbamidomethyl", 57.021464),
          7: ("Deamidated", 0.984016), 21: ("Phospho", 79.966331),
          28: ("Gln->pyro-Glu", -17.026549), 35: ("Oxidation", 15.994915),
          108: ("Nethylmaleimide", 125.047679), 6: ("Carboxymethyl", 58.005479),
          24: ("Propionamide", 71.037114), 39: ("Methylthio", 45.987721),
          31: ("Pyridylethyl", 105.057849)}


class OneTable(unittest.TestCase):
    def test_the_variable_mods_entries_are_built_from_var_mods(self):
        for word in ("ox", "deam"):
            with self.subTest(word=word):
                vm = ep.VAR_MODS[word]
                self.assertEqual(mm.UNIMOD[vm["unimod"]], (vm["mod"], vm["mass"]))
        with open(os.path.join(SCRIPTS, "make_methods.py"), encoding="utf-8") as fh:
            src = fh.read()
        for mass in ("15.994915", "0.984016"):
            self.assertNotIn(mass, src, f"{mass} is VAR_MODS', not make_methods' to write")

    def test_var_mods_engine_forms_are_unchanged(self):
        """(pin) the DIA-NN and Sage forms the search is built with."""
        self.assertEqual(ep.VAR_MODS["ox"]["diann"], "UniMod:35,15.994915,M")
        self.assertEqual(ep.VAR_MODS["ox"]["sage"], {"M": [15.994915]})
        self.assertEqual(ep.VAR_MODS["deam"]["diann"], "UniMod:7,0.984016,NQ")
        self.assertEqual(ep.VAR_MODS["deam"]["sage"], {"N": [0.984016], "Q": [0.984016]})

    def test_unimod_is_the_same_table(self):
        """(pin)"""
        self.assertEqual(mm.UNIMOD, BEFORE)


class SameText(unittest.TestCase):
    """(pin) the Methods' modification sentences, as before."""

    def sentences(self, engine, acquisition, out):
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, out)
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                                "--engine", engine, "--acquisition", acquisition, "--instrument",
                                "timsTOF HT", "--precursor-mz-range", "299.5", "1200.5",
                                "--var-mods", "ox,deam", "--out", cfg],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            prov = os.path.join(d, "search_provenance.json")
            with open(prov, "w") as fh:
                json.dump({"engine": engine, "version": "x", "params_file": cfg}, fh)
            text = mm.search_paragraph(mm.search_record(search_prov=prov))
        return [s for s in re.split(r"(?<=\.) ", text) if "modif" in s.lower()]

    def test_diann(self):
        self.assertEqual(self.sentences("diann", "DIA", "params.cfg"),
                         ["Fixed modifications: Carbamidomethyl (C).",
                          "Variable modifications: Oxidation (M); Deamidated (NQ) (at most 2 per "
                          "peptide)."])

    def test_sage(self):
        self.assertEqual(self.sentences("sage", "DDA", "sage.json"),
                         ["Fixed modifications: Carbamidomethyl (C).",
                          "Variable modifications: Oxidation (M); Deamidated (N); Deamidated (Q); "
                          "Acetyl (Protein N-term) (at most 2 per peptide)."])


if __name__ == "__main__":
    unittest.main()
