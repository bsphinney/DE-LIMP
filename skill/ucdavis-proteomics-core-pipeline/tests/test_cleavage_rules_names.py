#!/usr/bin/env python3
"""The common enzymes have names in the Methods and the SDRF (2.14.2 C2, E2E F5).

A chymotrypsin Sage search (cleave_at 'FYWLM', restrict 'P') reached the Methods as "the cleavage
rule cleave_at='FYWLM', restrict='P' [enzyme name not mapped — confirm]" and the SDRF's cleavage
agent as TO-FILL: cleavage_rules.SAGE_CUT knew trypsin only. Now one table, keyed by the rule
itself (the residues cut at, in any order; the no-cut residue; C- or N-terminal), names the
common agents with their PSI-MS terms, for Sage's cleave_at/restrict/c_terminal and MSFragger's
cut/nocut/sense alike. What these tests pin:
  * Sage: chymotrypsin (FWYL in any order, and FWYLM; not before P), Lys-C, Lys-C/P, Arg-C,
    Asp-N and Lys-N (c_terminal false), Glu-C in bicarbonate (E) and phosphate (DE), trypsin;
  * semi-specific keeps the name; non-specific and no-digestion are as before;
  * MSFragger: FragPipe 24.0's own list rows for the same enzymes, and FragPipe's way to enter a
    rule its list lacks (name `nonspecific`, ENZYMATIC, the cuts given);
  * an unmapped rule keeps "[enzyme name not mapped — confirm]" in the Methods;
  * a Sage cfg read by make_methods.search_record names Chymotrypsin MS:1001306.
The residue sets: Sage DOCS.md v0.14.7 ("FWYL" for chymotrypsin), FragPipe 24.0's
enzyme-defs-supported.txt, Skyline 26.1 (skyline_peaks.SKYLINE_ENZYMES) and psi-ms.obo 4.2.14.
Stdlib only.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import cleavage_rules as cr  # noqa: E402

SRC = "test cfg"

# Sage enzyme block -> (PSI-MS name, accession)
SAGE_FORMS = [
    ({"cleave_at": "FWYL", "restrict": "P"}, ("Chymotrypsin", "MS:1001306")),
    ({"cleave_at": "FYWL", "restrict": "P"}, ("Chymotrypsin", "MS:1001306")),
    ({"cleave_at": "FYWLM", "restrict": "P"}, ("Chymotrypsin", "MS:1001306")),
    ({"cleave_at": "K", "restrict": "P"}, ("Lys-C", "MS:1001309")),
    ({"cleave_at": "K"}, ("Lys-C/P", "MS:1001310")),
    ({"cleave_at": "R", "restrict": "P"}, ("Arg-C", "MS:1001303")),
    ({"cleave_at": "D", "c_terminal": False}, ("Asp-N", "MS:1001304")),
    ({"cleave_at": "K", "c_terminal": False}, ("Lys-N", "MS:1003093")),
    ({"cleave_at": "E", "restrict": "P"}, ("V8-E", "MS:1001315")),
    ({"cleave_at": "DE", "restrict": "P"}, ("V8-DE", "MS:1001314")),
    ({"cleave_at": "ED", "restrict": "P"}, ("V8-DE", "MS:1001314")),
    ({"cleave_at": "KR", "restrict": "P"}, ("Trypsin", "MS:1001251")),
    ({"cleave_at": "RK"}, ("Trypsin/P", "MS:1001313")),
]

# FragPipe 24.0 enzyme-defs-supported.txt rows (name, cut, nocut, sense) -> name, accession
FRAGPIPE_ROWS = [
    (("chymotrypsin", "FLWY", "P", "C"), ("Chymotrypsin", "MS:1001306")),
    (("lysc", "K", "P", "C"), ("Lys-C", "MS:1001309")),
    (("lysc-p", "K", "", "C"), ("Lys-C/P", "MS:1001310")),
    (("argc", "R", "P", "C"), ("Arg-C", "MS:1001303")),
    (("aspn", "D", "", "N"), ("Asp-N", "MS:1001304")),
    (("lysn", "K", "", "N"), ("Lys-N", "MS:1003093")),
    (("gluc", "DE", "P", "C"), ("V8-DE", "MS:1001314")),
    (("gluc_bicarb", "E", "P", "C"), ("V8-E", "MS:1001315")),
    (("trypsin/chymotrypsin", "FYWLKR", "P", "C"), ("TrypChymo", "MS:1001312")),
    (("trypsin", "KR", "P", "C"), ("Trypsin", "MS:1001251")),
    (("stricttrypsin", "KR", "", "C"), ("Trypsin/P", "MS:1001313")),
]


def fragger(name, cut, nocut, sense, termini=2):
    spec = {0: "none", 1: "semi", 2: "full"}[termini]
    return {"specificity": spec, "termini": termini,
            "enzymes": [{"name": name, "cut": cut, "nocut": nocut, "sense": sense, "missed": 2}]}


class SageForms(unittest.TestCase):
    def test_each_common_enzyme_is_named(self):
        for enz, want in SAGE_FORMS:
            with self.subTest(enz=enz):
                cl = cr.sage_cleavage(dict(enz, missed_cleavages=2), SRC)
                self.assertEqual((cl["name"], cl["ac"]), want)
                self.assertEqual(cl["specificity"], "full")
                self.assertIn(f"cleave_at={enz['cleave_at']!r}", cl["rule"], "the rule stays")

    def test_semi_specific_keeps_the_name(self):
        cl = cr.sage_cleavage({"cleave_at": "FWYL", "restrict": "P", "semi_enzymatic": True},
                              SRC)
        self.assertEqual((cl["name"], cl["specificity"]), ("Chymotrypsin", "semi"))
        import make_methods as mm
        self.assertEqual(mm.cleavage_phrase(cl), "semi-specific Chymotrypsin cleavage")

    def test_non_specific_and_no_digestion_are_unchanged(self):
        """(pin)"""
        self.assertEqual(cr.sage_cleavage({"cleave_at": ""}, SRC)["name"], "unspecific cleavage")
        self.assertEqual(cr.sage_cleavage({"cleave_at": "$"}, SRC)["name"], "no cleavage")

    def test_a_rule_cutting_the_other_side_is_not_named_as_the_c_terminal_one(self):
        """(pin) Lys-C cuts after K; a rule cutting before K is Lys-N, never Lys-C."""
        self.assertEqual(cr.sage_cleavage({"cleave_at": "K", "restrict": "P",
                                           "c_terminal": False}, SRC)["name"], None)
        self.assertEqual(cr.sage_cleavage({"cleave_at": "D"}, SRC)["name"], None,
                         "C-terminal of D is not Asp-N")

    def test_an_unmapped_rule_keeps_its_tag(self):
        """(pin)"""
        import make_methods as mm
        cl = cr.sage_cleavage({"cleave_at": "FYW", "restrict": "P"}, SRC)
        self.assertIsNone(cl["name"])
        self.assertIn("[enzyme name not mapped — confirm]", mm.cleavage_phrase(cl))
        self.assertIn("cleave_at='FYW', restrict='P'", mm.cleavage_phrase(cl))


class FragPipeForms(unittest.TestCase):
    def test_fragpipe_24_list_rows_are_named(self):
        for (name, cut, nocut, sense), want in FRAGPIPE_ROWS:
            with self.subTest(enzyme=name):
                cl = cr.fragger_cleavage(fragger(name, cut, nocut, sense), SRC)
                self.assertEqual((cl["name"], cl["ac"]), want)

    def test_fragpipes_custom_rule_route_is_read_by_its_cuts(self):
        """FragPipe 24.0 (TabMsfragger.java): a rule not in its list is entered with the name
        `nonspecific` and the cleavage ENZYMATIC -- the cuts are the rule."""
        cl = cr.fragger_cleavage(fragger("nonspecific", "FWYLM", "P", "C"), SRC)
        self.assertEqual((cl["name"], cl["specificity"]), ("Chymotrypsin", "full"))

    def test_a_non_enzymatic_search_is_still_non_specific(self):
        """(pin)"""
        cl = cr.fragger_cleavage(fragger("nonspecific", "@", "", "C", termini=0), SRC)
        self.assertEqual((cl["name"], cl["specificity"]), ("unspecific cleavage", "none"))

    def test_two_enzymes_are_not_named_as_one(self):
        """(pin)"""
        dg = fragger("trypsin", "KR", "P", "C")
        dg["enzymes"].append({"name": "lysc", "cut": "K", "nocut": "P", "sense": "C",
                              "missed": 2})
        self.assertIsNone(cr.fragger_cleavage(dg, SRC)["name"])


class MethodsAndSdrfRead(unittest.TestCase):
    def test_a_chymotrypsin_sage_cfg_names_the_agent(self):
        import make_methods as mm
        with tempfile.TemporaryDirectory() as d:
            cfg = os.path.join(d, "sage_config.json")
            with open(cfg, "w") as fh:
                json.dump({"database": {"enzyme": {"missed_cleavages": 2, "min_len": 7,
                                                   "max_len": 30, "cleave_at": "FYWLM",
                                                   "restrict": "P"},
                                        "static_mods": {"C": 57.0215}}, "quant": {"lfq": False}},
                          fh)
            rec = mm.search_record(params=cfg)
        cl = rec["cleavage"]
        self.assertEqual((cl["name"], cl["ac"]), ("Chymotrypsin", "MS:1001306"))
        self.assertEqual(mm.cleavage_phrase(cl), "Chymotrypsin specificity")
        self.assertNotIn("not mapped", mm.cleavage_phrase(cl))


if __name__ == "__main__":
    unittest.main()
