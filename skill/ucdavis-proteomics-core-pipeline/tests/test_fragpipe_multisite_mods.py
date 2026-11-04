#!/usr/bin/env python3
"""
A FragPipe variable modification on several sites is read site by site (2.14.2 C2 review:
proteomics L3).

fragger.params' "variable_mod_NN = -17.0265 nQnC 1" (FragPipe's pyro-glu line, off in the shipped
DIA presets) was read as one target string, so the Methods said "Gln->pyro-Glu (Any N-term QnC)".
fragpipe_params.site_parts splits a site at each terminus marker (n, c, [, ]) into position and
residues. Hermetic: in-process.
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import fragpipe_params as fpp  # noqa: E402


class SiteParts(unittest.TestCase):
    def test_a_multi_site_mod_is_split(self):
        self.assertEqual(fpp.mods({"variable_mod_01": "-17.0265 nQnC 1"}),
                         [("variable", -17.0265, "Any N-term", "QC", 1)])
        self.assertEqual(fpp.site_parts("STYn^"), [("Anywhere", "STY"), ("Any N-term", "")])

    def test_single_sites_read_as_before(self):
        """(pin) the nonspecific-peptidome fixture's nQ, nE and c^, a residue set, protein N-term"""
        for site, want in (("nQ", [("Any N-term", "Q")]), ("nE", [("Any N-term", "E")]),
                           ("c^", [("Any C-term", "")]), ("[^", [("Protein N-term", "")]),
                           ("STY", [("Anywhere", "STY")]), ("M", [("Anywhere", "M")]),
                           ("*", [("Anywhere", "")])):
            with self.subTest(site=site):
                self.assertEqual(fpp.site_parts(site), want)


if __name__ == "__main__":
    unittest.main()
