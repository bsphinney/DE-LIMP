#!/usr/bin/env python3
"""
A Sage identifications-only search says that no protein-level FDR filter was applied, and how
many protein groups pass Sage's own (2.14.2 C2 review: statistician H1 + proteomics L1).

adapt_sage_ids keeps rank-1 target PSMs at spectrum_q and peptide_q <= 0.01, and filters nothing
on protein_q. On the two-run DDA test the kept PSMs map to 198 protein groups, of which 143 pass
Sage's 1% picked-protein FDR, yet the Methods gave only the PSM and peptide rule: a reader counting
proteins in identifications.tsv got ~198 "at 1% FDR". Now:
  - sage_identifications.json `methods_sentence` (which make_methods quotes) says no protein-level
    filter was applied and gives N of M protein groups at protein_q <= 0.01, counted from the
    table (sage_ids_protein_groups); the record's totals carry both numbers;
  - its rule reads "over all runs together, listed per run" (the q-values are global);
  - the OUTPUT_FILES catalogue says the same, and each engine's row unit.
identifications.tsv already carries a per-row protein_q (SAGE_IDS_COLUMNS).
Hermetic: in-process, on a temp dir; the end-to-end part needs pyarrow.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_methods  # noqa: E402
import make_report  # noqa: E402
import run_search  # noqa: E402
import test_sage_dotd_ids as fx  # noqa: E402

PROTEIN_WORDS = "No protein-level FDR filter was applied"


def row(proteins, protein_q, run="r1", peptide="PEPTIDEK"):
    return (run, peptide, proteins, 2, 0.001, 0.001, protein_q)


class TheCount(unittest.TestCase):
    def test_groups_and_those_passing_protein_q(self):
        rows = [row("sp|A", 0.001), row("sp|A", 0.002, run="r2"), row("sp|B;sp|C", 0.04),
                row("sp|D", None), row("sp|E", 0.01)]
        self.assertEqual(run_search.sage_ids_protein_groups(rows), (4, 2))

    def test_the_methods_sentence_says_no_protein_filter_with_the_count(self):
        s = run_search.sage_ids_methods_sentence(198, 143)
        self.assertTrue(s.startswith(run_search.SAGE_IDS_METHODS), s)
        self.assertIn(PROTEIN_WORDS, s)
        self.assertIn("143 of the 198 protein groups", s)
        self.assertIn("over all runs together, listed per run", s)
        self.assertNotIn("Identifications in each run were", s)

    def test_the_catalogue_says_it(self):
        text = " ".join(desc for pat, _, desc in make_report.CATALOG
                        if "identifications" in pat)
        self.assertIn("No protein-level FDR filter was applied", text)
        self.assertIn("protein_q", text)
        self.assertIn("one row per run, modified peptide and charge", text)


@unittest.skipUnless(fx.HAVE_ARROW, "pyarrow is needed to write the Sage fixtures")
class TheRecordAndTheMethods(unittest.TestCase):
    def test_adapt_records_the_count_and_the_methods_quote_it(self):
        with tempfile.TemporaryDirectory() as out:
            fx.write_ids_results(out, ["runA.d", "runB.d"])
            with open(os.path.join(out, "search_provenance.json"), "w") as fh:
                json.dump({"engine": "sage", "version": "0.14.7", "quantified": False}, fh)
            run_search.adapt_sage_ids(out)
            with open(os.path.join(out, run_search.SAGE_IDS_RECORD)) as fh:
                rec = json.load(fh)
            # the fixture's kept rows map to AAA, BBB;CCC and FFF, every one at protein_q <= 0.01
            self.assertEqual((rec["totals"]["protein_groups"],
                              rec["totals"]["protein_groups_protein_q"]), (3, 3))
            self.assertIn(PROTEIN_WORDS, rec["methods_sentence"])
            self.assertIn("3 of the 3 protein groups", rec["methods_sentence"])
            self.assertIn("none:", rec["protein_filter"])
            cfg = os.path.join(out, "sage_config.json")
            with open(cfg, "w") as fh:
                json.dump({"database": {"enzyme": {"cleave_at": "KR", "restrict": "P"}},
                           "quant": {"lfq": False}}, fh)
            mrec = make_methods.search_record(params=cfg, search_prov=os.path.join(
                out, "search_provenance.json"))
            self.assertIn(PROTEIN_WORDS, make_methods.search_paragraph(mrec))


if __name__ == "__main__":
    unittest.main()
