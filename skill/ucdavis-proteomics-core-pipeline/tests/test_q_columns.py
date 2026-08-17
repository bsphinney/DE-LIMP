#!/usr/bin/env python3
"""
The q-value column definition must exist once, and the R mirror must not drift.

SKILL_OPEN_DEFECTS #2: the set was hand-written in five files with three
different answers. That is how run_de.R's dpc path drifted away from
build_maxlfq.R and applied a different identification FDR to the same report for
months (PR #48), and the same class of defect turned up again in the DE-LIMP app.

The definition now lives in scripts/diann_q_columns.py, with scripts/
diann_q_columns.R as a mirror -- because run_de.R and build_maxlfq.R are
standalone Rscripts that treat jsonlite as OPTIONAL, so reading a shared JSON
file would put a hard dependency in front of the identification filter.

A mirror is only safe if drift is impossible to miss. That is what this file is
for: it parses the R source and asserts the values equal the Python ones, so a
divergence is a CI failure rather than two subtly different FDR filters.
"""
import os
import re
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

from diann_q_columns import (  # noqa: E402
    FDR_REQUIRED, FDR_OPTIONAL, PROTEIN_Q_PREFERENCE,
    fdr_columns, protein_q_column,
)

R_FILE = os.path.join(SCRIPTS, "diann_q_columns.R")


def r_vector(name):
    """Pull a c("a", "b", ...) assignment out of the R source, spanning lines."""
    with open(R_FILE) as fh:
        src = fh.read()
    m = re.search(rf'^{re.escape(name)}\s*<-\s*c\((.*?)\)\s*$',
                  src, re.M | re.S)
    if not m:
        raise AssertionError(f"{name} not found in {R_FILE}")
    return re.findall(r'"([^"]+)"', m.group(1))


class TestRMirrorMatchesPython(unittest.TestCase):
    def test_fdr_required_matches(self):
        self.assertEqual(r_vector("DIANN_FDR_REQUIRED"), FDR_REQUIRED)

    def test_fdr_optional_matches(self):
        self.assertEqual(r_vector("DIANN_FDR_OPTIONAL"), FDR_OPTIONAL)

    def test_protein_preference_matches_including_ORDER(self):
        # Order is the whole point of a preference chain, so compare as a list.
        self.assertEqual(r_vector("DIANN_PROTEIN_Q_PREFERENCE"),
                         PROTEIN_Q_PREFERENCE)

    def test_r_functions_agree_with_python(self):
        """Run the R side for real -- parsing the literals is not enough if the
        accessor functions differ."""
        if not _have_rscript():
            self.skipTest("Rscript not available")
        available = ["Q.Value", "Lib.Q.Value", "Global.PG.Q.Value", "Junk"]
        r_expr = (
            f'source("{R_FILE}"); '
            f'a <- c({",".join(repr(c) for c in available)}); '
            'cat(paste(diann_fdr_columns(a), collapse=","), "|", '
            'diann_protein_q_column(a))'
        ).replace("'", '"')
        out = subprocess.run(["Rscript", "-e", r_expr],
                             capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stderr)
        r_fdr, r_pq = [x.strip() for x in out.stdout.split("|")]
        self.assertEqual(r_fdr.split(","), fdr_columns(available))
        self.assertEqual(r_pq, protein_q_column(available))


def _have_rscript():
    from shutil import which
    return which("Rscript") is not None


class TestFilterSetSemantics(unittest.TestCase):
    def test_fdr_columns_never_names_an_absent_column(self):
        # limpa::readDIANN() errors on an unknown name, and arrow silently
        # returns ZERO ROWS when filtering a column select() dropped -- the
        # defect that made MaxLFQ produce nothing in DE-LIMP v4.0.2-4.0.3.
        got = fdr_columns(["Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"])
        self.assertEqual(got, ["Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value"])
        self.assertNotIn("Global.PG.Q.Value", got)

    def test_all_six_when_the_report_has_them(self):
        every = FDR_REQUIRED + FDR_OPTIONAL
        self.assertEqual(fdr_columns(every), every)
        self.assertEqual(len(fdr_columns(None)), 6)

    def test_required_and_optional_do_not_overlap(self):
        self.assertEqual(set(FDR_REQUIRED) & set(FDR_OPTIONAL), set())

    def test_pg_q_value_is_in_the_filter_set(self):
        # Run-specific, but still applied: DE-LIMP's two --method values must
        # agree on FDR for the same report.
        self.assertIn("PG.Q.Value", FDR_REQUIRED + FDR_OPTIONAL)


class TestPreferenceSemantics(unittest.TestCase):
    def test_first_available_wins(self):
        self.assertEqual(protein_q_column(
            ["Q.Value", "Global.Q.Value", "Global.PG.Q.Value"]),
            "Global.PG.Q.Value")
        self.assertEqual(protein_q_column(["Q.Value", "Lib.PG.Q.Value"]),
                         "Lib.PG.Q.Value")

    def test_none_when_no_q_column_present(self):
        self.assertIsNone(protein_q_column(["Run", "Protein.Group"]))

    def test_preference_is_a_chain_not_a_filter_set(self):
        # Guard against a future "simplification" that merges the two concepts.
        # Applying the preference order as a filter would over-filter; filtering
        # on only the first available would under-filter.
        self.assertNotEqual(PROTEIN_Q_PREFERENCE, FDR_REQUIRED + FDR_OPTIONAL)
        self.assertEqual(len(PROTEIN_Q_PREFERENCE), 5)

    def test_every_preference_entry_is_a_real_diann_column(self):
        self.assertTrue(set(PROTEIN_Q_PREFERENCE)
                        <= set(FDR_REQUIRED + FDR_OPTIONAL))


class TestNoHandWrittenCopiesRemain(unittest.TestCase):
    """The point of the exercise: no call site may restate the set inline."""

    WIRED = {
        "run_de.R": "DIANN_FDR_REQUIRED",
        "build_maxlfq.R": "DIANN_FDR_REQUIRED",
        "compare_searches.py": "PROTEIN_Q_PREFERENCE",
        "run_search.py": "PROTEIN_Q_PREFERENCE",
    }

    def test_call_sites_reference_the_shared_definition(self):
        for fname, token in self.WIRED.items():
            with open(os.path.join(SCRIPTS, fname)) as fh:
                src = fh.read()
            self.assertIn(token, src, f"{fname} does not use the shared definition")

    def test_no_call_site_hardcodes_the_full_set(self):
        # An inline list of three-or-more q-columns in a call site means someone
        # started a fourth copy.
        pat = re.compile(r'"(?:Global\.)?(?:Lib\.)?(?:PG\.)?Q\.Value"')
        for fname in self.WIRED:
            path = os.path.join(SCRIPTS, fname)
            with open(path) as fh:
                lines = list(fh)
            for i, line in enumerate(lines, 1):
                if line.lstrip().startswith("#"):
                    continue
                if len(pat.findall(line)) >= 3:
                    self.fail(f"{fname}:{i} hardcodes a q-column set: {line.strip()}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
