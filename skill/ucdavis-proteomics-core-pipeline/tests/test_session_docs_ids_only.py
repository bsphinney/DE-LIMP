#!/usr/bin/env python3
"""
The session README and AGENTS.md name an identifications-only search's one result table, and say
what a no-DE study with conditions cannot show (2.14.2 C2 review: biologist M2 and L5,
statistician M3).

Before, the README never mentioned identifications.tsv ("Search output (report.parquet, logs)"),
and AGENTS.md said the search results were "not in this folder" next to eleven DE rows also "not
in this folder" -- "which peptides were found?" had no file to answer it. Now, for a search
sage_lfq_check.search_quantified calls identifications only:
  - README "Start here" names identifications.tsv and its *_identifications.json record; the
    locations table and the folder row say identifications only;
  - AGENTS.md has "Identifications, per run" and "How the identifications were filtered" rows,
    and no precursor-report row;
and for any session with no DE: AGENTS.md folds the DE-only rows into one "Differential
expression | none in this analysis" line, neither document offers the DE as R, and with two or
more conditions both carry session_docs.no_de_condition_caveat (the analysis brief's words).
Hermetic: a temp session, in-process.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import sage_lfq_check as slc  # noqa: E402
import session_docs  # noqa: E402


def write(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(obj if isinstance(obj, str) else json.dumps(obj))


class Base(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.s = os.path.join(self._t.name, "2026-01-01_study")
        os.makedirs(os.path.join(self.s, "output", "tables"))
        os.makedirs(os.path.join(self.s, "logs"))
        write(os.path.join(self.s, "input", "conditions.csv"),
              "File.Name,Group\nrunA,ctrl\nrunB,treat\n")

    def tearDown(self):
        self._t.cleanup()

    def ids_only(self, folder="search_sage"):
        d = os.path.join(self.s, "output", folder)
        write(os.path.join(d, "identifications.tsv"), "run\tpeptide\nrunA\tPEPTIDEK\n")
        write(os.path.join(d, "search_provenance.json"),
              {"engine": "sage", "version": "0.14.7", "quantified": False})
        write(os.path.join(d, slc.RECORD), slc.not_applicable())
        write(os.path.join(d, "sage_identifications.json"), {"per_run": {}})

    def docs(self):
        f = session_docs.gather(self.s)
        return session_docs.readme_md(f), session_docs.agents_md(f)


class IdentificationsOnly(Base):
    def test_the_readme_names_the_identifications(self):
        self.ids_only()
        readme, _ = self.docs()
        start = readme[readme.index("## Start here"):readme.index("## Summary")]
        self.assertIn("`output/search_sage/identifications.tsv`", start)
        self.assertIn("`output/search_sage/sage_identifications.json`", start)
        self.assertIn("identifications only: no quantities", start)
        self.assertIn("- Search: identifications only (no quantities)", readme)
        self.assertIn("Search output (identifications.tsv: identifications only, logs)", readme)
        self.assertNotIn("Search output (report.parquet", readme)
        self.assertIn("| `output/search_sage` | the search output (identifications.tsv", readme)

    def test_agents_md_answers_which_peptides_were_found(self):
        self.ids_only()
        _, agents = self.docs()
        self.assertIn("| Identifications, per run (identifications only: no quantities) | "
                      "`output/search_sage/identifications.tsv` |", agents)
        self.assertIn("| How the identifications were filtered | "
                      "`output/search_sage/sage_identifications.json` |", agents)
        self.assertNotIn("Precursor-level search results", agents)

    def test_the_default_search_folder_never_names_a_report(self):
        """The E2E session's shape (E2E low list, F10): Sage --ids-only in output/search/ --
        "Search output (report.parquet, logs)" named a file the search never wrote."""
        self.ids_only(folder="search")
        readme, agents = self.docs()
        for doc in (readme, agents):
            self.assertNotIn("report.parquet", doc)
            self.assertIn("output/search/identifications.tsv", doc)

    def test_a_quantifying_search_keeps_its_report_row(self):
        """(pin)"""
        d = os.path.join(self.s, "output", "search")
        write(os.path.join(d, "report.parquet"), "parquet")
        write(os.path.join(d, "search_provenance.json"), {"engine": "diann", "version": "2.7.0"})
        readme, agents = self.docs()
        self.assertIn("| Precursor-level search results | `output/search/report.parquet` |",
                      agents)
        self.assertNotIn("Identifications, per run", agents)
        self.assertIn("Search output (report.parquet, logs)", readme)


class NoDE(Base):
    def test_de_rows_fold_into_one_line(self):
        self.ids_only()
        readme, agents = self.docs()
        self.assertIn("| Differential expression | none in this analysis |", agents)
        for row in ("What the DE did", "DE results, one table per contrast",
                    "Protein abundance per sample", "The DE as runnable R"):
            self.assertNotIn(row, agents)
        self.assertNotIn("The DE as R", agents)
        self.assertNotIn("The analysis, as code", readme)
        self.assertIn("- Differential expression: none in this analysis", readme)

    def test_conditions_without_a_comparison_are_said(self):
        self.ids_only()
        caveat = session_docs.no_de_condition_caveat({"ctrl": 1, "treat": 1})
        self.assertIn("Each condition has one run", caveat)
        readme, agents = self.docs()
        for doc in (readme, agents):
            self.assertIn(f"- **Conditions:** {caveat}", doc)

    def test_one_condition_or_a_de_has_no_such_caveat(self):
        self.assertIsNone(session_docs.no_de_condition_caveat({"all": 3}))
        self.assertIn("No differential expression was run",
                      session_docs.no_de_condition_caveat({"a": 2, "b": 2}))
        self.assertIn("not shown to be absent",
                      session_docs.no_de_condition_caveat({"a": 1, "b": 1}, "DDA"))


class DeliveredAgents(Base):
    """C2 verify (biologist L5, part): session.NO_DE_REASON ends with a task for staff ("If one
    was wanted, run step 8, then finalize again"), and the delivered AGENTS.md quoted it in
    "Missing parts". The client's copy now says NO_DE_CLIENT; the MANIFEST and the session's own
    AGENTS.md keep the staff's words."""

    def test_no_staff_task_in_the_delivered_missing_parts(self):
        import session
        from make_deposit import Manifest
        self.ids_only()
        p = session.paths_for(self.s)
        man = Manifest()
        session.de_only_step(p, man)
        man.write(p["manifest_txt"], "MANIFEST")
        f = session_docs.gather(self.s)
        delivered = session_docs.agents_md(f, for_delivery=True)
        missing = next(ln for ln in delivered.splitlines() if "**Missing parts**" in ln)
        self.assertNotIn("run step 8", missing)
        self.assertNotIn("finalize again", missing)
        self.assertEqual(missing.count(session_docs.NO_DE_CLIENT), len(session.DE_ONLY_PARTS))
        self.assertIn(session.NO_DE_REASON, session_docs.agents_md(f))     # staff's copy
        with open(p["manifest_txt"]) as fh:
            self.assertIn("run step 8", fh.read())                             # the MANIFEST


if __name__ == "__main__":
    unittest.main()
