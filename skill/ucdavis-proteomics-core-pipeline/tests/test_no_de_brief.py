#!/usr/bin/env python3
"""
The analysis brief of a session with no differential expression is a brief for what it has
(2.14.2 C2 review: biologist H1, statistician M3, biologist M4, biologist L6).

Before, only the title and first sentence changed: the body still demanded volcanoes, a
significance rule, "Key Findings Per Comparison", a normalisation-check note and DE Methods, and
pointed the writer at tables/Expression_Matrix.csv, methods.txt and de_provenance.json -- files a
no-DE session does not have -- while never naming identifications.tsv or the search report.
Now, with no DE record and no DE table:
  - the body names the session's search files (session.search_folder): identifications.tsv and
    its *_identifications.json record for an identifications-only search, report.parquet for a
    quantifying one, and search_provenance.json; it says "identifications only: no quantities"
    when the search record says so (sage_lfq_check.search_quantified, the one predicate);
  - it makes no DE demand, and carries the first-pass caveat in its Methods form;
  - with two or more conditions it says differences between them are not tested (one run per
    condition: not separable from run-to-run variation) (statistician M3).
Every brief names the engine and version by session.search_engine (through the page's reader),
not the workflow manifest or --engine (biologist M4); a most-variable heatmap goes to QC (L6).
Hermetic: analysis_prompt.py on temp sessions.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import dda_first_pass  # noqa: E402
import run_search  # noqa: E402
import sage_lfq_check  # noqa: E402

DE_DEMANDS = ("Significance rule", "Key Findings Per Comparison", "Cross-Comparison Biomarkers",
              "High-Confidence Findings", "tables/Expression_Matrix.csv", "methods.txt",
              "de_provenance.json", "Normalisation check", "State clearly whether the data are "
              "fit for differential analysis", "adj.P.Val", "label-free quant")


def flat(text):
    return " ".join(text.split())


def write(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(obj if isinstance(obj, str) else json.dumps(obj))
    return path


class Base(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.s = os.path.join(self._t.name, "session")
        os.makedirs(os.path.join(self.s, "output", "tables"))
        self.conditions("runA,ctrl", "runB,treat")

    def tearDown(self):
        self._t.cleanup()

    def conditions(self, *rows):
        self.cond = write(os.path.join(self.s, "input", "conditions.csv"),
                          "File.Name,Group\n" + "\n".join(rows) + "\n")

    def ids_only_search(self, folder="search_sage"):
        d = os.path.join(self.s, "output", folder)
        write(os.path.join(d, "identifications.tsv"),
              "\t".join(run_search.SAGE_IDS_COLUMNS) + "\nrunA\tPEPTIDEK\tsp|P1|A\t2\t0.001\t"
              "0.001\t0.001\n")
        write(os.path.join(d, "sage_identifications.json"),
              {"methods_sentence": run_search.sage_ids_methods_sentence(198, 143),
               "row_unit": "one row per (run, modified peptide, charge)"})
        write(os.path.join(d, sage_lfq_check.RECORD), sage_lfq_check.not_applicable())
        write(os.path.join(d, "search_provenance.json"),
              {"engine": "sage", "version": "0.14.7", "quantified": False})
        return d

    def brief(self, *extra, de=None):
        if de:
            for name, obj in de.items():
                write(os.path.join(self.s, "output", "tables", name), obj)
        out = os.path.join(self._t.name, "ANALYSIS_PROMPT.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                            "--de-dir", os.path.join(self.s, "output", "tables"), "--out", out,
                            "--conditions", self.cond, "--session", self.s] + list(extra),
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out, encoding="utf-8") as fh:
            return flat(fh.read())


class TheNoDEBrief(Base):
    def test_an_ids_only_search_is_briefed_from_its_own_files(self):
        self.ids_only_search()
        b = self.brief("--acquisition", "DDA")
        self.assertIn("Identifications only: no quantities", b)
        for name in ("identifications.tsv", "sage_identifications.json",
                     "search_provenance.json"):
            self.assertIn(f"`output/search_sage/{name}`", b)
        self.assertIn("143 of the 198 protein groups", b)
        for demand in DE_DEMANDS:
            self.assertNotIn(demand, b)
        for section in ("### Identification depth per run", "### Notable proteins",
                        "### QC Assessment", "### Data Quality Notes", "### Next steps"):
            self.assertIn(section, b)

    def test_one_run_per_condition_gets_the_caveat(self):
        self.ids_only_search()
        b = self.brief("--acquisition", "DDA")
        self.assertIn("Each condition has one run, so differences between conditions", b)
        self.assertIn("do not present them as differences between conditions", b)
        self.assertIn("a protein missing from one run is not shown to be absent", b)

    def test_replicated_conditions_without_a_de_get_it_too(self):
        self.conditions("a1,ctrl", "a2,ctrl", "b1,treat", "b2,treat")
        self.ids_only_search()
        b = self.brief("--acquisition", "DIA")
        self.assertIn("No differential expression was run, so differences between conditions", b)
        self.assertNotIn("Each condition has one run", b)

    def test_one_condition_has_no_comparison_to_warn_about(self):
        self.conditions("a1,all", "a2,all")
        self.ids_only_search()
        self.assertNotIn("Conditions without a comparison", self.brief())

    def test_a_quantifying_search_names_its_report(self):
        d = os.path.join(self.s, "output", "search")
        write(os.path.join(d, "report.parquet"), "parquet")
        write(os.path.join(d, "report.stats.tsv"), "File.Name\tPrecursors.Identified\n")
        write(os.path.join(d, "search_provenance.json"), {"engine": "diann", "version": "2.7.0"})
        b = self.brief("--acquisition", "DIA")
        self.assertIn("`output/search/report.parquet`", b)
        self.assertIn("no differential expression was run", b)
        self.assertNotIn("Identifications only", b)
        for demand in DE_DEMANDS:
            self.assertNotIn(demand, b)

    def test_a_first_pass_report_carries_its_methods_caveat(self):
        d = os.path.join(self.s, "output", "search_dda")
        write(os.path.join(d, "report.parquet"), "parquet")
        prov = {"engine": "diann", "version": "2.7.0", "report_used": "first_pass",
                "dda_first_pass": {"report_used": "first_pass", "pg_maxlfq": "absent"}}
        write(os.path.join(d, "search_provenance.json"), prov)
        b = self.brief("--acquisition", "DDA")
        self.assertIn(flat(dda_first_pass.caution(prov, methods=True)), b)
        self.assertNotIn("DPC-Quant", b)


class TheE2EShape(Base):
    """The two-run timsTOF DDA end-to-end session (2.14.2 E2E F1): Sage --ids-only into the
    default output/search/, identifications.tsv and no report.parquet, output/tables/ and
    output/figures/ empty (step 8b is the DE's), the agent passing --report and --figures-dir as
    for a DE session. The brief named Expression_Matrix.csv, methods.txt and de_provenance.json,
    none of which existed, and never identifications.tsv."""

    def test_the_brief_names_what_the_session_has(self):
        d = self.ids_only_search(folder="search")
        figs = os.path.join(self.s, "output", "figures")
        os.makedirs(figs)
        b = self.brief("--acquisition", "DDA", "--engine", "sage", "--figures-dir", figs,
                       "--report", os.path.join(d, "report.parquet"))
        self.assertIn("`output/search/identifications.tsv`", b)
        self.assertIn("`output/search/sage_identifications.json`", b)
        self.assertIn("Identifications only: no quantities", b)
        for absent in DE_DEMANDS + ("report.parquet", "Expression_Matrix", "MaxLFQ",
                                    "figures.json"):
            self.assertNotIn(absent, b)
        self.assertIn("give its evidence from the identifications, per run", b)


class TheVerifyLeftovers(Base):
    """C2 verify (biologist, fixer A): with no search chosen the brief said "The search
    quantified precursors" (the predicate defaulted to True with no folder); and for a
    quantifying search it pointed at "the search's recorded FDR (above)", a rule it never stated."""

    def test_no_search_chosen_says_nothing_about_quantities(self):
        self.ids_only_search(folder="search_a")
        write(os.path.join(self.s, "output", "search_b", "search_provenance.json"),
              {"engine": "diann", "version": "2.7.0"})
        b = self.brief("--report", os.path.join(self.s, "output", "search_b", "report.parquet"))
        self.assertIn("Which search: not recorded", b)
        self.assertNotIn("The search quantified precursors", b)
        self.assertNotIn("Identifications only: no quantities", b)
        self.assertIn("Whether the search quantified anything is not recorded here", b)
        self.assertNotIn("search_b/report.parquet", b)      # never --report's folder instead

    def quantifying(self, cfg=None):
        d = os.path.join(self.s, "output", "search")
        write(os.path.join(d, "report.parquet"), "parquet")
        write(os.path.join(d, "search_provenance.json"), {"engine": "diann", "version": "2.7.0"})
        if cfg is not None:
            write(os.path.join(d, "params.resolved.cfg"), cfg)
        return self.brief("--acquisition", "DIA")

    def test_a_quantifying_search_states_its_fdr_rule(self):
        b = self.quantifying("--qvalue 0.01\n--threads 8\n")
        self.assertIn("- Identification rule: Precursor identifications were filtered at 1% FDR "
                      "(q ≤ 0.01).", b)
        self.assertIn("at the search's recorded FDR (the identification rule above)", b)

    def test_an_unrecorded_rule_is_said_never_referred_to(self):
        b = self.quantifying()
        self.assertIn("- Identification rule: **not recorded** in the search record", b)
        self.assertIn("(the FDR rule they rest on is not recorded: say so)", b)
        self.assertNotIn("(above)", b)


class TheEngine(Base):
    DE = {"DE_B_vs_A.csv": "Protein.Group,logFC,adj.P.Val\nP1,1.0,0.01\n",
          "de_provenance.json": {"contrasts": ["B-A"], "adjp": 0.05}}

    def intro(self, b):
        return b[b.index("(engine: "):b.index("acquisition:")]

    def test_the_search_record_wins_over_engine_and_the_manifest(self):
        write(os.path.join(self.s, "output", "search", "search_provenance.json"),
              {"engine": "diann", "version": "2.7.0"})
        write(os.path.join(self.s, "input", "wf", "workflow.manifest.json"),
              {"engine": {"name": "sage", "version": "0.14.7"}})
        self.assertIn("`diann 2.7.0`", self.intro(self.brief("--engine", "sage", de=self.DE)))

    def test_a_null_version_stays_unrecorded(self):
        """a conda Sage records no version: never the manifest's pin in its place"""
        write(os.path.join(self.s, "output", "search_sage", "search_provenance.json"),
              {"engine": "sage", "version": None})
        write(os.path.join(self.s, "input", "wf", "workflow.manifest.json"),
              {"engine": {"name": "sage", "version": "0.14.7"}})
        intro = self.intro(self.brief(de=self.DE))
        self.assertIn("sage (version not recorded)", intro)
        self.assertNotIn("0.14.7", intro)

    def test_two_records_none_chosen_name_neither(self):
        for f, eng in (("search_a", "sage"), ("search_b", "diann")):
            write(os.path.join(self.s, "output", f, "search_provenance.json"),
                  {"engine": eng, "version": "1"})
        intro = self.intro(self.brief("--engine", "sage", de=self.DE))
        self.assertIn("not recorded (2 search records, none chosen)", intro)
        self.assertNotIn("sage", intro)


class TheHeatmap(Base):
    def test_a_most_variable_heatmap_goes_to_qc(self):
        figs = os.path.join(self.s, "output", "figures")
        write(os.path.join(figs, "figures.json"), {"figures": [
            {"file": "heatmap_top.png", "type": "heatmap",
             "caption": "Heatmap of the 40 most variable proteins (only 2 passed adj.P < 0.05, "
                        "so the proteins shown were chosen by variance, not by the DE test)."}],
            "failed": []})
        write(os.path.join(figs, "heatmap_top.png"), "png")
        b = self.brief("--figures-dir", figs, de=TheEngine.DE)
        self.assertIn("the heatmap (the most variable proteins, chosen by variance, not by the DE "
                      "test) in **QC Assessment**", b)
        self.assertNotIn("Embed the heatmap and use it to show which proteins drive", b)


if __name__ == "__main__":
    unittest.main()
