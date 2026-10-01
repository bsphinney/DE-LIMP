#!/usr/bin/env python3
"""A Sage DE describes the quantification it actually did (architectural rule 1).

adapt_sage writes one row per PEPTIDE x run; build_maxlfq.R then takes max(PG.MaxLFQ) per protein
and run -- the single most intense peptide. Until 2.9 every Sage DE was nonetheless described
as "MaxLFQ + limma" / "DIA-NN PG.MaxLFQ" / "Quantification: DIA-NN MaxLFQ (Demichev et al. 2020)"
in methods.txt, de_provenance.json, the Methods, the AI brief and the reproducibility log.

Now the adapted report declares its quantity (parquet metadata delimp.quantity.*), build_maxlfq.R
builds the one descriptor from that declaration, and everything downstream reads the
descriptor. A DIA-NN report declares nothing and keeps its description.
"""
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

from test_run_de_contaminants import SYNTH_R, r_has, rscript   # noqa: E402

NEEDS = ("limma", "arrow", "dplyr", "tidyr", "jsonlite")
FORBIDDEN = ("DIA-NN", "MaxLFQ")
RUNS = [f"r{i:02d}" for i in range(1, 7)]
GROUPS = ["A"] * 3 + ["B"] * 3

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False


def write_sage_lfq(path, n_prot=60, n_pep=3, seed=7):
    """Sage-shaped lfq.parquet: n_prot proteins x n_pep peptides x 6 runs, each target with its
    decoy, a few targets failing q 0.05, the first 5 proteins 3x up in group B."""
    rnd = random.Random(seed)
    cols = {k: [] for k in ("peptide", "stripped_peptide", "charge", "proteins", "is_decoy",
                            "q_value", "filename", "intensity")}
    for i in range(n_prot):
        base = rnd.uniform(1e5, 1e7)
        for k in range(n_pep):
            pep = f"PEP{i}K{k}"
            q = 0.2 if (i + k) % 17 == 0 else 0.004
            for dec in (False, True):
                for r, g in zip(RUNS, GROUPS):
                    v = base * (0.5 + k) * rnd.uniform(0.8, 1.25)
                    if g == "B" and i < 5 and not dec:
                        v *= 3
                    cols["peptide"].append(pep)
                    cols["stripped_peptide"].append(pep)
                    cols["charge"].append(None)
                    cols["proteins"].append(f"sp|P{i:05d}|PROT{i}_HUMAN")
                    cols["is_decoy"].append(dec)
                    cols["q_value"].append(q)
                    cols["filename"].append(r + ".mzML")
                    cols["intensity"].append(v * (0.3 if dec else 1.0))
    pq.write_table(pa.table({**{k: v for k, v in cols.items()
                                if k not in ("charge", "q_value", "intensity")},
                             "charge": pa.array(cols["charge"], pa.int32()),
                             "q_value": pa.array(cols["q_value"], pa.float32()),
                             "intensity": pa.array(cols["intensity"], pa.float32())}), path)


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class SageDeDescribesItself(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.search = os.path.join(cls.tmp, "search")
        os.makedirs(cls.search)
        # the Sage v0.14.7 release writes its CRATE version, 0.14.6, into results.json
        with open(os.path.join(cls.search, "results.json"), "w") as fh:
            json.dump({"version": "0.14.6", "quant": {"lfq": True}}, fh)
        write_sage_lfq(os.path.join(cls.search, "lfq.parquet"))
        import run_search
        run_search.ENGINE_VERSION = None               # an --adapt-only run: no main()
        cls.report = run_search.adapt_sage(cls.search)
        with open(os.path.join(cls.tmp, "conditions.csv"), "w") as fh:
            fh.write("File.Name,Group\n" + "".join(f"{r},{g}\n" for r, g in zip(RUNS, GROUPS)))
        cls.out = os.path.join(cls.tmp, "de")
        cls.p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"), "--input", cls.report,
                                "--metadata", "conditions.csv", "--method", "maxlfq",
                                "--outdir", cls.out], capture_output=True, text=True, cwd=cls.tmp,
                               timeout=600)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def setUp(self):
        self.assertEqual(self.p.returncode, 0, self.p.stderr[-2000:])

    def prov(self):
        with open(os.path.join(self.out, "de_provenance.json")) as fh:
            return json.load(fh)

    def assertClean(self, text, where):
        for w in FORBIDDEN:
            self.assertNotIn(w, text, f"{where} says {w!r} for a Sage DE")

    def test_the_report_declares_its_quantity(self):
        meta = pq.read_schema(self.report).metadata
        self.assertEqual(meta[b"delimp.quantity.level"], b"peptide")
        self.assertIn(b"Sage 0.14.6 per its results.json (the v0.14.7 release also reports "
                      b"0.14.6)", meta[b"delimp.quantity.value"])
        self.assertEqual(meta[b"delimp.quantity.q_columns"], b"placeholder")

    def test_methods_txt_says_what_was_done(self):
        with open(os.path.join(self.out, "methods.txt")) as fh:
            text = fh.read()
        self.assertClean(text, "methods.txt")
        self.assertIn("highest Sage 0.14.6 per its results.json", text)
        # placeholder q-columns: no bare numeric cutoff, and the filter line says it kept all
        self.assertNotIn("ID FDR cutoff", text)
        self.assertIn("placeholder q-value columns", text)
        self.assertIn("Lazear 2023", text)
        self.assertIn("Caveat", text)
        self.assertIn("q-value columns are 0.0 placeholders", text)

    def test_the_provenance_record_says_what_was_done(self):
        prov = self.prov()
        self.assertEqual(prov["pipeline_id"], "peptide_max")
        self.assertIn("highest", prov["rollup_method"])
        self.assertTrue(prov["caveat"])
        self.assertEqual(prov["declared_quantity"]["level"], "peptide")
        # Everything the record says about the pipeline. The contaminant block is left out on
        # purpose: its `rule` names where the Cont_ rule comes from (DIA-NN's
        # --cont-quant-exclude, contaminants.R's one definition), not what quantified this run.
        rest = {k: v for k, v in prov.items() if k != "contaminants"}
        self.assertClean(json.dumps(rest), "de_provenance.json")

    def test_the_methods_paragraph_and_the_ai_brief(self):
        import make_methods
        prov = self.prov()
        para = make_methods.de_paragraph(prov)
        self.assertClean(para, "make_methods.de_paragraph")
        self.assertIn("Highest-peptide intensity + limma", para)
        self.assertIn("Sage's own FDR", para)
        self.assertTrue(prov["plain_language"].startswith("**Highest-peptide intensity"))
        self.assertClean(prov["plain_language"], "plain_language")

    def test_the_ai_brief_explains_the_pipeline_it_ran(self):
        brief = os.path.join(self.tmp, "ANALYSIS_PROMPT.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                            "--out", brief, "--de-dir", self.out, "--engine", "sage",
                            "--acquisition", "DDA", "--conditions", "conditions.csv"],
                           capture_output=True, text=True, cwd=self.tmp, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr[-1500:])
        text = open(brief, encoding="utf-8").read()
        self.assertIn("**Highest-peptide intensity + limma**", text)
        self.assertNotIn("MaxLFQ reconstructs protein quantities", text)
        self.assertClean(text, "the AI brief")

    def test_the_reproducibility_log_comments(self):
        path = os.path.join(self.out, "reproducibility_log.R")
        if not os.path.exists(path):
            self.skipTest("reproducibility_log.R not written in this environment")
        with open(path) as fh:
            comments = [ln for ln in fh if ln.lstrip().startswith("#")]
        text = "".join(comments).replace("PG.MaxLFQ", "<column>")   # the contract column's NAME
        self.assertClean(text, "reproducibility_log.R comments")
        self.assertIn("highest Sage", text)

    def test_audit_carries_the_caveat(self):
        audit = os.path.join(self.tmp, "AUDIT.md")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"), "--out", audit,
                        "--de-dir", self.out], capture_output=True, text=True, cwd=self.tmp,
                       timeout=120, check=True)
        text = open(audit).read()
        self.assertIn("⚠️ **quantification**", text)
        self.assertIn("single most intense peptide", text)


@unittest.skipUnless(r_has(*NEEDS), "needs R with " + "/".join(NEEDS))
class SageVersion(unittest.TestCase):
    """The version the methods name: what run_search.py recorded for the engine, else Sage's own
    results.json -- whose 0.14.6 from the v0.14.7 release is said as such -- else "not recorded"."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        import run_search
        self.rs = run_search
        self.addCleanup(setattr, run_search, "ENGINE_VERSION", run_search.ENGINE_VERSION)
        run_search.ENGINE_VERSION = None

    def value(self):
        return self.rs.sage_quantity_declaration(self.d, 0.05)["value"]

    def test_the_recorded_version_wins(self):
        with open(os.path.join(self.d, "results.json"), "w") as fh:
            json.dump({"version": "0.14.6"}, fh)
        with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "sage", "version": "v0.14.7"}, fh)
        self.assertTrue(self.value().startswith("Sage v0.14.7 label-free"))

    def test_nothing_recorded_is_said(self):
        self.assertTrue(self.value().startswith("Sage (version not recorded) label-free"))


@unittest.skipUnless(r_has(*NEEDS), "needs R with " + "/".join(NEEDS))
class DiannReportKeepsItsDescription(unittest.TestCase):
    def test_an_undeclared_report_is_still_diann_maxlfq(self):
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, True)
        rscript(f'OUT <- "{tmp}"\n' + SYNTH_R)
        out = rscript(f'setwd("{SCRIPTS}"); source("build_maxlfq.R"); '
                      f'm <- build_maxlfq(file.path("{tmp}", "report.parquet")); '
                      f'cat(m$descriptor$pipeline_id, m$descriptor$rollup_method, sep = "|"); '
                      f'cat("|", is.null(m$declared_quantity), sep = "")')
        self.assertEqual(out.strip(), "maxlfq|DIA-NN PG.MaxLFQ|TRUE")


if __name__ == "__main__":
    unittest.main()
