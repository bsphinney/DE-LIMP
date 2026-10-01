#!/usr/bin/env python3
"""Every adapted report says what its numbers are, and the DE says that -- not "DIA-NN".

adapt_fragpipe_dda (IonQuant combined_protein.tsv), adapt_alphadia (pg.matrix) and adapt_radiant
(Fulcrum) write one engine-made value per protein and run into the DE contract's PG.MaxLFQ
column. build_maxlfq.R used to describe every such DE as "MaxLFQ + limma" / "DIA-NN
PG.MaxLFQ" / "Quantification: DIA-NN MaxLFQ (Demichev et al. 2020)" -- Michelle's SET28
FragPipe 23 session (2026-09-25) carries exactly that. Each adapter now declares its quantity
(level "protein", the engine's own name and version, its FDR, the Crossref-verified citation),
and maxlfq_descriptor()'s protein branch passes it through with no re-rollup described.
test_sage_de_describes_itself.py covers Sage (level "peptide") and the undeclared DIA-NN control.
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

from test_run_de_contaminants import r_has   # noqa: E402

NEEDS = ("limma", "arrow", "dplyr", "tidyr", "jsonlite")
RUNS = [f"r{i:02d}" for i in range(1, 7)]
GROUPS = ["A"] * 3 + ["B"] * 3

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False


def protein_values(n_prot=60, seed=11):
    """{protein index: {run: intensity}} with the first 5 proteins 3x up in group B."""
    rnd = random.Random(seed)
    vals = {}
    for i in range(n_prot):
        base = rnd.uniform(1e5, 1e8)
        vals[i] = {r: base * rnd.uniform(0.8, 1.25) * (3 if g == "B" and i < 5 else 1)
                   for r, g in zip(RUNS, GROUPS)}
    return vals


def fragpipe_output(out):
    vals = protein_values()
    hdr = (["Protein", "Protein ID", "Entry Name", "Gene"]
           + [f"{r} Intensity" for r in RUNS] + [f"{r} MaxLFQ Intensity" for r in RUNS])
    with open(os.path.join(out, "combined_protein.tsv"), "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for i, per in vals.items():
            fh.write("\t".join([f"sp|P{i:05d}|PROT{i}_HUMAN", f"P{i:05d}", f"PROT{i}_HUMAN",
                                f"G{i}"] + [f"{per[r] * 1.1:.1f}" for r in RUNS]
                               + [f"{per[r]:.1f}" for r in RUNS]) + "\n")
    with open(os.path.join(out, "fragpipe.workflow"), "w") as fh:
        fh.write("# Version info:\n# FragPipe version 23.1\n# MSFragger version 4.4.1\n"
                 "# IonQuant version 1.11.20\n# DIA-NN version 1.8.2 beta 8\n"
                 "ionquant.maxlfq=1\nionquant.mbr=1\n"
                 "phi-report.filter=--sequential --prot 0.01 --picked\n")


# frozen_config.yaml / pg.matrix as AlphaDIA writes them. Shapes read from real HIVE runs
# (2026-09-29): a 2.0.1 service run (SERVICE/on_campus/PI_Example_B/...) -- `version: 2.0.1`, no
# search_output.normalization_method, pg.matrix.parquet; and 2.1.2 runs (brett/glendon/...) --
# `version: 2.1.2`, normalization_method: directlfq, file_format tsv -> pg.matrix.tsv whose id
# column is pg.name. Before 2.0.x the file carried the config SCHEMA version (`version: 1`,
# constants/default.yaml at v1.10.3 .. v2.0.0), which is never taken as the release.
FROZEN_1_10 = ("version: 1\nfdr:\n  fdr: 0.01\n  group_level: proteins\n"
               "search_output:\n  peptide_level_lfq: false\n  file_format: tsv\n")
FROZEN_2_0_1 = ("version: 2.0.1\nfdr:\n  fdr: 0.01\n  group_level: proteins\n"
                "search_output:\n  file_format: parquet\n  precursor_level_lfq: true\n")
FROZEN_2_1_2 = ("version: 2.1.2\nfdr:\n  fdr: 0.01\nsearch_output:\n  file_format: tsv\n"
                "  precursor_level_lfq: true\n  normalization_method: directlfq\n"
                "  normalize_directlfq: true\n")


def alphadia_output(out, frozen=FROZEN_2_0_1, tsv=False):
    vals = protein_values()
    if tsv:
        with open(os.path.join(out, "pg.matrix.tsv"), "w") as fh:
            fh.write("pg.name\t" + "\t".join(RUNS) + "\n")
            for i in vals:
                fh.write(f"P{i:05d}\t" + "\t".join(f"{vals[i][r]:.1f}" for r in RUNS) + "\n")
    else:
        pq.write_table(pa.table({"pg": [f"P{i:05d}" for i in vals],
                                 **{r: [vals[i][r] for i in vals] for r in RUNS}}),
                       os.path.join(out, "pg.matrix.parquet"))
    with open(os.path.join(out, "frozen_config.yaml"), "w") as fh:
        fh.write(frozen)


def radiant_output(out, column="PG.Normalised"):
    vals = protein_values()
    d = os.path.join(out, "radiant_results", "fulcrum-results")
    os.makedirs(d)
    rows = {"Run": [], "Protein.Group": [], column: [], "Global.PG.Q.Value": []}
    for i, per in vals.items():
        for r in RUNS:
            for _ in range(2):                          # precursor rows, PG value broadcast
                rows["Run"].append(f"file:///mnt/results/radiant-results/{r}.mzML.radiantDIA")
                rows["Protein.Group"].append(f"P{i:05d}")
                rows[column].append(per[r])
                rows["Global.PG.Q.Value"].append(0.001)
    pq.write_table(pa.table(rows), os.path.join(d, "part-00000.parquet"))
    with open(os.path.join(out, "search_provenance.json"), "w") as fh:
        json.dump({"engine": "radiant", "version": "2.3.3"}, fh)


class _Case(unittest.TestCase):
    """Adapt one engine's output, run the DE on it, keep the records."""
    make = adapt = None

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.search = os.path.join(cls.tmp, "search")
        os.makedirs(cls.search)
        cls.make(cls.search)
        import run_search
        run_search.ENGINE_VERSION = None               # an --adapt-only run: no main()
        cls.report = getattr(run_search, cls.adapt)(cls.search)
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

    def records(self):
        self.assertEqual(self.p.returncode, 0, self.p.stderr[-2000:])
        with open(os.path.join(self.out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        with open(os.path.join(self.out, "methods.txt")) as fh:
            methods = fh.read()
        import make_methods
        brief = os.path.join(self.tmp, "ANALYSIS_PROMPT.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"), "--out",
                            brief, "--de-dir", self.out, "--conditions", "conditions.csv"],
                           capture_output=True, text=True, cwd=self.tmp, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr[-1500:])
        with open(brief, encoding="utf-8") as fh:
            brief_text = fh.read()
        return prov, methods, make_methods.de_paragraph(prov), brief_text

    def check(self, label, must):
        prov, methods, para, brief = self.records()
        self.assertEqual(prov["pipeline_id"], "engine_protein")
        self.assertEqual(prov["display_label"], f"{label} + limma")
        self.assertEqual(prov["declared_quantity"]["level"], "protein")
        self.assertIn("no re-rollup here", prov["rollup_method"])
        self.assertIn(f"**{label} + limma**", brief)
        # the contaminant block is left out: its `rule` cites DIA-NN's --cont-quant-exclude as
        # where the Cont_ rule comes from, not what quantified this run
        described = json.dumps({k: v for k, v in prov.items() if k != "contaminants"})
        for where, text in (("methods.txt", methods), ("de_provenance.json", described),
                            ("make_methods", para), ("the AI brief", brief)):
            self.assertNotIn("DIA-NN", text, f"{where} names DIA-NN for a {label} DE")
            for m in must:
                if where != "the AI brief":
                    self.assertIn(m, text, f"{where} does not say {m!r}")
        return prov


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class FragPipeDda(_Case):
    make, adapt = staticmethod(fragpipe_output), "adapt_fragpipe"

    def test_the_de_names_ionquant_not_diann(self):
        prov = self.check("FragPipe IonQuant MaxLFQ",
                          ["FragPipe 23.1 / IonQuant 1.11.20 MaxLFQ intensity, "
                           "match-between-runs on",
                           "Philosopher filter --sequential --prot 0.01 --picked"])
        self.assertEqual(prov["q_columns_role"], "placeholder")
        _, methods, para, _ = self.records()
        # placeholder q-columns: no bare numeric cutoff anywhere, and the filter line says so
        self.assertNotIn("ID FDR cutoff", methods)
        self.assertIn("placeholder q-value columns, all 0.0", methods)
        self.assertNotIn("filtered at Q.Value", para)
        self.assertIn("Yu et al. 2021", prov["citation"])
        self.assertIn("Cox et al. 2014", prov["citation"])
        self.assertIsNone(prov.get("caveat"))


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class AlphaDia(_Case):
    make, adapt = staticmethod(alphadia_output), "adapt_alphadia"

    def test_a_2_0_1_run_is_directlfq_without_the_key(self):
        """2.0.1 records its version but no normalization_method: before 2.1 directLFQ is
        AlphaDIA's only protein LFQ, so it is named -- and the file is not called missing."""
        prov = self.check("AlphaDIA directLFQ",
                          ["AlphaDIA 2.0.1 protein-group LFQ intensity (directLFQ; "
                           "pg.matrix.parquet)", "precursors and protein groups at 0.01"])
        self.assertNotIn("no frozen_config.yaml", json.dumps(prov))
        self.assertIn("Wallmann et al. 2025", prov["citation"])
        self.assertIn("Ammar et al. 2023", prov["citation"])


def alphadia_2_1_output(out):
    alphadia_output(out, frozen=FROZEN_2_1_2, tsv=True)


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class AlphaDia21(_Case):
    make, adapt = staticmethod(alphadia_2_1_output), "adapt_alphadia"

    def test_a_2_1_2_tsv_run_names_its_own_version(self):
        prov = self.check("AlphaDIA directLFQ",
                          ["AlphaDIA 2.1.2 protein-group LFQ intensity (directLFQ; "
                           "pg.matrix.tsv)"])
        self.assertIn("Ammar et al. 2023", prov["citation"])


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class Radiant(_Case):
    make, adapt = staticmethod(radiant_output), "adapt_radiant"

    def test_the_de_names_radiant_fulcrum(self):
        prov = self.check("Radiant Fulcrum",
                          ["Radiant 2.3.3 Fulcrum PG.Normalised per protein group",
                           "Fulcrum's Global.PG.Q.Value"])
        self.assertIn("Just et al. 2026", prov["citation"])
        self.assertIsNone(prov.get("caveat"))
        # Fulcrum's q-values are real: the numeric cutoff stays in the methods
        self.assertEqual(prov["q_columns_role"], "real")
        _, methods, para, _ = self.records()
        self.assertIn("ID FDR cutoff : q <= 0.010", methods)
        self.assertIn("≤ 0.01", para)


@unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
class DeclarationDetails(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        import run_search
        self.rs = run_search
        self.addCleanup(setattr, run_search, "ENGINE_VERSION", run_search.ENGINE_VERSION)
        run_search.ENGINE_VERSION = None

    def test_a_precursor_level_fulcrum_column_is_declared_peptide_level(self):
        """The adapter kept the highest precursor per protein and run: that is the peptide-level
        rollup (build_maxlfq.R's peptide branch describes it and carries its caveat), not a
        protein quantity."""
        dec = self.rs.radiant_quantity_declaration(self.d, "Precursor.Normalised", None)
        self.assertEqual(dec["level"], "peptide")
        self.assertEqual(dec["value"], "Radiant (version not recorded) Fulcrum "
                                       "Precursor.Normalised precursor intensity")
        self.assertEqual(dec["q_columns"], "placeholder")
        self.assertIn("0.0 placeholders", dec["identification"])
        dec = self.rs.radiant_quantity_declaration(self.d, "PG.Quantity", "Global.PG.Q.Value")
        self.assertEqual((dec["level"], dec["q_columns"]), ("protein", "real"))

    def test_missing_version_sources_are_said_not_guessed(self):
        dec = self.rs.fragpipe_quantity_declaration(self.d, maxlfq_columns=False)
        self.assertIn("FragPipe (version not recorded) / IonQuant (version not recorded) "
                      "protein intensity (not MaxLFQ)", dec["value"])
        self.assertIn("not recorded: no fragpipe.workflow", dec["identification"])
        self.assertNotIn("Cox et al.", dec["citation"])
        dec = self.rs.alphadia_quantity_declaration(self.d, "pg.matrix.parquet")
        self.assertIn("method not recorded: no frozen_config.yaml", dec["value"])
        self.assertEqual(dec["label"], "AlphaDIA")

    def test_a_frozen_config_without_the_method_on_2_1_says_so(self):
        with open(os.path.join(self.d, "frozen_config.yaml"), "w") as fh:
            fh.write("version: 2.1.4\nfdr:\n  fdr: 0.01\n")
        dec = self.rs.alphadia_quantity_declaration(self.d, "pg.matrix.parquet")
        self.assertIn("AlphaDIA 2.1.4", dec["value"])
        self.assertIn("method not recorded in frozen_config.yaml", dec["value"])
        self.assertNotIn("directLFQ", dec["citation"])

    def test_a_schema_version_is_never_the_release(self):
        with open(os.path.join(self.d, "frozen_config.yaml"), "w") as fh:
            fh.write(FROZEN_1_10)
        dec = self.rs.alphadia_quantity_declaration(self.d, "pg.matrix.tsv")
        self.assertTrue(dec["value"].startswith("AlphaDIA (version not recorded) protein-group "
                                                "LFQ intensity (directLFQ;"), dec["value"])
        self.assertIn("directLFQ", dec["citation"])
        # what run_search.py recorded for the engine fills the gap
        with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "alphadia", "version": "1.10.3"}, fh)
        self.assertIn("AlphaDIA 1.10.3 protein-group", self.rs.alphadia_quantity_declaration(
            self.d, "pg.matrix.tsv")["value"])

    def test_fragpipe_mbr_off_and_unrecorded(self):
        with open(os.path.join(self.d, "fragpipe.workflow"), "w") as fh:
            fh.write("# FragPipe version 23.1\nionquant.mbr=0\n")
        self.assertIn("match-between-runs off",
                      self.rs.fragpipe_quantity_declaration(self.d, True)["value"])
        os.remove(os.path.join(self.d, "fragpipe.workflow"))
        self.assertIn("match-between-runs not recorded (no fragpipe.workflow)",
                      self.rs.fragpipe_quantity_declaration(self.d, True)["value"])

    def test_main_s_engine_version_is_used_when_the_output_has_none(self):
        self.rs.ENGINE_VERSION = "2.3.3"
        self.assertIn("Radiant 2.3.3", self.rs.radiant_quantity_declaration(
            self.d, "PG.Quantity", "Q.Value")["value"])

    def test_frozen_config_without_pyyaml(self):
        with open(os.path.join(self.d, "frozen_config.yaml"), "w") as fh:
            fh.write("version: 1.12.0\nfdr:\n  fdr: 0.01\nsearch_output:\n"
                     "  normalization_method: 'quantselect'\n")
        import builtins
        real = builtins.__import__

        def no_yaml(name, *a, **k):
            if name == "yaml":
                raise ImportError("no yaml")
            return real(name, *a, **k)
        builtins.__import__ = no_yaml
        try:
            got = self.rs._frozen_config(self.d)
        finally:
            builtins.__import__ = real
        self.assertEqual((got["version"], got["normalization_method"], got["fdr"]),
                         ("1.12.0", "quantselect", "0.01"))


if __name__ == "__main__":
    unittest.main()
