#!/usr/bin/env python3
"""The contaminant filter works on every engine's report, and says what each engine did.

contaminants.R's rule is DIA-NN's: an accession in the group STARTS with Cont_ (or FragPipe's
contam_). DIA-NN writes bare accessions (Cont_P02769). adapt_sage wrote Sage's `proteins` as
Sage gives it -- full FASTA IDs (sp|Cont_P02769|ALBU_BOVIN) -- so no Sage DE ever had a
contaminant removed: on gabrig's HeL50 UnvPe (2026-09-29) the Sage DE limma-tested 171 Cont_
groups (BSA, gelsolin, mouse keratins) while methods.txt said "Contaminants : none", and a
keratin sample's exemption could never apply. adapt_sage now writes bare accessions through
protein_ids.group_accessions (the one reading of an identifier). FragPipe's DDA adapter reads
"Protein ID" (Cont_A2I7N3 on the same data) and keeps FragPipe's contam_ tag.

And a keratin sample's kept keratin was described, for every engine, as DIA-NN's
--cont-quant-exclude having left its peptides out. Each pipeline's descriptor now says what its
engine did (kept_contaminant_quant): a Sage or FragPipe DE is never told DIA-NN did anything.
"""
import csv
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

import protein_ids   # noqa: E402
from test_run_de_contaminants import r_has   # noqa: E402
from test_keratin_sample import SIDECAR, REPORT_KERATIN, OWN_KERATIN   # noqa: E402

NEEDS = ("limma", "arrow", "dplyr", "tidyr", "jsonlite")
RUNS = [f"r{i:02d}" for i in range(1, 7)]
GROUPS = ["A"] * 3 + ["B"] * 3

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

# the contaminant and keratin entries every fixture carries, as their FASTA headers name them
BSA = "sp|Cont_P02769|ALBU_BOVIN"
TRYPSIN = "sp|Cont_P00761|TRYP_PIG"
FP_TRYPSIN = "contam_sp|P00761|TRYP_PIG"            # FragPipe/Philosopher --contam
OWN_KRTAP = "sp|Cont_Q9BYR0|KRA47_HUMAN"            # the sample species' own keratin
WOOL = "sp|Cont_P02438|KRB2A_SHEEP"                 # another species', naming no sample protein
SHARED_KRT = "sp|P31000|K1H1_HUMAN;sp|Cont_Q61765|K1H1_MOUSE"
# a sample protein (P00001, which has 3 peptides of its own) sharing one peptide with BSA: a
# MIXED Cont_ group. On gabrig's HeL50 Sage run 101 of the 171 Cont_ groups were mixed and 36 human
# proteins lost peptides to them, while methods.txt said "0 lost some" (sage-review B1).
MIXED_BSA = f"sp|P00001|PROT1_HUMAN;{BSA}"


def sample_proteins(n=60, seed=5):
    rnd = random.Random(seed)
    return [(f"sp|P{i:05d}|PROT{i}_HUMAN", rnd.uniform(1e5, 1e7), i < 5) for i in range(n)]


def intensity(base, g, up, rnd):
    return base * rnd.uniform(0.8, 1.25) * (3 if (up and g == "B") else 1)


class ProteinIdsReadAGroup(unittest.TestCase):
    def test_a_sage_group_becomes_bare_accessions(self):
        self.assertEqual(protein_ids.group_accessions(f"sp|P12345|ALBU_HUMAN;{BSA}"),
                         "P12345;Cont_P02769")
        self.assertEqual(protein_ids.group_accessions(FP_TRYPSIN), "contam_P00761")
        self.assertEqual(protein_ids.group_accessions("P1;Cont_P2;P1"), "P1;Cont_P2")   # as is
        self.assertEqual(protein_ids.group_entry_names(SHARED_KRT), "K1H1_HUMAN;K1H1_MOUSE")
        self.assertEqual(protein_ids.group_entry_names("P1;P2"), "")


def write_sage_output(out, seed=9):
    """Sage 0.14.7-shaped lfq.parquet: one row per peptide x run, proteins as FASTA IDs."""
    rnd = random.Random(seed)
    rows = [(p, base, up, k) for p, base, up in sample_proteins() for k in range(3)]
    rows += [(BSA, 5e6, False, k) for k in range(3)] + [(TRYPSIN, 4e6, False, k) for k in range(2)]
    rows += [(OWN_KRTAP, 3e6, True, k) for k in range(3)] + [(WOOL, 2e6, True, k) for k in range(3)]
    rows += [(SHARED_KRT, 6e6, True, 0), (MIXED_BSA, 5e6, False, 0)]
    cols = {k: [] for k in ("peptide", "stripped_peptide", "proteins", "is_decoy", "q_value",
                            "filename", "intensity")}
    for j, (prot, base, up, k) in enumerate(rows):
        for r, g in zip(RUNS, GROUPS):
            cols["peptide"].append(f"PEP{j}K")
            cols["stripped_peptide"].append(f"PEP{j}K")
            cols["proteins"].append(prot)
            cols["is_decoy"].append(False)
            cols["q_value"].append(0.001)
            cols["filename"].append(r + ".mzML")
            cols["intensity"].append(intensity(base * (0.5 + k), g, up, rnd))
    os.makedirs(out, exist_ok=True)
    pq.write_table(pa.table({**{k: v for k, v in cols.items() if k not in ("q_value", "intensity")},
                             "charge": pa.array([None] * len(cols["peptide"]), pa.int32()),
                             "q_value": pa.array(cols["q_value"], pa.float32()),
                             "intensity": pa.array(cols["intensity"], pa.float32())}),
                   os.path.join(out, "lfq.parquet"))
    with open(os.path.join(out, "results.json"), "w") as fh:
        json.dump({"version": "0.14.6", "quant": {"lfq": True}}, fh)


def write_fragpipe_output(out, seed=13):
    """combined_protein.tsv as FragPipe 23.1 writes it (Protein, Protein ID, Entry Name, Gene,
    '<sample> MaxLFQ Intensity'); Protein ID of a skill contaminant is Cont_... (gabrig's HeL50
    UnvPe MSFragger output, 2026-09-29), of a Philosopher --contam one the bare accession."""
    rnd = random.Random(seed)
    rows = [(p, p.split("|")[1], base, up) for p, base, up in sample_proteins()]
    rows += [(BSA, "Cont_P02769", 5e6, False), (FP_TRYPSIN, "P00761", 4e6, False),
             (OWN_KRTAP, "Cont_Q9BYR0", 3e6, True), (WOOL, "Cont_P02438", 2e6, True)]
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "combined_protein.tsv"), "w") as fh:
        fh.write("\t".join(["Protein", "Protein ID", "Entry Name", "Gene"]
                           + [f"{r} MaxLFQ Intensity" for r in RUNS]) + "\n")
        for prot, pid, base, up in rows:
            fh.write("\t".join([prot, pid, prot.split("|")[2], ""]
                               + [f"{intensity(base, g, up, rnd):.1f}"
                                  for r, g in zip(RUNS, GROUPS)]) + "\n")
    with open(os.path.join(out, "fragpipe.workflow"), "w") as fh:
        fh.write("# FragPipe version 23.1\n# IonQuant version 1.11.20\nionquant.mbr=1\n"
                 "phi-report.filter=--sequential --prot 0.01 --picked\n")


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class EveryEngineFiltersAndSaysWhatItDid(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = t = tempfile.mkdtemp()
        import run_search
        run_search.ENGINE_VERSION = None
        write_sage_output(os.path.join(t, "sage"))
        write_fragpipe_output(os.path.join(t, "fragpipe"))
        reports = {"sage": run_search.adapt_sage(os.path.join(t, "sage")),
                   "fragpipe": run_search.adapt_fragpipe(os.path.join(t, "fragpipe"))}
        with open(os.path.join(t, "conditions.csv"), "w") as fh:
            fh.write("File.Name,Group\n" + "".join(f"{r},{g}\n" for r, g in zip(RUNS, GROUPS)))
        # the searched database's keratin entries, as fetch_fasta.py's sidecar lists them (the user
        # said "not keratin" when it was built; the DE is told --keratin-sample); it holds both
        # tags, so no contaminant here is one the check never saw
        with open(os.path.join(t, "cont.json"), "w") as fh:
            json.dump(dict(SIDECAR, keratin_sample=False, keratin_sample_source="user",
                           keratin_contaminants_in_database=sorted(REPORT_KERATIN),
                           keratin_contaminants_same_species=sorted(OWN_KERATIN),
                           contaminant_tags_in_database=["Cont_", "contam_"]), fh)
        cls.runs = {}
        for engine, report in reports.items():
            for key, extra in (("plain", []),
                               ("keratin", ["--keratin-sample", "--fasta-meta",
                                            os.path.join(t, "cont.json")])):
                out = os.path.join(t, f"{engine}_{key}")
                p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"), "--input",
                                    report, "--metadata", "conditions.csv", "--method", "maxlfq",
                                    "--outdir", out, *extra], capture_output=True, text=True,
                                   cwd=t, timeout=600)
                cls.runs[(engine, key)] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, engine, key):
        p, out = self.runs[(engine, key)]
        self.assertEqual(p.returncode, 0, p.stderr[-2500:])
        return out

    def matrix(self, engine, key):
        with open(os.path.join(self.out(engine, key), "Expression_Matrix.csv"), newline="") as fh:
            return {r["Protein.Group"] for r in csv.DictReader(fh)}

    def prov(self, engine, key):
        with open(os.path.join(self.out(engine, key), "de_provenance.json")) as fh:
            return json.load(fh)

    def methods(self, engine, key):
        with open(os.path.join(self.out(engine, key), "methods.txt")) as fh:
            return fh.read()

    def test_a_sage_contaminant_is_removed(self):
        em = self.matrix("sage", "plain")
        for g in em:
            self.assertFalse(g.startswith("sp|"), f"{g}: Sage ids must be bare accessions")
        self.assertIn("P00000", em)
        self.assertFalse({"Cont_P02769", "Cont_P00761", "Cont_Q9BYR0"} & em)
        c = self.prov("sage", "plain")["contaminants"]
        self.assertEqual(c["policy"], "removed")
        self.assertGreaterEqual(c["n_protein_groups"], 4)
        self.assertIn("Contaminants  : REMOVED", self.methods("sage", "plain"))
        self.assertNotIn("Contaminants  : none", self.methods("sage", "plain"))

    def test_a_fragpipe_contaminant_is_removed_both_tags(self):
        em = self.matrix("fragpipe", "plain")
        self.assertIn("P00000", em)
        self.assertFalse({"Cont_P02769", "contam_P00761", "Cont_Q9BYR0", "P00761"} & em)
        c = self.prov("fragpipe", "plain")["contaminants"]
        self.assertEqual(c["policy"], "removed")
        self.assertEqual(c["n_protein_groups"], 4)

    def test_a_sage_keratin_sample_keeps_its_exempt_keratins(self):
        em = self.matrix("sage", "keratin")
        self.assertIn("Cont_Q9BYR0", em)                   # the species' own KRTAP
        self.assertIn("P31000;Cont_Q61765", em)             # a sample keratin sharing with mouse
        self.assertNotIn("Cont_P02438", em)                 # sheep wool alone: a contaminant
        self.assertFalse({"Cont_P02769", "Cont_P00761"} & em)
        k = self.prov("sage", "keratin")["contaminants"]["keratin_sample"]
        self.assertEqual(k["database"], "keratins_in_database")
        self.assertIs(k["requantified"], True)
        self.assertGreater(k["n_precursors_kept"], 0)

    def test_the_keratin_lines_say_what_each_engine_did(self):
        for engine, label, says in (
                ("sage", "Sage", "Each protein's value is taken here from the peptide rows the "
                                 "filter kept, so those keratin peptides count in it; Sage was "
                                 "given no contaminant exclusion."),
                ("fragpipe", "FragPipe IonQuant MaxLFQ",
                 "Their quantities are FragPipe IonQuant MaxLFQ's own protein quantities, used "
                 "as reported (not re-derived here); it was given no contaminant exclusion, and "
                 "whether its protein inference moved peptides shared with those entries to "
                 "another protein is not recorded.")):
            with self.subTest(engine=engine):
                m = self.methods(engine, "keratin")
                self.assertIn(says, " ".join(m.split()))
                self.assertNotIn("DIA-NN --cont-quant-exclude", m)
                self.assertNotIn("under-quantified", m)
                import make_methods
                sentence = make_methods.de_keratin_sentence(
                    self.prov(engine, "keratin")["contaminants"])
                self.assertIn(says, sentence)
                self.assertNotIn("DIA-NN", sentence)
                self.assertNotIn("resolve before publication", sentence)
        # sage-review N4: Philosopher's razor assignment can move shared peptides, so FragPipe's
        # is "not known", never a FALSE nobody checked; Sage's rollup is from the kept peptides
        self.assertIsNone(self.prov("fragpipe", "keratin")["contaminants"]["keratin_sample"]
                          ["under_quantified"])
        self.assertIn("[whether keratin is under-quantified is not known",
                      __import__("make_methods").de_keratin_sentence(
                          self.prov("fragpipe", "keratin")["contaminants"]))
        self.assertIs(self.prov("sage", "keratin")["contaminants"]["keratin_sample"]
                      ["under_quantified"], False)

    def test_counts_are_distinct_items_in_the_reports_own_unit(self):
        """sage-review B1: the census counted feature x run ROWS as "precursors" -- a Sage DE
        printed 2,744 for 686 peptides x 4 runs, a FragPipe DE 288 for 72 proteins x 4 runs. Now
        each counts DISTINCT items in the unit its declared quantity level names."""
        s = self.prov("sage", "plain")["contaminants"]
        # BSA 3 + trypsin 2 + KRTAP 3 + wool 3 + the shared keratin 1 + the mixed BSA peptide 1
        self.assertEqual((s["unit"], s["unit_level"], s["n_precursors"]), ("peptides", "peptide", 13))
        self.assertEqual(s["n_protein_groups"], 6)
        f = self.prov("fragpipe", "plain")["contaminants"]
        self.assertEqual((f["unit"], f["unit_level"], f["n_precursors"]),
                         ("protein groups", "protein", 4))
        m = " ".join(self.methods("sage", "plain").split())
        self.assertIn("REMOVED before quantification -- 13 peptides mapping to a Cont_ entry", m)
        self.assertIn("(a peptide is a contaminant when any accession in", m)

    def test_a_mixed_sage_group_is_counted_per_sample_protein(self):
        """P00001 lost one peptide to a mixed Cont_ group and keeps its own three;
        P31000 (human KRT31) was named only by its peptide shared with mouse Krt31, so it is not
        quantified at all. Never "0 lost some" by default."""
        c = self.prov("sage", "plain")["contaminants"]
        self.assertEqual(c["n_protein_groups_mixed"], 2)
        self.assertEqual((c["n_sample_groups_sharing"], c["n_sample_groups_all_shared"]), (1, 1))
        self.assertEqual(c["sample_loss"],
                         "Sample proteins named by a removed peptide (a peptide shared with a "
                         "Cont_ entry): 1 keep other peptides, 1 have no other peptide and are "
                         "not quantified.")
        em = self.matrix("sage", "plain")
        self.assertIn("P00001", em)
        self.assertNotIn("P31000", em)
        m = " ".join(self.methods("sage", "plain").split())
        self.assertIn("(2 of them also named a sample protein)", m)
        self.assertIn(c["sample_loss"], m)
        import make_methods
        para = make_methods._de_contaminant_sentence(self.prov("sage", "plain"))
        self.assertIn("13 peptides mapping to a Cont_-tagged contaminant entry", para)
        self.assertIn("(2 of them also naming a sample protein)", para)
        self.assertIn(c["sample_loss"], para)
        self.assertIn(c["rule"], para)          # the rule text is contaminants.R's (rule 3)

    def test_a_fragpipe_protein_report_never_says_precursors(self):
        """A protein-level report has no precursors: methods.txt, make_methods' paragraphs and the
        keratin lines count protein groups, and which sample proteins lost peptides inside the
        engine's own quantities is said to be not computed -- not 0."""
        import make_methods
        for key in ("plain", "keratin"):
            with self.subTest(run=key):
                prov = self.prov("fragpipe", key)
                c = prov["contaminants"]
                texts = {"methods.txt": self.methods("fragpipe", key),
                         "make_methods DE": make_methods._de_contaminant_sentence(prov),
                         "make_methods keratin": make_methods.de_keratin_sentence(c)}
                for name, text in texts.items():
                    self.assertNotIn("precursor", text.lower(), name)
                self.assertIsNone(c["n_sample_groups_sharing"])
                self.assertIn("is not computed (a protein-level report)", c["sample_loss"])
        m = " ".join(self.methods("fragpipe", "plain").split())
        self.assertIn("REMOVED before quantification -- 4 protein groups mapping to a", m)
        self.assertIn("(a protein group is a contaminant when any accession in", m)
        k = " ".join(self.methods("fragpipe", "keratin").split())
        # one kept group is "1 protein group ... was KEPT", not "1 protein groups" (sage-review)
        self.assertIn("1 protein group mapping only to keratin-family contaminant entries (3 such "
                      "entries in the database) was KEPT,", k)
        self.assertIn("1 protein group mapping only to keratin-family contaminant entries (3 in "
                      "the search database) was kept,", make_methods.de_keratin_sentence(
                          self.prov("fragpipe", "keratin")["contaminants"]))

    def test_the_rule_is_quoted_without_nested_parentheses(self):
        rule = ("a peptide is a contaminant when any accession in Protein.Group starts with "
                "'Cont_' or 'contam_' -- applied here, in the DE step, to the Sage report")
        c = self.prov("sage", "plain")["contaminants"]
        self.assertEqual(c["rule"], rule)
        self.assertIn(f"({rule})", " ".join(self.methods("sage", "plain").split()))
        self.assertNotIn("report))", " ".join(self.methods("sage", "plain").split()))

    def test_the_rules_origin_is_each_pipelines_never_diann_for_every_engine(self):
        """make_methods' contaminant sentence ended "the rule of DIA-NN's --cont-quant-exclude,
        applied here" for every engine (SKILL_OPEN_DEFECTS, 2.10), which a Sage or FragPipe reader
        could take to mean DIA-NN had run. The origin is now the pipeline descriptor's
        (contaminant_rule_origin): a DIA-NN report keeps those words (test_run_de_contaminants)."""
        import make_methods
        for engine, label in (("sage", "Sage"), ("fragpipe", "FragPipe IonQuant MaxLFQ")):
            with self.subTest(engine=engine):
                c = self.prov(engine, "plain")["contaminants"]
                self.assertEqual(c["rule_origin"],
                                 f"applied here, in the DE step, to the {label} report")
                self.assertTrue(c["rule"].endswith(c["rule_origin"]))
                para = make_methods._de_contaminant_sentence(self.prov(engine, "plain"))
                self.assertIn(c["rule"], para)
                self.assertNotIn("DIA-NN", para)
                self.assertNotIn("DIA-NN", self.methods(engine, "plain").split("Contaminants")[1]
                                 .split("Keratin")[0])


@unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
class FragPipeDiaContamWarning(unittest.TestCase):
    """An older FragPipe DIA output searched on a Philosopher --contam database is adapted with a
    WARNING: its DIA route drops the tag, so those contaminants reach the DE as sample proteins."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        os.makedirs(os.path.join(self.d, "dia-quant-output"))
        pq.write_table(pa.table({"Run": ["r01"], "Protein.Group": ["P1"], "PG.MaxLFQ": [1.0],
                                 **{c: [0.001] for c in ("Q.Value", "Lib.Q.Value",
                                                         "Lib.PG.Q.Value")}}),
                       os.path.join(self.d, "dia-quant-output", "report.parquet"))

    def adapt(self, db_text):
        db = os.path.join(self.d, "db.fasta")
        with open(db, "w") as fh:
            fh.write(db_text)
        with open(os.path.join(self.d, "fragpipe.workflow"), "w") as fh:
            fh.write(f"database.db-path={db}\n")
        return subprocess.run([sys.executable, "-c", "import sys; sys.path.insert(0, sys.argv[1]); "
                               "import run_search; run_search.adapt_fragpipe(sys.argv[2])",
                               SCRIPTS, self.d], capture_output=True, text=True, timeout=60)

    def test_a_contam_database_is_warned_about(self):
        r = self.adapt(">sp|P1|ONE_HUMAN\nPEPTIDEK\n>contam_sp|P00761|TRYP_PIG\nPEPTIDER\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("[adapt] WARNING", r.stderr)
        self.assertIn("contam_", r.stderr)
        self.assertIn(": 1)", r.stderr)

    def test_a_skill_database_is_not(self):
        r = self.adapt(">sp|P1|ONE_HUMAN\nPEPTIDEK\n>sp|Cont_P00761|TRYP_PIG\nPEPTIDER\n")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertNotIn("WARNING", r.stderr)


@unittest.skipUnless(r_has("jsonlite"), "needs R with jsonlite")
class UnitsPerLevel(unittest.TestCase):
    """contaminants.R contaminant_unit is the one definition of what a count counts, and the
    census counts distinct items -- a feature seen in 4 runs is one item."""

    def census(self, level, has_feature, feature):
        expr = (f'.script_dir <- "{SCRIPTS}"; source(file.path(.script_dir, "contaminants.R")); '
                f'u <- contaminant_unit({level}, has_feature = {has_feature}); '
                f'g <- rep(c("P1", "Cont_P2", "P3;Cont_P4"), each = 4); '
                f'cs <- contaminant_census(g, is_contaminant(g), feature = {feature}, unit = u); '
                f'cat(jsonlite::toJSON(list(unit = u, n = cs$n_precursors, '
                f'mixed = cs$n_groups_mixed, part = cs$n_groups_sample_part, '
                f'whole = cs$n_groups_sample_whole, cols = names(cs$groups), '
                f'loss = contaminant_sample_loss(cs, "Cont_")), auto_unbox = TRUE, na = "null"))')
        p = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        return json.loads(p.stdout)

    def test_each_level(self):
        # 3 groups x 4 runs = 12 rows; one feature per group and run for precursors/peptides
        feats = 'rep(c("a", "b", "c"), each = 4)'
        for level, has, feature, unit, n, part, whole in (
                ("NULL", "TRUE", feats, "precursors", 2, 0, 0),
                ('"peptide"', "TRUE", feats, "peptides", 2, 0, 1),
                ('"protein"', "TRUE", "g", "protein groups", 2, None, 1),
                ("NULL", "FALSE", "NULL", "report rows", 8, None, None)):
            with self.subTest(unit=unit):
                r = self.census(level, has, feature)
                self.assertEqual(r["unit"]["plural"], unit)
                self.assertEqual(r["n"], n)
                self.assertEqual(r["mixed"], 1)
                self.assertEqual((r["part"], r["whole"]), (part, whole))
                if unit in ("peptides", "protein groups"):
                    self.assertNotIn("precursor", r["loss"].lower())
                    self.assertNotIn("Precursors", r["cols"])
        self.assertIn("re-adapt it with run_search.py --adapt-only",
                      self.census("NULL", "FALSE", "NULL")["loss"])

    def test_one_item_is_singular(self):
        expr = (f'.script_dir <- "{SCRIPTS}"; source(file.path(.script_dir, "contaminants.R")); '
                f'u <- contaminant_unit("protein"); '
                f'cat(count_of(1, u), count_of(2, u), count_of(1234, u), count_of(NULL, u), '
                f'count_verb(1), count_verb(2), sep = "|")')
        p = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertEqual(p.stdout, "1 protein group|2 protein groups|1,234 protein groups|"
                                   "____ protein groups|was|were")
        import make_methods
        c = {"unit": "peptides", "unit_singular": "peptide"}
        self.assertEqual(make_methods._count_of(1, c), ("1 peptide", "was"))
        self.assertEqual(make_methods._count_of(686, c), ("686 peptides", "were"))
        self.assertEqual(make_methods._count_of(None, {}), ("____ precursors", "were"))


@unittest.skipUnless(r_has("jsonlite"), "needs R with jsonlite")
class DiannContQuantExcludeIsRead(unittest.TestCase):
    """sage-review N1: what DIA-NN's --cont-quant-exclude left out is said only when the run
    recorded the flag -- the command line in the DIA-NN log beside the report, else the parameters
    the search ran with (make_methods.diann_cont_quant_exclude, the one reader) -- and NOT
    RECORDED otherwise: FragPipe's DIA-NN step, or a DIA-NN run by hand, may not have had it."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.report = os.path.join(self.d, "report.parquet")
        open(self.report, "w").close()

    def log(self, flags):
        from test_keratin_sample import DIANN_LOG
        with open(os.path.join(self.d, "report.log.txt"), "w") as fh:
            fh.write(DIANN_LOG.format(out=self.report, flags=flags))

    def kept(self, requantified):
        expr = (f'.script_dir <- "{SCRIPTS}"; source(file.path(.script_dir, "contaminants.R")); '
                f'k <- diann_kept_quant(diann_cont_quant_exclude("{self.report}"), '
                f'{"TRUE" if requantified else "FALSE"}); '
                f'cat(jsonlite::toJSON(k, auto_unbox = TRUE, na = "null"))')
        p = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        return json.loads(p.stdout)

    def test_the_log_says_it_ran(self):
        import make_methods
        self.log(" --cont-quant-exclude Cont_")
        self.assertEqual(make_methods.diann_cont_quant_exclude(self.report)["value"], "Cont_")
        k = self.kept(False)
        self.assertIn("DIA-NN --cont-quant-exclude Cont_, per the DIA-NN command line in "
                      "report.log.txt", k["text"])
        self.assertIs(k["under_quantified"], True)
        self.assertIs(self.kept(True)["under_quantified"], False)

    def test_the_log_says_it_did_not(self):
        import make_methods
        self.log("")
        rec = make_methods.diann_cont_quant_exclude(self.report)
        self.assertIsNone(rec["value"])
        k = self.kept(False)
        self.assertIn("DIA-NN was given no --cont-quant-exclude", k["text"])
        self.assertIs(k["under_quantified"], False)
        self.assertIn("was not set", make_methods.diann_contaminant_sentence(
            {"engine": "diann", "cont_quant_exclude": rec}))

    def test_the_log_wins_over_the_parameters_file(self):
        import make_methods
        cfg = os.path.join(self.d, "params.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--cont-quant-exclude Cont_\n")
        with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "resolved_params_file": cfg}, fh)
        self.assertEqual(make_methods.diann_cont_quant_exclude(self.report),
                         {"value": "Cont_", "source": "params.cfg"})
        self.log("")                              # what ran had no flag: the log is what ran
        self.assertIsNone(make_methods.diann_cont_quant_exclude(self.report)["value"])

    def cfg(self, path, text):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w") as fh:
            fh.write(text)
        return path

    def raw_log(self, line, log=None):
        log = log or os.path.join(self.d, "report.log.txt")
        os.makedirs(os.path.dirname(log), exist_ok=True)
        with open(log, "w") as fh:
            fh.write("\nDIA-NN 2.7.0 Academia  (Data-Independent Acquisition by Neural Networks)\n"
                     "Compiled on Sep 15 2026 15:30:01\nLogical CPU cores: 64\n" + line + "\n\n")

    def test_a_cfg_on_the_command_line_is_read(self):
        """sage-review, ca9aff2: DIA-NN logs `--cfg <file>` UNEXPANDED (2.7.0 PROT_0002, 2.6.1
        taha_dog), so a --one-step search (`diann --cfg params.cfg ...`) read "not set" although
        params.cfg had the flag -- in the DE, the Methods and the run record."""
        import make_methods
        cfg = self.cfg(os.path.join(self.d, "wf", "params.cfg"),
                       "--qvalue 0.01\n--cont-quant-exclude Cont_\n--cut K*,R*\n")
        self.raw_log(f"/x/diann-linux --cfg {cfg} --f /d/a.raw --fasta /db.fasta "
                     f"--out {self.report} --threads 16")
        rec = make_methods.diann_cont_quant_exclude(self.report)
        self.assertEqual(rec, {"value": "Cont_", "source": "the DIA-NN command line in "
                               "report.log.txt and the --cfg file it names (params.cfg)"})
        k = self.kept(False)
        self.assertIn("DIA-NN --cont-quant-exclude Cont_, per the DIA-NN command line in "
                      "report.log.txt and the --cfg file it names (params.cfg)", k["text"])
        self.assertIs(k["under_quantified"], True)
        # a cfg without it, read: "not set" is then what ran
        self.cfg(cfg, "--qvalue 0.01\n")
        self.assertIsNone(make_methods.diann_cont_quant_exclude(self.report)["value"])
        # the flag later on the command line than the cfg still counts (DIA-NN: last wins)
        self.raw_log(f"/x/diann-linux --cfg {cfg} --out {self.report} --cont-quant-exclude Cont_")
        self.assertEqual(make_methods.diann_cont_quant_exclude(self.report)["value"], "Cont_")

    def test_a_relative_cfg_resolves_where_diann_ran(self):
        """FragPipe's DIA route runs DIA-NN in its workdir with --out dia-quant-output/report.tsv
        (HIVE, fragpipe24 DIA-NN 1.8.1: `--cfg <workdir>/filelist_diann.txt-- `); a relative --cfg
        resolves there, not beside the log."""
        import make_methods
        wd = os.path.join(self.d, "fp")
        report = os.path.join(wd, "dia-quant-output", "report.tsv")
        self.cfg(os.path.join(wd, "filelist_diann.txt"),
                 "--f /d/a.mzML\n--cont-quant-exclude contam_\n")
        self.raw_log("/fp/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 --lib library.tsv "
                     "--threads 32 --verbose 1 --out dia-quant-output/report.tsv --qvalue 0.01 "
                     "--matrices --no-prot-inf --cfg filelist_diann.txt-- ",
                     log=os.path.join(wd, "dia-quant-output", "report.log.txt"))
        self.assertEqual(make_methods.diann_cont_quant_exclude(report)["value"], "contam_")
        # the real FragPipe shape: an absolute --cfg with "--" glued on, no flag in it
        self.cfg(os.path.join(wd, "filelist_diann.txt"), "--f /d/a.mzML\n")
        self.raw_log(f"/fp/diann-1.8.1.8 --lib library.tsv --out dia-quant-output/report.tsv "
                     f"--report-lib-info --cfg {wd}/filelist_diann.txt-- ",
                     log=os.path.join(wd, "dia-quant-output", "report.log.txt"))
        rec = make_methods.diann_cont_quant_exclude(report)
        self.assertIsNone(rec["value"])
        self.assertIn("(filelist_diann.txt)", rec["source"])

    def test_a_cfg_that_cannot_be_read_is_never_not_set(self):
        """The flag may be in it: fall back to the parameters file the search ran with, else
        NOT RECORDED -- never "DIA-NN was given no --cont-quant-exclude"."""
        import make_methods
        self.raw_log(f"/x/diann-linux --cfg /gone/params.cfg --out {self.report}")
        self.assertIsNone(make_methods.diann_cont_quant_exclude(self.report))
        k = self.kept(False)
        self.assertIn("[not recorded -- confirm]", k["text"])
        self.assertIsNone(k["under_quantified"])
        cfg = self.cfg(os.path.join(self.d, "params.cfg"), "--cont-quant-exclude Cont_\n")
        with open(os.path.join(self.d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "params_file": cfg}, fh)
        self.assertEqual(make_methods.diann_cont_quant_exclude(self.report),
                         {"value": "Cont_", "source": "params.cfg"})

    def test_search_record_reads_the_cfg_the_log_names(self):
        """The database-search paragraph comes from the same reader (search_record)."""
        import make_methods
        cfg = self.cfg(os.path.join(self.d, "params.cfg"),
                       "--qvalue 0.01\n--cont-quant-exclude Cont_\n")
        self.raw_log(f"/x/diann-linux --cfg {cfg} --f /d/a.raw --out {self.report}")
        prov = os.path.join(self.d, "search_provenance.json")
        with open(prov, "w") as fh:
            json.dump({"engine": "diann", "params_file": cfg,
                       "result": {"report": self.report}}, fh)
        rec = make_methods.search_record(search_prov=prov)
        self.assertEqual(rec["cont_quant_exclude"]["value"], "Cont_")
        self.assertIn("were excluded from normalisation",
                      make_methods.diann_contaminant_sentence(rec))

    def test_nothing_recorded_is_said_so(self):
        import make_methods
        self.assertIsNone(make_methods.diann_cont_quant_exclude(self.report))
        k = self.kept(False)
        self.assertIn("[not recorded -- confirm]", k["text"])
        self.assertNotIn("DIA-NN --cont-quant-exclude", k["text"])
        self.assertIsNone(k["under_quantified"])


class KeratinRefusalNamesOnlyWhatTheEngineDoes(unittest.TestCase):
    """run_search.py refuses a keratin sample on a database still holding keratin Cont_ entries.
    Only DIA-NN is given --cont-quant-exclude, so only its refusal says so."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.db = os.path.join(self.d, "db.fasta")
        with open(self.db, "w") as fh:
            fh.write(">sp|P1|ONE_HUMAN One OS=Homo sapiens OX=9606 GN=ONE\nPEPTIDEK\n"
                     ">sp|Cont_P35527|K1C9_HUMAN Keratin, type I cytoskeletal 9 OS=Homo sapiens "
                     "OX=9606 GN=KRT9 PE=1 SV=3\nPEPTIDER\n")

    def refusal(self, engine):
        import run_search
        with self.assertRaises(SystemExit) as cm:
            run_search.keratin_sample_check(self.db, True, engine)
        return str(cm.exception)

    def test_only_diann_is_said_to_exclude(self):
        self.assertIn("DIA-NN --cont-quant-exclude", self.refusal("diann"))
        for engine in ("sage", "fragpipe", "alphadia", "radiant"):
            with self.subTest(engine=engine):
                msg = self.refusal(engine)
                self.assertNotIn("DIA-NN", msg)
                self.assertIn("run_de.R's contaminant filter would remove", msg)


if __name__ == "__main__":
    unittest.main()
