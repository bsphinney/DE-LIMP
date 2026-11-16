#!/usr/bin/env python3
"""make_protein_screen.py: find the protein, check the data, assemble the page data (offline)."""
import gzip
import re
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
FIXTURE = os.path.join(HERE, "fixtures", "screen_data_mini")
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("screen_data_mini_make", os.path.join(FIXTURE, "make.py"))
fixture_make = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fixture_make)

import make_protein_screen as mps  # noqa: E402
from make_protein_screen import (file_labels, page_data, parse_index, peptide_of, top_proteins,  # noqa: E402
                                 resolve_protein, verify)

SCRIPT = os.path.join(SCRIPTS, "make_protein_screen.py")
MAN = json.load(open(os.path.join(FIXTURE, "manifest.json")))
IDX = parse_index(os.path.join(FIXTURE, "index.tsv"))
REC = json.load(gzip.open(os.path.join(FIXTURE, "parts", "part-0001.json.gz"), "rt"))

SEQ = "MK" + "AEEATEAQEVVEATPEGACTEPR" + "GG" + "LGSEQGGGK" + "AAAA"
UNIPROT = {"acc": "P10001", "id": "AAA_HUMAN", "reviewed": "UniProtKB reviewed (Swiss-Prot)", "name": "Protein A",
           "gene": "GENEA", "organism": "Homo sapiens", "length": len(SEQ), "sequence": SEQ,
           "function": "Does a thing.", "location": ["Cytoplasm"], "uniprot_version": 3,
           "uniprot_updated": "2026-09-02", "release": "2026_03", "fetched_at": "2026-10-28"}


def run_cli(*args, screen_data=FIXTURE, out=None):
    out = out or os.path.join(tempfile.mkdtemp(), "page.html")
    return subprocess.run([sys.executable, SCRIPT, "--screen-data", screen_data, "--out", out, *args],
                          capture_output=True, text=True)


def copy_of_fixture():
    d = os.path.join(tempfile.mkdtemp(), "output", "screen_data")
    shutil.copytree(FIXTURE, d, ignore=shutil.ignore_patterns("make.py", "__pycache__"))
    return d


def corrupt_one_byte(screen_dir, rel):
    p = os.path.join(screen_dir, rel)
    b = bytearray(open(p, "rb").read())
    b[len(b) // 2] ^= 0xFF
    open(p, "wb").write(bytes(b))


class TestFind(unittest.TestCase):
    def test_resolve_by_gene_accession_group(self):
        self.assertEqual([r["Protein.Group"] for r in resolve_protein(IDX, "genea")], ["P10001"])
        self.assertEqual([r["Protein.Group"] for r in resolve_protein(IDX, "P10001")], ["P10001"])
        self.assertEqual([r["Protein.Group"] for r in resolve_protein(IDX, "p10004")], ["P10004"])

    def test_resolve_group_with_several_accessions(self):
        idx = [{"Protein.Group": "A0A804HLA8;P63218", "Genes": "GNG5;GNG5B"},
               {"Protein.Group": "P99999", "Genes": "OTHER"}]
        for q in ("GNG5B", "P63218", "a0a804hla8;p63218"):
            self.assertEqual([r["Protein.Group"] for r in resolve_protein(idx, q)], ["A0A804HLA8;P63218"], q)
        self.assertEqual(resolve_protein(idx, "GNG"), [])

    def test_unknown_protein(self):
        self.assertEqual(resolve_protein(IDX, "NOPE"), [])
        self.assertEqual(run_cli("--protein", "NOPE", "--offline").returncode, 2)

    def test_ambiguous_gene_lists_all(self):
        p = run_cli("--protein", "GENEX", "--offline")
        self.assertEqual(p.returncode, 2)
        self.assertIn("P10002", p.stderr)
        self.assertIn("P10003", p.stderr)

    def test_index_is_read_without_csv_quoting(self):
        # make_screen_data.R writes index.tsv with quote = FALSE: a double quote in a name is plain text, and must
        # not swallow the rows after it
        p = os.path.join(tempfile.mkdtemp(), "index.tsv")
        with open(p, "w") as fh:
            fh.write("Protein.Group\tGenes\tProtein.Names\tpart\tn_precursors\n")
            fh.write('P1\tG1\t"odd name\tpart-0001.json.gz\t2\n')
            fh.write("P2\tG2\tNAME2\tpart-0001.json.gz\t1\n")
        rows = parse_index(p)
        self.assertEqual([r["Protein.Group"] for r in rows], ["P1", "P2"])
        self.assertEqual(rows[0]["Protein.Names"], '"odd name')

    def test_peptide_of(self):
        self.assertEqual(peptide_of("AEEATEAQEVVEATPEGAC(UniMod:4)TEPR3"), "AEEATEAQEVVEATPEGACTEPR")
        self.assertEqual(peptide_of("(UniMod:1)MAAK12"), "MAAK")
        self.assertEqual(peptide_of("LGSEQGGGK2"), "LGSEQGGGK")

    def test_the_builder_has_no_network_code_of_its_own(self):
        import glob
        import protein_facts as pf
        self.assertIs(mps.lookup, pf.lookup)
        src = open(SCRIPT).read()
        self.assertNotIn("urllib", src)
        self.assertNotIn("def fetch_uniprot", src)
        defs = [f for f in glob.glob(os.path.join(SCRIPTS, "*.py"))
                if "def fetch_uniprot" in open(f).read()]
        self.assertEqual([os.path.basename(f) for f in defs], ["protein_facts.py"])      # one definition

    def test_the_builder_hashes_with_provenance_py(self):
        # repo rule 3: one sha256 helper, provenance.py's
        import provenance
        self.assertIs(mps.sha256_file, provenance.sha256_file)
        self.assertNotIn("def sha256_file", open(SCRIPT).read())


class TestFactsInPageData(unittest.TestCase):
    FACT = {"source": "NCBI Gene", "status": "ok", "url": "https://x", "release": None, "fetched_at": "2026-10-28",
            "licence": "l", "fields": [["Symbol", "A"]], "note": None}

    def D(self, uniprot, facts):
        return page_data(MAN, REC["P10001"], None, uniprot, index_rows=IDX, facts=facts)

    def test_facts_ride_on_info_when_uniprot_answered(self):
        D = self.D(UNIPROT, [self.FACT])
        self.assertEqual(D["info"]["P10001"]["facts"], [self.FACT])
        self.assertEqual(D["facts_no_uniprot"], {})

    def test_facts_still_reach_the_page_without_uniprot(self):
        D = self.D(None, [self.FACT])
        self.assertEqual(D["info"], {})
        self.assertEqual(D["facts_no_uniprot"], {"P10001": [self.FACT]})


class TestRefuse(unittest.TestCase):
    def test_checksum_mismatch_refuses(self):
        d = copy_of_fixture()
        corrupt_one_byte(d, "parts/part-0001.json.gz")
        p = run_cli("--protein", "GENEA", "--offline", screen_data=d)
        self.assertEqual(p.returncode, 2)
        self.assertIn("part-0001", p.stderr)
        self.assertTrue(verify(d, MAN))

    def test_clean_fixture_verifies(self):
        self.assertEqual(verify(FIXTURE, MAN), [])

    def test_missing_part_is_a_problem(self):
        d = copy_of_fixture()
        os.remove(os.path.join(d, "parts", "part-0001.json.gz"))
        self.assertTrue(any("part-0001" in m for m in verify(d, MAN)))

    def test_no_manifest_and_no_record_says_never_built_not_made_before(self):
        p = run_cli("--protein", "GENEA", "--offline", screen_data=tempfile.mkdtemp())
        self.assertEqual(p.returncode, 2)
        self.assertIn("no run record: the screen data was never built", p.stderr)
        self.assertNotIn("made before", p.stderr)

    def test_no_manifest_says_what_the_last_build_did(self):
        # the folder after a crashed rebuild: make_screen_data.R emptied it (no manifest.json) and wrote its
        # run record beside it; the page builder reports that, never "made before the protein screen"
        cases = [
            ({"state": "failed", "reason": "building the record of P1 failed: its worker returned nothing"},
             ["last build of the screen data failed (building the record of P1 failed: its worker returned nothing)",
              "Run make_screen_data.R on it again"]),
            ({"state": "running"}, ["started and did not finish", "killed", "job log"]),
            ({"state": "refused", "route": "maxlfq", "reason": "a MaxLFQ analysis"},
             ["was refused (maxlfq): a MaxLFQ analysis"]),
            ({"state": "built"}, ["says it was built", "changed or removed"]),
        ]
        for rec, needles in cases:
            out = tempfile.mkdtemp()
            sd = os.path.join(out, "screen_data")
            os.makedirs(os.path.join(sd, "parts"))
            with open(os.path.join(out, "screen_data_status.json"), "w") as fh:
                json.dump(dict({"schema": "screen_data_status/1"}, **rec), fh)
            p = run_cli("--protein", "GENEA", "--offline", screen_data=sd)
            self.assertEqual(p.returncode, 2, p.stderr)
            for n in needles:
                self.assertIn(n, p.stderr, rec)
            self.assertNotIn("made before", p.stderr)

    def test_corrupt_manifest_refuses_with_a_reason(self):
        d = copy_of_fixture()
        for text, why in (('{"schema": "screen_data/1", ', "not valid JSON"), ("[1, 2]", "not a JSON object")):
            open(os.path.join(d, "manifest.json"), "w").write(text)
            p = run_cli("--protein", "GENEA", "--offline", screen_data=d)
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("manifest.json", p.stderr)
            self.assertIn(why, p.stderr)
            self.assertIn("make_screen_data.R", p.stderr)
            self.assertNotIn("Traceback", p.stderr)

    def test_unsupported_refuses_with_reason(self):
        d = copy_of_fixture()
        m = json.load(open(os.path.join(d, "manifest.json")))
        m["supported"] = False
        m["unsupported_reason"] = "this analysis used MaxLFQ"
        json.dump(m, open(os.path.join(d, "manifest.json"), "w"))
        p = run_cli("--protein", "GENEA", "--offline", screen_data=d)
        self.assertEqual(p.returncode, 2)
        self.assertIn("this analysis used MaxLFQ", p.stderr)

    def test_unknown_contrast_lists_the_choices(self):
        p = run_cli("--protein", "GENEA", "--offline", "--contrast", "Z-Y")
        self.assertEqual(p.returncode, 2)
        self.assertIn("B-A", p.stderr)
        self.assertIn("C-A", p.stderr)


    def test_corrupt_index_refuses(self):
        d = copy_of_fixture()
        with open(os.path.join(d, "index.tsv"), "a") as fh:
            fh.write("P10009\tGENEZ\tZZZ_HUMAN\tpart-0001.json.gz\t1\n")
        p = run_cli("--protein", "GENEA", "--offline", screen_data=d)
        self.assertEqual(p.returncode, 2)
        self.assertIn("index.tsv", p.stderr)
        self.assertTrue(any("index.tsv" in m for m in verify(d, MAN)))

    def test_index_part_not_in_manifest_refuses(self):
        d = copy_of_fixture()
        shutil.copy(os.path.join(d, "parts", "part-0001.json.gz"), os.path.join(d, "parts", "part-0002.json.gz"))
        m = json.load(open(os.path.join(d, "manifest.json")))
        text = open(os.path.join(d, "index.tsv")).read().replace("part-0001", "part-0002")
        open(os.path.join(d, "index.tsv"), "w").write(text)
        m["parts"]["index.tsv"] = hashlib.sha256(text.encode()).hexdigest()
        json.dump(m, open(os.path.join(d, "manifest.json"), "w"))
        p = run_cli("--protein", "GENEA", "--offline", screen_data=d)
        self.assertEqual(p.returncode, 2)
        self.assertIn("part-0002", p.stderr)
        self.assertNotIn("Traceback", p.stderr)

    def test_protein_missing_from_its_part_refuses(self):
        d = copy_of_fixture()
        part = os.path.join(d, "parts", "part-0001.json.gz")
        recs = json.load(gzip.open(part, "rt"))
        del recs["P10001"]
        with gzip.open(part, "wt") as fh:
            json.dump(recs, fh)
        m = json.load(open(os.path.join(d, "manifest.json")))
        m["parts"]["parts/part-0001.json.gz"] = hashlib.sha256(open(part, "rb").read()).hexdigest()
        json.dump(m, open(os.path.join(d, "manifest.json"), "w"))
        p = run_cli("--protein", "GENEA", "--offline", screen_data=d)
        self.assertEqual(p.returncode, 2)
        self.assertIn("P10001", p.stderr)
        self.assertNotIn("Traceback", p.stderr)


    def test_a_record_filed_under_another_key_refuses(self):
        # a part whose key P10001 holds P10003's record (what a killed worker's misaligned records would file)
        d = copy_of_fixture()
        part = os.path.join(d, "parts", "part-0001.json.gz")
        recs = json.load(gzip.open(part, "rt"))
        recs["P10001"] = recs["P10003"]
        with gzip.open(part, "wt") as fh:
            json.dump(recs, fh)
        m = json.load(open(os.path.join(d, "manifest.json")))
        m["parts"]["parts/part-0001.json.gz"] = hashlib.sha256(open(part, "rb").read()).hexdigest()
        json.dump(m, open(os.path.join(d, "manifest.json"), "w"))
        p = run_cli("--protein", "GENEA", "--offline", screen_data=d)
        self.assertEqual(p.returncode, 2, p.stderr)
        self.assertIn("P10001", p.stderr)
        self.assertIn("P10003", p.stderr)
        self.assertNotIn("Traceback", p.stderr)


class TestPartial(unittest.TestCase):
    """A folder make_screen_data.R --proteins wrote (manifest partial: true) holds only some proteins."""

    def setUp(self):
        self.d = copy_of_fixture()
        m = json.load(open(os.path.join(self.d, "manifest.json")))
        m["partial"] = True
        json.dump(m, open(os.path.join(self.d, "manifest.json"), "w"))

    def test_top_refuses_a_partial_folder(self):
        p = run_top("--top", "5", "--offline", screen_data=self.d)
        self.assertEqual(p.returncode, 2, p.stderr)
        self.assertIn("partial", p.stderr)
        self.assertIn("--proteins", p.stderr)
        self.assertNotIn("Traceback", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(p.out_dir, "index.html")))

    def test_no_match_names_the_folder_as_partial(self):
        p = run_cli("--protein", "NOPE", "--offline", screen_data=self.d)
        self.assertEqual(p.returncode, 2, p.stderr)
        self.assertIn("no protein matches 'NOPE'", p.stderr)
        self.assertIn("partial", p.stderr)
        self.assertIn("4 proteins", p.stderr)
        full = run_cli("--protein", "NOPE", "--offline")
        self.assertNotIn("partial", full.stderr)

    def test_a_protein_in_a_partial_folder_still_builds(self):
        self.assertEqual(run_cli("--protein", "GENEA", "--offline", screen_data=self.d).returncode, 0)


class TestSources(unittest.TestCase):
    def setUp(self):
        self.d = copy_of_fixture()
        self.table = os.path.join(os.path.dirname(self.d), "tables", "DE_dpc_B.A.csv")
        os.makedirs(os.path.dirname(self.table))

    def test_absent_source_is_not_checked(self):
        verified, not_checked, bad = mps.check_sources(MAN, os.path.dirname(self.d))
        self.assertEqual((verified, not_checked, bad), ([], ["tables/DE_dpc_B.A.csv"], []))

    def test_matching_source_is_verified(self):
        open(self.table, "w").write(fixture_make.DE_TABLE_TEXT)
        verified, not_checked, bad = mps.check_sources(MAN, os.path.dirname(self.d))
        self.assertEqual((verified, not_checked, bad), (["tables/DE_dpc_B.A.csv"], [], []))
        self.assertEqual(verify(self.d, MAN), [])

    def test_each_source_is_hashed_once_per_build(self):
        # report.parquet can be gigabytes: one page, or one --top list with its pages, hashes it once
        open(self.table, "w").write(fixture_make.DE_TABLE_TEXT)
        seen, real = [], mps.sha256_file

        def counting(path):
            seen.append(os.path.abspath(path))
            return real(path)
        mps.sha256_file = counting
        try:
            out = os.path.join(tempfile.mkdtemp(), "p.html")
            self.assertEqual(mps.main(["--screen-data", self.d, "--protein", "GENEA", "--offline", "--out", out]), 0)
            self.assertEqual(seen.count(os.path.abspath(self.table)), 1)
            del seen[:]
            top = os.path.join(tempfile.mkdtemp(), "top")
            self.assertEqual(mps.main(["--screen-data", self.d, "--top", "5", "--offline", "--out-dir", top]), 0)
            self.assertEqual(seen.count(os.path.abspath(self.table)), 1)
        finally:
            mps.sha256_file = real

    def test_changed_source_refuses(self):
        open(self.table, "w").write(fixture_make.DE_TABLE_TEXT + "edited\n")
        p = run_cli("--protein", "GENEA", "--offline", screen_data=self.d)
        self.assertEqual(p.returncode, 2)
        self.assertIn("the analysis changed after its screen data was made; rebuild it with make_screen_data.R",
                      p.stderr)


class TestPageData(unittest.TestCase):
    def D(self, pg="P10001", contrast="B-A", uniprot=None, **kw):
        return page_data(MAN, REC[pg], contrast, uniprot, index_rows=IDX, **kw)

    def test_untested_contrast_page_data(self):
        D = page_data(MAN, REC["P10004"], "B-A", None)
        self.assertFalse(D["proteins"]["P10004"]["tested"])
        self.assertIsNone(D["proteins"]["P10004"]["logFC"])
        json.dumps(D, allow_nan=False)
        D = page_data(MAN, REC["P10004"], "C-A", None)
        self.assertTrue(D["proteins"]["P10004"]["tested"])

    def test_offline_marks_uniprot_missing(self):
        D = page_data(MAN, REC["P10001"], "B-A", None)
        self.assertIn("uniprot", D["missing"])
        self.assertNotIn("P10001", D["info"])

    def test_keys_and_values(self):
        D = self.D()
        for k in ("samples", "groups", "dpc", "prior", "n_prec", "n_tested", "n_sig", "proteins", "rows",
                  "contrib", "detprob", "routes", "stats", "info", "xic", "byrow", "contrast", "files",
                  "analysis", "missing", "provenance"):
            self.assertIn(k, D)
        self.assertEqual(D["contrast"], {"name": "B-A", "num": "B", "den": "A", "adjp": 0.05, "default": False})
        self.assertEqual((D["n_tested"], D["n_sig"], D["n_prec"]), (3, 2, 5))
        self.assertEqual(D["files"]["de_table"], "tables/DE_dpc_B.A.csv")
        self.assertEqual(D["files"]["session_rds"], "tables/DE-LIMP_session.rds")
        self.assertEqual(D["analysis"], MAN["analysis"])
        p = D["proteins"]["P10001"]
        self.assertEqual((p["gene"], p["name"], p["tested"], p["det"]), ("GENEA", "AAA_HUMAN", True,
                                                                         {"A": "2/2", "B": "2/2"}))
        self.assertEqual(p["CI"], [REC["P10001"]["de"]["B-A"]["CI.L"], REC["P10001"]["de"]["B-A"]["CI.R"]])
        self.assertEqual(p["prec_ids"], [x["id"] for x in REC["P10001"]["precursors"]])
        # limpa's count of the protein's precursors (genes$NPrec, NPeptides before limpa 1.4), under that name
        self.assertEqual(p["nprec"], REC["P10001"]["n_precursors"])
        self.assertEqual(p["nprec"], 2)
        self.assertNotIn("npep", p)
        self.assertEqual(len(D["rows"]), len(REC["P10001"]["rows"]))
        self.assertEqual({r["Protein.Group"] for r in D["rows"]}, {"P10001"})
        self.assertEqual(D["xic"], {})
        self.assertEqual(D["byrow"], {})

    def test_default_contrast_is_the_manifests_first(self):
        self.assertEqual(self.D(contrast=None)["contrast"]["name"], "B-A")

    def test_a_default_contrast_is_flagged_when_there_was_a_choice(self):
        # repo rule 2: the page tags a comparison nobody chose
        self.assertTrue(self.D(contrast=None)["contrast"]["default"])
        self.assertFalse(self.D(contrast="B-A")["contrast"]["default"])
        one = dict(MAN, contrasts=MAN["contrasts"][:1])
        self.assertFalse(page_data(one, REC["P10001"], None, None, index_rows=IDX)["contrast"]["default"])

    def test_top_pages_and_list_carry_the_default_flag(self):
        for extra, flag in (([], True), (["--contrast", "B-A"], False)):
            p = run_top("--top", "5", "--offline", *extra)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIs(list_data(read(p.out_dir, "index.html"))["contrast"]["default"], flag)
            page = read(p.out_dir, "GENEA.html")
            D = json.loads(re.search(r"const D = (\{.*?\});\n", page, re.S).group(1))
            self.assertIs(D["contrast"]["default"], flag)

    def test_blocks_map_from_the_record(self):
        D = self.D()
        pg = "P10001"
        self.assertEqual(D["contrib"][pg]["ids"], D["proteins"][pg]["prec_ids"])
        self.assertEqual(D["contrib"][pg]["own"], REC[pg]["shares"]["own"])
        self.assertEqual(D["contrib"][pg]["total"], REC[pg]["shares"]["total"])
        self.assertEqual(D["detprob"][pg]["pf_fit"], REC[pg]["detection"]["pf_fit"])
        self.assertEqual(D["detprob"][pg]["ids"], D["proteins"][pg]["prec_ids"])
        self.assertEqual(D["routes"][pg], {"w": REC[pg]["weights"]["w"],
                                           "wshare": REC[pg]["weights"]["share_in_group"], "df_res": 3,
                                           "df_prior": 8.5})
        self.assertNotIn("shares", D["missing"])
        self.assertNotIn("detection", D["missing"])

    def test_page_data_for_the_template(self):
        """What the page reads beyond Task 4's D: the comparison's cutoff, limma's prior df, the shared
        peptides (from the report rows, so an offline page has them too), how the shares were made, and
        the analysis's own q-value columns and cutoffs (de_provenance.json, via the manifest)."""
        D = self.D()
        self.assertEqual(D["contrast"]["adjp"], 0.05)
        self.assertEqual(D["routes"]["P10001"]["df_prior"], 8.5)
        self.assertEqual(D["shared"]["P10001"]["AEEATEAQEVVEATPEGACTEPR"]["others"],
                         [{"acc": "P10003", "gene": "GENEX"}])
        self.assertEqual(D["shared"]["P10001"], self.D(uniprot=UNIPROT)["info"]["P10001"]["shared"])
        sm = D["share_method"]
        self.assertEqual((sm["method"], sm["limpa_version"], sm["tolerance"]), ("derivative", "1.4.0", 1e-06))
        self.assertEqual(sm["check"], MAN["shares"]["nudge_check"])
        self.assertNotIn("withheld", sm)
        self.assertEqual(D["analysis"]["q_columns"], ["Q.Value", "Lib.Q.Value"])
        self.assertEqual(D["analysis"]["q_cutoffs"], [0.01, 0.01])
        # the page keys peptides by peptide_of (repo rule 3: one definition), never by its own regex
        self.assertEqual(D["peptides"], {"AEEATEAQEVVEATPEGAC(UniMod:4)TEPR3": "AEEATEAQEVVEATPEGACTEPR",
                                         "LGSEQGGGK2": "LGSEQGGGK"})

    def test_isoform_group_has_no_coverage_and_says_why(self):
        # UniProt's entry holds the canonical sequence: positions on it would not be the isoform's
        rec = json.loads(json.dumps(REC["P10001"]))
        rec["pg"] = "P10001-2"
        D = page_data(MAN, rec, "B-A", UNIPROT, index_rows=IDX)
        i = D["info"]["P10001-2"]
        self.assertIsNone(i["cov"])
        self.assertIsNone(i["pos"])
        self.assertIn("the isoform P10001-2", D["missing"]["coverage"])
        self.assertIn("canonical sequence of P10001", D["missing"]["coverage"])
        self.assertEqual(D["lookup"], {"accession": "P10001-2", "accessions": ["P10001-2"]})
        self.assertNotIn("coverage", self.D(uniprot=UNIPROT)["missing"])         # a canonical accession keeps it

    def test_multi_accession_group_records_which_accession_was_looked_up(self):
        rec = json.loads(json.dumps(REC["P10001"]))
        rec["pg"] = "P10001;P10009"
        D = page_data(MAN, rec, "B-A", UNIPROT, index_rows=IDX)
        self.assertEqual(D["lookup"], {"accession": "P10001", "accessions": ["P10001", "P10009"]})
        self.assertTrue(D["info"]["P10001;P10009"]["pos"])                       # the first is canonical
        self.assertNotIn("coverage", D["missing"])

    def test_a_dashed_accession_that_is_no_isoform_keeps_its_coverage(self):
        # protein_ids.py's isoform rule (-<digits>), not "anything after a dash"
        import protein_ids
        self.assertIs(mps.strip_isoform, protein_ids.strip_isoform)
        rec = json.loads(json.dumps(REC["P10001"]))
        rec["pg"] = "ABC123-ALT"
        D = page_data(MAN, rec, "B-A", UNIPROT, index_rows=IDX)
        self.assertTrue(D["info"]["ABC123-ALT"]["pos"])
        self.assertNotIn("coverage", D["missing"])
        self.assertIsNone(mps._fact_note("ABC123-ALT"))
        self.assertIsNotNone(mps._isoform_reason("P10001-2"))                    # an isoform still is one
        self.assertIn("isoform P10001-2", mps._fact_note("P10001-2"))

    def test_list_row_says_which_entry_names_a_group(self):
        man = json.load(open(os.path.join(FIXTURE, "manifest.json")))
        fact = [{"source": "UniProt", "status": "ok", "fields": [["Accession", "P10001"], ["Protein name", "Name"]]}]

        def row_for(pg):
            rec = json.loads(json.dumps(REC["P10001"]))
            rec["pg"] = pg
            r = dict(IDX[0], **{"Protein.Group": pg})
            return mps.list_data(man, man["contrasts"][0], 0.05, "both", 5, [r], 1, {pg: rec}, IDX, {pg: fact},
                                 {pg: "GENEA"}, {})["proteins"][0]
        self.assertIsNone(row_for("P10001")["fact_note"])
        self.assertEqual(row_for("P10001;P10009")["fact_note"], "UniProt entry of P10001, the group's first accession")
        self.assertEqual(row_for("P10001-2")["fact_note"],
                         "UniProt's canonical entry P10001; the group's accession is the isoform P10001-2")

    def test_withheld_shares_go_to_missing(self):
        D = self.D("P10004", "C-A")
        self.assertNotIn("P10004", D["contrib"])
        self.assertNotIn("P10004", D["detprob"])
        self.assertNotIn("P10004", D["routes"])
        self.assertIn("5.15e-05", D["missing"]["shares"])
        self.assertIn("detection", D["missing"])
        self.assertIn("weights", D["missing"])

    def test_manifest_skipped_blocks_are_missing(self):
        D = self.D()
        self.assertIn("does not exist", D["missing"]["run_stats"])
        self.assertIn("xic", D["missing"])
        self.assertEqual(D["stats"], {})

    def test_no_page_text_describes_a_nudge_of_0_01(self):
        text = json.dumps(self.D(uniprot=UNIPROT)).lower()
        self.assertNotIn("nudg", text)
        self.assertNotIn("0.01 ", text.replace('"0.01"', ""))

    def test_uniprot_info(self):
        D = self.D(uniprot=UNIPROT)
        self.assertNotIn("uniprot", D["missing"])
        i = D["info"]["P10001"]
        self.assertNotIn("fran", i)
        self.assertEqual((i["acc"], i["id"], i["length"], i["fetched_at"]), ("P10001", "AAA_HUMAN", len(SEQ),
                                                                             "2026-10-28"))
        pep = "AEEATEAQEVVEATPEGACTEPR"
        self.assertEqual(i["pos"][pep], [3, 3 + len(pep) - 1])
        self.assertEqual(i["pos"]["LGSEQGGGK"], [SEQ.find("LGSEQGGGK") + 1, SEQ.find("LGSEQGGGK") + 9])
        self.assertEqual(sum(i["cov"]), len(pep) + 9)
        self.assertEqual(len(i["cov"]), len(SEQ))

    def test_shared_peptides_named_from_the_index(self):
        i = self.D(uniprot=UNIPROT)["info"]["P10001"]
        self.assertEqual(i["shared"]["AEEATEAQEVVEATPEGACTEPR"]["others"], [{"acc": "P10003", "gene": "GENEX"}])
        self.assertEqual(i["shared"]["LGSEQGGGK"]["others"], [{"acc": "P10009", "gene": None}])
        self.assertEqual(i["shared"]["LGSEQGGGK"]["proteotypic"], 1)

    def test_peptide_not_in_sequence_has_no_position(self):
        u = dict(UNIPROT, sequence="MKKKKKKKKK", length=10)
        i = self.D(uniprot=u)["info"]["P10001"]
        self.assertIsNone(i["pos"]["LGSEQGGGK"])
        self.assertEqual(sum(i["cov"]), 0)

    def test_provenance(self):
        D = self.D(verified=["tables/x.csv"], not_checked=["search/report.parquet"])
        self.assertEqual(D["provenance"]["verified"], ["tables/x.csv"])
        self.assertEqual(D["provenance"]["not_checked"], ["search/report.parquet"])
        self.assertEqual(D["provenance"]["written_by"], MAN["written_by"])


class TestRender(unittest.TestCase):
    def test_embedded_data_survives_json_loads(self):
        D = {"info": {"function": "x </script><script>alert(1)</script> <!-- y --> \u2028 &"}, "n": 1}
        page = mps.render("<html><script>const D = __SCREEN_DATA__;</script></html>", D)
        self.assertEqual(page.count("</script>"), 1)
        self.assertNotIn("<!--", page)
        block = re.search(r"const D = (.*);</script>", page, re.S).group(1)
        self.assertEqual(json.loads(block), D)


class TestCli(unittest.TestCase):
    def test_dump_data_offline(self):
        dump = os.path.join(tempfile.mkdtemp(), "D.json")
        p = run_cli("--protein", "GENEA", "--offline", "--dump-data", dump)
        self.assertEqual(p.returncode, 0, p.stderr)
        D = json.load(open(dump))
        self.assertIn("P10001", D["proteins"])
        self.assertIn("uniprot", D["missing"])
        self.assertEqual(D["contrast"]["name"], "B-A")
        facts = D["facts_no_uniprot"]["P10001"]
        self.assertEqual([f["source"] for f in facts],
                         ["UniProt", "NCBI Gene", "HGNC", "Ensembl", "Human Protein Atlas"])
        self.assertTrue(all(f["status"] == "not_reached" and f["note"] == "--offline" for f in facts))

    def test_out_into_a_folder_not_made_yet(self):
        # SKILL.md's HIVE route writes <S>/output/protein_screens/<gene>.html before that folder exists
        out = os.path.join(tempfile.mkdtemp(), "protein_screens", "GENEA.html")
        p = run_cli("--protein", "GENEA", "--offline", out=out)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertTrue(os.path.isfile(out))

    def test_offline_never_calls_fetch(self):
        called = []
        saved = mps.lookup
        mps.lookup = lambda *a, **k: called.append(a)
        try:
            out = os.path.join(tempfile.mkdtemp(), "D.json")
            rc = mps.main(["--screen-data", FIXTURE, "--protein", "P10001", "--offline", "--out",
                           out + ".html", "--dump-data", out])
        finally:
            mps.lookup = saved
        self.assertEqual((rc, called), (0, []))

    def test_failed_fetch_is_missing_not_fatal(self):
        saved = mps.lookup

        def failing(acc, errors=None):
            errors.append("HTTP 500")
            return None, [{"source": "UniProt", "status": "not_reached", "url": None, "release": None,
                           "fetched_at": "2026-10-28", "licence": None, "fields": [], "note": "HTTP 500"}]
        mps.lookup = failing
        try:
            out = os.path.join(tempfile.mkdtemp(), "D.json")
            rc = mps.main(["--screen-data", FIXTURE, "--protein", "P10001", "--out", out + ".html",
                           "--dump-data", out])
        finally:
            mps.lookup = saved
        self.assertEqual(rc, 0)
        D = json.load(open(out))
        self.assertIn("HTTP 500", D["missing"]["uniprot"])

    def test_fetch_uses_first_accession_of_a_group(self):
        seen = []
        saved = mps.lookup
        mps.lookup = lambda acc, errors=None: seen.append(acc) or (None, [])
        try:
            out = os.path.join(tempfile.mkdtemp(), "D.json")
            mps.main(["--screen-data", FIXTURE, "--protein", "P10001", "--out", out + ".html", "--dump-data", out])
        finally:
            mps.lookup = saved
        self.assertEqual(seen, ["P10001"])


def run_top(*args, screen_data=FIXTURE, out_dir=None):
    out_dir = out_dir or os.path.join(tempfile.mkdtemp(), "top")
    p = subprocess.run([sys.executable, SCRIPT, "--screen-data", screen_data, "--out-dir", out_dir, *args],
                       capture_output=True, text=True)
    p.out_dir = out_dir
    return p


def read(d, name):
    with open(os.path.join(d, name), encoding="utf-8") as fh:
        return fh.read()


def list_data(html):
    """The D the list page was built with (the JSON at the data token), parsed back."""
    m = re.search(r"const D = (\{.*?\});\nconst CT", html, re.S)
    return json.loads(m.group(1))


class TestTop(unittest.TestCase):
    # fixture, B-A: P10001 p .001 adj .03 up; P10003 p .004 adj .04 down; P10002 p .01 adj .05 (NOT below
    # .05); P10004 untested. C-A: P10004 p .0005 adj .01 up; P10001 p .02 adj .04 down; others adj .9.
    def test_top_ranks_by_p_value_called_only(self):
        rows, k = top_proteins(IDX, "B-A", 10, "both", 0.05)
        self.assertEqual([r["Protein.Group"] for r in rows], ["P10001", "P10003"])
        self.assertEqual(k, 2)
        self.assertTrue(all(float(r["adj.P.Val__B-A"]) < 0.05 for r in rows))
        rows, k = top_proteins(IDX, "C-A", 10, "both", 0.05)
        self.assertEqual([r["Protein.Group"] for r in rows], ["P10004", "P10001"])

    def test_the_cut_is_strictly_below_adjp(self):
        # P10002 sits at adj.P.Val == 0.05: not called; at 0.051 cutoff it is
        self.assertNotIn("P10002", [r["Protein.Group"] for r in top_proteins(IDX, "B-A", 10, "both", 0.05)[0]])
        self.assertIn("P10002", [r["Protein.Group"] for r in top_proteins(IDX, "B-A", 10, "both", 0.0501)[0]])

    def test_never_ranks_by_fold_change_or_fills_with_uncalled(self):
        rows, k = top_proteins(IDX, "C-A", 10, "both", 0.05)
        self.assertEqual(k, 2)                                   # P10002/P10003 are tested but uncalled
        self.assertEqual(len(rows), 2)
        idx = [dict(r) for r in IDX]
        for r in idx:                                            # a huge fold change with a poor adjusted p
            if r["Protein.Group"] == "P10002":
                r["logFC__B-A"], r["P.Value__B-A"], r["adj.P.Val__B-A"] = "9", "0.0000001", "0.5"
        rows, _ = top_proteins(idx, "B-A", 10, "both", 0.05)
        self.assertEqual([r["Protein.Group"] for r in rows], ["P10001", "P10003"])

    def test_a_bigger_fold_change_does_not_rank_higher(self):
        idx = [dict(r) for r in IDX]
        for r in idx:
            if r["Protein.Group"] == "P10003":
                r["logFC__B-A"] = "-6"                           # 64-fold, but its p-value is the larger
        rows, _ = top_proteins(idx, "B-A", 10, "both", 0.05)
        self.assertEqual([r["Protein.Group"] for r in rows], ["P10001", "P10003"])

    def test_tied_p_values_have_a_stable_order(self):
        idx = [dict(r) for r in IDX]
        for r in idx:
            r["P.Value__B-A"] = "0.001"
        rows, _ = top_proteins(list(reversed(idx)), "B-A", 10, "both", 0.05)
        self.assertEqual([r["Protein.Group"] for r in rows], ["P10001", "P10003"])

    def test_adjusted_p_ties_break_by_p_value_then_accession(self):
        # BH gives these three one adjusted p; the order is then P.Value, then accession. P10003 has the
        # smallest p but the largest accession; P10001 and P10002 share a p, so accession decides.
        idx = [dict(r) for r in IDX]
        for r in idx:
            if r["Protein.Group"] in ("P10001", "P10002", "P10003"):
                r["adj.P.Val__B-A"] = "0.04"
                r["P.Value__B-A"] = {"P10001": "0.002", "P10002": "0.002", "P10003": "0.001"}[r["Protein.Group"]]
                r["logFC__B-A"] = "1"
        for order in (idx, list(reversed(idx))):
            rows, k = top_proteins(order, "B-A", 10, "both", 0.05)
            self.assertEqual([r["Protein.Group"] for r in rows], ["P10003", "P10001", "P10002"])
            self.assertEqual(k, 3)

    def test_top_direction(self):
        up, k = top_proteins(IDX, "B-A", 10, "up", 0.05)
        self.assertEqual([r["Protein.Group"] for r in up], ["P10001"])
        self.assertTrue(all(float(r["logFC__B-A"]) > 0 for r in up))
        down, _ = top_proteins(IDX, "B-A", 10, "down", 0.05)
        self.assertEqual([r["Protein.Group"] for r in down], ["P10003"])

    def test_n_caps_the_list_and_n_called_counts_all(self):
        rows, k = top_proteins(IDX, "B-A", 1, "both", 0.05)
        self.assertEqual(([r["Protein.Group"] for r in rows], k), (["P10001"], 2))

    def test_no_adjp_refuses_rather_than_assume(self):
        with self.assertRaises(mps.Refused) as c:
            top_proteins(IDX, "B-A", 10, "both", None)
        self.assertIn("adjp", str(c.exception))
        man = json.load(open(os.path.join(FIXTURE, "manifest.json")))
        d = copy_of_fixture()
        for co in man["contrasts"]:
            del co["adjp"]
        json.dump(man, open(os.path.join(d, "manifest.json"), "w"))
        p = run_top("--top", "5", "--offline", screen_data=d)
        self.assertEqual(p.returncode, 2)
        self.assertIn("adjp", p.stderr)

    def test_older_index_without_the_columns_refuses(self):
        old = [{k: v for k, v in r.items() if "__" not in k} for r in IDX]
        with self.assertRaises(mps.Refused) as c:
            top_proteins(old, "B-A", 10, "both", 0.05)
        self.assertIn("make_screen_data.R", str(c.exception))

    def test_fewer_called_than_asked(self):
        p = run_top("--top", "50", "--offline")
        self.assertEqual(p.returncode, 0, p.stderr)
        html = read(p.out_dir, "index.html")
        D = list_data(html)
        self.assertEqual((D["asked"], D["called"], len(D["proteins"])), (50, 2, 2))
        self.assertLess(D["called"], D["asked"])      # the K < N case; the Chrome test checks the rendered sentence
        self.assertEqual(len(D["proteins"]), D["called"])            # exactly the called ones, no filler
        self.assertEqual([r["pg"] for r in D["proteins"]], ["P10001", "P10003"])     # in P-value order
        self.assertEqual(sorted(os.listdir(p.out_dir)), ["GENEA.html", "GENEX.html", "index.html"])
        for bad in ("P10002", "P10004"):                              # tested-but-uncalled, untested
            self.assertNotIn(bad, [r["pg"] for r in D["proteins"]])

    def test_top_cap(self):
        p = run_top("--top", "101", "--offline")
        self.assertEqual(p.returncode, 2)
        self.assertIn("ask for 100 or fewer; the list page links each one", p.stderr)
        self.assertEqual(run_top("--top", "100", "--offline").returncode, 0)
        self.assertEqual(run_top("--top", "0", "--offline").returncode, 2)

    def test_top_and_protein_are_mutually_exclusive(self):
        p = run_top("--top", "3", "--protein", "GENEA", "--offline")
        self.assertEqual(p.returncode, 2)

    def test_top_needs_out_dir_and_refuses_out(self):
        p = subprocess.run([sys.executable, SCRIPT, "--screen-data", FIXTURE, "--top", "3", "--offline"],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 2)
        p = run_top("--top", "3", "--offline", "--out", "x.html")
        self.assertEqual(p.returncode, 2)

    def test_list_links_every_page(self):
        p = run_top("--top", "10", "--offline")
        D = list_data(read(p.out_dir, "index.html"))
        self.assertEqual([r["pg"] for r in D["proteins"]], ["P10001", "P10003"])
        for r in D["proteins"]:
            self.assertTrue(os.path.isfile(os.path.join(p.out_dir, r["link"])), r["link"])
        self.assertEqual(sorted(os.listdir(p.out_dir)), ["GENEA.html", "GENEX.html", "index.html"])

    def test_each_protein_page_is_the_protein_page(self):
        p = run_top("--top", "10", "--offline")
        one = os.path.join(tempfile.mkdtemp(), "one.html")
        self.assertEqual(run_cli("--protein", "P10001", "--offline", out=one).returncode, 0)
        self.assertEqual(read(p.out_dir, "GENEA.html"), read(os.path.dirname(one), "one.html"))

    def test_contrast_and_direction_reach_the_list(self):
        p = run_top("--top", "10", "--offline", "--contrast", "C-A", "--direction", "up")
        D = list_data(read(p.out_dir, "index.html"))
        self.assertEqual((D["contrast"]["name"], D["direction"]), ("C-A", "up"))
        self.assertEqual([r["pg"] for r in D["proteins"]], ["P10004"])
        self.assertEqual(run_top("--top", "3", "--offline", "--contrast", "NOPE").returncode, 2)

    def test_list_row_fields_come_from_the_record(self):
        D = list_data(read(run_top("--top", "10", "--offline").out_dir, "index.html"))
        a, x = D["proteins"]
        self.assertEqual((a["rank"], a["gene"], a["logFC"], a["adjP"], a["P"]), (1, "GENEA", 0.5, 0.03, 0.001))
        self.assertAlmostEqual(a["fold"], 2 ** 0.5)
        self.assertEqual(a["det"], [["A", "2/2"], ["B", "2/2"]])
        self.assertEqual((a["inferred"], a["shared"]), (False, True))       # shares peptides with P10003, P10009
        self.assertEqual((x["inferred"], x["shared"]), (False, False))
        self.assertIsNone(a["fact"])                                          # --offline: nothing was looked up
        self.assertEqual(a["fact_why"], "not reached")

    def test_inferred_means_a_sample_of_the_comparison_with_nothing_found(self):
        man = json.load(open(os.path.join(FIXTURE, "manifest.json")))
        rec = json.loads(json.dumps(REC["P10001"]))
        rec["values"]["n_obs"] = [2, 2, 0, 2, 0, 0]              # b1 (group B) is inferred; c1, c2 are outside B-A
        co = man["contrasts"][0]
        ld = mps.list_data(man, co, 0.05, "both", 5, [IDX[0]], 1, {"P10001": rec}, IDX,
                           {"P10001": mps.offline_facts()}, {"P10001": "GENEA"}, {})
        self.assertTrue(ld["proteins"][0]["inferred"])
        rec["values"]["n_obs"] = [2, 2, 2, 2, 0, 0]
        ld = mps.list_data(man, co, 0.05, "both", 5, [IDX[0]], 1, {"P10001": rec}, IDX,
                           {"P10001": mps.offline_facts()}, {"P10001": "GENEA"}, {})
        self.assertFalse(ld["proteins"][0]["inferred"])

    def test_fact_line_is_the_uniprot_recommended_name_verbatim(self):
        ok = {"source": "UniProt", "status": "ok", "fields": [["Accession", "P10001"],
                                                              ["Protein name", "Brain  acid <b>soluble</b> protein 1"]]}
        self.assertEqual(mps._recommended_name_of([ok, {"source": "HGNC", "status": "ok", "fields": [["Protein name", "x"]]}]),
                         ("Brain  acid <b>soluble</b> protein 1", None))
        self.assertEqual(mps._recommended_name_of([{"source": "UniProt", "status": "not_found"}]), (None, "not in UniProt"))
        self.assertEqual(mps._recommended_name_of([{"source": "UniProt", "status": "not_reached"}]), (None, "not reached"))
        self.assertEqual(mps._recommended_name_of([]), (None, "not reached"))


class TestFactsOncePerAccession(unittest.TestCase):
    def test_each_accession_is_looked_up_once_for_the_list_and_the_pages(self):
        calls, saved = [], mps.lookup

        def fake(acc, errors=None):
            calls.append(acc)
            return None, [{"source": "UniProt", "status": "ok", "url": "u", "release": "r", "fetched_at": "d",
                           "licence": "l", "note": None, "fields": [["Protein name", "Name of " + acc]]}]
        mps.lookup = fake
        try:
            out = os.path.join(tempfile.mkdtemp(), "top")
            rc = mps.main(["--screen-data", FIXTURE, "--top", "10", "--out-dir", out])
        finally:
            mps.lookup = saved
        self.assertEqual(rc, 0)
        self.assertEqual(sorted(calls), ["P10001", "P10003"])            # once each, not once for the list and again for the page
        D = list_data(read(out, "index.html"))
        self.assertEqual([r["fact"] for r in D["proteins"]], ["Name of P10001", "Name of P10003"])
        self.assertEqual([r["fact_src"] for r in D["proteins"]],
                         [{"source": "UniProt", "release": "r", "fetched_at": "d", "url": "u", "acc": None}] * 2)

    def test_a_failed_fetch_marks_that_protein_and_the_list_still_builds(self):
        saved = mps.lookup

        def fake(acc, errors=None):
            if acc == "P10001":
                errors.append("UniProt unreachable")
                return None, [{"source": "UniProt", "status": "not_reached", "url": None, "release": None,
                               "fetched_at": "d", "licence": None, "note": "boom", "fields": []}]
            return None, [{"source": "UniProt", "status": "ok", "url": "u", "release": "r", "fetched_at": "d",
                           "licence": "l", "note": None, "fields": [["Protein name", "Fine"]]}]
        mps.lookup = fake
        try:
            out = os.path.join(tempfile.mkdtemp(), "top")
            rc = mps.main(["--screen-data", FIXTURE, "--top", "10", "--out-dir", out])
        finally:
            mps.lookup = saved
        self.assertEqual(rc, 0)
        D = list_data(read(out, "index.html"))
        self.assertEqual([(r["fact"], r["fact_why"]) for r in D["proteins"]], [(None, "not reached"), ("Fine", None)])


class TestFileLabels(unittest.TestCase):
    def rows(self, *pairs):
        return [{"Protein.Group": pg, "Genes": g} for pg, g in pairs]

    def test_gene_else_leading_accession(self):
        self.assertEqual(file_labels(self.rows(("P1;P2", "ABC;DEF"), ("Q1", ""), ("R1", None))),
                         {"P1;P2": "ABC", "Q1": "Q1", "R1": "R1"})

    def test_safe_for_a_file_name(self):
        lab = file_labels(self.rows(("P1", "../../etc/passwd"), ("P2", "A/B C:D"), ("P3", "..")))
        for v in lab.values():
            self.assertRegex(v, r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
            self.assertNotIn("..", v)
        self.assertEqual(lab["P3"], "P3")

    def test_duplicates_take_the_accession_and_stay_unique(self):
        lab = file_labels(self.rows(("P1", "ABC"), ("P2", "abc"), ("P3", "XYZ"), ("P4", "index")))
        self.assertEqual(lab, {"P1": "ABC_P1", "P2": "abc_P2", "P3": "XYZ", "P4": "index_P4"})
        self.assertEqual(len({v.lower() for v in lab.values()}), 4)

    def test_still_unique_when_the_suffixed_names_collide(self):
        lab = file_labels(self.rows(("P1", "A"), ("P1;X", "A"), ("P9", "A_P1")))
        self.assertEqual(len({v.lower() for v in lab.values()}), 3)


if __name__ == "__main__":
    unittest.main()
