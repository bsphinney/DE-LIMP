#!/usr/bin/env python3
"""The contaminant competition check names the sample proteins the contaminant list took away.

The Core searches every sample with its universal contaminant list. Entries that share peptides
with abundant human and mouse proteins took them away: DIA-NN reported the peptides under the
Cont_ entry, and --cont-quant-exclude Cont_ (and run_de.R's filter) removed them. A mouse tissue
search had no Actb or Eef1a1 at all; a human IP lost 68 genes. The Core director chose (2026-10-20)
to keep the list and have every DIA-NN search report what it takes: contaminant_competition.py,
run by the search's last job, recorded in search_provenance.json and AUDIT.md, and read to the
client in the report's Audit & caveats table. It never blocks.

Stdlib only; the parquet reader's tests are skipped without pyarrow.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import audit_results as ar  # noqa: E402
import contaminant_competition as cc  # noqa: E402
import diann_q_columns as dq  # noqa: E402
import run_search  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:                       # the stdlib-only CI: those tests are skipped
    pa = pq = None

CHECK = os.path.join(SCRIPTS, "contaminant_competition.py")
AUDIT = os.path.join(SCRIPTS, "audit_results.py")

# Synthetic proteins, one tryptic stretch per letter (the search's digest: K*,R*, 7-30 aa).
# ACTB_BOVIN is the human protein residue for residue; ENOA_BOVIN shares two of human ENO1's
# peptides; BSA shares one of human ALB's and has two of its own; the horse entry differs from
# its target only by I/L in one peptide (test_contaminant_evidence's pair).
ACTB = "MEEEEEEEKGGGGGGGKHHHHHHHR"
ACTA1 = "MPPPPPPPKGGGGGGGKHHHHHHHR"          # a paralog sharing two of ACTB's peptides
ENO1 = "MAAAAAAAKWWWWWWWKYYYYYYYR"
ENO1_BOVIN = "MFFFFFFFKWWWWWWWKYYYYYYYR"
ALB = "MTTTTTTTKVVVVVVVKSSSSSSSR"
BSA = "MNNNNNNNKVVVVVVVKQQQQQQQR"
GAPDH = "MDDDDDDDKCCCCCCCR"
TGT = "MSDEQWHYKGHHEAEIKPLAQSHATKASEDLKKFNVWDTAGQEKYLEFISEAIIHVLHSR"
TGT_HORSE = "MCDEQWNYKGHHEAELKPLAQSHATKASEDLKKFNVWDTPGQEKYLEFISDAIIHVLHSR"
TRYP = "IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYK"


def entry(acc, name, gene, seq, cont=False):
    return (f">sp|{'Cont_' if cont else ''}{acc}|{name} x OS=x OX=1 GN={gene} PE=1 SV=1\n"
            f"{seq}\n")


FASTA = "".join([
    entry("P60709", "ACTB_HUMAN", "ACTB", ACTB), entry("P06733", "ENOA_HUMAN", "ENO1", ENO1),
    entry("P68133", "ACTS_HUMAN", "ACTA1", ACTA1),
    entry("P02768", "ALBU_HUMAN", "ALB", ALB), entry("P04406", "G3P_HUMAN", "GAPDH", GAPDH),
    entry("Q90001", "TGT_HUMAN", "TGT", TGT),
    entry("Q80001", "XG_HUMAN", "XG", "MKWQWQWQWQK"),
    entry("Q80002", "YG_HUMAN", "YG", "MRQWQWQWR"),
    entry("Q70001", "K1H1_HUMAN", "KRT31", "MCCCCDDDKFFFFGGGGR"),
    entry("Q70009", "K1H1_SHEEP", "KRT31", "MCCCCDDDKMMMMNNNNR", True),
    entry("P60712", "ACTB_BOVIN", "ACTB", ACTB, True),
    entry("Q9XSJ4", "ENOA_BOVIN", "ENO1", ENO1_BOVIN, True),
    entry("P02769", "ALBU_BOVIN", "ALB", BSA, True),
    entry("Q90002", "TGT_HORSE", "TGT", TGT_HORSE, True),
    entry("P00761", "TRYP_PIG", "", TRYP, True)])

COLS = ["Run", "Precursor.Id", "Protein.Ids", "Stripped.Sequence", "Precursor.Quantity",
        "Precursor.Normalised", "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.Q.Value"]

# (Protein.Ids, peptide, Precursor.Quantity per run)
CLEAN = [("P04406", "DDDDDDDK", 500), ("P04406", "CCCCCCCR", 400),
         ("Cont_P00761", "IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYK", 50)]
COMPETED = [("P60709;Cont_P60712", "GGGGGGGK", 900), ("P60709;Cont_P60712", "HHHHHHHR", 700)]
SIPHONED = [("P06733", "AAAAAAAK", 5), ("P06733;Cont_Q9XSJ4", "WWWWWWWK", 1000),
            ("P06733;Cont_Q9XSJ4", "YYYYYYYR", 800)]
MASKED = [("P02768;Cont_P02769", "VVVVVVVK", 300), ("Cont_P02769", "NNNNNNNK", 200),
          ("Cont_P02769", "QQQQQQQR", 200)]
PARALOGS = [("P60709;P68133;Cont_P60712", "GGGGGGGK", 900),
            ("P60709;P68133;Cont_P60712", "HHHHHHHR", 700)]
IL_TWIN = [("Cont_Q90002", "GHHEAELKPLAQSHATK", 600), ("Q90001;Cont_Q90002", "ASEDLKK", 100)]
# The review's case: XG's one taken precursor, one it shares with YG -- whose own three peptides
# name YG alone, in rows that name no seed accession.
CLOSURE = [("Q80001;Cont_P60712", "WQWQWQWQK", 900), ("Q80001;Q80002", "QWQWQWQK", 2000),
           ("Q80002", "QQWWQQWK", 50), ("Q80002", "WWQQWWQK", 50), ("Q80002", "QWWQQWWK", 50)]

# A hair sample's keratin, shared with the sheep wool keratin entry, in a database not built with
# --keratin-sample: the DE keeps it (contaminants.R is_contaminant(exempt=)).
KERATIN = [("Q70001;Cont_Q70009", "CCCCDDDK", 800)]
# ACTA1's one peptide of its own is shared with bovine actin; its others it shares with ACTB,
# which has one of its own.
SHARED_ELSEWHERE = [("P68133;Cont_P60712", "PPPPPPPK", 500), ("P60709;P68133", "GGGGGGGK", 900),
                    ("P60709;P68133", "HHHHHHHR", 700), ("P60709", "EEEEEEEK", 300)]

# The DE's contaminant record (de_provenance.json `contaminants`, contaminants.R)
DE_REMOVED = {"policy": "removed", "removed": True,
              "rule_origin": "the rule of DIA-NN's --cont-quant-exclude, applied here"}
DE_KEPT = {"policy": "kept", "removed": False}


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def write(path, text):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def report_rows(rows, runs=("r1", "r2"), failing=()):
    """Rows for each run; `failing` rows fail the q-value filter (Q.Value 0.5)."""
    out = []
    for run in runs:
        for ids, pep, q in rows:
            out.append([run, f"{pep}2", ids, pep, q, q * 3, 0.001, 0.001, 0.001, 0.01])
        for ids, pep, q in failing:
            out.append([run, f"{pep}2", ids, pep, q, q * 3, 0.5, 0.001, 0.001, 0.01])
    return out


def write_tsv(path, rows, failing=()):
    write(path, "\t".join(COLS) + "\n" + "".join(
        "\t".join(str(v) for v in r) + "\n" for r in report_rows(rows, failing=failing)))
    return path


class Fixture(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.fasta = os.path.join(self.d, "search.fasta")
        write(self.fasta, FASTA)
        self.meta = {"organism": "Homo sapiens", "fasta": self.fasta}
        write(self.fasta + ".meta.json", json.dumps(self.meta))

    def tearDown(self):
        self._tmp.cleanup()

    def run_check(self, rows, failing=(), fasta=True):
        rep = write_tsv(os.path.join(self.d, "report.tsv"), rows, failing)
        return cc.check(rep, self.fasta if fasta else None, self.meta if fasta else None)


class TheCheck(Fixture):
    def test_a_gene_whose_every_peptide_went_to_a_contaminant_is_competed_out(self):
        rec = self.run_check(CLEAN + COMPETED)
        self.assertEqual(rec["status"], "ok")
        (co,) = rec["competed_out"]
        self.assertEqual((co["gene"], co["contaminants"], co["n_removed_precursors"],
                          co["n_kept_precursors"]), ("ACTB", ["Cont_P60712"], 2, 0))
        self.assertEqual(rec["siphoned"], [])
        self.assertEqual(rec["drop_contaminant"], ["Cont_P60712"])
        self.assertIn("Cont_P60712", rec["contaminants_not_showing_themselves"])

    def test_a_gene_that_kept_the_minority_of_its_signal_is_siphoned(self):
        rec = self.run_check(CLEAN + SIPHONED)
        (si,) = rec["siphoned"]
        self.assertEqual((si["gene"], si["n_kept_precursors"], si["n_removed_precursors"]),
                         ("ENO1", 1, 2))
        self.assertGreater(si["removed_signal_pct"], 99)
        self.assertEqual(rec["competed_out"], [])
        self.assertEqual(rec["drop_contaminant"], ["Cont_Q9XSJ4"])

    def test_a_peptide_shared_by_paralogs_names_one_gene_not_each(self):
        """HeLa on HIVE: one actin peptide named ACTB, ACTG1, ACTA1, ACTC1, ACTA2, ACTG2 and
        POTEF, and each was listed as missing. Grouped as protein inference groups them."""
        rec = self.run_check(CLEAN + PARALOGS)
        (co,) = rec["competed_out"]
        self.assertEqual((co["gene"], co["indistinguishable"]), ("ACTB", ["ACTA1"]))
        self.assertIn("ACTB (= ACTA1; 2 precursors via Cont_P60712)", cc.staff_message(rec))
        self.assertTrue(cc.plain_competed(rec, DE_REMOVED).startswith("ACTB is missing"))

    def test_a_paralog_that_keeps_peptides_of_its_own_takes_the_shared_ones(self):
        rec = self.run_check(CLEAN + PARALOGS + [("P68133", "PPPPPPPK", 10)])
        self.assertEqual(rec["competed_out"], [])
        (si,) = rec["siphoned"]
        self.assertEqual((si["gene"], si["indistinguishable"], si["n_kept_precursors"]),
                         ("ACTA1", [], 1))

    def test_the_kept_evidence_is_read_to_closure(self):
        """Review of 2.14 (M1): only the kept rows naming a seed accession were read, so YG's own
        peptides were never seen, XG took the peptide it shares with YG, and came out SIPHONED --
        or, with that peptide's signal above the taken one, as a PASS while XG had nothing left."""
        rec = self.run_check(CLEAN + CLOSURE)
        self.assertEqual([c["gene"] for c in rec["competed_out"]], ["XG"])
        self.assertEqual(rec["siphoned"], [])
        self.assertEqual(rec["n_genes_lost_minority"], 0)
        kept = {ids: v for ids, v in cc.read_report(os.path.join(self.d, "report.tsv"))[0]
                ["kept"].items()}
        self.assertEqual(cc.connected(kept, {"Q80001"}, {}, {}), {"Q80001;Q80002", "Q80002"})

    def test_a_keratin_sample_s_exempt_keratin_is_not_taken(self):
        """Re-review of 2.14 (H1): the DE keeps a keratin sample's keratin-family precursors, and
        the client read "KRT31 is missing from the results"."""
        self.assertEqual([c["gene"] for c in self.run_check(CLEAN + KERATIN)["competed_out"]],
                         ["KRT31"])
        rep = os.path.join(self.d, "report.tsv")
        ex = {"shared": {"Cont_Q70009"}, "alone": set()}
        rec = cc.check(rep, self.fasta, self.meta, exempt=ex)
        self.assertEqual(rec["competed_out"], [])
        self.assertEqual(rec["keratin_exempt"]["shared"], ["Cont_Q70009"])
        self.assertIn("A keratin sample", cc.staff_message(rec))
        # the same exemption read from the search's records when there is no DE yet
        write(os.path.join(self.d, "search_provenance.json"), json.dumps({"keratin_sample": True}))
        write(self.fasta + ".meta.json", json.dumps(dict(
            self.meta, keratin_contaminants_in_database=["Cont_Q70009"])))
        self.assertEqual(cc.check(rep, self.fasta, self.meta)["competed_out"], [])

    def test_the_share_leaves_out_what_the_exemption_keeps(self):
        """Third review of 2.14 (L4): run_de.R's share leaves out exempted keratin."""
        self.assertEqual(self.run_check(CLEAN + KERATIN)["contaminant_share_pct"]["median"], 48.57)
        rec = cc.check(os.path.join(self.d, "report.tsv"), self.fasta, self.meta,
                       exempt={"shared": {"Cont_Q70009"}, "alone": set()})
        self.assertEqual(rec["contaminant_share_pct"]["median"], 2.86)

    def test_the_de_s_keratin_exemption_is_read_and_compared(self):
        k = {"value": True, "exempt_accessions": ["Cont_Q70009"], "exempt_alone_accessions": []}
        de = {"contaminants": dict(DE_REMOVED, keratin_sample=k)}
        self.assertEqual(cc.de_exemption(de), {"shared": {"Cont_Q70009"}, "alone": set()})
        self.assertIsNone(cc.de_exemption({"contaminants": DE_REMOVED}))
        # a MaxLFQ DE: DIA-NN's own quantities left them out -- nothing is exempt
        uq = {"contaminants": dict(DE_REMOVED, keratin_sample=dict(k, under_quantified=True))}
        self.assertEqual(cc.de_exemption(uq), {"shared": set(), "alone": set()})
        rec = self.run_check(CLEAN + KERATIN)
        self.assertFalse(cc.agrees_with_de(rec, de))             # made without the exemption
        ex = cc.check(rec["report"], self.fasta, self.meta, exempt=cc.de_exemption(de))
        self.assertTrue(cc.agrees_with_de(ex, de))
        alone = dict(k, exempt_alone_accessions=["Cont_Q70009"])
        self.assertFalse(cc.agrees_with_de(ex, {"contaminants": dict(DE_REMOVED,
                                                                     keratin_sample=alone)}))

    @unittest.skipUnless(shutil.which("Rscript"), "R is not installed")
    def test_the_exemption_is_contaminants_r_s_own(self):
        ids = ["P1;Cont_K2", "Cont_K1", "Cont_K2", "P1;Cont_K2;Cont_X", "P1;Cont_X", "P1"]
        expr = ('source("{}"); ex <- list(shared=c("Cont_K1","Cont_K2"), alone=c("Cont_K1")); '
                'cat(as.integer(is_contaminant(c({}), exempt=ex)))').format(
            os.path.join(SCRIPTS, "contaminants.R"), ", ".join(f'"{i}"' for i in ids))
        r = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        ex = {"shared": {"Cont_K1", "Cont_K2"}, "alone": {"Cont_K1"}}
        self.assertEqual([int(x) for x in r.stdout.split()],
                         [int(cc._is_cont(i) and not cc.exempted(i, ex)) for i in ids])

    def test_running_out_of_time_is_a_check_that_did_not_run(self):
        """Third review of 2.14 (M1): the stop was a TimeoutError, an OSError, which read_fasta's
        handler took for an unreadable FASTA -- an "ok" record with accessions for names."""
        import time
        rep = write_tsv(os.path.join(self.d, "report.tsv"), CLEAN + COMPETED)
        out = os.path.join(self.d, "rec.json")

        def slow(_fasta):
            time.sleep(5)
        with mock.patch.object(cc, "job_budget", return_value=(1, "a test's 1 s")), \
                mock.patch.object(cc, "MIN_BUDGET_S", 0), \
                mock.patch.object(cc, "read_fasta", side_effect=slow), \
                mock.patch("sys.stdout"), mock.patch("sys.stderr"):
            cc.main(["--report", rep, "--fasta", self.fasta, "--out", out, "--job-budget", "900"])
        rec = json.loads(read(out))
        self.assertEqual(rec["status"], "not_run")
        self.assertIn("stopped after 1 s", rec["reason"])

    def test_the_exemption_mirrors_contaminants_r(self):
        ex = {"shared": {"Cont_K1", "Cont_K2"}, "alone": {"Cont_K1"}}
        for ids, kept in (("P1;Cont_K2", True), ("Cont_K1", True), ("Cont_K2", False),
                          ("P1;Cont_K2;Cont_X", False), ("P1;Cont_X", False), ("P1", False)):
            with self.subTest(ids=ids):
                self.assertEqual(cc.exempted(ids, ex), kept)
        self.assertFalse(cc.exempted("P1;Cont_K2", None))

    def test_peptides_a_paralog_explains_are_not_called_the_gene_s(self):
        """Re-review of 2.14 (M1): ACTA1's other peptides are shared with ACTB, so "every peptide
        of it is shared with the contaminant list" was false."""
        rec = self.run_check(CLEAN + SHARED_ELSEWHERE)
        (co,) = rec["competed_out"]
        self.assertEqual((co["gene"], co["n_precursors_explained_by_other_genes"],
                          co["other_genes"]), ("ACTA1", 2, ["ACTB"]))
        self.assertIn("ACTA1 (1 precursor via Cont_P60712; 2 more precursors shared with ACTB, "
                      "counted there)", cc.staff_message(rec))
        plain = cc.plain_competed(rec, DE_REMOVED)
        self.assertIn("the peptides that tell it apart from other proteins are all shared", plain)
        self.assertNotIn("every peptide", plain)

    def test_a_gene_mostly_on_a_contaminant_that_shows_itself_is_said(self):
        """Re-review of 2.14 (L1): ALB with one weak peptide of its own and most of its signal on
        one shared with bovine albumin was reported nowhere."""
        rec = self.run_check(CLEAN + MASKED + [("P02768", "TTTTTTTK", 5)])
        self.assertEqual([(m["gene"], m["n_kept_precursors"]) for m in rec["masked"]],
                         [("ALB", 1)])
        (f,) = cc.audit_findings(rec)
        self.assertEqual((f[0], f[3]), ("PASS", cc.plain_masked(rec)))
        self.assertIn("cannot be told apart", f[3])

    def test_the_masked_words_follow_what_was_done_with_the_shared_peptides(self):
        """Third review of 2.14 (M2): with the shared peptides removed, a masked protein with no
        other peptide is not in the results, and the client read "no protein lost its signal"."""
        rec = self.run_check(CLEAN + MASKED)
        (f,) = cc.audit_findings(rec, DE_REMOVED)
        self.assertIn("one with none left is not in the results", f[3])
        self.assertNotIn("lost all or most", f[3])
        (f,) = cc.audit_findings(rec, DE_KEPT)
        self.assertNotIn("left out of the measurements", f[3])

    def test_only_a_true_il_twin_is_listed_as_one(self):
        rec = self.run_check(CLEAN + [("Cont_Q90002", "ASEDLKK", 100)])
        (co,) = rec["competed_out"]
        self.assertEqual(co["gene"], "TGT")                  # a TGT digest peptide, as written
        self.assertNotIn("il_twin_peptides", co)

    def test_a_clean_search_finds_nothing(self):
        rec = self.run_check(CLEAN)
        self.assertEqual((rec["status"], rec["competed_out"], rec["siphoned"], rec["masked"]),
                         ("ok", [], [], []))
        self.assertEqual(rec["drop_contaminant"], [])

    def test_a_contaminant_that_shows_itself_masks_rather_than_takes(self):
        """Bovine albumin in a cell-culture sample: its own peptides show it is there, so the one
        it shares with human ALB cannot be called ALB's -- listed for staff, never as lost."""
        rec = self.run_check(CLEAN + MASKED)
        self.assertEqual((rec["competed_out"], rec["siphoned"]), ([], []))
        self.assertEqual([m["gene"] for m in rec["masked"]], ["ALB"])
        self.assertIn("Cont_P02769", rec["contaminants_showing_themselves"])
        self.assertIn("Not counted as lost", cc.staff_message(rec))

    def test_one_own_peptide_is_not_enough_to_show_itself(self):
        rec = self.run_check(CLEAN + MASKED[:2])          # BSA with one peptide of its own
        self.assertEqual([c["gene"] for c in rec["competed_out"]], ["ALB"])

    def test_an_il_twin_is_the_sample_protein_s_peptide(self):
        """The horse entry's only lone peptide is the target's with I/L swapped: no evidence of
        the horse protein, so TGT lost both."""
        rec = self.run_check(CLEAN + IL_TWIN)
        (co,) = rec["competed_out"]
        self.assertEqual(co["gene"], "TGT")
        self.assertEqual(co["il_twin_peptides"], ["GHHEAELKPLAQSHATK"])
        self.assertEqual(co["n_removed_precursors"], 2)

    def test_the_q_value_filter_is_run_des(self):
        # an ACTB peptide of its own that fails Q.Value would otherwise keep ACTB
        rec = self.run_check(CLEAN + COMPETED, failing=[("P60709", "EEEEEEEK", 50)])
        self.assertEqual([c["gene"] for c in rec["competed_out"]], ["ACTB"])
        self.assertEqual(rec["q_columns"], ["Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
                                            "PG.Q.Value"])
        rec = self.run_check(CLEAN + COMPETED + [("P60709", "EEEEEEEK", 50)])
        self.assertEqual(rec["competed_out"], [])
        self.assertEqual([s["gene"] for s in rec["siphoned"]], ["ACTB"])

    def test_the_share_is_contaminants_r_s(self):
        rec = self.run_check(CLEAN + COMPETED)
        s = rec["contaminant_share_pct"]
        self.assertEqual(s["column"], "Precursor.Quantity")
        # per run: (900 + 700 + 50) / (500 + 400 + 50 + 900 + 700)
        self.assertEqual(s["per_run"], {"r1": 64.71, "r2": 64.71})
        self.assertEqual((s["median"], s["max"]), (64.71, 64.71))
        src = read(os.path.join(SCRIPTS, "contaminants.R"))
        m = re.search(r"CONTAMINANT_SHARE_COLUMNS <- c\(([^)]*)\)", src)
        self.assertEqual(tuple(re.findall(r'"([^"]+)"', m.group(1))), cc.SHARE_COLUMNS)

    def test_without_the_fasta_it_still_runs_and_says_what_it_could_not_do(self):
        rec = self.run_check(CLEAN + COMPETED, fasta=False)
        self.assertEqual([c["gene"] for c in rec["competed_out"]], ["P60709"])
        self.assertIn("no searched FASTA", rec["fasta_note"])

    def test_a_missing_report_is_not_run_with_the_reason(self):
        rec = cc.check(os.path.join(self.d, "nowhere", "report.parquet"), self.fasta, self.meta)
        self.assertEqual(rec["status"], "not_run")
        self.assertIn("not readable here", rec["reason"])
        (f,) = cc.audit_findings(rec)
        self.assertEqual(f[0], "INFO")
        self.assertTrue(f[1].startswith("Not run: the report"))

    def test_a_report_that_is_not_dianns_is_not_run(self):
        rep = os.path.join(self.d, "other.tsv")
        write(rep, "Run\tProtein\tIntensity\nr1\tP1\t5\n")
        rec = cc.check(rep, self.fasta, self.meta)
        self.assertEqual(rec["status"], "not_run")
        self.assertIn("not a DIA-NN precursor report", rec["reason"])

    def test_the_check_never_raises(self):
        """Review of 2.14 (M3): a parquet Q.Value stored as text made pyarrow raise
        ArrowNotImplementedError out of check(), and the audit died before AUDIT.md."""
        rep = write_tsv(os.path.join(self.d, "report.tsv"), CLEAN + COMPETED)
        for exc in (RuntimeError("boom"), KeyError("Run"), NotImplementedError("less_equal")):
            with self.subTest(exc=type(exc).__name__), \
                    mock.patch.object(cc, "read_report", side_effect=exc):
                rec = cc.check(rep, self.fasta, self.meta)
                self.assertEqual(rec["status"], "not_run")
                self.assertIn(f"the check failed on report.tsv ({type(exc).__name__}",
                              rec["reason"])
                self.assertEqual(cc.audit_findings(rec)[0][0], "INFO")

    def test_a_digestion_enzyme_is_never_offered_for_dropping(self):
        rec = self.run_check(CLEAN + [("P04406;Cont_P00761", "SSSSSSSSK", 9000)])
        self.assertEqual([s["gene"] for s in rec["siphoned"]], ["GAPDH"])
        self.assertEqual((rec["drop_contaminant"], rec["kept_digestion_enzymes"]),
                         ([], ["Cont_P00761"]))


class TheRecordDescribesOneReport(Fixture):
    def test_a_changed_report_is_not_described_by_the_old_record(self):
        rec = self.run_check(CLEAN + COMPETED)
        rep = rec["report"]
        self.assertTrue(cc.describes(rec, rep))
        write_tsv(rep, CLEAN)                         # a re-search into the same folder
        os.utime(rep, (1, 1))
        self.assertFalse(cc.describes(rec, rep))
        self.assertFalse(cc.describes(dict(rec, report_stat=None), rep))
        self.assertFalse(cc.describes(rec, os.path.join(self.d, "gone.tsv")))

    def test_a_record_made_with_other_cutoffs_does_not_stand_for_the_de(self):
        rec = self.run_check(CLEAN + COMPETED)
        self.assertTrue(cc.agrees_with_de(rec, {}))
        self.assertTrue(cc.agrees_with_de(rec, {"q_cutoff": 0.01}))
        self.assertFalse(cc.agrees_with_de(rec, {"q_cutoff": 0.05}))
        self.assertFalse(cc.agrees_with_de(rec, {"q_columns": ["Q.Value"]}))


class TheQFilterIsOne(unittest.TestCase):
    """Review of 2.14 (M6): three copies of run_de.R's q filter, and the parquet ones dropped a
    NaN that limpa (and the text reader) keep."""

    def test_a_row_is_dropped_only_when_a_value_exceeds_its_cutoff(self):
        cuts = {"Q.Value": 0.01, "PG.Q.Value": 0.05}
        for row, keep in (({"Q.Value": "0.01", "PG.Q.Value": "0.05"}, True),
                          ({"Q.Value": "nan", "PG.Q.Value": ""}, True),
                          ({"Q.Value": "NA"}, True),
                          ({"Q.Value": "0.02", "PG.Q.Value": "0.01"}, False),
                          ({"Q.Value": "x"}, False)):
            with self.subTest(row=row):
                self.assertEqual(dq.row_passes(row, cuts), keep)

    def test_every_reader_uses_it(self):
        for name in ("contaminant_competition.py", "audit_results.py", "pass_comparison.py"):
            src = read(os.path.join(SCRIPTS, name))
            with self.subTest(script=name):
                self.assertIn("keep_expression", src)
                self.assertNotIn("<= cut", src)


class TheWords(Fixture):
    TECHNICAL = ("Cont_", "precursor", "--drop", "q-value", "COMPETED", "SIPHONED", "P60712",
                 "fetch_fasta", "run_de", "DIA-NN")

    def test_competed_out_names_the_genes_says_why_and_offers_the_re_search(self):
        rec = self.run_check(CLEAN + COMPETED + IL_TWIN)
        (f,) = cc.audit_findings(rec, DE_REMOVED)
        status, msg, detail, plain = f
        self.assertEqual(status, "WARN")
        self.assertIn("COMPETED OUT", msg)
        self.assertIn("The DE removed every contaminant precursor before quantification (the rule "
                      "of DIA-NN's --cont-quant-exclude, applied here)", msg)
        self.assertIn("--drop-contaminant Cont_P60712 --drop-contaminant Cont_Q90002", msg)
        self.assertIn("search again", msg)
        self.assertEqual(detail["competed_out"], ["ACTB", "TGT"])
        self.assertTrue(plain.startswith("ACTB, TGT are missing from the results because the "
                                         "contaminant list absorbed their peptides"), plain)
        self.assertIn("The Core can search again without those contaminant entries", plain)
        for t in self.TECHNICAL:
            self.assertNotIn(t, plain)

    def test_the_words_follow_the_records_never_one_pipeline(self):
        """Review of 2.14 (M2, rule 1): "--cont-quant-exclude and run_de.R's contaminant filter
        removed them all" was said for a DE with --keep-contaminants, a cfg without the flag and
        any other DE pipeline alike."""
        rec = self.run_check(CLEAN + COMPETED + SIPHONED)
        cases = {
            "kept": (DE_KEPT, "--keep-contaminants", False),
            "flag recorded": (None, "DIA-NN's --cont-quant-exclude Cont_ (x.log.txt)", True),
            "flag absent": (None, "DIA-NN ran without --cont-quant-exclude", False),
            "not recorded": (None, "[not recorded -- confirm] (no DIA-NN log)", False)}
        cq = {"flag recorded": {"value": "Cont_", "source": "x.log.txt"},
              "flag absent": {"value": None, "source": "x.log.txt"},
              "not recorded": {"not_recorded": "no DIA-NN log"}}
        for name, (de, said, taken) in cases.items():
            with self.subTest(case=name):
                r = dict(rec, cont_quant_exclude=cq.get(name, rec.get("cont_quant_exclude")))
                (co, si) = cc.audit_findings(r, de)
                for f in (co, si):
                    self.assertIn(said, f[1])
                    self.assertNotIn("run_de.R's contaminant filter removed", f[1])
                    for t in self.TECHNICAL:
                        self.assertNotIn(t, f[3])
                self.assertEqual(co[3].startswith("ACTB is missing from the results"), taken)
                self.assertEqual("remaining peptides only" in si[3], taken)
                if not taken:
                    self.assertIn("cannot tell its signal from that entry's", co[3])
                    self.assertIn("cannot be told apart from the contaminant", si[3])

    def test_the_search_record_says_what_its_own_reader_found(self):
        rec = self.run_check(CLEAN + COMPETED)
        self.assertIn("not_recorded", rec["cont_quant_exclude"])      # no DIA-NN log here
        self.assertIn("[not recorded -- confirm]", cc.staff_message(rec))

    def test_a_database_built_by_an_older_rule_is_told_to_be_rebuilt(self):
        rep = write_tsv(os.path.join(self.d, "report.tsv"), CLEAN + COMPETED)
        old = dict(self.meta, contaminant_set="universal", n_contaminants_appended=381)
        msg = cc.audit_findings(cc.check(rep, self.fasta, old))[0][1]
        self.assertIn("older fetch_fasta.py rule (legacy)", msg)
        self.assertIn("--drop-contaminant Cont_P60712", msg)
        self.assertNotIn("older", cc.audit_findings(self.run_check(CLEAN + COMPETED))[0][1])

    def test_one_gene_reads_in_the_singular(self):
        plain = cc.audit_findings(self.run_check(CLEAN + COMPETED), DE_REMOVED)[0][3]
        self.assertTrue(plain.startswith("ACTB is missing from the results because the "
                                         "contaminant list absorbed its peptides"), plain)

    def test_siphoned_is_a_caveat_naming_the_genes(self):
        (f,) = cc.audit_findings(self.run_check(CLEAN + SIPHONED), DE_REMOVED)
        self.assertEqual(f[0], "WARN")
        self.assertIn("SIPHONED: 1 Homo sapiens gene kept precursors of its own", f[1])
        self.assertIn("ENO1", f[3])
        self.assertIn("remaining peptides only", f[3])
        for t in self.TECHNICAL:
            self.assertNotIn(t, f[3])

    def test_both_are_two_rows_of_findings_and_never_a_fail(self):
        fs = cc.audit_findings(self.run_check(CLEAN + COMPETED + SIPHONED), DE_REMOVED)
        self.assertEqual([f[0] for f in fs], ["WARN", "WARN"])
        self.assertEqual(sum("--drop-contaminant" in f[1] for f in fs), 1)   # the fix once
        md = ar.client_section_md([{"check": "contaminant_competition", "status": s,
                                    "message": m, "plain": p} for s, m, _d, p in fs])
        row = [ln for ln in md.splitlines() if ln.startswith("| Sample proteins")][0]
        self.assertIn("⚠️ Caution", row)
        self.assertIn("ACTB is missing", row)
        self.assertIn("ENO1", row)

    def test_clean_is_a_short_checked_none(self):
        (f,) = cc.audit_findings(self.run_check(CLEAN))
        self.assertEqual((f[0], f[3]), ("PASS", None))
        self.assertTrue(f[1].startswith("Checked: no Homo sapiens gene lost all or most"), f[1])
        label, meaning = ar.client_row({"check": "contaminant_competition", "status": "PASS",
                                        "message": f[1]})
        self.assertEqual(label, "Sample proteins taken by the contaminant list")
        self.assertEqual(meaning, "Checked: no protein of the sample lost all or most of its "
                                  "signal to the contaminant list.")

    def test_an_older_audit_md_reads_the_check_s_words(self):
        label, meaning = ar.client_row({"check": "contaminant_competition", "status": "WARN",
                                        "message": "COMPETED OUT: ..."})
        self.assertNotEqual(meaning, ar.GENERIC_MEANING["WARN"])
        self.assertIn("contaminant list", meaning)

    def test_not_run_reads_as_not_checked(self):
        label, meaning = ar.client_row({"check": "contaminant_competition", "status": "INFO",
                                        "message": "Not run: x."})
        self.assertEqual(meaning, ar.NOT_CHECKED)

    def test_another_engine_s_report_says_nothing(self):
        self.assertEqual(cc.audit_findings({"status": "not_applicable", "reason": "x"}), [])


class TheRecordAndTheAudit(Fixture):
    def search_out(self, rows, prov_extra=None, failing=()):
        out = os.path.join(self.d, "search_out")
        rep = write_tsv(os.path.join(out, "report.tsv"), rows, failing)
        write(os.path.join(out, "search_provenance.json"), json.dumps(
            dict({"engine": "diann", "fasta": self.fasta, "result": {"report": rep}},
                 **(prov_extra or {}))))
        return out, rep

    def cli(self, out):
        return subprocess.run([sys.executable, CHECK, "--search-out", out], capture_output=True,
                              text=True, timeout=120)

    def test_the_cli_writes_the_record_and_the_provenance_and_exits_0(self):
        out, rep = self.search_out(CLEAN + COMPETED)
        r = self.cli(out)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("[contaminant_competition] WARNING: COMPETED OUT", r.stderr)
        rec = json.loads(read(os.path.join(out, cc.RECORD)))
        prov = json.loads(read(os.path.join(out, "search_provenance.json")))
        self.assertEqual([c["gene"] for c in rec["competed_out"]], ["ACTB"])
        self.assertEqual(prov["contaminant_competition"]["competed_out"], rec["competed_out"])
        self.assertEqual(prov["engine"], "diann")                 # merged, never replaced
        self.assertEqual(json.loads(r.stdout)["competed_out"], ["ACTB"])
        self.assertTrue(cc.describes(rec, rep))

    def test_the_cli_on_a_missing_report_records_not_run_and_exits_0(self):
        r = subprocess.run([sys.executable, CHECK, "--report", os.path.join(self.d, "none.parquet"),
                            "--out", os.path.join(self.d, "rec.json")], capture_output=True,
                           text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("NOTE: Not run", r.stderr)
        self.assertEqual(json.loads(read(os.path.join(self.d, "rec.json")))["status"], "not_run")

    def test_an_unreadable_search_record_is_named_as_the_reason(self):
        out, rep = self.search_out(CLEAN)
        write(os.path.join(out, "search_provenance.json"), "{not json")
        _rep, fasta, _meta, why = cc.resolve_inputs(out)
        self.assertIsNone(fasta)
        self.assertIn("search_provenance.json in", why)

    def audit(self, de_prov=None, search_out=None):
        args = [sys.executable, AUDIT, "--out", os.path.join(self.d, "AUDIT.md"),
                "--fasta-meta", self.fasta + ".meta.json"]
        if de_prov is not None:
            de = os.path.join(self.d, "de")
            write(os.path.join(de, "de_provenance.json"), json.dumps(de_prov))
            args += ["--de-dir", de]
        if search_out:
            args += ["--search-out", search_out]
        r = subprocess.run(args, capture_output=True, text=True, timeout=180, cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        audit = json.loads(read(os.path.join(self.d, "AUDIT.json")))
        md = read(os.path.join(self.d, "AUDIT.md"))
        return [f for f in audit["findings"] if f["check"] == "contaminant_competition"], md

    def test_the_audit_quotes_the_search_s_record_in_the_de_s_words(self):
        out, rep = self.search_out(CLEAN + COMPETED)
        self.cli(out)
        fs, md = self.audit({"input": rep, "contaminants": DE_REMOVED}, out)
        (f,) = fs
        self.assertEqual(f["status"], "WARN")
        self.assertIn("ACTB is missing from the results", f["plain"])
        self.assertIn("**contaminant_competition** — COMPETED OUT", md)
        self.assertIn("The DE removed every contaminant precursor", f["message"])
        self.assertNotIn("checked_by", f["detail"])
        self.assertIn("ACTB is missing", ar.client_section_md([f]))

    def test_a_search_without_a_record_is_checked_by_the_audit(self):
        out, rep = self.search_out(CLEAN + SIPHONED)
        (f,), _md = self.audit({"input": rep, "q_cutoff": 0.01, "contaminants": DE_REMOVED})
        self.assertEqual(f["status"], "WARN")
        self.assertIn("ENO1", f["plain"])
        self.assertIn("checked_by", f["detail"])

    def test_a_record_of_another_report_is_not_used_for_this_one(self):
        out, rep = self.search_out(CLEAN)
        write(os.path.join(out, cc.RECORD), json.dumps(
            {"status": "ok", "report": os.path.join(out, "step3_assembly.parquet"),
             "competed_out": [{"gene": "WRONG"}]}))
        (f,), _md = self.audit({"input": rep}, out)
        self.assertEqual(f["status"], "PASS")

    def test_a_previous_search_s_record_is_not_read_for_a_re_search(self):
        """Review of 2.14 (M4): the step-5 job deletes the record but not its copy in
        search_provenance.json, and the audit matched the report by name alone."""
        out, rep = self.search_out(CLEAN + COMPETED)
        self.cli(out)                                 # ACTB taken in the first search
        write_tsv(rep, CLEAN)                         # the re-search's report, same name
        os.utime(rep, (1, 1))
        os.remove(os.path.join(out, cc.RECORD))       # step 5 cleared the file, not the copy
        (f,), _md = self.audit({"input": rep}, out)
        self.assertEqual(f["status"], "PASS")

    def test_a_record_made_with_other_cutoffs_is_checked_again(self):
        out, rep = self.search_out(CLEAN + COMPETED, failing=[("P60709", "EEEEEEEK", 50)])
        self.cli(out)                                 # at 1%: ACTB's own peptide fails
        (f,), _md = self.audit({"input": rep, "q_cutoff": 0.6, "contaminants": DE_REMOVED}, out)
        self.assertIn("SIPHONED", f["message"])       # at the DE's 60% it is kept
        self.assertIn("checked_by", f["detail"])

    def test_a_check_that_could_not_run_is_in_the_audit_with_its_reason(self):
        gone = os.path.join(self.d, "gone", "report.parquet")
        (f,), md = self.audit({"input": gone})
        self.assertEqual(f["status"], "INFO")
        self.assertIn("Not run: the report", f["message"])
        self.assertIn("ℹ️ **contaminant_competition** — Not run", md)

    def test_a_search_only_audit_still_says_a_check_did_not_run(self):
        """Review of 2.14 (M5): with no DE record a not_run record was dropped, so the audit
        said nothing at all."""
        out, rep = self.search_out(CLEAN)
        os.remove(rep)
        write(os.path.join(out, cc.RECORD), json.dumps(
            {"status": "not_run", "report": rep, "reason": "pyarrow is not installed"}))
        (f,), _md = self.audit(None, out)
        self.assertEqual(f["status"], "INFO")
        self.assertIn("Not run:", f["message"])

    def test_a_search_only_audit_reads_the_search_s_report(self):
        out, rep = self.search_out(CLEAN + COMPETED)
        (f,), _md = self.audit(None, out)
        self.assertEqual(f["status"], "WARN")
        self.assertIn("ACTB", f["plain"])

    def test_the_audit_is_written_whatever_the_check_does(self):
        out, rep = self.search_out(CLEAN + COMPETED)
        de = os.path.join(self.d, "de")
        write(os.path.join(de, "de_provenance.json"), json.dumps({"input": rep}))
        argv = ["audit_results.py", "--out", os.path.join(self.d, "AUDIT.md"), "--de-dir", de]
        with mock.patch.object(sys, "argv", argv), \
                mock.patch.object(cc, "audit_findings", side_effect=RuntimeError("boom")), \
                mock.patch("sys.stdout"):
            ar.main()
        audit = json.loads(read(os.path.join(self.d, "AUDIT.json")))
        (f,) = [x for x in audit["findings"] if x["check"] == "contaminant_competition"]
        self.assertEqual(f["status"], "INFO")
        self.assertIn("boom", f["message"])

    def test_another_engine_s_search_is_not_asked(self):
        out, rep = self.search_out(CLEAN + COMPETED, {"engine": "sage"})
        fs, _md = self.audit({"input": rep}, None)
        self.assertEqual(fs, [])

    def test_a_record_that_could_not_run_is_tried_again_by_the_audit(self):
        out, rep = self.search_out(CLEAN + COMPETED)
        rec = cc.check(rep, self.fasta, self.meta)
        rec.update(status="not_run", reason="pyarrow is not installed")
        cc.write(rec, os.path.join(out, cc.RECORD))
        (f,), _md = self.audit({"input": rep}, out)
        self.assertEqual(f["status"], "WARN")
        self.assertIn("ACTB", f["plain"])

    def test_only_a_record_with_a_status_goes_into_the_search_provenance(self):
        plan = {"contaminant_competition": {"file": "x", "produced": "runtime", "by": "s5"}}
        self.assertEqual(run_search.search_records(plan), {})
        rec = {"contaminant_competition": {"status": "ok"}, "cpu_sizing": {"x": 1}}
        self.assertEqual(run_search.search_records(rec), rec)


class TheJobsRunIt(unittest.TestCase):
    """Every DIA-NN search's last job ends with the check; it never fails the job."""

    def test_squeue_s_time_left_reads(self):
        for text, secs in (("1-02:03:04", 93784), ("02:03:04", 7384), ("05:00", 300), ("45", 45),
                           ("UNLIMITED", None), ("INVALID", None), ("", None)):
            with self.subTest(text=text):
                self.assertEqual(cc._slurm_seconds(text), secs)

    def test_the_check_stops_before_the_search_job_s_time_runs_out(self):
        """Re-review of 2.14 (M2): 15 minutes appended to a job that had used nearly all its time
        ended a complete search in TIMEOUT, whose fix is a resubmit."""
        margin = cc.job_margin()
        import notify_slack
        self.assertEqual(margin, notify_slack.JOB_END_HOOK_S + cc.JOB_MARGIN_EXTRA_S)
        done = subprocess.CompletedProcess([], 0, stdout="20:00\n", stderr="")
        with mock.patch.dict(os.environ, {"SLURM_JOB_ID": "42"}), \
                mock.patch.object(cc.subprocess, "run", return_value=done):
            self.assertEqual(cc.job_budget(900), (1200 - margin, f"1200 s were left in job 42, "
                                                                 f"{margin} s of them for its end"))
        done = subprocess.CompletedProcess([], 0, stdout="2-00:00:00\n", stderr="")
        with mock.patch.dict(os.environ, {"SLURM_JOB_ID": "42"}), \
                mock.patch.object(cc.subprocess, "run", return_value=done):
            self.assertEqual(cc.job_budget(900), (900, "the check's own limit of 900 s"))
        with mock.patch.dict(os.environ, {"SLURM_JOB_ID": "42"}), \
                mock.patch.object(cc.subprocess, "run", side_effect=OSError("no squeue")):
            self.assertEqual(cc.job_budget(900)[0], 900)
        with mock.patch.dict(os.environ, {}, clear=True):
            self.assertEqual(cc.job_budget(900), (900, "the check's own limit of 900 s"))
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "rec.json")
            left = (10, "130 s were left in job 7")
            with mock.patch.object(cc, "job_budget", return_value=left), \
                    mock.patch("sys.stdout"), mock.patch("sys.stderr"):
                self.assertEqual(cc.main(["--report", os.path.join(d, "report.parquet"),
                                          "--out", out, "--job-budget", "900"]), 0)
            rec = json.loads(read(out))
        self.assertEqual(rec["status"], "not_run")
        self.assertIn("too little of the search job's time was left", rec["reason"])

    def test_setting_a_record_aside_never_stops_a_generated_search(self):
        import diann_parallel
        with tempfile.TemporaryDirectory() as d, \
                mock.patch.object(diann_parallel, "set_aside", side_effect=OSError("read-only")), \
                mock.patch("sys.stderr"):
            self.assertIsNone(cc.set_aside_record(d))

    def test_the_chain_s_step5_runs_it_after_the_pass_comparison(self):
        import test_chain_pass_comparison as tcp
        with tempfile.TemporaryDirectory() as d:
            p, out = tcp.generate(d, tcp.PINNED)
            self.assertEqual(p.returncode, 0, p.stderr)
            s5 = read(os.path.join(out, "step5_report.sbatch"))
            info = json.loads(p.stdout)
            line = next(ln for ln in s5.splitlines() if "contaminant_competition.py" in ln)
            rec = os.path.join(out, cc.RECORD)
            self.assertIn(f"--report {os.path.join(out, 'report.parquet')}", line)
            self.assertIn(f"--fasta {os.path.join(d, 'db.fasta')}", line)
            self.assertIn(f"--out {rec}", line)
            self.assertIn(f"--provenance {os.path.join(out, 'search_provenance.json')}", line)
            self.assertIn(f'T="timeout {cc.JOB_TIME_LIMIT + cc.JOB_MARGIN_EXTRA_S}"', line)
            self.assertIn(f"--job-budget {cc.JOB_TIME_LIMIT}", line)       # bounded in the job
            self.assertIn('|| echo "WARNING', line)
            self.assertLess(s5.index("pass_comparison.py"), s5.index("contaminant_competition.py"))
            self.assertIn(f'"{rec}"', next(ln for ln in s5.splitlines() if ln.startswith("rm -f")))
            self.assertEqual(info["contaminant_competition"]["produced"], "runtime")
            # the line as generated, in bash: no report -> a not-run record, exit 0
            write(os.path.join(out, "search_provenance.json"), json.dumps({"engine": "diann"}))
            from job_env import job_env
            r = subprocess.run(["bash", "-c", line], capture_output=True, text=True,
                               env=job_env(d), timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(json.loads(read(rec))["status"], "not_run")
            self.assertEqual(json.loads(read(os.path.join(out, "search_provenance.json")))
                             ["contaminant_competition"]["status"], "not_run")
            # generating the chain again sets the record aside: it described another search
            p, out2 = tcp.generate(d, tcp.PINNED)
            self.assertEqual((p.returncode, out2), (0, out), p.stderr)
            self.assertFalse(os.path.exists(rec))
            self.assertTrue([f for f in os.listdir(out) if f.startswith(cc.RECORD + ".stale-")])
            self.assertIn("does not describe this search", p.stderr)

    def test_the_single_shot_search_job_runs_it_after_its_guard(self):
        import test_single_shot_sbatch as tss
        with tempfile.TemporaryDirectory() as d:
            w = tss._Workspace(d)
            w.generate()
            job = read(w.search_job)
            self.assertLess(job.index("check_report_runs.py"),
                            job.index("contaminant_competition.py"))
            self.assertIn(f"--fasta {w.fasta}", job)
            self.assertIn(os.path.join(w.out, cc.RECORD), next(
                ln for ln in job.splitlines() if ln.startswith("rm -f")))
            p = w.run_job(w.search_job, FAKE_REPORT_RUNS="sample_0,sample_1")
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("[contaminant_competition]", p.stderr)
            self.assertTrue(os.path.isfile(os.path.join(w.out, cc.RECORD)))
            w.generate()                              # a re-search into the same folder
            self.assertFalse(os.path.isfile(os.path.join(w.out, cc.RECORD)))


class TheWatcherSaysIt(Fixture):
    WATCH = os.path.join(SCRIPTS, "watch_run.sh")

    def watch(self, out, state="COMPLETED"):
        bindir = os.path.join(self.d, "bin")
        for name, body in (("sacct", f'echo "  {state}"'), ("squeue", "true")):
            write(os.path.join(bindir, name), f"#!/bin/bash\n{body}\n")
            os.chmod(os.path.join(bindir, name), 0o755)
        env = dict(os.environ, PATH=bindir + os.pathsep + os.environ.get("PATH", ""))
        r = subprocess.run(["bash", self.WATCH, "--slurm", "999", "--out", out],
                           capture_output=True, text=True, env=env, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        return json.loads(r.stdout)

    def recorded(self):
        out = os.path.join(self.d, "out")
        os.makedirs(out)
        rec = self.run_check(CLEAN + COMPETED)
        cc.write(rec, os.path.join(out, cc.RECORD))
        return out, rec

    def test_a_finished_search_names_the_genes(self):
        out, _rec = self.recorded()
        got = self.watch(out)
        self.assertEqual(got["contaminant_competition"]["competed_out"], ["ACTB"])
        self.assertEqual(got["contaminant_competition"]["drop_contaminant"], ["Cont_P60712"])
        self.assertTrue(any("COMPETED OUT" in w for w in got["warnings"]))
        self.assertFalse(got["failed"])

    def test_a_running_search_shows_no_record(self):
        out, _rec = self.recorded()
        self.assertNotIn("contaminant_competition", self.watch(out, "RUNNING"))

    def test_a_record_of_an_earlier_report_is_stale(self):
        out, rec = self.recorded()
        write_tsv(rec["report"], CLEAN)               # the report as a re-search left it
        os.utime(rec["report"], (1, 1))
        got = self.watch(out)
        self.assertEqual(got["contaminant_competition"]["status"], "stale")
        self.assertFalse(any("COMPETED OUT" in w for w in got.get("warnings", [])))

    def test_a_report_that_cannot_be_read_here_is_unchecked_not_stale(self):
        out, rec = self.recorded()
        os.rename(rec["report"], rec["report"] + ".moved")
        got = self.watch(out)
        self.assertEqual(got["contaminant_competition"]["status"], "unchecked")

    def test_no_record_no_key(self):
        os.makedirs(os.path.join(self.d, "out"))
        self.assertNotIn("contaminant_competition", self.watch(os.path.join(self.d, "out")))


@unittest.skipUnless(pq, "pyarrow is not installed")
class TheParquetReader(Fixture):
    def parquet(self, rows):
        par = os.path.join(self.d, "report.parquet")
        pq.write_table(pa.table({c: list(v) for c, v in zip(COLS, zip(*rows))}), par)
        return par

    def test_parquet_and_tsv_give_the_same_record(self):
        rows = CLEAN + COMPETED + SIPHONED + MASKED + IL_TWIN + CLOSURE
        failing = [("P60709", "EEEEEEEK", 50)]
        tsv = write_tsv(os.path.join(self.d, "report.tsv"), rows, failing)
        par = self.parquet(report_rows(rows, failing=failing))
        a, b = cc.check(tsv, self.fasta, self.meta), cc.check(par, self.fasta, self.meta)
        for k in ("competed_out", "siphoned", "masked", "drop_contaminant",
                  "contaminant_share_pct", "contaminants_showing_themselves",
                  "n_kept_protein_sets_read"):
            self.assertEqual(a[k], b[k], k)
        self.assertEqual([c["gene"] for c in b["competed_out"]], ["ACTB", "XG", "TGT"])

    def test_a_nan_q_value_is_kept_as_limpa_keeps_it(self):
        rows = report_rows(CLEAN + COMPETED + [("P60709", "EEEEEEEK", 50)])
        qi = COLS.index("Q.Value")
        for r in rows:
            if r[3] == "EEEEEEEK":
                r[qi] = float("nan")              # not above the cutoff: limpa keeps the row
        rec = cc.check(self.parquet(rows), self.fasta, self.meta)
        self.assertEqual(rec["competed_out"], [])
        self.assertEqual([s["gene"] for s in rec["siphoned"]], ["ACTB"])
        tsv = os.path.join(self.d, "report.tsv")
        write(tsv, "\t".join(COLS) + "\n" + "".join("\t".join(str(v) for v in r) + "\n"
                                                     for r in rows))
        self.assertEqual(cc.check(tsv, self.fasta, self.meta)["siphoned"], rec["siphoned"])

    def test_a_nan_quantity_counts_as_nothing_in_both_readers(self):
        """Re-review of 2.14 (L2): arrow's sum carried one NaN quantity into a run's whole total,
        and that run's share vanished."""
        rows = report_rows(CLEAN + COMPETED + SIPHONED)
        rows[1][COLS.index("Precursor.Quantity")] = float("nan")
        tsv = os.path.join(self.d, "report.tsv")
        write(tsv, "\t".join(COLS) + "\n" + "".join("\t".join(str(v) for v in r) + "\n"
                                                     for r in rows))
        a, b = cc.check(tsv, self.fasta, self.meta), cc.check(self.parquet(rows), self.fasta,
                                                               self.meta)
        self.assertEqual(sorted(a["contaminant_share_pct"]["per_run"]), ["r1", "r2"])
        self.assertEqual(a["contaminant_share_pct"], b["contaminant_share_pct"])
        self.assertEqual(a["siphoned"], b["siphoned"])

    def test_a_q_value_stored_as_text_is_a_check_that_did_not_run(self):
        rows = report_rows(CLEAN + COMPETED)
        t = pa.table({c: [str(x) if c == "Q.Value" else x for x in v]
                      for c, v in zip(COLS, zip(*rows))})
        par = os.path.join(self.d, "report.parquet")
        pq.write_table(t, par)
        rec = cc.check(par, self.fasta, self.meta)
        self.assertEqual(rec["status"], "not_run")
        self.assertIn("the check failed on report.parquet", rec["reason"])

    def test_an_adapted_report_is_not_dianns_question(self):
        par = os.path.join(self.d, "adapted.parquet")
        t = pa.table({c: [v] for c, v in zip(COLS, report_rows(CLEAN)[0])})
        pq.write_table(t.replace_schema_metadata({b"delimp.quantity.level": b"peptide"}), par)
        self.assertEqual(cc.check(par, self.fasta, self.meta)["status"], "not_applicable")


if __name__ == "__main__":
    unittest.main(verbosity=2)
