#!/usr/bin/env python3
"""The contaminant competition check names the sample proteins the contaminant list took away.

The Core searches every sample with its universal contaminant list. Entries that share peptides
with abundant human and mouse proteins took them away: DIA-NN reported the peptides under the
Cont_ entry, and --cont-quant-exclude Cont_ (and run_de.R's filter) removed them. A mouse tissue
search had no Actb or Eef1a1 at all; a human IP lost 68 genes. The Core director chose (2026-10-20)
to keep the list and have every DIA-NN search report what it takes: contaminant_competition.py,
run by the search's last job, recorded in search_provenance.json and AUDIT.md, and read to the
client in the report's Audit & caveats table. It never blocks.

Stdlib only; the parquet reader's test is skipped without pyarrow.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import audit_results as ar  # noqa: E402
import contaminant_competition as cc  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
except ImportError:                       # the stdlib-only CI: that test is skipped
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
        self.assertTrue(cc.plain_competed(rec).startswith("ACTB is missing"))

    def test_a_paralog_that_keeps_peptides_of_its_own_takes_the_shared_ones(self):
        rec = self.run_check(CLEAN + PARALOGS + [("P68133", "PPPPPPPK", 10)])
        self.assertEqual(rec["competed_out"], [])
        (si,) = rec["siphoned"]
        self.assertEqual((si["gene"], si["indistinguishable"], si["n_kept_precursors"]),
                         ("ACTA1", [], 1))

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

    def test_a_digestion_enzyme_is_never_offered_for_dropping(self):
        rec = self.run_check(CLEAN + [("P04406;Cont_P00761", "SSSSSSSSK", 9000)])
        self.assertEqual([s["gene"] for s in rec["siphoned"]], ["GAPDH"])
        self.assertEqual((rec["drop_contaminant"], rec["kept_digestion_enzymes"]),
                         ([], ["Cont_P00761"]))


class TheWords(Fixture):
    TECHNICAL = ("Cont_", "precursor", "--drop", "q-value", "COMPETED", "SIPHONED", "P60712",
                 "fetch_fasta", "run_de")

    def test_competed_out_names_the_genes_says_why_and_offers_the_re_search(self):
        rec = self.run_check(CLEAN + COMPETED + IL_TWIN)
        (f,) = cc.audit_findings(rec)
        status, msg, detail, plain = f
        self.assertEqual(status, "WARN")
        self.assertIn("COMPETED OUT", msg)
        self.assertIn("--drop-contaminant Cont_P60712 --drop-contaminant Cont_Q90002", msg)
        self.assertIn("search again", msg)
        self.assertEqual(detail["competed_out"], ["ACTB", "TGT"])
        self.assertTrue(plain.startswith("ACTB, TGT are missing from the results because the "
                                         "contaminant list absorbed their peptides"), plain)
        self.assertIn("The Core can search again without those contaminant entries", plain)
        for t in self.TECHNICAL:
            self.assertNotIn(t, plain)

    def test_a_database_built_by_an_older_rule_is_told_to_be_rebuilt(self):
        rep = write_tsv(os.path.join(self.d, "report.tsv"), CLEAN + COMPETED)
        old = dict(self.meta, contaminant_set="universal", n_contaminants_appended=381)
        msg = cc.audit_findings(cc.check(rep, self.fasta, old))[0][1]
        self.assertIn("older fetch_fasta.py rule (legacy)", msg)
        self.assertIn("--drop-contaminant Cont_P60712", msg)
        self.assertNotIn("older", cc.audit_findings(self.run_check(CLEAN + COMPETED))[0][1])

    def test_one_gene_reads_in_the_singular(self):
        plain = cc.audit_findings(self.run_check(CLEAN + COMPETED))[0][3]
        self.assertTrue(plain.startswith("ACTB is missing from the results because the "
                                         "contaminant list absorbed its peptides"), plain)

    def test_siphoned_is_a_caveat_naming_the_genes(self):
        (f,) = cc.audit_findings(self.run_check(CLEAN + SIPHONED))
        self.assertEqual(f[0], "WARN")
        self.assertIn("SIPHONED: 1 Homo sapiens gene kept", f[1])
        self.assertIn("ENO1", f[3])
        self.assertIn("remaining peptides only", f[3])
        for t in self.TECHNICAL:
            self.assertNotIn(t, f[3])

    def test_both_are_two_rows_of_findings_and_never_a_fail(self):
        fs = cc.audit_findings(self.run_check(CLEAN + COMPETED + SIPHONED))
        self.assertEqual([f[0] for f in fs], ["WARN", "WARN"])
        md = ar.client_section_md([{"check": "contaminant_competition", "status": s,
                                    "message": m, "plain": p} for s, m, _d, p in fs])
        row = [ln for ln in md.splitlines() if ln.startswith("| Sample proteins")][0]
        self.assertIn("⚠️ Caution", row)
        self.assertIn("ACTB is missing", row)
        self.assertIn("ENO1", row)

    def test_clean_is_a_short_checked_none(self):
        (f,) = cc.audit_findings(self.run_check(CLEAN))
        self.assertEqual((f[0], f[3]), ("PASS", None))
        self.assertTrue(f[1].startswith("Checked: no Homo sapiens gene lost"), f[1])
        label, meaning = ar.client_row({"check": "contaminant_competition", "status": "PASS",
                                        "message": f[1]})
        self.assertEqual(label, "Sample proteins taken by the contaminant list")
        self.assertTrue(meaning.startswith("Checked: no protein of the sample"), meaning)

    def test_not_run_reads_as_not_checked(self):
        label, meaning = ar.client_row({"check": "contaminant_competition", "status": "INFO",
                                        "message": "Not run: x."})
        self.assertEqual(meaning, ar.NOT_CHECKED)

    def test_another_engine_s_report_says_nothing(self):
        self.assertEqual(cc.audit_findings({"status": "not_applicable", "reason": "x"}), [])


class TheRecordAndTheAudit(Fixture):
    def search_out(self, rows, prov_extra=None):
        out = os.path.join(self.d, "search_out")
        rep = write_tsv(os.path.join(out, "report.tsv"), rows)
        write(os.path.join(out, "search_provenance.json"), json.dumps(
            dict({"engine": "diann", "fasta": self.fasta, "result": {"report": rep}},
                 **(prov_extra or {}))))
        return out, rep

    def test_the_cli_writes_the_record_and_the_provenance_and_exits_0(self):
        out, rep = self.search_out(CLEAN + COMPETED)
        r = subprocess.run([sys.executable, CHECK, "--search-out", out], capture_output=True,
                           text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("[contaminant_competition] WARNING: COMPETED OUT", r.stderr)
        rec = json.loads(read(os.path.join(out, cc.RECORD)))
        prov = json.loads(read(os.path.join(out, "search_provenance.json")))
        self.assertEqual([c["gene"] for c in rec["competed_out"]], ["ACTB"])
        self.assertEqual(prov["contaminant_competition"]["competed_out"], rec["competed_out"])
        self.assertEqual(prov["engine"], "diann")                 # merged, never replaced
        self.assertEqual(json.loads(r.stdout)["competed_out"], ["ACTB"])

    def test_the_cli_on_a_missing_report_records_not_run_and_exits_0(self):
        r = subprocess.run([sys.executable, CHECK, "--report", os.path.join(self.d, "none.parquet"),
                            "--out", os.path.join(self.d, "rec.json")], capture_output=True,
                           text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("NOTE: Not run", r.stderr)
        self.assertEqual(json.loads(read(os.path.join(self.d, "rec.json")))["status"], "not_run")

    def audit(self, de_prov, search_out=None):
        de = os.path.join(self.d, "de")
        write(os.path.join(de, "de_provenance.json"), json.dumps(de_prov))
        args = [sys.executable, AUDIT, "--out", os.path.join(self.d, "AUDIT.md"), "--de-dir", de,
                "--fasta-meta", self.fasta + ".meta.json"]
        if search_out:
            args += ["--search-out", search_out]
        r = subprocess.run(args, capture_output=True, text=True, timeout=180, cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stderr)
        audit = json.loads(read(os.path.join(self.d, "AUDIT.json")))
        md = read(os.path.join(self.d, "AUDIT.md"))
        return [f for f in audit["findings"] if f["check"] == "contaminant_competition"], md

    def test_the_audit_quotes_the_search_s_record(self):
        out, rep = self.search_out(CLEAN + COMPETED)
        subprocess.run([sys.executable, CHECK, "--search-out", out], capture_output=True,
                       timeout=120)
        fs, md = self.audit({"input": rep}, out)
        (f,) = fs
        self.assertEqual(f["status"], "WARN")
        self.assertIn("ACTB is missing from the results", f["plain"])
        self.assertIn("**contaminant_competition** — COMPETED OUT", md)
        self.assertNotIn("checked_by", json.dumps(f))
        self.assertIn("ACTB is missing", ar.client_section_md([f]))

    def test_a_search_without_a_record_is_checked_by_the_audit(self):
        out, rep = self.search_out(CLEAN + SIPHONED)
        (f,), _md = self.audit({"input": rep, "q_cutoff": 0.01})
        self.assertEqual(f["status"], "WARN")
        self.assertIn("ENO1", f["plain"])

    def test_a_record_of_another_report_is_not_used_for_this_one(self):
        out, rep = self.search_out(CLEAN)
        write(os.path.join(out, cc.RECORD), json.dumps(
            {"status": "ok", "report": os.path.join(out, "step3_assembly.parquet"),
             "competed_out": [{"gene": "WRONG"}]}))
        (f,), _md = self.audit({"input": rep}, out)
        self.assertEqual(f["status"], "PASS")

    def test_a_check_that_could_not_run_is_in_the_audit_with_its_reason(self):
        gone = os.path.join(self.d, "gone", "report.parquet")
        (f,), md = self.audit({"input": gone})
        self.assertEqual(f["status"], "INFO")
        self.assertIn("Not run: the report", f["message"])
        self.assertIn("ℹ️ **contaminant_competition** — Not run", md)

    def test_another_engine_s_search_is_not_asked(self):
        out, rep = self.search_out(CLEAN + COMPETED, {"engine": "sage"})
        fs, _md = self.audit({"input": rep}, None)
        self.assertEqual(fs, [])

    def test_a_record_that_could_not_run_is_tried_again_by_the_audit(self):
        out, rep = self.search_out(CLEAN + COMPETED)
        write(os.path.join(out, cc.RECORD), json.dumps(
            {"status": "not_run", "report": rep, "reason": "pyarrow is not installed"}))
        (f,), _md = self.audit({"input": rep}, out)
        self.assertEqual(f["status"], "WARN")
        self.assertIn("ACTB", f["plain"])


class TheJobsRunIt(unittest.TestCase):
    """Every DIA-NN search's last job ends with the check; it never fails the job."""

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


class TheWatcherSaysIt(Fixture):
    WATCH = os.path.join(SCRIPTS, "watch_run.sh")

    def watch(self, out):
        bindir = os.path.join(self.d, "bin")
        for name, body in (("sacct", 'echo "  COMPLETED"'), ("squeue", "true")):
            write(os.path.join(bindir, name), f"#!/bin/bash\n{body}\n")
            os.chmod(os.path.join(bindir, name), 0o755)
        env = dict(os.environ, PATH=bindir + os.pathsep + os.environ.get("PATH", ""))
        r = subprocess.run(["bash", self.WATCH, "--slurm", "999", "--out", out],
                           capture_output=True, text=True, env=env, timeout=60)
        self.assertEqual(r.returncode, 0, r.stderr)
        return json.loads(r.stdout)

    def test_a_finished_search_names_the_genes(self):
        out = os.path.join(self.d, "out")
        os.makedirs(out)
        cc.write(self.run_check(CLEAN + COMPETED), os.path.join(out, cc.RECORD))
        got = self.watch(out)
        self.assertEqual(got["contaminant_competition"]["competed_out"], ["ACTB"])
        self.assertEqual(got["contaminant_competition"]["drop_contaminant"], ["Cont_P60712"])
        self.assertTrue(any("COMPETED OUT" in w for w in got["warnings"]))
        self.assertFalse(got["failed"])

    def test_no_record_no_key(self):
        os.makedirs(os.path.join(self.d, "out"))
        self.assertNotIn("contaminant_competition", self.watch(os.path.join(self.d, "out")))


@unittest.skipUnless(pq, "pyarrow is not installed")
class TheParquetReader(Fixture):
    def test_parquet_and_tsv_give_the_same_record(self):
        rows = CLEAN + COMPETED + SIPHONED + MASKED + IL_TWIN
        failing = [("P60709", "EEEEEEEK", 50)]
        tsv = write_tsv(os.path.join(self.d, "report.tsv"), rows, failing)
        data = list(zip(*report_rows(rows, failing=failing)))
        par = os.path.join(self.d, "report.parquet")
        pq.write_table(pa.table({c: list(v) for c, v in zip(COLS, data)}), par)
        a, b = cc.check(tsv, self.fasta, self.meta), cc.check(par, self.fasta, self.meta)
        for k in ("competed_out", "siphoned", "masked", "drop_contaminant",
                  "contaminant_share_pct", "contaminants_showing_themselves"):
            self.assertEqual(a[k], b[k], k)
        self.assertEqual([c["gene"] for c in b["competed_out"]], ["ACTB", "TGT"])

    def test_an_adapted_report_is_not_dianns_question(self):
        par = os.path.join(self.d, "adapted.parquet")
        t = pa.table({c: [v] for c, v in zip(COLS, report_rows(CLEAN)[0])})
        pq.write_table(t.replace_schema_metadata({b"delimp.quantity.level": b"peptide"}), par)
        self.assertEqual(cc.check(par, self.fasta, self.meta)["status"], "not_applicable")


if __name__ == "__main__":
    unittest.main(verbosity=2)
