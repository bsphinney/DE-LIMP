#!/usr/bin/env python3
"""A contaminant protein group that rests only on the sample's own peptides is flagged.

Staff report (2026-10-14, 23 human heart dia-PASEF runs, fetch_fasta.py --proteome UP000005640
--contaminants universal, DIA-NN 2.7.0): the Universal set's horse myoglobin, Cont_P68082
(MYG_HORSE), was reported in every run, above human myoglobin, and its precursors were removed as
contaminant. Its group rested on GHHEAELKPLAQSHATK -- human MB's GHHEAEIKPLAQSHATK with I/L, the
same mass and fragments -- and ASEDLKK, which human MB (P02144) shares. A first-draft report read
it as an added horse-myoglobin standard.

What the build does, and does not, decide (fetch_fasta.py): its peptide rule already compares I and
L as equal (il_equal), and drops a contaminant entry with fewer than MIN_UNIQUE_PEPTIDES peptides of
its own. Horse myoglobin has 9 (checked on the UniProt sequences below), so it stays -- rightly for
a sample that does hold horse myoglobin. Which peptides a sample shows only the report can say:
audit_results.py's contaminant_evidence check reads the DE's input report (rows kept as run_de.R
keeps them) and flags each contaminant group whose every peptide is a target protein's -- shared,
or the same with I = L. "A target's peptide" is one the search's digest makes from it -- the build's
own rule (fetch_fasta.target_digest_peptides) -- never a string that merely occurs inside it
(2.12 statistics review: BSA's LCVLHEK sits in human ALB after Q).
"""
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fetch_fasta as ff  # noqa: E402

try:
    import pyarrow  # noqa: F401
    HAVE_PYARROW = True
except ImportError:
    HAVE_PYARROW = False

# UniProt REST (rest.uniprot.org/uniprotkb/<acc>.fasta), fetched 2026-10-15.
MYG_HUMAN = ("MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE"
             "DLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKH"
             "PGDFGADAQGAMNKALELFRKDMASNYKELGFQG")
MYG_HORSE = ("MGLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASE"
             "DLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKH"
             "PGDFGADAQGAMTKALELFRNDIAAKYKELGFQG")

# UniProt P02768 (human ALB) and P02769 (BSA, as the Universal set holds it), 70 residues around
# BSA's LCVLHEK: after Q in human ALB (not a trypsin site), after R in BSA.
ALB_WINDOW = "LGKVGSKCCKHPEAKRMPCAEDYLSVVLNQLCVLHEKTPVSDRVTKCCTESLVNRRPCFSALEVDETYVP"
BSA_WINDOW = "LGKVGTRCCTKPESERMPCTEDYLSLILNRLCVLHEKTPVSEKVTKCCTESLVNRRPCFSALTPDETYVP"

# Synthetic: a target protein, and a contaminant with peptides of its own (MCDEQWNYK,
# FNVWDTPGQEK, YLEFISDAIIHVLHSR), an I/L twin of the target's (GHHEAELK...) and a shared one.
TARGET = "MSDEQWHYKGHHEAEIKPLAQSHATKASEDLKKFNVWDTAGQEKYLEFISEAIIHVLHSR"
CONT = "MCDEQWNYKGHHEAELKPLAQSHATKASEDLKKFNVWDTPGQEKYLEFISDAIIHVLHSR"
TRYP = "IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYK"
FASTA = (f">sp|Q90001|TGT_HUMAN Target protein OS=Homo sapiens OX=9606 GN=TGT PE=1 SV=1\n{TARGET}\n"
         f">sp|Cont_Q90002|TGT_HORSE Target protein OS=Equus caballus OX=9796 GN=TGT PE=1 SV=1\n{CONT}\n"
         f">sp|Cont_P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1\n{TRYP}\n")


class IlEqual(unittest.TestCase):
    def test_one_reading_of_i_and_l(self):
        self.assertEqual(ff.il_equal("GHHEAEIKPLAQSHATK"), ff.il_equal("GHHEAELKPLAQSHATK"))
        self.assertEqual(ff.il_equal("aIk"), "ALK")

    def test_contaminant_accessions(self):
        for acc in ("Cont_P68082", "contam_P68082", " Cont_P00761"):
            self.assertTrue(ff.is_contaminant_accession(acc), acc)
        for acc in ("P02144", "XCont_P1", ""):
            self.assertFalse(ff.is_contaminant_accession(acc), acc)


class TheBuildKeepsAContaminantWithItsOwnPeptides(unittest.TestCase):
    """The rule is unchanged, and this is why the audit is the place for the check."""

    def drop(self, cont_seq, target_seq, name="X"):
        _, dropped, _ = ff.drop_target_contaminants(
            f">sp|Cont_{name}|C_HORSE c\n{cont_seq}\n", f">sp|T1|T_HUMAN t GN=T\n{target_seq}\n")
        return dropped

    def test_horse_myoglobin_has_nine_peptides_of_its_own_and_stays(self):
        (p,) = ff.peptide_overlap(ff._fasta_records(f">sp|Cont_P68082|MYG_HORSE\n{MYG_HORSE}\n"),
                                  ff._fasta_records(f">sp|P02144|MYG_HUMAN\n{MYG_HUMAN}\n")).values()
        self.assertEqual(p["n_unique_peptides"], 9)
        self.assertIn("GHHEAELKPLAQSHATK", p["shared_peptides"])      # shared, with I = L
        self.assertIn("ASEDLKK", p["shared_peptides"])
        self.assertEqual(self.drop(MYG_HORSE, MYG_HUMAN, "P68082"), [])

    def test_the_synthetic_contaminant_stays_and_its_il_only_twin_goes(self):
        self.assertEqual(self.drop(CONT, TARGET), [])
        (d,) = self.drop(ff.il_equal(TARGET), TARGET)          # differs from TARGET by I/L alone
        self.assertEqual((d["reason"], d["n_unique_peptides"]), ("shared_peptides", 0))


class GroupsOnTargetPeptides(unittest.TestCase):
    def flagged(self, groups, fasta=FASTA):
        return {g["group"]: g for g in ff.contaminant_groups_on_target_peptides(groups, fasta)}

    def test_an_il_twin_and_a_shared_peptide_flag_the_group(self):
        f = self.flagged({"Cont_Q90002": {"GHHEAELKPLAQSHATK": {"Cont_Q90002"},
                                          "ASEDLKK": {"Cont_Q90002", "Q90001"}}})
        g = f["Cont_Q90002"]
        self.assertEqual(g["cont_entries"], ["TGT_HORSE"])
        by = {p["peptide"]: p for p in g["peptides"]}
        self.assertEqual((by["GHHEAELKPLAQSHATK"]["how"], by["GHHEAELKPLAQSHATK"]["target_acc"],
                          by["GHHEAELKPLAQSHATK"]["target_peptide"], by["GHHEAELKPLAQSHATK"]["gene"]),
                         ("il_twin", "Q90001", "GHHEAEIKPLAQSHATK", "TGT"))
        self.assertEqual((by["ASEDLKK"]["how"], by["ASEDLKK"]["target_acc"]), ("shared", "Q90001"))

    def test_a_peptide_of_its_own_clears_the_group(self):
        self.assertEqual(self.flagged({"Cont_Q90002": {"GHHEAELKPLAQSHATK": {"Cont_Q90002"},
                                                       "FNVWDTPGQEK": {"Cont_Q90002"}}}), {})
        self.assertEqual(self.flagged({"Cont_P00761": {"IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYK":
                                                       {"Cont_P00761"}}}), {})

    def test_a_target_digest_peptide_counts_even_when_the_report_names_only_the_contaminant(self):
        # ASEDLKK is one of TARGET's digest peptides (after K, one missed cleavage)
        g = self.flagged({"Cont_Q90002": {"ASEDLKK": {"Cont_Q90002"}}})["Cont_Q90002"]
        self.assertEqual(g["peptides"][0]["how"], "shared")

    def test_a_string_inside_a_target_is_not_its_peptide(self):
        """2.12 statistics review: the audit matched raw substrings while the build matches the
        digest, so BSA's own LCVLHEK -- after R in BSA, after Q in human ALB -- was called ALB's,
        and the fix offered (--drop-contaminant BSA) would have moved bovine signal into human
        ALB. The windows are UniProt P02768 / P02769 around it, cleavage context intact."""
        fasta = (f">sp|P02768|ALBU_HUMAN Albumin OS=Homo sapiens OX=9606 GN=ALB\n{ALB_WINDOW}\n"
                 f">sp|Cont_P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB\n{BSA_WINDOW}\n")
        self.assertIn("LNQLCVLHEK", ALB_WINDOW)
        self.assertIn("LNRLCVLHEK", BSA_WINDOW)
        self.assertEqual(self.flagged({"Cont_P02769": {"LCVLHEK": {"Cont_P02769"}}}, fasta), {})
        # the same string after a K in the target IS its peptide
        alb_k = ALB_WINDOW.replace("LNQLCVLHEK", "LNKLCVLHEK")
        fasta_k = fasta.replace(ALB_WINDOW, alb_k)
        g = self.flagged({"Cont_P02769": {"LCVLHEK": {"Cont_P02769"}}}, fasta_k)["Cont_P02769"]
        self.assertEqual((g["peptides"][0]["how"], g["peptides"][0]["target_acc"]), ("shared", "P02768"))

    def test_no_own_peptide_of_the_build_is_called_a_target_s(self):
        """One definition: every peptide the build counts as the contaminant's own (peptide_overlap)
        is, to the audit, not a target's -- measured 0 of 7,464 on the Universal set vs human
        Swiss-Prot after the fix (99 of 221 kept entries had one called human's before)."""
        cont, targ = ff._fasta_records(FASTA)[1:2], ff._fasta_records(FASTA)[:1]
        shared = {ff.il_equal(p) for p in ff.peptide_overlap(cont, targ)[0]["shared_peptides"]}
        seq = ff._record_seq(cont[0][1])
        own = [seq[s:e] for s, e in ff.digest_spans(seq) if ff.il_equal(seq[s:e]) not in shared]
        self.assertTrue(own)
        for pep in own:
            self.assertEqual(self.flagged({"Cont_Q90002": {pep: {"Cont_Q90002"}}}), {}, pep)
        for pep in shared:
            self.assertIn("Cont_Q90002", self.flagged({"Cont_Q90002": {pep: {"Cont_Q90002"}}}), pep)

    def test_a_peptide_shared_only_with_another_contaminant_is_not_the_sample_s(self):
        fasta = FASTA + f">sp|Cont_Q90003|OTHER_BOVIN o\nMKKWWPEPTIDEWWK\n"
        self.assertEqual(self.flagged({"Cont_Q90002": {"WWPEPTIDEWWK": {"Cont_Q90002", "Cont_Q90003"}}},
                                      fasta), {})

    def test_a_mixed_group_is_flagged_too(self):
        f = self.flagged({"Cont_Q90002;Q90001": {"ASEDLKK": {"Cont_Q90002", "Q90001"}}})
        self.assertIn("Cont_Q90002;Q90001", f)

    def test_the_real_myoglobin_case(self):
        fasta = (f">sp|P02144|MYG_HUMAN Myoglobin OS=Homo sapiens OX=9606 GN=MB PE=1 SV=2\n{MYG_HUMAN}\n"
                 f">sp|Cont_P68082|MYG_HORSE Myoglobin OS=Equus caballus OX=9796 GN=MB PE=1 SV=2\n"
                 f"{MYG_HORSE}\n")
        g = self.flagged({"Cont_P68082": {"GHHEAELKPLAQSHATK": {"Cont_P68082"},
                                          "ASEDLKK": {"Cont_P68082", "P02144"}}}, fasta)["Cont_P68082"]
        self.assertEqual(sorted((p["peptide"], p["how"], p["target_acc"]) for p in g["peptides"]),
                         [("ASEDLKK", "shared", "P02144"),
                          ("GHHEAELKPLAQSHATK", "il_twin", "P02144")])
        # one horse-only peptide observed would have been evidence of horse myoglobin
        self.assertEqual(self.flagged({"Cont_P68082": {"GHHEAELKPLAQSHATK": {"Cont_P68082"},
                                                       "HGTVVLTALGGILK": {"Cont_P68082"}}}, fasta), {})


def report_rows():
    """Precursor rows as DIA-NN writes them: the twin group in 3 runs, trypsin in 3, and an own
    peptide of the twin group that fails the q-value filter (so the DE never kept it)."""
    rows = []
    for run in ("r1", "r2", "r3"):
        rows += [(run, "Cont_Q90002", "Cont_Q90002", "GHHEAELKPLAQSHATK", 0.001),
                 (run, "Cont_Q90002", "Cont_Q90002;Q90001", "ASEDLKK", 0.002),
                 (run, "Cont_P00761", "Cont_P00761", "IVGGYTCAANSIPYQVSLNSGSHFCGGSLINSQWVVSAAHCYK", 0.001),
                 (run, "Q90001", "Q90001", "FNVWDTAGQEK", 0.001)]
    rows.append(("r1", "Cont_Q90002", "Cont_Q90002", "FNVWDTPGQEK", 0.2))
    return rows


class AuditFlagsIt(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.fasta = os.path.join(self.root, "search.fasta")
        with open(self.fasta, "w") as fh:
            fh.write(FASTA)
        with open(self.fasta, "rb") as fh:
            sha = hashlib.sha256(fh.read()).hexdigest()
        self.meta = self.fasta + ".meta.json"
        with open(self.meta, "w") as fh:
            json.dump({"fasta": self.fasta, "sha256": sha, "organism": "Homo sapiens",
                       "taxid": 9606}, fh)
        self.de = os.path.join(self.root, "de")
        os.makedirs(self.de)

    def tearDown(self):
        self._td.cleanup()

    def write_report(self, kind):
        cols = ("Run", "Protein.Group", "Protein.Ids", "Stripped.Sequence", "Q.Value")
        rows = report_rows()
        if kind.startswith("parquet"):
            import pyarrow as pa
            import pyarrow.parquet as pq
            path = os.path.join(self.root, "report.parquet")
            t = pa.table({c: [r[i] for r in rows] for i, c in enumerate(cols)})
            if kind == "parquet_dict":          # another writer's dictionary-encoded strings
                t = t.set_column(2, "Protein.Ids", t["Protein.Ids"].dictionary_encode())
            pq.write_table(t, path)
        else:
            path = os.path.join(self.root, "report.tsv")
            with open(path, "w") as fh:
                fh.write("\t".join(cols) + "\n")
                fh.writelines("\t".join(str(x) for x in r) + "\n" for r in rows)
        return path

    def audit(self, report, removed=True, n_groups=2, meta=True):
        with open(os.path.join(self.de, "de_provenance.json"), "w") as fh:
            json.dump({"input": report, "q_cutoff": 0.01, "q_columns": ["Q.Value"],
                       "contaminants": {"n_protein_groups": n_groups, "removed": removed,
                                        "tag": "Cont_"}}, fh)
        out = os.path.join(self.root, "AUDIT.md")
        args = [sys.executable, os.path.join(SCRIPTS, "audit_results.py"), "--out", out,
                "--de-dir", self.de] + (["--fasta-meta", self.meta] if meta else [])
        r = subprocess.run(args, capture_output=True, text=True, cwd=self.root, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(self.root, "AUDIT.json")) as fh:
            found = [f for f in json.load(fh)["findings"] if f["check"] == "contaminant_evidence"]
        with open(out, encoding="utf-8") as fh:
            return found, fh.read()

    def check_warn(self, report):
        (f,), md = self.audit(report)
        self.assertEqual(f["status"], "WARN")
        msg = f["message"]
        self.assertIn("1 contaminant protein group(s) rest only on peptides of Homo sapiens "
                      "proteins, so the data cannot attribute them to the contaminant's species: "
                      "Cont_Q90002 (TGT_HORSE; 2 peptides in 3 run(s)): ASEDLKK is also Q90001 TGT's, "
                      "GHHEAELKPLAQSHATK differs from Q90001 TGT's GHHEAEIKPLAQSHATK only by I/L "
                      "(the same mass and fragments).", msg)
        self.assertIn("removed with the contaminants", msg)
        # the data cannot say whether it is the sample's own or handling contamination (a human
        # keratin behind a sheep keratin entry): neither is asserted (2.12 statistics review)
        self.assertIn("Report each as the Homo sapiens protein -- the sample's own, or handling "
                      "contamination (skin and hair keratins) -- never as a signal of the "
                      "contaminant species or an added standard.", msg)
        self.assertNotIn("probably the sample's own protein", msg)
        self.assertEqual(f["detail"]["digest_source"], "estimate_params.DIANN_DIGEST (the sidecar "
                                                       "records no contaminant_digest)")
        self.assertTrue(msg.endswith("rebuild the FASTA with the same fetch_fasta.py command plus "
                                     "--drop-contaminant Cont_Q90002 and search again."), msg)
        self.assertIn("To quantify it as the Homo sapiens protein, rebuild", msg)
        self.assertEqual(f["detail"]["drop_contaminant"], ["Cont_Q90002"])
        self.assertNotIn("Cont_P00761", msg)                 # trypsin has its own peptide
        self.assertNotIn(self.root, msg)
        self.assertEqual(f["detail"]["n_groups_checked"], 2)
        self.assertIn("contaminant_evidence", md)

    @unittest.skipUnless(HAVE_PYARROW, "needs pyarrow")
    def test_parquet_report(self):
        self.check_warn(self.write_report("parquet"))

    @unittest.skipUnless(HAVE_PYARROW, "needs pyarrow")
    def test_dictionary_encoded_parquet_report(self):
        self.check_warn(self.write_report("parquet_dict"))

    def test_tsv_report(self):
        self.check_warn(self.write_report("tsv"))

    @unittest.skipUnless(HAVE_PYARROW, "needs pyarrow")
    def test_kept_contaminants_are_said_to_be_quantified(self):
        (f,), _ = self.audit(self.write_report("parquet"), removed=False)
        self.assertIn("quantified and tested as contaminant groups (--keep-contaminants)", f["message"])

    def test_the_digest_is_the_one_the_sidecar_records(self):
        """GHHEAEIKPLAQSHATK is a TARGET peptide only with a missed cleavage (GHHEAEIK|PLAQSHATK):
        a build recorded with --missed-cleavages 0 never makes it, so the twin is the
        contaminant's own there and the group is not flagged."""
        with open(self.meta) as fh:
            m = json.load(fh)
        m["contaminant_digest"] = dict(ff.DIANN_DIGEST, missed_cleavages=0)
        with open(self.meta, "w") as fh:
            json.dump(m, fh)
        (f,), _ = self.audit(self.write_report("tsv"))
        self.assertEqual(f["status"], "PASS", f["message"])
        self.assertEqual(f["detail"]["digest_source"], "the sidecar's contaminant_digest")
        self.assertEqual(f["detail"]["digest"]["missed_cleavages"], 0)

    def test_an_own_peptide_that_passed_the_filter_clears_it(self):
        path = self.write_report("tsv")
        with open(path) as fh:
            text = fh.read()
        with open(path, "w") as fh:
            fh.write(text.replace("FNVWDTPGQEK\t0.2", "FNVWDTPGQEK\t0.001"))
        (f,), _ = self.audit(path)
        self.assertEqual(f["status"], "PASS", f["message"])

    def test_not_assessed_is_said(self):
        os.remove(self.fasta)
        (f,), _ = self.audit(self.write_report("tsv"))
        self.assertEqual(f["status"], "INFO")
        self.assertIn("Not assessed: the searched FASTA (search.fasta, sha256", f["message"])
        self.assertIn("is not readable here", f["message"])

    def test_no_contaminant_group_no_finding(self):
        found, _ = self.audit(self.write_report("tsv"), n_groups=0)
        self.assertEqual(found, [])


if __name__ == "__main__":
    unittest.main()
