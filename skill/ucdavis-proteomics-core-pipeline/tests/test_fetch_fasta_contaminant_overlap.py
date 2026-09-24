#!/usr/bin/env python3
"""A contaminant entry that IS a target protein must not take that protein out of quant.

HIVE, 2026-09-24: the Universal contaminant set (Frankenfield 2022, 381 entries) carries 153
sequences identical to human UP000005640 entries (human keratins/KRTAPs/FLG, and bovine ACTB =
human P60709, EEF1A1 = P68104, TUBB5 = P07437 ...; plus Cont_Q3LI67 as an exact substring) and
31 identical to mouse. In a real DIA-NN 2.7.0 HeLa search ACTB, EEF1A1 and KRT8 came out ONLY
as Cont_ protein groups, and --cont-quant-exclude Cont_ kept them out of quantification.

Guards, all offline on tiny synthetic FASTAs:
  * identical and contained contaminant entries are dropped and recorded; near-identical,
    too-short and I/L-only-different ones are kept; a target entry is never dropped;
  * the sidecar lists them, the counts stay truthful, the warning is present, and the
    methods text says what happened;
  * a database used as-is (contaminants already inside) cannot be fixed -- it warns loudly;
  * audit_results.py and sample_quality.py flag a matched protein as possible contamination
    KEPT in quantification, reading the list from the sidecar.
"""
import argparse
import contextlib
import csv
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fetch_fasta as ff  # noqa: E402

# Synthetic stand-ins, shaped like the real headers (>sp|Cont_P60712|ACTB_BOVIN ... GN=ACTB).
ACTB = "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEK"
KRT8 = "MSIRVTQKSYKVSTSGPRAFSSRSYTSGPGSRISSSSFSRVGSSNFRGGLGGGYGGASGMGGITAVTVNQSLLSPLVLEVDPNIQAVRTQEKEQIKTLNNKFASFIDKVRFLEQQNKMLETKWSLLQQQKTARSNMDNMFESYINNLRRQLETLGQEKLKLEAELGNMQGLVEDFKNKYEDEINKRTEMENEFVLIKKDVDEAYMNKVELESRLEGLTDEINFLRQLYEEEIRELQSQISDTSVVLSMDNSRSLDMDSIIAEVKAQYEDIANRSRAEAESMYQIKYEELQSLAGKHGDDLRRTKTEISEMNRNISRLQAEIEGLKGQRASLEAAIADAEQRGELAIKDANAKLSELEAALQRAKQDMARQLREYQELMNVKLALDIEIATYRKLLEGEESRLESGMQNMSIHTKTTSGYAGGLSSAYGGLTSPGLSYSLGSSFGSGAGSSSFSRTSSSRAVVVKKIETRDGKLVSESSDVLPK"
GAPDH = "MGKVKVGVNGFGRIGRLVTRAAFNSGKVDIVAINDPFIDLNYMVYMFQYDSTHGKFHGTVKAENGKLVINGNPITIFQERDPSKIKWGDAGAEYVVESTGVFTTMEKAGAHLQGGAKRVIISAPSADAPMFVMGVNHEKYDNSLKIISNASCTTNCLAPLAKVIHDNFGIVEGLMTTVHAITATQKTVDGPSGKLWRDGRGALQNIIPASTGAAKAVGKVIPELNGKLTGMAFRVPTANVSVVDLTCRLEKPAKYDDIKKVVKQASEGPLKGILGYTEDQVVSCDFNSDTHSSTFDAGAGIALNDHFVKLISWYDNEFGYSNRVVDLMAHMASKE"

TARGET = (f">sp|P60709|ACTB_HUMAN Actin, cytoplasmic 1 OS=Homo sapiens OX=9606 GN=ACTB PE=1 SV=1\n{ACTB}\n"
          f">sp|P05787|K2C8_HUMAN Keratin 8 OS=Homo sapiens OX=9606 GN=KRT8 PE=1 SV=7\n{KRT8}\n"
          f">sp|P04406|G3P_HUMAN GAPDH OS=Homo sapiens OX=9606 GN=GAPDH PE=1 SV=3\n{GAPDH}\n")

ACTB_ONE_OFF = ACTB[:10] + ("A" if ACTB[10] != "A" else "G") + ACTB[11:]
GAPDH_IL = GAPDH.replace("I", "L")         # an I/L-only difference: a different string to DIA-NN
assert GAPDH_IL != GAPDH

CONT = (f">sp|Cont_P60712|ACTB_BOVIN Actin, cytoplasmic 1 OS=Bos taurus OX=9913 GN=ACTB PE=1 SV=1\n"
        f"{ACTB[:40]}\n{ACTB[40:]}\n"                                  # identical, wrapped
        f">sp|Cont_Q3LI67|KRA63_HUMAN fragment OS=Homo sapiens OX=9606 GN=KRTAP6-3 PE=3 SV=3\n"
        f"{KRT8[50:120]}\n"                                              # exact substring of KRT8
        f">sp|Cont_Q00001|ACTB_OTHER one residue off OS=X GN=ACTB PE=1 SV=1\n{ACTB_ONE_OFF}\n"
        f">sp|Cont_Q00002|SHORT_X six residues OS=X GN=SHRT PE=1 SV=1\n{GAPDH[5:11]}\n"
        f">sp|Cont_Q00003|G3P_ILSWAP I->L only OS=X GN=GAPDH PE=1 SV=1\n{GAPDH_IL}\n"
        f">sp|Cont_P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1\nFPTDDDDKIVGGYTCAANSIPYQVSLNSG\n")


def headers(text):
    return [ln.split()[0] for ln in text.splitlines() if ln.startswith(">")]


class Matching(unittest.TestCase):
    def test_identical_and_substring_dropped_everything_else_kept(self):
        kept, dropped, enzymes = ff.drop_target_contaminants(CONT, TARGET)
        self.assertEqual(enzymes, [])      # the trypsin entry matches no target here
        by = {r["cont_acc"]: r for r in dropped}
        self.assertEqual(set(by), {"Cont_P60712", "Cont_Q3LI67"})
        self.assertEqual(by["Cont_P60712"]["reason"], "identical")
        self.assertEqual((by["Cont_P60712"]["target_acc"], by["Cont_P60712"]["gene"]),
                         ("P60709", "ACTB"))
        self.assertEqual(by["Cont_Q3LI67"]["reason"], "substring")
        # `gene` is the TARGET's -- the name the protein group now carries in the results.
        self.assertEqual((by["Cont_Q3LI67"]["target_acc"], by["Cont_Q3LI67"]["gene"]),
                         ("P05787", "KRT8"))
        self.assertEqual(by["Cont_Q3LI67"]["cont_gene"], "KRTAP6-3")
        # Near-identical, below the minimum length, I/L-only, and unrelated: all kept, and
        # written back byte-for-byte.
        self.assertEqual(headers(kept), [">sp|Cont_Q00001|ACTB_OTHER", ">sp|Cont_Q00002|SHORT_X",
                                         ">sp|Cont_Q00003|G3P_ILSWAP", ">sp|Cont_P00761|TRYP_PIG"])
        self.assertIn(f"{GAPDH_IL}\n", kept)

    def test_substring_minimum_is_seven_residues(self):
        self.assertEqual(ff.MIN_CONTAINED_LEN, 7)     # DIA-NN's default --min-pep-len
        tgt = ff._fasta_records(TARGET)
        six = ff._fasta_records(f">sp|Cont_A|X_Y\n{GAPDH[5:11]}\n")
        seven = ff._fasta_records(f">sp|Cont_A|X_Y\n{GAPDH[5:12]}\n")
        self.assertEqual(ff.contaminants_matching_targets(six, tgt), [])
        self.assertEqual(len(ff.contaminants_matching_targets(seven, tgt)), 1)

    def test_every_target_that_contains_it_is_listed(self):
        tgt = TARGET + f">sp|P99999|ACTB2_HUMAN dup OS=Homo sapiens GN=ACTB2 PE=1 SV=1\n{ACTB}\n"
        _, dropped, _ = ff.drop_target_contaminants(CONT, tgt)
        actb = next(r for r in dropped if r["cont_acc"] == "Cont_P60712")
        self.assertEqual(actb["n_targets"], 2)
        self.assertEqual(actb["target_accs"], ["P60709", "P99999"])


class FetchSidecar(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.target = os.path.join(self.root, "target.fasta")
        self.cont = os.path.join(self.root, "cont.fasta")
        with open(self.target, "w") as fh:
            fh.write(TARGET)
        with open(self.cont, "w") as fh:
            fh.write(CONT)

    def tearDown(self):
        self._td.cleanup()

    def fetch(self, *args):
        out = os.path.join(self.root, "out", "search.fasta")
        argv = ["fetch_fasta.py", "fetch", *args, "--out", out]
        err = io.StringIO()
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(err):
            self.assertEqual(ff.main(), 0, err.getvalue())
        with open(out + ".meta.json") as fh:
            meta = json.load(fh)
        with open(out) as fh:
            return meta, fh.read(), err.getvalue()

    def test_sidecar_records_the_drop_and_counts_stay_truthful(self):
        m, fasta, err = self.fetch("--path", self.target, "--contaminants", "universal",
                                   "--contaminants-path", self.cont)
        self.assertEqual(m["n_contaminants_in_set"], 6)
        self.assertEqual(m["n_contaminants_appended"], 4)
        self.assertEqual(m["n_contaminants_dropped_as_target"], 2)
        self.assertEqual({r["cont_acc"] for r in m["contaminants_dropped_as_target"]},
                         {"Cont_P60712", "Cont_Q3LI67"})
        self.assertEqual(m["n_sequences"], 3 + 4)
        self.assertEqual(m["n_entries"], 3 + 4)
        self.assertEqual(m["diann_cont_quant_exclude"], "Cont_")
        self.assertEqual(m["contaminants_identical_to_target_kept"], [])
        # The warning is there, and it IS the recorded note (make_methods keys on that).
        self.assertIn(m["contaminants_dropped_note"], m["warnings"])
        self.assertIn("removed 2 of the 6", m["contaminants_dropped_note"])
        self.assertIn("Cont_P60712 (ACTB_BOVIN) = P60709 ACTB", m["contaminants_dropped_note"])
        self.assertIn("removed 2 of the 6", err)
        # The searched database: the dropped contaminants are gone ...
        self.assertNotIn(">sp|Cont_P60712|", fasta)
        self.assertNotIn(">sp|Cont_Q3LI67|", fasta)

    def test_target_entries_are_never_dropped(self):
        _, fasta, _ = self.fetch("--path", self.target, "--contaminants", "universal",
                                 "--contaminants-path", self.cont)
        # ... and every target entry is still there, sequence intact.
        self.assertTrue(fasta.startswith(TARGET))
        self.assertEqual([h for h in headers(fasta) if "Cont_" not in h],
                         [">sp|P60709|ACTB_HUMAN", ">sp|P05787|K2C8_HUMAN", ">sp|P04406|G3P_HUMAN"])

    def test_staged_hive_source_names_the_organism(self):
        mrs = os.path.join(self.root, "MRS")
        os.makedirs(mrs)
        with open(os.path.join(mrs, "UP000005640_9606.fasta"), "w") as fh:
            fh.write(TARGET)
        data = {"taxonomy": {"scientificName": "Homo sapiens", "taxonId": 9606},
                "proteomeType": "Reference proteome", "geneCount": 3, "superkingdom": "eukaryota"}
        with mock.patch.object(ff, "HIVE_MRS", mrs), \
                mock.patch.object(ff, "_get_json", return_value=(data, {})):
            m, _, _ = self.fetch("--proteome", "UP000005640", "--hive", "--contaminants",
                                 "universal", "--contaminants-path", self.cont)
        self.assertEqual(m["n_contaminants_dropped_as_target"], 2)
        self.assertIn("is a Homo sapiens protein", m["contaminants_dropped_note"])

    def test_no_overlap_leaves_everything_as_before(self):
        clean = os.path.join(self.root, "clean.fasta")
        with open(clean, "w") as fh:
            fh.write(">sp|Cont_P00761|TRYP_PIG Trypsin\nFPTDDDDKIVGGYTCAANSIPYQVSLNSG\n")
        m, _, _ = self.fetch("--path", self.target, "--contaminants", "universal",
                             "--contaminants-path", clean)
        self.assertEqual(m["n_contaminants_dropped_as_target"], 0)
        self.assertEqual(m["contaminants_dropped_as_target"], [])
        self.assertIsNone(m["contaminants_dropped_note"])
        self.assertEqual(m["n_contaminants_appended"], 1)
        self.assertEqual(m["warnings"], [])

    def test_database_with_contaminants_inside_warns_it_cannot_be_fixed(self):
        # The staged ..._plus_universal_contam.fasta case: >= 20 Cont_ entries already inside.
        aa = "ACDEFGHKLM"
        filler = "".join(f">sp|Cont_Z{i:05d}|FILL_X filler\n"
                         f"WWWWWWWWHHHHHHHH{''.join(aa[int(d)] for d in f'{i:04d}')}\n"
                         for i in range(20))
        combined = os.path.join(self.root, "combined.fasta")
        with open(combined, "w") as fh:
            fh.write(TARGET + CONT + filler)
        m, fasta, _ = self.fetch("--path", combined, "--contaminants", "universal")
        self.assertEqual(m["n_contaminants_appended"], 0)       # used as-is ...
        self.assertEqual(m["n_contaminants_dropped_as_target"], 0)  # ... so nothing removable
        kept = {r["cont_acc"] for r in m["contaminants_identical_to_target_kept"]}
        self.assertEqual(kept, {"Cont_P60712", "Cont_Q3LI67"})
        w = next(x for x in m["warnings"] if "cannot be removed" in x)
        self.assertIn("Cont_P60712 (ACTB_BOVIN) = P60709 ACTB", w)
        self.assertIn("without --path", w)
        self.assertEqual(fasta.rstrip("\n"), (TARGET + CONT + filler).rstrip("\n"))

    def test_methods_text_says_the_proteins_are_quantified(self):
        m, _, _ = self.fetch("--path", self.target, "--contaminants", "universal",
                             "--contaminants-path", self.cont)
        m["organism"] = "Homo sapiens"          # a --path build records none; name it for the text
        meta_path = os.path.join(self.root, "meta.json")
        with open(meta_path, "w") as fh:
            json.dump(m, fh)
        md = os.path.join(self.root, "methods.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                            "--raw", os.path.join(self.root, "FL_hela.raw"),
                            "--fasta-meta", meta_path, "--out", md],
                           capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(md) as fh:
            text = fh.read()
        self.assertIn("2 contaminant entries identical to (or contained in) Homo sapiens "
                      "proteins were removed from the library first, so those proteins are "
                      "quantified under their own accessions.", text)
        # Described in the sentence, so NOT repeated as something to resolve before publication.
        self.assertNotIn("Database build warnings", text)


PIG_TRYPSIN = "FPTDDDDKIVGGYTCAANSIPYQVSLNSG"        # the fixture's Cont_P00761 sequence
PIG_TARGET = TARGET + f">sp|P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1\n{PIG_TRYPSIN}\n"
# A synthetic stand-in. The real case: S. aureus NCTC 8325's own SspA IS the Glu-C
# contaminant entry Q2FZL2 (reference proteome UP000008816; UniProt API, 2026-09-24).
SSPA = "VILPNNDRHQITDTTNGHYAPVTYIQVEAPTGTFIASGVVVGKDTLLTNKHVVDATHGDPHALKAFPSAINQDNYPNGGFTAEQITKYSGEGDLAIVKFSPNEQNKHIGEVVKPATMSNNAETQVNQNITVTGYPGDKPVATMWESKGKITYLKGEAMQYDLSTTGGNSGSPVFNEKNEVIGIHWGGVPNEFNGAVFINENVRNFLKQNIEDIHFANDDQPNNPDNPDNPNNPDNPNNPDEPNNPDNPNNPDNPDNGDNNNSDNPDAA"
CONT_GLUC = CONT + f">sp|Cont_Q2FZL2|SSPA_STAA8 Glutamyl endopeptidase OS=Staphylococcus aureus GN=sspA PE=1 SV=1\n{SSPA}\n"
SAUREUS_TARGET = TARGET + f">sp|Q2FZL2|SSPA_STAA8 Glutamyl endopeptidase GN=sspA PE=1 SV=1\n{SSPA}\n"


class DigestionEnzymes(unittest.TestCase):
    """The one exception, keyed on the enzyme(s) USED (--enzyme, default trypsin,lysc): that
    enzyme stays Cont_ even when it IS a target protein (porcine trypsin on a pig search), so
    its autolysis peptides stay out of normalisation. Every other protease -- and every
    non-enzyme, here ACTB -- follows the rule."""

    def fetch(self, root, *args):
        out = os.path.join(root, "out", "search.fasta")
        with mock.patch.object(sys, "argv", ["fetch_fasta.py", "fetch", *args, "--out", out]), \
                contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(ff.main(), 0)
        with open(out + ".meta.json") as fh, open(out) as fa:
            return json.load(fh), fa.read()

    def build(self, target, cont, *args):
        with tempfile.TemporaryDirectory() as root:
            tgt, cp = os.path.join(root, "target.fasta"), os.path.join(root, "cont.fasta")
            with open(tgt, "w") as fh:
                fh.write(target)
            with open(cp, "w") as fh:
                fh.write(cont)
            return self.fetch(root, "--path", tgt, "--contaminants", "universal",
                              "--contaminants-path", cp, *args)

    def test_default_keeps_porcine_trypsin_on_a_pig_identical_target(self):
        m, fasta = self.build(PIG_TARGET, CONT)          # no --enzyme: trypsin,lysc
        self.assertEqual(m["digestion_enzymes_used"], ["lysc", "trypsin"])
        kept = m["contaminants_kept_despite_target_match"]
        self.assertEqual([(r["cont_acc"], r["target_acc"], r["reason"], r["match"]) for r in kept],
                         [("Cont_P00761", "P00761",
                           "digestion enzyme used in this search (trypsin)", "identical")])
        self.assertIn(">sp|Cont_P00761|TRYP_PIG", headers(fasta))       # still Cont_, excluded
        self.assertEqual(m["diann_cont_quant_exclude"], "Cont_")
        dropped = {r["cont_acc"] for r in m["contaminants_dropped_as_target"]}
        self.assertIn("Cont_P60712", dropped)             # identical non-enzyme: dropped
        self.assertNotIn("Cont_P00761", dropped)
        w = next(x for x in m["warnings"] if x.startswith("kept 1 digestion-enzyme"))
        self.assertIn("Cont_P00761 (trypsin) = P00761", w)
        self.assertIn("--enzyme", w)                      # how to change it if the digest differed

    def test_gluc_digest_keeps_gluc_and_drops_trypsin(self):
        both = PIG_TARGET + f">sp|Q2FZL2|SSPA_STAA8 SspA GN=sspA\n{SSPA}\n"
        m, fasta = self.build(both, CONT_GLUC, "--enzyme", "gluc")
        self.assertEqual(m["digestion_enzymes_used"], ["gluc"])
        self.assertEqual([r["cont_acc"] for r in m["contaminants_kept_despite_target_match"]],
                         ["Cont_Q2FZL2"])
        self.assertIn(">sp|Cont_Q2FZL2|SSPA_STAA8", headers(fasta))
        drop = {r["cont_acc"]: r for r in m["contaminants_dropped_as_target"]}
        self.assertEqual(drop["Cont_P00761"]["enzyme_not_used"], "trypsin")   # not the digest
        self.assertNotIn(">sp|Cont_P00761|TRYP_PIG", headers(fasta))

    def test_default_drops_gluc_identical_to_the_s_aureus_protein(self):
        # Trypsin-digested S. aureus: Q2FZL2 is the organism's own SspA and must be quantified.
        m, fasta = self.build(SAUREUS_TARGET, CONT_GLUC)
        self.assertEqual(m["contaminants_kept_despite_target_match"], [])
        drop = {r["cont_acc"]: r for r in m["contaminants_dropped_as_target"]}
        self.assertEqual((drop["Cont_Q2FZL2"]["target_acc"], drop["Cont_Q2FZL2"]["enzyme_not_used"]),
                         ("Q2FZL2", "gluc"))
        self.assertNotIn(">sp|Cont_Q2FZL2|", fasta)
        self.assertIn(">sp|Q2FZL2|SSPA_STAA8", headers(fasta))      # the target is untouched
        w = next(x for x in m["warnings"] if "not in --enzyme" in x)
        self.assertIn("Cont_Q2FZL2 (gluc)", w)
        # Its own warning, not the drop note -- so the methods draft still lists it to confirm.
        self.assertNotEqual(w, m["contaminants_dropped_note"])

    def test_unknown_enzyme_is_a_clean_error(self):
        err = io.StringIO()
        with tempfile.TemporaryDirectory() as root, \
                mock.patch.object(sys, "argv", ["fetch_fasta.py", "fetch", "--proteome", "UP000005640",
                                                "--enzyme", "trypsin,papain",
                                                "--out", os.path.join(root, "x.fasta")]), \
                contextlib.redirect_stderr(err), self.assertRaises(SystemExit) as cm:
            ff.main()
        self.assertEqual(cm.exception.code, 2)                   # argparse usage error, no traceback
        self.assertIn("unknown enzyme 'papain'", err.getvalue())
        self.assertIn("trypsin", err.getvalue())                 # the choices are listed
        self.assertNotIn("Traceback", err.getvalue())

    def test_enzyme_names_are_normalised(self):
        self.assertEqual(ff.parse_enzymes("Trypsin, Lys-C"), ("lysc", "trypsin"))
        self.assertEqual(ff.parse_enzymes("glu_c"), ("gluc",))
        with self.assertRaises(argparse.ArgumentTypeError):
            ff.parse_enzymes("")

    def test_mature_enzyme_inside_a_proenzyme_target_is_kept_too(self):
        # A proteome usually carries the preproenzyme; the contaminant is the mature chain.
        tgt = TARGET + f">tr|F1SRS2|F1SRS2_PIG Trypsinogen GN=PRSS1\nMKTFIFLALLGAAVAF{PIG_TRYPSIN}\n"
        _, _, enzymes = ff.drop_target_contaminants(CONT, tgt)
        self.assertEqual([(e["cont_acc"], e["match"]) for e in enzymes], [("Cont_P00761", "substring")])

    def test_the_enzyme_list_is_one_constant(self):
        self.assertEqual(ff.DIGESTION_ENZYMES["P00761"], "trypsin")
        self.assertEqual(ff.ENZYME_FAMILIES, tuple(sorted(set(ff.DIGESTION_ENZYMES.values()))))
        self.assertEqual(set(ff.ENZYME_FAMILIES), {"trypsin", "lysc", "gluc", "chymotrypsin",
                                                   "aspn", "argc", "lysn", "pepsin"})
        self.assertEqual(ff._digestion_enzyme("Cont_Q2FZL2"), "gluc")
        self.assertIsNone(ff._digestion_enzyme("Cont_P00698"))   # LYSC_CHICK is lysozyme, not Lys-C
        self.assertIsNone(ff._digestion_enzyme("Cont_P02769"))   # BSA: sample protein, follows the rule

    def test_enzyme_inside_a_supplied_database_is_not_reported_as_a_lost_protein(self):
        aa = "ACDEFGHKLM"
        filler = "".join(f">sp|Cont_Z{i:05d}|FILL_X filler\n"
                         f"WWWWWWWWHHHHHHHH{''.join(aa[int(d)] for d in f'{i:04d}')}\n"
                         for i in range(20))
        with tempfile.TemporaryDirectory() as root:
            combined = os.path.join(root, "combined.fasta")
            with open(combined, "w") as fh:
                fh.write(PIG_TARGET + CONT + filler)
            m, _ = self.fetch(root, "--path", combined, "--contaminants", "universal")
        self.assertEqual({r["cont_acc"] for r in m["contaminants_identical_to_target_kept"]},
                         {"Cont_P60712", "Cont_Q3LI67"})
        self.assertEqual([r["cont_acc"] for r in m["contaminants_kept_despite_target_match"]],
                         ["Cont_P00761"])


RUN_DE_ID_COLUMNS = ["Protein.Group", "Genes", "Protein.Names"]   # run_de.R:447-451


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


class Auditors(unittest.TestCase):
    """The auditors read the list from the sidecar; neither keeps its own copy."""

    SAMPLES = ["S1", "S2", "S3", "S4", "S5", "S6"]

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.de = os.path.join(self.root, "de_results")
        os.makedirs(self.de)
        rec = lambda c, t, g, why="identical": {           # noqa: E731
            "cont_acc": c, "cont_entry": "", "cont_gene": g, "target_acc": t, "target_entry": "",
            "gene": g, "reason": why, "n_targets": 1, "target_accs": [t]}
        self.meta = {"organism": "Homo sapiens",
                     "contaminants_dropped_as_target": [
                         rec("Cont_P60712", "P60709", "ACTB"),
                         rec("Cont_P04264", "P04264", "KRT1"),
                         rec("Cont_XP1", "XP_000001.1", "")],   # NCBI: accession only
                     "contaminants_identical_to_target_kept": []}
        # KRT1 is 4 units up in S6 only: a per-sample excess for sample_quality to see.
        # run_de.R's REAL column order: Protein.Group, Genes, Protein.Names, then samples (R's
        # merge puts the key first). A fixture with Genes first hid a panel bug for months.
        write_csv(os.path.join(self.de, "Expression_Matrix.csv"),
                  RUN_DE_ID_COLUMNS + self.SAMPLES,
                  [["P60709", "ACTB", "ACTB_HUMAN"] + [20] * 6,
                   ["P04264;Q7Z794", "KRT1;KRT77", "K2C1_HUMAN;K2C1B_HUMAN"]
                   + [10, 10, 10, 10, 10, 14],
                   ["XP_000001.1", "", ""] + [12] * 6,
                   ["P04406", "GAPDH", "G3P_HUMAN"] + [18] * 6])
        write_csv(os.path.join(self.de, "DE_limpa_B.vs.A.csv"),
                  ["Protein.Group", "Genes", "logFC", "adj.P.Val"],
                  [["P04264;Q7Z794", "KRT1;KRT77", "2.1", "0.001"],
                   ["P60709", "ACTB", "0.1", "0.9"],
                   ["P04406", "GAPDH", "1.5", "0.002"]])

    def tearDown(self):
        self._td.cleanup()

    def write_meta(self, meta=None):
        p = os.path.join(self.root, "search.fasta.meta.json")
        with open(p, "w") as fh:
            json.dump(meta or self.meta, fh)
        return p

    def audit(self, *extra):
        out = os.path.join(self.root, "AUDIT.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                            "--out", out, "--de-dir", self.de, *extra],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(self.root, "AUDIT.json")) as fh:
            return json.load(fh)

    def finding(self, audit, check):
        return next(f for f in audit["findings"] if f["check"] == check)

    def test_audit_flags_matched_genes_as_kept_in_quantification(self):
        f = self.finding(self.audit("--fasta-meta", self.write_meta()), "target_contaminants")
        self.assertEqual(f["status"], "WARN")
        self.assertIn("KEPT in quantification", f["message"])
        self.assertEqual(sorted(f["detail"]["quantified"]), ["ACTB", "KRT1", "XP_000001.1"])
        self.assertNotIn("GAPDH", f["detail"]["quantified"])
        # A significant DE hit that may be contamination is called out by name.
        self.assertEqual(f["detail"]["significant"], {"DE_limpa_B.vs.A.csv": ["KRT1"]})

    def test_audit_reads_the_default_sidecar_in_cwd(self):
        self.write_meta()                        # ./search.fasta.meta.json, no --fasta-meta
        f = self.finding(self.audit(), "target_contaminants")
        self.assertEqual(f["status"], "WARN")

    def test_keratin_sample_does_not_flag_keratin(self):
        f = self.finding(self.audit("--fasta-meta", self.write_meta(), "--keratin-sample"),
                         "target_contaminants")
        self.assertNotIn("KRT1", f["detail"]["quantified"])
        self.assertIn("ACTB", f["detail"]["quantified"])
        self.assertEqual(f["detail"]["significant"], {})

    def test_audit_without_a_sidecar_says_not_assessed(self):
        a = self.audit()
        f = self.finding(a, "target_contaminants")
        self.assertEqual(f["status"], "INFO")
        self.assertIn("Not assessed", f["message"])

    def test_audit_warns_when_real_proteins_sit_only_as_cont(self):
        meta = dict(self.meta, contaminants_dropped_as_target=[],
                    contaminants_identical_to_target_kept=self.meta["contaminants_dropped_as_target"][:1])
        f = self.finding(self.audit("--fasta-meta", self.write_meta(meta)), "contaminant_overlap")
        self.assertEqual(f["status"], "WARN")
        self.assertIn("MISSING", f["message"])
        self.assertIn("ACTB", f["message"])

    def test_sample_quality_builds_the_panel_from_the_sidecar(self):
        out = os.path.join(self.root, "SAMPLE_QUALITY.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                            "--matrix", os.path.join(self.de, "Expression_Matrix.csv"),
                            "--fasta-meta", self.write_meta(), "--out", out],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out.replace(".md", ".json")) as fh:
            sq = json.load(fh)
        p = sq["panels"]["CONTAMINANT_IDENTICAL"]
        self.assertTrue(p["kept_in_quantification"])
        self.assertEqual(p["n_panel_proteins"], 3)          # ACTB, KRT1 and the NCBI accession
        self.assertEqual(p["elevated_samples"], ["S6"])
        self.assertTrue(any("CONTAMINANT_IDENTICAL: elevated in S6" in f
                            and "KEPT in quantification" in f for f in sq["flags"]), sq["flags"])
        with open(out) as fh:
            self.assertIn("kept in quantification", fh.read())


    # -- review 2026-09-24 --------------------------------------------------------------------
    def sample_quality(self, *extra):
        out = os.path.join(self.root, "SAMPLE_QUALITY.md")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                            "--matrix", os.path.join(self.de, "Expression_Matrix.csv"),
                            "--out", out, *extra],
                           capture_output=True, text=True, cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(out.replace(".md", ".json")) as fh:
            return json.load(fh)

    def test_curated_panels_match_genes_in_run_de_column_order(self):
        """read_matrix used the FIRST id column -- Protein.Group in run_de.R's output -- so the
        gene-symbol panels were compared against accessions: HBB/HBA1 +5 log2 in S1 gave
        HEMOLYSIS 0 proteins and no flag."""
        flat = [18.0] * 6
        write_csv(os.path.join(self.de, "Expression_Matrix.csv"), RUN_DE_ID_COLUMNS + self.SAMPLES,
                  [["P68871", "HBB", "HBB_HUMAN"] + [23, 18, 18, 18, 18, 18],
                   ["P69905", "HBA1;HBA2", "HBA_HUMAN"] + [23, 18, 18, 18, 18, 18],
                   ["P04406", "GAPDH", "G3P_HUMAN"] + flat,
                   ["P60709", "ACTB", "ACTB_HUMAN"] + flat])
        sq = self.sample_quality()
        h = sq["panels"]["HEMOLYSIS"]
        self.assertEqual(h["n_panel_proteins"], 2)
        self.assertEqual(h["matched_genes"], ["HBA1", "HBA2", "HBB"])
        self.assertEqual(h["elevated_samples"], ["S1"])
        self.assertTrue(any(f.startswith("HEMOLYSIS: elevated in S1") for f in sq["flags"]),
                        sq["flags"])

    def legacy_meta(self, fasta_text=None, sha=None, **extra):
        """A sidecar written BEFORE the overlap check: contaminants, no contaminant_target_rule."""
        meta = {"organism": "Homo sapiens", "taxid": 9606, "contaminant_set": "universal",
                "n_contaminants_appended": 6, "fasta": os.path.join(self.root, "search.fasta"),
                "sha256": sha or "0" * 64, **extra}
        if fasta_text is not None:
            with open(meta["fasta"], "w") as fh:
                fh.write(fasta_text)
            meta["sha256"] = sha or hashlib.sha256(fasta_text.encode()).hexdigest()
        return meta

    def cont_only_matrix(self):
        # What an old database's results look like: ACTB present ONLY as its Cont_ twin.
        write_csv(os.path.join(self.de, "Expression_Matrix.csv"), RUN_DE_ID_COLUMNS + self.SAMPLES,
                  [["Cont_P60712", "ACTB", "ACTB_BOVIN"] + [20] * 6,
                   ["P04406", "GAPDH", "G3P_HUMAN"] + [18] * 6])

    def test_old_sidecar_is_rechecked_against_its_fasta_and_made_concrete(self):
        self.cont_only_matrix()
        meta = self.legacy_meta(TARGET + CONT)
        f = self.finding(self.audit("--fasta-meta", self.write_meta(meta)), "contaminant_overlap")
        self.assertEqual(f["status"], "WARN")
        self.assertTrue(f["detail"]["legacy_database"])
        self.assertIn("built before fetch_fasta.py checked", f["message"])
        self.assertIn("2 of its Cont_ entries are Homo sapiens proteins", f["message"])
        self.assertIn("Affected: ACTB, KRT8", f["message"])
        self.assertEqual(f["detail"]["seen_only_as_cont"], ["ACTB (Cont_P60712)"])
        self.assertIn("Rebuild the FASTA", f["message"])

    def test_old_sidecar_without_its_fasta_warns_with_the_measured_size(self):
        f = self.finding(self.audit("--fasta-meta", self.write_meta(self.legacy_meta())),
                         "contaminant_overlap")
        self.assertEqual(f["status"], "WARN")
        self.assertIn("no longer readable", f["message"])
        self.assertIn("~153 target-identical", f["message"])
        self.assertIn("probably missing from quantification", f["message"])

    def test_old_sidecar_whose_fasta_changed_is_not_trusted(self):
        meta = self.legacy_meta(TARGET + CONT, sha="f" * 64)
        f = self.finding(self.audit("--fasta-meta", self.write_meta(meta)), "contaminant_overlap")
        self.assertIn("has changed since the search", f["message"])
        self.assertEqual(f["detail"]["genes"], [])            # no list from a file we can't trust

    def test_unreadable_fasta_is_named_not_the_sidecar(self):
        """An OSError on the searched FASTA used to escape to the auditors' sidecar handler and
        be reported as "could not read <meta.json>" -- the wrong file."""
        meta = self.legacy_meta(TARGET + CONT)
        with mock.patch.object(ff, "_read_fasta_text", side_effect=PermissionError(13, "denied")):
            tc = ff.target_contaminants(meta)
        self.assertEqual(tc["kept_as_contaminant"], [])
        self.assertIn(f"the searched FASTA {meta['fasta']} could not be read", tc["legacy_note"])
        self.assertIn("~153", tc["legacy_note"])

    @unittest.skipIf(hasattr(os, "geteuid") and os.geteuid() == 0, "root reads a mode-000 file")
    def test_audit_names_the_unreadable_fasta(self):
        meta = self.legacy_meta(TARGET + CONT)
        meta_path = self.write_meta(meta)
        os.chmod(meta["fasta"], 0)
        try:
            a = self.audit("--fasta-meta", meta_path)
        finally:
            os.chmod(meta["fasta"], 0o644)
        f = self.finding(a, "contaminant_overlap")
        self.assertEqual(f["status"], "WARN")
        self.assertIn(f"the searched FASTA {meta['fasta']} could not be read", f["message"])
        self.assertFalse([x for x in a["findings"] if "could not read" in x["message"]
                          and meta_path in x["message"]])

    def test_old_sidecar_without_contaminants_is_not_flagged(self):
        meta = {"organism": "Homo sapiens", "taxid": 9606, "contaminant_set": "none",
                "n_contaminants_appended": 0}
        a = self.audit("--fasta-meta", self.write_meta(meta))
        self.assertFalse([x for x in a["findings"] if x["check"] == "contaminant_overlap"])

    def test_sample_quality_carries_the_same_old_sidecar_note(self):
        self.cont_only_matrix()
        sq = self.sample_quality("--fasta-meta", self.write_meta(self.legacy_meta(TARGET + CONT)))
        self.assertTrue(any("built before fetch_fasta.py checked" in f
                            and "ACTB (Cont_P60712)" in f for f in sq["flags"]), sq["flags"])
        self.assertIn("built before", sq["lost_to_contaminants"])


class KeepTargetContaminants(unittest.TestCase):
    """fetch --keep-target-contaminants rebuilds an OLD database faithfully, and provenance.py's
    reproduce.sh uses it for a sidecar written before the overlap check -- so a replay does not
    quietly produce a different database (153 fewer human Cont_ entries)."""

    def test_flag_disables_the_drop_but_still_records_the_pairs(self):
        with tempfile.TemporaryDirectory() as root:
            tgt, cp = os.path.join(root, "target.fasta"), os.path.join(root, "cont.fasta")
            with open(tgt, "w") as fh:
                fh.write(TARGET)
            with open(cp, "w") as fh:
                fh.write(CONT)
            out = os.path.join(root, "out", "search.fasta")
            argv = ["fetch_fasta.py", "fetch", "--path", tgt, "--contaminants", "universal",
                    "--contaminants-path", cp, "--keep-target-contaminants", "--out", out]
            with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(ff.main(), 0)
            with open(out + ".meta.json") as fh:
                m = json.load(fh)
            with open(out) as fh:
                fasta = fh.read()
        self.assertEqual(m["contaminant_target_rule"], "disabled (--keep-target-contaminants)")
        self.assertEqual(m["n_contaminants_dropped_as_target"], 0)
        self.assertEqual(m["n_contaminants_appended"], 6)
        self.assertIn(">sp|Cont_P60712|ACTB_BOVIN", headers(fasta))
        self.assertEqual({r["cont_acc"] for r in m["contaminants_identical_to_target_kept"]},
                         {"Cont_P60712", "Cont_Q3LI67"})
        self.assertTrue(any(w.startswith("--keep-target-contaminants: kept 2") for w in m["warnings"]))

    def repro(self, fasta_info):
        with tempfile.TemporaryDirectory() as tmp:
            env = os.path.join(tmp, "env.json")
            with open(env, "w") as fh:
                subprocess.run(["bash", os.path.join(SCRIPTS, "detect_env.sh")], stdout=fh,
                               stderr=subprocess.DEVNULL, check=True)
            wf = os.path.join(tmp, "wf")
            subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--acquisition", "DIA", "--instrument", "timsTOF HT",
                            "--organism-taxid", "9606", "--env", env, "--dest", wf],
                           capture_output=True, text=True, check=True)
            out = os.path.join(tmp, "repro")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                                "--outdir", out, "--workflow-manifest",
                                os.path.join(wf, "workflow.manifest.json"), "--engine", "diann",
                                "--acquisition", "DIA", "--instrument", "timsTOF HT",
                                "--organism-taxid", "9606", "--fasta-info", json.dumps(fasta_info)],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(out, "reproduce.sh")) as fh:
                return fh.read()

    BASE = {"proteome": "UP000005640", "organism": "Homo sapiens", "taxid": 9606,
            "content_used": "one_per_gene", "contaminant_set": "universal",
            "n_contaminants_appended": 381, "n_proteome": 20663}

    def fetch_line(self, sh):
        return next(ln for ln in sh.splitlines() if "--contaminants " in ln and "--out" in ln)

    def test_old_sidecar_replays_with_the_flag_and_says_why(self):
        sh = self.repro(self.BASE)                                   # no contaminant_target_rule
        self.assertIn("--keep-target-contaminants", self.fetch_line(sh))
        self.assertIn("built before target-identical contaminants were removed; this replays "
                      "that faithfully -- drop --keep-target-contaminants to get the corrected "
                      "database", sh)

    def test_new_sidecar_replays_without_the_flag(self):
        sh = self.repro({**self.BASE, "contaminant_target_rule": ff.CONTAMINANT_TARGET_RULE,
                         "n_contaminants_appended": 228})
        self.assertNotIn("--keep-target-contaminants", self.fetch_line(sh))
        self.assertNotIn("replays that faithfully", sh)

    def test_a_replay_of_a_replay_keeps_the_flag(self):
        sh = self.repro({**self.BASE, "contaminant_target_rule": ff.KEEP_TARGET_CONTAMINANTS_RULE})
        self.assertIn("--keep-target-contaminants", self.fetch_line(sh))

    def test_old_sidecar_without_contaminants_needs_no_flag(self):
        sh = self.repro({**self.BASE, "contaminant_set": "none", "n_contaminants_appended": 0})
        self.assertNotIn("--keep-target-contaminants", self.fetch_line(sh))


class EnzymeNames(unittest.TestCase):
    def test_slash_forms(self):
        self.assertEqual(ff.parse_enzymes("Trypsin/P"), ("trypsin",))       # the engine's name
        self.assertEqual(ff.parse_enzymes("Trypsin/Lys-C"), ("lysc", "trypsin"))
        self.assertEqual(ff.parse_enzymes("trypsin/P, Lys-C"), ("lysc", "trypsin"))
        self.assertEqual(ff.parse_enzymes("trypsin/pepsin"), ("pepsin", "trypsin"))  # not /P
        # A trailing /P folds for ANY enzyme, not just trypsin.
        for text, want in (("LysC/P", ("lysc",)), ("Lys-C/P", ("lysc",)),
                           ("Asp-N/P", ("aspn",)), ("Trypsin/P, Lys-C/P", ("lysc", "trypsin")),
                           ("trypsin / P", ("trypsin",))):
            self.assertEqual(ff.parse_enzymes(text), want, text)
        with self.assertRaises(argparse.ArgumentTypeError):
            ff.parse_enzymes("trypsin/papain")


if __name__ == "__main__":
    unittest.main()
