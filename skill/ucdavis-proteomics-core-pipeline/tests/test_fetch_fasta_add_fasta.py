#!/usr/bin/env python3
"""--add-fasta: the user's own target sequences (a bait, a tag) win over contaminant entries.

A staff search (2026-10-01, a TurboID proximity-labelling study with an EGFP control) needed
EGFP and a TurboID construct in the database. fetch_fasta had no way to add target sequences
to a --proteome build, and the Universal contaminant set holds wild-type GFP (Cont_P42212). EGFP
differs from it by three substitutions and one insertion, so the ordinary rule (an entry with
2+ peptides of its own stays a contaminant) kept Cont_P42212 -- and every peptide the two share
would have been kept out of the bait's quantification by --cont-quant-exclude Cont_ and
run_de.R. The staff built the database by hand in three steps and patched the sidecar.

Guards, all offline:
  * added entries are written after the proteome, before the contaminants, and recorded in the
    sidecar (path, sha256, each entry's accession / name / length / sequence sha256);
  * a contaminant entry with at least ADDED_SEQUENCE_SHARED_FRACTION of its peptides in an added
    sequence (most of it IS the added protein) is removed and recorded -- from the contaminant
    set and from a supplied database's own Cont_ entries -- while one sharing fewer stays, its
    shared peptides flagged by name as ambiguous (sidecar, Methods, report), so real
    contamination with it stays visible; the digestion enzymes in use stay;
  * a tagged, duplicated or clashing entry is refused; an entry identical to a proteome entry
    is warned;
  * Methods names the added sequences and the removed contaminants; reproduce.sh replays it.
"""
import contextlib
import hashlib
import io
import json
import os
import random
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fetch_fasta as ff  # noqa: E402

# UniProt P42212 (GFP_AEQVI), as the Universal contaminant set holds it.
WT_GFP = ("MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTL"
          "VTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLV"
          "NRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLAD"
          "HYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK")
# EGFP: wild type with Val inserted after Met1, F64L, S65T and H231L (wild-type numbering).
EGFP = WT_GFP[:1] + "V" + WT_GFP[1:63] + "LT" + WT_GFP[65:230] + "L" + WT_GFP[231:]
PIG_TRYPSIN = "FPTDDDDKIVGGYTCAANSIPYQVSLNSG"

PROTEOME = (">sp|P60709|ACTB_HUMAN Actin OS=Homo sapiens OX=9606 GN=ACTB\n"
            "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLK\n"
            ">sp|P04406|G3P_HUMAN GAPDH OS=Homo sapiens OX=9606 GN=GAPDH\n"
            "MGKVKVGVNGFGRIGRLVTRAAFNSGKVDIVAINDPFIDLNYMVYMFQYDSTHGKFHGTVK\n")
CONT = (f">sp|Cont_P42212|GFP_AEQVI Green fluorescent protein OS=Aequorea victoria OX=6100 GN=GFP\n"
        f"{WT_GFP}\n"
        f">sp|Cont_P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB\n"
        f"MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVK\n"
        f">sp|Cont_P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823\n{PIG_TRYPSIN}\n")
ADDED = f">IN|1000001|EGFP enhanced GFP bait\n{EGFP}\n>IN|1000002|TAG synthetic tag\nMDYKDDDDKGGSAWSHPQFEK\n"


class AddFasta(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.proteome = self.write("proteome.fasta", PROTEOME)
        self.cont = self.write("cont.fasta", CONT)
        self.added = self.write("bait.fasta", ADDED)
        self.out = os.path.join(self.root, "out", "search.fasta")

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, text):
        p = os.path.join(self.root, name)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def fetch(self, *args, path=None, contaminants=True):
        cont = ["--contaminants", "universal", "--contaminants-path", self.cont] if contaminants \
            else ["--contaminants", "none"]
        argv = ["fetch_fasta.py", "fetch", "--path", path or self.proteome, "--organism",
                "Homo sapiens", "--taxid", "9606", *cont, *args, "--out", self.out]
        err = io.StringIO()
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(err):
            try:
                rc = ff.main()
            except SystemExit as e:
                rc = e.code
        if rc not in (0, None):
            return rc, None, None, err.getvalue()
        with open(self.out + ".meta.json") as fh, open(self.out) as fa:
            return rc, json.load(fh), fa.read(), err.getvalue()

    def test_the_ordinary_rule_alone_keeps_wild_type_gfp_beside_an_egfp_bait(self):
        """Why the added-sequence rule exists: measured with the search's digest."""
        cont = ff._fasta_records(CONT)[:1]
        egfp = ff._fasta_records(f">IN|1000001|EGFP\n{EGFP}\n")
        self.assertEqual(ff.contaminants_matching_targets(cont, egfp), [])
        p = ff.peptide_overlap(cont, egfp)[0]
        self.assertEqual((p["n_peptides"], p["n_shared_peptides"], p["n_unique_peptides"]),
                         (28, 21, 3))

    def test_added_entries_are_written_recorded_and_win_over_wild_type_gfp(self):
        rc, m, fasta, err = self.fetch("--add-fasta", self.added)
        self.assertIn(rc, (0, None), err)
        # order: proteome, added, contaminants
        self.assertLess(fasta.index(">sp|P04406|"), fasta.index(">IN|1000001|EGFP"))
        self.assertLess(fasta.index(">IN|1000002|TAG"), fasta.index(">sp|Cont_P02769|"))
        self.assertNotIn("Cont_P42212", fasta)
        self.assertIn("Cont_P00761", fasta)            # no shared peptide: untouched
        self.assertEqual(m["n_added_sequences"], 2)
        self.assertEqual(m["n_proteome"], 2)
        self.assertEqual(m["n_sequences"], 2 + 2 + 2)
        self.assertEqual(m["n_entries"], 6)
        (f,) = m["added_sequences"]
        self.assertEqual(f["file"], os.path.abspath(self.added))
        with open(self.added, "rb") as fh:
            self.assertEqual(f["sha256"], hashlib.sha256(fh.read()).hexdigest())
        e = f["entries"][0]
        self.assertEqual((e["accession"], e["name"], e["length"]), ("1000001", "IN|1000001|EGFP", 239))
        self.assertEqual(e["sequence_sha256"], hashlib.sha256(EGFP.encode()).hexdigest())
        (d,) = m["contaminants_dropped_for_added_sequences"]
        self.assertEqual((d["cont_acc"], d["target_acc"], d["source"]),
                         ("Cont_P42212", "1000001", "contaminant_set"))
        self.assertEqual((d["n_shared_peptides"], d["n_unique_peptides"]), (21, 3))
        self.assertEqual(d["shared_fraction"], 0.75)
        self.assertEqual(len(d["shared_peptides"]), 21)
        self.assertEqual(m["added_sequence_shared_fraction"], ff.ADDED_SEQUENCE_SHARED_FRACTION)
        self.assertEqual(m["contaminants_sharing_peptides_with_added_sequences"], [])
        self.assertEqual(m["n_contaminants_dropped_for_added_sequences"], 1)
        self.assertEqual(m["added_sequence_rule"], ff.ADDED_SEQUENCE_RULE)
        # not a target protein of the organism: the auditors' list stays clean
        self.assertNotIn("Cont_P42212", {r["cont_acc"] for r in m["contaminants_dropped_as_target"]})
        self.assertIn("Cont_P42212 (GFP_AEQVI) ~ 1000001: 21 of its 28 peptides shared, 3 of its own",
                      m["contaminants_dropped_for_added_sequences_note"])
        self.assertIn(m["contaminants_dropped_for_added_sequences_note"], m["warnings"])

    def test_a_supplied_databases_own_cont_entry_leaves_too(self):
        db = self.write("with_cont.fasta", PROTEOME + CONT)
        rc, m, fasta, err = self.fetch("--add-fasta", self.added, path=db, contaminants=False)
        self.assertIn(rc, (0, None), err)
        self.assertNotIn("Cont_P42212", fasta)
        self.assertEqual([(r["cont_acc"], r["source"]) for r in
                          m["contaminants_dropped_for_added_sequences"]],
                         [("Cont_P42212", "supplied_database")])
        self.assertEqual(m["n_contaminants_already_present"], 2)

    def test_a_digestion_enzyme_in_use_stays(self):
        added = self.write("tryp.fasta", f">IN|1000003|TRYPLIKE bait\nMAAAK{PIG_TRYPSIN}\n")
        rc, m, fasta, err = self.fetch("--add-fasta", added)
        self.assertIn(rc, (0, None), err)
        self.assertIn("Cont_P00761", fasta)
        self.assertEqual([r["cont_acc"] for r in m["contaminants_kept_near_added_sequences"]],
                         ["Cont_P00761"])
        self.assertTrue(any("digestion enzyme used in this search" in w for w in m["warnings"]))

    def test_bad_entries_are_refused(self):
        cases = {
            "tagged.fasta": (">sp|Cont_X1|BAIT tagged\nMPEPTIDEKAAAR\n", "contaminant tag"),
            "twice.fasta": (">IN|7|A\nMPEPTIDEK\n>IN|7|B\nMSTVWYKR\n", "used twice"),
            "clash.fasta": (">sp|P60709|ACTB_MINE\nMPEPTIDEKAAAR\n", "already in the database"),
            "empty.fasta": ("no header here\n", "no FASTA entries"),
            "noseq.fasta": (">IN|8|EMPTY\n", "has no sequence"),
        }
        for name, (text, why) in cases.items():
            rc, m, _, _ = self.fetch("--add-fasta", self.write(name, text))
            self.assertIn(why, str(rc), name)
            self.assertIsNone(m, name)

    def test_an_entry_identical_to_a_proteome_entry_is_warned(self):
        actb = PROTEOME.split("\n")[1]
        rc, m, _, err = self.fetch("--add-fasta", self.write("same.fasta", f">IN|9|MYACTB\n{actb}\n"))
        self.assertIn(rc, (0, None), err)
        self.assertTrue(any("identical to database entry P60709" in w for w in m["warnings"]))

    def test_without_add_fasta_nothing_changes(self):
        rc, m, fasta, err = self.fetch()
        self.assertIn(rc, (0, None), err)
        self.assertIn("Cont_P42212", fasta)
        self.assertEqual((m["n_added_sequences"], m["added_sequences"], m["added_sequence_rule"]),
                         (0, [], None))

    def test_a_contaminant_sharing_one_peptide_of_twenty_is_kept_and_flagged(self):
        """Dropping it would also hide real contamination with it: its 19 own peptides are the
        evidence. So it stays, and the one peptide it shares is named as ambiguous."""
        rng = random.Random(7)
        peps = ["".join(rng.choice("ACDEFGHNQSTVWY") for _ in range(15)) + "K" for _ in range(20)]
        cont = self.write("cont20.fasta", f">sp|Cont_Q00020|SYN20_BOVIN Synthetic OS=Bos taurus OX=9913\n"
                                          f"{''.join(peps)}\n")
        bait = self.write("bait1.fasta", f">IN|1000009|BAIT bait\nMAAAK{peps[4]}GGGGGGGK\n")
        self.cont = cont
        rc, m, fasta, err = self.fetch("--add-fasta", bait)
        self.assertIn(rc, (0, None), err)
        self.assertIn("Cont_Q00020", fasta)
        self.assertEqual(m["contaminants_dropped_for_added_sequences"], [])
        (k,) = m["contaminants_sharing_peptides_with_added_sequences"]
        self.assertEqual((k["cont_acc"], k["target_acc"], k["source"]),
                         ("Cont_Q00020", "1000009", "contaminant_set"))
        self.assertEqual((k["n_peptides"], k["n_shared_peptides"], k["shared_fraction"]), (20, 1, 0.05))
        self.assertEqual(k["shared_peptides"], [peps[4]])
        note = m["contaminants_sharing_peptides_with_added_sequences_note"]
        self.assertIn(f"Cont_Q00020 (SYN20_BOVIN) ~ 1000009: {peps[4]} (1 of its 20 peptides)", note)
        self.assertIn("AMBIGUOUS", note)
        import make_methods as mm
        sent = mm.added_sequences_sentence(m)
        self.assertIn(f"1 contaminant entry sharing fewer of its peptides with it was kept, and the "
                      f"peptides shared are ambiguous between the contaminant and the added sequence: "
                      f"Cont_Q00020 (SYN20_BOVIN) and 1000009: {peps[4]}.", sent)
        import make_analysis_html as mah
        session = os.path.join(self.root, "session")
        os.makedirs(os.path.join(session, "input"))
        with open(os.path.join(session, "input", "search.fasta.meta.json"), "w") as fh:
            json.dump(m, fh)
        callout = mah.added_sequences_note(session)
        self.assertEqual(callout["kind"], "warning")
        self.assertIn(peps[4], callout["text"])

    def test_a_contaminant_split_between_the_proteome_and_the_bait_leaves(self):
        """2.10 review: the ordinary rule counts own peptides against EVERY target -- the proteome
        and the added sequences -- so one with 1 own peptide in the whole database goes."""
        rng = random.Random(7)
        pep = ["".join(rng.choice("ACDEFGHLMNQSTVWY") for _ in range(16)) + "K" for _ in range(7)]
        proteome = self.write("split_proteome.fasta",
                              f">sp|Q00001|HUM_HUMAN human protein GN=HUM\n{''.join(pep[0:3])}\n")
        bait = self.write("split_bait.fasta", f">IN|1000010|BAIT bait\n{''.join(pep[3:6])}\n")
        self.cont = self.write("split_cont.fasta",
                               f">sp|Cont_Q99999|CON_BOVIN bovine protein\n{''.join(pep)}\n")
        rc, m, fasta, err = self.fetch("--add-fasta", bait, path=proteome)
        self.assertIn(rc, (0, None), err)
        self.assertNotIn("Cont_Q99999", fasta)
        (d,) = m["contaminants_dropped_as_target"]
        self.assertEqual((d["cont_acc"], d["reason"], d["n_unique_peptides"]),
                         ("Cont_Q99999", "shared_peptides", 1))
        self.assertIn(d["target_acc"], ("Q00001", "1000010"))

    def test_a_fragpipe_tagged_supplied_contaminant_follows_the_rule(self):
        db = self.write("fp.fasta", PROTEOME + f">contam_sp|P42212|GFP_AEQVI GFP\n{WT_GFP}\n")
        rc, m, fasta, err = self.fetch("--add-fasta", self.added, path=db, contaminants=False)
        self.assertIn(rc, (0, None), err)
        self.assertNotIn("P42212", fasta)
        self.assertEqual([(r["cont_acc"], r["source"]) for r in
                          m["contaminants_dropped_for_added_sequences"]],
                         [("contam_P42212", "supplied_database")])

    def test_text_before_the_first_header_is_refused(self):
        rc, m, _f, _e = self.fetch("--add-fasta", self.write("pre.fasta", f"from a plasmid map\n{ADDED}"))
        self.assertIn("text before the first '>' header", str(rc))
        self.assertIsNone(m)

    def test_methods_and_reproduce_carry_the_added_sequences(self):
        _, m, _, _ = self.fetch("--add-fasta", self.added)
        import make_methods as mm
        s = mm.added_sequences_sentence(m)
        self.assertIn("2 user-supplied sequences (IN|1000001|EGFP, IN|1000002|TAG) were added to "
                      "the database as target entries (bait.fasta, SHA-256 "
                      + m["added_sequences"][0]["sha256"] + ").", s)
        self.assertIn("1 contaminant entry sharing most of its peptides with them (Cont_P42212 "
                      "(GFP_AEQVI)) was removed from the library, so the added sequences keep those "
                      "peptides.", s)
        self.assertEqual(mm.added_sequences_sentence({}), "")

    def test_reproduce_sh_adds_the_same_file(self):
        _, m, _, _ = self.fetch("--add-fasta", self.added)
        wm = self.write("workflow.manifest.json", "{}")
        out = os.path.join(self.root, "repro")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir", out,
                            "--workflow-manifest", wm, "--engine", "diann", "--acquisition", "DIA",
                            "--instrument", "timsTOF HT", "--organism-taxid", "9606",
                            "--fasta-info", json.dumps(m)], capture_output=True, text=True,
                           cwd=self.root)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out, "reproduce.sh")) as fh:
            sh = fh.read()
        line = next(ln for ln in sh.splitlines() if "--contaminants " in ln and "--out" in ln)
        # X8 F9: the added file is bundled, and reproduce.sh reads the copy, not the path of the
        # machine that ran -- the same bytes, as the sidecar's sha256 says
        bundled = os.path.join(out, "inputs", "added_sequences", "1_" + os.path.basename(self.added))
        self.assertIn("--add-fasta ./inputs/added_sequences/1_" + os.path.basename(self.added), line)
        with open(bundled, "rb") as a, open(self.added, "rb") as b:
            self.assertEqual(a.read(), b.read())
        self.assertIn(f"Target sequences were added from {os.path.abspath(self.added)} (2 entries, "
                      f"sha256 {m['added_sequences'][0]['sha256'][:12]}...)", sh)


# --------------------------------------------------------------------------------------------
# Independent 2.10 review (wrong-results risks): fetch_fasta.py --add-fasta and --path organism.
#
#
#   * A supplied database made by FragPipe tags its contaminants contam_ (FRAGPIPE_CONT_TAG), and
#     the --add-fasta rule looks only at Cont_ entries: wild-type GFP stays beside an EGFP bait,
#     neither removed nor flagged, and nothing in the sidecar or Methods says so.
#   * Text before the first '>' of an --add-fasta file is written verbatim after the proteome, so it
#     becomes part of the LAST proteome entry's sequence in the searched database.
#   * The --add-fasta rule and the ordinary peptide rule (MIN_UNIQUE_PEPTIDES) are judged apart: a
#     contaminant whose peptides are split between a proteome protein and the bait keeps ONE peptide
#     of its own in the whole database and still stays, taking peptides from both.
#   * --path --organism / --taxid that contradict the curated organism table are recorded as given
#     (FRAN files the search under a taxid that is not the organism Methods names).
#
# Offline; no network.
# --------------------------------------------------------------------------------------------
import argparse  # noqa: E402



ACTB = (">sp|P60709|ACTB_HUMAN Actin OS=Homo sapiens OX=9606 GN=ACTB\n"
        "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPR\n")


class _Fetch(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.out = os.path.join(self.root, "out", "search.fasta")

    def tearDown(self):
        self._td.cleanup()

    def write(self, name, text):
        p = os.path.join(self.root, name)
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def fetch(self, proteome_text, added_text):
        argv = ["fetch_fasta.py", "fetch", "--path", self.write("p.fasta", proteome_text),
                "--organism", "Homo sapiens", "--taxid", "9606", "--contaminants", "none",
                "--add-fasta", self.write("a.fasta", added_text), "--out", self.out]
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            try:
                rc = ff.main()
            except SystemExit as e:
                rc = e.code
        if rc not in (0, None):
            return rc, None, None
        with open(self.out + ".meta.json") as fh, open(self.out) as fa:
            return rc, json.load(fh), fa.read()


class FragPipeTaggedSuppliedDatabase(_Fetch):
    def test_a_contam_tagged_wild_type_gfp_is_removed_or_flagged_beside_an_egfp_bait(self):
        # The same database with the Core's tag works (Cont_P42212 removed); FragPipe's tag is
        # ignored by the --add-fasta rule -- every other contaminant rule here reads both tags.
        db = ACTB + f">contam_sp|P42212|GFP_AEQVI Green fluorescent protein\n{WT_GFP}\n"
        rc, meta, fasta = self.fetch(db, f">sp|EGFP01|EGFP_SYN EGFP bait\n{EGFP}\n")
        self.assertIn(rc, (0, None))
        seen = ([r["cont_acc"] for r in meta["contaminants_dropped_for_added_sequences"]]
                + [r["cont_acc"] for r in
                   meta["contaminants_sharing_peptides_with_added_sequences"]])
        self.assertTrue(any("P42212" in acc for acc in seen),
                        "contam_sp|P42212 (wild-type GFP, 75% of its peptides in EGFP) was "
                        "neither removed nor flagged; the sidecar and Methods say nothing")


class TextBeforeTheFirstHeader(_Fetch):
    def test_a_preamble_line_never_joins_the_last_proteome_entry(self):
        rc, meta, fasta = self.fetch(
            ACTB, "EGFP construct from plasmid pX330\n>EGFP_construct\nMVSKGEELFTGVVPILVELDGDVNGHK\n")
        if rc not in (0, None):
            return                                    # refusing the file is a correct fix too
        recs = ff._fasta_records(fasta)
        actb = next(lines for h, lines in recs if "P60709" in h)
        self.assertEqual(ff._record_seq(actb), "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPR",
                         "the --add-fasta file's first line was appended to ACTB's sequence in "
                         "the searched database")


class RulesJudgedApart(unittest.TestCase):
    def test_a_contaminant_with_one_own_peptide_in_the_whole_database_leaves(self):
        random.seed(7)
        aa = "ACDEFGHLMNQSTVWY"                       # no K/R/P/I: one tryptic peptide each
        pep = ["".join(random.choice(aa) for _ in range(16)) + "K" for _ in range(7)]
        proteome = f">sp|Q00001|HUM_HUMAN human protein GN=HUM\n{''.join(pep[0:3])}\n"
        added = f">sp|BAIT01|BAIT_SYN bait\n{''.join(pep[3:6])}\n"
        contam = f">sp|Cont_Q99999|CON_BOVIN bovine protein\n{''.join(pep)}\n"
        used = ff.parse_enzymes(ff.DEFAULT_ENZYMES)
        # the order cmd_fetch applies them: the ordinary rule vs EVERY target -- the proteome and
        # the added sequences (the fix: judged together) -- then --add-fasta
        kept, _d, _e = ff.drop_target_contaminants(
            contam, proteome + ff.added_as_targets(ff._fasta_records(added)), used,
            ff.MIN_UNIQUE_PEPTIDES)
        kept, _d2, _flagged, _e2 = ff.drop_contaminants_near_added(
            kept, ff._fasta_records(added), used, "contaminant_set", only_tagged=False)
        own = ff.peptide_overlap(ff._fasta_records(contam),
                                 ff._fasta_records(proteome + added))[0]["n_unique_peptides"]
        self.assertEqual(own, 1)
        self.assertNotIn("Cont_Q99999", kept,
                         f"Cont_Q99999 has {own} peptide of its own among proteome + added "
                         f"sequences (MIN_UNIQUE_PEPTIDES = {ff.MIN_UNIQUE_PEPTIDES}) but stays, "
                         f"taking 3 peptides from Q00001 and 3 from the bait")


class PathOrganismContradiction(unittest.TestCase):
    def test_a_name_and_taxid_the_curated_table_contradicts_are_refused(self):
        a = argparse.Namespace(organism="Homo sapiens", taxid=10090, ncbi_organism="",
                               ncbi_taxid=0, path="/x.fasta", ncbi_accession="")
        with contextlib.redirect_stderr(io.StringIO()):
            try:
                got = ff.user_organism(a)
            except SystemExit:
                return
        self.fail(f"--organism 'Homo sapiens' --taxid 10090 (Mus musculus in ORGANISM_TAXIDS) "
                  f"was recorded as {got}")


if __name__ == "__main__":
    unittest.main()
