#!/usr/bin/env python3
"""A keratin sample (hair, wool, feather, skin, nail) keeps its keratin: it is the analyte.

msalemi's hair benchmark (SET28, skill 2.8.0): the Universal contaminant set carries ~190
keratin-family Cont_ entries. The target-identity rule removes the ones identical to human
proteins, but ~46 survive for human -- KRT34, KRTAPs, mouse hair and sheep wool keratins -- and
for a hair sample those are the analyte's paralogs. DIA-NN's --cont-quant-exclude Cont_ keeps
every peptide they share out of quantification, and run_de.R (2.8.0) removed every precursor
whose Protein.Ids names any Cont_ entry. --keratin-sample reached only the auditors, at step 8c,
after the search.

Guards:
  * the one keratin rule (fetch_fasta.is_keratin_gene) covers the sheep wool keratins that have
    no gene name, and not Krtcap2;
  * fetch_fasta.py --keratin-sample removes every keratin-family Cont_ entry -- from the appended
    set and from a supplied database's own -- and records keratin_sample + the entries; without
    the flag the database is exactly what it was;
  * run_search.py --keratin-sample refuses a database that still holds them, before anything is
    written; a keratin-built sidecar makes a search one, recorded in search_provenance.json;
  * run_de.R keeps a keratin sample's keratin precursors and still removes trypsin and BSA, on
    both methods, from every source that can say so; a non-keratin run is unchanged; when it
    cannot tell, it behaves as before and says so in its log and provenance;
  * methods text, the report and reproduce.sh say what was kept or dropped;
  * SKILL.md asks the sample type at step 3 and every documented --keratin-sample exists.
"""
import contextlib
import csv
import hashlib
import io
import json
import os
import random
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
REPO_CONTAMINANTS = os.path.join(os.path.dirname(os.path.dirname(SKILL)), "contaminants")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import fetch_fasta as ff        # noqa: E402
import make_methods as mm       # noqa: E402
import make_analysis_html as mah  # noqa: E402
import run_search as rs         # noqa: E402


def _organism(args):
    """A --path database needs its organism named (fetch_fasta.user_organism); these tests are
    about the contaminant rules, so they say "no single organism" unless they name one."""
    named = {"--organism", "--taxid", "--ncbi-organism", "--ncbi-taxid"}
    return ["--organism", "none"] if "--path" in args and not named & set(args) else []


def protein(seed, n=160):
    """A deterministic synthetic protein, tryptic sites every few residues."""
    r = random.Random(seed)
    aa = "ACDEFGHILMNPQSTVWY"
    return "M" + "".join(r.choice(aa) + ("K" if i % 9 == 8 else "") for i in range(n))


GAPDH, ACTB, KRT31 = protein("gapdh"), protein("actb"), protein("krt31")
TARGET = (f">sp|P04406|G3P_HUMAN GAPDH OS=Homo sapiens OX=9606 GN=GAPDH PE=1 SV=3\n{GAPDH}\n"
          f">sp|P60709|ACTB_HUMAN Actin OS=Homo sapiens OX=9606 GN=ACTB PE=1 SV=1\n{ACTB}\n"
          f">sp|Q15323|K1H1_HUMAN Keratin, type I cuticular Ha1 OS=Homo sapiens OX=9606 "
          f"GN=KRT31 PE=1 SV=3\n{KRT31}\n")
# Shaped like the Universal set's headers. Two are identical to a target (the identity rule takes
# them either way); three are keratin-family and target-unrelated (only --keratin-sample takes
# them): mouse Krt31, a sheep wool keratin with NO gene name, a human KRTAP. Trypsin, BSA and
# Krtcap2 (keratinocyte-associated, not a keratin) always stay Cont_.
CONT = (f">sp|Cont_P60712|ACTB_BOVIN Actin OS=Bos taurus OX=9913 GN=ACTB PE=1 SV=1\n{ACTB}\n"
        f">sp|Cont_Q15323|K1H1_HUMAN Keratin, type I cuticular Ha1 OS=Homo sapiens OX=9606 "
        f"GN=KRT31 PE=1 SV=3\n{KRT31}\n"
        f">sp|Cont_Q61765|K1H1_MOUSE Keratin, type I cuticular Ha1 OS=Mus musculus OX=10090 "
        f"GN=Krt31 PE=1 SV=2\n{protein('mkrt31')}\n"
        f">sp|Cont_P02438|KRB2A_SHEEP Keratin, high-sulfur matrix protein, B2A OS=Ovis aries "
        f"OX=9940 PE=1 SV=2\n{protein('krb2a', 90)}\n"
        f">sp|Cont_Q9BYR0|KRA47_HUMAN Keratin-associated protein 4-7 OS=Homo sapiens OX=9606 "
        f"GN=KRTAP4-7 PE=1 SV=2\n{protein('krtap47', 80)}\n"
        f">sp|Cont_P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1\n{protein('tryp')}\n"
        f">sp|Cont_P02769|ALBU_BOVIN Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4\n"
        f"{protein('bsa')}\n"
        f">sp|Cont_Q5RL79|KTAP2_MOUSE Keratinocyte-associated protein 2 OS=Mus musculus "
        f"OX=10090 GN=Krtcap2 PE=1 SV=2\n{protein('ktap2', 60)}\n")
KERATIN_CONT = {"Cont_Q61765", "Cont_P02438", "Cont_Q9BYR0"}
OTHER_CONT = {"Cont_P00761", "Cont_P02769", "Cont_Q5RL79"}


def headers(text):
    return [ln.split()[0][1:] for ln in text.splitlines() if ln.startswith(">")]


def accs(text):
    return {h.split("|")[1] if "|" in h else h for h in headers(text)}


def sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


class KeratinRule(unittest.TestCase):
    """fetch_fasta.is_keratin_gene is the one keratin rule: gene or protein name."""

    def test_genes_and_names(self):
        for g in ("KRT1", "Krt31", "KRTAP5-9", "KRT87P", "KRT222", "krt6a", "KRTHA1", "KRTHB6"):
            self.assertTrue(ff.is_keratin_gene(g), g)
        for g in ("KRTCAP2", "Krtcap2", "KRTDAP", "ALB", "", None):
            self.assertFalse(ff.is_keratin_gene(g), g)
        self.assertTrue(ff.is_keratin_gene("", "Keratin, high-sulfur matrix protein, B2A"))
        self.assertTrue(ff.is_keratin_gene("", "Putative keratin-87 protein"))
        self.assertTrue(ff.is_keratin_gene(None, "Keratin-associated protein 6-1"))
        self.assertFalse(ff.is_keratin_gene("Krtcap2", "Keratinocyte-associated protein 2"))

    def test_protein_name_is_read_from_the_header(self):
        self.assertEqual(ff._protein_name(">sp|Cont_P02438|KRB2A_SHEEP Keratin, high-sulfur "
                                          "matrix protein, B2A OS=Ovis aries OX=9940 PE=1 SV=2"),
                         "Keratin, high-sulfur matrix protein, B2A")
        self.assertEqual(ff._protein_name(">XP_012345.1 keratin 31 [Ovis aries]"),
                         "keratin 31 [Ovis aries]")
        self.assertEqual(ff._protein_name(">P1"), "")

    @unittest.skipUnless(os.path.isfile(os.path.join(REPO_CONTAMINANTS,
                                                     "Universal_Contaminants.fasta")),
                         "needs the repo's contaminants/ folder")
    def test_the_real_universal_set(self):
        """189 keratin-family entries, the 14 gene-less sheep wool keratins among them; every one
        says keratin in its gene or name, and no non-keratin does."""
        with open(os.path.join(REPO_CONTAMINANTS, "Universal_Contaminants.fasta")) as fh:
            recs = ff._fasta_records(fh.read())
        hits = ff.keratin_contaminants(recs)
        self.assertEqual(len(hits), 189)
        sheep = [r for _, r in hits if r["cont_entry"].endswith("_SHEEP") and not r["cont_gene"]]
        self.assertEqual(len(sheep), 14)
        hit_idx = {i for i, _ in hits}
        for i, (h, _l) in enumerate(recs):
            said = bool(re.search(r"\bkeratin\b|GN=KRT", h, re.I))
            self.assertEqual(i in hit_idx, said, h.strip())

    @unittest.skipUnless(os.path.isfile(os.path.join(REPO_CONTAMINANTS,
                                                     "Mouse_Tissue_Contaminants.fasta")),
                         "needs the repo's contaminants/ folder")
    def test_krtcap2_is_not_a_keratin(self):
        with open(os.path.join(REPO_CONTAMINANTS, "Mouse_Tissue_Contaminants.fasta")) as fh:
            hits = ff.keratin_contaminants(ff._fasta_records(fh.read()))
        self.assertNotIn("Cont_Q5RL79", {r["cont_acc"] for _, r in hits})


class FetchKeratinSample(unittest.TestCase):
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

    def fetch(self, *args, name="search.fasta"):
        out = os.path.join(self.root, "out", name)
        argv = ["fetch_fasta.py", "fetch", *args, *_organism(args), "--out", out]
        err = io.StringIO()
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(err):
            self.assertEqual(ff.main(), 0, err.getvalue())
        with open(out + ".meta.json") as fh:
            meta = json.load(fh)
        with open(out) as fh:
            return meta, fh.read(), err.getvalue(), out

    def test_keratin_sample_removes_every_keratin_family_entry(self):
        m, fasta, err, _ = self.fetch("--path", self.target, "--contaminants", "universal",
                                      "--contaminants-path", self.cont, "--keratin-sample")
        self.assertIs(m["keratin_sample"], True)
        # the identity rule still takes (and records) the target-identical ones ...
        self.assertEqual({r["cont_acc"] for r in m["contaminants_dropped_as_target"]},
                         {"Cont_P60712", "Cont_Q15323"})
        # ... and --keratin-sample the rest of the keratins, gene-less sheep one included
        dropped = {r["cont_acc"]: r for r in m["contaminants_dropped_keratin_sample"]}
        self.assertEqual(set(dropped), KERATIN_CONT)
        self.assertEqual({r["source"] for r in dropped.values()}, {"contaminant_set"})
        self.assertEqual(dropped["Cont_P02438"]["protein_name"],
                         "Keratin, high-sulfur matrix protein, B2A")
        self.assertEqual(m["n_contaminants_dropped_keratin_sample"], 3)
        self.assertEqual(m["keratin_contaminants_in_database"], [])
        self.assertEqual(m["keratin_family_rule"], ff.KERATIN_FAMILY_RULE)
        # counts stay truthful: 8 in the set - 2 target - 3 keratin
        self.assertEqual(m["n_contaminants_in_set"], 8)
        self.assertEqual(m["n_contaminants_appended"], 3)
        self.assertEqual(m["n_sequences"], 3 + 3)
        self.assertEqual(m["n_entries"], 3 + 3)
        # the database: targets intact, no keratin Cont_, the other contaminants still Cont_
        self.assertTrue(fasta.startswith(TARGET))
        cont_left = {a for a in accs(fasta) if a.startswith("Cont_")}
        self.assertEqual(cont_left, OTHER_CONT)
        self.assertEqual(ff.keratin_contaminants_in_fasta(os.path.join(self.root, "out",
                                                                       "search.fasta")), [])
        # said, as a record -- not a build warning to resolve before publication
        self.assertIn("keratin sample (--keratin-sample): removed 3 keratin-family", err)
        self.assertIn("besides 1 keratin entry already removed", m["keratin_sample_note"])
        self.assertNotIn(m["keratin_sample_note"], m["warnings"])
        self.assertEqual(m["diann_cont_quant_exclude"], "Cont_")

    def test_without_the_flag_the_database_is_as_before(self):
        """Non-keratin runs: the keratins stay Cont_ exactly as in 2.8.0; the sidecar only gains
        the record of them."""
        m, fasta, err, _ = self.fetch("--path", self.target, "--contaminants", "universal",
                                      "--contaminants-path", self.cont)
        expected = TARGET + "".join(
            h + "".join(ls) for h, ls in ff._fasta_records(CONT)
            if not h.startswith((">sp|Cont_P60712|", ">sp|Cont_Q15323|")))
        self.assertEqual(fasta, expected)
        self.assertIs(m["keratin_sample"], False)
        # nobody answered step 3: the "not keratin" is a default, and recorded as one (rule 2)
        self.assertEqual(m["keratin_sample_source"], "default")
        self.assertEqual(m["contaminants_dropped_keratin_sample"], [])
        self.assertEqual(m["n_contaminants_dropped_keratin_sample"], 0)
        self.assertIsNone(m["keratin_sample_note"])
        self.assertEqual(set(m["keratin_contaminants_in_database"]), KERATIN_CONT)
        self.assertEqual(m["contaminant_tags_in_database"], ["Cont_"])
        self.assertEqual(m["n_contaminants_appended"], 6)
        self.assertNotIn("keratin sample", err)
        # the user's "no" gives the same database, recorded as an answer
        m2, fasta2, _, _ = self.fetch("--path", self.target, "--contaminants", "universal",
                                      "--contaminants-path", self.cont, "--no-keratin-sample",
                                      name="no.fasta")
        self.assertEqual(fasta2, fasta)
        self.assertEqual((m2["keratin_sample"], m2["keratin_sample_source"]), (False, "user"))

    def test_the_searched_species_own_keratins_are_listed(self):
        """keratin_contaminants_same_species: the entries whose OX= is the build's taxid -- only
        they may stand alone as the sample's keratin (run_de.R)."""
        mrs = os.path.join(self.root, "MRS")
        os.makedirs(mrs)
        with open(os.path.join(mrs, "UP000005640_9606.fasta"), "w") as fh:
            fh.write(TARGET)
        data = {"taxonomy": {"scientificName": "Homo sapiens", "taxonId": 9606},
                "proteomeType": "Reference proteome", "geneCount": 3, "superkingdom": "eukaryota"}
        with mock.patch.object(ff, "HIVE_MRS", mrs), \
                mock.patch.object(ff, "_get_json", return_value=(data, {})):
            m, _, _, _ = self.fetch("--proteome", "UP000005640", "--hive", "--contaminants",
                                    "universal", "--contaminants-path", self.cont,
                                    "--no-keratin-sample")
        self.assertEqual(m["taxid"], 9606)
        self.assertEqual(set(m["keratin_contaminants_in_database"]), KERATIN_CONT)
        self.assertEqual(m["keratin_contaminants_same_species"], ["Cont_Q9BYR0"])

    def test_a_supplied_database_loses_its_own_keratin_entries(self):
        """The Core's staged human + contaminant FASTA, --path ... --keratin-sample: its keratin
        Cont_ entries go; everything else is written back byte for byte."""
        filler = "".join(f">sp|Cont_Z{i:05d}|FILL_X filler\n{protein(f'f{i}', 30)}\n"
                         for i in range(20))
        combined = os.path.join(self.root, "combined.fasta")
        with open(combined, "w") as fh:
            fh.write(TARGET + CONT + filler)
        m, fasta, _, _ = self.fetch("--path", combined, "--contaminants", "universal",
                                    "--keratin-sample")
        self.assertEqual(m["n_contaminants_appended"], 0)          # used as-is otherwise ...
        dropped = {r["cont_acc"]: r["source"] for r in m["contaminants_dropped_keratin_sample"]}
        self.assertEqual(set(dropped), KERATIN_CONT | {"Cont_Q15323"})
        self.assertEqual(set(dropped.values()), {"supplied_database"})
        expected = "".join(h + "".join(ls) for h, ls in ff._fasta_records(TARGET + CONT + filler)
                           if h.split("|")[1] not in KERATIN_CONT | {"Cont_Q15323"})
        self.assertEqual(fasta.rstrip("\n"), expected.rstrip("\n"))
        self.assertEqual(m["n_contaminants_already_present"], 8 + 20 - 4)
        self.assertEqual(m["keratin_contaminants_in_database"], [])
        # ... and without the flag it is untouched
        m2, fasta2, _, _ = self.fetch("--path", combined, "--contaminants", "universal",
                                      name="plain.fasta")
        self.assertEqual(fasta2.rstrip("\n"), (TARGET + CONT + filler).rstrip("\n"))
        self.assertEqual(set(m2["keratin_contaminants_in_database"]),
                         KERATIN_CONT | {"Cont_Q15323"})

    def test_the_auditors_list_follows_the_sidecar(self):
        """target_contaminants: a keratin-built sidecar makes keratin the analyte without the
        auditors' own --keratin-sample -- a gene-less wool keratin included (keratin-review F1:
        sheep wool has no gene on either side, and a gene-only test listed 17 of them as possible
        contamination of a wool sample)."""
        sheep = _write(os.path.join(self.root, "sheep.fasta"),
                       f">sp|P04406|G3P_SHEEP GAPDH OS=Ovis aries OX=9940 GN=GAPDH PE=1 SV=1\n"
                       f"{GAPDH}\n"
                       f">sp|P02438|KRB2A_SHEEP Keratin, high-sulfur matrix protein, B2A OS=Ovis "
                       f"aries OX=9940 PE=1 SV=2\n{protein('krb2a', 90)}\n")
        m, _, err, _ = self.fetch("--path", sheep, "--contaminants", "universal",
                                  "--contaminants-path", self.cont, "--keratin-sample")
        wool = next(r for r in m["contaminants_dropped_as_target"] if r["cont_acc"] == "Cont_P02438")
        self.assertEqual((wool["gene"], wool["cont_gene"], wool["keratin_family"]), ("", "", True))
        self.assertIn("besides 1 keratin entry already removed", m["keratin_sample_note"])
        self.assertEqual([r["cont_acc"] for r in ff.target_contaminants(m)["dropped"]], [])
        m["keratin_sample"] = False
        self.assertEqual({r["cont_acc"] for r in ff.target_contaminants(m)["dropped"]},
                         {"Cont_P02438"})
        self.assertEqual(ff.target_contaminants(m, True)["dropped"], [])
        # an older sidecar (no stamp) falls back to the genes
        old = {k: v for k, v in wool.items() if k != "keratin_family"}
        self.assertFalse(ff.record_is_keratin(old))
        self.assertTrue(ff.record_is_keratin(dict(old, cont_gene="KRTAP1-4")))
        # the human build: KRT31 identical to its target is keratin by gene too
        m, _, _, _ = self.fetch("--path", self.target, "--contaminants", "universal",
                                "--contaminants-path", self.cont, "--keratin-sample",
                                name="human.fasta")
        self.assertEqual({r["gene"] for r in ff.target_contaminants(m)["dropped"]}, {"ACTB"})
        m["keratin_sample"] = False
        self.assertEqual({r["gene"] for r in ff.target_contaminants(m)["dropped"]},
                         {"ACTB", "KRT31"})
        self.assertEqual({r["gene"] for r in ff.target_contaminants(m, True)["dropped"]}, {"ACTB"})
        self.assertIsNone(ff.keratin_sample_recorded({"organism": "x"}))
        self.assertIsNone(ff.keratin_sample_recorded({"keratin_sample": "yes"}))

    def test_keratin_database_check_every_source(self):
        m, _, _, out = self.fetch("--path", self.target, "--contaminants", "universal",
                                  "--contaminants-path", self.cont)
        meta = out + ".meta.json"
        # 1. the sidecar's own list
        r = ff.keratin_database_check(meta)
        self.assertTrue(r["checked"])
        self.assertEqual(set(r["accessions"]), KERATIN_CONT)
        self.assertEqual(r["source"], meta)
        # 2. an older sidecar (no list): its FASTA re-read, while the sha256 matches
        old = {k: v for k, v in m.items() if not k.startswith(("keratin", "n_contaminants_dropped_"
                                                                          "keratin",
                                                                "contaminants_dropped_keratin"))}
        old["taxid"] = 9606
        old_meta = os.path.join(self.root, "old.meta.json")
        with open(old_meta, "w") as fh:
            json.dump(old, fh)
        r = ff.keratin_database_check(old_meta)
        self.assertTrue(r["checked"], r["why"])
        self.assertEqual(set(r["accessions"]), KERATIN_CONT)
        self.assertEqual(r["own_species"], ["Cont_Q9BYR0"])       # the human KRTAP, by OX=
        self.assertEqual(r["source"], out)
        # ... not once the FASTA has changed
        with open(out, "a") as fh:
            fh.write(">sp|Cont_X|Y changed\nPEPTIDEK\n")
        r = ff.keratin_database_check(old_meta)
        self.assertFalse(r["checked"])
        self.assertIn("sha256 no longer matches", r["why"])
        # 3. no sidecar: the FASTA the search names
        sdir = os.path.join(self.root, "search_out")
        os.makedirs(sdir)
        both = _write(os.path.join(self.root, "both.fasta"), TARGET + CONT)
        with open(os.path.join(sdir, "search_provenance.json"), "w") as fh:
            json.dump({"fasta": both}, fh)
        r = ff.keratin_database_check(None, sdir)
        self.assertTrue(r["checked"], r["why"])
        self.assertEqual(set(r["accessions"]), KERATIN_CONT | {"Cont_Q15323"})
        # no sidecar: the organism is the one the FASTA's own targets carry (OX=9606)
        self.assertEqual((r["taxid"], set(r["own_species"])),
                         (9606, {"Cont_Q9BYR0", "Cont_Q15323"}))
        # 4. nothing to read: said, never guessed
        empty = os.path.join(self.root, "empty")
        os.makedirs(empty)
        r = ff.keratin_database_check(None, empty)
        self.assertFalse(r["checked"])
        self.assertIn("names no FASTA", r["why"])
        # the CLI run_de.R calls prints the same answer
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fetch_fasta.py"), "keratin-db",
                            "--fasta-meta", meta], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(set(json.loads(p.stdout)["accessions"]), KERATIN_CONT)


class RunSearchRefusesKeratinContaminants(unittest.TestCase):
    """run_search.py --keratin-sample on a database that still holds keratin-family Cont_
    entries: REFUSED before anything is written or submitted."""

    def setUp(self):
        from test_single_shot_sbatch import _Workspace
        self._td = tempfile.TemporaryDirectory()
        self.w = _Workspace(self._td.name)
        d = self._td.name
        for name, extra in (("plain.fasta", []), ("keratin.fasta", ["--keratin-sample"])):
            argv = ["fetch_fasta.py", "fetch", "--path", _write(os.path.join(d, "t.fasta"), TARGET),
                    "--contaminants", "universal", "--contaminants-path",
                    _write(os.path.join(d, "c.fasta"), CONT), *extra, "--organism", "none",
                    "--out", os.path.join(d, name)]
            with mock.patch.object(sys, "argv", argv), \
                    contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(ff.main(), 0)
        self.plain, self.keratin = os.path.join(d, "plain.fasta"), os.path.join(d, "keratin.fasta")

    def tearDown(self):
        self._td.cleanup()

    def search(self, fasta, *extra):
        self.w.fasta = fasta
        return self.w.generate(*extra, check=False)

    def prov(self):
        with open(os.path.join(self.w.out, "search_provenance.json")) as fh:
            return json.load(fh)

    def test_refused_before_anything_is_written(self):
        p = self.search(self.plain, "--keratin-sample")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("REFUSED -- nothing was converted, written or submitted", p.stderr)
        self.assertIn("still holds 3 keratin-family contaminant entries", p.stderr)
        self.assertIn("Cont_P02438 (KRB2A_SHEEP)", p.stderr)
        # the suggested rebuild names the organism the sidecar recorded (--path needs one)
        self.assertIn(f"--path {os.path.join(self._td.name, 't.fasta')} --organism none", p.stderr)
        self.assertIn("--keratin-sample --out", p.stderr)
        self.assertFalse(os.path.exists(self.w.out), "the search folder was created")
        self.assertFalse(any(n.startswith("diann_job") for n in os.listdir(self._td.name)))

    def test_a_keratin_database_is_searched_and_recorded(self):
        p = self.search(self.keratin, "--keratin-sample")
        self.assertEqual(p.returncode, 0, p.stderr)
        rec = self.prov()
        self.assertIs(rec["keratin_sample"], True)
        self.assertEqual(rec["keratin_sample_source"], "user")
        self.assertEqual(rec["keratin_sample_check"]["source"], "--keratin-sample")
        self.assertEqual(rec["keratin_sample_check"]["keratin_contaminants_in_database"], 0)

    def test_the_sidecar_makes_it_a_keratin_sample_without_the_flag(self):
        p = self.search(self.keratin)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("says the database was built with --keratin-sample", p.stderr)
        rec = self.prov()
        self.assertIs(rec["keratin_sample"], True)
        self.assertEqual(rec["keratin_sample_check"]["source"], self.keratin + ".meta.json")

    def test_a_non_keratin_search_is_unchanged(self):
        p = self.search(self.plain)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertNotIn("keratin", p.stderr.lower())
        rec = self.prov()
        self.assertIs(rec["keratin_sample"], False)
        # plain.fasta was built with neither flag: the sidecar's "default" is carried, not upgraded
        self.assertEqual(rec["keratin_sample_source"], "default")
        self.assertIsNone(rec["keratin_sample_check"]["keratin_contaminants_in_database"])

    def test_fragpipe_on_a_contam_database_is_refused(self):
        """keratin-review F2: FragPipe's DIA route drops Philosopher's contam_ tag (library.tsv), so
        those contaminants would reach the DE as plain accessions. A FragPipe search on such a
        database -- as --fasta, or as the workflow's database.db-path -- is refused."""
        d = self._td.name
        fp_db = _write(os.path.join(d, "fragpipe.fas"), FRAGPIPE_DB)
        p = self.search(fp_db, "--engine", "fragpipe")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("REFUSED -- nothing was converted, written or submitted", p.stderr)
        self.assertIn("3 contaminant entries tagged contam_", p.stderr)
        self.assertIn("fetch_fasta.py fetch", p.stderr)
        self.assertFalse(os.path.exists(self.w.out))
        wf = _write(os.path.join(d, "fp.workflow"),
                    f"database.decoy-tag=rev_\ndatabase.db-path={fp_db}\nworkflow.threads=4\n")
        with self.assertRaises(SystemExit) as cm:
            rs.refuse_fragpipe_contam_database(wf, self.keratin)
        self.assertIn(fp_db, str(cm.exception))
        # the skill's own database (Cont_ only) is searched
        clean = _write(os.path.join(d, "clean.workflow"), f"database.db-path={self.plain}\n")
        self.assertIsNone(rs.refuse_fragpipe_contam_database(clean, self.plain))

    def test_a_keratin_fasta_that_cannot_be_read_is_refused(self):
        with self.assertRaises(SystemExit) as cm:
            with contextlib.redirect_stderr(io.StringIO()):
                rs.keratin_sample_check(os.path.join(self._td.name, "missing.fasta"), True)
        self.assertIn("could not be read", str(cm.exception))


def _write(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    return path


def r_has(*pkgs):
    if not shutil.which("Rscript"):
        return False
    expr = ("quit(status = if (all(vapply(c(%s), requireNamespace, logical(1), quietly = TRUE))) "
            "0 else 1)" % ", ".join(f'"{p}"' for p in pkgs))
    return subprocess.run(["Rscript", "-e", expr], capture_output=True).returncode == 0


def rscript(expr):
    p = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError(p.stderr[-1500:])
    return p.stdout


# A hair-sample-shaped DIA-NN report (human): 80 sample proteins x 4 precursors over 6 runs (A x3,
# B x3), the hair keratin KRT31 (P31000) with one precursor it shares with MOUSE Krt31
# (Cont_Q61765), a SHEEP wool keratin group seen only as Cont_P02438, a HUMAN KRTAP group seen only
# as Cont_Q9BYR0, trypsin (Cont_P00761) and one sample precursor shared with it.
SYNTH_R = r'''
set.seed(7)
runs <- sprintf("run%02d", 1:6); grp <- rep(c("A", "B"), each = 3)
rows <- list()
add <- function(pg, ids, gene, prec, base, effB = 0) {
  for (r in seq_along(runs)) {
    lv <- base + rnorm(1, 0, 0.2) + (grp[r] == "B") * effB
    int <- 2^(lv + rnorm(1, 0, 0.3))
    if (runif(1) < plogis(-(lv - 13) * 2)) next
    rows[[length(rows) + 1]] <<- data.frame(Run = runs[r], Precursor.Id = prec,
      Protein.Group = pg, Protein.Ids = ids, Protein.Names = paste0(gene, "_X"), Genes = gene,
      Proteotypic = 1L, Precursor.Normalised = int, Precursor.Quantity = int,
      Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001,
      Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
  }
}
for (i in 1:80) {
  pg <- sprintf("P%05d", i); base <- runif(1, 11, 20)
  for (k in 1:4) add(pg, pg, sprintf("G%d", i), sprintf("PEP%dK%d2", i, k), base + rnorm(1, 0, 1),
                     effB = if (i <= 5) 1.5 else 0)
}
for (k in 1:3) add("P31000", "P31000", "KRT31", sprintf("KRTPEP%dK2", k), 19, effB = 1)
add("P31000", "P31000;Cont_Q61765", "KRT31", "KRTSHAREDK2", 19, effB = 1)
for (k in 1:3) add("Cont_P02438", "Cont_P02438", "", sprintf("WOOLPEP%dK2", k), 18, effB = 1)
for (k in 1:3) add("Cont_Q9BYR0", "Cont_Q9BYR0", "KRTAP4-7", sprintf("KRTAPPEP%dK2", k), 18, effB = 1)
for (k in 1:3) add("Cont_P00761", "Cont_P00761", "", sprintf("TRYPPEP%dK2", k), 18)
add("P00001", "P00001;Cont_P00761", "G1", "TRYPSHAREDK2", 17)
d <- do.call(rbind, rows)
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum)
names(pgm)[3] <- "PG.MaxLFQ"
d <- merge(d, pgm, by = c("Protein.Group", "Run"))
arrow::write_parquet(d, file.path(OUT, "report.parquet"))
write.csv(data.frame(File.Name = runs, Group = grp), file.path(OUT, "conditions.csv"),
          row.names = FALSE)
'''
# the report's keratin-family contaminant entries; only the human KRTAP is the sample species' own
REPORT_KERATIN = {"Cont_Q61765", "Cont_P02438", "Cont_Q9BYR0"}
OWN_KERATIN = {"Cont_Q9BYR0"}
# the searched database: human targets (OX=9606 -- the organism the FASTA's own targets carry)
# plus the contaminant entries the report names
DB_FASTA = (TARGET +
            f">sp|Cont_Q61765|K1H1_MOUSE Keratin, type I cuticular Ha1 OS=Mus musculus OX=10090 "
            f"GN=Krt31 PE=1 SV=2\n{protein('mkrt31')}\n"
            f">sp|Cont_P02438|KRB2A_SHEEP Keratin, high-sulfur matrix protein, B2A OS=Ovis aries "
            f"OX=9940 PE=1 SV=2\n{protein('krb2a', 90)}\n"
            f">sp|Cont_Q9BYR0|KRA47_HUMAN Keratin-associated protein 4-7 OS=Homo sapiens OX=9606 "
            f"GN=KRTAP4-7 PE=1 SV=2\n{protein('krtap47', 80)}\n"
            f">sp|Cont_P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1\n{protein('tryp')}\n")
# the head of a DIA-NN 2.7.0 log (HIVE, 2026-09-29): the command line is its sixth line. No
# --fasta: fetch_fasta.py reads the searched FASTA from it too, which these fixtures leave unnamed.
DIANN_LOG = ("\nDIA-NN 2.7.0 Academia  (Data-Independent Acquisition by Neural Networks)\n"
             "Compiled on Sep 15 2026 15:30:01\nCurrent date and time: Tue Sep 29 20:21:57 2026\n"
             "Logical CPU cores: 224\n/quobyte/proteomics-grp/dia-nn/diann-linux --f /d/r01.raw "
             "--out {out} --threads 16 --qvalue 0.01{flags} \n\n"
             "Existing .quant files will be used\n")
SIDECAR = {"organism": "Homo sapiens", "taxid": 9606, "n_contaminants_appended": 3,
           "contaminant_set": "universal", "contaminant_target_rule": ff.CONTAMINANT_TARGET_RULE,
           "min_unique_peptides": 2, "contaminants_dropped_as_target": [],
           "contaminants_identical_to_target_kept": [], "diann_cont_quant_exclude": "Cont_"}


@unittest.skipUnless(r_has("limpa", "limma", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limpa/limma/arrow/dplyr/tidyr/jsonlite")
class RunDeKeepsKeratinForAKeratinSample(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = t = tempfile.mkdtemp()
        rscript(f'OUT <- "{t}"\n' + SYNTH_R)
        db = _write(os.path.join(t, "db.fasta"), DB_FASTA)
        # the DIA-NN log that wrote the report, its command line as DIA-NN echoes it: what the
        # dpc/maxlfq descriptors say DIA-NN did with contaminant peptides is read from it
        # (make_methods.diann_cont_quant_exclude). sp_search / fp_search copy only the report.
        _write(os.path.join(t, "report.log.txt"), DIANN_LOG.format(
            out=os.path.join(t, "report.parquet"), flags=" --cont-quant-exclude Cont_"))
        listed = dict(SIDECAR, fasta=db, keratin_sample=False,
                      keratin_contaminants_in_database=sorted(REPORT_KERATIN),
                      keratin_contaminants_same_species=sorted(OWN_KERATIN),
                      contaminant_tags_in_database=["Cont_"])
        sidecars = {
            # skill 2.9.0, the user said "not keratin" (step 3): it lists the keratins it holds
            "notker.json": dict(listed, keratin_sample_source="user"),
            # ... built with neither flag: "not keratin" is an unconfirmed default
            "default.json": dict(listed, keratin_sample_source="default"),
            # built with --keratin-sample
            "built.json": dict(SIDECAR, fasta=db, keratin_sample=True,
                               keratin_sample_source="user",
                               n_contaminants_dropped_keratin_sample=40,
                               keratin_contaminants_in_database=[]),
            # skill 2.8.0: no keratin fields; its FASTA is still there and unchanged
            "old.json": dict(SIDECAR, fasta=db, sha256=sha256(db)),
        }
        for name, meta in sidecars.items():
            _write(os.path.join(t, name), json.dumps(meta))
        # a search folder whose search_provenance.json says keratin sample, no sidecar beside it
        sp_dir = os.path.join(t, "sp_search")
        os.makedirs(sp_dir)
        shutil.copy(os.path.join(t, "report.parquet"), sp_dir)
        sp_rec = {"fasta": db, "keratin_sample": True, "keratin_sample_source": "user",
                  "keratin_sample_check": {"keratin_contaminants_in_database": 3}}
        _write(os.path.join(sp_dir, "search_provenance.json"), json.dumps(sp_rec))
        # FragPipe's DIA route: the report in <workdir>/dia-quant-output/, the provenance above it
        fp_dir = os.path.join(t, "fp_search")
        os.makedirs(os.path.join(fp_dir, "dia-quant-output"))
        shutil.copy(os.path.join(t, "report.parquet"), os.path.join(fp_dir, "dia-quant-output"))
        _write(os.path.join(fp_dir, "search_provenance.json"), json.dumps(sp_rec))
        rep = os.path.join(t, "report.parquet")
        sp_rep = os.path.join(sp_dir, "report.parquet")
        fp_rep = os.path.join(fp_dir, "dia-quant-output", "report.parquet")
        fm = lambda n: ["--fasta-meta", os.path.join(t, n)]      # noqa: E731
        cls.runs = {}
        for key, report, args in (
                ("plain", rep, ["--method", "dpc"]),
                ("notker", rep, ["--method", "dpc", *fm("notker.json")]),
                ("default", rep, ["--method", "dpc", *fm("default.json")]),
                ("legacy", rep, ["--method", "dpc", *fm("old.json")]),
                ("built", rep, ["--method", "dpc", *fm("built.json")]),
                ("flag", rep, ["--method", "dpc", "--keratin-sample", *fm("notker.json")]),
                ("flag_ml", rep, ["--method", "maxlfq", "--keratin-sample", *fm("notker.json")]),
                ("plain_ml", rep, ["--method", "maxlfq"]),
                ("recheck", rep, ["--method", "dpc", "--keratin-sample", *fm("old.json")]),
                ("unknown", rep, ["--method", "dpc", "--keratin-sample"]),
                ("sp", sp_rep, ["--method", "dpc"]),
                ("fp_dia", fp_rep, ["--method", "dpc"])):
            out = os.path.join(t, key)
            p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"), "--input", report,
                                "--metadata", os.path.join(t, "conditions.csv"),
                                "--outdir", out, *args],
                               capture_output=True, text=True, cwd=t)
            cls.runs[key] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-2500:])
        return out

    def cont(self, key):
        with open(os.path.join(self.out(key), "de_provenance.json")) as fh:
            return json.load(fh)["contaminants"]

    def methods(self, key):
        with open(os.path.join(self.out(key), "methods.txt")) as fh:
            return fh.read()

    def matrix(self, key, name="Expression_Matrix.csv"):
        with open(os.path.join(self.out(key), name), newline="") as fh:
            return {r["Protein.Group"]: r for r in csv.DictReader(fh)}

    def precursors(self, key):
        rds = os.path.join(self.out(key), "DE-LIMP_session.rds")
        return set(rscript(f's <- readRDS("{rds}"); cat(rownames(s$raw_data$E), sep = "\\n")')
                   .split())

    def test_a_keratin_sample_keeps_its_keratin_and_removes_trypsin(self):
        for key in ("flag", "recheck", "sp", "fp_dia"):
            with self.subTest(run=key):
                em = self.matrix(key)
                self.assertIn("Cont_Q9BYR0", em)          # the sample species' own KRTAP stays
                self.assertIn("P31000", em)
                self.assertNotIn("Cont_P00761", em)       # trypsin does not
                pre = self.precursors(key)
                # mouse Krt31: kept, because the precursor also names the sample's own KRT31
                self.assertIn("KRTSHAREDK2", pre)
                self.assertNotIn("TRYPSHAREDK2", pre)     # shared with trypsin: removed
                k = self.cont(key)["keratin_sample"]
                self.assertIs(k["value"], True)
                self.assertEqual(k["database"], "keratins_in_database")
                self.assertEqual(set(k["exempt_accessions"]), REPORT_KERATIN)
                self.assertEqual(set(k["exempt_alone_accessions"]), OWN_KERATIN)
                self.assertIs(k["caution"], True)
                self.assertNotIn("unchecked_tags", k)       # every tag in the report was checked

    def test_another_species_keratin_alone_is_not_the_samples(self):
        """keratin-review F4: a sheep wool peptide in a human hair sample, naming no human protein,
        is not the sample's keratin -- it stays a contaminant."""
        for key in ("flag", "flag_ml"):
            with self.subTest(run=key):
                self.assertNotIn("Cont_P02438", self.matrix(key))
        self.assertFalse(any(p.startswith("WOOLPEP") for p in self.precursors("flag")))

    def test_the_counts_add_up(self):
        plain, flag = self.cont("plain"), self.cont("flag")
        kept = flag["keratin_sample"]["n_precursors_kept"]
        self.assertEqual(kept, 4)                           # KRTSHAREDK2 + 3 KRTAP precursors
        self.assertEqual(plain["n_precursors"] - flag["n_precursors"], kept)
        self.assertEqual(flag["n_protein_groups"], 2)       # trypsin + sheep wool
        self.assertEqual(plain["n_protein_groups"], 3)
        with open(os.path.join(self.out("flag"), "contaminants_removed.csv"), newline="") as fh:
            rem = {r["Protein.Group"] for r in csv.DictReader(fh)}
        self.assertIn("Cont_P02438", rem)
        self.assertNotIn("Cont_Q9BYR0", rem)
        self.assertNotIn("P31000", rem)

    def test_maxlfq_keeps_them_too_and_says_what_that_means(self):
        """keratin-review F3: on maxlfq the quantities are DIA-NN's PG.MaxLFQ, computed without
        the shared peptides; the record says so from the pipeline's own descriptor."""
        em, plain = self.matrix("flag_ml"), self.matrix("plain_ml")
        self.assertIn("Cont_Q9BYR0", em)
        self.assertNotIn("Cont_P00761", em)
        self.assertNotIn("Cont_Q9BYR0", plain)
        k = self.cont("flag_ml")["keratin_sample"]
        self.assertEqual((k["n_precursors_kept"], k["requantified"]), (4, False))
        self.assertIs(self.cont("flag")["keratin_sample"]["requantified"], True)
        m = " ".join(self.methods("flag_ml").split())
        self.assertIn("so their protein groups stay in the analysis;", m)
        # DIA-NN's own sentence, from the pipeline's descriptor (contaminants.R
        # DIANN_KEPT_AS_REPORTED) -- other engines say what they did (test_contaminants_every_engine)
        self.assertIn("But their quantities are the search engine's own (not re-derived from "
                      "precursors here), so what it left out (DIA-NN --cont-quant-exclude", m)
        self.assertIn("keratin is still under-quantified", m)
        self.assertIs(self.cont("flag_ml")["keratin_sample"]["under_quantified"], True)
        self.assertIn("--cont-quant-exclude Cont_, per the DIA-NN command line in report.log.txt",
                      m)
        self.assertNotIn("and quantified as sample proteins", m)
        self.assertIn("and quantified as sample proteins", self.methods("flag"))

    def test_diann_is_said_to_exclude_only_when_its_log_says_so(self):
        """sage-review N1: the dpc descriptor said DIA-NN's --cont-quant-exclude had left the kept
        keratin peptides out of its normalisation for every DIA-NN report. It is read from the
        command line in the DIA-NN log beside the report; FragPipe's DIA route here (fp_dia) and
        a search folder with no log (sp) are NOT RECORDED, never assumed."""
        k = self.cont("flag")["keratin_sample"]
        self.assertIn("DIA-NN --cont-quant-exclude Cont_, per the DIA-NN command line in "
                      "report.log.txt", k["kept_quant"])
        self.assertIs(k["under_quantified"], False)
        for key in ("fp_dia", "sp"):
            with self.subTest(run=key):
                k = self.cont(key)["keratin_sample"]
                self.assertIn("[not recorded -- confirm]", k["kept_quant"])
                self.assertNotIn("DIA-NN --cont-quant-exclude", k["kept_quant"])
                self.assertIs(k["under_quantified"], False)       # dpc re-derives from precursors

    def test_where_it_was_read_from(self):
        self.assertEqual(self.cont("flag")["keratin_sample"]["source"], "--keratin-sample")
        self.assertEqual(self.cont("recheck")["keratin_sample"]["checked_in"],
                         os.path.join(self.tmp, "db.fasta"))
        for key, folder in (("sp", "sp_search"), ("fp_dia", "fp_search")):
            with self.subTest(run=key):
                k = self.cont(key)["keratin_sample"]
                self.assertEqual(os.path.realpath(k["source"]),
                                 os.path.realpath(os.path.join(self.tmp, folder,
                                                               "search_provenance.json")))
                self.assertEqual(k["sample_source"], "user")
        b = self.cont("built")["keratin_sample"]
        self.assertEqual((b["value"], b["database"], b["n_removed_at_build"], b["exempt_accessions"]),
                         (True, "keratins_removed_at_build", 40, []))

    def test_methods_and_log_say_what_was_kept(self):
        m = self.methods("flag")
        self.assertIn("Keratin sample: yes (--keratin-sample) -- keratin is the analyte", m)
        self.assertIn("4 precursors mapping only to keratin-family contaminant entries (3 such", m)
        self.assertIn("were KEPT", m)
        self.assertIn("CAUTION: the search database was not built with fetch_fasta.py "
                      "--keratin-sample", m)
        self.assertIn("keratin sample:", self.runs["flag"][0].stderr)       # the R warning
        b = self.methods("built")
        self.assertIn("The FASTA was built with fetch_fasta.py --keratin-sample, which removed "
                      "its 40", b)
        self.assertNotIn("CAUTION", b)

    def test_cannot_tell_which_entries_behaves_as_before_and_says_so(self):
        self.assertEqual(set(self.matrix("unknown")), set(self.matrix("plain")))
        k = self.cont("unknown")["keratin_sample"]
        self.assertEqual((k["value"], k["database"]), (True, "unknown"))
        self.assertIn("could not be determined", k["note"])
        self.assertIn("REMOVED as contaminants", k["note"])
        self.assertIn("CAUTION: a keratin sample, but which contaminant (Cont_/contam_) entries",
                      self.methods("unknown"))
        self.assertNotIn("unchecked_tags", k)       # nothing was checked: said once, above
        self.assertIn("keratin sample:", self.runs["unknown"][0].stderr)

    def test_not_recorded_behaves_as_before_and_says_so(self):
        k = self.cont("plain")["keratin_sample"]
        self.assertIsNone(k["value"])
        self.assertIn("is not recorded", k["note"])
        self.assertIn("[run_de] keratin sample: whether these samples are keratin",
                      self.runs["plain"][0].stderr)
        # no database to learn the keratin entries from: nothing more to say in the methods
        self.assertNotIn("Keratin sample", self.methods("plain"))

    def test_an_unconfirmed_not_keratin_is_tagged_when_keratin_was_removed(self):
        """keratin-review F5/F6 (rule 2): "not keratin" as a default -- or not recorded at all --
        removed keratin precursors; ONE tagged methods line says so."""
        line = ("Keratin sample: not confirmed -- 7 precursors mapping to keratin-family "
                "contaminant entries were removed as contaminants; whether samples are "
                "keratinous was not recorded (DEFAULT — not user-confirmed).")
        for key in ("default", "legacy"):
            with self.subTest(run=key):
                self.assertIn(line, self.methods(key))
                self.assertEqual(self.methods(key).count("Keratin sample"), 1)
                k = self.cont(key)["keratin_sample"]
                self.assertEqual((k["default_removed"], k["n_keratin_precursors_removed"]),
                                 (True, 7))
                # the DE itself is the ordinary run's
                with open(os.path.join(self.out("plain"), "Expression_Matrix.csv")) as a, \
                        open(os.path.join(self.out(key), "Expression_Matrix.csv")) as b:
                    self.assertEqual(a.read(), b.read())
        self.assertEqual(self.cont("default")["keratin_sample"]["sample_source"], "default")
        self.assertIsNone(self.cont("legacy")["keratin_sample"]["value"])

    def test_a_non_keratin_run_is_unchanged(self):
        """The user said "not keratin": the same matrix, DE table, contaminant counts and methods
        as a run that records nothing -- only the keratin record differs."""
        for name in ("Expression_Matrix.csv", "DE_dpc_B.A.csv", "contaminants_removed.csv"):
            with self.subTest(table=name):
                with open(os.path.join(self.out("plain"), name)) as a, \
                        open(os.path.join(self.out("notker"), name)) as b:
                    self.assertEqual(a.read(), b.read())
        a, b = self.cont("plain"), self.cont("notker")
        for rec in (a, b):
            for k in ("keratin_sample", "fasta_meta", "database_checked", "database_risk",
                      "database_note"):
                rec.pop(k, None)
        self.assertEqual(a, b)
        k = self.cont("notker")["keratin_sample"]
        self.assertEqual((k["value"], k["sample_source"], k["exempt_accessions"]),
                         (False, "user", []))
        self.assertNotIn("note", k)
        self.assertNotIn("default_removed", k)
        self.assertNotIn("Keratin sample", self.methods("notker"))
        self.assertNotIn("keratin", self.runs["notker"][0].stderr.lower())

    def test_reproducibility_script_keeps_the_same_rows(self):
        for key, meth in (("flag", "dpc"), ("flag_ml", "maxlfq")):
            with self.subTest(method=meth):
                with open(os.path.join(self.out(key), "reproducibility_log.R")) as fh:
                    src = fh.read()
                self.assertIn("keratin_kept  <- c('Cont_", src)
                self.assertIn("keratin_alone <- c('Cont_Q9BYR0')", src)
                rerun = os.path.join(self.tmp, f"rerun_{key}")
                os.makedirs(rerun, exist_ok=True)
                shutil.copy(os.path.join(self.out(key), "reproducibility_log.R"), rerun)
                p = subprocess.run(["Rscript", "reproducibility_log.R"], cwd=rerun,
                                   capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, p.stderr[-1500:])
                de = f"DE_{meth}_B.A.csv"
                orig = read_lfc(os.path.join(self.out(key), de))
                again = read_lfc(os.path.join(rerun, "de_results_rerun", de))
                self.assertEqual(set(orig), set(again))
                self.assertIn("Cont_Q9BYR0", again)
                self.assertNotIn("Cont_P02438", again)
                self.assertEqual({k for k, v in orig.items() if v is None},
                                 {k for k, v in again.items() if v is None})
                self.assertLess(max(abs(orig[k] - again[k]) for k in orig
                                    if orig[k] is not None), 1e-9)


def read_lfc(path):
    with open(path, newline="") as fh:
        return {r["Protein.Group"]: (None if r["logFC"] in ("", "NA") else float(r["logFC"]))
                for r in csv.DictReader(fh)}


def de_prov(**keratin):
    c = dict(policy="removed", removed=True, tag="Cont_", id_column="Protein.Ids",
             n_precursors=12, n_protein_groups=3, n_sample_groups_sharing=1,
             n_sample_groups_all_shared=0)
    if keratin:
        c["keratin_sample"] = keratin
    return {"display_label": "DPC-Quant + limma (limpa)", "q_columns": ["Q.Value"],
            "q_cutoffs": [0.01], "design": "~ 0 + groups", "de_engine": "limpa::dpcDE",
            "adjp": 0.05, "logfc": 1, "logfc_role": "reference_line_only", "contaminants": c}


# `contam_` marks the contaminants of a database FragPipe/Philosopher builds with --contam
# (>contam_sp|P00761|TRYP_PIG, beside rev_ decoys). FragPipe adds none when it searches, and its
# DIA route drops the tag (library.tsv: bare accessions in report.tsv -- keratin-review, verified
# on HIVE), so the one report that carries it is the DDA adapter's: one row per protein x run,
# Protein.Group = contam_P00761 (run_search.py adapt_fragpipe + protein_ids.py). contam_P00761 is
# trypsin, contam_O43790 the human hair keratin KRT86, contam_P02438 a sheep wool keratin.
FRAG_R = r'''
set.seed(11)
runs <- sprintf("run%02d", 1:6); grp <- rep(c("A", "B"), each = 3)
rows <- list()
add <- function(pg, gene, base, effB = 0) {
  for (r in seq_along(runs))
    rows[[length(rows) + 1]] <<- data.frame(Run = runs[r], Protein.Group = pg,
      PG.MaxLFQ = 2^(base + rnorm(1, 0, 0.25) + (grp[r] == "B") * effB),
      Q.Value = 0, Lib.Q.Value = 0, Lib.PG.Q.Value = 0, Genes = gene,
      Protein.Names = paste0(gene, "_HUMAN"), stringsAsFactors = FALSE)
}
for (i in 1:60) add(sprintf("P%05d", i), sprintf("G%d", i), runif(1, 14, 22),
                    effB = if (i <= 5) 1.5 else 0)
add("contam_P00761", "", 18)
add("contam_O43790", "KRT86", 19, effB = 1)
add("contam_P02438", "", 18, effB = 1)
arrow::write_parquet(do.call(rbind, rows), file.path(OUT, "adapted.parquet"))
write.csv(data.frame(File.Name = runs, Group = grp), file.path(OUT, "conditions.csv"),
          row.names = FALSE)
'''
FRAGPIPE_DB = (f">sp|P04406|G3P_HUMAN GAPDH OS=Homo sapiens OX=9606 GN=GAPDH PE=1 SV=3\n{GAPDH}\n"
               f">contam_sp|P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1\n"
               f"{protein('tryp')}\n"
               f">contam_sp|O43790|KRT86_HUMAN Keratin, type II cuticular Hb6 OS=Homo sapiens "
               f"OX=9606 GN=KRT86 PE=1 SV=1\n{protein('krt86')}\n"
               f">contam_sp|P02438|KRB2A_SHEEP Keratin, high-sulfur matrix protein, B2A OS=Ovis "
               f"aries OX=9940 PE=1 SV=2\n{protein('krb2a', 90)}\n"
               f">rev_sp|P04406|G3P_HUMAN decoy\n{GAPDH[::-1]}\n"
               f">rev_contam_sp|O43790|KRT86_HUMAN decoy\n{protein('krt86')[::-1]}\n")

class FragPipeContamTag(unittest.TestCase):
    """contam_ is a contaminant tag exactly like Cont_: one set, in fetch_fasta.py and its R
    mirror; the keratin rule reads both."""

    def test_r_mirrors_the_tag_set(self):
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not available")
        out = rscript(f'source("{os.path.join(SCRIPTS, "contaminants.R")}"); '
                      f'cat(CONTAMINANT_TAGS, sep = "\\n")')
        self.assertEqual(tuple(out.split()), ff.CONTAMINANT_TAGS)
        self.assertEqual(ff.CONTAMINANT_TAGS, ("Cont_", "contam_"))

    def test_r_rule_reads_both_tags(self):
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not available")
        ids = ["P1", "Cont_P2", "contam_P3", "P4;contam_P5", "xcontam_P6", "P7;Xcontam_P8",
               "contam_O43790", "contam_O43790;Cont_Q61765", "contam_O43790;contam_P00761",
               "P9;Cont_Q61765", "NA"]
        vec = ", ".join(repr(i) for i in ids)
        both = 'c("contam_O43790", "Cont_Q61765")'
        out = rscript(f'source("{os.path.join(SCRIPTS, "contaminants.R")}"); x <- c({vec}); '
                      f'x[x == "NA"] <- NA; cat(is_contaminant(x), "\\n"); '
                      f'cat(is_contaminant(x, exempt = list(shared = {both}, alone = {both})), "\\n"); '
                      f'cat(is_contaminant(x, exempt = list(shared = {both}, alone = "contam_O43790")))'
                      .replace("'", '"'))
        plain, exempt, species = [ln.split() for ln in out.strip().splitlines()]
        T, F = "TRUE", "FALSE"
        self.assertEqual(plain, [F, T, T, T, F, F, T, T, T, T, F])
        # a keratin sample's exemption treats both tags alike: exempt only when EVERY tagged
        # accession is an exempt keratin
        self.assertEqual(exempt, [F, T, T, T, F, F, F, F, T, F, F])
        # ... and another species' keratin (Cont_Q61765, not in `alone`) only beside a sample
        # protein (P9): alone, or beside another tagged entry, it stays a contaminant
        self.assertEqual(species, [F, T, T, T, F, F, F, T, T, F, F])

    def test_keratin_rule_reads_fragpipe_headers(self):
        recs = ff._fasta_records(FRAGPIPE_DB)
        got = {r["cont_acc"]: r["tag"] for _, r in ff.keratin_contaminants(recs)}
        self.assertEqual(got, {"contam_O43790": "contam_", "contam_P02438": "contam_"})
        self.assertEqual(ff.contaminant_tags_in(recs), ["contam_"])
        self.assertEqual(ff.contaminant_tags_in(ff._fasta_records(CONT)), ["Cont_"])

    def fetch(self, root, *args):
        out = os.path.join(root, "out", "search.fasta")
        argv = ["fetch_fasta.py", "fetch", *args, *_organism(args), "--out", out]
        with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(ff.main(), 0)
        with open(out + ".meta.json") as fh, open(out) as fa:
            return json.load(fh), fa.read(), out

    def test_fetch_drops_fragpipe_keratins_from_a_supplied_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = _write(os.path.join(tmp, "fragpipe.fas"), FRAGPIPE_DB)
            m, fasta, _ = self.fetch(tmp, "--path", db, "--contaminants", "none",
                                     "--keratin-sample")
            self.assertEqual({(r["cont_acc"], r["source"]) for r in
                              m["contaminants_dropped_keratin_sample"]},
                             {("contam_O43790", "supplied_database"),
                              ("contam_P02438", "supplied_database")})
            self.assertIn(">contam_sp|P00761|TRYP_PIG", fasta)       # trypsin stays
            self.assertNotIn(">contam_sp|O43790|", fasta)
            self.assertIn(">rev_contam_sp|O43790|", fasta)          # a decoy is not an entry
            self.assertEqual(m["keratin_contaminants_in_database"], [])
            self.assertEqual(m["contaminant_tags_in_database"], ["contam_"])
            m2, fasta2, _ = self.fetch(tmp, "--path", db, "--contaminants", "none")
            self.assertEqual(fasta2.rstrip("\n"), FRAGPIPE_DB.rstrip("\n"))
            self.assertEqual(set(m2["keratin_contaminants_in_database"]),
                             {"contam_O43790", "contam_P02438"})

    def test_run_search_refuses_fragpipe_keratins_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            db = _write(os.path.join(tmp, "fragpipe.fas"), FRAGPIPE_DB)
            with self.assertRaises(SystemExit) as cm:
                with contextlib.redirect_stderr(io.StringIO()):
                    rs.keratin_sample_check(db, True)
        self.assertIn("contam_O43790 (KRT86_HUMAN)", str(cm.exception))
        self.assertIn("REFUSED", str(cm.exception))

    def test_report_flags_and_sentences(self):
        self.assertEqual(mah.background_flag("", "contam_P00761"), "contaminant")
        kept = frozenset({"contam_O43790"})
        self.assertIsNone(mah.background_flag("KRT86", "contam_O43790", None, kept))
        s = mm.de_paragraph(de_prov(value=True, database="keratins_removed_at_build",
                                    n_removed_at_build=46, exempt_accessions=[],
                                    unchecked_tags=["contam_"], n_precursors_unchecked=7,
                                    caution=True, note="7 precursors map to contam_ ..."))
        self.assertIn("7 precursors of contam_-tagged contaminant entries -- contaminants the "
                      "skill did not add, never checked for keratins -- were removed", s)
        self.assertNotIn("added itself", s)
        n = mah.keratin_note(de_prov(value=True, database="keratins_removed_at_build",
                                     exempt_accessions=[], unchecked_tags=["contam_"],
                                     caution=True, note="7 precursors map to contam_ ..."))
        self.assertEqual(n["kind"], "warning")


@unittest.skipUnless(r_has("limpa", "limma", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limpa/limma/arrow/dplyr/tidyr/jsonlite")
class RunDeFiltersFragPipeContam(unittest.TestCase):
    """The DDA adapter's contam_ rows (the one FragPipe report that keeps the tag) leave the DE
    like Cont_ ones; a keratin sample keeps its own species' contam_ keratins when the database
    check lists them, and says so when the database checked never held any contam_ entry."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = t = tempfile.mkdtemp()
        rscript(f'OUT <- "{t}"\n' + FRAG_R)
        # the Philosopher --contam database, checked by fetch_fasta.py: it lists both keratins,
        # one of them the searched organism's own (human KRT86), and holds contam_ entries only
        _write(os.path.join(t, "fp.json"), json.dumps(dict(
            SIDECAR, keratin_sample=False, keratin_sample_source="user",
            keratin_contaminants_in_database=["contam_O43790", "contam_P02438"],
            keratin_contaminants_same_species=["contam_O43790"],
            contaminant_tags_in_database=["contam_"])))
        # a keratin-sample build by fetch_fasta.py (Cont_ only) -- not the database these contam_
        # rows came from, so they were never checked
        _write(os.path.join(t, "built.json"), json.dumps(dict(
            SIDECAR, keratin_sample=True, keratin_sample_source="user",
            n_contaminants_dropped_keratin_sample=40, keratin_contaminants_in_database=[],
            contaminant_tags_in_database=["Cont_"])))
        ad = os.path.join(t, "adapted.parquet")
        fm = lambda n: ["--fasta-meta", os.path.join(t, n)]      # noqa: E731
        cls.runs = {}
        for key, args in (("plain", []),
                          ("ker", ["--keratin-sample", *fm("fp.json")]),
                          ("unchecked", fm("built.json"))):
            out = os.path.join(t, key)
            p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"), "--input", ad,
                                "--metadata", os.path.join(t, "conditions.csv"),
                                "--method", "maxlfq", "--outdir", out, *args],
                               capture_output=True, text=True, cwd=t)
            cls.runs[key] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-2500:])
        return out

    def cont(self, key):
        with open(os.path.join(self.out(key), "de_provenance.json")) as fh:
            return json.load(fh)["contaminants"]

    def matrix(self, key):
        with open(os.path.join(self.out(key), "Expression_Matrix.csv"), newline="") as fh:
            return {r["Protein.Group"] for r in csv.DictReader(fh)}

    def methods(self, key):
        with open(os.path.join(self.out(key), "methods.txt")) as fh:
            return fh.read()

    def test_contam_rows_are_removed_like_cont_rows(self):
        em = self.matrix("plain")
        for g in ("contam_P00761", "contam_O43790", "contam_P02438"):
            self.assertNotIn(g, em)
        self.assertIn("P00002", em)
        c = self.cont("plain")
        self.assertEqual(c["policy"], "removed")
        self.assertEqual(c["id_column"], "Protein.Group")
        self.assertEqual(c["tag"], "contam_")
        self.assertEqual(c["n_protein_groups"], 3)
        self.assertEqual(c["tags"], ["Cont_", "contam_"])
        self.assertIn("'Cont_' or 'contam_'", c["rule"])
        self.assertIn("mapping to a contam_ entry", self.methods("plain"))

    def test_a_keratin_sample_keeps_its_own_species_contam_keratin(self):
        em = self.matrix("ker")
        self.assertIn("contam_O43790", em)        # human KRT86 in a human sample
        self.assertNotIn("contam_P02438", em)     # sheep wool, naming no sample protein
        self.assertNotIn("contam_P00761", em)     # trypsin
        k = self.cont("ker")["keratin_sample"]
        self.assertEqual((k["database"], k["exempt_accessions"], k["exempt_alone_accessions"]),
                         ("keratins_in_database", ["contam_O43790", "contam_P02438"],
                          ["contam_O43790"]))
        self.assertNotIn("unchecked_tags", k)

    def test_contam_entries_the_database_checked_never_held_are_said(self):
        self.assertNotIn("contam_O43790", self.matrix("unchecked"))
        k = self.cont("unchecked")["keratin_sample"]
        self.assertEqual((k["database"], k["unchecked_tags"], k["caution"]),
                         ("keratins_removed_at_build", ["contam_"], True))
        self.assertEqual(k["n_precursors_unchecked"], self.cont("plain")["n_precursors"])
        self.assertIn("contaminants the skill did not add (a database built with "
                      "FragPipe/Philosopher's --contam)", k["note"])
        self.assertNotIn("adds its own contaminant list", k["note"])
        self.assertIn("CAUTION:", self.methods("unchecked"))

class MethodsReportAndReplay(unittest.TestCase):
    def test_de_sentence_per_database(self):
        built = mm.de_paragraph(de_prov(value=True, database="keratins_removed_at_build",
                                        n_removed_at_build=46, exempt_accessions=[]))
        self.assertIn("The samples are keratinous (a keratin sample), so keratin is the analyte",
                      built)
        self.assertIn("built without its 46 keratin-family contaminant entries", built)
        kept = mm.de_paragraph(de_prov(value=True, database="keratins_in_database",
                                       n_precursors_kept=4,
                                       exempt_accessions=["Cont_Q61765", "Cont_P02438"]))
        self.assertIn("4 precursors mapping only to keratin-family contaminant entries (2 in "
                      "the search database) were kept", kept)
        unknown = mm.de_paragraph(de_prov(value=True, database="unknown", note="x"))
        self.assertIn("could not be identified", unknown)
        self.assertIn("resolve before publication", unknown)
        for other in (de_prov(), de_prov(value=False, exempt_accessions=[]),
                      de_prov(value=None, note="not recorded")):
            self.assertNotIn("keratin", mm.de_paragraph(other).lower())

    def test_database_sentence(self):
        self.assertIn("the 46 keratin-family entries",
                      mm.keratin_database_sentence({"keratin_sample": True,
                                                    "n_contaminants_dropped_keratin_sample": 46}))
        self.assertEqual(mm.keratin_database_sentence({"keratin_sample": False}), "")
        self.assertEqual(mm.keratin_database_sentence({}), "")

    def test_methods_md(self):
        with tempfile.TemporaryDirectory() as tmp:
            meta = _write(os.path.join(tmp, "m.json"), json.dumps(dict(
                SIDECAR, proteome="UP000005640", uniprot_release="2026_03", organism="Homo sapiens",
                content_used="one_per_gene", n_proteome=20652, keratin_sample=True,
                n_contaminants_dropped_keratin_sample=46)))
            de = os.path.join(tmp, "de")
            os.makedirs(de)
            _write(os.path.join(de, "de_provenance.json"), json.dumps(de_prov(
                value=True, database="keratins_in_database", n_precursors_kept=4,
                exempt_accessions=["Cont_Q61765"], caution=True,
                note="the search database was not built ...")))
            out = os.path.join(tmp, "methods.md")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                                "--raw", os.path.join(tmp, "Ex_hair.raw"), "--fasta-meta", meta,
                                "--de-dir", de, "--out", out], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out) as fh:
                text = fh.read()
        self.assertIn("so the 46 keratin-family entries (keratins and keratin-associated "
                      "proteins) were also removed from the contaminant sequences", text)
        self.assertIn("> Keratin-sample caveat (resolve before publication): the search database "
                      "was not built", text)

    def test_report_notice_and_flags(self):
        self.assertIsNone(mah.keratin_note(de_prov()))
        self.assertIsNone(mah.keratin_note(de_prov(value=False, exempt_accessions=[])))
        n = mah.keratin_note(de_prov(value=True, database="keratins_removed_at_build",
                                     n_removed_at_build=46, exempt_accessions=[]))
        self.assertEqual(n["kind"], "info")
        self.assertIn("keratin is the analyte", n["text"])
        w = mah.keratin_note(de_prov(value=True, database="keratins_in_database",
                                     n_precursors_kept=4, exempt_accessions=["Cont_Q61765"],
                                     caution=True, note="rebuild the FASTA"))
        self.assertEqual(w["kind"], "warning")
        self.assertIn("Rebuild the FASTA", w["text"])
        kept = mah.keratin_kept_accessions(de_prov(value=True, exempt_accessions=["Cont_Q61765"]))
        self.assertIsNone(mah.background_flag("Krt31", "Cont_Q61765", None, kept))
        self.assertEqual(mah.background_flag("", "Cont_Q61765;Cont_P00761", None, kept),
                         "contaminant")
        self.assertEqual(mah.background_flag("", "Cont_Q61765"), "contaminant")

    def _repro(self, tmp, fasta_info):
        wf = os.path.join(tmp, "wf")
        env = os.path.join(tmp, "env.json")
        with open(env, "w") as fh:
            subprocess.run(["bash", os.path.join(SCRIPTS, "detect_env.sh")], stdout=fh,
                           stderr=subprocess.DEVNULL, check=True)
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

    def test_reproduce_sh_replays_the_keratin_build(self):
        base = {"proteome": "UP000005640", "content_used": "one_per_gene",
                "contaminant_set": "universal", "contaminant_target_rule":
                ff.CONTAMINANT_TARGET_RULE, "min_unique_peptides": 2}
        with tempfile.TemporaryDirectory() as tmp:
            sh = self._repro(tmp, dict(base, keratin_sample=True))
        fetch = next(ln for ln in sh.splitlines() if "--content" in ln and "--enzyme" in ln)
        self.assertIn("--keratin-sample", fetch)
        search = next(ln for ln in sh.splitlines() if "--fasta ./search.fasta" in ln)
        self.assertIn("--keratin-sample", search)
        with tempfile.TemporaryDirectory() as tmp:
            sh = self._repro(tmp, dict(base, keratin_sample=False, keratin_sample_source="user"))
        self.assertIn("--no-keratin-sample", sh)
        self.assertNotIn(" --keratin-sample", sh)
        with tempfile.TemporaryDirectory() as tmp:        # a default is replayed as a default
            sh = self._repro(tmp, dict(base, keratin_sample=False, keratin_sample_source="default"))
        self.assertNotIn("keratin-sample", sh)


class AuditorsReadTheSidecar(unittest.TestCase):
    def test_audit_does_not_flag_keratin_for_a_keratin_built_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            de = os.path.join(tmp, "de")
            os.makedirs(de)
            rows = "".join(f"P{i:05d},G{i},X_HUMAN,{20 - i * 0.01},{20 - i * 0.01}\n"
                           for i in range(600))
            _write(os.path.join(de, "Expression_Matrix.csv"),
                   "Protein.Group,Genes,Protein.Names,S1,S2\nQ15323,KRT31,K1H1_HUMAN,30,30\n" + rows)
            findings = {}
            for key, ker in (("keratin", True), ("plain", False)):
                meta = _write(os.path.join(tmp, f"{key}.json"),
                              json.dumps(dict(SIDECAR, keratin_sample=ker)))
                r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                                    "--out", os.path.join(tmp, f"AUDIT_{key}.md"), "--de-dir", de,
                                    "--fasta-meta", meta], capture_output=True, text=True, cwd=tmp)
                self.assertEqual(r.returncode, 0, r.stderr)
                with open(os.path.join(tmp, f"AUDIT_{key}.json")) as fh:
                    findings[key] = next(f for f in json.load(fh)["findings"]
                                         if f["check"] == "contamination")
        self.assertIn("Keratin-matrix sample — keratins are the analyte",
                      findings["keratin"]["message"])
        self.assertIn("(keratins/trypsin/casein)", findings["plain"]["message"])
        self.assertEqual(findings["plain"]["detail"]["examples"], ["KRT31"])


class SkillDocAsksAtStep3(unittest.TestCase):
    """SKILL.md decides a keratin sample BEFORE the database is built, and every script it
    tells to take --keratin-sample does."""

    @classmethod
    def setUpClass(cls):
        with open(os.path.join(SKILL, "SKILL.md")) as fh:
            cls.doc = fh.read()

    def section(self, head):
        start = self.doc.index(head)
        nxt = re.search(r"^### ", self.doc[start + len(head):], re.M)
        return self.doc[start:start + len(head) + (nxt.start() if nxt else len(self.doc))]

    def test_step3_asks_the_sample_type(self):
        s3 = self.section("### 3. Ask organism + experimental design")
        # asked for every analysis since 2.10 (tests/test_sample_type_asked.py); keratin is
        # decided from that answer
        self.assertIn("**Sample type — ask, for every analysis, and never infer it.**", s3)
        self.assertIn("**Is the sample keratin?** Decide it from the sample type the user gave",
                      s3)
        for tissue in ("hair", "wool", "feather", "skin", "nail"):
            self.assertIn(tissue, s3)
        self.assertIn("Deciding it after the search is too late", s3)
        self.assertIn("fetch_fasta.KERATIN_SAMPLE_TISSUES", s3)

    def test_each_step_carries_it(self):
        s6 = self.section("### 6. Build the FASTA")
        self.assertIn("--keratin-sample|--no-keratin-sample --out ./search.fasta", s6)
        self.assertIn("DEFAULT — not user-confirmed", s6)
        s7 = self.section("### 7. Run the search")
        self.assertIn("[--keratin-sample]", s7)
        self.assertIn("REFUSED", s7)
        self.assertIn("pass `--keratin-sample`", self.section("### 8. Differential expression"))
        self.assertIn("read it from `--fasta-meta`",
                      self.section("### 8c. Audit the results"))

    def test_the_limits_are_documented(self):
        s3 = self.section("### 3. Ask organism + experimental design")
        for part in ("KRT1/2/9/10", "FLG/FLG2/HRNR/DMKN", "Poorly annotated species"):
            self.assertIn(part, s3)
        s8 = self.section("### 8. Differential expression")
        self.assertIn("refuses a FragPipe search on such a database", s8)
        self.assertNotIn("FragPipe adds when it searches", s8)

    def test_every_documented_flag_exists(self):
        for script in ("fetch_fasta.py", "run_search.py", "audit_results.py", "sample_quality.py"):
            with self.subTest(script=script):
                with open(os.path.join(SCRIPTS, script)) as fh:
                    self.assertIn('"--keratin-sample"', fh.read())
        with open(os.path.join(SCRIPTS, "run_de.R")) as fh:
            self.assertIn('getarg("--keratin-sample"', fh.read())
        for t in ("hair", "wool", "feather", "nail"):
            self.assertIn(t, ff.KERATIN_SAMPLE_TISSUES)


if __name__ == "__main__":
    unittest.main(verbosity=2)
