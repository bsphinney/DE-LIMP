#!/usr/bin/env python3
"""
The FRAN database guard: no sequence that does not exist may reach the Core's corpus.

fasta_guard.py is the one definition of the rules; fran_deposit.check() applies them to the FASTA
files a search used (stage AND backfill), and a session can opt out wholesale. Everything here is
built from tiny synthetic FASTAs in a temp dir.
"""
import json
import os
import sys
import tempfile
import unittest
import unittest.mock
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))
sys.path.insert(0, HERE)

import fasta_fixture  # noqa: E402
import fasta_guard as fg  # noqa: E402
import fran_deposit as fd  # noqa: E402
from test_fran_deposit_gate import Args, search_dir  # noqa: E402

PROT = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV"


def fasta(path, entries, sidecar=True):
    with open(path, "w") as fh:
        for h, s in entries:
            fh.write(f">{h}\n{s}\n")
    if sidecar:
        with open(path + ".meta.json", "w") as fh:
            json.dump({"organism": "", "taxid": 0, "organism_source": "none"}, fh)
    return path


def real(i, seq=PROT, name="Some protein"):
    return (f"sp|P{i:05d}|G{i}_HUMAN {name} GN=G{i} PE=1 SV=1", seq)


def verdict(path):
    return fg.classify(path)[0]


def setUpModule():
    global _SAVED
    _SAVED = {k: os.environ.pop(k, None) for k in ("FRAN_DEPOSIT", "FRAN_DROP_DIR")}
    os.environ["FRAN_CORPUS_QUERY"] = "off"


def tearDownModule():
    for k, v in _SAVED.items():
        if v is not None:
            os.environ[k] = v


class RuleTests(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = self._t.name
        self.addCleanup(self._t.cleanup)

    def db(self, entries, **kw):
        return fasta(os.path.join(self.d, "x.fasta"), entries, **kw)

    def test_a_plain_database_passes(self):
        self.assertEqual(verdict(self.db([real(1), real(2)])), "ok")

    def test_the_test_header_patterns_are_hard(self):
        for head in ("ENTRAP_P00001", "sp|P1|shuffled_X", "sp|P1|X_ENT_3", "decoy_P1", "x_VAR_12",
                     "x_CON_1", "x_ENC_2", "Entrap_lower", "sp|P1|X_HUMAN protein ENTRAP_9"):
            with self.subTest(head=head):
                self.assertEqual(verdict(self.db([real(1), (head, PROT)])), "test_database")

    def test_ent_tag_is_case_sensitive_like_the_other_variant_tags(self):
        # a real custom proteome: genome-annotation scaffolds named scaffold_ent_12 (review M2)
        self.assertEqual(verdict(self.db([real(1), ("scaffold_ent_12 gene model", PROT)])), "ok")
        self.assertEqual(verdict(self.db([real(1), ("scaffold_ENT_12 gene model", PROT)])),
                         "test_database")
        self.assertFalse(fg.is_synthetic_accession("scaffold_ent_12"))
        self.assertTrue(fg.is_synthetic_accession("P1_ENT_12"))

    def test_the_residue_constants_are_what_the_rule_reads(self):
        with mock.patch.object(fg, "HARD_RESIDUES", "JW"):
            self.assertEqual(verdict(self.db([real(1), ("sp|P9|X", "MKWAK" + PROT)])),
                             "test_database")

    def test_added_sequences_ask_and_are_scanned_with_the_same_rules(self):
        base = self.db([real(1)])
        added = fasta(os.path.join(self.d, "bait.fasta"), [("bait1 EGFP", PROT[::-1])], sidecar=False)
        meta = {"organism": "", "taxid": 0, "organism_source": "none",
                "added_sequences": [{"file": added, "sha256": "x"}]}
        with open(base + ".meta.json", "w") as fh:
            json.dump(meta, fh)
        v, why, det = fg.classify(base)
        self.assertEqual(v, "needs_agent_check", why)       # an untagged, shuffled addition
        self.assertTrue(any("bait.fasta" in w for w in why), why)
        self.assertEqual(det["added_files"], [added])
        bad = fasta(os.path.join(self.d, "bad.fasta"), [("ENTRAP_X1", PROT), ("b2", "MKJA")],
                    sidecar=False)
        meta["added_sequences"] = [{"file": bad}]
        with open(base + ".meta.json", "w") as fh:
            json.dump(meta, fh)
        self.assertEqual(verdict(base), "test_database")
        meta["added_sequences"] = []                        # none added: back to ok
        with open(base + ".meta.json", "w") as fh:
            json.dump(meta, fh)
        self.assertEqual(verdict(base), "ok")

    def test_a_real_protein_NAMED_decoy_is_not_a_decoy(self):
        # Arabidopsis UniProt: a protein named "Decoy" with GN=DECOY (measured on the Core's cache)
        arath = "tr|F4HW90|F4HW90_ARATH Decoy OS=Arabidopsis thaliana OX=3702 GN=DECOY PE=1 SV=1"
        self.assertEqual(verdict(self.db([real(1), (arath, PROT)])), "ok")

    def test_real_lowercase_con_in_a_genome_annotation_is_not_a_variant_tag(self):
        # measured in the Lance corpus: a kelp annotation names contigs ...Genewise1Plus.C_con_440104
        kelp = "Mp-1236660-estExt_Genewise1Plus.C_con_440104 some protein"
        self.assertEqual(verdict(self.db([real(1), (kelp, PROT)])), "ok")
        self.assertFalse(fg.is_synthetic_accession("Mp-1237502-estExt_Genewise1Plus.C_con_630003"))
        self.assertTrue(fg.is_synthetic_accession("P12345_CON_7QE"))
        self.assertTrue(fg.is_synthetic_accession("entrap_p1"))

    def test_j_is_hard_but_u_o_b_z_only_ask(self):
        self.assertEqual(verdict(self.db([real(1), real(2, "MKJA")])), "test_database")
        for res in "UOBZ":
            with self.subTest(res=res):
                self.assertEqual(verdict(self.db([real(1), real(2, "MK" + res + "AK")])),
                                 "needs_agent_check")

    def test_a_stray_B_or_Z_in_a_big_real_proteome_is_tolerated_a_recoded_set_is_not(self):
        big = [real(i) for i in range(3, 2003)]              # 2,000 entries
        self.assertEqual(verdict(self.db(big + [real(1, "MKBAK"), real(2, "MKZAK")])), "ok")
        many = big + [real(i, "MKBAK") for i in range(5000, 5020)]   # 1%: machine-made
        self.assertEqual(verdict(self.db(many)), "needs_agent_check")

    def test_an_empty_database_asks(self):
        self.assertEqual(verdict(self.db([])), "needs_agent_check")

    def test_selenoproteins_with_nonobvious_names_are_explained(self):
        for name in ("Ethanolaminephosphotransferase 1", "Protein adenylyltransferase SelO, mitochondrial",
                     "Selenide, water dikinase 2", "Selenoprotein P", "Glutathione peroxidase 4",
                     "Thioredoxin reductase 1", "Iodothyronine deiodinase type 2"):
            with self.subTest(name=name):
                self.assertEqual(verdict(self.db([real(1), real(2, "MKUAK", name)])), "ok")

    def test_selenocysteine_in_a_selenoprotein_passes_in_another_entry_asks(self):
        sel = ("sp|P07203|GPX1_HUMAN Glutathione peroxidase 1 GN=GPX1 PE=1 SV=2", "MKUAK")
        self.assertEqual(verdict(self.db([real(1), sel])), "ok")
        r = fg.classify(self.db([real(1), sel, real(2, "MKUAK", "Kinase")]))
        self.assertEqual(r[0], "needs_agent_check")
        self.assertEqual(r[2]["u_unexplained"]["n"], 1)

    def test_pyrrolysine_header_explains_o(self):
        pyl = ("tr|Q1|MTMB_METBA Monomethylamine methyltransferase MtmB", "MKOAK")
        self.assertEqual(verdict(self.db([real(1), pyl])), "ok")

    def test_no_fetch_fasta_sidecar_asks(self):
        r = fg.classify(self.db([real(1)], sidecar=False))
        self.assertEqual(r[0], "needs_agent_check")
        self.assertTrue(r[2]["no_fetch_fasta_sidecar"])

    def test_the_nofran_sidecar_is_hard_with_or_without_content_evidence(self):
        p = self.db([real(1)])
        fg.mark_nofran(p, "recoded", "test", "pS recoded as U")
        self.assertEqual(verdict(p), "test_database")
        # unreadable marker reads as "no", never as permission
        with open(fg.sidecar_path(p), "w") as fh:
            fh.write("{not json")
        self.assertEqual(verdict(p), "test_database")
        # a corpus_eligible TRUE file is not a refusal
        with open(fg.sidecar_path(p), "w") as fh:
            json.dump({"corpus_eligible": True}, fh)
        self.assertEqual(verdict(p), "ok")

    def test_a_test_database_copied_without_its_marker_is_still_refused(self):
        """C2 review ux M7: fetch_fasta records the marker's path in the FASTA's .meta.json; a copy
        that kept the .meta.json and lost the .nofran.json must not read as a real database."""
        p = self.db([real(1), real(2)])
        with open(p + ".meta.json") as fh:
            meta = json.load(fh)
        meta["nofran_sidecar"] = "/elsewhere/x.fasta.nofran.json"
        with open(p + ".meta.json", "w") as fh:
            json.dump(meta, fh)
        v, reasons, detail = fg.classify(p)
        self.assertEqual(v, "test_database", reasons)
        self.assertEqual(detail["meta_nofran_sidecar"], "/elsewhere/x.fasta.nofran.json")
        # a marker that IS there and says eligible is a person's decision: not overruled
        with open(fg.sidecar_path(p), "w") as fh:
            json.dump({"corpus_eligible": True}, fh)
        self.assertEqual(verdict(p), "ok")
        # a record whose nofran_sidecar is null (a real fetch) is no refusal
        os.remove(fg.sidecar_path(p))
        meta["nofran_sidecar"] = None
        with open(p + ".meta.json", "w") as fh:
            json.dump(meta, fh)
        self.assertEqual(verdict(p), "ok")

    def test_the_marker_is_not_the_meta_json_and_bad_kinds_are_refused(self):
        p = self.db([real(1)])
        with open(p + ".meta.json") as fh:
            before = fh.read()
        out = fg.mark_nofran(p, "entrapment", "t", "w")
        self.assertTrue(out.endswith(".nofran.json"))
        with open(p + ".meta.json") as fh:
            self.assertEqual(fh.read(), before)
        with self.assertRaises(ValueError):
            fg.mark_nofran(p, "bogus", "t", "w")

    def test_a_complete_target_decoy_pair_set_is_fragpipes_ordinary_database(self):
        ents = [real(1), real(2)]
        ents += [("rev_" + h, s[::-1]) for h, s in ents]
        self.assertEqual(verdict(self.db(ents)), "ok")

    def test_decoys_without_targets_are_a_test_database(self):
        ents = [("rev_" + h, s[::-1]) for h, s in (real(1), real(2))]
        self.assertEqual(verdict(self.db(ents)), "test_database")
        # rev_rev_: the decoys of a decoy set
        full = [real(1)] + [("rev_" + real(1)[0], PROT[::-1])]
        full += [("rev_rev_" + real(1)[0], PROT)]
        self.assertEqual(verdict(self.db(full)), "test_database")

    def test_the_unreadable_database_asks_it_does_not_refuse(self):
        self.assertEqual(fg.classify(os.path.join(self.d, "missing.fasta"))[0], "needs_agent_check")


class CheckTests(unittest.TestCase):
    """fran_deposit.check() applies the guard to the FASTAs the search used."""
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = self._t.name
        self.addCleanup(self._t.cleanup)
        os.environ["FRAN_DROP_DIR"] = os.path.join(self.d, "incoming")
        self.addCleanup(os.environ.pop, "FRAN_DROP_DIR", None)

    def search(self, entries, sidecar=True, engine="diann"):
        out = search_dir(self.d, engine)
        db = fasta(os.path.join(self.d, "db.fasta"), entries, sidecar=sidecar)
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": engine, "version": "9.9", "fasta": db}, fh)
        return out, db

    def check(self, out, **kw):
        kw.setdefault("database_ok", None)
        return fd.check(Args(out, **kw))

    def test_a_real_database_is_eligible(self):
        out, _ = self.search([real(1)])
        c = self.check(out)
        self.assertTrue(c["eligible"], c.get("detail"))

    def test_an_unmarked_entrapment_database_is_refused_by_content(self):
        out, _ = self.search([real(1), ("ENTRAP_P00001", PROT[::-1])])
        c = self.check(out)
        self.assertEqual(c["reason"], "test_database")
        self.assertFalse(c["eligible"])

    def test_a_marked_database_is_refused_and_the_refusal_is_recorded_for_backfill(self):
        out, db = self.search([real(1)])
        fg.mark_nofran(db, "denovo_candidate", "t", "w")
        a = Args(out, database_ok=None)
        res = fd.do_stage(a)
        self.assertEqual(res["reason"], "test_database")
        self.assertFalse(res["staged"])
        self.assertEqual(fd.read_receipt(out)["status"], "test_database")
        # ...and a later flagless stage does not overturn it into a hand-over
        self.assertFalse(fd.do_stage(Args(out, database_ok=None))["staged"])

    def test_the_hand_built_database_asks_and_database_ok_is_recorded(self):
        out, _ = self.search([real(1)], sidecar=False)
        c = self.check(out)
        self.assertEqual(c["reason"], "needs_agent_check")
        self.assertIsNone(fd.read_receipt(out))               # nothing recorded
        res = fd.do_stage(Args(out, database_ok="client's own construct database, checked"))
        self.assertTrue(res["staged"], res.get("detail"))
        with open(os.path.join(res["entry"], fd.MANIFEST)) as fh:
            man = json.load(fh)
        self.assertEqual(man["database_ok"]["why"], "client's own construct database, checked")
        self.assertEqual(fd.read_receipt(out)["database_ok"]["why"],
                         "client's own construct database, checked")

    def test_a_skill_search_whose_sidecar_moved_is_still_a_skill_fetched_database(self):
        out, db = self.search([real(1)], sidecar=False)
        prov = os.path.join(out, "search_provenance.json")
        for fm, expect in (({"path": "/gone/search.fasta.meta.json", "rule": "own"}, True),
                           ({"path": None, "rule": "given"}, False), (None, False)):
            with open(prov, "w") as fh:
                json.dump({"engine": "diann", "version": "9.9", "fasta": db, "fasta_meta": fm}, fh)
            c = self.check(out)
            self.assertEqual(c["eligible"], expect, (fm, c.get("detail")))

    def test_the_provenance_route_never_clears_a_test_database(self):
        out, db = self.search([real(1), ("ENTRAP_P1", PROT)], sidecar=False)
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "fasta": db, "fasta_meta": {"path": "/x.meta.json"}}, fh)
        self.assertEqual(self.check(out)["reason"], "test_database")

    def test_database_ok_never_overrides_a_hard_refusal(self):
        out, _ = self.search([real(1, "MKJA")])
        self.assertEqual(self.check(out, database_ok="trust me")["reason"], "test_database")

    def test_a_wrong_hard_refusal_names_its_remedy(self):
        out, _ = self.search([real(1, "MKJA")])
        d = self.check(out, database_ok="trust me")["detail"]
        self.assertIn("does NOT clear", d)
        self.assertIn("NO_FRAN", d)

    def test_database_ok_needs_a_real_reason(self):
        out, _ = self.search([real(1)], sidecar=False)
        for blank in ("", "   "):
            c = self.check(out, database_ok=blank)
            self.assertEqual(c["reason"], "needs_agent_check", blank)
            self.assertNotIn("database_ok", c)
        self.assertTrue(self.check(out, database_ok="bait construct, looked at it")["eligible"])

    def test_database_ok_records_who_when_and_whether_unattended(self):
        out, _ = self.search([real(1)], sidecar=False)
        with mock.patch.dict(os.environ, {"SLURM_JOB_ID": "123"}):
            c = self.check(out, database_ok="own construct database, read it")
        rec = c["database_ok"]
        self.assertEqual((rec["by"], rec["unattended_job"]), (c["user"], "123"))
        self.assertRegex(rec["at"], r"^\d{4}-\d\d-\d\dT")

    def test_a_search_that_names_no_database_asks(self):
        out = search_dir(self.d)
        os.unlink(os.path.join(out, "search_provenance.json"))
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "version": "9.9"}, fh)
        self.assertEqual(self.check(out)["reason"], "needs_agent_check")

    def test_fragpipe_with_a_complete_decoy_set_is_eligible(self):
        ents = [real(1), real(2)]
        ents += [("rev_" + h, s[::-1]) for h, s in ents]
        out, _ = self.search(ents, engine="fragpipe")
        c = self.check(out)
        self.assertTrue(c["eligible"], c.get("detail"))

    def test_every_fasta_the_search_used_is_checked(self):
        out, _ = self.search([real(1)])
        bad = fasta(os.path.join(self.d, "second.fasta"), [real(2), ("x_VAR_1", PROT)])
        with open(os.path.join(out, "report.log.txt"), "w") as fh:
            fh.write(f"diann-linux --fasta {bad} --threads 4\n")
        self.assertEqual(self.check(out)["reason"], "test_database")


class SessionOptOutTests(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.d = self._t.name
        self.addCleanup(self._t.cleanup)
        os.environ["FRAN_DROP_DIR"] = os.path.join(self.d, "incoming")
        self.addCleanup(os.environ.pop, "FRAN_DROP_DIR", None)

    def session(self, **kw):
        sess = os.path.join(self.d, "service", "LabX", "sess")
        os.makedirs(os.path.join(sess, "output"), exist_ok=True)
        out = search_dir(os.path.join(sess, "output"), "search")
        return sess, out

    def test_a_NO_FRAN_file_at_any_ancestor_refuses(self):
        for where in ("sess", "."):
            with self.subTest(where=where):
                sess, out = self.session()
                target = sess if where == "sess" else os.path.join(self.d, "service")
                marker = os.path.join(target, "NO_FRAN")
                open(marker, "w").close()
                c = fd.check(Args(out))
                self.assertEqual(c["reason"], "opted_out_session")
                self.assertEqual(c["opt_out_file"], marker)
                os.unlink(marker)

    def test_session_json_fran_false_refuses_and_true_or_absent_does_not(self):
        sess, out = self.session()
        p = os.path.join(sess, "session.json")
        with open(p, "w") as fh:
            json.dump({"fran": False}, fh)
        self.assertEqual(fd.check(Args(out))["reason"], "opted_out_session")
        with open(p, "w") as fh:
            json.dump({"fran": True}, fh)
        self.assertNotEqual(fd.check(Args(out))["reason"], "opted_out_session")
        with open(p, "w") as fh:
            json.dump({"qc": False}, fh)
        self.assertNotEqual(fd.check(Args(out))["reason"], "opted_out_session")

    def test_a_malformed_session_json_is_said_and_asks(self):
        sess, out = self.session()
        with open(os.path.join(sess, "session.json"), "w") as fh:
            fh.write('{"fran": false,')
        import contextlib
        import io
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            c = fd.check(Args(out))
        self.assertEqual(c["reason"], "needs_agent_check")
        self.assertEqual(c["session_json_unreadable"], [os.path.join(sess, "session.json")])
        self.assertIn("cannot be read", err.getvalue())
        # a BOM is not malformed
        with open(os.path.join(sess, "session.json"), "w", encoding="utf-8-sig") as fh:
            json.dump({"fran": False}, fh)
        self.assertEqual(fd.check(Args(out))["reason"], "opted_out_session")

    def test_stage_records_it_and_backfill_honours_it(self):
        sess, out = self.session()
        open(os.path.join(sess, "NO_FRAN"), "w").close()
        res = fd.do_stage(Args(out))
        self.assertEqual(res["reason"], "opted_out_session")
        self.assertEqual(fd.read_receipt(out)["status"], "opted_out_session")
        rows = fd.plan_backfill([out], apply=True, prefixes=(os.path.realpath(self.d) + "/",))
        self.assertEqual(rows[0]["decision"], "skip")
        self.assertIn(rows[0]["reason"], ("opted_out_session",))
        self.assertFalse(os.path.isdir(os.path.join(self.d, "incoming", fd.entry_name(out))))

    def test_backfill_honours_it_even_with_no_receipt(self):
        sess, out = self.session()
        open(os.path.join(sess, "NO_FRAN"), "w").close()
        rows = fd.plan_backfill([out], apply=False, prefixes=(os.path.realpath(self.d) + "/",))
        self.assertEqual(rows[0]["reason"], "opted_out_session")


class WriterTests(unittest.TestCase):
    def test_fetch_fasta_writes_the_marker_with_test_database_and_not_otherwise(self):
        import subprocess
        scripts = os.path.join(os.path.dirname(HERE), "scripts")
        with tempfile.TemporaryDirectory() as d:
            src = fasta(os.path.join(d, "in.fasta"), [real(1)], sidecar=False)
            for flag, marked in ((["--test-database", "entrapment", "--test-why", "unit test"], True),
                                 ([], False)):
                out = os.path.join(d, "out_%s.fasta" % marked)
                r = subprocess.run([sys.executable, os.path.join(scripts, "fetch_fasta.py"), "fetch",
                                    "--path", src, "--organism", "none", "--contaminants", "none",
                                    "--out", out] + flag, capture_output=True, text=True)
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual(os.path.isfile(out + ".nofran.json"), marked)
                self.assertTrue(os.path.isfile(out + ".meta.json"))
            marker = fg.read_nofran(os.path.join(d, "out_True.fasta"))
            self.assertEqual((marker["kind"], marker["why"]), ("entrapment", "unit test"))


class AuditTests(unittest.TestCase):
    def test_the_audit_prints_read_only_sql_and_never_connects(self):
        import contextlib
        import io
        a = type("A", (), {"lance_dir": None, "report": None, "allow_login_node": False})()
        buf = io.StringIO()
        with unittest.mock.patch.dict(sys.modules, {"psycopg2": None}), \
                contextlib.redirect_stdout(buf), self.assertRaises(SystemExit):
            fd.audit(a)
        res = json.loads(buf.getvalue())
        self.assertIn("NOT CONTACTED", res["pg_farm"])
        for q in res["sql"]:
            self.assertRegex(q["sql"], r"^SELECT ")
            self.assertNotRegex(q["sql"], r"(?i)\b(insert|update|delete|drop|alter)\b")
        names = [q["name"] for q in res["sql"]]
        self.assertEqual(names[:4], ["peptides_with_J", "peptides_with_O_B_Z",
                                     "U_peptides_outside_selenoprotein_genes",
                                     "synthetic_protein_groups"])
        # the SQL carries the SAME accession pattern as the guard (one definition), PG-flavoured
        self.assertIn("ENTRAP_|shuffle", res["sql"][3]["sql"])
        self.assertIn("_VAR_|_CON_|_ENC_", res["sql"][3]["sql"])
        self.assertIn("how", res["lance_scan"])

    def test_the_lance_scan_refuses_a_login_node(self):
        import contextlib
        import io
        a = type("A", (), {"lance_dir": "/nonexistent", "report": None, "allow_login_node": False})()
        buf = io.StringIO()
        with unittest.mock.patch.dict(os.environ, {}, clear=False), contextlib.redirect_stdout(buf):
            os.environ.pop("SLURM_JOB_ID", None)
            with self.assertRaises(SystemExit):
                fd.audit(a)
        self.assertEqual(json.loads(buf.getvalue())["lance_scan"]["refused"], "login_node")

    def test_the_lance_scan_finds_planted_rows(self):
        try:
            import lance  # noqa: F401
            import pyarrow as pa
        except ImportError:
            self.skipTest("lance/pyarrow not installed here (HIVE: alphadia2 env)")
        with tempfile.TemporaryDirectory() as d:
            t = pa.table({"search_id": ["a"] * 4, "stripped_seq": ["PEPTIDEK", "PEPJK", "PEPOK", "AUK"],
                          "protein_group": ["P1", "P2", "ENTRAP_P3", "P4"],
                          "genes": ["G1", "G2", "G3", "XYZ"]})
            lance.write_dataset(t, os.path.join(d, "s1.lance"))
            r = fg.lance_scan(d)
        self.assertEqual((r["J"], r["OBZ"], r["U_outside"], r["synthetic_accessions"]), (1, 1, 1, 1))


class DeclaredPatternsTests(unittest.TestCase):
    def test_the_accession_pattern_is_pinned_for_frans_ingest(self):
        # tests/test_synthetic_guard.py in the FRAN repo pins the same literal; change both or neither
        self.assertEqual(fg.ACCESSION_TEST_CI, r"ENTRAP_|shuffle")
        self.assertEqual(fg.ACCESSION_TEST_CS, r"_VAR_|_CON_|_ENC_|_ENT_\d")
        self.assertEqual((fg.HARD_RESIDUES, fg.AMBIGUOUS_RESIDUES), ("J", "BZ"))

    def test_one_definition_of_the_header_patterns(self):
        # fran_deposit must not carry its own copy of the rules
        with open(os.path.join(os.path.dirname(HERE), "scripts", "fran_deposit.py")) as fh:
            src = fh.read()
        for needle in ("ENTRAP_|shuffle", "_VAR_|_CON_|_ENC_"):
            self.assertNotIn(needle, src)


if __name__ == "__main__":
    unittest.main()
