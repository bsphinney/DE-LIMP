#!/usr/bin/env python3
"""2.15 lane X8: FragPipe identifications-only / no-DE polish (F2, F4-F10, F13, F14 of the
PROT_9812 after-run). One class per finding."""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import fragpipe_params as fpp  # noqa: E402
import make_analysis_html as mah  # noqa: E402
import sage_lfq_check as slc  # noqa: E402
from test_fragpipe_dda_stage import FP_WORKFLOW  # noqa: E402


def write(path, obj):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh)


class F4DistinctPeptides(unittest.TestCase):
    def test_adapt_records_distinct_peptides_per_run(self):
        import run_search
        with tempfile.TemporaryDirectory() as out:
            with open(os.path.join(out, "psm.tsv"), "w") as fh:
                fh.write("Spectrum\tPeptide\tCharge\tProbability\tAssigned Modifications\t"
                         "Protein ID\tMapped Proteins\n"
                         "S1-a.00007.00007.2\tPEPTIDEK\t2\t0.99\t\tP02768\t\n"
                         "S1-a.00008.00008.2\tPEPTIDEK\t2\t0.99\t15.9949 M\tP02768\t\n"
                         "S1-a.00009.00009.2\tOTHERK\t2\t0.99\t\tP02768\t\n")
            with open(os.path.join(out, "peptide.tsv"), "w") as fh:
                fh.write("Peptide\tProbability\nPEPTIDEK\t0.99\nOTHERK\t0.99\n")
            q = fpp.workflow_file_quantification(FP_WORKFLOW)
            run_search.adapt_fragpipe_ids(out, q)
            with open(os.path.join(out, "fragpipe_identifications.json")) as fh:
                rec = json.load(fh)
        self.assertEqual(rec["per_run"], {"S1-a": 3})            # forms: unchanged
        self.assertEqual(rec["per_run_peptides"], {"S1-a": 2})   # distinct peptides

    def _note(self, rec):
        with tempfile.TemporaryDirectory() as s:
            d = os.path.join(s, "output", "search_fragpipe")
            write(os.path.join(d, "search_provenance.json"), {
                "engine": "fragpipe", "quantification": slc.not_applicable(engine="fragpipe")})
            write(os.path.join(d, "fragpipe_identifications.json"), rec)
            return mah.identifications_only_note({}, s)["text"]

    def test_the_page_shows_distinct_peptides_not_rows(self):
        t = self._note({"per_run": {"r1": 317}, "per_run_peptides": {"r1": 307}})
        self.assertIn("r1: 307 peptides", t)
        self.assertNotIn("317 peptides", t)
        self.assertIn("317 peptide forms", t)

    def test_an_older_record_says_forms_not_peptides(self):
        t = self._note({"per_run": {"r1": 40}})
        self.assertIn("r1: 40 peptide forms", t)
        self.assertNotIn("40 peptides", t)


class F6EnzymeName(unittest.TestCase):
    def test_list_entry_lookup(self):
        import cleavage_rules as cr
        self.assertEqual(cr.fragpipe_list_entry("FLWY", "P", "C"), "chymotrypsin")
        self.assertEqual(cr.fragpipe_list_entry("RK", "P", "c"), "trypsin")
        self.assertIsNone(cr.fragpipe_list_entry("FLMWY", "P", "C"))      # H5's own rule
        self.assertIsNone(cr.fragpipe_list_entry("@", "", "C"))           # never `nonspecific`

    def test_the_written_workflow_names_the_enzyme(self):
        import subprocess
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "c.workflow")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_presets.py"),
                                "--engine", "fragpipe", "--instrument", "timsTOF HT",
                                "--acquisition", "DDA", "--ids-only", "--enzyme", "chymotrypsin",
                                "--out", out], capture_output=True, text=True, timeout=120)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            kv, _ = fpp.read_workflow(out)
        self.assertEqual(kv["msfragger.search_enzyme_name_1"], "chymotrypsin")
        self.assertEqual(kv["msfragger.misc.fragger.enzyme-dropdown-1"], "chymotrypsin")
        self.assertEqual(kv["msfragger.search_enzyme_cut_1"], "FLWY")


class F5ManifestRecordsWhatRan(unittest.TestCase):
    def resolve(self, *args):
        import subprocess
        with tempfile.TemporaryDirectory() as d:
            dest = os.path.join(d, "wf")
            env = {k: v for k, v in os.environ.items() if k != "UCD_CORE_TOOLS"}
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                                "--dest", dest, "--engine", "fragpipe", "--acquisition", "DDA",
                                "--instrument", "timsTOF HT", *args],
                               capture_output=True, text=True, timeout=120, env=env)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(os.path.join(dest, "workflow.manifest.json")) as fh:
                man = json.load(fh)
            kv, _ = fpp.read_workflow(os.path.join(dest, "fragpipe.workflow"))
        return man, json.loads(p.stdout), kv

    def test_ids_only_has_no_de_block_and_the_ppm_that_ran(self):
        man, out, kv = self.resolve("--ids-only")
        self.assertIsNone(man["de"])
        self.assertIsNone(out["de"])
        ran = float(kv["msfragger.precursor_mass_upper"])
        self.assertEqual(man["search"]["ms1_ppm"], ran)
        self.assertEqual(man["search"]["ms2_ppm"], float(kv["msfragger.fragment_mass_tolerance"]))
        self.assertNotIn("DIA-NN README", man["search"]["ppm_source"])
        self.assertIn("FragPipe workflow", man["search"]["ppm_source"])

    def test_a_quantifying_search_keeps_its_de_block(self):
        man, _, _ = self.resolve()
        self.assertEqual(man["de"]["method"], "maxlfq")


class F7ContaminantDigest(unittest.TestCase):
    def test_search_digest_follows_the_enzyme(self):
        import fetch_fasta as ff
        self.assertEqual(ff.search_digest(("chymotrypsin",))["cut"], "F*,L*,W*,Y*,!*P")
        self.assertEqual(ff.search_digest(("gluc",))["cut"], "D*,E*,!*P")
        self.assertEqual(ff.search_digest(("lysc", "trypsin")), ff.DIANN_DIGEST)   # the default
        self.assertEqual(ff.search_digest(()), ff.DIANN_DIGEST)                    # stand-in
        self.assertEqual(ff.search_digest(("pepsin",)), ff.DIANN_DIGEST)           # no rule
        # lengths and missed cleavages are the search's, only the cut moves
        d = ff.search_digest(("chymotrypsin",))
        self.assertEqual({k: v for k, v in d.items() if k != "cut"},
                         {k: v for k, v in ff.DIANN_DIGEST.items() if k != "cut"})

    def test_the_sidecar_records_the_digest_it_counted_in(self):
        import contextlib
        import io
        from unittest import mock
        import fetch_fasta as ff
        import test_fetch_fasta_contaminant_overlap as t
        with tempfile.TemporaryDirectory() as root:
            tgt, cp = os.path.join(root, "t.fasta"), os.path.join(root, "c.fasta")
            with open(tgt, "w") as fh:
                fh.write(t.TARGET)
            with open(cp, "w") as fh:
                fh.write(t.CONT)
            out = os.path.join(root, "out", "search.fasta")
            argv = ["fetch_fasta.py", "fetch", "--path", tgt, "--contaminants", "universal",
                    "--contaminants-path", cp, "--organism", "none", "--enzyme", "chymotrypsin",
                    "--out", out]
            with mock.patch.object(sys, "argv", argv), contextlib.redirect_stdout(io.StringIO()), \
                    contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(ff.main(), 0)
            with open(out + ".meta.json") as fh:
                meta = json.load(fh)
        self.assertEqual(meta["contaminant_digest"]["cut"], "F*,L*,W*,Y*,!*P")
        self.assertIn("F*,L*,W*,Y*,!*P", meta["contaminant_target_rule"])
        self.assertNotIn("--cut K*,R*", meta["contaminant_target_rule"])


class F2WholeDatabaseAsAdded(unittest.TestCase):
    def read(self, added_text):
        import fetch_fasta as ff
        import test_fetch_fasta_add_fasta as t
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "staff.fasta")
            with open(p, "w") as fh:
                fh.write(added_text)
            with self.assertRaises(SystemExit) as cm:
                ff.read_added_sequences([p], t.PROTEOME)
        return str(cm.exception.code)

    def test_a_database_holding_the_proteome_gets_the_right_advice(self):
        import test_fetch_fasta_add_fasta as t
        msg = self.read(t.PROTEOME + ">IN|1|BAIT\nMPEPTIDEKAAAR\n"
                        ">sp|Cont_X1|KRT\nMKKKKRRRR\n")
        self.assertIn("2 of its 4 entries are already in the database", msg)
        self.assertIn("--path", msg)
        self.assertIn("only those", msg)
        self.assertNotIn("rename", msg)

    def test_a_real_accession_clash_still_says_rename(self):
        msg = self.read(">sp|P60709|ACTB_MINE\nMPEPTIDEKAAAR\n")
        self.assertIn("rename the added entry", msg)


class F8DeDirOptional(unittest.TestCase):
    def test_the_no_de_brief_needs_no_de_dir(self):
        import subprocess
        import test_no_de_brief as nb

        class T(nb.Base):
            def runTest(self):
                pass
        t = T()
        t.setUp()
        try:
            t.ids_only_search()
            out = os.path.join(t._t.name, "P.md")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                                "--out", out, "--conditions", t.cond, "--session", t.s,
                                "--acquisition", "DDA"], capture_output=True, text=True,
                               timeout=120)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out, encoding="utf-8") as fh:
                with_default = nb.flat(fh.read())
            self.assertEqual(with_default, t.brief("--acquisition", "DDA"))
            self.assertIn("no differential expression", with_default)
        finally:
            t.tearDown()


class F9ReproduceIdsOnly(unittest.TestCase):
    """reproduce.sh of an identifications-only search: --adapt-only, --sbatch on a cluster, and
    the added FASTA in the bundle."""

    def setUp(self):
        import test_reproduce_replays_what_ran as rw
        self.rw = rw

        class B(rw.SageIdentificationsOnly):
            def runTest(self):
                pass
        self.b = B()
        self.b.setUp()
        self.addCleanup(self.b.tearDown)

    def bundle(self, *extra):
        rw = self.rw
        cfg, rec = self.b.sage_ids_only()
        return self.b.provenance("--workflow-manifest", self.b.manifest(), "--params", cfg,
                                 "--search-record", rec, "--engine", "sage",
                                 "--acquisition", "DDA", "--instrument", "timsTOF HT", *extra)

    def test_off_a_cluster_the_search_is_followed_by_adapt_only(self):
        import test_reproduce_sh_runs as sr
        sh, _md, _man = self.bundle()
        self.assertEqual(sr.continuation_breaks(sh), [])
        calls, _ = sr.run_steps_5_and_6(self.b.d, self.b.out, sh)
        self.assertEqual([c["prog"] for c in calls], ["run_search.py", "run_search.py"])
        self.assertNotIn("--adapt-only", calls[0]["argv"])
        self.assertNotIn("--sbatch", calls[0]["argv"])
        self.assertEqual(calls[1]["argv"], calls[0]["argv"] + ["--adapt-only"])

    def test_on_a_cluster_it_is_a_job_and_then_adapt_only(self):
        import test_reproduce_sh_runs as sr
        sh, _md, _man = self.bundle()
        calls, _ = sr.run_steps_5_and_6(self.b.d, self.b.out, sh, REPRO_SBATCH="1")
        self.assertEqual([c["prog"] for c in calls], ["run_search.py", "sbatch", "run_search.py"])
        self.assertEqual(sr.arg(calls[0]["argv"], "--sbatch"), "./search_job.sh")
        self.assertEqual(calls[1]["argv"], ["--wait", "./search_job.sh"])
        self.assertIn("--adapt-only", calls[2]["argv"])

    def test_the_added_fasta_travels_with_the_bundle(self):
        import hashlib
        added = os.path.join(self.b.d, "added.fasta")
        with open(added, "w") as fh:
            fh.write(">IN|1|BAIT\nMPEPTIDEKAAAR\n")
        with open(added, "rb") as fh:
            sha = hashlib.sha256(fh.read()).hexdigest()
        info = os.path.join(self.b.d, "search.fasta.meta.json")
        with open(info, "w") as fh:
            json.dump({"proteome": "UP000005640", "added_sequences": [
                {"file": added, "sha256": sha, "n_entries": 1}]}, fh)
        sh, _md, man = self.bundle("--fasta-info-file", info)
        self.assertIn("--add-fasta ./inputs/added_sequences/1_added.fasta", sh)
        self.assertNotIn(f"--add-fasta {added}", sh)
        with open(os.path.join(self.b.out, "inputs", "added_sequences", "1_added.fasta")) as fh:
            self.assertIn("MPEPTIDEKAAAR", fh.read())

    def test_a_changed_added_fasta_is_referenced_by_its_hash(self):
        added = os.path.join(self.b.d, "added.fasta")
        with open(added, "w") as fh:
            fh.write(">IN|1|BAIT\nMPEPTIDEKAAAR\n")
        info = os.path.join(self.b.d, "search.fasta.meta.json")
        with open(info, "w") as fh:
            json.dump({"proteome": "UP000005640", "added_sequences": [
                {"file": added, "sha256": "ab" * 32, "n_entries": 1}]}, fh)
        sh, _md, man = self.bundle("--fasta-info-file", info)
        self.assertIn(f"--add-fasta {added}", sh)
        self.assertIn("ab" * 32, sh + man)
        self.assertFalse(os.path.exists(os.path.join(self.b.out, "inputs", "added_sequences")))


class F10OutputFiles(unittest.TestCase):
    def test_a_png_is_not_said_to_be_embedded_unless_it_is_a_report_figure(self):
        import make_report as mr
        cat, desc = mr.describe("RUN1_RT_calibration_curves.png")
        self.assertNotIn("Publication-quality", desc)
        self.assertIn("not embedded", desc)
        self.assertIn("MSBooster", desc)
        self.assertIn("figures.json are embedded", desc)

    def test_fragpipe_tables_are_recognized(self):
        import make_report as mr
        for name in ("psm.tsv", "peptide.tsv", "protein.tsv"):
            cat, desc = mr.describe(name)
            self.assertEqual(cat, "Search output", name)
            self.assertNotIn("unrecognized", desc, name)


class F13NoDeDeliveryHasTheProteinTable(unittest.TestCase):
    def plan(self, files):
        import core_submission
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "search"))
            for n in files:
                with open(os.path.join(d, "search", n), "w") as fh:
                    fh.write("x\n")
            return {i["name"]: i for i in core_submission.plan_delivery(d)}

    def test_a_pooled_basic_search_delivers_its_protein_and_peptide_tables(self):
        items = self.plan(["protein.tsv", "peptide.tsv", "psm.tsv"])
        self.assertEqual(items["search/protein.tsv"]["status"], "OK")
        self.assertEqual(items["search/peptide.tsv"]["status"], "OK")
        self.assertNotIn("search/psm.tsv", items, "one row per spectrum, with file paths")

    def test_without_them_the_plan_says_there_is_none(self):
        item = self.plan([])["search/protein.tsv"]
        self.assertEqual(item["status"], "SKIPPED")
        self.assertFalse(item["required"])


class F14AgentsNoDe(unittest.TestCase):
    def agents(self, de):
        import session_docs
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d, de=de)
            os.makedirs(p["repro_dir"], exist_ok=True)
            with open(os.path.join(p["repro_dir"], "REPRODUCE.md"), "w") as fh:
                fh.write("# Reproduce\n")
            return session_docs.agents_md(session_docs.gather(p["session_dir"]))

    def test_a_no_de_session_is_not_told_about_a_de(self):
        text = self.agents(False)
        for word in ("Search + DE", "`de_provenance.json` / `methods.txt`", "`output/tables/`?"):
            self.assertNotIn(word, text)
        self.assertIn("The search: ", text)
        self.assertIn("`search_provenance.json` and the Methods (`output/methods.md`)", text)

    def test_a_de_session_keeps_its_wording(self):
        text = self.agents(True)
        self.assertIn("Search + DE: ", text)
        self.assertIn("`de_provenance.json` / `methods.txt`", text)


if __name__ == "__main__":
    unittest.main()
