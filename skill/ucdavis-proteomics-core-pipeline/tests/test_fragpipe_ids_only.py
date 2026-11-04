#!/usr/bin/env python3
"""
A FragPipe workflow that runs no quantification is an identifications-only search, recorded and
checked as one (2.14.2, L11, Ruling 73a).

FragPipe 24.0's shipped Nonspecific-peptidome.workflow has quantitation.run-label-free-quant=false:
it writes each experiment's psm.tsv / peptide.tsv / protein.tsv and never combined_protein.tsv.
A 38-run hydrolysate search with it finished ("ALL JOBS DONE IN 112.6 MINUTES") and its job was
still reported FAILED, on a check for a quantification table the workflow never makes. Now:
  * fragpipe_params.workflow_quantification() is the one reading of the workflow's switches;
  * the job's success check follows it (identification tables always, combined_protein.tsv only
    when label-free quantification ran), so a finished search is never reported as failed;
  * the search is recorded as identifications only, with the quantification "not applicable"
    record (sage_lfq_check.not_applicable, the one definition) -- the file run_de.R asks the gate
    about, so a DE is refused with a message that says how to quantify;
  * --adapt-only writes identifications.tsv (and fragpipe_identifications.json), no
    report.parquet; the Methods say identifications only;
  * a workflow WITH label-free quantification keeps today's path (pin).
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import fragpipe_params as fpp  # noqa: E402
import sage_lfq_check  # noqa: E402
from test_fragpipe_dda_stage import FP_WORKFLOW, StageHarness  # noqa: E402

with open(FP_WORKFLOW, newline="") as _fh:
    SHIPPED = _fh.read()                                   # LFQ off, as FragPipe ships it
LFQ_ON = SHIPPED.replace("quantitation.run-label-free-quant=false",
                         "quantitation.run-label-free-quant=true")
assert LFQ_ON != SHIPPED


class TheReading(unittest.TestCase):
    def test_the_shipped_workflow_runs_no_quantification(self):
        q = fpp.workflow_file_quantification(FP_WORKFLOW)
        self.assertIs(q["quantified"], False)
        self.assertIn("no label-free quantification", q["why"])
        self.assertEqual(q["switches"][fpp.LFQ_KEY], "false")
        self.assertEqual(fpp.expected_outputs(q), fpp.ID_TABLES)

    def test_a_workflow_with_label_free_quantification(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "lfq.workflow")
            with open(p, "w", newline="") as fh:
                fh.write(LFQ_ON)
            q = fpp.workflow_file_quantification(p)
        self.assertIs(q["quantified"], True)
        self.assertTrue(q["label_free"])
        self.assertEqual(fpp.expected_outputs(q), fpp.ID_TABLES + (fpp.LFQ_TABLE,))

    def test_ionquant_s_radio_button_is_not_a_switch(self):
        """ionquant.run-ionquant=true in every shipped workflow, TMT and DIA ones included"""
        q = fpp.workflow_quantification({fpp.LFQ_KEY: "false", "ionquant.run-ionquant": "true"})
        self.assertIs(q["quantified"], False)
        self.assertIsNone(fpp.workflow_quantification({"ionquant.run-ionquant": "true"})
                          ["quantified"])


class IdentificationsOnly(StageHarness):
    def prov(self):
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            return json.load(fh)

    def test_the_search_is_recorded_as_identifications_only(self):
        p = self.generate()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        prov = self.prov()
        self.assertIs(prov["quantified"], False)
        q = prov["quantification"]
        self.assertEqual(q["status"], "not_applicable")
        self.assertEqual(q["engine"], "fragpipe")
        self.assertIn("the workflow ran no label-free quantification", q["why"])
        with open(os.path.join(self.out, sage_lfq_check.RECORD)) as fh:
            self.assertEqual(json.load(fh)["status"], "not_applicable")

    def test_a_finished_search_without_combined_protein_is_a_success(self):
        self.assertEqual(self.generate().returncode, 0)
        r = self.run_job()                        # writes the ID tables, no combined_protein.tsv
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("the search finished", r.stdout + r.stderr)
        self.assertNotIn("S1-a.raw", os.listdir(self.stage()))

    def test_missing_identification_tables_are_a_failure(self):
        self.assertEqual(self.generate().returncode, 0)
        r = self.run_job(FAKE_FP_IDS="0")
        self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("not written: psm.tsv peptide.tsv protein.tsv", r.stdout + r.stderr)
        self.assertIn("S1-a.raw", os.listdir(self.stage()))      # kept for a re-run

    def test_the_de_gate_refuses_and_says_how_to_quantify(self):
        self.assertEqual(self.generate().returncode, 0)
        g = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sage_lfq_check.py"), "gate",
                            "--out", self.out], capture_output=True, text=True, timeout=60)
        j = json.loads(g.stdout)
        self.assertEqual((j["checked"], j["refused"]), (True, True))
        for words in ("FragPipe search ran no quantification", "no report.parquet", "no DE",
                      "identifications.tsv", f"{fpp.LFQ_KEY}=true", "LFQ-MBR", "quant.lfq"):
            self.assertIn(words, j["message"])

    def test_adapt_only_writes_the_identifications_and_no_report(self):
        self.assertEqual(self.generate().returncode, 0)
        self.assertEqual(self.run_job().returncode, 0)
        with open(os.path.join(self.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(SHIPPED)                     # FragPipe saves the workflow it ran
        p = self.generate("--adapt-only")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "report.parquet")))
        with open(os.path.join(self.out, "identifications.tsv")) as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh]
        self.assertEqual(rows[0][:2], ["run", "peptide"])
        self.assertEqual(sorted({r[0] for r in rows[1:]}), ["S1-a", "S2-b", "S3"])   # run names
        with open(os.path.join(self.out, "fragpipe_identifications.json")) as fh:
            rec = json.load(fh)
        self.assertIn("--sequential --prot 1", rec["kept_rule"])
        self.assertIs(rec["quantification"]["quantified"], False)

    def test_the_methods_say_identifications_only(self):
        import make_methods as mm
        self.assertEqual(self.generate().returncode, 0)
        with open(os.path.join(self.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(SHIPPED)
        text = mm.search_paragraph(mm.search_record(
            search_prov=os.path.join(self.out, "search_provenance.json")))
        # publication words, no config key (C2 biologist L4)
        self.assertIn("The search was for identifications only, with no quantification: "
                      "FragPipe's workflow ran no label-free (IonQuant), isobaric or DIA-NN "
                      "quantification.", text)
        self.assertNotIn("quantitation.run-label-free-quant", text)


class ReviewFollowUps(StageHarness):
    """Review F6, F7, F8, F13 (of 9abc8390), and the per-run identifications of a search whose
    Philosopher filter ran over all runs together (its tables at the workdir root)."""

    def test_no_proteinprophet_means_no_protein_table_is_expected(self):
        q = fpp.workflow_quantification({fpp.LFQ_KEY: "false",
                                         "protein-prophet.run-protein-prophet": "false",
                                         "phi-report.run-report": "true"})
        self.assertEqual(fpp.expected_outputs(q), ("psm.tsv", "peptide.tsv"))
        q = fpp.workflow_quantification({fpp.LFQ_KEY: "false", "phi-report.run-report": "false"})
        self.assertEqual(fpp.expected_outputs(q), ())

    def test_a_missing_lfq_key_is_not_identifications_only(self):
        q = fpp.workflow_quantification({"tmtintegrator.run-tmtintegrator": "false",
                                         "diann.run-dia-nn": "false"})
        self.assertIsNone(q["quantified"])

    def test_adapt_only_follows_what_ran_and_clears_a_stale_record(self):
        """generated with LFQ off, but the workflow FragPipe saved ran LFQ and wrote
        combined_protein.tsv: --adapt-only builds report.parquet and the IDs-only record goes"""
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("pyarrow")
        self.assertEqual(self.generate().returncode, 0)
        self.assertTrue(os.path.exists(os.path.join(self.out, sage_lfq_check.RECORD)))
        self.assertEqual(self.run_job(FAKE_FP_LFQ="1").returncode, 0)
        with open(os.path.join(self.out, "combined_protein.tsv"), "w") as fh:
            fh.write("Protein\tProtein ID\tS1_a Intensity\nsp|P02768|ALBU_HUMAN\tP02768\t100\n")
        with open(os.path.join(self.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(LFQ_ON)
        p = self.generate("--adapt-only")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue(os.path.exists(os.path.join(self.out, "report.parquet")))
        self.assertFalse(os.path.exists(os.path.join(self.out, sage_lfq_check.RECORD)))
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            self.assertEqual(json.load(fh)["quantification"]["status"], "quantified")

    def test_a_missing_key_is_decided_from_the_outputs(self):
        self.assertEqual(self.generate().returncode, 0)
        self.assertEqual(self.run_job().returncode, 0)
        with open(os.path.join(self.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(SHIPPED.replace("quantitation.run-label-free-quant=false\r\n", ""))
        p = self.generate("--adapt-only")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue(os.path.exists(os.path.join(self.out, "identifications.tsv")))
        with open(os.path.join(self.out, "fragpipe_identifications.json")) as fh:
            self.assertIn("no combined_protein.tsv",
                          json.load(fh)["quantification"]["why"])

    def test_the_methods_read_the_route_from_the_search_record(self):
        """F8: no fragpipe.workflow left in the folder, --params fragger.params"""
        import make_methods as mm
        self.assertEqual(self.generate().returncode, 0)
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)
        prov["params_file"] = os.path.join(self.out, "fragger.params")
        with open(os.path.join(self.out, "search_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        os.remove(os.path.join(self.d, "search.workflow"))
        text = mm.search_paragraph(mm.search_record(
            search_prov=os.path.join(self.out, "search_provenance.json")))
        self.assertIn("The search was for identifications only", text)

    def test_the_job_s_return_says_what_adapt_only_writes(self):
        p = self.generate()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("--adapt-only to write identifications.tsv", p.stdout)

    def test_identifications_per_run_when_philosopher_filtered_all_runs_together(self):
        """the tables at the workdir root, as FragPipe left them for a 38-run LFQ-off search"""
        import run_search
        out = self.out
        os.makedirs(out)
        with open(os.path.join(out, "psm.tsv"), "w") as fh:
            fh.write("Spectrum\tPeptide\tCharge\tProbability\tAssigned Modifications\t"
                     "Protein ID\tMapped Proteins\n"
                     "S1-a.00007.00007.2\tPEPTIDEK\t2\t0.99\t\tP02768\t\n"
                     "S2-b.00009.00009.3\tPEPTIDEK\t3\t0.98\t\tP02768\t\n"
                     "S2-b.00011.00011.2\tWEAKONEK\t2\t0.95\t\tP60709\t\n")
        with open(os.path.join(out, "peptide.tsv"), "w") as fh:
            fh.write("Peptide\tProbability\nPEPTIDEK\t0.99\n")      # WEAKONEK failed peptide FDR
        q = fpp.workflow_file_quantification(FP_WORKFLOW)
        path = run_search.adapt_fragpipe_ids(out, q)
        with open(path) as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh][1:]
        self.assertEqual(sorted((r[0], r[1]) for r in rows),
                         [("S1-a", "PEPTIDEK"), ("S2-b", "PEPTIDEK")])


class ZeroIntensities(StageHarness):
    """Follow-up 3 (H7): with several experiments, Philosopher's abacus writes a
    combined_protein.tsv for an IDs-only search too -- spectral counts, every intensity 0. The F7
    fallback (no label-free key) must not call that "quantified"."""

    def _adapt_with(self, intensity):
        self.assertEqual(self.generate().returncode, 0)
        self.assertEqual(self.run_job().returncode, 0)
        with open(os.path.join(self.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(SHIPPED.replace("quantitation.run-label-free-quant=false\r\n", ""))
        with open(os.path.join(self.out, "combined_protein.tsv"), "w") as fh:
            fh.write("Protein\tProtein ID\tS1-a Spectral Count\tS1-a Intensity\n"
                     f"sp|P02768|ALBU_HUMAN\tP02768\t3\t{intensity}\n")
        return self.generate("--adapt-only")

    def test_all_zero_intensities_are_identifications_only(self):
        p = self._adapt_with(0)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue(os.path.exists(os.path.join(self.out, "identifications.tsv")))
        self.assertFalse(os.path.exists(os.path.join(self.out, "report.parquet")))
        with open(os.path.join(self.out, "fragpipe_identifications.json")) as fh:
            self.assertIn("no intensities", json.load(fh)["quantification"]["why"])

    def test_real_intensities_are_quantified(self):
        """(pin) the F7 fallback as before"""
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("pyarrow")
        p = self._adapt_with(1500.5)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue(os.path.exists(os.path.join(self.out, "report.parquet")))


class TheFilterAsRun(StageHarness):
    """Follow-up 3 (H7): fragpipe_identifications.json quoted the workflow's filter line
    (--sequential --prot 1), the Methods the filter as run from FragPipe's log (... --razor).
    One reading now: fragpipe_params.filter_as_recorded."""

    def test_the_ids_record_and_the_methods_read_the_same_filter(self):
        import make_methods as mm
        self.assertEqual(self.generate().returncode, 0)
        self.assertEqual(self.run_job().returncode, 0)
        with open(os.path.join(self.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(SHIPPED)
        with open(os.path.join(self.out, "log_2026-11-03_12-00-00.txt"), "w") as fh:
            fh.write(f"PhilosopherFilter [Work dir: {self.out}]\n/opt/philosopher-v5.1.3 filter "
                     f"--sequential --prot 1 --tag rev_ --pepxml {self.out} --protxml "
                     f"{self.out}/combined.prot.xml --razor\n")
        self.assertEqual(self.generate("--adapt-only").returncode, 0)
        with open(os.path.join(self.out, "fragpipe_identifications.json")) as fh:
            rec = json.load(fh)
        self.assertIn("--sequential --prot 1 --razor", rec["kept_rule"])
        self.assertIn("FragPipe's log", rec["kept_rule"])
        self.assertIn("Philosopher (filter --sequential --prot 1 --razor)", rec["methods_sentence"])
        text = mm.search_paragraph(mm.search_record(
            search_prov=os.path.join(self.out, "search_provenance.json")))
        self.assertIn(rec["methods_sentence"], text)


class FollowUp2(unittest.TestCase):
    def test_without_peptide_tsv_the_rule_says_psm_level_only(self):
        """verify R2: a folder with psm.tsv and no peptide.tsv keeps every PSM's peptide; the
        record must not claim a peptide-level filter"""
        import run_search
        with tempfile.TemporaryDirectory() as out:
            with open(os.path.join(out, "psm.tsv"), "w") as fh:
                fh.write("Spectrum\tPeptide\tCharge\tProbability\tAssigned Modifications\t"
                         "Protein ID\tMapped Proteins\n"
                         "S1.00007.00007.2\tPEPTIDEK\t2\t0.99\t\tP02768\t\n")
            run_search.adapt_fragpipe_ids(out, fpp.workflow_file_quantification(FP_WORKFLOW))
            with open(os.path.join(out, "fragpipe_identifications.json")) as fh:
                rec = json.load(fh)
        self.assertIs(rec["peptide_filter_applied"], False)
        self.assertNotIn("whose peptide is in the same folder's peptide.tsv", rec["kept_rule"])
        self.assertIn("no peptide.tsv", rec["kept_rule"])
        self.assertIn("PSM level only", rec["kept_rule"])


class OneReportDefault(unittest.TestCase):
    """verify R3: an absent phi-report.run-report was "requires nothing" for the job's table
    check and "ran" for the Methods. One definition (fragpipe_params.report_runs): absent counts
    as true, as all 78 of FragPipe 24.0's shipped workflows set it."""

    def test_absent_counts_as_true_for_the_job(self):
        q = fpp.workflow_quantification({fpp.LFQ_KEY: "false"})
        self.assertTrue(fpp.report_runs({}))
        self.assertEqual(fpp.expected_outputs(q), ("psm.tsv", "peptide.tsv"))

    def test_the_job_and_the_methods_read_the_one_definition(self):
        from unittest import mock
        import make_methods as mm
        with tempfile.TemporaryDirectory() as d:
            import shutil
            shutil.copy(FP_WORKFLOW, os.path.join(d, "fragpipe.workflow"))
            sp = os.path.join(d, "search_provenance.json")
            with open(sp, "w") as fh:
                json.dump({"engine": "fragpipe", "version": "24.0",
                           "params_file": os.path.join(d, "fragpipe.workflow")}, fh)
            with mock.patch.object(fpp, "report_runs", return_value=False):
                q = fpp.workflow_file_quantification(FP_WORKFLOW)
                text = mm.search_paragraph(mm.search_record(search_prov=sp))
        self.assertEqual(fpp.expected_outputs(q), ())
        self.assertNotIn("Identifications were filtered with Philosopher", text)


class LabelFreeKeepsTodaysPath(StageHarness):
    """(pin) a workflow with label-free quantification: quantified, combined_protein.tsv
    expected, no gate record, the Methods do not say identifications only"""

    def test_lfq_on(self):
        import make_methods as mm
        p = self.generate(workflow=LFQ_ON)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertIs(prov["quantified"], True)
        self.assertFalse(os.path.exists(os.path.join(self.out, sage_lfq_check.RECORD)))
        self.assertNotEqual(self.run_job().returncode, 0)        # no combined_protein.tsv
        r = self.run_job(FAKE_FP_LFQ="1")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(self.out, "fragpipe.workflow"), "w", newline="") as fh:
            fh.write(LFQ_ON)
        text = mm.search_paragraph(mm.search_record(
            search_prov=os.path.join(self.out, "search_provenance.json")))
        self.assertNotIn("identifications only", text)

    def test_an_earlier_identifications_only_record_is_set_aside(self):
        self.assertEqual(self.generate().returncode, 0)          # LFQ off first
        self.assertTrue(os.path.exists(os.path.join(self.out, sage_lfq_check.RECORD)))
        self.assertEqual(self.generate(workflow=LFQ_ON).returncode, 0)
        self.assertFalse(os.path.exists(os.path.join(self.out, sage_lfq_check.RECORD)))
        self.assertTrue([f for f in os.listdir(self.out)
                         if f.startswith(sage_lfq_check.RECORD + ".stale")])


if __name__ == "__main__":
    unittest.main()
