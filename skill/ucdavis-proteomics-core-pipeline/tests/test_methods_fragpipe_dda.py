#!/usr/bin/env python3
"""The Methods describe a FragPipe DDA search from its own records (2.15, F5; rule 1).

What a client reads in the Database-search paragraph of a FragPipe DDA search must come from what
ran, never from what the skill would have run:
  * the engine choice: the bundle's client sentence (engine_choice.choose), quoted as written, and
    nothing when it is null (the data type's default);
  * FragPipe's workflow (the preset's template_choice) and each step FragPipe's log says FINISHED,
    with the version the log states -- a tool in the log's version list that did not run is never
    named or cited (Percolator/MSBooster in a Basic-Search run with PeptideProphet);
  * citations only for the steps that ran (Crossref-checked DOIs: run_search.CITE);
  * the decoys the search's job added (fragpipe_db.py: Philosopher, rev_, no contaminants);
  * isotope errors from fragger.params (Ruling 90: read from the run, IDs-only searches too);
  * ion-trap / unit-resolution tolerances as "adapted from" FragPipe's own workflow (the preset's
    tolerances.mode; make_presets.TOLERANCE_METHODS, the one wording);
  * the quantity: directLFQ with IonQuant's Intensity beside it and why (fragpipe_quantity), or
    IonQuant's MaxLFQ;
  * MSBooster's RT bandwidth as run (Ruling 106: fragpipe_params.msbooster_rerun_sentence, the one
    sentence, which REPRODUCE.md uses too), "not stated in the log" for a run without one;
  * no path and no licence word.
Fixtures: invented logs in FragPipe 24.0's line shapes (fixtures/fragpipe24_methods_log/), F1's
fragger.params of the two-run reference without its paths, FragPipe 24.0's LFQ-MBR workflow.
"""
import json
import os
import re
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
FIX = os.path.join(HERE, "fixtures")
LOGS = os.path.join(FIX, "fragpipe24_methods_log")
sys.path.insert(0, SCRIPTS)

import fragpipe_params as fpp  # noqa: E402
import make_methods as mm  # noqa: E402
import make_presets  # noqa: E402

LOG_NAME = "log_2026-11-05_00-27-23.txt"
FALLBACK = ("FragPipe, the pipeline's default search engine for timsTOF DDA-PASEF data, was not "
            "available where this search ran, so the data were searched with Sage.")
LICENCE_WORDS = re.compile(r"\b(licen[cs]ed?|academic|commercial|for-profit)\b", re.I)


def _json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh)
    return path


def directlfq_quantity():
    return {"quantity": "directLFQ", "n_experiments": 2, "status": "ok",
            "why_not_maxlfq": ("IonQuant computes MaxLFQ only from 3 experiments on (FragPipe "
                               "issues #1389, #1770); this search has 2. IonQuant's log: \"There "
                               "are only 2 experiments. Will not calculate MaxLFQ intensity.\""),
            "companion_column": "PG.IonQuant.Intensity", "run_map": True, "caveat": None,
            "directlfq": {"status": "ok", "versions": {"directlfq": "0.3.3", "pandas": "2.2.3"},
                          "guard": {"passed": True}}}


class FragPipeSearch(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def search(self, log="log_lfq_mbr.txt", quantity=None, preset=None, prov_extra=None,
               database=True):
        d = self.d
        shutil.copy(os.path.join(FIX, "fragpipe24_dda_preset", "c0_r1.fragger.params"),
                    os.path.join(d, "fragger.params"))
        shutil.copy(os.path.join(SKILL, "presets", "fragpipe", "LFQ-MBR.workflow"),
                    os.path.join(d, "fragpipe.workflow"))
        if log:
            shutil.copy(os.path.join(LOGS, log), os.path.join(d, LOG_NAME))
        wf = os.path.join(d, "preset.workflow")
        shutil.copy(os.path.join(SKILL, "presets", "fragpipe", "LFQ-MBR.workflow"), wf)
        if preset is not None:
            _json(wf + ".provenance.json", dict({"engine": "fragpipe"}, **preset))
        prov = {"engine": "fragpipe", "version": "24.0", "params_file": wf,
                "engine_selection": {"engine": "fragpipe", "how": "default", "sentence": None},
                "fragpipe_workflow": {"given": wf, "run": os.path.join(d, "fragpipe_input.workflow")}}
        if database:
            params = fpp.read_fragger_params(os.path.join(d, "fragger.params"))
            prov["fragpipe_database"] = {
                "path": os.path.join(d, "fragpipe_db", os.path.basename(params["database_name"])),
                "decoys_built_in_job": True, "sha256": "ab" * 32, "n_entries": 4806,
                "n_decoys": 2403, "decoy_tag": "rev_", "decoy_tag_source":
                "default: FragPipe's own rev_ (the workflow sets no database.decoy-tag)",
                "decoy_build": {"tool": "Philosopher", "version": "5.1.3-RC9",
                                "contaminants_added": False}}
        if quantity is not None:
            prov["fragpipe_quantity"] = quantity
        prov.update(prov_extra or {})
        sp = _json(os.path.join(d, "search_provenance.json"), prov)
        rec = mm.search_record(search_prov=sp)
        return rec, mm.search_paragraph(rec)

    # -- what ran, and its citations ----------------------------------------------------------
    def test_a_requested_directlfq_is_said_to_be_requested_not_two_experiments(self):
        import fragpipe_quant as fq
        q = directlfq_quantity()
        q.update(n_experiments=6, why_not_maxlfq=fq.REQUESTED_WORDS)
        rec, text = self.search(quantity=q)
        self.assertIn("directLFQ was requested instead of IonQuant's MaxLFQ", text)
        self.assertNotIn("only from 3 experiments", text)
        self.assertIn("directLFQ", text)

    def test_calibration_is_stated_as_the_workflows_setting_not_a_warning(self):
        """C2 proteomics M5: LFQ-MBR (the timsTOF DDA default) ships calibrate_mass=2; the Methods
        say what MSFragger did, never "[... calibration was on ... — confirm]"."""
        rec, text = self.search()
        self.assertNotIn("calibration was on", text)
        self.assertFalse(any("calibrat" in w for w in rec.get("warnings") or []))
        self.assertIn("calibrate_mass=2, as the search's workflow sets it", text)
        self.assertIn("starting values", text)
        self.assertEqual(fpp.calibration_sentence("0", "x"), "")
        self.assertNotIn("starting values", fpp.calibration_sentence("1", "x"))

    def test_the_workflow_and_each_step_that_ran_with_its_version(self):
        _rec, text = self.search(preset={"template_choice": {"workflow": "LFQ-MBR"}})
        self.assertIn("Within FragPipe's LFQ-MBR workflow, MSFragger 4.4.1 searched the spectra; "
                      "MSBooster 1.4.14 added rescoring features", text)
        self.assertIn("Percolator 3.7.1 rescored the peptide-spectrum matches", text)
        self.assertIn("ProteinProphet inferred the proteins", text)
        self.assertIn("Philosopher 5.1.3-RC9 filtered the identifications", text)
        self.assertIn("IonQuant 1.11.20 quantified the ions, with match-between-runs", text)
        # a step that exited non-zero is said, never described as run
        self.assertIn("non-zero exit for: PTMShepherd", text)

    def test_citations_only_for_what_ran(self):
        _rec, text = self.search(quantity=directlfq_quantity())
        cites = text[text.index("Software citations:"):]
        for doi in ("10.1038/nmeth.4256", "10.1038/s41467-023-40129-9", "10.1038/nmeth1113",
                    "10.1007/s13361-016-1460-7",
                    "10.1038/s41592-020-0912-y", "10.1016/j.mcpro.2021.100077",
                    "10.1016/j.mcpro.2023.100581"):
            self.assertIn(doi, cites)
        # Basic-Search with PeptideProphet: MSBooster/Percolator/IonQuant are in the log's version
        # list but never ran, so they are neither named nor cited
        shutil.rmtree(self.d)
        os.makedirs(self.d)
        _rec, text = self.search(log="log_basic_peptideprophet.txt")
        self.assertIn("PeptideProphet (Philosopher 5.1.3-RC9) validated", text)
        for absent in ("MSBooster", "Percolator", "IonQuant 1.11.20", "10.1038/nmeth1113",
                       "10.1007/s13361-016-1460-7",
                       "10.1038/s41467-023-40129-9", "10.1016/j.mcpro.2021.100077",
                       "directLFQ"):
            self.assertNotIn(absent, text, absent)
        self.assertIn("10.1038/s41592-020-0912-y", text)      # Philosopher ran

    def test_every_cited_doi_is_in_the_one_table(self):
        import run_search
        for key in ("philosopher", "msbooster", "percolator", "percolator3"):
            self.assertRegex(run_search.CITE[key], r"\[10\.\d{4,}/\S+\]$")
        for step, keys in mm.FRAGPIPE_STEP_CITE.items():
            for key in keys:
                self.assertIn(key, run_search.CITE, step)

    def test_no_log_is_said_not_guessed(self):
        _rec, text = self.search(log=None)
        self.assertIn("the steps that ran: " + mm.NOT_RECORDED, text)
        self.assertNotIn("Software citations", text)

    # -- the database, tolerances, quantity, bandwidth -----------------------------------------
    def test_decoys_built_in_the_job(self):
        _rec, text = self.search()
        self.assertIn("Decoys (reversed protein sequences, tagged rev_, FragPipe's default) were "
                      "added to the target database by Philosopher 5.1.3-RC9 within the search "
                      "job; the job added no contaminant entries of its own (any contaminants are "
                      "those already in the target database), and MSFragger searched the "
                      "resulting 4,806 entries (2,403 decoys; sha256 abababababab…).", text)

    def test_isotope_errors_from_the_run(self):
        rec, text = self.search()
        self.assertEqual(rec["isotope_error"]["source"], "fragger.params")
        self.assertRegex(text, r"with precursor isotope errors [0-9, ]+ or \d\.")

    def test_ion_trap_and_unit_resolution_are_adapted_from_fragpipes_workflows(self):
        _rec, text = self.search(preset={"tolerances": {"mode": "ion_trap"}})
        self.assertIn("These tolerances are " + make_presets.TOLERANCE_METHODS["ion_trap"]
                      + " (FragPipe ships no label-free workflow with an ion-trap fragment "
                        "tolerance).", text)
        _rec, text = self.search(preset={"tolerances": {"mode": "unit_resolution"}})
        self.assertIn("adapted from its Stellar-DDA workflow", text)
        _rec, text = self.search(preset={"tolerances": {"mode": "vendor"}})
        self.assertNotIn("These tolerances are", text)

    def test_directlfq_with_ionquant_intensity_beside_it_and_why(self):
        _rec, text = self.search(quantity=directlfq_quantity())
        self.assertIn("IonQuant computes MaxLFQ only from 3 experiments on (FragPipe issues "
                      "#1389, #1770); this search has 2. Protein quantities were computed with "
                      "directLFQ 0.3.3 from IonQuant's ion intensities (combined_ion.tsv), with "
                      "directLFQ's own between-sample normalisation, checked to have taken "
                      "effect.", text)
        self.assertIn("IonQuant's protein Intensity is reported beside it (column "
                      "PG.IonQuant.Intensity) but not used to compare samples", text)
        self.assertNotIn("IonQuant's log:", text)            # the record's evidence, not Methods

    def test_one_experiment_says_not_normalised_in_plain_words(self):
        """C2 bio M1: IonQuant's plain Intensity for one experiment carries the caveat a reader
        needs, and "experiment" is glossed, never bare."""
        import fragpipe_quant as fq
        q = fq.plan(["Protein", "Protein ID", "r1 Spectral Count", "r1 Intensity"])
        self.assertEqual(q["quantity"], fq.QUANTITY_INTENSITY)
        rec, text = self.search(quantity=dict(q, status="ok"))
        self.assertIn("not normalised between runs; compare runs with care", text)
        self.assertIn("one experiment (FragPipe's grouping of runs)", text)
        self.assertNotIn("one experiment:", text)

    def test_a_semi_specific_search_carries_the_caveat_everywhere(self):
        """C2 bio M3: Methods, the report's callout and the brief quote the one caveat
        (cleavage_rules.SEMI_SPECIFIC_CAVEAT) for num_enzyme_termini 1, and only then."""
        import analysis_prompt
        import cleavage_rules
        import make_analysis_html as mah
        rec, text = self.search()
        self.assertNotIn(cleavage_rules.SEMI_SPECIFIC_CAVEAT, text)          # fully specific
        report = {"input": os.path.join(self.d, "report.parquet")}
        self.assertIsNone(mah.semi_specific_note(report))
        fp = os.path.join(self.d, "fragger.params")
        with open(fp) as fh:
            body = fh.read()
        with open(fp, "w") as fh:
            fh.write(body.replace("num_enzyme_termini = 2", "num_enzyme_termini = 1"))
        rec = mm.search_record(search_prov=os.path.join(self.d, "search_provenance.json"))
        self.assertEqual(rec["cleavage"]["specificity"], "semi")
        text = mm.search_paragraph(rec)
        self.assertIn(cleavage_rules.SEMI_SPECIFIC_CAVEAT, text)
        self.assertIn("in-source fragments", mah.semi_specific_note(report)["text"])
        lines = []
        analysis_prompt.semi_specific_brief(lines.append, self.d)
        self.assertTrue(any(cleavage_rules.SEMI_SPECIFIC_CAVEAT in ln for ln in lines))

    def test_maxlfq_from_three_experiments(self):
        _rec, text = self.search(quantity={"quantity": "IonQuant MaxLFQ Intensity",
                                           "n_experiments": 6, "status": "ok"})
        self.assertIn("Protein quantities are IonQuant's MaxLFQ intensities (6 experiments).",
                      text)

    def test_msbooster_bandwidth_as_run(self):
        _rec, text = self.search()
        self.assertIn("A re-run of this FragPipe search can differ slightly from it, because "
                      "MSBooster chooses its retention-time smoothing (the LOESS bandwidth) for "
                      "each run itself; this search ran with 0.1 for one run; the log does not "
                      "state it for 1 of the 2 runs (chosen by MSBooster, as FragPipe's log "
                      "states; each run's value is in the search record, msbooster_rt_bandwidth in "
                      "search_provenance.json).", text)
        # Ruling 110.2: the Methods never name a run file
        self.assertNotIn("S1-a", text)
        self.assertNotIn("S2-b", text)

    def test_the_records_bandwidth_wins_over_the_log(self):
        bw = {"runs": [{"run": "S1-a", "bandwidth": 0.11, "source": "chosen by MSBooster"}]}
        _rec, text = self.search(prov_extra={fpp.MSBOOSTER_KEY: bw})
        self.assertIn("this search ran with 0.11 for one run (chosen by MSBooster", text)

    def test_the_parameter_table_rows(self):
        import subprocess
        self.search(preset={"template_choice": {"workflow": "LFQ-MBR"}},
                    quantity=directlfq_quantity())
        out = os.path.join(self.d, "methods.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"), "--raw",
                            os.path.join(self.d, "none.d"), "--instrument", "timsTOF HT",
                            "--acquisition", "DDA", "--search-prov",
                            os.path.join(self.d, "search_provenance.json"), "--out", out],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(out, encoding="utf-8") as fh:
            md = fh.read()
        for row in ("| FragPipe workflow | LFQ-MBR |", "| Precursor isotope errors | 0/1/2 |",
                    "| Decoys | added in the search job by Philosopher (tag rev_) |",
                    "| Protein quantity | directLFQ (IonQuant Intensity beside it: "
                    "PG.IonQuant.Intensity) |"):
            self.assertIn(row, md)

    def test_the_de_paragraph_points_at_the_search_paragraphs_quantity(self):
        # Ruling 110.4: one description of the quantity in the Methods (rule 3)
        rec, _ = self.search(quantity=directlfq_quantity())
        de = {"display_label": "directLFQ (FragPipe IonQuant ions) + limma",
              "rollup_method": "protein intensity = directLFQ 0.3.3 protein intensity (...), as "
                               "the search engine reported it"}
        text = mm.de_paragraph(de, rec)
        self.assertIn("Protein quantities: the search's own (directLFQ), as described under "
                      "Database search; the DE used them as reported (no re-rollup).", text)
        self.assertNotIn("protein intensity = directLFQ", text)
        # any other search keeps the DE record's own description
        self.assertIn("Protein quantities: protein intensity = directLFQ",
                      mm.de_paragraph(de, {"engine": "sage"}))

    # -- the choice -------------------------------------------------------------------------
    def test_the_selection_sentence_is_quoted_and_nothing_for_the_default(self):
        rec, text = self.search()
        self.assertIsNone(rec["selection_sentence"])
        self.assertNotIn("was not available where this search ran", text)
        sp = _json(os.path.join(self.d, "sage_prov.json"),
                   {"engine": "sage", "version": "0.14.7",
                    "engine_selection": {"engine": "sage", "default": "fragpipe",
                                         "how": "fallback", "sentence": FALLBACK}})
        text = mm.search_paragraph(mm.search_record(search_prov=sp))
        self.assertIn("Raw data were processed with Sage 0.14.7. " + FALLBACK, text)

    def test_client_text_has_no_path_or_licence_word(self):
        _rec, text = self.search(preset={"template_choice": {"workflow": "LFQ-MBR"}},
                                 quantity=directlfq_quantity())
        self.assertNotIn(self.d, text)
        self.assertNotRegex(text, r"(?<![\w.])/(data|opt|quobyte|home)/")
        self.assertIsNone(LICENCE_WORDS.search(text), LICENCE_WORDS.search(text))


class TheSentences(unittest.TestCase):
    def test_many_runs_are_summed_up(self):
        runs = [{"run": f"r{i}", "bandwidth": b, "source": "chosen by MSBooster"}
                for i, b in enumerate((0.05, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1))]
        runs.append({"run": "r9", "bandwidth": "not stated in the log",
                     "source": "not stated in the log"})
        s = fpp.msbooster_rerun_sentence({"runs": runs})
        self.assertIn("this search ran with 0.05–0.2 across the 7 runs; the log does not state "
                      "it for 1 of the 8 runs (chosen by MSBooster", s)
        self.assertNotIn("r9", s)
        two = fpp.msbooster_rerun_sentence({"runs": [
            {"run": "S1-a", "bandwidth": 0.1, "source": "chosen by MSBooster"},
            {"run": "S2-b", "bandwidth": 0.17, "source": "chosen by MSBooster"}]})
        self.assertIn("this search ran with 0.1–0.17 across the 2 runs", two)
        self.assertNotIn("S1-a", two)
        self.assertIsNone(fpp.msbooster_rerun_sentence({"runs": []}))
        self.assertIsNone(fpp.msbooster_rerun_sentence(None))

    def test_tools_run_reads_only_finished_steps(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        self.assertIsNone(fpp.tools_run(d))
        shutil.copy(os.path.join(LOGS, "log_lfq_mbr.txt"), os.path.join(d, LOG_NAME))
        t = fpp.tools_run(d)
        self.assertEqual(t["log"], LOG_NAME)
        self.assertEqual(t["ran"]["Percolator"], 2)
        self.assertEqual(t["failed"], {"PTMShepherd": 1})
        self.assertEqual(t["versions"]["MSFragger"], "4.4.1")
        self.assertNotIn("Skyline", t["versions"])                 # "N/A" left out
        self.assertNotIn("/", json.dumps(t).replace("N/A", ""))     # file names only


if __name__ == "__main__":
    unittest.main()
