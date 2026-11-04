#!/usr/bin/env python3
"""
The Methods state the search as it ran, for a non-specific DDA search too (2.14.2, L11 N1/N4/N5).

Three staff reports, on a 38-run hydrolysate peptidomics search searched without an enzyme:

  * Sage (cleave_at ""): the Methods said "the cleavage rule cleave_at='', restrict='P' [enzyme
    name not mapped -- confirm], up to 0 missed cleavages" and "Identification FDR: ____". The
    config had no `restrict` key and Sage's results.json recorded "restrict": null. Sage 0.14.7
    (crates/sage/src/database.rs, From<EnzymeBuilder> for EnzymeParameters) passes an absent
    `restrict` through as None -- no proline rule -- and DOCS.md's "default 'P'" holds only when the
    WHOLE enzyme block is absent (EnzymeBuilder::default). cleave_at "" is a non-specific digest
    (enzyme.rs: Enzyme::new("") -> None, every peptide of min_len..max_len, missed cleavages
    forced to 0); "$" is no digestion (each FASTA entry as one peptide).
  * FragPipe (the Nonspecific-peptidome workflow): "No fixed modifications were set. No variable
    modifications were searched." when five were; tolerances and FDR blank. fragger.params and
    the .workflow were never read.

Fixtures: tests/fixtures/sage_nonspecific/ (results.json in Sage 0.14.7's form, the config it ran
from) and tests/fixtures/fragpipe24_nonspecific_peptidome/ (FragPipe 24.0's shipped
Nonspecific-peptidome.workflow, verbatim with its CRLF line ends; the fragger.params FragPipe 24.0
writes for it, in that file's own format).
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
FIX = os.path.join(HERE, "fixtures")
sys.path.insert(0, SCRIPTS)

import make_methods as mm  # noqa: E402
import sage_lfq_check  # noqa: E402

SAGE_FIX = os.path.join(FIX, "sage_nonspecific")
FP_FIX = os.path.join(FIX, "fragpipe24_nonspecific_peptidome")

# A tryptic Sage config of the shape estimate_params.py writes. Its Methods text below is FROZEN:
# it is what make_methods.py wrote for it at release/skill-2.14.2 ea2d819a, before this lane.
TRYPTIC = {"database": {"enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 30,
                                   "cleave_at": "KR", "restrict": "P"},
                        "static_mods": {"C": 57.0215},
                        "variable_mods": {"[": [42.0106], "M": [15.9949]},
                        "max_variable_mods": 2},
           "precursor_tol": {"ppm": [-10.0, 10.0]}, "fragment_tol": {"ppm": [-10.0, 10.0]},
           "quant": {"lfq": True}}
TRYPTIC_TEXT = (
    "Raw data were processed with Sage 0.14.7. In silico digestion used Trypsin specificity, up "
    "to 2 missed cleavages, peptide length 7–30 residues. Fixed modifications: Carbamidomethyl "
    "(C). Variable modifications: Acetyl (Protein N-term); Oxidation (M) (at most 2 per "
    "peptide). Mass tolerances were precursor (MS1) ±10 ppm and fragment (MS2) ±10 ppm. "
    "Identification FDR: ____ [not recorded — confirm].")


def _json(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh)
    return path


def _load(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


class Tmp(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def sage(self, cfg, results=None, **prov):
        """A Sage search folder: the config, search_provenance.json, and Sage's results.json when
        given -> (record, paragraph)."""
        cfgp = _json(os.path.join(self.d, "sage_config.json"), cfg)
        if results is not None:
            _json(os.path.join(self.d, "results.json"), results)
        sp = _json(os.path.join(self.d, "search_provenance.json"),
                   {"engine": "sage", "version": "0.14.7", "params_file": cfgp, **prov})
        rec = mm.search_record(search_prov=sp)
        return rec, mm.search_paragraph(rec)


class SageEnzymeAsRun(Tmp):
    def test_a_non_specific_search_from_results_json(self):
        rec, text = self.sage(_load(os.path.join(SAGE_FIX, "sage_config.json")),
                              _load(os.path.join(SAGE_FIX, "results.json")))
        self.assertIn("no enzyme specificity (non-specific search)", text)
        self.assertNotIn("restrict", text)
        self.assertNotIn("enzyme name not mapped", text)
        self.assertNotIn("missed cleavage", text)          # Sage forces 0 on a non-specific digest
        self.assertIn("peptide length 7–30 residues", text)
        self.assertIn("results.json", rec["cleavage"]["source"])
        # the SDRF's cleavage agent (make_deposit.py), a verified PSI-MS term
        self.assertEqual((rec["cleavage"]["name"], rec["cleavage"]["ac"]),
                         ("unspecific cleavage", "MS:1001956"))

    def test_a_non_specific_search_from_its_config_alone(self):
        rec, text = self.sage(_load(os.path.join(SAGE_FIX, "sage_config.json")))
        self.assertIn("no enzyme specificity (non-specific search)", text)
        self.assertNotIn("restrict", text)
        self.assertNotIn("enzyme name not mapped", text)
        self.assertEqual(os.path.basename(rec["cleavage"]["source"]), "sage_config.json")

    def test_results_json_is_what_ran(self):
        """The config says trypsin, Sage's record of the run says no enzyme: the record wins."""
        rec, text = self.sage(TRYPTIC, _load(os.path.join(SAGE_FIX, "results.json")))
        self.assertIn("no enzyme specificity (non-specific search)", text)
        self.assertNotIn("Trypsin", text)

    def test_an_absent_restrict_is_no_proline_rule(self):
        cfg = json.loads(json.dumps(TRYPTIC))
        del cfg["database"]["enzyme"]["restrict"]
        rec, text = self.sage(cfg)
        self.assertIn("Trypsin/P specificity", text)
        self.assertEqual(rec["cleavage"]["rule"], "cleave_at='KR'")
        self.assertEqual(rec["cleavage"]["ac"], "MS:1001313")

    def test_restrict_null_in_results_json_is_no_proline_rule(self):
        res = _load(os.path.join(SAGE_FIX, "results.json"))
        res["database"]["enzyme"].update(cleave_at="KR", missed_cleavages=2)
        rec, text = self.sage(TRYPTIC, res)
        self.assertIn("Trypsin/P specificity", text)
        self.assertNotIn("restrict='P'", rec["cleavage"]["rule"])

    def test_dollar_is_no_digestion(self):
        cfg = json.loads(json.dumps(TRYPTIC))
        cfg["database"]["enzyme"] = {"cleave_at": "$", "min_len": 5, "max_len": 50}
        rec, text = self.sage(cfg)
        self.assertIn("no in silico digestion (each database entry searched as one peptide)",
                      text)
        self.assertEqual((rec["cleavage"]["name"], rec["cleavage"]["ac"]),
                         ("no cleavage", "MS:1001955"))
        self.assertNotIn("enzyme name not mapped", text)

    def test_a_semi_specific_search_says_so(self):
        cfg = json.loads(json.dumps(TRYPTIC))
        cfg["database"]["enzyme"]["semi_enzymatic"] = True
        rec, text = self.sage(cfg)
        self.assertIn("semi-specific Trypsin cleavage", text)
        self.assertNotIn("Trypsin specificity", text)

    def test_a_tryptic_search_s_methods_are_unchanged(self):
        """(pin) the frozen text, from the config alone and with Sage's own record beside it"""
        _, text = self.sage(TRYPTIC)
        self.assertEqual(text, TRYPTIC_TEXT)
        res = _load(os.path.join(SAGE_FIX, "results.json"))
        res["database"].update(TRYPTIC["database"])            # the same search, as Sage ran it
        res["database"]["enzyme"] = dict(TRYPTIC["database"]["enzyme"], c_terminal=None,
                                         semi_enzymatic=None)
        res.update(precursor_tol=TRYPTIC["precursor_tol"], fragment_tol=TRYPTIC["fragment_tol"])
        _, text = self.sage(TRYPTIC, res)
        self.assertEqual(text, TRYPTIC_TEXT)

    def test_mods_and_tolerances_are_what_ran_too(self):
        """Review F16: the config given says C fixed and 10 ppm; Sage's own record of the search
        (another config ran) says no fixed mods and 20 ppm -- the Methods say what ran"""
        res = _load(os.path.join(SAGE_FIX, "results.json"))
        res["database"]["enzyme"] = dict(TRYPTIC["database"]["enzyme"])
        res.update(precursor_tol={"ppm": [-20.0, 20.0]}, fragment_tol={"ppm": [-20.0, 20.0]})
        rec, text = self.sage(TRYPTIC, res)
        self.assertIn("No fixed modifications were set.", text)
        self.assertIn("precursor (MS1) ±20 ppm", text)
        self.assertIn("results.json", rec["ms1_tol"]["source"])
        self.assertTrue(all("results.json" in m["source"] for m in rec["mods"]))


class SageLfqFdr(Tmp):
    """The FDR sentence of a quantifying Sage search: the skill's own rule (sage_lfq_check's
    PSM_Q and LFQ_Q_MAX, the one definition adapt_sage keeps rows with), never a blank."""
    ADAPT = {"q_value_max": 0.05, "report": "/data/search/report.parquet",
             "kept_rule": "is_decoy false and stored q_value <= 0.05 (the stored q is shared "
                          "with the precursor's decoy, so this is a conservative subset of "
                          "Sage's own count)"}

    def test_a_record_from_before_the_sentence_is_described_by_its_own_threshold(self):
        _, text = self.sage(TRYPTIC, sage_adapt=self.ADAPT)          # a 2.14.0 sage_adapt
        self.assertIn(sage_lfq_check.methods_sentence(0.05), text)
        self.assertNotIn("Identification FDR: ____", text)

    def test_the_sentence_is_the_one_rule(self):
        s = sage_lfq_check.methods_sentence(sage_lfq_check.LFQ_Q_MAX)
        self.assertIn(f"{sage_lfq_check.PSM_Q * 100:g}%", s)
        self.assertIn(f"{sage_lfq_check.LFQ_Q_MAX:g}", s)
        self.assertIn("no protein-level FDR", s)
        self.assertIn("over all runs together", s)
        self.assertIn("0.1", sage_lfq_check.methods_sentence(0.1))   # the record's own value

    def test_a_record_s_own_sentence_is_quoted(self):
        _, text = self.sage(TRYPTIC, sage_adapt=dict(self.ADAPT, methods_sentence="AS RECORDED."))
        self.assertIn("AS RECORDED.", text)

    def test_a_refused_search_states_no_rule(self):
        _, text = self.sage(TRYPTIC, sage_adapt={"q_value_max": 0.05, "report": None,
                                                 "status": "refused"})
        self.assertIn("Identification FDR: ____", text)

    def test_a_de_that_states_it_is_not_repeated(self):
        cfgp = _json(os.path.join(self.d, "sage_config.json"), TRYPTIC)
        sp = _json(os.path.join(self.d, "search_provenance.json"),
                   {"engine": "sage", "version": "0.14.7", "params_file": cfgp,
                    "sage_adapt": self.ADAPT})
        text = mm.search_paragraph(mm.search_record(search_prov=sp),
                                   de_prov={"q_columns": ["Q.Value"]})
        self.assertNotIn(sage_lfq_check.methods_sentence(0.05), text)

    def test_adapt_sage_records_the_sentence(self):
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("pyarrow is needed for the Sage fixtures")
        from test_sage_adapt import write_lfq
        import run_search
        write_lfq(os.path.join(self.d, "lfq.parquet"))
        run_search.adapt_sage(self.d)
        rec = _load(os.path.join(self.d, "sage_adapt.json"))
        self.assertEqual(rec["methods_sentence"],
                         sage_lfq_check.methods_sentence(sage_lfq_check.LFQ_Q_MAX))

    def test_the_identifications_route_is_unchanged(self):
        """(pin) sage_identifications.methods_sentence is quoted as before"""
        _, text = self.sage(dict(TRYPTIC, quant={"lfq": False}),
                            sage_identifications={"methods_sentence": "IDS RULE."})
        self.assertIn("IDS RULE.", text)


class FragPipeWorkflowSearch(Tmp):
    """A FragPipe search folder as FragPipe leaves it: fragger.params (the MSFragger parameters it
    wrote) and fragpipe.workflow (the workflow it saved) beside search_provenance.json."""

    def search(self, params=True, workflow=True, edit=None, given="fragger.params"):
        if params:
            shutil.copy(os.path.join(FP_FIX, "fragger.params"),
                        os.path.join(self.d, "fragger.params"))
            if edit:
                p = os.path.join(self.d, "fragger.params")
                with open(p, encoding="utf-8") as fh:
                    text = fh.read()
                with open(p, "w", encoding="utf-8") as fh:
                    fh.write(edit(text))
        wf = os.path.join(self.d, "Nonspecific-peptidome.workflow")
        shutil.copy(os.path.join(FP_FIX, "Nonspecific-peptidome.workflow"), wf)
        if workflow:
            shutil.copy(wf, os.path.join(self.d, "fragpipe.workflow"))
        sp = _json(os.path.join(self.d, "search_provenance.json"),
                   {"engine": "fragpipe", "version": "24.0",
                    "params_file": os.path.join(self.d, given) if given == "fragger.params"
                    else wf})
        rec = mm.search_record(search_prov=sp)
        return rec, mm.search_paragraph(rec)

    def test_the_nonspecific_peptidome_search_is_described_as_run(self):
        rec, text = self.search()
        self.assertIn("no enzyme specificity (non-specific search)", text)
        self.assertNotIn("missed cleavage", text)
        self.assertIn("peptide length 7–65 residues", text)
        self.assertIn("Fixed modifications: Carbamidomethyl (C).", text)
        self.assertIn("Variable modifications: Oxidation (M); Acetyl (Protein N-term); "
                      "Gln->pyro-Glu (Any N-term Q);", text)
        self.assertIn("(Any N-term E)", text)               # -18.0106 nE
        self.assertIn("(Any C-term)", text)                 # -0.98402 c^
        self.assertIn("(at most 3 per peptide)", text)
        self.assertEqual(len([m for m in rec["mods"] if m["type"] == "variable"]), 5)
        self.assertIn("Mass tolerances were precursor (MS1) ±20 ppm and fragment (MS2) ±25 ppm.",
                      text)
        self.assertIn("N-terminal methionine excision was enabled.", text)
        self.assertIn("Philosopher (filter --sequential --prot 1)", text)
        self.assertIn("no protein-level FDR filter (--prot 1)", text)
        self.assertIn("PSM 1%, peptide ion 1% and peptide 1% FDR (Philosopher's defaults; not "
                      "set in the filter)", text)
        for false in ("No fixed modifications were set", "No variable modifications",
                      "not parsed here", "Identification FDR: ____", "enzyme name not mapped"):
            self.assertNotIn(false, text)
        self.assertEqual((rec["cleavage"]["name"], rec["cleavage"]["ac"]),
                         ("unspecific cleavage", "MS:1001956"))

    def test_the_workflow_given_finds_fragger_params_beside_the_record(self):
        rec, text = self.search(given="workflow")
        self.assertIn("Fixed modifications: Carbamidomethyl (C).", text)
        self.assertIn("fragger.params", rec["cleavage"]["source"])

    def test_no_fixed_modifications_only_when_the_params_say_none(self):
        rec, text = self.search(edit=lambda t: t.replace("add_C_cysteine = 57.02146",
                                                         "add_C_cysteine = 0.0"))
        self.assertIn("No fixed modifications were set.", text)

    def test_a_disabled_variable_mod_is_not_searched(self):
        rec, _ = self.search()
        self.assertFalse([m for m in rec["mods"] if m["mass"] and abs(m["mass"] - 4.025107) < 1e-3])

    def test_without_fragger_params_nothing_is_claimed(self):
        rec, text = self.search(params=False, given="workflow")
        self.assertIn("modifications: ____ [not recorded — confirm]", text)
        self.assertNotIn("No fixed modifications were set", text)
        self.assertNotIn("No variable modifications", text)
        self.assertIn("fragger.params", text)                # the warning names what is missing
        self.assertIn("Philosopher (filter --sequential --prot 1)", text)   # from the workflow

    def test_without_any_workflow_the_fdr_is_not_recorded(self):
        rec, text = self.search(workflow=False)
        self.assertIn("Identification FDR: ____ [not recorded — confirm]", text)

    def test_a_tryptic_lfq_filter(self):
        rec, text = self.search(edit=lambda t: (
            t.replace("search_enzyme_name_1 = nonspecific", "search_enzyme_name_1 = trypsin")
             .replace("search_enzyme_cut_1 = @", "search_enzyme_cut_1 = KR")
             .replace("search_enzyme_nocut_1 = ", "search_enzyme_nocut_1 = P")
             .replace("num_enzyme_termini = 0", "num_enzyme_termini = 2")))
        self.assertIn("Trypsin specificity, up to 2 missed cleavages", text)
        semi = self.search(edit=lambda t: (
            t.replace("search_enzyme_cut_1 = @", "search_enzyme_cut_1 = KR")
             .replace("search_enzyme_name_1 = nonspecific", "search_enzyme_name_1 = stricttrypsin")
             .replace("num_enzyme_termini = 0", "num_enzyme_termini = 1")))[1]
        self.assertIn("semi-specific Trypsin/P cleavage", semi)

    def test_philosopher_filter_words(self):
        import fragpipe_params as fp
        s = fp.philosopher_methods_sentence("--sequential --prot 0.01 --picked --razor")
        self.assertIn("protein 1% FDR (picked protein FDR; razor peptides for protein scoring)", s)
        self.assertIn("PSM 1%, peptide ion 1% and peptide 1% FDR (Philosopher's defaults", s)
        s = fp.philosopher_methods_sentence("--psm 0.02 --prot 0.01")
        self.assertIn("PSM 2%", s)
        self.assertIn("peptide ion 1% and peptide 1% FDR (Philosopher's defaults", s)


# FragPipe's console log as FragPipe 24.0 writes it (log_<date>.txt in the workdir), cut to the
# Philosopher filter lines: one filter over the whole workdir (a 38-run search with LFQ off), with
# the options FragPipe added (--tag, --razor) beside the workflow's line
ONE_FILTER_LOG = """PhilosopherFilter [Work dir: {wd}]
/opt/fragpipe-24.0/tools/Philosopher/philosopher-v5.1.3-RC9 filter --sequential --prot 1 --tag rev_ --pepxml {wd} --protxml {wd}/combined.prot.xml --razor
time="21:36:08" level=info msg="Executing Filter  v5.1.3"
Process 'PhilosopherFilter' finished, exit code: 0
"""
PER_EXPERIMENT_LOG = "".join(
    f"PhilosopherFilter [Work dir: {{wd}}/{e}]\n/opt/fragpipe-24.0/tools/Philosopher/"
    f"philosopher-v5.1.3-RC9 filter --sequential --prot 0.01 --picked --tag rev_ --pepxml "
    f"{{wd}}/{e} --protxml {{wd}}/combined.prot.xml --razor\n" for e in ("S1", "S2"))


class FragPipePrecision(Tmp):
    """Review F3, F4, F9 (of 9abc8390)."""
    search = FragPipeWorkflowSearch.search              # the helper, not the tests

    def test_mass_offset_and_labile_settings_are_stated_and_tagged(self):
        rec, text = self.search(edit=lambda t: (
            t.replace("mass_offsets = 0.0", "mass_offsets = 0/79.96633/-97.9769")
             .replace("fragment_ion_series = b,y", "labile_search_mode = nglycan\n"
                      "mass_diff_to_variable_mod = 1\nrestrict_deltamass_to = STY\n"
                      "fragment_ion_series = b,y")))
        for words in ("mass offsets of 0, 79.96633 and -97.9769 Da", "labile search mode nglycan",
                      "mass_diff_to_variable_mod=1", "restricted to STY", "— confirm]"):
            self.assertIn(words, text)

    def test_a_plain_search_gets_no_such_sentence(self):
        """(pin) mass_offsets 0.0 and no labile mode: nothing added"""
        rec, text = self.search()
        self.assertNotIn("mass offset", text)
        self.assertNotIn("labile", text)

    def test_only_residue_and_terminus_add_keys_are_fixed_mods(self):
        rec, text = self.search(edit=lambda t: t.replace(
            "write_mzbin_all = 0", "write_mzbin_all = 0\nadd_topN_complementary = 3"))
        self.assertIn("Fixed modifications: Carbamidomethyl (C).", text)
        self.assertNotIn("topN", text)

    def write_log(self, text):
        with open(os.path.join(self.d, "log_2026-11-03_12-00-00.txt"), "w") as fh:
            fh.write(text.format(wd=self.d))

    def test_the_filter_as_run_over_all_runs_together(self):
        self.write_log(ONE_FILTER_LOG)
        rec, text = self.search()
        self.assertIn("Philosopher (filter --sequential --prot 1 --razor)", text)
        self.assertIn("over all runs together", text)
        self.assertNotIn("each run separately", text)

    def test_the_filter_as_run_per_experiment(self):
        self.write_log(PER_EXPERIMENT_LOG)
        rec, text = self.search()
        self.assertIn("--picked --razor", text)
        self.assertIn("each FragPipe experiment separately (2 experiments)", text)
        # Follow-up 3 (H7 bug 2): each experiment's filter read ONE combined.prot.xml, so the
        # protein-level FDR was over all runs together -- only PSM, ion and peptide were per run
        self.assertIn("the PSM, peptide-ion and peptide rates hold within each experiment", text)
        self.assertIn("the protein-level 1% FDR was estimated over all runs together", text)
        self.assertNotIn("so these rates hold within each", text)

    def test_an_apostrophe_in_a_logged_path_is_read(self):
        """Follow-up 2 (verify R1): FragPipe logs its argv unquoted, so a work dir like
        /data/it's/ made shlex.split raise and the Methods crash"""
        self.d = os.path.join(self.d, "it's here")
        os.makedirs(self.d)
        self.write_log(ONE_FILTER_LOG)
        rec, text = self.search()
        self.assertIn("Philosopher (filter --sequential --prot 1 --razor)", text)
        self.assertIn("over all runs together", text)

    def test_a_log_line_that_cannot_be_read_falls_back_tagged(self):
        from unittest import mock
        import fragpipe_params
        self.write_log(ONE_FILTER_LOG)
        with mock.patch.object(fragpipe_params, "filter_as_run",
                               side_effect=ValueError("No closing quotation")):
            rec, text = self.search()
        self.assertIn("Philosopher (filter --sequential --prot 1)", text)     # the workflow's
        self.assertIn("could not be read (ValueError: No closing quotation)", text)
        self.assertIn("— confirm", text)

    def test_the_database_searched_is_named(self):
        """Follow-up 3 (H7 bug 1): fragger.params database_name is the database MSFragger read"""
        rec, text = self.search()
        self.assertIn("MSFragger searched 2026-01-01-decoys-search.fasta.fas", text)
        self.assertIn("fragger.params database_name", rec["database"]["source"])

    def test_without_the_log_the_scope_is_tagged(self):
        rec, text = self.search()
        self.assertIn("Philosopher (filter --sequential --prot 1)", text)      # the workflow's
        self.assertIn("FragPipe may add options when it runs it", text)
        self.assertIn("per run or over all runs together", text)


class UnparsedEngineClaimsNothing(Tmp):
    def test_an_engine_whose_params_are_not_read_never_says_no_modifications(self):
        cfg = os.path.join(self.d, "alphadia.yaml")
        with open(cfg, "w") as fh:
            fh.write("search:\n  target_num_candidates: 2\n")
        sp = _json(os.path.join(self.d, "search_provenance.json"),
                   {"engine": "alphadia", "version": "1.10", "params_file": cfg})
        text = mm.search_paragraph(mm.search_record(search_prov=sp))
        self.assertNotIn("No fixed modifications were set", text)
        self.assertNotIn("No variable modifications", text)
        self.assertIn("modifications: ____ [not recorded — confirm]", text)


class MsParagraphReadsTheRecordedMode(Tmp):
    """(pin -- fixed in this release by 237d681e, S7) the session's input/acquisition.json is
    found from the search record alone, as `make_methods.py --params ... --search-prov
    output/search_x/search_provenance.json` was run for the Sage report."""

    def test_acquisition_json_from_the_search_record(self):
        sess = os.path.join(self.d, "session")
        os.makedirs(os.path.join(sess, "output", "search_lfq"))
        os.makedirs(os.path.join(sess, "input"))
        _json(os.path.join(sess, "input", "acquisition.json"),
              {"overall": "DDA", "instrument": "Orbitrap Exploris 480"})
        sp = _json(os.path.join(sess, "output", "search_lfq", "search_provenance.json"),
                   {"engine": "sage"})

        class A:
            acquisition = None
            workflow_manifest = None
            out = None
            de_dir = None
            params = None
            search_prov = sp
            raw = []
        self.assertEqual(mm.acquisition_from_records(A())[0], "DDA")


if __name__ == "__main__":
    unittest.main()
