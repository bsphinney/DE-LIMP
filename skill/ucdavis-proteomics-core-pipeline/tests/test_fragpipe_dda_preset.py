#!/usr/bin/env python3
"""
The FragPipe DDA preset (2.15 lane F1; review REVIEW_part3 F1, Rulings 86 and 90).

make_presets.py --engine fragpipe --acquisition DDA used to refuse ("these presets are the DIA
routes"), so a FragPipe DDA search ran a hand-edited workflow (H5) or none. Now it writes one from
FragPipe 24.0's own shipped workflows, chosen by the digest and by whether amounts are wanted:
LFQ-MBR (an enzyme rule, amounts), Basic-Search (--ids-only), Nonspecific-peptidome or -HLA
(--enzyme none; label-free quantification switched on only when amounts are wanted). What these
tests pin:
  * the vendor files are FragPipe 24.0's, byte for byte (presets/fragpipe/SOURCES.json);
  * H5's golden workflow: the skill's preset for H5's chymotrypsin-like rule is H5's hand-made
    workflow (h5_golden_lfqmbr.workflow = H5's fragpipe_dda_lfq.workflow, sha256 3f6469fc...,
    without its database.db-path and workflow.threads lines) MINUS its enzyme dropdown line: C0
    (2026-11-04) found that FragPipe's dropdown overwrites the rule (dropdown `custom` ran as
    argc, `nonspecific` ran a non-specific search), and that the rule reaches fragger.params only
    with no dropdown key. That dropdown-free file is C0's wf/v1a.workflow (cmp-identical without
    its db-path), which FragPipe 24.0 loaded and saved (c0_v1a_dryrun.fragpipe.workflow, a
    test job on HIVE) and ran for real (c0_r1.fragger.params, R1): name nonspecific, cut FLMWY, no-cut P,
    sense C, num_enzyme_termini 2;
  * the record's cleavage rule never quotes the enzyme NAME (FragPipe relabels it from the
    dropdown: H5's real fragger.params says argc);
  * template choice: IDs-only, non-specific, LFQ; no database.db-path (run_search.py writes it at
    search time); the vendor data-type keys are kept;
  * ion-trap MS2 = FragPipe's own TMT10-MS3 setting (0.6 Da, +/-20 ppm, isotope_error 0/1/2/3),
    unit resolution (Thermo Stellar) = its Stellar-DDA setting (0.3 Da, +/-1.5 Da);
  * --enzyme takes fetch_fasta.py's words; the name -> rule map is cleavage_rules', checked
    against FragPipe 24.0's enzyme list (enzyme-defs-supported.24.0.txt, sha256 c1c58a79...);
  * one DDA flag set for estimate_params.py and make_presets.py, with the shared constants;
  * Sage configs and the DIA presets are byte-identical to 2f9b229c for the same arguments;
  * fetch_workflows.py (nothing called it) is gone.
Hermetic: subprocesses on temp dirs; no network, no FragPipe.
"""
import hashlib
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
PRESETS = os.path.join(SKILL, "presets", "fragpipe")
FIX = os.path.join(HERE, "fixtures", "fragpipe24_dda_preset")
sys.path.insert(0, SCRIPTS)

DROPDOWN = "msfragger.misc.fragger.enzyme-dropdown-1"
TIMS = "timsTOF HT"


def props(text):
    out = {}
    for ln in text.splitlines():
        if ln and not ln.startswith("#") and "=" in ln:
            k, v = ln.split("=", 1)
            out[k.strip()] = v
    return out


def preset(d, *args, acquisition="DDA", instrument=TIMS, name="fp.workflow"):
    out = os.path.join(d, name)
    cmd = [sys.executable, os.path.join(SCRIPTS, "make_presets.py"), "--engine", "fragpipe",
           "--acquisition", acquisition, "--out", out, *args]
    if instrument is not None:
        cmd += ["--instrument", instrument]
    p = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
    text = prov = None
    if p.returncode == 0:
        with open(out) as fh:
            text = fh.read()
        with open(out + ".provenance.json") as fh:
            prov = json.load(fh)
    return p, out, text, prov


def ok(test, res):
    p, out, text, prov = res
    test.assertEqual(p.returncode, 0, p.stderr)
    return out, text, prov


def var_rows(text):
    return [r.strip() for r in props(text)["msfragger.table.var-mods"].split(";")]


class VendorFiles(unittest.TestCase):
    def test_each_preset_is_fragpipe_24s_own_file(self):
        with open(os.path.join(PRESETS, "SOURCES.json")) as fh:
            src = json.load(fh)
        want = {"DIA_SpecLib_Quant.workflow", "DIA_SpecLib_Quant_diaPASEF.workflow",
                "LFQ-MBR.workflow", "Basic-Search.workflow", "Nonspecific-peptidome.workflow",
                "Nonspecific-HLA.workflow"}
        self.assertEqual(set(src["presets"]), want)
        self.assertEqual({f for f in os.listdir(PRESETS) if f.endswith(".workflow")}, want)
        for name, rec in src["presets"].items():
            with self.subTest(preset=name):
                with open(os.path.join(PRESETS, name), "rb") as fh:
                    data = fh.read()
                self.assertEqual(hashlib.sha256(data).hexdigest(), rec["sha256"])
                self.assertNotIn(b"\r", data)

    def test_the_adapted_values_name_their_vendor_files(self):
        with open(os.path.join(PRESETS, "SOURCES.json")) as fh:
            src = json.load(fh)["values_adapted_from"]
        for name in ("TMT10-MS3.workflow", "Stellar-DDA.workflow",
                     "enzyme-defs-supported.txt"):
            fx = os.path.join(FIX, name if name.endswith(".workflow")
                              else "enzyme-defs-supported.24.0.txt")
            with open(fx, "rb") as fh:
                self.assertEqual(hashlib.sha256(fh.read()).hexdigest(), src[name]["sha256"])


class H5Golden(unittest.TestCase):
    """H5's chymotrypsin-like rule (cut FLMWY, not before P), peptides 6-30, Met-ox off, at most
    2 variable modifications: what H5 wrote by hand, as the skill's preset."""
    ARGS = ("--var-mods", "none", "--overrides", json.dumps({
        "msfragger.search_enzyme_cut_1": "FLMWY", "msfragger.search_enzyme_nocut_1": "P",
        "msfragger.digest_min_length": 6, "msfragger.digest_max_length": 30}))

    def test_the_preset_is_h5s_workflow_without_the_dropdown(self):
        with open(os.path.join(FIX, "h5_golden_lfqmbr.workflow")) as fh:
            golden = fh.read()
        self.assertIn(f"{DROPDOWN}=custom\n", golden)
        want = "".join(ln for ln in golden.splitlines(keepends=True)
                       if not ln.startswith(DROPDOWN + "="))
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, *self.ARGS))
        self.assertEqual(text, want)
        self.assertEqual(prov["template"], "fragpipe/LFQ-MBR.workflow")
        self.assertIn(DROPDOWN, prov["keys_deleted"])
        self.assertNotIn("database.db-path", text)

    def test_fragpipe_ran_that_rule(self):
        """C0's dry run saved the workflow FragPipe 24.0 loaded from that file, and R1's
        fragger.params is what MSFragger ran: the rule, under the name nonspecific, 2 termini."""
        import fragpipe_params as fpp
        saved, _ = fpp.read_workflow(os.path.join(FIX, "c0_v1a_dryrun.fragpipe.workflow"))
        ran = fpp.read_fragger_params(os.path.join(FIX, "c0_r1.fragger.params"))
        with tempfile.TemporaryDirectory() as d:
            _, text, _ = ok(self, preset(d, *self.ARGS))
        mine = props(text)
        for key in ("search_enzyme_name_1", "search_enzyme_cut_1", "search_enzyme_nocut_1",
                    "search_enzyme_sense_1", "num_enzyme_termini", "digest_min_length",
                    "digest_max_length", "max_variable_mods_per_peptide"):
            with self.subTest(key=key):
                self.assertEqual(mine["msfragger." + key], saved["msfragger." + key])
                self.assertEqual(mine["msfragger." + key], ran[key])
        self.assertEqual((ran["search_enzyme_name_1"], ran["search_enzyme_cut_1"],
                          ran["search_enzyme_nocut_1"], ran["num_enzyme_termini"]),
                         ("nonspecific", "FLMWY", "P", "2"))
        self.assertNotIn(DROPDOWN, mine)

    def test_the_ids_only_preset_is_what_c0_r2_ran(self):
        """C0's R2 (a test job on HIVE): vendor Basic-Search with the same rule keys and no dropdown,
        pooled, run for real -- its fragger.params (Basic-Search's own isotope errors 0/1/2/3)."""
        import fragpipe_params as fpp
        ran = fpp.read_fragger_params(os.path.join(FIX, "c0_r2_basic.fragger.params"))
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--ids-only", *self.ARGS))
        mine = props(text)
        self.assertEqual(prov["template"], "fragpipe/Basic-Search.workflow")
        for key in ("search_enzyme_name_1", "search_enzyme_cut_1", "search_enzyme_nocut_1",
                    "search_enzyme_sense_1", "num_enzyme_termini", "digest_min_length",
                    "digest_max_length", "max_variable_mods_per_peptide", "isotope_error",
                    "allowed_missed_cleavage_1"):
            with self.subTest(key=key):
                self.assertEqual(mine["msfragger." + key], ran[key])
        self.assertNotIn(DROPDOWN, mine)

    def test_the_record_reads_the_rule_not_the_name(self):
        import cleavage_rules as cr
        import fragpipe_params as fpp
        ran = fpp.read_fragger_params(os.path.join(FIX, "c0_r1.fragger.params"))
        cl = cr.fragger_cleavage(fpp.digest(ran), "fragger.params")
        self.assertEqual((cl["name"], cl["specificity"]), ("Chymotrypsin", "full"))
        self.assertNotIn("nonspecific", cl["rule"])
        self.assertIn("cut FLMWY", cl["rule"])
        with tempfile.TemporaryDirectory() as d:
            _, _, prov = ok(self, preset(d, *self.ARGS))
        self.assertEqual(prov["cleavage"]["name"], "Chymotrypsin")
        self.assertNotIn("nonspecific", prov["cleavage"]["rule"])

    def test_a_dropdown_relabelled_rule_never_names_the_dropdown(self):
        """H5's real fragger.params: its `custom` dropdown made FragPipe call the rule argc."""
        import cleavage_rules as cr
        dg = {"specificity": "full", "termini": 2, "enzymes": [
            {"name": "argc", "cut": "FLMWY", "nocut": "P", "sense": "C", "missed": 2}]}
        cl = cr.fragger_cleavage(dg, "fragger.params")
        self.assertNotIn("argc", cl["rule"])
        self.assertEqual(cl["name"], "Chymotrypsin")
        none = cr.fragger_cleavage({"specificity": "none", "termini": 0, "enzymes": [
            {"name": "argc", "cut": "@", "nocut": "", "sense": "C", "missed": 2}]}, "x")
        self.assertNotIn("argc", none["rule"])


class TemplateChoice(unittest.TestCase):
    def test_amounts_wanted_is_lfq_mbr(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d))
        self.assertEqual(prov["template"], "fragpipe/LFQ-MBR.workflow")
        self.assertEqual(props(text)["quantitation.run-label-free-quant"], "true")
        self.assertIs(prov["quantified"], True)
        self.assertIs(prov["ids_only"], False)

    def test_timstof_tolerances_are_said_to_be_calibration_starting_values(self):
        """C2 proteomics M5: LFQ-MBR's +/-20 ppm on timsTOF; the record says the workflow's
        calibrate_mass=2 makes them starting values, never that FragPipe "tunes" them."""
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d))
        self.assertEqual(props(text)["msfragger.calibrate_mass"], "2")
        src = prov["tolerances"]["source"]
        self.assertEqual(prov["tolerances"]["mode"], "vendor")
        self.assertIn("timsTOF too", src)
        self.assertIn("calibrate_mass=2", src)
        self.assertNotIn("tunes", src)

    def test_ids_only_is_basic_search(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--ids-only"))
        self.assertEqual(prov["template"], "fragpipe/Basic-Search.workflow")
        self.assertEqual(props(text)["quantitation.run-label-free-quant"], "false")
        self.assertIs(prov["quantified"], False)
        self.assertIn("--ids-only", prov["template_choice"]["why"])

    def test_no_enzyme_is_the_nonspecific_peptidome_with_lfq_only_for_amounts(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--enzyme", "none"))
            _, ids_text, ids_prov = ok(self, preset(d, "--enzyme", "none", "--ids-only",
                                                    name="ids.workflow"))
        self.assertEqual(prov["template"], "fragpipe/Nonspecific-peptidome.workflow")
        self.assertEqual(props(text)["quantitation.run-label-free-quant"], "true")
        self.assertIs(prov["quantified"], True)
        self.assertEqual(ids_prov["template"], "fragpipe/Nonspecific-peptidome.workflow")
        self.assertEqual(props(ids_text)["quantitation.run-label-free-quant"], "false")
        self.assertIs(ids_prov["quantified"], False)
        for t, pr in ((text, prov), (ids_text, ids_prov)):
            self.assertEqual(props(t)["msfragger.num_enzyme_termini"], "0")
            self.assertEqual(props(t)[DROPDOWN], "nonspecific", "the vendor's own, consistent")
            self.assertEqual(pr["cleavage"]["specificity"], "none")

    def test_hla_on_request(self):
        with tempfile.TemporaryDirectory() as d:
            _, _, prov = ok(self, preset(d, "--enzyme", "none", "--nonspecific-workflow", "hla"))
            p, *_ = preset(d, "--nonspecific-workflow", "hla", name="x.workflow")
        self.assertEqual(prov["template"], "fragpipe/Nonspecific-HLA.workflow")
        self.assertNotEqual(p.returncode, 0, "an enzyme rule has no non-specific workflow")
        self.assertIn("--enzyme none", p.stderr)

    def test_pepsin_is_searched_non_specifically(self):
        with tempfile.TemporaryDirectory() as d:
            _, _, prov = ok(self, preset(d, "--enzyme", "pepsin"))
        self.assertEqual(prov["template"], "fragpipe/Nonspecific-peptidome.workflow")
        self.assertIn("pepsin", prov["enzyme"]["why"])

    def test_no_database_path_or_threads_even_when_given(self):
        """run_search.py sets both at search time (F3): the database from --fasta (or the .fas
        it builds), the threads from the job's CPUs."""
        with tempfile.TemporaryDirectory() as d:
            for extra in ((), ("--ids-only",), ("--enzyme", "none")):
                _, text, prov = ok(self, preset(d, "--fasta", "/x/db.fasta", "--threads", "8",
                                                *extra))
                self.assertNotIn("database.db-path", props(text))
                self.assertNotIn("workflow.threads", props(text))
                self.assertIn("run_search", prov["database"])
                self.assertIn("--threads 8", prov["threads"])
                self.assertIn(props(text)["quantitation.run-label-free-quant"], ("true", "false"),
                              "every DDA preset says whether it quantifies")
                self.assertEqual(props(text)["database.decoy-tag"], "rev_")

    def test_the_vendor_data_type_keys_are_kept(self):
        """review M6: a .d forces IM-MS in FragPipe itself; a timsTOF mzML may have none."""
        with tempfile.TemporaryDirectory() as d:
            _, text, _ = ok(self, preset(d))
        p = props(text)
        self.assertEqual((p["workflow.input.data-type.im-ms"],
                          p["workflow.input.data-type.regular-ms"]), ("false", "true"))

    def test_dda_only_flags_on_a_dia_preset_are_refused(self):
        with tempfile.TemporaryDirectory() as d:
            for args in (("--ids-only",), ("--enzyme", "trypsin"), ("--var-mods", "ox"),
                         ("--overrides", "{}")):
                with self.subTest(args=args):
                    p, out, _, _ = preset(d, *args, acquisition="DIA")
                    self.assertEqual(p.returncode, 2, p.stderr)
                    self.assertIn("for the FragPipe DDA preset", p.stderr)
                    self.assertFalse(os.path.exists(out))


class IonTrapAndUnitResolution(unittest.TestCase):
    def vendor(self, name):
        with open(os.path.join(FIX, name)) as fh:
            return props(fh.read())

    def test_ion_trap_ms2_takes_fragpipes_tmt_ms3_setting(self):
        tmt = self.vendor("TMT10-MS3.workflow")
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--ms1-res", "120000", "--ms2-analyzer", "ITMS",
                                            "--resolution-source", "detected",
                                            instrument="Orbitrap Fusion Lumos"))
        p = props(text)
        for key in ("msfragger.fragment_mass_tolerance", "msfragger.fragment_mass_units",
                    "msfragger.isotope_error", "msfragger.precursor_mass_lower",
                    "msfragger.precursor_mass_upper", "msfragger.precursor_mass_units"):
            with self.subTest(key=key):
                self.assertEqual(p[key], tmt[key])
        self.assertEqual((p["msfragger.fragment_mass_tolerance"],
                          p["msfragger.fragment_mass_units"], p["msfragger.isotope_error"]),
                         ("0.6", "0", "0/1/2/3"))
        self.assertEqual(prov["tolerances"]["mode"], "ion_trap")
        self.assertIn("TMT10-MS3", prov["tolerances"]["source"])
        self.assertIn("adapted", prov["tolerances"]["source"])

    def test_mixed_analyzers_search_every_ms2_at_the_ion_trap_window(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--ms1-res", "120000", "--ms2-analyzer", "mixed",
                                            instrument="Orbitrap Fusion Lumos"))
        p = props(text)
        self.assertEqual(p["msfragger.fragment_mass_tolerance"], "0.6")
        self.assertEqual(p["msfragger.analyzer_types"], "all", "no analyzer filter")
        self.assertIn("analyzer", prov["tolerances"]["source"])
        self.assertIn("mixed", prov["tolerances"]["source"])

    def test_unit_resolution_takes_fragpipes_stellar_setting(self):
        st = self.vendor("Stellar-DDA.workflow")
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--ids-only", instrument="Stellar"))
        p = props(text)
        self.assertEqual((p["msfragger.fragment_mass_tolerance"],
                          p["msfragger.fragment_mass_units"]), ("0.3", "0"))
        self.assertEqual((p["msfragger.precursor_mass_lower"], p["msfragger.precursor_mass_upper"],
                          p["msfragger.precursor_mass_units"]), ("-1.5", "1.5", "0"))
        import make_presets as mp
        for key, val in mp.UNIT_RESOLUTION_KEYS.items():
            with self.subTest(key=key):
                self.assertEqual(val, st[key])
                self.assertEqual(p[key], val)
        self.assertEqual(prov["tolerances"]["mode"], "unit_resolution")
        self.assertIn("Stellar-DDA", prov["tolerances"]["source"])

    def test_unit_resolution_is_never_quantified(self):
        """Ruling 105: unit-resolution MS1 cannot support label-free MS1 quantities, and
        FragPipe's own Stellar-DDA workflow runs none -- identifications only, said why."""
        with tempfile.TemporaryDirectory() as d:
            for extra in ((), ("--enzyme", "none")):
                with self.subTest(extra=extra):
                    p, out, _, _ = preset(d, *extra, instrument="Stellar")
                    self.assertEqual(p.returncode, 1, p.stderr)
                    self.assertFalse(os.path.exists(out))
                    for words in ("unit resolution", "--ids-only", "PRM", "Skyline"):
                        self.assertIn(words, p.stderr)

    def test_the_ion_trap_keys_are_the_vendors(self):
        import make_presets as mp
        tmt = self.vendor("TMT10-MS3.workflow")
        for key, val in mp.ION_TRAP_KEYS.items():
            self.assertEqual(val, tmt[key], key)

    def test_a_ppm_fragment_override_on_an_ion_trap_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p, out, _, _ = preset(d, "--ms2-analyzer", "ITMS", "--ms2-ppm", "20",
                                  instrument="Orbitrap Fusion Lumos")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("ion trap", p.stderr)


class EnzymeVocabulary(unittest.TestCase):
    def test_fragpipes_enzyme_list_is_read_right(self):
        import cleavage_rules as cr
        got = {}
        with open(os.path.join(FIX, "enzyme-defs-supported.24.0.txt")) as fh:
            for ln in fh:
                m = re.match(r"^(\S+)\s*:\s*cut\(([^)]*)\)\s*(?:nocuts\(([^)]*)\))?\s*sense\((\w)\)",
                             ln.strip())
                if m:
                    got[m.group(1)] = (m.group(2), m.group(3) or "", m.group(4))
        self.assertEqual(cr.FRAGPIPE_ENZYMES, got)

    def test_the_vocabulary_is_fetch_fastas(self):
        import cleavage_rules as cr
        import fetch_fasta as ff
        self.assertEqual(set(cr.ENZYME_RULES), set(ff.ENZYME_FAMILIES))

    def test_each_word_is_fragpipes_rule(self):
        import cleavage_rules as cr
        want = {"trypsin": ("KR", "P", "C", "Trypsin"), "lysc": ("K", "P", "C", "Lys-C"),
                "argc": ("R", "P", "C", "Arg-C"), "aspn": ("D", "", "N", "Asp-N"),
                "lysn": ("K", "", "N", "Lys-N"), "chymotrypsin": ("FLWY", "P", "C", "Chymotrypsin"),
                "gluc": ("DE", "P", "C", "V8-DE")}
        for word, (cut, nocut, sense, name) in want.items():
            with self.subTest(word=word):
                r = cr.enzyme_rule(word)
                self.assertEqual((r["cut"], r["nocut"], r["sense"], r["name"], r["specificity"]),
                                 (cut, nocut, sense, name, "full"))

    def test_lists_slash_p_none_and_refusals(self):
        import cleavage_rules as cr
        r = cr.enzyme_rule("trypsin,lysc")
        self.assertEqual((r["cut"], r["nocut"], r["name"]), ("KR", "P", "Trypsin"))
        self.assertEqual(r["families"], ["lysc", "trypsin"])
        self.assertEqual((cr.enzyme_rule("Trypsin/P")["nocut"],
                          cr.enzyme_rule("Trypsin/P")["name"]), ("", "Trypsin/P"))
        self.assertEqual(cr.enzyme_rule("none")["specificity"], "none")
        self.assertEqual(cr.enzyme_rule("trypsin,pepsin")["specificity"], "none")
        for bad in ("trypsn", "trypsin,aspn", "none,trypsin", ""):
            with self.subTest(bad=bad):
                with self.assertRaises(cr.EnzymeError) as e:
                    cr.enzyme_rule(bad)
                if bad == "trypsn":
                    self.assertIn("chymotrypsin", str(e.exception), "names the vocabulary")

    def test_a_rule_that_is_a_list_entry_is_written_under_its_name(self):
        """X8 F6: chymotrypsin is a row of FragPipe's list: name and dropdown say so (FragPipe
        reloads the same cut/no-cut/sense from it), not `nonspecific` beside the cut rule."""
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--enzyme", "chymotrypsin"))
        p = props(text)
        self.assertEqual(p[DROPDOWN], "chymotrypsin")
        self.assertEqual((p["msfragger.search_enzyme_name_1"], p["msfragger.search_enzyme_cut_1"],
                          p["msfragger.search_enzyme_nocut_1"],
                          p["msfragger.search_enzyme_sense_1"], p["msfragger.num_enzyme_termini"]),
                         ("chymotrypsin", "FLWY", "P", "C", "2"))
        self.assertEqual(prov["cleavage"]["name"], "Chymotrypsin")
        self.assertIn("--enzyme chymotrypsin", prov["keys_changed"]["msfragger.search_enzyme_cut_1"])

    def test_no_enzyme_is_the_shared_dda_digest(self):
        """Ruling 57.1: "The FragPipe digest is the skill's trypsin (no cut before P), one digest
        shared with the Sage fallback" -- not the vendor's stricttrypsin (Trypsin/P), 7-50."""
        import estimate_params as ep
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d))
            _, ttext, _ = ok(self, preset(d, "--enzyme", "trypsin", name="t.workflow"))
        p = props(text)
        self.assertEqual(p[DROPDOWN], "trypsin")                  # FragPipe's `trypsin` (X8 F6)
        self.assertEqual((p["msfragger.search_enzyme_name_1"], p["msfragger.search_enzyme_cut_1"],
                          p["msfragger.search_enzyme_nocut_1"], p["msfragger.search_enzyme_sense_1"],
                          p["msfragger.num_enzyme_termini"]), ("trypsin", "KR", "P", "C", "2"))
        dg = ep.DDA_DIGEST
        self.assertEqual((p["msfragger.allowed_missed_cleavage_1"], p["msfragger.digest_min_length"],
                          p["msfragger.digest_max_length"]),
                         (str(dg["missed_cleavages"]), str(dg["min_len"]), str(dg["max_len"])))
        self.assertEqual(prov["cleavage"]["name"], "Trypsin")
        self.assertIn("DEFAULT", prov["cleavage"]["source"])
        self.assertIn("57.1", prov["cleavage"]["source"])
        self.assertEqual(text, ttext, "--enzyme trypsin is the default digest")
        self.assertIs(ep.SAGE_ENZYME, ep.DDA_DIGEST)

    def sage(self, d, *extra, name="sage.json"):
        out = os.path.join(d, name)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "sage", "--acquisition", "DDA", "--instrument", TIMS,
                            "--out", out, *extra], capture_output=True, text=True, timeout=60)
        return p, out

    def test_sage_takes_the_same_words(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self.sage(d, "--enzyme", "chymotrypsin")
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                enz = json.load(fh)["database"]["enzyme"]
            self.assertEqual((enz["cleave_at"], enz["restrict"]), ("FLWY", "P"))
            with open(out + ".rationale.json") as fh:
                rat = json.load(fh)
            self.assertIn("--enzyme chymotrypsin", rat["rationale"]["enzyme"]["source"])
            p, out = self.sage(d, "--enzyme", "none", name="n.json")
            with open(out) as fh:
                self.assertEqual(json.load(fh)["database"]["enzyme"]["cleave_at"], "")
            p, out = self.sage(d, "--enzyme", "aspn", name="a.json")
            with open(out) as fh:
                enz = json.load(fh)["database"]["enzyme"]
            self.assertEqual((enz["cleave_at"], enz.get("restrict"), enz["c_terminal"]),
                             ("D", None, False))

    def test_sage_trypsin_is_todays_cfg(self):
        with tempfile.TemporaryDirectory() as d:
            p1, a = self.sage(d)
            p2, b = self.sage(d, "--enzyme", "trypsin", name="t.json")
            with open(a, "rb") as fa, open(b, "rb") as fb:
                self.assertEqual(fa.read(), fb.read())

    def test_diann_refuses_enzyme(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "params.cfg")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                                "--engine", "diann", "--acquisition", "DIA", "--instrument", TIMS,
                                "--enzyme", "trypsin", "--out", out],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 2, p.stderr)
            self.assertIn("DIA-NN's digest is --cut K*,R*", p.stderr)
            self.assertFalse(os.path.exists(out))


class DirectLfqRequest(unittest.TestCase):
    """--directlfq: IonQuant's MaxLFQ off (ionquant.maxlfq=0), recorded as a request."""

    def test_the_flag_writes_maxlfq_off_and_records_the_request(self):
        with tempfile.TemporaryDirectory() as d:
            out, text, prov = ok(self, preset(d, "--directlfq"))
            self.assertEqual(props(text)["ionquant.maxlfq"], "0")
            self.assertEqual(prov["quantity_requested"], "directlfq")
            self.assertIn("--directlfq", prov["keys_changed"]["ionquant.maxlfq"])
            out, text, prov = ok(self, preset(d, name="plain.workflow"))
            self.assertEqual(props(text)["ionquant.maxlfq"], "1")
            self.assertIsNone(prov["quantity_requested"])

    def test_it_is_refused_with_ids_only_with_a_conflicting_override_and_for_sage(self):
        with tempfile.TemporaryDirectory() as d:
            p, out, _, _ = preset(d, "--directlfq", "--ids-only")
            self.assertNotEqual(p.returncode, 0)
            self.assertFalse(os.path.exists(out))
            p, out, _, _ = preset(d, "--directlfq", "--overrides", '{"ionquant.maxlfq": 1}',
                                  name="c.workflow")
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("give one or the other", p.stdout + p.stderr)
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                                "--engine", "sage", "--acquisition", "DDA", "--directlfq",
                                "--instrument", TIMS, "--out", os.path.join(d, "sage.json")],
                               capture_output=True, text=True, timeout=60)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("FragPipe DDA", p.stdout + p.stderr)

    def test_the_overrides_route_still_gives_the_same_request(self):
        with tempfile.TemporaryDirectory() as d:
            out, text, prov = ok(self, preset(d, "--overrides", '{"ionquant.maxlfq": 0}'))
            self.assertEqual(props(text)["ionquant.maxlfq"], "0")


class OneFlagSet(unittest.TestCase):
    FLAGS = ("--var-mods", "--ids-only", "--cys-alkylation", "--overrides", "--override-by",
             "--override-reason", "--enzyme", "--directlfq")

    def test_one_definition(self):
        import estimate_params as ep
        import make_presets as mp
        self.assertIs(mp.add_dda_args, ep.add_dda_args)
        for script in ("estimate_params.py", "make_presets.py"):
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, script), "--help"],
                               capture_output=True, text=True, timeout=60)
            for flag in self.FLAGS:
                self.assertIn(flag, p.stdout, f"{script} {flag}")

    def test_non_specific_workflows_keep_fragpipes_variable_mod_limit(self):
        """Ruling 105: the skill's limit of 2 is for enzymatic presets; Nonspecific-peptidome and
        -HLA keep FragPipe's own 3 (what H7 validated), recorded as FragPipe's default."""
        with tempfile.TemporaryDirectory() as d:
            for wf in ("peptidome", "hla"):
                with self.subTest(workflow=wf):
                    _, text, prov = ok(self, preset(d, "--enzyme", "none", "--nonspecific-workflow",
                                                    wf, name=f"{wf}.workflow"))
                    self.assertEqual(props(text)["msfragger.max_variable_mods_per_peptide"], "3")
                    self.assertEqual(prov["max_var_mods"]["value"], 3)
                    self.assertIn("FragPipe's own", prov["max_var_mods"]["source"])
                    self.assertIn("DEFAULT", prov["max_var_mods"]["source"])
                    self.assertNotIn("msfragger.max_variable_mods_per_peptide",
                                     prov["keys_changed"])
            _, text, prov = ok(self, preset(d, "--ids-only", name="ids.workflow"))
            self.assertEqual(props(text)["msfragger.max_variable_mods_per_peptide"], "2")
            self.assertEqual(prov["max_var_mods"]["value"], 2)

    def test_shared_constants(self):
        import estimate_params as ep
        self.assertEqual(ep.DDA_MAX_VAR_MODS, 2)
        self.assertEqual(ep.PROTEIN_NTERM_ACETYL, 42.0106)
        self.assertEqual(ep.VAR_MODS_DEFAULT[("fragpipe", "DDA")], ["ox"])
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d))
            out = os.path.join(d, "s.json")
            subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "sage", "--acquisition", "DDA", "--out", out],
                           capture_output=True, text=True, timeout=60, check=True)
            with open(out) as fh:
                sage = json.load(fh)["database"]
        self.assertEqual(int(props(text)["msfragger.max_variable_mods_per_peptide"]),
                         ep.DDA_MAX_VAR_MODS)
        self.assertEqual(sage["max_variable_mods"], ep.DDA_MAX_VAR_MODS)
        self.assertEqual(sage["variable_mods"]["["], [ep.PROTEIN_NTERM_ACETYL])
        self.assertIn(f"{ep.PROTEIN_NTERM_ACETYL},[^,true,1", var_rows(text))

    def test_met_ox_by_default_and_on_request(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d))
            self.assertIn("15.9949,M,true,3", var_rows(text))
            self.assertIn("DEFAULT", prov["var_mods"]["source"])
            _, text, prov = ok(self, preset(d, "--var-mods", "none", name="n.workflow"))
            self.assertIn("15.9949,M,false,3", var_rows(text))
            _, text, prov = ok(self, preset(d, "--var-mods", "ox,deam", name="od.workflow"))
            self.assertIn("15.9949,M,true,3", var_rows(text))
            self.assertIn("0.984016,NQ,true,2", var_rows(text))
            self.assertEqual(prov["var_mods"]["words"], ["ox", "deam"])
            p, out, _, _ = preset(d, "--var-mods", "phospho", name="bad.workflow")
            self.assertNotEqual(p.returncode, 0)
            self.assertFalse(os.path.exists(out))
            self.assertIn("ox, deam, none", p.stderr)

    def test_overrides_are_checked_and_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(
                d, "--overrides", json.dumps({"msfragger.allowed_missed_cleavage_1": 1}),
                "--override-reason", "the lab's SOP"))
            self.assertEqual(props(text)["msfragger.allowed_missed_cleavage_1"], "1")
            rec = prov["overrides"]["msfragger.allowed_missed_cleavage_1"]
            self.assertEqual(rec["reason"], "the lab's SOP")
            self.assertIn("user override", prov["keys_changed"]["msfragger.allowed_missed_cleavage_1"])
            for bad in ({"msfragger.no_such_key": 1}, {"database.db-path": "/x.fas"},
                        {DROPDOWN: "trypsin"}):
                with self.subTest(bad=bad):
                    p, out, _, _ = preset(d, "--overrides", json.dumps(bad), name="b.workflow")
                    self.assertNotEqual(p.returncode, 0)
                    self.assertFalse(os.path.exists(out))
            p, *_ = preset(d, "--ids-only", "--overrides",
                           json.dumps({"quantitation.run-label-free-quant": "true"}),
                           name="q.workflow")
            self.assertNotEqual(p.returncode, 0)


class Cysteine(unittest.TestCase):
    """--cys-alkylation on the DDA preset, and what the written workflow fixes on C."""

    def test_hla_without_an_answer_keeps_its_own_unmodified_cysteine(self):
        """Nonspecific-HLA ships its fixed C row OFF (HLA peptides are usually not alkylated):
        not given, that is the workflow's default -- never recorded as carbamidomethyl."""
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--enzyme", "none", "--nonspecific-workflow", "hla"))
        self.assertIn("57.02146,C (cysteine),false,-1", props(text)["msfragger.table.fix-mods"])
        self.assertIn("119.0041,C,true,1", var_rows(text))
        cys = prov["cys_alkylation"]
        self.assertEqual((cys["value"], cys["given"], cys["mod"]), ("none", False, None))
        self.assertIn("DEFAULT of the vendor workflow", cys["source"])

    def test_alkylated_hla_fixes_c_and_drops_cysteinylation(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--enzyme", "none", "--nonspecific-workflow", "hla",
                                            "--cys-alkylation", "cam"))
        self.assertIn("57.021464,C (cysteine),true,-1", props(text)["msfragger.table.fix-mods"])
        self.assertIn("119.0041,C,false,1", var_rows(text))
        self.assertEqual(prov["cys_alkylation"]["value"], "cam")

    def test_an_agent_replaces_carbamidomethyl(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = ok(self, preset(d, "--cys-alkylation", "nem"))
        self.assertIn("; 125.047679,C (cysteine),true,-1;", props(text)["msfragger.table.fix-mods"])
        self.assertEqual(prov["cys_alkylation"]["value"], "nem")

    def test_cys_searched_reads_fragpipe_files(self):
        """estimate_params.cys_searched, the one reader, read a .workflow as a DIA-NN cfg."""
        import estimate_params as ep
        self.assertEqual(ep.cys_searched(os.path.join(PRESETS, "LFQ-MBR.workflow"))["choice"], "cam")
        self.assertEqual(ep.cys_searched(os.path.join(PRESETS, "Nonspecific-HLA.workflow")),
                         {"fixed": False})
        self.assertEqual(ep.cys_searched(os.path.join(FIX, "c0_r1.fragger.params"))["choice"],
                         "cam")

    def test_a_fix_mods_override_is_recorded_as_written(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(PRESETS, "LFQ-MBR.workflow")) as fh:
                row = props(fh.read())["msfragger.table.fix-mods"]
            _, _, prov = ok(self, preset(d, "--overrides", json.dumps(
                {"msfragger.table.fix-mods": row.replace("57.02146,C (cysteine)",
                                                         "0.0,C (cysteine)")})))
        cys = prov["cys_alkylation"]
        self.assertEqual(cys["overridden"]["was"], "cam")
        self.assertIn("msfragger.table.fix-mods", cys["source"])
        self.assertNotIn("--fixed-mod", cys["source"])


class SharedDigest(unittest.TestCase):
    """review L7: run_search.warn_contaminant_digest skips FragPipe; the digest it needs."""

    def test_the_digest_of_each_preset(self):
        import cleavage_rules as cr
        with tempfile.TemporaryDirectory() as d:
            out, _, _ = ok(self, preset(d))
            self.assertEqual(cr.fragpipe_digest(out), {
                "cut": "K*,R*,!*P", "missed_cleavages": 2, "min_pep_len": 7, "max_pep_len": 30,
                "met_excision": True})
            out, _, _ = ok(self, preset(d, "--enzyme", "chymotrypsin", name="c.workflow"))
            self.assertEqual(cr.fragpipe_digest(out)["cut"], "F*,L*,W*,Y*,!*P")
            out, _, _ = ok(self, preset(d, "--enzyme", "none", name="n.workflow"))
            self.assertEqual(cr.fragpipe_digest(out)["cut"], "(none: a non-specific digest)")
        ran = cr.fragpipe_digest(os.path.join(FIX, "c0_r1.fragger.params"))
        self.assertEqual((ran["cut"], ran["min_pep_len"], ran["max_pep_len"]),
                         ("F*,L*,M*,W*,Y*,!*P", 6, 30))
        self.assertIsNone(cr.fragpipe_digest(os.path.join(FIX, "missing.workflow")))

    def test_one_cut_syntax_with_sage(self):
        import cleavage_rules as cr
        import run_search
        self.assertEqual(cr.diann_cut("KR", "P", True),
                         run_search.sage_digest({"cleave_at": "KR", "restrict": "P"})["cut"])

    def test_one_digest_with_sage(self):
        """Ruling 57.1: the FragPipe preset's default digest IS Sage's (estimate_params.DDA_DIGEST),
        so run_search knows one default difference for both."""
        import cleavage_rules as cr
        import estimate_params as ep
        import run_search
        with tempfile.TemporaryDirectory() as d:
            for extra in ((), ("--ids-only",)):
                out, _, _ = ok(self, preset(d, *extra, name=f"x{len(extra)}.workflow"))
                fp = cr.fragpipe_digest(out)
                fp.pop("met_excision")
                self.assertEqual(fp, run_search.sage_digest(ep.DDA_DIGEST))

    def test_a_dropdown_that_would_overwrite_the_rule_is_named(self):
        """C0 fact 1, and FragPipe's own saved workflow: it writes the dropdown back (here
        `trypsin`) beside the rule it ran, so a replay of THAT file would relabel the rule."""
        import cleavage_rules as cr
        import fragpipe_params as fpp
        saved, _ = fpp.read_workflow(os.path.join(FIX, "c0_v1a_dryrun.fragpipe.workflow"))
        self.assertIn("trypsin", cr.fragpipe_dropdown_conflict(saved))
        no_p = dict(saved, **{"msfragger.search_enzyme_nocut_1": ""})
        self.assertIn("no cut - -> P", cr.fragpipe_dropdown_conflict(no_p))
        for name in ("LFQ-MBR", "Basic-Search", "Nonspecific-peptidome", "Nonspecific-HLA"):
            vendor, _ = fpp.read_workflow(os.path.join(PRESETS, name + ".workflow"))
            self.assertIsNone(cr.fragpipe_dropdown_conflict(vendor), name)
        with tempfile.TemporaryDirectory() as d:
            out, _, _ = ok(self, preset(d, "--enzyme", "chymotrypsin"))
            mine, _ = fpp.read_workflow(out)
        self.assertIsNone(cr.fragpipe_dropdown_conflict(mine))

    def test_h5s_own_workflow_is_described_by_what_ran(self):
        """Ruling 104: a user's or vendor's workflow is never refused for its enzyme dropdown.
        H5's hand-made workflow (dropdown `custom` beside FLMWY / no-cut P) ran as `argc` (H5's
        real fragger.params, its database path replaced); the record reads the rule from that
        fragger.params, by its cuts, and nothing here refuses the workflow."""
        import cleavage_rules as cr
        import fragpipe_params as fpp
        with open(os.path.join(FIX, "h5_golden_lfqmbr.workflow")) as fh:
            h5 = props(fh.read())
        self.assertIn("argc", cr.fragpipe_dropdown_conflict(h5), "said, as information")
        self.assertFalse(hasattr(cr, "fragpipe_dropdown_changes_search"), "no refusal helper")
        ran = fpp.read_fragger_params(os.path.join(FIX, "h5.fragger.params"))
        self.assertEqual(ran["search_enzyme_name_1"], "argc")
        cl = cr.fragger_cleavage(fpp.digest(ran), "fragger.params")
        self.assertEqual((cl["name"], cl["ac"], cl["specificity"]),
                         ("Chymotrypsin", "MS:1001306", "full"))
        self.assertEqual(cl["rule"], "cut FLMWY, no cut P, sense C; num_enzyme_termini 2")
        dg = cr.fragpipe_digest(os.path.join(FIX, "h5.fragger.params"))
        self.assertEqual((dg["cut"], dg["min_pep_len"], dg["max_pep_len"]),
                         ("F*,L*,M*,W*,Y*,!*P", 6, 30))

    def test_the_skills_own_preset_is_refused_with_a_rule_under_a_dropdown(self):
        """Ruling 104's other half: the skill never writes a custom rule under a dropdown name."""
        import make_presets as mp
        with tempfile.TemporaryDirectory() as d:
            out, text, prov = ok(self, preset(d, "--enzyme", "chymotrypsin"))
            with open(out, "a") as fh:
                fh.write(f"{DROPDOWN}=trypsin\n")
            with self.assertRaises(ValueError) as e:
                mp.dda_record_as_written(dict(prov, warnings=[]), out,
                                         {"enzyme": None, "overrides": {}, "override_meta": {}})
        self.assertIn("would not search the rule it states", str(e.exception))


# sha256 of the outputs at 2f9b229c (computed there, 2026-11-04): cfg bytes, and the rationale /
# provenance JSON without its path (None: it holds an override's role, which depends on the host)
BASE_PINS = {
    # Sage re-pinned 2.15.1 (follow-up): the dda lane's 0e99e604 changed the cfg on purpose --
    # chimera off (it needs report_psms > 1), max_fragment_charge 1 for high-res MS2 (ion trap
    # stays 2), and the enzyme/charge/chimera rationale text. Only those fields differ from 2f9b229c.
    "sage_timstof_dda": ("2c1277a529797802a175a6b17e44488d8d3fe909c0c4a1551341ddef8c3ca7dd",
        "39fda83f8259a509b7e2de54f21b62cf0ce80d80cbfa083236227ef8b6a4ed24"),
    "sage_timstof_ids": ("d0eb6fc08fa1636523d3192822b87d24be415fe3883872ec09e7657f9d70c9bd",
        "bfb030469944e42bc988c37f6c8c5b43937b52a4d0da841b84ca5760125d6b4c"),
    "sage_lumos_itms": ("6402f5805d227fab3c29f72cedb0ad273407ccd96baafc84bd1162ef6ef31afb",
        "aa5f54ef731631bf9079a521a56fc44e0e6907fa5d085c187e8f1aff0ad6779e"),
    "sage_varmods_nocys": ("28ff76cc7755575017facc1217c5e8c598a69b9e36f3e07696016627465deb45",
        "4ab5b7f17b74a69695ff16e6e0cf09c5f01f975cc83f1f14ee120bdf37a6376c"),
    "sage_nonspecific_override": ("08bc41a84ae52fb6a31091740eecd948aa4bb6792776ba0524a787c6a61e63ce",
        None),
    "sage_astral_dia": ("a670507cda2f954a486a5172ccc54be016c6b6fded5c8c1debf661f867d9e93c",
        "80c531a287fe878d22dd3f05f75b702ebdca8f7f2b97e57fa8b281dbc131fb85"),
    "fp_timstof_dia": ("cf67a668b01dce364f677b8b936a0575f4eec2036df17fa4b0fced96c49b4c1b",
                       "41ba51aa37beb8253484cd553b78fac3e0c121bee0c8efb63d741a52fb519f8e"),
    "fp_orbitrap_dia": ("09e11fd2192ccf0b98ac514af096863dfcdc41df5d420726c320967233470793",
                        "f041f15579521e801feb7a5e7059ed00c4a3ae54c67707c1c4db4b49c3a26bee"),
    "fp_timstof_nem_db": ("3ac6f7e74c230183db3d2c051c0ab810a85df064061cfd686e2d828789cff4bb",
                          "710922462203510e0f40737fd17c497670e4a99f81546680596c4bb480f2be06"),
    "fp_orbitrap_ppm": ("5dd05f014375a9630d3df4500210ffa1ab661495631be9314c6266acdf4db660", None),
    "radiant_astral": ("3e0f0ff2b45b468287e50f7961d7bfeda4a6dd7d717e6db0582a2d58925ca343",
                       # sidecar re-pinned 2.15.1: + route_status, widths tagged UNVALIDATED
                       # (audit 2026-11-17 L1-3); the cfg bytes are unchanged
                       "28717ea5e98688f2db94292a16b190310e68ac9829552b2f4a2c01b925de463e"),
    "radiant_exploris_res": ("72a8bb62b07145f8fe512cb5a3daf7059c09c9999e69609aafd3d77695b78df2",
                             "28030a90234b4a34c86a4d633a8ba32989a1da57a6d95fb8807651af3d8b4f2f"),
}
SAGE = {
    "sage_timstof_dda": ["--engine", "sage", "--acquisition", "DDA", "--instrument", TIMS],
    "sage_timstof_ids": ["--engine", "sage", "--acquisition", "DDA", "--instrument", TIMS,
                         "--ids-only"],
    "sage_lumos_itms": ["--engine", "sage", "--acquisition", "DDA", "--instrument",
                        "Orbitrap Fusion Lumos", "--ms1-resolution", "120000", "--ms2-analyzer",
                        "ITMS", "--resolution-source", "detected"],
    "sage_varmods_nocys": ["--engine", "sage", "--acquisition", "DDA", "--instrument",
                           "Orbitrap Exploris 480", "--var-mods", "ox,deam",
                           "--cys-alkylation", "none"],
    "sage_nonspecific_override": ["--engine", "sage", "--acquisition", "DDA", "--instrument",
                                  "Orbitrap Exploris 480", "--overrides",
                                  '{"database": {"enzyme": {"cleave_at": ""}}}'],
    "sage_astral_dia": ["--engine", "sage", "--acquisition", "DIA", "--instrument",
                        "Orbitrap Astral"],
}
DIA = {
    "fp_timstof_dia": ["--engine", "fragpipe", "--instrument", TIMS],
    "fp_orbitrap_dia": ["--engine", "fragpipe", "--instrument", "Orbitrap Exploris 480"],
    "fp_timstof_nem_db": ["--engine", "fragpipe", "--instrument", TIMS, "--cys-alkylation",
                          "nem", "--fasta", "/x/db.fasta", "--threads", "8"],
    "fp_orbitrap_ppm": ["--engine", "fragpipe", "--instrument", "Orbitrap Exploris 480",
                        "--ms1-ppm", "7", "--ms2-ppm", "12"],
    "radiant_astral": ["--engine", "radiant", "--instrument", "Orbitrap Astral"],
    "radiant_exploris_res": ["--engine", "radiant", "--instrument", "Orbitrap Exploris 480",
                             "--ms1-res", "120000", "--ms2-res", "30000", "--resolution-source",
                             "detected", "--cys-alkylation", "none"],
}


class Unchanged(unittest.TestCase):
    """(pin) The Sage configs and the DIA presets are byte-identical to 2f9b229c."""

    def check(self, script, cases, sidecar, path_key):
        with tempfile.TemporaryDirectory() as d:
            for name, args in cases.items():
                with self.subTest(case=name):
                    ext = (".json" if script == "estimate_params.py" else
                           ".workflow" if "fragpipe" in args else ".radiantConfig")
                    out = os.path.join(d, name + ext)
                    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args,
                                        "--out", out], capture_output=True, text=True, timeout=60)
                    self.assertEqual(p.returncode, 0, p.stderr)
                    with open(out, "rb") as fh:
                        self.assertEqual(hashlib.sha256(fh.read()).hexdigest(),
                                         BASE_PINS[name][0])
                    if BASE_PINS[name][1]:
                        with open(out + sidecar) as fh:
                            rec = json.load(fh)
                        rec.pop(path_key)
                        self.assertEqual(hashlib.sha256(json.dumps(
                            rec, sort_keys=True).encode()).hexdigest(), BASE_PINS[name][1])

    def test_sage_cfgs(self):
        self.check("estimate_params.py", SAGE, ".rationale.json", "params_file")

    def test_dia_presets(self):
        self.check("make_presets.py", DIA, ".provenance.json", "output")


if __name__ == "__main__":
    unittest.main()
