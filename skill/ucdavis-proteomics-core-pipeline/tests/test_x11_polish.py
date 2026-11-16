#!/usr/bin/env python3
"""
2.15 acceptance polish (lane X11): the runs a Methods counts, the stray timsTOF column line, the
FragPipe DDA report description, --semi-specific, an mzML's instrument and MS2 analyser, the
contaminant-digest warning for a supplied FASTA, "identified" for an ids-only search, and the
charge-range default tag. Hermetic: temp dirs only.
"""
import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import make_methods as mm  # noqa: E402
import analysis_prompt as ap  # noqa: E402
import fragpipe_quant as fq  # noqa: E402
import detect_acquisition as da  # noqa: E402
import run_search as rs  # noqa: E402
import cleavage_rules as cr  # noqa: E402

TAG = "[facility default — confirm]"


def run(script, *args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script), *args],
                          capture_output=True, text=True, timeout=120)


class RunsCounted(unittest.TestCase):
    def test_search_record_carries_the_runs_the_search_used(self):
        with tempfile.TemporaryDirectory() as d:
            prov = os.path.join(d, "search_provenance.json")
            with open(prov, "w") as fh:
                json.dump({"engine": "sage", "version": "x", "n_files": 183}, fh)
            self.assertEqual(mm.search_record(None, prov)["n_runs"], 183)

    def test_methods_count_the_searched_runs_not_the_staged_files(self):
        with tempfile.TemporaryDirectory() as d:
            prov = os.path.join(d, "search_provenance.json")
            with open(prov, "w") as fh:
                json.dump({"engine": "sage", "version": "x", "n_files": 183}, fh)
            raws = [os.path.join(d, f"gone_{i}.raw") for i in range(300)]
            out = os.path.join(d, "methods.md")
            p = run("make_methods.py", "--raw", *raws, "--instrument", "Orbitrap Fusion",
                    "--acquisition", "DDA", "--search-prov", prov, "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            text = open(out).read()
            self.assertIn("| Files in series | 183 |", text)
            self.assertNotIn("| Files in series | 300 |", text)

    def test_no_search_record_counts_the_files(self):
        with tempfile.TemporaryDirectory() as d:
            raws = [os.path.join(d, f"gone_{i}.raw") for i in range(5)]
            out = os.path.join(d, "methods.md")
            p = run("make_methods.py", "--raw", *raws, "--instrument", "Orbitrap Fusion",
                    "--acquisition", "DDA", "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("| Files in series | 5 |", open(out).read())


class ColumnLine(unittest.TestCase):
    def test_a_non_timstof_run_is_not_told_the_standard_column_is_a_timstof_one(self):
        rec = mm.column_record(core=False)
        self.assertNotIn("timsTOF one", rec["source"])
        self.assertIn("Core timsTOF .d", rec["source"])          # why the default does not apply


class FragPipeReportDescription(unittest.TestCase):
    def test_the_column_list_is_the_adapters(self):
        self.assertEqual(fq.dda_report_columns("PG.IonQuant.Intensity"),
                         ["Run", "Protein.Group", "PG.MaxLFQ", "Genes", "Protein.Names",
                          "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "PG.IonQuant.Intensity"])
        self.assertEqual(fq.dda_report_columns()[-1], "Lib.PG.Q.Value")

    def test_fragpipe_dda_report_is_described_per_protein_and_run(self):
        sp = {"engine": "fragpipe", "fragpipe_quantity": {
            "status": "ok", "quantity": "directLFQ", "companion_column": "PG.IonQuant.Intensity"}}
        what = ap.report_parquet_what(sp)
        self.assertIn("one row per protein group and run", what)
        self.assertIn("NOT per precursor", what)
        self.assertIn("PG.IonQuant.Intensity", what)
        self.assertIn("placeholders", what)

    def test_other_engines_keep_the_precursor_description(self):
        self.assertIn("one row per precursor and run", ap.report_parquet_what({"engine": "diann"}))
        self.assertIn("one row per precursor and run", ap.report_parquet_what({}))
        self.assertIn("one row per precursor and run", ap.report_parquet_what(
            {"engine": "fragpipe", "fragpipe_quantity": {"status": "refused"}}))


class SemiSpecific(unittest.TestCase):
    def sage(self, d, *more):
        out = os.path.join(d, "sage.json")
        p = run("estimate_params.py", "--engine", "sage", "--acquisition", "DDA", "--instrument",
                "Orbitrap Fusion", "--precursor-mz-range", "375", "1500", "--out", out, *more)
        return p, out

    def test_sage_semi_specific_default_trypsin(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self.sage(d, "--semi-specific")
            self.assertEqual(p.returncode, 0, p.stderr)
            enz = json.load(open(out))["database"]["enzyme"]
            self.assertIs(enz["semi_enzymatic"], True)
            self.assertEqual((enz["cleave_at"], enz["restrict"]), ("KR", "P"))
            prov = os.path.join(d, "search_provenance.json")
            json.dump({"engine": "sage", "version": "x", "params_file": out, "result": {}},
                      open(prov, "w"))
            text = mm.search_paragraph(mm.search_record(out, prov))
            self.assertIn("semi-specific", text)

    def test_sage_without_the_flag_is_unchanged(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self.sage(d)
            self.assertNotIn("semi_enzymatic", json.load(open(out))["database"]["enzyme"])

    def test_semi_specific_with_enzyme(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = self.sage(d, "--enzyme", "lysc", "--semi-specific")
            self.assertEqual(p.returncode, 0, p.stderr)
            enz = json.load(open(out))["database"]["enzyme"]
            self.assertEqual((enz["cleave_at"], enz["semi_enzymatic"]), ("K", True))

    def test_refusals(self):
        with tempfile.TemporaryDirectory() as d:
            p, _ = self.sage(d, "--enzyme", "none", "--semi-specific")
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("non-specific", p.stderr + p.stdout)
            p, _ = self.sage(d, "--semi-specific", "--overrides",
                             '{"database": {"enzyme": {"semi_enzymatic": false}}}')
            self.assertNotEqual(p.returncode, 0)
            out = os.path.join(d, "p.cfg")
            p = run("estimate_params.py", "--engine", "diann", "--acquisition", "DIA",
                    "--instrument", "timsTOF HT", "--precursor-mz-range", "300", "1200",
                    "--semi-specific", "--out", out)
            self.assertEqual(p.returncode, 2)

    def test_fragpipe_semi_specific_sets_one_enzymatic_terminus(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "fp.workflow")
            p = run("make_presets.py", "--engine", "fragpipe", "--acquisition", "DDA",
                    "--instrument", "Orbitrap Fusion", "--ms2-analyzer", "FTMS",
                    "--semi-specific", "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            kv = dict(ln.split("=", 1) for ln in open(out).read().splitlines()
                      if ln and not ln.startswith("#") and "=" in ln)
            self.assertEqual(kv["msfragger.num_enzyme_termini"], "1")
            self.assertEqual(kv["msfragger.search_enzyme_cut_1"], "KR")
            prov = json.load(open(out + ".provenance.json"))
            self.assertEqual(prov["cleavage"]["specificity"], "semi")
            self.assertIn("semi-specific", mm.cleavage_phrase(prov["cleavage"]))

    def test_fragpipe_default_stays_two_termini(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "fp.workflow")
            p = run("make_presets.py", "--engine", "fragpipe", "--acquisition", "DDA",
                    "--instrument", "Orbitrap Fusion", "--ms2-analyzer", "FTMS", "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("msfragger.num_enzyme_termini=2", open(out).read())

    def test_resolve_defaults_passes_it_on(self):
        src = open(os.path.join(SCRIPTS, "resolve_defaults.py")).read()
        self.assertIn('cmd.append("--semi-specific")', src)

    def test_the_rule_is_semi_only_when_asked(self):
        r = cr.enzyme_rule("trypsin")
        self.assertEqual(r["specificity"], "full")
        self.assertEqual(cr.semi_specific(r)["specificity"], "semi")
        with self.assertRaises(cr.EnzymeError):
            cr.semi_specific(cr.enzyme_rule("none"))


MZML = """<?xml version="1.0" encoding="utf-8"?>
<mzML xmlns="http://psi.hupo.org/ms/mzml"><referenceableParamGroupList count="1">
<referenceableParamGroup id="CommonInstrumentParams">
<cvParam cvRef="MS" accession="MS:1002634" name="Q Exactive Plus" value=""/>
<cvParam cvRef="MS" accession="MS:1000529" name="instrument serial number" value="S1"/>
</referenceableParamGroup></referenceableParamGroupList>
<instrumentConfigurationList count="1"><instrumentConfiguration id="IC1">
<referenceableParamGroupRef ref="CommonInstrumentParams"/></instrumentConfiguration>
</instrumentConfigurationList><run><spectrumList count="3">
{spectra}</spectrumList></run></mzML>
"""


def spectrum(level, filt, target=None):
    iso = ""
    if level == 2:
        iso = (f'<precursorList count="1"><precursor><isolationWindow>'
               f'<cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" '
               f'value="{target}"/><cvParam cvRef="MS" accession="MS:1000828" name="lower offset" '
               f'value="0.8"/><cvParam cvRef="MS" accession="MS:1000829" name="upper offset" '
               f'value="0.8"/></isolationWindow></precursor></precursorList>')
    return (f'<spectrum index="0" id="x" defaultArrayLength="0"><cvParam cvRef="MS" '
            f'accession="MS:1000511" name="ms level" value="{level}"/><scanList count="1"><scan>'
            f'<cvParam cvRef="MS" accession="MS:1000512" name="filter string" value="{filt}"/>'
            f'</scan></scanList>{iso}</spectrum>')


class MzmlInstrument(unittest.TestCase):
    def write(self, d, filts, gz=False):
        sp = "".join(spectrum(1, "FTMS + p NSI Full ms [375.0000-1500.0000]") if f is None
                     else spectrum(2, f, 400 + i) for i, f in enumerate(filts))
        p = os.path.join(d, "a.mzML.gz" if gz else "a.mzML")
        data = MZML.format(spectra=sp).encode()
        with (gzip.open(p, "wb") if gz else open(p, "wb")) as fh:
            fh.write(data)
        return p

    def test_model_is_read_not_the_serial_number(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(da.mzml_instrument(self.write(d, [None])), "Q Exactive Plus")
            self.assertEqual(da.mzml_instrument(self.write(d, [None], gz=True)), "Q Exactive Plus")

    def test_classify_reports_instrument_and_ms2_analyzer(self):
        with tempfile.TemporaryDirectory() as d:
            r = da.classify(self.write(d, [None, "FTMS + c NSI d Full ms2 400.0@hcd30.00 [120-900]"]))
            self.assertEqual(r["instrument"], "Q Exactive Plus")
            self.assertEqual(r["ms2_analyzer"], "FTMS")
            r = da.classify(self.write(d, [None, "ITMS + c NSI r d Full ms2 400.0@cid35.00 [120-900]"]))
            self.assertEqual(r["ms2_analyzer"], "ITMS")
            r = da.classify(self.write(d, ["FTMS + c NSI d Full ms2 400.0@hcd30.00 [120-900]",
                                           "ITMS + c NSI r d Full ms2 401.0@cid35.00 [120-900]"]))
            self.assertEqual(r["ms2_analyzer"], "mixed")

    def test_no_filter_string_no_analyzer_claim(self):
        self.assertIsNone(da.mzml_ms2_analyzer([]))
        self.assertIsNone(da.mzml_ms2_analyzer({"TOF"}))


class ContaminantDigestWarning(unittest.TestCase):
    def meta(self, d, **extra):
        p = os.path.join(d, "search.fasta.meta.json")
        json.dump(dict({"min_unique_peptides": 2, "n_contaminants_in_set": 12,
                        "contaminant_digest": {"cut": "K*,R*", "missed_cleavages": 1,
                                               "min_pep_len": 7, "max_pep_len": 30}}, **extra),
                  open(p, "w"))
        return p

    def cfg(self, d, **enz):
        import estimate_params as ep
        p = os.path.join(d, "sage.json")
        json.dump({"database": {"enzyme": dict(ep.DDA_DIGEST, **enz)}}, open(p, "w"))
        return p

    def test_a_staff_fasta_with_its_own_contaminants_gets_no_warning(self):
        with tempfile.TemporaryDirectory() as d:
            side = self.meta(d, n_contaminants_in_set=0,
                             contaminant_set="already_in_supplied_database")
            self.assertIsNone(rs.warn_contaminant_digest("sage", self.cfg(d, cleave_at="KRE"),
                                                         os.path.join(d, "x.fasta"), side))

    def test_the_supplied_set_marker_alone_is_enough(self):
        with tempfile.TemporaryDirectory() as d:
            side = self.meta(d, contaminant_set="already_in_supplied_database")
            self.assertIsNone(rs.warn_contaminant_digest("sage", self.cfg(d, cleave_at="KRE"),
                                                         os.path.join(d, "x.fasta"), side))

    def test_a_judged_set_with_a_real_difference_still_warns(self):
        with tempfile.TemporaryDirectory() as d:
            side = self.meta(d)
            msg = rs.warn_contaminant_digest("sage", self.cfg(d, cleave_at="KRE"),
                                             os.path.join(d, "x.fasta"), side)
            self.assertIn("WARNING", msg or "")

    def test_semi_specific_default_trypsin_is_only_the_known_info(self):
        with tempfile.TemporaryDirectory() as d:
            side = self.meta(d)
            msg = rs.warn_contaminant_digest("sage", self.cfg(d, semi_enzymatic=True),
                                             os.path.join(d, "x.fasta"), side)
            self.assertIn("INFO", msg or "")
            self.assertNotIn("WARNING", msg or "")


class FragPipeSemiDigest(unittest.TestCase):
    def test_semi_specific_fragpipe_preset_is_only_the_known_info(self):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "fp.workflow")
            p = run("make_presets.py", "--engine", "fragpipe", "--acquisition", "DDA",
                    "--instrument", "Orbitrap Fusion", "--ms2-analyzer", "FTMS",
                    "--semi-specific", "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            side = ContaminantDigestWarning().meta(d)
            msg = rs.warn_contaminant_digest("fragpipe", out, os.path.join(d, "x.fasta"), side)
            self.assertIn("INFO", msg or "")
            self.assertNotIn("WARNING", msg or "")


class IdsOnlyWording(unittest.TestCase):
    def test_one_verb(self):
        import fetch_fasta as ff
        self.assertEqual(ff.own_accession_verb(True), "identified")
        self.assertEqual(ff.own_accession_verb(False), "quantified")
        self.assertEqual(ff.own_accession_verb(None), "quantified")

    def test_the_build_note_no_longer_claims_quantification(self):
        src = open(os.path.join(SCRIPTS, "fetch_fasta.py")).read()
        self.assertNotIn("They are now \"\n            f\"quantified", src)
        self.assertNotIn("quantified as {org} proteins -- and so count", src)


class ChargeRangeTag(unittest.TestCase):
    def search(self, d, *more):
        out = os.path.join(d, "params.cfg")
        p = run("estimate_params.py", "--engine", "diann", "--acquisition", "DIA", "--instrument",
                "timsTOF HT", "--precursor-mz-range", "300", "1200", "--out", out, *more)
        self.assertEqual(p.returncode, 0, p.stderr)
        prov = os.path.join(d, "search_provenance.json")
        json.dump({"engine": "diann", "version": "x", "params_file": out, "result": {}},
                  open(prov, "w"))
        return mm.search_record(out, prov)

    def test_default_charge_range_is_tagged_in_prose_and_it_is_one_place(self):
        with tempfile.TemporaryDirectory() as d:
            rec = self.search(d)
            self.assertIn(f"precursor charge 2–4 ({mm.DEF_INNER})", mm.search_paragraph(rec))
            self.assertEqual(rec["pr_charge"]["default"], mm.DEF_INNER)

    def test_a_charge_range_the_user_gave_is_not_tagged(self):
        with tempfile.TemporaryDirectory() as d:
            rec = self.search(d, "--overrides", '{"--min-pr-charge": 1, "--max-pr-charge": 5}')
            self.assertIsNone(rec["pr_charge"].get("default"))
            self.assertNotIn(f"precursor charge 1–5 ({mm.DEF_INNER})", mm.search_paragraph(rec))


if __name__ == "__main__":
    unittest.main()
