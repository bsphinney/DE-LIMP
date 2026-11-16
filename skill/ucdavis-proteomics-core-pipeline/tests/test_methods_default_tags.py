#!/usr/bin/env python3
"""
The Methods tag every default they print (DE-LIMP rule 2; 2.15 acceptance check 5, all five pairs).

make_methods.py printed fixed Carbamidomethyl (C) untagged, in the prose and in the parameter
table, although the params' own sidecar (estimate_params.py `cys_alkylation`) says "DEFAULT, not
user-confirmed". A DDA search also printed Ox(M), the protein N-terminal acetyl and the instrument
tolerances untagged. The tag now comes from the record that holds the value (the sidecar; the
FragPipe preset's provenance) -- make_methods.py decides nothing about what is a default, and
prints the same tag the LC defaults use. A value the user GAVE carries none.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
FIX = os.path.join(HERE, "fixtures", "fragpipe24_nonspecific_peptidome")
sys.path.insert(0, SCRIPTS)

import make_methods as mm  # noqa: E402

TAG = "[facility default — confirm]"


def estimate(d, engine, acq, *more):
    out = os.path.join(d, "params.cfg" if engine == "diann" else "sage_config.json")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                        "--engine", engine, "--acquisition", acq, "--instrument", "timsTOF HT",
                        "--precursor-mz-range", "300", "1200", "--out", out, *more],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError(p.stderr)
    return out


def record(d, engine, params):
    prov = os.path.join(d, "search_provenance.json")
    with open(prov, "w") as fh:
        json.dump({"engine": engine, "version": "x", "params_file": params, "result": {}}, fh)
    rec = mm.search_record(params, prov)
    return rec, mm.search_paragraph(rec)


def table_rows(rec):
    """The modification and tolerance rows of the parameter table, as make_methods writes them."""
    rows = [(f"{m['type'].capitalize()} modification", mm.mod_phrase(m)) for m in rec["mods"]]
    for key in ("ms1_tol", "ms2_tol"):
        t = rec.get(key)
        if t:
            rows.append((key, f"{mm._g(t['value'])} {t['unit']}"
                         + (f" ({t['default']})" if t.get("default") else "")))
    return rows


class DiannCarbamidomethyl(unittest.TestCase):
    def test_unasked_cam_is_tagged_in_prose_and_table(self):
        with tempfile.TemporaryDirectory() as d:
            rec, text = record(d, "diann", estimate(d, "diann", "DIA"))
            self.assertIn(f"Carbamidomethyl (C) {TAG}", text)
            self.assertIn(f"Carbamidomethyl (C) {TAG}",
                          [v for k, v in table_rows(rec) if k == "Fixed modification"])

    def test_given_cam_is_not_tagged(self):
        with tempfile.TemporaryDirectory() as d:
            _, text = record(d, "diann", estimate(d, "diann", "DIA", "--cys-alkylation", "cam"))
            self.assertIn("Carbamidomethyl (C)", text)
            self.assertNotIn(f"Carbamidomethyl (C) {TAG}", text)


class SageDda(unittest.TestCase):
    def test_default_ox_nterm_acetyl_cam_and_tolerances_are_tagged(self):
        with tempfile.TemporaryDirectory() as d:
            rec, text = record(d, "sage", estimate(d, "sage", "DDA"))
            for phrase in (f"Carbamidomethyl (C) {TAG}", f"Oxidation (M) {TAG}",
                           f"Acetyl (Protein N-term) {TAG}"):
                self.assertIn(phrase, text)
            self.assertRegex(text, r"precursor \(MS1\) ±[\d.]+ ppm \(facility default — confirm\)")
            self.assertRegex(text, r"fragment \(MS2\) ±[\d.]+ ppm \(facility default — confirm\)")
            rows = [v for _, v in table_rows(rec)]
            self.assertTrue(all(TAG in v or "ppm (facility default" in v for v in rows), rows)

    def test_given_choices_are_not_tagged(self):
        with tempfile.TemporaryDirectory() as d:
            _, text = record(d, "sage", estimate(d, "sage", "DDA", "--cys-alkylation", "cam",
                                                 "--var-mods", "ox"))
            self.assertIn("Carbamidomethyl (C)", text)
            self.assertNotIn(f"Carbamidomethyl (C) {TAG}", text)
            self.assertNotIn(f"Oxidation (M) {TAG}", text)
            # the protein N-terminal acetyl is added whatever was asked: still a default
            self.assertIn(f"Acetyl (Protein N-term) {TAG}", text)

    def test_no_record_no_claim(self):
        with tempfile.TemporaryDirectory() as d:
            params = estimate(d, "sage", "DDA")
            os.remove(params + ".rationale.json")
            _, text = record(d, "sage", params)
            self.assertNotIn("facility default", text)


class FragPipeDda(unittest.TestCase):
    DEFAULTS = {"cys_alkylation": {"value": "cam", "given": False},
                "var_mods": {"words": ["ox"], "source": "DEFAULT, not user-confirmed: Ox(M) on "
                             "by default for DDA; + protein N-term acetyl (always)"},
                "tolerances": {"mode": "vendor", "source": "the vendor workflow's values kept"}}

    def search(self, d, preset):
        shutil.copy(os.path.join(FIX, "fragger.params"), os.path.join(d, "fragger.params"))
        params = os.path.join(d, "fragger.params")
        with open(params + ".provenance.json", "w") as fh:
            json.dump(dict({"engine": "fragpipe"}, **preset), fh)
        return record(d, "fragpipe", params)

    def test_defaults_are_tagged(self):
        with tempfile.TemporaryDirectory() as d:
            _, text = self.search(d, self.DEFAULTS)
            for phrase in (f"Carbamidomethyl (C) {TAG}", f"Oxidation (M) {TAG}",
                           f"Acetyl (Protein N-term) {TAG}"):
                self.assertIn(phrase, text)
            self.assertIn("(facility default — confirm)", text)

    def test_given_values_are_not_tagged(self):
        with tempfile.TemporaryDirectory() as d:
            given = {"cys_alkylation": {"value": "cam", "given": True},
                     "var_mods": {"words": ["ox"], "source": "user requested (--var-mods ox); "
                                  "+ protein N-term acetyl (always)"},
                     "tolerances": {"mode": "caller", "source": "--ms2-ppm 10"}}
            _, text = self.search(d, given)
            self.assertNotIn(f"Carbamidomethyl (C) {TAG}", text)
            self.assertNotIn(f"Oxidation (M) {TAG}", text)
            self.assertNotIn("ppm (facility default", text)
            self.assertIn(f"Acetyl (Protein N-term) {TAG}", text)


class SessionDefaults(unittest.TestCase):
    """make_methods.py --session: where it writes, what it reads, what it says of an
    identifications-only search (2.15 acceptance: pep24-215fp F5, pep22 F6, pxd4 F8)."""

    def make(self, d, ids_only=False):
        for sub in ("input/wf", "output", "raw"):
            os.makedirs(os.path.join(d, sub))
        meta = {"fasta": "x.fasta", "organism": "Bos taurus", "taxid": 9913,
                "proteome": "UP000009136", "proteome_type": "Reference proteome",
                "content_used": "one_per_gene", "uniprot_release": "2026_03", "n_proteome": 22939,
                "n_contaminants_appended": 381, "contaminant_set": "universal",
                "n_contaminants_dropped_as_target": 3}
        for rel, obj in (("input/search.fasta.meta.json", meta),
                         ("input/acquisition.json", {
                             "overall": "DDA", "files": [
                                 {"file": "a.raw", "instrument": "Q Exactive Plus Orbitrap"},
                                 {"file": "b.raw", "instrument": "Q Exactive Plus Orbitrap"}]}),
                         ("input/wf/workflow.manifest.json", {"acquisition": "DDA"})):
            with open(os.path.join(d, rel), "w") as fh:
                json.dump(obj, fh)
        raws = []
        for n in ("PEP_a.raw", "PEP_b.raw"):
            r = os.path.join(d, "raw", n)
            with open(r, "wb") as fh:
                fh.write(b"\0" * 64)
            raws.append(r)
        params = estimate(d, "sage", "DDA", *(["--ids-only"] if ids_only else []))
        prov = os.path.join(d, "search_provenance.json")
        with open(prov, "w") as fh:
            json.dump({"engine": "sage", "version": "0.14.7", "params_file": params}, fh)
        return raws, prov

    def run_methods(self, d, raws, prov, *more):
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"), "--raw",
                            *raws, "--session", d, "--search-prov", prov, *more],
                           capture_output=True, text=True, cwd=tempfile.gettempdir())
        self.assertEqual(r.returncode, 0, r.stderr)
        return r

    def test_session_default_out_fasta_meta_and_instrument(self):
        with tempfile.TemporaryDirectory() as d:
            raws, prov = self.make(d)
            self.run_methods(d, raws, prov)
            out = os.path.join(d, "output", "methods.md")
            self.assertTrue(os.path.isfile(out), "methods.md is written into the session")
            text = open(out, encoding="utf-8").read()
            self.assertIn("Bos taurus", text)                    # --fasta-meta read from the session
            self.assertNotIn("Spectra were searched against ____", text)
            self.assertIn("Q Exactive Plus Orbitrap", text)      # the instrument the session knows
            self.assertNotIn("facility filename prefix", text)

    def test_ids_only_does_not_say_quantified(self):
        with tempfile.TemporaryDirectory() as d:
            raws, prov = self.make(d, ids_only=True)
            self.run_methods(d, raws, prov)
            text = open(os.path.join(d, "output", "methods.md"), encoding="utf-8").read()
            self.assertIn("removed from the library first", text)
            self.assertNotIn("are quantified under their own accessions", text)
            self.assertIn("are identified under their own accessions", text)

    def test_quantified_search_still_says_quantified(self):
        with tempfile.TemporaryDirectory() as d:
            raws, prov = self.make(d, ids_only=False)
            self.run_methods(d, raws, prov)
            text = open(os.path.join(d, "output", "methods.md"), encoding="utf-8").read()
            self.assertIn("are quantified under their own accessions", text)


class OrganismLine(unittest.TestCase):
    """The one organism line (session.database_organisms): every species the database holds."""

    def test_added_second_species_is_named(self):
        import session as sess
        fm = {"organism": "Homo sapiens", "taxid": 9606, "added_sequences": [{"entries": [
            {"accession": "P00330", "name": "ADH1_YEAST",
             "header": "sp|P00330|ADH1_YEAST Alcohol dehydrogenase 1 OS=Saccharomyces cerevisiae "
                       "(strain ATCC 204508 / S288c) OX=559292 GN=ADH1 PE=1 SV=2"},
            {"accession": "P00331", "name": "ADH2_YEAST",
             "header": "sp|P00331|ADH2_YEAST Alcohol dehydrogenase 2 OS=Saccharomyces cerevisiae "
                       "(strain ATCC 204508 / S288c) OX=559292 GN=ADH2 PE=1 SV=2"}]}]}
        line = sess.database_organisms(fm)
        self.assertIn("Homo sapiens (taxid 9606)", line)
        self.assertIn("Saccharomyces cerevisiae (strain ATCC 204508 / S288c) (taxid 559292", line)
        self.assertEqual(line.count("Saccharomyces"), 1)          # once, not per entry

    def test_one_species_and_none(self):
        import session as sess
        self.assertEqual(sess.database_organisms({"organism": "Homo sapiens", "taxid": 9606}),
                         "Homo sapiens (taxid 9606)")
        self.assertIsNone(sess.database_organisms({}))

    def test_readme_summary_and_report_header_use_it(self):
        import session_docs
        f = {"fm": {"organism": "Homo sapiens", "taxid": 9606, "added_sequences": [{"entries": [
            {"name": "ADH1_YEAST", "header": "sp|P00330|ADH1_YEAST X OS=Saccharomyces cerevisiae "
                                             "OX=559292 GN=ADH1"}]}]},
             "organism": "Homo sapiens", "taxid": 9606}
        import session as sess
        self.assertIn("Saccharomyces cerevisiae", sess.database_organisms(f["fm"], f["organism"],
                                                                        f["taxid"]))
        src = open(os.path.join(SCRIPTS, "session_docs.py"), encoding="utf-8").read()
        html = open(os.path.join(SCRIPTS, "make_analysis_html.py"), encoding="utf-8").read()
        self.assertIn("database_organisms(", src)
        self.assertIn("database_organisms(", html)


if __name__ == "__main__":
    unittest.main()
