#!/usr/bin/env python3
"""
Ruling 120's conditions for the non-specific Orbitrap DDA route (lane X3). The routing itself is
in test_engine_choice.py (RuleRows); here:

  1. Protein FDR. FragPipe's Nonspecific workflows ship Philosopher's filter `--sequential --prot 1`:
     no protein-level FDR. The skill KEEPS it, because Philosopher's protein FDR also removes the
     PSMs of every protein that fails it (philosopher lib/fil/fil.go extractPSMfromPepXML, in both
     the --sequential and the 2D filter), which would cut the endogenous peptides a peptidome
     search is for. So it is disclosed instead: the preset's record (the manifest) and
     fragpipe_identifications.json say `protein_fdr` applied false and why; the Methods sentence
     says the protein groups are not FDR-controlled; the analysis brief asks for no protein-group
     count. One reading: fragpipe_params.protein_fdr_record.
  3. Contaminants. Philosopher flags only its own contam_ entries; the skill's Cont_ entries
     (porcine pepsin was 23.5% of the PSMs of a peptidome search) passed unflagged. The FragPipe
     identifications adapter applies the skill's ONE classifier (fetch_fasta.
     is_contaminant_accession, on protein_ids' accessions; a PSM is a contaminant when any protein
     it maps to is one, as contaminants.R treats a group) and reports the share separately.
  4. Sage stays reachable: a Sage search's confident matches are checked against its precursor
     window (sage_lfq_check.precursor_window_check). Confident matches at a median +8.5 ppm in a
     +/-10 ppm window lost most true matches (pep22); the check warns with the number, writes the
     config with the window centred on the measured offset, and the re-run command.
Hermetic; stdlib + pyarrow (the Sage tests skip without it).
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
PRESETS = os.path.join(os.path.dirname(HERE), "presets")
sys.path.insert(0, SCRIPTS)

import fragpipe_params as fpp  # noqa: E402
import run_search  # noqa: E402
import sage_lfq_check as slc  # noqa: E402

PEPTIDOME = os.path.join(PRESETS, "fragpipe", "Nonspecific-peptidome.workflow")
LFQ_MBR = os.path.join(PRESETS, "fragpipe", "LFQ-MBR.workflow")
QE = "Q Exactive Plus"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text if isinstance(text, str) else json.dumps(text))
    return path


class ProteinFdrReading(unittest.TestCase):
    def test_prot_1_is_no_protein_fdr_and_says_why(self):
        r = fpp.protein_fdr_record("--sequential --prot 1")
        self.assertEqual((r["applied"], r["threshold"]), (False, None))
        self.assertIn("removes the PSMs", r["why"])
        self.assertIn("not FDR-controlled", r["protein_filter"])
        self.assertIn("no protein count", r["protein_filter"])

    def test_a_protein_threshold_is_applied(self):
        r = fpp.protein_fdr_record("--sequential --prot 0.01 --picked")
        self.assertEqual((r["applied"], r["threshold"]), (True, 0.01))
        self.assertIn("protein 1% FDR", r["protein_filter"])

    def test_an_unset_prot_is_philosophers_default(self):
        r = fpp.protein_fdr_record("--sequential --razor")
        self.assertEqual((r["applied"], r["threshold"]), (True, 0.01))
        self.assertIn("Philosopher's default", r["protein_filter"])

    def test_no_filter_recorded_is_unknown(self):
        r = fpp.protein_fdr_record(None)
        self.assertIsNone(r["applied"])
        self.assertIn("not recorded", r["protein_filter"])

    def test_the_methods_sentence_says_the_groups_are_not_controlled(self):
        s = fpp.philosopher_methods_sentence("--sequential --prot 1")
        self.assertIn("no protein-level FDR filter (--prot 1)", s)
        self.assertIn("not FDR-controlled", s)
        self.assertNotIn("not FDR-controlled",
                         fpp.philosopher_methods_sentence("--sequential --prot 0.01 --picked"))


class PresetRecordsTheChoice(unittest.TestCase):
    """make_presets.py records the protein-FDR choice of the workflow it wrote (the manifest's
    preset_provenance): self-describing, read from the written file (rule 1)."""

    def preset(self, *args):
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "fp.workflow")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_presets.py"),
                                "--engine", "fragpipe", "--acquisition", "DDA", "--instrument",
                                QE, "--out", out, *args], capture_output=True, text=True,
                               timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out + ".provenance.json") as fh:
                prov = json.load(fh)
            props, _ = fpp.read_workflow(out)
        return prov, props

    def test_the_peptidome_workflow_keeps_the_vendor_filter_and_records_no_protein_fdr(self):
        prov, props = self.preset("--enzyme", "none", "--ids-only")
        self.assertEqual(props[fpp.PHILOSOPHER_FILTER_KEY], "--sequential --prot 1")
        rec = prov["protein_fdr"]
        self.assertEqual((rec["applied"], rec["filter"]), (False, "--sequential --prot 1"))
        self.assertIn("fil.go", rec["why"])

    def test_the_lfq_workflow_records_its_protein_fdr(self):
        prov, props = self.preset()
        self.assertEqual((prov["protein_fdr"]["applied"], prov["protein_fdr"]["threshold"]),
                         (True, 0.01))


PSM_HEADER = ("Spectrum\tPeptide\tCharge\tProbability\tAssigned Modifications\tProtein\t"
              "Protein ID\tMapped Proteins\n")
PSMS = [
    # the skill's Cont_ entry (porcine pepsin), twice
    ("r1.00001.00001.2", "PEPSINPEP", "2", "sp|Cont_P00791|PEPA_PIG", "Cont_P00791", ""),
    ("r1.00002.00002.3", "PEPSINPEP", "3", "sp|Cont_P00791|PEPA_PIG", "Cont_P00791", ""),
    # a sample protein
    ("r1.00003.00003.2", "CASEINPEP", "2", "sp|P02662|CASA1_BOVIN", "P02662", ""),
    # a sample protein's peptide that a Cont_ entry also holds: a contaminant (any member)
    ("r2.00004.00004.2", "SHAREDPEP", "2", "sp|P02662|CASA1_BOVIN", "P02662",
     "sp|Cont_P02769|ALBU_BOVIN"),
    # FragPipe's own contaminant tag, which its Protein ID drops
    ("r2.00005.00005.2", "KERATINPEP", "2", "contam_sp|P04264|K2C1_HUMAN", "P04264", ""),
    ("r2.00006.00006.2", "LACTOPEP", "2", "sp|P02754|LACB_BOVIN", "P02754", ""),
]


class FragPipeIdentifications(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.out = self._t.name
        write(os.path.join(self.out, "psm.tsv"), PSM_HEADER + "".join(
            f"{s}\t{p}\t{c}\t0.99\t\t{prot}\t{pid}\t{m}\n" for s, p, c, prot, pid, m in PSMS))
        write(os.path.join(self.out, "peptide.tsv"),
              "Peptide\tProbability\n" + "".join(f"{p}\t0.99\n" for p in
                                                 sorted({x[1] for x in PSMS})))
        shutil.copyfile(PEPTIDOME, os.path.join(self.out, fpp.SAVED_WORKFLOW))
        self.path = run_search.adapt_fragpipe_ids(
            self.out, fpp.workflow_file_quantification(PEPTIDOME))
        with open(os.path.join(self.out, run_search.FRAGPIPE_IDS_RECORD)) as fh:
            self.rec = json.load(fh)

    def tearDown(self):
        self._t.cleanup()

    def rows(self):
        with open(self.path, newline="") as fh:
            return list(csv.DictReader(fh, delimiter="\t"))

    def test_each_row_says_whether_it_is_a_contaminant(self):
        flag = {r["peptide"]: r["contaminant"] for r in self.rows()}
        self.assertEqual(flag, {"PEPSINPEP": "yes", "CASEINPEP": "no", "SHAREDPEP": "yes",
                                "KERATINPEP": "yes", "LACTOPEP": "no"})

    def test_the_contaminant_share_is_reported_separately(self):
        c = self.rec["contaminants"]
        self.assertEqual((c["psms"], c["contaminant_psms"]), (6, 4))
        self.assertEqual((c["peptides"], c["contaminant_peptides"]), (5, 3))
        self.assertAlmostEqual(c["contaminant_psm_share"], 4 / 6, places=4)
        self.assertEqual(c["top"][0], {"accession": "Cont_P00791", "psms": 2})
        self.assertIn("fetch_fasta.is_contaminant_accession", c["rule"])
        s = self.rec["contaminant_summary"]
        self.assertIn("4 of 6 PSMs (66.7%)", s)
        self.assertIn("3 of 5 distinct peptides (60.0%)", s)
        self.assertIn("Cont_P00791 (2 PSMs)", s)
        self.assertEqual(c["per_run"], {"r1": {"psms": 3, "contaminant_psms": 2},
                                        "r2": {"psms": 3, "contaminant_psms": 2}})

    def test_the_record_says_the_protein_groups_have_no_protein_fdr(self):
        self.assertEqual(self.rec["protein_fdr"]["applied"], False)
        self.assertIn("not FDR-controlled", self.rec["protein_filter"])

    def test_without_contaminants_the_share_is_zero_not_missing(self):
        with tempfile.TemporaryDirectory() as out:
            write(os.path.join(out, "psm.tsv"), PSM_HEADER +
                  "r1.00003.00003.2\tCASEINPEP\t2\t0.99\t\tsp|P02662|CASA1_BOVIN\tP02662\t\n")
            run_search.adapt_fragpipe_ids(out, fpp.workflow_file_quantification(PEPTIDOME))
            with open(os.path.join(out, run_search.FRAGPIPE_IDS_RECORD)) as fh:
                rec = json.load(fh)
        self.assertEqual((rec["contaminants"]["contaminant_psms"],
                          rec["contaminants"]["contaminant_psm_share"]), (0, 0.0))
        self.assertIn("0 of 1 PSMs (0.0%)", rec["contaminant_summary"])


class TheBriefAndThePage(unittest.TestCase):
    """The analysis brief and the report page read the record: no protein-group count when the
    groups have no protein FDR, and the contaminant share as recorded."""

    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.s = os.path.join(self._t.name, "session")
        os.makedirs(os.path.join(self.s, "output", "tables"))
        self.cond = write(os.path.join(self.s, "input", "conditions.csv"),
                          "File.Name,Group\nr1,a\nr2,b\n")
        d = os.path.join(self.s, "output", "search_fragpipe")
        os.makedirs(d)
        for name, text in (("psm.tsv", PSM_HEADER + "".join(
                f"{s}\t{p}\t{c}\t0.99\t\t{prot}\t{pid}\t{m}\n" for s, p, c, prot, pid, m in PSMS)),):
            write(os.path.join(d, name), text)
        shutil.copyfile(PEPTIDOME, os.path.join(d, fpp.SAVED_WORKFLOW))
        write(os.path.join(d, "search_provenance.json"),
              {"engine": "fragpipe", "version": "24.0", "quantified": False})
        run_search.adapt_fragpipe_ids(d, fpp.workflow_file_quantification(PEPTIDOME))
        self.d = d

    def tearDown(self):
        self._t.cleanup()

    def brief(self):
        out = os.path.join(self._t.name, "ANALYSIS_PROMPT.md")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "analysis_prompt.py"),
                            "--de-dir", os.path.join(self.s, "output", "tables"), "--out", out,
                            "--conditions", self.cond, "--session", self.s, "--acquisition",
                            "DDA"], capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out, encoding="utf-8") as fh:
            return " ".join(fh.read().split())

    def test_the_brief_asks_for_no_protein_group_count(self):
        b = self.brief()
        self.assertIn("no protein-level FDR", b)
        self.assertIn("never a protein-group count", b)
        self.assertNotIn("distinct peptides and protein groups", b)

    def test_the_brief_carries_the_contaminant_share(self):
        b = self.brief()
        self.assertIn("4 of 6 PSMs (66.7%)", b)
        self.assertIn("`contaminant` column", b)

    def test_the_page_states_the_contaminant_share(self):
        import make_analysis_html as mah
        note = mah.identifications_only_note({}, self.s)
        self.assertIn("4 of 6 PSMs (66.7%)", note["text"])
        self.assertIn("not FDR-controlled", note["text"])


def sage_results(path, offsets, window_q=0.001):
    """A results.sage.parquet with one confident rank-1 target PSM per signed offset (ppm)."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    n = len(offsets)
    calc = [1000.0 + i for i in range(n)]
    exp = [c * (1 + o * 1e-6) for c, o in zip(calc, offsets)]
    pq.write_table(pa.table({
        "filename": [f"run{i % 3}.mzML" for i in range(n)],
        "peptide": [f"PEPTIDE{i}K" for i in range(n)],
        "proteins": ["sp|P02662|CASA1_BOVIN"] * n,
        "charge": [2] * n, "rank": [1] * n, "is_decoy": [False] * n,
        "spectrum_q": [window_q] * n, "peptide_q": [window_q] * n, "protein_q": [window_q] * n,
        "precursor_ppm": [abs(o) for o in offsets],
        "expmass": exp, "calcmass": calc, "isotope_error": [0.0] * n}), path)


class SagePrecursorWindow(unittest.TestCase):
    def setUp(self):
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("pyarrow")
        self._t = tempfile.TemporaryDirectory()
        self.out = self._t.name
        self.cfg = write(os.path.join(self.out, "sage_config.json"),
                         {"database": {"fasta": "x.fasta"},
                          "precursor_tol": {"ppm": [-10.0, 10.0]},
                          "fragment_tol": {"ppm": [-10.0, 10.0]}, "quant": {"lfq": False}})
        write(os.path.join(self.out, "search_provenance.json"),
              {"engine": "sage", "params_file": self.cfg, "tools": "/t/tools.json",
               "bundle": "/b/wf", "fasta": "/f/db.fasta", "files": ["/r/run0.mzML"]})

    def tearDown(self):
        self._t.cleanup()

    def offsets(self, centre, n=310):
        # spread evenly +/- 1.5 ppm around the centre (ten whole cycles), cut at the window edge
        return [o for o in (centre + (i % 31 - 15) / 10.0 for i in range(n)) if -10 <= o <= 10]

    def test_matches_at_the_window_edge_warn_with_the_number_and_a_shifted_window(self):
        sage_results(os.path.join(self.out, "results.sage.parquet"), self.offsets(8.5))
        rec = slc.precursor_window_check(self.out)
        self.assertEqual(rec["status"], "warn")
        self.assertEqual(rec["window_ppm"], [-10.0, 10.0])
        self.assertAlmostEqual(rec["median_offset_ppm"], 8.5, delta=0.3)
        # Sage's window is calculated minus observed: +8.5 measured sits at -8.5 in it, so the
        # centred window runs to -19 (PM check b: [-2, 19] gave 0 PSMs at 1%, [-19, 2] 4,034)
        self.assertEqual(rec["suggested_window_ppm"], [-19.0, 2.0])
        self.assertAlmostEqual(rec["median_in_window_ppm"], -8.5, delta=0.3)
        self.assertLess(rec["position_in_window"], -0.5)
        self.assertIn("lower edge", rec["message"])
        self.assertIn("[-19, 2] ppm", rec["message"])
        self.assertIn("precursor_-19_2ppm", rec["corrected_params"])
        self.assertIn("+8.", rec["message"])
        self.assertIn("+/-10 ppm", rec["message"])
        with open(rec["corrected_params"]) as fh:
            cfg = json.load(fh)
        self.assertEqual(cfg["precursor_tol"], {"ppm": [-19.0, 2.0]})
        self.assertEqual(cfg["fragment_tol"], {"ppm": [-10.0, 10.0]})       # nothing else changed
        self.assertIn("--params " + rec["corrected_params"], rec["rerun"]["command"])
        self.assertIn("--engine sage", rec["rerun"]["command"])

    def test_the_warning_never_says_fragpipe_calibrates(self):
        """C2 proteomics H2: the shipped Nonspecific-peptidome workflow has calibrate_mass=0, so
        the advice names FragPipe's wider window, never a calibration it does not do."""
        self.assertEqual(fpp.read_workflow(PEPTIDOME)[0].get("msfragger.calibrate_mass"), "0")
        sage_results(os.path.join(self.out, "results.sage.parquet"), self.offsets(8.5))
        msg = slc.precursor_window_check(self.out)["message"]
        self.assertNotIn("which calibrates", msg)
        self.assertIn("+/-20 ppm", msg)
        self.assertIn("no mass calibration either", msg)
        ref = os.path.join(os.path.dirname(HERE), "references", "search-engines.md")
        with open(ref, encoding="utf-8") as fh:
            self.assertNotIn("FragPipe, which calibrates", " ".join(fh.read().split()))

    def test_a_centred_offset_is_ok(self):
        sage_results(os.path.join(self.out, "results.sage.parquet"), self.offsets(1.0))
        rec = slc.precursor_window_check(self.out)
        self.assertEqual(rec["status"], "ok")
        self.assertNotIn("corrected_params", rec)

    def test_just_inside_the_outer_half_is_still_ok(self):
        """mutant guard on the edge rule: a median 4.5 ppm into a +/-10 window is inside half"""
        sage_results(os.path.join(self.out, "results.sage.parquet"), self.offsets(4.5))
        self.assertEqual(slc.precursor_window_check(self.out)["status"], "ok")

    def test_a_negative_offset_warns_too(self):
        sage_results(os.path.join(self.out, "results.sage.parquet"), self.offsets(-7.0))
        rec = slc.precursor_window_check(self.out)
        self.assertEqual(rec["status"], "warn")
        self.assertEqual(rec["suggested_window_ppm"], [-3.0, 17.0])

    def test_an_asymmetric_window_reads_the_position_in_sages_sign(self):
        """A window already shifted the right way ([-19, 2]) for +8.5 ppm measured fits: the
        matches sit at -8.5 in it, its centre. The old sign read them as far outside."""
        cfg = json.load(open(self.cfg))
        write(self.cfg, dict(cfg, precursor_tol={"ppm": [-19.0, 2.0]}))
        sage_results(os.path.join(self.out, "results.sage.parquet"),
                     [o for o in self.offsets(8.5) if o <= 19])
        rec = slc.precursor_window_check(self.out)
        self.assertEqual(rec["status"], "ok", rec.get("message"))
        self.assertLess(abs(rec["position_in_window"]), 0.2)

    def test_the_sign_conversion_is_one_function(self):
        self.assertEqual(slc.sage_window_ppm(8.4), -8.4)

    def test_too_few_matches_is_unchecked_and_says_so(self):
        sage_results(os.path.join(self.out, "results.sage.parquet"), [8.5] * 10)
        rec = slc.precursor_window_check(self.out)
        self.assertEqual(rec["status"], "unchecked")
        self.assertIn("10", rec["reason"])

    def test_the_ids_adapter_runs_it_and_records_it(self):
        sage_results(os.path.join(self.out, "results.sage.parquet"), self.offsets(8.5))
        p = subprocess.run([sys.executable, "-c",
                            "import sys; sys.path.insert(0, sys.argv[1]); import run_search; "
                            "run_search.adapt_sage_ids(sys.argv[2])", SCRIPTS, self.out],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("WARNING", p.stderr)
        self.assertIn("precursor window", p.stderr)
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["precursor_window_check"]["status"], "warn")
        self.assertTrue(os.path.exists(os.path.join(self.out, slc.WINDOW_RECORD)))


class JobSizedFromTheDatabase(unittest.TestCase):
    """Ruling 120 condition 2 (lane X2's helper): the FragPipe peptidome route's job memory is
    run_search.broad_digest_mem_gb("fragpipe", <the database's sequences>) and FragPipe's Java heap
    follows it (fragpipe_ram_gb). The route is sized where every single-job search is
    (single_job_sizing, which X2 changed in place), so this lane adds no second call; the test pins
    that the Ruling 120 preset reaches it. Guarded with hasattr: before X2's branch is merged the
    helper does not exist here and the test is skipped (said, never passed vacuously)."""

    @unittest.skipUnless(hasattr(run_search, "broad_digest_mem_gb"),
                         "lane X2's run_search.broad_digest_mem_gb is not on this branch yet")
    def test_the_peptidome_route_is_sized_from_the_database(self):
        sys.path.insert(0, HERE)
        from test_search_job_sizing import search, header, provenance
        n_seq = 22939
        with tempfile.TemporaryDirectory() as d:
            wf = os.path.join(d, "preset.workflow")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_presets.py"),
                                "--engine", "fragpipe", "--acquisition", "DDA", "--instrument",
                                QE, "--enzyme", "none", "--ids-only", "--out", wf],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            write(os.path.join(d, "db.fasta.meta.json"), {"n_sequences": n_seq})
            job = os.path.join(d, "fp_job.sh")
            with open(wf) as fh:
                p, out = search(d, "fragpipe", "--sbatch", job, cfg=fh.read())
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            want = run_search.broad_digest_mem_gb("fragpipe", n_seq)
            self.assertIsInstance(want, int)
            self.assertEqual(header(job, "mem"), f"{want}G")
            with open(job) as fh:
                self.assertIn(f"--ram {run_search.fragpipe_ram_gb(want)}", fh.read())
            self.assertEqual(provenance(out)["job_n_sequences"], n_seq)


if __name__ == "__main__":
    unittest.main()
