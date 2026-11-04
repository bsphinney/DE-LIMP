#!/usr/bin/env python3
"""
Sage searches a timsTOF .d directly for identifications, with the pinned 0.14.7 (I48,
I49; T19 of the 2.14.2 plan).

A Core staff session refused to run Sage on timsTOF DDA-PASEF .d (2026-10-29): mzml_plan() sent
every .d to msconvert, and bioconda's Linux msconvert has no Bruker reader, so the session left the
skill's path. On HIVE the pinned Sage 0.14.7 read the two .d itself: 896 peptides at 1% in 22 s,
against 753 from an 11.6 GB Windows mzML conversion (sage015_test/RESULTS.md, arm B). 0.14.7 has
no MS1 quantification on a .d, so this route is for identifications only, and it is chosen by the
cfg alone: estimate_params.py --ids-only writes "quant": {"lfq": false}.

  * an IDs-only job hands each .d to Sage as it is (no converter, no LFQ check); a .raw beside it
    is still converted with ThermoRawFileParser and an mzML goes straight in (Review Focus 3);
  * --adapt-only writes identifications.tsv from results.sage.parquet, never a report.parquet;
  * search_provenance.json says how each .d went in (sage_input), that nothing was quantified,
    and that the LFQ check does not apply;
  * a quantifying cfg still converts a .d first.

Sage, ThermoRawFileParser and msconvert are stand-ins (tests/test_sage_sbatch_conversion.py's
harness). The results.sage.parquet fixture has the columns Sage 0.14.7 wrote on HIVE for these .d
(C0 pre-read, 2026-11-01); its values are synthetic. Hermetic.
"""
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
RUN_SEARCH = os.path.join(SCRIPTS, "run_search.py")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import run_search  # noqa: E402
from test_sage_sbatch_conversion import _Harness, _exe  # noqa: E402
# self.env() is job_env() (see _Harness.env); imported so the guard can see that.
from job_env import job_env  # noqa: E402,F401

try:
    import pyarrow  # noqa: F401
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

# results.sage.parquet as Sage 0.14.7 (conda build, reports 0.14.6) wrote it for two timsTOF HT
# DDA-PASEF .d on HIVE: names and types, in order (C0_hive_reads.md, T19). "float" is float32.
SAGE_0147_RESULTS = [
    ("psm_id", "int64"), ("filename", "string"), ("scannr", "string"), ("peptide", "string"),
    ("stripped_peptide", "string"), ("proteins", "string"), ("num_proteins", "int32"),
    ("rank", "int32"), ("is_decoy", "bool"), ("expmass", "float"), ("calcmass", "float"),
    ("charge", "int32"), ("peptide_len", "int32"), ("missed_cleavages", "int32"),
    ("semi_enzymatic", "bool"), ("ms2_intensity", "float"), ("isotope_error", "float"),
    ("precursor_ppm", "float"), ("fragment_ppm", "float"), ("hyperscore", "float"),
    ("delta_next", "float"), ("delta_best", "float"), ("rt", "float"), ("aligned_rt", "float"),
    ("predicted_rt", "float"), ("delta_rt_model", "float"), ("ion_mobility", "float"),
    ("predicted_mobility", "float"), ("delta_mobility", "float"), ("matched_peaks", "int32"),
    ("longest_b", "int32"), ("longest_y", "int32"), ("longest_y_pct", "float"),
    ("matched_intensity_pct", "float"), ("scored_candidates", "int32"), ("poisson", "float"),
    ("sage_discriminant_score", "float"), ("posterior_error", "float"), ("spectrum_q", "float"),
    ("peptide_q", "float"), ("protein_q", "float"), ("reporter_ion_intensity", "list<float>"),
]

IDS_COLUMNS = ["run", "peptide", "proteins", "charge", "spectrum_q", "peptide_q", "protein_q"]


def _psm_rows(first):
    """(peptide, stripped, proteins, charge, rank, is_decoy, spectrum_q, peptide_q, protein_q)
    per PSM of one run (`first`: the first run): three confident targets -- one seen in two
    spectra, one at two charges -- a target above 1% peptide q, a decoy that passes every q, a
    chimeric rank-2 PSM that passes every q, and a peptide whose peptide_q (one value over all
    the runs) passes while its only spectrum in the first run is poor (spectrum_q 0.30; 0.0004 in
    the other run). The last two are not identifications in their run (L7R-F2)."""
    return [
        ("PEPT[+15.9949]IDEK", "PEPTIDEK", "sp|P11111|AAA_HUMAN", 2, 1, False, 0.004, 0.002,
         0.001),
        ("PEPT[+15.9949]IDEK", "PEPTIDEK", "sp|P11111|AAA_HUMAN", 2, 1, False, 0.001, 0.002,
         0.001),
        ("SAMPLER", "SAMPLER", "sp|P22222|BBB_HUMAN;sp|P33333|CCC_HUMAN", 2, 1, False, 0.003,
         0.003, 0.002),
        ("SAMPLER", "SAMPLER", "sp|P22222|BBB_HUMAN;sp|P33333|CCC_HUMAN", 3, 1, False, 0.006,
         0.003, 0.002),
        ("WEAKPEPK", "WEAKPEPK", "sp|P44444|DDD_HUMAN", 2, 1, False, 0.009, 0.05, 0.04),
        ("KEPTIDEP", "KEPTIDEP", "rev_sp|P11111|AAA_HUMAN", 2, 1, True, 0.001, 0.001, 0.001),
        ("CHIMERAK", "CHIMERAK", "sp|P55555|EEE_HUMAN", 2, 2, False, 0.001, 0.001, 0.001),
        ("SHAREDPEPK", "SHAREDPEPK", "sp|P66666|FFF_HUMAN", 2, 1, False,
         0.30 if first else 0.0004, 0.0019, 0.001),
    ]


def write_ids_results(out, runs, quant=None):
    """Sage 0.14.7-shaped output of an IDs-only search in `out`: results.sage.parquet (every
    column in SAGE_0147_RESULTS, `filename` = each run's file name as Sage records it) and
    results.json (the parameters as run). No lfq.parquet: Sage ran with quant.lfq false."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    types = {"int64": pa.int64(), "int32": pa.int32(), "string": pa.string(), "bool": pa.bool_(),
             "float": pa.float32(), "list<float>": pa.list_(pa.float32())}
    zero = {"int64": 0, "int32": 0, "string": "", "bool": False, "float": 0.0,
            "list<float>": None}
    cols = {n: [] for n, _ in SAGE_0147_RESULTS}
    psm_id = 0
    for k, run in enumerate(runs):
        for i, (pep, stripped, prots, z, rank, dec, sq, pq_, prq) in enumerate(
                _psm_rows(k == 0)):
            psm_id += 1
            row = {n: zero[t] for n, t in SAGE_0147_RESULTS}
            row.update(psm_id=psm_id, filename=run, scannr=f"index={1000 + i}", peptide=pep,
                       stripped_peptide=stripped, proteins=prots,
                       num_proteins=prots.count(";") + 1, rank=rank, is_decoy=dec, charge=z,
                       peptide_len=len(stripped), expmass=1000.0 + i, calcmass=1000.0 + i,
                       spectrum_q=sq, peptide_q=pq_, protein_q=prq, ion_mobility=0.9)
            for n in cols:
                cols[n].append(row[n])
    os.makedirs(out, exist_ok=True)
    pq.write_table(pa.table({n: pa.array(cols[n], types[t]) for n, t in SAGE_0147_RESULTS}),
                   os.path.join(out, "results.sage.parquet"))
    with open(os.path.join(out, "results.json"), "w") as fh:
        json.dump({"version": "0.14.6", "quant": quant or {"tmt": None, "lfq": False}}, fh)


# Plays Sage 0.14.7: logs its argv and its inputs (everything after the telemetry flag), then
# writes IDs-only output naming each input by its file name, as Sage does for a .d or an mzML.
FAKE_SAGE_IDS = r'''#!{py}
import json, os, sys
sys.path.insert(0, {here!r})
import test_sage_dotd_ids as fx
args = sys.argv[1:]
out = args[args.index("-o") + 1]
inputs = args[args.index("--disable-telemetry-i-dont-want-to-improve-sage") + 1:]
with open(os.environ["FAKE_CONV_LOG"], "a") as fh:
    fh.write(json.dumps({{"tool": "sage", "argv": args, "inputs": inputs,
                         "missing": [i for i in inputs if not os.path.exists(i)]}}) + "\n")
with open(args[0]) as fh:
    quant = json.load(fh).get("quant")
fx.write_ids_results(out, [os.path.basename(i.rstrip("/")) for i in inputs], quant=quant)
print("[INFO sage] finished", file=sys.stderr)
'''


class _IdsHarness(_Harness):
    def setUp(self):
        super().setUp()
        _exe(self.sage, FAKE_SAGE_IDS.format(py=sys.executable, here=HERE))
        tims = os.path.join(self.d, "tims")
        self.dots = []
        for n in ("runA_S1-A1_1_1001.d", "runB_S1-B1_1_1002.d"):
            p = os.path.join(tims, n)
            os.makedirs(p)
            with open(os.path.join(p, "analysis.tdf_bin"), "wb") as fh:
                fh.write(b"\0" * 16)                         # never read: no analysis.tdf here
            self.dots.append(p)
        self.ids_cfg = os.path.join(self.d, "sage_config_ids.json")
        with open(self.ids_cfg, "w") as fh:
            json.dump({"database": {"fasta": "REPLACED_AT_RUNTIME.fasta"},
                       "precursor_tol": {"ppm": [-20.0, 20.0]},
                       "fragment_tol": {"ppm": [-20.0, 20.0]}, "quant": {"lfq": False}}, fh)

    def run_ids(self, *extra, files=None, env=None):
        """run_search.py with the IDs-only cfg."""
        saved, self.cfg = self.cfg, self.ids_cfg
        try:
            return self.run_search(*extra, files=files or self.dots, env=env)
        finally:
            self.cfg = saved

    def generate(self, files=None, ids=True):
        job = os.path.join(self.d, "sage_job.sh")
        # a quantifying .d needs the ProteoWizard container (test_sage_dotd_pwiz_conversion.py)
        sif = os.path.join(self.d, "pwiz.sif")
        open(sif, "w").close()
        p = (self.run_ids if ids else self.run_search)("--sbatch", job, "--pwiz-sif", sif,
                                                       files=files or self.dots)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        with open(job) as fh:
            return job, fh.read()

    def provenance(self):
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            return json.load(fh)

    def adapt_only(self, files=None):
        p = self.run_ids("--adapt-only", files=files)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        return json.loads(p.stdout[p.stdout.index("{"):])


class TheIdentificationsAreDelivered(unittest.TestCase):
    """L7R-F4: deliver took only report.* and search_provenance.json from output/search, so an
    identifications-only delivery (--mode no-de, "an identification list") left the route's one
    product behind."""

    def test_the_delivery_plan_takes_them(self):
        import core_submission
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "search"))
            for n in (run_search.SAGE_IDS, run_search.SAGE_IDS_RECORD):
                with open(os.path.join(d, "search", n), "w") as fh:
                    fh.write("x\n")
            items = {i["name"]: i for i in core_submission.plan_delivery(d)}
        for n in (run_search.SAGE_IDS, run_search.SAGE_IDS_RECORD):
            self.assertEqual(items["search/" + n]["status"], "OK", items.get("search/" + n))

    def test_a_search_without_them_skips_them(self):
        import core_submission
        with tempfile.TemporaryDirectory() as d:
            items = {i["name"]: i for i in core_submission.plan_delivery(d)}
        item = items["search/" + run_search.SAGE_IDS]
        self.assertEqual(item["status"], "SKIPPED")
        self.assertFalse(item["required"])


class LfqSettingDecidesTheRoute(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def cfg(self, obj, raw=None):
        p = os.path.join(self.d, "c.json")
        with open(p, "w") as fh:
            fh.write(raw if raw is not None else json.dumps(obj))
        return p

    def test_lfq_false_is_ids_only(self):
        self.assertFalse(run_search.sage_lfq_enabled(self.cfg({"quant": {"lfq": False}})))

    def test_lfq_true_or_unset_quantifies(self):
        self.assertTrue(run_search.sage_lfq_enabled(self.cfg({"quant": {"lfq": True}})))
        self.assertTrue(run_search.sage_lfq_enabled(self.cfg({"quant": {}})))
        self.assertTrue(run_search.sage_lfq_enabled(self.cfg({"database": {}})))

    def test_an_unreadable_cfg_quantifies(self):
        """Converting a .d that did not need it costs a conversion; searching one directly that
        was meant to be quantified gives quantities Sage 0.14.7 cannot make."""
        self.assertTrue(run_search.sage_lfq_enabled(self.cfg(None, raw="{not json")))
        self.assertTrue(run_search.sage_lfq_enabled(os.path.join(self.d, "absent.json")))

    def test_the_default_is_tagged_where_it_is_said(self):
        _, how = run_search.sage_lfq_setting(self.cfg({"quant": {}}))
        self.assertIn("DEFAULT", how)
        _, how = run_search.sage_lfq_setting(self.cfg({"quant": {"lfq": False}}))
        self.assertNotIn("DEFAULT", how)


class MixedInputsIdsOnly(_IdsHarness):
    def test_mixed_inputs_ids_only(self):
        """Review Focus 3: only the .raw is converted (ThermoRawFileParser); the .d and the mzML
        are searched as given."""
        mz = os.path.join(self.d, "raw", "ready.mzML")
        open(mz, "w").close()
        trfp = os.path.join(self.bin, "ThermoRawFileParser")
        with mock.patch.dict(os.environ, {"THERMORAWFILEPARSER": trfp}):
            searched, steps = run_search.mzml_plan([self.dots[0], self.raws[0], mz], self.out,
                                                   dotd="direct")
        self.assertEqual([s["input"] for s in steps], [self.raws[0]])
        self.assertIn(trfp, steps[0]["argv"])
        self.assertEqual(searched, [self.dots[0],
                                    os.path.join(self.out, "mzml", "HeLa_1.mzML"), mz])

    def test_a_dotd_and_a_raw_of_one_run_name_are_refused_on_the_direct_route(self):
        """(pin) run_search.py stops it first (check_report_runs.repeated_names); the direct
        mzml_plan test below is the one L7R-F8 changes."""
        twin = os.path.join(self.d, "raw", os.path.basename(self.dots[0])[:-2] + ".raw")
        with open(twin, "wb") as fh:
            fh.write(b"\0")
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_ids("--sbatch", job, files=[self.dots[0], twin])
        self.assertNotEqual(p.returncode, 0, p.stdout)
        self.assertIn("share a run name", p.stderr)
        self.assertFalse(os.path.exists(job))

    def test_mzml_plan_keys_its_duplicate_check_on_the_run_name(self):
        """L7R-F8: run1.d is searched as run1.d and run1.raw as run1.mzML -- different file
        names -- but the adapters name both runs run1, so they would merge. mzml_plan's own
        check refuses it, as on the convert route (where both become run1.mzML)."""
        twin = os.path.join(self.d, "raw", os.path.basename(self.dots[0])[:-2] + ".raw")
        trfp = os.path.join(self.bin, "ThermoRawFileParser")
        with mock.patch.dict(os.environ, {"THERMORAWFILEPARSER": trfp}):
            with self.assertRaises(SystemExit) as cm:
                run_search.mzml_plan([self.dots[0], twin], self.out, dotd="direct")
        self.assertIn("share a Sage run name", str(cm.exception))

    def test_one_run_name_rule(self):
        for path in ("/x/runA.d", "/x/runA.d/", "/x/runA.mzML", "/x/runA.raw", "runA.d"):
            self.assertEqual(run_search.sage_run_name(path), "runA", path)

    def test_the_job_converts_only_the_raw(self):
        mz = os.path.join(self.d, "raw", "ready.mzML")
        open(mz, "w").close()
        _job, text = self.generate(files=[self.dots[0], self.raws[0], mz])
        self.assertIn("-i=" + self.raws[0], text)
        self.assertNotIn("msconvert", text)
        sage = text[text.index(self.sage):]
        self.assertIn(shlex.quote(self.dots[0]), sage)
        self.assertIn(shlex.quote(mz), sage)


class IdsOnlyJob(_IdsHarness):
    def test_ids_only_passes_dotd_straight_to_sage(self):
        _job, text = self.generate()
        self.assertNotIn("msconvert", text)
        self.assertNotIn("apptainer", text)
        self.assertNotIn("converting", text)
        for d in self.dots:
            self.assertIn(shlex.quote(d), text)
        self.assertNotIn("sage_lfq_check.py", text)
        self.assertEqual(self.calls(), [], "nothing runs at generation")

    def test_provenance_says_direct_and_not_quantified(self):
        self.generate()
        prov = self.provenance()
        self.assertEqual(prov["sage_input"]["dotd"], "direct")
        self.assertIn("identifications only", prov["sage_input"]["why"])
        self.assertEqual(prov["sage_input"]["files"], self.dots)
        self.assertIs(prov["quantified"], False)
        # L7R-F13: the one not-applicable record, sage_lfq_check's
        import sage_lfq_check as slc
        rec = {k: v for k, v in prov["sage_lfq_check"].items() if k != "written_to"}
        self.assertEqual(rec, slc.not_applicable())
        self.assertEqual(rec["why"], "identifications only (quant.lfq is false)")
        os.makedirs(self.out, exist_ok=True)
        with open(os.path.join(self.out, "results.json"), "w") as fh:
            json.dump({"quant": {"lfq": False}}, fh)
        checked = slc.check(self.out)
        self.assertEqual({k: checked[k] for k in ("status", "why", "message")},
                         {k: rec[k] for k in ("status", "why", "message")})
        self.assertIsNone(prov["result"]["mzml_conversion"])
        self.assertIn("identifications.tsv", prov["result"]["note"])

    def test_lfq_true_dotd_still_converts(self):
        """A quantifying cfg (the harness's, quant.lfq true) still converts each .d first."""
        _job, text = self.generate(ids=False)
        self.assertIn(os.path.join(self.out, "mzml", "runA_S1-A1_1_1001.mzML"), text)
        self.assertIn("wine msconvert", text)
        self.assertIn("sage_lfq_check.py", text)
        prov = self.provenance()
        self.assertEqual(prov["sage_input"]["dotd"], "converted")
        self.assertIs(prov["quantified"], True)
        self.assertNotIn("sage_lfq_check", prov, "the job writes the quantifying record itself")

    def test_no_dotd_input_says_so(self):
        mz = os.path.join(self.d, "raw", "ready.mzML")
        open(mz, "w").close()
        self.generate(files=[mz])
        rec = self.provenance()["sage_input"]
        self.assertIsNone(rec["dotd"])
        self.assertEqual(rec["files"], [])


class TheDotdRouteIsDecidedBeforeAnythingIsWritten(_IdsHarness):
    def cfg_with(self, text):
        p = os.path.join(self.d, "hand_written.json")
        with open(p, "w") as fh:
            fh.write(text)
        return p

    def refused(self, files=None, cfg=None, acquisition=None):
        if acquisition:
            with open(self.bundle, "w") as fh:
                json.dump({"acquisition": acquisition,
                           "engine": {"name": "sage", "version": "0.14.7"}}, fh)
        if cfg:
            self.cfg = cfg
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job, "--pwiz-sif", self.sif_file(),
                            files=files or self.dots)
        return job, p

    def sif_file(self):
        sif = os.path.join(self.d, "pwiz.sif")
        open(sif, "w").close()
        return sif

    def assert_nothing_written(self, job, p):
        self.assertNotEqual(p.returncode, 0, p.stdout)
        self.assertIn("Nothing was converted, written or submitted", p.stderr)
        self.assertFalse(os.path.exists(job))
        self.assertFalse(os.path.exists(os.path.join(self.out, "search_provenance.json")))

    def test_an_unset_lfq_with_a_dotd_is_refused(self):
        """L7R-F6: Sage's own default is no LFQ, so 'treated as quantifying' converted every .d
        for quantities Sage then did not make."""
        job, p = self.refused(cfg=self.cfg_with(json.dumps({"database": {}})))
        self.assert_nothing_written(job, p)
        self.assertIn("quant.lfq", p.stderr)
        self.assertIn("--ids-only", p.stderr)

    def test_an_unreadable_cfg_with_a_dotd_is_refused(self):
        job, p = self.refused(cfg=self.cfg_with("{not json"))
        self.assert_nothing_written(job, p)
        self.assertIn("quant.lfq", p.stderr)

    def test_an_unset_lfq_without_a_dotd_still_generates(self):
        """(pin) The rule is about a .d's route; a .raw is converted either way."""
        job, p = self.refused(files=self.raws, cfg=self.cfg_with(json.dumps({"database": {}})))
        self.assertEqual(p.returncode, 0, p.stderr)

    def test_quantifying_a_dia_dotd_is_refused(self):
        """L7R-F7: the conversion's scan summing is DDA-PASEF's; on diaPASEF it would sum each
        window's MS2 over about ten cycles."""
        job, p = self.refused(acquisition="DIA")
        self.assert_nothing_written(job, p)
        self.assertIn("DDA-PASEF", p.stderr)
        self.assertIn("DIA-NN", p.stderr)

    def test_identifications_only_on_a_dia_dotd_is_not_refused(self):
        """(pin) No conversion, so no summing: the IDs-only route is not touched by F7."""
        job, p = self.refused(acquisition="DIA", cfg=self.ids_cfg)
        self.assertEqual(p.returncode, 0, p.stderr)


class AnEarlierSearchsOutputsAreSetAsideAtGeneration(_IdsHarness):
    """L7R-F9: a search generated into a folder an earlier search of the OTHER kind used. Until
    --adapt-only, the watcher showed the earlier LFQ record, and run_de.R could gate on it and
    run the DE on the earlier report.parquet beside a provenance saying quantified: false."""

    def leave(self, files):
        os.makedirs(self.out, exist_ok=True)
        for name, text in files.items():
            with open(os.path.join(self.out, name), "w") as fh:
                fh.write(text)

    def stale(self, name):
        return [f for f in os.listdir(self.out) if f.startswith(name + ".stale-")]

    def test_an_ids_only_search_sets_aside_the_quantities_and_writes_its_record(self):
        self.leave({"report.parquet": "x", "lfq.parquet": "x", "sage_adapt.json": "{}",
                    "sage_lfq_check.json": json.dumps({"status": "warn",
                                                       "runs_outside_lfq_window": ["a"]})})
        self.generate()
        for name in ("report.parquet", "lfq.parquet", "sage_adapt.json"):
            self.assertFalse(os.path.exists(os.path.join(self.out, name)), name)
            self.assertEqual(len(self.stale(name)), 1, name)
        with open(os.path.join(self.out, "sage_lfq_check.json")) as fh:
            self.assertEqual(json.load(fh)["status"], "not_applicable")
        self.assertEqual(len(self.stale("sage_lfq_check.json")), 1)
        moved = self.provenance()["result"]["previous_artefacts_moved_to"]
        self.assertEqual(len(moved), 4, moved)

    def test_a_quantifying_search_sets_aside_identifications(self):
        self.leave({"identifications.tsv": "run\n", "sage_identifications.json": "{}",
                    "sage_lfq_check.json": json.dumps({"status": "not_applicable"})})
        self.generate(ids=False)
        for name in ("identifications.tsv", "sage_identifications.json", "sage_lfq_check.json"):
            self.assertFalse(os.path.exists(os.path.join(self.out, name)), name)
            self.assertEqual(len(self.stale(name)), 1, name)

    def test_a_quantifying_search_keeps_an_earlier_quantifying_record(self):
        """(pin) The job's check carries an earlier acceptance over (carry_acceptance)."""
        rec = {"status": "warn", "accepted": {"reason": "r", "by": "analyst1"}}
        self.leave({"sage_lfq_check.json": json.dumps(rec)})
        self.generate(ids=False)
        with open(os.path.join(self.out, "sage_lfq_check.json")) as fh:
            self.assertEqual(json.load(fh), rec)

    def test_a_refused_generation_moves_nothing(self):
        """(pin) The set-aside comes after every refusal, as DIA-NN's does."""
        self.leave({"identifications.tsv": "run\n"})
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job, "--pwiz-sif", os.path.join(self.d, "none.sif"),
                            files=self.dots)
        self.assertNotEqual(p.returncode, 0)
        self.assertTrue(os.path.exists(os.path.join(self.out, "identifications.tsv")))


@unittest.skipUnless(HAVE_ARROW, "pyarrow is needed to write the Sage fixtures")
class IdsOnlyResults(_IdsHarness):
    def read_ids(self, path):
        with open(path, newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        return rows[0], rows[1:]

    def test_adapt_only_ids_writes_identifications_not_report(self):
        write_ids_results(self.out, [os.path.basename(d) for d in self.dots])
        res = self.adapt_only()
        self.assertIsNone(res["report"])
        self.assertIs(res["quantified"], False)
        self.assertTrue(res["identifications"].endswith("identifications.tsv"))
        self.assertEqual(res["identifications"], os.path.join(self.out, run_search.SAGE_IDS))
        self.assertFalse(os.path.exists(os.path.join(self.out, "report.parquet")))
        header, rows = self.read_ids(res["identifications"])
        self.assertEqual(header, IDS_COLUMNS)
        runs = {r[0] for r in rows}
        self.assertEqual(runs, {"runA_S1-A1_1_1001", "runB_S1-B1_1_1002"})
        # targets at peptide_q <= 0.01 only; one row per (run, peptide, charge), at its best
        # spectrum_q
        lw1 = sorted((r[1], r[3], r[4]) for r in rows if r[0] == "runA_S1-A1_1_1001")
        self.assertEqual(lw1, [("PEPT[+15.9949]IDEK", "2", "0.001"), ("SAMPLER", "2", "0.003"),
                               ("SAMPLER", "3", "0.006")])
        # L7R-F2: a run's identification is a rank-1 PSM at spectrum_q AND peptide_q <= 0.01
        lw2 = {(r[1], r[4]) for r in rows if r[0] == "runB_S1-B1_1_1002"}
        self.assertIn(("SHAREDPEPK", "0.0004"), lw2)
        self.assertNotIn("SHAREDPEPK", {r[1] for r in rows if r[0] == "runA_S1-A1_1_1001"},
                         "spectrum_q 0.30 in this run: not identified here")
        self.assertNotIn("CHIMERAK", {r[1] for r in rows}, "rank 2")
        prots = {r[2] for r in rows}
        self.assertEqual(prots, {"sp|P11111|AAA_HUMAN", "sp|P22222|BBB_HUMAN;sp|P33333|CCC_HUMAN",
                                 "sp|P66666|FFF_HUMAN"})

    def test_the_identifications_record_says_what_was_kept(self):
        write_ids_results(self.out, [os.path.basename(d) for d in self.dots])
        self.generate()                       # the provenance the job's search would have
        self.adapt_only()
        prov = self.provenance()
        rec = prov["sage_identifications"]
        self.assertEqual(rec["identifications"], os.path.join(self.out, run_search.SAGE_IDS))
        self.assertEqual(rec["kept_rule"], run_search.SAGE_IDS_RULE)
        for part in ("rank 1", "spectrum_q <= 0.01", "peptide_q <= 0.01", "is_decoy false"):
            self.assertIn(part, rec["kept_rule"])
        one = rec["per_run"]["runA_S1-A1_1_1001"]
        self.assertEqual((one["psms"], one["decoy_dropped"], one["rank_dropped"],
                          one["spectrum_q_dropped"], one["peptide_q_dropped"],
                          one["psms_kept"], one["precursors"]), (8, 1, 1, 1, 1, 4, 3))
        self.assertEqual(rec["totals"]["peptides"], 3)
        self.assertEqual(prov["sage_lfq_check"]["status"], "not_applicable")
        with open(os.path.join(self.out, "sage_lfq_check.json")) as fh:
            self.assertEqual(json.load(fh)["status"], "not_applicable")

    def test_an_earlier_report_is_set_aside_not_left_beside_identifications(self):
        """A report.parquet left in the folder by an earlier, quantifying search would be read
        by the DE as this search's."""
        write_ids_results(self.out, [os.path.basename(d) for d in self.dots])
        stale = os.path.join(self.out, "report.parquet")
        with open(stale, "w") as fh:
            fh.write("an earlier search's report")
        p = self.run_ids("--adapt-only")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertFalse(os.path.exists(stale))
        moved = [f for f in os.listdir(self.out) if f.startswith("report.parquet.stale-")]
        self.assertEqual(len(moved), 1, os.listdir(self.out))
        self.assertIn(moved[0], p.stderr)

    def test_a_quantified_output_is_not_read_as_identifications(self):
        """The cfg handed to --adapt-only says IDs only, but Sage ran with LFQ (its results.json):
        refused, nothing written, nothing moved."""
        write_ids_results(self.out, [os.path.basename(d) for d in self.dots],
                          quant={"tmt": None, "lfq": True})
        report = os.path.join(self.out, "report.parquet")
        with open(report, "w") as fh:
            fh.write("this search's own report")
        # the quantified search's LFQ record -- a gate run_de.R reads -- must survive
        gate = os.path.join(self.out, "sage_lfq_check.json")
        with open(gate, "w") as fh:
            json.dump({"status": "warn", "runs_outside_lfq_window": ["runA_S1-A1_1_1001"]}, fh)
        p = self.run_ids("--adapt-only")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("quant.lfq", p.stderr)
        self.assertTrue(os.path.exists(report))
        self.assertFalse(os.path.exists(os.path.join(self.out, run_search.SAGE_IDS)))
        with open(gate) as fh:
            self.assertEqual(json.load(fh)["status"], "warn")

    def test_the_job_then_adapt_only(self):
        """The generated job hands Sage the .d paths; --adapt-only then writes the TSV."""
        job, _ = self.generate()
        r = subprocess.run(["bash", job], capture_output=True, text=True, timeout=180,
                           env=self.env(), cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        calls = self.calls()
        self.assertEqual([c["tool"] for c in calls], ["sage"], calls)
        self.assertEqual(calls[0]["inputs"], self.dots)
        self.assertEqual(calls[0]["missing"], [])
        self.assertNotIn("[sage_lfq_check]", r.stderr)
        res = self.adapt_only()
        _, rows = self.read_ids(res["identifications"])
        self.assertEqual(len({r[0] for r in rows}), 2)

    def test_inline_ids_only(self):
        p = self.run_ids()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual([c["tool"] for c in self.calls()], ["sage"])
        self.assertFalse(os.path.exists(os.path.join(self.out, "report.parquet")))
        self.assertTrue(os.path.isfile(os.path.join(self.out, run_search.SAGE_IDS)))
        prov = self.provenance()
        self.assertIs(prov["quantified"], False)
        self.assertEqual(prov["sage_input"]["dotd"], "direct")
        self.assertEqual(prov["sage_lfq_check"]["status"], "not_applicable")
        totals = prov["sage_identifications"]["totals"]
        self.assertEqual((totals["rows"], totals["precursors"], totals["peptides"]), (7, 4, 3))
        self.assertIsNone(prov["result"]["report"])

    def test_the_methods_quote_the_route_and_the_fdr_rule(self):
        """L7R-F3: make_methods quotes what run_search.py recorded (rule 1). An IDs-only search
        no longer says "Identification FDR: not recorded"."""
        import make_methods
        write_ids_results(self.out, [os.path.basename(d) for d in self.dots])
        self.generate()
        self.adapt_only()
        prov = self.provenance()
        route = prov["sage_input"]["methods_sentence"]
        fdr = prov["sage_identifications"]["methods_sentence"]
        self.assertIn("read each timsTOF .d file directly", route)
        self.assertIn("identifications only", route)
        self.assertIn("rank 1", fdr)
        self.assertIn("0.01", fdr)
        rec = make_methods.search_record(params=self.ids_cfg, search_prov=os.path.join(
            self.out, "search_provenance.json"))
        para = make_methods.search_paragraph(rec)
        self.assertIn(route, para)
        self.assertIn(fdr, para)
        self.assertNotIn("Identification FDR", para)
        self.assertIn("identifications only", rec["labelled"]["source"])

    def test_no_identification_is_refused_not_written_empty(self):
        """H2-2 (rule 4): an identifications.tsv with no rows looked like a result. Refused, exit
        2; an earlier identifications.tsv is set aside, so a delivery cannot ship it."""
        write_ids_results(self.out, [])
        stale = os.path.join(self.out, run_search.SAGE_IDS)
        with open(stale, "w") as fh:
            fh.write("an earlier search's identifications\n")
        p = self.run_ids("--adapt-only")
        self.assertEqual(p.returncode, 2, p.stdout + p.stderr)
        self.assertIn("REFUSED", p.stderr)
        self.assertFalse(os.path.exists(stale))
        with open(os.path.join(self.out, run_search.SAGE_IDS_RECORD)) as fh:
            rec = json.load(fh)
        self.assertEqual(rec["status"], "refused")
        self.assertIsNone(rec["identifications"])

    def test_a_quantifying_adapt_of_nothing_identified_is_refused(self):
        """H2-2, as it happened on HIVE: a quantifying search whose Sage found 0 peptides at 1%.
        --adapt-only exited 0 and wrote a 0-row report.parquet with no warning."""
        import test_sage_lfq_check as fx
        fx.write_sage_outputs(self.out, {"HeLa_1.mzML": 0.5, "HeLa_2.mzML": 0.5}, peaks=0,
                              target_q=0.5)
        p = self.run_search("--adapt-only", files=self.raws)      # the harness's lfq-true cfg
        self.assertEqual(p.returncode, 2, p.stdout + p.stderr)
        self.assertIn("no target peptide at 1% FDR", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "report.parquet")))
        with open(os.path.join(self.out, "sage_lfq_check.json")) as fh:
            self.assertEqual(json.load(fh)["status"], "unchecked")              # H2-3

    def test_a_quantifying_adapt_of_ids_output_says_why(self):
        """adapt_sage on an IDs-only output: no lfq.parquet, and the message says the search
        was identifications only rather than guessing at quant.lfq."""
        write_ids_results(self.out, [os.path.basename(d) for d in self.dots])
        p = self.run_search("--adapt-only", files=self.dots)        # the harness's lfq-true cfg
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("identifications only", p.stderr)
        self.assertIn(run_search.SAGE_IDS, p.stderr)


if __name__ == "__main__":
    unittest.main()
