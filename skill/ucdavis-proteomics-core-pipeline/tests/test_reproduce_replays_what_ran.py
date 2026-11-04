#!/usr/bin/env python3
"""The reproducibility bundle replays the analysis that RAN (2.14.2 C2, E2E F2 and F3;
LEDGER Ruling 84).

The E2E run (a two-run timsTOF DDA search with Sage, identifications only, chymotrypsin, no DE)
had finalize build a bundle whose reproduce.sh could not replay it:
  * step 5 searched with `./wf/$(ls wf | grep -vi manifest | head -n1)`, re-derived by
    resolve_defaults.py -- which writes no params file for Sage (or DIA-NN), so the replay was
    `--params ./wf/`, and even with one the chymotrypsin override and quant.lfq false were gone;
  * step 6 ran `run_de.R ... --method <method>` for a DE that never happened (`<method>` is an
    input redirection in bash), and the last line said "Compare ./de_results";
  * REPRODUCE.md described a DE ("The environment the DE ran in -- NOT RECORDED", DPC-Quant's CPU
    family, "methods.txt not found");
  * no search output was fingerprinted (F3, test_bundle_search_fingerprints.py).
What these tests pin:
  * a no-DE Sage identifications-only search with an enzyme override replays exactly its cfg
    (inputs/ holds it byte for byte), and nothing else runs after the search;
  * its REPRODUCE.md and MANIFEST.txt describe the search alone -- the engine, the .d read
    directly, identifications only, how the rows were kept -- and nothing of a DE;
  * a DIA-NN DE session keeps its DE step and DE sections (pin), and searches with its own cfg;
  * a FragPipe search replays the workflow it ran (fragpipe_input.workflow) without its
    database line, which run_search.py sets from --fasta;
  * finalize's command (make_deposit.bundle_command) gives provenance.py the search record and the
    params file the record names (the fingerprints, F3: test_bundle_search_fingerprints.py).
provenance.py runs for real against a fake setup.json (no conda, no R). Stdlib only.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import finish_fixture                                            # noqa: E402
from job_env import job_env                                      # noqa: E402
from test_reproduce_sh_runs import arg, continuation_breaks, run_steps_5_and_6  # noqa: E402

PY = sys.executable

# A Sage cfg as estimate_params.py --ids-only --overrides '{"database": {"enzyme": ...}}' writes
# it: chymotrypsin, no LFQ (identifications only)
SAGE_CFG = {
    "database": {"bucket_size": 8192,
                 "enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 30,
                            "cleave_at": "FYWLM", "restrict": "P"},
                 "static_mods": {"C": 57.0215}, "variable_mods": {"M": [15.994915]},
                 "decoy_tag": "rev_", "generate_decoys": True, "fasta": "REPLACED_AT_RUNTIME.fasta"},
    "precursor_tol": {"ppm": [-20.0, 20.0]}, "fragment_tol": {"ppm": [-20.0, 20.0]},
    "quant": {"lfq": False}, "output_directory": "sage_out"}
SAGE_ROUTE = ("Sage read each timsTOF .d file directly, without conversion to mzML; the search "
              "was for identifications only, with no MS1 quantification.")
SAGE_KEPT = ("a target PSM (is_decoy false) of rank 1 with Sage's spectrum_q <= 0.01 and "
             "peptide_q <= 0.01")
DIANN_CFG = "--qvalue 0.01\n--matrices\n--fasta-search\n--predictor\n--cut K*,R*\n"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


class Base(unittest.TestCase):
    maxDiff = None

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.home = finish_fixture.fake_pipeline_home(self.d)
        self.raws = []
        for r in ("S01_dda_A_1", "S02_dda_B_1"):
            os.makedirs(os.path.join(self.d, "raw", r + ".d"))
            self.raws.append(os.path.join(self.d, "raw", r + ".d"))
        self.out = os.path.join(self.d, "repro")

    def provenance(self, *more):
        r = subprocess.run([PY, os.path.join(SCRIPTS, "provenance.py"), "--outdir", self.out,
                            "--setup-json", os.path.join(self.home, "setup.json"), *more,
                            "--raw", *self.raws],
                           capture_output=True, text=True, timeout=300, env=job_env(self.d))
        self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
        return (read(os.path.join(self.out, "reproduce.sh")),
                read(os.path.join(self.out, "REPRODUCE.md")),
                read(os.path.join(self.out, "MANIFEST.txt")))

    def manifest(self, engine="sage", acq="DDA"):
        return write(os.path.join(self.d, "input", "wf", "workflow.manifest.json"), json.dumps({
            "id": f"{engine}_timstof_{acq.lower()}", "acquisition": acq,
            "instruments": ["timsTOF HT"], "engine": {"name": engine, "version": "0.14.7"
                                                      if engine == "sage" else "2.7.0"},
            "registry": {"source": "local skill defaults", "defaults_version": "2026-09-23"}}))

    def sage_ids_only(self):
        """input/wf/sage_config.json (+ its rationale), output/search/ as an identifications-only
        Sage search leaves it. -> (cfg, record)."""
        cfg = write(os.path.join(self.d, "input", "wf", "sage_config.json"),
                    json.dumps(SAGE_CFG, indent=2))
        write(cfg + ".rationale.json", json.dumps({"overrides": {"database": {
            "value": {"enzyme": {"cleave_at": "FYWLM", "restrict": "P"}},
            "set_by": "the user", "reason": "a chymotrypsin digest"}}}))
        search = os.path.join(self.d, "output", "search")
        ids = {"source": os.path.join(search, "results.sage.parquet"),
               "identifications": os.path.join(search, "identifications.tsv"),
               "kept_rule": SAGE_KEPT, "row_unit": "one row per (run, modified peptide, charge)",
               "quantified": False}
        write(ids["source"], "PAR1-sage-results")
        write(ids["identifications"], "run\tpeptide\tproteins\tcharge\n")
        write(os.path.join(search, "sage_identifications.json"), json.dumps(ids))
        rec = write(os.path.join(search, "search_provenance.json"), json.dumps({
            "engine": "sage", "version": "0.14.7", "params_file": cfg, "files": self.raws,
            "quantified": False, "sage_identifications": ids,
            "sage_input": {"dotd": "direct", "files": self.raws, "methods_sentence": SAGE_ROUTE},
            "result": {"engine": "sage", "out": search, "quantified": False}}))
        return cfg, rec


class SageIdentificationsOnly(Base):
    def bundle(self):
        cfg, rec = self.sage_ids_only()
        files = self.provenance("--workflow-manifest", self.manifest(), "--params", cfg,
                                "--search-record", rec, "--engine", "sage",
                                "--acquisition", "DDA", "--instrument", "timsTOF HT")
        return (cfg, rec) + files

    def test_step_5_searches_with_exactly_the_cfg_that_ran(self):
        cfg, _rec, sh, _md, _man = self.bundle()
        self.assertEqual(continuation_breaks(sh), [])
        self.assertEqual(read(os.path.join(self.out, "inputs", "sage_config.json")), read(cfg),
                         "the cfg that ran, byte for byte")
        calls, _err = run_steps_5_and_6(self.d, self.out, sh)
        self.assertEqual([c["prog"] for c in calls], ["run_search.py"], "nothing after the search")
        s = calls[0]["argv"]
        self.assertEqual(arg(s, "--params"), "./inputs/sage_config.json")
        self.assertEqual(s[s.index("--files") + 1:], self.raws)
        replayed = json.loads(read(os.path.join(self.out, arg(s, "--params"))))
        self.assertEqual(replayed["database"]["enzyme"]["cleave_at"], "FYWLM")
        self.assertEqual(replayed["database"]["enzyme"]["restrict"], "P")
        self.assertIs(replayed["quant"]["lfq"], False, "identifications only, as it ran")
        self.assertNotIn("grep -vi manifest", sh, "never re-derived")

    def test_no_de_step_is_written(self):
        _cfg, _rec, sh, _md, _man = self.bundle()
        for word in ("run_de.R", "<method>", "# 6.", "./de_results"):
            self.assertNotIn(word, sh)
        self.assertIn("No differential expression ran in this analysis, so there is no step 6.",
                      sh)
        self.assertIn('echo "Done. Compare ./search_out against checksums/checksums.json', sh)
        self.assertIn("Identifications only (its record says it did not quantify)", sh)
        self.assertIn("--adapt-only", sh, "how to finish an identifications-only job's replay")

    def test_reproduce_md_describes_the_search_that_ran(self):
        _cfg, _rec, _sh, md, man = self.bundle()
        self.assertIn("## What this analysis was", md)
        self.assertIn("- Search: Sage 0.14.7", md)
        self.assertIn(SAGE_ROUTE, md)
        self.assertIn("Identifications only: the search did not quantify", md)
        self.assertIn(SAGE_KEPT, md)
        self.assertIn("- DE: none ran", md)
        self.assertIn("re-runs the search with the parameters it ran with", md)
        self.assertIn("The publication Methods are the session's `output/methods.md`", md)
        self.assertIn("- registry: local skill defaults (defaults table `2026-09-23`)", md)
        for word in ("The environment the DE ran in", "DPC-Quant", "methods.txt not found",
                     "re-runs search + DE", "method=`None`", "None @ commit"):
            self.assertNotIn(word, md)
        self.assertIn("[NOTE]    no differential expression ran", man)
        for word in ("[SKIPPED] DE runtime", "[SKIPPED] reproducibility_log.R",
                     "default contrasts"):
            self.assertNotIn(word, man)
        self.assertIn("[OK]      step 5 searches with the parameters the search was given "
                      "(inputs/sage_config.json)", man)
        rm = json.loads(read(os.path.join(self.out, "run_manifest.json")))
        self.assertIs(rm["de"]["ran"], False)
        self.assertEqual(rm["search_record"], os.path.join(self.d, "output", "search",
                                                           "search_provenance.json"))


class DiannDeSession(Base):
    def bundle(self):
        cfg = write(os.path.join(self.d, "input", "wf", "diann.cfg"), DIANN_CFG)
        search = os.path.join(self.d, "output", "search")
        report = write(os.path.join(search, "report.parquet"), "PAR1-report")
        write(os.path.join(search, "search_provenance.json"), json.dumps({
            "engine": "diann", "version": "2.7.0", "params_file": cfg, "files": self.raws,
            "result": {"engine": "diann", "report": report}}))
        de = os.path.join(self.d, "output", "tables")
        write(os.path.join(de, "de_provenance.json"), json.dumps({
            "pipeline_id": "dpc", "method": "dpc", "contrasts": ["B-A"],
            "quantities": "normalised", "q_cutoff": 0.01, "logfc": 1, "adjp": 0.05,
            "compute": {"cpu_model": "AMD EPYC 9734", "cpu_family": "zen4",
                        "slurm_node": "node-1", "blas": "libopenblas", "lapack": "libopenblas"}}))
        write(os.path.join(de, "DE_B_vs_A.csv"), "Protein,logFC\nP1,1\n")
        conds = write(os.path.join(self.d, "input", "conditions.csv"),
                      "File.Name,Group\nS01_dda_A_1,A\nS02_dda_B_1,B\n")
        return cfg, self.provenance(
            "--workflow-manifest", self.manifest("diann", "DIA"), "--params", cfg,
            "--report", report, "--engine", "diann", "--de-dir", de, "--de-method", "dpc",
            "--conditions", conds, "--q-cutoff", "0.01", "--logfc", "1", "--adjp", "0.05",
            "--acquisition", "DIA", "--instrument", "timsTOF HT")

    def test_the_de_step_and_its_sections_are_as_before(self):
        """(pin) run_de.R with the DE's own method and contrasts; the DE's REPRODUCE.md sections."""
        _cfg, (sh, md, man) = self.bundle()
        self.assertEqual(continuation_breaks(sh), [])
        calls, _err = run_steps_5_and_6(self.d, self.out, sh)
        self.assertEqual([c["prog"] for c in calls], ["run_search.py", "Rscript"])
        de = calls[1]["argv"]
        self.assertTrue(de[0].endswith("/scripts/run_de.R"), de)
        self.assertEqual((arg(de, "--method"), arg(de, "--contrasts"), arg(de, "--metadata")),
                         ("dpc", "B-A", "inputs/conditions.csv"))
        self.assertIn("# 6. Re-run differential expression with identical settings.", sh)
        self.assertIn('echo "Done. Compare ./de_results against checksums/checksums.json', sh)
        for section in ("## The environment the DE ran in",
                        "## Exact numbers need the same kind of CPU",
                        "- DE: method=`dpc`, contrasts=`B-A`", "re-runs search + DE"):
            self.assertIn(section, md)
        self.assertNotIn("## What this analysis was", md)

    def test_the_search_replays_its_own_cfg(self):
        cfg, (sh, _md, _man) = self.bundle()
        calls, _err = run_steps_5_and_6(self.d, self.out, sh)
        self.assertEqual(arg(calls[0]["argv"], "--params"), "./inputs/diann.cfg")
        self.assertEqual(read(os.path.join(self.out, "inputs", "diann.cfg")), read(cfg))
        rm = json.loads(read(os.path.join(self.out, "run_manifest.json")))
        self.assertIs(rm["de"]["ran"], True)


class FragPipeSearch(Base):
    def test_the_workflow_fragpipe_ran_is_replayed_without_its_database(self):
        search = os.path.join(self.d, "output", "search")
        given = write(os.path.join(self.d, "input", "wf", "fragpipe.workflow"),
                      "msfragger.search_enzyme_name_1=chymotrypsin\n"
                      "msfragger.search_enzyme_cut_1=FLWY\nmsfragger.search_enzyme_nocut_1=P\n")
        ran = write(os.path.join(search, "fragpipe_input.workflow"),
                    "msfragger.search_enzyme_name_1=chymotrypsin\n"
                    "database.db-path=/orig/session/input/search.fasta\n"
                    "msfragger.search_enzyme_cut_1=FLWY\nmsfragger.search_enzyme_nocut_1=P\n")
        rec = write(os.path.join(search, "search_provenance.json"), json.dumps({
            "engine": "fragpipe", "version": "24.0", "params_file": given, "files": self.raws,
            "quantified": False,
            "fragpipe_workflow": {"given": given, "run": ran},
            "result": {"engine": "fragpipe", "fragpipe_workflow": {"given": given, "run": ran}}}))
        sh, _md, man = self.provenance("--workflow-manifest", self.manifest("fragpipe"),
                                       "--params", given, "--search-record", rec,
                                       "--engine", "fragpipe")
        self.assertEqual(continuation_breaks(sh), [])
        calls, _err = run_steps_5_and_6(self.d, self.out, sh)
        self.assertEqual(arg(calls[0]["argv"], "--params"), "./inputs/fragpipe_input.workflow")
        self.assertEqual(read(os.path.join(self.out, "inputs", "fragpipe_input.workflow")),
                         "msfragger.search_enzyme_name_1=chymotrypsin\n"
                         "msfragger.search_enzyme_cut_1=FLWY\nmsfragger.search_enzyme_nocut_1=P\n")
        self.assertTrue(os.path.isfile(os.path.join(self.out, "inputs", "fragpipe.workflow")),
                        "the workflow it was given is kept beside it")
        self.assertIn("copied the workflow FragPipe ran (fragpipe_input.workflow)", man)


class FinalizesCommand(Base):
    """make_deposit.bundle_command -- the command finalize runs -- for an identifications-only
    Sage session: the search record, and the params file it names."""

    def session(self):
        import test_deposit_package as tdp
        root = os.path.join(self.d, "s")                 # Base made <d>/raw already
        os.makedirs(root)
        p = tdp.dia_session(root, de=False)
        for n in ("report.parquet", "report.log.txt", "report.pg_matrix.tsv"):
            os.remove(os.path.join(p["search_out"], n))
        os.remove(os.path.join(p["workflow_dir"], "params.cfg"))
        cfg = write(os.path.join(p["workflow_dir"], "sage_config.json"), json.dumps(SAGE_CFG))
        # an older cfg the session's own guess (session.params_file) would take first
        write(os.path.join(p["workflow_dir"], "params.json"), json.dumps({"older": True}))
        search = p["search_out"]
        ids = {"source": os.path.join(search, "results.sage.parquet"),
               "identifications": os.path.join(search, "identifications.tsv"),
               "kept_rule": SAGE_KEPT, "row_unit": "one row per (run, modified peptide, charge)"}
        write(ids["source"], "PAR1")
        write(ids["identifications"], "run\tpeptide\n")
        write(os.path.join(search, "sage_identifications.json"), json.dumps(ids))
        write(p["search_prov"], json.dumps({
            "engine": "sage", "version": "0.14.7", "params_file": cfg, "quantified": False,
            "sage_identifications": ids, "sage_input": {"dotd": "direct",
                                                        "methods_sentence": SAGE_ROUTE},
            "result": {"engine": "sage", "quantified": False}}))
        write(p["workflow_manifest"], json.dumps({
            "acquisition": "DDA", "instruments": ["timsTOF HT"],
            "engine": {"name": "sage", "version": "0.14.7"}}))
        return p, cfg, ids

    def test_the_command_names_the_record_and_the_params_it_names(self):
        import make_deposit
        p, cfg, ids = self.session()
        cmd, missing, unnamed = make_deposit.bundle_command(make_deposit.gather(p["session_dir"]))
        self.assertEqual(missing, [])
        self.assertEqual(arg(cmd, "--params"), cfg, "the file the record says the search was given")
        self.assertEqual(arg(cmd, "--search-record"), p["search_prov"])
        self.assertNotIn("--report", cmd, "an identifications-only search has no report")

    def test_finalize_builds_a_bundle_that_replays_it(self):
        import make_deposit
        p, cfg, ids = self.session()
        env = job_env(self.d, PROTEOMICS_PIPELINE_HOME=self.home)
        with mock.patch.dict(os.environ, env, clear=True):
            man = make_deposit.Manifest()
            self.assertTrue(make_deposit.ensure_bundle(p["session_dir"], man), man.lines)
        repro = p["repro_dir"]
        sh = read(os.path.join(repro, "reproduce.sh"))
        self.assertIn("--params ./inputs/sage_config.json", sh)
        self.assertNotIn("run_de.R", sh)
        self.assertEqual(read(os.path.join(repro, "inputs", "sage_config.json")), read(cfg))
        md = read(os.path.join(repro, "REPRODUCE.md"))
        self.assertIn("Identifications only: the search did not quantify", md)
        self.assertNotIn("The environment the DE ran in", md)


if __name__ == "__main__":
    unittest.main()
