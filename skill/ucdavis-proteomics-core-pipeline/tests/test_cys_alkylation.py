#!/usr/bin/env python3
"""Cysteine that was never alkylated must not be searched as carbamidomethyl.

Staff reports (2026-10-09): estimate_params.py always wrote --unimod4 (fixed carbamidomethyl C),
and --overrides could only add or replace a flag, so on-bead / IP / pull-down digests with no
reduction and alkylation were searched with every cysteine peptide unsearchable -- 0.03-0.5% of
identified peptides held a C, against 7.8% on an S-Trap prep with DTT/IAA. The fix used was a
hand `sed` of params.cfg plus a hand-edited rationale sidecar.

  * estimate_params.py --cys-alkylation cam|none: none leaves --unimod4 out (Sage: no static C),
    and the sidecar records which, and on whose word -- cam unasked is tagged DEFAULT;
  * the Methods name carbamidomethylation only when it was searched, and say why it was not;
  * sample_quality.py --report warns when carbamidomethyl was fixed and under 1% of the identified
    precursors hold a C.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import estimate_params as ep  # noqa: E402
import make_methods as mm  # noqa: E402
import sample_quality as sq  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_PYARROW = True
except ImportError:
    HAVE_PYARROW = False


def estimate(d, engine="diann", *more):
    out = os.path.join(d, "params.cfg" if engine == "diann" else "sage_config.json")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                        "--engine", engine, "--acquisition", "DIA", "--instrument", "timsTOF HT",
                        "--precursor-mz-range", "300", "1200", "--out", out, *more],
                       capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError(p.stderr)
    with open(out) as fh, open(out + ".rationale.json") as side:
        return out, fh.read(), json.load(side)


class EstimateParams(unittest.TestCase):
    def test_the_default_is_fixed_carbamidomethyl_tagged_as_a_default(self):
        with tempfile.TemporaryDirectory() as d:
            _, cfg, side = estimate(d)
            self.assertIn("--unimod4\n", cfg)
            self.assertEqual(side["cys_alkylation"]["value"], "cam")
            self.assertFalse(side["cys_alkylation"]["given"])
            self.assertTrue(side["cys_alkylation"]["source"].startswith("DEFAULT, not user-confirmed"))
            self.assertIn("DEFAULT", side["rationale"]["--unimod4"]["source"])

    def test_none_leaves_cysteine_unmodified_and_says_who_said_so(self):
        with tempfile.TemporaryDirectory() as d:
            path, cfg, side = estimate(d, "diann", "--cys-alkylation", "none")
            self.assertNotIn("--unimod4", cfg)
            self.assertEqual(side["cys_alkylation"]["value"], "none")
            self.assertTrue(side["cys_alkylation"]["given"])
            self.assertIn("not reduced and alkylated", side["cys_alkylation"]["source"])
            self.assertIn("flag omitted", side["rationale"]["--unimod4"]["value"])
            self.assertFalse(ep.cam_searched(path))
            _, cfg, side = estimate(d, "diann", "--cys-alkylation", "cam")
            self.assertIn("--unimod4\n", cfg)
            self.assertTrue(side["cys_alkylation"]["given"])

    def test_a_yes_is_recorded_as_the_users_answer_not_a_default(self):
        """Staff review of 2.12: step 3's "yes, reduced and alkylated" must be passed as
        --cys-alkylation cam, and recorded as given -- only an unasked search is a DEFAULT."""
        with tempfile.TemporaryDirectory() as d:
            _, cfg, side = estimate(d, "diann", "--cys-alkylation", "cam")
            self.assertIn("--unimod4\n", cfg)
            self.assertEqual(side["cys_alkylation"]["value"], "cam")
            self.assertTrue(side["cys_alkylation"]["given"])
            self.assertTrue(side["cys_alkylation"]["source"].startswith("given (--cys-alkylation cam)"))
            self.assertNotIn("DEFAULT", side["rationale"]["--unimod4"]["source"])
            self.assertNotIn("DEFAULT", side["cys_alkylation"]["source"])

    def test_sage_none_has_no_static_cysteine(self):
        with tempfile.TemporaryDirectory() as d:
            path, text, side = estimate(d, "sage", "--cys-alkylation", "none")
            self.assertEqual(json.loads(text)["database"]["static_mods"], {})
            self.assertIn("unmodified", side["rationale"]["static_mods"]["source"])
            self.assertFalse(ep.cam_searched(path))
            path, text, _ = estimate(d, "sage")
            self.assertEqual(json.loads(text)["database"]["static_mods"], {"C": 57.0215})
            self.assertTrue(ep.cam_searched(path))

    def test_cam_searched_reads_every_way_a_cfg_can_fix_it(self):
        with tempfile.TemporaryDirectory() as d:
            def cfg(text):
                p = os.path.join(d, "x.cfg")
                with open(p, "w") as fh:
                    fh.write(text)
                return p
            self.assertTrue(ep.cam_searched(cfg("--qvalue 0.01\n--unimod4\n")))
            self.assertTrue(ep.cam_searched(cfg("--fixed-mod UniMod:4,57.021464,C\n")))
            self.assertTrue(ep.cam_searched(cfg("--fixed-mod Carbamidomethyl,57.021464,C\n")))
            self.assertFalse(ep.cam_searched(cfg("--qvalue 0.01\n# --unimod4\n")))
            self.assertIsNone(ep.cam_searched(os.path.join(d, "missing.cfg")))


class Methods(unittest.TestCase):
    def paragraph(self, d, *more):
        path, _, _ = estimate(d, "diann", *more)
        prov = os.path.join(d, "search_provenance.json")
        with open(prov, "w") as fh:
            json.dump({"engine": "diann", "version": "2.7.0", "params_file": path,
                       "result": {}}, fh)
        return mm.search_paragraph(mm.search_record(path, prov))

    def test_carbamidomethylation_is_named_only_when_it_was_searched(self):
        with tempfile.TemporaryDirectory() as d:
            text = self.paragraph(d)
            self.assertIn("Carbamidomethyl", text)
            self.assertNotIn("searched unmodified", text)
        with tempfile.TemporaryDirectory() as d:
            text = self.paragraph(d, "--cys-alkylation", "none")
            self.assertNotIn("Carbamidomethyl (", text)
            self.assertIn("No fixed modifications were set.", text)
            self.assertIn("Cysteine was searched unmodified (no carbamidomethylation), as the "
                          "samples were not reduced and alkylated.", text)


def write_report(path, seqs):
    """A DIA-NN-shaped report: one precursor per sequence, all at 1% FDR."""
    pq.write_table(pa.table({"Run": ["r1"] * len(seqs), "Protein.Group": ["P1"] * len(seqs),
                             "Precursor.Id": [f"{s}2" for s in seqs],
                             "Stripped.Sequence": seqs, "Q.Value": [0.001] * len(seqs)}), path)


def write_runs_report(path, runs):
    """{run: [sequences]} -> a DIA-NN-shaped report, every precursor at 1% FDR."""
    rows = [(r, s) for r, seqs in runs.items() for s in seqs]
    pq.write_table(pa.table({"Run": [r for r, _ in rows], "Protein.Group": ["P1"] * len(rows),
                             "Precursor.Id": [f"{s}2" for _, s in rows],
                             "Stripped.Sequence": [s for _, s in rows],
                             "Q.Value": [0.001] * len(rows)}), path)


@unittest.skipUnless(HAVE_PYARROW, "the cysteine check reads report.parquet with pyarrow")
class SampleQualityPerRun(unittest.TestCase):
    """Statistics review of 2.12: a cohort half of whose samples were never alkylated stays above
    1% pooled, so each run is judged too."""

    def search(self, d, runs):
        params, _, _ = estimate(d, "diann", "--cys-alkylation", "cam")
        report = os.path.join(d, "report.parquet")
        write_runs_report(report, runs)
        with open(os.path.join(d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "params_file": params}, fh)
        return report

    def test_one_unalkylated_run_warns_though_the_pool_is_above_1pc(self):
        alkylated = [f"PEPCTIDE{i}" for i in range(16)] + [f"PEPTIDE{i}" for i in range(184)]
        not_alkylated = [f"PEPTIDEK{i}" for i in range(200)]
        few = [f"PEPTIDER{i}" for i in range(20)]                  # no C, but only 20 IDs
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, {"r1": alkylated, "r2": not_alkylated,
                                               "blank": few}))
            self.assertGreater(rec["fraction"], 0.01, "the pooled share hides it")
            self.assertEqual(rec["status"], "WARN")
            self.assertEqual(rec["low_runs"], ["r2"])
            self.assertEqual(rec["too_few_to_judge"], ["blank"])
            self.assertEqual(rec["per_run"]["r1"], {"n": 200, "n_cys": 16, "fraction": 0.08})
            self.assertIn("1 of 3 runs have under 1%", rec["message"])
            self.assertIn("r2 0.00% (0 of 200)", rec["message"])
            self.assertIn("mixed prep", rec["message"])
            self.assertIn("Not judged alone (under 100 identified precursors): blank", rec["message"])

    def test_the_warning_names_other_agents_and_cysteine_poor_samples(self):
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, {"r1": [f"PEPTIDEK{i}" for i in range(200)]}))
            self.assertEqual(rec["status"], "WARN")
            for words in ("the search fixed Carbamidomethyl (+57.0215) on every C",
                          "NOT reduced and alkylated", "--cys-alkylation none",
                          "alkylated with another agent than the one searched",
                          "nem (N-ethylmaleimide (NEM))", "propionamide (acrylamide)",
                          "carboxymethyl (iodoacetic acid)", "pyridylethyl (4-vinylpyridine)",
                          "mmts (S-methyl methanethiosulfonate (MMTS))", "histone or collagen"):
                self.assertIn(words, rec["message"])
            self.assertNotIn("cam (", rec["message"], "the searched agent is not offered")

    def test_every_run_alkylated_is_info(self):
        alk = [f"PEPCTIDE{i}" for i in range(16)] + [f"PEPTIDE{i}" for i in range(184)]
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, {"r1": alk, "r2": [x + "Z" for x in alk]}))
            self.assertEqual(rec["status"], "INFO")
            self.assertEqual(rec["low_runs"], [])


@unittest.skipUnless(HAVE_PYARROW, "the cysteine check reads report.parquet with pyarrow")
class SampleQuality(unittest.TestCase):
    # 1 of 200 precursors with a C: 0.5%, the top of what the on-bead searches gave
    LOW = ["PEPTCDEK"] + [f"PEPTIDE{'K' * (i % 5 + 1)}{i}" for i in range(199)]
    NORMAL = [f"PEPCTIDE{i}" for i in range(16)] + [f"PEPTIDE{i}" for i in range(184)]

    def search(self, d, seqs, cys="cam"):
        params, _, _ = estimate(d, "diann", "--cys-alkylation", cys)
        report = os.path.join(d, "report.parquet")
        write_report(report, seqs)
        with open(os.path.join(d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "params_file": params}, fh)
        return report

    def test_fixed_cam_with_almost_no_cysteine_warns_and_says_the_fix(self):
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, self.LOW))
            self.assertEqual(rec["status"], "WARN")
            self.assertTrue(rec["cam_searched"])
            self.assertEqual((rec["n_cys"], rec["n"]), (1, 200))
            self.assertIn("NOT reduced and alkylated", rec["message"])
            self.assertIn("--cys-alkylation none", rec["message"])

    def test_a_normal_share_or_an_unmodified_search_does_not_warn(self):
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, self.NORMAL))
            self.assertEqual(rec["status"], "INFO")
            self.assertAlmostEqual(rec["fraction"], 0.08)
        with tempfile.TemporaryDirectory() as d:
            rec = sq.cys_check(self.search(d, self.LOW, cys="none"))
            self.assertEqual(rec["status"], "INFO")
            self.assertFalse(rec["cam_searched"])
            self.assertIn("left cysteine unmodified", rec["message"])

    def test_no_params_found_is_said_not_guessed(self):
        with tempfile.TemporaryDirectory() as d:
            report = os.path.join(d, "report.parquet")
            write_report(report, self.LOW)
            rec = sq.cys_check(report)
            self.assertEqual(rec["status"], "INFO")
            self.assertIsNone(rec["cam_searched"])
            self.assertIn("pass --params", rec["message"])

    def test_the_warning_reaches_sample_quality_md(self):
        with tempfile.TemporaryDirectory() as d:
            report = self.search(d, self.LOW)
            matrix = os.path.join(d, "Expression_Matrix.csv")
            with open(matrix, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["Protein.Group", "Genes", "Protein.Names", "S1", "S2", "S3"])
                w.writerow(["P1", "ACTB", "ACTB_HUMAN", 16, 16, 16])
            out = os.path.join(d, "SAMPLE_QUALITY.md")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "sample_quality.py"),
                                "--matrix", matrix, "--report", report, "--out", out],
                               capture_output=True, text=True, cwd=d)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                md = fh.read()
            self.assertIn("## Cysteine", md)
            self.assertIn("- CYSTEINE: only 0.50% of the identified precursors", md)
            self.assertEqual(json.loads(p.stdout)["cysteine"]["status"], "WARN")


class Reagents(unittest.TestCase):
    """Review of 2.12 (proteomics, P3): one map of alkylating agents, each with its own fixed C
    modification (UniMod, read 2026-10-15) -- the Core has seen acrylamide (propionamide)."""

    def test_the_map_is_unimods(self):
        want = {"cam": (4, 57.021464), "nem": (108, 125.047679), "carboxymethyl": (6, 58.005479),
                "propionamide": (24, 71.037114), "mmts": (39, 45.987721),
                "pyridylethyl": (31, 105.057849)}
        self.assertEqual({k: (r["unimod"], r["mass"]) for k, r in ep.CYS_REAGENTS.items() if r},
                         want)
        self.assertIsNone(ep.CYS_REAGENTS["none"])

    def test_each_engine_gets_the_agents_modification(self):
        with tempfile.TemporaryDirectory() as d:
            path, cfg, side = estimate(d, "diann", "--cys-alkylation", "propionamide")
            self.assertIn("--fixed-mod UniMod:24,71.037114,C\n", cfg)
            self.assertNotIn("--unimod4", cfg)
            self.assertEqual(side["cys_alkylation"]["reagent"], "acrylamide")
            self.assertIn("alkylated with acrylamide", side["cys_alkylation"]["source"])
            self.assertEqual(ep.cys_searched(path), {"fixed": True, "choice": "propionamide",
                                                     "mass": 71.037114, "mod": "Propionamide"})
            self.assertFalse(ep.cam_searched(path), "another agent is not carbamidomethyl")
            path, text, _ = estimate(d, "sage", "--cys-alkylation", "nem")
            self.assertEqual(json.loads(text)["database"]["static_mods"], {"C": 125.047679})
            self.assertEqual(ep.cys_searched(path)["choice"], "nem")
            path, _, _ = estimate(d, "diann", "--cys-alkylation", "none")
            self.assertEqual(ep.cys_searched(path), {"fixed": False})

    def test_the_methods_name_the_agents_modification(self):
        with tempfile.TemporaryDirectory() as d:
            text = Methods().paragraph(d, "--cys-alkylation", "propionamide")
            self.assertIn("Fixed modifications: Propionamide (C).", text)
            self.assertIn("Cysteines were alkylated with acrylamide.", text)
            self.assertNotIn("Carbamidomethyl", text)
            self.assertNotIn("name not verified", text)

    @unittest.skipUnless(HAVE_PYARROW, "the cysteine check reads report.parquet with pyarrow")
    def test_a_search_for_the_wrong_agent_warns_with_the_others_offered(self):
        with tempfile.TemporaryDirectory() as d:
            params, _, _ = estimate(d, "diann", "--cys-alkylation", "nem")
            report = os.path.join(d, "report.parquet")
            write_report(report, [f"PEPTIDEK{i}" for i in range(200)])
            rec = sq.cys_check(report, params)
            self.assertEqual(rec["status"], "WARN")
            self.assertIn("fixed Nethylmaleimide (+125.048) on every C", rec["message"])
            self.assertIn("cam (iodoacetamide or chloroacetamide)", rec["message"])
            self.assertNotIn("nem (", rec["message"])


class UnreadableReport(unittest.TestCase):
    """Error audit of 2.12: a damaged report.parquet made cys_check() and detected_depth() raise
    with a traceback. Each now says, in one line, that it could not be read."""

    @unittest.skipUnless(HAVE_PYARROW, "the cysteine check reads report.parquet with pyarrow")
    def test_a_damaged_report_is_one_info_line(self):
        with tempfile.TemporaryDirectory() as d:
            rep = os.path.join(d, "report.parquet")
            with open(rep, "wb") as fh:
                fh.write(b"PAR1 this is not parquet")
            rec = sq.cys_check(rep)
            self.assertEqual(rec["status"], "INFO")
            self.assertRegex(rec["message"], r"^cysteine check not run: report\.parquet could not "
                                             r"be read \(\w+\)$")
            depth, note = sq.detected_depth(rep)
            self.assertIsNone(depth)
            self.assertRegex(note, r"^report\.parquet could not be read \(\w+\)$")


class OtherRoutes(unittest.TestCase):
    """Review of 2.12 (proteomics, P2): Radiant never searched carbamidomethyl-C -- its DIA-NN
    library was predicted without --unimod4 and its config has no modifications (HIVE, DIA-NN
    2.7.0, a test job: 0 of 6,044 C-peptides carried UniMod:4) -- and FragPipe's presets fixed
    carbamidomethyl whatever --cys-alkylation said."""

    def presets(self, d, engine, *more):
        out = os.path.join(d, "wf", f"x.{'workflow' if engine == 'fragpipe' else 'radiantConfig'}")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_presets.py"), "--engine",
                            engine, "--instrument", "timsTOF HT" if engine == "fragpipe" else
                            "Orbitrap Astral", "--out", out, *more],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out) as fh:
            return out, fh.read(), json.loads(p.stdout)

    def test_fragpipe_fixes_the_agents_mass_on_c(self):
        with tempfile.TemporaryDirectory() as d:
            _, text, prov = self.presets(d, "fragpipe", "--cys-alkylation", "propionamide")
            row = next(ln for ln in text.splitlines() if ln.startswith("msfragger.table.fix-mods="))
            self.assertIn("; 71.037114,C (cysteine),true,-1;", row)
            self.assertNotIn("57.02146,C (cysteine)", row)
            self.assertEqual(prov["cys_alkylation"]["value"], "propionamide")
            _, text, _ = self.presets(d, "fragpipe", "--cys-alkylation", "none")
            self.assertIn("; 0.0,C (cysteine),true,-1;", text)
            _, text, prov = self.presets(d, "fragpipe")
            self.assertIn("57.02146,C (cysteine)", text, "cam keeps the vendor row as shipped")
            self.assertFalse(prov["cys_alkylation"]["given"])

    def test_radiant_libraries_carry_the_fixed_c(self):
        import run_search
        with tempfile.TemporaryDirectory() as d:
            cfg, _, prov = self.presets(d, "radiant", "--cys-alkylation", "nem")
            self.assertIn("--fixed-mod UniMod:108,125.047679,C",
                          prov["keys_changed"]["cysteine (DIA-NN library prediction)"])
            self.assertEqual(run_search.radiant_cys(cfg), "nem")
            self.assertEqual(run_search.radiant_cys(os.path.join(d, "no_record.radiantConfig")),
                             "cam")
            lib, lines = run_search.radiant_library_job({"diann": "/opt/diann"}, "/x/db.fasta",
                                                        os.path.join(d, "out"), 4, cfg)
            self.assertIn("--cys-alkylation nem", " ".join(lines))
            # the array chain's step 1 predicts with the flags
            runs = [os.path.join(d, f"r{i}.mzML") for i in range(2)]
            for r in runs:
                open(r, "w").close()
            listing = os.path.join(d, "list.txt")
            with open(listing, "w") as fh:
                fh.write("\n".join(runs) + "\n")
            fasta = os.path.join(d, "db.fasta")
            with open(fasta, "w") as fh:
                fh.write(">sp|P1|X\nPEPTIDEK\n")
            out = os.path.join(d, "rp")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "radiant_parallel.py"),
                                "--runtime", "apptainer", "--image", "/x/r.sif", "--raw-list",
                                listing, "--fasta", fasta, "--config", cfg, "--out", out,
                                "--diann", "/opt/diann", "--cys-alkylation", "cam",
                                "--partition", "high", "--account", "genome-center-grp"],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(os.path.join(out, "radiant_parallel", "step1_libpred.sbatch")) as fh:
                self.assertIn("--threads 16 --unimod4", fh.read())

    def test_a_radiant_library_is_reused_only_for_the_same_choice(self):
        with tempfile.TemporaryDirectory() as d:
            tsv = os.path.join(d, "radiant_library.tsv")
            with open(tsv, "w") as fh:
                fh.write("x" * 1_100_000)
            argv = [sys.executable, os.path.join(SCRIPTS, "make_radiant_library.py"),
                    "--out-dir", d, "--diann", "/nonexistent/diann", "--fasta", "/x.fasta"]
            # no record: a pre-2.12 library, predicted with nothing fixed on C -- rebuilt
            p = subprocess.run(argv, capture_output=True, text=True)
            self.assertNotIn('"reused": true', p.stdout)
            with open(tsv + ".meta.json", "w") as fh:
                json.dump({"cys_alkylation": "cam"}, fh)
            p = subprocess.run(argv, capture_output=True, text=True)
            self.assertEqual(json.loads(p.stdout)["reused"], True)
            p = subprocess.run(argv + ["--cys-alkylation", "none"], capture_output=True, text=True)
            self.assertNotIn('"reused": true', p.stdout)


if __name__ == "__main__":
    unittest.main()
